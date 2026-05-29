"""Commit-check tool for the PHYEX offline testprogs.

Compiles the PHYEX package using a specific commit, runs the different
test programs and checks whether results are identical to a given version.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

from pyphyextools.testing.check_common import (
    CheckCommitBase, CheckCommitError, escape_commit, PHYEXCONF_DEFAULT, run_tool)
from pyphyextools.testing.compare import comp_testprogs


class CheckCommitTestprogs(CheckCommitBase):
    """Check a commit against the offline testprogs reference."""

    default_expand = True
    default_buildSystem = 'fcm'
    alternative_buildSystem = 'ecbuild'

    @staticmethod
    def _extrapolation_configs():
        """Return ``(conf_extra_tag, conf_extra_opts)`` dicts for extrapolation configurations."""
        conf_extra_tag = {
            0: '',
            1: '_Z120_NPRO32_BLK1024',
            2: '_Z120_NPRO32_BLK256_TIMES4',
            3: '_Z120_NPRO32_BLK64_TIMES16',
            4: '_Z120_NPRO${NPROMA}_BLK${NBLOCKS}_TIMES${NTIMES}',
        }
        conf_extra_opts = {
            0: '',
            1: '--nflevg 120 --nproma 32 --blocks 1024',
            2: '--nflevg 120 --nproma 32 --blocks 256 --times 4',
            3: '--nflevg 120 --nproma 32 --blocks 64 --times 16',
            4: '--nflevg 120 --nproma ${NPROMA} --blocks ${NBLOCKS} --times ${NTIMES}',
        }
        return conf_extra_tag, conf_extra_opts

    def __init__(self, archfile=None, refarchfile=None, extrapolation=0,
                 perf=True, checkOpt='--check', buildSys=None, **kwargs):
        super().__init__(**kwargs)
        self.defaultTest = ["ice_adjust", "rain_ice", "rain_ice_old",
                            "turb", "shallow", "lima_adjust", "lima"]
        self.allowedTests = ["ice_adjust", "rain_ice", "rain_ice_old",
                             "turb", "shallow", "lima_adjust", "lima"]
        self.buildSys = buildSys or self.default_buildSystem
        self.perf = perf
        self.checkOpt = checkOpt
        self.conf_extra_tag, self.conf_extra_opts = self._extrapolation_configs()
        self.archfile = archfile
        self.refarchfile = refarchfile
        self.extrapolation = extrapolation
        self.TESTDIR = os.environ.get('TESTPROGSDIR',
                                       os.path.expanduser('~/TESTPROGS'))
        self.dirdata = os.path.join(PHYEXCONF_DEFAULT, 'testprogs_data')
        self.defaultarchfile = 'gnu'
        self.submit_method = ''
        self.varToExport = ("NPROMA,NBLOCKS,OMP_NUM_THREADS,DR_HOOK_OPT,"
                            "DR_HOOK,DR_HOOK_IGNORE_SIGNALS,"
                            "NVCOMPILER_ACC_GANGLIMIT")
        self.extrapolation_tag = ''
        self.extrapolation_opts = ''

        hostname = os.uname().nodename
        if hostname[:7] in ('belenos', 'taranis'):
            self.defaultarchfile = 'MIMPIIFC1805.EPONA'
            self.submit_method = 'slurm_belenos'
        elif hostname == 'aurora01':
            self.defaultarchfile = 'ECMWF_NEC440MPI225SP.AU.x'
            self.submit_method = ''
        else:
            self.defaultarchfile = 'gnu'
            self.submit_method = ''

        self.archfile = archfile or self.defaultarchfile
        self.refarchfile = refarchfile or self.defaultarchfile

        tag = self.conf_extra_tag.get(self.extrapolation, '')
        self.extrapolation_tag = os.path.expandvars(tag) if tag else ''
        opts = self.conf_extra_opts.get(self.extrapolation, '')
        self.extrapolation_opts = os.path.expandvars(opts) if opts else ''

        if os.path.isdir(os.path.join(self.commit, 'src')):
            self.model_ready = False
            json_path = os.path.join(self.commit, 'src', 'offline', 'testprogs_version.json')
        else:
            self.model_ready = True
            json_path = os.path.join(self.commit, 'testprogs_version.json')
        if os.path.isfile(json_path):
            with open(json_path, encoding='utf-8') as f:
                self.json_content = json.load(f)
        else:
            self.json_content = {}

        testing = self.json_content.get('testing', '')
        all_tests = []
        if testing:
            for t in self.allowedTests:
                ref = testing.get(t, '')
                if ref:
                    all_tests.append(t)
        if not self.tests:
            self.tests = list(self.defaultTest)
        elif 'ALL' in self.tests:
            expanded = []
            for t in self.tests:
                expanded.extend(all_tests if t == 'ALL' else [t])
            self.tests = expanded
        if not self.name:
            self.name = escape_commit(self.commit)

    @classmethod
    def parse_arguments(cls):
        """Build the argument parser, add testprogs-specific args, parse and return the dict."""
        parser = super().parse_arguments()
        parser.description = (
            "Compile the PHYEX package using a specific commit, "
            "run the different test progs and check if results "
            "are identical to a given version")
        parser.add_argument('--no-perf', action='store_false', dest='perf',
                            help='deactivate DR_HOOK')
        parser.add_argument('--no-check', action='store_const', dest='checkOpt',
                            const='', default='--check',
                            help='suppress value printing (comparison will be impossible)')
        parser.add_argument(
            f'--no{cls.default_buildSystem}',
            action='store_const', dest='buildSys',
            const=cls.alternative_buildSystem,
            help=f'Do not use {cls.default_buildSystem}')
        parser.add_argument('-a', dest='archfile', default=None,
                            help='architecture name to build and run the commit')
        parser.add_argument('-A', dest='refarchfile', default=None,
                            help='architecture name for the reference simulation')
        parser.add_argument('-e', dest='extrapolation', type=int, default=0,
                            help='extrapolation configuration index\n'
                             + '\n'.join(f'  {k}: {v}'
                                         for k, v in cls._extrapolation_configs()[0].items()))
        args = parser.parse_args()
        return vars(args), cls._build_commitcmd(parser)

    def submit(self, output, error, cmd_args, cwd=None):
        """Run a command, either through SLURM or directly."""
        if self.submit_method == 'slurm_belenos':
            with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
                f.write('#!/bin/bash\n')
                f.write('#SBATCH -n 1\n')
                f.write('#SBATCH -N 1\n')
                f.write('#SBATCH -t 10\n')
                f.write(f'#SBATCH --export={self.varToExport}\n')
                if subprocess.run(['ldd', cmd_args[0]],
                                  capture_output=True, check=True).stdout.count(b'libcuda') > 0:
                    f.write('#SBATCH -p ndl\n')
                f.write(f'\ncd {cwd or os.getcwd()}\n')
                f.write(' '.join(cmd_args) + '\n')
                f.flush()
                os.chmod(f.name, 0o755)

            outtmp = tempfile.mktemp()
            subprocess.run(['sbatch', '--wait', '-o', outtmp, '-e', error, f.name],
                           check=True)
            with open(outtmp, encoding='utf-8') as fh:
                content = fh.read()
            sep = '#' * 41
            idx = content.find(sep)
            if idx >= 0:
                with open(error, 'a', encoding='utf-8') as ef:
                    ef.write(content[idx-1:])
                with open(output, 'w', encoding='utf-8') as of:
                    of.write(content[:idx-2])
            else:
                with open(output, 'w', encoding='utf-8') as of:
                    of.write(content)
            os.unlink(f.name)
            os.unlink(outtmp)
        else:
            with open(output, 'w', encoding='utf-8') as out, open(error, 'w', encoding='utf-8') as err:
                subprocess.run(cmd_args, stdout=out, stderr=err, cwd=cwd, check=False)

    def fill_dirdata(self):
        """Download and extract test data if not already present."""
        os.makedirs(self.dirdata, exist_ok=True)
        files = [
            'https://github.com/UMR-CNRM/PHYEX/files/12783926/ice_adjust.tar.gz',
            'https://github.com/UMR-CNRM/PHYEX/files/12783935/rain_ice.tar.gz',
            'https://github.com/UMR-CNRM/PHYEX/files/12783942/rain_ice_old.tar.gz',
            'https://github.com/UMR-CNRM/PHYEX/files/12783945/shallow.tar.gz',
            'https://github.com/UMR-CNRM/PHYEX/files/12783952/turb.tar.gz',
            'https://github.com/user-attachments/files/17030876/lima_adjust.tar.gz',
            'https://github.com/user-attachments/files/17330177/lima.tar.gz',
        ]
        for url in files:
            basefile = os.path.basename(url)
            basename_noext = os.path.splitext(os.path.splitext(basefile)[0])[0]
            data_file = os.path.join(self.dirdata, basename_noext, '00000000.dat')
            if not os.path.isfile(data_file):
                subprocess.run(['wget', '--no-check-certificate', url, '-O', basefile],
                               cwd=self.dirdata, check=True)
                subprocess.run(['tar', 'xf', basefile], cwd=self.dirdata, check=True)
                os.remove(os.path.join(self.dirdata, basefile))

    def pack_creation(self):
        testdir = os.path.join(self.TESTDIR, self.name)
        build_dir = os.path.join(testdir, 'build', f'with_{self.buildSys}',
                                 f'arch_{self.archfile}')

        if self.suppress and os.path.isdir(testdir):
            shutil.rmtree(testdir)

        if os.path.isdir(build_dir) and self.onlyIfNeeded:
            self.packcreation = False  # To prevent pack_update from being executed
            return

        if os.path.isdir(build_dir):
            raise CheckCommitError(
                f"Directory already exists ({build_dir}),\n"
                "suppress it to be able to compile it again "
                "(or use the -s option to automatically suppress it)", 5)

        print(f"### Pack creation for commit {self.commit}")
        os.makedirs(testdir, exist_ok=True)
        build_base = os.path.join(testdir, 'build')
        if not os.path.isdir(build_base):
            toolsdir = os.path.abspath(
                os.path.join(os.path.dirname(__file__), '..', '..'))
            build_src = os.path.join(toolsdir, '..', '..', 'build')
            shutil.copytree(build_src, build_base, symlinks=True)
            arch_base = os.path.join(build_base, f'with_{self.buildSys}')
            for d in os.listdir(arch_base):
                if d.startswith('arch_'):
                    shutil.rmtree(os.path.join(arch_base, d))
        else:
            print("WARNING: the compilation system is already there, "
                  "we use it but it could be outdated")

    def pack_update(self):
        testdir = os.path.join(self.TESTDIR, self.name)
        build_dir = os.path.join(testdir, 'build', f'with_{self.buildSys}')
        if not os.path.isdir(build_dir):
            raise CheckCommitError(
                f"Compilation directory must exist ({build_dir})", 9)

        makeargs = ['./make_' + self.buildSys + '.sh']
        if self.packupdate:
            makeargs.append('-u')
        if self.packcreation:
            makeargs.append('-p')
        makeargs.extend(['--commit', self.commit, '--arch', self.archfile])
        out_path = os.path.join(build_dir, 'Output_compilation_step1')
        self._run_with_tee(makeargs, build_dir, out_path)

    def compilation_step(self):
        testdir = os.path.join(self.TESTDIR, self.name)
        build_dir = os.path.join(testdir, 'build', f'with_{self.buildSys}')
        libphyex = os.path.join(build_dir, f'arch_{self.archfile}', 'build', 'lib',
                                'libphyex.so')
        if not self.onlyIfNeeded or not os.path.isfile(libphyex):
            print(f"### Compilation of commit {self.commit}")
            makeargs = ['./make_' + self.buildSys + '.sh', '-c', '--jobs=10',
                        '--commit', self.commit, '--arch', self.archfile]
            if not self.useexpand:
                makeargs.append('--noexpand')
            out_path = os.path.join(build_dir, 'Output_compilation_step2')
            self._run_with_tee(makeargs, build_dir, out_path)

    def execution(self):
        testdir = os.path.join(self.TESTDIR, self.name)

        if not self.onlyIfNeeded:
            for t in self.tests:
                test_out = os.path.join(testdir, 'tests', f'with_{self.buildSys}',
                                        f'arch_{self.archfile}',
                                        f'{t}{self.extrapolation_tag}')
                if os.path.isdir(test_out):
                    shutil.rmtree(test_out)

        firstrun = 1
        for t in self.tests:
            if t not in self.allowedTests:
                print(f"The test {t} is not allowed")
                continue

            test_out = os.path.join(testdir, 'tests', f'with_{self.buildSys}',
                                    f'arch_{self.archfile}',
                                    f'{t}{self.extrapolation_tag}')
            if os.path.isdir(test_out) and self.onlyIfNeeded:
                continue

            if firstrun:
                print(f"### Running of commit {self.commit}")
                firstrun = 0

            prec = ''
            for p in ('', '_dp', '_sp'):
                exe = os.path.join(testdir, 'build', f'with_{self.buildSys}',
                                   f'arch_{self.archfile}', 'build', 'bin',
                                   f'main_{t}{p}.exe')
                if os.path.isfile(exe):
                    prec = p
                    break
            exe_main = os.path.join(testdir, 'build', f'with_{self.buildSys}',
                                    f'arch_{self.archfile}', 'build', 'bin',
                                    f'main_{t}{prec}.exe')
            if not os.path.isfile(exe_main):
                raise CheckCommitError(
                    f"Directory does not exist ({testdir}) or compilation has failed, "
                    "please check\n"
                    f"Run '-p -c {self.commit}' to compile.", 6)

            os.makedirs(test_out, exist_ok=True)
            self.fill_dirdata()
            os.symlink(os.path.join(self.dirdata, t), os.path.join(test_out, 'data'))

            if self.perf:
                env = os.environ.copy()
                env['DR_HOOK_OPT'] = 'prof'
                env['DR_HOOK'] = '1'
                env['DR_HOOK_IGNORE_SIGNALS'] = '-1'
            else:
                env = os.environ.copy()

            arch_env = os.path.join(testdir, 'build', f'with_{self.buildSys}',
                                    f'arch_{self.archfile}', 'arch.env')
            if os.path.isfile(arch_env):
                with open(arch_env, encoding='utf-8') as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('export '):
                            parts = line[7:].split('=', 1)
                            if len(parts) == 2:
                                env[parts[0].strip()] = os.path.expandvars(
                                    parts[1].strip())

            cmd = [exe_main]
            if self.checkOpt:
                cmd.append(self.checkOpt)
            if self.extrapolation_opts:
                cmd.extend(self.extrapolation_opts.split())

            self.submit(os.path.join(test_out, 'Output_run'),
                        os.path.join(test_out, 'Stderr_run'), cmd, cwd=test_out)
            stderr_run = os.path.join(test_out, 'Stderr_run')
            if os.path.isfile(stderr_run):
                with open(stderr_run, encoding='utf-8') as f:
                    err_content = f.read()
                if err_content.strip():
                    print(err_content, file=sys.stderr)

            if self.perf:
                prof0 = os.path.join(test_out, 'drhook.prof.0')
                if os.path.isfile(prof0):
                    with open(prof0, encoding='utf-8') as f:
                        content = f.read()
                    first_line = None
                    for i, line in enumerate(content.split('\n')):
                        if re.match(r'^\s*1', line):
                            first_line = i + 1
                            break
                    if first_line is not None:
                        subprocess.run([sys.executable, '-c', """
import numpy, pandas, re
d = {'time': ('<f4', ('mean',)), 'self': ('<f4', ('mean', 'max', 'min', 'std', 'sum')),
     'total': ('<f4', ('mean', 'max', 'min', 'std', 'sum')), 'calls': ('<i4', ('sum',)),
     'self_per_call': ('<f4', ('mean',)), 'total_per_call': ('<f4', ('mean',)), 'routine': ('U256', '')}
first_line = """ + str(first_line) + """
arraynp = numpy.loadtxt('drhook.prof.0',
    dtype=[(k, v[0]) for (k, v) in d.items()],
    converters={8: lambda s: s.split(b'@')[0].lstrip(b'*')},
    skiprows=first_line - 1, usecols=[1, 3, 4, 5, 6, 7, 8], encoding='bytes')
df = pandas.DataFrame(arraynp).groupby('routine').agg(
    **{k + '_' + agg: pandas.NamedAgg(column=k, aggfunc=agg)
       for k in d.keys() for agg in d[k][1]
       if k != 'routine'}).sort_values('self_sum', ascending=False)
df.index.name += ' ordered by self_sum'
with open('drhook.prof.agg', 'w', encoding='utf-8') as f:
    f.write(df.to_string())
"""], env=env, cwd=test_out, check=True)

    def performance_evaluation(self):
        testdir = os.path.join(self.TESTDIR, self.name)
        if not self.run_tests or not self.perffile:
            return

        print("### Evaluate performance for commit", self.commit)

        ZTD_sum = 0.0
        ZTC_sum = 0.0
        firstrun = 1
        for t in self.tests:
            if t not in self.allowedTests:
                continue

            exe_main = os.path.join(testdir, 'build', f'with_{self.buildSys}',
                                    f'arch_{self.archfile}', 'build', 'bin',
                                    f'main_{t}.exe')
            if not os.path.isfile(exe_main):
                raise CheckCommitError(
                    f"Directory does not exist ({testdir}) or compilation has failed\n"
                    f"Run '-p -c {self.commit}' to compile.", 7)

            if firstrun:
                firstrun = 0
                arch_env = os.path.join(testdir, 'build', f'with_{self.buildSys}',
                                        f'arch_{self.archfile}', 'arch.env')
                env = os.environ.copy()
                if os.path.isfile(arch_env):
                    with open(arch_env, encoding='utf-8') as f:
                        for line in f:
                            line = line.strip()
                            if line.startswith('export '):
                                parts = line[7:].split('=', 1)
                                if len(parts) == 2:
                                    env[parts[0].strip()] = os.path.expandvars(
                                        parts[1].strip())

                NPOINTS = 100000
                NPOINTS_perf = int(env.get('NPOINTS_perf', str(NPOINTS)))
                if NPOINTS_perf < NPOINTS:
                    NTIMES = round(NPOINTS / NPOINTS_perf)
                    NPOINTS = NPOINTS_perf
                else:
                    NTIMES = 1
                NPROMA = int(env.get('NPROMA_perf', '32'))
                OMP_NUM_THREADS = int(env.get('OMP_NUM_THREADS_perf', '8'))
                NBLOCKS = NPOINTS // NPROMA // 8 * 8

                perf_tag_template = self.conf_extra_tag[4]
                perf_tag = (perf_tag_template
                            .replace('${NPROMA}', str(NPROMA))
                            .replace('${NBLOCKS}', str(NBLOCKS))
                            .replace('${NTIMES}', str(NTIMES)))

                if not self.onlyIfNeeded:
                    for t2 in self.tests:
                        perf_dir = os.path.join(
                            testdir, 'tests', f'with_{self.buildSys}',
                            f'arch_{self.archfile}', f'{t2}{perf_tag}')
                        if os.path.isdir(perf_dir):
                            shutil.rmtree(perf_dir)

            perf_tag_template = self.conf_extra_tag[4]
            perf_tag = (perf_tag_template
                        .replace('${NPROMA}', str(NPROMA))
                        .replace('${NBLOCKS}', str(NBLOCKS))
                        .replace('${NTIMES}', str(NTIMES)))

            env_perf = os.environ.copy()
            env_perf['NPROMA'] = str(NPROMA)
            env_perf['NBLOCKS'] = str(NBLOCKS)
            env_perf['NTIMES'] = str(NTIMES)
            env_perf['OMP_NUM_THREADS'] = str(OMP_NUM_THREADS)

            subprocess.run(
                [sys.executable, '-m', 'pyphyextools.testing.check_testprogs',
                 '-r', '-t', t, '-a', self.archfile, '--no-check', '--no-perf',
                 '-e', '4', '--name', self.name, self.commit],
                env=env_perf, check=True)

            output_file = os.path.join(
                testdir, 'tests', f'with_{self.buildSys}',
                f'arch_{self.archfile}', f'{t}{perf_tag}', 'Output_run')
            if os.path.isfile(output_file):
                with open(output_file, encoding='utf-8') as f:
                    content = f.read()

                ZTD = None
                for line in content.split('\n'):
                    m = re.match(r'.*\bZTD\s*=\s*(\S+)', line)
                    if m:
                        ZTD = float(m.group(1))
                        break
                if ZTD is not None:
                    ZTD_sum = ZTD_sum + ZTD if ZTD_sum >= 0 else ZTD
                else:
                    ZTD = -999.0
                    ZTD_sum = -999.0

                ZTC = None
                for line in content.split('\n'):
                    m = re.match(r'.*\bZTC\s*=\s*(\S+)', line)
                    if m:
                        ZTC = float(m.group(1))
                        break
                if ZTC is not None:
                    ZTC_sum = ZTC_sum + ZTC if ZTC_sum >= 0 else ZTC
                else:
                    ZTC = -999.0
                    ZTC_sum = -999.0
            else:
                ZTD = -999.0
                ZTD_sum = -999.0
                ZTC = -999.0
                ZTC_sum = -999.0

            with open(self.perffile, 'a', encoding='utf-8') as pf:
                pf.write(f"{self.commit} testprogs {t} {ZTD} {ZTC}\n")

        with open(self.perffile, 'a', encoding='utf-8') as pf:
            pf.write(f"{self.commit} testprogs ALL {ZTD_sum} {ZTC_sum}\n")

    def comparison(self):
        testdir = os.path.join(self.TESTDIR, self.name)
        print(f"### Check commit {self.commit} against commit {self.reference}")
        refByTest = self.resolve_ref_per_test(self.json_content)

        # First part: build test info dict, collect missing refs by caseref
        test_info = {}
        missing_by_ref = {}
        for t in self.tests:
            if t not in self.allowedTests:
                continue

            caseref = refByTest.get(t, '')
            if '/' in caseref:
                refname = escape_commit(self.reference)
            else:
                refname = caseref

            file1 = os.path.join(testdir, 'tests', f'with_{self.buildSys}',
                                 f'arch_{self.archfile}',
                                 f'{t}{self.extrapolation_tag}', 'Output_run')
            file2 = os.path.join(self.TESTDIR, refname, 'tests',
                                 f'with_{self.buildSys}',
                                 f'arch_{self.refarchfile}',
                                 f'{t}{self.extrapolation_tag}', 'Output_run')

            test_info[t] = dict(refname=refname, file1=file1, file2=file2,
                                caseref=caseref)
            if not os.path.isfile(file2):
                missing_by_ref.setdefault(caseref, []).append(t)

        # Between parts: compute missing references in batch (one per caseref)
        if self.computeRefIfNeeded:
            build_sys_arg = ''
            if self.default_buildSystem != self.buildSys:
                build_sys_arg = f'--no{self.default_buildSystem}'
            for caseref, tests_list in missing_by_ref.items():
                # We cannot use directly a python instance because clone_and_run
                # use arguments received on sys.argv to build the command to execute
                subprocess.run(
                    [sys.executable, '-m', 'pyphyextools.testing.check_testprogs',
                     '-p', '-c', '-r', '-t', ','.join(tests_list),
                     '-a', self.refarchfile, '--onlyIfNeeded',
                     '-e', str(self.extrapolation), '--no-perf', caseref] +
                    ([build_sys_arg] if build_sys_arg else []),
                    check=True)

        # Second part: compare each test
        alltests = 0
        message = ""
        for t, info in test_info.items():
            file1 = info['file1']
            file2 = info['file2']
            mess = ""
            te = 0
            if not os.path.isfile(file1):
                mess = f"Result ({file1}) for commit {self.commit} does not exist, " \
                       "please run the simulation"
                te = 1
            if not os.path.isfile(file2):
                mess2 = f"Result ({file2}) for commit {info['caseref']} does not exist, " \
                        "please run the simulation"
                te = 1
                mess = mess2 if not mess else f"{mess} and {mess2}"

            if te == 0:
                te = comp_testprogs(file1, file2)
                mess = ""

            if te != 0:
                message += f" {mess}\n"
            alltests += te

        if alltests == 0:
            print("SUCCESS, files are identical")
            return 0
        else:
            print("*************** Files are different *******************")
            print(message)
            return 50

    def cleaning(self):
        testdir = os.path.join(self.TESTDIR, self.name)
        print(f"### Remove model directory for commit {self.commit}")
        if os.path.isdir(testdir):
            shutil.rmtree(testdir)


def main():
    """Entry point for the ``phyex-check_testprogs`` command."""
    run_tool(CheckCommitTestprogs, module_name='pyphyextools.testing.check_testprogs')


if __name__ == '__main__':
    main()
