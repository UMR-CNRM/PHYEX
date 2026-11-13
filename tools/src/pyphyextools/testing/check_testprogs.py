"""Commit-check tool for the PHYEX offline testprogs.

Compiles the PHYEX package using a specific commit, runs the different
test programs and checks whether results are identical to a given version.
"""

import json
import os
import re
import shlex
import shutil
import string
import sys
import tempfile
from pathlib import Path

from pyphyextools.testing.check_common import (
    CheckCommitBase, CheckCommitError, aggregate_drhook_profile, escape_commit, PHYEXCONF_DEFAULT,
    run_tool, SEPARATOR, sourced_env)
from pyphyextools import download_and_extract, run_command
from pyphyextools.testing.compare import comp_testprogs


class CheckCommitTestprogs(CheckCommitBase):
    """Check a commit against the offline testprogs reference."""

    default_expand = False
    default_shuman_trans = False
    default_precision = 'dp'

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
                 perf=True, checkOpt='--check', precision=None, **kwargs):
        super().__init__(**kwargs)
        self.defaultTest = ["ice_adjust", "rain_ice", "rain_ice_old",
                            "turb", "shallow", "lima_adjust", "lima"]
        self.allowedTests = ["ice_adjust", "rain_ice", "rain_ice_old",
                             "turb", "shallow", "lima_adjust", "lima"]
        self.perf = perf
        self.checkOpt = checkOpt
        self.conf_extra_tag, self.conf_extra_opts = self._extrapolation_configs()
        self.archfile = archfile
        self.refarchfile = refarchfile
        self.extrapolation = extrapolation
        self.precision = precision or self.default_precision
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
        if opts:
            self.extrapolation_opts = os.path.expandvars(opts)
            self.partition = 'normal256'  # job can be memory-intensive
        else:
            self.extrapolation_opts = ''
            self.partition = 'shared'  # job is small

        if os.path.isdir(os.path.join(self.commit, 'src')):
            self.model_ready = False
            json_path = os.path.join(self.commit, 'src', 'phyex', 'testprogs_version.json')
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
        parser.add_argument('-a', dest='archfile', default=None,
                            help='architecture name to build and run the commit')
        parser.add_argument('-A', dest='refarchfile', default=None,
                            help='architecture name for the reference simulation')
        parser.add_argument('-e', dest='extrapolation', type=int, default=0,
                            help='extrapolation configuration index\n'
                             + '\n'.join(f'  {k}: {v}'
                                         for k, v in cls._extrapolation_configs()[0].items()))
        parser.add_argument('--prec', dest='precision', type=str, default=cls.default_precision,
                            help="precision choice")
        args = parser.parse_args()
        return vars(args), cls._build_commitcmd(parser)

    def _arch_env(self):
        """Environment file of the architecture (copied in the build directory)"""
        return os.path.join(self.TESTDIR, self.name, 'build', f'arch_{self.archfile}',
                            'arch.env')

    def submit(self, output, error, cmd_args, cwd=None, env=None):
        """Run a command, either through SLURM or directly. With SLURM, only
        some variables are exported to the job (varToExport), the environment
        file is then sourced in the job script."""
        if self.submit_method == 'slurm_belenos':
            with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
                f.write('#!/bin/bash\n')
                f.write('#SBATCH -n 1\n')
                f.write('#SBATCH -N 1\n')
                f.write('#SBATCH -t 10\n')
                f.write(f'#SBATCH -p {self.partition}\n')
                f.write(f'#SBATCH --export={self.varToExport}\n')
                if run_command(['ldd', cmd_args[0]]).stdout.count('libcuda') > 0:
                    f.write('#SBATCH -p ndl\n')
                f.write(f'\n. {shlex.quote(self._arch_env())}\n')
                f.write(f'cd {cwd or os.getcwd()}\n')
                f.write(' '.join(cmd_args) + '\n')
                f.flush()
                os.chmod(f.name, 0o755)

            fd, outtmp = tempfile.mkstemp()
            os.close(fd)
            run_command(['sbatch', '--wait', '-o', outtmp, '-e', error, f.name], env=env)
            content = Path(outtmp).read_text(encoding='utf-8')
            # The SLURM report (after a line of '#') goes in the error file
            idx = content.find('#' * 41)
            if idx >= 0:
                with open(error, 'a', encoding='utf-8') as ef:
                    ef.write(content[idx-1:])
                Path(output).write_text(content[:idx-2], encoding='utf-8')
            else:
                Path(output).write_text(content, encoding='utf-8')
            os.unlink(f.name)
            os.unlink(outtmp)
        else:
            run_command(cmd_args, cwd=cwd, out_path=output, err_path=error, env=env,
                        check=False, display='never')

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
                download_and_extract(url, self.dirdata)

    def pack_creation(self):
        testdir = os.path.join(self.TESTDIR, self.name)
        build_dir = os.path.join(testdir, 'build', f'arch_{self.archfile}')

        if self.suppress and os.path.isdir(testdir):
            shutil.rmtree(testdir)

        if os.path.isdir(build_dir) and self.onlyIfNeeded:
            if any(os.path.isfile(lib) for lib in self._libphyex()):
                self.packcreation = False  # To prevent pack_update from being executed
                return
            # The library has not been built: the directory has been left in an
            # unknown state by a previous attempt (e.g. failed CMake configuration),
            # it is created again
            print(f"WARNING: the compilation directory ({build_dir}) is incomplete, "
                  "it is suppressed to be created again")
            shutil.rmtree(build_dir)

        if os.path.isdir(build_dir):
            raise CheckCommitError(
                f"Directory already exists ({build_dir}),\n"
                "suppress it to be able to compile it again "
                "(or use the -s option to automatically suppress it)", 5)

        if os.path.isdir(os.path.join(testdir, 'build')):
            print("WARNING: the compilation system is already there, "
                  "we use it but it could be outdated")
        # Directory needed by pack_update
        os.makedirs(os.path.join(testdir, 'build'), exist_ok=True)

    def pack_update(self):
        testdir = os.path.join(self.TESTDIR, self.name)
        build_dir = os.path.join(testdir, 'build')
        if not os.path.isdir(build_dir):
            raise CheckCommitError(
                f"Compilation directory must exist ({build_dir})", 9)

        # Source code, in the PHYEX subdirectory of testdir
        phyex_dir = os.path.join(testdir, 'PHYEX')
        if self.commit.startswith('testprogs' + SEPARATOR) or \
           (os.path.isdir(self.commit) and self.model_ready):
            raise NotImplementedError(
                "Commits ready for testprogs (testprogs branches) cannot be compiled "
                "with the CMake build system")
        if os.path.isdir(phyex_dir):
            shutil.rmtree(phyex_dir)
        if os.path.isdir(self.commit):
            print(f"Copy {self.commit}")
            self._copy_source(self.commit, phyex_dir)
        else:
            # Commit hash or tag
            if self.repo_protocol == 'https':
                repo_url = f'https://github.com/{self.repo_user}/PHYEX.git'
            else:
                repo_url = f'git@github.com:{self.repo_user}/PHYEX.git'
            print(f"Clone repository, and checkout commit {self.commit}")
            run_command(['git', 'clone', repo_url, phyex_dir])
            run_command(['git', 'checkout', self.commit], cwd=phyex_dir)

        # Architecture files, taken in the arch subdirectory of the source code.
        # They are copied in the architecture specific build directory (arch.env,
        # used by the next steps, and arch.ecbuild, used as CMake toolchain file)
        arch_dir = os.path.join(phyex_dir, 'arch')
        envfile = os.path.join(arch_dir, f'arch-{self.archfile}.env')
        toolchain = os.path.join(arch_dir, f'arch-{self.archfile}.ecbuild')
        for f in (envfile, toolchain):
            if not os.path.isfile(f):
                raise CheckCommitError(f"Architecture file not found: {f}", 6)
        arch_build_dir = os.path.join(build_dir, f'arch_{self.archfile}')
        os.makedirs(arch_build_dir, exist_ok=True)
        arch_env = os.path.join(arch_build_dir, 'arch.env')
        shutil.copy(envfile, arch_env)
        arch_toolchain = os.path.join(arch_build_dir, 'arch.ecbuild')
        shutil.copy(toolchain, arch_toolchain)

        # CMake configuration only (the compilation is done by the next step),
        # executed in the environment defined by the environment file
        if self.useexpand and self.useshumantrans:
            transformation = 'IAL'
        elif self.useexpand:
            transformation = 'IAL_NOSHUGRAD'
        elif self.useshumantrans:
            transformation = 'IAL_NOEXPAND'
        else:
            transformation = 'NONE'
        cmake_args = ['cmake', '-S', phyex_dir, '-B', arch_build_dir,
                      '-DENABLE_PHYEX_BUILD_PROGS=ON',
                      '-DENABLE_SINGLE_PRECISION=ON',
                      '-DCMAKE_BUILD_TYPE=Release',
                      '-DBUILD_SHARED_LIBS=OFF',
                      f'-DPHYEX_TRANSFORMATION={transformation}',
                      f'-DCMAKE_TOOLCHAIN_FILE={arch_toolchain}',
                      '-DENABLE_DR_HOOK_NVTX=OFF',
                      '-DCMAKE_CXX_FLAGS=-lstdc++']
        out_path = os.path.join(arch_build_dir, 'Output_compilation_step1')
        run_command(cmake_args, build_dir, out_path, env=sourced_env(arch_env),
                    display='always')

    @staticmethod
    def _copy_source(srcdir, destdir):
        """Copy the source code of a PHYEX directory: in a git repository,
        tracked files and untracked files not ignored (local modifications are
        kept, build directories and virtual environments are not copied),
        otherwise the whole directory (except .git)."""
        if os.path.exists(os.path.join(srcdir, '.git')):
            result = run_command(
                ['git', 'ls-files', '-z', '--cached', '--others', '--exclude-standard'],
                cwd=srcdir, banner=False)
            for rel in sorted(set(result.stdout.split('\0'))):
                src = os.path.join(srcdir, rel)
                if rel == '' or not os.path.lexists(src):
                    continue  # deleted but still tracked file
                dst = os.path.join(destdir, rel)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy2(src, dst, follow_symlinks=False)
        else:
            shutil.copytree(srcdir, destdir, symlinks=True,
                            ignore=shutil.ignore_patterns('.git'))

    def _libphyex(self):
        """Possible paths of the library built by CMake (static with
        BUILD_SHARED_LIBS=OFF, see pack_update)"""
        arch_build_dir = os.path.join(self.TESTDIR, self.name, 'build', f'arch_{self.archfile}')
        return [os.path.join(arch_build_dir, 'lib', f'libphyex_{self.precision}{ext}')
                for ext in ('.a', '.so')]

    def compilation_step(self):
        testdir = os.path.join(self.TESTDIR, self.name)
        build_dir = os.path.join(testdir, 'build')
        arch_build_dir = os.path.join(build_dir, f'arch_{self.archfile}')
        if not self.onlyIfNeeded or not any(os.path.isfile(lib) for lib in self._libphyex()):
            print(f"### Compilation of commit {self.commit}")
            # Build (the CMake configuration is done by pack_update), executed
            # in the environment defined by the environment file
            cmake_args = ['cmake', '--build', arch_build_dir, '--parallel', '10']
            out_path = os.path.join(arch_build_dir, 'Output_compilation_step2')
            run_command(cmake_args, build_dir, out_path, env=sourced_env(self._arch_env()),
                        display='always')

    def execution(self):
        testdir = os.path.join(self.TESTDIR, self.name)

        if not self.onlyIfNeeded:
            for t in self.tests:
                test_out = os.path.join(testdir, 'tests',
                                        f'arch_{self.archfile}',
                                        f'{t}{self.extrapolation_tag}')
                if os.path.isdir(test_out):
                    shutil.rmtree(test_out)

        firstrun = True
        for t in self.tests:
            if t not in self.allowedTests:
                print(f"The test {t} is not allowed")
                continue

            test_out = os.path.join(testdir, 'tests',
                                    f'arch_{self.archfile}',
                                    f'{t}{self.extrapolation_tag}')
            if os.path.isdir(test_out) and self.onlyIfNeeded:
                continue

            if firstrun:
                print(f"### Running of commit {self.commit}")
                firstrun = False

            exe_main = os.path.join(testdir, 'build',
                                    f'arch_{self.archfile}', 'bin',
                                    f'main_{t}_{self.precision}.exe')
            if not os.path.isfile(exe_main):
                raise CheckCommitError(
                    f"Directory does not exist ({testdir}) or compilation has failed, "
                    "please check\n"
                    f"Run '-p -c {self.commit}' to compile.", 6)

            os.makedirs(test_out, exist_ok=True)
            self.fill_dirdata()
            os.symlink(os.path.join(self.dirdata, t), os.path.join(test_out, 'data'))

            env = sourced_env(self._arch_env())
            if self.perf:
                env.update({'DR_HOOK_OPT': 'prof', 'DR_HOOK': '1',
                            'DR_HOOK_IGNORE_SIGNALS': '-1'})

            cmd = [exe_main]
            if self.checkOpt:
                cmd.append(self.checkOpt)
            if self.extrapolation_opts:
                cmd.extend(self.extrapolation_opts.split())

            self.submit(os.path.join(test_out, 'Output_run'),
                        os.path.join(test_out, 'Stderr_run'), cmd,
                        cwd=test_out, env=env)
            stderr_run = Path(test_out, 'Stderr_run')
            if stderr_run.is_file():
                err_content = stderr_run.read_text(encoding='utf-8')
                if err_content.strip():
                    print(err_content, file=sys.stderr)

            if self.perf:
                aggregate_drhook_profile(os.path.join(test_out, 'drhook.prof.1'),
                                         os.path.join(test_out, 'drhook.prof.agg'))

    @staticmethod
    def _add_timing(total, value):
        """Sum of timings: -999 as soon as a value is not available (None)"""
        if value is None or total < 0:
            return -999.0
        return total + value

    def performance_evaluation(self):
        testdir = os.path.join(self.TESTDIR, self.name)
        if not self.run_tests or not self.perffile:
            return

        print("### Evaluate performance for commit", self.commit)

        ZTD_sum = 0.0
        ZTC_sum = 0.0
        firstrun = True
        for t in self.tests:
            if t not in self.allowedTests:
                continue

            exe_main = os.path.join(testdir, 'build',
                                    f'arch_{self.archfile}', 'bin',
                                    f'main_{t}_{self.precision}.exe')
            if not os.path.isfile(exe_main):
                raise CheckCommitError(
                    f"Directory does not exist ({testdir}) or compilation has failed\n"
                    f"Run '-p -c {self.commit}' to compile.", 7)

            if firstrun:
                firstrun = False
                env = sourced_env(self._arch_env())

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

                perf_tag = string.Template(self.conf_extra_tag[4]).substitute(
                    NPROMA=NPROMA, NBLOCKS=NBLOCKS, NTIMES=NTIMES)

                if not self.onlyIfNeeded:
                    for t2 in self.tests:
                        perf_dir = os.path.join(
                            testdir, 'tests',
                            f'arch_{self.archfile}', f'{t2}{perf_tag}')
                        if os.path.isdir(perf_dir):
                            shutil.rmtree(perf_dir)

            env_perf = os.environ.copy()
            env_perf['NPROMA'] = str(NPROMA)
            env_perf['NBLOCKS'] = str(NBLOCKS)
            env_perf['NTIMES'] = str(NTIMES)
            env_perf['OMP_NUM_THREADS'] = str(OMP_NUM_THREADS)

            run_command(
                [sys.executable, '-m', 'pyphyextools.testing.check_testprogs',
                 '-r', '-t', t, '-a', self.archfile, '--no-check', '--no-perf',
                 '-e', '4', '--name', self.name, '--prec', self.precision, self.commit],
                env=env_perf)

            # Timings (ZTD and ZTC values printed by the program), -999 if not found
            output_file = Path(testdir, 'tests', f'arch_{self.archfile}', f'{t}{perf_tag}',
                               'Output_run')
            content = output_file.read_text(encoding='utf-8') if output_file.is_file() else ''
            timings = {}
            for name in ('ZTD', 'ZTC'):
                m = re.search(rf'\b{name}[ \t]*=[ \t]*(\S+)', content)
                timings[name] = float(m.group(1)) if m else None
            ZTD, ZTC = (-999.0 if timings[name] is None else timings[name]
                        for name in ('ZTD', 'ZTC'))
            ZTD_sum = self._add_timing(ZTD_sum, timings['ZTD'])
            ZTC_sum = self._add_timing(ZTC_sum, timings['ZTC'])

            with open(self.perffile, 'a', encoding='utf-8') as pf:
                pf.write(f"{self.name} testprogs {t} {ZTD} {ZTC}\n")

        with open(self.perffile, 'a', encoding='utf-8') as pf:
            pf.write(f"{self.name} testprogs ALL {ZTD_sum} {ZTC_sum}\n")

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

            file1 = os.path.join(testdir, 'tests',
                                 f'arch_{self.archfile}',
                                 f'{t}{self.extrapolation_tag}', 'Output_run')
            file2 = os.path.join(self.TESTDIR, refname, 'tests',
                                 f'arch_{self.refarchfile}',
                                 f'{t}{self.extrapolation_tag}', 'Output_run')

            test_info[t] = dict(refname=refname, file1=file1, file2=file2,
                                caseref=caseref)
            if not os.path.isfile(file2):
                missing_by_ref.setdefault(caseref, []).append(t)

        # Between parts: compute missing references in batch (one per caseref)
        if self.computeRefIfNeeded:
            for caseref, tests_list in missing_by_ref.items():
                # We cannot use directly a python instance because clone_and_run
                # use arguments received on sys.argv to build the command to execute
                run_command(
                    [sys.executable, '-m', 'pyphyextools.testing.check_testprogs',
                     '-p', '-c', '-r', '-t', ','.join(tests_list),
                     '-a', self.refarchfile, '--onlyIfNeeded',
                     '-e', str(self.extrapolation), '--no-perf',
                     '--prec', self.precision, caseref])

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
