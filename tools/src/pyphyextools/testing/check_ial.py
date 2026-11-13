"""Commit-check tool for the IAL / AROME model.

Compiles the AROME model using a specific commit for the externalised
physics, runs a small 3D case and checks whether results are identical
to a given reference.
"""

import glob
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

from pyphyextools.testing.check_common import (
    CheckCommitBase, CheckCommitError, aggregate_drhook_profile, escape_commit, mvdiff,
    run_tool, timed)
from pyphyextools import REPO_DIR, run_command
from pyphyextools.prep_code import prep_code
from pyphyextools.testing.compare import comp_binary, comp_NODE


class CheckCommitIAL(CheckCommitBase):
    """Check a commit against the IAL / AROME reference."""

    default_expand = True
    default_shuman_trans = True

    def __init__(self, prepCodeOpts='', fullcompilation=False, **kwargs):
        super().__init__(**kwargs)
        self.defaultTest = ["small_3D"]
        self.allowedTests = ["small_3D", "small_3D_np2", "small_3D_alt1",
                             "small_3D_alt2", "small_3D_alt3", "small_3D_alt4",
                             "small_3D_alt5", "small_3D_alt6", "small_3D_alt7",
                             "small_3D_alt8", "small_3D_alt9", "small_3D_alt10",
                             "small_3D_alt11", "small_3D_lima",
                             "small_3D_xfrmin", "arp_t31"]
        self.prepCodeOpts = prepCodeOpts
        self.fullcompilation = fullcompilation
        self.HOMEPACK = os.environ.get('HOMEPACK', os.path.expanduser('~/pack'))
        self.HPC = 0
        self.gmkpack_l = 'BOMPI5GFORT141DP50T2'
        self.gmkpack_o = 'dp'
        self.ialdir = None
        self.cycle = None
        self.phytoolsdir = os.path.join(REPO_DIR, 'tools')

        hostname = os.uname().nodename
        if hostname[:7] in ('belenos', 'taranis'):
            self.HPC = 1
            self.gmkpack_l = 'BIMPIIFCI2302SP50T2'
            self.gmkpack_o = 'sp'
            self.allowedTests.append('big_3D')

        if os.path.isdir(os.path.join(self.commit, 'src')):
            self.model_ready = False
            json_path = os.path.join(self.commit, 'src', 'phyex', 'ial_version.json')
        else:
            self.model_ready = True
            json_path = os.path.join(self.commit, 'ial_version.json')
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

        self.cycle = self.json_content.get('cycle', '')
        self.ialdir = f"PHYEX/{self.cycle}_{self.name}.01.{self.gmkpack_l}.{self.gmkpack_o}"

    @classmethod
    def parse_arguments(cls):
        """Build the argument parser, add IAL-specific args, parse and return the dict."""
        parser = super().parse_arguments()
        extra_doc = ("The PHYEXROOTPACK environment variable, if set, is used as the argument\n"
                     "of the --rootpack option of ial-git2pack/ial-to_pack, for incremental packs.")
        parser.description = (
            "Compile the AROME model using a specific commit for the externalised physics, "
            "run a small 3D case and check if results are identical to a given version")
        parser.epilog = extra_doc
        parser.add_argument('--prep_code-opts', dest='prepCodeOpts', default='',
                            help='options added to prep_code call')
        parser.add_argument('-f', action='store_true', dest='fullcompilation',
                            help='full compilation (no pre-compiled pack)')
        args = parser.parse_args()
        return vars(args), cls._build_commitcmd(parser)

    def submit(self, output, script_args, cwd=None, env=None):
        """Run a script, either through SLURM or directly."""
        if self.HPC:
            run_command(['sbatch', '--wait', '-o', output, ' '.join(script_args)],
                        cwd=cwd, env=env)
            with open(output, encoding='utf-8') as f:
                print(f.read(), end='')
        else:
            run_command(script_args, cwd, output, env=env, display='always')

    def pack_creation(self):
        """Create the IAL pack via *ial-to_pack*."""
        ialdir_full = os.path.join(self.HOMEPACK, self.ialdir)

        if self.suppress and os.path.isdir(ialdir_full):
            shutil.rmtree(ialdir_full)

        if os.path.isdir(ialdir_full) and self.onlyIfNeeded:
            self.packcreation = False  # To prevent pack_update from being executed
            return

        if os.path.isdir(ialdir_full):
            raise CheckCommitError(
                f"Pack already exists ({ialdir_full}),\n"
                "suppress it to be able to compile it again "
                "(or use the -s option to automatically suppress it)", 5)

        print(f"### Pack creation for commit {self.commit}")
        os.environ['GMKTMP'] = '/dev/shm'

        if shutil.which('ial-to_pack') is None:
            raise CheckCommitError("ial-to_pack not found, please install it", 7)

        tmpbuilddir = tempfile.mkdtemp()
        try:
            repo_url = self.json_content.get('IALrepo', 'git@github.com:ACCORD-NWP/IAL.git')
            run_command(['git', 'clone', repo_url], cwd=tmpbuilddir)
            ial_repo = os.path.join(tmpbuilddir, 'IAL')
            ial_commit = self.json_content.get('IALcommit', 'XXXXXXX')
            run_command(['git', 'checkout', ial_commit], cwd=ial_repo)

            IALbundle_tag = self.json_content.get('IALbundle_tag', '')
            bundle_args = ['--hub_bundle_tag', IALbundle_tag] if IALbundle_tag else []

            if self.fullcompilation:
                kind = 'main'
                rootpack_opt = []
            else:
                kind = 'incr'
                phyexrootpack = os.environ.get('PHYEXROOTPACK', '')
                rootpack_opt = ['--rootpack', phyexrootpack] if phyexrootpack else []

            pack_dir = os.path.join(tmpbuilddir, 'pack')
            cmd = ['ial-to_pack', '-l', self.gmkpack_l, '-o', self.gmkpack_o,
                   '-t', kind, '-n', '10', '-r', ial_repo, '-p', 'masterodb',
                   '--homepack', pack_dir] + rootpack_opt + bundle_args
            env = os.environ.copy()
            env['LANG'] = 'C'
            run_command(cmd, cwd=tmpbuilddir, input_string='y\n', env=env)

            oldname = os.listdir(pack_dir)[0]
            old_path = os.path.join(pack_dir, oldname)
            shutil.move(old_path, ialdir_full)

            for fname in os.listdir(ialdir_full):
                if fname.startswith('ics_'):
                    fpath = Path(ialdir_full, fname)
                    content = fpath.read_text(encoding='utf-8')
                    fpath.write_text(content.replace(old_path, ialdir_full), encoding='utf-8')

            phyex_dir = os.path.join(ialdir_full, 'hub', 'local', 'src', 'PHYEX', 'phyex')
            shutil.rmtree(phyex_dir, ignore_errors=True)
            os.makedirs(phyex_dir)
        finally:
            shutil.rmtree(tmpbuilddir, ignore_errors=True)

    def pack_update(self):
        """Copy and prepare source files with *prep_code*."""
        ialdir_full = os.path.join(self.HOMEPACK, self.ialdir)
        phyex_dir = os.path.join(ialdir_full, 'hub', 'local', 'src', 'PHYEX', 'phyex')

        print(f"Copy {self.commit}")
        os.makedirs(os.path.join(phyex_dir, 'PHYEX'), exist_ok=True)

        if not self.model_ready:
            run_command(['scp', '-q', '-r', os.path.join(self.commit, 'src'), 'PHYEX/'],
                        cwd=phyex_dir)
            prep_kwargs = self._parse_prep_code_opts(self.prepCodeOpts)
            # We do not use the mnh_expand option of prep_code because expansion
            # must be done after all other transformations
            pyfortool_options = ['--removeACC']
            if self.useshumantrans:
                pyfortool_options.append('--shugradTransfo')
            if self.useexpand:
                # --mnhExpand after --shugradTransfo to deal with the directives
                # added during the transformation
                pyfortool_options.append('--mnhExpand')
            if self.useshumantrans or self.useexpand:
                pyfortool_options.append('--splitLongLines')
            prep_code(
                directory=os.path.join(phyex_dir, 'PHYEX'),
                model='arome',
                subs=['gmkpack_ignored_files', 'turb', 'micro', 'aux', 'conv',
                      'CMakeLists.txt', 'VERSION', 'cmake'],
                pyfortool_options=pyfortool_options,
                **prep_kwargs)
        else:
            print("model ready")
            # '*' expanded as the shell would do (the command is not run by a shell)
            run_command(['scp', '-q', '-r'] + sorted(glob.glob(os.path.join(self.commit, '*'))) +
                        ['PHYEX/'], cwd=phyex_dir)
            prep_kwargs = self._parse_prep_code_opts(self.prepCodeOpts)
            prep_code(directory=os.path.join(phyex_dir, 'PHYEX'), **prep_kwargs)

        for root, dirs, files in os.walk(os.path.join(phyex_dir, 'PHYEX')):
            for fname in files:
                os.utime(os.path.join(root, fname), None)

        if self.packupdate:
            find_root = os.path.join(phyex_dir, 'PHYEX')
            for subdir in ('turb', 'micro', 'conv', 'aux', 'cmake'):
                sub_path = os.path.join(find_root, subdir)
                if not os.path.isdir(sub_path):
                    continue
                for root, dirs, files in os.walk(sub_path):
                    for fname in files:
                        rel = os.path.relpath(os.path.join(root, fname), find_root)
                        mvdiff(os.path.join(find_root, rel), os.path.join(phyex_dir, rel))
            for cmake_file in ('CMakeLists.txt', 'VERSION'):
                if os.path.isfile(os.path.join(phyex_dir, 'PHYEX', cmake_file)):
                    mvdiff(os.path.join(phyex_dir, 'PHYEX', cmake_file),
                           os.path.join(phyex_dir, cmake_file))
            shutil.rmtree(os.path.join(phyex_dir, 'PHYEX'))
        else:
            for rep in ('turb', 'micro', 'conv', 'aux', 'cmake'):
                src = os.path.join(phyex_dir, 'PHYEX', rep)
                if os.path.isdir(src):
                    shutil.move(src, phyex_dir)
            for cmake_file in ('CMakeLists.txt', 'VERSION'):
                cmake_src = os.path.join(phyex_dir, 'PHYEX', cmake_file)
                if os.path.isfile(cmake_src):
                    shutil.move(cmake_src, phyex_dir)

        ignored_file = os.path.join(phyex_dir, 'PHYEX', 'gmkpack_ignored_files')
        if os.path.isfile(ignored_file):
            if not self.fullcompilation:
                if not self.packupdate:
                    with open(ignored_file, encoding='utf-8') as fh:
                        lines = [l.strip() for l in fh if l.strip()]
                    if lines:
                        ics_file = Path(ialdir_full, 'ics_masterodb')
                        insert = '\n'.join(lines)
                        ics_file.write_text(ics_file.read_text(encoding='utf-8').replace(
                            'end_of_ignored_files', f'end_of_ignored_files\n{insert}', 1),
                            encoding='utf-8')
            else:
                with open(ignored_file, encoding='utf-8') as fh:
                    for line in fh:
                        fname = line.strip()
                        fpath = os.path.join(ialdir_full, 'src', 'local', fname)
                        if os.path.isfile(fpath):
                            os.remove(fpath)
                dummy_dir = os.path.join(ialdir_full, 'src', 'local', 'mpa', 'dummy')
                if os.path.isdir(dummy_dir) and not os.listdir(dummy_dir):
                    os.rmdir(dummy_dir)

        shutil.rmtree(os.path.join(phyex_dir, 'PHYEX'), ignore_errors=True)

        # Incremental pack with updated Hub, we must put, in src/local, the files using
        # PHYEX to force a re-compilation
        src = os.path.join(ialdir_full, 'src')
        phyex_hub = os.path.join(ialdir_full, 'hub', 'local', 'src', 'PHYEX', 'phyex')
        src_main = os.path.join(ialdir_full, 'src', 'main')
        if os.path.isdir(phyex_hub) and os.path.isdir(src_main):
            grep_files = []
            for subdir in ('main/mpa', 'main/arpifs'):
                sub_path = os.path.join(src, subdir)
                if not os.path.isdir(sub_path):
                    continue
                for root, dirs, files in os.walk(sub_path):
                    for fname in files:
                        if not fname.endswith(('.F90', '.h')):
                            continue
                        full = os.path.join(root, fname)
                        with open(full, 'rb') as fh:
                            content = fh.read().lower()
                            if b'phyex_t' in content or b'modd_nsv' in content:
                                grep_files.append(os.path.relpath(full, src))
            for fname in ('main/mpa/conv/externals/aro_conv_mnh.F90',
                          'main/arpifs/phys_dmn/suphmpa.F90'):
                if os.path.isfile(os.path.join(src_main, fname)):
                    if fname not in grep_files:
                        grep_files.append(fname)
            for fpath in grep_files:
                localfile = re.sub(r'^main/', 'local/', fpath)
                localpath = os.path.join(src, localfile)
                if not os.path.isfile(localpath):
                    shutil.copy2(os.path.join(src, fpath), localpath)

    def compilation_step(self):
        """Compile the pack with *ics_masterodb* and check for errors."""
        ialdir_full = os.path.join(self.HOMEPACK, self.ialdir)
        if not self.onlyIfNeeded or not os.path.isfile(os.path.join(ialdir_full, 'bin', 'MASTERODB')):
            print(f"### Compilation of commit {self.commit}")

            # Compilation with 10 threads ('$': end of line, re.M)
            ics_file = Path(ialdir_full, 'ics_masterodb')
            ics_file.write_text(re.sub(r'GMK_THREADS=1$', 'GMK_THREADS=10',
                                       ics_file.read_text(encoding='utf-8'), flags=re.M),
                                encoding='utf-8')

            if os.path.isfile(os.path.join(ialdir_full, 'ics_packages')):
                self.submit(os.path.join(ialdir_full, 'Output_compilation_hub'),
                            ['./ics_packages'], cwd=ialdir_full)
            self.submit(os.path.join(ialdir_full, 'Output_compilation'),
                        ['./ics_masterodb'], cwd=ialdir_full)

    def execution(self):
        """Run each test case and optionally collect profiling data."""
        ialdir_full = os.path.join(self.HOMEPACK, self.ialdir)
        dirconf = os.path.join(self.phytoolsdir, 'conf_tests')

        if not self.onlyIfNeeded:
            for t in self.tests:
                test_dir = os.path.join(ialdir_full, 'conf_tests', t)
                if os.path.isdir(test_dir):
                    shutil.rmtree(test_dir)

        firstrun = True
        for t in self.tests:
            if t not in self.allowedTests:
                print(f"The test {t} is not allowed")
                continue

            test_dir = os.path.join(ialdir_full, 'conf_tests', t)
            if os.path.isdir(test_dir) and self.onlyIfNeeded:
                continue

            if firstrun:
                print(f"### Running of commit {self.commit}")
                firstrun = False

            if not os.path.isfile(os.path.join(ialdir_full, 'bin', 'MASTERODB')):
                raise CheckCommitError(
                    "Pack does not exist or compilation has failed, please check", 6)

            os.makedirs(test_dir, exist_ok=True)

            script_dir = os.path.join(dirconf, t)
            script = next((os.path.join(script_dir, name) for name in ('aro.sh', 'arp.sh')
                           if os.path.isfile(os.path.join(script_dir, name))), None)
            if script is None:
                print(f"No script found in {script_dir}")
                continue

            env = os.environ.copy()
            env['MYLIB'] = self.ialdir
            env['TESTDIR'] = script_dir
            with timed(self.perffile, self.name, 'ial', t):
                self.submit(os.path.join(test_dir, 'Output_run'), [script], cwd=test_dir,
                            env=env)

            # DrHook profiles of the different MPI tasks (drhook.prof.<n>): concatenation
            # (header of the first profile, data lines of all of them) and aggregation
            prof_files = sorted(f for f in os.listdir(test_dir)
                                if re.fullmatch(r'drhook\.prof\.\d+', f))
            if prof_files:
                concat_file = os.path.join(test_dir, 'drhook.prof.concat')
                with open(concat_file, 'w', encoding='utf-8') as dst:
                    for i, fname in enumerate(prof_files):
                        with open(os.path.join(test_dir, fname), encoding='utf-8') as src:
                            dst.writelines(line for line in src
                                           if i == 0 or re.match(r'^\s*\d', line))
                aggregate_drhook_profile(concat_file, os.path.join(test_dir, 'drhook.prof.agg'))

    def comparison(self):
        """Compare output files against the reference commit."""
        ialdir_full = os.path.join(self.HOMEPACK, self.ialdir)
        print(f"### Check commit {self.commit} against commit {self.reference}")
        refByTest = self.resolve_ref_per_test(self.json_content)

        # First part: build test info dict, collect missing refs by caseref
        test_info = {}
        missing_by_ref = {}
        for t in self.tests:
            if t in self.allowedTests:
                if 'small' in t:
                    filepath = f"conf_tests/{t}/ICMSHFPOS+0002:00"
                else:
                    filepath = f"conf_tests/{t}/NODE.001_01"

                caseref = refByTest.get(t, '')
                if '/' in caseref:
                    refname = f"PHYEX/*_{escape_commit(caseref)}.01.{self.gmkpack_l}.{self.gmkpack_o}"
                else:
                    refname = f"PHYEX/*_{caseref}.01.{self.gmkpack_l}.{self.gmkpack_o}"

                file1 = os.path.join(self.HOMEPACK, self.ialdir, filepath)
                ref_glob = os.path.join(self.HOMEPACK, refname, filepath)
                matched = glob.glob(ref_glob)
                file2 = matched[0] if matched else ref_glob

                test_info[t] = dict(filepath=filepath, file1=file1, file2=file2,
                                    caseref=caseref, refname=refname, ref_glob=ref_glob)
                if not os.path.isfile(file2):
                    missing_by_ref.setdefault(caseref, []).append(t)

        # Between parts: compute missing references in batch (one per caseref)
        if self.computeRefIfNeeded:
            for caseref, tests_list in missing_by_ref.items():
                run_command(
                    [sys.executable, '-m', 'pyphyextools.testing.check_ial',
                     '-p', '-c', '-r', '-t', ','.join(tests_list),
                     '--onlyIfNeeded', caseref])
                # Re-resolve file2 for the newly computed refs
                for t in tests_list:
                    info = test_info[t]
                    matched = glob.glob(info['ref_glob'])
                    if matched:
                        info['file2'] = matched[0]

        # Second part: compare each test
        allt = 0
        message = ""
        for t, info in test_info.items():
            file1 = info['file1']
            file2 = info['file2']
            filepath = info['filepath']
            mess = ""
            t_status = 0
            if not os.path.isfile(file1):
                mess = f"Result ({file1}) for commit {self.commit} does not exist, " \
                       "please run the simulation"
                t_status = 1
            if not os.path.isfile(file2):
                mess2 = f"Reference result ({file2}) for commit {info['caseref']} does not exist, " \
                        "please run the simulation"
                t_status = 1
                mess = mess2 if not mess else f"{mess} and {mess2}"

            if t_status == 0:
                basename = os.path.basename(filepath)
                if basename == 'ICMSHFPOS+0002:00':
                    t_status = comp_binary(file1, file2, 256)
                elif basename == 'NODE.001_01':
                    t_status = comp_NODE(file1, file2)
                else:
                    t_status = comp_binary(file1, file2, 0)
                mess = ""

            if t_status != 0:
                message += f" {filepath} : {mess}\n"
            allt += t_status

        if allt == 0:
            print("SUCCESS, files are (nearly) identical")
            return 0
        else:
            print("*************** Files are different *******************")
            print(message)
            return 50

    def cleaning(self):
        """Remove the pack directory."""
        ialdir_full = os.path.join(self.HOMEPACK, self.ialdir)
        print(f"### Remove model directory for commit {self.commit}")
        if os.path.isdir(ialdir_full):
            shutil.rmtree(ialdir_full)


def main():
    """Entry point for the ``phyex-check_ial`` command."""
    run_tool(CheckCommitIAL, module_name='pyphyextools.testing.check_ial')


if __name__ == '__main__':
    main()
