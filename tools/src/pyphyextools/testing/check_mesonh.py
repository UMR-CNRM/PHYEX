"""Commit-check tool for the MesoNH model.

Compiles MesoNH using a specific commit for the externalised physics,
runs tests and checks if results are identical to a given version.
"""

import glob
import json
import os
import shutil
import subprocess
import sys
import time

from pyphyextools.testing.check_common import (
    CheckCommitBase, CheckCommitError, escape_commit, mvdiff, run_tool)
from pyphyextools.prep_code import prep_code
from pyphyextools.testing.compare import compareBACKUPFiles, compareTSERIESFiles, comp_ncdump


class CheckCommitMesonh(CheckCommitBase):
    """Check a commit against the MesoNH reference."""

    default_expand = False

    def __init__(self, prepCodeOpts='', **kwargs):
        super().__init__(**kwargs)
        self.defaultTest = ["007_16janvier/008_run2"]
        self.allowedTests = ["007_16janvier/008_run2", "007_16janvier/008_run2_turb3D",
                             "007_16janvier/008_run2_lredf",
                             "COLD_BUBBLE/002_mesonh", "ARMLES/RUN",
                             "COLD_BUBBLE_3D/002_mesonh", "OCEAN_LES/004_run2",
                             "014_LIMA/002_mesonh"]
        self.prepCodeOpts = prepCodeOpts
        self.MNHPACK = os.environ.get('MNHPACK', os.path.expanduser('~/MesoNH/PHYEX'))
        self.TARGZDIR = os.environ.get(
            'TARGZDIR',
            os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'pack')))
        self.refversion = None
        self.mnhdir = None

        if os.path.isdir(os.path.join(self.commit, 'src')):
            self.model_ready = False
            json_path = os.path.join(self.commit, 'src', 'mesonh', 'mesonh_version.json')
        else:
            self.model_ready = True
            json_path = os.path.join(self.commit, 'mesonh_version.json')
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

        self.refversion = self.json_content.get('refversion', '')
        self.mnhdir = f"{self.refversion}-{self.name}"

    @classmethod
    def parse_arguments(cls):
        """Build the argument parser, add MesoNH-specific args, parse and return the dict."""
        parser = super().parse_arguments()
        parser.description = (
            "Compile the MESONH model using a specific commit for the "
            "externalised physics, run tests and check if results are "
            "identical to a given version")
        parser.add_argument('--prep_code-opts', dest='prepCodeOpts', default='',
                            help='options added to prep_code call')
        args = parser.parse_args()
        return vars(args), cls._build_commitcmd(parser)

    def pack_creation(self):
        """Create the MesoNH pack by extracting a tar.gz archive."""
        mnhdir_full = os.path.join(self.MNHPACK, self.mnhdir)

        if self.suppress and os.path.isdir(mnhdir_full):
            shutil.rmtree(mnhdir_full)

        if os.path.isdir(mnhdir_full) and self.onlyIfNeeded:
            self.packcreation = False  # To prevent pack_update from being executed
            return

        if os.path.isdir(mnhdir_full):
            raise CheckCommitError(
                f"Pack already exists ({mnhdir_full}),\n"
                "suppress it to be able to compile it again "
                "(or use the -s option to automatically suppress it)", 5)

        print(f"### Pack creation for commit {self.commit}")

        tmpdir = os.path.join(self.MNHPACK, f"{self.mnhdir}_$$")
        os.makedirs(tmpdir, exist_ok=True)
        try:
            tarball = os.path.join(self.TARGZDIR, f"{self.refversion}.tar.gz")
            if not os.path.isfile(tarball):
                raise CheckCommitError(f"Tarball not found: {tarball}", 8)
            subprocess.run(['cp', tarball, tmpdir], check=True)
            subprocess.run(['tar', 'xfz', os.path.join(tmpdir, f"{self.refversion}.tar.gz")],
                           cwd=tmpdir, check=True)
            os.remove(os.path.join(tmpdir, f"{self.refversion}.tar.gz"))
            shutil.move(os.path.join(tmpdir, self.refversion), mnhdir_full)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def pack_update(self):
        """Copy and prepare source files with *prep_code*."""
        mnhdir_full = os.path.join(self.MNHPACK, self.mnhdir)
        src_phyex = os.path.join(mnhdir_full, 'src', 'PHYEX')

        if not os.path.isdir(src_phyex):
            raise CheckCommitError(
                f"PHYEX directory doesn't exist in pack ({mnhdir_full})", 9)

        src_base = os.path.join(mnhdir_full, 'src')

        if self.packupdate:
            shutil.move(os.path.join(src_base, 'PHYEX'),
                        os.path.join(src_base, 'PHYEXori'))
        else:
            shutil.rmtree(os.path.join(src_base, 'PHYEX'), ignore_errors=True)

        print(f"Copy {self.commit}")
        os.makedirs(os.path.join(src_base, 'PHYEX'), exist_ok=True)

        if not self.model_ready:
            subprocess.run(['scp', '-q', '-r', os.path.join(self.commit, 'src'),
                            os.path.join(src_base, 'PHYEX') + '/'],
                           check=True)
            prep_kwargs = self._parse_prep_code_opts(self.prepCodeOpts)
            prep_code(
                directory=os.path.join(src_base, 'PHYEX'),
                model='mesonh',
                mnh_expand=self.useexpand,
                subs=['turb', 'micro', 'aux', 'ext', 'conv'],
                rename_Ff_flag=True,
                ilooprm=True,
                no_raise_on_coding_norms=True,
                pyfortool_options=['--removeExtraDOinMnhDoConcurrent'],
                **prep_kwargs)
        else:
            prep_kwargs = self._parse_prep_code_opts(self.prepCodeOpts)
            prep_code(directory=os.path.join(src_base, 'PHYEX'), **prep_kwargs)
            git_dir = os.path.join(src_base, 'PHYEX', '.git')
            if os.path.isdir(git_dir):
                shutil.rmtree(git_dir)

        subprocess.run(['find', 'PHYEX', '-type', 'f', '-exec', 'touch', '{}', ';'],
                       cwd=src_base, check=True)

        phyex_dir = os.path.join(src_base, 'PHYEX')
        if os.path.isdir(os.path.join(phyex_dir, 'ext')):
            for fname in os.listdir(os.path.join(phyex_dir, 'ext')):
                fpath = os.path.join(phyex_dir, 'ext', fname)
                if not os.path.isfile(fpath):
                    continue
                if fname == 'modd_salt.f90':
                    aclib_dir = os.path.join(mnhdir_full, 'src', 'ACLIB', 'aux')
                    os.makedirs(aclib_dir, exist_ok=True)
                    mvdiff(fpath, aclib_dir)
                else:
                    mvdiff(fpath, os.path.join(mnhdir_full, 'src', 'MNH'))
            if self.packupdate:
                shutil.rmtree(os.path.join(phyex_dir, 'ext'))
            else:
                os.rmdir(os.path.join(phyex_dir, 'ext'))

        if self.packupdate:
            phyexori_dir = os.path.join(src_base, 'PHYEXori')
            for subdir in ('turb', 'micro', 'conv', 'aux'):
                subdir_path = os.path.join(phyex_dir, subdir)
                if not os.path.isdir(subdir_path):
                    continue
                for fname in os.listdir(subdir_path):
                    fpath = os.path.join(subdir_path, fname)
                    if os.path.isfile(fpath):
                        mvdiff(fpath, os.path.join(phyexori_dir, subdir, fname))
            shutil.rmtree(phyex_dir)
            shutil.move(phyexori_dir, phyex_dir)
        else:
            for subdir in ('turb', 'micro', 'conv', 'aux'):
                subdir_path = os.path.join(phyex_dir, subdir)
                for f in os.listdir(subdir_path):
                    if f.endswith('.f90'):
                        target_path = os.path.join(mnhdir_full, 'src', 'MNH', f)
                        if os.path.isfile(target_path):
                            os.remove(target_path)

            for subdir in ('turb', 'micro', 'conv'):
                for f in os.listdir(os.path.join(phyex_dir, subdir)):
                    if f.startswith('mode_'):
                        filepath = os.path.join(mnhdir_full, 'src', 'MNH', f)
                        if os.path.exists(filepath):
                            os.remove(filepath)

        exe_dir = os.path.join(mnhdir_full, 'exe')
        if os.path.isdir(exe_dir):
            for f in os.listdir(exe_dir):
                os.remove(os.path.join(exe_dir, f))

        for t in self.allowedTests:
            case1 = t.split('/')[0]
            case2 = t.split('/')[1]
            casedir = os.path.join(mnhdir_full, 'MY_RUN', 'KTEST', case1, case2)
            if os.path.isdir(casedir):
                shutil.rmtree(casedir)

    def compilation_step(self):
        """Compile MesoNH with *make*."""
        mnhdir_full = os.path.join(self.MNHPACK, self.mnhdir)
        exe_dir = os.path.join(mnhdir_full, 'exe')
        need_comp = not self.onlyIfNeeded
        if self.onlyIfNeeded:
            need_comp = not glob.glob(os.path.join(exe_dir, 'MESONH*'))

        if need_comp:
            print(f"### Compilation of commit {self.commit}")
            src_dir = os.path.join(mnhdir_full, 'src')

            subprocess.run(['./configure'], cwd=src_dir, check=True)

            for f in os.listdir(exe_dir):
                os.remove(os.path.join(exe_dir, f))
            shell_cmd = (
                f". {mnhdir_full}/conf/profile_mesonh-* 2>/dev/null; "
                f"make -j 8 && make installmaster"
            )
            self._run_with_tee(['bash', '-c', shell_cmd],
                               os.path.join(mnhdir_full, 'src'),
                               os.path.join(mnhdir_full, 'Output_compilation'))

    def execution(self):
        """Run each test case."""
        mnhdir_full = os.path.join(self.MNHPACK, self.mnhdir)

        if not self.onlyIfNeeded:
            for t in self.tests:
                if t not in self.allowedTests:
                    continue
                case1 = t.split('/')[0]
                case2 = t.split('/')[1]
                casedir = os.path.join(mnhdir_full, 'MY_RUN', 'KTEST', case1, case2)
                if os.path.isdir(casedir):
                    shutil.rmtree(casedir)

        firstrun = 1
        for t in self.tests:
            if t not in self.allowedTests:
                print(f"The test {t} is not allowed")
                continue

            case1 = t.split('/')[0]
            case2 = t.split('/')[1]
            casedir = os.path.join(mnhdir_full, 'MY_RUN', 'KTEST', case1)
            full_casedir = os.path.join(casedir, case2)

            if os.path.isdir(full_casedir) and self.onlyIfNeeded:
                continue

            if firstrun:
                print(f"### Running of commit {self.commit}")
                firstrun = 0

            exe_glob = os.path.join(mnhdir_full, 'exe', 'MESONH*')
            if not glob.glob(exe_glob):
                raise CheckCommitError(
                    f"Pack does not exist ({mnhdir_full}) or compilation has failed, "
                    "please check", 6)

            if not os.path.isdir(casedir):
                ref_casedir = os.path.join(self.MNHPACK, f"{self.refversion}",
                                           'MY_RUN', 'KTEST', case1)
                if os.path.isdir(ref_casedir):
                    shutil.copytree(ref_casedir, casedir, symlinks=True)
                for newt in self.allowedTests:
                    nc1 = newt.split('/')[0]
                    nc2 = newt.split('/')[1]
                    if case1 == nc1:
                        d = os.path.join(casedir, nc2)
                        if os.path.isdir(d):
                            shutil.rmtree(d)

            for d in os.listdir(casedir):
                d_path = os.path.join(casedir, d)
                if os.path.isdir(d_path) or (os.path.islink(d_path) and not os.path.exists(d_path)):
                    if f"{case1}/{d}" not in self.allowedTests:
                        if os.path.islink(d_path):
                            os.unlink(d_path)
                        else:
                            shutil.rmtree(d_path, ignore_errors=True)
                        ref_link = os.path.join(
                            self.MNHPACK, self.refversion, 'MY_RUN', 'KTEST', case1, d)
                        if os.path.exists(ref_link):
                            os.symlink(ref_link, d_path)

            shutil.copytree(os.path.join(self.MNHPACK, self.refversion,
                                         'MY_RUN', 'KTEST', case1, case2),
                            full_casedir, symlinks=True)

            subprocess.run(['./clean_mesonh_xyz'], cwd=full_casedir, check=False)
            t1_us = time.time_ns() // 1000
            shell_cmd = (
                f". {mnhdir_full}/conf/profile_mesonh-* 2>/dev/null; "
                f"./run_mesonh_xyz"
            )
            self._run_with_tee(['bash', '-c', shell_cmd], full_casedir,
                               os.path.join(full_casedir, 'Output_run'))
            t2_us = time.time_ns() // 1000
            elapsed_ms = (t2_us - t1_us) // 1000

            if self.perffile:
                with open(self.perffile, 'a', encoding='utf-8') as pf:
                    pf.write(f"{self.commit} mesonh {t} {elapsed_ms}\n")

    def comparison(self):
        """Compare NetCDF output files against the reference."""
        mnhdir_full = os.path.join(self.MNHPACK, self.mnhdir)
        print(f"### Check commit {self.commit} against commit {self.reference}")
        refByTest = self.resolve_ref_per_test(self.json_content)

        result = subprocess.run(
            ['bash', '-c',
             f'. {mnhdir_full}/conf/profile_mesonh-* 2>/dev/null; echo $HDF5_PLUGIN_PATH'],
            capture_output=True, text=True, check=False)
        hdf5_path = result.stdout.strip()
        if hdf5_path:
            os.environ['HDF5_PLUGIN_PATH'] = hdf5_path

        # First part: build test info dict, collect missing refs by caseref
        test_info = {}
        missing_by_ref = {}
        for t in self.tests:
            if t not in self.allowedTests:
                continue

            caseref = refByTest.get(t, '')
            if '/' in caseref:
                refname = f"{self.refversion}-{escape_commit(caseref)}"
            else:
                refname = f"{self.refversion}-{caseref}"

            case1 = t.split('/')[0]
            case2 = t.split('/')[1]
            path_user = os.path.join(mnhdir_full, 'MY_RUN', 'KTEST', case1, case2)
            path_ref = os.path.join(self.MNHPACK, refname, 'MY_RUN', 'KTEST', case1, case2)

            file1 = file2 = file3 = file4 = None
            bit_diff = 0

            if case1 == '007_16janvier':
                file1 = os.path.join(path_user, '16JAN.1.12B18.001.nc')
                file2 = os.path.join(path_ref, '16JAN.1.12B18.001.nc')
                file3 = os.path.join(path_user, '16JAN.1.12B18.000.nc')
                file4 = os.path.join(path_ref, '16JAN.1.12B18.000.nc')
                bit_diff = 40500
            elif case1 == 'COLD_BUBBLE':
                file1 = os.path.join(path_user, 'BUBBL.1.CEN4T.001.nc')
                file2 = os.path.join(path_ref, 'BUBBL.1.CEN4T.001.nc')
                bit_diff = 27300
            elif case1 == 'OCEAN_LES':
                file1 = os.path.join(path_user, 'SPWAN.2.25m00.001.nc')
                file2 = os.path.join(path_ref, 'SPWAN.2.25m00.001.nc')
                bit_diff = 18400
            elif case1 == 'COLD_BUBBLE_3D':
                file1 = os.path.join(path_user, 'BUBBL.1.CEN4T.001.nc')
                file2 = os.path.join(path_ref, 'BUBBL.1.CEN4T.001.nc')
                file3 = os.path.join(path_user, 'BUBBL.1.CEN4T.000.nc')
                file4 = os.path.join(path_ref, 'BUBBL.1.CEN4T.000.nc')
                bit_diff = 27300
            elif case1 == 'ARMLES':
                file1 = os.path.join(path_user, 'ARM__.1.CEN4T.001.nc')
                file2 = os.path.join(path_ref, 'ARM__.1.CEN4T.001.nc')
                file3 = os.path.join(path_user, 'ARM__.1.CEN4T.000.nc')
                file4 = os.path.join(path_ref, 'ARM__.1.CEN4T.000.nc')
                bit_diff = 71100
            elif case1 == '014_LIMA':
                file1 = os.path.join(path_user, 'XPREF.1.SEG01.002.nc')
                file2 = os.path.join(path_ref, 'XPREF.1.SEG01.002.nc')
                file3 = os.path.join(path_user, 'XPREF.1.SEG01.000.nc')
                file4 = os.path.join(path_ref, 'XPREF.1.SEG01.000.nc')
                bit_diff = 32200
            else:
                print(f"Unknown case: {t}")
                continue

            test_info[t] = dict(refname=refname, case1=case1, case2=case2,
                                path_user=path_user, path_ref=path_ref,
                                file1=file1, file2=file2, file3=file3,
                                file4=file4, bit_diff=bit_diff, caseref=caseref)
            if not os.path.isfile(file2):
                missing_by_ref.setdefault(caseref, []).append(t)

        # Between parts: compute missing references in batch (one per caseref)
        if self.computeRefIfNeeded:
            for caseref, tests_list in missing_by_ref.items():
                env = os.environ.copy()
                env['MNHPACK'] = self.MNHPACK
                env['TARGZDIR'] = self.TARGZDIR
                env['PHYEXREPOuser'] = self.repo_user
                env['PHYEXREPOprotocol'] = self.repo_protocol
                subprocess.run(
                    [sys.executable, '-m', 'pyphyextools.testing.check_mesonh',
                     '-p', '-c', '-r', '-t', ','.join(tests_list),
                     '--onlyIfNeeded', caseref],
                    env=env, check=True)

        # Second part: compare each test
        allt = 0
        for t, info in test_info.items():
            path_user = info['path_user']
            path_ref = info['path_ref']
            file1 = info['file1']
            file2 = info['file2']
            file3 = info['file3']
            file4 = info['file4']
            bit_diff = info['bit_diff']

            if not os.path.isdir(path_user):
                raise CheckCommitError(
                    f"{path_user} is missing, please run the simulation", 7)

            if not os.path.isdir(path_ref):
                raise CheckCommitError(
                    f"{path_ref} is missing, please run the reference simulation", 8)

            if os.path.isfile(file1) and os.path.isfile(file2):
                print(f"Comparison for case {t}...")
                if file3 is None:
                    r = compareBACKUPFiles(file1, file2)
                else:
                    r = compareBACKUPFiles(file1, file2)
                    r += compareTSERIESFiles(file3, file4)
                allt += r

                r = comp_ncdump(file1, file2, str(bit_diff))
                allt += r
            else:
                if not os.path.isfile(file1):
                    print(f"  {file1} is missing")
                if not os.path.isfile(file2):
                    print(f"  {file2} is missing")
                allt += 1

        if allt == 0:
            status = "OK"
            print(f"...comparison done: {status}")
            return 0
        else:
            status = "Files are different"
            print(f"...comparison done: {status}")
            return 50

    def cleaning(self):
        """Remove the MesoNH pack directory."""
        mnhdir_full = os.path.join(self.MNHPACK, self.mnhdir)
        print(f"### Remove model directory for commit {self.commit}")
        if os.path.isdir(mnhdir_full):
            shutil.rmtree(mnhdir_full)


def main():
    """Entry point for the ``phyex-check_mesonh`` command."""
    run_tool(CheckCommitMesonh, module_name='pyphyextools.testing.check_mesonh')


if __name__ == '__main__':
    main()
