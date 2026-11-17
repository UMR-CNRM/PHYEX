"""Commit-check tool for the phyex1d 1D model.

Installs the pyphyex python package (python binding and 1D model) of a specific
commit in a virtual environment, runs the 1D model on DEPHY-SCM cases and checks
whether results are identical to a given version.
"""

import json
import os
import shutil
import sys
from pathlib import Path

from pyphyextools.testing.check_common import (
    CheckCommitBase, CheckCommitError, arch_file, copy_source, default_arch, escape_commit,
    PHYEXCONF, run_tool, SEPARATOR, sourced_env, timed)
from pyphyextools import run_command
from pyphyextools.testing.compare import comp_phyex1d

# Experiment run on each case (options of phyex1d-exp_run) and result file name
# (relative to the output directory of the case)
EXPERIMENT = ['--name', 'phyex1d', '--exp', 'dt=60', 'name=default']
RESULT_FILE = os.path.join('phyex1d', 'default.hdf5')


def case_tag(case):
    """Directory name corresponding to a case (CASE/SUBCASE)"""
    return case.strip('/').replace('/', SEPARATOR)


class CheckCommitPhyex1d(CheckCommitBase):
    """Check a commit against the phyex1d reference."""

    default_expand = None  # no source code transformation option
    default_shuman_trans = None

    def __init__(self, archfile=None, refarchfile=None, **kwargs):
        super().__init__(**kwargs)
        self.defaultTest = ['ARMCU/REF']
        self.TESTDIR = os.environ.get('PHYEX1DDIR', os.path.expanduser('~/PHYEX1D'))
        self.archfile = archfile or default_arch()
        self.refarchfile = refarchfile or default_arch()

        if os.path.isdir(os.path.join(self.commit, 'src')):
            json_path = os.path.join(self.commit, 'src', 'phyex', 'phyex1d_version.json')
        else:
            json_path = os.path.join(self.commit, 'phyex1d_version.json')
        if os.path.isfile(json_path):
            with open(json_path, encoding='utf-8') as f:
                self.json_content = json.load(f)
        else:
            self.json_content = {}

        # Any installed DEPHY-SCM case can be tested, ALL means all the cases having
        # a reference in the json file
        all_tests = list(self.json_content.get('testing', {}))
        if not self.tests:
            self.tests = list(self.defaultTest)
        elif 'ALL' in self.tests:
            expanded = []
            for t in self.tests:
                expanded.extend(all_tests if t == 'ALL' else [t])
            self.tests = expanded
        self.tests = [t.strip('/') for t in self.tests]
        self.allowedTests = list(dict.fromkeys(all_tests + self.tests))

    @classmethod
    def parse_arguments(cls):
        """Build the argument parser, add phyex1d-specific args, parse and return the dict."""
        parser = super().parse_arguments()
        parser.description = (
            "Install the pyphyex python package of a specific commit in a virtual "
            "environment, run the phyex1d 1D model on DEPHY-SCM cases (CASE/SUBCASE, "
            "ARMCU/REF by default, see phyex1d-cases_list) and check if results are "
            "identical to a given version. The cases are installed automatically in "
            f"{os.path.join(PHYEXCONF, 'phyex1d', 'DEPHY-SCM')} if needed.")
        parser.add_argument('-a', dest='archfile', default=None,
                            help='architecture name to build and run the commit (the files '
                                 'arch-<ARCH>.env and arch-<ARCH>.ecbuild are searched in '
                                 '${PHYEXCONF}/arch, then in the arch directory of the '
                                 'source code)')
        parser.add_argument('-A', dest='refarchfile', default=None,
                            help='architecture name for the reference simulation (files '
                                 'searched as for -a)')
        args = parser.parse_args()
        return vars(args), cls._build_commitcmd(parser)

    # Directories and files
    def _testdir(self, name=None):
        return os.path.join(self.TESTDIR, name or self.name)

    def _venv(self, name=None, archfile=None):
        """Virtual environment (next to the PHYEX source code), one per architecture; it
        also contains the architecture files, the build directory and the logs"""
        return os.path.join(self._testdir(name), f'venv_{archfile or self.archfile}')

    def _arch_env(self):
        return os.path.join(self._venv(), 'arch.env')

    def _bin(self, command):
        return os.path.join(self._venv(), 'bin', command)

    def _installed(self):
        return os.path.isfile(self._bin('phyex1d-exp_run'))

    def _test_out(self, case, name=None, archfile=None):
        return os.path.join(self._testdir(name), 'tests', f'arch_{archfile or self.archfile}',
                            case_tag(case))

    def pack_creation(self):
        testdir = self._testdir()
        venv = self._venv()

        if self.suppress and os.path.isdir(testdir):
            shutil.rmtree(testdir)

        if os.path.isdir(venv) and self.onlyIfNeeded:
            if self._installed():
                self.packcreation = False  # To prevent pack_update from being executed
                return
            # The package has not been installed: the directory has been left in an
            # unknown state by a previous attempt, it is created again
            print(f"WARNING: the virtual environment ({venv}) is incomplete, "
                  "it is suppressed to be created again")
            shutil.rmtree(venv)

        if os.path.isdir(venv):
            raise CheckCommitError(
                f"Directory already exists ({venv}),\n"
                "suppress it to be able to compile it again "
                "(or use the -s option to automatically suppress it)", 5)
        os.makedirs(testdir, exist_ok=True)

    def pack_update(self):
        testdir = self._testdir()
        if not os.path.isdir(testdir):
            raise CheckCommitError(f"Directory must exist ({testdir})", 9)

        # Source code, in the PHYEX subdirectory of testdir
        phyex_dir = os.path.join(testdir, 'PHYEX')
        if os.path.isdir(phyex_dir):
            shutil.rmtree(phyex_dir)
        if os.path.isdir(self.commit):
            print(f"Copy {self.commit}")
            copy_source(self.commit, phyex_dir)
        else:
            # Commit hash or tag
            if self.repo_protocol == 'https':
                repo_url = f'https://github.com/{self.repo_user}/PHYEX.git'
            else:
                repo_url = f'git@github.com:{self.repo_user}/PHYEX.git'
            print(f"Clone repository, and checkout commit {self.commit}")
            run_command(['git', 'clone', repo_url, phyex_dir])
            run_command(['git', 'checkout', self.commit], cwd=phyex_dir)
        if not os.path.isdir(os.path.join(phyex_dir, 'src', 'pyphyex')):
            raise CheckCommitError(f"The commit {self.commit} does not contain the pyphyex "
                                   "package (src/pyphyex)", 6)

        # Virtual environment, next to the source code, with the architecture files
        # (searched in ${PHYEXCONF}/arch then in the arch directory of the source code)
        venv = self._venv()
        if not os.path.isdir(venv):
            run_command([sys.executable, '-m', 'venv', venv])
        for ext in ('env', 'ecbuild'):
            shutil.copy(arch_file(self.archfile, phyex_dir, ext),
                        os.path.join(venv, f'arch.{ext}'))
        run_command([self._bin('pip'), 'install', '--upgrade', 'pip'],
                    out_path=os.path.join(venv, 'Output_pip_upgrade'))

    def compilation_step(self):
        venv = self._venv()
        if not os.path.isdir(venv):
            raise CheckCommitError(f"Virtual environment must exist ({venv}), "
                                   f"run '-p {self.commit}' to create it", 9)
        if not self.onlyIfNeeded or not self._installed():
            print(f"### Compilation of commit {self.commit}")
            # Installation of the package (compilation of the library and of the python
            # binding), executed in the environment defined by the environment file; the
            # build directory is kept for the following compilations
            cmd = [self._bin('pip'), 'install', '--prefer-binary',
                   '-C', 'build-dir=' + os.path.join(venv, 'build'),
                   '-C', 'cmake.define.CMAKE_TOOLCHAIN_FILE=' + os.path.join(venv, 'arch.ecbuild'),
                   os.path.join(self._testdir(), 'PHYEX')]
            run_command(cmd, venv, os.path.join(venv, 'Output_compilation'),
                        env=sourced_env(self._arch_env()), display='always')

    def _install_cases(self, env):
        """Installation (or update) of the DEPHY-SCM cases if some are missing"""
        result = run_command([self._bin('phyex1d-cases_list'), '-n'], env=env,
                             check=False, banner=False)
        available = result.stdout.split() if result.returncode == 0 else []
        if any(case not in available for case in self.tests):
            run_command([self._bin('phyex1d-cases_install'), 'download'], env=env,
                        display='always')

    def execution(self):
        if not self._installed():
            raise CheckCommitError(
                f"The pyphyex package is not installed in {self._venv()}, "
                f"run '-p -c {self.commit}' to compile.", 6)
        env = sourced_env(self._arch_env())
        cases_checked = False
        firstrun = True
        for case in self.tests:
            test_out = self._test_out(case)
            if os.path.isdir(test_out):
                if self.onlyIfNeeded and os.path.isfile(os.path.join(test_out, RESULT_FILE)):
                    continue
                shutil.rmtree(test_out)

            if firstrun:
                print(f"### Running of commit {self.commit}")
                firstrun = False
            if not cases_checked:
                self._install_cases(env)
                cases_checked = True

            os.makedirs(test_out)
            with timed(self.perffile, self.name, 'phyex1d', case):
                result = run_command([self._bin('phyex1d-exp_run'), case, '-o', test_out]
                                     + EXPERIMENT,
                                     cwd=test_out, out_path=os.path.join(test_out, 'Output_run'),
                                     err_path=os.path.join(test_out, 'Stderr_run'),
                                     env=env, check=False, display='never')
                stderr_run = Path(test_out, 'Stderr_run').read_text(encoding='utf-8')
                if stderr_run.strip():
                    print(stderr_run, file=sys.stderr)
                if result.returncode != 0:  # no timing written for a failed execution
                    print(Path(test_out, 'Output_run').read_text(encoding='utf-8'))
                    raise CheckCommitError(f"Execution of case {case} failed "
                                           f"(see {test_out})", 7)

    def comparison(self):
        print(f"### Check commit {self.commit} against commit {self.reference}")
        refByTest = self.resolve_ref_per_test(self.json_content)

        # Result files, missing references collected by reference commit
        test_info = {}
        missing_by_ref = {}
        for case in self.tests:
            caseref = refByTest.get(case, '')
            refname = escape_commit(self.reference) if '/' in caseref else caseref
            file1 = os.path.join(self._test_out(case), RESULT_FILE)
            file2 = os.path.join(self._test_out(case, refname, self.refarchfile), RESULT_FILE)
            test_info[case] = dict(file1=file1, file2=file2, caseref=caseref)
            if not os.path.isfile(file2) and caseref:
                missing_by_ref.setdefault(caseref, []).append(case)

        # Computation of the missing references (one call per reference commit)
        if self.computeRefIfNeeded:
            for caseref, cases in missing_by_ref.items():
                run_command(
                    [sys.executable, '-m', 'pyphyextools.testing.check_phyex1d',
                     '-p', '-c', '-r', '-t', ','.join(cases),
                     '-a', self.refarchfile, '--onlyIfNeeded', caseref])

        alltests = 0
        message = ""
        for case, info in test_info.items():
            mess = []
            if not info['caseref']:
                mess.append(f"No reference for case {case}")
            if not os.path.isfile(info['file1']):
                mess.append(f"Result ({info['file1']}) for commit {self.commit} does not "
                            "exist, please run the simulation")
            if info['caseref'] and not os.path.isfile(info['file2']):
                mess.append(f"Result ({info['file2']}) for commit {info['caseref']} does "
                            "not exist, please run the simulation")
            if mess:
                message += ' ' + ' and '.join(mess) + '\n'
                alltests += 1
            elif comp_phyex1d(info['file1'], info['file2']) != 0:
                message += f" Results differ for case {case}\n"
                alltests += 1

        if alltests == 0:
            print("SUCCESS, files are identical")
            return 0
        print("*************** Files are different *******************")
        print(message)
        return 50

    def cleaning(self):
        testdir = self._testdir()
        print(f"### Remove model directory for commit {self.commit}")
        if os.path.isdir(testdir):
            shutil.rmtree(testdir)


def main():
    """Entry point for the ``phyex-check_phyex1d`` command."""
    run_tool(CheckCommitPhyex1d, module_name='pyphyextools.testing.check_phyex1d')


if __name__ == '__main__':
    main()
