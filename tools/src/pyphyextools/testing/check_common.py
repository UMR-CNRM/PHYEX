"""Base class for PHYEX commit-checking tools.

Provides the shared infrastructure (argument parsing, cloning, property
resolution, step orchestration) that each model-specific tool specialises.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

SEPARATOR = '_'
PHYEXCONF_DEFAULT = os.path.expanduser('~/.phyex')


def _parse_test_list(s):
    """Parse a comma-separated test list from CLI into a list of stripped strings."""
    return [t.strip() for t in s.split(',') if t.strip()]


class CheckCommitError(Exception):
    """Exception raised by commit-check tools to signal an error with a specific exit status."""

    def __init__(self, message, exitcode=1):
        super().__init__(message)
        self.exitcode = exitcode


def escape_commit(commit):
    """Replace ``/``, ``:`` and ``.`` with *SEPARATOR* for use in directory names."""
    for c in ('/', ':', '.'):
        commit = commit.replace(c, SEPARATOR)
    return commit


def mvdiff(src, dst):
    """Move *src* to *dst* only if the files differ; otherwise delete *src*."""
    if os.path.isdir(dst):
        dst = os.path.join(dst, os.path.basename(src))
    if not os.path.exists(dst):
        shutil.move(src, dst)
    else:
        with open(src, 'rb') as f_src, open(dst, 'rb') as f_dst:
            if f_src.read() != f_dst.read():
                shutil.move(src, dst)
            else:
                os.remove(src)


class CheckCommitBase:
    """Base class for commit-checking tools.

    Subclasses override :meth:`pack_creation`, :meth:`pack_update`,
    :meth:`compilation_step`, :meth:`execution`, :meth:`comparison` and
    :meth:`cleaning` to provide model-specific behaviour.
    """

    default_expand = None

    @staticmethod
    def _run_with_tee(cmd, cwd, out_path, mode='w', env=None,
                      input_string=None):
        """Run *cmd* in *cwd*, teeing output to both terminal and *out_path*."""
        with open(out_path, mode, encoding='utf-8') as out:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, bufsize=1, text=True,
                                    cwd=cwd, env=env,
                                    stdin=subprocess.PIPE if input_string is not None else None)
            if input_string is not None:
                proc.stdin.write(input_string)
                proc.stdin.close()
            for line in proc.stdout:
                print(line, end='')
                out.write(line)
            proc.wait()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)

    def __init__(self, packcreation, packupdate, compilation, run_tests, check,
                 remove, suppress, onlyIfNeeded, computeRefIfNeeded, perffile,
                 name, repo_user, repo_protocol, tests, commit, reference,
                 useexpand=None, module_name=None, commitcmd=None, **kwargs):
        self.packcreation = packcreation
        self.packupdate = packupdate
        self.compilation = compilation
        self.run_tests = run_tests
        self.check = check
        self.remove = remove
        self.suppress = suppress
        self.onlyIfNeeded = onlyIfNeeded
        self.computeRefIfNeeded = computeRefIfNeeded
        self.perffile = perffile
        self.name = name
        self.repo_user = repo_user
        self.repo_protocol = repo_protocol
        self.tests = tests
        self.commit = commit
        self.reference = reference

        if self.default_expand is True:
            self.useexpand = useexpand if useexpand is not None else True
        elif self.default_expand is False:
            self.useexpand = useexpand if useexpand is not None else False
        else:
            self.useexpand = True

        if not any([self.packcreation, self.packupdate,
                     self.compilation, self.run_tests, self.check, self.remove]):
            self.packcreation = True
            self.compilation = True
            self.run_tests = True
            self.check = True

        if self.commit is None:
            raise CheckCommitError(
                "At least one commit hash must be provided on command line", 2)

        if self.check and self.reference is None:
            raise CheckCommitError(
                "To perform a comparison two commit hashes are mandatory "
                "on the command line", 3)

        self.commitcmd = commitcmd if commitcmd is not None else []
        self._module_name = module_name
        self.json_content = ""
        self.model_ready = False

    @classmethod
    def parse_arguments(cls):
        """Build the argument parser, add common arguments and return it.

        Subclasses override to add model-specific arguments and call
        ``super().parse_arguments()`` to get the parser, then call
        ``parse_args()`` and return ``vars(...)``.
        """
        parser = argparse.ArgumentParser()
        cls._add_common_arguments(parser)
        return parser

    @classmethod
    def _add_common_arguments(cls, parser):
        parser.add_argument('-s', action='store_true', dest='suppress',
                            help='suppress compilation directory')
        parser.add_argument('-p', action='store_true', dest='packcreation',
                            help='create compilation directory')
        parser.add_argument('-u', action='store_true', dest='packupdate',
                            help='update compilation directory (experimental)')
        parser.add_argument('-c', action='store_true', dest='compilation',
                            help='performs compilation')
        parser.add_argument('-r', action='store_true', dest='run_tests',
                            help='runs the tests')
        parser.add_argument('-C', action='store_true', dest='check',
                            help='checks the result against the reference')
        parser.add_argument('--remove', action='store_true', dest='remove',
                            help='removes the compilation directory')
        parser.add_argument('-t', dest='tests', default=None, type=_parse_test_list,
                            help='comma separated list of tests to execute, or ALL')
        parser.add_argument('--onlyIfNeeded', action='store_true', dest='onlyIfNeeded',
                            help='perform step only if not already done')
        parser.add_argument('--computeRefIfNeeded', action='store_true',
                            dest='computeRefIfNeeded',
                            help='compute missing references')
        parser.add_argument('--perf', dest='perffile', default=None,
                            help='add performance statistics in file FILE')
        parser.add_argument('--name', dest='name', default=None,
                            help='specific name for the directory')
        parser.add_argument('--repo-user', dest='repo_user',
                            default=os.environ.get('PHYEXREPOuser', 'UMR-CNRM'),
                            help='GitHub user (default: PHYEXREPOuser env)')
        parser.add_argument('--repo-protocol', dest='repo_protocol',
                            default=os.environ.get('PHYEXREPOprotocol', 'https'),
                            choices=['https', 'ssh'],
                            help='protocol (https or ssh)')
        if cls.default_expand is True:
            parser.add_argument('--noexpand', action='store_false', dest='useexpand',
                                help='do not expand mnh_expand blocks')
        elif cls.default_expand is False:
            parser.add_argument('--expand', action='store_true', dest='useexpand',
                                help='expand mnh_expand blocks')
        parser.add_argument('commit', nargs='?', default=None,
                            help='commit hash or a directory to test')
        parser.add_argument('reference', nargs='?', default=None,
                            help='commit hash, directory, or REF')

    @classmethod
    def _build_commitcmd(cls, parser):
        """Build the list of flag arguments from *sys.argv* for subprocess re-invocation.

        Uses argparse action types to correctly identify which flags take
        a value (``store``, ``append``) versus boolean flags (``store_true``,
        ``store_false``, ``store_const``).
        """
        commitcmd = []
        i = 1
        while i < len(sys.argv):
            arg = sys.argv[i]
            if not arg.startswith('-'):
                break
            commitcmd.append(arg)
            action = parser._option_string_actions.get(arg)
            if action is not None:
                takes_value = not isinstance(action, (
                    argparse._StoreTrueAction,
                    argparse._StoreFalseAction,
                    argparse._StoreConstAction))
                if takes_value and i + 1 < len(sys.argv):
                    commitcmd.append(sys.argv[i + 1])
                    i += 1
            elif i + 1 < len(sys.argv) and not sys.argv[i + 1].startswith('-'):
                commitcmd.append(sys.argv[i + 1])
                i += 1
            i += 1
        return commitcmd

    def clone_and_run(self):
        """If *commit* is a hash (no ``/`` outside ``tags/``), clone the
        repository, checkout the commit, install the tools package and
        re-invoke the same command with the cloned directory as argument."""
        needs_clone = not ('/' in self.commit and not self.commit.startswith('tags/'))
        if not needs_clone:
            return

        tmpdir = tempfile.mkdtemp()
        try:
            if self.repo_protocol == 'https':
                repo_url = f'https://github.com/{self.repo_user}/PHYEX.git'
            else:
                repo_url = f'git@github.com:{self.repo_user}/PHYEX.git'

            subprocess.run(['git', 'clone', repo_url], cwd=tmpdir, check=True)
            phyex_dir = os.path.join(tmpdir, 'PHYEX')
            subprocess.run(['git', 'checkout', self.commit], cwd=phyex_dir, check=True)

            if not os.path.isdir(os.path.join(phyex_dir, 'docs')):
                result = subprocess.run(
                    ['git', 'log', '--all', '--full-history', '--', 'tools'],
                    cwd=phyex_dir, capture_output=True, text=True, check=True)
                toolscommit = result.stdout.strip().split('\n')[0].split()[1]
                parents = subprocess.run(
                    ['git', 'rev-parse', f'{toolscommit}^@'],
                    cwd=phyex_dir, capture_output=True, text=True, check=True)
                for parent in parents.stdout.strip().split():
                    nfiles = subprocess.run(
                        ['git', 'ls-tree', '-r', parent, '--', 'tools'],
                        cwd=phyex_dir, capture_output=True, text=True, check=True)
                    if nfiles.stdout.count('\n') != 0:
                        toolscommit = parent
                        break
                for path in ('tools', 'pyproject.toml', 'src/pyphyex'):
                    subprocess.run(['git', 'checkout', toolscommit, '--', path],
                                   cwd=phyex_dir, check=True)

            subprocess.run([sys.executable, '-m', 'venv', 'venv'],
                           cwd=phyex_dir, check=True)
            pip_path = os.path.join(phyex_dir, 'venv', 'bin', 'pip')
            subprocess.run([pip_path, 'install', '-e', 'tools/'],
                           cwd=phyex_dir, check=True)

            python_path = os.path.join(phyex_dir, 'venv', 'bin', 'python')
            # Use the explicit module name stored by run_tool, or fall back
            # to the class's __module__ (which is correct under entry-point
            # invocation but is '__main__' under python -m).
            module = getattr(self, '_module_name', None) or self.__class__.__module__
            cmd = [python_path, '-m', module]
            cmd.extend(self.commitcmd)
            cmd.append(phyex_dir)
            if self.reference:
                cmd.append(self.reference)
            if '--name' not in cmd:
                cmd.extend(['--name', self.commit])

            result = subprocess.run(cmd, check=False)
            if result.returncode != 0:
                raise CheckCommitError(
                    f"Subprocess exited with status {result.returncode}",
                    result.returncode)
            sys.exit(0)

        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def _parse_prep_code_opts(self, opts_str):
        """Convert a *prep_code* option string into keyword arguments."""
        kwargs = {}
        if not opts_str:
            return kwargs
        for part in opts_str.split():
            if part == '--noRaiseOnCodingNorms':
                kwargs['no_raise_on_coding_norms'] = True
            elif part == '--renameFf':
                kwargs['rename_Ff_flag'] = True
            elif part == '--ilooprm':
                kwargs['ilooprm'] = True
            else:
                print(f"Warning: unknown prep_code option '{part}' ignored")
        return kwargs

    def resolve_ref_per_test(self, json_content):
        """Build and return a dict mapping each test name to its reference commit."""
        testing = json_content.get('testing', '')
        refByTest = {}
        if testing:
            for t in self.allowedTests:
                ref = testing.get(t, '')
                if ref:
                    refByTest[t] = ref
        if self.reference:
            tests = self.tests or self.defaultTest
            for t in tests:
                if self.reference == 'REF':
                    caseref = refByTest.get(t, '')
                else:
                    caseref = self.reference
                refByTest[t] = caseref
        return refByTest

    def run(self):
        """Drive the workflow: clone, pack, compile, run, compare, clean."""
        self.clone_and_run()

        if self.packcreation:
            self.pack_creation()

        if self.packcreation or self.packupdate:
            self.pack_update()

        if self.compilation:
            self.compilation_step()

        if self.run_tests:
            self.execution()

        if self.run_tests and self.perffile:
            self.performance_evaluation()

        if self.check:
            cmpstatus = self.comparison()
            if cmpstatus:
                raise CheckCommitError("Comparison step reported differences", cmpstatus)

        if self.remove:
            self.cleaning()

    def pack_creation(self):
        """Create the compilation directory (model-specific)."""
        raise NotImplementedError

    def pack_update(self):
        """Update the source code in the pack (model-specific)."""
        raise NotImplementedError

    def compilation_step(self):
        """Compile the pack (model-specific)."""
        raise NotImplementedError

    def execution(self):
        """Run the test cases (model-specific)."""
        raise NotImplementedError

    def comparison(self):
        """Compare results against the reference (model-specific)."""
        raise NotImplementedError

    def cleaning(self):
        """Remove the pack directory (model-specific)."""
        raise NotImplementedError

    def performance_evaluation(self):
        """Optional performance evaluation."""


def run_tool(tool_class, module_name=None):
    """Create a tool instance, parse arguments, set up and run.

    Parameters
    ----------
    tool_class : type
        A subclass of ``CheckCommitBase``.
    module_name : str or None
        Fully qualified module name for subprocess re-invocation
        (e.g. ``'pyphyextools.testing.check_ial'``).
    """
    kwargs, commitcmd = tool_class.parse_arguments()
    tool = tool_class(module_name=module_name, commitcmd=commitcmd, **kwargs)
    try:
        tool.run()
    except CheckCommitError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(e.exitcode)
