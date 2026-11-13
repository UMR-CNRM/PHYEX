"""Base class for PHYEX commit-checking tools.

Provides the shared infrastructure (argument parsing, cloning, property
resolution, step orchestration) that each model-specific tool specialises.
"""

import argparse
import contextlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy
import pandas

from pyphyextools import run_command

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

    default_expand = True
    default_shuman_trans = None

    def __init__(self, packcreation, packupdate, compilation, run_tests, check,
                 remove, suppress, onlyIfNeeded, computeRefIfNeeded, perffile,
                 name, repo_user, repo_protocol, tests, commit, reference,
                 useexpand=None, useshumantrans=None, module_name=None, commitcmd=None, **kwargs):
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
        self.repo_user = repo_user
        self.repo_protocol = repo_protocol
        self.tests = tests
        self.commit = commit
        self.name = escape_commit(self.commit) if not name else name
        self.reference = reference
        self.useexpand = useexpand if useexpand is not None else self.default_expand
        self.useshumantrans = useshumantrans if useshumantrans is not None \
                              else self.default_shuman_trans

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
        if cls.default_shuman_trans is True:
            parser.add_argument('--noshumantrans', action='store_false', dest='useshumantrans',
                                help='do not transform shuman functions into subroutines')
        elif cls.default_shuman_trans is False:
            parser.add_argument('--shumantrans', action='store_true', dest='useshumantrans',
                                help='transform shuman functions into subroutines')
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

            run_command(['git', 'clone', repo_url], cwd=tmpdir)
            phyex_dir = os.path.join(tmpdir, 'PHYEX')
            run_command(['git', 'checkout', self.commit], cwd=phyex_dir)

            if not os.path.isdir(os.path.join(phyex_dir, 'docs')):
                toolscommit = run_command(
                    ['git', 'log', '-1', '--format=%H', '--all', '--full-history', '--', 'tools'],
                    cwd=phyex_dir).stdout.strip()
                parents = run_command(
                    ['git', 'rev-parse', f'{toolscommit}^@'],
                    cwd=phyex_dir)
                for parent in parents.stdout.strip().split():
                    if run_command(['git', 'ls-tree', '-r', parent, '--', 'tools'],
                                   cwd=phyex_dir).stdout.strip():
                        toolscommit = parent
                        break
                for path in ('tools', 'pyproject.toml', 'src/pyphyex'):
                    run_command(['git', 'checkout', toolscommit, '--', path],
                                cwd=phyex_dir)

            run_command([sys.executable, '-m', 'venv', 'venv'], cwd=phyex_dir)
            pip_path = os.path.join(phyex_dir, 'venv', 'bin', 'pip')
            run_command([pip_path, 'install', '--prefer-binary', '-e', 'tools/'],
                        cwd=phyex_dir)

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
            elif part == '--genModi':
                kwargs['gen_modi'] = True
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


_SOURCED_ENV_MARKER = '#### PHYEX ENVIRONMENT AFTER SOURCING ####'
_SOURCED_ENV_CACHE = {}


def sourced_env(script, cwd=None):
    """Return the environment (dict) obtained by sourcing the shell file *script*.

    The file is sourced by bash with ``set -a``: all the variables it sets are
    exported, even if they are not explicitly exported (as when the file is
    sourced in a shell script using these variables). The output of the
    sourcing is displayed; failing commands in the file only produce a
    warning (environment files are known to produce harmless errors). The
    result is cached (for a given file, modification time and *cwd*), a copy
    is returned.

    Parameters
    ----------
    script : str
        Environment file to source (it must exist).
    cwd : str or None
        Directory in which the file is sourced.

    Returns
    -------
    dict
        Environment to give to the commands (``env`` argument of
        :func:`~pyphyextools.run_command`).
    """
    script = os.path.abspath(script)
    if not os.path.isfile(script):
        raise CheckCommitError(f"Environment file not found: {script}", 9)
    key = (script, os.stat(script).st_mtime_ns, cwd)
    if key not in _SOURCED_ENV_CACHE:
        # The environment is printed (after a marker, following the output of the
        # sourcing) by the python interpreter executed in the sourced environment
        dump = ('import json, os, sys; '
                f'print({_SOURCED_ENV_MARKER!r}); print(sys.argv[1]); '
                'print(json.dumps(dict(os.environ)))')
        # The failing commands of the file are counted with an ERR trap. The
        # arguments are saved in variables (__phyex_*) and removed before
        # sourcing (some environment files modify the positional parameters)
        result = run_command(['bash', '-c', '__phyex_file=$1; __phyex_python=$2; __phyex_dump=$3; '
                                            'set --; __phyex_errors=0; '
                                            'trap \'__phyex_errors=$((__phyex_errors+1))\' ERR; '
                                            'set -a; . "$__phyex_file"; set +a; trap - ERR; '
                                            'exec "$__phyex_python" -c "$__phyex_dump" '
                                            '"$__phyex_errors"',
                              'bash', script, sys.executable, dump],
                             cwd=cwd, display='never', banner=False, check=False)
        output, marker, env_part = result.stdout.rpartition(_SOURCED_ENV_MARKER + '\n')
        print(f'### . {script}')
        print(output if marker else result.stdout, end='')
        if not marker or result.returncode != 0:
            raise CheckCommitError(f"Environment file {script} cannot be sourced", 9)
        errors, env = env_part.split('\n', 1)
        if errors != '0':
            print(f"WARNING: {errors} command(s) failed during the sourcing of {script} "
                  "(ignored)")
        env = {var: value for var, value in json.loads(env).items()
               if not var.startswith('__phyex_')}
        # Variables set by bash itself
        for var in ('_', 'PWD', 'OLDPWD', 'SHLVL'):
            env.pop(var, None)
            if var in os.environ:
                env[var] = os.environ[var]
        _SOURCED_ENV_CACHE[key] = env
    return dict(_SOURCED_ENV_CACHE[key])


@contextlib.contextmanager
def timed(perffile, *fields):
    """Context manager timing the block (elapsed time, measured with
    time.perf_counter) and, if *perffile* is given, appending to it a line
    made of *fields* followed by the elapsed time in ms. Nothing is written
    if the block raises an exception.

    Example::

        with timed(self.perffile, self.name, 'ial', test):
            ...
    """
    start = time.perf_counter()
    yield
    elapsed_ms = int((time.perf_counter() - start) * 1000)
    if perffile:
        with open(perffile, 'a', encoding='utf-8') as pf:
            pf.write(' '.join(map(str, fields)) + f' {elapsed_ms}\n')


def aggregate_drhook_profile(prof_file, agg_file):
    """Aggregate, by routine, a DrHook profile (statistics on the time spent
    in each routine) and write the result, ordered by self time, in *agg_file*.

    Parameters
    ----------
    prof_file : str
        DrHook profile (drhook.prof.* file, possibly the concatenation of the
        profiles of several MPI tasks).
    agg_file : str
        Output file.
    """
    # The data start at the first line beginning with '1' (first routine)
    lines = Path(prof_file).read_text(encoding='utf-8').splitlines()
    first_line = next((i for i, line in enumerate(lines) if re.match(r'^\s*1', line)), None)
    if first_line is None:
        print(f"No data found in the DrHook profile {prof_file}")
        return
    columns = {'time': ('<f4', ('mean',)),
               'self': ('<f4', ('mean', 'max', 'min', 'std', 'sum')),
               'total': ('<f4', ('mean', 'max', 'min', 'std', 'sum')),
               'calls': ('<i4', ('sum',)),
               'self_per_call': ('<f4', ('mean',)),
               'total_per_call': ('<f4', ('mean',)),
               'routine': ('U256', '')}
    arraynp = numpy.loadtxt(prof_file,
                            dtype=[(k, v[0]) for (k, v) in columns.items()],
                            converters={8: lambda s: s.split(b'@')[0].lstrip(b'*')},
                            skiprows=first_line, usecols=[1, 3, 4, 5, 6, 7, 8],
                            encoding='bytes')
    df = pandas.DataFrame(arraynp).groupby('routine').agg(
        **{k + '_' + agg: pandas.NamedAgg(column=k, aggfunc=agg)
           for k, (_, aggs) in columns.items() for agg in aggs
           if k != 'routine'}).sort_values('self_sum', ascending=False)
    df.index.name += ' ordered by self_sum'
    Path(agg_file).write_text(df.to_string(), encoding='utf-8')


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
