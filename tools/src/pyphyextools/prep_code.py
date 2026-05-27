#!/usr/bin/env python3
"""
Extract a tag/commit from the PHYEX repository, merge common and model-specific
code, apply pyfortool transformations, and optionally push the result.
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import glob

import pyfortool.scripting
from pyphyextools.testing.coding_norms import coding_norms

try:
    from pyphyex import __version__
except ImportError:
    __version__ = "UNKNOWN"


SEPARATOR = '_'


def parse_args(argv=None):
    """Parse and return command-line arguments."""
    if argv is None:
        argv = sys.argv[1:]

    # Split on '--': everything after is for pyfortool
    pyfortool_options = []
    if '--' in argv:
        idx = argv.index('--')
        pyfortool_options = argv[idx + 1:]
        argv = argv[:idx]

    parser = argparse.ArgumentParser(
        description="Prepare PHYEX code: extract, merge, transform, push.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
If -c is not provided, DIRECTORY must already contain files and directories
as if it was the result of a git checkout.
If -m is used, directory tree is modified, only relevant code is kept.
If none of --mnhExpand, --removeACC, --pyfortool_opts_env or PYFORTOOL_OPTIONS
is used, pyfortool is not called at all.
-s options are mandatory for -m, -D and -p options.
-p option is allowed only if -c and -m options are provided.
Everything after '--' is passed to pyfortool for source-to-source transformation.
""",
    )
    parser.add_argument('directory', nargs='?', help='directory containing the script result')
    parser.add_argument('-c', dest='checkout_point',
                        help='git object to checkout (commit or tags/TAG)')
    parser.add_argument('-m', dest='model', help='merge common code with code specific to MODEL')
    parser.add_argument('--mnhExpand', action='store_true', help='option passed to pyfortool')
    parser.add_argument('--removeACC', action='store_true', help='option passed to pyfortool')
    parser.add_argument('-s', dest='subs', action='append', default=[],
                        help='subdirectory or file (under src) to consider')
    parser.add_argument('-p', dest='push', action='store_true',
                        help='push the result as a new branch')
    parser.add_argument('--renameFf', action='store_true', help='rename .F90 into .f90')
    parser.add_argument('--ilooprm', action='store_true',
                        help='replace indexes in do loop (and mnh_expand) by :')
    parser.add_argument('--repo', help='repository URL (overrides env variables)')
    parser.add_argument('--noRaiseOnCodingNorms', action='store_true',
                        help='deactivate blocking check on coding norms')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('-v', dest='verbose', action='count', default=0,
                        help='increase verbosity: -v = INFO, -vv = DEBUG')
    parser.add_argument('--useParallelPyForTool', action='store_true',
                        help='use parallel pyfortool')
    parser.add_argument('--pyfortool_opts_env',
                        help='env var containing per-file pyfortool options')

    args = parser.parse_args(argv)
    args.pyfortool_options = pyfortool_options
    return args


def run(cmd, **kwargs):
    """Execute *cmd* via subprocess, forwarding extra kwargs.  Exit on failure."""
    logging.debug("Running: %s", ' '.join(cmd) if isinstance(cmd, list) else cmd)
    return subprocess.run(cmd, check=True, **kwargs)


def get_repository():
    """Build the PHYEX repository URL from ``PHYEXREPOprotocol`` and ``PHYEXREPOuser`` env vars."""
    protocol = os.environ.get('PHYEXREPOprotocol', '')
    user = os.environ.get('PHYEXREPOuser', '')
    if protocol == 'https':
        return f'https://github.com/{user}/PHYEX.git'
    if protocol == 'ssh':
        return f'git@github.com:{user}/PHYEX.git'
    return ''


def rename_Ff(directory, mv_func):
    """Rename every ``.F90`` file under *directory* to ``.f90`` using *mv_func*."""
    for root, _, files in os.walk(directory):
        for f in files:
            if f.endswith('.F90'):
                src = os.path.join(root, f)
                dst = os.path.join(root, f[:-4] + '.f90')
                mv_func(src, dst)


def merge_code(directory, model, subs, mv_func, rm_func):
    """Merge ``src/common`` with model-specific sources, then remove unneeded files."""
    src_dir_candidates = [f'src/{model}']
    if model == 'offline' and not os.path.isdir(os.path.join(directory, 'src', model)):
        src_dir_candidates.append('src/testprogs')
    src_dir = None
    for cand in src_dir_candidates:
        full = os.path.join(directory, cand)
        if os.path.isdir(full):
            src_dir = cand
            break
    if src_dir is None:
        raise RuntimeError(f"src/{model} directory does not exist")
    if not subs:
        raise RuntimeError("It is not possible to merge without -s options")

    logging.info("Merge common code and %s specific code", model)

    all_files = [f for f in os.listdir(directory) if not f.startswith('.')]

    for sub in subs:
        logging.debug("Merging %s directory/file", sub)
        sub_path = os.path.join(directory, sub)
        if os.path.exists(sub_path):
            raise RuntimeError(f"{sub} must not exist in the repository root, "
                              "this is a limitation of the script")

        common_sub = os.path.join(directory, 'src', 'common', sub)
        if os.path.exists(common_sub):
            mv_func(common_sub, sub_path)

        model_sub = os.path.join(directory, src_dir, sub)
        if os.path.exists(model_sub):
            if os.path.isfile(model_sub):
                mv_func(model_sub, sub_path)
            else:
                for root, _, files in os.walk(model_sub):
                    for f in files:
                        rel_path = os.path.relpath(os.path.join(root, f), model_sub)
                        dst_path = os.path.join(directory, sub, rel_path)
                        os.makedirs(os.path.dirname(dst_path), exist_ok=True)
                        mv_func(os.path.join(root, f), dst_path)
                try:
                    os.removedirs(model_sub)
                except OSError:
                    pass

    # Suppress unwanted files
    suppress_file = os.path.join(directory, src_dir, 'filesToSuppress.txt')
    if os.path.isfile(suppress_file):
        with open(suppress_file, encoding='UTF-8') as fh:
            for line in fh:
                filename = line.strip()
                full_path = os.path.join(directory, filename)
                if os.path.isfile(full_path):
                    rm_func(full_path)

    # Clean unrelevant files
    logging.info("Cleaning unrelevant files")
    for f in all_files:
        if f in ('.git', '.git/'):
            continue
        full_path = os.path.join(directory, f)
        if os.path.exists(full_path):
            logging.debug("Suppression of %s", f)
            rm_func(full_path)


def apply_ilooprm(directory):
    """Replace Fortran loop-index ranges (e.g. ``1:IKT``) with ``:`` in all source files."""
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        if not os.path.isdir(item_path):
            continue
        for root, _, files in os.walk(item_path):
            for f in files:
                fpath = os.path.join(root, f)
                fname = os.path.basename(f)
                if fname == 'minpack':
                    continue
                if fname.startswith('gradient_m'):
                    continue

                with open(fpath, encoding='UTF-8') as fh:
                    content = fh.read()

                # Protection for turb files
                if fname.startswith('turb'):
                    content = re.sub(
                        r'PLM(IIJB:IIJE,IKTB:IKTE) = ' + 
                            r'PZZ(IIJB:IIJE,IKTB\+IKL:IKTE\+IKL) - PZZ(IIJB:IIJE,IKTB:IKTE)',
                        r'PLM(IIJB:IIJE,IKTB : IKTE) = ' +
                            r'PZZ(IIJB:IIJE,IKTB+IKL:IKTE+IKL) - PZZ(IIJB:IIJE,IKTB : IKTE)',
                        content
                    )

                protections = {
                    'JK=IKTB:IKTE': 'transJKIKTB',
                    'JK=1:IKT': 'transIKT',
                    'JIJ=IIJB:IIJE': 'transJIJ',
                    'IKTB+1:IKTE': 'IKTB1IKTE',
                    'IKB+1:IKT': 'IKB1IKT',
                }
                for orig, repl in protections.items():
                    content = content.replace(orig, repl)

                # Apply transformations
                content = content.replace('1:IKT', ':')
                content = content.replace('IKTB:IKTE', ':')
                content = content.replace('IIJB:IIJE', ':')

                # Restore protections
                reversals = {v: k for k, v in protections.items()}
                for orig, repl in reversals.items():
                    content = content.replace(orig, repl)

                if fname.startswith('turb'):
                    content = content.replace('IKTB : IKTE', 'IKTB:IKTE')

                with open(fpath, 'w', encoding='UTF-8') as fh:
                    fh.write(content)


def prep_code(
    directory,
    checkout_point=None,
    model=None,
    mnh_expand=False,
    remove_acc=False,
    subs=None,
    push=False,
    rename_Ff_flag=False,
    ilooprm=False,
    repo=None,
    no_raise_on_coding_norms=False,
    verbose=0,
    use_parallel_pyfortool=False,
    pyfortool_opts_env=None,
    pyfortool_options=None,
):
    """Orchestrate checkout, renaming, merging, transformation, and push.

    Parameters
    ----------
    directory : str
        Directory containing the script result (or where to clone).
    checkout_point : str, optional
        Git object to checkout (commit hash or ``tags/TAG``).
    model : str, optional
        Model name for merging common + model-specific code.
    mnh_expand : bool
        Pass ``--mnhExpand`` to pyfortool.
    remove_acc : bool
        Pass ``--removeACC`` to pyfortool.
    subs : list of str, optional
        Subdirectories or files (under ``src``) to consider.
    push : bool
        Push the result as a new branch.
    rename_Ff_flag : bool
        Rename ``.F90`` files to ``.f90``.
    ilooprm : bool
        Replace loop-index ranges with ``:``.
    repo : str, optional
        Repository URL (overrides env variables).
    no_raise_on_coding_norms : bool
        Deactivate blocking check on coding norms.
    verbose : int
        Verbosity level (0-3).
    use_parallel_pyfortool : bool
        Use parallel version of pyfortool.
    pyfortool_opts_env : str, optional
        Environment variable name containing per-file pyfortool options.
    pyfortool_options : list of str, optional
        Extra options passed to pyfortool.
    """
    subs = subs or []
    pyfortool_options = pyfortool_options or []

    if verbose <= 0:
        logging.basicConfig(level=logging.WARNING, format="%(message)s")
    elif verbose == 1:
        logging.basicConfig(level=logging.INFO, format="%(message)s")
    else:
        logging.basicConfig(level=logging.DEBUG, format="%(message)s")

    # Get full directory path
    directory = os.path.abspath(directory)

    # Build pyfortool options list
    pyfortool_opts = []
    if mnh_expand:
        pyfortool_opts.append('--mnhExpand')
    if remove_acc:
        pyfortool_opts.append('--removeACC')
    pyfortool_opts.extend(pyfortool_options)

    # Determine repository URL
    repository = repo if repo else get_repository()

    # Build branch name for push
    branch = None
    if checkout_point and model and push:
        branch = f'{model}{SEPARATOR}{checkout_point}'

    # -------- WORKING DIRECTORY --------
    if not checkout_point:
        logging.info("No checkout point provided, we use the content of %s directory", directory)
        if not (os.path.isdir(os.path.join(directory, 'src')) or
                (os.path.isdir(os.path.join(directory, 'turb')) and
                 os.path.isdir(os.path.join(directory, 'micro')) and
                 os.path.isdir(os.path.join(directory, 'aux')))):
            raise RuntimeError(f"{directory} must be filled with files and directories "
                              "as if it was obtained through a checkout")

        def mv_func(src, dst):
            shutil.move(src, dst)
        def rm_func(path):
            if os.path.isdir(path):
                shutil.rmtree(path)
            else:
                os.remove(path)
    else:
        logging.info("Clone and checkout %s into %s directory", checkout_point, directory)
        if os.path.exists(directory):
            raise RuntimeError(f"{directory} already exists, suppress it before executing "
                              "the script (or remove the -c option)")
        if not repository:
            raise RuntimeError("A repository must be set (use -h option to get help)")

        run(['git', 'clone', repository, directory])
        if branch:
            result = subprocess.run(
                ['git', 'ls-remote', '--heads', 'origin', 'SR_GPU'],
                capture_output=True, text=True, check=True,
                cwd=directory
            )
            if result.stdout.strip():
                raise RuntimeError(f"{branch} branch already exists on remote")

        checkout_args = ['git', 'checkout']
        if branch:
            checkout_args.extend(['-b', branch])
        checkout_args.append(checkout_point)
        run(checkout_args, cwd=directory)

        def mv_func(src, dst):
            run(['git', 'mv', '-f', src, dst], cwd=directory)
        def rm_func(path):
            run(['git', 'rm', '-q', '-f', '-r', path], cwd=directory)

    # -------- RENAME .F90 into .f90 --------
    if rename_Ff_flag:
        rename_Ff(directory, mv_func)

    # -------- MERGE --------
    if model:
        merge_code(directory, model, subs, mv_func, rm_func)

    # -------- ILOOPRM --------
    if ilooprm:
        apply_ilooprm(directory)

    # -------- PyForTool --------
    if pyfortool_opts_env or pyfortool_opts:
        logging.info("Applying pyfortool")

        if model:
            reps = list(subs)
        else:
            reps = []
            for sub in subs:
                reps.extend(glob.glob(os.path.join('src', '*', sub)))

        extra_opts = []
        if pyfortool_opts_env:
            extra_opts = ['--optsByEnv', pyfortool_opts_env]

        if extra_opts or pyfortool_opts:
            orig_cwd = os.getcwd()
            os.chdir(directory)
            try:
                for rep in reps:
                    gitkeep = os.path.join(rep, '.gitkeep')
                    if os.path.isfile(gitkeep):
                        os.remove(gitkeep)

                if use_parallel_pyfortool:
                    cmd = (['pyfortool_parallel', '--wrapH'] + pyfortool_opts +
                           extra_opts + ['--nbPar', '8'])
                    logging.debug(' '.join(cmd))
                    pyfortool.scripting.mainParallel(cmd)
                else:
                    for rep in reps:
                        if not os.path.isdir(rep):
                            continue
                        for root, _, files in os.walk(rep):
                            for f in files:
                                if f.startswith('.'):
                                    continue
                                fpath = os.path.join(root, f)
                                ext = os.path.splitext(f)[1].lstrip('.')
                                if ext in ('fypp', 'yaml', 'in', 'cmake'):
                                    continue
                                if not ext:
                                    continue
                                cmd = (['pyfortool', '--wrapH'] + pyfortool_opts +
                                       extra_opts + [fpath])
                                logging.debug(' '.join(cmd))
                                pyfortool.scripting.main(cmd)
            finally:
                os.chdir(orig_cwd)

    # -------- Check coding conventions --------
    if not coding_norms(directory, verbose=True):
        if not no_raise_on_coding_norms:
            raise RuntimeError("Coding norms verification failed, see messages above")

    # -------- PUSH --------
    if branch:
        logging.info("commit and push")
        run(['git', 'add', '-A'], cwd=directory)
        commit_msg = (
            f"Version '{checkout_point}' of source code ready for "
            f"inclusion into {model} source tree"
        )
        run(['git', 'commit', '-m', commit_msg], cwd=directory)
        run(['git', 'push', '-u', 'origin', 'HEAD'], cwd=directory)

    logging.info("Finished!")


def main():
    """Parse command-line arguments and delegate to :func:`prep_code`."""
    args = parse_args()

    if not args.directory:
        raise RuntimeError("A directory must be provided on command line "
                           "(use -h option to get help)")

    prep_code(
        directory=args.directory,
        checkout_point=args.checkout_point,
        model=args.model,
        mnh_expand=args.mnhExpand,
        remove_acc=args.removeACC,
        subs=args.subs,
        push=args.push,
        rename_Ff_flag=args.renameFf,
        ilooprm=args.ilooprm,
        repo=args.repo,
        no_raise_on_coding_norms=args.noRaiseOnCodingNorms,
        verbose=args.verbose,
        use_parallel_pyfortool=args.useParallelPyForTool,
        pyfortool_opts_env=args.pyfortool_opts_env,
        pyfortool_options=args.pyfortool_options,
    )


if __name__ == '__main__':
    main()
