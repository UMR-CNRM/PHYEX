"""
Consistency checks of the PHYEX repository, also run by the GitHub workflow
.github/workflows/consistency.yaml:
  - version: the PHYEX version (src/phyex/VERSION) vs the release tags
  - pyfortool: same minimal pyfortool version everywhere
The checks are the scripts of .github/scripts (they only need the python
standard library, so that the workflow does not need to install the tools).
"""

import argparse
import os
import sys

from pyphyextools import REPO_DIR, run_command

CHECKS = {'version': 'check_version.py',
          'pyfortool': 'check_pyfortool_version.py'}


def check_consistency(checks=None, repo_dir=REPO_DIR):
    """
    Runs the consistency checks
    :param checks: list of checks to run (keys of CHECKS), all by default
    :param repo_dir: PHYEX repository to check (the one of the tools by default)
    :return: list of the checks that failed
    """
    failed = []
    for check in checks or CHECKS:
        script = os.path.join(repo_dir, '.github', 'scripts', CHECKS[check])
        print(f'### Check {check}')
        sys.stdout.flush()
        if run_command([sys.executable, script], cwd=repo_dir, check=False,
                       display='always', banner=False).returncode != 0:
            failed.append(check)
    return failed


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Consistency checks of the PHYEX repository')
    parser.add_argument('checks', nargs='*', metavar='CHECK',
                        help='checks to run among ' + ', '.join(CHECKS) + ' (default: all)')
    parser.add_argument('--repo-dir', default=REPO_DIR,
                        help='PHYEX repository to check (default: the repository of the tools, '
                             f'{REPO_DIR})')
    args = parser.parse_args()
    for check in args.checks:
        if check not in CHECKS:
            parser.error(f"invalid check '{check}' (choose from {', '.join(CHECKS)})")
    failed = check_consistency(args.checks, os.path.abspath(args.repo_dir))
    if failed:
        print('Failed checks: ' + ', '.join(failed))
        sys.exit(1)
    print('All checks passed')


if __name__ == '__main__':
    main()
