#!/usr/bin/env python3

"""Check the reference hashes stored in a json file (--refHashFile option of the
phyex-check_testprogs and phyex-check_phyex1d tools) against the reference results present
on disk (computed during the workflow run, for instance by the detailed comparison done after
a failure).

A stored hash which differs from the hash of a freshly computed reference means that something
changing the results (compiler, system libraries...) is not taken into account in the cache key
(see .github/actions/refhash/action.yml).

Usage: check_refhash.py testprogs REFHASHFILE ARCH EXTRAPOLATION
       check_refhash.py phyex1d REFHASHFILE ARCH
"""

import json
import os
import sys

from pyphyextools.testing.check_phyex1d import RESULT_FILE, case_tag
from pyphyextools.testing.check_testprogs import CheckCommitTestprogs
from pyphyextools.testing.compare import HASH_FORMAT, hash_phyex1d, hash_testprogs


def result_file(tool, arch, extrap, commit, case):
    """Result file of a reference simulation (same paths as in the tools)"""
    if tool == 'testprogs':
        testdir = os.environ.get('TESTPROGSDIR', os.path.expanduser('~/TESTPROGS'))
        tag = CheckCommitTestprogs._extrapolation_configs()[0][extrap]  # pylint: disable=protected-access
        return os.path.join(testdir, commit, 'tests', f'arch_{arch}', f'{case}{tag}',
                            'Output_run')
    testdir = os.environ.get('PHYEX1DDIR', os.path.expanduser('~/PHYEX1D'))
    return os.path.join(testdir, commit, 'tests', f'arch_{arch}', case_tag(case), RESULT_FILE)


def main():
    """Compare the stored hashes with the hashes of the reference results present on disk"""
    tool, refhashfile, arch = sys.argv[1:4]
    extrap = int(sys.argv[4]) if tool == 'testprogs' else None
    hash_function = {'testprogs': hash_testprogs, 'phyex1d': hash_phyex1d}[tool]
    if not os.path.isfile(refhashfile):
        print(f"No reference hash file ({refhashfile})")
        return 0
    with open(refhashfile, encoding='utf-8') as f:
        content = json.load(f)
    if content.get('format') != HASH_FORMAT:
        print(f"Hash format of {refhashfile} is not the current one, nothing to check")
        return 0
    status = 0
    for case, entry in sorted(content.get('cases', {}).items()):
        filename = result_file(tool, arch, extrap, entry['commit'], case)
        if not os.path.isfile(filename):
            print(f"{case}: no reference result on disk for commit {entry['commit']}")
        elif hash_function(filename) == entry['hash']:
            print(f"{case}: the stored hash is the hash of the reference result "
                  f"(commit {entry['commit']})")
        else:
            print(f"ERROR {case}: the stored hash differs from the hash of the reference result "
                  f"(commit {entry['commit']}, {filename})")
            status = 1
    return status


if __name__ == '__main__':
    sys.exit(main())
