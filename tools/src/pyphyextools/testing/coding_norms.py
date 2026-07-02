#!/usr/bin/env python3
"""
Check coding norms in PHYEX Fortran sources:
- presence of IMPLICIT NONE
- INTENT attribute on all dummy arguments
- no (:) in subroutine calls
- no unused local variables
- matching interfaces for mode_ice4_pack / mode_ice4_stepping
- ONLY clause in USE statements
- no operations in CALL statements
"""

import argparse
import os
import subprocess
import sys
import tempfile

import pyfortool


def _args_from_file(filename):
    """Return sorted dummy-argument names for *filename* (empty list if missing)."""
    if not os.path.exists(filename):
        return []
    vars_ = [v for v in pyfortool.PYFT(filename).varList if v["arg"]]
    vars_ = sorted(vars_, key=lambda v: v["argorder"])
    return [v["n"] for v in vars_]


def coding_norms(sourcedir, verbose=False):
    """Run all coding-norm checks on *sourcedir*, return ``True`` if clean."""
    check = True
    results = []

    with tempfile.TemporaryDirectory() as tmpdir:
        tree_json = os.path.join(tmpdir, "tree.json")

        # Check interfaces of mode_ice4_pack and mode_ice4_stepping
        pack_file = os.path.join(sourcedir, "micro", "mode_ice4_pack.F90")
        step_file = os.path.join(sourcedir, "micro", "mode_ice4_stepping.F90")

        dummy_pack = _args_from_file(pack_file)
        dummy_pack = [v for v in dummy_pack if v not in ("D", "KSIZE")]
        dummy_step = _args_from_file(step_file)
        dummy_step = [v for v in dummy_step if v not in ("KMICRO",)]

        if dummy_pack != dummy_step:
            msg = "mode_ice4_pack and mode_ice4_stepping have different interfaces!"
            results.append(msg)
            check = False

        # Run pyfortool_parallel checks
        cmd = [
            "pyfortool_parallel",
            "--tree", sourcedir,
            "--descTree", tree_json,
            "--wrapH",
            "--enableCache",
            "--checkIMPLICIT", "Warn",
            "--checkINTENT", "Warn",
            "--checkONLY", "Warn",
            "--checkPHYEXUnusedLocalVar", "Warn",
            "--checkOpInCall", "Warn",
            "--checkEmptyParensInCall", "Warn",
        ]

        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
            out, err = proc.communicate()
            captured = (out + err).decode()
        if captured.strip():
            results.append(captured.rstrip())
            check = False

    if not check:
        if verbose:
            for r in results:
                print(r, end="")
        return False
    return True


def main():
    """Parse command-line arguments and delegate to :func:`coding_norms`."""
    parser = argparse.ArgumentParser(description="Check PHYEX Fortran coding norms")
    parser.add_argument("-v", action="store_true", help="print non-conforming messages")
    parser.add_argument("source", help="source directory to check")
    args = parser.parse_args()

    if not coding_norms(args.source, verbose=args.v):
        sys.exit(2)


if __name__ == "__main__":
    main()
