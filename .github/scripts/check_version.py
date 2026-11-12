#!/usr/bin/env python3

"""
Checks that the PHYEX version (src/phyex/VERSION) is the version of the last
release tag (vX.Y.Z) reachable from the checked commit (when the checked commit
is a release tag, this is this tag).
The repository must be cloned with its history and its tags.
Must be run from the repository root.
"""

import re
import subprocess
import sys

VERSION_FILE = 'src/phyex/VERSION'
RELEASE_TAG = re.compile(r'v([0-9]+(?:\.[0-9]+)*)')


def to_tuple(version):
    """'0.9' -> (0, 9, 0) (padded to compare versions of different lengths)"""
    values = [int(value) for value in version.split('.')]
    return tuple(values + [0] * (3 - len(values)))


def error(message, filename=VERSION_FILE):
    """Prints an error (GitHub annotation) and exits"""
    print(f'::error file={filename}::{message}')
    sys.exit(1)


def main():
    """Main function"""
    with open(VERSION_FILE, encoding='utf-8') as version_file:
        version = version_file.read().strip()
    if not re.fullmatch(r'[0-9]+(?:\.[0-9]+)*', version):
        error(f'Invalid version "{version}" in {VERSION_FILE} (x.y.z expected)')
    print(f'PHYEX version in {VERSION_FILE}: {version}')

    # Last release tag reachable from the checked commit
    tags = subprocess.run(['git', 'tag', '--merged', 'HEAD'], check=True,
                          capture_output=True, text=True).stdout.split()
    releases = [(to_tuple(match.group(1)), tag) for tag in tags
                if (match := RELEASE_TAG.fullmatch(tag))]
    if not releases:
        error('No release tag (vX.Y.Z) found, the repository must be cloned '
              'with its history and its tags')
    last_version, last_tag = max(releases)
    print(f'Last release tag: {last_tag}')
    if to_tuple(version) != last_version:
        error(f'The PHYEX version {version} is not the version of the last release tag {last_tag}')
    print('The PHYEX version is the version of the last release tag')


if __name__ == '__main__':
    main()
