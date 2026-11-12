#!/usr/bin/env python3

"""
Checks that the minimal pyfortool version is the same in:
  - pyproject.toml (build requirements of the python wheel)
  - tools/pyproject.toml (dependencies of the python tools)
  - src/phyex/cmake/phyex_pyfortool.cmake (pyfortool installed by CMake)
Must be run from the repository root (python >= 3.11).
"""

import re
import sys
import tomllib

MIN_VERSION = re.compile(r'>=\s*([0-9]+(?:\.[0-9]+)*)')


def find_spec(requirements, filename, section):
    """Returns the pyfortool requirement of a list of requirements"""
    for requirement in requirements:
        if re.match(r'pyfortool\b', requirement.strip(), re.IGNORECASE):
            return requirement
    error(f'pyfortool not found in {section}', filename)
    return None


def error(message, filename):
    """Prints an error (GitHub annotation) and exits"""
    print(f'::error file={filename}::{message}')
    sys.exit(1)


def main():
    """Main function"""
    specs = {}

    filename = 'pyproject.toml'
    with open(filename, 'rb') as toml_file:
        requirements = tomllib.load(toml_file)['build-system']['requires']
    specs[filename] = find_spec(requirements, filename, 'build-system.requires')

    filename = 'tools/pyproject.toml'
    with open(filename, 'rb') as toml_file:
        requirements = tomllib.load(toml_file)['project']['dependencies']
    specs[filename] = find_spec(requirements, filename, 'project.dependencies')

    filename = 'src/phyex/cmake/phyex_pyfortool.cmake'
    with open(filename, encoding='utf-8') as cmake_file:
        match = re.search(r'set\(\s*_default\s+"(pyfortool[^"]*)"\s*\)', cmake_file.read())
    if not match:
        error('Default value of PHYEX_PYFORTOOL_SPEC not found', filename)
    specs[filename] = match.group(1)

    versions = {}
    for filename, spec in specs.items():
        match = MIN_VERSION.search(spec)
        if not match:
            error(f'No minimal version (>=) in the pyfortool requirement "{spec}"', filename)
        versions[filename] = match.group(1)
        print(f'{filename}: {spec}')

    if len(set(versions.values())) != 1:
        details = ', '.join(f'{version} in {filename}' for filename, version in versions.items())
        for filename in versions:
            print(f'::error file={filename}::Minimal pyfortool versions differ: {details}')
        sys.exit(1)
    print(f'The minimal pyfortool version is the same everywhere: {versions[filename]}')


if __name__ == '__main__':
    main()
