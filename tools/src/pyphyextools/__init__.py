"""
pyphyextools python package
"""

import os

__version__ = "0.9.0"

# Check if current package is installed in editable mode
pyproject = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '..', '..', 'pyproject.toml')
if not os.path.exists(pyproject):
    raise RuntimeError('pyphyextools package must be installed in editable mode')
