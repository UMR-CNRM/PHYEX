"""
pyphyex: python interface to the PHYEX atmospheric physics package

Modules:
  - lowlevel: python binding of some PHYEX routines (PYINI_PHYEX, PYICE_ADJUST,
    PYRAIN_ICE, PYTURB...); this module, generated during the build, loads the
    PHYEX shared library when imported
  - phyex1d: 1D model using PHYEX
"""

try:
    # File written during the build (see cmake/phyex_python_wheel.cmake)
    from ._version import __version__
except ImportError:  # source tree
    __version__ = 'unknown'
