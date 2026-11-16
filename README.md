# PHYEX
PHYsique EXternalisée

[![tests](https://github.com/UMR-CNRM/PHYEX/actions/workflows/tests.yaml/badge.svg)](https://github.com/UMR-CNRM/PHYEX/actions/workflows/tests.yaml)
[![consistency](https://github.com/UMR-CNRM/PHYEX/actions/workflows/consistency.yaml/badge.svg)](https://github.com/UMR-CNRM/PHYEX/actions/workflows/consistency.yaml)
[![Doxygen](https://github.com/UMR-CNRM/PHYEX/actions/workflows/doxygen.yaml/badge.svg?branch=master)](https://github.com/UMR-CNRM/PHYEX/actions/workflows/doxygen.yaml)
[![AROME branch](https://github.com/UMR-CNRM/PHYEX/actions/workflows/arome_branch.yaml/badge.svg?branch=master)](https://github.com/UMR-CNRM/PHYEX/actions/workflows/arome_branch.yaml)

PHYEX is the atmospheric physics package (turbulence, shallow convection and microphysics
schemes) shared by the Meso-NH, AROME (and HARMONIE-AROME) and LMDZ models.

The documentation can be found in the [docs directory](./docs/PHYEX.md). Several
presentations were done, the materials can be found on the
[wiki](https://github.com/UMR-CNRM/PHYEX/wiki).

## FORTRAN users: compiling the library

The library (and, optionally, the offline test programs) is compiled with CMake directly from
the main CMakeLists.txt file.

Requirements:
  - CMake >= 3.18
  - FORTRAN and C compilers
  - git and an internet connection if [ecbuild](https://github.com/ecmwf/ecbuild) and
    [fiat](https://github.com/ecmwf-ifs/fiat) are not already installed (they are then
    downloaded and built automatically)
  - python >= 3.8 (with the venv module) only if the source code is transformed at build time
    (PHYEX\_TRANSFORMATION option, the pyfortool package is then installed automatically if
    needed)

Example:
```
git clone https://github.com/UMR-CNRM/PHYEX.git
cmake -S PHYEX -B build -DCMAKE_TOOLCHAIN_FILE=PHYEX/arch/arch-gnu.ecbuild
cmake --build build -j
```

The build options (precision, test programs, transformations...) are described in the
[Offline documentation](./docs/Offline.md). The way to plug the library in a model is
described in the [Plugging documentation](./docs/Plugging.md).

## Python users: the pyphyex package

The pyphyex package provides a python binding of the main PHYEX routines (pyphyex.lowlevel)
and a 1D model (pyphyex.phyex1d, with the phyex1d-\* commands). It can be installed from this
repository (same requirements as for the FORTRAN library, CMake >= 3.18 and the compilers,
python >= 3.10):
```
pip install ./PHYEX
```

See the [Offline documentation](./docs/Offline.md) for the python binding and the
[phyex1d documentation](./docs/phyex1d.md) for the 1D model (for instance:
`phyex1d-cases_install download; phyex1d-exp_run ARMCU/REF --exp dt=60`).

## Developers

The contribution workflows are described in the [Developer documentation](./docs/Developer.md)
and the coding norms in the [CodingNorms documentation](./docs/CodingNorms.md).

The tools (source code preparation and test scripts, phyex-check\_\* commands) are provided as a
python package to install, in editable mode, from the repository (python >= 3.10 recommended):
```
git clone https://github.com/UMR-CNRM/PHYEX.git
python3 -m venv phyex.env
. phyex.env/bin/activate
pip install -e PHYEX/tools
```

For instance, the offline test programs (testprogs) can be compiled, executed and compared to
the reference results with the source code of a local repository (input data are downloaded
automatically):
```
phyex-check_testprogs -p -c -r -C --computeRefIfNeeded $PWD/PHYEX REF
```

All the checks needed before a contribution (coding norms, consistency, documentation, python
wheel and the test cases of the different models) can be run with `phyex-check_commit`.
The tools are described in the [Tools documentation](./docs/Tools.md).

## Main references
  - PHYEX:
    - [code](https://github.com/UMR-CNRM/PHYEX)
  - MesoNH:
    - [reference](https://doi.org/10.5194/gmd-11-1929-2018)
    - [documentation](https://mesonh.cnrs.fr/)
    - [code](https://src.koda.cnrs.fr/mesonh/mesonh-code)
  - AROME:
    - [reference](https://doi.org/10.1175/2010MWR3425.1)
  - Common format:
    - [documentation and cases](https://github.com/GdR-DEPHY/DEPHY-SCM)

## License

PHYEX is distributed under the CeCILL-C license (see the [LICENSE](./LICENSE) file).
