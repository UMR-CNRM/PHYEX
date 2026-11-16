# PHYEX OFFLINE DOCUMENTATION

## About this section

This document is intended for persons who want to use PHYEX in an offline mode.
Some offline test programs are provided with the package and a python package (python binding and 1D model) is also available.

## Compilation

The library and the test programs are compiled with CMake (version 3.18 or later), from the main CMakeLists.txt file:

```
. PHYEX/arch/arch-<architecture>.env
cmake -S PHYEX -B build -DCMAKE_TOOLCHAIN_FILE=PHYEX/arch/arch-<architecture>.ecbuild -DENABLE_PHYEX_BUILD_PROGS=ON
cmake --build build -j
```

The PHYEX source code needs the [fiat](https://github.com/ecmwf-ifs/fiat) package and the build system needs [ecbuild](https://github.com/ecmwf/ecbuild).
If they are not found (installed versions can be given with the usual CMake mechanisms, e.g. CMAKE\_PREFIX\_PATH), they are downloaded and built automatically (an internet connection is then needed).
Their versions and repositories can be changed with the PHYEX\_ECBUILD\_VERSION, PHYEX\_ECBUILD\_REPO, PHYEX\_FIAT\_VERSION and PHYEX\_FIAT\_REPO variables (environment or -D options).

The arch directory contains the architecture files: arch-\<architecture\>.ecbuild is a CMake toolchain file (compilers and compilation options) and arch-\<architecture\>.env is an environment file to source before compiling and executing (modules to load...).

The main options (-D\<option\>=\<value\>) are:

  - ENABLE\_DOUBLE\_PRECISION (ON by default) and ENABLE\_SINGLE\_PRECISION (OFF by default): precisions to compile (the libraries are named libphyex\_dp and libphyex\_sp)
  - ENABLE\_PHYEX\_BUILD\_PROGS (OFF by default): compilation of the test programs (main\_\<test\>\_\<precision\>.exe)
  - ENABLE\_OMP, ENABLE\_ACC: OpenMP and OpenACC support
  - PHYEX\_TRANSFORMATION (also read from the environment variable of the same name): transformation of the source code at build time by [pyfortool](https://github.com/UMR-CNRM/pyfortool), among IAL (transformations applied for the AROME model: mnh\_expand directives expanded in DO loops, shuman functions transformed...), IAL\_NOSHUGRAD and IAL\_NOEXPAND (empty or NONE, the default, means no transformation). If pyfortool is not found, or is too old, it is installed automatically in a python virtual environment in the build directory
  - ENABLE\_CHECK\_NORMS (ON by default): check of the coding norms of the transformed source code
  - ENABLE\_PYTHON\_WHEEL (OFF by default): generation of the python binding (used to build the python package, see [below](#python-package))

The src/phyex directory can also be compiled alone, using its own CMakeLists.txt file (for instance when included in a bigger project providing ecbuild and fiat).

## Test programs

### Usage through the phyex-check\_testprogs tool

The phyex-check\_testprogs tool (see the [Tools documentation](./Tools.md)) compiles the test programs (with the CMake build system), executes them and compares the results with those of a reference version. The source code can be a commit of a remote repository or the content of a local repository (with its uncommitted modifications), for instance:

```
phyex-check_testprogs -p -c -r -C --computeRefIfNeeded $HOME/PHYEX REF
```

This command creates a directory (in $TESTPROGSDIR, $HOME/TESTPROGS by default) with a copy of the source code, builds and executes the testprogs, then compares the results with the reference results (REF: the reference commits given in the 'testing' entry of the src/phyex/testprogs\_version.json file, computed if needed).

### Data

The test programs need input data, stored in ${PHYEXCONF}/testprogs\_data (the PHYEXCONF environment variable defaults to ${HOME}/.phyex).
If this directory is empty, a reduced dataset is automatically downloaded and installed when phyex-check\_testprogs is used.

The data can also be generated:

  - The branch testprogs\_data contains modified source code for the AROME model to enable the generation of data samples for the turb, shallow, rain\_ice and ice\_adjust testprogs.
  - The branch testprogs\_data2 contains modified source code for the AROME model to enable the generation of data samples for the rain\_ice\_old testprog.
  - The branch testprogs\_data3 contains modified source code for the AROME model to enable the generation of data samples for the lima and lima\_adjust testprogs.

Using these branches, in the drivers of the different parametrisations (aro\_\* files), output can be enable for the AROME model.
Running the AROME model with these modifications outputs files in the running directory.
This must be done once by parametrisation (note that the phyex-check\_ial tool can be used to execute an AROME simulation).

These files should be renamed with the following command:
i=0; for file in ????\_??\_????????.dat; do mv $file `printf %08d $i`.dat; i=$((i+1)); done

and put in the corresponding directories under ${PHYEXCONF}/testprogs\_data.

### Usage directly with the testprogs executables

The different main\_\*.exe programs obtained by the compilation can be run. Each of these executables is expecting the presence of a 'data' directory in their working directory containing the different files.

## Python package

The pyphyex python package contains:

- pyphyex.lowlevel: python binding of the main subroutines of PHYEX (PYINI\_PHYEX, PYICE\_ADJUST,
  PYRAIN\_ICE, PYTURB, PYSHALLOW\_MF...). This module loads the PHYEX shared library when imported.
- pyphyex.phyex1d: a 1D model using PHYEX, also available through the phyex1d-\* commands
  (see the [phyex1d documentation](./phyex1d.md)).

The package is compiled and installed with pip (python >= 3.10, a Fortran compiler and CMake >= 3.18 are needed) from the root of the PHYEX repository:

```
pip install .          # or, to build a wheel: pip wheel . -w dist
```

CMake options (e.g. toolchain file) can be given with
'-C cmake.define.CMAKE\_TOOLCHAIN\_FILE=\<arch\>.ecbuild' or with the CMAKE\_ARGS environment variable.

### Development installation

The package can be installed in editable mode (pip install -e .):

- the python modules (src/pyphyex, among them phyex1d) are used directly from the source tree, while
  the files produced by CMake (pyphyex/lowlevel.py, pyphyex/\_version.py, pyphyex/lib/) are installed in
  site-packages. A modification of the Fortran code (or of the bound routines) needs a new installation.
- by default the CMake build is done in a temporary directory (full rebuild at each installation). With
  '-C build-dir=build/{wheel\_tag}' the build tree is kept and a new installation is incremental.
- with '-C editable.rebuild=true' (experimental), the CMake build is run at each import of the package, so
  Fortran modifications are taken into account automatically. It needs build-dir and --no-build-isolation
  (the build requirements, scikit-build-core and pyfortool, must then be installed in the environment).
  Each import is slower and needs the compilers.

For instance:
```
pip install scikit-build-core pyfortool
pip install --no-build-isolation -e . -C build-dir=build/{wheel_tag} -C editable.rebuild=true
```
or, without editable.rebuild, run the same command again after a Fortran modification.

### Python binding

The low-level module can be used like in this example for ice\_adjust:

```
#!/usr/bin/env python3

"""
Offline python binding example
"""

import tempfile
from pyphyex.lowlevel import PYICE_ADJUST, PYINI_PHYEX, close
import numpy
from matplotlib import pyplot as plt
import f90nml

#Config
NIT, NKT = 150, 100
VSIGQSAT = 0.02
PSIGS = numpy.ones((NKT, NIT)) * 0.0005
PPABST = numpy.ones((NKT, NIT)) * 101325.
PRC = numpy.zeros((NKT, NIT))
PRI = numpy.zeros((NKT, NIT))
PRV = numpy.linspace(0.003, 0.01, num=NKT)
PTH = numpy.linspace(280., 295., num=NIT)
PTH, PRV = numpy.meshgrid(PTH, PRV)

#Other values
PTSTEP = 60.
PROGRAM = 'AROME'
HBUNAME = 'DEPO'
NKL, NVEXT, KRR = 1, 0, 6
PSIGQSAT = numpy.ones((NIT,)) * VSIGQSAT
PRHODJ = numpy.ones((NKT, NIT))
PEXNREF = numpy.ones((NKT, NIT))
PEXN = numpy.ones((NKT, NIT))
PRHODREF = numpy.ones((NKT, NIT))
LMFCONV = False
PMFCONV = numpy.zeros((NKT, NIT))
PZZ = numpy.zeros((NKT, NIT))
PCF_MF = numpy.zeros((NKT, NIT))
PRC_MF = numpy.zeros((NKT, NIT))
PRI_MF = numpy.zeros((NKT, NIT))
PWEIGHT_MF_CLOUD = numpy.zeros((NKT, NIT))
OCOMPUTE_SRC = True
PRR = numpy.zeros((NKT, NIT))
PRS = numpy.zeros((NKT, NIT))
PRG = numpy.zeros((NKT, NIT))
KBUDGETS = 0

with tempfile.NamedTemporaryFile() as f:
    nml = f90nml.read(f.name)
    nml['NAM_NEBn'] = {}
    nml['NAM_NEBn']['CCONDENS'] = 'CB02'
    nml['NAM_NEBn']['LSUBG_COND'] = True
    nml['NAM_NEBn']['LSIGMAS'] = True
    nml.write(f.name, force=True)
    PYINI_PHYEX(PROGRAM, 33, f.name, False, 20, 0, 1, PTSTEP, 20., 'ICE3', 'EDKF', 'TKEL',
                LDCHANGEMODEL=True, LDDEFAULTVAL=True, LDREADNAM=True, LDCHECK=True,
                KPRINT=0, LDINIT=True)
PRVS = PRV / PTSTEP
PRCS = PRC / PTSTEP
PRIS = PRI / PTSTEP
PTHS = PTH / PTSTEP
result = PYICE_ADJUST(NIT, NKT, NKL, NVEXT, KRR, HBUNAME, PTSTEP, PSIGQSAT, PRHODJ,
                      PEXNREF, PRHODREF, PSIGS, LMFCONV, PMFCONV, PPABST, PZZ, PEXN,
                      PCF_MF, PRC_MF, PRI_MF, PWEIGHT_MF_CLOUD, PRV, PRC, PRVS, PRCS, PTH, PTHS,
                      OCOMPUTE_SRC, PRR, PRI, PRIS, PRS, PRG, KBUDGETS)
(_, _, _, _, _, PRVS, PRCS, PTHS, PSRCS, PCLDFR, PRIS,
 POUT_RV, POUT_RC, POUT_RI, POUT_TH, PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF) = result
PRC = PRCS / PTSTEP

close()

fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True)

cf = ax[0].contourf(PRV, PTH, PRC, levels=numpy.linspace(0., 5.6E-7, 15))
fig.colorbar(cf, ax=ax[0])
ax[0].set_title('Cloud mixing ratio (kg/kg)')
ax[0].set_ylabel('Potential temperature (K)')
ax[0].set_xlabel('Total water mixing ratio (kg/kg)')

cf = ax[1].contourf(PRV, PTH, PCLDFR, levels=numpy.linspace(0., 1., 15))
fig.colorbar(cf, ax=ax[1])
ax[1].set_title('Cloud fraction (0-1)')
ax[1].set_xlabel('Total water mixing ratio (kg/kg)')

plt.show()
```
