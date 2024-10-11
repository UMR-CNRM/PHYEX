# PHYEX OFFLINE DOCUMENTATION

## About this section

This document is intended for persons who want to use PHYEX in an offline mode.
Some offline test programs are provided with the package and a library suitable for use with python is also available.

## Compilation

The build/with\_fcm directory in the master branch contains a build system.
This build system has two dependencies (installation is done automatically by the compilation script):

  - [fcm](https://metomi.github.io/fcm/doc/user_guide/)
  - [fiat](https://github.com/ecmwf-ifs/fiat)

The script build/with\_fcm/make\_fcm.sh uses configuration files and build the library and test programs.
These executables can be found in the build/bin subdirectory in the architecture specific directory arch\_\<architecture name\>.

Some configuration files are stored in build/with\_fcm/arch but other can be maintained by the end users in their
${HOME}/.phyex/fcm\_arch.

Some more details on the build system can be found in [build/with\_fcm/README.md file](../build/with_fcm/README.md).

### Compilation directly in the repository without execution (or manual execution)

When on a master commit, the build/with\_fcm/make\_fcm.sh script can be used to compile the offline tools.

### Compilation and execution

When on a master commit, the tools/check\_commit\_testprogs.sh script can be used to compile and execute the testprogs (offline tests).
The check\_commit\_testprogs.sh script uses the PHYEX source code:

  - of a specific commit on the master branch available on a remote repository
  - or, the last commit of a offline\_\<commit\_hash\> branch available on a remote repository
  - or, the content of a local repository.

In the latter case, it can be interesting to clone the PHYEX repository twice.
A first one to have the build tools on the master branch, and a second one to checkout the source code version to use.
This solution is especially useful when working on a offline\_\<commit\_hash\> branch (because these branches does not
contain the build tools).

Something like this can be used:

- cd $HOME; git clone \<PHYEX url\> PHYEXtools
- cd PHYEXtools; git checkout master
- cd $HOME; git clone \<PHYEX url\> PHYEX
- cd PHYEX; git checkout arome\_\<commit\_hash\>; source code moddifications...
- . PHYEXtools/tools/env.sh; check\_commit\_testprogs.sh $HOME/PHYEX REF

The last step will create a directory (in $HOME/TESTPROGS) with a copy of your source code and the build system, builds the testprogs and executes them.

## Test program

### Data generation

The branch testprogs\_data contains modified source code for the AROME model to enable the generation of data samples for the turb, shallow, rain\_ice and ice\_adjust testprogs.
The branch testprogs\_data2 contains modified source code for the AROME model to enable the generation of data samples for the rain\_ice\_old testprog.
The branch testprogs\_data3 contains modified source code for the AROME model to enable the generation of data samples for the lima and lima\_adjust testprogs.
Using these branches, in the drivers of the different parametrisations (aro\_\* files), output can be enable for the AROME model.
Running the AROME model with these modifications outputs files in the running directory.
This must be done once by parametrisation (note that the check\_commit\_ial.sh script can be used to execute an AROME simulation).

These files should be renamed with the following command:
i=0; for file in ????_??_????????.dat; do mv $file `printf %08d $i`.dat; i=$((i+1)); done

### Usage directly with the testprogs executables

The different main\_\*.exe programs obtained by the compilation can be run. Each of these executables is expecting the presence of a 'data' directory in their working directory containing the different files.

### Usage through the check\_commit\_testprogs.sh script

As described in [COMPILATION](#compilation).

## Python bindings

The package includes a python binding to call the main subroutines of PHYEX.
Using the default gfortran compiler, the needed compilation is achieved with:

- . PHYEX/tools/env.sh
- cd PHYEX/build/with\_fcm
- ./make\_fcm.sh

This builds a shared library and writes a python wrapper.
This script needs the ctypesForFortran package (available on PyPI).
And it can be used like in this example for ice\_adjust (the sourcing of PHYEX/tools/env.sh is also needed at this step):

```
#!/usr/bin/env python3

"""
Offline python binding example, valid with a June 2024 version of PHYEX
"""

import tempfile
from pyphyex import PYICE_ADJUST, PYINI_PHYEX, close
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
                      PCF_MF, PRC_MF, PRI_MF, PRV, PRC, PRVS, PRCS, PTH, PTHS,
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
