# PHYEX OFFLINE DOCUMENTATION

## About this section

This document is intended for persons who want to use PHYEX in an offline mode.
Some offline test programs are provided with the package and a library suitable for use with python is also available.

## Compilation

The build/with\_fcm directory contains a build system.
This build system has two dependencies (installation is done automatically by the compilation script):

  - [fcm](https://metomi.github.io/fcm/doc/user_guide/)
  - [fiat](https://github.com/ecmwf-ifs/fiat)

The script build/with\_fcm/make\_fcm.sh uses a configuration file and build the library and test programs.
They can be found in the build/bin subdirectory in the architecture specific directory arch\_\<architecture name\>.

Some more details on the build system can be found in [build/with\_fcm/README.md file](../build/with_fcm/README.md).

### Compilation in the repository

When on a master commit, the build/with\_fcm/make\_fcm.sh script can be used to compile the offline tools.

### Compilation of a testprogs specific commit

When on a master commit, the tools/check\_commit\_testprogs.sh script can be used to compile a testprogs\_\<commit\_hash\> commit already available in the PHYEX repository.

### Modification of a testprogs specific commit

The idea is to use two clones of PHYEX: one to provide the building tools, the another one to provide the source code.
Something like this can be used:

- cd $HOME; git clone \<PHYEX url\> PHYEXtools
- cd PHYEXtools; git checkout master
- cd $HOME; git clone \<PHYEX url\> PHYEX
- cd PHYEX; git checkout arome\_\<commit\_hash\>; source code moddifications...
- . PHYEXtools/tools/env.sh; check\_commit\_testprogs.sh $HOME/PHYEX REF

The last step will creates a directory (in $HOME/TESTPROGS) with a copy of your source code and the build system, builds the testprogs and executes them.

## Test program

### Data generation

The branch testprogs\_data contains modified source code for the AROME model to enable the generation of data samples.
Using this branch, in the drivers of the different parametrisations (aro\_\* files), output can be enable for the AROME model.
Running the AROME model with these modifications outputs files in the running directory.
This must be done once by parametrisation (note that the check\_commit\_ial.sh script can be used to execute an AROME simulation).

These files should be renamed with the following command:
i=0; for file in ????_??_????????.dat; do mv $file `printf %08d $i`.dat; i=$((i+1)); done

### Usage directly with the testprogs executables

The different main\_\*.exe programs obtained by the compilation can be run. Each of these executables is expecting the presence of a 'data' directory in their working directory containing the different files.

### Usage through the check\_commit\_testprogs.sh script

As described in [COMPILATION](#compilation).

## Python bindings

**TODO** This section (and code) must be written. Key ideas are:

  - ctypesforfortran
  - example
