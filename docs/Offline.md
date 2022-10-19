# OFFLINE documentation

## ABOUT THIS DOCUMENT

This document is intended for persons who want to use PHYEX in an offline mode.
Some offline test programs are provided with the package and a library suitable for use with python is also available.

This document is written using the markdown language. With pandoc, it can be converted to HTML (pandoc -s \<filename\>.md -o \<filename\>.html) or PDF (pandoc -s \<filename\>.md -o \<filename\>.pdf).

## COMPILATION

The build/with\_fcm directory contains a build system.
This build system has two dependencies (installation is done automatically by the compilation script):

  - [fcm](https://metomi.github.io/fcm/doc/user_guide/)
  - [fiat](https://github.com/ecmwf-ifs/fiat)

The script build/with\_fcm/make\_fcm.sh uses a configuration file and build the library and test programs.
They can be found in the build/bin subdirectory in the architecture specific directory arch\_\<architecture name\>.

Some more details can be found in [build/with\_fcm/README.md file](../build/with_fcm/README.md).

## TEST PROGRAM

### Data generation

The branch testprogs\_data contains modified source code for the AROME model to enable the generation of data samples.
Using this branch, in the drivers of the different parametrisations (aro\_\* files), output can be enable for the AROME model.
Running the AROME model with these modifications outputs files in the running directory.
This must be done once by parametrisation (note that the check\_commit\_ial.sh script can be used to execute an AROME simulation).

These files should be renamed with the following command:
i=0; for file in ????_??_????????.dat; do mv $file `printf %08d $i`.dat; i=$((i+1)); done

### Usage directly with the testprogs executables

The different main\_\*.exe programs obtained by the compilation can be run. Each of these executables is expecting the presence of a 'data' directory in their working directory containing the different files.

## PYTHON BINDING

**TODO** This section must be written. Key ideas are:

  - ctypesforfortran
  - example
