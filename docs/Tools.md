# TOOLS DOCUMENTATION

## ABOUT THIS DOCUMENT

This document is intended for persons who want to use the prep\_code.sh or the check\_commit\_\*.sh scripts.

This document is written using the markdown language. With pandoc, it can be converted to HTML (pandoc -s \<filename\>.md -o \<filename\>.html) or PDF (pandoc -s \<filename\>.md -o \<filename\>.pdf).

## INSTALLATION, PATH...

Installation is covered in the [INSTALL documentation](../tools/INSTALL.md).

Environment variables can be set with:

```
. <git repository>/tools/env.sh
```

## TOOLS

### check\_commit\_ial.sh

The check\_commit\_ial script compiles, executes IAL test cases and compare the results against a reference simulations.

Script options can be displayed with the -h option.

Before being usable, the AROME model must be installed following the [INSTALL\_pack\_ial documentation](../tools/INSTALL_pack\_ial.md).

### check\_commit\_mesonh.sh

The check\_commit\_mesonhsh script compiles, runs a test case of the Meso-NH model and compares the results against a reference simulation.

Script options can be displayed with the -h option.

Before being usable, the mesonh model must be installed following the [INSTALL\_pack\_mesonh documentation](../tools/INSTALL_pack_mesonh.md).

For check\_commit\_mesonh.sh the following environment variables can be set:

  - MNHPACK: directory in which MNH pack will be created (default is $HOME/MesoNH/PHYEX)
  - REFDIR: directory in which reference pack can be found (default is the pack directory near the check\_commit\_mesonh.sh file)
  - TARGZDIR: directory in which the tar.gz file can be found (default is the pack directory near the check\_commit\_mesonh.sh file)

### check\_commit\_testprogs.sh

The check\_commit\_testprogs.sh script runs offline simulations in the directory given
by the environment variable TESTPROGSDIR ($HOME/TESTPROGS will be used if the variable is not set).
This directory must exist.

Script options can be displayed with the -h option.

To be usable the check\_commit\_testprogs.sh script needs input data. The generation and installation of these data are described in the [INSTALL\_testprogs documentation](../tools/INSTALL_testprogs.md).

The goal of the script is to compare outputs between two simulations (to check if bit-reproducibilty is achieved or not).
A reference simulation must be performed and save. This reference simulation is run the same way as the
test experiment but cannot be compared to something else:
check\_commit\_testprogs.sh -c -r <reference_commit>

If this reference simulation must become the 'absolute' reference (used when invoking the check\_commit\_testprogs.sh
script with the 'REF' argument), the reference simulation directory (under $TESTPROGSDIR) must be renamed 'ref'.

### prep\_code.sh

This script is used by the different check\_commit\_\* scripts and can be used directly to pre-process the source code.
 The installation is described in the [INSTALL\_mnh\_expand documentation](../tools/INSTALL_mnh_expand.md)

### others

Other scripts are:

  - comp\_DDH.py: compare DDH outputs (used by check\_commit\_ial.sh)
  - compare.py: compare MESO-NH outputs (used by check\_commit\_mesonh.sh)
  - correct\_indent.py: correct source code indentation in mnh\_expand blocs
  - diffNODE.001\_01: compare NODE.0001\_01 output files
  - verify\_mnh\_expand.py: check syntax in mnh\_expand blocs
