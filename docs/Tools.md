# PHYEX TOOLS DOCUMENATTION

## About this section

This document is intended for persons who want to use the prep\_code.sh or the check\_commit\_\*.sh scripts.

## Installation, path...

Installation is covered in the [tools/INSTALL.md file](../tools/INSTALL.md).

Environment variables can be set with:

```
. <git repository>/tools/env.sh
```

## Tools

### check\_commit\_ial.sh

The check\_commit\_ial script compiles, executes IAL test cases and compare the results against a reference simulations.

Script options can be displayed with the -h option.

Before being usable, some packages must be installed following the [tools/INSTALL.md file](../tools/INSTALL.md).

### check\_commit\_mesonh.sh

The check\_commit\_mesonh.sh script compiles, runs a test case of the Meso-NH model and compares the results against a reference simulation.

Script options can be displayed with the -h option.

Before being usable, the mesonh model must be installed following the [tools/INSTALL\_pack\_mesonh.md file](../tools/INSTALL_pack_mesonh.md).

For check\_commit\_mesonh.sh the following environment variables can be set:

  - MNHPACK: directory in which MNH pack will be created (default is $HOME/MesoNH/PHYEX)
  - REFDIR: directory in which reference pack can be found (default is the pack directory near the check\_commit\_mesonh.sh file)
  - TARGZDIR: directory in which the tar.gz file can be found (default is the pack directory near the check\_commit\_mesonh.sh file)

### check\_commit\_testprogs.sh

The check\_commit\_testprogs.sh script runs offline simulations in the directory given
by the environment variable TESTPROGSDIR ($HOME/TESTPROGS will be used if the variable is not set).
This directory must exist.

Script options can be displayed with the -h option.

To be usable the check\_commit\_testprogs.sh script needs input data. The generation and installation of these data are described in the [tools/INSTALL\_testprogs.md file](../tools/INSTALL_testprogs.md).

The goal of the script is to compare outputs between two simulations (to check if bit-reproducibility is achieved or not).
A reference simulation must be performed and save. This reference simulation is run the same way as the
test experiment but cannot be compared to something else:
check\_commit\_testprogs.sh -c -r \<reference\_commit\>

If this reference simulation must become the 'absolute' reference (used when invoking the check\_commit\_testprogs.sh
script with the 'REF' argument), the reference simulation directory (under $TESTPROGSDIR) must be renamed 'ref'.

### prep\_code.sh

This script is used by the different check\_commit\_\* scripts and can be used directly to pre-process the source code.

### testing.sh

This script is designed to be run periodically by cron. It searches for the last commit on a github repository,
use the different check\_commit\_\* scripts to run the test cases and add a comments attached to the gihub commit
to report success or unsecess of these tests. In case of an error, the script can also send an email.

### others

Other scripts are:

  - comp\_DDH.py: compare DDH outputs (used by check\_commit\_ial.sh)
  - compare.py: compare MESO-NH outputs (used by check\_commit\_mesonh.sh)
  - diffNODE.001\_01: compare NODE.0001\_01 output files
  - generate\_standalone\_doc.sh: to generate a standalone doc from the different md files
