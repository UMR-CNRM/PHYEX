# PHYEX
PHYsique EXternalisée

Documentation can be found in the [docs directory](./docs/PHYEX.md).

Several presentations were done, the materials can be found on the [wiki](https://github.com/UMR-CNRM/PHYEX/wiki).

Prerequisites:
  - an internet connexion (with access to the github servers) is needed only for the installation
  - python > 3.7 (but only tested with version 3.10)
  - a C compiler (tested with cc 11.4.0)
  - a FORTRAN compiler (tested with ifort and gfortran, but automatic installation only works with gfortran)
  - some classical unix tools: wget, tar, make and git

Quick Start Guide:
  - open a terminal on a system satisfying the prerequisites and enter the following commands
  - if you don't have a github ssh key or don't know what it is:
    > git@github.com:UMR-CNRM/PHYEX.git  
    > ./PHYEX/tools/INSTALL.sh --ALL --test 
  - if you have a github ssh key:
    > git clone git@github.com:UMR-CNRM/PHYEX.git  
    > ./PHYEX/tools/INSTALL.sh --ALL --test --ssh
  - If all goes well, the last line should be "SUCCESS, files are identical"
