# PHYEX
PHYsique EXternalisée

Documentation can be found in the [docs directory](./docs/PHYEX.md).

Several presentations were done, the materials can be found on the [wiki](https://github.com/UMR-CNRM/PHYEX/wiki).

Prerequisites:
  - an internet connexion (with access to the github servers) is needed only for the installation
    and, for offline tests (testprogs), when the fiat version to use change
  - python > 3.10
  - some python packages (available on PyPI). All these modules can be installed by ```pip install -r requirements.txt```:
    - pyfortool to transform the source code
    - numpy and pandas for the testprogs
    - epygram and matplotlib for AROME
    - xarray for Meso-NH
    - ctypesForFortran to use the python binding
  - a C compiler (tested with cc 11.4.0)
  - a FORTRAN compiler (tested with ifort and gfortran, but automatic installation only works with gfortran >= 10)
  - some classical unix tools: wget, tar, make and git

Quick Start Guide:
  - open a terminal on a system satisfying the prerequisites and enter the following commands
  - if you don't have a github ssh key or don't know what it is:
    > git clone https://github.com/UMR-CNRM/PHYEX.git  
    > ./PHYEX/tools/INSTALL.sh --ALL --test
  - if you have a github ssh key:
    > git clone git@github.com:UMR-CNRM/PHYEX.git  
    > ./PHYEX/tools/INSTALL.sh --ALL --test --ssh
  - If all goes well, the last line should be "SUCCESS, files are identical"
