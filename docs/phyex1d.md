# PHYEX1D: 1D MODEL USING PHYEX

## About this section

The phyex1d module of the pyphyex python package runs 1D atmospheric simulations using ideal
cases described in a "common format" file (or built from xarray data). The physics is provided
by the PHYEX package (through the python binding, pyphyex.lowlevel), the same as the one used by
the Meso-NH and AROME models. A simplified physics, not using PHYEX, is also available.

## Installation

phyex1d is part of the pyphyex package (see the [Offline documentation](./Offline.md)):
```
pip install ./PHYEX       # from the PHYEX repository
```

The ecRad radiation scheme is used through the pyecrad package, installed as a dependency from
PyPI (as pyecrad-SR, built from the [pyecrad-main branch](https://github.com/SebastienRietteMTO/ecrad/tree/pyecrad-main)
of ecRad, until pyecrad is published by ECMWF). Prebuilt wheels are available for Linux (x86\_64 and
aarch64); elsewhere, the installation compiles ecRad from the source distribution: a Fortran compiler,
CMake and NetCDF-Fortran are then needed (see the pyecrad documentation for the build options).
A development version of pyecrad can be used
instead with the 'pyecrad' option of the PHYEX physics (path to the pyecrad package directory,
the libecrad4py library being given by the PYECRAD\_LIBRARY environment variable).

This installs the following commands (each of them has a -h option):

  - phyex1d-cases\_install: install the DEPHY-SCM cases (download or link)
  - phyex1d-cases\_list: list the installed DEPHY-SCM cases
  - phyex1d-exp\_run: run experiments on a case
  - phyex1d-exp\_list: list the experiments having results
  - phyex1d-exp\_plot: plot the results of experiments

## Cases

The cases in the common format are provided by the [DEPHY-SCM](https://github.com/GdR-DEPHY/DEPHY-SCM) repository.
They are installed in ${PHYEXCONF}/phyex1d/DEPHY-SCM (the PHYEXCONF environment variable defaults to ${HOME}/.phyex):
```
phyex1d-cases_install download       # clone (or update) the DEPHY-SCM repository (about 1 GB)
phyex1d-cases_install link PATH      # or use an existing DEPHY-SCM repository (symbolic link)
phyex1d-cases_list                   # list of the cases (CASE/SUBCASE) and of their driver files
```

## Experiments

An experiment is run on a case given by the name of an installed case (CASE/SUBCASE) or by the path to a
netCDF driver file. For instance, the ARMCU case can be run with the PHYEX physics (default) and with the
simplified physics:
```
phyex1d-exp_run ARMCU/REF -o out --name armcu --exp dt=60 name=phyex \
                                              --exp dt=60 name=simple class=SimplePhysicsTQ
phyex1d-exp_list -o out --name armcu
phyex1d-exp_plot -o out --name armcu --exp phyex --plot rc.png y_var=P rc
```
And a rc plot is available in the rc.png file.

Each --exp option describes an experiment with key=value strings: dt (timestep), name, grid, class (physics
class: PhysicsAromeTQ, PhysicsAromeThetaR, SimplePhysicsTQ or SimplePhysicsThetaR) and the options of the
physics class (namelist for the PHYEX physics, activation of the different processes for the simplified
physics...). The list of the keys, with their default values, is given by phyex1d-exp\_run -h.

The results are stored in OUTPUT\_DIR/\<comparison name\>/\<experiment name\>.hdf5 (OUTPUT\_DIR is given by
the -o option and the comparison name by the --name option). An experiment is not run again if its result
file exists (use --force).

Grids (.grid files) and namelists (.namel files) are searched in ${PHYEXCONF}/phyex1d/grids and
${PHYEXCONF}/phyex1d/namelists, then in the package. The description of the .grid file format is given by
the docstring of the pyphyex.phyex1d.Grid class.

## Usage from python

The experiments can also be run from python, with the ExperimentRunner class (pyphyex.phyex1d module).
Instead of using a netCDF driver file, a case can be built programmatically from an xarray Dataset using
the CaseXarray class. Complete examples with profiles, forcing and plots are available in the
src/pyphyex/phyex1d/examples directory (case\_xarray\_demo.py with the PHYEX physics, and
case\_xarray\_simple.py with the simplified physics).
