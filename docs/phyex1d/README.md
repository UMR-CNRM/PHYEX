# phyex1d: 1D model using PHYEX

## Description
The phyex1D python module runs 1D atmospheric simulations using ideal
cases described in a "common format" file. The physics is provided by
the PHYEX package, also used by the Meso-NH and AROME models.

Online documentation is [here](https://sebastienriettemto.github.io/phyex1d).

## Usage

### Installation

The installation procedure is as follows:
```
# phyex1d and DEPHY-SCM installation
pip install phyex1d
git clone https://github.com/GdR-DEPHY/DEPHY-SCM.git

# PHYEX compilation
git clone https://github.com/UMR-CNRM/PHYEX.git
pip install -r PHYEX/requirements.txt
. PHYEX/tools/env.sh
cd PHYEX/build/with_ecbuild
./make_ecbuild.sh
cd ../../..
```

An easy install procedure also exists:
```
wget https://github.com/SebastienRietteMTO/phyex1d/raw/refs/heads/main/easy_install.sh
chmod +x easy_install.sh
./easy_install.sh INSTALL_DIR MODE
with MODE being 'maxidev', 'dev' or 'user'.
```

### Experiment

The ARMCU case can be run with the following command:
```
. PHYEX/tools/env.sh
phyex1d ./DEPHY-SCM/ARMCU/REF/ARMCU_REF_SCM_driver.nc --exp dt=60 --plot rc.png rc y_var=P
```
And a rc plot is available in the ```rc.png``` file.

### Creating a case from xarray data

Instead of using a netCDF driver file, you can build a case programmatically
from an xarray Dataset using `CaseXarray`. A complete example with profiles,
forcing and plots is available at
[`examples/case_xarray_demo.py`](examples/case_xarray_demo.py).

## Main references
  - PHYEX:
    - [code](https://github.com/UMR-CNRM/PHYEX)
  - MesoNH:
    - [reference](https://doi.org/10.5194%2Fgmd-11-1929-2018)
    - [documentation](http://mesonh.aero.obs-mip.fr/)
    - [code](https://src.koda.cnrs.fr/mesonh/mesonh-code)
  - AROME:
    - [reference](https://doi.org/10.1175/2010MWR3425.1)
  - Common format:
    - [documentation and cases](https://github.com/GdR-DEPHY/DEPHY-SCM)

## TODO (to be integrated when the documentation is rewritten): development installation

Notes on the editable installation of the pyphyex package (scikit-build-core):

- ```pip install -e .``` works with the compiled library: the python modules (src/pyphyex, among
  them phyex1d) are used directly from the source tree, while the files produced by CMake
  (pyphyex/lowlevel.py, pyphyex/_version.py, pyphyex/lib/) are installed in site-packages.
  A modification of the Fortran code (or of the bound routines) needs a new ```pip install -e .```.
- ```build-dir```: by default the CMake build is done in a temporary directory (full rebuild at each
  installation). With ```-C build-dir=build/{wheel_tag}``` the build tree is kept and a new
  installation is incremental (the directory must be git-ignored).
- ```editable.rebuild=true``` (experimental): the CMake build is run at each import of the package,
  so Fortran modifications are taken into account automatically. It needs ```build-dir``` and
  ```--no-build-isolation``` (the build requirements, scikit-build-core and pyfortool, must be
  installed in the environment: the binding is generated at configure time with pyfortool).
  Each import is slower and needs the compilers.
- Proposed usage, on the command line rather than in pyproject.toml (which would also affect
  the normal wheel builds):
  ```
  pip install scikit-build-core pyfortool
  pip install --no-build-isolation -e . -C build-dir=build/{wheel_tag} -C editable.rebuild=true
  ```
  or, without editable.rebuild, run the same ```pip install --no-build-isolation -e . -C build-dir=...```
  again after a Fortran modification.
