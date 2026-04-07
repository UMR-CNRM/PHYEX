# PHYEX - Physics Parameterization Package

PHYEX is a Fortran library for atmospheric physics parameterizations with Python bindings.

## Features

- Mixed-phase cloud microphysics (ICE_ADJUST, RAIN_ICE)
- Shallow convection schemes
- Python bindings via Cython for easy integration

## Building

### Fortran Library Only

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
cmake --install .
```

### Python Package with scikit-build-core

The project uses [scikit-build-core](https://scikit-build-core.readthedocs.io/) for building Python wheels.

#### Requirements

- Python >= 3.11
- gfortran (GNU Fortran compiler)
- CMake >= 3.15
- NumPy < 2.0
- Cython

#### Installation from source

```bash
# Install in editable mode for development
pip install -e .

# Or build a wheel
pip install build
python -m build
```

#### Installation from wheel

```bash
pip install dist/phyex-0.1.0-*.whl
```
## Installation from spack

The package is installable from a spack repo 

- To install the python extensions : 

```bash
    spack install phyex +python +double
```

- To install fortran package :

```bash
    spack install phyex ~python +double
```

## Usage

### Python

```python
import numpy as np
from phyex import ice_adjust, rain_ice, init_rain_ice

# Initialize microphysics
init_rain_ice(timestep=60.0, dzmin=20.0, krr=6, hcloud="AROME")

# Prepare input arrays (Fortran-contiguous, float32)
nlon, nlev = 10, 50
pabs = np.ones((nlon, nlev), dtype=np.float32, order='F') * 85000.0
# ... initialize other arrays ...

# Call ice_adjust
ice_adjust(
    timestep=60.0, krr=6,
    sigqsat=sigqsat, pabs=pabs, sigs=sigs, th=th,
    exn=exn, exn_ref=exn_ref, rho_dry_ref=rho_dry_ref,
    rv=rv, rc=rc, ri=ri, rr=rr, rs=rs, rg=rg,
    cf_mf=cf_mf, rc_mf=rc_mf, ri_mf=ri_mf,
    rvs=rvs, rcs=rcs, ris=ris, ths=ths,
    cldfr=cldfr, icldfr=icldfr, wcldfr=wcldfr
)
```

## Project Structure

```
PHYEX/
├── aux/              # Auxiliary routines
├── turb/             # Turbulence parameterizations
├── micro/            # Microphysics schemes
├── conv/             # Convection schemes
├── bridge/           # Python-Fortran bridge
│   ├── phyex_bridge.F90       # Fortran C-interop layer
│   ├── _phyex_wrapper.pyx     # Cython wrapper (CPU)
│   └── _phyex_wrapper_acc.pyx # Cython wrapper (GPU/OpenACC)
├── CMakeLists.txt    # Build configuration
└── pyproject.toml    # Python package configuration
```

## Development

### Testing

```bash
pip install -e ".[dev]"
pytest tests/
```

### Code Formatting

```bash
black phyex/
ruff check phyex/
```

## References

Based on the PHYEX physics package from Météo-France.

## License

Apache-2.0
