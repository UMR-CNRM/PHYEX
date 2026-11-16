# PHYEX PLUGGING DOCUMENTATION

## About this section

The PHYEX parametrisations can be called from the Meso-NH, AROME and LMDZ models, from
test programs and from a driver.
This document is intended for developers who want to plug in the physics in a new model or program.

## Interfaces

The following routines are identified as the interfaces of the physics:

  - lima\_adjust\_split
  - ice\_adjust
  - shallow\_mf
  - turb
  - lima
  - rain\_ice, rain\_ice\_old
  - ini\_phyex

These interfaces are declared in the corresponding modi\_\* files that are generated at compile time.

The code provided in the src/phyex directory is independent, it can be compiled and used without
dependency except the [fiat package](https://github.com/ecmwf-ifs/fiat) (see the
[Offline documentation](./Offline.md) for the compilation of the library).

