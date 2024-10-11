# PHYEX PLUGGING DOCUMENTATION

## About this section

The PHYEX parametrisations can be called from the Meso-NH and AROME models, from
test programs and from a driver.
This document is intended for developers who want to plug in the physics in a new model or program.

## Interfaces

The folowing routines are identified as the interfaces of the physics:

  - lima\_adjust\_split
  - ice\_adjust
  - shallow\_mf
  - turb
  - lima
  - rain\_ice, rain\_ice\_old
  - ini\_\phyex

These interfaces are declared in the corresponding modi\_\* files.

## Hooks

The code provided in the common directory is independent, it can be compiled and used without
dependency except the [fiat package](https://github.com/ecmwf-ifs/fiat).
For more interaction with the hosting model, some subroutine can receive a specific implementation.
The following codes already have specific implementations for the Meso-NH and AROME models and are therefore
quite likely to receive a new implementation before plugging into another host model.

  - mode\_budget: to store and/or compute statistics on variable tendencies
  - mode\_msg: to print messages and abort on error

