# PHYEX plugging documentation

# ABOUT THIS DOCUMENT

The PHYEX parameterizations can be called from the Meso-NH and AROME models, from
test programs and from a driver.
This document is intended for developers who want to plug in the physics in a new model or program.

This document is written using the markdown "language". With pandoc, it can be converted to HTML (pandoc -s \<filename\>.md -o \<filename\>.html) or PDF (pandoc -s \<filename\>.md -o \<filename\>.pdf).

# INTERFACES

The folowing routines are identified as the interface of the physics:

  - lima\_adjust
  - ice\_adjust
  - shallow\_mf
  - turb
  - lima, lima\_warm, lima\_cold and lima\_mixed
  - rain\_ice, rain\_ice\_old
  - ini\_\* **TODO: list the different ini subroutine needed**

This interface is declared in the corresponding modi\_\* files.

# HOOKS

The code provided in the common directory is independent, it can be compiled and used without
dependency except the [fiat package](https://github.com/ecmwf-ifs/fiat).
For more interaction with the hosting model, some subroutine can receive a specific implementation.
The following codes already have specific implementations for the Meso-NH and AROME models and are therefore
quite likely to receive a new implementation before plugging into another host model.

  - mode\_budget: to store and/or compute statistics on variable tendencies
  - mode\_msg: to print messages and abort on error

