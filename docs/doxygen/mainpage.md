Doxygen is used to document the source code, general documentation is available
directly in the [docs directory](../..) as Makdown files and are linked here:

  - [Introduction](../PHYEX.md)
  - [Developer](../Developer.md): package organisation, how to contribute, coding norms
  - [Coding norms](../CodingNorms.md): coding norms
  - [Integrator](../Integrator.md): how to merge contributions
  - [Offline](../Offline.md): how to compile the library and the test programs, how to use the library with python and how to use the test programs
  - [Plugging](../Plugging.md) : how to plug the physics package in a model
  - [Tools](../Tools.md): description of the check\_commit\_\*.sh scripts (to check bit reproducibility between two commits) and of the prep\_code.sh script

The main entry points of the package are:

  - [ini\_phyex](@ref ini\_phyex): to setup the different schemes
  - [rain\_ice](@ref rain\_ice): to run the ICE3/ICE4 microphysics scheme
  - [rain\_ice\_old](@ref rain\_ice\_old): to run the old version of ICE3/ICE4 microphysics scheme
  - [ice\_adjust](@ref ice\_adjust): to run the adjustment process linked to the ICE3/ICE4 scheme
  - [turb](@ref turb): to run the turbulence scheme
  - [shallow\_mf](@ref shallow\_mf): to run the shallow convection scheme

In addition, some offline tools are available:

  - [main\_rain\_ice](@ref main\_rain\_ice): to drive the [rain\_ice](@ref rain\_ice) subroutine
  - [main\_rain\_ice\_old](@ref main\_rain\_ice\_old): to drive the [rain\_ice\_old](@ref rain\_ice\_old) subroutine
  - [main\_turb](@ref main\_turb): to drive the [turb](@ref turb) subroutine
  - [main\_ice\_adjust](@ref main\_ice\_adjust): to drive the [ice\_adjust](@ref ice\_adjust) subroutine
  - [main\_shallow](@ref main\_shallow): to drive the [shallow\_mf](@ref shallow\_mf) subroutine
