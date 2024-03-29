# PHYEX INTEGRATOR DOCUMENTATION

## About this section

This document is intended for integrators who are in charge of assembling contributions received through pull requests.

## Branches and norms

Regarding array-syntax, the [applicable norm](./CodingNorms.md) depends on the branch:

  - The legacy main branch of PHYEX (and all branches based on main) was written using array-syntax
  - The master branch is written using array-syntax with mnh\_expand directives
  - arome specific branches based on the master branch are written using DO loops
  - mesonh specific branches based on the master branch are written using array-syntax

Pull requests can be received on all these kind of branches and must be merged into the master branch with according norm.

## Normal workflow for a contribution developed in AROME and HARMONIE-AROME

### Scientific contributions until cycle 49t1

![](./AROMEworkflow1.svg)

The pull request comes from the IAL integrator. It must be based on an arome specific branch.

Details for point 6, the PHYEX administrator:

  - validates (see [below](#tests)) the contribution
  - integrates the contribution in the arome branch and merges it in the master branch
  - regularly, he tags a new (minor) version of the master branch
  - when asked by the IAL integrator, he builds a new arome specific branch
  - when an arome specific branch is used in an official cycle, the arome specific branch is tagged accordingly

### Refactoring contributions from now on, and scientific contributions after cycle 49t1

![](./AROMEworkflow2.svg)

The pull request comes directly from a developer. It must be based on an arome specific branch.

Details for point 6:

  - The PHYEX administrator checks the pull requests in the other applications (see [below](#tests))
  - The IAL integrator integrates the contribution on the arome specific branch
  - The PHYEX administrator
    - integrates the modifications in the master branch
    - regularly, tags a new (minor) version of the master branch
    - when asked by the IAL integrator, builds a new arome specific branch (see [below](#code-preparation))
    - when an arome specific branch is used in an official cycle, the arome specific branch is tagged accordingly

### PHYEX - IAL compatibility

Because the interfaces between the physics and the rest of the model can change, the right version of IAL to use with PHYEX must be referenced in the 'src/arome/ial\_version.json' file.

If no IAL version is able to receive directly the PHYEX physics, an 'ext' directory can be added to hold the files that must be modified on top of the IAL version described in the 'src/arome/ial\_version.json' file.

In addition, if some files must be excluded from the IAL compilation, they must be listed in the 'src/arome/gmkpack\_ignored\_files' file.

## Normal workflow for a contribution developed in MESONH

The developer sends its pull request on the Méso-NH repository (the physics source code is embedded in the model source code).

Integration details:

  - The Meso-NH integrator extracts, from the different pull requests, what concern the PHYEX repository and send a pull request on PHYEX based on a mesonh specific branch
  - The PHYEX administrator:
    - validates (see [below](#tests)) the contribution
    - integrates the contribution in the mesonh branch and merges it in the master branch
    - regularly, he tags a new (minor) version of the master branch
    - when asked by the Meso-NH integrator, he builds a new mesonh specific branch (see [below](#code-preparation))
    - when a mesonh specific branch is used in an official release, the mesonh specific branch is tagged accordingly

## Normal workflow for another contribution

Pull requests must be based on the master branch.

The PHYEX administrator:

  - validates (see [below](#tests)) the contribution
  - integrates the contribution in the master branch
  - regularly, he tags a new (minor) version of the master branch

## Tests

The source code must follow strict mnh\_expand directives (described in the [Coding Norms documentation](./CodingNorms.md)). The script verify\_mnh\_expand.py must be used to give an additional check.

In addition to the scientific validation, the following tests must give the same results (with bit-reproducibility) in each of the model (arome, mesonh and testprogs):

  - compilation transforming the mnh\_expand directives in DO loop
  - compilation keeping the array-syntax
  - execution with a different number of processors

When possible, the new version of PHYEX must reproduce the old results (scientific modifications must be activated with namelist keys).

## Code preparation

The source code stored in the master branches must be usable by all the models. But these models can have contradictory constraints. To bypass this difficulty, the source code is preprocessed before being included in the compilation environment of each model.

This preprocessing step can be done on the fly (in this case the preprocessing tools must be available aside of the compilation tools), or the result of the preprocessing can be stored in the PHYEX package (in this case, the preprocessing is done once and can be used by several users).
This second possibility is useful to historize the source code really used during the model compilation and enables contributions to the PHYEX package without the need of the preprocessing tools.

The preprocessed versions of the source code are put in branches named \<model\>\_\<commit\> where \<model\> is the name of the model for which the source code have been preprocessed and \<commit\> is the commit hash used as a basis.

The preprocessing tools are described in the [Tools documentation](./Tools.md).

## Code validation

The testing.h script can be automatically run to help for validation. This script is described
in the [Tools documentation](./Tools.md).

