# PHYEX DEVELOPER DOCUMENTATION

## About this section

This section is intended for developers who want to contribute to the PHYEX package.
Developer who is interested in plugging the physics in a new model can refer to the Plugging documentation.
The topics covered are as follows:

  - [Package organisation](#package-organisation)
  - [Contribution workflow for AROME and HARMONIE-AROME developers](#contribution-workflow-for-arome-and-harmonie-arome-developers)
  - [Contribution workflow for Méso-NH developers](#contribution-workflow-for-mesonh-developers)
  - [Contribution workflow for other developers](#contribution-workflow-for-other-developers)

## Package organisation

The package contains two kinds of branches:

  - generic branches which contain codes for all the models and applications (eg: master branch)
  - model specific branches which are automatically derived from generic branches (eg: arome\_\<commit\_hash\>, mesonh\_\<commit\_hash\>, offline\_\<commit\_hash\>)

The directories found in the package are different depending on the branches (generic or model specific).

For model specific branches, only the source code adapted for a given model is present (one directory per parametrisation and an aux directory). No compilation engine or scripts are present in these branches. They are intended to be included directly in the compilation system of the hosting model.

The generic branches contains the following directories:

  - docs: for documentation
  - build: an autonomous build system is included in the package. Its usage is covered in the [Offline documentation](./Offline.md)
  - src/common: the main source code which is the basis for all models
  - src/\<model\>: the source code specific to one model that must replace or complement the source code found in the common directory
  - tools: scripts to build model specific branches and run test cases (described in the [Integrator](./Integrator.md) documentation).

Here is a short description of the different generic branches:

  - master: source code adapted for GPU transformations
  - testHUGE: modified source code to check if incomplete NPROMA blocs are working well (only useful for testing)
  - testprogs\_data: modified source code used to generate samples for the test programs (more on this topic in the [Offline documentation](./Offline.md))

## Contribution workflow for AROME and HARMONIE-AROME developers

The build systems are evolving.
Until cycle 49t1 (included), the physics source code is directly included in the source code tree.
After cycle 49t1, the physics source code (as well as other model parts such as ectrans, fiat...) will be available as a _component_ of a _bundle_.

This evolution will impact the way to contribute to the PHYEX repository.

Whatever is the cycle, the AROME and HARMONIE-AROME developers only see codes coming from arome (or offline) specific branches (branches named arome\_\<commit\_hash\> or offline\_\<commit\_hash\>). This code is ready for inclusion (array-syntax already transformed into DO loops for instance).

Said differently, developers do not need to manipulate code transformation tools.

The workflow was chosen so that the developers would not have to change their working methods several times:

- Developers who have a scientific contribution will submit their pull request on the IAL repository until the ecbundle mechanism is active (49t2 or 50t1), afterwards they will submit the pull request directly on PHYEX.
- Developers who work on the refactoring cannot use the IAL repository as a starting point and must use the source codes in the PHYEX repository. They will directly use the PHYEX repository.

### Scientific contributions until cycle 49t1

Who: developers with scientific contributions based on cycles 48t1, 48t2, 48t3, 49 and 49t1 (as long as PHYEX is not a bundle component).

Workflow summary: because the physics source code is still included in the IAL source code, pull requests concerning the physics continue to follow the same path as before (ie pull requests are submitted to the IAL repository). Afterwards, the IAL integrator will submit a pull request to the PHYEX repository with only the relevant files.

![](./AROMEworkflow1.svg)

Workflow details (getting the source code in blue, pull request in red, integration in green):

  - 1: PHYEX administrator sends (pull request) the content of a specific arome branch to the IAL Integrator. The IAL integrator tags a new release of IAL.
  - 2: The AROME or HARMONIE-AROME developer forks the IAL repository
  - 3: The AROME or HARMONIE-AROME developer compiles, executes, modifies the source code in its environment
  - 4: The AROME or HARMONIE-AROME developer sends a pull request to the IAL repository
  - 5: The IAL integrator extracts the physics source files and sends a pull request to the PHYEX repository
  - 6: The PHYEX administrator checks and integrates the modifications in the master branch and, eventually, produce a new arome specific branch for future integration in IAL

### Refactoring contributions from now on, and scientific contributions after cycle 49t1

Who: developers with scientific contributions based on cycles 50 and following (as soon as PHYEX is a bundle component); and developers with GPU-refactoring contributions.

Workflow summary: after the cycle 49t1 (probably starting from cycle 50), AROME and HARMONIE-AROME will become a bundle. Il will be built with source codes coming from various places. One of these places will be the PHYEX repository. Pull requests must be sent to each modified components of the bundle.

Developer must use a model specific branch (arome\_\<commit\_hash\> when working with the model, or offline\_\<commit\_hash\> when working with the offline tools).
These branches receive tags based on the master branch version. For example the commit, in the master branch, corresponding to the version 1.0.0 of PHYEX will receive the tag "v1.0.0".
The arome specific commit corresponding to this version will be tagged "v1.0.0\_arome".

![](./AROMEworkflow2.svg)

Workflow details (getting the source code in blue, pull request in red, integration in green):

  - 1 and 2: AROME or HARMONIE-AROME developer forks the different repositories needed to build the model
  - 3: AROME or HARMONIE-AROME developer compiles, executes, modifies the source code in its environment
  - 4 and 5: AROME or HARMONIE-AROME developer sends pull requests to the different repositories where files have been modified
  - 6: The PHYEX administrator checks the pull requests in the other applications, the IAL integrator integrates on the arome specific branch; then the PHYEX administrator integrates the modifications in the master branch and, eventually, produce a new arome specific branch for future integration in IAL

### Special notes for building the AROME or HARMONIE-AROME model from PHYEX until cycle 49t1 included

Because the interfaces between the physics and the rest of the model can change, one have to choose the right version of IAL to use with PHYEX.
The file 'src/arome/ial\_version.json' contains a description of this IAL version.

If no IAL version suits correctly, this json file is accompanied by the 'ext' directory and/or by the 'src/arome/gmkpack\_ignored\_files' file.

To build the model from PHYEX, you must:

- checkout the IAL source code using the version described in the file src/arome/ial\_version.json
- remove the directories 'mpa/\*/internals' and 'mpa/\*/modules' (if they still exist in IAL, eg: 48t3)
- put the PHYEX directories 'aux', 'conv', 'micro' and 'turb' into a directory (at the same level as 'mpa') named 'phyex'
- if the 'ext' directory exists, dispatch its content into the subdirectories of IAL
- remove from the source tree the files listed in the 'src/arome/gmkpack\_ignored\_files' file

However, for scientific contributions to 49t1, we suggest scientists to use the physics code present in IAL, rather than from PHYEX.

## Contribution workflow for MESO-NH developers

The physics source code is embedded in the Méso-NH source code.

The physics source code comes directly from a mesonh specific branch (these branches are named mesonh\_\<commit\_hash\>) which contain code ready for use in the Méso-NH model (array-syntax...).

Pull requests concerning the physics continue to follow the same path as before (ie pull requests are submitted to the Meso-NH repository). The Meso-NH integrator will submit a pull request to the PHYEX repository with only the relevant files.

## Contribution workflow for other developers

Other developers must work with source code coming directly from the master branch. They issue pull requests directly on this branch as usual with git repositories.

## Code documentation

The code must contain comments.
The documentation rules to folllow are described in the [CodingNorms](./CodingNorms.md) document.
