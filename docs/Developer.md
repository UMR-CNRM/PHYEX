# PHYEX developer documentation

# ABOUT THIS DOCUMENT

This document is intended for developers who want to contribute to the PHYEX package.
Developer who is interested in plugging the physics in a new model can refere to the Plugging documentation.
The topics covered are as follows:
 - [Package organisation](#PACKAGE-ORGANISATION)
 - [Code preparation](#CODE-PREPARATION)
 - [Coding norms](#CODING-NORMS)
 - [Pull requests](#PULL-REQUESTS)

This document is written using the markdown "language". With pandoc, it can be converted to HTML (pandoc \<filename\>.md -o \<filename\>.html) or PDF (pandoc \<filename\>.md -o \<filename\>.pdf).

# PACKAGE ORGANISATION

The package contains the folowing directories:
 - docs: for documentation
 - build: an autonomous build system is included in the package. Its usage is covered in the Offline documentation
 - src/common: the main source code which is the basis for all models
 - src/\<model\>: the source code specific to one model that must replace source code found in the common directory

In addition to this organisation, the package uses git branches. The main branches are as follows:
 - main: source code without rewriting for GPU transformation (used for official Meso-NH source code)
 - GPU: source code adapted for GPU transformations (used for official AROME source code, starting from the 48t3 cycle)
 - arome\_\<commit\>: source code ready for inclusion in the AROME compilation environment (the generation of such a branch is described in [Code preparation](#CODE-PREPARATION))
 - testprogs\_data: modified source code used to generate samples for the test programs (more on this topic in the Offline documentation)

# CODE PREPARATION

The source code stored in the main and GPU branches must be usable by all the models. But these models can have contradictory constraints. To bypass this difficulty, the source code is preprocessed before being included in the compilation environment of each model.

This preprocessing step can be done on the fly (in this case the preprocessing tools must be available aside of the compilation tools), or the result of the preprocessing can be stored in the PHYEX package (in this case, the proeprocessing is done once and can be used by several users).
This second possibility is usefull to historize the source code really used during the model compilation and enables contributions to the PHYEX package without the need of the preprocessing tools.

The preprocessed versions of the source code are put in branches named \<model\>\_\<commit\> where \<model\> is the name of the model for which the source code have been preprocessed and \<commit\> is the commit hash used as a basis.

During the initial development step, it was found easier to store the preprocessing tools outside of the PHYEX repository in the [PHYEX\_tools](https://github.com/SebastienRietteMTO/PHYEX_tools) repository. In the future, they will be moved in the PHYEX repository.

Installation and usage of the preprocessing tools are described in the PHYEX\_tools package.

# CODING NORMS

## File names
The fortran file names use a capital F letter (eg: foo.F90) except if working a branch (mesonh\_\<commit\>) or in the folder (src/mesonh) specifci to the Meso-NH model.

Names for the module:
 - modd\_ for module containing only variable declaration (eg: tuning parameters)
 - modi\_ for module containing only interface declaration
 - modn\_ for namelist declaration
 - mode\_ for module containing executable source code (subroutine or function)

## When using mode\_ or modi\_?
When writing a new subroutine, should we put it in a module (in a mode\_ file) or should we write the subroutine in a file and write the interface bloc in another file (modi\_ file)?

The answer depends on whether the routine is the 'main' routine of the parameterisation or not. If it is the 'main' routine, the interface bloc is declared apart, if not we can use a module.
The idea behind is to break compilation dependency at the parameterisation level, and to isolate the interface declaration of the different routines that must be pluged in the hosting model.

## Norm
Several constraints are imposed:
 - The code must be written with up to 132 characters per line.
 - CODE IS IN CAPITAL LETTERS! comments in small letters
 - All variables must be declared: IMPLICIT NONE
 - except in rare cases, use automatic arrays, no allocatable
 - dimensions of dummy argument arrays are explicit (no (:,:))
 - use parenthesis when manipulating arrays (eg: A(:,:)=B(:,:)+C(:,:) instead of A=B+C)

The variables are named according to the doctor norm:
|Type / Status | INTEGER    | REAL       | LOGICAL     | CHARACTER      | TYPE             |
|--------------|------------|------------|-------------|----------------|------------------|
|Global        | N          | X          | L (not LP)  | C              | T (not TP,TS,TZ) |
|Dummy argument| K          | P (not PP) | O           | H              | TP               |
|Local         | I (not IS) | Z (not ZS) | G (not GS)  | Y (not YS, YP) | TZ               |
|Loop control  | J (not JP) | -          | -           | -              | -                |

Regarding array-syntax, code is written using array-syntax in the main branch and in mesonh specific branches based on the GPU branch, using array-syntax with mnh\_expand directives in the GPU branch, using DO loops in arome specific branches based on the GPU branch. If in doublt, check what is done in other routines in the brach you are working in.
Be carrefull when using the mnh\_expand directives, code must respect some constraints:
 - parenthesis after array variables are mandatory (no A=B+C, but A(:,:)=B(:,:)+C(:,:))
 - no space between array variables and the opening parenthesis (no A (:)=B (:), but A(:)=B(:))
 - same bounds as declared in the mnh\_expand directive should be used in the array-syntax (A(D%NIB;D%NIE)=...)
A tool (verify\_mnh\_expand.py) can help at checking the validity of the written code.

For the GPU branch (and branches on GPU, including model specific branches):
 - except variables declared with the PARAMETER attribute, no variable from modules can be used in the physics. Varaibles must be put in a type received by interface.
 - subroutines or functions must not be called from within a loop on horizontal or vertical dimensions (see below for exception)
 - functions returning arrays must be rewritten as subroutine

Call to external subroutine in loop on horizontal or vertical dimensions must be suppressed in the GPU version. If possible, the call must be put outside of the loop (acting on the full array as a whole) or the subroutine must be put in the CONTAINS part but, in this case, the included subroutine cannot use local array. There are 3 cases:
 - the subroutine does't use local array: subroutine is put in an include file (with the .h extension) and included with the fortran INCLUDE statement.
 - the subroutine use local arrays but it is called from only one place in the code: the source code of the subroutine is moved (no INCLUDE) in the CONTAINS part and the array declarations are moved in the main subroutine.
 - the subroutine use local arrays and is called from several places: the previous technique is not recommended. The source code is put in an include file (with the .h extension) and an extra argument is provided to the subroutine and is used as a buffer so there is no more need to declare local arrays in the called subroutine.

# PULL REQUESTS
This section deals with the pull request procedure from the developer point of view. The integrator point of view is described in the Intergator documentation.

To contribute to the PHYEX repository, developer must fork the repository, contribute on the main or on the GPU branch and send a pull request. Alternatively, a contribution on a model specific branch is also possible (especially for minor modifications).

If a modification must be applied to the main and to the GPU branches, the pull request must be made on the main branch (and will be merged into the GPU branch).

