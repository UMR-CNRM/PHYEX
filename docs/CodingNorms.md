# PHYEX coding norms documentation

## ABOUT THIS DOCUMENT

This document is intended for developers and integrators and describes the coding norms to use.

This document is written using the markdown language. With pandoc, it can be converted to HTML (pandoc -s \<filename\>.md -o \<filename\>.html) or PDF (pandoc -s \<filename\>.md -o \<filename\>.pdf).

## CODING NORMS

### Namelists
We must be able to reproduce (binary comparison of the output files) the model results before and after code modifications. It means that every modification must be controlled by a namelist key (with the exception of bug corrections).

### File names
The fortran file names use a capital F letter (eg: foo.F90) except if working in a Meso-NH branch (mesonh\_\<commit\>) or in the folder (src/mesonh) specific to the Meso-NH model.

Names for the module:

  - modd\_ for module containing only variable declaration (eg: tuning parameters)
  - modi\_ for module containing only interface declaration
  - modn\_ for namelist declaration
  - mode\_ for module containing executable source code (subroutine or function)

### When using mode\_ or modi\_?
When writing a new subroutine, should we put it in a module (in a mode\_ file) or should we write the subroutine in a file and write the interface bloc in another file (modi\_ file)?

The answer depends on whether the routine is the 'main' routine of the parametrisation or not. If it is the 'main' routine, the interface bloc is declared apart, if not we can use a module.
The idea behind is to break compilation dependency at the parametrisation level, and to isolate the interface declaration of the different routines that must be plugged in the hosting model.

### Norm
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

Regarding array-syntax, code is written using array-syntax in the legacy main branch, using array-syntax with mnh\_expand directives in the master branch and in mesonh specific branches based on the master branch, using DO loops in arome specific branches based on the master branch. If in doubt, check what is done in other routines in the branch you are working in.
Be carrefull when using the mnh\_expand directives, code must respect some constraints:

  - parenthesis after array variables are mandatory (no A=B+C, but A(:,:)=B(:,:)+C(:,:))
  - no space between array variables and the opening parenthesis (no A (:)=B (:), but A(:)=B(:))
  - same bounds as declared in the mnh\_expand directive should be used in the array-syntax (A(D%NIB;D%NIE)=...)

A tool (verify\_mnh\_expand.py) can help at checking the validity of the written code.

For the master branch (and branches on master, including model specific branches):

  - except variables declared with the PARAMETER attribute, no variable from modules can be used in the physics. Variables must be put in a type received by interface.
  - subroutines or functions must not be called from within a loop on horizontal or vertical dimensions (see below for exception)
  - functions returning arrays must be rewritten as subroutine

Call to external subroutine in loop on horizontal or vertical dimensions must be suppressed in the master version. If possible, the call must be put outside of the loop (acting on the full array as a whole) or the subroutine must be put in the CONTAINS part but, in this case, the included subroutine cannot use local array. There are 3 cases:

  - the subroutine doesn't use local array: subroutine is put in an include file (with the .h extension) and included with the fortran INCLUDE statement.
  - the subroutine use local arrays but it is called from only one place in the code: the source code of the subroutine is moved (no INCLUDE) in the CONTAINS part and the array declarations are moved in the main subroutine.
  - the subroutine use local arrays and is called from several places: the previous technique is not recommended. The source code is put in an include file (with the .h extension) and an extra argument is provided to the subroutine and is used as a buffer so there is no more need to declare local arrays in the called subroutine.

### Budgets

In Meso-NH, the budget can be used in two ways:

  - by giving to the budget machinery the tendency due to a given process
  - by giving to the budget machinery the total tendency (S variable) before and after a given process. The budget mechanism recomputes by difference the tendency only due to the given process.

In AROME, we cannot provide the total tendency (S variable) before the process. This total tendency is stored internally by the machinery but cannot be set to a different value before doing a computation.

The physics package must be usable from AROME and Meso-NH, several examples are given:

Invalid for AROME:
```
budget_store_init(tempo_s)
modification of tempo_s
budget_store_end(tempo_s)
```

Valid:
```
budget_store_init(pronostic_s) #useless for AROME, but needed for Meso-NH
modification of pronostic_s
budget_store_end(pronostic_s)
```

Valid:
```
computation of delta tempo_s
budget_store_add(delta tempo_s)
```
