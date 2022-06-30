# PHYEX integrator documentation

# ABOUT THIS DOCUMENT

This document is intended for integrators who are in charge of assembling contributions received through pull requests.

his document is written using the markdown "language". With pandoc, it can be converted to HTML (pandoc \<filename\>.md -o \<filename\>.html) or PDF (pandoc \<filename\>.md -o \<filename\>.pdf).

# BRANCHES AND NORMS

Regarding array-syntax, the applicalble norm depends on the branch:
 - The main branch of PHYEX (and all branches based on main) is written using array-syntax
 - The GPU branch is written using array-syntax with mnh\_expand directives
 - arome specific branches based on the GPU branch are written using DO loops
 - mesonh specific branches based on the GPU branch are written using array-syntax withour mnh\_expand directives

Pull requests can be received on all these kind of branches and must be merged into the main or the GPU branch with according norm.

# TESTS

The source code must follow strict mnh\_expand directives (described in the Developer documentation). The script verify\_mnh\_expand.py must be used to give an additional check.

In addition to the scientific validation, the folowing tests must give the same results (with bit-reproducibility) in each of the model:
 - compilation transforming the mnh\_expand directives in DO loop
 - compilation keeping the array-syntax
 - execution with a different umber of processors

