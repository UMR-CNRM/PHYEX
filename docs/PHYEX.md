# PHYEX
PHYsique EXternalisée

## ABOUT THIS DOCUMENT

This document is a general presentation of the PHYEX package.
More specific documentation can be found:

 - [Developer](./Developer.md): package organisation, how to contribute, coding norms
 - [Coding norms](./CodingNorms.md): coding norms
 - [Integrator](./Integrator.md): how to merge contributions
 - [Offline](./Offline.md): how to compile the library and the test programs, how to use the library with python and how to use the test programs
 - [Plugging](./Plugging.md) : how to plug the physics package in a model
 - [Tools](./Tools.md): description of the check\_commit\_\*.sh scripts (to check bit reproducibility between two commits) and of the prep\_code.sh script

This document is written using the markdown language. With pandoc, it can be converted to HTML (pandoc -s \<filename\>.md -o \<filename\>.html) or PDF (pandoc -s \<filename\>.md -o \<filename\>.pdf).

## HISTORY

The physics was first developed for the [Meso-NH model](http://mesonh.aero.obs-mip.fr/).

Then, a part of the physics have been used to build the AROME model ([Seity et al, 2011](http://dx.doi.org/10.1175/2010MWR3425.1)).

The last evolution was to extract the physics to build the independent PHYEX package.

## CONTENT

The following parametrisations are included in the PHYEX package (see the MesoNH documentation for references):

 - turbulence scheme
 - shallow convection scheme
 - microphysics schemes

In addition to the parametrisations source code, test programs and a library for python binding are also provided.
