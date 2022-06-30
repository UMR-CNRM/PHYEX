# PHYEX
PHYsique EXternalisée

# ABOUT THIS DOCUMENT

This document is a general presentation of the PHYEX package.
More specific documentation can be found:
 - Developer: package organisation, how to contribute, coding norms
 - Integrator: how to merge contributions
 - Offline: how to compile the library and the test programs, how to use the library with python, how to use the test programs
 - Plugging : how to plug the physics package in a model

This document is written using the markdown "language". With pandoc, it can be converted to HTML (pandoc \<filename\>.md -o \<filename\>.html) or PDF (pandoc \<filename\>.md -o \<filename\>.pdf).

# HISTORY

The physics was first developed for the Meso-NH model (http://mesonh.aero.obs-mip.fr/).
Then, a part of the physics have been used to build the AROME model ([Seity et al, 2011](http://dx.doi.org/10.1175/2010MWR3425.1)).
The last step was to extract the physics to build the independent PHYEX package.

# CONTENT

The folowing paramerisations are included in the PHYEX package (see the MesoNH documentation for references):
 - turbulence scheme
 - shallow convection scheme
 - microphysics scheme

In addition to the parametrisatin source code, test programs and a library for python binding are also provided.
