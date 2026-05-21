# INSTALLATION NEEDED FOR THE TEST PROGRAMS

## ABOUT THIS DOCUMENT 

This document is intended for persons who want to use the testprogs programs.
These programs need data. This document describes how to generate the data, the the programs can be used directly or through the check\_commit\_testprogs.sh script.

This document is written using the markdown language. With pandoc, it can be converted to HTML (pandoc -s \<filename\>.md -o \<filename\>.html) or PDF (pandoc -s \<filename\>.md -o \<filename\>.pdf). 

## DATA

There are two options.

On one hand, the data can be generated (please, refer to the [Offline](../docs/Offline.md) documentation) and, once produced, be put in the corresponding directories under ${PHYEXCONF}/testprogs\_data.

On the other hand, if {PHYEXCONF}/testprogs\_data is empty, a reduced dataset will be automatically dowloaded and installed in ${PHYEXCONF}/testprogs\_data when check\_commit\_testprogs.sh is used.
