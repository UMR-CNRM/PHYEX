# TOOLS INSTALLATION

## ABOUT THIS DOCUMENT 

This document is intended for persons who want to install part or all of the tools provided with the PHYEX package.
   
This document is written using the markdown language. With pandoc, it can be converted to HTML (pandoc -s \<filename\>.md -o \<filename\>.html) or PDF (pandoc -s \<filename\>.md -o \<filename\>.pdf).

## INSTALLATIONS

1. IAL
   Some packages must be installed:

   - gmkpack
   - fypp python module: pip3 install fypp
   - ial-build python module: pip3 install ial-build
   - ecbundle: pip install git+https://github.com/ecmwf/ecbundle

   Note that some of these tools can be pre-installed (eg. on belenos, module use ~mary/public/modulefiles; module load IAL-build)

2. MESONH REFERENCE PACK
   The reference pack for Meso-NH must be installed. Instructions can be found in
   [INSTALL\_pack\_mesonh](./INSTALL_pack_mesonh.md)

3. MNH\_EXPAND
   This tool is composed of two parts: the filepp software, and the MNH\_Expand\_Array package.
   It is needed by the different check\_commit\_\*.sh scripts, depending on options used at execution.
   Installation is described in [INSTALL\_mnh\_expand](./INSTALL_mnh_expand.md)

4. TESTPROGS REFERENCE DATA
   The script check\_commit\_testprogs.sh needs reference data.
   The data generation is explained in [INSTALL\_testprogs](./INSTALL_testprogs.md)
