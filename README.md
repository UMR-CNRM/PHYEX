# PHYEX
PHYsique EXternalisée

## Build

The build systems can be found in the `build` directory.

The PHYEX compilation depends on the fiat (https://github.com/ecmwf-ifs/fiat) package.

### Build with FCM

In the `with_fcm` subdirectory, a build system based on the FCM (https://github.com/metomi/fcm) tool is available.
The command `make_fcm.sh` (call it with the '-h' option to get help):
  - clone the fcm tool
  - clone the fiat package
  - compile the PHYEX and (part of) the fiat package
  - create a shared library (.so)
The resulting shared library (libphyex.so) is under the architecture specific directory created by the script.
