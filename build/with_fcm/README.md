Directory organisation:
- the fcm subdirectory contains the fcm tool. In the PHYEX git repository, this directory is empty.
  It will be populated on first call to the make\_fcm.sh script and content is not tracked by git.
- the fiat subdirectory contains the fiat package from the ECMWF. In the PHYEX git repository, this directory is empty.
  It will be populated on first call to the make\_fcm.sh script and content is not tracked by git.
- the arch subdirectory contains architecture specific files. An alternative arch directory can be
  provided on the command line when calling the make\_fcm.sh script
- arch\_\* subdirectories are automatically created by the make\_fcm.sh script and are tracked by git.
- the make\_fcm.sh script will:
  - populate the fcm and fiat directories on first call
  - create the arch\_$ARCH directory, poulate it with arch specific files and a compilation script
  - execute the newly created compilation script

Note: full cleaning is achieved by removing the arch\_\* subdirectories.
Note: documentation of the make\_fcm.sh script can be obtained with the -h option

