# OFFLINE COMPILATION

## Directory organisation

The different directories are:
  - the fcm subdirectory contains the fcm tool. In the PHYEX git repository, this directory is empty.
    It will be populated on first call to the make\_fcm.sh script and content is not tracked by git.
  - the fiat subdirectory contains the fiat package from the ECMWF. In the PHYEX git repository, this directory is empty.
    It will be populated on first call to the make\_fcm.sh script and content is not tracked by git.
  - the arch subdirectory contains architecture specific files. An alternative arch directory can be
    provided on the command line when calling the make\_fcm.sh script. And, arch files are also (and in
    first place) looked for in the ${HOME}/.phyex/fcm\_arch directory.
  - arch\_\* subdirectories are automatically created by the make\_fcm.sh script and are not tracked by git.

## Compilation

The make\_fcm.sh script will:
  - populate the fcm and fiat directories on first call
  - create the arch\_$ARCH directory, poulate it with arch specific files and a compilation script
  - populate the compilation directory with source code (using the prep\_code.sh script). Source code
    transformation can occur during this step (in particular using the PYFT_OPTS variable set
    in the architecture env file).
  - execute the newly created compilation script

Note: full cleaning is achieved by removing the arch\_\* subdirectories.
Note: documentation of the make\_fcm.sh script can be obtained with the -h option

## Architecture files

By default, the code is compiled using the _gnu_ configuration. This configuration is stored in the
arch-gnu.env and arch-gnu.fcm files in the arch directory.
Other acthitectures exist in this directory.
In addition, other architectures can be added in the ${HOME}/.phyex/fcm\_arch directory.

The fcm file contains the compilation options to be used by the fcm tool.
Such a file can be build using a gmkfile with the command:
  make\_fcm.sh --gmkfile \<GMKFILE\> --arch \<new arch name\>
The env file is sourced several times:
 - before populating the build directory with the source code.
   Indeed, during this second step, the source code is copied or cloned and transformed by the pyft\_tool.py script.
   The active transformations are controlled by the --noexpand/--expand options given to the
   different check\_commit\_\* scripts and by the PYFT\_OPTS that can be set in the env file (only for testprogs).
   The syntax is given below.
 - just before the compilation step for loading modules or to defined (LIBS variable)
   the list of system libraries to link with (defaults to 'rt ld' to link with librt and libdl).
 - just before execution to set environment variable specific to this architecture
   needed during execution
   
The PYFT\_OPTS environment variable can contain a multi-lines string.
For each file, the PYFT\_OPTS is read line by line and the last applicable line is used.
A line can take one of these two forms:
  - FILE\_DESCRIPTOR:=:OPTIONS
    where FILE\_DESCRIPTOR is a regular expression to test against the filename. If there
    is a match, the OPTIONS can be used for the file. The regular expression is
    tested using 'grep -e'.
  - OPTIONS
    If the line doesn't contain the FILE\_DESCRIPTOR part, it applies to all source code.

For example, to transform all source code in lower case:
> export OPTS='--lowerCase'; pyft\_tool.py --pyft\_opts\_env OPTS ...

To transform all source code in lower case, except routines in turb directory which must be
in upper case but keeping the turb.F90 in lower case:
> export OPTS='--lowerCase 
> ^turb/:=:--upperCase 
> ^turb/turb\..90:=:--lowerCase'; pyft\_tool.py --pyft\_opts\_env OPTS ...
