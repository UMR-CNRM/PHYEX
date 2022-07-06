# INSTALLATION NEEDED FOR MESONH COMPILATION WITH PHYEX

## ABOUT THIS DOCUMENT 

This document is intended for persons who want to compile the MESONH model using the PHYEX package.
The installation described here is also needed in order to use the check\_commit\_mesonh.sh script which enable to check if a given commit reproduces a reference version.

This document is written using the markdown language. With pandoc, it can be converted to HTML (pandoc -s \<filename\>.md -o \<filename\>.html) or PDF (pandoc -V geometry:landscape -s \<filename\>.md -o \<filename\>.pdf).

The directory in which the repository lies is designated by the TRUNK variable.
The $TRUNK dir can be put on a shared directory to share this installation among several users.

## COMPILATION OF THE MASTER

An official version of Meso-Nh is installed in the pack directory:

```
cd $TRUNK/tools/pack
scp sxphynh.cnrm.meteo.fr:/home/rodierq/MNH-V5-5-0.tar.gz . #For MF users, can be retrieve on the Meso-NH website
tar xvfz MNH-V5-5-0.tar.gz
cd MNH-V5-5-0/src
./configure
. ../conf/profile_mesonh-*
make -j 8
make installmaster
```

## PREPROCESSING STEPS FOR THE TEST CASE

The preprocessing steps must be done at least once on the master pack:

```
cd ../MY_RUN/KTEST/007_16janvier
#The official namelist can be modified to enable more options
rm -Rf 008_run2; scp -r sxphynh.cnrm.meteo.fr:/home/rodierq/MNH-V5-5-0/MY_RUN/KTEST/007_16janvier/008_run2 .

make clean
make #the step #10 succeed only if a isplay (X11) is available
     #after the error, the next steps can be executed by "make E011_ncl E012_spectre"
```

## GET A MODIFIED PACK SUITABLE FOR PHYEX

The master pack cannot be used directly to compile Meso-NH using PHYEX, the Makefile must be modified
```
cd $TRUNK/pack
scp sxphynh.cnrm.meteo.fr:/home/rodierq/MNH-V5-5-0_PHYEX.tar.gz .
```

## MESONH COMPILATION WITH PHYEX

To compile MESONH using the PHYEX package, some manipulation is needed.
The easiest way is to use the check\_commit\_mesonh.sh script.

## COMPARISON OF SOME COMMITS
Some commits doesn't reproduce the reference commit but are comparable to the b1e20 commit.
If someone is interested in compiling this commit, the argument order of the functions DZM\_MF, MZM\_MF, GZ\_M\_W\_MF must be reversed in compute\_bl89\_ml, tridiag\_massflux.f90 and shuman\_mf.f90
