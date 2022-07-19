# INSTALLATION NEEDED FOR AROME COMPILATION WITH PHYEX

## ABOUT THIS DOCUMENT

This document is intended for persons who want to compile the AROME model using the PHYEX package.
The installation described here is also needed in order to use the check\_commit\_ial.sh script which enable to check if a given commit reproduces a reference version.

This document is written using the markdown language. With pandoc, it can be converted to HTML (pandoc -s \<filename\>.md -o \<filename\>.html) or PDF (pandoc -V geometry:landscape -s \<filename\>.md -o \<filename\>.pdf).

This document describes how to build reference packs on sxphynh (ubuntu) and belenos.
The following instructions focuse on the 'phyex' pack in which the physics
have been put in a new gmkpack project, named phyex.

Another pack ('main') have been produced as a super-reference without code
move. Only a recompilation has been done.

This document ends with instruction on how to use these packs. This step
is now automatically performed with the prep\_code.sh script.

The same installation guide applies to sxphynh and belenos except for some commands.
The directory in which the repository lies is designated by the TRUNK variable.
The $TRUNK dir can be put on a shared directory to share this installation among several users.

## REFERENCE PACK CREATION

### Create the pack

```
version=01
compiler=MPIGFORTRAN920DBL on ubuntu, MIMPIIFC1805 on belenos
gmkfile=${compiler}.GMAP on ubuntu, ${compiler}.EPONA on belenos
option=xfftw on ubuntu, 2y on belenos
getpack 48t1_main.01.${compiler}.${option} #get source code on ubuntu, is it really necessary?
export GMKTMP=/dev/shm
(. berootpack)
gmkpack -a -r 48t1 -b phyex -n $version -l ${compiler} -o ${option} -p masterodb -h $TRUNK/tools/pack/ #create main pack
```

### Populate main pack with source code

```
cd $TRUNK/tools/pack/48t1_phyex.${version}.${compiler}.${option}/src/local
if sxphynh; then
  wget http://anonymous:mto@webdav.cnrm.meteo.fr/public/algo/khatib/src/48t1_main.01.tgz #only available at MF but equivalent must exist elsewhere
else
  ssh sxphynh.cnrm.meteo.fr "wget http://anonymous:mto@webdav.cnrm.meteo.fr/public/algo/khatib/src/48t1_main.01.tgz -O -" > 48t1_main.01.tgz
fi
tar xf 48t1_main.01.tgz
rm -f 48t1_main.01.tgz
for rep in turb micro conv; do
  mkdir -p phyex/$rep
  mv mpa/$rep/internals/* phyex/$rep/
  mv mpa/$rep/module/* phyex/$rep/
  rmdir mpa/$rep/internals mpa/$rep/module
done
tar xf /cnrm/algo/khatib/drhook.c_for_ubuntu.tar #only on ubuntu
```

### Apply some bug corrections

```
sed -i 's/IF (LBUDGET_RH)/IF (LBUDGET_RH .AND. KRR==7)/' mpa/micro/externals/aro_rain_ice.F90
```

Edition of arpifs/phys\_dmn/apl\_arome.F90 to modift (line 1573):

  ```
  IF (LMFSHAL .AND. (CMF_CLOUD=='DIRE'.OR.CMF_CLOUD=='BIGA')) THEN
    IOFF_MFSHAL=IOFF_MFSHAL+3
    ...
  ENDIF
  ```

into:

  ```
  IF (LMFSHAL .AND. (CMF_CLOUD=='DIRE'.OR.CMF_CLOUD=='BIGA')) THEN
    IOFF_MFSHAL=IOFF_MFSHAL+3
    ...
  ELSE
    DO JLEV = 1, KLEV
      ZRC_MF_(KIDIA:KFDIA,JLEV)=0._JPRB
      ZRI_MF_(KIDIA:KFDIA,JLEV)=0._JPRB
      ZCF_MF_(KIDIA:KFDIA,JLEV)=0._JPRB
    ENDDO
  ENDIF
  ```

Edition of apl\_arome.F90 to modify (line 1406):

  ```
  IF ( LKFBCONV.AND.LOSUBG_COND.AND..NOT.LOSIGMAS) THEN
    DO JLEV = 1, KLEV
      ZMFM_(KIDIA:KFDIA,JLEV)=PSIGM(KIDIA:KFDIA,JLEV)
    ENDDO
  ENDIF
  ```

into:

  ```
  IF (LOSUBG_COND.AND..NOT.LOSIGMAS) THEN
    IF (LKFBCONV) THEN
      DO JLEV = 1, KLEV
        ZMFM_(KIDIA:KFDIA,JLEV)=PSIGM(KIDIA:KFDIA,JLEV)
      ENDDO
    ELSE
      DO JLEV = 1, KLEV
        ZMFM_(KIDIA:KFDIA,JLEV)=0._JPRB
      ENDDO
    ENDIF
  ENDIF
  ```

### Compilation

```
cd $TRUNK/tools/pack/48t1_phyex.${version}.${compiler}.${option}
#Not needed anymore: grep MPA .gmkfile/${compiler}.GMAP | sed 's/MPA/PHYEX/g' >> .gmkfile/${gmkfile}
Edition of .gmkfile/${gmkfile} to add -DREPRO48 to the MACROS_FRT variable in order to suppress bug corrections and be able to reproduce the original cy48
#Not needed anymore: Edition of .gmkfile/${gmkfile} to suppress on ubuntu -ftree-vectorize
sed -i 's/GMK_THREADS=1/GMK_THREADS=10/' ics_masterodb
cleanpack -f
resetpack -f
./ics_masterodb
```

## USER\'S PACK CREATION

```
version=01
compiler=MPIGFORTRAN920DBL on ubuntu, MIMPIIFC1805 on belenos
gmkfile=${compiler}.GMAP on ubuntu, ${compiler}.EPONA on belenos
option=xfftw on ubuntu, 2y on belenos
getpack 48t1_main.01.${compiler}.${option} #get source code on ubuntu, is it really necessary?

commit=9ce8119430dd603d35308d8ae94cf18636157473 #exemple of commit to test against the reference pack

gmkpack -r 48t1 -b phyex -v ${version} -l ${compiler} -o ${option} -p masterodb -f $TRUNK/tools/pack/ -u PHYEX/$commit

cd $HOMEPACK/PHYEX/$commit/src/local/phyex
git clone git@github.com:QuentinRodier/PHYEX.git
cd PHYEX
git checkout $commit
#The exact manipulation to perform depends on the commit to test. For a full description, please see the check\_commit\_ial.sh script
for rep in turb micro conv; do
  mv -f src/common/$rep/* ../$rep/
  mv -f src/arome/$rep/* ../$rep/
  touch ../$rep/*
done
cd ..
rm -rf PHYEX

cd $HOMEPACK/PHYEX/$commit
sed -i 's/GMK_THREADS=1/GMK_THREADS=10/' ics_masterodb
cleanpack -f
./ics_masterodb
```
