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
Tools are designed and tested with TRUNK=\<git repository\>/tools/pack/

## REFERENCE PACK CREATION

### Prerequiste
gmkpack must be installed
fypp python module must be instaled (pip3 install --user fypp)

### Create the pack

```
version=01
cycle=48t1 or cy48t3 (after commit XXX on 22 September 2022)
compiler=MPIGFORTRAN920DBL on ubuntu, MIMPIIFC1805 on belenos
gmkfile=${compiler}.GMAP on ubuntu, ${compiler}.EPONA on belenos
option=xfftw on ubuntu, 2y on belenos
export GMKTMP=/dev/shm
gmkpack -a -r ${cycle} -b phyex -n $version -l ${compiler} -o ${option} -p masterodb -h $TRUNK/tools/pack/ #create main pack
```

### Populate main pack with source code

```
cd $TRUNK${cycle}_phyex.${version}.${compiler}.${option}/src/local
if sxphynh; then
  wget http://anonymous:mto@webdav.cnrm.meteo.fr/public/algo/khatib/src/${cycle}_main.01.tgz #only available at MF but equivalent must exist elsewhere
else
  ssh sxphynh.cnrm.meteo.fr "wget http://anonymous:mto@webdav.cnrm.meteo.fr/public/algo/khatib/src/${cycle}_main.01.tgz -O -" > ${cycle}_main.01.tgz
fi
tar xf ${cycle}_main.01.tgz
rm -f ${cycle}_main.01.tgz
for rep in turb micro conv; do
  mkdir -p phyex/$rep
  mv -f mpa/$rep/internals/* phyex/$rep/
  mv -f mpa/$rep/module/* phyex/$rep/
  rmdir mpa/$rep/internals mpa/$rep/module
done
[ $cycle == 48t1 ] && tar xf /cnrm/algo/khatib/drhook.c_for_ubuntu.tar #only on ubuntu
```

### Apply some bug corrections

```
sed -i 's/IF (LBUDGET_RH)/IF (LBUDGET_RH .AND. KRR==7)/' mpa/micro/externals/aro_rain_ice.F90
```

Edition of arpifs/phys\_dmn/apl\_arome.F90 to modify (line 1573 in 48t1, 1496 in 48t3):

  ```
  IF (LMFSHAL .AND. (CMF_CLOUD=='DIRE'.OR.CMF_CLOUD=='BIGA')) THEN
    IOFF_MFSHAL=IOFF_MFSHAL+3
    ...
  ENDIF
  ```

into (48t1):

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

or (48t3):

  ```
  IF (LMFSHAL .AND. (CMF_CLOUD=='DIRE'.OR.CMF_CLOUD=='BIGA')) THEN
    IOFF_MFSHAL=IOFF_MFSHAL+3
    ...
  ELSE
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      ZRC_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB                         
      ZRI_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB                         
      ZCF_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB                         
    ENDDO
  ENDIF
  ```

Edition of apl\_arome.F90 to modify (line 1406 in 48t1, 1329 in 48t3):

  ```
  IF ( LKFBCONV.AND.LOSUBG_COND.AND..NOT.LOSIGMAS) THEN
    DO JLEV = 1, KLEV
      ZMFM_(...
    ENDDO
  ENDIF
  ```

into:

  ```
  IF (LOSUBG_COND.AND..NOT.LOSIGMAS) THEN
    IF (LKFBCONV) THEN
      DO JLEV = 1, KLEV
        ZMFM_(...
      ENDDO
    ELSE
      DO JLEV = 1, KLEV
        ZMFM_(...)=0._JPRB
      ENDDO
    ENDIF
  ENDIF
  ```

If cycle is 48t3, edition of apl\_arome.F90 to modify (line 3616)

  ```
      ZTMP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=ZFPR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,2)+ZFPR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,3)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP2(:,:),'FQTPRECISTL',YDDDH)
      ZTMP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=ZFPR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,4)+ZFPR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,5)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP2(:,:),'FQTPRECISTN',YDDDH)
  ```

into:

  ```
      ZTMP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,0)=0._JPRB
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        ZTMP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZPFPR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,2)+ZPFPR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,4)
      ENDDO
      !ZTMP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=ZFPR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,2)+ZFPR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,3)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP2(:,:),'FQTPRECISTL',YDDDH)
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        ZTMP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZPFPR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,4)+ZPFPR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,5)
      ENDDO
      !ZTMP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=ZFPR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,4)+ZFPR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,5)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP2(:,:),'FQTPRECISTN',YDDDH)
  ```

If cycle is 48t3, edition of apl\_arome.F90 to:

add the folowing lines after the line 236 (YOMTRAJ use statement):

  ```
  #ifdef REPRO48
  !To compensate a bug introduced in 48t3
  !Must be suppressed as soon as the bug is corrected
  USE MODD_BUDGET
  #endif
  ```

and add the folowing lines at line 912 (to be one of the first execution statement but exact emplacement is not sensitive):

  ```
  #ifdef REPRO48
  !see comment associated to the MODD_BUDGET use statement
  LBU_ENABLE = YDMODEL%YRML_DIAG%YRLDDH%LSDDH
  LBUDGET_U =LBU_ENABLE
  LBUDGET_V =LBU_ENABLE
  LBUDGET_W =LBU_ENABLE
  LBUDGET_TH=LBU_ENABLE
  LBUDGET_TKE=LBU_ENABLE
  LBUDGET_RV=LBU_ENABLE
  LBUDGET_RC=LBU_ENABLE
  LBUDGET_RR =LBU_ENABLE
  LBUDGET_RI =LBU_ENABLE
  LBUDGET_RS =LBU_ENABLE
  LBUDGET_RG =LBU_ENABLE
  LBUDGET_RH =LBU_ENABLE
  LBUDGET_SV=LBU_ENABLE
  #endif
  ```

Edition of phyex/turb/compute\_mf\_cloud\_bigaus.F90 to modify (line 120):

  ```
  DO JK=KKB,KKE,KKL
  ```

into:

  ```
  DO JK=KKB,KKE-KKL,KKL
  ```

### Compilation

```
cd $TRUNK/${cycle}_phyex.${version}.${compiler}.${option}
Edition of .gmkfile/${gmkfile} to add -DREPRO48 to the MACROS_FRT variable in order to suppress bug corrections and be able to reproduce the original cy48
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

commit=9ce8119430dd603d35308d8ae94cf18636157473 #exemple of commit to test against the reference pack

gmkpack -r ${cycle} -b phyex -v ${version} -l ${compiler} -o ${option} -p masterodb -f $TRUNK -u PHYEX/$commit

cd $HOMEPACK/PHYEX/$commit/src/local/phyex
git clone git@github.com:UMR-CNRM/PHYEX.git
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
