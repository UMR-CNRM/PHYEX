#!/bin/bash

#set -x
set -e

# Repertoire où Mesonh MNH-V5-5-0 officiel est installe
#REFDIR=$HOME
# Repertoire où MNH-V5-5-0_PHYEX.tar.gz modifie pour accueillir PHYEX se trouve
#TARGZDIR=$HOME

MNHPACK=${MNHPACK:=$HOME/MesoNH/PHYEX}
PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

function usage {
  echo "Usage: $0 [-h] [-c] [-r] [-C] commit reference"
  echo "commit          commit hash"
  echo "reference       commit hash or nothing for ref"
  echo "-s              suppress compilation pack"
  echo "-c              performs compilation"
  echo "-r              runs the tests, if option appears twice, run is also executed on only 2 procs (instead of 4 procs)"
  echo "-C              checks the result against the reference"
  echo ""
  echo "If nothing is asked (compilation, running, check) everything is done"
}

compilation=0
run=0
check=0
commit=""
reference=""
suppress=0

while [ -n "$1" ]; do
  case "$1" in
    '-h') usage;;
    '-s') suppress=1;;
    '-c') compilation=1;;
    '-r') run=$(($run+1));;
    '-C') check=1;;
    #-b) param="$2"; shift ;;
    #--) shift; break ;;
     *) if [ -z "${commit-}" ]; then
          commit=$1
        else
          if [ -z "${reference-}" ]; then
            reference=$1
          else
            echo "Only two commit hash allowed on command line"
            exit 1
          fi
        fi;;
  esac
  shift
done

if [ $compilation -eq 0 -a \
     $run -eq 0 -a \
     $check -eq 0 ]; then
  compilation=1
  run=1
  check=1
fi

if [ -z "${commit-}" ]; then
  echo "At least one commit hash must be provided on command line"
  exit 2
fi

#if [ $check -eq 1 -a -z "${reference-}" ]; then
#  echo "To perform a comparison two commit hashes are mandatory on the command line"
#  exit 3
#fi

if echo $commit | grep '/' > /dev/null; then
  name=MNH-V5-5-0-$(echo $commit | sed 's/\//_/g')
  [ $suppress -eq 1 -a -d $MNHPACK/$name ] && rm -rf $MNHPACK/$name
else
  name=MNH-V5-5-0-$commit
  [ $suppress -eq 1 -a -d $MNHPACK/$name ] && rm -rf $MNHPACK/$name
fi

#Two possibilities are supported for the simulations
# - they can be done in the the pack we are currently checking
# - they can be done in the reference pack
#They are done in the current pack except if the reference pack
#already contains a tested simulation
run_in_ref=$(ls -d $REFDIR/MNH-V5-5-0/MY_RUN/KTEST/007_16janvier/008_run2_* | tail -1 |wc -l)
if [ $run_in_ref -eq 1 ]; then
  path_user_beg=$REFDIR/MNH-V5-5-0
  path_user_end=_$commit
  path_ref_beg=$REFDIR/MNH-V5-5-0
  if [ "$reference" == "" ]; then
    path_ref_end=
  else
    path_ref_end=_$reference
  fi
else
  path_user_beg=$MNHPACK/$name
  path_user_end=
  path_ref_end=
  if [ "$reference" == "" ]; then
    path_ref_beg=$REFDIR/MNH-V5-5-0
  else
    path_ref_beg=$MNHPACK/MNH-V5-5-0-$reference
  fi
fi

if [ $compilation -eq 1 ]; then
  echo "### Compilation of commit $commit"

  if [ -d $MNHPACK/$name ]; then
    echo "Pack already exists ($MNHPACK/$name), suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
    exit 5
  fi

  # Prepare the pack
  cd $MNHPACK
  cp $TARGZDIR/MNH-V5-5-0_PHYEX.tar.gz .
  tar xvfz MNH-V5-5-0_PHYEX.tar.gz 
  rm MNH-V5-5-0_PHYEX.tar.gz
  mv MNH-V5-5-0 $name
  cd $name/src

  cd $MNHPACK
  echo "Clone repository, and checkout commit $commit"
  git clone https://github.com/QuentinRodier/PHYEX.git
  cd PHYEX
  git checkout $commit
  
  cd src/common/turb
  # Rename all *.F90 to *.f90
  for rep in turb micro aux; do
    cd ../$rep  
    for f in *.F90; do 
      mv -- "$f" "${f%.F90}.f90"
    done
  done
  cd ../../../

  for rep in turb micro conv ext aux; do
    [ ! -d ../$rep ] && mkdir ../$rep
    [ -d src/common/$rep ] && mv -f src/common/$rep/* ../$rep/
    [ -d src/mesonh/$rep ] && mv -f src/mesonh/$rep/* ../$rep/
    touch ../$rep/*
  done
  cd ..
  # Move PHYEX files inside MNH/src/PHYEX
  for rep in turb micro conv aux; do
    mv $rep/* $name/src/PHYEX/$rep/.
    rmdir $rep
  done

  # Move manually ext/ files in src/MNH
  mv -f ext/* $name/src/MNH/. 

  # Clean folder
  rmdir ext
  rm -Rf PHYEX

  cd $name/src/PHYEX/turb
  # Delete files of MNH-V5-5-0/src/MNH and MNH/src/LIB/SURCOUCHE/src with same name
  for rep in turb micro conv aux ; do
    cd ../$rep
    for f in *.f90; do
      echo $f
      rm -f ../../MNH/$f
      rm -f ../../LIB/SURCOUCHE/src/$f
    done
  done
  cd ..
  
  # Delete old files of MNH-V5-5-0/src/MNH that is now called by mode_... NO /aux NEEDED!
  find turb micro conv -name 'mode_*' > remove_non_mode.sh
  sed -i 's/turb\/mode_/rm -f MNH\//g' remove_non_mode.sh
  sed -i 's/micro\/mode_/rm -f MNH\//g' remove_non_mode.sh
  sed -i 's/conv\/mode_/rm -f MNH\//g' remove_non_mode.sh
  chmod +x remove_non_mode.sh
  mv remove_non_mode.sh ../.
  cd ../
  ./remove_non_mode.sh
  # nettoyage, routine non appellee : 
  rm -f MNH/mf_turb_greyzone.f90
  rm -f MNH/compute_frac_ice.f90

  #Configure and compilation
  ./configure
  set +e #file ends with a test that can return false
  . ../conf/profile_mesonh-LXgfortran-R8I4-MNH-V5-5-0-MPIAUTO-DEBUG
  set -e
  make -j 8
  make installmaster
fi

if [ $run -ge 1 ]; then
  echo "### Running of commit $commit"
  echo $commit
  if [ ! -f $MNHPACK/$name/exe/MESONH* ]; then
    echo "Pack does not exist ($MNHPACK/$name) or compilation has failed, please check"
    exit 6
  fi

  if [ $run_in_ref -eq 1 ]; then
    cd $REFDIR/MNH-V5-5-0/MY_RUN/KTEST/007_16janvier/
    [ ! -d 008_run2_$commit ] && cp -R 008_run2 008_run2_$commit
    cd $REFDIR/MNH-V5-5-0/MY_RUN/KTEST/007_16janvier/008_run2_$commit
  else
    cd $MNHPACK/$name/MY_RUN/KTEST/007_16janvier/
    for rep in ???_*; do
      if [ $rep != 008_run2 ]; then
        rm -rf $rep
        ln -s $REFDIR/MNH-V5-5-0/MY_RUN/KTEST/007_16janvier/$rep .
      fi
    done
    [ -d 008_run2 ] && rm -rf 008_run2
    cp -R $REFDIR/MNH-V5-5-0/MY_RUN/KTEST/007_16janvier/008_run2 .
    cd 008_run2
  fi

  set +e #file ends with a test that can return false
  [ $compilation -eq 0 ] && . $MNHPACK/$name/conf/profile_mesonh-LXgfortran-R8I4-MNH-V5-5-0-MPIAUTO-DEBUG
  set -e
  ./clean_mesonh_xyz
  ./run_mesonh_xyz
fi
  
if [ $check -eq 1 ]; then
  allt=0
  echo "### Check commit $commit against commit $reference"
  path_user=$path_user_beg/MY_RUN/KTEST/007_16janvier/008_run2$path_user_end
  path_ref=$path_ref_beg/MY_RUN/KTEST/007_16janvier/008_run2$path_ref_end

  if [ ! -d $path_user ]; then
    echo "$path_user is missing, please run the simulation"
    exit 7
  fi
  if [ ! -d $path_ref ]; then
    echo "$path_ref is missing, please run the reference simulation"
    exit 8
  fi

  echo "Compare with python..."
  # Compare variable of both Synchronous and Diachronic files with printing difference
  set +e
  $PHYEXTOOLSDIR/compare.py $path_user $path_ref
  t=$?
  set -e
  allt=$(($allt+$t))
  
  #Check bit-repro before date of creation of Synchronous file from ncdump of all values (pb with direct .nc file checks)
  file1=$path_user/16JAN.1.12B18.001.nc 
  file2=$path_ref/16JAN.1.12B18.001.nc
  echo "Compare with ncdump..."
  set +e
  diff <(ncdump $file1 | head -c 62889) <(ncdump $file2 | head -c 62889)
  t=$?
  set -e
  allt=$(($allt+$t))

  if [ $allt -eq 0 ]; then
    status="OK"
  else
    status="Files are different"
  fi
  echo "...comparison done: $status"
fi
