#!/bin/bash

set -x
set -e

#This script:
# - compiles the PHYEX package using a specific commit
# - runs the different test progs and checks if results are identical to a given version

#ice_adjust: the ice adjust test case

specialName="ref"
availTests="ice_adjust,rain_ice"
defaultTest='ALL'
separator='_' #- seprator must be in sync with prep_code.sh separator

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

dirdata=$PHYEXTOOLSDIR/testprogs_data
if [ $(hostname | cut -c 1-7) == 'belenos' -o $(hostname | cut -c 1-7) == 'taranis' ]; then
  HPC=1
  archfile=MIMPIIFC1805.EPONA
else
  HPC=0
  archfile=gnu
fi
defaultRef=ref

function usage {
  echo "Usage: $0 [-h] [-c] [-r] [-C] [-s] [-f] [--noexpand] [-t test] commit reference"
  echo "commit          commit hash (or a directory, or among $specialName) to test"
  echo "reference       commit hash (or a directory, or among $specialName) REF to use as a reference"
  echo "-s              suppress compilation pack"
  echo "-c              performs compilation"
  echo "-r              runs the tests"
  echo "-C              checks the result against the reference"
  echo "-t              comma separated list of tests to execute"
  echo "                or ALL to execute all tests"
  echo "--noexpand      do not use mnh_expand (code will be in array-syntax)"
  echo ""
  echo "If nothing is asked (compilation, running, check) everything is done"
  echo
  echo "With the special reference REF commit, a suitable reference is guessed"
  echo
  echo "If no test is aked for, the default one ($defaultTest) is executed"
  echo
  echo "The directory (for commit only, not ref) can take the form server:directory"
  echo
  echo "If using a directory (for commit or reference) it must contain at least one '/'"
}

compilation=0
run=0
check=0
commit=""
reference=""
tests=""
suppress=0
useexpand=1

while [ -n "$1" ]; do
  case "$1" in
    '-h') usage;;
    '-s') suppress=1;;
    '-c') compilation=1;;
    '-r') run=$(($run+1));;
    '-C') check=1;;
    '-t') tests="$2"; shift;;
    '--noexpand') useexpand=0;;
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

TESTDIR=${TESTPROGSDIR:=$HOME/TESTPROGS}

function exescript () {
  #usage: exescript <output file> <script> [arg [arg ...]]
  output=$1
  shift
  if [ $HPC -eq 1 ]; then
    sbatch --wait -o $output $@
    cat $output
  else
    $@ 2>&1 | tee $output
  fi
}

if [ -z "${tests-}" ]; then
  tests=$defaultTest
fi
if [ $tests == 'ALL' ]; then
  tests=$availTests
fi

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

if [ $check -eq 1 -a -z "${reference-}" ]; then
  echo "To perform a comparison two commit hashes are mandatory on the command line"
  exit 3
fi

#Name is choosen such as it can be produced with a main pack: PHYEX/48t1_XXXXXXXXX.01.${gmkpack_l}.${gmkpack_o}
fromdir=''
if echo $commit | grep '/' > /dev/null; then
  fromdir=$commit
  name=$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
  [ $suppress -eq 1 -a -d $TESTDIR/$name ] && rm -rf $TESTDIR/$name
elif echo $specialName | grep -w $commit > /dev/null; then
  name="$commit"
else
  name="COMMIT$commit"
  [ $suppress -eq 1 -a -d $TESTDIR/$name ] && rm -rf $TESTDIR/$name
fi
if [ ! -z "${reference-}" ]; then
  [ $reference == 'REF' ] && reference=$defaultRef
  reffromdir=''
  if echo $reference | grep '/' > /dev/null; then
    reffromdir=$reference
    refname=$(echo $reference | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
  elif echo $specialName | grep -w $reference > /dev/null; then
    refname="$reference"
  else
    refname="COMMIT${reference}"
  fi
fi

if [ $compilation -eq 1 ]; then
  echo "### Compilation of commit $commit"

  if echo $specialName | grep -w $commit > /dev/null; then
    echo "Special commit '$commit' cannot be compiled with this script"
    exit 4
  fi

  if [ -d $TESTDIR/$name ]; then
    echo "Pack already exists ($TESTDIR/$name), suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
    exit 5
  fi
  mkdir $TESTDIR/$name
  cd $TESTDIR/$name/

  MNH_EXPAND_DIR=$PHYEXTOOLSDIR/mnh_expand
  export PATH=$MNH_EXPAND_DIR/filepp:$MNH_EXPAND_DIR/MNH_Expand_Array:$PATH

  if [ $useexpand == 1 ]; then
    expand_options="-D MNH_EXPAND -D MNH_EXPAND_LOOP"
  else
    expand_options=""
  fi
  subs="$subs -s turb -s micro -s aux -s ice_adjust"
  prep_code=$PHYEXTOOLSDIR/prep_code.sh
  if [ "$fromdir" == '' ]; then
    echo "Clone repository, and checkout commit $commit (using prep_code.sh)"
    $prep_code -c $commit $expand_options $subs src
  else
    echo "Copy $fromdir"
    scp -q -r $fromdir src
    $prep_code $expand_options $subs src
  fi

  cd $TESTDIR/$name/src/build/with_fcm/
  rm -rf arch_*
  ./make_fcm.sh --arch $archfile 2>&1 | tee Output_compilation
fi

if [ $run -ge 1 ]; then
  echo "### Running of commit $commit"

  for t in $(echo $tests | sed 's/,/ /g'); do
    if [ ! -f $TESTDIR/$name/src/build/with_fcm/arch_${archfile}/build/bin/main_${t}.exe ]; then
      echo "Pack does not exist ($TESTDIR/$name) or compilation has failed, please check"
      exit 6
    fi
  done

  #Cleaning to suppress old results that may be confusing in case of a crash during the run
  for t in $(echo $tests | sed 's/,/ /g'); do
    cd $TESTDIR/$name
    if [ -d tests/$t ]; then
      rm -rf tests/$t
    fi
  done

  #Run the tests one after the other
  for t in $(echo $tests | sed 's/,/ /g'); do
    cd $TESTDIR/$name
    mkdir -p tests/$t
    cd tests/$t
    ln -s $dirdata/$t data
    $TESTDIR/$name/src/build/with_fcm/arch_${archfile}/build/bin/main_${t}.exe --check 2>&1 > Output_run
  done
fi

if [ $check -eq 1 ]; then
  echo "### Check commit $commit against commit $reference"

  alltests=0
  message=""
  for t in $(echo $tests | sed 's/,/ /g'); do
    file1=$TESTDIR/$name/tests/$t/Output_run
    file2=$TESTDIR/$refname/tests/$t/Output_run
    mess=""
    te=0
    if [ ! -f "$file1" ]; then
      mess="Result ($file1) for commit $commit does not exist, please run the simulation"
      te=1
    fi
    if [ ! -f "$file2" ]; then
      mess2="Result ($file2) for commit $reference does not exist, please run the simulation"
      te=1
      if [ "$mess" = "" ]; then
        mess=$mess2
      else
        mess="$mess and $mess2"
      fi
    fi
    if [ $te -eq 0 ]; then
      set +e
      mess=$(cmp $file1 $file2 246 246 2>&1)
      te=$?
      set -e
    fi
    [ $te -ne 0 ] && message="$message $mess \n"
    alltests=$(($alltests+$te))
  done
  if [ $alltests -eq 0 ]; then
    echo "SUCCESS, files are identical"
  else
    echo "*************** Files are different *******************"
    echo -e "$message"
  fi
fi
