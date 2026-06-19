#!/bin/bash

#This script:
# - compiles the MESONH model using a specific commit for the externalised physics
# - runs tests and checks if results are identical to a given version

#The folowing environment variables can be defined:
# MNHPACK: directory where tests are build


PHYEXTOOLSDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Mutualised functions and definitions
. ${PHYEXTOOLSDIR}/check_commit_common.sh

#######################
#### CONFIGURATION ####
#######################

#About the tests:
# - fileFromCase is an array defining allowed tests which can depend on platform, if we ask to perform an action
#   with a test not in the fileFromCase array, the action is ignored. First argument is the output file of the
#   run step and the second one is an optional diachronic file resulting from the run step.
# - defaultTest is the list of tests to perform when no '-t' option is provided on the command line.
declare -A fileFromCase=(
  ["KTEST/007_16janvier"]="008_run2/16JAN.1.12B18.001.nc 008_run2/16JAN.1.12B18.000.nc"
  ["INTEGRATION_CASES/LOCAL/ARMCU_1D_CONDSAMP"]="002_mesonh/ARM__.1.CEN4T.001.nc 002_mesonh/ARM__.1.CEN4T.001.nc"
  ["INTEGRATION_CASES/HPC/ARMCU_LES/DEAR"]="ARM__.1.CEN4T.001.nc ARM__.1.CEN4T.000.nc"
  ["KTEST/014_LIMA"]="002_mesonh/XPREF.1.SEG01.002.nc 002_mesonh/XPREF.1.SEG01.000.nc"
  ["INTEGRATION_CASES/HPC/OCEAN_LES"]="004_run2/SPWAN.2.25m00.001.nc"
)
allowedTests=$(IFS=','; echo "${!fileFromCase[*]}")
defaultTest="KTEST/007_16janvier"

export MNHPACK=${MNHPACK:=$HOME/MesoNH/PHYEX}

################################
#### COMMAND LINE ARGUMENTS ####
################################

default_expand=false
enable_prepCodeOpts=true
command_line $@

#######################
#### FROM A COMMIT ####
#######################

# In this case, we clone the commit and run the check_commit script
# of this commit using the cloned directory.
hash_prefix=''
tools_install=''
clone_and_run

###########################
#### COMMIT PROPERTIES ####
###########################

if [ -d $commit/src ]; then
  model_ready=false
  json_content=$(scp $commit/src/mesonh/mesonh_version.json /dev/stdout 2>/dev/null || echo "")
else
  model_ready=true
  json_content=$(scp $commit/mesonh_version.json /dev/stdout 2>/dev/null || echo "")          
fi

declare -A refByTest
get_properties

refversion=$(json_dictkey2value "$json_content" 'refversion' '')
mnhdir=${refversion}-$name

#######################
#### PACK CREATION ####
#######################

[ $suppress == true -a -d $MNHPACK/$mnhdir ] && rm -rf $MNHPACK/$mnhdir
if [ $packcreation == true -a -d $MNHPACK/$mnhdir -a $onlyIfNeeded == true ]; then
  packcreation=false
fi
if [ $packcreation == true ]; then
  if [ -d $MNHPACK/$mnhdir ]; then
    echo "Pack already exists ($MNHPACK/$mnhdir),"
    echo "suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
    exit 5
  else
    echo "### Pack creation for commit $commit"

    # Prepare the pack (we use tar.gz instead of 'git clone' to avoid depending on 'git lfs')
    cd $MNHPACK
    refcommit=$(json_dictkey2value "$json_content" 'MESONHcommit' 'XXXXXXX')
    url=$(json_dictkey2value "$json_content" 'MESONHrepo' 'https://src.koda.cnrs.fr/mesonh/mesonh-code')
    url="${url}/-/archive/${refcommit}/mesonh-code-${refcommit}.tar.gz"
    wget --no-check-certificate ${url} -O mesonh-code-${refcommit}.tar.gz
    tar xf mesonh-code-${refcommit}.tar.gz
    rm -f mesonh-code-${refcommit}.tar.gz
    mv mesonh-code-${refcommit} ${mnhdir}
  fi
fi
if [ $packupdate == true -o $packcreation == true ]; then
  if [ !  -d $MNHPACK/$mnhdir/src/PHYEX ]; then
      echo "PHYEX directory doesn't exist in pack ($MNHPACK/$mnhdir)"
      exit 9
  else
    cd $MNHPACK/$mnhdir/src
    if [ $packupdate == true ]; then
      mv PHYEX PHYEXori
    else
      rm -rf PHYEX
    fi

    if [ $useexpand == true ]; then
      expand_options="--mnhExpand"
    else
      expand_options=""
    fi
    prep_code=$PHYEXTOOLSDIR/prep_code.sh
    echo "Copy $commit"
    mkdir PHYEX
    if [ $model_ready == false ]; then
      scp -q -r $commit/src PHYEX/
      subs="-s turb -s micro -s aux -s ext -s conv"
      $prep_code $prepCodeOpts --renameFf --ilooprm --noRaiseOnCodingNorms $expand_options $subs -m mesonh PHYEX -- --removeExtraDOinMnhDoConcurrent
    else
      $prep_code $prepCodeOpts PHYEX #Ready for inclusion
      rm -rf PHYEX/.git
    fi
    find PHYEX -type f -exec touch {} \; #to be sure a recompilation occurs

    if [ $packupdate == true ]; then
      #Update only modified files
      cd PHYEX
      for file in $(find turb micro conv aux -type f); do
        mvdiff $file ../PHYEXori/$file
      done
      cd ..
      rm -rf PHYEX
      mv PHYEXori PHYEX
    fi

    # Remove binaries
    rm -f $MNHPACK/$mnhdir/exe/*
  fi
fi

#####################
#### COMPILATION ####
#####################

profile_sourced=0
if [ $compilation == true ]; then
  if [ $onlyIfNeeded == false -o ! -f $MNHPACK/$mnhdir/exe/MESONH* ]; then
    echo "### Compilation of commit $commit"
    cd $MNHPACK/$mnhdir/src
    #Configure and compilation
    command -v module && modulelist=$(module -t list 2>&1 | tail -n +2) #save loaded modules
    ./configure
    set +e #file ends with a test that can return false
    . ../conf/profile_mesonh-* #This lines modifies the list of loaded modules
    set -e
    profile_sourced=1
    rm -f ../exe/* #Suppress old executables, if any
    make -j 8 2>&1 | tee ../Output_compilation
    make installmaster 2>&1 | tee -a ../Output_compilation
    command -v module && module load $modulelist #restore loaded modules
  fi
fi

###################
#### EXECUTION ####
###################

if [ $run == true ]; then
  #Cleaning to suppress old results that may be confusing in case of a crash during the run
  if [ $onlyIfNeeded == false ]; then
    for t in $(echo $tests | sed 's/,/ /g'); do
      (cd $MNHPACK/$mnhdir/MY_RUN/$t; make clean)
    done
  fi

  #Run the tests one after the other
  firstrun=1
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then #test is allowed on this plateform
      casedir=$MNHPACK/$mnhdir/MY_RUN/$t
      read -r file1 file2 <<< "${fileFromCase[$t]}"
      if [ ! -d $casedir/$file1 ]; then #We do not enter systematically this part if onlyIfNeeded=true
        if [ $firstrun -eq 1 ]; then
          echo "### Running of commit $commit"
          firstrun=0
        fi

        if [ ! -f $MNHPACK/$mnhdir/exe/MESONH* ]; then
          echo "Pack does not exist ($MNHPACK/$mnhdir) or compilation has failed, please check"
          exit 6
        fi

        #execution
        cd ${casedir}
        if [ $profile_sourced -eq 0 ]; then
          set +e #file ends with a test that can return false
          . $MNHPACK/$mnhdir/conf/profile_mesonh-*
          set -e
          profile_sourced=1
        fi
        set +o pipefail #We want to go through all tests
        t1=$(($(date +%s%N)/1000)) #current time in milliseconds
        export POSTRUN=echo #to disable plotting
        export PHYEX_REDUCED_GRID=yes #to reduce the grid size for some tests
        #echo yes to accept the dowloading of PGD files
        echo yes | make all 2>&1 | tee Output_run
        t2=$(($(date +%s%N)/1000))
        set -o pipefail
        if [ "$perffile" != "" ]; then
          echo "$commit mesonh $t $(($t2-$t1))" >> "$perffile"
        fi
      fi
    fi
  done
fi

####################
#### COMPARISON ####
####################

if [ $check == true ]; then
  echo "### Check commit $commit against commit $reference"

  allt=0
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then
      #Files to compare
      caseref=${refByTest[$t]}
      if echo $caseref | grep '/' > /dev/null; then
        refname=${refversion}-$(escape_commit $caseref)
      else
        refname=${refversion}-$caseref
      fi
      read -r file1 file2 <<< "${fileFromCase[$t]}"
      path_user=$MNHPACK/$mnhdir/MY_RUN/$t
      path_ref=$MNHPACK/$refname/MY_RUN/$t
      file1u=$path_user/$file1
      file1r=$path_ref/$file1
      if [ "$file2" != "" ]; then
        file2u=$path_user/$file2
        file2r=$path_ref/$file2
      else
        file2u=""
        file2r=""
      fi

      if [ ! -d $path_user ]; then
        echo "$path_user is missing, please run the simulation"
        exit 7
      fi

      #Run the reference if needed
      if [ ! -f $file1r -a $computeRefIfNeeded == true ]; then
        $0 -p -c -r -t $t --onlyIfNeeded ${caseref}
      fi

      if [ ! -d $path_ref ]; then
        echo "$path_ref is missing, please run the reference simulation"
        exit 8
      fi

      #Comparison
      if [ -f $file1u -a -f $file1r ]; then
        # Compare variable of both Synchronous and Diachronic files with printing difference
        echo "Comparison for case $t..."
        set +e
        if [ "$file2u" == "" ]; then
          $PHYEXTOOLSDIR/compare.py --backup $file1u $file1r
          r=$?
        else
          $PHYEXTOOLSDIR/compare.py --backup $file1u $file1r --diac $file2u $file2r
          r=$?
        fi
        set -e
        allt=$(($allt+$r))

        #Check bit-repro except date of creation of Synchronous file from ncdump of all values
        #(pb with direct .nc file checks)
        set +e
        $PHYEXTOOLSDIR/compare.py --ncdump $file1u $file1r
        r=$?
        set -e
        allt=$(($allt+$r))
      else
        [ ! -f $file1u ] && echo "  $file1u is missing"
        [ ! -f $file1r ] && echo "  $file1r is missing"
        allt=$(($allt+1))
      fi
    fi
  done

  if [ $allt -eq 0 ]; then
    status="OK"
  else
    status="Files are different"
    cmpstatus=50
  fi
  echo "...comparison done: $status"
fi

##################
#### CLEANING ####
##################

if [ $remove == true ]; then
  echo "### Remove model directory for commit $commit"
  [ -d $MNHPACK/$mnhdir ] && rm -rf $MNHPACK/$mnhdir
fi

exit $cmpstatus
