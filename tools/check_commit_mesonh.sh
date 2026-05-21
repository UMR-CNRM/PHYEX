#!/bin/bash

#This script:
# - compiles the MESONH model using a specific commit for the externalised physics
# - runs tests and checks if results are identical to a given version

#The folowing environment variables can be defined:
# TARGZDIR: directory where tar.gz files are searched for
# MNHPACK: directory where tests are build


PHYEXTOOLSDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Mutualised functions and definitions
. ${PHYEXTOOLSDIR}/check_commit_common.sh

#######################
#### CONFIGURATION ####
#######################

#About the tests:
# - allowedTests is the list of allowed tests which can depend on platform, if we ask to perform an action
#   with a test not in the allowedTests list, the action is ignored
# - defaultTest is the list of tests to perform when no '-t' option is provided on the command line.
allowedTests="007_16janvier/008_run2, 007_16janvier/008_run2_turb3D, 007_16janvier/008_run2_lredf, 
              COLD_BUBBLE/002_mesonh, ARMLES/RUN, COLD_BUBBLE_3D/002_mesonh, OCEAN_LES/004_run2,
              014_LIMA/002_mesonh"
defaultTest="007_16janvier/008_run2"

export MNHPACK=${MNHPACK:=$HOME/MesoNH/PHYEX}
export TARGZDIR=${TARGZDIR:=$PHYEXTOOLSDIR/pack/}

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

    # Prepare the pack
    cd $MNHPACK
    mkdir ${mnhdir}_$$
    cd ${mnhdir}_$$
    cp $TARGZDIR/${refversion}.tar.gz .
    tar xfz ${refversion}.tar.gz 
    rm ${refversion}.tar.gz
    mv ${refversion} ../$mnhdir
    cd ..
    rmdir ${mnhdir}_$$
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

    # Move manually ext/ files in src/MNH
    if [ -d PHYEX/ext ]; then
      for file in PHYEX/ext/*; do
        [ $file = PHYEX/ext/modd_salt.f90 ] && mvdiff $file ACLIB/aux/
        [ -f $file ] && mvdiff $file MNH/ # Not an ACLIB file
      done
      if [ $packupdate == true ]; then
        #Only modified files are moved
        rm -rf PHYEX/ext
      else
        #All files must have been moved
        rmdir PHYEX/ext
      fi
    fi

    if [ $packupdate == true ]; then
      #Update only modified files
      cd PHYEX
      for file in $(find turb micro conv aux -type f); do
        mvdiff $file ../PHYEXori/$file
      done
      cd ..
      rm -rf PHYEX
      mv PHYEXori PHYEX
    else
      cd $MNHPACK/$mnhdir/src/PHYEX/turb
      # Delete files of ${refversion}/src/MNH and MNH/src/LIB/SURCOUCHE/src with same name
      for rep in turb micro conv aux ; do
        cd ../$rep
        for f in *.f90; do
          echo $f
          rm -f ../../MNH/$f
          rm -f ../../LIB/SURCOUCHE/src/$f
        done
      done
      cd ..
    
      # Delete old files of ${refversion}/src/MNH that is now called by mode_... NO /aux NEEDED!
      find turb micro conv -name 'mode_*' > remove_non_mode.sh
      sed -i 's/turb\/mode_/rm -f MNH\//g' remove_non_mode.sh
      sed -i 's/micro\/mode_/rm -f MNH\//g' remove_non_mode.sh
      sed -i 's/conv\/mode_/rm -f MNH\//g' remove_non_mode.sh
      chmod +x remove_non_mode.sh
      mv remove_non_mode.sh ../.
      cd ../
      ./remove_non_mode.sh
    fi

    # Remove binaries
    rm -f $MNHPACK/$mnhdir/exe/*

    # Remove execution results
    for t in $(echo $allowedTests | sed 's/,/ /g'); do
      case1=$(echo $t | cut -d / -f 1)
      case2=$(echo $t | cut -d / -f 2)
      casedir=$MNHPACK/$mnhdir/MY_RUN/KTEST/$case1
      [ -d $casedir/$case2 ] && rm -rf $casedir/$case2
    done
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
      case1=$(echo $t | cut -d / -f 1)
      case2=$(echo $t | cut -d / -f 2)
      casedir=$MNHPACK/$mnhdir/MY_RUN/KTEST/$case1
      [ -d $casedir/$case2 ] && rm -rf $casedir/$case2
    done
  fi

  #Run the tests one after the other
  firstrun=1
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then #test is allowed on this plateform
      case1=$(echo $t | cut -d / -f 1)
      case2=$(echo $t | cut -d / -f 2)
      casedir=$MNHPACK/$mnhdir/MY_RUN/KTEST/$case1
      if [ ! -d $casedir/$case2 ]; then #We do not enter systematically this part if onlyIfNeeded=true
        if [ $firstrun -eq 1 ]; then
          echo "### Running of commit $commit"
          firstrun=0
        fi

        if [ ! -f $MNHPACK/$mnhdir/exe/MESONH* ]; then
          echo "Pack does not exist ($MNHPACK/$mnhdir) or compilation has failed, please check"
          exit 6
        fi

        #If the test case didn't exist in the tar.gz, we copy it from from the reference version
        #and we suppress all the test directories for this case
        if [ ! -d $casedir ]; then
          cp -r $MNHPACK/${refversion}/MY_RUN/KTEST/$case1 $casedir/
          for newt in $(echo $allowedTests | sed 's/,/ /g'); do
            newcase1=$(echo $newt | cut -d / -f 1)
            newcase2=$(echo $newt | cut -d / -f 2)
            if [ $case1 == $newcase1 ]; then
              [ -d $casedir/$newcase2 ] && rm -rf $casedir/$newcase2
            fi
          done
        fi

        #Loop on the subdirectories to replace them by links to their reference version
        cd $casedir
        for d in *; do
          if [[ -d "$d" || ( -L "$d" && ! -e "$d" ) ]]; then #directory (or a link to a directory) or a broken link
            if ! echo $allowedTests | grep ${case1}/$d > /dev/null; then
              #This directory is not a test case but might be needed to run the test case,
              #we take the reference version
              rm -rf $d
              ln -s $MNHPACK/${refversion}/MY_RUN/KTEST/$case1/$d 
            fi
          fi
        done

        #We create the test case directory
        cp -r $MNHPACK/${refversion}/MY_RUN/KTEST/$case1/${case2} .

        #execution
        cd ${case2}
        if [ $profile_sourced -eq 0 ]; then
          set +e #file ends with a test that can return false
          . $MNHPACK/$mnhdir/conf/profile_mesonh-*
          set -e
          profile_sourced=1
        fi
        ./clean_mesonh_xyz
        set +o pipefail #We want to go through all tests
        t1=$(($(date +%s%N)/1000)) #current time in milliseconds
        ./run_mesonh_xyz | tee Output_run
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
      case1=$(echo $t | cut -d / -f 1)
      case2=$(echo $t | cut -d / -f 2)
      path_user=$MNHPACK/$mnhdir/MY_RUN/KTEST/$case1/$case2
      path_ref=$MNHPACK/$refname/MY_RUN/KTEST/$case1/$case2
      file3=""
      file4=""
      if [ $case1 == 007_16janvier ]; then
        file1=$path_user/16JAN.1.12B18.001.nc 
        file2=$path_ref/16JAN.1.12B18.001.nc
        file3=$path_user/16JAN.1.12B18.000.nc 
        file4=$path_ref/16JAN.1.12B18.000.nc
        bit_diff=40500
      elif [ $case1 == COLD_BUBBLE ]; then
        file1=$path_user/BUBBL.1.CEN4T.001.nc
        file2=$path_ref/BUBBL.1.CEN4T.001.nc
        bit_diff=27300
      elif [ $case1 == OCEAN_LES ]; then
        file1=$path_user/SPWAN.2.25m00.001.nc
        file2=$path_ref/SPWAN.2.25m00.001.nc
        bit_diff=18400
      elif [ $case1 == COLD_BUBBLE_3D ]; then
        file1=$path_user/BUBBL.1.CEN4T.001.nc
        file2=$path_ref/BUBBL.1.CEN4T.001.nc
        file3=$path_user/BUBBL.1.CEN4T.000.nc
        file4=$path_ref/BUBBL.1.CEN4T.000.nc
        bit_diff=27300
      elif [ $case1 == ARMLES ]; then
        file1=$path_user/ARM__.1.CEN4T.001.nc
        file2=$path_ref/ARM__.1.CEN4T.001.nc
        file3=$path_user/ARM__.1.CEN4T.000.nc
        file4=$path_ref/ARM__.1.CEN4T.000.nc
        bit_diff=71100
      elif [ $case1 == 014_LIMA ]; then
        file1=$path_user/XPREF.1.SEG01.002.nc
        file2=$path_ref/XPREF.1.SEG01.002.nc
        file3=$path_user/XPREF.1.SEG01.000.nc
        file4=$path_ref/XPREF.1.SEG01.000.nc
        bit_diff=32200
      else
        echo "cas $t non reconnu"
      fi

      if [ ! -d $path_user ]; then
        echo "$path_user is missing, please run the simulation"
        exit 7
      fi

      #Run the reference if needed
      if [ ! -f $file2 -a $computeRefIfNeeded == true ]; then
          #We must call it in another shell because of the potentially loaded MesoNH profile
          #because we cannot load two MesoNH profiles in the same shell
          env -i $SHELL -l -c "MNHPACK=${MNHPACK} TARGZDIR=${TARGZDIR} \
                               PHYEXREPOuser=${PHYEXREPOuser} PHYEXREPOprotocol=${PHYEXREPOprotocol} \
                               $0 -p -c -r -t $t --onlyIfNeeded ${caseref}"
      fi

      if [ ! -d $path_ref ]; then
        echo "$path_ref is missing, please run the reference simulation"
        exit 8
      fi

      #Comparison
      if [ -f $file1 -a -f $file2 ]; then
        # Compare variable of both Synchronous and Diachronic files with printing difference
        echo "Comparison for case $t..."
        set +e
        if [ "$file3" == "" ]; then
          $PHYEXTOOLSDIR/compare.py --backup $file1 $file2
          r=$?
        else
          $PHYEXTOOLSDIR/compare.py --backup $file1 $file2 --diac $file3 $file4
          r=$?
        fi
        set -e
        allt=$(($allt+$r))

        #Check bit-repro before date of creation of Synchronous file from ncdump of all values
        #(pb with direct .nc file checks)
        set +e
        $PHYEXTOOLSDIR/compare.py --ncdump $file1 $file2 $bit_diff
        r=$?
        set -e
        allt=$(($allt+$r))
      else
        [ ! -f $file1 ] && echo "  $file1 is missing"
        [ ! -f $file2 ] && echo "  $file2 is missing"
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
