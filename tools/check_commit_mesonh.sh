#!/bin/bash

#set -x
set -e
set -o pipefail #abort if left command on a pipe fails

#Notes for v5.5.0
#For the OCEAN_LES/004_run2 case, results obtained are different from those obtained with the original version
#of Meso-NH because of new developments and bug correction. The reference version is given by commit e053c59.
#In this commit two modifications must be done in turb/mode_tke_eps_sources.f90 to change twice LOCEAN into OOCEAN.

#######################
#### CONFIGURATION ####
#######################

#The folowing environment variables can be defined:
# TARGZDIR: directory where tar.gz files are searched for
# MNHPACK: directory where tests are build

#About the tests:
# - ALLTests is a list of tests to be done when '-t ALL' is used. This list is filled here
#   in case there is no mesonh_version.json file containig a 'testing' section. If this 'testing'
#   section exists, this list is overridden.
# - allowedTests is the list of allowed tests which can depend on platform, if we ask to perform an action
#   with a test not in the allowedTests list, the action is ignored
# - defaultTest is the list of tests to perform when no '-t' option is provided on the command line.
ALLTests="007_16janvier/008_run2, 007_16janvier/008_run2_turb3D, 007_16janvier/008_run2_lredf, 
          COLD_BUBBLE/002_mesonh, ARMLES/RUN, COLD_BUBBLE_3D/002_mesonh,OCEAN_LES/004_run2,014_LIMA/002_mesonh"
defaultTest="007_16janvier/008_run2"
allowedTests=$ALLTests

separator='_' #- be carrefull, gmkpack (at least on belenos) has multiple allergies (':', '.', '@')
              #- seprator must be in sync with prep_code.sh separator

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

MNHPACK=${MNHPACK:=$HOME/MesoNH/PHYEX}
TARGZDIR=${TARGZDIR:=$PHYEXTOOLSDIR/pack/}

################################
#### COMMAND LINE ARGUMENTS ####
################################

function usage {
  echo "Usage: $0 [-h] [-p] [-u] [-c] [-r] [-C] [-s] [--expand] [-t TEST] [--repo-user USER] [--repo-protocol PROTOCOL] [--remove] [--onlyIfNeeded] [--computeRefIfNeeded] [--prep_code-opts 'OPTS'] [--per FILE] commit [reference]"
  echo "commit          commit hash (or a directory) to test"
  echo "reference       commit hash or a directory or nothing for ref"
  echo "-s              suppress compilation pack"
  echo "-p              creates pack"
  echo "-u              updates pack (experimental, only for 'small' updates where"
  echo "                              the list of files does not change)"
  echo "-c              performs compilation"
  echo "-r              runs the tests"
  echo "-C              checks the result against the reference"
  echo "-t TEST         comma separated list of tests to execute"
  echo "                or ALL to execute all tests"
  echo "--expand        expand mnh_expand blocks (code will use do loops)"
  echo "--repo-user USER"
  echo "                user hosting the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOuser (=$PHYEXREOuser)"
  echo "--repo-protocol PROTOCOL"
  echo "                protocol (https or ssh) to reach the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOprotocol (=$PHYEXREOprotocol)"
  echo "--remove        removes the pack"
  echo "--onlyIfNeeded  performs the pack creation and/or the compilation and/or the execution"
  echo "                only if the step has not already been done"
  echo "--computeRefIfNeeded"
  echo "                computes the missing references"
  echo "--prep_code-opts 'OPTS'"
  echo "                OPTS is added to the call to prep_code (e.g. --prep_code_opts '--lowerCase'"
  echo "                to transfor all source codes in lower case). Help on prep_code.sh options"
  echo "                can be found with 'prep_code.sh -h'. Note: don't forget to enclose OPTS in ' or \""
  echo "--perf FILE     add performance statistics in file FILE"
  echo ""
  echo "If nothing is asked (pack creation, compilation, running, check, removing) everything"
  echo "except the removing is done"
  echo
  echo "If no test is aked for, the default one ($defaultTest) is executed"
  echo
  echo "With the special reference REF commit, a suitable reference is guessed"
  echo
  echo "The directory (for commit only, not ref) can take the form server:directory"
  echo
  echo "If using a directory (for commit or reference) it must contain at least one '/'"
  echo "The commit can be a tag, written with syntagx tags/<TAG>"
}

packcreation=0
packupdate=0
compilation=0
run=0
check=0
commit=""
reference=""
tests=""
suppress=0
useexpand=0
remove=0
onlyIfNeeded=0
computeRefIfNeeded=0
prepCodeOpts=""
perffile=""

while [ -n "$1" ]; do
  case "$1" in
    '-h') usage; exit;;
    '-s') suppress=1;;
    '-p') packcreation=1;;
    '-u') packupdate=1;;
    '-c') compilation=1;;
    '-r') run=$(($run+1));;
    '-C') check=1;;
    '-t') tests="$2"; shift;;
    '--expand') useexpand=1;;
    '--repo-user') export PHYEXREPOuser=$2; shift;;
    '--repo-protocol') export PHYEXREPOprotocol=$2; shift;;
    '--remove') remove=1;;
    '--onlyIfNeeded') onlyIfNeeded=1;;
    '--computeRefIfNeeded') computeRefIfNeeded=1;;
    '--prep_code-opts') prepCodeOpts=$2; shift;;
    '--perf') perffile="$(realpath $2)"; shift;;
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

if [ $packcreation -eq 0 -a \
     $packupdate -eq 0 -a \
     $compilation -eq 0 -a \
     $run -eq 0 -a \
     $check -eq 0 -a \
     $remove -eq 0 ]; then
  packcreation=1
  compilation=1
  run=1
  check=1
fi

if [ -z "${commit-}" ]; then
  echo "At least one commit hash must be provided on command line"
  exit 2
fi

##############################
#### FUNCTION DEFINITIONS ####
##############################

function json_dictkey2value {
  # $1 must contain the json string
  # $2 must be the key name
  # $3 is the default value
  json_content="$1" python3 -c "import json; import os; result=json.loads(os.environ['json_content']).get('$2', '$3'); print(json.dumps(result) if isinstance(result, dict) else result)"
}

function mvdiff {
  # $1 is the file to move
  # $2 is the destination
  if [ -d $2 ]; then
    #$2 is a directory and file must be moved in this directory, keeping its name
    dest=$2/$(basename $1)
  else
    dest=$2
  fi
  if [ $packupdate -eq 0 -o ! -f $dest ]; then
    #When creating the pack or if destination file doesn't exist
    mv $1 $2
  else
    #Destination file exists and we are in the update mode
    #File is moved only if different
    if ! cmp $1 $dest > /dev/null; then
      mv $1 $2
    else
      rm $1
    fi
  fi
}

###################################
#### VERSION/COMMIT ADAPTATION ####
###################################

#Name and directory for compiling and executing user pack
declare -A refByTest
if echo $commit | grep '/' | grep -v '^tags/' > /dev/null; then
  fromdir=$commit
  content_mesonh_version=$(scp $commit/src/mesonh/mesonh_version.json /dev/stdout 2>/dev/null || echo "")
else
  fromdir=''
  if [[ $commit == mesonh${separator}* ]]; then
    mesonh_version_file="mesonh_version.json"
  else
    mesonh_version_file="src/mesonh/mesonh_version.json"
  fi
  if echo $commit | grep '^tags/' > /dev/null; then
    urlcommit=$(echo $commit | cut -d / -f 2-)
  else
    urlcommit=$commit
  fi
  content_mesonh_version=$(wget --no-check-certificate https://raw.githubusercontent.com/$PHYEXREPOuser/PHYEX/${urlcommit}/$mesonh_version_file -O - 2>/dev/null || echo "")
fi
if [ ! "${content_mesonh_version}" == "" ]; then
  testing=$(json_dictkey2value "$content_mesonh_version" 'testing' '')
  if [ ! "$testing" == "" ]; then
    ALLTests='' #We reset the list of tests
    for t in $(echo $allowedTests | sed 's/,/ /g'); do
      ref=$(json_dictkey2value "$testing" "$t" '')
      if [ ! "$ref" == "" ]; then
        ALLTests="${ALLTests},$t"
        refByTest[$t]=$ref
      fi
    done
    ALLTests="${ALLTests:1}" #Remove first character (',')
  fi
fi
[ "${content_mesonh_version}" == "" ] && content_mesonh_version='{}'
refversion=$(json_dictkey2value "$content_mesonh_version" 'refversion' 'MNH-V5-5-0')
if [ $refversion == "MNH-V5-5-0" ]; then
  targzsuffix="_PHYEX"
else
  targzsuffix=""
fi
tag=$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
name=${refversion}-$tag
[ $suppress -eq 1 -a -d $MNHPACK/$name ] && rm -rf $MNHPACK/$name

#Name and directory for the reference version
declare -A refnameByTest
#Reference to use for each test
for t in $(echo $ALLTests | sed 's/,/ /g'); do
  #Name of the reference
  if [ "$reference" == "" -o "$reference" == "REF" ]; then
    if [[ ! -z "${refByTest[$t]+unset}" ]]; then #the -v test is valid only with bash > 4.3
      #The json file contained the references to use on a per test case basis
      reftag=${refByTest[$t]}
    else
      reftag=""
    fi
    refByTest[$t]=$reftag
  else
    if echo $reference | grep '/' > /dev/null; then
      reftag=$(echo $reference | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
    else
      reftag=$reference
    fi
    refByTest[$t]=$reference
  fi
  #Conversion into directory name
  if [ "$reftag" == "" ]; then
    refname=${refversion}
  else
    refname=${refversion}-$reftag
  fi
  refnameByTest[$t]=$refname
done

if [ -z "${tests-}" ]; then
  tests=$defaultTest
elif echo "$tests" | grep -w 'ALL' > /dev/null; then
  tests=$(echo "$tests" | sed "s:\bALL\b:$ALLTests:g")
fi

#######################
#### PACK CREATION ####
#######################

if [ $packcreation -eq 1 -a -d $MNHPACK/$name -a $onlyIfNeeded -eq 1 ]; then
  packcreation=0
fi
if [ $packcreation -eq 1 ]; then
  if [ -d $MNHPACK/$name ]; then
    echo "Pack already exists ($MNHPACK/$name), suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
    exit 5
  else
    echo "### Pack creation for commit $commit"

    # Prepare the pack
    cd $MNHPACK
    mkdir ${name}_$$
    cd ${name}_$$
    cp $TARGZDIR/${refversion}${targzsuffix}.tar.gz .
    tar xfz ${refversion}${targzsuffix}.tar.gz 
    rm ${refversion}${targzsuffix}.tar.gz
    mv ${refversion} ../$name
    cd ..
    rmdir ${name}_$$
    cd $name/src

    # Routine that changed names
    [ -f PHYEX/turb/modd_diag_in_run.f90 ] && mv -f PHYEX/turb/modd_diag_in_run.f90 MNH/. #To be removed once, this is done in MNH-git-lfs repo before inclusion of last version of PHYEX
  fi
fi
if [ $packupdate -eq 1 -o $packcreation -eq 1 ]; then
  if [ !  -d $MNHPACK/$name/src/PHYEX ]; then
      echo "PHYEX directory doesn't exist in pack ($MNHPACK/$name)"
      exit 9
  else
    cd $MNHPACK/$name/src
    if [ $packupdate -eq 1 ]; then
      mv PHYEX PHYEXori
    else
      rm -rf PHYEX
    fi

    if [ $useexpand == 1 ]; then
      expand_options="--mnhExpand"
    else
      expand_options=""
    fi
    subs="-s turb -s micro -s aux -s ext -s conv"
    prep_code=$PHYEXTOOLSDIR/prep_code.sh
    if [ "$fromdir" == '' ]; then
      echo "Clone repository, and checkout commit $commit (using prep_code.sh)"
      if [[ $commit == mesonh${separator}* ]]; then
        $prep_code $prepCodeOpts --renameFf --ilooprm -c $commit PHYEX #This commit is ready for inclusion
      else
        $prep_code $prepCodeOpts --renameFf --ilooprm -c $commit $expand_options $subs -m mesonh PHYEX
      fi
    else
      echo "Copy $fromdir"
      mkdir PHYEX
      scp -q -r $fromdir/src PHYEX/
      $prep_code $prepCodeOpts --renameFf --ilooprm $expand_options $subs -m mesonh PHYEX
    fi
    rm -rf PHYEX/.git
    find PHYEX -type f -exec touch {} \; #to be sure a recompilation occurs

    # Move manually ext/ files in src/MNH
    [ -f PHYEX/ext/yomhook.f90 ] && mv PHYEX/ext/yomhook.f90 PHYEX/ext/yomhook.F90
    if [ -d PHYEX/ext ]; then
      for file in PHYEX/ext/*; do
        mvdiff $file MNH/
      done
      if [ $packupdate -eq 1 ]; then
        #Only modified files are moved
        rm -rf PHYEX/ext
      else
        #All files must have been moved
        rmdir PHYEX/ext
      fi
    fi

    if [ $packupdate -eq 1 ]; then
      #Update only modified files
      cd PHYEX
      for file in $(find turb micro conv aux -type f); do
        if ! cmp $file ../PHYEXori/$file > /dev/null; then
          mv $file ../PHYEXori/$file
        fi
      done
      cd ..
      rm -rf PHYEX
      mv PHYEXori PHYEX
    fi

    if [ $packupdate -eq 0 ]; then
      cd $MNHPACK/$name/src/PHYEX/turb
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
      # Supress some files if they are not used anymore
      ! grep -i MODI_COMPUTE_ENTR_DETR $(ls MNH/*compute_updraft* PHYEX/turb/*compute_updraft* 2>/dev/null) && rm -f MNH/compute_entr_detr.f90
      ! grep -i MODI_TH_R_FROM_THL_RT_ $(ls MNH/compute_entr_detr.f90 MNH/compute_entr_detr.f90 PHYEX/turb/mode_compute_updraft*.f90 MNH/ice_adjust_bis.f90 MNH/prep_ideal_case.f90 MNH/set_rsou.f90 2>/dev/null)  > /dev/null && rm -f MNH/th_r_from_thl_rt_1d.f90 MNH/th_r_from_thl_rt_2d.f90 MNH/th_r_from_thl_rt_3d.f90
 
      # Routine that changed names
      #To be removed once, this is done in MNH-git-lfs repo before inclusion of last version of PHYEX
      rm -f PHYEX/micro/ini_rain_ice.f90
      rm -f PHYEX/micro/lima_nucleation_procs.f90
    fi

    # Remove binaries
    rm -f $MNHPACK/$name/exe/*

    # Remove execution results
    for t in $(echo $ALLTests | sed 's/,/ /g'); do
      case1=$(echo $t | cut -d / -f 1)
      case2=$(echo $t | cut -d / -f 2)
      casedir=$MNHPACK/$name/MY_RUN/KTEST/$case1
      [ -d $casedir/$case2 ] && rm -rf $casedir/$case2
    done
  fi
fi

#####################
#### COMPILATION ####
#####################

profile_sourced=0
if [ $compilation -eq 1 ]; then
  if [ $onlyIfNeeded -eq 0 -o ! -f $MNHPACK/$name/exe/MESONH* ]; then
    echo "### Compilation of commit $commit"
    cd $MNHPACK/$name/src
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

if [ $run -ge 1 ]; then
  #Cleaning to suppress old results that may be confusing in case of a crash during the run
  if [ $onlyIfNeeded -eq 0 ]; then
    for t in $(echo $tests | sed 's/,/ /g'); do
      case1=$(echo $t | cut -d / -f 1)
      case2=$(echo $t | cut -d / -f 2)
      casedir=$MNHPACK/$name/MY_RUN/KTEST/$case1
      [ -d $casedir/$case2 ] && rm -rf $casedir/$case2
    done
  fi

  #Run the tests one after the other
  firstrun=1
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then #test is allowed on this plateform
      case1=$(echo $t | cut -d / -f 1)
      case2=$(echo $t | cut -d / -f 2)
      casedir=$MNHPACK/$name/MY_RUN/KTEST/$case1
      if [ ! -d $casedir/$case2 ]; then #We do not enter systematically this part if onlyIfNeeded=1
        if [ $firstrun -eq 1 ]; then
          echo "### Running of commit $commit"
          firstrun=0
        fi

        if [ ! -f $MNHPACK/$name/exe/MESONH* ]; then
          echo "Pack does not exist ($MNHPACK/$name) or compilation has failed, please check"
          exit 6
        fi

        #If the test case didn't exist in the tar.gz, we copy it from from the reference version
        #and we suppress all the test directories for this case
        if [ ! -d $casedir ]; then
          cp -r $MNHPACK/${refversion}/MY_RUN/KTEST/$case1 $casedir/
          for newt in $(echo $ALLTests | sed 's/,/ /g'); do
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
            if ! echo $ALLTests | grep ${case1}/$d > /dev/null; then
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
          . $MNHPACK/$name/conf/profile_mesonh-*
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

if [ $check -eq 1 ]; then
  echo "### Check commit $commit against commit $reference"

  allt=0
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then
      #Run the reference if needed
      if [ $computeRefIfNeeded -eq 1 ]; then
        if [ "${refByTest[$t]}" == "" ]; then
          echo "Don't know how to compile/run $refversion"
          exit 3
        else
          #We must call it in another shell because of the potentially loaded MesoNH profile
          #because we cannot load two MesoNH profiles in the same shell
          env -i $SHELL -l -c "MNHPACK=${MNHPACK} TARGZDIR=${TARGZDIR} \
                               PHYEXREPOuser=${PHYEXREPOuser} PHYEXREPOprotocol=${PHYEXREPOprotocol} \
                               $0 -p -c -r -t $t --onlyIfNeeded ${refByTest[$t]}"
        fi
      fi

      #Files to compare
      refname=${refnameByTest[$t]}
      case1=$(echo $t | cut -d / -f 1)
      case2=$(echo $t | cut -d / -f 2)
      path_user=$MNHPACK/$name/MY_RUN/KTEST/$case1/$case2
      path_ref=$MNHPACK/$refname/MY_RUN/KTEST/$case1/$case2
      file3=""
      file4=""
      if [ $case1 == 007_16janvier ]; then
        file1=$path_user/16JAN.1.12B18.001.nc 
        file2=$path_ref/16JAN.1.12B18.001.nc
        file3=$path_user/16JAN.1.12B18.000.nc 
        file4=$path_ref/16JAN.1.12B18.000.nc
        bit_diff=57100
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
        bit_diff=76300
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

if [ $remove -eq 1 ]; then
  echo "### Remove model directory for commit $commit"
  [ -d $MNHPACK/$name ] && rm -rf $MNHPACK/$name
fi

exit $cmpstatus
