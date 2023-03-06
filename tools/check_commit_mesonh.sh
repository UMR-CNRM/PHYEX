#!/bin/bash

#set -x
set -e

#The folowing environment variables can be defined:
# REFDIR: directory in which the reference compilation directory can be found
# TARGZDIR: directory where tar.gz files are searched for
# MNHPACK: directory where tests are build

availTests="007_16janvier/008_run2, 007_16janvier/008_run2_turb3D, 007_16janvier/008_run2_lredf, COLD_BUBBLE/002_mesonh, 
            ARMLES/RUN, COLD_BUBBLE_3D/002_mesonh,OCEAN_LES/004_run2"
defaultTest="007_16janvier/008_run2"
separator='_' #- be carrefull, gmkpack (at least on belenos) has multiple allergies (':', '.', '@')
              #- seprator must be in sync with prep_code.sh separator

#Notes for v5.5.0
#For the OCEAN_LES/004_run2 case, results obtained are different from those obtained with the original version
#of Meso-NH because of new developments and bug correction. The reference version is given by commit e053c59.
#In this commit two modifications must be done in turb/mode_tke_eps_sources.f90 to change twice LOCEAN into OOCEAN.

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
function usage {
  echo "Usage: $0 [-h] [-c] [-r] [-C] [-s] [--expand]  [-t test] commit reference"
  echo "commit          commit hash (or a directory)"
  echo "reference       commit hash or a directory or nothing for ref"
  echo "-s              suppress compilation pack"
  echo "-c              performs compilation"
  echo "-r              runs the tests"
  echo "-C              checks the result against the reference"
  echo "-t              comma separated list of tests to execute"
  echo "                or ALL to execute all tests"
  echo "--expand        use mnh_expand (code will use do loops)"
  echo "--repo-user     user hosting the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOuser (=$PHYEXREOuser)"
  echo "--repo-protocol protocol (https or ssh) to reach the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOprotocol (=$PHYEXREOprotocol)"
  echo ""
  echo "If nothing is asked (compilation, running, check) everything is done"
  echo 
  echo "If no test is aked for, the default on ($defaultTest) is executed"
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
useexpand=0

while [ -n "$1" ]; do
  case "$1" in
    '-h') usage;;
    '-s') suppress=1;;
    '-c') compilation=1;;
    '-r') run=$(($run+1));;
    '-C') check=1;;
    '-t') tests="$2"; shift;;
    '--expand') useexpand=1;;
    '--repo-user') export PHYEXREPOuser=$2; shift;;
    '--repo-protocol') export PHYEXREPOprotocol=$2; shift;;
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

[ "$reference" == 'REF' ] && reference="" #Compatibility with check_commit_arome.sh

MNHPACK=${MNHPACK:=$HOME/MesoNH/PHYEX}
REFDIR=${REFDIR:=$PHYEXTOOLSDIR/pack/}
TARGZDIR=${TARGZDIR:=$PHYEXTOOLSDIR/pack/}
if [ -z "${tests-}" ]; then
  tests=$defaultTest
elif [ $tests == 'ALL' ]; then
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

#Name, directory and reference for compiling and executing user pack
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
refversion=$(content_mesonh_version=$content_mesonh_version python3 -c "import json, os; v=os.environ['content_mesonh_version']; print(json.loads(v if len(v)!=0 else '{}').get('refversion', 'MNH-V5-5-0'))")
if [ $refversion == "MNH-V5-5-0" ]; then
  targzsuffix="_PHYEX"
else
  targzsuffix=""
fi
tag=$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
name=${refversion}-$tag
[ $suppress -eq 1 -a -d $MNHPACK/$name ] && rm -rf $MNHPACK/$name

#Two possibilities are supported for the simulations
# - they can be done in the the pack we are currently checking
# - they can be done in the reference pack
#They are done in the current pack except if the reference pack
#already contains a tested simulation
#To check this, we use the case 007_16janvier/008_run2_turb3D
run_in_ref=$(ls -d $REFDIR/${refversion}/MY_RUN/KTEST/007_16janvier/008_run2_turb3D_* 2> /dev/null | tail -1 |wc -l)
if [ $run_in_ref -eq 1 ]; then
  path_user_beg=$REFDIR/${refversion} #pack directory containing the simulation
  path_user_end=_$tag #to be appended to the 'run' simulation directory
else
  path_user_beg=$MNHPACK/$name #pack directory containing the simulation
  path_user_end= #to be appended to the 'run' simulation directory
fi

#Name and directory for the reference
reffromdir=''
if echo $reference | grep '/' > /dev/null; then
  reffromdir=$reference
  reftag=$(echo $reference | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
else
  reftag=$reference
fi
refname=${refversion}-$reftag
if [ $run_in_ref -eq 1 ]; then
  path_ref_beg=$REFDIR/${refversion}
  if [ "$reference" == "" ]; then
    path_ref_end=
  else
    path_ref_end=_$reftag
  fi
else
  path_ref_end=
  if [ "$reference" == "" ]; then
    path_ref_beg=$REFDIR/${refversion}
  else
    path_ref_beg=$MNHPACK/${refversion}-$reftag
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
  cp $TARGZDIR/${refversion}${targzsuffix}.tar.gz .
  tar xfz ${refversion}${targzsuffix}.tar.gz 
  rm ${refversion}${targzsuffix}.tar.gz
  mv ${refversion} $name
  cd $name/src
  rm -rf PHYEX

  MNH_EXPAND_DIR=$PHYEXTOOLSDIR/mnh_expand
  export PATH=$MNH_EXPAND_DIR/filepp:$MNH_EXPAND_DIR/MNH_Expand_Array:$PATH

  if [ $useexpand == 1 ]; then
    expand_options="-D MNH_EXPAND -D MNH_EXPAND_LOOP"
  else
    expand_options=""
  fi
  subs="-s turb -s micro -s aux -s ext -s conv"
  prep_code=$PHYEXTOOLSDIR/prep_code.sh
  if [ "$fromdir" == '' ]; then
    echo "Clone repository, and checkout commit $commit (using prep_code.sh)"
    if [[ $commit == mesonh${separator}* ]]; then
      $prep_code --renameFf -c $commit PHYEX #This commit is ready for inclusion
    else
      $prep_code --renameFf -c $commit $expand_options $subs -m mesonh PHYEX
    fi
  else
    echo "Copy $fromdir"
    mkdir PHYEX
    scp -q -r $fromdir/src PHYEX/
    $prep_code --renameFf $expand_options $subs -m mesonh PHYEX
  fi
  rm -rf PHYEX/.git
  find PHYEX -type f -exec touch {} \; #to be sure a recompilation occurs

  # Move manually ext/ files in src/MNH
  if [ -d PHYEX/ext ]; then
    mv -f PHYEX/ext/* MNH/
    rmdir PHYEX/ext
  fi

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
  # nettoyage, routines non appellees : 
  rm -f MNH/mf_turb_greyzone.f90
  rm -f MNH/compute_frac_ice.f90
  rm -f MNH/rain_ice_red.f90
  # Supress some files if they are not used anymore
  ! grep -i MODI_COMPUTE_ENTR_DETR $(ls MNH/*compute_updraft* PHYEX/turb/*compute_updraft* 2>/dev/null) && rm -f MNH/compute_entr_detr.f90
  ! grep -i MODI_TH_R_FROM_THL_RT_ $(ls MNH/compute_entr_detr.f90 MNH/compute_entr_detr.f90 PHYEX/turb/mode_compute_updraft*.f90 MNH/ice_adjust_bis.f90 MNH/prep_ideal_case.f90 MNH/set_rsou.f90 2>/dev/null)  > /dev/null && rm -f MNH/th_r_from_thl_rt_1d.f90 MNH/th_r_from_thl_rt_2d.f90 MNH/th_r_from_thl_rt_3d.f90
  # Routine that changed names (if mode_budget.f90 is present)
  set +e
  mv -f PHYEX/aux/mode_budget.f90 MNH/budget.f90
  set -e

  #Configure and compilation
  command -v module && modulelist=$(module -t list 2>&1 | tail -n +2) #save loaded modules
  ./configure
  set +e #file ends with a test that can return false
  . ../conf/profile_mesonh-* #This lines modifies the list of loaded modules
  set -e
  make -j 8 | tee ../Output_compilation
  make installmaster | tee -a ../Output_compilation
  command -v module && module load $modulelist #restore loaded modules
fi

if [ $run -ge 1 ]; then
  echo "### Running of commit $commit"

  if [ ! -f $MNHPACK/$name/exe/MESONH* ]; then
    echo "Pack does not exist ($MNHPACK/$name) or compilation has failed, please check"
    exit 6
  fi

  for t in $(echo $tests | sed 's/,/ /g'); do
    case=$(echo $t | cut -d / -f 1)
    exedir=$(echo $t | cut -d / -f 2)
    if [ $run_in_ref -eq 1 ]; then
      cd $REFDIR/${refversion}/MY_RUN/KTEST/$case/
      [ ! -d ${exedir}_$commit ] && cp -R ${exedir} ${exedir}_$commit
      cd $REFDIR/${refversion}/MY_RUN/KTEST/$case/${exedir}_$commit
    else
      #If the test case didn't exist in the tar.gz, we copy it from from the reference version
      rep=$MNHPACK/$name/MY_RUN/KTEST/$case
      [ ! -d $rep ] && cp -r $REFDIR/${refversion}/MY_RUN/KTEST/$case $rep
      cd $rep

      #Loop on the directories
      for rep in *; do
        if [ -d "$rep" ]; then
          if echo $availTests | grep ${case}/$rep > /dev/null; then
            #This directory is a test case
            if [ $rep == ${exedir} ]; then
              #this is the case we want to run
              rm -rf $rep
              cp -r $REFDIR/${refversion}/MY_RUN/KTEST/$case/$rep .
            fi
          else
            #This directory might be neede to run the test case, we take the reference version
            rm -rf $rep
            ln -s $REFDIR/${refversion}/MY_RUN/KTEST/$case/$rep 
          fi
        fi
      done

      #In case subcase does not exist we create it
      [ ! -d ${exedir} ] && cp -r $REFDIR/${refversion}/MY_RUN/KTEST/$case/${exedir} .
      cd ${exedir}
    fi

    set +e #file ends with a test that can return false
    [ $compilation -eq 0 ] && . $MNHPACK/$name/conf/profile_mesonh-*
    set -e
    ./clean_mesonh_xyz
    ./run_mesonh_xyz | tee Output_run
  done
fi

if [ $check -eq 1 ]; then
  echo "### Check commit $commit against commit $reference"

  allt=0
  for t in $(echo $tests | sed 's/,/ /g'); do
    case=$(echo $t | cut -d / -f 1)
    exedir=$(echo $t | cut -d / -f 2)
    if [ $t == 007_16janvier/008_run2 ]; then
      path_user=$path_user_beg/MY_RUN/KTEST/007_16janvier/008_run2$path_user_end
      path_ref=$path_ref_beg/MY_RUN/KTEST/007_16janvier/008_run2$path_ref_end
    elif  [ $t == 007_16janvier/008_run2_turb3D ]; then
      path_user=$path_user_beg/MY_RUN/KTEST/007_16janvier/008_run2_turb3D$path_user_end
      path_ref=$path_ref_beg/MY_RUN/KTEST/007_16janvier/008_run2_turb3D$path_ref_end
    elif  [ $t == 007_16janvier/008_run2_lredf ]; then
      path_user=$path_user_beg/MY_RUN/KTEST/007_16janvier/008_run2_lredf$path_user_end
      path_ref=$path_ref_beg/MY_RUN/KTEST/007_16janvier/008_run2_lredf$path_ref_end
    elif   [ $t == COLD_BUBBLE/002_mesonh ]; then
      path_user=$path_user_beg/MY_RUN/KTEST/COLD_BUBBLE/002_mesonh$path_user_end
      path_ref=$path_ref_beg/MY_RUN/KTEST/COLD_BUBBLE/002_mesonh$path_ref_end
    elif   [ $t == COLD_BUBBLE_3D/002_mesonh ]; then
      path_user=$path_user_beg/MY_RUN/KTEST/COLD_BUBBLE_3D/002_mesonh$path_user_end
      path_ref=$path_ref_beg/MY_RUN/KTEST/COLD_BUBBLE_3D/002_mesonh$path_ref_end
    elif   [ $t == ARMLES/RUN ]; then
      path_user=$path_user_beg/MY_RUN/KTEST/ARMLES/RUN$path_user_end
      path_ref=$path_ref_beg/MY_RUN/KTEST/ARMLES/RUN$path_ref_end
    elif   [ $t == OCEAN_LES/004_run2 ]; then
      path_user=$path_user_beg/MY_RUN/KTEST/OCEAN_LES/004_run2$path_user_end
      path_ref=$path_ref_beg/MY_RUN/KTEST/OCEAN_LES/004_run2$path_ref_end
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

    if [ $case == 007_16janvier ]; then
      echo "Compare with python..."
      # Compare variable of both Synchronous and Diachronic files with printing difference
      file1=$path_user/16JAN.1.12B18.001.nc 
      file2=$path_ref/16JAN.1.12B18.001.nc
      file3=$path_user/16JAN.1.12B18.000.nc 
      file4=$path_ref/16JAN.1.12B18.000.nc
      set +e
      $PHYEXTOOLSDIR/compare.py --f1 $file1 --f2 $file2 --f3 $file3 --f4 $file4
      t=$?
      set -e
      allt=$(($allt+$t))
      
      #Check bit-repro before date of creation of Synchronous file from ncdump of all values (pb with direct .nc file checks)
      echo "Compare with ncdump..."
      if [ -f $file1 -a -f $file2 ]; then
        set +e
        diff <(ncdump $file1 | head -c 62889) <(ncdump $file2 | head -c 62889)
        t=$?
        set -e
        allt=$(($allt+$t))
      else
        [ ! -f $file1 ] && echo "  $file1 is missing"
        [ ! -f $file2 ] && echo "  $file2 is missing"
        allt=$(($allt+1))
      fi
    fi

    if [ $case == COLD_BUBBLE ]; then
      echo "Compare with python..."
      # Compare variable of both Synchronous files with printing difference
      file1=$path_user/BUBBL.1.CEN4T.001.nc
      file2=$path_ref/BUBBL.1.CEN4T.001.nc
      set +e
      $PHYEXTOOLSDIR/compare.py --f1 $file1 --f2 $file2
      t=$?
      set -e
      allt=$(($allt+$t))
      
      #Check bit-repro before date of creation of Synchronous file from ncdump of all values (pb with direct .nc file checks)
      echo "Compare with ncdump..."
      if [ -f $file1 -a -f $file2 ]; then
        set +e
        diff <(ncdump $file1 | head -c 27300) <(ncdump $file2 | head -c 27300)
        t=$?
        set -e
        allt=$(($allt+$t))
      else
        [ ! -f $file1 ] && echo "  $file1 is missing"
        [ ! -f $file2 ] && echo "  $file2 is missing"
        allt=$(($allt+1))
      fi
    fi

   if [ $case == OCEAN_LES ]; then
        echo "Compare with python..."
        # Compare variable of both Synchronous files with printing difference
        file1=$path_user/SPWAN.2.25m00.001.nc
        file2=$path_ref/SPWAN.2.25m00.001.nc
        set +e
        $PHYEXTOOLSDIR/compare.py --f1 $file1 --f2 $file2
        t=$?
        set -e
        allt=$(($allt+$t))
  
        #Check bit-repro before date of creation of Synchronous file from ncdump of all values (pb with direct .nc file checks)
        echo "Compare with ncdump..."
        if [ -f $file1 -a -f $file2 ]; then
          set +e
          diff <(ncdump $file1 | head -c 18400) <(ncdump $file2 | head -c 18400)
          t=$?
          set -e
          allt=$(($allt+$t))
        else
          [ ! -f $file1 ] && echo "  $file1 is missing"
          [ ! -f $file2 ] && echo "  $file2 is missing"
          allt=$(($allt+1))
        fi
      fi

    if [ $case == COLD_BUBBLE_3D ]; then
      echo "Compare with python..."
      # Compare variable of both Synchronous and Diachronic files with printing difference
      file1=$path_user/BUBBL.1.CEN4T.001.nc
      file2=$path_ref/BUBBL.1.CEN4T.001.nc
      file3=$path_user/BUBBL.1.CEN4T.000.nc
      file4=$path_ref/BUBBL.1.CEN4T.000.nc
      set +e
      $PHYEXTOOLSDIR/compare.py --f1 $file1 --f2 $file2 --f3 $file3 --f4 $file4
      t=$?
      set -e
      allt=$(($allt+$t))

      #Check bit-repro before date of creation of Synchronous file from ncdump of all values (pb with direct .nc file checks)
      echo "Compare with ncdump..."
      if [ -f $file1 -a -f $file2 ]; then
        set +e
        diff <(ncdump $file1 | head -c 27300) <(ncdump $file2 | head -c 27300)
        t=$?
        set -e
        allt=$(($allt+$t))
      else
        [ ! -f $file1 ] && echo "  $file1 is missing"
        [ ! -f $file2 ] && echo "  $file2 is missing"
        allt=$(($allt+1))
      fi
    fi

    if [ $case == ARMLES ]; then
      echo "Compare with python..."
      # Compare variable of both Synchronous and Diachronic files with printing difference
      file1=$path_user/ARM__.1.CEN4T.001.nc
      file2=$path_ref/ARM__.1.CEN4T.001.nc
      file3=$path_user/ARM__.1.CEN4T.000.nc
      file4=$path_ref/ARM__.1.CEN4T.000.nc
      set +e
      $PHYEXTOOLSDIR/compare.py --f1 $file1 --f2 $file2 --f3 $file3 --f4 $file4
      t=$?
      set -e
      allt=$(($allt+$t))

      #Check bit-repro before date of creation of Synchronous file from ncdump of all values (pb with direct .nc file checks)
      echo "Compare with ncdump..."
      if [ -f $file1 -a -f $file2 ]; then
        set +e
        diff <(ncdump $file1 | head -c 62889) <(ncdump $file2 | head -c 62889)
        t=$?
        set -e
        allt=$(($allt+$t))
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
  fi
  echo "...comparison done: $status"
fi
