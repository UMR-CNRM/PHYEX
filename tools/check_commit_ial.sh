#!/bin/bash

#set -x
set -e
set -o pipefail #abort if left command on a pipe fails

#This script:
# - compiles the AROME model using a specific commit for the externalised physics
# - runs a small 3D case and checks if results are identical to a given version

#small_3D_np2: on only 2 procs
#small_3D_alt1: options around time-step dependency, CFRAC_ICE_*='T', CSEDIM='SPLI', LSEDIM_AFTER=.T.
#small_3D_alt2: CCLOUD='OLD3'
#small_3D_alt3: PRFR
#small_3D_alt4: small_3D_alt1 + CSNOWRIMING='OLD'
#small_3D_alt5: CCLOUD='ICE4'
#small_3D_alt6: CMF_UPDRAFT='RAHA', CMF_CLOUD='BIGA'
#small_3D_alt7: CMF_CLOUD='STAT', LOSIGMAS=.FALSE. #Needs 2 corrections in original cycle 48
#small_3D_alt8: CMF_UPDRAFT='RHCJ'
#small_3D_alt9: CCLOUD='OLD3', OCND2=.T.
#small_3D_alt10: LCRIAUTI=F
#small_3D_lima: LIMA scheme
#small_3D_alt11: same as small_3D but with a different value for NPROMICRO (must give exactly the same results)
#small_3D_alt12: same as small_3D but with LPACK_MICRO=.F. (must give exactly the same results)
#small_3D_xfrmin: same as small_3D_alt2 but with specified values for XFRMIN(16:17)
#arp_t31: Ph. Marguinaud's ARPEGE toy

#When running in 49t0 after the f065e64 commit (23 June 2023) all configurations must be compared to this same commit.
#79fe47e (previous commit) is identical to the different references for all the test cases.
#When running in 49t0 after the 00148b1 commit (27 June 2023) all configurations must be compared to this same commit.

#The small_3D_alt7 needed a correction in apl_arome which has been introduced in d37dd1f. But the reference pack has been modified
#                  afterwards to enable this test case to be run (documented in INSTALL_pack_ial.md). In consequence, the reference
#                  to use is the same as for the other test cases and this case cannot be run for commit before d37dd1f (20 April 2022).

#The small_3D_alt8 is not included in the list of available tests because it needs to be compared against a special commit.
#                  Indeed, on 3 February 2022 (commit 907e906) the mesonh version of compute_updraft_rhcj.F90 has been put in the common directory.
#                  The reference is
#                       the commit 907e906 when running in 48t1
#                       the commit d10ed48 when running in 48t3 (edc3f88 (last commit in 48t1) is identical to 907e906)
#                       the commit 7e55649 when running in 49t0 (9164c67 (last commit in 48t3) is identical to d10ed48)

#The small_3D_alt9 is not included in the list of available tests because it needs to be compared against a special commit.
#                  Indeed, some pieces are missing in the reference pack.
#                  Theses pieces have been added in commit edc3f88 during phasing with 48t3.
#                  The reference is
#                       the commit edc3f88 (21 September 2022) when running in 48t1
#                       the commit d10ed48 in 48t3 (29 september 2022) when running in 48t3
#                       the commit 110a5aa in 49t0 (13 June 2023) when running in 49t0 (bd44ba7 (patch on the last commint in 48t3) is identical to d10ed48)

#The small_3D_alt10 is not included in the list because it is not sufficiently different from other tests
#                  Be careful that namelists were wrong before commit 3c01df4 (8 June 2023)

#The small_3D_lima is not included in the list of available tests because it needs to be compared against a special commit.
#                  Indeed, the lima version in arome has been changed.
#                  The reference commit is
#                       the commit d095d11 (20 March 2023) when running in 48t3
#                       the commit 7e55649 when running in 49t0 (9164c67 (last commit in 48t3) is identical to d095d11)

#######################
#### CONFIGURATION ####
#######################

#Special pack names:
# - recompil: original source code (everything under mpa)
# - split_48t1: original 48t1 source code but with physics source code under phyex directory
# - split_48t3: same as split_48t1 but for the 48t3 cycle
# - split: symbolic link to split_48t1 (backward compatibility)
# - split_49t0: same as split_48t1 but for the 49t0 cycle
specialPack="ori split split_48t1 split_48t3 recompil split_49t0"

#About the tests:
# - ALLTests is a list of tests to be done when '-t ALL' is used. This list is filled here
#   in case there is no ial_version.json file containig a 'testing' section. If this 'testing'
#   section exists, this list is overridden.
# - allowedTests is the list of allowed tests which can depend on platform, if we ask to perform an action
#   with a test not in the allowedTests list, the action is ignored
# - defaultTest is the list of tests to perform when no '-t' option is provided on the command line.
ALLTests="small_3D,small_3D_np2,small_3D_alt1,small_3D_alt2,small_3D_alt3,small_3D_alt4,small_3D_alt5,small_3D_alt6,small_3D_alt7"
defaultTest="small_3D"
allowedTests="small_3D,small_3D_np2,small_3D_alt1,small_3D_alt2,small_3D_alt3,small_3D_alt4,small_3D_alt5,small_3D_alt6,small_3D_alt7,small_3D_alt8,small_3D_alt9,small_3D_alt10,small_3D_alt11,small_3D_alt12,small_3D_lima,small_3D_xfrmin,arp_t31"

separator='_' #- be carrefull, gmkpack (at least on belenos) has multiple allergies (':', '.', '@')
              #- seprator must be in sync with prep_code.sh separator

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

HOMEPACK=${HOMEPACK:=$HOME/pack}

dirpack=$PHYEXTOOLSDIR/pack
dirconf=$PHYEXTOOLSDIR/conf_tests
declare -A gmkpack_l
declare -A gmkpack_o
if [ $(hostname | cut -c 1-7) == 'belenos' -o $(hostname | cut -c 1-7) == 'taranis' ]; then
  HPC=1
  gmkpack_l[default]=MIMPIIFC1805
  gmkpack_l[49t0]=IMPIIFC2018
  gmkpack_l[49t2]=IMPIIFCI2302DP
  gmkpack_l[50t1]=IMPIIFCI2302REPRODP
  gmkpack_o[default]=2y
  gmkpack_o[49t0]=x
  gmkpack_o[49t2]=y
  gmkpack_o[50t1]=y
  defaultMainPackVersion=01
  defaultRef='split_${cycle}'
  ALLTests="${ALLTests},big_3D"
  allowedTests="${allowedTests},big_3D"
else
  HPC=0
  gmkpack_l[default]=MPIGFORTRAN920DBL
  gmkpack_l[49t0]=OMPIGFORT920DBL
  gmkpack_l[49t2]=OMPIGFORT920DBL
  gmkpack_l[50t1]=OMPI5GFORT141DP50
  gmkpack_o[default]=xfftw
  gmkpack_o[49t0]=x
  gmkpack_o[49t2]=x
  gmkpack_o[50t1]=x
  defaultMainPackVersion=01
  defaultRef='split_${cycle}'
fi
mainPackVersion=${mainPackVersion:-${defaultMainPackVersion}}

################################
#### COMMAND LINE ARGUMENTS ####
################################

function usage {
  echo "Usage: $0 [-h] [-p] [-u] [-c] [-r] [-C] [-s] [-f] [--noexpand] [-t TEST] [--cycle CYCLE] [--scripttag TAG] [--repo-user USER] [--repo-protocol PROTOCOL] [--remove] [--onlyIfNeeded] [--computeRefIfNeeded] [--prep_code-opts 'OPTS'] [--perf FILE] commit [reference]"
  echo "commit          commit hash (or a directory, or among $specialPack) to test"
  echo "reference       commit hash (or a directory, or among $specialPack) REF to use as a reference"
  echo "-s              suppress compilation pack"
  echo "-p              creates pack"
  echo "-u              updates pack (experimental, only for 'small' updates where"
  echo "                              the list of files does not change)"
  echo "-c              performs compilation"
  echo "-r              runs the tests"
  echo "-C              checks the result against the reference"
  echo "-t TEST         comma separated list of tests to execute"
  echo "                or ALL to execute all tests"
  echo "--noexpand      do not expand mnh_expand blocks (code will be in array-syntax)"
  echo "-f              full compilation (do not use pre-compiled pack)"
  echo "--cycle CYCLE   to force using CYCLE"
  echo "--scripttag TAG script tag to use in case --cycle is provided"
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
  echo
  echo "The cycle will be guessed from the source code"
  echo
  echo "The -f flag (full recompilation) is active only at pack creation"
  echo
  echo "The PHYEXROOTPACK environment variable, if set, is used as the argument"
  echo "of the --rootpack option of ial-git2pack, for incremental packs."
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
useexpand=1
fullcompilation=0
cycle=""
scripttag=''
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
    '--noexpand') useexpand=0;;
    '-f') fullcompilation=1;;
    '--cycle') cycle="$2"; shift;;
    '--scripttag') scripttag="$2"; shift;;
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

if [ $check -eq 1 -a -z "${reference-}" ]; then
  echo "To perform a comparison two commit hashes are mandatory on the command line"
  exit 3
fi

##############################
#### FUNCTION DEFINITIONS ####
##############################

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

function apl_arome_content2cycle {
  # variable content_apl_arome must contain the source code of apl_arome.F90
  if grep CPG_DYN_TYPE <(echo $content_apl_arome) > /dev/null; then
    echo 48t3
  else
    echo 48t1
  fi
}

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

#################################
#### CYCLE/COMMIT ADAPTATION ####
#################################

is_directory=0
is_special=0
is_commit=0
if echo $commit | grep '/' | grep -v '^tags/' > /dev/null; then
  is_directory=1
elif echo $specialPack | grep -w $commit > /dev/null; then
  is_special=1
else
  is_commit=1
fi

#Name is choosen such as it can be produced with a main pack: PHYEX/${cycle}_XXXXXXXXX.01.${gmkpack_l}.${gmkpack_o}
declare -A refByTest
fromdir=''
if [ $is_directory -eq 1 -o $is_commit -eq 1 ]; then
  if [ $is_directory -eq 1 ]; then
    #The git repository is a directory
    fromdir=$commit
    if [ -d $commit/src ]; then
      content_ial_version=$(scp $commit/src/arome/ial_version.json /dev/stdout 2>/dev/null || echo "")
      cmd_apl_arome="scp $commit/src/arome/ext/apl_arome.F90 /dev/stdout"
    else
      content_ial_version=$(scp $commit/ial_version.json /dev/stdout 2>/dev/null || echo "")
      cmd_apl_arome="scp $commit/ext/apl_arome.F90 /dev/stdout"
    fi
    packBranch=$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
  else
    #The git repository is on github
    if [[ $commit == arome${separator}* ]]; then
      apl_arome_file="ext/apl_arome.F90"
      ial_version_file="ial_version.json"
    else
      apl_arome_file="src/arome/ext/apl_arome.F90"
      ial_version_file="src/arome/ial_version.json"
    fi
    if echo $commit | grep '^tags/' > /dev/null; then
      urlcommit=$(echo $commit | cut -d / -f 2-)
    else
      urlcommit=$commit
    fi
    content_ial_version=$(wget --no-check-certificate https://raw.githubusercontent.com/$PHYEXREPOuser/PHYEX/${urlcommit}/$ial_version_file -O - 2>/dev/null || echo "")
    cmd_apl_arome="wget --no-check-certificate https://raw.githubusercontent.com/$PHYEXREPOuser/PHYEX/${urlcommit}/$apl_arome_file -O - 2>/dev/null"
    packBranch="COMMIT$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')"
  fi
  if [ ! "$content_ial_version" == "" ]; then
    testing=$(json_dictkey2value "$content_ial_version" 'testing' '')
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
  if [ "$cycle" == "" ]; then
    if [ "$content_ial_version" == "" ]; then
      content_apl_arome=$($cmd_apl_arome)
      cycle=$(apl_arome_content2cycle)
    else
      cycle=$(json_dictkey2value "$content_ial_version" 'cycle' '')
      scripttag=$(json_dictkey2value "$content_ial_version" 'scripttag' '')
    fi
  fi
  if [[ ! -z "${gmkpack_o[$cycle]+unset}" ]]; then #the -v approach is valid only with bash > 4.3:  if [[ -v gmkpack_o[$cycle] ]]; then
    gmkpack_o=${gmkpack_o[$cycle]}
  else
    gmkpack_o=${gmkpack_o[default]}
  fi
  if [[ ! -z "${gmkpack_l[$cycle]+unset}" ]]; then
    gmkpack_l=${gmkpack_l[$cycle]}
  else
    gmkpack_l=${gmkpack_l[default]}
  fi
  name="PHYEX/${cycle}_${packBranch}.01.${gmkpack_l}.${gmkpack_o}"
  [ $suppress -eq 1 -a -d $HOMEPACK/$name ] && rm -rf $HOMEPACK/$name
elif [ $is_special -eq 1 ]; then
  name="PHYEX/$commit"
  if [ $commit == split_49t0 ]; then
    cycle=49t0
  elif [ $commit == split_48t3 ]; then
    cycle=48t3
  else
    cycle=48t1
  fi
fi
if [ ! -z "${reference-}" ]; then
  declare -A refnameByTest
  #Reference to use for each test
  for t in $(echo $allowedTests | sed 's/,/ /g'); do
    #Name of the reference
    if [ "$reference" == 'REF' ]; then
      if [[ ! -z "${refByTest[$t]+unset}" ]]; then #the -v approach is valid only with bash > 4.3:  if [[ -v gmkpack_o[$cycle] ]]; then
        #The json file contained the references to use on a per test case basis
        caseref=${refByTest[$t]}
      else
        #No json file, we use the global default reference
        caseref=$(eval echo $defaultRef) #echo to replace ${cycle} by value
      fi
    else
      #The exact reference to use was given on the command line
      caseref=$reference
    fi
    refByTest[$t]=$caseref
    #Conversion into directory name
    if echo $caseref | grep '/' > /dev/null; then
      refname="PHYEX/*_$(echo $caseref | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g').01.${gmkpack_l}.${gmkpack_o}"
    elif echo $specialPack | grep -w $caseref > /dev/null; then
      refname="PHYEX/$caseref"
    else
      refname="PHYEX/*_COMMIT${caseref}.01.${gmkpack_l}.${gmkpack_o}"
    fi
    refnameByTest[$t]=$refname
  done
fi

if [ -z "${tests-}" ]; then
  tests=$defaultTest
elif echo "$tests" | grep -w 'ALL' > /dev/null; then
  tests=$(echo "$tests" | sed "s/\bALL\b/$ALLTests/g")
fi

#######################
#### PACK CREATION ####
#######################

if [ $packcreation -eq 1 -a -d $HOMEPACK/$name -a $onlyIfNeeded -eq 1 ]; then
  packcreation=0
fi
if [ $packcreation -eq 1 ]; then
  if [ -d $HOMEPACK/$name ]; then
    echo "Pack already exists ($HOMEPACK/$name), suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
    exit 5
  else
    echo "### Pack creation for commit $commit"

    if echo $specialPack | grep -w $commit > /dev/null; then
      echo "Special commit '$commit' cannot be compiled with this script"
      exit 4
    fi

    export GMKTMP=/dev/shm

    if [ $cycle == '49t2' -o $cycle == '50t1' ]; then
      if ! which ial-git2pack > /dev/null 2>&1; then
        echo "ial-git2pack not found, please install it"
        exit 7
      fi
      #Pack creation using ial-git2pack
      tmpbuilddir=$(mktemp -d)
      trap "\rm -rf $tmpbuilddir" EXIT
      cd $tmpbuilddir
      git clone $(json_dictkey2value "$content_ial_version" 'IALrepo' 'git@github.com:ACCORD-NWP/IAL.git')
      cd IAL
      git checkout $(json_dictkey2value "$content_ial_version" 'IALcommit' 'XXXXXXX')
      IALbundle_tag=$(json_dictkey2value "$content_ial_version" 'IALbundle_tag' '')
      [ "$IALbundle_tag" != "" ] && IALbundle_tag="--hub_bundle_tag $IALbundle_tag"
      ROOTPACKopt=""
      if [ $fullcompilation == 0 ]; then
        kind=incr
        if [ "$PHYEXROOTPACK" != "" ]; then
          ROOTPACKopt="--rootpack $PHYEXROOTPACK"
        fi
      else
        kind=main
      fi
      echo y | ial-git2pack -l ${gmkpack_l} -o ${gmkpack_o} -t $kind -n 10 \
                            -r $tmpbuilddir/IAL -p masterodb \
                            $ROOTPACKopt \
                            --homepack $tmpbuilddir/pack $IALbundle_tag

      #Moving
      oldname=$(echo $tmpbuilddir/pack/*)
      mv $oldname $HOMEPACK/$name
      for file in $HOMEPACK/$name/ics_*; do
        sed -i "s*$oldname*$HOMEPACK/$name*g" $file
      done
      rm -rf $tmpbuilddir

      cd $HOMEPACK/$name
      if [ $cycle == '49t2' ]; then
        #Workarounds for 49t2 compilation
        rm -rf src/local/obstat src/local/oopsifs
        rm -rf hub/local/src/Atlas hub/local/src/OOPS hub/local/src/ecSDK/?ckit
        if [ $fullcompilation != 0 ]; then
          gmkfile=$HOMEPACK/$name/.gmkfile/${gmkpack_l}.*
          for key in REPRO48 WITH_ODB WITHOUT_ECFLOW; do
            if ! grep "^MACROS_FRT" $gmkfile | grep -e -D$key; then
              sed -i "/^MACROS_FRT/s/$/ -D$key/" $gmkfile
            fi
          done
          if ! grep "^GMK_CMAKE_ectrans" $gmkfile | grep -e -DENABLE_ETRANS=ON; then
            sed -i "/^GMK_CMAKE_ectrans/s/$/ -DENABLE_ETRANS=ON/" $gmkfile
          fi
        fi
      fi
      if [ $cycle == '50t1' ]; then
        #Workarounds for 50t1 compilation
        rm -rf src/local/oopsifs
        rm -rf hub/local/src/Atlas hub/local/src/OOPS
      fi

      #Prepare PHYEX inclusion
      rm -rf src/local/phyex
      mkdir src/local/phyex

    elif [ $fullcompilation == 0 ]; then
      #Incremental compilation old style
      basepack=${cycle}_main.01.${gmkpack_l}.${gmkpack_o}
      #[ $HPC -eq 0 -a ! -d $ROOTPACK/$basepack ] &&  getpack $basepack
      gmkpack -r ${cycle} -b phyex -v $mainPackVersion -l ${gmkpack_l} -o ${gmkpack_o} -p masterodb \
              -f $dirpack/ \
              -u $name
    else
      #Main pack creation old style
      if [ $(echo $cycle | cut -c 1-2) -ne 48 ]; then
        hub='-K'
      else
        hub=''
      fi
      #Create main pack
      gmkpack -a $hub -r ${cycle} -b ${packBranch} -n 01 -l ${gmkpack_l} -o ${gmkpack_o} -p masterodb -h $HOMEPACK/PHYEX
      #Populate hub
      if [ -d $HOMEPACK/$name/hub ]; then
        cd $HOMEPACK/$name/hub/local/src
        if [ $HPC -eq 1 ]; then
          ssh sxphynh.cnrm.meteo.fr "wget http://anonymous:mto@webdav.cnrm.meteo.fr/public/algo/khatib/src/hub49.tgz -O -" > hub49.tgz
        else
          wget http://anonymous:mto@webdav.cnrm.meteo.fr/public/algo/khatib/src/hub49.tgz
        fi
        tar xf hub49.tgz
        rm -f hub49.tgz
      fi
      #Populate
      cd $HOMEPACK/$name/src/local/
      if [ $HPC -eq 1 ]; then
        ssh sxphynh.cnrm.meteo.fr "wget http://anonymous:mto@webdav.cnrm.meteo.fr/public/algo/khatib/src/${cycle}_main.01.tgz -O -" > ${cycle}_main.01.tgz
      else
        wget http://anonymous:mto@webdav.cnrm.meteo.fr/public/algo/khatib/src/${cycle}_main.01.tgz
      fi
      tar xf ${cycle}_main.01.tgz
      rm -f ${cycle}_main.01.tgz
      #Cleaning and moving
      if [ "$cycle" == '48t3' -o "$cycle" == '49t0' ]; then
        #extracting budgets from micro
        mkdir mpa/budgets
        for file in mpa/micro/module/moddb_intbudget.F90 mpa/micro/externals/aro_suintbudget_omp.F90 \
                    mpa/micro/interface/aro_convbu.h mpa/micro/externals/aro_convbu.F90 \
                    mpa/micro/interface/aro_startbu.h mpa/micro/externals/aro_startbu.F90 \
                    mpa/micro/externals/aro_suintbudget.F90 mpa/micro/externals/aro_suintbudget_omp.F90 \
                    mpa/micro/interface/aroini_budget.h mpa/micro/externals/aroini_budget.F90; do
          [ -f $file ] && mv $file mpa/budgets/
        done 
        mkdir mpa/aux
        for file in mpa/micro/interface/aroini_frommpa.h mpa/micro/externals/aroini_frommpa.F90 \
                    mpa/micro/externals/modd_spp_type.F90 mpa/micro/externals/spp_mod_type.F90 \
                    mpa/micro/interface/aroini_cstmnh.h mpa/micro/externals/aroini_cstmnh.F90; do
          [ -f $file ] && mv $file mpa/aux/
        done
        [ -f mpa/micro/externals/add_bounds.F90 ] && rm -f mpa/micro/externals/add_bounds.F90
        [ -f mpa/micro/externals/aroini_wet_dep.F90 ] && mv mpa/micro/externals/aroini_wet_dep.F90 mpa/chem/externals/aroini_wet_dep.F90
        [ -f mpa/micro/interface/aroini_wet_dep.h ] && mv mpa/micro/interface/aroini_wet_dep.h mpa/chem/interface/aroini_wet_dep.h
      fi
      if [ "$cycle" == '49t0' ]; then
        rm -rf oopsifs
      fi
      #we keep everything from the official source code except internals and module subdirectories of mpa
      #and except some files of mpa/conv/module
      for file in modi_shallow_convection.F90 modi_shallow_convection_part1.F90 \
                  modi_shallow_convection_part2.F90 modi_shallow_convection_part2_select.F90; do
        if [ -f mpa/conv/module/$file ]; then
          [ ! -d mpa/conv/module_save ] && mkdir mpa/conv/module_save
          mv mpa/conv/module/$file mpa/conv/module_save/
        fi
      done
      for rep in turb micro conv; do
        mkdir -p phyex/$rep
        rm -rf mpa/$rep/internals mpa/$rep/module
      done
      [ -d mpa/conv/module_save ] && mv mpa/conv/module_save mpa/conv/module
      if [ -f /cnrm/algo/khatib/drhook.c_for_ubuntu.tar -a $(echo $cycle | cut -c 1-2) -eq 48 ]; then
        #If file exists it means that we are running on a CTI computer, so we are using ubuntu
        tar xf /cnrm/algo/khatib/drhook.c_for_ubuntu.tar
      fi
      #Special modification of the compilation configuration file and script
      sed -i 's/-ftree-vectorize//' $HOMEPACK/$name/.gmkfile/${gmkpack_l}.*
      sed -i "/^MACROS_FRT/s/$/ -DREPRO48/" $HOMEPACK/$name/.gmkfile/${gmkpack_l}.*
      #sed -i "s/PHYEX\/${cycle}_$$.01.${gmkpack_l}.${gmkpack_o}/$(echo $name | sed 's/\//\\\//')/" $HOMEPACK/$name/ics_masterodb #this line could be used if pack was renamed before compilation but it does not work on belenos

      cleanpack -f #Is it really useful?
      resetpack -f #Is it really useful?
    fi
  fi
fi
if [ $packupdate -eq 1 -o $packcreation -eq 1 ]; then
  if [ ! -d $HOMEPACK/$name/src/local/phyex ]; then
    echo "phyex directory doesn't exist in pack ($HOMEPACK/$name)"
    exit 7
  else
    cd $HOMEPACK/$name/src/local/phyex

    if [ $useexpand == 1 ]; then
      expand_options="--mnhExpand"
    else
      expand_options=""
    fi
    subs="-s gmkpack_ignored_files -s turb -s micro -s aux -s ext -s conv -s externals" #externals is the old name for aux/ext
    prep_code=$PHYEXTOOLSDIR/prep_code.sh
    if [ "$fromdir" == '' ]; then
      echo "Clone repository, and checkout commit $commit (using prep_code.sh)"
      if [[ $commit == arome${separator}* ]]; then
        $prep_code $prepCodeOpts -c $commit PHYEX #This commit is ready for inclusion
      else
        $prep_code $prepCodeOpts -c $commit $subs -m arome PHYEX -- --shumanFUNCtoCALL --removeACC $expand_options
      fi
    else
      echo "Copy $fromdir"
      mkdir PHYEX
      if [ -d $fromdir/src ]; then
        scp -q -r $fromdir/src PHYEX/
        $prep_code $prepCodeOpts $subs -m arome PHYEX -- --shumanFUNCtoCALL --removeACC $expand_options
      else
        scp -q -r $fromdir/* PHYEX/
        $prep_code $prepCodeOpts PHYEX #Ready for inclusion
      fi
    fi
    find PHYEX -type f -exec touch {} \; #to be sure a recompilation occurs
    if [ $packupdate -eq 1 ]; then
      #Update only modified files
      cd PHYEX
      for file in $(find turb micro conv aux -type f); do
        mvdiff $file ../$file
      done
      cd ..
      rm -rf PHYEX
    else
      #Move PHYEX source files
      for rep in turb micro conv aux; do
        [ -d PHYEX/$rep ] && mv PHYEX/$rep .
      done
    fi
    #modd_nsv.F90 has been moved and gmkpack is lost in case a different version exists in main/.../micro
    if [ -f ../../main/phyex/micro/modd_nsv.F90 -a -f aux/modd_nsv.F90 ]; then
      mvdiff aux/modd_nsv.F90 micro/
      if [ -f PHYEX/gmkpack_ignored_files ]; then
        grep -v micro/modd_nsv.F90 PHYEX/gmkpack_ignored_files > PHYEX/gmkpack_ignored_files_new
        mv PHYEX/gmkpack_ignored_files_new PHYEX/gmkpack_ignored_files
      fi
    fi
    #shuman_phy.F90 has been renamed and gmkpack is lost in case a different version exists in main
    if [ -f ../../main/phyex/aux/shuman_phy.F90 -a -f aux/mode_shuman_phy.F90 ]; then
      mv aux/mode_shuman_phy.F90 aux/shuman_phy.F90
    fi
    if [ -f PHYEX/gmkpack_ignored_files ]; then
      #gmkpack_ignored_files contains a list of file, present in the reference pack, that is not used anymore
      #and must be excluded from compilation (in case of a full comilation) or from re-compilation (in case of a non-full
      #compilation).
      if [ $fullcompilation == 0 ]; then
        #Content is added in the ics_masterodb script
        if [ $packupdate -eq 0 ]; then
          if [ $(cat PHYEX/gmkpack_ignored_files | wc -l) != 0 ]; then
            sed -i "/^end_of_ignored_files/i $(first=1; for line in $(cat PHYEX/gmkpack_ignored_files); do echo -n $(test $first -ne 1 && echo \\n)${line}; first=0; done)" $HOMEPACK/$name/ics_masterodb
          fi
        fi
      else
        #Files must be suppressed (non phyex files)
        for file in $(cat PHYEX/gmkpack_ignored_files); do
          [ -f $HOMEPACK/$name/src/local/$file ] && rm -f $HOMEPACK/$name/src/local/$file
        done
        if [ -d $HOMEPACK/$name/src/local/mpa/dummy ]; then
          [ ! "$(ls -A $HOMEPACK/$name/src/local/mpa/dummy)" ] && rmdir $HOMEPACK/$name/src/local/mpa/dummy
        fi
      fi
    fi

    EXT=PHYEX/ext
    [ ! -d $EXT ] && EXT=PHYEX/externals #old name for ext/aux
    if [ -d $EXT ]; then
      #Move manually files outside of mpa (a find on the whole repository would take too much a long time)
      for file in yomparar.F90 cpg_opts_type_mod.fypp field_variables_mod.fypp cpg_type_mod.fypp \
                  field_registry_mod.fypp mf_phys_next_state_type_mod.fypp yemlbc_model.F90 \
                  field_config.yaml field_definitions.F90 field_gfl_wrapper.F90 \
                  yomfa.F90 yom_ygfl.F90 type_gflflds.F90 mf_phys_base_state_type_mod.fypp; do
        [ -f $EXT/$file ] && mvdiff $EXT/$file ../arpifs/module/
      done
      for file in namparar.nam.h namlima.nam.h namgfl.nam.h; do
        [ -f $EXT/$file ] && mvdiff $EXT/$file ../arpifs/namelist/
      done
      for file in aplpar.F90 acvppkf.F90 writemusc.F90 vdfhghtnhl.F90 suphmse.F90 suphmpa.F90 \
                  suparar.F90 apl_arome.F90 apl_arome_adjust.F90 apl_arome_micro.F90 \
                  apl_arome_shallow.F90 apl_arome_turbulence.F90; do
        [ -f $EXT/$file ] && mvdiff $EXT/$file ../arpifs/phys_dmn/
      done
      for file in cpg_pt_ulp_expl.fypp cpg_gp.F90; do
        [ -f $EXT/$file ] && mvdiff $EXT/$file ../arpifs/adiab/
      done
      for file in su0yomb.F90 suctrl_gflattr.F90 sudefo_gflattr.F90 sufa.F90 sugfl1.F90 sugfl2.F90 \
                  sugfl3.F90; do
        [ -f $EXT/$file ] && mvdiff $EXT/$file ../arpifs/setup/
      done
      for file in vpos_prep.F90; do
        [ -f $EXT/$file ] && mvdiff $EXT/$file ../arpifs/fullpos/
      done
      if [ $cycle == '49t2' ]; then
        file=$EXT/arp_shallow_mf.F90; [ -f $file ] && mvdiff $file ../arpifs/phys_dmn/
      fi
      #Special mpa case
      [ -f $EXT/modd_ch_aerosol.F90 ] && mvdiff $EXT/modd_ch_aerosol.F90 ../mpa/chem/module/modd_ch_aerosol.F90
      for file in modd_spp_type.F90 spp_mod_type.F90 aroini_conf.h aroini_conf.F90; do
        if [ -f $EXT/$file ]; then
          [ ! -d ../mpa/aux ] && mkdir ../mpa/aux
          mvdiff $EXT/$file ../mpa/aux/
        fi
      done
      [ -d $EXT/dead_code ] && rm -rf $EXT/dead_code/
      if [ $EXT == "PHYEX/externals" ]; then
        mv $EXT .
      else
        #Move automatically all codes under mpa
        if [ -d $HOMEPACK/$name/src/main ]; then
          #if main exists this is not a main pack (because it references a main pack)
          reftree='main'
        else
          reftree='local'
        fi
        if [ ! -z "$(\ls $EXT)" ]; then
          #$EXT is not empty
          for file in $EXT/*; do
            extname=`basename $file`
            loc=`find ../../$reftree/mpa/ -name $extname | sed "s/\/$reftree\//\/local\//g"`
            nb=`echo $loc | wc -w`
            if [ $nb -ne 1 ]; then
              echo "Don't know where $file must be moved, none or several places found!"
              exit 9
            fi
            mvdiff $file $loc
          done
        fi
      fi
    fi
    rm -rf PHYEX
  fi
fi

#####################
#### COMPILATION ####
#####################

if [ $compilation -eq 1 ]; then
  if [ $onlyIfNeeded -eq 0 -o ! -f $HOMEPACK/$name/bin/MASTERODB ]; then
    echo "### Compilation of commit $commit"

    cd $HOMEPACK/$name
    sed -i 's/GMK_THREADS=1$/GMK_THREADS=10/' ics_masterodb
  
    [ -f ics_packages ] && exescript Output_compilation_hub ics_packages
    exescript Output_compilation ics_masterodb
    if [ -f bin/MASTERODB \
         -a $(grep Error Output_compilation | \
              grep -v TestErrorHandler | \
              grep -v "'Error" | \
              grep -v "'CPLNG: Error" | \
              grep -v '"Error' | \
              grep -v "'*** Error" | \
              grep -v "\-\- Up-to-date:" | wc -l) -ne 0 ]; then
      echo "MASTERODB was produced but errors occured during compilation:"
      grep Error Output_compilation | \
            grep -v TestErrorHandler | \
            grep -v "'Error" | \
            grep -v "'CPLNG: Error" | \
            grep -v '"Error' | \
            grep -v "'*** Error" | \
            grep -v "\-\- Up-to-date:"
      echo "MASTERODB suppressed!"
      rm -f bin/MASTERODB
      exit 12
    fi
  fi
fi

###################
#### EXECUTION ####
###################

if [ $run -ge 1 ]; then
  #Cleaning to suppress old results that may be confusing in case of a crash during the run
  if [ $onlyIfNeeded -eq 0 ]; then
    for t in $(echo $tests | sed 's/,/ /g'); do
      cd $HOMEPACK/$name
      if [ -d conf_tests/$t ]; then
        rm -rf conf_tests/$t
      fi
    done
  fi

  #Run the tests one after the other
  firstrun=1
  for t in $(echo $tests | sed 's/,/ /g'); do #loop on tests
    if echo $allowedTests | grep -w $t > /dev/null; then #test is allowed on this plateform
      cd $HOMEPACK/$name
      if [ ! -d conf_tests/$t ]; then #We do not enter systematically this part if onlyIfNeeded=1
        if [ $firstrun -eq 1 ]; then
          echo "### Running of commit $commit"
          firstrun=0
        fi

        if [ ! -f $HOMEPACK/$name/bin/MASTERODB ]; then
          echo "Pack does not exist ($HOMEPACK/$name) or compilation has failed, please check"
          exit 6
        fi

        mkdir -p conf_tests/$t
        cd conf_tests/$t
        t1=$(($(date +%s%N)/1000)) #current time in milliseconds
        MYLIB=$name TESTDIR=$dirconf/$t exescript Output_run $dirconf/$t/ar?${cycle}${scripttag}.sh
        t2=$(($(date +%s%N)/1000))
        if [ "$perffile" != "" ]; then
          #The elapsed time is not relevant when the model runs with a queuing system (HPC)
          echo "$commit ial $t $(($t2-$t1))" >> "$perffile"
        fi

        #Profiling
        if ls drhook.prof.* > /dev/null 2>&1; then
          firstfile=1
          for file in drhook.prof.*; do
            if [ $firstfile -eq 1 ]; then
              cp $file drhook.prof.concat
              firstfile=0
            else
              #append only relevant lines
              grep -e '^ *[0-9]' $file >> drhook.prof.concat
            fi
          done
          firstLine=$(grep -m 1 -n "^ *1" drhook.prof.concat | cut -d: -f1)
          python3 -c "import numpy, pandas
d = {'time': ('<f4', ('mean', )), 'self': ('<f4', ('mean', 'max', 'min', 'std', 'sum')),
     'total': ('<f4', ('mean', 'max', 'min', 'std', 'sum')), 'calls': ('<i4', ('sum', )),
     'self_per_call': ('<f4', ('mean', )), 'total_per_call': ('<f4', ('mean', )), 'routine': ('U256', '')}
arraynp = numpy.loadtxt('drhook.prof.concat', dtype=[(k, v[0]) for (k, v) in d.items()],
                        converters={8: lambda s: s.split(b'@')[0].lstrip(b'*')},
                        skiprows=$firstLine - 1, usecols=[1, 3, 4, 5, 6, 7, 8], encoding='bytes')
df = pandas.DataFrame(arraynp).groupby('routine').agg(
      **{k + '_' + agg:pandas.NamedAgg(column=k, aggfunc=agg)
         for (k, agg) in [(k, agg) for k in d.keys() for agg in d[k][1]]
         if k != 'routine'}).sort_values('self_sum', ascending=False)
df.index.name += ' ordered by self_sum'
with open('drhook.prof.agg', 'w') as f: f.write(df.to_string())
"
        fi
      fi
    else
      echo "The test $t is not allowed"
    fi
  done
fi

####################
#### COMPARISON ####
####################
if [ $check -eq 1 ]; then
  echo "### Check commit $commit against commit $reference"

  allt=0
  message=""
  filestocheck=""
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then
      #Run the reference if needed
      if [ $computeRefIfNeeded -eq 1 ]; then
        $0 -p -c -r -t $t --onlyIfNeeded ${refByTest[$t]}
      fi

      #Files to compare
      if echo $t | grep 'small' > /dev/null; then
        if [ $cycle == '49t2' -o $cycle == '50t1' ]; then
          filestocheck="$filestocheck ${t},conf_tests/$t/ICMSHFPOS+0002:00"
        else
          filestocheck="$filestocheck ${t},conf_tests/$t/ICMSHFPOS+0002:00 ${t},conf_tests/$t/DHFDLFPOS+0002"
        fi
      else
        filestocheck="$filestocheck ${t},conf_tests/$t/NODE.001_01"
      fi
    else
      echo "The test $t is not allowed"
    fi
  done

  for tag_file in $filestocheck; do
      tag=$(echo $tag_file | cut -d, -f1)
      file=$(echo $tag_file | cut -d, -f2)
      refname=${refnameByTest[$tag]}
      ref=${refByTest[$tag]}
      file1=$HOMEPACK/$name/$file
      file2=$(echo $HOMEPACK/$refname/$file) #echo to enable shell substitution

      mess=""
      t=0
      if [ ! -f "$file1" ]; then
        mess="Result ($file1) for commit $commit does not exist, please run the simulation"
        t=1
      fi
      if [ ! -f "$file2" ]; then
        mess2="Reference result ($file2) for commit $ref does not exist, please run the simulation"
        t=1
        if [ "$mess" = "" ]; then
          mess=$mess2
        else
          mess="$mess and $mess2"
        fi
      fi
      if [ $t -eq 0 ]; then
        cmd="$PHYEXTOOLSDIR/compare.py"
        if [ ! -x $cmd ]; then
          echo "Command not found: \"$cmd\""
          exit 10
        fi
        if [ $(basename $file) == ICMSHFPOS+0002:00 ]; then
          #historic files
          cmd="$cmd --binary $file1 $file2 256"
        elif [ $(basename $file) == DHFDLFPOS+0002 ]; then
          #DDH files
          ddh_images="$HOMEPACK/$name/ddh_diff_${tag}.png"
          if [ `hostname` == 'sxphynh' ]; then
            [ ! -d /d0/images/$USER ] && mkdir /d0/images/$USER
            ddh_images="$ddh_images /d0/images/$USER/ddh_diff_${tag}.png"
          fi
          cmd="$cmd --ddh $file1 $file2 --ddhplots $ddh_images"
        elif [ $(basename $file) == NODE.001_01 ]; then
          #Output listing
          cmd="$cmd --node $file1 $file2"
        else
          cmd="$cmd --binary $file1 $file2 0"
        fi
        set +e
        mess=$($cmd)
        t=$?
        set -e
      fi
      [ $t -ne 0 ] && message="$message $file : $mess \n"
      allt=$(($allt+$t))
  done
  if [ $allt -eq 0 ]; then
    echo "SUCCESS, files are (nearly) identical"
  else
    echo "*************** Files are different *******************"
    echo -e "$message"
    cmpstatus=50
  fi
fi

##################
#### CLEANING ####
##################

if [ $remove -eq 1 ]; then
  echo "### Remove model directory for commit $commit"
  [ -d $HOMEPACK/$name ] && rm -rf $HOMEPACK/$name
fi

exit $cmpstatus
