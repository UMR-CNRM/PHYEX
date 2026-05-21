#!/bin/bash

#set -x
set -e
set -o pipefail #abort if left command on a pipe fails

#This script:
# - compiles the LMDZ model using a specific commit for the externalised physics
# - runs RICO and ARM-CU 1D cases

#######################
#### CONFIGURATION ####
#######################

#About the tests:
# - allowedTests is the list of allowed tests which can depend on platform, if we ask to perform an action
#   with a test not in the allowedTests list, the action is ignored
# - defaultTest is the list of tests to perform when no '-t' option is provided on the command line.
allowedTests="rico arm_cu"
defaultTest="rico"

separator='_' #- be carrefull, gmkpack (at least on belenos) has multiple allergies (':', '.', '@')
              #- seprator must be in sync with prep_code.sh separator

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

LMDZPACK=${LMDZPACK:=$HOME/LMDZ/PHYEX}

L=79

################################
#### COMMAND LINE ARGUMENTS ####
################################

function usage {
  echo "Usage: $0 [-h] [-p] [-u] [-c] [-C] [-r] [-s] [--expand] [-t test] [--version VERSION] [--repo-user] [--repo-protocol] [--remove] [--perf FILE] [--name NAME] commit [reference]"
  echo "commit          commit hash (or a directory) to test"
  echo "reference       commit hash (or a directory) REF to use as a reference"
  echo "-s              suppress compilation pack"
  echo "-p              creates pack"
  echo "-u              updates pack (experimental, only for 'small' updates where"
  echo "                              the list of files does not change)"
  echo "-c              performs compilation"
  echo "-r              runs the tests"
  echo "-C              checks the result against the reference"
  echo "-t              comma separated list of tests to execute"
  echo "                or ALL to execute all tests"
  echo "--nofcm         don't use fcm (be carreful, with fcm compilation exits with a (false) error"
  echo "--expand        expand mnh_expand blocks (code will be in do loops)"
  echo "--version VERSION   to force using lmdz VERSION"
  echo "--repo-user     user hosting the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOuser (=$PHYEXREOuser)"
  echo "--repo-protocol protocol (https or ssh) to reach the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOprotocol (=$PHYEXREOprotocol)"
  echo "--remove        removes the pack"
  echo "--onlyIfNeeded  performs the pack creation and/or the compilation and/or the execution"
  echo "                only if the step has not already been done"
  echo "--computeRefIfNeeded"
  echo "                computes the missing references"
  echo "--perf FILE     add performance statistics in file FILE"
  echo "--name NAME     to give a name to the pack"
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

fcm=1
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
version=""
link=0 #Not yet put in command line argument because this option has not been tested here
remove=0
onlyIfNeeded=0
computeRefIfNeeded=0
perffile=""
name=""

commitcmd="" # To execute again check_commit with the same options (except commit and ref)
while [ -n "$1" ]; do
  toadd="$1"
  case "$1" in
    '-h') usage; exit;;
    '-s') suppress=1;;
    '-p') packcreation=1;;
    '-u') packupdate=1;;
    '-c') compilation=1;;
    '-r') run=$(($run+1));;
    '-C') check=1;;
    '-t') tests="$2"; toadd="$toadd $2"; shift;;
    '--nofcm') fcm=0;;
    '--expand') useexpand=1;;
    '--version') version="$2"; toadd="$toadd $2"; shift;;
    '--repo-user') export PHYEXREPOuser=$2; toadd="$toadd $2"; shift;;
    '--repo-protocol') export PHYEXREPOprotocol=$2; toadd="$toadd $2"; shift;;
    '--remove') remove=1;;
    '--onlyIfNeeded') onlyIfNeeded=1;;
    '--computeRefIfNeeded') computeRefIfNeeded=1;;
    '--perf') perffile="$(realpath $2)"; toadd="$toadd $2"; shift;;
    '--name') name=$2; toadd="$toadd $2"; shift;;
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
        fi
        toadd="";;
  esac
  commitcmd="$commitcmd $toadd"
  shift
done

if [ $packcreation -eq 0 -a \
     $compilation -eq 0 -a \
     $packupdate -eq 0 -a \
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

#######################
#### FROM A COMMIT ####
#######################

# In this case, we clone the commit and run the check_commit script
# of this commit using the cloned directory.
#
if ! echo $commit | grep '/' | grep -v '^tags/' > /dev/null; then
  # We must run check_commit on a commit
  # We use a sub-shell for isolation
  (
    # Temporary directory
    export TMP_LOC=$(mktemp -d) # Workind directory
    trap "\rm -rf $TMP_LOC" EXIT # Automatic removing
    cd $TMP_LOC

    # Clone
    if [ $PHYEXREPOprotocol == 'https' ]; then
      repository=https://github.com/$PHYEXREPOuser/PHYEX.git
    elif [ $PHYEXREPOprotocol == 'ssh' ]; then
      repository=git@github.com:$PHYEXREPOuser/PHYEX.git
    fi
    git clone $repository
    cd PHYEX
    git checkout $commit
    if [ ! -d docs ]; then
      # LMDZ adapted commit, we must fill the repository with
      # the last version (in history) of tools and requirements
      toolscommit=$(git log --all --full-history -- tools | head -1 | cut -d\  -f2) # deleted in this commit
      for parent in $(git rev-parse ${toolscommit}^@); do
        if [ $(git ls-tree -r $parent -- tools | wc -l) -ne 0 ]; then
          toolscommit=$parent # last commit with tools
          break
        fi
      done
      git checkout $toolscommit -- "tools"
      git checkout $toolscommit -- "requirements.txt"
    fi

    # Requirements
    python3 -m venv venv
    . venv/bin/activate
    python3 -m pip install -r requirements.txt

    # Running commit
    . tools/env.sh
    if [ "$packBranch" == "" ]; then
      mypackname="--name $commit"
    else
      mypackname=""
    fi
    check_commit_lmdz.sh $commitcmd $TMP_LOC/PHYEX $reference $mypackname
  )
  exit $?
fi

###########################
#### COMMIT PROPERTIES ####
###########################

if [ -d $commit/src ]; then
  lmdz_ready=false
  content_lmdz_version=$(scp $commit/src/lmdz/lmdz_version.json /dev/stdout 2>/dev/null || echo "")
else
  lmdz_ready=true
  content_lmdz_version=$(scp $commit/lmdz_version.json /dev/stdout 2>/dev/null || echo "")
fi

[ "$version" == "" ] && version=$(json_dictkey2value "$content_lmdz_version" version '')
rad=$(json_dictkey2value "$content_lmdz_version" rad '')
install_arg=$(json_dictkey2value "$content_lmdz_version" install_arg '')

if [ -z "${tests-}" ]; then
  tests=$defaultTest
elif echo "$tests" | grep -w 'ALL' > /dev/null; then
  tests=$allowedTests
fi

if [ "$name" == "" ]; then
  name=$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
fi
lmdzdir=$LMDZPACK/$name/LMDZ
phyexdir=$LMDZPACK/$name/PHYEX/
main=lmdz1d
##-debug is the default value and there's a bug (in the current script) if we try to specify it here
##compilecmd="./compile -L $L -rad $rad -cosp 0 -opt \"-debug \" -main $main"
compilecmd="./compile -L $L -rad $rad -cosp 0 -main $main"

if [ ! -z "${reference-}" ]; then
  declare -A refnameByTest
  #Reference to use for each test
  for t in $(echo $allowedTests | sed 's/,/ /g'); do
    #Name of the reference
    if [ "$reference" == 'REF' ]; then
      if [[ ! -z "${refByTest[$t]+unset}" ]]; then
        #The json file contained the references to use on a per test case basis
        caseref=${refByTest[$t]}
      fi
    else
      #The exact reference to use was given on the command line
      caseref=$reference
    fi
    refByTest[$t]=$caseref
    #Conversion into directory name
    if echo $caseref | grep '/' > /dev/null; then
      refname="PHYEX/*_$(echo $caseref | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g').01.${gmkpack_l}.${gmkpack_o}"
    else
      refname="PHYEX/*_${caseref}.01.${gmkpack_l}.${gmkpack_o}"
    fi
    refnameByTest[$t]=$refname
  done
fi

#######################
#### PACK CREATION ####
#######################

[ $suppress -eq 1 -a -d $LMDZPACK/$name ] && rm -rf $LMDZPACK/$name
if [ $packcreation -eq 1 -a -d $LMDZPACK/$name -a $onlyIfNeeded -eq 1 ]; then
  packcreation=0
fi
if [ $packcreation -eq 1 ]; then
  echo "### Compilation of commit $commit"

  if [ -d $LMDZPACK/$name ]; then
    echo "Pack already exists ($LMDZPACK/$name), suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
    exit 5
  fi

  #Create directory
  mkdir -p $LMDZPACK/$name
  cd $LMDZPACK/$name

  #Populate with arch files
  touch arch-mylocal.env
  cat - <<EOF > arch-mylocal.fcm
%COMPILER            gfortran
%LINK                gfortran
%FPP                 cpp
%AR                  ar
%ARFLAGS             rU
%MAKE                make
%FPP_FLAGS           -P -traditional
%FPP_DEF             NC_DOUBLE
%BASE_FFLAGS          -cpp -ffree-line-length-0 -fdefault-real-8 -DNC_DOUBLE
%PROD_FFLAGS         -O3 -fallow-argument-mismatch
%DEV_FFLAGS          -Wall -fbounds-check  -fallow-argument-mismatch
%DEBUG_FFLAGS        -g3 -Wall -fbounds-check -ffpe-trap=invalid,zero,overflow -O0 -fstack-protector-all -fbacktrace -finit-real=snan  -fallow-argument-mismatch
%MPI_FFLAGS
%OMP_FFLAGS          
%BASE_LD              
%MPI_LD
%OMP_LD   
EOF

  cat - <<EOF > arch-mylocal.path
NETCDF_INCDIR="-I/usr/include"
NETCDF_LIBDIR=""
NETCDF_LIB="-lnetcdff -lnetcdf"

NETCDF95_INCDIR=-I\$LMDGCM/../../include/
NETCDF95_LIBDIR=-L\$LMDGCM/../../lib
NETCDF95_LIB=-lnetcdf95

IOIPSL_INCDIR="-I\$LMDGCM/../../lib -I\$LMDGCM/../IOIPSL/inc"
IOIPSL_LIBDIR="-L\$LMDGCM/../../lib -L\$LMDGCM/../IOIPSL/lib"
IOIPSL_LIB="-lioipsl"

XIOS_INCDIR="-I\$LMDGCM/../XIOS/inc"
XIOS_LIBDIR="-L\$LMDGCM/../XIOS/lib"
XIOS_LIB="-lxios -lstdc++"

ORCH_INCDIR="-I\$LMDGCM/../../lib"
ORCH_LIBDIR="-L\$LMDGCM/../../lib"

OASIS_INCDIR="-I\$LMDGCM/../../oasis3-mct/BLD/build/lib/psmile.MPI1"
OASIS_LIBDIR="-L\$LMDGCM/../../oasis3-mct/BLD/lib"
OASIS_LIB="-lpsmile.MPI1 -lscrip -lmct -lmpeu"

INCA_INCDIR="-I\$LMDGCM/../INCA/build/inc"
INCA_LIBDIR="-L\$LMDGCM/../INCA/build/lib"
INCA_LIB="-lchimie"
EOF

  #Compilation
  wget https://lmdz.lmd.jussieu.fr/pub/install_lmdz.sh -O install_lmdz.sh
  bash install_lmdz.sh -v $version $install_arg -bench 0 -rad $rad -name LMDZ -arch_dir $PWD -arch mylocal 2>&1 | tee Install.log

  #Populate with test cases (1D directory needed for compilation)
  cd $lmdzdir
  wget https://lmdz.lmd.jussieu.fr/pub/1D/1D.tar.gz
  tar xf 1D.tar.gz
fi
if [ $packupdate -eq 1 -o $packcreation -eq 1 ]; then
  #PHYEX code
  if [ $link -eq 1 ]; then
    if [ $packupdate -eq 1 ]; then
      echo "link option not compatible with the update option"
      exit 10
    fi
    #Special case when a PHYEX repository exist locally
    #We can link the LMDZ source tree with the PHYEX repository
    #This can be useful for debuging or developping
    cd ${lmdzdir}/modipsl/modeles/LMDZ/libf/phylmd/
    ln -s ~/PHYEX/src/common/*/* .
    ln -sf ~/PHYEX/src/lmdz/*/* .
  else
    #Checkout PHYEX
    if [ ! -d $LMDZPACK/$name ]; then
      echo "Pack directory doesn't exist ($LMDZPACK/$name)"
      exit 9
    fi
    cd $LMDZPACK/$name
    if [ $packupdate -eq 1 ]; then
      mv PHYEX PHYEXori
    fi
  
    if [ $useexpand == 1 ]; then
      expand_options="--mnhExpand"
    else
      expand_options=""
    fi
    subs="-s turb -s micro -s aux -s ext"
    prep_code=$PHYEXTOOLSDIR/prep_code.sh
    echo "Copy $commit"
    mkdir PHYEX
    scp -q -r $commit/src PHYEX/
    $prep_code $expand_options $subs -m lmdz PHYEX -- --shumanFUNCtoCALL --removeACC
    if [ $packupdate -eq 1 ]; then
      #Update only modified files
      cd PHYEX
      for file in $(find turb micro aux ext -type f); do
        mvdiff $file ../PHYEXori/$file
      done
      cd ..
      rm -rf PHYEX
      mv PHYEXori PHYEX
    else
      #Put PHYEX source code in the LMDZ source tree
      cd $lmdzdir/modipsl/modeles/LMDZ/libf/phylmd/
      ln -sf $phyexdir/*/* .
    fi
  fi

  #Update code
  cd $lmdzdir/modipsl/modeles/LMDZ/libf/phylmd/
  [ $packupdate -eq 0 ] && cp -r . ../phylmdorig
  if [ $fcm -eq 0 ]; then
    mvdiff modd_dimphyexn.F90 modd_dimphyex.F90
    for name in `grep -i 'END MODULE' modd*n.F90 | cut -d: -f1 | sed -e 's/n.F90//'` ; do
      mvdiff ${name}n.F90 ${name}_n.F90
    done
    mvdiff hypgeo.F90 modi_hypgeo.F90
    mvdiff hypser.F90 modi_hypser.F90
    mvdiff momg.F90 modi_momg.F90
    mvdiff tools.F90 mode_tools.F90
    [ -f shuman_mf.F90 ] && mvdiff shuman_mf.F90 modi_shuman_mf.F90
    [ -f shuman_phy.F90 ] && mvdiff shuman_phy.F90 mode_shuman_phy.F90
  fi

  #Missing files in case ecrad is not used
  if [ $packupdate -eq 0 ]; then
    if [ "$rad" != "ecrad" ] ; then
      for file in ecrad/yom* ecrad/abor1.F90 ecrad/abor1.intfb.h ecrad/parkind1.F90 \
                  ecrad/*/yom* ecrad/*/abor1.F90 ecrad/*/abor1.intfb.h ecrad/*/parkind1.F90; do
        [ -f "$file" -a ! -f "$(basename $file)" ] && ln -s $file .
      done
    fi
  fi
fi

#####################
#### COMPILATION ####
#####################

if [ $compilation -eq 1 ]; then
  if [ $onlyIfNeeded -eq 0 -o ! -f $HOMEPACK/$name/bin/MASTERODB ]; then
    echo "### Compilation of commit $commit"
    cd $lmdzdir/1D/bin
    rm -f $lmdzdir/1D/bin/${main}.e
    rm -f $lmdzdir/modipsl/modeles/LMDZ/bin/lmdz1d_${L}_phylmd_${rad}_seq.e
    if [ $fcm -eq 1 ]; then
      sed -i "s/fcm=0/fcm=1/g" compile
    fi
    $compilecmd 2>&1 | tee $lmdzdir/compilation.log
    if [ $fcm -eq 1 ]; then
      echo "Using fcm, compilation exits with error even if everything is OK"
    fi
  fi
fi

###################
#### EXECUTION ####
###################

if [ $run -eq 1 ]; then
  echo "### Execution of commit $commit"
  cd $lmdzdir/1D/INPUT/PHYS
  sed '1 i\iflag_physiq=1\n' physiq.def_6A > physiq.def_PHYLMD
  sed '1 i\iflag_physiq=2\n' physiq.def_6A > physiq.def_PHYEX
  for t in $tests; do
    for DEF in PHYEX PHYLMD; do
      d=${lmdzdir}/1D/EXEC/${DEF}L$L/$t
      [ ! -d $d ] && mkdir -p $d
      cd $d
      ln -sf ${lmdzdir}/1D/OLDCASES/$t/* .
      cp -f ${lmdzdir}/1D/INPUT/DEF/*.def .
      cp -f ${lmdzdir}/1D/INPUT/PHYS/physiq.def_$DEF physiq.def
      if [ $rad = oldrad ] ; then
        sed -i'' -e 's/iflag_rrtm=.*$/iflag_rrtm=0/' -e 's/NSW=.*$/NSW=2/' physiq.def
      fi
      if [ $rad = ecrad ] ; then
        cp -f $lmdzdir/modipsl/modeles/LMDZ/DefLists/namelist_ecrad .
        cp -rf $lmdzdir/modipsl/modeles/LMDZ/libf/phylmd/ecrad/data .
        sed -e 's@iflag_rrtm=1@iflag_rrtm=2@' physiq.def > tmp
        \mv tmp physiq.def
      fi
      cp -f ${lmdzdir}/1D/INPUT/VERT/L$L/* .
      ln -sf L$L.def vert.def
      set +e
      cp -f $lmdzdir/1D/OLDCASES/$t/*.d[ae]* .
      set -e
      cat <<......eod>| compile.sh
       cd $lmdzdir/1D/bin
       $compilecmd
......eod
      chmod +x compile.sh
  
      if [ $fcm -eq 0 ]; then
        ln -sf $lmdzdir/1D/bin/${main}.e ${main}.e
      else
        execname=$lmdzdir/modipsl/modeles/LMDZ/bin/lmdz1d_${L}_phylmd_${rad}_seq.e
        if [ -f $execname ]; then
          ln -sf $execname ${main}.e
        else
          ln -sf $lmdzdir/modipsl/modeles/LMDZ/lmdz1d_${rad}_L${L}.e ${main}.e
        fi
      fi
  
      if [ $DEF == PHYEX ]; then
        sed -i -e 's/day_step=144$/day_step=1440/' gcm1d.def
      fi
      t1=$(($(date +%s%N)/1000)) #current time in milliseconds
      ./lmdz1d.e 2>&1 | tee execution.log
      t2=$(($(date +%s%N)/1000))
      if [ "$perffile" != "" ]; then
        echo "$commit lmdz $t $(($t2-$t1))" >> "$perffile"
      fi
    done
  done
fi

if [ $check -eq 1 ]; then
  echo "### Check commit $commit against commit $reference"
  echo "This functionnality is not yet implemented because:"
  echo "  1) the PHYEX interface will evolve and bit-reproducibility is not guaranted"
  echo "  2) there is no surface scheme plugged in LMZD-PHYEX, a recompilation"
  echo "     must be done to change the surface fluxes"
  exit 6
fi

if [ $remove -eq 1 ]; then
  echo "### Remove model directory for commit $commit"
  [ -d $LMDZPACK/$name ] && rm -rf $LMDZPACK/$name
fi
