#!/bin/bash

set -x
set -e
set -o pipefail #abort if left command on a pipe fails

#This script:
# - compiles the LMDZ model using a specific commit for the externalised physics
# - runs RICO and ARM-CU 1D cases

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
specialPack="ref"
availTests="rico arm_cu"
defaultTest="rico"
defaultRef='ref'
LMDZPACK=${LMDZPACK:=$HOME/LMDZ/PHYEX}
separator='_' #- be carrefull, gmkpack (at least on belenos) has multiple allergies (':', '.', '@')
              #- seprator must be in sync with prep_code.sh separator

function usage {
  echo "Usage: $0 [-h] [-p] [-c] [-C] [-r] [-s] [--expand] [-t test] [--version VERSION] [--repo-user] [--repo-protocol] commit [reference]"
  echo "commit          commit hash (or a directory, or among $specialPack) to test"
  echo "reference       commit hash (or a directory, or among $specialPack) REF to use as a reference"
  echo "-s              suppress compilation pack"
  echo "-p              creates pack"
  echo "-c              performs compilation"
  echo "-r              runs the tests"
  echo "-C              checks the result against the reference"
  echo "-t              comma separated list of tests to execute"
  echo "                or ALL to execute all tests"
  echo "--nofcm         don't use fcm (be carreful, with fcm compilation exits with a (false) error"
  echo "--expand        use mnh_expand (code will be in do loops)"
  echo "--version VERSION   to force using lmdz VERSION"
  echo "--repo-user     user hosting the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOuser (=$PHYEXREOuser)"
  echo "--repo-protocol protocol (https or ssh) to reach the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOprotocol (=$PHYEXREOprotocol)"
  echo ""
  echo "If nothing is asked (pack creation, compilation, running, check) everything is done"
  echo "If no test is aked for, the default one ($defaultTest) is executed"
  echo
  echo "With the special reference REF commit, a suitable reference is guessed"
  echo "The directory (for commit only, not ref) can take the form server:directory"
  echo "If using a directory (for commit or reference) it must contain at least one '/'"
  echo "The commit can be a tag, written with syntagx tags/<TAG>"
}

fcm=1
packcreation=0
compilation=0
run=0
check=0
commit=""
reference=""
tests=""
suppress=0
useexpand=0
version=""

while [ -n "$1" ]; do
  case "$1" in
    '-h') usage;;
    '-s') suppress=1;;
    '-p') packcreation=1;;
    '-c') compilation=1;;
    '-r') run=$(($run+1));;
    '-C') check=1;;
    '-t') tests="$2"; shift;;
    '--nofcm') fcm=0;;
    '--expand') useexpand=1;;
    '--version') version="$2"; shift;;
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

if [ -z "${tests-}" ]; then
  tests=$defaultTest
elif [ $tests == 'ALL' ]; then
  tests=$availTests
fi

if [ $packcreation -eq 0 -a \
     $compilation -eq 0 -a \
     $run -eq 0 -a \
     $check -eq 0 ]; then
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

function jsonarg {
  #$1 is the file containing the json dictionnary
  #$2 is the dictionnary key to return
  python3 -c "import json; print(json.load(open('$1', 'r'))['$2'])"
}

fromdir=''
if echo $commit | grep '/' | grep -v '^tags/' > /dev/null; then
  fromdir=$commit
  if [ "$version" == "" ]; then
    content_lmdz_version=$(scp $commit/src/lmdz/lmdz_version.json /dev/stdout 2>/dev/null)
    version=$(jsonarg <(echo $content_lmdz_version) version)
    rad=$(jsonarg <(echo $content_lmdz_version) rad)
    install_arg=$(jsonarg <(echo $content_lmdz_version) install_arg)
  fi
  name=$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
  [ $suppress -eq 1 -a -d $LMDZPACK/$name ] && rm -rf $LMDZPACK/$name
elif echo $specialPack | grep -w $commit > /dev/null; then
  name="PHYEX/$commit"
else
  if [ "$version" == "" ]; then
    if [[ $commit == lmdz${separator}* ]]; then
      lmdz_version_file="lmdz_version.json"
    else
      lmdz_version_file="src/lmdz/lmdz_version.json"
    fi
    if echo $commit | grep '^tags/' > /dev/null; then
      urlcommit=$(echo $commit | cut -d / -f 2-)
    else
      urlcommit=$commit
    fi
    content_lmdz_version=$(wget --no-check-certificate https://raw.githubusercontent.com/$PHYEXREPOuser/PHYEX/${urlcommit}/$lmdz_version_file -O - 2>/dev/null || echo "")
    version=$(jsonarg <(echo $content_lmdz_version) version)
    rad=$(jsonarg <(echo $content_lmdz_version) rad)
    install_arg=$(jsonarg <(echo $content_lmdz_version) install_arg)
  fi
  name="COMMIT$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')"
  [ $suppress -eq 1 -a -d $LMDZPACK/$name ] && rm -rf $LMDZPACK/$name
fi
if [ ! -z "${reference-}" ]; then
  [ $reference == 'REF' ] && reference=$defaultRef
  reffromdir=''
  if echo $reference | grep '/' > /dev/null; then
    reffromdir=$reference
    refname="PHYEX/$(echo $reference | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')"
  elif echo $specialPack | grep -w $reference > /dev/null; then
    refname="PHYEX/$reference"
  else
    refname="PHYEX/COMMIT${reference}"
  fi
fi

if [ $packcreation -eq 1 ]; then
  echo "### Compilation of commit $commit"

  if echo $specialPack | grep -w $commit > /dev/null; then
    echo "Special commit '$commit' cannot be compiled with this script"
    exit 4
  fi

  if [ -d $LMDZPACK/$name ]; then
    echo "Pack already exists ($LMDZPACK/$name), suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
    exit 5
  fi

  #Create directory
  cd $LMDZPACK
  mkdir $name
  cd $name
  base=$PWD
  wget https://lmdz.lmd.jussieu.fr/pub/install_lmdz.sh
  bash install_lmdz.sh -v $version $install_arg -bench 0 -rad $rad -name LMDZ > Install.log
  lmdzdir=$PWD/LMDZ

  #Populate with test cases
  cd $lmdzdir
  wget https://lmdz.lmd.jussieu.fr/pub/1D/1D.tar.gz
  tar xf 1D.tar.gz
  cd 1D
  rad=$(echo $rad) #to suppress spaces
  sed -i'' -e 's:^listecas=.*$:listecas="arm_cu rico":' -e "s/^rad=.*$/rad='$rad'/" run.sh
  
  cd INPUT/PHYS
  cp physiq.def_6A physiq.def_PHYLMD
  echo "iflag_physiq=1" >> physiq.def_PHYLMD
  sed -e "s/iflag_physiq=1/iflag_physiq=2/" physiq.def_PHYLMD > physiq.def_PHYEX

  #Update compilation cript
  phylmd=${lmdzdir}/modipsl/modeles/LMDZ/libf/phylmd/
  if [ $fcm -eq 1 ]; then
    sed -i "s/fcm=0/fcm=1/g" $lmdzdir/1D/bin/compile
  fi

  #Checkout PHYEX
  cd $base
  phyex=$base/PHYEX/
  MNH_EXPAND_DIR=$PHYEXTOOLSDIR/mnh_expand
  export PATH=$MNH_EXPAND_DIR/filepp:$MNH_EXPAND_DIR/MNH_Expand_Array:$PATH

  if [ $useexpand == 1 ]; then
    expand_options="-D MNH_EXPAND -D MNH_EXPAND_LOOP"
  else
    expand_options=""
  fi
  subs="-s turb -s micro -s aux -s ext"
  prep_code=$PHYEXTOOLSDIR/prep_code.sh
  if [ "$fromdir" == '' ]; then
    echo "Clone repository, and checkout commit $commit (using prep_code.sh)"
    if [[ $commit == lmdz${separator}* ]]; then
      $prep_code -c $commit PHYEX #This commit is ready for inclusion
    else
      $prep_code -c $commit $expand_options $subs -m arome PHYEX
    fi
  else
    echo "Copy $fromdir"
    mkdir PHYEX
    scp -q -r $fromdir/src PHYEX/
    $prep_code $expand_options $subs -m lmdz PHYEX
  fi

  #Update code
  cd $phylmd
  cp -r . ../phylmdorig
  ln -sf $phyex/*/* .
  if [ $fcm -eq 0 ]; then
    mv modd_dimphyexn.F90 modd_dimphyex.F90
    for name in `grep -i 'END MODULE' modd*n.F90 | cut -d: -f1 | sed -e 's/n.F90//'` ; do mv ${name}n.F90 ${name}_n.F90 ; done
    mv hypgeo.F90 modi_hypgeo.F90
    mv hypser.f90 modi_hypser.F90
    mv tools.F90 mode_tools.F90
    mv shuman_mf.F90 modi_shuman_mf.F90
    mv shuman_phy.F90 mode_shuman_phy.F90
  fi

  #Missing files in case ecrad is not used
  if [ "$rad" != "ecrad" ] ; then
    ln -s ecrad/yom* ecrad/abor1.F90 ecrad/abor1.intfb.h ecrad/parkind1.F90 .
  fi
fi

if [ $compilation -eq 1 ]; then
  echo "### Compilation of commit $commit"
  cd $LMDZPACK/$name/LMDZ/1D
  sed -i'' -e 's/^listedef=.*$/listedef="PHYLMD PHYEX"/' run.sh
  ./run.sh > log.$$ 2>&1
fi

if [ $run -eq 1 ]; then
  echo "### Execution of commit $commit"
  cd $LMDZPACK/$name/LMDZ/1D
  ./run.sh -r > log.$$ 2>&1
  wget https://www.lmd.jussieu.fr/~hourdin/phyex/compare.sh
  bash compare.sh
fi

