#!/bin/bash

#set -x
set -e
set -o pipefail #abort if left command on a pipe fails

pyft_version=50c0efa5877feec73cff2b1ce73f197d52417b51

#This script installs PHYEX
#Call the script with the -h option to get more information.

################################
#### COMMAND LINE ARGUMENTS ####
################################

function usage {
  echo "Usage: $0 [-h] [--ALL] [--dataset] [--pyft] [--clean]"
  echo "  --ALL              Install or clean everything"
  echo "  --dataset          Install or clean a reduced dataset for the test programs"
  echo "  --pyft             Install or clean the pyft tool"
  echo "  --fiatfcm          Install or clean the fiat and fcm tools"
  echo "  --clean            Clean instead of installing"
  echo "  --test             Perform a test"
  echo "  --ssh              Use the ssh protocol to clone the pyft, fxtran, fiat and fcm"
  echo "                     repositories instead of https"
  echo ""
  echo "If the installation has already been done, calling again this script will update"
  echo "the different tools."
}

ALL=0
dataset=0
pyft=0
clean=0
dotest=0
ssh=0
fiatfcm=0

while [ -n "$1" ]; do
  case "$1" in
    '--ALL') ALL=1;;
    '--dataset') dataset=1;;
    '--pyft') pyft=1;;
    '--fiatfcm') fiatfcm=1;;
    '--clean') clean=1;;
    '--test') dotest=1;;
    '--ssh') ssh=1;;
    '-h') usage; exit;;
    *) echo "Unknown option $1"; exit 1;;
  esac
  shift
done

if [ $ALL == 1 ]; then
  dataset=1
  pyft=1
  fiatfcm=1
fi

#################################
#### INSTALLATION / CLEANING ####
#################################

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $dataset -eq 1 ]; then
  cd $PHYEXTOOLSDIR/testprogs_data/
  for file in https://github.com/UMR-CNRM/PHYEX/files/12783926/ice_adjust.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783935/rain_ice.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783942/rain_ice_old.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783945/shallow.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783952/turb.tar.gz; do
    basefile=$(basename $file)
    if [ $clean -eq 1 ]; then
      find $(basename $basefile .tar.gz) -name \*.dat -type f -delete
    elif [ ! -f $(basename $basefile .tar.gz)/00000000.dat ]; then
      wget --no-check-certificate $file -O $basefile
      tar xf $basefile
      rm -f $basefile
    fi
  done
fi

if [ $pyft -eq 1 ]; then
  cd $PHYEXTOOLSDIR/site
  if [ $clean -eq 1 ]; then
    rm -rf pyft
  else
    if [ ! -d pyft ]; then
      #Install
      if [ $ssh -eq 1 ]; then
        git clone git@github.com:UMR-CNRM/pyft.git
      else
        git clone https://github.com/UMR-CNRM/pyft.git
      fi
    else
      #Update
      cd pyft
      if ! git rev-parse ${pyft_version}^{commit} > /dev/null 2>&1; then
        if [ $ssh -eq 1 ]; then
          git fetch git@github.com:UMR-CNRM/pyft.git
        else
          git fetch https://github.com/UMR-CNRM/pyft.git
        fi
      fi
      if [ $(git rev-parse HEAD^{commit}) != $(git rev-parse ${pyft_version}^{commit}) ]; then
        git checkout ${pyft_version}
      fi
      cd ..
    fi
    #Install/update
    if [ $ssh -eq 1 ]; then
      ./pyft/bin/INSTALL.sh --ssh
    else
      ./pyft/bin/INSTALL.sh
    fi
    echo ""
    echo "To use PHYEX, you only need to source the file $PHYEXTOOLSDIR/env.sh"
  fi
fi

if [ $fiatfcm -eq 1 ]; then
  if [ $clean -eq 1 ]; then
    $PHYEXTOOLSDIR/../build/with_fcm/make_fcm.sh --inplace-clean
  else
    #Install/update
    if [ $ssh -eq 1 ]; then
      $PHYEXTOOLSDIR/../build/with_fcm/make_fcm.sh --inplace-install --ssh
    else
      $PHYEXTOOLSDIR/../build/with_fcm/make_fcm.sh --inplace-install
    fi
  fi
fi

if [ $dotest -eq 1 ]; then
  . $PHYEXTOOLSDIR/env.sh
  if [ $ssh -eq 1 ]; then
    protocol="ssh"
  else
    protocol="https"
  fi
  check_commit_testprogs.sh $(cd $PHYEXTOOLSDIR/.. && pwd) REF \
                            --computeRefIfNeeded \
                            --repo-user UMR-CNRM --repo-protocol $protocol
fi
