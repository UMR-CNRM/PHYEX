#!/bin/bash

#set -x
set -e
set -o pipefail #abort if left command on a pipe fails

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

function usage {
  echo "Usage: $0 [-h] [--ALL] [--dataset] [--pyft] [--clean]"
  echo "  --ALL              Install or clean everything"
  echo "  --dataset          Install or clean a reduced dataset for the test programs"
  echo "  --pyft             Install or clean the pyft tool"
  echo "  --clean            Clean instead of installing"
}

ALL=0
dataset=0
pyft=0
clean=0

while [ -n "$1" ]; do
  case "$1" in
    '--ALL') ALL=1;;
    '--dataset') dataset=1;;
    '--pyft') pyft=1;;
    '--clean') clean=1;;
    '-h') usage;;
  esac
  shift
done

if [ $ALL == 1 ]; then
  dataset=1
  pyft=1
fi

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
    else
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
    git clone git@github.com:UMR-CNRM/pyft.git
    pyft/bin/INSTALL.sh
    echo "Instead of sourcing the previous file, you can source $PHYEXTOOLSDIR/env.sh"
    echo "This file, among other things, source the previous displayed file."
  fi
fi
    
