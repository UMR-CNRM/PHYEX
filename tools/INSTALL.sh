#!/bin/bash

#set -x
set -e
set -o pipefail #abort if left command on a pipe fails

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

function usage {
  echo "Usage: $0 [-h] [--install-dataset]"
  echo "--install-dataset          Install a reduced dataset for the test programs"
}

install_dataset=0

while [ -n "$1" ]; do
  case "$1" in
    '--install-dataset') install_dataset=1;;
  esac
  shift
done

if [ $install_dataset -eq 1 ]; then
  cd $PHYEXTOOLSDIR/testprogs_data/
  for file in https://github.com/UMR-CNRM/PHYEX/files/12783926/ice_adjust.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783935/rain_ice.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783942/rain_ice_old.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783945/shallow.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783952/turb.tar.gz; do
    basefile=$(basename $file)
    wget --no-check-certificate $file -O $basefile
    tar xf $basefile
    rm -f $basefile
  done
fi
