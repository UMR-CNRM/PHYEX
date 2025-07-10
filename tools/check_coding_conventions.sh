#!/bin/bash

#set -x
set -e
set -o pipefail #abort if left command on a pipe fails

#This script checks:
# - the presence of IMPLICIT NONE
# - that all dummy arguments are declared using the INTENT attribute

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

SOURCEDIR=$PHYEXTOOLSDIR/../src/common

################################
#### COMMAND LINE ARGUMENTS ####
################################

verbose=0
function usage {
  echo "Usage: $0 [-h] [-v]"
  echo "-v  print non conforming messages"
}
while [ -n "$1" ]; do
  case "$1" in
    '-h') usage; exit;;
    '-v') verbose=1;;
    *) echo "Too many arguments"; exit 3;;
  esac
  shift
done

################
#### CHECKS ####
################

JSONDIR=$(mktemp -d)
trap "\rm -rf $JSONDIR" EXIT
result=$(pyfortool_parallel --tree $SOURCEDIR --descTree $JSONDIR/tree.json \
                            --wrapH --enableCache \
                            --checkIMPLICIT Warn --checkINTENT Warn \
                            --checkPHYEXUnusedLocalVar Warn 2>&1)
if [ "$result" == "" ]; then
  exit 0
else
  if [ $verbose -eq 1 ]; then
    echo "$result"
  fi
  exit 2
fi
