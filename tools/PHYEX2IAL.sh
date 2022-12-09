#!/bin/bash

set -e
set -x

################################
#Command line arguments and help

full_command="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}") $@"

function usage {
  echo "Usage: $0 [-h] [--phyex-repo-user PHYEXREPOuser] [--phyex-repo-protocol PHYEXREPOprotocol]"
  echo "       IALDIRECTORY:IALVERSION PHYEXVERSION BRANCH"
  echo
  echo "--phyex-repo-user PHYEXREPOuser         user hosting the PHYEX repository on github,"
  echo "                                        defaults to the env variable PHYEXREOuser (=$PHYEXREPOuser)"
  echo "--phyex-repo-protocol PHYEXREPOprotocol protocol (https or ssh) to reach the PHYEX repository on github,"
  echo "                                        defaults to the env variable PHYEXREOprotocol (=$PHYEXREPOprotocol)"
  echo "IALDIRECTORY                            local directory containing the IAL repository"
  echo "IALVERSION                              version to checkout in the IAL repository"
  echo "PHYEXVERSION                            commit, tag (as tags/<tag>) of PHYEX to use"
  echo
  echo "The scripts builds a pack using PHYEX (with the help of the check_commit_ial.sh script)"
  echo "and puts the content of the pack in the IAL repository."
  echo "It is important that the PHYEXVERSION is based on the same version as the selected IALVERSION."
}

IALDIRECTORY=''
IALVERSION=''
PHYEXVERSION=''
BRANCH=''

positional=0
while [ -n "$1" ]; do
  case "$1" in
    '-h') usage; exit;;
    '--phyex-repo-user') export PHYEXREPOuser="$2"; shift;;
    '--phyex-repo-protocol') export PHYEXREPOprotocol="$2"; shift;;
    *) positional=$(($positional + 1))
       case $positional in
         1) if echo "$1" | grep ':'; then
              IALDIRECTORY=$(echo "$1" | cut -d: -f1)
              IALVERSION=$(echo "$1" | cut -d: -f2-)
            else
              echo "First mandatory argument must take the form IALDIRECTORY:IALVERSION with ':' as separator"
              exit 1
            fi;;
         2) PHYEXVERSION="$1";;
         3) BRANCH="$1";;
         *) echo "Only three positional arguments are allowed"
            exit 2;;
       esac
  esac
  shift
done

if [ -z "${BRANCH-}" ]; then
  echo "This script needs positional arguments, you can use the -h option to get help"
  exit 3
fi

##################################################
#Create a gmkpack's pack and filling it with PHYEX

# Create temporary directory and set up its automatic destruction
export TMP_LOC=$(mktemp -d --tmpdir=$TMP XXXXXX)
trap "echo Removing now temporary directory $TMP_LOC ; \rm -rf $TMP_LOC" EXIT

# Creates a pack using check_commit_ial.sh script
echo "Creating pack in $TMP_LOC using $PHYEXVERSION PHYEX version with"
echo "PHYEXREPOuser=$PHYEXREPOuser and PHYEXREPOprotocol=$PHYEXREPOprotocol"
HOMEPACK=$TMP_LOC check_commit_ial.sh -p -f "${PHYEXVERSION}"

#########################################
#Create branch in IAL, fill it and commit

# Create the branch in the IAL repository
cd "${IALDIRECTORY}"
if [ ! -z "$(git status --porcelain)" ]; then
  echo "The IAL repository ($IALDIRECTORY) cannot be used as it is not clean"
  exit 4
fi
git checkout -b "${BRANCH}" "${IALVERSION}"

# copy the pack created with PHYEX into the IAL branch
cd $TMP_LOC/PHYEX/*/src/local/
for rep in *; do
  rsync -r --delete $rep /home/riette/IAL/
done

# commit
cd "${IALDIRECTORY}"
git add .
git commit -m "Integration of PHYEX version $PHYEXVERSION in IAL" \
           -m "PHYEXREPOuser=$PHYEXREPOuser PHYEXREPOprotocol=$PHYEXREPOprotocol $full_command"

