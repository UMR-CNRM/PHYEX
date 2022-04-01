#!/bin/bash

set -e

#This script can:
# - extract a tag or a commit from the PHYEX repository
# - merge code from common and model specific directories
# - apply mnh_expand tool
# - push the result in a new branch of the repository


###### CONFIGURATION
repository_https=https://github.com/QuentinRodier/PHYEX.git
repository_ssh=git@github.com:QuentinRodier/PHYEX.git

###### COMMAND LINE ARGUMENTS
function usage {
  echo "Usage: $0 [-h] [-c CHECKOUT_POINT] [-m MODEL] [-D OPTION [-D OPTION [...]]]] \\"
  echo "          [-s SUBDIR [-s SUBDIR [...]]]  [-v [-v [-v]]] DIRECTORY"
  echo "DIRECTORY             directory containing the script result"
  echo "-c CHECKOUT_POINT     git object to checkout, can be a specific commit"
  echo "                      or a tag with the following syntax: tags/TAG where TAG is the tag name"
  echo "-m MODEL              merge the code under the common directory with the code specific to MODEL model"
  echo "-D OPTION             option to use with mnh_expand"
  echo "                      BE CARREFULL, a space between -D and the option is required here"
  echo "-p                    push the result as a new branch"
  echo "-s SUB                subdiretory or file (under src) to consider when merging and applying mnh_expand"
  echo "--renameFf            rename .F90 into .f90"
  echo "--ssh                 use ssh instead of https for git cloning"
  echo "-v                    add verbosity (up to 3 -v)"
  echo ""
  echo "* If the -c option is not provided, DIRECTORY must already contain files and directory as if"
  echo "  it was the result of a git checkout"
  echo "* If the -m option is used, directory tree is modified, only relevant code is kept"
  echo "* If no -D options are used, mnh_expand is not called at all"
  echo "* -s options are mandatory for -m, -D and -p options"
  echo "* -p option is allowed only if -c and -m options are provided"
  echo ""
  echo "To use mnh_expand... it must be installed (alongside the filepp tool)"
}

full_command="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) $@"
separator='_' #- be carrefull, gmkpack (at least on belenos) has multiple allergies (':', '.', '@')
              #- seprator must be in sync with prep_code.sh separator

directory=""
checkout_point=""
mnh_expand_options=""
model=""
push=0
subs=""
renameFf=0
verbose=0
repository=$repository_https

while [ -n "$1" ]; do
  case "$1" in
    '-h') usage;;
    '-c') checkout_point="$2"; shift;;
    '-m') model="$2"; shift;;
    '-D') mnh_expand_options="$mnh_expand_options -D$2"; shift;;
    '-s') subs="$subs $2"; shift;;
    '-p') push=1;;
    '--renameFf') renameFf=1;;
    '--ssh') repository=$repository_ssh;;
    '-v') verbose=$(($verbose+1));;
     *) directory="$1";;
  esac
  shift
done

if [ $verbose -ge 3 ]; then
  set -x
fi

###### BRANCH OR NOT BRANCH
if [ -n "${checkout_point-}" -a -n "${model-}" -a $push == 1 ]; then
  branch=${model}${separator}${checkout_point}
fi

###### WORKING DIRECTORY
if [ -z "${directory-}" ]; then
  echo "A directory must be provided on command line (use -h option to get help)"
  exit 1
fi

if [ -z "${checkout_point-}" ]; then
  [ $verbose -gt 0 ] && echo "No checkout point provided, we use the content of $directory directory"
  if [ ! -d $directory/src ]; then
    echo "$directory must be filled with files and directories as if it was obtained through a checkout"
    exit 2
  fi
  cd $directory
  from='dir'
else
  [ $verbose -gt 0 ] && echo "Clone and checkout $checkout_point into $directory directory"
  if [ -d $directory ]; then
    echo "$directory already exists, suppress it before executing the script (or remove the -c option)"
    exit 3
  fi
  git clone $repository $directory
  cd $directory
  if [ -n "${branch-}" ]; then
    if [ $(git ls-remote --heads origin SR_GPU | wc -l) -eq 1 ]; then
      echo "$branch branch already exists on remote"
      exit 4
    fi
    branch="-b $branch"
  fi
  git checkout $branch $checkout_point
  from='git'
fi

###### MERGE
if [ -n "${model-}" ]; then
  if [ ! -d src/$model ]; then
    echo "src/$model directory does not exist"
    exit 5
  fi
  if [ -z "${subs-}" ]; then
    echo "It is not possible to merge common and model specific codes if no subs are provided"
    exit 6
  fi
  [ $verbose -gt 0 ] && echo "Merge common code and $model specific code"

  files=$(\ls -A) #files to suppress at the end

  #Merge
  for sub in $subs; do
    [ $verbose -gt 1 ] && echo "Merging $sub directory/file"
    if [ -e $sub ]; then
      echo "$sub must not exist in the repository, this is a limitation of the script"
       exit 7
    fi
    [ -e src/common/$sub ] && mv src/common/$sub .
    [ -e src/$model/$sub ] && cp -rlf src/$model/$sub . && rm -rf src/$model/$sub
  done

  #cleaning
  [ $verbose -gt 0 ] && echo "Cleaning unrelevant files"
  #multiple checks to prevent error
  if [ $from == 'git' -a ! $(git config --get remote.origin.url) == "$repository" ]; then
    echo "Not inside the right git!!!!!!!!!!!!!!!!"
    exit 8
  fi
  for file in $files; do
    if [ $from == 'dir' -o $(git ls-files --error-unmatch $file 2>/dev/null | wc -l) -gt 0 ] ; then
      [ $verbose -gt 1 ] && echo "Suppression of $file"
      rm -rf $file
    fi
  done
fi

###### MNH_EXPAND
if [ -n "${mnh_expand_options-}" ]; then
  [ $verbose -gt 0 ] && echo "Applying mnh_expand"
  function apply_mnh_expand () {
    if grep mnh_expand $1 > /dev/null 2>&1 ; then
      [ $verbose -gt 1 ] && echo "Applying mnh_expand on $1"
      mnh_expand -DMNH_EXPAND_NOCPP $mnh_expand_options $1 > tempo_mnh_expand
      mv tempo_mnh_expand $1
    fi
  }
  if [ -n "${model-}" ]; then
    reps=$subs
  else
    reps=""
    for sub in $subs; do
      reps="$reps src/*/$sub"
    done
  fi
  for rep in $reps; do
    if [ -d $rep ]; then
      #find $rep -type f | while read file; do
      find $rep -type f -not -name '.*.swp' | while read file; do
        apply_mnh_expand "$file"
      done
    fi
  done
fi

###### RENAME .F90 into .f90
if [ $renameFf -eq 1 ]; then
  find . -type f  -name \*.F90 -print0 | \
    while IFS= read -r -d '' file; do
      mv -- "$file" "${file%.F90}.f90"
    done
fi

###### PUSH
if [ -n "${branch-}" ]; then
  [ $verbose -gt 0 ] && echo "commit and push"
  git add -A
  git commit -m "Version '$checkout_point' of source code ready for inclusion into $model source tree" -m "$full_command"
  git push -u origin HEAD
fi

[ $verbose -gt 0 ] && echo "Finished!"
exit 0
