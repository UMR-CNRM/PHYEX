#!/bin/bash

set -e
#set -x

#This script can:
# - extract a tag or a commit from the PHYEX repository
# - merge code from common and model specific directories
# - apply the pyft tool
# - push the result in a new branch of the repository


###### CONFIGURATION
PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

###### COMMAND LINE ARGUMENTS
function usage {
  echo "Usage: $0 [-h] [-c CHECKOUT_POINT] [-m MODEL] [-D OPTION [-D OPTION [...]]]] \\"
  echo "          [-s SUBDIR [-s SUBDIR [...]]] [--pyft_opts_env VAR] [-v [-v [-v]]] \\"
  echo "          DIRECTORY -- PYFT_OPTIONS"
  echo "DIRECTORY             directory containing the script result"
  echo "-c CHECKOUT_POINT     git object to checkout, can be a specific commit"
  echo "                      or a tag with the following syntax: tags/TAG where TAG is the tag name"
  echo "-m MODEL              merge the code under the common directory with the code specific to MODEL model"
  echo "--mnhExpand           option passed to the pyft tool"
  echo "-p                    push the result as a new branch"
  echo "-s SUB                subdiretory or file (under src) to consider when merging and applying pyft"
  echo "--renameFf            rename .F90 into .f90"
  echo "--ilooprm             replace indexes in do loop (and mnh_expand) by :"
  echo "--repo                use this repository instead of the one derived (if any) from the env variables"
  echo "                      PHYEXREPOuser (=$PHYEXREPOuser) and PHYEXREPOprotocol (=$PHYEXREPOprotocol)"
  echo "-v                    add verbosity (up to 3 -v)"
  echo "--pyft_opts_env VAR   name of an environment variable containing options to use to call"
  echo "                      the pyft_tool.py script"
  echo "-- PYFT_OPTIONS       everything after '--' are used as options for pyft_tool.py"
  echo "                      These options are used for all the files."
  echo ""
  echo "* If the -c option is not provided, DIRECTORY must already contain files and directory as if"
  echo "  it was the result of a git checkout"
  echo "* If the -m option is used, directory tree is modified, only relevant code is kept"
  echo "* If --mnhExpand is not used, pyft is not called at all"
  echo "* -s options are mandatory for -m, -D and -p options"
  echo "* -p option is allowed only if -c and -m options are provided"
  echo ""
  echo "Everything after the '--' is passed to pyft for source-to-source transformation"
  echo ""
  echo "To use the pyft tool... it must be installed"
  echo ""
  echo "The variable name sent with --pyft_opts_env must correspond to an exported environement"
  echo "variable.  The variable can contain a multi-lines string."
  echo "The variable is read line by line and the last applicable line is used."
  echo "A line can take one of these two forms:"
  echo "  - FILE_DESCRIPTOR:OPTIONS"
  echo "    where FILE_DESCRIPTOR is a regular expression to test against the filename. If there"
  echo "    is a match, the OPTIONS can be used for the file. The regular expression is"
  echo "    tested using 'grep -e'."
  echo "  - OPTIONS"
  echo "    If the line doesn't contain the FILE_DESCRIPTOR part, it applies to all source code."
  echo ""
  echo "For example, to transform all source code in lower case:"
  echo "> export OPTS='--lowerCase'; $0 --pyft\_opts\_env OPTS ..."
  echo ""
  echo "To transform all source code in lower case, except routines in turb directory which must be"
  echo "in upper case but keeping the turb.F90 in lower case:"
  echo "> export OPTS='--lowerCase "
  echo "> ^turb/:--upperCase "
  echo "> ^turb/turb\..90:--lowerCase'; $0 --pyft\_opts\_env OPTS ..."
}

full_command="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}") $@"
separator='_' #- be carrefull, gmkpack (at least on belenos) has multiple allergies (':', '.', '@')
              #- seprator must be in sync with prep_code.sh separator

directory=""
checkout_point=""
pyft_options=""
model=""
push=0
subs=""
renameFf=0
verbose=0
ilooprm=0
forpyft=0
pyft_opts_env=""

if [ -z "${PHYEXREPOprotocol-}" ]; then
  repository=""
else
  if [ $PHYEXREPOprotocol == 'https' ]; then
    repository=https://github.com/$PHYEXREPOuser/PHYEX.git
  elif [ $PHYEXREPOprotocol == 'ssh' ]; then
    repository=git@github.com:$PHYEXREPOuser/PHYEX.git
  else
    repository=""
  fi
fi

while [ -n "$1" ]; do
  case "$1" in
    '-h') usage; exit;;
    '-c') checkout_point="$2"; shift;;
    '-m') model="$2"; shift;;
    '--mnhExpand') pyft_options="$pyft_options $1";;
    '-s') subs="$subs $2"; shift;;
    '-p') push=1;;
    '--renameFf') renameFf=1;;
    '--ilooprm') ilooprm=1;;
    '--repo') repository=$2; shift;;
    '-v') verbose=$(($verbose+1));;
    '--pyft_opts_env') pyft_opts_env=$2; shift;;
    '--') forpyft=1;;
     *) if [ $forpyft -eq 0 ]; then
          directory="$1"
        else
          pyft_options="$pyft_options $1"
        fi;;
  esac
  shift
done

if [ "$pyft_opts_env" != "" ]; then
  #pyft_opts_env contains the name of the environment variable to use
  pyft_opts_env=${!pyft_opts_env}
  #now, pyft_opts_env contains the configuration to use
fi
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
  mv='mv -f'
  rm='rm -f'
else
  [ $verbose -gt 0 ] && echo "Clone and checkout $checkout_point into $directory directory"
  if [ -d $directory ]; then
    echo "$directory already exists, suppress it before executing the script (or remove the -c option)"
    exit 3
  fi
  if [ -z "${repository-}" ]; then
    echo "A repository must be set (use -h option to get help)"
    exit 1
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
  mv='git mv -f'
  rm='git rm -q -f'
fi

###### RENAME .F90 into .f90
#This step could also be achieved by pyft_tool.py but it must be done *before* the call to pyft_tool
if [ $renameFf -eq 1 ]; then
  #we use find/while/read in case the number of files is too big to be hold on a single shell line
  find . -type f  -name \*.F90 -print0 | \
    while IFS= read -r -d '' file; do
      $mv "$file" "${file%.F90}.f90"
    done
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
      echo "$sub must not exist in the repository root, this is a limitation of the script"
      exit 7
    fi
    [ -e src/common/$sub ] && $mv src/common/$sub . #sub doesn't exist, we can move it directly
    if [ -e src/$model/$sub ]; then
      if [ -f src/$model/$sub ]; then
        #$sub is a file, it can be overwritten
        $mv "src/$model/$sub" $sub
      else
        #directory can exist, we must move files one by one
        #we use find/while/read in case the number of files is too big to be hold on a single shell line
        (cd src/$model/$sub; find . -type f -print0) | \
          while IFS= read -r -d '' file; do
            dname=$(dirname $file)
            [ ! -d $sub/$dname ] && mkdir -p $sub/$dname
            $mv "src/$model/$sub/$file" "$sub/$file"
          done
        rmdir --ignore-fail-on-non-empty -p "src/$model/$sub" #suppress tree if empty
      fi
    fi
  done

  #Supression of unwanted files
  if [ -f src/$model/filesToSuppress.txt ]; then
    #Some files can be present in the common directory but are not wanted for a model export
    #because these files are already existing elsewhere in the model source code
    while read -r line; do
      filename=$(echo $line | sed -e 's/^[[:space:]]*//' | sed -e 's/[[:space:]]*$//') #trim
      [ -f "$filename" ] && $rm "$filename"
    done < src/$model/filesToSuppress.txt
  fi

  #Cleaning
  [ $verbose -gt 0 ] && echo "Cleaning unrelevant files"
  #multiple checks to prevent error
  if [ $from == 'git' -a ! "$(git config --get remote.origin.url)" == "$repository" ]; then
    echo "Not inside the right git!!!!!!!!!!!!!!!!"
    exit 8
  fi
  for file in $files; do
    [ $verbose -gt 1 ] && echo "Suppression of $file"
    if [ -e "$file" -a "$file" != '.git' -a "$file" != '.git/' ]; then
      $rm -r "$file"
    fi
  done
fi

##### Replace index in do loop and mnh_expand directives by :
if [ $ilooprm -eq 1 ]; then
  subs=$(\ls)
  for sub in $subs; do
    cd $sub
    files=$(\ls -A)
    for file in $files; do
      if [[ "$file" != "gradient_m"* ]]; then
        # Protection only for one line in turb.f90/.F90
        if [[ "$file" == "turb"* ]]; then
          sed -i 's/PLM(IIJB:IIJE,IKTB:IKTE) = PZZ(IIJB:IIJE,IKTB+IKL:IKTE+IKL) - PZZ(IIJB:IIJE,IKTB:IKTE)/PLM(IIJB:IIJE,IKTB : IKTE) = PZZ(IIJB:IIJE,IKTB+IKL:IKTE+IKL) - PZZ(IIJB:IIJE,IKTB : IKTE)/g' $file
        fi
        # Protection
        sed -i 's/JK=IKTB:IKTE/transJKIKTB/g' $file
        sed -i 's/JK=1:IKT/transIKT/g' $file
        sed -i 's/JIJ=IIJB:IIJE/transJIJ/g' $file
        sed -i 's/IKTB+1:IKTE/IKTB1IKTE/g' $file
        # Apply transformation
        sed -i 's/1:IKT/:/g' $file
        sed -i 's/IKTB:IKTE/:/g' $file
        sed -i 's/IIJB:IIJE/:/g' $file
        # Supression protection
        sed -i 's/transJKIKTB/JK=IKTB:IKTE/g' $file
        sed -i 's/transIKT/JK=1:IKT/g' $file
        sed -i 's/transJIJ/JIJ=IIJB:IIJE/g' $file
        sed -i 's/IKTB1IKTE/IKTB+1:IKTE/g' $file
        if [[ "$file" == "turb"* ]]; then
          sed -i 's/IKTB : IKTE/IKTB:IKTE/g' $file
        fi
      fi
    done
    cd ..
  done
fi

###### PYFT
if [ "$pyft_opts_env" != "" -o -n "${pyft_options-}" ]; then
  [ $verbose -gt 0 ] && echo "Applying pyft_tool"

  #Update PATH and PYTHONPATH if needed
  which pyft_tool.py > /dev/null || . $PHYEXTOOLSDIR/site/pyft/bin/env.sh

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
      find $rep -type f -not -name '.*.swp'  -not -name '.*.swo' | while read file; do
        if [ "$(echo $file | grep '\.')" != '' -a $(echo $file | rev | cut -d. -f1 | rev) != 'fypp' ]; then
          #Files without extension are certainly not source code files
          #.fypp files cannot be read by pyft_tool.py
          extra_opts=""
          if [ "$pyft_opts_env" != "" ]; then
            while read line; do
              if echo $line | grep ':' > /dev/null; then
                #This line has the form FILE_DESCRIPTOR:OPTIONS
                fd=$(echo $line | cut -d: -f1)
                if echo $file | grep -e $fd > /dev/null; then
                  extra_opts=$(echo $line | cut -d: -f2-)
                fi
              else
                extra_opts=$line
              fi
            done < <(echo "$pyft_opts_env")
          fi
          if [ "$extra_opts" != "" -o -n "${pyft_options-}" ]; then
            cmd="pyft_tool.py --wrapH $pyft_options $extra_opts" #--wrapH allows to deal with h files
            [ $verbose -gt 1 ] && echo $cmd "$file"
            $cmd "$file"
          fi
        fi
      done
    fi
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
