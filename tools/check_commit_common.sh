# Common codes for check_commit_ scripts

#set -x
set -e
set -o pipefail #abort if left command on a pipe fails

###################
#### VARIABLES ####
###################

separator='_' #- be carrefull, gmkpack (at least on belenos) has multiple allergies (':', '.', '@')
              #- separator must be in sync with prep_code separator

export PHYEXCONF=${PHYEXCONF:-${HOME}/.phyex}

########################
#### MISC FUNCTIONS ####
########################

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
    # $2 is a directory and file must be moved in this directory, keeping its name
    dest=$2/$(basename $1)
  else
    dest=$2
  fi
  if [ ! -f $dest ]; then
    mv $1 $2
  else
    # File is moved only if different
    if ! cmp $1 $dest > /dev/null; then
      mv $1 $2
    else
      rm $1
    fi
  fi
}

function escape_commit {
  echo $1 | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g'
}

######################
#### COMMAND LINE ####
######################

function command_line {
  # Some variables can be set to control the argument parsing
  # default_expand=true|false: the defaut value for the expand option
  # default_buildSystem=''|option1 '' to disable the build system choice
  #                                option1 to set the default build system; in this case
  #                                         alternative_buildSystem must contain option2
  # enable_full=true|false to enable the -f option
  # enable_testprogs_opt=true|false
  # enable_prepCodeOpts=true|false
  # extra_doc can contain lines of documentation to print at the end of the message issued by -h

  # Default values for function options
  enable_full=${enable_full:-false}
  enable_testprogs_opt=${enable_testprogs_opt:-false}
  enable_prepCodeOpts=${enable_prepCodeOpts:-false}
  if [ "$default_buildSystem" == "" ]; then
    enable_buildSystem=false
  else
    enable_buildSystem=true
  fi

  # Default values for options
  packcreation=false
  packupdate=false
  compilation=false
  run=false
  check=false
  remove=false
  tests=""
  suppress=false
  [ $enable_full == true ] && fullcompilation=false
  useexpand=$default_expand
  commit=""
  reference=""
  onlyIfNeeded=false
  computeRefIfNeeded=false
  [ $enable_prepCodeOpts == true ] && prepCodeOpts=""
  perffile=""
  name=""
  [ $enable_buildSystem == true ] && buildSys=$default_buildSystem
  if [ $enable_testprogs_opt == true ]; then
    archfile=$defaultarchfile
    refarchfile=$defaultarchfile
    perf=true
    extrapolation=0
    checkOpt="--check"
  fi

  function usage {
    cmd="Usage: $0 [-h] [-s] [-p] [-u] [-c] [-r] [-C] [--repo-user USER] [--repo-protocol PROTOCOL] [--remove] [--onlyIfNeeded] [--computeRefIfNeeded] [--perf FILE] [--name NAME]"
    [ $enable_full == true ] && cmd="$cmd [-f]"
    [ $enable_buildSystem == true ] && cmd="$cmd [--no$default_buildSystem]"
    if [ $default_expand == true ]; then
       echo "--noexpand     do not expand mnh_expand blocks (code will be in array-syntax)"
    else
       echo "--expand       expand mnh_expand blocks"
    fi
    cmd="$cmd commit [reference]"
    echo "commit          commit hash or a directory to test"
    echo "                The directory can take the form server:directory"
    echo "reference       commit hash or a directory or 'REF' to use as a reference"
    echo "                With the special reference REF commit, a suitable reference is guessed"
    echo "-s              suppress compilation directory"
    echo "-p              create compilation directory"
    echo "-u              update compilation directory (experimental, only for 'small' updates where"
    echo "                the list of files does not change)"
    echo "-c              performs compilation"
    echo "-r              runs the tests"
    echo "-C              checks the result against the reference"
    echo "-t TEST         comma separated list of tests to execute"
    echo "                or ALL to execute all tests"
    [ $enable_full == true ] && echo "-f              full compilation (do not use pre-compiled pack)"
    [ $enable_buildSystem == true ] && echo "--no$default_buildSystem      Do not use $default_buildSystem"
    if [ $default_expand == true ]; then
       echo "--noexpand     do not expand mnh_expand blocks (code will be in array-syntax)"
    else
       echo "--expand       expand mnh_expand blocks"
    fi
    echo "--repo-user USER"
    echo "                user hosting the PHYEX repository on github,"
    echo "                defaults to the env variable PHYEXREOuser (=$PHYEXREOuser)"
    echo "--repo-protocol PROTOCOL"
    echo "                protocol (https or ssh) to reach the PHYEX repository on github,"
    echo "                defaults to the env variable PHYEXREOprotocol (=$PHYEXREOprotocol)"
    echo "--remove        removes the compilation directory"
    echo "--onlyIfNeeded  performs the directory creation and/or the compilation and/or the execution"
    echo "                only if the step has not already been done"
    echo "--computeRefIfNeeded"
    echo "                computes the missing references"
    echo "--perf FILE     add performance statistics in file FILE"
    echo "--name NAME     to give a specifiq name to the directory"
    echo "--prep_code-opts 'OPTS'"
    echo "                OPTS is added to the call to prep_code (e.g. --prep_code_opts '--lowerCase'"
    echo "                to transfor all source codes in lower case). Help on prep_code.sh options"
    echo "                can be found with 'prep_code.sh -h'. Note: don't forget to enclose OPTS in ' or \""
    if [ $enable_testprogs_opt == true ]; then
      echo "--no-perf       deactivate DR_HOOK"
      echo "--no-check      suppress value printing (comparison will be impossible)"
      echo "                this option can reduce drastically the running time but only allow"
      echo "                to access performance statistics."
      echo "-a arch ARCH    architecture name to use to build and run the commit (=$defaultarchfile)"
      echo "-A arch ARCH    architecture name to use for the reference simulation (=$defaultarchfile)"
      echo "-e EXTRAPOLATION"
      echo "                extrapolate data. EXTRAPOLATION corresponds to a configuration:"
      for i in $(seq 1 $((${#conf_extra_tag[@]}-1))); do
        echo "                  - '$i': ${conf_extra_opts[$i]} (${conf_extra_tag[$i]})"
      done
    fi
    echo ""
    echo "If nothing is asked (directory creation, compilation, running, check, removing) everything"
    echo "except the removing is done"
    echo ""
    echo "If no test is aked for, the default one ($defaultTest) is executed"
    echo ""
    echo "If using a directory (for commit or reference) it must contain at least one '/'."
    echo ""
    echo "Architecture and/or data can be added in the \${PHYEXCONF} (=${PHYEXCONF}) directory."
    echo ""
    echo "The commit can be a tag, written with syntagx tags/<TAG>"
    if [ $enable_full == true ]; then
      echo ""
      echo "The -f flag (full recompilation) is active only at pack creation"
    fi
    echo $extra_doc
  }

  function check_allowed {
    # print error message if $2 is false
    if [ $2 == false ]; then
      echo "$1 option forbidden for $0"
      exit 20
    fi
  }

  commitcmd="" # To execute again check_commit with the same options (except commit and ref)
  while [ -n "$1" ]; do
    toadd="$1"
    case "$1" in
      '-h') usage; exit;;
      '-p') packcreation=true;;
      '-u') packupdate=true;;
      '-c') compilation=true;;
      '-r') run=true;;
      '-C') check=true;;
      '--remove') remove=true;;
      '-t') tests="$2"; toadd="$toadd $2"; shift;;
      '-s') suppress=true;;
      '-f') if [ $enable_full == false ]; then
              echo "-f option forbidden for $0"
              exit 20
            fi
            fullcompilation=true;;
      '--noexpand') useexpand=false;;
      '--expand') useexpand=true;;
      '--repo-user') export PHYEXREPOuser=$2; toadd="$toadd $2"; shift;;
      '--repo-protocol') export PHYEXREPOprotocol=$2; toadd="$toadd $2"; shift;;
      '--onlyIfNeeded') onlyIfNeeded=true;;
      '--computeRefIfNeeded') computeRefIfNeeded=true;;
      '--prep_code-opts') check_allowed $1 $enable_prepCodeOpts 
                          prepCodeOpts=$2; toadd="$toadd $2"; shift;;
      '--perf') perffile="$(realpath $2)"; toadd="$toadd $2"; shift;;
      '--name') name=$2; toadd="$toadd $2"; shift;;
      "--no$default_buildSystem") check_allowed $1 $enable_buildSystem
                                  buildSys=$alternative_buildSystem;;
      '-a') archfile="$2"; check_allowed $1 $enable_testprogs_opt
                           toadd="$toadd $2"; shift;;
      '-A') refarchfile="$2"; check_allowed $1 $enable_testprogs_opt
                              toadd="$toadd $2"; shift;;
      '--no-perf') check_allowed $1 $enable_testprogs_opt; perf=false;;
      '--no-check') check_allowed $1 $enable_testprogs_opt; checkOpt="";;
      '-e') check_allowed $1 $enable_testprogs_opt;
            extrapolation=$2; toadd="$toadd $2"; shift;;
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

  if [ $packcreation == false -a \
       $packupdate == false -a \
       $compilation == false -a \
       $run == false -a \
       $check == false -a \
       $remove == false ]; then
    packcreation=true
    compilation=true
    run=true
    check=true
  fi

  if [ -z "${commit-}" ]; then
    echo "At least one commit hash must be provided on command line"
    exit 2
  fi
  
  if [ $check == true -a -z "${reference-}" ]; then
    echo "To perform a comparison two commit hashes are mandatory on the command line"
    exit 3
  fi
}

#######################
#### CLONE AND RUN ####
#######################

function clone_and_run {
  # If commit is hash or a tag, clone the repo and run the script using
  # the version found inside
  # Variables used:
  #   - commit
  #   - PHYEXREPOprotocol, PHYEXREPOuser
  #   - reference (if set)
  #   - commitcmd: list of arguments to use with the check_commit script
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
        # model adapted commit, we must fill the repository with
        # the last version (in history) of tools
        toolscommit=$(git log --all --full-history -- tools | head -1 | cut -d\  -f2) # deleted in this commit
        for parent in $(git rev-parse ${toolscommit}^@); do
          if [ $(git ls-tree -r $parent -- tools | wc -l) -ne 0 ]; then
            toolscommit=$parent # last commit with tools
            break
          fi
        done
        git checkout $toolscommit -- "tools"
        git checkout $toolscommit -- "pyproject.toml"
        git checkout $toolscommit -- "src/pyphyex"
      fi
  
      # Requirements
      python3 -m venv venv
      . venv/bin/activate
      python3 -m pip install -e .
  
      # Running commit
      . tools/env.sh
      if [ "$name" == "" ]; then
        mypackname="--name $commit"
      else
        mypackname=""
      fi
      $(basename $0) $commitcmd $TMP_LOC/PHYEX $reference $mypackname
    )
    exit $?
  fi
}

#############################
#### TESTS, REF AND NAME ####
#############################

function get_properties {
  testing=$(json_dictkey2value "$json_content" 'testing' '')
  ALLTests=''
  if [ "$testing" != "" ]; then
    for t in $(echo $allowedTests | sed 's/,/ /g'); do
      ref=$(json_dictkey2value "$testing" "$t" '')
      if [ ! "$ref" == "" ]; then
        ALLTests="${ALLTests},$t"
        refByTest[$t]=$ref
      fi
    done
    ALLTests="${ALLTests:1}" #Remove first character (',')
  fi
  
  if [ -z "${tests-}" ]; then
    tests=$defaultTest
  elif echo "$tests" | grep -w 'ALL' > /dev/null; then
    tests=$(echo "$tests" | sed "s:\bALL\b:$ALLTests:g")
  fi
  
  if [ ! -z "${reference-}" ]; then
    #Reference to use for each test
    for t in $(echo $tests | sed 's/,/ /g'); do
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
    done
  fi

  if [ "$name" == "" ]; then
    name=$(escape_commit $commit)
  fi
}
