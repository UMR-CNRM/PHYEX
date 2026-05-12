#!/bin/bash

#set -x
set -e
set -o pipefail #abort if left command on a pipe fails

#This script:
# - compiles the AROME model using a specific commit for the externalised physics
# - runs a small 3D case and checks if results are identical to a given version

#small_3D_np2: on only 2 procs
#small_3D_alt1: options around time-step dependency, CFRAC_ICE_*='T', CSEDIM='SPLI', LSEDIM_AFTER=.T.
#small_3D_alt2: CCLOUD='OLD3'
#small_3D_alt3: PRFR
#small_3D_alt4: small_3D_alt1 + CSNOWRIMING='OLD'
#small_3D_alt5: CCLOUD='ICE4'
#small_3D_alt6: CMF_UPDRAFT='RAHA', CMF_CLOUD='BIGA'
#small_3D_alt7: CMF_CLOUD='STAT', LOSIGMAS=.FALSE.
#small_3D_alt8: CMF_UPDRAFT='RHCJ'
#small_3D_alt9: CCLOUD='OLD3', OCND2=.T.
#small_3D_alt10: LCRIAUTI=F, not included in the list because it is not sufficiently different from other tests
#small_3D_lima: LIMA scheme
#small_3D_alt11: same as small_3D but with a different value for NPROMICRO (must give exactly the same results)
#small_3D_alt12: same as small_3D but with LPACK_MICRO=.F. (must give exactly the same results)
#small_3D_xfrmin: same as small_3D_alt2 but with specified values for XFRMIN(16:17)
#arp_t31: Ph. Marguinaud's ARPEGE toy

#######################
#### CONFIGURATION ####
#######################

cycle='50t2' # Used to buil pack names

#About the tests:
# - allowedTests is the list of allowed tests which can depend on platform, if we ask to perform an action
#   with a test not in the allowedTests list, the action is ignored
# - defaultTest is the list of tests to perform when no '-t' option is provided on the command line.
defaultTest="small_3D"
allowedTests="small_3D,small_3D_np2,small_3D_alt1,small_3D_alt2,small_3D_alt3,small_3D_alt4,small_3D_alt5,small_3D_alt6,small_3D_alt7,small_3D_alt8,small_3D_alt9,small_3D_alt10,small_3D_alt11,small_3D_alt12,small_3D_lima,small_3D_xfrmin,arp_t31"

separator='_' #- be carrefull, gmkpack (at least on belenos) has multiple allergies (':', '.', '@')
              #- seprator must be in sync with prep_code.sh separator

PHYEXTOOLSDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

HOMEPACK=${HOMEPACK:=$HOME/pack}

dirconf=$PHYEXTOOLSDIR/conf_tests
if [ $(hostname | cut -c 1-7) == 'belenos' -o $(hostname | cut -c 1-7) == 'taranis' ]; then
  HPC=1
  gmkpack_l=PHYEX50T2ifort
  gmkpack_o=sp
  allowedTests="${allowedTests},big_3D"
else
  HPC=0
  gmkpack_l=PHYEX50T2gfort
  gmkpack_o=dp
fi

################################
#### COMMAND LINE ARGUMENTS ####
################################

function usage {
  echo "Usage: $0 [-h] [-p] [-u] [-c] [-r] [-C] [-s] [-f] [--noexpand] [-t TEST] [--repo-user USER] [--repo-protocol PROTOCOL] [--remove] [--onlyIfNeeded] [--computeRefIfNeeded] [--prep_code-opts 'OPTS'] [--perf FILE] [--name NAME] commit [reference]"
  echo "commit          commit hash (or a directory) to test"
  echo "reference       commit hash (or a directory) REF to use as a reference"
  echo "-s              suppress compilation pack"
  echo "-p              creates pack"
  echo "-u              updates pack (experimental, only for 'small' updates where"
  echo "                              the list of files does not change)"
  echo "-c              performs compilation"
  echo "-r              runs the tests"
  echo "-C              checks the result against the reference"
  echo "-t TEST         comma separated list of tests to execute"
  echo "                or ALL to execute all tests"
  echo "--noexpand      do not expand mnh_expand blocks (code will be in array-syntax)"
  echo "-f              full compilation (do not use pre-compiled pack)"
  echo "--repo-user USER"
  echo "                user hosting the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOuser (=$PHYEXREOuser)"
  echo "--repo-protocol PROTOCOL"
  echo "                protocol (https or ssh) to reach the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOprotocol (=$PHYEXREOprotocol)"
  echo "--remove        removes the pack"
  echo "--onlyIfNeeded  performs the pack creation and/or the compilation and/or the execution"
  echo "                only if the step has not already been done"
  echo "--computeRefIfNeeded"
  echo "                computes the missing references"
  echo "--prep_code-opts 'OPTS'"
  echo "                OPTS is added to the call to prep_code (e.g. --prep_code_opts '--lowerCase'"
  echo "                to transfor all source codes in lower case). Help on prep_code.sh options"
  echo "                can be found with 'prep_code.sh -h'. Note: don't forget to enclose OPTS in ' or \""
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
  echo
  echo "The -f flag (full recompilation) is active only at pack creation"
  echo
  echo "The PHYEXROOTPACK environment variable, if set, is used as the argument"
  echo "of the --rootpack option of ial-git2pack/ial-to_pack, for incremental packs."
}

packcreation=0
packupdate=0
compilation=0
run=0
check=0
commit=""
reference=""
tests=""
suppress=0
useexpand=1
fullcompilation=0
remove=0
onlyIfNeeded=0
computeRefIfNeeded=0
prepCodeOpts=""
perffile=""
packBranch=""

commitcmd="" # To execute again check_commit with the same options (except commit and ref)
while [ -n "$1" ]; do
  toadd="$1"
  case "$1" in
    '-h') usage; exit;;
    '-s') suppress=1;;
    '-p') packcreation=1;;
    '-u') packupdate=1;;
    '-c') compilation=1;;
    '-r') run=1;;
    '-C') check=1;;
    '-t') tests="$2"; toadd="$toadd $2"; shift;;
    '--noexpand') useexpand=0;;
    '-f') fullcompilation=1;;
    '--repo-user') export PHYEXREPOuser=$2; toadd="$toadd $2"; shift;;
    '--repo-protocol') export PHYEXREPOprotocol=$2; toadd="$toadd $2"; shift;;
    '--remove') remove=1;;
    '--onlyIfNeeded') onlyIfNeeded=1;;
    '--computeRefIfNeeded') computeRefIfNeeded=1;;
    '--prep_code-opts') prepCodeOpts=$2; toadd="$toadd $2"; shift;;
    '--perf') perffile="$(realpath $2)"; toadd="$toadd $2"; shift;;
    '--name') packBranch=$2; toadd="$toadd $2"; shift;;
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
     $packupdate -eq 0 -a \
     $compilation -eq 0 -a \
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

function exescript () {
  #usage: exescript <output file> <script> [arg [arg ...]]
  output=$1
  shift
  if [ $HPC -eq 1 ]; then
    sbatch --wait -o $output $@
    cat $output
  else
    $@ 2>&1 | tee $output
  fi
}

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
      # AROME adapted commit, we must fill the repository with
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
      mypackname="--name COMMIT$commit"
    else
      mypackname=""
    fi
    check_commit_ial.sh $commitcmd $TMP_LOC/PHYEX $reference $mypackname
  )
  exit $?
fi

###########################
#### COMMIT PROPERTIES ####
###########################

if [ -d $commit/src ]; then
  arome_ready=false
  content_ial_version=$(scp $commit/src/arome/ial_version.json /dev/stdout 2>/dev/null || echo "")
else
  arome_ready=true
  content_ial_version=$(scp $commit/ial_version.json /dev/stdout 2>/dev/null || echo "")
fi

declare -A refByTest
testing=$(json_dictkey2value "$content_ial_version" 'testing' '')
ALLTests=''
for t in $(echo $allowedTests | sed 's/,/ /g'); do
  ref=$(json_dictkey2value "$testing" "$t" '')
  if [ ! "$ref" == "" ]; then
    ALLTests="${ALLTests},$t"
    refByTest[$t]=$ref
  fi
done
ALLTests="${ALLTests:1}" #Remove first character (',')

if [ -z "${tests-}" ]; then
  tests=$defaultTest
elif echo "$tests" | grep -w 'ALL' > /dev/null; then
  tests=$(echo "$tests" | sed "s/\bALL\b/$ALLTests/g")
fi

#Name is choosen such as it can be produced with a main pack: PHYEX/${cycle}_XXXXXXXXX.01.${gmkpack_l}.${gmkpack_o}
if [ "$packBranch" == "" ]; then
  packBranch=$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
fi
name="PHYEX/${cycle}_${packBranch}.01.${gmkpack_l}.${gmkpack_o}"

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
      refname="PHYEX/*_COMMIT${caseref}.01.${gmkpack_l}.${gmkpack_o}"
    fi
    refnameByTest[$t]=$refname
  done
fi

#######################
#### PACK CREATION ####
#######################

[ $suppress -eq 1 -a -d $HOMEPACK/$name ] && rm -rf $HOMEPACK/$name
if [ $packcreation -eq 1 -a -d $HOMEPACK/$name -a $onlyIfNeeded -eq 1 ]; then
  packcreation=0
fi
if [ $packcreation -eq 1 ]; then
  if [ -d $HOMEPACK/$name ]; then
    echo "Pack already exists ($HOMEPACK/$name), suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
    exit 5
  else
    echo "### Pack creation for commit $commit"

    export GMKTMP=/dev/shm

    ialcmd=ial-to_pack
    if ! which $ialcmd > /dev/null 2>&1; then
      echo "$ialcmd not found, please install it"
      exit 7
    fi
    #Pack creation using ial-git2pack
    tmpbuilddir=$(mktemp -d)
    trap "\rm -rf $tmpbuilddir" EXIT
    cd $tmpbuilddir
    git clone $(json_dictkey2value "$content_ial_version" 'IALrepo' 'git@github.com:ACCORD-NWP/IAL.git')
    cd IAL
    git checkout $(json_dictkey2value "$content_ial_version" 'IALcommit' 'XXXXXXX')
    IALbundle_tag=$(json_dictkey2value "$content_ial_version" 'IALbundle_tag' '')
    [ "$IALbundle_tag" != "" ] && IALbundle_tag="--hub_bundle_tag $IALbundle_tag"
    ROOTPACKopt=""
    if [ $fullcompilation == 0 ]; then
      kind=incr
      if [ "$PHYEXROOTPACK" != "" ]; then
        ROOTPACKopt="--rootpack $PHYEXROOTPACK"
      fi
    else
      kind=main
    fi
    echo y | $ialcmd -l ${gmkpack_l} -o ${gmkpack_o} -t $kind -n 10 \
                     -r $tmpbuilddir/IAL -p masterodb \
                     $ROOTPACKopt \
                     --homepack $tmpbuilddir/pack $IALbundle_tag

    #Moving
    oldname=$(echo $tmpbuilddir/pack/*)
    mv $oldname $HOMEPACK/$name
    for file in $HOMEPACK/$name/ics_*; do
      sed -i "s*$oldname*$HOMEPACK/$name*g" $file
    done
    rm -rf $tmpbuilddir

    cd $HOMEPACK/$name
    if [ "$kind" == 'main' ]; then
      #Workarounds for 50t2 compilation: update falfilfa to relax eccodes version check
      #Only needed on SIRES computer, latest eccode version is available on HPC
      cd hub/local/src/FALFILFA/falfilfa
      git cherry-pick 15359c1
      cd $HOMEPACK/$name
    fi

    #Prepare PHYEX inclusion
    rm -rf hub/local/src/PHYEX/phyex
    mkdir -p hub/local/src/PHYEX/phyex
  fi
fi
if [ $packupdate -eq 1 -o $packcreation -eq 1 ]; then
  cd $HOMEPACK/$name/hub/local/src/PHYEX/phyex

  if [ $useexpand == 1 ]; then
    expand_options="--mnhExpand"
  else
    expand_options=""
  fi
  prep_code=$PHYEXTOOLSDIR/prep_code.sh
  echo "Copy $commit"
  mkdir PHYEX
  if [ $arome_ready == false ]; then
    scp -q -r $commit/src PHYEX/
    subs="-s gmkpack_ignored_files -s turb -s micro -s aux -s conv  -s CMakeLists.txt -s cmake"
    $prep_code $prepCodeOpts $subs -m arome PHYEX -- --shumanFUNCtoCALL --removeACC $expand_options
  else
    scp -q -r $commit/* PHYEX/
    $prep_code $prepCodeOpts PHYEX #Ready for inclusion
  fi
  find PHYEX -type f -exec touch {} \; #to be sure a recompilation occurs
  if [ $packupdate -eq 1 ]; then
    #Update only modified files
    cd PHYEX
    for file in $(find turb micro conv aux cmake -type f); do
      mvdiff $file ../$file
    done
    file=CMakeLists.txt; [ -f $file ] && mvdiff $file ../$file
    cd ..
    rm -rf PHYEX
  else
    #Move PHYEX source files
    for rep in turb micro conv aux cmake; do
      [ -d PHYEX/$rep ] && mv PHYEX/$rep .
    done
    file=PHYEX/CMakeLists.txt; [ -f $file ] && mv $file .
  fi
  if [ -f PHYEX/gmkpack_ignored_files ]; then
    #gmkpack_ignored_files contains a list of file, present in the reference pack, that is not used anymore
    #and must be excluded from compilation (in case of a full comilation) or from re-compilation (in case of a non-full
    #compilation).
    if [ $fullcompilation == 0 ]; then
      #Content is added in the ics_masterodb script
      if [ $packupdate -eq 0 ]; then
        if [ $(cat PHYEX/gmkpack_ignored_files | wc -l) != 0 ]; then
          sed -i "/^end_of_ignored_files/i $(first=1; for line in $(cat PHYEX/gmkpack_ignored_files); do echo -n $(test $first -ne 1 && echo \\n)${line}; first=0; done)" $HOMEPACK/$name/ics_masterodb
        fi
      fi
    else
      #Files must be suppressed (non phyex files)
      for file in $(cat PHYEX/gmkpack_ignored_files); do
        [ -f $HOMEPACK/$name/src/local/$file ] && rm -f $HOMEPACK/$name/src/local/$file
      done
      if [ -d $HOMEPACK/$name/src/local/mpa/dummy ]; then
        [ ! "$(ls -A $HOMEPACK/$name/src/local/mpa/dummy)" ] && rmdir $HOMEPACK/$name/src/local/mpa/dummy
      fi
    fi
  fi

  rm -rf PHYEX

  if [ -d $HOMEPACK/$name/hub/local/src/PHYEX/phyex -a -d $HOMEPACK/$name/src/main ]; then
    # Incremental pack with updated Hub, we must put, in src/local, the files using
    # PHYEX to force a re-compilation
    cd $HOMEPACK/$name/src

    for file in main/mpa/conv/externals/aro_conv_mnh.F90 \
                main/arpifs/phys_dmn/suphmpa.F90 \
                $(grep -l -i -e phyex_t -e modd_nsv $(find main/mpa main/arpifs -name \*.F90 -o -name \*.h)); do
      localfile=$(echo $file | sed 's;^main/;local/;')
      if [ ! -f $localfile ]; then
        cp $file $localfile
      fi
    done
  fi
fi

#####################
#### COMPILATION ####
#####################

if [ $compilation -eq 1 ]; then
  if [ $onlyIfNeeded -eq 0 -o ! -f $HOMEPACK/$name/bin/MASTERODB ]; then
    echo "### Compilation of commit $commit"

    cd $HOMEPACK/$name
    sed -i 's/GMK_THREADS=1$/GMK_THREADS=10/' ics_masterodb
  
    [ -f ics_packages ] && exescript Output_compilation_hub ics_packages
    exescript Output_compilation ics_masterodb
    if [ -f bin/MASTERODB \
         -a $(grep Error Output_compilation | \
              grep -v ErrorCovariance3D | \
              grep -v TestErrorHandler | \
              grep -v instantiateObsErrorFactory | \
              grep -v "'Error" | \
              grep -v "'CPLNG: Error" | \
              grep -v '"Error' | \
              grep -v "'*** Error" | \
              grep -v "\-\- Up-to-date:" | wc -l) -ne 0 ]; then
      echo "check_commit_ial: MASTERODB was produced but errors occured during compilation:"
      grep Error Output_compilation | \
            grep -v ErrorCovariance3D | \
            grep -v TestErrorHandler | \
            grep -v instantiateObsErrorFactory | \
            grep -v "'Error" | \
            grep -v "'CPLNG: Error" | \
            grep -v '"Error' | \
            grep -v "'*** Error" | \
            grep -v "\-\- Up-to-date:"
      echo "check_commit_ial: MASTERODB suppressed!"
      rm -f bin/MASTERODB
      exit 12
    fi
  fi
fi

###################
#### EXECUTION ####
###################

if [ $run -eq 1 ]; then
  #Cleaning to suppress old results that may be confusing in case of a crash during the run
  if [ $onlyIfNeeded -eq 0 ]; then
    for t in $(echo $tests | sed 's/,/ /g'); do
      cd $HOMEPACK/$name
      if [ -d conf_tests/$t ]; then
        rm -rf conf_tests/$t
      fi
    done
  fi

  #Run the tests one after the other
  firstrun=1
  for t in $(echo $tests | sed 's/,/ /g'); do #loop on tests
    if echo $allowedTests | grep -w $t > /dev/null; then #test is allowed on this plateform
      cd $HOMEPACK/$name
      if [ ! -d conf_tests/$t ]; then #We do not enter systematically this part if onlyIfNeeded=1
        if [ $firstrun -eq 1 ]; then
          echo "### Running of commit $commit"
          firstrun=0
        fi

        if [ ! -f $HOMEPACK/$name/bin/MASTERODB ]; then
          echo "Pack does not exist ($HOMEPACK/$name) or compilation has failed, please check"
          exit 6
        fi

        mkdir -p conf_tests/$t
        cd conf_tests/$t
        t1=$(($(date +%s%N)/1000)) #current time in milliseconds
        MYLIB=$name TESTDIR=$dirconf/$t exescript Output_run $dirconf/$t/ar?.sh
        t2=$(($(date +%s%N)/1000))
        if [ "$perffile" != "" ]; then
          #The elapsed time is not relevant when the model runs with a queuing system (HPC)
          echo "$commit ial $t $(($t2-$t1))" >> "$perffile"
        fi

        #Profiling
        if ls drhook.prof.* > /dev/null 2>&1; then
          firstfile=1
          for file in drhook.prof.*; do
            if [ $firstfile -eq 1 ]; then
              cp $file drhook.prof.concat
              firstfile=0
            else
              #append only relevant lines
              grep -e '^ *[0-9]' $file >> drhook.prof.concat
            fi
          done
          firstLine=$(grep -m 1 -n "^ *1" drhook.prof.concat | cut -d: -f1)
          python3 -c "import numpy, pandas
d = {'time': ('<f4', ('mean', )), 'self': ('<f4', ('mean', 'max', 'min', 'std', 'sum')),
     'total': ('<f4', ('mean', 'max', 'min', 'std', 'sum')), 'calls': ('<i4', ('sum', )),
     'self_per_call': ('<f4', ('mean', )), 'total_per_call': ('<f4', ('mean', )), 'routine': ('U256', '')}
arraynp = numpy.loadtxt('drhook.prof.concat', dtype=[(k, v[0]) for (k, v) in d.items()],
                        converters={8: lambda s: s.split(b'@')[0].lstrip(b'*')},
                        skiprows=$firstLine - 1, usecols=[1, 3, 4, 5, 6, 7, 8], encoding='bytes')
df = pandas.DataFrame(arraynp).groupby('routine').agg(
      **{k + '_' + agg:pandas.NamedAgg(column=k, aggfunc=agg)
         for (k, agg) in [(k, agg) for k in d.keys() for agg in d[k][1]]
         if k != 'routine'}).sort_values('self_sum', ascending=False)
df.index.name += ' ordered by self_sum'
with open('drhook.prof.agg', 'w') as f: f.write(df.to_string())
"
        fi
      fi
    else
      echo "The test $t is not allowed"
    fi
  done
fi

####################
#### COMPARISON ####
####################
if [ $check -eq 1 ]; then
  echo "### Check commit $commit against commit $reference"

  allt=0
  message=""
  filestocheck=""
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then
      #Files to compare
      if echo $t | grep 'small' > /dev/null; then
        filestocheck="$filestocheck ${t},conf_tests/$t/ICMSHFPOS+0002:00"
        #filestocheck="$filestocheck ${t},conf_tests/$t/ICMSHFPOS+0002:00 ${t},conf_tests/$t/DHFDLFPOS+0002"
      else
        filestocheck="$filestocheck ${t},conf_tests/$t/NODE.001_01"
      fi
    else
      echo "The test $t is not allowed"
    fi
  done

  for tag_file in $filestocheck; do
      tag=$(echo $tag_file | cut -d, -f1)
      file=$(echo $tag_file | cut -d, -f2)
      refname=${refnameByTest[$tag]}
      ref=${refByTest[$tag]}
      file1=$HOMEPACK/$name/$file
      file2=$(echo $HOMEPACK/$refname/$file) #echo to enable shell substitution

      if [ ! -f $file2 -a $computeRefIfNeeded -eq 1 ]; then
        # The reference has not been run yet, we run it
        $0 -p -c -r -t $t --onlyIfNeeded ${ref}
        file2=$(echo $HOMEPACK/$refname/$file) # reperform shell substitution
      fi

      mess=""
      t=0
      if [ ! -f "$file1" ]; then
        mess="Result ($file1) for commit $commit does not exist, please run the simulation"
        t=1
      fi
      if [ ! -f "$file2" ]; then
        mess2="Reference result ($file2) for commit $ref does not exist, please run the simulation"
        t=1
        if [ "$mess" = "" ]; then
          mess=$mess2
        else
          mess="$mess and $mess2"
        fi
      fi
      if [ $t -eq 0 ]; then
        cmd="$PHYEXTOOLSDIR/compare.py"
        if [ ! -x $cmd ]; then
          echo "Command not found: \"$cmd\""
          exit 10
        fi
        if [ $(basename $file) == ICMSHFPOS+0002:00 ]; then
          #historic files
          cmd="$cmd --binary $file1 $file2 256"
        elif [ $(basename $file) == DHFDLFPOS+0002 ]; then
          #DDH files
          ddh_images="$HOMEPACK/$name/ddh_diff_${tag}.png"
          if [ `hostname` == 'sxphynh' ]; then
            [ ! -d /d0/images/$USER ] && mkdir /d0/images/$USER
            ddh_images="$ddh_images /d0/images/$USER/ddh_diff_${tag}.png"
          fi
          cmd="$cmd --ddh $file1 $file2 --ddhplots $ddh_images"
        elif [ $(basename $file) == NODE.001_01 ]; then
          #Output listing
          cmd="$cmd --node $file1 $file2"
        else
          cmd="$cmd --binary $file1 $file2 0"
        fi
        set +e
        mess=$($cmd)
        t=$?
        set -e
      fi
      [ $t -ne 0 ] && message="$message $file : $mess \n"
      allt=$(($allt+$t))
  done
  if [ $allt -eq 0 ]; then
    echo "SUCCESS, files are (nearly) identical"
  else
    echo "*************** Files are different *******************"
    echo -e "$message"
    cmpstatus=50
  fi
fi

##################
#### CLEANING ####
##################

if [ $remove -eq 1 ]; then
  echo "### Remove model directory for commit $commit"
  [ -d $HOMEPACK/$name ] && rm -rf $HOMEPACK/$name
fi

exit $cmpstatus
