#!/bin/bash

#set -x
set -e
set -o pipefail #abort if left command on a pipe fails

#This script:
# - compiles the PHYEX package using a specific commit
# - runs the different test progs and checks if results are identical to a given version

#ice_adjust: the ice adjust test case

#ref is commit 855b8f8 for ice_adjust, rain_ice
#ref is commit ??????? for turb
#ref is commit 7e44ab1 for shallow
#ref is commit e070d16 for rain_ice_old

#Commit e070d16 can be used for rain_ice_old (ref commit for this testprogs), and for
#turb, shallow, rain_ice and ice_adjust (as it gives the same results for these test cases).

#Some modifications have been introduced and new reference commit is 00148b1

#Data generation:
# - The last commit of the testprogs_data branch (based on 46t1) is able to produce the data
#   for the turb, shallow, rain_ice and ice_adjust testprogs. The code is present but must be
#   activated in the corresponding aro_* routine (as only one set of data can be produced during
#   a single execution).
# - The last commit of the testprogs_data2 branch (based on 48t3) is able to produce the data
#   for the rain_ice_old testprog.

#######################
#### CONFIGURATION ####
#######################

#Special pack names:
# - ref: symbolic name to the commit to use as a reference
#        useless for the commits containing a json file
specialName="ref"

#About the tests:
# - ALLTests is a list of tests to be done when '-t ALL' is used. This list is filled here
#   in case there is no ial_version.json file containig a 'testing' section. If this 'testing'
#   section exists, this list is overridden.
# - allowedTests is the list of allowed tests which can depend on platform, if we ask to perform an action
#   with a test not in the allowedTests list, the action is ignored
# - defaultTest is the list of tests to perform when no '-t' option is provided on the command line.
ALLTests="ice_adjust,rain_ice,rain_ice_old,turb,shallow,lima_adjust,lima"
defaultTest=${ALLTests}
allowedTests=${ALLTests}

separator='_' #- seprator must be in sync with prep_code.sh separator

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

TESTDIR=${TESTPROGSDIR:=$HOME/TESTPROGS}

dirdata=$PHYEXTOOLSDIR/testprogs_data
if [ $(hostname | cut -c 1-7) == 'belenos' -o $(hostname | cut -c 1-7) == 'taranis' ]; then
  defaultarchfile=MIMPIIFC1805.EPONA
  submit_method=slurm_belenos
elif [ $(hostname) == 'aurora01' ]; then
  defaultarchfile=ECMWF_NEC440MPI225SP.AU.x
  submit_method=''
else
  defaultarchfile=gnu
  submit_method=''
fi
defaultRef=ref

#Comma separated list of variables that must be set if job is executed on other node
varToExport="NPROMA,NBLOCKS,OMP_NUM_THREADS,DR_HOOK_OPT,DR_HOOK,DR_HOOK_IGNORE_SIGNALS,NVCOMPILER_ACC_GANGLIMIT"

#Options to have longer simulations, tag is used to build the directory name of the result
declare -A conf_extra_tag
declare -A conf_extra_opts
i=-1
i=$((i+1)); conf_extra_tag[$i]=""
            conf_extra_opts[$i]=""
i=$((i+1)); conf_extra_tag[$i]="_Z120_NPRO32_BLK1024"
            conf_extra_opts[$i]="--nflevg 120 --nproma 32 --blocks 1024"
i=$((i+1)); conf_extra_tag[$i]="_Z120_NPRO32_BLK256_TIMES4"
            conf_extra_opts[$i]="--nflevg 120 --nproma 32 --blocks 256 --times 4"
i=$((i+1)); conf_extra_tag[$i]="_Z120_NPRO32_BLK64_TIMES16"
            conf_extra_opts[$i]="--nflevg 120 --nproma 32 --blocks 64 --times 16"
#The following case is the one used for performance evaluation, it must remains the 4th one
i=$((i+1)); conf_extra_tag[$i]='_Z120_NPRO${NPROMA}_BLK${NBLOCKS}_TIMES${NTIMES}'
            conf_extra_opts[$i]='--nflevg 120 --nproma ${NPROMA} --blocks ${NBLOCKS} --times ${NTIMES}'


################################
#### COMMAND LINE ARGUMENTS ####
################################

function usage {
  echo "Usage: $0 [-h] [-p] [-u] [-c] [-r] [-C] [-s] [--noexpand] [-t TEST] [--repo-user USER] [--repo-protocol PROTOCOL] [-a ARCH] [-A ARCH] [--remove] [--onlyIfNeeded] [--computeRefIfNeeded] [--no-perf] [--no-check] [-e EXTRAPOLATION] [--perf FILE] commit [reference]"
  echo "commit          commit hash (or a directory, or among $specialName) to test"
  echo "reference       commit hash (or a directory, or among $specialName) REF to use as a reference"
  echo "-s              suppress compilation directory"
  echo "-p              creates pack"
  echo "-u              updates pack"
  echo "-c              performs compilation"
  echo "-r              runs the tests"
  echo "-C              checks the result against the reference"
  echo "-t TEST         comma separated list of tests to execute"
  echo "                or ALL to execute all tests"
  echo "--noexpand      do not expand mnh_expand blocks (code will be in array-syntax)"
  echo "--repo-user USER"
  echo "                user hosting the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOuser (=$PHYEXREOuser)"
  echo "--repo-protocol PROTOCOL"
  echo "                protocol (https or ssh) to reach the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOprotocol (=$PHYEXREOprotocol)"
  echo "--remove        removes the pack"
  echo "--onlyIfNeeded  do not rerun already run steps"
  echo "--computeRefIfNeeded"
  echo "                compute the reference if not already present"
  echo "--no-perf       deactivate DR_HOOK"
  echo "--no-check      suppress value printing (comparison will be impossible)"
  echo "                this option can reduce drastically the running time but only allow"
  echo "                to access performance statistics."
  echo "-a arch ARCH    architecture name to use to build and run the commit (=$defaultarchfile)"
  echo "-A arch ARCH    architecture name to use for the reference simulation (=$defaultarchfile)"
  echo "--perf FILE     add performance statistics in file FILE"
  echo "-e EXTRAPOLATION"
  echo "                extrapolate data. EXTRAPOLATION corresponds to a configuration:"
  for i in $(seq 1 $((${#conf_extra_tag[@]}-1))); do
    echo "                  - '$i': ${conf_extra_opts[$i]} (${conf_extra_tag[$i]})"
  done
  echo ""
  echo "If nothing is asked (pack creation compilation, running, check, removing) everything"
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
useexpand=""
archfile=$defaultarchfile
refarchfile=$defaultarchfile
remove=0
onlyIfNeeded=0
computeRefIfNeeded=0
perf=1
extrapolation=0
checkOpt="--check"
perffile=""

while [ -n "$1" ]; do
  case "$1" in
    '-h') usage; exit;;
    '-s') suppress=1;;
    '-p') packcreation=1;;
    '-u') packupdate=1;;
    '-c') compilation=1;;
    '-r') run=$(($run+1));;
    '-C') check=1;;
    '-t') tests="$2"; shift;;
    '--noexpand') useexpand=$1;;
    '--repo-user') export PHYEXREPOuser=$2; shift;;
    '--repo-protocol') export PHYEXREPOprotocol=$2; shift;;
    '--remove') remove=1;;
    '-a') archfile="$2"; shift;;
    '-A') refarchfile="$2"; shift;;
    '--onlyIfNeeded') onlyIfNeeded=1;;
    '--computeRefIfNeeded') computeRefIfNeeded=1;;
    '--no-perf') perf=0;;
    '--no-check') checkOpt="";;
    '--perf') perffile="$(realpath $2)"; shift;;
    '-e') extrapolation=$2; shift;;

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

if [[ ! -z "${conf_extra_tag[$extrapolation]+unset}" ]]; then
  extrapolation_tag=$(eval echo ${conf_extra_tag[$extrapolation]})
else
  echo "The extrapolation option ($extrapolation) doesn't have associated tag"
fi
if [[ ! -z "${conf_extra_opts[$extrapolation]+unset}" ]]; then
  extrapolation_opts=$(eval echo ${conf_extra_opts[$extrapolation]})
else
  echo "The extrapolation option ($extrapolation) doesn't have associated options"
fi

##############################
#### FUNCTION DEFINITIONS ####
##############################

function json_dictkey2value {
  # $1 must contain the json string
  # $2 must be the key name
  # $3 is the default value
  json_content="$1" python3 -c "import json; import os; result=json.loads(os.environ['json_content']).get('$2', '$3'); print(json.dumps(result) if isinstance(result, dict) else result)"
}

function submit {
  #usage: submit <output file> <error file> <command> [arg [arg ...]]
  output=$(realpath $1); shift
  error=$(realpath $1); shift
  if [ "$submit_method" == 'slurm_belenos' ]; then
    myscript=$TMP/testprogs$$
    if ldd $1 | grep libcuda > /dev/null; then
      #We need GPU node
      GPU="#SBATCH -p ndl"
    else
      GPU=""
    fi
    cat - << EOF > $myscript
#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 10
#SBATCH --export=$varToExport
$GPU

cd $PWD
$@
EOF
    chmod +x $myscript
    outtmp=$(mktemp)
    sbatch --wait -o $outtmp -e $error $myscript
    #Move job accounting in Stderrr (if present)
    str='#########################################'
    num=$(grep -n -m 1 $str $outtmp)
    if [ "$num" != "" ]; then
      #Acounting is present in this run
      num=$(echo $num | cut -d : -f 1)
      tail -n +$((${num}-1)) $outtmp >> $error
      head -n $((${num}-2)) $outtmp > $output
    else
      cp $outtmp $output
    fi
    rm -f $myscript $outtmp
  else
    $@ > $output 2> $error
  fi
}

###########################
#### COMMIT ADAPTATION ####
###########################

#Name and directory for compiling and executing user pack
declare -A refByTest
if echo $commit | grep '/' | grep -v '^tags/' > /dev/null; then
  #The git repository is a directory
  name=$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
  content_testprogs_version=$(scp $commit/src/testprogs/testprogs_version.json /dev/stdout 2>/dev/null || echo "")
  if [ "${content_testprogs_version}" == "" ]; then
    content_testprogs_version=$(scp $commit/src/offline/testprogs_version.json /dev/stdout 2>/dev/null || echo "")
  fi
  [ $suppress -eq 1 -a -d $TESTDIR/$name ] && rm -rf $TESTDIR/$name
elif echo $specialName | grep -w $commit > /dev/null; then
  name="$commit"
else
  #The git repository is on github
  if [[ $commit == testprogs${separator}*  || $commit == offline${separator}* ]]; then
    testprogs_version_file="testprogs_version.json"
    testprogs_version_file_alt=""
  else
    testprogs_version_file="src/offline/testprogs_version.json"
    testprogs_version_file_alt="src/testprogs/testprogs_version.json"
  fi
  if echo $commit | grep '^tags/' > /dev/null; then
    urlcommit=$(echo $commit | cut -d / -f 2-)
  else
    urlcommit=$commit
  fi
  content_testprogs_version=$(wget --no-check-certificate https://raw.githubusercontent.com/$PHYEXREPOuser/PHYEX/${urlcommit}/$testprogs_version_file -O - 2>/dev/null || echo "")
  if [ "${content_testprogs_version}" == "" -a "${testprogs_version_file_alt}" != "" ]; then
    content_testprogs_version=$(wget --no-check-certificate https://raw.githubusercontent.com/$PHYEXREPOuser/PHYEX/${urlcommit}/$testprogs_version_file_alt -O - 2>/dev/null || echo "")
  fi
  name="COMMIT$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')"
  [ $suppress -eq 1 -a -d $TESTDIR/$name ] && rm -rf $TESTDIR/$name
fi
if [ ! "${content_testprogs_version}" == "" ]; then
  testing=$(json_dictkey2value "$content_testprogs_version" 'testing' '')
  refALL=$(json_dictkey2value "$testing" "ALL" '')
  if [ ! "$testing" == "" ]; then
    ALLTests='' #We reset the list of tests
    for t in $(echo $allowedTests | sed 's/,/ /g'); do
      ref=$(json_dictkey2value "$testing" "$t" "$refALL")
      if [ ! "$ref" == "" ]; then
        ALLTests="${ALLTests},$t"
        refByTest[$t]=$ref
      fi
    done
    ALLTests="${ALLTests:1}" #Remove first character (',')
  fi
fi

#Name and directory for the reference version
if [ ! -z "${reference-}" ]; then
  declare -A refnameByTest
  #Reference to use for each test
  for t in $(echo $ALLTests | sed 's/,/ /g'); do
    #Name of the reference
    if [ "$reference" == "REF" ]; then
      if [[ ! -z "${refByTest[$t]+unset}" ]]; then #the -v test is valid only with bash > 4.3
        #The json file contained the references to use on a per test case basis
        caseref=${refByTest[$t]}
      else
        caseref=$defaultRef
      fi
      refByTest[$t]=$caseref
    else
      #The exact reference to use was given on the command line
      caseref=$reference
    fi
    refByTest[$t]=$caseref
  
    #Conversion into directory name
    if echo $caseref | grep '/' > /dev/null; then
      refname=$(echo $reference | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
    elif echo $specialName | grep -w $caseref > /dev/null; then
      refname="$caseref"
    else
      refname="COMMIT${caseref}"
    fi
    refnameByTest[$t]=$refname
  done
fi

if [ -z "${tests-}" ]; then
  tests=$defaultTest
elif echo "$tests" | grep -w 'ALL' > /dev/null; then
  tests=$(echo "$tests" | sed "s/\bALL\b/$ALLTests/g")
fi

#######################
#### PACK CREATION ####
#######################

if [ $packcreation -eq 1 -a -d $TESTDIR/$name/build/with_fcm/arch_${archfile} -a $onlyIfNeeded -eq 1 ]; then
  packcreation=0
fi
if [ $packcreation -eq 1 ]; then
  if [ -d $TESTDIR/$name/build/with_fcm/arch_${archfile} ]; then
    echo "Directory already exists ($TESTDIR/$name/build/with_fcm/arch_${archfile}),"
    echo "suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
    exit 5
  else
    echo "### Pack creation for commit $commit"

    if echo $specialName | grep -w $commit > /dev/null; then
      echo "Special commit '$commit' cannot be compiled with this script"
      exit 4
    fi

    mkdir -p $TESTDIR/$name
    cd $TESTDIR/$name/
    if [ ! -d build ]; then
      cp -r $PHYEXTOOLSDIR/../build . #We use the compilation system from the same commit as the current script
      rm -rf build/with_fcm/arch_*
    else
      echo "WARNING: the compilation system is already there, we use it but it could be outdated"
    fi
  fi
fi
if [ $packupdate -eq 1 -o $packcreation -eq 1 ]; then
  if [ ! -d $TESTDIR/$name/build/with_fcm/ ]; then
    echo "Compilation directory must exist ($TESTDIR/$name/build/with_fcm)"
    exit 9
  else
    cd $TESTDIR/$name/build/with_fcm/
    packbuild=""
    [ $packupdate -eq 1 ] && packbuild='-u'
    [ $packcreation -eq 1 ] && packbuild="$packbuild -p"
    ./make_fcm.sh $packbuild $useexpand --commit $commit --arch $archfile 2>&1 | tee Output_compilation_step1
  fi
fi

#####################
#### COMPILATION ####
#####################

if [ $compilation -eq 1 ]; then
  if [ $onlyIfNeeded -eq 0 -o ! -f $TESTDIR/$name/build/with_fcm/arch_${archfile}/build/bin/libphyex.so ]; then
    echo "### Compilation of commit $commit"

    cd $TESTDIR/$name/build/with_fcm/
    ./make_fcm.sh -c --jobs=10 $useexpand --commit $commit --arch $archfile 2>&1 | tee Output_compilation_step2
  fi
fi

###################
#### EXECUTION ####
###################
if [ $run -ge 1 ]; then
  cd $TESTDIR/$name

  #Cleaning to suppress old results that may be confusing in case of a crash during the run
  if [ $onlyIfNeeded -eq 0 ]; then
    for t in $(echo $tests | sed 's/,/ /g'); do
      if [ -d tests/with_fcm/arch_${archfile}/${t}${extrapolation_tag} ]; then
        rm -rf tests/with_fcm/arch_${archfile}/${t}${extrapolation_tag}
      fi
    done
  fi

  #Run the tests one after the other
  firstrun=1
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then #test is allowed on this plateform
      if  [ ! -d tests/with_fcm/arch_${archfile}/${t}${extrapolation_tag} ]; then #We do not enter systematically this part if onlyIfNeeded=1
        if [ $firstrun -eq 1 ]; then
          echo "### Running of commit $commit"
          firstrun=0
        fi

        if [ ! -f $TESTDIR/$name/build/with_fcm/arch_${archfile}/build/bin/main_${t}.exe ]; then
          echo "Directory does not exist ($TESTDIR/$name) or compilation has failed, please check"
          echo "Run '$0 -p -c $commit' to compile."
          exit 6
        fi

        #execution
        cd $TESTDIR/$name
        mkdir -p tests/with_fcm/arch_${archfile}/${t}${extrapolation_tag}
        cd tests/with_fcm/arch_${archfile}/${t}${extrapolation_tag}
        ln -s $dirdata/$t data
        if [ $perf -eq 1 ]; then
            export DR_HOOK_OPT=prof
            export DR_HOOK=1
            export DR_HOOK_IGNORE_SIGNALS=-1
        fi
        . $TESTDIR/$name/build/with_fcm/arch_${archfile}/arch.env
        set +e
        submit Output_run Stderr_run $TESTDIR/$name/build/with_fcm/arch_${archfile}/build/bin/main_${t}.exe $checkOpt $extrapolation_opts
        stat=$?
        set -e
        if [ $stat -ne 0 ]; then
          cat Stderr_run
          exit $stat
        fi
        if [ $perf -eq 1 ]; then
            firstLine=$(grep -m 1 -n "^ *1" drhook.prof.0 | cut -d: -f1)
            python3 -c "import numpy, pandas
d = {'time': ('<f4', ('mean', )), 'self': ('<f4', ('mean', 'max', 'min', 'std', 'sum')),
     'total': ('<f4', ('mean', 'max', 'min', 'std', 'sum')), 'calls': ('<i4', ('sum', )),
     'self_per_call': ('<f4', ('mean', )), 'total_per_call': ('<f4', ('mean', )), 'routine': ('U256', '')}
arraynp = numpy.loadtxt('drhook.prof.0', dtype=[(k, v[0]) for (k, v) in d.items()],
                        converters={8: lambda s: s.split(b'@')[0].lstrip(b'*')},
                        skiprows=$firstLine - 1, usecols=[1, 3, 4, 5, 6, 7, 8])
df = pandas.DataFrame(arraynp).groupby('routine').agg(
      **{k + '_' + agg:pandas.NamedAgg(column=k, aggfunc=agg)
         for (k, agg) in [(k, agg) for k in d.keys() for agg in d[k][1]]
         if k != 'routine'}).sort_values('self_sum', ascending=False)
df.index.name += ' ordered by self_sum'
with open('drhook.prof.agg', 'w') as f: f.write(df.to_string())
"
        fi
      fi
    fi
  done
fi

#####################
#### PERFORMANCE ####
#####################

if [ $run -ge 1 -a "$perffile" != "" ]; then
  echo "### Evaluate performance for commit $commit"

  ZTD_sum=0
  ZTC_sum=0
  firstrun=1
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then
      if [ ! -f $TESTDIR/$name/build/with_fcm/arch_${archfile}/build/bin/main_${t}.exe ]; then
        echo "Directory does not exist ($TESTDIR/$name) or compilation has failed, please check"
        echo "Run '$0 -p -c $commit' to compile."
        exit 7
      fi

      if [ $firstrun -eq 1 ]; then
        firstrun=0
        #Read prefered NPROMA, maximum number of points and number of openMP threads
        #for performance evaluation
        . $TESTDIR/$name/build/with_fcm/arch_${archfile}/arch.env

        #Experiement size 
        NPOINTS=100000
        if [ ${NPOINTS_perf-${NPOINTS}} -lt ${NPOINTS} ]; then
          NTIMES=$(python3 -c "print(round(${NPOINTS}/${NPOINTS_perf}))")
          NPOINTS=${NPOINTS_perf}
        else
          NTIMES=1
        fi
        NPROMA=${NPROMA_perf-32}
        OMP_NUM_THREADS=${OMP_NUM_THREADS_perf-8}
        NBLOCKS=$(($NPOINTS/$NPROMA/8*8)) #must be divisible by 8
        perf_extrapolation_tag=$(NPROMA=$NPROMA; NBLOCKS=$NBLOCKS; NTIMES=$NTIMES; eval echo ${conf_extra_tag[4]})

        #Cleaning to suppress old results that may be confusing in case of a crash during the run
        if [ $onlyIfNeeded -eq 0 ]; then
          for t2 in $(echo $tests | sed 's/,/ /g'); do
            if [ -d tests/with_fcm/arch_${archfile}/${t2}${perf_extrapolation_tag} ]; then
              rm -rf tests/with_fcm/arch_${archfile}/${t2}${perf_extrapolation_tag}
            fi
          done
        fi
      fi

      NPROMA=$NPROMA NBLOCKS=$NBLOCKS NTIMES=$NTIMES OMP_NUM_THREADS=${OMP_NUM_THREADS} $0 -r -t $t -a ${archfile} --no-check --no-perf -e 4 ${commit}
      file=$TESTDIR/$name/tests/with_fcm/arch_${archfile}/${t}${perf_extrapolation_tag}/Output_run
      if [ -f $file ]; then
        ZTD=$(grep -m 1 "ZTD =" $file | awk '{print $4}')
        if [ "$ZTD" != "" ]; then
          ZTD_sum=$(python3 -c "print(${ZTD_sum} if ${ZTD_sum} < 0. else (${ZTD_sum} + ${ZTD}))")
        else
          ZTD=-999
          ZTD_sum=-999
        fi
        ZTC=$(grep -m 1 "ZTC =" $file | awk '{print $4}')
        if [ "$ZTC" != "" ]; then
          ZTC_sum=$(python3 -c "print(${ZTC_sum} if ${ZTC_sum} < 0. else (${ZTC_sum} + ${ZTC}))")
        else
          ZTC=-999
          ZTC_sum=-999
        fi
      else
        ZTD=-999
        ZTD_sum=-999
        ZTC=-999
        ZTC_sum=-999
      fi
      echo "$commit testprogs $t $ZTD $ZTC" >> "$perffile"
    fi
  done
  echo "$commit testprogs ALL $ZTD_sum $ZTC_sum" >> "$perffile"
fi

####################
#### COMPARISON ####
####################

if [ $check -eq 1 ]; then
  echo "### Check commit $commit against commit $reference"

  alltests=0
  message=""
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then
      #Run the reference if needed
      if [ $computeRefIfNeeded -eq 1 ]; then
        $0 -p -c -r -t $t -a ${refarchfile} --onlyIfNeeded -e $extrapolation --no-perf ${refByTest[$t]}
      fi

      #File comparison
      file1=$TESTDIR/$name/tests/with_fcm/arch_${archfile}/${t}${extrapolation_tag}/Output_run
      file2=$TESTDIR/${refnameByTest[$t]}/tests/with_fcm/arch_${refarchfile}/${t}${extrapolation_tag}/Output_run
      mess=""
      te=0
      if [ ! -f "$file1" ]; then
        mess="Result ($file1) for commit $commit does not exist, please run the simulation"
        te=1
      fi
      if [ ! -f "$file2" ]; then
        mess2="Result ($file2) for commit ${refByTest[$t]} does not exist, please run the simulation"
        te=1
        if [ "$mess" = "" ]; then
          mess=$mess2
        else
          mess="$mess and $mess2"
        fi
      fi
      if [ $te -eq 0 ]; then
        set +e
        mess=$($PHYEXTOOLSDIR/compare.py --testprogs $file1 $file2)
        te=$?
        set -e
      fi
      [ $te -ne 0 ] && message="$message $mess \n"
      alltests=$(($alltests+$te))
    fi
  done
  if [ $alltests -eq 0 ]; then
    echo "SUCCESS, files are identical"
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
  [ -d $TESTDIR/$name ] && rm -rf $TESTDIR/$name
fi

exit $cmpstatus
