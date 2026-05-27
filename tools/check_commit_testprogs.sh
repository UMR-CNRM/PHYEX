#!/bin/bash

#This script:
# - compiles the PHYEX package using a specific commit
# - runs the different test progs and checks if results are identical to a given version

#Data generation:
# - The last commit of the testprogs_data branch (based on 46t1) is able to produce the data
#   for the turb, shallow, rain_ice and ice_adjust testprogs. The code is present but must be
#   activated in the corresponding aro_* routine (as only one set of data can be produced during
#   a single execution).
# - The last commit of the testprogs_data2 branch (based on 48t3) is able to produce the data
#   for the rain_ice_old testprog.


PHYEXTOOLSDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Mutualised functions and definitions
. ${PHYEXTOOLSDIR}/check_commit_common.sh

#######################
#### CONFIGURATION ####
#######################

#About the tests:
# - allowedTests is the list of allowed tests which can depend on platform, if we ask to perform an action
#   with a test not in the allowedTests list, the action is ignored
# - defaultTest is the list of tests to perform when no '-t' option is provided on the command line.
allowedTests="ice_adjust,rain_ice,rain_ice_old,turb,shallow,lima_adjust,lima"
defaultTest=${allowedTests}

TESTDIR=${TESTPROGSDIR:=$HOME/TESTPROGS}

dirdata=$PHYEXCONF/testprogs_data
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

#Build system
default_buildSystem='fcm'
alternative_buildSystem='ecbuild'

################################
#### COMMAND LINE ARGUMENTS ####
################################

enable_testprogs_opt=true
enable_buildSystem=true
default_expand=true
command_line $@

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

if [ $useexpand == true ]; then
  expand_options=""
else
  expand_options="--noexpand"
fi

##############################
#### FUNCTION DEFINITIONS ####
##############################

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

function fill_dirdata {
  [ ! -d $dirdata ] && mkdir -p $dirdata
  CWD=$PWD
  cd $dirdata
  for file in https://github.com/UMR-CNRM/PHYEX/files/12783926/ice_adjust.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783935/rain_ice.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783942/rain_ice_old.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783945/shallow.tar.gz \
              https://github.com/UMR-CNRM/PHYEX/files/12783952/turb.tar.gz \
              https://github.com/user-attachments/files/17030876/lima_adjust.tar.gz \
              https://github.com/user-attachments/files/17330177/lima.tar.gz; do
    basefile=$(basename $file)
    if [ ! -f $(basename $basefile .tar.gz)/00000000.dat ]; then
      wget --no-check-certificate $file -O $basefile
      tar xf $basefile
      rm -f $basefile
    fi
  done
  cd $CWD
}

#######################
#### FROM A COMMIT ####
#######################

# In this case, we clone the commit and run the check_commit script
# of this commit using the cloned directory.
tools_install='--dataset'
clone_and_run

###########################
#### COMMIT PROPERTIES ####
###########################

if [ -d $commit/src ]; then
  model_ready=false
  json_content=$(scp $commit/src/common/testprogs_version.json /dev/stdout 2>/dev/null || echo "")
else
  model_ready=true
  json_content=$(scp $commit/testprogs_version.json /dev/stdout 2>/dev/null || echo "")
fi

declare -A refByTest
get_properties

#######################
#### PACK CREATION ####
#######################

[ $suppress == true -a -d $TESTDIR/$name ] && rm -rf $TESTDIR/$name
if [ $packcreation == true -a -d $TESTDIR/$name/build/with_${buildSys}/arch_${archfile} -a $onlyIfNeeded == true ]; then
  packcreation=false
fi
if [ $packcreation == true ]; then
  if [ -d $TESTDIR/$name/build/with_${buildSys}/arch_${archfile} ]; then
    echo "Directory already exists ($TESTDIR/$name/build/with_${buildSys}/arch_${archfile}),"
    echo "suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
    exit 5
  else
    echo "### Pack creation for commit $commit"

    mkdir -p $TESTDIR/$name
    cd $TESTDIR/$name/
    if [ ! -d build ]; then
      cp -r $PHYEXTOOLSDIR/../build . #We use the compilation system from the same commit as the current script
      rm -rf build/with_${buildSys}/arch_*
    else
      echo "WARNING: the compilation system is already there, we use it but it could be outdated"
    fi
  fi
fi
if [ $packupdate == true -o $packcreation == true ]; then
  if [ ! -d $TESTDIR/$name/build/with_${buildSys}/ ]; then
    echo "Compilation directory must exist ($TESTDIR/$name/build/with_${buildSys})"
    exit 9
  else
    cd $TESTDIR/$name/build/with_${buildSys}/
    packbuild=""
    [ $packupdate == true ] && packbuild='-u'
    [ $packcreation == true ] && packbuild="$packbuild -p"
    ./make_${buildSys}.sh $packbuild $expand_options --commit $commit --arch $archfile 2>&1 | tee Output_compilation_step1
  fi
fi

#####################
#### COMPILATION ####
#####################

if [ $compilation == true ]; then
  if [ $onlyIfNeeded == false -o ! -f $TESTDIR/$name/build/with_${buildSys}/arch_${archfile}/build/lib/libphyex.so ]; then
    echo "### Compilation of commit $commit"

    cd $TESTDIR/$name/build/with_${buildSys}/
    ./make_${buildSys}.sh -c --jobs=10 $expand_options --commit $commit --arch $archfile 2>&1 | tee Output_compilation_step2
  fi
fi

###################
#### EXECUTION ####
###################
if [ $run == true ]; then
  cd $TESTDIR/$name

  #Cleaning to suppress old results that may be confusing in case of a crash during the run
  if [ $onlyIfNeeded == false ]; then
    for t in $(echo $tests | sed 's/,/ /g'); do
      if [ -d tests/with_${buildSys}/arch_${archfile}/${t}${extrapolation_tag} ]; then
        rm -rf tests/with_${buildSys}/arch_${archfile}/${t}${extrapolation_tag}
      fi
    done
  fi

  #Run the tests one after the other
  firstrun=1
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then #test is allowed on this plateform
      if  [ ! -d tests/with_${buildSys}/arch_${archfile}/${t}${extrapolation_tag} ]; then #We do not enter systematically this part if onlyIfNeeded=true
        if [ $firstrun -eq 1 ]; then
          echo "### Running of commit $commit"
          firstrun=0
        fi

        for prec in "" _dp _sp; do
          if [ -f $TESTDIR/$name/build/with_${buildSys}/arch_${archfile}/build/bin/main_${t}${prec}.exe ]; then
            break
          fi
        done
        if [ ! -f $TESTDIR/$name/build/with_${buildSys}/arch_${archfile}/build/bin/main_${t}${prec}.exe ]; then
          echo "Directory does not exist ($TESTDIR/$name) or compilation has failed, please check"
          echo "Run '$0 -p -c $commit' to compile."
          exit 6
        fi

        #execution
        cd $TESTDIR/$name
        mkdir -p tests/with_${buildSys}/arch_${archfile}/${t}${extrapolation_tag}
        cd tests/with_${buildSys}/arch_${archfile}/${t}${extrapolation_tag}
        fill_dirdata
        ln -s $dirdata/$t data
        if [ $perf == true ]; then
            export DR_HOOK_OPT=prof
            export DR_HOOK=1
            export DR_HOOK_IGNORE_SIGNALS=-1
        fi
        . $TESTDIR/$name/build/with_${buildSys}/arch_${archfile}/arch.env
        set +e
        submit Output_run Stderr_run $TESTDIR/$name/build/with_${buildSys}/arch_${archfile}/build/bin/main_${t}${prec}.exe $checkOpt $extrapolation_opts
        stat=$?
        set -e
        if [ $stat -ne 0 ]; then
          cat Stderr_run
          exit $stat
        fi
        if [ $perf == true ]; then
            firstLine=$(grep -m 1 -n "^ *1" drhook.prof.0 | cut -d: -f1)
            python3 -c "import numpy, pandas
d = {'time': ('<f4', ('mean', )), 'self': ('<f4', ('mean', 'max', 'min', 'std', 'sum')),
     'total': ('<f4', ('mean', 'max', 'min', 'std', 'sum')), 'calls': ('<i4', ('sum', )),
     'self_per_call': ('<f4', ('mean', )), 'total_per_call': ('<f4', ('mean', )), 'routine': ('U256', '')}
arraynp = numpy.loadtxt('drhook.prof.0', dtype=[(k, v[0]) for (k, v) in d.items()],
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
    fi
  done
fi

#####################
#### PERFORMANCE ####
#####################

if [ $run == true -a "$perffile" != "" ]; then
  echo "### Evaluate performance for commit $commit"

  ZTD_sum=0
  ZTC_sum=0
  firstrun=1
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then
      if [ ! -f $TESTDIR/$name/build/with_${buildSys}/arch_${archfile}/build/bin/main_${t}.exe ]; then
        echo "Directory does not exist ($TESTDIR/$name) or compilation has failed, please check"
        echo "Run '$0 -p -c $commit' to compile."
        exit 7
      fi

      if [ $firstrun -eq 1 ]; then
        firstrun=0
        #Read prefered NPROMA, maximum number of points and number of openMP threads
        #for performance evaluation
        . $TESTDIR/$name/build/with_${buildSys}/arch_${archfile}/arch.env

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
        if [ $onlyIfNeeded == false ]; then
          for t2 in $(echo $tests | sed 's/,/ /g'); do
            if [ -d tests/with_${buildSys}/arch_${archfile}/${t2}${perf_extrapolation_tag} ]; then
              rm -rf tests/with_${buildSys}/arch_${archfile}/${t2}${perf_extrapolation_tag}
            fi
          done
        fi
      fi

      NPROMA=$NPROMA NBLOCKS=$NBLOCKS NTIMES=$NTIMES OMP_NUM_THREADS=${OMP_NUM_THREADS} \
              $0 -r -t $t -a ${archfile} --no-check --no-perf -e 4 --name $name ${commit}
      file=$TESTDIR/$name/tests/with_${buildSys}/arch_${archfile}/${t}${perf_extrapolation_tag}/Output_run
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

if [ $check == true ]; then
  echo "### Check commit $commit against commit $reference"

  alltests=0
  message=""
  for t in $(echo $tests | sed 's/,/ /g'); do
    if echo $allowedTests | grep -w $t > /dev/null; then
      caseref=${refByTest[$t]}
      #Conversion into directory name
      if echo $caseref | grep '/' > /dev/null; then
        refname=$(escape_commit $reference)
      else
        refname="${caseref}"
      fi
      #File comparison
      file1=$TESTDIR/$name/tests/with_${buildSys}/arch_${archfile}/${t}${extrapolation_tag}/Output_run
      file2=$TESTDIR/${refname}/tests/with_${buildSys}/arch_${refarchfile}/${t}${extrapolation_tag}/Output_run
      if [ ! -f $file2 -a $computeRefIfNeeded == true ]; then
        # The reference has not been run yet, we run it
        if [ $default_buildSystem == ${buildSys} ]; then
          buildSysArg=""
        else
          buildSysArg="--no${default_buildSystem}"
        fi
        $0 -p -c -r -t $t -a ${refarchfile} --onlyIfNeeded -e $extrapolation --no-perf ${caseref} ${buildSysArg}
      fi
      mess=""
      te=0
      if [ ! -f "$file1" ]; then
        mess="Result ($file1) for commit $commit does not exist, please run the simulation"
        te=1
      fi
      if [ ! -f "$file2" ]; then
        mess2="Result ($file2) for commit ${caseref} does not exist, please run the simulation"
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

if [ $remove == true ]; then
  echo "### Remove model directory for commit $commit"
  [ -d $TESTDIR/$name ] && rm -rf $TESTDIR/$name
fi

exit $cmpstatus
