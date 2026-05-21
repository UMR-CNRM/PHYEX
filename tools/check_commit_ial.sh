#!/bin/bash

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
defaultTest="small_3D"
allowedTests="small_3D,small_3D_np2,small_3D_alt1,small_3D_alt2,small_3D_alt3,small_3D_alt4,small_3D_alt5,small_3D_alt6,small_3D_alt7,small_3D_alt8,small_3D_alt9,small_3D_alt10,small_3D_alt11,small_3D_alt12,small_3D_lima,small_3D_xfrmin,arp_t31"

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

default_expand=true
extra_doc="The PHYEXROOTPACK environment variable, if set, is used as the argument
of the --rootpack option of ial-git2pack/ial-to_pack, for incremental packs."
enable_prepCodeOpts=true
command_line $@

##############################
#### FUNCTION DEFINITIONS ####
##############################

function submit () {
  #usage: submit <output file> <script> [arg [arg ...]]
  output=$1
  shift
  if [ $HPC -eq 1 ]; then
    sbatch --wait -o $output $@
    cat $output
  else
    $@ 2>&1 | tee $output
  fi
}

#######################
#### FROM A COMMIT ####
#######################

# In this case, we clone the commit and run the check_commit script
# of this commit using the cloned directory.
tools_install=''
clone_and_run

###########################
#### COMMIT PROPERTIES ####
###########################

if [ -d $commit/src ]; then
  model_ready=false
  json_content=$(scp $commit/src/arome/ial_version.json /dev/stdout 2>/dev/null || echo "")
else
  model_ready=true
  json_content=$(scp $commit/ial_version.json /dev/stdout 2>/dev/null || echo "")
fi

declare -A refByTest
get_properties

cycle=$(json_dictkey2value "$json_content" 'cycle' '')
ialdir="PHYEX/${cycle}_${name}.01.${gmkpack_l}.${gmkpack_o}"

#######################
#### PACK CREATION ####
#######################

[ $suppress == true -a -d $HOMEPACK/$ialdir ] && rm -rf $HOMEPACK/$ialdir
if [ $packcreation == true -a -d $HOMEPACK/$ialdir -a $onlyIfNeeded == true ]; then
  packcreation=false
fi
if [ $packcreation == true ]; then
  if [ -d $HOMEPACK/$ialdir ]; then
    echo "Pack already exists ($HOMEPACK/$ialdir),"
    echo "suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
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
    git clone $(json_dictkey2value "$json_content" 'IALrepo' 'git@github.com:ACCORD-NWP/IAL.git')
    cd IAL
    git checkout $(json_dictkey2value "$json_content" 'IALcommit' 'XXXXXXX')
    IALbundle_tag=$(json_dictkey2value "$json_content" 'IALbundle_tag' '')
    [ "$IALbundle_tag" != "" ] && IALbundle_tag="--hub_bundle_tag $IALbundle_tag"
    ROOTPACKopt=""
    if [ $fullcompilation == false ]; then
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
    mv $oldname $HOMEPACK/$ialdir
    for file in $HOMEPACK/$ialdir/ics_*; do
      sed -i "s*$oldname*$HOMEPACK/$ialdir*g" $file
    done
    rm -rf $tmpbuilddir

    cd $HOMEPACK/$ialdir
    if [ "$kind" == 'main' ]; then
      #Workarounds for 50t2 compilation: update falfilfa to relax eccodes version check
      #Only needed on SIRES computer, latest eccode version is available on HPC
      cd hub/local/src/FALFILFA/falfilfa
      git cherry-pick 15359c1
      cd $HOMEPACK/$ialdir
    fi

    #Prepare PHYEX inclusion
    rm -rf hub/local/src/PHYEX/phyex
    mkdir -p hub/local/src/PHYEX/phyex
  fi
fi
if [ $packupdate == true -o $packcreation == true ]; then
  cd $HOMEPACK/$ialdir/hub/local/src/PHYEX/phyex

  if [ $useexpand == true ]; then
    expand_options="--mnhExpand"
  else
    expand_options=""
  fi
  prep_code=$PHYEXTOOLSDIR/prep_code.sh
  echo "Copy $commit"
  mkdir PHYEX
  if [ $model_ready == false ]; then
    scp -q -r $commit/src PHYEX/
    subs="-s gmkpack_ignored_files -s turb -s micro -s aux -s conv  -s CMakeLists.txt -s cmake"
    $prep_code $prepCodeOpts $subs -m arome PHYEX -- --shumanFUNCtoCALL --removeACC $expand_options
  else
    scp -q -r $commit/* PHYEX/
    $prep_code $prepCodeOpts PHYEX #Ready for inclusion
  fi
  find PHYEX -type f -exec touch {} \; #to be sure a recompilation occurs
  if [ $packupdate == true ]; then
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
    if [ $fullcompilation == false ]; then
      #Content is added in the ics_masterodb script
      if [ $packupdate == false ]; then
        if [ $(cat PHYEX/gmkpack_ignored_files | wc -l) != 0 ]; then
          sed -i "/^end_of_ignored_files/i $(first=1; for line in $(cat PHYEX/gmkpack_ignored_files); do echo -n $(test $first -ne 1 && echo \\n)${line}; first=0; done)" $HOMEPACK/$ialdir/ics_masterodb
        fi
      fi
    else
      #Files must be suppressed (non phyex files)
      for file in $(cat PHYEX/gmkpack_ignored_files); do
        [ -f $HOMEPACK/$ialdir/src/local/$file ] && rm -f $HOMEPACK/$ialdir/src/local/$file
      done
      if [ -d $HOMEPACK/$ialdir/src/local/mpa/dummy ]; then
        [ ! "$(ls -A $HOMEPACK/$ialdir/src/local/mpa/dummy)" ] && rmdir $HOMEPACK/$ialdir/src/local/mpa/dummy
      fi
    fi
  fi

  rm -rf PHYEX

  if [ -d $HOMEPACK/$ialdir/hub/local/src/PHYEX/phyex -a -d $HOMEPACK/$ialdir/src/main ]; then
    # Incremental pack with updated Hub, we must put, in src/local, the files using
    # PHYEX to force a re-compilation
    cd $HOMEPACK/$ialdir/src

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

if [ $compilation == true ]; then
  if [ $onlyIfNeeded == false -o ! -f $HOMEPACK/$ialdir/bin/MASTERODB ]; then
    echo "### Compilation of commit $commit"

    cd $HOMEPACK/$ialdir
    sed -i 's/GMK_THREADS=1$/GMK_THREADS=10/' ics_masterodb
  
    [ -f ics_packages ] && submit Output_compilation_hub ics_packages
    submit Output_compilation ics_masterodb
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

if [ $run == true ]; then
  #Cleaning to suppress old results that may be confusing in case of a crash during the run
  if [ $onlyIfNeeded == false ]; then
    for t in $(echo $tests | sed 's/,/ /g'); do
      cd $HOMEPACK/$ialdir
      if [ -d conf_tests/$t ]; then
        rm -rf conf_tests/$t
      fi
    done
  fi

  #Run the tests one after the other
  firstrun=1
  for t in $(echo $tests | sed 's/,/ /g'); do #loop on tests
    if echo $allowedTests | grep -w $t > /dev/null; then #test is allowed on this plateform
      cd $HOMEPACK/$ialdir
      if [ ! -d conf_tests/$t ]; then #We do not enter systematically this part if onlyIfNeeded=true
        if [ $firstrun -eq 1 ]; then
          echo "### Running of commit $commit"
          firstrun=0
        fi

        if [ ! -f $HOMEPACK/$ialdir/bin/MASTERODB ]; then
          echo "Pack does not exist ($HOMEPACK/$ialdir) or compilation has failed, please check"
          exit 6
        fi

        mkdir -p conf_tests/$t
        cd conf_tests/$t
        t1=$(($(date +%s%N)/1000)) #current time in milliseconds
        MYLIB=$ialdir TESTDIR=$dirconf/$t submit Output_run $dirconf/$t/ar?.sh
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
if [ $check == true ]; then
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
      ref=${refByTest[$tag]}

      #Conversion into directory name
      if echo $ref | grep '/' > /dev/null; then
        refname="PHYEX/*_$(escape_commit $ref).01.${gmkpack_l}.${gmkpack_o}"
      else
        refname="PHYEX/*_${ref}.01.${gmkpack_l}.${gmkpack_o}"
      fi

      file1=$HOMEPACK/$ialdir/$file
      file2=$(echo $HOMEPACK/$refname/$file) #echo to enable shell substitution

      if [ ! -f $file2 -a $computeRefIfNeeded == true ]; then
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
          ddh_images="$HOMEPACK/$ialdir/ddh_diff_${tag}.png"
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

if [ $remove == true ]; then
  echo "### Remove model directory for commit $commit"
  [ -d $HOMEPACK/$ialdir ] && rm -rf $HOMEPACK/$ialdir
fi

exit $cmpstatus
