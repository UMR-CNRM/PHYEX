#!/bin/bash

set -e
set -o pipefail #abort if left command on a pipe fails

#This script:
# - compiles the LMDZ model using a specific commit for the externalised physics
# - runs RICO and ARM-CU 1D cases

PHYEXTOOLSDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
specialPack="ref"
availTests="rico arm_cu"
defaultTest="rico"
defaultRef='ref'
LMDZPACK=${LMDZPACK:=$HOME/LMDZ/PHYEX}
separator='_' #- be carrefull, gmkpack (at least on belenos) has multiple allergies (':', '.', '@')
              #- seprator must be in sync with prep_code.sh separator

function usage {
  echo "Usage: $0 [-h] [-p] [-u] [-c] [-C] [-r] [-s] [--expand] [-t test] [--version VERSION] [--repo-user] [--repo-protocol] [--remove] [--perf FILE] commit [reference]"
  echo "commit          commit hash (or a directory, or among $specialPack) to test"
  echo "reference       commit hash (or a directory, or among $specialPack) REF to use as a reference"
  echo "-s              suppress compilation pack"
  echo "-p              creates pack"
  echo "-u              updates pack (experimental, only for 'small' updates where"
  echo "                              the list of files does not change)"
  echo "-c              performs compilation"
  echo "-r              runs the tests"
  echo "-C              checks the result against the reference"
  echo "-t              comma separated list of tests to execute"
  echo "                or ALL to execute all tests"
  echo "--nofcm         don't use fcm (be carreful, with fcm compilation exits with a (false) error"
  echo "--expand        expand mnh_expand blocks (code will be in do loops)"
  echo "--version VERSION   to force using lmdz VERSION"
  echo "--repo-user     user hosting the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOuser (=$PHYEXREOuser)"
  echo "--repo-protocol protocol (https or ssh) to reach the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREOprotocol (=$PHYEXREOprotocol)"
  echo "--remove        removes the pack"
  echo "--perf FILE     add performance statistics in file FILE"
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
}

fcm=1
packcreation=0
packupdate=0
compilation=0
run=0
check=0
commit=""
reference=""
tests=""
suppress=0
useexpand=0
version=""
link=0 #Not yet put in command line argument because this option has not been tested here
remove=0
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
    '--nofcm') fcm=0;;
    '--expand') useexpand=1;;
    '--version') version="$2"; shift;;
    '--repo-user') export PHYEXREPOuser=$2; shift;;
    '--repo-protocol') export PHYEXREPOprotocol=$2; shift;;
    '--remove') remove=1;;
    '--perf') perffile="$(realpath $2)"; shift;;
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

if [ -z "${tests-}" ]; then
  tests=$defaultTest
elif [ $tests == 'ALL' ]; then
  tests=$availTests
fi

if [ $packcreation -eq 0 -a \
     $compilation -eq 0 -a \
     $packupdate -eq 0 -a \
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

function jsonarg {
  #$1 is the file containing the json dictionnary
  #$2 is the dictionnary key to return
  python3 -c "import json; print(json.load(open('$1', 'r'))['$2'])"
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

fromdir=''
if echo $commit | grep '/' | grep -v '^tags/' > /dev/null; then
  fromdir=$commit
  if [ "$version" == "" ]; then
    content_lmdz_version=$(scp $commit/src/lmdz/lmdz_version.json /dev/stdout 2>/dev/null)
    version=$(jsonarg <(echo $content_lmdz_version) version)
    rad=$(jsonarg <(echo $content_lmdz_version) rad)
    install_arg=$(jsonarg <(echo $content_lmdz_version) install_arg)
  fi
  name=$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')
  [ $suppress -eq 1 -a -d $LMDZPACK/$name ] && rm -rf $LMDZPACK/$name
elif echo $specialPack | grep -w $commit > /dev/null; then
  name="PHYEX/$commit"
else
  if [ "$version" == "" ]; then
    if [[ $commit == lmdz${separator}* ]]; then
      lmdz_version_file="lmdz_version.json"
    else
      lmdz_version_file="src/lmdz/lmdz_version.json"
    fi
    if echo $commit | grep '^tags/' > /dev/null; then
      urlcommit=$(echo $commit | cut -d / -f 2-)
    else
      urlcommit=$commit
    fi
    content_lmdz_version=$(wget --no-check-certificate https://raw.githubusercontent.com/$PHYEXREPOuser/PHYEX/${urlcommit}/$lmdz_version_file -O - 2>/dev/null || echo "")
    version=$(jsonarg <(echo $content_lmdz_version) version)
    rad=$(jsonarg <(echo $content_lmdz_version) rad)
    install_arg=$(jsonarg <(echo $content_lmdz_version) install_arg)
  fi
  name="COMMIT$(echo $commit | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')"
  [ $suppress -eq 1 -a -d $LMDZPACK/$name ] && rm -rf $LMDZPACK/$name
fi
if [ ! -z "${reference-}" ]; then
  [ $reference == 'REF' ] && reference=$defaultRef
  reffromdir=''
  if echo $reference | grep '/' > /dev/null; then
    reffromdir=$reference
    refname="PHYEX/$(echo $reference | sed 's/\//'${separator}'/g' | sed 's/:/'${separator}'/g' | sed 's/\./'${separator}'/g')"
  elif echo $specialPack | grep -w $reference > /dev/null; then
    refname="PHYEX/$reference"
  else
    refname="PHYEX/COMMIT${reference}"
  fi
fi

lmdzdir=$LMDZPACK/$name/LMDZ
phyexdir=$LMDZPACK/$name/PHYEX/
main=lmdz1d
L=79
##-debug is the default value and there's a bug (in the current script) if we try to specify it here
##compilecmd="./compile -L $L -rad $rad -cosp 0 -opt \"-debug \" -main $main"
compilecmd="./compile -L $L -rad $rad -cosp 0 -main $main"

if [ $packcreation -eq 1 ]; then
  echo "### Compilation of commit $commit"

  if echo $specialPack | grep -w $commit > /dev/null; then
    echo "Special commit '$commit' cannot be compiled with this script"
    exit 4
  fi

  if [ -d $LMDZPACK/$name ]; then
    echo "Pack already exists ($LMDZPACK/$name), suppress it to be able to compile it again (or use the -s option to automatically suppress it)"
    exit 5
  fi

  #Create directory
  mkdir -p $LMDZPACK/$name
  cd $LMDZPACK/$name

  #Populate with arch files
  touch arch-mylocal.env
  cat - <<EOF > arch-mylocal.fcm
%COMPILER            gfortran
%LINK                gfortran
%AR                  ar
%ARFLAGS             rU
%MAKE                make
%FPP_FLAGS           -P -traditional
%FPP_DEF             NC_DOUBLE
%BASE_FFLAGS          -cpp -ffree-line-length-0 -fdefault-real-8 -DNC_DOUBLE
%PROD_FFLAGS         -O3 -fallow-argument-mismatch
%DEV_FFLAGS          -Wall -fbounds-check  -fallow-argument-mismatch
%DEBUG_FFLAGS        -g3 -Wall -fbounds-check -ffpe-trap=invalid,zero,overflow -O0 -fstack-protector-all -fbacktrace -finit-real=snan  -fallow-argument-mismatch
%MPI_FFLAGS
%OMP_FFLAGS          
%BASE_LD              
%MPI_LD
%OMP_LD   
EOF

  cat - <<EOF > arch-mylocal.path
NETCDF_INCDIR="-I/usr/include"
NETCDF_LIBDIR=""
NETCDF_LIB="-lnetcdff -lnetcdf"

NETCDF95_INCDIR=-I\$LMDGCM/../../include/
NETCDF95_LIBDIR=-L\$LMDGCM/../../lib
NETCDF95_LIB=-lnetcdf95

IOIPSL_INCDIR="-I\$LMDGCM/../../lib -I\$LMDGCM/../IOIPSL/inc"
IOIPSL_LIBDIR="-L\$LMDGCM/../../lib -L\$LMDGCM/../IOIPSL/lib"
IOIPSL_LIB="-lioipsl"

XIOS_INCDIR="-I\$LMDGCM/../XIOS/inc"
XIOS_LIBDIR="-L\$LMDGCM/../XIOS/lib"
XIOS_LIB="-lxios -lstdc++"

ORCH_INCDIR="-I\$LMDGCM/../../lib"
ORCH_LIBDIR="-L\$LMDGCM/../../lib"

OASIS_INCDIR="-I\$LMDGCM/../../oasis3-mct/BLD/build/lib/psmile.MPI1"
OASIS_LIBDIR="-L\$LMDGCM/../../oasis3-mct/BLD/lib"
OASIS_LIB="-lpsmile.MPI1 -lscrip -lmct -lmpeu"

INCA_INCDIR="-I\$LMDGCM/../INCA/build/inc"
INCA_LIBDIR="-L\$LMDGCM/../INCA/build/lib"
INCA_LIB="-lchimie"
EOF

  #Compilation
  wget https://lmdz.lmd.jussieu.fr/pub/install_lmdz.sh -O install_lmdz.sh
  bash install_lmdz.sh -v $version $install_arg -bench 0 -rad $rad -name LMDZ -arch_dir $PWD -arch mylocal 2>&1 | tee Install.log

  #Populate with test cases (1D directory needed for compilation)
  cd $lmdzdir
  wget https://lmdz.lmd.jussieu.fr/pub/1D/1D.tar.gz
  tar xf 1D.tar.gz
fi
if [ $packupdate -eq 1 -o $packcreation -eq 1 ]; then
  #PHYEX code
  if [ $link -eq 1 ]; then
    if [ $packupdate -eq 1 ]; then
      echo "link option not compatible with the update option"
      exit 10
    fi
    #Special case when a PHYEX repository exist locally
    #We can link the LMDZ source tree with the PHYEX repository
    #This can be useful for debuging or developping
    cd ${lmdzdir}/modipsl/modeles/LMDZ/libf/phylmd/
    ln -s ~/PHYEX/src/common/*/* .
    ln -sf ~/PHYEX/src/lmdz/*/* .
  else
    #Checkout PHYEX
    if [ ! -d $LMDZPACK/$name ]; then
      echo "Pack directory doesn't exist ($LMDZPACK/$name)"
      exit 9
    fi
    cd $LMDZPACK/$name
    if [ $packupdate -eq 1 ]; then
      mv PHYEX PHYEXori
    fi
  
    if [ $useexpand == 1 ]; then
      expand_options="--mnhExpand"
    else
      expand_options=""
    fi
    subs="-s turb -s micro -s aux -s ext"
    prep_code=$PHYEXTOOLSDIR/prep_code.sh
    if [ "$fromdir" == '' ]; then
      echo "Clone repository, and checkout commit $commit (using prep_code.sh)"
      if [[ $commit == lmdz${separator}* ]]; then
        $prep_code -c $commit PHYEX #This commit is ready for inclusion
      else
        $prep_code -c $commit $expand_options $subs -m lmdz PHYEX -- --shumanFUNCtoCALL --removeACC
      fi
    else
      echo "Copy $fromdir"
      mkdir PHYEX
      scp -q -r $fromdir/src PHYEX/
      $prep_code $expand_options $subs -m lmdz PHYEX -- --shumanFUNCtoCALL --removeACC
    fi
    if [ $packupdate -eq 1 ]; then
      #Update only modified files
      cd PHYEX
      for file in $(find turb micro aux ext -type f); do
        mvdiff $file ../PHYEXori/$file
      done
      cd ..
      rm -rf PHYEX
      mv PHYEXori PHYEX
    else
      #Put PHYEX source code in the LMDZ source tree
      cd $lmdzdir/modipsl/modeles/LMDZ/libf/phylmd/
      ln -sf $phyexdir/*/* .
    fi
  fi

  #Update code
  cd $lmdzdir/modipsl/modeles/LMDZ/libf/phylmd/
  [ $packupdate -eq 0 ] && cp -r . ../phylmdorig
  if [ $fcm -eq 0 ]; then
    mvdiff modd_dimphyexn.F90 modd_dimphyex.F90
    for name in `grep -i 'END MODULE' modd*n.F90 | cut -d: -f1 | sed -e 's/n.F90//'` ; do
      mvdiff ${name}n.F90 ${name}_n.F90
    done
    mvdiff hypgeo.F90 modi_hypgeo.F90
    mvdiff hypser.f90 modi_hypser.F90
    mvdiff momg.F90 modi_momg.F90
    mvdiff tools.F90 mode_tools.F90
    [ -f shuman_mf.F90 ] && mvdiff shuman_mf.F90 modi_shuman_mf.F90
    [ -f shuman_phy.F90 ] && mvdiff shuman_phy.F90 mode_shuman_phy.F90
  fi

  #Missing files in case ecrad is not used
  if [ $packupdate -eq 0 ]; then
    if [ "$rad" != "ecrad" ] ; then
      for file in ecrad/yom* ecrad/abor1.F90 ecrad/abor1.intfb.h ecrad/parkind1.F90; do
        [ ! -f $(basename $file) ] && ln -s $file .
      done
    fi
  fi
fi

if [ $compilation -eq 1 ]; then
  echo "### Compilation of commit $commit"
  cd $lmdzdir/1D/bin
  rm -f $lmdzdir/1D/bin/${main}.e
  rm -f $lmdzdir/modipsl/modeles/LMDZ/bin/lmdz1d_${L}_phylmd_${rad}_seq.e
  if [ $fcm -eq 1 ]; then
    sed -i "s/fcm=0/fcm=1/g" compile
  fi
  $compilecmd 2>&1 | tee $lmdzdir/compilation.log
  if [ $fcm -eq 1 ]; then
    echo "Using fcm, compilation exits with error even if everything is OK"
  fi
fi

if [ $run -eq 1 ]; then
  echo "### Execution of commit $commit"
  cd $lmdzdir/1D/INPUT/PHYS
  sed '1 i\iflag_physiq=1\n' physiq.def_6A > physiq.def_PHYLMD
  sed '1 i\iflag_physiq=2\n' physiq.def_6A > physiq.def_PHYEX
  for t in $tests; do
    for DEF in PHYEX PHYLMD; do
      d=${lmdzdir}/1D/EXEC/${DEF}L$L/$t
      [ ! -d $d ] && mkdir -p $d
      cd $d
      ln -sf ${lmdzdir}/1D/OLDCASES/$t/* .
      cp -f ${lmdzdir}/1D/INPUT/DEF/*.def .
      cp -f ${lmdzdir}/1D/INPUT/PHYS/physiq.def_$DEF physiq.def
      if [ $rad = oldrad ] ; then
        sed -i'' -e 's/iflag_rrtm=.*$/iflag_rrtm=0/' -e 's/NSW=.*$/NSW=2/' physiq.def
      fi
      if [ $rad = ecrad ] ; then
        cp -f $lmdzdir/modipsl/modeles/LMDZ/DefLists/namelist_ecrad .
        cp -rf $lmdzdir/modipsl/modeles/LMDZ/libf/phylmd/ecrad/data .
        sed -e 's@iflag_rrtm=1@iflag_rrtm=2@' physiq.def > tmp
        \mv tmp physiq.def
      fi
      cp -f ${lmdzdir}/1D/INPUT/VERT/L$L/* .
      ln -sf L$L.def vert.def
      set +e
      cp -f $lmdzdir/1D/OLDCASES/$t/*.d[ae]* .
      set -e
      cat <<......eod>| compile.sh
       cd $lmdzdir/1D/bin
       $compilecmd
......eod
      chmod +x compile.sh
  
      if [ $fcm -eq 0 ]; then
        ln -sf $lmdzdir/1D/bin/${main}.e ${main}.e
      else
        execname=$lmdzdir/modipsl/modeles/LMDZ/bin/lmdz1d_${L}_phylmd_${rad}_seq.e
        if [ -f $execname ]; then
          ln -sf $execname ${main}.e
        else
          ln -sf $lmdzdir/modipsl/modeles/LMDZ/lmdz1d_${rad}_L${L}.e ${main}.e
        fi
      fi
  
      if [ $DEF == PHYEX ]; then
        sed -i -e 's/day_step=144$/day_step=1440/' gcm1d.def
      fi
      t1=$(($(date +%s%N)/1000)) #current time in milliseconds
      ./lmdz1d.e 2>&1 | tee execution.log
      t2=$(($(date +%s%N)/1000))
      if [ "$perffile" != "" ]; then
        echo "$commit lmdz $t $(($t2-$t1))" >> "$perffile"
      fi
    done
  done
fi

if [ $check -eq 1 ]; then
  echo "### Check commit $commit against commit $reference"
  echo "This functionnality is not yet implemented because:"
  echo "  1) the PHYEX interface will evolve and bit-reproducibility is not guaranted"
  echo "  2) there is no surface scheme plugged in LMZD-PHYEX, a recompilation"
  echo "     must be done to change the surface fluxes"
  exit 6
fi

if [ $remove -eq 1 ]; then
  echo "### Remove model directory for commit $commit"
  [ -d $LMDZPACK/$name ] && rm -rf $LMDZPACK/$name
fi
