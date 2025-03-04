# This file contains the generic functions needed to compile PHYEX
# It cannot be used directly, instead you can use script under one of the with_* subdirectories

fiat_version=5eef5552c3002aa962caef56c6bdc88932739e77 #this specific version is needed for NEC
fiat_gh_user=ACCORD-NWP #for official repo, use ecmwf-ifs

function parse_args() {
  # default values
  ARCH_PATH=$PWD/arch
  ARCH=
  GMKFILE=
  MESONHPROFILE=
  useexpand=1
  commit=""
  packcreation=0
  packupdate=0
  compilation=0
  inplaceClean=0
  inplaceInstall=0
  ssh=0
  # pass unrecognized arguments to the build system
  BS_ARGS=""
  
  while (($# > 0)); do
    OPTION="$1" ; shift
    case "$OPTION" in
      "-h") cat <<EOF
	    Usage :
$0 [options]
--help -h   help
--arch-path ARCH_PATH directory for architecture specific files (see below) [./arch]
                      note that arch files are first looked for in ${HOME}/.phyex/<bs>_arch where <bs> is the build system
--arch ARCH  	        build using arch file [gnu]
--gmkfile FILE        build using a gmkpack configuration file (--arch must be used to give a name to the build dir)
--mesonhprofile FILE  build using Méso-NH profile and rules (--arch must be used to give a name to the build dir)
--noexpand            do not use mnh_expand (code will be in array-syntax)"
--commit              commit hash (or a directory) to test; do not use this option from within a repository
-p                    creates 'pack' (compilation directory)
-u                    updates 'pack'
-c                    performs compilation
--inplace-install     install or update, if needed, fiat and the build system in the directory where the current script is
--inplace-clean       remove the fiat and the build system installation present in the directory where the current script is
--ssh                 use the ssh protocol to clone the pyfortool and fxtran repositories instead of https"

Unrecognized options are passed to the the build system. Useful options for FCM include:
--new                   clean build tree before building
--jobs=N                parallel build, similar to make -j N
--ignore-lock           ignore lock indicating another build is ongoing, useful after an interrupted build

For details on FCM, see 
    http://metomi.github.io/fcm/doc/user_guide/build.html
    http://metomi.github.io/fcm/doc/user_guide/command_ref.html#fcm-build

If neither creation nor execution is requested, both steps are performed.
EOF
        exit;;
      "--arch")
        ARCH=$1 ; shift ;; 
      "--arch-path")
        ARCH_PATH=$1 ; shift ;; 
      "--gmkfile")
        GMKFILE=$1 ; shift ;;
      "--mesonhprofile")
        MESONHPROFILE=$1 ; shift ;;
      '--noexpand') useexpand=0;;
      '--commit') commit=$1; shift;;
      '-p') packcreation=1;;
      '-u') packupdate=1;;
      '-c') compilation=1;;
      '--inplace-install') inplaceInstall=1;;
      '--inplace-clean') inplaceClean=1;;
      '--ssh') ssh=1;;
      *)
        BS_ARGS="$BS_ARGS $OPTION" ;;
    esac
  done
  [ "$GMKFILE" == "" -a "$MESONHPROFILE" == "" -a "$ARCH" == "" ] && ARCH=gnu
  if [ "$GMKFILE" != "" -a "$ARCH" == "" ]; then
    echo "--arch option is mandatory if --gmkfile option is used"
    exit 2
  fi
  if [ "$MESONHPROFILE" != "" -a "$ARCH" == "" ]; then
    echo "--arch option is mandatory if --mesonhprofile option is used"
    exit 3
  fi
  if [ $inplaceInstall -eq 0 -a \
       $inplaceClean -eq 0 -a \
       $packcreation -eq 0 -a \
       $packupdate -eq 0 -a \
       $compilation -eq 0 ]; then
    packcreation=1
    compilation=1
  fi
}

function check_install_fiat() {
  if [ $ssh -eq 1 ]; then
    repo=git@github.com:$fiat_gh_user/fiat.git
    remote=ssh_$fiat_gh_user
  else
    repo=https://github.com/$fiat_gh_user/fiat.git
    remote=https_$fiat_gh_user
  fi
  cd fiat
  if [ ! -d src ]; then
    echo "Performing fiat cloning..."
    rm -f .gitkeep
    git clone $repo .
    git remote rename origin $remote
    touch .gitkeep
    echo "...fiat cloning done"
  fi
  if [ $(git remote -v | grep -e "^$remote\s" | wc -l) -eq 0 ]; then
    #the repository and/or the protocol is new
    git remote add $remote $repo
  fi
  #Checkout the right version
  set +e
  #try to directly checkout to reduce network need
  git checkout $fiat_version 2> /dev/null
  stat=$?
  set -e
  if [ $stat -ne 0 ]; then
    echo "Performing fiat fetching..."
    git fetch $remote
    git checkout $fiat_version
    echo "...fiat fetching done"
  fi
  cd ..
  echo
}

####################################

function main() {
  BS=$1 # build system
  shift
  DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/with_$BS

  # Parse command line arguments
  parse_args $*

  if [ $inplaceClean -eq 1 ]; then
    # Change current working dir
    cd $DIR
  
    # Restore build system
    remove_build_system
  
    # Restore fiat
    rm -rf fiat
  fi
  
  if [ $inplaceInstall -eq 1 ]; then
    # Change current working dir
    cd $DIR
  
    # Check the build system installation
    check_install_build_system
  
    # Check the fiat installation
    check_install_fiat
  fi
  
  builddir=arch_$ARCH
  if [ $packcreation -eq 1 ]; then
    # Change current working dir
    cd -P $DIR
    
    # Check the build system installation
    check_install_build_system
    
    # Check the fiat installation
    check_install_fiat
    
    # Create the build directory and set up the build system
    if [ -d $builddir ]; then
      echo "$builddir already exists. To rerun compilation, please enter this directory and use the compilation.sh script."
      echo "Otherwise, you can remove the $builddir directory and execute again this script."
      exit 1
    fi
    mkdir $builddir
    if [ "$GMKFILE" != "" ]; then
      touch $builddir/arch.env
      gmkfile2arch $GMKFILE $builddir
    elif [ "$MESONHPROFILE" != "" ]; then
      touch $builddir/arch.env
      mesonhprofile2arch $MESONHPROFILE $builddir
    else
      setuparch ${ARCH} $builddir
    fi
  fi
  if [ $packupdate -eq 1 -o $packcreation -eq 1 ]; then
    if [ ! -d $builddir ]; then
      echo "$builddir doesn't exist"
      exit 4
    fi
    cd $builddir
    if [ $packupdate -eq 1 ]; then
      rm -rf src
    fi
    . arch.env
    [ -z ${PYFT_OPTS+x} ] && PYFT_OPTS=''
    export PYFT_OPTS
    
    # Populate the source directory with (modified) PHYEX source code
    [ "$commit" == "" ] && commit=$PWD/../../.. #Current script run from within a PHYEX repository
    if echo $commit | grep '/' | grep -v '^tags/' > /dev/null; then
      # We get the source code directly from a directory
      fromdir=$commit
    else
      # We use a commit to checkout
      fromdir=''
    fi
    #Expand options
    if [ $useexpand == 1 ]; then
      expand_options="--mnhExpand"
    else
      expand_options=""
    fi
    PHYEXTOOLSDIR="$DIR/../../../tools" #if run from within a PHYEX repository
    UPDATEDPATH=$PATH
  
    #Temporary file for the description tree
    descTree=${TMPDIR:-/tmp}/descTree_$$
    trap "\rm -f $descTree" EXIT
  
    which prep_code.sh > /dev/null || export UPDATEDPATH=$PHYEXTOOLSDIR:$PATH
    subs="$subs -s turb -s shallow -s turb_mnh -s micro -s aux -s ice_adjust -s rain_ice -s rain_ice_old -s support -s progs $EXTRASUBS"
    if [ "$fromdir" == '' ]; then
      echo "Clone repository, and checkout commit $commit (using prep_code.sh)"
      if [[ $commit == testprogs${separator}* || $commit == offline${separator}* ]]; then
        #This commit is ready for inclusion
        PATH=$UPDATEDPATH prep_code.sh -c $commit src
      else
        PATH=$UPDATEDPATH prep_code.sh --pyfortool_opts_env PYFT_OPTS -c $commit $expand_options $subs \
                                       -m offline src --useParallelPyForTool -- --tree . --descTree $descTree
      fi
    else
      echo "Copy $fromdir"
      mkdir src
      scp -q -r $fromdir/src src/
      PATH=$UPDATEDPATH prep_code.sh --pyfortool_opts_env PYFT_OPTS $expand_options $subs \
                                     -m offline src --useParallelPyForTool -- --tree . --descTree $descTree
    fi
    
    # Add some code
    mkdir -p build/bin
    cd src
    if [ ${PYBINDING-yes} == 'yes' ]; then
      pybinding.py micro/ice_adjust.F90 sub:ICE_ADJUST pyphyex.F90 ../build/bin/pyphyex.py ./../lib/libphyex.so
      pybinding.py micro/rain_ice.F90 sub:RAIN_ICE pyphyex.F90 ../build/bin/pyphyex.py ./../lib/libphyex.so
      pybinding.py micro/mode_ice4_sedimentation.F90 \
                   module:MODE_ICE4_SEDIMENTATION/sub:ICE4_SEDIMENTATION \
                   pyphyex.F90 ../build/bin/pyphyex.py ./../lib/libphyex.so
      pybinding.py micro/rain_ice_old.F90 sub:RAIN_ICE_OLD pyphyex.F90 ../build/bin/pyphyex.py ./../lib/libphyex.so
      pybinding.py turb/shallow_mf.F90 sub:SHALLOW_MF pyphyex.F90 ../build/bin/pyphyex.py ./../lib/libphyex.so
      pybinding.py turb/turb.F90 sub:TURB pyphyex.F90 ../build/bin/pyphyex.py ./../lib/libphyex.so
      pybinding.py aux/ini_phyex.F90 sub:INI_PHYEX pyphyex.F90 ../build/bin/pyphyex.py ./../lib/libphyex.so
      pybinding.py micro/lima_adjust_split.F90 sub:LIMA_ADJUST_SPLIT pyphyex.F90 \
                   ../build/bin/pyphyex.py ./../lib/libphyex.so
      pybinding.py micro/lima.F90 sub:LIMA pyphyex.F90 ../build/bin/pyphyex.py ./../lib/libphyex.so
    else
      cat <<......EOF > pyphyex.F90
      SUBROUTINE PYPHYEXSUB
      END SUBROUTINE PYPHYEXSUB
......EOF
    fi
    ln -s ../../fiat fiat
    cat <<....EOF > dummyprog.F90
    PROGRAM DUMMYPROG
      PRINT*, "CREATED TO FORCE FCM TO LINK SOMETHING"
    END PROGRAM DUMMYPROG
....EOF
    #needed with commit 6b9b61b3d17228fcb5c0186e38d72aea987acd10 by P. Marguinaud
    #due to a weakness of fcm
    cat <<....EOF > nvtx_dummy.F90
    MODULE NVTX
      !Unused module but wrongly detected as dependency by fcm
      !whereas it is within an ifdef directive
    END MODULE NVTX
....EOF
    #needed with commit 6dc52c45fef39e1202a281eea21f7bdda047db28 due to the same weakness of fcm
    if [ ! -f aux/modd_util_phyex_t.F90 ]; then
      cat <<......EOF > aux/modd_util_phyex_t.F90
      MODULE MODD_UTIL_PHYEX_T
      !Unused module but wrongly detected as dependency by fcm
      !whereas it is within an ifdef directive
      USE MODD_PHYEX, ONLY: PHYEX_T
      CONTAINS
        SUBROUTINE COPY_PHYEX_T (YD, LDCREATED)
          TYPE (PHYEX_T), INTENT(IN), TARGET :: YD
          LOGICAL, OPTIONAL, INTENT(IN) :: LDCREATED
        END SUBROUTINE COPY_PHYEX_T
        SUBROUTINE WIPE_PHYEX_T (YD, LDDELETED)
          TYPE (PHYEX_T), INTENT(IN), TARGET :: YD
          LOGICAL, OPTIONAL, INTENT(IN) :: LDDELETED
        END SUBROUTINE WIPE_PHYEX_T
      END MODULE MODD_UTIL_PHYEX_T
......EOF
    fi
  fi
  
  # Build the compilation script and run it
  if [ $compilation -eq 1 ]; then
    cd -P $DIR/arch_$ARCH
    . arch.env
    build_compilation_script src
    ./compilation.sh
    
    # Check if python can open the resulting shared lib
    set +e
    #On NEC, the shared library cannot be loaded as other lib if it was compiled for the Vector Engine
    python3 -c "from ctypes import cdll; cdll.LoadLibrary('./build/lib/libphyex.so')"
    if [ $? -ne 0 ]; then
      echo "On some systems (cross-compilation) it's normal to obtain an error here"
      echo "when python tries to open the shared lib."
    fi
    set -e
    
    # ldd -r ./build/lib/libphyex.so should also give interesting results
  fi
}
