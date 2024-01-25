#!/bin/bash

set -e
#set -x

fcm_version=tags/2021.05.0
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
  # pass unrecognized arguments to fcm
  FCM_ARGS=""
  
  while (($# > 0)); do
    OPTION="$1" ; shift
    case "$OPTION" in
      "-h") cat <<EOF
	    Usage :
$0 [options]
--help -h   help
--arch-path ARCH_PATH directory for architecture specific files (see below) [./arch]
                      note that arch files are first looked for in ${HOME}/.phyex/fcm_arch
--arch ARCH  	        build using arch file arch-ARCH.fcm [gnu]
--gmkfile FILE        build using a gmkpack configuration file (--arch must be used to give a name to the build dir)
--mesonhprofile FILE  build using Méso-NH profile and rules (--arch must be used to give a name to the build dir)
--noexpand            do not use mnh_expand (code will be in array-syntax)"
--commit              commit hash (or a directory) to test; do not use this option from within a repository
-p                    creates 'pack' (compilation directory)
-u                    updates 'pack'
-c                    performs compilation
--inplace-install     install, if needed, fiat and fcm in the directory where the current script is
--inplace-clean       remove the fiat and fcm installation present in the directory where the current script is
--ssh                 use the ssh protocol to clone the pyft and fxtran repositories instead of https"

Unrecognized options are passed to the fcm build command. Useful options include :
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
        FCM_ARGS="$FCM_ARGS $OPTION" ;;
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

function check_install_fcm() {
  if [ ! -f fcm/bin/fcm ]; then
    echo "Performing FCM installation..."
    cd fcm
    rm -f .gitkeep
    if [ $ssh -eq 1 ]; then
      git clone git@github.com:metomi/fcm.git .
    else
      git clone https://github.com/metomi/fcm.git .
    fi
    git checkout $fcm_version
    touch .gitkeep
    cd ..
    echo "...FCM installation done"
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

function gmkfile2arch() {
  GMKFILE=$1
  ARCHFILE=$2
cat <<EOF > $ARCHFILE
# Compilation
\$FCOMPILER     =     $(grep "^FRTNAME =" $GMKFILE | cut -d = -f 2)
\$BASE_FFLAGS   =     $(grep "^FRTFLAGS =" $GMKFILE | cut -d = -f 2-) $(grep "^GMK_FCFLAGS_PHYEX =" $GMKFILE | cut -d = -f 2-)
\$PROD_FFLAGS   =     $(grep "^OPT_FRTFLAGS =" $GMKFILE | cut -d = -f 2-)
\$DEV_FFLAGS    =     $(grep "^DBG_FRTFLAGS =" $GMKFILE | cut -d = -f 2-)
\$DEBUG_FFLAGS  =     $(grep "^DBG_FRTFLAGS =" $GMKFILE | cut -d = -f 2-) $(grep "^BCD_FRTFLAGS =" $GMKFILE | cut -d = -f 2-) $(grep "^NAN_FRTFLAGS =" $GMKFILE | cut -d = -f 2-)
\$CCOMPILER     =     $(grep "^VCCNAME =" $GMKFILE | cut -d = -f 2)
\$BASE_CFLAGS   =     $(grep "^VCCFLAGS =" $GMKFILE | cut -d = -f 2-)
\$PROD_CFLAGS   =     $(grep "^OPT_VCCFLAGS =" $GMKFILE | cut -d = -f 2-)
\$DEV_CFLAGS    =     
\$DEBUG_CFLAGS  =     
\$OMP_FFLAGS    =

# Preprocessor
\$FPP_FLAGS     =     $(grep "^MACROS_FRT =" $GMKFILE | cut -d = -f 2- | sed 's/-D//g')
\$CPP_FLAGS     =     $(grep "^MACROS_CC =" $GMKFILE | cut -d = -f 2- | sed 's/-D//g')

# Linker
\$LINK          =     $(grep "^LNK_MPI =" $GMKFILE | cut -d = -f 2-)
\$BASE_LD       =     $(grep "^LNK_FLAGS =" $GMKFILE | cut -d = -f 2-)
\$OMP_LD        =
\$LD_EXE_TO_SHARED = $(grep "^LNK_SOLIB =" $GMKFILE | cut -d = -f 2-  | sed 's/-o a.out//')

# Other
\$AR            =     $(grep "^AR =" $GMKFILE | cut -d = -f 2-)
EOF
}

function mesonhprofile2archenv() {
  MESONHPROFILE=$1
  ARCHFILE=$2
  ENVFILE=$3

  echo "
   You are trying to produce a configuration file for fcm from a Meso-NH configuration.
   The resulting file is certainly incomplete and must be modified as follows:
      Optimisation level:
        The opt level is set in the mesonh profile file; as a consequence, the BASE_FFLAGS contains
        the base *and* the opt flags.
        To compile with other opt level, the profile file must be modified before executing this function.
      Long lines:
        Meso-NH rules does not allow the compilation of long lines. Depending on compilers, it might be needed to
        manually add an option to allow long lines.
        For gfortran: add '-ffree-line-length-none' to BASE_FFLAGS
      OpenMP:
        Meso-NH does not use OpenMP but testprogs do; as a consequence, openmp flags are not included in the
        Meso-NH rules, they must be manually added.
        For gfortran: add '-fopenmp' to BASE_FFLAGS and to BASE_LD
      Position Independent Code:
        Meso-NH does not need to build position independent code, flags must be set manually.
        For gfortran ('-fPIC' already in BASE_FFLAGS): add '-fPIC' to BASE_CFLAGS
      Shared lib:
        Flags needed to build shared lib are not defined in Meso-NH rules, only hard coded in Makefile to build a
        specific lib. The flags to set for building a shared lib, in addition to flags used to build an object, must
        be manually set.
        For gfortran: add '-shared' to LD_EXE_TO_SHARED
      Swap:
        Meso-NH rules does not swap IO byte order (litle-/big-endian). Depending on your endianess, the
        corresponding flag may have to be set manually.
        For gfortran: add '-fconvert=swap' to BASE_FFLAGS"
  tac $MESONHPROFILE | grep -m1 '#' -B $(cat $MESONHPROFILE | wc -l) | tac | grep -v '#' > $ENVFILE
  MAKEFILE='
include Rules.$(ARCH)$(F).mk

archfile :
	echo "# Compilation"
	echo "\$$FCOMPILER     =     $(F90)"
	echo "\$$BASE_FFLAGS   =     -c $(F90FLAGS)"
	echo "\$$PROD_FFLAGS   =     "
	echo "\$$DEV_FFLAGS    =     "
	echo "\$$DEBUG_FFLAGS  =     "
	echo "\$$CCOMPILER     =     $(CC)"
	echo "\$$BASE_CFLAGS   =     -c $(CFLAGS)"
	echo "\$$PROD_CFLAGS   =     "
	echo "\$$DEV_CFLAGS    =     "
	echo "\$$DEBUG_CFLAGS  =     "
	echo "\$$OMP_FFLAGS    ="
	echo ""
	echo "# Preprocessor"
	echo "\$$FPP_FLAGS     =     $(CPPFLAGS)"
	echo "\$$CPP_FLAGS     =     $(CPPFLAGS)" 
	echo ""
	echo "# Linker"
	echo "\$$LINK          =     $(FC)"
	echo "\$$BASE_LD       =     $(LDFLAGS)"
	echo "\$$OMP_LD        ="
	echo "\$$LD_EXE_TO_SHARED =  "
	echo ""
	echo "# Other" 
	echo "\$$AR            =     $(AR)"

'
  (. $MESONHPROFILE; make -f <(echo -e "$MAKEFILE") -s -I $(dirname $MESONHPROFILE)/../src archfile) | sed 's/-D//g' > $ARCHFILE
}

function build_compilation_script() {
srcdir=$1

#fcm doesn't like if a source directory doesn't exist.
#To be able to compile an old commit, we must filter the source directories
TESTPROGS_DIR=""
#support is not a testprog but is needed
for testprog in ice_adjust rain_ice turb_mnh shallow rain_ice_old support; do
  [ -d $srcdir/$testprog ] && TESTPROGS_DIR+="src/$testprog "
done

cat <<EOF > compilation.sh
#!/bin/bash

. arch.env

level=PROD #PROD DEV or DEBUG

#fcm variables begin with a dollar sign

COMPIL_FFLAGS="\\\$\${level}_FFLAGS"
COMPIL_FFLAGS="\$COMPIL_FFLAGS \\\$OMP_FFLAGS"

COMPIL_CFLAGS="\\\$\${level}_CFLAGS"

LD_FLAGS="\\\$BASE_LD"
LD_FLAGS="\$LD_FLAGS \$OMP_LD"

LIBS="${LIBS:-rt dl}"

ENTRYPOINTS="rain_ice.o shallow_mf.o turb.o ice_adjust.o"

FCM_ARGS="$FCM_ARGS"

echo "\\\$COMPIL_FFLAGS = \$COMPIL_FFLAGS" > config.fcm
echo "\\\$COMPIL_CFLAGS = \$COMPIL_CFLAGS" >> config.fcm
echo "\\\$LD_FLAGS = \$LD_FLAGS" >> config.fcm
echo "\\\$ENTRYPOINTS = \$ENTRYPOINTS" >> config.fcm
echo "\\\$LIBS = \$LIBS" >> config.fcm
echo "\\\$TESTPROGS_DIR=$TESTPROGS_DIR" >> config.fcm

export PATH=$PWD/../fcm/bin/:\$PATH

echo "This script has generated config.fcm which is included by fcm-make.cfg, the FCM configuration file."
echo "Running : fcm make \$FCM_ARGS"

fcm make \$FCM_ARGS
EOF
chmod +x compilation.sh
}

####################################

# Where we are
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Parse command line arguments
parse_args $*

if [ $inplaceClean -eq 1 ]; then
  # Change current working dir
  cd $DIR

  # Restore fcm
  rm -rf fcm
  git restore fcm

  # Restore fiat
  rm -rf fiat
  git restore fiat
fi

if [ $inplaceInstall -eq 1 ]; then
  # Change current working dir
  cd $DIR

  # Check the fcm installation
  check_install_fcm

  # Check the fiat installation
  check_install_fiat
fi

builddir=arch_$ARCH
if [ $packcreation -eq 1 ]; then
  # Change current working dir
  cd -P $DIR
  
  # Check the fcm installation
  check_install_fcm
  
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
    gmkfile2arch $GMKFILE $builddir/arch.fcm
  elif [ "$MESONHPROFILE" != "" ]; then
    touch $builddir/arch.env
    mesonhprofile2archenv $MESONHPROFILE $builddir/arch.fcm $builddir/arch.env
  else
    if [ -f ${HOME}/.phyex/fcm_arch/arch-${ARCH}.fcm ]; then
      cp ${HOME}/.phyex/fcm_arch/arch-${ARCH}.env $builddir/arch.env
      cp ${HOME}/.phyex/fcm_arch/arch-${ARCH}.fcm $builddir/arch.fcm
    else 
      cp ${ARCH_PATH}/arch-${ARCH}.env $builddir/arch.env
      cp ${ARCH_PATH}/arch-${ARCH}.fcm $builddir/arch.fcm 
    fi
  fi
  cp fcm-make.cfg $builddir
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
  which prep_code.sh > /dev/null || export UPDATEDPATH=$PHYEXTOOLSDIR:$PATH
  subs="$subs -s turb -s shallow -s turb_mnh -s micro -s aux -s ice_adjust -s rain_ice -s rain_ice_old -s support"
  if [ "$fromdir" == '' ]; then
    echo "Clone repository, and checkout commit $commit (using prep_code.sh)"
    if [[ $commit == testprogs${separator}* ]]; then
      PATH=$UPDATEDPATH prep_code.sh --pyft_opts_env PYFT_OPTS -c $commit src #This commit is ready for inclusion
    else
      PATH=$UPDATEDPATH prep_code.sh --pyft_opts_env PYFT_OPTS -c $commit $expand_options $subs -m testprogs src
    fi
  else
    echo "Copy $fromdir"
    mkdir src
    scp -q -r $fromdir/src src/
    PATH=$UPDATEDPATH prep_code.sh --pyft_opts_env PYFT_OPTS $expand_options $subs -m testprogs src
  fi
  
  # Add some code
  cd src
  ln -s ../../fiat/src fiat
  cat <<..EOF > dummyprog.F90
  PROGRAM DUMMYPROG
    PRINT*, "CREATED TO FORCE FCM TO LINK SOMETHING"
  END PROGRAM DUMMYPROG
..EOF
  #needed with commit 6b9b61b3d17228fcb5c0186e38d72aea987acd10 by P. Marguinaud
  #due to a weakness of fcm
  cat <<..EOF > nvtx_dummy.F90
  MODULE NVTX
    !Unused module but wrongly detected as dependency by fcm
    !whereas it is within an ifdef directive
  END MODULE NVTX
..EOF
fi

# Build the compilation script and run it
if [ $compilation -eq 1 ]; then
  cd -P $DIR/arch_$ARCH
  . arch.env
  build_compilation_script src
  ./compilation.sh
  [ ! -f libphyex.so ] && ln -s build/bin/libphyex.so .
  
  # Check if python can open the resulting shared lib
  set +e
  #On NEC, the shared library cannot be loaded as other lib if it was compiled for the Vector Engine
  python3 -c "from ctypes import cdll; cdll.LoadLibrary('./libphyex.so')"
  if [ $? -ne 0 ]; then
    echo "On some systems (cross-compilation) it's normal to obtain an error here"
    echo "when python tries to open the shared lib."
  fi
  set -e
  
  # ldd -r ./libphyex.so should also give interesting results
fi
