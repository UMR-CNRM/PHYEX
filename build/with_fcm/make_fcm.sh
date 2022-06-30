#!/bin/bash

set -e

fcm_version=tags/2021.05.0
fiat_version=1295120464c3905e5edcbb887e4921686653eab8

function parse_args() {
  # default values
  ARCH_PATH=$PWD/arch
  ARCH=
  GMKFILE=
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
--arch ARCH  	        build using arch files $ARCH_PATH/arch-ARCH.* [gnu]
--gmkfile FILE        build using a gmkpack configuration file (--arch must be used to give a name to the build dir)

Unrecognized options are passed to the fcm build command. Useful options include :
--new                   clean build tree before building
--jobs=N                parallel build, similar to make -j N
--ignore-lock           ignore lock indicating another build is ongoing, useful after an interrupted build

For details on FCM, see 
    http://metomi.github.io/fcm/doc/user_guide/build.html
    http://metomi.github.io/fcm/doc/user_guide/command_ref.html#fcm-build
EOF
        exit;;
      "--arch")
        ARCH=$1 ; shift ;; 
      "--arch-path")
        ARCH_PATH=$1 ; shift ;; 
      "--gmkfile")
        GMKFILE=$1 ; shift ;;
      *)
        FCM_ARGS="$FCM_ARGS $OPTION" ;;
    esac
  done
  [ "$GMKFILE" == "" -a "$ARCH" == "" ] && ARCH=gnu
  if [ "$GMKFILE" != "" -a "$ARCH" == "" ]; then
    echo "--arch option is mandatory if --gmkfile option is used"
    exit 2
  fi
}

function check_install_fcm() {
  if [ ! -f fcm/bin/fcm ]; then
    echo "Performing FCM installation..."
    cd fcm
    rm -f .gitkeep
    git clone https://github.com/metomi/fcm.git .
    git checkout $fcm_version
    touch .gitkeep
    cd ..
    echo "...FCM installation done"
  fi
}

function check_install_fiat() {
  if [ ! -d fiat/src ]; then
    echo "Performing fiat cloning..."
    cd fiat
    rm -f .gitkeep
    git clone https://github.com/ecmwf-ifs/fiat.git .
    git checkout $fiat_version
    touch .gitkeep
    cd ..
    echo "...fiat cloning done"
  fi
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

function build_compilation_script() {
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

LIBS="rt dl"

ENTRYPOINTS="rain_ice.o shallow_mf.o turb.o ice_adjust.o ini_neb.o"

FCM_ARGS="$FCM_ARGS"

echo "\\\$COMPIL_FFLAGS = \$COMPIL_FFLAGS" > config.fcm
echo "\\\$COMPIL_CFLAGS = \$COMPIL_CFLAGS" >> config.fcm
echo "\\\$LD_FLAGS = \$LD_FLAGS" >> config.fcm
echo "\\\$ENTRYPOINTS = \$ENTRYPOINTS" >> config.fcm
echo "\\\$LIBS = \$LIBS" >> config.fcm

export PATH=$PWD/../fcm/bin/:\$PATH

echo "This script has generated config.fcm which is included by fcm-make.cfg, the FCM configuration file."
echo "Running : fcm make \$FCM_ARGS"

fcm make \$FCM_ARGS
EOF
chmod +x compilation.sh
}

####################################

# Parse command line arguments
parse_args $*

# Change current working dir
cd -P $(dirname $0)

# Check the fcm installation
check_install_fcm

# Check the fiat installation
check_install_fiat

# Create the build directory and populate it
builddir=arch_$ARCH
if [ -d $builddir ]; then
  echo "$builddir already exists. To rerun compilation, please enter this directory and use the compilation.sh script."
  echo "Otherwise, you can remove the $builddir directory and execute again this script."
  exit 1
fi
mkdir $builddir
if [ "$GMKFILE" != "" ]; then
  touch $builddir/arch.env
  gmkfile2arch $GMKFILE $builddir/arch.fcm
else
  cp ${ARCH_PATH}/arch-${ARCH}.env $builddir/arch.env
  cp ${ARCH_PATH}/arch-${ARCH}.fcm $builddir/arch.fcm 
fi
cp fcm-make.cfg $builddir
cd $builddir
mkdir src
cd src
ln -s ../../../../src/common .
ln -s ../../../../src/testprogs .
ln -s ../../fiat/src fiat
cat <<EOF > dummyprog.F90
PROGRAM DUMMYPROG
  PRINT*, "CREATED TO FORCE FCM TO LINK SOMETHING"
END PROGRAM DUMMYPROG
EOF
cd ..
build_compilation_script

# Run the compilation
./compilation.sh
ln -s build/bin/libphyex.so .

# Check if python can open the resulting shared lib
python3 -c "from ctypes import cdll; cdll.LoadLibrary('./libphyex.so')"

# ldd -r ./libphyex.so should also give interesting results
