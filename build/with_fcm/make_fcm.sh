#!/bin/bash

set -e
#set -x

. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"/../make_util.sh

fcm_version=tags/2021.05.0

function check_install_build_system() {
  if [ ! -f fcm/bin/fcm ]; then
    echo "Performing FCM installation..."
    [ ! -d fcm ] && mkdir fcm
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
  cd fcm
  if [ $(git rev-parse HEAD^{commit}) != $(git rev-parse ${fcm_version}^{commit}) ]; then
    echo "Updating FCM installation..."
    cd ..
    echo "Cleaning FCM installation..."
    remove_build_system
    check_install_build_system
  else
    cd ..
  fi
}

function remove_build_system() {
  rm -rf fcm
}

function gmkfile2arch() {
  GMKFILE=$1
  ARCHFILE=$2/arch.fcm
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
cp fcm-make.cfg $builddir
}

function mesonhprofile2arch() {
  MESONHPROFILE=$1
  ARCHFILE=$2/arch.fcm
  ENVFILE=$2/arch.env

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
  (set +e; . $MESONHPROFILE; set -e; make -f <(echo -e "$MAKEFILE") -s -I $(dirname $MESONHPROFILE)/../src archfile) | sed 's/-D//g' > $ARCHFILE
  cp fcm-make.cfg $builddir
}

function setuparch() {
  ARCH=$1
  ARCHFILE=$2/arch.fcm
  ENVFILE=$2/arch.env
  if [ -f ${HOME}/.phyex/fcm_arch/arch-${ARCH}.fcm ]; then
    cp ${HOME}/.phyex/fcm_arch/arch-${ARCH}.env $ENVFILE
    cp ${HOME}/.phyex/fcm_arch/arch-${ARCH}.fcm $ARCHFILE
  else 
    cp ${ARCH_PATH}/arch-${ARCH}.env $ENVFILE
    cp ${ARCH_PATH}/arch-${ARCH}.fcm $ARCHFILE
  fi
  cp fcm-make.cfg $builddir
}

function build_compilation_script() {
srcdir=$1

#fcm doesn't like if a source directory doesn't exist.
#To be able to compile an old commit, we must filter the source directories
TESTPROGS_DIR=""
#support is not a testprog but is needed
for testprog in ice_adjust rain_ice turb_mnh shallow rain_ice_old support progs; do
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

ENTRYPOINTS="rain_ice.o shallow_mf.o turb.o ice_adjust.o lima.o lima_adjust_split.o pyphyex.o"

FCM_ARGS="$BS_ARGS"

echo "\\\$COMPIL_FFLAGS = \$COMPIL_FFLAGS" > config.fcm
echo "\\\$COMPIL_CFLAGS = \$COMPIL_CFLAGS" >> config.fcm
echo "\\\$LD_FLAGS = \$LD_FLAGS" >> config.fcm
echo "\\\$ENTRYPOINTS = \$ENTRYPOINTS" >> config.fcm
echo "\\\$LIBS = \$LIBS" >> config.fcm
echo "\\\$TESTPROGS_DIR=$TESTPROGS_DIR" >> config.fcm

export PATH=$PWD/../fcm/bin/:\$PATH

if [ \$(head -1 src/pyphyex.F90 | grep MODULE | wc -l) == 0 ]; then
  content=\$(cat src/pyphyex.F90)
  echo -e "MODULE PYPHYEX\nCONTAINS\n\${content}\nEND MODULE PYPHYEX" > src/pyphyex.F90
fi

echo "This script has generated config.fcm which is included by fcm-make.cfg, the FCM configuration file."
echo "Running : fcm make \$FCM_ARGS"

fcm make \$FCM_ARGS
[ ! -d build/lib/ ] && mkdir build/lib
mv build/bin/libphyex.so build/lib/
EOF
chmod +x compilation.sh
}

####################################

# Do everything
main fcm $*
