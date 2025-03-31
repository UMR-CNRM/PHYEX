#!/bin/bash

set -e
#set -x

. ../make_util.sh

ecbuild_version=tags/3.9.1

function check_install_build_system() {
  if [ ! -f ecbuild/bin/ecbuild ]; then
    echo "Performing ecbuild installation..."
    [ ! -d ecbuild ] && mkdir ecbuild
    cd ecbuild
    rm -f .gitkeep
    if [ $ssh -eq 1 ]; then
      git clone git@github.com:ecmwf/ecbuild.git .
    else
      git clone https://github.com/ecmwf/ecbuild .
    fi
    git checkout $ecbuild_version
    touch .gitkeep
    cd ..
    echo "...ecbuild installation done"
  fi
  cd ecbuild
  if [ $(git rev-parse HEAD^{commit}) != $(git rev-parse ${ecbuild_version}^{commit}) ]; then
    echo "Updating ecbuild installation..."
    cd ..
    echo "Cleaning ecbuild installation..."
    remove_build_system
    check_install_build_system
  else
    cd ..
  fi
}

function remove_build_system() {
  rm -rf ecbuild
}

function gmkfile2arch() {
  GMKFILE=$1
  ARCHFILE=$2/arch.ecbuild

  echo "
    Your are building an architecture file from a gmkfile.
    Please note that linker (command and and arguments) are ignored; ar command is also ignored"

cat <<EOF > $ARCHFILE
# Compilation
set(CMAKE_Fortran_COMPILER $(grep "^FRTNAME =" $GMKFILE | cut -d = -f 2))
set(CMAKE_Fortran_FLAGS "$(grep "^FRTFLAGS =" $GMKFILE | cut -d = -f 2-) $(grep "^GMK_FCFLAGS_PHYEX =" $GMKFILE | cut -d = -f 2-)")
set(CMAKE_Fortran_FLAGS_RELEASE "$(grep "^OPT_FRTFLAGS =" $GMKFILE | cut -d = -f 2-)")
set(CMAKE_Fortran_FLAGS_DEBUG "$(grep "^DBG_FRTFLAGS =" $GMKFILE | cut -d = -f 2-) $(grep "^BCD_FRTFLAGS =" $GMKFILE | cut -d = -f 2-) $(grep "^NAN_FRTFLAGS =" $GMKFILE | cut -d = -f 2-)")
set(CMAKE_C_COMPILER $(grep "^VCCNAME =" $GMKFILE | cut -d = -f 2))
set(CMAKE_C_FLAGS "$(grep "^VCCFLAGS =" $GMKFILE | cut -d = -f 2-)")
set(CMAKE_C_FLAGS_RELEASE "$(grep "^OPT_VCCFLAGS =" $GMKFILE | cut -d = -f 2-)")

# Preprocessor
add_compile_definitions("$<$<COMPILE_LANGUAGE:Fortran>:$(grep "^MACROS_FRT =" $GMKFILE | cut -d = -f 2-)>")
add_compile_definitions("$<$<COMPILE_LANGUAGE:C,CXX>:$(grep "^MACROS_CC =" $GMKFILE | cut -d = -f 2-)>")
EOF
}

function mesonhprofile2arch() {
  MESONHPROFILE=$1
  ARCHFILE=$2/arch.ecbuild
  ENVFILE=$2/arch.env

  echo "
   You are trying to produce a configuration file for fcm from a Meso-NH configuration.
   Please note that linker (command and and arguments) are ignored; ar command is also ignored
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
	echo "set(CMAKE_Fortran_COMPILER $(F90))"
	echo "set(CMAKE_Fortran_FLAGS \"-c $(F90FLAGS)\")"
	echo "set(CMAKE_Fortran_FLAGS_RELEASE \"\")"
	echo "set(CMAKE_Fortran_FLAGS_DEBUG \"\")"
	echo "set(CMAKE_C_COMPILER $(CC))"
	echo "set(CMAKE_C_FLAGS \"-c $(CFLAGS)\")"
	echo "set(CMAKE_C_FLAGS_RELEASE \"\")"
	echo "set(CMAKE_C_FLAGS_DEBUG \"\")"
	echo ""
	echo "# Preprocessor"
	echo "add_compile_definitions(\"$<$<COMPILE_LANGUAGE:Fortran>:$(CPPFLAGS)>\")"
	echo "add_compile_definitions(\"$<$<COMPILE_LANGUAGE:C,CXX>:$(CPPFLAGS)>\")"

'
  (set +e; . $MESONHPROFILE; set -e; make -f <(echo -e "$MAKEFILE") -s -I $(dirname $MESONHPROFILE)/../src archfile) > $ARCHFILE
}

function setuparch() {
  ARCH=$1
  ARCHFILE=$2/arch.ecbuild
  ENVFILE=$2/arch.env
  if [ -f ${HOME}/.phyex/ecbuild_arch/arch-${ARCH}.ecbuild ]; then
    cp ${HOME}/.phyex/ecbuild_arch/arch-${ARCH}.env $ENVFILE
    cp ${HOME}/.phyex/ecbuild_arch/arch-${ARCH}.ecbuild $ARCHFILE
  else 
    cp ${ARCH_PATH}/arch-${ARCH}.env $ENVFILE
    cp ${ARCH_PATH}/arch-${ARCH}.ecbuild $ARCHFILE
  fi
  touch $ENVFILE
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

level=Release #Release or Debug

export PATH=$PWD/../ecbuild/bin/:\$PATH

if [ \$(head -1 src/pyphyex.F90 | grep MODULE | wc -l) == 0 ]; then
  content=\$(cat src/pyphyex.F90)
  echo -e "MODULE PYPHYEX\nCONTAINS\n\${content}\nEND MODULE PYPHYEX" > src/pyphyex.F90
fi

BUILDDIR=$PWD/build
mkdir \$BUILDDIR
cd \$BUILDDIR

#fiat compilation
mkdir build_fiat
cd build_fiat
cmake ../../src/fiat -DCMAKE_INSTALL_PREFIX=\$BUILDDIR \
                     -DCMAKE_CXX_FLAGS=-lstdc++ \
                     -DCMAKE_BUILD_TYPE=\$level \
                     -DBUILD_SHARED_LIBS=OFF \
                     --toolchain \$BUILDDIR/../arch.ecbuild
make -j
make install
cd ..

#PHYEX common compilation
mkdir build_PHYEX
cd build_PHYEX
cmake ../../src -DCMAKE_INSTALL_PREFIX=\$BUILDDIR \
                -DFIAT_MODULE_DIR=\$BUILDDIR/module \
                -DFIAT_INCLUDE_DIR=\$BUILDDIR/include \
                -DFIAT_LIB_DIR=\$BUILDDIR/lib \
                -DCMAKE_BUILD_TYPE=\$level \
                --toolchain \$BUILDDIR/../arch.ecbuild
make -j VERBOSE=1
make install
cd ..
EOF
chmod +x compilation.sh
}

####################################

# Do everything
EXTRASUBS="-s CMakeLists.txt"
main ecbuild $*
