#!/bin/bash

set -e
#set -x

. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"/../make_util.sh

ecbuild_version=tags/3.12.0

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

function setuparch() {
  ARCH=$1
  ARCHFILE=$2/arch.ecbuild
  ENVFILE=$2/arch.env
  if [ -f ${PHYEXCONF}/ecbuild_arch/arch-${ARCH}.ecbuild ]; then
    cp ${PHYEXCONF}/ecbuild_arch/arch-${ARCH}.env $ENVFILE
    cp ${PHYEXCONF}/ecbuild_arch/arch-${ARCH}.ecbuild $ARCHFILE
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

# ecbuild
export PATH=$PWD/../ecbuild/bin/:\$PATH
[ ! -f src/ecbuild ] && ln -s $PWD/../ecbuild src/

if [ \$(head -1 src/pyphyex.F90 | grep MODULE | wc -l) == 0 ]; then
  content=\$(cat src/pyphyex.F90)
  echo -e "MODULE PYPHYEX\nCONTAINS\n\${content}\nEND MODULE PYPHYEX" > src/pyphyex.F90
fi

BUILDDIR=$PWD/build
[ ! -d \$BUILDDIR ] && mkdir \$BUILDDIR
cd \$BUILDDIR

#fiat compilation
mkdir build_fiat
cd build_fiat
cmake ../../src/fiat -DCMAKE_INSTALL_PREFIX=\$BUILDDIR \\
                     -DCMAKE_CXX_FLAGS=-lstdc++ \\
                     -DCMAKE_BUILD_TYPE=\$level \\
                     -DBUILD_SHARED_LIBS=OFF \\
                     -DCMAKE_TOOLCHAIN_FILE=\$BUILDDIR/../arch.ecbuild
make -j
make install
cd ..

#PHYEX compilation
mkdir build_PHYEX
cd build_PHYEX
cmake ../../src -DCMAKE_INSTALL_PREFIX=\$BUILDDIR \\
                -DCMAKE_BUILD_TYPE=\$level \\
                -DCMAKE_TOOLCHAIN_FILE=\$BUILDDIR/../arch.ecbuild \\
                -DENABLE_PHYEX_BUILD_PROGS=ON \\
                -DENABLE_SINGLE_PRECISION=ON
make -j
make install
cd ..
EOF
chmod +x compilation.sh
}

####################################

# Do everything
EXTRASUBS="-s CMakeLists.txt -s cmake"
main ecbuild $*
