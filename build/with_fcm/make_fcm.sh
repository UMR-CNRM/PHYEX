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

function setuparch() {
  ARCH=$1
  ARCHFILE=$2/arch.fcm
  ENVFILE=$2/arch.env
  if [ -f ${PHYEXCONF}/fcm_arch/arch-${ARCH}.fcm ]; then
    cp ${PHYEXCONF}/fcm_arch/arch-${ARCH}.env $ENVFILE
    cp ${PHYEXCONF}/fcm_arch/arch-${ARCH}.fcm $ARCHFILE
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
