#!/usr/bin/env bash
#
# Build the pyphyex wheel INSIDE a manylinux container (see the release
# workflow), with the repository bind-mounted at /src. The wheel is written in
# /src/wheelhouse, repaired by auditwheel: the Fortran and OpenMP runtime
# libraries are bundled.
#
# Usage: build_wheel.sh HOST_UID HOST_GID
#   HOST_UID, HOST_GID: user and group of the runner (the container runs as
#   root: the outputs are given back to the runner user)

set -e

HOST_UID="$1"
HOST_GID="$2"
WHEELHOUSE=/src/wheelhouse

cd /src
trap 'chown -R "${HOST_UID}:${HOST_GID}" "${WHEELHOUSE}" 2>/dev/null || true' EXIT

# The repository belongs to the runner user, not to root
git config --global --add safe.directory '*'

if [ -x /usr/bin/dnf ]; then
    PKG=dnf   # AlmaLinux 8 based (manylinux_2_28, glibc 2.28)
else
    PKG=yum   # CentOS 7 based (manylinux2014, glibc 2.17)
    # CentOS 7 is EOL: the SCLo repositories (devtoolset packages) still point
    # at the dead mirrorlist.centos.org. Disable them and provide vault-based
    # repositories.
    enable_vault_sclo() {
        yum -y install centos-release-scl 2>/dev/null || true
        sed -i "s/^enabled=1/enabled=0/" /etc/yum.repos.d/CentOS-SCLo-*.repo 2>/dev/null || true
        printf "%s\n" \
            "[vault-sclo-rh]" \
            "name=Vault CentOS-7 SCLo rh (glibc 2.17)" \
            "baseurl=http://vault.centos.org/7.9.2009/sclo/\$basearch/rh/" \
            "enabled=1" \
            "gpgcheck=0" \
            > /etc/yum.repos.d/vault-sclo-rh.repo
    }
fi

# Python (>= 3.10, needed by pyphyex): the manylinux images ship CPythons in
# /opt/python, enough to build a wheel. Otherwise, install the one of the
# distribution (only possible on AlmaLinux 8: CentOS 7 provides python 3.8 at
# most, through SCL).
PYEXE=""
for py in /opt/python/cp312-cp312/bin/python /opt/python/cp311-cp311/bin/python \
          /opt/python/cp310-cp310/bin/python; do
    if [ -x "${py}" ]; then
        PYEXE="${py}"
        break
    fi
done
if [ -z "${PYEXE}" ]; then
    if [ "${PKG}" = dnf ]; then
        dnf -y install python3.11 python3.11-pip
        PYEXE=/usr/bin/python3.11
    else
        echo "No python >= 3.10 found in /opt/python, and none can be installed on CentOS 7"
        exit 1
    fi
fi
"${PYEXE}" --version

# Build tools (cmake and ninja, also used by the build of fiat), auditwheel
"${PYEXE}" -m venv /tmp/buildtools
/tmp/buildtools/bin/pip install --upgrade pip
/tmp/buildtools/bin/pip install cmake ninja
export PATH="/tmp/buildtools/bin:${PATH}"
command -v auditwheel > /dev/null || /tmp/buildtools/bin/pip install auditwheel

# Compilers: gcc of the image (a gcc-toolset/devtoolset), with its gfortran
if ! command -v gfortran > /dev/null; then
    gccdir="$(dirname "$(command -v gcc)")"
    if [[ "${gccdir}" =~ ^/opt/rh/([^/]+)/ ]]; then
        fortran_pkg="${BASH_REMATCH[1]}-gcc-gfortran"
    else
        fortran_pkg=gcc-gfortran
    fi
    [ "${PKG}" = yum ] && enable_vault_sclo
    ${PKG} -y install "${fortran_pkg}"
fi
export CC=gcc CXX=g++ FC=gfortran
gcc --version | head -1
gfortran --version | head -1

# Wheel (ecbuild and fiat are downloaded and built by the top-level
# CMakeLists.txt; pyfortool is installed in the build environment), then
# bundling of the external libraries by auditwheel
rm -rf /tmp/rawwheel "${WHEELHOUSE}"
"${PYEXE}" -m pip wheel --no-deps -v -w /tmp/rawwheel .
auditwheel show /tmp/rawwheel/*.whl
auditwheel repair -w "${WHEELHOUSE}" /tmp/rawwheel/*.whl
ls -l "${WHEELHOUSE}"
