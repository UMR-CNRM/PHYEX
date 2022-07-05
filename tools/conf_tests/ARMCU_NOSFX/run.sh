#!/bin/sh

set -ex

export OMP_NUM_THREADS=1
export DR_HOOK_IGNORE_SIGNALS=-1
export DR_HOOK=0
ulimit -s unlimited
unset LD_LIBRARY_PATH

tar xf rrtm.tgz

./MASTER >lola 2>&1

