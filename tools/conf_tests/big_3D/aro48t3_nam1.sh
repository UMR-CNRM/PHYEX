#!/bin/bash
#SBATCH -p normal256
#SBATCH --export=MYLIB,HOMEPACK,TESTDIR
#SBATCH -n 1280
#SBATCH -c 4
#SBATCH -N 40
#SBATCH -t 00:40:00
#SBATCH --mem=247000
#SBATCH --exclusiv

# Job management :
# --------------
JOB_INITDIR=$SLURM_SUBMIT_DIR
export JOB_NAME=arome_e700
export JOB_ID=$SLURM_JOB_ID

echo JOB_INITDIR=$JOB_INITDIR
echo JOB_NAME=$JOB_NAME
echo JOB_ID=$JOB_ID

# =============================================================================

#                               RESOURCES ALLOCATIONS
#                               =====================

# Number of nodes/mpi-tasks/omp-threads:
# -------------------------------------
NNODES=$SLURM_JOB_NUM_NODES
# Number of MPI tasks per node:
MPITASKS_PER_NODE=$((SLURM_NTASKS/SLURM_JOB_NUM_NODES))
# Number of OPEN-MP threads per MPI task:
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Total number of MPI tasks:
MPI_TASKS=$SLURM_NTASKS
# Number of tasks reserved for the I/O server : 2 (hyperthreaded) nodes
NTASKS_IO=$(($(grep processor /proc/cpuinfo | wc -l)/1/$OMP_NUM_THREADS))

echo NNODES=$NNODES
echo MPITASKS_PER_NODE=$MPITASKS_PER_NODE
echo
# Number of MPI tasks and OMP threads used in the application :
echo MPI_TASKS=$MPI_TASKS
echo OMP_NUM_THREADS=$OMP_NUM_THREADS

# =============================================================================

#                               SYSTEM PREFERENCES
#                               ==================

# OMP/MPI submission management :
# -----------------------------
# LOCAL_MPI_WRAPPER : could be "mpiauto", "mpdrun", "mpiexec" ... or empty string
# LOCAL_STACK_LIMIT : could be "unlimited" or empty string

set -x
#LOCAL_MPI_WRAPPER="/opt/softs/mpiauto/mpiauto --wrap --wrap-stdeo --wrap-stdeo-pack"
LOCAL_MPI_WRAPPER="/opt/softs/mpiauto/mpiauto"
LOCAL_STACK_LIMIT=unlimited
ulimit -l unlimited
set +x

# Specific environment variables :
# ------------------------------
set -x
export OMP_STACKSIZE=4G
export KMP_STACKSIZE=4G
export KMP_MONITOR_STACKSIZE=4G
export I_MPI_HARD_FINALIZE=1
export I_MPI_SCALABLE_OPTIMIZATION=0
export I_MPI_DAPL_UD_RNDV_EP_NUM=4
export I_MPI_SHM_SPIN_COUNT=10
export I_MPI_SPIN_COUNT=10
set +x

# File systems :
# ------------
# Global file system:
export TMPGFS=$TMPDIR
export WORKGFS=$WORKDIR/benchmarks
#MTOOL export TMPGFS=$MTOOL_STEP_WORKSPACE
# Local file system (if preferred):
export TMPLOC=$TMPGFS
echo TMPGFS=$TMPGFS
echo TMPLOC=$TMPLOC

# Local disks synchronization :
# ---------------------------
export ISYNC=0
if [ "$MTOOL_IS" = "ON" ] ; then
# synchronization is needed anyway between the steps
  export ISYNC=1
elif [ $NNODES -gt 1 ] && [ "$TMPLOC" != "$TMPGFS" ] ; then
# Local disk synchronization needed:
  export ISYNC=1
fi
echo ISYNC=$ISYNC

# Miscellaneous :
# -------------

# =============================================================================

#                               USER PREFERENCES
#                               ================

#export NAMELDIR=/home/gmap/mrpm/khatib/pack/48t1_main.01#myref/run/cy47.forecast_arome_e700/Namelists
export NAMELDIR=$TESTDIR/Namelists

HOMEPACK=${HOMEPACK:=$HOME/pack}
export BINDIR=$HOMEPACK/$MYLIB/bin
#export BINDIR=/home/gmap/mrpm/khatib/pack/48t1_main.01#myref/bin
OUTPUTDIR=${OUTPUTDIR:-$PWD} #No cd command have been done before this line



export DATADIR=/scratch/work/khatib/data/cy47.forecast_arome_e700
#export REFDIR=/home/gmap/mrpm/khatib/benchmarks/apps/modules/cy47.forecast_arome_e700/References
export TOOLSDIR=/home/gmap/mrpm/khatib/benchmarks/tools
#export ROOTDIR_ODB=/home/gmap/mrpm/khatib/odbpools/36t1_bench/cy47.forecast_arome_e700

# Check reliability of auxilary directories :
# -----------------------------------------
ierr=0
#for var in NAMELDIR BINDIR DATADIR REFDIR TOOLSDIR ; do
for var in NAMELDIR BINDIR DATADIR TOOLSDIR ; do
  eval "dir=\$$var"
  if [ ! "$dir" ] ; then
    echo "$var is not set."
    ierr=1
  fi
  if [ $ierr -ne 0 ] ; then
    exit 1
  fi
done
ierr=0
for dir in $NAMELDIR $BINDIR $REFDIR $TOOLSDIR ; do
  if [ ! -d $dir ] ; then
    echo "$dir does not exists."
    ierr=1
  fi
  if [ $ierr -ne 0 ] ; then
    exit 1
  fi
done

echo TOOLSDIR=$TOOLSDIR
echo NAMELDIR=$NAMELDIR
echo DATADIR=$DATADIR
#echo REFDIR=$REFDIR
echo BINDIR=$BINDIR
#echo ROOTDIR_ODB=$ROOTDIR_ODB

export PATH=$TOOLSDIR:$PATH
export TOOLSDIR
export DATADIR

# Software default environment variables :
# --------------------------------------
set -x
export DR_HOOK=0
export DR_HOOK_IGNORE_SIGNALS=-1
export DR_HOOK_SILENT=1
export DR_HOOK_SHOW_PROCESS_OPTIONS=0
export MPL_MBX_SIZE=2048000000
export EC_PROFILE_HEAP=0
export EC_PROFILE_MEM=0
export EC_MPI_ATEXIT=0
export EC_MEMINFO=0
export OPENBLAS_NUM_THREADS=1
export MKL_CBWR="AUTO,STRICT"
export MKL_NUM_THREADS=1
export MKL_DEBUG_CPU_TYPE=5
set +x

# Profilers management :
# --------------------
# FTRACE_JOB : profiler switch
#              =0 : no profiler
#              =1 : integrated DrHook profiler
#              =2 : specific profiler
#
export FTRACE_JOB=1
echo "FTRACE_JOB=$FTRACE_JOB"

if [ $FTRACE_JOB -ne 0 ] ; then
# profilings main directory:
  if [ -d $JOB_INITDIR ] ; then
#   Use the initial job's dir
    FTRACE_DIR=$JOB_INITDIR
  else
#   Set one :
    FTRACEDIR=
    if [ ! "$FTRACEDIR" ] ; then
      echo "FTRACEDIR is not set."
      exit 1
    fi
    if [ ! -d $FTRACEDIR ] ; then
      mkdir -p $FTRACEDIR
      if [ $? -ne 0 ] ; then
        echo "Can't make directory $FTRACEDIR"
        exit 1
      fi
    fi
    FTRACE_DIR=$FTRACEDIR
  fi
   SCRATCH_FTRACE_DIR=$TMPGFS
#  SCRATCH_FTRACE_DIR=$FTRACE_DIR
  if [ $FTRACE_JOB -eq 1 ] ; then
    set -x
    export DR_HOOK=1
    export DR_HOOK_OPT=prof
#   Directory of individual profiles :
    export PROFDIR=$SCRATCH_FTRACE_DIR/${JOB_NAME}.d${JOB_ID}
#   Merged profiles report filename :
    export PROFMRG=$FTRACE_DIR/${JOB_NAME}.h${JOB_ID}
    set +x
  elif [ $FTRACE_JOB -ge 2 ] ; then
    set -x
#   Directory of individual profiles :
    export PROFDIR=$SCRATCH_FTRACE_DIR/${JOB_NAME}.f${JOB_ID}
#   Merged profiles report filename :
    export PROFMRG=$FTRACE_DIR/${JOB_NAME}.t${JOB_ID}
    set +x
  fi
fi

set +x

# Listings :
# --------
set -x
export ECHO_MPSH=OFF
export OUTPUT_LISTING=YES
export LOGDIR=$JOB_INITDIR/${JOB_NAME}.l${JOB_ID}
set +x

# ODB archives :
# ------------
# if set to 1, archived Odbs will be gzipped :
set -x
export ARCHIVE_AND_ZIP_ODB=0
set +x

# Directory for application output data files :
# -------------------------------------------
echo
OUTDIR=
OUTDIR=${OUTDIR:=$TMPGFS}
if [ "$TMPGFS" != "$TMPLOC" ] ; then
  if [ "$OUTDIR" = "$TMPLOC" ] ; then
    echo "Output files on LOCAL file system"
  elif [ "$OUTDIR" = "$TMPGFS" ] ; then
    echo "Output files on GLOBAL file system"
  else
    echo "Output files on directory : $OUTDIR"
  fi
else
  echo "Output files on directory : $OUTDIR"
fi

# NFS temporary directory for small I/Os
echo
TMPNFS=$(mktemp -d --tmpdir=/tmp/$LOGNAME)
if [ -d $TMPNFS ] ; then
  echo "temporary directory on NFS for small I/Os : $TMPNFS"
else
  TMPNFS="."
fi

# =============================================================================

#                               APPLICATION TUNING
#                               ==================

# ARPEGE : Forecast
# *****************

mkdir -p $TMPLOC
if [ $ISYNC -gt 0 ] ; then
  mkdir -p $TMPGFS
  cd $TMPGFS
else
  cd $TMPLOC
fi

# Driver-specific environment variables :
# -------------------------------------
set -x
NAMELIST=namel_previ.48
CTRLLIST=extra_namelists48.list
LINKS=links_inline48.scpt
EXECUTABLE=MASTERODB
#REFLIST=$REFDIR/forecast.out
EXPLIST=./NODE.001_01
set +x

#MTOOL common join=step_2

# Namelists modifications :
# -----------------------
 
set -x

# Number of MPI tasks for the I/O server :
NPROC_IO=$NTASKS_IO

# Remaining number of MPI tasks :
NPROC=$((MPI_TASKS-NPROC_IO))

# Memory cache optimisation:
NPROMA=-16
NFPROMA=-24

# Overall scalar optimisation:
LOPT_SCALAR=.TRUE.

# Output packing distribution:
NSTROUT=${NPROC}
NSTRIN=${NPROC}

NPRGPEW=16
#NPRGPNS=((NPROC/NPRGPEW))
NPRTRV=16
#NPRTRW=((NPROC/NPRTRV))

set +x

cat > namelist_mods2 <<EOF
 &NAM_PARAM_ICE
   CFRAC_ICE_ADJUST='S',
   CFRAC_ICE_SHALLOW_MF='S',
   CSEDIM='STAT',
   CSNOWRIMING='M90',
   LCONVHG=.TRUE.,
   LCRFLIMIT=.TRUE.,
   LEVLIMIT=.TRUE.,
   LFEEDBACKT=.TRUE.,
   LNULLWETG=.TRUE.,
   LNULLWETH=.TRUE.,
   LSEDIM_AFTER=.FALSE.,
   LWETGPOST=.TRUE.,
   LWETHPOST=.TRUE.,
   NMAXITER_MICRO=1,
   XFRACM90=0.1,
   XMRSTEP=0.00005,
   XSPLIT_MAXCFL=0.8,
   XTSTEP_TS=0.,
   LCRIAUTI=.TRUE.,
   XCRIAUTC_NAM=0.001,
   XCRIAUTI_NAM=0.0002,
   XT0CRIAUTI_NAM=-5.,
   LRED=.TRUE.,
   LSEDIC=.TRUE.,
 /
 &NAMPARAR
   LOSEDIC=-,
   CFRAC_ICE_ADJUST=-,
   CFRAC_ICE_SHALLOW_MF=-,
   CSEDIM=-,
   CSNOWRIMING=-,
   LCONVHG=-,
   LCRFLIMIT=-,
   LEVLIMIT=-,
   LFEEDBACKT=-,
   LNULLWETG=-,
   LNULLWETH=-,
   LSEDIM_AFTER=-,
   LWETGPOST=-,
   LWETHPOST=-,
   NMAXITER_MICRO=-,
   XFRACM90=-,
   XMRSTEP=-,
   XSPLIT_MAXCFL=-,
   XTSTEP_TS=-,
   LCRIAUTI=-,
   RCRIAUTC=-,
   RCRIAUTI=-,
   RT0CRIAUTI=-,
 /
 &NAMTRANS
   LFFTW=.TRUE.,
 /
 &NAMPAR0
   NPRINTLEV=1,
   LOPT_SCALAR=${LOPT_SCALAR},
   MBX_SIZE=2048000000,
   NPROC=${NPROC},
   NPRGPNS=-,
   NPRGPEW=-,
   NPRTRW=-,
   NPRTRV=-,
 /
 &NAMDIM
   NPROMA=$NPROMA,
 /
 &NAMFPSC2
   NFPROMA=$NFPROMA,
 /
 &NAMFPSC2_DEP
   NFPROMA_DEP=$NFPROMA,
 /
 &NAMPAR1
   LSPLIT=.TRUE.,
   NSTRIN=${NSTRIN},
   NSTROUT=${NSTROUT},
 /
 &NAMFA
   CMODEL=' ',
 /
 &NAMIAU
   LIAU=.FALSE.,
 /
 &NAMARG
   CNMEXP='0000',
 /
 &NAMCT0
   CSCRIPT_LAMRTC=' ',
   CSCRIPT_PPSERVER=' ',
   CFPNCF='ECHFP',
   NSDITS(0)=0,
   NFRSDI=4,
   NFPOS=1,
 /
 &NAMCT1
   N1POS=1,
 /
 &NAMFPC
   CFPDIR='${OUTDIR}/PF',
 /
 &NAMOPH
   CFNHWF='${OUTDIR}/ECHIS',
   CFPATH='${OUTDIR}/',
 /
 &NAMIO_SERV 
   NPROC_IO=${NPROC_IO}, 
   NMSG_LEVEL_SERVER=1, 
   NMSG_LEVEL_CLIENT=1, 
   NPROCESS_LEVEL=5,
 /
 &NAMRIP
    CSTOP='h24',
    TSTEP=50.,
 /
EOF
cat namelist_mods2 > namelist_modset
\rm -f namelist_mods2
echo
echo Namelists adaptations :
cat namelist_modset
echo

set +x
cp $NAMELDIR/$NAMELIST namelist
perl -w $TOOLSDIR/xpnam namelist --dfile=namelist_modset
set -x
echo
/bin/cat namelist.new
set +x
\rm -f namelist_modset namelist
\mv namelist.new fort.4
set -x

#MTOOL common

# =============================================================================

#                               DRIVER
#                               ======

#      ******************************
#      *  fetch initial data files  *
#      ******************************

#MTOOL common join=step_1

set -x
$TOOLSDIR/getdata.sh
set +x

#MTOOL common

#MTOOL common join=step_2

for file in $(cat $NAMELDIR/$CTRLLIST) ; do
  set -x
  cp $NAMELDIR/$file .
  set +x
done
if [ -s $NAMELDIR/$LINKS ] ; then
  set -x
  cp $NAMELDIR/$LINKS .
  chmod 755 $LINKS
  . ./$LINKS
  \rm $LINKS
  set +x
fi

#MTOOL common

#      ***************
#      *  Executable *
#      ***************

#MTOOL common join=step_1

echo
set -x
cp $BINDIR/$EXECUTABLE .
set +x
if [ ! -f $EXECUTABLE ] ; then
  echo "executable $BINDIR/$EXECUTABLE could not be copied."
  exit 1
fi

#MTOOL common

#      ********************************
#      *  Prepare parallel executions *
#      ********************************

#MTOOL common join=step_2

if [ "$LOCAL_STACK_LIMIT" ] ; then
  set -x
  ulimit -s $LOCAL_STACK_LIMIT
  set +x
fi
# for mpsh :
export MPSH_NPES=$NNODES

# grib_api environment variables may be determined by the executable : 
. grib_api_profile $EXECUTABLE

# Intel mpi fabric setup depending on what is found in the executable :
. intel_mpi_fabric $EXECUTABLE

#MTOOL common

set -x
cd $TMPLOC
set +x
. rttov_profile

#      *******************************************************
#      *  Unarchive datasets and local disks Synchronisation *
#      *******************************************************

if [ $ISYNC -eq 0 ] ; then
  set -x
#MTOOL common join=step_1
  $TOOLSDIR/input_sync.sh
#MTOOL common
  set +x
else
  set -x
#MTOOL common join=step_2
  $TOOLSDIR/input_sync.sh
#MTOOL common
  set +x
fi

#      ***************
#      *  Execution  *
#      ***************

#MTOOL common join=step_2

mkdir -p $OUTDIR
echo
if [ $(echo $LOCAL_MPI_WRAPPER | grep -c mpiauto) -ne 0 ] ; then
  set -x
  time $LOCAL_MPI_WRAPPER -np $MPI_TASKS -nnp $MPITASKS_PER_NODE -- ./$EXECUTABLE </dev/null \
  errorcode=$?
  2>&1 | grep -v "FA[DC]GR[AM]: Field .* is not declared in \`faFieldName.def'"
  set +x
elif [ "$LOCAL_MPI_WRAPPER" = "srun" ] ; then
  set -x
  time $LOCAL_MPI_WRAPPER ./$EXECUTABLE </dev/null \
  errorcode=$?
  2>&1 | grep -v "FA[DC]GR[AM]: Field .* is not declared in \`faFieldName.def'"
  set +x
elif [ "$LOCAL_MPI_WRAPPER" ] ; then
  set -x
  time $LOCAL_MPI_WRAPPER -np $MPI_TASKS ./$EXECUTABLE </dev/null \
  errorcode=$?
  2>&1 | grep -v "FA[DC]GR[AM]: Field .* is not declared in \`faFieldName.def'"
  set +x
else
  set -x
  time ./$EXECUTABLE \
  errorcode=$?
  2>&1 | grep -v "FA[DC]GR[AM]: Field .* is not declared in \`faFieldName.def'"
  set +x
fi

#      **********************
#      *  Post-processings  *
#      **********************

echo
if [ "$OUTPUT_LISTING" = "YES" ] ; then
  set -x
  $TOOLSDIR/outsync.sh
  set +x
fi

if [ $FTRACE_JOB -gt 0 ] ; then
  set -x
  $TOOLSDIR/profsync.sh
  set +x
fi

set -x
ls -l $OUTDIR
set +x

set -x
#errorcode returned by executable is not reliable (always different from 0)
if grep " NSTEP =  1728 CNT0" NODE.001_01 > /dev/null; then
  cp $EXPLIST $OUTPUTDIR/
else
  mkdir $OUTPUTDIR/error
  cp $EXPLIST $OUTPUTDIR/error/
fi
#if [ -f $REFLIST ] && [ -f $EXPLIST ] ; then $TOOLSDIR/diffNODE.001_01 $EXPLIST $REFLIST ; fi
set +x
#      ****************
#      *  Cleanups    *
#      ****************

set -x
cd $TMPGFS
$TOOLSDIR/cleansync.sh
set +x

#MTOOL common

#      ****************
#      *  Epilogue    *
#       ****************

set -x
$TOOLSDIR/epilog.sh
set +x
if [ "$MTOOL_IS" != "ON" ] && [ "$AUTO_CLEAN" = "ON" ] ; then
  cd $HOME
  \rm -rf $TMPGFS
fi

#MTOOL step id=step_1 target=FRONTEND
#MTOOL step id=step_2 target=SUPERCOMPUTER
