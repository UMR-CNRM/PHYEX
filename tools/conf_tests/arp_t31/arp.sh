#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH --export=MYLIB,HOME,HOMEPACK,TMPDIR,OUTPUTDIR,TESTDIR
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -p normal256

#The MYLIB variable must contain the gmkpack pack name
#The TESTDIR variable must contain the test directory
#Results will be stored in the local directory

#Other environment variables that can be set:
#OUTPUTDIR

date
set -e

case=t31

OUTPUTDIR=${OUTPUTDIR:-$PWD}

NPROC=4
NSTRIN=$NPROC
NSTROUT=12
NPRTRW_NPRTRV="  NPRTRW=$NPROC,
  NPRTRV=1,
  NPRGPNS=$NPROC,
  NPRGPEW=1,"
export OMP_NUM_THREADS=1

export DR_HOOK=1
#export DR_HOOK_IGNORE_SIGNALS=-1
export DR_HOOK_NOT_MPI=1
export DR_HOOK_SILENT=1
export DR_HOOK_OPT=prof

export EC_PROFILE_HEAP=0
export EC_PROFILE_MEM=0
export EC_MPI_ATEXIT=0
export DR_HOOK_SHOW_PROCESS_OPTIONS=0
export EC_MEMINFO=0
export TVSEARCHPATH=$SOURCE

HOMEPACK=${HOMEPACK:=$HOME/pack}
SOURCE=$HOMEPACK/$MYLIB/src/local
LOADIR=$HOMEPACK/$MYLIB/bin

TMPDIR=${TMPDIR:=$HOME/tmp}
TMPLOC=$TMPDIR/rundir.$$
TMPWAIT=$TMPDIR/wait_queue.$$
mkdir $TMPWAIT
mkdir $TMPLOC
cd $TMPLOC

export RTTOV_COEFDIR=$PWD

#      **************************
#      *  Saisie des NAMELISTS  *
#      **************************

CNMEXP='FCST'

echo
/bin/cat <<FIN > fort.4
 &NACIETEO
 /
 &NACOBS
 /
 &NACTAN
 /
 &NACTEX
 /
 &NACVEG
 /
 &NADOCK
 /
 &NAEAEM7
 /
 &NAEAER
 /
 &NAECOAPHY
 /
 &NAEPHLI
 /
 &NAEPHY
 /
 &NAERAD
   LRRTM=.TRUE.,
   LSRTM=.FALSE.,
   NMCICA=0,
   NRADFR=-1,
   NSW=6,
   RLWINHF=0.9,
 /
 &NAERCLI
 /
 &NAETLDIAG
 /
 &NAEVOL
 /
 &NAIMPO
 /
 &NALORI
 /
 &NAMACV
 /
 &NAMAERDET
 /
 &NAMAFN
 /
 &NAMARG
   CNMEXP='${CNMEXP}',
   NCONF=1,
   LELAM=.FALSE.,
   LECMWF=.FALSE.,
   NSUPERSEDE=1,
 /
 &NAMARPHY
 /
 &NAMCA
 /
 &NAMCAPE
 /
 &NAMCFU
 /
 &NAMCHEM
 /
 &NAMCHET
 /
 &NAMCHK
 /
 &NAMCLA
 /
 &NAMCLDP
 /
 &NAMCLI
 /
 &NAMCLRADLID
 /
 &NAMCLTC
 /
 &NAMCOK
 /
 &NAMCOM
 /
 &NAMCOMPO
 /
 &NAMCOSJO
 /
 &NAMCT0
   LFBDAP=.TRUE.,
   LFDBOP=.FALSE.,
   NFRHIS=1,
   NFRPOS=10000,
   NFRSDI=1,
   NPOSTS(0)=1,
   NPOSTS(1)=1,
   NSPPR=0,
   LALLOPR=.FALSE.,
   LGRIB_API=.FALSE.,
 /
 &NAMCT1
   LRFILAF=.FALSE.,
   N1RES=0,
 /
 &NAMCUMF
 /
 &NAMCUMFS
 /
 &NAMCVER
   LVERTFE=.TRUE.,
   CVFE_ETAH='CHORDAL',
   NVFE_TYPE=3,
   LREGETA=.FALSE.,
   LAPRXPK=.TRUE.,
 /
 &NAMCVMNH
   LSMOOTH=.FALSE.,
   OTADJS=10800.,
   XATPERT=500.,
   XBTPERT=0.5,
   XAW=1.,
   XBW=0.5,
   XCDEPTH=1.,
   XCDEPTH_D=4000.,
   XDTPERT=0.3,
   XENTR=0.01,
 /
 &NAMDDH
 /
 &NAMDFI
 /
 &NAMDIM
   NPROMA=-16,
 /
 &NAMDIMO
 /
 &NAMDIM_TRAJ
 /
 &NAMDPHY
 /
 &NAMDPRECIPS
 /
 &NAMDVISI
 /
 &NAMDWET
 /
 &NAMDYN
   BETADT=1.,
   LADVF=.TRUE.,
   NTLAG=3,
   NVLAG=3,
   NWLAG=3,
   RCMSLP0=1.,
   SIPR=100000.,
   SITR=350.,
   VESL=0.,
   XIDT=0.,
   NITMP=5,
   RW2TLFF=1.,
   RATIO_HDI_TOP=20.,
   SLEVDH1=0.15,
   SLEVDH2=0.08,
   RDAMPDIV=4.,
   RDAMPVOR=50.,
   RDAMPT=0.,
   RDAMPQ=0.,
 /
 &NAMDYNA
   LNESC=.FALSE.,
   LNESCT=.FALSE.,
   LNESCV=.FALSE.,
   LSETTLS=.TRUE.,
   LSETTLST=.TRUE.,
   LSETTLSV=.TRUE.,
   LSPRT=.TRUE.,
   LTWOTL=.TRUE.,
   LSLAG=.TRUE.,
 /
 &NAMDYNA_STATIC
   LDRY_ECMWF=.FALSE.,
 /
 &NAMDYNCORE
 /
 &NAMECV
 /
 &NAMECVDESC
 /
 &NAMECVGRB
 /
 &NAMEMIS_CONF
 /
 &NAMENKF
 /
 &NAMENSCOV
 /
 &NAMFA
   NBITCS=30,
   NBITPG=30,
   NSTRON=-1,
   YFAL%NBITS=16,
   YFAI%NBITS=16,
   YFAR%NBITS=16,
   YFAS%NBITS=16,
   YFALRAD%NBITS=16,
   YFAIRAD%NBITS=16,
   YFACLF%NBITS=6,
   YFATKE%NBITS=16,
   NVGRIB=123,
 /
 &NAMFAINIT
 /
 &NAMFPC
 /
 &NAMFPD
 /
 &NAMFPDY2
 /
 &NAMFPDYF
 /
 &NAMFPDYH
 /
 &NAMFPDYI
 /
 &NAMFPDYP
 /
 &NAMFPDYS
 /
 &NAMFPDYT
 /
 &NAMFPDYV
 /
 &NAMFPF
 /
 &NAMFPG
 /
 &NAMFPIOS
 /
 &NAMFPMOVE
 /
 &NAMFPOBJ
 /
 &NAMFPPHY
 /
 &NAMFPSC2
 /
 &NAMFPSC2_DEP
 /
 &NAMGEM
   LNONHYD_GEOM=.FALSE.,
   LNHX_GEOM=.FALSE.,
 /
 &NAMGFL
   YI_NL%LGPINGP=.TRUE.,
   YI_NL%LGP=.TRUE.,
   YI_NL%LT1=.TRUE.,
   YI_NL%LPHY=.FALSE.,
   YI_NL%NREQIN=1,
   YI_NL%LREQOUT=.FALSE.,
   YI_NL%LADV=.TRUE.,
   YI_NL%LQM=.TRUE.,
   YL_NL%LGPINGP=.TRUE.,
   YL_NL%LGP=.TRUE.,
   YL_NL%LT1=.TRUE.,
   YL_NL%LPHY=.FALSE.,
   YL_NL%NREQIN=1,
   YL_NL%LREQOUT=.FALSE.,
   YL_NL%LADV=.TRUE.,
   YL_NL%LQM=.TRUE.,
   YR_NL%LGPINGP=.TRUE.,
   YR_NL%LGP=.TRUE.,
   YR_NL%LT1=.TRUE.,
   YR_NL%LPHY=.FALSE.,
   YR_NL%NREQIN=1,
   YR_NL%LREQOUT=.FALSE.,
   YR_NL%LADV=.TRUE.,
   YR_NL%LQM=.TRUE.,
   YS_NL%LGPINGP=.TRUE.,
   YS_NL%LGP=.TRUE.,
   YS_NL%LT1=.TRUE.,
   YS_NL%LPHY=.FALSE.,
   YS_NL%NREQIN=1,
   YS_NL%LREQOUT=.FALSE.,
   YS_NL%LADV=.TRUE.,
   YS_NL%LQM=.TRUE.,
   YTKE_NL%LGPINGP=.TRUE.,
   YTKE_NL%LGP=.TRUE.,
   YTKE_NL%LT1=.TRUE.,
   YTKE_NL%NREQIN=1,
   YTKE_NL%LREQOUT=.TRUE.,
   YTKE_NL%LADV=.TRUE.,
   YTKE_NL%LQM=.TRUE.,
   YIRAD_NL%LGP=.TRUE.,
   YIRAD_NL%NREQIN=0,
   YIRAD_NL%LREQOUT=.TRUE.,
   YLRAD_NL%LGP=.TRUE.,
   YLRAD_NL%NREQIN=0,
   YLRAD_NL%LREQOUT=.TRUE.,
   YA_NL%LGP=.TRUE.,
   YA_NL%NREQIN=0,
   YA_NL%LREQOUT=.TRUE.,
 /
 &NAMGRIB
 /
 &NAMGWD
 /
 &NAMGWDIAG
 /
 &NAMGWWMS
 /
 &NAMIAU
 /
 &NAMICE
 /
 &NAMINI
   LDFI=.FALSE.,
 /
 &NAMINTFLEX
 /
 &NAMIOMI
 /
 &NAMIOS
 /
 &NAMIO_SERV
 /
 &NAMJBALPHACV
 /
 &NAMJBCODES
 /
 &NAMJBECPHYSECV
 /
 &NAMJBHYBACV
 /
 &NAMJBSKTECV
 /
 &NAMJG
 /
 &NAMLCZ
 /
 &NAMLSFORC
 /
 &NAMMARS
 /
 &NAMMCC
 /
 &NAMMCUF
 /
 &NAMMETHOX
 /
 &NAMMKODB
 /
 &NAMMODERR
 /
 &NAMMODERRCONF
 /
 &NAMMODERRCOV
 /
 &NAMMODERRINCRCONF
 /
 &NAMMODERRMOD
 /
 &NAMMTS
 /
 &NAMMWAVE
 /
 &NAMNORGWD
 /
 &NAMNPROF
 /
 &NAMNRTAER
 /
 &NAMNUD
 /
 &NAMNUDGLH
 /
 &NAMOBS
 /
 &NAMOOPS
 /
 &NAMOPH
   LINC=.TRUE.,
 /
 &NAMOPTCMEM
 /
 &NAMPAR0
   MBX_SIZE=2048000000,
   NOUTPUT=1,
   NPROC=$NPROC,
   MP_TYPE=2,
   NPRINTLEV=1,
 /
 &NAMPAR1
   NSTRIN=$NSTRIN,
   NSTROUT=$NSTROUT,
   NCOMBFLEN=1638400,
   LEQ_REGIONS=.FALSE.,
 /
 &NAMPARAR
 /
 &NAMPARECV
 /
 &NAMPERTOBS
 /
 &NAMPERTPAR
 /
 &NAMPHMSE
 /
 &NAMPHY
   CGMIXLEN='AY',
   LAERODES=.TRUE.,
   LAEROLAN=.TRUE.,
   LAEROSEA=.TRUE.,
   LAEROSOO=.TRUE.,
   LCONDWT=.TRUE.,
   LDIFCONS=.TRUE.,
   LFPCOR=.TRUE.,
   LNEWD=.TRUE.,
   LNOIAS=.TRUE.,
   LO3ABC=.TRUE.,
   LPROCLD=.TRUE.,
   LRAY=.FALSE.,
   LRAYFM=.TRUE.,
   LRAYLU=.TRUE.,
   LRNUMX=.TRUE.,
   LSSD=.TRUE.,
   LSTRA=.FALSE.,
   LVGSN=.TRUE.,
   LCVPPKF=.TRUE.,
   LECDEEP=.TRUE.,
   LECSHAL=.TRUE.,
   LECT=.TRUE.,
   LFLUSO=.TRUE.,
   LNEBECT=.FALSE.,
   LO3FL=.TRUE.,
   LECTFL=.TRUE.,
   LECTFL0=.TRUE.,
   LZ0HSREL=.TRUE.,
   LADJCLD=.TRUE.,
   LSMITH_CDEV=.TRUE.,
   NCALLRAD=2,
   LGLACIERS=.TRUE.,
   NDPSFI=0,
   LNEBN=.FALSE.,
 /
 &NAMPHY0
   EDD=1.,
   EDK=1.,
   GCVNU=0.00005,
   GCVPSI=1.,
   GCVPSIE=1.,
   GWDCD=5.4,
   GWDSE=0.005,
   GWDVALI=0.5,
   QSNEBC=-1.,
   QSSUSC=5.,
   RCVEVAP=0.25,
   REVASX=2.0D-7,
   RICRLM=0.5,
   TDDGP=0.6,
   TUDGP=0.6,
   USURIC=0.175,
   USURICL=1.,
   USURID=0.1,
   VZ0CM=0.0001,
   XBLM=8.5,
   XMAXLM=5000.,
   XMINLM=10.,
   ALMAVX=1000.,
   GCVHMIN=30000.,
   RFACNSM=1.2,
   RFLCHCE=0.25,
   RKFBTAU=7200.,
   RPRTH=1.,
   RQICRT2=0.,
   RQICVMIN=0.00001,
   SXNBCO=1.,
   ALMAVE=0.,
   RQCRNS=0.,
   RQICRSN=1.,
   TFVI=0.08,
   TFVL=0.02,
   TFVS=1.5,
   GCVOMGE=1.,
   ECTMAX=35.,
   AECLS4=3.,
   GDDEVA=0.2,
   GCVOMGQ=0.,
   REFLKUO=5000.,
   GCVRHMIN=0.5,
   GCVOMDPS=2000.,
   GCVOMCA=-1.918,
   GCVRHMAX=1.,
   GCVOMGS=1.,
   GCVOMCC=1.35,
   GCVOMGSX=0.5,
 /
 &NAMPHY1
   ALBMIN=0.65,
   ALCRIN=0.81,
   EMMGLA=0.98,
   EMMMER=0.99,
   ALRCN1=0.01,
 /
 &NAMPHY2
   FACRAF=4.,
   XMULAF=0.,
   LRAFTKE=.TRUE.,
   HTKERAF=20.,
 /
 &NAMPHY3
 /
 &NAMPHYDS
 /
 &NAMPPC
 /
 &NAMPPVI
 /
 &NAMPRE
 /
 &NAMRADCMEM
 /
 &NAMRCF
 /
 &NAMRCOEF
 /
 &NAMRES
 /
 &NAMRGRI
 /
 &NAMRINC
 /
 &NAMRIP
   CSTOP='h6',
   TSTEP=600,
 /
 &NAMRIP0
 /
 &NAMRLX
 /
 &NAMRPP
 /
 &NAMRSTRHBIAS
 /
 &NAMSATS
 /
 &NAMSATSIM
 /
 &NAMSCC
 /
 &NAMSCEN
 /
 &NAMSCM
 /
 &NAMSEKF
 /
 &NAMSENS
 /
 &NAMSIMPHL
 /
 &NAMSPNG
 /
 &NAMSPP
 /
 &NAMSPSDT
 /
 &NAMSTA
 /
 &NAMSTOPH
 /
 &NAMSWE
 /
 &NAMTESTVAR
 /
 &NAMTHLIM
 /
 &NAMTOPH
   ETCVIM=5000.,
   ETNEBU=5000.,
   ETPLUI=5000.,
   XDRMTK=2.0D-8,
   XDRMTP=800.,
   XDRMTX=4.0D-7,
   XDRMUK=1.0D-7,
   XDRMUP=800.,
   XDRMUX=0.000002,
 /
 &NAMTRAJ
 /
 &NAMTRAJP
 /
 &NAMTRANS
   LUSEFLT=.FALSE.,
   LUSERPNM=.TRUE.,
   LKEEPRPNM=.TRUE.,
   LFFTW=.TRUE.,
 /
 &NAMTRANS0
 /
 &NAMVAR
 /
 &NAMVARBC
 /
 &NAMVARBC_AIREP
 /
 &NAMVARBC_ALLSKY
 /
 &NAMVARBC_GBRAD
 /
 &NAMVARBC_MODES
 /
 &NAMVARBC_RAD
 /
 &NAMVARBC_SFCOBS
 /
 &NAMVARBC_TCWV
 /
 &NAMVARBC_TO3
 /
 &NAMVAREPS
 /
 &NAMVDF
 /
 &NAMVDOZ
 /
 &NAMVOLCANO
 /
 &NAMVRTL
 /
 &NAMVV0
 /
 &NAMVV1
 /
 &NAMVWRK
 /
 &NAMWAVELETJB
 /
 &NAMXFU
 /
 &NAM_CANAPE
 /
 &NAM_DISTRIBUTED_VECTORS
 /
 &NAM_NEBN
 /
 &NAM_PARAM_ICEN
 /
 &NAM_PARAM_LIMA
 /
 &NAM_PARAM_MFSHALLN
 /
 &NAM_SPP_ACTIVE
 /
 &NAM_SPP_MODIF
 /
 &NAM_TURBN
 /
 &NAPHLC
 /
 &NEMCT0
 /
 &NEMDIM
 /
 &NEMDYN
 /
 &NEMELBC0A
 /
 &NEMELBC0B
 /
 &NEMFPEZO
 /
 &NEMGEO
 /
 &NEMJK
 /
 &NEMVAR
 /
 &NEMWAVELET
 /
FIN
/bin/cat fort.4

#      *****************************************
#      *  Acquisition du fichier de demarrage  *
#      *****************************************

echo
set -x
ln -s $TESTDIR/data/arp/$case/ICMSHFCSTINIT ICMSH${CNMEXP}INIT
set +x
tar xfz $TESTDIR/data/rtm/rrtm.const.04.tgz

#      ***************
#      *  Chargement *
#      ***************

echo
set -x
\ln -s $LOADIR/MASTERODB MASTER
set +x
if ldd MASTER | grep openmpi > /dev/null; then
  #On est sur PC
  set +e
  MPILIB=$(ldd MASTER | grep openmpi | tail -1 | awk '{print $3}' | awk -F "/" '{print $(NF-2)}')
  MPIRUN="$(echo $(dirname $(dirname $(ldd MASTER | grep openmpi| tail -1 | awk '{print $3}'))))/bin/mpirun --oversubscribe -np $NPROC"
  GRIB_API=$(dirname $(dirname $(ldd MASTER | grep grib_api | head -1 | awk '{print $3}') 2>/dev/null) 2>/dev/null)
  ECCODES=$(dirname $(dirname $(ldd MASTER | grep eccodes | head -1 | awk '{print $3}') 2>/dev/null) 2>/dev/null)
  set -e
  export GRIB_SAMPLES_PATH=$GRIB_API/share/grib_api/ifs_samples/grib1
  export GRIB_DEFINITION_PATH=$GRIB_API/share/grib_api/definitions
  export ECCODES_SAMPLES_PATH=$ECCODES/share/eccodes/ifs_samples/grib1
  export ECCODES_DEFINITION_PATH=$TESTDIR/data/eccodes_extras_definitions:$ECCODES/share/eccodes/definitions
else
  #On est sur HPC
  #MPIRUN="$(echo $(dirname $(dirname $(ldd MASTER | grep libmpi| tail -1 | awk '{print $3}'))))/bin/mpirun -wdir $PWD"
  NNODES=$SLURM_JOB_NUM_NODES
  MPITASKS_PER_NODE=$((SLURM_NTASKS/SLURM_JOB_NUM_NODES))
  MPI_TASKS=$SLURM_NTASKS
  MPIRUN="/home/gmap/mrpm/marguina/SAVE/mpiauto/mpiauto -np $MPI_TASKS -nnp $MPITASKS_PER_NODE --"
  export OMP_STACKSIZE=4G
  export KMP_STACKSIZE=4G
  export KMP_MONITOR_STACKSIZE=4G
  export DR_HOOK=1
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
  #export ECCODES_SAMPLES_PATH=/opt/softs/libraries/ICC_2018.5.274/eccodes-2.17.0/share/eccodes/ifs_samples/grib1
  #export ECCODES_DEFINITION_PATH=/opt/softs/libraries/ICC_2018.5.274/eccodes-2.17.0/share/eccodes/definitions
fi
echo $MPIRUN
set +x
if [ ! -f MASTER ] ; then echo No executable MASTER;exit 1;fi

#      ***************
#      *  Execution  *
#      ***************

echo
echo OMP_NUM_THREADS=$OMP_NUM_THREADS
set -x
ulimit -s unlimited
$MPIRUN $PWD/MASTER >lola
set +x
echo
##if [ -f lola ] ; then
##  echo;echo Standard output :;echo;cat lola
##fi
##if [ -f stderr.* ] ; then
##  for file in stderr.* ; do
##    echo;echo $file :;cat $file
##  done
##fi
##if [ -f stdout.* ] ; then
##echo;echo stdout :;echo;cat stdout.*
##fi
##if [ -a NODE.001_01 ] ; then
##  for file in NODE* ; do
##    echo;echo Listing $file;echo
##    cat $file
##  done
##fi
##if [ $(find . -name "drhook.prof.*" | wc -l) -ne 0 ] ; then
### Top 25 for each MPI task :
##  for file in drhook.prof.* ; do
##    echo;echo $file :;head -38 $file
##  done
##fi
#cat drhook.prof.* | perl -w $HOME/bin/drhook_merge_walltime_max.pl

#      *******************
#      *  Sauvegardes    *
#      *******************

ls
#if [ -f PFFPOS000+0000 ] ; then
#  cp PFFPOS000+0000 $WAIT_QUEUE/PFFPOS000+0000.$PBS_JOBID
#fi
set +e
cp drhook.prof.* $OUTPUTDIR/
cp lola NODE.001_01 ICMSHFPOS+00* DHFDLFPOS+00* $OUTPUTDIR/
/bin/rm fort.4 lola ICMSHFPOS+0000* PFFPOSFRANGP0025+0000* ICMSHFPOS+0001*
/bin/rm PFFPOSFRANGP0025+0001* ICMSHFPOS+0002* ECHIS PFFPOSFRANGP0025+0002* DHFDLFPOS+00*
/bin/rm ECHFP NODE.001_01 ifs.stat $(tar tfz $TESTDIR/data/rtm/rrtm.const.04.tgz)

#      ****************
#      *  Epilogue    *
#      ****************

ls -ltr | grep -v "\->"
echo Wait_queue :
ls -ltr $TMPWAIT
cd $TMPDIR
\rm -rf rundir.$$
\rm -rf wait_queue.$$
date
set +x
