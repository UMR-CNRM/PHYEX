PROGRAM MAIN_TURB

USE XRD_GETOPTIONS
USE GETDATA_TURB_MOD
USE COMPUTE_DIFF
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_CST,        ONLY: CST
USE MODD_CTURB
USE MODD_LES,        ONLY: TLES
USE MODD_TURB_n,     ONLY: TURBN
USE MODD_IO,         ONLY: TFILEDATA
USE MODI_TURB
USE MODI_INI_CST
USE MODD_BUDGET!, ONLY: TBUCONF_ASSOCIATE, TBUDGETDATA, NBUDGET_RH, TBUCONF
USE STACK_MOD
USE OMP_LIB
USE YOMHOOK, ONLY : LHOOK, DR_HOOK
USE PARKIND1, ONLY : JPRB, JPIM


IMPLICIT NONE

INTEGER      :: KLON 
INTEGER      :: KLEV
INTEGER      :: KRR, KRRL, KRRI
INTEGER      :: KSV

!IN and INOUTS
REAL, ALLOCATABLE   :: ZDXX               (:,:,:,:)
REAL, ALLOCATABLE   :: ZDYY               (:,:,:,:)
REAL, ALLOCATABLE   :: ZDZZ               (:,:,:,:)
REAL, ALLOCATABLE   :: ZDZX               (:,:,:,:)
REAL, ALLOCATABLE   :: ZDZY               (:,:,:,:)
REAL, ALLOCATABLE   :: ZZZ                (:,:,:,:)
REAL, ALLOCATABLE   :: ZDIRCOSXW          (:,:,:)
REAL, ALLOCATABLE   :: ZDIRCOSYW          (:,:,:)
REAL, ALLOCATABLE   :: ZDIRCOSZW          (:,:,:)
REAL, ALLOCATABLE   :: ZCOSSLOPE          (:,:,:)
REAL, ALLOCATABLE   :: ZSINSLOPE          (:,:,:)
REAL, ALLOCATABLE   :: PRHODJ             (:,:,:,:)
REAL, ALLOCATABLE   :: PTHVREF            (:,:,:,:)
REAL, ALLOCATABLE   :: PSFTH              (:,:,:)
REAL, ALLOCATABLE   :: PSFRV              (:,:,:)
REAL, ALLOCATABLE   :: PSFU               (:,:,:)
REAL, ALLOCATABLE   :: PSFV               (:,:,:)
REAL, ALLOCATABLE   :: PSFSV              (:,:,:,:) !(KLON, 1, KSV)
REAL, ALLOCATABLE   :: PPABSM             (:,:,:,:)
REAL, ALLOCATABLE   :: PUM                (:,:,:,:)
REAL, ALLOCATABLE   :: PVM                (:,:,:,:)
REAL, ALLOCATABLE   :: PWM                (:,:,:,:)
REAL, ALLOCATABLE   :: PTKEM              (:,:,:,:)
REAL, ALLOCATABLE   :: ZSVM               (:,:,:,:,:) !(KLON,1,KLEV+2,KSV)
REAL, ALLOCATABLE   :: PSRCM              (:,:,:,:)
REAL, ALLOCATABLE   :: PLENGTHM           (:,:,:,:)
REAL, ALLOCATABLE   :: PLENGTHH           (:,:,:,:)
REAL, ALLOCATABLE   :: MFMOIST            (:,:,:,:)
REAL, ALLOCATABLE   :: ZBL_DEPTH          (:,:,:)
REAL, ALLOCATABLE   :: ZSBL_DEPTH         (:,:,:)
REAL, ALLOCATABLE   :: ZCEI               (:,:,:,:)
REAL, ALLOCATABLE   :: PTHM               (:,:,:,:)
REAL, ALLOCATABLE   :: ZRM                (:,:,:,:,:) !(KLON,1,KLEV+2,KRR)
REAL, ALLOCATABLE   :: PRUS               (:,:,:,:)
REAL, ALLOCATABLE   :: PRVS               (:,:,:,:)
REAL, ALLOCATABLE   :: PRWS               (:,:,:,:)
REAL, ALLOCATABLE   :: PRTHS              (:,:,:,:)
REAL, ALLOCATABLE   :: ZRRS               (:,:,:,:,:) !(KLON,1,KLEV+2,KRR)
REAL, ALLOCATABLE   :: ZRSVS              (:,:,:,:,:) !(KLON,1,KLEV+2,KSV)
REAL, ALLOCATABLE   :: PRTKES_OUT         (:,:,:,:)
REAL, ALLOCATABLE   :: PFLXZTHVMF         (:,:,:,:)
REAL, ALLOCATABLE   :: PHGRAD             (:,:,:,:,:) !(KLON,1,KLEV+2,KGRADIENTS)
REAL, ALLOCATABLE   :: PZS                (:,:,:)

!OUT
REAL, ALLOCATABLE   :: PSIGS              (:,:,:,:)
REAL, ALLOCATABLE   :: ZWTH               (:,:,:,:)
REAL, ALLOCATABLE   :: ZWRC               (:,:,:,:)
REAL, ALLOCATABLE   :: ZWSV               (:,:,:,:,:) !(KLON,1,KLEV+2,KSV)
REAL, ALLOCATABLE   :: PDP                (:,:,:,:)
REAL, ALLOCATABLE   :: PTP                (:,:,:,:)
REAL, ALLOCATABLE   :: PTDIFF             (:,:,:,:)
REAL, ALLOCATABLE   :: PTDISS             (:,:,:,:)
REAL, ALLOCATABLE   :: PEDR               (:,:,:,:)
REAL, ALLOCATABLE   :: PTPMF              (:,:,:,:)
REAL, ALLOCATABLE   :: PDRUS_TURB         (:,:,:,:)
REAL, ALLOCATABLE   :: PDRVS_TURB         (:,:,:,:)
REAL, ALLOCATABLE   :: PDRTHLS_TURB       (:,:,:,:)
REAL, ALLOCATABLE   :: PDRRTS_TURB        (:,:,:,:)
REAL, ALLOCATABLE   :: ZDRSVS_TURB        (:,:,:,:,:) !(KLON,1,KLEV+2,KSV)

!Expected values
REAL, ALLOCATABLE   :: ZBL_DEPTH_OUT      (:,:,:)
REAL, ALLOCATABLE   :: ZSBL_DEPTH_OUT     (:,:,:)
REAL, ALLOCATABLE   :: PTHM_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: ZRM_OUT            (:,:,:,:,:) !(KLON,1,KLEV+2,KRR)
REAL, ALLOCATABLE   :: PRUS_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: PRVS_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: PRWS_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: PRTHS_OUT          (:,:,:,:)
REAL, ALLOCATABLE   :: ZRRS_OUT           (:,:,:,:,:) !(KLON,1,KLEV+2,KRR)
REAL, ALLOCATABLE   :: ZRSVS_OUT          (:,:,:,:,:) !(KLON,1,KLEV+2,KSV)
REAL, ALLOCATABLE   :: PRTKES_OUT_OUT     (:,:,:,:)
REAL, ALLOCATABLE   :: PSIGS_OUT          (:,:,:,:)
REAL, ALLOCATABLE   :: ZWTH_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: ZWRC_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: ZWSV_OUT           (:,:,:,:,:) !(KLON,1,KLEV+2,KSV)
REAL, ALLOCATABLE   :: PDP_OUT            (:,:,:,:)
REAL, ALLOCATABLE   :: PTP_OUT            (:,:,:,:)
REAL, ALLOCATABLE   :: PTDIFF_OUT         (:,:,:,:)
REAL, ALLOCATABLE   :: PTDISS_OUT         (:,:,:,:)
REAL, ALLOCATABLE   :: PEDR_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: PTPMF_OUT          (:,:,:,:)
REAL, ALLOCATABLE   :: PDRUS_TURB_OUT     (:,:,:,:)
REAL, ALLOCATABLE   :: PDRVS_TURB_OUT     (:,:,:,:)
REAL, ALLOCATABLE   :: PDRTHLS_TURB_OUT   (:,:,:,:)
REAL, ALLOCATABLE   :: PDRRTS_TURB_OUT    (:,:,:,:)
REAL, ALLOCATABLE   :: ZDRSVS_TURB_OUT    (:,:,:,:,:) !(KLON,1,KLEV+2,KSV)

INTEGER :: NPROMA, NGPBLKS, NFLEVG
INTEGER :: IBL, JLON, JLEV

TYPE(DIMPHYEX_t)         :: D, D0
INTEGER                  :: IMI, ISPLIT, KSV_LGBEG, KSV_LGEND, KGRADIENTS
INTEGER                  :: KSV_LIMA_NR, KSV_LIMA_NS, KSV_LIMA_NG, KSV_LIMA_NH
CHARACTER(LEN=4),DIMENSION(2)  :: HLBCX, HLBCY
CHARACTER(LEN=6)   ::  HPROGRAM
LOGICAL            :: O2D, ONOMIXLG, OFLAT, OCOUPLES, OBLOWSNOW, OCOMPUTE_SRC, OOCEAN, ODEEPOC
TYPE(TFILEDATA)    :: ZTFILE
REAL :: ZCEI_MAX, ZCEI_MIN, ZCOEF_AMPL_SAT
CHARACTER (LEN=4)   :: CMICRO  
REAL                :: PTSTEP 
TYPE(TBUDGETDATA), DIMENSION(NBUDGET_RH) :: YLBUDGET
LOGICAL                  :: LLCHECK
LOGICAL                  :: LLCHECKDIFF
LOGICAL                  :: LLDIFF
INTEGER                  :: IBLOCK1, IBLOCK2
INTEGER                  :: ISTSZ, JBLK1, JBLK2
INTEGER                  :: NTID, ITID
INTEGER                  :: JRR


REAL, ALLOCATABLE :: PSTACK(:,:)
TYPE (STACK) :: YLSTACK

REAL(KIND=8) :: TS,TE
REAL(KIND=8) :: TSC, TEC, TSD, TED, ZTC, ZTD 
INTEGER :: ITIME, NTIME
INTEGER :: IRANK, ISIZE
LOGICAL :: LLVERBOSE, LLSTAT, LLBIND
REAL (KIND=JPRB) :: ZHOOK_HANDLE
CHARACTER(LEN=32) :: CLTEXT

CALL INITOPTIONS ()
NGPBLKS = 150
CALL GETOPTION ("--blocks", NGPBLKS)
NPROMA = 32
CALL GETOPTION ("--nproma", NPROMA)
NFLEVG = -1
CALL GETOPTION ("--nflevg", NFLEVG)
CALL GETOPTION ("--check",  LLCHECK)
CALL GETOPTION ("--checkdiff",  LLCHECKDIFF)
IBLOCK1 = 1
CALL GETOPTION ("--check-block-1", IBLOCK1)
IBLOCK2 = NGPBLKS
CALL GETOPTION ("--check-block-2", IBLOCK2)
CALL GETOPTION ("--stat", LLSTAT)
NTIME = 1
CALL GETOPTION ("--times", NTIME)
CALL GETOPTION ("--verbose", LLVERBOSE)
CALL GETOPTION ("--bind", LLBIND)
CALL CHECKOPTIONS ()

LLDIFF = .FALSE.

IRANK = 0
ISIZE = 1
IF (LLBIND) THEN
  CALL LINUX_BIND      (IRANK, ISIZE)
  CALL LINUX_BIND_DUMP (IRANK, ISIZE)
ENDIF

CALL GETDATA_TURB (NPROMA, NGPBLKS, NFLEVG, KRR, KRRL, KRRI, KSV, KLEV, &
                  !IN and INOUT (allocation and values are needed for the call)
                  &ZDXX, ZDYY, ZDZZ, ZDZX, ZDZY, ZZZ, &
                  &ZDIRCOSXW, ZDIRCOSYW, ZDIRCOSZW, ZCOSSLOPE, ZSINSLOPE, &
                  &PRHODJ, PTHVREF, &
                  &PSFTH, PSFRV, PSFU, PSFV, PSFSV, &
                  &PPABSM, PUM, PVM, PWM, PTKEM, ZSVM, PSRCM, &
                  &PLENGTHM, PLENGTHH, MFMOIST, &
                  &ZBL_DEPTH, ZSBL_DEPTH, &
                  &ZCEI, &
                  &PTHM, ZRM, &
                  &PRUS, PRVS, PRWS, PRTHS, ZRRS, ZRSVS, PRTKES_OUT, &
                  &PFLXZTHVMF, &
                  &PHGRAD, PZS, &
                  !OUT only (needed to allocate the array to be passed to the subroutine)
                  &PSIGS, &
                  &ZWTH,ZWRC,ZWSV,PDP,PTP,PTDIFF,PTDISS, &
                  &PEDR,PTPMF, &
                  &PDRUS_TURB,PDRVS_TURB, &
                  &PDRTHLS_TURB,PDRRTS_TURB,ZDRSVS_TURB, &
                  !OUT and INOUT (expected values)
                  &ZBL_DEPTH_OUT, ZSBL_DEPTH_OUT, &
                  &PTHM_OUT, ZRM_OUT, &
                  &PRUS_OUT, PRVS_OUT, PRWS_OUT, PRTHS_OUT, ZRRS_OUT, ZRSVS_OUT, PRTKES_OUT_OUT, &
                  &PSIGS_OUT, &
                  &ZWTH_OUT, ZWRC_OUT, ZWSV_OUT, PDP_OUT, PTP_OUT, PTDIFF_OUT, PTDISS_OUT, &
                  &PEDR_OUT, PTPMF_OUT, &
                  &PDRUS_TURB_OUT, PDRVS_TURB_OUT, &
                  &PDRTHLS_TURB_OUT, PDRRTS_TURB_OUT, ZDRSVS_TURB_OUT)

IF (LLVERBOSE) PRINT *, " KLEV = ", KLEV, " KRR = ", KRR

PRINT *, " NPROMA = ", NPROMA, " KLEV = ", KLEV, " NGPBLKS = ", NGPBLKS

IMI = 1
HLBCX(:)='CYCLCYCL'
HLBCY(:)='CYCLCYCL'
ISPLIT = 1
KSV_LGBEG = 0
KSV_LGEND = 0
HPROGRAM='AROME '
O2D=.FALSE.
ONOMIXLG=.FALSE.
OFLAT=.FALSE.
OCOUPLES=.FALSE.
OBLOWSNOW=.FALSE.
OCOMPUTE_SRC=SIZE(PSIGS, 3)/=0
OOCEAN=.FALSE.
ODEEPOC=.FALSE.
CMICRO='ICE3'
ZTFILE%LOPENED=.FALSE.
ZCEI_MAX=1.0
ZCEI_MIN=0.0
ZCOEF_AMPL_SAT=0.0
KGRADIENTS=0
KSV_LIMA_NR=0
KSV_LIMA_NS=0
KSV_LIMA_NG=0
KSV_LIMA_NH=0
TLES%LLES=.FALSE.
!
PTSTEP = 25.0000000000000

CALL INIT_PHYEX ()

DO JRR=1, NBUDGET_RH
  YLBUDGET(JRR)%NBUDGET=JRR
ENDDO

D0%NIT  = NPROMA
D0%NIB  = 1
D0%NIE  = NPROMA
D0%NJT  = 1
D0%NJB  = 1
D0%NJE  = 1
D0%NIJT = D0%NIT * D0%NJT
D0%NIJB = 1
D0%NIJE = NPROMA
D0%NKL  = -1
D0%NKT  = KLEV
D0%NKA  = KLEV
D0%NKU  = 1
D0%NKB  = KLEV-1
D0%NKE  = 2
D0%NKTB = 2
D0%NKTE = KLEV-1
D0%NIBC = 1
D0%NJBC = 1
D0%NIEC = D0%NIE
D0%NJEC = D0%NJT

ISTSZ = NPROMA * 20 * KLEV
ALLOCATE (PSTACK (ISTSZ, NGPBLKS))

TS = OMP_GET_WTIME ()

ZTD = 0.
ZTC = 0.

IF (LHOOK) CALL DR_HOOK ('MAIN',0,ZHOOK_HANDLE)

DO ITIME = 1, NTIME

  TSD = OMP_GET_WTIME ()

!!!              !directives pas a jour !$acc data &
!!!              !directives pas a jour !$acc      & copyin  (D0, CST, ICEP, NEB, KRR, HFRAC_ICE, HCONDENS, HLAMBDA3, HBUNAME, OSUBG_COND, OSIGMAS, OCND2, HSUBG_MF_PDF, PTSTEP, LMFCONV, &
!!!              !directives pas a jour !$acc      &          ZSIGQSAT, PRHODJ, PEXNREF, PRHODREF, PSIGS, PMFCONV, PPABSM, ZZZ, PCF_MF, PRC_MF, PRI_MF, ZRS, ZICE_CLD_WGT) &
!!!              !directives pas a jour !$acc      & copy    (PRS, PTHS), &
!!!              !directives pas a jour !$acc      & copyout (PSRCS, PCLDFR, PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF) &
!!!              !directives pas a jour !$acc      & create  (PSTACK) 

  TSC = OMP_GET_WTIME ()

#ifdef USE_OPENMP
!$OMP PARALLEL PRIVATE (D, ITID, JBLK1, JBLK2)
#endif

#ifdef _OPENACC
JBLK1 = 1 
JBLK2 = NGPBLKS
#endif

#ifdef USE_OPENMP
NTID = OMP_GET_MAX_THREADS ()
ITID = OMP_GET_THREAD_NUM ()
JBLK1 = 1 +  (NGPBLKS * (ITID+0)) / NTID
JBLK2 =      (NGPBLKS * (ITID+1)) / NTID


!PRINT *, ITID, JBLK1, JBLK2

#endif

!$acc parallel loop gang vector private (YLSTACK, IBL, JLON, D) collapse (2)

  DO IBL = JBLK1, JBLK2


#ifdef _OPENACC
  DO JLON = 1, NPROMA
    D = D0
    D%NIB = JLON
    D%NIE = JLON
    D%NIJB = JLON
    D%NIJE = JLON
    D%NIBC = JLON
    D%NIEC = JLON
#endif

#ifdef USE_OPENMP
    D = D0
#endif

#ifdef USE_STACK
    YLSTACK%L = LOC (PSTACK (1, IBL))
    YLSTACK%U = YLSTACK%L + ISTSZ * KIND (PSTACK)
#else
    YLSTACK%L = 0
    YLSTACK%U = 0
#endif

CALL TURB (CST,CSTURB,TBUCONF,TURBN, D, TLES,&
   & IMI, KRR, KRRL, KRRI, HLBCX, HLBCY, KGRADIENTS, 1,&
   & ISPLIT,IMI, KSV, KSV_LGBEG, KSV_LGEND, &
   & HPROGRAM, &
   & KSV_LIMA_NR, KSV_LIMA_NS, KSV_LIMA_NG, KSV_LIMA_NH,&
   & O2D, ONOMIXLG, OFLAT, OCOUPLES,OBLOWSNOW,.FALSE.,.TRUE.,& 
   & OCOMPUTE_SRC, 1.0, &
   & OOCEAN,ODEEPOC, .FALSE.,   &
   & 'NONE',CMICRO,           &
   & 2*PTSTEP,ZTFILE,                                      &
   & ZDXX(:,:,:,IBL),ZDYY(:,:,:,IBL),ZDZZ(:,:,:,IBL),ZDZX(:,:,:,IBL),ZDZY(:,:,:,IBL),ZZZ(:,:,:,IBL),          &
   & ZDIRCOSXW,ZDIRCOSYW,ZDIRCOSZW,ZCOSSLOPE,ZSINSLOPE,    &
   & PRHODJ(:,:,:,IBL),PTHVREF(:,:,:,IBL), PHGRAD, PZS,                             &
   & PSFTH(:,:,IBL),PSFRV(:,:,IBL),PSFSV(:,:,:,IBL),PSFU(:,:,IBL),PSFV(:,:,IBL),                          &
   & PPABSM(:,:,:,IBL),PUM(:,:,:,IBL),PVM(:,:,:,IBL),PWM(:,:,:,IBL),PTKEM(:,:,:,IBL),ZSVM(:,:,:,:,IBL),PSRCM(:,:,:,IBL),                  &
   & PLENGTHM(:,:,:,IBL),PLENGTHH(:,:,:,IBL),MFMOIST(:,:,:,IBL),                            &
   & ZBL_DEPTH(:,:,IBL),ZSBL_DEPTH(:,:,IBL),                                 &
   & ZCEI(:,:,:,IBL),ZCEI_MIN,ZCEI_MAX,ZCOEF_AMPL_SAT,    &
   & PTHM(:,:,:,IBL),ZRM(:,:,:,:,IBL), &
   & PRUS(:,:,:,IBL),PRVS(:,:,:,IBL),PRWS(:,:,:,IBL),PRTHS(:,:,:,IBL),ZRRS(:,:,:,:,IBL),ZRSVS(:,:,:,:,IBL),PRTKES_OUT(:,:,:,IBL),         &
   & PSIGS(:,:,:,IBL),                                         &
   & PFLXZTHVMF(:,:,:,IBL),ZWTH(:,:,:,IBL),ZWRC(:,:,:,IBL),ZWSV(:,:,:,:,IBL),PDP(:,:,:,IBL),PTP(:,:,:,IBL),PTDIFF(:,:,:,IBL),PTDISS(:,:,:,IBL),&
   & YLBUDGET, KBUDGETS=SIZE(YLBUDGET),PEDR=PEDR(:,:,:,IBL),PTPMF=PTPMF(:,:,:,IBL),&
   & PDRUS_TURB=PDRUS_TURB(:,:,:,IBL),PDRVS_TURB=PDRVS_TURB(:,:,:,IBL),          &
   & PDRTHLS_TURB=PDRTHLS_TURB(:,:,:,IBL),PDRRTS_TURB=PDRRTS_TURB(:,:,:,IBL),PDRSVS_TURB=ZDRSVS_TURB(:,:,:,:,IBL))

#ifdef _OPENACC
    ENDDO
#endif

  ENDDO

#ifdef USE_OPENMP
!$OMP END PARALLEL
#endif

!$acc end parallel loop

  TEC = OMP_GET_WTIME ()

!$acc end data

  TED = OMP_GET_WTIME ()

  ZTC = ZTC + (TEC - TSC)
  ZTD = ZTD + (TED - TSD)

ENDDO

IF (LHOOK) CALL DR_HOOK ('MAIN',1,ZHOOK_HANDLE)

TE = OMP_GET_WTIME()

WRITE (*,'(A,F8.2,A)') 'elapsed time : ',TE-TS,' s'
WRITE (*,'(A,F8.4,A)') '          i.e. ',1000.*(TE-TS)/(NPROMA*NGPBLKS)/NTIME,' ms/gp'

PRINT *, " ZTD = ", ZTD, ZTD / REAL (NPROMA*NGPBLKS*NTIME)
PRINT *, " ZTC = ", ZTC, ZTC / REAL (NPROMA*NGPBLKS*NTIME)

IF (LLCHECK .OR. LLSTAT .OR. LLCHECKDIFF) THEN
  DO IBL = IBLOCK1, IBLOCK2
    PRINT *, " IBL = ", IBL
    DO JRR=1, KRR
      WRITE (CLTEXT, '("ZRM JRR=",I3.3)') JRR
      CALL DIFF3 (CLTEXT,      ZRM_OUT       (:,:,:,JRR,IBL), ZRM      (:,:,:,JRR,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
      WRITE (CLTEXT, '("ZRRS JRR=",I3.3)') JRR
      CALL DIFF3 (CLTEXT,      ZRRS_OUT       (:,:,:,JRR,IBL), ZRRS      (:,:,:,JRR,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    ENDDO
    CALL DIFF2 ("ZBL_DEPTH   ", ZBL_DEPTH_OUT    (:,:,IBL)  , ZBL_DEPTH   (:,:,IBL)  , LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF2 ("ZSBL_DEPTH  ", ZSBL_DEPTH_OUT   (:,:,IBL)  , ZSBL_DEPTH  (:,:,IBL)  , LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PTHM        ", PTHM_OUT         (:,:,:,IBL), PTHM        (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PRUS        ", PRUS_OUT         (:,:,:,IBL), PRUS        (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PRVS        ", PRVS_OUT         (:,:,:,IBL), PRVS        (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PRWS        ", PRWS_OUT         (:,:,:,IBL), PRWS        (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PRTHS       ", PRTHS_OUT        (:,:,:,IBL), PRTHS       (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PRTKES_OUT  ", PRTKES_OUT_OUT   (:,:,:,IBL), PRTKES_OUT  (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PSIGS       ", PSIGS_OUT        (:,:,:,IBL), PSIGS       (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("ZWTH        ", ZWTH_OUT         (:,:,:,IBL), ZWTH        (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("ZWRC        ", ZWRC_OUT         (:,:,:,IBL), ZWRC        (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PDP         ", PDP_OUT          (:,:,:,IBL), PDP         (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PTP         ", PTP_OUT          (:,:,:,IBL), PTP         (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PTDIFF      ", PTDIFF_OUT       (:,:,:,IBL), PTDIFF      (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PTDISS      ", PTDISS_OUT       (:,:,:,IBL), PTDISS      (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PEDR        ", PEDR_OUT         (:,:,:,IBL), PEDR        (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PTPMF       ", PTPMF_OUT        (:,:,:,IBL), PTPMF       (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PDRUS_TURB  ", PDRUS_TURB_OUT   (:,:,:,IBL), PDRUS_TURB  (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PDRVS_TURB  ", PDRVS_TURB_OUT   (:,:,:,IBL), PDRVS_TURB  (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PDRTHLS_TURB", PDRTHLS_TURB_OUT (:,:,:,IBL), PDRTHLS_TURB(:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PDRRTS_TURB ", PDRRTS_TURB_OUT  (:,:,:,IBL), PDRRTS_TURB (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
  ENDDO
ENDIF

IF (LLCHECKDIFF) THEN
  IF (LLDIFF) THEN
    PRINT*, "THERE ARE DIFF SOMEWHERE"
  ELSE
    PRINT*, "THERE IS NO DIFF AT ALL"
  ENDIF
ENDIF

STOP

CONTAINS

SUBROUTINE INIT_PHYEX()
USE MODD_TURB_N, ONLY: TURB_GOTO_MODEL
IMPLICIT NONE
!
CALL INI_CST
CALL TURB_GOTO_MODEL(1,1)
CALL CTURB_ASSOCIATE()
CALL TBUCONF_ASSOCIATE
LBU_ENABLE=.FALSE.                                                                                                       
LBUDGET_U=.FALSE.
LBUDGET_V=.FALSE.
LBUDGET_W=.FALSE.
LBUDGET_TH=.FALSE.
LBUDGET_TKE=.FALSE.
LBUDGET_RV=.FALSE.
LBUDGET_RC=.FALSE.
LBUDGET_RR=.FALSE.
LBUDGET_RI=.FALSE.
LBUDGET_RS=.FALSE.
LBUDGET_RG=.FALSE.
LBUDGET_RH=.FALSE.
LBUDGET_SV=.FALSE.

TURBN%LHARAT=.FALSE.
TURBN%CTURBDIM = '1DIM'
TURBN%XIMPL=1.
TURBN%CTURBLEN='BL89'
TURBN%LSTATNW=.FALSE.
TURBN%LTURB_DIAG=.FALSE.
TURBN%LTURB_FLX=.FALSE.
TURBN%LSUBG_COND=.TRUE.
TURBN%LRMC01=.FALSE.
TURBN%CTOM='NONE'
TURBN%LLEONARD=.FALSE.

XCED  = 0.85
XCEP  = 2.11 
XA0   = 0.6
XA2   = 1.
XA3   = 0.
XCTD  = 1.2
IF (TURBN%LSTATNW) THEN
    XCTP  = 4.0
  ELSE
    XCTP  = 4.65
ENDIF
XA5   = 1./3.
XCET  = 0.40
XALPSBL = 4.63
XRM17 = 0.5  ! Rodier et al 2017
XCMFS= 2./3./XCEP*(1.-XA0)   !Constant for the momentum flux due to shear (RS)
XCSHF= 2./3./XCTP            !Constant for the sensible heat flux(RS)
XCHF= XCSHF                  !Constant for the humidity flux(RS)
XCTV= 2./3./XCTP/XCTD        !Constant for the temperature variance(RS)
XCHV=  XCTV                  !Constant for the humidity variance(RS)
XCHT1= XCTV/2.      !Constants for the temperature-humidity correlation(RS)
XCHT2= XCTV/2.
XCPR1= XCTV         !Constants for the turbulent Prandtl and Schmidt numbers
XCPR2= XCHT1
XCPR3= XCPR2        ! used only for the Schmidt number for scalar variables
XCPR4= XCPR2
XCPR5= XCPR2
XTKEMIN=1.E-6
!XLINI=10.   ! BL mixing length
XLINI=0.1   ! BL mixing length
XLINF=1.E-10! to prevent division by zero
XPHI_LIM = 3.
XCDP  =  1.46
XCDD  =  1.83
XCDT  =  0.42
XSBL_O_BL     = 0.05 ! SBL height / BL height ratio
XFTOP_O_FSURF = 0.05 ! Fraction of surface (heat or momentum) flux used to define top of BL

!
END SUBROUTINE INIT_PHYEX

END PROGRAM

