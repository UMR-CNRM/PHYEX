PROGRAM MAIN_SHALLOW

USE XRD_GETOPTIONS
USE GETDATA_SHALLOW_MOD
USE COMPUTE_DIFF
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_CST,        ONLY: CST_t
USE MODD_NEB,        ONLY: NEB
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALLN, PARAM_MFSHALL_GOTO_MODEL
USE MODD_CTURB
USE MODD_TURB_n,     ONLY: TURBN
USE MODI_SHALLOW_MF
USE MODI_INI_NEB
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
REAL, ALLOCATABLE   :: PDZZF              (:,:,:,:)
REAL, ALLOCATABLE   :: PZZ               (:,:,:,:)
REAL, ALLOCATABLE   :: PRHODJ              (:,:,:,:)
REAL, ALLOCATABLE   :: PRHODREF               (:,:,:,:)
REAL, ALLOCATABLE   :: PPABSM               (:,:,:,:)
REAL, ALLOCATABLE   :: PEXNM               (:,:,:,:)
REAL, ALLOCATABLE   :: PSFTH          (:,:,:)
REAL, ALLOCATABLE   :: PSFRV          (:,:,:)
REAL, ALLOCATABLE   :: PTHM             (:,:,:,:)
REAL, ALLOCATABLE   :: PRM              (:,:,:,:,:) !(KLON, 1, KLEV, KRR)

REAL, ALLOCATABLE   :: PUM                (:,:,:,:)
REAL, ALLOCATABLE   :: PVM                (:,:,:,:)
REAL, ALLOCATABLE   :: PTKEM              (:,:,:,:)
REAL, ALLOCATABLE   :: PSVM               (:,:,:,:,:) !(KLON,1,KLEV,KSV)
REAL, ALLOCATABLE   :: PTHL_UP              (:,:,:,:)
REAL, ALLOCATABLE   :: PRT_UP           (:,:,:,:)
REAL, ALLOCATABLE   :: PRV_UP           (:,:,:,:)
REAL, ALLOCATABLE   :: PRC_UP            (:,:,:,:)
REAL, ALLOCATABLE   :: PRI_UP               (:,:,:,:)
REAL, ALLOCATABLE   :: PU_UP               (:,:,:,:)
REAL, ALLOCATABLE   :: PV_UP               (:,:,:,:)
REAL, ALLOCATABLE   :: PTHV_UP               (:,:,:,:)
REAL, ALLOCATABLE   :: PW_UP               (:,:,:,:)
REAL, ALLOCATABLE   :: PFRAC_UP              (:,:,:,:)
REAL, ALLOCATABLE   :: PEMF         (:,:,:,:)

!OUT
REAL, ALLOCATABLE   :: PDUDT_MF              (:,:,:,:)
REAL, ALLOCATABLE   :: PDVDT_MF               (:,:,:,:)
REAL, ALLOCATABLE   :: PDTHLDT_MF               (:,:,:,:)
REAL, ALLOCATABLE   :: PDRTDT_MF               (:,:,:,:)
REAL, ALLOCATABLE   :: PDSVDT_MF               (:,:,:,:,:) !(KLON,1,KLEV,KSV)
REAL, ALLOCATABLE   :: PSIGMF                (:,:,:,:)
REAL, ALLOCATABLE   :: PRC_MF                (:,:,:,:)
REAL, ALLOCATABLE   :: PRI_MF             (:,:,:,:)
REAL, ALLOCATABLE   :: PCF_MF             (:,:,:,:)
REAL, ALLOCATABLE   :: PFLXZTHVMF               (:,:,:,:)
REAL, ALLOCATABLE   :: ZFLXZTHMF              (:,:,:,:)
REAL, ALLOCATABLE   :: ZFLXZRMF         (:,:,:,:)
REAL, ALLOCATABLE   :: ZFLXZUMF         (:,:,:,:)
REAL, ALLOCATABLE   :: ZFLXZVMF       (:,:,:,:)
REAL, ALLOCATABLE   :: ZDETR        (:,:,:,:)
REAL, ALLOCATABLE   :: ZENTR        (:,:,:,:)
INTEGER, ALLOCATABLE:: IKLCL  (:,:,:)
INTEGER, ALLOCATABLE:: IKETL  (:,:,:)
INTEGER, ALLOCATABLE:: IKCTL (:,:,:)

!Expected values
REAL, ALLOCATABLE   :: PDUDT_MF_OUT      (:,:,:,:)
REAL, ALLOCATABLE   :: PDVDT_MF_OUT     (:,:,:,:)
REAL, ALLOCATABLE   :: PDTHLDT_MF_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: PDRTDT_MF_OUT            (:,:,:,:)
REAL, ALLOCATABLE   :: PDSVDT_MF_OUT           (:,:,:,:,:) !(KLON,1,KLEV,KSV)
REAL, ALLOCATABLE   :: PSIGMF_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: PRC_MF_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: PRI_MF_OUT          (:,:,:,:)
REAL, ALLOCATABLE   :: PCF_MF_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: PFLXZTHVMF_OUT          (:,:,:,:)
REAL, ALLOCATABLE   :: ZFLXZTHMF_OUT     (:,:,:,:)
REAL, ALLOCATABLE   :: ZFLXZRMF_OUT          (:,:,:,:)
REAL, ALLOCATABLE   :: ZFLXZUMF_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: ZFLXZVMF_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: PTHL_UP_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: PRT_UP_OUT            (:,:,:,:)
REAL, ALLOCATABLE   :: PRV_UP_OUT            (:,:,:,:)
REAL, ALLOCATABLE   :: PRC_UP_OUT         (:,:,:,:)
REAL, ALLOCATABLE   :: PRI_UP_OUT         (:,:,:,:)
REAL, ALLOCATABLE   :: PU_UP_OUT           (:,:,:,:)
REAL, ALLOCATABLE   :: PV_UP_OUT          (:,:,:,:)
REAL, ALLOCATABLE   :: PTHV_UP_OUT     (:,:,:,:)
REAL, ALLOCATABLE   :: PW_UP_OUT     (:,:,:,:)
REAL, ALLOCATABLE   :: PFRAC_UP_OUT   (:,:,:,:)
REAL, ALLOCATABLE   :: PEMF_OUT    (:,:,:,:)
REAL, ALLOCATABLE   :: ZDETR_OUT    (:,:,:,:)
REAL, ALLOCATABLE   :: ZENTR_OUT    (:,:,:,:)
INTEGER, ALLOCATABLE:: IKLCL_OUT  (:,:,:)
INTEGER, ALLOCATABLE:: IKETL_OUT  (:,:,:)
INTEGER, ALLOCATABLE:: IKCTL_OUT (:,:,:)

INTEGER :: NPROMA, NGPBLKS, NFLEVG
INTEGER :: IBL, JLON, JLEV

TYPE(DIMPHYEX_t)         :: D, D0
TYPE(CST_t)              :: CST
CHARACTER (LEN=4)        :: HMF_CLOUD, HMF_UPDRAFT
CHARACTER (LEN=1)        :: HFRAC_ICE
LOGICAL                  :: OMIXUV, ONOMIXLG, OSTATNW
REAL                     :: ZIMPL
INTEGER                  :: KSV_LGBEG, KSV_LGEND
REAL                     :: PTSTEP 
LOGICAL                  :: LLCHECK
LOGICAL                  :: LLCHECKDIFF
LOGICAL                  :: LLDIFF
INTEGER                  :: IBLOCK1, IBLOCK2
INTEGER                  :: ISTSZ, JBLK1, JBLK2
INTEGER                  :: NTID, ITID


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

CALL GETDATA_SHALLOW (NPROMA, NGPBLKS, NFLEVG, KRR, KRRL, KRRI, KSV, KLEV, &
                  !IN and INOUT (allocation and values are needed for the call)
                  &PDZZF, PZZ, PRHODJ, PRHODREF, PPABSM, PEXNM, &
                  &PSFTH, PSFRV, &
                  &PTHM, PRM, &
                  &PUM, PVM, PTKEM, PSVM, PTHL_UP, &
                  &PRT_UP, PRV_UP, PRC_UP, &
                  &PRI_UP, &
                  &PU_UP, &
                  &PV_UP, PTHV_UP, PW_UP, PFRAC_UP, PEMF, &
                  !OUT only (needed to allocate the array to be passed to the subroutine)
                  &PDUDT_MF, &
                  &PDVDT_MF,PDTHLDT_MF,PDRTDT_MF,PDSVDT_MF,PSIGMF,PRC_MF,PRI_MF,PCF_MF, &
                  &PFLXZTHVMF,ZFLXZTHMF, &
                  &ZFLXZRMF,ZFLXZUMF, &
                  &ZFLXZVMF,ZDETR,ZENTR, IKLCL, IKETL, IKCTL, &
                  !OUT and INOUT (expected values)
                  &PDUDT_MF_OUT, PDVDT_MF_OUT, &
                  &PDTHLDT_MF_OUT, PDRTDT_MF_OUT, &
                  &PDSVDT_MF_OUT, PSIGMF_OUT, PRC_MF_OUT, PRI_MF_OUT, PCF_MF_OUT, PFLXZTHVMF_OUT, ZFLXZTHMF_OUT, &
                  &ZFLXZRMF_OUT, &
                  &ZFLXZUMF_OUT, ZFLXZVMF_OUT, PTHL_UP_OUT, PRT_UP_OUT, PRV_UP_OUT, PRC_UP_OUT, PRI_UP_OUT, &
                  &PU_UP_OUT, PV_UP_OUT, &
                  &PTHV_UP_OUT, PW_UP_OUT, &
                  &PFRAC_UP_OUT, PEMF_OUT, ZDETR_OUT, ZENTR_OUT, IKLCL_OUT, IKETL_OUT, IKCTL_OUT)

IF (LLVERBOSE) PRINT *, " KLEV = ", KLEV, " KRR = ", KRR

PRINT *, " NPROMA = ", NPROMA, " KLEV = ", KLEV, " NGPBLKS = ", NGPBLKS

KSV_LGBEG = 0
KSV_LGEND = 0
HMF_CLOUD='DIRE'
HMF_UPDRAFT='EDKF'
HFRAC_ICE='S'
OMIXUV=.TRUE.
ONOMIXLG=.FALSE.
ZIMPL=1.
OSTATNW=.FALSE.
!
PTSTEP = 25.0000000000000

CALL INIT_PHYEX (20, &
                 CST)

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
!!!              !directives pas a jour !$acc      &          ZSIGQSAT, PTHM, PEXNREF, PRHODREF, PSIGS, PMFCONV, PPABSM, ZZZ, PCF_MF, PRC_MF, PRI_MF, ZRS, ZICE_CLD_WGT) &
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

  CALL SHALLOW_MF(D, CST, NEB, PARAM_MFSHALLN, TURBN, CSTURB,                    &
     &KRR=KRR, KRRL=KRRL, KRRI=KRRI, KSV=KSV,                                             &
     &HFRAC_ICE=HFRAC_ICE,ONOMIXLG=ONOMIXLG,KSV_LGBEG=KSV_LGBEG,KSV_LGEND=KSV_LGEND,      &
     &PIMPL_MF=ZIMPL, PTSTEP=PTSTEP,                                                      &
     &PDZZ=PDZZF(:,:,:,IBL),PZZ=PZZ(:,:,:,IBL),                                                                 &
     &PRHODJ=PRHODJ(:,:,:,IBL),PRHODREF=PRHODREF(:,:,:,IBL),                                                    &
     &PPABSM=PPABSM(:,:,:,IBL),PEXNM=PEXNM(:,:,:,IBL),                                                          &
     &PSFTH=PSFTH(:,:,IBL),PSFRV=PSFRV(:,:,IBL),                                                            &
     &PTHM=PTHM(:,:,:,IBL),PRM=PRM(:,:,:,:,IBL),PUM=PUM(:,:,:,IBL),PVM=PVM(:,:,:,IBL),&
     &PTKEM=PTKEM(:,:,:,IBL),PSVM=PSVM(:,:,:,:,IBL),                            &
     &PDUDT_MF=PDUDT_MF(:,:,:,IBL),PDVDT_MF=PDVDT_MF(:,:,:,IBL),                                                &
     &PDTHLDT_MF=PDTHLDT_MF(:,:,:,IBL),PDRTDT_MF=PDRTDT_MF(:,:,:,IBL),PDSVDT_MF=PDSVDT_MF(:,:,:,:,IBL),                      &
     &PSIGMF=PSIGMF(:,:,:,IBL),PRC_MF=PRC_MF(:,:,:,IBL),PRI_MF=PRI_MF(:,:,:,IBL),PCF_MF=PCF_MF(:,:,:,IBL),&
     &PFLXZTHVMF=PFLXZTHVMF(:,:,:,IBL),      &
     &PFLXZTHMF=ZFLXZTHMF(:,:,:,IBL),PFLXZRMF=ZFLXZRMF(:,:,:,IBL),PFLXZUMF=ZFLXZUMF(:,:,:,IBL),PFLXZVMF=ZFLXZVMF(:,:,:,IBL),     &
     &PTHL_UP=PTHL_UP(:,:,:,IBL),PRT_UP=PRT_UP(:,:,:,IBL),PRV_UP=PRV_UP(:,:,:,IBL),&
     &PRC_UP=PRC_UP(:,:,:,IBL),PRI_UP=PRI_UP(:,:,:,IBL),            &
     &PU_UP=PU_UP(:,:,:,IBL), PV_UP=PV_UP(:,:,:,IBL), PTHV_UP=PTHV_UP(:,:,:,IBL), PW_UP=PW_UP(:,:,:,IBL),                        &
     &PFRAC_UP=PFRAC_UP(:,:,:,IBL),PEMF=PEMF(:,:,:,IBL),PDETR=ZDETR(:,:,:,IBL),PENTR=ZENTR(:,:,:,IBL),                           &
     &KKLCL=IKLCL(:,:,IBL),KKETL=IKETL(:,:,IBL),KKCTL=IKCTL(:,:,IBL),PDX=0.,PDY=0.,KBUDGETS=0                                  )

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
    CALL DIFF3 ("PDUDT_MF   ", PDUDT_MF_OUT   (:,:,:,IBL), PDUDT_MF   (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PDVDT_MF   ", PDVDT_MF_OUT   (:,:,:,IBL), PDVDT_MF   (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PDTHLDT_MF ", PDTHLDT_MF_OUT (:,:,:,IBL), PDTHLDT_MF (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PDRTDT_MF  ", PDRTDT_MF_OUT  (:,:,:,IBL), PDRTDT_MF  (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PSIGMF     ", PSIGMF_OUT     (:,:,:,IBL), PSIGMF     (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PRC_MF     ", PRC_MF_OUT     (:,:,:,IBL), PRC_MF     (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PRI_MF     ", PRI_MF_OUT     (:,:,:,IBL), PRI_MF     (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PCF_MF     ", PCF_MF_OUT     (:,:,:,IBL), PCF_MF     (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PFLXZTHVMF ", PFLXZTHVMF_OUT (:,:,:,IBL), PFLXZTHVMF (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("ZFLXZTHMF  ", ZFLXZTHMF_OUT  (:,:,:,IBL), ZFLXZTHMF  (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("ZFLXZRMF   ", ZFLXZRMF_OUT   (:,:,:,IBL), ZFLXZRMF   (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("ZFLXZUMF   ", ZFLXZUMF_OUT   (:,:,:,IBL), ZFLXZUMF   (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("ZFLXZVMF   ", ZFLXZVMF_OUT   (:,:,:,IBL), ZFLXZVMF   (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PTHL_UP    ", PTHL_UP_OUT    (:,:,:,IBL), PTHL_UP    (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PRT_UP     ", PRT_UP_OUT     (:,:,:,IBL), PRT_UP     (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PRV_UP     ", PRV_UP_OUT     (:,:,:,IBL), PRV_UP     (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PRC_UP     ", PRC_UP_OUT     (:,:,:,IBL), PRC_UP     (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PRI_UP     ", PRI_UP_OUT     (:,:,:,IBL), PRI_UP     (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PU_UP      ", PU_UP_OUT      (:,:,:,IBL), PU_UP      (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PV_UP      ", PV_UP_OUT      (:,:,:,IBL), PV_UP      (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PTHV_UP    ", PTHV_UP_OUT    (:,:,:,IBL), PTHV_UP    (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PW_UP      ", PW_UP_OUT      (:,:,:,IBL), PW_UP      (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PFRAC_UP   ", PFRAC_UP_OUT   (:,:,:,IBL), PFRAC_UP   (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("PEMF       ", PEMF_OUT       (:,:,:,IBL), PEMF       (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("ZDETR      ", ZDETR_OUT      (:,:,:,IBL), ZDETR      (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF3 ("ZENTR      ", ZENTR_OUT      (:,:,:,IBL), ZENTR      (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
!    CALL DIFF2 ("IKLCL      ", IKLCL_OUT      (:,:,IBL),   IKLCL      (:,:,IBL),   LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
!    CALL DIFF2 ("IKETL      ", IKETL_OUT      (:,:,IBL),   IKETL      (:,:,IBL),   LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
!    CALL DIFF2 ("IKCTL      ", IKCTL_OUT      (:,:,IBL),   IKCTL      (:,:,IBL),   LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
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

SUBROUTINE INIT_PHYEX(KULOUT,&
                      CST)

USE MODD_CST, ONLY: CST_t
USE MODD_TURB_N, ONLY: TURB_GOTO_MODEL
USE MODI_INI_PHYEX, ONLY: INI_PHYEX
IMPLICIT NONE

! -----------------------------------------------------------------------
!     DUMMY VARIABLES
INTEGER, INTENT (IN) :: KULOUT
TYPE(CST_t),            INTENT(OUT) :: CST
!-----------------------------------------------------------------------
!    LOCAL VARIABLES
REAL :: ZDZMIN, ZTSTEP
CHARACTER(LEN=6) :: CPROGRAM
CHARACTER(LEN=4) :: CMICRO
! -----------------------------------------------------------------------

CPROGRAM='AROME'
ZDZMIN=999.
ZTSTEP=999
CMICRO='NONE'

!Default values
CALL INI_PHYEX(CPROGRAM, 0, .TRUE., KULOUT, 0, 1, &
              &ZTSTEP, ZDZMIN, &
              &CMICRO, &
              &LDDEFAULTVAL=.TRUE., LDREADNAM=.FALSE., LDCHECK=.FALSE., LDPRINT=.FALSE., LDINIT=.FALSE.)

!Param initialisation
CALL INI_PHYEX(CPROGRAM, 0, .TRUE., KULOUT, 0, 1, &
              &ZTSTEP, ZDZMIN, &
              &CMICRO, &
              &LDDEFAULTVAL=.FALSE., LDREADNAM=.FALSE., LDCHECK=.TRUE., LDPRINT=.TRUE., LDINIT=.TRUE., &
              &CST_INOUT=CST)




!
CALL INI_NEB
CALL TURB_GOTO_MODEL(1,1)
CALL CTURB_ASSOCIATE()
!CALL TBUCONF_ASSOCIATE
CALL PARAM_MFSHALL_GOTO_MODEL(1,1)
!LBU_ENABLE=.FALSE.                                                                                                       
!LBUDGET_U=.FALSE.
!LBUDGET_V=.FALSE.
!LBUDGET_W=.FALSE.
!LBUDGET_TH=.FALSE.
!LBUDGET_TKE=.FALSE.
!LBUDGET_RV=.FALSE.
!LBUDGET_RC=.FALSE.
!LBUDGET_RR=.FALSE.
!LBUDGET_RI=.FALSE.
!LBUDGET_RS=.FALSE.
!LBUDGET_RG=.FALSE.
!LBUDGET_RH=.FALSE.
!LBUDGET_SV=.FALSE.

PARAM_MFSHALLN%XALP_PERT   = 0.3
PARAM_MFSHALLN%XABUO       = 1.
PARAM_MFSHALLN%XBENTR      = 1.
PARAM_MFSHALLN%XBDETR      = 0.
PARAM_MFSHALLN%XCMF        = 0.065
PARAM_MFSHALLN%XENTR_MF    = 0.035
PARAM_MFSHALLN%XCRAD_MF    = 50.
PARAM_MFSHALLN%XENTR_DRY   = 0.55
PARAM_MFSHALLN%XDETR_DRY   = 10.
PARAM_MFSHALLN%XDETR_LUP   = 1.
PARAM_MFSHALLN%XKCF_MF     = 2.75
PARAM_MFSHALLN%XKRC_MF     = 1.
PARAM_MFSHALLN%XTAUSIGMF   = 600.
PARAM_MFSHALLN%XPRES_UV    = 0.5
PARAM_MFSHALLN%XFRAC_UP_MAX= 0.33
PARAM_MFSHALLN%XALPHA_MF = 2.
PARAM_MFSHALLN%XSIGMA_MF = 20.
PARAM_MFSHALLN%XA1    =  2.
PARAM_MFSHALLN%XB     =  0.002
PARAM_MFSHALLN%XC     =  0.012
PARAM_MFSHALLN%XBETA1 =  0.9
PARAM_MFSHALLN%XR     =  2.
PARAM_MFSHALLN%XLAMBDA_MF  = 0.
PARAM_MFSHALLN%LGZ = .FALSE.
PARAM_MFSHALLN%XGZ=1.
PARAM_MFSHALLN%CMF_UPDRAFT=HMF_UPDRAFT
PARAM_MFSHALLN%CMF_CLOUD=HMF_CLOUD
PARAM_MFSHALLN%LMIXUV=OMIXUV
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

