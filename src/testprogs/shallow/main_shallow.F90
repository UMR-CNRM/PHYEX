PROGRAM MAIN_SHALLOW

USE XRD_GETOPTIONS
USE GETDATA_SHALLOW_MOD
USE COMPUTE_DIFF
USE MODI_SHALLOW_MF
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_PHYEX,      ONLY: PHYEX_t
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
TYPE(PHYEX_t)            :: PHYEX
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
                  &PFRAC_UP_OUT, PEMF_OUT, ZDETR_OUT, ZENTR_OUT, IKLCL_OUT, IKETL_OUT, IKCTL_OUT, LLVERBOSE)

IF (LLVERBOSE) PRINT *, " KLEV = ", KLEV, " KRR = ", KRR

PRINT *, " NPROMA = ", NPROMA, " KLEV = ", KLEV, " NGPBLKS = ", NGPBLKS

CALL INIT_PHYEX(KRR, KRRL, KRRI, KSV, PHYEX)

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
!!!              !directives pas a jour !$acc      & copyin  (D0, CST, ICEP, NEBN, KRR, HCONDENS, HLAMBDA3, HBUNAME, OSIGMAS, OCND2, PTSTEP, LMFCONV, &
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
  CALL SHALLOW_MF(D, PHYEX%CST, PHYEX%NEBN, PHYEX%PARAM_MFSHALLN, PHYEX%TURBN, PHYEX%CSTURB,                    &
     &KRR=PHYEX%MISC%KRR, KRRL=PHYEX%MISC%KRRL, KRRI=PHYEX%MISC%KRRI, KSV=PHYEX%MISC%KSV,                                             &
     &ONOMIXLG=PHYEX%MISC%ONOMIXLG,KSV_LGBEG=PHYEX%MISC%KSV_LGBEG,KSV_LGEND=PHYEX%MISC%KSV_LGEND,      &
     &PTSTEP=PHYEX%MISC%PTSTEP, &
     &PDZZ=PDZZF(:,:,:,IBL),PZZ=PZZ(:,:,:,IBL),                                                                 &
     &PRHODJ=PRHODJ(:,:,:,IBL),PRHODREF=PRHODREF(:,:,:,IBL),                                                    &
     &PPABSM=PPABSM(:,:,:,IBL),PEXNM=PEXNM(:,:,:,IBL),                                                          &
     &PSFTH=PSFTH(JLON,1,IBL),PSFRV=PSFRV(JLON,1,IBL),                                                            &
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
     &KKLCL=IKLCL(JLON,1,IBL),KKETL=IKETL(JLON,1,IBL),KKCTL=IKCTL(JLON,1,IBL),PDX=PHYEX%MISC%PDX,PDY=PHYEX%MISC%PDY,KBUDGETS=PHYEX%MISC%NBUDGET )

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

SUBROUTINE INIT_PHYEX(KRR, KRRL, KRRI, KSV, PHYEX)

USE MODD_BUDGET, ONLY: TBUCONF_ASSOCIATE, NBUDGET_RI, TBUCONF, LBU_ENABLE, LBUDGET_U, LBUDGET_V, LBUDGET_W, LBUDGET_TH, &
                       LBUDGET_TKE, LBUDGET_RV, LBUDGET_RC, LBUDGET_RR, LBUDGET_RI, LBUDGET_RS, LBUDGET_RG, LBUDGET_RH, LBUDGET_SV
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODI_INI_PHYEX, ONLY: INI_PHYEX

IMPLICIT NONE

! -----------------------------------------------------------------------
!     DUMMY VARIABLES
INTEGER,          INTENT(IN)  :: KRR, KRRL, KRRI, KSV
TYPE(PHYEX_t),    INTENT(OUT) :: PHYEX

!-----------------------------------------------------------------------
!    LOCAL VARIABLES
INTEGER :: IULOUT, JRR
REAL :: ZDZMIN
CHARACTER(LEN=6) :: CPROGRAM
CHARACTER(LEN=4) :: CMICRO, CSCONV, CTURB
REAL             :: PTSTEP
! -----------------------------------------------------------------------

IULOUT=20
CPROGRAM='AROME'
ZDZMIN=999.
CMICRO='NONE'
CSCONV='EDKF'
CTURB='TKEL'
PTSTEP = 25.0000000000000

!Default values
CALL INI_PHYEX(CPROGRAM, 0, .TRUE., IULOUT, 0, 1, &
              &PTSTEP, ZDZMIN, &
              &CMICRO, CSCONV, CTURB, &
              &LDDEFAULTVAL=.TRUE., LDREADNAM=.FALSE., LDCHECK=.FALSE., KPRINT=0, LDINIT=.FALSE., &
              &PHYEX_OUT=PHYEX)

!Control parameters
PHYEX%MISC%PTSTEP       = PTSTEP
PHYEX%MISC%KSV_LGBEG = 0
PHYEX%MISC%KSV_LGEND = 0
PHYEX%MISC%ONOMIXLG=.FALSE.
PHYEX%MISC%KRR          = KRR
PHYEX%MISC%KRRL         = KRRL
PHYEX%MISC%KRRI         = KRRI
PHYEX%MISC%KSV          = KSV

!Emulate the namelist reading
PHYEX%NEBN%LSUBG_COND=.TRUE.
PHYEX%NEBN%CFRAC_ICE_SHALLOW_MF='S'

!Param initialisation
CALL INI_PHYEX(CPROGRAM, 0, .TRUE., IULOUT, 0, 1, &
              &PTSTEP, ZDZMIN, &
              &CMICRO, CSCONV, CTURB, &
              &LDDEFAULTVAL=.FALSE., LDREADNAM=.FALSE., LDCHECK=.TRUE., KPRINT=2, LDINIT=.TRUE., &
              &PHYEX_IN=PHYEX, PHYEX_OUT=PHYEX)

!Budgets
CALL TBUCONF_ASSOCIATE
PHYEX%MISC%NBUDGET=0
DO JRR=1, PHYEX%MISC%NBUDGET
  PHYEX%MISC%YLBUDGET(JRR)%NBUDGET=JRR
ENDDO
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
PHYEX%MISC%TBUCONF=TBUCONF

END SUBROUTINE INIT_PHYEX

END PROGRAM

