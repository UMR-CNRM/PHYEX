PROGRAM MAIN_ICE_ADJUST

USE XRD_GETOPTIONS
USE GETDATA_ICE_ADJUST_MOD
USE COMPUTE_DIFF
USE MODI_ICE_ADJUST
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_CST,        ONLY: CST_t
USE MODD_NEB,        ONLY: NEB
USE MODD_TURB_n,     ONLY: TURB_t
USE MODD_RAIN_ICE_PARAM, ONLY : RAIN_ICE_PARAM_t
USE MODD_PARAM_ICE,      ONLY : PARAM_ICE_t
USE MODI_INI_NEB
USE MODD_BUDGET !, ONLY: TBUCONF_ASSOCIATE, TBUDGETDATA, NBUDGET_RI, TBUCONF
USE STACK_MOD
USE OMP_LIB
USE YOMHOOK, ONLY : LHOOK, DR_HOOK
USE PARKIND1, ONLY : JPRB, JPIM


IMPLICIT NONE

INTEGER      :: KLON 
INTEGER      :: KLEV
INTEGER      :: KRR

REAL,    ALLOCATABLE   :: PRHODJ         (:,:,:,:)   
REAL,    ALLOCATABLE   :: PEXNREF        (:,:,:,:)   
REAL,    ALLOCATABLE   :: PRHODREF       (:,:,:,:)   
REAL,    ALLOCATABLE   :: PPABSM         (:,:,:,:)   
REAL,    ALLOCATABLE   :: PTHT           (:,:,:,:)   
REAL,    ALLOCATABLE   :: PSIGS          (:,:,:,:)   
REAL,    ALLOCATABLE   :: PMFCONV        (:,:,:,:)   
REAL,    ALLOCATABLE   :: PRC_MF         (:,:,:,:)   
REAL,    ALLOCATABLE   :: PRI_MF         (:,:,:,:)   
REAL,    ALLOCATABLE   :: PCF_MF         (:,:,:,:)   
REAL,    ALLOCATABLE   :: PTHS           (:,:,:,:)   
REAL,    ALLOCATABLE   :: PRS            (:,:,:,:,:) 
REAL,    ALLOCATABLE   :: PSRCS          (:,:,:,:)   
REAL,    ALLOCATABLE   :: PCLDFR         (:,:,:,:)   
REAL,    ALLOCATABLE   :: PHLC_HRC       (:,:,:,:)   
REAL,    ALLOCATABLE   :: PHLC_HCF       (:,:,:,:)   
REAL,    ALLOCATABLE   :: PHLI_HRI       (:,:,:,:)   
REAL,    ALLOCATABLE   :: PHLI_HCF       (:,:,:,:)   
REAL,    ALLOCATABLE   :: ZRS            (:,:,:,:,:) 
REAL,    ALLOCATABLE   :: ZZZ            (:,:,:,:)   
REAL,    ALLOCATABLE   :: ZSIGQSAT       (:,:,:)   
REAL,    ALLOCATABLE   :: ZICE_CLD_WGT   (:,:,:)   
REAL,    ALLOCATABLE   :: ZDUM1          (:,:,:,:)
REAL,    ALLOCATABLE   :: ZDUM2          (:,:,:,:)
REAL,    ALLOCATABLE   :: ZDUM3          (:,:,:,:)
REAL,    ALLOCATABLE   :: ZDUM4          (:,:,:,:)
REAL,    ALLOCATABLE   :: ZDUM5          (:,:,:,:)

REAL,    ALLOCATABLE   :: PRS_OUT        (:,:,:,:,:) 
REAL,    ALLOCATABLE   :: PSRCS_OUT      (:,:,:,:)   
REAL,    ALLOCATABLE   :: PCLDFR_OUT     (:,:,:,:)   
REAL,    ALLOCATABLE   :: PHLC_HRC_OUT   (:,:,:,:)   
REAL,    ALLOCATABLE   :: PHLC_HCF_OUT   (:,:,:,:)   
REAL,    ALLOCATABLE   :: PHLI_HRI_OUT   (:,:,:,:)   
REAL,    ALLOCATABLE   :: PHLI_HCF_OUT   (:,:,:,:)   

INTEGER :: NPROMA, NGPBLKS, NFLEVG
INTEGER :: IBL, JLON, JLEV

TYPE(DIMPHYEX_t)         :: D, D0
TYPE(CST_t)              :: CST
TYPE(PARAM_ICE_t)        :: PARAMI
TYPE(RAIN_ICE_PARAM_t)   :: ICEP
TYPE(TURB_t)             :: TURB
CHARACTER(LEN=4)         :: HBUNAME  
LOGICAL                  :: LHGT_QS
LOGICAL                  :: LMFCONV
CHARACTER (LEN=4)   :: CMICRO
REAL                :: PTSTEP
TYPE(TBUDGETDATA), DIMENSION(NBUDGET_RI) :: YLBUDGET
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

CALL INITOPTIONS ()
NGPBLKS = 296
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

CALL GETDATA_ICE_ADJUST (NPROMA, NGPBLKS, NFLEVG, PRHODJ, PEXNREF, PRHODREF, PPABSM, PTHT, ZICE_CLD_WGT,     &
& ZSIGQSAT, PSIGS, PMFCONV, PRC_MF, PRI_MF, PCF_MF, ZDUM1, ZDUM2, ZDUM3, ZDUM4, ZDUM5, PTHS, PRS, PSRCS, PCLDFR, PHLC_HRC, PHLC_HCF, &
& PHLI_HRI, PHLI_HCF, ZRS, ZZZ, PRS_OUT, PSRCS_OUT, PCLDFR_OUT, PHLC_HRC_OUT, PHLC_HCF_OUT,       &
& PHLI_HRI_OUT, PHLI_HCF_OUT, LLVERBOSE)

KLEV = SIZE (PRS, 3)
KRR  = SIZE (PRS, 4)

IF (LLVERBOSE) PRINT *, " KLEV = ", KLEV, " KRR = ", KRR

PRINT *, " NPROMA = ", NPROMA, " KLEV = ", KLEV, " NGPBLKS = ", NGPBLKS

CMICRO='ICE3'
PTSTEP       = 50.000000000000000
HBUNAME      = 'DEPI'
LMFCONV      = .TRUE.
LHGT_QS      = .FALSE.
CALL INIT_PHYEX (20, CMICRO, PTSTEP, &
                 CST, &
                 PARAMI, ICEP, TURB)

DO JRR=1, NBUDGET_RI
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
D0%NKB  = KLEV 
D0%NKE  = 1
D0%NKTB = 1
D0%NKTE = KLEV

ISTSZ = NPROMA * 20 * KLEV
ALLOCATE (PSTACK (ISTSZ, NGPBLKS))

TS = OMP_GET_WTIME ()

ZTD = 0.
ZTC = 0.

IF (LHOOK) CALL DR_HOOK ('MAIN',0,ZHOOK_HANDLE)

DO ITIME = 1, NTIME

  TSD = OMP_GET_WTIME ()

!$acc data &
!$acc      & copyin  (D0, CST, ICEP, NEB, TURB, KRR, HCONDENS, HLAMBDA3, HBUNAME, OSUBG_COND, OSIGMAS, LHGT_QS, HSUBG_MF_PDF, PTSTEP, LMFCONV, &
!$acc      &          ZSIGQSAT, PRHODJ, PEXNREF, PRHODREF, PSIGS, PMFCONV, PPABSM, ZZZ, PCF_MF, PRC_MF, PRI_MF, ZDUM1, ZDUM2, ZDUM3, ZDUM4, ZDUM5, ZRS, ZICE_CLD_WGT) &
!$acc      & copy    (PRS, PTHS), &
!$acc      & copyout (PSRCS, PCLDFR, PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF) &
!$acc      & create  (PSTACK) 

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

    CALL ICE_ADJUST (D, CST, ICEP, NEB, TURB, TBUCONF, KRR, PARAMI%CFRAC_ICE_ADJUST, HBUNAME,                                &
    & PARAMI%LOCND2, LHGT_QS, PTSTEP, ZSIGQSAT (:, :, IBL), PRHODJ=PRHODJ (:, :, :, IBL),                                               &
    & PEXNREF=PEXNREF (:, :, :, IBL),                                                                                           &
    & PRHODREF=PRHODREF (:, :, :, IBL), PSIGS=PSIGS (:, :, :, IBL), LMFCONV=LMFCONV, PMFCONV=PMFCONV (:, :, :, IBL),            &
    & PPABST=PPABSM (:, :, :, IBL), PZZ=ZZZ (:, :, :, IBL), PEXN=PEXNREF (:, :, :, IBL), PCF_MF=PCF_MF (:, :, :, IBL),          &
    & PRC_MF=PRC_MF (:, :, :, IBL), PRI_MF=PRI_MF  (:, :, :, IBL),                                                              &
    & PICLDFR=ZDUM1(:, :, :, IBL), PWCLDFR=ZDUM2(:, :, :, IBL), PSSIO=ZDUM3(:, :, :, IBL),                                      &
    & PSSIU=ZDUM4(:, :, :, IBL), PIFR=ZDUM5(:, :, :, IBL),                                                                      &
    & PRV=ZRS(:, :, :, 1, IBL), PRC=ZRS(:, :, :, 2, IBL),                                                                       &
    & PRVS=PRS(:, :, :, 1, IBL), PRCS=PRS(:, :, :, 2, IBL), PTH=ZRS(:, :, :, 0, IBL), PTHS=PTHS (:, :, :, IBL),                 &
    & OCOMPUTE_SRC=.TRUE.,                                                                                                      &
    & PSRCS=PSRCS (:, :, :, IBL), PCLDFR=PCLDFR (:, :, :, IBL), PRR=ZRS(:, :, :, 3, IBL), PRI=ZRS(:, :, :, 4, IBL),             &
    & PRIS=PRS(:, :, :, 4, IBL), PRS=ZRS(:, :, :, 5, IBL), PRG=ZRS(:, :, :, 6, IBL), TBUDGETS=YLBUDGET, KBUDGETS=NBUDGET_RI,    &
    & PICE_CLD_WGT=ZICE_CLD_WGT(:, :, IBL),                                                                                     &
    & PHLC_HRC=PHLC_HRC(:, :, :, IBL), PHLC_HCF=PHLC_HCF(:, :, :, IBL),                                                         &
    & PHLI_HRI=PHLI_HRI(:, :, :, IBL), PHLI_HCF=PHLI_HCF(:, :, :, IBL)                                                          &
#ifdef USE_STACK
    & , YDSTACK=YLSTACK &
#endif
    & )

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
    CALL DIFF ("PSRCS",    PSRCS_OUT    (:,:,:,IBL), PSRCS    (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF ("PCLDFR",   PCLDFR_OUT   (:,:,:,IBL), PCLDFR   (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF ("PHLC_HRC", PHLC_HRC_OUT (:,:,:,IBL), PHLC_HRC (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF ("PHLC_HCF", PHLC_HCF_OUT (:,:,:,IBL), PHLC_HCF (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF ("PHLI_HRI", PHLI_HRI_OUT (:,:,:,IBL), PHLI_HRI (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
    CALL DIFF ("PHLI_HCF", PHLI_HCF_OUT (:,:,:,IBL), PHLI_HCF (:,:,:,IBL), LLSTAT, LLCHECK, NPROMA, LLCHECKDIFF, LLDIFF)
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

SUBROUTINE INIT_PHYEX(KULOUT,CMICRO,PTSTEP, &
                      CST, &
                      PARAM_ICE, RAIN_ICE_PARAM, TURB)

USE MODD_CST, ONLY: CST_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE MODD_PARAM_ICE, ONLY: PARAM_ICE_t
USE MODD_TURB_N, ONLY: TURB_GOTO_MODEL, TURB_t, TURBN
USE MODI_INI_PHYEX, ONLY: INI_PHYEX

IMPLICIT NONE

! -----------------------------------------------------------------------
!     DUMMY VARIABLES
INTEGER, INTENT (IN) :: KULOUT
CHARACTER(4), INTENT (IN) :: CMICRO 
REAL, INTENT(IN) :: PTSTEP
TYPE(CST_t),            INTENT(OUT) :: CST
TYPE(PARAM_ICE_t)     , INTENT(OUT) :: PARAM_ICE
TYPE(RAIN_ICE_PARAM_t), INTENT(OUT) :: RAIN_ICE_PARAM
TYPE(TURB_t),           INTENT(OUT) :: TURB

!-----------------------------------------------------------------------
!    LOCAL VARIABLES
REAL :: ZDZMIN
CHARACTER(LEN=6) :: CPROGRAM
! -----------------------------------------------------------------------

CPROGRAM='AROME'
ZDZMIN=20.

!Default values
CALL INI_PHYEX(CPROGRAM, 0, .TRUE., KULOUT, 0, 1, &
              &PTSTEP, ZDZMIN, &
              &CMICRO, &
              &LDDEFAULTVAL=.TRUE., LDREADNAM=.FALSE., LDCHECK=.FALSE., LDPRINT=.FALSE., LDINIT=.FALSE., &
              &PARAM_ICE_INOUT=PARAM_ICE)

!Emulate the namelist reading
PARAM_ICE%LCRIAUTI=.TRUE.
PARAM_ICE%XCRIAUTI_NAM=0.2E-3
PARAM_ICE%XT0CRIAUTI_NAM=-5.
PARAM_ICE%XCRIAUTC_NAM=0.1E-2
PARAM_ICE%LOCND2=.FALSE.
PARAM_ICE%CSEDIM='STAT'
PARAM_ICE%LWARM=.TRUE.
PARAM_ICE%LSEDIC=.TRUE.
PARAM_ICE%CFRAC_ICE_ADJUST='S' ! Ice/liquid partition rule to use in adjustment
PARAM_ICE%CFRAC_ICE_SHALLOW_MF='S' ! Ice/liquid partition rule to use in shallow_mf
PARAM_ICE%CSNOWRIMING='M90 '
PARAM_ICE%XFRACM90=0.1 ! Fraction used for the Murakami 1990 formulation
PARAM_ICE%LCONVHG=.TRUE. ! TRUE to allow the conversion from hail to graupel
PARAM_ICE%LCRFLIMIT=.TRUE. !True to limit rain contact freezing to possible heat exchange
PARAM_ICE%LFEEDBACKT=.TRUE. ! When .TRUE. feed back on temperature is taken into account
PARAM_ICE%LEVLIMIT=.TRUE.   ! When .TRUE. water vapour pressure is limited by saturation
PARAM_ICE%LNULLWETG=.TRUE.  ! When .TRUE. graupel wet growth is activated with null rate (to allow water shedding)
PARAM_ICE%LWETGPOST=.TRUE.  ! When .TRUE. graupel wet growth is activated with positive temperature (to allow water shedding)
PARAM_ICE%LNULLWETH=.TRUE.  ! Same as LNULLWETG but for hail
PARAM_ICE%LWETHPOST=.TRUE.  ! Same as LWETGPOST but for hail
PARAM_ICE%LSEDIM_AFTER=.FALSE. ! Sedimentation done after microphysics
PARAM_ICE%XSPLIT_MAXCFL=0.8
PARAM_ICE%LDEPOSC=.FALSE.  ! water deposition on vegetation
PARAM_ICE%XVDEPOSC=0.02    ! deposition speed (2 cm.s-1)
PARAM_ICE%CSUBG_RC_RR_ACCR='NONE'
PARAM_ICE%CSUBG_RR_EVAP='NONE'
PARAM_ICE%CSUBG_PR_PDF='SIGM'

!Param initialisation
CALL INI_PHYEX(CPROGRAM, 0, .TRUE., KULOUT, 0, 1, &
              &PTSTEP, ZDZMIN, &
              &CMICRO, &
              &LDDEFAULTVAL=.FALSE., LDREADNAM=.FALSE., LDCHECK=.TRUE., LDPRINT=.TRUE., LDINIT=.TRUE., &
              &CST_INOUT=CST, &
              &PARAM_ICE_INOUT=PARAM_ICE, RAIN_ICE_PARAM_INOUT=RAIN_ICE_PARAM)

CALL TURB_GOTO_MODEL(1,1)
TURBN%LSTATNW=.FALSE.
TURBN%CSUBG_MF_PDF='TRIANGLE'
TURBN%LSTATNW=.FALSE.
TURBN%CCONDENS     = 'CB02'
TURBN%CLAMBDA3     = 'CB'
TURBN%LSUBG_COND   = .TRUE.
TURBN%LSIGMAS      = .TRUE.
TURBN%CSUBG_MF_PDF = 'TRIANGLE'
TURB=TURBN
CALL INI_NEB

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

END SUBROUTINE INIT_PHYEX

END PROGRAM

