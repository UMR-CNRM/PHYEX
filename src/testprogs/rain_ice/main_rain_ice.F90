PROGRAM MAIN_RAIN_ICE

USE MODD_CONF
USE MODD_CST, ONLY: CST
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM
USE MODD_PARAM_ICE,      ONLY: PARAM_ICE
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t

USE MODD_BUDGET, ONLY: TBUDGETDATA, NBUDGET_RH, TBUCONF

USE MODI_RAIN_ICE

IMPLICIT NONE

INTEGER   :: KPROMA  
INTEGER   :: KKA  
INTEGER   :: KKU  
INTEGER   :: KKL  
INTEGER   :: KLON    
INTEGER   :: KLEV     
INTEGER   :: KRR      
INTEGER   :: KDUM

CHARACTER (LEN=4)   :: CSUBG_AUCV_RC
CHARACTER (LEN=80)  :: CSUBG_AUCV_RI
LOGICAL             :: OSEDIC 
CHARACTER (LEN=4)   :: CSEDIM  
CHARACTER (LEN=4)   :: CMICRO  
REAL                :: PTSTEP 


REAL,    ALLOCATABLE, DIMENSION(:,:,:) :: PRS, PRS_1
REAL,    ALLOCATABLE, DIMENSION(:,:,:) :: PFPR, PFPR_1
REAL,    ALLOCATABLE, DIMENSION(:,:,:) :: PRT   
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PDZZ     
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PRHODJ  
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PRHODREF
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PEXNREF 
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PPABSM  
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PHLC_HRC
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PHLC_HCF
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PHLI_HRI
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PHLI_HCF
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PTHT    
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PSIGS   
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PCLDFR  
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PTHS    
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PEVAP, PEVAP_1
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: PCIT, PCIT_1
REAL,    ALLOCATABLE, DIMENSION(:)     :: PSEA  
REAL,    ALLOCATABLE, DIMENSION(:)     :: PTOWN  
REAL,    ALLOCATABLE, DIMENSION(:)     :: PINPRR, PINPRR_1
REAL,    ALLOCATABLE, DIMENSION(:)     :: PINPRS, PINPRS_1
REAL,    ALLOCATABLE, DIMENSION(:)     :: PINPRG, PINPRG_1
REAL,    ALLOCATABLE, DIMENSION(:)     :: ZINDEP, ZINDEP_1
REAL,    ALLOCATABLE, DIMENSION(:,:)   :: ZRAINFR, ZRAINFR_1
REAL,    ALLOCATABLE, DIMENSION(:)     :: ZINPRC, ZINPRC_1
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LLMICRO 

LOGICAL :: OWARM 
LOGICAL :: OCND2 

INTEGER :: IFILE

LOGICAL :: LCRIAUTI
REAL :: ZCRIAUTI, ZT0CRIAUTI, ZCRIAUTC

TYPE(TBUDGETDATA), DIMENSION(NBUDGET_RH) :: YLBUDGET 
TYPE(DIMPHYEX_t) :: D

INTEGER :: IPROMA, ISIZE
INTEGER :: JLON, JLEV

CHARACTER(LEN=128) :: CLFILE



CMICRO='ICE3'

CALL GETARG (1, CLFILE)

PTSTEP = 25.0000000000000
KRR = 6
OSEDIC = .TRUE.
OCND2 = .FALSE.
CSEDIM = 'STAT'
CSUBG_AUCV_RC = 'PDF'
CSUBG_AUCV_RI = 'NONE'
OWARM = .TRUE.

LCRIAUTI=.TRUE.
ZCRIAUTI=0.2E-3
ZT0CRIAUTI=-5.
ZCRIAUTC=0.1E-2

CALL INIT_PHYEX (20, OWARM, CMICRO, CSEDIM, &
            & LCRIAUTI, ZCRIAUTI, ZT0CRIAUTI, ZCRIAUTC)

IFILE = 77
OPEN (IFILE, FILE=TRIM (CLFILE), FORM='UNFORMATTED') 
READ (IFILE) IPROMA, ISIZE
READ (IFILE) KLON, KDUM, KLEV, KRR

PRINT *, KLON, KDUM, KLEV, KRR

ALLOCATE (PRT       (KLON,KLEV,KRR))
ALLOCATE (PRS       (KLON,KLEV,KRR), PRS_1  (KLON,KLEV,KRR))
ALLOCATE (PFPR      (KLON,KLEV,KRR), PFPR_1 (KLON,KLEV,KRR))
ALLOCATE (PDZZ      (KLON,KLEV))
ALLOCATE (PRHODJ    (KLON,KLEV))
ALLOCATE (PRHODREF  (KLON,KLEV))
ALLOCATE (PEXNREF   (KLON,KLEV))
ALLOCATE (PPABSM    (KLON,KLEV))
ALLOCATE (PHLC_HRC  (KLON,KLEV))
ALLOCATE (PHLC_HCF  (KLON,KLEV))
ALLOCATE (PHLI_HRI  (KLON,KLEV))
ALLOCATE (PHLI_HCF  (KLON,KLEV))
ALLOCATE (PTHT      (KLON,KLEV))
ALLOCATE (PSIGS     (KLON,KLEV))
ALLOCATE (PCLDFR    (KLON,KLEV))
ALLOCATE (PTHS      (KLON,KLEV))
ALLOCATE (PEVAP     (KLON,KLEV),     PEVAP_1 (KLON,KLEV))
ALLOCATE (PCIT      (KLON,KLEV),     PCIT_1 (KLON,KLEV))
ALLOCATE (PSEA      (KLON))
ALLOCATE (PTOWN     (KLON))
ALLOCATE (PINPRR    (KLON))
ALLOCATE (PINPRS    (KLON),          PINPRS_1 (KLON))
ALLOCATE (PINPRG    (KLON),          PINPRG_1 (KLON))
ALLOCATE (LLMICRO   (KLON,KLEV))
ALLOCATE (ZINDEP    (KLON),          ZINDEP_1 (KLON))
ALLOCATE (ZRAINFR   (KLON,KLEV),     ZRAINFR_1 (KLON,KLEV))
ALLOCATE (ZINPRC    (KLON),          ZINPRC_1 (KLON))

READ (IFILE) LLMICRO
READ (IFILE) PEXNREF
READ (IFILE) PDZZ
READ (IFILE) PRHODJ
READ (IFILE) PRHODREF
READ (IFILE) PEXNREF
READ (IFILE) PPABSM
READ (IFILE) PCIT
READ (IFILE) PCLDFR
READ (IFILE) PHLC_HRC
READ (IFILE) PHLC_HCF
READ (IFILE) PHLI_HRI
READ (IFILE) PHLI_HCF
READ (IFILE) PTHT
READ (IFILE) PRT
READ (IFILE) PTHS
READ (IFILE) PRS
READ (IFILE) PSIGS
READ (IFILE) PSEA
READ (IFILE) PTOWN

READ (IFILE) PCIT_1
READ (IFILE) PRS_1
READ (IFILE) ZINPRC_1
READ (IFILE) PINPRR_1
READ (IFILE) PEVAP_1
READ (IFILE) PINPRS_1
READ (IFILE) PINPRG_1
READ (IFILE) ZINDEP_1
READ (IFILE) ZRAINFR_1
READ (IFILE) PFPR_1
CLOSE (IFILE)

D%NIT  = KLON
D%NIB  = 1                                                                                                                         
D%NIE  = KLON
D%NJT  = 1                                                                                                                         
D%NJB  = 1                                                                                                                         
D%NJE  = 1                                                                                                                         
D%NKL  = -1                                                                                                                        
D%NKT  = KLEV                                                                                                                      
D%NKA  = KLEV                                                                                                                      
D%NKU  = 1                                                                                                                         
D%NKB  = KLEV                                                                                                                      
D%NKE  = 1                                                                                                                         
D%NKTB = 1                                                                                                                         
D%NKTE = KLEV   

CALL RAIN_ICE (D, CST, PARAM_ICE, RAIN_ICE_PARAM, &
             & RAIN_ICE_DESCR, TBUCONF, &
             & IPROMA, ISIZE, &
             & OSEDIC=OSEDIC, OCND2=OCND2, HSEDIM=CSEDIM, &
             & HSUBG_AUCV_RC=CSUBG_AUCV_RC, HSUBG_AUCV_RI=CSUBG_AUCV_RI,&
             & OWARM=OWARM, &
             & PTSTEP=2*PTSTEP, &
             & KRR=KRR, ODMICRO=LLMICRO, PEXN=PEXNREF,            &
             & PDZZ=PDZZ, PRHODJ=PRHODJ, PRHODREF=PRHODREF,PEXNREF=PEXNREF,&
             & PPABST=PPABSM, PCIT=PCIT, PCLDFR=PCLDFR,  &
             & PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF, &
             & PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF, &
             & PTHT=PTHT,PRVT=PRT(1,1,1),PRCT=PRT(1,1,2), &
             & PRRT=PRT(1,1,3), &
             & PRIT=PRT(1,1,4), PRST=PRT(1,1,5), &
             & PRGT=PRT(1,1,6),       &
             & PTHS=PTHS, PRVS=PRS(1,1,1),PRCS=PRS(1,1,2),&
             & PRRS=PRS(1,1,3),&
             & PRIS=PRS(1,1,4),PRSS= PRS(1,1,5),PRGS= PRS(1,1,6),&
             & PINPRC=ZINPRC,PINPRR=PINPRR,PEVAP3D=PEVAP,&
             & PINPRS=PINPRS, PINPRG=PINPRG, PINDEP=ZINDEP, PRAINFR=ZRAINFR, &
             & PSIGS=PSIGS, &
             & TBUDGETS=YLBUDGET, KBUDGETS=SIZE(YLBUDGET), &
             & PSEA=PSEA, PTOWN=PTOWN, PFPR=PFPR)

DO JLEV = 1, KLEV
  WRITE (*, '(I4," | ")', ADVANCE='NO') JLEV
  DO JLON = 1, KLON
    WRITE (*, '(E12.5)', ADVANCE='NO') ZRAINFR_1 (JLON, JLEV)
  ENDDO
  WRITE (*, *)
ENDDO


CONTAINS

SUBROUTINE INIT_PHYEX(KULOUT,LDWARM,CMICRO,CCSEDIM,LDCRIAUTI,&
                   PCRIAUTI,PT0CRIAUTI,PCRIAUTC)

USE MODD_RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM
USE MODD_PARAM_ICE
USE MODD_TURB_N, ONLY: TURB_GOTO_MODEL, CSUBG_MF_PDF

USE MODD_REF
USE MODI_INI_RAIN_ICE

IMPLICIT NONE
! -----------------------------------------------------------------------
!     DUMMY INTEGER SCALARS
INTEGER, INTENT (IN) :: KULOUT
LOGICAL, INTENT (IN) :: LDWARM
CHARACTER(4), INTENT (IN) :: CMICRO 
CHARACTER(4), INTENT (IN) :: CCSEDIM
LOGICAL, INTENT (IN) :: LDCRIAUTI
REAL, INTENT (IN) :: PCRIAUTI
REAL, INTENT (IN) :: PT0CRIAUTI
REAL, INTENT (IN) :: PCRIAUTC
!-----------------------------------------------------------------------
!    LOCAL VARIABLES
REAL :: ZCRI0, ZTCRI0    
! -----------------------------------------------------------------------

CALL INI_CST
CALL TURB_GOTO_MODEL(1,1)
CALL PARAM_ICE_ASSOCIATE

!        1. Set implicit default values for MODD_PARAM_ICE
LWARM=LDWARM
CPRISTINE_ICE='PLAT'
CSEDIM=CCSEDIM
CSUBG_AUCV_RC='PDF'
CSUBG_AUCV_RI='NONE'
CSUBG_RC_RR_ACCR='NONE'
CSUBG_RR_EVAP='NONE'
CSUBG_PR_PDF='SIGM'
CSUBG_MF_PDF='TRIANGLE'
! Snow riming                                                                                                                       
CSNOWRIMING='M90 '
XFRACM90=0.1 ! Fraction used for the Murakami 1990 formulation
!                                                                                                                                   
LFEEDBACKT=.TRUE. ! When .TRUE. feed back on temperature is taken into account
LEVLIMIT=.TRUE.   ! When .TRUE. water vapour pressure is limited by saturation
LNULLWETG=.TRUE.  ! When .TRUE. graupel wet growth is activated with null rate (to allow water shedding)
LWETGPOST=.TRUE.  ! When .TRUE. graupel wet growth is activated with positive temperature (to allow water shedding)
LNULLWETH=.TRUE.  ! Same as LNULLWETG but for hail
LWETHPOST=.TRUE.  ! Same as LWETGPOST but for hail
LCONVHG=.TRUE. ! TRUE to allow the conversion from hail to graupel
LCRFLIMIT=.TRUE. !True to limit rain contact freezing to possible heat exchange
CFRAC_ICE_ADJUST='S' ! Ice/liquid partition rule to use in adjustment
CFRAC_ICE_SHALLOW_MF='S' ! Ice/liquid partition rule to use in shallow_mf
LSEDIM_AFTER=.FALSE. ! Sedimentation done after microphysics
XSPLIT_MAXCFL=0.8
LDEPOSC=.FALSE.  ! water deposition on vegetation
XVDEPOSC=0.02    ! deposition speed (2 cm.s-1)
!
!        2. Set implicit default values for MODD_RAIN_ICE_DESCR 
!                     et MODD_RAIN_ICE_PARAM
XTHVREFZ=300.
!
CALL INI_RAIN_ICE (KULOUT, CMICRO)
!update values from namparar
IF (LDCRIAUTI) THEN

  XCRIAUTI=PCRIAUTI
  XCRIAUTC=PCRIAUTC
  XT0CRIAUTI=PT0CRIAUTI
  !second point to determine 10**(aT+b) law
  ZTCRI0=-40.0
  ZCRI0=1.25E-6
  
  XBCRIAUTI=-( LOG10(XCRIAUTI) - LOG10(ZCRI0)*PT0CRIAUTI/ZTCRI0 )&
                   *ZTCRI0/(XT0CRIAUTI-ZTCRI0)
  XACRIAUTI=(LOG10(ZCRI0)-XBCRIAUTI)/ZTCRI0
  
ENDIF
! -----------------------------------------------------------------------

END SUBROUTINE INIT_PHYEX

END PROGRAM

