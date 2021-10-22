!     ######spl
SUBROUTINE CH_ORILAM(PAERO, PCHEM, PM, PSIG0, PRG0, PN0, PCTOTG, PCTOTA,&
                           PCCTOT, PDTACT, PSEDA,&
                           PMU, PLAMBDA, PRHOP0, POM, PSO4RAT,  &
                           PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PFRAC, PMI,&
                           PTIME, GSCHEME, PSOLORG)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!! #####################################################################################
!!
!!    PURPOSE
!!    -------
!!
!!    ORILAM aerosol Code
!!    
!!
!!    Inputs:
!!    PCHEM : Chemical (gaseous and aerosol) species (in molec./cm3)  
!!    PSEDA : Moments 
!!
!!    Outputs:
!!
!!
!!
!!    REFERENCE
!!    ---------
!!    P. Tulet, V. Crassier, F. Cousin, K. Suhre, R. Rosset, jgr
!!    ORILAM, A three moment lognormal aerosol scheme for mesoscale atmospheric 
!!    model.
!!    On-line coupling into the Meso-NH-C model and validation  on the Escompte 
!!    campaign.
!!
!!    AUTHOR
!!    ------
!!    Pierre Tulet (GMEI) and Vincent Crassier (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_AER_TRANS
USE MODI_CH_AER_DRIVER

!!
!!  IMPLICIT ARGUMENTS
!!  ------------------
!!
USE MODD_CH_AEROSOL
!
!-------------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
REAL,                                   INTENT(IN)    :: PDTACT, PTIME
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PRHOP0, POM
REAL,                 DIMENSION(:),     INTENT(INOUT) :: PLAMBDA, PMU, PSO4RAT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSEDA
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCHEM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PAERO
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PFRAC
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PMI
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSIG0, PRG0, PN0
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG
REAL,                 DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG
REAL,                 DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC
CHARACTER(LEN=10),                      INTENT(IN)    :: GSCHEME

REAL, DIMENSION(SIZE(PAERO,1),JPMODE)                 :: ZMASK
!
!-------------------------------------------------------------------------------
! transfer gas phase variables into aerosol variables
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_ORILAM',0,ZHOOK_HANDLE)
CALL CH_AER_TRANS(0, PM, PSIG0, PRG0, PN0, PRHOP0,PAERO, PCHEM, PCTOTG, PCTOTA, PCCTOT,&
                         PFRAC, PMI, ZMASK,GSCHEME)

! integrate aerosol variables
CALL CH_AER_DRIVER(PM,PSIG0, PRG0, PN0,  PCTOTG, PCTOTA, PCCTOT,         &
                      PDTACT, PSEDA, PMU, PLAMBDA, PRHOP0, POM, PSO4RAT, &
                      PRV, PDENAIR, PPRESSURE, PTEMP, PRC, ZMASK, PTIME, &
                      PSOLORG)
!
! transfer aerosol variables back into gas phase variables
 CALL CH_AER_TRANS(1, PM, PSIG0, PRG0, PN0, PRHOP0, PAERO, PCHEM, PCTOTG, PCTOTA, PCCTOT,&
                      PFRAC, PMI, ZMASK,GSCHEME)
!
!
IF (LHOOK) CALL DR_HOOK('CH_ORILAM',1,ZHOOK_HANDLE)
END SUBROUTINE CH_ORILAM
