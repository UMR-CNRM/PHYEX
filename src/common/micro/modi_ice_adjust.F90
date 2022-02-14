!     ######################
      MODULE MODI_ICE_ADJUST
!     ######################
!
INTERFACE
!
      SUBROUTINE ICE_ADJUST (KKA, KKU, KKL, KRR, HFRAC_ICE, HCONDENS, HLAMBDA3,&
                             HBUNAME, OSUBG_COND, OSIGMAS, OCND2, HSUBG_MF_PDF,&
                             PTSTEP, PSIGQSAT,                                 &
                             PRHODJ, PEXNREF, PRHODREF, PSIGS, PMFCONV,        &
                             PPABST, PZZ,                                      &
                             PEXN, PCF_MF, PRC_MF, PRI_MF,                     &
                             PRV, PRC, PRVS, PRCS, PTH, PTHS, PSRCS, PCLDFR,   &
                             PRR, PRI, PRIS, PRS, PRG, PRH,                    &
                             POUT_RV, POUT_RC, POUT_RI, POUT_TH,               &
                             PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,           &
                             TBUDGETS, KBUDGETS,                               &
                             PICE_CLD_WGT)
USE MODD_BUDGET,         ONLY: TBUDGETDATA
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN)    :: KKA  !near ground array index
INTEGER,                  INTENT(IN)    :: KKU  !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL  !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER(LEN=1),         INTENT(IN)    :: HFRAC_ICE
CHARACTER(LEN=80),        INTENT(IN)    :: HCONDENS
CHARACTER(LEN=4),         INTENT(IN)    :: HLAMBDA3 ! formulation for lambda3 coeff
CHARACTER(LEN=4),         INTENT(IN)    :: HBUNAME  ! Name of the budget
LOGICAL,                  INTENT(IN)    :: OSUBG_COND ! Switch for Subgrid
                                                    ! Condensation
LOGICAL                                 :: OSIGMAS  ! Switch for Sigma_s:
                                                    ! use values computed in CONDENSATION
                                                    ! or that from turbulence scheme
LOGICAL                                 :: OCND2    ! logical switch to sparate liquid
                                                    ! and ice
                                                    ! more rigid (DEFALT value : .FALSE.)
CHARACTER(LEN=80),        INTENT(IN)    :: HSUBG_MF_PDF
REAL,                     INTENT(IN)   :: PTSTEP    ! Double Time step
                                                    ! (single if cold start)
REAL, DIMENSION(:,:),     INTENT(IN)   :: PSIGQSAT  ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODREF
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PSIGS   ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST  ! Absolute Pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PZZ     ! height of model layer
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXN    ! Exner function
!
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRI_MF! Convective Mass Flux ice mixing ratio
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRV     ! Water vapor m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRC     ! Cloud water m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTH     ! Theta to adjust
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                   ! s'rc'/2Sigma_s2 at time t+1
                                                   ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PCLDFR  ! Cloud fraction
!
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PRIS ! Cloud ice  m.r. at t+1
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PRR  ! Rain water m.r. to adjust
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PRI  ! Cloud ice  m.r. to adjust
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PRS  ! Aggregate  m.r. to adjust
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PRG  ! Graupel    m.r. to adjust
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRH  ! Hail       m.r. to adjust
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RV ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RC ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RI ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_TH ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  PHLC_HRC
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  PHLC_HCF
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  PHLI_HRI
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  PHLI_HCF
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(IN)   :: PICE_CLD_WGT
!
END SUBROUTINE ICE_ADJUST
!
END INTERFACE
!
END MODULE MODI_ICE_ADJUST

