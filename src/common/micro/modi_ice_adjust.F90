!     ######################
      MODULE MODI_ICE_ADJUST
!     ######################
!
INTERFACE
!
      SUBROUTINE ICE_ADJUST (KKA, KKU, KKL, KRR, HFRAC_ICE, HCONDENS, HLAMBDA3,&
                             HBUNAME, OSUBG_COND, OSIGMAS, OCND2, LHGT_QS,     &
                             HSUBG_MF_PDF, PTSTEP, PSIGQSAT,                   &
                             PRHODJ, PEXNREF, PRHODREF, PSIGS, PMFCONV,        &
                             PPABST, PZZ,                                      &
                             PEXN, PCF_MF, PRC_MF, PRI_MF,                     &
                             PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR,             &
                             PRV, PRC, PRVS, PRCS, PTH, PTHS, PSRCS, PCLDFR,   &
                             PRR, PRI, PRIS, PRS, PRG, TBUDGETS, KBUDGETS,     &
                             PICE_CLD_WGT,                                     &
                             PRH,                                              &
                             POUT_RV, POUT_RC, POUT_RI, POUT_TH,               &
                             PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF)
!
USE MODD_BUDGET,         ONLY: TBUDGETDATA
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN)    :: KKA      !near ground array index
INTEGER,                  INTENT(IN)    :: KKU      !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL      !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER(LEN=1),         INTENT(IN)    :: HFRAC_ICE
CHARACTER(LEN=80),        INTENT(IN)    :: HCONDENS
CHARACTER(LEN=4),         INTENT(IN)    :: HLAMBDA3 ! formulation for lambda3 coeff
CHARACTER(LEN=4),         INTENT(IN)    :: HBUNAME  ! Name of the budget
LOGICAL,                  INTENT(IN)    :: OSUBG_COND ! Switch for Subgrid
                                                    ! Condensation
LOGICAL,                  INTENT(IN)    :: OSIGMAS  ! Switch for Sigma_s:
                                                    ! use values computed in CONDENSATION
                                                    ! or that from turbulence scheme
LOGICAL,                  INTENT(IN)    :: OCND2    ! logical switch to separate liquid
                                                    ! and ice
                                                    ! more rigid (DEFAULT value : .FALSE.)
LOGICAL,                  INTENT(IN)   :: LHGT_QS   ! logical switch for height dependent VQSIGSAT
CHARACTER(LEN=80),        INTENT(IN)   :: HSUBG_MF_PDF
REAL,                     INTENT(IN)   :: PTSTEP    ! Double Time step
                                                    ! (single if cold start)
REAL, DIMENSION(:,:),     INTENT(IN)   :: PSIGQSAT  ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PRHODREF
!
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PSIGS   ! Sigma_s at time t
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PPABST  ! Absolute Pressure at t
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PZZ     ! height of model layer
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PEXN    ! Exner function
!
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)    :: PCF_MF  ! Convective Mass Flux Cloud fraction
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)    :: PRC_MF  ! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)    :: PRI_MF  ! Convective Mass Flux ice mixing ratio
!
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)    :: PRV     ! Water vapor m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)    :: PRC     ! Cloud water m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)    :: PTH     ! Theta to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                   ! s'rc'/2Sigma_s2 at time t+1
                                                   ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(OUT)   :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(OUT)   :: PICLDFR ! ice cloud fraction
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(OUT)   :: PWCLDFR ! water or mixed-phase cloud fraction
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(OUT)   :: PSSIO   ! Super-saturation with respect to ice in the  
                                                              ! supersaturated fraction
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(OUT)   :: PSSIU   ! Sub-saturation with respect to ice in the  
                                                              ! subsaturated fraction 
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(OUT)   :: PIFR    ! Ratio cloud ice moist part to dry part
!
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(INOUT) :: PRIS    ! Cloud ice  m.r. at t+1
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(IN)    :: PRR     ! Rain water m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(IN)    :: PRI     ! Cloud ice  m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(IN)    :: PRS     ! Aggregate  m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(IN)    :: PRG     ! Graupel    m.r. to adjust
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER,                             INTENT(IN)    :: KBUDGETS
REAL, DIMENSION(:,:),   CONTIGUOUS, OPTIONAL, INTENT(IN)   ::  PICE_CLD_WGT
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(IN)   ::  PRH  ! Hail       m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  POUT_RV ! Adjusted value
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  POUT_RC ! Adjusted value
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  POUT_RI ! Adjusted value
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  POUT_TH ! Adjusted value
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  PHLC_HRC
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  PHLC_HCF
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  PHLI_HRI
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  PHLI_HCF
!
END SUBROUTINE ICE_ADJUST
!
END INTERFACE
!
END MODULE MODI_ICE_ADJUST

