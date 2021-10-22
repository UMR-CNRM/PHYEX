!     ######################
      MODULE MODI_ICE_ADJUST
!     ######################
!
INTERFACE
!
      SUBROUTINE ICE_ADJUST (KKA, KKU, KKL, KRR, HFRAC_ICE,                    &
                             HBUNAME, OSUBG_COND, OSIGMAS, OCND2,              &
                             PTSTEP, PSIGQSAT,                                 &
                             PRHODJ, PEXNREF, PSIGS, PMFCONV, PPABST, PZZ,     &
                             PEXN, PCF_MF, PRC_MF, PRI_MF,                     &
                             PRV, PRC, PRVS, PRCS, PTH, PTHS, PSRCS, PCLDFR ,  &
                             PRR, PRI, PRIS, PRS, PRG,                         &
                             PRH, POUT_RV, POUT_RC, POUT_RI, POUT_TH,          &
                             YSPP_PSIGQSAT, YSPP_ICE_CLD_WGT,                  &
                             YDDDH, YDLDDH, YDMDDH                 )
!
USE MODD_SPP_TYPE
!
USE DDH_MIX, ONLY : TYP_DDH
USE YOMLDDH, ONLY : TLDDH
USE YOMMDDH, ONLY : TMDDH
!
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER*1,              INTENT(IN)    :: HFRAC_ICE
CHARACTER*4,              INTENT(IN)    :: HBUNAME  ! Name of the budget
LOGICAL,                  INTENT(IN)    :: OSUBG_COND ! Switch for Subgrid 
                                                    ! Condensation
LOGICAL                                 :: OSIGMAS  ! Switch for Sigma_s: 
                                                    ! use values computed in CONDENSATION
                                                    ! or that from turbulence scheme
LOGICAL,                  INTENT(IN)    :: OCND2    ! logical switch to separate liquid
                                                    ! and ice
REAL,                     INTENT(IN)   :: PTSTEP    ! Double Time step
                                                    ! (single if cold start)
REAL,                     INTENT(IN)   :: PSIGQSAT  ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF ! Reference Exner function
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PSIGS   ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST  ! Absolute Pressure at t        
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PZZ     ! height of model layer
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXN    ! Exner function
!
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRI_MF! Convective Mass Flux ice mixing ratio
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
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
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRIS   ! Cloud ice  m.r. at t+1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRR    ! Rain water m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRI    ! Cloud ice  m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRS    ! Aggregate  m.r. to adjust
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRG    ! Graupel    m.r. to adjust
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)   ::  PRH  ! Hail       m.r. to adjust
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RV ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RC ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_RI ! Adjusted value
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)  ::  POUT_TH ! Adjusted value
TYPE(TSPP_CONFIG_MPA),    INTENT(IN)   :: YSPP_PSIGQSAT,YSPP_ICE_CLD_WGT
TYPE(TYP_DDH)                   , INTENT(INOUT) :: YDDDH
TYPE(TLDDH)                   , INTENT(IN) :: YDLDDH
TYPE(TMDDH)                   , INTENT(IN) :: YDMDDH
!
!
END SUBROUTINE ICE_ADJUST
!
END INTERFACE
!
END MODULE MODI_ICE_ADJUST

