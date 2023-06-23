!     #############################
      MODULE MODI_LIMA_ADJUST_SPLIT
!     #############################
!
IMPLICIT NONE
INTERFACE
!
   SUBROUTINE LIMA_ADJUST_SPLIT(D, CST, BUCONF, TBUDGETS, KBUDGETS, &
                             KRR, KMI, HCONDENS, HLAMBDA3,        &
                             OSUBG_COND, OSIGMAS, PTSTEP, PSIGQSAT,             &
                             PRHODREF, PRHODJ, PEXNREF, PSIGS, PMFCONV,         &
                             PPABST, PPABSTT, PZZ, PDTHRAD, PW_NU,              &
                             PRT, PRS, PSVT, PSVS,                              &
                             PTHS, PSRCS, PCLDFR, PICEFR, PRC_MF, PRI_MF, PCF_MF)
!
!USE MODD_IO,    ONLY: TFILEDATA
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_BUDGET,   ONLY: TBUDGETDATA, TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),         INTENT(IN)   :: D
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
INTEGER,                  INTENT(IN)   :: KRR        ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KMI        ! Model index 
CHARACTER(len=80),        INTENT(IN)   :: HCONDENS
CHARACTER(len=4),         INTENT(IN)   :: HLAMBDA3   ! formulation for lambda3 coeff
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid
                                                     ! Condensation
LOGICAL,                  INTENT(IN)   :: OSIGMAS    ! Switch for Sigma_s:
                                                     ! use values computed in CONDENSATION
                                                     ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PTSTEP     ! Time step
REAL,                     INTENT(IN)   :: PSIGQSAT   ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODREF  ! Dry density of the 
                                                     ! reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ    ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF   ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PSIGS     ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PMFCONV   ! 
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST    ! Absolute Pressure at t     
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABSTT   ! Absolute Pressure at t+dt     
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PZZ       !     
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PDTHRAD   ! Radiative temperature tendency
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PW_NU     ! updraft velocity used for
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT       ! m.r. at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS       ! m.r. source
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT ! Concentrations at time t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS ! Concentration sources
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS      ! Theta source
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS     ! Second-order flux
                                                     ! s'rc'/2Sigma_s2 at time t+1
                                                     ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   :: PCLDFR    ! Cloud fraction          
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   :: PICEFR    ! Cloud fraction          
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRI_MF! Convective Mass Flux ice mixing ratio
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
!
END SUBROUTINE LIMA_ADJUST_SPLIT
END INTERFACE
END MODULE MODI_LIMA_ADJUST_SPLIT
