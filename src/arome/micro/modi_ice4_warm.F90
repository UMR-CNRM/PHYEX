MODULE MODI_ICE4_WARM
INTERFACE
SUBROUTINE ICE4_WARM(KSIZE, LDSOFT, PCOMPUTE, HSUBG_RC_RR_ACCR, HSUBG_RR_EVAP, &
                    &PRHODREF, PLVFACT, PT, PPRES, PTHT, &
                    &PLBDAR, PLBDAR_RF, PKA, PDV, PCJ, &
                    &PHLC_LCF, PHLC_HCF, PHLC_LRC, PHLC_HRC, &
                    &PCF, PRF, &
                    &PRVT, PRCT, PRRT, &
                    &PRCAUTR, PRCACCR, PRREVAV, &
                    &PA_TH, PA_RV, PA_RC, PA_RR)
IMPLICIT NONE
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
CHARACTER*80,                 INTENT(IN)    :: HSUBG_RC_RR_ACCR ! subgrid rc-rr accretion
CHARACTER*80,                 INTENT(IN)    :: HSUBG_RR_EVAP ! subgrid rr evaporation
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PTHT     ! Theta at time t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR_RF!like PLBDAR but for the Rain Fraction part
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_HCF ! HLCLOUDS : fraction of High Cloud Fraction in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_LCF ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_HRC ! HLCLOUDS : LWC that is High LWC in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_LRC ! HLCLOUDS : LWC that is Low  LWC in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCF      ! Cloud fraction
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRF      ! Rain fraction
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCAUTR   ! Autoconversion of r_c for r_r production
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCACCR  ! Accretion of r_c for r_r production
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRREVAV  ! Evaporation of r_r
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RV
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
END SUBROUTINE ICE4_WARM
END INTERFACE
END MODULE MODI_ICE4_WARM
