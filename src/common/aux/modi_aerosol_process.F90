!     ######spl
       MODULE MODI_AEROSOL_PROCESS
!      ####################
!
INTERFACE
 SUBROUTINE AEROSOL_PROCESS(PHYEX, KKA,KKU,KKL,KNAERO,PTSTEP, &
        &  PPABST,PTHT,PDZZ,PRHODREF,PRCT,PRIT, &
        &  PCLDFR,PFPR, &
        &  PLSM,PUGST,PVGST, &
        &  PAERMR,PAERTEND)
!
USE MODD_PHYEX, ONLY: PHYEX_t
!
TYPE(PHYEX_t),            INTENT(IN) :: PHYEX
INTEGER,                  INTENT(IN)    :: KKA
INTEGER,                  INTENT(IN)    :: KKU
INTEGER,                  INTENT(IN)    :: KKL
INTEGER,                  INTENT(IN)    :: KNAERO
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold
                                                   ! start)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCLDFR  ! Convective Mass Flux Cloud fraction
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PFPR    ! upper-air precipitation fluxes
REAL, DIMENSION(:),       INTENT(IN)    :: PLSM    ! Land/Sea Mask
REAL, DIMENSION(:),       INTENT(IN)    :: PUGST   ! 10m u-wind gust
REAL, DIMENSION(:),       INTENT(IN)    :: PVGST   ! 10m v-wind gust
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PAERMR  ! Aerosol mixing ratio
REAL, DIMENSION(:,:,:,:), INTENT(  OUT) :: PAERTEND! Aerosol tendencies: scavenging + emission
END SUBROUTINE AEROSOL_PROCESS
END INTERFACE
END MODULE MODI_AEROSOL_PROCESS

