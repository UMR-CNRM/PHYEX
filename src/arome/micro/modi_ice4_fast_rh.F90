MODULE MODI_ICE4_FAST_RH
INTERFACE
SUBROUTINE ICE4_FAST_RH(KSIZE, LDSOFT, PCOMPUTE, PWETG, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, &
                       &PLBDAS, PLBDAG, PLBDAR, PLBDAH, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT, &
                       &PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH, &
                       &PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH, PRDRYHG, PRHMLTR, &
                       &PRH_TEND, &
                       &PA_TH, PA_RC, PA_RR, PA_RI, PA_RS, PA_RG, PA_RH)
IMPLICIT NONE
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PWETG    ! 1. where graupel grows in wet mode, 0. elsewhere
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the rain      distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAH   ! Slope parameter of the hail      distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHT     ! Hail m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRGWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRGDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRDRYHG  ! Conversion of hailstone into graupel
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRHMLTR  ! Melting of the hailstones
REAL, DIMENSION(KSIZE, 10),   INTENT(INOUT) :: PRH_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RH
END SUBROUTINE ICE4_FAST_RH
END INTERFACE
END MODULE MODI_ICE4_FAST_RH
