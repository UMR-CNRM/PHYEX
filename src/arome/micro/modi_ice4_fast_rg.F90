MODULE MODI_ICE4_FAST_RG
INTERFACE
SUBROUTINE ICE4_FAST_RG(KSIZE, LDSOFT, PCOMPUTE, KRR, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, PCIT, &
                       &PLBDAR, PLBDAS, PLBDAG, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                       &PRGSI, PRGSI_MR, &
                       &PWETG, &
                       &PRICFRRG, PRRCFRIG, PRICFRR, PRCWETG, PRIWETG, PRRWETG, PRSWETG, &
                       &PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG, PRWETGH, PRWETGH_MR, PRGMLTR, &
                       &PRG_TEND, &
                       &PA_TH, PA_RC, PA_RR, PA_RI, PA_RS, PA_RG, PA_RH, PB_RG, PB_RH)
IMPLICIT NONE
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
INTEGER,                      INTENT(IN)    :: KRR      ! Number of moist variable
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGSI    ! Graupel tendency by other processes
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGSI_MR ! Graupel mr change by other processes
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PWETG    ! 1. where graupel grows in wet mode, 0. elsewhere
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRICFRRG ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRCFRIG ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRICFRR  ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRWETGH  ! Conversion of graupel into hail
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRWETGH_MR ! Conversion of graupel into hail, mr change
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRGMLTR  ! Melting of the graupel
REAL, DIMENSION(KSIZE, 8),    INTENT(INOUT) :: PRG_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RH
END SUBROUTINE ICE4_FAST_RG
END INTERFACE
END MODULE MODI_ICE4_FAST_RG
