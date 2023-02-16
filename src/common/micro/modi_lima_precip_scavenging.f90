!#################################
MODULE MODI_LIMA_PRECIP_SCAVENGING
!#################################
!
  INTERFACE
!
     SUBROUTINE LIMA_PRECIP_SCAVENGING (D, CST, BUCONF, TBUDGETS, KBUDGETS, &
                                        HCLOUD, KLUOUT, KTCOUNT, PTSTEP,    &
                                        PRRT, PRHODREF, PRHODJ, PZZ,        &
                                        PPABST, PTHT, PSVT, PRSVS, PINPAP )
       USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
       use modd_budget,          only: TBUDGETDATA,TBUDGETCONF_t
       USE MODD_CST,             ONLY: CST_t
!
       TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
       TYPE(CST_t),              INTENT(IN)    :: CST
       TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
       TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
       INTEGER,                  INTENT(IN)    :: KBUDGETS
!
       CHARACTER(LEN=4),       INTENT(IN)      :: HCLOUD   ! cloud paramerization
       INTEGER,                INTENT(IN)      :: KLUOUT   ! unit for output listing
       INTEGER,                INTENT(IN)      :: KTCOUNT  ! iteration count
       REAL,                   INTENT(IN)      :: PTSTEP   ! Double timestep except 
                                                           ! for the first time step
!
       REAL, DIMENSION(:,:,:), INTENT(IN)      :: PRRT     ! Rain mixing ratio at t
       REAL, DIMENSION(:,:,:), INTENT(IN)      :: PRHODREF ! Air Density [kg/m**3]
       REAL, DIMENSION(:,:,:), INTENT(IN)      :: PRHODJ   ! Dry Density [kg]
       REAL, DIMENSION(:,:,:), INTENT(IN)      :: PZZ      ! Altitude
       REAL, DIMENSION(:,:,:), INTENT(IN)      :: PPABST   ! Absolute pressure at t
       REAL, DIMENSION(:,:,:), INTENT(IN)      :: PTHT     ! Theta at time t 
!
       REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT   ! Particle Concentration [/m**3]
       REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS  ! Total Number Scavenging Rate
!
       REAL, DIMENSION(:,:),   INTENT(INOUT)   :: PINPAP
     END SUBROUTINE LIMA_PRECIP_SCAVENGING
  END INTERFACE
END MODULE MODI_LIMA_PRECIP_SCAVENGING
