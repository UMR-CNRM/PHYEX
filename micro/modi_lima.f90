MODULE MODI_LIMA
!
IMPLICIT NONE
INTERFACE
!
   SUBROUTINE LIMA ( D, CST, BUCONF, TBUDGETS, KBUDGETS,                     &
                     PTSTEP, OELEC,                                          &
                     PRHODREF, PEXNREF, PDZZ,                                &
                     PRHODJ, PPABST,                                         &
                     NCCN, NIFN, NIMM,                                       &
                     PDTHRAD, PTHT, PRT, PSVT, PW_NU,                        &
                     PTHS, PRS, PSVS,                                        &
                     PINPRC, PINDEP, PINPRR, PINPRI, PINPRS, PINPRG, PINPRH, &
                     PEVAP3D, PCLDFR, PICEFR, PPRCFR, PFPR,                  &
                     PLATHAM_IAGGS, PEFIELDW, PSV_ELEC_T, PSV_ELEC_S         )
!
USE MODD_IO,  ONLY: TFILEDATA
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_BUDGET,   ONLY: TBUDGETDATA, TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step
!
LOGICAL,                  INTENT(IN)    :: OELEC      ! if true, cloud electrification is activated
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ       ! Layer thikness (m)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST     ! absolute pressure at t
!
INTEGER,                  INTENT(IN)    :: NCCN       ! for array size declarations
INTEGER,                  INTENT(IN)    :: NIFN       ! for array size declarations
INTEGER,                  INTENT(IN)    :: NIMM       ! for array size declarations
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDTHRAD    ! dT/dt due to radiation
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT       ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT        ! Mixing ratios at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT       ! Concentrations at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU      ! w for CCN activation
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS       ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS        ! Mixing ratios sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS       ! Concentration sources
!
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRC     ! Cloud instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINDEP     ! Cloud droplets deposition
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRR     ! Rain instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRI     ! Rain instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRS     ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRG     ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRH     ! Rain instant precip
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PEVAP3D    ! Rain evap profile
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCLDFR     ! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR     ! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PPRCFR     ! Cloud fraction
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PFPR       ! Precipitation fluxes in altitude
!
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PLATHAM_IAGGS  ! Factor for IAGGS modification due to Efield
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(IN)    :: PEFIELDW   ! Vertical component of the electric field
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN)    :: PSV_ELEC_T ! Charge density at time t
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(INOUT) :: PSV_ELEC_S ! Charge density sources
!
END SUBROUTINE LIMA
END INTERFACE
END MODULE MODI_LIMA
