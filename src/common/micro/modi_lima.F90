MODULE MODI_LIMA
!
IMPLICIT NONE
INTERFACE
!
   SUBROUTINE LIMA ( LIMAP, LIMAW, LIMAC, LIMAM,                             &
                  D, CST, ICED, ICEP, ELECD, ELECP, BUCONF, TBUDGETS, KBUDGETS, KRR, &
                  PTSTEP, OELEC,                                          &
                  PRHODREF, PEXNREF, PDZZ, PTHVREFZIKB,                   &
                  PRHODJ, PPABST,                                         &
                  KCCN, KIFN, KIMM,                                       &
                  ODTHRAD, PDTHRAD, PTHT, PRT, PSVT, PW_NU,               &
                  PTHS, PRS, PSVS,                                        &
                  PINPRC, PINDEP, PINPRR, PINPRI, PINPRS, PINPRG, PINPRH, &
                  PEVAP3D, PCLDFR, PICEFR, PPRCFR, PFPR,                  &
                  PLATHAM_IAGGS, PEFIELDW, PSV_ELEC_T, PSV_ELEC_S         )
!
USE MODD_IO,  ONLY: TFILEDATA
USE MODD_PARAM_LIMA_MIXED, ONLY:PARAM_LIMA_MIXED_t
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_t
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_t
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_t
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_RAIN_ICE_DESCR_n,ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n,ONLY: RAIN_ICE_PARAM_t
USE MODD_ELEC_PARAM,      ONLY: ELEC_PARAM_t
USE MODD_ELEC_DESCR,      ONLY: ELEC_DESCR_t
USE MODD_BUDGET,   ONLY: TBUDGETDATA, TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
USE MODD_NSV, ONLY: NSV
IMPLICIT NONE
!
TYPE(PARAM_LIMA_MIXED_t),INTENT(IN)::LIMAM
TYPE(PARAM_LIMA_COLD_t),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_WARM_t),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_t),INTENT(IN)::LIMAP
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(ELEC_PARAM_t),       INTENT(IN)    :: ELECP   ! electrical parameters
TYPE(ELEC_DESCR_t),       INTENT(IN)    :: ELECD   ! electrical descriptive csts
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER,                  INTENT(IN)    :: KBUDGETS
INTEGER,                  INTENT(IN)    :: KRR
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step
!
LOGICAL,                  INTENT(IN)    :: OELEC      ! if true, cloud electrification is activated
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)    :: PDZZ       ! Layer thikness (m)
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)    :: PRHODJ     ! Dry density * Jacobian
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)    :: PPABST     ! absolute pressure at t
!
INTEGER,                  INTENT(IN)    :: KCCN       ! for array size declarations
INTEGER,                  INTENT(IN)    :: KIFN       ! for array size declarations
INTEGER,                  INTENT(IN)    :: KIMM       ! for array size declarations
!
LOGICAL,                                 INTENT(IN)   :: ODTHRAD    ! Use radiative temperature tendency
REAL, DIMENSION(MERGE(D%NIT,0,ODTHRAD), &
                MERGE(D%NJT,0,ODTHRAD), &
                MERGE(D%NKT,0,ODTHRAD)),   INTENT(IN) :: PDTHRAD   ! dT/dt due to radiation
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)    :: PTHT       ! Theta at time t
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, KRR), INTENT(IN) :: PRT        ! Mixing ratios at time t
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, NSV), INTENT(IN) :: PSVT       ! Concentrations at time t
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)    :: PW_NU      ! w for CCN activation
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(INOUT)    :: PTHS       ! Theta source
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, KRR), INTENT(INOUT) :: PRS        ! Mixing ratios sources
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, NSV), INTENT(INOUT) :: PSVS       ! Concentration sources
!
REAL, DIMENSION(D%NIT, D%NJT),     INTENT(OUT)        :: PINPRC     ! Cloud instant precip
REAL, DIMENSION(D%NIT, D%NJT),     INTENT(OUT)        :: PINDEP     ! Cloud droplets deposition
REAL, DIMENSION(D%NIT, D%NJT),     INTENT(OUT)        :: PINPRR     ! Rain instant precip
REAL, DIMENSION(D%NIT, D%NJT),     INTENT(OUT)        :: PINPRI     ! Rain instant precip
REAL, DIMENSION(D%NIT, D%NJT),     INTENT(OUT)        :: PINPRS     ! Snow instant precip
REAL, DIMENSION(D%NIT, D%NJT),     INTENT(OUT)        :: PINPRG     ! Graupel instant precip
REAL, DIMENSION(D%NIT, D%NJT),     INTENT(OUT)        :: PINPRH     ! Rain instant precip
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(OUT)   :: PEVAP3D    ! Rain evap profile
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(INOUT) :: PCLDFR     ! Cloud fraction
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(INOUT) :: PICEFR     ! Cloud fraction
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(INOUT) :: PPRCFR     ! Cloud fraction
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, KRR), INTENT(INOUT) :: PFPR    ! Precipitation fluxes in altitude
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   OPTIONAL, INTENT(IN)       :: PLATHAM_IAGGS  ! Factor for IAGGS modification due to Efield
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   OPTIONAL, INTENT(IN)       :: PEFIELDW   ! Vertical component of the electric field
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, NSV), OPTIONAL, INTENT(IN)    :: PSV_ELEC_T ! Charge density at time t
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, NSV), OPTIONAL, INTENT(INOUT) :: PSV_ELEC_S ! Charge density sources
!
REAL, INTENT(IN)                :: PTHVREFZIKB ! Reference thv at IKB for electricity
END SUBROUTINE LIMA
END INTERFACE
END MODULE MODI_LIMA
