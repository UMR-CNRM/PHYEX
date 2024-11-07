MODULE MODI_LIMA
!
IMPLICIT NONE
INTERFACE
!
   SUBROUTINE LIMA ( LIMAP, LIMAW, LIMAC, LIMAM, TNSV,                    &
                  D, CST, NEBN, ICED, ICEP, ELECD, ELECP, BUCONF, TBUDGETS, KBUDGETS, KRR, &
                  PTSTEP, OELEC,                                          &
                  PRHODREF, PEXNREF, PDZZ, PTHVREFZIKB,                   &
                  PRHODJ, PPABST,                                         &
                  ODTHRAD, PDTHRAD, PTHT, PRT, PSVT, PW_NU,               &
                  PTHS, PRS, PSVS,                                        &
                  PINPRC, PINDEP, PINPRR, PINPRI, PINPRS, PINPRG, PINPRH, &
                  PEVAP3D, PCLDFR, PICEFR, PPRCFR, PFPR,                  &
                  PLATHAM_IAGGS, PEFIELDW, PSV_ELEC_T, PSV_ELEC_S         )
!
USE MODD_IO,  ONLY: TFILEDATA
USE MODD_PARAM_LIMA_MIXED, ONLY:PARAM_LIMA_MIXED_T
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_RAIN_ICE_DESCR_N,ONLY: RAIN_ICE_DESCR_T
USE MODD_RAIN_ICE_PARAM_N,ONLY: RAIN_ICE_PARAM_T
USE MODD_ELEC_PARAM,      ONLY: ELEC_PARAM_T
USE MODD_ELEC_DESCR,      ONLY: ELEC_DESCR_T
USE MODD_BUDGET,   ONLY: TBUDGETDATA, TBUDGETCONF_T
USE MODD_CST,            ONLY: CST_T
USE MODD_NEB_N,          ONLY: NEB_T
USE MODD_NSV, ONLY: NSV_T
IMPLICIT NONE
!
TYPE(PARAM_LIMA_MIXED_T),INTENT(IN)::LIMAM
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(NSV_T),              INTENT(IN)    :: TNSV
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
TYPE(CST_T),              INTENT(IN)    :: CST
TYPE(NEB_T),              INTENT(IN)    :: NEBN
TYPE(RAIN_ICE_DESCR_T),   INTENT(IN)    :: ICED
TYPE(RAIN_ICE_PARAM_T),   INTENT(IN)    :: ICEP
TYPE(ELEC_PARAM_T),       INTENT(IN)    :: ELECP   ! electrical parameters
TYPE(ELEC_DESCR_T),       INTENT(IN)    :: ELECD   ! electrical descriptive csts
TYPE(TBUDGETCONF_T),      INTENT(IN)    :: BUCONF
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
LOGICAL,                                 INTENT(IN)   :: ODTHRAD    ! Use radiative temperature tendency
REAL, DIMENSION(MERGE(D%NIT,0,ODTHRAD), &
                MERGE(D%NJT,0,ODTHRAD), &
                MERGE(D%NKT,0,ODTHRAD)),   INTENT(IN) :: PDTHRAD   ! dT/dt due to radiation
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)    :: PTHT       ! Theta at time t
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, KRR), INTENT(IN) :: PRT        ! Mixing ratios at time t
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, TNSV%NSV), INTENT(IN) :: PSVT       ! Concentrations at time t
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)    :: PW_NU      ! w for CCN activation
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(INOUT)    :: PTHS       ! Theta source
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, KRR), INTENT(INOUT) :: PRS        ! Mixing ratios sources
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, TNSV%NSV), INTENT(INOUT) :: PSVS       ! Concentration sources
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
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, TNSV%NSV), OPTIONAL, INTENT(IN)    :: PSV_ELEC_T ! Charge density at time t
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, TNSV%NSV), OPTIONAL, INTENT(INOUT) :: PSV_ELEC_S ! Charge density sources
!
REAL, INTENT(IN)                :: PTHVREFZIKB ! Reference thv at IKB for electricity
END SUBROUTINE LIMA
END INTERFACE
END MODULE MODI_LIMA
