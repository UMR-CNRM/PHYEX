MODULE MODI_LIMA
!
IMPLICIT NONE
INTERFACE
!
SUBROUTINE LIMA ( LIMAP, LIMAW, LIMAC, LIMAM, TNSV, D, CST, NEBN,         &
                  ICED, ICEP, ELECD, ELECP,                               &
                  BUCONF, TBUDGETS, HACTCCN, KBUDGETS, KRR,               &
                  PTSTEP, OELEC,                                          &
                  PRHODREF, PEXNREF, PDZZ, PTHVREFZIKB,                   &
                  PRHODJ, PPABST,                                         &
                  KCARB, KSOA, KSP, ODUST, OSALT, OORILAM,                &
                  ODTHRAD, PDTHRAD, PTHT, PRT, PSVT, PW_NU,               &
                  PAERO,PSOLORG, PMI, PTHS, PRS, PSVS,                    &
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
USE MODD_BUDGET,   ONLY: TBUDGETDATA_PTR, TBUDGETCONF_T
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
CHARACTER(LEN=4),         INTENT(IN)    :: HACTCCN  ! kind of CCN activation
TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER,                  INTENT(IN)    :: KBUDGETS
INTEGER,                  INTENT(IN)    :: KRR
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step
!
LOGICAL,                  INTENT(IN)    :: OELEC      ! if true, cloud electrification is activated
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PDZZ       ! Layer thikness (m)
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PRHODJ     ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PPABST     ! absolute pressure at t
INTEGER,                  INTENT(IN)    :: KCARB, KSOA, KSP ! for array size declarations
LOGICAL,                  INTENT(IN)    :: ODUST, OSALT, OORILAM
!
LOGICAL,                                 INTENT(IN)   :: ODTHRAD    ! Use radiative temperature tendency
REAL, DIMENSION(MERGE(D%NIJT,0,ODTHRAD), &
                MERGE(D%NKT,0,ODTHRAD)),   INTENT(IN) :: PDTHRAD   ! dT/dt due to radiation
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PTHT       ! Theta at time t
REAL, DIMENSION(D%NIJT, D%NKT, KRR), INTENT(IN) :: PRT        ! Mixing ratios at time t
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), INTENT(IN) :: PSVT       ! Concentrations at time t
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PW_NU      ! w for CCN activation
REAL, DIMENSION(D%NIJT, D%NKT ,TNSV%NSV), INTENT(INOUT) :: PAERO    ! Aerosol concentration
REAL, DIMENSION(D%NIJT, D%NKT, 10),  INTENT(IN)    :: PSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(D%NIJT, D%NKT, KSP+KCARB+KSOA), INTENT(IN)    :: PMI
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(INOUT)    :: PTHS       ! Theta source
REAL, DIMENSION(D%NIJT, D%NKT, KRR), INTENT(INOUT) :: PRS        ! Mixing ratios sources
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), INTENT(INOUT) :: PSVS       ! Concentration sources
!
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRC     ! Cloud instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINDEP     ! Cloud droplets deposition
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRR     ! Rain instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRI     ! Rain instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRS     ! Snow instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRG     ! Graupel instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRH     ! Rain instant precip
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(OUT)   :: PEVAP3D    ! Rain evap profile
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(INOUT) :: PCLDFR     ! Cloud fraction
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(INOUT) :: PICEFR     ! Cloud fraction
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(INOUT) :: PPRCFR     ! Cloud fraction
REAL, DIMENSION(D%NIJT, D%NKT, KRR), INTENT(OUT) :: PFPR    ! Precipitation fluxes in altitude
!
REAL, DIMENSION(D%NIJT, D%NKT),   OPTIONAL, INTENT(IN)       :: PLATHAM_IAGGS  ! Factor for IAGGS modification due to Efield
REAL, DIMENSION(D%NIJT, D%NKT),   OPTIONAL, INTENT(IN)       :: PEFIELDW   ! Vertical component of the electric field
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), OPTIONAL, INTENT(IN)    :: PSV_ELEC_T ! Charge density at time t
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), OPTIONAL, INTENT(INOUT) :: PSV_ELEC_S ! Charge density sources
!
REAL, INTENT(IN)                :: PTHVREFZIKB ! Reference thv at IKB for electricity
END SUBROUTINE LIMA
END INTERFACE
END MODULE MODI_LIMA
