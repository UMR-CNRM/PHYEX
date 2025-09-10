!     #############################
      MODULE MODI_LIMA_ADJUST_SPLIT
!     #############################
!
IMPLICIT NONE
INTERFACE
!
   SUBROUTINE LIMA_ADJUST_SPLIT(LIMAP, LIMAW, TNSV, D, CST, NEBN, TURBN, BUCONF, TBUDGETS, KBUDGETS, &
                             KRR, HCONDENS, HLAMBDA3,                           &
                             KCARB, KSOA, KSP, ODUST, OSALT, OORILAM,           &
                             OSUBG_COND, OSIGMAS, PTSTEP, PSIGQSAT,             &
                             PRHODREF, PRHODJ, PEXNREF, PSIGS, OMFCONV, PMFCONV,&
                             PPABST, PZZ, ODTHRAD, PDTHRAD, PW_NU,              &
                             PRT, PRS, PSVT, PSVS,                              &
                             HACTCCN, PAERO,PSOLORG, PMI,                       &
                             PTHS, OCOMPUTE_SRC, PSRCS, PCLDFR, PICEFR,         &
                             PRC_MF, PRI_MF, PCF_MF)
!
!USE MODD_IO,    ONLY: TFILEDATA
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_T
USE MODD_BUDGET,         ONLY: TBUDGETDATA_PTR, TBUDGETCONF_T
USE MODD_CST,            ONLY: CST_T
USE MODD_NSV, ONLY: NSV_T
USE MODD_NEB_N,            ONLY: NEB_T
USE MODD_TURB_N,            ONLY: TURB_T
USE MODD_PARAM_LIMA, ONLY: PARAM_LIMA_T
USE MODD_PARAM_LIMA_WARM, ONLY: PARAM_LIMA_WARM_T
IMPLICIT NONE
!
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(NSV_T),              INTENT(IN)    :: TNSV
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
TYPE(CST_T),              INTENT(IN)    :: CST
TYPE(NEB_T),              INTENT(IN)    :: NEBN
TYPE(TURB_T),             INTENT(IN)    :: TURBN
TYPE(TBUDGETCONF_T),      INTENT(IN)    :: BUCONF
TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
INTEGER,                  INTENT(IN)   :: KRR        ! Number of moist variables
CHARACTER(LEN=80),        INTENT(IN)   :: HCONDENS
CHARACTER(LEN=4),         INTENT(IN)   :: HLAMBDA3   ! formulation for lambda3 coeff
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid
                                                     ! Condensation
LOGICAL,                  INTENT(IN)   :: OSIGMAS    ! Switch for Sigma_s:
                                                     ! use values computed in CONDENSATION
                                                     ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PTSTEP     ! Time step
REAL, DIMENSION(D%NIJT),  INTENT(IN)   :: PSIGQSAT   ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)   ::  PRHODREF  ! Dry density of the 
                                                             ! reference state
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)   ::  PRHODJ    ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)   ::  PEXNREF   ! Reference Exner function
REAL, DIMENSION(MERGE(D%NIJT,0,NEBN%LSUBG_COND), &
                MERGE(D%NKT,0,NEBN%LSUBG_COND)),   INTENT(IN)   ::  PSIGS     ! Sigma_s at time t
LOGICAL,                                  INTENT(IN)    ::  OMFCONV ! T to use PMFCONV
REAL, DIMENSION(MERGE(D%NIJT,0,OMFCONV), &
                MERGE(D%NKT,0,OMFCONV)),   INTENT(IN)   ::  PMFCONV   ! 
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)   ::  PPABST    ! Absolute Pressure at t     
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)   ::  PZZ       !     
LOGICAL,                                INTENT(IN)   :: ODTHRAD    ! Use radiative temperature tendency
REAL, DIMENSION(MERGE(D%NIJT,0,ODTHRAD), &
                MERGE(D%NKT,0,ODTHRAD)),   INTENT(IN) :: PDTHRAD   ! Radiative temperature tendency
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PW_NU     ! updraft velocity used for
INTEGER,                  INTENT(IN)    :: KCARB, KSOA, KSP ! for array size declarations
LOGICAL,                  INTENT(IN)    :: ODUST, OSALT, OORILAM
REAL, DIMENSION(D%NIJT, D%NKT ,TNSV%NSV), INTENT(INOUT) :: PAERO    ! Aerosol concentration
REAL, DIMENSION(D%NIJT, D%NKT, 10),  INTENT(IN)    :: PSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(D%NIJT, D%NKT, KSP+KCARB+KSOA), INTENT(IN)    :: PMI
CHARACTER(LEN=4),         INTENT(IN)    :: HACTCCN  ! kind of CCN activation
!
REAL, DIMENSION(D%NIJT, D%NKT, KRR), INTENT(IN)    :: PRT       ! m.r. at t
!
REAL, DIMENSION(D%NIJT, D%NKT, KRR), INTENT(INOUT) :: PRS       ! m.r. source
!
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), INTENT(IN)    :: PSVT ! Concentrations at time t
!
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), INTENT(INOUT) :: PSVS ! Concentration sources
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(INOUT) :: PTHS      ! Theta source
!
LOGICAL,                                      INTENT(IN)    :: OCOMPUTE_SRC ! T to comput PSRCS
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC), &
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(OUT)   :: PSRCS     ! Second-order flux
                                                                         ! s'rc'/2Sigma_s2 at time t+1
                                                                         ! multiplied by Lambda_3
REAL, DIMENSION(D%NIJT, D%NKT),     INTENT(OUT) :: PCLDFR    ! Cloud fraction          
REAL, DIMENSION(D%NIJT, D%NKT),     INTENT(OUT) :: PICEFR    ! Cloud fraction          
REAL, DIMENSION(D%NIJT, D%NKT),     INTENT(IN)  :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(D%NIJT, D%NKT),     INTENT(IN)  :: PRI_MF! Convective Mass Flux ice mixing ratio
REAL, DIMENSION(D%NIJT, D%NKT),     INTENT(IN)  :: PCF_MF! Convective Mass Flux Cloud fraction 
!
END SUBROUTINE LIMA_ADJUST_SPLIT
END INTERFACE
END MODULE MODI_LIMA_ADJUST_SPLIT
