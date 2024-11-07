!     #############################
      MODULE MODI_LIMA_ADJUST_SPLIT
!     #############################
!
IMPLICIT NONE
INTERFACE
!
   SUBROUTINE LIMA_ADJUST_SPLIT(LIMAP, LIMAW, TNSV, D, CST, NEBN, TURBN, BUCONF, TBUDGETS, KBUDGETS, &
                             KRR, KMI, HCONDENS, HLAMBDA3,                      &
                             OSUBG_COND, OSIGMAS, PTSTEP, PSIGQSAT,             &
                             PRHODREF, PRHODJ, PEXNREF, PSIGS, OMFCONV, PMFCONV,&
                             PPABST, PPABSTT, PZZ, ODTHRAD, PDTHRAD, PW_NU,     &
                             PRT, PRS, PSVT, PSVS,                              &
                             PTHS, OCOMPUTE_SRC, PSRCS, PCLDFR, PICEFR,         &
                             PRC_MF, PRI_MF, PCF_MF)
!
!USE MODD_IO,    ONLY: TFILEDATA
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_BUDGET,   ONLY: TBUDGETDATA, TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
USE MODD_NSV, ONLY: NSV_t
USE MODD_NEB_n,            ONLY: NEB_t
USE MODD_TURB_n,            ONLY: TURB_t
USE MODD_PARAM_LIMA, ONLY: PARAM_LIMA_t
USE MODD_PARAM_LIMA_WARM, ONLY: PARAM_LIMA_WARM_t
IMPLICIT NONE
!
TYPE(PARAM_LIMA_t),INTENT(IN)::LIMAP
TYPE(PARAM_LIMA_WARM_t),INTENT(IN)::LIMAW
TYPE(NSV_t),              INTENT(IN)    :: TNSV
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(NEB_t),              INTENT(IN)    :: NEBN
TYPE(TURB_t),              INTENT(IN)   :: TURBN
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
INTEGER,                  INTENT(IN)   :: KRR        ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KMI        ! Model index 
CHARACTER(len=80),        INTENT(IN)   :: HCONDENS
CHARACTER(len=4),         INTENT(IN)   :: HLAMBDA3   ! formulation for lambda3 coeff
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid
                                                     ! Condensation
LOGICAL,                  INTENT(IN)   :: OSIGMAS    ! Switch for Sigma_s:
                                                     ! use values computed in CONDENSATION
                                                     ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PTSTEP     ! Time step
REAL,                     INTENT(IN)   :: PSIGQSAT   ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PRHODREF  ! Dry density of the 
                                                                   ! reference state
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PRHODJ    ! Dry density * Jacobian
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PEXNREF   ! Reference Exner function
REAL, DIMENSION(MERGE(D%NIT,0,NEBN%LSUBG_COND), &
                MERGE(D%NJT,0,NEBN%LSUBG_COND), &
                MERGE(D%NKT,0,NEBN%LSUBG_COND)),   INTENT(IN)   ::  PSIGS     ! Sigma_s at time t
LOGICAL,                                  INTENT(IN)    ::  OMFCONV ! T to use PMFCONV
REAL, DIMENSION(MERGE(D%NIT,0,OMFCONV), &
                MERGE(D%NJT,0,OMFCONV), &
                MERGE(D%NKT,0,OMFCONV)),   INTENT(IN)   ::  PMFCONV   ! 
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PPABST    ! Absolute Pressure at t     
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PPABSTT   ! Absolute Pressure at t+dt     
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PZZ       !     
LOGICAL,                                INTENT(IN)   :: ODTHRAD    ! Use radiative temperature tendency
REAL, DIMENSION(MERGE(D%NIT,0,ODTHRAD), &
                MERGE(D%NJT,0,ODTHRAD), &
                MERGE(D%NKT,0,ODTHRAD)),   INTENT(IN) :: PDTHRAD   ! Radiative temperature tendency
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)    :: PW_NU     ! updraft velocity used for
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, KRR), INTENT(IN)    :: PRT       ! m.r. at t
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, KRR), INTENT(INOUT) :: PRS       ! m.r. source
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, TNSV%NSV), INTENT(IN)    :: PSVT ! Concentrations at time t
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, TNSV%NSV), INTENT(INOUT) :: PSVS ! Concentration sources
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(INOUT) :: PTHS      ! Theta source
!
LOGICAL,                                      INTENT(IN)    :: OCOMPUTE_SRC ! T to comput PSRCS
REAL, DIMENSION(MERGE(D%NIT,0,OCOMPUTE_SRC), &
                MERGE(D%NJT,0,OCOMPUTE_SRC), &
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(OUT)   :: PSRCS     ! Second-order flux
                                                                         ! s'rc'/2Sigma_s2 at time t+1
                                                                         ! multiplied by Lambda_3
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(INOUT)   :: PCLDFR    ! Cloud fraction          
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(INOUT)   :: PICEFR    ! Cloud fraction          
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),     INTENT(IN)    :: PRI_MF! Convective Mass Flux ice mixing ratio
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
!
END SUBROUTINE LIMA_ADJUST_SPLIT
END INTERFACE
END MODULE MODI_LIMA_ADJUST_SPLIT
