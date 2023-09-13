!     ######spl
     MODULE MODI_SHALLOW_MF
!    ######################
!
IMPLICIT NONE
INTERFACE
!     #################################################################
      SUBROUTINE SHALLOW_MF(D, CST, NEBN, PARAMMF, TURBN, CSTURB,     &
                KRR, KRRL, KRRI, KSV,                                 &
                ONOMIXLG,KSV_LGBEG,KSV_LGEND,                         &
                PTSTEP,                                               &
                PDZZ, PZZ,                                            &
                PRHODJ, PRHODREF,                                     &
                PPABSM, PEXNM,                                        &
                PSFTH,PSFRV,                                          &
                PTHM,PRM,PUM,PVM,PTKEM,PSVM,                          &
                PDUDT_MF,PDVDT_MF,                                    &
                PDTHLDT_MF,PDRTDT_MF,PDSVDT_MF,                       &
                PSIGMF,PRC_MF,PRI_MF,PCF_MF,PFLXZTHVMF,               &
                PFLXZTHMF,PFLXZRMF,PFLXZUMF,PFLXZVMF,                 &
                PTHL_UP,PRT_UP,PRV_UP,PRC_UP,PRI_UP,                  &
                PU_UP, PV_UP, PTHV_UP, PW_UP,                         &
                PFRAC_UP,PEMF,PDETR,PENTR,                            &
                KKLCL,KKETL,KKCTL,PDX,PDY,PRSVS,PSVMIN,               &
                BUCONF, TBUDGETS, KBUDGETS                            )
!     #################################################################
!!
USE MODD_BUDGET,          ONLY: TBUDGETCONF_t, TBUDGETDATA
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST,             ONLY: CST_t
USE MODD_NEB_n,           ONLY: NEB_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
USE MODD_TURB_n,          ONLY: TURB_t
USE MODD_CTURB,           ONLY: CSTURB_t
USE MODD_PARAMETERS,      ONLY: JPSVMAX
IMPLICIT NONE
!               
!*               1.1  Declaration of Arguments
!                
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D            ! PHYEX variables dimensions structure
TYPE(CST_t),            INTENT(IN)   :: CST          ! modd_cst general constant structure
TYPE(NEB_t),            INTENT(IN)   :: NEBN
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF
TYPE(TURB_t),           INTENT(IN)   :: TURBN        ! modn_turbn (turb namelist) structure
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB       ! modd_csturb turb constant structure
INTEGER,                INTENT(IN)   :: KRR          ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL         ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI         ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV          ! number of scalar var.
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL,                   INTENT(IN)   :: PTSTEP    ! Dynamical timestep 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PZZ         ! Height of flux point
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PDZZ        ! Metric coefficients
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PRHODJ      ! dry density * Grid size
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PRHODREF    ! dry density of the
                                                     ! reference state
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PPABSM      ! Pressure at time t-1
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PEXNM       ! Exner function at t-dt

REAL, DIMENSION(D%NIJT),         INTENT(IN) ::  PSFTH,PSFRV ! normal surface fluxes of theta and Rv 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PTHM        ! Theta at t-dt
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN) ::  PRM         ! water var. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PUM,PVM     ! wind components at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PTKEM       ! tke at t-dt
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVM        ! scalar variable a t-dt

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)::  PDUDT_MF     ! tendency of U   by massflux scheme
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)::  PDVDT_MF     ! tendency of V   by massflux scheme
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)::  PDTHLDT_MF   ! tendency of thl by massflux scheme
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)::  PDRTDT_MF    ! tendency of rt  by massflux scheme
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT)::  PDSVDT_MF    ! tendency of Sv  by massflux scheme

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PSIGMF,PRC_MF,PRI_MF,PCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PFLXZTHVMF           ! Thermal production for TKE scheme
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PFLXZTHMF
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PFLXZRMF
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PFLXZUMF
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PFLXZVMF
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PTHL_UP   ! Thl updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PRT_UP    ! Rt  updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PRV_UP    ! Vapor updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PU_UP     ! U wind updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PV_UP     ! V wind updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PRC_UP    ! cloud content updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PRI_UP    ! ice content   updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PTHV_UP   ! Thv   updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PW_UP     ! vertical speed updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PFRAC_UP  ! updraft fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PEMF      ! updraft mass flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) ::  PDETR     ! updraft detrainment
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) ::  PENTR     ! updraft entrainment
INTEGER,DIMENSION(D%NIJT),     INTENT(OUT) :: KKLCL,KKETL,KKCTL ! level of LCL,ETL and CTL
REAL,                          INTENT(IN)  :: PDX, PDY
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN),OPTIONAL ::  PRSVS ! sources of sv (for Budgets with lagrangian tracer)
TYPE(TBUDGETCONF_t),    INTENT(IN),OPTIONAL   :: BUCONF       ! budget structure
INTEGER,                INTENT(IN)   :: KBUDGETS     ! option. because not used in arpifs
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT),OPTIONAL :: TBUDGETS
REAL,DIMENSION(JPSVMAX),INTENT(IN),OPTIONAL   :: PSVMIN       ! minimum value for SV variables (for Budgets)

END SUBROUTINE SHALLOW_MF

END INTERFACE
!
END MODULE MODI_SHALLOW_MF
