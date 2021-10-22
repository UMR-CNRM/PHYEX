!     ######spl
     MODULE MODI_MF_TURB
!    ######################
!
INTERFACE
!     #################################################################
      SUBROUTINE MF_TURB(KKA,KKB,KKE,KKU,KKL,OMIXUV,                  &
                ONOMIXLG,KSV_LGBEG,KSV_LGEND,                         &
                PIMPL, PTSTEP, PTSTEP_MET, PTSTEP_SV,                 &
                PDZZ,                                                 &
                PRHODJ,                                               &
                PTHLM,PTHVM,PRTM,PUM,PVM,PSVM,                        &
                PTHLDT,PRTDT,PUDT,PVDT,PSVDT,                         &
                PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,PSV_UP,       &
                PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF,      &
                PFLXZSVMF                                             )

!     #################################################################
!
!               
!*               1.1  Declaration of Arguments
!
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise

LOGICAL,                INTENT(IN)   :: OMIXUV      ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL,                   INTENT(IN)   :: PIMPL       ! degree of implicitness
REAL,                 INTENT(IN)     ::  PTSTEP   ! Dynamical timestep 
REAL,                 INTENT(IN)     ::  PTSTEP_MET! Timestep for meteorological variables                        
REAL,                 INTENT(IN)     ::  PTSTEP_SV! Timestep for tracer variables
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PDZZ        ! metric coefficients

REAL, DIMENSION(:,:), INTENT(IN)   ::  PRHODJ    ! dry density * Grid size

!   Conservative var. at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHLM       ! conservative pot. temp.
REAL, DIMENSION(:,:), INTENT(IN) ::  PRTM         ! water var.  where 
!  Virtual potential temperature at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHVM 
!  Momentum at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PUM
REAL, DIMENSION(:,:), INTENT(IN) ::  PVM
!  scalar variables at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PSVM
!
! Tendencies of conservative variables
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PTHLDT

REAL, DIMENSION(:,:),   INTENT(OUT) ::  PRTDT 
! Tendencies of momentum
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PUDT
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PVDT
! Tendencies of scalar variables
REAL, DIMENSION(:,:,:), INTENT(OUT) ::  PSVDT


! Updraft characteritics
REAL, DIMENSION(:,:), INTENT(IN)   ::  PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PSV_UP
! Fluxes
REAL, DIMENSION(:,:), INTENT(OUT)  ::  PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF

REAL, DIMENSION(:,:,:), INTENT(OUT)::  PFLXZSVMF

END SUBROUTINE MF_TURB

END INTERFACE
!
END MODULE MODI_MF_TURB
