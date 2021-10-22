!     ######spl
     MODULE MODI_MF_TURB_EXPL
!    ######################
!
INTERFACE
!     #################################################################
      SUBROUTINE MF_TURB_EXPL(KKA,KKB,KKE,KKU,KKL,OMIXUV,             &
                PRHODJ,                                               &
                PTHLM,PTHVM,PRTM,PUM,PVM,                             &
                PTHLDT,PRTDT,PUDT,PVDT,                               &
                PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,              &
                PFLXZTHLMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF)

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

REAL, DIMENSION(:,:), INTENT(IN)   :: PRHODJ      ! dry density * Grid size

!   Conservative var. at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHLM        ! conservative pot. temp.
REAL, DIMENSION(:,:), INTENT(IN) ::  PRTM         ! water var.  where 
!REAL, DIMENSION(:,:), INTENT(IN) ::  PRVM  
!  Virtual potential temperature at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHVM 
!  Potential temperature at t-dt
!REAL, DIMENSION(:,:), INTENT(IN) ::  PTHM 
!  Momentum at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PUM
REAL, DIMENSION(:,:), INTENT(IN) ::  PVM
!
! Tendencies of conservative variables
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PTHLDT
!REAL, DIMENSION(:,:),   INTENT(OUT) ::  PTHVDT
!REAL, DIMENSION(:,:),   INTENT(OUT) ::  PTHDT

REAL, DIMENSION(:,:),   INTENT(OUT) ::  PRTDT 
!REAL, DIMENSION(:,:),   INTENT(OUT) ::  PRVDT 
! Tendencies of momentum
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PUDT
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PVDT

! Updraft characteritics
REAL, DIMENSION(:,:), INTENT(IN)   ::  PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP
!REAL, DIMENSION(:,:), INTENT(IN)   ::  PRV_UP
! Fluxes
REAL, DIMENSION(:,:), INTENT(OUT)  ::  PFLXZTHLMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF
!REAL, DIMENSION(:,:), INTENT(OUT)  ::  PFLXZTHMF

END SUBROUTINE MF_TURB_EXPL

END INTERFACE
!
END MODULE MODI_MF_TURB_EXPL
