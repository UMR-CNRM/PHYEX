!     ######spl
     MODULE MODI_SHALLOW_MF
!    ######################
!
INTERFACE
!     #################################################################
      SUBROUTINE SHALLOW_MF(KKA,KKU,KKL,KRR,KRRL,KRRI,                &
                HMF_UPDRAFT, HMF_CLOUD, HFRAC_ICE, OMIXUV,            &
                ONOMIXLG,KSV_LGBEG,KSV_LGEND,                         &
                PIMPL_MF, PTSTEP, PTSTEP_MET, PTSTEP_SV,              &
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
                KKLCL,KKETL,KKCTL                                     )
!     #################################################################
!!
!               
!*               1.1  Declaration of Arguments
!                
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER,                INTENT(IN)   :: KRR          ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL         ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI         ! number of ice water var.
CHARACTER (LEN=4),      INTENT(IN)   :: HMF_UPDRAFT  ! Type of Mass Flux Scheme
                                     ! 'NONE' if no parameterization 
CHARACTER (LEN=4),      INTENT(IN)   :: HMF_CLOUD    ! Type of statistical cloud
                                                     ! scheme
CHARACTER*1,            INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
LOGICAL,                INTENT(IN)   :: OMIXUV    ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL,                   INTENT(IN)   :: PIMPL_MF     ! degre of implicitness
REAL,              INTENT(IN)     ::  PTSTEP   ! Dynamical timestep 
REAL,              INTENT(IN)     ::  PTSTEP_MET! Timestep for meteorological variables                        
REAL,              INTENT(IN)     ::  PTSTEP_SV! Timestep for tracer variables

REAL, DIMENSION(:,:),   INTENT(IN) ::  PZZ         ! Height of flux point
REAL, DIMENSION(:,:),   INTENT(IN) ::  PDZZ        ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN) ::  PRHODJ      ! dry density * Grid size
REAL, DIMENSION(:,:),   INTENT(IN) ::  PRHODREF    ! dry density of the
                                                     ! reference state
REAL, DIMENSION(:,:),   INTENT(IN) ::  PPABSM      ! Pressure at time t-1
REAL, DIMENSION(:,:),   INTENT(IN) ::  PEXNM       ! Exner function at t-dt

REAL, DIMENSION(:),     INTENT(IN) ::  PSFTH,PSFRV ! normal surface fluxes of theta and Rv 
REAL, DIMENSION(:,:),   INTENT(IN) ::  PTHM        ! Theta at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PRM         ! water var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN) ::  PUM,PVM     ! wind components at t-dt
REAL, DIMENSION(:,:),   INTENT(IN) ::  PTKEM       ! tke at t-dt

REAL, DIMENSION(:,:,:), INTENT(IN) ::  PSVM        ! scalar variable a t-dt

REAL, DIMENSION(:,:),   INTENT(OUT)::  PDUDT_MF     ! tendency of U   by massflux scheme
REAL, DIMENSION(:,:),   INTENT(OUT)::  PDVDT_MF     ! tendency of V   by massflux scheme
REAL, DIMENSION(:,:),   INTENT(OUT)::  PDTHLDT_MF   ! tendency of thl by massflux scheme
REAL, DIMENSION(:,:),   INTENT(OUT)::  PDRTDT_MF    ! tendency of rt  by massflux scheme
REAL, DIMENSION(:,:,:), INTENT(OUT)::  PDSVDT_MF    ! tendency of Sv  by massflux scheme

REAL, DIMENSION(:,:),   INTENT(OUT)   ::  PSIGMF,PRC_MF,PRI_MF,PCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(:,:), INTENT(OUT)     ::  PFLXZTHVMF           ! Thermal production for TKE scheme
REAL, DIMENSION(:,:), INTENT(OUT)     ::  PFLXZTHMF
REAL, DIMENSION(:,:), INTENT(OUT)     ::  PFLXZRMF
REAL, DIMENSION(:,:), INTENT(OUT)     ::  PFLXZUMF
REAL, DIMENSION(:,:), INTENT(OUT)     ::  PFLXZVMF
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PTHL_UP   ! Thl updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PRT_UP    ! Rt  updraft characteristics
REAL, DIMENSION(:,:), INTENT(OUT) ::  PRV_UP    ! Vapor updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PU_UP     ! U wind updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PV_UP     ! V wind updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PRC_UP    ! cloud content updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PRI_UP    ! ice content   updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PTHV_UP   ! Thv   updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PW_UP     ! vertical speed updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PFRAC_UP  ! updraft fraction
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PEMF      ! updraft mass flux
REAL, DIMENSION(:,:), INTENT(OUT) ::  PDETR     ! updraft detrainment
REAL, DIMENSION(:,:), INTENT(OUT) ::  PENTR     ! updraft entrainment
INTEGER,DIMENSION(:), INTENT(OUT) :: KKLCL,KKETL,KKCTL ! level of LCL,ETL and CTL


END SUBROUTINE SHALLOW_MF

END INTERFACE
!
END MODULE MODI_SHALLOW_MF
