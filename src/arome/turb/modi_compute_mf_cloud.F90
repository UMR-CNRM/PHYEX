!     ######spl
     MODULE MODI_COMPUTE_MF_CLOUD
!    ############################
!
INTERFACE
!     #################################################################
      SUBROUTINE COMPUTE_MF_CLOUD(KKA,KKB,KKE,KKU,KKL,KRR,KRRL,KRRI,HMF_CLOUD,&
                                  PFRAC_ICE,                                &
                                  PRC_UP,PRI_UP,PEMF,                       &
                                  PTHL_UP, PRT_UP, PFRAC_UP,                &
                                  PTHV_UP, PFRAC_ICE_UP, PRSAT_UP,          &
                                  PEXNM, PTHLM, PRTM, PTHM, PTHVM, PRM,     &
                                  PDZZ, PZZ, KKLCL,                         &
                                  PPABSM, PRHODREF,                         &
                                  PRC_MF, PRI_MF, PCF_MF, PSIGMF, PDEPTH    )
!     #################################################################
!!
!               
!*               1.1  Declaration of Arguments
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   ::  KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
INTEGER,                INTENT(IN)   ::  KRRL         ! number of liquid water var.
INTEGER,                INTENT(IN)   ::  KRRI         ! number of ice water var.
CHARACTER (LEN=4),      INTENT(IN)   ::  HMF_CLOUD    ! Type of statistical cloud scheme
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PFRAC_ICE    ! liquid/ice fraction
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PRC_UP,PRI_UP,PEMF ! updraft characteritics
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHL_UP, PRT_UP  
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PFRAC_UP
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHV_UP           ! updraft thetaV
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PFRAC_ICE_UP      ! liquid/solid fraction in updraft
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PRSAT_UP          ! Rsat in updraft
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PEXNM             ! exner function
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHLM, PRTM ! cons. var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHM, PTHVM ! theta and thetaV
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRM         ! water var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDZZ, PZZ
INTEGER, DIMENSION(:),  INTENT(IN)   ::  KKLCL       ! index of updraft condensation level
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PPABSM, PRHODREF ! environement
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PRC_MF, PRI_MF   ! cloud content and
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PCF_MF           ! cloud fraction for MF scheme
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PSIGMF ! SQRT(variance) for statistical cloud scheme
REAL, DIMENSION(:),     INTENT(IN)   ::  PDEPTH ! Deepness of cloud

END SUBROUTINE COMPUTE_MF_CLOUD

END INTERFACE
!
END MODULE MODI_COMPUTE_MF_CLOUD
