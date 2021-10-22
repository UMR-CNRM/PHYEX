!     ######spl
     MODULE MODI_COMPUTE_MF_CLOUD_BIGAUS
!    ###################################
!
INTERFACE
!     #################################################################
      SUBROUTINE COMPUTE_MF_CLOUD_BIGAUS(KKA, KKB, KKE, KKU, KKL,&
                                  PRC_UP, PRI_UP, PEMF, PDEPTH,&
                                  PRT_UP, PTHV_UP, PFRAC_ICE_UP, PRSAT_UP,&
                                  PRTM, PTHM, PTHVM,&
                                  PDZZ, PZZ, PRHODREF,&
                                  PRC_MF, PRI_MF, PCF_MF)
!     #################################################################
!!
!
!*               1.1  Declaration of Arguments
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL                     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:),   INTENT(IN)   :: PRC_UP,PRI_UP,PEMF      ! updraft characteritics
REAL, DIMENSION(:),     INTENT(IN)   :: PDEPTH                  ! Deepness of cloud
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHV_UP, PRSAT_UP, PRT_UP ! updraft characteritics
REAL, DIMENSION(:,:),   INTENT(IN)   :: PFRAC_ICE_UP            ! liquid/ice fraction in updraft
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHM, PRTM, PTHVM       ! env. var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   :: PDZZ, PZZ
REAL, DIMENSION(:,:),   INTENT(IN)   :: PRHODREF
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PRC_MF, PRI_MF          ! cloud content
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PCF_MF                  ! and cloud fraction for MF scheme

END SUBROUTINE COMPUTE_MF_CLOUD_BIGAUS

END INTERFACE
!
END MODULE MODI_COMPUTE_MF_CLOUD_BIGAUS
