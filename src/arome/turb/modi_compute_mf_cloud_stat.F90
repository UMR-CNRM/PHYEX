!     ######spl
     MODULE MODI_COMPUTE_MF_CLOUD_STAT
!    ############################
!
INTERFACE
!     #################################################################
      SUBROUTINE COMPUTE_MF_CLOUD_STAT(KKA, KKB, KKE, KKU, KKL, KRR, KRRL, KRRI,&
                            &PFRAC_ICE,&
                            &PTHLM, PRTM, PPABSM, PRM,&
                            &PDZZ, PTHM, PEXNM,&
                            &PEMF, PTHL_UP, PRT_UP,&
                            &PSIGMF)
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
INTEGER,                INTENT(IN)   :: KRR                     ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL                    ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI                    ! number of ice water var.
REAL, DIMENSION(:,:),   INTENT(IN)   :: PFRAC_ICE               ! liquid/ice fraction
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHLM, PRTM             ! cons. var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   :: PPABSM                  ! Pressure at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   :: PRM                     ! water var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   :: PDZZ
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHM                    ! environement
REAL, DIMENSION(:,:),   INTENT(IN)   :: PEXNM
REAL, DIMENSION(:,:),   INTENT(IN)   :: PEMF                    ! updraft characteritics
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHL_UP, PRT_UP         ! rc,w,Mass Flux,Thetal,rt
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PSIGMF                  ! SQRT(variance) for statistical cloud scheme


END SUBROUTINE COMPUTE_MF_CLOUD_STAT

END INTERFACE
!
END MODULE MODI_COMPUTE_MF_CLOUD_STAT
