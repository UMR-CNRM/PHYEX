!     ######spl
     MODULE MODI_COMPUTE_MF_CLOUD_DIRECT
!    ###################################
!
INTERFACE
!     #################################################################
      SUBROUTINE COMPUTE_MF_CLOUD_DIRECT(KKB, KKE, KKL, &
                                        &KKLCL, PFRAC_UP, PRC_UP, PRI_UP,&
                                        &PRC_MF, PRI_MF, PCF_MF)
!     #################################################################
!!
!
!*               1.1  Declaration of Arguments
!
INTEGER,                INTENT(IN)   :: KKB            ! near groud physical index
INTEGER,                INTENT(IN)   :: KKE            ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKL            ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER, DIMENSION(:),  INTENT(IN)   :: KKLCL          ! index of updraft condensation level
REAL, DIMENSION(:,:),   INTENT(IN)   :: PFRAC_UP       ! Updraft Fraction
REAL, DIMENSION(:,:),   INTENT(IN)   :: PRC_UP,PRI_UP  ! updraft characteritics
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PRC_MF, PRI_MF ! cloud content (INPUT=environment, OUTPUT=conv. cloud)
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PCF_MF         ! and cloud fraction for MF scheme

END SUBROUTINE COMPUTE_MF_CLOUD_DIRECT

END INTERFACE
!
END MODULE MODI_COMPUTE_MF_CLOUD_DIRECT
