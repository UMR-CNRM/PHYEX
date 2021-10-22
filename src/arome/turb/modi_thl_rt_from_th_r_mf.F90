!     ######spl
     MODULE MODI_THL_RT_FROM_TH_R_MF
!    ###############################
!
INTERFACE
!     #################################################################
      SUBROUTINE THL_RT_FROM_TH_R_MF( KRR,KRRL,KRRI,                  &
                                      PTH, PR, PEXN, &
                                      PTHL, PRT                      )
!     #################################################################
!
!               
!*               1.1  Declaration of Arguments
!                
!
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.

REAL, DIMENSION(:,:), INTENT(IN)   :: PTH      ! theta
REAL, DIMENSION(:,:,:), INTENT(IN) :: PR       ! water species
REAL, DIMENSION(:,:), INTENT(IN)   :: PEXN    ! exner function

REAL, DIMENSION(:,:), INTENT(OUT)  :: PTHL     ! th_l
REAL, DIMENSION(:,:), INTENT(OUT)  :: PRT      ! total non precip. water
!
END SUBROUTINE THL_RT_FROM_TH_R_MF

END INTERFACE
!
END MODULE MODI_THL_RT_FROM_TH_R_MF
