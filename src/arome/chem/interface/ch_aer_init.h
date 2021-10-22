INTERFACE
     SUBROUTINE CH_AER_INIT(PCHEM,PAERO, PRHODREF)
!!   ############################################################
!!
!!
!!    EXTERNAL
!!    --------
!!    None
!!
USE PARKIND1  ,ONLY : JPIM     ,JPRB
!!
IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL(KIND=JPRB),   DIMENSION(:,:,:,:),    INTENT(INOUT)   :: PCHEM, PAERO
REAL(KIND=JPRB),   DIMENSION(:,:,:),      INTENT(IN)      :: PRHODREF
!
END SUBROUTINE CH_AER_INIT
END INTERFACE
