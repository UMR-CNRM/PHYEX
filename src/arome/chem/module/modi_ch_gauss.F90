!     ######spl
      MODULE MODI_CH_GAUSS
!!    ############################# 

INTERFACE

SUBROUTINE CH_GAUSS(PIN,POUT,KDIM,KFAIL)
IMPLICIT NONE
INTEGER, INTENT(IN)                     :: KDIM   ! dimension of the matrix
INTEGER, INTENT(INOUT)                  :: KFAIL  ! error flag
REAL, DIMENSION(KDIM,KDIM), INTENT(IN)  :: PIN    ! matrix to be inverted
REAL, DIMENSION(KDIM,KDIM), INTENT(OUT) :: POUT   ! inverse matrix
END SUBROUTINE CH_GAUSS

END INTERFACE

END MODULE MODI_CH_GAUSS
