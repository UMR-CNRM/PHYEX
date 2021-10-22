!     ######spl
      MODULE MODI_CH_CRANCK
!!    ##################### 
!!
!
INTERFACE
!!
SUBROUTINE CH_CRANCK(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                       PALPHA)
IMPLICIT NONE
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC ! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN) :: KMI      ! model number
REAL, INTENT(IN) :: PALPHA
END SUBROUTINE CH_CRANCK
!!
END INTERFACE
!!
END MODULE MODI_CH_CRANCK
