!     ######spl
      MODULE MODI_CH_SIS
!!    ################## 
!!
!
INTERFACE
!!
SUBROUTINE CH_SIS(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI)
IMPLICIT NONE
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT  
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC ! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN)  :: KMI      ! model number
END SUBROUTINE CH_SIS
!!
END INTERFACE
!!
END MODULE MODI_CH_SIS
