!     ######spl
    MODULE MODI_CH_EXQSSA
!!    ##################### 
!!
!
INTERFACE
!!
SUBROUTINE CH_EXQSSA(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                       PSLOW, PFAST, PDTMAX)
IMPLICIT NONE
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    
                              ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC 
                              ! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN) :: KMI    ! model number
                              ! reac. rates, auxiliary variables
REAL,    INTENT(IN) :: PSLOW, PFAST, PDTMAX
END SUBROUTINE CH_EXQSSA
!!
END INTERFACE
!!
END MODULE MODI_CH_EXQSSA
