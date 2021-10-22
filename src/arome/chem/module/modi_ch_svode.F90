!     ######spl
      MODULE MODI_CH_SVODE
!!    #############################

!
INTERFACE

SUBROUTINE CH_SVODE(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                    PRTOL, PATOL, KPED )
IMPLICIT NONE
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC
                 ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC
                 ! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN) :: KMI      ! model index
REAL,    INTENT(IN) :: PRTOL, PATOL
INTEGER, INTENT(IN) :: KPED
END SUBROUTINE CH_SVODE

END INTERFACE

END MODULE MODI_CH_SVODE
