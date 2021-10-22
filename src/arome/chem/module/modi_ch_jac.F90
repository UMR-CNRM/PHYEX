!     ######spl
      MODULE MODI_CH_JAC
!!    ##################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
INTERFACE
SUBROUTINE CH_JAC(PTIME,PCONC,PJAC,KMI,KVECNPT,KEQ)
IMPLICIT NONE
REAL,    INTENT(IN)                              :: PTIME
INTEGER, INTENT(IN)                              :: KVECNPT
INTEGER, INTENT(IN)                              :: KEQ
REAL,    INTENT(IN),  DIMENSION(KVECNPT,KEQ)     :: PCONC
REAL,    INTENT(OUT), DIMENSION(KVECNPT,KEQ,KEQ) :: PJAC
INTEGER, INTENT(IN)                              :: KMI
END SUBROUTINE CH_JAC
END INTERFACE
END MODULE MODI_CH_JAC
