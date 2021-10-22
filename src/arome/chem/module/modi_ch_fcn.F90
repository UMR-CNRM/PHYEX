!     ######spl
      MODULE MODI_CH_FCN
!!    ##################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
INTERFACE
SUBROUTINE CH_FCN(PTIME,PCONC,PDCDT,KMI,KVECNPT,KEQ)
IMPLICIT NONE
REAL,    INTENT(IN)                          :: PTIME
INTEGER, INTENT(IN)                          :: KVECNPT
INTEGER, INTENT(IN)                          :: KEQ
REAL,    INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC
REAL,    INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PDCDT
INTEGER, INTENT(IN)                          :: KMI
!!
END SUBROUTINE CH_FCN
END INTERFACE
END MODULE MODI_CH_FCN
