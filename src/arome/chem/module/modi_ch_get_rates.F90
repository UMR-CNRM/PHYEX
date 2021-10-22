!     ######spl
      MODULE MODI_CH_GET_RATES
!!    ########################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
INTERFACE
SUBROUTINE CH_GET_RATES(PRATE,KMI,KVECNPT,KREAC)
IMPLICIT NONE
INTEGER, INTENT(IN)                         :: KVECNPT
INTEGER, INTENT(IN)                         :: KREAC
REAL, INTENT(OUT), DIMENSION(KVECNPT,KREAC) :: PRATE
INTEGER, INTENT(IN)                         :: KMI
END SUBROUTINE CH_GET_RATES
END INTERFACE
END MODULE MODI_CH_GET_RATES
