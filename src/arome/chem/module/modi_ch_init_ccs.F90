!     ######spl
      MODULE MODI_CH_INIT_CCS
!     #######################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
INTERFACE
SUBROUTINE CH_INIT_CCS(KEQ,KOUT,KVERB)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: KEQ         ! number of scalar variables
INTEGER, INTENT(IN)  :: KOUT, KVERB ! stdout output, verbosity level 
END SUBROUTINE CH_INIT_CCS
END INTERFACE
END MODULE MODI_CH_INIT_CCS
