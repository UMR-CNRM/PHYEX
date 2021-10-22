!     ######spl
      MODULE MODI_CH_SET_RATES
!!    ########################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
INTERFACE
SUBROUTINE CH_SET_RATES(PTIME,PCONC,TPM,KMI,KOUT,KVERB,KVECNPT,KEQ)
USE MODD_CH_M9, ONLY: METEOTRANSTYPE
IMPLICIT NONE
REAL,    INTENT(IN)                      :: PTIME
INTEGER, INTENT(IN)                      :: KVECNPT
INTEGER, INTENT(IN)                      :: KEQ
REAL,    INTENT(IN),  DIMENSION(KVECNPT,KEQ)        :: PCONC
TYPE(METEOTRANSTYPE), DIMENSION(KVECNPT), INTENT(IN):: TPM
INTEGER, INTENT(IN)                      :: KMI
INTEGER, INTENT(IN)                      :: KOUT,KVERB
END SUBROUTINE CH_SET_RATES
END INTERFACE
END MODULE MODI_CH_SET_RATES
