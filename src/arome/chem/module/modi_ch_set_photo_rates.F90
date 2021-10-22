!     ######spl
      MODULE MODI_CH_SET_PHOTO_RATES
!!    ##############################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
INTERFACE
SUBROUTINE CH_SET_PHOTO_RATES(PTIME,PCONC,KL,TPM,KMI,KOUT,KVERB,KVECNPT,KVECMASK,KEQ,PJVALUES)
USE MODD_CH_M9, ONLY: METEOTRANSTYPE
IMPLICIT NONE
REAL,    INTENT(IN)                              :: PTIME
INTEGER, INTENT(IN)                              :: KVECNPT,KL,KEQ,KMI
INTEGER, DIMENSION(:,:), INTENT(IN)              :: KVECMASK
REAL,    INTENT(IN),  DIMENSION(KVECNPT,KEQ)     :: PCONC
TYPE(METEOTRANSTYPE), DIMENSION(KVECNPT), INTENT(IN) :: TPM
INTEGER, INTENT(IN)                              :: KOUT,KVERB
REAL,DIMENSION(:,:,:,:), INTENT(IN) :: PJVALUES    ! Tuv coefficient
END SUBROUTINE CH_SET_PHOTO_RATES
END INTERFACE
END MODULE MODI_CH_SET_PHOTO_RATES
