!     ######spl
      SUBROUTINE CH_GET_CNAMES(HNAMES)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #################################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
!!
!!*** *CH_GET_CNAMES*
!!
!!    PURPOSE
!!    -------
!       return the names for the chemical species in HNAMES
!!
!!**  METHOD
!!    ------
!!      simple
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Karsten Suhre (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 26/07/96
!!    Modified 05/05/98: Vectorization (Vincent Crassier & KS)
!!    Modified 31/10/03: New interface for better MesoNH compilation (D. Gazen)
!!
!!----------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_M9_SCHEME
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
!!
!!    LOCAL VARIABLES
!!    ---------------
CHARACTER(LEN=32), DIMENSION(:), INTENT(OUT) :: HNAMES
!!
!!----------------------------------------------------------------------
!!
! copy the names of the chemical species into HNAMES
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_GET_CNAMES',0,ZHOOK_HANDLE)
HNAMES = CNAMES
IF (LHOOK) CALL DR_HOOK('CH_GET_CNAMES',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE CH_GET_CNAMES
