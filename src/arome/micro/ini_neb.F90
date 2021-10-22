!     ######spl
        SUBROUTINE INI_NEB
!       #######################
!
!!****     *INI_NEB*  - routine to initialize the nebulosity computation
!!                        constants.
!!
!!      PURPOSE
!!      -------
!         The purpose of this routine is to initialize
!       constants used for nebulosity computation
!
!!      METHOD
!!      ------
!!        The constants are set to their numerical values
!!
!!      EXTERNAL
!!      --------
!!        NONE
!!
!!      IMPLICIT ARGUMENTS
!!      ------------------
!!        Module MODD_NEB
!!
!!      REFERENCE
!!      ---------
!!
!!      AUTHOR
!!      ------
!!        S. Riette (Meteo France)
!!
!!      MODIFICATIONS
!!      -------------
!!        Original 24 Aug 2011
!! --------------------------------------------------------------------------
!
!*        0. DECLARATIONS
!            ------------
!
USE MODD_NEB
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!  ---------------------------------------------------------------------------
!
!         1. SETTING THE NUMERICAL VALUES
!
IF (LHOOK) CALL DR_HOOK('INI_NEB',0,ZHOOK_HANDLE)
!Freezing between 0 and -20. Other possibilities are 0/-40 or -5/-25
XTMAXMIX    = 273.16
XTMINMIX    = 253.16
IF (LHOOK) CALL DR_HOOK('INI_NEB',1,ZHOOK_HANDLE)
END SUBROUTINE INI_NEB
