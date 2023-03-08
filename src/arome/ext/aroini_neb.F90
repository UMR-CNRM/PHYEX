SUBROUTINE AROINI_NEB
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!**** *AROINI_NEB*   -
!     Purpose.
!     --------
!           Call Meso-NH routine INI_NEB
!              (setup of constants for nebulosity computation)
!
!**   Interface.
!     ----------
!        *CALL* *AROINI_NEB

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation AROME 

!     Author.
!     -------
!        S. Riette

!     Modifications.
!     --------------
!        Original : 24 Aug 2011
!     ------------------------------------------------------------------

USE MODI_INI_NEB

IMPLICIT NONE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

!         1. Set implicit default values for MODD_NEB

IF (LHOOK) CALL DR_HOOK('AROINI_NEB',0,ZHOOK_HANDLE)
CALL INI_NEB
IF (LHOOK) CALL DR_HOOK('AROINI_NEB',1,ZHOOK_HANDLE)

END SUBROUTINE AROINI_NEB
