!     ######spl
SUBROUTINE AROINI_TURB(PLINI,OHARATU)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!**** *INI_TURB*   - Initialize common meso_NH MODD_ used in Turbulence scheme
!                    for AROME
!     Purpose.
!     --------
!           Initialize MODD_LES and  MODD_TKE
!           parameters used in AROME turbulence scheme 

!**   Interface.
!     ----------
!        *CALL* *INI_TURB

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
!        Y. Seity 

!     Modifications.
!     --------------
!        Original : 03-12-12
!     ------------------------------------------------------------------

USE MODD_LES,   ONLY : LLES, LLES_CALL
USE MODD_CTURB, ONLY : XLINI
USE MODD_TURB_n, ONLY: LHARAT, TURB_GOTO_MODEL
USE MODI_INI_CTURB

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL,   INTENT(IN) :: PLINI ! minimum bl89 mixing length
LOGICAL,INTENT(IN) :: OHARATU ! switch HARATU
!
!     ------------------------------------------------------------------

!         1. Set implicit default values for MODD_CTURB

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AROINI_TURB',0,ZHOOK_HANDLE)
!
CALL TURB_GOTO_MODEL(1,1)
!
CALL INI_CTURB

!         1bis. Modification of MODD_CTURB values
XLINI=PLINI
LHARAT=OHARATU

!         2. Set implicit default values for MODD_LES

LLES=.FALSE.
LLES_CALL=.FALSE.


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AROINI_TURB',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE AROINI_TURB
