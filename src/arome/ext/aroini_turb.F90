!     ######spl
SUBROUTINE AROINI_TURB(PLINI,OHARATU,OSTATNW,OSUBG_COND,HCONDENS,HLAMBDA3,HSUBG_MF_PDF,OSIGMAS)
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

USE MODD_LES,   ONLY : TLES, LES_ASSOCIATE
USE MODD_CTURB, ONLY : XLINI
USE MODD_TURB_n, ONLY: LHARAT, LSTATNW, CTURBLEN, TURB_GOTO_MODEL, LTURB_FLX, LTURB_DIAG, &
                       LSUBG_COND, LRMC01, CTURBDIM, XIMPL, CTOM, CCONDENS, CLAMBDA3,  &
                       CSUBG_MF_PDF, LSIGMAS
USE MODI_INI_CTURB

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL,   INTENT(IN) :: PLINI ! minimum bl89 mixing length
LOGICAL,INTENT(IN) :: OHARATU ! switch HARATU
LOGICAL,INTENT(IN) :: OSTATNW ! switch LSTATNW
LOGICAL,INTENT(IN) :: OSUBG_COND ! switch of subgrid condensation
CHARACTER(LEN=80),INTENT(IN)   :: HCONDENS ! subrgrid condensation PDF
CHARACTER(LEN=4),INTENT(IN)    :: HLAMBDA3 ! lambda3 choice for subgrid cloud scheme
CHARACTER(LEN=80),INTENT(IN)  :: HSUBG_MF_PDF ! PDF to use for MF cloud autoconversions
LOGICAL, INTENT(IN) :: OSIGMAS
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
LSTATNW=OSTATNW

!         2. Set implicit default values for MODD_LES

CALL LES_ASSOCIATE()
TLES%LLES=.FALSE.
TLES%LLES_CALL=.FALSE.

!         3. Set implicit default values for MODD_TURB_n

CTURBLEN  = 'BL89'
CTURBDIM  = '1DIM'
LTURB_FLX = .FALSE.
LTURB_DIAG = .FALSE.
LSIGMAS=OSIGMAS
XIMPL = 1.
LSUBG_COND = OSUBG_COND
CCONDENS=HCONDENS
CLAMBDA3=HLAMBDA3
CSUBG_MF_PDF=HSUBG_MF_PDF
LRMC01 = .FALSE.
CTOM = 'NONE'

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AROINI_TURB',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE AROINI_TURB
