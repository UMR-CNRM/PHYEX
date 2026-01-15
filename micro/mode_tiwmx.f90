!@no_insert_drhook
!     ######spl
      MODULE MODE_TIWMX
!     ###############
!
!!****  *MODE_TIWMX* - 
!!
!!    PURPOSE
!!    -------
!       The purpose of this  ...
!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (ha ha)
!!          
!!    AUTHOR
!!    ------
!!      K. I. Ivarsson   *SMHI*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/11/14  
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
IMPLICIT NONE

CONTAINS

  REAL FUNCTION ESATW(TIWMX, TT)
    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    ESATW = TIWMX%ESTABW(NINT(XNDEGR*TT))
  END FUNCTION ESATW

  REAL FUNCTION DESDTW(TIWMX, TT)
    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    DESDTW = TIWMX%DESTABW(NINT(XNDEGR*TT))
  END FUNCTION

  REAL FUNCTION ESATI(TIWMX, TT)
    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    ESATI = TIWMX%ESTABI(NINT(XNDEGR*TT))
  END FUNCTION

  REAL FUNCTION DESDTI(TIWMX, TT)
    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    DESDTI = TIWMX%DESTABI(NINT(XNDEGR*TT))
  END FUNCTION

! Water droplet function:
  REAL FUNCTION AA2W(TIWMX, TT)
    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    AA2W = TIWMX%A2WTAB(NINT(XNDEGR*TT))
  END FUNCTION

! Ice crystal function
  PURE REAL FUNCTION AA2(TIWMX, TT)
    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    AA2 = TIWMX%A2TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Meyers IN concentration function:
  PURE REAL FUNCTION AM3(TIWMX, TT)
    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    AM3 = TIWMX%AM3TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Fletchers IN concentration function:
  PURE REAL FUNCTION AF3(TIWMX, TT)
    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    AF3 = TIWMX%AF3TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Ice crystal function
  PURE REAL FUNCTION BB3(TIWMX, TT)
    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    BB3 = TIWMX%BB3TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Water droplet function:
  REAL FUNCTION BB3W(TIWMX, TT)
    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    BB3W = TIWMX%BB3WTAB(NINT(XNDEGR*TT))
  END FUNCTION

! Function for IN concentration reduction between 0 and -25 C:
  PURE REAL FUNCTION REDIN(TIWMX, TT)
    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    REDIN = TIWMX%REDINTAB(NINT(XNDEGR*TT))
  END FUNCTION
END MODULE MODE_TIWMX
