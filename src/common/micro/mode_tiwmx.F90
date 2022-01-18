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
IMPLICIT NONE

REAL, PARAMETER :: XNDEGR = 100.0
INTEGER, PARAMETER :: NSTART = 10000
INTEGER, PARAMETER :: NSTOP = 37316

! Saturation tables and derivatives
REAL ::  ESTABW(NSTART:NSTOP)
REAL :: DESTABW(NSTART:NSTOP)
REAL ::  ESTABI(NSTART:NSTOP)
REAL :: DESTABI(NSTART:NSTOP)

! Ice crystal- or water droplet tables
REAL ::   A2TAB(NSTART:NSTOP)
REAL ::  BB3TAB(NSTART:NSTOP)
REAL ::  AM3TAB(NSTART:NSTOP)
REAL ::  AF3TAB(NSTART:NSTOP)
REAL ::  A2WTAB(NSTART:NSTOP)
REAL :: BB3WTAB(NSTART:NSTOP)
REAL :: REDINTAB(NSTART:NSTOP)

CONTAINS

  REAL FUNCTION ESATW(TT)
    REAL,INTENT(IN) :: TT 
    ESATW = ESTABW(NINT(XNDEGR*TT))
  END FUNCTION ESATW

  REAL FUNCTION DESDTW(TT)
    REAL,INTENT(IN) :: TT 
    DESDTW = DESTABW(NINT(XNDEGR*TT))
  END FUNCTION

  REAL FUNCTION ESATI(TT)
    REAL,INTENT(IN) :: TT 
    ESATI = ESTABI(NINT(XNDEGR*TT))
  END FUNCTION

  REAL FUNCTION DESDTI(TT)
    REAL,INTENT(IN) :: TT 
    DESDTI = DESTABI(NINT(XNDEGR*TT))
  END FUNCTION

! Water droplet function:
  REAL FUNCTION AA2W(TT)
    REAL,INTENT(IN) :: TT 
    AA2W = A2WTAB(NINT(XNDEGR*TT))
  END FUNCTION

! Ice crystal function
  REAL FUNCTION AA2(TT)
    REAL,INTENT(IN) :: TT 
    AA2 = A2TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Meyers IN concentration function:
  REAL FUNCTION AM3(TT)
    REAL,INTENT(IN) :: TT 
    AM3 = AM3TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Fletchers IN concentration function:
  REAL FUNCTION AF3(TT)
    REAL,INTENT(IN) :: TT 
    AF3 = AF3TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Ice crystal function
  REAL FUNCTION BB3(TT)
    REAL,INTENT(IN) :: TT 
    BB3 = BB3TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Water droplet function:
  REAL FUNCTION BB3W(TT)
    REAL,INTENT(IN) :: TT 
    BB3W = BB3WTAB(NINT(XNDEGR*TT))
  END FUNCTION

! Function for IN concentration reduction between 0 and -25 C:
  REAL FUNCTION REDIN(TT)
    REAL,INTENT(IN) :: TT 
    REDIN = REDINTAB(NINT(XNDEGR*TT))
  END FUNCTION
END MODULE MODE_TIWMX
