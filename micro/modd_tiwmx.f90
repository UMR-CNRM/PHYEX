!@no_insert_drhook
!     ######spl
      MODULE MODD_TIWMX
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
INTEGER, PARAMETER :: NSTART = 13200 ! A too small value may result into a FPE in single precision mode. REK.
INTEGER, PARAMETER :: NSTOP = 37316

TYPE TIWMX_t
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
END TYPE TIWMX_t

TYPE(TIWMX_t), SAVE, TARGET :: TIWMX

END MODULE MODD_TIWMX
