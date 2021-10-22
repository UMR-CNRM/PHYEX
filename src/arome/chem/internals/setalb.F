!     ######spl
      SUBROUTINE setalb(nw,wl,albedo)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE
!_________________________________________________
!ks:   INCLUDE 'params' has been inserted here
!_________________________________________________
! BROADLY USED PARAMETERS:
!_________________________________________________
! i/o file unit numbers
      INTEGER kout, kin
! output
      PARAMETER(kout=6)
! input
      PARAMETER(kin=78)
!_________________________________________________
! altitude, wavelength grids
      INTEGER kz, kw
! altitude
      PARAMETER(kz=151)
! wavelength
      PARAMETER(kw=150)
!_________________________________________________
! number of weighting functions
      INTEGER ks, kj
!  wavelength dependent
      PARAMETER(ks=40)
!  wavelength and altitude dependent
      PARAMETER(kj=50)

! delta for adding points at beginning or end of data grids
      REAL deltax
      PARAMETER (deltax = 1.E-4)

! some constants...

! pi:
      REAL pi
      PARAMETER(pi=3.1415926535898)

! radius of the earth:
      REAL radius
      PARAMETER(radius=6.371E+3)

! largest number of the machine:
      REAL largest
      PARAMETER(largest=1.E+36)
!_________________________________________________

! input: (wavelength working grid data)
      INTEGER nw
      REAL wl(kw)

! output:
      REAL albedo(kw)

! local:
      INTEGER iw
      REAL alb
!_______________________________________________________________________

! set

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SETALB',0,ZHOOK_HANDLE)
      alb = 0.10
      WRITE(kout,*)'wavelength-independent albedo = ', alb
      DO 10, iw = 1, nw - 1

         albedo(iw) = alb

   10 CONTINUE
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('SETALB',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE SETALB
