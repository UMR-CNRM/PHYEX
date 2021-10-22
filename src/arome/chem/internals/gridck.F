!     ######spl
      SUBROUTINE gridck(k,n,x,ok)
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

! input:
      INTEGER k, n
      REAL x(k)

! output:
      LOGICAL ok

! local:
      INTEGER i
!_______________________________________________________________________

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('GRIDCK',0,ZHOOK_HANDLE)
      ok = .TRUE.

! check if dimension meaningful and within bounds

      IF (n .GT. k) THEN
         ok = .false.
         WRITE(kout,100)
         IF (LHOOK) CALL DR_HOOK('GRIDCK',1,ZHOOK_HANDLE)
         RETURN
      ENDIF
  100 FORMAT('Number of data exceeds dimension')

      IF (n .LT. 2) THEN
         ok = .FALSE.
         WRITE(kout,101)
         IF (LHOOK) CALL DR_HOOK('GRIDCK',1,ZHOOK_HANDLE)
         RETURN
      ENDIF
  101 FORMAT('Too few data, number of data points must be >= 2')

! disallow negative grid values

      IF(x(1) .LT. 0.) THEN
         ok = .FALSE.
         WRITE(kout,105)
         IF (LHOOK) CALL DR_HOOK('GRIDCK',1,ZHOOK_HANDLE)
         RETURN
      ENDIF
  105 FORMAT('Grid cannot start below zero')

! check sorting

      DO 10, i = 2, n
         IF( x(i) .LE. x(i-1)) THEN
            ok = .FALSE.
            WRITE(kout,110)
            IF (LHOOK) CALL DR_HOOK('GRIDCK',1,ZHOOK_HANDLE)
            RETURN
         ENDIF
   10 CONTINUE
  110 FORMAT('Grid is not sorted or contains multiple values')
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('GRIDCK',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE GRIDCK
