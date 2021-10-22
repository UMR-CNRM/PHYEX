!     ######spl
      SUBROUTINE tridag(a,b,c,r,u,n)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!_______________________________________________________________________
! solves tridiagonal system.  From Numerical Recipies, p. 40
!_______________________________________________________________________

      IMPLICIT NONE

! input:
      INTEGER n
      REAL a, b, c, r
      DIMENSION a(n),b(n),c(n),r(n)

! output:
      REAL u
      DIMENSION u(n)

! local:
      INTEGER j

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
      REAL bet, gam
      DIMENSION gam(2*kz)
!_______________________________________________________________________

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('TRIDAG',0,ZHOOK_HANDLE)
      IF (b(1) .EQ. 0.) STOP 1001
      bet   = b(1)
      u(1) = r(1)/bet
      DO 11, j = 2, n
         gam(j) = c(j - 1)/bet
         bet = b(j) - a(j)*gam(j)
         IF (bet .EQ. 0.) STOP 2002
         u(j) = (r(j) - a(j)*u(j - 1))/bet
   11 CONTINUE
      DO 12, j = n - 1, 1, -1
         u(j) = u(j) - gam(j + 1)*u(j + 1)
   12 CONTINUE
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('TRIDAG',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE tridag
