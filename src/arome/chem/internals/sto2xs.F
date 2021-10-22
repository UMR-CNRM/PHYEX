!     ######spl
       SUBROUTINE sto2xs(nz,nw,xso2,nj,sq)
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

! calling parameters

      INTEGER nz, nw
      REAL xso2(kz,kw)
      INTEGER nj
      REAL sq(kj,kz,kw)

! local

      REAL qy
      INTEGER iw, iz
!______________________________________________________________________________

! quantum yield assumed to be unity
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('STO2XS',0,ZHOOK_HANDLE)
      qy = 1.

! assign cross section values at all wavelengths and at all altitudes
      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(nj,iz,iw) = qy * xso2(iz,iw)
        ENDDO
      ENDDO
!______________________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('STO2XS',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE STO2XS
