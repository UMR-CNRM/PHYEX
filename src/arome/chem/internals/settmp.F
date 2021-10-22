!     ######spl
      SUBROUTINE settmp(nz,z,tlev,tlay)
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

      INTEGER kdata
      PARAMETER(kdata=150)

! input: (altitude working grid)
      REAL z(kz)
      INTEGER nz

! output:
      REAL tlev(kz), tlay(kz)

! local:
      REAL zd(kdata), td(kdata)
      INTEGER i, nd
!_______________________________________________________________________


! read in temperature profile

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SETTMP',0,ZHOOK_HANDLE)
      WRITE(kout,*) 'air temperature: USSA, 1976'

      OPEN(kin,FILE='DATAE1/ATM/ussa.temp',STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      nd = 1
 4    CONTINUE
         READ(kin,*,END=5) zd(nd), td(nd)
         nd = nd+1
         GOTO 4
 5    CONTINUE
      CLOSE(kin)
      nd = nd-1

! use constant temperature to infinity:

      zd(nd) = 1.E10

! alternative input temperature data could include, e.g., a read file here:

!**********
!********** end data input.

! interpolate onto z-grid

      CALL ch_inter1(nz,z,tlev,nd,zd,td)

! compute layer-averages

      DO 20, i = 1, nz - 1
         tlay(i) = (tlev(i+1) + tlev(i))/2.
 20   CONTINUE
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('SETTMP',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE SETTMP
