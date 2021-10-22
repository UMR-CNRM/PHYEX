!     ######spl
      SUBROUTINE gridz(nz,z,az,nlevel)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!XXX      SUBROUTINE gridz(nz,z)
!KS:8==========================================================================

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

!KS:9==========================================================================
! input: altitudes of calling model (CAUTION: AZ is given in meters,
!                                             Z will be in kilometers!!)
      REAL az(*)
      INTEGER nlevel
!KS:9==========================================================================

! output: altitude working grid:

      REAL z(kz)
      INTEGER nz

! local:

      REAL zincr
      INTEGER i
      LOGICAL ok
!_______________________________________________________________________

! set vertical grid of the atmosphere.  All values should be in km.
! User specifies upright grid (surface at lowest km value, increasing
! upwards:
!     -  NZ = total number of user levels
!     -  Z(I) = altitude in km for each level.
! Note "levels" are vertical points
!      "layers" are vertical distances between levels

! set atmospheric level altitudes (in real km), including
! top-most level.
! non-uniform spacing is possible

!KS:10=========================================================================
!XXX      nz = 51
!XXX      z(1) = 0.
!XXX      zincr = 1.
!XXX      DO 10, i = 2, nz
!XXX         z(i) = z(1) + FLOAT(i-1) * zincr
!XXX    10 CONTINUE
!
!      use model levels for vertical grid where available
       REAL(KIND=JPRB) :: ZHOOK_HANDLE
       IF (LHOOK) CALL DR_HOOK('GRIDZ',0,ZHOOK_HANDLE)
       do 10, i = 1, nlevel
         z(i) = az(i) *1E-3
10     continue
!
       nz = nlevel
!
!      fill up between model top and 50km with 1km grid spacing
20     continue
       if (z(nz) .ge. 50.) goto 30
         nz = nz + 1
         if (nz .gt. kz) stop "GRIDZ: not enough memory, increase kz"
         z(nz) = z(nz-1) + 1.
       goto 20
!
30     continue
!
!KS:10=========================================================================

! write to record:

      WRITE(kout,*)'z-grid:',nz,z(1),z(nz)

! check grid for assorted improprieties:

      CALL gridck(kz,nz,z,ok)

      IF (.NOT. ok) THEN
         WRITE(kout,*)'STOP in GRIDZ:  The z-grid does not make sense'
         STOP
      ENDIF
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('GRIDZ',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE GRIDZ
