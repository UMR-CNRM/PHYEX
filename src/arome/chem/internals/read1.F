!     ######spl
      SUBROUTINE read1(nw,wl,f)
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

! input: (wavelength grid)
      INTEGER nw
      REAL wl(kw)
      INTEGER iw

! output: (extra terrestrial solar flux)
      REAL f(kw)

! local:

      REAL lambda_hi(10000),irrad_hi(10000)
      REAL lambda
      INTEGER ierr
      INTEGER i, j, n
      CHARACTER*40 FIL

!_______________________________________________________________________

!****** SUSIM irradiance
!_______________________________________________________________________
! VanHoosier, M. E., J.-D. F. Bartoe, G. E. Brueckner, and
! D. K. Prinz, Absolute solar spectral irradiance 120 nm -
! 400 nm (Results from the Solar Ultraviolet Spectral Irradiance
! Monitor - SUSIM- Experiment on board Spacelab 2),
! Astro. Lett. and Communications, 1988, vol. 27, pp. 163-168.
!     SUSIM SL2 high resolution (0.15nm) Solar Irridance data.
!     Irradiance values are given in milliwatts/m^2/nanomenters
!     and are listed at 0.05nm intervals.  The wavelength given is
!     the center wavelength of the 0.15nm triangular bandpass.
!     Normalized to 1 astronomical unit.
!  DATA for wavelengths > 350 nm are unreliable
! (Van Hoosier, personal communication, 1994).
!_______________________________________________________________________

!* high resolution

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('READ1',0,ZHOOK_HANDLE)
      fil = 'DATAE1/SUN/susim_hi.flx'
      write(kout,*) fil
      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      DO 11, i = 1, 7
         READ(kin,*)
   11 CONTINUE
      DO 12, i = 1, 559
         READ(kin,*)lambda,(irrad_hi(10*(i-1)+j), j=1, 10)
   12 CONTINUE
      CLOSE (kin)

! compute wavelengths, convert from mW to W

      n = 559*10
      DO 13, i = 1, n
         lambda_hi(i)=120.5 + FLOAT(i-1)*.05
         irrad_hi(i) = irrad_hi(i)  /  1000.
   13 CONTINUE
!_______________________________________________________________________

      CALL addpnt(lambda_hi,irrad_hi,10000,n,                           &
     &            lambda_hi(1)*(1.-deltax),0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,                 0.,0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,                           &
     &            lambda_hi(n)*(1.+deltax),0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,              1.e38,0.)
      CALL ch_inter2(nw,wl,f,n,lambda_hi,irrad_hi,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('READ1',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE READ1
