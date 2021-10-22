!     ######spl
      SUBROUTINE read2(nw,wl,f)
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
      REAL yg(kw)

!
      INTEGER iw

! output: (extra terrestrial solar flux)
      REAL f(kw)

! local:

      REAL x1(1000), y1(1000)
      REAL x2(1000), y2(1000)
      REAL x3(1000), y3(1000)
      INTEGER i, n
      REAL DUM
      INTEGER IDUM

!_______________________________________________________________________

!********WMO 85 irradiance

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('READ2',0,ZHOOK_HANDLE)
      WRITE(kout,*) 'DATAE1/SUN/wmo85.flx'
      OPEN(UNIT=kin,FILE='DATAE1/SUN/wmo85.flx',STATUS='old')
      DO 11, i = 1, 3
         READ(kin,*)
   11 CONTINUE
      n = 158
      DO 12, i = 1, n
         READ(kin,*) idum, x1(i),x2(i),y1(i), dum, dum, dum
         x3(i) = 0.5 * (x1(i) + x2(i))

! average value needs to be calculated only if inter2 is
! used to interpolate onto wavelength grid (see below)
!        y1(i) =  y1(i) / (x2(i) - x1(i))

   12 CONTINUE
      CLOSE (kin)

      x1(n+1) = x2(n)

! inter2: INPUT : average value in each bin
!         OUTPUT: average value in each bin
! inter3: INPUT : total area in each bin
!         OUTPUT: total area in each bin

      CALL inter3(nw,wl,yg, n+1,x1,y1,0)
!      CALL ch_inter2(nw,wl,yg,n,x3,y1,ierr)

      DO 10,  iw = 1, nw-1
! from quanta s-1 cm-2 bin-1 to  watts m-2 nm-1
! 1.e4 * ([hc =] 6.62E-34 * 2.998E8)/(wc*1e-9)

! the scaling by bin width needs to be done only if
! inter3 is used for interpolation

         yg(iw) = yg(iw) / (wl(iw+1)-wl(iw))
         f(iw) = yg(iw) * 1.e4 * (6.62E-34 * 2.998E8) /                 &
     &        ( 0.5 * (wl(iw+1)+wl(iw)) * 1.e-9)

   10 CONTINUE

!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('READ2',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE READ2
