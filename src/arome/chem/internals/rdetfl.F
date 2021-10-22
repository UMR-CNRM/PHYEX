!     ######spl
      SUBROUTINE rdetfl(nw,wl,f)
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

      integer kdata
      parameter(kdata=20000)

! input: (wavelength grid)
      INTEGER nw
      REAL wl(kw)
      INTEGER iw

! output: (extra terrestrial solar flux)
      REAL f(kw)

! INTERNAL:

! work arrays for input data files:

      CHARACTER*40 fil
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)
      INTEGER nhead, n, i, ierr

! data gridded onto wl(kw) grid:

      REAL yg1(kw)
      REAL yg2(kw)

      REAL hc
      PARAMETER(hc = 6.62E-34 * 2.998E8)

!_______________________________________________________________________
! simple files are read and interpolated here in-line. Reading of
! more complex files may be done with longer code in a read#.f subroutine.

!****** SUSIM irradiance

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('RDETFL',0,ZHOOK_HANDLE)
      fil = 'DATAE1/SUN/susim_hr.flx'
      CALL read1(nw,wl,yg1)

!      fil = 'DATAE1/SUN/wmo85.flx'
!      CALL read2(nw,wl,yg1)

!      fil = 'DATAE1/SUN/extsol.flx'
!      nhead = 3
!      n =121

!      fil = 'DATAE1/SUN/lowsun3.flx'
!      nhead = 3
!      n = 4327

!      fil = 'DATAE1/SUN/modtran1.flx'
!      nhead = 6
!      n = 14980

!      fil = 'DATAE1/SUN/nicolarv.flx'
!      nhead = 8
!      n = 1260

!      fil = 'DATAE2/SUN/solstice.flx'
! unofficial - do not use
!      nhead = 11
!      n = 2047

!      fil = 'DATAE2/SUN/suntoms.flx'
! unofficial - do not use
!      nhead = 3
!      n = 1200
!* suntoms.flx: original data is in milli-Watts m-2 nm-1
!*         y1(i) = y1(i) / 1000.

      fil = 'DATAE1/SUN/neckel.flx'
      nhead = 11
      n = 496
! neckel.flx:  convert from quanta cm-2 s-1 to watts m-2 nm-1
!         y1(i) = y1(i) * 1.E4 * hc / (x1(i) * 1.E-9)

      write(kout,*) fil
      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      DO i = 1, nhead
         READ(kin,*)
      ENDDO
      DO 12, i = 1, n

! do file-specific read, and apply conversions here

         READ(kin,*) x1(i), y1(i)

! neckel.flx:  convert from quanta cm-2 s-1 to watts m-2 nm-1
         y1(i) = y1(i) * 1.E4 * hc / (x1(i) * 1.E-9)

   12 CONTINUE
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
      CALL ch_inter2(nw,wl,yg2,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

! Adjoin susim and neckel

      DO 10,  iw = 1, nw-1
         IF (wl(iw) .GT. 350.) THEN
            f(iw) = yg2(iw)
         ELSE
           f(iw) = yg1(iw)
         ENDIF
   10 CONTINUE


!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('RDETFL',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE RDETFL
