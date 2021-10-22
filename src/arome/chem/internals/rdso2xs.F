!     ######spl
      SUBROUTINE rdso2xs(nw,wl,xsso2)
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
      PARAMETER(kdata=1000)

! input: (altitude working grid)
      INTEGER nw
      REAL wl(kw)

! output:

      REAL xsso2(kw)

! local:
      REAL x1(kdata)
      REAL y1(kdata),y2(kdata),y3(kdata)
      REAL yg(kw)
      REAL a1, a2, dum
      INTEGER ierr
      INTEGER i, l, n, idum
      CHARACTER*40 fil
!_______________________________________________________________________

!************ absorption cross sections:
! SO2 absorption cross sections from J. Quant. Spectrosc. Radiat. Transfer
! 37, 165-182, 1987, T. J. McGee and J. Burris Jr.
! Angstrom vs. cm2/molecule, value at 221 K

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('RDSO2XS',0,ZHOOK_HANDLE)
      fil = 'DATA/McGee87'
      OPEN(UNIT=kin,FILE='DATAE1/SO2/SO2xs.all',STATUS='old')
      DO 11, i = 1,3
         read(kin,*)
   11 CONTINUE
!      n = 681
      n = 704
      DO 12, i = 1, n
         READ(kin,*) x1(i), y1(i)
         x1(i) = x1(i)/10.
   12 CONTINUE
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      DO 13, l = 1, nw-1
         xsso2(l) = yg(l)
   13 CONTINUE

!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('RDSO2XS',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE RDSO2XS
