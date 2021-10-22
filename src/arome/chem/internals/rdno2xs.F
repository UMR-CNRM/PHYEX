!     ######spl
      SUBROUTINE rdno2xs(nw,wl,xsno2)
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

      REAL xsno2(kw)

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
!     measurements of Davidson et al. (198x) at 273K
!     from 263.8 to 648.8 nm in approximately 0.5 nm intervals

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('RDNO2XS',0,ZHOOK_HANDLE)
      fil = 'DATAE1/NO2/NO2_ncar_00.abs'
      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      n = 750
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), dum, dum, idum
      ENDDO
      CLOSE(kin)

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
         xsno2(l) = yg(l)
   13 CONTINUE

!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('RDNO2XS',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE RDNO2XS
