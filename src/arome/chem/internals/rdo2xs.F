!     ######spl
      SUBROUTINE rdo2xs(nw,wl,xso2)
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
! O2 absorption cross section. Units are cm2 molecule-1

      REAL xso2(kw)

! local:
      REAL x1(kdata), x
      REAL y1(kdata), y
      REAL yg(kw)
      INTEGER i, iw, n, ifirst
      INTEGER icount
      INTEGER ierr
      CHARACTER*40 fil

!_______________________________________________________________________

!************ O2 absorption cross sections:
! from 116 nm to 245 nm, including Schumann-Runge continumm
! from Brasseur and Solomon 1986.

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('RDO2XS',0,ZHOOK_HANDLE)
      fil =  'DATAE1/O2/O2_src.abs'
      OPEN(UNIT=kin,FILE='DATAE1/O2/O2_src.abs')
      READ(kin,*) ifirst, n
      DO i = 1, ifirst-2
        READ(kin,*)
      ENDDO
      icount = 0
      DO i = 1, n
        READ(kin,*) x, y
        IF (x .LT. 204.) THEN
          icount = icount+1
          x1(icount) = x
          y1(icount) = y
        ENDIF
      ENDDO
      CLOSE(kin)

! overwrite from 204 to 241 nm (Herzberg continuum)

      OPEN(UNIT=kin,FILE='DATAE1/O2/O2_jpl94.abs',STATUS='old')
      read(kin,*)
      n = 40
      DO i = 1, n
         icount = icount+1
         READ(kin,*) y
         y1(icount) = y*1E-24
         x1(icount) = 204. + FLOAT(i-1)
      END DO
      CLOSE (kin)

         CALL addpnt(x1,y1,kdata,icount,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,icount,               0.,0.)
         CALL addpnt(x1,y1,kdata,icount,x1(icount)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,icount,            1.e38,0.)

      CALL ch_inter2(nw,wl,yg,icount,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      DO iw = 1, nw-1
         xso2(iw) = yg(iw)
      END DO

      IF (LHOOK) CALL DR_HOOK('RDO2XS',1,ZHOOK_HANDLE)
      ENDSUBROUTINE RDO2XS
