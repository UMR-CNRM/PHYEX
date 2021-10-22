!     ######spl
      SUBROUTINE setno2(no2new,                                         &
     &     nz,z,nw,wl,                                                  &
     &     xsno2, tlay,                                                 &
     &     dtno2)
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
      PARAMETER(kdata=51)

!*******
! input:
!*******

! grids:

      REAL wl(kw)
      REAL z(kz)
      INTEGER nw
      INTEGER nz
      REAL no2new

! mid-layer temperature:

      REAL tlay(kz)

!*******
! output:
!*******

      REAL dtno2(kz,kw)

!*******
! local:
!*******

! absorption cross sections

      REAL xsno2(kw)
      REAL cz(kz)

! nitrogen dioxide profile data:

      REAL zd(kdata), no2(kdata)
      REAL cd(kdata)
      REAL hscale
      REAL colold, colnew
      REAL scale
      REAL sno2

! other:

      INTEGER i, l, nd

!*******
! External functions:
!*******
      REAL fsum
      EXTERNAL fsum

!_______________________________________________________________________
! Data input:

! Example:  set to 1 ppb in lowest 1 km, set to zero above that.
! - do by specifying concentration at 3 altitudes.

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SETNO2',0,ZHOOK_HANDLE)
      write(kout,*) 'NO2:  1 ppb in lowest 1 km, 0 above'

      nd = 3
      zd(1) = 0.
      no2(1) = 1. * 2.69e10

      zd(2) = 1.
      no2(2) = 1. * 2.69e10

      zd(3) = zd(2)* 1.000001
      no2(3) = 0.

!     zd(4) = zd(3)*1.1
!     no2(4) = 0.

! compute column increments (alternatively, can specify these directly)

      DO 11, i = 1, nd - 1
         cd(i) = (no2(i+1)+no2(i)) * 1.E5 * (zd(i+1)-zd(i)) / 2.
   11 CONTINUE

! Include exponential tail integral from top level to infinity.
! fold tail integral into top layer
! specify scale height near top of data (use ozone value)

      hscale = 4.50e5
      cd(nd-1) = cd(nd-1) + hscale * no2(nd)

!**********
!********** end data input.

! Compute column increments on standard z-grid.

      CALL inter3(nz,z,cz, nd,zd,cd, 1)

! scale values of cz(i)

      colold = fsum(nz-1,cz)
      WRITE(kout,100) colold, colold/2.687E16
  100 FORMAT(5x,'old NO2 Column = ', 1pe11.4,1x,'# cm-2  = ',           &
     &     0pf8.2, '  Dobson Units ')

      if ( (no2new .lt. 0.)  .or.  (colold .le. 0.) ) then
         scale = 1.
      else
         scale =  no2new/colold
      endif

      do i = 1, nz-1
         cz(i) = cz(i) * scale
      enddo
      colnew = fsum(nz-1,cz)
      WRITE(kout,105) colnew, colnew/2.687E16
  105 format(5x,'new NO2 Column = ', 1pe11.4,1x,'# cm-2  = ',           &
     &     0pf8.2, '  Dobson Units ')

!***********************************
! calculate optical depth for each layer, with temperature
! correction.  Output, dtno2(kz,kw)

      DO 20, l = 1, nw-1
         sno2 = xsno2(l)
         DO 10, i = 1, nz-1
            dtno2(i,l) = cz(i)*sno2
   10    CONTINUE
   20 CONTINUE
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('SETNO2',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE SETNO2
