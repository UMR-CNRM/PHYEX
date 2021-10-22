!     ######spl
      SUBROUTINE setozo(dobnew,                                         &
     &     nz,z,nw,wl,                                                  &
     &     xso3,s226,s263,s298,tlay,                                    &
     &     dto3)
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

!*******
! input:
!*******

! grids:
      REAL wl(kw)
      REAL z(kz)
      INTEGER nw
      INTEGER nz
! ozone absorption cross sections at 226, 263, and 298 K:
      REAL xso3(kw), s226(kw),s263(kw),s298(kw)
      REAL dobnew
! mid-layer temperature:
      REAL tlay(kz)

!*******
! output:
!*******
      REAL dto3(kz,kw)


!*******
! local:
!*******

      REAL cz(kz)

! ozone profile data:

      REAL zd(kdata), o3(kdata)
      REAL cd(kdata)
      REAL hscale
      REAL dobold, scale
      REAL colold, colnew
      REAL so3

! other:
      INTEGER i, iw, nd

!*******
! External functions:
!*******
      REAL fsum
      EXTERNAL fsum
!_______________________________________________________________________


! read in ozone profile

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SETOZO',0,ZHOOK_HANDLE)
      WRITE(kout,*) 'ozone profile: USSA, 1976'

      OPEN(kin,FILE='DATAE1/ATM/ussa.ozone',STATUS='old')
      DO i = 1, 7
        READ(kin,*)
      ENDDO
      nd = 1
 4    CONTINUE
         READ(kin,*,END=5) zd(nd), o3(nd)
         nd = nd+1
         GOTO 4
 5    CONTINUE
      CLOSE(kin)
      nd = nd-1

! compute column increments

      DO 11, i = 1, nd - 1
         cd(i) = (o3(i+1)+o3(i)) * 1.E5 * (zd(i+1)-zd(i)) / 2.
   11 CONTINUE

! Include exponential tail integral from infinity to 50 km,
! fold tail integral into top layer
! specify scale height near top of data.

      hscale = 4.50e5
      cd(nd-1) = cd(nd-1) + hscale * o3(nd)

! alternative input ozone concentration data could include, e.g.,
! a read file here:

!**********
!********** end data input.

! Compute column increments on standard z-grid.

      CALL inter3(nz,z,cz, nd,zd,cd, 1)

! scale values of cz(i) by any dobson unit

      colold = fsum(nz-1,cz)
      dobold = colold/2.687e16
      WRITE(kout,100) colold, dobold
  100 FORMAT(5x,'old O3 Column = ', 1pe11.4,1x,'# cm-2  = ',            &
     &     0pf8.2, '  Dobson Units ')

      if (dobnew .lt. 0.) then
         scale = 1.
      else
         scale = dobnew/dobold
      endif

      do i = 1, nz-1
         cz(i) = cz(i) * scale
      enddo
      colnew = fsum(nz-1,cz)
      WRITE(kout,105) colnew, colnew/2.687E16
  105 format(5x,'new O3 Column = ', 1pe11.4,1x,'# cm-2  = ',            &
     &     0pf8.2, '  Dobson Units ')

!***********************************
! calculate ozone optical depth for each layer, with temperature
! correction.  Output, dto3(kz,kw)

      DO 20, iw = 1, nw-1
         so3 = xso3(iw)
         DO 10, i = 1, nz - 1

            IF ( wl(iw) .GT. 240.5  .AND. wl(iw+1) .LT. 350. ) THEN
               IF (tlay(i) .LT. 263.) THEN
                  so3 = s226(iw) + (s263(iw)-s226(iw)) / (263.-226.) *  &
     &                 (tlay(i)-226.)
               ELSE
                  so3 = s263(iw) + (s298(iw)-s263(iw)) / (298.-263.) *  &
     &              (tlay(i)-263.)
               ENDIF
            ENDIF

            dto3(i,iw) = cz(i)*so3

   10    CONTINUE
   20 CONTINUE
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('SETOZO',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE setozo
