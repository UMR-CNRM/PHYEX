!     ######spl
      SUBROUTINE gridw(nw,wl,wc,wu)
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
      PARAMETER(kout=20)
! input
      PARAMETER(kin=85)
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

! output:

      REAL wl(kw), wc(kw), wu(kw)
      INTEGER nw

! local:

      REAL wincr
      INTEGER iw
      LOGICAL ok
      INTEGER idum
      REAL dum
      INTEGER mopt

!_______________________________________________________________________

!*** chose wavelengths

! some pre-set options
!     mopt = 1    equal spacing
!     mopt = 2    Isaksen's grid
!     mopt = 3    combined Kockarts/Isaksen grid
!     mopt = 4    user-defined

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('GRIDW',0,ZHOOK_HANDLE)
      mopt = 2
      IF (mopt .EQ. 1) GO TO 1
      IF (mopt .EQ. 2) GO TO 2
      IF (mopt .EQ. 3) GO TO 3
      IF (mopt .EQ. 4) GO TO 4

 1    CONTINUE
      nw = 140 + 1
      wincr = 1.0
      DO 10, iw = 1, nw-1
         wl(iw) = 280. + wincr*FLOAT(iw-1)
         wu(iw) = wl(iw) + wincr
         wc(iw) = ( wl(iw) + wu(iw) )/2.
   10 CONTINUE
      wl(nw) = wu(nw-1)
      GO TO 9

 2    CONTINUE
      nw = 0
      OPEN(unit=kin,file='DATAE1/GRIDS/isaksen.grid',status='old')
      DO iw = 1, 2
         READ(kin,*)
      ENDDO
      DO iw = 1, 130
         nw = nw + 1
         READ(kin,*) idum, dum, wc(nw), wl(nw), wu(nw)
      ENDDO
      CLOSE(kin)
      nw = nw + 1
      wl(nw) = wu(nw-1)
      GO TO 9

!** grid for strat photolysis calculations, extended at short wavelengths

 3    CONTINUE
      nw = 0
      OPEN(unit=kin,file='DATAE1/GRIDS/kockarts.grid',status='old')
      DO iw = 1, 16
         nw = nw + 1
         READ(kin,*) wl(nw), wu(nw)
         wc(iw) = ( wl(nw) + wu(nw) ) / 2.
      ENDDO
      CLOSE(kin)

      OPEN(unit=kin,file='DATAE1/GRIDS/isaksen.grid',status='old')
      DO iw = 1, 2 + 10
         READ(kin,*)
      ENDDO
      DO iw = 11, 130
         nw = nw + 1
         READ(kin,*) idum, dum, wc(nw), wl(nw), wu(nw)
      ENDDO
      CLOSE(kin)
      nw = nw + 1
      wl(nw) = wu(nw-1)
      GO TO 9

 4    CONTINUE
! define wavelength intervals of width 1 nm from 150 - 420 nm:
      nw = 1
      wl(1) = 150.
      DO iw = 151, 420
        wu(nw) = Float(iw)
        wc(nw) = (wl(nw) + wu(nw))/2.
        nw = nw+1
        wl(nw) = Float(iw)
      ENDDO
! define wavelength intervals of width 10 nm from 420 - 700 nm:
      DO iw = 430, 700, 10
        wu(nw) = Float(iw)
        wc(nw) = (wl(nw) + wu(nw))/2.
        nw = nw+1
        wl(nw) = Float(iw)
      ENDDO


 9    CONTINUE

!**
! write to record

      WRITE(kout,*)'w-grid:',nw,wl(1),wl(nw)

! check grid for assorted improprieties:

      CALL gridck(kw,nw,wl,ok)

      IF (.NOT. ok) THEN
         WRITE(kout,*)'STOP in GRIDW:  The w-grid does not make sense'
         STOP
      ENDIF

!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('GRIDW',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE GRIDW
