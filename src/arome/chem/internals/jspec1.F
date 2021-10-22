!     ######spl
      SUBROUTINE jspec1(nw,wl,nz,tlev,airlev,                           &
     &     j,sq,jlabel)
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
      PARAMETER(kdata=16000)

! input

      INTEGER nw
      REAL wl(kw), wc(kw)

      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

! weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

! input/output:
      INTEGER j

! data arrays

      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

! local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL xso3(kw),s226(kw),s263(kw),s298(kw)
      REAL xsno2(kz,kw)
      REAL qy, qy1d, qy3p
      REAL p298, p230
      REAL phi1, phi2, phi20
      REAL tau, tau2, tau3
      REAL ak300, akt, temp
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL xl, xl0
      REAL so3
      REAL sigma
      REAL xs
      REAL t
      REAL dum
      INTEGER iyld
      INTEGER kmich, kjpl87, kjpl92, kshet
      INTEGER irow, icol
      INTEGER i, irev, iw, n, idum
      INTEGER ierr
      INTEGER iz, k
      REAL slope,lambda, sum
      REAL sumA, sumB, chi
      CHARACTER*120 inline

      REAL coeff(4,3), TBar, LBar

!_______________________________________________________________________

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('JSPEC1',0,ZHOOK_HANDLE)
      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

!***************************************************************
!*************** NO2 photodissociation

      j = j + 1
      jlabel(j) = 'NO2 -> NO + O(3P)'

! cross section
!     measurements of Davidson et al. (198x) at 273K
!     from 263.8 to 648.8 nm in approximately 0.5 nm intervals

!     OPEN(UNIT=kin,FILE='DATAE1/NO2/NO2_ncar_00.abs',STATUS='old')
!     n = 750
!     DO i = 1, n
!        READ(kin,*) x1(i), y1(i), dum, dum, idum
!     ENDDO
!     CLOSE(kin)

!     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
!     CALL addpnt(x1,y1,kdata,n,               0.,0.)
!     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
!     CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
!     CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
!     IF (ierr .NE. 0) THEN
!        WRITE(*,*) ierr, jlabel(j)
!        STOP
!     ENDIF

! cross section data from JPL 94 recommendation

      OPEN(UNIT=kin,FILE='DATAE1/NO2/NO2_jpl94.abs',STATUS='old')
      READ(kin,*) idum, n
      DO i = 1, idum-2
         READ(kin,*)
      ENDDO
! read in wavelength bins, cross section at T0 and temperature correction
! coefficient a;  see input file for details.
! data need to be scaled to total area per bin so that they can be used with
! inter3
      DO i = 1, n
         READ(kin,*) x1(i), x3(i), y1(i), dum, y2(i)
         y1(i) = (x3(i)-x1(i)) * y1(i)*1E-20
         y2(i) = (x3(i)-x1(i)) * y2(i)*1E-22
         x2(i) = x1(i)
      ENDDO
      CLOSE(kin)

      x1(n+1) = x3(n)
      x2(n+1) = x3(n)
      n = n+1
      n1 = n

      CALL inter3(nw,wl,yg1,n,x1,y1,0)

      CALL inter3(nw,wl,yg2,n1,x2,y2,0)

! yg1, yg2 are per nm, so rescale by bin widths
      DO iw = 1, nw-1
        yg1(iw) = yg1(iw)/(wl(iw+1)-wl(iw))
        yg2(iw) = yg2(iw)/(wl(iw+1)-wl(iw))
      ENDDO

      DO iw = 1, nw-1
        DO i = 1, nz
           xsno2(i,iw) = yg1(iw) + yg2(iw)*(tlev(i)-273.15)
        ENDDO
      ENDDO

! quantum yield
! from Gardiner, Sperry, and Calvert

      OPEN(UNIT=kin,FILE='DATAJ1/YLD/NO2_calvert.yld',STATUS='old')
      DO i = 1, 8
         READ(kin,*)
      ENDDO
      n = 66
      DO i = 1, n
         READ(kin,*) x1(i),y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,   0.)
      CALL ch_inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

! combine

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = xsno2(i,iw)*yg1(iw)
         ENDDO
      ENDDO


!************       jlabel(j) = 'O3 -> O2 + O(1D)'
!************       jlabel(j) = 'O3 -> O2 + O(3P)'

      j = j + 1
      jlabel(j) = 'O3 -> O2 + O(1D)'

      j = j + 1
      jlabel(j) = 'O3 -> O2 + O(3P)'

! cross sections:
! from WMO 1985 Ozone Assessment
! from 175.439 to 847.500 nm
! use value at 273 K

      OPEN(UNIT=kin,FILE='DATAE1/wmo85',STATUS='old')
      DO i = 1, 3
         read(kin,*)
      ENDDO
      n = 158
      DO i = 1, n
         READ(kin,*) idum, a1, a2, dum, dum, dum, dum, y1(i)
         x1(i) = (a1+a2)/2.
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw-1
         xso3(iw) = yg(iw)
      ENDDO

! For Hartley and Huggins bands, use temperature-dependent values from
! Molina, L. T., and M. J. Molina, Absolute absorption cross sections
! of ozone in the 185- to 350-nm wavelength range,
! J. Geophys. Res., vol. 91, 14501-14508, 1986.

      OPEN(UNIT=kin,FILE='DATAE1/O3/O3.molina.abs',STATUS='old')
      DO i = 1, 5
         READ(kin,*)
      ENDDO
      n1 = 220
      n2 = 220
      n3 = 220
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i), y2(i), y3(i)
         x2(i) = x1(i)
         x3(i) = x1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,            1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s226(iw) = yg(iw)*1.E-20
      ENDDO

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,            1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s263(iw) = yg(iw)*1.E-20
      ENDDO

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,               0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,            1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s298(iw) = yg(iw)*1.E-20
      ENDDO

! quantum yield:

      kmich = 1
      kjpl87 = 2
      kjpl92 = 3
      kshet = 4

! choose quantum yield recommendation:
!    kjpl87:  JPL recommendation 1987                - JPL 87, 90, 92 do not
!    kjpl92:  JPL recommendation 1990/92 (identical) - include "tail" region
!    kmich :  Michelsen et al., 1994
!    kshet :  Shetter et al., 1996

      iyld = kjpl87
      iyld = kjpl92
      iyld = kmich
      iyld = kshet

! read parameters from Michelsen, H. A., R.J. Salawitch, P. O. Wennber,
! and J. G. Anderson, Geophys. Res. Lett., 21, 2227-2230, 1994.

      IF (iyld .EQ. kmich) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/YLD/O3.param.yld',STATUS='old')
        READ(kin,*)
        READ(kin,*)
        READ(kin,*)
        n1 = 21
        n2 = n1
        DO i = 1, n1
           READ(kin,*) x1(i), y1(i), y2(i)
           x2(i) = x1(i)
        ENDDO
        CLOSE(kin)

        CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
        CALL addpnt(x1,y1,kdata,n1,               0.,y1(1))
        CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),y1(n1))
        CALL addpnt(x1,y1,kdata,n1,            1.e+38,y1(n1))
        CALL ch_inter2(nw,wl,yg1,n1,x1,y1,ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF

        CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
        CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
        CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
        CALL addpnt(x2,y2,kdata,n2,            1.e+38,y2(n2))
        CALL ch_inter2(nw,wl,yg2,n2,x2,y2,ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF
      ENDIF

! quantum yield data from
! Shetter et al, J.Geophys.Res., v 101 (D9), pg. 14,631-14,641, June 20, 1996

      IF (iyld .EQ. kshet) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/YLD/O3_shetter.yld',STATUS='OLD')
        READ(kin,*) idum, n
        DO i = 1, idum-2
          READ(kin,*)
        ENDDO
        n = n-2
        DO i = 1, n
          READ(kin,*) x1(i),y3(i),y4(i),y1(i),y2(i)
          x2(i) = x1(i)
          x3(i) = x1(i)
          x4(i) = x1(i)
        ENDDO
        DO i = n+1, n+2
           READ(kin,*) x3(i),y3(i),y4(i)
           x4(i) = x3(i)
        ENDDO
        CLOSE(kin)

        n1 = n
        n2 = n
        n3 = n+2
        n4 = n+2

! coefficients for exponential fit:

        CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax), y1(1))
        CALL addpnt(x1,y1,kdata,n1,                0., y1(1))
        CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
        CALL addpnt(x1,y1,kdata,n1,              1E38,0.)

        CALL ch_inter2(nw,wl,yg1, n1,x1,y1, ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF

        CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
        CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
        CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
        CALL addpnt(x2,y2,kdata,n2,              1E38,0.)

        CALL ch_inter2(nw,wl,yg2, n2,x2,y2, ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF

! phi data at 298 and 230 K

        CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),y3(1))
        CALL addpnt(x3,y3,kdata,n3,               0.,y3(1))
        CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
        CALL addpnt(x3,y3,kdata,n3,              1E38,0.)

        CALL ch_inter2(nw,wl,yg3, n3,x3,y3, ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr,jlabel(j)
           STOP
        ENDIF

        CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),y4(1))
        CALL addpnt(x4,y4,kdata,n4,               0.,y4(1))
        CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
        CALL addpnt(x4,y4,kdata,n4,              1E38,0.)

        CALL ch_inter2(nw,wl,yg4, n4,x4,y4, ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr,jlabel(j)
           STOP
        ENDIF
      ENDIF

! compute cross sections and yields at different wavelengths, altitudes:

      DO iw = 1, nw-1

         so3 = xso3(iw)

         DO i = 1, nz

            IF ( wl(iw) .GT. 240.5  .AND. wl(iw+1) .LT. 350. ) THEN
               IF (tlev(i) .LT. 263.) THEN
                  so3 = s226(iw) + (s263(iw)-s226(iw)) / (263.-226.) *  &
     &                 (tlev(i)-226.)
               ELSE
                  so3 = s263(iw) + (s298(iw)-s263(iw)) / (298.-263.) *  &
     &              (tlev(i)-263.)
               ENDIF
            ENDIF

! quantum yields
! coefficients from jpl 87:

             IF (iyld .EQ. kjpl87) THEN
               tau = tlev(i) - 230.
               tau2 = tau*tau
               tau3 = tau2*tau
               xl = wc(iw)
               xl0 = 308.2 + 4.4871e-2*tau + 6.938e-5*tau2 -            &
     &               2.5452e-6*tau3
               a = 0.9*(0.369 + 2.85e-4*tau + 1.28e-5*tau2 +            &
     &                  2.57e-8*tau3)
               b     = -0.575 + 5.59e-3*tau - 1.439e-5*tau2 -           &
     &                  3.27e-8*tau3
               c = 0.9*(0.518 + 9.87e-4*tau - 3.94e-5*tau2 +            &
     &                  3.91e-7*tau3)
               qy1d = a*atan(b*(xl-xl0)) + c
               qy1d = max(0.,qy1d)
               qy1d = min(0.9,qy1d)
             ENDIF

! from jpl90, jpl92:
! (caution: error in JPL92 for first term of a3)

             IF (iyld .EQ. kjpl92) THEN
               tau = 298. - tlev(i)
               tau2 = tau*tau
               xl0 = wc(iw) - 305.
               a0 =   .94932   - 1.7039e-4*tau + 1.4072E-6*tau2
               a1 = -2.4052e-2 + 1.0479e-3*tau - 1.0655e-5*tau2
               a2 =  1.8771e-2 - 3.6401e-4*tau - 1.8587e-5*tau2
               a3 = -1.4540e-2 - 4.7787e-5*tau + 8.1277e-6*tau2
               a4 =  2.3287e-3 + 1.9891e-5*tau - 1.1801e-6*tau2
               a5 = -1.4471e-4 - 1.7188e-6*tau + 7.2661e-8*tau2
               a6 =  3.1830e-6 + 4.6209e-8*tau - 1.6266e-9*tau2
               qy1d = a0 + a1*xl0 + a2*(xl0)**2 + a3*(xl0)**3 +         &
     &                a4*(xl0)**4 + a5*(xl0)**5 + a6*(xl0)**6
               IF (wc(iw) .LT. 305.) qy1d = 0.95
               IF (wc(iw) .GT. 320.) qy1d = 0.
               IF (qy1d .LT. 0.02) qy1d = 0.
             ENDIF

! from Michelsen, H. A., R.J. Salawitch, P. O. Wennber, and J. G. Anderson
! Geophys. Res. Lett., 21, 2227-2230, 1994.

           IF (iyld .EQ. kmich) THEN
             IF (wc(iw) .LT. 271.) THEN
                qy1d = 0.87
             ELSE IF (wc(iw) .GE. 271. .AND. wc(iw) .LT. 305.) THEN
                qy1d = 1.98 - 301./wc(iw)
             ELSE IF (wc(iw) .GE. 305. .AND. wc(iw) .LE. 325.) THEN
                qy1d = yg1(iw) * EXP (-yg2(iw) /(0.6951*tlev(i)))
             ELSE
                qy1d = 0.
             ENDIF
           ENDIF

! Shetter et al.:
! phi = A * exp(-B/T), A and B are based on meas. at 298 and 230 K
! do linear interpolation between phi(298) and phi(230) for wavelengths > 321
! as phi(230)=0. for those wavelengths, so there are no A and B factors

           IF (iyld .EQ. kshet) THEN
             IF (wl(iw+1) .LE. 321.) THEN
               qy1d = yg1(iw) * EXP(-1. * yg2(iw)/tlev(i))
             ELSE
               qy1d = (yg3(iw) - yg4(iw))/(298.-230.) * (tlev(i)-230.) +&
     &                 yg4(iw)
             ENDIF
           ENDIF

           sq(j-1,i,iw) = qy1d*so3
           qy3p = 1.0 - qy1d
           sq(j,i,iw) = qy3p*so3

         ENDDO
      ENDDO

!***************      jlabel(j) = 'NO3 -> NO2 + O(3P)'
!***************      jlabel(j) = 'NO3 -> NO + O2'

! cross section
!     measurements of Graham and Johnston 1978

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/NO3_gj78.abs',STATUS='old')
      DO i = 1, 9
         READ(kin,*)
      ENDDO
      n = 305
      DO irow = 1, 30
         READ(kin,*) ( y1(10*(irow-1) + icol), icol =  1, 10 )
      ENDDO
      READ(kin,*) ( y1(300 + icol), icol = 1, 5 )
      CLOSE (kin)
      DO i = 1, n
         y1(i) =  y1(i) * 1.E-19
         x1(i) = 400. + 1.*FLOAT(i-1)
      ENDDO

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!     cross section from JPL94:

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/NO3_jpl94.abs',STATUS='old')
      READ(kin,*) idum, n
      DO i = 1, idum-2
         READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i)*1E-20
      ENDDO
      CLOSE (kin)
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

! use JPL94 for wavelengths longer than 600 nm

      DO iw = 1, nw-1
         IF(wl(iw) .GT. 600.) yg(iw) = yg1(iw)
      ENDDO

! quantum yield:
! from Madronich (1988) see CEC NO3 book.

! for   NO3 ->NO+O2

      j = j + 1
      jlabel(j) = 'NO3 -> NO + O2'
      DO iw = 1, nw - 1
         IF (wc(iw).LT.584.) THEN
            qy = 0.
         ELSEIF (wc(iw).GE.640.) THEN
            qy = 0.
         ELSEIF (wc(iw).GE.595.) THEN
            qy = 0.35*(1.-(wc(iw)-595.)/45.)
         ELSE
            qy = 0.35*(wc(iw)-584.)/11.
         ENDIF
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

! for  NO3 ->NO2+O
      j = j + 1
      jlabel(j) = 'NO3 -> NO2 + O(3P)'
      DO iw = 1, nw - 1
         IF (wc(iw).LT.584.) THEN
            qy = 1.
         ELSEIF (wc(iw).GT.640.) THEN
            qy = 0.
         ELSEIF (wc(iw).GT.595.) THEN
            qy = 0.65*(1-(wc(iw)-595.)/45.)
         ELSE
            qy = 1.-0.35*(wc(iw)-584.)/11.
         ENDIF
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

!*************** N2O5 photodissociation

      j = j + 1
      jlabel(j) = 'N2O5 -> NO3 + NO + O(3P)'

      j = j + 1
      jlabel(j) = 'N2O5 -> NO3 + NO2'

! cross section from jpl94, table up to 280 nm

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/N2O5.abs',STATUS='old')
      DO i = 1, 21
         READ(kin,*)
      ENDDO
      n = 17
      DO i = 1,  n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE(kin)

! fill in from 281 nm trough 400 nm with value at 262.5 K:
! this will be overwritten in most cases.

      DO i = 1, 120
         n = n + 1
         x1(n) = 281. + 1.*FLOAT(i-1)
         y1(n) = 1.E-20*EXP( 2.735 + (4728.5-17.127*x1(n)) / 262.5 )
      ENDDO

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,0.               ,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

! quantum yield : see DATAJ1/YLD/N2O5.qy for explanation
! correct for T-dependence of cross section

      DO iw = 1, nw - 1

         qy = MIN( 1., 3.832441 - 0.012809638 * wc(iw) )
         qy = MAX( 0., qy )

         DO i = 1, nz
            IF (wc(iw) .GE. 285. .AND. wc(iw) .LE. 380.) THEN
              t = MAX(225.,MIN(tlev(i),300.))
              sq(j-1,i,iw) = qy *                                       &
     &            1.E-20*EXP( 2.735 + (4728.5-17.127*wc(iw)) / t )
               sq(j,i,iw) = (1.-qy) *                                   &
     &            1.E-20*EXP( 2.735 + (4728.5-17.127*wc(iw)) / t )
            ELSE
               sq(j-1,i,iw) = qy * yg(iw)
               sq(j,i,iw) = (1.-qy) * yg(iw)
            ENDIF
         ENDDO
      ENDDO

!*************** HNO2 photodissociation

! cross section from JPL92
! (from Bongartz et al., identical to JPL94 recommendation)

      j = j + 1
      jlabel(j) = 'HNO2 -> OH + NO'
!      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO2_jpl85.abs',STATUS='old')
      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO2_jpl92.abs',STATUS='old')
      DO i = 1, 13
         READ(kin,*)
      ENDDO
!      n = 86
      n = 91
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

! quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

!*************** HNO3 photodissociation

       j = j + 1
       jlabel(j) = 'HNO3 -> OH + NO2'

!* cross section from JPL85
!
!      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO3.abs',STATUS='old')
!      DO i = 1, 9
!         READ(kin,*)
!      ENDDO
!      n = 29
!      DO i = 1, n
!         READ(kin,*) x1(i), y1(i)
!         y1(i) = y1(i) * 1.E-20
!      ENDDO
!      CLOSE (kin)
!
!      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,               0.,0.)
!      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
!      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
!      IF (ierr .NE. 0) THEN
!         WRITE(*,*) ierr, jlabel(j)
!         STOP
!      ENDIF
!
!* quantum yield = 1
!
!      qy = 1.
!      DO iw = 1, nw - 1
!         DO i = 1, nz
!            sq(j,i,iw) = yg(iw)*qy
!         ENDDO
!      ENDDO


! HNO3 cross section parameters from Burkholder et al. 1993

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO3_burk.abs',STATUS='old')
      DO i = 1, 6
         READ(kin,*)
      END DO
      n1 =  83
      n2 = n1
      DO i = 1, n1
         READ(kin,*) y1(i), y2(i)
         x1(i) = 184. + i*2.
         x2(i) = x1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,            1.e+38,0.)
      CALL ch_inter2(nw,wl,yg1,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF


      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
      CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
      CALL addpnt(x2,y2,kdata,n2,            1.e+38,y2(n2))
      CALL ch_inter2(nw,wl,yg2,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

! quantum yield = 1
! correct for temperature dependence

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg1(iw) * 1.E-20                               &
     &           * exp( yg2(iw)/1.e3*(tlev(i)-298.) )
         END DO
      END DO


!*************** HNO4 photodissociation

! cross section from JPL85 (identical to JPL92 and JPL94)

      j = j + 1
      jlabel(j) = 'HNO4 -> HO2 + NO2'
!      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO4.abs',STATUS='old')
      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO4_jpl92.abs',STATUS='old')
      DO i = 1, 4
         READ(kin,*)
      ENDDO
      n = 31
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

! quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

!*************** H2O2 photodissociation

! cross section from Lin et al. 1978

      j = j + 1
      jlabel(j) = 'H2O2 -> 2 OH'
!     OPEN(UNIT=kin,FILE='DATAJ1/ABS/H2O2_lin.abs',STATUS='old')
!     DO i = 1, 7
!        READ(kin,*)
!     ENDDO
!     n = 32
!     DO i = 1, n
!        READ(kin,*) x1(i), y1(i)
!        y1(i) = y1(i) * 1.E-20
!     ENDDO
!     CLOSE (kin)
!
!      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,               0.,0.)
!      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
!      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
!      IF (ierr .NE. 0) THEN
!         WRITE(*,*) ierr, jlabel(j)
!         STOP
!      ENDIF

! cross section from JPL94
! tabulated data up to 260 nm
      OPEN(UNIT=kin,FILE='DATAJ1/ABS/H2O2_jpl94.abs',STATUS='old')
      READ(kin,*) idum,n
      DO i = 1, idum-2
         READ(kin,*)
      ENDDO
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE (kin)


      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF


      A0 = 6.4761E+04
      A1 = -9.2170972E+02
      A2 = 4.535649
      A3 = -4.4589016E-03
      A4 = -4.035101E-05
      A5 = 1.6878206E-07
      A6 = -2.652014E-10
      A7 = 1.5534675E-13

      B0 = 6.8123E+03
      B1 = -5.1351E+01
      B2 = 1.1522E-01
      B3 = -3.0493E-05
      B4 = -1.0924E-07

! quantum yield = 1
      qy = 1.

      DO iw = 1, nw - 1

! Parameterization (JPL94)
! Range 260-350 nm; 200-400 K

         IF ((wl(iw) .GE. 260.) .AND. (wl(iw) .LT. 350.)) THEN

           lambda = wc(iw)
           sumA = ((((((A7*lambda + A6)*lambda + A5)*lambda +           &
     &                  A4)*lambda +A3)*lambda + A2)*lambda +           &
     &                  A1)*lambda + A0
           sumB = (((B4*lambda + B3)*lambda + B2)*lambda +              &
     &               B1)*lambda + B0

           DO i = 1, nz
              t = MIN(MAX(tlev(i),200.),400.)
              chi = 1./(1.+EXP(-1265./t))
              xs = (chi * sumA + (1.-chi)*sumB)*1E-21
              sq(j,i,iw) = xs*qy
           ENDDO
         ELSE
           DO i = 1, nz
              sq(j,i,iw) = yg(iw)*qy
           ENDDO
         ENDIF

      ENDDO

!*************** CH2O photodissociatation

      j = j+1
      jlabel(j) = 'CH2O -> H + HCO'

      j = j+1
      jlabel(j) = 'CH2O -> H2 + CO'

! cross section data from Cantrell et al., as recommended in JPL 94.
! JPL recommended values are averaged over 2.5 nm intervals from 301.25-356.25
! nm;  a temperature dependence is given between -50 and +20 degrees C.


!     OPEN(UNIT=kin,FILE='DATAJ1/ABS/CH2O_jpl94.abs',STATUS='OLD')
!     READ(kin,*) idum, n
!     DO i = 1, idum-2
!       READ(kin,*)
!     ENDDO
!     DO i = 1, n
!       READ(kin,*) x1(i),dum,dum,y1(i),y2(i)
!       x1(i) = x1(i)-1.25
!       x2(i) = x1(i)
!       y1(i) = y1(i)*2.5
!       y2(i) = y2(i)*2.5
!     ENDDO
!     x1(n+1) = x1(n)+1.25
!     x2(n+1) = x1(n+1)
!     n = n+1
!     n1 = n
!     n2 = n
!     CLOSE(kin)

!     CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax), 0.)
!     CALL addpnt(x1,y1,kdata,n1,               0., 0.)
!     CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
!     CALL addpnt(x1,y1,kdata,n1,           1.e+38, 0.)
!     CALL inter3(nw,wl,yg1,n1,x1,y1,ierr)
!     IF (ierr .NE. 0) THEN
!        WRITE(*,*) ierr, jlabel(j)
!        STOP
!     ENDIF

!     CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax), 0.)
!     CALL addpnt(x2,y2,kdata,n2,               0., 0.)
!     CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
!     CALL addpnt(x2,y2,kdata,n2,            1.e+38,0.)
!     CALL inter3(nw,wl,yg2,n2,x2,y2,ierr)
!     IF (ierr .NE. 0) THEN
!        WRITE(*,*) ierr, jlabel(j)
!        STOP
!     ENDIF

!     DO iw = 1, nw-1
!        yg1(iw) = yg1(iw)/(wl(iw+1)-wl(iw))
!        yg2(iw) = yg2(iw)/(wl(iw+1)-wl(iw))
!     ENDDO


! cross section (Cantrell et al., 199?, original data set)

      OPEN(kin,FILE='DATAJ1/ABS/CH2O_can.abs',STATUS='old')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i)
! data are on wavenumber grid (cm-1), so convert to wavelength in nm:
        x1(i) = 1./x1(i) * 1E7
! cross section assumed to be zero for wavelengths longer than 360 nm
        IF (x1(i) .GT. 360.) THEN
          y1(i) = 0.
          y2(i) = 0.
        ENDIF
      ENDDO
      CLOSE(kin)


! grid was on increasing wavenumbers, so need to reverse to get increasing
! wavelengths

      DO i = 1, n/2

         irev = n+1-i

         dum = x1(i)
         x1(i) = x1(irev)
         x1(irev) = dum

         dum = y1(i)
         y1(i) = y1(irev)
         y1(irev) = dum

         dum = y2(i)
         y2(i) = y2(irev)
         y2(irev) = dum

      ENDDO

      DO i = 1, n
        x2(i) = x1(i)
      ENDDO

      n1 = n
      n2 = n

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,             1E38,0.)

      CALL ch_inter2(nw,wl,yg3,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,              1E38,0.)

      CALL ch_inter2(nw,wl,yg4,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF


! quantum yields

      OPEN(UNIT=kin,FILE='DATAJ1/YLD/CH2O_i_mad.yld',STATUS='old')
      DO i = 1, 11
         READ(kin,*)
      ENDDO
      n = 20
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j-1)
         STOP
      ENDIF
      DO iw = 1, nw - 1
         DO i = 1, nz
!           sq(j-1,i,iw) = yg(iw)*yg1(iw)
            t = Max(223.15,MIN(tlev(i),293.15))-273.15
            xs = Max(yg3(iw)+yg4(iw)*t,0.)
            sq(j-1,i,iw) = xs * yg1(iw)
         ENDDO
      ENDDO

      OPEN(UNIT=kin,FILE='DATAJ1/YLD/CH2O_ii_mad.yld',STATUS='old')
      DO i = 1, 9
         READ(kin,*)
      ENDDO
      n = 33
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg2,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

! pressure and temmperature dependent for w > 330.
      DO iw = 1, nw - 1
         DO i = 1, nz
            IF ( (wc(iw) .GE. 330.) .AND. (yg2(iw) .GT. 0.) ) THEN
               phi1 = yg1(iw)
               phi2 = yg2(iw)
               phi20 = 1. - phi1
               ak300=((1./phi2)-(1./phi20))/2.54E+19
               akt=ak300*(1.+61.69*(1.-tlev(i)/300.)*(wc(iw)/329.-1.))
               qy = 1. / ( (1./phi20) + airlev(i)*akt)
            ELSE
               qy = yg2(iw)
            ENDIF
            qy = MAX(0.,qy)
            qy = MIN(1.,qy)
!           sq(j,i,iw) = yg(iw)*qy
            t = Max(223.15,MIN(tlev(i),293.15))-273.15
            xs = Max(yg3(iw)+yg4(iw)*t,0.)
            sq(j,i,iw) = xs * qy
         ENDDO
      ENDDO

!*************** CH3OOH photodissociation

! cross section from JPL 94 (identical to JPL92)

      j = j + 1
      jlabel(j) = 'CH3OOH -> CH3O + OH'
!      OPEN(UNIT=kin,FILE='DATAJ1/ABS/h011_jpl85.abs',STATUS='old')
!     OPEN(UNIT=kin,FILE='DATAJ1/ABS/h011_jpl92.abs',STATUS='old')
      OPEN(UNIT=kin,FILE='DATAJ1/ABS/CH3OOH_jpl94.abs',STATUS='old')
      READ(kin,*) idum, n
      DO i = 1, idum-2
         READ(kin,*)
      ENDDO
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

! quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

!*************** CH3ONO2 photodissociation

! cross section from Calvert and  Pitts

      j = j + 1
      jlabel(j) = 'CH3ONO2 -> CH3O+NO2'
      OPEN(UNIT=kin,FILE='DATAJ1/ABS/n011_cp.abs',STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      n = 15
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

! quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

!*************** CH3CHO photodissociation

! cross section from Calvert and  Pitts

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/d021_cp.abs',STATUS='old')
      DO i = 1, 14
         READ(kin,*)
      ENDDO
      n = 54
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         x1(i) = x1(i)/10.
         y1(i) = y1(i) * 3.82E-21
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

! quantum yield

      j = j + 1
      jlabel(j) = 'CH3CHO -> CH3+CHO'
      OPEN(UNIT=kin,FILE='DATAJ1/YLD/d021_i.yld',STATUS='old')
      DO i = 1, 18
         READ(kin,*)
      ENDDO
      n = 10
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            IF (yg1(iw).NE.0.) THEN
               qy = 1./(1.+ (1./yg1(iw) -1.)*airlev(i)/2.465E19)
            ELSE
               qy = 0.
            ENDIF
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

      j = j + 1
      jlabel(j) = 'CH3CHO -> ii'
      OPEN(UNIT=kin,FILE='DATAJ1/YLD/d021_ii.yld',STATUS='old')
      DO i = 1, 10
         READ(kin,*)
      ENDDO
      n = 9
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*yg1(iw)
         ENDDO
      ENDDO

      j = j + 1
      jlabel(j) = 'CH3CHO -> iii'
      OPEN(UNIT=kin,FILE='DATAJ1/YLD/d021_iii.yld',STATUS='old')
      DO i = 1, 10
         READ(kin,*)
      ENDDO
      n = 9
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*yg1(iw)
         ENDDO
      ENDDO

!*************** CH3COCH3 photodissociation

! cross section from Calvert and  Pitts

      j = j + 1
      jlabel(j) = 'CH3COCH3 '
      OPEN(UNIT=kin,FILE='DATAJ1/ABS/k031_cp.abs',STATUS='old')
      DO i = 1, 6
         READ(kin,*)
      ENDDO
      n = 35
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 3.82E-21
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF
! quantum yield from Gardiner et al. 1984:

      DO iw = 1, nw - 1
         DO i = 1, nz
            qy = 0.0766 + 0.09415*EXP(-airlev(i)/3.222e18)
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

!*************** PAN photodissociation

      j = j+1
      jlabel(j) = 'PAN + hv -> products'

! cross section from Senum et al., 1984, J.Phys.Chem. 88/7, 1269-1270

!     OPEN(UNIT=kin,FILE='DATAJ1/ABS/PAN_senum.abs',STATUS='OLD')
!     DO i = 1, 14
!        READ(kin,*)
!     ENDDO
!     n = 21
!     DO i = 1, n
!        READ(kin,*) x1(i), y1(i)
!        y1(i) = y1(i) * 1.E-20
!     ENDDO
!     CLOSE(kin)

!      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,               0.,0.)
!      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
!      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
!      IF (ierr .NE. 0) THEN
!         WRITE(*,*) ierr, jlabel(j)
!         STOP
!      ENDIF

! cross section from
!      Talukdar et al., 1995, J.Geophys.Res. 100/D7, 14163-14174

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/PAN_talukdar.abs',STATUS='OLD')
      DO i = 1, 14
         READ(kin,*)
      ENDDO
      n = 78
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), y2(i)
         y1(i) = y1(i) * 1.E-20
         y2(i) = y2(i) * 1E-3
         x2(i) = x1(i)
      ENDDO
      n2 = n
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,          0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,      1.e+38,0.)
      CALL ch_inter2(nw,wl,yg2,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

! quantum yield
! yet unknown, but assumed to be 1.0 (Talukdar et al., 1995)

      qy = 1.0

      DO iw = 1, nw-1
        DO i = 1, nz
          xs = yg(iw) * EXP(yg2(iw)*(tlev(i)-298.))
          sq(j,i,iw) = qy * xs
        ENDDO
      ENDDO

!*************************************************************
!************ CCl2O photodissociation

      j = j+1
      jlabel(j) = 'CCl2O + hv -> Products'

!** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CCl2O_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!** quantum yield unity (Calvert and Pitts)
      qy = 1.
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CCl4 photodissociation

      j = j+1
      jlabel(j) = 'CCl4 + hv -> Products'

!** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CCl4_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!** quantum yield assumed to be unity
      qy = 1.
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CClFO photodissociation

      j = j+1
      jlabel(j) = 'CClFO + hv -> Products'

!** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CClFO_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!** quantum yield unity
      qy = 1.
      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CF2O photodissociation

      j = j+1
      jlabel(j) = 'CF2O + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CF2O_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!** quantum yield unity (Nolle et al.)
      qy = 1.
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CF2ClCFCl2 (CFC-113) photodissociation

      j = j+1
      jlabel(j) = 'CF2ClCFCl2 (CFC-113) + hv -> Products'

!** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CFC-113_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        x2(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n

!* sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,         1E38,0.)

      CALL ch_inter2(nw,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

! sigma @ 210 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,         1E38,0.)

      CALL ch_inter2(nw,wl,yg2,n2,x2,y2,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = MAX(210.,MIN(tlev(iz),295.))
        slope = (t-210.)/(295.-210.)
        DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg2(iw) + slope*(yg1(iw)-yg2(iw)))
        ENDDO
      ENDDO

!*************************************************************
!************ CF2ClCF2Cl (CFC-114) photodissociation

      j = j+1
      jlabel(j) = 'CF2ClCF2Cl (CFC-114) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CFC-114_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        x2(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n

!* sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,         1E38,0.)

      CALL ch_inter2(nw,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

! sigma @ 210 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,         1E38,0.)

      CALL ch_inter2(nw,wl,yg2,n2,x2,y2,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = MAX(210.,MIN(tlev(iz),295.))
        slope = (t-210.)/(295.-210.)
        DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg2(iw) + slope*(yg1(iw)-yg2(iw)))
        ENDDO
      ENDDO

!*************************************************************
!************ CF3CF2Cl (CFC-115) photodissociation

      j = j+1
      jlabel(j) = 'CF3CF2Cl (CFC-115) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CFC-115_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO


!*************************************************************
!************ CCl3F (CFC-11) photodissociation

      j = j+1
      jlabel(j) = 'CCl3F (CFC-11) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CFC-11_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

!* sigma @ 298 K

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = 1E-04 * (tlev(iz)-298.)
        DO iw = 1, nw-1
          sq(j,iz,iw) = qy * yg(iw) * EXP((wc(iw)-184.9) * t)
        ENDDO
      ENDDO

!*************************************************************
!************ CCl2F2 (CFC-12) photodissociation

      j = j+1
      jlabel(j) = 'CCl2F2 (CFC-12) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CFC-12_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

!* sigma @ 298 K

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = 1E-04 * (tlev(i)-298.)
        DO iw = 1, nw-1
          sq(j,iz,iw) = qy * yg(iw) * EXP((wc(iw)-184.9) * t)
        ENDDO
      ENDDO

!*************************************************************
!************ CH3Br photodissociation

      j = j+1
      jlabel(j) = 'CH3Br + hv -> Products'
      OPEN(kin,FILE='DATAJ1/ABS/CH3Br_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CH3CCl3 photodissociation

      j = j+1
      jlabel(j) = 'CH3CCl3 + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CH3CCl3_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        y3(i) = y3(i) * 1E-20
        x2(i) = x1(i)
        x3(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n

!* sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,         1E38,0.)

      CALL ch_inter2(nw,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!* sigma @ 250 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,         1E38,0.)

      CALL ch_inter2(nw,wl,yg2,n2,x2,y2,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!* sigma @ 210 K

      CALL addpnt(x3,y3,kdata,n3, x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,           0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,         1E38,0.)

      CALL ch_inter2(nw,wl,yg3,n3,x3,y3,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = MIN(295.,MAX(tlev(iz),210.))
        IF (t .LE. 250.) THEN
          slope = (t-210.)/(250.-210.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg3(iw) + slope*(yg2(iw)-yg3(iw)))
          ENDDO
        ELSE
          slope = (t-250.)/(295.-250.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg2(iw) + slope*(yg1(iw)-yg2(iw)))
          ENDDO
        ENDIF
      ENDDO

!*************************************************************
!************ CH3Cl photodissociation

      j = j+1
      jlabel(j) = 'CH3Cl + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CH3Cl_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        y3(i) = y3(i) * 1E-20
        x2(i) = x1(i)
        x3(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n

!* sigma @ 296 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,         1E38,0.)

      CALL ch_inter2(nw,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!* sigma @ 279 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,         1E38,0.)

      CALL ch_inter2(nw,wl,yg2,n2,x2,y2,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!* sigma @ 255 K

      CALL addpnt(x3,y3,kdata,n3, x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,           0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,         1E38,0.)

      CALL ch_inter2(nw,wl,yg3,n3,x3,y3,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = MAX(255.,MIN(tlev(i),296.))
        IF (t .LE. 279.) THEN
          slope = (t-255.)/(279.-255.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg3(iw)+slope*(yg2(iw)-yg3(iw)))
          ENDDO
        ELSE
          slope = (t-279.)/(296.-279.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg2(iw)+slope*(yg1(iw)-yg2(iw)))
          ENDDO
        ENDIF
      ENDDO

!*************************************************************
!************ ClOO photodissociation

      j = j+1
      jlabel(j) = 'ClOO + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/ClOO_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CF3CHCl2 (HCFC-123) photodissociation

      j = j+1
      jlabel(j) = 'CF3CHCl2 (HCFC-123) + hv -> Products'

!*** cross sections from JPL94 recommendation

!     OPEN(kin,FILE='DATAJ1/ABS/HCFC-123_jpl94.abs',STATUS='OLD')
!     READ(kin,*) idum, n
!     DO i = 1, idum-2
!       READ(kin,*)
!     ENDDO
!     DO i = 1, n
!       READ(kin,*) x1(i), y1(i)
!       y1(i) = y1(i) * 1E-20
!     ENDDO
!     CLOSE(kin)

!     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
!     CALL addpnt(x1,y1,kdata,n,          0.,0.)
!     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
!     CALL addpnt(x1,y1,kdata,n,        1E38,0.)

!     CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

!     IF (ierr .NE. 0) THEN
!        WRITE(*,*) ierr, jlabel(j)
!        STOP
!     ENDIF

!*** quantum yield assumed to be unity
!     qy = 1.

!     DO iw = 1, nw-1
!       DO iz = 1, nz
!         sq(j,iz,iw) = qy * yg(iw)
!       ENDDO
!     ENDDO


!*** cross section from Orlando et al., 1991

      OPEN(kin,FILE='DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,('(A)')) inline
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      LBar = 206.214

!*** quantum yield assumed to be unity

      qy = 1.

      DO iw = 1, nw-1
        lambda = wc(iw)
        IF (lambda .GE. 190. .AND. lambda .LE. 230.) THEN
          DO iz = 1, nz
             t = MIN(295.,MAX(tlev(i),203.))-TBar
             sum = 0.
             DO i = 1, 4
                sum = (coeff(i,1)+t*(coeff(i,2)+t*coeff(i,3))) *        &
     &                (lambda-LBar)**(i-1) + sum
             ENDDO
             sq(j,iz,iw) = qy * EXP(sum)
          ENDDO
        ELSE
          DO iz = 1, nz
            sq(j,iz,iw) = 0.
          ENDDO
        ENDIF
      ENDDO

!*************************************************************
!************ CF3CHFCl (HCFC-124) photodissociation

      j = j+1
      jlabel(j) = 'CF3CHFCl (HCFC-124) + hv -> Products'

!*** cross sections from JPL94 recommendation

!     OPEN(kin,FILE='DATAJ1/ABS/HCFC-124_jpl94.abs',STATUS='OLD')
!     READ(kin,*) idum, n
!     DO i = 1, idum-2
!       READ(kin,*)
!     ENDDO
!     DO i = 1, n
!       READ(kin,*) x1(i), y1(i)
!       y1(i) = y1(i) * 1E-20
!     ENDDO
!     CLOSE(kin)

!     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
!     CALL addpnt(x1,y1,kdata,n,          0.,0.)
!     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
!     CALL addpnt(x1,y1,kdata,n,        1E38,0.)

!     CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

!     IF (ierr .NE. 0) THEN
!       WRITE(*,*) ierr, jlabel(j)
!       STOP
!     ENDIF

!*** quantum yield assumed to be unity
!     qy = 1.

!     DO iw = 1, nw-1
!       DO iz = 1, nz
!         sq(j,iz,iw) = qy * yg(iw)
!       ENDDO
!     ENDDO

!*** cross section from Orlando et al., 1991

      OPEN(kin,FILE='DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      idum = idum+5
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,('(A)')) inline
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      LBar = 206.214

!*** quantum yield assumed to be unity

      qy = 1.

      DO iw = 1, nw-1
        lambda = wc(iw)
        IF (lambda .GE. 190. .AND. lambda .LE. 230.) THEN
          DO iz = 1, nz
             t = MIN(295.,MAX(tlev(i),203.))-TBar
             sum = 0.
             DO i = 1, 4
                sum = (coeff(i,1)+t*(coeff(i,2)+t*coeff(i,3))) *        &
     &                (lambda-LBar)**(i-1) + sum
             ENDDO
             sq(j,iz,iw) = qy * EXP(sum)
          ENDDO
        ELSE
          DO iz = 1, nz
            sq(j,iz,iw) = 0.
          ENDDO
        ENDIF
      ENDDO

!*************************************************************
!************ CH3CFCl2 (HCFC-141b) photodissociation

      j = j+1
      jlabel(j) = 'CH3CFCl2 (HCFC-141b) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-141b_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CH3CF2Cl (HCFC-142b) photodissociation

      j = j+1
      jlabel(j) = 'CH3CF2Cl (HCFC-142b) + hv -> Products'

!*** cross sections from JPL94 recommendation

!     OPEN(kin,FILE='DATAJ1/ABS/HCFC-142b_jpl94.abs',STATUS='OLD')
!     READ(kin,*) idum, n
!     DO i = 1, idum-2
!       READ(kin,*)
!     ENDDO
!     DO i = 1, n
!       READ(kin,*) x1(i), y1(i)
!       y1(i) = y1(i) * 1E-20
!     ENDDO
!     CLOSE(kin)

!     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
!     CALL addpnt(x1,y1,kdata,n,          0.,0.)
!     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
!     CALL addpnt(x1,y1,kdata,n,        1E38,0.)

!     CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

!     IF (ierr .NE. 0) THEN
!       WRITE(*,*) ierr, jlabel(j)
!       STOP
!     ENDIF

!*** quantum yield assumed to be unity
!     qy = 1.

!     DO iw = 1, nw-1
!       DO iz = 1, nz
!         sq(j,iz,iw) = qy * yg(iw)
!       ENDDO
!     ENDDO

!*** cross section from Orlando et al., 1991

      OPEN(kin,FILE='DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      idum = idum+10
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,('(A)')) inline
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      LBar = 206.214

!*** quantum yield assumed to be unity

      qy = 1.

      DO iw = 1, nw-1
        lambda = wc(iw)
        IF (lambda .GE. 190. .AND. lambda .LE. 230.) THEN
          DO iz = 1, nz
             t = MIN(295.,MAX(tlev(i),203.))-TBar
             sum = 0.
             DO i = 1, 4
                sum = (coeff(i,1)+t*(coeff(i,2)+t*coeff(i,3))) *        &
     &                (lambda-LBar)**(i-1) + sum
             ENDDO
             sq(j,iz,iw) = qy * EXP(sum)
          ENDDO
        ELSE
          DO iz = 1, nz
            sq(j,iz,iw) = 0.
          ENDDO
        ENDIF
      ENDDO

!*************************************************************
!************ CF3CF2CHCl2 (HCFC-225ca) photodissociation

      j = j+1
      jlabel(j) = 'CF3CF2CHCl2 (HCFC-225ca) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-225ca_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CF2ClCF2CHFCl (HCFC-225cb) photodissociation

      j = j+1
      jlabel(j) = 'CF2ClCF2CHFCl (HCFC-225cb) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-225cb_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CHClF2 (HCFC-22) photodissociation

      j = j+1
      jlabel(j) = 'CHClF2 (HCFC-22) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-22_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i), y5(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        y3(i) = y3(i) * 1E-20
        y4(i) = y4(i) * 1E-20
        y5(i) = y5(i) * 1E-20
        x2(i) = x1(i)
        x3(i) = x1(i)
        x4(i) = x1(i)
        x5(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n
      n4 = n
      n5 = n

!* sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,         1E38,0.)

      CALL ch_inter2(nw,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!* sigma @ 270 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,         1E38,0.)

      CALL ch_inter2(nw,wl,yg2,n2,x2,y2,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!* sigma @ 250 K

      CALL addpnt(x3,y3,kdata,n3, x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,           0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,         1E38,0.)

      CALL ch_inter2(nw,wl,yg3,n3,x3,y3,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!* sigma @ 230 K

      CALL addpnt(x4,y4,kdata,n4, x4(1)*(1.-deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,           0.,0.)
      CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,         1E38,0.)

      CALL ch_inter2(nw,wl,yg4,n4,x4,y4,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!* sigma @ 210 K

      CALL addpnt(x5,y5,kdata,n5, x5(1)*(1.-deltax),0.)
      CALL addpnt(x5,y5,kdata,n5,           0.,0.)
      CALL addpnt(x5,y5,kdata,n5,x5(n5)*(1.+deltax),0.)
      CALL addpnt(x5,y5,kdata,n5,         1E38,0.)

      CALL ch_inter2(nw,wl,yg5,n5,x5,y5,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
         t = MIN(295.,MAX(tlev(iz),210.))
         IF (t .LE. 230.) THEN
            slope = (t-210.)/(230.-210.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = qy * (yg5(iw)+slope*(yg4(iw)-yg5(iw)))
            ENDDO
         ELSEIF (t .LE. 250.) THEN
            slope = (t-230.)/(250.-230.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = qy * (yg4(iw)+slope*(yg3(iw)-yg4(iw)))
            ENDDO
         ELSEIF (t .LE. 270.) THEN
            slope = (t-250.)/(270.-250.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = qy * (yg3(iw)+slope*(yg2(iw)-yg3(iw)))
            ENDDO
         ELSE
            slope = (t-270.)/(295.-270.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = qy * (yg2(iw)+slope*(yg1(iw)-yg2(iw)))
            ENDDO
         ENDIF
      ENDDO

!*************************************************************
!************ HO2 photodissociation

      j = j+1
      jlabel(j) = 'HO2 + hv -> OH + O'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/HO2_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield:  absolute quantum yield has not been reported yet, but
!***                 Lee measured a quantum yield for O(1D) production at 248
!***                 nm that was 15 time larger than at 193 nm
!*** here:  a quantum yield of unity is assumed at 248 nm and beyond, for
!***        shorter wavelengths a linear decrease with lambda is assumed

      DO iw = 1, nw-1
         IF (wc(iw) .GE. 248.) THEN
            qy = 1.
         ELSE
            qy = 1./15. + (wc(iw)-193.)*(14./15.)/(248.-193.)
            qy = MAX(qy,0.)
         ENDIF
         DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
         ENDDO
      ENDDO

!*************************************************************
!************ CF2Br2 (Halon-1202) photodissociation

      j = j+1
      jlabel(j) = 'CF2Br2 (Halon-1202) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/Halon-1202_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield unity (Molina and Molina)
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CF2BrCl (Halon-1211) photodissociation

      j = j+1
      jlabel(j) = 'CF2BrCl (Halon-1211) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/Halon-1211_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CF3Br (Halon-1301) photodissociation

      j = j+1
      jlabel(j) = 'CF3Br (Halon-1301) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/Halon-1301_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ CF2BrCF2Br (Halon-2402) photodissociation

      j = j+1
      jlabel(j) = 'CF2BrCF2Br (Halon-2402) + hv -> Products'

!*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/Halon-2402_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)+deltax,0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL ch_inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

!*** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

!*************************************************************
!************ N2O photodissociation

      j = j+1
      jlabel(j) = 'N2O + hv -> N2 + O(1D)'

!*** cross sections according to JPL94 recommendation
!*** see file DATAJ1/ABS/N2O_jpl94.abs for detail

      A0 = 68.21023
      A1 = -4.071805
      A2 = 4.301146E-02
      A3 = -1.777846E-04
      A4 = 2.520672E-07

      B0 = 123.4014
      B1 = -2.116255
      B2 = 1.111572E-02
      B3 = -1.881058E-05

!*** quantum yield of N(4s) and NO(2Pi) is less than 1% (Greenblatt and
!*** Ravishankara), so quantum yield of O(1D) is assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
         lambda = wc(iw)
         IF (lambda .GE. 173. .AND. lambda .LE. 240.) THEN
           DO iz = 1, nz
             t = MAX(194.,MIN(tlev(iz),320.))
             A = (((A4*lambda+A3)*lambda+A2)*lambda+A1)*lambda+A0
             B = (((B3*lambda+B2)*lambda+B1)*lambda+B0)
             B = (t-300.)*EXP(B)
             sq(j,iz,iw) = qy * EXP(A+B)
           ENDDO
         ELSE
           DO iz = 1, nz
             sq(j,iz,iw) = 0.
           ENDDO
         ENDIF
      ENDDO




!***************************************************************
!***************************************************************

      IF (j .GT. kj) STOP '1002'
      IF (LHOOK) CALL DR_HOOK('JSPEC1',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE JSPEC1
