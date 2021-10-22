!     ######spl
      SUBROUTINE sphers(nz, z, zen, dsdh, nid)
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

! input
      INTEGER nz
      REAL zen, z(kz)

! output
      INTEGER nid(0:kz)
      REAL dsdh(0:kz,kz)

! more program constants
      REAL re, ze(kz)
      REAL  dr
      PARAMETER ( dr = pi/180.)

! local

      REAL zenrad, rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm
      INTEGER i, j, k
      INTEGER id

      INTEGER nlayer
      REAL zd(0:kz-1)


      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SPHERS',0,ZHOOK_HANDLE)
      zenrad = zen*dr

! number of layers:
      nlayer = nz - 1

! include the elevation above sea level to the radius of the earth:
      re = radius + z(1)
! correspondingly z changed to the elevation above earth surface:
      DO k = 1, nz
         ze(k) = z(k) - z(1)
      END DO

! inverse coordinate of z
      zd(0) = ze(nz)
      DO k = 1, nlayer
        zd(k) = ze(nz - k)
      END DO

! initialize dsdh(i,j), nid(i)
      DO i = 0, kz
       nid(i) = 0
       DO j = 1, kz
        dsdh(i,j) = 0.
       END DO
      END DO

! calculate ds/dh of every layer
      DO 100 i = 0, nlayer

        rpsinz = (re + zd(i)) * SIN(zenrad)

        IF ( (zen .GT. 90.0) .AND. (rpsinz .LT. re) ) THEN
           nid(i) = -1
        ELSE

!
! Find index of layer in which the screening height lies
!
           id = i
           IF( zen .GT. 90.0 ) THEN
              DO 10 j = 1, nlayer
                 IF( (rpsinz .LT. ( zd(j-1) + re ) ) .AND.              &
     &               (rpsinz .GE. ( zd(j) + re )) ) id = j
 10           CONTINUE
           END IF

           DO 20 j = 1, id

             sm = 1.0
             IF(j .EQ. id .AND. id .EQ. i .AND. zen .GT. 90.0)          &
     &          sm = -1.0

             rj = re + zd(j-1)
             rjp1 = re + zd(j)

             dhj = zd(j-1) - zd(j)

             ga = rj*rj - rpsinz*rpsinz
             gb = rjp1*rjp1 - rpsinz*rpsinz
             IF (ga .LT. 0.0) ga = 0.0
             IF (gb .LT. 0.0) gb = 0.0

             IF(id.GT.i .AND. j.EQ.id) THEN
                dsj = SQRT( ga )
             ELSE
                dsj = SQRT( ga ) - sm*SQRT( gb )
             END IF
             dsdh(i,j) = dsj / dhj

 20        CONTINUE

           nid(i) = id

        END IF

 100  CONTINUE


      IF (LHOOK) CALL DR_HOOK('SPHERS',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE SPHERS
