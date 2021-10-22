!     ######spl
      SUBROUTINE sundis(idate,esrm2)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE

! input:
      INTEGER idate

! output:
      REAL esrm2

! internal:
      INTEGER iyear, imonth, iday, mday, month, jday
      REAL dayn, thet0
      REAL sinth, costh, sin2th, cos2th
      INTEGER imn(12)

      REAL pi
      PARAMETER(pi=3.1415926535898)
!_______________________________________________________________________

      DATA imn/31,28,31,30,31,30,31,31,30,31,30,31/
!_______________________________________________________________________

! parse date to find day number (Julian day)

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SUNDIS',0,ZHOOK_HANDLE)
      iyear = int(idate/10000)
      imonth = int( (idate-10000*iyear)/100 )
      iday = idate - (10000*iyear + 100*imonth)
      IF ( MOD(iyear,4) .EQ. 0) THEN
         imn(2) = 29
      ELSE
         imn(2) = 28
      ENDIF

      mday = 0
      DO 12, month = 1, imonth-1
         mday = mday + imn(month)
   12 CONTINUE
      jday = mday + iday
      dayn = FLOAT(jday - 1) + 0.5

! define angular day number and compute esrm2:

      thet0 = 2.*pi*dayn/365.

! calculate SIN(2*thet0), COS(2*thet0) from
! addition theoremes for trig functions for better
! performance;  the computation of sin2th, cos2th
! is about 5-6 times faster than the evaluation
! of the intrinsic functions SIN and COS
!
      sinth = SIN(thet0)
      costh = COS(thet0)
      sin2th = 2.*sinth*costh
      cos2th = costh*costh - sinth*sinth
      esrm2  = 1.000110 +                                               &
     &         0.034221*costh  +  0.001280*sinth +                      &
     &         0.000719*cos2th +  0.000077*sin2th
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('SUNDIS',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE SUNDIS
