FUNCTION ARO_TIWMX(P,T,QR,FICE,QRSN,RS,EPS)

!     Purpose:  (*)
!     The fuction arotiwmx returns the wet bulb temperature, but also the
!     corresponding saturation specific humidity as an output parameter.
!     -------------------------------------------------------------------
!     Computation of wet bulb temperature.
!     TIW is found iteratively. Note that Q is used instead of WV.
!     Converges VERY quickly. EPS is the threshold value for
!     TIW(n)-TIW(n-1) at wich the iteration is interupted. (n is
!     iteration number).Modified in July 1988 by Stefan Gollvik
!     Convereted to new fortran standard in Dec. 2006 by Karl-Ivar Ivarsson.
!     Convereted to AROME in 2014 by Karl-Ivar Ivarsson.
!     INTERFACE  :  the function is intended to be used everywhere.
!     INPUT  arguments  (arguments d'entree)
!     -----------------------------------------------------
!     P        : pressure  (Pa)
!     T        : temperature (K)
!     QR       : mixing ratio humidity (kg/Kg)
!     FICE     : fraction of ice (0 to 1)
!     EPS      : The value determens the accuracy of the output value. 0.1
!              : is suffient in most cases. Low value means high accuracy but
!              : also high computational cost.

!     OUTPUT  arguments  (arguments d'sortie)
!     -----------------------------------------------------
!     RS        : saturation mixing ratio (fice determines if it is over ice or water)
!                
!     QRSN     : saturation mixing ratio for the wet bulb temperature.
!                (Kg/Kg)
!     ( the function itself is the wet bulb temperature (K) )

!     Work  variables  :
!     -----------------------------------------------------
!      f        : temporary variable , temperature (K)
!      dfdt     : temporary variable   (K/K)
!      t2       : temperature used in iteration (K)
!      thigh    : temperature ( t2 + EPS K )
!      tlow     : temperature ( t2 - EPS K )
!      dt       : temperature residual (K)
!      dqsdt    : approximative d(qsat)/d(T) in iteration (1/K)
!      b        : (latent heat)/(heat capacity for dry air)  (K)
!      iter     : iteration number

!     1. Declarations.
!     ==================================================================
!     1.1 MODULES USED
  USE MODD_CST,ONLY : XCPD,XCPV,XLVTT,XLSTT,XRD,XTT,XEPSILO, &
       &    XALPW,XBETAW,XGAMW,XALPI,XBETAI,XGAMI
  USE PARKIND1, ONLY : JPRB 
  USE YOMHOOK , ONLY : LHOOK, DR_HOOK
  IMPLICIT NONE

!  Function name :
  REAL                        :: ARO_TIWMX

!  Input Arguments
  REAL,  INTENT(IN)           ::  P,T,QR,FICE,EPS

!  Output Arguments
  REAL,  INTENT(OUT)          :: QRSN,RS

!  Work variables :
  REAL :: F,DFDT,T2,DT,QSN,DQSDT,B
  REAL :: ZB,ZG,ZES,ZDESDT
  INTEGER :: ITER
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('ARO_TIWMX',0,ZHOOK_HANDLE)

  T2 = T
!!  Q = QR !approximation

  B = ( XLVTT*(1.-FICE) + FICE*XLSTT )/XCPD

  ARO_TIWMX = T2
  DO ITER=1,10
     ZB = (XBETAI*FICE + XBETAW*(1.-FICE))/T2
     ZG = XGAMI*FICE + XGAMW*(1.-FICE)
     ZES = EXP(XALPI*FICE + XALPW*(1.-FICE) - ZB - ZG*ALOG(T2))
     ZDESDT = ZES*(ZB - ZG)/T2
     IF(ZES >= P)THEN
        QSN = 1.
        DQSDT = 0.
     ELSE
        QSN = XEPSILO*ZES/(P-ZES)
        DQSDT = QSN*ZDESDT*( 1.0/ZES + 1.0/(P-ZES) )
     ENDIF
     IF ( ITER == 1 ) RS = QSN
     F = T2 - T + B*(QSN - QR)
     DFDT = 1. + B*DQSDT
     DT = -F / DFDT
     T2 = T2 + DT
     IF(ABS(DT) <= EPS)THEN
        ARO_TIWMX = T2
        QSN = MIN(1.,QSN+DT*DQSDT)
        QRSN = QSN ! approximation
        IF (LHOOK) CALL DR_HOOK('ARO_TIWMX',1,ZHOOK_HANDLE)
        RETURN
     ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('ARO_TIWMX',1,ZHOOK_HANDLE)
END FUNCTION ARO_TIWMX
!!!!!!!!!!!!
FUNCTION AROQSATMX(P,T,FICE)
  USE MODD_CST,ONLY : XALPW,XBETAW,XGAMW,XALPI,XBETAI,XGAMI,XTT,XEPSILO 
  IMPLICIT NONE
  REAL ZES,ZQS,P,T,FICE,&
       &AROQSATMX

  ZES = EXP(XALPI*FICE + XALPW*(1.-FICE) - &
       & (XBETAI*FICE + XBETAW*(1.-FICE))/T -  &
       & (XGAMI*FICE + XGAMW*(1.-FICE))*ALOG(T))

  IF(ZES >= P)THEN ! temp > boiling point, condensation not possible.
                   ! Then this function lacks physical meaning, 
                   ! here set to one
     ZQS=1.
  ELSE
!     ZQS=0.622*ZES/(P-0.378*ZES) !q
     ZQS=XEPSILO*ZES/(P-ZES) !r
  ENDIF
  AROQSATMX=ZQS
END FUNCTION AROQSATMX
