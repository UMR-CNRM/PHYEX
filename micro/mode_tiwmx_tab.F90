MODULE MODE_TIWMX_TAB
IMPLICIT NONE
CONTAINS
FUNCTION TIWMX_TAB(P,T,QR,FICE,QRSN,RS,EPS)

!     Purpose:  (*)
!     The fuction tiwmx_tab returns the wet bulb temperature, but also the
!     corresponding saturation specific humidity as an output parameter.
!     -------------------------------------------------------------------
!     Computation of wet bulb temperature.
!     TIW is found iteratively. Note that Q is used instead of WV.
!     Converges VERY quickly. EPS is the threshold value for
!     TIW(n)-TIW(n-1) at wich the iteration is interupted. (n is
!     iteration number).
!     Modified in July 1988 by Stefan Gollvik
!     Converted to new fortran standard in Dec. 2006 by Karl-Ivar Ivarsson.

!     INTERFACE  :  the function is intended to be used everywhere.
!     INPUT  arguments  (arguments d'entree)
!     -----------------------------------------------------
!     P        : pressure  (Pa)
!     T        : temperature (K)
!     QR       : mixing ratio humidity (kg/Kg)
!     FICE     : fraction of ice (0 to 1)
!     EPS      : The value determines the accuracy of the output value. 0.1
!              : is suffient in most cases. Low value means high accuracy but
!              : also high computational cost.

!     OUTPUT  arguments  (arguments d'sortie)
!     -----------------------------------------------------
!     RS       : saturation mixing ratio (fice determines if it is over ice or water)
!                
!     QRSN     : saturation mixing ratio for the wet bulb temperature.
!                (Kg/Kg)
!     ( the function itself is the wet bulb temperature (K) )

!     Work  variables  :
!     -----------------------------------------------------
!      f        : temporary variable , temperature (K)
!      dfdt     : temporary variable   (K/K)
!      t2       : temperature used in iteration (K)
!      dt       : temperature residual (K)
!      dqsdt    : d(qsat)/d(T) in iteration (1/K)
!      b        : (latent heat)/(heat capacity for dry air)  (K)
!      iter     : iteration number

!     1. Declarations.
!     ==================================================================
!     1.1 MODULES USED
 
  USE MODD_CST, ONLY  : XEPSILO, XCPD, XLSTT, XLVTT
  USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
  USE MODE_TIWMX, ONLY : ESATI, ESATW, DESDTI, DESDTW
  
  IMPLICIT NONE

!  Function name :
  REAL :: TIWMX_TAB

!  Input Arguments
  REAL, INTENT(IN) :: P,T,QR,FICE,EPS

!  Output Arguments
  REAL, INTENT(OUT) :: QRSN,RS

!  Work variables :
  REAL :: F,DFDT,T2,DT,QSN,DQSDT,B
  REAL :: ZES,ZDESDT
  INTEGER :: ITER
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('TIWMX_TAB',0,ZHOOK_HANDLE)

  T2 = T

  B = ( XLVTT*(1.-FICE) + FICE*XLSTT )/XCPD

  TIWMX_TAB = T2
  DO ITER=1,10
     ZES = ESATI(T2)*FICE + ESATW(T2)*(1.-FICE)
     ZDESDT= DESDTI(T2)*FICE + DESDTW(T2)*(1.-FICE)
     IF(ZES >= P*0.61)THEN ! Do not to compute when condensation 
        QRSN = 1.          ! not possible and avoid mixing ratio > 1  
        TIWMX_TAB = T2
        IF (LHOOK) CALL DR_HOOK('TIWMX_TAB',1,ZHOOK_HANDLE)
        RETURN 
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
        TIWMX_TAB = T2
        QSN = MIN(1.,QSN+DT*DQSDT)
        QRSN = QSN ! approximation
        IF (LHOOK) CALL DR_HOOK('TIWMX_TAB',1,ZHOOK_HANDLE)
        RETURN
     ENDIF
  ENDDO

  IF (LHOOK) CALL DR_HOOK('TIWMX_TAB',1,ZHOOK_HANDLE)

END FUNCTION TIWMX_TAB
END MODULE MODE_TIWMX_TAB
