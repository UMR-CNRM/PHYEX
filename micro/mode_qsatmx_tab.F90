MODULE MODE_QSATMX_TAB
IMPLICIT NONE
CONTAINS
FUNCTION QSATMX_TAB(P,T,FICE)

  USE PARKIND1, ONLY : JPRB
  USE MODD_CST ,ONLY : XEPSILO
  USE MODE_TIWMX, ONLY : ESATI,ESATW

  IMPLICIT NONE

  REAL :: QSATMX_TAB
  REAL, INTENT(IN) :: P,T,FICE

  REAL :: ZES

  ZES = ESATI(T)*FICE + ESATW(T)*(1.-FICE)
  IF(ZES >= P)THEN ! temp > boiling point, condensation not possible.
                   ! Then this function lacks physical meaning, 
                   ! here set to one
     QSATMX_TAB = 1.
  ELSE
     QSATMX_TAB = XEPSILO*ZES/(P-ZES) !r
  ENDIF

END FUNCTION QSATMX_TAB
END MODULE MODE_QSATMX_TAB
