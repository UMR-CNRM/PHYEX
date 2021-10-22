!     ######spl
      FUNCTION GAMMA_INC(PA,PX)  RESULT(PGAMMA_INC)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #############################################
!
!
!!****  *GAMMA_INC * -  Generalized gamma  function
!!
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the generalized
!!   incomplete Gamma function of its argument.
!!
!!                             /X
!!                       1     |
!!    GAMMA_INC(A,X)= -------- | Z**(A-1) EXP(-Z) dZ
!!                    GAMMA(A) |
!!                             /0
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODULE MODI_GAMMA : computation of the Gamma function
!!
!!    REFERENCE
!!    ---------
!!      Press, Teukolsky, Vetterling and Flannery: Numerical Recipes, 209-213
!!
!!
!!    AUTHOR
!!    ------
!!     Jean-Pierre Pinty *LA/OMP*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     7/12/95
!
!*       0. DECLARATIONS
!           ------------
!
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
REAL, INTENT(IN)                     :: PA
REAL, INTENT(IN)                     :: PX
REAL                                 :: PGAMMA_INC
!
!*       0.2 declarations of local variables
!
INTEGER                              :: JN
INTEGER                              :: ITMAX=100
REAL                                 :: ZEPS=3.E-7
REAL                                 :: ZFPMIN=1.E-30
REAL                                 :: ZAP,ZDEL,ZSUM
REAL                                 :: ZAN,ZB,ZC,ZD,ZH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GAMMA_INC',0,ZHOOK_HANDLE)
IF( (PX.LT.0.0).OR.(PA.LE.0.0) ) THEN
  PRINT *,' BAD ARGUMENTS IN GAMMA_INC'
!callabortstop
CALL ABORT
  STOP
END IF
!
IF( (PX.LT.PA+1.0) ) THEN
  ZAP = PA
  ZSUM = 1.0/PA
  ZDEL = ZSUM
  JN = 1
!
  LOOP_SERIES: DO
    ZAP = ZAP +1.0
    ZDEL = ZDEL*PX/ZAP
    ZSUM = ZSUM + ZDEL
    IF( ABS(ZDEL).LT.ABS(ZSUM)*ZEPS ) EXIT LOOP_SERIES
    JN = JN + 1
    IF( JN.GT.ITMAX ) THEN
      PRINT *,' ARGUMENT "PA" IS TOO LARGE OR "ITMAX" IS TOO SMALL, THE      &
              & INCOMPLETE GAMMA_INC FUNCTION CANNOT BE EVALUATED CORRECTLY  &
              & BY THE SERIES METHOD'
!callabortstop
CALL ABORT
      STOP
    END IF
  END DO LOOP_SERIES
  PGAMMA_INC = ZSUM * EXP( -PX+PA*ALOG(PX)-ALOG(GAMMA(PA)) )
!
  ELSE
!
  ZB = PX + 1.0 - PA
  ZC = 1.0/TINY(PX)
  ZD = 1.0/ZB
  ZH = ZD
  JN = 1
!
  LOOP_FRACTION: DO
    ZAN = -FLOAT(JN)*(FLOAT(JN)-PA)
    ZB = ZB + 2.0
    ZD = ZAN*ZD + ZB
    IF( ABS(ZD).LT.TINY(PX) ) THEN
      ZD = ZFPMIN
    END IF
    ZC = ZB + ZAN/ZC
    IF( ABS(ZC).LT.TINY(PX) ) THEN
      ZC = ZFPMIN
    END IF
    ZD = 1.0/ZD
    ZDEL = ZD*ZC
    ZH = ZH*ZDEL
    IF( ABS(ZDEL-1.0).LT.ZEPS ) EXIT LOOP_FRACTION
    JN = JN + 1
    IF( JN.GT.ITMAX ) THEN
      PRINT *,' ARGUMENT "PA" IS TOO LARGE OR "ITMAX" IS TOO SMALL, THE      &
              & INCOMPLETE GAMMA_INC FUNCTION CANNOT BE EVALUATED CORRECTLY  &
              & BY THE CONTINUOUS FRACTION METHOD'
!callabortstop
CALL ABORT
      STOP
    END IF
  END DO LOOP_FRACTION
  PGAMMA_INC = 1.0 - ZH*EXP( -PX+PA*ALOG(PX)-ALOG(GAMMA(PA)) )
!
END IF
!
IF (LHOOK) CALL DR_HOOK('GAMMA_INC',1,ZHOOK_HANDLE)
RETURN
!
END FUNCTION GAMMA_INC
