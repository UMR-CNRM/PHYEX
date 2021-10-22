!     ######spl
      FUNCTION GY_M_V(KKA,KKU,KL,PY,PDYY,PDZZ,PDZY) RESULT(PGY_M_V)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##################################################
!
!!****  *GY_M_V * - Compute the gradient along y for a variable localized at
!!                  a mass point
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute a gradient along y
!     direction for a field PY localized at a mass point. The result PGY_M_V
!     is localized at a y-flux point (v point).
!
!                    (           ____________z )
!                    (               ________y )
!                 1  (                dzm(PY)  )
!   PGY_M_V =   ---- (dym(PY) - d*zy --------  )
!               d*yy (                 d*zz    )
!
!
!
!
!!**  METHOD
!!    ------
!!      We employ the Shuman operators to compute the derivatives and the
!!    averages. The metric coefficients PDYY,PDZY,PDZZ are dummy arguments.
!!
!!
!!    EXTERNAL
!!    --------
!!      FUNCTION DYM: compute a finite difference along the y direction for
!!      a variable at a mass localization
!!      FUNCTION DZM: compute a finite difference along the y direction for
!!      a variable at a mass localization
!!      FUNCTION MYM: compute an average in the x direction for a variable
!!      at a mass localization
!!      FUNCTION MZF: compute an average in the z direction for a variable
!!      at a flux side
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODD_CONF : LFLAT
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (function GY_M_V)
!!
!!
!!    AUTHOR
!!    ------
!!      P. Hereil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94
!!      Modification       16/03/95  change the order of the arguments
!!                         19/07/00  add the LFLAT switch + inlining(J. Stein)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_SHUMAN
USE MODD_CONF
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!              -------------------------------------
!
INTEGER,                INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                INTENT(IN)  :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY                   !d*yy
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZY                   !d*zy
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGY_M_V  ! result at flux
                                                              ! side
INTEGER  IJU,IKU,JJ,JK
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG Y
!              ----------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GY_M_V',0,ZHOOK_HANDLE)
IJU=SIZE(PY,2)
IKU=SIZE(PY,3)
IF (.NOT. LFLAT) THEN
!  PGY_M_V = (   DYM(PY)  -  MZF (   MYM(  DZM(PY) /PDZZ  ) * PDZY   )   )/PDYY
  DO JK=1+JPVEXT_TURB,IKU-JPVEXT_TURB
    DO JJ=1+JPHEXT,IJU
        PGY_M_V(:,JJ,JK)=                                                 &
           (  PY(:,JJ,JK)-PY(:,JJ-1,JK)                                   &
             -(  (PY(:,JJ,JK)-PY(:,JJ,JK-KL))     / PDZZ(:,JJ,JK)          &
                +(PY(:,JJ-1,JK)-PY(:,JJ-KL,JK-KL)) / PDZZ(:,JJ-1,JK)        &
              ) * PDZY(:,JJ,JK)* 0.25                                     &
             -(  (PY(:,JJ,JK+KL)-PY(:,JJ,JK))     / PDZZ(:,JJ,JK+KL)        &
                +(PY(:,JJ-1,JK+KL)-PY(:,JJ-1,JK)) / PDZZ(:,JJ-1,JK+KL)      &
              ) * PDZY(:,JJ,JK+KL)* 0.25                                   &
            )  / PDYY(:,JJ,JK)
    END DO
  END DO
!
  DO JJ=1+JPHEXT,IJU
    PGY_M_V(:,JJ,KKU)=  ( PY(:,JJ,KKU)-PY(:,JJ-1,KKU)  )  / PDYY(:,JJ,KKU)
    PGY_M_V(:,JJ,KKA)=  -999.
  END DO
!
  PGY_M_V(:,1,:)=PGY_M_V(:,IJU-2*JPHEXT+1,:)
ELSE
!  PGY_M_V = DYM(PY)/PDYY
  PGY_M_V(:,1+JPHEXT:IJU,:) = ( PY(:,1+JPHEXT:IJU,:)-PY(:,JPHEXT:IJU-1,:) ) &
                               / PDYY(:,1+JPHEXT:IJU,:)
!
  PGY_M_V(:,1,:)=PGY_M_V(:,IJU-2*JPHEXT+1,:)
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GY_M_V',1,ZHOOK_HANDLE)
END FUNCTION GY_M_V
