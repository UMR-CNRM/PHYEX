!     ######spl
      FUNCTION GX_M_U(KKA,KKU,KL,PY,PDXX,PDZZ,PDZX) RESULT(PGX_M_U)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##################################################
!
!!****  *GX_M_U * - Compute the gradient along x for a variable localized at
!!                  a mass point
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute a gradient along x
!     direction for a field PY localized at a mass point. The result PGX_M_U
!     is localized at a x-flux point (u point).
!
!                    (           ____________z )
!                    (               ________x )
!                 1  (                dzm(PY)  )
!   PGX_M_U =   ---- (dxm(PY) - d*zx --------  )
!               d*xx (                 d*zz    )
!
!
!
!!**  METHOD
!!    ------
!!      We employ the Shuman operators to compute the derivatives and the
!!    averages. The metric coefficients PDXX,PDZX,PDZZ are dummy arguments.
!!
!!
!!    EXTERNAL
!!    --------
!!      FUNCTION DXM: compute a finite difference along the x direction for
!!      a variable at a mass localization
!!      FUNCTION DZM: compute a finite difference along the y direction for
!!      a variable at a mass localization
!!      FUNCTION MXM: compute an average in the x direction for a variable
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
!!      Book2 of documentation (function GX_M_U)
!!
!!
!!    AUTHOR
!!    ------
!!      P. Hereil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original           05/07/94
!!      Modification       16/03/95  change the order of the arguments
!!                         19/07/00  add the LFLAT switch  + inlining(J. Stein)
!!                         20/08/00  optimization (J. Escobar)
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
!*       0.1   Declarations of arguments and result
!              ------------------------------------
!
INTEGER,                INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                INTENT(IN)  :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX                   ! d*xx
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZX                   ! d*zx
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   ! d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGX_M_U  ! result at flux
                                                              ! side
INTEGER  IIU,IKU,JI,JK
!
INTEGER :: JJK,IJU
INTEGER :: JIJK,JIJKOR,JIJKEND
INTEGER :: JI_1JK, JIJK_1, JI_1JK_1, JIJKP1, JI_1JKP1
!
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG X
!              -----------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GX_M_U',0,ZHOOK_HANDLE)
IIU=SIZE(PY,1)
IJU=SIZE(PY,2)
IKU=SIZE(PY,3)
IF (.NOT. LFLAT) THEN
! PGX_M_U = (   DXM(PY)  -  MZF (   MXM(  DZM(PY) /PDZZ  ) * PDZX   )   )/PDXX
!!  DO JK=1+JPVEXT_TURB,IKU-JPVEXT_TURB
!!    DO JI=1+JPHEXT,IIU
!!        PGX_M_U(JI,:,JK)=                                                 &
!!           (  PY(JI,:,JK)-PY(JI-1,:,JK)                                   &
!!             -(  (PY(JI,:,JK)-PY(JI,:,JK-1))     / PDZZ(JI,:,JK)          &
!!                +(PY(JI-1,:,JK)-PY(JI-1,:,JK-1)) / PDZZ(JI-1,:,JK)        &
!!              ) * PDZX(JI,:,JK)* 0.25                                     &
!!             -(  (PY(JI,:,JK+1)-PY(JI,:,JK))     / PDZZ(JI,:,JK+1)        &
!!                +(PY(JI-1,:,JK+1)-PY(JI-1,:,JK)) / PDZZ(JI-1,:,JK+1)      &
!!              ) * PDZX(JI,:,JK+1)* 0.25                                   &
!!            )  / PDXX(JI,:,JK)
!!    END DO
!!  END DO
  JIJKOR  = 1 + JPHEXT + IIU*IJU*(JPVEXT_TURB+1 - 1)
  JIJKEND = IIU*IJU*(IKU-JPVEXT_TURB)
!CDIR NODEP
!OCL NOVREC
  DO JIJK=JIJKOR , JIJKEND
! indexation
    JI_1JK   = JIJK - 1
    JIJK_1   = JIJK     - IIU*IJU*KL
    JI_1JK_1 = JIJK - 1 - IIU*IJU*KL
    JIJKP1   = JIJK     + IIU*IJU*KL
    JI_1JKP1 = JIJK - 1 + IIU*IJU*KL
!
    PGX_M_U(JIJK,1,1)=                                              &
       (  PY(JIJK,1,1)-PY(JI_1JK,1,1)                               &
       -(  (PY(JIJK,1,1)-PY(JIJK_1,1,1))     / PDZZ(JIJK,1,1)       &
       +(PY(JI_1JK,1,1)-PY(JI_1JK_1,1,1)) / PDZZ(JI_1JK,1,1)        &
       ) * PDZX(JIJK,1,1)* 0.25                                     &
       -(  (PY(JIJKP1,1,1)-PY(JIJK,1,1))     / PDZZ(JIJKP1,1,1)     &
       +(PY(JI_1JKP1,1,1)-PY(JI_1JK,1,1)) / PDZZ(JI_1JKP1,1,1)      &
       ) * PDZX(JIJKP1,1,1)* 0.25                                   &
        )  / PDXX(JIJK,1,1)
  END DO

!
  DO JI=1+JPHEXT,IIU
    PGX_M_U(JI,:,KKU)=  ( PY(JI,:,KKU)-PY(JI-1,:,KKU)  )  / PDXX(JI,:,KKU)
    PGX_M_U(JI,:,KKA)=  -999.
  END DO
!
  PGX_M_U(1,:,:)=PGX_M_U(IIU-2*JPHEXT+1,:,:)
ELSE
!  PGX_M_U = DXM(PY) / PDXX
  PGX_M_U(1+JPHEXT:IIU,:,:) = ( PY(1+JPHEXT:IIU,:,:)-PY(JPHEXT:IIU-1,:,:) ) &
                             / PDXX(1+JPHEXT:IIU,:,:)
!
  PGX_M_U(1,:,:)=PGX_M_U(IIU-2*JPHEXT+1,:,:)
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GX_M_U',1,ZHOOK_HANDLE)
END FUNCTION GX_M_U
