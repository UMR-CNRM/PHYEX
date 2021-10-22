!     ######spl
      FUNCTION GZ_M_W(KKA,KKU,KL,PY,PDZZ) RESULT(PGZ_M_W)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #########################################
!
!!****  *GZ_M_W * - Compute the gradient along z direction for a
!!       variable localized at a mass point
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute a gradient along x,y,z
!     directions for a field PY localized at a mass point. The result PGZ_M_W
!     is localized at a z-flux point (w point)
!
!
!                    dzm(PY)
!       PGZ_M_W =    -------
!                     d*zz
!
!!**  METHOD
!!    ------
!!      We employ the Shuman operators to compute the derivatives and the
!!    averages. The metric coefficients PDZZ are dummy arguments.
!!
!!
!!    EXTERNAL
!!    --------
!!      FUNCTION DZM : compute a finite difference along the z
!!    direction for a variable at a mass localization
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODI_SHUMAN : interface for the Shuman functions
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (function GZ_M_W)
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
!!                         19/07/00  inlining(J. Stein)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_SHUMAN
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!              -------------------------------------
!
INTEGER,           INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,           INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise

                                                          ! Metric coefficient:
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGZ_M_W  ! result at flux
                                                              ! side
!
INTEGER :: IKT,IKTB,IKTE
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG Z
!              -----------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GZ_M_W',0,ZHOOK_HANDLE)
IKT=SIZE(PY,3)
IKTB=1+JPVEXT_TURB
IKTE=IKT-JPVEXT_TURB

PGZ_M_W(:,:,IKTB:IKTE) =  (PY(:,:,IKTB:IKTE)-PY(:,:,IKTB-KL:IKTE-KL))  &
                           / PDZZ(:,:,IKTB:IKTE)
PGZ_M_W(:,:,KKU)=  (PY(:,:,KKU)-PY(:,:,KKU-KL))  &
                           / PDZZ(:,:,KKU)
PGZ_M_W(:,:,KKA)=-999.
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GZ_M_W',1,ZHOOK_HANDLE)
END FUNCTION GZ_M_W
