!     ######spl
     SUBROUTINE CONVECT_CLOSURE_ADJUST_SHAL( KLON, KLEV, PADJ,                      &
                                             PUMF, PZUMF, PUER, PZUER, PUDR, PZUDR  )
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!    ################################################################################
!
!!**** Uses closure adjustment factor to adjust mass flux and to modify
!!     precipitation efficiency  when necessary. The computations are
!!     similar to routine CONVECT_PRECIP_ADJUST.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to adjust the mass flux using the
!!      factor PADJ computed in CONVECT_CLOSURE
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!
!!
!!    EXTERNAL
!!    --------
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    None
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CLOSURE_ADJUST)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  15/11/96
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAREXT, ONLY : JCVEXB, JCVEXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
REAL, DIMENSION(KLON),      INTENT(IN) :: PADJ     ! mass adjustment factor
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUMF ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUER ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUDR ! initial value of  "
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB, IKE                 ! vert. loop bounds
INTEGER :: JK                       ! vertical loop index
!
!
!-------------------------------------------------------------------------------
!
!*       0.3   Compute loop bounds
!              -------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CONVECT_CLOSURE_ADJUST_SHAL',0,ZHOOK_HANDLE)
IKB  = 1 + JCVEXB
IKE  = KLEV - JCVEXT
!
!
!*       1.     Adjust mass flux by the factor PADJ to converge to
!               specified degree of stabilization
!               ----------------------------------------------------
!
     DO JK = IKB + 1, IKE
          PUMF(:,JK)  = PZUMF(:,JK)   * PADJ(:)
          PUER(:,JK)  = PZUER(:,JK)   * PADJ(:)
          PUDR(:,JK)  = PZUDR(:,JK)   * PADJ(:)
     END DO
!
IF (LHOOK) CALL DR_HOOK('CONVECT_CLOSURE_ADJUST_SHAL',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECT_CLOSURE_ADJUST_SHAL
