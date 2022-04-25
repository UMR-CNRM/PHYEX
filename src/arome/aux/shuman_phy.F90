MODULE SHUMAN_PHY
IMPLICIT NONE
CONTAINS
!     ###############################
      SUBROUTINE MZM_PHY(D,PA,PMZM)
!     ###############################
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *MZM* -  Shuman operator : mean operator in z direction for a
!!                                 mass variable
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean
!     along the z direction (K index) for a field PA localized at a mass
!     point. The result is localized at a z-flux point (w point).
!
!!**  METHOD
!!    ------
!!        The result PMZM(:,:,k) is defined by 0.5*(PA(:,:,k)+PA(:,:,k-1))
!!        At k=1, PMZM(:,:,1) is defined by -999.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PMZM   ! result at flux localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK             ! Loop index in z direction
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZM
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MZM',0,ZHOOK_HANDLE)
DO JK=2,SIZE(PA,3)-1
  PMZM(:,:,JK) = 0.5*( PA(:,:,JK)+PA(:,:,JK-D%NKL) )
END DO
PMZM(:,:,D%NKA)    = -999.
PMZM(:,:,D%NKU) = 0.5*( PA(:,:,D%NKU)+PA(:,:,D%NKU-D%NKL) )
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MZM',1,ZHOOK_HANDLE)
END SUBROUTINE MZM_PHY
END MODULE SHUMAN_PHY
