MODULE MODI_SHUMAN_DEVICE
! This module is temporary existing to handle turb_hor specificities in GPU porting
! Merge in _PHY shuman and gradient should be coded in mode_tridiag_w and mode_turb_hor_uv to leviate discrepancies and remove this module
! MZM_DEVICE has another implementation in mesonh
! It is not used in AROME (turb_3D)
IMPLICIT NONE
CONTAINS
SUBROUTINE MZM_DEVICE(PA,PMZM)
IMPLICIT NONE
!
REAL, DIMENSION(:,:), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(:,:), INTENT(OUT) :: PMZM   ! result at flux localization
!
CALL ABORT
END SUBROUTINE MZM_DEVICE
END MODULE MODI_SHUMAN_DEVICE
