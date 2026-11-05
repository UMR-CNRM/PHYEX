MODULE MODE_INI_MFSHALL
IMPLICIT NONE
CONTAINS
SUBROUTINE INI_MFSHALL()
!     ###########################################################
!
!!****  *INI_MFSHALL * - initialize the constants necessary for the
!!                       shallow convection scheme.
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the constants used by
!!    the shallow convection scheme.
!!
!!**  METHOD
!!    ------
!!      The constants are initialized to their numerical values.
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_MFSHALL_n, ONLY: LTHETAS_MF, XLAMBDA_MF
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!

!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INI_MFSHALL',0,ZHOOK_HANDLE)
!
IF(LTHETAS_MF) THEN
  XLAMBDA_MF=5.87
ELSE
  XLAMBDA_MF=0.
ENDIF
!
IF (LHOOK) CALL DR_HOOK('INI_MFSHALL',1,ZHOOK_HANDLE)
!
END SUBROUTINE INI_MFSHALL
!
END MODULE MODE_INI_MFSHALL
