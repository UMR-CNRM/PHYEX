!     ######spl
 MODULE MODI_BUDGET
!##################
!
INTERFACE
!
SUBROUTINE BUDGET(PVARS,KBUDN,HBUVAR,YDDDH, YDLDDH, YDMDDH)
!
USE DDH_MIX, ONLY : TYP_DDH
USE YOMLDDH, ONLY : TLDDH
USE YOMMDDH, ONLY : TMDDH
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARS    ! Source 
INTEGER               , INTENT(IN) :: KBUDN    ! variable number
CHARACTER (LEN=*)    , INTENT(IN) :: HBUVAR   ! Identifier of the Budget of the
                                               ! variable that is considered 
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
!
END SUBROUTINE BUDGET
!
END INTERFACE
!
END MODULE MODI_BUDGET
