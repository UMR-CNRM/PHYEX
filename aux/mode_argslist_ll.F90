MODULE MODE_ARGSLIST_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
CONTAINS

!
  SUBROUTINE CLEANLIST_ll(TPLIST)
IMPLICIT NONE
    TYPE(LIST_ll),  POINTER :: TPLIST ! List of fields
    CALL ABORT
  END SUBROUTINE CLEANLIST_ll
!
  SUBROUTINE ADD2DFIELD_ll(TPLIST, PFIELD, HNAME)
IMPLICIT NONE
          
    TYPE(LIST_ll), POINTER         :: TPLIST   ! list of fields
    REAL, DIMENSION(:,:), TARGET :: PFIELD   ! field to be added to the list
  !                                              of fields
    character(len=*), intent(in) :: HNAME ! Name of the field to be added
  !
   CALL ABORT  
END SUBROUTINE ADD2DFIELD_ll
!
  SUBROUTINE ADD3DFIELD_ll(TPLIST, PFIELD, HNAME)
IMPLICIT NONE
          
    TYPE(LIST_ll), POINTER         :: TPLIST   ! list of fields
    REAL, DIMENSION(:,:,:), TARGET :: PFIELD   ! field to be added to the list
  !                                              of fields
    character(len=*), intent(in) :: HNAME ! Name of the field to be added
  !
   CALL ABORT  
END SUBROUTINE ADD3DFIELD_ll
!
  SUBROUTINE ADD4DFIELD_ll(TPLIST, PFIELD, HNAME)
IMPLICIT NONE

    TYPE(LIST_ll), POINTER         :: TPLIST   ! list of fields
    REAL, DIMENSION(:,:,:,:), TARGET :: PFIELD   ! field to be added to the list
  !                                              of fields
    character(len=*), intent(in) :: HNAME ! Name of the field to be added
  !
   CALL ABORT
END SUBROUTINE ADD4DFIELD_ll
END MODULE MODE_ARGSLIST_ll
