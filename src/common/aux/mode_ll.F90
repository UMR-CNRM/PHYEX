MODULE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODE_TOOLS
IMPLICIT NONE
CONTAINS
  SUBROUTINE GET_INDICE_ll(KXOR, KYOR, KXEND, KYEND, KSIZE1, KSIZE2)
  USE MODD_PARAMETERS, ONLY : JPHEXT
  IMPLICIT NONE
  INTEGER, INTENT(IN),OPTIONAL :: KSIZE1, KSIZE2
  INTEGER, INTENT(OUT) :: KXOR, KYOR, KXEND, KYEND
  KXOR=1+JPHEXT
  KYOR=1+JPHEXT
  KXEND=KSIZE1-JPHEXT
  KYEND=KSIZE2-JPHEXT
  END SUBROUTINE GET_INDICE_ll

  SUBROUTINE UPDATE_HALO_ll(TPLIST, KINFO)
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll  
  TYPE(LIST_ll), POINTER :: TPLIST ! pointer to the list of fields to be updated
  INTEGER                :: KINFO  ! return status
  CALL ABORT
  END SUBROUTINE UPDATE_HALO_ll

  SUBROUTINE GET_DIM_EXT_ll(CBORD,IIU,IJU)
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN) :: CBORD
  INTEGER, INTENT(IN) :: IIU,IJU
  END SUBROUTINE GET_DIM_EXT_ll
LOGICAL FUNCTION LNORTH_ll()
  LNORTH_ll=.FALSE.
END FUNCTION LNORTH_ll
!
LOGICAL FUNCTION LEAST_ll()
  LEAST_ll=.FALSE.
END FUNCTION LEAST_ll
!
LOGICAL FUNCTION LWEST_ll()
  LWEST_ll=.FALSE.
END FUNCTION LWEST_ll
!
LOGICAL FUNCTION LSOUTH_ll()
  LSOUTH_ll=.FALSE.
END FUNCTION LSOUTH_ll
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

END MODULE MODE_ll
