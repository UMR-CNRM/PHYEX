!     ######spl
      MODULE MODI_CH_READ_CHEM
!!    ######################## 
!!
!
INTERFACE
SUBROUTINE CH_READ_CHEM(PCONC, PAERO, PRHODREF, HFILE)
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(INOUT) :: PCONC ! gas concentration vector to be read
REAL, DIMENSION(:), INTENT(INOUT) :: PAERO ! aerosol concentration vector to be read
REAL, DIMENSION(:), INTENT(IN) :: PRHODREF ! air density
CHARACTER(LEN=*), INTENT(IN)      :: HFILE ! name of the file to be read from
END SUBROUTINE CH_READ_CHEM
END INTERFACE
!!
END MODULE MODI_CH_READ_CHEM
