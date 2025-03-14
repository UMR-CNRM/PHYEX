MODULE MODE_POSNAM_PHY
IMPLICIT NONE
CONTAINS
SUBROUTINE POSNAM_PHY(TFILENAM, CDNAML, LDNEEDNAM, LDFOUND)
!Wrapper to call the AROME version of posnam

USE MODD_IO, ONLY: TFILEDATA

IMPLICIT NONE

TYPE(TFILEDATA),  INTENT(IN)    :: TFILENAM  !< Namelist file
CHARACTER(LEN=*), INTENT(IN)    :: CDNAML    !< Namelist name
LOGICAL,          INTENT(IN)    :: LDNEEDNAM !< True to abort if namelist is absent
LOGICAL,          INTENT(OUT)   :: LDFOUND   !< True if namelist has been found

#include "posnam.intfb.h"
CALL POSNAM(TFILENAM%NLU, CDNAML)
LDFOUND=.TRUE. !Posnam aborts if not found

END SUBROUTINE POSNAM_PHY
END MODULE MODE_POSNAM_PHY
