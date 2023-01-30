MODULE MODE_POSNAM_PHY
CONTAINS
SUBROUTINE POSNAM_PHY(KUNITNML, CDNAML, LDNEEDNAM, LDFOUND, KLUOUT)
!Wrapper to call the AROME version of posnam

IMPLICIT NONE

INTEGER,          INTENT(IN)    :: KUNITNML  !< Logical unit to access the namelist
CHARACTER(LEN=*), INTENT(IN)    :: CDNAML    !< Namelist name
LOGICAL,          INTENT(IN)    :: LDNEEDNAM !< True to abort if namelist is absent
LOGICAL,          INTENT(OUT)   :: LDFOUND   !< True if namelist has been found
INTEGER,          INTENT(IN)    :: KLUOUT    !< Logical unit for output

#include "posnam.intfb.h"
CALL POSNAM(KUNITNML, CDNAML)
LDFOUND=.TRUE. !Posnam aborts if not found

END SUBROUTINE POSNAM_PHY
END MODULE MODE_POSNAM_PHY
