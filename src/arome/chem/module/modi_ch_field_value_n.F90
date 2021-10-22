!     ######spl
      MODULE MODI_CH_FIELD_VALUE_n
!!    ############################
!!
INTERFACE
!!
FUNCTION CH_FIELD_VALUE_n(PZ, HGRIDTYPE, HNAME, &
                          HUNIT, KLUOUT, KVERB)
IMPLICIT NONE
REAL                         :: CH_FIELD_VALUE_n! the function name
REAL,             INTENT(IN) :: PZ              ! x-y-z coo. to initialize
CHARACTER(LEN=*), INTENT(IN) :: HGRIDTYPE       ! grid type passed as x-y-z coo.
CHARACTER(LEN=*), INTENT(IN) :: HNAME           ! name of species to initialize
CHARACTER(LEN=*), INTENT(OUT):: HUNIT           ! unit ("CON" or "MIX")
INTEGER,          INTENT(IN) :: KLUOUT          ! output listing channel
INTEGER,          INTENT(IN) :: KVERB           ! verbosity level
!!
END FUNCTION CH_FIELD_VALUE_n
!!
END INTERFACE
!!
END MODULE MODI_CH_FIELD_VALUE_n
