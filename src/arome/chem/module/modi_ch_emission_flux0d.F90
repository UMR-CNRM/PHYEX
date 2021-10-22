!     ######spl
      MODULE MODI_CH_EMISSION_FLUX0D
!!    ############################
!!    
INTERFACE
SUBROUTINE CH_EMISSION_FLUX0D(PTIME, PFLUX, HINPUTFILE, KLUOUT, KVERB)
USE MODD_CH_M9,      ONLY: NEQ
IMPLICIT NONE
REAL,                 INTENT(IN)  :: PTIME      ! time of simulation in sec UTC
                                                ! (counting from midnight)
REAL, DIMENSION(NEQ), INTENT(OUT) :: PFLUX      ! emission flux in ppp*m/s
CHARACTER*(*),        INTENT(IN)  :: HINPUTFILE ! name of the input file
INTEGER,              INTENT(IN)  :: KLUOUT     ! output listing channel
INTEGER,              INTENT(IN)  :: KVERB      ! verbosity level
END SUBROUTINE CH_EMISSION_FLUX0D
END INTERFACE
!!
END MODULE MODI_CH_EMISSION_FLUX0D
