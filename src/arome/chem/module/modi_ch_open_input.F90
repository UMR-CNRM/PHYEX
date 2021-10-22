!     ######spl
      MODULE MODI_CH_OPEN_INPUT
!!    ######################### 
!!
INTERFACE
!!
SUBROUTINE CH_OPEN_INPUT(HCHEM_INPUT_FILE,HKEYWORD,KCHANNEL,KLUOUT,KVERB)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: HCHEM_INPUT_FILE ! general purpose input file
CHARACTER(LEN=*), INTENT(IN) :: HKEYWORD         ! keyword for positioning
INTEGER         , INTENT(OUT):: KCHANNEL         ! I/O channel to choose
INTEGER,          INTENT(IN) :: KLUOUT           ! output listing channel
INTEGER,          INTENT(IN) :: KVERB            ! verbosity level
END SUBROUTINE CH_OPEN_INPUT
!!
END INTERFACE
!!
END MODULE MODI_CH_OPEN_INPUT
