!     ######spl
      MODULE MODI_CH_READ_VECTOR
!!    ##########################
!!
INTERFACE
!!
SUBROUTINE CH_READ_VECTOR(KEQ, HNAMES, PVAR, PDEFAULT, KIN, KOUT, KVERB)
IMPLICIT NONE
!!
INTEGER,                       INTENT(IN) :: KEQ
                        ! number of variables to be defined
CHARACTER*(*), DIMENSION(KEQ), INTENT(IN) :: HNAMES
                        ! names of the variables to be defined
REAL,          DIMENSION(KEQ), INTENT(OUT):: PVAR
                        ! value of the variable to be read
REAL,                          INTENT(IN) :: PDEFAULT
                        ! default value
INTEGER,                       INTENT(IN) :: KIN
                        ! I/O channel for file input
INTEGER,                       INTENT(IN) :: KOUT
                        ! I/O channel for printing
INTEGER,                       INTENT(IN) :: KVERB
                        ! verbosity level
END SUBROUTINE CH_READ_VECTOR
!!
END INTERFACE
!!
END MODULE MODI_CH_READ_VECTOR
