!     ######spl
      MODULE MODD_CH_M9
!!    #################
!!
!! This code is the MESONH interface to constant and variables defined in module
!! MODD_CH_M9_SCHEME that is proper to one chemical scheme. This interface should
!! reduce the source dependances and then improve the compilation time.
!!
!!*** *MODD_CH_M9*
!!
!!    PURPOSE
!!    -------
!     definition of variables and constant for the chemical core system
!!
!!**  METHOD
!!    ------
!!    The constants NEQ and NREAC are duplicated here in order to avoid
!!    decouple the CCS from the other modules of MNHC.
!!
!!    BEWARE : you must call the procedure 'CH_INIT_SCHEME' before using any
!!             variables from this module.
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Didier Gazen (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 19/10/03
!!
!!----------------------------------------------------------------------
!!    DECLARATIONS
!!    ------------
IMPLICIT NONE
!
INTEGER :: NEQ           ! number of prognostic chemical species
INTEGER :: NREAC         ! number of chemical reactions
INTEGER :: NMETEOVARS    ! number of meteorological variables
INTEGER :: NNONZEROTERMS ! number of non-zero terms returned by CH_TERMS
!
CHARACTER(LEN=32),  DIMENSION(:), POINTER  :: CNAMES=>NULL() ! names of the species
CHARACTER(LEN=32),  DIMENSION(:), POINTER  :: CREACS=>NULL() ! the reaction rate names
CHARACTER(LEN=256), DIMENSION(:), POINTER  :: CFULLREACS=>NULL() ! the full reactions
!
TYPE METEOTRANSTYPE ! variables from the meteorological part
  REAL,              DIMENSION(20) :: XMETEOVAR  ! the meteorological variables
  CHARACTER(LEN=32), DIMENSION(20) :: CMETEOVAR  ! their names
END TYPE METEOTRANSTYPE
!
END MODULE MODD_CH_M9
