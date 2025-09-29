MODULE MODI_INI_PHYEX
IMPLICIT NONE
INTERFACE
SUBROUTINE INI_PHYEX(HPROGRAM, TPFILE, LDNEEDNAM, KLUOUT, KFROM, KTO, &
                    &PTSTEP, PDZMIN, &
                    &CMICRO, CSCONV, CTURB, &
                    &LDCHANGEMODEL, LDDEFAULTVAL, LDREADNAM, LDCHECK,&
                    &KPRINT, LDINIT, &
                    &PHYEX_IN, PHYEX_OUT)
!
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODD_IO,             ONLY: TFILEDATA
!
IMPLICIT NONE

CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM     !< Current program
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE       !< Logical unit to access the namelist
LOGICAL,           INTENT(IN) :: LDNEEDNAM    !< True to abort if namelist is absent
INTEGER,           INTENT(IN) :: KLUOUT       !< Logical unit for outputs
INTEGER,           INTENT(IN) :: KFROM        !< Old model number
INTEGER,           INTENT(IN) :: KTO          !< New model number
REAL,              INTENT(IN) :: PTSTEP       !< Timestep
REAL,              INTENT(IN) :: PDZMIN       !< Minimum thickness
CHARACTER(LEN=4),  INTENT(IN) :: CMICRO       !< Microphysical scheme to use
CHARACTER(LEN=4),  INTENT(IN) :: CTURB        !< Turbulence scheme to use
CHARACTER(LEN=4),  INTENT(IN) :: CSCONV       !< Shallow convection scheme to use
LOGICAL, OPTIONAL, INTENT(IN) :: LDCHANGEMODEL!< Must we change the active model
LOGICAL, OPTIONAL, INTENT(IN) :: LDDEFAULTVAL !< Must we initialize variables with default values (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDREADNAM    !< Must we read the namelist (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDCHECK      !< Must we perform some checks on values (defaults to .TRUE.)
INTEGER, OPTIONAL, INTENT(IN) :: KPRINT       !< Print level (defaults to 0): 0 for no print, 1 to safely print namelist,
                                              !! 2 to print informative messages
LOGICAL, OPTIONAL, INTENT(IN) :: LDINIT       !< Must we call the init routines
TYPE(PHYEX_t), OPTIONAL, INTENT(IN)    :: PHYEX_IN    !< Structure for constants (IN)
TYPE(PHYEX_t), OPTIONAL, INTENT(INOUT) :: PHYEX_OUT   !< Structure for constants (OUT)

!IMPORTANT NOTE on PHYEX_OUT arguments.
!Logically this argument should be declared with INTENT(OUT) but in this case ifort (at least) breaks the
!execution when the same structure is given for the PHYEX_IN and the PHYEX_OUT argument.
!When INITENT(INOUT) is used, execution is OK on ifort.




END SUBROUTINE INI_PHYEX
END INTERFACE
END MODULE MODI_INI_PHYEX
