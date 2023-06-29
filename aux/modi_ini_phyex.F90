MODULE MODI_INI_PHYEX
IMPLICIT NONE
INTERFACE
SUBROUTINE INI_PHYEX(HPROGRAM, KUNITNML, LDNEEDNAM, KLUOUT, KFROM, KTO, &
                    &PTSTEP, PDZMIN, &
                    &CMICRO, CSCONV, CTURB, &
                    &LDCHANGEMODEL, LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT, LDINIT, &
                    &PHYEX_IN, PHYEX_OUT)
!
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODD_CST, ONLY: CST_t
USE MODD_PARAM_ICE_n, ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_CLOUDPAR_N,     ONLY: CLOUDPAR_t
USE MODD_PARAM_MFSHALL_N,ONLY: PARAM_MFSHALL_t
USE MODD_TURB_N,         ONLY: TURB_t
USE MODD_CTURB,          ONLY: CSTURB_t
USE MODD_NEB_N,          ONLY: NEB_t
!
IMPLICIT NONE

CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM     !< Current program
INTEGER,           INTENT(IN) :: KUNITNML     !< Logical unit to access the namelist
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
