MODULE MODI_INI_PHYEX
INTERFACE
SUBROUTINE INI_PHYEX(HPROGRAM, KUNITNML, LDNEEDNAM, KLUOUT, KFROM, KTO, &
                    &PTSTEP, PDZMIN, &
                    &CMICRO, &
                    &LDCHANGEMODEL, LDDEFAULTVAL, LDREADNAM, LDCHECK, LDPRINT, LDINIT, &
                    &CST_INOUT, &
                    &PARAM_ICE_INOUT, RAIN_ICE_DESCR_INOUT, RAIN_ICE_PARAM_INOUT, CLOUDPARN_INOUT)
!
USE MODD_CST, ONLY: CST_t
USE MODD_PARAM_ICE, ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE MODD_CLOUDPAR_N,     ONLY: CLOUDPAR_t
!
CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM     !< Current program
INTEGER,           INTENT(IN) :: KUNITNML     !< Logical unit to access the namelist
LOGICAL,           INTENT(IN) :: LDNEEDNAM    !< True to abort if namelist is absent
INTEGER,           INTENT(IN) :: KLUOUT       !< Logical unit for outputs
INTEGER,           INTENT(IN) :: KFROM        !< Old model number
INTEGER,           INTENT(IN) :: KTO          !< New model number
REAL,              INTENT(IN) :: PTSTEP       !< Timestep
REAL,              INTENT(IN) :: PDZMIN       !< Minimum thickness
CHARACTER(4),      INTENT(IN) :: CMICRO       !< Microphysical scheme to use
LOGICAL, OPTIONAL, INTENT(IN) :: LDCHANGEMODEL!< Must we change the active model
LOGICAL, OPTIONAL, INTENT(IN) :: LDDEFAULTVAL !< Must we initialize variables with default values (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDREADNAM    !< Must we read the namelist (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDCHECK      !< Must we perform some checks on values (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDPRINT      !< Must we print the effective values (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDINIT       !< Must we call the init routines
TYPE(CST_t),             OPTIONAL, INTENT(INOUT) :: CST_INOUT
TYPE(PARAM_ICE_t),       OPTIONAL, INTENT(INOUT) :: PARAM_ICE_INOUT
TYPE(RAIN_ICE_DESCR_t) , OPTIONAL, INTENT(INOUT) :: RAIN_ICE_DESCR_INOUT
TYPE(RAIN_ICE_PARAM_t) , OPTIONAL, INTENT(INOUT) :: RAIN_ICE_PARAM_INOUT
TYPE(CLOUDPAR_t),        OPTIONAL, INTENT(INOUT) :: CLOUDPARN_INOUT

END SUBROUTINE INI_PHYEX
END INTERFACE
END MODULE MODI_INI_PHYEX
