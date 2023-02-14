MODULE MODI_INI_PHYEX
INTERFACE
SUBROUTINE INI_PHYEX(HPROGRAM, KUNITNML, LDNEEDNAM, KLUOUT, KFROM, KTO, &
                    &PTSTEP, PDZMIN, &
                    &CMICRO, CSCONV, CTURB, &
                    &LDCHANGEMODEL, LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT, LDINIT, &
                    &CST_IN, CST_OUT, &
                    &PARAM_ICE_IN, PARAM_ICE_OUT, RAIN_ICE_DESCR_IN, RAIN_ICE_DESCR_OUT, RAIN_ICE_PARAM_IN, RAIN_ICE_PARAM_OUT, &
                    &CLOUDPARN_IN, CLOUDPARN_OUT, &
                    &PARAM_MFSHALLN_IN, PARAM_MFSHALLN_OUT, &
                    &TURBN_IN, TURBN_OUT, CSTURB_IN, CSTURB_OUT)
!
USE MODD_CST, ONLY: CST_t
USE MODD_PARAM_ICE, ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE MODD_CLOUDPAR_N,     ONLY: CLOUDPAR_t
USE MODD_PARAM_MFSHALL_N,ONLY: PARAM_MFSHALL_t
USE MODD_TURB_N,         ONLY: TURB_t
USE MODD_CTURB,          ONLY: CSTURB_t
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
TYPE(CST_t),             OPTIONAL, INTENT(IN)  :: CST_IN
TYPE(CST_t),             OPTIONAL, INTENT(OUT) :: CST_OUT
TYPE(PARAM_ICE_t),       OPTIONAL, INTENT(IN)  :: PARAM_ICE_IN
TYPE(PARAM_ICE_t),       OPTIONAL, INTENT(OUT) :: PARAM_ICE_OUT
TYPE(RAIN_ICE_DESCR_t) , OPTIONAL, INTENT(IN)  :: RAIN_ICE_DESCR_IN
TYPE(RAIN_ICE_DESCR_t) , OPTIONAL, INTENT(OUT) :: RAIN_ICE_DESCR_OUT
TYPE(RAIN_ICE_PARAM_t) , OPTIONAL, INTENT(IN)  :: RAIN_ICE_PARAM_IN
TYPE(RAIN_ICE_PARAM_t) , OPTIONAL, INTENT(OUT) :: RAIN_ICE_PARAM_OUT
TYPE(CLOUDPAR_t),        OPTIONAL, INTENT(IN)  :: CLOUDPARN_IN
TYPE(CLOUDPAR_t),        OPTIONAL, INTENT(OUT) :: CLOUDPARN_OUT
TYPE(PARAM_MFSHALL_t),   OPTIONAL, INTENT(IN)  :: PARAM_MFSHALLN_IN
TYPE(PARAM_MFSHALL_t),   OPTIONAL, INTENT(OUT) :: PARAM_MFSHALLN_OUT
TYPE(TURB_t),            OPTIONAL, INTENT(IN)  :: TURBN_IN
TYPE(TURB_t),            OPTIONAL, INTENT(OUT) :: TURBN_OUT
TYPE(CSTURB_t),          OPTIONAL, INTENT(IN)  :: CSTURB_IN
TYPE(CSTURB_t),          OPTIONAL, INTENT(OUT) :: CSTURB_OUT



END SUBROUTINE INI_PHYEX
END INTERFACE
END MODULE MODI_INI_PHYEX
