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
TYPE(CST_t),             OPTIONAL, INTENT(IN)  :: CST_IN               !< Structure for constants (IN)
TYPE(CST_t),             OPTIONAL, INTENT(INOUT) :: CST_OUT            !< Structure for constants (OUT)
TYPE(PARAM_ICE_t),       OPTIONAL, INTENT(IN)  :: PARAM_ICE_IN         !< Structure for controling ICE3/ICE4 (IN)
TYPE(PARAM_ICE_t),       OPTIONAL, INTENT(INOUT) :: PARAM_ICE_OUT      !< Structure for controling ICE3/ICE4 (OUT)
TYPE(RAIN_ICE_DESCR_t) , OPTIONAL, INTENT(IN)  :: RAIN_ICE_DESCR_IN    !< Structure for describing hydrometeors (IN)
TYPE(RAIN_ICE_DESCR_t) , OPTIONAL, INTENT(INOUT) :: RAIN_ICE_DESCR_OUT !< Structure for describing hydrometeors (OUT)
TYPE(RAIN_ICE_PARAM_t) , OPTIONAL, INTENT(IN)  :: RAIN_ICE_PARAM_IN    !< Structure for ICE3/ICE4 precomputed values (IN)
TYPE(RAIN_ICE_PARAM_t) , OPTIONAL, INTENT(INOUT) :: RAIN_ICE_PARAM_OUT !< Structure for ICE3/ICE4 precomputed values (OUT)
TYPE(CLOUDPAR_t),        OPTIONAL, INTENT(IN)  :: CLOUDPARN_IN         !< Structure for model dependant microphysics variables (IN)
TYPE(CLOUDPAR_t),        OPTIONAL, INTENT(INOUT) :: CLOUDPARN_OUT      !< Structure for model dependant microphysics variables (IN
TYPE(PARAM_MFSHALL_t),   OPTIONAL, INTENT(IN)  :: PARAM_MFSHALLN_IN    !< Structure for controling shallow convection scheme (IN)
TYPE(PARAM_MFSHALL_t),   OPTIONAL, INTENT(INOUT) :: PARAM_MFSHALLN_OUT !< Structure for controling shallow convection scheme (OUT)
TYPE(TURB_t),            OPTIONAL, INTENT(IN)  :: TURBN_IN             !< Structure for controling the turbulence scheme (IN)
TYPE(TURB_t),            OPTIONAL, INTENT(INOUT) :: TURBN_OUT          !< Structure for controling the turbulence scheme (IN)
TYPE(CSTURB_t),          OPTIONAL, INTENT(IN)  :: CSTURB_IN            !< Structure for the turbulence scheme constants (IN)
TYPE(CSTURB_t),          OPTIONAL, INTENT(INOUT) :: CSTURB_OUT         !< Structure for the turbulence scheme constants (IN)

!IMPORTANT NOTE on *_OUT arguments.
!Logically those arguments should be declared with INTENT(OUT) but in this case ifort (at least) breaks the
!execution when same the structure is given for the _IN and the _OUT argument.
!When INITENT(INOUT) is used, execution is OK on ifort.




END SUBROUTINE INI_PHYEX
END INTERFACE
END MODULE MODI_INI_PHYEX
