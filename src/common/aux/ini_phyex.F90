SUBROUTINE INI_PHYEX(HPROGRAM, KUNITNML, LDNEEDNAM, KLUOUT, KFROM, KTO, &
                    &PTSTEP, PDZMIN, &
                    &CMICRO, CSCONV, CTURB, &
                    &LDCHANGEMODEL, LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT, LDINIT, &
                    &CST_IN, CST_OUT, &
                    &PARAM_ICEN_IN, PARAM_ICEN_OUT, &
                    &RAIN_ICE_DESCRN_IN, RAIN_ICE_DESCRN_OUT, RAIN_ICE_PARAMN_IN, RAIN_ICE_PARAMN_OUT, &
                    &CLOUDPARN_IN, CLOUDPARN_OUT, &
                    &PARAM_MFSHALLN_IN, PARAM_MFSHALLN_OUT, &
                    &TURBN_IN, TURBN_OUT, CSTURB_IN, CSTURB_OUT)
!
USE MODD_CST, ONLY: CST, CST_t, PRINT_CST
USE MODD_PARAM_ICE_n, ONLY: PARAM_ICE_GOTO_MODEL, PARAM_ICEN_INIT, PARAM_ICEN, PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_GOTO_MODEL, RAIN_ICE_DESCRN, RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_GOTO_MODEL, RAIN_ICE_PARAMN, RAIN_ICE_PARAM_t
USE MODD_CLOUDPAR_N,     ONLY: CLOUDPAR_GOTO_MODEL, CLOUDPARN, CLOUDPAR_t
USE MODD_PARAM_MFSHALL_N,ONLY: PARAM_MFSHALLN, PARAM_MFSHALL_t, PARAM_MFSHALLN_INIT, PARAM_MFSHALL_GOTO_MODEL
USE MODD_TURB_N,         ONLY: TURBN, TURB_t, TURBN_INIT, TURB_GOTO_MODEL
USE MODD_CTURB,          ONLY: CSTURB, CSTURB_t, CTURB_ASSOCIATE
!
USE MODE_INI_CST, ONLY: INI_CST
USE MODE_INI_RAIN_ICE, ONLY: INI_RAIN_ICE
USE MODE_INI_TIWMX, ONLY: INI_TIWMX
USE MODE_INI_SNOW, ONLY: INI_SNOW
USE MODE_INI_MFSHALL, ONLY: INI_MFSHALL
USE MODE_INI_TURB, ONLY: INI_TURB
!
!!
!!      *INI_PHYEX* - PHYEX initialisation routine
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialise (default values, namelist, checks, prints)
!!      the different modules used by the physics
!!
!!
!!    AUTHOR
!!    ------
!!     S. Riette
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    -  Original      Fev 2023
!!
!-------------------------------------------------------------------------------
!
!**       DECLARATIONS
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
TYPE(CST_t),             OPTIONAL, INTENT(IN)    :: CST_IN              !< Structure for constants (IN)
TYPE(CST_t),             OPTIONAL, INTENT(INOUT) :: CST_OUT             !< Structure for constants (OUT)
TYPE(PARAM_ICE_t),       OPTIONAL, INTENT(IN)    :: PARAM_ICEN_IN       !< Structure for controling ICE3/ICE4 (IN)
TYPE(PARAM_ICE_t),       OPTIONAL, INTENT(INOUT) :: PARAM_ICEN_OUT      !< Structure for controling ICE3/ICE4 (OUT)
TYPE(RAIN_ICE_DESCR_t) , OPTIONAL, INTENT(IN)    :: RAIN_ICE_DESCRN_IN  !< Structure for describing hydrometeors (IN)
TYPE(RAIN_ICE_DESCR_t) , OPTIONAL, INTENT(INOUT) :: RAIN_ICE_DESCRN_OUT !< Structure for describing hydrometeors (OUT)
TYPE(RAIN_ICE_PARAM_t) , OPTIONAL, INTENT(IN)    :: RAIN_ICE_PARAMN_IN  !< Structure for ICE3/ICE4 precomputed values (IN)
TYPE(RAIN_ICE_PARAM_t) , OPTIONAL, INTENT(INOUT) :: RAIN_ICE_PARAMN_OUT !< Structure for ICE3/ICE4 precomputed values (OUT)
TYPE(CLOUDPAR_t),        OPTIONAL, INTENT(IN)    :: CLOUDPARN_IN        !< Structure for model dependant microphysics variables (IN)
TYPE(CLOUDPAR_t),        OPTIONAL, INTENT(INOUT) :: CLOUDPARN_OUT       !< Structure for model dependant microphysics variables (IN
TYPE(PARAM_MFSHALL_t),   OPTIONAL, INTENT(IN)    :: PARAM_MFSHALLN_IN   !< Structure for controling shallow convection scheme (IN)
TYPE(PARAM_MFSHALL_t),   OPTIONAL, INTENT(INOUT) :: PARAM_MFSHALLN_OUT  !< Structure for controling shallow convection scheme (OUT)
TYPE(TURB_t),            OPTIONAL, INTENT(IN)    :: TURBN_IN            !< Structure for controling the turbulence scheme (IN)
TYPE(TURB_t),            OPTIONAL, INTENT(INOUT) :: TURBN_OUT           !< Structure for controling the turbulence scheme (IN)
TYPE(CSTURB_t),          OPTIONAL, INTENT(IN)    :: CSTURB_IN           !< Structure for the turbulence scheme constants (IN)
TYPE(CSTURB_t),          OPTIONAL, INTENT(INOUT) :: CSTURB_OUT          !< Structure for the turbulence scheme constants (IN)

!IMPORTANT NOTE on *_OUT arguments.
!Logically those arguments should be declared with INTENT(OUT) but in this case ifort (at least) breaks the
!execution when same the structure is given for the _IN and the _OUT argument.
!When INITENT(INOUT) is used, execution is OK on ifort.

LOGICAL :: LLINIT, LLCHANGEMODEL
INTEGER :: IPRINT
!
!**       ARGUMENTS
!
LLINIT=.TRUE.
IF(PRESENT(LDINIT)) LLINIT=LDINIT
LLCHANGEMODEL=.TRUE.
IF(PRESENT(LDCHANGEMODEL)) LLCHANGEMODEL=LDCHANGEMODEL
IPRINT=0
IF(PRESENT(KPRINT)) IPRINT=KPRINT
!
!**       CST
!
IF(LLINIT) THEN
  IF(IPRINT==2) WRITE(UNIT=KLUOUT,FMT='('' MODD_CST '')')
  IF(PRESENT(CST_IN)) CST=CST_IN
  CALL INI_CST()
  IF(IPRINT==2) CALL PRINT_CST(KLUOUT)
  IF(PRESENT(CST_OUT)) CST_OUT=CST
ENDIF
!
!**       MICROPHYSICS SCHEME
!
IF(CMICRO=='ICE3' .OR. CMICRO=='ICE4' .OR. CMICRO=='LIMA') THEN
  IF(IPRINT==2) WRITE(UNIT=KLUOUT,FMT='('' MODD_PARAM_ICEN, MODD_RAIN_ICE_DESCRN, MODD_RAIN_ICE_PARAMN, MODD_CLOUDPARN '')')
  IF(LLCHANGEMODEL) THEN
    CALL CLOUDPAR_GOTO_MODEL(KFROM, KTO)
    CALL PARAM_ICE_GOTO_MODEL(KFROM, KTO)
    CALL RAIN_ICE_DESCR_GOTO_MODEL(KFROM, KTO)
    CALL RAIN_ICE_PARAM_GOTO_MODEL(KFROM, KTO)
  ENDIF
  IF(PRESENT(PARAM_ICEN_IN)) PARAM_ICEN=PARAM_ICEN_IN
  IF(PRESENT(RAIN_ICE_DESCRN_IN)) RAIN_ICE_DESCRN=RAIN_ICE_DESCRN_IN
  IF(PRESENT(RAIN_ICE_PARAMN_IN)) RAIN_ICE_PARAMN=RAIN_ICE_PARAMN_IN
  IF(PRESENT(CLOUDPARN_IN)) CLOUDPARN=CLOUDPARN_IN

  CALL PARAM_ICEN_INIT(HPROGRAM, KUNITNML, LDNEEDNAM, KLUOUT, &
                      &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
  IF(LLINIT) THEN
    CALL INI_RAIN_ICE(KLUOUT, PTSTEP, PDZMIN, CLOUDPARN%NSPLITR, CMICRO)
    CALL INI_TIWMX
  
    IF(RAIN_ICE_PARAMN%XFRMIN(16) > 0.) THEN
       CALL INI_SNOW(KLUOUT) ! Recalculate snow parameters :  XCCS = XFRMIN(16),XCXS = XFRMIN(17)
    ENDIF
  ENDIF

  IF(PRESENT(PARAM_ICEN_OUT)) PARAM_ICEN_OUT=PARAM_ICEN
  IF(PRESENT(RAIN_ICE_DESCRN_OUT)) RAIN_ICE_DESCRN_OUT=RAIN_ICE_DESCRN
  IF(PRESENT(RAIN_ICE_PARAMN_OUT)) RAIN_ICE_PARAMN_OUT=RAIN_ICE_PARAMN
  IF(PRESENT(CLOUDPARN_OUT)) CLOUDPARN_OUT=CLOUDPARN
ENDIF
!
!**       SHALLOW CONVECTION SCHEME
!
IF(CSCONV=='EDKF') THEN
  IF(IPRINT==2) WRITE(UNIT=KLUOUT,FMT='('' MODD_PARAM_MFSHALL_n '')')
  IF(LLCHANGEMODEL) CALL PARAM_MFSHALL_GOTO_MODEL(KFROM, KTO)
  IF(PRESENT(PARAM_MFSHALLN_IN)) PARAM_MFSHALLN=PARAM_MFSHALLN_IN

  CALL PARAM_MFSHALLN_INIT(HPROGRAM, KUNITNML, LDNEEDNAM, KLUOUT, &
                          &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
  IF(LLINIT) THEN
    CALL INI_MFSHALL()
  ENDIF

  IF(PRESENT(PARAM_MFSHALLN_OUT)) PARAM_MFSHALLN_OUT=PARAM_MFSHALLN
ENDIF
!
!**       TURBULENCE SCHEME
!
IF(CTURB=='TKEL') THEN
  IF(IPRINT==2) WRITE(UNIT=KLUOUT,FMT='('' MODD_TURB_n MODD_CTURB '')')
  IF(LLCHANGEMODEL) CALL TURB_GOTO_MODEL(KFROM, KTO)
  IF(PRESENT(TURBN_IN)) TURBN=TURBN_IN
  IF(PRESENT(CSTURB_IN)) CSTURB=CSTURB_IN

  CALL CTURB_ASSOCIATE()
  CALL TURBN_INIT(HPROGRAM, KUNITNML, LDNEEDNAM, KLUOUT, &
                 &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
  IF(LLINIT) THEN
    CALL INI_TURB(HPROGRAM)
  ENDIF

  IF(PRESENT(TURBN_OUT)) TURBN_OUT=TURBN
  IF(PRESENT(CSTURB_OUT)) CSTURB_OUT=CSTURB
ENDIF

END SUBROUTINE INI_PHYEX
