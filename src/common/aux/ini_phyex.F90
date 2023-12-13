SUBROUTINE INI_PHYEX(HPROGRAM, TPFILE, LDNEEDNAM, KLUOUT, KFROM, KTO, &
                    &PTSTEP, PDZMIN, &
                    &CMICRO, CSCONV, CTURB, &
                    &LDCHANGEMODEL, LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT, LDINIT, &
                    &PHYEX_IN, PHYEX_OUT)
!
!IMPORTANT NOTE ON HOW TO DEAL WITH KEYS USED IN THE PHYSICS *AND* IN THE HOST MODEL
!
!When this situation occurs:
!For example, if a physics option (eg logical switch) has an impact on the model setup, this
!option must be accessible in the host model *and* in the physics.
!
!The prefered solution (when it is possible):
!The option is implemented in a module and in a namelist of PHYEX. And, after the initialisation
!is done, the option is available in the PHYEX_OUT structure (in addition to the PHYEX module) for use
!by the host model.
!
!When the physics option is needed in the host model before the PHYEX initialisation, or if the physics
!option must be computed from other options known by the host model, there are several solutions:
!
!  If this situation is specific to only one host model:
!  The option is declared (module) and initialised (namelist) twice: once in the host model, once in PHYEX;
!  and a consistency check is added after the physics initialisation.
!  This solution is preferred to overwriting the variable stored in the PHYEX module, because between
!  initialization and overwriting, the variable may have been used to initialize a parameterization.
!
!  If this issue is common to all the host models:
!  The option is declared and initialised in the host model. Then, the value is given, through INI_PHYEX, to
!  the parametrisation INIT subroutine (eg TURBN_INIT). In the INIT subroutine the value is assigned
!  to a module variable specific to PHYEX (in order to be accessible from inside the parametrisations).
!  In this case the variable *must not* be added in a PHYEX namelist.
!  An alternative to this solution is not to copy the variable inside a PHYEX module but add it directly
!  to the parametrisation call arguments.
!
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODD_CST, ONLY: CST, PRINT_CST
USE MODD_PARAM_ICE_n, ONLY: PARAM_ICE_GOTO_MODEL, PARAM_ICEN_INIT, PARAM_ICEN, &
                          & CSUBG_AUCV_RC, CSUBG_AUCV_RI
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_GOTO_MODEL, RAIN_ICE_DESCRN
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_GOTO_MODEL, RAIN_ICE_PARAMN
USE MODD_CLOUDPAR_N,     ONLY: CLOUDPAR_GOTO_MODEL, CLOUDPARN
USE MODD_PARAM_MFSHALL_N,ONLY: PARAM_MFSHALLN, PARAM_MFSHALLN_INIT, PARAM_MFSHALL_GOTO_MODEL
USE MODD_TURB_N,         ONLY: TURBN, TURBN_INIT, TURB_GOTO_MODEL, LHARAT
USE MODD_CTURB,          ONLY: CSTURB, CTURB_ASSOCIATE
USE MODD_NEB_N,          ONLY: NEBN, NEBN_INIT, NEB_GOTO_MODEL, CCONDENS, LSTATNW, LSUBG_COND
USE MODD_PARAM_LIMA,     ONLY: PARAM_LIMA, PARAM_LIMA_INIT, PARAM_LIMA_ASSOCIATE, LPTSPLIT, LADJ
USE MODD_PARAM_LIMA_WARM, ONLY: PARAM_LIMA_WARM
USE MODD_PARAM_LIMA_COLD, ONLY: PARAM_LIMA_COLD
USE MODD_PARAM_LIMA_MIXED, ONLY: PARAM_LIMA_MIXED
USE MODD_NSV,            ONLY: TNSV, NSV_ASSOCIATE
USE MODD_IO,             ONLY: TFILEDATA
!
USE MODE_INI_CST, ONLY: INI_CST
USE MODE_INI_RAIN_ICE, ONLY: INI_RAIN_ICE
USE MODE_INI_TIWMX, ONLY: INI_TIWMX
USE MODE_INI_SNOW, ONLY: INI_SNOW
USE MODE_INI_MFSHALL, ONLY: INI_MFSHALL
USE MODE_INI_TURB, ONLY: INI_TURB
USE MODE_LIMA_UPDATE_NSV, ONLY: LIMA_UPDATE_NSV
USE MODE_INIT_AEROSOL_PROPERTIES, ONLY: INIT_AEROSOL_PROPERTIES
USE MODE_INI_LIMA, ONLY: INI_LIMA
!
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
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
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE       !< Namelist file
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

LOGICAL :: LLINIT, LLCHANGEMODEL, LLCHECK
INTEGER :: IPRINT
!
!**       ARGUMENTS
!
LLINIT=.TRUE.
IF(PRESENT(LDINIT)) LLINIT=LDINIT
LLCHECK=.TRUE.
IF(PRESENT(LDCHECK)) LLCHECK=LDCHECK
LLCHANGEMODEL=.TRUE.
IF(PRESENT(LDCHANGEMODEL)) LLCHANGEMODEL=LDCHANGEMODEL
IPRINT=0
IF(PRESENT(KPRINT)) IPRINT=KPRINT
!
!**       CST
!
IF(LLINIT) THEN
  IF(IPRINT==2) WRITE(UNIT=KLUOUT,FMT='('' MODD_CST '')')
  IF(PRESENT(PHYEX_IN)) CST=PHYEX_IN%CST
  CALL INI_CST()
  IF(IPRINT==2) CALL PRINT_CST(KLUOUT)
  IF(PRESENT(PHYEX_OUT)) PHYEX_OUT%CST=CST
ENDIF
!
!**       MICROPHYSICS SCHEME
!
IF(CMICRO=='ICE3' .OR. CMICRO=='ICE4' .OR. CMICRO=='LIMA') THEN
  !The LIMA scheme makes use of the condensation routine of the ICE3/ICE4 scheme
  !It is why the ICE3/ICE4 initialisation is needed in for the LIMA scheme
  IF(IPRINT==2) WRITE(UNIT=KLUOUT,FMT='('' MODD_PARAM_ICEN, MODD_RAIN_ICE_DESCRN, MODD_RAIN_ICE_PARAMN, MODD_CLOUDPARN '')')
  IF(LLCHANGEMODEL) THEN
    CALL CLOUDPAR_GOTO_MODEL(KFROM, KTO)
    CALL PARAM_ICE_GOTO_MODEL(KFROM, KTO)
    CALL RAIN_ICE_DESCR_GOTO_MODEL(KFROM, KTO)
    CALL RAIN_ICE_PARAM_GOTO_MODEL(KFROM, KTO)
    IF(CMICRO=='LIMA') THEN
      CALL PARAM_LIMA_ASSOCIATE()
    ENDIF
  ENDIF
  IF(PRESENT(PHYEX_IN)) THEN
    PARAM_ICEN=PHYEX_IN%PARAM_ICEN
    RAIN_ICE_DESCRN=PHYEX_IN%RAIN_ICE_DESCRN
    RAIN_ICE_PARAMN=PHYEX_IN%RAIN_ICE_PARAMN
    CLOUDPARN=PHYEX_IN%CLOUDPARN
    IF(CMICRO=='LIMA') THEN
      PARAM_LIMA=PHYEX_IN%PARAM_LIMA
      PARAM_LIMA_WARM=PHYEX_IN%PARAM_LIMA_WARM
      PARAM_LIMA_COLD=PHYEX_IN%PARAM_LIMA_COLD
      PARAM_LIMA_MIXED=PHYEX_IN%PARAM_LIMA_MIXED
    ENDIF
  ENDIF

  CALL PARAM_ICEN_INIT(HPROGRAM, TPFILE, LDNEEDNAM, KLUOUT, &
                      &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
  IF(CMICRO=='LIMA') THEN
    CALL PARAM_LIMA_INIT(HPROGRAM, TPFILE, LDNEEDNAM, KLUOUT, &
                        &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
  ENDIF

  IF(LLINIT) THEN
    CALL INI_RAIN_ICE(KLUOUT, PTSTEP, PDZMIN, CLOUDPARN%NSPLITR, CMICRO)
    CALL INI_TIWMX
  
    IF(RAIN_ICE_PARAMN%XFRMIN(16) > 0.) THEN
       CALL INI_SNOW(KLUOUT) ! Recalculate snow parameters :  XCCS = XFRMIN(16),XCXS = XFRMIN(17)
    ENDIF

    IF(CMICRO=='LIMA') THEN
      CALL INIT_AEROSOL_PROPERTIES
      CALL INI_LIMA(PTSTEP, PDZMIN, CLOUDPARN%NSPLITR, CLOUDPARN%NSPLITG)
    ENDIF
  ENDIF

  IF(PRESENT(PHYEX_OUT)) THEN
    PHYEX_OUT%PARAM_ICEN=PARAM_ICEN
    PHYEX_OUT%RAIN_ICE_DESCRN=RAIN_ICE_DESCRN
    PHYEX_OUT%RAIN_ICE_PARAMN=RAIN_ICE_PARAMN
    PHYEX_OUT%CLOUDPARN=CLOUDPARN
    IF(CMICRO=='LIMA') THEN
      PHYEX_OUT%PARAM_LIMA=PARAM_LIMA
      PHYEX_OUT%PARAM_LIMA_WARM=PARAM_LIMA_WARM
      PHYEX_OUT%PARAM_LIMA_COLD=PARAM_LIMA_COLD
      PHYEX_OUT%PARAM_LIMA_MIXED=PARAM_LIMA_MIXED
    ENDIF
  ENDIF
ENDIF
!
!**       NSV init (must be after LIMA init)
!
IF(LLINIT) THEN
  IF(LLCHANGEMODEL) CALL NSV_ASSOCIATE()
  IF(PRESENT(PHYEX_IN)) TNSV=PHYEX_IN%TNSV
  TNSV%NSV=0
  CALL LIMA_UPDATE_NSV(LDINIT=.TRUE., KMI=KTO, KSV=TNSV%NSV, &
                       &CDCLOUD=CMICRO, LDUPDATE=.TRUE.)
  TNSV%NSV=TNSV%NSV_LIMA
  IF(PRESENT(PHYEX_OUT)) PHYEX_OUT%TNSV=TNSV
ENDIF
!
!**       SHALLOW CONVECTION SCHEME
!
IF(CSCONV=='EDKF') THEN
  IF(IPRINT==2) WRITE(UNIT=KLUOUT,FMT='('' MODD_PARAM_MFSHALL_n '')')
  IF(LLCHANGEMODEL) CALL PARAM_MFSHALL_GOTO_MODEL(KFROM, KTO)
  IF(PRESENT(PHYEX_IN)) PARAM_MFSHALLN=PHYEX_IN%PARAM_MFSHALLN

  CALL PARAM_MFSHALLN_INIT(HPROGRAM, TPFILE, LDNEEDNAM, KLUOUT, &
                          &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
  IF(LLINIT) THEN
    CALL INI_MFSHALL()
  ENDIF

  IF(PRESENT(PHYEX_OUT)) PHYEX_OUT%PARAM_MFSHALLN=PARAM_MFSHALLN
ENDIF
!
!**       CLOUD SCHEME
!
IF(.TRUE.) THEN !Placeholder for configuration without cloud scheme or a different one
  IF(IPRINT==2) WRITE(UNIT=KLUOUT,FMT='('' MODD_NEB_n '')')
  IF(LLCHANGEMODEL) CALL NEB_GOTO_MODEL(KFROM, KTO)
  IF(PRESENT(PHYEX_IN)) NEBN=PHYEX_IN%NEBN

  CALL NEBN_INIT(HPROGRAM, TPFILE, LDNEEDNAM, KLUOUT, &
                &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
  IF(LLINIT) THEN
    !Nothing to do, everything is read from namelist
  ENDIF

  IF(PRESENT(PHYEX_OUT)) PHYEX_OUT%NEBN=NEBN
ENDIF
!
!**       TURBULENCE SCHEME
!
IF(CTURB=='TKEL') THEN
  IF(IPRINT==2) WRITE(UNIT=KLUOUT,FMT='('' MODD_TURB_n MODD_CTURB '')')
  IF(LLCHANGEMODEL) THEN
    CALL TURB_GOTO_MODEL(KFROM, KTO)
    CALL CTURB_ASSOCIATE()
  ENDIF
  IF(PRESENT(PHYEX_IN)) THEN
    TURBN=PHYEX_IN%TURBN
    CSTURB=PHYEX_IN%CSTURB
  ENDIF

  CALL TURBN_INIT(HPROGRAM, TPFILE, LDNEEDNAM, KLUOUT, &
                 &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
  IF(LLINIT) THEN
    CALL INI_TURB(HPROGRAM)
  ENDIF

  IF(PRESENT(PHYEX_OUT)) THEN
    PHYEX_OUT%TURBN=TURBN
    PHYEX_OUT%CSTURB=CSTURB
  ENDIF
ENDIF
!
!**       GLOBAL CONSISTENCY TESTS
!
IF(LLCHECK) THEN
  IF((CMICRO=='ICE3' .OR. CMICRO=='ICE4' .OR. CMICRO=='LIMA') .AND. CTURB=='TKEL') THEN
    IF ((CSUBG_AUCV_RC == 'ADJU' .OR. CSUBG_AUCV_RI == 'ADJU') .AND. CCONDENS /= 'GAUS') THEN
      CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'INI_PHYEX', &
                    &"CSUBG_AUCV_RC and/or CSUBG_AUCV_RI cannot be 'ADJU' if CCONDENS is not 'GAUS'")
    ENDIF
    IF (.NOT. LHARAT .AND. LSTATNW) THEN
      CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'INI_PHYEX', &
                    &'LSTATNW only tested in combination with HARATU and EDMFm!')
    ENDIF
  ENDIF
  !
  IF(CMICRO=='LIMA') THEN
   IF (LSUBG_COND .AND. (.NOT. LPTSPLIT)) THEN
      CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'INI_PHYEX', &                                                                             
                    &"YOU MUST USE LPTSPLIT=T WITH CMICRO=LIMA AND LSUBG_COND")
    END IF
    IF (LSUBG_COND .AND. (.NOT. LADJ)) THEN
      CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'INI_PHYEX', &
                    &"YOU MUST USE LADJ=t WITH CMICRO=LIMA AND LSUBG_COND")
    END IF
  ENDIF
ENDIF

END SUBROUTINE INI_PHYEX
