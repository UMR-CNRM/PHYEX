!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
      MODULE MODD_PARAM_ICE_n
!     #####################
!> @file
!!      *MODD_PARAM_ICE_n* - declaration of the control parameters for the
!!                           mixed phase cloud parameterization
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to define the set of space
!!    and time control parameters for the microphysics.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_PARAM_ICE)
!!
!!    AUTHOR
!!    ------
!!     J.-P. Pinty   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    -  Original      14/12/95
!!    -  Jan 2015 S. Riette: new ICE3/ICE4 parameters
!!    -  01/10/16 (C.Lac)  Add droplet deposition for fog
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE
!
TYPE PARAM_ICE_t
LOGICAL :: LWARM       !< When .TRUE. activates the formation of rain by
                       !! the warm microphysical processes
LOGICAL :: LSEDIC      !< TRUE to enable the droplet sedimentation
LOGICAL :: LDEPOSC     !< TRUE to enable cloud droplet deposition
REAL    :: XVDEPOSC    !< Droplet deposition velocity
!
CHARACTER(LEN=4) :: CPRISTINE_ICE !< Pristine ice type PLAT, COLU or BURO
CHARACTER(LEN=4) :: CSEDIM        !< Sedimentation calculation mode      
!
LOGICAL :: LRED       !< To use modified ICE3/ICE4 to reduce time step dependency
LOGICAL :: LFEEDBACKT !< When .TRUE. feed back on temperature is taken into account
LOGICAL :: LEVLIMIT   !< When .TRUE. water vapour pressure is limited by saturation
LOGICAL :: LNULLWETG  !< When .TRUE. graupel wet growth is activated with null rate (to allow water shedding)
LOGICAL :: LWETGPOST  !< When .TRUE. graupel wet growth is activated with positive temperature (to allow water shedding)
LOGICAL :: LNULLWETH  !< Same as LNULLWETG but for hail
LOGICAL :: LWETHPOST  !< Same as LWETGPOST but for hail
CHARACTER(LEN=4) :: CSNOWRIMING !< OLD or M90 for Murakami 1990 formulation
REAL :: XFRACM90      !< Fraction used for the Murakami 1990 formulation
INTEGER :: NMAXITER_MICRO   !< Maximum number of iterations for mixing ratio or time splitting
REAL :: XMRSTEP       !< maximum mixing ratio step for mixing ratio splitting
LOGICAL :: LCONVHG    !< TRUE to allow the conversion from hail to graupel
LOGICAL :: LCRFLIMIT  !< True to limit rain contact freezing to possible heat exchange
!
REAL :: XTSTEP_TS     !< Approximative time step for time-splitting (0 for no time-splitting)
!
CHARACTER(LEN=80) :: CSUBG_RC_RR_ACCR !< subgrid rc-rr accretion
CHARACTER(LEN=80) :: CSUBG_RR_EVAP    !< subgrid rr evaporation
CHARACTER(LEN=80) :: CSUBG_PR_PDF     !< pdf for subgrid precipitation
CHARACTER(LEN=4)  :: CSUBG_AUCV_RC    !< type of subgrid rc->rr autoconv. method
CHARACTER(LEN=80) :: CSUBG_AUCV_RI    !< type of subgrid ri->rs autoconv. method
CHARACTER(LEN=80) :: CSUBG_MF_PDF     !< PDF to use for MF cloud autoconversions
!
LOGICAL :: LADJ_BEFORE !< must we perform an adjustment before rain_ice call
LOGICAL :: LADJ_AFTER  !< must we perform an adjustment after rain_ice call
LOGICAL :: LSEDIM_AFTER !< sedimentation done before (.FALSE.) or after (.TRUE.) microphysics
!
REAL :: XSPLIT_MAXCFL   !< Maximum CFL number allowed for SPLIT scheme
LOGICAL :: LSNOW_T      !< Snow parameterization from Wurtz (2021)
!
LOGICAL :: LPACK_INTERP !< To pack arrays before computing the different interpolations (kernels and other)
LOGICAL :: LPACK_MICRO  !< To pack arrays before computing the process tendencies
!
INTEGER :: NPROMICRO    !< Size of cache-blocking bloc (0 to deactivate)
!
LOGICAL :: LCRIAUTI     !< .T. to compute XACRIAUTI and XBCRIAUTI (from XCRIAUTI and XT0CRIAUTI);
                        !! .F. to compute XT0CRIAUTI (from XCRIAUTI and XBCRIAUTI)
REAL :: XCRIAUTI_NAM    !< Minimum value for the ice->snow autoconversion threshold
REAL :: XT0CRIAUTI_NAM  !< Threshold temperature (???) for the ice->snow autoconversion threshold
REAL :: XBCRIAUTI_NAM   !< B barameter for the ice->snow autoconversion 10**(aT+b) law
REAL :: XACRIAUTI_NAM   !< A barameter for the ice->snow autoconversion 10**(aT+b) law
REAL :: XCRIAUTC_NAM    !< Threshold for liquid cloud -> rain autoconversion (kg/m**3)
REAL :: XRDEPSRED_NAM   !< Tuning factor of sublimation of snow
REAL :: XRDEPGRED_NAM   !< Tuning factor of sublimation of graupel
!
LOGICAL :: LOCND2       !< Logical switch to separate liquid and ice
REAL, DIMENSION(40) :: XFRMIN_NAM
!
END TYPE PARAM_ICE_t
!
TYPE(PARAM_ICE_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: PARAM_ICE_MODEL
TYPE(PARAM_ICE_t), POINTER, SAVE :: PARAM_ICEN => NULL()
!
LOGICAL, POINTER :: LWARM => NULL(), &
                    LSEDIC => NULL(), &
                    LDEPOSC => NULL(), &
                    LRED => NULL(), &
                    LFEEDBACKT => NULL(), &
                    LEVLIMIT => NULL(), &
                    LNULLWETG => NULL(), &
                    LWETGPOST => NULL(), &
                    LNULLWETH => NULL(), &
                    LWETHPOST => NULL(), &
                    LCONVHG => NULL(), &
                    LCRFLIMIT => NULL(), &
                    LADJ_BEFORE => NULL(), &
                    LADJ_AFTER => NULL(), &
                    LSEDIM_AFTER => NULL(), &
                    LSNOW_T => NULL(), &
                    LPACK_INTERP => NULL(), &
                    LPACK_MICRO => NULL(), &
                    LCRIAUTI => NULL(), &
                    LOCND2 => NULL()

REAL, POINTER :: XVDEPOSC => NULL(), &
                 XFRACM90 => NULL(), &
                 XMRSTEP => NULL(), &
                 XTSTEP_TS => NULL(), &
                 XSPLIT_MAXCFL => NULL(), &
                 XCRIAUTI_NAM => NULL(), &
                 XT0CRIAUTI_NAM => NULL(), &
                 XBCRIAUTI_NAM => NULL(), &
                 XACRIAUTI_NAM => NULL(), &
                 XCRIAUTC_NAM => NULL(), &
                 XRDEPSRED_NAM => NULL(), &
                 XRDEPGRED_NAM => NULL()
REAL, DIMENSION(:), POINTER :: XFRMIN_NAM => NULL()

INTEGER, POINTER :: NMAXITER_MICRO => NULL(), &
                    NPROMICRO => NULL()

CHARACTER(LEN=4), POINTER :: CPRISTINE_ICE => NULL()
CHARACTER(LEN=4), POINTER :: CSEDIM => NULL()
CHARACTER(LEN=4), POINTER :: CSNOWRIMING => NULL()
CHARACTER(LEN=80),POINTER :: CSUBG_RC_RR_ACCR => NULL()
CHARACTER(LEN=80),POINTER :: CSUBG_RR_EVAP => NULL()
CHARACTER(LEN=80),POINTER :: CSUBG_PR_PDF => NULL()
CHARACTER(LEN=4),POINTER :: CSUBG_AUCV_RC=>NULL()
CHARACTER(LEN=80),POINTER :: CSUBG_AUCV_RI=>NULL()
CHARACTER(LEN=80),POINTER :: CSUBG_MF_PDF=>NULL()
!
NAMELIST/NAM_PARAM_ICEn/LWARM,LSEDIC,LCONVHG,CPRISTINE_ICE,CSEDIM,LDEPOSC,XVDEPOSC, &
                       LRED, LFEEDBACKT, &
                       LEVLIMIT,LNULLWETG,LWETGPOST,LNULLWETH,LWETHPOST, &
                       CSNOWRIMING,XFRACM90,NMAXITER_MICRO,XMRSTEP,XTSTEP_TS, &
                       LADJ_BEFORE, LADJ_AFTER, LCRFLIMIT, &
                       XSPLIT_MAXCFL, LSEDIM_AFTER, LSNOW_T, &
                       LPACK_INTERP, LPACK_MICRO, NPROMICRO, CSUBG_RC_RR_ACCR, &
                       CSUBG_RR_EVAP, CSUBG_PR_PDF, CSUBG_AUCV_RC, CSUBG_AUCV_RI, &
                       LCRIAUTI, XCRIAUTI_NAM, XT0CRIAUTI_NAM, XBCRIAUTI_NAM, &
                       XACRIAUTI_NAM, XCRIAUTC_NAM, XRDEPSRED_NAM, XRDEPGRED_NAM, &
                       LOCND2, XFRMIN_NAM, CSUBG_MF_PDF
!
!-------------------------------------------------------------------------------
!
CONTAINS
SUBROUTINE PARAM_ICE_GOTO_MODEL(KFROM, KTO)
!! This subroutine associate all the pointers to the right component of
!! the right strucuture. A value can be accessed through the structure PARAM_ICEN
!! or through the strucuture PARAM_ICE_MODEL(KTO) or directly through these pointers.
IMPLICIT NONE
INTEGER, INTENT(IN) :: KFROM, KTO
!
IF(.NOT. ASSOCIATED(PARAM_ICEN, PARAM_ICE_MODEL(KTO))) THEN
  !
  PARAM_ICEN => PARAM_ICE_MODEL(KTO)
  !
  LWARM => PARAM_ICEN%LWARM
  LSEDIC => PARAM_ICEN%LSEDIC
  LDEPOSC => PARAM_ICEN%LDEPOSC
  LRED => PARAM_ICEN%LRED
  LFEEDBACKT => PARAM_ICEN%LFEEDBACKT
  LEVLIMIT => PARAM_ICEN%LEVLIMIT
  LNULLWETG => PARAM_ICEN%LNULLWETG
  LWETGPOST => PARAM_ICEN%LWETGPOST
  LNULLWETH => PARAM_ICEN%LNULLWETH
  LWETHPOST => PARAM_ICEN%LWETHPOST
  LCONVHG => PARAM_ICEN%LCONVHG
  LCRFLIMIT => PARAM_ICEN%LCRFLIMIT
  LADJ_BEFORE => PARAM_ICEN%LADJ_BEFORE
  LADJ_AFTER => PARAM_ICEN%LADJ_AFTER
  LSEDIM_AFTER => PARAM_ICEN%LSEDIM_AFTER
  LSNOW_T => PARAM_ICEN%LSNOW_T
  LPACK_INTERP => PARAM_ICEN%LPACK_INTERP
  LPACK_MICRO => PARAM_ICEN%LPACK_MICRO
  LCRIAUTI => PARAM_ICEN%LCRIAUTI
  LOCND2 => PARAM_ICEN%LOCND2
  !
  XVDEPOSC => PARAM_ICEN%XVDEPOSC
  XFRACM90 => PARAM_ICEN%XFRACM90
  XMRSTEP => PARAM_ICEN%XMRSTEP
  XTSTEP_TS => PARAM_ICEN%XTSTEP_TS
  XSPLIT_MAXCFL => PARAM_ICEN%XSPLIT_MAXCFL
  XCRIAUTI_NAM => PARAM_ICEN%XCRIAUTI_NAM
  XT0CRIAUTI_NAM => PARAM_ICEN%XT0CRIAUTI_NAM
  XBCRIAUTI_NAM => PARAM_ICEN%XBCRIAUTI_NAM
  XACRIAUTI_NAM => PARAM_ICEN%XACRIAUTI_NAM
  XCRIAUTC_NAM => PARAM_ICEN%XCRIAUTC_NAM
  XRDEPSRED_NAM => PARAM_ICEN%XRDEPSRED_NAM
  XRDEPGRED_NAM => PARAM_ICEN%XRDEPGRED_NAM
  XFRMIN_NAM => PARAM_ICEN%XFRMIN_NAM
  !
  NMAXITER_MICRO => PARAM_ICEN%NMAXITER_MICRO
  NPROMICRO => PARAM_ICEN%NPROMICRO
  !
  CPRISTINE_ICE => PARAM_ICEN%CPRISTINE_ICE
  CSEDIM => PARAM_ICEN%CSEDIM
  CSNOWRIMING => PARAM_ICEN%CSNOWRIMING
  CSUBG_RC_RR_ACCR => PARAM_ICEN%CSUBG_RC_RR_ACCR
  CSUBG_RR_EVAP => PARAM_ICEN%CSUBG_RR_EVAP
  CSUBG_PR_PDF => PARAM_ICEN%CSUBG_PR_PDF
  CSUBG_AUCV_RC=>PARAM_ICEN%CSUBG_AUCV_RC
  CSUBG_AUCV_RI=>PARAM_ICEN%CSUBG_AUCV_RI
  CSUBG_MF_PDF=>PARAM_ICEN%CSUBG_MF_PDF
ENDIF
END SUBROUTINE PARAM_ICE_GOTO_MODEL
!
SUBROUTINE PARAM_ICEN_INIT(HPROGRAM, KUNITNML, LDNEEDNAM, KLUOUT, &
                          &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
!!*** *PARAM_ICEN_INIT* - Code needed to initialize the MODD_PARAM_ICE_n module
!!
!!*   PURPOSE
!!    -------
!!    Sets the default values, reads the namelist, performs the checks and prints
!!
!!*   METHOD
!!    ------
!!    0. Declarations
!!       1. Declaration of arguments
!!       2. Declaration of local variables
!!    1. Default values
!!    2. Namelist
!!    3. Checks
!!    4. Prints
!!
!!    AUTHOR
!!    ------
!!    S. Riette
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    Feb 2023
!-------------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!       ---------------
!
USE MODE_POSNAM_PHY, ONLY: POSNAM_PHY
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
USE MODE_CHECK_NAM_VAL, ONLY: CHECK_NAM_VAL_CHAR, CHECK_NAM_VAL_REAL, CHECK_NAM_VAL_INT
!
IMPLICIT NONE
!
!* 0.1. Declaration of arguments
!       ------------------------
!
CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM     !< Name of the calling program
INTEGER,           INTENT(IN) :: KUNITNML     !< Logical unit to access the namelist
LOGICAL,           INTENT(IN) :: LDNEEDNAM    !< True to abort if namelist is absent
INTEGER,           INTENT(IN) :: KLUOUT       !< Logical unit for outputs
LOGICAL, OPTIONAL, INTENT(IN) :: LDDEFAULTVAL !< Must we initialize variables with default values (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDREADNAM    !< Must we read the namelist (defaults to .TRUE.)
LOGICAL, OPTIONAL, INTENT(IN) :: LDCHECK      !< Must we perform some checks on values (defaults to .TRUE.)
INTEGER, OPTIONAL, INTENT(IN) :: KPRINT       !< Print level (defaults to 0): 0 for no print, 1 to safely print namelist,
                                              !! 2 to print informative messages
!
!* 0.2 Declaration of local variables
!      ------------------------------
!
LOGICAL :: LLDEFAULTVAL, LLREADNAM, LLCHECK, LLFOUND
INTEGER :: IPRINT

LLDEFAULTVAL=.TRUE.
LLREADNAM=.TRUE.
LLCHECK=.TRUE.
IPRINT=0
IF(PRESENT(LDDEFAULTVAL)) LLDEFAULTVAL=LDDEFAULTVAL
IF(PRESENT(LDREADNAM   )) LLREADNAM   =LDREADNAM
IF(PRESENT(LDCHECK     )) LLCHECK     =LDCHECK
IF(PRESENT(KPRINT      )) IPRINT      =KPRINT
!
!*      1. DEFAULT VALUES
!       -----------------
!
IF(LLDEFAULTVAL) THEN
  !NOTES ON GENERAL DEFAULTS AND MODEL-SPECIFIC DEFAULTS :
  !- General default values *MUST* remain unchanged.
  !- To change the default value for a given application,                                                 
  !  an "IF(HPROGRAM=='...')" condition must be used.

  LWARM=.TRUE.
  LSEDIC=.TRUE.
  LDEPOSC=.FALSE.
  XVDEPOSC= 0.02 ! 2 cm/s
  CPRISTINE_ICE='PLAT'
  CSEDIM='SPLI'
  LRED=.TRUE.
  LFEEDBACKT=.TRUE.
  LEVLIMIT=.TRUE.
  LNULLWETG=.TRUE.
  LWETGPOST=.TRUE.
  LNULLWETH=.TRUE.
  LWETHPOST=.TRUE.
  CSNOWRIMING='M90'
  XFRACM90=0.1
  NMAXITER_MICRO=5
  XMRSTEP=0.00005
  LCONVHG=.FALSE.
  LCRFLIMIT=.TRUE.
  XTSTEP_TS=0.
  CSUBG_RC_RR_ACCR='NONE'
  CSUBG_RR_EVAP='NONE'
  CSUBG_PR_PDF='SIGM'
  CSUBG_AUCV_RC='NONE'
  CSUBG_AUCV_RI='NONE'
  CSUBG_MF_PDF='TRIANGLE'
  LADJ_BEFORE=.TRUE.
  LADJ_AFTER=.TRUE.
  LSEDIM_AFTER=.FALSE.
  XSPLIT_MAXCFL=0.8
  LSNOW_T=.FALSE.
  LPACK_INTERP=.TRUE.
  LPACK_MICRO=.TRUE.
  NPROMICRO=0
  LCRIAUTI=.FALSE.
  !!XCRIAUTIi_NAM = 0.25E-3 !  Critical ice content for the autoconversion to occur
  XCRIAUTI_NAM = 0.2E-4 !  Revised value by Chaboureau et al. (2001)
  XACRIAUTI_NAM=0.06
  XBCRIAUTI_NAM=-3.5
  XT0CRIAUTI_NAM=(LOG10(XCRIAUTI_NAM)-XBCRIAUTI_NAM)/0.06
  XCRIAUTC_NAM=0.5E-3
  XRDEPSRED_NAM=1.
  XRDEPGRED_NAM=1.
  LOCND2=.FALSE.
  ! Tuning and modication of graupeln etc:
  XFRMIN_NAM(1:6)=0.
  XFRMIN_NAM(7:9)=1.
  XFRMIN_NAM(10) =10.
  XFRMIN_NAM(11) =1.
  XFRMIN_NAM(12) =0.
  XFRMIN_NAM(13) =1.0E-15
  XFRMIN_NAM(14) =120.
  XFRMIN_NAM(15) =1.0E-4
  XFRMIN_NAM(16:20)=0.
  XFRMIN_NAM(21:22)=1.
  XFRMIN_NAM(23)=0.5
  XFRMIN_NAM(24)=1.5
  XFRMIN_NAM(25)=30.
  XFRMIN_NAM(26:38)=0.
  XFRMIN_NAM(39)=0.25
  XFRMIN_NAM(40)=0.15

  IF(HPROGRAM=='AROME') THEN
    LCONVHG=.TRUE.
    LADJ_BEFORE=.TRUE.
    LADJ_AFTER=.FALSE.
    LRED=.FALSE.
    CSEDIM='STAT'
    LSEDIC=.FALSE.
    XMRSTEP=0.
    CSUBG_AUCV_RC='PDF'
  ELSEIF(HPROGRAM=='LMDZ') THEN
    CSUBG_AUCV_RC='PDF'
    CSEDIM='STAT'
    NMAXITER_MICRO=1
    LCRIAUTI=.TRUE.
    XCRIAUTC_NAM=0.001
    XCRIAUTI_NAM=0.0002
    XT0CRIAUTI_NAM=-5.
    LRED=.TRUE.
    LCONVHG=.TRUE.
    LADJ_BEFORE=.TRUE.
    LADJ_AFTER=.FALSE.
  ENDIF
ENDIF
!
!*      2. NAMELIST
!       -----------
!
IF(LLREADNAM) THEN
  CALL POSNAM_PHY(KUNITNML, 'NAM_PARAM_ICEN', LDNEEDNAM, LLFOUND, KLUOUT)
  IF(LLFOUND) READ(UNIT=KUNITNML, NML=NAM_PARAM_ICEn)
ENDIF
!
!*      3. CHECKS
!       ---------
!
IF(LLCHECK) THEN
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CPRISTINE_ICE', CPRISTINE_ICE, 'PLAT', 'COLU', 'BURO')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CSEDIM', CSEDIM, 'SPLI', 'STAT', 'NONE')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CSUBG_RC_RR_ACCR', CSUBG_RC_RR_ACCR, 'NONE', 'PRFR')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CSUBG_RR_EVAP', CSUBG_RR_EVAP, 'NONE', 'CLFR', 'PRFR')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CSUBG_PR_PDF', CSUBG_PR_PDF, 'SIGM', 'HLCRECTPD', 'HLCTRIANGPDF', &
                                                                'HLCQUADRAPDF', 'HLCISOTRIPDF')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CSUBG_AUCV_RC', CSUBG_AUCV_RC, 'PDF ', 'CLFR', 'NONE', 'ADJU', 'SIGM')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CSUBG_AUCV_RI', CSUBG_AUCV_RI, 'NONE', 'CLFR', 'ADJU')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CSUBG_MF_PDF', CSUBG_MF_PDF, 'NONE', 'TRIANGLE')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CSNOWRIMING', CSNOWRIMING, 'OLD ', 'M90 ')
  CALL CHECK_NAM_VAL_REAL(KLUOUT, 'XTSTEP_TS', XTSTEP_TS, '>=', 0.)
  CALL CHECK_NAM_VAL_REAL(KLUOUT, 'XSPLIT_MAXCFL', XSPLIT_MAXCFL, '>', 0., '<=', 1.)
  CALL CHECK_NAM_VAL_REAL(KLUOUT, 'XVDEPOSC', XVDEPOSC, '>=', 0.)
  CALL CHECK_NAM_VAL_REAL(KLUOUT, 'XFRACM90', XFRACM90, '>=', 0., '<=', 1.)
  CALL CHECK_NAM_VAL_REAL(KLUOUT, 'XMRSTEP', XMRSTEP, '>=', 0.)
  CALL CHECK_NAM_VAL_INT(KLUOUT, 'NPROMICRO', NPROMICRO, '>=', 0)

  IF (LOCND2 .AND. (XRDEPSRED_NAM /= 1 .OR. XRDEPGRED_NAM /= 1)) THEN
    CALL ABOR1 ("XRDESRED_NAM and XRDEGRED_NAM must not be activated together with LOCND2")
  ENDIF

  IF(HPROGRAM=='AROME' .OR. HPROGRAM=='LMDZ') THEN
    IF(.NOT. (LADJ_BEFORE .AND. .NOT. LADJ_AFTER)) THEN
      CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODD_PARAM_ICE_n', 'With AROME/LMDZ, LADJ_BEFORE must be .T. and LADJ_AFTER must be .F.')
    ENDIF
  ENDIF
ENDIF
!
!*      3. PRINTS
!       ---------
!
IF(IPRINT>=1) THEN
  WRITE(UNIT=KLUOUT,NML=NAM_PARAM_ICEn)
ENDIF
!
END SUBROUTINE PARAM_ICEN_INIT
!
END MODULE MODD_PARAM_ICE_n
