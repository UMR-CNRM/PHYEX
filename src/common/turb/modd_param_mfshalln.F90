!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     #############################
      MODULE MODD_PARAM_MFSHALL_n
!     #############################
!> @file
!!      *MODD_PARAM_MFSHALL_n* - Declaration of Mass flux scheme free parameters
!!
!!    PURPOSE
!!    -------
!!    The purpose of this declarative module is to declare the
!!    variables that may be set by namelist for the mass flux scheme
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      
!!          
!!    AUTHOR
!!    ------
!!       S. Malardel, J. Pergaud (Meteo France)      
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/02/07
!!      10/16 R.Honnert Update with AROME
!!      01/2019 R.Honnert add parameters for the reduction of mass-flux surface closure with resolution
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE PARAM_MFSHALL_t

REAL               :: XIMPL_MF     !< degre of implicitness
      
CHARACTER (LEN=4)  :: CMF_UPDRAFT  !< Type of Mass Flux Scheme
                                   !! 'NONE' if no parameterization 
CHARACTER (LEN=4)  :: CMF_CLOUD    !< Type of cloud scheme associated

                                     
LOGICAL       :: LMIXUV      !< True if mixing of momentum
LOGICAL       :: LMF_FLX     !< logical switch for the storage of the mass flux fluxes
REAL          :: XALP_PERT   !< coefficient for the perturbation of
                             !! theta_l and r_t at the first level of the updraft
REAL          ::    XABUO    !< coefficient of the buoyancy term in the w_up equation
REAL          ::    XBENTR   !< coefficient of the entrainment term in the w_up equation
REAL          ::    XBDETR   !< coefficient of the detrainment term in the w_up equation
REAL          ::    XCMF     !< coefficient for the mass flux at the first level of the updraft (closure)
REAL          :: XENTR_MF    !< entrainment constant (m/Pa) = 0.2 (m) 
REAL          :: XCRAD_MF    !< cloud radius in cloudy part
REAL          :: XENTR_DRY   !< coefficient for entrainment in dry part 
REAL          :: XDETR_DRY   !< coefficient for detrainment in dry part
REAL          :: XDETR_LUP   !< coefficient for detrainment in dry part
REAL          :: XKCF_MF     !< coefficient for cloud fraction
REAL          :: XKRC_MF     !< coefficient for convective rc
REAL          :: XTAUSIGMF
REAL          :: XPRES_UV    !< coefficient for pressure term in wind mixing
REAL          :: XALPHA_MF   !< coefficient for cloudy fraction
REAL          :: XSIGMA_MF   !< coefficient for sigma computation
REAL          :: XFRAC_UP_MAX!< maximum Updraft fraction
!
REAL          :: XA1         !< Parameter for Rio et al (2010) formulation for entrainment and detrainment (RHCJ10)
REAL          :: XB          !!
REAL          :: XC          !!
REAL          :: XBETA1      !!
!
REAL          :: XR          !< Parameter for closure assumption of Hourdin et al (2002): aspect ratio of updraft
!
LOGICAL       :: LGZ         !< Grey Zone Surface Closure
REAL          :: XGZ         !< Tuning of the surface initialisation for Grey Zone
!
LOGICAL       :: LTHETAS_MF      !< .TRUE. to use ThetaS1 instead of ThetaL
REAL          :: XLAMBDA_MF      !< Thermodynamic parameter: Lambda to compute ThetaS1 from ThetaL

END TYPE PARAM_MFSHALL_t

TYPE(PARAM_MFSHALL_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: PARAM_MFSHALL_MODEL
TYPE(PARAM_MFSHALL_t), POINTER, SAVE :: PARAM_MFSHALLN => NULL()
                                     
REAL             , POINTER :: XIMPL_MF=>NULL()
CHARACTER (LEN=4), POINTER :: CMF_UPDRAFT=>NULL() 
CHARACTER (LEN=4), POINTER :: CMF_CLOUD=>NULL()
LOGICAL          , POINTER :: LMIXUV=>NULL() 
LOGICAL          , POINTER :: LMF_FLX=>NULL() 
!
REAL, POINTER          :: XALP_PERT=>NULL()   
REAL, POINTER          :: XABUO=>NULL()   
REAL, POINTER          :: XBENTR=>NULL()   
REAL, POINTER          :: XBDETR=>NULL()   
REAL, POINTER          :: XCMF=>NULL()    
REAL, POINTER          :: XENTR_MF=>NULL()   
REAL, POINTER          :: XCRAD_MF=>NULL()   
REAL, POINTER          :: XENTR_DRY=>NULL()    
REAL, POINTER          :: XDETR_DRY=>NULL()  
REAL, POINTER          :: XDETR_LUP=>NULL()  
REAL, POINTER          :: XKCF_MF=>NULL()    
REAL, POINTER          :: XKRC_MF=>NULL()    
REAL, POINTER          :: XTAUSIGMF=>NULL()
REAL, POINTER          :: XPRES_UV=>NULL()   
REAL, POINTER          :: XALPHA_MF=>NULL()   
REAL, POINTER          :: XSIGMA_MF=>NULL()  
REAL, POINTER          :: XFRAC_UP_MAX=>NULL()
REAL, POINTER          :: XA1=>NULL() 
REAL, POINTER          :: XB=>NULL()
REAL, POINTER          :: XC=>NULL()  
REAL, POINTER          :: XBETA1=>NULL()
REAL, POINTER          :: XR=>NULL() 
LOGICAL, POINTER       :: LTHETAS_MF=>NULL()
REAL, POINTER          :: XLAMBDA_MF=>NULL() 
LOGICAL, POINTER       :: LGZ=>NULL() 
REAL, POINTER          :: XGZ=>NULL() 
!
NAMELIST/NAM_PARAM_MFSHALLn/XIMPL_MF,CMF_UPDRAFT,CMF_CLOUD,LMIXUV,LMF_FLX,&
                            XALP_PERT,XABUO,XBENTR,XBDETR,XCMF,XENTR_MF,&
                            XCRAD_MF,XENTR_DRY,XDETR_DRY,XDETR_LUP,XKCF_MF,&
                            XKRC_MF,XTAUSIGMF,XPRES_UV,XALPHA_MF,XSIGMA_MF,&
                            XFRAC_UP_MAX,XA1,XB,XC,XBETA1,XR,LTHETAS_MF,LGZ,XGZ
!
!-------------------------------------------------------------------------------
!
CONTAINS

SUBROUTINE PARAM_MFSHALL_GOTO_MODEL(KFROM, KTO)
!! This subroutine associate all the pointers to the right component of
!! the right strucuture. A value can be accessed through the structure PARAM_MFSHALLN
!! or through the strucuture PARAM_MFSHALL_MODEL(KTO) or directly through these pointers.
INTEGER, INTENT(IN) :: KFROM, KTO
!
IF(.NOT. ASSOCIATED(PARAM_MFSHALLN, PARAM_MFSHALL_MODEL(KTO))) THEN
!
PARAM_MFSHALLN => PARAM_MFSHALL_MODEL(KTO)
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
XIMPL_MF=>PARAM_MFSHALL_MODEL(KTO)%XIMPL_MF
CMF_UPDRAFT=>PARAM_MFSHALL_MODEL(KTO)%CMF_UPDRAFT
CMF_CLOUD=>PARAM_MFSHALL_MODEL(KTO)%CMF_CLOUD
LMIXUV=>PARAM_MFSHALL_MODEL(KTO)%LMIXUV
LMF_FLX=>PARAM_MFSHALL_MODEL(KTO)%LMF_FLX
!
XALP_PERT=>PARAM_MFSHALL_MODEL(KTO)%XALP_PERT
XABUO=>PARAM_MFSHALL_MODEL(KTO)%XABUO
XBENTR=>PARAM_MFSHALL_MODEL(KTO)%XBENTR
XBDETR=>PARAM_MFSHALL_MODEL(KTO)%XBDETR
XCMF=>PARAM_MFSHALL_MODEL(KTO)%XCMF
XENTR_MF=>PARAM_MFSHALL_MODEL(KTO)%XENTR_MF
XCRAD_MF=>PARAM_MFSHALL_MODEL(KTO)%XCRAD_MF
XENTR_DRY=>PARAM_MFSHALL_MODEL(KTO)%XENTR_DRY
XDETR_DRY=>PARAM_MFSHALL_MODEL(KTO)%XDETR_DRY
XDETR_LUP=>PARAM_MFSHALL_MODEL(KTO)%XDETR_LUP
XKCF_MF=>PARAM_MFSHALL_MODEL(KTO)%XKCF_MF
XKRC_MF=>PARAM_MFSHALL_MODEL(KTO)%XKRC_MF
XTAUSIGMF=>PARAM_MFSHALL_MODEL(KTO)%XTAUSIGMF
XPRES_UV=>PARAM_MFSHALL_MODEL(KTO)%XPRES_UV
XALPHA_MF=>PARAM_MFSHALL_MODEL(KTO)%XALPHA_MF
XSIGMA_MF=>PARAM_MFSHALL_MODEL(KTO)%XSIGMA_MF
XFRAC_UP_MAX=>PARAM_MFSHALL_MODEL(KTO)%XFRAC_UP_MAX
XA1=>PARAM_MFSHALL_MODEL(KTO)%XA1
XB=>PARAM_MFSHALL_MODEL(KTO)%XB
XC=>PARAM_MFSHALL_MODEL(KTO)%XC
XBETA1=>PARAM_MFSHALL_MODEL(KTO)%XBETA1
XR=>PARAM_MFSHALL_MODEL(KTO)%XR
LTHETAS_MF=>PARAM_MFSHALL_MODEL(KTO)%LTHETAS_MF
XLAMBDA_MF=>PARAM_MFSHALL_MODEL(KTO)%XLAMBDA_MF
LGZ=>PARAM_MFSHALL_MODEL(KTO)%LGZ
XGZ=>PARAM_MFSHALL_MODEL(KTO)%XGZ
!
ENDIF
!
END SUBROUTINE PARAM_MFSHALL_GOTO_MODEL

SUBROUTINE PARAM_MFSHALLN_INIT(HPROGRAM, TFILENAM, LDNEEDNAM, KLUOUT, &
                         &LDDEFAULTVAL, LDREADNAM, LDCHECK, KPRINT)
!!*** *PARAM_MFSHALLN* - Code needed to initialize the MODD_PARAM_MFSHALL_n module
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
USE MODE_CHECK_NAM_VAL, ONLY: CHECK_NAM_VAL_CHAR
USE MODD_IO,  ONLY: TFILEDATA
!
IMPLICIT NONE
!
!* 0.1. Declaration of arguments
!       ------------------------
!
CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM     !< Name of the calling program
TYPE(TFILEDATA),   INTENT(IN) :: TFILENAM     !< Namelist file
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

  XIMPL_MF=1.
  CMF_UPDRAFT='EDKF'
  CMF_CLOUD='DIRE'
  LMIXUV=.TRUE.
  LMF_FLX=.FALSE.
  XALP_PERT=0.3
  XABUO=1.
  XBENTR=1.
  XBDETR=0.
  XCMF=0.065
  XENTR_MF=0.035
  XCRAD_MF=50.
  XENTR_DRY=0.55
  XDETR_DRY=10.
  XDETR_LUP=1.
  XKCF_MF=2.75
  XKRC_MF=1.
  XTAUSIGMF=600.
  XPRES_UV=0.5
  XALPHA_MF=2.
  XSIGMA_MF=20.
  XFRAC_UP_MAX=0.33
  XA1=2./3.
  XB=0.002
  XC=0.012
  XBETA1=0.9
  XR=2.
  LTHETAS_MF=.FALSE.
  XLAMBDA_MF=0.
  LGZ=.FALSE.
  XGZ=1.83 ! between 1.83 and 1.33
ENDIF
!
!*      2. NAMELIST
!       -----------
!
IF(LLREADNAM) THEN
  CALL POSNAM_PHY(TFILENAM, 'NAM_PARAM_MFSHALLN', LDNEEDNAM, LLFOUND)
  IF(LLFOUND) READ(UNIT=TFILENAM%NLU, NML=NAM_PARAM_MFSHALLn)
ENDIF
!
!*      3. CHECKS
!       ---------
!
IF(LLCHECK) THEN
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CMF_CLOUD', CMF_CLOUD, 'NONE', 'STAT', 'DIRE', 'BIGA')
  CALL CHECK_NAM_VAL_CHAR(KLUOUT, 'CMF_UPDRAFT', CMF_UPDRAFT, 'NONE', 'EDKF', 'RHCJ', 'RAHA')
ENDIF
!
!*      3. PRINTS
!       ---------
!
IF(IPRINT>=1) THEN
  WRITE(UNIT=KLUOUT,NML=NAM_PARAM_MFSHALLn)
ENDIF
!
END SUBROUTINE PARAM_MFSHALLN_INIT
!
END MODULE MODD_PARAM_MFSHALL_n
