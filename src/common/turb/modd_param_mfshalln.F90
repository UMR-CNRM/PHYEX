!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     #############################
      MODULE MODD_PARAM_MFSHALL_n
!     #############################
!
!!****  *MODD_PARAM_MFSHALL_n* - Declaration of Mass flux scheme free parameters
!!
!!    PURPOSE
!!    -------
!!    The purpose of this declarative module is to declare the
!!    variables that may be set by namelist for the mass flux scheme
!!
!!**  IMPLICIT ARGUMENTS
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

REAL               :: XIMPL_MF     ! degre of implicitness
      
CHARACTER (LEN=4)  :: CMF_UPDRAFT  ! Type of Mass Flux Scheme
                                     ! 'NONE' if no parameterization 
CHARACTER (LEN=4)  :: CMF_CLOUD

                                     
LOGICAL       :: LMIXUV      ! True if mixing of momentum
LOGICAL       :: LMF_FLX     ! logical switch for the storage of
                             ! the mass flux fluxes
REAL          :: XALP_PERT   ! coefficient for the perturbation of
                             ! theta_l and r_t at the first level of 
                             ! the updraft
REAL          ::    XABUO    ! coefficient of the buoyancy term in the w_up equation
REAL          ::    XBENTR   ! coefficient of the entrainment term in the w_up equation
REAL          ::    XBDETR   ! coefficient of the detrainment term in the w_up equation
REAL          ::    XCMF     ! coefficient for the mass flux at the first level 
                             ! of the updraft (closure)
REAL          :: XENTR_MF    ! entrainment constant (m/Pa) = 0.2 (m) 
REAL          :: XCRAD_MF    ! cloud radius in cloudy part
REAL          :: XENTR_DRY   ! coefficient for entrainment in dry part 
REAL          :: XDETR_DRY   ! coefficient for detrainment in dry part
REAL          :: XDETR_LUP   ! coefficient for detrainment in dry part
REAL          :: XKCF_MF     ! coefficient for cloud fraction
REAL          :: XKRC_MF     ! coefficient for convective rc
REAL          :: XTAUSIGMF
REAL          :: XPRES_UV    ! coefficient for pressure term in wind
                             ! mixing
REAL          :: XALPHA_MF   ! coefficient for cloudy fraction
REAL          :: XSIGMA_MF   ! coefficient for sigma computation
REAL          :: XFRAC_UP_MAX! maximum Updraft fraction
!
!  Parameter for Rio et al (2010) formulation for entrainment and detrainment (RHCJ10)
REAL          :: XA1
REAL          :: XB
REAL          :: XC
REAL          :: XBETA1
!
!  Parameters for closure assumption of Hourdin et al 2002

REAL          :: XR      ! Aspect ratio of updraft
!
!  Grey Zone 
LOGICAL       :: LGZ      ! Grey Zone Surface Closure
REAL          :: XGZ      ! Tuning of the surface initialisation
!
!  Thermodynamic parameter
REAL          :: XLAMBDA_MF      ! Lambda to compute ThetaS1 from ThetaL

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
REAL, POINTER          :: XLAMBDA_MF=>NULL() 
LOGICAL, POINTER       :: LGZ=>NULL() 
REAL, POINTER          :: XGZ=>NULL() 
CONTAINS

SUBROUTINE PARAM_MFSHALL_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
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
XLAMBDA_MF=>PARAM_MFSHALL_MODEL(KTO)%XLAMBDA_MF
LGZ=>PARAM_MFSHALL_MODEL(KTO)%LGZ
XGZ=>PARAM_MFSHALL_MODEL(KTO)%XGZ
!
END SUBROUTINE PARAM_MFSHALL_GOTO_MODEL

END MODULE MODD_PARAM_MFSHALL_n
