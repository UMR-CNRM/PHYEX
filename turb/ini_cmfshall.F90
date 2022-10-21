!     ######spl
        SUBROUTINE INI_CMFSHALL(PALP_PERT,PABUO,PBENTR,PBDETR,PCMF,PENTR_MF,PCRAD_MF,PENTR_DRY,&
 &          PDETR_DRY,PDETR_LUP,PKCF_MF,PKRC_MF,PTAUSIGMF,PPRES_UV,PFRAC_UP_MAX,&
 &          PALPHA_MF,PSIGMA_MF,PA1,PB,PC,PBETA1,PR,PLAMBDA)
        USE PARKIND1, ONLY : JPRB
        USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!       #######################
!
!!****     *INI_CMFSHALL*  - routine to initialize the mass flux scheme 
!!                        constants.
!!
!!      PURPOSE
!!      -------
!         The purpose of this routine is to initialize the mass flux  
!       scheme constants that are stored in module MODD_PARAM_MFSHALL_n
!
!!      METHOD
!!      ------
!!        The constants are set to their numerical values
!!
!!      EXTERNAL
!!      --------
!!        NONE
!!
!!      IMPLICIT ARGUMENTS
!!      ------------------
!!        Module MODD_CTURB
!!
!!      REFERENCE
!!      ---------
!!
!!      AUTHOR
!!      ------
!!        S. Malardel, J. Pergaud (Meteo France)
!!
!!      MODIFICATIONS
!!      -------------
!!        S. Riette april 2011 : XALPHA and XSIGMA added
!!        S. Riette Jan 2022: Merge with Méso-NH: MODD_MCFSHALL -> MODD_PARAM_MFSHALL_n
!! --------------------------------------------------------------------------
!
!*        0. DECLARATIONS
!            ------------
!
USE MODD_PARAM_MFSHALL_n
!
IMPLICIT NONE

REAL,   INTENT(IN)   :: PALP_PERT
REAL,   INTENT(IN)   :: PABUO
REAL,   INTENT(IN)   :: PBENTR
REAL,   INTENT(IN)   :: PBDETR
REAL,   INTENT(IN)   :: PCMF
REAL,   INTENT(IN)   :: PENTR_MF
REAL,   INTENT(IN)   :: PCRAD_MF
REAL,   INTENT(IN)   :: PENTR_DRY
REAL,   INTENT(IN)   :: PDETR_DRY
REAL,   INTENT(IN)   :: PDETR_LUP
REAL,   INTENT(IN)   :: PKCF_MF
REAL,   INTENT(IN)   :: PKRC_MF
REAL,   INTENT(IN)   :: PTAUSIGMF
REAL,   INTENT(IN)   :: PPRES_UV
REAL,   INTENT(IN)   :: PFRAC_UP_MAX
REAL,   INTENT(IN)   :: PALPHA_MF
REAL,   INTENT(IN)   :: PSIGMA_MF
REAL,   INTENT(IN)   :: PA1
REAL,   INTENT(IN)   :: PB
REAL,   INTENT(IN)   :: PC
REAL,   INTENT(IN)   :: PBETA1   
REAL,   INTENT(IN)   :: PR   
REAL,   INTENT(IN)   :: PLAMBDA   

!
!  ---------------------------------------------------------------------------
!
!         1. SETTING THE NUMERICAL VALUES
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('INI_CMFSHALL',0,ZHOOK_HANDLE)

CALL PARAM_MFSHALL_GOTO_MODEL(1, 1)

XALP_PERT   = PALP_PERT  ! coefficient for the perturbation of
                         ! theta_l and r_t at the first level of 
                         ! the updraft
XABUO       = PABUO      ! coefficient of the buoyancy term in the w_up equation  
XBENTR      = PBENTR     ! coefficient of the entrainment term in the w_up equation
XBDETR      = PBDETR     ! coefficient of the detrainment term in the w_up equation
XCMF        = PCMF       ! coefficient for the mass flux at the first level 
                         ! of the updraft (closure) 
XENTR_MF    = PENTR_MF   ! entrainment constant (m/Pa) = 0.2 (m)
XCRAD_MF    = PCRAD_MF   ! cloud radius in cloudy part
XENTR_DRY   = PENTR_DRY  ! coefficient for entrainment in dry part 
XDETR_DRY   = PDETR_DRY  ! coefficient for detrainment in dry part
XDETR_LUP   = PDETR_LUP  !  coefficient for detrainment in dry part
XKCF_MF     = PKCF_MF    ! coefficient for cloud fraction
XKRC_MF     = PKRC_MF    ! coefficient for convective rc
XTAUSIGMF   = PTAUSIGMF  
XPRES_UV    = PPRES_UV   ! coefficient for pressure term in wind mixing
!
XFRAC_UP_MAX= PFRAC_UP_MAX ! maximum Updraft fraction
!
XALPHA_MF = PALPHA_MF    ! coefficient for updraft fraction in STA2 cloud scheme
XSIGMA_MF = PSIGMA_MF    ! coefficient for sigma in STA2 cloud scheme

!  tuning variables for compute_updraft_rhcj10
XA1    =  PA1       !  Value Rio et al 2010
XB     =  PB        !  Value Rio et al 2010
XC     =  PC        !  Value Rio et al 2010
XBETA1 =  PBETA1    !  Value Rio et al 2010

!  Parameters for closure assumption of Hourdin et al 2002

XR = PR            ! Aspect ratio of updraft

!  Thermodynamic parameter

XLAMBDA_MF = PLAMBDA         ! Lambda to compute ThetaS1 from ThetaL

IF (LHOOK) CALL DR_HOOK('INI_CMFSHALL',1,ZHOOK_HANDLE)
END SUBROUTINE INI_CMFSHALL
