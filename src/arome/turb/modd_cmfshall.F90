!     ######spl
      MODULE MODD_CMFSHALL
!     #############################
!
!!****  *MODD_CMFSHALL* - Declaration of Mass flux scheme constants 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to declare some
!!      constants for Mass Flux Scheme in the shallow convection 
!!      parameterization.  
!
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
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------

IMPLICIT NONE

REAL,SAVE          :: XALP_PERT   ! coefficient for the perturbation of
                                ! theta_l and r_t at the first level of 
                                ! the updraft
REAL,SAVE          ::    XABUO    ! coefficient of the buoyancy term in the w_up equation
REAL,SAVE          ::    XBENTR   ! coefficient of the entrainment term in the w_up equation
REAL,SAVE          ::    XBDETR   ! coefficient of the detrainment term in the w_up equation
REAL,SAVE          ::    XCMF     ! coefficient for the mass flux at the first level 
                                ! of the updraft (closure)
REAL,SAVE          :: XENTR_MF    ! entrainment constant (m/Pa) = 0.2 (m) 
REAL,SAVE          :: XCRAD_MF    ! cloud radius in cloudy part
REAL,SAVE          :: XENTR_DRY   ! coefficient for entrainment in dry part 
REAL,SAVE          :: XDETR_DRY   ! coefficient for detrainment in dry part
REAL,SAVE          :: XDETR_LUP   ! coefficient for detrainment in dry part
REAL,SAVE          :: XKCF_MF     ! coefficient for cloud fraction
REAL,SAVE          :: XKRC_MF     ! coefficient for convective rc
REAL,SAVE          :: XTAUSIGMF
REAL,SAVE          :: XPRES_UV    ! coefficient for pressure term in wind
                                  ! mixing

REAL,SAVE          :: XALPHA_MF   ! coefficient for cloudy fraction
REAL,SAVE          :: XSIGMA_MF   ! coefficient for sigma computation

REAL,SAVE          :: XFRAC_UP_MAX! maximum Updraft fraction


!  Parameter for Rio et al (2010) formulation for entrainment and detrainment

REAL,SAVE          :: XA1      ! a1
REAL,SAVE          :: XB       ! b
REAL,SAVE          :: XC       ! c
REAL,SAVE          :: XBETA1   ! beta1

!  Parameters for closure assumption of Hourdin et al 2002

REAL,SAVE          :: XR      ! Aspect ratio of updraft

!  Thermodynamic parameter

REAL,SAVE          :: XLAMBDA      ! Lambda to compute ThetaS1 from ThetaL

END MODULE MODD_CMFSHALL
