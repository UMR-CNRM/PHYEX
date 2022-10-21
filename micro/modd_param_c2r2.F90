!     ######spl
      MODULE MODD_PARAM_C2R2
!     ######################
!
!!****  *MODD_PARAM_C2R2* - declaration of the control parameters
!!                               for use in the warm scheme.
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the microphysical
!     constants. This includes the descriptive parameters for the raindrop
!     and the parameters relevant of the dimensional distributions.
!
!       The four constants used to set the actiavtion spectrum are also
!     defined
!
!          N_activated = (C*S**k)*F_hypgeo(mu,k/2,k/2+1;-beta*S**2)
!         DN_activated/DS = k*(C*S**(k-1))*(1+beta*S**2)**-mu
!
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_PARAM_C2R2)
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty  *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/11/2000
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
REAL,SAVE :: XALPHAR,XNUR,           & ! Raindrop       distribution parameters
             XALPHAC,XNUC              ! Cloud droplet  distribution parameters
!
LOGICAL, SAVE :: LRAIN                 ! TRUE to enable the formation of rain
LOGICAL, SAVE :: LSEDC                 ! TRUE to enable the droplet sedimentation
LOGICAL, SAVE :: LACTIT                ! TRUE to enable the usage of
                                       ! dT/dt in CCN activation (twomey and CPB98)
!
REAL,SAVE ::  XCHEN,XKHEN,           & ! Parameters used to define the CCN
              XMUHEN,XBETAHEN          ! activation spectra (CPB or TWO)
!
CHARACTER(LEN=3),SAVE :: HPARAM_CCN    ! Parameterization used for the CCN activation
CHARACTER(LEN=3),SAVE :: HINI_CCN      ! Initialization type of the CCN activation
CHARACTER(LEN=1),SAVE :: HTYPE_CCN     ! 'M' or 'C' standard type of CCN
REAL,SAVE ::  XCONC_CCN,             & ! Concentration of the CCN
              XR_MEAN_CCN,           & ! Geometric mean radius of the CCN
              XLOGSIG_CCN,           & ! Log of geometric dispersion of the CCN
              XFSOLUB_CCN,           & ! Fractionnal solubility of the CCN
              XACTEMP_CCN,           & ! Expected temperature of CCN activation
              XAERDIFF, XAERHEIGHT     ! For the vertical gradient of
                                       ! aerosol distribution
END MODULE MODD_PARAM_C2R2
