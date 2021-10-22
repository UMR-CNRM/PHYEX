!     ######spl
      MODULE MODD_PARAM_ICE
!     #####################
!
!!****  *MODD_PARAM_ICE* - declaration of the control parameters for the
!!                           mixed phase cloud parameterization
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to define the set of space
!!    and time control parameters for the microphysics.
!!
!!
!!**  IMPLICIT ARGUMENTS
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
!!      Original      14/12/95
!!      Jan 2015 S. Riette: new ICE3/ICE4 parameters
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
LOGICAL, SAVE :: LWARM       ! When .TRUE. activates the formation of rain by
                             ! the warm microphysical processes
LOGICAL, SAVE :: LSEDIC                 ! TRUE to enable the droplet sedimentation
!
CHARACTER(LEN=4), SAVE :: CPRISTINE_ICE ! Pristine ice type PLAT, COLU or BURO
CHARACTER(LEN=4), SAVE :: CSEDIM        ! Sedimentation calculation mode      
!
LOGICAL, SAVE :: LFEEDBACKT ! When .TRUE. feed back on temperature is taken into account
LOGICAL, SAVE :: LEVLIMIT   ! When .TRUE. water vapour pressure is limited by saturation
LOGICAL, SAVE :: LNULLWETG  ! When .TRUE. graupel wet growth is activated with null rate (to allow water shedding)
LOGICAL, SAVE :: LWETGPOST  ! When .TRUE. graupel wet growth is activated with positive temperature (to allow water shedding)
LOGICAL, SAVE :: LNULLWETH  ! Same as LNULLWETG but for hail
LOGICAL, SAVE :: LWETHPOST  ! Same as LWETGPOST but for hail
CHARACTER(LEN=4), SAVE :: CSNOWRIMING ! OLD or M90 for Murakami 1990 formulation
REAL, SAVE :: XFRACM90 ! Fraction used for the Murakami 1990 formulation
INTEGER, SAVE :: NMAXITER ! Maximum number of iterations for mixing ratio or time splitting
REAL, SAVE :: XMRSTEP ! maximum mixing ratio step for mixing ratio splitting
LOGICAL, SAVE :: LCONVHG ! TRUE to allow the conversion from hail to graupel
LOGICAL, SAVE :: LCRFLIMIT !True to limit rain contact freezing to possible heat exchange
!
REAL, SAVE :: XTSTEP_TS ! Approximative time step for time-splitting (0 for no time-splitting)
!
CHARACTER*80, SAVE :: CSUBG_RC_RR_ACCR ! subgrid rc-rr accretion
CHARACTER*80, SAVE :: CSUBG_RR_EVAP ! subgrid rr evaporation
CHARACTER*80, SAVE :: CSUBG_PR_PDF ! pdf for subgrid precipitation
!
LOGICAL, SAVE :: LADJ_BEFORE ! must we perform an adjustment before rain_ice call
LOGICAL, SAVE :: LADJ_AFTER ! must we perform an adjustment after rain_ice call
CHARACTER*1, SAVE :: CFRAC_ICE_ADJUST ! ice fraction for adjustments
CHARACTER*1, SAVE :: CFRAC_ICE_SHALLOW_MF ! ice fraction for shallow_mf
LOGICAL, SAVE :: LSEDIM_AFTER ! sedimentation done before (.FALSE.) or after (.TRUE.) microphysics
!
LOGICAL, SAVE :: LDEPOSC ! switch on water deposition on vegetation
REAL, SAVE    :: XVDEPOSC ! water deposition speed

REAL, SAVE :: XSPLIT_MAXCFL ! Maximum CFL number allowed for SPLIT scheme
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_PARAM_ICE
