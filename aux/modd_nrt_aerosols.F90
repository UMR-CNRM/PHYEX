!     ######spl
      MODULE MODD_NRT_AEROSOLS
!     #####################
!
!!****  *MODD_NRT_AEROSOLS* - declaration of the control parameters for 
!!                           the use of near real time aerosols
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to define the set of 
!!    parameters for the use of nrt aerosols in the microphysics.
!!
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
!!       D. Martin-Perez
!!
!!    MODIFICATIONS
!!    -------------
!!       Original 05/08/2022
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
LOGICAL, SAVE :: LAEIFN       ! switch to activate ice nuclei
LOGICAL, SAVE :: LAECCN2CLDR  ! switch to convert CCN into Cloud
                              ! droplets
LOGICAL, SAVE :: LAERDRDEP    ! Switch for dry deposition activation
LOGICAL, SAVE :: LAERSSEM     ! Switch for sea salt emission.
REAL,    SAVE :: XSSMINLO     ! Minimum supersaturation at lowest level
REAL,    SAVE :: XSSMINUP     ! Minimum supersaturation at high levels
REAL,    SAVE :: XSSMAX       ! Maximum supersaturation
REAL,    SAVE :: XSSHEIGHT    ! Height limit for SSMIN separation
REAL,    SAVE :: XSSFACVV     ! Factor for SS dependence with vert. velo.
REAL,    SAVE :: XSSFACSS     ! Factor for SS dependence with coarse sea salt.
REAL,    SAVE :: XPANMIN      ! Minimum particle concentration number
                              ! (m-3)
REAL,    SAVE :: XCCNMIN      ! Minimum cloud condensation number
                              ! concentration (m-3)
REAL,    SAVE :: XCLDROPMIN   ! Minimum cloud droplet number
                              ! concentration (m-3)
REAL,    SAVE :: XIFNMINSIZE  ! Minimum size of aerosol to be considered
                              ! ice nuclei (micrometers)

!
!-------------------------------------------------------------------------------
!
END MODULE MODD_NRT_AEROSOLS

