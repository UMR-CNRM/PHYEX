!     ######spl
      MODULE MODD_AEROSOL_PROP
!     ##########################
!
!!****  *MODD_AEROSOL_PROP* - declaration of aerosol properties.
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare ...
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/12/17
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE

! CAMS Aerosols   : 14 species
! MOCAGE Aerosols : ???

INTEGER, PARAMETER :: NAEROMAX=16

REAL,SAVE :: XCONC_MIN   ! Minimun droplet concentration

INTEGER, SAVE :: NCCN
INTEGER, SAVE :: NIFN
INTEGER, SAVE :: NCOAR
INTEGER, SAVE :: NSEASALT
INTEGER, SAVE :: NCOARSEASALT
INTEGER, SAVE :: NDUST
INTEGER, SAVE :: NDRYDEP

! Cloud condensation nuclei
INTEGER, DIMENSION(NAEROMAX), SAVE :: XACCNI
! Ice freezing Nuclei
INTEGER, DIMENSION(NAEROMAX), SAVE :: XIFNI
! Coarse modes for dry sedimentation
INTEGER, DIMENSION(NAEROMAX), SAVE :: XCOARSE
! Sea salt species
INTEGER, DIMENSION(NAEROMAX), SAVE :: XSEASALT
INTEGER, DIMENSION(NAEROMAX), SAVE :: XCOARSEASALT
! Desert dust
INTEGER, DIMENSION(NAEROMAX), SAVE :: XDUST
! Dry deposition active
INTEGER, DIMENSION(NAEROMAX), SAVE :: XDRYDEP

! Aerosol distribution: log normal size distribution 
!     number mode radius(rm)
!     geometric standard deviation(s)
REAL, DIMENSION(NAEROMAX), SAVE :: XMRAE, XSDAE
! Variables related with the bin limits
REAL, DIMENSION(NAEROMAX), SAVE :: XBINDOWN,XBINUP
REAL, DIMENSION(NAEROMAX), SAVE :: XERFDOWN,XERFUP,XR3CORR

! Physical properties
REAL, DIMENSION(NAEROMAX), SAVE :: XKHYGROS    ! hygroscopic factor
REAL, DIMENSION(NAEROMAX), SAVE :: XRHOAE      ! aerosol density (Kg m-3)
! Optical properties
!     Mass extinction
REAL, DIMENSION(NAEROMAX,12), SAVE :: XMEXT    ! mass extinction

REAL, DIMENSION(NAEROMAX,2), SAVE :: XVDRY     ! Dry Deposition Velocity
                                               ! over sea and land
REAL, DIMENSION(NAEROMAX), SAVE :: XFRICS      ! Fracion of aerosols 
                                               ! in aquaeous phase
                                               ! for in-cloud scavenging

REAL, DIMENSION(NAEROMAX), SAVE :: XINTSS      ! Sea salt emision

END MODULE MODD_AEROSOL_PROP
