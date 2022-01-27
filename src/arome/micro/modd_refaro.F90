!     ######spl
      MODULE MODD_REF
!     ###############
!
!!****  *MODD_REF* - declaration of reference state  profile
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the vertical
!     profile of  the reference state, used for the anelastic
!     approximation.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_REF)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   07/06/94
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!  temperature parameter used in INI_RAIN_ICE,
!  set to constant value 300k for AROME
REAL,SAVE, DIMENSION(2) :: XTHVREFZ=300.     ! Thetav(z) for reference
                                             ! state without orography
!
END MODULE MODD_REF
