!     ######spl
      MODULE MODD_CONVPAR
!     ###################
!
!!****  *MODD_CONVPAR* - Declaration of convection constants
!!
!!    PURPOSE
!!    -------
!      The purpose of this declarative module is to declare  the
!      constants in the deep convection parameterization.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_CONVPAR)
!!
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  15/11/96
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
REAL, SAVE :: XA25        ! 25 km x 25 km reference grid area
!
REAL, SAVE :: XCRAD       ! cloud radius
REAL, SAVE :: XCDEPTH     ! minimum necessary cloud depth
REAL, SAVE :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)
!
REAL, SAVE :: XZLCL       ! maximum allowed allowed height
                          ! difference between departure level and surface
REAL, SAVE :: XZPBL       ! minimum mixed layer depth to sustain convection
REAL, SAVE :: XWTRIG      ! constant in vertical velocity trigger
!
!
REAL, SAVE :: XNHGAM      ! accounts for non-hydrost. pressure
                          ! in buoyancy term of w equation
                          ! = 2 / (1+gamma)
REAL, SAVE :: XTFRZ1      ! begin of freezing interval
REAL, SAVE :: XTFRZ2      ! end of freezing interval
!
REAL, SAVE :: XRHDBC      ! relative humidity below cloud in downdraft
!
REAL, SAVE :: XRCONV      ! constant in precipitation conversion
REAL, SAVE :: XSTABT      ! factor to assure stability in  fractional time
                          ! integration, routine CONVECT_CLOSURE
REAL, SAVE :: XSTABC      ! factor to assure stability in CAPE adjustment,
                          !  routine CONVECT_CLOSURE
REAL, SAVE :: XUSRDPTH    ! pressure thickness used to compute updraft
                          ! moisture supply rate for downdraft
REAL, SAVE :: XMELDPTH    ! layer (Pa) through which precipitation melt is
                          ! allowed below  melting level
REAL, SAVE :: XUVDP       ! constant for pressure perturb in momentum transport
!
END MODULE MODD_CONVPAR
