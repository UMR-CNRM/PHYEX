!     ######spl
      MODULE MODD_CONVPAR_SHAL
!     ########################
!
!!****  *MODD_CONVPAR_SHAL* - Declaration of convection constants
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to declare  the
!!      constants in the deep convection parameterization.
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_CONVPAR_SHAL)
!!
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  04/10/98
!!      E. Bazile   05/05/09
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
TYPE CONVPAR_SHAL
REAL :: XA25        ! 25 km x 25 km reference grid area
!
REAL :: XCRAD       ! cloud radius
REAL :: XCTIME_SHAL ! convective adjustment time
REAL :: XCDEPTH     ! minimum necessary cloud depth
REAL :: XCDEPTH_D   ! maximum allowed cloud thickness
REAL :: XDTPERT     ! add small Temp perturb. at LCL
REAL :: XATPERT     ! Parameter for temp Perturb
REAL :: XBTPERT     ! Parameter for temp Perturb
                    ! (XATPERT* TKE/Cp + XBTPERT) * XDTPERT
REAL :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)
!
REAL :: XZLCL       ! maximum allowed allowed height
                    ! difference between departure level and surface
REAL :: XZPBL       ! minimum mixed layer depth to sustain convection
REAL :: XWTRIG      ! constant in vertical velocity trigger
!
!
REAL :: XNHGAM      ! accounts for non-hydrost. pressure
                    ! in buoyancy term of w equation
                    ! = 2 / (1+gamma)
REAL :: XTFRZ1      ! begin of freezing interval
REAL :: XTFRZ2      ! end of freezing interval
!
!
REAL :: XSTABT      ! factor to assure stability in  fractional time
                    ! integration, routine CONVECT_CLOSURE
REAL :: XSTABC      ! factor to assure stability in CAPE adjustment,
                    !  routine CONVECT_CLOSURE
REAL :: XAW,XBW     ! Parameters for WLCL = XAW * W + XBW
LOGICAL :: LLSMOOTH ! Default=TRUE but not necessary
END TYPE CONVPAR_SHAL
!Keep global variables for parts of the code not ported to the type ye
REAL :: XA25        ! 25 km x 25 km reference grid area
!
REAL :: XCRAD       ! cloud radius
REAL :: XCTIME_SHAL ! convective adjustment time
REAL :: XCDEPTH     ! minimum necessary cloud depth
REAL :: XCDEPTH_D   ! maximum allowed cloud thickness
REAL :: XDTPERT     ! add small Temp perturb. at LCL
REAL :: XATPERT     ! Parameter for temp Perturb
REAL :: XBTPERT     ! Parameter for temp Perturb
                    ! (XATPERT* TKE/Cp + XBTPERT) * XDTPERT
REAL :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)
!
REAL :: XZLCL       ! maximum allowed allowed height
                    ! difference between departure level and surface
REAL :: XZPBL       ! minimum mixed layer depth to sustain convection
REAL :: XWTRIG      ! constant in vertical velocity trigger
!
!
REAL :: XNHGAM      ! accounts for non-hydrost. pressure
                    ! in buoyancy term of w equation
                    ! = 2 / (1+gamma)
REAL :: XTFRZ1      ! begin of freezing interval
REAL :: XTFRZ2      ! end of freezing interval
!
!
REAL :: XSTABT      ! factor to assure stability in  fractional time
                    ! integration, routine CONVECT_CLOSURE
REAL :: XSTABC      ! factor to assure stability in CAPE adjustment,
                    !  routine CONVECT_CLOSURE
REAL :: XAW,XBW     ! Parameters for WLCL = XAW * W + XBW
LOGICAL :: LLSMOOTH ! Default=TRUE but not necessary
!
END MODULE MODD_CONVPAR_SHAL
