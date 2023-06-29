!     ######spl
      SUBROUTINE INI_CONVPAR(CONVPAR)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ######################
!
!!****  *INI_CONVPAR * - routine to initialize the constants modules 
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize  the constants
!     stored in  modules MODD_CONVPAR, MODD_CST, MODD_CONVPAREXT.
!      
!
!!**  METHOD
!!    ------
!!      The deep convection constants are set to their numerical values 
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONVPAR   : contains deep convection constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module MODD_CONVPAR, routine INI_CONVPAR)
!!      
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Last modified  15/04/98 adapted for ARPEGE
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAR, ONLY : CONVPAR_T
IMPLICIT NONE
TYPE(CONVPAR_T), INTENT(INOUT) :: CONVPAR

!
!  
!-------------------------------------------------------------------------------
!
!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------
!
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('INI_CONVPAR',0,ZHOOK_HANDLE)
CONVPAR%XA25     = 625.E6    ! 25 km x 25 km reference grid area
!
CONVPAR%XCRAD    = 1500.     ! cloud radius 
CONVPAR%XCDEPTH  = 2.5E3      ! minimum necessary cloud depth
CONVPAR%XENTR    = 0.03      ! entrainment constant (m/Pa) = 0.2 (m)  
!
CONVPAR%XZLCL    = 3.5E3     ! maximum allowed allowed height 
                     ! difference between the surface and the LCL
CONVPAR%XZPBL    = 60.E2     ! minimum mixed layer depth to sustain convection
CONVPAR%XWTRIG   = 6.00      ! constant in vertical velocity trigger
!
!
CONVPAR%XNHGAM   = 1.3333    ! accounts for non-hydrost. pressure 
                     ! in buoyancy term of w equation
                     ! = 2 / (1+gamma)
CONVPAR%XTFRZ1   = 268.16    ! begin of freezing interval
CONVPAR%XTFRZ2   = 248.16    ! end of freezing interval
!
CONVPAR%XRHDBC   = 0.9       ! relative humidity below cloud in downdraft

CONVPAR%XRCONV   = 0.015     ! constant in precipitation conversion 
CONVPAR%XSTABT   = 0.75      ! factor to assure stability in  fractional time
                     ! integration, routine CONVECT_CLOSURE
CONVPAR%XSTABC   = 0.95      ! factor to assure stability in CAPE adjustment,
                     !  routine CONVECT_CLOSURE
CONVPAR%XUSRDPTH = 165.E2    ! pressure thickness used to compute updraft
                     ! moisture supply rate for downdraft
CONVPAR%XMELDPTH = 100.E2    ! layer (Pa) through which precipitation melt is
                     ! allowed below downdraft
CONVPAR%XUVDP    = 0.7       ! constant for pressure perturb in momentum transport
!
!
IF (LHOOK) CALL DR_HOOK('INI_CONVPAR',1,ZHOOK_HANDLE)
END SUBROUTINE INI_CONVPAR 
