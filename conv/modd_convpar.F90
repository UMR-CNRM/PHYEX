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
TYPE CONVPAR_T
REAL :: XA25        ! 25 km x 25 km reference grid area
!
REAL :: XCRAD       ! cloud radius
REAL :: XCDEPTH     ! minimum necessary cloud depth
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
REAL :: XRHDBC      ! relative humidity below cloud in downdraft
!
REAL :: XRCONV      ! constant in precipitation conversion
REAL :: XSTABT      ! factor to assure stability in  fractional time
                    ! integration, routine CONVECT_CLOSURE
REAL :: XSTABC      ! factor to assure stability in CAPE adjustment,
                    !  routine CONVECT_CLOSURE
REAL :: XUSRDPTH    ! pressure thickness used to compute updraft
                    ! moisture supply rate for downdraft
REAL :: XMELDPTH    ! layer (Pa) through which precipitation melt is
                    ! allowed below  melting level
REAL :: XUVDP       ! constant for pressure perturb in momentum transport
END TYPE CONVPAR_T

!Keep global variables for parts of the code not ported to the type yet
REAL :: XA25        ! 25 km x 25 km reference grid area
!
REAL :: XCRAD       ! cloud radius
REAL :: XCDEPTH     ! minimum necessary cloud depth
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
REAL :: XRHDBC      ! relative humidity below cloud in downdraft
!
REAL :: XRCONV      ! constant in precipitation conversion
REAL :: XSTABT      ! factor to assure stability in  fractional time
                    ! integration, routine CONVECT_CLOSURE
REAL :: XSTABC      ! factor to assure stability in CAPE adjustment,
                    !  routine CONVECT_CLOSURE
REAL :: XUSRDPTH    ! pressure thickness used to compute updraft
                    ! moisture supply rate for downdraft
REAL :: XMELDPTH    ! layer (Pa) through which precipitation melt is
                    ! allowed below  melting level
REAL :: XUVDP       ! constant for pressure perturb in momentum transport

INTERFACE INI_CONVPAR
  MODULE PROCEDURE INI_CONVPAR0
  MODULE PROCEDURE INI_CONVPAR1
END INTERFACE



CONTAINS


SUBROUTINE INI_CONVPAR1 (CONVPAR)

USE YOMHOOK , ONLY : LHOOK, JPHOOK, DR_HOOK

TYPE (CONVPAR_T), INTENT (OUT) :: CONVPAR

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('INI_CONVPAR1',0,ZHOOK_HANDLE)

CONVPAR%XA25     = XA25        
CONVPAR%XCRAD    = XCRAD       
CONVPAR%XCDEPTH  = XCDEPTH     
CONVPAR%XENTR    = XENTR       
CONVPAR%XZLCL    = XZLCL       
CONVPAR%XZPBL    = XZPBL       
CONVPAR%XWTRIG   = XWTRIG      
CONVPAR%XNHGAM   = XNHGAM      
CONVPAR%XTFRZ1   = XTFRZ1      
CONVPAR%XTFRZ2   = XTFRZ2      
CONVPAR%XRHDBC   = XRHDBC      
CONVPAR%XRCONV   = XRCONV      
CONVPAR%XSTABT   = XSTABT      
CONVPAR%XSTABC   = XSTABC      
CONVPAR%XUSRDPTH = XUSRDPTH    
CONVPAR%XMELDPTH = XMELDPTH    
CONVPAR%XUVDP    = XUVDP       

IF (LHOOK) CALL DR_HOOK('INI_CONVPAR1',1,ZHOOK_HANDLE)

END SUBROUTINE INI_CONVPAR1

SUBROUTINE INI_CONVPAR0
USE YOMHOOK , ONLY : LHOOK, JPHOOK, DR_HOOK
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
!
IMPLICIT NONE
!  
!-------------------------------------------------------------------------------
!
!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------
!
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('INI_CONVPAR0',0,ZHOOK_HANDLE)
XA25     = 625.E6    ! 25 km x 25 km reference grid area
!
XCRAD    = 1500.     ! cloud radius 
XCDEPTH  = 2.5E3      ! minimum necessary cloud depth
XENTR    = 0.03      ! entrainment constant (m/Pa) = 0.2 (m)  
!
XZLCL    = 3.5E3     ! maximum allowed allowed height 
                     ! difference between the surface and the LCL
XZPBL    = 60.E2     ! minimum mixed layer depth to sustain convection
XWTRIG   = 6.00      ! constant in vertical velocity trigger
!
!
XNHGAM   = 1.3333    ! accounts for non-hydrost. pressure 
                     ! in buoyancy term of w equation
                     ! = 2 / (1+gamma)
XTFRZ1   = 268.16    ! begin of freezing interval
XTFRZ2   = 248.16    ! end of freezing interval
!
XRHDBC   = 0.9       ! relative humidity below cloud in downdraft

XRCONV   = 0.015     ! constant in precipitation conversion 
XSTABT   = 0.75      ! factor to assure stability in  fractional time
                     ! integration, routine CONVECT_CLOSURE
XSTABC   = 0.95      ! factor to assure stability in CAPE adjustment,
                     !  routine CONVECT_CLOSURE
XUSRDPTH = 165.E2    ! pressure thickness used to compute updraft
                     ! moisture supply rate for downdraft
XMELDPTH = 100.E2    ! layer (Pa) through which precipitation melt is
                     ! allowed below downdraft
XUVDP    = 0.7       ! constant for pressure perturb in momentum transport
!
!
IF (LHOOK) CALL DR_HOOK('INI_CONVPAR0',1,ZHOOK_HANDLE)
END SUBROUTINE INI_CONVPAR0
!
END MODULE MODD_CONVPAR
