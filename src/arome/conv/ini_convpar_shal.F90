!     ######spl
      SUBROUTINE INI_CONVPAR_SHAL
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ###########################
!
!!****  *INI_CONVPAR * - routine to initialize the constants modules 
!!
!!    PURPOSE
!!    -------
!!       The purpose of this routine is to initialize  the constants
!!     stored in  modules MODD_CONVPAR_SHAL
!!      
!!
!!**  METHOD
!!    ------
!!      The shallow convection constants are set to their numerical values 
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONVPAR_SHAL   : contains deep convection constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module MODD_CONVPAR_SHAL, routine INI_CONVPAR)
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
!!                  05/05/09 E. Bazile
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAR_SHAL
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
IF (LHOOK) CALL DR_HOOK('INI_CONVPAR_SHAL',0,ZHOOK_HANDLE)
XA25     = 625.E6    ! 25 km x 25 km reference grid area
!
XCRAD       = 50.    ! cloud radius 
XCTIME_SHAL = 10800. ! convective adjustment time
XCDEPTH     = 0.5E3  ! minimum necessary shallow cloud depth
XCDEPTH_D   = 2.5E3  ! maximum allowed shallow cloud depth
XDTPERT     = .2     ! add small Temp perturbation at LCL
XATPERT     = 0.     ! 0.=original scheme , recommended = 1000. 
XBTPERT     = 1.     ! 1.=original scheme , recommended = 0.
!
XENTR    = 0.02      ! entrainment constant (m/Pa) = 0.2 (m)  
!
XZLCL    = 0.5E3     ! maximum allowed allowed height 
                     ! difference between the DPL and the surface
XZPBL    = 40.E2     ! minimum mixed layer depth to sustain convection
!
!
XNHGAM   = 1.3333    ! accounts for non-hydrost. pressure 
                     ! in buoyancy term of w equation
                     ! = 2 / (1+gamma)
XTFRZ1   = 268.16    ! begin of freezing interval
XTFRZ2   = 248.16    ! end of freezing interval
!

XSTABT   = 0.75      ! factor to assure stability in  fractional time
                     ! integration, routine CONVECT_CLOSURE
XSTABC   = 0.95      ! factor to assure stability in CAPE adjustment,
                     !  routine CONVECT_CLOSURE
XAW      = 0.        ! 0.= Original scheme , 1 = recommended 
XBW      = 1.        ! 1.= Original scheme , 0 = recommended
LLSMOOTH = .TRUE.
!
!
IF (LHOOK) CALL DR_HOOK('INI_CONVPAR_SHAL',1,ZHOOK_HANDLE)
END SUBROUTINE INI_CONVPAR_SHAL 
