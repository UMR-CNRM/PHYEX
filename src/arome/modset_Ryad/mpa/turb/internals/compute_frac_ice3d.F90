!     ######spl
      SUBROUTINE COMPUTE_FRAC_ICE3D(HFRAC_ICE,PFRAC_ICE,PT)
!     #################################################################
!
!
!!****  *COMPUTE_FRAC_ICE* - computes ice fraction
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Julien PERGAUD      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         13/03/06
!!      S. Riette        April 2011 optimisation
!!      R. El Khatib 24-Aug-2021 Optimize by loop collapsing  + assume data is contiguous
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODI_COMPUTE_FRAC_ICE1D
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
CHARACTER*1           , INTENT(IN)    :: HFRAC_ICE ! scheme to use
REAL, TARGET, CONTIGUOUS, DIMENSION(:,:,:), INTENT(IN)    :: PT        ! Temperature
REAL, TARGET, CONTIGUOUS, DIMENSION(:,:,:), INTENT(INOUT) :: PFRAC_ICE ! Ice fraction (1 for ice only, 0 for liquid only)
!-------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
REAL, POINTER, DIMENSION(:) :: ZT, ZFRAC_ICE

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------
!
!       0.3  Initialisation
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE3D',0,ZHOOK_HANDLE)
!
!----------------------------------------------------------------------------
!
!       1 Compute FRAC_ICE
!         ----------------
!
ZT(1:SIZE(PT))=>PT
ZFRAC_ICE(1:SIZE(PFRAC_ICE))=>PFRAC_ICE
CALL COMPUTE_FRAC_ICE1D(HFRAC_ICE,ZFRAC_ICE,ZT)

IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE3D',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_FRAC_ICE3D
