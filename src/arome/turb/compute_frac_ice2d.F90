!     ######spl
      SUBROUTINE COMPUTE_FRAC_ICE2D(HFRAC_ICE,PFRAC_ICE,PT)
!    ##########################################################
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
!!
!! --------------------------------------------------------------------------
!       0. DECLARATIONS
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
CHARACTER*1         , INTENT(IN)    :: HFRAC_ICE ! scheme to use
REAL, DIMENSION(:,:), INTENT(IN)    :: PT        ! Temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PFRAC_ICE ! Ice fraction (1 for ice only, 0 for liquid only)
!-------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
INTEGER :: JK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------
!
!       0.3  Initialisation
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE2D',0,ZHOOK_HANDLE)
!
!----------------------------------------------------------------------------
!
!       1 Compute FRAC_ICE
!         ----------------
!
DO JK=1, SIZE(PT,2)
  CALL COMPUTE_FRAC_ICE1D(HFRAC_ICE,PFRAC_ICE(:,JK),PT(:,JK))
ENDDO

IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE2D',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_FRAC_ICE2D
