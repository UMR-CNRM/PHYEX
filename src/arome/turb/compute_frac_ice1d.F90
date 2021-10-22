!     ######spl
      SUBROUTINE COMPUTE_FRAC_ICE1D(HFRAC_ICE,PFRAC_ICE,PT)
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
!!      S. Riette        08/2016 add option O
!!
!! --------------------------------------------------------------------------
!       0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_NEB, ONLY : XTMINMIX, XTMAXMIX
USE MODD_CST, ONLY : XTT
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
CHARACTER*1       , INTENT(IN)    :: HFRAC_ICE  ! scheme to use
REAL, DIMENSION(:), INTENT(IN)    :: PT         ! temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PFRAC_ICE  ! Ice fraction (1 for ice only, 0 for liquid only)
!
!               0.2  declaration of local variables
! 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!               0.2  initialisation
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE1D',0,ZHOOK_HANDLE)
!
!------------------------------------------------------------------------
!                1. Compute FRAC_ICE
!
IF (HFRAC_ICE=='T') THEN !using Temperature
  PFRAC_ICE(:) = ( XTMAXMIX - PT(:) ) / ( XTMAXMIX - XTMINMIX ) ! freezing interval
ELSEIF (HFRAC_ICE=='O') THEN !using Temperature with old formulae
  PFRAC_ICE(:) = ( XTT - PT(:) ) / 40. ! freezing interval
ELSEIF (HFRAC_ICE=='N') THEN !No ice
  PFRAC_ICE(:) = 0.
ELSEIF (HFRAC_ICE=='S') THEN !Same as previous
  !nothing to do
ELSE
  WRITE(*,*) ' STOP'
  WRITE(*,*) ' INVALID OPTION IN COMPUTE_FRAC_ICE, HFRAC_ICE=',HFRAC_ICE
  CALL ABORT
  STOP
ENDIF

PFRAC_ICE(:) = MAX( 0., MIN(1., PFRAC_ICE(:) ) )

IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE1D',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_FRAC_ICE1D
