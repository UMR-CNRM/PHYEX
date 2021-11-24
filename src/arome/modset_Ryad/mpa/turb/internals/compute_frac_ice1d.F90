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
!!      R. El Khatib 24-Aug-2021 Optimization by cache re-use + assume data is contiguous
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
REAL, CONTIGUOUS, DIMENSION(:), INTENT(IN)    :: PT         ! temperature
REAL, CONTIGUOUS, DIMENSION(:), INTENT(INOUT) :: PFRAC_ICE  ! Ice fraction (1 for ice only, 0 for liquid only)
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
SELECT CASE(HFRAC_ICE)
  CASE ('T') !using Temperature
    PFRAC_ICE(:) = MAX( 0., MIN(1., (( XTMAXMIX - PT(:) ) / ( XTMAXMIX - XTMINMIX )) ) ) ! freezing interval
  CASE ('O') !using Temperature with old formulae
    PFRAC_ICE(:) = MAX( 0., MIN(1., (( XTT - PT(:) ) / 40.) ) ) ! freezing interval
  CASE ('N') !No ice
    PFRAC_ICE(:) = 0.
  CASE ('S') !Same as previous
    ! (almost) nothing to do
    PFRAC_ICE(:) = MAX( 0., MIN(1., PFRAC_ICE(:) ) )
  CASE DEFAULT
    WRITE(*,*) ' STOP'
    WRITE(*,*) ' INVALID OPTION IN COMPUTE_FRAC_ICE, HFRAC_ICE=',HFRAC_ICE
    CALL ABORT
    STOP
END SELECT

IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE1D',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_FRAC_ICE1D
