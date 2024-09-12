      ELEMENTAL SUBROUTINE CONVECT_SATMIXRATIO(PPRES, PT, PEPS, PEW, PLV, PLS, PCPH)

! ******* TO BE INCLUDED IN THE *CONTAINS* OF A SUBROUTINE, IN ORDER TO EASE AUTOMATIC INLINING ******
! => Don't use drHook !!!

!     ################################################################
!
!!**** Compute vapor saturation mixing ratio over liquid water
!!
!!
!!    PDRPOSE
!!    -------
!!     The purpose of this routine is to determine saturation mixing ratio
!!     and to return values for L_v L_s and C_ph
!!
!!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!     None
!!     
!!
!!    IMPLICIT ARGUMENTS    
!!    ------------------
!!      Module MODD_CST
!!          XALPW, XBETAW, XGAMW ! constants for water saturation pressure
!!          XRD, XRV             ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV           ! specific heat for dry air and water vapor
!!          XCL, XCI             ! specific heat for liquid water and ice
!!          XTT                  ! triple point temperature
!!          XLVTT, XLSTT         ! vaporization, sublimation heat constant
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_SATMIXRATIO)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!!      R. El Khatib     01-Jun-2023 written as a include file
!------------------------- ------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL, INTENT(IN) :: PPRES   ! pressure
REAL, INTENT(IN) :: PT      ! temperature   
REAL, INTENT(IN) :: PEPS    ! CST%XRD / CST%XRV (ideally pre-computed in CST%)
!
REAL, INTENT(OUT):: PEW     ! vapor saturation mixing ratio
REAL, INTENT(OUT):: PLV     ! latent heat L_v    
REAL, INTENT(OUT):: PLS     ! latent heat L_s  
REAL, INTENT(OUT):: PCPH    ! specific heat C_ph   
!
!*       0.2   Declarations of local variables :
!
REAL :: ZT      ! temperature   
!
!-------------------------------------------------------------------------------
!
    ZT     = MIN( 400., MAX( PT, 10. ) ) ! overflow bound
    PEW    = EXP( CST%XALPW - CST%XBETAW / ZT - CST%XGAMW * ALOG( ZT ) )
    PEW    = PEPS * PEW / ( PPRES - PEW )
!
    PLV    = CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZT - CST%XTT ) ! compute L_v
    PLS    = CST%XLSTT + ( CST%XCPV - CST%XCI ) * ( ZT - CST%XTT ) ! compute L_i
!    
    PCPH   = CST%XCPD + CST%XCPV * PEW                     ! compute C_ph 
!
END SUBROUTINE CONVECT_SATMIXRATIO
