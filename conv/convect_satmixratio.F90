!     ######spl
      SUBROUTINE CONVECT_SATMIXRATIO(CST, D, PPRES, PT, PEW, PLV, PLS, PCPH)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
!------------------------- ------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : CST_T
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
TYPE(CST_T),            INTENT(IN) :: CST
TYPE(DIMPHYEX_T),       INTENT(IN) :: D
REAL, DIMENSION(D%NIT),  INTENT(IN) :: PPRES   ! pressure
REAL, DIMENSION(D%NIT),  INTENT(IN) :: PT      ! temperature   
!
REAL, DIMENSION(D%NIT),  INTENT(OUT):: PEW     ! vapor saturation mixing ratio
REAL, DIMENSION(D%NIT),  INTENT(OUT):: PLV     ! latent heat L_v    
REAL, DIMENSION(D%NIT),  INTENT(OUT):: PLS     ! latent heat L_s  
REAL, DIMENSION(D%NIT),  INTENT(OUT):: PCPH    ! specific heat C_ph   
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(D%NIT)              :: ZT      ! temperature   
REAL    :: ZEPS           ! R_d / R_v
!
!
!-------------------------------------------------------------------------------
!
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('CONVECT_SATMIXRATIO',0,ZHOOK_HANDLE)
    ZEPS      = CST%XRD / CST%XRV
!
    ZT(D%NIB:D%NIE)     = MIN( 400., MAX( PT(D%NIB:D%NIE), 10. ) ) ! overflow bound
    PEW(D%NIB:D%NIE)    = EXP( CST%XALPW - CST%XBETAW / ZT(D%NIB:D%NIE) - CST%XGAMW * ALOG( ZT(D%NIB:D%NIE) ) )
    PEW(D%NIB:D%NIE)    = ZEPS * PEW(D%NIB:D%NIE) / ( PPRES(D%NIB:D%NIE) - PEW(D%NIB:D%NIE) )
!
    PLV(D%NIB:D%NIE)    = CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZT(D%NIB:D%NIE) - CST%XTT ) ! compute L_v
    PLS(D%NIB:D%NIE)    = CST%XLSTT + ( CST%XCPV - CST%XCI ) * ( ZT(D%NIB:D%NIE) - CST%XTT ) ! compute L_i
!    
    PCPH(D%NIB:D%NIE)   = CST%XCPD + CST%XCPV * PEW(D%NIB:D%NIE)                     ! compute C_ph 
!
IF (LHOOK) CALL DR_HOOK('CONVECT_SATMIXRATIO',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECT_SATMIXRATIO
