!     ######spl
      SUBROUTINE CONVECT_CONDENS( KLON,                                           &
                                  KICE, PPRES, PTHL, PRW, PRCO, PRIO, PZ, OWORK1, &
                                  PT, PEW, PRC, PRI, PLV, PLS, PCPH   )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #############################################################################
!
!!**** Compute temperature cloud and ice water content from enthalpy and r_w 
!!
!!
!!    PURPOSE
!!    -------
!!     The purpose of this routine is to determine cloud condensate
!!     and to return values for L_v, L_s and C_ph
!!
!!
!!**  METHOD
!!    ------
!!     Condensate is extracted iteratively 
!!     
!!
!!    EXTERNAL
!!    --------
!!     None
!!     
!!
!!    IMPLICIT ARGUMENTS     
!!    ------------------
!!
!!      Module MODD_CST
!!          XG                   ! gravity constant
!!          XALPW, XBETAW, XGAMW ! constants for water saturation pressure
!!          XALPI, XBETAI, XGAMI ! constants for ice saturation pressure
!!          XP00                 ! reference pressure
!!          XRD, XRV             ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV           ! specific heat for dry air and water vapor
!!          XCL, XCI             ! specific heat for liquid water and ice
!!          XTT                  ! triple point temperature
!!          XLVTT, XLSTT         ! vaporization, sublimation heat constant
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONVPAR
!!          XTFRZ1               ! begin of freezing interval
!!          XTFRZ2               ! end of freezing interval
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CONDENS)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : XALPI, XALPW, XBETAI, XBETAW, XCI, XCL, XCPD, XCPV, XG, XGAMI, XGAMW, XLSTT, XLVTT, XRD, XRV, XTT
USE MODD_CONVPAR, ONLY : XTFRZ1, XTFRZ2
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN)                :: KLON    ! horizontal loop index
INTEGER, INTENT(IN)                :: KICE    ! flag for ice ( 1 = yes,
                                              !                0 = no ice )
REAL, DIMENSION(KLON),   INTENT(IN) :: PPRES  ! pressure
REAL, DIMENSION(KLON),   INTENT(IN) :: PTHL   ! enthalpy (J/kg)
REAL, DIMENSION(KLON),   INTENT(IN) :: PRW    ! total water mixing ratio  
REAL, DIMENSION(KLON),   INTENT(IN) :: PRCO   ! cloud water estimate (kg/kg)
REAL, DIMENSION(KLON),   INTENT(IN) :: PRIO   ! cloud ice   estimate (kg/kg)
REAL, DIMENSION(KLON),   INTENT(IN) :: PZ     ! level height (m)
LOGICAL, DIMENSION(KLON),INTENT(IN) :: OWORK1 ! logical mask         
!
!
REAL, DIMENSION(KLON),   INTENT(OUT):: PT     ! temperature   
REAL, DIMENSION(KLON),   INTENT(OUT):: PRC    ! cloud water mixing ratio(kg/kg)
REAL, DIMENSION(KLON),   INTENT(OUT):: PRI    ! cloud ice mixing ratio  (kg/kg)
REAL, DIMENSION(KLON),   INTENT(OUT):: PLV    ! latent heat L_v    
REAL, DIMENSION(KLON),   INTENT(OUT):: PLS    ! latent heat L_s  
REAL, DIMENSION(KLON),   INTENT(OUT):: PCPH   ! specific heat C_ph   
REAL, DIMENSION(KLON),   INTENT(OUT):: PEW    ! water saturation mixing ratio  
!
!*       0.2   Declarations of local variables KLON
!
INTEGER :: JITER          ! iteration index
REAL    :: ZEPS           ! R_d / R_v
!
REAL, DIMENSION(KLON)    :: ZEI           ! ice saturation mixing ratio
REAL, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3, ZT ! work arrays
!
!
!-------------------------------------------------------------------------------
!
!*       1.     Initialize temperature and Exner function
!               -----------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CONVECT_CONDENS',0,ZHOOK_HANDLE)
ZEPS        = XRD / XRV
!
!
    ! Make a first temperature estimate, based e.g. on values of
    !  r_c and r_i at lower level
!
      !! Note that the definition of ZCPH is not the same as used in
      !! routine CONVECT_SATMIXRATIO
     PCPH(:)   = XCPD + XCPV * PRW(:)
     ZWORK1(:) = ( 1. + PRW(:) ) * XG * PZ(:)
     PT(:)     = ( PTHL(:) + PRCO(:) * XLVTT + PRIO(:) * XLSTT - ZWORK1(:) )   &
                 / PCPH(:)
     PT(:)     = MAX(180., MIN( 330., PT(:) ) ) ! set overflow bounds in
                                                    ! case that PTHL=0     
!
!
!*       2.     Enter the iteration loop
!               ------------------------
!    
DO JITER = 1,6
     PEW(:) = EXP( XALPW - XBETAW / PT(:) - XGAMW * ALOG( PT(:) ) )
     ZEI(:) = EXP( XALPI - XBETAI / PT(:) - XGAMI * ALOG( PT(:) ) )
     PEW(:) = ZEPS * PEW(:) / ( PPRES(:) - PEW(:) )
     ZEI(:) = ZEPS * ZEI(:) / ( PPRES(:) - ZEI(:) )    
!
     PLV(:)    = XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT ) ! compute L_v
     PLS(:)    = XLSTT + ( XCPV - XCI ) * ( PT(:) - XTT ) ! compute L_i
!    
     ZWORK2(:) = ( XTFRZ1 - PT(:) ) / ( XTFRZ1 - XTFRZ2 ) ! freezing interval
     ZWORK2(:) = MAX( 0., MIN(1., ZWORK2(:) ) ) * REAL( KICE )
     ZWORK3(:) = ( 1. - ZWORK2(:) ) * PEW(:) + ZWORK2(:) * ZEI(:)
     PRC(:)    = MAX( 0., ( 1. - ZWORK2(:) ) * ( PRW(:) - ZWORK3(:) ) )
     PRI(:)    = MAX( 0.,  ZWORK2(:) * ( PRW(:) - ZWORK3(:) ) )
     ZT(:)     = ( PTHL(:) + PRC(:) * PLV(:) + PRI(:) * PLS(:) - ZWORK1(:) )   &
                 / PCPH(:)
     PT(:) = PT(:) + ( ZT(:) - PT(:) ) * 0.4  ! force convergence
     PT(:) = MAX( 175., MIN( 330., PT(:) ) )
END DO
!
!
IF (LHOOK) CALL DR_HOOK('CONVECT_CONDENS',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECT_CONDENS
