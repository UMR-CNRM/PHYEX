!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_CONVECT_CONDENS
IMPLICIT NONE
CONTAINS
      SUBROUTINE CONVECT_CONDENS( CST, D, CONVPAR,                        &
                                  PICE, PEPS0, PPRES, PTHL, PRW, PRCO, PRIO, PZ, &
                                  PT, PEW, PRC, PRI, PLV, PLS, PCPH   )

!$ACDC singlecolumn

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
!!     R. El Khatib 16-Jul-2025 Optimization
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODD_CST, ONLY : CST_t
USE MODD_CONVPAR, ONLY : CONVPAR_t
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN) :: CST
TYPE(DIMPHYEX_t),         INTENT(IN) :: D
TYPE(CONVPAR_t),          INTENT(IN) :: CONVPAR
REAL,                     INTENT(IN) :: PICE   ! flag for ice ( 1 = yes, 0 = no ice )
REAL,                     INTENT(IN) :: PEPS0  ! R_d / R_v
REAL, DIMENSION(D%NIJT),   INTENT(IN) :: PPRES  ! pressure
REAL, DIMENSION(D%NIJT),   INTENT(IN) :: PTHL   ! enthalpy (J/kg)
REAL, DIMENSION(D%NIJT),   INTENT(IN) :: PRW    ! total water mixing ratio  
REAL, DIMENSION(D%NIJT),   INTENT(IN) :: PRCO   ! cloud water estimate (kg/kg)
REAL, DIMENSION(D%NIJT),   INTENT(IN) :: PRIO   ! cloud ice   estimate (kg/kg)
REAL, DIMENSION(D%NIJT),   INTENT(IN) :: PZ     ! level height (m)
!
!
REAL, DIMENSION(D%NIJT),   INTENT(OUT):: PT     ! temperature   
REAL, DIMENSION(D%NIJT),   INTENT(OUT):: PEW    ! water saturation mixing ratio  
REAL, DIMENSION(D%NIJT),   INTENT(OUT):: PRC    ! cloud water mixing ratio(kg/kg)
REAL, DIMENSION(D%NIJT),   INTENT(OUT):: PRI    ! cloud ice mixing ratio  (kg/kg)
REAL, DIMENSION(D%NIJT),   INTENT(OUT):: PLV    ! latent heat L_v    
REAL, DIMENSION(D%NIJT),   INTENT(OUT):: PLS    ! latent heat L_s  
REAL, DIMENSION(D%NIJT),   INTENT(OUT):: PCPH   ! specific heat C_ph   
!
!*       0.2   Declarations of local variables D%NIJT
!
INTEGER :: JITER, JI      ! iteration index
!
REAL :: ZEI           ! ice saturation mixing ratio
REAL, DIMENSION(D%NIJT)    :: ZWORK1 ! work arrays
REAL :: ZWORK2, ZTMP, ZXTFRDZ, ZCPVMCL, ZCPVMCI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!
!-------------------------------------------------------------------------------
!
!*       1.     Initialize temperature and Exner function
!               -----------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CONVECT_CONDENS',0,ZHOOK_HANDLE)
!
!
    ! Make a first temperature estimate, based e.g. on values of
    !  r_c and r_i at lower level
!
      !! Note that the definition of ZCPH is not the same as used in
      !! routine CONVECT_SATMIXRATIO
  DO JI=D%NIJB,D%NIJE
     PCPH(JI)   = CST%XCPD + CST%XCPV * PRW(JI)
     ZWORK1(JI) = ( 1. + PRW(JI) ) * CST%XG * PZ(JI)
     PT(JI)     = MAX(180., MIN( 330., &
                 & ( PTHL(JI) + PRCO(JI) * CST%XLVTT + PRIO(JI) * CST%XLSTT - ZWORK1(JI) ) / PCPH(JI) &
                 & ) ) ! set overflow bounds in case that PTHL=0     
  ENDDO
!
!
!*       2.     Enter the iteration loop
!               ------------------------
!
ZXTFRDZ=CONVPAR%XTFRZ1 - CONVPAR%XTFRZ2
ZCPVMCL=CST%XCPV-CST%XCL
ZCPVMCI=CST%XCPV-CST%XCI
DO JITER = 1,6
  DO JI=D%NIJB,D%NIJE
    ZTMP = LOG( PT(JI) )
    PEW(JI) = EXP( CST%XALPW - CST%XBETAW / PT(JI) - CST%XGAMW * ZTMP )
    PEW(JI) = PEPS0 * PEW(JI) / ( PPRES(JI) - PEW(JI) )
    ZEI = EXP( CST%XALPI - CST%XBETAI / PT(JI) - CST%XGAMI * ZTMP )
    ZEI = PEPS0 * ZEI / ( PPRES(JI) - ZEI )    
    !
    ZWORK2 = MAX( 0., MIN(1., ( CONVPAR%XTFRZ1 - PT(JI) ) / ZXTFRDZ ) ) * PICE ! freezing interval
    ZTMP = PRW(JI) - ( ( 1. - ZWORK2 ) * PEW(JI) + ZWORK2 * ZEI )
    PRC(JI)    = MAX( 0., ( 1. - ZWORK2 ) * ZTMP )
    PRI(JI)    = MAX( 0.,  ZWORK2 * ZTMP )

    ZTMP = PT(JI) - CST%XTT
    PLV(JI)    = CST%XLVTT + ZCPVMCL * ZTMP ! compute L_v
    PLS(JI)    = CST%XLSTT + ZCPVMCI * ZTMP ! compute L_i

    ZTMP     = ( PTHL(JI) + PRC(JI) * PLV(JI) + PRI(JI) * PLS(JI) - ZWORK1(JI) ) / PCPH(JI)
    PT(JI) = MAX( 175., MIN( 330., PT(JI) + ( ZTMP - PT(JI) ) * 0.4 ) ) ! force convergence
  END DO
END DO
!
!
IF (LHOOK) CALL DR_HOOK('CONVECT_CONDENS',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECT_CONDENS
END MODULE MODE_CONVECT_CONDENS
