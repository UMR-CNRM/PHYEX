!MNH_LIC Copyright 1995-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 25/02/2019: split rain_ice (cleaner and easier to maintain/debug)
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet 29/05/2019: remove PACK/UNPACK intrinsics (to get more performance and better OpenACC support)
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_NUCLEATION

  IMPLICIT NONE

  PRIVATE

  PUBLIC RAIN_ICE_NUCLEATION

CONTAINS

SUBROUTINE RAIN_ICE_NUCLEATION(KIB, KIE, KJB, KJE, KKTB, KKTE,KRR,PTSTEP,&
     PTHT,PPABST,PRHODJ,PRHODREF,PRVT,PRCT,PRRT,PRIT,PRST,PRGT,&
     PCIT,PEXNREF,PTHS,PRVS,PRIS,PT,PRHT)
!
!*      0. DECLARATIONS
!          ------------
!
use modd_budget,          only: lbudget_th, lbudget_rv, lbudget_ri, &
                                NBUDGET_TH, NBUDGET_RV, NBUDGET_RI, &
                                tbudgets
use MODD_CST,             only: XALPI, XALPW, XBETAI, XBETAW, XCI, XCL, XCPD, XCPV, XGAMI, XGAMW, &
                                XLSTT, XMD, XMV, XP00, XRD, XTT
use MODD_RAIN_ICE_PARAM,  only: XALPHA1, XALPHA2, XBETA1, XBETA2, XMNU0, XNU10, XNU20

use mode_budget,          only: Budget_store_init, Budget_store_end
use mode_tools,           only: Countjv

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KIB, KIE, KJB, KJE, KKTB, KKTE
INTEGER,                          INTENT(IN)    :: KRR     ! Number of moist variable
REAL,                             INTENT(IN)    :: PTSTEP  ! Double Time step
                                                           ! (single if cold start)
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PTHT    ! Theta at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),           INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),           INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),           INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),           INTENT(OUT)   :: PT      ! Temperature
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
!
!*       0.2  declaration of local variables
!
INTEGER                            :: INEGT
INTEGER                            :: JL       ! and PACK intrinsics
INTEGER, DIMENSION(SIZE(PEXNREF))  :: I1,I2,I3 ! Used to replace the COUNT
LOGICAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) &
                                   :: GNEGT    ! Test where to compute the HEN process
REAL,    DIMENSION(:), ALLOCATABLE :: ZRVT     ! Water vapor m.r. at t
REAL,    DIMENSION(:), ALLOCATABLE :: ZCIT     ! Pristine ice conc. at t
REAL,    DIMENSION(:), ALLOCATABLE :: ZZT,   & ! Temperature
                                      ZPRES, & ! Pressure
                                      ZZW,   & ! Work array
                                      ZUSW,  & ! Undersaturation over water
                                      ZSSI     ! Supersaturation over ice
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
                                    :: ZW      ! work array
!
!-------------------------------------------------------------------------------

if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HENU', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'HENU', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'HENU', pris(:, :, :) * prhodj(:, :, :) )
!
!  compute the temperature and the pressure
!
PT(:,:,:) = PTHT(:,:,:) * ( PPABST(:,:,:) / XP00 ) ** (XRD/XCPD)
!
!  optimization by looking for locations where
!  the temperature is negative only !!!
!
GNEGT(:,:,:) = .FALSE.
GNEGT(KIB:KIE,KJB:KJE,KKTB:KKTE) = PT(KIB:KIE,KJB:KJE,KKTB:KKTE)<XTT
INEGT = COUNTJV( GNEGT(:,:,:),I1(:),I2(:),I3(:))
IF( INEGT >= 1 ) THEN
  ALLOCATE(ZRVT(INEGT)) ;
  ALLOCATE(ZCIT(INEGT)) ;
  ALLOCATE(ZZT(INEGT))  ;
  ALLOCATE(ZPRES(INEGT));
  DO JL=1,INEGT
    ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
    ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
    ZZT(JL) = PT(I1(JL),I2(JL),I3(JL))
    ZPRES(JL) = PPABST(I1(JL),I2(JL),I3(JL))
  ENDDO
  ALLOCATE(ZZW(INEGT))
  ALLOCATE(ZUSW(INEGT))
  ALLOCATE(ZSSI(INEGT))
    ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) )           ! es_i
    ZZW(:) = MIN(ZPRES(:)/2., ZZW(:))             ! safety limitation
    ZSSI(:) = ZRVT(:)*( ZPRES(:)-ZZW(:) ) / ( (XMV/XMD) * ZZW(:) ) - 1.0
                                                      ! Supersaturation over ice
    ZUSW(:) = EXP( XALPW - XBETAW/ZZT(:) - XGAMW*ALOG(ZZT(:) ) )          ! es_w
    ZUSW(:) = MIN(ZPRES(:)/2.,ZUSW(:))            ! safety limitation
    ZUSW(:) = ( ZUSW(:)/ZZW(:) )*( (ZPRES(:)-ZZW(:))/(ZPRES(:)-ZUSW(:)) ) - 1.0
                             ! Supersaturation of saturated water vapor over ice
!
!*       3.1     compute the heterogeneous nucleation source: RVHENI
!
!*       3.1.1   compute the cloud ice concentration
!
  ZZW(:) = 0.0
  ZSSI(:) = MIN( ZSSI(:), ZUSW(:) ) ! limitation of SSi according to SSw=0
  WHERE( (ZZT(:)<XTT-5.0) .AND. (ZSSI(:)>0.0) )
    ZZW(:) = XNU20 * EXP( XALPHA2*ZSSI(:)-XBETA2 )
  END WHERE
  WHERE( (ZZT(:)<=XTT-2.0) .AND. (ZZT(:)>=XTT-5.0) .AND. (ZSSI(:)>0.0) )
    ZZW(:) = MAX( XNU20 * EXP( -XBETA2 ),XNU10 * EXP( -XBETA1*(ZZT(:)-XTT) ) * &
                               ( ZSSI(:)/ZUSW(:) )**XALPHA1 )
  END WHERE
  ZZW(:) = ZZW(:) - ZCIT(:)
  IF( MAXVAL(ZZW(:)) > 0.0 ) THEN
!
!*       3.1.2   update the r_i and r_v mixing ratios
!
    ZZW(:) = MIN( ZZW(:),50.E3 ) ! limitation provisoire a 50 l^-1
    ZW(:,:,:) = 0.0
    DO JL=1, INEGT
      ZW(I1(JL), I2(JL), I3(JL)) = ZZW( JL )
    END DO
    ZW(:,:,:) = MAX( ZW(:,:,:) ,0.0 ) *XMNU0/(PRHODREF(:,:,:)*PTSTEP)
    PRIS(:,:,:) = PRIS(:,:,:) + ZW(:,:,:)
    PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
    IF ( KRR == 7 ) THEN
        PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(XLSTT+(XCPV-XCI)*(PT(:,:,:)-XTT))   &
                 /( (XCPD + XCPV*PRVT(:,:,:) + XCL*(PRCT(:,:,:)+PRRT(:,:,:))   &
       + XCI*(PRIT(:,:,:)+PRST(:,:,:)+PRGT(:,:,:)+PRHT(:,:,:)))*PEXNREF(:,:,:) )
      ELSE IF( KRR == 6 ) THEN
        PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(XLSTT+(XCPV-XCI)*(PT(:,:,:)-XTT))   &
                 /( (XCPD + XCPV*PRVT(:,:,:) + XCL*(PRCT(:,:,:)+PRRT(:,:,:))   &
                   + XCI*(PRIT(:,:,:)+PRST(:,:,:)+PRGT(:,:,:)))*PEXNREF(:,:,:) )
    END IF
                                 ! f(L_s*(RVHENI))
    ZZW(:) = MAX( ZZW(:)+ZCIT(:),ZCIT(:) )
    PCIT(:,:,:) = MAX( PCIT(:,:,:), 0.0 )
    DO JL=1, INEGT
      PCIT(I1(JL), I2(JL), I3(JL)) = MAX( ZZW( JL ), PCIT(I1(JL), I2(JL), I3(JL)), 0.0 )
    END DO
  END IF
  DEALLOCATE(ZSSI)
  DEALLOCATE(ZUSW)
  DEALLOCATE(ZZW)
  DEALLOCATE(ZPRES)
  DEALLOCATE(ZZT)
  DEALLOCATE(ZCIT)
  DEALLOCATE(ZRVT)
END IF
!
!*       3.1.3   budget storage
!
if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HENU', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'HENU', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'HENU', pris(:, :, :) * prhodj(:, :, :) )

END SUBROUTINE RAIN_ICE_NUCLEATION

END MODULE MODE_RAIN_ICE_NUCLEATION
