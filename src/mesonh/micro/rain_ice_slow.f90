!MNH_LIC Copyright 1995-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 25/02/2019: split rain_ice (cleaner and easier to maintain/debug)
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_SLOW

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: RAIN_ICE_SLOW

CONTAINS

SUBROUTINE RAIN_ICE_SLOW(OMICRO, PINVTSTEP, PRHODREF,                                      &
                         PRCT, PRRT, PRIT, PRST, PRGT, PRHODJ, PZT, PPRES,                 &
                         PLSFACT, PLVFACT, PSSI,                                           &
                         PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, PTHS,                         &
                         PAI, PCJ, PKA, PDV, PLBDAS, PLBDAG)
!
!*      0. DECLARATIONS
!          ------------
!
use modd_budget,         only: lbudget_th, lbudget_rv, lbudget_rc, lbudget_rr, lbudget_ri, lbudget_rs, lbudget_rg, &
                               NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, NBUDGET_RS, NBUDGET_RG, &
                               tbudgets
USE MODD_CST,            only: XALPI, XBETAI, XCI, XCPV, XGAMI, XLSTT, XMNH_HUGE_12_LOG, XP00, XRV, XTT
USE MODD_RAIN_ICE_DESCR_n, only: XCEXVT, XLBDAS_MAX, XLBEXG, XLBEXS, XLBG, XLBS, XRTMIN, XBS
USE MODD_RAIN_ICE_PARAM_n, only: X0DEPG, X0DEPS, X1DEPG, X1DEPS, XACRIAUTI, XALPHA3, XBCRIAUTI, XBETA3, XCOLEXIS, XCRIAUTI, &
                               XEX0DEPG, XEX0DEPS, XEX1DEPG, XEX1DEPS, XEXIAGGS, XFIAGGS, XHON, XSCFAC, XTEXAUTI, XTIMAUTI

use mode_budget,         only: Budget_store_add
#ifdef MNH_OPENACC
USE MODE_MNH_ZWORK,      ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
#endif
use mode_mppdb
#if defined(MNH_BITREP) || defined(MNH_BITREP_OMP)
USE MODI_BITREP
#endif
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,  DIMENSION(:,:,:), intent(in)    :: OMICRO   ! Test where to compute all processes
REAL,                       intent(in)    :: PINVTSTEP
REAL,     DIMENSION(:),     intent(in)    :: PRHODREF ! RHO Dry REFerence
REAL,     DIMENSION(:),     intent(in)    :: PRCT     ! Cloud water m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRRT     ! Rain water m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRIT     ! Pristine ice m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRST     ! Snow/aggregate m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRGT     ! Graupel m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRHODJ   ! RHO times Jacobian
REAL,     DIMENSION(:),     intent(in)    :: PZT      ! Temperature
REAL,     DIMENSION(:),     intent(in)    :: PPRES    ! Pressure
REAL,     DIMENSION(:),     intent(in)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL,     DIMENSION(:),     intent(in)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL,     DIMENSION(:),     intent(in)    :: PSSI     ! Supersaturation over ice
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRVS     ! Water vapor m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRCS     ! Cloud water m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRRS     ! Rain water m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRIS     ! Pristine ice m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRSS     ! Snow/aggregate m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRGS     ! Graupel m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PTHS     ! Theta source
REAL,     DIMENSION(:),     intent(OUT)   :: PAI      ! Thermodynamical function
REAL,     DIMENSION(:),     intent(OUT)   :: PCJ      ! Function to compute the ventilation coefficient
REAL,     DIMENSION(:),     intent(OUT)   :: PKA      ! Thermal conductivity of the air
REAL,     DIMENSION(:),     intent(OUT)   :: PDV      ! Diffusivity of water vapor in the air
REAL,     DIMENSION(:),     intent(OUT)   :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL,     DIMENSION(:),     intent(OUT)   :: PLBDAG   ! Slope parameter of the graupel   distribution
!
!*       0.2  declaration of local variables
!
#ifndef MNH_OPENACC
LOGICAL, DIMENSION(:), ALLOCATABLE :: GWORK
REAL,    DIMENSION(:), ALLOCATABLE :: ZCRIAUTI ! Snow-to-ice autoconversion thres.
REAL,    DIMENSION(:), ALLOCATABLE :: ZZW      ! Work array
real,    dimension(:), ALLOCATABLE :: zz_diff
#else
LOGICAL, DIMENSION(:), pointer, contiguous :: GWORK
REAL,    DIMENSION(:), pointer, contiguous :: ZCRIAUTI ! Snow-to-ice autoconversion thres.
REAL,    DIMENSION(:), pointer, contiguous :: ZZW      ! Work array
real,    dimension(:), pointer, contiguous :: zz_diff
#endif
!
INTEGER                            :: JL,JLU
!-------------------------------------------------------------------------------
!
! IN variables
!
!$acc data present( OMICRO, PRHODREF, PRCT, PRRT, PRIT, PRST, PRGT, &
!$acc &             PRHODJ, PZT, PPRES, PLSFACT, PLVFACT, PSSI,     &
!
! INOUT variables
!
!$acc &             PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, PTHS,       &
!
! OUT variables
!
!$acc &             PAI, PCJ, PKA, PDV, PLBDAS, PLBDAG )

IF (MPPDB_INITIALIZED) THEN
  !Check all IN arrays
  CALL MPPDB_CHECK(OMICRO,"RAIN_ICE_SLOW beg:OMICRO")
  CALL MPPDB_CHECK(PRHODREF,"RAIN_ICE_SLOW beg:PRHODREF")
  CALL MPPDB_CHECK(PRCT,"RAIN_ICE_SLOW beg:PRCT")
  CALL MPPDB_CHECK(PRRT,"RAIN_ICE_SLOW beg:PRRT")
  CALL MPPDB_CHECK(PRIT,"RAIN_ICE_SLOW beg:PRIT")
  CALL MPPDB_CHECK(PRST,"RAIN_ICE_SLOW beg:PRST")
  CALL MPPDB_CHECK(PRGT,"RAIN_ICE_SLOW beg:PRGT")
  CALL MPPDB_CHECK(PRHODJ,"RAIN_ICE_SLOW beg:PRHODJ")
  CALL MPPDB_CHECK(PZT,"RAIN_ICE_SLOW beg:PZT")
  CALL MPPDB_CHECK(PPRES,"RAIN_ICE_SLOW beg:PPRES")
  CALL MPPDB_CHECK(PLSFACT,"RAIN_ICE_SLOW beg:PLSFACT")
  CALL MPPDB_CHECK(PLVFACT,"RAIN_ICE_SLOW beg:PLVFACT")
  CALL MPPDB_CHECK(PSSI,"RAIN_ICE_SLOW beg:PSSI")
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PRVS,"RAIN_ICE_SLOW beg:PRVS")
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_SLOW beg:PRCS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_SLOW beg:PRIS")
  CALL MPPDB_CHECK(PRSS,"RAIN_ICE_SLOW beg:PRSS")
  CALL MPPDB_CHECK(PRGS,"RAIN_ICE_SLOW beg:PRGS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_SLOW beg:PTHS")
END IF
!
JLU = size(PRHODREF)
!
#ifndef MNH_OPENACC
ALLOCATE( GWORK   (JLU) )
ALLOCATE( ZZW     (JLU) )
ALLOCATE( ZCRIAUTI(JLU) )
ALLOCATE( zz_diff (size(PLSFACT))  )
#else
!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_SLOW' )

CALL MNH_MEM_GET( GWORK,    JLU )
CALL MNH_MEM_GET( ZZW,      JLU )
CALL MNH_MEM_GET( ZCRIAUTI, JLU )
CALL MNH_MEM_GET( zz_diff,  SIZE(PLSFACT) )

!$acc data present( gwork, zzw, zcriauti, zz_diff )
#endif
!
!*       3.2     compute the homogeneous nucleation source: RCHONI
!
!$acc kernels
  zz_diff(:) = plsfact(:) - plvfact(:)

  ZZW(:) = 0.0
  GWORK(:) = PZT(:)<XTT-35.0 .AND. PRCT(:)>XRTMIN(2) .AND. PRCS(:)>0.
  !$mnh_do_concurrent( JL=1:JLU )
  IF ( GWORK(JL) ) THEN
    ZZW(JL) = MIN( PRCS(JL),XHON*PRHODREF(JL)*PRCT(JL)       &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
                                 *EXP( MIN(XMNH_HUGE_12_LOG,XALPHA3*(PZT(JL)-XTT)-XBETA3) ) )
#else
                                 *BR_EXP(MIN(XMNH_HUGE_12_LOG, XALPHA3*(PZT(JL)-XTT)-XBETA3) ) )
#endif
    PRIS(JL) = PRIS(JL) + ZZW(JL)
    PRCS(JL) = PRCS(JL) - ZZW(JL)
    PTHS(JL) = PTHS(JL) + ZZW(JL) * zz_diff(JL) ! f(L_f*(RCHONI))
  END IF
  !$mnh_end_do()
!$acc end kernels

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'HON', &
                                           Unpack(  zzw(:) * zz_diff(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
  if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'HON', &
                                           Unpack( -zzw(:)              * prhodj(:), mask = omicro(:,:,:), field = 0.) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'HON', &
                                           Unpack(  zzw(:)              * prhodj(:), mask = omicro(:,:,:), field = 0.) )
!
!*       3.3     compute the spontaneous freezing source: RRHONG
!
!$acc kernels
!$mnh_do_concurrent (JL=1:JLU)  
  ZZW(JL) = 0.0
  GWORK(JL) = PZT(JL)<XTT-35.0 .AND. PRRT(JL)>XRTMIN(3) .AND. PRRS(JL)>0.
  IF( GWORK(JL) )THEN
    ZZW(JL) = MIN( PRRS(JL),PRRT(JL)* PINVTSTEP )
    PRGS(JL) = PRGS(JL) + ZZW(JL)
    PRRS(JL) = PRRS(JL) - ZZW(JL)
    PTHS(JL) = PTHS(JL) + ZZW(JL) * zz_diff(JL) ! f(L_f*(RRHONG))
  ENDIF
!$mnh_end_do()
!$acc end kernels

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'SFR', &
                                           Unpack(  zzw(:) * zz_diff(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'SFR', &
                                           Unpack( -zzw(:)              * prhodj(:), mask = omicro(:,:,:), field = 0.) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'SFR', &
                                           Unpack(  zzw(:)              * prhodj(:), mask = omicro(:,:,:), field = 0.) )
!
!*       3.4    compute the deposition, aggregation and autoconversion sources
!
!$acc kernels
  PKA(:) = 2.38E-2 + 0.0071E-2 * ( PZT(:) - XTT )          ! k_a
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
  PDV(:) = 0.211E-4 * (PZT(:)/XTT)**1.94 * (XP00/PPRES(:)) ! D_v
#else
!$mnh_do_concurrent ( JL=1:JLU )
   PDV(JL) = 0.211E-4 * BR_POW(PZT(JL)/XTT,1.94) * (XP00/PPRES(JL)) ! D_v
!$mnh_end_do()   
#endif
!
!*       3.4.1  compute the thermodynamical function A_i(T,P)
!*              and the c^prime_j (in the ventilation factor)
!
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
  PAI(:) = EXP( XALPI - XBETAI/PZT(:) - XGAMI*ALOG(PZT(:) ) ) ! es_i
  PAI(:) = ( XLSTT + (XCPV-XCI)*(PZT(:)-XTT) )**2 / (PKA(:)*XRV*PZT(:)**2) &
                                 + ( XRV*PZT(:) ) / (PDV(:)*PAI(:))
  PCJ(:) = XSCFAC * PRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(PZT(:)-XTT) )
#else
!$mnh_do_concurrent ( JL=1:JLU )
  PAI(JL) = BR_EXP( XALPI - XBETAI/PZT(JL) - XGAMI*BR_LOG(PZT(JL) ) ) ! es_i
  PAI(JL) = BR_P2( XLSTT + (XCPV-XCI)*(PZT(JL)-XTT) ) / (PKA(JL)*XRV*BR_P2(PZT(JL))) &
                                 + ( XRV*PZT(JL) ) / (PDV(JL)*PAI(JL))
  PCJ(JL) = XSCFAC * BR_POW(PRHODREF(JL),0.3) / BR_POW( 1.718E-5+0.0049E-5*(PZT(JL)-XTT) , 0.5)
!$mnh_end_do()
#endif
!
!*       3.4.2  compute the riming-conversion of r_c for r_i production: RCAUTI
!
!  ZZW(:) = 0.0
!  ZTIMAUTIC = SQRT( XTIMAUTI*XTIMAUTC )
!  WHERE ( (PRCT(:)>0.0) .AND. (PRIT(:)>0.0) .AND. (PRCS(:)>0.0) )
!    ZZW(:) = MIN( PRCS(:),ZTIMAUTIC * MAX( SQRT( PRIT(:)*PRCT(:) ),0.0 ) )
!    PRIS(:) = PRIS(:) + ZZW(:)
!    PRCS(:) = PRCS(:) - ZZW(:)
!    PTHS(:) = PTHS(:) + ZZW(:) * zz_diff(:) ! f(L_f*(RCAUTI))
!  END WHERE
!
!*       3.4.3  compute the deposition on r_s: RVDEPS
!
  GWORK(:) = PRST(:)>0.0
  !$mnh_do_concurrent( JL=1:JLU )
    IF ( GWORK(JL) ) THEN
      PLBDAS(JL)  = MIN( XLBDAS_MAX,                                           &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
                         XLBS*( PRHODREF(JL)*MAX( PRST(JL),XRTMIN(5) ) )**XLBEXS )
#else
                         XLBS*BR_POW( PRHODREF(JL)*MAX( PRST(JL),XRTMIN(5) ),XLBEXS ) )
#endif
    ELSE
      PLBDAS(JL) = 0.
    END IF
  !$mnh_end_do() ! CONCURRENT
  ZZW(:) = 0.0
  GWORK(:) = (PRST(:)>XRTMIN(5)) .AND. (PRSS(:)>0.0)
 !$mnh_do_concurrent ( JL=1:JLU )
    IF ( GWORK(JL) ) THEN
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
       ZZW(JL) = ( PRST(JL) * PLBDAS(JL)**XBS * PSSI(JL)/PAI(JL) ) *                                &
            ( X0DEPS*PLBDAS(JL)**XEX0DEPS + X1DEPS*PCJ(JL)*PLBDAS(JL)**XEX1DEPS )
#else
       ZZW(JL) = ( PRST(JL) * BR_POW(PLBDAS(JL),XBS) * PSSI(JL)/PAI(JL) ) *                                &
            ( X0DEPS*BR_POW(PLBDAS(JL),XEX0DEPS) + X1DEPS*PCJ(JL)*BR_POW(PLBDAS(JL),XEX1DEPS) )
#endif
       ZZW(JL) =         MIN( PRVS(JL),ZZW(JL)      )*(0.5+SIGN(0.5,ZZW(JL))) &
            - MIN( PRSS(JL),ABS(ZZW(JL)) )*(0.5-SIGN(0.5,ZZW(JL)))
       PRSS(JL) = PRSS(JL) + ZZW(JL)
       PRVS(JL) = PRVS(JL) - ZZW(JL)
       PTHS(JL) = PTHS(JL) + ZZW(JL)*PLSFACT(JL)
    END IF
 !$mnh_end_do()
!$acc end kernels

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'DEPS', &
                                           Unpack(  zzw(:) * plsfact(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
  if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'DEPS', &
                                           Unpack( -zzw(:)              * prhodj(:), mask = omicro(:,:,:), field = 0.) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'DEPS', &
                                           Unpack(  zzw(:)              * prhodj(:), mask = omicro(:,:,:), field = 0.) )
!
!*       3.4.4  compute the aggregation on r_s: RIAGGS
!
!$acc kernels
  ZZW(:) = 0.0
  GWORK(:) = PRIT(:)>XRTMIN(4) .AND. PRST(:)>XRTMIN(5) .AND. PRIS(:)>0.0
  !$mnh_do_concurrent ( JL=1:JLU )
    IF ( GWORK(JL) ) THEN
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)  
      ZZW(JL) = MIN( PRIS(JL),XFIAGGS * EXP( XCOLEXIS*(PZT(JL)-XTT) ) &
                                      * PRIT(JL)                      &
                                      * PRST(JL) * PLBDAS(JL)**(XBS+XEXIAGGS)          &
                                      * PRHODREF(JL)**(-XCEXVT+1)       )
#else
      ZZW(JL) = MIN( PRIS(JL),XFIAGGS * BR_EXP( XCOLEXIS*(PZT(JL)-XTT) ) &
                                      * PRIT(JL)                         &
                                      * PRST(JL) * BR_POW(PLBDAS(JL),XBS+XEXIAGGS)      &
                                      * BR_POW(PRHODREF(JL),-XCEXVT+1)     )
#endif
      PRSS(JL)  = PRSS(JL)  + ZZW(JL)
      PRIS(JL)  = PRIS(JL)  - ZZW(JL)
   END IF
  !$mnh_end_do() ! CONCURRENT
!$acc end kernels

  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'AGGS', &
                                           Unpack( -zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'AGGS', &
                                           Unpack(  zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
!
!*       3.4.5  compute the autoconversion of r_i for r_s production: RIAUTS
!
!$acc kernels
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
  ZCRIAUTI(:)=MIN(XCRIAUTI,10**(XACRIAUTI*(PZT(:)-XTT)+XBCRIAUTI))
#else
!$mnh_do_concurrent ( JL=1:JLU )
  ZCRIAUTI(JL)=MIN(XCRIAUTI, BR_POW(10.,XACRIAUTI*(PZT(JL)-XTT)+XBCRIAUTI) )
!$mnh_end_do()
#endif
  ZZW(:) = 0.0
  GWORK(:) = PRIT(:)>XRTMIN(4) .AND. PRIS(:)>0.0
!$mnh_do_concurrent( JL=1:JLU )
  IF ( GWORK(JL) ) THEN
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
    ZZW(JL) = MIN( PRIS(JL),XTIMAUTI * EXP( XTEXAUTI*(PZT(JL)-XTT) ) &
                            * MAX( PRIT(JL)-ZCRIAUTI(JL),0.0 ) )
#else
    ZZW(JL) = MIN( PRIS(JL),XTIMAUTI * BR_EXP( XTEXAUTI*(PZT(JL)-XTT) ) &
                            * MAX( PRIT(JL)-ZCRIAUTI(JL),0.0 ) )
#endif
    PRSS(JL)  = PRSS(JL)  + ZZW(JL)
    PRIS(JL)  = PRIS(JL)  - ZZW(JL)
  !!END WHERE
  END IF
!$mnh_end_do()
!$acc end kernels

  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'AUTS', &
                                           Unpack( -zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'AUTS', &
                                           Unpack(  zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
!
!*       3.4.6  compute the deposition on r_g: RVDEPG
!
!
!$acc kernels
  GWORK(:) = PRGT(:)>0.0
  !$mnh_do_concurrent ( JL=1:JLU )
    IF ( GWORK(JL) ) THEN
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)  
      PLBDAG(JL)  = XLBG*( PRHODREF(JL)*MAX( PRGT(JL),XRTMIN(6) ) )**XLBEXG
#else
      PLBDAG(JL)  = XLBG*BR_POW( PRHODREF(JL)*MAX( PRGT(JL),XRTMIN(6) ), XLBEXG)
#endif
    ELSE
      PLBDAG(JL) = 0.
    END IF
  !$mnh_end_do() ! CONCURRENT
  ZZW(:) = 0.0
  GWORK(:) = PRGT(:)>XRTMIN(6) .AND. PRGS(:)>0.0
  !$mnh_do_concurrent ( JL=1:JLU )
    IF ( GWORK(JL) ) THEN
      ZZW(JL) = ( PSSI(JL)/(PRHODREF(JL)*PAI(JL)) ) *                               &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
                ( X0DEPG*PLBDAG(JL)**XEX0DEPG + X1DEPG*PCJ(JL)*PLBDAG(JL)**XEX1DEPG )
#else
                ( X0DEPG*BR_POW(PLBDAG(JL),XEX0DEPG) + X1DEPG*PCJ(JL)*BR_POW(PLBDAG(JL),XEX1DEPG) )
#endif
      ZZW(JL) =   MIN( PRVS(JL),ZZW(JL)      )*(0.5+SIGN(0.5,ZZW(JL))) &
                - MIN( PRGS(JL),ABS(ZZW(JL)) )*(0.5-SIGN(0.5,ZZW(JL)))
      PRGS(JL) = PRGS(JL) + ZZW(JL)
      PRVS(JL) = PRVS(JL) - ZZW(JL)
      PTHS(JL) = PTHS(JL) + ZZW(JL)*PLSFACT(JL)
    END IF
  !$mnh_end_do() ! CONCURRENT
!$acc end kernels

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'DEPG', &
                                           Unpack(  zzw(:) * plsfact(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
  if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'DEPG', &
                                           Unpack( -zzw(:)              * prhodj(:), mask = omicro(:,:,:), field = 0.) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'DEPG', &
                                           Unpack(  zzw(:)              * prhodj(:), mask = omicro(:,:,:), field = 0.) )
IF (MPPDB_INITIALIZED) THEN
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PRVS,"RAIN_ICE_SLOW end:PRVS")
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_SLOW end:PRCS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_SLOW end:PRIS")
  CALL MPPDB_CHECK(PRSS,"RAIN_ICE_SLOW end:PRSS")
  CALL MPPDB_CHECK(PRGS,"RAIN_ICE_SLOW end:PRGS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_SLOW end:PTHS")
  !Check all OUT arrays
  CALL MPPDB_CHECK(PAI,"RAIN_ICE_SLOW end:PAI")
  CALL MPPDB_CHECK(PCJ,"RAIN_ICE_SLOW end:PCJ")
  CALL MPPDB_CHECK(PKA,"RAIN_ICE_SLOW end:PKA")
  CALL MPPDB_CHECK(PDV,"RAIN_ICE_SLOW end:PDV")
  CALL MPPDB_CHECK(PLBDAS,"RAIN_ICE_SLOW end:PLBDAS")
  CALL MPPDB_CHECK(PLBDAG,"RAIN_ICE_SLOW end:PLBDAG")
END IF

!$acc end data

#ifdef MNH_OPENACC
!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'RAIN_ICE_SLOW' )
#endif

!$acc end data

END SUBROUTINE RAIN_ICE_SLOW

END MODULE MODE_RAIN_ICE_SLOW
