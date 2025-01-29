!MNH_LIC Copyright 1995-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 25/02/2019: split rain_ice (cleaner and easier to maintain/debug)
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 03/06/2019: remove PACK/UNPACK intrinsics (to get more performance and better OpenACC support)
!  P. Wautelet 05/06/2019: optimisations
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!  J. Escobar  11/08/2020: Bypass PGI/NVHPC OPENACC BUG, error 700: Illegal address during kernel execution => DO CONCURRENT
!  J. Escobar  12/08/2020: Bypass PGI/NVHPC OPENACC BUG data partially present => enter data in ini_rain_ce & DO CONCURRENT
!  J. Escobar  13/08/2020: Openacc PB , missing enter/update data rain_ice_fast_rs/g
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_FAST_RG

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: RAIN_ICE_FAST_RG

CONTAINS

SUBROUTINE RAIN_ICE_FAST_RG(KRR, OMICRO, PRHODREF, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PCIT, &
                            PRHODJ, PPRES, PZT, PLBDAR, PLBDAS, PLBDAG, PLSFACT, PLVFACT, &
                            PCJ, PKA, PDV, &
                            PRCS, PRRS, PRIS, PRSS, PRGS, PRHS, PTHS, &
                            PUSW, PRDRYG, PRWETG)

!
!*      0. DECLARATIONS
!          ------------
!
use modd_budget,         only: lbudget_th, lbudget_rc, lbudget_rr, lbudget_ri, lbudget_rs, lbudget_rg, lbudget_rh, &
                               NBUDGET_TH, NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, &
                               tbudgets
USE MODD_CST,            only: XCI, XCL, XCPV, XESTT, XLMTT, XLVTT, XMD, XMV, XRV, XTT
USE MODD_RAIN_ICE_DESCR_n, only: XBS, XCEXVT, XCXG, XCXS, XDG, XRTMIN
USE MODD_RAIN_ICE_PARAM_n, only: NDRYLBDAG, NDRYLBDAR, NDRYLBDAS, X0DEPG, X1DEPG, XCOLEXIG, XCOLEXSG, XCOLIG, XCOLSG, XDRYINTP1G, &
                               XDRYINTP1R, XDRYINTP1S, XDRYINTP2G, XDRYINTP2R, XDRYINTP2S, XEX0DEPG, XEX1DEPG, XEXICFRR,        &
                               XEXRCFRI, XFCDRYG, XFIDRYG, XFRDRYG, XFSDRYG, XICFRR, XKER_RDRYG, XKER_SDRYG, XLBRDRYG1,         &
                               XLBRDRYG2, XLBRDRYG3, XLBSDRYG1, XLBSDRYG2, XLBSDRYG3, XRCFRI

use mode_budget,         only: Budget_store_add, Budget_store_end, Budget_store_init
#ifdef MNH_OPENACC
USE MODE_MNH_ZWORK,      ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
#endif
use mode_mppdb
#ifndef MNH_OPENACC
use mode_tools,                        only: Countjv
#else
use mode_tools,                        only: Countjv_device
#endif

#if defined(MNH_BITREP) || defined(MNH_BITREP_OMP)
USE MODI_BITREP
#endif
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
! mnh_undef(OPENACC)
#endif

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                    INTENT(IN)    :: KRR      ! Number of moist variables
LOGICAL,  DIMENSION(:,:,:), intent(in)    :: OMICRO   ! Test where to compute all processes
REAL,     DIMENSION(:),     intent(in)    :: PRHODREF ! RHO Dry REFerence
REAL,     DIMENSION(:),     intent(in)    :: PRVT     ! Water vapor m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRCT     ! Cloud water m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRRT     ! Rain water m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRIT     ! Pristine ice m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRST     ! Snow/aggregate m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRGT     ! Graupel m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PCIT     ! Pristine ice conc. at t
REAL,     DIMENSION(:),     intent(in)    :: PRHODJ   ! RHO times Jacobian
REAL,     DIMENSION(:),     intent(in)    :: PPRES    ! Pressure
REAL,     DIMENSION(:),     intent(in)    :: PZT      ! Temperature
REAL,     DIMENSION(:),     intent(in)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL,     DIMENSION(:),     intent(in)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL,     DIMENSION(:),     intent(in)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL,     DIMENSION(:),     intent(in)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL,     DIMENSION(:),     intent(in)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL,     DIMENSION(:),     intent(in)    :: PCJ      ! Function to compute the ventilation coefficient
REAL,     DIMENSION(:),     intent(in)    :: PKA      ! Thermal conductivity of the air
REAL,     DIMENSION(:),     intent(in)    :: PDV      ! Diffusivity of water vapor in the air
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRCS     ! Cloud water m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRRS     ! Rain water m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRIS     ! Pristine ice m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRSS     ! Snow/aggregate m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRGS     ! Graupel m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRHS     ! Hail m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PTHS     ! Theta source
REAL,     DIMENSION(:),     intent(inout) :: PUSW     ! Undersaturation over water
REAL,     DIMENSION(:),     intent(out)   :: PRDRYG   ! Dry growth rate of the graupeln
REAL,     DIMENSION(:),     intent(out)   :: PRWETG   ! Wet growth rate of the graupeln
!
!*       0.2  declaration of local variables
!
INTEGER                              :: IGDRY
INTEGER                              :: JJ, JL
#ifndef MNH_OPENACC
INTEGER, DIMENSION(:),   ALLOCATABLE :: I1
INTEGER, DIMENSION(:),   ALLOCATABLE :: IVEC1, IVEC2      ! Vectors of indices for interpolations
LOGICAL, DIMENSION(:),   ALLOCATABLE :: GWORK
REAL,    DIMENSION(:),   ALLOCATABLE :: ZZW               ! Work array
REAL,    DIMENSION(:),   ALLOCATABLE :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for interpolations
REAL,    DIMENSION(:),   ALLOCATABLE :: ZVECLBDAG, ZVECLBDAR, ZVECLBDAS
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZZW1              ! Work arrays
#else
INTEGER, DIMENSION(:),   pointer, contiguous :: I1
INTEGER, DIMENSION(:),   pointer, contiguous :: IVEC1, IVEC2      ! Vectors of indices for interpolations
LOGICAL, DIMENSION(:),   pointer, contiguous :: GWORK
REAL,    DIMENSION(:),   pointer, contiguous :: ZZW               ! Work array
REAL,    DIMENSION(:),   pointer, contiguous :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for interpolations
REAL,    DIMENSION(:),   pointer, contiguous :: ZVECLBDAG, ZVECLBDAR, ZVECLBDAS
REAL,    DIMENSION(:,:), pointer, contiguous :: ZZW1              ! Work arrays
#endif
!
INTEGER                              :: JLU
!-------------------------------------------------------------------------------
!
! IN variables
!
!$acc data present( OMICRO, PRHODREF, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PCIT,   &
!$acc &             PRHODJ, PPRES, PZT, PLBDAR, PLBDAS, PLBDAG, PLSFACT, PLVFACT, &
!$acc &             PCJ, PKA, PDV,                                                &
!
! INOUT variables
!
!$acc &             PRCS, PRRS, PRIS, PRSS, PRGS, PRHS, PTHS, PUSW,               &
!
! OUT variables
!
!$acc &             PRDRYG, PRWETG, &
!
! use variables
!
!$acc &              XKER_SDRYG,XKER_RDRYG )

IF (MPPDB_INITIALIZED) THEN
  !Check all IN arrays
  CALL MPPDB_CHECK(OMICRO,"RAIN_ICE_FAST_RG beg:OMICRO")
  CALL MPPDB_CHECK(PRHODREF,"RAIN_ICE_FAST_RG beg:PRHODREF")
  CALL MPPDB_CHECK(PRVT,"RAIN_ICE_FAST_RG beg:PRVT")
  CALL MPPDB_CHECK(PRCT,"RAIN_ICE_FAST_RG beg:PRCT")
  CALL MPPDB_CHECK(PRRT,"RAIN_ICE_FAST_RG beg:PRRT")
  CALL MPPDB_CHECK(PRIT,"RAIN_ICE_FAST_RG beg:PRIT")
  CALL MPPDB_CHECK(PRST,"RAIN_ICE_FAST_RG beg:PRST")
  CALL MPPDB_CHECK(PRGT,"RAIN_ICE_FAST_RG beg:PRGT")
  CALL MPPDB_CHECK(PCIT,"RAIN_ICE_FAST_RG beg:PCIT")
  CALL MPPDB_CHECK(PRHODJ,"RAIN_ICE_FAST_RG beg:PRHODJ")
  CALL MPPDB_CHECK(PPRES,"RAIN_ICE_FAST_RG beg:PPRES")
  CALL MPPDB_CHECK(PZT,"RAIN_ICE_FAST_RG beg:PZT")
  CALL MPPDB_CHECK(PLBDAR,"RAIN_ICE_FAST_RG beg:PLBDAR")
  CALL MPPDB_CHECK(PLBDAS,"RAIN_ICE_FAST_RG beg:PLBDAS")
  CALL MPPDB_CHECK(PLBDAG,"RAIN_ICE_FAST_RG beg:PLBDAG")
  CALL MPPDB_CHECK(PLSFACT,"RAIN_ICE_FAST_RG beg:PLSFACT")
  CALL MPPDB_CHECK(PLVFACT,"RAIN_ICE_FAST_RG beg:PLVFACT")
  CALL MPPDB_CHECK(PCJ,"RAIN_ICE_FAST_RG beg:PCJ")
  CALL MPPDB_CHECK(PKA,"RAIN_ICE_FAST_RG beg:PKA")
  CALL MPPDB_CHECK(PDV,"RAIN_ICE_FAST_RG beg:PDV")
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_FAST_RG beg:PRCS")
  CALL MPPDB_CHECK(PRRS,"RAIN_ICE_FAST_RG beg:PRRS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_FAST_RG beg:PRIS")
  CALL MPPDB_CHECK(PRSS,"RAIN_ICE_FAST_RG beg:PRSS")
  CALL MPPDB_CHECK(PRGS,"RAIN_ICE_FAST_RG beg:PRGS")
  CALL MPPDB_CHECK(PRHS,"RAIN_ICE_FAST_RG beg:PRHS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_FAST_RG beg:PTHS")
  CALL MPPDB_CHECK(PUSW,"RAIN_ICE_FAST_RG beg:PUSW")
  !Check use variable
  CALL MPPDB_CHECK(XKER_SDRYG,"RAIN_ICE_FAST_RG beg:XKER_SDRYG")
  CALL MPPDB_CHECK(XKER_SDRYG,"RAIN_ICE_FAST_RG beg:XKER_RDRYG")
END IF
!
JLU = size(PRHODREF) 
!
#ifndef MNH_OPENACC
ALLOCATE( I1   (size(PRHODREF)) )
ALLOCATE( GWORK(size(PRHODREF)) )
ALLOCATE( ZZW  (size(PRHODREF)) )
ALLOCATE( ZZW1 (size(PRHODREF),7) )
#else
!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_FAST_RG 1' )

CALL MNH_MEM_GET( I1,    SIZE(PRHODREF) )
CALL MNH_MEM_GET( GWORK, SIZE(PRHODREF) )
CALL MNH_MEM_GET( ZZW,   SIZE(PRHODREF) )
CALL MNH_MEM_GET( ZZW1,  SIZE(PRHODREF), 7 )

!$acc data present( I1, GWORK, ZZW, ZZW1 )
#endif

!
!*       6.1    rain contact freezing
!
!$acc kernels present_cr(ZZW1,GWORK)
  ZZW1(:,:) = 0.0
  GWORK(:) = PRIT(:)>XRTMIN(4) .AND. PRRT(:)>XRTMIN(3) .AND. PRIS(:)>0.0 .AND. PRRS(:)>0.0
!$acc end kernels
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
!$acc parallel present_cr(ZZW1,GWORK)
  !$mnh_do_concurrent ( JL=1:JLU )
    IF ( GWORK(JL) ) THEN
      ZZW1(JL,3) = MIN( PRIS(JL),XICFRR * PRIT(JL)                & ! RICFRRG
                                      * PLBDAR(JL)**XEXICFRR      &
                                      * PRHODREF(JL)**(-XCEXVT) )
      ZZW1(JL,4) = MIN( PRRS(JL),XRCFRI * PCIT(JL)                & ! RRCFRIG
                                      * PLBDAR(JL)**XEXRCFRI      &
                                      * PRHODREF(JL)**(-XCEXVT-1.) )
      PRIS(JL) = PRIS(JL) - ZZW1(JL,3)
      PRRS(JL) = PRRS(JL) - ZZW1(JL,4)
      PRGS(JL) = PRGS(JL) + ZZW1(JL,3)+ZZW1(JL,4)
      PTHS(JL) = PTHS(JL) + ZZW1(JL,4)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*RRCFRIG)
    END IF
  !$mnh_end_do() ! CONCURRENT
!$acc end parallel
#else
!$acc parallel present_cr(ZZW1,GWORK)
  !$mnh_do_concurrent ( JL=1:JLU )
    IF ( GWORK(JL) ) THEN
       ZZW1(JL,3) = MIN( PRIS(JL),XICFRR * PRIT(JL)                     & ! RICFRRG
            * BR_POW(PLBDAR(JL),XEXICFRR)  &
            * BR_POW(PRHODREF(JL),-XCEXVT) )
       ZZW1(JL,4) = MIN( PRRS(JL),XRCFRI * PCIT(JL)                     & ! RRCFRIG
            * BR_POW(PLBDAR(JL),XEXRCFRI)  &
            * BR_POW(PRHODREF(JL),-XCEXVT-1.) )
       PRIS(JL) = PRIS(JL) - ZZW1(JL,3)
       PRRS(JL) = PRRS(JL) - ZZW1(JL,4)
       PRGS(JL) = PRGS(JL) + ZZW1(JL,3)+ZZW1(JL,4)
       PTHS(JL) = PTHS(JL) + ZZW1(JL,4)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*RRCFRIG)
    END IF
  !$mnh_end_do() ! CONCURRENT
!$acc end parallel
#endif 


IF (MPPDB_INITIALIZED) THEN
   CALL MPPDB_CHECK(PRRS,"RAIN_ICE_FAST_RG 6.1:PRRS")
END IF

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'CFRZ', Unpack (  zzw1(:, 4) * ( plsfact(:) - plvfact(:) ) &
                                                                     * prhodj(:), mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'CFRZ', Unpack ( -zzw1(:, 4)                * prhodj(:), &
                                                                                  mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'CFRZ', Unpack ( -zzw1(:, 3)                * prhodj(:), &
                                                                                  mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'CFRZ', Unpack ( (  zzw1(:,3) + zzw1(:,4) ) * prhodj(:), &
                                                                                  mask = omicro(:,:,:), field = 0. ) )

!PW:used init/end instead of add because zzw1 is produced with a where(...) and is used with other where(...)
! => can not use directly zzw1 in Budget_store_add
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'WETG', Unpack ( pths(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'WETG', Unpack ( prcs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'WETG', Unpack ( prrs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'WETG', Unpack ( pris(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'WETG', Unpack ( prss(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'WETG', Unpack ( prgs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'WETG', Unpack ( prhs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
!
!*       6.2    compute the Dry growth case
!
!$acc kernels present_cr(GWORK)
  ZZW1(:,:) = 0.0
!$acc end kernels
!$acc parallel present_cr(GWORK)
  !$mnh_expand_where(JL=1:JLU)
  GWORK(:) = PRGT(:)>XRTMIN(6) .AND. PRCT(:)>XRTMIN(2) .AND. PRCS(:)>0.0
  WHERE( GWORK(:) )
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
    ZZW(:) = PLBDAG(:)**(XCXG-XDG-2.0) * PRHODREF(:)**(-XCEXVT)
#else
    ZZW(:) = BR_POW(PLBDAG(:),XCXG-XDG-2.0) * BR_POW(PRHODREF(:),-XCEXVT)
#endif
    ZZW1(:,1) = MIN( PRCS(:),XFCDRYG * PRCT(:) * ZZW(:) )             ! RCDRYG
  END WHERE
  !$mnh_end_expand_where()
!$acc end parallel
!$acc parallel present_cr(GWORK)
  !$mnh_expand_where(JL=1:JLU)
  GWORK(:) = PRGT(:)>XRTMIN(6) .AND. PRIT(:)>XRTMIN(4) .AND. PRIS(:)>0.0
  WHERE( GWORK(:) )
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
    ZZW(:) = PLBDAG(:)**(XCXG-XDG-2.0) * PRHODREF(:)**(-XCEXVT)
    ZZW1(:,2) = MIN( PRIS(:),XFIDRYG * EXP( XCOLEXIG*(PZT(:)-XTT) ) &
                                     * PRIT(:) * ZZW(:) )             ! RIDRYG
#else
    ZZW(:) = BR_POW(PLBDAG(:),XCXG-XDG-2.0) * BR_POW(PRHODREF(:),-XCEXVT)
    ZZW1(:,2) = MIN( PRIS(:),XFIDRYG * BR_EXP( XCOLEXIG*(PZT(:)-XTT) ) &
                                     * PRIT(:) * ZZW(:) )             ! RIDRYG
#endif
  END WHERE
  !$mnh_end_expand_where()
!$acc end parallel
IF (MPPDB_INITIALIZED) THEN
   CALL MPPDB_CHECK(ZZW1,"RAIN_ICE_FAST_RG 6.2:ZZW1")
END IF   
!
!*       6.2.1  accretion of aggregates on the graupeln
!
!$acc kernels present_cr(GWORK)
  GWORK(:) = PRST(:)>XRTMIN(5) .AND. PRGT(:)>XRTMIN(6) .AND. PRSS(:)>0.0
!$acc end kernels
#ifndef MNH_OPENACC
  IGDRY = COUNTJV( GWORK(:), I1(:) )
#else
  CALL COUNTJV_DEVICE( GWORK(:), I1(:), IGDRY )
#endif

  IF( IGDRY>0 ) THEN
!
!*       6.2.2  allocations
!
#ifndef MNH_OPENACC
    ALLOCATE(ZVECLBDAG(IGDRY))
    ALLOCATE(ZVECLBDAS(IGDRY))
    ALLOCATE(ZVEC1(IGDRY))
    ALLOCATE(ZVEC2(IGDRY))
    ALLOCATE(ZVEC3(IGDRY))
    ALLOCATE(IVEC1(IGDRY))
    ALLOCATE(IVEC2(IGDRY))
#else
    !Pin positions in the pools of MNH memory
    CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_FAST_RG 2' )

    CALL MNH_MEM_GET( ZVECLBDAG, IGDRY )
    CALL MNH_MEM_GET( ZVECLBDAS, IGDRY )
    CALL MNH_MEM_GET( ZVEC1,     IGDRY )
    CALL MNH_MEM_GET( ZVEC2,     IGDRY )
    CALL MNH_MEM_GET( ZVEC3,     IGDRY )
    CALL MNH_MEM_GET( IVEC1,     IGDRY )
    CALL MNH_MEM_GET( IVEC2,     IGDRY )

!$acc data present( ZVECLBDAG, ZVECLBDAS, ZVEC1, ZVEC2, ZVEC3, IVEC1, IVEC2 )
#endif
!
!*       6.2.3  select the (PLBDAG,PLBDAS) couplet
!
!$acc parallel present_cr(ZVECLBDAG,ZVECLBDAS)
 !$mnh_expand_where(JL=1:IGDRY)
    ZVECLBDAG(1:IGDRY) = PLBDAG(I1(1:IGDRY))
    ZVECLBDAS(1:IGDRY) = PLBDAS(I1(1:IGDRY))
!
!*       6.2.4  find the next lower indice for the PLBDAG and for the PLBDAS
!               in the geometrical set of (Lbda_g,Lbda_s) couplet use to
!               tabulate the SDRYG-kernel
!
    ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAG)-0.00001,           &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
                          XDRYINTP1G * LOG( ZVECLBDAG(1:IGDRY) ) + XDRYINTP2G ) )
#else
                          XDRYINTP1G * BR_LOG( ZVECLBDAG(1:IGDRY) ) + XDRYINTP2G ) )
#endif
    IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
    ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - REAL( IVEC1(1:IGDRY) )
!
    ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAS)-0.00001,           &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
                          XDRYINTP1S * LOG( ZVECLBDAS(1:IGDRY) ) + XDRYINTP2S ) )
#else
                          XDRYINTP1S * BR_LOG( ZVECLBDAS(1:IGDRY) ) + XDRYINTP2S ) )
#endif
    IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
    ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
  !$mnh_end_expand_where()
!$acc end parallel
!$acc kernels ! present_cr(ZVECLBDAG,ZVECLBDAS)
    !
!*       6.2.5  perform the bilinear interpolation of the normalized
!               SDRYG-kernel
!
    !$mnh_do_concurrent ( JJ = 1:IGDRY )
      ZVEC3(JJ) =  (  XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * ZVEC1(JJ) &
                 - (  XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
    !$mnh_end_do() ! CONCURRENT
!
    ! acc loop independent , private (JL) 
    !$mnh_do_concurrent (JJ=1:IGDRY)
      JL = I1(JJ)
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
      ZZW1(JL,3) = MIN( PRSS(JL),XFSDRYG*ZVEC3(JJ)                           & ! RSDRYG
                                      * EXP( XCOLEXSG*(PZT(JL)-XTT) )        &
                    * PRST(JJ) * ZVECLBDAG(JJ)**XCXG        &
                    * PRHODREF(JL)**(-XCEXVT)                             &
                         * ( XLBSDRYG1 /   ZVECLBDAG(JJ)**2                  &
                           + XLBSDRYG2 / ( ZVECLBDAG(JJ)   * ZVECLBDAS(JJ) ) &
                           + XLBSDRYG3 /   ZVECLBDAS(JJ)**2 ) )
#else
      ZZW1(JL,3) = MIN( PRSS(JL),XFSDRYG*ZVEC3(JJ)                                & ! RSDRYG
                                      * BR_EXP( XCOLEXSG*(PZT(JL)-XTT) )          &
                    * PRST(JL) * BR_POW(ZVECLBDAG(JJ),XCXG) &
                    * BR_POW(PRHODREF(JL),-XCEXVT)                             &
                         * ( XLBSDRYG1 /   BR_P2(ZVECLBDAG(JJ))                   &
                           + XLBSDRYG2 / ( ZVECLBDAG(JJ)   * ZVECLBDAS(JJ) )      &
                           + XLBSDRYG3 /   BR_P2(ZVECLBDAS(JJ)) ) )
#endif
    !$mnh_end_do()
!$acc end kernels
IF (MPPDB_INITIALIZED) THEN
    CALL MPPDB_CHECK(ZZW1,"RAIN_ICE_FAST_RG 6.2.5:ZZW1")
    CALL MPPDB_CHECK(ZVEC3,"RAIN_ICE_FAST_RG 6.2.5:ZVEC3")
END IF    
!$acc end data
#ifndef MNH_OPENACC
    DEALLOCATE(ZVECLBDAS)
    DEALLOCATE(ZVECLBDAG)
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
#else
    !Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
    CALL MNH_MEM_RELEASE( 'RAIN_ICE_FAST_RG 2' )
#endif
  END IF
!
!*       6.2.6  accretion of raindrops on the graupeln
!
!$acc kernels present_cr(GWORK)
  GWORK(:) = PRRT(:)>XRTMIN(3) .AND. PRGT(:)>XRTMIN(6) .AND. PRSS(:)>0.0
!$acc end kernels
#ifndef MNH_OPENACC
  IGDRY = COUNTJV( GWORK(:), I1(:) )
#else
  CALL COUNTJV_DEVICE( GWORK(:), I1(:), IGDRY )
#endif
!
  IF( IGDRY>0 ) THEN
!
!*       6.2.7  allocations
!
#ifndef MNH_OPENACC
    ALLOCATE(ZVECLBDAG(IGDRY))
    ALLOCATE(ZVECLBDAR(IGDRY))
    ALLOCATE(ZVEC1(IGDRY))
    ALLOCATE(ZVEC2(IGDRY))
    ALLOCATE(ZVEC3(IGDRY))
    ALLOCATE(IVEC1(IGDRY))
    ALLOCATE(IVEC2(IGDRY))
#else
    !Pin positions in the pools of MNH memory
    CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_FAST_RG 3' )

    CALL MNH_MEM_GET( ZVECLBDAG, IGDRY )
    CALL MNH_MEM_GET( ZVECLBDAR, IGDRY )
    CALL MNH_MEM_GET( ZVEC1,     IGDRY )
    CALL MNH_MEM_GET( ZVEC2,     IGDRY )
    CALL MNH_MEM_GET( ZVEC3,     IGDRY )
    CALL MNH_MEM_GET( IVEC1,     IGDRY )
    CALL MNH_MEM_GET( IVEC2,     IGDRY )

!$acc data present( ZVECLBDAG, ZVECLBDAR, ZVEC1, ZVEC2, ZVEC3, IVEC1, IVEC2 )
#endif
!
!*       6.2.8  select the (PLBDAG,PLBDAR) couplet
!
!$acc parallel present_cr(ZVECLBDAG,ZVECLBDAR)
 !$mnh_expand_where(JL=1:IGDRY)    
    ZVECLBDAG(1:IGDRY) = PLBDAG(I1(1:IGDRY))
    ZVECLBDAR(1:IGDRY) = PLBDAR(I1(1:IGDRY))
!
!*       6.2.9  find the next lower indice for the PLBDAG and for the PLBDAR
!               in the geometrical set of (Lbda_g,Lbda_r) couplet use to
!               tabulate the RDRYG-kernel
!
    ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAG)-0.00001,           &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
                          XDRYINTP1G * LOG( ZVECLBDAG(1:IGDRY) ) + XDRYINTP2G ) )
#else
                          XDRYINTP1G * BR_LOG( ZVECLBDAG(1:IGDRY) ) + XDRYINTP2G ) )
#endif
    IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
    ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - REAL( IVEC1(1:IGDRY) )
!
    ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAR)-0.00001,           &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
                          XDRYINTP1R * LOG( ZVECLBDAR(1:IGDRY) ) + XDRYINTP2R ) )
#else
                          XDRYINTP1R * BR_LOG( ZVECLBDAR(1:IGDRY) ) + XDRYINTP2R ) )
#endif
    IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
    ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
  !$mnh_end_expand_where()
!$acc end parallel 
!
!*       6.2.10 perform the bilinear interpolation of the normalized
!               RDRYG-kernel
!
!$acc kernels ! present_cr(ZVECLBDAG,ZVECLBDAR,ZZW1)    
    !$mnh_do_concurrent (JJ = 1:IGDRY )
      ZVEC3(JJ) =  (  XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                  * ZVEC1(JJ) &
                 - (  XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
    !$mnh_end_do() ! CONCURRENT
    ! acc loop independent , private (JL)
    !$mnh_do_concurrent (JJ=1:IGDRY)
      JL = I1(JJ)
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
      ZZW1(JL,4) = MIN( PRRS(JL),XFRDRYG*ZVEC3(JJ)                  & ! RRDRYG
                        * ZVECLBDAR(JJ)**(-4) * ZVECLBDAG(JJ)**XCXG &
                               * PRHODREF(JL)**(-XCEXVT-1.)         &
                    * ( XLBRDRYG1/  ZVECLBDAG(JJ)**2                &
                      + XLBRDRYG2/( ZVECLBDAG(JJ) * ZVECLBDAR(JJ) ) &
                      + XLBRDRYG3/  ZVECLBDAR(JJ)**2 ) )
#else
      ZZW1(JL,4) = MIN( PRRS(JL),XFRDRYG*ZVEC3(JJ)                               & ! RRDRYG
                        * BR_POW(ZVECLBDAR(JJ),-4.) * BR_POW(ZVECLBDAG(JJ),XCXG) &
                               * BR_POW(PRHODREF(JL),-XCEXVT-1.)                 &
                    * ( XLBRDRYG1/  BR_P2(ZVECLBDAG(JJ))                         &
                      + XLBRDRYG2/( ZVECLBDAG(JJ) * ZVECLBDAR(JJ) )              &
                      + XLBRDRYG3/  BR_P2(ZVECLBDAR(JJ)) ) )
#endif
    !$mnh_end_do()
!$acc end kernels
IF (MPPDB_INITIALIZED) THEN
    CALL MPPDB_CHECK(ZZW1,"RAIN_ICE_FAST_RG 6.2.10:ZZW1")
    CALL MPPDB_CHECK(ZVEC3,"RAIN_ICE_FAST_RG 6.2.10:ZVEC3")
END IF    
!$acc end data
#ifndef MNH_OPENACC
    DEALLOCATE(ZVECLBDAR)
    DEALLOCATE(ZVECLBDAG)
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
#else
    !Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
    CALL MNH_MEM_RELEASE( 'RAIN_ICE_FAST_RG 3' )
#endif
  END IF
!
!$acc kernels present_cr(GWORK)
  PRDRYG(:) = ZZW1(:,1) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,4)
!
!*       6.3    compute the Wet growth case
!
  PRWETG(:) = 0.0
  GWORK(:) = PRGT(:)>XRTMIN(6)
!$acc end kernels  
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
!$acc parallel present_cr(GWORK)
  !$mnh_do_concurrent ( JL=1:JLU )
    IF ( GWORK(JL) ) THEN
       ZZW1(JL,5) = MIN( PRIS(JL),                                    &
            ZZW1(JL,2) / (XCOLIG*EXP(XCOLEXIG*(PZT(JL)-XTT)) ) ) ! RIWETG
       ZZW1(JL,6) = MIN( PRSS(JL),                                    &
            ZZW1(JL,3) / (XCOLSG*EXP(XCOLEXSG*(PZT(JL)-XTT)) ) ) ! RSWETG
       !
       ZZW(JL) = PRVT(JL)*PPRES(JL)/((XMV/XMD)+PRVT(JL)) ! Vapor pressure
       ZZW(JL) =   PKA(JL)*(XTT-PZT(JL)) +                              &
            ( PDV(JL)*(XLVTT + ( XCPV - XCL ) * ( PZT(JL) - XTT )) &
            *(XESTT-ZZW(JL))/(XRV*PZT(JL))           )
       !
       ! compute RWETG
       !
       PRWETG(JL)=MAX( 0.0,                                        &
            ( ZZW(JL) * ( X0DEPG*        PLBDAG(JL)**XEX0DEPG +    &
                          X1DEPG*PCJ(JL)*PLBDAG(JL)**XEX1DEPG ) +  &
            ( ZZW1(JL,5)+ZZW1(JL,6) ) *                            &
            ( PRHODREF(JL)*(XLMTT+(XCI-XCL)*(XTT-PZT(JL)))   ) ) / &
            ( PRHODREF(JL)*(XLMTT-XCL*(XTT-PZT(JL))) )   )
    END IF
  !$mnh_end_do() ! CONCURRENT
  !$acc end parallel
#else
!$acc parallel present_cr(GWORK)
  !$mnh_do_concurrent ( JL=1:JLU )
    IF ( GWORK(JL) ) THEN
       ZZW1(JL,5) = MIN( PRIS(JL),                                    &
            ZZW1(JL,2) / (XCOLIG*BR_EXP(XCOLEXIG*(PZT(JL)-XTT)) ) ) ! RIWETG
       ZZW1(JL,6) = MIN( PRSS(JL),                                    &
            ZZW1(JL,3) / (XCOLSG*BR_EXP(XCOLEXSG*(PZT(JL)-XTT)) ) ) ! RSWETG
       !
       ZZW(JL) = PRVT(JL)*PPRES(JL)/((XMV/XMD)+PRVT(JL)) ! Vapor pressure
       ZZW(JL) =   PKA(JL)*(XTT-PZT(JL)) +                              &
            ( PDV(JL)*(XLVTT + ( XCPV - XCL ) * ( PZT(JL) - XTT )) &
            *(XESTT-ZZW(JL))/(XRV*PZT(JL))           )
       !
       ! compute RWETG
       !
       PRWETG(JL)=MAX( 0.0,                                               &
            ( ZZW(JL) * ( X0DEPG*       BR_POW(PLBDAG(JL),XEX0DEPG) +     &
            X1DEPG*PCJ(JL)*BR_POW(PLBDAG(JL),XEX1DEPG) ) +   &
            ( ZZW1(JL,5)+ZZW1(JL,6) ) *                            &
            ( PRHODREF(JL)*(XLMTT+(XCI-XCL)*(XTT-PZT(JL)))   ) ) / &
            ( PRHODREF(JL)*(XLMTT-XCL*(XTT-PZT(JL))) )   )
    END IF
  !$mnh_end_do() ! CONCURRENT
  !$acc end parallel
#endif 
!
!*       6.4    Select Wet or Dry case
!
!$acc kernels present_cr(GWORK)
  IF     ( KRR == 7 ) THEN
!$mnh_do_concurrent (JL=1:JLU)     
   GWORK(JL) = PRGT(JL)>XRTMIN(6) .AND. PZT(JL)<XTT .and. PRDRYG(JL)>=PRWETG(JL) .AND. PRWETG(JL)>0.0 ! Wet case
   IF( GWORK(JL) )THEN
     ZZW(JL) = PRWETG(JL) - ZZW1(JL,5) - ZZW1(JL,6) ! RCWETG+RRWETG
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
    ZZW1(JL,7) = MAX( 0.0,MIN( ZZW(JL),PRRS(JL)+ZZW1(JL,1) ) )
    PUSW(JL)   = ZZW1(JL,7) / ZZW(JL)
    ZZW1(JL,5) = ZZW1(JL,5)*PUSW(JL)
    ZZW1(JL,6) = ZZW1(JL,6)*PUSW(JL)
    PRWETG(JL) = ZZW1(JL,7) + ZZW1(JL,5) + ZZW1(JL,6)
!
    PRCS(JL) = PRCS(JL) - ZZW1(JL,1)
    PRIS(JL) = PRIS(JL) - ZZW1(JL,5)
    PRSS(JL) = PRSS(JL) - ZZW1(JL,6)
!
! assume a linear percent of conversion of graupel into hail
!
    PRGS(JL) = PRGS(JL) + PRWETG(JL)                     !     Wet growth
    ZZW(JL)  = PRGS(JL)*PRDRYG(JL)/(PRWETG(JL)+PRDRYG(JL)) !        and
    PRGS(JL) = PRGS(JL) - ZZW(JL)                        !   partial conversion
    PRHS(JL) = PRHS(JL) + ZZW(JL)                        ! of the graupel into hail
!
    PRRS(JL) = MAX( 0.0,PRRS(JL) - ZZW1(JL,7) + ZZW1(JL,1) )
    PTHS(JL) = PTHS(JL) + ZZW1(JL,7)*(PLSFACT(JL)-PLVFACT(JL))
                                                 ! f(L_f*(RCWETG+RRWETG))
   ENDIF
!$mnh_end_do()
ELSE IF( KRR == 6 ) THEN
 !$mnh_expand_where(JL=1:JLU)
   GWORK(:) = PRGT(:)>XRTMIN(6) .AND. PZT(:)<XTT .AND. PRDRYG(:)>=PRWETG(:) .AND. PRWETG(:)>0.0 ! Wet case
    WHERE( GWORK(:) )
    PRCS(:) = PRCS(:) - ZZW1(:,1)
    PRIS(:) = PRIS(:) - ZZW1(:,5)
    PRSS(:) = PRSS(:) - ZZW1(:,6)
    PRGS(:) = PRGS(:) + PRWETG(:)
!
    PRRS(:) = PRRS(:) - PRWETG(:) + ZZW1(:,5) + ZZW1(:,6) + ZZW1(:,1)
    PTHS(:) = PTHS(:) + (PRWETG(:)-ZZW1(:,5)-ZZW1(:,6))*(PLSFACT(:)-PLVFACT(:))
                                                 ! f(L_f*(RCWETG+RRWETG))
   END WHERE
  !$mnh_end_expand_where()
END IF
!$acc end kernels

IF (MPPDB_INITIALIZED) THEN
   CALL MPPDB_CHECK(PRRS,"RAIN_ICE_FAST_RG 6.4:PRRS")
END IF

  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'WETG', Unpack ( pths(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'WETG', Unpack ( prcs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'WETG', Unpack ( prrs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'WETG', Unpack ( pris(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'WETG', Unpack ( prss(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'WETG', Unpack ( prgs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'WETG', Unpack ( prhs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )

!PW:used init/end instead of add because zzw1 is produced with a where(...) and is used with other where(...)
! => can not use directly zzw1 in Budget_store_add
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'DRYG', Unpack ( pths(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'DRYG', Unpack ( prcs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'DRYG', Unpack ( prrs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'DRYG', Unpack ( pris(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'DRYG', Unpack ( prss(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'DRYG', Unpack ( prgs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )

!$acc kernels present_cr(GWORK)
!$mnh_expand_where(JL=1:JLU)  
  GWORK(:) = PRGT(:)>XRTMIN(6) .AND. PZT(:)<XTT .AND. PRDRYG(:)<PRWETG(:) .AND. PRDRYG(:)>0.0 ! Dry case

  WHERE( GWORK(:) )
    PRCS(:) = PRCS(:) - ZZW1(:,1)
    PRIS(:) = PRIS(:) - ZZW1(:,2)
    PRSS(:) = PRSS(:) - ZZW1(:,3)
    PRRS(:) = PRRS(:) - ZZW1(:,4)
    PRGS(:) = PRGS(:) + PRDRYG(:)
    PTHS(:) = PTHS(:) + (ZZW1(:,1)+ZZW1(:,4))*(PLSFACT(:)-PLVFACT(:)) !
                      ! f(L_f*(RCDRYG+RRDRYG))
  END WHERE
  !$mnh_end_expand_where()
!$acc end kernels

IF (MPPDB_INITIALIZED) THEN
   CALL MPPDB_CHECK(PRRS,"RAIN_ICE_FAST_RG 6.4b:PRRS")
END IF

  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'DRYG', Unpack ( pths(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'DRYG', Unpack ( prcs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'DRYG', Unpack ( prrs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'DRYG', Unpack ( pris(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'DRYG', Unpack ( prss(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'DRYG', Unpack ( prgs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
!
!      WHERE ( PZT(:) > XTT ) ! RSWETG case only
!        PRSS(:) = PRSS(:) - ZZW1(:,6)
!        PRGS(:) = PRGS(:) + ZZW1(:,6)
!      END WHERE
!
!*       6.5    Melting of the graupeln
!
!$acc kernels present_cr(GWORK)
  GWORK(:) = PRGT(:)>XRTMIN(6) .AND. PRGS(:)>0.0 .AND. PZT(:)>XTT
!$acc end kernels
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
!$acc parallel present_cr(GWORK)
  !$mnh_do_concurrent ( JL=1:JLU )
    IF ( GWORK(JL) ) THEN
      ZZW(JL) = PRVT(JL)*PPRES(JL)/((XMV/XMD)+PRVT(JL)) ! Vapor pressure
      ZZW(JL) =  PKA(JL)*(XTT-PZT(JL)) +                                &
                 ( PDV(JL)*(XLVTT + ( XCPV - XCL ) * ( PZT(JL) - XTT )) &
                             *(XESTT-ZZW(JL))/(XRV*PZT(JL))             )
!
! compute RGMLTR
!
      ZZW(JL)  = MIN( PRGS(JL), MAX( 0.0,( -ZZW(JL) *                   &
                             ( X0DEPG*       PLBDAG(JL)**XEX0DEPG +     &
                               X1DEPG*PCJ(JL)*PLBDAG(JL)**XEX1DEPG ) -  &
                                       ( ZZW1(JL,1)+ZZW1(JL,4) ) *      &
                                ( PRHODREF(JL)*XCL*(XTT-PZT(JL))) ) /   &
                                               ( PRHODREF(JL)*XLMTT ) ) )
      PRRS(JL) = PRRS(JL) + ZZW(JL)
      PRGS(JL) = PRGS(JL) - ZZW(JL)
      PTHS(JL) = PTHS(JL) - ZZW(JL)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(-RGMLTR))
    END IF
  !$mnh_end_do() ! CONCURRENT
!$acc end parallel
#else
!$acc parallel present_cr(GWORK)
  !$mnh_do_concurrent ( JL=1:JLU )
    IF ( GWORK(JL) ) THEN
       ZZW(JL) = PRVT(JL)*PPRES(JL)/((XMV/XMD)+PRVT(JL)) ! Vapor pressure
       ZZW(JL) =  PKA(JL)*(XTT-PZT(JL)) +                                 &
            ( PDV(JL)*(XLVTT + ( XCPV - XCL ) * ( PZT(JL) - XTT )) &
            *(XESTT-ZZW(JL))/(XRV*PZT(JL))             )
       !
       ! compute RGMLTR
       !
       ZZW(JL)  = MIN( PRGS(JL), MAX( 0.0,( -ZZW(JL) *                     &
            ( X0DEPG*       BR_POW(PLBDAG(JL),XEX0DEPG) +     &
            X1DEPG*PCJ(JL)*BR_POW(PLBDAG(JL),XEX1DEPG) ) -   &
            ( ZZW1(JL,1)+ZZW1(JL,4) ) *       &
            ( PRHODREF(JL)*XCL*(XTT-PZT(JL))) ) /    &
            ( PRHODREF(JL)*XLMTT ) ) )
       PRRS(JL) = PRRS(JL) + ZZW(JL)
       PRGS(JL) = PRGS(JL) - ZZW(JL)
       PTHS(JL) = PTHS(JL) - ZZW(JL)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(-RGMLTR))
    END IF
  !$mnh_end_do() ! CONCURRENT
!$acc end parallel
#endif 

IF (MPPDB_INITIALIZED) THEN
   CALL MPPDB_CHECK(PRRS,"RAIN_ICE_FAST_RG 6.5:PRRS")
END IF

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'GMLT', Unpack ( -zzw(:) * ( plsfact(:) - plvfact(:) ) &
                                                                     * prhodj(:), mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'GMLT', Unpack (  zzw(:) * prhodj(:), &
                                                                                  mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'GMLT', Unpack ( -zzw(:) * prhodj(:), &
                                                                                  mask = omicro(:,:,:), field = 0. ) )
IF (MPPDB_INITIALIZED) THEN
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_FAST_RG end:PRCS")
  CALL MPPDB_CHECK(PRRS,"RAIN_ICE_FAST_RG end:PRRS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_FAST_RG end:PRIS")
  CALL MPPDB_CHECK(PRSS,"RAIN_ICE_FAST_RG end:PRSS")
  CALL MPPDB_CHECK(PRGS,"RAIN_ICE_FAST_RG end:PRGS")
  CALL MPPDB_CHECK(PRHS,"RAIN_ICE_FAST_RG end:PRHS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_FAST_RG end:PTHS")
  CALL MPPDB_CHECK(PUSW,"RAIN_ICE_FAST_RG end:PUSW")
  !Check all OUT arrays
  CALL MPPDB_CHECK(PRDRYG,"RAIN_ICE_FAST_RG end:PRDRYG")
  CALL MPPDB_CHECK(PRWETG,"RAIN_ICE_FAST_RG end:PRWETG")
END IF

!$acc end data

#ifdef MNH_OPENACC
!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'RAIN_ICE_FAST_RG 1' )
#endif

!$acc end data

END SUBROUTINE RAIN_ICE_FAST_RG

END MODULE MODE_RAIN_ICE_FAST_RG
