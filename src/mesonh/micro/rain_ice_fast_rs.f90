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
!  P. Wautelet 19/02/2021: bugfix: RIM and ACC terms for budgets are now correctly stored
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_FAST_RS

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: RAIN_ICE_FAST_RS

CONTAINS

SUBROUTINE RAIN_ICE_FAST_RS(PTSTEP, OMICRO, PRHODREF, PRVT, PRCT, PRRT, PRST, PRHODJ, PPRES, PZT, &
                            PLBDAR, PLBDAS, PLSFACT, PLVFACT, PCJ, PKA, PDV, &
                            PRCS, PRRS, PRSS, PRGS, PTHS)
!
!*      0. DECLARATIONS
!          ------------
!
use modd_budget,         only: lbudget_th, lbudget_rc, lbudget_rr, lbudget_rs, lbudget_rg, &
                               NBUDGET_TH, NBUDGET_RC, NBUDGET_RR, NBUDGET_RS, NBUDGET_RG, &
                               tbudgets
USE MODD_CST,            only: XCL, XCPV, XESTT, XLMTT, XLVTT, XMD, XMV, XRV, XTT
USE MODD_RAIN_ICE_DESCR_n, only: XBS, XCEXVT, XCXS, XRTMIN
USE MODD_RAIN_ICE_PARAM_n, only: NACCLBDAR, NACCLBDAS, NGAMINC, X0DEPS, X1DEPS, XACCINTP1R, XACCINTP1S, XACCINTP2R, XACCINTP2S, &
                               XCRIMSG, XCRIMSS, XEX0DEPS, XEX1DEPS, XEXCRIMSG, XEXCRIMSS, XEXSRIMCG, XFRACCSS,               &
                               XFSACCRG, XFSCVMG, XGAMINC_RIM1, XGAMINC_RIM1, XGAMINC_RIM2, XKER_RACCS,                       &
                               XKER_RACCSS, XKER_SACCRG, XLBRACCS1, XLBRACCS2, XLBRACCS3, XLBSACCR1, XLBSACCR2, XLBSACCR3,    &
                               XRIMINTP1, XRIMINTP2, XSRIMCG

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

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                       intent(in)    :: PTSTEP   ! Double Time step
                                                      ! (single if cold start)
LOGICAL,  DIMENSION(:,:,:), intent(in)    :: OMICRO   ! Test where to compute all processes
REAL,     DIMENSION(:),     intent(in)    :: PRHODREF ! RHO Dry REFerence
REAL,     DIMENSION(:),     intent(in)    :: PRVT     ! Water vapor m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRCT     ! Cloud water m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRRT     ! Rain water m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRST     ! Snow/aggregate m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRHODJ   ! RHO times Jacobian
REAL,     DIMENSION(:),     intent(in)    :: PPRES    ! Pressure
REAL,     DIMENSION(:),     intent(in)    :: PZT      ! Temperature
REAL,     DIMENSION(:),     intent(in)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL,     DIMENSION(:),     intent(in)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL,     DIMENSION(:),     intent(in)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL,     DIMENSION(:),     intent(in)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL,     DIMENSION(:),     intent(in)    :: PCJ      ! Function to compute the ventilation coefficient
REAL,     DIMENSION(:),     intent(in)    :: PKA      ! Thermal conductivity of the air
REAL,     DIMENSION(:),     intent(in)    :: PDV      ! Diffusivity of water vapor in the air
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRCS     ! Cloud water m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRRS     ! Rain water m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRSS     ! Snow/aggregate m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRGS     ! Graupel m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PTHS     ! Theta source
!
!*       0.2  declaration of local variables
!
INTEGER                            :: IGRIM, IGACC
INTEGER                            :: JJ, JL
#ifndef MNH_OPENACC
INTEGER, DIMENSION(:), ALLOCATABLE :: I1
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1, IVEC2      ! Vectors of indices for interpolations
LOGICAL, DIMENSION(:), ALLOCATABLE :: GWORK
REAL,    DIMENSION(:), ALLOCATABLE :: ZZW               ! Work array
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for interpolations
REAL,    DIMENSION(:), ALLOCATABLE :: ZVECLBDAR, ZVECLBDAS
REAL,    DIMENSION(:), ALLOCATABLE :: ZZW1, ZZW2, ZZW3, ZZW4 ! Work arrays
#else
INTEGER, DIMENSION(:), pointer, contiguous :: I1
INTEGER, DIMENSION(:), pointer, contiguous :: IVEC1, IVEC2      ! Vectors of indices for interpolations
LOGICAL, DIMENSION(:), pointer, contiguous :: GWORK
REAL,    DIMENSION(:), pointer, contiguous :: ZZW               ! Work array
REAL,    DIMENSION(:), pointer, contiguous :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for interpolations
REAL,    DIMENSION(:), pointer, contiguous :: ZVECLBDAR, ZVECLBDAS
REAL,    DIMENSION(:), pointer, contiguous :: ZZW1, ZZW2, ZZW3, ZZW4 ! Work arrays
#endif
!
INTEGER                            :: JJU
!
!-------------------------------------------------------------------------------
!
! IN variables
!
!$acc data present( OMICRO, PRHODREF, PRVT, PRCT, PRRT, PRST, PRHODJ, &
!$acc &             PPRES, PZT, PLBDAR, PLBDAS, PLSFACT, PLVFACT,     &
!$acc &             PCJ, PKA, PDV,                                    &
!
! INOUT variables
!
!$acc &             PRCS, PRRS, PRSS, PRGS, PTHS , &
!
! use variable
!$acc &             XGAMINC_RIM1,XGAMINC_RIM2, &
!$acc &             XKER_RACCSS,XKER_RACCS,XKER_SACCRG )
!
! OUT variables
!
!NONE
!
!
IF (MPPDB_INITIALIZED) THEN
  !Check all IN arrays
  CALL MPPDB_CHECK(OMICRO,"RAIN_ICE_FAST_RS beg:OMICRO")
  CALL MPPDB_CHECK(PRHODREF,"RAIN_ICE_FAST_RS beg:PRHODREF")
  CALL MPPDB_CHECK(PRVT,"RAIN_ICE_FAST_RS beg:PRVT")
  CALL MPPDB_CHECK(PRCT,"RAIN_ICE_FAST_RS beg:PRCT")
  CALL MPPDB_CHECK(PRRT,"RAIN_ICE_FAST_RS beg:PRRT")
  CALL MPPDB_CHECK(PRST,"RAIN_ICE_FAST_RS beg:PRST")
  CALL MPPDB_CHECK(PRHODJ,"RAIN_ICE_FAST_RS beg:PRHODJ")
  CALL MPPDB_CHECK(PPRES,"RAIN_ICE_FAST_RS beg:PPRES")
  CALL MPPDB_CHECK(PZT,"RAIN_ICE_FAST_RS beg:PZT")
  CALL MPPDB_CHECK(PLBDAR,"RAIN_ICE_FAST_RS beg:PLBDAR")
  CALL MPPDB_CHECK(PLBDAS,"RAIN_ICE_FAST_RS beg:PLBDAS")
  CALL MPPDB_CHECK(PLSFACT,"RAIN_ICE_FAST_RS beg:PLSFACT")
  CALL MPPDB_CHECK(PLVFACT,"RAIN_ICE_FAST_RS beg:PLVFACT")
  CALL MPPDB_CHECK(PCJ,"RAIN_ICE_FAST_RS beg:PCJ")
  CALL MPPDB_CHECK(PKA,"RAIN_ICE_FAST_RS beg:PKA")
  CALL MPPDB_CHECK(PDV,"RAIN_ICE_FAST_RS beg:PDV")
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_FAST_RS beg:PRCS")
  CALL MPPDB_CHECK(PRRS,"RAIN_ICE_FAST_RS beg:PRRS")
  CALL MPPDB_CHECK(PRSS,"RAIN_ICE_FAST_RS beg:PRSS")
  CALL MPPDB_CHECK(PRGS,"RAIN_ICE_FAST_RS beg:PRGS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_FAST_RS beg:PTHS")
  !Check use variables
  CALL MPPDB_CHECK(XKER_RACCSS,"RAIN_ICE_FAST_RS beg:XKER_RACCSS")
END IF
!
#ifndef MNH_OPENACC
ALLOCATE( I1   (size(PRHODREF)) )
ALLOCATE( GWORK(size(PRHODREF)) )
ALLOCATE( ZZW  (size(PRHODREF)) )
#else
!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_FAST_RS 1' )

CALL MNH_MEM_GET( I1,    SIZE(PRHODREF) )
CALL MNH_MEM_GET( GWORK, SIZE(PRHODREF) )
CALL MNH_MEM_GET( ZZW,   SIZE(PRHODREF) )

!$acc data present( I1, GWORK, ZZW )
#endif

JJU = size(PRHODREF)
!
!*       5.1    cloud droplet riming of the aggregates
!
!$acc kernels present_cr(GWORK)
GWORK(:) = PRCT(:)>XRTMIN(2) .AND. PRST(:)>XRTMIN(5) .AND. PRCS(:)>0.0 .AND. PZT(:)<XTT
!$acc end kernels
#ifndef MNH_OPENACC
IGRIM = COUNTJV( GWORK(:), I1(:) )
#else
CALL COUNTJV_DEVICE( GWORK(:), I1(:), IGRIM )
#endif
  !
  IF( IGRIM>0 ) THEN
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'RIM', Unpack ( pths(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'RIM', Unpack ( prcs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'RIM', Unpack ( prss(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'RIM', Unpack ( prgs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )

!
!        5.1.0  allocations
!
#ifndef MNH_OPENACC
    ALLOCATE(ZVECLBDAS(IGRIM))
    ALLOCATE(ZVEC1(IGRIM))
    ALLOCATE(ZVEC2(IGRIM))
    ALLOCATE(IVEC2(IGRIM))
    ALLOCATE(ZZW1(IGRIM))
    ALLOCATE(ZZW2(IGRIM))
    ALLOCATE(ZZW3(IGRIM))
#else
    !Pin positions in the pools of MNH memory
    CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_FAST_RS 2' )

    CALL MNH_MEM_GET( ZVECLBDAS, IGRIM )
    CALL MNH_MEM_GET( ZVEC1,     IGRIM )
    CALL MNH_MEM_GET( ZVEC2,     IGRIM )
    CALL MNH_MEM_GET( IVEC2,     IGRIM )
    CALL MNH_MEM_GET( ZZW1,      IGRIM )
    CALL MNH_MEM_GET( ZZW2,      IGRIM )
    CALL MNH_MEM_GET( ZZW3,      IGRIM )

!$acc data present( ZVECLBDAS, ZVEC1, ZVEC2, IVEC2, ZZW1, ZZW2, ZZW3 )
#endif
!
!        5.1.1  select the PLBDAS
!
#ifdef MNH_COMPILER_CCE
! loop with indirection don't parallelize with Cray -> Keep do concurrent    
!$mnh_undef(LOOP)
#endif    
!$acc kernels
!$mnh_expand_array(JJ = 1:IGRIM)
    ZVECLBDAS(1:IGRIM) = PLBDAS(I1(1:IGRIM))
!
!        5.1.2  find the next lower indice for the PLBDAS in the geometrical
!               set of Lbda_s used to tabulate some moments of the incomplete
!               gamma function
!
    ZVEC2(1:IGRIM) = MAX( 1.00001, MIN( REAL(NGAMINC)-0.00001,           &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
                          XRIMINTP1 * LOG( ZVECLBDAS(1:IGRIM) ) + XRIMINTP2 ) )
#else
                          XRIMINTP1 * BR_LOG( ZVECLBDAS(1:IGRIM) ) + XRIMINTP2 ) )
#endif
    IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
    ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - REAL( IVEC2(1:IGRIM) )
!$mnh_end_expand_array()
!
!        5.1.3  perform the linear interpolation of the normalized
!               "2+XDS"-moment of the incomplete gamma function
!
    !$mnh_do_concurrent (JJ=1:IGRIM)
       ZVEC1(JJ) =   XGAMINC_RIM1( IVEC2(JJ)+1 )* ZVEC2(JJ)      &
                   - XGAMINC_RIM1( IVEC2(JJ)   )*(ZVEC2(JJ) - 1.0)
    !$mnh_end_do() ! CONCURRENT
!
!        5.1.4  riming of the small sized aggregates
!
    ! acc loop independent , private (JL)
    !$mnh_do_concurrent ( JJ = 1:IGRIM )
      JL = I1(JJ)
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
      ZZW1(JJ) = MIN( PRCS(JL),                           &
                     XCRIMSS * ZVEC1(JJ) * PRCT(JL) * PRST(JL)      & ! RCRIMSS
                             *   ZVECLBDAS(JJ)**(XBS+XEXCRIMSS) &
                             * PRHODREF(JL)**(-XCEXVT+1) )
#else
      ZZW1(JJ) = MIN( PRCS(JL),                                &
                     XCRIMSS * ZVEC1(JJ) * PRCT(JL) * PRST(JL)      & ! RCRIMSS
                             * BR_POW(ZVECLBDAS(JJ),XBS+XEXCRIMSS) &
                             * BR_POW(PRHODREF(JL),-XCEXVT+1) )
#endif
      PRCS(JL) = PRCS(JL) - ZZW1(JJ)
      PRSS(JL) = PRSS(JL) + ZZW1(JJ)
      PTHS(JL) = PTHS(JL) + ZZW1(JJ)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RCRIMSS))
   !$mnh_end_do() ! CONCURRENT
   !
!$acc end kernels   
IF (MPPDB_INITIALIZED) THEN
   CALL MPPDB_CHECK(ZZW1,"RAIN_ICE_FAST_RS 5.1.4:ZZW1")   
   CALL MPPDB_CHECK(PRCS,"RAIN_ICE_FAST_RS 5.1.4:PRCS")
   CALL MPPDB_CHECK(PRSS,"RAIN_ICE_FAST_RS 5.1.4:PRSS")
   CALL MPPDB_CHECK(PTHS,"RAIN_ICE_FAST_RS 5.1.4:PTHS")
END IF
!$acc kernels
!
!        5.1.5  perform the linear interpolation of the normalized
!               "XBS"-moment of the incomplete gamma function
!
    !$mnh_do_concurrent (JJ=1:IGRIM)
       ZVEC1(JJ) =  XGAMINC_RIM2( IVEC2(JJ)+1 )* ZVEC2(JJ)      &
                  - XGAMINC_RIM2( IVEC2(JJ)   )*(ZVEC2(JJ) - 1.0)
    !$mnh_end_do() ! CONCURRENT
!
!        5.1.6  riming-conversion of the large sized aggregates into graupeln
!
    !
    ! acc loop independent , private (JL)
    !$mnh_do_concurrent (JJ = 1:IGRIM )
      JL = I1(JJ)
      IF ( PRSS(JL) > 0.0 ) THEN
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
        ZZW2(JJ) = MIN( PRCS(JL),                       &
                    XCRIMSG * PRCT(JL) *PRST(JL)                 & ! RCRIMSG
                            * ZVECLBDAS(JJ)**(XBS+XEXCRIMSG)  &
                            * PRHODREF(JL)**(-XCEXVT+1)   &
                            - ZZW1(JJ) )
        ZZW3(JJ) = MIN( PRSS(JL),                                          &
                        PRST(JL) * PRHODREF(JL) * XSRIMCG * ZVECLBDAS(JJ)**(XBS+XEXSRIMCG)                 & ! RSRIMCG
                                * (1.0 - ZVEC1(JJ) )/(PTSTEP*PRHODREF(JL)) )
#else
        ZZW2(JJ) = MIN( PRCS(JL),                             &
                    XCRIMSG * PRCT(JL) *PRST(JL)                        & ! RCRIMSG
                            * BR_POW(ZVECLBDAS(JJ),XBS+XEXCRIMSG) &
                            * BR_POW(PRHODREF(JL),-XCEXVT+1)    &
                            - ZZW1(JJ) )
        ZZW3(JJ) = MIN( PRSS(JL),                                          &
                        PRST(JL) * PRHODREF(JL) * XSRIMCG * BR_POW(ZVECLBDAS(JJ),XBS+XEXSRIMCG)          & ! RSRIMCG
                                * (1.0 - ZVEC1(JJ) )/(PTSTEP*PRHODREF(JL)) )
#endif
        PRCS(JL) = PRCS(JL) - ZZW2(JJ)
        PRSS(JL) = PRSS(JL) - ZZW3(JJ)
        PRGS(JL) = PRGS(JL) + ZZW2(JJ)+ZZW3(JJ)
        PTHS(JL) = PTHS(JL) + ZZW2(JJ)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RCRIMSG))
      END IF
   !$mnh_end_do() ! CONCURRENT
!$acc end kernels

IF (MPPDB_INITIALIZED) THEN
   CALL MPPDB_CHECK(ZZW1,"RAIN_ICE_FAST_RS 5.1.5:ZZW1")
   CALL MPPDB_CHECK(ZZW2,"RAIN_ICE_FAST_RS 5.1.5:ZZW2")
   CALL MPPDB_CHECK(ZZW3,"RAIN_ICE_FAST_RS 5.1.5:ZZW3")
   CALL MPPDB_CHECK(PRCS,"RAIN_ICE_FAST_RS 5.1.5:PRCS")
   CALL MPPDB_CHECK(PRSS,"RAIN_ICE_FAST_RS 5.1.5:PRSS")
   CALL MPPDB_CHECK(PRGS,"RAIN_ICE_FAST_RS 5.1.5:PRGS")
   CALL MPPDB_CHECK(PTHS,"RAIN_ICE_FAST_RS 5.1.5:PTHS")
END IF
!$acc end data

    !Remark: not possible to use Budget_store_add here
    !        because variables modified a second time but with a if on prss + jl/=jj
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'RIM', Unpack ( pths(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'RIM', Unpack ( prcs(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'RIM', Unpack ( prss(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'RIM', Unpack ( prgs(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )

#ifndef MNH_OPENACC
    DEALLOCATE(ZZW3)
    DEALLOCATE(ZZW2)
    DEALLOCATE(ZZW1)
    DEALLOCATE(IVEC2)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
    DEALLOCATE(ZVECLBDAS)
#else
    !Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
    CALL MNH_MEM_RELEASE( 'RAIN_ICE_FAST_RS 2' )
#endif
  END IF
!
!*       5.2    rain accretion onto the aggregates
!
!$acc kernels present_cr(GWORK)
  GWORK(:) = PRRT(:)>XRTMIN(3) .AND. PRST(:)>XRTMIN(5) .AND. PRRS(:)>0.0 .AND. PZT(:)<XTT
!$acc end kernels
#ifndef MNH_OPENACC
  IGACC = COUNTJV( GWORK(:), I1(:) )
#else
  CALL COUNTJV_DEVICE( GWORK(:), I1(:), IGACC )
#endif
  !
  IF( IGACC>0 ) THEN
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'ACC', Unpack ( pths(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'ACC', Unpack ( prrs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'ACC', Unpack ( prss(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'ACC', Unpack ( prgs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
!
!        5.2.0  allocations
!
#ifndef MNH_OPENACC
    ALLOCATE(ZVECLBDAR(IGACC))
    ALLOCATE(ZVECLBDAS(IGACC))
    ALLOCATE(ZVEC1(IGACC))
    ALLOCATE(ZVEC2(IGACC))
    ALLOCATE(ZVEC3(IGACC))
    ALLOCATE(IVEC1(IGACC))
    ALLOCATE(IVEC2(IGACC))
    ALLOCATE(ZZW2(IGACC))
    ALLOCATE(ZZW3(IGACC))
    ALLOCATE(ZZW4(IGACC))
#else
    !Pin positions in the pools of MNH memory
    CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_FAST_RS 3' )

    CALL MNH_MEM_GET( ZVECLBDAR, IGACC )
    CALL MNH_MEM_GET( ZVECLBDAS, IGACC )
    CALL MNH_MEM_GET( ZVEC1,     IGACC )
    CALL MNH_MEM_GET( ZVEC2,     IGACC )
    CALL MNH_MEM_GET( ZVEC3,     IGACC )
    CALL MNH_MEM_GET( IVEC1,     IGACC )
    CALL MNH_MEM_GET( IVEC2,     IGACC )
    CALL MNH_MEM_GET( ZZW2,      IGACC )
    CALL MNH_MEM_GET( ZZW3,      IGACC )
    CALL MNH_MEM_GET( ZZW4,      IGACC )

!$acc data present( ZVECLBDAR, ZVECLBDAS, ZVEC1, ZVEC2, ZVEC3, IVEC1, IVEC2, ZZW2, ZZW3, ZZW4 )
#endif
!
!        5.2.1  select the (PLBDAS,PLBDAR) couplet
!
!$acc kernels
!$mnh_expand_array(JJ = 1:IGACC)
    ZVECLBDAS(1:IGACC) = PLBDAS(I1(1:IGACC))
    ZVECLBDAR(1:IGACC) = PLBDAR(I1(1:IGACC))
!
!        5.2.2  find the next lower indice for the PLBDAS and for the PLBDAR
!               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
!               tabulate the RACCSS-kernel
!
    ZVEC1(1:IGACC) = MAX( 1.00001, MIN( REAL(NACCLBDAS)-0.00001,           &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
                          XACCINTP1S * LOG( ZVECLBDAS(1:IGACC) ) + XACCINTP2S ) )
#else
                          XACCINTP1S * BR_LOG( ZVECLBDAS(1:IGACC) ) + XACCINTP2S ) )
#endif
    IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
    ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - REAL( IVEC1(1:IGACC) )
!
    ZVEC2(1:IGACC) = MAX( 1.00001, MIN( REAL(NACCLBDAR)-0.00001,           &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
                          XACCINTP1R * LOG( ZVECLBDAR(1:IGACC) ) + XACCINTP2R ) )
#else
                          XACCINTP1R * BR_LOG( ZVECLBDAR(1:IGACC) ) + XACCINTP2R ) )
#endif
    IVEC2(1:IGACC) = INT( ZVEC2(1:IGACC) )
    ZVEC2(1:IGACC) = ZVEC2(1:IGACC) - REAL( IVEC2(1:IGACC) )
!$mnh_end_expand_array()
!
!        5.2.3  perform the bilinear interpolation of the normalized
!               RACCSS-kernel
!
    !$mnh_do_concurrent ( JJ = 1:IGACC )
      ZVEC3(JJ) =  (  XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * ZVEC1(JJ) &
                 - (  XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
    !$mnh_end_do() ! CONCURRENT
!
!        5.2.4  raindrop accretion on the small sized aggregates
!
! acc loop independent , private (JL)
    !$mnh_do_concurrent ( JJ = 1:IGACC )
      JL = I1(JJ)
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
      ZZW2(JJ) =                                                          & !! coef of RRACCS
              XFRACCSS * ( PRST(JL)*ZVECLBDAS(JJ)**XBS )*( PRHODREF(JL)**(-XCEXVT) )  &
         *( XLBRACCS1 / ZVECLBDAS(JJ)**2                                  &
          + XLBRACCS2 / ( ZVECLBDAS(JJ)    * ZVECLBDAR(JJ)    )           &
          + XLBRACCS3 / ZVECLBDAR(JJ)**2 ) / ZVECLBDAR(JJ)**4
#else
      ZZW2(JJ) =                                                                      & !! coef of RRACCS
              XFRACCSS * ( PRST(JL)*BR_POW(ZVECLBDAS(JJ),XBS) * BR_POW(PRHODREF(JL),-XCEXVT)) &
         *( XLBRACCS1 / BR_P2(ZVECLBDAS(JJ))                                          &
          + XLBRACCS2 / ( ZVECLBDAS(JJ)    * ZVECLBDAR(JJ)    )                       &
          + XLBRACCS3 / BR_P2(ZVECLBDAR(JJ)) ) / BR_POW(ZVECLBDAR(JJ),4.0)
#endif
      ZZW4(JJ) = MIN( PRRS(JL),ZZW2(JJ)*ZVEC3(JJ) )           ! RRACCSS
      PRRS(JL) = PRRS(JL) - ZZW4(JJ)
      PRSS(JL) = PRSS(JL) + ZZW4(JJ)
      PTHS(JL) = PTHS(JL) + ZZW4(JJ)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RRACCSS))
   !$mnh_end_do() ! CONCURRENT
!$acc end kernels
IF (MPPDB_INITIALIZED) THEN
   CALL MPPDB_CHECK(ZVEC1,"RAIN_ICE_FAST_RS 5.2.4:IVEC1")
   CALL MPPDB_CHECK(ZVEC2,"RAIN_ICE_FAST_RS 5.2.4:IVEC2")
   CALL MPPDB_CHECK(ZVEC1,"RAIN_ICE_FAST_RS 5.2.4:ZVEC1")
   CALL MPPDB_CHECK(ZVEC2,"RAIN_ICE_FAST_RS 5.2.4:ZVEC2")
   CALL MPPDB_CHECK(ZVEC3,"RAIN_ICE_FAST_RS 5.2.4:ZVEC3")   
   CALL MPPDB_CHECK(ZZW2,"RAIN_ICE_FAST_RS 5.2.4:ZZW2")   
   CALL MPPDB_CHECK(ZZW4,"RAIN_ICE_FAST_RS 5.2.4:ZZW4")   
   CALL MPPDB_CHECK(PRRS,"RAIN_ICE_FAST_RS 5.2.4:PRRS")
END IF   
!$acc kernels    
!
!        5.2.4b perform the bilinear interpolation of the normalized
!               RACCS-kernel
!
    !$mnh_do_concurrent (JJ = 1:IGACC )
      ZVEC3(JJ) =  (   XKER_RACCS(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    -  XKER_RACCS(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                                   * ZVEC2(JJ) &
                 - (   XKER_RACCS(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    -  XKER_RACCS(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                           * (ZVEC2(JJ) - 1.0)
      ZZW2(JJ) = ZZW2(JJ) * ZVEC3(JJ)
    !$mnh_end_do() ! CONCURRENT
                                                                       !! RRACCS!
!        5.2.5  perform the bilinear interpolation of the normalized
!               SACCRG-kernel
    !
    !$mnh_do_concurrent (JJ = 1:IGACC )
      ZVEC3(JJ) =  (  XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    - XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                          * ZVEC2(JJ) &
                 - (  XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    - XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                          * (ZVEC2(JJ) - 1.0)
    !$mnh_end_do() ! CONCURRENT
!
!        5.2.6  raindrop accretion-conversion of the large sized aggregates
!               into graupeln
    !
    ! acc loop independent , private (JL)
    !$mnh_do_concurrent ( JJ = 1:IGACC )
      JL = I1(JJ)
      IF ( PRSS(JL) > 0.0 ) THEN
        ZZW2(JJ) = MAX( MIN( PRRS(JL),ZZW2(JJ)-ZZW4(JJ) ),0.0 )       ! RRACCSG
        IF ( ZZW2(JJ) > 0.0 ) THEN
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
          ZZW3(JJ) = MIN( PRSS(JL),XFSACCRG*ZVEC3(JJ)*                  & ! RSACCRG
                 PRST(JL) * PRHODREF(JL)**(-XCEXVT) &
             * ( XLBSACCR1 / ZVECLBDAR(JJ)**2                           &
               + XLBSACCR2 /( ZVECLBDAR(JJ) * ZVECLBDAS(JJ) )           &
               + XLBSACCR3 / ZVECLBDAS(JJ)**2 ) / ZVECLBDAR(JJ) )
#else
          ZZW3(JJ) = MIN( PRSS(JL),XFSACCRG*ZVEC3(JJ)*                  & ! RSACCRG
                 PRST(JL) * BR_POW(PRHODREF(JL),-XCEXVT) &
             * ( XLBSACCR1 / BR_P2(ZVECLBDAR(JJ))                           &
               + XLBSACCR2 /( ZVECLBDAR(JJ) * ZVECLBDAS(JJ) )           &
               + XLBSACCR3 / BR_P2(ZVECLBDAS(JJ)) ) / ZVECLBDAR(JJ) )
#endif
          PRRS(JL) = PRRS(JL) - ZZW2(JJ)
          PRSS(JL) = PRSS(JL) - ZZW3(JJ)
          PRGS(JL) = PRGS(JL) + ZZW2(JJ)+ZZW3(JJ)
          PTHS(JL) = PTHS(JL) + ZZW2(JJ)*(PLSFACT(JL)-PLVFACT(JL)) !
                                ! f(L_f*(RRACCSG))
        END IF
      END IF
    !$mnh_end_do() ! CONCURRENT
!$acc end kernels

!$acc end data

    !Remark: not possible to use Budget_store_add here
    !        because variables modified a second time but with a if on prss + jl/=jj
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'ACC', Unpack ( pths(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'ACC', Unpack ( prrs(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'ACC', Unpack ( prss(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'ACC', Unpack ( prgs(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )

#ifndef MNH_OPENACC
    DEALLOCATE(ZZW4)
    DEALLOCATE(ZZW3)
    DEALLOCATE(ZZW2)
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
    DEALLOCATE(ZVECLBDAS)
    DEALLOCATE(ZVECLBDAR)
#else
    !Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
    CALL MNH_MEM_RELEASE( 'RAIN_ICE_FAST_RS 3' )
#endif
  END IF
!
!*       5.3    Conversion-Melting of the aggregates
!
!$acc kernels present_cr(GWORK,zzw)
  zzw(:) = 0.
  GWORK(:) = PRST(:)>XRTMIN(5) .AND. PRSS(:)>0.0 .AND. PZT(:)>XTT
!$acc end kernels
!$acc parallel present_cr(GWORK,zzw)  
  !$mnh_do_concurrent (JJ=1:JJU)
   IF ( GWORK(JJ) ) THEN
    ZZW(JJ) = PRVT(JJ)*PPRES(JJ)/((XMV/XMD)+PRVT(JJ)) ! Vapor pressure
    ZZW(JJ) =  PKA(JJ)*(XTT-PZT(JJ)) +                                 &
               ( PDV(JJ)*(XLVTT + ( XCPV - XCL ) * ( PZT(JJ) - XTT )) &
                           *(XESTT-ZZW(JJ))/(XRV*PZT(JJ))             )
!
! compute RSMLT
!
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
    ZZW(JJ)  = MIN( PRSS(JJ), XFSCVMG*MAX( 0.0,( -ZZW(JJ) * PRST(JJ) * PRHODREF(JJ) *     &
                           ( X0DEPS*       PLBDAS(JJ)**(XBS+XEX0DEPS) +     &
                             X1DEPS*PCJ(JJ)*PLBDAS(JJ)**(XBS+XEX1DEPS) ) ) / &
                                             ( PRHODREF(JJ)*XLMTT ) ) )
#else
    ZZW(JJ)  = MIN( PRSS(JJ), XFSCVMG*MAX( 0.0,( -ZZW(JJ) * PRST(JJ) * PRHODREF(JJ) * &
                           ( X0DEPS*       BR_POW(PLBDAS(JJ),XBS+XEX0DEPS) +     &
                             X1DEPS*PCJ(JJ)*BR_POW(PLBDAS(JJ),XBS+XEX1DEPS) ) ) / &
                                             ( PRHODREF(JJ)*XLMTT ) ) )
#endif
!
! note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
! because the graupeln produced by this process are still icy!!!
!
    PRSS(JJ) = PRSS(JJ) - ZZW(JJ)
    PRGS(JJ) = PRGS(JJ) + ZZW(JJ)
  END IF
 !$mnh_end_do() ! CONCURRENT
!$acc end parallel

  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'CMEL', &
                                           Unpack ( -zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'CMEL', &
                                           Unpack (  zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0. ) )

IF (MPPDB_INITIALIZED) THEN
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_FAST_RS end:PRCS")
  CALL MPPDB_CHECK(PRRS,"RAIN_ICE_FAST_RS end:PRRS")
  CALL MPPDB_CHECK(PRSS,"RAIN_ICE_FAST_RS end:PRSS")
  CALL MPPDB_CHECK(PRGS,"RAIN_ICE_FAST_RS end:PRGS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_FAST_RS end:PTHS")
END IF

!$acc end data

#ifdef MNH_OPENACC
!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'RAIN_ICE_FAST_RS 1' )
#endif

!$acc end data

END SUBROUTINE RAIN_ICE_FAST_RS

END MODULE MODE_RAIN_ICE_FAST_RS
