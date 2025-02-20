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
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_FAST_RH

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: RAIN_ICE_FAST_RH

CONTAINS

SUBROUTINE RAIN_ICE_FAST_RH(OMICRO, PRHODREF, PRVT, PRCT, PRIT, PRST, PRGT, PRHT, PRHODJ, PPRES, &
                            PZT, PLBDAS, PLBDAG, PLBDAH, PLSFACT, PLVFACT, PCJ, PKA, PDV, &
                            PRCS, PRRS, PRIS, PRSS, PRGS, PRHS, PTHS, PUSW)
!
!*      0. DECLARATIONS
!          ------------
!
use modd_budget,         only: lbudget_th, lbudget_rc, lbudget_rr, lbudget_ri, lbudget_rs, lbudget_rg, lbudget_rh, &
                               NBUDGET_TH, NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, &
                               tbudgets
USE MODD_CST,            only: XCI, XCL, XCPV, XESTT, XLMTT, XLVTT, XMD, XMV, XRV, XTT
USE MODD_RAIN_ICE_DESCR_n, only: XBG, XBS, XCEXVT, XCXG, XCXH, XCXS, XDH, XLBEXH, XLBH, XRTMIN
USE MODD_RAIN_ICE_PARAM_n, only: NWETLBDAG, NWETLBDAH, NWETLBDAS, X0DEPH, X1DEPH, &
                               XEX0DEPH, XEX1DEPH, XFGWETH, XFSWETH, XFWETH, XKER_GWETH, XKER_SWETH, &
                               XLBGWETH1, XLBGWETH2, XLBGWETH3, XLBSWETH1, XLBSWETH2, XLBSWETH3,     &
                               XWETINTP1G, XWETINTP1H, XWETINTP1S, XWETINTP2G, XWETINTP2H, XWETINTP2S

use mode_budget,         only: Budget_store_add, Budget_store_end, Budget_store_init
#ifdef MNH_OPENACC
USE MODE_MNH_ZWORK, ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
#endif
use mode_mppdb
use mode_msg
#ifndef MNH_OPENACC
use mode_tools,                        only: Countjv
#else
use mode_tools,                        only: Countjv_device
#endif

#ifdef MNH_BITREP
USE MODI_BITREP
#endif

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,  DIMENSION(:,:,:), intent(in)    :: OMICRO   ! Test where to compute all processes
REAL,     DIMENSION(:),     intent(in)    :: PRHODREF ! RHO Dry REFerence
REAL,     DIMENSION(:),     intent(in)    :: PRVT     ! Water vapor m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRCT     ! Cloud water m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRIT     ! Pristine ice m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRST     ! Snow/aggregate m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRGT     ! Graupel m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRHT    ! Hail m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRHODJ   ! RHO times Jacobian
REAL,     DIMENSION(:),     intent(in)    :: PPRES    ! Pressure
REAL,     DIMENSION(:),     intent(in)    :: PZT      ! Temperature
REAL,     DIMENSION(:),     intent(in)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL,     DIMENSION(:),     intent(in)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL,     DIMENSION(:),     intent(inout) :: PLBDAH   ! Slope parameter of the hail      distribution
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
!
!*       0.2  declaration of local variables
!
INTEGER                              :: IHAIL, IGWET
INTEGER                              :: JJ, JL
#ifndef MNH_OPENACC
INTEGER, DIMENSION(:),   ALLOCATABLE :: I1H, I1W
INTEGER, DIMENSION(:),   ALLOCATABLE :: IVEC1, IVEC2      ! Vectors of indices for interpolations
LOGICAL, DIMENSION(:),   ALLOCATABLE :: GWORK
REAL,    DIMENSION(:),   ALLOCATABLE :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for interpolations
REAL,    DIMENSION(:),   ALLOCATABLE :: ZVECLBDAG, ZVECLBDAH, ZVECLBDAS
REAL,    DIMENSION(:),   ALLOCATABLE :: ZZW               ! Work array
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZZW1              ! Work arrays
#else
INTEGER, DIMENSION(:),   pointer, contiguous :: I1H, I1W
INTEGER, DIMENSION(:),   pointer, contiguous :: IVEC1, IVEC2      ! Vectors of indices for interpolations
LOGICAL, DIMENSION(:),   pointer, contiguous :: GWORK
REAL,    DIMENSION(:),   pointer, contiguous :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for interpolations
REAL,    DIMENSION(:),   pointer, contiguous :: ZVECLBDAG, ZVECLBDAH, ZVECLBDAS
REAL,    DIMENSION(:),   pointer, contiguous :: ZZW               ! Work array
REAL,    DIMENSION(:,:), pointer, contiguous :: ZZW1              ! Work arrays
#endif
!
!-------------------------------------------------------------------------------
!
! IN variables
!
!$acc data present( OMICRO, PRHODREF, PRVT, PRCT, PRIT, PRST, PRGT, PRHT, &
!$acc &             PRHODJ, PPRES, PZT, PLBDAS, PLBDAG, PLSFACT, PLVFACT, &
!$acc &             PCJ, PKA, PDV,                                        &
!
! INOUT variables
!
!$acc &             PLBDAH, PRCS, PRRS, PRIS, PRSS, PRGS, PRHS, PTHS, PUSW )
!
! OUT variables
!
!NONE

#ifdef MNH_OPENACC
CALL PRINT_MSG(NVERB_WARNING,'GEN','RAIN_ICE_FAST_RH','OPENACC: not yet tested')
#endif
!
IF (MPPDB_INITIALIZED) THEN
  !Check all IN arrays
  CALL MPPDB_CHECK(OMICRO,"RAIN_ICE_FAST_RH beg:OMICRO")
  CALL MPPDB_CHECK(PRHODREF,"RAIN_ICE_FAST_RH beg:PRHODREF")
  CALL MPPDB_CHECK(PRVT,"RAIN_ICE_FAST_RH beg:PRVT")
  CALL MPPDB_CHECK(PRCT,"RAIN_ICE_FAST_RH beg:PRCT")
  CALL MPPDB_CHECK(PRIT,"RAIN_ICE_FAST_RH beg:PRIT")
  CALL MPPDB_CHECK(PRST,"RAIN_ICE_FAST_RH beg:PRST")
  CALL MPPDB_CHECK(PRGT,"RAIN_ICE_FAST_RH beg:PRGT")
  CALL MPPDB_CHECK(PRHT,"RAIN_ICE_FAST_RH beg:PRHT")
  CALL MPPDB_CHECK(PRHODJ,"RAIN_ICE_FAST_RH beg:PRHODJ")
  CALL MPPDB_CHECK(PPRES,"RAIN_ICE_FAST_RH beg:PPRES")
  CALL MPPDB_CHECK(PZT,"RAIN_ICE_FAST_RH beg:PZT")
  CALL MPPDB_CHECK(PLBDAS,"RAIN_ICE_FAST_RH beg:PLBDAS")
  CALL MPPDB_CHECK(PLBDAG,"RAIN_ICE_FAST_RH beg:PLBDAG")
  CALL MPPDB_CHECK(PLSFACT,"RAIN_ICE_FAST_RH beg:PLSFACT")
  CALL MPPDB_CHECK(PLVFACT,"RAIN_ICE_FAST_RH beg:PLVFACT")
  CALL MPPDB_CHECK(PCJ,"RAIN_ICE_FAST_RH beg:PCJ")
  CALL MPPDB_CHECK(PKA,"RAIN_ICE_FAST_RH beg:PKA")
  CALL MPPDB_CHECK(PDV,"RAIN_ICE_FAST_RH beg:PDV")
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PLBDAH,"RAIN_ICE_FAST_RH beg:PLBDAH")
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_FAST_RH beg:PRCS")
  CALL MPPDB_CHECK(PRRS,"RAIN_ICE_FAST_RH beg:PRRS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_FAST_RH beg:PRIS")
  CALL MPPDB_CHECK(PRSS,"RAIN_ICE_FAST_RH beg:PRSS")
  CALL MPPDB_CHECK(PRGS,"RAIN_ICE_FAST_RH beg:PRGS")
  CALL MPPDB_CHECK(PRHS,"RAIN_ICE_FAST_RH beg:PRHS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_FAST_RH beg:PTHS")
  CALL MPPDB_CHECK(PUSW,"RAIN_ICE_FAST_RH beg:PUSW")
END IF
!
#ifndef MNH_OPENACC
ALLOCATE( I1H  (size(PRHODREF)) )
ALLOCATE( I1W  (size(PRHODREF)) )
ALLOCATE( GWORK(size(PRHODREF)) )
ALLOCATE( ZZW  (size(PRHODREF)) )
ALLOCATE( ZZW1 (size(PRHODREF),7) )
#else
!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_FAST_RH 1' )

CALL MNH_MEM_GET( I1H,   SIZE(PRHODREF) )
CALL MNH_MEM_GET( I1W,   SIZE(PRHODREF) )
CALL MNH_MEM_GET( GWORK, SIZE(PRHODREF) )
CALL MNH_MEM_GET( ZZW,   SIZE(PRHODREF) )
CALL MNH_MEM_GET( ZZW1,  SIZE(PRHODREF), 7 )

!$acc data present( I1H, I1W, GWORK, ZZW, ZZW1 )
#endif

!$acc kernels present_cr(GWORK)
  GWORK(:) = PRHT(:)>XRTMIN(7)
!$acc end kernels
#ifndef MNH_OPENACC
  IHAIL = COUNTJV( GWORK(:), I1H(:) )
#else
  CALL COUNTJV_DEVICE( GWORK(:), I1H(:), IHAIL )
#endif
!
  IF( IHAIL>0 ) THEN
!PW:used init/end instead of add because zzw1 is produced and used with different conditions
! => can not use directly zzw1 in Budget_store_add
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'WETH', Unpack ( pths(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'WETH', Unpack ( prcs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'WETH', Unpack ( prrs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'WETH', Unpack ( pris(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'WETH', Unpack ( prss(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'WETH', Unpack ( prgs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'WETH', Unpack ( prhs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
!
!*       7.2    compute the Wet growth of hail
!
!$acc kernels
    ZZW1(:,:) = 0.0
!
!$acc loop independent
    DO JJ = 1, IHAIL
      JL = I1H(JJ)
      PLBDAH(JL)  = XLBH * ( PRHODREF(JL) * MAX( PRHT(JL), XRTMIN(7) ) )**XLBEXH

      IF ( PRCT(JL)>XRTMIN(2) .AND. PRCS(JL)>0.0 ) THEN
        ZZW(JL) = PLBDAH(JL)**(XCXH-XDH-2.0) * PRHODREF(JL)**(-XCEXVT)
        ZZW1(JL,1) = MIN( PRCS(JL),XFWETH * PRCT(JL) * ZZW(JL) )             ! RCWETH
      END IF

      IF ( PRIT(JL)>XRTMIN(4) .AND. PRIS(JL)>0.0 ) THEN
        ZZW(JL) = PLBDAH(JL)**(XCXH-XDH-2.0) * PRHODREF(JL)**(-XCEXVT)
        ZZW1(JL,2) = MIN( PRIS(JL),XFWETH * PRIT(JL) * ZZW(JL) )             ! RIWETH
      END IF
    END DO
!$acc end kernels
!
!*       7.2.1  accretion of aggregates on the hailstones
!
!$acc kernels present_cr(GWORK)
    GWORK(1:IHAIL) = PRST(I1H(1:IHAIL))>XRTMIN(5) .AND. PRSS(I1H(1:IHAIL))>0.0
!$acc end kernels
#ifndef MNH_OPENACC
    IGWET = COUNTJV( GWORK(1:IHAIL), I1W(:) )
#else
    CALL COUNTJV_DEVICE( GWORK(1:IHAIL), I1W(:), IGWET )
#endif
!
    IF( IGWET>0 ) THEN
!
!*       7.2.2  allocations
!
#ifndef MNH_OPENACC
      ALLOCATE(ZVECLBDAH(IGWET))
      ALLOCATE(ZVECLBDAS(IGWET))
      ALLOCATE(ZVEC1(IGWET))
      ALLOCATE(ZVEC2(IGWET))
      ALLOCATE(ZVEC3(IGWET))
      ALLOCATE(IVEC1(IGWET))
      ALLOCATE(IVEC2(IGWET))
#else
      !Pin positions in the pools of MNH memory
      CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_FAST_RH 2' )

      CALL MNH_MEM_GET( ZVECLBDAH, IGWET )
      CALL MNH_MEM_GET( ZVECLBDAS, IGWET )
      CALL MNH_MEM_GET( ZVEC1,     IGWET )
      CALL MNH_MEM_GET( ZVEC2,     IGWET )
      CALL MNH_MEM_GET( ZVEC3,     IGWET )
      CALL MNH_MEM_GET( IVEC1,     IGWET )
      CALL MNH_MEM_GET( IVEC2,     IGWET )

!$acc data present( ZVECLBDAH, ZVECLBDAS, ZVEC1, ZVEC2, ZVEC3, IVEC1, IVEC2 )
#endif
!
!*       7.2.3  select the (PLBDAH,PLBDAS) couplet
!
!$acc kernels
      ZVECLBDAH(1:IGWET) = PLBDAH(I1W(1:IGWET))
      ZVECLBDAS(1:IGWET) = PLBDAS(I1W(1:IGWET))
!
!*       7.2.4  find the next lower indice for the PLBDAG and for the PLBDAS
!               in the geometrical set of (Lbda_h,Lbda_s) couplet use to
!               tabulate the SWETH-kernel
!
      ZVEC1(1:IGWET) = MAX( 1.00001, MIN( REAL(NWETLBDAH)-0.00001,           &
                            XWETINTP1H * LOG( ZVECLBDAH(1:IGWET) ) + XWETINTP2H ) )
      IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
      ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - REAL( IVEC1(1:IGWET) )
!
      ZVEC2(1:IGWET) = MAX( 1.00001, MIN( REAL(NWETLBDAS)-0.00001,           &
                            XWETINTP1S * LOG( ZVECLBDAS(1:IGWET) ) + XWETINTP2S ) )
      IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
      ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - REAL( IVEC2(1:IGWET) )
!
!*       7.2.5  perform the bilinear interpolation of the normalized
!               SWETH-kernel
!
!$acc loop independent
      DO JJ = 1,IGWET
        ZVEC3(JJ) = (  XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                   * ZVEC1(JJ) &
                   - ( XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
      END DO
!
!$acc loop independent
      DO JJ = 1, IGWET
        JL = I1W(JJ)
        ZZW1(JL,3) = MIN( PRSS(JL),XFSWETH*ZVEC3(JJ)                       & ! RSWETH
                      *PRST(JL)*( ZVECLBDAH(JJ)**XCXH )  &
                         *( PRHODREF(JL)**(-XCEXVT) )               &
                         *( XLBSWETH1/( ZVECLBDAH(JJ)**2              ) + &
                            XLBSWETH2/( ZVECLBDAH(JJ)   * ZVECLBDAS(JJ)   ) + &
                            XLBSWETH3/(               ZVECLBDAS(JJ)**2) ) )
      END DO
!$acc end kernels

!$acc end data
#ifndef MNH_OPENACC
      DEALLOCATE(ZVECLBDAS)
      DEALLOCATE(ZVECLBDAH)
      DEALLOCATE(IVEC2)
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC3)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(ZVEC1)
#else
      !Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
      CALL MNH_MEM_RELEASE( 'RAIN_ICE_FAST_RH 2' )
#endif
    END IF
!
!*       7.2.6  accretion of graupeln on the hailstones
!
!$acc kernels present_cr(GWORK)
    GWORK(1:IHAIL) = PRGT(I1H(1:IHAIL))>XRTMIN(6) .AND. PRGS(I1H(1:IHAIL))>0.0
!$acc end kernels
#ifndef MNH_OPENACC
    IGWET = COUNTJV( GWORK(1:IHAIL), I1W(:) )
#else
    CALL COUNTJV_DEVICE( GWORK(1:IHAIL), I1W(:), IGWET )
#endif
!
    IF( IGWET>0 ) THEN
!
!*       7.2.7  allocations
!
#ifndef MNH_OPENACC
      ALLOCATE(ZVECLBDAG(IGWET))
      ALLOCATE(ZVECLBDAH(IGWET))
      ALLOCATE(ZVEC1(IGWET))
      ALLOCATE(ZVEC2(IGWET))
      ALLOCATE(ZVEC3(IGWET))
      ALLOCATE(IVEC1(IGWET))
      ALLOCATE(IVEC2(IGWET))
#else
      !Pin positions in the pools of MNH memory
      CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_FAST_RH 3' )

      CALL MNH_MEM_GET( ZVECLBDAG, IGWET )
      CALL MNH_MEM_GET( ZVECLBDAH, IGWET )
      CALL MNH_MEM_GET( ZVEC1,     IGWET )
      CALL MNH_MEM_GET( ZVEC2,     IGWET )
      CALL MNH_MEM_GET( ZVEC3,     IGWET )
      CALL MNH_MEM_GET( IVEC1,     IGWET )
      CALL MNH_MEM_GET( IVEC2,     IGWET )

!$acc data present( ZVECLBDAG, ZVECLBDAH, ZVEC1, ZVEC2, ZVEC3, IVEC1, IVEC2 )
#endif
!
!*       7.2.8  select the (PLBDAH,PLBDAG) couplet
!
!$acc kernels
      ZVECLBDAG(1:IGWET) = PLBDAG(I1W(1:IGWET))
      ZVECLBDAH(1:IGWET) = PLBDAH(I1W(1:IGWET))
!
!*       7.2.9  find the next lower indice for the PLBDAH and for the PLBDAG
!               in the geometrical set of (Lbda_h,Lbda_g) couplet use to
!               tabulate the GWETH-kernel
!
      ZVEC1(1:IGWET) = MAX( 1.00001, MIN( REAL(NWETLBDAG)-0.00001,           &
                            XWETINTP1H * LOG( ZVECLBDAH(1:IGWET) ) + XWETINTP2H ) )
      IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
      ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - REAL( IVEC1(1:IGWET) )
!
      ZVEC2(1:IGWET) = MAX( 1.00001, MIN( REAL(NWETLBDAG)-0.00001,           &
                            XWETINTP1G * LOG( ZVECLBDAG(1:IGWET) ) + XWETINTP2G ) )
      IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
      ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - REAL( IVEC2(1:IGWET) )
!
!*       7.2.10 perform the bilinear interpolation of the normalized
!               GWETH-kernel
!
!$acc loop independent
      DO JJ = 1,IGWET
        ZVEC3(JJ) = (  XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                   * ZVEC1(JJ) &
                  - (  XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
      END DO
!
!$acc loop independent
      DO JJ = 1, IGWET
        JL = I1W(JJ)
        ZZW1(JL,5) = MAX(MIN( PRGS(JL),XFGWETH*ZVEC3(JJ)                       & ! RGWETH
                      *( ZVECLBDAG(JJ)**(XCXG-XBG) )*( ZVECLBDAH(JJ)**XCXH )  &
                         *( PRHODREF(JL)**(-XCEXVT-1.) )               &
                         *( XLBGWETH1/( ZVECLBDAH(JJ)**2              ) + &
                            XLBGWETH2/( ZVECLBDAH(JJ)   * ZVECLBDAG(JJ)   ) + &
                            XLBGWETH3/(               ZVECLBDAG(JJ)**2) ) ),0. )
      END DO
!$acc end kernels

!$acc end data
#ifndef MNH_OPENACC
      DEALLOCATE(ZVECLBDAH)
      DEALLOCATE(ZVECLBDAG)
      DEALLOCATE(IVEC2)
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC3)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(ZVEC1)
#else
      !Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
      CALL MNH_MEM_RELEASE( 'RAIN_ICE_FAST_RH 3' )
#endif
    END IF
!
!*       7.3    compute the Wet growth of hail
!
!$acc kernels
!$acc loop independent
    DO JJ = 1, IHAIL
      JL = I1H(JJ)
      IF ( PZT(JL)<XTT ) THEN
        ZZW(JL) = PRVT(JL)*PPRES(JL)/((XMV/XMD)+PRVT(JL)) ! Vapor pressure
        ZZW(JL) = PKA(JL)*(XTT-PZT(JL)) +                                 &
                  ( PDV(JL)*(XLVTT + ( XCPV - XCL ) * ( PZT(JL) - XTT )) &
                              *(XESTT-ZZW(JL))/(XRV*PZT(JL))             )
!
! compute RWETH
!
        ZZW(JL)  =  MAX(0.,  ( ZZW(JL) * ( X0DEPH*       PLBDAH(JL)**XEX0DEPH +     &
                                  X1DEPH*PCJ(JL)*PLBDAH(JL)**XEX1DEPH ) +   &
                    ( ZZW1(JL,2)+ZZW1(JL,3)+ZZW1(JL,5) ) *                  &
                    ( PRHODREF(JL)*(XLMTT+(XCI-XCL)*(XTT-PZT(JL)))   ) ) / &
                          ( PRHODREF(JL)*(XLMTT-XCL*(XTT-PZT(JL))) ) )

!
        ZZW1(JL,6) = MAX( ZZW(JL) - ZZW1(JL,2) - ZZW1(JL,3) - ZZW1(JL,5),0.) ! RCWETH+RRWETH
        IF ( ZZW1(JL,6)/=0.) THEN
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
          ZZW1(JL,4) = MAX( 0.0,MIN( ZZW1(JL,6),PRRS(JL)+ZZW1(JL,1) ) )
          PUSW(JL)   = ZZW1(JL,4) / ZZW1(JL,6)
          ZZW1(JL,2) = ZZW1(JL,2)*PUSW(JL)
          ZZW1(JL,3) = ZZW1(JL,3)*PUSW(JL)
          ZZW1(JL,5) = ZZW1(JL,5)*PUSW(JL)
          ZZW(JL)    = ZZW1(JL,4) + ZZW1(JL,2) + ZZW1(JL,3) + ZZW1(JL,5)
!
!*       7.1.6  integrate the Wet growth of hail
!
          PRCS(JL) = PRCS(JL) - ZZW1(JL,1)
          PRIS(JL) = PRIS(JL) - ZZW1(JL,2)
          PRSS(JL) = PRSS(JL) - ZZW1(JL,3)
          PRGS(JL) = PRGS(JL) - ZZW1(JL,5)
          PRHS(JL) = PRHS(JL) + ZZW(JL)
          PRRS(JL) = MAX( 0.0,PRRS(JL) - ZZW1(JL,4) + ZZW1(JL,1) )
          PTHS(JL) = PTHS(JL) + ZZW1(JL,4)*(PLSFACT(JL)-PLVFACT(JL))
                           ! f(L_f*(RCWETH+RRWETH))
        END IF
      END IF
    END DO
!$acc end kernels

    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'WETH', Unpack ( pths(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'WETH', Unpack ( prcs(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'WETH', Unpack ( prrs(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'WETH', Unpack ( pris(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'WETH', Unpack ( prss(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'WETH', Unpack ( prgs(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'WETH', Unpack ( prhs(:) * prhodj(:), &
                                             mask = omicro(:,:,:), field = 0. ) )
!
!
! ici LRECONVH et un flag pour autoriser une reconversion partielle de
!la grele en gresil
!
!  IF( IHAIL>0  ) THEN
!
!UPG_CD
!
!
!*       7.45   Conversion of the hailstones into graupel
!
!    XDUMMY6=0.01E-3
!    XDUMMY7=0.001E-3
!    WHERE( PRHT(:)<XDUMMY6 .AND. PRCT(:)<XDUMMY7 .AND. PZT(:)<XTT )
!      ZZW(:) = MIN( 1.0,MAX( 0.0,1.0-(PRCT(:)/XDUMMY7) ) )
!
! assume a linear percent conversion rate of hail into graupel
!
!      ZZW(:)  = PRHS(:)*ZZW(:)
!      PRGS(:) = PRGS(:) + ZZW(:)                      !   partial conversion
!      PRHS(:) = PRHS(:) - ZZW(:)                      ! of hail into graupel
!
!    END WHERE
!  END IF




!
!*       7.5    Melting of the hailstones
!
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HMLT', Unpack ( pths(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'HMLT', Unpack ( prrs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'HMLT', Unpack ( prhs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )

!$acc kernels
!$acc loop independent
    DO JJ = 1, IHAIL
      JL = I1H(JJ)
      IF( PRHS(JL)>0.0 .AND. PZT(JL)>XTT ) THEN
        ZZW(JL) = PRVT(JL)*PPRES(JL)/((XMV/XMD)+PRVT(JL)) ! Vapor pressure
        ZZW(JL) = PKA(JL)*(XTT-PZT(JL)) +                              &
              ( PDV(JL)*(XLVTT + ( XCPV - XCL ) * ( PZT(JL) - XTT )) &
                              *(XESTT-ZZW(JL))/(XRV*PZT(JL))         )
!
! compute RHMLTR
!
        ZZW(JL)  = MIN( PRHS(JL), MAX( 0.0,( -ZZW(JL) *                     &
                              ( X0DEPH*       PLBDAH(JL)**XEX0DEPH +     &
                                X1DEPH*PCJ(JL)*PLBDAH(JL)**XEX1DEPH ) ) / &
                                                ( PRHODREF(JL)*XLMTT ) ) )
        PRRS(JL) = PRRS(JL) + ZZW(JL)
        PRHS(JL) = PRHS(JL) - ZZW(JL)
        PTHS(JL) = PTHS(JL) - ZZW(JL)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(-RHMLTR))
      END IF
    END DO
!$acc end kernels

    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HMLT', Unpack ( pths(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'HMLT', Unpack ( prrs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'HMLT', Unpack ( prhs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
  END IF
!
IF (MPPDB_INITIALIZED) THEN
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PLBDAH,"RAIN_ICE_FAST_RH end:PLBDAH")
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_FAST_RH end:PRCS")
  CALL MPPDB_CHECK(PRRS,"RAIN_ICE_FAST_RH end:PRRS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_FAST_RH end:PRIS")
  CALL MPPDB_CHECK(PRSS,"RAIN_ICE_FAST_RH end:PRSS")
  CALL MPPDB_CHECK(PRGS,"RAIN_ICE_FAST_RH end:PRGS")
  CALL MPPDB_CHECK(PRHS,"RAIN_ICE_FAST_RH end:PRHS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_FAST_RH end:PTHS")
  CALL MPPDB_CHECK(PUSW,"RAIN_ICE_FAST_RH end:PUSW")
END IF

!$acc end data

#ifdef MNH_OPENACC
!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'RAIN_ICE_FAST_RH 1' )
#endif

!$acc end data

END SUBROUTINE RAIN_ICE_FAST_RH

END MODULE MODE_RAIN_ICE_FAST_RH
