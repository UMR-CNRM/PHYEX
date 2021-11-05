!MNH_LIC Copyright 1995-2020 CNRS, Meteo-France and Universite Paul Sabatier
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
use MODD_CST,            only: XCI, XCL, XCPV, XESTT, XLMTT, XLVTT, XMD, XMV, XRV, XTT
use MODD_RAIN_ICE_DESCR, only: XBG, XBS, XCEXVT, XCXG, XCXH, XCXS, XDH, XLBEXH, XLBH, XRTMIN
use MODD_RAIN_ICE_PARAM, only: NWETLBDAG, NWETLBDAH, NWETLBDAS, X0DEPH, X1DEPH, &
                               XEX0DEPH, XEX1DEPH, XFGWETH, XFSWETH, XFWETH, XKER_GWETH, XKER_SWETH, &
                               XLBGWETH1, XLBGWETH2, XLBGWETH3, XLBSWETH1, XLBSWETH2, XLBSWETH3,     &
                               XWETINTP1G, XWETINTP1H, XWETINTP1S, XWETINTP2G, XWETINTP2H, XWETINTP2S

use mode_budget,         only: Budget_store_add, Budget_store_end, Budget_store_init

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
INTEGER, DIMENSION(size(PRHODREF))   :: I1H, I1W
INTEGER, DIMENSION(:), ALLOCATABLE   :: IVEC1, IVEC2      ! Vectors of indices for interpolations
REAL,    DIMENSION(:), ALLOCATABLE   :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for interpolations
REAL,    DIMENSION(:), ALLOCATABLE   :: ZVECLBDAG, ZVECLBDAH, ZVECLBDAS
REAL,    DIMENSION(size(PRHODREF))   :: ZZW               ! Work array
REAL,    DIMENSION(size(PRHODREF),6) :: ZZW1              ! Work arrays
!
!-------------------------------------------------------------------------------
!
  IHAIL = 0
  DO JJ = 1, SIZE(PRHT)
    IF ( PRHT(JJ)>XRTMIN(7) ) THEN
      IHAIL = IHAIL + 1
      I1H(IHAIL) = JJ
    END IF
  END DO
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
    ZZW1(:,:) = 0.0
!
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
!
!*       7.2.1  accretion of aggregates on the hailstones
!
    IGWET = 0
    DO JJ = 1, IHAIL
      JL = I1H(JJ)
      IF ( PRST(JL)>XRTMIN(5) .AND. PRSS(JL)>0.0 ) THEN
        IGWET = IGWET + 1
        I1W(IGWET) = JL
      END IF
    END DO
!
    IF( IGWET>0 ) THEN
!
!*       7.2.2  allocations
!
      ALLOCATE(ZVECLBDAH(IGWET))
      ALLOCATE(ZVECLBDAS(IGWET))
      ALLOCATE(ZVEC1(IGWET))
      ALLOCATE(ZVEC2(IGWET))
      ALLOCATE(ZVEC3(IGWET))
      ALLOCATE(IVEC1(IGWET))
      ALLOCATE(IVEC2(IGWET))
!
!*       7.2.3  select the (PLBDAH,PLBDAS) couplet
!
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
      DO JJ = 1,IGWET
        ZVEC3(JJ) = (  XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                   * ZVEC1(JJ) &
                   - ( XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
      END DO
!
      DO JJ = 1, IGWET
        JL = I1W(JJ)
        ZZW1(JL,3) = MIN( PRSS(JL),XFSWETH*ZVEC3(JJ)                       & ! RSWETH
                      *( ZVECLBDAS(JJ)**(XCXS-XBS) )*( ZVECLBDAH(JJ)**XCXH )  &
                         *( PRHODREF(JL)**(-XCEXVT-1.) )               &
                         *( XLBSWETH1/( ZVECLBDAH(JJ)**2              ) + &
                            XLBSWETH2/( ZVECLBDAH(JJ)   * ZVECLBDAS(JJ)   ) + &
                            XLBSWETH3/(               ZVECLBDAS(JJ)**2) ) )
      END DO
      DEALLOCATE(ZVECLBDAS)
      DEALLOCATE(ZVECLBDAH)
      DEALLOCATE(IVEC2)
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC3)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(ZVEC1)
    END IF
!
!*       7.2.6  accretion of graupeln on the hailstones
!
    IGWET = 0
    DO JJ = 1, IHAIL
      JL = I1H(JJ)
      IF ( PRGT(JL)>XRTMIN(6) .AND. PRGS(JL)>0.0 ) THEN
        IGWET = IGWET + 1
        I1W(IGWET) = JL
      END IF
    END DO
!
    IF( IGWET>0 ) THEN
!
!*       7.2.7  allocations
!
      ALLOCATE(ZVECLBDAG(IGWET))
      ALLOCATE(ZVECLBDAH(IGWET))
      ALLOCATE(ZVEC1(IGWET))
      ALLOCATE(ZVEC2(IGWET))
      ALLOCATE(ZVEC3(IGWET))
      ALLOCATE(IVEC1(IGWET))
      ALLOCATE(IVEC2(IGWET))
!
!*       7.2.8  select the (PLBDAH,PLBDAG) couplet
!
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
      DO JJ = 1,IGWET
        ZVEC3(JJ) = (  XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                   * ZVEC1(JJ) &
                  - (  XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
      END DO
!
      DO JJ = 1, IGWET
        JL = I1W(JJ)
        ZZW1(JL,5) = MAX(MIN( PRGS(JL),XFGWETH*ZVEC3(JJ)                       & ! RGWETH
                      *( ZVECLBDAG(JJ)**(XCXG-XBG) )*( ZVECLBDAH(JJ)**XCXH )  &
                         *( PRHODREF(JL)**(-XCEXVT-1.) )               &
                         *( XLBGWETH1/( ZVECLBDAH(JJ)**2              ) + &
                            XLBGWETH2/( ZVECLBDAH(JJ)   * ZVECLBDAG(JJ)   ) + &
                            XLBGWETH3/(               ZVECLBDAG(JJ)**2) ) ),0. )
      END DO
      DEALLOCATE(ZVECLBDAH)
      DEALLOCATE(ZVECLBDAG)
      DEALLOCATE(IVEC2)
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC3)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(ZVEC1)
    END IF
!
!*       7.3    compute the Wet growth of hail
!
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

    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HMLT', Unpack ( pths(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'HMLT', Unpack ( prrs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
    if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'HMLT', Unpack ( prhs(:) * prhodj(:), &
                                              mask = omicro(:,:,:), field = 0. ) )
  END IF
!
END SUBROUTINE RAIN_ICE_FAST_RH

END MODULE MODE_RAIN_ICE_FAST_RH
