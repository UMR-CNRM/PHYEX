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
use MODD_CST,            only: XCI, XCL, XCPV, XESTT, XLMTT, XLVTT, XMD, XMV, XRV, XTT
use MODD_RAIN_ICE_DESCR, only: XBS, XCEXVT, XCXG, XCXS, XDG, XRTMIN
use MODD_RAIN_ICE_PARAM, only: NDRYLBDAG, NDRYLBDAR, NDRYLBDAS, X0DEPG, X1DEPG, XCOLEXIG, XCOLEXSG, XCOLIG, XCOLSG, XDRYINTP1G, &
                               XDRYINTP1R, XDRYINTP1S, XDRYINTP2G, XDRYINTP2R, XDRYINTP2S, XEX0DEPG, XEX1DEPG, XEXICFRR,        &
                               XEXRCFRI, XFCDRYG, XFIDRYG, XFRDRYG, XFSDRYG, XICFRR, XKER_RDRYG, XKER_SDRYG, XLBRDRYG1,         &
                               XLBRDRYG2, XLBRDRYG3, XLBSDRYG1, XLBSDRYG2, XLBSDRYG3, XRCFRI

use mode_budget,         only: Budget_store_add, Budget_store_end, Budget_store_init

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
INTEGER, DIMENSION(size(PRHODREF))   :: I1
INTEGER, DIMENSION(:), ALLOCATABLE   :: IVEC1, IVEC2      ! Vectors of indices for interpolations
REAL,    DIMENSION(size(PRHODREF))   :: ZZW               ! Work array
REAL,    DIMENSION(:), ALLOCATABLE   :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for interpolations
REAL,    DIMENSION(:), ALLOCATABLE   :: ZVECLBDAG, ZVECLBDAR, ZVECLBDAS
REAL,    DIMENSION(size(PRHODREF),7) :: ZZW1              ! Work arrays
!
!-------------------------------------------------------------------------------
!
!*       6.1    rain contact freezing
!
  ZZW1(:,:) = 0.0
  WHERE( (PRIT(:)>XRTMIN(4)) .AND. (PRRT(:)>XRTMIN(3)) .AND.  &
                             (PRIS(:)>0.0) .AND. (PRRS(:)>0.0) )
    ZZW1(:,3) = MIN( PRIS(:),XICFRR * PRIT(:)                & ! RICFRRG
                                    * PLBDAR(:)**XEXICFRR    &
                                    * PRHODREF(:)**(-XCEXVT) )
    ZZW1(:,4) = MIN( PRRS(:),XRCFRI * PCIT(:)                & ! RRCFRIG
                                    * PLBDAR(:)**XEXRCFRI    &
                                    * PRHODREF(:)**(-XCEXVT-1.) )
    PRIS(:) = PRIS(:) - ZZW1(:,3)
    PRRS(:) = PRRS(:) - ZZW1(:,4)
    PRGS(:) = PRGS(:) + ZZW1(:,3)+ZZW1(:,4)
    PTHS(:) = PTHS(:) + ZZW1(:,4)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*RRCFRIG)
  END WHERE

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
  ZZW1(:,:) = 0.0
  WHERE( PRGT(:)>XRTMIN(6) .AND. PRCT(:)>XRTMIN(2) .AND. PRCS(:)>0.0 )
    ZZW(:) = PLBDAG(:)**(XCXG-XDG-2.0) * PRHODREF(:)**(-XCEXVT)
    ZZW1(:,1) = MIN( PRCS(:),XFCDRYG * PRCT(:) * ZZW(:) )             ! RCDRYG
  END WHERE
  WHERE( (PRGT(:)>XRTMIN(6)) .AND. PRIT(:)>XRTMIN(4) .AND. PRIS(:)>0.0 )
    ZZW(:) = PLBDAG(:)**(XCXG-XDG-2.0) * PRHODREF(:)**(-XCEXVT)
    ZZW1(:,2) = MIN( PRIS(:),XFIDRYG * EXP( XCOLEXIG*(PZT(:)-XTT) ) &
                                     * PRIT(:) * ZZW(:) )             ! RIDRYG
  END WHERE
!
!*       6.2.1  accretion of aggregates on the graupeln
!
  IGDRY = 0
  DO JJ = 1, SIZE(PRST)
    IF ( PRST(JJ)>XRTMIN(5) .AND. PRGT(JJ)>XRTMIN(6) .AND. PRSS(JJ)>0.0 ) THEN
      IGDRY = IGDRY + 1
      I1(IGDRY) = JJ
    END IF
  END DO

  IF( IGDRY>0 ) THEN
!
!*       6.2.2  allocations
!
    ALLOCATE(ZVECLBDAG(IGDRY))
    ALLOCATE(ZVECLBDAS(IGDRY))
    ALLOCATE(ZVEC1(IGDRY))
    ALLOCATE(ZVEC2(IGDRY))
    ALLOCATE(ZVEC3(IGDRY))
    ALLOCATE(IVEC1(IGDRY))
    ALLOCATE(IVEC2(IGDRY))
!
!*       6.2.3  select the (PLBDAG,PLBDAS) couplet
!
    ZVECLBDAG(1:IGDRY) = PLBDAG(I1(1:IGDRY))
    ZVECLBDAS(1:IGDRY) = PLBDAS(I1(1:IGDRY))
!
!*       6.2.4  find the next lower indice for the PLBDAG and for the PLBDAS
!               in the geometrical set of (Lbda_g,Lbda_s) couplet use to
!               tabulate the SDRYG-kernel
!
    ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAG)-0.00001,           &
                          XDRYINTP1G * LOG( ZVECLBDAG(1:IGDRY) ) + XDRYINTP2G ) )
    IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
    ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - REAL( IVEC1(1:IGDRY) )
!
    ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAS)-0.00001,           &
                          XDRYINTP1S * LOG( ZVECLBDAS(1:IGDRY) ) + XDRYINTP2S ) )
    IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
    ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
!
!*       6.2.5  perform the bilinear interpolation of the normalized
!               SDRYG-kernel
!
    DO JJ = 1,IGDRY
      ZVEC3(JJ) =  (  XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * ZVEC1(JJ) &
                 - (  XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
    END DO
!
    DO JJ = 1, IGDRY
      JL = I1(JJ)
      ZZW1(JL,3) = MIN( PRSS(JL),XFSDRYG*ZVEC3(JJ)                         & ! RSDRYG
                                      * EXP( XCOLEXSG*(PZT(JL)-XTT) )  &
                    *( ZVECLBDAS(JJ)**(XCXS-XBS) )*( ZVECLBDAG(JJ)**XCXG )    &
                    *( PRHODREF(JL)**(-XCEXVT-1.) )                    &
                         *( XLBSDRYG1/( ZVECLBDAG(JJ)**2              ) + &
                            XLBSDRYG2/( ZVECLBDAG(JJ)   * ZVECLBDAS(JJ)   ) + &
                            XLBSDRYG3/(               ZVECLBDAS(JJ)**2) ) )
    END DO
    DEALLOCATE(ZVECLBDAS)
    DEALLOCATE(ZVECLBDAG)
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
  END IF
!
!*       6.2.6  accretion of raindrops on the graupeln
!
  IGDRY = 0
  DO JJ = 1, SIZE(PRRT)
    IF ( PRRT(JJ)>XRTMIN(3) .AND. PRGT(JJ)>XRTMIN(6) .AND. PRRS(JJ)>0.0 ) THEN
      IGDRY = IGDRY + 1
      I1(IGDRY) = JJ
    END IF
  END DO
!
  IF( IGDRY>0 ) THEN
!
!*       6.2.7  allocations
!
    ALLOCATE(ZVECLBDAG(IGDRY))
    ALLOCATE(ZVECLBDAR(IGDRY))
    ALLOCATE(ZVEC1(IGDRY))
    ALLOCATE(ZVEC2(IGDRY))
    ALLOCATE(ZVEC3(IGDRY))
    ALLOCATE(IVEC1(IGDRY))
    ALLOCATE(IVEC2(IGDRY))
!
!*       6.2.8  select the (PLBDAG,PLBDAR) couplet
!
    ZVECLBDAG(1:IGDRY) = PLBDAG(I1(1:IGDRY))
    ZVECLBDAR(1:IGDRY) = PLBDAR(I1(1:IGDRY))
!
!*       6.2.9  find the next lower indice for the PLBDAG and for the PLBDAR
!               in the geometrical set of (Lbda_g,Lbda_r) couplet use to
!               tabulate the RDRYG-kernel
!
    ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAG)-0.00001,           &
                          XDRYINTP1G * LOG( ZVECLBDAG(1:IGDRY) ) + XDRYINTP2G ) )
    IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
    ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - REAL( IVEC1(1:IGDRY) )
!
    ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAR)-0.00001,           &
                          XDRYINTP1R * LOG( ZVECLBDAR(1:IGDRY) ) + XDRYINTP2R ) )
    IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
    ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
!
!*       6.2.10 perform the bilinear interpolation of the normalized
!               RDRYG-kernel
!
    DO JJ = 1,IGDRY
      ZVEC3(JJ) =  (  XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                  * ZVEC1(JJ) &
                 - (  XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
    END DO
!
    DO JJ = 1, IGDRY
      JL = I1(JJ)
      ZZW1(JL,4) = MIN( PRRS(JL),XFRDRYG*ZVEC3(JJ)                    & ! RRDRYG
                        *( ZVECLBDAR(JJ)**(-4) )*( ZVECLBDAG(JJ)**XCXG ) &
                               *( PRHODREF(JL)**(-XCEXVT-1.) )   &
                    *( XLBRDRYG1/( ZVECLBDAG(JJ)**2              ) + &
                       XLBRDRYG2/( ZVECLBDAG(JJ)   * ZVECLBDAR(JJ)   ) + &
                       XLBRDRYG3/(               ZVECLBDAR(JJ)**2) ) )
    END DO
    DEALLOCATE(ZVECLBDAR)
    DEALLOCATE(ZVECLBDAG)
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
  END IF
!
  PRDRYG(:) = ZZW1(:,1) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,4)
!
!*       6.3    compute the Wet growth case
!
  PRWETG(:) = 0.0
  WHERE( PRGT(:)>XRTMIN(6) )
    ZZW1(:,5) = MIN( PRIS(:),                                    &
                ZZW1(:,2) / (XCOLIG*EXP(XCOLEXIG*(PZT(:)-XTT)) ) ) ! RIWETG
    ZZW1(:,6) = MIN( PRSS(:),                                    &
                ZZW1(:,3) / (XCOLSG*EXP(XCOLEXSG*(PZT(:)-XTT)) ) ) ! RSWETG
!
    ZZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
    ZZW(:) =   PKA(:)*(XTT-PZT(:)) +                              &
             ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT )) &
                           *(XESTT-ZZW(:))/(XRV*PZT(:))           )
!
! compute RWETG
!
    PRWETG(:)=MAX( 0.0,                                               &
                 ( ZZW(:) * ( X0DEPG*       PLBDAG(:)**XEX0DEPG +     &
                              X1DEPG*PCJ(:)*PLBDAG(:)**XEX1DEPG ) +   &
                 ( ZZW1(:,5)+ZZW1(:,6) ) *                            &
                 ( PRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-PZT(:)))   ) ) / &
                            ( PRHODREF(:)*(XLMTT-XCL*(XTT-PZT(:))) )   )
  END WHERE
!
!*       6.4    Select Wet or Dry case
!
  IF     ( KRR == 7 ) THEN
   WHERE( PRGT(:)>XRTMIN(6) .AND. PZT(:)<XTT                            &
                                        .AND.                          & ! Wet
                              PRDRYG(:)>=PRWETG(:) .AND. PRWETG(:)>0.0 ) ! case
     ZZW(:) = PRWETG(:) - ZZW1(:,5) - ZZW1(:,6) ! RCWETG+RRWETG
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
    ZZW1(:,7) = MAX( 0.0,MIN( ZZW(:),PRRS(:)+ZZW1(:,1) ) )
    PUSW(:)   = ZZW1(:,7) / ZZW(:)
    ZZW1(:,5) = ZZW1(:,5)*PUSW(:)
    ZZW1(:,6) = ZZW1(:,6)*PUSW(:)
    PRWETG(:) = ZZW1(:,7) + ZZW1(:,5) + ZZW1(:,6)
!
    PRCS(:) = PRCS(:) - ZZW1(:,1)
    PRIS(:) = PRIS(:) - ZZW1(:,5)
    PRSS(:) = PRSS(:) - ZZW1(:,6)
!
! assume a linear percent of conversion of graupel into hail
!
    PRGS(:) = PRGS(:) + PRWETG(:)                     !     Wet growth
    ZZW(:)  = PRGS(:)*PRDRYG(:)/(PRWETG(:)+PRDRYG(:)) !        and
    PRGS(:) = PRGS(:) - ZZW(:)                        !   partial conversion
    PRHS(:) = PRHS(:) + ZZW(:)                        ! of the graupel into hail
!
    PRRS(:) = MAX( 0.0,PRRS(:) - ZZW1(:,7) + ZZW1(:,1) )
    PTHS(:) = PTHS(:) + ZZW1(:,7)*(PLSFACT(:)-PLVFACT(:))
                                                 ! f(L_f*(RCWETG+RRWETG))
   END WHERE
  ELSE IF( KRR == 6 ) THEN
   WHERE( PRGT(:)>XRTMIN(6) .AND. PZT(:)<XTT                            &
                                        .AND.                          & ! Wet
                              PRDRYG(:)>=PRWETG(:) .AND. PRWETG(:)>0.0 ) ! case
    PRCS(:) = PRCS(:) - ZZW1(:,1)
    PRIS(:) = PRIS(:) - ZZW1(:,5)
    PRSS(:) = PRSS(:) - ZZW1(:,6)
    PRGS(:) = PRGS(:) + PRWETG(:)
!
    PRRS(:) = PRRS(:) - PRWETG(:) + ZZW1(:,5) + ZZW1(:,6) + ZZW1(:,1)
    PTHS(:) = PTHS(:) + (PRWETG(:)-ZZW1(:,5)-ZZW1(:,6))*(PLSFACT(:)-PLVFACT(:))
                                                 ! f(L_f*(RCWETG+RRWETG))
   END WHERE
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

  WHERE( PRGT(:)>XRTMIN(6) .AND. PZT(:)<XTT                            &
                                        .AND.                          &
                               PRDRYG(:)<PRWETG(:) .AND. PRDRYG(:)>0.0 ) ! Dry
    PRCS(:) = PRCS(:) - ZZW1(:,1)
    PRIS(:) = PRIS(:) - ZZW1(:,2)
    PRSS(:) = PRSS(:) - ZZW1(:,3)
    PRRS(:) = PRRS(:) - ZZW1(:,4)
    PRGS(:) = PRGS(:) + PRDRYG(:)
    PTHS(:) = PTHS(:) + (ZZW1(:,1)+ZZW1(:,4))*(PLSFACT(:)-PLVFACT(:)) !
                      ! f(L_f*(RCDRYG+RRDRYG))
  END WHERE

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
  zzw(:) = 0. !initialization necessary (for budgets)
  WHERE( PRGT(:)>XRTMIN(6) .AND. PRGS(:)>0.0 .AND. PZT(:)>XTT )
    ZZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
    ZZW(:) =  PKA(:)*(XTT-PZT(:)) +                                 &
               ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT )) &
                           *(XESTT-ZZW(:))/(XRV*PZT(:))             )
!
! compute RGMLTR
!
    ZZW(:)  = MIN( PRGS(:), MAX( 0.0,( -ZZW(:) *                     &
                           ( X0DEPG*       PLBDAG(:)**XEX0DEPG +     &
                             X1DEPG*PCJ(:)*PLBDAG(:)**XEX1DEPG ) -   &
                                     ( ZZW1(:,1)+ZZW1(:,4) ) *       &
                              ( PRHODREF(:)*XCL*(XTT-PZT(:))) ) /    &
                                             ( PRHODREF(:)*XLMTT ) ) )
    PRRS(:) = PRRS(:) + ZZW(:)
    PRGS(:) = PRGS(:) - ZZW(:)
    PTHS(:) = PTHS(:) - ZZW(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(-RGMLTR))
  END WHERE

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'GMLT', Unpack ( -zzw(:) * ( plsfact(:) - plvfact(:) ) &
                                                                     * prhodj(:), mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'GMLT', Unpack (  zzw(:) * prhodj(:), &
                                                                                  mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'GMLT', Unpack ( -zzw(:) * prhodj(:), &
                                                                                  mask = omicro(:,:,:), field = 0. ) )
END SUBROUTINE RAIN_ICE_FAST_RG

END MODULE MODE_RAIN_ICE_FAST_RG
