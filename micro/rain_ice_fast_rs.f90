!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
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
INTEGER                              :: IGRIM, IGACC
INTEGER                              :: JJ, JL
INTEGER, DIMENSION(size(PRHODREF))   :: I1
INTEGER, DIMENSION(:), ALLOCATABLE   :: IVEC1, IVEC2      ! Vectors of indices for interpolations
REAL,    DIMENSION(size(PRHODREF))   :: ZZW               ! Work array
REAL,    DIMENSION(:), ALLOCATABLE   :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for interpolations
REAL,    DIMENSION(:), ALLOCATABLE   :: ZVECLBDAR, ZVECLBDAS
REAL,    DIMENSION(:), ALLOCATABLE :: ZZW1, ZZW2, ZZW3, ZZW4 ! Work arrays
!-------------------------------------------------------------------------------
!
!*       5.1    cloud droplet riming of the aggregates
!
  IGRIM = 0
  DO JJ = 1, SIZE(PRCT)
    IF ( PRCT(JJ)>XRTMIN(2) .AND. PRST(JJ)>XRTMIN(5) .AND. PRCS(JJ)>0.0 .AND. PZT(JJ)<XTT ) THEN
      IGRIM = IGRIM + 1
      I1(IGRIM) = JJ
    END IF
  END DO
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
    ALLOCATE(ZVECLBDAS(IGRIM))
    ALLOCATE(ZVEC1(IGRIM))
    ALLOCATE(ZVEC2(IGRIM))
    ALLOCATE(IVEC2(IGRIM))
    ALLOCATE(ZZW1(IGRIM))
    ALLOCATE(ZZW2(IGRIM))
    ALLOCATE(ZZW3(IGRIM))
!
!        5.1.1  select the PLBDAS
!
    ZVECLBDAS(1:IGRIM) = PLBDAS(I1(1:IGRIM))
!
!        5.1.2  find the next lower indice for the PLBDAS in the geometrical
!               set of Lbda_s used to tabulate some moments of the incomplete
!               gamma function
!
    ZVEC2(1:IGRIM) = MAX( 1.00001, MIN( REAL(NGAMINC)-0.00001,           &
                          XRIMINTP1 * LOG( ZVECLBDAS(1:IGRIM) ) + XRIMINTP2 ) )
    IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
    ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - REAL( IVEC2(1:IGRIM) )
!
!        5.1.3  perform the linear interpolation of the normalized
!               "2+XDS"-moment of the incomplete gamma function
!
    ZVEC1(1:IGRIM) =   XGAMINC_RIM1( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                     - XGAMINC_RIM1( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
!
!        5.1.4  riming of the small sized aggregates
!
    DO JJ = 1, IGRIM
      JL = I1(JJ)
      ZZW1(JJ) = MIN( PRCS(JL),                                &
                     XCRIMSS * ZVEC1(JJ) * PRCT(JL) * PRST(JL)      & ! RCRIMSS
                                       *   ZVECLBDAS(JJ)**(XBS+XEXCRIMSS) &
                                       * PRHODREF(JL)**(-XCEXVT+1) )
      PRCS(JL) = PRCS(JL) - ZZW1(JJ)
      PRSS(JL) = PRSS(JL) + ZZW1(JJ)
      PTHS(JL) = PTHS(JL) + ZZW1(JJ)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RCRIMSS))
    END DO
!
!        5.1.5  perform the linear interpolation of the normalized
!               "XBS"-moment of the incomplete gamma function
!
    ZVEC1(1:IGRIM) =  XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
!
!        5.1.6  riming-conversion of the large sized aggregates into graupeln
!
!
    DO JJ = 1, IGRIM
      JL = I1(JJ)
      IF ( PRSS(JL) > 0.0 ) THEN
        ZZW2(JJ) = MIN( PRCS(JL),                     &
                    XCRIMSG * PRCT(JL) *PRST(JL)         & ! RCRIMSG
                            *  ZVECLBDAS(JJ)**(XBS+XEXCRIMSG)  &
                            * PRHODREF(JL)**(-XCEXVT+1) &
                            - ZZW1(JJ)              )
        ZZW3(JJ) = MIN( PRSS(JL),                         &
                        PRST(JL) * PRHODREF(JL) * XSRIMCG * ZVECLBDAS(JJ)**(XBS+XEXSRIMCG)   & ! RSRIMCG
                                * (1.0 - ZVEC1(JJ) )/(PTSTEP*PRHODREF(JL)) )
        PRCS(JL) = PRCS(JL) - ZZW2(JJ)
        PRSS(JL) = PRSS(JL) - ZZW3(JJ)
        PRGS(JL) = PRGS(JL) + ZZW2(JJ)+ZZW3(JJ)
        PTHS(JL) = PTHS(JL) + ZZW2(JJ)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RCRIMSG))
      END IF
    END DO

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

    DEALLOCATE(ZZW3)
    DEALLOCATE(ZZW2)
    DEALLOCATE(ZZW1)
    DEALLOCATE(IVEC2)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
    DEALLOCATE(ZVECLBDAS)
  END IF
!
!*       5.2    rain accretion onto the aggregates
!
  IGACC = 0
  DO JJ = 1, SIZE(PRRT)
    IF ( PRRT(JJ)>XRTMIN(3) .AND. PRST(JJ)>XRTMIN(5) .AND. PRRS(JJ)>0.0 .AND. PZT(JJ)<XTT ) THEN
      IGACC = IGACC + 1
      I1(IGACC) = JJ
    END IF
  END DO
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
!
!        5.2.1  select the (PLBDAS,PLBDAR) couplet
!
    ZVECLBDAS(1:IGACC) = PLBDAS(I1(1:IGACC))
    ZVECLBDAR(1:IGACC) = PLBDAR(I1(1:IGACC))
!
!        5.2.2  find the next lower indice for the PLBDAS and for the PLBDAR
!               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
!               tabulate the RACCSS-kernel
!
    ZVEC1(1:IGACC) = MAX( 1.00001, MIN( REAL(NACCLBDAS)-0.00001,           &
                          XACCINTP1S * LOG( ZVECLBDAS(1:IGACC) ) + XACCINTP2S ) )
    IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
    ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - REAL( IVEC1(1:IGACC) )
!
    ZVEC2(1:IGACC) = MAX( 1.00001, MIN( REAL(NACCLBDAR)-0.00001,           &
                          XACCINTP1R * LOG( ZVECLBDAR(1:IGACC) ) + XACCINTP2R ) )
    IVEC2(1:IGACC) = INT( ZVEC2(1:IGACC) )
    ZVEC2(1:IGACC) = ZVEC2(1:IGACC) - REAL( IVEC2(1:IGACC) )
!
!        5.2.3  perform the bilinear interpolation of the normalized
!               RACCSS-kernel
!
    DO JJ = 1,IGACC
      ZVEC3(JJ) =  (  XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * ZVEC1(JJ) &
                 - (  XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
    END DO
!
!        5.2.4  raindrop accretion on the small sized aggregates
!
    DO JJ = 1, IGACC
      JL = I1(JJ)
      ZZW2(JJ) =                                            & !! coef of RRACCS
              XFRACCSS*( PRST(JL)*ZVECLBDAS(JJ)**XBS )*( PRHODREF(JL)**(-XCEXVT) ) &
         *( XLBRACCS1/((ZVECLBDAS(JJ)**2)               ) +                  &
            XLBRACCS2/( ZVECLBDAS(JJ)    * ZVECLBDAR(JJ)    ) +                  &
            XLBRACCS3/(               (ZVECLBDAR(JJ)**2)) )/ZVECLBDAR(JJ)**4
      ZZW4(JJ) = MIN( PRRS(JL),ZZW2(JJ)*ZVEC3(JJ) )           ! RRACCSS
      PRRS(JL) = PRRS(JL) - ZZW4(JJ)
      PRSS(JL) = PRSS(JL) + ZZW4(JJ)
      PTHS(JL) = PTHS(JL) + ZZW4(JJ)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RRACCSS))
    END DO
!
!        5.2.4b perform the bilinear interpolation of the normalized
!               RACCS-kernel
!
    DO JJ = 1,IGACC
      ZVEC3(JJ) =  (   XKER_RACCS(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    -  XKER_RACCS(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                                   * ZVEC2(JJ) &
                 - (   XKER_RACCS(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    -  XKER_RACCS(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                           * (ZVEC2(JJ) - 1.0)
    END DO
    DO JJ = 1, IGACC
      ZZW2(JJ) = ZZW2(JJ) * ZVEC3(JJ)
    END DO
                                                                       !! RRACCS!
!        5.2.5  perform the bilinear interpolation of the normalized
!               SACCRG-kernel
!
    DO JJ = 1,IGACC
      ZVEC3(JJ) =  (  XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    - XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                          * ZVEC2(JJ) &
                 - (  XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    - XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                          * (ZVEC2(JJ) - 1.0)
    END DO
!
!        5.2.6  raindrop accretion-conversion of the large sized aggregates
!               into graupeln
!
    DO JJ = 1, IGACC
      JL = I1(JJ)
      IF ( PRSS(JL) > 0.0 ) THEN
        ZZW2(JJ) = MAX( MIN( PRRS(JL),ZZW2(JJ)-ZZW4(JJ) ),0.0 )       ! RRACCSG
        IF ( ZZW2(JJ) > 0.0 ) THEN
          ZZW3(JJ) = MIN( PRSS(JL),XFSACCRG*ZVEC3(JJ)*                     & ! RSACCRG
                PRST(JL)*( PRHODREF(JL)**(-XCEXVT) ) &
              *( XLBSACCR1/((ZVECLBDAR(JJ)**2)               ) +           &
                  XLBSACCR2/( ZVECLBDAR(JJ)    * ZVECLBDAS(JJ)    ) +           &
                  XLBSACCR3/(               (ZVECLBDAS(JJ)**2)) )/ZVECLBDAR(JJ) )
          PRRS(JL) = PRRS(JL) - ZZW2(JJ)
          PRSS(JL) = PRSS(JL) - ZZW3(JJ)
          PRGS(JL) = PRGS(JL) + ZZW2(JJ)+ZZW3(JJ)
          PTHS(JL) = PTHS(JL) + ZZW2(JJ)*(PLSFACT(JL)-PLVFACT(JL)) !
                                ! f(L_f*(RRACCSG))
        END IF
      END IF
    END DO

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
  END IF
!
!*       5.3    Conversion-Melting of the aggregates
!
  zzw(:) = 0.
  WHERE( PRST(:)>XRTMIN(5) .AND. PRSS(:)>0.0 .AND. PZT(:)>XTT )
    ZZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
    ZZW(:) =  PKA(:)*(XTT-PZT(:)) +                                 &
               ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT )) &
                           *(XESTT-ZZW(:))/(XRV*PZT(:))             )
!
! compute RSMLT
!
    ZZW(:)  = MIN( PRSS(:), XFSCVMG*MAX( 0.0,( -ZZW(:) * PRST(:) * PRHODREF(:) * &
                           ( X0DEPS*       PLBDAS(:)**(XBS+XEX0DEPS) +     &
                             X1DEPS*PCJ(:)*PLBDAS(:)**(XBS+XEX1DEPS) ) ) / &
                                             ( PRHODREF(:)*XLMTT ) ) )
!
! note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
! because the graupeln produced by this process are still icy!!!
!
    PRSS(:) = PRSS(:) - ZZW(:)
    PRGS(:) = PRGS(:) + ZZW(:)
  END WHERE

  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'CMEL', &
                                           Unpack ( -zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'CMEL', &
                                           Unpack (  zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0. ) )

END SUBROUTINE RAIN_ICE_FAST_RS

END MODULE MODE_RAIN_ICE_FAST_RS
