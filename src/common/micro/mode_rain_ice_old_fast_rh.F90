!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_FAST_RH

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_FAST_RH(D, CST, ICEP, ICED, BUCONF, &
                                  KSIZE, KRR, &
                                  GMICRO, &
                                  PTHS, PRHODJ, &
                                  PRVT, PRCT, PRIT, PRST, PRGT, PRHT, &
                                  PRIS, PRRS, PRCS, PRSS, PRGS, PRHS, PZTHS, &
                                  PRHODREF, PZRHODJ, PLSFACT, PLVFACT, &
                                  PLBDAS, PLBDAG, PLBDAH, &
                                  PCJ, PKA, PDV, &
                                  PZT, PRES, &
                                  TBUDGETS, KBUDGETS)

    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_T
    USE MODD_CST,            ONLY: CST_T
    USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_T

    USE MODE_BUDGET_PHY,     ONLY: BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY
    USE MODD_BUDGET,         ONLY: TBUDGETDATA, TBUDGETCONF_t, NBUDGET_TH, NBUDGET_RG, NBUDGET_RR, NBUDGET_RC, &
                                   NBUDGET_RI, NBUDGET_RS, NBUDGET_RH

    IMPLICIT NONE

    TYPE(DIMPHYEX_T),       INTENT(IN) :: D
    TYPE(CST_T),            INTENT(IN) :: CST
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED
    TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF

    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, INTENT(IN) :: KRR

    LOGICAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODJ ! Dry density * Jacobian

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PTHS   ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRGT     ! Graupel m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRHT     ! Hail m.r. at t

    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRIS     ! Pristine ice m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRRS     ! Rain water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRCS     ! Cloud water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRSS     ! Snow/aggregate m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRGS     ! Graupel m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRHS     ! Hail m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PZTHS    ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRHODREF ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PZRHODJ  ! RHO times Jacobian
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: PLBDAH   ! Slope parameter of the hail      distribution

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PKA      ! Thermal conductivity of the air
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PZT      ! Temperature
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRES     ! Pressure

    TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
    INTEGER, INTENT(IN) :: KBUDGETS

    LOGICAL, DIMENSION(KSIZE) :: GWET  ! Test where to compute wet growth
    LOGICAL, DIMENSION(KSIZE) :: GHAIL ! Test where to compute hail growth

    INTEGER, DIMENSION(KSIZE) :: IVEC1 ! Vectors of indices for
    INTEGER, DIMENSION(KSIZE) :: IVEC2 ! Vectors of indices for

    REAL, DIMENSION(KSIZE) :: ZVEC1 ! Work vectors for interpolations
    REAL, DIMENSION(KSIZE) :: ZVEC2 ! Work vectors for interpolations
    REAL, DIMENSION(KSIZE) :: ZVEC3 ! Work vectors for interpolations

    REAL, DIMENSION(KSIZE) :: ZUSW  ! Undersaturation over water

    REAL, DIMENSION(KSIZE)      :: ZZW  ! Work array
    REAL, DIMENSION(KSIZE, KRR) :: ZZW1 ! Work array

    INTEGER :: IGWET, IHAIL
    INTEGER :: JJ, JK

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RH',0,ZHOOK_HANDLE)

    GHAIL(:) = PRHT(:)>ICED%XRTMIN(7)
    IHAIL = COUNT(GHAIL(:))

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'WETH', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'WETH', UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'WETH', UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'WETH', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'WETH', UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'WETH', UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RH), 'WETH', UNPACK(PRHS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
!
    IF( IHAIL>0 ) THEN
!
!*       7.2    compute the Wet growth of hail
!
      DO JK = 1, KSIZE
        IF (GHAIL(JK)) THEN
          PLBDAH(JK)  = ICED%XLBH*( PRHODREF(JK)*MAX( PRHT(JK),ICED%XRTMIN(7) ) )**ICED%XLBEXH
        END IF
      END DO

      ZZW1(:,:) = 0.0
      DO JK = 1, KSIZE
        IF (GHAIL(JK) .AND. ((PRCT(JK)>ICED%XRTMIN(2) .AND. PRCS(JK)>0.0))) THEN
          ZZW(JK) = PLBDAH(JK)**(ICED%XCXH-ICED%XDH-2.0) * PRHODREF(JK)**(-ICED%XCEXVT)
          ZZW1(JK,1) = MIN( PRCS(JK),ICEP%XFWETH * PRCT(JK) * ZZW(JK) )             ! RCWETH
        END IF
      END DO
      DO JK = 1, KSIZE
        IF (GHAIL(JK) .AND. ((PRIT(JK)>ICED%XRTMIN(4) .AND. PRIS(JK)>0.0))) THEN
          ZZW(JK) = PLBDAH(JK)**(ICED%XCXH-ICED%XDH-2.0) * PRHODREF(JK)**(-ICED%XCEXVT)
          ZZW1(JK,2) = MIN( PRIS(JK),ICEP%XFWETH * PRIT(JK) * ZZW(JK) )             ! RIWETH
        END IF
      END DO
!
!*       7.2.1  accretion of aggregates on the hailstones
!
      GWET(:) = GHAIL(:) .AND. (PRST(:)>ICED%XRTMIN(5) .AND. PRSS(:)>0.0)
      IGWET = COUNT( GWET(:) )
!
      IF( IGWET>0 ) THEN
!
!*       7.2.3  select the (PLBDAH,PLBDAS) couplet
!
        ZVEC1(1:IGWET) = PACK( PLBDAH(:),MASK=GWET(:) )
        ZVEC2(1:IGWET) = PACK( PLBDAS(:),MASK=GWET(:) )
!
!*       7.2.4  find the next lower indice for the PLBDAG and for the PLBDAS
!               in the geometrical set of (Lbda_h,Lbda_s) couplet use to
!               tabulate the SWETH-kernel
!
        ZVEC1(1:IGWET) = MAX( 1.00001, MIN( FLOAT(ICEP%NWETLBDAH)-0.00001,           &
                              ICEP%XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + ICEP%XWETINTP2H ) )
        IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
        ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - FLOAT( IVEC1(1:IGWET) )

        ZVEC2(1:IGWET) = MAX( 1.00001, MIN( FLOAT(ICEP%NWETLBDAS)-0.00001,           &
                              ICEP%XWETINTP1S * LOG( ZVEC2(1:IGWET) ) + ICEP%XWETINTP2S ) )
        IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
        ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - FLOAT( IVEC2(1:IGWET) )
!
!*       7.2.5  perform the bilinear interpolation of the normalized
!               SWETH-kernel
!
        DO JJ = 1,IGWET
          ZVEC3(JJ) = (  ICEP%XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                       - ICEP%XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                 * ZVEC1(JJ) &
                     - ( ICEP%XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                       - ICEP%XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                * (ZVEC1(JJ) - 1.0)
        END DO
        ZZW(:) = UNPACK( VECTOR=ZVEC3(1:IGWET),MASK=GWET,FIELD=0.0 )

        DO JK = 1, KSIZE
          IF (GWET(JK)) THEN
            ZZW1(JK,3) = MIN(PRSS(JK),ICEP%XFSWETH*ZZW(JK)                       & ! RSWETH
                          *(PLBDAS(JK)**(ICED%XCXS-ICED%XBS) )*( PLBDAH(JK)**ICED%XCXH )  &
                          *(PRHODREF(JK)**(-ICED%XCEXVT-1.) )               &
                          *(ICEP%XLBSWETH1/( PLBDAH(JK)**2              ) + &
                            ICEP%XLBSWETH2/( PLBDAH(JK)   * PLBDAS(JK)   ) + &
                            ICEP%XLBSWETH3/(               PLBDAS(JK)**2)))
          END IF
        END DO
      END IF
!
!*       7.2.6  accretion of graupeln on the hailstones
!
      GWET(:) = GHAIL(:) .AND. (PRGT(:)>ICED%XRTMIN(6) .AND. PRGS(:)>0.0)
      IGWET = COUNT( GWET(:) )
!
      IF( IGWET>0 ) THEN
!
!*       7.2.8  select the (PLBDAH,PLBDAG) couplet
!
        ZVEC1(1:IGWET) = PACK( PLBDAH(:),MASK=GWET(:) )
        ZVEC2(1:IGWET) = PACK( PLBDAG(:),MASK=GWET(:) )
!
!*       7.2.9  find the next lower indice for the PLBDAH and for the PLBDAG
!               in the geometrical set of (Lbda_h,Lbda_g) couplet use to
!               tabulate the GWETH-kernel
!
        ZVEC1(1:IGWET) = MAX( 1.00001, MIN( FLOAT(ICEP%NWETLBDAG)-0.00001,           &
                              ICEP%XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + ICEP%XWETINTP2H ) )
        IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
        ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - FLOAT( IVEC1(1:IGWET) )

        ZVEC2(1:IGWET) = MAX( 1.00001, MIN( FLOAT(ICEP%NWETLBDAG)-0.00001,           &
                              ICEP%XWETINTP1G * LOG( ZVEC2(1:IGWET) ) + ICEP%XWETINTP2G ) )
        IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
        ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - FLOAT( IVEC2(1:IGWET) )
!
!*       7.2.10 perform the bilinear interpolation of the normalized
!               GWETH-kernel
!
        DO JJ = 1,IGWET
          ZVEC3(JJ) = (  ICEP%XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                       - ICEP%XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                 * ZVEC1(JJ) &
                    - (  ICEP%XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                       - ICEP%XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                * (ZVEC1(JJ) - 1.0)
        END DO
        ZZW(:) = UNPACK( VECTOR=ZVEC3(1:IGWET),MASK=GWET,FIELD=0.0 )

        DO JK = 1, KSIZE
          IF (GWET(JK)) THEN
            ZZW1(JK,5) = MAX(MIN( PRGS(JK),ICEP%XFGWETH*ZZW(JK)                       & ! RGWETH
                          *( PLBDAG(JK)**(ICED%XCXG-ICED%XBG) )*( PLBDAH(JK)**ICED%XCXH )  &
                             *( PRHODREF(JK)**(-ICED%XCEXVT-1.) )               &
                             *( ICEP%XLBGWETH1/( PLBDAH(JK)**2              ) + &
                                ICEP%XLBGWETH2/( PLBDAH(JK)   * PLBDAG(JK)   ) + &
                                ICEP%XLBGWETH3/(               PLBDAG(JK)**2) ) ),0. )
          END IF
        END DO
      END IF
!
!*       7.3    compute the Wet growth of hail
!
      ZZW(:) = 0.0
      DO JK = 1, KSIZE
        IF (GHAIL(JK) .AND. PZT(JK)<CST%XTT) THEN
          ZZW(JK) = PRVT(JK)*PRES(JK)/(CST%XEPSILO+PRVT(JK)) ! Vapor pressure
          ZZW(JK) = PKA(JK)*(CST%XTT-PZT(JK)) +                                 &
                    ( PDV(JK)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PZT(JK) - CST%XTT )) &
                                *(CST%XESTT-ZZW(JK))/(CST%XRV*PZT(JK)))
!
!         compute RWETH
!
          ZZW(JK)  =  MAX(0.,  ( ZZW(JK) * ( ICEP%X0DEPH*       PLBDAH(JK)**ICEP%XEX0DEPH + &
                                    ICEP%X1DEPH*PCJ(JK)*PLBDAH(JK)**ICEP%XEX1DEPH ) + &
                       ( ZZW1(JK,2)+ZZW1(JK,3)+ZZW1(JK,5) ) *                  &
                       ( PRHODREF(JK)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-PZT(JK))))) / &
                             ( PRHODREF(JK)*(CST%XLMTT-CST%XCL*(CST%XTT-PZT(JK))) ) )
!
          ZZW1(JK,6) = MAX(ZZW(JK) - ZZW1(JK,2) - ZZW1(JK,3) - ZZW1(JK,5), 0.) ! RCWETH+RRWETH
        END IF
      END DO

      DO JK = 1, KSIZE
        IF (GHAIL(JK) .AND. PZT(JK)<CST%XTT  .AND. ZZW1(JK,6)/=0.) THEN
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
          ZZW1(JK,4) = MAX( 0.0,MIN( ZZW1(JK,6),PRRS(JK)+ZZW1(JK,1) ) )
          ZUSW(JK)   = ZZW1(JK,4) / ZZW1(JK,6)
          ZZW1(JK,2) = ZZW1(JK,2)*ZUSW(JK)
          ZZW1(JK,3) = ZZW1(JK,3)*ZUSW(JK)
          ZZW1(JK,5) = ZZW1(JK,5)*ZUSW(JK)
          ZZW(JK)    = ZZW1(JK,4) + ZZW1(JK,2) + ZZW1(JK,3) + ZZW1(JK,5)
!
!*       7.1.6  integrate the Wet growth of hail
!
          PRCS(JK) = PRCS(JK) - ZZW1(JK,1)
          PRIS(JK) = PRIS(JK) - ZZW1(JK,2)
          PRSS(JK) = PRSS(JK) - ZZW1(JK,3)
          PRGS(JK) = PRGS(JK) - ZZW1(JK,5)
          PRHS(JK) = PRHS(JK) + ZZW(JK)
          PRRS(JK) = MAX( 0.0,PRRS(JK) - ZZW1(JK,4) + ZZW1(JK,1) )
          PZTHS(JK) = PZTHS(JK) + ZZW1(JK,4)*(PLSFACT(JK)-PLVFACT(JK))
                           ! f(L_f*(RCWETH+RRWETH))
        END IF
      END DO
    END IF

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'WETH', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'WETH', UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'WETH', UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'WETH', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'WETH', UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'WETH', UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RH), 'WETH', UNPACK(PRHS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

!*       7.45   Conversion of the hailstones into graupel

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'HMLT', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'HMLT', UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RH), 'HMLT', UNPACK(PRHS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    IF( IHAIL>0 ) THEN
!
!*       7.5    Melting of the hailstones
!
      ZZW(:) = 0.0
      DO JK = 1, KSIZE
        IF (GHAIL(JK) .AND. (PRHS(JK)>0.0) .AND. (PZT(JK)>CST%XTT)) THEN
          ZZW(JK) = PRVT(JK)*PRES(JK)/(CST%XEPSILO+PRVT(JK)) ! Vapor pressure
          ZZW(JK) = PKA(JK)*(CST%XTT-PZT(JK)) +                              &
                 ( PDV(JK)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PZT(JK) - CST%XTT )) &
                                 *(CST%XESTT-ZZW(JK))/(CST%XRV*PZT(JK)))
!
! compute RHMLTR
!
          ZZW(JK)  = MIN( PRHS(JK), MAX( 0.0,( -ZZW(JK) *                     &
                                 ( ICEP%X0DEPH*       PLBDAH(JK)**ICEP%XEX0DEPH +     &
                                   ICEP%X1DEPH*PCJ(JK)*PLBDAH(JK)**ICEP%XEX1DEPH ) -   &
                          ZZW1(JK,6)*( PRHODREF(JK)*CST%XCL*(CST%XTT-PZT(JK))) ) /    &
                                                   ( PRHODREF(JK)*CST%XLMTT)))
          PRRS(JK) = PRRS(JK) + ZZW(JK)
          PRHS(JK) = PRHS(JK) - ZZW(JK)
          PZTHS(JK) = PZTHS(JK) - ZZW(JK)*(PLSFACT(JK)-PLVFACT(JK)) ! f(L_f*(-RHMLTR))
        END IF
      END DO
    END IF

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'HMLT', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'HMLT', UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RH), 'HMLT', UNPACK(PRHS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RH',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_FAST_RH

END MODULE MODE_RAIN_ICE_OLD_FAST_RH
