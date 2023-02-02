!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_FAST_RH

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_FAST_RH(D, CST, ICEP, ICED, &
                                  KSIZE, KRR, &
                                  GMICRO, &
                                  PTHS, PRHODJ, &
                                  ZRVT, ZRCT, ZRIT, ZRST, ZRGT, ZRHT, &
                                  ZRIS, ZRRS, ZRCS, ZRSS, ZRGS, ZRHS, ZTHS, &
                                  ZRHODREF, ZRHODJ, ZLSFACT, ZLVFACT, &
                                  ZLBDAS, ZLBDAG, ZLBDAH, &
                                  ZCJ, ZKA, ZDV, &
                                  ZZT, ZPRES, &
                                  YDDDH, YDLDDH, YDMDDH)

    USE PARKIND1,            ONLY: JPRB
    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
    USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_T
    USE MODD_CST,            ONLY: CST_T
    USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_T

    USE MODE_BUDGET,         ONLY: BUDGET_DDH
    USE MODD_BUDGET,         ONLY: LBUDGET_TH, LBUDGET_RG, LBUDGET_RR, LBUDGET_RC, &
                                   LBUDGET_RI, LBUDGET_RS, LBUDGET_RH
    USE DDH_MIX,             ONLY: TYP_DDH
    USE YOMLDDH,             ONLY: TLDDH
    USE YOMMDDH,             ONLY: TMDDH

    IMPLICIT NONE

    TYPE(DIMPHYEX_T),       INTENT(IN) :: D
    TYPE(CST_T),            INTENT(IN) :: CST
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED

    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, INTENT(IN) :: KRR

    LOGICAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PTHS    ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRVT     ! Water vapor m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRCT     ! Cloud water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRIT     ! Pristine ice m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRST     ! Snow/aggregate m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRGT     ! Graupel m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRHT     ! Hail m.r. at t

    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRIS     ! Pristine ice m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRRS     ! Rain water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRCS     ! Cloud water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRSS     ! Snow/aggregate m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRGS     ! Graupel m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRHS     ! Hail m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZTHS     ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRHODREF ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRHODJ   ! RHO times Jacobian
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLSFACT  ! L_s/(Pi_ref*C_ph)
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLVFACT  ! L_v/(Pi_ref*C_ph)

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLBDAS   ! Slope parameter of the aggregate distribution
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLBDAG   ! Slope parameter of the graupel   distribution
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZLBDAH   ! Slope parameter of the hail      distribution

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZCJ      ! Function to compute the ventilation coefficient
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZKA      ! Thermal conductivity of the air
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZDV      ! Diffusivity of water vapor in the air

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZZT      ! Temperature
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZPRES    ! Pressure

    TYPE(TYP_DDH),          INTENT(INOUT) :: YDDDH
    TYPE(TLDDH),            INTENT(IN)    :: YDLDDH
    TYPE(TMDDH),            INTENT(IN)    :: YDMDDH

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

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RH',0,ZHOOK_HANDLE)

    GHAIL(:) = ZRHT(:)>ICED%XRTMIN(7)
    IHAIL = COUNT(GHAIL(:))
!
    IF( IHAIL>0 ) THEN
!
!*       7.2    compute the Wet growth of hail
!
      DO JK = 1, KSIZE
        IF (GHAIL(JK)) THEN
          ZLBDAH(JK)  = ICED%XLBH*( ZRHODREF(JK)*MAX( ZRHT(JK),ICED%XRTMIN(7) ) )**ICED%XLBEXH
        END IF
      END DO

      ZZW1(:,:) = 0.0
      DO JK = 1, KSIZE
        IF (GHAIL(JK) .AND. ((ZRCT(JK)>ICED%XRTMIN(2) .AND. ZRCS(JK)>0.0))) THEN
          ZZW(JK) = ZLBDAH(JK)**(ICED%XCXH-ICED%XDH-2.0) * ZRHODREF(JK)**(-ICED%XCEXVT)
          ZZW1(JK,1) = MIN( ZRCS(JK),ICEP%XFWETH * ZRCT(JK) * ZZW(JK) )             ! RCWETH
        END IF
      END DO
      DO JK = 1, KSIZE
        IF (GHAIL(JK) .AND. ((ZRIT(JK)>ICED%XRTMIN(4) .AND. ZRIS(JK)>0.0))) THEN
          ZZW(JK) = ZLBDAH(JK)**(ICED%XCXH-ICED%XDH-2.0) * ZRHODREF(JK)**(-ICED%XCEXVT)
          ZZW1(JK,2) = MIN( ZRIS(JK),ICEP%XFWETH * ZRIT(JK) * ZZW(JK) )             ! RIWETH
        END IF
      END DO
!
!*       7.2.1  accretion of aggregates on the hailstones
!
      GWET(:) = GHAIL(:) .AND. (ZRST(:)>ICED%XRTMIN(5) .AND. ZRSS(:)>0.0)
      IGWET = COUNT( GWET(:) )
!
      IF( IGWET>0 ) THEN
!
!*       7.2.3  select the (ZLBDAH,ZLBDAS) couplet
!
        ZVEC1(:) = PACK( ZLBDAH(:),MASK=GWET(:) )
        ZVEC2(:) = PACK( ZLBDAS(:),MASK=GWET(:) )
!
!*       7.2.4  find the next lower indice for the ZLBDAG and for the ZLBDAS
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
        ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )

        DO JK = 1, KSIZE
          IF (GWET(JK)) THEN
            ZZW1(JK,3) = MIN(ZRSS(JK),ICEP%XFSWETH*ZZW(JK)                       & ! RSWETH
                          *(ZLBDAS(JK)**(ICED%XCXS-ICED%XBS) )*( ZLBDAH(JK)**ICED%XCXH )  &
                          *(ZRHODREF(JK)**(-ICED%XCEXVT-1.) )               &
                          *(ICEP%XLBSWETH1/( ZLBDAH(JK)**2              ) + &
                            ICEP%XLBSWETH2/( ZLBDAH(JK)   * ZLBDAS(JK)   ) + &
                            ICEP%XLBSWETH3/(               ZLBDAS(JK)**2)))
          END IF
        END DO
      END IF
!
!*       7.2.6  accretion of graupeln on the hailstones
!
      GWET(:) = GHAIL(:) .AND. (ZRGT(:)>ICED%XRTMIN(6) .AND. ZRGS(:)>0.0)
      IGWET = COUNT( GWET(:) )
!
      IF( IGWET>0 ) THEN
!
!*       7.2.8  select the (ZLBDAH,ZLBDAG) couplet
!
        ZVEC1(:) = PACK( ZLBDAH(:),MASK=GWET(:) )
        ZVEC2(:) = PACK( ZLBDAG(:),MASK=GWET(:) )
!
!*       7.2.9  find the next lower indice for the ZLBDAH and for the ZLBDAG
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
        ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )

        DO JK = 1, KSIZE
          IF (GWET(JK)) THEN
            ZZW1(JK,5) = MAX(MIN( ZRGS(JK),ICEP%XFGWETH*ZZW(JK)                       & ! RGWETH
                          *( ZLBDAG(JK)**(ICED%XCXG-ICED%XBG) )*( ZLBDAH(JK)**ICED%XCXH )  &
                             *( ZRHODREF(JK)**(-ICED%XCEXVT-1.) )               &
                             *( ICEP%XLBGWETH1/( ZLBDAH(JK)**2              ) + &
                                ICEP%XLBGWETH2/( ZLBDAH(JK)   * ZLBDAG(JK)   ) + &
                                ICEP%XLBGWETH3/(               ZLBDAG(JK)**2) ) ),0. )
          END IF
        END DO
      END IF
!
!*       7.3    compute the Wet growth of hail
!
      ZZW(:) = 0.0
      DO JK = 1, KSIZE
        IF (GHAIL(JK) .AND. ZZT(JK)<CST%XTT) THEN
          ZZW(JK) = ZRVT(JK)*ZPRES(JK)/(CST%XEPSILO+ZRVT(JK)) ! Vapor pressure
          ZZW(JK) = ZKA(JK)*(CST%XTT-ZZT(JK)) +                                 &
                    ( ZDV(JK)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZZT(JK) - CST%XTT )) &
                                *(CST%XESTT-ZZW(JK))/(CST%XRV*ZZT(JK)))
!
!         compute RWETH
!
          ZZW(JK)  =  MAX(0.,  ( ZZW(JK) * ( ICEP%X0DEPH*       ZLBDAH(JK)**ICEP%XEX0DEPH + &
                                    ICEP%X1DEPH*ZCJ(JK)*ZLBDAH(JK)**ICEP%XEX1DEPH ) + &
                       ( ZZW1(JK,2)+ZZW1(JK,3)+ZZW1(JK,5) ) *                  &
                       ( ZRHODREF(JK)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-ZZT(JK))))) / &
                             ( ZRHODREF(JK)*(CST%XLMTT-CST%XCL*(CST%XTT-ZZT(JK))) ) )
!
          ZZW1(JK,6) = MAX(ZZW(JK) - ZZW1(JK,2) - ZZW1(JK,3) - ZZW1(JK,5), 0.) ! RCWETH+RRWETH
        END IF
      END DO

      DO JK = 1, KSIZE
        IF (GHAIL(JK) .AND. ZZT(JK)<CST%XTT  .AND. ZZW1(JK,6)/=0.) THEN
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
          ZZW1(JK,4) = MAX( 0.0,MIN( ZZW1(JK,6),ZRRS(JK)+ZZW1(JK,1) ) )
          ZUSW(JK)   = ZZW1(JK,4) / ZZW1(JK,6)
          ZZW1(JK,2) = ZZW1(JK,2)*ZUSW(JK)
          ZZW1(JK,3) = ZZW1(JK,3)*ZUSW(JK)
          ZZW1(JK,5) = ZZW1(JK,5)*ZUSW(JK)
          ZZW(JK)    = ZZW1(JK,4) + ZZW1(JK,2) + ZZW1(JK,3) + ZZW1(JK,5)
!
!*       7.1.6  integrate the Wet growth of hail
!
          ZRCS(JK) = ZRCS(JK) - ZZW1(JK,1)
          ZRIS(JK) = ZRIS(JK) - ZZW1(JK,2)
          ZRSS(JK) = ZRSS(JK) - ZZW1(JK,3)
          ZRGS(JK) = ZRGS(JK) - ZZW1(JK,5)
          ZRHS(JK) = ZRHS(JK) + ZZW(JK)
          ZRRS(JK) = MAX( 0.0,ZRRS(JK) - ZZW1(JK,4) + ZZW1(JK,1) )
          ZTHS(JK) = ZTHS(JK) + ZZW1(JK,4)*(ZLSFACT(JK)-ZLVFACT(JK))
                           ! f(L_f*(RCWETH+RRWETH))
        END IF
      END DO
    END IF

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:), &
                                     4,'WETH_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     7,'WETH_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     8,'WETH_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     9,'WETH_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    10,'WETH_BU_RRS',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    11,'WETH_BU_RRG',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RH) CALL BUDGET_DDH(UNPACK(ZRHS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    12,'WETH_BU_RRH',YDDDH, YDLDDH, YDMDDH)

!*       7.45   Conversion of the hailstones into graupel

    IF( IHAIL>0 ) THEN
!
!*       7.5    Melting of the hailstones
!
      ZZW(:) = 0.0
      DO JK = 1, KSIZE
        IF (GHAIL(JK) .AND. (ZRHS(JK)>0.0) .AND. (ZZT(JK)>CST%XTT)) THEN
          ZZW(JK) = ZRVT(JK)*ZPRES(JK)/(CST%XEPSILO+ZRVT(JK)) ! Vapor pressure
          ZZW(JK) = ZKA(JK)*(CST%XTT-ZZT(JK)) +                              &
                 ( ZDV(JK)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZZT(JK) - CST%XTT )) &
                                 *(CST%XESTT-ZZW(JK))/(CST%XRV*ZZT(JK)))
!
! compute RHMLTR
!
          ZZW(JK)  = MIN( ZRHS(JK), MAX( 0.0,( -ZZW(JK) *                     &
                                 ( ICEP%X0DEPH*       ZLBDAH(JK)**ICEP%XEX0DEPH +     &
                                   ICEP%X1DEPH*ZCJ(JK)*ZLBDAH(JK)**ICEP%XEX1DEPH ) -   &
                          ZZW1(JK,6)*( ZRHODREF(JK)*CST%XCL*(CST%XTT-ZZT(JK))) ) /    &
                                                   ( ZRHODREF(JK)*CST%XLMTT)))
          ZRRS(JK) = ZRRS(JK) + ZZW(JK)
          ZRHS(JK) = ZRHS(JK) - ZZW(JK)
          ZTHS(JK) = ZTHS(JK) - ZZW(JK)*(ZLSFACT(JK)-ZLVFACT(JK)) ! f(L_f*(-RHMLTR))
        END IF
      END DO
    END IF

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),&
                                     4,'HMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                     8,'HMLT_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RH) CALL BUDGET_DDH(UNPACK(ZRHS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                    12,'HMLT_BU_RRH',YDDDH, YDLDDH, YDMDDH)

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RH',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_FAST_RH

END MODULE MODE_RAIN_ICE_OLD_FAST_RH
