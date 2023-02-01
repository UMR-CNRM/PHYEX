!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_FAST_RG

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_FAST_RG(D, CST, ICEP, ICED,                       &
                                  PTSTEP, KSIZE, KRR,                       &
                                  OCND2, LTIW, GMICRO,                      &
                                  PRHODJ, PTHS,                             &
                                  ZRVT, ZRCT, ZRIT, ZRRT, ZRST, ZRGT, ZCIT, &
                                  ZRIS, ZRRS, ZRCS, ZRSS, ZRGS, ZRHS, ZTHS, &
                                  ZRHODREF, ZRHODJ, ZLSFACT, ZLVFACT,       &
                                  ZCJ, ZKA, ZDV,                            &
                                  ZLBDAR, ZLBDAG, ZLBDAS,                   &
                                  ZTIW, ZZT, ZPRES,                         &
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

    REAL,    INTENT(IN) :: PTSTEP  ! Double Time step
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, INTENT(IN) :: KRR

    LOGICAL, INTENT(IN) :: OCND2
    LOGICAL, INTENT(IN) :: LTIW

    LOGICAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PTHS    ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRVT     ! Water vapor m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRCT     ! Cloud water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRRT     ! Rain water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRIT     ! Pristine ice m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRST     ! Snow/aggregate m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRGT     ! Graupel m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZCIT     ! Pristine ice conc. at t

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

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLBDAR   ! Slope parameter of the raindrop  distribution
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLBDAG   ! Slope parameter of the graupel   distribution
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLBDAS   ! Slope parameter of the aggregate distribution

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZCJ      ! Function to compute the ventilation coefficient
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZKA      ! Thermal conductivity of the air
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZDV      ! Diffusivity of water vapor in the air

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZTIW     ! Wet bulb temperature

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZZT      ! Temperature
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZPRES    ! Pressure

    TYPE(TYP_DDH),          INTENT(INOUT) :: YDDDH
    TYPE(TLDDH),            INTENT(IN)    :: YDLDDH
    TYPE(TMDDH),            INTENT(IN)    :: YDMDDH

    LOGICAL, DIMENSION(KSIZE) :: GDRY ! Test where to compute dry growth

    INTEGER, DIMENSION(KSIZE) :: IVEC1 ! Vectors of indices for
    INTEGER, DIMENSION(KSIZE) :: IVEC2 ! Vectors of indices for

    REAL, DIMENSION(KSIZE) :: ZVEC1 ! Work vectors for interpolations
    REAL, DIMENSION(KSIZE) :: ZVEC2 ! Work vectors for interpolations
    REAL, DIMENSION(KSIZE) :: ZVEC3 ! Work vectors for interpolations

    REAL, DIMENSION(KSIZE) :: ZRDRYG   ! Dry growth rate of the graupeln
    REAL, DIMENSION(KSIZE) :: ZRWETG   ! Wet growth rate of the graupeln

    REAL, DIMENSION(KSIZE) :: ZUSW  ! Undersaturation over water

    REAL, DIMENSION(KSIZE)      :: ZZW  ! Work array
    REAL, DIMENSION(KSIZE, KRR) :: ZZW1 ! Work array

    INTEGER :: IGDRY
    INTEGER :: JJ, JK

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       6.1    rain contact freezing
!
    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RG',0,ZHOOK_HANDLE)

    ZZW1(:,3:4) = 0.0
    DO JK = 1, KSIZE
      IF ((ZRIT(JK) > ICED%XRTMIN(4) .AND. ZRIT(JK)>ICEP%XFRMIN(2)) .AND. &
          (ZRRT(JK) > ICED%XRTMIN(3)) .AND. &
          (ZRIS(JK) > 0.0) .AND. &
          (ZRRS(JK) > 0.0)) THEN
        ZZW1(JK,3) = MIN( ZRIS(JK),ICEP%XICFRR * ZRIT(JK)                & ! RICFRRG
                                        * ZLBDAR(JK)**ICEP%XEXICFRR    &
                                        * ZRHODREF(JK)**(-ICED%XCEXVT) )
        ZZW1(JK,4) = MIN( ZRRS(JK),ICEP%XRCFRI * ZCIT(JK)                & ! RRCFRIG
                                        * ZLBDAR(JK)**ICEP%XEXRCFRI    &
                                        * ZRHODREF(JK)**(-ICED%XCEXVT-1.) )
        ZRIS(JK) = ZRIS(JK) - ZZW1(JK,3)
        ZRRS(JK) = ZRRS(JK) - ZZW1(JK,4)
        ZRGS(JK) = ZRGS(JK) + ZZW1(JK,3)+ZZW1(JK,4)
        ZTHS(JK) = ZTHS(JK) + ZZW1(JK,4)*(ZLSFACT(JK)-ZLVFACT(JK)) ! f(L_f*RRCFRIG)
      END IF
    END DO

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:), &
                                     4,'CFRZ_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                     8,'CFRZ_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                     9,'CFRZ_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                    11,'CFRZ_BU_RRG',YDDDH, YDLDDH, YDMDDH)

!*       6.2    compute the Dry growth case

    ZZW1(:,:) = 0.0
    DO JK = 1, KSIZE
      IF ((ZRGT(JK) > ICED%XRTMIN(6)) .AND. &
         ((ZRCT(JK) > ICED%XRTMIN(2)  .AND. ZRCS(JK) > 0.0))) THEN
        ZZW(JK) = ZLBDAG(JK)**(ICED%XCXG-ICED%XDG-2.0) * ZRHODREF(JK)**(-ICED%XCEXVT)
        ZZW1(JK,1) = MIN( ZRCS(JK),ICEP%XFCDRYG * ZRCT(JK) * ZZW(JK) )             ! RCDRYG
      END IF
    END DO

    DO JK = 1, KSIZE
      IF ((ZRGT(JK) > ICED%XRTMIN(6)) .AND. &
         ((ZRIT(JK) > ICED%XRTMIN(4)  .AND. ZRIS(JK)>0.0))) THEN

        ZZW(JK) = ZLBDAG(JK)**(ICED%XCXG-ICED%XDG-2.0) * ZRHODREF(JK)**(-ICED%XCEXVT)
        ZZW1(JK,2) = MIN(ZRIS(JK), ICEP%XFIDRYG * EXP(ICEP%XCOLEXIG*(ZZT(JK) - CST%XTT)) &
                                                * ZRIT(JK) * ZZW(JK) )             ! RIDRYG
      END IF
    END DO
!
!*       6.2.1  accretion of aggregates on the graupeln
!
    GDRY(:) = (ZRST(:)>ICED%XRTMIN(5)) .AND. (ZRGT(:)>ICED%XRTMIN(6)) .AND. (ZRSS(:)>0.0)
    IGDRY = COUNT( GDRY(:) )
!
    IF( IGDRY>0 ) THEN
!
!*       6.2.3  select the (ZLBDAG,ZLBDAS) couplet
!
      ZVEC1(:) = PACK( ZLBDAG(:),MASK=GDRY(:) )
      ZVEC2(:) = PACK( ZLBDAS(:),MASK=GDRY(:) )
!
!*       6.2.4  find the next lower indice for the ZLBDAG and for the ZLBDAS
!               in the geometrical set of (Lbda_g,Lbda_s) couplet use to
!               tabulate the SDRYG-kernel
!
      ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( FLOAT(ICEP%NDRYLBDAG)-0.00001,           &
                          ICEP%XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + ICEP%XDRYINTP2G ) )
      IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
      ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - FLOAT( IVEC1(1:IGDRY) )
!
      ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( FLOAT(ICEP%NDRYLBDAS)-0.00001,           &
                          ICEP%XDRYINTP1S * LOG( ZVEC2(1:IGDRY) ) + ICEP%XDRYINTP2S ) )
      IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
      ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - FLOAT( IVEC2(1:IGDRY) )
!
!*       6.2.5  perform the bilinear interpolation of the normalized
!               SDRYG-kernel
!
      DO JJ = 1,IGDRY
        ZVEC3(JJ) =  (  ICEP%XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                      - ICEP%XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                           * ZVEC1(JJ) &
                   - (  ICEP%XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                      - ICEP%XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                           * (ZVEC1(JJ) - 1.0)
      END DO
      ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
      IF (OCND2) THEN
        ZZW1(:,3) = 0.
      ELSE
        DO JK = 1, KSIZE
          IF (GDRY(JK)) THEN
            ZZW1(JK,3) = MIN(ZRSS(JK), ICEP%XFSDRYG*ZZW(JK)            & ! RSDRYG
                               * EXP(ICEP%XCOLEXSG*(ZZT(JK)-CST%XTT))  &
                               *(ZLBDAS(JK)**(ICED%XCXS-ICED%XBS))     &
                               *(ZLBDAG(JK)**ICED%XCXG)                &
                               *(ZRHODREF(JK)**(-ICED%XCEXVT-1.))      &
                               *(ICEP%XLBSDRYG1/( ZLBDAG(JK)**2               ) + &
                                 ICEP%XLBSDRYG2/( ZLBDAG(JK)   * ZLBDAS(JK)   ) + &
                                 ICEP%XLBSDRYG3/(                ZLBDAS(JK)**2) ) )
          END IF
        END DO
      ENDIF
    END IF
!
!*       6.2.6  accretion of raindrops on the graupeln
!
    GDRY(:) = (ZRRT(:)>ICED%XRTMIN(3)) .AND. (ZRGT(:)>ICED%XRTMIN(6)) .AND. (ZRRS(:)>0.0)
    IGDRY = COUNT( GDRY(:) )
!
    IF (IGDRY>0) THEN
!
!*       6.2.8  select the (ZLBDAG,ZLBDAR) couplet
!
      ZVEC1(:) = PACK( ZLBDAG(:),MASK=GDRY(:) )
      ZVEC2(:) = PACK( ZLBDAR(:),MASK=GDRY(:) )
!
!*       6.2.9  find the next lower indice for the ZLBDAG and for the ZLBDAR
!               in the geometrical set of (Lbda_g,Lbda_r) couplet use to
!               tabulate the RDRYG-kernel
!
      ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( FLOAT(ICEP%NDRYLBDAG)-0.00001,           &
                            ICEP%XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + ICEP%XDRYINTP2G ) )
      IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
      ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - FLOAT( IVEC1(1:IGDRY) )

      ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( FLOAT(ICEP%NDRYLBDAR)-0.00001,           &
                          ICEP%XDRYINTP1R * LOG( ZVEC2(1:IGDRY) ) + ICEP%XDRYINTP2R ) )
      IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
      ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - FLOAT( IVEC2(1:IGDRY) )
!
!*       6.2.10 perform the bilinear interpolation of the normalized
!               RDRYG-kernel
!
      DO JJ = 1,IGDRY
        ZVEC3(JJ) =  (  ICEP%XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)         &
                      - ICEP%XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0)) &
                                                                * ZVEC1(JJ)         &
                   - (  ICEP%XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)         &
                      - ICEP%XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0)) &
                                                                *(ZVEC1(JJ) - 1.0)
      END DO
      ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )

      DO JK = 1, KSIZE
        IF (GDRY(JK)) THEN
          ZZW1(JK,4) = MIN(ZRRS(JK),ICEP%XFRDRYG*ZZW(JK)                     & ! RRDRYG
                           *(ZLBDAR(JK)**(-4) )*(ZLBDAG(JK)**ICED%XCXG)      &
                           *(ZRHODREF(JK)**(-ICED%XCEXVT-1.))                &
                           *(ICEP%XLBRDRYG1/(ZLBDAG(JK)**2               ) + &
                             ICEP%XLBRDRYG2/(ZLBDAG(JK)   * ZLBDAR(JK)   ) + &
                             ICEP%XLBRDRYG3/(               ZLBDAR(JK)**2) ) )
        END IF
      END DO
    END IF

    ZRDRYG(:) = ZZW1(:,1) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,4)
!
!*       6.3    compute the Wet growth case
!
    ZZW(:) = 0.0
    ZRWETG(:) = 0.0

    DO JK = 1, JK
      IF (ZRGT(JK)>ICED%XRTMIN(6)) THEN
        ZZW1(JK,5) = MIN( ZRIS(JK),                                    &
                    ZZW1(JK,2) / (ICEP%XCOLIG*EXP(ICEP%XCOLEXIG*(ZZT(JK)-CST%XTT)) ) ) ! RIWETG
        ZZW1(JK,6) = MIN( ZRSS(JK),                                    &
                    ZZW1(JK,3) / (ICEP%XCOLSG*EXP(ICEP%XCOLEXSG*(ZZT(JK)-CST%XTT)) ) ) ! RSWETG
!
        ZZW(JK) = ZRVT(JK)*ZPRES(JK)/(CST%XEPSILO+ZRVT(JK)) ! Vapor pressure
        ZZW(JK) = ZKA(JK)*(CST%XTT-ZZT(JK)) +                              &
                (ZDV(JK)*(CST%XLVTT + (CST%XCPV - CST%XCL) * (ZZT(JK) - CST%XTT)) &
                       *(CST%XESTT-ZZW(JK))/(CST%XRV*ZZT(JK))           )
!
! compute RWETG
!
        ZRWETG(JK)=MAX(0.0,                                               &
                     (ZZW(JK) * ( ICEP%X0DEPG*       ZLBDAG(JK)**ICEP%XEX0DEPG  + &
                                 ICEP%X1DEPG*ZCJ(JK)*ZLBDAG(JK)**ICEP%XEX1DEPG) + &
                     (ZZW1(JK,5)+ZZW1(JK,6)) *                            &
                     (ZRHODREF(JK)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-ZZT(JK))))) / &
                               (ZRHODREF(JK)*(CST%XLMTT-CST%XCL*(CST%XTT-ZZT(JK)))))
      END IF
    END DO
!
!*       6.4    Select Wet or Dry case
!
    ZZW(:) = 0.0
    IF (KRR == 7) THEN
      DO JK = 1, KSIZE
        IF (ZRGT(JK) > ICED%XRTMIN(6) .AND. &
            ZZT(JK) < CST%XTT .AND.         &
            ZRDRYG(JK) >= ZRWETG(JK) .AND.  &
            ZRWETG(JK) > 0.0) THEN

          ZZW(JK) = ZRWETG(JK) - ZZW1(JK,5) - ZZW1(JK,6) ! RCWETG+RRWETG
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
          ZZW1(JK,7) = MAX( 0.0,MIN( ZZW(JK),ZRRS(JK)+ZZW1(JK,1) ) )
          ZUSW(JK)   = ZZW1(JK,7) / ZZW(JK)
          ZZW1(JK,5) = ZZW1(JK,5)*ZUSW(JK)
          ZZW1(JK,6) = ZZW1(JK,6)*ZUSW(JK)
          ZRWETG(JK) = ZZW1(JK,7) + ZZW1(JK,5) + ZZW1(JK,6)

          ZRCS(JK) = ZRCS(JK) - ZZW1(JK,1)
          ZRIS(JK) = ZRIS(JK) - ZZW1(JK,5)
          ZRSS(JK) = ZRSS(JK) - ZZW1(JK,6)
!
! assume a linear percent of conversion of graupel into hail
!
          ZRGS(JK) = ZRGS(JK) + ZRWETG(JK)                     !     Wet growth
          ZZW(JK)  = ZRGS(JK)*ZRDRYG(JK)/(ZRWETG(JK)+ZRDRYG(JK)) !        and
          ZRGS(JK) = ZRGS(JK) - ZZW(JK)                        !   partial conversion
          ZRHS(JK) = ZRHS(JK) + ZZW(JK)                        ! of the graupel into hail
!
          ZRRS(JK) = MAX( 0.0,ZRRS(JK) - ZZW1(JK,7) + ZZW1(JK,1) )
          ZTHS(JK) = ZTHS(JK) + ZZW1(JK,7)*(ZLSFACT(JK)-ZLVFACT(JK)) ! f(L_f*(RCWETG+RRWETG)
        END IF
      END DO

    ELSE IF( KRR == 6 ) THEN
      DO JK = 1, KSIZE
        IF (ZRGT(JK) > ICED%XRTMIN(6) .AND. &
            ZRGT(JK) > ICEP%XFRMIN(3) .AND. &
            ZRIS(JK)*PTSTEP > ICEP%XFRMIN(3) .AND. &
            ZZT(JK) < CST%XTT .AND. &
            ZRDRYG(JK) >= ZRWETG(JK) .AND. &
            ZRWETG(JK) > 0.0) THEN
          ZZW(JK)  = ZRWETG(JK)
          ZRCS(JK) = ZRCS(JK) - ZZW1(JK,1)
          ZRIS(JK) = ZRIS(JK) - ZZW1(JK,5)
          ZRSS(JK) = ZRSS(JK) - ZZW1(JK,6)
          ZRGS(JK) = ZRGS(JK) + ZZW(JK)

          ZRRS(JK) = ZRRS(JK) - ZZW(JK) + ZZW1(JK,5) + ZZW1(JK,6) + ZZW1(JK,1)
          ZTHS(JK) = ZTHS(JK) + (ZZW(JK)-ZZW1(JK,5)-ZZW1(JK,6))*(ZLSFACT(JK)-ZLVFACT(JK))
                                                 ! f(L_f*(RCWETG+RRWETG))
        END IF
      END DO
    END IF

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:), &
                                     4,'WETG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     7,'WETG_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     8,'WETG_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     9,'WETG_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    10,'WETG_BU_RRS',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    11,'WETG_BU_RRG',YDDDH, YDLDDH, YDMDDH)
    IF ( KRR == 7 ) THEN
      IF (LBUDGET_RH) CALL BUDGET_DDH(UNPACK(ZRHS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                      12,'WETG_BU_RRH',YDDDH, YDLDDH, YDMDDH)
    END IF

    DO JK = 1, KSIZE
      IF (ZRGT(JK) > ICED%XRTMIN(6) .AND. &
          ZRGT(JK) > ICEP%XFRMIN(4) .AND. &
          ZRIS(JK)*PTSTEP > ICEP%XFRMIN(4) .AND. &
          ZZT(JK) < CST%XTT .AND. &
          ZRDRYG(JK) < ZRWETG(JK) .AND. &
          ZRDRYG(JK) > 0.0) THEN
        ZRCS(JK) = ZRCS(JK) - ZZW1(JK,1)
        ZRIS(JK) = ZRIS(JK) - ZZW1(JK,2)
        ZRSS(JK) = ZRSS(JK) - ZZW1(JK,3)
        ZRRS(JK) = ZRRS(JK) - ZZW1(JK,4)
        ZRGS(JK) = ZRGS(JK) + ZRDRYG(JK)
        ZTHS(JK) = ZTHS(JK) + (ZZW1(JK,1)+ZZW1(JK,4))*(ZLSFACT(JK)-ZLVFACT(JK))
                          ! f(L_f*(RCDRYG+RRDRYG))
      END IF
    END DO

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),   &
                                     4,'DRYG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     7,'DRYG_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     8,'DRYG_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     9,'DRYG_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    10,'DRYG_BU_RRS',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    11,'DRYG_BU_RRG',YDDDH, YDLDDH, YDMDDH)

!*       6.5    Melting of the graupeln

    ZZW(:) = 0.0
    IF (LTIW) THEN

      DO JK = 1, KSIZE
        IF ((ZRGT(JK) > ICED%XRTMIN(6)) .AND. &
            (ZRGS(JK) > 0.0) .AND. &
            (ZTIW(JK) > CST%XTT)) THEN
          ZZW(JK) = ZRVT(JK)*ZPRES(JK)/(CST%XEPSILO+ZRVT(JK)) ! Vapor pressure
          ZZW(JK) = ZKA(JK)*(CST%XTT-ZTIW(JK)) +                                 &
                  (ZDV(JK)*(CST%XLVTT + (CST%XCPV - CST%XCL) * (ZTIW(JK) - CST%XTT)) &
                         *(CST%XESTT-ZZW(JK))/(CST%XRV*ZTIW(JK))             )
!
! compute RGMLTR
!
          ZZW(JK)  = ICEP%XFRMIN(8)*MIN( ZRGS(JK), MAX( 0.0,( -ZZW(JK) *           &
                                 ( ICEP%X0DEPG*       ZLBDAG(JK)**ICEP%XEX0DEPG + &
                                   ICEP%X1DEPG*ZCJ(JK)*ZLBDAG(JK)**ICEP%XEX1DEPG ) - &
                                           ( ZZW1(JK,1)+ZZW1(JK,4) ) *       &
                                    ( ZRHODREF(JK)*CST%XCL*(CST%XTT-ZTIW(JK))) ) /   &
                                                   ( ZRHODREF(JK)*CST%XLMTT)))


          ZRRS(JK) = ZRRS(JK) + ZZW(JK)
          ZRGS(JK) = ZRGS(JK) - ZZW(JK)
          ZTHS(JK) = ZTHS(JK) - ZZW(JK)*(ZLSFACT(JK)-ZLVFACT(JK)) ! f(L_f*(-RGMLTR))
        END IF
      END DO
    ELSE

      DO JK = 1, KSIZE
        IF ((ZRGT(JK)>ICED%XRTMIN(6)) .AND. &
            (ZRGS(JK)>0.0) .AND. &
            (ZZT(JK)>CST%XTT)) THEN
          ZZW(JK) = ZRVT(JK)*ZPRES(JK)/(CST%XEPSILO+ZRVT(JK)) ! Vapor pressure
          ZZW(JK) = ZKA(JK)*(CST%XTT-ZZT(JK)) +                                 &
                  (ZDV(JK)*(CST%XLVTT + (CST%XCPV - CST%XCL) * (ZZT(JK) - CST%XTT)) &
                         *(CST%XESTT-ZZW(JK))/(CST%XRV*ZZT(JK)))
!
! compute RGMLTR
!
          ZZW(JK)  = ICEP%XFRMIN(8)*MIN(ZRGS(JK), MAX(0.0, (-ZZW(JK) *         &
                             (ICEP%X0DEPG*       ZLBDAG(JK)**ICEP%XEX0DEPG +   &
                              ICEP%X1DEPG*ZCJ(JK)*ZLBDAG(JK)**ICEP%XEX1DEPG) - &
                             (ZZW1(JK,1)+ZZW1(JK,4)) *                         &
                                (ZRHODREF(JK)*CST%XCL*(CST%XTT-ZZT(JK)))) /    &
                                               (ZRHODREF(JK)*CST%XLMTT)))
          ZRRS(JK) = ZRRS(JK) + ZZW(JK)
          ZRGS(JK) = ZRGS(JK) - ZZW(JK)
          ZTHS(JK) = ZTHS(JK) - ZZW(JK)*(ZLSFACT(JK)-ZLVFACT(JK)) ! f(L_f*(-RGMLTR))
        END IF
      END DO
    ENDIF

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:), &
                                     4,'GMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     8,'GMLT_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    11,'GMLT_BU_RRG',YDDDH, YDLDDH, YDMDDH)
    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RG',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_FAST_RG

END MODULE MODE_RAIN_ICE_OLD_FAST_RG
