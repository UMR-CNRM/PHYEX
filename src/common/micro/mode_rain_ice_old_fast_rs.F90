!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_FAST_RS

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_FAST_RS(D, CST, ICEP, ICED, BUCONF,         &
                                  PTSTEP, KSIZE, KRR, GMICRO,         &
                                  PRHODJ, PTHS,                       &
                                  ZRVT, ZRCT, ZRRT, ZRST,             &
                                  ZRRS, ZRCS, ZRSS, ZRGS, ZTHS,       &
                                  ZRHODREF, ZRHODJ, ZLSFACT, ZLVFACT, &
                                  ZCJ, ZKA, ZDV,                      &
                                  ZLBDAR, ZLBDAS, ZCOLF, ZPRES, ZZT,  &
                                  TBUDGETS, KBUDGETS)

    USE PARKIND1,            ONLY: JPRB
    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
    USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_T
    USE MODD_CST,            ONLY: CST_T
    USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_T

    USE MODE_BUDGET_PHY, ONLY: BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY
    USE MODD_BUDGET,     ONLY: TBUDGETDATA, TBUDGETCONF_t, &
                               LBU_ENABLE, NBUDGET_TH, NBUDGET_RG, NBUDGET_RR, NBUDGET_RC, NBUDGET_RS

    IMPLICIT NONE

    TYPE(DIMPHYEX_T),       INTENT(IN) :: D
    TYPE(CST_T),            INTENT(IN) :: CST 
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED
    TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF

    REAL,    INTENT(IN) :: PTSTEP  ! Double Time step
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, INTENT(IN) :: KRR

    LOGICAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PTHS    ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRVT   ! Water vapor m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRCT   ! Cloud water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRRT   ! Rain water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRST   ! Snow/aggregate m.r. at t

    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRRS   ! Rain water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRCS   ! Cloud water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRSS   ! Snow/aggregate m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRGS   ! Graupel m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZTHS   ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRHODREF ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRHODJ   ! RHO times Jacobian
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLSFACT  ! L_s/(Pi_ref*C_ph)
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLVFACT  ! L_v/(Pi_ref*C_ph)

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZCJ      ! Function to compute the ventilation coefficient
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZKA      ! Thermal conductivity of the air
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZDV      ! Diffusivity of water vapor in the air

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLBDAR   ! Slope parameter of the raindrop  distribution
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLBDAS   ! Slope parameter of the aggregate distribution

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZCOLF    ! collision factor cloud liquid to snow / graupel
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZZT      ! Temperature
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZPRES    ! Pressure

    TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
    INTEGER, INTENT(IN) :: KBUDGETS

    LOGICAL, DIMENSION(KSIZE) :: GMASK ! Test where to compute riming/accretion

    INTEGER, DIMENSION(KSIZE) :: IVEC1 ! Vectors of indices for
    INTEGER, DIMENSION(KSIZE) :: IVEC2 ! Vectors of indices for

    REAL, DIMENSION(KSIZE) :: ZVEC1 ! Work vectors for interpolations
    REAL, DIMENSION(KSIZE) :: ZVEC2 ! Work vectors for interpolations
    REAL, DIMENSION(KSIZE) :: ZVEC3 ! Work vectors for interpolations

    REAL, DIMENSION(KSIZE)      :: ZZW      ! Work array
    REAL, DIMENSION(KSIZE, KRR) :: ZZW1     ! Work array

    INTEGER :: IGRIM, IGACC
    INTEGER :: JJ, JK

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       5.1    cloud droplet riming of the aggregates
!
    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RS',0,ZHOOK_HANDLE)

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'RIM', UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'RIM', UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'RIM', UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'RIM', UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    ZZW1(:,:) = 0.0
!
    GMASK(:) = (ZRCT(:)>ICED%XRTMIN(2)) .AND. (ZRST(:)>ICED%XRTMIN(5)) .AND. &
                                  (ZRCS(:)>0.0) .AND. (ZZT(:)<CST%XTT)
    IGRIM = COUNT( GMASK(:) )
!
    IF( IGRIM>0 ) THEN
!
!        5.1.1  select the ZLBDAS
!
      ZVEC1(:) = PACK( ZLBDAS(:),MASK=GMASK(:) )
!
!        5.1.2  find the next lower indice for the ZLBDAS in the geometrical
!               set of Lbda_s used to tabulate some moments of the incomplete
!               gamma function
!
      ZVEC2(1:IGRIM) = MAX(1.00001, MIN(FLOAT(ICEP%NGAMINC) - 0.00001,           &
                           ICEP%XRIMINTP1 * LOG(ZVEC1(1:IGRIM)) + ICEP%XRIMINTP2))
      IVEC2(1:IGRIM) = INT(ZVEC2(1:IGRIM))
      ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - FLOAT(IVEC2(1:IGRIM))
!
!        5.1.3  perform the linear interpolation of the normalized
!               "2+XDS"-moment of the incomplete gamma function
!
      ZVEC1(1:IGRIM) = ICEP%XGAMINC_RIM1(IVEC2(1:IGRIM)+1)* ZVEC2(1:IGRIM)      &
                   - ICEP%XGAMINC_RIM1(IVEC2(1:IGRIM)  )*(ZVEC2(1:IGRIM) - 1.0)
      ZZW(:) = UNPACK(VECTOR=ZVEC1(:), MASK=GMASK, FIELD=0.0)
!
!        5.1.4  riming of the small sized aggregates
!
      DO JK = 1, KSIZE
        IF (GMASK(JK)) THEN
          ZZW1(JK,1) = MIN( ZRCS(JK),                                 &
                         ICEP%XCRIMSS * ZZW(JK) * ZRCT(JK)*ZCOLF(JK)   & ! RCRIMSS
                                      *   ZLBDAS(JK)**ICEP%XEXCRIMSS &
                                      * ZRHODREF(JK)**(-ICED%XCEXVT) )
          ZRCS(JK) = ZRCS(JK) - ZZW1(JK,1)
          ZRSS(JK) = ZRSS(JK) + ZZW1(JK,1)
          ZTHS(JK) = ZTHS(JK) + ZZW1(JK,1)*(ZLSFACT(JK)-ZLVFACT(JK)) ! f(L_f*(RCRIMSS))
        END IF
      END DO
!
!        5.1.5  perform the linear interpolation of the normalized
!               "XBS"-moment of the incomplete gamma function
!
      ZVEC1(1:IGRIM) = ICEP%XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                     - ICEP%XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
      ZZW(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GMASK,FIELD=0.0 )
!
!        5.1.6  riming-conversion of the large sized aggregates into graupeln
!
!
      DO JK = 1, KSIZE
        IF (GMASK(JK) .AND. (ZRSS(JK) > 0.0)) THEN
          ZZW1(JK,2) = MIN(ZRCS(JK),                                   &
                          ICEP%XCRIMSG * ZRCT(JK)*ZCOLF(JK)            & ! RCRIMSG
                                       * ZLBDAS(JK)**ICEP%XEXCRIMSG   &
                                       * ZRHODREF(JK)**(-ICED%XCEXVT) &
                                       - ZZW1(JK,1))

          ZZW1(JK,3) = MIN(ZRSS(JK),                                 &
                          ICEP%XSRIMCG * ZLBDAS(JK)**ICEP%XEXSRIMCG & ! RSRIMCG
                                       * (1.0 - ZZW(JK))/(PTSTEP*ZRHODREF(JK)))

          ZRCS(JK) = ZRCS(JK) - ZZW1(JK,2)
          ZRSS(JK) = ZRSS(JK) - ZZW1(JK,3)
          ZRGS(JK) = ZRGS(JK) + ZZW1(JK,2)+ZZW1(JK,3)
          ZTHS(JK) = ZTHS(JK) + ZZW1(JK,2)*(ZLSFACT(JK)-ZLVFACT(JK)) ! f(L_f*(RCRIMSG))
        END IF
      END DO
    END IF

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'RIM', UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'RIM', UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'RIM', UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'RIM', UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

!*       5.2    rain accretion onto the aggregates

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'ACC', UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'ACC', UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'ACC', UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'ACC', UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    ZZW1(:,2:3) = 0.0
    GMASK(:) = (ZRRT(:) > ICED%XRTMIN(3)) .AND. &
              (ZRST(:) > ICED%XRTMIN(5)) .AND. &
              (ZRRS(:) > 0.0)            .AND. &
              (ZZT(:) < CST%XTT)

    IGACC = COUNT(GMASK(:))
!
    IF( IGACC>0 ) THEN
!
!        5.2.1  select the (ZLBDAS,ZLBDAR) couplet
!
      ZVEC1(:) = PACK( ZLBDAS(:),MASK=GMASK(:) )
      ZVEC2(:) = PACK( ZLBDAR(:),MASK=GMASK(:) )
!
!        5.2.2  find the next lower indice for the ZLBDAS and for the ZLBDAR
!               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
!               tabulate the RACCSS-kernel
!
      ZVEC1(1:IGACC) = MAX( 1.00001, MIN( FLOAT(ICEP%NACCLBDAS)-0.00001,           &
                            ICEP%XACCINTP1S * LOG( ZVEC1(1:IGACC) ) + ICEP%XACCINTP2S ) )
      IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
      ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - FLOAT( IVEC1(1:IGACC) )

      ZVEC2(1:IGACC) = MAX( 1.00001, MIN( FLOAT(ICEP%NACCLBDAR)-0.00001,           &
                            ICEP%XACCINTP1R * LOG( ZVEC2(1:IGACC) ) + ICEP%XACCINTP2R ) )
      IVEC2(1:IGACC) = INT( ZVEC2(1:IGACC) )
      ZVEC2(1:IGACC) = ZVEC2(1:IGACC) - FLOAT( IVEC2(1:IGACC) )

!        5.2.3  perform the bilinear interpolation of the normalized
!               RACCSS-kernel

      DO JJ = 1,IGACC
        ZVEC3(JJ) = ( ICEP%XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)         &
                    - ICEP%XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0)) &
                                                               * ZVEC1(JJ)         &
                  - ( ICEP%XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)         &
                    - ICEP%XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0)) &
                                                               *(ZVEC1(JJ) - 1.0)
      END DO
      ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GMASK,FIELD=0.0 )
!
!        5.2.4  raindrop accretion on the small sized aggregates
!
      DO JK = 1, KSIZE
        IF (GMASK(JK)) THEN
          ZZW1(JK,2) =                                            & !! coef of RRACCS
                  ICEP%XFRACCSS*( ZLBDAS(JK)**ICED%XCXS )*( ZRHODREF(JK)**(-ICED%XCEXVT-1.) ) &
           *(ICEP%XLBRACCS1/((ZLBDAS(JK)**2)               ) +                  &
             ICEP%XLBRACCS2/( ZLBDAS(JK)   * ZLBDAR(JK)    ) +                  &
             ICEP%XLBRACCS3/(               (ZLBDAR(JK)**2)) )/ZLBDAR(JK)**4
          ZZW1(JK,4) = MIN( ZRRS(JK),ZZW1(JK,2)*ZZW(JK) )           ! RRACCSS
          ZRRS(JK) = ZRRS(JK) - ZZW1(JK,4)*ICEP%XFRMIN(7)
          ZRSS(JK) = ZRSS(JK) + ZZW1(JK,4)*ICEP%XFRMIN(7)
          ZTHS(JK) = ZTHS(JK) + ZZW1(JK,4)*(ZLSFACT(JK)-ZLVFACT(JK))*ICEP%XFRMIN(7) ! f(L_f*(RRACCSS))
        END IF
      END DO
!
!        5.2.4b perform the bilinear interpolation of the normalized
!               RACCS-kernel
!
      DO JJ = 1,IGACC
        ZVEC3(JJ) =  (ICEP%XKER_RACCS(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    - ICEP%XKER_RACCS(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                              * ZVEC2(JJ)          &
                   - (ICEP%XKER_RACCS(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                   -  ICEP%XKER_RACCS(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                              *(ZVEC2(JJ) - 1.0)
      END DO
      ZZW1(:,2) = ZZW1(:,2)*UNPACK( VECTOR=ZVEC3(:),MASK=GMASK(:),FIELD=0.0 )
                                                                       !! RRACCS!
!        5.2.5  perform the bilinear interpolation of the normalized
!               SACCRG-kernel
!
      DO JJ = 1,IGACC
        ZVEC3(JJ) =  (  ICEP%XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                      - ICEP%XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                            * ZVEC2(JJ) &
                   - (  ICEP%XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                      - ICEP%XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                          * (ZVEC2(JJ) - 1.0)
      END DO
      ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GMASK,FIELD=0.0 )
!
!        5.2.6  raindrop accretion-conversion of the large sized aggregates
!               into graupeln
!
      DO JK = 1, KSIZE
        IF (GMASK(JK) .AND. (ZRSS(JK) > 0.0)) THEN
          ZZW1(JK,2) = MAX( MIN( ZRRS(JK),ZZW1(JK,2)-ZZW1(JK,4) ),0.0 )       ! RRACCSG
        END IF
      END DO

      DO JK = 1, KSIZE
        IF (GMASK(JK) .AND. (ZRSS(JK)>0.0) .AND. ZZW1(JK,2) > 0.0 .AND. ZRSS(JK) > ICEP%XFRMIN(1)/PTSTEP) THEN
          ZZW1(JK,3) = MIN( ZRSS(JK),ICEP%XFSACCRG*ZZW(JK)*                     & ! RSACCRG
                ( ZLBDAS(JK)**(ICED%XCXS-ICED%XBS) )*( ZRHODREF(JK)**(-ICED%XCEXVT-1.) ) &
               *( ICEP%XLBSACCR1/((ZLBDAR(JK)**2)               ) +           &
                  ICEP%XLBSACCR2/( ZLBDAR(JK)    * ZLBDAS(JK)    ) +           &
                  ICEP%XLBSACCR3/(               (ZLBDAS(JK)**2)) )/ZLBDAR(JK) )
          ZRRS(JK) = ZRRS(JK) - ZZW1(JK,2)
          ZRSS(JK) = ZRSS(JK) - ZZW1(JK,3)
          ZRGS(JK) = ZRGS(JK) + ZZW1(JK,2)+ZZW1(JK,3)
          ZTHS(JK) = ZTHS(JK) + ZZW1(JK,2)*(ZLSFACT(JK)-ZLVFACT(JK)) !
                                 ! f(L_f*(RRACCSG))
        END IF
      END DO
    END IF

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'ACC', UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'ACC', UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'ACC', UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'ACC', UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

!*       5.3    Conversion-Melting of the aggregates

    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'CMEL', UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'CMEL', UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    ZZW(:) = 0.0
    DO JK = 1, KSIZE
      IF ((ZRST(JK)>ICED%XRTMIN(5)) .AND. (ZRSS(JK)>0.0) .AND. (ZZT(JK)>CST%XTT)) THEN
        ZZW(JK) = ZRVT(JK)*ZPRES(JK)/(CST%XEPSILO+ZRVT(JK)) ! Vapor pressure
        ZZW(JK) =  ZKA(JK)*(CST%XTT-ZZT(JK)) +                                 &
                 ( ZDV(JK)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZZT(JK) - CST%XTT )) &
                             *(CST%XESTT-ZZW(JK))/(CST%XRV*ZZT(JK))             )
!
! compute RSMLT
!
        ZZW(JK)  = MIN( ZRSS(JK), ICEP%XFSCVMG*MAX( 0.0,( -ZZW(JK) *             &
                             ( ICEP%X0DEPS*       ZLBDAS(JK)**ICEP%XEX0DEPS + &
                               ICEP%X1DEPS*ZCJ(JK)*ZLBDAS(JK)**ICEP%XEX1DEPS ) -   &
                                       ( ZZW1(JK,1)+ZZW1(JK,4) ) *       &
                                (ZRHODREF(JK)*CST%XCL*(CST%XTT-ZZT(JK)))) /    &
                                               ( ZRHODREF(JK)*CST%XLMTT ) ) )
!
! note that RSCVMG = RSMLT*ICEP%XFSCVMG but no heat is exchanged (at the rate RSMLT)
! because the graupeln produced by this process are still icy!!!
!
        ZRSS(JK) = ZRSS(JK) - ZZW(JK)
        ZRGS(JK) = ZRGS(JK) + ZZW(JK)
      END IF
    END DO

    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'CMEL', UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'CMEL', UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RS',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_FAST_RS

END MODULE MODE_RAIN_ICE_OLD_FAST_RS
