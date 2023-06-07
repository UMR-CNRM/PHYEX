!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_WARM

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_WARM(D, CST, PARAMI, ICEP, ICED, BUCONF,           &
                               KSIZE, OCND2, LKOGAN, GMICRO,                 &
                               PRHODJ, PEVAP3D, PTHS, PRVS,                  &
                               ZRVT, ZRCT, ZRRT, ZRCS, ZRRS, ZTHS,           &
                               ZRVS, ZTHT, ZTHLT,                            &
                               ZCJ, ZKA, ZCF, ZDV, ZRF,                      &
                               ZACRF, ZCONCM,                                &
                               ZRHODREF, ZRHODJ, ZLVFACT, ZLBDAR, ZLBDAR_RF, &
                               ZZKGN_ACON, ZZKGN_SBGR,                       &
                               ZHLC_HCF, ZHLC_LCF, ZHLC_HRC, ZHLC_LRC,       &
                               ZAA2W, ZBB3W,                                 &
                               ZZT, ZPRES, ZESW,                             &
                               TBUDGETS, KBUDGETS)

    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_T
    USE MODD_CST,            ONLY: CST_T
    USE MODD_PARAM_ICE_n,    ONLY: PARAM_ICE_t
    USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_T

    USE MODE_TIWMX,          ONLY: ESATW, AA2W, BB3W

    USE MODE_BUDGET_PHY, ONLY: BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY
    USE MODD_BUDGET,     ONLY: TBUDGETDATA, TBUDGETCONF_t, &
                               NBUDGET_TH, NBUDGET_RR, NBUDGET_RC, NBUDGET_RV

    IMPLICIT NONE

    TYPE(DIMPHYEX_T),       INTENT(IN) :: D
    TYPE(CST_T),            INTENT(IN) :: CST
    TYPE(PARAM_ICE_t),      INTENT(IN) :: PARAMI
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED
    TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF

    INTEGER, INTENT(IN) :: KSIZE

    LOGICAL, INTENT(IN) :: OCND2  ! Logical switch to separate liquid and ice
    LOGICAL, INTENT(IN) :: LKOGAN ! Logical switch for using Kogan autoconversion of liquid.
    LOGICAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PEVAP3D ! Rain evap profile

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PTHS    ! Theta source
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PRVS    ! Water vapor m.r. source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRVT   ! Water vapor m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRCT   ! Cloud water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRRT   ! Rain water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRCS   ! Cloud water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRRS   ! Rain water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZTHS   ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRVS      ! Water vapor m.r. source
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZTHT      ! Potential temperature
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZTHLT     ! Liquid potential temperature

    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZRHODREF  ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZRHODJ    ! RHO times Jacobian
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZLVFACT   ! L_v/(Pi_ref*C_ph)
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZLBDAR    ! Slope parameter of the raindrop  distribution
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZLBDAR_RF ! Slope parameter of the raindrop  distribution

    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZCJ       ! Function to compute the ventilation coefficient
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZKA       ! Thermal conductivity of the air
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZCF       ! Cloud fraction
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZDV       ! Diffusivity of water vapor in the air
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZRF       ! Rain fraction

    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZACRF     ! collision factor cloud liquid to rain
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZCONCM    ! Same as ZCONC3D but GMICRO-array only and cm-3 instead of m-3

    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZHLC_HCF  ! HLCLOUDS : fraction of High Cloud Fraction in grid
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZHLC_LCF  ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZHLC_HRC  ! HLCLOUDS : LWC that is High LWC in grid
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZHLC_LRC  ! HLCLOUDS : LWC that is Low  LWC in grid

    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZAA2W     ! as ZAA2 but for liquid
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZBB3W     ! as ZBB3 but for liquid

    REAL, DIMENSION(KSIZE), INTENT(IN)  :: ZZT       ! Temperature
    REAL, DIMENSION(KSIZE), INTENT(IN)  :: ZPRES     ! Pressure

    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZESW      ! saturation pressure over water

!   SPP arrays
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZZKGN_ACON
    REAL, DIMENSION(KSIZE), INTENT(IN) :: ZZKGN_SBGR

    TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
    INTEGER, INTENT(IN) :: KBUDGETS

    REAL, DIMENSION(KSIZE) :: ZUSW     ! Undersaturation over water

    REAL, DIMENSION(KSIZE) :: ZZW      ! Work array
    REAL, DIMENSION(KSIZE) :: ZZW2     ! Work array
    REAL, DIMENSION(KSIZE) :: ZZW3     ! Work array
    REAL, DIMENSION(KSIZE) :: ZZW4     ! Work array
    REAL, DIMENSION(KSIZE) :: ZARTMP   ! temporary work array

    REAL, DIMENSION(D%NIT,D%NKT) :: ZW ! work array

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    INTEGER :: JK
    LOGICAL :: LTEST ! Only for test !
!
!-------------------------------------------------------------------------------
!
!*       4.2    compute the autoconversion of r_c for r_r production: RCAUTR
!
    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_WARM',0,ZHOOK_HANDLE)

    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'AUTO', UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'AUTO', UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    IF (LKOGAN) THEN
      DO JK = 1, KSIZE
        IF (ZRCT(JK) >  1.0E-8) THEN ! Closely following Kogan autoconversion
          ZZW(JK) = 1350.0*ZZKGN_ACON(JK)* ZCONCM(JK)**(-1.79) * &
                &  (ZRCT(JK)/(MAX(ZZKGN_SBGR(JK),ZCF(JK))))**2.47
          ZZW(JK) = ZZW(JK)*MAX(ZZKGN_SBGR(JK),ZCF(JK))
          ZZW(JK) = MIN( ZRCS(JK),ZZW(JK))
          ZRCS(JK) = ZRCS(JK) - ZZW(JK)
          ZRRS(JK) = ZRRS(JK) + ZZW(JK)
        END IF
      END DO
    ELSE
      DO JK = 1, KSIZE
        IF (ZRCS(JK) > 0.0 .AND. ZHLC_HCF(JK) > 0.0) THEN
          ZZW(JK) = ICEP%XTIMAUTC*MAX(ZHLC_HRC(JK)/ZHLC_HCF(JK) - ICEP%XCRIAUTC/ZRHODREF(JK),0.0)
          ZZW(JK) = MIN(ZRCS(JK),ZHLC_HCF(JK)*ZZW(JK))
          ZRCS(JK) = ZRCS(JK) - ZZW(JK)
          ZRRS(JK) = ZRRS(JK) + ZZW(JK)
        END IF
      END DO
    ENDIF

    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'AUTO', UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'AUTO', UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

!*       4.3    compute the accretion of r_c for r_r production: RCACCR

    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'ACCR', UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'ACCR', UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    IF (PARAMI%CSUBG_RC_RR_ACCR == 'NONE') THEN
      !CLoud water and rain are diluted over the grid box
      DO JK = 1, KSIZE
        IF (ZRCT(JK)>ICED%XRTMIN(2) .AND. ZRRT(JK)>ICED%XRTMIN(3) .AND. ZRCS(JK)>0.0) THEN
          ZZW(JK) = MIN( ZRCS(JK), ICEP%XFCACCR * ZRCT(JK)*ZACRF(JK) &
                 * ZLBDAR(JK)**ICEP%XEXCACCR    &
                 * ZRHODREF(JK)**(-ICED%XCEXVT) )
          ZRCS(JK) = ZRCS(JK) - ZZW(JK)
          ZRRS(JK) = ZRRS(JK) + ZZW(JK)
        END IF
      END DO

    ELSEIF (PARAMI%CSUBG_RC_RR_ACCR=='PRFR') THEN
      !Cloud water is concentrated over its fraction with possibly to parts with high and low content as set for autoconversion
      !Rain is concnetrated over its fraction
      !Rain in high content area fraction: ZHLC_HCF
      !Rain in low content area fraction:
      ! if ZRF<ZCF (rain is entirely falling in cloud): ZRF-ZHLC_HCF
      ! if ZRF>ZCF (rain is falling in cloud and in clear sky): ZCF-ZHLC_HCF
      ! => min(ZCF, ZRF)-ZHLC_HCF
      ZZW(:) = 0.
      DO JK = 1, KSIZE
        IF (ZHLC_HRC(JK) > ICED%XRTMIN(2) .AND. &
            ZRRT(JK) > ICED%XRTMIN(3) .AND.     &
            ZRCS(JK) > 0.0 .AND.                &
            ZHLC_HCF(JK) > 0) THEN

            !Accretion due to rain falling in high cloud content
            ZZW(JK) = ICEP%XFCACCR * (ZHLC_HRC(JK)/ZHLC_HCF(JK)) &
                    * ZLBDAR_RF(JK)**ICEP%XEXCACCR &
                    * ZRHODREF(JK)**(-ICED%XCEXVT) &
                    * ZHLC_HCF(JK)
        END IF
      END DO

      DO JK = 1, KSIZE
        IF (ZHLC_LRC(JK) > ICED%XRTMIN(2) .AND. &
            ZRRT(JK) > ICED%XRTMIN(3) .AND.     &
            ZRCS(JK) > 0.0 .AND.                &
            ZHLC_LCF(JK) > 0) THEN

          !We add acrretion due to rain falling in low cloud content
          ZZW(JK) = ZZW(JK) + ICEP%XFCACCR * ( ZHLC_LRC(JK)/ZHLC_LCF(JK) ) &
                            * ZLBDAR_RF(JK)**ICEP%XEXCACCR &
                            * ZRHODREF(JK)**(-ICED%XCEXVT) &
                            * (MIN(ZCF(JK), ZRF(JK))-ZHLC_HCF(JK))
        END IF
      END DO
      ZZW(:)=MIN(ZRCS(:), ZZW(:))
      ZRCS(:) = ZRCS(:) - ZZW(:)
      ZRRS(:) = ZRRS(:) + ZZW(:)

    ELSE
      !wrong CSUBG_RC_RR_ACCR case
      CALL ABORT
      STOP 'wrong CSUBG_RC_RR_ACCR case'
    ENDIF

    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'ACCR', UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'ACCR', UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

!*       4.4    compute the evaporation of r_r: RREVAV

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'REVA', UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'REVA', UNPACK(ZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'REVA', UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    ZZW(:) = 0.0

    IF (PARAMI%CSUBG_RR_EVAP=='NONE') THEN

      !Evaporation only when there's no cloud (RC must be 0)
      IF (OCND2) THEN
        DO JK = 1, KSIZE
          IF ((ZRRT(JK)>ICED%XRTMIN(3)) .AND. (ZRCT(JK)<=ICED%XRTMIN(2))) THEN
             ZZW(JK) = ZAA2W(JK) + ZBB3W(JK)*ZPRES(JK)
             ZUSW(JK) = 1.0 - ZRVT(JK)*( ZPRES(JK)-ZESW(JK) ) / ( CST%XEPSILO * ZESW(JK) )
                                                    ! Undersaturation over water
             ZZW(JK) = MIN(ZRRS(JK), (MAX(0.0, ZUSW(JK))/(ZRHODREF(JK)*ZZW(JK))) *      &
                  (ICEP%X0EVAR*ZLBDAR(JK)**ICEP%XEX0EVAR + ICEP%X1EVAR*ZCJ(JK)*ZLBDAR(JK)**ICEP%XEX1EVAR))
             ZRRS(JK) = ZRRS(JK) - ZZW(JK)
             ZRVS(JK) = ZRVS(JK) + ZZW(JK)
             ZTHS(JK) = ZTHS(JK) - ZZW(JK)*ZLVFACT(JK)
          END IF
        END DO
      ELSE
        DO JK = 1, KSIZE
          IF ((ZRRT(JK) > ICED%XRTMIN(3)) .AND. (ZRCT(JK) <= ICED%XRTMIN(2))) THEN

            ZZW(JK)  = EXP( CST%XALPW - CST%XBETAW/ZZT(JK) - CST%XGAMW*ALOG(ZZT(JK) ) ) ! es_w
            ZUSW(JK) = 1.0 - ZRVT(JK)*( ZPRES(JK)-ZZW(JK) ) / ( CST%XEPSILO * ZZW(JK) )
                                                        ! Undersaturation over water
            ZZW(JK) = (CST%XLVTT+(CST%XCPV-CST%XCL)*(ZZT(JK)-CST%XTT))**2 / (ZKA(JK)*CST%XRV*ZZT(JK)**2) &
                   + (CST%XRV*ZZT(JK) ) / ( ZDV(JK)*ZZW(JK))
            ZZW(JK) = MIN(ZRRS(JK), (MAX(0.0, ZUSW(JK))/(ZRHODREF(JK)*ZZW(JK))) *      &
                        (ICEP%X0EVAR*ZLBDAR(JK)**ICEP%XEX0EVAR + ICEP%X1EVAR*ZCJ(JK)*ZLBDAR(JK)**ICEP%XEX1EVAR))
            ZRRS(JK) = ZRRS(JK) - ZZW(JK)
            ZRVS(JK) = ZRVS(JK) + ZZW(JK)
            ZTHS(JK) = ZTHS(JK) - ZZW(JK)*ZLVFACT(JK)
          END IF
        END DO
      ENDIF

    ELSEIF (PARAMI%CSUBG_RR_EVAP=='CLFR' .OR. PARAMI%CSUBG_RR_EVAP=='PRFR') THEN
      !Evaporation in clear sky part
      !With CLFR, rain is diluted over the grid box
      !With PRFR, rain is concentrated in its fraction
      !Use temperature and humidity in clear sky part like Bechtold et al. (1993)
      IF (PARAMI%CSUBG_RR_EVAP=='CLFR') THEN
        ZZW4(:)=1. !Precipitation fraction
        ZZW3(:)=ZLBDAR(:)
      ELSE
        ZZW4(:)=ZRF(:) !Precipitation fraction
        ZZW3(:)=ZLBDAR_RF(:)
      ENDIF

      !ATTENTION
      !Il faudrait recalculer les variables ZKA, ZDV, ZCJ en tenant compte de la température T^u
      !Ces variables devraient être sorties de rain_ice_slow et on mettrait le calcul de T^u, T^s
      !et plusieurs versions (comme actuellement, en ciel clair, en ciel nuageux) de ZKA, ZDV, ZCJ dans rain_ice
      !On utiliserait la bonne version suivant l'option NONE, CLFR... dans l'évaporation et ailleurs

      LTEST = .FALSE.
      IF(OCND2.AND.LTEST) THEN
        DO JK= 1,KSIZE
          IF ((ZRRT(JK) > ICED%XRTMIN(3)) .AND. (ZZW4(JK) > ZCF(JK))) THEN
            ! outside the cloud (environment) the use of T^u (unsaturated) instead of T
            ! Bechtold et al. 1993
            !
            ! T^u = T_l = theta_l * (T/theta)
            ZZW2(JK) =  ZTHLT(JK) * ZZT(JK) / ZTHT(JK) ! ZZW2 = Temperature
            ZZW(JK) = AA2W(ZZW2(JK)) + BB3W(ZZW2(JK))*ZPRES(JK) ! ZZW = Droplet function
            ZARTMP(JK)= ESATW(ZZW2(JK)) ! saturation pressure, water
          ENDIF
        ENDDO

        DO JK= 1,KSIZE
          IF ((ZRRT(JK) > ICED%XRTMIN(3)) .AND. (ZZW4(JK) > ZCF(JK))) THEN

            ! S, Undersaturation over water (with new theta^u)
            ZUSW(JK) = 1.0 - ZRVT(JK)*(ZPRES(JK) - ZARTMP(JK)) / (CST%XEPSILO * ZARTMP(JK))
            ! New ZCJ(JK) for  T^u
            ZARTMP(JK) = ICEP%XSCFAC * ZRHODREF(JK)**0.3 / SQRT(1.718E-5+0.0049E-5*(ZZW2(JK)-CST%XTT))
            ZZW(JK) = MAX( 0.0,ZUSW(JK) )/(ZRHODREF(JK)*ZZW(JK))  *      &
                    (ICEP%X0EVAR*ZZW3(JK)**ICEP%XEX0EVAR + ICEP%X1EVAR*ZARTMP(JK)*ZZW3(JK)**ICEP%XEX1EVAR)
            ZZW(JK) = MIN( ZRRS(JK),  ZZW(JK) *( ZZW4(JK) - ZCF(JK) ) )
            ZRRS(JK) = ZRRS(JK) - ZZW(JK)
            ZRVS(JK) = ZRVS(JK) + ZZW(JK)
            ZTHS(JK) = ZTHS(JK) - ZZW(JK)*ZLVFACT(JK)
          ENDIF
        ENDDO

      ELSE
        DO JK = 1, KSIZE
          IF ((ZRRT(JK)>ICED%XRTMIN(3)) .AND. (ZZW4(JK) > ZCF(JK))) THEN
            ! outside the cloud (environment) the use of T^u (unsaturated) instead of T
            ! Bechtold et al. 1993

            ! T^u = T_l = theta_l * (T/theta)
            ZZW2(JK) =  ZTHLT(JK) * ZZT(JK) / ZTHT(JK)

            ! es_w with new T^u
            ZZW(JK)  = EXP(CST%XALPW - CST%XBETAW/ZZW2(JK) - CST%XGAMW*ALOG(ZZW2(JK)))

            ! S, Undersaturation over water (with new theta^u)
            ZUSW(JK) = 1.0 - ZRVT(JK)*(ZPRES(JK) - ZZW(JK)) / (CST%XEPSILO * ZZW(JK))

            ZZW(JK) = (CST%XLVTT+(CST%XCPV-CST%XCL)*(ZZW2(JK)-CST%XTT))**2 / (ZKA(JK)*CST%XRV*ZZW2(JK)**2) &
                   + (CST%XRV*ZZW2(JK)) / (ZDV(JK)*ZZW(JK))

            ZZW(JK) = MAX(0.0, ZUSW(JK))/(ZRHODREF(JK)*ZZW(JK)) * &
                   (ICEP%X0EVAR*ZZW3(JK)**ICEP%XEX0EVAR + ICEP%X1EVAR*ZCJ(JK)*ZZW3(JK)**ICEP%XEX1EVAR)

            ZZW(JK) = MIN( ZRRS(JK),  ZZW(JK) *( ZZW4(JK) - ZCF(JK) ) )

            ZRRS(JK) = ZRRS(JK) - ZZW(JK)
            ZRVS(JK) = ZRVS(JK) + ZZW(JK)
            ZTHS(JK) = ZTHS(JK) - ZZW(JK)*ZLVFACT(JK)
          END IF
        END DO
      ENDIF
    ELSE
      !wrong CSUBG_RR_EVAP case
      CALL ABORT
      STOP 'wrong CSUBG_RR_EVAP case'
    END IF

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'REVA', UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'REVA', UNPACK(ZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'REVA', UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    ZW(:,:) = PEVAP3D(:,:)
    PEVAP3D(:,:) = UNPACK(ZZW(:), MASK=GMICRO(:,:), FIELD=ZW(:,:))
!
  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_WARM',1,ZHOOK_HANDLE)
  END SUBROUTINE RAIN_ICE_OLD_WARM

END MODULE MODE_RAIN_ICE_OLD_WARM
