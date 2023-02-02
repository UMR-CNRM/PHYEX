!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_FAST_RI

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_FAST_RI(D, CST, ICEP, ICED,           &
                                  PTSTEP, KSIZE,                &
                                  OCND2, LMODICEDEP, GMICRO,    &
                                  PRHODJ, PTHS,                 &
                                  ZRIT, ZCIT,                   &
                                  ZRVS, ZRCS, ZRIS, ZRSS, ZTHS, &
                                  ZRHODREF, ZRHODJ,             &
                                  ZLSFACT, ZLVFACT,             &
                                  ZAI, ZCJ,                     &
                                  ZSSIO, ZSSIU, ZW2D, ZXW2D13,  &
                                  ZZT, ZPRES, ZSSI,             &
                                  ZSIFRC, ZESI,                 &
                                  ZCITRED, ZCITRED23, ZDICRIT,  &
                                  YDDDH, YDLDDH, YDMDDH)

    USE PARKIND1,            ONLY: JPRB
    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
    USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_T
    USE MODD_CST,            ONLY: CST_T
    USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_T

    USE MODE_RAIN_ICE_OLD_ICENUMBER2, ONLY: ICENUMBER2

    USE MODE_BUDGET,         ONLY: BUDGET_DDH
    USE MODD_BUDGET,         ONLY: LBUDGET_TH, LBUDGET_RC, LBUDGET_RI
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

    LOGICAL, INTENT(IN) :: OCND2
    LOGICAL, INTENT(IN) :: LMODICEDEP

    LOGICAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PTHS    ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRIT     ! Pristine ice m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZCIT     ! Pristine ice conc. at t

    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRVS    ! Water vapor m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRCS     ! Cloud water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRIS     ! Pristine ice m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRSS     ! Snow/aggregate m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZTHS     ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRHODREF ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRHODJ   ! RHO times Jacobian
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLSFACT  ! L_s/(Pi_ref*C_ph)
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLVFACT  ! L_v/(Pi_ref*C_ph)

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZAI      ! Thermodynamical function
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZCJ      ! Function to compute the ventilation coefficient

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZSSIO    ! Super-saturation with respect to ice in the
                                                      ! supersaturated fraction of gridbox
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZSSIU    ! Sub-saturation with respect to ice in the
                                                      ! sub-saturated fraction of gridbox
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZW2D     ! Factor for subgridscale calculations
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZXW2D13  ! ZXW2D**0.333 or other expression for LMODICEDEP=T

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZZT      ! Temperature
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZPRES    ! Pressure
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZSSI     ! Supersaturation over ice

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZSIFRC   ! subgridscale fraction with supersaturation with
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZESI     ! saturation pressure over ice

    REAL, INTENT(IN) :: ZCITRED
    REAL, INTENT(IN) :: ZCITRED23
    REAL, INTENT(IN) :: ZDICRIT

    TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
    TYPE(TLDDH),   INTENT(IN)    :: YDLDDH
    TYPE(TMDDH),   INTENT(IN)    :: YDMDDH

    REAL, DIMENSION(KSIZE) :: ZCRYSHA ! Ice crystal habit factor
    REAL, DIMENSION(KSIZE) :: ZCI2S   ! factor to turn cloud ice with few lagre crystals into snow
    REAL, DIMENSION(KSIZE) :: ZZW     ! Work array
    REAL, DIMENSION(KSIZE) :: ZZWC    ! Work array
    REAL, DIMENSION(KSIZE) :: ZZW2    ! Work array

    REAL :: ZTC
    REAL :: ZQIMAX
    REAL :: ZHU

    INTEGER :: JJ, JK, JL

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       7.1    cloud ice melting
!

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RI',0,ZHOOK_HANDLE)

    ZZW(:) = 0.0
    DO JK = 1, KSIZE
      IF ((ZRIS(JK)>0.0) .AND. (ZZT(JK)>CST%XTT)) THEN
        ZZW(JK)  = ZRIS(JK)
        ZRCS(JK) = ZRCS(JK) + ZRIS(JK)
        ZTHS(JK) = ZTHS(JK) - ZRIS(JK)*(ZLSFACT(JK)-ZLVFACT(JK)) ! f(L_f*(-RIMLTC))
        ZRIS(JK) = 0.0
        ZCIT(JK) = 0.0
      END IF
    END DO

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),   &
                                     4,'IMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     7,'IMLT_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     9,'IMLT_BU_RRI',YDDDH, YDLDDH, YDMDDH)

!*       7.2    Bergeron-Findeisen effect: RCBERI

    ZZW(:) = 0.0
    IF(OCND2)THEN

      ! Sub gridscale decomposition into a supersaturation part of the gridbox,
      ! ZSIFRC with a superaturation ZSSIO and a subsaturated part (1.- ZSIFRC)
      ! with a (negative) superaturation of ZSSIU

      IF (LMODICEDEP) THEN

        DO JL=1,KSIZE
          ZZW2(JL) = MAX(ZCIT(JL),ICENUMBER2(ZRIS(JL)*PTSTEP,ZZT(JL))* &
          ZRHODREF(JL))
        ENDDO

        DO JK = 1, KSIZE
          IF (ZZW2(JK)>0.0 .AND. ZESI(JK) < ZPRES(JK)*0.5) THEN
            ZZW(JK)= ICEP%X0DEPI/(ICED%XLBI*ZAI(JK)) *(ZZW2(JK)/ZRHODREF(JK))**(1.+ICED%XLBEXI) * &
               & (PTSTEP*MAX(ICED%XRTMIN(4)/PTSTEP,ZRIS(JK))*ZW2D(JK) )**(-ICED%XLBEXI)
            ZZW(JK)=  MAX(-ZRIS(JK)*ZW2D(JK)*(1.-ZSIFRC(JK))+ZZW(JK)*ZSSIO(JK)* ZSIFRC(JK)* ZXW2D13(JK), &
          &  ZZW(JK)* ( ZSSIO(JK)* ZSIFRC(JK)* ZXW2D13(JK)  + ZCITRED23*ZSSIU(JK)* (1.-ZSIFRC(JK)) ))

            ZRIS(JK) = ZRIS(JK) + ZZW(JK)
            ZRVS(JK) = ZRVS(JK) - ZZW(JK)  ! Budget here: ! cloud ice + vapor = const
            ZTHS(JK) = ZTHS(JK) + ZZW(JK)*ZLSFACT(JK) ! f(L_f*(RCBERI))
          END IF
        END DO
      ELSE

        DO JK=1,KSIZE

          ZTC =  MAX(-18.,MIN(-1.,ZZT(JK)-CST%XTT))
          ZHU =  MIN(0.15,MAX(0.,ZSSI(JK)))
          ZCRYSHA(JK)=1.1+ 3.*ZHU*(1.+ SIN(0.64*ZTC -1.3))
!         icedensity*4/3 *pi /8. =366.5 ; icedensity=700 kg/m3
          ZQIMAX = 366.5 * ZDICRIT**3 * ZCIT(JK)*ZCITRED/ZRHODREF(JK)
          ZCI2S(JK) = 0.
          IF(ZRIS(JK)*PTSTEP > 1.0e-12)THEN
              ZCI2S(JK)  =  ZRIS(JK)*(1. - MIN(1., 0.5*ZQIMAX /ZRIS(JK)/PTSTEP))* &
                  &  (1.-ZSIFRC(JK))*ZW2D(JK)
          ENDIF
        ENDDO

        DO JK = 1, KSIZE
          IF (ZCIT(JK)>0.0 .AND. ZESI(JK) < ZPRES(JK)*0.5) THEN
            ZZWC(JK)=ZCRYSHA(JK)*0.878/ZAI(JK)*(ZCIT(JK)/ZRHODREF(JK))**0.667 &
                 &*(MAX(ICED%XRTMIN(4)/PTSTEP,ZRIS(JK))*PTSTEP*ZW2D(JK))**0.333
!     Ice supersaturated part of grid box:
            IF (ZSSIO(JK)>0. .AND. ZSIFRC(JK) > 0.02_JPRB) THEN
              ZZW(JK)  = ZZWC(JK)*ZXW2D13(JK)*ZSSIO(JK)
              ZRIS(JK) = ZRIS(JK) + ZZW(JK)*ZSIFRC(JK)
              ZRVS(JK) = ZRVS(JK) - ZZW(JK)*ZSIFRC(JK)  ! Budget here: ! cloud ice + vapor = const
              ZTHS(JK) = ZTHS(JK) + ZZW(JK)*ZLSFACT(JK)*ZSIFRC(JK) ! f(L_f*(RCBERI))
            END IF
!    Ice subsaturated part of grid box:
            IF (ZSSIU(JK)<0. .AND. ZSIFRC(JK) <0.98_JPRB) THEN
              ZRIS(JK) = ZRIS(JK) - ZCI2S(JK)
              ZRSS(JK) = ZRSS(JK) + ZCI2S(JK)
              ZZW(JK)  = ZZWC(JK)*ZCITRED23*ZSSIU(JK)
              ZRIS(JK) = ZRIS(JK) + ZZW(JK)*(1.-ZSIFRC(JK))
              ZRVS(JK) = ZRVS(JK) - ZZW(JK)*(1.-ZSIFRC(JK))
              ZTHS(JK) = ZTHS(JK) + ZZW(JK)*ZLSFACT(JK)*(1.-ZSIFRC(JK))
            END IF
          END IF
        END DO
      ENDIF

    ELSE ! End OCND2

      DO JK = 1, KSIZE
        IF ((ZRCS(JK) > 0.0) .AND. &
            (ZSSI(JK) > 0.0) .AND. &
            (ZRIT(JK) > ICED%XRTMIN(4)) .AND. &
            (ZCIT(JK) > 0.0)) THEN
          ZZW(JK) = MIN(1.E8,ICED%XLBI*( ZRHODREF(JK)*ZRIT(JK)/ZCIT(JK) )**ICED%XLBEXI) ! Lbda_i
          ZZW(JK) = MIN( ZRCS(JK),( ZSSI(JK) / (ZRHODREF(JK)*ZAI(JK)) ) * ZCIT(JK) * &
                        ( ICEP%X0DEPI/ZZW(JK) + ICEP%X2DEPI*ZCJ(JK)*ZCJ(JK)/ZZW(JK)**(ICED%XDI+2.0) ) )
          ZRCS(JK) = ZRCS(JK) - ZZW(JK)
          ZRIS(JK) = ZRIS(JK) + ZZW(JK)
          ZTHS(JK) = ZTHS(JK) + ZZW(JK)*(ZLSFACT(JK)-ZLVFACT(JK))
        END IF
      END DO
    ENDIF

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),   &
                                     4,'BERFI_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     7,'BERFI_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     9,'BERFI_BU_RRI',YDDDH, YDLDDH, YDMDDH)

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RI',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_FAST_RI

END MODULE MODE_RAIN_ICE_OLD_FAST_RI
