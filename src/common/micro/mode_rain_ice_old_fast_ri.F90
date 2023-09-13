!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_FAST_RI

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_FAST_RI(D, CST, ICEP, ICED, BUCONF,   &
                                  PTSTEP, KSIZE,                &
                                  OCND2, LMODICEDEP, GMICRO,    &
                                  PRHODJ, PTHS,                 &
                                  PRIT, PCIT,                   &
                                  PRVS, PRCS, PRIS, PRSS, PZTHS, &
                                  PRHODREF, PZRHODJ,             &
                                  PLSFACT, PLVFACT,             &
                                  PAI, PCJ,                     &
                                  PSSIO, PSSIU, PW2D, PXW2D13,  &
                                  PZT, PRES, PSSI,             &
                                  PSIFRC, PESI,                 &
                                  PCITRED, PCITRED23, PDICRIT,  &
                                  TBUDGETS, KBUDGETS)

    USE PARKIND1,            ONLY: JPRB
    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_T
    USE MODD_CST,            ONLY: CST_T
    USE MODD_RAIN_ICE_PARAM_N, ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR_N, ONLY: RAIN_ICE_DESCR_T

    USE MODE_RAIN_ICE_OLD_ICENUMBER2, ONLY: ICENUMBER2

    USE MODE_BUDGET_PHY,     ONLY: BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY
    USE MODD_BUDGET,         ONLY: TBUDGETDATA, TBUDGETCONF_t, &
                                   NBUDGET_TH, NBUDGET_RC, NBUDGET_RI

    IMPLICIT NONE

    TYPE(DIMPHYEX_T),       INTENT(IN) :: D
    TYPE(CST_T),            INTENT(IN) :: CST
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED
    TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF

    REAL,    INTENT(IN) :: PTSTEP  ! Double Time step
    INTEGER, INTENT(IN) :: KSIZE

    LOGICAL, INTENT(IN) :: OCND2
    LOGICAL, INTENT(IN) :: LMODICEDEP

    LOGICAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODJ ! Dry density * Jacobian
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PTHS   ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PCIT     ! Pristine ice conc. at t

    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRVS     ! Water vapor m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRCS     ! Cloud water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRIS     ! Pristine ice m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRSS     ! Snow/aggregate m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PZTHS    ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRHODREF ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PZRHODJ  ! RHO times Jacobian
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PAI      ! Thermodynamical function
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PSSIO    ! Super-saturation with respect to ice in the
                                                      ! supersaturated fraction of gridbox
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PSSIU    ! Sub-saturation with respect to ice in the
                                                      ! sub-saturated fraction of gridbox
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PW2D     ! Factor for subgridscale calculations
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PXW2D13  ! ZXW2D**0.333 or other expression for LMODICEDEP=T

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PZT      ! Temperature
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRES     ! Pressure
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PSSI     ! Supersaturation over ice

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PSIFRC   ! subgridscale fraction with supersaturation with
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PESI     ! saturation pressure over ice

    REAL, INTENT(IN) :: PCITRED
    REAL, INTENT(IN) :: PCITRED23
    REAL, INTENT(IN) :: PDICRIT

    TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
    INTEGER, INTENT(IN) :: KBUDGETS

    REAL, DIMENSION(KSIZE) :: ZCRYSHA ! Ice crystal habit factor
    REAL, DIMENSION(KSIZE) :: ZCI2S   ! factor to turn cloud ice with few lagre crystals into snow
    REAL, DIMENSION(KSIZE) :: ZZW     ! Work array
    REAL, DIMENSION(KSIZE) :: ZZWC    ! Work array
    REAL, DIMENSION(KSIZE) :: ZZW2    ! Work array

    REAL :: ZTC
    REAL :: ZQIMAX
    REAL :: ZHU

    INTEGER :: JK, JL

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!*       7.1    cloud ice melting
!

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RI',0,ZHOOK_HANDLE)

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'IMLT', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'IMLT', UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'IMLT', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    ZZW(:) = 0.0
    DO JK = 1, KSIZE
      IF ((PRIS(JK)>0.0) .AND. (PZT(JK)>CST%XTT)) THEN
        ZZW(JK)  = PRIS(JK)
        PRCS(JK) = PRCS(JK) + PRIS(JK)
        PZTHS(JK) = PZTHS(JK) - PRIS(JK)*(PLSFACT(JK)-PLVFACT(JK)) ! f(L_f*(-RIMLTC))
        PRIS(JK) = 0.0
        PCIT(JK) = 0.0
      END IF
    END DO

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'IMLT', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'IMLT', UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'IMLT', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

!*       7.2    Bergeron-Findeisen effect: RCBERI

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'BERFI', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'BERFI', UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'BERFI', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    ZZW(:) = 0.0
    IF(OCND2)THEN

      ! Sub gridscale decomposition into a supersaturation part of the gridbox,
      ! PSIFRC with a superaturation PSSIO and a subsaturated part (1.- PSIFRC)
      ! with a (negative) superaturation of PSSIU

      IF (LMODICEDEP) THEN

        DO JL=1,KSIZE
          ZZW2(JL) = MAX(PCIT(JL),ICENUMBER2(PRIS(JL)*PTSTEP,PZT(JL))* &
          PRHODREF(JL))
        ENDDO

        DO JK = 1, KSIZE
          IF (ZZW2(JK)>0.0 .AND. PESI(JK) < PRES(JK)*0.5) THEN
            ZZW(JK)= ICEP%X0DEPI/(ICED%XLBI*PAI(JK)) *(ZZW2(JK)/PRHODREF(JK))**(1.+ICED%XLBEXI) * &
               & (PTSTEP*MAX(ICED%XRTMIN(4)/PTSTEP,PRIS(JK))*PW2D(JK) )**(-ICED%XLBEXI)
            ZZW(JK)=  MAX(-PRIS(JK)*PW2D(JK)*(1.-PSIFRC(JK))+ZZW(JK)*PSSIO(JK)* PSIFRC(JK)* PXW2D13(JK), &
          &  ZZW(JK)* ( PSSIO(JK)* PSIFRC(JK)* PXW2D13(JK)  + PCITRED23*PSSIU(JK)* (1.-PSIFRC(JK)) ))

            PRIS(JK) = PRIS(JK) + ZZW(JK)
            PRVS(JK) = PRVS(JK) - ZZW(JK)  ! Budget here: ! cloud ice + vapor = const
            PZTHS(JK) = PZTHS(JK) + ZZW(JK)*PLSFACT(JK) ! f(L_f*(RCBERI))
          END IF
        END DO
      ELSE

        DO JK=1,KSIZE

          ZTC =  MAX(-18.,MIN(-1.,PZT(JK)-CST%XTT))
          ZHU =  MIN(0.15,MAX(0.,PSSI(JK)))
          ZCRYSHA(JK)=1.1+ 3.*ZHU*(1.+ SIN(0.64*ZTC -1.3))
!         icedensity*4/3 *pi /8. =366.5 ; icedensity=700 kg/m3
          ZQIMAX = 366.5 * PDICRIT**3 * PCIT(JK)*PCITRED/PRHODREF(JK)
          ZCI2S(JK) = 0.
          IF(PRIS(JK)*PTSTEP > 1.0e-12)THEN
              ZCI2S(JK)  =  PRIS(JK)*(1. - MIN(1., 0.5*ZQIMAX /PRIS(JK)/PTSTEP))* &
                  &  (1.-PSIFRC(JK))*PW2D(JK)
          ENDIF
        ENDDO

        DO JK = 1, KSIZE
          IF (PCIT(JK)>0.0 .AND. PESI(JK) < PRES(JK)*0.5) THEN
            ZZWC(JK)=ZCRYSHA(JK)*0.878/PAI(JK)*(PCIT(JK)/PRHODREF(JK))**0.667 &
                 &*(MAX(ICED%XRTMIN(4)/PTSTEP,PRIS(JK))*PTSTEP*PW2D(JK))**0.333
!     Ice supersaturated part of grid box:
            IF (PSSIO(JK)>0. .AND. PSIFRC(JK) > 0.02_JPRB) THEN
              ZZW(JK)  = ZZWC(JK)*PXW2D13(JK)*PSSIO(JK)
              PRIS(JK) = PRIS(JK) + ZZW(JK)*PSIFRC(JK)
              PRVS(JK) = PRVS(JK) - ZZW(JK)*PSIFRC(JK)  ! Budget here: ! cloud ice + vapor = const
              PZTHS(JK) = PZTHS(JK) + ZZW(JK)*PLSFACT(JK)*PSIFRC(JK) ! f(L_f*(RCBERI))
            END IF
!    Ice subsaturated part of grid box:
            IF (PSSIU(JK)<0. .AND. PSIFRC(JK) <0.98_JPRB) THEN
              PRIS(JK) = PRIS(JK) - ZCI2S(JK)
              PRSS(JK) = PRSS(JK) + ZCI2S(JK)
              ZZW(JK)  = ZZWC(JK)*PCITRED23*PSSIU(JK)
              PRIS(JK) = PRIS(JK) + ZZW(JK)*(1.-PSIFRC(JK))
              PRVS(JK) = PRVS(JK) - ZZW(JK)*(1.-PSIFRC(JK))
              PZTHS(JK) = PZTHS(JK) + ZZW(JK)*PLSFACT(JK)*(1.-PSIFRC(JK))
            END IF
          END IF
        END DO
      ENDIF

    ELSE ! End OCND2

      DO JK = 1, KSIZE
        IF ((PRCS(JK) > 0.0) .AND. &
            (PSSI(JK) > 0.0) .AND. &
            (PRIT(JK) > ICED%XRTMIN(4)) .AND. &
            (PCIT(JK) > 0.0)) THEN
          ZZW(JK) = MIN(1.E8,ICED%XLBI*( PRHODREF(JK)*PRIT(JK)/PCIT(JK) )**ICED%XLBEXI) ! Lbda_i
          ZZW(JK) = MIN( PRCS(JK),( PSSI(JK) / (PRHODREF(JK)*PAI(JK)) ) * PCIT(JK) * &
                        ( ICEP%X0DEPI/ZZW(JK) + ICEP%X2DEPI*PCJ(JK)*PCJ(JK)/ZZW(JK)**(ICED%XDI+2.0) ) )
          PRCS(JK) = PRCS(JK) - ZZW(JK)
          PRIS(JK) = PRIS(JK) + ZZW(JK)
          PZTHS(JK) = PZTHS(JK) + ZZW(JK)*(PLSFACT(JK)-PLVFACT(JK))
        END IF
      END DO
    ENDIF

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'BERFI', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'BERFI', UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'BERFI', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RI',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_FAST_RI

END MODULE MODE_RAIN_ICE_OLD_FAST_RI
