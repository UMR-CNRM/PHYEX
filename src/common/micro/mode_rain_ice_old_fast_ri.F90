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
  WHERE( (ZRIS(:)>0.0) .AND. (ZZT(:)>CST%XTT) )
    ZZW(:)  = ZRIS(:)
    ZRCS(:) = ZRCS(:) + ZRIS(:)
    ZTHS(:) = ZTHS(:) - ZRIS(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(-RIMLTC))
    ZRIS(:) = 0.0
    ZCIT(:) = 0.0
  END WHERE

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

       WHERE( ZZW2(:)>0.0 .AND. ZESI(:) < ZPRES(:)*0.5)
          ZZW(:)= ICEP%X0DEPI/(ICED%XLBI*ZAI(:)) *(ZZW2(:)/ZRHODREF(:))**(1.+ICED%XLBEXI) * &
             & (PTSTEP*MAX(ICED%XRTMIN(4)/PTSTEP,ZRIS(:))*ZW2D(:) )**(-ICED%XLBEXI)
          ZZW(:)=  MAX(-ZRIS(:)*ZW2D(:)*(1.-ZSIFRC(:))+ZZW(:)*ZSSIO(:)* ZSIFRC(:)* ZXW2D13(:), &
        &  ZZW(:)* ( ZSSIO(:)* ZSIFRC(:)* ZXW2D13(:)  + ZCITRED23*ZSSIU(:)* (1.-ZSIFRC(:)) ))

          ZRIS(:) = ZRIS(:) + ZZW(:)
          ZRVS(:) = ZRVS(:) - ZZW(:)  ! Budget here: ! cloud ice + vapor = const
          ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:) ! f(L_f*(RCBERI))

       END WHERE

     ELSE

      DO JK=1,KSIZE

        ZTC =  MAX(-18.,MIN(-1.,ZZT(JK)-CST%XTT))
        ZHU =  MIN(0.15,MAX(0.,ZSSI(JK)))
        ZCRYSHA(JK)=1.1+ 3.*ZHU*(1.+ SIN(0.64*ZTC -1.3))
!       icedensity*4/3 *pi /8. =366.5 ; icedensity=700 kg/m3
        ZQIMAX = 366.5 * ZDICRIT**3 * ZCIT(JK)*ZCITRED/ZRHODREF(JK)
        ZCI2S(JK) = 0.
        IF(ZRIS(JK)*PTSTEP > 1.0e-12)THEN
            ZCI2S(JK)  =  ZRIS(JK)*(1. - MIN(1., 0.5*ZQIMAX /ZRIS(JK)/PTSTEP))* &
                &  (1.-ZSIFRC(JK))*ZW2D(JK)
        ENDIF

      ENDDO

      WHERE( ZCIT(:)>0.0 .AND. ZESI(:) < ZPRES(:)*0.5)
        ZZWC(:)=ZCRYSHA(:)*0.878/ZAI(:)*(ZCIT(:)/ZRHODREF(:))**0.667 &
             &*(MAX(ICED%XRTMIN(4)/PTSTEP,ZRIS(:))*PTSTEP*ZW2D(:))**0.333
!     Ice supersaturated part of grid box:
        WHERE( ZSSIO(:)>0. .AND. ZSIFRC(:) > 0.02_JPRB )
           ZZW(:)  = ZZWC(:)*ZXW2D13(:)*ZSSIO(:)
           ZRIS(:) = ZRIS(:) + ZZW(:)*ZSIFRC(:)
           ZRVS(:) = ZRVS(:) - ZZW(:)*ZSIFRC(:)  ! Budget here: ! cloud ice + vapor = const
           ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)*ZSIFRC(:) ! f(L_f*(RCBERI))
        END WHERE
!    Ice subsaturated part of grid box:
        WHERE(  ZSSIU(:)<0. .AND. ZSIFRC(:) <0.98_JPRB )
           ZRIS(:) = ZRIS(:) - ZCI2S(:)
           ZRSS(:) = ZRSS(:) + ZCI2S(:)
           ZZW(:)  = ZZWC(:)*ZCITRED23*ZSSIU(:)
           ZRIS(:) = ZRIS(:) + ZZW(:)*(1.-ZSIFRC(:))
           ZRVS(:) = ZRVS(:) - ZZW(:)*(1.-ZSIFRC(:))
           ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)*(1.-ZSIFRC(:))
        END WHERE
      END WHERE

     ENDIF

  ELSE ! End OCND2
  WHERE( (ZRCS(:)>0.0) .AND. (ZSSI(:)>0.0) .AND. &
         (ZRIT(:)>ICED%XRTMIN(4)) .AND. (ZCIT(:)>0.0)       )
    ZZW(:) = MIN(1.E8,ICED%XLBI*( ZRHODREF(:)*ZRIT(:)/ZCIT(:) )**ICED%XLBEXI) ! Lbda_i
    ZZW(:) = MIN( ZRCS(:),( ZSSI(:) / (ZRHODREF(:)*ZAI(:)) ) * ZCIT(:) * &
                  ( ICEP%X0DEPI/ZZW(:) + ICEP%X2DEPI*ZCJ(:)*ZCJ(:)/ZZW(:)**(ICED%XDI+2.0) ) )
    ZRCS(:) = ZRCS(:) - ZZW(:)
    ZRIS(:) = ZRIS(:) + ZZW(:)
    ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:))
  END WHERE
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
