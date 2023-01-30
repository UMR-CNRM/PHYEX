!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_SLOW

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_SLOW(D, CST, ICED, ICEP,                 &
                               KSIZE, OCND2, LMODICEDEP,           &
                               PTSTEP, ZREDSN,                     &
                               GMICRO, PRHODJ, PTHS, PRVS,         &
                               ZRCT, ZRRT, ZRIT, ZRRS,             &
                               ZRGS, ZRST, ZRGT, ZCIT,             &
                               ZRHODREF, ZRHODJ, ZZW2, ZLBDAS,     &
                               ZZT, ZLSFACT, ZLVFACT, ZPRES, ZSSI, &
                               ZRVS, ZRCS, ZRIS, ZRSS, ZTHS,       &
                               ZLBDAG, ZKA, ZDV,                   &
                               ZAI, ZCJ, ZAA2, ZBB3,               &
                               ZRDEPSRED, ZRDEPGRED,               &
                               ZDICRIT, ZREDGR, ZKVO,              &
                               YDDDH, YDLDDH, YDMDDH)

    USE PARKIND1,             ONLY: JPRB
    USE YOMHOOK,              ONLY: LHOOK, DR_HOOK
    USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
    USE MODD_CST,             ONLY: CST_T
    USE MODD_RAIN_ICE_PARAM,  ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR,  ONLY: RAIN_ICE_DESCR_T

    USE MODE_BUDGET,          ONLY: BUDGET_DDH
    USE DDH_MIX,              ONLY: TYP_DDH
    USE YOMLDDH,              ONLY: TLDDH
    USE YOMMDDH,              ONLY: TMDDH

    USE MODD_BUDGET,     ONLY: LBUDGET_TH, LBUDGET_RG, LBUDGET_RR, LBUDGET_RC, &
                               LBUDGET_RI, LBUDGET_RS, LBUDGET_RV

    USE MODE_RAIN_ICE_OLD_ICENUMBER2, ONLY: ICENUMBER2

    IMPLICIT NONE

    TYPE(DIMPHYEX_T), INTENT(IN)       :: D
    TYPE(CST_T), INTENT(IN)            :: CST
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_T), INTENT(IN) :: ICED

    INTEGER, INTENT(IN) :: KSIZE
    LOGICAL, INTENT(IN) :: OCND2
    LOGICAL, INTENT(IN) :: LMODICEDEP ! Logical switch for alternative dep/evap of ice

    REAL, INTENT(IN) :: PTSTEP ! Double Time step (single if cold start)

    REAL, INTENT(IN) :: ZREDSN

    LOGICAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PRHODJ  ! Dry density * Jacobian
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PTHS    ! Theta source
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PRVS    ! Water vapor m.r. source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRCT  ! Cloud water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRRT  ! Rain water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRIT  ! Pristine ice m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRRS  ! Rain water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRGS  ! Graupel m.r. source
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRST  ! Snow/aggregate m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRGT  ! Graupel m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZCIT  ! Pristine ice conc. at t

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRHODREF ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZRHODJ   ! RHO times Jacobian
    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: ZZW2     ! Work array
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZLBDAS   ! Slope parameter of the aggregate distribution

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZZT      ! Temperature
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLSFACT  ! L_s/(Pi_ref*C_ph)
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZLVFACT  ! L_v/(Pi_ref*C_ph)
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZPRES    ! Pressure
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZSSI     ! Supersaturation over ice

    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRVS ! Water vapor m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRCS ! Cloud water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRIS ! Pristine ice m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZRSS ! Snow/aggregate m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZTHS ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: ZLBDAG ! Slope parameter of the graupel distribution
    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: ZKA    ! Thermal conductivity of the air
    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: ZDV    ! Diffusivity of water vapor in the air

    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: ZAI  ! Thermodynamical function
    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: ZCJ  ! Function to compute the ventilation coefficient
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZAA2 ! Part of ZAI used for optimized code
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: ZBB3 ! Part of ZAI used for optimized code

    REAL, INTENT(IN) :: ZRDEPSRED, ZRDEPGRED ! Tuning of sublimation (neg. sublimation)
    REAL, INTENT(IN) :: ZDICRIT, ZREDGR ! Possible reduction of the rate of graupel,snow growth
    REAL, INTENT(IN) :: ZKVO ! factor used for caluclate maximum mass in the ice distubution.

    TYPE(TYP_DDH),        INTENT(INOUT)   :: YDDDH
    TYPE(TLDDH),          INTENT(IN)      :: YDLDDH
    TYPE(TMDDH),          INTENT(IN)      :: YDMDDH
!
!*       3.2     compute the homogeneous nucleation source: RCHONI
!
    REAL, DIMENSION(KSIZE) :: ZBFT ! Mean time for a pristine ice crystal to reach
                                   ! size of an snow/graupel particle (ZDICRIT)
    REAL, DIMENSION(KSIZE) :: ZCRIAUTI ! Snow-to-ice autoconversion thres.
    REAL, DIMENSION(KSIZE) :: ZZW      ! Work array

    INTEGER :: JL
    REAL    :: ZINVTSTEP

    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SLOW',0,ZHOOK_HANDLE)

    ZINVTSTEP=1./PTSTEP
    ZZW(:) = 0.0

    DO JL = 1, KSIZE
      IF ((ZZT(JL)<CST%XTT-35.0) .AND. (ZRCT(JL)>ICED%XRTMIN(2)) .AND. (ZRCS(JL)>0.)) THEN
        ZZW(JL) = MIN( ZRCS(JL),ICEP%XHON*ZRHODREF(JL)*ZRCT(JL)       &
                                     *EXP(ICEP%XALPHA3*(ZZT(JL)-CST%XTT)-ICEP%XBETA3))
        ZRIS(JL) = ZRIS(JL) + ZZW(JL)
        ZRCS(JL) = ZRCS(JL) - ZZW(JL)
        ZTHS(JL) = ZTHS(JL) + ZZW(JL)*(ZLSFACT(JL)-ZLVFACT(JL)) ! f(L_f*(RCHONI))
      END IF
    END DO

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),  &
                                     4,'HON_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),   &
                                     7,'HON_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),   &
                                     9,'HON_BU_RRI',YDDDH, YDLDDH, YDMDDH)

!*       3.3     compute the spontaneous freezing source: RRHONG

    ZZW(:) = 0.0
    DO JL = 1, KSIZE
      IF ((ZZT(JL)<CST%XTT-35.0) .AND. (ZRRT(JL)>ICED%XRTMIN(3)) .AND. (ZRRS(JL)>0.)) THEN
        ZZW(JL) = MIN( ZRRS(JL),ZRRT(JL)* ZINVTSTEP )
        ZRGS(JL) = ZRGS(JL) + ZZW(JL)
        ZRRS(JL) = ZRRS(JL) - ZZW(JL)
        ZTHS(JL) = ZTHS(JL) + ZZW(JL)*(ZLSFACT(JL)-ZLVFACT(JL)) ! f(L_f*(RRHONG))
      END IF
    END DO

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),  &
                                     4,'SFR_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),   &
                                     8,'SFR_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),   &
                                    11,'SFR_BU_RRG',YDDDH, YDLDDH, YDMDDH)

!*       3.4    compute the deposition, aggregation and autoconversion sources

    ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - CST%XTT )              ! k_a
    ZDV(:) = 0.211E-4 * (ZZT(:)/CST%XTT)**1.94 * (CST%XP00/ZPRES(:)) ! D_v
!
!*       3.4.1  compute the thermodynamical function A_i(T,P)
!*              and the c^prime_j (in the ventilation factor)
!
    IF(OCND2)THEN
       ZAI(:) = ZAA2(:) + ZBB3(:)*ZPRES(:)
    ELSE
       ZAI(:) = EXP( CST%XALPI - CST%XBETAI/ZZT(:) - CST%XGAMI*ALOG(ZZT(:) ) ) ! es_i
       ZAI(:) = ( CST%XLSTT + (CST%XCPV-CST%XCI)*(ZZT(:)-CST%XTT) )**2 / (ZKA(:)*CST%XRV*ZZT(:)**2) &
                                   + ( CST%XRV*ZZT(:) ) / (ZDV(:)*ZAI(:))
    ENDIF
    ZCJ(:) = ICEP%XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-CST%XTT) )
!
!*       3.4.3  compute the deposition on r_s: RVDEPS
!
    DO JL = 1, KSIZE
      IF (ZRST(JL)>0.0) THEN
        ZLBDAS(JL)  = MIN( ICED%XLBDAS_MAX,                                           &
                          ICED%XLBS*( ZRHODREF(JL)*MAX( ZRST(JL),ICED%XRTMIN(5) ) )**ICED%XLBEXS )
      END IF
    END DO
    ZZW(:) = 0.0

    IF(OCND2)THEN
      WHERE ( (ZRST(:)>ICED%XRTMIN(5)) .AND. (ZRSS(:)>0.0) )
        ZZW(:) = ( ZSSI(:)/(ZRHODREF(:)*ZAI(:)) ) *  &
               ( ICEP%X0DEPS*ZLBDAS(:)**ICEP%XEX0DEPS + ICEP%X1DEPS*ZCJ(:)*ZLBDAS(:)**ICEP%XEX1DEPS )
        ZZW(:) = MIN( ZRVS(:),MAX(-ZRSS(:),ZZW(:)))  ! Simpler
        ZZW(:) = ZZW(:)*ZREDSN ! Possible tuning by using ZREDSN /=  1
        ZRSS(:) = ZRSS(:) + ZZW(:)
        ZRVS(:) = ZRVS(:) - ZZW(:)
        ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
      END WHERE
    ELSE
      WHERE ( (ZRST(:)>ICED%XRTMIN(5)) .AND. (ZRSS(:)>0.0) )
        ZZW(:) = ( ZSSI(:)/(ZRHODREF(:)*ZAI(:)) ) *          &
             ( ICEP%X0DEPS*ZLBDAS(:)**ICEP%XEX0DEPS + ICEP%X1DEPS*ZCJ(:)*ZLBDAS(:)**ICEP%XEX1DEPS )
        ZZW(:) =         MIN( ZRVS(:),ZZW(:)      )*(0.5+SIGN(0.5,ZZW(:))) &
             - MIN( ZRSS(:),ABS(ZZW(:)) )*(0.5-SIGN(0.5,ZZW(:)))
        WHERE (ZZW(:) < 0.0 )
          ZZW(:) = ZZW(:) * ZRDEPSRED
        END WHERE
        ZRSS(:) = ZRSS(:) + ZZW(:)
        ZRVS(:) = ZRVS(:) - ZZW(:)
        ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
      END WHERE
    ENDIF

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),  &
                                     4,'DEPS_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RV) CALL BUDGET_DDH(UNPACK(ZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:),  &
                                     6,'DEPS_BU_RRV',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),   &
                                    10,'DEPS_BU_RRS',YDDDH, YDLDDH, YDMDDH)

!*       3.4.4  compute the aggregation on r_s: RIAGGS

    ZZW(:) = 0.0
    WHERE ( (ZRIT(:)>ICED%XRTMIN(4)) .AND. (ZRST(:)>ICED%XRTMIN(5)) .AND. (ZRIS(:)>0.0) )
      ZZW(:) = MIN( ZRIS(:),ICEP%XFIAGGS * EXP( ICEP%XCOLEXIS*(ZZT(:)-CST%XTT) ) &
                                    * ZRIT(:)                      &
                                    * ZLBDAS(:)**ICEP%XEXIAGGS     &
                                    * ZRHODREF(:)**(-ICED%XCEXVT)       )
      ZRSS(:)  = ZRSS(:)  + ZZW(:)
      ZRIS(:)  = ZRIS(:)  - ZZW(:)
    END WHERE

    IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     9,'AGGS_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    10,'AGGS_BU_RRS',YDDDH, YDLDDH, YDMDDH)

!*       3.4.5  compute the autoconversion of r_i for r_s production: RIAUTS

    ZCRIAUTI(:)=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(ZZT(:)-CST%XTT)+ICEP%XBCRIAUTI))
    ZZW(:) = 0.0
    WHERE ( (ZRIT(:)>ICED%XRTMIN(4)) .AND. (ZRIS(:)>0.0) )
      ZZW(:) = MIN( ZRIS(:),ICEP%XTIMAUTI * EXP( ICEP%XTEXAUTI*(ZZT(:)-CST%XTT) ) &
                              * MAX( ZRIT(:)-ZCRIAUTI(:),0.0 ) )
      ZRSS(:)  = ZRSS(:)  + ZZW(:)
      ZRIS(:)  = ZRIS(:)  - ZZW(:)
    END WHERE

    IF (OCND2 .AND. .NOT. LMODICEDEP) THEN ! 3.4.5 B:

      ! Turn ice crystals lagrer than a precribed size into snow:
      ! (For the moment sperical ice crystals are assumed)

      WHERE (  (ZRIS(:)>0.0_JPRB) .AND.(ZSSI(:)>0.001_JPRB) )
        ZBFT(:) =   0.5_JPRB*87.5_JPRB*(ZDICRIT)**2*ZAI(:)/ ZSSI(:)
        ZBFT(:) =   PTSTEP/ MAX(PTSTEP,ZBFT(:)*2._JPRB)
        ZRSS(:) =   ZRSS(:)  + ZBFT(:)*ZRIS(:)
        ZRIS(:) =   ZRIS(:)  - ZBFT(:)*ZRIS(:)
      END WHERE
    ENDIF

    IF (OCND2 .AND. LMODICEDEP) THEN ! 3.4.5 B:

      ! Turn ice to snow if ice crystal distrubution is such that
      ! the ice crystal diameter for the (mass x N_i) maximum
      ! is lagrer than a precribed size.
      ! (ZDICRIT) The general gamma function is assumed

      DO JL=1,KSIZE
        ZZW2(JL) = &
        MAX(ZCIT(JL),ICENUMBER2(ZRIS(JL)*PTSTEP,ZZT(JL))*ZRHODREF(JL))
      ENDDO

      WHERE (  ZRIS(:)>ICEP%XFRMIN(13) .AND.ZCIT(:) > 0. )
        ! LAMBDA for ICE
        ZZW2(:) = MIN(1.E8,ICED%XLBI*( ZRHODREF(:)*ZRIS(:)* PTSTEP/ZZW2(:) )**ICED%XLBEXI)
        ZBFT(:) = 1. - 0.5**( ZKVO /ZZW2(:))
        ZBFT(:) = MIN(0.9*ZRIS(:)*PTSTEP, ZBFT(:)*ZRIS(:)*PTSTEP)
        ZRSS(:) =   ZRSS(:)  + ZBFT(:)
        ZRIS(:) =   ZRIS(:)  - ZBFT(:)
      END WHERE
    ENDIF

    IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     9,'AUTS_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    10,'AUTS_BU_RRS',YDDDH, YDLDDH, YDMDDH)

!*       3.4.6  compute the deposition on r_g: RVDEPG

    ZZW2(:) = 0.0
    IF (ICEP%XFRMIN(5)> 1.0E-12 .AND. ICEP%XFRMIN(6) > 0.01) &
          &        ZZW2(:) = MAX(0., MIN(1., (ICEP%XFRMIN(5) - ZRGS(:))/ICEP%XFRMIN(5)))* &
          & MAX(0.,MIN(1.,ZSSI(:)/ICEP%XFRMIN(6)))


    WHERE ( ZRGT(:)>0.0 )
      ZLBDAG(:)  = ICED%XLBG*( ZRHODREF(:)*MAX( ZRGT(:),ICED%XRTMIN(6) ) )**ICED%XLBEXG
    END WHERE
    ZZW(:) = 0.0
    WHERE ( (ZRGT(:)>ICED%XRTMIN(6)) .AND. (ZRGS(:)>0.0) )
      ZZW(:) = ( ZSSI(:)/(ZRHODREF(:)*ZAI(:)) ) *                               &
               ( ICEP%X0DEPG*ZLBDAG(:)**ICEP%XEX0DEPG + ICEP%X1DEPG*ZCJ(:)*ZLBDAG(:)**ICEP%XEX1DEPG )
      ZZW(:) =         MIN( ZRVS(:),ZZW(:)      )*(0.5+SIGN(0.5,ZZW(:))) &
                     - MIN( ZRGS(:),ABS(ZZW(:)) )*(0.5-SIGN(0.5,ZZW(:)))
      ZZW(:)  = ZZW(:)*ZREDGR
      WHERE (ZZW(:) < 0.0 )
        ZZW(:)  = ZZW(:) * ZRDEPGRED
      END WHERE
      ZRSS(:) = (ZZW(:) + ZRGS(:))* ZZW2(:) + ZRSS(:)
      ZRGS(:) = (ZZW(:) + ZRGS(:))*(1. - ZZW2(:))
      ZRVS(:) = ZRVS(:) - ZZW(:)
      ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
    END WHERE
    WHERE (ZZW(:) < 0.0 )
      ZZW(:)  = ZZW(:) * ZRDEPGRED
    END WHERE

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),   &
                                     4,'DEPG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RV) CALL BUDGET_DDH(UNPACK(ZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:),   &
                                     6,'DEPG_BU_RRV',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    11,'DEPG_BU_RRG',YDDDH, YDLDDH, YDMDDH)

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SLOW',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_SLOW

END MODULE MODE_RAIN_ICE_OLD_SLOW
