!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_SLOW

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_SLOW(D, CST, ICED, ICEP, BUCONF,         &
                               KSIZE, OCND2, LMODICEDEP,           &
                               PTSTEP, ZREDSN,                     &
                               GMICRO, PRHODJ, PTHS, PRVS,         &
                               PRCT, PRRT, PRIT, PRRS,             &
                               PRGS, PRST, PRGT, PCIT,             &
                               PRHODREF, PZRHODJ, PLDBAS,           &
                               PZT, PLSFACT, PLVFACT, PPRES, PSSI, &
                               PZRVS, PRCS, PRIS, PRSS, PZTHS,       &
                               PLBDAG, PKA, PDV,                   &
                               PAI, PCJ, PAA2, PBB3,               &
                               ZDICRIT, ZREDGR, ZKVO,              &
                               TBUDGETS, KBUDGETS)

    USE PARKIND1,             ONLY: JPRB
    USE YOMHOOK,              ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
    USE MODD_CST,             ONLY: CST_T
    USE MODD_RAIN_ICE_PARAM_n,  ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR_n,  ONLY: RAIN_ICE_DESCR_T

    USE MODE_BUDGET_PHY,      ONLY: BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY

    USE MODD_BUDGET,     ONLY: TBUDGETDATA, TBUDGETCONF_t, NBUDGET_TH, NBUDGET_RG, NBUDGET_RR, NBUDGET_RC, &
                               NBUDGET_RI, NBUDGET_RS, NBUDGET_RV

    USE MODE_RAIN_ICE_OLD_ICENUMBER2, ONLY: ICENUMBER2

    IMPLICIT NONE

    TYPE(DIMPHYEX_T), INTENT(IN)       :: D
    TYPE(CST_T), INTENT(IN)            :: CST
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_T), INTENT(IN) :: ICED
    TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF

    INTEGER, INTENT(IN) :: KSIZE
    LOGICAL, INTENT(IN) :: OCND2
    LOGICAL, INTENT(IN) :: LMODICEDEP ! Logical switch for alternative dep/evap of ice

    REAL, INTENT(IN) :: PTSTEP ! Double Time step (single if cold start)

    REAL, INTENT(IN) :: ZREDSN

    LOGICAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PRHODJ  ! Dry density * Jacobian
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PTHS    ! Theta source
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: PRVS    ! Water vapor m.r. source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRCT  ! Cloud water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRRT  ! Rain water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRIT  ! Pristine ice m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRRS  ! Rain water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRGS  ! Graupel m.r. source
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRST  ! Snow/aggregate m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRGT  ! Graupel m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PCIT  ! Pristine ice conc. at t

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRHODREF ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PZRHODJ  ! RHO times Jacobian
    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: PLDBAS   ! Slope parameter of the aggregate distribution

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PZT      ! Temperature
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PPRES    ! Pressure
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PSSI     ! Supersaturation over ice

    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PZRVS ! Water vapor m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRCS ! Cloud water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRIS ! Pristine ice m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRSS ! Snow/aggregate m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PZTHS ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: PLBDAG ! Slope parameter of the graupel distribution
    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: PKA    ! Thermal conductivity of the air
    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: PDV    ! Diffusivity of water vapor in the air

    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: PAI  ! Thermodynamical function
    REAL, DIMENSION(KSIZE), INTENT(OUT)   :: PCJ  ! Function to compute the ventilation coefficient
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PAA2 ! Part of PAI used for optimized code
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PBB3 ! Part of PAI used for optimized code

    REAL, INTENT(IN) :: ZDICRIT, ZREDGR ! Possible reduction of the rate of graupel,snow growth
    REAL, INTENT(IN) :: ZKVO ! factor used for caluclate maximum mass in the ice distubution.

    TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
    INTEGER, INTENT(IN) :: KBUDGETS
!
!*       3.2     compute the homogeneous nucleation source: RCHONI
!
    REAL, DIMENSION(KSIZE) :: ZBFT ! Mean time for a pristine ice crystal to reach
                                   ! size of an snow/graupel particle (ZDICRIT)
    REAL, DIMENSION(KSIZE) :: ZCRIAUTI ! Snow-to-ice autoconversion thres.
    REAL, DIMENSION(KSIZE) :: ZZW      ! Work array
    REAL, DIMENSION(KSIZE) :: ZZW2     ! Work array

    INTEGER :: JL
    REAL    :: ZINVTSTEP

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SLOW',0,ZHOOK_HANDLE)

    ZINVTSTEP=1./PTSTEP
    ZZW(:) = 0.0

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'HON', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'HON', UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'HON', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    DO JL = 1, KSIZE
      IF ((PZT(JL)<CST%XTT-35.0) .AND. (PRCT(JL)>ICED%XRTMIN(2)) .AND. (PRCS(JL)>0.)) THEN
        ZZW(JL) = MIN( PRCS(JL),ICEP%XHON*PRHODREF(JL)*PRCT(JL)       &
                                     *EXP(ICEP%XALPHA3*(PZT(JL)-CST%XTT)-ICEP%XBETA3))
        PRIS(JL) = PRIS(JL) + ZZW(JL)
        PRCS(JL) = PRCS(JL) - ZZW(JL)
        PZTHS(JL) = PZTHS(JL) + ZZW(JL)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RCHONI))
      END IF
    END DO

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'HON', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'HON', UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'HON', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

!*       3.3     compute the spontaneous freezing source: RRHONG

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'SFR', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'SFR', UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'SFR', UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    ZZW(:) = 0.0
    DO JL = 1, KSIZE
      IF ((PZT(JL)<CST%XTT-35.0) .AND. (PRRT(JL)>ICED%XRTMIN(3)) .AND. (PRRS(JL)>0.)) THEN
        ZZW(JL) = MIN( PRRS(JL),PRRT(JL)* ZINVTSTEP )
        PRGS(JL) = PRGS(JL) + ZZW(JL)
        PRRS(JL) = PRRS(JL) - ZZW(JL)
        PZTHS(JL) = PZTHS(JL) + ZZW(JL)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RRHONG))
      END IF
    END DO

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'SFR', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'SFR', UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'SFR', UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

!*       3.4    compute the deposition, aggregation and autoconversion sources

    PKA(:) = 2.38E-2 + 0.0071E-2 * ( PZT(:) - CST%XTT )              ! k_a
    PDV(:) = 0.211E-4 * (PZT(:)/CST%XTT)**1.94 * (CST%XP00/PPRES(:)) ! D_v
!
!*       3.4.1  compute the thermodynamical function A_i(T,P)
!*              and the c^prime_j (in the ventilation factor)
!
    IF(OCND2)THEN
       PAI(:) = PAA2(:) + PBB3(:)*PPRES(:)
    ELSE
       PAI(:) = EXP( CST%XALPI - CST%XBETAI/PZT(:) - CST%XGAMI*ALOG(PZT(:) ) ) ! es_i
       PAI(:) = ( CST%XLSTT + (CST%XCPV-CST%XCI)*(PZT(:)-CST%XTT) )**2 / (PKA(:)*CST%XRV*PZT(:)**2) &
                                   + ( CST%XRV*PZT(:) ) / (PDV(:)*PAI(:))
    ENDIF
    PCJ(:) = ICEP%XSCFAC * PRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(PZT(:)-CST%XTT) )
!
!*       3.4.3  compute the deposition on r_s: RVDEPS
!
    DO JL = 1, KSIZE
      IF (PRST(JL)>0.0) THEN
        PLDBAS(JL)  = MIN( ICED%XLBDAS_MAX,                                           &
                          ICED%XLBS*( PRHODREF(JL)*MAX( PRST(JL),ICED%XRTMIN(5) ) )**ICED%XLBEXS )
      END IF
    END DO
    ZZW(:) = 0.0

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'DEPS', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'DEPS', UNPACK(PZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'DEPS', UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    IF(OCND2)THEN
      DO JL = 1, KSIZE
        IF ((PRST(JL)>ICED%XRTMIN(5)) .AND. (PRSS(JL)>0.0)) THEN
          ZZW(JL) = (PSSI(JL)/(PRHODREF(JL)*PAI(JL))) *  &
                    (ICEP%X0DEPS*PLDBAS(JL)**ICEP%XEX0DEPS + ICEP%X1DEPS*PCJ(JL)*PLDBAS(JL)**ICEP%XEX1DEPS)
          ZZW(JL) = MIN( PZRVS(JL),MAX(-PRSS(JL),ZZW(JL)))  ! Simpler
          ZZW(JL) = ZZW(JL)*ZREDSN ! Possible tuning by using ZREDSN /=  1
          PRSS(JL) = PRSS(JL) + ZZW(JL)
          PZRVS(JL) = PZRVS(JL) - ZZW(JL)
          PZTHS(JL) = PZTHS(JL) + ZZW(JL)*PLSFACT(JL)
        END IF
      END DO
    ELSE
      DO JL = 1, KSIZE
        IF ((PRST(JL)>ICED%XRTMIN(5)) .AND. (PRSS(JL)>0.0)) THEN
          ZZW(JL) = ( PSSI(JL)/(PRHODREF(JL)*PAI(JL)) ) *          &
               ( ICEP%X0DEPS*PLDBAS(JL)**ICEP%XEX0DEPS + ICEP%X1DEPS*PCJ(JL)*PLDBAS(JL)**ICEP%XEX1DEPS )
          ZZW(JL) = MIN(PZRVS(JL),ZZW(JL)     )*(0.5+SIGN(0.5,ZZW(JL))) &
                  - MIN(PRSS(JL),ABS(ZZW(JL)))*(0.5-SIGN(0.5,ZZW(JL)))

          IF (ZZW(JL) < 0.0) THEN
            ZZW(JL) = ZZW(JL) * ICEP%XRDEPSRED
          END IF

          PRSS(JL) = PRSS(JL) + ZZW(JL)
          PZRVS(JL) = PZRVS(JL) - ZZW(JL)
          PZTHS(JL) = PZTHS(JL) + ZZW(JL)*PLSFACT(JL)
        END IF
      END DO
    ENDIF

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'DEPS', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'DEPS', UNPACK(PZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'DEPS', UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

!*       3.4.4  compute the aggregation on r_s: RIAGGS

    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'AGGS', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'AGGS', UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    ZZW(:) = 0.0
    DO JL = 1, KSIZE
      IF ((PRIT(JL)>ICED%XRTMIN(4)) .AND. (PRST(JL)>ICED%XRTMIN(5)) .AND. (PRIS(JL)>0.0)) THEN
        ZZW(JL) = MIN(PRIS(JL),ICEP%XFIAGGS * EXP( ICEP%XCOLEXIS*(PZT(JL)-CST%XTT)) &
                                            * PRIT(JL)                              &
                                            * PLDBAS(JL)**ICEP%XEXIAGGS             &
                                            * PRHODREF(JL)**(-ICED%XCEXVT))
        PRSS(JL)  = PRSS(JL)  + ZZW(JL)
        PRIS(JL)  = PRIS(JL)  - ZZW(JL)
      END IF
    END DO

    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'AGGS', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'AGGS', UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

!*       3.4.5  compute the autoconversion of r_i for r_s production: RIAUTS

    ZCRIAUTI(:)=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(PZT(:)-CST%XTT)+ICEP%XBCRIAUTI))
    ZZW(:) = 0.0
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'AUTS', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'AUTS', UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    DO JL = 1, KSIZE
      IF ((PRIT(JL)>ICED%XRTMIN(4)) .AND. (PRIS(JL)>0.0)) THEN
        ZZW(JL) = MIN(PRIS(JL),ICEP%XTIMAUTI * EXP(ICEP%XTEXAUTI*(PZT(JL)-CST%XTT)) &
                                             * MAX(PRIT(JL)-ZCRIAUTI(JL),0.0 ))
        PRSS(JL) = PRSS(JL) + ZZW(JL)
        PRIS(JL) = PRIS(JL) - ZZW(JL)
      END IF
    END DO

    IF (OCND2 .AND. .NOT. LMODICEDEP) THEN ! 3.4.5 B:

      ! Turn ice crystals lagrer than a precribed size into snow:
      ! (For the moment sperical ice crystals are assumed)

      DO JL = 1, KSIZE
        IF ((PRIS(JL)>0.0_JPRB) .AND. (PSSI(JL)>0.001_JPRB)) THEN
          ZBFT(JL) = 0.5_JPRB*87.5_JPRB*(ZDICRIT)**2*PAI(JL)/ PSSI(JL)
          ZBFT(JL) = PTSTEP/ MAX(PTSTEP,ZBFT(JL)*2._JPRB)
          PRSS(JL) = PRSS(JL) + ZBFT(JL)*PRIS(JL)
          PRIS(JL) = PRIS(JL) - ZBFT(JL)*PRIS(JL)
        END IF
      END DO
    ENDIF

    IF (OCND2 .AND. LMODICEDEP) THEN ! 3.4.5 B:

      ! Turn ice to snow if ice crystal distrubution is such that
      ! the ice crystal diameter for the (mass x N_i) maximum
      ! is lagrer than a precribed size.
      ! (ZDICRIT) The general gamma function is assumed

      DO JL=1,KSIZE
        ZZW2(JL) = &
        MAX(PCIT(JL),ICENUMBER2(PRIS(JL)*PTSTEP,PZT(JL))*PRHODREF(JL))
      ENDDO

      DO JL = 1, KSIZE
        IF (PRIS(JL)>ICEP%XFRMIN(13) .AND.PCIT(JL) > 0.) THEN
          ! LAMBDA for ICE
          ZZW2(JL) = MIN(1.E8,ICED%XLBI*(PRHODREF(JL)*PRIS(JL)* PTSTEP/ZZW2(JL))**ICED%XLBEXI)
          ZBFT(JL) = 1. - 0.5**(ZKVO /ZZW2(JL))
          ZBFT(JL) = MIN(0.9*PRIS(JL)*PTSTEP, ZBFT(JL)*PRIS(JL)*PTSTEP)
          PRSS(JL) = PRSS(JL) + ZBFT(JL)
          PRIS(JL) = PRIS(JL) - ZBFT(JL)
        END IF
      END DO
    ENDIF

    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'AUTS', UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'AUTS', UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

!*       3.4.6  compute the deposition on r_g: RVDEPG

    ZZW2(:) = 0.0
    IF (ICEP%XFRMIN(5)> 1.0E-12 .AND. ICEP%XFRMIN(6) > 0.01) THEN
      ZZW2(:) = MAX(0., MIN(1., (ICEP%XFRMIN(5) - PRGS(:))/ICEP%XFRMIN(5)))* &
              & MAX(0., MIN(1., PSSI(:)/ICEP%XFRMIN(6)))
    ENDIF


    DO JL = 1, KSIZE
      IF (PRGT(JL)>0.0) THEN
        PLBDAG(JL)  = ICED%XLBG*(PRHODREF(JL)*MAX(PRGT(JL), ICED%XRTMIN(6)))**ICED%XLBEXG
      END IF
    END DO

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'DEPG', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'DEPG', UNPACK(PZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'DEPG', UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    ZZW(:) = 0.0
    DO JL = 1, KSIZE
      IF ((PRGT(JL)>ICED%XRTMIN(6)) .AND. (PRGS(JL)>0.0)) THEN
        ZZW(JL) = (PSSI(JL)/(PRHODREF(JL)*PAI(JL))) * &
                  (ICEP%X0DEPG*PLBDAG(JL)**ICEP%XEX0DEPG + &
                   ICEP%X1DEPG*PCJ(JL)*PLBDAG(JL)**ICEP%XEX1DEPG)

        ZZW(JL) = MIN(PZRVS(JL),ZZW(JL)      )*(0.5+SIGN(0.5,ZZW(JL))) &
                - MIN(PRGS(JL),ABS(ZZW(JL)) )*(0.5-SIGN(0.5,ZZW(JL)))
        ZZW(JL) = ZZW(JL)*ZREDGR

        IF (ZZW(JL) < 0.0 ) THEN
          ZZW(JL)  = ZZW(JL) * ICEP%XRDEPGRED
        END IF

        PRSS(JL) = (ZZW(JL) + PRGS(JL))* ZZW2(JL) + PRSS(JL)
        PRGS(JL) = (ZZW(JL) + PRGS(JL))*(1. - ZZW2(JL))
        PZRVS(JL) = PZRVS(JL) - ZZW(JL)
        PZTHS(JL) = PZTHS(JL) + ZZW(JL)*PLSFACT(JL)
      END IF
    END DO

    DO JL = 1, KSIZE
      IF (ZZW(JL) < 0.0) THEN
        ZZW(JL)  = ZZW(JL) * ICEP%XRDEPGRED
      END IF
    END DO

    IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'DEPG', UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'DEPG', UNPACK(PZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'DEPG', UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0))

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SLOW',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_SLOW

END MODULE MODE_RAIN_ICE_OLD_SLOW
