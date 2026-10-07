!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_SLOW

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_SLOW(D, CST, ICED, ICEP, ICE_T_PARAMETERS, BUCONF, &
                               KSIZE, OCND2, OICE_T, LMODICEDEP,             &
                               PTSTEP, ZREDSN,                               &
                               GMICRO, PRHODJ, PTHS, PRVS,                   &
                               PRCT, PRRT, PRIT, PRRS,                       &
                               PRGS, PRST, PRGT, PCIT,                       &
                               PRHODREF, PZRHODJ, PLDBAS,                    &
                               PZT, PLSFACT, PLVFACT, PPRES, PSSI,           &
                               PZRVS, PRCS, PRIS, PRSS, PZTHS,               &
                               PLBDAG, PKA, PDV,                             &
                               PAI, PCJ, PAA2, PBB3,                         &
                               ZDICRIT, ZREDGR, ZKVO, PNT_C, PPRS_SDE,       &
                               TBUDGETS, KBUDGETS)

    USE YOMHOOK,              ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
    USE MODD_CST,             ONLY: CST_T
    USE MODD_CST,             ONLY: XPI, XRHOLW
    USE MODD_PRECISION,        ONLY: MNHREAL64
    USE MODD_ICET_PARAM,       ONLY: XNT_IN, ICET_PARAM_t
    USE MODD_ICET_PARAM,       ONLY: NTB_IN, NBC, NTB_C, NTB_R, NTB_R1, XR_R, XR_C
    USE MODD_RAIN_ICE_PARAM_n,  ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR_n,  ONLY: RAIN_ICE_DESCR_T

    USE MODD_BUDGET,     ONLY: TBUDGETDATA_PTR, TBUDGETCONF_t, NBUDGET_TH, NBUDGET_RG, NBUDGET_RR, NBUDGET_RC, &
                               NBUDGET_RI, NBUDGET_RS, NBUDGET_RV

    USE MODE_RAIN_ICE_OLD_ICENUMBER2, ONLY: ICENUMBER2

    IMPLICIT NONE

    TYPE(DIMPHYEX_T), INTENT(IN)       :: D
    TYPE(CST_T), INTENT(IN)            :: CST
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_T), INTENT(IN) :: ICED
    TYPE(ICET_PARAM_t),     INTENT(IN) :: ICE_T_PARAMETERS
    TYPE(TBUDGETCONF_t),    INTENT(IN) :: BUCONF

    INTEGER, INTENT(IN) :: KSIZE
    LOGICAL, INTENT(IN) :: OCND2
    LOGICAL, INTENT(IN) :: OICE_T
    LOGICAL, INTENT(IN) :: LMODICEDEP ! Logical switch for alternative dep/evap of ice

    REAL, INTENT(IN) :: PTSTEP ! Double Time step (single if cold start)

    REAL, INTENT(IN) :: ZREDSN

    LOGICAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRHODJ  ! Dry density * Jacobian
    REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTHS    ! Theta source
    REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRVS    ! Water vapor m.r. source

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
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PNT_C

    REAL(KIND=MNHREAL64), DIMENSION(KSIZE), INTENT(OUT) :: PPRS_SDE

    REAL, INTENT(IN) :: ZDICRIT, ZREDGR ! Possible reduction of the rate of graupel,snow growth
    REAL, INTENT(IN) :: ZKVO ! factor used for caluclate maximum mass in the ice distubution.

    TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
    INTEGER, INTENT(IN) :: KBUDGETS
    REAL, DIMENSION(KSIZE) :: ZBFT ! Mean time for a pristine ice crystal to reach
                                   ! size of an snow/graupel particle (ZDICRIT)
    REAL, DIMENSION(KSIZE) :: ZCRIAUTI ! Snow-to-ice autoconversion thres.
    REAL, DIMENSION(KSIZE) :: ZZW      ! Work array
    REAL, DIMENSION(KSIZE) :: ZZW2     ! Work array
    REAL, DIMENSION(D%NIJT,D%NKT) :: ZWKBUD

    REAL(KIND=MNHREAL64), DIMENSION(KSIZE) :: ZPRG_RFZ, ZPRI_RFZ, ZPRI_WFZ
    REAL(KIND=MNHREAL64) :: ZLAM_R, ZLAM_EXP, ZN0_EXP
    REAL :: ZHOMFRZ
    REAL :: ZTEMPC
    REAL :: ZNI
    INTEGER :: IDX_TC, IDX_IN, IDX_C, IDX_N, IDX_R, IDX_R1
    INTEGER :: II, IC, IR, IX
    INTEGER :: J

    INTEGER :: JL
    REAL    :: ZINVTSTEP

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!*       3.2     compute the homogeneous nucleation source: RCHONI

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SLOW',0,ZHOOK_HANDLE)

    ZINVTSTEP=1./PTSTEP
    ZZW(:) = 0.0

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'HON', ZWKBUD)
    IF (BUCONF%LBUDGET_RC) ZWKBUD(:,:) = UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'HON', ZWKBUD)
    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'HON', ZWKBUD)

!*********************BJKE begin**********************************
! Added Bigg (1953) freezing as in Thompson et al. (2008), yet with
! diagnosed number concentration for rain. This is a new process in addition
! to the ones we already have.
    IF (OICE_T) THEN
      ZHOMFRZ=-38.
      DO JL=1, KSIZE
        ZTEMPC = PZT(JL) - 273.15
        IDX_TC = MAX(1, MIN(NINT(-ZTEMPC), 45) )

        ZNI = 1.0 *1000.
        !..Ice nuclei lookup table index.
        IF (ZNI .GT. XNT_IN(1)) THEN
          II = NINT(ALOG10(ZNI))
          DO J = II-1, II+1
            IX = J
            IF ((ZNI/10.**J) .GE. 1.0 .AND. (ZNI/10.**J) .LT. 10.0) EXIT
          ENDDO
          IDX_IN = INT(ZNI/10.**IX) + 10*(IX - ICE_T_PARAMETERS%NIIN2) - (IX-ICE_T_PARAMETERS%NIIN2)
          IDX_IN = MAX(1, MIN(IDX_IN, NTB_IN))
        ELSE
          IDX_IN = 1
        ENDIF
        IF( (PRRT(JL) .GT. XR_R(1)) .AND. (PRRS(JL)>0.)) THEN
          ! Calculate rain drop number concentration
          ZLAM_R = SQRT(SQRT(XPI*XRHOLW*ICED%XCCR/(PRRT(JL)*PRHODREF(JL))))
          IR = NINT(ALOG10(PRRT(JL)))
          DO J = IR-1, IR+1
            IX = J
            IF ((PRRT(JL)/10.**J) .GE. 1.0 .AND. (PRRT(JL)/10.**J) .LT. 10.0) EXIT
          ENDDO
          IDX_R = INT(PRRT(JL)/10.**IX) + 10*(IX - ICE_T_PARAMETERS%NIR2) - (IX - ICE_T_PARAMETERS%NIR2)
          IDX_R = MAX(1, MIN(IDX_R, NTB_R))

          ZLAM_EXP = ZLAM_R * (ICE_T_PARAMETERS%XCR_GM(3)*ICE_T_PARAMETERS%XORG2*ICE_T_PARAMETERS%XORG1)**ICED%XBR
          ZN0_EXP = ICE_T_PARAMETERS%XORG1*PRRT(JL)/ICED%XAR * ZLAM_EXP**ICE_T_PARAMETERS%XCR_EX(1)
          IR = NINT(DLOG10(ZN0_EXP))
          DO J = IR-1, IR+1
            IX = J
            IF ((ZN0_EXP/10.**J) .GE. 1.0 .AND. (ZN0_EXP/10.**J) .LT. 10.0) EXIT
          ENDDO
          IDX_R1 = INT(ZN0_EXP/10.**IX) + 10*(IX - ICE_T_PARAMETERS%NIR3) - (IX - ICE_T_PARAMETERS%NIR3)
          IDX_R1 = MAX(1, MIN(IDX_R1, NTB_R1))
        ELSE
          IDX_R = 1
          IDX_R1 = NTB_R1
        ENDIF

        IF( (PRRT(JL) .GT. XR_R(1)) .AND. (PRRS(JL)>0.)) THEN
          ZPRG_RFZ(JL) = ICE_T_PARAMETERS%XTPG_QRFZ(IDX_R,IDX_R1, IDX_TC, IDX_IN)*ICE_T_PARAMETERS%XODTS
          ZPRI_RFZ(JL) = ICE_T_PARAMETERS%XTPI_QRFZ(IDX_R,IDX_R1, IDX_TC, IDX_IN)*ICE_T_PARAMETERS%XODTS
          ! Budget: Ice is created from rain
          PRIS(JL) = PRIS(JL) + ZPRI_RFZ(JL)
          PRRS(JL) = PRRS(JL) - ZPRI_RFZ(JL)
          ! Budget: Graupel is created from rain
          PRGS(JL) = PRGS(JL) + ZPRG_RFZ(JL)
          PRRS(JL) = PRRS(JL) - ZPRG_RFZ(JL)
          ! Budget: Latent heat release
          PZTHS(JL) = PZTHS(JL) + ZPRI_RFZ(JL)*(PLSFACT(JL) - PLVFACT(JL))
          PZTHS(JL) = PZTHS(JL) + ZPRG_RFZ(JL)*(PLSFACT(JL) - PLVFACT(JL))
        ENDIF
        !..Cloud water lookup table index.
        IF (PRCT(JL) .GT. XR_C(1).AND. (PRCS(JL)>0.)) THEN
          IC = NINT(ALOG10(PRCT(JL)))
          DO J = IC-1, IC+1
            IX = J
            IF ( (PRCT(JL)/10.**J).GE.1.0 .AND. &
                 (PRCT(JL)/10.**J).LT.10.0) EXIT
          ENDDO
          IDX_C = INT(PRCT(JL)/10.**IX) + 10*(IX - ICE_T_PARAMETERS%NIC2) - (IX - ICE_T_PARAMETERS%NIC2)
          IDX_C = MAX(1, MIN(IDX_C, NTB_C))
        ELSE
          IDX_C = 1
        ENDIF
        !..Cloud droplet number lookup table index.
        IDX_N = NINT(1.0 + FLOAT(NBC) * DLOG(PNT_C(JL)/ICE_T_PARAMETERS%XT_NC(1)) / ICE_T_PARAMETERS%NIC1)
        IDX_N = MAX(1, MIN(IDX_N, NBC))

        IF((PRCT(JL) .GT. XR_C(1)) .AND. (PRCS(JL) > 0.)) THEN
          ZPRI_WFZ(JL) = ICE_T_PARAMETERS%XTPI_QCFZ(IDX_C, IDX_N, IDX_TC, IDX_IN)*ICE_T_PARAMETERS%XODTS
          ZPRI_WFZ(JL) = MIN(DBLE(PRCT(JL)*ICE_T_PARAMETERS%XODTS), ZPRI_WFZ(JL))
          ! Budget: Ice is created by cloud water
          PRIS(JL) = PRIS(JL) + ZPRI_WFZ(JL)
          PRCS(JL) = PRCS(JL) - ZPRI_WFZ(JL)
          PZTHS(JL) = PZTHS(JL) + ZPRI_WFZ(JL)*(PLSFACT(JL) - PLVFACT(JL))
        ENDIF
      ENDDO
    ENDIF
!*************************BJKE out*********************************

    DO JL = 1, KSIZE
      IF ((PZT(JL)<CST%XTT-35.0) .AND. (PRCT(JL)>ICED%XRTMIN(2)) .AND. (PRCS(JL)>0.)) THEN
        ZZW(JL) = MIN( PRCS(JL),ICEP%XHON*PRHODREF(JL)*PRCT(JL)       &
                                     *EXP(ICEP%XALPHA3*(PZT(JL)-CST%XTT)-ICEP%XBETA3))
        PRIS(JL) = PRIS(JL) + ZZW(JL)
        PRCS(JL) = PRCS(JL) - ZZW(JL)
        PZTHS(JL) = PZTHS(JL) + ZZW(JL)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RCHONI))
      END IF
    END DO

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'HON', ZWKBUD)
    IF (BUCONF%LBUDGET_RC) ZWKBUD(:,:) = UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'HON', ZWKBUD)
    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'HON', ZWKBUD)

!*       3.3     compute the spontaneous freezing source: RRHONG

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'SFR', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'SFR', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'SFR', ZWKBUD)

    ZZW(:) = 0.0
    DO JL = 1, KSIZE
      IF ((PZT(JL)<CST%XTT-35.0) .AND. (PRRT(JL)>ICED%XRTMIN(3)) .AND. (PRRS(JL)>0.)) THEN
        ZZW(JL) = MIN( PRRS(JL),PRRT(JL)* ZINVTSTEP )
        PRGS(JL) = PRGS(JL) + ZZW(JL)
        PRRS(JL) = PRRS(JL) - ZZW(JL)
        PZTHS(JL) = PZTHS(JL) + ZZW(JL)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RRHONG))
      END IF
    END DO

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'SFR', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'SFR', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'SFR', ZWKBUD)

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
       PAI(:) = EXP( CST%XALPI - CST%XBETAI/PZT(:) - CST%XGAMI*LOG(PZT(:) ) ) ! es_i
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

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'DEPS', ZWKBUD)
    IF (BUCONF%LBUDGET_RV) ZWKBUD(:,:) = UNPACK(PZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_RV) CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'DEPS', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'DEPS', ZWKBUD)

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

    IF(OICE_T)THEN
      DO JL = 1, KSIZE
        ! needed for Thompson snow collecting cloud water
        PPRS_SDE(JL) = ZZW(JL)
      ENDDO
    ENDIF

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'DEPS', ZWKBUD)
    IF (BUCONF%LBUDGET_RV) ZWKBUD(:,:) = UNPACK(PZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_RV) CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'DEPS', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'DEPS', ZWKBUD)

!*       3.4.4  compute the aggregation on r_s: RIAGGS

    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'AGGS', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'AGGS', ZWKBUD)

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

    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'AGGS', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'AGGS', ZWKBUD)

!*       3.4.5  compute the autoconversion of r_i for r_s production: RIAUTS

    ZCRIAUTI(:)=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(PZT(:)-CST%XTT)+ICEP%XBCRIAUTI))
    ZZW(:) = 0.0
    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'AUTS', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'AUTS', ZWKBUD)
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
        IF ((PRIS(JL)>0.0) .AND. (PSSI(JL)>0.001)) THEN
          ZBFT(JL) = 0.5*87.5*(ZDICRIT)**2*PAI(JL)/ PSSI(JL)
          ZBFT(JL) = PTSTEP/ MAX(PTSTEP,ZBFT(JL)*2.)
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

    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'AUTS', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'AUTS', ZWKBUD)

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

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'DEPG', ZWKBUD)
    IF (BUCONF%LBUDGET_RV) ZWKBUD(:,:) = UNPACK(PZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_RV) CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'DEPG', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'DEPG', ZWKBUD)

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

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'DEPG', ZWKBUD)
    IF (BUCONF%LBUDGET_RV) ZWKBUD(:,:) = UNPACK(PZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_RV) CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'DEPG', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'DEPG', ZWKBUD)

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SLOW',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_SLOW

END MODULE MODE_RAIN_ICE_OLD_SLOW
