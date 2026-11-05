!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_FAST_RS

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_FAST_RS(D, CST, ICEP, ICED, ICE_T_PARAMETERS, BUCONF,  &
                                  OICE_T, PTSTEP, KSIZE, KRR, GMICRO,            &
                                  PRHODJ, PTHS,                                  &
                                  PRVT, PRCT, PRRT, PRST,                        &
                                  PRRS, PRCS, PRSS, PRGS, PZTHS,                 &
                                  ZRHODREF, ZRHODJ, ZLSFACT, ZLVFACT,            &
                                  ZCJ, ZKA, ZDV,                                 &
                                  ZLBDAR, ZLBDAS, ZCOLF, ZPRES, ZZT,             &
                                  PMVD_C, PMVD_R, PPRS_SDE, PVTR, PCCR_V, PRHOF, &
                                  TBUDGETS, KBUDGETS)

    USE YOMHOOK,               ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_PRECISION,        ONLY: MNHREAL64
    USE MODD_PARAMETERS,       ONLY: JPVEXT
    USE MODD_ICET_PARAM,       ONLY: XR_THOM, XRHO_NOT
    USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_T
    USE MODD_CST,              ONLY: CST_T, XPI
    USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_T
    USE MODD_ICET_PARAM,       ONLY: ICET_PARAM_t
    USE MODD_ICET_PARAM,       ONLY: XSA, XSB, XD0S, XD0C, NBS, XICET_EPS, XTHVREFZ

    USE MODD_BUDGET,     ONLY: TBUDGETDATA_PTR, TBUDGETCONF_t, &
                               NBUDGET_TH, NBUDGET_RG, NBUDGET_RR, NBUDGET_RC, NBUDGET_RS

    IMPLICIT NONE

    TYPE(DIMPHYEX_T),       INTENT(IN) :: D
    TYPE(CST_T),            INTENT(IN) :: CST
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED
    TYPE(ICET_PARAM_t),     INTENT(IN) :: ICE_T_PARAMETERS
    TYPE(TBUDGETCONF_t),    INTENT(IN) :: BUCONF

    LOGICAL, INTENT(IN) :: OICE_T
    REAL,    INTENT(IN) :: PTSTEP  ! Double Time step
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, INTENT(IN) :: KRR

    LOGICAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian

    REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PTHS    ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRVT   ! Water vapor m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRCT   ! Cloud water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRRT   ! Rain water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRST   ! Snow/aggregate m.r. at t

    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRRS   ! Rain water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRCS   ! Cloud water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRSS   ! Snow/aggregate m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRGS   ! Graupel m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PZTHS   ! Theta source

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

    REAL, DIMENSION(KSIZE), INTENT(IN)                  :: PMVD_C
    REAL, DIMENSION(KSIZE), INTENT(IN)                  :: PMVD_R
    REAL(KIND=MNHREAL64), DIMENSION(KSIZE), INTENT(IN)  :: PPRS_SDE
    REAL, DIMENSION(KSIZE), INTENT(OUT)                 :: PVTR
    REAL, DIMENSION(KSIZE), INTENT(IN)                  :: PCCR_V
    REAL(KIND=MNHREAL64), DIMENSION(KSIZE), INTENT(OUT) :: PRHOF
    
    TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
    INTEGER, INTENT(IN) :: KBUDGETS

    LOGICAL, DIMENSION(KSIZE) :: GMASK ! Test where to compute riming/accretion

    INTEGER, DIMENSION(KSIZE) :: IVEC1 ! Vectors of indices for
    INTEGER, DIMENSION(KSIZE) :: IVEC2 ! Vectors of indices for

    REAL, DIMENSION(KSIZE) :: ZVEC1 ! Work vectors for interpolations
    REAL, DIMENSION(KSIZE) :: ZVEC2 ! Work vectors for interpolations
    REAL, DIMENSION(KSIZE) :: ZVEC3 ! Work vectors for interpolations

    REAL, DIMENSION(KSIZE)      :: ZZW      ! Work array
    REAL, DIMENSION(KSIZE, KRR) :: ZZW1     ! Work array
    REAL, DIMENSION(D%NIJT,D%NKT) :: ZWKBUD

    REAL :: ZSMOB, ZSMO2, ZSMOC, ZSMOE
    REAL :: ZTC0, ZLOGA_A, ZA_A, ZB_B
    REAL :: ZR_FRAC, ZG_FRAC
    REAL :: ZRATIO, ZFRACCSS_V
    REAL :: ZRHO00
    REAL :: ZFSACCRG
    REAL, DIMENSION(KSIZE)      :: ZX_DS
    REAL, DIMENSION(KSIZE)      :: ZVTS
    REAL, DIMENSION(KSIZE)      :: ZEF_SR
    REAL(KIND=MNHREAL64) :: ZRHO
    REAL(KIND=MNHREAL64) :: ZEF_SW
    REAL(KIND=MNHREAL64) :: ZPRG_SCW
    REAL(KIND=MNHREAL64) :: ZPRS_SCW

    INTEGER :: IGRIM, IGACC
    INTEGER, DIMENSION(KSIZE) :: I1
    INTEGER :: JL, JK
    INTEGER :: IDX

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       5.1    cloud droplet riming of the aggregates
!
    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RS',0,ZHOOK_HANDLE)

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'RIM', ZWKBUD)
    IF (BUCONF%LBUDGET_RC) ZWKBUD(:,:) = UNPACK(PRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'RIM', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'RIM', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'RIM', ZWKBUD)

    ZZW1(:,:) = 0.0

    IF (OICE_T) THEN
      DO JK=1,KSIZE
        !Important: These are necessary also for graupel collecting cloud water
        ZRHO = 0.622*ZPRES(JK)/(XR_THOM*ZZT(JK)*(PRVT(JK)+0.622)) !This might already exist
        PRHOF(JK) = SQRT(XRHO_NOT/ZRHO)
        ZPRG_SCW = 0.0
        ZPRS_SCW=0.0
        ZX_DS(JK) = 0.0
        ZTC0 = MIN(-0.1, ZZT(JK)-273.15)
        ZSMOB = PRST(JK)*ICE_T_PARAMETERS%XOAMS

        ! All other moments based on reference, 2nd moment.  If XBS.ne.2,
        ! then we must compute actual 2nd moment and use as reference.
        IF (ICED%XBS .GT. (2.0-1.e-3) .AND. ICED%XBS .LT. (2.0+1.e-3)) THEN
          ZSMO2 = ZSMOB
        ELSE
          ZLOGA_A = XSA(1) + XSA(2)*ZTC0 + XSA(3)*ICED%XBS &
                  + XSA(4)*ZTC0*ICED%XBS + XSA(5)*ZTC0*ZTC0 &
                  + XSA(6)*ICED%XBS*ICED%XBS + XSA(7)*ZTC0*ZTC0*ICED%XBS &
                  + XSA(8)*ZTC0*ICED%XBS*ICED%XBS + XSA(9)*ZTC0*ZTC0*ZTC0 &
                  + XSA(10)*ICED%XBS*ICED%XBS*ICED%XBS
          ZA_A = 10.0**ZLOGA_A
          ZB_B = XSB(1) + XSB(2)*ZTC0 + XSB(3)*ICED%XBS &
               + XSB(4)*ZTC0*ICED%XBS + XSB(5)*ZTC0*ZTC0 &
               + XSB(6)*ICED%XBS*ICED%XBS + XSB(7)*ZTC0*ZTC0*ICED%XBS &
               + XSB(8)*ZTC0*ICED%XBS*ICED%XBS + XSB(9)*ZTC0*ZTC0*ZTC0 &
               + XSB(10)*ICED%XBS*ICED%XBS*ICED%XBS
          ZSMO2 = (ZSMOB/ZA_A)**(1./ZB_B)
        ENDIF

        ! Calculate XBS+1 (th) moment.  Useful for diameter calcs.
        ZLOGA_A = XSA(1) + XSA(2)*ZTC0 + XSA(3)*ICE_T_PARAMETERS%XCS_EX(1) &
              & + XSA(4)*ZTC0*ICE_T_PARAMETERS%XCS_EX(1) + XSA(5)*ZTC0*ZTC0 &
              & + XSA(6)*ICE_T_PARAMETERS%XCS_EX(1)*ICE_T_PARAMETERS%XCS_EX(1) &
              & + XSA(7)*ZTC0*ZTC0*ICE_T_PARAMETERS%XCS_EX(1) &
              & + XSA(8)*ZTC0*ICE_T_PARAMETERS%XCS_EX(1)*ICE_T_PARAMETERS%XCS_EX(1) + XSA(9)*ZTC0*ZTC0*ZTC0 &
              & + XSA(10)*ICE_T_PARAMETERS%XCS_EX(1)*ICE_T_PARAMETERS%XCS_EX(1)*ICE_T_PARAMETERS%XCS_EX(1)
        ZA_A = 10.0**ZLOGA_A
        ZB_B = XSB(1)+ XSB(2)*ZTC0 + XSB(3)*ICE_T_PARAMETERS%XCS_EX(1) + XSB(4)*ZTC0*ICE_T_PARAMETERS%XCS_EX(1) &
           & + XSB(5)*ZTC0*ZTC0 + XSB(6)*ICE_T_PARAMETERS%XCS_EX(1)*ICE_T_PARAMETERS%XCS_EX(1) &
           & + XSB(7)*ZTC0*ZTC0*ICE_T_PARAMETERS%XCS_EX(1) &
           & + XSB(8)*ZTC0*ICE_T_PARAMETERS%XCS_EX(1)*ICE_T_PARAMETERS%XCS_EX(1) &
           & + XSB(9)*ZTC0*ZTC0*ZTC0 &
           & + XSB(10)*ICE_T_PARAMETERS%XCS_EX(1)*ICE_T_PARAMETERS%XCS_EX(1)*ICE_T_PARAMETERS%XCS_EX(1)
        ZSMOC = ZA_A * ZSMO2**ZB_B

        ! Calculate XDS+2 (th) moment.  Useful for riming.
        ZLOGA_A = XSA(1) + XSA(2)*ZTC0 + XSA(3)*ICE_T_PARAMETERS%XCS_EX(13) &
              & + XSA(4)*ZTC0*ICE_T_PARAMETERS%XCS_EX(13) + XSA(5)*ZTC0*ZTC0 &
              & + XSA(6)*ICE_T_PARAMETERS%XCS_EX(13)*ICE_T_PARAMETERS%XCS_EX(13) &
              & + XSA(7)*ZTC0*ZTC0*ICE_T_PARAMETERS%XCS_EX(13) &
              & + XSA(8)*ZTC0*ICE_T_PARAMETERS%XCS_EX(13)*ICE_T_PARAMETERS%XCS_EX(13) + XSA(9)*ZTC0*ZTC0*ZTC0 &
              & + XSA(10)*ICE_T_PARAMETERS%XCS_EX(13)*ICE_T_PARAMETERS%XCS_EX(13)*ICE_T_PARAMETERS%XCS_EX(13)
        ZA_A = 10.0**ZLOGA_A
        ZB_B = XSB(1)+ XSB(2)*ZTC0 + XSB(3)*ICE_T_PARAMETERS%XCS_EX(13) + XSB(4)*ZTC0*ICE_T_PARAMETERS%XCS_EX(13) &
           & + XSB(5)*ZTC0*ZTC0 + XSB(6)*ICE_T_PARAMETERS%XCS_EX(13)*ICE_T_PARAMETERS%XCS_EX(13) &
           & + XSB(7)*ZTC0*ZTC0*ICE_T_PARAMETERS%XCS_EX(13) &
           & + XSB(8)*ZTC0*ICE_T_PARAMETERS%XCS_EX(13)*ICE_T_PARAMETERS%XCS_EX(13) &
           & + XSB(9)*ZTC0*ZTC0*ZTC0 &
           & + XSB(10)*ICE_T_PARAMETERS%XCS_EX(13)*ICE_T_PARAMETERS%XCS_EX(13)*ICE_T_PARAMETERS%XCS_EX(13)
        ZSMOE = ZA_A * ZSMO2**ZB_B

        ! Snow collecting cloud water.  In CE, assume XITDC<<Ds and vtc=~0.
        IF (PRST(JK)>0.0) ZX_DS(JK) = ZSMOC / ZSMOB
!       Add conditions for snow collecting cloud water
        IF (ZX_DS(JK) > XD0S .AND. PRCT(JK) > 0.0 .AND. PMVD_C(JK) > XD0C .AND. PRCS(JK) > 0.0) THEN
          IDX = 1 + INT(NBS*LOG(ZX_DS(JK) &
                    & / ICE_T_PARAMETERS%XITDS(1)) &
                    & / LOG(ICE_T_PARAMETERS%XITDS(NBS) &
                    & / ICE_T_PARAMETERS%XITDS(1)))
          IDX = MIN(IDX, NBS)
          ZEF_SW = ICE_T_PARAMETERS%XT_EFSW(IDX, INT(PMVD_C(JK)*1.E6))
          ZPRS_SCW = PRHOF(JK)*ICE_T_PARAMETERS%XT1_QS_QC*ZEF_SW*PRCT(JK)*ZSMOE

          ! A portion of rimed snow converts to graupel but some remains snow.
          ! Interp from 5 to 75% as riming factor increases from 5.0 to 30.0
          ! 0.028 came from (.75-.05)/(30.-5.).  This remains ad-hoc and should
          ! be revisited.
          IF (ZPRS_SCW .GT. 5.0*PPRS_SDE(JK) .AND. PPRS_SDE(JK) .GT. XICET_EPS) THEN
            ZR_FRAC = MIN(30.0D0, ZPRS_SCW/PPRS_SDE(JK))
            ZG_FRAC = MIN(0.75, 0.05 + (ZR_FRAC-5.)*.028)
            ZPRG_SCW = ZG_FRAC*ZPRS_SCW
            ZPRS_SCW = (1. - ZG_FRAC)*ZPRS_SCW
          ENDIF

          PRCS(JK) = PRCS(JK) - ZPRS_SCW - ZPRG_SCW
          PRSS(JK) = PRSS(JK) + ZPRS_SCW
          PRGS(JK) = PRGS(JK) + ZPRG_SCW
          PZTHS(JK) = PZTHS(JK) + (ZPRS_SCW + ZPRG_SCW) *(ZLSFACT(JK)-ZLVFACT(JK))
        ENDIF
      ENDDO

    ELSE  ! NOT ICE_T
!
      IGRIM=0
      DO JK=1, KSIZE
        IF((PRCT(JK)>ICED%XRTMIN(2)) .AND. (PRST(JK)>ICED%XRTMIN(5)) .AND. &
                                    (PRCS(JK)>0.0) .AND. (ZZT(JK)<CST%XTT)) THEN
          IGRIM=IGRIM+1
          GMASK(JK)=.TRUE.
          ! 5.1.1  select the ZLBDAS
          I1(IGRIM)=JK
          ZVEC1(IGRIM)=ZLBDAS(JK)
        ELSE
          GMASK(JK)=.FALSE.
        ENDIF
      ENDDO

      IF( IGRIM>0 ) THEN
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
        ZZW(:) = 0.
        DO JK=1, IGRIM
          ZZW(I1(JK))=ZVEC1(JK)
        ENDDO
!
!        5.1.4  riming of the small sized aggregates
!
        DO JK = 1, KSIZE
          IF (GMASK(JK)) THEN
            ZZW1(JK,1) = MIN( PRCS(JK),                                 &
                          ICEP%XCRIMSS * ZZW(JK) * PRCT(JK)*ZCOLF(JK)   & ! RCRIMSS
                                        *   ZLBDAS(JK)**ICEP%XEXCRIMSS &
                                        * ZRHODREF(JK)**(-ICED%XCEXVT) )
            PRCS(JK) = PRCS(JK) - ZZW1(JK,1)
            PRSS(JK) = PRSS(JK) + ZZW1(JK,1)
            PZTHS(JK) = PZTHS(JK) + ZZW1(JK,1)*(ZLSFACT(JK)-ZLVFACT(JK)) ! f(L_f*(RCRIMSS))
          END IF
        END DO
!
!        5.1.5  perform the linear interpolation of the normalized
!               "XBS"-moment of the incomplete gamma function
!
        ZVEC1(1:IGRIM) = ICEP%XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                      - ICEP%XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
        ZZW(:) = 0.
        DO JK=1, IGRIM
          ZZW(I1(JK))=ZVEC1(JK)
        ENDDO
!
!        5.1.6  riming-conversion of the large sized aggregates into graupeln
!
!
        DO JK = 1, KSIZE
          IF (GMASK(JK) .AND. (PRSS(JK) > 0.0)) THEN
            ZZW1(JK,2) = MIN(PRCS(JK),                                   &
                            ICEP%XCRIMSG * PRCT(JK)*ZCOLF(JK)            & ! RCRIMSG
                                        * ZLBDAS(JK)**ICEP%XEXCRIMSG   &
                                        * ZRHODREF(JK)**(-ICED%XCEXVT) &
                                        - ZZW1(JK,1))

            ZZW1(JK,3) = MIN(PRSS(JK),                                 &
                            ICEP%XSRIMCG * ZLBDAS(JK)**ICEP%XEXSRIMCG & ! RSRIMCG
                                        * (1.0 - ZZW(JK))/(PTSTEP*ZRHODREF(JK)))

            PRCS(JK) = PRCS(JK) - ZZW1(JK,2)
            PRSS(JK) = PRSS(JK) - ZZW1(JK,3)
            PRGS(JK) = PRGS(JK) + ZZW1(JK,2)+ZZW1(JK,3)
            PZTHS(JK) = PZTHS(JK) + ZZW1(JK,2)*(ZLSFACT(JK)-ZLVFACT(JK)) ! f(L_f*(RCRIMSG))
          END IF
        END DO
      END IF
    ENDIF !END ICE_T

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'RIM', ZWKBUD)
    IF (BUCONF%LBUDGET_RC) ZWKBUD(:,:) = UNPACK(PRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'RIM', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'RIM', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'RIM', ZWKBUD)

!*       5.2    rain accretion onto the aggregates

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'ACC', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'ACC', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'ACC', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'ACC', ZWKBUD)

    ZZW1(:,2:3) = 0.0
    IGACC=0
    DO JK=1, KSIZE
      IF((PRRT(JK) > ICED%XRTMIN(3)) .AND. &
         (PRST(JK) > ICED%XRTMIN(5)) .AND. &
         (PRRS(JK) > 0.0)            .AND. &
         (ZZT(JK) < CST%XTT)) THEN
        IGACC=IGACC+1
        GMASK(JK)=.TRUE.
        ! 5.1.1  select the (ZLBDAS,ZLBDAR) couplet
        I1(IGACC)=JK
        ZVEC1(IGACC)=ZLBDAS(JK)
        ZVEC2(IGACC)=ZLBDAR(JK)
      ELSE
        GMASK(JK)=.FALSE.
      ENDIF
    ENDDO

    IF (OICE_T) THEN
      PVTR = 0.0
      ZVTS = 0.0
      ! added stricter terms for rain accreting snow
      DO JK=1, KSIZE
        GMASK(JK) = GMASK(JK) .AND. &
                  & (PMVD_R(JK) > 50e-6) .AND. &
                  & (ZX_DS(JK) > 100e-6)
      ENDDO
    ENDIF

    IF( IGACC>0 ) THEN
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

      DO JL = 1,IGACC
        ZVEC3(JL) = ( ICEP%XKER_RACCSS(IVEC1(JL)+1,IVEC2(JL)+1)* ZVEC2(JL)         &
                    - ICEP%XKER_RACCSS(IVEC1(JL)+1,IVEC2(JL)  )*(ZVEC2(JL) - 1.0)) &
                                                               * ZVEC1(JL)         &
                  - ( ICEP%XKER_RACCSS(IVEC1(JL)  ,IVEC2(JL)+1)* ZVEC2(JL)         &
                    - ICEP%XKER_RACCSS(IVEC1(JL)  ,IVEC2(JL)  )*(ZVEC2(JL) - 1.0)) &
                                                               *(ZVEC1(JL) - 1.0)
      END DO
      ZZW(:) = 0.
      DO JK=1, IGACC
        ZZW(I1(JK))=ZVEC3(JK)
      ENDDO
!
!        5.2.4  raindrop accretion on the small sized aggregates
!
    IF (OICE_T) THEN
      ZRHO00 = CST%XP00/(CST%XRD*XTHVREFZ(1+JPVEXT))
      DO JK = 1, KSIZE
        ZEF_SR(JK) = 1.0
        ZVTS(JK) = ICED%XCS*ZX_DS(JK)**ICED%XDS*ZRHODREF(JK)**(-ICED%XCEXVT)
        PVTR(JK) = ICED%XCR*PMVD_R(JK)**ICED%XDR*ZRHODREF(JK)**(-ICED%XCEXVT)
        ZRATIO = 0.0

        IF (GMASK(JK) .AND. PVTR(JK) > 0.0000001) THEN
          ZRATIO = ZVTS(JK)/PVTR(JK)
        ENDIF

        IF (ZRATIO < 0.25) THEN
          ZEF_SR(JK)=0.9
        ELSEIF (ZRATIO > 1.75) THEN
          ZEF_SR(JK)=0.9
        ELSE
          ZEF_SR(JK) = (COS(((ZRATIO - 0.25)/0.75)*XPI)+1)*0.40 + 0.1
        ENDIF

        IF (GMASK(JK)) THEN
          ZFRACCSS_V = ((XPI**2)/24.0)*ICED%XCCS*PCCR_V(JK)*CST%XRHOLW*(ZRHO00**ICED%XCEXVT)
          ! coef of RRACCS
          ZZW1(JK,2) = ZFRACCSS_V*(ZLBDAS(JK)**ICED%XCXS)*(ZRHODREF(JK)**(-ICED%XCEXVT-1.)) &
                   & * (ICEP%XLBRACCS1/((ZLBDAS(JK)**2)                )               &
                   & +  ICEP%XLBRACCS2/( ZLBDAS(JK)    * ZLBDAR(JK)    )               &
                   & +  ICEP%XLBRACCS3/(                (ZLBDAR(JK)**2)))/ZLBDAR(JK)**4
                    !BJKE: added variable collection efficiency
          ZZW1(JK,4) = MIN(PRRS(JK),ZZW1(JK,2)*ZZW(JK)*ZEF_SR(JK))           ! RRACCSS
          PRRS(JK) = PRRS(JK) - ZZW1(JK,4)*ICEP%XFRMIN(7)
          PRSS(JK) = PRSS(JK) + ZZW1(JK,4)*ICEP%XFRMIN(7)
          PZTHS(JK) = PZTHS(JK) + ZZW1(JK,4)*(ZLSFACT(JK)-ZLVFACT(JK))*ICEP%XFRMIN(7) ! f(L_f*(RRACCSS))
        ENDIF
      ENDDO
    ELSE ! End ICE-T
      DO JK = 1, KSIZE
        IF (GMASK(JK)) THEN
          ZZW1(JK,2) =                                            & !! coef of RRACCS
                  ICEP%XFRACCSS*( ZLBDAS(JK)**ICED%XCXS )*( ZRHODREF(JK)**(-ICED%XCEXVT-1.) ) &
           *(ICEP%XLBRACCS1/((ZLBDAS(JK)**2)               ) +                  &
             ICEP%XLBRACCS2/( ZLBDAS(JK)   * ZLBDAR(JK)    ) +                  &
             ICEP%XLBRACCS3/(               (ZLBDAR(JK)**2)) )/ZLBDAR(JK)**4
          ZZW1(JK,4) = MIN( PRRS(JK),ZZW1(JK,2)*ZZW(JK) )           ! RRACCSS
          PRRS(JK) = PRRS(JK) - ZZW1(JK,4)*ICEP%XFRMIN(7)
          PRSS(JK) = PRSS(JK) + ZZW1(JK,4)*ICEP%XFRMIN(7)
          PZTHS(JK) = PZTHS(JK) + ZZW1(JK,4)*(ZLSFACT(JK)-ZLVFACT(JK))*ICEP%XFRMIN(7) ! f(L_f*(RRACCSS))
        END IF
      END DO
    ENDIF
!
!        5.2.4b perform the bilinear interpolation of the normalized
!               RACCS-kernel
!
      DO JL = 1,IGACC
        ZVEC3(JL) =  (ICEP%XKER_RACCS(IVEC2(JL)+1,IVEC1(JL)+1)* ZVEC1(JL)          &
                    - ICEP%XKER_RACCS(IVEC2(JL)+1,IVEC1(JL)  )*(ZVEC1(JL) - 1.0) ) &
                                                              * ZVEC2(JL)          &
                   - (ICEP%XKER_RACCS(IVEC2(JL)  ,IVEC1(JL)+1)* ZVEC1(JL)          &
                   -  ICEP%XKER_RACCS(IVEC2(JL)  ,IVEC1(JL)  )*(ZVEC1(JL) - 1.0) ) &
                                                              *(ZVEC2(JL) - 1.0)
      END DO
      ZZW(:) = 0.
      DO JK=1, IGACC
        ZZW(I1(JK)) = ZVEC3(JK)
      ENDDO
      ZZW1(:,2) = ZZW1(:,2)*ZZW(:)

                                                                       !! RRACCS!
!        5.2.5  perform the bilinear interpolation of the normalized
!               SACCRG-kernel
!
      DO JL = 1,IGACC
        ZVEC3(JL) =  (  ICEP%XKER_SACCRG(IVEC2(JL)+1,IVEC1(JL)+1)* ZVEC1(JL)          &
                      - ICEP%XKER_SACCRG(IVEC2(JL)+1,IVEC1(JL)  )*(ZVEC1(JL) - 1.0) ) &
                                                            * ZVEC2(JL) &
                   - (  ICEP%XKER_SACCRG(IVEC2(JL)  ,IVEC1(JL)+1)* ZVEC1(JL)          &
                      - ICEP%XKER_SACCRG(IVEC2(JL)  ,IVEC1(JL)  )*(ZVEC1(JL) - 1.0) ) &
                                                          * (ZVEC2(JL) - 1.0)
      END DO
      ZZW(:) = 0.
      DO JK=1, IGACC
        ZZW(I1(JK)) = ZVEC3(JK)
      ENDDO
!
!        5.2.6  raindrop accretion-conversion of the large sized aggregates
!               into graupeln
!
      DO JK = 1, KSIZE
        IF (GMASK(JK) .AND. (PRSS(JK) > 0.0)) THEN
          ZZW1(JK,2) = MAX( MIN( PRRS(JK),ZZW1(JK,2)-ZZW1(JK,4) ),0.0 )       ! RRACCSG
        END IF
      END DO

      DO JK = 1, KSIZE
        IF (GMASK(JK) .AND. (PRSS(JK)>0.0) .AND. ZZW1(JK,2) > 0.0 .AND. PRSS(JK) > ICEP%XFRMIN(1)/PTSTEP) THEN

          IF (OICE_T) THEN
            ZFSACCRG = (XPI/4.0)*ICED%XAS*ICED%XCCS*PCCR_V(JK)*(ZRHO00**ICED%XCEXVT)
          ELSE
            ZFSACCRG = ICEP%XFSACCRG
          ENDIF

          ZZW1(JK,3) = MIN( PRSS(JK), ZFSACCRG*ZZW(JK)*                     & ! RSACCRG
                ( ZLBDAS(JK)**(ICED%XCXS-ICED%XBS) )*( ZRHODREF(JK)**(-ICED%XCEXVT-1.) ) &
               *( ICEP%XLBSACCR1/((ZLBDAR(JK)**2)               ) +           &
                  ICEP%XLBSACCR2/( ZLBDAR(JK)    * ZLBDAS(JK)    ) +           &
                  ICEP%XLBSACCR3/(               (ZLBDAS(JK)**2)) )/ZLBDAR(JK) )
          PRRS(JK) = PRRS(JK) - ZZW1(JK,2)
          PRSS(JK) = PRSS(JK) - ZZW1(JK,3)
          PRGS(JK) = PRGS(JK) + ZZW1(JK,2)+ZZW1(JK,3)
          PZTHS(JK) = PZTHS(JK) + ZZW1(JK,2)*(ZLSFACT(JK)-ZLVFACT(JK)) !
                                 ! f(L_f*(RRACCSG))
        END IF
      END DO
    END IF

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'ACC', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'ACC', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'ACC', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'ACC', ZWKBUD)

!*       5.3    Conversion-Melting of the aggregates

    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'CMEL', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'CMEL', ZWKBUD)

    ZZW(:) = 0.0
    DO JK = 1, KSIZE
      IF ((PRST(JK)>ICED%XRTMIN(5)) .AND. (PRSS(JK)>0.0) .AND. (ZZT(JK)>CST%XTT)) THEN
        ZZW(JK) = PRVT(JK)*ZPRES(JK)/(CST%XEPSILO+PRVT(JK)) ! Vapor pressure
        ZZW(JK) =  ZKA(JK)*(CST%XTT-ZZT(JK)) +                                 &
                 ( ZDV(JK)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZZT(JK) - CST%XTT )) &
                             *(CST%XESTT-ZZW(JK))/(CST%XRV*ZZT(JK))             )
!
! compute RSMLT
!
        ZZW(JK)  = MIN( PRSS(JK), ICEP%XFSCVMG*MAX( 0.0,( -ZZW(JK) *             &
                             ( ICEP%X0DEPS*       ZLBDAS(JK)**ICEP%XEX0DEPS + &
                               ICEP%X1DEPS*ZCJ(JK)*ZLBDAS(JK)**ICEP%XEX1DEPS ) -   &
                                       ( ZZW1(JK,1)+ZZW1(JK,4) ) *       &
                                (ZRHODREF(JK)*CST%XCL*(CST%XTT-ZZT(JK)))) /    &
                                               ( ZRHODREF(JK)*CST%XLMTT ) ) )
!
! note that RSCVMG = RSMLT*ICEP%XFSCVMG but no heat is exchanged (at the rate RSMLT)
! because the graupeln produced by this process are still icy!!!
!
        PRSS(JK) = PRSS(JK) - ZZW(JK)
        PRGS(JK) = PRGS(JK) + ZZW(JK)
      END IF
    END DO

    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'CMEL', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'CMEL', ZWKBUD)

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RS',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_FAST_RS

END MODULE MODE_RAIN_ICE_OLD_FAST_RS
