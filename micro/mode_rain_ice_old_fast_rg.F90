!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_FAST_RG

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_FAST_RG(D, CST, ICEP, ICED, ICE_T_PARAMETERS, BUCONF, &
                                & PTSTEP, KSIZE, KRR,                           &
                                & OCND2, OICE_T, LTIW, GMICRO,                  &
                                & PRHODJ, PTHS,                                 &
                                & PRVT, PRCT, PRIT, PRRT, PRST, PRGT, PCIT,     &
                                & PRIS, PRRS, PRCS, PRSS, PRGS, PRHS, PZTHS,    &
                                & PRHODREF, PZRHODJ, PLSFACT, PLVFACT,          &
                                & PCJ, PKA, PDV,                                &
                                & PLBDAR, PLBDAG, PLBDAS,                       &
                                & PTIW, PZT, PRES,                              &
                                & PMVD_C, PMVD_R, PRHOF, PVTR, PCCR_V,          &
                                & TBUDGETS, KBUDGETS)

    USE MODD_PRECISION,        ONLY: MNHREAL64
    USE MODD_PARAMETERS,       ONLY: JPVEXT
    USE YOMHOOK,               ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_T
    USE MODD_CST,              ONLY: CST_T, XRHOLW, XPI
    USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_T
    USE MODD_ICET_PARAM,       ONLY: ICET_PARAM_t
    USE MODD_ICET_PARAM,       ONLY: XGONV_MIN, XGONV_MAX, XBM_G, XMU_G, XR_G, XD0C, XD0G, XTHVREFZ

    USE MODD_BUDGET,           ONLY: TBUDGETDATA_PTR, TBUDGETCONF_t, NBUDGET_TH, NBUDGET_RG, NBUDGET_RR, NBUDGET_RC, &
                                     NBUDGET_RI, NBUDGET_RS, NBUDGET_RH

    IMPLICIT NONE

    TYPE(DIMPHYEX_T),       INTENT(IN) :: D
    TYPE(CST_T),            INTENT(IN) :: CST
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED
    TYPE(ICET_PARAM_t),     INTENT(IN) :: ICE_T_PARAMETERS
    TYPE(TBUDGETCONF_t),    INTENT(IN) :: BUCONF

    REAL,    INTENT(IN) :: PTSTEP  ! Double Time step
    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, INTENT(IN) :: KRR

    LOGICAL, INTENT(IN) :: OCND2
    LOGICAL, INTENT(IN) :: OICE_T
    LOGICAL, INTENT(IN) :: LTIW

    LOGICAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: GMICRO ! Layer thickness (m)

    REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian

    REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PTHS    ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRRT     ! Rain water m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRGT     ! Graupel m.r. at t
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PCIT     ! Pristine ice conc. at t

    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRIS     ! Pristine ice m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRRS     ! Rain water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRCS     ! Cloud water m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRSS     ! Snow/aggregate m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRGS     ! Graupel m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PRHS     ! Hail m.r. source
    REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PZTHS     ! Theta source

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRHODREF ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PZRHODJ   ! RHO times Jacobian
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PKA      ! Thermal conductivity of the air
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PTIW     ! Wet bulb temperature

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PZT      ! Temperature
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRES    ! Pressure

    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PMVD_C
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PMVD_R
    REAL(KIND=MNHREAL64), DIMENSION(KSIZE), INTENT(IN) :: PRHOF
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PVTR
    REAL, DIMENSION(KSIZE), INTENT(IN)    :: PCCR_V

    TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
    INTEGER, INTENT(IN) :: KBUDGETS

    LOGICAL, DIMENSION(KSIZE) :: GDRY ! Test where to compute dry growth

    INTEGER, DIMENSION(KSIZE) :: IVEC1 ! Vectors of indices for
    INTEGER, DIMENSION(KSIZE) :: IVEC2 ! Vectors of indices for

    REAL, DIMENSION(KSIZE) :: ZVEC1 ! Work vectors for interpolations
    REAL, DIMENSION(KSIZE) :: ZVEC2 ! Work vectors for interpolations
    REAL, DIMENSION(KSIZE) :: ZVEC3 ! Work vectors for interpolations

    REAL, DIMENSION(KSIZE) :: ZRDRYG   ! Dry growth rate of the graupeln
    REAL, DIMENSION(KSIZE) :: ZRWETG   ! Wet growth rate of the graupeln

    REAL, DIMENSION(KSIZE) :: ZUSW  ! Undersaturation over water

    REAL :: ZEF_GW, ZEF_GR, ZTEMPC, ZVISCO
    REAL :: ZSLW1, ZYGRA1, ZZANS1, ZX_DG, ZSTOKE_G, ZFRDRYG_V
    REAL(KIND=MNHREAL64) :: ZN0_MIN, ZN0_EXP, ZLAMG, ZLAM_EXP, ZILAMG, ZN0_G
    REAL, DIMENSION(KSIZE) :: ZVTG
    REAL :: ZRHO00

    REAL, DIMENSION(KSIZE)      :: ZZW  ! Work array
    REAL, DIMENSION(KSIZE, KRR) :: ZZW1 ! Work array
    REAL, DIMENSION(D%NIJT,D%NKT) :: ZWKBUD

    INTEGER :: IGDRY
    INTEGER, DIMENSION(KSIZE) :: I1
    INTEGER :: JL, JK

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       6.1    rain contact freezing
!
    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RG',0,ZHOOK_HANDLE)

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'CFRZ', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'CFRZ', ZWKBUD)
    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'CFRZ', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'CFRZ', ZWKBUD)

    ZZW1(:,3:4) = 0.0
    DO JK = 1, KSIZE
      IF ((PRIT(JK) > ICED%XRTMIN(4) .AND. PRIT(JK)>ICEP%XFRMIN(2)) .AND. &
          (PRRT(JK) > ICED%XRTMIN(3)) .AND. &
          (PRIS(JK) > 0.0) .AND. &
          (PRRS(JK) > 0.0)) THEN
        ZZW1(JK,3) = MIN( PRIS(JK),ICEP%XICFRR * PRIT(JK)                & ! RICFRRG
                                        * PLBDAR(JK)**ICEP%XEXICFRR    &
                                        * PRHODREF(JK)**(-ICED%XCEXVT) )
        ZZW1(JK,4) = MIN( PRRS(JK),ICEP%XRCFRI * PCIT(JK)                & ! RRCFRIG
                                        * PLBDAR(JK)**ICEP%XEXRCFRI    &
                                        * PRHODREF(JK)**(-ICED%XCEXVT-1.) )
        PRIS(JK) = PRIS(JK) - ZZW1(JK,3)
        PRRS(JK) = PRRS(JK) - ZZW1(JK,4)
        PRGS(JK) = PRGS(JK) + ZZW1(JK,3)+ZZW1(JK,4)
        PZTHS(JK) = PZTHS(JK) + ZZW1(JK,4)*(PLSFACT(JK)-PLVFACT(JK)) ! f(L_f*RRCFRIG)
      END IF
    END DO

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'CFRZ', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'CFRZ', ZWKBUD)
    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'CFRZ', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'CFRZ', ZWKBUD)

!*       6.2    compute the Dry growth case

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'WETG', ZWKBUD)
    IF (BUCONF%LBUDGET_RC) ZWKBUD(:,:) = UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'WETG', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'WETG', ZWKBUD)
    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'WETG', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'WETG', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'WETG', ZWKBUD)
    IF ( KRR == 7 ) THEN
      IF (BUCONF%LBUDGET_RH) ZWKBUD(:,:) = UNPACK(PRHS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
      IF (BUCONF%LBUDGET_RH) CALL TBUDGETS(NBUDGET_RH)%PTR%INIT_PHY(D, 'WETG', ZWKBUD)
    END IF

    ZZW1(:,:) = 0.0
    IF (OICE_T) THEN
      ! replaced the graupel collecting cloud water formulation with Thompson et al. (2008)
      ZVTG(:) = 0.0
      ZN0_MIN = XGONV_MAX
      DO JK = 1, KSIZE

        IF ((PRGT(JK) > ICED%XRTMIN(6))) THEN

          ZTEMPC = PZT(JK) - CST%XTT
          IF (ZTEMPC > 0.0) THEN
            ZVISCO = (1.718 + 0.0049*ZTEMPC)*1.0E-5
          ELSE
            ZVISCO = (1.718 + 0.0049*ZTEMPC - 1.2E-5*ZTEMPC*ZTEMPC)*1.0E-5
          ENDIF
          IF (PZT(JK) < 270.65 .AND. PRRT(JK) > 0.0 .AND. PMVD_R(JK) > 100.E-6) THEN
            ZSLW1 = 4.01 + ALOG10(PMVD_R(JK))
          ELSE
            ZSLW1 = 0.01
          ENDIF

          ZYGRA1 = 4.31 + ALOG10(max(5.E-5, PRGT(JK)))
          ZZANS1 = 3.1 + (100./(300.*ZSLW1*ZYGRA1/(10./ZSLW1 + 1. + 0.25*ZYGRA1)+30. + 10.*ZYGRA1))
          ZN0_EXP = 10.**(ZZANS1)
          ZN0_EXP = MAX(DBLE(XGONV_MIN), MIN(ZN0_EXP, DBLE(XGONV_MAX)))
          ZN0_MIN = MIN(ZN0_EXP, ZN0_MIN)
          ZN0_EXP = ZN0_MIN
          ZLAM_EXP = (ZN0_EXP*ICE_T_PARAMETERS%XAM_G*ICE_T_PARAMETERS%XCG_G(1)/PRGT(JK))**ICE_T_PARAMETERS%XOGE1
          ZLAMG = ZLAM_EXP &
              & * (ICE_T_PARAMETERS%XCG_G(3)*ICE_T_PARAMETERS%XOGG2*ICE_T_PARAMETERS%XOGG1)**ICE_T_PARAMETERS%XOBMG
          ZILAMG = 1./ZLAMG
          ZN0_G = ZN0_EXP/(ICE_T_PARAMETERS%XCG_G(2)*ZLAM_EXP) * ZLAMG**ICE_T_PARAMETERS%XCG_E(2)
          ZX_DG = (XBM_G + XMU_G + 1.) * ZILAMG
          ZVTG(JK) = ICED%XCG*ZX_DG**ICED%XDG*PRHODREF(JK)**(-ICED%XCEXVT)

          ! Graupel collecting cloud water.  In CE, assume XITDC<<Dg and vtc=~0.
          IF (PRGT(JK) > XR_G(1) .AND. PMVD_C(JK) > XD0C .AND. PRCS(JK) > 0.0) THEN
            ZSTOKE_G = PMVD_C(JK)*PMVD_C(JK)*ZVTG(JK)*XRHOLW/(9.*ZVISCO*ZX_DG)
            IF (ZX_DG > XD0G) THEN
                IF (ZSTOKE_G>0.4 .AND. ZSTOKE_G<10.) THEN
                  ZEF_GW = 0.55*ALOG10(2.51*ZSTOKE_G)
                ELSEIF (ZSTOKE_G<0.4) THEN
                  ZEF_GW = 0.0
                ELSEIF (ZSTOKE_G>10) THEN
                  ZEF_GW = 0.77
                ENDIF
                ZZW1(JK,1) = PRHOF(JK)*ICE_T_PARAMETERS%XT1_QG_QC*ZEF_GW*PRCT(JK)*ZN0_G*ZILAMG**ICE_T_PARAMETERS%XCG_E(9)
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      !Graupel collecting cloud ice is set to zero
      ZZW1(:,2)=0.0
    ELSE
      DO JK = 1, KSIZE
        IF ((PRGT(JK) > ICED%XRTMIN(6)) .AND. &
           ((PRCT(JK) > ICED%XRTMIN(2)  .AND. PRCS(JK) > 0.0))) THEN
          ZZW(JK) = PLBDAG(JK)**(ICED%XCXG-ICED%XDG-2.0) * PRHODREF(JK)**(-ICED%XCEXVT)
          ZZW1(JK,1) = MIN( PRCS(JK),ICEP%XFCDRYG * PRCT(JK) * ZZW(JK) )             ! RCDRYG
        END IF
      END DO

      DO JK = 1, KSIZE
        IF ((PRGT(JK) > ICED%XRTMIN(6)) .AND. &
           ((PRIT(JK) > ICED%XRTMIN(4)  .AND. PRIS(JK)>0.0))) THEN

          ZZW(JK) = PLBDAG(JK)**(ICED%XCXG-ICED%XDG-2.0) * PRHODREF(JK)**(-ICED%XCEXVT)
          ZZW1(JK,2) = MIN(PRIS(JK), ICEP%XFIDRYG * EXP(ICEP%XCOLEXIG*(PZT(JK) - CST%XTT)) &
                                                  * PRIT(JK) * ZZW(JK) )             ! RIDRYG
        END IF
      END DO
    ENDIF
!
!*       6.2.1  accretion of aggregates on the graupeln
!
    IGDRY=0
    DO JK=1, KSIZE
      IF ((PRST(JK) > ICED%XRTMIN(5)) .AND. &
          (PRGT(JK) > ICED%XRTMIN(6)) .AND. &
          (PRSS(JK) > 0.0)            .AND. &
          (.NOT. OICE_T .OR. ((ZVTG(JK) > 50E-5) .AND. & ! added stricter conditions for ICE-T
                              (PVTR(JK) > 50E-5)))) THEN   
        IGDRY=IGDRY + 1
        GDRY(JK)=.TRUE.

        ! 6.2.3 select the (PLBDAG,PLBDAS) couplet
        I1(IGDRY)=JK
        ZVEC1(IGDRY)=PLBDAG(JK)
        ZVEC2(IGDRY)=PLBDAS(JK)
      ELSE
        GDRY(JK)=.FALSE.
      ENDIF
    ENDDO

    IF( IGDRY>0 ) THEN
!
!*       6.2.4  find the next lower indice for the PLBDAG and for the PLBDAS
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
      DO JL = 1,IGDRY
        ZVEC3(JL) =  (  ICEP%XKER_SDRYG(IVEC1(JL)+1,IVEC2(JL)+1)* ZVEC2(JL)          &
                      - ICEP%XKER_SDRYG(IVEC1(JL)+1,IVEC2(JL)  )*(ZVEC2(JL) - 1.0) ) &
                                                           * ZVEC1(JL) &
                   - (  ICEP%XKER_SDRYG(IVEC1(JL)  ,IVEC2(JL)+1)* ZVEC2(JL)          &
                      - ICEP%XKER_SDRYG(IVEC1(JL)  ,IVEC2(JL)  )*(ZVEC2(JL) - 1.0) ) &
                                                           * (ZVEC1(JL) - 1.0)
      END DO
      ZZW(:) = 0.
      DO JK=1, IGDRY
        ZZW(I1(JK))=ZVEC3(JK)
      ENDDO
!
      IF (OCND2) THEN
        ZZW1(:,3) = 0.
      ELSE
        DO JK = 1, KSIZE
          IF (GDRY(JK)) THEN
            ZZW1(JK,3) = MIN(PRSS(JK), ICEP%XFSDRYG*ZZW(JK)            & ! RSDRYG
                               * EXP(ICEP%XCOLEXSG*(PZT(JK)-CST%XTT))  &
                               *(PLBDAS(JK)**(ICED%XCXS-ICED%XBS))     &
                               *(PLBDAG(JK)**ICED%XCXG)                &
                               *(PRHODREF(JK)**(-ICED%XCEXVT-1.))      &
                               *(ICEP%XLBSDRYG1/( PLBDAG(JK)**2               ) + &
                                 ICEP%XLBSDRYG2/( PLBDAG(JK)   * PLBDAS(JK)   ) + &
                                 ICEP%XLBSDRYG3/(                PLBDAS(JK)**2) ) )
          END IF
        END DO
      ENDIF
    END IF
!
!*       6.2.6  accretion of raindrops on the graupeln
!
    IGDRY=0
    DO JK=1, KSIZE
      IF((PRRT(JK)>ICED%XRTMIN(3)) .AND. (PRGT(JK)>ICED%XRTMIN(6)) .AND. (PRRS(JK)>0.0)) THEN
        IGDRY=IGDRY+1
        GDRY(JK)=.TRUE.
        ! 6.2.8  select the (PLBDAG,PLBDAR) couplet
        I1(IGDRY)=JK
        ZVEC1(IGDRY)=PLBDAG(JK)
        ZVEC2(IGDRY)=PLBDAR(JK)
      ELSE
        GDRY(JK)=.FALSE.
      ENDIF
    ENDDO

    IF (IGDRY>0) THEN
!
!*       6.2.9  find the next lower indice for the PLBDAG and for the PLBDAR
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
      DO JL = 1,IGDRY
        ZVEC3(JL) =  (  ICEP%XKER_RDRYG(IVEC1(JL)+1,IVEC2(JL)+1)* ZVEC2(JL)         &
                      - ICEP%XKER_RDRYG(IVEC1(JL)+1,IVEC2(JL)  )*(ZVEC2(JL) - 1.0)) &
                                                                * ZVEC1(JL)         &
                   - (  ICEP%XKER_RDRYG(IVEC1(JL)  ,IVEC2(JL)+1)* ZVEC2(JL)         &
                      - ICEP%XKER_RDRYG(IVEC1(JL)  ,IVEC2(JL)  )*(ZVEC2(JL) - 1.0)) &
                                                                *(ZVEC1(JL) - 1.0)
      END DO
      ZZW(:) = 0.
      DO JK=1, IGDRY
        ZZW(I1(JK))=ZVEC3(JK)
      ENDDO

      IF (OICE_T) THEN
        ZRHO00 = CST%XP00/(CST%XRD*XTHVREFZ(1+JPVEXT))
        DO JK = 1, KSIZE
          IF (GDRY(JK)) THEN

            IF (ZVTG(JK)/PVTR(JK) < 0.25) THEN
              ZEF_GR = 0.9
            ELSEIF (ZVTG(JK)/PVTR(JK) > 1.75) THEN
              ZEF_GR = 0.9
            ELSE
              ZEF_GR = (COS(((ZVTG(JK)/PVTR(JK) - 0.25)/0.75)*XPI) + 1)*0.40 + 0.1
            ENDIF

            !BJKE: variable rain intercept parameter
            ZFRDRYG_V = ((XPI**2)/24.0)*ICED%XCCG*PCCR_V(JK)*XRHOLW*(ZRHO00**ICED%XCEXVT)

            !BJKE: added variable collection efficiency
            ZZW1(JK,4) = MIN(PRRS(JK), ZFRDRYG_V*ZZW(JK)*ZEF_GR          & ! RRDRYG
                     & *  (PLBDAR(JK)**(-4))*(PLBDAG(JK)**ICED%XCXG)     &
                     & *  (PRHODREF(JK)**(-ICED%XCEXVT-1.))              &
                     & *  (ICEP%XLBRDRYG1/(PLBDAG(JK)**2)                &
                     & +   ICEP%XLBRDRYG2/(PLBDAG(JK) * PLBDAR(JK))      &
                     & +   ICEP%XLBRDRYG3/(             PLBDAR(JK)**2)))
          ENDIF
        ENDDO
      ELSE
        DO JK = 1, KSIZE
          IF (GDRY(JK)) THEN
            ZZW1(JK,4) = MIN(PRRS(JK),ICEP%XFRDRYG*ZZW(JK)                     & ! RRDRYG
                             *(PLBDAR(JK)**(-4) )*(PLBDAG(JK)**ICED%XCXG)      &
                             *(PRHODREF(JK)**(-ICED%XCEXVT-1.))                &
                             *(ICEP%XLBRDRYG1/(PLBDAG(JK)**2               ) + &
                               ICEP%XLBRDRYG2/(PLBDAG(JK)   * PLBDAR(JK)   ) + &
                               ICEP%XLBRDRYG3/(               PLBDAR(JK)**2) ) )
          END IF
        END DO
      ENDIF
    END IF

    ZRDRYG(:) = ZZW1(:,1) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,4)
!
!*       6.3    compute the Wet growth case
!
    ZZW(:) = 0.0
    ZRWETG(:) = 0.0

    DO JK = 1, KSIZE
      IF (PRGT(JK)>ICED%XRTMIN(6)) THEN
        ZZW1(JK,5) = MIN( PRIS(JK),                                    &
                    ZZW1(JK,2) / (ICEP%XCOLIG*EXP(ICEP%XCOLEXIG*(PZT(JK)-CST%XTT)) ) ) ! RIWETG
        ZZW1(JK,6) = MIN( PRSS(JK),                                    &
                    ZZW1(JK,3) / (ICEP%XCOLSG*EXP(ICEP%XCOLEXSG*(PZT(JK)-CST%XTT)) ) ) ! RSWETG
!
        ZZW(JK) = PRVT(JK)*PRES(JK)/(CST%XEPSILO+PRVT(JK)) ! Vapor pressure
        ZZW(JK) = PKA(JK)*(CST%XTT-PZT(JK)) +                              &
                (PDV(JK)*(CST%XLVTT + (CST%XCPV - CST%XCL) * (PZT(JK) - CST%XTT)) &
                       *(CST%XESTT-ZZW(JK))/(CST%XRV*PZT(JK))           )
!
! compute RWETG
!
        ZRWETG(JK)=MAX(0.0,                                               &
                     (ZZW(JK) * ( ICEP%X0DEPG*       PLBDAG(JK)**ICEP%XEX0DEPG  + &
                                 ICEP%X1DEPG*PCJ(JK)*PLBDAG(JK)**ICEP%XEX1DEPG) + &
                     (ZZW1(JK,5)+ZZW1(JK,6)) *                            &
                     (PRHODREF(JK)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-PZT(JK))))) / &
                               (PRHODREF(JK)*(CST%XLMTT-CST%XCL*(CST%XTT-PZT(JK)))))
      END IF
    END DO
!
!*       6.4    Select Wet or Dry case
!
    ZZW(:) = 0.0
    IF (KRR == 7) THEN
      DO JK = 1, KSIZE
        IF (PRGT(JK) > ICED%XRTMIN(6) .AND. &
            PZT(JK) < CST%XTT .AND.         &
            ZRDRYG(JK) >= ZRWETG(JK) .AND.  &
            ZRWETG(JK) > 0.0) THEN

          ZZW(JK) = ZRWETG(JK) - ZZW1(JK,5) - ZZW1(JK,6) ! RCWETG+RRWETG
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
          ZZW1(JK,7) = MAX( 0.0,MIN( ZZW(JK),PRRS(JK)+ZZW1(JK,1) ) )
          ZUSW(JK)   = ZZW1(JK,7) / ZZW(JK)
          ZZW1(JK,5) = ZZW1(JK,5)*ZUSW(JK)
          ZZW1(JK,6) = ZZW1(JK,6)*ZUSW(JK)
          ZRWETG(JK) = ZZW1(JK,7) + ZZW1(JK,5) + ZZW1(JK,6)

          PRCS(JK) = PRCS(JK) - ZZW1(JK,1)
          PRIS(JK) = PRIS(JK) - ZZW1(JK,5)
          PRSS(JK) = PRSS(JK) - ZZW1(JK,6)
!
! assume a linear percent of conversion of graupel into hail
!
          PRGS(JK) = PRGS(JK) + ZRWETG(JK)                     !     Wet growth
          ZZW(JK)  = PRGS(JK)*ZRDRYG(JK)/(ZRWETG(JK)+ZRDRYG(JK)) !        and
          PRGS(JK) = PRGS(JK) - ZZW(JK)                        !   partial conversion
          PRHS(JK) = PRHS(JK) + ZZW(JK)                        ! of the graupel into hail
!
          PRRS(JK) = MAX( 0.0,PRRS(JK) - ZZW1(JK,7) + ZZW1(JK,1) )
          PZTHS(JK) = PZTHS(JK) + ZZW1(JK,7)*(PLSFACT(JK)-PLVFACT(JK)) ! f(L_f*(RCWETG+RRWETG)
        END IF
      END DO

    ELSE IF( KRR == 6 ) THEN
      DO JK = 1, KSIZE
        IF (PRGT(JK) > ICED%XRTMIN(6) .AND. &
            PRGT(JK) > ICEP%XFRMIN(3) .AND. &
            PRIS(JK)*PTSTEP > ICEP%XFRMIN(3) .AND. &
            PZT(JK) < CST%XTT .AND. &
            ZRDRYG(JK) >= ZRWETG(JK) .AND. &
            ZRWETG(JK) > 0.0) THEN
          ZZW(JK)  = ZRWETG(JK)
          PRCS(JK) = PRCS(JK) - ZZW1(JK,1)
          PRIS(JK) = PRIS(JK) - ZZW1(JK,5)
          PRSS(JK) = PRSS(JK) - ZZW1(JK,6)
          PRGS(JK) = PRGS(JK) + ZZW(JK)

          PRRS(JK) = PRRS(JK) - ZZW(JK) + ZZW1(JK,5) + ZZW1(JK,6) + ZZW1(JK,1)
          PZTHS(JK) = PZTHS(JK) + (ZZW(JK)-ZZW1(JK,5)-ZZW1(JK,6))*(PLSFACT(JK)-PLVFACT(JK))
                                                 ! f(L_f*(RCWETG+RRWETG))
        END IF
      END DO
    END IF

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'WETG', ZWKBUD)
    IF (BUCONF%LBUDGET_RC) ZWKBUD(:,:) = UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'WETG', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'WETG', ZWKBUD)
    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'WETG', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'WETG', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'WETG', ZWKBUD)
    IF ( KRR == 7 ) THEN
      IF (BUCONF%LBUDGET_RH) ZWKBUD(:,:) = UNPACK(PRHS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
      IF (BUCONF%LBUDGET_RH) CALL TBUDGETS(NBUDGET_RH)%PTR%END_PHY(D, 'WETG', ZWKBUD)
    END IF

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'DRYG', ZWKBUD)
    IF (BUCONF%LBUDGET_RC) ZWKBUD(:,:) = UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'DRYG', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'DRYG', ZWKBUD)
    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'DRYG', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'DRYG', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'DRYG', ZWKBUD)

    DO JK = 1, KSIZE
      IF (PRGT(JK) > ICED%XRTMIN(6) .AND. &
          PRGT(JK) > ICEP%XFRMIN(4) .AND. &
          PRIS(JK)*PTSTEP > ICEP%XFRMIN(4) .AND. &
          PZT(JK) < CST%XTT .AND. &
          ZRDRYG(JK) < ZRWETG(JK) .AND. &
          ZRDRYG(JK) > 0.0) THEN
        PRCS(JK) = PRCS(JK) - ZZW1(JK,1)
        PRIS(JK) = PRIS(JK) - ZZW1(JK,2)
        PRSS(JK) = PRSS(JK) - ZZW1(JK,3)
        PRRS(JK) = PRRS(JK) - ZZW1(JK,4)
        PRGS(JK) = PRGS(JK) + ZRDRYG(JK)
        PZTHS(JK) = PZTHS(JK) + (ZZW1(JK,1)+ZZW1(JK,4))*(PLSFACT(JK)-PLVFACT(JK))
                          ! f(L_f*(RCDRYG+RRDRYG))
      END IF
    END DO

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'DRYG', ZWKBUD)
    IF (BUCONF%LBUDGET_RC) ZWKBUD(:,:) = UNPACK(PRCS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'DRYG', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'DRYG', ZWKBUD)
    IF (BUCONF%LBUDGET_RI) ZWKBUD(:,:) = UNPACK(PRIS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'DRYG', ZWKBUD)
    IF (BUCONF%LBUDGET_RS) ZWKBUD(:,:) = UNPACK(PRSS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'DRYG', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'DRYG', ZWKBUD)

!*       6.5    Melting of the graupeln

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'GMLT', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'GMLT', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'GMLT', ZWKBUD)

    ZZW(:) = 0.0
    IF (LTIW) THEN

      DO JK = 1, KSIZE
        IF ((PRGT(JK) > ICED%XRTMIN(6)) .AND. &
            (PRGS(JK) > 0.0) .AND. &
            (PTIW(JK) > CST%XTT)) THEN
          ZZW(JK) = PRVT(JK)*PRES(JK)/(CST%XEPSILO+PRVT(JK)) ! Vapor pressure
          ZZW(JK) = PKA(JK)*(CST%XTT-PTIW(JK)) +                                 &
                  (PDV(JK)*(CST%XLVTT + (CST%XCPV - CST%XCL) * (PTIW(JK) - CST%XTT)) &
                         *(CST%XESTT-ZZW(JK))/(CST%XRV*PTIW(JK))             )
!
! compute RGMLTR
!
          ZZW(JK)  = ICEP%XFRMIN(8)*MIN( PRGS(JK), MAX( 0.0,( -ZZW(JK) *           &
                                 ( ICEP%X0DEPG*       PLBDAG(JK)**ICEP%XEX0DEPG + &
                                   ICEP%X1DEPG*PCJ(JK)*PLBDAG(JK)**ICEP%XEX1DEPG ) - &
                                           ( ZZW1(JK,1)+ZZW1(JK,4) ) *       &
                                    ( PRHODREF(JK)*CST%XCL*(CST%XTT-PTIW(JK))) ) /   &
                                                   ( PRHODREF(JK)*CST%XLMTT)))


          PRRS(JK) = PRRS(JK) + ZZW(JK)
          PRGS(JK) = PRGS(JK) - ZZW(JK)
          PZTHS(JK) = PZTHS(JK) - ZZW(JK)*(PLSFACT(JK)-PLVFACT(JK)) ! f(L_f*(-RGMLTR))
        END IF
      END DO
    ELSE

      DO JK = 1, KSIZE
        IF ((PRGT(JK)>ICED%XRTMIN(6)) .AND. &
            (PRGS(JK)>0.0) .AND. &
            (PZT(JK)>CST%XTT)) THEN
          ZZW(JK) = PRVT(JK)*PRES(JK)/(CST%XEPSILO+PRVT(JK)) ! Vapor pressure
          ZZW(JK) = PKA(JK)*(CST%XTT-PZT(JK)) +                                 &
                  (PDV(JK)*(CST%XLVTT + (CST%XCPV - CST%XCL) * (PZT(JK) - CST%XTT)) &
                         *(CST%XESTT-ZZW(JK))/(CST%XRV*PZT(JK)))
!
! compute RGMLTR
!
          ZZW(JK)  = ICEP%XFRMIN(8)*MIN(PRGS(JK), MAX(0.0, (-ZZW(JK) *         &
                             (ICEP%X0DEPG*       PLBDAG(JK)**ICEP%XEX0DEPG +   &
                              ICEP%X1DEPG*PCJ(JK)*PLBDAG(JK)**ICEP%XEX1DEPG) - &
                             (ZZW1(JK,1)+ZZW1(JK,4)) *                         &
                                (PRHODREF(JK)*CST%XCL*(CST%XTT-PZT(JK)))) /    &
                                               (PRHODREF(JK)*CST%XLMTT)))
          PRRS(JK) = PRRS(JK) + ZZW(JK)
          PRGS(JK) = PRGS(JK) - ZZW(JK)
          PZTHS(JK) = PZTHS(JK) - ZZW(JK)*(PLSFACT(JK)-PLVFACT(JK)) ! f(L_f*(-RGMLTR))
        END IF
      END DO
    ENDIF

    IF (BUCONF%LBUDGET_TH) ZWKBUD(:,:) = UNPACK(PZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:)
    IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'GMLT', ZWKBUD)
    IF (BUCONF%LBUDGET_RR) ZWKBUD(:,:) = UNPACK(PRRS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'GMLT', ZWKBUD)
    IF (BUCONF%LBUDGET_RG) ZWKBUD(:,:) = UNPACK(PRGS(:)*PZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0)
    IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'GMLT', ZWKBUD)

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RG',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_FAST_RG

END MODULE MODE_RAIN_ICE_OLD_FAST_RG
