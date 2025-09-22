!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_SEDIMENTATION_STAT

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_SEDIMENTATION_STAT(D, CST, ICEP, ICED,                 &
                                             KRR, OSEDIC, PTSTEP, KKL, IKB, IKE, &
                                             PDZZ, PRHODJ, PRHODREF, PPABST,     &
                                             PTHT, PRCT, PRRT, PRST, PRGT,       &
                                             PRCS, PRRS, PRIS, PRSS, PRGS,       &
                                             PINPRC, PINPRR, PINPRS, PINPRG,     &
                                             ZRAY, ZLBC, ZFSEDC, ZCONC3D,        &
                                             PRHT, PRHS, PINPRH, PFPR)

    USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
    USE MODD_CST,             ONLY: CST_T
    USE MODD_RAIN_ICE_PARAM_n,  ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR_n,  ONLY: RAIN_ICE_DESCR_T
    USE YOMHOOK,              ONLY: LHOOK, DR_HOOK, JPHOOK

!*      0. DECLARATIONS
!          ------------
!
    IMPLICIT NONE

    TYPE(DIMPHYEX_T),       INTENT(IN) :: D
    TYPE(CST_T),            INTENT(IN) :: CST
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED

    INTEGER, INTENT(IN) :: KRR
    LOGICAL, INTENT(IN) :: OSEDIC ! Switch for droplet sedim.

    REAL,    INTENT(IN) :: PTSTEP  ! Double Time step
    INTEGER, INTENT(IN) :: KKL !vert. levels type 1=MNH -1=ARO
    INTEGER, INTENT(IN) :: IKB
    INTEGER, INTENT(IN) :: IKE

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PDZZ     ! Layer thickness (m)
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODJ   ! Dry density * Jacobian
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODREF ! Reference density
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PPABST   ! absolute pressure at t

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PTHT ! Theta at time t
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRCT ! Cloud water m.r. at t
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRRT ! Rain water m.r. at t
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRST ! Snow/aggregate m.r. at t
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRGT ! Graupel/hail m.r. at t

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRCS ! Cloud water m.r. source
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRRS ! Rain water m.r. source
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRIS ! Pristine ice m.r. source
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRSS ! Snow/aggregate m.r. source
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRGS ! Graupel m.r. source

    REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRC ! Cloud instant precip
    REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRR ! Rain instant precip
    REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRS ! Snow instant precip
    REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRG ! Graupel instant precip

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: ZRAY    ! Cloud Mean radius
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: ZLBC    ! XLBC weighted by sea fraction
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: ZFSEDC
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: ZCONC3D !  droplet concentration m-3

    REAL, DIMENSION(D%NIT,D%NKT),     OPTIONAL, INTENT(IN)    :: PRHT   ! Hail m.r. at t
    REAL, DIMENSION(D%NIT,D%NKT),     OPTIONAL, INTENT(INOUT) :: PRHS   ! Hail m.r. source
    REAL, DIMENSION(D%NIT),           OPTIONAL, INTENT(OUT)   :: PINPRH ! Hail instant precip
    REAL, DIMENSION(D%NIT,D%NKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR   ! upper-air precipitation fluxes

    REAL, DIMENSION(D%NIT, 0:D%NKT+1) :: ZWSED   ! sedimentation fluxes
    REAL, DIMENSION(D%NIT, 0:D%NKT+1) :: ZWSEDW1 ! sedimentation speed
    REAL, DIMENSION(D%NIT, 0:D%NKT+1) :: ZWSEDW2 ! sedimentation speed

    REAL, DIMENSION(D%NIT, D%NKT)     :: ZW ! work array

    REAL :: ZP1,ZP2,ZH,ZZWLBDA,ZZWLBDC,ZZCC
    REAL, DIMENSION(D%NIT) :: ZQP
    INTEGER :: JI,JK
    INTEGER :: JCOUNT, JL
    INTEGER, DIMENSION(D%NIT) :: I1
    LOGICAL, DIMENSION(D%NIT) :: GMASK

    REAL, DIMENSION(D%NIT,D%NKT) :: ZPRCS, ZPRRS, ZPRSS, ZPRGS, ZPRHS ! Mixing ratios created during the time step

    REAL, DIMENSION(SIZE(ICED%XRTMIN)) :: ZRTMIN

    REAL :: ZINVTSTEP

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SEDIMENTATION_STAT',0,ZHOOK_HANDLE)
!
!*       2.    compute the fluxes
!
    ZINVTSTEP = 1./PTSTEP
    ZRTMIN(:) = ICED%XRTMIN(:) * ZINVTSTEP
!
    IF (OSEDIC) THEN
      ZPRCS(D%NIB:D%NIE,:) = 0.0
      ZPRCS(D%NIB:D%NIE,:) = PRCS(D%NIB:D%NIE,:) - PRCT(D%NIB:D%NIE,:) * ZINVTSTEP
      PRCS(D%NIB:D%NIE,:)  = PRCT(D%NIB:D%NIE,:) * ZINVTSTEP
    END IF
    ZPRRS(D%NIB:D%NIE,:) = 0.0
    ZPRSS(D%NIB:D%NIE,:) = 0.0
    ZPRGS(D%NIB:D%NIE,:) = 0.0
    IF (KRR == 7) ZPRHS(D%NIB:D%NIE,:) = 0.0
!
    ZPRRS(D%NIB:D%NIE,:) = PRRS(D%NIB:D%NIE,:)-PRRT(D%NIB:D%NIE,:)* ZINVTSTEP
    ZPRSS(D%NIB:D%NIE,:) = PRSS(D%NIB:D%NIE,:)-PRST(D%NIB:D%NIE,:)* ZINVTSTEP
    ZPRGS(D%NIB:D%NIE,:) = PRGS(D%NIB:D%NIE,:)-PRGT(D%NIB:D%NIE,:)* ZINVTSTEP
    IF (KRR == 7) ZPRHS(D%NIB:D%NIE,:) = PRHS(D%NIB:D%NIE,:)-PRHT(D%NIB:D%NIE,:)* ZINVTSTEP
    PRRS(D%NIB:D%NIE,:)  = PRRT(D%NIB:D%NIE,:)* ZINVTSTEP
    PRSS(D%NIB:D%NIE,:)  = PRST(D%NIB:D%NIE,:)* ZINVTSTEP
    PRGS(D%NIB:D%NIE,:)  = PRGT(D%NIB:D%NIE,:)* ZINVTSTEP
    IF (KRR == 7) PRHS(D%NIB:D%NIE,:)  = PRHT(D%NIB:D%NIE,:)* ZINVTSTEP
!
    IF (OSEDIC) PRCS(D%NIB:D%NIE,:) = PRCS(D%NIB:D%NIE,:) + ZPRCS(D%NIB:D%NIE,:)
    PRRS(D%NIB:D%NIE,:) = PRRS(D%NIB:D%NIE,:) + ZPRRS(D%NIB:D%NIE,:)
    PRSS(D%NIB:D%NIE,:) = PRSS(D%NIB:D%NIE,:) + ZPRSS(D%NIB:D%NIE,:)
    PRGS(D%NIB:D%NIE,:) = PRGS(D%NIB:D%NIE,:) + ZPRGS(D%NIB:D%NIE,:)
    IF ( KRR == 7 ) PRHS(D%NIB:D%NIE,:) = PRHS(D%NIB:D%NIE,:) + ZPRHS(D%NIB:D%NIE,:)
    DO JK = D%NKTB , D%NKTE
      DO JI = D%NIB , D%NIE
        ZW(JI,JK) =PTSTEP/(PRHODREF(JI,JK)* PDZZ(JI,JK) )
      END DO
    END DO
!
!*       2.1   for cloud
!
    IF (OSEDIC) THEN
      PRCS(D%NIB:D%NIE,:) = PRCS(D%NIB:D%NIE,:) * PTSTEP
      ZWSED(D%NIB:D%NIE,:) = 0.
      ZWSEDW1(D%NIB:D%NIE,:) = 0.
      ZWSEDW2(D%NIB:D%NIE,:) = 0.

      ! calculation of P1, P2 and sedimentation flux
      DO JK = IKE , IKB, -1*KKL
        !estimation of q' taking into account incomming ZWSED
        DO JI = D%NIB , D%NIE
          ZQP(JI)=ZWSED(JI,JK+KKL)*ZW(JI,JK)
        END DO

        GMASK(:)=.FALSE.
        GMASK(D%NIB:D%NIE)=(PRCS(D%NIB:D%NIE,JK) > ZRTMIN(2) .AND. PRCT(D%NIB:D%NIE,JK) > ZRTMIN(2)) .OR. &
                          &ZQP(D%NIB:D%NIE) > ZRTMIN(2)
        CALL COUNTJV2(D, JCOUNT, GMASK(:), I1(:))

        DO JL=1, JCOUNT
          JI=I1(JL)
          !calculation of w
          ! mars 2009 : ajout d'un test
          IF(PRCS(JI,JK) > ZRTMIN(2) .AND. PRCT(JI,JK) > ZRTMIN(2)) THEN
            ZZWLBDA=6.6E-8*(101325./PPABST(JI,JK))*(PTHT(JI,JK)/293.15)
            ZZWLBDC=(ZLBC(JI,JK)*ZCONC3D(JI,JK)  &
                  &/(PRHODREF(JI,JK)*PRCT(JI,JK)))**ICED%XLBEXC
            ZZCC=ICED%XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY(JI,JK)) !! ZCC  : Fall speed
            ZWSEDW1(JI,JK)=PRHODREF(JI,JK)**(-ICED%XCEXVT ) *   &
               &  ZZWLBDC**(-ICED%XDC)*ZZCC*ZFSEDC(JI,JK)
          ENDIF
          IF ( ZQP(JI) > ZRTMIN(2) ) THEN
            ZZWLBDA=6.6E-8*(101325./PPABST(JI,JK))*(PTHT(JI,JK)/293.15)
            ZZWLBDC=(ZLBC(JI,JK)*ZCONC3D(JI,JK)  &
                  &/(PRHODREF(JI,JK)*ZQP(JI)))**ICED%XLBEXC
            ZZCC=ICED%XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY(JI,JK)) !! ZCC  : Fall speed
            ZWSEDW2(JI,JK)=PRHODREF(JI,JK)**(-ICED%XCEXVT ) *   &
               &  ZZWLBDC**(-ICED%XDC)*ZZCC*ZFSEDC(JI,JK)
          ENDIF
        ENDDO

        DO JI = D%NIB, D%NIE
          ZH=PDZZ(JI,JK)
          ZP1 = MIN(1., ZWSEDW1(JI,JK) * PTSTEP / ZH)
          ! mars 2009 : correction : ZWSEDW1 =>  ZWSEDW2
          IF (ZWSEDW2(JI,JK) /= 0.) THEN
            ZP2 = MAX(0.,1 -  ZH &
          &  / (PTSTEP*ZWSEDW2(JI,JK)) )
          ELSE
            ZP2 = 0.
          ENDIF
          ZWSED(JI,JK)=ZP1*PRHODREF(JI,JK)*&
          &ZH*PRCS(JI,JK)&
          &* ZINVTSTEP+ ZP2 * ZWSED(JI,JK+KKL)
        ENDDO
      ENDDO

      DO JK = D%NKTB , D%NKTE
        PRCS(D%NIB:D%NIE,JK) = PRCS(D%NIB:D%NIE,JK) + ZW(D%NIB:D%NIE,JK)*(ZWSED(D%NIB:D%NIE,JK+KKL)-ZWSED(D%NIB:D%NIE,JK))
      END DO

      IF (PRESENT(PFPR)) THEN
        DO JK = D%NKTB , D%NKTE
          PFPR(D%NIB:D%NIE,JK,2)=ZWSED(D%NIB:D%NIE,JK)
        ENDDO
      ENDIF

      PINPRC(D%NIB:D%NIE) = ZWSED(D%NIB:D%NIE,IKB)/CST%XRHOLW ! in m/s
      PRCS(D%NIB:D%NIE,:) = PRCS(D%NIB:D%NIE,:) * ZINVTSTEP
    ENDIF
!
!*       2.2   for rain
!
    PRRS(D%NIB:D%NIE,:) = PRRS(D%NIB:D%NIE,:) * PTSTEP
    ZWSED(D%NIB:D%NIE,:) = 0.
    ZWSEDW1(D%NIB:D%NIE,:) = 0.
    ZWSEDW2(D%NIB:D%NIE,:) = 0.

    ! calculation of ZP1, ZP2 and sedimentation flux
    DO JK = IKE , IKB, -1*KKL
      !estimation of q' taking into account incomming ZWSED
      DO JI = D%NIB, D%NIE
        ZQP(JI)=ZWSED(JI,JK+KKL)*ZW(JI,JK)
      END DO

      GMASK(:)=.FALSE.
      GMASK(D%NIB:D%NIE)=PRRS(D%NIB:D%NIE,JK) > ZRTMIN(3) .OR. ZQP(D%NIB:D%NIE) > ZRTMIN(3)
      CALL COUNTJV2(D, JCOUNT, GMASK(:), I1(:))
      DO JL=1, JCOUNT
        JI=I1(JL)

        !calculation of w
        IF ( PRRS(JI,JK) > ZRTMIN(3) ) THEN
          ZWSEDW1(JI,JK)= ICEP%XFSEDR *PRRS(JI,JK)**(ICEP%XEXSEDR-1)* &
          PRHODREF(JI,JK)**(ICEP%XEXSEDR-ICED%XCEXVT-1)
        ENDIF

        IF ( ZQP(JI) > ZRTMIN(3) ) THEN
          ZWSEDW2(JI,JK)= ICEP%XFSEDR *(ZQP(JI))**(ICEP%XEXSEDR-1)* &
          PRHODREF(JI,JK)**(ICEP%XEXSEDR-ICED%XCEXVT-1)
        ENDIF
      ENDDO

      DO JI = D%NIB, D%NIE
        ZH=PDZZ(JI,JK)
        ZP1 = MIN(1., ZWSEDW1(JI,JK) * PTSTEP / ZH )
        IF (ZWSEDW2(JI,JK) /= 0.) THEN
          ZP2 = MAX(0.,1 -  ZH/(PTSTEP*ZWSEDW2(JI,JK)) )
        ELSE
          ZP2 = 0.
        ENDIF
        ZWSED(JI,JK)=ZP1*PRHODREF(JI,JK)*ZH*PRRS(JI,JK)*ZINVTSTEP+ZP2*ZWSED(JI,JK+KKL)
      ENDDO
    ENDDO

    DO JK = D%NKTB , D%NKTE
      PRRS(D%NIB:D%NIE,JK) = PRRS(D%NIB:D%NIE,JK) + ZW(D%NIB:D%NIE,JK)*(ZWSED(D%NIB:D%NIE,JK+KKL)-ZWSED(D%NIB:D%NIE,JK))
    ENDDO
    IF (PRESENT(PFPR)) THEN
      DO JK = D%NKTB , D%NKTE
        PFPR(D%NIB:D%NIE,JK,3)=ZWSED(D%NIB:D%NIE,JK)
      ENDDO
    ENDIF
    PINPRR(D%NIB:D%NIE) = ZWSED(D%NIB:D%NIE,IKB)/CST%XRHOLW ! in m/s
    PRRS(D%NIB:D%NIE,:) = PRRS(D%NIB:D%NIE,:) * ZINVTSTEP
!
!*       2.3   for pristine ice
!
    PRIS(D%NIB:D%NIE,:) = PRIS(D%NIB:D%NIE,:) * PTSTEP
    ZWSED(D%NIB:D%NIE,:) = 0.
    ZWSEDW1(D%NIB:D%NIE,:) = 0.
    ZWSEDW2(D%NIB:D%NIE,:) = 0.

    ! calculation of ZP1, ZP2 and sedimentation flux
    DO JK = IKE , IKB, -1*KKL
      !estimation of q' taking into account incomming ZWSED
      DO JI = D%NIB, D%NIE
        ZQP(JI)=ZWSED(JI,JK+KKL)*ZW(JI,JK)
      ENDDO

      GMASK(:)=.FALSE.
      GMASK(D%NIB:D%NIE)=PRIS(D%NIB:D%NIE,JK) > MAX(ZRTMIN(4), 1.0E-7) .OR. ZQP(D%NIB:D%NIE) > MAX(ZRTMIN(4), 1.0E-7)
      CALL COUNTJV2(D, JCOUNT, GMASK(:), I1(:))

      DO JL=1, JCOUNT
        JI=I1(JL)

        !calculation of w
        IF ( PRIS(JI,JK) > MAX(ZRTMIN(4),1.0E-7 ) ) THEN
          ZWSEDW1(JI,JK)= ICEP%XFSEDI *  &
          &  PRHODREF(JI,JK)**(ICED%XCEXVT) * & !    McF&H
          &  MAX( 0.05E6,-0.15319E6-0.021454E6* &
          &  ALOG(PRHODREF(JI,JK)*PRIS(JI,JK)) )**ICEP%XEXCSEDI
        ENDIF

        IF ( ZQP(JI) > MAX(ZRTMIN(4),1.0E-7 ) ) THEN
          ZWSEDW2(JI,JK)= ICEP%XFSEDI *  &
          &  PRHODREF(JI,JK)**(ICED%XCEXVT) * & !    McF&H
          &  MAX( 0.05E6,-0.15319E6-0.021454E6* &
          &  ALOG(PRHODREF(JI,JK)*ZQP(JI)) )**ICEP%XEXCSEDI
        ENDIF
      ENDDO

      DO JI = D%NIB, D%NIE
        ZH=PDZZ(JI,JK)
        ZP1 = MIN(1., ZWSEDW1(JI,JK) * PTSTEP / ZH )

        IF (ZWSEDW2(JI,JK) /= 0.) THEN
          ZP2 = MAX(0.,1 - ZH/(PTSTEP*ZWSEDW2(JI,JK)))
        ELSE
          ZP2 = 0.
        ENDIF

        ZWSED(JI,JK)=ZP1*PRHODREF(JI,JK)*ZH*PRIS(JI,JK)*ZINVTSTEP+ZP2*ZWSED(JI,JK+KKL)
      ENDDO
    ENDDO

    DO JK = D%NKTB , D%NKTE
      PRIS(D%NIB:D%NIE,JK) = PRIS(D%NIB:D%NIE,JK) + ZW(D%NIB:D%NIE,JK)*(ZWSED(D%NIB:D%NIE,JK+KKL)-ZWSED(D%NIB:D%NIE,JK))
    ENDDO

    IF (PRESENT(PFPR)) THEN
      DO JK = D%NKTB , D%NKTE
        PFPR(D%NIB:D%NIE,JK,4)=ZWSED(D%NIB:D%NIE,JK)
      ENDDO
    ENDIF

    PRIS(D%NIB:D%NIE,:) = PRIS(D%NIB:D%NIE,:) * ZINVTSTEP

    PINPRS(D%NIB:D%NIE) = ZWSED(D%NIB:D%NIE,IKB)/CST%XRHOLW
!
!*       2.4   for aggregates/snow
!
    PRSS(D%NIB:D%NIE,:) = PRSS(D%NIB:D%NIE,:) * PTSTEP
    ZWSED(D%NIB:D%NIE,:) = 0.
    ZWSEDW1(D%NIB:D%NIE,:) = 0.
    ZWSEDW2(D%NIB:D%NIE,:) = 0.

    ! calculation of ZP1, ZP2 and sedimentation flux
    DO JK = IKE , IKB, -1*KKL
      !estimation of q' taking into account incomming ZWSED
      ZQP(D%NIB:D%NIE)=ZWSED(D%NIB:D%NIE,JK+KKL)*ZW(D%NIB:D%NIE,JK)

      GMASK(:)=.FALSE.
      GMASK(D%NIB:D%NIE)=PRSS(D%NIB:D%NIE,JK) > ZRTMIN(5) .OR. ZQP(D%NIB:D%NIE) > ZRTMIN(5)
      CALL COUNTJV2(D, JCOUNT, GMASK(:), I1(:))
      DO JL=1, JCOUNT
        JI=I1(JL)

        !calculation of w
        IF (PRSS(JI,JK) > ZRTMIN(5) ) THEN
          ZWSEDW1(JI,JK)=ICEP%XFSEDS*(PRSS(JI,JK))**(ICEP%XEXSEDS-1)*&
          PRHODREF(JI,JK)**(ICEP%XEXSEDS-ICED%XCEXVT-1)
        ENDIF

        IF ( ZQP(JI) > ZRTMIN(5) ) THEN
          ZWSEDW2(JI,JK)=ICEP%XFSEDS*(ZQP(JI))**(ICEP%XEXSEDS-1)*&
          PRHODREF(JI,JK)**(ICEP%XEXSEDS-ICED%XCEXVT-1)
        ENDIF
      ENDDO

      DO JI = D%NIB, D%NIE
        ZH=PDZZ(JI,JK)
        ZP1 = MIN(1., ZWSEDW1(JI,JK) * PTSTEP / ZH )
        IF (ZWSEDW2(JI,JK) /= 0.) THEN
          ZP2 = MAX(0.,1 - ZH/(PTSTEP*ZWSEDW2(JI,JK)))
        ELSE
          ZP2 = 0.
        ENDIF

        ZWSED(JI,JK)=ZP1*PRHODREF(JI,JK)*ZH*PRSS(JI,JK)* ZINVTSTEP+ZP2*ZWSED(JI,JK+KKL)
      ENDDO
    ENDDO

    DO JK = D%NKTB , D%NKTE
      PRSS(D%NIB:D%NIE,JK) = PRSS(D%NIB:D%NIE,JK) + ZW(D%NIB:D%NIE,JK)*(ZWSED(D%NIB:D%NIE,JK+KKL)-ZWSED(D%NIB:D%NIE,JK))
    ENDDO

    IF (PRESENT(PFPR)) THEN
      DO JK = D%NKTB , D%NKTE
        PFPR(D%NIB:D%NIE,JK,5)=ZWSED(D%NIB:D%NIE,JK)
      ENDDO
    ENDIF

    PINPRS(D%NIB:D%NIE) = ZWSED(D%NIB:D%NIE,IKB)/CST%XRHOLW + PINPRS(D%NIB:D%NIE)    ! in m/s (add ice fall)

    PRSS(D%NIB:D%NIE,:) = PRSS(D%NIB:D%NIE,:) * ZINVTSTEP
!
!*       2.5   for graupeln
!
    PRGS(D%NIB:D%NIE,:) = PRGS(D%NIB:D%NIE,:) * PTSTEP
    ZWSED(D%NIB:D%NIE,:) = 0.
    ZWSEDW1(D%NIB:D%NIE,:) = 0.
    ZWSEDW2(D%NIB:D%NIE,:) = 0.

    ! calculation of ZP1, ZP2 and sedimentation flux
    DO JK = IKE,  IKB, -1*KKL
      !estimation of q' taking into account incomming ZWSED
      ZQP(D%NIB:D%NIE)=ZWSED(D%NIB:D%NIE,JK+KKL)*ZW(D%NIB:D%NIE,JK)

      GMASK(:)=.FALSE.
      GMASK(D%NIB:D%NIE)=PRGS(D%NIB:D%NIE,JK) > ZRTMIN(6) .OR. ZQP(D%NIB:D%NIE) > ZRTMIN(6)
      CALL COUNTJV2(D, JCOUNT, GMASK(:), I1(:))

      DO JL = 1,JCOUNT
        JI=I1(JL)

        !calculation of w
        IF ( PRGS(JI,JK) > ZRTMIN(6) ) THEN
          ZWSEDW1(JI,JK)= ICEP%XFSEDG*(PRGS(JI,JK))**(ICEP%XEXSEDG-1) * &
                                  PRHODREF(JI,JK)**(ICEP%XEXSEDG-ICED%XCEXVT-1)
        ENDIF

        IF ( ZQP(JI) > ZRTMIN(6) ) THEN
          ZWSEDW2(JI,JK)= ICEP%XFSEDG*(ZQP(JI))**(ICEP%XEXSEDG-1) * &
                                  PRHODREF(JI,JK)**(ICEP%XEXSEDG-ICED%XCEXVT-1)
        ENDIF
      ENDDO

      DO JI = D%NIB, D%NIE
        ZH=PDZZ(JI,JK)
        ZP1 = MIN(1., ZWSEDW1(JI,JK) * PTSTEP / ZH )
        IF (ZWSEDW2(JI,JK) /= 0.) THEN
          ZP2 = MAX(0.,1 - ZH/(PTSTEP*ZWSEDW2(JI,JK)))
        ELSE
          ZP2 = 0.
        ENDIF
        ZWSED(JI,JK) = ZP1*PRHODREF(JI,JK)*ZH*PRGS(JI,JK)*ZINVTSTEP+ZP2*ZWSED(JI,JK+KKL)
      ENDDO
    ENDDO

    DO JK = D%NKTB , D%NKTE
      PRGS(D%NIB:D%NIE,JK) = PRGS(D%NIB:D%NIE,JK) + ZW(D%NIB:D%NIE,JK)*(ZWSED(D%NIB:D%NIE,JK+KKL)-ZWSED(D%NIB:D%NIE,JK))
    ENDDO

    IF (PRESENT(PFPR)) THEN
      DO JK = D%NKTB , D%NKTE
        PFPR(D%NIB:D%NIE,JK,6)=ZWSED(D%NIB:D%NIE,JK)
      ENDDO
    ENDIF

    PINPRG(D%NIB:D%NIE) = ZWSED(D%NIB:D%NIE,IKB)/CST%XRHOLW ! in m/s

    PRGS(D%NIB:D%NIE,:) = PRGS(D%NIB:D%NIE,:) * ZINVTSTEP
!
!*       2.6   for hail
!
    IF ( KRR == 7 ) THEN
      PRHS(D%NIB:D%NIE,:) = PRHS(D%NIB:D%NIE,:) * PTSTEP
      ZWSED(D%NIB:D%NIE,:) = 0.
      ZWSEDW1(D%NIB:D%NIE,:) = 0.
      ZWSEDW2(D%NIB:D%NIE,:) = 0.

      ! calculation of ZP1, ZP2 and sedimentation flux
      DO JK = IKE, IKB, -1*KKL
        !estimation of q' taking into account incomming ZWSED
        ZQP(D%NIB:D%NIE)=ZWSED(D%NIB:D%NIE,JK+KKL)*ZW(D%NIB:D%NIE,JK)

        GMASK(:)=.FALSE.
        GMASK(D%NIB:D%NIE)=PRHS(D%NIB:D%NIE,JK)+ZQP(:) > ZRTMIN(7) .OR. ZQP(D%NIB:D%NIE) > ZRTMIN(7)
        CALL COUNTJV2(D, JCOUNT, GMASK(:), I1(:))

        DO JL=1, JCOUNT
          JI=I1(JL)

          !calculation of w
          IF ((PRHS(JI,JK)+ZQP(JI)) > ZRTMIN(7)) THEN
            ZWSEDW1(JI,JK)= ICEP%XFSEDH  * (PRHS(JI,JK))**(ICEP%XEXSEDH-1) *   &
                                      PRHODREF(JI,JK)**(ICEP%XEXSEDH-ICED%XCEXVT-1)
          ENDIF

          IF ( ZQP(JI) > ZRTMIN(7) ) THEN
            ZWSEDW2(JI,JK) = ICEP%XFSEDH * ZQP(JI)**(ICEP%XEXSEDH-1) *   &
                                    PRHODREF(JI,JK)**(ICEP%XEXSEDH-ICED%XCEXVT-1)
          ENDIF
        ENDDO

        DO JI = D%NIB, D%NIE
          ZH=PDZZ(JI,JK)
          ZP1 = MIN(1., ZWSEDW1(JI,JK) * PTSTEP/ZH)
          IF (ZWSEDW2(JI,JK) /= 0.) THEN
            ZP2 = MAX(0.,1 - ZH/(PTSTEP*ZWSEDW2(JI,JK)))
          ELSE
            ZP2 = 0.
          ENDIF

          ZWSED(JI,JK) = ZP1*PRHODREF(JI,JK)*ZH*PRHS(JI,JK)*ZINVTSTEP + ZP2*ZWSED(JI,JK+KKL)
        ENDDO
      ENDDO

      DO JK = D%NKTB , D%NKTE
        PRHS(D%NIB:D%NIE,JK) = PRHS(D%NIB:D%NIE,JK) + ZW(D%NIB:D%NIE,JK)*(ZWSED(D%NIB:D%NIE,JK+KKL)-ZWSED(D%NIB:D%NIE,JK))
      ENDDO

      IF (PRESENT(PFPR)) THEN
        DO JK = D%NKTB , D%NKTE
          PFPR(D%NIB:D%NIE,JK,7)=ZWSED(D%NIB:D%NIE,JK)
        ENDDO
      ENDIF

      PINPRH(D%NIB:D%NIE) = ZWSED(D%NIB:D%NIE,IKB)/CST%XRHOLW ! in m/s

      PRHS(D%NIB:D%NIE,:) = PRHS(D%NIB:D%NIE,:) * ZINVTSTEP

    ENDIF

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SEDIMENTATION_STAT',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_SEDIMENTATION_STAT

  SUBROUTINE COUNTJV2(D, IC, LTAB, I1)
    USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
    IMPLICIT NONE

    TYPE(DIMPHYEX_T),       INTENT(IN) :: D
    INTEGER, INTENT(OUT) :: IC
    LOGICAL, DIMENSION(D%NIT), INTENT(IN)  :: LTAB ! Mask
    INTEGER, DIMENSION(D%NIT), INTENT(OUT) :: I1   ! Used to replace the COUNT and PACK
    INTEGER :: JI

    IC = 0
    DO JI = D%NIB, D%NIE
      IF(LTAB(JI)) THEN
        IC = IC +1
        I1(IC) = JI
      END IF
    END DO

  END SUBROUTINE COUNTJV2

END MODULE MODE_RAIN_ICE_OLD_SEDIMENTATION_STAT
