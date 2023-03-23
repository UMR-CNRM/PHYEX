!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_NUCLEATION

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_NUCLEATION(D, CST, ICEP, KSIZE, OCND2, LMODICEDEP, KRR, PTSTEP, &
                                     PTHT, PPABST, PEXNREF, PICLDFR, PRHODJ, PRHODREF, &
                                     PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                     PTHS, PRVS, PRIS, PCIT, &
                                     PICENU, PT, PZZZ, &
                                     PRHT)
!
    USE PARKIND1,             ONLY: JPRB
    USE YOMHOOK ,             ONLY: LHOOK, DR_HOOK
    USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
    USE MODD_CST,             ONLY: CST_T
    USE MODE_TIWMX,           ONLY: ESATI, ESATW, AM3, REDIN
    USE MODD_RAIN_ICE_PARAM_N,ONLY: RAIN_ICE_PARAM_T
!
!*      0. DECLARATIONS
!          ------------
!
    IMPLICIT NONE

    TYPE(DIMPHYEX_T), INTENT(IN)       :: D
    TYPE(CST_T), INTENT(IN)            :: CST 
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP

    INTEGER, INTENT(IN) :: KSIZE

    LOGICAL, INTENT(IN) :: OCND2      ! Logical switch to separate liquid and ice
    LOGICAL, INTENT(IN) :: LMODICEDEP ! Logical switch for alternative dep/evap of ice

    INTEGER, INTENT(IN) :: KRR        ! Number of moist variable

    REAL, INTENT(IN) :: PTSTEP  ! Double Time step

    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PTHT     ! Theta at time t
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PPABST   ! absolute pressure at t
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PEXNREF  ! Reference Exner function
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PICLDFR  ! ice cloud fraction
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PRHODJ   ! Dry density * Jacobian
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PRHODREF ! Reference density
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PRVT     ! Water vapor m.r. at t
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PRCT     ! Cloud water m.r. at t
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PRRT     ! Rain water m.r. at t
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PRIT     ! Pristine ice m.r. at t
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PRST     ! Snow/aggregate m.r. at t
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PRGT     ! Graupel/hail m.r. at t
!
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(INOUT) :: PTHS ! Theta source
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(INOUT) :: PRVS ! Water vapor m.r. source
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(INOUT) :: PRIS ! Pristine ice m.r. source
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(INOUT) :: PCIT ! Pristine ice n.c. at t
!
    REAL, DIMENSION(D%NIT), INTENT(IN) :: PICENU
!
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PT   ! Temperature
    REAL, DIMENSION(D%NIT, D%NKT), INTENT(IN) :: PZZZ ! Temperature

    REAL, DIMENSION(D%NIT, D%NKT), OPTIONAL, INTENT(IN) :: PRHT ! Hail m.r. at t
!
!*       0.2  declaration of local variables
!
    REAL, DIMENSION(KSIZE) :: ZRVT    ! Water vapor m.r. at t

    INTEGER , DIMENSION(KSIZE) :: I1,I2 ! Used to replace the COUNT
    INTEGER                    :: JI, JK, JL ! and PACK intrinsics
!
!-------------------------------------------------------------------------------
!
    REAL, DIMENSION(KSIZE) :: ZZZ     ! height of model layer (m)
    REAL, DIMENSION(KSIZE) :: ZREDIN  ! Reduction of IN concentration between 0 and -25 C
    REAL, DIMENSION(KSIZE) :: ZCIT    ! Pristine ice conc. at t
    REAL, DIMENSION(KSIZE) :: ZZT     ! Temperatur
    REAL, DIMENSION(KSIZE) :: ZPRES   ! Pressure
    REAL, DIMENSION(KSIZE) :: ZZICENU ! Pressure
    REAL, DIMENSION(KSIZE) :: ZAM3    ! Meyers IN concentration function
    REAL, DIMENSION(KSIZE) :: ZESI    ! saturation pressure over ice
    REAL, DIMENSION(KSIZE) :: ZESW    ! saturation pressure over water
    REAL, DIMENSION(KSIZE) :: ZUSW    ! Undersaturation over water
    REAL, DIMENSION(KSIZE) :: ZSSI    ! Supersaturation over ice
    REAL, DIMENSION(KSIZE) :: ZSIFRC  ! subgridscale fraction with supersaturation with
                                               ! respect to ice.
    REAL, DIMENSION(KSIZE) :: ZZW     ! Work array
    REAL, DIMENSION(D%NIT,D%NKT) :: ZW      ! work array
!
!   compute the temperature and the pressure
!
    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_OLD_NUCLEATION',0,ZHOOK_HANDLE)

    IF( KSIZE >= 1 ) THEN
!
!     optimization by looking for locations where
!     the temperature is negative only !!!
!
      JL = 0
      DO JK = 1, D%NKT
        DO JI = 1, D%NIT
          IF (PT(JI, JK) < CST%XTT) THEN
            JL = JL + 1
            I1(JL) = JI
            I2(JL) = JK
          ENDIF
        ENDDO
      ENDDO


      DO JL=1, KSIZE
        ZRVT(JL) = PRVT(I1(JL),I2(JL))
        ZCIT(JL) = PCIT(I1(JL),I2(JL))
        ZZT(JL) = PT(I1(JL),I2(JL))
        ZPRES(JL) = PPABST(I1(JL),I2(JL))
        ZZICENU(JL) = PICENU(I1(JL))

        IF (OCND2) THEN
          ZZZ(JL) = PZZZ(I1(JL),I2(JL))
          ZESI(JL) = ESATI(ZZT(JL))
          ZESW(JL) = ESATW(ZZT(JL))
          ZAM3(JL) = AM3(MAX(ICEP%XFRMIN(27),ZZT(JL))) ! Avoid too high IN for very low temp.
          ZREDIN(JL) = REDIN(ZZT(JL))
          ZSIFRC(JL) = PICLDFR(I1(JL),I2(JL))
        ENDIF

      ENDDO

      IF(OCND2)THEN ! try to do some optimazation :

        DO JL = 1, KSIZE
          ZZW(JL) = MIN(ZPRES(JL)/2., ZESI(JL))             ! safety limitation  es_i
          ZZW(JL) = MIN(ZPRES(JL)/2., ZZW(JL))             ! safety limitation
          ZSSI(JL) = ZRVT(JL)*(ZPRES(JL)-ZZW(JL)) / (CST%XEPSILO * ZZW(JL)) - 1.0
        ENDDO
                                                      ! Supersaturation over ice
        DO JL = 1, KSIZE
          ZUSW(JL) = MIN(ZPRES(JL)/2.,ZESW(JL))            ! safety limitation   es_w
          ZUSW(JL) = (ZUSW(JL)/ZZW(JL))*((ZPRES(JL)-ZZW(JL))/(ZPRES(JL)-ZUSW(JL))) - 1.0
        ENDDO
                             ! Supersaturation of saturated water vapor over ice
      ELSE

        DO JL = 1, KSIZE
          ZZW(JL)  = EXP( CST%XALPI - CST%XBETAI/ZZT(JL) - CST%XGAMI*ALOG(ZZT(JL))) ! es_i
          ZZW(JL)  = MIN(ZPRES(JL)/2., ZZW(JL))             ! safety limitation
          ZSSI(JL) = ZRVT(JL)*( ZPRES(JL)-ZZW(JL) ) / ( CST%XEPSILO * ZZW(JL) ) - 1.0
        ENDDO
                                                      ! Supersaturation over ice
        DO JL = 1, KSIZE
          ZUSW(JL) = EXP( CST%XALPW - CST%XBETAW/ZZT(JL) - CST%XGAMW*ALOG(ZZT(JL)))  ! es_w
          ZUSW(JL) = MIN(ZPRES(JL)/2.,ZUSW(JL))            ! safety limitation
                             ! Supersaturation of saturated water vapor over ice
          ZUSW(JL) = (ZUSW(JL)/ZZW(JL))*((ZPRES(JL)-ZZW(JL))/(ZPRES(JL)-ZUSW(JL))) - 1.0
        ENDDO

      ENDIF
!
!*       3.1     compute the heterogeneous nucleation source: RVHENI
!
!*       3.1.1   compute the cloud ice concentration
!
      ZZW(:) = 0.0

      ZSSI(:) = MIN( ZSSI(:), ZUSW(:) ) ! limitation of SSi according to SSw=0

      IF(OCND2)THEN

        IF (LMODICEDEP) THEN
          DO JL = 1, KSIZE
            ZZW(JL) = 5.*EXP(0.304*(CST%XTT-ZZT(JL)))
            ZZW(JL) = MIN(1.,MAX(ZSSI(JL)*10.,0.01))*ZZW(JL)
          ENDDO
        ELSE

          DO JL = 1, KSIZE
            ZZW(JL) = ZREDIN(JL)* MAX(0.1,((20000.- MIN(20000.,ZZZ(JL)))/20000.)**4) &
                &   *ZAM3(JL)*(0.0001 + 0.9999*ZSIFRC(JL))
          ENDDO

        ENDIF

      ELSE

        DO JL = 1, KSIZE
          IF ((ZZT(JL)<CST%XTT-5.0) .AND. (ZSSI(JL)>0.0)) THEN
            ZZW(JL) = ICEP%XNU20 * EXP( ICEP%XALPHA2*ZSSI(JL)-ICEP%XBETA2 )
          END IF
        ENDDO

        DO JL = 1, KSIZE
          IF ((ZZT(JL)<=CST%XTT-2.0) .AND. (ZZT(JL)>=CST%XTT-5.0) .AND. (ZSSI(JL)>0.0)) THEN
            ZZW(JL) = MAX(ICEP%XNU20 * EXP(-ICEP%XBETA2),ICEP%XNU10 * EXP(-ICEP%XBETA1*(ZZT(JL)-CST%XTT)) * &
                     (ZSSI(JL)/ZUSW(JL))**ICEP%XALPHA1)
          END IF
        ENDDO
      ENDIF

      ZZW(:) = ZZW(:)*ZZICENU(:) - ZCIT(:)

      IF( MAXVAL(ZZW(:)) > 0.0 ) THEN
!
!*       3.1.2   update the r_i and r_v mixing ratios
!
        ZZW(:) = MIN( ZZW(:),50.E3 ) ! limitation provisoire a 50 l^-1

        IF(.NOT.OCND2)THEN
          ZW(:,:) = UNPACK( ZZW(:),MASK=PT(D%NIB:D%NIE,D%NKTB:D%NKTE) < CST%XTT, FIELD=0.0 )
          ZW(:,:) = MAX( ZW(:,:) ,0.0 ) *ICEP%XMNU0/(PRHODREF(:,:)*PTSTEP)
          PRIS(:,:) = PRIS(:,:) + ZW(:,:)
          PRVS(:,:) = PRVS(:,:) - ZW(:,:)

          IF ( KRR == 7 ) THEN
            PTHS(:,:) = PTHS(:,:) + ZW(:,:)*(CST%XLSTT+(CST%XCPV-CST%XCI)*(PT(:,:)-CST%XTT))   &
                 /( (CST%XCPD + CST%XCPV*PRVT(:,:) + CST%XCL*(PRCT(:,:)+PRRT(:,:))   &
                 + CST%XCI*(PRIT(:,:)+PRST(:,:)+PRGT(:,:)+PRHT(:,:)))*PEXNREF(:,:) )
          ELSE IF( KRR == 6 ) THEN
            PTHS(:,:) = PTHS(:,:) + ZW(:,:)*(CST%XLSTT+(CST%XCPV-CST%XCI)*(PT(:,:)-CST%XTT))   &
                 /( (CST%XCPD + CST%XCPV*PRVT(:,:) + CST%XCL*(PRCT(:,:)+PRRT(:,:))   &
                 + CST%XCI*(PRIT(:,:)+PRST(:,:)+PRGT(:,:)))*PEXNREF(:,:) )
          END IF
        ENDIF
                                 ! f(L_s*(RVHENI))
        ZZW(:) = MAX(ZZW(:)+ZCIT(:),ZCIT(:))

        PCIT(:,:) = MAX(UNPACK(ZZW(:), MASK=PT(D%NIB:D%NIE,D%NKTB:D%NKTE) < CST%XTT, FIELD=0.0), PCIT(:,:))
      END IF

    END IF

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_OLD_NUCLEATION',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_NUCLEATION

END MODULE MODE_RAIN_ICE_OLD_NUCLEATION
