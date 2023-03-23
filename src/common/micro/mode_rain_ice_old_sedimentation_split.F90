!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_OLD_SEDIMENTATION_SPLIT

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RAIN_ICE_OLD_SEDIMENTATION_SPLIT(D, CST, ICEP, ICED, KSIZE, &
                                              KRR, OSEDIC, PTSTEP, KKL, IKB, KSPLITR, &
                                              PDZZ, PRHODJ, PRHODREF, PPABST, &
                                              PTHT, PRCT, PRRT, PRST, PRGT, &
                                              PRCS, PRRS, PRIS, PRSS, PRGS, &
                                              PINPRC, PINPRR, PINPRS, PINPRG, &
                                              ZRAY, ZLBC, ZFSEDC, ZCONC3D,  &
                                              PRHT, PRHS, PINPRH, PFPR)

    USE PARKIND1,             ONLY: JPRB
    USE YOMHOOK ,             ONLY: LHOOK, DR_HOOK
    USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
    USE MODD_CST,             ONLY: CST_T
    USE MODD_RAIN_ICE_PARAM_n,  ONLY: RAIN_ICE_PARAM_T
    USE MODD_RAIN_ICE_DESCR_n,  ONLY: RAIN_ICE_DESCR_T
!
!*      0. DECLARATIONS
!          ------------
    IMPLICIT NONE

    TYPE(DIMPHYEX_T),       INTENT(IN) :: D
    TYPE(CST_T),            INTENT(IN) :: CST
    TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
    TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED

    INTEGER, INTENT(IN) :: KSIZE
    INTEGER, INTENT(IN) :: KRR
    LOGICAL, INTENT(IN) :: OSEDIC ! Switch for droplet sedim.

    REAL,    INTENT(IN) :: PTSTEP  ! Double Time step
    INTEGER, INTENT(IN) :: KKL     !vert. levels type 1=MNH -1=ARO
    INTEGER, INTENT(IN) :: IKB
    INTEGER, INTENT(IN) :: KSPLITR ! Number of small time step

    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PDZZ     ! Layer thickness (m)
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODJ   ! Dry density * Jacobian
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODREF ! Reference density
    REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PPABST  ! absolute pressure at t

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

    REAL, DIMENSION(D%NIT,D%NKT) :: ZPRCS,ZPRRS,ZPRSS,ZPRGS,ZPRHS ! Mixing ratios created during the time step
    INTEGER :: ISEDIMR, ISEDIMC, ISEDIMI, ISEDIMS, ISEDIMG, ISEDIMH

    REAL, DIMENSION(D%NIT, 0:D%NKT+1) :: ZWSED ! sedimentation fluxes

    LOGICAL, DIMENSION(D%NIT,D%NKT) :: GSEDIMR,GSEDIMC, GSEDIMI, GSEDIMS, GSEDIMG, GSEDIMH ! Test where to compute the SED processes
    INTEGER, DIMENSION(D%NIT*D%NKT) :: IC1, IC2
    INTEGER, DIMENSION(D%NIT*D%NKT) :: IR1, IR2
    INTEGER, DIMENSION(D%NIT*D%NKT) :: IS1, IS2
    INTEGER, DIMENSION(D%NIT*D%NKT) :: II1, II2
    INTEGER, DIMENSION(D%NIT*D%NKT) :: IG1, IG2
    INTEGER, DIMENSION(D%NIT*D%NKT) :: IH1, IH2

    INTEGER :: ILISTLENC, ILISTLENR, ILISTLENI, ILISTLENS, ILISTLENG, ILISTLENH

    INTEGER, DIMENSION(KSIZE) :: ILISTR
    INTEGER, DIMENSION(KSIZE) :: ILISTC
    INTEGER, DIMENSION(KSIZE) :: ILISTI
    INTEGER, DIMENSION(KSIZE) :: ILISTS
    INTEGER, DIMENSION(KSIZE) :: ILISTG
    INTEGER, DIMENSION(KSIZE) :: ILISTH

    REAL, DIMENSION(KSIZE) :: ZRCT    ! Cloud water m.r. at t

    REAL, DIMENSION(KSIZE) :: ZRCS    ! Cloud water m.r. source
    REAL, DIMENSION(KSIZE) :: ZRRS    ! Rain water m.r. source
    REAL, DIMENSION(KSIZE) :: ZRIS    ! Pristine ice m.r. source
    REAL, DIMENSION(KSIZE) :: ZRSS    ! Snow/aggregate m.r. source
    REAL, DIMENSION(KSIZE) :: ZRGS    ! Graupel m.r. source
    REAL, DIMENSION(KSIZE) :: ZRHS    ! Hail m.r. source

    REAL, DIMENSION(KSIZE) :: ZRHODREFC ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE) :: ZRHODREFR ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE) :: ZRHODREFI ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE) :: ZRHODREFS ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE) :: ZRHODREFG ! RHO Dry REFerence
    REAL, DIMENSION(KSIZE) :: ZRHODREFH ! RHO Dry REFerence

    REAL, DIMENSION(KSIZE) :: ZCC       ! terminal velocity
    REAL, DIMENSION(KSIZE) :: ZFSEDC1D  ! For cloud sedimentation
    REAL, DIMENSION(KSIZE) :: ZCONC     ! Concentration des aerosols
    REAL, DIMENSION(KSIZE) :: ZRAY1D    ! Mean radius
    REAL, DIMENSION(KSIZE) :: ZWLBDA    ! Libre parcours moyen

    REAL, DIMENSION(D%NIT, D%NKT) :: ZW ! work array

    REAL, DIMENSION(KSIZE) :: ZZT       ! Temperature
    REAL, DIMENSION(KSIZE) :: ZPRES     ! Pressure
    REAL, DIMENSION(KSIZE) :: ZWLBDC    ! Slope parameter of the droplet  distribution

    REAL, DIMENSION(SIZE(ICED%XRTMIN)) :: ZRTMIN

    REAL    :: ZTSPLITR ! Small time step for rain sedimentation
    REAL    :: ZINVTSTEP

    INTEGER :: JN, JL, JK, JI, JJ

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!        O. Initialization of for sedimentation
!
    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SEDIMENTATION_SPLIT',0,ZHOOK_HANDLE)

    ZTSPLITR  = PTSTEP / FLOAT(KSPLITR)
    ZINVTSTEP = 1./PTSTEP
    ZRTMIN(:) = ICED%XRTMIN(:) * ZINVTSTEP

    IF (OSEDIC) PINPRC (:) = 0.
    PINPRR (:) = 0.
    PINPRS (:) = 0.
    PINPRG (:) = 0.
    IF ( KRR == 7 ) PINPRH(:) = 0.
!
!*       1. Parameters for cloud sedimentation
!        Computation moved to beginning of rain_ice
!
!*       2.    compute the fluxes
!
!   optimization by looking for locations where
!   the precipitating fields are larger than a minimal value only !!!
!   For optimization we consider each variable separately

    ZRTMIN(:)    = ICED%XRTMIN(:) * ZINVTSTEP
    IF (OSEDIC) GSEDIMC(:,:) = .FALSE.
    GSEDIMR(:,:) = .FALSE.
    GSEDIMI(:,:) = .FALSE.
    GSEDIMS(:,:) = .FALSE.
    GSEDIMG(:,:) = .FALSE.
    IF ( KRR == 7 ) GSEDIMH(:,:) = .FALSE.
!
! ZPiS = Specie i source creating during the current time step
! PRiS = Source of the previous time step
!
    IF (OSEDIC) THEN
      ZPRCS(:,:) = 0.0
      ZPRCS(:,:) = PRCS(:,:)-PRCT(:,:)* ZINVTSTEP
      PRCS(:,:)  = PRCT(:,:)* ZINVTSTEP
    END IF

    ZPRRS(:,:) = 0.0
    ZPRSS(:,:) = 0.0
    ZPRGS(:,:) = 0.0
    IF ( KRR == 7 ) ZPRHS(:,:) = 0.0
!
    ZPRRS(:,:) = PRRS(:,:)-PRRT(:,:)* ZINVTSTEP
    ZPRSS(:,:) = PRSS(:,:)-PRST(:,:)* ZINVTSTEP
    ZPRGS(:,:) = PRGS(:,:)-PRGT(:,:)* ZINVTSTEP
    IF ( KRR == 7 ) ZPRHS(:,:) = PRHS(:,:)-PRHT(:,:)* ZINVTSTEP
    PRRS(:,:)  = PRRT(:,:)* ZINVTSTEP
    PRSS(:,:)  = PRST(:,:)* ZINVTSTEP
    PRGS(:,:)  = PRGT(:,:)* ZINVTSTEP
    IF ( KRR == 7 ) PRHS(:,:)  = PRHT(:,:)* ZINVTSTEP
!
! PRiS = Source of the previous time step + source created during the subtime
! step
!
    DO JN = 1, KSPLITR
      IF( JN==1 ) THEN
        IF (OSEDIC) PRCS(:,:) = PRCS(:,:) + ZPRCS(:,:)/KSPLITR
          PRRS(:,:) = PRRS(:,:) + ZPRRS(:,:)/KSPLITR
          PRSS(:,:) = PRSS(:,:) + ZPRSS(:,:)/KSPLITR
          PRGS(:,:) = PRGS(:,:) + ZPRGS(:,:)/KSPLITR
        IF ( KRR == 7 ) PRHS(:,:) = PRHS(:,:) + ZPRHS(:,:)/KSPLITR
        DO JK = D%NKTB , D%NKTE
          DO JI = D%NIB , D%NIE
            ZW(JI,JK) =ZTSPLITR/(PRHODREF(JI,JK)* PDZZ(JI,JK))
          END DO
        END DO
      ELSE
        IF (OSEDIC) PRCS(:,:) = PRCS(:,:) + ZPRCS(:,:)*ZTSPLITR
        PRRS(:,:) = PRRS(:,:) + ZPRRS(:,:)*ZTSPLITR
        PRSS(:,:) = PRSS(:,:) + ZPRSS(:,:)*ZTSPLITR
        PRGS(:,:) = PRGS(:,:) + ZPRGS(:,:)*ZTSPLITR
        IF ( KRR == 7 ) PRHS(:,:) = PRHS(:,:) + ZPRHS(:,:)*ZTSPLITR
      END IF
 !
      IF (OSEDIC) GSEDIMC(D%NIB:D%NIE,D%NKTB:D%NKTE) =                &
                      PRCS(D%NIB:D%NIE,D%NKTB:D%NKTE)>ZRTMIN(2)
      GSEDIMR(D%NIB:D%NIE,D%NKTB:D%NKTE) =                            &
          PRRS(D%NIB:D%NIE,D%NKTB:D%NKTE)>ZRTMIN(3)
      GSEDIMI(D%NIB:D%NIE,D%NKTB:D%NKTE) =                            &
          PRIS(D%NIB:D%NIE,D%NKTB:D%NKTE)>ZRTMIN(4)
      GSEDIMS(D%NIB:D%NIE,D%NKTB:D%NKTE) =                            &
          PRSS(D%NIB:D%NIE,D%NKTB:D%NKTE)>ZRTMIN(5)
      GSEDIMG(D%NIB:D%NIE,D%NKTB:D%NKTE) =                            &
          PRGS(D%NIB:D%NIE,D%NKTB:D%NKTE)>ZRTMIN(6)
      IF ( KRR == 7 ) GSEDIMH(D%NIB:D%NIE,D%NKTB:D%NKTE) =            &
                          PRHS(D%NIB:D%NIE,D%NKTB:D%NKTE)>ZRTMIN(7)
!
      IF (OSEDIC) CALL  COUNTJV(ISEDIMC, GSEDIMC(:,:),IC1(:),IC2(:))
      CALL COUNTJV(ISEDIMR, GSEDIMR(:,:),IR1(:),IR2(:))
      CALL COUNTJV(ISEDIMI, GSEDIMI(:,:),II1(:),II2(:))
      CALL COUNTJV(ISEDIMS, GSEDIMS(:,:),IS1(:),IS2(:))
      CALL COUNTJV(ISEDIMG, GSEDIMG(:,:),IG1(:),IG2(:))
      IF ( KRR == 7 ) CALL COUNTJV(ISEDIMH, GSEDIMH(:,:),IH1(:),IH2(:))
!
!*       2.1   for cloud
!
      IF (OSEDIC) THEN

        ZWSED(:,:) = 0.
        IF( JN==1 ) PRCS(:,:) = PRCS(:,:) * PTSTEP

        IF( ISEDIMC >= 1 ) THEN

          DO JL=1,ISEDIMC
            ZRCS(JL) = PRCS(IC1(JL),IC2(JL))
            ZRHODREFC(JL) =  PRHODREF(IC1(JL),IC2(JL))
            ZWLBDC(JL) = ZLBC(IC1(JL),IC2(JL))
            ZCONC(JL) = ZCONC3D(IC1(JL),IC2(JL))
            ZRCT(JL) = PRCT(IC1(JL),IC2(JL))
            ZZT(JL) = PTHT(IC1(JL),IC2(JL))
            ZPRES(JL) = PPABST(IC1(JL),IC2(JL))
            ZRAY1D(JL) = ZRAY(IC1(JL),IC2(JL))
            ZFSEDC1D(JL) = ZFSEDC(IC1(JL),IC2(JL))
          END DO

          ILISTLENC = 0
          DO JL=1,ISEDIMC
            IF( ZRCS(JL) .GT. ZRTMIN(2) ) THEN
              ILISTLENC = ILISTLENC + 1
              ILISTC(ILISTLENC) = JL
            END IF
          END DO

          DO JJ = 1, ILISTLENC
            JL = ILISTC(JJ)
            IF (ZRCS(JL) .GT. ZRTMIN(2) .AND. ZRCT(JL) .GT. ICED%XRTMIN(2)) THEN
              ZWLBDC(JL) = ZWLBDC(JL) * ZCONC(JL) / (ZRHODREFC(JL) * ZRCT(JL))
              ZWLBDC(JL) = ZWLBDC(JL)**ICED%XLBEXC
              ZRAY1D(JL) = ZRAY1D(JL) / ZWLBDC(JL) !! ZRAY : mean diameter=M(1)/2
              ZZT(JL)    = ZZT(JL) * (ZPRES(JL)/CST%XP00)**(CST%XRD/CST%XCPD)
              ZWLBDA(JL) = 6.6E-8*(101325./ZPRES(JL))*(ZZT(JL)/293.15)
              ZCC(JL)    = ICED%XCC*(1.+1.26*ZWLBDA(JL)/ZRAY1D(JL)) !! XCC modified for cloud
              ZWSED(IC1(JL),IC2(JL))= ZRHODREFC(JL)**(-ICED%XCEXVT +1 ) *   &
              ZWLBDC(JL)**(-ICED%XDC)*ZCC(JL)*ZFSEDC1D(JL) * ZRCS(JL)
            END IF
          END DO

        END IF

        DO JK = D%NKTB , D%NKTE
          PRCS(:,JK) = PRCS(:,JK) + ZW(:,JK)*(ZWSED(:,JK+KKL)-ZWSED(:,JK))
        END DO

        IF (PRESENT(PFPR)) THEN
          DO JK = D%NKTB , D%NKTE
            PFPR(:,JK,2)=ZWSED(:,JK)
          ENDDO
        ENDIF

        PINPRC(:) = PINPRC(:) + ZWSED(:,IKB) / CST%XRHOLW / KSPLITR

        IF( JN==KSPLITR ) THEN
          PRCS(:,:) = PRCS(:,:) * ZINVTSTEP
        END IF

      END IF !OSEDIC
!
!*       2.2   for rain
!
      IF( JN==1 ) PRRS(:,:) = PRRS(:,:) * PTSTEP
      ZWSED(:,:) = 0.

      IF( ISEDIMR >= 1 ) THEN
!
        DO JL=1,ISEDIMR
          ZRRS(JL) = PRRS(IR1(JL),IR2(JL))
          ZRHODREFR(JL) =  PRHODREF(IR1(JL),IR2(JL))
        END DO
!
        ILISTLENR = 0
        DO JL=1,ISEDIMR
          IF( ZRRS(JL) .GT. ZRTMIN(3) ) THEN
            ILISTLENR = ILISTLENR + 1
            ILISTR(ILISTLENR) = JL
          END IF
        END DO

        DO JJ = 1, ILISTLENR
          JL = ILISTR(JJ)
          ZWSED(IR1(JL),IR2(JL))= ICEP%XFSEDR  * ZRRS(JL)**ICEP%XEXSEDR *   &
                                   ZRHODREFR(JL)**(ICEP%XEXSEDR-ICED%XCEXVT)
        END DO
      END IF ! ISEDIMR

      DO JK = D%NKTB , D%NKTE
        PRRS(:,JK) = PRRS(:,JK) + ZW(:,JK)*(ZWSED(:,JK+KKL)-ZWSED(:,JK))
      END DO

      IF (PRESENT(PFPR)) THEN
        DO JK = D%NKTB , D%NKTE
          PFPR(:,JK,3)=ZWSED(:,JK)
        ENDDO
      ENDIF

      PINPRR(:) = PINPRR(:) + ZWSED(:,IKB)/CST%XRHOLW/KSPLITR
      IF( JN==KSPLITR ) THEN
        PRRS(:,:) = PRRS(:,:) * ZINVTSTEP
      END IF
!
!*       2.3   for pristine ice
!

      IF( JN==1 ) PRIS(:,:) = PRIS(:,:) * PTSTEP

      ZWSED(:,:) = 0.
      IF( ISEDIMI >= 1 ) THEN

        DO JL=1,ISEDIMI
          ZRIS(JL) = PRIS(II1(JL),II2(JL))
          ZRHODREFI(JL) =  PRHODREF(II1(JL),II2(JL))
        END DO

        ILISTLENI = 0
        DO JL=1,ISEDIMI
          IF( ZRIS(JL) .GT.  MAX(ZRTMIN(4),1.0E-7 )) THEN ! limitation of the McF&H formula
            ILISTLENI = ILISTLENI + 1
            ILISTI(ILISTLENI) = JL
          END IF
        END DO

        DO JJ = 1, ILISTLENI
          JL = ILISTI(JJ)
          ZWSED(II1(JL),II2(JL))= ICEP%XFSEDI * ZRIS(JL) *  &
                          ZRHODREFI(JL)**(1.0-ICED%XCEXVT) * & !    McF&H
                          MAX( 0.05E6,-0.15319E6-0.021454E6* &
                          ALOG(ZRHODREFI(JL)*ZRIS(JL)) )**ICEP%XEXCSEDI
        END DO
      END IF !ISEDIMI

      DO JK = D%NKTB , D%NKTE
        PRIS(:,JK) = PRIS(:,JK) + ZW(:,JK)*(ZWSED(:,JK+KKL)-ZWSED(:,JK))
      END DO

      IF (PRESENT(PFPR)) THEN
        DO JK = D%NKTB , D%NKTE
          PFPR(:,JK,4)=ZWSED(:,JK)
        ENDDO
      ENDIF

      IF( JN==KSPLITR ) THEN
        PRIS(:,:) = PRIS(:,:) * ZINVTSTEP
      END IF
!
!*       2.4   for aggregates/snow
!
      IF( JN==1 ) PRSS(:,:) = PRSS(:,:) * PTSTEP

      ZWSED(:,:) = 0.
      IF( ISEDIMS >= 1 ) THEN
!
        DO JL=1,ISEDIMS
          ZRSS(JL) = PRSS(IS1(JL),IS2(JL))
          ZRHODREFS(JL) =  PRHODREF(IS1(JL),IS2(JL))
        END DO
!
        ILISTLENS = 0
        DO JL=1,ISEDIMS
          IF( ZRSS(JL) .GT. ZRTMIN(5) ) THEN
            ILISTLENS = ILISTLENS + 1
            ILISTS(ILISTLENS) = JL
          END IF
        END DO

        DO JJ = 1, ILISTLENS
          JL = ILISTS(JJ)
          ZWSED(IS1(JL),IS2(JL))= ICEP%XFSEDS * ZRSS(JL)**ICEP%XEXSEDS *  &
                                   ZRHODREFS(JL)**(ICEP%XEXSEDS-ICED%XCEXVT)
        END DO
      END IF !ISEDIMS

      DO JK = D%NKTB , D%NKTE
        PRSS(:,JK) = PRSS(:,JK) + ZW(:,JK)*(ZWSED(:,JK+KKL)-ZWSED(:,JK))
      END DO

      IF (PRESENT(PFPR)) THEN
        DO JK = D%NKTB , D%NKTE
          PFPR(:,JK,5)=ZWSED(:,JK)
        ENDDO
      ENDIF

      PINPRS(:) = PINPRS(:) + ZWSED(:,IKB)/CST%XRHOLW/KSPLITR
      IF( JN==KSPLITR ) THEN
        PRSS(:,:) = PRSS(:,:) * ZINVTSTEP
      END IF
!
!*       2.5   for graupeln
!
      ZWSED(:,:) = 0.
      IF( JN==1 ) PRGS(:,:) = PRGS(:,:) * PTSTEP

      IF( ISEDIMG >= 1 ) THEN
!
        DO JL=1,ISEDIMG
          ZRGS(JL) = PRGS(IG1(JL),IG2(JL))
          ZRHODREFG(JL) =  PRHODREF(IG1(JL),IG2(JL))
        END DO
!
        ILISTLENG = 0
        DO JL=1,ISEDIMG
          IF( ZRGS(JL) .GT. ZRTMIN(6) ) THEN
            ILISTLENG = ILISTLENG + 1
            ILISTG(ILISTLENG) = JL
          END IF
        END DO

        DO JJ = 1, ILISTLENG
          JL = ILISTG(JJ)
          ZWSED(IG1(JL),IG2(JL)) = ICEP%XFSEDG  * ZRGS(JL)**ICEP%XEXSEDG *   &
                                  ZRHODREFG(JL)**(ICEP%XEXSEDG-ICED%XCEXVT)
        END DO
      END IF !ISEDIMG

      DO JK = D%NKTB , D%NKTE
        PRGS(:,JK) = PRGS(:,JK) + ZW(:,JK)*(ZWSED(:,JK+KKL)-ZWSED(:,JK))
      END DO

      IF (PRESENT(PFPR)) THEN
        DO JK = D%NKTB , D%NKTE
          PFPR(:,JK,6)=ZWSED(:,JK)
        ENDDO
      ENDIF

      PINPRG(:) = PINPRG(:) + ZWSED(:,IKB)/CST%XRHOLW/KSPLITR
      IF( JN==KSPLITR ) THEN
        PRGS(:,:) = PRGS(:,:) * ZINVTSTEP
      END IF
!
!*       2.6   for hail
!
      IF ( KRR == 7 ) THEN
        IF( JN==1 ) PRHS(:,:) = PRHS(:,:) * PTSTEP
        ZWSED(:,:) = 0.

        IF( ISEDIMH >= 1 ) THEN

          DO JL=1,ISEDIMH
            ZRHS(JL) = PRHS(IH1(JL),IH2(JL))
            ZRHODREFH(JL) =  PRHODREF(IH1(JL),IH2(JL))
          END DO

          ILISTLENH = 0
          DO JL=1,ISEDIMH
            IF( ZRHS(JL) .GT. ZRTMIN(7) ) THEN
              ILISTLENH = ILISTLENH + 1
              ILISTH(ILISTLENH) = JL
            END IF
          END DO

          DO JJ = 1, ILISTLENH
            JL = ILISTH(JJ)
            ZWSED(IH1(JL),IH2(JL))= ICEP%XFSEDH  * ZRHS(JL)**ICEP%XEXSEDH *   &
                                     ZRHODREFH(JL)**(ICEP%XEXSEDH-ICED%XCEXVT)
          END DO

        END IF !ISEDIMH

        DO JK = D%NKTB , D%NKTE
          PRHS(:,JK) = PRHS(:,JK) + ZW(:,JK)*(ZWSED(:,JK+KKL)-ZWSED(:,JK))
        END DO

        IF (PRESENT(PFPR)) THEN
          DO JK = D%NKTB , D%NKTE
            PFPR(:,JK,7)=ZWSED(:,JK)
          ENDDO
        ENDIF

        PINPRH(:) = PINPRH(:) + ZWSED(:,IKB)/CST%XRHOLW/KSPLITR
        IF( JN==KSPLITR ) THEN
          PRHS(:,:) = PRHS(:,:) * ZINVTSTEP
        END IF
      END IF !KRR == 7

    END DO !JN = 1, KSPLITR

    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SEDIMENTATION_SPLIT',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_OLD_SEDIMENTATION_SPLIT

  SUBROUTINE COUNTJV(IC, LTAB, I1, I2)

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: IC
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: LTAB ! Mask
    INTEGER, DIMENSION(:), INTENT(OUT) :: I1,I2 ! Used to replace the COUNT and PACK
    INTEGER :: JI,JK

    IC = 0
    DO JK = 1,SIZE(LTAB,2)
      DO JI = 1,SIZE(LTAB,1)
        IF(LTAB(JI,JK)) THEN
          IC = IC +1
          I1(IC) = JI
          I2(IC) = JK
        END IF
      END DO
    END DO

  END SUBROUTINE COUNTJV

END MODULE MODE_RAIN_ICE_OLD_SEDIMENTATION_SPLIT
