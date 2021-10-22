SUBROUTINE ICE4_SEDIMENTATION_SPLIT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, &
                                   &PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PSEA, PTOWN, PINPRH, PRHT, PRHS, PFPR)
!!
!!**  PURPOSE
!!    -------
!!      Computes the sedimentation
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the plitting of rain_ice source code (nov. 2014)
!!                and modified for optimisation
!!
!!    MODIFICATIONS
!!    -------------
!!
!
!
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_CST
USE MODD_RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM
USE MODD_PARAM_ICE
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT
INTEGER,                      INTENT(IN)              :: KKL     !vert. levels type 1=MNH -1=ARO
REAL,                         INTENT(IN)              :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                      INTENT(IN)              :: KRR     ! Number of moist variable
LOGICAL,                      INTENT(IN)              :: OSEDIC  ! Switch for droplet sedim.
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRHODREF! Reference density
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PPABST  ! absolute pressure at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PTHT    ! Theta at time t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRGS    ! Graupel m.r. source
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRR  ! Rain instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRI  ! Pristine ice instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRS  ! Snow instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(KIT,KJT),         OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(KIT,KJT),         OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(KIT,KJT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
!
LOGICAL, DIMENSION(KIT, KJT, KKT) &
    :: GSEDIM ! Test where to compute the SED processes
INTEGER , DIMENSION(SIZE(GSEDIM)) :: I1,I2,I3 ! Used to replace the COUNT

REAL, DIMENSION(KIT, KJT, KKT) :: ZCONC3D, & !  droplet condensation
                                & ZRAY,   & ! Cloud Mean radius
                                & ZLBC,   & ! XLBC weighted by sea fraction
                                & ZFSEDC, &
                                & ZPRCS,ZPRRS,ZPRIS,ZPRSS,ZPRGS,ZPRHS, &   ! Mixing ratios created during the time step
                                & ZW, & ! work array
                                & ZRCT, &
                                & ZRRT, &
                                & ZRIT, &
                                & ZRST, &
                                & ZRGT, &
                                & ZRHT
REAL,    DIMENSION(KIT, KJT,0:KKT+1) :: ZWSED        ! sedimentation fluxes
REAL,    DIMENSION(KIT, KJT) :: ZCONC_TMP    ! Weighted concentration
REAL,    DIMENSION(KIT, KJT) :: ZREMAINT ! Remaining time until the timestep end
REAL :: ZINVTSTEP
INTEGER :: ISEDIM ! ! Case number of sedimentation
INTEGER :: JK
REAL, DIMENSION(SIZE(XRTMIN)) :: ZRSMIN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT', 0, ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!
!        O. Initialization of for sedimentation
!
ZINVTSTEP=1./PTSTEP
ZRSMIN(:) = XRTMIN(:) * ZINVTSTEP
IF (OSEDIC) PINPRC (:,:) = 0.
PINPRR (:,:) = 0.
PINPRI (:,:) = 0.
PINPRS (:,:) = 0.
PINPRG (:,:) = 0.
IF ( KRR == 7 ) PINPRH (:,:) = 0.
IF (PRESENT(PFPR)) PFPR(:,:,:,:) = 0.
!
!*       1. Parameters for cloud sedimentation
!
IF (OSEDIC) THEN
  IF(PRESENT(PSEA) .AND. PRESENT(PTOWN)) THEN
    ZRAY(:,:,:)   = 0.
    ZCONC_TMP(:,:)=PSEA(:,:)*XCONC_SEA+(1.-PSEA(:,:))*XCONC_LAND
    DO JK=KKTB, KKTE
      ZLBC(:,:,JK)   = PSEA(:,:)*XLBC(2)+(1.-PSEA(:,:))*XLBC(1)
      ZFSEDC(:,:,JK) = (PSEA(:,:)*XFSEDC(2)+(1.-PSEA(:,:))*XFSEDC(1))
      ZFSEDC(:,:,JK) = MAX(MIN(XFSEDC(1),XFSEDC(2)),ZFSEDC(:,:,JK))
      ZCONC3D(:,:,JK)= (1.-PTOWN(:,:))*ZCONC_TMP(:,:)+PTOWN(:,:)*XCONC_URBAN
      ZRAY(:,:,JK)   = 0.5*((1.-PSEA(:,:))*GAMMA(XNUC+1.0/XALPHAC)/(GAMMA(XNUC)) + &
              PSEA(:,:)*GAMMA(XNUC2+1.0/XALPHAC2)/(GAMMA(XNUC2)))
    END DO
  ELSE
    ZLBC(:,:,:)   = XLBC(1)
    ZFSEDC(:,:,:) = XFSEDC(1)
    ZCONC3D(:,:,:)= XCONC_LAND
    ZRAY(:,:,:)  = 0.5*(GAMMA(XNUC+1.0/XALPHAC)/(GAMMA(XNUC)))
  ENDIF
  ZRAY(:,:,:)      = MAX(1.,ZRAY(:,:,:))
  ZLBC(:,:,:)      = MAX(MIN(XLBC(1),XLBC(2)),ZLBC(:,:,:))
ENDIF
!
!*       2.    compute the fluxes
!
!  optimization by looking for locations where
!  the precipitating fields are larger than a minimal value only !!!
!  For optimization we consider each variable separately
!
! External tendecies
IF (OSEDIC) ZPRCS(:,:,:) = PRCS(:,:,:)-PRCT(:,:,:)*ZINVTSTEP
ZPRRS(:,:,:) = PRRS(:,:,:)-PRRT(:,:,:)*ZINVTSTEP
ZPRIS(:,:,:) = PRIS(:,:,:)-PRIT(:,:,:)*ZINVTSTEP
ZPRSS(:,:,:) = PRSS(:,:,:)-PRST(:,:,:)*ZINVTSTEP
ZPRGS(:,:,:) = PRGS(:,:,:)-PRGT(:,:,:)*ZINVTSTEP
IF ( KRR == 7 ) ZPRHS(:,:,:) = PRHS(:,:,:)-PRHT(:,:,:)*ZINVTSTEP
!
! mr values inside the time-splitting loop
ZRCT(:,:,:) = PRCT(:,:,:)
ZRRT(:,:,:) = PRRT(:,:,:)
ZRIT(:,:,:) = PRIT(:,:,:)
ZRST(:,:,:) = PRST(:,:,:)
ZRGT(:,:,:) = PRGT(:,:,:)
IF (KRR==7) ZRHT(:,:,:) = PRHT(:,:,:)
!
DO JK = KKTB , KKTE
  ZW(:,:,JK) =1./(PRHODREF(:,:,JK)* PDZZ(:,:,JK))
END DO
!
!
!*       2.1   for cloud
!
IF (OSEDIC) THEN
  ZREMAINT(:,:) = PTSTEP
  DO WHILE (ANY(ZREMAINT>0.))
    GSEDIM(:,:,:)=.FALSE.
    DO JK = KKTB , KKTE
      GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                    (ZRCT(KIB:KIE,KJB:KJE,JK)>XRTMIN(2) .OR.    &
                     ZPRCS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(2)) .AND. &
                    ZREMAINT(KIB:KIE,KJB:KJE)>0.
    ENDDO
    ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                             &SIZE(I1),I1(:),I2(:),I3(:))
    CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &2, &
                          &ZRCT, PRCS, ZWSED, PINPRC, ZPRCS, &
                          &ZRAY, ZLBC, ZFSEDC, ZCONC3D, PFPR=PFPR)
  ENDDO
ENDIF
!
!*       2.2   for rain
!
ZREMAINT(:,:) = PTSTEP
DO WHILE (ANY(ZREMAINT>0.))
  GSEDIM(:,:,:)=.FALSE.
  DO JK = KKTB , KKTE
    GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                  (ZRRT(KIB:KIE,KJB:KJE,JK)>XRTMIN(3) .OR.    &
                   ZPRRS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(3)) .AND. &
                  ZREMAINT(KIB:KIE,KJB:KJE)>0.
  ENDDO
  ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                           &SIZE(I1),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &3, &
                          &ZRRT, PRRS, ZWSED, PINPRR, ZPRRS, &
                          &PFPR=PFPR)
ENDDO
!
!*       2.3   for pristine ice
!
ZREMAINT(:,:) = PTSTEP
DO WHILE (ANY(ZREMAINT>0.))
  GSEDIM(:,:,:)=.FALSE.
  DO JK = KKTB , KKTE
    GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                  (ZRIT(KIB:KIE,KJB:KJE,JK)>XRTMIN(4) .OR.    &
                   ZPRIS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(4)) .AND. &
                  ZREMAINT(KIB:KIE,KJB:KJE)>0.
  ENDDO
  ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                           &SIZE(I1),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &4, &
                          &ZRIT, PRIS, ZWSED, PINPRI, ZPRIS, PFPR=PFPR)
ENDDO
!
!*       2.4   for aggregates/snow
!
ZREMAINT(:,:) = PTSTEP
DO WHILE (ANY(ZREMAINT>0.))
  GSEDIM(:,:,:)=.FALSE.
  DO JK = KKTB , KKTE
    GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                  (ZRST(KIB:KIE,KJB:KJE,JK)>XRTMIN(5) .OR.    &
                   ZPRSS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(5)) .AND. &
                  ZREMAINT(KIB:KIE,KJB:KJE)>0.
  ENDDO
  ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                           &SIZE(I1),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &5, &
                          &ZRST, PRSS, ZWSED, PINPRS, ZPRSS, PFPR=PFPR)
ENDDO
!
!*       2.5   for graupeln
!
ZREMAINT(:,:) = PTSTEP
DO WHILE (ANY(ZREMAINT>0.))
  GSEDIM(:,:,:)=.FALSE.
  DO JK = KKTB , KKTE
    GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                  (ZRGT(KIB:KIE,KJB:KJE,JK)>XRTMIN(6) .OR.    &
                   ZPRGS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(6)) .AND. &
                  ZREMAINT(KIB:KIE,KJB:KJE)>0.
  ENDDO
  ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                           &SIZE(I1),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &6, &
                          &ZRGT, PRGS, ZWSED, PINPRG, ZPRGS, PFPR=PFPR)
ENDDO
!
!*       2.6   for hail
!
IF (KRR==7) THEN
  ZREMAINT(:,:) = PTSTEP
  DO WHILE (ANY(ZREMAINT>0.))
    GSEDIM(:,:,:)=.FALSE.
    DO JK = KKTB , KKTE
      GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                    (ZRHT(KIB:KIE,KJB:KJE,JK)>XRTMIN(7) .OR.    &
                     ZPRHS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(7)) .AND. &
                    ZREMAINT(KIB:KIE,KJB:KJE)>0.
    ENDDO
    ISEDIM = ICE4_SEDIMENTATION_SPLIT_COUNTJV(GSEDIM(:,:,:),KIT,KJT,KKT,&
                                             &SIZE(I1),I1(:),I2(:),I3(:))
    CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                            &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                            &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                            &7, &
                            &ZRHT, PRHS, ZWSED, PINPRH, ZPRHS, PFPR=PFPR)
  END DO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT', 1, ZHOOK_HANDLE)
!
CONTAINS
!
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                                &KSEDIM, LDSEDIM, I1, I2, I3, PMAXCFL, PREMAINT, &
                                &PRHODREF, POORHODZ, PDZZ, PPABST, PTHT, PTSTEP, &
                                &KSPE, &
                                &PRXT, PRXS, PWSED, PINPRX, PPRXS, &
                                &PRAY, PLBC, PFSEDC, PCONC3D, PFPR)
    !
    !*      0. DECLARATIONS
    !          ------------
    !
    USE MODD_RAIN_ICE_DESCR
    USE MODD_RAIN_ICE_PARAM
    !
    IMPLICIT NONE
    !
    !*       0.1   Declarations of dummy arguments :
    !
    INTEGER, INTENT(IN) :: KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR
    INTEGER, INTENT(IN) :: KSEDIM
    LOGICAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: LDSEDIM
    INTEGER, DIMENSION(KSEDIM), INTENT(IN) :: I1, I2, I3
    REAL,                         INTENT(IN)              :: PMAXCFL ! maximum CFL allowed
    REAL, DIMENSION(KIT,KJT),     INTENT(INOUT)           :: PREMAINT ! Time remaining until the end of the timestep
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRHODREF ! Reference density
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: POORHODZ ! One Over (Rhodref times delta Z)
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PDZZ ! layer thikness (m)
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PPABST
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PTHT
    REAL,                         INTENT(IN)              :: PTSTEP  ! total timestep
    INTEGER,                      INTENT(IN)              :: KSPE ! 1 for rc, 2 for rr...
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRXT ! mr of specy X
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRXS !Tendency of the specy KSPE
    REAL, DIMENSION(KIT,KJT,0:KKT+1), INTENT(OUT)         :: PWSED ! sedimentation flux
    REAL, DIMENSION(KIT,KJT),     INTENT(INOUT)           :: PINPRX ! instant precip
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PPRXS ! external tendencie
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN), OPTIONAL    :: PRAY, PLBC, PFSEDC, PCONC3D
    REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(INOUT)   :: PFPR    ! upper-air precipitation fluxes
    !
    !*       0.2  declaration of local variables
    !
    !
    INTEGER :: JK, JL, JI, JJ
    REAL :: ZINVTSTEP
    REAL :: ZZWLBDC, ZRAY, ZZT, ZZWLBDA, ZZCC
    REAL :: ZFSED, ZEXSED
    REAL, DIMENSION(KIT, KJT) :: ZMRCHANGE
    REAL, DIMENSION(KIT, KJT) :: ZMAX_TSTEP ! Maximum CFL in column
    REAL, DIMENSION(SIZE(XRTMIN)) :: ZRSMIN
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT:INTERNAL_SEDIM_SPLIT', 0, ZHOOK_HANDLE)
    !
    !-------------------------------------------------------------------------------
    !
    !
    !*       1. Parameters for cloud sedimentation
    !
    !
    !*       2.    compute the fluxes
    !
    !
    ZINVTSTEP = 1./PTSTEP
    ZRSMIN(:) = XRTMIN(:) * ZINVTSTEP
    IF(KSPE==2) THEN
      !******* for cloud
      PWSED(:,:,:) = 0.
      DO JL=1, KSEDIM
        JI=I1(JL)
        JJ=I2(JL)
        JK=I3(JL)
        IF(PRXT(JI,JJ,JK)>XRTMIN(KSPE)) THEN
          ZZWLBDC = PLBC(JI,JJ,JK) * PCONC3D(JI,JJ,JK) / &
                   (PRHODREF(JI,JJ,JK) * PRXT(JI,JJ,JK))
          ZZWLBDC = ZZWLBDC**XLBEXC
          ZRAY = PRAY(JI,JJ,JK) / ZZWLBDC
          ZZT = PTHT(JI,JJ,JK) * (PPABST(JI,JJ,JK)/XP00)**(XRD/XCPD)
          ZZWLBDA = 6.6E-8*(101325./PPABST(JI,JJ,JK))*(ZZT/293.15)
          ZZCC = XCC*(1.+1.26*ZZWLBDA/ZRAY)
          PWSED(JI, JJ, JK) = PRHODREF(JI,JJ,JK)**(-XCEXVT +1 ) *   &
                   ZZWLBDC**(-XDC)*ZZCC*PFSEDC(JI,JJ,JK) * PRXT(JI,JJ,JK)
        ENDIF
      ENDDO
    ELSEIF(KSPE==4) THEN
      ! ******* for pristine ice
      PWSED(:,:,:) = 0.
      DO JL=1, KSEDIM
        JI=I1(JL)
        JJ=I2(JL)
        JK=I3(JL)
        IF(PRXT(JI, JJ, JK) .GT. MAX(XRTMIN(4), 1.0E-7)) THEN
          PWSED(JI, JJ, JK) =  XFSEDI * PRXT(JI, JJ, JK) *  &
                              & PRHODREF(JI,JJ,JK)**(1.-XCEXVT) * & !    McF&H
                              & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                              &      ALOG(PRHODREF(JI,JJ,JK)*PRXT(JI,JJ,JK)) )**XEXCSEDI
        ENDIF
      ENDDO
    ELSE
      ! ******* for other species
      IF(KSPE==3) THEN
        ZFSED=XFSEDR
        ZEXSED=XEXSEDR
      ELSEIF(KSPE==5) THEN
        ZFSED=XFSEDS
        ZEXSED=XEXSEDS
      ELSEIF(KSPE==6) THEN
        ZFSED=XFSEDG
        ZEXSED=XEXSEDG
      ELSEIF(KSPE==7) THEN
        ZFSED=XFSEDH
        ZEXSED=XEXSEDH
      ELSE
        WRITE(*,*) ' STOP'
        WRITE(*,*) ' NO SEDIMENTATION PARAMETER FOR KSPE==', KSPE
        CALL ABORT
        STOP
      ENDIF
      PWSED(:,:,:) = 0.
      DO JL=1, KSEDIM
        JI=I1(JL)
        JJ=I2(JL)
        JK=I3(JL)
        IF(PRXT(JI,JJ,JK)>XRTMIN(KSPE)) THEN
          PWSED(JI, JJ, JK) = ZFSED  * PRXT(JI, JJ, JK)**ZEXSED *   &
                                       PRHODREF(JI, JJ, JK)**(ZEXSED-XCEXVT)
        ENDIF
      ENDDO
    ENDIF
    ZMAX_TSTEP(:,:) = PREMAINT(:,:)
    DO JL=1, KSEDIM
      JI=I1(JL)
      JJ=I2(JL)
      JK=I3(JL)
      IF(PRXT(JI,JJ,JK)>XRTMIN(KSPE) .AND. PWSED(JI, JJ, JK)>1.E-20) THEN
        ZMAX_TSTEP(JI, JJ) = MIN(ZMAX_TSTEP(JI, JJ), PMAXCFL * PRHODREF(JI, JJ, JK) * &
                      PRXT(JI, JJ, JK) * PDZZ(JI, JJ, JK) / PWSED(JI, JJ, JK))
      ENDIF
    ENDDO
    ZMRCHANGE(:,:) = 0.
    PREMAINT(:,:) = PREMAINT(:,:) - ZMAX_TSTEP(:,:)
    DO JK = KKTB , KKTE
      ZMRCHANGE(:,:) = ZMAX_TSTEP(:,:) * POORHODZ(:,:,JK)*(PWSED(:,:,JK+KKL)-PWSED(:,:,JK))
      PRXT(:,:,JK) = PRXT(:,:,JK) + ZMRCHANGE(:,:) + PPRXS(:,:,JK) * ZMAX_TSTEP(:,:)
      PRXS(:,:,JK) = PRXS(:,:,JK) + ZMRCHANGE(:,:) * ZINVTSTEP
    ENDDO
    PINPRX(:,:) = PINPRX(:,:) + ZWSED(:,:,KKB) / XRHOLW * (ZMAX_TSTEP(:,:) * ZINVTSTEP)
    IF (PRESENT(PFPR)) THEN
      DO JK = KKTB , KKTE
        PFPR(:,:,JK,KSPE) = PFPR(:,:,JK,KSPE) + ZWSED(:,:,JK) * (ZMAX_TSTEP(:,:) * ZINVTSTEP)
      ENDDO
    ENDIF
    !
    IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT:INTERNAL_SEDIM_SPLIT', 1, ZHOOK_HANDLE)
  END SUBROUTINE INTERNAL_SEDIM_SPLI
  !
  FUNCTION ICE4_SEDIMENTATION_SPLIT_COUNTJV(LTAB,KIT,KJT,KKT,KSIZE,I1,I2,I3) RESULT(IC)
  !
  !*      0. DECLARATIONS
  !          ------------
  !
  IMPLICIT NONE
  !
  !*       0.2  declaration of local variables
  !
  INTEGER, INTENT(IN) :: KIT,KJT,KKT,KSIZE
  LOGICAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)  :: LTAB ! Mask
  INTEGER, DIMENSION(KSIZE), INTENT(OUT) :: I1,I2,I3 ! Used to replace the COUNT and PACK
  INTEGER :: JI,JJ,JK,IC
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  !
  !-------------------------------------------------------------------------------
  !
  IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT:ICE4_SEDIMENTATION_SPLIT_COUNTJV', 0, ZHOOK_HANDLE)
  IC = 0
  DO JK = 1,SIZE(LTAB,3)
    DO JJ = 1,SIZE(LTAB,2)
      DO JI = 1,SIZE(LTAB,1)
        IF( LTAB(JI,JJ,JK) ) THEN
          IC = IC +1
          I1(IC) = JI
          I2(IC) = JJ
          I3(IC) = JK
        END IF
      END DO
    END DO
  END DO
  !
  IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT:ICE4_SEDIMENTATION_SPLIT_COUNTJV', 1, ZHOOK_HANDLE)
  END FUNCTION ICE4_SEDIMENTATION_SPLIT_COUNTJV
  !
END SUBROUTINE ICE4_SEDIMENTATION_SPLIT
