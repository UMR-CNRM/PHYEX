SUBROUTINE ICE4_SEDIMENTATION_STAT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, &
                                  &PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, &
                                  &PRSS, PRST, PRGS, PRGT,&
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
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Ryad El Khatib 09-Oct-2019 Substantial re-write for optimization
!!       (outerunrolling, vectorization, memory cache saving, unrolling)
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
LOGICAL :: LLSEA_AND_TOWN
INTEGER :: JRR, JI, JJ, JK
INTEGER :: ISHIFT, IK, IKPLUS
REAL :: ZQP, ZP1, ZINVTSTEP, ZGAC, ZGC, ZGAC2, ZGC2, ZRAYDEFO
REAL, DIMENSION(KIT)     :: ZWSEDW1, ZWSEDW2 ! sedimentation speed
REAL, DIMENSION(KIT,KJT) :: ZTSORHODZ        ! TimeStep Over (Rhodref times delta Z)
REAL, DIMENSION(KIT,KJT,0:1,2:KRR) :: ZSED   ! sedimentation flux array for each species and for above and current levels
REAL :: FWSED1, FWSED2, PWSEDW, PWSEDWSUP, PINVTSTEP, PTSTEP1, PDZZ1, PRHODREF1, PRXT1

REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZHOOK_HANDLE_PRXS
REAL(KIND=JPRB) :: ZHOOK_HANDLE_PINPRX
!
#ifndef REK
! 5 multiplications + 1 division => cost = 7X
FWSED1(PWSEDW,PTSTEP1,PDZZ1,PRHODREF1,PRXT1,PINVTSTEP)=MIN(1.,PWSEDW*PTSTEP1/PDZZ1 )*PRHODREF1*PDZZ1*PRXT1*PINVTSTEP
#else
! 5 multiplications only => cost = 5X
FWSED1(PWSEDW,PTSTEP1,PDZZ1,PRHODREF1,PRXT1,PINVTSTEP)=MIN(PRHODREF1*PDZZ1*PRXT1*PINVTSTEP,PWSEDW*PRHODREF1*PRXT1)
#endif

FWSED2(PWSEDW,PTSTEP1,PDZZ1,PWSEDWSUP)=MAX(0.,1.-PDZZ1/(PTSTEP1*PWSEDW))*PWSEDWSUP

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.    compute the fluxes
!
ZINVTSTEP = 1./PTSTEP
ZGAC=GAMMA(XNUC+1.0/XALPHAC)
ZGC=GAMMA(XNUC)
ZGAC2=GAMMA(XNUC2+1.0/XALPHAC2)
ZGC2=GAMMA(XNUC2)
ZRAYDEFO=MAX(1.,0.5*(ZGAC/ZGC))
LLSEA_AND_TOWN=PRESENT(PSEA).AND.PRESENT(PTOWN)

!
!*       2.    compute the fluxes
!
! Start shift mechanism:
ISHIFT=0
CALL SHIFT

! Initialize vertical loop
DO JRR=2,KRR
  ZSED(:,:,IKPLUS,JRR) = 0.
ENDDO

! calculation sedimentation flux
DO JK = KKE , KKB, -1*KKL

  DO JJ = KJB, KJE
    DO JI = KIB, KIE
      ZTSORHODZ(JI,JJ) =PTSTEP/(PRHODREF(JI,JJ,JK)*PDZZ(JI,JJ,JK))
    ENDDO
  ENDDO
!
  DO JRR=2,KRR

    IF (JRR==2) THEN

      !******* for cloud
      IF (OSEDIC) THEN
        CALL CLOUD(PRCT(:,:,JK))
      ELSE
        ZSED(:,:,IK,JRR)=0.
      ENDIF

    ELSEIF (JRR==3) THEN

      !*       2.2   for rain
      CALL OTHER_SPECIES(XFSEDR,XEXSEDR,PRRT(:,:,JK))

    ELSEIF (JRR==4) THEN

      CALL PRISTINE_ICE(PRIT(:,:,JK))

    ELSEIF (JRR==5) THEN

      !*       2.4   for aggregates/snow
      CALL OTHER_SPECIES(XFSEDS,XEXSEDS,PRST(:,:,JK))

    ELSEIF (JRR==6) THEN

      !*       2.5   for graupeln
      CALL OTHER_SPECIES(XFSEDG,XEXSEDG,PRGT(:,:,JK))

    ELSEIF (JRR==7) THEN

      !*       2.6   for hail
      IF (PRESENT(PRHT))  THEN
        CALL OTHER_SPECIES(XFSEDH,XEXSEDH,PRHT(:,:,JK))
      ENDIF

    ENDIF

  ENDDO ! JRR

  ! Wrap-up

  IF(PRESENT(PFPR)) THEN
    DO JRR=2,KRR
      PFPR(:,:,JK,JRR)=ZSED(:,:,IK,JRR)
    ENDDO
  ENDIF

  !IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:UPDATE_PRXS',0,ZHOOK_HANDLE_PRXS)
    DO JJ = KJB, KJE
      DO JI = KIB, KIE
        PRCS(JI,JJ,JK) = PRCS(JI,JJ,JK)+ZTSORHODZ(JI,JJ)*(ZSED(JI,JJ,IKPLUS,2)-ZSED(JI,JJ,IK,2))*ZINVTSTEP
        PRRS(JI,JJ,JK) = PRRS(JI,JJ,JK)+ZTSORHODZ(JI,JJ)*(ZSED(JI,JJ,IKPLUS,3)-ZSED(JI,JJ,IK,3))*ZINVTSTEP
        PRIS(JI,JJ,JK) = PRIS(JI,JJ,JK)+ZTSORHODZ(JI,JJ)*(ZSED(JI,JJ,IKPLUS,4)-ZSED(JI,JJ,IK,4))*ZINVTSTEP
        PRSS(JI,JJ,JK) = PRSS(JI,JJ,JK)+ZTSORHODZ(JI,JJ)*(ZSED(JI,JJ,IKPLUS,5)-ZSED(JI,JJ,IK,5))*ZINVTSTEP
        PRGS(JI,JJ,JK) = PRGS(JI,JJ,JK)+ZTSORHODZ(JI,JJ)*(ZSED(JI,JJ,IKPLUS,6)-ZSED(JI,JJ,IK,6))*ZINVTSTEP
        IF (PRESENT(PRHS))  THEN
          PRHS(JI,JJ,JK) = PRHS(JI,JJ,JK)+ZTSORHODZ(JI,JJ)*(ZSED(JI,JJ,IKPLUS,7)-ZSED(JI,JJ,IK,7))*ZINVTSTEP
        ENDIF
      ENDDO
    ENDDO
  !IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:UPDATE_PRXS',1,ZHOOK_HANDLE_PRXS)

  IF (JK==KKB) THEN
    !IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:CP_INSTANT_PRECIP',0,ZHOOK_HANDLE_PINPRX)
    DO JJ = KJB, KJE
      DO JI = KIB, KIE
        PINPRC(JI,JJ) = ZSED(JI,JJ,IK,2)/XRHOLW
        PINPRR(JI,JJ) = ZSED(JI,JJ,IK,3)/XRHOLW
        PINPRI(JI,JJ) = ZSED(JI,JJ,IK,4)/XRHOLW
        PINPRS(JI,JJ) = ZSED(JI,JJ,IK,5)/XRHOLW
        PINPRG(JI,JJ) = ZSED(JI,JJ,IK,6)/XRHOLW
        IF (PRESENT(PINPRH)) THEN
          PINPRH(JI,JJ) = ZSED(JI,JJ,IK,7)/XRHOLW
        ENDIF
      ENDDO
    ENDDO
    !IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:CP_INSTANT_PRECIP',1,ZHOOK_HANDLE_PINPRX)
  ENDIF

  ! shift mechanism : current level now takes the place of previous one
  CALL SHIFT

ENDDO ! JK

IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT',1,ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE CLOUD(PRXT)

    REAL, INTENT(IN)    :: PRXT(KIT,KJT) ! mr of specy X

    REAL :: ZLBC    ! XLBC weighted by sea fraction
    REAL :: ZFSEDC
    REAL :: ZCONC3D ! droplet condensation
    REAL :: ZRAY    ! Cloud Mean radius
    REAL :: ZZWLBDA, ZZWLBDC, ZZCC

    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:CLOUD',0,ZHOOK_HANDLE)

    DO JJ = KJB, KJE
      DO JI = KIB, KIE
        !estimation of q' taking into account incoming ZWSED from previous vertical level
        ZQP=ZSED(JI,JJ,IKPLUS,JRR)*ZTSORHODZ(JI,JJ)
        IF ((PRXT(JI,JJ) > XRTMIN(JRR)) .OR. (ZQP > XRTMIN(JRR))) THEN
          IF (LLSEA_AND_TOWN) THEN
            ZRAY   = MAX(1.,0.5*((1.-PSEA(JI,JJ))*ZGAC/ZGC+PSEA(JI,JJ)*ZGAC2/ZGC2))
            ZLBC   = MAX(MIN(XLBC(1),XLBC(2)),(PSEA(JI,JJ)*XLBC(2)+(1.-PSEA(JI,JJ))*XLBC(1)) )
            ZFSEDC = MAX(MIN(XFSEDC(1),XFSEDC(2)), (PSEA(JI,JJ)*XFSEDC(2)+(1.-PSEA(JI,JJ))*XFSEDC(1)) )
            ZCONC3D= (1.-PTOWN(JI,JJ))*(PSEA(JI,JJ)*XCONC_SEA+(1.-PSEA(JI,JJ))*XCONC_LAND) + PTOWN(JI,JJ)  *XCONC_URBAN
          ELSE
            ZRAY   = ZRAYDEFO
            ZLBC   = XLBC(1)
            ZFSEDC = XFSEDC(1)
            ZCONC3D= XCONC_LAND
          ENDIF
          !calculation of w
          IF(PRXT(JI,JJ) > XRTMIN(JRR)) THEN
            ZZWLBDA=6.6E-8*(101325./PPABST(JI,JJ,JK))*(PTHT(JI,JJ,JK)/293.15)
            ZZWLBDC=(ZLBC*ZCONC3D/(PRHODREF(JI,JJ,JK)*PRXT(JI,JJ)))**XLBEXC
            ZZCC=XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY) !! ZCC  : Fall speed
            ZWSEDW1(JI)=PRHODREF(JI,JJ,JK)**(-XCEXVT ) * ZZWLBDC**(-XDC)*ZZCC*ZFSEDC
          ELSE
            ZWSEDW1(JI)=0.
          ENDIF
          IF ( ZQP > XRTMIN(JRR) ) THEN
            ZZWLBDA=6.6E-8*(101325./PPABST(JI,JJ,JK))*(PTHT(JI,JJ,JK)/293.15)
            ZZWLBDC=(ZLBC*ZCONC3D/(PRHODREF(JI,JJ,JK)*ZQP))**XLBEXC
            ZZCC=XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY) !! ZCC  : Fall speed
            ZWSEDW2(JI)=PRHODREF(JI,JJ,JK)**(-XCEXVT ) * ZZWLBDC**(-XDC)*ZZCC*ZFSEDC
          ELSE
            ZWSEDW2(JI)=0.
          ENDIF
        ELSE
          ZWSEDW1(JI)=0.
          ZWSEDW2(JI)=0.
        ENDIF
!- duplicated code -------------------------------------------------------------------------
        IF (ZWSEDW2(JI) /= 0.) THEN
          ZSED(JI,JJ,IK,JRR)=FWSED1(ZWSEDW1(JI),PTSTEP,PDZZ(JI,JJ,JK),PRHODREF(JI,JJ,JK),PRXT(JI,JJ),ZINVTSTEP) &
           & + FWSED2(ZWSEDW2(JI),PTSTEP,PDZZ(JI,JJ,JK),ZSED(JI,JJ,IKPLUS,JRR))
        ELSE
          ZSED(JI,JJ,IK,JRR)=FWSED1(ZWSEDW1(JI),PTSTEP,PDZZ(JI,JJ,JK),PRHODREF(JI,JJ,JK),PRXT(JI,JJ),ZINVTSTEP)
        ENDIF
!-------------------------------------------------------------------------------------------
      ENDDO
    ENDDO

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:CLOUD',1,ZHOOK_HANDLE)

  END SUBROUTINE CLOUD

  SUBROUTINE PRISTINE_ICE(PRXT)

    REAL, INTENT(IN)    :: PRXT(KIT,KJT) ! mr of specy X

    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:PRISTINE_ICE',0,ZHOOK_HANDLE)

    ! ******* for pristine ice
    DO JJ = KJB, KJE
      DO JI = KIB, KIE
        ZQP=ZSED(JI,JJ,IKPLUS,JRR)*ZTSORHODZ(JI,JJ)
        IF ((PRXT(JI,JJ) > XRTMIN(JRR)) .OR. (ZQP > XRTMIN(JRR))) THEN
          !calculation of w
          IF ( PRXT(JI,JJ) > MAX(XRTMIN(JRR),1.0E-7 ) ) THEN
            ZWSEDW1(JI)= XFSEDI *  &
                              & PRHODREF(JI,JJ,JK)**(-XCEXVT) * & !    McF&H
                              & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                              &      ALOG(PRHODREF(JI,JJ,JK)*PRXT(JI,JJ)) )**XEXCSEDI
          ELSE
            ZWSEDW1(JI)=0.
          ENDIF
          IF ( ZQP > MAX(XRTMIN(JRR),1.0E-7 ) ) THEN
            ZWSEDW2(JI)= XFSEDI *  &
                              & PRHODREF(JI,JJ,JK)**(-XCEXVT) * & !    McF&H
                              & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                              &      ALOG(PRHODREF(JI,JJ,JK)*ZQP) )**XEXCSEDI
          ELSE
            ZWSEDW2(JI)=0.
          ENDIF
        ELSE
          ZWSEDW1(JI)=0.
          ZWSEDW2(JI)=0.
        ENDIF
!- duplicated code -------------------------------------------------------------------------
        IF (ZWSEDW2(JI) /= 0.) THEN
          ZSED(JI,JJ,IK,JRR)=FWSED1(ZWSEDW1(JI),PTSTEP,PDZZ(JI,JJ,JK),PRHODREF(JI,JJ,JK),PRXT(JI,JJ),ZINVTSTEP) &
           & + FWSED2(ZWSEDW2(JI),PTSTEP,PDZZ(JI,JJ,JK),ZSED(JI,JJ,IKPLUS,JRR))
        ELSE
          ZSED(JI,JJ,IK,JRR)=FWSED1(ZWSEDW1(JI),PTSTEP,PDZZ(JI,JJ,JK),PRHODREF(JI,JJ,JK),PRXT(JI,JJ),ZINVTSTEP)
        ENDIF
!-------------------------------------------------------------------------------------------
      ENDDO
    ENDDO

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:PRISTINE_ICE',1,ZHOOK_HANDLE)

  END SUBROUTINE PRISTINE_ICE

  SUBROUTINE OTHER_SPECIES(PFSED,PEXSED,PRXT)

    REAL, INTENT(IN)    :: PFSED
    REAL, INTENT(IN)    :: PEXSED
    REAL, INTENT(IN)    :: PRXT(KIT,KJT) ! mr of specy X

    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:OTHER_SPECIES',0,ZHOOK_HANDLE)

    ! for all but cloud and pristine ice :
    DO JJ = KJB, KJE
      DO JI = KIB, KIE
        ZQP=ZSED(JI,JJ,IKPLUS,JRR)*ZTSORHODZ(JI,JJ)
        IF ((PRXT(JI,JJ) > XRTMIN(JRR)) .OR. (ZQP > XRTMIN(JRR))) THEN
          !calculation of w
          IF ( PRXT(JI,JJ) > XRTMIN(JRR) ) THEN
            ZWSEDW1(JI)= PFSED *PRXT(JI,JJ)**(PEXSED-1)*PRHODREF(JI,JJ,JK)**(PEXSED-XCEXVT-1)
          ELSE
            ZWSEDW1(JI)=0.
          ENDIF
          IF ( ZQP > XRTMIN(JRR) ) THEN
            ZWSEDW2(JI)= PFSED *ZQP**(PEXSED-1)*PRHODREF(JI,JJ,JK)**(PEXSED-XCEXVT-1)
          ELSE
            ZWSEDW2(JI)=0.
          ENDIF
        ELSE
          ZWSEDW1(JI)=0.
          ZWSEDW2(JI)=0.
        ENDIF
!- duplicated code -------------------------------------------------------------------------
        IF (ZWSEDW2(JI) /= 0.) THEN
          ZSED(JI,JJ,IK,JRR)=FWSED1(ZWSEDW1(JI),PTSTEP,PDZZ(JI,JJ,JK),PRHODREF(JI,JJ,JK),PRXT(JI,JJ),ZINVTSTEP) &
           & + FWSED2(ZWSEDW2(JI),PTSTEP,PDZZ(JI,JJ,JK),ZSED(JI,JJ,IKPLUS,JRR))
        ELSE
          ZSED(JI,JJ,IK,JRR)=FWSED1(ZWSEDW1(JI),PTSTEP,PDZZ(JI,JJ,JK),PRHODREF(JI,JJ,JK),PRXT(JI,JJ),ZINVTSTEP)
        ENDIF
!-------------------------------------------------------------------------------------------
      ENDDO
    ENDDO

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:OTHER_SPECIES',1,ZHOOK_HANDLE)

  END SUBROUTINE OTHER_SPECIES

  SUBROUTINE SHIFT

    IKPLUS=ISHIFT
    IK=1-ISHIFT
    ISHIFT=1-ISHIFT

  END SUBROUTINE SHIFT
END SUBROUTINE ICE4_SEDIMENTATION_STAT
