!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_SEDIMENTATION_STAT
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_SEDIMENTATION_STAT(D, CST, ICEP, ICED, &
                                  &PTSTEP, KRR, OSEDIC, PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PLBDAS, &
                                  &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, &
                                  &PRSS, PRST, PRGS, PRGT,&
                                  &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, &
                                  &PSEA, PTOWN, &
                                  &PINPRH, PRHT, PRHS, PFPR)

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
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!!      Ryad El Khatib 09-Oct-2019 Substantial re-write for optimization
!!       (outerunrolling, vectorization, memory cache saving, unrolling)
!  P. Wautelet 21/01/2021: initialize untouched part of PFPR
!  J. Wurtz       03/2022: New snow characteristics with LSNOW_T
!
!
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_CST, ONLY: CST_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE MODI_GAMMA, ONLY: GAMMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t),             INTENT(IN)              :: D       !array dimensions
TYPE(CST_t),                  INTENT(IN)              :: CST
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)              :: ICEP
TYPE(RAIN_ICE_DESCR_t),       INTENT(IN)              :: ICED
REAL,                         INTENT(IN)              :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                      INTENT(IN)              :: KRR     ! Number of moist variable
LOGICAL,                      INTENT(IN)              :: OSEDIC  ! Switch for droplet sedim.
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PRHODREF! Reference density
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PPABST  ! absolute pressure at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PLBDAS  ! lambda parameter for snow
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT)           :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT)           :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT)           :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT)           :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT)           :: PRGS    ! Graupel m.r. source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(OUT)             :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(OUT)             :: PINPRR  ! Rain instant precip
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(OUT)             :: PINPRI  ! Pristine ice instant precip
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(OUT)             :: PINPRS  ! Snow instant precip
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(OUT)             :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(D%NIT,D%NJT),     OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(D%NIT,D%NJT),     OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(D%NIT,D%NJT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
LOGICAL :: LLSEA_AND_TOWN
INTEGER :: JRR, JI, JJ, JK, IKB, IKE,IKL, IIE, IIB, IJB, IJE, IKTB, IKTE
INTEGER :: ISHIFT, IK, IKPLUS
REAL :: ZQP, ZP1, ZINVTSTEP, ZGAC, ZGC, ZGAC2, ZGC2, ZRAYDEFO
REAL, DIMENSION(D%NIT)     :: ZWSEDW1, ZWSEDW2 ! sedimentation speed
REAL, DIMENSION(D%NIT,D%NJT) :: ZTSORHODZ        ! TimeStep Over (Rhodref times delta Z)
REAL, DIMENSION(D%NIT,D%NJT,0:1,2:KRR) :: ZSED   ! sedimentation flux array for each species and for above and current levels
REAL :: FWSED1, FWSED2, PWSEDW, PWSEDWSUP, PINVTSTEP, PTSTEP1, PDZZ1, PRHODREF1, PRXT1

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
#if defined(REPRO48) || defined(REPRO55)
! 5 multiplications + 1 division => cost = 7X
FWSED1(PWSEDW,PTSTEP1,PDZZ1,PRHODREF1,PRXT1,PINVTSTEP)=MIN(1.,PWSEDW*PTSTEP1/PDZZ1 )*PRHODREF1*PDZZ1*PRXT1*PINVTSTEP
#else
! 5 multiplications only => cost = 5X
FWSED1(PWSEDW,PTSTEP1,PDZZ1,PRHODREF1,PRXT1,PINVTSTEP)=MIN(PRHODREF1*PDZZ1*PRXT1*PINVTSTEP,PWSEDW*PRHODREF1*PRXT1)
#endif

FWSED2(PWSEDW,PTSTEP1,PDZZ1,PWSEDWSUP)=MAX(0.,1.-PDZZ1/(PTSTEP1*PWSEDW))*PWSEDWSUP

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT',0,ZHOOK_HANDLE)
!
IKB=D%NKB
IKE=D%NKE
IKL=D%NKL
IIB=D%NIB
IIE=D%NIE
IJB=D%NJB
IJE=D%NJE
IKTB=D%NKTB
IKTE=D%NKTE
!
IF ( PRESENT( PFPR ) ) THEN
 !Set to 0. to avoid undefined values (in files)
 PFPR(:, :, : IKTB, :) = 0.
 PFPR(:, :, IKTE :, :) = 0.
END IF

!-------------------------------------------------------------------------------
!
!*       1.    compute the fluxes
!
ZINVTSTEP = 1./PTSTEP
ZGAC=GAMMA(ICED%XNUC+1.0/ICED%XALPHAC)
ZGC=GAMMA(ICED%XNUC)
ZGAC2=GAMMA(ICED%XNUC2+1.0/ICED%XALPHAC2)
ZGC2=GAMMA(ICED%XNUC2)
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
DO JK = IKE , IKB, -1*IKL

  DO JJ = IJB, IJE
    DO JI = IIB, IIE
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
      CALL OTHER_SPECIES(ICEP%XFSEDR,ICEP%XEXSEDR,PRRT(:,:,JK))

    ELSEIF (JRR==4) THEN

      CALL PRISTINE_ICE(PRIT(:,:,JK))

    ELSEIF (JRR==5) THEN

      !*       2.4   for aggregates/snow
#ifdef REPRO48
      CALL OTHER_SPECIES(ICEP%XFSEDS,ICEP%XEXSEDS,PRST(:,:,JK))
#else
      CALL SNOW(PRST(:,:,JK))
#endif

    ELSEIF (JRR==6) THEN

      !*       2.5   for graupeln
      CALL OTHER_SPECIES(ICEP%XFSEDG,ICEP%XEXSEDG,PRGT(:,:,JK))

    ELSEIF (JRR==7) THEN

      !*       2.6   for hail
      IF (PRESENT(PRHT))  THEN
        CALL OTHER_SPECIES(ICEP%XFSEDH,ICEP%XEXSEDH,PRHT(:,:,JK))
      ENDIF

    ENDIF

  ENDDO ! JRR

  ! Wrap-up

  IF(PRESENT(PFPR)) THEN
    DO JRR=2,KRR
      PFPR(:,:,JK,JRR)=ZSED(:,:,IK,JRR)
    ENDDO
  ENDIF

    DO JJ = IJB, IJE
      DO JI = IIB, IIE
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

  IF (JK==IKB) THEN
    DO JJ = IJB, IJE
      DO JI = IIB, IIE
        IF(OSEDIC) PINPRC(JI,JJ) = ZSED(JI,JJ,IK,2)/CST%XRHOLW
        PINPRR(JI,JJ) = ZSED(JI,JJ,IK,3)/CST%XRHOLW
        PINPRI(JI,JJ) = ZSED(JI,JJ,IK,4)/CST%XRHOLW
        PINPRS(JI,JJ) = ZSED(JI,JJ,IK,5)/CST%XRHOLW
        PINPRG(JI,JJ) = ZSED(JI,JJ,IK,6)/CST%XRHOLW
        IF (PRESENT(PINPRH)) THEN
          PINPRH(JI,JJ) = ZSED(JI,JJ,IK,7)/CST%XRHOLW
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  ! shift mechanism : current level now takes the place of previous one
  CALL SHIFT

ENDDO ! JK

IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT',1,ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE CLOUD(PRXT)

    REAL, INTENT(IN)    :: PRXT(D%NIT,D%NJT) ! mr of specy X

    REAL :: ZLBC    ! XLBC weighted by sea fraction
    REAL :: ZFSEDC
    REAL :: ZCONC3D ! droplet condensation
    REAL :: ZRAY    ! Cloud Mean radius
    REAL :: ZZWLBDA, ZZWLBDC, ZZCC

    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:CLOUD',0,ZHOOK_HANDLE)

    DO JJ = IJB, IJE
      DO JI = IIB, IIE
        !estimation of q' taking into account incoming ZWSED from previous vertical level
        ZQP=ZSED(JI,JJ,IKPLUS,JRR)*ZTSORHODZ(JI,JJ)
        IF ((PRXT(JI,JJ) > ICED%XRTMIN(JRR)) .OR. (ZQP > ICED%XRTMIN(JRR))) THEN
          IF (LLSEA_AND_TOWN) THEN
            ZRAY   = MAX(1.,0.5*((1.-PSEA(JI,JJ))*ZGAC/ZGC+PSEA(JI,JJ)*ZGAC2/ZGC2))
            ZLBC   = MAX(MIN(ICED%XLBC(1),ICED%XLBC(2)),(PSEA(JI,JJ)*ICED%XLBC(2)+(1.-PSEA(JI,JJ))*ICED%XLBC(1)) )
            ZFSEDC = MAX(MIN(ICEP%XFSEDC(1),ICEP%XFSEDC(2)), (PSEA(JI,JJ)*ICEP%XFSEDC(2)+(1.-PSEA(JI,JJ))*ICEP%XFSEDC(1)) )
            ZCONC3D= (1.-PTOWN(JI,JJ))*(PSEA(JI,JJ)*ICED%XCONC_SEA+(1.-PSEA(JI,JJ))*ICED%XCONC_LAND) + &
                      PTOWN(JI,JJ)  *ICED%XCONC_URBAN
          ELSE
            ZRAY   = ZRAYDEFO
            ZLBC   = ICED%XLBC(1)
            ZFSEDC = ICEP%XFSEDC(1)
            ZCONC3D= ICED%XCONC_LAND
          ENDIF
          !calculation of w
          IF(PRXT(JI,JJ) > ICED%XRTMIN(JRR)) THEN
            ZZWLBDA=6.6E-8*(101325./PPABST(JI,JJ,JK))*(PTHT(JI,JJ,JK)/293.15)
            ZZWLBDC=(ZLBC*ZCONC3D/(PRHODREF(JI,JJ,JK)*PRXT(JI,JJ)))**ICED%XLBEXC
            ZZCC=ICED%XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY) !! ZCC  : Fall speed
            ZWSEDW1(JI)=PRHODREF(JI,JJ,JK)**(-ICED%XCEXVT ) * ZZWLBDC**(-ICED%XDC)*ZZCC*ZFSEDC
          ELSE
            ZWSEDW1(JI)=0.
          ENDIF
          IF ( ZQP > ICED%XRTMIN(JRR) ) THEN
            ZZWLBDA=6.6E-8*(101325./PPABST(JI,JJ,JK))*(PTHT(JI,JJ,JK)/293.15)
            ZZWLBDC=(ZLBC*ZCONC3D/(PRHODREF(JI,JJ,JK)*ZQP))**ICED%XLBEXC
            ZZCC=ICED%XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY) !! ZCC  : Fall speed
            ZWSEDW2(JI)=PRHODREF(JI,JJ,JK)**(-ICED%XCEXVT ) * ZZWLBDC**(-ICED%XDC)*ZZCC*ZFSEDC
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

    REAL, INTENT(IN)    :: PRXT(D%NIT,D%NJT) ! mr of specy X

    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:PRISTINE_ICE',0,ZHOOK_HANDLE)

    ! ******* for pristine ice
    DO JJ = IJB, IJE
      DO JI = IIB, IIE
        ZQP=ZSED(JI,JJ,IKPLUS,JRR)*ZTSORHODZ(JI,JJ)
        IF ((PRXT(JI,JJ) > ICED%XRTMIN(JRR)) .OR. (ZQP > ICED%XRTMIN(JRR))) THEN
          !calculation of w
          IF ( PRXT(JI,JJ) > MAX(ICED%XRTMIN(JRR),1.0E-7 ) ) THEN
            ZWSEDW1(JI)= ICEP%XFSEDI *  &
                              & PRHODREF(JI,JJ,JK)**(-ICED%XCEXVT) * & !    McF&H
                              & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                              &      ALOG(PRHODREF(JI,JJ,JK)*PRXT(JI,JJ)) )**ICEP%XEXCSEDI
          ELSE
            ZWSEDW1(JI)=0.
          ENDIF
          IF ( ZQP > MAX(ICED%XRTMIN(JRR),1.0E-7 ) ) THEN
            ZWSEDW2(JI)= ICEP%XFSEDI *  &
                              & PRHODREF(JI,JJ,JK)**(-ICED%XCEXVT) * & !    McF&H
                              & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                              &      ALOG(PRHODREF(JI,JJ,JK)*ZQP) )**ICEP%XEXCSEDI
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

  SUBROUTINE SNOW(PRXT)

    REAL, INTENT(IN)    :: PRXT(D%NIT,D%NJT) ! mr of specy X

    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:SNOW',0,ZHOOK_HANDLE)

    ! ******* for snow
    DO JJ = IJB, IJE
      DO JI = IIB, IIE
        ZQP=ZSED(JI,JJ,IKPLUS,JRR)*ZTSORHODZ(JI,JJ)
        IF ((PRXT(JI,JJ) > ICED%XRTMIN(JRR)) .OR. (ZQP > ICED%XRTMIN(JRR))) THEN
          !calculation of w
          IF ( PRXT(JI,JJ) > ICED%XRTMIN(JRR) ) THEN
            ZWSEDW1(JI)= ICEP%XFSEDS *  &
                          & PRHODREF(JI,JJ,JK)**(-ICED%XCEXVT) * &
                          & (1+(ICED%XFVELOS/PLBDAS(JI,JJ,JK))**ICED%XALPHAS)**(-ICED%XNUS+ICEP%XEXSEDS/ICED%XALPHAS)* &
			   & PLBDAS(JI,JJ,JK)**(ICED%XBS+ICEP%XEXSEDS) 
          ELSE
            ZWSEDW1(JI)=0.
          ENDIF
          IF ( ZQP > ICED%XRTMIN(JRR) ) THEN
            ZWSEDW2(JI)= ICEP%XFSEDS *  &
                          & PRHODREF(JI,JJ,JK)**(-ICED%XCEXVT) * &
                          & (1+(ICED%XFVELOS/PLBDAS(JI,JJ,JK))**ICED%XALPHAS)**(-ICED%XNUS+ICEP%XEXSEDS/ICED%XALPHAS)* &
			   & PLBDAS(JI,JJ,JK)**(ICED%XBS+ICEP%XEXSEDS) 
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

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:SNOW',1,ZHOOK_HANDLE)

  END SUBROUTINE SNOW

  SUBROUTINE OTHER_SPECIES(PFSED,PEXSED,PRXT)

    REAL, INTENT(IN)    :: PFSED
    REAL, INTENT(IN)    :: PEXSED
    REAL, INTENT(IN)    :: PRXT(D%NIT,D%NJT) ! mr of specy X

    REAL(KIND=JPRB) :: ZHOOK_HANDLE

    !!IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_STAT:OTHER_SPECIES',0,ZHOOK_HANDLE)

    ! for all but cloud and pristine ice :
    DO JJ = IJB, IJE
      DO JI = IIB, IIE
        ZQP=ZSED(JI,JJ,IKPLUS,JRR)*ZTSORHODZ(JI,JJ)
        IF ((PRXT(JI,JJ) > ICED%XRTMIN(JRR)) .OR. (ZQP > ICED%XRTMIN(JRR))) THEN
          !calculation of w
          IF ( PRXT(JI,JJ) > ICED%XRTMIN(JRR) ) THEN
            ZWSEDW1(JI)= PFSED *PRXT(JI,JJ)**(PEXSED-1)*PRHODREF(JI,JJ,JK)**(PEXSED-ICED%XCEXVT-1)
          ELSE
            ZWSEDW1(JI)=0.
          ENDIF
          IF ( ZQP > ICED%XRTMIN(JRR) ) THEN
            ZWSEDW2(JI)= PFSED *ZQP**(PEXSED-1)*PRHODREF(JI,JJ,JK)**(PEXSED-ICED%XCEXVT-1)
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
END MODULE MODE_ICE4_SEDIMENTATION_STAT
