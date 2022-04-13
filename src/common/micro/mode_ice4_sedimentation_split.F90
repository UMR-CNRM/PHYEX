!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_SEDIMENTATION_SPLIT
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_SEDIMENTATION_SPLIT(D, CST, ICEP, ICED, PARAMI, &
                                   &PTSTEP, KRR, OSEDIC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PSEA, PTOWN,  &
                                   &PINPRH, PRHT, PRHS, PFPR)
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
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
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
USE MODD_PARAM_ICE,      ONLY: PARAM_ICE_t
!
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
!
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
TYPE(PARAM_ICE_t),            INTENT(IN)              :: PARAMI
REAL,                         INTENT(IN)              :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                      INTENT(IN)              :: KRR     ! Number of moist variable
LOGICAL,                      INTENT(IN)              :: OSEDIC  ! Switch for droplet sedim.
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PRHODREF! Reference density
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PPABST  ! absolute pressure at t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PRHODJ  ! Dry density * Jacobian
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
!
INTEGER                                                             :: JI, JJ, JK
INTEGER                                                             :: IRR !Workaround of PGI bug with OpenACC (at least up to 18.10 version)
LOGICAL                                                             :: GSEDIC !Workaround of PGI bug with OpenACC (at least up to 18.10 version)
LOGICAL                                                             :: GPRESENT_PFPR, GPRESENT_PSEA
REAL                                                                :: ZINVTSTEP
REAL, DIMENSION(D%NIT, D%NJT)                                           :: ZCONC_TMP    ! Weighted concentration
REAL, DIMENSION(D%NIT,D%NJT,D%NKTB:D%NKTE)                                  :: ZW ! work array
REAL, DIMENSION(D%NIT, D%NJT, D%NKT)                                      :: ZCONC3D, & !  droplet condensation
                                                                     & ZRAY,   & ! Cloud Mean radius
                                                                     & ZLBC,   & ! XLBC weighted by sea fraction
                                                                     & ZFSEDC, &
                                                                     & ZPRCS,ZPRRS,ZPRIS,ZPRSS,ZPRGS,ZPRHS, &   ! Mixing ratios created during the time step
                                                                     & ZRCT, &
                                                                     & ZRRT, &
                                                                     & ZRIT, &
                                                                     & ZRST, &
                                                                     & ZRGT, &
                                                                     & ZRHT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT', 0, ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!
GSEDIC = OSEDIC
IRR    = KRR
!
IF (PRESENT(PFPR)) THEN
  GPRESENT_PFPR = .TRUE.
ELSE
  GPRESENT_PFPR = .FALSE.
END IF
!
IF (PRESENT(PSEA)) THEN
  GPRESENT_PSEA = .TRUE.
ELSE
  GPRESENT_PSEA = .FALSE.
END IF
!
!        O. Initialization of for sedimentation
!
ZINVTSTEP=1./PTSTEP
IF (GPRESENT_PFPR) THEN
  PFPR(:,:,:,:) = 0.
END IF
!
!*       1. Parameters for cloud sedimentation
!
IF (GSEDIC) THEN
  ZRAY(:,:,:)   = 0.
  ZLBC(:,:,:)   = ICED%XLBC(1)
  ZFSEDC(:,:,:) = ICEP%XFSEDC(1)
  ZCONC3D(:,:,:)= ICED%XCONC_LAND
  ZCONC_TMP(:,:)= ICED%XCONC_LAND
  IF (GPRESENT_PSEA) THEN
    DO JJ = D%NJB, D%NJE
      DO JI = D%NIB, D%NIE
        ZCONC_TMP(JI,JJ)=PSEA(JI,JJ)*ICED%XCONC_SEA+(1.-PSEA(JI,JJ))*ICED%XCONC_LAND
      ENDDO
    ENDDO
    DO JK=D%NKTB, D%NKTE
      DO JJ = D%NJB, D%NJE
        DO JI = D%NIB, D%NIE
          ZLBC(JI,JJ,JK)   = PSEA(JI,JJ)*ICED%XLBC(2)+(1.-PSEA(JI,JJ))*ICED%XLBC(1)
          ZFSEDC(JI,JJ,JK) = (PSEA(JI,JJ)*ICEP%XFSEDC(2)+(1.-PSEA(JI,JJ))*ICEP%XFSEDC(1))
          ZFSEDC(JI,JJ,JK) = MAX(MIN(ICEP%XFSEDC(1),ICEP%XFSEDC(2)),ZFSEDC(JI,JJ,JK))
          ZCONC3D(JI,JJ,JK)= (1.-PTOWN(JI,JJ))*ZCONC_TMP(JI,JJ)+PTOWN(JI,JJ)*ICED%XCONC_URBAN
          ZRAY(JI,JJ,JK)   = 0.5*((1.-PSEA(JI,JJ))*GAMMA(ICED%XNUC+1.0/ICED%XALPHAC)/(GAMMA(ICED%XNUC)) + &
                         & PSEA(JI,JJ)*GAMMA(ICED%XNUC2+1.0/ICED%XALPHAC2)/(GAMMA(ICED%XNUC2)))
        ENDDO
      ENDDO
    END DO
  ELSE
    ZCONC3D(:,:,:) = ICED%XCONC_LAND
    ZRAY(:,:,:)  = 0.5*(GAMMA(ICED%XNUC+1.0/ICED%XALPHAC)/(GAMMA(ICED%XNUC)))
  END IF
  ZRAY(:,:,:)      = MAX(1.,ZRAY(:,:,:))
  ZLBC(:,:,:)      = MAX(MIN(ICED%XLBC(1),ICED%XLBC(2)),ZLBC(:,:,:))
ENDIF
!
!*       2.    compute the fluxes
!
!  optimization by looking for locations where
!  the precipitating fields are larger than a minimal value only !!!
!  For optimization we consider each variable separately
!
DO JK=D%NKTB, D%NKTE
  DO JJ = D%NJB, D%NJE
    DO JI = D%NIB, D%NIE
      ! External tendecies
      IF (GSEDIC) THEN
        ZPRCS(JI,JJ,JK) = PRCS(JI,JJ,JK)-PRCT(JI,JJ,JK)*ZINVTSTEP
      ENDIF
      ZPRRS(JI,JJ,JK) = PRRS(JI,JJ,JK)-PRRT(JI,JJ,JK)*ZINVTSTEP
      ZPRIS(JI,JJ,JK) = PRIS(JI,JJ,JK)-PRIT(JI,JJ,JK)*ZINVTSTEP
      ZPRSS(JI,JJ,JK) = PRSS(JI,JJ,JK)-PRST(JI,JJ,JK)*ZINVTSTEP
      ZPRGS(JI,JJ,JK) = PRGS(JI,JJ,JK)-PRGT(JI,JJ,JK)*ZINVTSTEP
      IF ( IRR == 7 ) THEN
        ZPRHS(JI,JJ,JK) = PRHS(JI,JJ,JK)-PRHT(JI,JJ,JK)*ZINVTSTEP
      END IF
      !
      ! mr values inside the time-splitting loop
      ZRCT(JI,JJ,JK) = PRCT(JI,JJ,JK)
      ZRRT(JI,JJ,JK) = PRRT(JI,JJ,JK)
      ZRIT(JI,JJ,JK) = PRIT(JI,JJ,JK)
      ZRST(JI,JJ,JK) = PRST(JI,JJ,JK)
      ZRGT(JI,JJ,JK) = PRGT(JI,JJ,JK)
      IF (IRR==7) THEN
        ZRHT(JI,JJ,JK) = PRHT(JI,JJ,JK)
      END IF
      !
      ZW(JI,JJ,JK) =1./(PRHODREF(JI,JJ,JK)* PDZZ(JI,JJ,JK))
    ENDDO
  ENDDO
ENDDO
!
!
!*       2.1   for cloud
!
IF (GSEDIC) THEN
    CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, KRR, &
                          &PARAMI%XSPLIT_MAXCFL, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &2, &
                          &ZRCT, PRCS, PINPRC, ZPRCS, &
                          &ZRAY, ZLBC, ZFSEDC, ZCONC3D, PFPR=PFPR)
ENDIF
!
!*       2.2   for rain
!
  CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, KRR, &
                          &PARAMI%XSPLIT_MAXCFL, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &3, &
                          &ZRRT, PRRS, PINPRR, ZPRRS, &
                          &PFPR=PFPR)
!
!*       2.3   for pristine ice
!
  CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, KRR, &
                          &PARAMI%XSPLIT_MAXCFL, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &4, &
                          &ZRIT, PRIS, PINPRI, ZPRIS, &
                          &PFPR=PFPR)
!
!*       2.4   for aggregates/snow
!
  CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, KRR, &
                        &PARAMI%XSPLIT_MAXCFL, &
                        &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                        &5, &
                        &ZRST, PRSS, PINPRS, ZPRSS, &
                        &PFPR=PFPR)
!
!*       2.5   for graupeln
!
  CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, KRR, &
                        &PARAMI%XSPLIT_MAXCFL, &
                        &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                        &6, &
                        &ZRGT, PRGS, PINPRG, ZPRGS, &
                        &PFPR=PFPR)
!
!*       2.6   for hail
!
IF (IRR==7) THEN
    CALL INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, KRR, &
                          &PARAMI%XSPLIT_MAXCFL, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &7, &
                          &ZRHT, PRHS, PINPRH, ZPRHS, &
                          &PFPR=PFPR)
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
SUBROUTINE INTERNAL_SEDIM_SPLI(D, CST, ICEP, ICED, KRR, &
                              &PMAXCFL, PRHODREF, POORHODZ, PDZZ, PPABST, PTHT, PTSTEP, &
                              &KSPE, PRXT, PRXS, PINPRX, PPRXS, &
                              &PRAY, PLBC, PFSEDC, PCONC3D, PFPR)
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: CST_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t),             INTENT(IN)              :: D
TYPE(CST_t),                  INTENT(IN)              :: CST
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)              :: ICEP
TYPE(RAIN_ICE_DESCR_t),       INTENT(IN)              :: ICED
INTEGER,                      INTENT(IN)              :: KRR
REAL,                         INTENT(IN)              :: PMAXCFL ! maximum CFL allowed
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PRHODREF ! Reference density
REAL, DIMENSION(D%NIT,D%NJT,D%NKTB:D%NKTE), INTENT(IN)        :: POORHODZ ! One Over (Rhodref times delta Z)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PDZZ ! layer thikness (m)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PPABST
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PTHT
REAL,                         INTENT(IN)              :: PTSTEP  ! total timestep
INTEGER,                      INTENT(IN)              :: KSPE ! 1 for rc, 2 for rr...
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT)           :: PRXT ! mr of specy X
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT)           :: PRXS !Tendency of the specy KSPE
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(OUT)             :: PINPRX ! instant precip
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)              :: PPRXS ! external tendencie
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN), OPTIONAL    :: PRAY, PLBC, PFSEDC, PCONC3D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(INOUT), OPTIONAL :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
CHARACTER(LEN=10) :: YSPE ! String for error message
INTEGER                         :: IDX, ISEDIM
INTEGER                         :: JI, JJ, JK, JL
INTEGER, DIMENSION(D%NIT*D%NJT*D%NKT) :: I1,I2,I3   ! Used to replace the COUNT
LOGICAL                         :: GPRESENT_PFPR
REAL                            :: ZINVTSTEP
REAL                            :: ZZWLBDC, ZRAY, ZZT, ZZWLBDA, ZZCC
REAL                            :: ZFSED, ZEXSED
REAL                                :: ZMRCHANGE
REAL, DIMENSION(D%NIT, D%NJT)       :: ZMAX_TSTEP ! Maximum CFL in column
REAL, DIMENSION(SIZE(ICED%XRTMIN))   :: ZRSMIN
REAL, DIMENSION(D%NIT, D%NJT)       :: ZREMAINT   ! Remaining time until the timestep end
REAL, DIMENSION(D%NIT, D%NJT, 0:D%NKT+1) :: ZWSED   ! Sedimentation fluxes
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT:INTERNAL_SEDIM_SPLIT', 0, ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
IF (KSPE<2 .OR. KSPE>7) CALL PRINT_MSG(NVERB_FATAL,'GEN','INTERNAL_SEDIM_SPLIT','invalid species (KSPE variable)')
!
IF (PRESENT(PFPR)) THEN
  GPRESENT_PFPR = .TRUE.
ELSE
  GPRESENT_PFPR = .FALSE.
END IF
!
PINPRX(:,:) = 0.
ZINVTSTEP=1./PTSTEP
ZRSMIN(:) = ICED%XRTMIN(:) * ZINVTSTEP
ZREMAINT(:,:) = 0.
ZREMAINT(D%NIB:D%NIE,D%NJB:D%NJE) = PTSTEP
!
DO WHILE (ANY(ZREMAINT>0.))
  ISEDIM = 0
  DO JK = D%NKTB,D%NKTE
    DO JJ = D%NJB,D%NJE
      DO JI = D%NIB,D%NIE
        IF( (PRXT (JI,JJ,JK)>ICED%XRTMIN(KSPE) .OR.    &
             PPRXS(JI,JJ,JK)>ZRSMIN(KSPE)) .AND. &
             ZREMAINT(JI,JJ)>0. ) THEN
          ISEDIM = ISEDIM + 1
          IDX = ISEDIM
          I1(IDX) = JI
          I2(IDX) = JJ
          I3(IDX) = JK
        END IF
      END DO
    END DO
  END DO
  !
  !
  !*       1. Parameters for cloud sedimentation
  !
  !
  !*       2.    compute the fluxes
  !
  !
  IF(KSPE==2) THEN
    !******* for cloud
    ZWSED(:,:,:) = 0.
    DO JL=1, ISEDIM
      JI=I1(JL)
      JJ=I2(JL)
      JK=I3(JL)
      IF(PRXT(JI,JJ,JK)>ICED%XRTMIN(KSPE)) THEN
        ZZWLBDC = PLBC(JI,JJ,JK) * PCONC3D(JI,JJ,JK) / &
                 &(PRHODREF(JI,JJ,JK) * PRXT(JI,JJ,JK))
        ZZWLBDC = ZZWLBDC**ICED%XLBEXC
        ZRAY = PRAY(JI,JJ,JK) / ZZWLBDC
        ZZT = PTHT(JI,JJ,JK) * (PPABST(JI,JJ,JK)/CST%XP00)**(CST%XRD/CST%XCPD)
        ZZWLBDA = 6.6E-8*(101325./PPABST(JI,JJ,JK))*(ZZT/293.15)
        ZZCC = ICED%XCC*(1.+1.26*ZZWLBDA/ZRAY)
        ZWSED(JI, JJ, JK) = PRHODREF(JI,JJ,JK)**(-ICED%XCEXVT +1 ) *   &
                           &ZZWLBDC**(-ICED%XDC)*ZZCC*PFSEDC(JI,JJ,JK) * PRXT(JI,JJ,JK)
      ENDIF
    ENDDO
  ELSEIF(KSPE==4) THEN
    ! ******* for pristine ice
    ZWSED(:,:,:) = 0.
    DO JL=1, ISEDIM
      JI=I1(JL)
      JJ=I2(JL)
      JK=I3(JL)
      IF(PRXT(JI, JJ, JK) .GT. MAX(ICED%XRTMIN(4), 1.0E-7)) THEN
        ZWSED(JI, JJ, JK) =  ICEP%XFSEDI * PRXT(JI, JJ, JK) *  &
                            & PRHODREF(JI,JJ,JK)**(1.-ICED%XCEXVT) * & !    McF&H
                            & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                            &      ALOG(PRHODREF(JI,JJ,JK)*PRXT(JI,JJ,JK)) )**ICEP%XEXCSEDI
      ENDIF
    ENDDO
  ELSE
    ! ******* for other species
    SELECT CASE(KSPE)
      CASE(3)
        ZFSED=ICEP%XFSEDR
        ZEXSED=ICEP%XEXSEDR
      CASE(5)
        ZFSED=ICEP%XFSEDS
        ZEXSED=ICEP%XEXSEDS
      CASE(6)
        ZFSED=ICEP%XFSEDG
        ZEXSED=ICEP%XEXSEDG
      CASE(7)
        ZFSED=ICEP%XFSEDH
        ZEXSED=ICEP%XEXSEDH
      CASE DEFAULT
        WRITE(YSPE, '( I10 )' ) KSPE
        CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'ICE4_SEDIMENTATION_SPLIT', 'no sedimentation parameter for KSPE='//TRIM(YSPE) )
    END SELECT
    !
    ZWSED(:,:,:) = 0.
    DO JL=1, ISEDIM
      JI=I1(JL)
      JJ=I2(JL)
      JK=I3(JL)
      IF(PRXT(JI,JJ,JK)>ICED%XRTMIN(KSPE)) THEN
        ZWSED(JI, JJ, JK) = ZFSED  * PRXT(JI, JJ, JK)**ZEXSED            &
                          &        * PRHODREF(JI, JJ, JK)**(ZEXSED-ICED%XCEXVT)
      ENDIF
    ENDDO
  ENDIF
  ZMAX_TSTEP(:,:) = ZREMAINT(:,:)
  DO JL=1, ISEDIM
    JI=I1(JL)
    JJ=I2(JL)
    JK=I3(JL)
    IF(PRXT(JI,JJ,JK)>ICED%XRTMIN(KSPE) .AND. ZWSED(JI, JJ, JK)>1.E-20) THEN
      ZMAX_TSTEP(JI, JJ) = MIN(ZMAX_TSTEP(JI, JJ), PMAXCFL * PRHODREF(JI, JJ, JK) * &
                         & PRXT(JI, JJ, JK) * PDZZ(JI, JJ, JK) / ZWSED(JI, JJ, JK))
    ENDIF
  ENDDO

  DO JJ = D%NJB, D%NJE
    DO JI = D%NIB, D%NIE
      ZREMAINT(JI,JJ) = ZREMAINT(JI,JJ) - ZMAX_TSTEP(JI,JJ)
      PINPRX(JI,JJ) = PINPRX(JI,JJ) + ZWSED(JI,JJ,D%NKB) / CST%XRHOLW * (ZMAX_TSTEP(JI,JJ) * ZINVTSTEP)
    ENDDO
  ENDDO

  DO JK = D%NKTB , D%NKTE
    DO JJ = D%NJB, D%NJE
      DO JI = D%NIB, D%NIE
        ZMRCHANGE = ZMAX_TSTEP(JI,JJ) * POORHODZ(JI,JJ,JK)*(ZWSED(JI,JJ,JK+D%NKL)-ZWSED(JI,JJ,JK))
        PRXT(JI,JJ,JK) = PRXT(JI,JJ,JK) + ZMRCHANGE + PPRXS(JI,JJ,JK) * ZMAX_TSTEP(JI,JJ)
        PRXS(JI,JJ,JK) = PRXS(JI,JJ,JK) + ZMRCHANGE * ZINVTSTEP
        IF (GPRESENT_PFPR) THEN
          PFPR(JI,JJ,JK,KSPE) = PFPR(JI,JJ,JK,KSPE) + ZWSED(JI,JJ,JK) * (ZMAX_TSTEP(JI,JJ) * ZINVTSTEP)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
!
END DO
!
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION_SPLIT:INTERNAL_SEDIM_SPLIT', 1, ZHOOK_HANDLE)
END SUBROUTINE INTERNAL_SEDIM_SPLI
!
END SUBROUTINE ICE4_SEDIMENTATION_SPLIT
END MODULE MODE_ICE4_SEDIMENTATION_SPLIT
