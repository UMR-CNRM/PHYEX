!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_ICE4_SEDIMENTATION_SPLIT
INTERFACE
SUBROUTINE ICE4_SEDIMENTATION_SPLIT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, ODEPOSC, PVDEPOSC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINDEP, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PSEA, PTOWN,  &
                                   &PINPRH, PRHT, PRHS, PFPR)
IMPLICIT NONE
INTEGER, INTENT(IN) :: KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT
INTEGER,                      INTENT(IN)              :: KKL     !vert. levels type 1=MNH -1=ARO
REAL,                         INTENT(IN)              :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                      INTENT(IN)              :: KRR     ! Number of moist variable
LOGICAL,                      INTENT(IN)              :: OSEDIC  ! Switch for droplet sedim.
LOGICAL,                      INTENT(IN)              :: ODEPOSC ! Switch for droplet depos.
REAL,                         INTENT(IN)              :: PVDEPOSC! Droplet deposition velocity
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
REAL, DIMENSION(:,:),         INTENT(OUT)             :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(:,:),         INTENT(OUT)             :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRR  ! Rain instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRI  ! Pristine ice instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRS  ! Snow instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(KIT,KJT), OPTIONAL,INTENT(IN)         :: PSEA    ! Sea Mask
REAL, DIMENSION(KIT,KJT), OPTIONAL,INTENT(IN)         :: PTOWN   ! Fraction that is town
REAL, DIMENSION(KIT,KJT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
END SUBROUTINE ICE4_SEDIMENTATION_SPLIT
END INTERFACE
END MODULE MODI_ICE4_SEDIMENTATION_SPLIT
SUBROUTINE ICE4_SEDIMENTATION_SPLIT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, ODEPOSC, PVDEPOSC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINDEP, PINPRR, PINPRI, PINPRS, PINPRG, &
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
!  P. Wautelet 11/02/2019: dimensions of PINPRC and PINDEP not necessarily KIT,KJT
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: XRHOLW
USE MODD_PARAM_ICE,      ONLY: XSPLIT_MAXCFL
USE MODD_RAIN_ICE_DESCR, ONLY: XALPHAC,XALPHAC2,XCONC_LAND,XCONC_SEA,XCONC_URBAN,XLBC,XNUC,XNUC2
USE MODD_RAIN_ICE_PARAM, ONLY: XFSEDC
!
USE MODE_MSG
!
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
LOGICAL,                      INTENT(IN)              :: ODEPOSC ! Switch for droplet depos.
REAL,                         INTENT(IN)              :: PVDEPOSC! Droplet deposition velocity
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
REAL, DIMENSION(:,:),         INTENT(OUT)             :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(:,:),         INTENT(OUT)             :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRR  ! Rain instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRI  ! Pristine ice instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRS  ! Snow instant precip
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(KIT,KJT), OPTIONAL,INTENT(IN)         :: PSEA    ! Sea Mask
REAL, DIMENSION(KIT,KJT), OPTIONAL,INTENT(IN)         :: PTOWN   ! Fraction that is town
REAL, DIMENSION(KIT,KJT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
!
INTEGER                                                             :: JI,JJ,JK
INTEGER                                                             :: IRR !Workaround of PGI bug with OpenACC (at least up to 18.10 version)
LOGICAL                                                             :: GDEPOSC, GSEDIC !Workaround of PGI bug with OpenACC (at least up to 18.10 version)
LOGICAL                                                             :: GPRESENT_PFPR, GPRESENT_PSEA
REAL                                                                :: ZINVTSTEP
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2))                  :: ZCONC_TMP    ! Weighted concentration
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),KKTB:KKTE)        :: ZW ! work array
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZCONC3D, & !  droplet condensation
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
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!
GDEPOSC = ODEPOSC
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
  ZLBC(:,:,:)   = XLBC(1)
  ZFSEDC(:,:,:) = XFSEDC(1)
  ZCONC3D(:,:,:)= XCONC_LAND
  ZCONC_TMP(:,:)= XCONC_LAND
  IF (GPRESENT_PSEA) THEN
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
        ZCONC3D(:,:,:) = XCONC_LAND
        ZRAY(:,:,:)  = 0.5*(GAMMA(XNUC+1.0/XALPHAC)/(GAMMA(XNUC)))
   END IF
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
IF (GSEDIC) THEN
  ZPRCS(:,:,:) = PRCS(:,:,:)-PRCT(:,:,:)*ZINVTSTEP
ENDIF
ZPRRS(:,:,:) = PRRS(:,:,:)-PRRT(:,:,:)*ZINVTSTEP
ZPRIS(:,:,:) = PRIS(:,:,:)-PRIT(:,:,:)*ZINVTSTEP
ZPRSS(:,:,:) = PRSS(:,:,:)-PRST(:,:,:)*ZINVTSTEP
ZPRGS(:,:,:) = PRGS(:,:,:)-PRGT(:,:,:)*ZINVTSTEP
IF ( IRR == 7 ) THEN
  ZPRHS(:,:,:) = PRHS(:,:,:)-PRHT(:,:,:)*ZINVTSTEP
END IF
!
! mr values inside the time-splitting loop
ZRCT(:,:,:) = PRCT(:,:,:)
ZRRT(:,:,:) = PRRT(:,:,:)
ZRIT(:,:,:) = PRIT(:,:,:)
ZRST(:,:,:) = PRST(:,:,:)
ZRGT(:,:,:) = PRGT(:,:,:)
IF (IRR==7) THEN
  ZRHT(:,:,:) = PRHT(:,:,:)
END IF
!
ZW(:,:,KKTB:KKTE) =1./(PRHODREF(:,:,KKTB:KKTE)* PDZZ(:,:,KKTB:KKTE))
!
!
!*       2.1   for cloud
!
IF (GSEDIC) THEN
    CALL INTERNAL_SEDIM_SPLI(KIB,KIE,KIT,KJB,KJE,KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &XSPLIT_MAXCFL, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &2, &
                          &ZRCT, PRCS, PINPRC, ZPRCS, &
                          &ZRAY, ZLBC, ZFSEDC, ZCONC3D, PFPR=PFPR)
ENDIF
!
!
!*       2.1bis  DROPLET DEPOSITION AT THE 1ST LEVEL ABOVE GROUND
!
IF (GDEPOSC) THEN
  PINDEP (:,:) = 0.
  DO JJ=KJB,KJE
    DO JI=KIB,KIE
      IF (PRCS(JI,JJ,KKB)>0.) THEN
        PRCS(JI,JJ,KKB) = PRCS(JI,JJ,KKB) - PVDEPOSC * PRCT(JI,JJ,KKB) / PDZZ(JI,JJ,KKB)
        PINPRC(JI,JJ) = PINPRC(JI,JJ) + PVDEPOSC * PRCT(JI,JJ,KKB) * PRHODREF(JI,JJ,KKB) /XRHOLW
        PINDEP(JI,JJ) = PVDEPOSC * PRCT(JI,JJ,KKB) * PRHODREF(JI,JJ,KKB) /XRHOLW
      END IF
    END DO
  END DO
END IF
!
!*       2.2   for rain
!
  CALL INTERNAL_SEDIM_SPLI(KIB,KIE,KIT,KJB,KJE,KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &XSPLIT_MAXCFL, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &3, &
                          &ZRRT, PRRS, PINPRR, ZPRRS, &
                          PFPR=PFPR)
!
!*       2.3   for pristine ice
!
  CALL INTERNAL_SEDIM_SPLI(KIB,KIE,KIT,KJB,KJE,KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &XSPLIT_MAXCFL, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &4, &
                          &ZRIT, PRIS, PINPRI, ZPRIS, &
                          PFPR=PFPR)
!
!*       2.4   for aggregates/snow
!
  CALL INTERNAL_SEDIM_SPLI(KIB,KIE,KIT,KJB,KJE,KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                        &XSPLIT_MAXCFL, &
                        &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                        &5, &
                        &ZRST, PRSS, PINPRS, ZPRSS, &
                        PFPR=PFPR)
!
!*       2.5   for graupeln
!
  CALL INTERNAL_SEDIM_SPLI(KIB,KIE,KIT,KJB,KJE,KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                        &XSPLIT_MAXCFL, &
                        &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                        &6, &
                        &ZRGT, PRGS, PINPRG, ZPRGS, &
                        PFPR=PFPR)
!
!*       2.6   for hail
!
IF (IRR==7) THEN
    CALL INTERNAL_SEDIM_SPLI(KIB,KIE,KIT,KJB,KJE,KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &XSPLIT_MAXCFL, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PTSTEP, &
                          &7, &
                          &ZRHT, PRHS, PINPRH, ZPRHS, &
                          PFPR=PFPR)
ENDIF
!
!
CONTAINS
!
!
!-------------------------------------------------------------------------------
!
!
SUBROUTINE INTERNAL_SEDIM_SPLI(KIB,KIE,KIT,KJB,KJE,KJT,KKB,KKTB,KKTE,KKT,KKL,KRR, &
                              &PMAXCFL,PRHODREF,POORHODZ,PDZZ,PPABST,PTHT,PTSTEP, &
                              &KSPE,PRXT,PRXS,PINPRX,PPRXS,                       &
                              &PRAY,PLBC,PFSEDC,PCONC3D,PFPR)
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: XCPD,XP00,XRD
USE MODD_RAIN_ICE_DESCR, ONLY: XCC,XCEXVT,XDC,XLBEXC,XRTMIN
USE MODD_RAIN_ICE_PARAM, ONLY: XEXCSEDI,XEXSEDG,XEXSEDH,XEXSEDR,XEXSEDS,XFSEDG,XFSEDH,XFSEDI,XFSEDR,XFSEDS
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KIB,KIE,KIT, KJB,KJE,KJT, KKB, KKTB, KKTE, KKT, KKL, KRR
REAL,                         INTENT(IN)              :: PMAXCFL ! maximum CFL allowed
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRHODREF ! Reference density
REAL, DIMENSION(KIT,KJT,KKTB:KKTE), INTENT(IN)        :: POORHODZ ! One Over (Rhodref times delta Z)
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PDZZ ! layer thikness (m)
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PPABST
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PTHT
REAL,                         INTENT(IN)              :: PTSTEP  ! total timestep
INTEGER,                      INTENT(IN)              :: KSPE ! 1 for rc, 2 for rr...
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRXT ! mr of specy X
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRXS !Tendency of the specy KSPE
REAL, DIMENSION(KIT,KJT),     INTENT(OUT)             :: PINPRX ! instant precip
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PPRXS ! external tendencie
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN), OPTIONAL    :: PRAY, PLBC, PFSEDC, PCONC3D
REAL, DIMENSION(KIT,KJT,KKT,KRR), INTENT(INOUT), OPTIONAL :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
character(len=10) :: yspe ! String for error message
INTEGER                         :: IDX, ISEDIM
INTEGER                         :: JI, JJ, JK, JL
INTEGER, DIMENSION(KIT*KJT*KKT) :: I1,I2,I3 ! Used to replace the COUNT
LOGICAL                         :: GPRESENT_PFPR
REAL                            :: ZINVTSTEP
REAL                            :: ZZWLBDC, ZRAY, ZZT, ZZWLBDA, ZZCC
REAL                            :: ZFSED, ZEXSED
REAL, DIMENSION(KIT, KJT)       :: ZMRCHANGE
REAL, DIMENSION(KIT, KJT)       :: ZMAX_TSTEP ! Maximum CFL in column
REAL, DIMENSION(SIZE(XRTMIN))   :: ZRSMIN
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)) :: ZREMAINT                   ! Remaining time until the timestep end
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),0:SIZE(PRHODREF,3)+1) :: ZWSED ! Sedimentation fluxes
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
ZRSMIN(:) = XRTMIN(:) * ZINVTSTEP
ZREMAINT(:,:) = PTSTEP
!
DO WHILE (ANY(ZREMAINT>0.))
  ISEDIM = 0
  DO JK = KKTB,KKTE
    DO JJ = KJB,KJE
      DO JI = KIB,KIE
        IF( (PRXT (JI,JJ,JK)>XRTMIN(KSPE) .OR.    &
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
      IF(PRXT(JI,JJ,JK)>XRTMIN(KSPE)) THEN
        ZZWLBDC = PLBC(JI,JJ,JK) * PCONC3D(JI,JJ,JK) / &
                  (PRHODREF(JI,JJ,JK) * PRXT(JI,JJ,JK))
        ZZWLBDC = ZZWLBDC**XLBEXC
        ZRAY = PRAY(JI,JJ,JK) / ZZWLBDC
        ZZT = PTHT(JI,JJ,JK) * (PPABST(JI,JJ,JK)/XP00)**(XRD/XCPD)
        ZZWLBDA = 6.6E-8*(101325./PPABST(JI,JJ,JK))*(ZZT/293.15)
        ZZCC = XCC*(1.+1.26*ZZWLBDA/ZRAY)
        ZWSED(JI, JJ, JK) = PRHODREF(JI,JJ,JK)**(-XCEXVT +1 ) *   &
                  ZZWLBDC**(-XDC)*ZZCC*PFSEDC(JI,JJ,JK) * PRXT(JI,JJ,JK)
      ENDIF
    ENDDO
  ELSEIF(KSPE==4) THEN
    ! ******* for pristine ice
    ZWSED(:,:,:) = 0.
    DO JL=1, ISEDIM
      JI=I1(JL)
      JJ=I2(JL)
      JK=I3(JL)
      IF(PRXT(JI, JJ, JK) .GT. MAX(XRTMIN(4), 1.0E-7)) THEN
        ZWSED(JI, JJ, JK) =  XFSEDI * PRXT(JI, JJ, JK) *  &
                            & PRHODREF(JI,JJ,JK)**(1.-XCEXVT) * & !    McF&H
                            & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                            &      ALOG(PRHODREF(JI,JJ,JK)*PRXT(JI,JJ,JK)) )**XEXCSEDI
      ENDIF
    ENDDO
  ELSE
    ! ******* for other species
    SELECT CASE(KSPE)
      CASE(3)
        ZFSED=XFSEDR
        ZEXSED=XEXSEDR
      CASE(5)
        ZFSED=XFSEDS
        ZEXSED=XEXSEDS
      CASE(6)
        ZFSED=XFSEDG
        ZEXSED=XEXSEDG
      CASE(7)
        ZFSED=XFSEDH
        ZEXSED=XEXSEDH
      CASE DEFAULT
        write( yspe, '( I10 )' ) kspe
        call Print_msg( NVERB_FATAL, 'GEN', 'ICE4_SEDIMENTATION_SPLIT', 'no sedimentation parameter for KSPE='//trim(yspe) )
    END SELECT
    !
    ZWSED(:,:,:) = 0.
    DO JL=1, ISEDIM
      JI=I1(JL)
      JJ=I2(JL)
      JK=I3(JL)
      IF(PRXT(JI,JJ,JK)>XRTMIN(KSPE)) THEN
        ZWSED(JI, JJ, JK) = ZFSED  * PRXT(JI, JJ, JK)**ZEXSED            &
                                   * PRHODREF(JI, JJ, JK)**(ZEXSED-XCEXVT)
      ENDIF
    ENDDO
  ENDIF
  ZMAX_TSTEP(:,:) = ZREMAINT(:,:)
  DO JL=1, ISEDIM
    JI=I1(JL)
    JJ=I2(JL)
    JK=I3(JL)
    IF(PRXT(JI,JJ,JK)>XRTMIN(KSPE) .AND. ZWSED(JI, JJ, JK)>1.E-20) THEN
      ZMAX_TSTEP(JI, JJ) = MIN(ZMAX_TSTEP(JI, JJ), PMAXCFL * PRHODREF(JI, JJ, JK) * &
                           PRXT(JI, JJ, JK) * PDZZ(JI, JJ, JK) / ZWSED(JI, JJ, JK))
    ENDIF
  ENDDO
  ZREMAINT(:,:) = ZREMAINT(:,:) - ZMAX_TSTEP(:,:)
  DO JK = KKTB , KKTE
    ZMRCHANGE(:,:) = ZMAX_TSTEP(:,:) * POORHODZ(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
    PRXT(:,:,JK) = PRXT(:,:,JK) + ZMRCHANGE(:,:) + PPRXS(:,:,JK) * ZMAX_TSTEP(:,:)
    PRXS(:,:,JK) = PRXS(:,:,JK) + ZMRCHANGE(:,:) * ZINVTSTEP
  ENDDO
  PINPRX(:,:) = PINPRX(:,:) + ZWSED(:,:,KKB) / XRHOLW * (ZMAX_TSTEP(:,:) * ZINVTSTEP)
  IF (GPRESENT_PFPR) THEN
    DO JK = KKTB , KKTE
      PFPR(:,:,JK,KSPE) = PFPR(:,:,JK,KSPE) + ZWSED(:,:,JK) * (ZMAX_TSTEP(:,:) * ZINVTSTEP)
    ENDDO
  ENDIF
!
END DO
!
END SUBROUTINE INTERNAL_SEDIM_SPLI
!
END SUBROUTINE ICE4_SEDIMENTATION_SPLIT
