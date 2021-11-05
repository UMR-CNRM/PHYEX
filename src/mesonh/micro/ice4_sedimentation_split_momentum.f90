!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_ICE4_SEDIMENTATION_SPLIT_MOMENTUM
INTERFACE
SUBROUTINE ICE4_SEDIMENTATION_SPLIT_MOMENTUM(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, OMOMENTUM, &
                                   &PSEA, PTOWN, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PINPRH, PRHT, PRHS, PFPR)
IMPLICIT NONE
INTEGER, INTENT(IN) :: KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT
INTEGER,                      INTENT(IN)              :: KKL     !vert. levels type 1=MNH -1=ARO
REAL,                         INTENT(IN)              :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                      INTENT(IN)              :: KRR     ! Number of moist variable
LOGICAL,                      INTENT(IN)              :: OSEDIC  ! Switch for droplet sedim.
LOGICAL,                      INTENT(IN)              :: OMOMENTUM  ! Switch to use momentum flux
REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PSEA    ! Sea Mask
REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PTOWN   ! Fraction that is town
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
REAL, DIMENSION(KIT,KJT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
END SUBROUTINE ICE4_SEDIMENTATION_SPLIT_MOMENTUM
END INTERFACE
END MODULE MODI_ICE4_SEDIMENTATION_SPLIT_MOMENTUM
SUBROUTINE ICE4_SEDIMENTATION_SPLIT_MOMENTUM(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, OMOMENTUM, &
                                   &PSEA, PTOWN, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PINPRH, PRHT, PRHS, PFPR)
!!
!!**  PURPOSE
!!    -------
!!      Computes the sedimentation
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the plitting of rain_ice source code (nov. 2014)
!!                and modified to use momentum
!!
!!    MODIFICATIONS
!!    -------------
!!
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_PARAM_ICE
USE MODD_RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM

USE MODE_MSG
use mode_tools,           only: Countjv

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
LOGICAL,                      INTENT(IN)              :: OMOMENTUM  ! Switch to use momentum flux
REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PSEA    ! Sea Mask
REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PTOWN   ! Fraction that is town
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
REAL, DIMENSION(KIT,KJT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
!
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
    :: GSEDIM ! Test where to compute the SED processes
INTEGER , DIMENSION(SIZE(GSEDIM)) :: I1,I2,I3 ! Used to replace the COUNT

REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZCONC3D, & !  droplet condensation
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
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZMOMC, ZMOMR, ZMOMI, ZMOMS, ZMOMG, ZMOMH ! momentum
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZMOMC_EXT, ZMOMR_EXT, ZMOMI_EXT, &
                                                                       ZMOMS_EXT, ZMOMG_EXT, ZMOMH_EXT ! momentum tendencies
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),0:SIZE(PRHODREF,3)+1) :: ZWSED        ! sedimentation fluxes
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)) :: ZCONC_TMP    ! Weighted concentration
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)) :: ZREMAINT ! Remaining time until the timestep end
REAL :: ZINVTSTEP
INTEGER :: ISEDIM ! ! Case number of sedimentation
INTEGER :: JK
LOGICAL :: FIRST
REAL, DIMENSION(SIZE(XRTMIN)) :: ZRSMIN
!-------------------------------------------------------------------------------
!
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
  FIRST = .TRUE.
  DO WHILE (ANY(ZREMAINT>0.))
    GSEDIM(:,:,:)=.FALSE.
    DO JK = KKTB , KKTE
      GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                    (ZRCT(KIB:KIE,KJB:KJE,JK)>XRTMIN(2) .OR.    &
                     ZPRCS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(2)) .AND. &
                    ZREMAINT(KIB:KIE,KJB:KJE)>0.
    ENDDO
    ISEDIM = COUNTJV(GSEDIM(:,:,:),I1(:),I2(:),I3(:))
    CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &OMOMENTUM, FIRST .AND. OMOMENTUM, &
                          &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PSEA, PTOWN, PTSTEP, &
                          &2, &
                          &ZRCT, PRCS, ZWSED, PINPRC, ZPRCS, ZMOMC, ZMOMC_EXT, &
                          &ZRAY, ZLBC, ZFSEDC, ZCONC3D, PFPR=PFPR)
    FIRST = .FALSE.
  ENDDO
ENDIF
!
!*       2.2   for rain
!
ZREMAINT(:,:) = PTSTEP
FIRST = .TRUE.
DO WHILE (ANY(ZREMAINT>0.))
  GSEDIM(:,:,:)=.FALSE.
  DO JK = KKTB , KKTE
    GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                  (ZRRT(KIB:KIE,KJB:KJE,JK)>XRTMIN(3) .OR.    &
                   ZPRRS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(3)) .AND. &
                  ZREMAINT(KIB:KIE,KJB:KJE)>0.
  ENDDO
  ISEDIM = COUNTJV(GSEDIM(:,:,:),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &OMOMENTUM, FIRST .AND. OMOMENTUM, &
                          &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PSEA, PTOWN, PTSTEP, &
                          &3, &
                          &ZRRT, PRRS, ZWSED, PINPRR, ZPRRS, ZMOMR, ZMOMR_EXT, &
                          &PFPR=PFPR)
  FIRST = .FALSE.
ENDDO
!
!*       2.3   for pristine ice
!
ZREMAINT(:,:) = PTSTEP
FIRST = .TRUE.
DO WHILE (ANY(ZREMAINT>0.))
  GSEDIM(:,:,:)=.FALSE.
  DO JK = KKTB , KKTE
    GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                  (ZRIT(KIB:KIE,KJB:KJE,JK)>XRTMIN(4) .OR.    &
                   ZPRIS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(4)) .AND. &
                  ZREMAINT(KIB:KIE,KJB:KJE)>0.
  ENDDO
  ISEDIM = COUNTJV(GSEDIM(:,:,:),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &OMOMENTUM, FIRST .AND. OMOMENTUM, &
                          &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PSEA, PTOWN, PTSTEP, &
                          &4, &
                          &ZRIT, PRIS, ZWSED, PINPRI, ZPRIS, ZMOMI, ZMOMI_EXT, PFPR=PFPR)
  FIRST = .FALSE.
ENDDO
!
!*       2.4   for aggregates/snow
!
ZREMAINT(:,:) = PTSTEP
FIRST = .TRUE.
DO WHILE (ANY(ZREMAINT>0.))
  GSEDIM(:,:,:)=.FALSE.
  DO JK = KKTB , KKTE
    GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                  (ZRST(KIB:KIE,KJB:KJE,JK)>XRTMIN(5) .OR.    &
                   ZPRSS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(5)) .AND. &
                  ZREMAINT(KIB:KIE,KJB:KJE)>0.
  ENDDO
  ISEDIM = COUNTJV(GSEDIM(:,:,:),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &OMOMENTUM, FIRST .AND. OMOMENTUM, &
                          &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PSEA, PTOWN, PTSTEP, &
                          &5, &
                          &ZRST, PRSS, ZWSED, PINPRS, ZPRSS, ZMOMS, ZMOMS_EXT, PFPR=PFPR)
  FIRST = .FALSE.
ENDDO
!
!*       2.5   for graupeln
!
ZREMAINT(:,:) = PTSTEP
FIRST = .TRUE.
DO WHILE (ANY(ZREMAINT>0.))
  GSEDIM(:,:,:)=.FALSE.
  DO JK = KKTB , KKTE
    GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                  (ZRGT(KIB:KIE,KJB:KJE,JK)>XRTMIN(6) .OR.    &
                   ZPRGS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(6)) .AND. &
                  ZREMAINT(KIB:KIE,KJB:KJE)>0.
  ENDDO
  ISEDIM = COUNTJV(GSEDIM(:,:,:),I1(:),I2(:),I3(:))
  CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                          &OMOMENTUM, FIRST .AND. OMOMENTUM, &
                          &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                          &PRHODREF, ZW, PDZZ, PPABST, PTHT, PSEA, PTOWN, PTSTEP, &
                          &6, &
                          &ZRGT, PRGS, ZWSED, PINPRG, ZPRGS, ZMOMG, ZMOMG_EXT, PFPR=PFPR)
  FIRST = .FALSE.
ENDDO
!
!*       2.6   for hail
!
IF (KRR==7) THEN
  ZREMAINT(:,:) = PTSTEP
  FIRST = .TRUE.
  DO WHILE (ANY(ZREMAINT>0.))
    GSEDIM(:,:,:)=.FALSE.
    DO JK = KKTB , KKTE
      GSEDIM(KIB:KIE,KJB:KJE,JK) =                              &
                    (ZRHT(KIB:KIE,KJB:KJE,JK)>XRTMIN(7) .OR.    &
                     ZPRHS(KIB:KIE,KJB:KJE,JK)>ZRSMIN(7)) .AND. &
                    ZREMAINT(KIB:KIE,KJB:KJE)>0.
    ENDDO
    ISEDIM = COUNTJV(GSEDIM(:,:,:),I1(:),I2(:),I3(:))
    CALL INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                            &OMOMENTUM, FIRST .AND. OMOMENTUM, &
                            &ISEDIM, GSEDIM, I1, I2, I3, XSPLIT_MAXCFL, ZREMAINT, &
                            &PRHODREF, ZW, PDZZ, PPABST, PTHT, PSEA, PTOWN, PTSTEP, &
                            &7, &
                            &ZRHT, PRHS, ZWSED, PINPRH, ZPRHS, ZMOMH, ZMOMH_EXT, PFPR=PFPR)
    FIRST = .FALSE.
  END DO
ENDIF
!
!
CONTAINS
!
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE INTERNAL_SEDIM_SPLI(KIT, KJT, KKB, KKTB, KKTE, KKT, KKL, KRR, &
                                &OMOMENTUM, OCOMPUTE_MOM, &
                                &KSEDIM, LDSEDIM, I1, I2, I3, PMAXCFL, PREMAINT, &
                                &PRHODREF, POORHODZ, PDZZ, PPABST, PTHT, PSEA, PTOWN, PTSTEP, &
                                &KSPE, &
                                &PRXT, PRXS, PWSED, PINPRX, PPRXS, PMOM, PMOM_EXT, &
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
    LOGICAL, INTENT(IN) :: OMOMENTUM, OCOMPUTE_MOM
    INTEGER, INTENT(IN) :: KSEDIM
    LOGICAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: LDSEDIM
    INTEGER, DIMENSION(KSEDIM), INTENT(IN) :: I1, I2, I3
    REAL,                         INTENT(IN)              :: PMAXCFL ! maximum CFL allowed
    REAL, DIMENSION(KIT,KJT),     INTENT(INOUT)           :: PREMAINT ! Time remaining until the end of the timestep
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRHODREF ! Reference density
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: POORHODZ ! One Over (Rhodref times delta Z)
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PDZZ ! layer thikness (m)
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PPABST
    REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PSEA    ! Sea Mask
    REAL, DIMENSION(KIT,KJT),     INTENT(IN)              :: PTOWN   ! Fraction that is town
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PTHT
    REAL,                         INTENT(IN)              :: PTSTEP  ! total timestep
    INTEGER,                      INTENT(IN)              :: KSPE ! 1 for rc, 2 for rr...
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRXT ! mr of specy X
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRXS !Tendency of the specy KSPE
    REAL, DIMENSION(KIT,KJT,0:KKT+1), INTENT(OUT)         :: PWSED ! sedimentation flux
    REAL, DIMENSION(KIT,KJT),     INTENT(INOUT)           :: PINPRX ! instant precip
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PPRXS ! external tendencie
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PMOM ! momentum associated to PRXT
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PMOM_EXT ! momentum tendency associated to PPRXS
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN), OPTIONAL    :: PRAY, PLBC, PFSEDC, PCONC3D
    REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(INOUT)   :: PFPR    ! upper-air precipitation fluxes
    !
    !*       0.2  declaration of local variables
    !
    !
    character(len=10) :: yspe ! String for error message
    INTEGER :: JK, JL, JI, JJ
    REAL :: ZINVTSTEP
    REAL :: ZZWLBDC, ZRAY, ZZT, ZZWLBDA, ZZCC
    REAL :: ZFSED, ZEXSED
    REAL, DIMENSION(KIT, KJT) :: ZMRCHANGE
    REAL, DIMENSION(KIT, KJT) :: ZMAX_TSTEP ! Maximum CFL in column
    REAL, DIMENSION(KIT,KJT,0:KKT+1) :: ZWSED_MOM ! Momentum flux
    REAL, DIMENSION(SIZE(XRTMIN)) :: ZRSMIN
    !
    !-------------------------------------------------------------------------------
    !
    !
    !*       1. Parameters for cloud sedimentation
    !
    !
    IF(OCOMPUTE_MOM .AND. .NOT. OMOMENTUM) THEN
      call Print_msg( NVERB_FATAL, 'GEN', 'ICE4_SEDIMENTATION_SPLIT_MOMENTUM',  &
                      'OCOMPUTE_MOM cannot be .TRUE. if we do not use momentum' )
    ENDIF
    !*       2.    compute the fluxes
    !
    !
    ZINVTSTEP = 1./PTSTEP
    ZRSMIN(:) = XRTMIN(:) * ZINVTSTEP
    IF(KSPE==2) THEN
      !******* for cloud
      IF(OCOMPUTE_MOM .OR. .NOT. OMOMENTUM) THEN
        PWSED(:,:,:) = 0.
        PMOM_EXT(:,:,:) = 0.
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
          IF(PPRXS(JI,JJ,JK)>ZRSMIN(KSPE) .AND. OCOMPUTE_MOM) THEN
            ZZWLBDC = PLBC(JI,JJ,JK) * PCONC3D(JI,JJ,JK) / &
                     (PRHODREF(JI,JJ,JK) * PPRXS(JI,JJ,JK) * PTSTEP)
            ZZWLBDC = ZZWLBDC**XLBEXC
            ZRAY = PRAY(JI,JJ,JK) / ZZWLBDC
            ZZT = PTHT(JI,JJ,JK) * (PPABST(JI,JJ,JK)/XP00)**(XRD/XCPD)
            ZZWLBDA = 6.6E-8*(101325./PPABST(JI,JJ,JK))*(ZZT/293.15)
            ZZCC = XCC*(1.+1.26*ZZWLBDA/ZRAY)
            PMOM_EXT(JI, JJ, JK) = PRHODREF(JI,JJ,JK)**(-XCEXVT +1 -1) *   &
                     ZZWLBDC**(-XDC)*ZZCC*PFSEDC(JI,JJ,JK) * PPRXS(JI,JJ,JK)
          ENDIF
        ENDDO
        IF(OCOMPUTE_MOM) PMOM(:, :, :)=PWSED(:, :, 1:KKT)
      ENDIF
    ELSEIF(KSPE==4) THEN
      ! ******* for pristine ice
      IF(OCOMPUTE_MOM .OR. .NOT. OMOMENTUM) THEN
        PWSED(:,:,:) = 0.
        PMOM_EXT(:,:,:) = 0.
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
          IF(PPRXS(JI,JJ,JK)>MAX(ZRSMIN(4), 1.0E-7/PTSTEP) .AND. OCOMPUTE_MOM) THEN
            PMOM_EXT(JI, JJ, JK) = XFSEDI * PPRXS(JI, JJ, JK) *  &
                                & PRHODREF(JI,JJ,JK)**(1.-XCEXVT-1) * & !    McF&H
                                & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                                &      ALOG(PRHODREF(JI,JJ,JK)*PPRXS(JI,JJ,JK)*PTSTEP) )**XEXCSEDI
          ENDIF
        ENDDO
        IF(OCOMPUTE_MOM) PMOM(:, :, :)=PWSED(:, :, 1:KKT)
      ENDIF
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
        write( yspe, '( I10 )' ) kspe
        call Print_msg( NVERB_FATAL, 'GEN', 'ICE4_SEDIMENTATION_SPLIT_MOMENTUM', &
                        'no sedimentation parameter for KSPE='//trim(yspe) )
      ENDIF
      IF(OCOMPUTE_MOM .OR. .NOT. OMOMENTUM) THEN
        !Momentum (per m3) and mass flux are given by the same formulae
        PWSED(:,:,:) = 0.
        PMOM_EXT(:,:,:) = 0.
        DO JL=1, KSEDIM
          JI=I1(JL)
          JJ=I2(JL)
          JK=I3(JL)
          IF(PRXT(JI,JJ,JK)>XRTMIN(KSPE)) THEN
            PWSED(JI, JJ, JK) = ZFSED  * PRXT(JI, JJ, JK)**ZEXSED *   &
                                         PRHODREF(JI, JJ, JK)**(ZEXSED-XCEXVT)
          ENDIF
          IF(PPRXS(JI,JJ,JK)>ZRSMIN(KSPE) .AND. OCOMPUTE_MOM) THEN
            PMOM_EXT(JI, JJ, JK) = ZFSED  * (PPRXS(JI, JJ, JK)*PTSTEP)**ZEXSED *   &
                                         PRHODREF(JI, JJ, JK)**(ZEXSED-XCEXVT-1) * ZINVTSTEP
          ENDIF
        ENDDO
        IF(OCOMPUTE_MOM) PMOM(:, :, :)=PWSED(:, :, 1:KKT) / PRHODREF(:, :, :) ! momentum per kg of dry air
      ENDIF
    ENDIF
    IF(OMOMENTUM) THEN
      PWSED(:,:,:) = 0.
      ZWSED_MOM(:,:,:) = 0.
      DO JL=1, KSEDIM
        JI=I1(JL)
        JJ=I2(JL)
        JK=I3(JL)
        IF(PRXT(JI,JJ,JK)>XRTMIN(KSPE)) THEN
          ZWSED_MOM(JI, JJ, JK) = PMOM(JI, JJ, JK)**2 / PRXT(JI, JJ, JK) * PRHODREF(JI, JJ, JK) ! (kg*m/s)/(s*m**2)
        ENDIF
      ENDDO
      PWSED(:, :, 1:KKT) = PMOM(:, :, :)*PRHODREF(:, :, :) !PMOM divided by r to get speed and multiply by rho*r to get flux
    ENDIF
    ZMAX_TSTEP(:,:) = PREMAINT(:,:)
    DO JL=1, KSEDIM
      JI=I1(JL)
      JJ=I2(JL)
      JK=I3(JL)
      IF(PRXT(JI,JJ,JK)>XRTMIN(KSPE)) THEN
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
    IF(OMOMENTUM) THEN
      DO JK = KKTB , KKTE
        PMOM(:,:,JK) = PMOM(:,:,JK) + ZMAX_TSTEP(:,:) * POORHODZ(:,:,JK) * (ZWSED_MOM(:,:,JK+KKL)-ZWSED_MOM(:,:,JK))
        PMOM(:,:,JK) = PMOM(:,:,JK) + ZMAX_TSTEP(:,:) * PMOM_EXT(:,:,JK)
        PMOM(:,:,JK) = MAX(0., PMOM(:,:,JK))
      ENDDO
    ENDIF
    PINPRX(:,:) = PINPRX(:,:) + ZWSED(:,:,KKB) / XRHOLW * (ZMAX_TSTEP(:,:) * ZINVTSTEP)
    IF (PRESENT(PFPR)) THEN
      DO JK = KKTB , KKTE
        PFPR(:,:,JK,KSPE) = PFPR(:,:,JK,KSPE) + ZWSED(:,:,JK) * (ZMAX_TSTEP(:,:) * ZINVTSTEP)
      ENDDO
    ENDIF
    !
  END SUBROUTINE INTERNAL_SEDIM_SPLI
  !
END SUBROUTINE ICE4_SEDIMENTATION_SPLIT_MOMENTUM
