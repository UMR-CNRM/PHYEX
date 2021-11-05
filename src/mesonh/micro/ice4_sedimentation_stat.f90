!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_ICE4_SEDIMENTATION_STAT
INTERFACE
SUBROUTINE ICE4_SEDIMENTATION_STAT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, ODEPOSC, PVDEPOSC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT,&
                                   &PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINDEP, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PSEA, PTOWN, &
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
REAL, DIMENSION(KIT,KJT),     OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(KIT,KJT),     OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(KIT,KJT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
END SUBROUTINE ICE4_SEDIMENTATION_STAT
END INTERFACE
END MODULE MODI_ICE4_SEDIMENTATION_STAT
SUBROUTINE ICE4_SEDIMENTATION_STAT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKTB, KKTE, KKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, ODEPOSC, PVDEPOSC, PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, &
                                  &PRSS, PRST, PRGS, PRGT,&
                                  &PINPRC, PINDEP, PINPRR, PINPRI, PINPRS, PINPRG, &
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
!  P. Wautelet 21/01/2021: initialize untouched part of PFPR
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST

USE MODE_MSG

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
REAL, DIMENSION(KIT,KJT),     OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(KIT,KJT),     OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(KIT,KJT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
!
INTEGER :: JK
!
REAL            :: ZINVTSTEP
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZW ! work array
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),0:SIZE(PRHODREF,3)+1)   &
                                  :: ZWSED        ! sedimentation fluxes
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)):: GDEP
!
!
!-------------------------------------------------------------------------------
!
ZINVTSTEP=1./PTSTEP

IF ( PRESENT( PFPR ) ) THEN
 !Set to 0. to avoid undefined values (in files)
 PFPR(:, :, : KKTB - 1, :) = 0.
 PFPR(:, :, KKTE + 1 :, :) = 0.
END IF
!-------------------------------------------------------------------------------
!
!*       1.    compute the fluxes
!
!
DO JK = KKTB , KKTE
  ZW(:,:,JK) =PTSTEP/(PRHODREF(:,:,JK)* PDZZ(:,:,JK) )
END DO
!
!*       2.1   for cloud
!
IF (OSEDIC) THEN
  CALL INTERNAL_SEDIM_STAT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKTB, KKTE, KKL, &
                          &PRHODREF, PDZZ, ZW, PPABST, PTHT, PTSTEP, &
                          &2, &
                          &PRCT, PRCS, ZWSED, PSEA, PTOWN)
  IF (PRESENT(PFPR)) THEN
    DO JK = KKTB , KKTE
      PFPR(:,:,JK,2)=ZWSED(:,:,JK)
    ENDDO
  ENDIF
  PINPRC(:,:) = ZWSED(:,:,KKB)/XRHOLW                        ! in m/s
ENDIF
!
!
!*       2.1bis  DROPLET DEPOSITION AT THE 1ST LEVEL ABOVE GROUND
!
IF (ODEPOSC) THEN
  GDEP(:,:) = .FALSE.
  GDEP(KIB:KIE,KJB:KJE) =    PRCS(KIB:KIE,KJB:KJE,KKB) >0 
  WHERE (GDEP)
     PRCS(:,:,KKB) = PRCS(:,:,KKB) - PVDEPOSC * PRCT(:,:,KKB) / PDZZ(:,:,KKB)
     PINPRC(:,:) = PINPRC(:,:) + PVDEPOSC * PRCT(:,:,KKB) * PRHODREF(:,:,KKB) /XRHOLW 
     PINDEP(:,:) = PVDEPOSC * PRCT(:,:,KKB) * PRHODREF(:,:,KKB) /XRHOLW 
  END WHERE
END IF
!
!*       2.2   for rain
!
CALL INTERNAL_SEDIM_STAT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKTB, KKTE, KKL, &
                        &PRHODREF, PDZZ, ZW, PPABST, PTHT, PTSTEP, &
                        &3, &
                        &PRRT, PRRS, ZWSED)
IF (PRESENT(PFPR)) THEN
  DO JK = KKTB , KKTE
    PFPR(:,:,JK,3)=ZWSED(:,:,JK)
  ENDDO
ENDIF
PINPRR(:,:) = ZWSED(:,:,KKB)/XRHOLW                        ! in m/s
!
!*       2.3   for pristine ice
!
CALL INTERNAL_SEDIM_STAT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKTB, KKTE, KKL, &
                        &PRHODREF, PDZZ, ZW, PPABST, PTHT, PTSTEP, &
                        &4, &
                        &PRIT, PRIS, ZWSED)
IF (PRESENT(PFPR)) THEN
  DO JK = KKTB , KKTE
    PFPR(:,:,JK,4)=ZWSED(:,:,JK)
  ENDDO
ENDIF
PINPRI(:,:) = ZWSED(:,:,KKB)/XRHOLW                        ! in m/s
!
!*       2.4   for aggregates/snow
!
CALL INTERNAL_SEDIM_STAT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKTB, KKTE, KKL, &
                        &PRHODREF, PDZZ, ZW, PPABST, PTHT, PTSTEP, &
                        &5, &
                        &PRST, PRSS, ZWSED)
IF (PRESENT(PFPR)) THEN
  DO JK = KKTB , KKTE
    PFPR(:,:,JK,5)=ZWSED(:,:,JK)
  ENDDO
ENDIF
PINPRS(:,:) = ZWSED(:,:,KKB)/XRHOLW                        ! in m/s
!
!*       2.5   for graupeln
!
CALL INTERNAL_SEDIM_STAT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKTB, KKTE, KKL, &
                        &PRHODREF, PDZZ, ZW, PPABST, PTHT, PTSTEP, &
                        &6, &
                        &PRGT, PRGS, ZWSED)
IF (PRESENT(PFPR)) THEN
  DO JK = KKTB , KKTE
    PFPR(:,:,JK,6)=ZWSED(:,:,JK)
  ENDDO
ENDIF
PINPRG(:,:) = ZWSED(:,:,KKB)/XRHOLW                        ! in m/s
!
!*       2.6   for hail
!
IF ( KRR == 7 ) THEN
  CALL INTERNAL_SEDIM_STAT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKTB, KKTE, KKL, &
                          &PRHODREF, PDZZ, ZW, PPABST, PTHT, PTSTEP, &
                          &7, &
                          &PRHT, PRHS, ZWSED)
  IF (PRESENT(PFPR)) THEN
    DO JK = KKTB , KKTE
      PFPR(:,:,JK,7)=ZWSED(:,:,JK)
    ENDDO
  ENDIF
  PINPRH(:,:) = ZWSED(:,:,KKB)/XRHOLW                        ! in m/s
ENDIF
!
!
CONTAINS
  SUBROUTINE INTERNAL_SEDIM_STAT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKTB, KKTE, KKL, &
                                &PRHODREF, PDZZ, PTSORHODZ, PPABST, PTHT, PTSTEP, &
                                &KSPE, &
                                &PRXT, PRXS, PWSED, PSEA, PTOWN)
    !
    !*      0. DECLARATIONS
    !          ------------
    !
    use mode_tools,  only: Countjv

    USE MODD_RAIN_ICE_DESCR
    USE MODD_RAIN_ICE_PARAM

    USE MODI_GAMMA

    IMPLICIT NONE
    !
    !*       0.1   Declarations of dummy arguments :
    !
    INTEGER, INTENT(IN) :: KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKT, KKE, KKTB, KKTE, KKL
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRHODREF ! Reference density
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PDZZ    ! Layer thikness (m)
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PTSORHODZ ! TimeStep Over (Rhodref times delta Z)
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PPABST
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PTHT
    REAL,                         INTENT(IN)              :: PTSTEP
    INTEGER,                      INTENT(IN)              :: KSPE ! 1 for rc, 2 for rr...
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)              :: PRXT ! mr of specy X
    REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)           :: PRXS !Tendency of the specy KSPE
    REAL, DIMENSION(KIT,KJT,0:KKT+1), INTENT(OUT)         :: PWSED ! sedimentation flux
    REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(IN)        :: PSEA    ! Sea Mask
    REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(IN)        :: PTOWN   ! Fraction that is town
    !
    !*       0.2  declaration of local variables
    !
    !
    character(len=10) :: yspe ! String for error message
    INTEGER :: JK, JCOUNT, JL, JI, JJ
    INTEGER, DIMENSION(SIZE(PRHODREF,1)*SIZE(PRHODREF,2)) :: I1, I2
    REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),0:SIZE(PRHODREF,3)+1)   &
                                  :: ZWSEDW1, ZWSEDW2    ! sedimentation speed
    REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)) :: ZQP
    REAL :: ZINVTSTEP, ZH, ZP1, ZP2, ZZWLBDA, ZZWLBDC, ZZCC
    REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZCONC3D !  droplet condensation
    REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) ::  &
                                     ZRAY,   & ! Cloud Mean radius
                                     ZLBC,   & ! XLBC weighted by sea fraction
                                     ZFSEDC
    REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2))                   &
                                  :: ZCONC_TMP    ! Weighted concentration
    REAL :: ZFSED, ZEXSED
    !
    !-------------------------------------------------------------------------------
    !
    !
    !*       1. Parameters for cloud sedimentation
    !
    IF(KSPE==2) THEN
     ZRAY(:,:,:)   = 0.
     ZLBC(:,:,:)   = XLBC(1)
     ZFSEDC(:,:,:) = XFSEDC(1)
     ZCONC3D(:,:,:)= XCONC_LAND
     ZCONC_TMP(:,:)= XCONC_LAND
     IF (PRESENT(PSEA)) THEN
      ZCONC_TMP(:,:)=PSEA(:,:)*XCONC_SEA+(1.-PSEA(:,:))*XCONC_LAND
      DO JK=KKTB,KKTE
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
    !
    ZINVTSTEP = 1./PTSTEP
    PWSED(:,:,:) = 0.
    ZWSEDW1(:,:,:) = 0.
    ZWSEDW2(:,:,:) = 0.
    ! calculation of ZP1, ZP2 and sedimentation flux
    DO JK = KKE , KKB, -1*KKL
      !estimation of q' taking into account incomming PWSED
      ZQP(:,:)=PWSED(:,:,JK+KKL)*PTSORHODZ(:,:,JK)
      JCOUNT=COUNTJV( (PRXT(:,:,JK) > XRTMIN(KSPE)) .OR. (ZQP(:,:) > XRTMIN(KSPE)) ,I1(:),I2(:))
      IF(KSPE==2) THEN
        !******* for cloud
        DO JL=1, JCOUNT
          JI=I1(JL)
          JJ=I2(JL)
          !calculation of w
          IF(PRXT(JI,JJ,JK) > XRTMIN(KSPE)) THEN
            ZZWLBDA=6.6E-8*(101325./PPABST(JI,JJ,JK))*(PTHT(JI,JJ,JK)/293.15)
            ZZWLBDC=(ZLBC(JI,JJ,JK)*ZCONC3D(JI,JJ,JK)  &
                  &/(PRHODREF(JI,JJ,JK)*PRXT(JI,JJ,JK)))**XLBEXC
            ZZCC=XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY(JI,JJ,JK)) !! ZCC  : Fall speed
            ZWSEDW1 (JI,JJ,JK)=PRHODREF(JI,JJ,JK)**(-XCEXVT ) *   &
                &  ZZWLBDC**(-XDC)*ZZCC*ZFSEDC(JI,JJ,JK)
          ENDIF
          IF ( ZQP(JI,JJ) > XRTMIN(KSPE) ) THEN
            ZZWLBDA=6.6E-8*(101325./PPABST(JI,JJ,JK))*(PTHT(JI,JJ,JK)/293.15)
            ZZWLBDC=(ZLBC(JI,JJ,JK)*ZCONC3D(JI,JJ,JK)  &
                 &/(PRHODREF(JI,JJ,JK)*ZQP(JI,JJ)))**XLBEXC
            ZZCC=XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY(JI,JJ,JK)) !! ZCC  : Fall speed
            ZWSEDW2 (JI,JJ,JK)=PRHODREF(JI,JJ,JK)**(-XCEXVT ) *   &
               &  ZZWLBDC**(-XDC)*ZZCC*ZFSEDC(JI,JJ,JK)
          ENDIF
        ENDDO
      ELSEIF(KSPE==4) THEN
        ! ******* for pristine ice
        DO JL=1, JCOUNT
          JI=I1(JL)
          JJ=I2(JL)
          !calculation of w
          IF ( PRXT(JI,JJ,JK) > MAX(XRTMIN(KSPE),1.0E-7 ) ) THEN
            ZWSEDW1 (JI,JJ,JK)= XFSEDI *  &
                              & PRHODREF(JI,JJ,JK)**(-XCEXVT) * & !    McF&H
                              & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                              &      ALOG(PRHODREF(JI,JJ,JK)*PRXT(JI,JJ,JK)) )**XEXCSEDI
          ENDIF
          IF ( ZQP(JI,JJ) > MAX(XRTMIN(KSPE),1.0E-7 ) ) THEN
            ZWSEDW2 (JI,JJ,JK)= XFSEDI *  &
                              & PRHODREF(JI,JJ,JK)**(-XCEXVT) * & !    McF&H
                              & MAX( 0.05E6,-0.15319E6-0.021454E6* &
                              &      ALOG(PRHODREF(JI,JJ,JK)*ZQP(JI,JJ)) )**XEXCSEDI
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
          write( yspe, '( I10 )' ) kspe
          call Print_msg( NVERB_FATAL, 'GEN', 'ICE4_SEDIMENTATION_STAT', &
                          'no sedimentation parameter for KSPE='//trim(yspe) )
        ENDIF
        DO JL=1, JCOUNT
          JI=I1(JL)
          JJ=I2(JL)
          !calculation of w
          IF ( PRXT(JI,JJ,JK) > XRTMIN(KSPE) ) THEN
            ZWSEDW1 (JI,JJ,JK)= ZFSED *PRXT(JI,JJ,JK)**(ZEXSED-1)* &
                                PRHODREF(JI,JJ,JK)**(ZEXSED-XCEXVT-1)
          ENDIF
          IF ( ZQP(JI,JJ) > XRTMIN(KSPE) ) THEN
            ZWSEDW2 (JI,JJ,JK)= ZFSED *ZQP(JI,JJ)**(ZEXSED-1)* &
                                PRHODREF(JI,JJ,JK)**(ZEXSED-XCEXVT-1)
          ENDIF
        ENDDO
      ENDIF
      DO JJ = KJB, KJE
        DO JI = KIB, KIE
          ZH=PDZZ(JI,JJ,JK)
          ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH )
          IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
            ZP2 = MAX(0.,1 -  ZH &
                & / (PTSTEP*ZWSEDW2(JI,JJ,JK)) )
          ELSE
            ZP2 = 0.
          ENDIF
          PWSED (JI,JJ,JK)=ZP1*PRHODREF(JI,JJ,JK)*&
                          &ZH*PRXT(JI,JJ,JK)&
                          &* ZINVTSTEP+ ZP2 * PWSED (JI,JJ,JK+KKL)
        ENDDO
      ENDDO
    ENDDO
    DO JK = KKTB , KKTE
      PRXS(:,:,JK) = PRXS(:,:,JK) + &
                   & PTSORHODZ(:,:,JK)*(PWSED(:,:,JK+KKL)-PWSED(:,:,JK))*ZINVTSTEP
    ENDDO
  END SUBROUTINE INTERNAL_SEDIM_STAT
  !
END SUBROUTINE ICE4_SEDIMENTATION_STAT
