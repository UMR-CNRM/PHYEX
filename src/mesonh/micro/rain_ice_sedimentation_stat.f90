!MNH_LIC Copyright 1995-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 25/02/2019: split rain_ice (cleaner and easier to maintain/debug)
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_SEDIMENTATION_STAT

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: RAIN_ICE_SEDIMENTATION_STAT

CONTAINS

SUBROUTINE RAIN_ICE_SEDIMENTATION_STAT( KIB, KIE, KJB, KJE, KKB, KKE, KKTB, KKTE, KKT, KKL, KRR,                &
                                        PTSTEP, OSEDIC, PINPRC, PINDEP,                                         &
                                        PINPRR, PINPRS, PINPRG, PDZZ, PRHODREF, PPABST, PTHT, PRHODJ, PINPRR3D, &
                                        PRCS, PRCT, PRRS, PRRT, PRIS, PRSS, PRST, PRGS, PRGT,                   &
                                        PSEA, PTOWN, PINPRH, PRHS, PRHT, PFPR )
!
!*      0. DECLARATIONS
!          ------------
!
use modd_budget,         only: lbudget_rc, lbudget_rr, lbudget_ri, lbudget_rs, lbudget_rg, lbudget_rh, &
                               NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, &
                               tbudgets
use MODD_CST,            only: XRHOLW
use MODD_PARAM_ICE,      only: LDEPOSC, XVDEPOSC
use MODD_RAIN_ICE_PARAM, only: XEXSEDG, XEXSEDH, XEXCSEDI, XEXSEDR, XEXSEDS, &
                               XFSEDC, XFSEDG, XFSEDH, XFSEDI, XFSEDR, XFSEDS
use MODD_RAIN_ICE_DESCR, only: XALPHAC, XALPHAC2, XCC, XCEXVT, XCONC_LAND, XCONC_SEA, XCONC_URBAN, &
                               XDC, XLBC, XLBEXC, XNUC, XNUC2, XRTMIN

use mode_budget,         only: Budget_store_init, Budget_store_end
use mode_tools,          only: Countjv

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                            INTENT(IN)    :: KIB, KIE, KJB, KJE, KKB, KKE, KKTB, KKTE, KKT
INTEGER,                            INTENT(IN)    :: KKL     ! vert. levels type 1=MNH -1=ARO
INTEGER,                            INTENT(IN)    :: KRR     ! Number of moist variable
REAL,                               INTENT(IN)    :: PTSTEP  ! Double Time step
                                                             ! (single if cold start)
LOGICAL,                            INTENT(IN)    :: OSEDIC  ! Switch for droplet sedim.
REAL, DIMENSION(:,:),               INTENT(INOUT) :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(:,:),               INTENT(INOUT) :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(:,:),               INTENT(INOUT) :: PINPRR  ! Rain instant precip
REAL, DIMENSION(:,:),               INTENT(INOUT) :: PINPRS  ! Snow instant precip
REAL, DIMENSION(:,:),               INTENT(INOUT) :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),             INTENT(OUT)   :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(:,:,:),             INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:),             INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),             INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),             INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),             INTENT(INOUT) :: PRGS    ! Graupel m.r. source
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(INOUT) :: PINPRH  ! Hail instant precip
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
INTEGER                           :: JI,JJ,JK
INTEGER                           :: JCOUNT, JL
INTEGER, DIMENSION(SIZE(PRHODREF,1)*SIZE(PRHODREF,2)) :: I1, I2
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)) :: GDEP
REAL                              :: ZINVTSTEP
REAL                              :: ZP1,ZP2,ZH,ZZWLBDA,ZZWLBDC,ZZCC
REAL,    DIMENSION(SIZE(XRTMIN))     :: ZRTMIN
! XRTMIN = Minimum value for the mixing ratio
! ZRTMIN = Minimum value for the source (tendency)
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)) :: ZQP
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2))                   &
                                  :: ZCONC_TMP    ! Weighted concentration
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZPRCS,ZPRRS,ZPRSS,ZPRGS,ZPRHS   ! Mixing ratios created during the time step
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZCONC3D !  droplet condensation
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) ::  &
                                     ZRAY,   & ! Cloud Mean radius
                                     ZLBC,   & ! XLBC weighted by sea fraction
                                     ZFSEDC
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZW ! work array
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),0:SIZE(PRHODREF,3)+1)   &
                                  :: ZWSED        ! sedimentation fluxes
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),0:SIZE(PRHODREF,3)+1)   &
                                  :: ZWSEDW1       ! sedimentation speed
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),0:SIZE(PRHODREF,3)+1)   &
                                  :: ZWSEDW2       ! sedimentation speed
!-------------------------------------------------------------------------------

if ( lbudget_rc .and. osedic ) call Budget_store_init( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rr )              call Budget_store_init( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_ri )              call Budget_store_init( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rs )              call Budget_store_init( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rg )              call Budget_store_init( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rh )              call Budget_store_init( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )

ZINVTSTEP=1./PTSTEP
!
!*       1. Parameters for cloud sedimentation
!
  IF (OSEDIC) THEN
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
  IF (LDEPOSC) PINDEP (:,:) = 0.
!
!*       2.    compute the fluxes
!


ZRTMIN(:)    = XRTMIN(:) * ZINVTSTEP
!
IF (OSEDIC) THEN
  ZPRCS(:,:,:) = 0.0
  ZPRCS(:,:,:) = PRCS(:,:,:)-PRCT(:,:,:)* ZINVTSTEP
  PRCS(:,:,:)  = PRCT(:,:,:)* ZINVTSTEP
END IF
ZPRRS(:,:,:) = 0.0
ZPRSS(:,:,:) = 0.0
ZPRGS(:,:,:) = 0.0
IF ( KRR == 7 ) ZPRHS(:,:,:) = 0.0
!
ZPRRS(:,:,:) = PRRS(:,:,:)-PRRT(:,:,:)* ZINVTSTEP
ZPRSS(:,:,:) = PRSS(:,:,:)-PRST(:,:,:)* ZINVTSTEP
ZPRGS(:,:,:) = PRGS(:,:,:)-PRGT(:,:,:)* ZINVTSTEP
IF ( KRR == 7 ) ZPRHS(:,:,:) = PRHS(:,:,:)-PRHT(:,:,:)* ZINVTSTEP
PRRS(:,:,:)  = PRRT(:,:,:)* ZINVTSTEP
PRSS(:,:,:)  = PRST(:,:,:)* ZINVTSTEP
PRGS(:,:,:)  = PRGT(:,:,:)* ZINVTSTEP
IF ( KRR == 7 ) PRHS(:,:,:)  = PRHT(:,:,:)* ZINVTSTEP
!
IF (OSEDIC) PRCS(:,:,:) = PRCS(:,:,:) + ZPRCS(:,:,:)
PRRS(:,:,:) = PRRS(:,:,:) + ZPRRS(:,:,:)
PRSS(:,:,:) = PRSS(:,:,:) + ZPRSS(:,:,:)
PRGS(:,:,:) = PRGS(:,:,:) + ZPRGS(:,:,:)
IF ( KRR == 7 ) PRHS(:,:,:) = PRHS(:,:,:) + ZPRHS(:,:,:)
IF (PRESENT(PFPR)) PFPR(:,:,:,:) = 0.
DO JK = KKTB , KKTE
  ZW(:,:,JK) =PTSTEP/(PRHODREF(:,:,JK)* PDZZ(:,:,JK) )
END DO
PINPRR3D (:,:,:) = 0.

!
!*       2.1   for cloud
!
 IF (OSEDIC) THEN
     PRCS(:,:,:) = PRCS(:,:,:) * PTSTEP
     ZWSED(:,:,:) = 0.
     ZWSEDW1(:,:,:) = 0.
     ZWSEDW2(:,:,:) = 0.

! calculation of P1, P2 and sedimentation flux
     DO JK = KKE , KKB, -1*KKL
       !estimation of q' taking into account incomming ZWSED
       ZQP(:,:)=ZWSED(:,:,JK+KKL)*ZW(:,:,JK)

       JCOUNT=COUNTJV((PRCS(:,:,JK) > ZRTMIN(2) .AND. PRCT(:,:,JK) > ZRTMIN(2)) .OR. &
                      (ZQP(:,:) > ZRTMIN(2)),I1(:),I2(:))
       DO JL=1, JCOUNT
         JI=I1(JL)
         JJ=I2(JL)
         !calculation of w
         ! mars 2009 : ajout d'un test
         !IF ( PRCS(JI,JJ,JK) > ZRTMIN(2) ) THEN
         IF(PRCS(JI,JJ,JK) > ZRTMIN(2) .AND. PRCT(JI,JJ,JK) > ZRTMIN(2)) THEN
           ZZWLBDA=6.6E-8*(101325./PPABST(JI,JJ,JK))*(PTHT(JI,JJ,JK)/293.15)
           ZZWLBDC=(ZLBC(JI,JJ,JK)*ZCONC3D(JI,JJ,JK)  &
                &/(PRHODREF(JI,JJ,JK)*PRCT(JI,JJ,JK)))**XLBEXC
           ZZCC=XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY(JI,JJ,JK)) !! ZCC  : Fall speed
           ZWSEDW1 (JI,JJ,JK)=PRHODREF(JI,JJ,JK)**(-XCEXVT ) *   &
             &  ZZWLBDC**(-XDC)*ZZCC*ZFSEDC(JI,JJ,JK)
         ENDIF
         IF ( ZQP(JI,JJ) > ZRTMIN(2) ) THEN
           ZZWLBDA=6.6E-8*(101325./PPABST(JI,JJ,JK))*(PTHT(JI,JJ,JK)/293.15)
           ZZWLBDC=(ZLBC(JI,JJ,JK)*ZCONC3D(JI,JJ,JK)  &
                &/(PRHODREF(JI,JJ,JK)*ZQP(JI,JJ)))**XLBEXC
           ZZCC=XCC*(1.+1.26*ZZWLBDA*ZZWLBDC/ZRAY(JI,JJ,JK)) !! ZCC  : Fall speed
           ZWSEDW2 (JI,JJ,JK)=PRHODREF(JI,JJ,JK)**(-XCEXVT ) *   &
             &  ZZWLBDC**(-XDC)*ZZCC*ZFSEDC(JI,JJ,JK)
         ENDIF
       ENDDO

       DO JJ = KJB, KJE
         DO JI = KIB, KIE
           ZH=PDZZ(JI,JJ,JK)
           ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH)
           ! mars 2009 : correction : ZWSEDW1 =>  ZWSEDW2
           !IF (ZWSEDW1(JI,JJ,JK) /= 0.) THEN
           IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
             ZP2 = MAX(0.,1 -  ZH &
           &  / (PTSTEP*ZWSEDW2(JI,JJ,JK)) )
           ELSE
             ZP2 = 0.
           ENDIF
           ZWSED (JI,JJ,JK)=ZP1*PRHODREF(JI,JJ,JK)*&
           &ZH*PRCS(JI,JJ,JK)&
           &* ZINVTSTEP+ ZP2 * ZWSED (JI,JJ,JK+KKL)
         ENDDO
       ENDDO
     ENDDO

     DO JK = KKTB , KKTE
       PRCS(:,:,JK) = PRCS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
     END DO
     IF (PRESENT(PFPR)) THEN
       DO JK = KKTB , KKTE
         PFPR(:,:,JK,2)=ZWSED(:,:,JK)
       ENDDO
     ENDIF

     PINPRC(:,:) = ZWSED(:,:,KKB)/XRHOLW                        ! in m/s
     PRCS(:,:,:) = PRCS(:,:,:) * ZINVTSTEP
 ENDIF

!
!*       2.2   for rain
!

   PRRS(:,:,:) = PRRS(:,:,:) * PTSTEP
   ZWSED(:,:,:) = 0.
   ZWSEDW1(:,:,:) = 0.
   ZWSEDW2(:,:,:) = 0.

! calculation of ZP1, ZP2 and sedimentation flux
   DO JK = KKE , KKB, -1*KKL
     !estimation of q' taking into account incomming ZWSED
     ZQP(:,:)=ZWSED(:,:,JK+KKL)*ZW(:,:,JK)

     JCOUNT=COUNTJV((PRRS(:,:,JK) > ZRTMIN(3)) .OR. &
                    (ZQP(:,:) > ZRTMIN(3)),I1(:),I2(:))
     DO JL=1, JCOUNT
       JI=I1(JL)
       JJ=I2(JL)
       !calculation of w
       IF ( PRRS(JI,JJ,JK) > ZRTMIN(3) ) THEN
         ZWSEDW1 (JI,JJ,JK)= XFSEDR *PRRS(JI,JJ,JK)**(XEXSEDR-1)* &
         PRHODREF(JI,JJ,JK)**(XEXSEDR-XCEXVT-1)
       ENDIF
       IF ( ZQP(JI,JJ) > ZRTMIN(3) ) THEN
         ZWSEDW2 (JI,JJ,JK)= XFSEDR *(ZQP(JI,JJ))**(XEXSEDR-1)* &
         PRHODREF(JI,JJ,JK)**(XEXSEDR-XCEXVT-1)
       ENDIF
     ENDDO
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
         ZWSED (JI,JJ,JK)=ZP1*PRHODREF(JI,JJ,JK)*&
         &ZH*PRRS(JI,JJ,JK)&
         &* ZINVTSTEP+ ZP2 * ZWSED (JI,JJ,JK+KKL)
       ENDDO
     ENDDO
   ENDDO

   DO JK = KKTB , KKTE
     PRRS(:,:,JK) = PRRS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
   ENDDO
   IF (PRESENT(PFPR)) THEN
     DO JK = KKTB , KKTE
       PFPR(:,:,JK,3)=ZWSED(:,:,JK)
     ENDDO
   ENDIF
   PINPRR(:,:) = ZWSED(:,:,KKB)/XRHOLW                        ! in m/s
   PINPRR3D(:,:,:) = ZWSED(:,:,1:KKT)/XRHOLW                        ! in m/s
   PRRS(:,:,:) = PRRS(:,:,:) * ZINVTSTEP

!
!*       2.3   for pristine ice
!

   PRIS(:,:,:) = PRIS(:,:,:) * PTSTEP
   ZWSED(:,:,:) = 0.
   ZWSEDW1(:,:,:) = 0.
   ZWSEDW2(:,:,:) = 0.
! calculation of ZP1, ZP2 and sedimentation flux
   DO JK = KKE , KKB, -1*KKL
     !estimation of q' taking into account incomming ZWSED
     ZQP(:,:)=ZWSED(:,:,JK+KKL)*ZW(:,:,JK)

     JCOUNT=COUNTJV((PRIS(:,:,JK) > MAX(ZRTMIN(4),1.0E-7 )) .OR. &
                    (ZQP(:,:) > MAX(ZRTMIN(4),1.0E-7 )),I1(:),I2(:))
     DO JL=1, JCOUNT
       JI=I1(JL)
       JJ=I2(JL)
       !calculation of w
       IF ( PRIS(JI,JJ,JK) > MAX(ZRTMIN(4),1.0E-7 ) ) THEN
         ZWSEDW1 (JI,JJ,JK)= XFSEDI *  &
         &  PRHODREF(JI,JJ,JK)**(XCEXVT) * & !    McF&H
         &  MAX( 0.05E6,-0.15319E6-0.021454E6* &
         &  ALOG(PRHODREF(JI,JJ,JK)*PRIS(JI,JJ,JK)) )**XEXCSEDI
       ENDIF
       IF ( ZQP(JI,JJ) > MAX(ZRTMIN(4),1.0E-7 ) ) THEN
         ZWSEDW2 (JI,JJ,JK)= XFSEDI *  &
         &  PRHODREF(JI,JJ,JK)**(XCEXVT) * & !    McF&H
         &  MAX( 0.05E6,-0.15319E6-0.021454E6* &
         &  ALOG(PRHODREF(JI,JJ,JK)*ZQP(JI,JJ)) )**XEXCSEDI
       ENDIF
     ENDDO
     DO JJ = KJB, KJE
       DO JI = KIB, KIE
         ZH=PDZZ(JI,JJ,JK)
         ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH )
         IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
           ZP2 = MAX(0.,1 - ZH  &
           &  / (PTSTEP*ZWSEDW2(JI,JJ,JK)) )
         ELSE
           ZP2 = 0.
         ENDIF
         ZWSED (JI,JJ,JK)=ZP1*PRHODREF(JI,JJ,JK)*&
         &ZH*PRIS(JI,JJ,JK)&
         &* ZINVTSTEP+ ZP2 * ZWSED (JI,JJ,JK+KKL)
       ENDDO
     ENDDO
   ENDDO

   DO JK = KKTB , KKTE
     PRIS(:,:,JK) = PRIS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
   ENDDO
   IF (PRESENT(PFPR)) THEN
     DO JK = KKTB , KKTE
       PFPR(:,:,JK,4)=ZWSED(:,:,JK)
     ENDDO
   ENDIF

   PRIS(:,:,:) = PRIS(:,:,:) * ZINVTSTEP


!
!*       2.4   for aggregates/snow
!

   PRSS(:,:,:) = PRSS(:,:,:) * PTSTEP
   ZWSED(:,:,:) = 0.
   ZWSEDW1(:,:,:) = 0.
   ZWSEDW2(:,:,:) = 0.

! calculation of ZP1, ZP2 and sedimentation flux
   DO JK = KKE , KKB, -1*KKL
     !estimation of q' taking into account incomming ZWSED
     ZQP(:,:)=ZWSED(:,:,JK+KKL)*ZW(:,:,JK)

     JCOUNT=COUNTJV((PRSS(:,:,JK) > ZRTMIN(5)) .OR. &
                    (ZQP(:,:) > ZRTMIN(5)),I1(:),I2(:))
     DO JL=1, JCOUNT
       JI=I1(JL)
       JJ=I2(JL)
       !calculation of w
       IF (PRSS(JI,JJ,JK) > ZRTMIN(5) ) THEN
         ZWSEDW1(JI,JJ,JK)=XFSEDS*(PRSS(JI,JJ,JK))**(XEXSEDS-1)*&
         PRHODREF(JI,JJ,JK)**(XEXSEDS-XCEXVT-1)
       ENDIF
       IF ( ZQP(JI,JJ) > ZRTMIN(5) ) THEN
         ZWSEDW2(JI,JJ,JK)=XFSEDS*(ZQP(JI,JJ))**(XEXSEDS-1)*&
         PRHODREF(JI,JJ,JK)**(XEXSEDS-XCEXVT-1)
       ENDIF
     ENDDO
     DO JJ = KJB, KJE
       DO JI = KIB, KIE
         ZH=PDZZ(JI,JJ,JK)
         ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH )
         IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
           ZP2 = MAX(0.,1 - ZH&
          / (PTSTEP*ZWSEDW2(JI,JJ,JK)) )
         ELSE
           ZP2 = 0.
         ENDIF
         ZWSED (JI,JJ,JK)=ZP1*PRHODREF(JI,JJ,JK)*&
         &ZH*PRSS(JI,JJ,JK)&
         &* ZINVTSTEP+ ZP2 * ZWSED (JI,JJ,JK+KKL)
       ENDDO
     ENDDO
   ENDDO

   DO JK = KKTB , KKTE
     PRSS(:,:,JK) = PRSS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
   ENDDO
   IF (PRESENT(PFPR)) THEN
     DO JK = KKTB , KKTE
       PFPR(:,:,JK,5)=ZWSED(:,:,JK)
     ENDDO
   ENDIF

   PINPRS(:,:) = ZWSED(:,:,KKB)/XRHOLW                        ! in m/s

   PRSS(:,:,:) = PRSS(:,:,:) * ZINVTSTEP


!
!*       2.5   for graupeln
!

   PRGS(:,:,:) = PRGS(:,:,:) * PTSTEP
   ZWSED(:,:,:) = 0.
   ZWSEDW1(:,:,:) = 0.
   ZWSEDW2(:,:,:) = 0.

! calculation of ZP1, ZP2 and sedimentation flux
   DO JK = KKE,  KKB, -1*KKL
     !estimation of q' taking into account incomming ZWSED
     ZQP(:,:)=ZWSED(:,:,JK+KKL)*ZW(:,:,JK)

     JCOUNT=COUNTJV((PRGS(:,:,JK) > ZRTMIN(6)) .OR. &
                    (ZQP(:,:) > ZRTMIN(6)),I1(:),I2(:))
     DO JL=1, JCOUNT
       JI=I1(JL)
       JJ=I2(JL)
       !calculation of w
       IF ( PRGS(JI,JJ,JK) > ZRTMIN(6) ) THEN
         ZWSEDW1 (JI,JJ,JK)= XFSEDG*(PRGS(JI,JJ,JK))**(XEXSEDG-1) * &
                                  PRHODREF(JI,JJ,JK)**(XEXSEDG-XCEXVT-1)
       ENDIF
       IF ( ZQP(JI,JJ) > ZRTMIN(6) ) THEN
         ZWSEDW2 (JI,JJ,JK)= XFSEDG*(ZQP(JI,JJ))**(XEXSEDG-1) * &
                                  PRHODREF(JI,JJ,JK)**(XEXSEDG-XCEXVT-1)
       ENDIF
     ENDDO
     DO JJ = KJB, KJE
       DO JI = KIB, KIE
         ZH=PDZZ(JI,JJ,JK)
         ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH )
         IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
           ZP2 = MAX(0.,1 - ZH &
         & / (PTSTEP*ZWSEDW2(JI,JJ,JK)) )
         ELSE
           ZP2 = 0.
         ENDIF
         ZWSED (JI,JJ,JK)=ZP1*PRHODREF(JI,JJ,JK)*&
         &ZH*PRGS(JI,JJ,JK)&
         &* ZINVTSTEP+ ZP2 * ZWSED (JI,JJ,JK+KKL)
       ENDDO
     ENDDO
   ENDDO

   DO JK = KKTB , KKTE
         PRGS(:,:,JK) = PRGS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
   ENDDO
   IF (PRESENT(PFPR)) THEN
     DO JK = KKTB , KKTE
       PFPR(:,:,JK,6)=ZWSED(:,:,JK)
     ENDDO
   ENDIF

   PINPRG(:,:) = ZWSED(:,:,KKB)/XRHOLW                        ! in m/s

   PRGS(:,:,:) = PRGS(:,:,:) * ZINVTSTEP

!
!*       2.6   for hail
!
 IF ( KRR == 7 ) THEN
     PRHS(:,:,:) = PRHS(:,:,:) * PTSTEP
     ZWSED(:,:,:) = 0.
     ZWSEDW1(:,:,:) = 0.
     ZWSEDW2(:,:,:) = 0.
! calculation of ZP1, ZP2 and sedimentation flux
     DO JK = KKE , KKB, -1*KKL
     !estimation of q' taking into account incomming ZWSED
     ZQP(:,:)=ZWSED(:,:,JK+KKL)*ZW(:,:,JK)

     JCOUNT=COUNTJV((PRHS(:,:,JK)+ZQP(JI,JJ) > ZRTMIN(7)) .OR. &
                    (ZQP(:,:) > ZRTMIN(7)),I1(:),I2(:))
     DO JL=1, JCOUNT
       JI=I1(JL)
       JJ=I2(JL)
         !calculation of w
         IF ((PRHS(JI,JJ,JK)+ZQP(JI,JJ)) > ZRTMIN(7) ) THEN
           ZWSEDW1 (JI,JJ,JK)= XFSEDH  * (PRHS(JI,JJ,JK))**(XEXSEDH-1) *   &
                                    PRHODREF(JI,JJ,JK)**(XEXSEDH-XCEXVT-1)
         ENDIF
         IF ( ZQP(JI,JJ) > ZRTMIN(7) ) THEN
           ZWSEDW2 (JI,JJ,JK)= XFSEDH  * ZQP(JI,JJ)**(XEXSEDH-1) *   &
                                    PRHODREF(JI,JJ,JK)**(XEXSEDH-XCEXVT-1)
         ENDIF
       ENDDO
       DO JJ = KJB, KJE
         DO JI = KIB, KIE
           ZH=PDZZ(JI,JJ,JK)
           ZP1 = MIN(1., ZWSEDW1(JI,JJ,JK) * PTSTEP / ZH)
           IF (ZWSEDW2(JI,JJ,JK) /= 0.) THEN
             ZP2 = MAX(0.,1 - ZH &
           &  / (PTSTEP*ZWSEDW2(JI,JJ,JK)) )
           ELSE
             ZP2 = 0.
           ENDIF
           ZWSED (JI,JJ,JK)=ZP1*PRHODREF(JI,JJ,JK)*&
           &ZH*PRHS(JI,JJ,JK)&
           &* ZINVTSTEP+ ZP2 * ZWSED (JI,JJ,JK+KKL)
         ENDDO
       ENDDO
     ENDDO

     DO JK = KKTB , KKTE
       PRHS(:,:,JK) = PRHS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
     ENDDO
     IF (PRESENT(PFPR)) THEN
       DO JK = KKTB , KKTE
         PFPR(:,:,JK,7)=ZWSED(:,:,JK)
       ENDDO
     ENDIF

     PINPRH(:,:) = ZWSED(:,:,KKB)/XRHOLW                        ! in m/s

     PRHS(:,:,:) = PRHS(:,:,:) * ZINVTSTEP

 ENDIF
!
!*       2.3     budget storage
!
if ( lbudget_rc .and. osedic ) call Budget_store_end( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rr )              call Budget_store_end( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_ri )              call Budget_store_end( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rs )              call Budget_store_end( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rg )              call Budget_store_end( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rh )              call Budget_store_end( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )
!
!
!*       2.4  DROPLET DEPOSITION AT THE 1ST LEVEL ABOVE GROUND
!
IF (LDEPOSC) THEN
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'DEPO', prcs(:, :, :) * prhodj(:, :, :) )

  GDEP(:,:) = .FALSE.
  GDEP(KIB:KIE,KJB:KJE) =    PRCS(KIB:KIE,KJB:KJE,KKB) >0
  WHERE (GDEP)
     PRCS(:,:,KKB) = PRCS(:,:,KKB) - XVDEPOSC * PRCT(:,:,KKB) / PDZZ(:,:,KKB)
     PINPRC(:,:) = PINPRC(:,:) + XVDEPOSC * PRCT(:,:,KKB) * PRHODREF(:,:,KKB) /XRHOLW
     PINDEP(:,:) = XVDEPOSC * PRCT(:,:,KKB) * PRHODREF(:,:,KKB) /XRHOLW
  END WHERE

  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'DEPO', prcs(:, :, :) * prhodj(:, :, :) )
END IF

END SUBROUTINE RAIN_ICE_SEDIMENTATION_STAT

END MODULE MODE_RAIN_ICE_SEDIMENTATION_STAT
