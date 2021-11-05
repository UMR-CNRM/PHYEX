!MNH_LIC Copyright 1995-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 25/02/2019: split rain_ice (cleaner and easier to maintain/debug)
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_SEDIMENTATION_SPLIT

  IMPLICIT NONE

  PRIVATE

  PUBLIC RAIN_ICE_SEDIMENTATION_SPLIT

CONTAINS

SUBROUTINE RAIN_ICE_SEDIMENTATION_SPLIT(KIB, KIE, KJB, KJE, KKB, KKE, KKTB, KKTE, KKT, KKL,&
  KSPLITR,PTSTEP, &
  KRR,OSEDIC,ODEPOSC,PINPRC,PINDEP,PINPRR,PINPRS,PINPRG,PDZZ,PRHODREF,PPABST,PTHT,PRHODJ,&
      PINPRR3D,PRCS,PRCT,PRRS,PRRT,PRIS,PRIT,PRSS,PRST,PRGS,PRGT,PSEA,PTOWN,PINPRH,PRHS,PRHT,PFPR)
!
!*      0. DECLARATIONS
!          ------------
!
use modd_budget,         only: lbudget_rc, lbudget_rr, lbudget_ri, lbudget_rs, lbudget_rg, lbudget_rh, &
                               NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, &
                               tbudgets
use MODD_CST,            only: XCPD, XP00, XRD, XRHOLW
use MODD_PARAM_ICE,      only: XVDEPOSC
use MODD_RAIN_ICE_DESCR, only: XCC, XCONC_LAND, xconc_sea, xconc_urban, XDC, XCEXVT, &
                               XALPHAC, XNUC, XALPHAC2, XNUC2, XLBEXC, XRTMIN, XLBEXC, XLBC
use MODD_RAIN_ICE_PARAM, only: XEXSEDG, XEXSEDH, XEXCSEDI, XEXSEDR, XEXSEDS, &
                               XFSEDG, XFSEDH, XFSEDI, XFSEDR, XFSEDS, XFSEDC

use mode_budget,         only: Budget_store_init, Budget_store_end
use mode_tools,          only: Countjv

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                            INTENT(IN)    :: KIB, KIE, KJB, KJE, KKB, KKE, KKTB, KKTE, KKT
INTEGER,                            INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
INTEGER,                            INTENT(IN)    :: KSPLITR ! Number of small time step
                                                             ! integration for  rain sedimendation
REAL,                               INTENT(IN)    :: PTSTEP  ! Double Time step
                                                            ! (single if cold start)
INTEGER,                            INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL,                            INTENT(IN)    :: OSEDIC  ! Switch for droplet sedim.
LOGICAL,                            INTENT(IN)    :: ODEPOSC ! Switch for droplet depos.
REAL, DIMENSION(:,:),               INTENT(INOUT) :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(:,:),               INTENT(INOUT) :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(:,:),               INTENT(OUT)   :: PINPRR  ! Rain instant precip
REAL, DIMENSION(:,:),               INTENT(OUT)   :: PINPRS  ! Snow instant precip
REAL, DIMENSION(:,:),               INTENT(OUT)   :: PINPRG  ! Graupel instant precip
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
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),             INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),             INTENT(INOUT) :: PRGS    ! Graupel m.r. source
REAL, DIMENSION(:,:,:),             INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
!
!*       0.2  declaration of local variables
!
!
INTEGER, SAVE                      :: IOLDALLOCC = 6000
INTEGER, SAVE                      :: IOLDALLOCR = 6000
INTEGER, SAVE                      :: IOLDALLOCI = 6000
INTEGER, SAVE                      :: IOLDALLOCS = 6000
INTEGER, SAVE                      :: IOLDALLOCG = 6000
INTEGER, SAVE                      :: IOLDALLOCH = 6000
INTEGER                            :: ILENALLOCC,ILENALLOCR,ILENALLOCI,ILENALLOCS,ILENALLOCG,ILENALLOCH
INTEGER                            :: ILISTLENC,ILISTLENR,ILISTLENI,ILISTLENS,ILISTLENG,ILISTLENH
INTEGER                            :: ISEDIMR,ISEDIMC, ISEDIMI, ISEDIMS, ISEDIMG, ISEDIMH
INTEGER                            :: JK            ! Vertical loop index for the rain sedimentation
INTEGER                            :: JN            ! Temporal loop index for the rain sedimentation
INTEGER                            :: JJ            ! Loop index for the interpolation
INTEGER                            :: JL
INTEGER, DIMENSION(SIZE(PRCS))     :: IC1,IC2,IC3 ! Used to replace the COUNT
INTEGER, DIMENSION(SIZE(PRCS))     :: IR1,IR2,IR3 ! Used to replace the COUNT
INTEGER, DIMENSION(SIZE(PRCS))     :: IS1,IS2,IS3 ! Used to replace the COUNT
INTEGER, DIMENSION(SIZE(PRCS))     :: II1,II2,II3 ! Used to replace the COUNT
INTEGER, DIMENSION(SIZE(PRCS))     :: IG1,IG2,IG3 ! Used to replace the COUNT
INTEGER, DIMENSION(SIZE(PRCS))     :: IH1,IH2,IH3 ! Used to replace the COUNT
INTEGER, DIMENSION(:), ALLOCATABLE :: ILISTR,ILISTC,ILISTI,ILISTS,ILISTG,ILISTH
LOGICAL, DIMENSION(SIZE(PRCS,1),SIZE(PRCS,2)):: GDEP
LOGICAL, DIMENSION(SIZE(PRCS,1),SIZE(PRCS,2),SIZE(PRCS,3)) &
                                   :: GSEDIMR,GSEDIMC, GSEDIMI, GSEDIMS, GSEDIMG, GSEDIMH ! Test where to compute the SED processes
REAL                               :: ZINVTSTEP
REAL                               :: ZTSPLITR ! Small time step for rain sedimentation
REAL,    DIMENSION(SIZE(XRTMIN))   :: ZRTMIN
! XRTMIN = Minimum value for the mixing ratio
! ZRTMIN = Minimum value for the source (tendency)
REAL,    DIMENSION(:), ALLOCATABLE :: ZRCS    ! Cloud water m.r. source
REAL,    DIMENSION(:), ALLOCATABLE :: ZRRS    ! Rain water m.r. source
REAL,    DIMENSION(:), ALLOCATABLE :: ZRIS    ! Pristine ice m.r. source
REAL,    DIMENSION(:), ALLOCATABLE :: ZRSS    ! Snow/aggregate m.r. source
REAL,    DIMENSION(:), ALLOCATABLE :: ZRGS    ! Graupel m.r. source
REAL,    DIMENSION(:), ALLOCATABLE :: ZRHS    ! Hail m.r. source
REAL,    DIMENSION(:), ALLOCATABLE :: ZRCT    ! Cloud water m.r. at t
REAL,    DIMENSION(:), ALLOCATABLE :: ZRHODREFC,& ! RHO Dry REFerence
                                      ZRHODREFR,& ! RHO Dry REFerence
                                      ZRHODREFI,& ! RHO Dry REFerence
                                      ZRHODREFS,& ! RHO Dry REFerence
                                      ZRHODREFG,& ! RHO Dry REFerence
                                      ZRHODREFH,& ! RHO Dry REFerence
                                      ZCC,      & ! terminal velocity
                                      ZFSEDC1D, & ! For cloud sedimentation
                                      ZWLBDC,   & ! Slope parameter of the droplet  distribution
                                      ZCONC,    & ! Concentration des aerosols
                                      ZRAY1D,   & ! Mean radius
                                      ZWLBDA,   & ! Libre parcours moyen
                                      ZZT,      & ! Temperature
                                      ZPRES       ! Pressure
REAL,    DIMENSION(SIZE(PRCS,1),SIZE(PRCS,2))                   &
                                   :: ZCONC_TMP    ! Weighted concentration
REAL,    DIMENSION(SIZE(PRCS,1),SIZE(PRCS,2),SIZE(PRCS,3)) :: ZCONC3D !  droplet condensation
REAL,    DIMENSION(SIZE(PRCS,1),SIZE(PRCS,2),SIZE(PRCS,3)) ::  &
                                      ZRAY,   & ! Cloud Mean radius
                                      ZLBC,   & ! XLBC weighted by sea fraction
                                      ZFSEDC
REAL,    DIMENSION(SIZE(PRCS,1),SIZE(PRCS,2),SIZE(PRCS,3))   &
                                   :: ZPRCS,ZPRRS,ZPRSS,ZPRGS,ZPRHS   ! Mixing ratios created during the time step
REAL,    DIMENSION(SIZE(PRCS,1),SIZE(PRCS,2),SIZE(PRCS,3))   &
                                   :: ZW ! work array
REAL,    DIMENSION(SIZE(PRCS,1),SIZE(PRCS,2),0:SIZE(PRCS,3)+1)   &
                                   :: ZWSED        ! sedimentation fluxes
!-------------------------------------------------------------------------------

if ( lbudget_rc .and. osedic ) call Budget_store_init( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rr )              call Budget_store_init( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_ri )              call Budget_store_init( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rs )              call Budget_store_init( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rg )              call Budget_store_init( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rh )              call Budget_store_init( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )
!
!        O. Initialization of for sedimentation
!
ZINVTSTEP=1./PTSTEP
ZTSPLITR= PTSTEP / REAL(KSPLITR)
!
IF (OSEDIC) PINPRC (:,:) = 0.
IF (ODEPOSC) PINDEP (:,:) = 0.
PINPRR (:,:) = 0.
PINPRR3D (:,:,:) = 0.
PINPRS (:,:) = 0.
PINPRG (:,:) = 0.
IF ( KRR == 7 ) PINPRH (:,:) = 0.
IF (PRESENT(PFPR)) PFPR(:,:,:,:) = 0.
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
!
!*       2.    compute the fluxes
!
!  optimization by looking for locations where
!  the precipitating fields are larger than a minimal value only !!!
!  For optimization we consider each variable separately

ZRTMIN(:)    = XRTMIN(:) * ZINVTSTEP
IF (OSEDIC) GSEDIMC(:,:,:) = .FALSE.
GSEDIMR(:,:,:) = .FALSE.
GSEDIMI(:,:,:) = .FALSE.
GSEDIMS(:,:,:) = .FALSE.
GSEDIMG(:,:,:) = .FALSE.
IF ( KRR == 7 ) GSEDIMH(:,:,:) = .FALSE.
!
ILENALLOCR = 0
IF (OSEDIC) ILENALLOCC = 0
ILENALLOCI = 0
ILENALLOCS = 0
ILENALLOCG = 0
IF ( KRR == 7 ) ILENALLOCH = 0
!
! ZPiS = Specie i source creating during the current time step
! PRiS = Source of the previous time step
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
! PRiS = Source of the previous time step + source created during the subtime
! step
!
DO JN = 1 , KSPLITR
  IF( JN==1 ) THEN
   IF (OSEDIC) PRCS(:,:,:) = PRCS(:,:,:) + ZPRCS(:,:,:)/KSPLITR
   PRRS(:,:,:) = PRRS(:,:,:) + ZPRRS(:,:,:)/KSPLITR
   PRSS(:,:,:) = PRSS(:,:,:) + ZPRSS(:,:,:)/KSPLITR
   PRGS(:,:,:) = PRGS(:,:,:) + ZPRGS(:,:,:)/KSPLITR
   IF ( KRR == 7 ) PRHS(:,:,:) = PRHS(:,:,:) + ZPRHS(:,:,:)/KSPLITR
   DO JK = KKTB , KKTE
     ZW(:,:,JK) =ZTSPLITR/(PRHODREF(:,:,JK)* PDZZ(:,:,JK))
   END DO
 ELSE
   IF (OSEDIC) PRCS(:,:,:) = PRCS(:,:,:) + ZPRCS(:,:,:)*ZTSPLITR
   PRRS(:,:,:) = PRRS(:,:,:) + ZPRRS(:,:,:)*ZTSPLITR
   PRSS(:,:,:) = PRSS(:,:,:) + ZPRSS(:,:,:)*ZTSPLITR
   PRGS(:,:,:) = PRGS(:,:,:) + ZPRGS(:,:,:)*ZTSPLITR
   IF ( KRR == 7 ) PRHS(:,:,:) = PRHS(:,:,:) + ZPRHS(:,:,:)*ZTSPLITR
 END IF
 !
 IF (OSEDIC) GSEDIMC(KIB:KIE,KJB:KJE,KKTB:KKTE) =                &
                  PRCS(KIB:KIE,KJB:KJE,KKTB:KKTE)>ZRTMIN(2)
 GSEDIMR(KIB:KIE,KJB:KJE,KKTB:KKTE) =                            &
                  PRRS(KIB:KIE,KJB:KJE,KKTB:KKTE)>ZRTMIN(3)
 GSEDIMI(KIB:KIE,KJB:KJE,KKTB:KKTE) =                            &
                  PRIS(KIB:KIE,KJB:KJE,KKTB:KKTE)>ZRTMIN(4)
 GSEDIMS(KIB:KIE,KJB:KJE,KKTB:KKTE) =                            &
                  PRSS(KIB:KIE,KJB:KJE,KKTB:KKTE)>ZRTMIN(5)
 GSEDIMG(KIB:KIE,KJB:KJE,KKTB:KKTE) =                            &
                  PRGS(KIB:KIE,KJB:KJE,KKTB:KKTE)>ZRTMIN(6)
 IF ( KRR == 7 ) GSEDIMH(KIB:KIE,KJB:KJE,KKTB:KKTE) =            &
                  PRHS(KIB:KIE,KJB:KJE,KKTB:KKTE)>ZRTMIN(7)
!
 IF (OSEDIC) ISEDIMC = COUNTJV( GSEDIMC(:,:,:),IC1(:),IC2(:),IC3(:))
 ISEDIMR = COUNTJV( GSEDIMR(:,:,:),IR1(:),IR2(:),IR3(:))
 ISEDIMI = COUNTJV( GSEDIMI(:,:,:),II1(:),II2(:),II3(:))
 ISEDIMS = COUNTJV( GSEDIMS(:,:,:),IS1(:),IS2(:),IS3(:))
 ISEDIMG = COUNTJV( GSEDIMG(:,:,:),IG1(:),IG2(:),IG3(:))
 IF ( KRR == 7 ) ISEDIMH = COUNTJV( GSEDIMH(:,:,:),IH1(:),IH2(:),IH3(:))
!
!*       2.1   for cloud
!
 IF (OSEDIC) THEN
  ZWSED(:,:,:) = 0.
  IF( JN==1 ) PRCS(:,:,:) = PRCS(:,:,:) * PTSTEP
  IF( ISEDIMC >= 1 ) THEN
    IF ( ISEDIMC .GT. ILENALLOCC ) THEN
      IF ( ILENALLOCC .GT. 0 ) THEN
        DEALLOCATE (ZRCS, ZRHODREFC, ILISTC,ZWLBDC,ZCONC,ZRCT,  &
                    ZZT,ZPRES,ZRAY1D,ZFSEDC1D,ZWLBDA,ZCC )
      END IF
      ILENALLOCC = MAX (IOLDALLOCC, 2*ISEDIMC )
      IOLDALLOCC = ILENALLOCC
      ALLOCATE(ZRCS(ILENALLOCC), ZRHODREFC(ILENALLOCC), ILISTC(ILENALLOCC), &
        ZWLBDC(ILENALLOCC), ZCONC(ILENALLOCC), ZRCT(ILENALLOCC), ZZT(ILENALLOCC), &
        ZPRES(ILENALLOCC), ZRAY1D(ILENALLOCC), ZFSEDC1D(ILENALLOCC), &
        ZWLBDA(ILENALLOCC), ZCC(ILENALLOCC)  )
    END IF
!
    DO JL=1,ISEDIMC
      ZRCS(JL) = PRCS(IC1(JL),IC2(JL),IC3(JL))
      ZRHODREFC(JL) =  PRHODREF(IC1(JL),IC2(JL),IC3(JL))
      ZWLBDC(JL) = ZLBC(IC1(JL),IC2(JL),IC3(JL))
      ZCONC(JL) = ZCONC3D(IC1(JL),IC2(JL),IC3(JL))
      ZRCT(JL) = PRCT(IC1(JL),IC2(JL),IC3(JL))
      ZZT(JL) = PTHT(IC1(JL),IC2(JL),IC3(JL))
      ZPRES(JL) = PPABST(IC1(JL),IC2(JL),IC3(JL))
      ZRAY1D(JL) = ZRAY(IC1(JL),IC2(JL),IC3(JL))
      ZFSEDC1D(JL) = ZFSEDC(IC1(JL),IC2(JL),IC3(JL))
    END DO
!
    ILISTLENC = 0
    DO JL=1,ISEDIMC
     IF( ZRCS(JL) .GT. ZRTMIN(2) ) THEN
       ILISTLENC = ILISTLENC + 1
       ILISTC(ILISTLENC) = JL
     END IF
    END DO
       DO JJ = 1, ILISTLENC
          JL = ILISTC(JJ)
          IF (ZRCS(JL) .GT. ZRTMIN(2) .AND. ZRCT(JL) .GT. XRTMIN(2)) THEN
            ZWLBDC(JL) = ZWLBDC(JL) * ZCONC(JL) / (ZRHODREFC(JL) * ZRCT(JL))
            ZWLBDC(JL) = ZWLBDC(JL)**XLBEXC
            ZRAY1D(JL) = ZRAY1D(JL) / ZWLBDC(JL) !! ZRAY : mean diameter=M(1)/2
            ZZT(JL)    = ZZT(JL) * (ZPRES(JL)/XP00)**(XRD/XCPD)
            ZWLBDA(JL) = 6.6E-8*(101325./ZPRES(JL))*(ZZT(JL)/293.15)
            ZCC(JL)    = XCC*(1.+1.26*ZWLBDA(JL)/ZRAY1D(JL)) !! XCC modified for cloud
            ZWSED (IC1(JL),IC2(JL),IC3(JL))= ZRHODREFC(JL)**(-XCEXVT +1 ) *   &
              ZWLBDC(JL)**(-XDC)*ZCC(JL)*ZFSEDC1D(JL) * ZRCS(JL)
          END IF
       END DO
  END IF
       DO JK = KKTB , KKTE
         PRCS(:,:,JK) = PRCS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
       END DO
       IF (PRESENT(PFPR)) THEN
         DO JK = KKTB , KKTE
           PFPR(:,:,JK,2)=ZWSED(:,:,JK)
         ENDDO
       ENDIF
      PINPRC(:,:) = PINPRC(:,:) + ZWSED(:,:,KKB) / XRHOLW / KSPLITR
      IF( JN==KSPLITR ) THEN
        PRCS(:,:,:) = PRCS(:,:,:) * ZINVTSTEP
      END IF
 END IF
!
!*       2.2   for rain
!
  IF( JN==1 ) PRRS(:,:,:) = PRRS(:,:,:) * PTSTEP
  ZWSED(:,:,:) = 0.
  IF( ISEDIMR >= 1 ) THEN
    IF ( ISEDIMR .GT. ILENALLOCR ) THEN
      IF ( ILENALLOCR .GT. 0 ) THEN
        DEALLOCATE (ZRRS, ZRHODREFR, ILISTR)
      END IF
      ILENALLOCR = MAX (IOLDALLOCR, 2*ISEDIMR )
      IOLDALLOCR = ILENALLOCR
      ALLOCATE(ZRRS(ILENALLOCR), ZRHODREFR(ILENALLOCR), ILISTR(ILENALLOCR))
    END IF
!
    DO JL=1,ISEDIMR
      ZRRS(JL) = PRRS(IR1(JL),IR2(JL),IR3(JL))
      ZRHODREFR(JL) =  PRHODREF(IR1(JL),IR2(JL),IR3(JL))
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
           ZWSED (IR1(JL),IR2(JL),IR3(JL))= XFSEDR  * ZRRS(JL)**XEXSEDR *   &
                                        ZRHODREFR(JL)**(XEXSEDR-XCEXVT)
       END DO
  END IF
       DO JK = KKTB , KKTE
         PRRS(:,:,JK) = PRRS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
       END DO
       IF (PRESENT(PFPR)) THEN
         DO JK = KKTB , KKTE
           PFPR(:,:,JK,3)=ZWSED(:,:,JK)
         ENDDO
       ENDIF
       PINPRR(:,:) = PINPRR(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
       PINPRR3D(:,:,:) = PINPRR3D(:,:,:) + ZWSED(:,:,1:KKT)/XRHOLW/KSPLITR
      IF( JN==KSPLITR ) THEN
        PRRS(:,:,:) = PRRS(:,:,:) * ZINVTSTEP
      END IF
!
!*       2.3   for pristine ice
!
  IF( JN==1 ) PRIS(:,:,:) = PRIS(:,:,:) * PTSTEP
  ZWSED(:,:,:) = 0.
  IF( ISEDIMI >= 1 ) THEN
    IF ( ISEDIMI .GT. ILENALLOCI ) THEN
      IF ( ILENALLOCI .GT. 0 ) THEN
        DEALLOCATE (ZRIS, ZRHODREFI, ILISTI)
      END IF
      ILENALLOCI = MAX (IOLDALLOCI, 2*ISEDIMI )
      IOLDALLOCI = ILENALLOCI
      ALLOCATE(ZRIS(ILENALLOCI), ZRHODREFI(ILENALLOCI), ILISTI(ILENALLOCI))
    END IF
!
    DO JL=1,ISEDIMI
      ZRIS(JL) = PRIS(II1(JL),II2(JL),II3(JL))
      ZRHODREFI(JL) =  PRHODREF(II1(JL),II2(JL),II3(JL))
    END DO
!
    ILISTLENI = 0
    DO JL=1,ISEDIMI
     IF( ZRIS(JL) .GT.  MAX(ZRTMIN(4),1.0E-7 )) THEN ! limitation of the McF&H formula
       ILISTLENI = ILISTLENI + 1
       ILISTI(ILISTLENI) = JL
     END IF
    END DO
       DO JJ = 1, ILISTLENI
          JL = ILISTI(JJ)
              ZWSED (II1(JL),II2(JL),II3(JL))= XFSEDI * ZRIS(JL) *  &
                               ZRHODREFI(JL)**(1.0-XCEXVT) * & !    McF&H
                               MAX( 0.05E6,-0.15319E6-0.021454E6* &
                               ALOG(ZRHODREFI(JL)*ZRIS(JL)) )**XEXCSEDI
       END DO
  END IF
       DO JK = KKTB , KKTE
         PRIS(:,:,JK) = PRIS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
       END DO
       IF (PRESENT(PFPR)) THEN
         DO JK = KKTB , KKTE
           PFPR(:,:,JK,4)=ZWSED(:,:,JK)
         ENDDO
       ENDIF
      IF( JN==KSPLITR ) THEN
        PRIS(:,:,:) = PRIS(:,:,:) * ZINVTSTEP
      END IF
!
!*       2.4   for aggregates/snow
!
  IF( JN==1 ) PRSS(:,:,:) = PRSS(:,:,:) * PTSTEP
  ZWSED(:,:,:) = 0.
  IF( ISEDIMS >= 1 ) THEN
    IF ( ISEDIMS .GT. ILENALLOCS ) THEN
      IF ( ILENALLOCS .GT. 0 ) THEN
        DEALLOCATE (ZRSS, ZRHODREFS, ILISTS)
      END IF
      ILENALLOCS = MAX (IOLDALLOCS, 2*ISEDIMS )
      IOLDALLOCS = ILENALLOCS
      ALLOCATE(ZRSS(ILENALLOCS), ZRHODREFS(ILENALLOCS), ILISTS(ILENALLOCS))
    END IF
!
    DO JL=1,ISEDIMS
      ZRSS(JL) = PRSS(IS1(JL),IS2(JL),IS3(JL))
      ZRHODREFS(JL) =  PRHODREF(IS1(JL),IS2(JL),IS3(JL))
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
             ZWSED (IS1(JL),IS2(JL),IS3(JL))= XFSEDS * ZRSS(JL)**XEXSEDS *  &
                                        ZRHODREFS(JL)**(XEXSEDS-XCEXVT)
       END DO
  END IF
       DO JK = KKTB , KKTE
         PRSS(:,:,JK) = PRSS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
       END DO
       IF (PRESENT(PFPR)) THEN
         DO JK = KKTB , KKTE
           PFPR(:,:,JK,5)=ZWSED(:,:,JK)
         ENDDO
       ENDIF
      PINPRS(:,:) = PINPRS(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
      IF( JN==KSPLITR ) THEN
        PRSS(:,:,:) = PRSS(:,:,:) * ZINVTSTEP
      END IF
!
!*       2.5   for graupeln
!
  ZWSED(:,:,:) = 0.
  IF( JN==1 ) PRGS(:,:,:) = PRGS(:,:,:) * PTSTEP
  IF( ISEDIMG >= 1 ) THEN
    IF ( ISEDIMG .GT. ILENALLOCG ) THEN
      IF ( ILENALLOCG .GT. 0 ) THEN
        DEALLOCATE (ZRGS, ZRHODREFG, ILISTG)
      END IF
      ILENALLOCG = MAX (IOLDALLOCG, 2*ISEDIMG )
      IOLDALLOCG = ILENALLOCG
      ALLOCATE(ZRGS(ILENALLOCG), ZRHODREFG(ILENALLOCG), ILISTG(ILENALLOCG))
    END IF
!
    DO JL=1,ISEDIMG
      ZRGS(JL) = PRGS(IG1(JL),IG2(JL),IG3(JL))
      ZRHODREFG(JL) =  PRHODREF(IG1(JL),IG2(JL),IG3(JL))
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
             ZWSED (IG1(JL),IG2(JL),IG3(JL))= XFSEDG  * ZRGS(JL)**XEXSEDG *   &
                                        ZRHODREFG(JL)**(XEXSEDG-XCEXVT)
       END DO
END IF
       DO JK = KKTB , KKTE
         PRGS(:,:,JK) = PRGS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
       END DO
       IF (PRESENT(PFPR)) THEN
         DO JK = KKTB , KKTE
           PFPR(:,:,JK,6)=ZWSED(:,:,JK)
         ENDDO
       ENDIF
       PINPRG(:,:) = PINPRG(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
      IF( JN==KSPLITR ) THEN
        PRGS(:,:,:) = PRGS(:,:,:) * ZINVTSTEP
      END IF
!
!*       2.6   for hail
!
 IF ( KRR == 7 ) THEN
  IF( JN==1 ) PRHS(:,:,:) = PRHS(:,:,:) * PTSTEP
  ZWSED(:,:,:) = 0.
  IF( ISEDIMH >= 1 ) THEN
    IF ( ISEDIMH .GT. ILENALLOCH ) THEN
      IF ( ILENALLOCH .GT. 0 ) THEN
        DEALLOCATE (ZRHS, ZRHODREFH, ILISTH)
      END IF
      ILENALLOCH = MAX (IOLDALLOCH, 2*ISEDIMH )
      IOLDALLOCH = ILENALLOCH
      ALLOCATE(ZRHS(ILENALLOCH), ZRHODREFH(ILENALLOCH), ILISTH(ILENALLOCH))
    END IF
!
    DO JL=1,ISEDIMH
      ZRHS(JL) = PRHS(IH1(JL),IH2(JL),IH3(JL))
      ZRHODREFH(JL) =  PRHODREF(IH1(JL),IH2(JL),IH3(JL))
    END DO
!
    ILISTLENH = 0
    DO JL=1,ISEDIMH
     IF( ZRHS(JL) .GT. ZRTMIN(7) ) THEN
       ILISTLENH = ILISTLENH + 1
       ILISTH(ILISTLENH) = JL
     END IF
    END DO
       DO JJ = 1, ILISTLENH
          JL = ILISTH(JJ)
             ZWSED (IH1(JL),IH2(JL),IH3(JL))= XFSEDH  * ZRHS(JL)**XEXSEDH *   &
                                        ZRHODREFH(JL)**(XEXSEDH-XCEXVT)
       END DO
  END IF
       DO JK = KKTB , KKTE
         PRHS(:,:,JK) = PRHS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
       END DO
       IF (PRESENT(PFPR)) THEN
         DO JK = KKTB , KKTE
           PFPR(:,:,JK,7)=ZWSED(:,:,JK)
         ENDDO
       ENDIF
       PINPRH(:,:) = PINPRH(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
      IF( JN==KSPLITR ) THEN
        PRHS(:,:,:) = PRHS(:,:,:) * ZINVTSTEP
      END IF
 END IF
!
END DO
!
IF (OSEDIC) THEN
   IF (ILENALLOCC .GT. 0) DEALLOCATE (ZRCS, ZRHODREFC,  &
  ILISTC,ZWLBDC,ZCONC,ZRCT, ZZT,ZPRES,ZRAY1D,ZFSEDC1D, ZWLBDA,ZCC)
END IF
IF (ILENALLOCR .GT. 0 ) DEALLOCATE(ZRHODREFR,ZRRS,ILISTR)
IF (ILENALLOCI .GT. 0 ) DEALLOCATE(ZRHODREFI,ZRIS,ILISTI)
IF (ILENALLOCS .GT. 0 ) DEALLOCATE(ZRHODREFS,ZRSS,ILISTS)
IF (ILENALLOCG .GT. 0 ) DEALLOCATE(ZRHODREFG,ZRGS,ILISTG)
IF (KRR == 7 .AND. (ILENALLOCH .GT. 0 )) DEALLOCATE(ZRHODREFH,ZRHS,ILISTH)
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
!*       2.4  DROPLET DEPOSITION AT THE 1ST LEVEL ABOVE GROUND
!
IF (ODEPOSC) THEN
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

END SUBROUTINE RAIN_ICE_SEDIMENTATION_SPLIT

END MODULE MODE_RAIN_ICE_SEDIMENTATION_SPLIT
