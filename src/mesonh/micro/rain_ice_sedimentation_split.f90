!MNH_LIC Copyright 1995-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 25/02/2019: split rain_ice (cleaner and easier to maintain/debug)
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet 25/06/2019: OpenACC: optimisation of rain_ice_sedimentation_split
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_SEDIMENTATION_SPLIT

  IMPLICIT NONE

  PRIVATE

  PUBLIC RAIN_ICE_SEDIMENTATION_SPLIT

CONTAINS

SUBROUTINE RAIN_ICE_SEDIMENTATION_SPLIT( KIB, KIE, KJB, KJE, KKB, KKE, KKTB, KKTE, KKT, KKL,                           &
                                         KSPLITR, PTSTEP, KRR, OSEDIC, ODEPOSC,                                        &
                                         PINPRC, PINDEP, PINPRR, PINPRS, PINPRG, PDZZ, PRHODREF, PPABST, PTHT, PRHODJ, &
                                         PINPRR3D, PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,         &
                                         PSEA, PTOWN, PINPRH, PRHS, PRHT, PFPR                                         )
!
!*      0. DECLARATIONS
!          ------------
!
use modd_budget,         only: lbudget_rc, lbudget_rr, lbudget_ri, lbudget_rs, lbudget_rg, lbudget_rh, &
                               NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, &
                               tbudgets
USE MODD_CST,            only: XCPD, XP00, XRD, XRHOLW
USE MODD_PARAM_ICE_n,      only: XVDEPOSC
USE MODD_RAIN_ICE_DESCR_n, only: XCC, XCONC_LAND, xconc_sea, xconc_urban, XDC, XCEXVT, &
                               XALPHAC, XNUC, XALPHAC2, XNUC2, XLBEXC, XRTMIN, XLBEXC, XLBC
USE MODD_RAIN_ICE_PARAM_n, only: XEXSEDG, XEXSEDH, XEXCSEDI, XEXSEDR, XEXSEDS, &
                               XFSEDG, XFSEDH, XFSEDI, XFSEDR, XFSEDS, XFSEDC

use mode_budget,         only: Budget_store_init, Budget_store_end
use mode_mppdb
#ifndef MNH_OPENACC
use mode_tools,           only: Countjv
#else
use mode_tools,           only: Countjv_device
#endif

#if defined(MNH_BITREP) || defined(MNH_BITREP_OMP)
USE MODI_BITREP
#endif

#ifdef MNH_OPENACC
USE MODE_MNH_ZWORK,      ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
#endif

#if defined(MNH_COMPILER_CCE) && defined(MNH_BITREP_OMP)
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif

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
INTEGER                                :: ISEDIMR, ISEDIMC, ISEDIMI, ISEDIMS, ISEDIMG, ISEDIMH
INTEGER                                :: JI, JJ, JK    ! Loop indices on grid
INTEGER                                :: JN            ! Temporal loop index for the rain sedimentation
INTEGER                                :: JL
INTEGER, DIMENSION(:),     POINTER, CONTIGUOUS :: IC1, IC2, IC3 ! Used to replace the COUNT
INTEGER, DIMENSION(:),     POINTER, CONTIGUOUS :: IR1, IR2, IR3 ! Used to replace the COUNT
INTEGER, DIMENSION(:),     POINTER, CONTIGUOUS :: IS1, IS2, IS3 ! Used to replace the COUNT
INTEGER, DIMENSION(:),     POINTER, CONTIGUOUS :: II1, II2, II3 ! Used to replace the COUNT
INTEGER, DIMENSION(:),     POINTER, CONTIGUOUS :: IG1, IG2, IG3 ! Used to replace the COUNT
INTEGER, DIMENSION(:),     POINTER, CONTIGUOUS :: IH1, IH2, IH3 ! Used to replace the COUNT

LOGICAL                                :: GPRESENT_PFPR, GPRESENT_PSEA
LOGICAL, DIMENSION(:,:),   POINTER, CONTIGUOUS :: GDEP
LOGICAL, DIMENSION(:,:,:), POINTER, CONTIGUOUS :: GSEDIMR, GSEDIMC, GSEDIMI, GSEDIMS, GSEDIMG, GSEDIMH ! Where to compute the SED processes
REAL                                   :: ZINVTSTEP
REAL                                   :: ZTSPLITR            ! Small time step for rain sedimentation
REAL                                   :: ZTMP1, ZTMP2, ZTMP3 ! Intermediate variables
REAL                                   :: ZRHODREFLOC         ! RHO Dry REFerence
REAL                                   :: ZRSLOC, ZRTLOC      ! Intermediate variables
REAL,    DIMENSION(:), POINTER, CONTIGUOUS     :: ZRTMIN

! XRTMIN = Minimum value for the mixing ratio
! ZRTMIN = Minimum value for the source (tendency)
REAL                                   :: ZCC,      & ! terminal velocity
                                          ZFSEDC1D, & ! For cloud sedimentation
                                          ZWLBDC,   & ! Slope parameter of the droplet  distribution
                                          ZCONC,    & ! Concentration des aerosols
                                          ZRAY1D,   & ! Mean radius
                                          ZWLBDA,   & ! Free mean path
                                          ZZT,      & ! Temperature
                                          ZPRES       ! Pressure
REAL,    DIMENSION(:,:),   POINTER, CONTIGUOUS :: ZCONC_TMP   ! Weighted concentration
REAL,    DIMENSION(:,:,:), POINTER, CONTIGUOUS :: ZCONC3D     ! Doplet condensation
REAL,    DIMENSION(:,:),   POINTER, CONTIGUOUS :: ZOMPSEA,ZTMP1_2D,ZTMP2_2D,ZTMP3_2D,ZTMP4_2D !Work arrays
REAL,    DIMENSION(:,:,:), POINTER, CONTIGUOUS :: ZRAY,   &   ! Cloud Mean radius
                                          ZLBC,   &   ! XLBC weighted by sea fraction
                                          ZFSEDC

REAL,    DIMENSION(:,:,:), POINTER, CONTIGUOUS :: ZPRCS,ZPRRS,ZPRSS,ZPRGS,ZPRHS   ! Mixing ratios created during the time step
REAL,    DIMENSION(:,:,:), POINTER, CONTIGUOUS :: ZW          ! Work array
REAL,    DIMENSION(:,:,:), POINTER, CONTIGUOUS :: ZWSED       ! sedimentation fluxes

INTEGER :: IIU,IJU,IKU, IIJKU
LOGICAL :: GKRR_7,GSEDIC
!
!-------------------------------------------------------------------------------
!
! IN variables
!
!$acc data present( PDZZ, PRHODREF, PPABST, PTHT, PRHODJ, PRCT, PRRT, PRIT, PRST, PRGT, &
!$acc &             PSEA, PTOWN, PRHT,                                                  &
!
! INOUT variables
!
!$acc &             PINPRC, PINDEP, PRCS, PRRS, PRIS, PRSS, PRGS, PRHS,                 &
!
! OUT variables
!
!$acc &             PINPRR, PINPRS, PINPRG, PINPRR3D, PINPRH, PFPR )                    &

!$acc      present( XRTMIN )

! !$acc &     copyin( XFSEDC, XLBC, XLBEXC )

IF (MPPDB_INITIALIZED) THEN
  !Check all IN arrays
  CALL MPPDB_CHECK(PDZZ,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PDZZ")
  CALL MPPDB_CHECK(PRHODREF,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRHODREF")
  CALL MPPDB_CHECK(PPABST,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PPABST")
  CALL MPPDB_CHECK(PTHT,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PTHT")
  CALL MPPDB_CHECK(PRHODJ,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRHODJ")
  CALL MPPDB_CHECK(PRCT,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRCT")
  CALL MPPDB_CHECK(PRRT,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRRT")
  CALL MPPDB_CHECK(PRIT,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRIT")
  CALL MPPDB_CHECK(PRST,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRST")
  CALL MPPDB_CHECK(PRGT,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRGT")
  IF (PRESENT(PSEA))  CALL MPPDB_CHECK(PSEA,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PSEA")
  IF (PRESENT(PTOWN)) CALL MPPDB_CHECK(PTOWN,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PTOWN")
  IF (PRESENT(PRHT))  CALL MPPDB_CHECK(PRHT,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRHT")
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PINPRC,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PINPRC")
  CALL MPPDB_CHECK(PINDEP,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PINDEP")
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRCS")
  CALL MPPDB_CHECK(PRRS,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRRS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRIS")
  CALL MPPDB_CHECK(PRSS,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRSS")
  CALL MPPDB_CHECK(PRGS,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRGS")
  IF (PRESENT(PRHS)) CALL MPPDB_CHECK(PRHS,"RAIN_ICE_SEDIMENTATION_SPLIT beg:PRHS")
END IF

IIU =  size(PRCS, 1 )
IJU =  size(PRCS, 2 )
IKU =  size(PRCS, 3 )
IIJKU = IIU * IJU * IKU

!
#ifndef MNH_OPENACC
ALLOCATE( IC1(IIJKU), IC2(IIJKU), IC3(IIJKU) )
ALLOCATE( IR1(IIJKU), IR2(IIJKU), IR3(IIJKU) )
ALLOCATE( IS1(IIJKU), IS2(IIJKU), IS3(IIJKU) )
ALLOCATE( II1(IIJKU), II2(IIJKU), II3(IIJKU) )
ALLOCATE( IG1(IIJKU), IG2(IIJKU), IG3(IIJKU) )
ALLOCATE( IH1(IIJKU), IH2(IIJKU), IH3(IIJKU) )
#else
!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_SEDIMENTATION_SPLIT' )

CALL MNH_MEM_GET(IC1,IIJKU)
CALL MNH_MEM_GET(IC2,IIJKU)
CALL MNH_MEM_GET(IC3,IIJKU)
CALL MNH_MEM_GET(IR1,IIJKU)
CALL MNH_MEM_GET(IR2,IIJKU)
CALL MNH_MEM_GET(IR3,IIJKU)
CALL MNH_MEM_GET(IS1,IIJKU)
CALL MNH_MEM_GET(IS2,IIJKU)
CALL MNH_MEM_GET(IS3,IIJKU)
CALL MNH_MEM_GET(II1,IIJKU)
CALL MNH_MEM_GET(II2,IIJKU)
CALL MNH_MEM_GET(II3,IIJKU)
CALL MNH_MEM_GET(IG1,IIJKU)
CALL MNH_MEM_GET(IG2,IIJKU)
CALL MNH_MEM_GET(IG3,IIJKU)
CALL MNH_MEM_GET(IH1,IIJKU)
CALL MNH_MEM_GET(IH2,IIJKU)
CALL MNH_MEM_GET(IH3,IIJKU)
#endif
#ifndef MNH_OPENACC
ALLOCATE( GDEP(IIU,IJU) )
ALLOCATE( GSEDIMR(IIU,IJU,IKU) )
ALLOCATE( GSEDIMC(IIU,IJU,IKU) )
ALLOCATE( GSEDIMI(IIU,IJU,IKU) )
ALLOCATE( GSEDIMS(IIU,IJU,IKU) )
ALLOCATE( GSEDIMG(IIU,IJU,IKU) )
ALLOCATE( GSEDIMH(IIU,IJU,IKU) )
#else
CALL MNH_MEM_GET( GDEP,IIU,IJU )
CALL MNH_MEM_GET( GSEDIMR,IIU,IJU,IKU )
CALL MNH_MEM_GET( GSEDIMC,IIU,IJU,IKU )
CALL MNH_MEM_GET( GSEDIMI,IIU,IJU,IKU )
CALL MNH_MEM_GET( GSEDIMS,IIU,IJU,IKU )
CALL MNH_MEM_GET( GSEDIMG,IIU,IJU,IKU )
CALL MNH_MEM_GET( GSEDIMH,IIU,IJU,IKU )
#endif
#ifndef MNH_OPENACC
ALLOCATE( ZCONC_TMP(IIU,IJU) )
ALLOCATE( ZOMPSEA  (IIU,IJU) )
ALLOCATE( ZTMP1_2D (IIU,IJU) )
ALLOCATE( ZTMP2_2D (IIU,IJU) )
ALLOCATE( ZTMP3_2D (IIU,IJU) )
ALLOCATE( ZTMP4_2D (IIU,IJU) )
ALLOCATE( ZCONC3D(IIU,IJU,IKU) )
ALLOCATE( ZRAY   (IIU,IJU,IKU) )
ALLOCATE( ZLBC   (IIU,IJU,IKU) )
ALLOCATE( ZFSEDC (IIU,IJU,IKU) )
ALLOCATE( ZPRCS  (IIU,IJU,IKU) )
ALLOCATE( ZPRRS  (IIU,IJU,IKU) )
ALLOCATE( ZPRSS  (IIU,IJU,IKU) )
ALLOCATE( ZPRGS  (IIU,IJU,IKU) )
ALLOCATE( ZPRHS  (IIU,IJU,IKU) )
ALLOCATE( ZW     (IIU,IJU,IKU) )
ALLOCATE( ZWSED  (IIU,IJU,0:IKU+1) )
ALLOCATE( ZRTMIN(SIZE(XRTMIN)) )
#else
CALL MNH_MEM_GET( ZCONC_TMP,IIU,IJU )
CALL MNH_MEM_GET( ZOMPSEA  ,IIU,IJU )
CALL MNH_MEM_GET( ZTMP1_2D ,IIU,IJU )
CALL MNH_MEM_GET( ZTMP2_2D ,IIU,IJU )
CALL MNH_MEM_GET( ZTMP3_2D ,IIU,IJU )
CALL MNH_MEM_GET( ZTMP4_2D ,IIU,IJU )
CALL MNH_MEM_GET( ZCONC3D,IIU,IJU,IKU )
CALL MNH_MEM_GET( ZRAY   ,IIU,IJU,IKU )
CALL MNH_MEM_GET( ZLBC   ,IIU,IJU,IKU )
CALL MNH_MEM_GET( ZFSEDC ,IIU,IJU,IKU )
CALL MNH_MEM_GET( ZPRCS  ,IIU,IJU,IKU )
CALL MNH_MEM_GET( ZPRRS  ,IIU,IJU,IKU )
CALL MNH_MEM_GET( ZPRSS  ,IIU,IJU,IKU )
CALL MNH_MEM_GET( ZPRGS  ,IIU,IJU,IKU )
CALL MNH_MEM_GET( ZPRHS  ,IIU,IJU,IKU )
CALL MNH_MEM_GET( ZW     ,IIU,IJU,IKU )
CALL MNH_MEM_GET( ZWSED, 1, IIU, 1, IJU, 0, IKU+1 )
CALL MNH_MEM_GET( ZRTMIN , SIZE(XRTMIN) )
#endif

!$acc data present( IC1, IC2, IC3, IR1, IR2, IR3, IS1, IS2, IS3, II1, II2, II3, IG1, IG2, IG3, IH1, IH2, IH3,&
!$acc &            GDEP, GSEDIMR, GSEDIMC,  GSEDIMI,  GSEDIMS,  GSEDIMG,  GSEDIMH )                          &
!$acc &    present( ZRTMIN )                                                                                 &
!$acc &    present(ZCONC_TMP,                                                                                &
!$acc &            ZOMPSEA, ZTMP1_2D, ZTMP2_2D, ZTMP3_2D, ZTMP4_2D, ZCONC3D,                                 &
!$acc &            ZRAY, ZLBC, ZFSEDC,                                                                       &
!$acc &            ZPRCS, ZPRRS, ZPRSS, ZPRGS, ZPRHS, ZW, ZWSED                                              )

IF ( PRESENT( PFPR ) ) THEN
  GPRESENT_PFPR = .TRUE.
ELSE
  GPRESENT_PFPR = .FALSE.
END IF
IF ( PRESENT( PSEA ) ) THEN
  GPRESENT_PSEA = .TRUE.
ELSE
  GPRESENT_PSEA = .FALSE.
END IF
IF ( KRR == 7 ) THEN
   GKRR_7 = .TRUE.
ELSE
   GKRR_7 = .FALSE.
END IF
GSEDIC = OSEDIC

if ( lbudget_rc .and. osedic ) call Budget_store_init( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rr )              call Budget_store_init( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_ri )              call Budget_store_init( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rs )              call Budget_store_init( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rg )              call Budget_store_init( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rh )              call Budget_store_init( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )
!
!        O. Initialization of for sedimentation
!
!$acc kernels present_cr(ZOMPSEA,ZTMP1_2D,zconc_tmp,ztmp3_2d,ztmp2_2d,ztmp4_2d,ZLBC,ZFSEDC) &
!$acc & present_cr(zconc3d,zray,zprrs,zprss)
ZINVTSTEP=1./PTSTEP
ZTSPLITR= PTSTEP / REAL(KSPLITR)
!
IF ( OSEDIC )  PINPRC (:,:) = 0.
IF ( ODEPOSC ) PINDEP (:,:) = 0.
PINPRR (:,:)     = 0.
PINPRR3D (:,:,:) = 0.
PINPRS (:,:)     = 0.
PINPRG (:,:)     = 0.
IF ( KRR == 7 ) PINPRH (:,:) = 0.
IF ( GPRESENT_PFPR ) PFPR(:,:,:,:) = 0.
!
!*       1. Parameters for cloud sedimentation
!
IF ( OSEDIC ) THEN
  ZTMP1 = 0.5 * GAMMA( XNUC  + 1.0 / XALPHAC  ) / ( GAMMA( XNUC  ) )
  ZTMP2 = 0.5 * GAMMA( XNUC2 + 1.0 / XALPHAC2 ) / ( GAMMA( XNUC2 ) )

  IF ( GPRESENT_PSEA ) THEN
     !$mnh_do_concurrent( JI=1:IIU , JJ=1:IJU )
        ZOMPSEA  (JI,JJ) = 1.-PSEA(JI,JJ)
        ZTMP1_2D (JI,JJ) = PSEA(JI,JJ)*XLBC(2)  +ZOMPSEA(JI,JJ)*XLBC(1)
        ZTMP2_2D (JI,JJ) = PSEA(JI,JJ)*XFSEDC(2)+ZOMPSEA(JI,JJ)*XFSEDC(1)
        ZCONC_TMP(JI,JJ) = PSEA(JI,JJ)*XCONC_SEA+ZOMPSEA(JI,JJ)*XCONC_LAND
        ZTMP3_2D (JI,JJ) = (1.-PTOWN(JI,JJ))*ZCONC_TMP(JI,JJ)+PTOWN(JI,JJ)*XCONC_URBAN
        ZTMP4_2D (JI,JJ) = MAX( 1. , ZOMPSEA(JI,JJ)*ZTMP1 + PSEA(JI,JJ)*ZTMP2 )
      !$mnh_end_do()
      !$mnh_do_concurrent( JI=1:IIU ,  JJ=1:IJU , JK=KKTB:KKTE )  
          ZLBC   (JI,JJ,JK) = ZTMP1_2D(JI,JJ)
          ZFSEDC (JI,JJ,JK) = ZTMP2_2D(JI,JJ)
          ZCONC3D(JI,JJ,JK) = ZTMP3_2D(JI,JJ)
          ZRAY   (JI,JJ,JK) = ZTMP4_2D(JI,JJ)
      !$mnh_end_do()
  ELSE
    ZLBC   (:,:,:) = XLBC(1)
    ZFSEDC (:,:,:) = XFSEDC(1)
    ZCONC3D(:,:,:) = XCONC_LAND
    ZTMP3 = MAX(1.,ZTMP1)
    ZRAY   (:,:,:) = ZTMP3
  END IF
END IF
!
!*       2.    compute the fluxes
!
!  optimization by looking for locations where
!  the precipitating fields are larger than a minimal value only !!!
!  For optimization we consider each variable separately

ZRTMIN(:) = XRTMIN(:) * ZINVTSTEP
IF ( OSEDIC ) GSEDIMC(:,:,:) = .FALSE.
GSEDIMR(:,:,:) = .FALSE.
GSEDIMI(:,:,:) = .FALSE.
GSEDIMS(:,:,:) = .FALSE.
GSEDIMG(:,:,:) = .FALSE.
IF ( KRR == 7 ) GSEDIMH(:,:,:) = .FALSE.
!
! ZPiS = Specie i source creating during the current time step
! PRiS = Source of the previous time step
!
IF ( OSEDIC ) THEN
  ZPRCS(:,:,:) = 0.0
  ZPRCS(:,:,:) = PRCS(:,:,:) - PRCT(:,:,:) * ZINVTSTEP
  PRCS (:,:,:) = PRCT(:,:,:)* ZINVTSTEP
END IF
#if 0
ZPRRS(:,:,:) = PRRS(:,:,:)-PRRT(:,:,:)* ZINVTSTEP
ZPRSS(:,:,:) = PRSS(:,:,:)-PRST(:,:,:)* ZINVTSTEP
ZPRGS(:,:,:) = PRGS(:,:,:)-PRGT(:,:,:)* ZINVTSTEP
IF ( KRR == 7 ) ZPRHS(:,:,:) = PRHS(:,:,:)-PRHT(:,:,:)* ZINVTSTEP
PRRS(:,:,:)  = PRRT(:,:,:)* ZINVTSTEP
PRSS(:,:,:)  = PRST(:,:,:)* ZINVTSTEP
PRGS(:,:,:)  = PRGT(:,:,:)* ZINVTSTEP
IF ( KRR == 7 ) PRHS(:,:,:)  = PRHT(:,:,:)* ZINVTSTEP
#else
!$mnh_do_concurrent( JI=1:IIU , JJ=1:IJU , JK=1:IKU )
      ZPRRS(JI,JJ,JK) = PRRS(JI,JJ,JK) - PRRT(JI,JJ,JK) * ZINVTSTEP
      ZPRSS(JI,JJ,JK) = PRSS(JI,JJ,JK) - PRST(JI,JJ,JK) * ZINVTSTEP
      ZPRGS(JI,JJ,JK) = PRGS(JI,JJ,JK) - PRGT(JI,JJ,JK) * ZINVTSTEP
      PRRS (JI,JJ,JK) =                  PRRT(JI,JJ,JK) * ZINVTSTEP
      PRSS (JI,JJ,JK) =                  PRST(JI,JJ,JK) * ZINVTSTEP
      PRGS (JI,JJ,JK) =                  PRGT(JI,JJ,JK) * ZINVTSTEP
!$mnh_end_do()
IF ( KRR == 7 ) THEN
!$mnh_do_concurrent( JI=1:IIU , JJ=1:IJU , JK=1:IKU )   
        ZPRHS(JI,JJ,JK) = PRHS(JI,JJ,JK) - PRHT(JI,JJ,JK) * ZINVTSTEP
        PRHS (JI,JJ,JK) =                  PRHT(JI,JJ,JK) * ZINVTSTEP
!$mnh_end_do()
END IF
#endif
!$acc end kernels
!
! PRiS = Source of the previous time step + source created during the subtime
! step
!
DO JN = 1 , KSPLITR
!$acc kernels present_cr(gsedimc,gsedimr,gsedimi,gsedims,gsedimg,gsedimh)
  IF( JN == 1 ) THEN
    IF ( OSEDIC ) PRCS(:,:,:) = PRCS(:,:,:) + ZPRCS(:,:,:) / KSPLITR
    PRRS(:,:,:) = PRRS(:,:,:) + ZPRRS(:,:,:) / KSPLITR
    PRSS(:,:,:) = PRSS(:,:,:) + ZPRSS(:,:,:) / KSPLITR
    PRGS(:,:,:) = PRGS(:,:,:) + ZPRGS(:,:,:) / KSPLITR
    IF ( KRR == 7 ) PRHS(:,:,:) = PRHS(:,:,:) + ZPRHS(:,:,:) / KSPLITR
    DO JK = KKTB , KKTE
      ZW(:,:,JK) = ZTSPLITR / ( PRHODREF(:,:,JK) * PDZZ(:,:,JK) )
    END DO
  ELSE
    IF ( OSEDIC ) PRCS(:,:,:) = PRCS(:,:,:) + ZPRCS(:,:,:) * ZTSPLITR
    PRRS(:,:,:) = PRRS(:,:,:) + ZPRRS(:,:,:) * ZTSPLITR
    PRSS(:,:,:) = PRSS(:,:,:) + ZPRSS(:,:,:) * ZTSPLITR
    PRGS(:,:,:) = PRGS(:,:,:) + ZPRGS(:,:,:) * ZTSPLITR
    IF ( KRR == 7 ) PRHS(:,:,:) = PRHS(:,:,:) + ZPRHS(:,:,:) * ZTSPLITR
  END IF
  !
#ifndef MNH_COMPILER_CCE  
  !$mnh_expand_array( JI=KIB:KIE,JJ=KJB:KJE,JK=KKTB:KKTE )
#endif
  IF ( GSEDIC ) GSEDIMC(:,:,:) =                &
                   PRCS(:,:,:) > ZRTMIN(2)
  GSEDIMR(:,:,:) =                            &
                   PRRS(:,:,:) > ZRTMIN(3)
  GSEDIMI(:,:,:) =                            &
                   PRIS(:,:,:) > ZRTMIN(4)
  GSEDIMS(:,:,:) =                            &
                   PRSS(:,:,:) > ZRTMIN(5)
  GSEDIMG(:,:,:) =                            &
                   PRGS(:,:,:) > ZRTMIN(6)
  IF ( GKRR_7 ) GSEDIMH(:,:,:) =            &
                   PRHS(:,:,:) > ZRTMIN(7)
#ifndef MNH_COMPILER_CCE  
  !$mnh_end_expand_array() ! CONCURRENT
#endif 
  !$acc end kernels
!
#ifndef MNH_OPENACC
  IF ( OSEDIC ) ISEDIMC = COUNTJV( GSEDIMC(:,:,:), IC1(:), IC2(:), IC3(:) )
  ISEDIMR = COUNTJV( GSEDIMR(:,:,:), IR1(:), IR2(:), IR3(:) )
  ISEDIMI = COUNTJV( GSEDIMI(:,:,:), II1(:), II2(:), II3(:) )
  ISEDIMS = COUNTJV( GSEDIMS(:,:,:), IS1(:), IS2(:), IS3(:) )
  ISEDIMG = COUNTJV( GSEDIMG(:,:,:), IG1(:), IG2(:), IG3(:) )
  IF ( KRR == 7 ) ISEDIMH = COUNTJV( GSEDIMH(:,:,:), IH1(:), IH2(:), IH3(:) )
#else
  IF ( OSEDIC ) CALL COUNTJV_DEVICE( GSEDIMC, IC1, IC2, IC3, ISEDIMC )
  CALL COUNTJV_DEVICE( GSEDIMR, IR1, IR2, IR3, ISEDIMR )
  CALL COUNTJV_DEVICE( GSEDIMI, II1, II2, II3, ISEDIMI )
  CALL COUNTJV_DEVICE( GSEDIMS, IS1, IS2, IS3, ISEDIMS )
  CALL COUNTJV_DEVICE( GSEDIMG, IG1, IG2, IG3, ISEDIMG )
  IF ( KRR == 7 ) CALL COUNTJV_DEVICE( GSEDIMH, IH1, IH2, IH3, ISEDIMH )
#endif
!
!*       2.1   for cloud
!
!$acc kernels
  IF ( OSEDIC ) THEN
    ZWSED(:,:,:) = 0.
    IF( JN==1 ) PRCS(:,:,:) = PRCS(:,:,:) * PTSTEP
     !$mnh_do_concurrent (JL=1:ISEDIMC)
       ZRSLOC = PRCS(IC1(JL),IC2(JL),IC3(JL))
       ZRTLOC = PRCT(IC1(JL),IC2(JL),IC3(JL))
       IF (ZRSLOC > ZRTMIN(2) .AND. ZRTLOC > XRTMIN(2)) THEN
         ZRHODREFLOC =  PRHODREF(IC1(JL),IC2(JL),IC3(JL))
         ZFSEDC1D = ZFSEDC(IC1(JL),IC2(JL),IC3(JL))
         ZWLBDC = ZLBC(IC1(JL),IC2(JL),IC3(JL))
         ZCONC = ZCONC3D(IC1(JL),IC2(JL),IC3(JL))
         ZRAY1D = ZRAY(IC1(JL),IC2(JL),IC3(JL))
         ZZT = PTHT(IC1(JL),IC2(JL),IC3(JL))
         ZPRES = PPABST(IC1(JL),IC2(JL),IC3(JL))
!Problems with PGI (18.10). OK if 2 lines are merged!
!          ZWLBDC = ZWLBDC * ZCONC / (ZRHODREFLOC * ZRTLOC)
! #ifndef MNH_BITREP
!          ZWLBDC = ZWLBDC**XLBEXC
! #else
!          ZWLBDC = BR_POW(ZWLBDC,XLBEXC)
! #endif
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
        ZWLBDC = (ZWLBDC * ZCONC / (ZRHODREFLOC * ZRTLOC))**XLBEXC
#else
        ZWLBDC = BR_POW(ZWLBDC * ZCONC / (ZRHODREFLOC * ZRTLOC),XLBEXC)
#endif
        ZRAY1D = ZRAY1D / ZWLBDC !! ZRAY : mean diameter=M(1)/2
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
        ZZT    = ZZT * (ZPRES/XP00)**(XRD/XCPD)
#else
        ZZT    = ZZT * BR_POW(ZPRES/XP00,XRD/XCPD)
#endif
!!$        ZWLBDA = 6.6E-8*(101325./ZPRES)*(ZZT/293.15)
        ZWLBDA = 2281.238e-8*(ZZT/ZPRES)        
        ZCC    = XCC*(1.+1.26*ZWLBDA/ZRAY1D) !! XCC modified for cloud
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
        ZWSED (IC1(JL),IC2(JL),IC3(JL)) = ZRHODREFLOC**(-XCEXVT +1.0 )           &
                                          * ZWLBDC**(-XDC)*ZCC*ZFSEDC1D * ZRSLOC
#else
        ZWSED (IC1(JL),IC2(JL),IC3(JL)) = BR_POW(ZRHODREFLOC,-XCEXVT +1.0 )           &
                                          * BR_POW(ZWLBDC,-XDC)*ZCC*ZFSEDC1D * ZRSLOC
#endif
      END IF
    !$mnh_end_do()
    DO JK = KKTB , KKTE
      PRCS(:,:,JK) = PRCS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
    END DO
    IF ( GPRESENT_PFPR ) THEN
      DO JK = KKTB , KKTE
        PFPR(:,:,JK,2)=ZWSED(:,:,JK)
      ENDDO
    ENDIF
    PINPRC(:,:) = PINPRC(:,:) + ZWSED(:,:,KKB) / XRHOLW / KSPLITR
    IF( JN == KSPLITR ) THEN
      PRCS(:,:,:) = PRCS(:,:,:) * ZINVTSTEP
    END IF
  END IF
!
!*       2.2   for rain
!
  IF( JN==1 ) PRRS(:,:,:) = PRRS(:,:,:) * PTSTEP
  ZWSED(:,:,:) = 0.
  !$mnh_do_concurrent (JL=1:ISEDIMR)
    ZRSLOC = PRRS(IR1(JL),IR2(JL),IR3(JL))
    IF( ZRSLOC > ZRTMIN(3) ) THEN
      ZRHODREFLOC =  PRHODREF(IR1(JL),IR2(JL),IR3(JL))

#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
      ZWSED (IR1(JL),IR2(JL),IR3(JL))= XFSEDR  * ZRSLOC**XEXSEDR *   &
                                       ZRHODREFLOC**(XEXSEDR-XCEXVT)
#else
      ZWSED (IR1(JL),IR2(JL),IR3(JL))= XFSEDR  * BR_POW(ZRSLOC,XEXSEDR) *   &
                                       BR_POW(ZRHODREFLOC,XEXSEDR-XCEXVT)
#endif
    END IF
  !$mnh_end_do()
  DO JK = KKTB , KKTE
    PRRS(:,:,JK) = PRRS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
  END DO
  IF ( GPRESENT_PFPR ) THEN
    DO JK = KKTB , KKTE
      PFPR(:,:,JK,3)=ZWSED(:,:,JK)
    ENDDO
  ENDIF
  PINPRR(:,:) = PINPRR(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
  PINPRR3D(:,:,:) = PINPRR3D(:,:,:) + ZWSED(:,:,1:KKT)/XRHOLW/KSPLITR
  IF( JN == KSPLITR ) THEN
    PRRS(:,:,:) = PRRS(:,:,:) * ZINVTSTEP
  END IF
!
!*       2.3   for pristine ice
!
  IF( JN==1 ) PRIS(:,:,:) = PRIS(:,:,:) * PTSTEP
  ZWSED(:,:,:) = 0.
  !$mnh_do_concurrent (JL=1:ISEDIMI)
    ZRSLOC = PRIS(II1(JL),II2(JL),II3(JL))
    IF( ZRSLOC >  MAX(ZRTMIN(4),1.0E-7 )) THEN ! limitation of the McF&H formula
      ZRHODREFLOC =  PRHODREF(II1(JL),II2(JL),II3(JL))
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
      ZWSED (II1(JL),II2(JL),II3(JL)) = XFSEDI * ZRSLOC *                  &
                                        ZRHODREFLOC**(1.0-XCEXVT) *        & !    McF&H
                                        MAX( 0.05E6,-0.15319E6-0.021454E6* &
                                        ALOG(ZRHODREFLOC*ZRSLOC) )**XEXCSEDI
#else
      ZWSED (II1(JL),II2(JL),II3(JL)) = XFSEDI * ZRSLOC *                          &
                                        BR_POW(ZRHODREFLOC,1.0-XCEXVT) *           & !    McF&H
                                        BR_POW( MAX( 0.05E6,-0.15319E6-0.021454E6* &
                                        BR_LOG(ZRHODREFLOC*ZRSLOC) ), XEXCSEDI)
#endif
    END IF
  !$mnh_end_do()
  DO JK = KKTB , KKTE
    PRIS(:,:,JK) = PRIS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
  END DO
  IF ( GPRESENT_PFPR ) THEN
    DO JK = KKTB , KKTE
      PFPR(:,:,JK,4)=ZWSED(:,:,JK)
    ENDDO
  ENDIF
  IF( JN == KSPLITR ) THEN
    PRIS(:,:,:) = PRIS(:,:,:) * ZINVTSTEP
  END IF
!
!*       2.4   for aggregates/snow
!
  IF( JN==1 ) PRSS(:,:,:) = PRSS(:,:,:) * PTSTEP
  ZWSED(:,:,:) = 0.
  !$mnh_do_concurrent (JL=1:ISEDIMS)
    ZRSLOC = PRSS(IS1(JL),IS2(JL),IS3(JL))
    IF( ZRSLOC > ZRTMIN(5) ) THEN
      ZRHODREFLOC =  PRHODREF(IS1(JL),IS2(JL),IS3(JL))
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
      ZWSED (IS1(JL),IS2(JL),IS3(JL)) = XFSEDS * ZRSLOC**XEXSEDS *  &
                                        ZRHODREFLOC**(XEXSEDS-XCEXVT)
#else
      ZWSED (IS1(JL),IS2(JL),IS3(JL)) = XFSEDS * BR_POW(ZRSLOC,XEXSEDS) *  &
                                        BR_POW(ZRHODREFLOC,XEXSEDS-XCEXVT)
#endif
    END IF
  !$mnh_end_do()
  DO JK = KKTB , KKTE
    PRSS(:,:,JK) = PRSS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
  END DO
  IF ( GPRESENT_PFPR ) THEN
    DO JK = KKTB , KKTE
      PFPR(:,:,JK,5)=ZWSED(:,:,JK)
    ENDDO
  ENDIF
  PINPRS(:,:) = PINPRS(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
  IF( JN == KSPLITR ) THEN
    PRSS(:,:,:) = PRSS(:,:,:) * ZINVTSTEP
  END IF
!
!*       2.5   for graupeln
!
  ZWSED(:,:,:) = 0.
  IF( JN==1 ) PRGS(:,:,:) = PRGS(:,:,:) * PTSTEP
  !$mnh_do_concurrent (JL=1:ISEDIMG)
    ZRSLOC = PRGS(IG1(JL),IG2(JL),IG3(JL))
    IF( ZRSLOC > ZRTMIN(6) ) THEN
      ZRHODREFLOC =  PRHODREF(IG1(JL),IG2(JL),IG3(JL))
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
      ZWSED (IG1(JL),IG2(JL),IG3(JL)) = XFSEDG  * ZRSLOC**XEXSEDG *   &
                                        ZRHODREFLOC**(XEXSEDG-XCEXVT)
#else
      ZWSED (IG1(JL),IG2(JL),IG3(JL)) = XFSEDG  * BR_POW(ZRSLOC,XEXSEDG) *   &
                                        BR_POW(ZRHODREFLOC,XEXSEDG-XCEXVT)
#endif
    END IF
  !$mnh_end_do()
  DO JK = KKTB , KKTE
    PRGS(:,:,JK) = PRGS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
  END DO
  IF ( GPRESENT_PFPR ) THEN
    DO JK = KKTB , KKTE
      PFPR(:,:,JK,6)=ZWSED(:,:,JK)
    ENDDO
  ENDIF
  PINPRG(:,:) = PINPRG(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
  IF( JN == KSPLITR ) THEN
    PRGS(:,:,:) = PRGS(:,:,:) * ZINVTSTEP
  END IF
!
!*       2.6   for hail
!
  IF ( KRR == 7 ) THEN
    IF( JN==1 ) PRHS(:,:,:) = PRHS(:,:,:) * PTSTEP
    ZWSED(:,:,:) = 0.
    !$mnh_do_concurrent (JL=1:ISEDIMH)
      ZRSLOC = PRHS(IH1(JL),IH2(JL),IH3(JL))
      IF( ZRSLOC > ZRTMIN(7) ) THEN
        ZRHODREFLOC =  PRHODREF(IH1(JL),IH2(JL),IH3(JL))
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
        ZWSED (IH1(JL),IH2(JL),IH3(JL)) = XFSEDH  * ZRSLOC**XEXSEDH *   &
                                         ZRHODREFLOC**(XEXSEDH-XCEXVT)
#else
        ZWSED (IH1(JL),IH2(JL),IH3(JL)) = XFSEDH  * BR_POW(ZRSLOC,XEXSEDH) *   &
                                          BR_POW(ZRHODREFLOC,XEXSEDH-XCEXVT)
#endif
      END IF
    !$mnh_end_do()
    DO JK = KKTB , KKTE
      PRHS(:,:,JK) = PRHS(:,:,JK) + ZW(:,:,JK)*(ZWSED(:,:,JK+KKL)-ZWSED(:,:,JK))
    END DO
    IF ( GPRESENT_PFPR ) THEN
      DO JK = KKTB , KKTE
        PFPR(:,:,JK,7)=ZWSED(:,:,JK)
      ENDDO
    ENDIF
    PINPRH(:,:) = PINPRH(:,:) + ZWSED(:,:,KKB)/XRHOLW/KSPLITR
    IF( JN == KSPLITR ) THEN
      PRHS(:,:,:) = PRHS(:,:,:) * ZINVTSTEP
    END IF
  END IF
!$acc end kernels
!
END DO
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

!$acc kernels
  GDEP(:,:) = .FALSE.
  GDEP(KIB:KIE,KJB:KJE) =    PRCS(KIB:KIE,KJB:KJE,KKB) >0
  WHERE (GDEP)
    PRCS(:,:,KKB) = PRCS(:,:,KKB) - XVDEPOSC * PRCT(:,:,KKB) / PDZZ(:,:,KKB)
    PINPRC(:,:) = PINPRC(:,:) + XVDEPOSC * PRCT(:,:,KKB) * PRHODREF(:,:,KKB) /XRHOLW
    PINDEP(:,:) = XVDEPOSC * PRCT(:,:,KKB) * PRHODREF(:,:,KKB) /XRHOLW
  END WHERE
!$acc end kernels

  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'DEPO', prcs(:, :, :) * prhodj(:, :, :) )
END IF

IF (MPPDB_INITIALIZED) THEN
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PINPRC,"RAIN_ICE_SEDIMENTATION_SPLIT end:PINPRC")
  CALL MPPDB_CHECK(PINDEP,"RAIN_ICE_SEDIMENTATION_SPLIT end:PINDEP")
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_SEDIMENTATION_SPLIT end:PRCS")
  CALL MPPDB_CHECK(PRRS,"RAIN_ICE_SEDIMENTATION_SPLIT end:PRRS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_SEDIMENTATION_SPLIT end:PRIS")
  CALL MPPDB_CHECK(PRSS,"RAIN_ICE_SEDIMENTATION_SPLIT end:PRSS")
  CALL MPPDB_CHECK(PRGS,"RAIN_ICE_SEDIMENTATION_SPLIT end:PRGS")
  IF (PRESENT(PRHS)) CALL MPPDB_CHECK(PRHS,"RAIN_ICE_SEDIMENTATION_SPLIT end:PRHS")
  !Check all OUT arrays
  CALL MPPDB_CHECK(PINPRR,"RAIN_ICE_SEDIMENTATION_SPLIT end:PINPRR")
  CALL MPPDB_CHECK(PINPRS,"RAIN_ICE_SEDIMENTATION_SPLIT end:PINPRS")
  CALL MPPDB_CHECK(PINPRG,"RAIN_ICE_SEDIMENTATION_SPLIT end:PINPRG")
  CALL MPPDB_CHECK(PINPRR3D,"RAIN_ICE_SEDIMENTATION_SPLIT end:PINPRR3D")
  IF (PRESENT(PINPRH)) CALL MPPDB_CHECK(PRHS,"RAIN_ICE_SEDIMENTATION_SPLIT end:PINPRH")
  IF (PRESENT(PFPR))   CALL MPPDB_CHECK(PFPR,"RAIN_ICE_SEDIMENTATION_SPLIT end:PFPR")
END IF

!$acc end data

#ifndef MNH_OPENACC
DEALLOCATE(IC1, IC2, IC3, IR1, IR2, IR3, IS1, IS2, IS3, II1, II2, II3, IG1, IG2, IG3, IH1, IH2, IH3)
DEALLOCATE(GDEP, GSEDIMR, GSEDIMC,  GSEDIMI,  GSEDIMS,  GSEDIMG,  GSEDIMH)
DEALLOCATE(ZRTMIN)
DEALLOCATE(ZCONC_TMP,ZOMPSEA, ZTMP1_2D, ZTMP2_2D, ZTMP3_2D, ZTMP4_2D, ZCONC3D)
DEALLOCATE(ZRAY, ZLBC, ZFSEDC,ZPRCS, ZPRRS, ZPRSS, ZPRGS, ZPRHS, ZW)
DEALLOCATE(ZWSED)
#else
!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'RAIN_ICE_SEDIMENTATION_SPLIT' )
#endif

!$acc end data

END SUBROUTINE RAIN_ICE_SEDIMENTATION_SPLIT

END MODULE MODE_RAIN_ICE_SEDIMENTATION_SPLIT
