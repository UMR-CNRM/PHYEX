!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_SEDIMENTATION
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_SEDIMENTATION(D, CST, ICEP, ICED, PARAMI, ELECP, ELECD, BUCONF, &
                             &OELEC, OSEDIM_BEARD, HCLOUD, PTSTEP, KRR, PDZZ, PTHVREFZIKB, &
                             &PLVFACT, PLSFACT, PRHODREF, PPABST, PTHT, PT, PRHODJ, &
                             &PTHS, PRVS, PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                             &PINPRC, PINPRR, PINPRS, PINPRG, &
                             &PQCT, PQRT, PQIT, PQST, PQGT, PQCS, PQRS, PQIS, PQSS, PQGS, PEFIELDW, &
                             &TBUDGETS, KBUDGETS, &
                             &PSEA, PTOWN,  &
                             &PINPRH, PRHT, PRHS, PFPR, &
                             &PQHT, PQHS)
!!
!!**  PURPOSE
!!    -------
!!      Computes the sedimentation
!!
!!    AUTHOR
!!    ------
!!      S. Riette code extracted from rain_ice
!!
!!    MODIFICATIONS
!!    -------------
!!
!
!
!*      0. DECLARATIONS
!          ------------
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_BUDGET,         ONLY: TBUDGETDATA, TBUDGETCONF_t, NBUDGET_RC, NBUDGET_SV1, &
                               NBUDGET_RI, NBUDGET_RR, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH
USE MODD_CST, ONLY: CST_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_ELEC_DESCR,     ONLY: ELEC_DESCR_t
USE MODD_ELEC_PARAM,     ONLY: ELEC_PARAM_t
USE MODD_NSV,            ONLY: NSV_ELECBEG
!
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
USE MODE_BUDGET_PHY,         ONLY: BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY
!
USE MODE_ICE4_SEDIMENTATION_STAT, ONLY: ICE4_SEDIMENTATION_STAT                                                                     
USE MODE_ICE4_SEDIMENTATION_SPLIT, ONLY: ICE4_SEDIMENTATION_SPLIT                                                                   
USE MODE_ICE4_CORRECT_NEGATIVITIES, ONLY: ICE4_CORRECT_NEGATIVITIES
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t),                            INTENT(IN)    :: D       !array dimensions
TYPE(CST_t),                                 INTENT(IN)    :: CST
TYPE(RAIN_ICE_PARAM_t),                      INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),                      INTENT(IN)    :: ICED
TYPE(PARAM_ICE_t),                           INTENT(IN)    :: PARAMI
TYPE(ELEC_PARAM_t),                          INTENT(IN)    :: ELECP   ! electrical parameters
TYPE(ELEC_DESCR_t),                          INTENT(IN)    :: ELECD   ! electrical descriptive csts
TYPE(TBUDGETCONF_t),                         INTENT(IN)    :: BUCONF
LOGICAL,                                     INTENT(IN)    :: OELEC   ! switch to activate cloud electrification
LOGICAL,                                     INTENT(IN)    :: OSEDIM_BEARD  ! Switch for effect of electrical forces on sedim.
CHARACTER (LEN=4),                           INTENT(IN)    :: HCLOUD  ! Kind of microphysical scheme
REAL,                                        INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                                     INTENT(IN)    :: KRR     ! Number of moist variable
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PLVFACT
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PLSFACT
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PT      ! Temperature at time t
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(INOUT) :: PTHS
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(INOUT) :: PRVS
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(INOUT) :: PRGS    ! Graupel m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),               INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(D%NIJT),                     INTENT(OUT)   :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(D%NIJT),                     INTENT(OUT)   :: PINPRR  ! Rain instant precip
REAL, DIMENSION(D%NIJT),                     INTENT(OUT)   :: PINPRS  ! Snow instant precip
REAL, DIMENSION(D%NIJT),                     INTENT(OUT)   :: PINPRG  ! Graupel instant precip
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS),      INTENT(INOUT) :: TBUDGETS                                                                   
INTEGER,                                     INTENT(IN)    :: KBUDGETS
REAL, INTENT(IN)                :: PTHVREFZIKB ! Reference thv at IKB for electricity
!
! variables for cloud electricity
!
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(IN) :: PQCT   ! Cloud droplet | 
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(IN) :: PQRT   ! Rain          | electric
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(IN) :: PQIT   ! Ice crystals  |  charge 
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(IN) :: PQST   ! Snow          |   at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(IN) :: PQGT   ! Graupel       |
!
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQCS   ! Cloud droplet | 
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQRS   ! Rain          | electric
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQIS   ! Ice crystals  |  charge 
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQSS   ! Snow          |  source
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), INTENT(INOUT) :: PQGS   ! Graupel       |
!
REAL, DIMENSION(MERGE(D%NIJT,0,OSEDIM_BEARD),MERGE(D%NKT,0,OSEDIM_BEARD)), INTENT(IN) :: PEFIELDW ! vert. E field
!
! optional variables
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(D%NIJT,D%NKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(D%NIJT,D%NKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(IN)    :: PQHT ! Hail electric charge at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQHS ! Hail electric charge source
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRCT, ZRRT, ZRIT, ZRST, ZRGT, ZRHT
REAL, DIMENSION(D%NIJT) :: ZINPRI
INTEGER :: JK, JIJ, IKTB, IKTE, IIJB, IIJE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION', 0, ZHOOK_HANDLE)
IKTB=D%NKTB
IKTE=D%NKTE
IIJB=D%NIJB
IIJE=D%NIJE
!
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!               -------------------------------------
!
!
!
IF (BUCONF%LBUDGET_RC .AND. PARAMI%LSEDIC) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'SEDI', PRCS(:, :) * PRHODJ(:, :))
IF (BUCONF%LBUDGET_RR)              CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'SEDI', PRRS(:, :) * PRHODJ(:, :))
IF (BUCONF%LBUDGET_RI)              CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'SEDI', PRIS(:, :) * PRHODJ(:, :))
IF (BUCONF%LBUDGET_RS)              CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'SEDI', PRSS(:, :) * PRHODJ(:, :))
IF (BUCONF%LBUDGET_RG)              CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'SEDI', PRGS(:, :) * PRHODJ(:, :))
IF (BUCONF%LBUDGET_RH .AND. KRR==7) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RH), 'SEDI', PRHS(:, :) * PRHODJ(:, :))
!
! budget of electric charges
IF (BUCONF%LBUDGET_SV .AND. OELEC) THEN
  IF (PARAMI%LSEDIC) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+1), 'SEDI', PQCS(:, :) * PRHODJ(:, :))
  CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+2), 'SEDI', PQRS(:, :) * PRHODJ(:, :))
  CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+3), 'SEDI', PQIS(:, :) * PRHODJ(:, :))
  CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+4), 'SEDI', PQSS(:, :) * PRHODJ(:, :))
  CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+5), 'SEDI', PQGS(:, :) * PRHODJ(:, :))
  IF (KRR == 7)      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+6), 'SEDI', PQHS(:, :) * PRHODJ(:, :))
END IF
!
IF(PARAMI%CSEDIM=='STAT') THEN
  DO JK = IKTB,IKTE
    DO JIJ = IIJB,IIJE
      ZRCT(JIJ,JK)=PRCS(JIJ,JK)*PTSTEP
      ZRRT(JIJ,JK)=PRRS(JIJ,JK)*PTSTEP
      ZRIT(JIJ,JK)=PRIS(JIJ,JK)*PTSTEP
      ZRST(JIJ,JK)=PRSS(JIJ,JK)*PTSTEP
      ZRGT(JIJ,JK)=PRGS(JIJ,JK)*PTSTEP
      IF (KRR==7) ZRHT(JIJ,JK)=PRHS(JIJ,JK)*PTSTEP
    ENDDO
  ENDDO
  CALL ICE4_SEDIMENTATION_STAT(D, CST, ICEP, ICED, PARAMI, &
                              &PTSTEP, KRR, PDZZ, &
                              &PRHODREF, PPABST, PTHT, PT, PRHODJ, &
                              &PRCS, ZRCT, PRRS, ZRRT, PRIS, ZRIT,&
                              &PRSS, ZRST, PRGS, ZRGT,&
                              &PINPRC, PINPRR, ZINPRI, PINPRS, PINPRG, &
                              &PSEA=PSEA, PTOWN=PTOWN, &
                              &PINPRH=PINPRH, PRHT=ZRHT, PRHS=PRHS, PFPR=PFPR)
  PINPRS(:) = PINPRS(:) + ZINPRI(:)
  !No negativity correction here as we apply sedimentation on PR.S*PTSTEP variables
ELSEIF(PARAMI%CSEDIM=='SPLI') THEN
  CALL ICE4_SEDIMENTATION_SPLIT(D, CST, ICEP, ICED, PARAMI, ELECP, ELECD, &
                               &OELEC, OSEDIM_BEARD, PTHVREFZIKB, HCLOUD, PTSTEP, KRR, PDZZ, &
                               &PRHODREF, PPABST, PTHT, PT, PRHODJ, &
                               &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                               &PINPRC, PINPRR, ZINPRI, PINPRS, PINPRG, &
                               &PQCT, PQRT, PQIT, PQST, PQGT, PQCS, PQRS, PQIS, PQSS, PQGS, &
                               &PEFIELDW, &
                               &PSEA=PSEA, PTOWN=PTOWN, &
                               &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR, &
                               &PQHT=PQHT, PQHS=PQHS)
  PINPRS(:) = PINPRS(:) + ZINPRI(:)
  !We correct negativities with conservation
  !SPLI algorith uses a time-splitting. Inside the loop a temporary m.r. is used.
  !   It is initialized with the m.r. at T and is modified by two tendencies:
  !   sedimentation tendency and an external tendency which represents all other
  !   processes (mainly advection and microphysical processes). If both tendencies
  !   are negative, sedimentation can remove a species at a given sub-timestep. From
  !   this point sedimentation stops for the remaining sub-timesteps but the other tendency
  !   will be still active and will lead to negative values.
  !   We could prevent the algorithm to not consume too much a species, instead we apply
  !   a correction here.
!++cb-- il faudrait faire la correction correspondante sur les charges electriques pour eviter de se retrouver
! avec des points ou il y a de la charge mais pas de masse !
  CALL ICE4_CORRECT_NEGATIVITIES(D, ICED, KRR, PRVS, PRCS, PRRS, &
                                &PRIS, PRSS, PRGS, &
                                &PTHS, PLVFACT, PLSFACT, PRHS)
!  CALL ICE4_CORRECT_NEGATIVITIES(D, ICED, KRR, PRVS, PRCS, PRRS, &
!                               &PRIS, PRSS, PRGS, &
!                                &PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS, &
!                                &PTHS, PLVFACT, PLSFACT, PRHS, PQHS)
ELSEIF(PARAMI%CSEDIM=='NONE') THEN
ELSE
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'RAIN_ICE', 'no sedimentation scheme for PARAMI%CSEDIM='//PARAMI%CSEDIM)
END IF
!
!
IF (BUCONF%LBUDGET_RC .AND. PARAMI%LSEDIC) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'SEDI', PRCS(:, :) * PRHODJ(:, :))
IF (BUCONF%LBUDGET_RR)              CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'SEDI', PRRS(:, :) * PRHODJ(:, :))
IF (BUCONF%LBUDGET_RI)              CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'SEDI', PRIS(:, :) * PRHODJ(:, :))
IF (BUCONF%LBUDGET_RS)              CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'SEDI', PRSS(:, :) * PRHODJ(:, :))
IF (BUCONF%LBUDGET_RG)              CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'SEDI', PRGS(:, :) * PRHODJ(:, :))
IF (BUCONF%LBUDGET_RH .AND. KRR==7) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RH), 'SEDI', PRHS(:, :) * PRHODJ(:, :))
!
! Budget for electric charges
IF (BUCONF%LBUDGET_SV .AND. OELEC) THEN
  IF (PARAMI%LSEDIC) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+1), 'SEDI', PQCS(:, :) * PRHODJ(:, :))
  CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+2), 'SEDI', PQRS(:, :) * PRHODJ(:, :))
  CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+3), 'SEDI', PQIS(:, :) * PRHODJ(:, :))
  CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+4), 'SEDI', PQSS(:, :) * PRHODJ(:, :))
  CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+5), 'SEDI', PQGS(:, :) * PRHODJ(:, :))
  IF (KRR == 7)      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1-1+NSV_ELECBEG+6), 'SEDI', PQHS(:, :) * PRHODJ(:, :))
END IF
!
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_SEDIMENTATION
END MODULE MODE_ICE4_SEDIMENTATION
