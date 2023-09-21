!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_SEDIMENTATION
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_SEDIMENTATION(D, CST, ICEP, ICED, PARAMI, BUCONF, &
                             &PTSTEP, KRR, PDZZ, &
                             &PLVFACT, PLSFACT, PRHODREF, PPABST, PTHT, PT, PRHODJ, &
                             &PTHS, PRVS, PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                             &PINPRC, PINPRR, PINPRS, PINPRG, &
                             &TBUDGETS, KBUDGETS, &
                             &PSEA, PTOWN,  &
                             &PINPRH, PRHT, PRHS, PFPR)
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
USE MODD_BUDGET,         ONLY: TBUDGETDATA, TBUDGETCONF_t, NBUDGET_RC, &
                               NBUDGET_RI, NBUDGET_RR, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH
USE MODD_CST, ONLY: CST_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
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
TYPE(TBUDGETCONF_t),                         INTENT(IN)    :: BUCONF
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
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(D%NIJT,D%NKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(D%NIJT,D%NKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
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
  CALL ICE4_SEDIMENTATION_SPLIT(D, CST, ICEP, ICED, PARAMI, &
                               &PTSTEP, KRR, PDZZ, &
                               &PRHODREF, PPABST, PTHT, PT, PRHODJ, &
                               &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                               &PINPRC, PINPRR, ZINPRI, PINPRS, PINPRG, &
                               &PSEA=PSEA, PTOWN=PTOWN, &
                               &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR)
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
  CALL ICE4_CORRECT_NEGATIVITIES(D, ICED, KRR, PRVS, PRCS, PRRS, &
                                &PRIS, PRSS, PRGS, &
                                &PTHS, PLVFACT, PLSFACT, PRHS)
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
IF (LHOOK) CALL DR_HOOK('ICE4_SEDIMENTATION', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_SEDIMENTATION
END MODULE MODE_ICE4_SEDIMENTATION
