!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_NUCLEATION_PROCS
  IMPLICIT NONE
CONTAINS
!     ###############################################################################
  SUBROUTINE LIMA_NUCLEATION_PROCS (LIMAP, LIMAW, LIMAC, TNSV, D, CST, NEBN, BUCONF, TBUDGETS, KBUDGETS, &
                                    KCARB, KSOA, KSP, ODUST, OSALT, OORILAM,        &
                                    PTSTEP, PRHODJ,                                 &
                                    PRHODREF, PEXNREF, PPABST, PT, PDTHRAD, PW_NU,  &
                                    PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT, &
                                    PCCT, PCRT, PCIT, PCIT_SHAPE, PAERO,PSOLORG, PMI, HACTCCN,           &
                                    PNFT, PNAT, PIFT, PINT, PNIT, PNHT,             &
                                    PCLDFR, PICEFR, PPRCFR,                         &
                                    PTOT_RV_HENU, PTOT_RC_HINC, PTOT_RI_HIND,       &
                                    PTOT_RV_HONH                                    )
!     ###############################################################################
!
!!    PURPOSE
!!    -------
!!      Compute nucleation processes for the time-split version of LIMA
!!
!!    AUTHOR
!!    ------
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!  M. Leriche     06/2019: missing update of PNFT after CCN hom. ncl.
!  P. Wautelet 27/02/2020: bugfix: PNFT was not updated after LIMA_CCN_HOM_FREEZING
!  P. Wautelet 27/02/2020: add Z_TH_HINC variable (for budgets)
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!  B. Vie      03/03/2020: use DTHRAD instead of dT/dt in Smax diagnostic computation
!  B. Vie         03/2022: Add option for 1-moment pristine ice
!  C. Barthe      06/2022: add dummy arguments (mass transfer rates) for cloud electrication 
!-------------------------------------------------------------------------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_CST,            ONLY: CST_T
USE MODE_BUDGET_IAL, ONLY: TBUDGETDATA_IAL
USE MODD_BUDGET,     ONLY: NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RI, NBUDGET_SV1, TBUDGETDATA_PTR, TBUDGETCONF_T
USE MODD_NSV,        ONLY : NSV_T
USE MODD_NEB_N,      ONLY : NEB_T


USE MODE_LIMA_CCN_ACTIVATION, ONLY: LIMA_CCN_ACTIVATION
USE MODE_LIMA_CCN_HOM_FREEZING, ONLY: LIMA_CCN_HOM_FREEZING
USE MODE_LIMA_MEYERS_NUCLEATION, ONLY: LIMA_MEYERS_NUCLEATION
USE MODE_LIMA_PHILLIPS_IFN_NUCLEATION, ONLY: LIMA_PHILLIPS_IFN_NUCLEATION
USE MODE_LIMA_ICE4_NUCLEATION, ONLY: LIMA_ICE4_NUCLEATION
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
!-------------------------------------------------------------------------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(NSV_T),              INTENT(IN)    :: TNSV
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
TYPE(CST_T),              INTENT(IN)    :: CST
TYPE(NEB_T),              INTENT(IN)    :: NEBN
TYPE(TBUDGETCONF_T),      INTENT(IN)    :: BUCONF
TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS 
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step 
!TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE     ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODJ     ! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PT         ! Temperature
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PDTHRAD    ! Radiative temperature tendency
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PW_NU      ! updraft velocity used for
REAL, DIMENSION(D%NIJT, D%NKT ,TNSV%NSV), INTENT(INOUT) :: PAERO  ! Aerosol concentration
REAL, DIMENSION(D%NIJT, D%NKT, 10), INTENT(IN)    :: PSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(D%NIJT, D%NKT, KSP+KCARB+KSOA), INTENT(IN)    :: PMI
CHARACTER(LEN=4),         INTENT(IN)    :: HACTCCN  ! kind of CCN activation
INTEGER,                  INTENT(IN)    :: KCARB, KSOA, KSP ! for array size declarations
LOGICAL,                  INTENT(IN)    :: ODUST, OSALT, OORILAM
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHT       ! Theta at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRIT       ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRST       ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRGT       ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHT       ! Hail m.r. at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCCT       ! Cloud water conc. at t 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCRT       ! Rain water conc. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCIT       ! Prinstine ice conc. at t
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(INOUT) :: PCIT_SHAPE ! Ice crystal conc. at t for different shapes !++cb--
!
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_CCN), INTENT(INOUT) :: PNFT       ! CCN C. available at t
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_CCN), INTENT(INOUT) :: PNAT       ! CCN C. activated at t
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_IFN), INTENT(INOUT) :: PIFT       ! IFN C. available at t
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_IFN), INTENT(INOUT) :: PINT       ! IFN C. activated at t
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_IMM), INTENT(INOUT) :: PNIT       ! Coated IFN activated at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PNHT       ! CCN hom. freezing
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCLDFR     ! Cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PICEFR     ! Ice fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PPRCFR     ! Precipitation fraction
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTOT_RV_HENU  ! Mixing ratio change due to HENU
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTOT_RC_HINC  ! Mixing ratio change due to HINC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTOT_RI_HIND  ! Mixing ratio change due to HIND
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTOT_RV_HONH  ! Mixing ratio change due to HONH
!
!-------------------------------------------------------------------------------
!
!REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3)) :: Z_TH_HIND, Z_RI_HIND, Z_CI_HIND, Z_TH_HINC, Z_RC_HINC, Z_CC_HINC
REAL, DIMENSION(D%NIJT,D%NKT) :: Z_TH_HIND, Z_CI_HIND, Z_TH_HINC, Z_CC_HINC, &
                                 ZLSFACT, ZRVHENIMR
REAL, DIMENSION(:,:,:), ALLOCATABLE :: Z_SHCI_HIND, Z_SHCI_HINC
!
INTEGER :: IDX, IL
INTEGER :: II,IJ,ISH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_NUCLEATION_PROCS', 0, ZHOOK_HANDLE)
IF ( LIMAP%LACTI .AND. LIMAP%NMOD_CCN >=1 .AND. LIMAP%NMOM_C.GE.2) THEN

  IF (.NOT.NEBN%LSUBG_COND .AND. .NOT.LIMAP%LSPRO) THEN

    IF ( BUCONF%LBU_ENABLE ) then
       IF ( BUCONF%LBUDGET_TH ) &
            CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'HENU', PTHT(:,:) * PRHODJ(:,:) / PTSTEP )
       IF ( BUCONF%LBUDGET_RV ) &
            CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'HENU', PRVT(:,:) * PRHODJ(:,:) / PTSTEP )
       IF ( BUCONF%LBUDGET_RC ) &
            CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'HENU', PRCT(:,:) * PRHODJ(:,:) / PTSTEP )
      IF ( BUCONF%LBUDGET_SV ) then
        CALL TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NC)%PTR%INIT_PHY(D, 'HENU', PCCT(:,:) * PRHODJ(:,:) / PTSTEP )
        DO IL = 1, LIMAP%NMOD_CCN
          IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_FREE - 1 + IL
          CALL TBUDGETS(IDX)%PTR%INIT_PHY(D, 'HENU', PNFT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
          IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_ACTI - 1 + IL
          CALL TBUDGETS(IDX)%PTR%INIT_PHY(D, 'HENU', PNAT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
        END DO
      END IF
    END IF

    CALL LIMA_CCN_ACTIVATION( LIMAP, LIMAW, TNSV, D, CST, NEBN,                 &
                              KCARB, KSOA, KSP, ODUST, OSALT, OORILAM,          &
                              PRHODREF, PEXNREF, PPABST, PT, PDTHRAD, PW_NU,    &
                              PAERO, PSOLORG, PMI, HACTCCN,               & 
                              PTHT, PRVT, PRCT, PCCT, PRRT, PNFT, PNAT, PCLDFR, &
                              PTOT_RV_HENU )

    IF ( BUCONF%LBU_ENABLE ) then
       IF ( BUCONF%LBUDGET_TH ) &
            CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'HENU', PTHT(:,:) * PRHODJ(:,:) / PTSTEP )
       IF ( BUCONF%LBUDGET_RV ) &
            CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'HENU', PRVT(:,:) * PRHODJ(:,:) / PTSTEP )
       IF ( BUCONF%LBUDGET_RC ) &
            CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'HENU', PRCT(:,:) * PRHODJ(:,:) / PTSTEP )
      IF ( BUCONF%LBUDGET_SV ) then
        CALL TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NC)%PTR%END_PHY(D, 'HENU', PCCT(:,:) * PRHODJ(:,:) / PTSTEP )
        DO IL = 1, LIMAP%NMOD_CCN
          IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_FREE - 1 + IL
          CALL TBUDGETS(IDX)%PTR%END_PHY(D, 'HENU', PNFT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
          IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_ACTI - 1 + IL
          CALL TBUDGETS(IDX)%PTR%END_PHY(D, 'HENU', PNAT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
        END DO
      END IF
    END IF

  END IF

  WHERE(PCLDFR(:,:)<1.E-10 .AND. PRCT(:,:)>LIMAP%XRTMIN(2) .AND. PCCT(:,:)>LIMAP%XCTMIN(2)) PCLDFR(:,:)=1.

END IF
!
!-------------------------------------------------------------------------------
!
IF ( LIMAP%LNUCL .AND. LIMAP%NMOM_I>=2 .AND. .NOT.LIMAP%LMEYERS .AND. LIMAP%NMOD_IFN >= 1 ) THEN
  IF ( BUCONF%LBU_ENABLE ) then
    IF ( BUCONF%LBUDGET_SV ) then
      DO IL = 1, LIMAP%NMOD_IFN
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_IFN_FREE - 1 + IL
        CALL TBUDGETS(IDX)%PTR%INIT_PHY(D, 'HIND', PIFT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_IFN_NUCL - 1 + IL
        CALL TBUDGETS(IDX)%PTR%INIT_PHY(D, 'HIND', PINT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
      END DO

      DO IL = 1, LIMAP%NMOD_CCN
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_ACTI - 1 + IL
        CALL TBUDGETS(IDX)%PTR%INIT_PHY(D, 'HINC', PNAT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
      END DO
      DO IL = 1, LIMAP%NMOD_IMM
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_IMM_NUCL - 1 + IL
        CALL TBUDGETS(IDX)%PTR%INIT_PHY(D, 'HINC', PNIT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
      END DO
    END IF
  END IF

  IF (LIMAP%LCRYSTAL_SHAPE) THEN
    ALLOCATE(Z_SHCI_HIND(SIZE(PT,1),SIZE(PT,2),LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_HIND(:,:,:) = 0.
    ALLOCATE(Z_SHCI_HINC(SIZE(PT,1),SIZE(PT,2),LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_HINC(:,:,:) = 0.
  ELSE
    ALLOCATE(Z_SHCI_HIND(0,0,0))
    ALLOCATE(Z_SHCI_HINC(0,0,0))
  END IF

   CALL LIMA_PHILLIPS_IFN_NUCLEATION (LIMAP, LIMAC, D, CST, PTSTEP,                                   &
                                      PRHODREF, PEXNREF, PPABST,                        &
                                      PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,         &
                                      PCCT, PCIT, PCIT_SHAPE, PNAT, PIFT, PINT, PNIT,   &
                                      Z_TH_HIND, PTOT_RI_HIND, Z_CI_HIND, Z_SHCI_HIND,  &
                                      Z_TH_HINC, PTOT_RC_HINC, Z_CC_HINC, Z_SHCI_HINC,  &
                                      PICEFR                                            )
  WHERE(PICEFR(:,:)<1.E-10 .AND. PRIT(:,:)>LIMAP%XRTMIN(4) .AND. PCIT(:,:)>LIMAP%XCTMIN(4)) PICEFR(:,:)=1.
!
  IF ( BUCONF%LBU_ENABLE ) then
     IF ( BUCONF%LBUDGET_TH ) &
          CALL  TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'HIND',     Z_TH_HIND(:,:) * PRHODJ(:,:) / PTSTEP )
     IF ( BUCONF%LBUDGET_RV ) &
          CALL  TBUDGETS(NBUDGET_RV)%PTR%ADD_PHY(D, 'HIND', -PTOT_RI_HIND(:,:) * PRHODJ(:,:) / PTSTEP )
     IF ( BUCONF%LBUDGET_RI ) &
          CALL  TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'HIND',  PTOT_RI_HIND(:,:) * PRHODJ(:,:) / PTSTEP )
    IF ( BUCONF%LBUDGET_SV ) THEN
      IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN      
        CALL  TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NI)%PTR%ADD_PHY(D, 'HIND', Z_CI_HIND(:,:) * PRHODJ(:,:) / PTSTEP )
      ELSE
        DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
          CALL  TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NI+ISH-1)%PTR%ADD_PHY(D, 'HIND', &
                                    Z_SHCI_HIND(:, :, ISH) * PRHODJ(:, :) / PTSTEP )
        END DO
      END IF
      DO IL = 1, LIMAP%NMOD_IFN
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_IFN_FREE - 1 + IL
        CALL TBUDGETS(IDX)%PTR%END_PHY(D, 'HIND', PIFT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_IFN_NUCL - 1 + IL
        CALL TBUDGETS(IDX)%PTR%END_PHY(D, 'HIND', PINT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
      END DO
    END IF

    IF ( BUCONF%LBUDGET_TH ) &
         CALL  TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'HINC',     Z_TH_HINC(:,:) * PRHODJ(:,:) / PTSTEP )
    IF ( BUCONF%LBUDGET_RC ) &
         CALL  TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'HINC',  PTOT_RC_HINC(:,:) * PRHODJ(:,:) / PTSTEP )
    IF ( BUCONF%LBUDGET_RI ) &
         CALL  TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'HINC', -PTOT_RC_HINC(:,:) * PRHODJ(:,:) / PTSTEP )
    IF ( BUCONF%LBUDGET_SV ) then
       IF (LIMAP%NMOM_C.GE.2) then
          CALL  TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NC)%PTR%ADD_PHY(D, 'HINC',  Z_CC_HINC(:,:) * PRHODJ(:,:) / PTSTEP )
       END IF
      IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
          CALL  TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NI)%PTR%ADD_PHY(D, 'HINC', -Z_CC_HINC(:,:) * PRHODJ(:,:) / PTSTEP )
      ELSE
        DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
          CALL  TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NI+ISH-1)%PTR%ADD_PHY(D, 'HINC', &
                                    Z_SHCI_HIND(:, :, ISH) * PRHODJ(:, :) / PTSTEP )
        END DO
      END IF
      DO IL = 1, LIMAP%NMOD_CCN
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_ACTI - 1 + IL
        CALL TBUDGETS(IDX)%PTR%END_PHY(D, 'HINC', PNAT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
      END DO
      DO IL = 1, LIMAP%NMOD_IMM
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_IMM_NUCL - 1 + IL
        CALL TBUDGETS(IDX)%PTR%END_PHY(D, 'HINC', PNIT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
      END DO
    END IF
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!++cb-- 26/01/24 : pour l'instant Meyers n'est pas traite avec les formes des cristaux
IF (LIMAP%LNUCL .AND. LIMAP%NMOM_I>=2 .AND. LIMAP%LMEYERS) THEN
   CALL LIMA_MEYERS_NUCLEATION (LIMAP, LIMAC, D, CST, PTSTEP,               &
                                PRHODREF, PEXNREF, PPABST,                  &
                                PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,   &
                                PCCT, PCIT, PINT,                           &
                                Z_TH_HIND, PTOT_RI_HIND, Z_CI_HIND,         &
                                Z_TH_HINC, PTOT_RC_HINC, Z_CC_HINC,         &
                                PICEFR                                      )
  WHERE(PICEFR(:,:)<1.E-10 .AND. PRIT(:,:)>LIMAP%XRTMIN(4) .AND. PCIT(:,:)>LIMAP%XCTMIN(4)) PICEFR(:,:)=1.
  !
  IF ( BUCONF%LBU_ENABLE ) then
     IF ( BUCONF%LBUDGET_TH ) &
          CALL  TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'HIND',     Z_TH_HIND(:,:) * PRHODJ(:,:) / PTSTEP )
     IF ( BUCONF%LBUDGET_RV ) &
          CALL  TBUDGETS(NBUDGET_RV)%PTR%ADD_PHY(D, 'HIND', -PTOT_RI_HIND(:,:) * PRHODJ(:,:) / PTSTEP )
     IF ( BUCONF%LBUDGET_RI ) &
          CALL  TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'HIND',  PTOT_RI_HIND(:,:) * PRHODJ(:,:) / PTSTEP )
    IF ( BUCONF%LBUDGET_SV ) then
      CALL  TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NI)%PTR%ADD_PHY(D, 'HIND', Z_CI_HIND(:,:) * PRHODJ(:,:) / PTSTEP )
      IF (LIMAP%NMOD_IFN > 0 ) &
        CALL  TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_IFN_NUCL)%PTR%ADD_PHY(D, 'HIND', &
                               Z_CI_HIND(:,:) * PRHODJ(:,:) / PTSTEP )
    END IF

    IF ( BUCONF%LBUDGET_TH ) &
         CALL  TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'HINC',     Z_TH_HINC(:,:) * PRHODJ(:,:) / PTSTEP )
    IF ( BUCONF%LBUDGET_RC ) &
         CALL  TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'HINC',  PTOT_RC_HINC(:,:) * PRHODJ(:,:) / PTSTEP )
    IF ( BUCONF%LBUDGET_RI ) &
         CALL  TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'HINC', -PTOT_RC_HINC(:,:) * PRHODJ(:,:) / PTSTEP )
    IF ( BUCONF%LBUDGET_SV ) then
       IF (LIMAP%NMOM_C.GE.2) then
          CALL  TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NC)%PTR%ADD_PHY(D, 'HINC',  Z_CC_HINC(:,:) * PRHODJ(:,:) / PTSTEP )
       END IF
      CALL  TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NI)%PTR%ADD_PHY(D, 'HINC', -Z_CC_HINC(:,:) * PRHODJ(:,:) / PTSTEP )
      IF (LIMAP%NMOD_IFN > 0 ) &
        CALL  TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_IFN_NUCL)%PTR%ADD_PHY(D, 'HINC', &
                               -Z_CC_HINC(:,:) * PRHODJ(:,:) / PTSTEP )
    END IF
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!++cb-- pour l'instant, on ne recupere pas cette tendance
! actuellement, les echanges vapeur-->glace/eau lies a la nucleation ne sont pas traites dans l'electrisation
IF (LIMAP%LNUCL .AND. LIMAP%NMOM_I.EQ.1) THEN
  WHERE(PICEFR(:,:)<1.E-10 .AND. PRIT(:,:)>LIMAP%XRTMIN(4) .AND. PCIT(:,:)>LIMAP%XCTMIN(4)) PICEFR(:,:)=1.
  !
  ZLSFACT(:,:)=(CST%XLSTT+(CST%XCPV-CST%XCI)*(PT(:,:)-CST%XTT)) / &
     ( ( CST%XCPD +                                  &
         CST%XCPV*PRVT(:,:) +                      &
         CST%XCL*(PRCT(:,:)+PRRT(:,:)) +         &
         CST%XCI*(PRIT(:,:)+PRST(:,:)+PRGT(:,:)+PRHT(:,:)) ) * PEXNREF(:,:) ) 
  DO II = 1, D%NIJT
     CALL LIMA_ICE4_NUCLEATION(LIMAP, LIMAC, CST, D%NKT, &
          PTHT(II,:), PPABST(II,:), PRHODREF(II,:), PEXNREF(II,:), ZLSFACT(II,:), PT(II,:), &
          PRVT(II,:), &
          PCIT(II,:), ZRVHENIMR(II,:) )
  END DO
  !
  PRIT(:,:)=PRIT(:,:)+ZRVHENIMR(:,:)
  PTHT(:,:)=PTHT(:,:)+ZRVHENIMR(:,:)*ZLSFACT(:,:)
  PRVT(:,:)=PRVT(:,:)-ZRVHENIMR(:,:)
  !
  IF ( BUCONF%LBU_ENABLE ) then
     IF ( BUCONF%LBUDGET_TH ) &
          CALL  TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'HIND',  ZRVHENIMR(:,:)*ZLSFACT(:,:) * PRHODJ(:,:) / PTSTEP )
     IF ( BUCONF%LBUDGET_RV ) &
          CALL  TBUDGETS(NBUDGET_RV)%PTR%ADD_PHY(D, 'HIND', -ZRVHENIMR(:,:) * PRHODJ(:,:) / PTSTEP )
     IF ( BUCONF%LBUDGET_RI ) &
          CALL  TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'HIND',  ZRVHENIMR(:,:) * PRHODJ(:,:) / PTSTEP )
  END IF
!  Z_TH_HINC=0.
!  Z_RC_HINC=0.
!  Z_CC_HINC=0.
!  !
!  if ( BUCONF%lbu_enable ) then
!    if ( BUCONF%lbudget_th ) call  TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'HIND',  z_th_hind(:,:) * prhodj(:,:) / ptstep )
!    if ( BUCONF%lbudget_rv ) call  TBUDGETS(NBUDGET_RV)%PTR%ADD_PHY(D, 'HIND', -z_ri_hind(:,:) * prhodj(:,:) / ptstep )
!    if ( BUCONF%lbudget_ri ) call  TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'HIND',  z_ri_hind(:,:) * prhodj(:,:) / ptstep )
!    if ( BUCONF%lbudget_sv ) then
!      call  TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ni)%PTR%ADD_PHY(D, 'HIND', z_ci_hind(:,:) * prhodj(:,:) / ptstep )
!      if (nmod_ifn > 0 ) &
!        call  TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl)%PTR%ADD_PHY(D, 'HIND', &
!                               z_ci_hind(:,:) * prhodj(:,:) / ptstep )
!    end if
!
!    if ( BUCONF%lbudget_th ) call  TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'HINC',  z_th_hinc(:,:) * prhodj(:,:) / ptstep )
!    if ( BUCONF%lbudget_rc ) call  TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'HINC',  z_rc_hinc(:,:) * prhodj(:,:) / ptstep )
!    if ( BUCONF%lbudget_ri ) call  TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'HINC', -z_rc_hinc(:,:) * prhodj(:,:) / ptstep )
!    if ( BUCONF%lbudget_sv ) then
!      call  TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nc)%PTR%ADD_PHY(D, 'HINC',  z_cc_hinc(:,:) * prhodj(:,:) / ptstep )
!      call  TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ni)%PTR%ADD_PHY(D, 'HINC', -z_cc_hinc(:,:) * prhodj(:,:) / ptstep )
!      if (nmod_ifn > 0 ) &
!        call  TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl)%PTR%ADD_PHY(D, 'HINC', &
!                               -z_cc_hinc(:,:) * prhodj(:,:) / ptstep )
!    end if
!  end if
END IF
!
!-------------------------------------------------------------------------------
!
IF ( LIMAP%LNUCL .AND. LIMAP%LHHONI .AND. LIMAP%NMOD_CCN >= 1 .AND. LIMAP%NMOM_I.GE.2) THEN
  IF ( BUCONF%LBU_ENABLE ) then
     IF ( BUCONF%LBUDGET_TH ) &
          CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'HONH', PTHT(:,:) * PRHODJ(:,:) / PTSTEP )
     IF ( BUCONF%LBUDGET_RV ) &
          CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'HONH', PRVT(:,:) * PRHODJ(:,:) / PTSTEP )
     IF ( BUCONF%LBUDGET_RI ) &
          CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'HONH', PRIT(:,:) * PRHODJ(:,:) / PTSTEP )
    IF ( BUCONF%LBUDGET_SV ) THEN
      IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
        CALL TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NI)%PTR%INIT_PHY(D,'HONH',PCIT(:,:)*PRHODJ(:,:)/PTSTEP )
      ELSE
        DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
          CALL TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NI+ISH-1)%PTR%INIT_PHY(D,'HONH', &
                                     PCIT_SHAPE(:, :, ISH)*PRHODJ(:, :)/PTSTEP )
        END DO
      END IF
      DO IL = 1, LIMAP%NMOD_CCN
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_FREE - 1 + IL
        CALL TBUDGETS(IDX)%PTR%INIT_PHY(D, 'HONH', PNFT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
      END DO
      CALL TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_HOM_HAZE)%PTR%INIT_PHY(D,'HONH',PNHT(:,:)*PRHODJ(:,:)/PTSTEP)
    END IF
  END IF

  CALL LIMA_CCN_HOM_FREEZING (LIMAP, LIMAC, D, CST, PRHODREF, PEXNREF, PPABST, PW_NU, &
                              PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                              PCCT, PCRT, PCIT, PCIT_SHAPE, PNFT, PNHT, &
                              PICEFR, PTOT_RV_HONH                      )
  WHERE(PICEFR(:,:)<1.E-10 .AND. PRIT(:,:)>LIMAP%XRTMIN(4) .AND. PCIT(:,:)>LIMAP%XCTMIN(4)) PICEFR(:,:)=1.
!
  IF ( BUCONF%LBU_ENABLE ) then
     IF ( BUCONF%LBUDGET_TH ) &
          CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'HONH', PTHT(:,:) * PRHODJ(:,:) / PTSTEP )
     IF ( BUCONF%LBUDGET_RV ) &
          CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'HONH', PRVT(:,:) * PRHODJ(:,:) / PTSTEP )
     IF ( BUCONF%LBUDGET_RI ) &
          CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'HONH', PRIT(:,:) * PRHODJ(:,:) / PTSTEP )
    IF ( BUCONF%LBUDGET_SV ) then
      IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
        CALL TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NI)%PTR%END_PHY(D,'HONH',PCIT(:,:)*PRHODJ(:,:)/PTSTEP )
      ELSE
        DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
          CALL TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_NI+ISH-1)%PTR%END_PHY(D,'HONH', &
                                    PCIT_SHAPE(:, :, ISH)*PRHODJ(:, :)/PTSTEP )
        END DO
      END IF
      DO IL = 1, LIMAP%NMOD_CCN
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_FREE - 1 + IL
        CALL TBUDGETS(IDX)%PTR%END_PHY(D, 'HONH', PNFT(:,:, IL) * PRHODJ(:,:) / PTSTEP )
      END DO
      CALL TBUDGETS(NBUDGET_SV1-1+TNSV%NSV_LIMA_HOM_HAZE)%PTR%END_PHY(D,'HONH',PNHT(:,:)*PRHODJ(:,:)/PTSTEP)
    END IF
  END IF
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_NUCLEATION_PROCS', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_NUCLEATION_PROCS
END MODULE MODE_LIMA_NUCLEATION_PROCS
