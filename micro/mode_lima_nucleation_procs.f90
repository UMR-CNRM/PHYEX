!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_NUCLEATION_PROCS
  IMPLICIT NONE
CONTAINS
!     #############################################################################
  SUBROUTINE LIMA_NUCLEATION_PROCS (D, CST, BUCONF, TBUDGETS, KBUDGETS,             &
                                    PTSTEP, PRHODJ,                                 &
                                    PRHODREF, PEXNREF, PPABST, PT, PDTHRAD, PW_NU,  &
                                    PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT, &
                                    PCCT, PCRT, PCIT,                               &
                                    PNFT, PNAT, PIFT, PINT, PNIT, PNHT,             &
                                    PCLDFR, PICEFR, PPRCFR                          )
!     #############################################################################
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
!-------------------------------------------------------------------------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_BUDGET,   ONLY: TBUDGETDATA, TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
use modd_budget,     only: NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RI, NBUDGET_SV1
!USE MODD_IO,         ONLY: TFILEDATA
USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT
USE MODD_NSV,        ONLY : NSV_LIMA_NC, NSV_LIMA_NR, NSV_LIMA_CCN_FREE, NSV_LIMA_CCN_ACTI, &
                            NSV_LIMA_NI, NSV_LIMA_IFN_FREE, NSV_LIMA_IFN_NUCL, NSV_LIMA_IMM_NUCL, NSV_LIMA_HOM_HAZE
USE MODD_PARAM_LIMA, ONLY : LNUCL, LMEYERS, LACTI, LHHONI,  &
                            NMOD_CCN, NMOD_IFN, NMOD_IMM, XCTMIN, XRTMIN, LSPRO, NMOM_I, NMOM_C
USE MODD_TURB_n,     ONLY : LSUBG_COND

USE MODE_BUDGET_PHY,     ONLY: BUDGET_STORE_ADD_PHY, BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY

USE MODE_LIMA_CCN_ACTIVATION, ONLY: LIMA_CCN_ACTIVATION
USE MODE_LIMA_CCN_HOM_FREEZING, ONLY: LIMA_CCN_HOM_FREEZING
USE MODE_LIMA_MEYERS_NUCLEATION, ONLY: LIMA_MEYERS_NUCLEATION
USE MODE_LIMA_PHILLIPS_IFN_NUCLEATION, ONLY: LIMA_PHILLIPS_IFN_NUCLEATION
USE MODE_LIMA_ICE4_NUCLEATION, ONLY: LIMA_ICE4_NUCLEATION
!
!-------------------------------------------------------------------------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS 
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step 
!TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE     ! Output file
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PT         ! Temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDTHRAD    ! Radiative temperature tendency
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU      ! updraft velocity used for
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHT       ! Theta at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIT       ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST       ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT       ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHT       ! Hail m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCT       ! Cloud water conc. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT       ! Rain water conc. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT       ! Prinstine ice conc. at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNFT       ! CCN C. available at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNAT       ! CCN C. activated at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PIFT       ! IFN C. available at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PINT       ! IFN C. activated at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PNIT       ! Coated IFN activated at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNHT       ! CCN hom. freezing
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCLDFR     ! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR     ! Ice fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PPRCFR     ! Precipitation fraction
!
!-------------------------------------------------------------------------------
!
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3))          :: Z_TH_HIND, Z_RI_HIND, Z_CI_HIND, Z_TH_HINC, Z_RC_HINC, Z_CC_HINC
REAL, DIMENSION(SIZE(PT,1),SIZE(PT,2),SIZE(PT,3))          :: ZLSFACT, ZRVHENIMR
!
integer :: idx, jl
INTEGER :: JI,JJ
!
!-------------------------------------------------------------------------------
!
IF ( LACTI .AND. NMOD_CCN >=1 .AND. NMOM_C.GE.2) THEN

  IF (.NOT.LSUBG_COND .AND. .NOT.LSPRO) THEN

    if ( BUCONF%lbu_enable ) then
       if ( BUCONF%lbudget_th ) &
            call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'HENU', ptht(:, :, :) * prhodj(:, :, :) / ptstep )
       if ( BUCONF%lbudget_rv ) &
            call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'HENU', prvt(:, :, :) * prhodj(:, :, :) / ptstep )
       if ( BUCONF%lbudget_rc ) &
            call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'HENU', prct(:, :, :) * prhodj(:, :, :) / ptstep )
      if ( BUCONF%lbudget_sv ) then
        call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1-1+nsv_lima_nc), 'HENU', pcct(:, :, :) * prhodj(:, :, :) / ptstep )
        do jl = 1, nmod_ccn
          idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_free - 1 + jl
          call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'HENU', pnft(:, :, :, jl) * prhodj(:, :, :) / ptstep )
          idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_acti - 1 + jl
          call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'HENU', pnat(:, :, :, jl) * prhodj(:, :, :) / ptstep )
        end do
      end if
    end if

    CALL LIMA_CCN_ACTIVATION( CST,                                      &
                              PRHODREF, PEXNREF, PPABST, PT, PDTHRAD, PW_NU,    &
                              PTHT, PRVT, PRCT, PCCT, PRRT, PNFT, PNAT, PCLDFR  )
    if ( BUCONF%lbu_enable ) then
       if ( BUCONF%lbudget_th ) &
            call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'HENU', ptht(:, :, :) * prhodj(:, :, :) / ptstep )
       if ( BUCONF%lbudget_rv ) &
            call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'HENU', prvt(:, :, :) * prhodj(:, :, :) / ptstep )
       if ( BUCONF%lbudget_rc ) &
            call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'HENU', prct(:, :, :) * prhodj(:, :, :) / ptstep )
      if ( BUCONF%lbudget_sv ) then
        call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1-1+nsv_lima_nc), 'HENU', pcct(:, :, :) * prhodj(:, :, :) / ptstep )
        do jl = 1, nmod_ccn
          idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_free - 1 + jl
          call BUDGET_STORE_END_PHY(D, TBUDGETS(idx), 'HENU', pnft(:, :, :, jl) * prhodj(:, :, :) / ptstep )
          idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_acti - 1 + jl
          call BUDGET_STORE_END_PHY(D, TBUDGETS(idx), 'HENU', pnat(:, :, :, jl) * prhodj(:, :, :) / ptstep )
        end do
      end if
    end if

  END IF

  WHERE(PCLDFR(:,:,:)<1.E-10 .AND. PRCT(:,:,:)>XRTMIN(2) .AND. PCCT(:,:,:)>XCTMIN(2)) PCLDFR(:,:,:)=1.

END IF
!
!-------------------------------------------------------------------------------
!
IF ( LNUCL .AND. NMOM_I>=2 .AND. .NOT.LMEYERS .AND. NMOD_IFN >= 1 ) THEN
  if ( BUCONF%lbu_enable ) then
    if ( BUCONF%lbudget_sv ) then
      do jl = 1, nmod_ifn
        idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_free - 1 + jl
        call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'HIND', pift(:, :, :, jl) * prhodj(:, :, :) / ptstep )
        idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl - 1 + jl
        call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'HIND', pint(:, :, :, jl) * prhodj(:, :, :) / ptstep )
      end do

      do jl = 1, nmod_ccn
        idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_acti - 1 + jl
        call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'HINC', pnat(:, :, :, jl) * prhodj(:, :, :) / ptstep )
      end do
      do jl = 1, nmod_imm
        idx = NBUDGET_SV1 - 1 + nsv_lima_imm_nucl - 1 + jl
        call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'HINC', pnit(:, :, :, jl) * prhodj(:, :, :) / ptstep )
      end do
    end if
  end if

   CALL LIMA_PHILLIPS_IFN_NUCLEATION (CST, PTSTEP,                                      &
                                      PRHODREF, PEXNREF, PPABST,                        &
                                      PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,         &
                                      PCCT, PCIT, PNAT, PIFT, PINT, PNIT,               &
                                      Z_TH_HIND, Z_RI_HIND, Z_CI_HIND,                  &
                                      Z_TH_HINC, Z_RC_HINC, Z_CC_HINC,                  &
                                      PICEFR                                            )
  WHERE(PICEFR(:,:,:)<1.E-10 .AND. PRIT(:,:,:)>XRTMIN(4) .AND. PCIT(:,:,:)>XCTMIN(4)) PICEFR(:,:,:)=1.
!
  if ( BUCONF%lbu_enable ) then
     if ( BUCONF%lbudget_th ) &
          call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'HIND',  z_th_hind(:, :, :) * prhodj(:, :, :) / ptstep )
     if ( BUCONF%lbudget_rv ) &
          call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RV), 'HIND', -z_ri_hind(:, :, :) * prhodj(:, :, :) / ptstep )
     if ( BUCONF%lbudget_ri ) &
          call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'HIND',  z_ri_hind(:, :, :) * prhodj(:, :, :) / ptstep )
    if ( BUCONF%lbudget_sv ) then
      call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1-1+nsv_lima_ni), 'HIND', z_ci_hind(:, :, :) * prhodj(:, :, :) / ptstep )
      do jl = 1, nmod_ifn
        idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_free - 1 + jl
        call BUDGET_STORE_END_PHY(D, TBUDGETS(idx), 'HIND', pift(:, :, :, jl) * prhodj(:, :, :) / ptstep )
        idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl - 1 + jl
        call BUDGET_STORE_END_PHY(D, TBUDGETS(idx), 'HIND', pint(:, :, :, jl) * prhodj(:, :, :) / ptstep )
      end do
    end if

    if ( BUCONF%lbudget_th ) &
         call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'HINC',  z_th_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
    if ( BUCONF%lbudget_rc ) &
         call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'HINC',  z_rc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
    if ( BUCONF%lbudget_ri ) &
         call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'HINC', -z_rc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
    if ( BUCONF%lbudget_sv ) then
       if (nmom_c.ge.2) then
          call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1-1+nsv_lima_nc), 'HINC',  z_cc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
       end if
      call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1-1+nsv_lima_ni), 'HINC', -z_cc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
      do jl = 1, nmod_ccn
        idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_acti - 1 + jl
        call BUDGET_STORE_END_PHY(D, TBUDGETS(idx), 'HINC', pnat(:, :, :, jl) * prhodj(:, :, :) / ptstep )
      end do
      do jl = 1, nmod_imm
        idx = NBUDGET_SV1 - 1 + nsv_lima_imm_nucl - 1 + jl
        call BUDGET_STORE_END_PHY(D, TBUDGETS(idx), 'HINC', pnit(:, :, :, jl) * prhodj(:, :, :) / ptstep )
      end do
    end if
  end if
END IF
!
!-------------------------------------------------------------------------------
!
IF (LNUCL .AND. NMOM_I>=2 .AND. LMEYERS) THEN
   CALL LIMA_MEYERS_NUCLEATION (CST, PTSTEP,                                &
                                PRHODREF, PEXNREF, PPABST,                  &
                                PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,   &
                                PCCT, PCIT, PINT,                           &
                                Z_TH_HIND, Z_RI_HIND, Z_CI_HIND,            &
                                Z_TH_HINC, Z_RC_HINC, Z_CC_HINC,            &
                                PICEFR                                      )
  WHERE(PICEFR(:,:,:)<1.E-10 .AND. PRIT(:,:,:)>XRTMIN(4) .AND. PCIT(:,:,:)>XCTMIN(4)) PICEFR(:,:,:)=1.
  !
  if ( BUCONF%lbu_enable ) then
     if ( BUCONF%lbudget_th ) &
          call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'HIND',  z_th_hind(:, :, :) * prhodj(:, :, :) / ptstep )
     if ( BUCONF%lbudget_rv ) &
          call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RV), 'HIND', -z_ri_hind(:, :, :) * prhodj(:, :, :) / ptstep )
     if ( BUCONF%lbudget_ri ) &
          call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'HIND',  z_ri_hind(:, :, :) * prhodj(:, :, :) / ptstep )
    if ( BUCONF%lbudget_sv ) then
      call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1-1+nsv_lima_ni), 'HIND', z_ci_hind(:, :, :) * prhodj(:, :, :) / ptstep )
      if (nmod_ifn > 0 ) &
        call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1-1+nsv_lima_ifn_nucl), 'HIND', &
                               z_ci_hind(:, :, :) * prhodj(:, :, :) / ptstep )
    end if

    if ( BUCONF%lbudget_th ) &
         call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'HINC',  z_th_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
    if ( BUCONF%lbudget_rc ) &
         call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'HINC',  z_rc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
    if ( BUCONF%lbudget_ri ) &
         call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'HINC', -z_rc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
    if ( BUCONF%lbudget_sv ) then
       if (nmom_c.ge.2) then
          call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1-1+nsv_lima_nc), 'HINC',  z_cc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
       end if
      call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1-1+nsv_lima_ni), 'HINC', -z_cc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
      if (nmod_ifn > 0 ) &
        call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1-1+nsv_lima_ifn_nucl), 'HINC', &
                               -z_cc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
    end if
  end if
END IF
!
!-------------------------------------------------------------------------------
!
IF (LNUCL .AND. NMOM_I.EQ.1) THEN
  WHERE(PICEFR(:,:,:)<1.E-10 .AND. PRIT(:,:,:)>XRTMIN(4) .AND. PCIT(:,:,:)>XCTMIN(4)) PICEFR(:,:,:)=1.
  !
  ZLSFACT(:,:,:)=(CST%XLSTT+(CST%XCPV-CST%XCI)*(PT(:,:,:)-CST%XTT)) / &
     ( ( CST%XCPD +                                  &
         CST%XCPV*PRVT(:,:,:) +                      &
         CST%XCL*(PRCT(:,:,:)+PRRT(:,:,:)) +         &
         CST%XCI*(PRIT(:,:,:)+PRST(:,:,:)+PRGT(:,:,:)+PRHT(:,:,:)) ) * PEXNREF(:,:,:) ) 
  DO JI = 1, SIZE(PTHT,1)
     DO JJ = 1, SIZE(PTHT,2)
        CALL LIMA_ICE4_NUCLEATION(CST, SIZE(PTHT,3), &
             PTHT(JI,JJ,:), PPABST(JI,JJ,:), PRHODREF(JI,JJ,:), PEXNREF(JI,JJ,:), ZLSFACT(JI,JJ,:), PT(JI,JJ,:), &
             PRVT(JI,JJ,:), &
             PCIT(JI,JJ,:), ZRVHENIMR(JI,JJ,:) )
     END DO
  END DO
  !
  PRIT(:,:,:)=PRIT(:,:,:)+ZRVHENIMR(:,:,:)
  PTHT(:,:,:)=PTHT(:,:,:)+ZRVHENIMR(:,:,:)*ZLSFACT(:,:,:)
  PRVT(:,:,:)=PRVT(:,:,:)-ZRVHENIMR(:,:,:)
  !
  if ( BUCONF%lbu_enable ) then
     if ( BUCONF%lbudget_th ) &
          call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'HIND',  ZRVHENIMR(:,:,:)*ZLSFACT(:,:,:) * prhodj(:, :, :) / ptstep )
     if ( BUCONF%lbudget_rv ) &
          call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RV), 'HIND', -ZRVHENIMR(:, :, :) * prhodj(:, :, :) / ptstep )
     if ( BUCONF%lbudget_ri ) &
          call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'HIND',  ZRVHENIMR(:, :, :) * prhodj(:, :, :) / ptstep )
  end if
!  Z_TH_HINC=0.
!  Z_RC_HINC=0.
!  Z_CC_HINC=0.
!  !
!  if ( BUCONF%lbu_enable ) then
!    if ( BUCONF%lbudget_th ) call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'HIND',  z_th_hind(:, :, :) * prhodj(:, :, :) / ptstep )
!    if ( BUCONF%lbudget_rv ) call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RV), 'HIND', -z_ri_hind(:, :, :) * prhodj(:, :, :) / ptstep )
!    if ( BUCONF%lbudget_ri ) call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'HIND',  z_ri_hind(:, :, :) * prhodj(:, :, :) / ptstep )
!    if ( BUCONF%lbudget_sv ) then
!      call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HIND', z_ci_hind(:, :, :) * prhodj(:, :, :) / ptstep )
!      if (nmod_ifn > 0 ) &
!        call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl), 'HIND', &
!                               z_ci_hind(:, :, :) * prhodj(:, :, :) / ptstep )
!    end if
!
!    if ( BUCONF%lbudget_th ) call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'HINC',  z_th_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
!    if ( BUCONF%lbudget_rc ) call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'HINC',  z_rc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
!    if ( BUCONF%lbudget_ri ) call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'HINC', -z_rc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
!    if ( BUCONF%lbudget_sv ) then
!      call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nc), 'HINC',  z_cc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
!      call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HINC', -z_cc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
!      if (nmod_ifn > 0 ) &
!        call BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl), 'HINC', &
!                               -z_cc_hinc(:, :, :) * prhodj(:, :, :) / ptstep )
!    end if
!  end if
END IF
!
!-------------------------------------------------------------------------------
!
IF ( LNUCL .AND. LHHONI .AND. NMOD_CCN >= 1 .AND. NMOM_I.GE.2) THEN
  if ( BUCONF%lbu_enable ) then
     if ( BUCONF%lbudget_th ) &
          call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'HONH', PTHT(:, :, :) * prhodj(:, :, :) / ptstep )
     if ( BUCONF%lbudget_rv ) &
          call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'HONH', PRVT(:, :, :) * prhodj(:, :, :) / ptstep )
     if ( BUCONF%lbudget_ri ) &
          call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'HONH', PRIT(:, :, :) * prhodj(:, :, :) / ptstep )
    if ( BUCONF%lbudget_sv ) then
      call BUDGET_STORE_INIT_PHY(D,TBUDGETS(NBUDGET_SV1-1+nsv_lima_ni),'HONH',PCIT(:, :, :)*prhodj(:, :, :)/ptstep )
      do jl = 1, nmod_ccn
        idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_free - 1 + jl
        call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'HONH', PNFT(:, :, :, jl) * prhodj(:, :, :) / ptstep )
      end do
      call BUDGET_STORE_INIT_PHY(D,TBUDGETS(NBUDGET_SV1-1+nsv_lima_hom_haze),'HONH',PNHT(:, :, :)*prhodj(:, :, :)/ptstep)
    end if
  end if

  CALL LIMA_CCN_HOM_FREEZING (CST, PRHODREF, PEXNREF, PPABST, PW_NU,    &
                              PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                              PCCT, PCRT, PCIT, PNFT, PNHT,             &
                              PICEFR                                    )
  WHERE(PICEFR(:,:,:)<1.E-10 .AND. PRIT(:,:,:)>XRTMIN(4) .AND. PCIT(:,:,:)>XCTMIN(4)) PICEFR(:,:,:)=1.
!
  if ( BUCONF%lbu_enable ) then
     if ( BUCONF%lbudget_th ) &
          call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'HONH', PTHT(:, :, :) * prhodj(:, :, :) / ptstep )
     if ( BUCONF%lbudget_rv ) &
          call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'HONH', PRVT(:, :, :) * prhodj(:, :, :) / ptstep )
     if ( BUCONF%lbudget_ri ) &
          call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'HONH', PRIT(:, :, :) * prhodj(:, :, :) / ptstep )
    if ( BUCONF%lbudget_sv ) then
      call BUDGET_STORE_END_PHY(D,TBUDGETS(NBUDGET_SV1-1+nsv_lima_ni),'HONH',PCIT(:, :, :)*prhodj(:, :, :)/ptstep )
      do jl = 1, nmod_ccn
        idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_free - 1 + jl
        call BUDGET_STORE_END_PHY(D, TBUDGETS(idx), 'HONH', PNFT(:, :, :, jl) * prhodj(:, :, :) / ptstep )
      end do
      call BUDGET_STORE_END_PHY(D,TBUDGETS(NBUDGET_SV1-1+nsv_lima_hom_haze),'HONH',PNHT(:, :, :)*prhodj(:, :, :)/ptstep)
    end if
  end if
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_NUCLEATION_PROCS
END MODULE MODE_LIMA_NUCLEATION_PROCS
