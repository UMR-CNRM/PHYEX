!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ######spl
MODULE MODI_LIMA
!      ####################
!
INTERFACE
!
   SUBROUTINE LIMA ( KKA, KKU, KKL,                                          &
                     PTSTEP, TPFILE,                                         &
                     PRHODREF, PEXNREF, PDZZ,                                &
                     PRHODJ, PPABSM, PPABST,                                 &
                     NCCN, NIFN, NIMM,                                       &
                     PDTHRAD, PTHT, PRT, PSVT, PW_NU,                        &
                     PTHS, PRS, PSVS,                                        &
                     PINPRC, PINDEP, PINPRR, PINPRI, PINPRS, PINPRG, PINPRH, &
                     PEVAP3D, PCLDFR, PICEFR, PPRCFR                         )
!
USE MODD_IO,  ONLY: TFILEDATA
USE MODD_NSV, only: NSV_LIMA_BEG
!
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE     ! Output file
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ       ! Layer thikness (m)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABSM     ! absolute pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST     ! absolute pressure at t
!
INTEGER,                  INTENT(IN)    :: NCCN       ! for array size declarations
INTEGER,                  INTENT(IN)    :: NIFN       ! for array size declarations
INTEGER,                  INTENT(IN)    :: NIMM       ! for array size declarations
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDTHRAD    ! dT/dt due to radiation
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT       ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT        ! Mixing ratios at time t
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(IN)    :: PSVT ! Concentrations at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU      ! w for CCN activation
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS       ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS        ! Mixing ratios sources
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(INOUT) :: PSVS ! Concentration sources
!
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRC     ! Cloud instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINDEP     ! Cloud droplets deposition
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRR     ! Rain instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRI     ! Rain instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRS     ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRG     ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRH     ! Rain instant precip
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PEVAP3D    ! Rain evap profile
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCLDFR     ! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR     ! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PPRCFR     ! Cloud fraction
!
END SUBROUTINE LIMA
END INTERFACE
END MODULE MODI_LIMA
!
!
!     ######spl
      SUBROUTINE LIMA ( KKA, KKU, KKL,                                          &
                        PTSTEP, TPFILE,                                         &
                        PRHODREF, PEXNREF, PDZZ,                                &
                        PRHODJ, PPABSM, PPABST,                                 &
                        NCCN, NIFN, NIMM,                                       &
                        PDTHRAD, PTHT, PRT, PSVT, PW_NU,                        &
                        PTHS, PRS, PSVS,                                        &
                        PINPRC, PINDEP, PINPRR, PINPRI, PINPRS, PINPRG, PINPRH, &
                        PEVAP3D, PCLDFR, PICEFR, PPRCFR                         )
!     ######################################################################
!
!!    PURPOSE
!!    -------
!!      Compute explicit microphysical sources using the 2-moment scheme LIMA     
!!    using the time-splitting method
!!
!!    REFERENCE
!!    ---------
!!      Vié et al. (GMD, 2016)
!!      Meso-NH scientific documentation
!!
!!    AUTHOR
!!    ------
!!      S. Riette  * CNRM *
!!      B. Vié     * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   15/03/2018
!!
!  B. Vie         02/2019: minor correction on budget
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets (no more budget calls in this subroutine)
!  P. Wautelet 26/02/2020: bugfix: corrected condition to write budget CORR_BU_RRS
!  B. Vie      03/03/2020: use DTHRAD instead of dT/dt in Smax diagnostic computation
!  P. Wautelet 28/05/2020: bugfix: correct array start for PSVT and PSVS
!  P. Wautelet 03/02/2021: budgets: add new source if LIMA splitting: CORR2
!  B. Vie         06/2021: add subgrid condensation with LIMA
!-----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
use modd_budget,          only: lbu_enable,                                                  &
                                lbudget_th, lbudget_rv, lbudget_rc, lbudget_rr, lbudget_ri,  &
                                lbudget_rs, lbudget_rg, lbudget_rh, lbudget_sv,              &
                                NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RR, NBUDGET_RI,  &
                                NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, NBUDGET_SV1,             &
                                tbudgets
USE MODD_CST,             ONLY: XCI, XCL, XCPD, XCPV, XLSTT, XLVTT, XTT, XRHOLW, XP00, XRD
USE MODD_IO,              ONLY: TFILEDATA
USE MODD_NSV,             ONLY: NSV_LIMA_BEG,                                                   &
                                NSV_LIMA_NC, NSV_LIMA_NR, NSV_LIMA_CCN_FREE, NSV_LIMA_CCN_ACTI, &
                                NSV_LIMA_NI, NSV_LIMA_IFN_FREE,                                 &
                                NSV_LIMA_IFN_NUCL, NSV_LIMA_IMM_NUCL, NSV_LIMA_HOM_HAZE
USE MODD_PARAMETERS,      ONLY: JPHEXT, JPVEXT
USE MODD_PARAM_LIMA,      ONLY: LCOLD, LRAIN, LWARM, NMOD_CCN, NMOD_IFN, NMOD_IMM, LHHONI,      &
                                LACTIT, LFEEDBACKT, NMAXITER, XMRSTEP, XTSTEP_TS,               &
                                LSEDC, LSEDI, XRTMIN, XCTMIN, LDEPOC, XVDEPOC,                  &
                                LHAIL, LSNOW
USE MODD_PARAM_LIMA_COLD, ONLY: XAI, XBI
USE MODD_PARAM_LIMA_WARM, ONLY: XLBC, XLBEXC, XAC, XBC, XAR, XBR
USE MODD_TURB_n,          ONLY: LSUBG_COND

use mode_budget,          only: Budget_store_add, Budget_store_init, Budget_store_end
use mode_tools,           only: Countjv

USE MODI_LIMA_COMPUTE_CLOUD_FRACTIONS
USE MODI_LIMA_DROPS_TO_DROPLETS_CONV
USE MODI_LIMA_INST_PROCS
USE MODI_LIMA_NUCLEATION_PROCS
USE MODI_LIMA_SEDIMENTATION
USE MODI_LIMA_TENDENCIES
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE     ! Output file
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ       ! Layer thikness (m)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ     ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABSM     ! absolute pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST     ! absolute pressure at t
!
INTEGER,                  INTENT(IN)    :: NCCN       ! for array size declarations
INTEGER,                  INTENT(IN)    :: NIFN       ! for array size declarations
INTEGER,                  INTENT(IN)    :: NIMM       ! for array size declarations
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDTHRAD    ! dT/dt due to radiation
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT       ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT        ! Mixing ratios at time t
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(IN)    :: PSVT ! Concentrations at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU      ! w for CCN activation
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS       ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS        ! Mixing ratios sources
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(INOUT) :: PSVS ! Concentration sources
!
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRC     ! Cloud instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINDEP     ! Cloud droplets deposition
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRR     ! Rain instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRI     ! Rain instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRS     ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRG     ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRH     ! Rain instant precip
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PEVAP3D    ! Rain evap profile
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCLDFR     ! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PICEFR     ! Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PPRCFR     ! Cloud fraction
!
!*       0.2   Declarations of local variables :
!
!
! Prognostic variables and sources
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3))      :: ZTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT, ZRHT
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3))      :: ZCCT, ZCRT, ZCIT
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3))      :: ZTHS, ZRVS, ZRCS, ZRRS, ZRIS, ZRSS, ZRGS, ZRHS
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3))      :: ZCCS, ZCRS, ZCIS
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),NCCN) :: ZCCNFT, ZCCNAT
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),NCCN) :: ZCCNFS, ZCCNAS
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),NIFN) :: ZIFNFT, ZIFNNT
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),NIFN) :: ZIFNFS, ZIFNNS
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),NIMM) :: ZIMMNT
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),NIMM) :: ZIMMNS
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3))      :: ZHOMFT
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3))      :: ZHOMFS

!
! Other 3D thermodynamical variables
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3))      :: ZEXN, ZT

!
! Packed prognostic & thermo variables
REAL, DIMENSION(:),   ALLOCATABLE ::                               &
     ZP1D, ZRHODREF1D, ZEXNREF1D, ZEXN1D,                          &
     ZTHT1D,                                                       &
     ZRVT1D, ZRCT1D, ZRRT1D, ZRIT1D, ZRST1D, ZRGT1D, ZRHT1D,       &
     ZCCT1D, ZCRT1D, ZCIT1D,                                       &
     ZEVAP1D
REAL, DIMENSION(:,:), ALLOCATABLE :: ZIFNN1D

!
! for each process & species inside the loop, we need 1D packed variables to store instant tendencies for hydrometeors
REAL, DIMENSION(:), ALLOCATABLE ::                          &
! mixing ratio & concentration changes by instantaneous processes (kg/kg and #/kg) :
     Z_CR_BRKU,                                             & ! spontaneous break up of drops (BRKU) : Nr
     Z_TH_HONR, Z_RR_HONR, Z_CR_HONR,                       & ! rain drops homogeneous freezing (HONR) : rr, Nr, rg=-rr, th
     Z_TH_IMLT, Z_RC_IMLT, Z_CC_IMLT,                       & ! ice melting (IMLT) : rc, Nc, ri=-rc, Ni=-Nc, th, IFNF, IFNA
! mixing ratio & concentration tendencies by continuous processes (kg/kg/s and #/kg/s) :
     Z_TH_HONC, Z_RC_HONC, Z_CC_HONC,                       & ! droplets homogeneous freezing (HONC) : rc, Nc, ri=-rc, Ni=-Nc, th
     Z_CC_SELF,                                             & ! self collection of droplets (SELF) : Nc
     Z_RC_AUTO, Z_CC_AUTO, Z_CR_AUTO,                       & ! autoconversion of cloud droplets (AUTO) : rc, Nc, rr=-rc, Nr
     Z_RC_ACCR, Z_CC_ACCR,                                  & ! accretion of droplets by rain drops (ACCR) : rc, Nc, rr=-rr
     Z_CR_SCBU,                                             & ! self collectio break up of drops (SCBU) : Nr
!      Z_TH_EVAP, Z_RC_EVAP, Z_CC_EVAP, Z_RR_EVAP, Z_CR_EVAP, & ! evaporation of rain drops (EVAP) : rv=-rr-rc, rc, Nc, rr, Nr, th
     Z_TH_EVAP, Z_RR_EVAP, & ! evaporation of rain drops (EVAP) : rv=-rr-rc, rc, Nc, rr, Nr, th
     Z_RI_CNVI, Z_CI_CNVI,                                  & ! conversion snow -> ice (CNVI) : ri, Ni, rs=-ri
     Z_TH_DEPS, Z_RS_DEPS,                                  & ! deposition of vapor on snow (DEPS) : rv=-rs, rs, th
     Z_TH_DEPI, Z_RI_DEPI,                                  & ! deposition of vapor on ice (DEPI) : rv=-ri, ri, th
     Z_RI_CNVS, Z_CI_CNVS,                                  & ! conversion ice -> snow (CNVS) : ri, Ni, rs=-ri
     Z_RI_AGGS, Z_CI_AGGS,                                  & ! aggregation of ice on snow (AGGS) : ri, Ni, rs=-ri
     Z_TH_DEPG, Z_RG_DEPG,                                  & ! deposition of vapor on graupel (DEPG) : rv=-rg, rg, th
     Z_TH_BERFI, Z_RC_BERFI,                                & ! Bergeron (BERFI) : rc, ri=-rc, th
     Z_TH_RIM, Z_RC_RIM, Z_CC_RIM, Z_RS_RIM, Z_RG_RIM,      & ! cloud droplet riming (RIM) : rc, Nc, rs, rg, th
     Z_RI_HMS, Z_CI_HMS, Z_RS_HMS,                          & ! hallett mossop snow (HMS) : ri, Ni, rs
     Z_TH_ACC, Z_RR_ACC, Z_CR_ACC, Z_RS_ACC, Z_RG_ACC,      & ! rain accretion on aggregates (ACC) : rr, Nr, rs, rg, th
     Z_RS_CMEL,                                             & ! conversion-melting (CMEL) : rs, rg=-rs
     Z_TH_CFRZ, Z_RR_CFRZ, Z_CR_CFRZ, Z_RI_CFRZ, Z_CI_CFRZ, & ! rain freezing (CFRZ) : rr, Nr, ri, Ni, rg=-rr-ri, th
     Z_TH_WETG, Z_RC_WETG, Z_CC_WETG, Z_RR_WETG, Z_CR_WETG, & ! wet growth of graupel (WETG) : rc, NC, rr, Nr, ri, Ni, rs, rg, rh, th
     Z_RI_WETG, Z_CI_WETG, Z_RS_WETG, Z_RG_WETG, Z_RH_WETG, & ! wet growth of graupel (WETG) : rc, NC, rr, Nr, ri, Ni, rs, rg, rh, th
     Z_TH_DRYG, Z_RC_DRYG, Z_CC_DRYG, Z_RR_DRYG, Z_CR_DRYG, & ! dry growth of graupel (DRYG) : rc, Nc, rr, Nr, ri, Ni, rs, rg, th
     Z_RI_DRYG, Z_CI_DRYG, Z_RS_DRYG, Z_RG_DRYG,            & ! dry growth of graupel (DRYG) : rc, Nc, rr, Nr, ri, Ni, rs, rg, th
     Z_RI_HMG, Z_CI_HMG, Z_RG_HMG,                          & ! hallett mossop graupel (HMG) : ri, Ni, rg
     Z_TH_GMLT, Z_RR_GMLT, Z_CR_GMLT,                       & ! graupel melting (GMLT) : rr, Nr, rg=-rr, th
!      Z_RC_WETH, Z_CC_WETH, Z_RR_WETH, Z_CR_WETH,            & ! wet growth of hail (WETH) : rc, Nc, rr, Nr, ri, Ni, rs, rg, rh, th
!      Z_RI_WETH, Z_CI_WETH, Z_RS_WETH, Z_RG_WETH, Z_RH_WETH, & ! wet growth of hail (WETH) : rc, Nc, rr, Nr, ri, Ni, rs, rg, rh, th
!      Z_RG_COHG,                                             & ! conversion of hail into graupel (COHG) : rg, rh
!      Z_RR_HMLT, Z_CR_HMLT                                     ! hail melting (HMLT) : rr, Nr, rh=-rr, th
     Z_RV_CORR2, Z_RC_CORR2, Z_RR_CORR2, Z_RI_CORR2,        &
     Z_CC_CORR2, Z_CR_CORR2, Z_CI_CORR2
!
! for the conversion from rain to cloud, we need a 3D variable instead of a 1D packed variable
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) ::  &
     Z_RR_CVRC, Z_CR_CVRC                                     ! conversion of rain into cloud droplets (CVRC)

!
! Packed variables for total tendencies
REAL, DIMENSION(:),   ALLOCATABLE ::                                              &
     ZA_TH, ZA_RV, ZA_RC, ZA_CC, ZA_RR, ZA_CR, ZA_RI, ZA_CI, ZA_RS, ZA_RG, ZA_RH, & ! ZA = continuous tendencies (kg/kg/s = S variable)
     ZB_TH, ZB_RV, ZB_RC, ZB_CC, ZB_RR, ZB_CR, ZB_RI, ZB_CI, ZB_RS, ZB_RG, ZB_RH    ! ZB = instant mixing ratio change (kg/kg = T variable)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZB_IFNN

!
! for each process & species, we need 3D variables to store total mmr and conc change (kg/kg and #/kg and theta)
REAL, DIMENSION(:,:,:), ALLOCATABLE ::                                     &
! instantaneous processes :
     ZTOT_CR_BRKU,                                                         & ! spontaneous break up of drops (BRKU)
     ZTOT_TH_HONR, ZTOT_RR_HONR, ZTOT_CR_HONR,                             & ! rain drops homogeneous freezing (HONR)
     ZTOT_TH_IMLT, ZTOT_RC_IMLT, ZTOT_CC_IMLT,                             & ! ice melting (IMLT)
! continuous processes :
     ZTOT_TH_HONC, ZTOT_RC_HONC, ZTOT_CC_HONC,                             & ! droplets homogeneous freezing (HONC)
     ZTOT_CC_SELF,                                                         & ! self collection of droplets (SELF)
     ZTOT_RC_AUTO, ZTOT_CC_AUTO, ZTOT_CR_AUTO,                             & ! autoconversion of cloud droplets (AUTO)
     ZTOT_RC_ACCR, ZTOT_CC_ACCR,                                           & ! accretion of droplets by rain drops (ACCR)
     ZTOT_CR_SCBU,                                                         & ! self collectio break up of drops (SCBU)
!      ZTOT_TH_EVAP, ZTOT_RC_EVAP, ZTOT_CC_EVAP, ZTOT_RR_EVAP, ZTOT_CR_EVAP, & ! evaporation of rain drops (EVAP)
     ZTOT_TH_EVAP, ZTOT_RR_EVAP, & ! evaporation of rain drops (EVAP)
     ZTOT_RI_CNVI, ZTOT_CI_CNVI,                                           & ! conversion snow -> ice (CNVI)
     ZTOT_TH_DEPS, ZTOT_RS_DEPS,                                           & ! deposition of vapor on snow (DEPS)
     ZTOT_TH_DEPI, ZTOT_RI_DEPI,                                           & ! deposition of vapor on ice (DEPI)
     ZTOT_RI_CNVS, ZTOT_CI_CNVS,                                           & ! conversion ice -> snow (CNVS)
     ZTOT_RI_AGGS, ZTOT_CI_AGGS,                                           & ! aggregation of ice on snow (AGGS)
     ZTOT_TH_DEPG, ZTOT_RG_DEPG,                                           & ! deposition of vapor on graupel (DEPG)
     ZTOT_TH_BERFI, ZTOT_RC_BERFI,                                         & ! Bergeron (BERFI)
     ZTOT_TH_RIM, ZTOT_RC_RIM, ZTOT_CC_RIM, ZTOT_RS_RIM, ZTOT_RG_RIM,      & ! cloud droplet riming (RIM)
     ZTOT_RI_HMS, ZTOT_CI_HMS, ZTOT_RS_HMS,                                & ! hallett mossop snow (HMS)
     ZTOT_TH_ACC, ZTOT_RR_ACC, ZTOT_CR_ACC, ZTOT_RS_ACC, ZTOT_RG_ACC,      & ! rain accretion on aggregates (ACC)
     ZTOT_RS_CMEL,                                                         & ! conversion-melting (CMEL)
     ZTOT_TH_CFRZ, ZTOT_RR_CFRZ, ZTOT_CR_CFRZ, ZTOT_RI_CFRZ, ZTOT_CI_CFRZ, & ! rain freezing (CFRZ)
     ZTOT_TH_WETG, ZTOT_RC_WETG, ZTOT_CC_WETG, ZTOT_RR_WETG, ZTOT_CR_WETG, & ! wet growth of graupel (WETG)
     ZTOT_RI_WETG, ZTOT_CI_WETG, ZTOT_RS_WETG, ZTOT_RG_WETG, ZTOT_RH_WETG, & ! wet growth of graupel (WETG)
     ZTOT_TH_DRYG, ZTOT_RC_DRYG, ZTOT_CC_DRYG, ZTOT_RR_DRYG, ZTOT_CR_DRYG, & ! dry growth of graupel (DRYG)
     ZTOT_RI_DRYG, ZTOT_CI_DRYG, ZTOT_RS_DRYG, ZTOT_RG_DRYG,               & ! dry growth of graupel (DRYG)
     ZTOT_RI_HMG, ZTOT_CI_HMG, ZTOT_RG_HMG,                                & ! hallett mossop graupel (HMG)
     ZTOT_TH_GMLT, ZTOT_RR_GMLT, ZTOT_CR_GMLT,                             & ! graupel melting (GMLT)
!      ZTOT_RC_WETH, ZTOT_CC_WETH, ZTOT_RR_WETH, ZTOT_CR_WETH,               & ! wet growth of hail (WETH)
!      ZTOT_RI_WETH, ZTOT_CI_WETH, ZTOT_RS_WETH, ZTOT_RG_WETH, ZTOT_RH_WETH, & ! wet growth of hail (WETH)
!      ZTOT_RG_COHG,                                                         & ! conversion of hail into graupel (COHG)
!      ZTOT_RR_HMLT, ZTOT_CR_HMLT,                                           & ! hail melting (HMLT)
     ZTOT_RR_CVRC, ZTOT_CR_CVRC,                                           & ! conversion of rain into cloud droplets if diameter too small
     ZTOT_RV_CORR2, ZTOT_RC_CORR2, ZTOT_RR_CORR2, ZTOT_RI_CORR2,           &
     ZTOT_CC_CORR2, ZTOT_CR_CORR2, ZTOT_CI_CORR2
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZTOT_IFNN_IMLT

!
!For mixing-ratio splitting
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)) :: Z0RVT,   Z0RCT,   Z0RRT,   Z0RIT,   Z0RST,   Z0RGT,   Z0RHT
REAL, DIMENSION(:), ALLOCATABLE                      :: Z0RVT1D, Z0RCT1D, Z0RRT1D, Z0RIT1D, Z0RST1D, Z0RGT1D, Z0RHT1D 

!
! Loop control variables
REAL,    DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)) :: ZTIME,   ZTIME_LASTCALL,   IITER
REAL,    DIMENSION(:), ALLOCATABLE                      :: ZTIME1D, ZTIME_LASTCALL1D, IITER1D, ZMAXTIME, ZTIME_THRESHOLD
LOGICAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)) :: LLCOMPUTE
LOGICAL, DIMENSION(:), ALLOCATABLE                      :: LLCOMPUTE1D
REAL                                                    :: ZTSTEP
INTEGER                                                 :: INB_ITER_MAX
!
!For subgrid clouds
REAL, DIMENSION(:), ALLOCATABLE                      :: ZCF1D, ZIF1D, ZPF1D     ! 1D packed cloud, ice and precip. frac.

!
! Various parameters
! domain size and levels (AROME compatibility)
INTEGER :: KRR
INTEGER :: IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKT, IKTB, IKTE
! loops and packing
INTEGER :: II, IPACK, JI, JJ, JK
integer :: idx
INTEGER, DIMENSION(:), ALLOCATABLE :: I1, I2, I3
! Inverse ov PTSTEP
REAL :: ZINV_TSTEP
! Work arrays
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)) :: ZW3D
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2))             :: ZW2D
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)) :: ZRT_SUM ! Total condensed water mr
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)) :: ZCPT    ! Total condensed water mr
LOGICAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2))          :: GDEP
real, dimension(:,:,:), allocatable :: zrhodjontstep
!
!-------------------------------------------------------------------------------
!
!*       0.     Init
!               ----
!
!
IIB=1+JPHEXT      ! first physical point in i
IIT=SIZE(PDZZ,1)  ! total number of points in i
IIE=IIT - JPHEXT  ! last physical point in i
!
IJB=1+JPHEXT      ! first physical point in j
IJT=SIZE(PDZZ,2)  ! total number of points in j
IJE=IJT - JPHEXT  ! last physical point in j
!
IKB=KKA+JPVEXT*KKL ! near ground physical point
IKE=KKU-JPVEXT*KKL ! near TOA physical point
IKT=SIZE(PDZZ,3)   ! total number of points in k
!
IKTB=1+JPVEXT      ! first index for a physical point in k
IKTE=IKT-JPVEXT    ! last index for a physical point in k
!
ZTHS(:,:,:) = PTHS(:,:,:)
ZTHT(:,:,:) = PTHS(:,:,:) * PTSTEP
ZRVT(:,:,:) = 0.
ZRVS(:,:,:) = 0.
ZRCT(:,:,:) = 0.
ZRCS(:,:,:) = 0.
ZRRT(:,:,:) = 0.
ZRRS(:,:,:) = 0.
ZRIT(:,:,:) = 0.
ZRIS(:,:,:) = 0.
ZRST(:,:,:) = 0.
ZRSS(:,:,:) = 0.
ZRGT(:,:,:) = 0.
ZRGS(:,:,:) = 0.
ZRHT(:,:,:) = 0.
ZRHS(:,:,:) = 0.
ZRT_SUM(:,:,:) = 0.
ZCCT(:,:,:)   = 0.
ZCCS(:,:,:)   = 0.
ZCRT(:,:,:)   = 0.
ZCRS(:,:,:)   = 0.
ZCIT(:,:,:)   = 0.
ZCIS(:,:,:)   = 0.
ZCCNFT(:,:,:,:) = 0.
ZCCNAT(:,:,:,:) = 0.
ZCCNFS(:,:,:,:) = 0.
ZCCNAS(:,:,:,:) = 0.
ZIFNFT(:,:,:,:) = 0.
ZIFNNT(:,:,:,:) = 0.
ZIFNFS(:,:,:,:) = 0.
ZIFNNS(:,:,:,:) = 0.
ZIMMNT(:,:,:,:) = 0.
ZIMMNS(:,:,:,:) = 0.
ZHOMFT(:,:,:)   = 0.
ZHOMFS(:,:,:)   = 0.

if ( lbu_enable ) then
  allocate( ZTOT_CR_BRKU (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_BRKU(:,:,:) = 0.
  allocate( ZTOT_TH_HONR (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_HONR(:,:,:) = 0.
  allocate( ZTOT_RR_HONR (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RR_HONR(:,:,:) = 0.
  allocate( ZTOT_CR_HONR (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_HONR(:,:,:) = 0.
  allocate( ZTOT_TH_IMLT (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_IMLT(:,:,:) = 0.
  allocate( ZTOT_RC_IMLT (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RC_IMLT(:,:,:) = 0.
  allocate( ZTOT_CC_IMLT (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CC_IMLT(:,:,:) = 0.
  allocate( ZTOT_IFNN_IMLT (size( ptht, 1), size( ptht, 2), size( ptht, 3), nmod_ifn ) ); ZTOT_IFNN_IMLT(:,:,:,:) = 0.
  allocate( ZTOT_TH_HONC (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_HONC(:,:,:) = 0.
  allocate( ZTOT_RC_HONC (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RC_HONC(:,:,:) = 0.
  allocate( ZTOT_CC_HONC (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CC_HONC(:,:,:) = 0.
  allocate( ZTOT_CC_SELF (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CC_SELF(:,:,:) = 0.
  allocate( ZTOT_RC_AUTO (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RC_AUTO(:,:,:) = 0.
  allocate( ZTOT_CC_AUTO (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CC_AUTO(:,:,:) = 0.
  allocate( ZTOT_CR_AUTO (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_AUTO(:,:,:) = 0.
  allocate( ZTOT_RC_ACCR (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RC_ACCR(:,:,:) = 0.
  allocate( ZTOT_CC_ACCR (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CC_ACCR(:,:,:) = 0.
  allocate( ZTOT_CR_SCBU (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_SCBU(:,:,:) = 0.
  allocate( ZTOT_TH_EVAP (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_EVAP(:,:,:) = 0.
!   allocate( ZTOT_RC_EVAP (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RC_EVAP(:,:,:) = 0.
!   allocate( ZTOT_CC_EVAP (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CC_EVAP(:,:,:) = 0.
  allocate( ZTOT_RR_EVAP (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RR_EVAP(:,:,:) = 0.
!   allocate( ZTOT_CR_EVAP (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_EVAP(:,:,:) = 0.
  allocate( ZTOT_RI_CNVI (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RI_CNVI(:,:,:) = 0.
  allocate( ZTOT_CI_CNVI (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CI_CNVI(:,:,:) = 0.
  allocate( ZTOT_TH_DEPS (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_DEPS(:,:,:) = 0.
  allocate( ZTOT_RS_DEPS (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RS_DEPS(:,:,:) = 0.
  allocate( ZTOT_TH_DEPI (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_DEPI(:,:,:) = 0.
  allocate( ZTOT_RI_DEPI (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RI_DEPI(:,:,:) = 0.
  allocate( ZTOT_RI_CNVS (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RI_CNVS(:,:,:) = 0.
  allocate( ZTOT_CI_CNVS (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CI_CNVS(:,:,:) = 0.
  allocate( ZTOT_RI_AGGS (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RI_AGGS(:,:,:) = 0.
  allocate( ZTOT_CI_AGGS (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CI_AGGS(:,:,:) = 0.
  allocate( ZTOT_TH_DEPG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_DEPG(:,:,:) = 0.
  allocate( ZTOT_RG_DEPG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RG_DEPG(:,:,:) = 0.
  allocate( ZTOT_TH_BERFI(size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_BERFI(:,:,:) = 0.
  allocate( ZTOT_RC_BERFI(size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RC_BERFI(:,:,:) = 0.
  allocate( ZTOT_TH_RIM  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_RIM(:,:,:) = 0.
  allocate( ZTOT_RC_RIM  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RC_RIM(:,:,:) = 0.
  allocate( ZTOT_CC_RIM  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CC_RIM(:,:,:) = 0.
  allocate( ZTOT_RS_RIM  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RS_RIM(:,:,:) = 0.
  allocate( ZTOT_RG_RIM  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RG_RIM(:,:,:) = 0.
  allocate( ZTOT_RI_HMS  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RI_HMS(:,:,:) = 0.
  allocate( ZTOT_CI_HMS  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CI_HMS(:,:,:) = 0.
  allocate( ZTOT_RS_HMS  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RS_HMS(:,:,:) = 0.
  allocate( ZTOT_TH_ACC  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_ACC(:,:,:) = 0.
  allocate( ZTOT_RR_ACC  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RR_ACC(:,:,:) = 0.
  allocate( ZTOT_CR_ACC  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_ACC(:,:,:) = 0.
  allocate( ZTOT_RS_ACC  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RS_ACC(:,:,:) = 0.
  allocate( ZTOT_RG_ACC  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RG_ACC(:,:,:) = 0.
  allocate( ZTOT_RS_CMEL (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RS_CMEL(:,:,:) = 0.
  allocate( ZTOT_TH_CFRZ (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_CFRZ(:,:,:) = 0.
  allocate( ZTOT_RR_CFRZ (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RR_CFRZ(:,:,:) = 0.
  allocate( ZTOT_CR_CFRZ (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_CFRZ(:,:,:) = 0.
  allocate( ZTOT_RI_CFRZ (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RI_CFRZ(:,:,:) = 0.
  allocate( ZTOT_CI_CFRZ (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CI_CFRZ(:,:,:) = 0.
  allocate( ZTOT_TH_WETG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_WETG(:,:,:) = 0.
  allocate( ZTOT_RC_WETG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RC_WETG(:,:,:) = 0.
  allocate( ZTOT_CC_WETG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CC_WETG(:,:,:) = 0.
  allocate( ZTOT_RR_WETG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RR_WETG(:,:,:) = 0.
  allocate( ZTOT_CR_WETG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_WETG(:,:,:) = 0.
  allocate( ZTOT_RI_WETG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RI_WETG(:,:,:) = 0.
  allocate( ZTOT_CI_WETG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CI_WETG(:,:,:) = 0.
  allocate( ZTOT_RS_WETG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RS_WETG(:,:,:) = 0.
  allocate( ZTOT_RG_WETG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RG_WETG(:,:,:) = 0.
  allocate( ZTOT_RH_WETG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RH_WETG(:,:,:) = 0.
  allocate( ZTOT_TH_DRYG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_DRYG(:,:,:) = 0.
  allocate( ZTOT_RC_DRYG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RC_DRYG(:,:,:) = 0.
  allocate( ZTOT_CC_DRYG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CC_DRYG(:,:,:) = 0.
  allocate( ZTOT_RR_DRYG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RR_DRYG(:,:,:) = 0.
  allocate( ZTOT_CR_DRYG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_DRYG(:,:,:) = 0.
  allocate( ZTOT_RI_DRYG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RI_DRYG(:,:,:) = 0.
  allocate( ZTOT_CI_DRYG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CI_DRYG(:,:,:) = 0.
  allocate( ZTOT_RS_DRYG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RS_DRYG(:,:,:) = 0.
  allocate( ZTOT_RG_DRYG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RG_DRYG(:,:,:) = 0.
  allocate( ZTOT_RI_HMG  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RI_HMG(:,:,:) = 0.
  allocate( ZTOT_CI_HMG  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CI_HMG(:,:,:) = 0.
  allocate( ZTOT_RG_HMG  (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RG_HMG(:,:,:) = 0.
  allocate( ZTOT_TH_GMLT (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_TH_GMLT(:,:,:) = 0.
  allocate( ZTOT_RR_GMLT (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RR_GMLT(:,:,:) = 0.
  allocate( ZTOT_CR_GMLT (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_GMLT(:,:,:) = 0.
!   allocate( ZTOT_RC_WETH (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RC_WETH(:,:,:) = 0.
!   allocate( ZTOT_CC_WETH (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CC_WETH(:,:,:) = 0.
!   allocate( ZTOT_RR_WETH (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RR_WETH(:,:,:) = 0.
!   allocate( ZTOT_CR_WETH (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_WETH(:,:,:) = 0.
!   allocate( ZTOT_RI_WETH (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RI_WETH(:,:,:) = 0.
!   allocate( ZTOT_CI_WETH (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CI_WETH(:,:,:) = 0.
!   allocate( ZTOT_RS_WETH (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RS_WETH(:,:,:) = 0.
!   allocate( ZTOT_RG_WETH (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RG_WETH(:,:,:) = 0.
!   allocate( ZTOT_RH_WETH (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RH_WETH(:,:,:) = 0.
!   allocate( ZTOT_RG_COHG (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RG_COHG(:,:,:) = 0.
!   allocate( ZTOT_RR_HMLT (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RR_HMLT(:,:,:) = 0.
!   allocate( ZTOT_CR_HMLT (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_HMLT(:,:,:) = 0.
  allocate( ZTOT_RR_CVRC (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RR_CVRC(:,:,:) = 0.
  allocate( ZTOT_CR_CVRC (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_CVRC(:,:,:) = 0.

  allocate( ZTOT_RV_CORR2 (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RV_CORR2(:,:,:) = 0.
  allocate( ZTOT_RC_CORR2 (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RC_CORR2(:,:,:) = 0.
  allocate( ZTOT_RR_CORR2 (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RR_CORR2(:,:,:) = 0.
  allocate( ZTOT_RI_CORR2 (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_RI_CORR2(:,:,:) = 0.
  allocate( ZTOT_CC_CORR2 (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CC_CORR2(:,:,:) = 0.
  allocate( ZTOT_CR_CORR2 (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CR_CORR2(:,:,:) = 0.
  allocate( ZTOT_CI_CORR2 (size( ptht, 1), size( ptht, 2), size( ptht, 3) ) ); ZTOT_CI_CORR2(:,:,:) = 0.
END IF
!
! Initial values computed as source * PTSTEP
!
! Mixing ratios
!
KRR=SIZE(PRT,4)
ZRVT(:,:,:) = PRS(:,:,:,1) * PTSTEP
ZRVS(:,:,:) = PRS(:,:,:,1)
IF ( KRR .GE. 2 ) ZRCT(:,:,:) = PRS(:,:,:,2) * PTSTEP
IF ( KRR .GE. 2 ) ZRCS(:,:,:) = PRS(:,:,:,2)
IF ( KRR .GE. 3 ) ZRRT(:,:,:) = PRS(:,:,:,3) * PTSTEP
IF ( KRR .GE. 3 ) ZRRS(:,:,:) = PRS(:,:,:,3)
IF ( KRR .GE. 4 ) ZRIT(:,:,:) = PRS(:,:,:,4) * PTSTEP
IF ( KRR .GE. 4 ) ZRIS(:,:,:) = PRS(:,:,:,4)
IF ( KRR .GE. 5 ) ZRST(:,:,:) = PRS(:,:,:,5) * PTSTEP
IF ( KRR .GE. 5 ) ZRSS(:,:,:) = PRS(:,:,:,5)
IF ( KRR .GE. 6 ) ZRGT(:,:,:) = PRS(:,:,:,6) * PTSTEP
IF ( KRR .GE. 6 ) ZRGS(:,:,:) = PRS(:,:,:,6)
IF ( KRR .GE. 7 ) ZRHT(:,:,:) = PRS(:,:,:,7) * PTSTEP
IF ( KRR .GE. 7 ) ZRHS(:,:,:) = PRS(:,:,:,7)
!
! Concentrations
!
IF ( LWARM )             ZCCT(:,:,:)   = PSVS(:,:,:,NSV_LIMA_NC) * PTSTEP
IF ( LWARM )             ZCCS(:,:,:)   = PSVS(:,:,:,NSV_LIMA_NC)
IF ( LWARM .AND. LRAIN ) ZCRT(:,:,:)   = PSVS(:,:,:,NSV_LIMA_NR) * PTSTEP
IF ( LWARM .AND. LRAIN ) ZCRS(:,:,:)   = PSVS(:,:,:,NSV_LIMA_NR)
IF ( LCOLD )             ZCIT(:,:,:)   = PSVS(:,:,:,NSV_LIMA_NI) * PTSTEP
IF ( LCOLD )             ZCIS(:,:,:)   = PSVS(:,:,:,NSV_LIMA_NI)
!
IF ( NMOD_CCN .GE. 1 )   ZCCNFT(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1) * PTSTEP
IF ( NMOD_CCN .GE. 1 )   ZCCNAT(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1) * PTSTEP
IF ( NMOD_CCN .GE. 1 )   ZCCNFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1)
IF ( NMOD_CCN .GE. 1 )   ZCCNAS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1)
!
IF ( NMOD_IFN .GE. 1 )   ZIFNFT(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1) * PTSTEP
IF ( NMOD_IFN .GE. 1 )   ZIFNNT(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1) * PTSTEP
IF ( NMOD_IFN .GE. 1 )   ZIFNFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1)
IF ( NMOD_IFN .GE. 1 )   ZIFNNS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1)
!
IF ( NMOD_IMM .GE. 1 )   ZIMMNT(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1) * PTSTEP
IF ( NMOD_IMM .GE. 1 )   ZIMMNS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1)
!
IF ( LCOLD .AND. LHHONI ) ZHOMFT(:,:,:)  = PSVS(:,:,:,NSV_LIMA_HOM_HAZE) * PTSTEP
IF ( LCOLD .AND. LHHONI ) ZHOMFS(:,:,:)  = PSVS(:,:,:,NSV_LIMA_HOM_HAZE)
!
ZINV_TSTEP  = 1./PTSTEP
ZEXN(:,:,:) = (PPABST(:,:,:)/XP00)**(XRD/XCPD)
ZT(:,:,:)   = ZTHT(:,:,:) * ZEXN(:,:,:)
!
!-------------------------------------------------------------------------------
!
!*       0.     Check mean diameter for cloud, rain  and ice
!               --------------------------------------------
! if ( lbu_enable ) then
!   if ( lbudget_rc .and. lwarm .and. lrain ) call Budget_store_init( tbudgets(NBUDGET_RC), 'CORR', zrcs(:, :, :) * prhodj(:, :, :) )
!   if ( lbudget_rr .and. lwarm .and. lrain ) call Budget_store_init( tbudgets(NBUDGET_RR), 'CORR', zrrs(:, :, :) * prhodj(:, :, :) )
!   if ( lbudget_ri .and. lcold .and. lsnow ) call Budget_store_init( tbudgets(NBUDGET_RI), 'CORR', zris(:, :, :) * prhodj(:, :, :) )
!   if ( lbudget_rs .and. lcold .and. lsnow ) call Budget_store_init( tbudgets(NBUDGET_RS), 'CORR', zrss(:, :, :) * prhodj(:, :, :) )
!   if ( lbudget_sv ) then
!     if ( lwarm .and. lrain ) &
!       call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'CORR', zccs(:, :, :) * prhodj(:, :, :) )
!     if ( lwarm .and. lrain ) &
!       call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'CORR', zcrs(:, :, :) * prhodj(:, :, :) )
!     if ( lcold .and. lsnow ) &
!       call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CORR', zcis(:, :, :) * prhodj(:, :, :) )
!   end if
! end if
!!$IF (LWARM .AND. LRAIN) THEN
!!$   WHERE( ZRCT>XRTMIN(2) .AND. ZCCT>XCTMIN(2) .AND. ZRCT>XAC*ZCCT*(100.E-6)**XBC )
!!$      ZRRT=ZRRT+ZRCT
!!$      ZRRS=ZRRS+ZRCS
!!$      ZCRT=ZCRT+ZCCT
!!$      ZCRS=ZCRS+ZCCS
!!$      ZRCT=0.
!!$      ZCCT=0.
!!$      ZRCS=0.
!!$      ZCCS=0.
!!$   END WHERE
!!$END IF
!!$!
!!$IF (LWARM .AND. LRAIN) THEN
!!$   WHERE( ZRRT>XRTMIN(3) .AND. ZCRT>XCTMIN(3) .AND. ZRRT<XAR*ZCRT*(60.E-6)**XBR )
!!$      ZRCT=ZRCT+ZRRT
!!$      ZRCS=ZRCS+ZRRS
!!$      ZCCT=ZCCT+ZCRT
!!$      ZCCS=ZCCS+ZCRS
!!$      ZRRT=0.
!!$      ZCRT=0.
!!$      ZRRS=0.
!!$      ZCRS=0.
!!$   END WHERE
!!$END IF
!!$!
!!$IF (LCOLD .AND. LSNOW) THEN
!!$   WHERE( ZRIT>XRTMIN(4) .AND. ZCIT>XCTMIN(4) .AND. ZRIT>XAI*ZCIT*(250.E-6)**XBI )
!!$      ZRST=ZRST+ZRIT
!!$      ZRSS=ZRSS+ZRIS
!!$      ZRIT=0.
!!$      ZCIT=0.
!!$      ZRIS=0.
!!$      ZCIS=0.
!!$   END WHERE
!!$END IF
!
! if ( lbu_enable ) then
!   if ( lbudget_rc .and. lwarm .and. lrain ) call Budget_store_end( tbudgets(NBUDGET_RC), 'CORR', zrcs(:, :, :) * prhodj(:, :, :) )
!   if ( lbudget_rr .and. lwarm .and. lrain ) call Budget_store_end( tbudgets(NBUDGET_RR), 'CORR', zrrs(:, :, :) * prhodj(:, :, :) )
!   if ( lbudget_ri .and. lcold .and. lsnow ) call Budget_store_end( tbudgets(NBUDGET_RI), 'CORR', zris(:, :, :) * prhodj(:, :, :) )
!   if ( lbudget_rs .and. lcold .and. lsnow ) call Budget_store_end( tbudgets(NBUDGET_RS), 'CORR', zrss(:, :, :) * prhodj(:, :, :) )
!   if ( lbudget_sv ) then
!     if ( lwarm .and. lrain ) &
!       call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'CORR', zccs(:, :, :) * prhodj(:, :, :) )
!     if ( lwarm .and. lrain ) &
!       call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'CORR', zcrs(:, :, :) * prhodj(:, :, :) )
!     if ( lcold .and. lsnow ) &
!       call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CORR', zcis(:, :, :) * prhodj(:, :, :) )
!   end if
! end if
!-------------------------------------------------------------------------------
!
!*       1.     Sedimentation
!               -------------
!
!
if ( lbu_enable ) then
  if ( lbudget_th )                         call Budget_store_init( tbudgets(NBUDGET_TH), 'SEDI', zths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc .and. lwarm .and. lsedc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'SEDI', zrcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr .and. lwarm .and. lrain ) call Budget_store_init( tbudgets(NBUDGET_RR), 'SEDI', zrrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri .and. lcold .and. lsedi ) call Budget_store_init( tbudgets(NBUDGET_RI), 'SEDI', zris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs .and. lcold .and. lsnow ) call Budget_store_init( tbudgets(NBUDGET_RS), 'SEDI', zrss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg .and. lcold .and. lsnow ) call Budget_store_init( tbudgets(NBUDGET_RG), 'SEDI', zrgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh .and. lcold .and. lhail ) call Budget_store_init( tbudgets(NBUDGET_RH), 'SEDI', zrhs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if ( lwarm .and. lsedc ) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'SEDI', zccs(:, :, :) * prhodj(:, :, :) )
    if ( lwarm .and. lrain ) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'SEDI', zcrs(:, :, :) * prhodj(:, :, :) )
    if ( lcold .and. lsedi ) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'SEDI', zcis(:, :, :) * prhodj(:, :, :) )
  end if
end if

ZRT_SUM = (ZRVS + ZRCS + ZRRS + ZRIS + ZRSS + ZRGS + ZRHS)*PTSTEP
ZCPT    = XCPD + (XCPV * ZRVS + XCL * (ZRCS + ZRRS) + XCI * (ZRIS + ZRSS + ZRGS + ZRHS))*PTSTEP
IF (LWARM .AND. LSEDC) CALL LIMA_SEDIMENTATION(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                               'L', 2, 2, 1, PTSTEP, PDZZ, PRHODREF, PPABST, ZT, ZRT_SUM, ZCPT, ZRCS, ZCCS, PINPRC)
ZRT_SUM = (ZRVS + ZRCS + ZRRS + ZRIS + ZRSS + ZRGS + ZRHS)*PTSTEP
ZCPT    = XCPD + (XCPV * ZRVS + XCL * (ZRCS + ZRRS) + XCI * (ZRIS + ZRSS + ZRGS + ZRHS))*PTSTEP
IF (LWARM .AND. LRAIN) CALL LIMA_SEDIMENTATION(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                               'L', 2, 3, 1, PTSTEP, PDZZ, PRHODREF, PPABST, ZT, ZRT_SUM, ZCPT, ZRRS, ZCRS, PINPRR)
ZRT_SUM = (ZRVS + ZRCS + ZRRS + ZRIS + ZRSS + ZRGS + ZRHS)*PTSTEP
ZCPT    = XCPD + (XCPV * ZRVS + XCL * (ZRCS + ZRRS) + XCI * (ZRIS + ZRSS + ZRGS + ZRHS))*PTSTEP
IF (LCOLD .AND. LSEDI) CALL LIMA_SEDIMENTATION(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                               'I', 2, 4, 1, PTSTEP, PDZZ, PRHODREF, PPABST, ZT, ZRT_SUM, ZCPT, ZRIS, ZCIS, ZW2D)
ZRT_SUM = (ZRVS + ZRCS + ZRRS + ZRIS + ZRSS + ZRGS + ZRHS)*PTSTEP
ZCPT    = XCPD + (XCPV * ZRVS + XCL * (ZRCS + ZRRS) + XCI * (ZRIS + ZRSS + ZRGS + ZRHS))*PTSTEP
IF (LCOLD .AND. LSNOW) CALL LIMA_SEDIMENTATION(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                               'I', 1, 5, 1, PTSTEP, PDZZ, PRHODREF, PPABST, ZT, ZRT_SUM, ZCPT, ZRSS, ZW3D, PINPRS)
ZRT_SUM = (ZRVS + ZRCS + ZRRS + ZRIS + ZRSS + ZRGS + ZRHS)*PTSTEP
ZCPT    = XCPD + (XCPV * ZRVS + XCL * (ZRCS + ZRRS) + XCI * (ZRIS + ZRSS + ZRGS + ZRHS))*PTSTEP
IF (LCOLD .AND. LSNOW) CALL LIMA_SEDIMENTATION(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                               'I', 1, 6, 1, PTSTEP, PDZZ, PRHODREF, PPABST, ZT, ZRT_SUM, ZCPT, ZRGS, ZW3D, PINPRG)
ZRT_SUM = (ZRVS + ZRCS + ZRRS + ZRIS + ZRSS + ZRGS + ZRHS)*PTSTEP
ZCPT    = XCPD + (XCPV * ZRVS + XCL * (ZRCS + ZRRS) + XCI * (ZRIS + ZRSS + ZRGS + ZRHS))*PTSTEP
IF (LCOLD .AND. LHAIL) CALL LIMA_SEDIMENTATION(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                               'I', 1, 7, 1, PTSTEP, PDZZ, PRHODREF, PPABST, ZT, ZRT_SUM, ZCPT, ZRHS, ZW3D, PINPRH)
!
ZTHS(:,:,:) = ZT(:,:,:) / ZEXN(:,:,:) * ZINV_TSTEP
!
! Call budgets
!
if ( lbu_enable ) then
  if ( lbudget_th )                         call Budget_store_end( tbudgets(NBUDGET_TH), 'SEDI', zths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc .and. lwarm .and. lsedc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'SEDI', zrcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr .and. lwarm .and. lrain ) call Budget_store_end( tbudgets(NBUDGET_RR), 'SEDI', zrrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri .and. lcold .and. lsedi ) call Budget_store_end( tbudgets(NBUDGET_RI), 'SEDI', zris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs .and. lcold .and. lsnow ) call Budget_store_end( tbudgets(NBUDGET_RS), 'SEDI', zrss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg .and. lcold .and. lsnow ) call Budget_store_end( tbudgets(NBUDGET_RG), 'SEDI', zrgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh .and. lcold .and. lhail ) call Budget_store_end( tbudgets(NBUDGET_RH), 'SEDI', zrhs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if ( lwarm .and. lsedc ) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'SEDI', zccs(:, :, :) * prhodj(:, :, :) )
    if ( lwarm .and. lrain ) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'SEDI', zcrs(:, :, :) * prhodj(:, :, :) )
    if ( lcold .and. lsedi ) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'SEDI', zcis(:, :, :) * prhodj(:, :, :) )
  end if
end if
!
! 1.bis Deposition at 1st level above ground
!
IF (LWARM .AND. LDEPOC) THEN
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC),                    'DEPO', zrcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'DEPO', zccs(:, :, :) * prhodj(:, :, :) )

  PINDEP(:,:)=0.
  GDEP(:,:) = .FALSE.
  GDEP(:,:) = ZRCS(:,:,IKB) >0 .AND. ZCCS(:,:,IKB) >0 .AND. ZRCT(:,:,IKB) >0 .AND. ZCCT(:,:,IKB) >0
  WHERE (GDEP)
     ZRCS(:,:,IKB) = ZRCS(:,:,IKB) - XVDEPOC * ZRCT(:,:,IKB) / PDZZ(:,:,IKB)
     ZCCS(:,:,IKB) = ZCCS(:,:,IKB) - XVDEPOC * ZCCT(:,:,IKB) / PDZZ(:,:,IKB)
     PINPRC(:,:) = PINPRC(:,:) + XVDEPOC * ZRCT(:,:,IKB) * PRHODREF(:,:,IKB) /XRHOLW
     PINDEP(:,:) = XVDEPOC * ZRCT(:,:,IKB) *  PRHODREF(:,:,IKB) /XRHOLW
  END WHERE

  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC),                    'DEPO', zrcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'DEPO', zccs(:, :, :) * prhodj(:, :, :) )
END IF
!
!
Z_RR_CVRC(:,:,:) = 0.
Z_CR_CVRC(:,:,:) = 0.
IF (LWARM .AND. LRAIN) THEN
   if( lbu_enable ) then
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC),                    'R2C1', zrcs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR),                    'R2C1', zrrs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'R2C1', zccs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'R2C1', zcrs(:, :, :) * prhodj(:, :, :) )
   end if

   CALL LIMA_DROPS_TO_DROPLETS_CONV(PRHODREF, ZRCS*PTSTEP, ZRRS*PTSTEP, ZCCS*PTSTEP, ZCRS*PTSTEP, &
                                    Z_RR_CVRC, Z_CR_CVRC)
   !
   ZRCS(:,:,:) = ZRCS(:,:,:) - Z_RR_CVRC(:,:,:)/PTSTEP
   ZRRS(:,:,:) = ZRRS(:,:,:) + Z_RR_CVRC(:,:,:)/PTSTEP
   ZCCS(:,:,:) = ZCCS(:,:,:) - Z_CR_CVRC(:,:,:)/PTSTEP
   ZCRS(:,:,:) = ZCRS(:,:,:) + Z_CR_CVRC(:,:,:)/PTSTEP

   if( lbu_enable ) then
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC),                    'R2C1', zrcs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR),                    'R2C1', zrrs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'R2C1', zccs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'R2C1', zcrs(:, :, :) * prhodj(:, :, :) )
   end if
END IF
!
! Update variables
!
ZTHT(:,:,:) = ZTHS(:,:,:) * PTSTEP
ZT(:,:,:)   = ZTHT(:,:,:) * ZEXN(:,:,:)
!
IF ( KRR .GE. 2 ) ZRCT(:,:,:) = ZRCS(:,:,:) * PTSTEP
IF ( KRR .GE. 3 ) ZRRT(:,:,:) = ZRRS(:,:,:) * PTSTEP
IF ( KRR .GE. 4 ) ZRIT(:,:,:) = ZRIS(:,:,:) * PTSTEP
IF ( KRR .GE. 5 ) ZRST(:,:,:) = ZRSS(:,:,:) * PTSTEP
IF ( KRR .GE. 6 ) ZRGT(:,:,:) = ZRGS(:,:,:) * PTSTEP
IF ( KRR .GE. 7 ) ZRHT(:,:,:) = ZRHS(:,:,:) * PTSTEP
!
IF ( LWARM )             ZCCT(:,:,:)   = ZCCS(:,:,:) * PTSTEP
IF ( LWARM .AND. LRAIN ) ZCRT(:,:,:)   = ZCRS(:,:,:) * PTSTEP
IF ( LCOLD )             ZCIT(:,:,:)   = ZCIS(:,:,:) * PTSTEP
! 
!-------------------------------------------------------------------------------
!
!*       2.     Compute cloud, ice and precipitation fractions
!               ----------------------------------------------
!
IF (LSUBG_COND) THEN
   CALL LIMA_COMPUTE_CLOUD_FRACTIONS (IIB, IIE, IJB, IJE, IKB, IKE, KKL, &
                                      ZCCT, ZRCT,                        &
                                      ZCRT, ZRRT,                        &
                                      ZCIT, ZRIT,                        &
                                      ZRST, ZRGT, ZRHT,                  &
                                      PCLDFR, PICEFR, PPRCFR             )
ELSE
   PCLDFR(:,:,:)=1.
   PICEFR(:,:,:)=1.
   PPRCFR(:,:,:)=1.
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.     Nucleation processes
!               --------------------
!
CALL LIMA_NUCLEATION_PROCS (PTSTEP, TPFILE, PRHODJ,                             &
                            PRHODREF, ZEXN, PPABST, ZT, PDTHRAD, PW_NU,         &
                            ZTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT,           &
                            ZCCT, ZCRT, ZCIT,                                   &
                            ZCCNFT, ZCCNAT, ZIFNFT, ZIFNNT, ZIMMNT, ZHOMFT,     &
                            PCLDFR, PICEFR, PPRCFR                              )
!
! Saving sources before microphysics time-splitting loop
!
ZRVS(:,:,:) = ZRVT(:,:,:) *ZINV_TSTEP
ZRCS(:,:,:) = ZRCT(:,:,:) *ZINV_TSTEP
ZRRS(:,:,:) = ZRRT(:,:,:) *ZINV_TSTEP
ZRIS(:,:,:) = ZRIT(:,:,:) *ZINV_TSTEP
ZRSS(:,:,:) = ZRST(:,:,:) *ZINV_TSTEP
ZRGS(:,:,:) = ZRGT(:,:,:) *ZINV_TSTEP
ZRHS(:,:,:) = ZRHT(:,:,:) *ZINV_TSTEP
!
ZCCS(:,:,:) = ZCCT(:,:,:) *ZINV_TSTEP
ZCRS(:,:,:) = ZCRT(:,:,:) *ZINV_TSTEP
ZCIS(:,:,:) = ZCIT(:,:,:) *ZINV_TSTEP
!
ZCCNFS(:,:,:,:) = ZCCNFT(:,:,:,:) *ZINV_TSTEP
ZCCNAS(:,:,:,:) = ZCCNAT(:,:,:,:) *ZINV_TSTEP
ZIFNFS(:,:,:,:) = ZIFNFT(:,:,:,:) *ZINV_TSTEP
ZIFNNS(:,:,:,:) = ZIFNNT(:,:,:,:) *ZINV_TSTEP
ZIMMNS(:,:,:,:) = ZIMMNT(:,:,:,:) *ZINV_TSTEP
ZHOMFS(:,:,:)   = ZHOMFT(:,:,:)   *ZINV_TSTEP
!
ZTHS(:,:,:) = ZTHT(:,:,:) *ZINV_TSTEP
ZT(:,:,:)   = ZTHT(:,:,:) * ZEXN(:,:,:)
!
!
!-------------------------------------------------------------------------------
!
!*       2.     LOOP
!               ----
!
!
! Maximum number of iterations
INB_ITER_MAX=NMAXITER
IF(XTSTEP_TS/=0.)THEN
  INB_ITER_MAX=MAX(1, INT(PTSTEP/XTSTEP_TS)) !At least the number of iterations needed for the time-splitting
  ZTSTEP=PTSTEP/INB_ITER_MAX
  INB_ITER_MAX=MAX(NMAXITER, INB_ITER_MAX) !Fot the case XMRSTEP/=0. at the same time
ENDIF
IITER(:,:,:)=0
ZTIME(:,:,:)=0. ! Current integration time (all points may have a different integration time)
!
! Begin the huge time splitting loop
!
ZRT_SUM(:,:,:) = ZRCT(:,:,:) + ZRRT(:,:,:) + ZRIT(:,:,:) + ZRST(:,:,:) + ZRGT(:,:,:) + ZRHT(:,:,:)
WHERE (ZRT_SUM(:,:,:)<XRTMIN(2)) ZTIME(:,:,:)=PTSTEP ! no need to treat hydrometeor-free point
!
DO WHILE(ANY(ZTIME(IIB:IIE,IJB:IJE,IKTB:IKTE)<PTSTEP))
   !
   IF(XMRSTEP/=0.) THEN
      ! In this case we need to remember the mixing ratios used to compute the tendencies
      ! because when mixing ratio has evolved more than a threshold, we must re-compute tendecies
      Z0RVT(:,:,:)=ZRVT(:,:,:)
      Z0RCT(:,:,:)=ZRCT(:,:,:)
      Z0RRT(:,:,:)=ZRRT(:,:,:)
      Z0RIT(:,:,:)=ZRIT(:,:,:)
      Z0RST(:,:,:)=ZRST(:,:,:)
      Z0RGT(:,:,:)=ZRGT(:,:,:)
      Z0RHT(:,:,:)=ZRHT(:,:,:)
   ENDIF
   !
   IF(XTSTEP_TS/=0.) THEN
      ! In this case we need to remember the time when tendencies were computed
      ! because when time has evolved more than a limit, we must re-compute tendecies
      ZTIME_LASTCALL(:,:,:)=ZTIME(:,:,:)
   ENDIF
   !
   LLCOMPUTE(:,:,:)=.FALSE.
   LLCOMPUTE(IIB:IIE,IJB:IJE,IKTB:IKTE) = ZTIME(IIB:IIE,IJB:IJE,IKTB:IKTE)<PTSTEP ! Compuation only for points for which integration time has not reached the timestep
   WHERE(LLCOMPUTE(:,:,:))
      IITER(:,:,:)=IITER(:,:,:)+1
   END WHERE
   ! 
   DO WHILE(ANY(LLCOMPUTE(:,:,:))) ! Loop to adjust tendencies when we cross the 0°C or when a species disappears

      !
      ! Packing variables to run computations only where necessary
      !
      IPACK = COUNT(LLCOMPUTE)
      ALLOCATE(I1(IPACK))
      ALLOCATE(I2(IPACK))
      ALLOCATE(I3(IPACK))
      ALLOCATE(ZRHODREF1D(IPACK))
      ALLOCATE(ZEXNREF1D(IPACK))
      ALLOCATE(ZEXN1D(IPACK))
      ALLOCATE(ZP1D(IPACK))     
      ALLOCATE(ZTHT1D(IPACK))
      ALLOCATE(ZRVT1D(IPACK))
      ALLOCATE(ZRCT1D(IPACK))
      ALLOCATE(ZRRT1D(IPACK))
      ALLOCATE(ZRIT1D(IPACK))
      ALLOCATE(ZRST1D(IPACK))
      ALLOCATE(ZRGT1D(IPACK))
      ALLOCATE(ZRHT1D(IPACK))
      ALLOCATE(ZCCT1D(IPACK))
      ALLOCATE(ZCRT1D(IPACK))
      ALLOCATE(ZCIT1D(IPACK))
      ALLOCATE(ZIFNN1D(IPACK,NMOD_IFN))
      ALLOCATE(ZEVAP1D(IPACK))
      ALLOCATE(ZTIME1D(IPACK))
      ALLOCATE(LLCOMPUTE1D(IPACK))
      ALLOCATE(IITER1D(IPACK))
      ALLOCATE(ZTIME_LASTCALL1D(IPACK))
      ALLOCATE(Z0RVT1D(IPACK))
      ALLOCATE(Z0RCT1D(IPACK))
      ALLOCATE(Z0RRT1D(IPACK))
      ALLOCATE(Z0RIT1D(IPACK))
      ALLOCATE(Z0RST1D(IPACK))
      ALLOCATE(Z0RGT1D(IPACK))
      ALLOCATE(Z0RHT1D(IPACK))
      ALLOCATE(ZCF1D(IPACK))
      ALLOCATE(ZIF1D(IPACK))
      ALLOCATE(ZPF1D(IPACK))
      IPACK = COUNTJV(LLCOMPUTE,I1,I2,I3)
      DO II=1,IPACK
         ZRHODREF1D(II)       = PRHODREF(I1(II),I2(II),I3(II))
         ZEXNREF1D(II)        = PEXNREF(I1(II),I2(II),I3(II))
         ZEXN1D(II)           = ZEXN(I1(II),I2(II),I3(II))
         ZP1D(II)             = PPABST(I1(II),I2(II),I3(II))
         ZTHT1D(II)           = ZTHT(I1(II),I2(II),I3(II))
         ZRVT1D(II)           = ZRVT(I1(II),I2(II),I3(II))
         ZRCT1D(II)           = ZRCT(I1(II),I2(II),I3(II))
         ZRRT1D(II)           = ZRRT(I1(II),I2(II),I3(II))
         ZRIT1D(II)           = ZRIT(I1(II),I2(II),I3(II))
         ZRST1D(II)           = ZRST(I1(II),I2(II),I3(II))
         ZRGT1D(II)           = ZRGT(I1(II),I2(II),I3(II))
         ZRHT1D(II)           = ZRHT(I1(II),I2(II),I3(II))
         ZCCT1D(II)           = ZCCT(I1(II),I2(II),I3(II))
         ZCRT1D(II)           = ZCRT(I1(II),I2(II),I3(II))
         ZCIT1D(II)           = ZCIT(I1(II),I2(II),I3(II))
         ZIFNN1D(II,:)        = ZIFNNT(I1(II),I2(II),I3(II),:)
         ZEVAP1D(II)          = PEVAP3D(I1(II),I2(II),I3(II))
         ZTIME1D(II)          = ZTIME(I1(II),I2(II),I3(II))         
         LLCOMPUTE1D(II)      = LLCOMPUTE(I1(II),I2(II),I3(II))         
         IITER1D(II)          = IITER(I1(II),I2(II),I3(II))         
         ZTIME_LASTCALL1D(II) = ZTIME_LASTCALL(I1(II),I2(II),I3(II))         
         Z0RVT1D(II)          = Z0RVT(I1(II),I2(II),I3(II))
         Z0RCT1D(II)          = Z0RCT(I1(II),I2(II),I3(II))
         Z0RRT1D(II)          = Z0RRT(I1(II),I2(II),I3(II))
         Z0RIT1D(II)          = Z0RIT(I1(II),I2(II),I3(II))
         Z0RST1D(II)          = Z0RST(I1(II),I2(II),I3(II))
         Z0RGT1D(II)          = Z0RGT(I1(II),I2(II),I3(II))
         Z0RHT1D(II)          = Z0RHT(I1(II),I2(II),I3(II))
         ZCF1D(II)            = PCLDFR(I1(II),I2(II),I3(II))
         ZIF1D(II)            = PICEFR(I1(II),I2(II),I3(II))
         ZPF1D(II)            = PPRCFR(I1(II),I2(II),I3(II))
      END DO
      !
      WHERE(ZCF1D(:)<1.E-10 .AND. ZRCT1D(:)>XRTMIN(2) .AND. ZCCT1D(:)>XCTMIN(2)) ZCF1D(:)=1.
      WHERE(ZIF1D(:)<1.E-10 .AND. ZRIT1D(:)>XRTMIN(4) .AND. ZCIT1D(:)>XCTMIN(4)) ZIF1D(:)=1.
      WHERE(ZPF1D(:)<1.E-10 .AND. (ZRRT1D(:)>XRTMIN(3) .OR. ZRST1D(:)>XRTMIN(5) &
                              .OR. ZRGT1D(:)>XRTMIN(6) .OR. ZRHT1D(:)>XRTMIN(7) ) ) ZPF1D(:)=1.
      !
      ! Allocating 1D variables
      !
      ALLOCATE(ZMAXTIME(IPACK))           ;  ZMAXTIME(:) = 0.
      ALLOCATE(ZTIME_THRESHOLD(IPACK))    ;  ZTIME_THRESHOLD(:) = 0.
      !
      ALLOCATE(ZA_TH(IPACK))              ;  ZA_TH(:) = 0.
      ALLOCATE(ZA_RV(IPACK))              ;  ZA_RV(:) = 0.
      ALLOCATE(ZA_RC(IPACK))              ;  ZA_RC(:) = 0.
      ALLOCATE(ZA_RR(IPACK))              ;  ZA_RR(:) = 0.
      ALLOCATE(ZA_RI(IPACK))              ;  ZA_RI(:) = 0.
      ALLOCATE(ZA_RS(IPACK))              ;  ZA_RS(:) = 0.
      ALLOCATE(ZA_RG(IPACK))              ;  ZA_RG(:) = 0.
      ALLOCATE(ZA_RH(IPACK))              ;  ZA_RH(:) = 0.
      ALLOCATE(ZA_CC(IPACK))              ;  ZA_CC(:) = 0.
      ALLOCATE(ZA_CR(IPACK))              ;  ZA_CR(:) = 0.
      ALLOCATE(ZA_CI(IPACK))              ;  ZA_CI(:) = 0.
      !
      ALLOCATE(ZB_TH(IPACK))              ;  ZB_TH(:) = 0.
      ALLOCATE(ZB_RV(IPACK))              ;  ZB_RV(:) = 0.
      ALLOCATE(ZB_RC(IPACK))              ;  ZB_RC(:) = 0.
      ALLOCATE(ZB_RR(IPACK))              ;  ZB_RR(:) = 0.
      ALLOCATE(ZB_RI(IPACK))              ;  ZB_RI(:) = 0.
      ALLOCATE(ZB_RS(IPACK))              ;  ZB_RS(:) = 0.
      ALLOCATE(ZB_RG(IPACK))              ;  ZB_RG(:) = 0.
      ALLOCATE(ZB_RH(IPACK))              ;  ZB_RH(:) = 0.
      ALLOCATE(ZB_CC(IPACK))              ;  ZB_CC(:) = 0.
      ALLOCATE(ZB_CR(IPACK))              ;  ZB_CR(:) = 0.
      ALLOCATE(ZB_CI(IPACK))              ;  ZB_CI(:) = 0.
      ALLOCATE(ZB_IFNN(IPACK,NMOD_IFN))   ;  ZB_IFNN(:,:) = 0.
      !
      ALLOCATE(Z_CR_BRKU(IPACK))          ; Z_CR_BRKU(:) = 0.
      ALLOCATE(Z_TH_HONR(IPACK))          ; Z_TH_HONR(:) = 0.
      ALLOCATE(Z_RR_HONR(IPACK))          ; Z_RR_HONR(:) = 0.
      ALLOCATE(Z_CR_HONR(IPACK))          ; Z_CR_HONR(:) = 0.
      ALLOCATE(Z_TH_IMLT(IPACK))          ; Z_TH_IMLT(:) = 0.
      ALLOCATE(Z_RC_IMLT(IPACK))          ; Z_RC_IMLT(:) = 0.
      ALLOCATE(Z_CC_IMLT(IPACK))          ; Z_CC_IMLT(:) = 0.
      ALLOCATE(Z_TH_HONC(IPACK))          ; Z_TH_HONC(:) = 0.
      ALLOCATE(Z_RC_HONC(IPACK))          ; Z_RC_HONC(:) = 0.
      ALLOCATE(Z_CC_HONC(IPACK))          ; Z_CC_HONC(:) = 0.
      ALLOCATE(Z_CC_SELF(IPACK))          ; Z_CC_SELF(:) = 0.
      ALLOCATE(Z_RC_AUTO(IPACK))          ; Z_RC_AUTO(:) = 0.
      ALLOCATE(Z_CC_AUTO(IPACK))          ; Z_CC_AUTO(:) = 0.
      ALLOCATE(Z_CR_AUTO(IPACK))          ; Z_CR_AUTO(:) = 0.
      ALLOCATE(Z_RC_ACCR(IPACK))          ; Z_RC_ACCR(:) = 0.
      ALLOCATE(Z_CC_ACCR(IPACK))          ; Z_CC_ACCR(:) = 0.
      ALLOCATE(Z_CR_SCBU(IPACK))          ; Z_CR_SCBU(:) = 0.
      ALLOCATE(Z_TH_EVAP(IPACK))          ; Z_TH_EVAP(:) = 0.
      ALLOCATE(Z_RR_EVAP(IPACK))          ; Z_RR_EVAP(:) = 0.
      ALLOCATE(Z_RI_CNVI(IPACK))          ; Z_RI_CNVI(:) = 0.
      ALLOCATE(Z_CI_CNVI(IPACK))          ; Z_CI_CNVI(:) = 0.
      ALLOCATE(Z_TH_DEPS(IPACK))          ; Z_TH_DEPS(:) = 0.
      ALLOCATE(Z_RS_DEPS(IPACK))          ; Z_RS_DEPS(:) = 0.
      ALLOCATE(Z_TH_DEPI(IPACK))          ; Z_TH_DEPI(:) = 0.
      ALLOCATE(Z_RI_DEPI(IPACK))          ; Z_RI_DEPI(:) = 0.
      ALLOCATE(Z_RI_CNVS(IPACK))          ; Z_RI_CNVS(:) = 0.
      ALLOCATE(Z_CI_CNVS(IPACK))          ; Z_CI_CNVS(:) = 0.
      ALLOCATE(Z_RI_AGGS(IPACK))          ; Z_RI_AGGS(:) = 0.
      ALLOCATE(Z_CI_AGGS(IPACK))          ; Z_CI_AGGS(:) = 0.
      ALLOCATE(Z_TH_DEPG(IPACK))          ; Z_TH_DEPG(:) = 0.
      ALLOCATE(Z_RG_DEPG(IPACK))          ; Z_RG_DEPG(:) = 0.
      ALLOCATE(Z_TH_BERFI(IPACK))         ; Z_TH_BERFI(:) = 0.
      ALLOCATE(Z_RC_BERFI(IPACK))         ; Z_RC_BERFI(:) = 0.
      ALLOCATE(Z_TH_RIM(IPACK))           ; Z_TH_RIM = 0.
      ALLOCATE(Z_RC_RIM(IPACK))           ; Z_RC_RIM = 0.
      ALLOCATE(Z_CC_RIM(IPACK))           ; Z_CC_RIM = 0.
      ALLOCATE(Z_RS_RIM(IPACK))           ; Z_RS_RIM = 0.
      ALLOCATE(Z_RG_RIM(IPACK))           ; Z_RG_RIM = 0.
      ALLOCATE(Z_RI_HMS(IPACK))           ; Z_RI_HMS = 0.
      ALLOCATE(Z_CI_HMS(IPACK))           ; Z_CI_HMS = 0.
      ALLOCATE(Z_RS_HMS(IPACK))           ; Z_RS_HMS = 0.
      ALLOCATE(Z_TH_ACC(IPACK))           ; Z_TH_ACC = 0.
      ALLOCATE(Z_RR_ACC(IPACK))           ; Z_RR_ACC = 0.
      ALLOCATE(Z_CR_ACC(IPACK))           ; Z_CR_ACC = 0.
      ALLOCATE(Z_RS_ACC(IPACK))           ; Z_RS_ACC = 0.
      ALLOCATE(Z_RG_ACC(IPACK))           ; Z_RG_ACC = 0.
      ALLOCATE(Z_RS_CMEL(IPACK))          ; Z_RS_CMEL = 0.
      ALLOCATE(Z_TH_CFRZ(IPACK))          ; Z_TH_CFRZ = 0.
      ALLOCATE(Z_RR_CFRZ(IPACK))          ; Z_RR_CFRZ = 0.
      ALLOCATE(Z_CR_CFRZ(IPACK))          ; Z_CR_CFRZ = 0.
      ALLOCATE(Z_RI_CFRZ(IPACK))          ; Z_RI_CFRZ = 0.
      ALLOCATE(Z_CI_CFRZ(IPACK))          ; Z_CI_CFRZ = 0.
      ALLOCATE(Z_TH_WETG(IPACK))          ; Z_TH_WETG = 0.
      ALLOCATE(Z_RC_WETG(IPACK))          ; Z_RC_WETG = 0.
      ALLOCATE(Z_CC_WETG(IPACK))          ; Z_CC_WETG = 0.
      ALLOCATE(Z_RR_WETG(IPACK))          ; Z_RR_WETG = 0.
      ALLOCATE(Z_CR_WETG(IPACK))          ; Z_CR_WETG = 0.
      ALLOCATE(Z_RI_WETG(IPACK))          ; Z_RI_WETG = 0.
      ALLOCATE(Z_CI_WETG(IPACK))          ; Z_CI_WETG = 0.
      ALLOCATE(Z_RS_WETG(IPACK))          ; Z_RS_WETG = 0.
      ALLOCATE(Z_RG_WETG(IPACK))          ; Z_RG_WETG = 0.
      ALLOCATE(Z_RH_WETG(IPACK))          ; Z_RH_WETG = 0.
      ALLOCATE(Z_TH_DRYG(IPACK))          ; Z_TH_DRYG = 0.
      ALLOCATE(Z_RC_DRYG(IPACK))          ; Z_RC_DRYG = 0.
      ALLOCATE(Z_CC_DRYG(IPACK))          ; Z_CC_DRYG = 0.
      ALLOCATE(Z_RR_DRYG(IPACK))          ; Z_RR_DRYG = 0.
      ALLOCATE(Z_CR_DRYG(IPACK))          ; Z_CR_DRYG = 0.
      ALLOCATE(Z_RI_DRYG(IPACK))          ; Z_RI_DRYG = 0.
      ALLOCATE(Z_CI_DRYG(IPACK))          ; Z_CI_DRYG = 0.
      ALLOCATE(Z_RS_DRYG(IPACK))          ; Z_RS_DRYG = 0.
      ALLOCATE(Z_RG_DRYG(IPACK))          ; Z_RG_DRYG = 0.
      ALLOCATE(Z_RI_HMG(IPACK))           ; Z_RI_HMG = 0.
      ALLOCATE(Z_CI_HMG(IPACK))           ; Z_CI_HMG = 0.
      ALLOCATE(Z_RG_HMG(IPACK))           ; Z_RG_HMG = 0.
      ALLOCATE(Z_TH_GMLT(IPACK))          ; Z_TH_GMLT = 0.
      ALLOCATE(Z_RR_GMLT(IPACK))          ; Z_RR_GMLT = 0.
      ALLOCATE(Z_CR_GMLT(IPACK))          ; Z_CR_GMLT = 0.

      ALLOCATE(Z_RV_CORR2(IPACK))         ; Z_RV_CORR2 = 0.
      ALLOCATE(Z_RC_CORR2(IPACK))         ; Z_RC_CORR2 = 0.
      ALLOCATE(Z_RR_CORR2(IPACK))         ; Z_RR_CORR2 = 0.
      ALLOCATE(Z_RI_CORR2(IPACK))         ; Z_RI_CORR2 = 0.
      ALLOCATE(Z_CC_CORR2(IPACK))         ; Z_CC_CORR2 = 0.
      ALLOCATE(Z_CR_CORR2(IPACK))         ; Z_CR_CORR2 = 0.
      ALLOCATE(Z_CI_CORR2(IPACK))         ; Z_CI_CORR2 = 0.
      !
      !***       4.1 Tendecies computation
      !
      
      CALL LIMA_INST_PROCS (PTSTEP, LLCOMPUTE1D,                                &
                            ZEXNREF1D, ZP1D,                                    &
                            ZTHT1D, ZRVT1D, ZRCT1D, ZRRT1D, ZRIT1D, ZRST1D, ZRGT1D, &
                            ZCCT1D, ZCRT1D, ZCIT1D,                             &
                            ZIFNN1D,                                            &
                            Z_CR_BRKU,                                          & ! spontaneous break up of drops (BRKU) : Nr
                            Z_TH_HONR, Z_RR_HONR, Z_CR_HONR,                    & ! rain drops homogeneous freezing (HONR) : rr, Nr, rg=-rr, th
                            Z_TH_IMLT, Z_RC_IMLT, Z_CC_IMLT,                    & ! ice melting (IMLT) : rc, Nc, ri=-rc, Ni=-Nc, th, IFNF, IFNA
                            ZB_TH, ZB_RV, ZB_RC, ZB_RR, ZB_RI, ZB_RG,           &
                            ZB_CC, ZB_CR, ZB_CI,                                &
                            ZB_IFNN,                                            &
                            ZCF1D, ZIF1D, ZPF1D                                 )
      
      CALL LIMA_TENDENCIES (PTSTEP, LLCOMPUTE1D,                                   &
                            ZEXNREF1D, ZRHODREF1D, ZP1D, ZTHT1D,                   &
                            ZRVT1D, ZRCT1D, ZRRT1D, ZRIT1D, ZRST1D, ZRGT1D, ZRHT1D,&
                            ZCCT1D, ZCRT1D, ZCIT1D,                                &
                            Z_TH_HONC, Z_RC_HONC, Z_CC_HONC,                       & 
                            Z_CC_SELF,                                             & 
                            Z_RC_AUTO, Z_CC_AUTO, Z_CR_AUTO,                       & 
                            Z_RC_ACCR, Z_CC_ACCR,                                  & 
                            Z_CR_SCBU,                                             & 
                            Z_TH_EVAP, Z_RR_EVAP,                                  & 
                            Z_RI_CNVI, Z_CI_CNVI,                                  & 
                            Z_TH_DEPS, Z_RS_DEPS,                                  & 
                            Z_TH_DEPI, Z_RI_DEPI,                                  & 
                            Z_RI_CNVS, Z_CI_CNVS,                                  & 
                            Z_RI_AGGS, Z_CI_AGGS,                                  & 
                            Z_TH_DEPG, Z_RG_DEPG,                                  & 
                            Z_TH_BERFI, Z_RC_BERFI,                                & 
                            Z_TH_RIM, Z_RC_RIM, Z_CC_RIM, Z_RS_RIM, Z_RG_RIM,      & 
                            Z_RI_HMS, Z_CI_HMS, Z_RS_HMS,                          & 
                            Z_TH_ACC, Z_RR_ACC, Z_CR_ACC, Z_RS_ACC, Z_RG_ACC,      & 
                            Z_RS_CMEL,                                             & 
                            Z_TH_CFRZ, Z_RR_CFRZ, Z_CR_CFRZ, Z_RI_CFRZ, Z_CI_CFRZ, & 
                            Z_TH_WETG, Z_RC_WETG, Z_CC_WETG, Z_RR_WETG, Z_CR_WETG, & 
                            Z_RI_WETG, Z_CI_WETG, Z_RS_WETG, Z_RG_WETG, Z_RH_WETG, & 
                            Z_TH_DRYG, Z_RC_DRYG, Z_CC_DRYG, Z_RR_DRYG, Z_CR_DRYG, & 
                            Z_RI_DRYG, Z_CI_DRYG, Z_RS_DRYG, Z_RG_DRYG,            & 
                            Z_RI_HMG, Z_CI_HMG, Z_RG_HMG,                          & 
                            Z_TH_GMLT, Z_RR_GMLT, Z_CR_GMLT,                       & 
!!!     Z_RC_WETH, Z_CC_WETH, Z_RR_WETH, Z_CR_WETH,  &           ! wet growth of hail (WETH) : rc, Nc, rr, Nr, ri, Ni, rs, rg, rh, th
!!!     Z_RI_WETH, Z_CI_WETH, Z_RS_WETH, Z_RG_WETH, Z_RH_WETH, & ! wet growth of hail (WETH) : rc, Nc, rr, Nr, ri, Ni, rs, rg, rh, th
!!!     Z_RG_COHG, &                                             ! conversion of hail into graupel (COHG) : rg, rh
!!!     Z_RR_HMLT, Z_CR_HMLT                                     ! hail melting (HMLT) : rr, Nr, rh=-rr, th
                            ZA_TH, ZA_RV, ZA_RC, ZA_CC, ZA_RR, ZA_CR,              &
                            ZA_RI, ZA_CI, ZA_RS, ZA_RG, ZA_RH,                     &
                            ZEVAP1D,                                               &
                            ZCF1D, ZIF1D, ZPF1D                                    )

      !
      !***       4.2 Integration time
      !
      ! If we can, we will use these tendecies until the end of the timestep
      ZMAXTIME(:)=PTSTEP-ZTIME1D(:) ! Remaining time until the end of the timestep

      ! We need to adjust tendencies when temperature reaches 0
      IF(LFEEDBACKT) THEN
         !Is ZB_TH enough to change temperature sign?
         WHERE( ((ZTHT1D(:) - XTT/ZEXN1D(:)) * (ZTHT1D(:) + ZB_TH(:) - XTT/ZEXN1D(:))) < 0. )
            ZMAXTIME(:)=0.
         ENDWHERE
         !Can ZA_TH make temperature change of sign?
         ZTIME_THRESHOLD(:)=-1.
         WHERE(ABS(ZA_TH(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(XTT/ZEXN1D(:) - ZB_TH(:) - ZTHT1D(:))/ZA_TH(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>0.)
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
         ENDWHERE
      ENDIF

      ! We need to adjust tendencies when a species disappears
      ! When a species is missing, only the external tendencies can be negative (and we must keep track of it)
      WHERE(ZA_RV(:)<-1.E-20 .AND. ZRVT1D(:)>XRTMIN(1))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RV(:)+ZRVT1D(:))/ZA_RV(:))
      END WHERE
      WHERE(ZA_RC(:)<-1.E-20 .AND. ZRCT1D(:)>XRTMIN(2))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RC(:)+ZRCT1D(:))/ZA_RC(:))
      END WHERE
      WHERE(ZA_RR(:)<-1.E-20 .AND. ZRRT1D(:)>XRTMIN(3))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RR(:)+ZRRT1D(:))/ZA_RR(:))
      END WHERE
      WHERE(ZA_RI(:)<-1.E-20 .AND. ZRIT1D(:)>XRTMIN(4))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RI(:)+ZRIT1D(:))/ZA_RI(:))
      END WHERE
      WHERE(ZA_RS(:)<-1.E-20 .AND. ZRST1D(:)>XRTMIN(5))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RS(:)+ZRST1D(:))/ZA_RS(:))
      END WHERE
      WHERE(ZA_RG(:)<-1.E-20 .AND. ZRGT1D(:)>XRTMIN(6))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RG(:)+ZRGT1D(:))/ZA_RG(:))
      END WHERE
      WHERE(ZA_RH(:)<-1.E-20 .AND. ZRHT1D(:)>XRTMIN(7))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RH(:)+ZRHT1D(:))/ZA_RH(:))
      END WHERE

      ! We stop when the end of the timestep is reached
      WHERE(PTSTEP-ZTIME1D(:)-ZMAXTIME(:)<=0.)
         LLCOMPUTE1D(:)=.FALSE.
      ENDWHERE

      ! We must recompute tendencies when the end of the sub-timestep is reached
      IF(XTSTEP_TS/=0.) THEN
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ZTIME1D(:)+ZMAXTIME(:)>ZTIME_LASTCALL1D(:)+ZTSTEP)
            ZMAXTIME(:)=ZTIME_LASTCALL1D(:)-ZTIME1D(:)+ZTSTEP
            LLCOMPUTE1D(:)=.FALSE.
         ENDWHERE
      ENDIF

      ! We must recompute tendencies when the maximum allowed change is reached
      ! When a species is missing, only the external tendencies can be active and we do not want to recompute
      ! the microphysical tendencies when external tendencies are negative (results won't change because species was already missing)
      IF(XMRSTEP/=0.) THEN
         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RV(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RV(:))*XMRSTEP+Z0RVT1D(:)-ZRVT1D(:)-ZB_RV(:))/ZA_RV(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRVT1D(:)>XRTMIN(1) .OR. ZA_RV(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            LLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RC(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RC(:))*XMRSTEP+Z0RCT1D(:)-ZRCT1D(:)-ZB_RC(:))/ZA_RC(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRCT1D(:)>XRTMIN(2) .OR. ZA_RC(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            LLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RR(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RR(:))*XMRSTEP+Z0RRT1D(:)-ZRRT1D(:)-ZB_RR(:))/ZA_RR(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRRT1D(:)>XRTMIN(3) .OR. ZA_RR(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            LLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RI(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RI(:))*XMRSTEP+Z0RIT1D(:)-ZRIT1D(:)-ZB_RI(:))/ZA_RI(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRIT1D(:)>XRTMIN(4) .OR. ZA_RI(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            LLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RS(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RS(:))*XMRSTEP+Z0RST1D(:)-ZRST1D(:)-ZB_RS(:))/ZA_RS(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRST1D(:)>XRTMIN(5) .OR. ZA_RS(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            LLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RG(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RG(:))*XMRSTEP+Z0RGT1D(:)-ZRGT1D(:)-ZB_RG(:))/ZA_RG(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRGT1D(:)>XRTMIN(6) .OR. ZA_RG(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            LLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RH(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RH(:))*XMRSTEP+Z0RHT1D(:)-ZRHT1D(:)-ZB_RH(:))/ZA_RH(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRHT1D(:)>XRTMIN(7) .OR. ZA_RH(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            LLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         WHERE(IITER1D(:)<INB_ITER_MAX .AND. MAX(ABS(ZB_RV(:)),  &
              ABS(ZB_RC(:)), ABS(ZB_RR(:)), ABS(ZB_RI(:)), &
              ABS(ZB_RS(:)), ABS(ZB_RG(:)), ABS(ZB_RH(:)))>XMRSTEP)
            ZMAXTIME(:)=0.
            LLCOMPUTE1D(:)=.FALSE.
         ENDWHERE
      ENDIF
      !
      !***       4.3 New values of variables for next iteration
      !
      ZTHT1D = ZTHT1D + ZA_TH(:) * ZMAXTIME(:) + ZB_TH(:)
      ZRVT1D = ZRVT1D + ZA_RV(:) * ZMAXTIME(:) + ZB_RV(:)
      ZRCT1D = ZRCT1D + ZA_RC(:) * ZMAXTIME(:) + ZB_RC(:)
      ZCCT1D = ZCCT1D + ZA_CC(:) * ZMAXTIME(:) + ZB_CC(:)
      ZRRT1D = ZRRT1D + ZA_RR(:) * ZMAXTIME(:) + ZB_RR(:)
      ZCRT1D = ZCRT1D + ZA_CR(:) * ZMAXTIME(:) + ZB_CR(:)
      ZRIT1D = ZRIT1D + ZA_RI(:) * ZMAXTIME(:) + ZB_RI(:)
      ZCIT1D = ZCIT1D + ZA_CI(:) * ZMAXTIME(:) + ZB_CI(:)
      ZRST1D = ZRST1D + ZA_RS(:) * ZMAXTIME(:) + ZB_RS(:)
      ZRGT1D = ZRGT1D + ZA_RG(:) * ZMAXTIME(:) + ZB_RG(:)
      ZRHT1D = ZRHT1D + ZA_RH(:) * ZMAXTIME(:) + ZB_RH(:)
      !
      DO II=1,NMOD_IFN
         ZIFNN1D(:,II) = ZIFNN1D(:,II) + ZB_IFNN(:,II)
      END DO
      !
      !***       4.5 
      !
      WHERE (ZRCT1D .LE. XRTMIN(2))
         Z_RV_CORR2(:) = ZRCT1D(:)
         Z_RC_CORR2(:) = -ZRCT1D(:)
         Z_CC_CORR2(:) = -ZCCT1D(:)

         ZRVT1D = ZRVT1D + ZRCT1D
         ZRCT1D = 0.
         ZCCT1D = 0.
      END WHERE
      WHERE (ZRRT1D .LE. XRTMIN(3))
         Z_RV_CORR2(:) = Z_RV_CORR2(:) + ZRRT1D(:)
         Z_RR_CORR2(:) = -ZRRT1D(:)
         Z_CR_CORR2(:) = -ZCRT1D(:)

         ZRVT1D = ZRVT1D + ZRRT1D
         ZRRT1D = 0.
         ZCRT1D = 0.
      END WHERE
      WHERE (ZRIT1D .LE. XRTMIN(4))
         Z_RV_CORR2(:) = Z_RV_CORR2(:) + ZRIT1D(:)
         Z_RI_CORR2(:) = -ZRIT1D(:)
         Z_CI_CORR2(:) = -ZCIT1D(:)

         ZRVT1D = ZRVT1D + ZRIT1D
         ZRIT1D = 0.
         ZCIT1D = 0.
      END WHERE
      
      !
      !***       4.5 Next loop
      !
      ZTIME1D(:)=ZTIME1D(:)+ZMAXTIME(:)
      !
      !***       4.4 Unpacking
      !
      DO II=1,IPACK
         ZTHT(I1(II),I2(II),I3(II))      = ZTHT1D(II)
         ZRVT(I1(II),I2(II),I3(II))      = ZRVT1D(II)
         ZRCT(I1(II),I2(II),I3(II))      = ZRCT1D(II)
         ZRRT(I1(II),I2(II),I3(II))      = ZRRT1D(II)
         ZRIT(I1(II),I2(II),I3(II))      = ZRIT1D(II)
         ZRST(I1(II),I2(II),I3(II))      = ZRST1D(II)
         ZRGT(I1(II),I2(II),I3(II))      = ZRGT1D(II)
         ZRHT(I1(II),I2(II),I3(II))      = ZRHT1D(II)
         ZCCT(I1(II),I2(II),I3(II))      = ZCCT1D(II)
         ZCRT(I1(II),I2(II),I3(II))      = ZCRT1D(II)
         ZCIT(I1(II),I2(II),I3(II))      = ZCIT1D(II)
         ZIFNNT(I1(II),I2(II),I3(II),:)  = ZIFNN1D(II,:)
         PEVAP3D(I1(II),I2(II),I3(II))   = ZEVAP1D(II)
         ZTIME(I1(II),I2(II),I3(II))     = ZTIME1D(II)
         LLCOMPUTE(I1(II),I2(II),I3(II)) = LLCOMPUTE1D(II)
         IITER(I1(II),I2(II),I3(II))     = IITER1D(II)
      END DO
      !
      CALL LIMA_DROPS_TO_DROPLETS_CONV(PRHODREF, ZRCT, ZRRT, ZCCT, ZCRT, &
                                       Z_RR_CVRC, Z_CR_CVRC    )
      ZRCT(:,:,:) = ZRCT(:,:,:) - Z_RR_CVRC(:,:,:)
      ZRRT(:,:,:) = ZRRT(:,:,:) + Z_RR_CVRC(:,:,:)
      ZCCT(:,:,:) = ZCCT(:,:,:) - Z_CR_CVRC(:,:,:)
      ZCRT(:,:,:) = ZCRT(:,:,:) + Z_CR_CVRC(:,:,:)
      !
      !***       4.4 Unpacking for budgets
      !
      IF(LBU_ENABLE) THEN
        ZTOT_RR_CVRC(:,:,:) = ZTOT_RR_CVRC(:,:,:) + Z_RR_CVRC(:,:,:)
        ZTOT_CR_CVRC(:,:,:) = ZTOT_CR_CVRC(:,:,:) + Z_CR_CVRC(:,:,:)

         DO II=1,IPACK
            ! Instantaneous processes
            ZTOT_CR_BRKU(I1(II),I2(II),I3(II)) =   ZTOT_CR_BRKU(I1(II),I2(II),I3(II))   + Z_CR_BRKU(II)
            ZTOT_TH_HONR(I1(II),I2(II),I3(II)) =   ZTOT_TH_HONR(I1(II),I2(II),I3(II))   + Z_TH_HONR(II)
            ZTOT_RR_HONR(I1(II),I2(II),I3(II)) =   ZTOT_RR_HONR(I1(II),I2(II),I3(II))   + Z_RR_HONR(II)
            ZTOT_CR_HONR(I1(II),I2(II),I3(II)) =   ZTOT_CR_HONR(I1(II),I2(II),I3(II))   + Z_CR_HONR(II)
            ZTOT_TH_IMLT(I1(II),I2(II),I3(II)) =   ZTOT_TH_IMLT(I1(II),I2(II),I3(II))   + Z_TH_IMLT(II)
            ZTOT_RC_IMLT(I1(II),I2(II),I3(II)) =   ZTOT_RC_IMLT(I1(II),I2(II),I3(II))   + Z_RC_IMLT(II)
            ZTOT_CC_IMLT(I1(II),I2(II),I3(II)) =   ZTOT_CC_IMLT(I1(II),I2(II),I3(II))   + Z_CC_IMLT(II)
            DO JI = 1, NMOD_IFN
              ZTOT_IFNN_IMLT(I1(II),I2(II),I3(II),JI) = ZTOT_IFNN_IMLT(I1(II),I2(II),I3(II),JI) + ZB_IFNN(II,JI)
            END DO

            ! Tendencies
            ZTOT_TH_HONC(I1(II),I2(II),I3(II)) =   ZTOT_TH_HONC(I1(II),I2(II),I3(II))   + Z_TH_HONC(II)  * ZMAXTIME(II)
            ZTOT_RC_HONC(I1(II),I2(II),I3(II)) =   ZTOT_RC_HONC(I1(II),I2(II),I3(II))   + Z_RC_HONC(II)  * ZMAXTIME(II)
            ZTOT_CC_HONC(I1(II),I2(II),I3(II)) =   ZTOT_CC_HONC(I1(II),I2(II),I3(II))   + Z_CC_HONC(II)  * ZMAXTIME(II)
            ZTOT_CC_SELF(I1(II),I2(II),I3(II)) =   ZTOT_CC_SELF(I1(II),I2(II),I3(II))   + Z_CC_SELF(II)  * ZMAXTIME(II)
            ZTOT_RC_AUTO(I1(II),I2(II),I3(II)) =   ZTOT_RC_AUTO(I1(II),I2(II),I3(II))   + Z_RC_AUTO(II)  * ZMAXTIME(II)
            ZTOT_CC_AUTO(I1(II),I2(II),I3(II)) =   ZTOT_CC_AUTO(I1(II),I2(II),I3(II))   + Z_CC_AUTO(II)  * ZMAXTIME(II)
            ZTOT_CR_AUTO(I1(II),I2(II),I3(II)) =   ZTOT_CR_AUTO(I1(II),I2(II),I3(II))   + Z_CR_AUTO(II)  * ZMAXTIME(II)
            ZTOT_RC_ACCR(I1(II),I2(II),I3(II)) =   ZTOT_RC_ACCR(I1(II),I2(II),I3(II))   + Z_RC_ACCR(II)  * ZMAXTIME(II)
            ZTOT_CC_ACCR(I1(II),I2(II),I3(II)) =   ZTOT_CC_ACCR(I1(II),I2(II),I3(II))   + Z_CC_ACCR(II)  * ZMAXTIME(II)
            ZTOT_CR_SCBU(I1(II),I2(II),I3(II)) =   ZTOT_CR_SCBU(I1(II),I2(II),I3(II))   + Z_CR_SCBU(II)  * ZMAXTIME(II)
            ZTOT_TH_EVAP(I1(II),I2(II),I3(II)) =   ZTOT_TH_EVAP(I1(II),I2(II),I3(II))   + Z_TH_EVAP(II)  * ZMAXTIME(II)
!!$            ZTOT_RC_EVAP(I1(II),I2(II),I3(II)) =   ZTOT_RC_EVAP(I1(II),I2(II),I3(II))   + Z_RC_EVAP(II)  * ZMAXTIME(II)
!!$            ZTOT_CC_EVAP(I1(II),I2(II),I3(II)) =   ZTOT_CC_EVAP(I1(II),I2(II),I3(II))   + Z_CC_EVAP(II)  * ZMAXTIME(II)
            ZTOT_RR_EVAP(I1(II),I2(II),I3(II)) =   ZTOT_RR_EVAP(I1(II),I2(II),I3(II))   + Z_RR_EVAP(II)  * ZMAXTIME(II)
!!$            ZTOT_CR_EVAP(I1(II),I2(II),I3(II)) =   ZTOT_CR_EVAP(I1(II),I2(II),I3(II))   + Z_CR_EVAP(II)  * ZMAXTIME(II)
            ZTOT_RI_CNVI(I1(II),I2(II),I3(II)) =   ZTOT_RI_CNVI(I1(II),I2(II),I3(II))   + Z_RI_CNVI(II)  * ZMAXTIME(II)
            ZTOT_CI_CNVI(I1(II),I2(II),I3(II)) =   ZTOT_CI_CNVI(I1(II),I2(II),I3(II))   + Z_CI_CNVI(II)  * ZMAXTIME(II)
            ZTOT_TH_DEPS(I1(II),I2(II),I3(II)) =   ZTOT_TH_DEPS(I1(II),I2(II),I3(II))   + Z_TH_DEPS(II)  * ZMAXTIME(II)
            ZTOT_RS_DEPS(I1(II),I2(II),I3(II)) =   ZTOT_RS_DEPS(I1(II),I2(II),I3(II))   + Z_RS_DEPS(II)  * ZMAXTIME(II)
            ZTOT_TH_DEPI(I1(II),I2(II),I3(II)) =   ZTOT_TH_DEPI(I1(II),I2(II),I3(II))   + Z_TH_DEPI(II)  * ZMAXTIME(II)
            ZTOT_RI_DEPI(I1(II),I2(II),I3(II)) =   ZTOT_RI_DEPI(I1(II),I2(II),I3(II))   + Z_RI_DEPI(II)  * ZMAXTIME(II)
            ZTOT_RI_CNVS(I1(II),I2(II),I3(II)) =   ZTOT_RI_CNVS(I1(II),I2(II),I3(II))   + Z_RI_CNVS(II)  * ZMAXTIME(II)
            ZTOT_CI_CNVS(I1(II),I2(II),I3(II)) =   ZTOT_CI_CNVS(I1(II),I2(II),I3(II))   + Z_CI_CNVS(II)  * ZMAXTIME(II)
            ZTOT_RI_AGGS(I1(II),I2(II),I3(II)) =   ZTOT_RI_AGGS(I1(II),I2(II),I3(II))   + Z_RI_AGGS(II)  * ZMAXTIME(II)
            ZTOT_CI_AGGS(I1(II),I2(II),I3(II)) =   ZTOT_CI_AGGS(I1(II),I2(II),I3(II))   + Z_CI_AGGS(II)  * ZMAXTIME(II)
            ZTOT_TH_DEPG(I1(II),I2(II),I3(II)) =   ZTOT_TH_DEPG(I1(II),I2(II),I3(II))   + Z_TH_DEPG(II)  * ZMAXTIME(II)
            ZTOT_RG_DEPG(I1(II),I2(II),I3(II)) =   ZTOT_RG_DEPG(I1(II),I2(II),I3(II))   + Z_RG_DEPG(II)  * ZMAXTIME(II)
            ZTOT_TH_BERFI(I1(II),I2(II),I3(II))=   ZTOT_TH_BERFI(I1(II),I2(II),I3(II))  + Z_TH_BERFI(II) * ZMAXTIME(II)
            ZTOT_RC_BERFI(I1(II),I2(II),I3(II))=   ZTOT_RC_BERFI(I1(II),I2(II),I3(II))  + Z_RC_BERFI(II) * ZMAXTIME(II)
            ZTOT_TH_RIM(I1(II),I2(II),I3(II))  =   ZTOT_TH_RIM(I1(II),I2(II),I3(II))    + Z_TH_RIM(II)   * ZMAXTIME(II)
            ZTOT_RC_RIM(I1(II),I2(II),I3(II))  =   ZTOT_RC_RIM(I1(II),I2(II),I3(II))    + Z_RC_RIM(II)   * ZMAXTIME(II)
            ZTOT_CC_RIM(I1(II),I2(II),I3(II))  =   ZTOT_CC_RIM(I1(II),I2(II),I3(II))    + Z_CC_RIM(II)   * ZMAXTIME(II)
            ZTOT_RS_RIM(I1(II),I2(II),I3(II))  =   ZTOT_RS_RIM(I1(II),I2(II),I3(II))    + Z_RS_RIM(II)   * ZMAXTIME(II)
            ZTOT_RG_RIM(I1(II),I2(II),I3(II))  =   ZTOT_RG_RIM(I1(II),I2(II),I3(II))    + Z_RG_RIM(II)   * ZMAXTIME(II)
            ZTOT_RI_HMS(I1(II),I2(II),I3(II))  =   ZTOT_RI_HMS(I1(II),I2(II),I3(II))    + Z_RI_HMS(II)   * ZMAXTIME(II)
            ZTOT_CI_HMS(I1(II),I2(II),I3(II))  =   ZTOT_CI_HMS(I1(II),I2(II),I3(II))    + Z_CI_HMS(II)   * ZMAXTIME(II)
            ZTOT_RS_HMS(I1(II),I2(II),I3(II))  =   ZTOT_RS_HMS(I1(II),I2(II),I3(II))    + Z_RS_HMS(II)   * ZMAXTIME(II)
            ZTOT_TH_ACC(I1(II),I2(II),I3(II))  =   ZTOT_TH_ACC(I1(II),I2(II),I3(II))    + Z_TH_ACC(II)   * ZMAXTIME(II)
            ZTOT_RR_ACC(I1(II),I2(II),I3(II))  =   ZTOT_RR_ACC(I1(II),I2(II),I3(II))    + Z_RR_ACC(II)   * ZMAXTIME(II)
            ZTOT_CR_ACC(I1(II),I2(II),I3(II))  =   ZTOT_CR_ACC(I1(II),I2(II),I3(II))    + Z_CR_ACC(II)   * ZMAXTIME(II)
            ZTOT_RS_ACC(I1(II),I2(II),I3(II))  =   ZTOT_RS_ACC(I1(II),I2(II),I3(II))    + Z_RS_ACC(II)   * ZMAXTIME(II)
            ZTOT_RG_ACC(I1(II),I2(II),I3(II))  =   ZTOT_RG_ACC(I1(II),I2(II),I3(II))    + Z_RG_ACC(II)   * ZMAXTIME(II)
            ZTOT_RS_CMEL(I1(II),I2(II),I3(II)) =   ZTOT_RS_CMEL(I1(II),I2(II),I3(II))   + Z_RS_CMEL(II)  * ZMAXTIME(II)
            ZTOT_TH_CFRZ(I1(II),I2(II),I3(II)) =   ZTOT_TH_CFRZ(I1(II),I2(II),I3(II))   + Z_TH_CFRZ(II)  * ZMAXTIME(II)
            ZTOT_RR_CFRZ(I1(II),I2(II),I3(II)) =   ZTOT_RR_CFRZ(I1(II),I2(II),I3(II))   + Z_RR_CFRZ(II)  * ZMAXTIME(II)
            ZTOT_CR_CFRZ(I1(II),I2(II),I3(II)) =   ZTOT_CR_CFRZ(I1(II),I2(II),I3(II))   + Z_CR_CFRZ(II)  * ZMAXTIME(II)
            ZTOT_RI_CFRZ(I1(II),I2(II),I3(II)) =   ZTOT_RI_CFRZ(I1(II),I2(II),I3(II))   + Z_RI_CFRZ(II)  * ZMAXTIME(II)
            ZTOT_CI_CFRZ(I1(II),I2(II),I3(II)) =   ZTOT_CI_CFRZ(I1(II),I2(II),I3(II))   + Z_CI_CFRZ(II)  * ZMAXTIME(II)
            ZTOT_TH_WETG(I1(II),I2(II),I3(II)) =   ZTOT_TH_WETG(I1(II),I2(II),I3(II))   + Z_TH_WETG(II)  * ZMAXTIME(II)
            ZTOT_RC_WETG(I1(II),I2(II),I3(II)) =   ZTOT_RC_WETG(I1(II),I2(II),I3(II))   + Z_RC_WETG(II)  * ZMAXTIME(II)
            ZTOT_CC_WETG(I1(II),I2(II),I3(II)) =   ZTOT_CC_WETG(I1(II),I2(II),I3(II))   + Z_CC_WETG(II)  * ZMAXTIME(II)
            ZTOT_RR_WETG(I1(II),I2(II),I3(II)) =   ZTOT_RR_WETG(I1(II),I2(II),I3(II))   + Z_RR_WETG(II)  * ZMAXTIME(II)
            ZTOT_CR_WETG(I1(II),I2(II),I3(II)) =   ZTOT_CR_WETG(I1(II),I2(II),I3(II))   + Z_CR_WETG(II)  * ZMAXTIME(II)
            ZTOT_RI_WETG(I1(II),I2(II),I3(II)) =   ZTOT_RI_WETG(I1(II),I2(II),I3(II))   + Z_RI_WETG(II)  * ZMAXTIME(II)
            ZTOT_CI_WETG(I1(II),I2(II),I3(II)) =   ZTOT_CI_WETG(I1(II),I2(II),I3(II))   + Z_CI_WETG(II)  * ZMAXTIME(II)
            ZTOT_RS_WETG(I1(II),I2(II),I3(II)) =   ZTOT_RS_WETG(I1(II),I2(II),I3(II))   + Z_RS_WETG(II)  * ZMAXTIME(II)
            ZTOT_RG_WETG(I1(II),I2(II),I3(II)) =   ZTOT_RG_WETG(I1(II),I2(II),I3(II))   + Z_RG_WETG(II)  * ZMAXTIME(II)
            ZTOT_RH_WETG(I1(II),I2(II),I3(II)) =   ZTOT_RH_WETG(I1(II),I2(II),I3(II))   + Z_RH_WETG(II)  * ZMAXTIME(II)
            ZTOT_TH_DRYG(I1(II),I2(II),I3(II)) =   ZTOT_TH_DRYG(I1(II),I2(II),I3(II))   + Z_TH_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RC_DRYG(I1(II),I2(II),I3(II)) =   ZTOT_RC_DRYG(I1(II),I2(II),I3(II))   + Z_RC_DRYG(II)  * ZMAXTIME(II)
            ZTOT_CC_DRYG(I1(II),I2(II),I3(II)) =   ZTOT_CC_DRYG(I1(II),I2(II),I3(II))   + Z_CC_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RR_DRYG(I1(II),I2(II),I3(II)) =   ZTOT_RR_DRYG(I1(II),I2(II),I3(II))   + Z_RR_DRYG(II)  * ZMAXTIME(II)
            ZTOT_CR_DRYG(I1(II),I2(II),I3(II)) =   ZTOT_CR_DRYG(I1(II),I2(II),I3(II))   + Z_CR_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RI_DRYG(I1(II),I2(II),I3(II)) =   ZTOT_RI_DRYG(I1(II),I2(II),I3(II))   + Z_RI_DRYG(II)  * ZMAXTIME(II)
            ZTOT_CI_DRYG(I1(II),I2(II),I3(II)) =   ZTOT_CI_DRYG(I1(II),I2(II),I3(II))   + Z_CI_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RS_DRYG(I1(II),I2(II),I3(II)) =   ZTOT_RS_DRYG(I1(II),I2(II),I3(II))   + Z_RS_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RG_DRYG(I1(II),I2(II),I3(II)) =   ZTOT_RG_DRYG(I1(II),I2(II),I3(II))   + Z_RG_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RI_HMG(I1(II),I2(II),I3(II))  =   ZTOT_RI_HMG(I1(II),I2(II),I3(II))    + Z_RI_HMG(II)   * ZMAXTIME(II)
            ZTOT_CI_HMG(I1(II),I2(II),I3(II))  =   ZTOT_CI_HMG(I1(II),I2(II),I3(II))    + Z_CI_HMG(II)   * ZMAXTIME(II)
            ZTOT_RG_HMG(I1(II),I2(II),I3(II))  =   ZTOT_RG_HMG(I1(II),I2(II),I3(II))    + Z_RG_HMG(II)   * ZMAXTIME(II)
            ZTOT_TH_GMLT(I1(II),I2(II),I3(II)) =   ZTOT_TH_GMLT(I1(II),I2(II),I3(II))   + Z_TH_GMLT(II)  * ZMAXTIME(II)
            ZTOT_RR_GMLT(I1(II),I2(II),I3(II)) =   ZTOT_RR_GMLT(I1(II),I2(II),I3(II))   + Z_RR_GMLT(II)  * ZMAXTIME(II)
            ZTOT_CR_GMLT(I1(II),I2(II),I3(II)) =   ZTOT_CR_GMLT(I1(II),I2(II),I3(II))   + Z_CR_GMLT(II)  * ZMAXTIME(II)
!!$            ZTOT_RC_WETH(I1(II),I2(II),I3(II)) =   ZTOT_RC_WETH(I1(II),I2(II),I3(II))   + Z_RC_WETH(II)  * ZMAXTIME(II)
!!$            ZTOT_CC_WETH(I1(II),I2(II),I3(II)) =   ZTOT_CC_WETH(I1(II),I2(II),I3(II))   + Z_CC_WETH(II)  * ZMAXTIME(II)
!!$            ZTOT_RR_WETH(I1(II),I2(II),I3(II)) =   ZTOT_RR_WETH(I1(II),I2(II),I3(II))   + Z_RR_WETH(II)  * ZMAXTIME(II)
!!$            ZTOT_CR_WETH(I1(II),I2(II),I3(II)) =   ZTOT_CR_WETH(I1(II),I2(II),I3(II))   + Z_CR_WETH(II)  * ZMAXTIME(II)
!!$            ZTOT_RI_WETH(I1(II),I2(II),I3(II)) =   ZTOT_RI_WETH(I1(II),I2(II),I3(II))   + Z_RI_WETH(II)  * ZMAXTIME(II)
!!$            ZTOT_CI_WETH(I1(II),I2(II),I3(II)) =   ZTOT_CI_WETH(I1(II),I2(II),I3(II))   + Z_CI_WETH(II)  * ZMAXTIME(II)
!!$            ZTOT_RS_WETH(I1(II),I2(II),I3(II)) =   ZTOT_RS_WETH(I1(II),I2(II),I3(II))   + Z_RS_WETH(II)  * ZMAXTIME(II)
!!$            ZTOT_RG_WETH(I1(II),I2(II),I3(II)) =   ZTOT_RG_WETH(I1(II),I2(II),I3(II))   + Z_RG_WETH(II)  * ZMAXTIME(II)
!!$            ZTOT_RH_WETH(I1(II),I2(II),I3(II)) =   ZTOT_RH_WETH(I1(II),I2(II),I3(II))   + Z_RH_WETH(II)  * ZMAXTIME(II)
!!$            ZTOT_RG_COHG(I1(II),I2(II),I3(II)) =   ZTOT_RG_COHG(I1(II),I2(II),I3(II))   + Z_RG_COHG(II)  * ZMAXTIME(II)
!!$            ZTOT_RR_HMLT(I1(II),I2(II),I3(II)) =   ZTOT_RR_HMLT(I1(II),I2(II),I3(II))   + Z_RR_HMLT(II)  * ZMAXTIME(II)
!!$            ZTOT_CR_HMLT(I1(II),I2(II),I3(II)) =   ZTOT_CR_HMLT(I1(II),I2(II),I3(II))   + Z_CR_HMLT(II)  * ZMAXTIME(II)

          !Correction term
          ZTOT_RV_CORR2(I1(II),I2(II),I3(II)) =   ZTOT_RV_CORR2(I1(II),I2(II),I3(II)) + Z_RV_CORR2(II)
          ZTOT_RC_CORR2(I1(II),I2(II),I3(II)) =   ZTOT_RC_CORR2(I1(II),I2(II),I3(II)) + Z_RC_CORR2(II)
          ZTOT_RR_CORR2(I1(II),I2(II),I3(II)) =   ZTOT_RR_CORR2(I1(II),I2(II),I3(II)) + Z_RR_CORR2(II)
          ZTOT_RI_CORR2(I1(II),I2(II),I3(II)) =   ZTOT_RI_CORR2(I1(II),I2(II),I3(II)) + Z_RI_CORR2(II)
          ZTOT_CC_CORR2(I1(II),I2(II),I3(II)) =   ZTOT_CC_CORR2(I1(II),I2(II),I3(II)) + Z_CC_CORR2(II)
          ZTOT_CR_CORR2(I1(II),I2(II),I3(II)) =   ZTOT_CR_CORR2(I1(II),I2(II),I3(II)) + Z_CR_CORR2(II)
          ZTOT_CI_CORR2(I1(II),I2(II),I3(II)) =   ZTOT_CI_CORR2(I1(II),I2(II),I3(II)) + Z_CI_CORR2(II)
         END DO
      ENDIF
      !
      ! Deallocating variables
      !
      DEALLOCATE(I1)
      DEALLOCATE(I2)
      DEALLOCATE(I3)
      DEALLOCATE(ZRHODREF1D)
      DEALLOCATE(ZEXNREF1D)
      DEALLOCATE(ZEXN1D)
      DEALLOCATE(ZP1D)     
      DEALLOCATE(ZTHT1D)
      DEALLOCATE(ZRVT1D)
      DEALLOCATE(ZRCT1D)
      DEALLOCATE(ZRRT1D)
      DEALLOCATE(ZRIT1D)
      DEALLOCATE(ZRST1D)
      DEALLOCATE(ZRGT1D)
      DEALLOCATE(ZRHT1D)
      DEALLOCATE(ZCCT1D)
      DEALLOCATE(ZCRT1D)
      DEALLOCATE(ZCIT1D)
      DEALLOCATE(ZIFNN1D)
      DEALLOCATE(ZEVAP1D)
      DEALLOCATE(ZTIME1D)
      DEALLOCATE(LLCOMPUTE1D)
      DEALLOCATE(IITER1D)
      DEALLOCATE(ZTIME_LASTCALL1D)
      DEALLOCATE(Z0RVT1D)
      DEALLOCATE(Z0RCT1D)
      DEALLOCATE(Z0RRT1D)
      DEALLOCATE(Z0RIT1D)
      DEALLOCATE(Z0RST1D)
      DEALLOCATE(Z0RGT1D)
      DEALLOCATE(Z0RHT1D)
      DEALLOCATE(ZCF1D)
      DEALLOCATE(ZIF1D)
      DEALLOCATE(ZPF1D)
      !
      DEALLOCATE(ZMAXTIME)
      DEALLOCATE(ZTIME_THRESHOLD)
      !
      DEALLOCATE(ZA_TH)
      DEALLOCATE(ZA_RV)
      DEALLOCATE(ZA_RC)
      DEALLOCATE(ZA_RR)
      DEALLOCATE(ZA_RI)
      DEALLOCATE(ZA_RS)
      DEALLOCATE(ZA_RG)
      DEALLOCATE(ZA_RH)
      DEALLOCATE(ZA_CC)
      DEALLOCATE(ZA_CR)
      DEALLOCATE(ZA_CI)
      !
      DEALLOCATE(ZB_TH)
      DEALLOCATE(ZB_RV)   
      DEALLOCATE(ZB_RC)   
      DEALLOCATE(ZB_RR)   
      DEALLOCATE(ZB_RI)   
      DEALLOCATE(ZB_RS)   
      DEALLOCATE(ZB_RG)   
      DEALLOCATE(ZB_RH)   
      DEALLOCATE(ZB_CC)   
      DEALLOCATE(ZB_CR)  
      DEALLOCATE(ZB_CI)  
      DEALLOCATE(ZB_IFNN)
      !
      DEALLOCATE(Z_CR_BRKU)
      DEALLOCATE(Z_TH_HONR)
      DEALLOCATE(Z_RR_HONR)
      DEALLOCATE(Z_CR_HONR)
      DEALLOCATE(Z_TH_IMLT)
      DEALLOCATE(Z_RC_IMLT)
      DEALLOCATE(Z_CC_IMLT)
      DEALLOCATE(Z_TH_HONC)
      DEALLOCATE(Z_RC_HONC)
      DEALLOCATE(Z_CC_HONC)
      DEALLOCATE(Z_CC_SELF) 
      DEALLOCATE(Z_RC_AUTO) 
      DEALLOCATE(Z_CC_AUTO)
      DEALLOCATE(Z_CR_AUTO) 
      DEALLOCATE(Z_RC_ACCR) 
      DEALLOCATE(Z_CC_ACCR) 
      DEALLOCATE(Z_CR_SCBU)
      DEALLOCATE(Z_TH_EVAP) 
      DEALLOCATE(Z_RR_EVAP) 
      DEALLOCATE(Z_RI_CNVI)
      DEALLOCATE(Z_CI_CNVI)
      DEALLOCATE(Z_TH_DEPS)
      DEALLOCATE(Z_RS_DEPS)
      DEALLOCATE(Z_TH_DEPI)
      DEALLOCATE(Z_RI_DEPI)
      DEALLOCATE(Z_RI_CNVS)
      DEALLOCATE(Z_CI_CNVS)
      DEALLOCATE(Z_RI_AGGS) 
      DEALLOCATE(Z_CI_AGGS) 
      DEALLOCATE(Z_TH_DEPG) 
      DEALLOCATE(Z_RG_DEPG) 
      DEALLOCATE(Z_TH_BERFI)
      DEALLOCATE(Z_RC_BERFI)
      DEALLOCATE(Z_TH_RIM) 
      DEALLOCATE(Z_RC_RIM) 
      DEALLOCATE(Z_CC_RIM)  
      DEALLOCATE(Z_RS_RIM) 
      DEALLOCATE(Z_RG_RIM) 
      DEALLOCATE(Z_RI_HMS) 
      DEALLOCATE(Z_CI_HMS) 
      DEALLOCATE(Z_RS_HMS) 
      DEALLOCATE(Z_TH_ACC) 
      DEALLOCATE(Z_RR_ACC) 
      DEALLOCATE(Z_CR_ACC) 
      DEALLOCATE(Z_RS_ACC) 
      DEALLOCATE(Z_RG_ACC)  
      DEALLOCATE(Z_RS_CMEL) 
      DEALLOCATE(Z_TH_CFRZ)
      DEALLOCATE(Z_RR_CFRZ)
      DEALLOCATE(Z_CR_CFRZ)
      DEALLOCATE(Z_RI_CFRZ)
      DEALLOCATE(Z_CI_CFRZ)
      DEALLOCATE(Z_TH_WETG)
      DEALLOCATE(Z_RC_WETG)
      DEALLOCATE(Z_CC_WETG)
      DEALLOCATE(Z_RR_WETG) 
      DEALLOCATE(Z_CR_WETG) 
      DEALLOCATE(Z_RI_WETG)
      DEALLOCATE(Z_CI_WETG)
      DEALLOCATE(Z_RS_WETG)
      DEALLOCATE(Z_RG_WETG)
      DEALLOCATE(Z_RH_WETG) 
      DEALLOCATE(Z_TH_DRYG) 
      DEALLOCATE(Z_RC_DRYG) 
      DEALLOCATE(Z_CC_DRYG)
      DEALLOCATE(Z_RR_DRYG)
      DEALLOCATE(Z_CR_DRYG)
      DEALLOCATE(Z_RI_DRYG)
      DEALLOCATE(Z_CI_DRYG)
      DEALLOCATE(Z_RS_DRYG) 
      DEALLOCATE(Z_RG_DRYG)
      DEALLOCATE(Z_RI_HMG) 
      DEALLOCATE(Z_CI_HMG) 
      DEALLOCATE(Z_RG_HMG) 
      DEALLOCATE(Z_TH_GMLT)
      DEALLOCATE(Z_RR_GMLT)
      DEALLOCATE(Z_CR_GMLT)

      DEALLOCATE(Z_RV_CORR2)
      DEALLOCATE(Z_RC_CORR2)
      DEALLOCATE(Z_RR_CORR2)
      DEALLOCATE(Z_RI_CORR2)
      DEALLOCATE(Z_CC_CORR2)
      DEALLOCATE(Z_CR_CORR2)
      DEALLOCATE(Z_CI_CORR2)
      !
   ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       7.     TOTAL TENDENCIES
!               ----------------
!
! Source at the end of microphysics = new state / PTSTEP
!
PTHS(:,:,:) = ZTHT(:,:,:) * ZINV_TSTEP
!
PRS(:,:,:,1) = ZRVT(:,:,:) *ZINV_TSTEP
IF ( KRR .GE. 2 ) PRS(:,:,:,2) = ZRCT(:,:,:) *ZINV_TSTEP
IF ( KRR .GE. 3 ) PRS(:,:,:,3) = ZRRT(:,:,:) *ZINV_TSTEP
IF ( KRR .GE. 4 ) PRS(:,:,:,4) = ZRIT(:,:,:) *ZINV_TSTEP
IF ( KRR .GE. 5 ) PRS(:,:,:,5) = ZRST(:,:,:) *ZINV_TSTEP
IF ( KRR .GE. 6 ) PRS(:,:,:,6) = ZRGT(:,:,:) *ZINV_TSTEP
IF ( KRR .GE. 7 ) PRS(:,:,:,7) = ZRHT(:,:,:) *ZINV_TSTEP
!
IF ( LWARM )             PSVS(:,:,:,NSV_LIMA_NC) = ZCCT(:,:,:) *ZINV_TSTEP
IF ( LWARM .AND. LRAIN ) PSVS(:,:,:,NSV_LIMA_NR) = ZCRT(:,:,:) *ZINV_TSTEP
IF ( LCOLD )             PSVS(:,:,:,NSV_LIMA_NI) = ZCIT(:,:,:) *ZINV_TSTEP
!
IF ( NMOD_CCN .GE. 1 )   PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1) = ZCCNFT(:,:,:,:) *ZINV_TSTEP
IF ( NMOD_CCN .GE. 1 )   PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1) = ZCCNAT(:,:,:,:) *ZINV_TSTEP
IF ( NMOD_IFN .GE. 1 )   PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1) = ZIFNFT(:,:,:,:) *ZINV_TSTEP
IF ( NMOD_IFN .GE. 1 )   PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1) = ZIFNNT(:,:,:,:) *ZINV_TSTEP
IF ( NMOD_IMM .GE. 1 )   PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1) = ZIMMNT(:,:,:,:) *ZINV_TSTEP
IF ( LCOLD .AND. LHHONI) PSVS(:,:,:,NSV_LIMA_HOM_HAZE) = ZHOMFT(:,:,:) *ZINV_TSTEP
!
!
!
! Call budgets
!
if ( lbu_enable ) then
  allocate( zrhodjontstep(size( prhodj, 1), size( prhodj, 2), size( prhodj, 3) ) )
  zrhodjontstep(:, :, :) = zinv_tstep * prhodj(:, :, :)

  if ( lbudget_th ) then
    call Budget_store_add( tbudgets(NBUDGET_TH), 'REVA',  ztot_th_evap (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'HONC',  ztot_th_honc (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'HONR',  ztot_th_honr (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'DEPS',  ztot_th_deps (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'DEPI',  ztot_th_depi (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'DEPG',  ztot_th_depg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'IMLT',  ztot_th_imlt (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'BERFI', ztot_th_berfi(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'RIM',   ztot_th_rim  (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'ACC',   ztot_th_acc  (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'CFRZ',  ztot_th_cfrz (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'WETG',  ztot_th_wetg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'DRYG',  ztot_th_dryg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_TH), 'GMLT',  ztot_th_gmlt (:, :, :) * zrhodjontstep(:, :, :) )
  end if

  if ( lbudget_rv ) then
    call Budget_store_add( tbudgets(NBUDGET_RV), 'REVA', -ztot_rr_evap (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RV), 'DEPS', -ztot_rs_deps (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RV), 'DEPI', -ztot_ri_depi (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RV), 'DEPG', -ztot_rg_depg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RV), 'CORR2', ztot_rv_corr2(:, :, :) * zrhodjontstep(:, :, :) )
  end if

  if ( lbudget_rc ) then
    call Budget_store_add( tbudgets(NBUDGET_RC), 'AUTO',  ztot_rc_auto (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RC), 'ACCR',  ztot_rc_accr (:, :, :) * zrhodjontstep(:, :, :) )
    !call Budget_store_add( tbudgets(NBUDGET_RC), 'REVA',  0. )
    call Budget_store_add( tbudgets(NBUDGET_RC), 'HONC',  ztot_rc_honc (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RC), 'IMLT',  ztot_rc_imlt (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RC), 'BERFI', ztot_rc_berfi(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RC), 'RIM',   ztot_rc_rim  (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RC), 'WETG',  ztot_rc_wetg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RC), 'DRYG',  ztot_rc_dryg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RC), 'CVRC', -ztot_rr_cvrc (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RC), 'CORR2', ztot_rc_corr2(:, :, :) * zrhodjontstep(:, :, :) )
  end if

  if ( lbudget_rr ) then
    call Budget_store_add( tbudgets(NBUDGET_RR), 'AUTO', -ztot_rc_auto(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RR), 'ACCR', -ztot_rc_accr(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RR), 'REVA',  ztot_rr_evap(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RR), 'HONR',  ztot_rr_honr(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RR), 'ACC',   ztot_rr_acc (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RR), 'CFRZ',  ztot_rr_cfrz(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RR), 'WETG',  ztot_rr_wetg(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RR), 'DRYG',  ztot_rr_dryg(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RR), 'GMLT',  ztot_rr_gmlt(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RR), 'CVRC',  ztot_rr_cvrc(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RR), 'CORR2', ztot_rr_corr2(:, :, :) * zrhodjontstep(:, :, :) )
  end if

  if ( lbudget_ri ) then
    call Budget_store_add( tbudgets(NBUDGET_RI), 'HONC',  -ztot_rc_honc (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'CNVI',   ztot_ri_cnvi (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'CNVS',   ztot_ri_cnvs (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'AGGS',   ztot_ri_aggs (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'IMLT',  -ztot_rc_imlt (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'BERFI', -ztot_rc_berfi(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'HMS',    ztot_ri_hms  (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'CFRZ',   ztot_ri_cfrz (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'DEPI',   ztot_ri_depi (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'WETG',   ztot_ri_wetg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'DRYG',   ztot_ri_dryg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'HMG',    ztot_ri_hmg  (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RI), 'CORR2',  ztot_ri_corr2(:, :, :) * zrhodjontstep(:, :, :) )
  end if

  if ( lbudget_rs ) then
    call Budget_store_add( tbudgets(NBUDGET_RS), 'CNVI', -ztot_ri_cnvi(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RS), 'DEPS',  ztot_rs_deps(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RS), 'CNVS', -ztot_ri_cnvs(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RS), 'AGGS', -ztot_ri_aggs(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RS), 'RIM',   ztot_rs_rim (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RS), 'HMS',   ztot_rs_hms (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RS), 'ACC',   ztot_rs_acc (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RS), 'CMEL',  ztot_rs_cmel(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RS), 'WETG',  ztot_rs_wetg(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RS), 'DRYG',  ztot_rs_dryg(:, :, :) * zrhodjontstep(:, :, :) )
  end if

  if ( lbudget_rg ) then
    call Budget_store_add( tbudgets(NBUDGET_RG), 'HONR', -ztot_rr_honr(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RG), 'DEPG',  ztot_rg_depg(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RG), 'RIM',   ztot_rg_rim (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RG), 'ACC',   ztot_rg_acc (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RG), 'CMEL', -ztot_rs_cmel(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RG), 'CFRZ', ( -ztot_rr_cfrz(:, :, :) - ztot_ri_cfrz(:, :, :) ) &
                                                         * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RG), 'WETG',  ztot_rg_wetg(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RG), 'DRYG',  ztot_rg_dryg(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RG), 'HMG',   ztot_rg_hmg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(NBUDGET_RG), 'GMLT', -ztot_rr_gmlt(:, :, :) * zrhodjontstep(:, :, :) )
  end if

  if ( lbudget_rh ) then
    call Budget_store_add( tbudgets(NBUDGET_RH), 'WETG', ztot_rh_wetg(:, :, :) * zrhodjontstep(:, :, :) )
  end if

  if ( lbudget_sv ) then
    !
    ! Cloud droplets
    !
    idx = NBUDGET_SV1 - 1 + nsv_lima_nc
    call Budget_store_add( tbudgets(idx), 'SELF',  ztot_cc_self (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'AUTO',  ztot_cc_auto (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'ACCR',  ztot_cc_accr (:, :, :) * zrhodjontstep(:, :, :) )
    !call Budget_store_add( tbudgets(idx), 'REVA',  0. )
    call Budget_store_add( tbudgets(idx), 'HONC',  ztot_cc_honc (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'IMLT',  ztot_cc_imlt (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'RIM',   ztot_cc_rim  (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'WETG',  ztot_cc_wetg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'DRYG',  ztot_cc_dryg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'CVRC', -ztot_cr_cvrc (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'CORR2', ztot_cc_corr2(:, :, :) * zrhodjontstep(:, :, :) )
    !
    ! Rain drops
    !
    idx = NBUDGET_SV1 - 1 + nsv_lima_nr
    call Budget_store_add( tbudgets(idx), 'AUTO',  ztot_cr_auto(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'SCBU',  ztot_cr_scbu(:, :, :) * zrhodjontstep(:, :, :) )
    !call Budget_store_add( tbudgets(idx), 'REVA',  0. )
    call Budget_store_add( tbudgets(idx), 'BRKU',  ztot_cr_brku(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'HONR',  ztot_cr_honr(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'ACC',   ztot_cr_acc (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'CFRZ',  ztot_cr_cfrz(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'WETG',  ztot_cr_wetg(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'DRYG',  ztot_cr_dryg(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'GMLT',  ztot_cr_gmlt(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'CVRC',  ztot_cr_cvrc(:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'CORR2', ztot_cr_corr2(:, :, :) * zrhodjontstep(:, :, :) )
    !
    ! Ice crystals
    !
    idx = NBUDGET_SV1 - 1 + nsv_lima_ni
    call Budget_store_add( tbudgets(idx), 'HONC',  -ztot_cc_honc (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'CNVI',   ztot_ci_cnvi (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'CNVS',   ztot_ci_cnvs (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'AGGS',   ztot_ci_aggs (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'IMLT',  -ztot_cc_imlt (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'HMS',    ztot_ci_hms  (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'CFRZ',   ztot_ci_cfrz (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'WETG',   ztot_ci_wetg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'DRYG',   ztot_ci_dryg (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'HMG',    ztot_ci_hmg  (:, :, :) * zrhodjontstep(:, :, :) )
    call Budget_store_add( tbudgets(idx), 'CORR2',  ztot_ci_corr2(:, :, :) * zrhodjontstep(:, :, :) )

    do ii = 1, nmod_ifn
      idx = nsv_lima_ifn_nucl + ii - 1
      call Budget_store_add( tbudgets(idx), 'IMLT', ztot_ifnn_imlt(:, :, :, ii) * zrhodjontstep(:, :, :) )
    end do
  end if

  deallocate( zrhodjontstep )
end if
!
END SUBROUTINE LIMA
