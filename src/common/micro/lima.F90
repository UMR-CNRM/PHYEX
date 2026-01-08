!MNH_LIC Copyright 2013-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #####################################################################
SUBROUTINE LIMA ( LIMAP, LIMAW, LIMAC, LIMAM, TNSV, D, CST, NEBN,         &
                  ICED, ICEP, ELECD, ELECP,                               &
                  BUCONF, TBUDGETS, HACTCCN, KBUDGETS, KRR,               &
                  PTSTEP, OELEC,                                          &
                  PRHODREF, PEXNREF, PDZZ, PTHVREFZIKB,                   &
                  PRHODJ, PPABST,                                         &
                  KCARB, KSOA, KSP, ODUST, OSALT, OORILAM,                &
                  ODTHRAD, PDTHRAD, PTHT, PRT, PSVT, PCIT, PW_NU,         &
                  PAERO,PSOLORG, PMI, PTHS, PRS, PSVS,                    &
                  PINPRC, PINDEP, PINPRR, PINPRI, PINPRS, PINPRG, PINPRH, &
                  PEVAP3D, PCLDFR, PICEFR, PPRCFR, PFPR,                  &
                  PHLC_HCF, PHLC_HRC,                                     &
                  PHLI_HCF, PHLI_HRI,                                     &
                  PLATHAM_IAGGS, PEFIELDW, PSV_ELEC_T, PSV_ELEC_S         )
!     #####################################################################
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
!  C. Barthe      04/2022: add cloud electrification
!  C. Barthe      03/2023: add CIBU, RDSF and 2 moments for s, g and h in cloud electrification
!
!-----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_PARAM_LIMA_MIXED,ONLY:PARAM_LIMA_MIXED_T
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
USE MODD_PARAM_LIMA,      ONLY:PARAM_LIMA_T
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
USE MODD_RAIN_ICE_DESCR_N,ONLY: RAIN_ICE_DESCR_T
USE MODD_RAIN_ICE_PARAM_N,ONLY: RAIN_ICE_PARAM_T
USE MODD_ELEC_PARAM,      ONLY: ELEC_PARAM_T
USE MODD_ELEC_DESCR,      ONLY: ELEC_DESCR_T
USE MODD_BUDGET,          ONLY: TBUDGETDATA_PTR, TBUDGETCONF_T, NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, &
                                NBUDGET_RI, NBUDGET_RR, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, NBUDGET_SV1
USE MODD_CST,             ONLY: CST_T
USE MODD_NSV,             ONLY: NSV_T
USE MODD_NEB_N,           ONLY: NEB_T
USE MODE_TOOLS,           only: COUNTJV

USE MODE_LIMA_COMPUTE_PDF, ONLY: LIMA_COMPUTE_PDF
USE MODE_LIMA_RAINFR_VERT, ONLY: LIMA_RAINFR_VERT
USE MODE_LIMA_COMPUTE_CLOUD_FRACTIONS, ONLY: LIMA_COMPUTE_CLOUD_FRACTIONS
USE MODE_LIMA_SHAPE_COMPUTE_LBDA, ONLY : LIMA_SHAPE_COMPUTE_LBDA_3D
USE MODE_LIMA_INST_PROCS, ONLY: LIMA_INST_PROCS
USE MODE_LIMA_NUCLEATION_PROCS, ONLY: LIMA_NUCLEATION_PROCS
USE MODE_LIMA_SEDIMENTATION, ONLY: LIMA_SEDIMENTATION
USE MODE_LIMA_TENDENCIES, ONLY: LIMA_TENDENCIES
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
USE MODE_ELEC_TENDENCIES, ONLY : ELEC_TENDENCIES
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_MIXED_T),INTENT(IN)::LIMAM
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(NSV_T),              INTENT(IN)    :: TNSV
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
TYPE(CST_T),              INTENT(IN)    :: CST
TYPE(NEB_T),              INTENT(IN)    :: NEBN
TYPE(RAIN_ICE_DESCR_T),   INTENT(IN)    :: ICED
TYPE(RAIN_ICE_PARAM_T),   INTENT(IN)    :: ICEP
TYPE(ELEC_PARAM_T),       INTENT(IN)    :: ELECP   ! electrical parameters
TYPE(ELEC_DESCR_T),       INTENT(IN)    :: ELECD   ! electrical descriptive csts
TYPE(TBUDGETCONF_T),      INTENT(IN)    :: BUCONF
CHARACTER(LEN=4),         INTENT(IN)    :: HACTCCN  ! kind of CCN activation
TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER,                  INTENT(IN)    :: KBUDGETS
INTEGER,                  INTENT(IN)    :: KRR
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Time step
!
LOGICAL,                  INTENT(IN)    :: OELEC      ! if true, cloud electrification is activated
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PDZZ       ! Layer thikness (m)
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PRHODJ     ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PPABST     ! absolute pressure at t
INTEGER,                  INTENT(IN)    :: KCARB, KSOA, KSP ! for array size declarations
LOGICAL,                  INTENT(IN)    :: ODUST, OSALT, OORILAM
!
LOGICAL,                                 INTENT(IN)   :: ODTHRAD    ! Use radiative temperature tendency
REAL, DIMENSION(MERGE(D%NIJT,0,ODTHRAD), &
                MERGE(D%NKT,0,ODTHRAD)),   INTENT(IN) :: PDTHRAD   ! dT/dt due to radiation
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PTHT       ! Theta at time t
REAL, DIMENSION(D%NIJT, D%NKT, KRR), INTENT(IN) :: PRT        ! Mixing ratios at time t
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), INTENT(IN) :: PSVT       ! Concentrations at time t
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(INOUT)    :: PCIT       ! Theta at time t
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PW_NU      ! w for CCN activation
REAL, DIMENSION(D%NIJT, D%NKT ,TNSV%NSV), INTENT(INOUT) :: PAERO    ! Aerosol concentration
REAL, DIMENSION(D%NIJT, D%NKT, 10),  INTENT(IN)    :: PSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(D%NIJT, D%NKT, KSP+KCARB+KSOA), INTENT(IN)    :: PMI
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(INOUT)    :: PTHS       ! Theta source
REAL, DIMENSION(D%NIJT, D%NKT, KRR), INTENT(INOUT) :: PRS        ! Mixing ratios sources
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), INTENT(INOUT) :: PSVS       ! Concentration sources
!
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRC     ! Cloud instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINDEP     ! Cloud droplets deposition
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRR     ! Rain instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRI     ! Rain instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRS     ! Snow instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRG     ! Graupel instant precip
REAL, DIMENSION(D%NIJT),     INTENT(OUT)        :: PINPRH     ! Rain instant precip
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(OUT)   :: PEVAP3D    ! Rain evap profile
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(INOUT) :: PCLDFR     ! Cloud fraction
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(INOUT) :: PICEFR     ! Cloud fraction
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(INOUT) :: PPRCFR     ! Cloud fraction
REAL, DIMENSION(D%NIJT, D%NKT, KRR), INTENT(OUT) :: PFPR    ! Precipitation fluxes in altitude
!
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL,   INTENT(INOUT) :: PHLC_HCF
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL,   INTENT(INOUT) :: PHLC_HRC
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL,   INTENT(INOUT) :: PHLI_HCF
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL,   INTENT(INOUT) :: PHLI_HRI
!
REAL, DIMENSION(D%NIJT, D%NKT),   OPTIONAL, INTENT(IN)       :: PLATHAM_IAGGS  ! Factor for IAGGS modification due to Efield
REAL, DIMENSION(D%NIJT, D%NKT),   OPTIONAL, INTENT(IN)       :: PEFIELDW   ! Vertical component of the electric field
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), OPTIONAL, INTENT(IN)    :: PSV_ELEC_T ! Charge density at time t
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), OPTIONAL, INTENT(INOUT) :: PSV_ELEC_S ! Charge density sources
!
REAL, INTENT(IN)                :: PTHVREFZIKB ! Reference thv at IKB for electricity
!*       0.2   Declarations of local variables :
!
!
! Prognostic variables and sources
REAL, DIMENSION(D%NIJT,D%NKT)      :: ZTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT, ZRHT
REAL, DIMENSION(D%NIJT,D%NKT)      :: ZCCT, ZCRT, ZCIT, ZCST, ZCGT, ZCHT
REAL, DIMENSION(D%NIJT,D%NKT)      :: ZTHS, ZRVS, ZRCS, ZRRS, ZRIS, ZRSS, ZRGS, ZRHS
REAL, DIMENSION(D%NIJT,D%NKT)      :: ZCCS, ZCRS, ZCIS, ZCSS, ZCGS, ZCHS
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_CCN) :: ZCCNFT, ZCCNAT
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_CCN) :: ZCCNFS, ZCCNAS
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_IFN) :: ZIFNFT, ZIFNNT
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_IFN) :: ZIFNFS, ZIFNNS
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_IMM) :: ZIMMNT
REAL, DIMENSION(D%NIJT,D%NKT,LIMAP%NMOD_IMM) :: ZIMMNS
REAL, DIMENSION(D%NIJT,D%NKT)      :: ZHOMFT
REAL, DIMENSION(D%NIJT,D%NKT)      :: ZHOMFS

REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCIS_SHAPE  ! Nb concentration for each ice habit (source)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZRIS_SHAPE  ! Mixing ratio for each ice habit (source)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCIT_SHAPE  ! Nb concentration for each ice habit (at t)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZRIT_SHAPE  ! Mixing ratio for each ice habit (at t)
!
! Other 3D thermodynamical variables
REAL, DIMENSION(D%NIJT,D%NKT)      :: ZEXN, ZT, ZLSFACT, ZW, ZTMP

!
! Packed prognostic & thermo variables
REAL, DIMENSION(:),   ALLOCATABLE ::                               &
     ZP1D, ZRHODREF1D, ZEXNREF1D, ZEXN1D,                          &
     ZTHT1D,                                                       &
     ZRVT1D, ZRCT1D, ZRRT1D, ZRIT1D, ZRST1D, ZRGT1D, ZRHT1D,       &
     ZCCT1D, ZCRT1D, ZCIT1D, ZCST1D, ZCGT1D, ZCHT1D,               &
     ZEVAP1D,                                                      &
     ZHLC_LCF1D, ZHLC_LRC1D, ZHLI_LCF1D, ZHLI_LRI1D,               &
     ZHLC_HCF1D, ZHLC_HRC1D, ZHLI_HCF1D, ZHLI_HRI1D

REAL, DIMENSION(:,:), ALLOCATABLE :: ZIFNN1D
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCIT1D_SHAPE

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
     Z_TH_EVAP, Z_RR_EVAP, Z_CR_EVAP,                       & ! evaporation of rain drops (EVAP) : rv=-rr-rc, rc, Nc, rr, Nr, th
     Z_RI_CNVI, Z_CI_CNVI,                                  & ! conversion snow -> ice (CNVI) : ri, Ni, rs=-ri
     Z_TH_DEPS, Z_RS_DEPS,                                  & ! deposition of vapor on snow (DEPS) : rv=-rs, rs, th
     Z_TH_DEPI, Z_RI_DEPI,                                  & ! deposition of vapor on ice (DEPI) : rv=-ri, ri, th
     Z_RI_CNVS, Z_CI_CNVS,                                  & ! conversion ice -> snow (CNVS) : ri, Ni, rs=-ri
     Z_CS_SSC,                                              & ! self collection of snow (SSC) : Ns
     Z_CI_ISC,                                              & ! self collection of ice (ISC) : Ni      ++mt--
     Z_RI_AGGS, Z_CI_AGGS,                                  & ! aggregation of ice on snow (AGGS) : ri, Ni, rs=-ri
     Z_TH_DEPG, Z_RG_DEPG,                                  & ! deposition of vapor on graupel (DEPG) : rv=-rg, rg, th
     Z_TH_BERFI, Z_RC_BERFI,                                & ! Bergeron (BERFI) : rc, ri=-rc, th
     Z_TH_RIM, Z_CC_RIM, Z_CS_RIM, Z_RC_RIMSS, Z_RC_RIMSG, Z_RS_RIMCG, & ! cloud droplet riming (RIM) : rc, Nc, rs, rg, th
     Z_RI_HMS, Z_CI_HMS, Z_RS_HMS,                          & ! hallett mossop snow (HMS) : ri, Ni, rs
     Z_TH_ACC, Z_CR_ACC, Z_CS_ACC, Z_RR_ACCSS, Z_RR_ACCSG, Z_RS_ACCRG, & ! rain accretion on aggregates (ACC) : rr, Nr, rs, rg, th
     Z_RS_CMEL, Z_CS_CMEL,                                  & ! conversion-melting (CMEL) : rs, rg=-rs
     Z_TH_CFRZ, Z_RR_CFRZ, Z_CR_CFRZ, Z_RI_CFRZ, Z_CI_CFRZ, & ! rain freezing (CFRZ) : rr, Nr, ri, Ni, rg=-rr-ri, th
     Z_RI_CIBU, Z_CI_CIBU,                                  & ! collisional ice break-up (CIBU) : ri, Ni, rs=-ri
     Z_RI_RDSF, Z_CI_RDSF,                                  & ! rain drops freezing shattering (RDSF) : ri, Ni, rg=-ri
     Z_TH_WETG, Z_RC_WETG, Z_CC_WETG, Z_RR_WETG, Z_CR_WETG, & ! wet growth of graupel (WETG) : rc, NC, rr, Nr, ri, Ni, rs, rg, rh, th
     Z_RI_WETG, Z_CI_WETG, Z_RS_WETG, Z_CS_WETG, Z_RG_WETG, Z_CG_WETG, Z_RH_WETG, & ! wet growth of graupel (WETG) : rc, NC, rr, Nr, ri, Ni, rs, Ns, rg, Ng, rh, Nh=-Ng, th
     Z_TH_DRYG, Z_RC_DRYG, Z_CC_DRYG, Z_RR_DRYG, Z_CR_DRYG, & ! dry growth of graupel (DRYG) : rc, Nc, rr, Nr, ri, Ni, rs, rg, th
     Z_RI_DRYG, Z_CI_DRYG, Z_RS_DRYG, Z_CS_DRYG, Z_RG_DRYG, & ! dry growth of graupel (DRYG) : rc, Nc, rr, Nr, ri, Ni, rs, Ns, rg, th
     Z_RI_HMG, Z_CI_HMG, Z_RG_HMG,                          & ! hallett mossop graupel (HMG) : ri, Ni, rg
     Z_TH_GMLT, Z_RR_GMLT, Z_CR_GMLT, Z_CG_GMLT,            & ! graupel melting (GMLT) : rr, Nr, rg=-rr, th
     Z_TH_DEPH, Z_RH_DEPH,                                  & ! deposition of vapor on hail (DEPH) : rv=-rh, rh, th
     Z_TH_WETH, Z_RC_WETH, Z_CC_WETH, Z_RR_WETH, Z_CR_WETH, & ! wet growth of hail (WETH) : rc, Nc, rr, Nr, ri, Ni, rs, rg, rh, th
     Z_RI_WETH, Z_CI_WETH, Z_RS_WETH, Z_CS_WETH, Z_RG_WETH, Z_CG_WETH, Z_RH_WETH, & ! wet growth of hail (WETH) : rc, Nc, rr, Nr, ri, Ni, rs, rg, rh, th
     Z_RG_COHG, Z_CG_COHG,                                  & ! conversion of hail into graupel (COHG) : rg, rh
     Z_TH_HMLT, Z_RR_HMLT, Z_CR_HMLT, Z_CH_HMLT,            & ! hail melting (HMLT) : rr, Nr, rh=-rr, th
     Z_RV_CORR2, Z_RC_CORR2, Z_RR_CORR2, Z_RI_CORR2,        &
     Z_CC_CORR2, Z_CR_CORR2, Z_CI_CORR2
!--cb--
!
! for the conversion from rain to cloud, we need a 3D variable instead of a 1D packed variable
REAL, DIMENSION(D%NIJT,D%NKT) ::  Z_RR_CVRC, Z_CR_CVRC        ! conversion of rain into cloud droplets (CVRC)

!
! Packed variables for total tendencies
REAL, DIMENSION(:),   ALLOCATABLE ::                                              &
     ZA_TH, ZA_RV, ZA_RC, ZA_CC, ZA_RR, ZA_CR, ZA_RI, ZA_CI, ZA_RS, ZA_CS, ZA_RG, ZA_CG, ZA_RH, ZA_CH, & ! ZA = continuous tendencies (kg/kg/s = S variable)
     ZB_TH, ZB_RV, ZB_RC, ZB_CC, ZB_RR, ZB_CR, ZB_RI, ZB_CI, ZB_RS, ZB_CS, ZB_RG, ZB_CG, ZB_RH, ZB_CH    ! ZB = instant mixing ratio change (kg/kg = T variable)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZB_IFNN

REAL, DIMENSION(:,:), ALLOCATABLE :: ZA_CI_SHAPE ! continuous tendencies (kg/kg/s = S variable) for ice shapes C.
REAL, DIMENSION(:,:), ALLOCATABLE :: ZB_CI_SHAPE ! instant change (kg/kg = T variable) for ice shapes C.
!
! for each process & species, we need 3D variables to store total mmr and conc change (kg/kg and #/kg and theta)
REAL, DIMENSION(:,:), ALLOCATABLE ::                                       &
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
     ZTOT_TH_EVAP, ZTOT_RR_EVAP, ZTOT_CR_EVAP,                             & ! evaporation of rain drops (EVAP)
     ZTOT_RI_CNVI, ZTOT_CI_CNVI,                                           & ! conversion snow -> ice (CNVI)
     ZTOT_TH_DEPS, ZTOT_RS_DEPS,                                           & ! deposition of vapor on snow (DEPS)
     ZTOT_TH_DEPI, ZTOT_RI_DEPI,                                           & ! deposition of vapor on ice (DEPI)
     ZTOT_RI_CNVS, ZTOT_CI_CNVS,                                           & ! conversion ice -> snow (CNVS)
     ZTOT_CS_SSC,                                                          & ! self collection of snow (SSC)
     ZTOT_CI_ISC,                                                          & ! self collection of ice (ISC)
     ZTOT_RI_AGGS, ZTOT_CI_AGGS,                                           & ! aggregation of ice on snow (AGGS)
     ZTOT_TH_DEPG, ZTOT_RG_DEPG,                                           & ! deposition of vapor on graupel (DEPG)
     ZTOT_TH_BERFI, ZTOT_RC_BERFI,                                         & ! Bergeron (BERFI)
     ZTOT_TH_RIM, ZTOT_CC_RIM, ZTOT_CS_RIM, ZTOT_RC_RIMSS, ZTOT_RC_RIMSG, ZTOT_RS_RIMCG, & ! cloud droplet riming (RIM)
     ZTOT_RI_HMS, ZTOT_CI_HMS, ZTOT_RS_HMS,                                & ! hallett mossop snow (HMS)
     ZTOT_TH_ACC, ZTOT_CR_ACC, ZTOT_CS_ACC, ZTOT_RR_ACCSS, ZTOT_RR_ACCSG, ZTOT_RS_ACCRG, & ! rain accretion on aggregates (ACC)
     ZTOT_RS_CMEL, ZTOT_CS_CMEL,                                                        & ! conversion-melting (CMEL)
     ZTOT_TH_CFRZ, ZTOT_RR_CFRZ, ZTOT_CR_CFRZ, ZTOT_RI_CFRZ, ZTOT_CI_CFRZ, & ! rain freezing (CFRZ)
     ZTOT_RI_CIBU, ZTOT_CI_CIBU,                                           & ! collisional ice break-up (CIBU)
     ZTOT_RI_RDSF, ZTOT_CI_RDSF,                                           & ! rain drops freezing shattering (RDSF)
     ZTOT_TH_WETG, ZTOT_RC_WETG, ZTOT_CC_WETG, ZTOT_RR_WETG, ZTOT_CR_WETG, & ! wet growth of graupel (WETG)
     ZTOT_RI_WETG, ZTOT_CI_WETG, ZTOT_RS_WETG, ZTOT_CS_WETG, ZTOT_RG_WETG, ZTOT_CG_WETG, ZTOT_RH_WETG, & ! wet growth of graupel (WETG)
     ZTOT_TH_DRYG, ZTOT_RC_DRYG, ZTOT_CC_DRYG, ZTOT_RR_DRYG, ZTOT_CR_DRYG, & ! dry growth of graupel (DRYG)
     ZTOT_RI_DRYG, ZTOT_CI_DRYG, ZTOT_RS_DRYG, ZTOT_CS_DRYG, ZTOT_RG_DRYG,               & ! dry growth of graupel (DRYG)
     ZTOT_RI_HMG, ZTOT_CI_HMG, ZTOT_RG_HMG,                                & ! hallett mossop graupel (HMG)
     ZTOT_TH_GMLT, ZTOT_RR_GMLT, ZTOT_CR_GMLT, ZTOT_CG_GMLT,               & ! graupel melting (GMLT)
     ZTOT_TH_DEPH, ZTOT_RH_DEPH,                                           & ! deposition of vapor on hail (DEPH)
     ZTOT_TH_WETH, ZTOT_RC_WETH, ZTOT_CC_WETH, ZTOT_RR_WETH, ZTOT_CR_WETH, & ! wet growth of hail (WETH)
     ZTOT_RI_WETH, ZTOT_CI_WETH, ZTOT_RS_WETH, ZTOT_CS_WETH, ZTOT_RG_WETH, ZTOT_CG_WETH, ZTOT_RH_WETH, & ! wet growth of hail (WETH)
     ZTOT_RG_COHG, ZTOT_CG_COHG,                                           & ! conversion of hail into graupel (COHG)
     ZTOT_TH_HMLT, ZTOT_RR_HMLT, ZTOT_CR_HMLT, ZTOT_CH_HMLT,               & ! hail melting (HMLT)
     ZTOT_RR_CVRC, ZTOT_CR_CVRC,                                           & ! conversion of rain into cloud droplets if diameter too small
     ZTOT_RV_CORR2, ZTOT_RC_CORR2, ZTOT_RR_CORR2, ZTOT_RI_CORR2,           &
     ZTOT_CC_CORR2, ZTOT_CR_CORR2, ZTOT_CI_CORR2,                          &
     ZTOT_RI_HIND, ZTOT_RC_HINC, ZTOT_RV_HENU, ZTOT_RV_HONH

REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTOT_IFNN_IMLT

REAL, DIMENSION(:,:),     ALLOCATABLE :: Z_SHCI_HONC, Z_SHCI_CNVI, &
                                         Z_SHCI_HMS,  Z_SHCI_HMG,  &
                                         Z_SHCI_CFRZ, Z_SHCI_CIBU, &
                                         Z_SHCI_RDSF, Z_SHCI_WETG, &
                                         Z_SHCI_DRYG, Z_SHCI_CORR2,&
                                         Z_SHCI_IMLT, Z_SHCI_HACH, &
                                         Z_SHCI_CNVS, Z_SHCI_AGGS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTOT_SHCI_HONC, ZTOT_SHCI_CNVI, &
                                         ZTOT_SHCI_HMS,  ZTOT_SHCI_HMG,  &
                                         ZTOT_SHCI_CFRZ, ZTOT_SHCI_CIBU, &
                                         ZTOT_SHCI_RDSF, ZTOT_SHCI_WETG, &
                                         ZTOT_SHCI_DRYG, ZTOT_SHCI_CORR2,&
                                         ZTOT_SHCI_IMLT, ZTOT_SHCI_HACH, &
                                         ZTOT_SHCI_CNVS, ZTOT_SHCI_AGGS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLBDAI_SHAPE
REAL, DIMENSION(:,:),     ALLOCATABLE :: Z_SHCI_ISC  ! ISC  : Ni COL/PLA/DRO + COL/PLA/DRO --> Ni BUR
REAL, DIMENSION(:,:),     ALLOCATABLE :: Z_SHCI_ISCS ! ISCS : Ni BUR + COL/PLA/DRO/BUR --> Ns
REAL, DIMENSION(:),       ALLOCATABLE :: Z_SHRI_ISCS ! ISCS : ri --> rs
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTOT_SHCI_ISC
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTOT_SHCI_ISCS ! perte pour chaque Ni, gain pour Ns = somme des pertes pour les Ni
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZTOT_SHRI_ISCS
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZAUX_ISCS  ! for snow budget: to take into account contribution from all shapes
!
!For mixing-ratio splitting
REAL, DIMENSION(D%NIJT,D%NKT) :: Z0RVT,   Z0RCT,   Z0RRT,   Z0RIT,   Z0RST,   Z0RGT,   Z0RHT
REAL, DIMENSION(:), ALLOCATABLE                      :: Z0RVT1D, Z0RCT1D, Z0RRT1D, Z0RIT1D, Z0RST1D, Z0RGT1D, Z0RHT1D 

!
! Loop control variables
REAL,    DIMENSION(D%NIJT,D%NKT)   :: ZTIME,   ZTIME_LASTCALL
INTEGER, DIMENSION(D%NIJT,D%NKT)   :: IITER
REAL,    DIMENSION(:), ALLOCATABLE :: ZTIME1D, ZTIME_LASTCALL1D, ZMAXTIME, ZTIME_THRESHOLD
INTEGER, DIMENSION(:), ALLOCATABLE :: IITER1D
LOGICAL, DIMENSION(D%NIJT,D%NKT)   :: GLCOMPUTE
LOGICAL, DIMENSION(:), ALLOCATABLE :: GLCOMPUTE1D
REAL                               :: ZTSTEP
INTEGER                            :: INB_ITER_MAX
!
!For subgrid clouds
REAL, DIMENSION(:), ALLOCATABLE                      :: ZCF1D, ZIF1D, ZPF1D     ! 1D packed cloud, ice and precip. frac.

!
! Various parameters
! domain size and levels (AROME compatibility)
! loops and packing
INTEGER :: II, IPACK, IN, IK
INTEGER :: IDX
INTEGER, DIMENSION(:), ALLOCATABLE :: I1, I3
! Inverse ov PTSTEP
REAL :: ZINV_TSTEP
! Work arrays
REAL, DIMENSION(D%NIJT)       :: ZW2D
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRT_SUM ! Total condensed water mr
REAL, DIMENSION(D%NIJT,D%NKT) :: ZCPT    ! Total condensed water mr
LOGICAL, DIMENSION(D%NIJT)          :: GDEP
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRHODJONTSTEP
!
INTEGER :: ISV_LIMA_NC
INTEGER :: ISV_LIMA_NR
INTEGER :: ISV_LIMA_CCN_FREE
INTEGER :: ISV_LIMA_CCN_ACTI
INTEGER :: ISV_LIMA_NI
INTEGER :: ISV_LIMA_NS
INTEGER :: ISV_LIMA_NG
INTEGER :: ISV_LIMA_NH
INTEGER :: ISV_LIMA_IFN_FREE
INTEGER :: ISV_LIMA_IFN_NUCL
INTEGER :: ISV_LIMA_IMM_NUCL
INTEGER :: ISV_LIMA_HOM_HAZE
INTEGER :: ISH !++cb--
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
! Variables for the electrification scheme
LOGICAL, DIMENSION(:,:), ALLOCATABLE :: GMASK_ELEC
INTEGER :: IELEC ! nb of points where the electrification scheme may apply
REAL, DIMENSION(:,:), ALLOCATABLE :: ZQPIT, ZQNIT, ZQCT, ZQRT, ZQIT, ZQST, ZQGT, ZQHT
REAL, DIMENSION(:,:), ALLOCATABLE :: ZQPIS, ZQNIS, ZQCS, ZQRS, ZQIS, ZQSS, ZQGS, ZQHS
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRVT_ELEC, ZRCT_ELEC, ZRRT_ELEC, ZRIT_ELEC, ZRST_ELEC, ZRGT_ELEC, ZRHT_ELEC
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCCT_ELEC, ZCRT_ELEC, ZCIT_ELEC, ZCST_ELEC, ZCGT_ELEC, ZCHT_ELEC
REAL, DIMENSION(:),     ALLOCATABLE :: ZLATHAM_IAGGS
!
! ICE4 COMPUTE PDF
LOGICAL, DIMENSION(D%NIJT,D%NKT) :: LLMICRO
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLC_LCF
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLC_LRC
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLI_LCF
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLI_LRI
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLC_HCF
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLC_HRC
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLI_HCF
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLI_HRI
REAL, DIMENSION(D%NIJT, D%NKT) :: ZSIGMA_RC
!-------------------------------------------------------------------------------
!
!*       0.     Init
!               ----
!
IF (LHOOK) CALL DR_HOOK('LIMA', 0, ZHOOK_HANDLE)
!
ISV_LIMA_NC       = TNSV%NSV_LIMA_NC       - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_NR       = TNSV%NSV_LIMA_NR       - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_CCN_FREE = TNSV%NSV_LIMA_CCN_FREE - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_CCN_ACTI = TNSV%NSV_LIMA_CCN_ACTI - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_NI       = TNSV%NSV_LIMA_NI       - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_NS       = TNSV%NSV_LIMA_NS       - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_NG       = TNSV%NSV_LIMA_NG       - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_NH       = TNSV%NSV_LIMA_NH       - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_IFN_FREE = TNSV%NSV_LIMA_IFN_FREE - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_IFN_NUCL = TNSV%NSV_LIMA_IFN_NUCL - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_IMM_NUCL = TNSV%NSV_LIMA_IMM_NUCL - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_HOM_HAZE = TNSV%NSV_LIMA_HOM_HAZE - TNSV%NSV_LIMA_BEG + 1
!
ZTHS(D%NIJB:D%NIJE,:) = PTHS(D%NIJB:D%NIJE,:)
ZTHT(D%NIJB:D%NIJE,:) = PTHS(D%NIJB:D%NIJE,:) * PTSTEP
ZRVT(:,:) = 0.
ZRVS(:,:) = 0.
ZRCT(:,:) = 0.
ZRCS(:,:) = 0.
ZRRT(:,:) = 0.
ZRRS(:,:) = 0.
ZRIT(:,:) = 0.
ZRIS(:,:) = 0.
ZRST(:,:) = 0.
ZRSS(:,:) = 0.
ZRGT(:,:) = 0.
ZRGS(:,:) = 0.
ZRHT(:,:) = 0.
ZRHS(:,:) = 0.
ZRT_SUM(:,:) = 0.
ZCCT(:,:)   = 0.
ZCCS(:,:)   = 0.
ZCRT(:,:)   = 0.
ZCRS(:,:)   = 0.
ZCIT(:,:)   = 0.
ZCIS(:,:)   = 0.
ZCST(:,:)   = 0.
ZCSS(:,:)   = 0.
ZCGT(:,:)   = 0.
ZCGS(:,:)   = 0.
ZCHT(:,:)   = 0.
ZCHS(:,:)   = 0.
ZCCNFT(:,:,:) = 0.
ZCCNAT(:,:,:) = 0.
ZCCNFS(:,:,:) = 0.
ZCCNAS(:,:,:) = 0.
ZIFNFT(:,:,:) = 0.
ZIFNNT(:,:,:) = 0.
ZIFNFS(:,:,:) = 0.
ZIFNNS(:,:,:) = 0.
ZIMMNT(:,:,:) = 0.
ZIMMNS(:,:,:) = 0.
ZHOMFT(:,:)   = 0.
ZHOMFS(:,:)   = 0.

IF ( BUCONF%LBU_ENABLE .OR. OELEC) THEN
  Z_RR_CVRC(:,:) = 0.
  Z_CR_CVRC(:,:) = 0.
  ALLOCATE( ZTOT_CR_BRKU (D%NIJT, D%NKT ) ); ZTOT_CR_BRKU(:,:) = 0.
  ALLOCATE( ZTOT_TH_HONR (D%NIJT, D%NKT ) ); ZTOT_TH_HONR(:,:) = 0.
  ALLOCATE( ZTOT_RR_HONR (D%NIJT, D%NKT ) ); ZTOT_RR_HONR(:,:) = 0.
  ALLOCATE( ZTOT_CR_HONR (D%NIJT, D%NKT ) ); ZTOT_CR_HONR(:,:) = 0.
  ALLOCATE( ZTOT_TH_IMLT (D%NIJT, D%NKT ) ); ZTOT_TH_IMLT(:,:) = 0.
  ALLOCATE( ZTOT_RC_IMLT (D%NIJT, D%NKT ) ); ZTOT_RC_IMLT(:,:) = 0.
  ALLOCATE( ZTOT_CC_IMLT (D%NIJT, D%NKT ) ); ZTOT_CC_IMLT(:,:) = 0.
  ALLOCATE( ZTOT_IFNN_IMLT (D%NIJT, D%NKT, LIMAP%NMOD_IFN ) ); ZTOT_IFNN_IMLT(:,:,:) = 0.
  ALLOCATE( ZTOT_TH_HONC (D%NIJT, D%NKT ) ); ZTOT_TH_HONC(:,:) = 0.
  ALLOCATE( ZTOT_RC_HONC (D%NIJT, D%NKT ) ); ZTOT_RC_HONC(:,:) = 0.
  ALLOCATE( ZTOT_CC_HONC (D%NIJT, D%NKT ) ); ZTOT_CC_HONC(:,:) = 0.
  ALLOCATE( ZTOT_CC_SELF (D%NIJT, D%NKT ) ); ZTOT_CC_SELF(:,:) = 0.
  ALLOCATE( ZTOT_RC_AUTO (D%NIJT, D%NKT ) ); ZTOT_RC_AUTO(:,:) = 0.
  ALLOCATE( ZTOT_CC_AUTO (D%NIJT, D%NKT ) ); ZTOT_CC_AUTO(:,:) = 0.
  ALLOCATE( ZTOT_CR_AUTO (D%NIJT, D%NKT ) ); ZTOT_CR_AUTO(:,:) = 0.
  ALLOCATE( ZTOT_RC_ACCR (D%NIJT, D%NKT ) ); ZTOT_RC_ACCR(:,:) = 0.
  ALLOCATE( ZTOT_CC_ACCR (D%NIJT, D%NKT ) ); ZTOT_CC_ACCR(:,:) = 0.
  ALLOCATE( ZTOT_CR_SCBU (D%NIJT, D%NKT ) ); ZTOT_CR_SCBU(:,:) = 0.
  ALLOCATE( ZTOT_TH_EVAP (D%NIJT, D%NKT ) ); ZTOT_TH_EVAP(:,:) = 0.
!   allocate( ZTOT_RC_EVAP (D%NIJT, D%NKT ) ); ZTOT_RC_EVAP(:,:) = 0.
!   allocate( ZTOT_CC_EVAP (D%NIJT, D%NKT ) ); ZTOT_CC_EVAP(:,:) = 0.
  ALLOCATE( ZTOT_RR_EVAP (D%NIJT, D%NKT ) ); ZTOT_RR_EVAP(:,:) = 0.
  ALLOCATE( ZTOT_CR_EVAP (D%NIJT, D%NKT ) ); ZTOT_CR_EVAP(:,:) = 0.
  ALLOCATE( ZTOT_RI_CNVI (D%NIJT, D%NKT ) ); ZTOT_RI_CNVI(:,:) = 0.
  ALLOCATE( ZTOT_CI_CNVI (D%NIJT, D%NKT ) ); ZTOT_CI_CNVI(:,:) = 0.
  ALLOCATE( ZTOT_TH_DEPS (D%NIJT, D%NKT ) ); ZTOT_TH_DEPS(:,:) = 0.
  ALLOCATE( ZTOT_RS_DEPS (D%NIJT, D%NKT ) ); ZTOT_RS_DEPS(:,:) = 0.
  ALLOCATE( ZTOT_TH_DEPI (D%NIJT, D%NKT ) ); ZTOT_TH_DEPI(:,:) = 0.
  ALLOCATE( ZTOT_RI_DEPI (D%NIJT, D%NKT ) ); ZTOT_RI_DEPI(:,:) = 0.
  ALLOCATE( ZTOT_RI_CNVS (D%NIJT, D%NKT ) ); ZTOT_RI_CNVS(:,:) = 0.
  ALLOCATE( ZTOT_CI_CNVS (D%NIJT, D%NKT ) ); ZTOT_CI_CNVS(:,:) = 0.
  ALLOCATE( ZTOT_CI_ISC  (D%NIJT, D%NKT ) ); ZTOT_CI_ISC(:,:)  = 0.
  ALLOCATE( ZTOT_CS_SSC  (D%NIJT, D%NKT ) ); ZTOT_CS_SSC(:,:) = 0.
  ALLOCATE( ZTOT_RI_AGGS (D%NIJT, D%NKT ) ); ZTOT_RI_AGGS(:,:) = 0.
  ALLOCATE( ZTOT_CI_AGGS (D%NIJT, D%NKT ) ); ZTOT_CI_AGGS(:,:) = 0.
  ALLOCATE( ZTOT_TH_DEPG (D%NIJT, D%NKT ) ); ZTOT_TH_DEPG(:,:) = 0.
  ALLOCATE( ZTOT_RG_DEPG (D%NIJT, D%NKT ) ); ZTOT_RG_DEPG(:,:) = 0.
  ALLOCATE( ZTOT_TH_BERFI(D%NIJT, D%NKT ) ); ZTOT_TH_BERFI(:,:) = 0.
  ALLOCATE( ZTOT_RC_BERFI(D%NIJT, D%NKT ) ); ZTOT_RC_BERFI(:,:) = 0.
!++cb++ need rcrimss, rcrimsg and rsrimcg to be consistent with ice3
  ALLOCATE( ZTOT_TH_RIM  (D%NIJT, D%NKT ) ); ZTOT_TH_RIM(:,:) = 0.
!  allocate( ZTOT_RC_RIM  (D%NIJT, D%NKT ) ); ZTOT_RC_RIM(:,:) = 0.
  ALLOCATE( ZTOT_CC_RIM  (D%NIJT, D%NKT ) ); ZTOT_CC_RIM(:,:) = 0.
!  allocate( ZTOT_RS_RIM  (D%NIJT, D%NKT ) ); ZTOT_RS_RIM(:,:) = 0.
  ALLOCATE( ZTOT_CS_RIM  (D%NIJT, D%NKT ) ); ZTOT_CS_RIM(:,:) = 0.
!  allocate( ZTOT_RG_RIM  (D%NIJT, D%NKT ) ); ZTOT_RG_RIM(:,:) = 0.
  ALLOCATE( ZTOT_RC_RIMSS (D%NIJT, D%NKT ) ); ZTOT_RC_RIMSS(:,:) = 0.
  ALLOCATE( ZTOT_RC_RIMSG (D%NIJT, D%NKT ) ); ZTOT_RC_RIMSG(:,:) = 0.
  ALLOCATE( ZTOT_RS_RIMCG (D%NIJT, D%NKT ) ); ZTOT_RS_RIMCG(:,:) = 0.
!--cb--
  ALLOCATE( ZTOT_RI_HMS  (D%NIJT, D%NKT ) ); ZTOT_RI_HMS(:,:) = 0.
  ALLOCATE( ZTOT_CI_HMS  (D%NIJT, D%NKT ) ); ZTOT_CI_HMS(:,:) = 0.
  ALLOCATE( ZTOT_RS_HMS  (D%NIJT, D%NKT ) ); ZTOT_RS_HMS(:,:) = 0.
!++cb++ need rraccss, rraccsg and rsaccrg to be consistent with ice3
  ALLOCATE( ZTOT_TH_ACC  (D%NIJT, D%NKT ) ); ZTOT_TH_ACC(:,:) = 0.
!  allocate( ZTOT_RR_ACC  (D%NIJT, D%NKT ) ); ZTOT_RR_ACC(:,:) = 0.
  ALLOCATE( ZTOT_CR_ACC  (D%NIJT, D%NKT ) ); ZTOT_CR_ACC(:,:) = 0.
!  allocate( ZTOT_RS_ACC  (D%NIJT, D%NKT ) ); ZTOT_RS_ACC(:,:) = 0.
  ALLOCATE( ZTOT_CS_ACC  (D%NIJT, D%NKT ) ); ZTOT_CS_ACC(:,:) = 0.
!  allocate( ZTOT_RG_ACC  (D%NIJT, D%NKT ) ); ZTOT_RG_ACC(:,:) = 0.
  ALLOCATE( ZTOT_RR_ACCSS(D%NIJT, D%NKT ) ); ZTOT_RR_ACCSS(:,:) = 0.
  ALLOCATE( ZTOT_RR_ACCSG(D%NIJT, D%NKT ) ); ZTOT_RR_ACCSG(:,:) = 0.
  ALLOCATE( ZTOT_RS_ACCRG(D%NIJT, D%NKT ) ); ZTOT_RS_ACCRG(:,:) = 0.
!--cb--
  ALLOCATE( ZTOT_RS_CMEL (D%NIJT, D%NKT ) ); ZTOT_RS_CMEL(:,:) = 0.
  ALLOCATE( ZTOT_CS_CMEL (D%NIJT, D%NKT ) ); ZTOT_CS_CMEL(:,:) = 0.
  ALLOCATE( ZTOT_TH_CFRZ (D%NIJT, D%NKT ) ); ZTOT_TH_CFRZ(:,:) = 0.
  ALLOCATE( ZTOT_RR_CFRZ (D%NIJT, D%NKT ) ); ZTOT_RR_CFRZ(:,:) = 0.
  ALLOCATE( ZTOT_CR_CFRZ (D%NIJT, D%NKT ) ); ZTOT_CR_CFRZ(:,:) = 0.
  ALLOCATE( ZTOT_RI_CFRZ (D%NIJT, D%NKT ) ); ZTOT_RI_CFRZ(:,:) = 0.
  ALLOCATE( ZTOT_CI_CFRZ (D%NIJT, D%NKT ) ); ZTOT_CI_CFRZ(:,:) = 0.
  ALLOCATE( ZTOT_RI_CIBU (D%NIJT, D%NKT ) ); ZTOT_RI_CIBU(:,:) = 0.
  ALLOCATE( ZTOT_CI_CIBU (D%NIJT, D%NKT ) ); ZTOT_CI_CIBU(:,:) = 0.
  ALLOCATE( ZTOT_RI_RDSF (D%NIJT, D%NKT ) ); ZTOT_RI_RDSF(:,:) = 0.
  ALLOCATE( ZTOT_CI_RDSF (D%NIJT, D%NKT ) ); ZTOT_CI_RDSF(:,:) = 0.
  ALLOCATE( ZTOT_TH_WETG (D%NIJT, D%NKT ) ); ZTOT_TH_WETG(:,:) = 0.
  ALLOCATE( ZTOT_RC_WETG (D%NIJT, D%NKT ) ); ZTOT_RC_WETG(:,:) = 0.
  ALLOCATE( ZTOT_CC_WETG (D%NIJT, D%NKT ) ); ZTOT_CC_WETG(:,:) = 0.
  ALLOCATE( ZTOT_RR_WETG (D%NIJT, D%NKT ) ); ZTOT_RR_WETG(:,:) = 0.
  ALLOCATE( ZTOT_CR_WETG (D%NIJT, D%NKT ) ); ZTOT_CR_WETG(:,:) = 0.
  ALLOCATE( ZTOT_RI_WETG (D%NIJT, D%NKT ) ); ZTOT_RI_WETG(:,:) = 0.
  ALLOCATE( ZTOT_CI_WETG (D%NIJT, D%NKT ) ); ZTOT_CI_WETG(:,:) = 0.
  ALLOCATE( ZTOT_RS_WETG (D%NIJT, D%NKT ) ); ZTOT_RS_WETG(:,:) = 0.
  ALLOCATE( ZTOT_CS_WETG (D%NIJT, D%NKT ) ); ZTOT_CS_WETG(:,:) = 0.
  ALLOCATE( ZTOT_RG_WETG (D%NIJT, D%NKT ) ); ZTOT_RG_WETG(:,:) = 0.
  ALLOCATE( ZTOT_CG_WETG (D%NIJT, D%NKT ) ); ZTOT_CG_WETG(:,:) = 0.
  ALLOCATE( ZTOT_RH_WETG (D%NIJT, D%NKT ) ); ZTOT_RH_WETG(:,:) = 0.
  ALLOCATE( ZTOT_TH_DRYG (D%NIJT, D%NKT ) ); ZTOT_TH_DRYG(:,:) = 0.
  ALLOCATE( ZTOT_RC_DRYG (D%NIJT, D%NKT ) ); ZTOT_RC_DRYG(:,:) = 0.
  ALLOCATE( ZTOT_CC_DRYG (D%NIJT, D%NKT ) ); ZTOT_CC_DRYG(:,:) = 0.
  ALLOCATE( ZTOT_RR_DRYG (D%NIJT, D%NKT ) ); ZTOT_RR_DRYG(:,:) = 0.
  ALLOCATE( ZTOT_CR_DRYG (D%NIJT, D%NKT ) ); ZTOT_CR_DRYG(:,:) = 0.
  ALLOCATE( ZTOT_RI_DRYG (D%NIJT, D%NKT ) ); ZTOT_RI_DRYG(:,:) = 0.
  ALLOCATE( ZTOT_CI_DRYG (D%NIJT, D%NKT ) ); ZTOT_CI_DRYG(:,:) = 0.
  ALLOCATE( ZTOT_RS_DRYG (D%NIJT, D%NKT ) ); ZTOT_RS_DRYG(:,:) = 0.
  ALLOCATE( ZTOT_CS_DRYG (D%NIJT, D%NKT ) ); ZTOT_CS_DRYG(:,:) = 0.
  ALLOCATE( ZTOT_RG_DRYG (D%NIJT, D%NKT ) ); ZTOT_RG_DRYG(:,:) = 0.
  ALLOCATE( ZTOT_RI_HMG  (D%NIJT, D%NKT ) ); ZTOT_RI_HMG(:,:) = 0.
  ALLOCATE( ZTOT_CI_HMG  (D%NIJT, D%NKT ) ); ZTOT_CI_HMG(:,:) = 0.
  ALLOCATE( ZTOT_RG_HMG  (D%NIJT, D%NKT ) ); ZTOT_RG_HMG(:,:) = 0.
  ALLOCATE( ZTOT_TH_GMLT (D%NIJT, D%NKT ) ); ZTOT_TH_GMLT(:,:) = 0.
  ALLOCATE( ZTOT_RR_GMLT (D%NIJT, D%NKT ) ); ZTOT_RR_GMLT(:,:) = 0.
  ALLOCATE( ZTOT_CR_GMLT (D%NIJT, D%NKT ) ); ZTOT_CR_GMLT(:,:) = 0.
  ALLOCATE( ZTOT_CG_GMLT (D%NIJT, D%NKT ) ); ZTOT_CG_GMLT(:,:) = 0.
  ALLOCATE( ZTOT_TH_DEPH (D%NIJT, D%NKT ) ); ZTOT_TH_DEPH(:,:) = 0.
  ALLOCATE( ZTOT_RH_DEPH (D%NIJT, D%NKT ) ); ZTOT_RH_DEPH(:,:) = 0.
  ALLOCATE( ZTOT_TH_WETH (D%NIJT, D%NKT ) ); ZTOT_TH_WETH(:,:) = 0.
  ALLOCATE( ZTOT_RC_WETH (D%NIJT, D%NKT ) ); ZTOT_RC_WETH(:,:) = 0.
  ALLOCATE( ZTOT_CC_WETH (D%NIJT, D%NKT ) ); ZTOT_CC_WETH(:,:) = 0.
  ALLOCATE( ZTOT_RR_WETH (D%NIJT, D%NKT ) ); ZTOT_RR_WETH(:,:) = 0.
  ALLOCATE( ZTOT_CR_WETH (D%NIJT, D%NKT ) ); ZTOT_CR_WETH(:,:) = 0.
  ALLOCATE( ZTOT_RI_WETH (D%NIJT, D%NKT ) ); ZTOT_RI_WETH(:,:) = 0.
  ALLOCATE( ZTOT_CI_WETH (D%NIJT, D%NKT ) ); ZTOT_CI_WETH(:,:) = 0.
  ALLOCATE( ZTOT_RS_WETH (D%NIJT, D%NKT ) ); ZTOT_RS_WETH(:,:) = 0.
  ALLOCATE( ZTOT_CS_WETH (D%NIJT, D%NKT ) ); ZTOT_CS_WETH(:,:) = 0.
  ALLOCATE( ZTOT_RG_WETH (D%NIJT, D%NKT ) ); ZTOT_RG_WETH(:,:) = 0.
  ALLOCATE( ZTOT_CG_WETH (D%NIJT, D%NKT ) ); ZTOT_CG_WETH(:,:) = 0.
  ALLOCATE( ZTOT_RH_WETH (D%NIJT, D%NKT ) ); ZTOT_RH_WETH(:,:) = 0.
  ALLOCATE( ZTOT_RG_COHG (D%NIJT, D%NKT ) ); ZTOT_RG_COHG(:,:) = 0.
  ALLOCATE( ZTOT_CG_COHG (D%NIJT, D%NKT ) ); ZTOT_CG_COHG(:,:) = 0.
  ALLOCATE( ZTOT_TH_HMLT (D%NIJT, D%NKT ) ); ZTOT_TH_HMLT(:,:) = 0.
  ALLOCATE( ZTOT_RR_HMLT (D%NIJT, D%NKT ) ); ZTOT_RR_HMLT(:,:) = 0.
  ALLOCATE( ZTOT_CR_HMLT (D%NIJT, D%NKT ) ); ZTOT_CR_HMLT(:,:) = 0.
  ALLOCATE( ZTOT_CH_HMLT (D%NIJT, D%NKT ) ); ZTOT_CH_HMLT(:,:) = 0.
  ALLOCATE( ZTOT_RR_CVRC (D%NIJT, D%NKT ) ); ZTOT_RR_CVRC(:,:) = 0.
  ALLOCATE( ZTOT_CR_CVRC (D%NIJT, D%NKT ) ); ZTOT_CR_CVRC(:,:) = 0.

  ALLOCATE( ZTOT_RV_CORR2 (D%NIJT, D%NKT ) ); ZTOT_RV_CORR2(:,:) = 0.
  ALLOCATE( ZTOT_RC_CORR2 (D%NIJT, D%NKT ) ); ZTOT_RC_CORR2(:,:) = 0.
  ALLOCATE( ZTOT_RR_CORR2 (D%NIJT, D%NKT ) ); ZTOT_RR_CORR2(:,:) = 0.
  ALLOCATE( ZTOT_RI_CORR2 (D%NIJT, D%NKT ) ); ZTOT_RI_CORR2(:,:) = 0.
  ALLOCATE( ZTOT_CC_CORR2 (D%NIJT, D%NKT ) ); ZTOT_CC_CORR2(:,:) = 0.
  ALLOCATE( ZTOT_CR_CORR2 (D%NIJT, D%NKT ) ); ZTOT_CR_CORR2(:,:) = 0.
  ALLOCATE( ZTOT_CI_CORR2 (D%NIJT, D%NKT ) ); ZTOT_CI_CORR2(:,:) = 0.
END IF
!
IF ( BUCONF%LBU_ENABLE .AND. LIMAP%LCRYSTAL_SHAPE ) THEN
  ALLOCATE( ZTOT_SHCI_IMLT(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_IMLT(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHCI_CORR2(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ); ZTOT_SHCI_CORR2(:,:,:) = 0.
  ALLOCATE( ZTOT_SHCI_HONC(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_HONC(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHCI_AGGS(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_AGGS(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHCI_CNVI(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_CNVI(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHCI_HMS(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) )  ; ZTOT_SHCI_HMS(:,:,:)   = 0.
  ALLOCATE( ZTOT_SHCI_CFRZ(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_CFRZ(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHCI_CIBU(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_CIBU(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHCI_RDSF(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_RDSF(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHCI_WETG(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_WETG(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHCI_DRYG(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_DRYG(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHCI_HMG(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) )  ; ZTOT_SHCI_HMG(:,:,:)   = 0.
  ALLOCATE( ZTOT_SHCI_CNVS(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_CNVS(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHCI_HACH(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_HACH(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHCI_ISC(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) )  ; ZTOT_SHCI_ISC(:,:,:)   = 0.
  ALLOCATE( ZTOT_SHCI_ISCS(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE) ) ; ZTOT_SHCI_ISCS(:,:,:)  = 0.
  ALLOCATE( ZTOT_SHRI_ISCS(D%NIJT,D%NKT) )                         ; ZTOT_SHRI_ISCS(:,:)    = 0.
END IF!
!
ALLOCATE (ZTOT_RI_HIND(D%NIJT,D%NKT)) ; ZTOT_RI_HIND(:,:) = 0.
ALLOCATE (ZTOT_RC_HINC(D%NIJT,D%NKT)) ; ZTOT_RC_HINC(:,:) = 0.
ALLOCATE (ZTOT_RV_HENU(D%NIJT,D%NKT)) ; ZTOT_RV_HENU(:,:) = 0.
ALLOCATE (ZTOT_RV_HONH(D%NIJT,D%NKT)) ; ZTOT_RV_HONH(:,:) = 0.
!
ZINV_TSTEP  = 1./PTSTEP
ZEXN(D%NIJB:D%NIJE,:) = (PPABST(D%NIJB:D%NIJE,:)/CST%XP00)**(CST%XRD/CST%XCPD)
ZT(D%NIJB:D%NIJE,:)   = ZTHT(D%NIJB:D%NIJE,:) * ZEXN(D%NIJB:D%NIJE,:)
!
! Initial values computed as source * PTSTEP
!
! Mixing ratios
!
ZRVT(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,1) * PTSTEP
ZRVS(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,1)
IF ( KRR .GE. 2 ) ZRCT(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,2) * PTSTEP
IF ( KRR .GE. 2 ) ZRCS(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,2)
IF ( KRR .GE. 3 ) ZRRT(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,3) * PTSTEP
IF ( KRR .GE. 3 ) ZRRS(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,3)
IF ( KRR .GE. 4 ) ZRIT(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,4) * PTSTEP
IF ( KRR .GE. 4 ) ZRIS(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,4)
IF ( KRR .GE. 5 ) ZRST(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,5) * PTSTEP
IF ( KRR .GE. 5 ) ZRSS(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,5)
IF ( KRR .GE. 6 ) ZRGT(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,6) * PTSTEP
IF ( KRR .GE. 6 ) ZRGS(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,6)
IF ( KRR .GE. 7 ) ZRHT(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,7) * PTSTEP
IF ( KRR .GE. 7 ) ZRHS(D%NIJB:D%NIJE,:) = PRS(D%NIJB:D%NIJE,:,7)
!
! Concentrations
!
IF ( LIMAP%NMOM_C.GE.2) THEN
   ZCCT(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NC) * PTSTEP
   ZCCS(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NC)
ELSE
   IF (LIMAP%LICE3) THEN
      ZCCT(D%NIJB:D%NIJE,:)   = 300.E6 / PRHODREF(D%NIJB:D%NIJE,:)   
   ELSE
      ZCCT(D%NIJB:D%NIJE,:)   = 100.E6 / PRHODREF(D%NIJB:D%NIJE,:)
   END IF
   ZCCS(D%NIJB:D%NIJE,:)   = ZCCT(D%NIJB:D%NIJE,:) / PTSTEP
END IF
IF ( LIMAP%NMOM_I.GE.2) THEN
  IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
    ZCIT(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NI) * PTSTEP
    ZCIS(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NI)
    ALLOCATE(ZCIT_SHAPE (0,0,0))
    ALLOCATE(ZCIS_SHAPE (0,0,0))
    ALLOCATE(ZRIT_SHAPE (0,0,0))
    ALLOCATE(ZRIS_SHAPE (0,0,0))
  ELSE
    ALLOCATE(ZCIT_SHAPE (D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE)); ZCIT_SHAPE(:,:,:) = 0.
    ALLOCATE(ZCIS_SHAPE (D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE)); ZCIS_SHAPE(:,:,:) = 0.
    ALLOCATE(ZRIT_SHAPE (D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE)); ZRIT_SHAPE(:,:,:) = 0.
    ALLOCATE(ZRIS_SHAPE (D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE)); ZRIS_SHAPE(:,:,:) = 0.
    DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
      ZCIT_SHAPE(D%NIJB:D%NIJE,:,ISH) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NI+ISH-1) * PTSTEP
      ZCIS_SHAPE(D%NIJB:D%NIJE,:,ISH) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NI+ISH-1)
    END DO
  END IF
ELSE
   ! ICE3 uses PCIT in m-3, but LIMA is in kg-1
   ZCIT(:,:)=PCIT(:,:)/PRHODREF(:,:)
   ZCIS(:,:)=PCIT(:,:)/PRHODREF(:,:)*ZINV_TSTEP
   ALLOCATE(ZCIT_SHAPE (0,0,0))
   ALLOCATE(ZCIS_SHAPE (0,0,0))
   ALLOCATE(ZRIT_SHAPE (0,0,0))
   ALLOCATE(ZRIS_SHAPE (0,0,0))
END IF
IF ( LIMAP%NMOM_R.GE.2) ZCRT(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NR) * PTSTEP
IF ( LIMAP%NMOM_R.GE.2) ZCRS(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NR)
IF ( LIMAP%NMOM_I.GE.2) ZCIT(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NI) * PTSTEP
IF ( LIMAP%NMOM_I.GE.2) ZCIS(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NI)
IF ( LIMAP%NMOM_S.GE.2) ZCST(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NS) * PTSTEP
IF ( LIMAP%NMOM_S.GE.2) ZCSS(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NS)
IF ( LIMAP%NMOM_G.GE.2) ZCGT(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NG) * PTSTEP
IF ( LIMAP%NMOM_G.GE.2) ZCGS(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NG)
IF ( LIMAP%NMOM_H.GE.2) ZCHT(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NH) * PTSTEP
IF ( LIMAP%NMOM_H.GE.2) ZCHS(D%NIJB:D%NIJE,:)   = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NH)
!
IF ( LIMAP%NMOD_CCN .GE. 1 ) ZCCNFT(D%NIJB:D%NIJE,:,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_CCN_FREE:ISV_LIMA_CCN_FREE+LIMAP%NMOD_CCN-1) * PTSTEP
IF ( LIMAP%NMOD_CCN .GE. 1 ) ZCCNAT(D%NIJB:D%NIJE,:,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_CCN_ACTI:ISV_LIMA_CCN_ACTI+LIMAP%NMOD_CCN-1) * PTSTEP
IF ( LIMAP%NMOD_CCN .GE. 1 ) ZCCNFS(D%NIJB:D%NIJE,:,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_CCN_FREE:ISV_LIMA_CCN_FREE+LIMAP%NMOD_CCN-1)
IF ( LIMAP%NMOD_CCN .GE. 1 ) ZCCNAS(D%NIJB:D%NIJE,:,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_CCN_ACTI:ISV_LIMA_CCN_ACTI+LIMAP%NMOD_CCN-1)
!
IF ( LIMAP%NMOD_IFN .GE. 1 ) ZIFNFT(D%NIJB:D%NIJE,:,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_IFN_FREE:ISV_LIMA_IFN_FREE+LIMAP%NMOD_IFN-1) * PTSTEP
IF ( LIMAP%NMOD_IFN .GE. 1 ) ZIFNNT(D%NIJB:D%NIJE,:,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_IFN_NUCL:ISV_LIMA_IFN_NUCL+LIMAP%NMOD_IFN-1) * PTSTEP
IF ( LIMAP%NMOD_IFN .GE. 1 ) ZIFNFS(D%NIJB:D%NIJE,:,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_IFN_FREE:ISV_LIMA_IFN_FREE+LIMAP%NMOD_IFN-1)
IF ( LIMAP%NMOD_IFN .GE. 1 ) ZIFNNS(D%NIJB:D%NIJE,:,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_IFN_NUCL:ISV_LIMA_IFN_NUCL+LIMAP%NMOD_IFN-1)
!
IF ( LIMAP%NMOD_IMM .GE. 1 ) ZIMMNT(D%NIJB:D%NIJE,:,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_IMM_NUCL:ISV_LIMA_IMM_NUCL+LIMAP%NMOD_IMM-1) * PTSTEP
IF ( LIMAP%NMOD_IMM .GE. 1 ) ZIMMNS(D%NIJB:D%NIJE,:,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_IMM_NUCL:ISV_LIMA_IMM_NUCL+LIMAP%NMOD_IMM-1)
!
IF ( LIMAP%LHHONI ) ZHOMFT(D%NIJB:D%NIJE,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_HOM_HAZE) * PTSTEP
IF ( LIMAP%LHHONI ) ZHOMFS(D%NIJB:D%NIJE,:) = PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_HOM_HAZE)
!
! Electric charge density
!
IF (OELEC) THEN
  ALLOCATE(ZQPIT(D%NIJT,D%NKT))
  ALLOCATE(ZQCT(D%NIJT,D%NKT))
  ALLOCATE(ZQRT(D%NIJT,D%NKT))
  ALLOCATE(ZQIT(D%NIJT,D%NKT))
  ALLOCATE(ZQST(D%NIJT,D%NKT))
  ALLOCATE(ZQGT(D%NIJT,D%NKT))
  ALLOCATE(ZQNIT(D%NIJT,D%NKT))
  IF (KRR == 7) ALLOCATE(ZQHT(D%NIJT,D%NKT))
  !
  ALLOCATE(ZQPIS(D%NIJT,D%NKT))
  ALLOCATE(ZQCS(D%NIJT,D%NKT))
  ALLOCATE(ZQRS(D%NIJT,D%NKT))
  ALLOCATE(ZQIS(D%NIJT,D%NKT))
  ALLOCATE(ZQSS(D%NIJT,D%NKT))
  ALLOCATE(ZQGS(D%NIJT,D%NKT))
  ALLOCATE(ZQNIS(D%NIJT,D%NKT))
  IF (KRR == 7) ALLOCATE(ZQHS(D%NIJT,D%NKT))
  !
  ALLOCATE(ZRVT_ELEC(D%NIJT,D%NKT))
  ALLOCATE(ZRCT_ELEC(D%NIJT,D%NKT))
  ALLOCATE(ZRRT_ELEC(D%NIJT,D%NKT))
  ALLOCATE(ZRIT_ELEC(D%NIJT,D%NKT))
  ALLOCATE(ZRST_ELEC(D%NIJT,D%NKT))
  ALLOCATE(ZRGT_ELEC(D%NIJT,D%NKT))
  IF (KRR == 7) ALLOCATE(ZRHT_ELEC(D%NIJT,D%NKT))
  ALLOCATE(ZCCT_ELEC(D%NIJT,D%NKT))
  ALLOCATE(ZCRT_ELEC(D%NIJT,D%NKT))
  ALLOCATE(ZCIT_ELEC(D%NIJT,D%NKT))
  ALLOCATE(ZCST_ELEC(D%NIJT,D%NKT))
  ALLOCATE(ZCGT_ELEC(D%NIJT,D%NKT))
  IF (KRR == 7) ALLOCATE(ZCHT_ELEC(D%NIJT,D%NKT))
!
!++cb++ 21/04/23 source * ptstep
  ZQPIT(D%NIJB:D%NIJE,:) = PSV_ELEC_S(D%NIJB:D%NIJE,:,1) * PTSTEP
  ZQCT(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,2) * PTSTEP
  ZQRT(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,3) * PTSTEP
  ZQIT(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,4) * PTSTEP
  ZQST(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,5) * PTSTEP
  ZQGT(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,6) * PTSTEP
  IF (KRR == 6) THEN
    ZQNIT(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,7) * PTSTEP
  ELSE IF (KRR == 7) THEN
    ZQHT(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,7) * PTSTEP
    ZQNIT(D%NIJB:D%NIJE,:) = PSV_ELEC_S(D%NIJB:D%NIJE,:,8) * PTSTEP
  END IF
  !
  ZQPIS(D%NIJB:D%NIJE,:) = PSV_ELEC_S(D%NIJB:D%NIJE,:,1)
  ZQCS(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,2)
  ZQRS(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,3)
  ZQIS(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,4)
  ZQSS(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,5)
  ZQGS(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,6)
  IF (KRR == 6) THEN
    ZQNIS(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,7)
  ELSE IF (KRR == 7) THEN
    ZQHS(D%NIJB:D%NIJE,:)  = PSV_ELEC_S(D%NIJB:D%NIJE,:,7)
    ZQNIS(D%NIJB:D%NIJE,:) = PSV_ELEC_S(D%NIJB:D%NIJE,:,8) 
  END IF
END IF
!
PINPRC=0.
PINDEP=0.
PINPRR=0.
PINPRI=0.
PINPRS=0.
PINPRG=0.
PINPRH=0.
PEVAP3D(:,:)=0.
!
!-------------------------------------------------------------------------------
!
!*       0.     Check mean diameter for cloud, rain  and ice
!               --------------------------------------------
! if ( BUCONF%LBU_ENABLE ) then
!   if ( BUCONF%lbudget_rc .and. lwarm .and. lrain ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'CORR', zrcs(:,:) * prhodj(:,:) )
!   if ( BUCONF%lbudget_rr .and. lwarm .and. lrain ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'CORR', zrrs(:,:) * prhodj(:,:) )
!   if ( BUCONF%lbudget_ri .and. lcold .and. lsnow ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'CORR', zris(:,:) * prhodj(:,:) )
!   if ( BUCONF%lbudget_rs .and. lcold .and. lsnow ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'CORR', zrss(:,:) * prhodj(:,:) )
!   if ( BUCONF%lbudget_sv ) then
!     if ( lwarm .and. lrain .and. nmom_c.ge.2) &
!       call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nc), 'CORR', zccs(:,:) * prhodj(:,:) )
!     if ( lwarm .and. lrain .and. nmom_r.ge.2) &
!       call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nr), 'CORR', zcrs(:,:) * prhodj(:,:) )
!     if ( lcold .and. lsnow .and. nmom_i.ge.2) &
!       call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CORR', zcis(:,:) * prhodj(:,:) )
!   end if
! end if
!!$IF (LIMAP%NMOM_R.GE.2) THEN
!!$   WHERE( ZRCT>LIMAP%XRTMIN(2) .AND. ZCCT>LIMAP%XCTMIN(2) .AND. ZRCT>XAC*ZCCT*(100.E-6)**LIMAW%XBC )
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
!!$IF (LIMAP%NMOM_R.GE.2) THEN
!!$   WHERE( ZRRT>LIMAP%XRTMIN(3) .AND. ZCRT>LIMAP%XCTMIN(3) .AND. ZRRT<LIMAW%XAR*ZCRT*(60.E-6)**LIMAW%XBR )
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
!!$IF (LIMAP%NMOM_S.GE.2) THEN
!!$   WHERE( ZRIT>LIMAP%XRTMIN(4) .AND. ZCIT>LIMAP%XCTMIN(4) .AND. ZRIT>LIMAC%XAI*ZCIT*(250.E-6)**LIMAC%XBI )
!!$      ZRST=ZRST+ZRIT
!!$      ZRSS=ZRSS+ZRIS
!!$      ZRIT=0.
!!$      ZCIT=0.
!!$      ZRIS=0.
!!$      ZCIS=0.
!!$   END WHERE
!!$END IF
!
! if ( BUCONF%LBU_ENABLE ) then
!   if ( BUCONF%lbudget_rc .and. lwarm .and. lrain ) call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'CORR', zrcs(:,:) * prhodj(:,:) )
!   if ( BUCONF%lbudget_rr .and. lwarm .and. lrain ) call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'CORR', zrrs(:,:) * prhodj(:,:) )
!   if ( BUCONF%lbudget_ri .and. lcold .and. lsnow ) call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'CORR', zris(:,:) * prhodj(:,:) )
!   if ( BUCONF%lbudget_rs .and. lcold .and. lsnow ) call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'CORR', zrss(:,:) * prhodj(:,:) )
!   if ( BUCONF%lbudget_sv ) then
!     if ( lwarm .and. lrain .and. nmom_c.ge.2) &
!       call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nc), 'CORR', zccs(:,:) * prhodj(:,:) )
!     if ( lwarm .and. lrain .and. nmom_r.ge.2) &
!       call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nr), 'CORR', zcrs(:,:) * prhodj(:,:) )
!     if ( lcold .and. lsnow .and. nmom_i.ge.2) &
!       call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CORR', zcis(:,:) * prhodj(:,:) )
!   end if
! end if
!-------------------------------------------------------------------------------
!
!*       1.     Sedimentation
!               -------------
!
!
IF ( BUCONF%LBU_ENABLE ) THEN
  IF ( BUCONF%LBUDGET_TH ) &
       CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'SEDI', ZTHS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RC .AND. LIMAP%NMOM_C.GE.1 .AND. LIMAP%LSEDC ) &
       CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'SEDI', ZRCS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RR .AND. LIMAP%NMOM_R.GE.1 ) &
       CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'SEDI', ZRRS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RI .AND. LIMAP%NMOM_I.GE.1 .AND. LIMAP%LSEDI ) &
       CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'SEDI', ZRIS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RS .AND. LIMAP%NMOM_S.GE.1 ) &
       CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'SEDI', ZRSS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RG .AND. LIMAP%NMOM_G.GE.1 ) &
       CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'SEDI', ZRGS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RH .AND. LIMAP%NMOM_H.GE.1 ) &
       CALL TBUDGETS(NBUDGET_RH)%PTR%INIT_PHY(D, 'SEDI', ZRHS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_SV ) THEN
    IF ( LIMAP%LSEDC .AND. LIMAP%NMOM_C.GE.2) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NC)%PTR%INIT_PHY(D, 'SEDI', ZCCS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%NMOM_R.GE.2) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NR)%PTR%INIT_PHY(D, 'SEDI', ZCRS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%LSEDI .AND. LIMAP%NMOM_I.GE.2) THEN
      IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NI)%PTR%INIT_PHY(D, 'SEDI', ZCIS(:,:) * PRHODJ(:,:) )
      ELSE
        DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
          CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NI + ISH - 1)%PTR%INIT_PHY(D, 'SEDI', ZCIS_SHAPE( :, :, ISH) * PRHODJ( :, :) )
        END DO
      END IF
    END IF
    IF ( LIMAP%NMOM_S.GE.2) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NS)%PTR%INIT_PHY(D, 'SEDI', ZCSS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%NMOM_G.GE.2) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NG)%PTR%INIT_PHY(D, 'SEDI', ZCGS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%NMOM_H.GE.2) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NH)%PTR%INIT_PHY(D, 'SEDI', ZCHS(:,:) * PRHODJ(:,:) )
    !
    IF (OELEC) THEN
      IF ( LIMAP%LSEDC ) &
        CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 1)%PTR%INIT_PHY(D, 'SEDI', ZQCS(:,:) * PRHODJ(:,:) )
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 2)%PTR%INIT_PHY(D, 'SEDI', ZQRS(:,:) * PRHODJ(:,:) )
      IF ( LIMAP%LSEDI ) &
        CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 3)%PTR%INIT_PHY(D, 'SEDI', ZQIS(:,:) * PRHODJ(:,:) )
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 4)%PTR%INIT_PHY(D, 'SEDI', ZQSS(:,:) * PRHODJ(:,:) )
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 5)%PTR%INIT_PHY(D, 'SEDI', ZQGS(:,:) * PRHODJ(:,:) )
      IF (LIMAP%NMOM_H .GE. 1) &
        CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 6)%PTR%INIT_PHY(D, 'SEDI', ZQHS(:,:) * PRHODJ(:,:) )
    END IF
  END IF
END IF
!
PFPR(:,:,:)=0.
!
! sedimentation of cloud droplets
ZRT_SUM(D%NIJB:D%NIJE,:) = (ZRVS(D%NIJB:D%NIJE,:) + ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:) + ZRIS(D%NIJB:D%NIJE,:) + &
                         &ZRSS(D%NIJB:D%NIJE,:) + ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:))*PTSTEP
ZCPT(D%NIJB:D%NIJE,:)    = CST%XCPD + (CST%XCPV * ZRVS(D%NIJB:D%NIJE,:) + &
                                     CST%XCL * (ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:)) + &
                                     CST%XCI * (ZRIS(D%NIJB:D%NIJE,:) + ZRSS(D%NIJB:D%NIJE,:) + &
                                                ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:)))*PTSTEP
IF (LIMAP%NMOM_C.GE.1 .AND. LIMAP%LSEDC) THEN
  IF (OELEC) THEN
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'L', 2, 2, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF, PTHVREFZIKB, PPABST, ZT, ZRT_SUM, ZCPT, &
                            ZRCS, ZCCS, PINPRC, PFPR(:,:,2), PEFIELDW, ZQCS)
  ELSE
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'L', 2, 2, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF, PTHVREFZIKB, PPABST, ZT, ZRT_SUM, ZCPT, &
                            ZRCS, ZCCS, PINPRC, PFPR(:,:,2))
  END IF
END IF
!
! sedimentation of raindrops
ZRT_SUM(D%NIJB:D%NIJE,:) = (ZRVS(D%NIJB:D%NIJE,:) + ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:) + ZRIS(D%NIJB:D%NIJE,:) + &
                          ZRSS(D%NIJB:D%NIJE,:) + ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:))*PTSTEP
ZCPT(D%NIJB:D%NIJE,:)    = CST%XCPD + (CST%XCPV * ZRVS(D%NIJB:D%NIJE,:) + &
                                     CST%XCL * (ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:)) + &
                                     CST%XCI * (ZRIS(D%NIJB:D%NIJE,:) + ZRSS(D%NIJB:D%NIJE,:) + &
                                                ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:)))*PTSTEP
IF (LIMAP%NMOM_R.GE.1) THEN
  IF (OELEC) THEN
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'L', LIMAP%NMOM_R, 3, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF, PTHVREFZIKB,PPABST, ZT, ZRT_SUM, ZCPT, &
                            ZRRS, ZCRS, PINPRR, PFPR(:,:,3), PEFIELDW, ZQRS)
  ELSE
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'L', LIMAP%NMOM_R, 3, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF,PTHVREFZIKB, PPABST, ZT, ZRT_SUM, ZCPT, &
                            ZRRS, ZCRS, PINPRR, PFPR(:,:,3))
  END IF
END IF
!
! sedimentation of ice crystals
ZRT_SUM(D%NIJB:D%NIJE,:) = (ZRVS(D%NIJB:D%NIJE,:) + ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:) + ZRIS(D%NIJB:D%NIJE,:) + &
                          ZRSS(D%NIJB:D%NIJE,:) + ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:))*PTSTEP
ZCPT(D%NIJB:D%NIJE,:)    = CST%XCPD + (CST%XCPV * ZRVS(D%NIJB:D%NIJE,:) + &
                                     CST%XCL * (ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:)) + &
                                     CST%XCI * (ZRIS(D%NIJB:D%NIJE,:) + ZRSS(D%NIJB:D%NIJE,:) + &
                                                ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:)))*PTSTEP
IF (LIMAP%NMOM_I.GE.1 .AND. LIMAP%LSEDI) THEN
  IF (OELEC) THEN
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'I', LIMAP%NMOM_I, 4, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF,PTHVREFZIKB, PPABST, ZT, ZRT_SUM, ZCPT, &
                            ZRIS, ZCIS, ZW2D, PFPR(:,:,4), PEFIELDW, ZQIS)
  ELSE
    IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'I', LIMAP%NMOM_I, 4, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF,PTHVREFZIKB, PPABST, ZT, ZRT_SUM, ZCPT, &
                            ZRIS, ZCIS, ZW2D, PFPR(:,:,4))
    ELSE
      ALLOCATE(ZLBDAI_SHAPE(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE))
      CALL LIMA_SHAPE_COMPUTE_LBDA_3D(D, LIMAP, LIMAC, ZRIS, ZCIS_SHAPE, ZRIS_SHAPE, ZLBDAI_SHAPE)
      DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
        CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
             'I', LIMAP%NMOM_I, 4, ISH, 1, PTSTEP, OELEC, ELECP, ELECD, &
             PDZZ, PRHODREF,PTHVREFZIKB, PPABST, ZT, ZRT_SUM, ZCPT, &
             ZRIS_SHAPE(:,:,ISH), ZCIS_SHAPE(:,:,ISH), ZW2D, PFPR(:,:,4), &
             PLBDAI_SHAPE=ZLBDAI_SHAPE(:,:,ISH))
      END DO
      ZCIS(D%NIJB:D%NIJE,:) = SUM(ZCIS_SHAPE(D%NIJB:D%NIJE,:,:),DIM=3)
      ZRIS(D%NIJB:D%NIJE,:) = SUM(ZRIS_SHAPE(D%NIJB:D%NIJE,:,:),DIM=3)
      DEALLOCATE(ZLBDAI_SHAPE)
    END IF
  END IF
END IF
!
! sedimentation of snow/aggregates
ZRT_SUM(D%NIJB:D%NIJE,:) = (ZRVS(D%NIJB:D%NIJE,:) + ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:) + ZRIS(D%NIJB:D%NIJE,:) + &
                          ZRSS(D%NIJB:D%NIJE,:) + ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:))*PTSTEP
ZCPT(D%NIJB:D%NIJE,:)    = CST%XCPD + (CST%XCPV * ZRVS(D%NIJB:D%NIJE,:) + &
                                     CST%XCL * (ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:)) + &
                                     CST%XCI * (ZRIS(D%NIJB:D%NIJE,:) + ZRSS(D%NIJB:D%NIJE,:) + &
                                                ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:)))*PTSTEP
IF (LIMAP%NMOM_S.GE.1) THEN
  IF (OELEC) THEN
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'I', LIMAP%NMOM_S, 5, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF, PTHVREFZIKB, PPABST, ZT, ZRT_SUM, ZCPT, &
                            ZRSS, ZCSS, PINPRS, PFPR(:,:,5), PEFIELDW, ZQSS)
  ELSE
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'I', LIMAP%NMOM_S, 5, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF, PTHVREFZIKB, PPABST,ZT, ZRT_SUM, ZCPT, &
                            ZRSS, ZCSS, PINPRS, PFPR(:,:,5))
  END IF
END IF
!
! sedimentation of graupel
ZRT_SUM(D%NIJB:D%NIJE,:) = (ZRVS(D%NIJB:D%NIJE,:) + ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:) + ZRIS(D%NIJB:D%NIJE,:) + &
                          ZRSS(D%NIJB:D%NIJE,:) + ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:))*PTSTEP
ZCPT(D%NIJB:D%NIJE,:)    = CST%XCPD + (CST%XCPV * ZRVS(D%NIJB:D%NIJE,:) + &
                                     CST%XCL * (ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:)) + &
                                     CST%XCI * (ZRIS(D%NIJB:D%NIJE,:) + ZRSS(D%NIJB:D%NIJE,:) + &
                                                ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:)))*PTSTEP
IF (LIMAP%NMOM_G.GE.1) THEN
  IF (OELEC) THEN
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'I', LIMAP%NMOM_G, 6, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF, PTHVREFZIKB, PPABST, ZT, ZRT_SUM, ZCPT, &
                            ZRGS, ZCGS, PINPRG, PFPR(:,:,6), PEFIELDW, ZQGS)
  ELSE
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'I', LIMAP%NMOM_G, 6, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF, PTHVREFZIKB, PPABST, ZT, ZRT_SUM, ZCPT, &
                            ZRGS, ZCGS, PINPRG, PFPR(:,:,6))
  END IF
END IF
!
! sedimentation of hail
ZRT_SUM(D%NIJB:D%NIJE,:) = (ZRVS(D%NIJB:D%NIJE,:) + ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:) + ZRIS(D%NIJB:D%NIJE,:) + &
                          ZRSS(D%NIJB:D%NIJE,:) + ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:))*PTSTEP
ZCPT(D%NIJB:D%NIJE,:)    = CST%XCPD + (CST%XCPV * ZRVS(D%NIJB:D%NIJE,:) + &
                                     CST%XCL * (ZRCS(D%NIJB:D%NIJE,:) + ZRRS(D%NIJB:D%NIJE,:)) + &
                                     CST%XCI * (ZRIS(D%NIJB:D%NIJE,:) + ZRSS(D%NIJB:D%NIJE,:) + &
                                                ZRGS(D%NIJB:D%NIJE,:) + ZRHS(D%NIJB:D%NIJE,:)))*PTSTEP
IF (LIMAP%NMOM_H.GE.1) THEN
  IF (OELEC) THEN
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'I', LIMAP%NMOM_H, 7, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF, PTHVREFZIKB, PPABST, ZT, ZRT_SUM, ZCPT, &
                            ZRHS, ZCHS, PINPRH, PFPR(:,:,7), PEFIELDW, ZQHS)
  ELSE
    CALL LIMA_SEDIMENTATION(LIMAP, LIMAW, LIMAC, LIMAM, D, CST, &
         'I', LIMAP%NMOM_H, 7, 0, 1, PTSTEP, OELEC, ELECP, ELECD, &
         PDZZ, PRHODREF, PTHVREFZIKB, PPABST, ZT, ZRT_SUM, ZCPT, &
                            ZRHS, ZCHS, PINPRH, PFPR(:,:,7))
  END IF
END IF
!
ZTHS(D%NIJB:D%NIJE,:) = ZT(D%NIJB:D%NIJE,:) / ZEXN(D%NIJB:D%NIJE,:) * ZINV_TSTEP
!
! Call budgets
!
IF ( BUCONF%LBU_ENABLE ) then
  IF ( BUCONF%LBUDGET_TH ) &
       CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'SEDI', ZTHS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RC .AND. LIMAP%NMOM_C.GE.1 .AND. LIMAP%LSEDC ) &
       CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'SEDI', ZRCS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RR .AND. LIMAP%NMOM_R.GE.1 ) &
       CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'SEDI', ZRRS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RI .AND. LIMAP%NMOM_I.GE.1 .AND. LIMAP%LSEDI ) &
       CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'SEDI', ZRIS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RS .AND. LIMAP%NMOM_S.GE.1 ) &
       CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'SEDI', ZRSS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RG .AND. LIMAP%NMOM_G.GE.1 ) &
       CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'SEDI', ZRGS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RH .AND. LIMAP%NMOM_H.GE.1 ) &
       CALL TBUDGETS(NBUDGET_RH)%PTR%END_PHY(D, 'SEDI', ZRHS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_SV ) then
    IF ( LIMAP%LSEDC .AND. LIMAP%NMOM_C.GE.2 ) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NC)%PTR%END_PHY(D, 'SEDI', ZCCS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%NMOM_R.GE.2 ) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NR)%PTR%END_PHY(D, 'SEDI', ZCRS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%LSEDI .AND. LIMAP%NMOM_I.GE.2 ) THEN
      IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NI)%PTR%END_PHY(D, 'SEDI', ZCIS(:,:) * PRHODJ(:,:) )
      ELSE
        DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
          CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NI + ISH - 1)%PTR%END_PHY(D, 'SEDI', ZCIS_SHAPE( :, :, ISH) * PRHODJ( :, :) )
        END DO
      END IF
    END IF
    IF ( LIMAP%NMOM_S.GE.2 ) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NS)%PTR%END_PHY(D, 'SEDI', ZCSS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%NMOM_G.GE.2 ) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NG)%PTR%END_PHY(D, 'SEDI', ZCGS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%NMOM_H.GE.2 ) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NH)%PTR%END_PHY(D, 'SEDI', ZCHS(:,:) * PRHODJ(:,:) )
!
    IF (OELEC) then
      IF ( LIMAP%LSEDC ) &
        CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 1)%PTR%END_PHY(D, 'SEDI', ZQCS(:,:) * PRHODJ(:,:) )
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 2)%PTR%END_PHY(D, 'SEDI', ZQRS(:,:) * PRHODJ(:,:) )
      IF ( LIMAP%LSEDI ) &
        CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 3)%PTR%END_PHY(D, 'SEDI', ZQIS(:,:) * PRHODJ(:,:) )
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 4)%PTR%END_PHY(D, 'SEDI', ZQSS(:,:) * PRHODJ(:,:) )
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 5)%PTR%END_PHY(D, 'SEDI', ZQGS(:,:) * PRHODJ(:,:) )
      IF (LIMAP%NMOM_H .GE. 1) &
        CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_ELECBEG + 6)%PTR%END_PHY(D, 'SEDI', ZQHS(:,:) * PRHODJ(:,:) )
    END IF
  END IF
END IF
!
! 1.bis Deposition at 1st level above ground
!
IF (LIMAP%NMOM_C.GE.1 .AND. LIMAP%LDEPOC) THEN
  IF ( BUCONF%LBUDGET_RC ) CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'DEPO', ZRCS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_SV .AND. LIMAP%NMOM_C.GE.2) &
       CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NC)%PTR%INIT_PHY(D, 'DEPO', ZCCS(:,:) * PRHODJ(:,:) )

  GDEP(:) = .FALSE.
  GDEP(D%NIJB:D%NIJE) = ZRCS(D%NIJB:D%NIJE,D%NKB) >0 .AND. ZCCS(D%NIJB:D%NIJE,D%NKB) >0 .AND. &
                      ZRCT(D%NIJB:D%NIJE,D%NKB) >0 .AND. ZCCT(D%NIJB:D%NIJE,D%NKB) >0
  WHERE (GDEP)
     ZRCS(D%NIJB:D%NIJE,D%NKB) = ZRCS(D%NIJB:D%NIJE,D%NKB) - LIMAP%XVDEPOC * ZRCT(D%NIJB:D%NIJE,D%NKB) / PDZZ(D%NIJB:D%NIJE,D%NKB)
     ZCCS(D%NIJB:D%NIJE,D%NKB) = ZCCS(D%NIJB:D%NIJE,D%NKB) - LIMAP%XVDEPOC * ZCCT(D%NIJB:D%NIJE,D%NKB) / PDZZ(D%NIJB:D%NIJE,D%NKB)
     PINPRC(D%NIJB:D%NIJE) = PINPRC(D%NIJB:D%NIJE) + LIMAP%XVDEPOC * ZRCT(D%NIJB:D%NIJE,D%NKB) * PRHODREF(D%NIJB:D%NIJE,D%NKB) /CST%XRHOLW
     PINDEP(D%NIJB:D%NIJE) = LIMAP%XVDEPOC * ZRCT(D%NIJB:D%NIJE,D%NKB) *  PRHODREF(D%NIJB:D%NIJE,D%NKB) /CST%XRHOLW
  END WHERE

  IF ( BUCONF%LBUDGET_RC ) CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'DEPO', ZRCS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_SV .AND. LIMAP%NMOM_C.GE.2) &
       CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NC)%PTR%END_PHY(D, 'DEPO', ZCCS(:,:) * PRHODJ(:,:) )
END IF
!
!
!!$Z_RR_CVRC(:,:) = 0.
!!$Z_CR_CVRC(:,:) = 0.
!!$IF (LIMAP%NMOM_R.GE.2) THEN
!!$   if( BUCONF%LBU_ENABLE ) then
!!$    if ( BUCONF%lbudget_rc ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC),                    'R2C1', zrcs(:,:) * prhodj(:,:) )
!!$    if ( BUCONF%lbudget_rr ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR),                    'R2C1', zrrs(:,:) * prhodj(:,:) )
!!$    if ( BUCONF%lbudget_sv .and. nmom_c.ge.2) &
!!$         call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nc), 'R2C1', zccs(:,:) * prhodj(:,:) )
!!$    if ( BUCONF%lbudget_sv .and. nmom_r.ge.2) &
!!$         call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nr), 'R2C1', zcrs(:,:) * prhodj(:,:) )
!!$   end if
!!$
!!$   CALL LIMA_DROPS_TO_DROPLETS_CONV(LIMAP, D, PRHODREF, ZRCS*PTSTEP, ZRRS*PTSTEP, ZCCS*PTSTEP, ZCRS*PTSTEP, &
!!$                                    Z_RR_CVRC, Z_CR_CVRC)
!!$   !
!!$   ZRCS(:,:) = ZRCS(:,:) - Z_RR_CVRC(:,:)/PTSTEP
!!$   ZRRS(:,:) = ZRRS(:,:) + Z_RR_CVRC(:,:)/PTSTEP
!!$   ZCCS(:,:) = ZCCS(:,:) - Z_CR_CVRC(:,:)/PTSTEP
!!$   ZCRS(:,:) = ZCRS(:,:) + Z_CR_CVRC(:,:)/PTSTEP
!!$
!!$   if( BUCONF%LBU_ENABLE ) then
!!$    if ( BUCONF%lbudget_rc ) call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC),                    'R2C1', zrcs(:,:) * prhodj(:,:) )
!!$    if ( BUCONF%lbudget_rr ) call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR),                    'R2C1', zrrs(:,:) * prhodj(:,:) )
!!$    if ( BUCONF%lbudget_sv .and. nmom_c.ge.2) &
!!$         call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nc), 'R2C1', zccs(:,:) * prhodj(:,:) )
!!$    if ( BUCONF%lbudget_sv .and. nmom_r.ge.2) &
!!$         call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nr), 'R2C1', zcrs(:,:) * prhodj(:,:) )
!!$   end if
!!$END IF
!
! Update variables
!
ZTHT(D%NIJB:D%NIJE,:) = ZTHS(D%NIJB:D%NIJE,:) * PTSTEP
ZT(D%NIJB:D%NIJE,:)   = ZTHT(D%NIJB:D%NIJE,:) * ZEXN(D%NIJB:D%NIJE,:)
!
IF ( KRR .GE. 2 ) ZRCT(D%NIJB:D%NIJE,:) = ZRCS(D%NIJB:D%NIJE,:) * PTSTEP
IF ( KRR .GE. 3 ) ZRRT(D%NIJB:D%NIJE,:) = ZRRS(D%NIJB:D%NIJE,:) * PTSTEP
IF ( KRR .GE. 4 ) ZRIT(D%NIJB:D%NIJE,:) = ZRIS(D%NIJB:D%NIJE,:) * PTSTEP
IF ( KRR .GE. 5 ) ZRST(D%NIJB:D%NIJE,:) = ZRSS(D%NIJB:D%NIJE,:) * PTSTEP
IF ( KRR .GE. 6 ) ZRGT(D%NIJB:D%NIJE,:) = ZRGS(D%NIJB:D%NIJE,:) * PTSTEP
IF ( KRR .GE. 7 ) ZRHT(D%NIJB:D%NIJE,:) = ZRHS(D%NIJB:D%NIJE,:) * PTSTEP
!
IF ( LIMAP%NMOM_C.GE.2 ) ZCCT(D%NIJB:D%NIJE,:) = ZCCS(D%NIJB:D%NIJE,:) * PTSTEP
IF ( LIMAP%NMOM_R.GE.2 ) ZCRT(D%NIJB:D%NIJE,:) = ZCRS(D%NIJB:D%NIJE,:) * PTSTEP
IF ( LIMAP%NMOM_I.GE.2 ) ZCIT(D%NIJB:D%NIJE,:) = ZCIS(D%NIJB:D%NIJE,:) * PTSTEP
IF ( LIMAP%NMOM_S.GE.2 ) ZCST(D%NIJB:D%NIJE,:) = ZCSS(D%NIJB:D%NIJE,:) * PTSTEP
IF ( LIMAP%NMOM_G.GE.2 ) ZCGT(D%NIJB:D%NIJE,:) = ZCGS(D%NIJB:D%NIJE,:) * PTSTEP
IF ( LIMAP%NMOM_H.GE.2 ) ZCHT(D%NIJB:D%NIJE,:) = ZCHS(D%NIJB:D%NIJE,:) * PTSTEP
!
IF (LIMAP%LCRYSTAL_SHAPE) THEN
  ZCIT_SHAPE(D%NIJB:D%NIJE,:,:) = ZCIS_SHAPE(D%NIJB:D%NIJE,:,:) * PTSTEP
  ZCIT(D%NIJB:D%NIJE,:) = SUM(ZCIT_SHAPE(D%NIJB:D%NIJE,:,:),DIM=3)
END IF
!
IF (OELEC) THEN
  ZQCT(D%NIJB:D%NIJE,:) = ZQCS(D%NIJB:D%NIJE,:) * PTSTEP
  ZQRT(D%NIJB:D%NIJE,:) = ZQRS(D%NIJB:D%NIJE,:) * PTSTEP
  ZQIT(D%NIJB:D%NIJE,:) = ZQIS(D%NIJB:D%NIJE,:) * PTSTEP
  ZQST(D%NIJB:D%NIJE,:) = ZQSS(D%NIJB:D%NIJE,:) * PTSTEP
  ZQGT(D%NIJB:D%NIJE,:) = ZQGS(D%NIJB:D%NIJE,:) * PTSTEP
  IF (LIMAP%NMOM_H .GE. 1) ZQHT(D%NIJB:D%NIJE,:) = ZQHS(D%NIJB:D%NIJE,:) * PTSTEP
END IF
! 
!-------------------------------------------------------------------------------
!
!*       2.     Compute cloud, ice and precipitation fractions
!               ----------------------------------------------
!
CALL LIMA_COMPUTE_CLOUD_FRACTIONS (LIMAP, D,                          &
                                   ZCCT, ZRCT,                        &
                                   ZCRT, ZRRT,                        &
                                   ZCIT, ZRIT,                        &
                                   ZCST, ZRST,                        &
                                   ZCGT, ZRGT,                        &
                                   ZCHT, ZRHT,                        &
                                   PCLDFR, PICEFR, PPRCFR             )
!
!
!-------------------------------------------------------------------------------
!
!*       2.     Nucleation processes
!               --------------------
!
CALL LIMA_NUCLEATION_PROCS (LIMAP, LIMAW, LIMAC, TNSV, D, CST, NEBN, BUCONF, TBUDGETS, KBUDGETS,  &
                            KCARB, KSOA, KSP, ODUST, OSALT, OORILAM,            &
                            PTSTEP, PRHODJ,                                     &
                            PRHODREF, ZEXN, PPABST, ZT, PDTHRAD, PW_NU,         &
                            ZTHT, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT, ZRHT,     &
                            ZCCT, ZCRT, ZCIT, ZCIT_SHAPE, PAERO,PSOLORG, PMI, HACTCCN, &
                            ZCCNFT, ZCCNAT, ZIFNFT, ZIFNNT, ZIMMNT, ZHOMFT,     &
                            PCLDFR, PICEFR, PPRCFR,                             &
                            ZTOT_RV_HENU, ZTOT_RC_HINC, ZTOT_RI_HIND, ZTOT_RV_HONH)
! dans lima_nucleation_proc, pcit a ete mis a jour comme la somme des pcit_shape si lcrystal_shape=t 
!
! Saving sources before microphysics time-splitting loop
!
ZRVS(D%NIJB:D%NIJE,:) = ZRVT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
ZRCS(D%NIJB:D%NIJE,:) = ZRCT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
ZRRS(D%NIJB:D%NIJE,:) = ZRRT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
ZRIS(D%NIJB:D%NIJE,:) = ZRIT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
ZRSS(D%NIJB:D%NIJE,:) = ZRST(D%NIJB:D%NIJE,:) *ZINV_TSTEP
ZRGS(D%NIJB:D%NIJE,:) = ZRGT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
ZRHS(D%NIJB:D%NIJE,:) = ZRHT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
!
IF (LIMAP%NMOM_C.GE.2) ZCCS(D%NIJB:D%NIJE,:) = ZCCT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF (LIMAP%NMOM_R.GE.2) ZCRS(D%NIJB:D%NIJE,:) = ZCRT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF (LIMAP%NMOM_I.GE.2) ZCIS(D%NIJB:D%NIJE,:) = ZCIT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF (LIMAP%NMOM_S.GE.2) ZCSS(D%NIJB:D%NIJE,:) = ZCST(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF (LIMAP%NMOM_G.GE.2) ZCGS(D%NIJB:D%NIJE,:) = ZCGT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF (LIMAP%NMOM_H.GE.2) ZCHS(D%NIJB:D%NIJE,:) = ZCHT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
!
IF (LIMAP%LCRYSTAL_SHAPE) THEN
  ZCIS_SHAPE(D%NIJB:D%NIJE,:,:) = ZCIT_SHAPE(D%NIJB:D%NIJE,:,:) * ZINV_TSTEP
  ZCIS(D%NIJB:D%NIJE,:) = SUM(ZCIS_SHAPE(D%NIJB:D%NIJE,:,:),DIM=3)
END IF
!
ZCCNFS(D%NIJB:D%NIJE,:,:) = ZCCNFT(D%NIJB:D%NIJE,:,:) *ZINV_TSTEP
ZCCNAS(D%NIJB:D%NIJE,:,:) = ZCCNAT(D%NIJB:D%NIJE,:,:) *ZINV_TSTEP
ZIFNFS(D%NIJB:D%NIJE,:,:) = ZIFNFT(D%NIJB:D%NIJE,:,:) *ZINV_TSTEP
ZIFNNS(D%NIJB:D%NIJE,:,:) = ZIFNNT(D%NIJB:D%NIJE,:,:) *ZINV_TSTEP
ZIMMNS(D%NIJB:D%NIJE,:,:) = ZIMMNT(D%NIJB:D%NIJE,:,:) *ZINV_TSTEP
ZHOMFS(D%NIJB:D%NIJE,:)   = ZHOMFT(D%NIJB:D%NIJE,:)   *ZINV_TSTEP
!
ZTHS(D%NIJB:D%NIJE,:) = ZTHT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
ZT(D%NIJB:D%NIJE,:)   = ZTHT(D%NIJB:D%NIJE,:) * ZEXN(D%NIJB:D%NIJE,:)
!
IF (OELEC) THEN
  ZRVT_ELEC(D%NIJB:D%NIJE,:) = ZRVT(D%NIJB:D%NIJE,:)
  ZRCT_ELEC(D%NIJB:D%NIJE,:) = ZRCT(D%NIJB:D%NIJE,:)
  ZRRT_ELEC(D%NIJB:D%NIJE,:) = ZRRT(D%NIJB:D%NIJE,:)
  ZRIT_ELEC(D%NIJB:D%NIJE,:) = ZRIT(D%NIJB:D%NIJE,:)
  ZRST_ELEC(D%NIJB:D%NIJE,:) = ZRST(D%NIJB:D%NIJE,:)
  ZRGT_ELEC(D%NIJB:D%NIJE,:) = ZRGT(D%NIJB:D%NIJE,:)
  IF (LIMAP%NMOM_H .GE. 1) ZRHT_ELEC(D%NIJB:D%NIJE,:) = ZRHT(D%NIJB:D%NIJE,:)
  IF (LIMAP%NMOM_C .GE. 2) ZCCT_ELEC(D%NIJB:D%NIJE,:) = ZCCT(D%NIJB:D%NIJE,:)
  IF (LIMAP%NMOM_R .GE. 2) ZCRT_ELEC(D%NIJB:D%NIJE,:) = ZCRT(D%NIJB:D%NIJE,:)
  IF (LIMAP%NMOM_I .GE. 2) ZCIT_ELEC(D%NIJB:D%NIJE,:) = ZCIT(D%NIJB:D%NIJE,:)
  IF (LIMAP%NMOM_S .GE. 2) ZCST_ELEC(D%NIJB:D%NIJE,:) = ZCST(D%NIJB:D%NIJE,:)
  IF (LIMAP%NMOM_G .GE. 2) ZCGT_ELEC(D%NIJB:D%NIJE,:) = ZCGT(D%NIJB:D%NIJE,:)
  IF (LIMAP%NMOM_H .GE. 2) ZCHT_ELEC(D%NIJB:D%NIJE,:) = ZCHT(D%NIJB:D%NIJE,:)
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE PRECIPITATION FRACTION
!               ------------------------------
!
LLMICRO(:,:)=.TRUE.
IF (KRR==7) THEN
   LLMICRO(:,:)=ZRCT(:,:)>LIMAP%XRTMIN(2) .OR. &
        ZRRT(:,:)>LIMAP%XRTMIN(3) .OR. &
        ZRIT(:,:)>LIMAP%XRTMIN(4) .OR. &
        ZRST(:,:)>LIMAP%XRTMIN(5) .OR. &
        ZRGT(:,:)>LIMAP%XRTMIN(6) .OR. &
        ZRHT(:,:)>LIMAP%XRTMIN(7)
ELSE
   LLMICRO(:,:)=ZRCT(:,:)>LIMAP%XRTMIN(2) .OR. &
        ZRRT(:,:)>LIMAP%XRTMIN(3) .OR. &
        ZRIT(:,:)>LIMAP%XRTMIN(4) .OR. &
        ZRST(:,:)>LIMAP%XRTMIN(5) .OR. &
        ZRGT(:,:)>LIMAP%XRTMIN(6)
ENDIF
! PHLC_HRC should never be larger than PRCT
!PCLDFR(:,:)=MAX(PCLDFR(:,:),PHLC_HCF(:,:))
!PHLC_HRC(:,:)=MIN(PHLC_HRC(:,:),ZRCT(:,:))
!PICEFR(:,:)=MAX(PICEFR(:,:),PHLI_HCF(:,:))
!PHLI_HRI(:,:)=MIN(PHLI_HRI(:,:),ZRIT(:,:))
IF ( NEBN%LSUBG_COND ) THEN
   IF (PRESENT(PHLC_HRC)) THEN
      ZHLC_HRC(:,:)=PHLC_HRC(:,:)
      ZHLC_HCF(:,:)=PHLC_HCF(:,:)
      ZHLI_HRI(:,:)=PHLI_HRI(:,:)
      ZHLI_HCF(:,:)=PHLI_HCF(:,:)
      DO IK = D%NKTB, D%NKTE     
         DO II=D%NIJB, D%NIJE
            ZHLC_LRC(II, IK) = ZRCT(II, IK) - ZHLC_HRC(II, IK)
            ZHLI_LRI(II, IK) = ZRIT(II, IK) - ZHLI_HRI(II, IK)
            IF(ZRCT(II, IK)>0.) THEN
               ZHLC_LCF(II, IK) = PCLDFR(II, IK)- ZHLC_HCF(II, IK)
            ELSE
               ZHLC_LCF(II, IK)=0.
            ENDIF
            IF(ZRIT(II, IK)>0.) THEN
               ZHLI_LCF(II, IK) = PICEFR(II, IK)- ZHLI_HCF(II, IK)
            ELSE
               ZHLI_LCF(II, IK)=0.
            ENDIF
         ENDDO
      ENDDO
   ELSE
      ZHLC_LRC(:,:)=0.
      ZHLI_LRI(:,:)=0.
      ZHLC_LCF(:,:)=0.
      ZHLI_LCF(:,:)=0.
      ZHLC_HRC(:,:)=ZRCT(:,:)
      ZHLC_HCF(:,:)=PCLDFR(:,:)
      ZHLI_HRI(:,:)=ZRIT(:,:)
      ZHLI_HCF(:,:)=PICEFR(:,:)
   END IF
  CALL LIMA_COMPUTE_PDF(CST, LIMAP, D%NIJT*(D%NKTE-D%NKTB+1), LIMAP%CSUBG_AUCV_RC, LIMAP%CSUBG_AUCV_RI, LIMAP%CSUBG_PR_PDF,&
                        LLMICRO(:,D%NKTB:D%NKTE), PRHODREF(:,D%NKTB:D%NKTE), ZRCT(:,D%NKTB:D%NKTE), ZRIT(:,D%NKTB:D%NKTE), &
                        PCLDFR(:,D%NKTB:D%NKTE), ZT(:,D%NKTB:D%NKTE), ZSIGMA_RC(:,D%NKTB:D%NKTE), &
                        ZHLC_HCF(:,D%NKTB:D%NKTE), ZHLC_LCF(:,D%NKTB:D%NKTE), ZHLC_HRC(:,D%NKTB:D%NKTE), ZHLC_LRC(:,D%NKTB:D%NKTE), &
                        ZHLI_HCF(:,D%NKTB:D%NKTE), ZHLI_LCF(:,D%NKTB:D%NKTE), ZHLI_HRI(:,D%NKTB:D%NKTE), ZHLI_LRI(:,D%NKTB:D%NKTE), &
                        PPRCFR(:,D%NKTB:D%NKTE))
  CALL LIMA_RAINFR_VERT(D, LIMAP, PPRCFR, ZRRT, ZRST, ZRGT, ZRHT)
ELSE
   PPRCFR(:,:)=1.
   ZHLC_LRC(:,:)=0.
   ZHLI_LRI(:,:)=0.
   ZHLC_LCF(:,:)=0.
   ZHLI_LCF(:,:)=0.
   ZHLC_HRC(:,:)=ZRCT(:,:)
   ZHLI_HRI(:,:)=ZRIT(:,:)
   ZHLC_HCF(:,:)=1.
   ZHLI_HCF(:,:)=1.
ENDIF

!-------------------------------------------------------------------------------
!
!*       2.     LOOP
!               ----
!
!
! Maximum number of iterations
INB_ITER_MAX=LIMAP%NMAXITER
IF(LIMAP%XTSTEP_TS/=0.)THEN
  INB_ITER_MAX=MAX(1, INT(PTSTEP/LIMAP%XTSTEP_TS)) !At least the number of iterations needed for the time-splitting
  ZTSTEP=PTSTEP/INB_ITER_MAX
  INB_ITER_MAX=MAX(LIMAP%NMAXITER, INB_ITER_MAX) !Fot the case LIMAP%XMRSTEP/=0. at the same time
ENDIF
IITER(:,:)=0
ZTIME(:,:)=0. ! Current integration time (all points may have a different integration time)
!
! Begin the huge time splitting loop
!
ZRT_SUM(D%NIJB:D%NIJE,:) = ZRCT(D%NIJB:D%NIJE,:) + ZRRT(D%NIJB:D%NIJE,:) + ZRIT(D%NIJB:D%NIJE,:) + &
                         ZRST(D%NIJB:D%NIJE,:) + ZRGT(D%NIJB:D%NIJE,:) + ZRHT(D%NIJB:D%NIJE,:)
WHERE (ZRT_SUM(D%NIJB:D%NIJE,:)<LIMAP%XRTMIN(2)) ZTIME(D%NIJB:D%NIJE,:)=PTSTEP ! no need to treat hydrometeor-free point
!
DO WHILE(ANY(ZTIME(D%NIJB:D%NIJE,D%NKTB:D%NKTE)<PTSTEP))
   !
   IF(LIMAP%XMRSTEP/=0.) THEN
      ! In this case we need to remember the mixing ratios used to compute the tendencies
      ! because when mixing ratio has evolved more than a threshold, we must re-compute tendecies
      Z0RVT(D%NIJB:D%NIJE,:)=ZRVT(D%NIJB:D%NIJE,:)
      Z0RCT(D%NIJB:D%NIJE,:)=ZRCT(D%NIJB:D%NIJE,:)
      Z0RRT(D%NIJB:D%NIJE,:)=ZRRT(D%NIJB:D%NIJE,:)
      Z0RIT(D%NIJB:D%NIJE,:)=ZRIT(D%NIJB:D%NIJE,:)
      Z0RST(D%NIJB:D%NIJE,:)=ZRST(D%NIJB:D%NIJE,:)
      Z0RGT(D%NIJB:D%NIJE,:)=ZRGT(D%NIJB:D%NIJE,:)
      Z0RHT(D%NIJB:D%NIJE,:)=ZRHT(D%NIJB:D%NIJE,:)
   ENDIF
   !
   IF(LIMAP%XTSTEP_TS/=0.) THEN
      ! In this case we need to remember the time when tendencies were computed
      ! because when time has evolved more than a limit, we must re-compute tendecies
      ZTIME_LASTCALL(D%NIJB:D%NIJE,:)=ZTIME(D%NIJB:D%NIJE,:)
   ENDIF
   !
   GLCOMPUTE(:,:)=.FALSE.
   GLCOMPUTE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = ZTIME(D%NIJB:D%NIJE,D%NKTB:D%NKTE)<PTSTEP ! Compuation only for points for which integration time has not reached the timestep
   WHERE(GLCOMPUTE(:,:))
      IITER(:,:)=IITER(:,:)+1
   END WHERE
   ! 
   DO WHILE(ANY(GLCOMPUTE(:,:))) ! Loop to adjust tendencies when we cross the 0°C or when a species disappears

      !
      ! Packing variables to run computations only where necessary
      !
      IPACK = COUNT(GLCOMPUTE)
      ALLOCATE(I1(IPACK))
      ALLOCATE(I3(IPACK))
      IPACK = COUNTJV(GLCOMPUTE,I1,I3)
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
      ALLOCATE(ZCST1D(IPACK))
      ALLOCATE(ZCGT1D(IPACK))
      ALLOCATE(ZCHT1D(IPACK))
      ALLOCATE(ZIFNN1D(IPACK,LIMAP%NMOD_IFN))
      ALLOCATE(ZEVAP1D(IPACK))
      ALLOCATE(ZTIME1D(IPACK))
      ALLOCATE(GLCOMPUTE1D(IPACK))
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
      ALLOCATE(ZHLC_HRC1D(IPACK))
      ALLOCATE(ZHLC_HCF1D(IPACK))
      ALLOCATE(ZHLI_HRI1D(IPACK))
      ALLOCATE(ZHLI_HCF1D(IPACK))
      ALLOCATE(ZHLC_LRC1D(IPACK))
      ALLOCATE(ZHLC_LCF1D(IPACK))
      ALLOCATE(ZHLI_LRI1D(IPACK))
      ALLOCATE(ZHLI_LCF1D(IPACK))
      ALLOCATE(ZLATHAM_IAGGS(IPACK))
      IF (LIMAP%LCRYSTAL_SHAPE) THEN
        ALLOCATE(ZCIT1D_SHAPE(IPACK,LIMAP%NNB_CRYSTAL_SHAPE))
      ELSE
        ALLOCATE(ZCIT1D_SHAPE(0,0))
      END IF
      DO II=1,IPACK
         ZRHODREF1D(II)       = PRHODREF(I1(II),I3(II))
         ZEXNREF1D(II)        = PEXNREF(I1(II),I3(II))
         ZEXN1D(II)           = ZEXN(I1(II),I3(II))
         ZP1D(II)             = PPABST(I1(II),I3(II))
         ZTHT1D(II)           = ZTHT(I1(II),I3(II))
         ZRVT1D(II)           = ZRVT(I1(II),I3(II))
         ZRCT1D(II)           = ZRCT(I1(II),I3(II))
         ZRRT1D(II)           = ZRRT(I1(II),I3(II))
         ZRIT1D(II)           = ZRIT(I1(II),I3(II))
         ZRST1D(II)           = ZRST(I1(II),I3(II))
         ZRGT1D(II)           = ZRGT(I1(II),I3(II))
         ZRHT1D(II)           = ZRHT(I1(II),I3(II))
         ZCCT1D(II)           = ZCCT(I1(II),I3(II))
         ZCRT1D(II)           = ZCRT(I1(II),I3(II))
         ZCIT1D(II)           = ZCIT(I1(II),I3(II))
         ZCST1D(II)           = ZCST(I1(II),I3(II))
         ZCGT1D(II)           = ZCGT(I1(II),I3(II))
         ZCHT1D(II)           = ZCHT(I1(II),I3(II))
         ZIFNN1D(II,:)        = ZIFNNT(I1(II),I3(II),:)
         ZEVAP1D(II)          = PEVAP3D(I1(II),I3(II))
         ZTIME1D(II)          = ZTIME(I1(II),I3(II))         
         GLCOMPUTE1D(II)      = GLCOMPUTE(I1(II),I3(II))         
         IITER1D(II)          = IITER(I1(II),I3(II))         
         ZTIME_LASTCALL1D(II) = ZTIME_LASTCALL(I1(II),I3(II))         
         Z0RVT1D(II)          = Z0RVT(I1(II),I3(II))
         Z0RCT1D(II)          = Z0RCT(I1(II),I3(II))
         Z0RRT1D(II)          = Z0RRT(I1(II),I3(II))
         Z0RIT1D(II)          = Z0RIT(I1(II),I3(II))
         Z0RST1D(II)          = Z0RST(I1(II),I3(II))
         Z0RGT1D(II)          = Z0RGT(I1(II),I3(II))
         Z0RHT1D(II)          = Z0RHT(I1(II),I3(II))
         ZCF1D(II)            = PCLDFR(I1(II),I3(II))
         ZIF1D(II)            = PICEFR(I1(II),I3(II))
         ZPF1D(II)            = PPRCFR(I1(II),I3(II))
         ZHLC_HRC1D(II)       = ZHLC_HRC(I1(II),I3(II))
         ZHLC_HCF1D(II)       = ZHLC_HCF(I1(II),I3(II))
         ZHLI_HRI1D(II)       = ZHLI_HRI(I1(II),I3(II))
         ZHLI_HCF1D(II)       = ZHLI_HCF(I1(II),I3(II))
         ZHLC_LRC1D(II)       = ZHLC_LRC(I1(II),I3(II))
         ZHLC_LCF1D(II)       = ZHLC_LCF(I1(II),I3(II))
         ZHLI_LRI1D(II)       = ZHLI_LRI(I1(II),I3(II))
         ZHLI_LCF1D(II)       = ZHLI_LCF(I1(II),I3(II))
         IF (OELEC) THEN
           ZLATHAM_IAGGS(II) = PLATHAM_IAGGS(I1(II),I3(II))
         ELSE
           ZLATHAM_IAGGS(II) = 1.0 
         END IF
         IF (LIMAP%LCRYSTAL_SHAPE) THEN
           DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
             ZCIT1D_SHAPE(II,ISH) = ZCIT_SHAPE(I1(II),I3(II),ISH)
           END DO
         END IF
      END DO
      !
      WHERE(ZCF1D(:)<1.E-10 .AND. ZRCT1D(:)>LIMAP%XRTMIN(2) .AND. ZCCT1D(:)>LIMAP%XCTMIN(2)) ZCF1D(:)=1.
      WHERE(ZIF1D(:)<1.E-10 .AND. ZRIT1D(:)>LIMAP%XRTMIN(4) .AND. ZCIT1D(:)>LIMAP%XCTMIN(4)) ZIF1D(:)=1.
      WHERE(ZPF1D(:)<1.E-10 .AND. ZRRT1D(:)>LIMAP%XRTMIN(3) .AND. ZCRT1D(:)>LIMAP%XCTMIN(3)) ZPF1D(:)=1.
      WHERE(ZPF1D(:)<1.E-10 .AND. ZRST1D(:)>LIMAP%XRTMIN(5) .AND. ZCST1D(:)>LIMAP%XCTMIN(5)) ZPF1D(:)=1.
      WHERE(ZPF1D(:)<1.E-10 .AND. ZRGT1D(:)>LIMAP%XRTMIN(6) .AND. ZCGT1D(:)>LIMAP%XCTMIN(6)) ZPF1D(:)=1.
      WHERE(ZPF1D(:)<1.E-10 .AND. ZRHT1D(:)>LIMAP%XRTMIN(7) .AND. ZCHT1D(:)>LIMAP%XCTMIN(7)) ZPF1D(:)=1.
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
      ALLOCATE(ZA_CS(IPACK))              ;  ZA_CS(:) = 0.
      ALLOCATE(ZA_CG(IPACK))              ;  ZA_CG(:) = 0.
      ALLOCATE(ZA_CH(IPACK))              ;  ZA_CH(:) = 0.
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
      ALLOCATE(ZB_CS(IPACK))              ;  ZB_CS(:) = 0.
      ALLOCATE(ZB_CG(IPACK))              ;  ZB_CG(:) = 0.
      ALLOCATE(ZB_CH(IPACK))              ;  ZB_CH(:) = 0.
      ALLOCATE(ZB_IFNN(IPACK,LIMAP%NMOD_IFN))   ;  ZB_IFNN(:,:) = 0.
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
      ALLOCATE(Z_CR_EVAP(IPACK))          ; Z_CR_EVAP(:) = 0.
      ALLOCATE(Z_RI_CNVI(IPACK))          ; Z_RI_CNVI(:) = 0.
      ALLOCATE(Z_CI_CNVI(IPACK))          ; Z_CI_CNVI(:) = 0.
      ALLOCATE(Z_TH_DEPS(IPACK))          ; Z_TH_DEPS(:) = 0.
      ALLOCATE(Z_RS_DEPS(IPACK))          ; Z_RS_DEPS(:) = 0.
      ALLOCATE(Z_TH_DEPI(IPACK))          ; Z_TH_DEPI(:) = 0.
      ALLOCATE(Z_RI_DEPI(IPACK))          ; Z_RI_DEPI(:) = 0.
      ALLOCATE(Z_RI_CNVS(IPACK))          ; Z_RI_CNVS(:) = 0.
      ALLOCATE(Z_CI_CNVS(IPACK))          ; Z_CI_CNVS(:) = 0.
      ALLOCATE(Z_CS_SSC(IPACK))           ; Z_CS_SSC(:) = 0.
      ALLOCATE(Z_CI_ISC(IPACK))           ; Z_CI_ISC(:) = 0.
      ALLOCATE(Z_RI_AGGS(IPACK))          ; Z_RI_AGGS(:) = 0.
      ALLOCATE(Z_CI_AGGS(IPACK))          ; Z_CI_AGGS(:) = 0.
      ALLOCATE(Z_TH_DEPG(IPACK))          ; Z_TH_DEPG(:) = 0.
      ALLOCATE(Z_RG_DEPG(IPACK))          ; Z_RG_DEPG(:) = 0.
      ALLOCATE(Z_TH_BERFI(IPACK))         ; Z_TH_BERFI(:) = 0.
      ALLOCATE(Z_RC_BERFI(IPACK))         ; Z_RC_BERFI(:) = 0.
      ALLOCATE(Z_TH_RIM(IPACK))           ; Z_TH_RIM(:) = 0.
      ALLOCATE(Z_CC_RIM(IPACK))           ; Z_CC_RIM(:) = 0.
      ALLOCATE(Z_CS_RIM(IPACK))           ; Z_CS_RIM(:) = 0.
      ALLOCATE(Z_RC_RIMSS(IPACK))         ; Z_RC_RIMSS = 0.
      ALLOCATE(Z_RC_RIMSG(IPACK))         ; Z_RC_RIMSG = 0.
      ALLOCATE(Z_RS_RIMCG(IPACK))         ; Z_RS_RIMCG = 0.
      ALLOCATE(Z_RI_HMS(IPACK))           ; Z_RI_HMS(:) = 0.
      ALLOCATE(Z_CI_HMS(IPACK))           ; Z_CI_HMS(:) = 0.
      ALLOCATE(Z_RS_HMS(IPACK))           ; Z_RS_HMS(:) = 0.
      ALLOCATE(Z_TH_ACC(IPACK))           ; Z_TH_ACC(:) = 0.
      ALLOCATE(Z_CR_ACC(IPACK))           ; Z_CR_ACC(:) = 0.
      ALLOCATE(Z_CS_ACC(IPACK))           ; Z_CS_ACC(:) = 0.
      ALLOCATE(Z_RR_ACCSS(IPACK))         ; Z_RR_ACCSS = 0.
      ALLOCATE(Z_RR_ACCSG(IPACK))         ; Z_RR_ACCSG = 0.
      ALLOCATE(Z_RS_ACCRG(IPACK))         ; Z_RS_ACCRG = 0.
      ALLOCATE(Z_RS_CMEL(IPACK))          ; Z_RS_CMEL(:) = 0.
      ALLOCATE(Z_CS_CMEL(IPACK))          ; Z_CS_CMEL(:) = 0.
      ALLOCATE(Z_TH_CFRZ(IPACK))          ; Z_TH_CFRZ(:) = 0.
      ALLOCATE(Z_RR_CFRZ(IPACK))          ; Z_RR_CFRZ(:) = 0.
      ALLOCATE(Z_CR_CFRZ(IPACK))          ; Z_CR_CFRZ(:) = 0.
      ALLOCATE(Z_RI_CFRZ(IPACK))          ; Z_RI_CFRZ(:) = 0.
      ALLOCATE(Z_CI_CFRZ(IPACK))          ; Z_CI_CFRZ(:) = 0.
      ALLOCATE(Z_RI_CIBU(IPACK))          ; Z_RI_CIBU(:) = 0.
      ALLOCATE(Z_CI_CIBU(IPACK))          ; Z_CI_CIBU(:) = 0.
      ALLOCATE(Z_RI_RDSF(IPACK))          ; Z_RI_RDSF(:) = 0.
      ALLOCATE(Z_CI_RDSF(IPACK))          ; Z_CI_RDSF(:) = 0.
      ALLOCATE(Z_TH_WETG(IPACK))          ; Z_TH_WETG(:) = 0.
      ALLOCATE(Z_RC_WETG(IPACK))          ; Z_RC_WETG(:) = 0.
      ALLOCATE(Z_CC_WETG(IPACK))          ; Z_CC_WETG(:) = 0.
      ALLOCATE(Z_RR_WETG(IPACK))          ; Z_RR_WETG(:) = 0.
      ALLOCATE(Z_CR_WETG(IPACK))          ; Z_CR_WETG(:) = 0.
      ALLOCATE(Z_RI_WETG(IPACK))          ; Z_RI_WETG(:) = 0.
      ALLOCATE(Z_CI_WETG(IPACK))          ; Z_CI_WETG(:) = 0.
      ALLOCATE(Z_RS_WETG(IPACK))          ; Z_RS_WETG(:) = 0.
      ALLOCATE(Z_CS_WETG(IPACK))          ; Z_CS_WETG(:) = 0.
      ALLOCATE(Z_RG_WETG(IPACK))          ; Z_RG_WETG(:) = 0.
      ALLOCATE(Z_CG_WETG(IPACK))          ; Z_CG_WETG(:) = 0.
      ALLOCATE(Z_RH_WETG(IPACK))          ; Z_RH_WETG(:) = 0.
      ALLOCATE(Z_TH_DRYG(IPACK))          ; Z_TH_DRYG(:) = 0.
      ALLOCATE(Z_RC_DRYG(IPACK))          ; Z_RC_DRYG(:) = 0.
      ALLOCATE(Z_CC_DRYG(IPACK))          ; Z_CC_DRYG(:) = 0.
      ALLOCATE(Z_RR_DRYG(IPACK))          ; Z_RR_DRYG(:) = 0.
      ALLOCATE(Z_CR_DRYG(IPACK))          ; Z_CR_DRYG(:) = 0.
      ALLOCATE(Z_RI_DRYG(IPACK))          ; Z_RI_DRYG(:) = 0.
      ALLOCATE(Z_CI_DRYG(IPACK))          ; Z_CI_DRYG(:) = 0.
      ALLOCATE(Z_RS_DRYG(IPACK))          ; Z_RS_DRYG(:) = 0.
      ALLOCATE(Z_CS_DRYG(IPACK))          ; Z_CS_DRYG(:) = 0.
      ALLOCATE(Z_RG_DRYG(IPACK))          ; Z_RG_DRYG(:) = 0.
      ALLOCATE(Z_RI_HMG(IPACK))           ; Z_RI_HMG(:) = 0.
      ALLOCATE(Z_CI_HMG(IPACK))           ; Z_CI_HMG(:) = 0.
      ALLOCATE(Z_RG_HMG(IPACK))           ; Z_RG_HMG(:) = 0.
      ALLOCATE(Z_TH_GMLT(IPACK))          ; Z_TH_GMLT(:) = 0.
      ALLOCATE(Z_RR_GMLT(IPACK))          ; Z_RR_GMLT(:) = 0.
      ALLOCATE(Z_CR_GMLT(IPACK))          ; Z_CR_GMLT(:) = 0.
      ALLOCATE(Z_CG_GMLT(IPACK))          ; Z_CG_GMLT(:) = 0.
      ALLOCATE(Z_TH_DEPH(IPACK))          ; Z_TH_DEPH(:) = 0.
      ALLOCATE(Z_RH_DEPH(IPACK))          ; Z_RH_DEPH(:) = 0.
      ALLOCATE(Z_TH_WETH(IPACK))          ; Z_TH_WETH(:) = 0.
      ALLOCATE(Z_RC_WETH(IPACK))          ; Z_RC_WETH(:) = 0.
      ALLOCATE(Z_CC_WETH(IPACK))          ; Z_CC_WETH(:) = 0.
      ALLOCATE(Z_RR_WETH(IPACK))          ; Z_RR_WETH(:) = 0.
      ALLOCATE(Z_CR_WETH(IPACK))          ; Z_CR_WETH(:) = 0.
      ALLOCATE(Z_RI_WETH(IPACK))          ; Z_RI_WETH(:) = 0.
      ALLOCATE(Z_CI_WETH(IPACK))          ; Z_CI_WETH(:) = 0.
      ALLOCATE(Z_RS_WETH(IPACK))          ; Z_RS_WETH(:) = 0.
      ALLOCATE(Z_CS_WETH(IPACK))          ; Z_CS_WETH(:) = 0.
      ALLOCATE(Z_RG_WETH(IPACK))          ; Z_RG_WETH(:) = 0.
      ALLOCATE(Z_CG_WETH(IPACK))          ; Z_CG_WETH(:) = 0.
      ALLOCATE(Z_RH_WETH(IPACK))          ; Z_RH_WETH(:) = 0.
      ALLOCATE(Z_RG_COHG(IPACK))          ; Z_RG_COHG(:) = 0.
      ALLOCATE(Z_CG_COHG(IPACK))          ; Z_CG_COHG(:) = 0.
      ALLOCATE(Z_TH_HMLT(IPACK))          ; Z_TH_HMLT(:) = 0.
      ALLOCATE(Z_RR_HMLT(IPACK))          ; Z_RR_HMLT(:) = 0.
      ALLOCATE(Z_CR_HMLT(IPACK))          ; Z_CR_HMLT(:) = 0.
      ALLOCATE(Z_CH_HMLT(IPACK))          ; Z_CH_HMLT(:) = 0.

      ALLOCATE(Z_RV_CORR2(IPACK))         ; Z_RV_CORR2(:) = 0.
      ALLOCATE(Z_RC_CORR2(IPACK))         ; Z_RC_CORR2(:) = 0.
      ALLOCATE(Z_RR_CORR2(IPACK))         ; Z_RR_CORR2(:) = 0.
      ALLOCATE(Z_RI_CORR2(IPACK))         ; Z_RI_CORR2(:) = 0.
      ALLOCATE(Z_CC_CORR2(IPACK))         ; Z_CC_CORR2(:) = 0.
      ALLOCATE(Z_CR_CORR2(IPACK))         ; Z_CR_CORR2(:) = 0.
      ALLOCATE(Z_CI_CORR2(IPACK))         ; Z_CI_CORR2(:) = 0.
!
      IF (LIMAP%LCRYSTAL_SHAPE) THEN
        ALLOCATE(ZA_CI_SHAPE(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; ZA_CI_SHAPE(:,:) = 0.
        ALLOCATE(ZB_CI_SHAPE(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; ZB_CI_SHAPE(:,:) = 0.
        ALLOCATE(Z_SHCI_IMLT(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_IMLT(:,:) = 0.
        ALLOCATE(Z_SHCI_HONC(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_HONC(:,:) = 0.
        ALLOCATE(Z_SHCI_CNVI(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_CNVI(:,:) = 0.
        ALLOCATE(Z_SHCI_HACH(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_HACH(:,:) = 0.
        ALLOCATE(Z_SHCI_CNVS(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_CNVS(:,:) = 0.
        ALLOCATE(Z_SHCI_AGGS(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_AGGS(:,:) = 0.
        ALLOCATE(Z_SHCI_HMS(IPACK,LIMAP%NNB_CRYSTAL_SHAPE))  ; Z_SHCI_HMS(:,:)  = 0.
        ALLOCATE(Z_SHCI_CFRZ(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_CFRZ(:,:) = 0.
        ALLOCATE(Z_SHCI_CIBU(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_CIBU(:,:) = 0.
        ALLOCATE(Z_SHCI_RDSF(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_RDSF(:,:) = 0.
        ALLOCATE(Z_SHCI_WETG(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_WETG(:,:) = 0.
        ALLOCATE(Z_SHCI_DRYG(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_DRYG(:,:) = 0.
        ALLOCATE(Z_SHCI_HMG(IPACK,LIMAP%NNB_CRYSTAL_SHAPE))  ; Z_SHCI_HMG(:,:)  = 0.
        ALLOCATE(Z_SHCI_CORR2(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)); Z_SHCI_CORR2(:,:) = 0.
        ALLOCATE(Z_SHCI_ISC(IPACK,LIMAP%NNB_CRYSTAL_SHAPE))  ; Z_SHCI_ISC(:,:)  = 0.
        ALLOCATE(Z_SHCI_ISCS(IPACK,LIMAP%NNB_CRYSTAL_SHAPE)) ; Z_SHCI_ISCS(:,:) = 0.
        ALLOCATE(Z_SHRI_ISCS(IPACK))                   ; Z_SHRI_ISCS(:) = 0.
      ELSE
        ALLOCATE(ZA_CI_SHAPE(0,0))
        ALLOCATE(ZB_CI_SHAPE(0,0))
        ALLOCATE(Z_SHCI_IMLT(0,0))
        ALLOCATE(Z_SHCI_HONC(0,0))
        ALLOCATE(Z_SHCI_CNVI(0,0))
        ALLOCATE(Z_SHCI_HACH(0,0))
        ALLOCATE(Z_SHCI_CNVS(0,0))
        ALLOCATE(Z_SHCI_AGGS(0,0))
        ALLOCATE(Z_SHCI_HMS(0,0))
        ALLOCATE(Z_SHCI_CFRZ(0,0))
        ALLOCATE(Z_SHCI_CIBU(0,0))
        ALLOCATE(Z_SHCI_RDSF(0,0))
        ALLOCATE(Z_SHCI_WETG(0,0))
        ALLOCATE(Z_SHCI_DRYG(0,0))
        ALLOCATE(Z_SHCI_HMG(0,0))
        ALLOCATE(Z_SHCI_CORR2(0,0))
        ALLOCATE(Z_SHCI_ISC(0,0))
        ALLOCATE(Z_SHCI_ISCS(0,0))
        ALLOCATE(Z_SHRI_ISCS(0))
      END IF
      !
      !***       4.1 Tendencies computation
      !
      CALL LIMA_INST_PROCS (CST, LIMAP, LIMAW, IPACK, PTSTEP, GLCOMPUTE1D,      &
                            ZEXNREF1D, ZP1D,                                    &
                            ZTHT1D, ZRVT1D, ZRCT1D, ZRRT1D, ZRIT1D, ZRST1D, ZRGT1D, &
                            ZCCT1D, ZCRT1D, ZCIT1D, ZCIT1D_SHAPE,               &
                            ZIFNN1D,                                            &
                            Z_CR_BRKU,                                          & ! spontaneous break up of drops (BRKU) : Nr
                            Z_TH_HONR, Z_RR_HONR, Z_CR_HONR,                    & ! rain drops homogeneous freezing (HONR) : rr, Nr, rg=-rr, th
                            Z_TH_IMLT, Z_RC_IMLT, Z_CC_IMLT,                    & ! ice melting (IMLT) : rc, Nc, ri=-rc, Ni=-Nc, th, IFNF, IFNA
                            Z_SHCI_IMLT,                                        &
                            ZB_TH, ZB_RV, ZB_RC, ZB_RR, ZB_RI, ZB_RG,           &
                            ZB_CC, ZB_CR, ZB_CI, ZB_CI_SHAPE, ZB_CG,            &
                            ZB_IFNN,                                            &
                            ZCF1D, ZIF1D, ZPF1D                                 )
      
      CALL LIMA_TENDENCIES (CST, LIMAP, LIMAW, LIMAC, LIMAM, IPACK, PTSTEP, GLCOMPUTE1D, &
                            ZEXNREF1D, ZRHODREF1D, ZP1D, ZTHT1D,                    &
                            ZRVT1D, ZRCT1D, ZRRT1D, ZRIT1D, ZRST1D, ZRGT1D, ZRHT1D, &
                            ZCCT1D, ZCRT1D, ZCIT1D, ZCIT1D_SHAPE, ZCST1D, ZCGT1D, ZCHT1D, &
                            Z_TH_HONC, Z_RC_HONC, Z_CC_HONC, Z_SHCI_HONC,           & 
                            Z_CC_SELF,                                              & 
                            Z_RC_AUTO, Z_CC_AUTO, Z_CR_AUTO,                        & 
                            Z_RC_ACCR, Z_CC_ACCR,                                   & 
                            Z_CR_SCBU,                                              & 
                            Z_TH_EVAP, Z_RR_EVAP, Z_CR_EVAP,                        & 
                            Z_RI_CNVI, Z_CI_CNVI, Z_SHCI_CNVI,                      & 
                            Z_TH_DEPS, Z_RS_DEPS,                                   & 
                            Z_TH_DEPI, Z_RI_DEPI, Z_SHCI_HACH,                      & 
                            Z_RI_CNVS, Z_CI_CNVS, Z_SHCI_CNVS,                      &
                            Z_CS_SSC,  Z_CI_ISC, Z_SHCI_ISC, Z_SHRI_ISCS, Z_SHCI_ISCS, &
                            Z_RI_AGGS, Z_CI_AGGS, Z_SHCI_AGGS,                      &
                            Z_TH_DEPG, Z_RG_DEPG,                                   & 
                            Z_TH_BERFI, Z_RC_BERFI,                                 & 
                            Z_TH_RIM, Z_CC_RIM, Z_CS_RIM, Z_RC_RIMSS, Z_RC_RIMSG, Z_RS_RIMCG, &
                            Z_RI_HMS, Z_CI_HMS, Z_RS_HMS, Z_SHCI_HMS,                         &
                            Z_TH_ACC, Z_CR_ACC, Z_CS_ACC, Z_RR_ACCSS, Z_RR_ACCSG, Z_RS_ACCRG, &
                            Z_RS_CMEL, Z_CS_CMEL,                                   & 
                            Z_TH_CFRZ, Z_RR_CFRZ, Z_CR_CFRZ, Z_RI_CFRZ, Z_CI_CFRZ, Z_SHCI_CFRZ,& 
                            Z_RI_CIBU, Z_CI_CIBU, Z_SHCI_CIBU,                      & 
                            Z_RI_RDSF, Z_CI_RDSF, Z_SHCI_RDSF,                      & 
                            Z_TH_WETG, Z_RC_WETG, Z_CC_WETG, Z_RR_WETG, Z_CR_WETG,  & 
                            Z_RI_WETG, Z_CI_WETG, Z_RS_WETG, Z_CS_WETG, Z_RG_WETG, Z_CG_WETG,  &
                            Z_RH_WETG, Z_SHCI_WETG,                                 & 
                            Z_TH_DRYG, Z_RC_DRYG, Z_CC_DRYG, Z_RR_DRYG, Z_CR_DRYG,  & 
                            Z_RI_DRYG, Z_CI_DRYG, Z_RS_DRYG, Z_CS_DRYG, Z_RG_DRYG, Z_SHCI_DRYG,& 
                            Z_RI_HMG, Z_CI_HMG, Z_RG_HMG, Z_SHCI_HMG,               & 
                            Z_TH_GMLT, Z_RR_GMLT, Z_CR_GMLT, Z_CG_GMLT,             &
                            Z_TH_DEPH, Z_RH_DEPH,                                   &
                            Z_TH_WETH, Z_RC_WETH, Z_CC_WETH, Z_RR_WETH, Z_CR_WETH,  &
                            Z_RI_WETH, Z_CI_WETH, Z_RS_WETH, Z_CS_WETH, Z_RG_WETH,  &
                            Z_CG_WETH, Z_RH_WETH, &
                            Z_RG_COHG, Z_CG_COHG,                                   &
                            Z_TH_HMLT, Z_RR_HMLT, Z_CR_HMLT, Z_CH_HMLT,             &
                            ZA_TH, ZA_RV, ZA_RC, ZA_CC, ZA_RR, ZA_CR,               &
                            ZA_RI, ZA_CI, ZA_CI_SHAPE, ZA_RS, ZA_CS, ZA_RG, ZA_CG,  &
                            ZA_RH, ZA_CH,                                           &
                            ZEVAP1D,                                                &
                            ZCF1D, ZIF1D, ZPF1D,                                    &
                            ZHLC_HCF1D, ZHLC_LCF1D, ZHLC_HRC1D, ZHLC_LRC1D,         &
                            ZHLI_HCF1D, ZHLI_LCF1D, ZHLI_HRI1D, ZHLI_LRI1D,         &
                            ZLATHAM_IAGGS                                           )

      !
      !***       4.2 Integration time
      !
      ! If we can, we will use these tendecies until the end of the timestep
      ZMAXTIME(:)=PTSTEP-ZTIME1D(:) ! Remaining time until the end of the timestep

      ! We need to adjust tendencies when temperature reaches 0
      IF(LIMAP%LFEEDBACKT) THEN
         !Is ZB_TH enough to change temperature sign?
         WHERE( ((ZTHT1D(:) - CST%XTT/ZEXN1D(:)) * (ZTHT1D(:) + ZB_TH(:) - CST%XTT/ZEXN1D(:))) < 0. )
            ZMAXTIME(:)=0.
         ENDWHERE
         !Can ZA_TH make temperature change of sign?
         ZTIME_THRESHOLD(:)=-1.
         WHERE(ABS(ZA_TH(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(CST%XTT/ZEXN1D(:) - ZB_TH(:) - ZTHT1D(:))/ZA_TH(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>0.)
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
         ENDWHERE
      ENDIF

      ! We need to adjust tendencies when a species disappears
      ! When a species is missing, only the external tendencies can be negative (and we must keep track of it)
      WHERE(ZA_RV(:)<-1.E-20 .AND. ZRVT1D(:)>LIMAP%XRTMIN(1))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RV(:)+ZRVT1D(:))/ZA_RV(:))
      END WHERE
      WHERE(ZA_RC(:)<-1.E-20 .AND. ZRCT1D(:)>LIMAP%XRTMIN(2))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RC(:)+ZRCT1D(:))/ZA_RC(:))
      END WHERE
      IF (LIMAP%NMOM_C.GE.2) THEN
         WHERE(ZA_CC(:)<-1.E-20 .AND. ZCCT1D(:)>LIMAP%XCTMIN(2))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_CC(:)+ZCCT1D(:))/ZA_CC(:))
         END WHERE
      ENDIF
      WHERE(ZA_RR(:)<-1.E-20 .AND. ZRRT1D(:)>LIMAP%XRTMIN(3))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RR(:)+ZRRT1D(:))/ZA_RR(:))
      END WHERE
      IF (LIMAP%NMOM_R.GE.2) THEN
         WHERE(ZA_CR(:)<-1.E-20 .AND. ZCRT1D(:)>LIMAP%XCTMIN(3))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_CR(:)+ZCRT1D(:))/ZA_CR(:))
         END WHERE
      ENDIF
      WHERE(ZA_RI(:)<-1.E-20 .AND. ZRIT1D(:)>LIMAP%XRTMIN(4))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RI(:)+ZRIT1D(:))/ZA_RI(:))
      END WHERE
      IF (LIMAP%NMOM_I.GE.2) THEN
         WHERE(ZA_CI(:)<-1.E-20 .AND. ZCIT1D(:)>LIMAP%XCTMIN(4))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_CI(:)+ZCIT1D(:))/ZA_CI(:))
         END WHERE
      ENDIF
      WHERE(ZA_RS(:)<-1.E-20 .AND. ZRST1D(:)>LIMAP%XRTMIN(5))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RS(:)+ZRST1D(:))/ZA_RS(:))
      END WHERE
      IF (LIMAP%NMOM_S.GE.2) THEN
         WHERE(ZA_CS(:)<-1.E-20 .AND. ZCST1D(:)>LIMAP%XCTMIN(5))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_CS(:)+ZCST1D(:))/ZA_CS(:))
         END WHERE
      ENDIF
      WHERE(ZA_RG(:)<-1.E-20 .AND. ZRGT1D(:)>LIMAP%XRTMIN(6))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RG(:)+ZRGT1D(:))/ZA_RG(:))
      END WHERE
      IF (LIMAP%NMOM_G.GE.2) THEN
         WHERE(ZA_CG(:)<-1.E-20 .AND. ZCGT1D(:)>LIMAP%XCTMIN(6))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_CG(:)+ZCGT1D(:))/ZA_CG(:))
         END WHERE
      ENDIF
      WHERE(ZA_RH(:)<-1.E-20 .AND. ZRHT1D(:)>LIMAP%XRTMIN(7))
         ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RH(:)+ZRHT1D(:))/ZA_RH(:))
      END WHERE
      IF (LIMAP%NMOM_H.GE.2) THEN
         WHERE(ZA_CH(:)<-1.E-20 .AND. ZCHT1D(:)>LIMAP%XCTMIN(7))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_CH(:)+ZCHT1D(:))/ZA_CH(:))
         END WHERE
      ENDIF

      ! We stop when the end of the timestep is reached
      WHERE(PTSTEP-ZTIME1D(:)-ZMAXTIME(:)<=0.)
         GLCOMPUTE1D(:)=.FALSE.
      ENDWHERE

      ! We must recompute tendencies when the end of the sub-timestep is reached
      IF(LIMAP%XTSTEP_TS/=0.) THEN
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ZTIME1D(:)+ZMAXTIME(:)>ZTIME_LASTCALL1D(:)+ZTSTEP)
            ZMAXTIME(:)=ZTIME_LASTCALL1D(:)-ZTIME1D(:)+ZTSTEP
            GLCOMPUTE1D(:)=.FALSE.
         ENDWHERE
      ENDIF

      ! We must recompute tendencies when the maximum allowed change is reached
      ! When a species is missing, only the external tendencies can be active and we do not want to recompute
      ! the microphysical tendencies when external tendencies are negative (results won't change because species was already missing)
      IF(LIMAP%XMRSTEP/=0.) THEN
         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RV(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RV(:))*LIMAP%XMRSTEP+Z0RVT1D(:)-ZRVT1D(:)-ZB_RV(:))/ZA_RV(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRVT1D(:)>LIMAP%XRTMIN(1) .OR. ZA_RV(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            GLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RC(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RC(:))*LIMAP%XMRSTEP+Z0RCT1D(:)-ZRCT1D(:)-ZB_RC(:))/ZA_RC(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRCT1D(:)>LIMAP%XRTMIN(2) .OR. ZA_RC(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            GLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RR(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RR(:))*LIMAP%XMRSTEP+Z0RRT1D(:)-ZRRT1D(:)-ZB_RR(:))/ZA_RR(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRRT1D(:)>LIMAP%XRTMIN(3) .OR. ZA_RR(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            GLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RI(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RI(:))*LIMAP%XMRSTEP+Z0RIT1D(:)-ZRIT1D(:)-ZB_RI(:))/ZA_RI(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRIT1D(:)>LIMAP%XRTMIN(4) .OR. ZA_RI(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            GLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RS(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RS(:))*LIMAP%XMRSTEP+Z0RST1D(:)-ZRST1D(:)-ZB_RS(:))/ZA_RS(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRST1D(:)>LIMAP%XRTMIN(5) .OR. ZA_RS(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            GLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RG(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RG(:))*LIMAP%XMRSTEP+Z0RGT1D(:)-ZRGT1D(:)-ZB_RG(:))/ZA_RG(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRGT1D(:)>LIMAP%XRTMIN(6) .OR. ZA_RG(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            GLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         ZTIME_THRESHOLD(:)=-1.
         WHERE(IITER1D(:)<INB_ITER_MAX .AND. ABS(ZA_RH(:))>1.E-20)
            ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RH(:))*LIMAP%XMRSTEP+Z0RHT1D(:)-ZRHT1D(:)-ZB_RH(:))/ZA_RH(:)
         ENDWHERE
         WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
              &(ZRHT1D(:)>LIMAP%XRTMIN(7) .OR. ZA_RH(:)>0.))
            ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
            GLCOMPUTE1D(:)=.FALSE.
         ENDWHERE

         WHERE(IITER1D(:)<INB_ITER_MAX .AND. MAX(ABS(ZB_RV(:)),  &
              ABS(ZB_RC(:)), ABS(ZB_RR(:)), ABS(ZB_RI(:)), &
              ABS(ZB_RS(:)), ABS(ZB_RG(:)), ABS(ZB_RH(:)))>LIMAP%XMRSTEP)
            ZMAXTIME(:)=0.
            GLCOMPUTE1D(:)=.FALSE.
         ENDWHERE
      ENDIF
      !
      !***       4.3 New values of variables for next iteration
      !
      ZTHT1D = ZTHT1D + ZA_TH(:) * ZMAXTIME(:) + ZB_TH(:)
      ZRVT1D = ZRVT1D + ZA_RV(:) * ZMAXTIME(:) + ZB_RV(:)
      ZRCT1D = ZRCT1D + ZA_RC(:) * ZMAXTIME(:) + ZB_RC(:)
      IF (LIMAP%NMOM_C.GE.2) ZCCT1D = ZCCT1D + ZA_CC(:) * ZMAXTIME(:) + ZB_CC(:)
      ZRRT1D = ZRRT1D + ZA_RR(:) * ZMAXTIME(:) + ZB_RR(:)
      IF (LIMAP%NMOM_R.GE.2) ZCRT1D = ZCRT1D + ZA_CR(:) * ZMAXTIME(:) + ZB_CR(:)
      ZRIT1D = ZRIT1D + ZA_RI(:) * ZMAXTIME(:) + ZB_RI(:)
      IF (LIMAP%NMOM_I.GE.2) ZCIT1D = ZCIT1D + ZA_CI(:) * ZMAXTIME(:) + ZB_CI(:)
      ZRST1D = ZRST1D + ZA_RS(:) * ZMAXTIME(:) + ZB_RS(:)
      IF (LIMAP%NMOM_S.GE.2) ZCST1D = ZCST1D + ZA_CS(:) * ZMAXTIME(:) + ZB_CS(:)
      ZRGT1D = ZRGT1D + ZA_RG(:) * ZMAXTIME(:) + ZB_RG(:)
      IF (LIMAP%NMOM_G.GE.2) ZCGT1D = ZCGT1D + ZA_CG(:) * ZMAXTIME(:) + ZB_CG(:)
      ZRHT1D = ZRHT1D + ZA_RH(:) * ZMAXTIME(:) + ZB_RH(:)
      IF (LIMAP%NMOM_H.GE.2) ZCHT1D = ZCHT1D + ZA_CH(:) * ZMAXTIME(:) + ZB_CH(:)
      IF (LIMAP%LCRYSTAL_SHAPE) THEN
        DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
          ZCIT1D_SHAPE(:,ISH) = ZCIT1D_SHAPE(:,ISH) + &
                                ZA_CI_SHAPE(:,ISH) * ZMAXTIME(:) + ZB_CI_SHAPE(:,ISH)
        END DO
        ZCIT1D(:) = SUM(ZCIT1D_SHAPE,DIM=2)
      END IF
      !
      DO II=1,LIMAP%NMOD_IFN
         ZIFNN1D(:,II) = ZIFNN1D(:,II) + ZB_IFNN(:,II)
      END DO
      !
      !***       4.5 
      !
      WHERE (ZRCT1D .LE. LIMAP%XRTMIN(2))
         Z_RV_CORR2(:) = ZRCT1D(:)
         Z_RC_CORR2(:) = -ZRCT1D(:)
         Z_CC_CORR2(:) = -ZCCT1D(:)

         ZRVT1D = ZRVT1D + ZRCT1D
         ZRCT1D = 0.
         ZCCT1D = 0.
      END WHERE
      WHERE (ZRRT1D .LE. LIMAP%XRTMIN(3))
         Z_RV_CORR2(:) = Z_RV_CORR2(:) + ZRRT1D(:)
         Z_RR_CORR2(:) = -ZRRT1D(:)
         Z_CR_CORR2(:) = -ZCRT1D(:)

         ZRVT1D = ZRVT1D + ZRRT1D
         ZRRT1D = 0.
         ZCRT1D = 0.
      END WHERE
! on traite les formes avant pour eviter que zrit1d ne soit modifie (where)
      IF (LIMAP%LCRYSTAL_SHAPE) THEN
        DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
          WHERE (ZRIT1D .LE. LIMAP%XRTMIN(4))
            Z_SHCI_CORR2(:,ISH) = -ZCIT1D_SHAPE(:,ISH)
            ZCIT1D_SHAPE(:,ISH) = 0.
          END WHERE
        END DO
        ZCIT1D(:) = SUM(ZCIT1D_SHAPE,DIM=2) !++cb-- 08/03/24
      END IF
      WHERE (ZRIT1D .LE. LIMAP%XRTMIN(4))
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
         ZTHT(I1(II),I3(II))      = ZTHT1D(II)
         ZRVT(I1(II),I3(II))      = ZRVT1D(II)
         ZRCT(I1(II),I3(II))      = ZRCT1D(II)
         ZRRT(I1(II),I3(II))      = ZRRT1D(II)
         ZRIT(I1(II),I3(II))      = ZRIT1D(II)
         ZRST(I1(II),I3(II))      = ZRST1D(II)
         ZRGT(I1(II),I3(II))      = ZRGT1D(II)
         ZRHT(I1(II),I3(II))      = ZRHT1D(II)
         IF (LIMAP%NMOM_C.GE.2) ZCCT(I1(II),I3(II))      = ZCCT1D(II)
         IF (LIMAP%NMOM_R.GE.2) ZCRT(I1(II),I3(II))      = ZCRT1D(II)
         IF (LIMAP%NMOM_I.GE.1) ZCIT(I1(II),I3(II))      = ZCIT1D(II)
         IF (LIMAP%NMOM_S.GE.2) ZCST(I1(II),I3(II))      = ZCST1D(II)
         IF (LIMAP%NMOM_G.GE.2) ZCGT(I1(II),I3(II))      = ZCGT1D(II)
         IF (LIMAP%NMOM_H.GE.2) ZCHT(I1(II),I3(II))      = ZCHT1D(II)
         IF (LIMAP%LCRYSTAL_SHAPE) THEN
           DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
             ZCIT_SHAPE(I1(II),I3(II),ISH) = ZCIT1D_SHAPE(II,ISH)
           END DO
         END IF
         ZIFNNT(I1(II),I3(II),:)  = ZIFNN1D(II,:)
         PEVAP3D(I1(II),I3(II))   = ZEVAP1D(II)
         ZTIME(I1(II),I3(II))     = ZTIME1D(II)
         GLCOMPUTE(I1(II),I3(II)) = GLCOMPUTE1D(II)
         IITER(I1(II),I3(II))     = IITER1D(II)
      END DO
      !
!!$      IF (LIMAP%NMOM_C.GE.2 .AND. LIMAP%NMOM_R.GE.2) THEN
!!$         CALL LIMA_DROPS_TO_DROPLETS_CONV(LIMAP, D, PRHODREF, ZRCT, ZRRT, ZCCT, ZCRT, &
!!$              Z_RR_CVRC, Z_CR_CVRC    )
!!$         ZRCT(:,:) = ZRCT(:,:) - Z_RR_CVRC(:,:)
!!$         ZRRT(:,:) = ZRRT(:,:) + Z_RR_CVRC(:,:)
!!$         ZCCT(:,:) = ZCCT(:,:) - Z_CR_CVRC(:,:)
!!$         ZCRT(:,:) = ZCRT(:,:) + Z_CR_CVRC(:,:)
!!$      END IF
      !
      !***       4.4 Unpacking for budgets
      !
      IF(BUCONF%LBU_ENABLE .OR. OELEC) THEN
        ZTOT_RR_CVRC(D%NIJB:D%NIJE,:) = ZTOT_RR_CVRC(D%NIJB:D%NIJE,:) + Z_RR_CVRC(D%NIJB:D%NIJE,:)
        ZTOT_CR_CVRC(D%NIJB:D%NIJE,:) = ZTOT_CR_CVRC(D%NIJB:D%NIJE,:) + Z_CR_CVRC(D%NIJB:D%NIJE,:)

         DO II=1,IPACK
            ! Instantaneous processes
            ZTOT_CR_BRKU(I1(II),I3(II)) =   ZTOT_CR_BRKU(I1(II),I3(II))   + Z_CR_BRKU(II)
            ZTOT_TH_HONR(I1(II),I3(II)) =   ZTOT_TH_HONR(I1(II),I3(II))   + Z_TH_HONR(II)
            ZTOT_RR_HONR(I1(II),I3(II)) =   ZTOT_RR_HONR(I1(II),I3(II))   + Z_RR_HONR(II)
            ZTOT_CR_HONR(I1(II),I3(II)) =   ZTOT_CR_HONR(I1(II),I3(II))   + Z_CR_HONR(II)
            ZTOT_TH_IMLT(I1(II),I3(II)) =   ZTOT_TH_IMLT(I1(II),I3(II))   + Z_TH_IMLT(II)
            ZTOT_RC_IMLT(I1(II),I3(II)) =   ZTOT_RC_IMLT(I1(II),I3(II))   + Z_RC_IMLT(II)
            ZTOT_CC_IMLT(I1(II),I3(II)) =   ZTOT_CC_IMLT(I1(II),I3(II))   + Z_CC_IMLT(II)
            DO IN = 1, LIMAP%NMOD_IFN
              ZTOT_IFNN_IMLT(I1(II),I3(II),IN) = ZTOT_IFNN_IMLT(I1(II),I3(II),IN) + ZB_IFNN(II,IN)
            END DO

            ! Tendencies
            ZTOT_TH_HONC(I1(II),I3(II)) =   ZTOT_TH_HONC(I1(II),I3(II))   + Z_TH_HONC(II)  * ZMAXTIME(II)
            ZTOT_RC_HONC(I1(II),I3(II)) =   ZTOT_RC_HONC(I1(II),I3(II))   + Z_RC_HONC(II)  * ZMAXTIME(II)
            ZTOT_CC_HONC(I1(II),I3(II)) =   ZTOT_CC_HONC(I1(II),I3(II))   + Z_CC_HONC(II)  * ZMAXTIME(II)
            ZTOT_CC_SELF(I1(II),I3(II)) =   ZTOT_CC_SELF(I1(II),I3(II))   + Z_CC_SELF(II)  * ZMAXTIME(II)
            ZTOT_RC_AUTO(I1(II),I3(II)) =   ZTOT_RC_AUTO(I1(II),I3(II))   + Z_RC_AUTO(II)  * ZMAXTIME(II)
            ZTOT_CC_AUTO(I1(II),I3(II)) =   ZTOT_CC_AUTO(I1(II),I3(II))   + Z_CC_AUTO(II)  * ZMAXTIME(II)
            ZTOT_CR_AUTO(I1(II),I3(II)) =   ZTOT_CR_AUTO(I1(II),I3(II))   + Z_CR_AUTO(II)  * ZMAXTIME(II)
            ZTOT_RC_ACCR(I1(II),I3(II)) =   ZTOT_RC_ACCR(I1(II),I3(II))   + Z_RC_ACCR(II)  * ZMAXTIME(II)
            ZTOT_CC_ACCR(I1(II),I3(II)) =   ZTOT_CC_ACCR(I1(II),I3(II))   + Z_CC_ACCR(II)  * ZMAXTIME(II)
            ZTOT_CR_SCBU(I1(II),I3(II)) =   ZTOT_CR_SCBU(I1(II),I3(II))   + Z_CR_SCBU(II)  * ZMAXTIME(II)
            ZTOT_TH_EVAP(I1(II),I3(II)) =   ZTOT_TH_EVAP(I1(II),I3(II))   + Z_TH_EVAP(II)  * ZMAXTIME(II)
!!$            ZTOT_RC_EVAP(I1(II),I3(II)) =   ZTOT_RC_EVAP(I1(II),I3(II))   + Z_RC_EVAP(II)  * ZMAXTIME(II)
!!$            ZTOT_CC_EVAP(I1(II),I3(II)) =   ZTOT_CC_EVAP(I1(II),I3(II))   + Z_CC_EVAP(II)  * ZMAXTIME(II)
            ZTOT_RR_EVAP(I1(II),I3(II)) =   ZTOT_RR_EVAP(I1(II),I3(II))   + Z_RR_EVAP(II)  * ZMAXTIME(II)
            ZTOT_CR_EVAP(I1(II),I3(II)) =   ZTOT_CR_EVAP(I1(II),I3(II))   + Z_CR_EVAP(II)  * ZMAXTIME(II)
            ZTOT_RI_CNVI(I1(II),I3(II)) =   ZTOT_RI_CNVI(I1(II),I3(II))   + Z_RI_CNVI(II)  * ZMAXTIME(II)
            ZTOT_CI_CNVI(I1(II),I3(II)) =   ZTOT_CI_CNVI(I1(II),I3(II))   + Z_CI_CNVI(II)  * ZMAXTIME(II)
            ZTOT_TH_DEPS(I1(II),I3(II)) =   ZTOT_TH_DEPS(I1(II),I3(II))   + Z_TH_DEPS(II)  * ZMAXTIME(II)
            ZTOT_RS_DEPS(I1(II),I3(II)) =   ZTOT_RS_DEPS(I1(II),I3(II))   + Z_RS_DEPS(II)  * ZMAXTIME(II)
            ZTOT_TH_DEPI(I1(II),I3(II)) =   ZTOT_TH_DEPI(I1(II),I3(II))   + Z_TH_DEPI(II)  * ZMAXTIME(II)
            ZTOT_RI_DEPI(I1(II),I3(II)) =   ZTOT_RI_DEPI(I1(II),I3(II))   + Z_RI_DEPI(II)  * ZMAXTIME(II)
            ZTOT_RI_CNVS(I1(II),I3(II)) =   ZTOT_RI_CNVS(I1(II),I3(II))   + Z_RI_CNVS(II)  * ZMAXTIME(II)
            ZTOT_CI_CNVS(I1(II),I3(II)) =   ZTOT_CI_CNVS(I1(II),I3(II))   + Z_CI_CNVS(II)  * ZMAXTIME(II)
            ZTOT_CS_SSC(I1(II),I3(II))  =   ZTOT_CS_SSC(I1(II),I3(II))    + Z_CS_SSC(II)   * ZMAXTIME(II)
            ZTOT_CI_ISC(I1(II),I3(II))  =   ZTOT_CI_ISC(I1(II),I3(II))    + Z_CI_ISC(II)   * ZMAXTIME(II)
            ZTOT_RI_AGGS(I1(II),I3(II)) =   ZTOT_RI_AGGS(I1(II),I3(II))   + Z_RI_AGGS(II)  * ZMAXTIME(II)
            ZTOT_CI_AGGS(I1(II),I3(II)) =   ZTOT_CI_AGGS(I1(II),I3(II))   + Z_CI_AGGS(II)  * ZMAXTIME(II)
            ZTOT_TH_DEPG(I1(II),I3(II)) =   ZTOT_TH_DEPG(I1(II),I3(II))   + Z_TH_DEPG(II)  * ZMAXTIME(II)
            ZTOT_RG_DEPG(I1(II),I3(II)) =   ZTOT_RG_DEPG(I1(II),I3(II))   + Z_RG_DEPG(II)  * ZMAXTIME(II)
            ZTOT_TH_BERFI(I1(II),I3(II))=   ZTOT_TH_BERFI(I1(II),I3(II))  + Z_TH_BERFI(II) * ZMAXTIME(II)
            ZTOT_RC_BERFI(I1(II),I3(II))=   ZTOT_RC_BERFI(I1(II),I3(II))  + Z_RC_BERFI(II) * ZMAXTIME(II)
            ZTOT_TH_RIM(I1(II),I3(II))  =   ZTOT_TH_RIM(I1(II),I3(II))    + Z_TH_RIM(II)   * ZMAXTIME(II)
            ZTOT_CC_RIM(I1(II),I3(II))  =   ZTOT_CC_RIM(I1(II),I3(II))    + Z_CC_RIM(II)   * ZMAXTIME(II)
            ZTOT_CS_RIM(I1(II),I3(II))  =   ZTOT_CS_RIM(I1(II),I3(II))    + Z_CS_RIM(II)   * ZMAXTIME(II)
            ZTOT_RC_RIMSS(I1(II),I3(II))=   ZTOT_RC_RIMSS(I1(II),I3(II))  + Z_RC_RIMSS(II) * ZMAXTIME(II)
            ZTOT_RC_RIMSG(I1(II),I3(II))=   ZTOT_RC_RIMSG(I1(II),I3(II))  + Z_RC_RIMSG(II) * ZMAXTIME(II)
            ZTOT_RS_RIMCG(I1(II),I3(II))=   ZTOT_RS_RIMCG(I1(II),I3(II))  + Z_RS_RIMCG(II) * ZMAXTIME(II)
            ZTOT_RI_HMS(I1(II),I3(II))  =   ZTOT_RI_HMS(I1(II),I3(II))    + Z_RI_HMS(II)   * ZMAXTIME(II)
            ZTOT_CI_HMS(I1(II),I3(II))  =   ZTOT_CI_HMS(I1(II),I3(II))    + Z_CI_HMS(II)   * ZMAXTIME(II)
            ZTOT_RS_HMS(I1(II),I3(II))  =   ZTOT_RS_HMS(I1(II),I3(II))    + Z_RS_HMS(II)   * ZMAXTIME(II)
            ZTOT_TH_ACC(I1(II),I3(II))  =   ZTOT_TH_ACC(I1(II),I3(II))    + Z_TH_ACC(II)   * ZMAXTIME(II)
            ZTOT_CR_ACC(I1(II),I3(II))  =   ZTOT_CR_ACC(I1(II),I3(II))    + Z_CR_ACC(II)   * ZMAXTIME(II)
            ZTOT_CS_ACC(I1(II),I3(II))  =   ZTOT_CS_ACC(I1(II),I3(II))    + Z_CS_ACC(II)   * ZMAXTIME(II)
            ZTOT_RR_ACCSS(I1(II),I3(II))=   ZTOT_RR_ACCSS(I1(II),I3(II))  + Z_RR_ACCSS(II) * ZMAXTIME(II)
            ZTOT_RR_ACCSG(I1(II),I3(II))=   ZTOT_RR_ACCSG(I1(II),I3(II))  + Z_RR_ACCSG(II) * ZMAXTIME(II)
            ZTOT_RS_ACCRG(I1(II),I3(II))=   ZTOT_RS_ACCRG(I1(II),I3(II))  + Z_RS_ACCRG(II) * ZMAXTIME(II)
            ZTOT_CS_CMEL(I1(II),I3(II)) =   ZTOT_CS_CMEL(I1(II),I3(II))   + Z_CS_CMEL(II)  * ZMAXTIME(II)
            ZTOT_RS_CMEL(I1(II),I3(II)) =   ZTOT_RS_CMEL(I1(II),I3(II))   + Z_RS_CMEL(II)  * ZMAXTIME(II)
            ZTOT_TH_CFRZ(I1(II),I3(II)) =   ZTOT_TH_CFRZ(I1(II),I3(II))   + Z_TH_CFRZ(II)  * ZMAXTIME(II)
            ZTOT_RR_CFRZ(I1(II),I3(II)) =   ZTOT_RR_CFRZ(I1(II),I3(II))   + Z_RR_CFRZ(II)  * ZMAXTIME(II)
            ZTOT_CR_CFRZ(I1(II),I3(II)) =   ZTOT_CR_CFRZ(I1(II),I3(II))   + Z_CR_CFRZ(II)  * ZMAXTIME(II)
            ZTOT_RI_CFRZ(I1(II),I3(II)) =   ZTOT_RI_CFRZ(I1(II),I3(II))   + Z_RI_CFRZ(II)  * ZMAXTIME(II)
            ZTOT_CI_CFRZ(I1(II),I3(II)) =   ZTOT_CI_CFRZ(I1(II),I3(II))   + Z_CI_CFRZ(II)  * ZMAXTIME(II)
            ZTOT_RI_CIBU(I1(II),I3(II)) =   ZTOT_RI_CIBU(I1(II),I3(II))   + Z_RI_CIBU(II)  * ZMAXTIME(II)
            ZTOT_CI_CIBU(I1(II),I3(II)) =   ZTOT_CI_CIBU(I1(II),I3(II))   + Z_CI_CIBU(II)  * ZMAXTIME(II)
            ZTOT_RI_RDSF(I1(II),I3(II)) =   ZTOT_RI_RDSF(I1(II),I3(II))   + Z_RI_RDSF(II)  * ZMAXTIME(II)
            ZTOT_CI_RDSF(I1(II),I3(II)) =   ZTOT_CI_RDSF(I1(II),I3(II))   + Z_CI_RDSF(II)  * ZMAXTIME(II)
            ZTOT_TH_WETG(I1(II),I3(II)) =   ZTOT_TH_WETG(I1(II),I3(II))   + Z_TH_WETG(II)  * ZMAXTIME(II)
            ZTOT_RC_WETG(I1(II),I3(II)) =   ZTOT_RC_WETG(I1(II),I3(II))   + Z_RC_WETG(II)  * ZMAXTIME(II)
            ZTOT_CC_WETG(I1(II),I3(II)) =   ZTOT_CC_WETG(I1(II),I3(II))   + Z_CC_WETG(II)  * ZMAXTIME(II)
            ZTOT_RR_WETG(I1(II),I3(II)) =   ZTOT_RR_WETG(I1(II),I3(II))   + Z_RR_WETG(II)  * ZMAXTIME(II)
            ZTOT_CR_WETG(I1(II),I3(II)) =   ZTOT_CR_WETG(I1(II),I3(II))   + Z_CR_WETG(II)  * ZMAXTIME(II)
            ZTOT_RI_WETG(I1(II),I3(II)) =   ZTOT_RI_WETG(I1(II),I3(II))   + Z_RI_WETG(II)  * ZMAXTIME(II)
            ZTOT_CI_WETG(I1(II),I3(II)) =   ZTOT_CI_WETG(I1(II),I3(II))   + Z_CI_WETG(II)  * ZMAXTIME(II)
            ZTOT_RS_WETG(I1(II),I3(II)) =   ZTOT_RS_WETG(I1(II),I3(II))   + Z_RS_WETG(II)  * ZMAXTIME(II)
            ZTOT_CS_WETG(I1(II),I3(II)) =   ZTOT_CS_WETG(I1(II),I3(II))   + Z_CS_WETG(II)  * ZMAXTIME(II)
            ZTOT_RG_WETG(I1(II),I3(II)) =   ZTOT_RG_WETG(I1(II),I3(II))   + Z_RG_WETG(II)  * ZMAXTIME(II)
            ZTOT_CG_WETG(I1(II),I3(II)) =   ZTOT_CG_WETG(I1(II),I3(II))   + Z_CG_WETG(II)  * ZMAXTIME(II)
            ZTOT_RH_WETG(I1(II),I3(II)) =   ZTOT_RH_WETG(I1(II),I3(II))   + Z_RH_WETG(II)  * ZMAXTIME(II)
            ZTOT_TH_DRYG(I1(II),I3(II)) =   ZTOT_TH_DRYG(I1(II),I3(II))   + Z_TH_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RC_DRYG(I1(II),I3(II)) =   ZTOT_RC_DRYG(I1(II),I3(II))   + Z_RC_DRYG(II)  * ZMAXTIME(II)
            ZTOT_CC_DRYG(I1(II),I3(II)) =   ZTOT_CC_DRYG(I1(II),I3(II))   + Z_CC_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RR_DRYG(I1(II),I3(II)) =   ZTOT_RR_DRYG(I1(II),I3(II))   + Z_RR_DRYG(II)  * ZMAXTIME(II)
            ZTOT_CR_DRYG(I1(II),I3(II)) =   ZTOT_CR_DRYG(I1(II),I3(II))   + Z_CR_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RI_DRYG(I1(II),I3(II)) =   ZTOT_RI_DRYG(I1(II),I3(II))   + Z_RI_DRYG(II)  * ZMAXTIME(II)
            ZTOT_CI_DRYG(I1(II),I3(II)) =   ZTOT_CI_DRYG(I1(II),I3(II))   + Z_CI_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RS_DRYG(I1(II),I3(II)) =   ZTOT_RS_DRYG(I1(II),I3(II))   + Z_RS_DRYG(II)  * ZMAXTIME(II)
            ZTOT_CS_DRYG(I1(II),I3(II)) =   ZTOT_CS_DRYG(I1(II),I3(II))   + Z_CS_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RG_DRYG(I1(II),I3(II)) =   ZTOT_RG_DRYG(I1(II),I3(II))   + Z_RG_DRYG(II)  * ZMAXTIME(II)
            ZTOT_RI_HMG(I1(II),I3(II))  =   ZTOT_RI_HMG(I1(II),I3(II))    + Z_RI_HMG(II)   * ZMAXTIME(II)
            ZTOT_CI_HMG(I1(II),I3(II))  =   ZTOT_CI_HMG(I1(II),I3(II))    + Z_CI_HMG(II)   * ZMAXTIME(II)
            ZTOT_RG_HMG(I1(II),I3(II))  =   ZTOT_RG_HMG(I1(II),I3(II))    + Z_RG_HMG(II)   * ZMAXTIME(II)
            ZTOT_TH_GMLT(I1(II),I3(II)) =   ZTOT_TH_GMLT(I1(II),I3(II))   + Z_TH_GMLT(II)  * ZMAXTIME(II)
            ZTOT_RR_GMLT(I1(II),I3(II)) =   ZTOT_RR_GMLT(I1(II),I3(II))   + Z_RR_GMLT(II)  * ZMAXTIME(II)
            ZTOT_CR_GMLT(I1(II),I3(II)) =   ZTOT_CR_GMLT(I1(II),I3(II))   + Z_CR_GMLT(II)  * ZMAXTIME(II)
            ZTOT_CG_GMLT(I1(II),I3(II)) =   ZTOT_CG_GMLT(I1(II),I3(II))   + Z_CG_GMLT(II)  * ZMAXTIME(II)
            ZTOT_TH_DEPH(I1(II),I3(II)) =   ZTOT_TH_DEPH(I1(II),I3(II))   + Z_TH_DEPH(II)  * ZMAXTIME(II)
            ZTOT_RH_DEPH(I1(II),I3(II)) =   ZTOT_RH_DEPH(I1(II),I3(II))   + Z_RH_DEPH(II)  * ZMAXTIME(II)
            ZTOT_TH_WETH(I1(II),I3(II)) =   ZTOT_TH_WETH(I1(II),I3(II))   + Z_TH_WETH(II)  * ZMAXTIME(II)
            ZTOT_RC_WETH(I1(II),I3(II)) =   ZTOT_RC_WETH(I1(II),I3(II))   + Z_RC_WETH(II)  * ZMAXTIME(II)
            ZTOT_CC_WETH(I1(II),I3(II)) =   ZTOT_CC_WETH(I1(II),I3(II))   + Z_CC_WETH(II)  * ZMAXTIME(II)
            ZTOT_RR_WETH(I1(II),I3(II)) =   ZTOT_RR_WETH(I1(II),I3(II))   + Z_RR_WETH(II)  * ZMAXTIME(II)
            ZTOT_CR_WETH(I1(II),I3(II)) =   ZTOT_CR_WETH(I1(II),I3(II))   + Z_CR_WETH(II)  * ZMAXTIME(II)
            ZTOT_RI_WETH(I1(II),I3(II)) =   ZTOT_RI_WETH(I1(II),I3(II))   + Z_RI_WETH(II)  * ZMAXTIME(II)
            ZTOT_CI_WETH(I1(II),I3(II)) =   ZTOT_CI_WETH(I1(II),I3(II))   + Z_CI_WETH(II)  * ZMAXTIME(II)
            ZTOT_RS_WETH(I1(II),I3(II)) =   ZTOT_RS_WETH(I1(II),I3(II))   + Z_RS_WETH(II)  * ZMAXTIME(II)
            ZTOT_CS_WETH(I1(II),I3(II)) =   ZTOT_CS_WETH(I1(II),I3(II))   + Z_CS_WETH(II)  * ZMAXTIME(II)
            ZTOT_RG_WETH(I1(II),I3(II)) =   ZTOT_RG_WETH(I1(II),I3(II))   + Z_RG_WETH(II)  * ZMAXTIME(II)
            ZTOT_CG_WETH(I1(II),I3(II)) =   ZTOT_CG_WETH(I1(II),I3(II))   + Z_CG_WETH(II)  * ZMAXTIME(II)
            ZTOT_RH_WETH(I1(II),I3(II)) =   ZTOT_RH_WETH(I1(II),I3(II))   + Z_RH_WETH(II)  * ZMAXTIME(II)
            ZTOT_RG_COHG(I1(II),I3(II)) =   ZTOT_RG_COHG(I1(II),I3(II))   + Z_RG_COHG(II)  * ZMAXTIME(II)
            ZTOT_CG_COHG(I1(II),I3(II)) =   ZTOT_CG_COHG(I1(II),I3(II))   + Z_CG_COHG(II)  * ZMAXTIME(II)
            ZTOT_TH_HMLT(I1(II),I3(II)) =   ZTOT_RR_HMLT(I1(II),I3(II))   + Z_RR_HMLT(II)  * ZMAXTIME(II)
            ZTOT_RR_HMLT(I1(II),I3(II)) =   ZTOT_RR_HMLT(I1(II),I3(II))   + Z_RR_HMLT(II)  * ZMAXTIME(II)
            ZTOT_CR_HMLT(I1(II),I3(II)) =   ZTOT_CR_HMLT(I1(II),I3(II))   + Z_CR_HMLT(II)  * ZMAXTIME(II)
            ZTOT_CH_HMLT(I1(II),I3(II)) =   ZTOT_CH_HMLT(I1(II),I3(II))   + Z_CH_HMLT(II)  * ZMAXTIME(II)

            ! Correction term
            ZTOT_RV_CORR2(I1(II),I3(II)) =   ZTOT_RV_CORR2(I1(II),I3(II)) + Z_RV_CORR2(II)
            ZTOT_RC_CORR2(I1(II),I3(II)) =   ZTOT_RC_CORR2(I1(II),I3(II)) + Z_RC_CORR2(II)
            ZTOT_RR_CORR2(I1(II),I3(II)) =   ZTOT_RR_CORR2(I1(II),I3(II)) + Z_RR_CORR2(II)
            ZTOT_RI_CORR2(I1(II),I3(II)) =   ZTOT_RI_CORR2(I1(II),I3(II)) + Z_RI_CORR2(II)
            ZTOT_CC_CORR2(I1(II),I3(II)) =   ZTOT_CC_CORR2(I1(II),I3(II)) + Z_CC_CORR2(II)
            ZTOT_CR_CORR2(I1(II),I3(II)) =   ZTOT_CR_CORR2(I1(II),I3(II)) + Z_CR_CORR2(II)
            ZTOT_CI_CORR2(I1(II),I3(II)) =   ZTOT_CI_CORR2(I1(II),I3(II)) + Z_CI_CORR2(II)

            ! Terms for ice crystal shapes
            IF (LIMAP%LCRYSTAL_SHAPE) THEN
              DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
                ! instantaneous process
                ZTOT_SHCI_IMLT(I1(II),I3(II),ISH) = ZTOT_SHCI_IMLT(I1(II),I3(II),ISH) + Z_SHCI_IMLT(II,ISH)
                ! tendencies
                ZTOT_SHCI_HONC(I1(II),I3(II),ISH) = ZTOT_SHCI_HONC(I1(II),I3(II),ISH) + &
                                                           Z_SHCI_HONC(II,ISH) * ZMAXTIME(II)
                ZTOT_SHCI_CNVI(I1(II),I3(II),ISH) = ZTOT_SHCI_CNVI(I1(II),I3(II),ISH) + &
                                                           Z_SHCI_CNVI(II,ISH) * ZMAXTIME(II)
                ZTOT_SHCI_HACH(I1(II),I3(II),ISH) = ZTOT_SHCI_HACH(I1(II),I3(II),ISH) + &
                                                           Z_SHCI_HACH(II,ISH) * ZMAXTIME(II)
                ZTOT_SHCI_CNVS(I1(II),I3(II),ISH) = ZTOT_SHCI_CNVS(I1(II),I3(II),ISH) + &
                                                           Z_SHCI_CNVS(II,ISH) * ZMAXTIME(II)
                ZTOT_SHCI_AGGS(I1(II),I3(II),ISH) = ZTOT_SHCI_AGGS(I1(II),I3(II),ISH) + &
                                                           Z_SHCI_AGGS(II,ISH) * ZMAXTIME(II)
                ZTOT_SHCI_HMS(I1(II),I3(II),ISH)  = ZTOT_SHCI_HMS(I1(II),I3(II),ISH)  + &
                                                           Z_SHCI_HMS(II,ISH)  * ZMAXTIME(II)
                ZTOT_SHCI_CFRZ(I1(II),I3(II),ISH) = ZTOT_SHCI_CFRZ(I1(II),I3(II),ISH) + &
                                                           Z_SHCI_CFRZ(II,ISH) * ZMAXTIME(II)
                ZTOT_SHCI_CIBU(I1(II),I3(II),ISH) = ZTOT_SHCI_CIBU(I1(II),I3(II),ISH) + &
                                                           Z_SHCI_CIBU(II,ISH) * ZMAXTIME(II)
                ZTOT_SHCI_RDSF(I1(II),I3(II),ISH) = ZTOT_SHCI_RDSF(I1(II),I3(II),ISH) + &
                                                           Z_SHCI_RDSF(II,ISH) * ZMAXTIME(II)
                ZTOT_SHCI_WETG(I1(II),I3(II),ISH) = ZTOT_SHCI_WETG(I1(II),I3(II),ISH) + &
                                                           Z_SHCI_WETG(II,ISH) * ZMAXTIME(II)
                ZTOT_SHCI_DRYG(I1(II),I3(II),ISH) = ZTOT_SHCI_DRYG(I1(II),I3(II),ISH) + &
                                                           Z_SHCI_DRYG(II,ISH) * ZMAXTIME(II)
                ZTOT_SHCI_HMG(I1(II),I3(II),ISH)  = ZTOT_SHCI_HMG(I1(II),I3(II),ISH)  + &
                                                           Z_SHCI_HMG(II,ISH)  * ZMAXTIME(II)
                ZTOT_SHCI_ISC(I1(II),I3(II),ISH)  = ZTOT_SHCI_ISC(I1(II),I3(II),ISH)  + &
                                                           Z_SHCI_ISC(II,ISH) * ZMAXTIME(II)
                ZTOT_SHCI_ISCS(I1(II),I3(II),ISH) = ZTOT_SHCI_ISCS(I1(II),I3(II),ISH) + &    
                                                           Z_SHCI_ISCS(II,ISH) * ZMAXTIME(II)
                ! Correction term
                ZTOT_SHCI_CORR2(I1(II),I3(II),ISH)  = ZTOT_SHCI_CORR2(I1(II),I3(II),ISH) + &
                                                             Z_SHCI_CORR2(II,ISH)
              END DO
              ZTOT_SHRI_ISCS(I1(II),I3(II)) = ZTOT_SHRI_ISCS(I1(II),I3(II)) + &    
                                                        Z_SHRI_ISCS(II) * ZMAXTIME(II) 
            END IF

         END DO
      ENDIF
      !
      ! Deallocating variables
      !
      DEALLOCATE(I1)
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
      DEALLOCATE(ZCST1D)
      DEALLOCATE(ZCGT1D)
      DEALLOCATE(ZCHT1D)
      DEALLOCATE(ZIFNN1D)
      DEALLOCATE(ZEVAP1D)
      DEALLOCATE(ZTIME1D)
      DEALLOCATE(GLCOMPUTE1D)
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
      DEALLOCATE(ZLATHAM_IAGGS)
      DEALLOCATE(ZHLC_HRC1D)
      DEALLOCATE(ZHLC_HCF1D)
      DEALLOCATE(ZHLI_HRI1D)
      DEALLOCATE(ZHLI_HCF1D)
      DEALLOCATE(ZHLC_LRC1D)
      DEALLOCATE(ZHLC_LCF1D)
      DEALLOCATE(ZHLI_LRI1D)
      DEALLOCATE(ZHLI_LCF1D)
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
      DEALLOCATE(ZA_CS)
      DEALLOCATE(ZA_CG)
      DEALLOCATE(ZA_CH)
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
      DEALLOCATE(ZB_CS)  
      DEALLOCATE(ZB_CG)  
      DEALLOCATE(ZB_CH)  
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
      DEALLOCATE(Z_CR_EVAP) 
      DEALLOCATE(Z_RI_CNVI)
      DEALLOCATE(Z_CI_CNVI)
      DEALLOCATE(Z_TH_DEPS)
      DEALLOCATE(Z_RS_DEPS)
      DEALLOCATE(Z_TH_DEPI)
      DEALLOCATE(Z_RI_DEPI)
      DEALLOCATE(Z_RI_CNVS)
      DEALLOCATE(Z_CI_CNVS)
      DEALLOCATE(Z_CS_SSC) 
      DEALLOCATE(Z_CI_ISC)
      DEALLOCATE(Z_RI_AGGS) 
      DEALLOCATE(Z_CI_AGGS) 
      DEALLOCATE(Z_TH_DEPG) 
      DEALLOCATE(Z_RG_DEPG) 
      DEALLOCATE(Z_TH_BERFI)
      DEALLOCATE(Z_RC_BERFI)
      DEALLOCATE(Z_TH_RIM) 
      DEALLOCATE(Z_CC_RIM)  
      DEALLOCATE(Z_CS_RIM) 
      DEALLOCATE(Z_RC_RIMSS)
      DEALLOCATE(Z_RC_RIMSG)
      DEALLOCATE(Z_RS_RIMCG)    
      DEALLOCATE(Z_RI_HMS) 
      DEALLOCATE(Z_CI_HMS) 
      DEALLOCATE(Z_RS_HMS)
      DEALLOCATE(Z_TH_ACC) 
      DEALLOCATE(Z_CR_ACC) 
      DEALLOCATE(Z_CS_ACC) 
      DEALLOCATE(Z_RR_ACCSS)
      DEALLOCATE(Z_RR_ACCSG)
      DEALLOCATE(Z_RS_ACCRG)
      DEALLOCATE(Z_CS_CMEL) 
      DEALLOCATE(Z_RS_CMEL) 
      DEALLOCATE(Z_TH_CFRZ)
      DEALLOCATE(Z_RR_CFRZ)
      DEALLOCATE(Z_CR_CFRZ)
      DEALLOCATE(Z_RI_CFRZ)
      DEALLOCATE(Z_CI_CFRZ)
      DEALLOCATE(Z_RI_CIBU) 
      DEALLOCATE(Z_CI_CIBU) 
      DEALLOCATE(Z_RI_RDSF) 
      DEALLOCATE(Z_CI_RDSF) 
      DEALLOCATE(Z_TH_WETG)
      DEALLOCATE(Z_RC_WETG)
      DEALLOCATE(Z_CC_WETG)
      DEALLOCATE(Z_RR_WETG) 
      DEALLOCATE(Z_CR_WETG) 
      DEALLOCATE(Z_RI_WETG)
      DEALLOCATE(Z_CI_WETG)
      DEALLOCATE(Z_RS_WETG)
      DEALLOCATE(Z_CS_WETG)
      DEALLOCATE(Z_RG_WETG)
      DEALLOCATE(Z_CG_WETG)
      DEALLOCATE(Z_RH_WETG) 
      DEALLOCATE(Z_TH_DRYG) 
      DEALLOCATE(Z_RC_DRYG) 
      DEALLOCATE(Z_CC_DRYG)
      DEALLOCATE(Z_RR_DRYG)
      DEALLOCATE(Z_CR_DRYG)
      DEALLOCATE(Z_RI_DRYG)
      DEALLOCATE(Z_CI_DRYG)
      DEALLOCATE(Z_RS_DRYG) 
      DEALLOCATE(Z_CS_DRYG) 
      DEALLOCATE(Z_RG_DRYG)
      DEALLOCATE(Z_RI_HMG) 
      DEALLOCATE(Z_CI_HMG) 
      DEALLOCATE(Z_RG_HMG) 
      DEALLOCATE(Z_TH_GMLT)
      DEALLOCATE(Z_RR_GMLT)
      DEALLOCATE(Z_CR_GMLT)
      DEALLOCATE(Z_CG_GMLT)
      DEALLOCATE(Z_TH_DEPH) 
      DEALLOCATE(Z_RH_DEPH) 
      DEALLOCATE(Z_TH_WETH)
      DEALLOCATE(Z_RC_WETH)
      DEALLOCATE(Z_CC_WETH)
      DEALLOCATE(Z_RR_WETH) 
      DEALLOCATE(Z_CR_WETH) 
      DEALLOCATE(Z_RI_WETH)
      DEALLOCATE(Z_CI_WETH)
      DEALLOCATE(Z_RS_WETH)
      DEALLOCATE(Z_CS_WETH)
      DEALLOCATE(Z_RG_WETH)
      DEALLOCATE(Z_CG_WETH)
      DEALLOCATE(Z_RH_WETH) 
      DEALLOCATE(Z_RG_COHG)
      DEALLOCATE(Z_CG_COHG)
      DEALLOCATE(Z_TH_HMLT) 
      DEALLOCATE(Z_RR_HMLT) 
      DEALLOCATE(Z_CR_HMLT) 
      DEALLOCATE(Z_CH_HMLT) 
      !
      DEALLOCATE(Z_RV_CORR2)
      DEALLOCATE(Z_RC_CORR2)
      DEALLOCATE(Z_RR_CORR2)
      DEALLOCATE(Z_RI_CORR2)
      DEALLOCATE(Z_CC_CORR2)
      DEALLOCATE(Z_CR_CORR2)
      DEALLOCATE(Z_CI_CORR2)
      !
      DEALLOCATE(ZCIT1D_SHAPE)
      DEALLOCATE(ZA_CI_SHAPE)
      DEALLOCATE(ZB_CI_SHAPE)
      DEALLOCATE(Z_SHCI_IMLT)
      DEALLOCATE(Z_SHCI_ISC)
      DEALLOCATE(Z_SHCI_ISCS)
      DEALLOCATE(Z_SHRI_ISCS)
      DEALLOCATE(Z_SHCI_HONC)
      DEALLOCATE(Z_SHCI_CNVI)
      DEALLOCATE(Z_SHCI_HACH)
      DEALLOCATE(Z_SHCI_CNVS)
      DEALLOCATE(Z_SHCI_AGGS)
      DEALLOCATE(Z_SHCI_HMS)
      DEALLOCATE(Z_SHCI_CFRZ)
      DEALLOCATE(Z_SHCI_CIBU)
      DEALLOCATE(Z_SHCI_RDSF)
      DEALLOCATE(Z_SHCI_WETG)
      DEALLOCATE(Z_SHCI_DRYG)
      DEALLOCATE(Z_SHCI_HMG)      
      DEALLOCATE(Z_SHCI_CORR2)
      !
   ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       7.     CLOUD ELECTRIFICATION
!               ---------------------
!
!*       7.1    Packing variables
!               -----------------
!
IF (OELEC) THEN
  ALLOCATE(GMASK_ELEC(D%NIJT,D%NKT))
  GMASK_ELEC(:,:) = .FALSE.
  GMASK_ELEC(D%NIJB:D%NIJE,:) = ZTOT_RI_HIND (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RR_HONR (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RC_IMLT (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RC_HONC (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RS_DEPS (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RI_AGGS (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RI_CNVS (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RG_DEPG (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RC_AUTO (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RC_ACCR (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RR_EVAP (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RC_RIMSS(D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RC_RIMSG(D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RS_RIMCG(D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RR_ACCSS(D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RR_ACCSG(D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RS_ACCRG(D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RS_CMEL (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RR_CFRZ (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RI_CFRZ (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RI_CIBU (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RI_RDSF (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RC_WETG (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RI_WETG (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RR_WETG (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RS_WETG (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RC_DRYG (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RI_DRYG (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RR_DRYG (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RS_DRYG (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RH_WETG (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RR_GMLT (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RC_BERFI(D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RV_HENU (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RC_HINC (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RV_HONH (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RR_CVRC (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RI_CNVI (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RI_DEPI (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RI_HMS  (D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RI_HMG  (D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RC_CORR2(D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                              ZTOT_RR_CORR2(D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RI_CORR2(D%NIJB:D%NIJE,:) .NE. 0.
  IF (LIMAP%NMOM_H .GE. 1) &
    GMASK_ELEC(D%NIJB:D%NIJE,:) = GMASK_ELEC(D%NIJB:D%NIJE,:)           .OR.                                  &
                        ZTOT_RC_WETH(D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RI_WETH(D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                        ZTOT_RS_WETH(D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RG_WETH(D%NIJB:D%NIJE,:) .NE. 0. .OR. &
                        ZTOT_RR_WETH(D%NIJB:D%NIJE,:) .NE. 0. .OR.                                  &
                        !ZTOT_RC_DRYH(:,:) .NE. 0. .OR. ZTOT_RI_DRYH(:,:) .NE. 0. .OR. &
                        !ZTOT_RS_DRYH(:,:) .NE. 0. .OR. ZTOT_RR_DRYH(:,:) .NE. 0. .OR. &
                        !ZTOT_RG_DRYH(:,:) .NE. 0. .OR.                                  &
                        ZTOT_RG_COHG(D%NIJB:D%NIJE,:) .NE. 0. .OR. ZTOT_RR_HMLT(D%NIJB:D%NIJE,:) .NE. 0.
  !
  IELEC = COUNT(GMASK_ELEC)
  !
!
!*       7.2    Cloud electrification:
!               ---------------------
!
! Attention, les signes des tendances ne sont pas traites de la meme facon dans ice3 et lima
! On se cale sur la facon de faire dans ice3 => on fait en sorte que les tendances soient positives
  IF (LIMAP%NMOM_H .GE. 1) THEN
    CALL ELEC_TENDENCIES(D, CST, ICED, ICEP, ELECD, ELECP,                                              &
                         KRR, IELEC, PTSTEP, GMASK_ELEC,                                                &
                         BUCONF, TBUDGETS, KBUDGETS,                                                    &
                         'LIMA', PTHVREFZIKB,                                                           &
                         PRHODREF,  PRHODJ, ZT, ZCIT_ELEC,                                              &
                         ZRVT_ELEC, ZRCT_ELEC, ZRRT_ELEC, ZRIT_ELEC, ZRST_ELEC, ZRGT_ELEC,              &
                         ZQPIT, ZQCT, ZQRT, ZQIT, ZQST, ZQGT, ZQNIT,                                    &
                         ZQPIS, ZQCS, ZQRS, ZQIS, ZQSS, ZQGS, ZQNIS,                                    &
                         ZTOT_RI_HIND*ZINV_TSTEP,  -ZTOT_RR_HONR*ZINV_TSTEP,   ZTOT_RC_IMLT*ZINV_TSTEP, &
                        -ZTOT_RC_HONC*ZINV_TSTEP,   ZTOT_RS_DEPS*ZINV_TSTEP,  -ZTOT_RI_AGGS*ZINV_TSTEP, &
                        -ZTOT_RI_CNVS*ZINV_TSTEP,   ZTOT_RG_DEPG*ZINV_TSTEP,  -ZTOT_RC_AUTO*ZINV_TSTEP, &
                        -ZTOT_RC_ACCR*ZINV_TSTEP,  -ZTOT_RR_EVAP*ZINV_TSTEP,                            &
                         ZTOT_RC_RIMSS*ZINV_TSTEP,  ZTOT_RC_RIMSG*ZINV_TSTEP,  ZTOT_RS_RIMCG*ZINV_TSTEP,&
                         ZTOT_RR_ACCSS*ZINV_TSTEP,  ZTOT_RR_ACCSG*ZINV_TSTEP,  ZTOT_RS_ACCRG*ZINV_TSTEP,&
                        -ZTOT_RS_CMEL*ZINV_TSTEP,  -ZTOT_RI_CFRZ*ZINV_TSTEP,  -ZTOT_RR_CFRZ*ZINV_TSTEP, &
                        -ZTOT_RC_WETG*ZINV_TSTEP,  -ZTOT_RI_WETG*ZINV_TSTEP,  -ZTOT_RR_WETG*ZINV_TSTEP, &
                        -ZTOT_RS_WETG*ZINV_TSTEP,                                                       &
                        -ZTOT_RC_DRYG*ZINV_TSTEP,  -ZTOT_RI_DRYG*ZINV_TSTEP,  -ZTOT_RR_DRYG*ZINV_TSTEP, &
                        -ZTOT_RS_DRYG*ZINV_TSTEP,                                                       &
                         ZTOT_RR_GMLT*ZINV_TSTEP,  -ZTOT_RC_BERFI*ZINV_TSTEP,                           &
! variables et processus optionnels propres a la grele : pas encore teste
                         PRWETGH=ZTOT_RH_WETG*ZINV_TSTEP,                                               &
                         PRCWETH=ZTOT_RC_WETH, PRIWETH=ZTOT_RI_WETH, PRSWETH=ZTOT_RS_WETH,              &
                         PRGWETH=ZTOT_RG_WETH, PRRWETH=ZTOT_RR_WETH,                                    &
!                         PRCDRYH=ZTOT_RC_DRYH, PRIDRYH=ZTOT_RI_DRYH, PRSDRYH=ZTOT_RS_DRYH,              &
!                         PRRDRYH=ZTOT_RR_DRYH, PRGDRYH=ZTOT_RG_DRYH,                                    &
                         PRDRYHG=ZTOT_RG_COHG, PRHMLTR=ZTOT_RR_HMLT,                                    &
                         PRHT=ZRHT, PQHT=ZQHT, PQHS=ZQHS, PCHT=ZCHT,                         &
! variables et processus optionnels propres a lima
                         PCCT=ZCCT_ELEC, PCRT=ZCRT_ELEC, PCST=ZCST_ELEC, PCGT=ZCGT_ELEC,                &
                         PRVHENC=ZTOT_RV_HENU*ZINV_TSTEP,    PRCHINC=-ZTOT_RC_HINC*ZINV_TSTEP,          &
                         PRVHONH=-ZTOT_RV_HONH*ZINV_TSTEP,   PRRCVRC=-ZTOT_RR_CVRC*ZINV_TSTEP,          &
                         PRICNVI=ZTOT_RI_CNVI*ZINV_TSTEP,    PRVDEPI=ZTOT_RI_DEPI*ZINV_TSTEP,           &
                         PRSHMSI=ZTOT_RI_HMS*ZINV_TSTEP,     PRGHMGI=ZTOT_RI_HMG*ZINV_TSTEP,            &
                         PRICIBU=ZTOT_RI_CIBU*ZINV_TSTEP,    PRIRDSF=ZTOT_RI_RDSF*ZINV_TSTEP,           &
                         PRCCORR2=-ZTOT_RC_CORR2*ZINV_TSTEP, PRRCORR2=-ZTOT_RR_CORR2*ZINV_TSTEP,        &
                         PRICORR2=-ZTOT_RI_CORR2*ZINV_TSTEP)
  ELSE
    CALL ELEC_TENDENCIES(D, CST, ICED, ICEP, ELECD, ELECP,                                              &
                         KRR, IELEC, PTSTEP, GMASK_ELEC,                                                &
                         BUCONF, TBUDGETS, KBUDGETS,                                                    &
                         'LIMA', PTHVREFZIKB,                                                           &
                         PRHODREF, PRHODJ, ZT, ZCIT_ELEC,                                               &
                         ZRVT_ELEC, ZRCT_ELEC, ZRRT_ELEC, ZRIT_ELEC, ZRST_ELEC, ZRGT_ELEC,              &
                         ZQPIT, ZQCT, ZQRT, ZQIT, ZQST, ZQGT, ZQNIT,                                    &
                         ZQPIS, ZQCS, ZQRS, ZQIS, ZQSS, ZQGS, ZQNIS,                                    &
                         ZTOT_RI_HIND*ZINV_TSTEP,  -ZTOT_RR_HONR*ZINV_TSTEP,   ZTOT_RC_IMLT*ZINV_TSTEP, &
                        -ZTOT_RC_HONC*ZINV_TSTEP,   ZTOT_RS_DEPS*ZINV_TSTEP,  -ZTOT_RI_AGGS*ZINV_TSTEP, &
                        -ZTOT_RI_CNVS*ZINV_TSTEP,   ZTOT_RG_DEPG*ZINV_TSTEP,  -ZTOT_RC_AUTO*ZINV_TSTEP, &
                        -ZTOT_RC_ACCR*ZINV_TSTEP,  -ZTOT_RR_EVAP*ZINV_TSTEP,                            &
                         ZTOT_RC_RIMSS*ZINV_TSTEP,  ZTOT_RC_RIMSG*ZINV_TSTEP,  ZTOT_RS_RIMCG*ZINV_TSTEP,&
                         ZTOT_RR_ACCSS*ZINV_TSTEP,  ZTOT_RR_ACCSG*ZINV_TSTEP,  ZTOT_RS_ACCRG*ZINV_TSTEP,&
                        -ZTOT_RS_CMEL*ZINV_TSTEP,  -ZTOT_RI_CFRZ*ZINV_TSTEP,  -ZTOT_RR_CFRZ*ZINV_TSTEP, &
                        -ZTOT_RC_WETG*ZINV_TSTEP,  -ZTOT_RI_WETG*ZINV_TSTEP,  -ZTOT_RR_WETG*ZINV_TSTEP, &
                        -ZTOT_RS_WETG*ZINV_TSTEP,                                                       &
                        -ZTOT_RC_DRYG*ZINV_TSTEP,  -ZTOT_RI_DRYG*ZINV_TSTEP,  -ZTOT_RR_DRYG*ZINV_TSTEP, &
                        -ZTOT_RS_DRYG*ZINV_TSTEP,                                                       &
                         ZTOT_RR_GMLT*ZINV_TSTEP,  -ZTOT_RC_BERFI*ZINV_TSTEP,                           &
! variables et processus optionnels propres a lima
                         PCCT=ZCCT, PCRT=ZCRT, PCST=ZCST, PCGT=ZCGT,                                    &
                         PRVHENC=ZTOT_RV_HENU*ZINV_TSTEP,    PRCHINC=-ZTOT_RC_HINC*ZINV_TSTEP,          &
                         PRVHONH=-ZTOT_RV_HONH*ZINV_TSTEP,   PRRCVRC=-ZTOT_RR_CVRC*ZINV_TSTEP,          &
                         PRICNVI=ZTOT_RI_CNVI*ZINV_TSTEP,    PRVDEPI=ZTOT_RI_DEPI*ZINV_TSTEP,           &
                         PRSHMSI=ZTOT_RI_HMS*ZINV_TSTEP,     PRGHMGI=ZTOT_RI_HMG*ZINV_TSTEP,            &
                         PRICIBU=ZTOT_RI_CIBU*ZINV_TSTEP,    PRIRDSF=ZTOT_RI_RDSF*ZINV_TSTEP,           &
                         PRCCORR2=-ZTOT_RC_CORR2*ZINV_TSTEP, PRRCORR2=-ZTOT_RR_CORR2*ZINV_TSTEP,        &
                         PRICORR2=-ZTOT_RI_CORR2*ZINV_TSTEP)            
  END IF
  !
  ! update the source variables
  PSV_ELEC_S(D%NIJB:D%NIJE,:,1) = ZQPIS(D%NIJB:D%NIJE,:)
  PSV_ELEC_S(D%NIJB:D%NIJE,:,2) = ZQCS(D%NIJB:D%NIJE,:)
  PSV_ELEC_S(D%NIJB:D%NIJE,:,3) = ZQRS(D%NIJB:D%NIJE,:)
  PSV_ELEC_S(D%NIJB:D%NIJE,:,4) = ZQIS(D%NIJB:D%NIJE,:)
  PSV_ELEC_S(D%NIJB:D%NIJE,:,5) = ZQSS(D%NIJB:D%NIJE,:)
  PSV_ELEC_S(D%NIJB:D%NIJE,:,6) = ZQGS(D%NIJB:D%NIJE,:)
  IF (KRR == 6) THEN
    PSV_ELEC_S(D%NIJB:D%NIJE,:,7) = ZQNIS(D%NIJB:D%NIJE,:)
  ELSE IF (KRR == 7) THEN
    PSV_ELEC_S(D%NIJB:D%NIJE,:,7) = ZQHS(D%NIJB:D%NIJE,:)
    PSV_ELEC_S(D%NIJB:D%NIJE,:,8) = ZQNIS(D%NIJB:D%NIJE,:)
  END IF
  !
  DEALLOCATE(GMASK_ELEC)
  !
  DEALLOCATE(ZQPIT)
  DEALLOCATE(ZQNIT)
  DEALLOCATE(ZQCT)
  DEALLOCATE(ZQRT)
  DEALLOCATE(ZQIT)
  DEALLOCATE(ZQST)
  DEALLOCATE(ZQGT)
  IF (ALLOCATED(ZQHT)) DEALLOCATE(ZQHT)
  DEALLOCATE(ZQPIS)
  DEALLOCATE(ZQNIS)
  DEALLOCATE(ZQCS)
  DEALLOCATE(ZQRS)
  DEALLOCATE(ZQIS)
  DEALLOCATE(ZQSS)
  DEALLOCATE(ZQGS)
  IF (ALLOCATED(ZQHS)) DEALLOCATE(ZQHS)
  !
  DEALLOCATE(ZRVT_ELEC)
  DEALLOCATE(ZRCT_ELEC)
  DEALLOCATE(ZRRT_ELEC)
  DEALLOCATE(ZRIT_ELEC)
  DEALLOCATE(ZRST_ELEC)
  DEALLOCATE(ZRGT_ELEC)
  IF (ALLOCATED(ZRHT_ELEC)) DEALLOCATE(ZRHT_ELEC)
  IF (ALLOCATED(ZCCT_ELEC)) DEALLOCATE(ZCCT_ELEC)
  IF (ALLOCATED(ZCRT_ELEC)) DEALLOCATE(ZCRT_ELEC)
  IF (ALLOCATED(ZCIT_ELEC)) DEALLOCATE(ZCIT_ELEC)
  IF (ALLOCATED(ZCST_ELEC)) DEALLOCATE(ZCST_ELEC)
  IF (ALLOCATED(ZCGT_ELEC)) DEALLOCATE(ZCGT_ELEC)
  IF (ALLOCATED(ZCHT_ELEC)) DEALLOCATE(ZCHT_ELEC)
  !
END IF
!
DEALLOCATE(ZTOT_RI_HIND)
DEALLOCATE(ZTOT_RC_HINC)
DEALLOCATE(ZTOT_RV_HENU)
DEALLOCATE(ZTOT_RV_HONH)
!
!
!*       7.3    Unpacking variables
!               -------------------
!
! not necessary! the only variables needed in the following (PQxS) are already 3D
!
!
!-------------------------------------------------------------------------------
!
!*       7.     TOTAL TENDENCIES
!               ----------------
!
! Source at the end of microphysics = new state / PTSTEP
!
PTHS(D%NIJB:D%NIJE,:) = ZTHT(D%NIJB:D%NIJE,:) * ZINV_TSTEP
!
PRS(D%NIJB:D%NIJE,:,1) = ZRVT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF ( KRR .GE. 2 ) PRS(D%NIJB:D%NIJE,:,2) = ZRCT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF ( KRR .GE. 3 ) PRS(D%NIJB:D%NIJE,:,3) = ZRRT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF ( KRR .GE. 4 ) PRS(D%NIJB:D%NIJE,:,4) = ZRIT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF ( KRR .GE. 5 ) PRS(D%NIJB:D%NIJE,:,5) = ZRST(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF ( KRR .GE. 6 ) PRS(D%NIJB:D%NIJE,:,6) = ZRGT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF ( KRR .GE. 7 ) PRS(D%NIJB:D%NIJE,:,7) = ZRHT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
!
IF ( LIMAP%NMOM_C.GE.2 ) PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NC) = ZCCT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF ( LIMAP%NMOM_R.GE.2 ) PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NR) = ZCRT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF ( LIMAP%NMOM_I.GE.2 ) THEN
  IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
    PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NI) = ZCIT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
  ELSE
    DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
      PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NI+ISH-1) = ZCIT_SHAPE(D%NIJB:D%NIJE,:,ISH) *ZINV_TSTEP
    END DO
  END IF
END IF
IF ( LIMAP%NMOM_S.GE.2 ) PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NS) = ZCST(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF ( LIMAP%NMOM_G.GE.2 ) PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NG) = ZCGT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
IF ( LIMAP%NMOM_H.GE.2 ) PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_NH) = ZCHT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
!
IF ( LIMAP%NMOD_CCN .GE. 1 )   PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_CCN_FREE:ISV_LIMA_CCN_FREE+LIMAP%NMOD_CCN-1) = ZCCNFT(D%NIJB:D%NIJE,:,:) *ZINV_TSTEP
IF ( LIMAP%NMOD_CCN .GE. 1 )   PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_CCN_ACTI:ISV_LIMA_CCN_ACTI+LIMAP%NMOD_CCN-1) = ZCCNAT(D%NIJB:D%NIJE,:,:) *ZINV_TSTEP
IF ( LIMAP%NMOD_IFN .GE. 1 )   PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_IFN_FREE:ISV_LIMA_IFN_FREE+LIMAP%NMOD_IFN-1) = ZIFNFT(D%NIJB:D%NIJE,:,:) *ZINV_TSTEP
IF ( LIMAP%NMOD_IFN .GE. 1 )   PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_IFN_NUCL:ISV_LIMA_IFN_NUCL+LIMAP%NMOD_IFN-1) = ZIFNNT(D%NIJB:D%NIJE,:,:) *ZINV_TSTEP
IF ( LIMAP%NMOD_IMM .GE. 1 )   PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_IMM_NUCL:ISV_LIMA_IMM_NUCL+LIMAP%NMOD_IMM-1) = ZIMMNT(D%NIJB:D%NIJE,:,:) *ZINV_TSTEP
IF ( LIMAP%LHHONI) PSVS(D%NIJB:D%NIJE,:,ISV_LIMA_HOM_HAZE) = ZHOMFT(D%NIJB:D%NIJE,:) *ZINV_TSTEP
!
IF ( LIMAP%NMOM_I.EQ.1 ) PCIT(:,:) = ZCIT(:,:)*PRHODREF(:,:)
!
!
! Call budgets
!
IF ( BUCONF%LBU_ENABLE ) then
  ALLOCATE(ZRHODJONTSTEP(D%NIJT,D%NKT))
  ZRHODJONTSTEP(D%NIJB:D%NIJE,:) = ZINV_TSTEP * PRHODJ(D%NIJB:D%NIJE,:)

  IF ( BUCONF%LBUDGET_TH ) then
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'REVA',  ZTOT_TH_EVAP (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'HONC',  ZTOT_TH_HONC (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'HONR',  ZTOT_TH_HONR (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'DEPS',  ZTOT_TH_DEPS (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'DEPI',  ZTOT_TH_DEPI (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'DEPG',  ZTOT_TH_DEPG (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'IMLT',  ZTOT_TH_IMLT (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'BERFI', ZTOT_TH_BERFI(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'RIM',   ZTOT_TH_RIM  (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'ACC',   ZTOT_TH_ACC  (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'CFRZ',  ZTOT_TH_CFRZ (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'WETG',  ZTOT_TH_WETG (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'DRYG',  ZTOT_TH_DRYG (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'GMLT',  ZTOT_TH_GMLT (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'DEPH',  ZTOT_TH_DEPH (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'WETH',  ZTOT_TH_WETH (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_TH)%PTR%ADD_PHY(D, 'HMLT',  ZTOT_TH_HMLT (:,:) * ZRHODJONTSTEP(:,:) )
  END IF

  IF ( BUCONF%LBUDGET_RV ) then
    CALL TBUDGETS(NBUDGET_RV)%PTR%ADD_PHY(D, 'REVA', -ZTOT_RR_EVAP (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RV)%PTR%ADD_PHY(D, 'DEPS', -ZTOT_RS_DEPS (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RV)%PTR%ADD_PHY(D, 'DEPI', -ZTOT_RI_DEPI (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RV)%PTR%ADD_PHY(D, 'DEPG', -ZTOT_RG_DEPG (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RV)%PTR%ADD_PHY(D, 'CORR2', ZTOT_RV_CORR2(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RV)%PTR%ADD_PHY(D, 'DEPH', -ZTOT_RH_DEPH (:,:) * ZRHODJONTSTEP(:,:) )
  END IF

  IF ( BUCONF%LBUDGET_RC ) then
    CALL TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'AUTO',  ZTOT_RC_AUTO (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'ACCR',  ZTOT_RC_ACCR (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'HONC',  ZTOT_RC_HONC (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'IMLT',  ZTOT_RC_IMLT (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'BERFI', ZTOT_RC_BERFI(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'RIM',  (-ZTOT_RR_ACCSS(:,:) - ZTOT_RR_ACCSG(:,:)) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'WETG',  ZTOT_RC_WETG (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'DRYG',  ZTOT_RC_DRYG (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'CVRC', -ZTOT_RR_CVRC (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'CORR2', ZTOT_RC_CORR2(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RC)%PTR%ADD_PHY(D, 'WETH',  ZTOT_RC_WETH(:,:)  * ZRHODJONTSTEP(:,:) )
  END IF

  IF ( BUCONF%LBUDGET_RR ) then
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'AUTO', -ZTOT_RC_AUTO(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'ACCR', -ZTOT_RC_ACCR(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'REVA',  ZTOT_RR_EVAP(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'HONR',  ZTOT_RR_HONR(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'ACC',  (ZTOT_RC_RIMSS(:,:) + ZTOT_RC_RIMSG(:,:)) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'CFRZ',  ZTOT_RR_CFRZ(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'WETG',  ZTOT_RR_WETG(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'DRYG',  ZTOT_RR_DRYG(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'GMLT',  ZTOT_RR_GMLT(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'CVRC',  ZTOT_RR_CVRC(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'CORR2', ZTOT_RR_CORR2(:,:)* ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'WETH',  ZTOT_RR_WETH(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RR)%PTR%ADD_PHY(D, 'HMLT',  ZTOT_RR_HMLT(:,:) * ZRHODJONTSTEP(:,:) )
  END IF

  IF ( BUCONF%LBUDGET_RI ) then
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'HONC',  -ZTOT_RC_HONC (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'CNVI',   ZTOT_RI_CNVI (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'CNVS',   ZTOT_RI_CNVS (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'AGGS',   ZTOT_RI_AGGS (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'IMLT',  -ZTOT_RC_IMLT (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'BERFI', -ZTOT_RC_BERFI(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'HMS',    ZTOT_RI_HMS  (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'CFRZ',   ZTOT_RI_CFRZ (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'DEPI',   ZTOT_RI_DEPI (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'CIBU',   ZTOT_RI_CIBU (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'RDSF',   ZTOT_RI_RDSF (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'WETG',   ZTOT_RI_WETG (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'DRYG',   ZTOT_RI_DRYG (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'HMG',    ZTOT_RI_HMG  (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'CORR2',  ZTOT_RI_CORR2(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'WETH',   ZTOT_RI_WETH (:,:) * ZRHODJONTSTEP(:,:) )
    IF (LIMAP%LCRYSTAL_SHAPE) THEN
      CALL TBUDGETS(NBUDGET_RI)%PTR%ADD_PHY(D, 'ISCS', ZTOT_SHRI_ISCS (:,:) * ZRHODJONTSTEP(:,:) )
    END IF
  END IF

  IF ( BUCONF%LBUDGET_RS ) then
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'CNVI', -ZTOT_RI_CNVI(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'DEPS',  ZTOT_RS_DEPS(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'CNVS', -ZTOT_RI_CNVS(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'AGGS', -ZTOT_RI_AGGS(:,:) * ZRHODJONTSTEP(:,:) )
!    call TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'RIM',   ztot_rs_rim (:,:) * zrhodjontstep(:,:) )
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'RIM', (-ZTOT_RC_RIMSS(:,:) -  ZTOT_RS_RIMCG(:,:)) & 
                                                          * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'HMS',   ZTOT_RS_HMS (:,:) * ZRHODJONTSTEP(:,:) )
!    call TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'ACC',   ztot_rs_acc (:,:) * zrhodjontstep(:,:) )
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'ACC',  (ZTOT_RR_ACCSS(:,:) - ZTOT_RS_ACCRG (:,:)) &
                                                          * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'CMEL',  ZTOT_RS_CMEL(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'CIBU', -ZTOT_RI_CIBU(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'WETG',  ZTOT_RS_WETG(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'DRYG',  ZTOT_RS_DRYG(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'WETH',  ZTOT_RS_WETH(:,:) * ZRHODJONTSTEP(:,:) )
    IF (LIMAP%LCRYSTAL_SHAPE) THEN
      CALL TBUDGETS(NBUDGET_RS)%PTR%ADD_PHY(D, 'ISCS', -ZTOT_SHRI_ISCS (:,:) * ZRHODJONTSTEP(:,:) )
    END IF
  END IF

  IF ( BUCONF%LBUDGET_RG ) then
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'HONR', -ZTOT_RR_HONR(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'DEPG',  ZTOT_RG_DEPG(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'RIM',  (-ZTOT_RC_RIMSG(:,:) +  ZTOT_RS_RIMCG(:,:) ) &
                                                          * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'ACC',  (ZTOT_RR_ACCSG(:,:) + ZTOT_RS_ACCRG (:,:)) & 
                                                          * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'CMEL', -ZTOT_RS_CMEL(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'CFRZ', ( -ZTOT_RR_CFRZ(:,:) - ZTOT_RI_CFRZ(:,:) ) &
                                                         * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'RDSF', -ZTOT_RI_RDSF(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'WETG',  ZTOT_RG_WETG(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'DRYG',  ZTOT_RG_DRYG(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'HMG',   ZTOT_RG_HMG (:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'GMLT', -ZTOT_RR_GMLT(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'WETH',  ZTOT_RG_WETH(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RG)%PTR%ADD_PHY(D, 'COHG',  ZTOT_RG_COHG(:,:) * ZRHODJONTSTEP(:,:) )
  END IF

  IF ( BUCONF%LBUDGET_RH ) then
    CALL TBUDGETS(NBUDGET_RH)%PTR%ADD_PHY(D, 'WETG',  ZTOT_RH_WETG(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RH)%PTR%ADD_PHY(D, 'DEPH',  ZTOT_RH_DEPH(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RH)%PTR%ADD_PHY(D, 'WETH',  ZTOT_RH_WETH(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RH)%PTR%ADD_PHY(D, 'COHG', -ZTOT_RG_COHG(:,:) * ZRHODJONTSTEP(:,:) )
    CALL TBUDGETS(NBUDGET_RH)%PTR%ADD_PHY(D, 'HMLT', -ZTOT_RR_HMLT(:,:) * ZRHODJONTSTEP(:,:) )
  END IF

  IF ( BUCONF%LBUDGET_SV ) then
    !
    ! Cloud droplets
    !
    IF (LIMAP%NMOM_C.GE.2) then 
       IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NC
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'SELF',  ZTOT_CC_SELF (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'AUTO',  ZTOT_CC_AUTO (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'ACCR',  ZTOT_CC_ACCR (:,:) * ZRHODJONTSTEP(:,:) )
       !call TBUDGETS(IDX)%PTR%ADD_PHY(D, 'REVA',  0. )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'HONC',  ZTOT_CC_HONC (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'IMLT',  ZTOT_CC_IMLT (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'RIM',   ZTOT_CC_RIM  (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETG',  ZTOT_CC_WETG (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'DRYG',  ZTOT_CC_DRYG (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CVRC', -ZTOT_CR_CVRC (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CORR2', ZTOT_CC_CORR2(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETH',  ZTOT_CC_WETH (:,:) * ZRHODJONTSTEP(:,:) )
    END IF
    !
    ! Rain drops
    !
    IF (LIMAP%NMOM_R.GE.2) then
       IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NR
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'AUTO',  ZTOT_CR_AUTO(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'SCBU',  ZTOT_CR_SCBU(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'REVA',  ZTOT_CR_EVAP(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'BRKU',  ZTOT_CR_BRKU(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'HONR',  ZTOT_CR_HONR(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'ACC',   ZTOT_CR_ACC (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CFRZ',  ZTOT_CR_CFRZ(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETG',  ZTOT_CR_WETG(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'DRYG',  ZTOT_CR_DRYG(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'GMLT',  ZTOT_CR_GMLT(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CVRC',  ZTOT_CR_CVRC(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CORR2', ZTOT_CR_CORR2(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETH',  ZTOT_CR_WETH(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'HMLT',  ZTOT_CR_HMLT(:,:) * ZRHODJONTSTEP(:,:) )
    END IF
    !
    ! Ice crystals
    !
    IF (LIMAP%NMOM_I.GE.2) then
      IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
       IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NI
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'HONC',  -ZTOT_CC_HONC (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CNVI',   ZTOT_CI_CNVI (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CNVS',   ZTOT_CI_CNVS (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'AGGS',   ZTOT_CI_AGGS (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'IMLT',  -ZTOT_CC_IMLT (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'HMS',    ZTOT_CI_HMS  (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CFRZ',   ZTOT_CI_CFRZ (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CIBU',   ZTOT_CI_CIBU (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'RDSF',   ZTOT_CI_RDSF (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETG',   ZTOT_CI_WETG (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'DRYG',   ZTOT_CI_DRYG (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'HMG',    ZTOT_CI_HMG  (:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CORR2',  ZTOT_CI_CORR2(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETH',   ZTOT_CI_WETH (:,:) * ZRHODJONTSTEP(:,:) )
      ELSE
         DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
           IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NI + ISH - 1
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'HONC', ztot_shci_honc (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CNVI', ztot_shci_cnvi (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'HACH', ztot_shci_hach (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CNVS', ztot_shci_cnvs (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'AGGS', ztot_shci_aggs (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'IMLT', ztot_shci_imlt (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'HMS',  ztot_shci_hms  (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CFRZ', ztot_shci_cfrz (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CIBU', ztot_shci_cibu (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'RDSF', ztot_shci_rdsf (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETG', ztot_shci_wetg (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'DRYG', ztot_shci_dryg (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'HMG',  ztot_shci_hmg  (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CORR2',ztot_shci_corr2(:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'ISC',  ztot_shci_isc  (:,:,ISH) * ZRHODJONTSTEP(:,:) )
           CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'ISCS', ztot_shci_iscs (:,:,ISH) * ZRHODJONTSTEP(:,:) )
         END DO
      END IF
    END IF
    !
    ! Snow
    !
    IF (LIMAP%NMOM_S.GE.2) then
       IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NS
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CNVI',  -ZTOT_CI_CNVI(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CNVS',  -ZTOT_CI_CNVS(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'RIM',    ZTOT_CS_RIM(:,:)  * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'ACC',    ZTOT_CS_ACC(:,:)  * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CMEL',   ZTOT_CS_CMEL(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'SSC',    ZTOT_CS_SSC(:,:)  * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETG',   ZTOT_CS_WETG(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'DRYG',   ZTOT_CS_DRYG(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETH',   ZTOT_CS_WETH(:,:) * ZRHODJONTSTEP(:,:) )
       IF (LIMAP%LCRYSTAL_SHAPE) THEN
         ALLOCATE(ZAUX_ISCS(D%NIJT,D%NKT)) ; ZAUX_ISCS(:,:) = 0.    
         DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
           ZAUX_ISCS(:,:) = SUM(ZTOT_SHCI_ISCS,DIM=3)
         END DO
         CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'ISCS',  -ZAUX_ISCS(:,:) * ZRHODJONTSTEP(:,:) ) 
         DEALLOCATE(ZAUX_ISCS)
       END IF
    END IF
    !
    ! Graupel
    !
    IF (LIMAP%NMOM_G.GE.2) then
       IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NG
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'RIM',   -ZTOT_CS_RIM(:,:)  * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'ACC',   -ZTOT_CS_ACC(:,:)  * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CMEL',  -ZTOT_CS_CMEL(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'CFRZ',  -ZTOT_CR_CFRZ(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETG',   ZTOT_CG_WETG(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'GMLT',   ZTOT_CG_GMLT(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETH',   ZTOT_CG_WETH(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'COHG',   ZTOT_CG_COHG(:,:) * ZRHODJONTSTEP(:,:) )
    END IF
    !
    ! Hail
    !
    IF (LIMAP%NMOM_H.GE.2) then
       IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NH
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'WETG',  -ZTOT_CG_WETG(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'COHG',  -ZTOT_CG_COHG(:,:) * ZRHODJONTSTEP(:,:) )
       CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'HMLT',   ZTOT_CH_HMLT(:,:) * ZRHODJONTSTEP(:,:) )
    END IF

    DO II = 1, LIMAP%NMOD_IFN
      IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_IFN_NUCL + II - 1
      CALL TBUDGETS(IDX)%PTR%ADD_PHY(D, 'IMLT', ZTOT_IFNN_IMLT(:, :, II) * ZRHODJONTSTEP(:,:) )
    END DO
  END IF

  DEALLOCATE( ZRHODJONTSTEP )
END IF
!
IF (ALLOCATED(ZCIT_SHAPE)) DEALLOCATE(ZCIT_SHAPE)
IF (ALLOCATED(ZCIS_SHAPE)) DEALLOCATE(ZCIS_SHAPE)
IF (ALLOCATED(ZRIT_SHAPE)) DEALLOCATE(ZRIT_SHAPE)
IF (ALLOCATED(ZRIS_SHAPE)) DEALLOCATE(ZRIS_SHAPE)
IF (ALLOCATED(ZTOT_CR_BRKU)) DEALLOCATE(ZTOT_CR_BRKU)
IF (ALLOCATED(ZTOT_TH_HONR)) DEALLOCATE(ZTOT_TH_HONR)
IF (ALLOCATED(ZTOT_RR_HONR)) DEALLOCATE(ZTOT_RR_HONR)
IF (ALLOCATED(ZTOT_CR_HONR)) DEALLOCATE(ZTOT_CR_HONR)
IF (ALLOCATED(ZTOT_TH_IMLT)) DEALLOCATE(ZTOT_TH_IMLT)
IF (ALLOCATED(ZTOT_RC_IMLT)) DEALLOCATE(ZTOT_RC_IMLT)
IF (ALLOCATED(ZTOT_CC_IMLT)) DEALLOCATE(ZTOT_CC_IMLT)
IF (ALLOCATED(ZTOT_TH_HONC)) DEALLOCATE(ZTOT_TH_HONC)
IF (ALLOCATED(ZTOT_RC_HONC)) DEALLOCATE(ZTOT_RC_HONC)
IF (ALLOCATED(ZTOT_CC_HONC)) DEALLOCATE(ZTOT_CC_HONC)
IF (ALLOCATED(ZTOT_CC_SELF)) DEALLOCATE(ZTOT_CC_SELF)
IF (ALLOCATED(ZTOT_RC_AUTO)) DEALLOCATE(ZTOT_RC_AUTO)
IF (ALLOCATED(ZTOT_CC_AUTO)) DEALLOCATE(ZTOT_CC_AUTO)
IF (ALLOCATED(ZTOT_CR_AUTO)) DEALLOCATE(ZTOT_CR_AUTO)
IF (ALLOCATED(ZTOT_RC_ACCR)) DEALLOCATE(ZTOT_RC_ACCR)
IF (ALLOCATED(ZTOT_CC_ACCR)) DEALLOCATE(ZTOT_CC_ACCR)
IF (ALLOCATED(ZTOT_CR_SCBU)) DEALLOCATE(ZTOT_CR_SCBU)
IF (ALLOCATED(ZTOT_TH_EVAP)) DEALLOCATE(ZTOT_TH_EVAP)
IF (ALLOCATED(ZTOT_RR_EVAP)) DEALLOCATE(ZTOT_RR_EVAP)
IF (ALLOCATED(ZTOT_CR_EVAP)) DEALLOCATE(ZTOT_CR_EVAP)
IF (ALLOCATED(ZTOT_RI_CNVI)) DEALLOCATE(ZTOT_RI_CNVI)
IF (ALLOCATED(ZTOT_CI_CNVI)) DEALLOCATE(ZTOT_CI_CNVI)
IF (ALLOCATED(ZTOT_TH_DEPS)) DEALLOCATE(ZTOT_TH_DEPS)
IF (ALLOCATED(ZTOT_RS_DEPS)) DEALLOCATE(ZTOT_RS_DEPS)
IF (ALLOCATED(ZTOT_TH_DEPI)) DEALLOCATE(ZTOT_TH_DEPI)
IF (ALLOCATED(ZTOT_RI_DEPI)) DEALLOCATE(ZTOT_RI_DEPI)
IF (ALLOCATED(ZTOT_RI_CNVS)) DEALLOCATE(ZTOT_RI_CNVS)
IF (ALLOCATED(ZTOT_CI_CNVS)) DEALLOCATE(ZTOT_CI_CNVS)
IF (ALLOCATED(ZTOT_CS_SSC)) DEALLOCATE(ZTOT_CS_SSC)
IF (ALLOCATED(ZTOT_CI_ISC)) DEALLOCATE(ZTOT_CI_ISC)
IF (ALLOCATED(ZTOT_RI_AGGS)) DEALLOCATE(ZTOT_RI_AGGS)
IF (ALLOCATED(ZTOT_CI_AGGS)) DEALLOCATE(ZTOT_CI_AGGS)
IF (ALLOCATED(ZTOT_TH_DEPG)) DEALLOCATE(ZTOT_TH_DEPG)
IF (ALLOCATED(ZTOT_RG_DEPG)) DEALLOCATE(ZTOT_RG_DEPG)
IF (ALLOCATED(ZTOT_TH_BERFI)) DEALLOCATE(ZTOT_TH_BERFI)
IF (ALLOCATED(ZTOT_RC_BERFI)) DEALLOCATE(ZTOT_RC_BERFI)
IF (ALLOCATED(ZTOT_TH_RIM)) DEALLOCATE(ZTOT_TH_RIM)
IF (ALLOCATED(ZTOT_CC_RIM)) DEALLOCATE(ZTOT_CC_RIM)
IF (ALLOCATED(ZTOT_CS_RIM)) DEALLOCATE(ZTOT_CS_RIM)
IF (ALLOCATED(ZTOT_RC_RIMSS)) DEALLOCATE(ZTOT_RC_RIMSS)
IF (ALLOCATED(ZTOT_RC_RIMSG)) DEALLOCATE(ZTOT_RC_RIMSG)
IF (ALLOCATED(ZTOT_RS_RIMCG)) DEALLOCATE(ZTOT_RS_RIMCG)
IF (ALLOCATED(ZTOT_RI_HMS)) DEALLOCATE(ZTOT_RI_HMS)
IF (ALLOCATED(ZTOT_CI_HMS)) DEALLOCATE(ZTOT_CI_HMS)
IF (ALLOCATED(ZTOT_RS_HMS)) DEALLOCATE(ZTOT_RS_HMS)
IF (ALLOCATED(ZTOT_TH_ACC)) DEALLOCATE(ZTOT_TH_ACC)
IF (ALLOCATED(ZTOT_CR_ACC)) DEALLOCATE(ZTOT_CR_ACC)
IF (ALLOCATED(ZTOT_CS_ACC)) DEALLOCATE(ZTOT_CS_ACC)
IF (ALLOCATED(ZTOT_RR_ACCSS)) DEALLOCATE(ZTOT_RR_ACCSS)
IF (ALLOCATED(ZTOT_RR_ACCSG)) DEALLOCATE(ZTOT_RR_ACCSG)
IF (ALLOCATED(ZTOT_RS_ACCRG)) DEALLOCATE(ZTOT_RS_ACCRG)
IF (ALLOCATED(ZTOT_RS_CMEL)) DEALLOCATE(ZTOT_RS_CMEL)
IF (ALLOCATED(ZTOT_CS_CMEL)) DEALLOCATE(ZTOT_CS_CMEL)
IF (ALLOCATED(ZTOT_TH_CFRZ)) DEALLOCATE(ZTOT_TH_CFRZ)
IF (ALLOCATED(ZTOT_RR_CFRZ)) DEALLOCATE(ZTOT_RR_CFRZ)
IF (ALLOCATED(ZTOT_CR_CFRZ)) DEALLOCATE(ZTOT_CR_CFRZ)
IF (ALLOCATED(ZTOT_RI_CFRZ)) DEALLOCATE(ZTOT_RI_CFRZ)
IF (ALLOCATED(ZTOT_CI_CFRZ)) DEALLOCATE(ZTOT_CI_CFRZ)
IF (ALLOCATED(ZTOT_RI_CIBU)) DEALLOCATE(ZTOT_RI_CIBU)
IF (ALLOCATED(ZTOT_CI_CIBU)) DEALLOCATE(ZTOT_CI_CIBU)
IF (ALLOCATED(ZTOT_RI_RDSF)) DEALLOCATE(ZTOT_RI_RDSF)
IF (ALLOCATED(ZTOT_CI_RDSF)) DEALLOCATE(ZTOT_CI_RDSF)
IF (ALLOCATED(ZTOT_TH_WETG)) DEALLOCATE(ZTOT_TH_WETG)
IF (ALLOCATED(ZTOT_RC_WETG)) DEALLOCATE(ZTOT_RC_WETG)
IF (ALLOCATED(ZTOT_CC_WETG)) DEALLOCATE(ZTOT_CC_WETG)
IF (ALLOCATED(ZTOT_RR_WETG)) DEALLOCATE(ZTOT_RR_WETG)
IF (ALLOCATED(ZTOT_CR_WETG)) DEALLOCATE(ZTOT_CR_WETG)
IF (ALLOCATED(ZTOT_RI_WETG)) DEALLOCATE(ZTOT_RI_WETG)
IF (ALLOCATED(ZTOT_CI_WETG)) DEALLOCATE(ZTOT_CI_WETG)
IF (ALLOCATED(ZTOT_RS_WETG)) DEALLOCATE(ZTOT_RS_WETG)
IF (ALLOCATED(ZTOT_CS_WETG)) DEALLOCATE(ZTOT_CS_WETG)
IF (ALLOCATED(ZTOT_RG_WETG)) DEALLOCATE(ZTOT_RG_WETG)
IF (ALLOCATED(ZTOT_CG_WETG)) DEALLOCATE(ZTOT_CG_WETG)
IF (ALLOCATED(ZTOT_RH_WETG)) DEALLOCATE(ZTOT_RH_WETG)
IF (ALLOCATED(ZTOT_TH_DRYG)) DEALLOCATE(ZTOT_TH_DRYG)
IF (ALLOCATED(ZTOT_RC_DRYG)) DEALLOCATE(ZTOT_RC_DRYG)
IF (ALLOCATED(ZTOT_CC_DRYG)) DEALLOCATE(ZTOT_CC_DRYG)
IF (ALLOCATED(ZTOT_RR_DRYG)) DEALLOCATE(ZTOT_RR_DRYG)
IF (ALLOCATED(ZTOT_CR_DRYG)) DEALLOCATE(ZTOT_CR_DRYG)
IF (ALLOCATED(ZTOT_RI_DRYG)) DEALLOCATE(ZTOT_RI_DRYG)
IF (ALLOCATED(ZTOT_CI_DRYG)) DEALLOCATE(ZTOT_CI_DRYG)
IF (ALLOCATED(ZTOT_RS_DRYG)) DEALLOCATE(ZTOT_RS_DRYG)
IF (ALLOCATED(ZTOT_CS_DRYG)) DEALLOCATE(ZTOT_CS_DRYG)
IF (ALLOCATED(ZTOT_RG_DRYG)) DEALLOCATE(ZTOT_RG_DRYG)
IF (ALLOCATED(ZTOT_RI_HMG)) DEALLOCATE(ZTOT_RI_HMG)
IF (ALLOCATED(ZTOT_CI_HMG)) DEALLOCATE(ZTOT_CI_HMG)
IF (ALLOCATED(ZTOT_RG_HMG)) DEALLOCATE(ZTOT_RG_HMG)
IF (ALLOCATED(ZTOT_TH_GMLT)) DEALLOCATE(ZTOT_TH_GMLT)
IF (ALLOCATED(ZTOT_RR_GMLT)) DEALLOCATE(ZTOT_RR_GMLT)
IF (ALLOCATED(ZTOT_CR_GMLT)) DEALLOCATE(ZTOT_CR_GMLT)
IF (ALLOCATED(ZTOT_CG_GMLT)) DEALLOCATE(ZTOT_CG_GMLT)
IF (ALLOCATED(ZTOT_TH_DEPH)) DEALLOCATE(ZTOT_TH_DEPH)
IF (ALLOCATED(ZTOT_RH_DEPH)) DEALLOCATE(ZTOT_RH_DEPH)
IF (ALLOCATED(ZTOT_TH_WETH)) DEALLOCATE(ZTOT_TH_WETH)
IF (ALLOCATED(ZTOT_RC_WETH)) DEALLOCATE(ZTOT_RC_WETH)
IF (ALLOCATED(ZTOT_CC_WETH)) DEALLOCATE(ZTOT_CC_WETH)
IF (ALLOCATED(ZTOT_RR_WETH)) DEALLOCATE(ZTOT_RR_WETH)
IF (ALLOCATED(ZTOT_CR_WETH)) DEALLOCATE(ZTOT_CR_WETH)
IF (ALLOCATED(ZTOT_RI_WETH)) DEALLOCATE(ZTOT_RI_WETH)
IF (ALLOCATED(ZTOT_CI_WETH)) DEALLOCATE(ZTOT_CI_WETH)
IF (ALLOCATED(ZTOT_RS_WETH)) DEALLOCATE(ZTOT_RS_WETH)
IF (ALLOCATED(ZTOT_CS_WETH)) DEALLOCATE(ZTOT_CS_WETH)
IF (ALLOCATED(ZTOT_RG_WETH)) DEALLOCATE(ZTOT_RG_WETH)
IF (ALLOCATED(ZTOT_CG_WETH)) DEALLOCATE(ZTOT_CG_WETH)
IF (ALLOCATED(ZTOT_RH_WETH)) DEALLOCATE(ZTOT_RH_WETH)
IF (ALLOCATED(ZTOT_RG_COHG)) DEALLOCATE(ZTOT_RG_COHG)
IF (ALLOCATED(ZTOT_CG_COHG)) DEALLOCATE(ZTOT_CG_COHG)
IF (ALLOCATED(ZTOT_TH_HMLT)) DEALLOCATE(ZTOT_TH_HMLT)
IF (ALLOCATED(ZTOT_RR_HMLT)) DEALLOCATE(ZTOT_RR_HMLT)
IF (ALLOCATED(ZTOT_CR_HMLT)) DEALLOCATE(ZTOT_CR_HMLT)
IF (ALLOCATED(ZTOT_CH_HMLT)) DEALLOCATE(ZTOT_CH_HMLT)
IF (ALLOCATED(ZTOT_RR_CVRC)) DEALLOCATE(ZTOT_RR_CVRC)
IF (ALLOCATED(ZTOT_CR_CVRC)) DEALLOCATE(ZTOT_CR_CVRC)
IF (ALLOCATED(ZTOT_RV_CORR2)) DEALLOCATE(ZTOT_RV_CORR2)
IF (ALLOCATED(ZTOT_RC_CORR2)) DEALLOCATE(ZTOT_RC_CORR2)
IF (ALLOCATED(ZTOT_RR_CORR2)) DEALLOCATE(ZTOT_RR_CORR2)
IF (ALLOCATED(ZTOT_RI_CORR2)) DEALLOCATE(ZTOT_RI_CORR2)
IF (ALLOCATED(ZTOT_CC_CORR2)) DEALLOCATE(ZTOT_CC_CORR2)
IF (ALLOCATED(ZTOT_CR_CORR2)) DEALLOCATE(ZTOT_CR_CORR2)
IF (ALLOCATED(ZTOT_CI_CORR2)) DEALLOCATE(ZTOT_CI_CORR2)
IF (ALLOCATED(ZTOT_IFNN_IMLT)) DEALLOCATE(ZTOT_IFNN_IMLT)
IF (ALLOCATED(ZTOT_SHCI_HONC)) DEALLOCATE(ZTOT_SHCI_HONC)
IF (ALLOCATED(ZTOT_SHCI_CNVI)) DEALLOCATE(ZTOT_SHCI_CNVI)
IF (ALLOCATED(ZTOT_SHCI_HMS)) DEALLOCATE(ZTOT_SHCI_HMS)
IF (ALLOCATED(ZTOT_SHCI_HMG)) DEALLOCATE(ZTOT_SHCI_HMG)
IF (ALLOCATED(ZTOT_SHCI_CFRZ)) DEALLOCATE(ZTOT_SHCI_CFRZ)
IF (ALLOCATED(ZTOT_SHCI_CIBU)) DEALLOCATE(ZTOT_SHCI_CIBU)
IF (ALLOCATED(ZTOT_SHCI_RDSF)) DEALLOCATE(ZTOT_SHCI_RDSF)
IF (ALLOCATED(ZTOT_SHCI_WETG)) DEALLOCATE(ZTOT_SHCI_WETG)
IF (ALLOCATED(ZTOT_SHCI_DRYG)) DEALLOCATE(ZTOT_SHCI_DRYG)
IF (ALLOCATED(ZTOT_SHCI_CORR2)) DEALLOCATE(ZTOT_SHCI_CORR2)
IF (ALLOCATED(ZTOT_SHCI_IMLT)) DEALLOCATE(ZTOT_SHCI_IMLT)
IF (ALLOCATED(ZTOT_SHCI_HACH)) DEALLOCATE(ZTOT_SHCI_HACH)
IF (ALLOCATED(ZTOT_SHCI_CNVS)) DEALLOCATE(ZTOT_SHCI_CNVS)
IF (ALLOCATED(ZTOT_SHCI_AGGS)) DEALLOCATE(ZTOT_SHCI_AGGS)
IF (ALLOCATED(ZTOT_SHCI_ISC)) DEALLOCATE(ZTOT_SHCI_ISC)
IF (ALLOCATED(ZTOT_SHCI_ISCS)) DEALLOCATE(ZTOT_SHCI_ISCS)
IF (ALLOCATED(ZTOT_SHRI_ISCS)) DEALLOCATE(ZTOT_SHRI_ISCS)
!
IF (LHOOK) CALL DR_HOOK('LIMA', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA
