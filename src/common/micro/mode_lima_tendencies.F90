!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_TENDENCIES
  IMPLICIT NONE
CONTAINS
!#####################################################################
  SUBROUTINE LIMA_TENDENCIES (PTSTEP, LDCOMPUTE,                                                &
                              PEXNREF, PRHODREF, PPABST, PTHT,                                  &
                              PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT,                         &
                              PCCT, PCRT, PCIT, PCST, PCGT, PCHT,                               &
                              P_TH_HONC, P_RC_HONC, P_CC_HONC,                                  & 
                              P_CC_SELF,                                                        & 
                              P_RC_AUTO, P_CC_AUTO, P_CR_AUTO,                                  & 
                              P_RC_ACCR, P_CC_ACCR,                                             & 
                              P_CR_SCBU,                                                        & 
                              P_TH_EVAP, P_RR_EVAP, P_CR_EVAP,                                  & 
                              P_RI_CNVI, P_CI_CNVI,                                             & 
                              P_TH_DEPS, P_RS_DEPS,                                             & 
                              P_TH_DEPI, P_RI_DEPI,                                             & 
                              P_RI_CNVS, P_CI_CNVS,                                             &
                              P_CS_SSC,                                                         &
                              P_RI_AGGS, P_CI_AGGS,                                             & 
                              P_TH_DEPG, P_RG_DEPG,                                             & 
                              P_TH_BERFI, P_RC_BERFI,                                           & 
!++cb++
!                              P_TH_RIM, P_RC_RIM, P_CC_RIM, P_RS_RIM, P_CS_RIM, P_RG_RIM,       & 
                              P_TH_RIM, P_CC_RIM, P_CS_RIM, P_RC_RIMSS, P_RC_RIMSG, P_RS_RIMCG, &
                              P_RI_HMS, P_CI_HMS, P_RS_HMS,                                     & 
!                              P_TH_ACC, P_RR_ACC, P_CR_ACC, P_RS_ACC, P_CS_ACC, P_RG_ACC,       & 
                              P_TH_ACC, P_CR_ACC, P_CS_ACC, P_RR_ACCSS, P_RR_ACCSG, P_RS_ACCRG, &
!--cb--
                              P_RS_CMEL, P_CS_CMEL,                                             & 
                              P_TH_CFRZ, P_RR_CFRZ, P_CR_CFRZ, P_RI_CFRZ, P_CI_CFRZ,            & 
                              P_RI_CIBU, P_CI_CIBU,                                             & 
                              P_RI_RDSF, P_CI_RDSF,                                             & 
                              P_TH_WETG, P_RC_WETG, P_CC_WETG, P_RR_WETG, P_CR_WETG,            & 
                              P_RI_WETG, P_CI_WETG, P_RS_WETG, P_CS_WETG, P_RG_WETG, P_CG_WETG, P_RH_WETG, & 
                              P_TH_DRYG, P_RC_DRYG, P_CC_DRYG, P_RR_DRYG, P_CR_DRYG,            & 
                              P_RI_DRYG, P_CI_DRYG, P_RS_DRYG, P_CS_DRYG, P_RG_DRYG,            & 
                              P_RI_HMG, P_CI_HMG, P_RG_HMG,                                     & 
                              P_TH_GMLT, P_RR_GMLT, P_CR_GMLT, P_CG_GMLT,                       &
                              P_TH_DEPH, P_RH_DEPH,                                             &
                              P_TH_WETH, P_RC_WETH, P_CC_WETH, P_RR_WETH, P_CR_WETH,            &
                              P_RI_WETH, P_CI_WETH, P_RS_WETH, P_CS_WETH, P_RG_WETH, P_CG_WETH, P_RH_WETH, &
                              P_RG_COHG, P_CG_COHG,                                             &
                              P_TH_HMLT, P_RR_HMLT, P_CR_HMLT, P_CH_HMLT,                       &
                              PA_TH, PA_RV, PA_RC, PA_CC, PA_RR, PA_CR,                         &
                              PA_RI, PA_CI, PA_RS, PA_CS, PA_RG, PA_CG, PA_RH, PA_CH,           &
                              PEVAP3D,                                                          &
                              PCF1D, PIF1D, PPF1D,                                              &
                              PLATHAM_IAGGS                                                     )
!     ######################################################################
!!
!!    PURPOSE
!!    -------
!!      Compute sources of non-instantaneous microphysical processes for the
!!    time-split version of LIMA
!!
!!    AUTHOR
!!    ------
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!!
!       Delbeke/Vie     03/2022 : KHKO option
!       J. Wurtz        03/2022 : new snow characteristics
!       B. Vie          03/2022: Add option for 1-moment pristine ice
!       C. Barthe       06/2022: change some mass transfer rates to be consistent with ICE3, for cloud electrification
!       C. Barthe       06/2023: add Latham effet for IAGGS
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY : XP00, XRD, XRV, XMD, XMV, XCPD, XCPV, XCL, XCI, XLVTT, XLSTT, XTT, &
                                  XALPW, XBETAW, XGAMW, XALPI, XBETAI, XGAMI
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCTMIN,  LCIBU, LRDSF, LKHKO, LSNOW_T, &
                                  NMOM_C, NMOM_R, NMOM_I, NMOM_S, NMOM_G, NMOM_H
USE MODD_PARAM_LIMA_WARM,  ONLY : XLBC, XLBEXC, XLBR, XLBEXR, XCCR, XCXR
USE MODD_PARAM_LIMA_MIXED, ONLY : XLBG, XLBEXG, XCCG, XCXG, XLBH, XLBEXH, XCCH, XCXH
USE MODD_PARAM_LIMA_COLD,  ONLY : XSCFAC, XLBI, XLBEXI, XLBS, XLBEXS, XLBDAS_MAX, XTRANS_MP_GAMMAS,  &
                                   XLBDAS_MIN, XCCS, XCXS, XBS, XNS
!
USE MODE_LIMA_DROPLETS_HOM_FREEZING, ONLY: LIMA_DROPLETS_HOM_FREEZING
USE MODE_LIMA_DROPLETS_SELF_COLLECTION, ONLY: LIMA_DROPLETS_SELF_COLLECTION
USE MODE_LIMA_DROPLETS_AUTOCONVERSION, ONLY: LIMA_DROPLETS_AUTOCONVERSION
USE MODE_LIMA_DROPLETS_ACCRETION, ONLY: LIMA_DROPLETS_ACCRETION
USE MODE_LIMA_DROPS_SELF_COLLECTION, ONLY: LIMA_DROPS_SELF_COLLECTION
USE MODE_LIMA_RAIN_EVAPORATION, ONLY: LIMA_RAIN_EVAPORATION
USE MODE_LIMA_ICE_DEPOSITION, ONLY: LIMA_ICE_DEPOSITION
USE MODE_LIMA_SNOW_DEPOSITION, ONLY: LIMA_SNOW_DEPOSITION
USE MODE_LIMA_SNOW_SELF_COLLECTION, ONLY: LIMA_SNOW_SELF_COLLECTION
USE MODE_LIMA_ICE_AGGREGATION_SNOW, ONLY: LIMA_ICE_AGGREGATION_SNOW
USE MODE_LIMA_GRAUPEL_DEPOSITION, ONLY: LIMA_GRAUPEL_DEPOSITION
USE MODE_LIMA_DROPLETS_RIMING_SNOW, ONLY: LIMA_DROPLETS_RIMING_SNOW
USE MODE_LIMA_RAIN_ACCR_SNOW, ONLY: LIMA_RAIN_ACCR_SNOW
USE MODE_LIMA_CONVERSION_MELTING_SNOW, ONLY: LIMA_CONVERSION_MELTING_SNOW
USE MODE_LIMA_RAIN_FREEZING, ONLY: LIMA_RAIN_FREEZING
USE MODE_LIMA_COLLISIONAL_ICE_BREAKUP, ONLY: LIMA_COLLISIONAL_ICE_BREAKUP
USE MODE_LIMA_RAINDROP_SHATTERING_FREEZING, ONLY: LIMA_RAINDROP_SHATTERING_FREEZING
USE MODE_LIMA_GRAUPEL, ONLY: LIMA_GRAUPEL
USE MODE_LIMA_HAIL_DEPOSITION, ONLY: LIMA_HAIL_DEPOSITION
USE MODE_LIMA_HAIL, ONLY: LIMA_HAIL
!
USE MODE_LIMA_BERGERON, ONLY: LIMA_BERGERON
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                 INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PEXNREF   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF  ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PPABST    ! Pressure
REAL, DIMENSION(:),   INTENT(IN)    :: PTHT      ! Potential temperature
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT      ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT      ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT      ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT      ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRST      ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT      !
REAL, DIMENSION(:),   INTENT(IN)    :: PRHT      ! Mixing ratios (kg/kg)
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT      !
REAL, DIMENSION(:),   INTENT(INOUT) :: PCRT      !
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT      ! 
REAL, DIMENSION(:),   INTENT(INOUT) :: PCST      ! 
REAL, DIMENSION(:),   INTENT(INOUT) :: PCGT      !
REAL, DIMENSION(:),   INTENT(INOUT) :: PCHT      ! Number concentrations (/kg)
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_HONC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_HONC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_HONC ! droplets homogeneous freezing (HONC) : rc, Nc, ri=-rc, Ni=-Nc, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_SELF ! self collection of droplets (SELF) : Nc
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_AUTO
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_AUTO
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_AUTO ! autoconversion of cloud droplets (AUTO) : rc, Nc, rr=-rc, Nr
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_ACCR
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_ACCR ! accretion of droplets by rain drops (ACCR) : rc, Nc, rr=-rr
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_SCBU ! self collectio break up of drops (SCBU) : Nr
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_EVAP
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_EVAP
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_EVAP ! evaporation of rain drops (EVAP) : rr, Nr, rv=-rr
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CNVI
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CNVI  ! conversion snow -> ice (CNVI) : ri, Ni, rs=-ri
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_DEPS  ! deposition of vapor on snow (DEPS) : rv=-rs, rs, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPI
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_DEPI  ! deposition of vapor on ice (DEPI) : rv=-ri, ri, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CNVS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CNVS  ! conversion ice -> snow (CNVS) : ri, Ni, rs=-ri
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_SSC   ! self collection of snow (SSC) : Ns
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_AGGS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_AGGS  ! aggregation of ice on snow (AGGS) : ri, Ni, rs=-ri
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_DEPG  ! deposition of vapor on graupel (DEPG) : rv=-rg, rg, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_BERFI
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_BERFI ! Bergeron (BERFI) : rc, ri=-rc, th
!
!++cb++
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_RIM
!REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_RIM
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_RIM
!REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_RIM
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_RIM
!REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_RIM   ! cloud droplet riming (RIM) : rc, Nc, rs, rg, th
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_RIMSS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_RIMSG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_RIMCG   ! cloud droplet riming (RIM) : rc, Nc, rs, rg, th
!--cb--
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_HMS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_HMS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_HMS   ! hallett mossop snow (HMS) : ri, Ni, rs
!
!++cb++
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_ACC
!REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_ACC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_ACC
!REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_ACC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_ACC
!REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_ACC   ! rain accretion on aggregates (ACC) : rr, Nr, rs, Ns, rg, Ng=-Ns, th
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_ACCSS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_ACCSG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_ACCRG ! rain accretion on aggregates (ACC) : rr, Nr, rs, rg, th
!--cb--
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_CMEL
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_CMEL  ! conversion-melting (CMEL) : rs, Ns, rg=-rs, Ng=-Ns
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CFRZ  ! rain freezing (CFRZ) : rr, Nr, ri, Ni, rg=-rr-ri, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CIBU
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CIBU  ! collisional ice break-up (CIBU) : ri, Ni, rs=-ri
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_RDSF
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_RDSF  ! rain drops freezing shattering (RDSF) : ri, Ni, rg=-ri
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CG_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RH_WETG ! wet growth of graupel (WETG) : rc, NC, rr, Nr, ri, Ni, rs, Ns, rg, Ng, rh, Nh=-Ng, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_DRYG ! dry growth of graupel (DRYG) : rc, Nc, rr, Nr, ri, Ni, rs, rg, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_HMG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_HMG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_HMG  ! hallett mossop graupel (HMG) : ri, Ni, rg
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CG_GMLT ! graupel melting (GMLT) : rr, Nr, rg=-rr, Ng, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RH_DEPH  ! deposition of vapor on hail (DEPH) : rv=-rh, rh, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CS_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CG_WETH
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RH_WETH ! wet growth of hail (WETH) : rc, NC, rr, Nr, ri, Ni, rs, Ns, rg, Ng, rh, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_COHG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CG_COHG ! conversion hail -> graupel (COHG) : rg, Ng, rh=-rg; Nh=-Ng
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_HMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_HMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_HMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CH_HMLT ! hail melting (HMLT) : rr, Nr, rh=-rr, Nh, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RV
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RR
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CR
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RS
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CS
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RG
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CG
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CH
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PEVAP3D
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCF1D
REAL, DIMENSION(:),   INTENT(IN)    :: PIF1D
REAL, DIMENSION(:),   INTENT(IN)    :: PPF1D
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLATHAM_IAGGS ! factor to account for the effect of Efield on IAGGS
!
!*       0.2   Declarations of local variables :
!
REAL,    DIMENSION(SIZE(PRCT))  :: ZT

REAL,    DIMENSION(SIZE(PRCT))  :: ZLBDC
REAL,    DIMENSION(SIZE(PRCT))  :: ZLBDC3
REAL,    DIMENSION(SIZE(PRCT))  :: ZLBDR
REAL,    DIMENSION(SIZE(PRCT))  :: ZLBDR3
REAL,    DIMENSION(SIZE(PRCT))  :: ZLBDI
REAL,    DIMENSION(SIZE(PRCT))  :: ZLBDS
REAL,    DIMENSION(SIZE(PRCT))  :: ZLBDS3
REAL,    DIMENSION(SIZE(PRCT))  :: ZLBDG
REAL,    DIMENSION(SIZE(PRCT))  :: ZLBDH

REAL,    DIMENSION(SIZE(PRCT))  :: ZAI
REAL,    DIMENSION(SIZE(PRCT))  :: ZKA
REAL,    DIMENSION(SIZE(PRCT))  :: ZDV
REAL,    DIMENSION(SIZE(PRCT))  :: ZCJ

REAL,    DIMENSION(SIZE(PRCT))  :: ZEPS
REAL,    DIMENSION(SIZE(PRCT))  :: ZEVSAT
REAL,    DIMENSION(SIZE(PRCT))  :: ZEISAT
REAL,    DIMENSION(SIZE(PRCT))  :: ZRVSAT
REAL,    DIMENSION(SIZE(PRCT))  :: ZRISAT
!
REAL,    DIMENSION(SIZE(PRCT))  :: ZSSI
REAL,    DIMENSION(SIZE(PRCT))  :: ZSSIW

REAL,    DIMENSION(SIZE(PRCT))  :: ZLV
REAL,    DIMENSION(SIZE(PRCT))  :: ZLS
REAL,    DIMENSION(SIZE(PRCT))  :: ZLVFACT
REAL,    DIMENSION(SIZE(PRCT))  :: ZLSFACT
!
REAL,    DIMENSION(SIZE(PRCT))  :: ZW
!
REAL,    DIMENSION(SIZE(PRCT))  :: ZCF1D
REAL,    DIMENSION(SIZE(PRCT))  :: ZIF1D
REAL,    DIMENSION(SIZE(PRCT))  :: ZPF1D
!
REAL,    DIMENSION(SIZE(PRCT))  :: ZRCT
REAL,    DIMENSION(SIZE(PRCT))  :: ZRRT
REAL,    DIMENSION(SIZE(PRCT))  :: ZRIT
REAL,    DIMENSION(SIZE(PRCT))  :: ZRST
REAL,    DIMENSION(SIZE(PRCT))  :: ZRGT
REAL,    DIMENSION(SIZE(PRCT))  :: ZRHT
!-------------------------------------------------------------------------------
! Pre-compute quantities
!
! Prevent fractions to reach 0 (divide by 0)
!
ZCF1D(:) = MAX(PCF1D(:),0.01)
ZIF1D(:) = MAX(PIF1D(:),0.01)
ZPF1D(:) = MAX(PPF1D(:),0.01)
!
!
!
WHERE (PRCT(:).GE.XRTMIN(2))
   ZRCT(:)=PRCT(:)
ELSEWHERE
   ZRCT(:)=0.
END WHERE
WHERE (PRRT(:).GE.XRTMIN(3))
   ZRRT(:)=PRRT(:)
ELSEWHERE
   ZRRT(:)=0.
END WHERE
WHERE (PRIT(:).GE.XRTMIN(4))
   ZRIT(:)=PRIT(:)
ELSEWHERE
   ZRIT(:)=0.
END WHERE
WHERE (PRST(:).GE.XRTMIN(5))
   ZRST(:)=PRST(:)
ELSEWHERE
   ZRST(:)=0.
END WHERE
WHERE (PRGT(:).GE.XRTMIN(6))
   ZRGT(:)=PRGT(:)
ELSEWHERE
   ZRGT(:)=0.
END WHERE
WHERE (PRHT(:).GE.XRTMIN(7))
   ZRHT(:)=PRHT(:)
ELSEWHERE
   ZRHT(:)=0.
END WHERE
!
! Is it necessary to compute the following quantities
! accounting for subgrig cloud fraction ?
! lambda does not depend on cloud fraction for 2-m species
! lambda depends on CF for 1-m species ?
!
!
! Is it necessary to change water vapour in cloudy / non cloudy parts ?
!
!
WHERE (LDCOMPUTE(:))
   ZT(:) = PTHT(:) * ( PPABST(:)/XP00 ) ** (XRD/XCPD)
!
   ZW(:) = PEXNREF(:)*( XCPD &
                               +XCPV*PRVT(:) &
                               +XCL*(ZRCT(:)+ZRRT(:)) &
                               +XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)+ZRHT(:)) )
!
   ZLV(:) = XLVTT + (XCPV-XCL)*(ZT(:)-XTT)
   ZLVFACT(:) = ZLV(:)/ZW(:)               ! L_v/(Pi_ref*C_ph)
   ZLS(:) = XLSTT + (XCPV-XCI)*(ZT(:)-XTT)
   ZLSFACT(:) = ZLS(:)/ZW(:)               ! L_s/(Pi_ref*C_ph)
!
   ZEVSAT(:)  = EXP( XALPW - XBETAW/ZT(:) - XGAMW*ALOG(ZT(:) ) )
   ZEISAT(:)  = EXP( XALPI - XBETAI/ZT(:) - XGAMI*ALOG(ZT(:) ) )
   !
   ZEPS= XMV / XMD
   ZRVSAT(:) = ZEPS * ZEVSAT(:) / (PPABST(:) - ZEVSAT(:))
   ZRISAT(:) = ZEPS * ZEISAT(:) / (PPABST(:) - ZEISAT(:))
   !
   ZSSI(:)  = PRVT(:)/ZRISAT(:) - 1.0             ! Si  =  rv/rsi - 1
   ZSSIW(:) = ZRVSAT(:)/ZRISAT(:) - 1.0 ! Siw = rsw/rsi - 1
!
   ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZT(:) - XTT )
!
   ZDV(:) = 0.211E-4 * (ZT(:)/XTT)**1.94 * (XP00/PPABST(:))
!
   ZAI(:) =   ( XLSTT + (XCPV-XCI)*(ZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZT(:)**2) &
                + ( XRV*ZT(:) ) / (ZDV(:)*ZEISAT(:))
!
   ZCJ(:) = XSCFAC * PRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZT(:)-XTT) )
!
END WHERE
!
! Cloud droplets : same formula for 1 and 2 moments, but using real or fixed Nc value
ZLBDC(:)  = 1.E10
ZLBDC3(:) = 1.E30
WHERE (ZRCT(:)>XRTMIN(2) .AND. PCCT(:)>XCTMIN(2) .AND. LDCOMPUTE(:))
   ZLBDC3(:) = XLBC*PCCT(:) / ZRCT(:)
   ZLBDC(:)  = ZLBDC3(:)**XLBEXC
END WHERE
!
! Rain drops
ZLBDR(:)  = 1.E10
ZLBDR3(:) = 1.E30
IF (NMOM_R.EQ.1) THEN
   WHERE (ZRRT(:)>XRTMIN(3) .AND. LDCOMPUTE(:) )
      ZLBDR(:) = XLBR*(PRHODREF(:)*ZRRT(:) )**XLBEXR
      ZLBDR3(:) = ZLBDR(:)**3.
   END WHERE
   PCRT(:) = XCCR * ZLBDR(:)**XCXR / PRHODREF(:)
ELSE
   WHERE (ZRRT(:)>XRTMIN(3) .AND. PCRT(:)>XCTMIN(3) .AND. LDCOMPUTE(:))
      ZLBDR3(:) = XLBR*PCRT(:) / ZRRT(:)
      ZLBDR(:)  = ZLBDR3(:)**XLBEXR
   END WHERE
END IF
!
! Pristine ice : same formula for 1 and 2 moments, using real or diagnosed Ni
ZLBDI(:)  = 1.E10
WHERE (ZRIT(:)>XRTMIN(4) .AND. PCIT(:)>XCTMIN(4) .AND. LDCOMPUTE(:))
   ZLBDI(:) = ( XLBI*PCIT(:) / ZRIT(:) )**XLBEXI
END WHERE
!
! Snow : additional option for LSNOW_T if NMOM_S=1
ZLBDS(:)  = 1.E10
IF (NMOM_S.EQ.1) THEN
   IF (LSNOW_T) THEN
      WHERE (ZRST(:)>XRTMIN(5) .AND. LDCOMPUTE(:) )
         WHERE(ZT(:)>263.15)
            ZLBDS(:) = MAX(MIN(XLBDAS_MAX, 10**(14.554-0.0423*ZT(:))),XLBDAS_MIN)
         END WHERE
         WHERE(ZT(:)<=263.15)
            ZLBDS(:) = MAX(MIN(XLBDAS_MAX, 10**(6.226-0.0106*ZT(:))),XLBDAS_MIN)
         END WHERE
      END WHERE
      ZLBDS(:) =  ZLBDS(:) * XTRANS_MP_GAMMAS
      PCST(:) = XNS * ZRST(:) * ZLBDS(:)**XBS
   ELSE
      WHERE (ZRST(:)>XRTMIN(5) .AND. LDCOMPUTE(:) )
         ZLBDS(:) = XLBS*( PRHODREF(:)*ZRST(:) )**XLBEXS
      END WHERE
      PCST(:) = XCCS * ZLBDS(:)**XCXS / PRHODREF(:)
   END IF
ELSE
   WHERE (ZRST(:)>XRTMIN(5) .AND. PCST(:)>XCTMIN(5) .AND. LDCOMPUTE(:))
      ZLBDS(:) = (XLBS*PCST(:)/ZRST(:))**XLBEXS
   END WHERE
END IF
ZLBDS3(:) = ZLBDS(:)**3.
!
! Graupel
ZLBDG(:)  = 1.E10
IF (NMOM_G.EQ.1) THEN
   WHERE (ZRGT(:)>XRTMIN(6) .AND. LDCOMPUTE(:) )
      ZLBDG(:) = XLBG*( PRHODREF(:)*ZRGT(:) )**XLBEXG
   END WHERE
      PCGT(:) = XCCG * ZLBDG(:)**XCXG / PRHODREF(:)
ELSE
   WHERE (ZRGT(:)>XRTMIN(6) .AND. PCGT(:)>XCTMIN(6) .AND. LDCOMPUTE(:))
      ZLBDG(:) = (XLBG*PCGT(:)/ZRGT(:))**XLBEXG
   END WHERE
END IF
!
! Hail
ZLBDH(:)  = 1.E10
IF (NMOM_H.EQ.1) THEN
   WHERE (ZRHT(:)>XRTMIN(7) .AND. LDCOMPUTE(:) )
      ZLBDH(:) = XLBH*( PRHODREF(:)*ZRHT(:) )**XLBEXH
      PCHT(:) = XCCH * ZLBDH(:)**XCXH / PRHODREF(:)
   END WHERE
ELSE
   WHERE (ZRHT(:)>XRTMIN(7) .AND. PCHT(:)>XCTMIN(7) .AND. LDCOMPUTE(:))
      ZLBDH(:) = (XLBH*PCHT(:)/ZRHT(:))**XLBEXH
   END WHERE
END IF
!
!-------------------------------------------------------------------------------
! Call microphysical processes   
!
IF (NMOM_C.GE.1 .AND. NMOM_I.GE.1) THEN
   CALL LIMA_DROPLETS_HOM_FREEZING (PTSTEP, LDCOMPUTE,                 & ! independent from CF,IF,PF
                                    ZT, ZLVFACT, ZLSFACT,              &
                                    ZRCT, PCCT, ZLBDC,                 &
                                    P_TH_HONC, P_RC_HONC, P_CC_HONC    )
   PA_RC(:) = PA_RC(:) + P_RC_HONC(:)
   IF (NMOM_C.GE.2) PA_CC(:) = PA_CC(:) + P_CC_HONC(:)
   PA_RI(:) = PA_RI(:) - P_RC_HONC(:)
   IF (NMOM_I.GE.2) PA_CI(:) = PA_CI(:) - P_CC_HONC(:) 
   P_TH_ACC(:) = - P_RC_HONC(:) * (ZLSFACT(:)-ZLVFACT(:))
   PA_TH(:) = PA_TH(:) + P_TH_HONC(:)
END IF
!
IF ((.NOT. LKHKO) .AND. NMOM_C.GE.2) THEN
   CALL LIMA_DROPLETS_SELF_COLLECTION (LDCOMPUTE,          & ! depends on CF
                                       PRHODREF,           &
                                       PCCT/ZCF1D, ZLBDC3, &
                                       P_CC_SELF           )
   P_CC_SELF(:) = P_CC_SELF(:) * ZCF1D(:)
   PA_CC(:) = PA_CC(:) + P_CC_SELF(:)
END IF
!
IF (NMOM_C.GE.1 .AND. NMOM_R.GE.1) THEN
   CALL LIMA_DROPLETS_AUTOCONVERSION (LDCOMPUTE,                             & ! depends on CF
                                      PRHODREF,                              &
                                      ZRCT/ZCF1D, PCCT/ZCF1D, ZLBDC, ZLBDR,  &
                                      P_RC_AUTO, P_CC_AUTO, P_CR_AUTO        )
   P_RC_AUTO(:) = P_RC_AUTO(:) * ZCF1D(:)
   P_CC_AUTO(:) = P_CC_AUTO(:) * ZCF1D(:)
   P_CR_AUTO(:) = P_CR_AUTO(:) * ZCF1D(:)
   !
   PA_RC(:) = PA_RC(:) + P_RC_AUTO(:)
   IF (NMOM_C.GE.2) PA_CC(:) = PA_CC(:) + P_CC_AUTO(:)
   PA_RR(:) = PA_RR(:) - P_RC_AUTO(:)
   IF (NMOM_R.GE.2) PA_CR(:) = PA_CR(:) + P_CR_AUTO(:)
END IF
!
IF (NMOM_C.GE.1 .AND. NMOM_R.GE.1) THEN
   CALL LIMA_DROPLETS_ACCRETION (LDCOMPUTE,                                     & ! depends on CF, PF
                                 PRHODREF,                                      &
                                 ZRCT/ZCF1D, ZRRT/ZPF1D, PCCT/ZCF1D, PCRT/ZPF1D,&
                                 ZLBDC, ZLBDC3, ZLBDR, ZLBDR3,                  &
                                 P_RC_ACCR, P_CC_ACCR                           )
   !
   P_CC_ACCR(:) = P_CC_ACCR(:) * ZCF1D(:)
   P_RC_ACCR(:) = P_RC_ACCR(:) * ZCF1D(:)
   !
   PA_RC(:) = PA_RC(:) + P_RC_ACCR(:)
   IF (NMOM_C.GE.2) PA_CC(:) = PA_CC(:) + P_CC_ACCR(:)
   PA_RR(:) = PA_RR(:) - P_RC_ACCR(:)
END IF
!
IF ((.NOT. LKHKO) .AND. NMOM_R.GE.2) THEN 
   CALL LIMA_DROPS_SELF_COLLECTION (LDCOMPUTE,           & ! depends on PF
                                    PRHODREF,            &
                                    PCRT/ZPF1D(:), ZLBDR, ZLBDR3, &
                                    P_CR_SCBU            )
   !
   P_CR_SCBU(:) = P_CR_SCBU(:) * ZPF1D(:)
   !
   PA_CR(:) = PA_CR(:) + P_CR_SCBU(:)
END IF
!
IF (NMOM_R.GE.1) THEN
   CALL LIMA_RAIN_EVAPORATION (PTSTEP, LDCOMPUTE,                               & ! depends on PF > CF 
                               PRHODREF, ZT, ZLV, ZLVFACT, ZEVSAT, ZRVSAT,      &
                               PRVT, ZRCT/ZPF1D, ZRRT/ZPF1D, PCRT/ZPF1D, ZLBDR, &
                               P_TH_EVAP, P_RR_EVAP, P_CR_EVAP,                 &
                               PEVAP3D                                          )
   P_RR_EVAP(:) = P_RR_EVAP(:) * MAX((ZPF1D(:) - ZCF1D(:)),0.)
   P_CR_EVAP(:) = P_RR_EVAP(:) * MAX((ZPF1D(:) - ZCF1D(:)),0.)
   P_TH_EVAP(:) = P_RR_EVAP(:) * ZLVFACT(:)
   PEVAP3D(:) = - P_RR_EVAP(:)
   !
   PA_TH(:) = PA_TH(:) + P_TH_EVAP(:)
   PA_RV(:) = PA_RV(:) - P_RR_EVAP(:)
   PA_RR(:) = PA_RR(:) + P_RR_EVAP(:)
   IF (NMOM_R.GE.2) PA_CR(:) = PA_CR(:) + P_CR_EVAP(:)
END IF
!
IF (NMOM_I.GE.1) THEN
   !
   ! Includes vapour deposition on ice, ice -> snow conversion
   !
   CALL LIMA_ICE_DEPOSITION (PTSTEP, LDCOMPUTE,                 & ! depends on IF, PF
                             PRHODREF, ZT, ZSSI, ZAI, ZCJ, ZLSFACT, &
                             ZRIT/ZIF1D, PCIT/ZIF1D, ZLBDI,     &
                             P_TH_DEPI, P_RI_DEPI,              &
                             P_RI_CNVS, P_CI_CNVS               )
   !
   P_RI_DEPI(:) = P_RI_DEPI(:) * ZIF1D(:)
   P_RI_CNVS(:) = P_RI_CNVS(:) * ZIF1D(:)
   P_CI_CNVS(:) = P_CI_CNVS(:) * ZIF1D(:)
   P_TH_DEPI(:) = P_RI_DEPI(:) * ZLSFACT(:)
   !
   PA_TH(:) = PA_TH(:) + P_TH_DEPI(:)
   PA_RV(:) = PA_RV(:) - P_RI_DEPI(:) 
   PA_RI(:) = PA_RI(:) + P_RI_DEPI(:) + P_RI_CNVS(:)
   IF (NMOM_I.GE.2) PA_CI(:) = PA_CI(:)                + P_CI_CNVS(:)
   PA_RS(:) = PA_RS(:)                - P_RI_CNVS(:)
   IF (NMOM_S.GE.2) PA_CS(:) = PA_CS(:)                - P_CI_CNVS(:)

END IF
!
IF (NMOM_S.GE.1) THEN
   !
   ! Includes vapour deposition on snow, snow -> ice conversion
   !
   CALL LIMA_SNOW_DEPOSITION (LDCOMPUTE,                         & ! depends on IF, PF
                              PRHODREF, ZSSI, ZAI, ZCJ, ZLSFACT, &
                              ZRST/ZPF1D, PCST/ZPF1D, ZLBDS,     &
                              P_RI_CNVI, P_CI_CNVI,              &
                              P_TH_DEPS, P_RS_DEPS               )
   !
   P_RI_CNVI(:) = P_RI_CNVI(:) * ZPF1D(:)
   P_CI_CNVI(:) = P_CI_CNVI(:) * ZPF1D(:)
   P_RS_DEPS(:) = P_RS_DEPS(:) * ZPF1D(:)
   P_TH_DEPS(:) = P_RS_DEPS(:) * ZLSFACT(:)
   !
   PA_RI(:) = PA_RI(:) + P_RI_CNVI(:)
   IF (NMOM_I.GE.2) PA_CI(:) = PA_CI(:) + P_CI_CNVI(:)
   PA_RS(:) = PA_RS(:) - P_RI_CNVI(:) + P_RS_DEPS(:) 
   IF (NMOM_S.GE.2) PA_CS(:) = PA_CS(:) - P_CI_CNVI(:)
   PA_TH(:) = PA_TH(:)                + P_TH_DEPS(:)
   PA_RV(:) = PA_RV(:)                - P_RS_DEPS(:) 

END IF
!
IF (NMOM_S.GE.2) THEN 
   CALL LIMA_SNOW_SELF_COLLECTION (LDCOMPUTE,           & ! depends on PF
                                   PRHODREF,            &
                                   ZRST(:)/ZPF1D(:), PCST/ZPF1D(:), ZLBDS, ZLBDS3, &
                                   P_CS_SSC             )
   !
   P_CS_SSC(:) = P_CS_SSC(:) * ZPF1D(:)
   !
   PA_CS(:) = PA_CS(:) + P_CS_SSC(:)
END IF
!
! Lambda_s limited for collection processes to prevent too high concentrations
! must be changed or removed if C and x modified
!
!ZLBDS(:) = MIN( XLBDAS_MAX, ZLBDS(:))
!
!
IF (NMOM_I.GE.1 .AND. NMOM_S.GE.1) THEN
   CALL LIMA_ICE_AGGREGATION_SNOW (LDCOMPUTE,                                                    & ! depends on IF, PF
                                   ZT, PRHODREF,                                                 &
                                   ZRIT/ZIF1D, ZRST/ZPF1D, PCIT/ZIF1D, PCST/ZPF1D, ZLBDI, ZLBDS, &
                                   PLATHAM_IAGGS,                                                &
                                   P_RI_AGGS, P_CI_AGGS                                          )
   P_CI_AGGS(:) = P_CI_AGGS(:) * ZIF1D(:)
   P_RI_AGGS(:) = P_RI_AGGS(:) * ZIF1D(:)
   !
   PA_RI(:) = PA_RI(:) + P_RI_AGGS(:)
   IF (NMOM_I.GE.2) PA_CI(:) = PA_CI(:) + P_CI_AGGS(:)
   PA_RS(:) = PA_RS(:) - P_RI_AGGS(:)
END IF
!
IF (NMOM_G.GE.1) THEN
   CALL LIMA_GRAUPEL_DEPOSITION (LDCOMPUTE, PRHODREF,                                    & ! depends on PF ?
                                 ZRGT/ZPF1D, PCGT/ZPF1D, ZSSI, ZLBDG, ZAI, ZCJ, ZLSFACT, &
                                 P_TH_DEPG, P_RG_DEPG                                    )
   P_RG_DEPG(:) = P_RG_DEPG(:) * ZPF1D(:)
   P_TH_DEPG(:) = P_RG_DEPG(:) * ZLSFACT(:)
   !
   PA_RV(:) = PA_RV(:) - P_RG_DEPG(:)
   PA_RG(:) = PA_RG(:) + P_RG_DEPG(:)
   PA_TH(:) = PA_TH(:) + P_TH_DEPG(:)
END IF
!
IF (NMOM_C.GE.1 .AND. NMOM_I.EQ.1) THEN
   CALL LIMA_BERGERON (LDCOMPUTE,                                 & ! depends on CF, IF
                       ZRCT/ZCF1D, ZRIT/ZIF1D, PCIT/ZIF1D, ZLBDI, &
                       ZSSIW, ZAI, ZCJ, ZLVFACT, ZLSFACT,         &
                       P_TH_BERFI, P_RC_BERFI                     )
   P_TH_BERFI(:) = P_TH_BERFI(:) * MIN(ZCF1D,ZIF1D)
   P_RC_BERFI(:) = P_RC_BERFI(:) * MIN(ZCF1D,ZIF1D)
!
   PA_RC(:) = PA_RC(:) + P_RC_BERFI(:)
   PA_RI(:) = PA_RI(:) - P_RC_BERFI(:)
   PA_TH(:) = PA_TH(:) + P_TH_BERFI(:)
END IF
!
IF (NMOM_C.GE.1 .AND. NMOM_S.GE.1) THEN
     !
     ! Graupel production as tendency (or should be tendency + instant to stick to the previous version ?)
     ! Includes the Hallett Mossop process for riming of droplets by snow (HMS)
     !
!++cb++
   CALL LIMA_DROPLETS_RIMING_SNOW (PTSTEP, LDCOMPUTE,                                & ! depends on CF
                                   PRHODREF, ZT,                                     &
                                   ZRCT/ZCF1D, PCCT/ZCF1D, ZRST/ZPF1D, PCST/ZPF1D, ZLBDC, ZLBDS, ZLVFACT, ZLSFACT, &
!                                   P_TH_RIM, P_RC_RIM, P_CC_RIM, P_RS_RIM, P_CS_RIM, P_RG_RIM, &
                                   P_TH_RIM, P_CC_RIM, P_CS_RIM, P_RC_RIMSS, P_RC_RIMSG, P_RS_RIMCG, &
                                   P_RI_HMS, P_CI_HMS, P_RS_HMS                      )
!   P_RC_RIM(:) = P_RC_RIM(:) * ZCF1D(:)
   P_CC_RIM(:) = P_CC_RIM(:) * ZCF1D(:)
!   P_RS_RIM(:) = P_RS_RIM(:) * ZCF1D(:)
   P_CS_RIM(:) = P_CS_RIM(:) * ZCF1D(:)
!   P_RG_RIM(:) = P_RG_RIM(:) * ZCF1D(:)
   P_RC_RIMSS(:) = P_RC_RIMSS(:) * ZCF1D(:)
   P_RC_RIMSG(:) = P_RC_RIMSG(:) * ZCF1D(:)
   P_RS_RIMCG(:) = P_RS_RIMCG(:) * ZCF1D(:)
!   P_TH_RIM(:) = - P_RC_RIM(:) * (ZLSFACT(:)-ZLVFACT(:))
   P_TH_RIM(:) = - (P_RC_RIMSS(:) + P_RC_RIMSG(:)) * (ZLSFACT(:)-ZLVFACT(:))
   P_RI_HMS(:) = P_RI_HMS(:) * ZCF1D(:)
   P_CI_HMS(:) = P_CI_HMS(:) * ZCF1D(:)
   P_RS_HMS(:) = P_RS_HMS(:) * ZCF1D(:)
   !
!   PA_RC(:) = PA_RC(:) + P_RC_RIM(:)
   PA_RC(:) = PA_RC(:) + P_RC_RIMSS(:) + P_RC_RIMSG(:) ! RCRIMSS < 0 and RCRIMSG < 0 (both loss for rc)
   IF (NMOM_C.GE.2) PA_CC(:) = PA_CC(:) + P_CC_RIM(:) 
   PA_RI(:) = PA_RI(:)               + P_RI_HMS(:)
   IF (NMOM_I.GE.2) PA_CI(:) = PA_CI(:)               + P_CI_HMS(:)
!   PA_RS(:) = PA_RS(:) + P_RS_RIM(:) + P_RS_HMS(:)
   PA_RS(:) = PA_RS(:) - P_RC_RIMSS(:) - P_RS_RIMCG(:) + P_RS_HMS(:) ! RCRIMSS < 0 (gain for rs), RSRIMCG > 0 (loss for rs)
   IF (NMOM_S.GE.2) PA_CS(:) = PA_CS(:) + P_CS_RIM(:)
!   PA_RG(:) = PA_RG(:) + P_RG_RIM(:)
   PA_RG(:) = PA_RG(:) - P_RC_RIMSG(:) + P_RS_RIMCG(:) ! RCRIMSG < 0 (gain for rg), RSRIMCG > 0 (gain for rg)
   IF (NMOM_G.GE.2) PA_CG(:) = PA_CG(:) - P_CS_RIM(:)
   PA_TH(:) = PA_TH(:) + P_TH_RIM(:)
!--cb--
END IF
!
IF (NMOM_R.GE.1 .AND. NMOM_S.GE.1) THEN
!++cb++
   CALL LIMA_RAIN_ACCR_SNOW (PTSTEP, LDCOMPUTE,                                & ! depends on PF
                             PRHODREF, ZT,                                     &
                             ZRRT/ZPF1D, PCRT/ZPF1D, ZRST/ZPF1D, PCST/ZPF1D, ZLBDR, ZLBDS, ZLVFACT, ZLSFACT, &
!                             P_TH_ACC, P_RR_ACC, P_CR_ACC, P_RS_ACC, P_CS_ACC, P_RG_ACC )
                             P_TH_ACC, P_CR_ACC, P_CS_ACC, P_RR_ACCSS, P_RR_ACCSG, P_RS_ACCRG )
!   P_RR_ACC(:) = P_RR_ACC(:) * ZPF1D(:)
   P_CR_ACC(:) = P_CR_ACC(:) * ZPF1D(:)
!   P_RS_ACC(:) = P_RS_ACC(:) * ZPF1D(:)
   P_CS_ACC(:) = P_CS_ACC(:) * ZPF1D(:)
!   P_RG_ACC(:) = P_RG_ACC(:) * ZPF1D(:)
!   P_TH_ACC(:) = - P_RR_ACC(:) * (ZLSFACT(:)-ZLVFACT(:))
   P_RR_ACCSS(:) = P_RR_ACCSS(:) * ZPF1D(:)
   P_RR_ACCSG(:) = P_RR_ACCSG(:) * ZPF1D(:)
   P_RS_ACCRG(:) = P_RS_ACCRG(:) * ZPF1D(:)
   P_TH_ACC(:)   = (P_RR_ACCSS(:) + P_RR_ACCSG(:)) * (ZLSFACT(:)-ZLVFACT(:))
   !
!   PA_RR(:) = PA_RR(:) + P_RR_ACC(:)
   PA_RR(:) = PA_RR(:) - P_RR_ACCSS(:) - P_RR_ACCSG(:)
   IF (NMOM_R.GE.2) PA_CR(:) = PA_CR(:) + P_CR_ACC(:)
!   PA_RS(:) = PA_RS(:) + P_RS_ACC(:)
   PA_RS(:) = PA_RS(:) + P_RR_ACCSS(:) - P_RS_ACCRG(:)
   IF (NMOM_S.GE.2) PA_CS(:) = PA_CS(:) + P_CS_ACC(:)
!   PA_RG(:) = PA_RG(:) + P_RG_ACC(:)
   PA_RG(:) = PA_RG(:) + P_RR_ACCSG(:) + P_RS_ACCRG(:)
   IF (NMOM_G.GE.2) PA_CG(:) = PA_CG(:) - P_CS_ACC(:)
   PA_TH(:) = PA_TH(:) + P_TH_ACC(:)
!--cb--
END IF
!
IF (NMOM_S.GE.1) THEN
   !
   ! Conversion melting of snow should account for collected droplets and drops where T>0C, but does not !
   ! Some thermodynamical computations inside, to externalize ?
   !
   CALL LIMA_CONVERSION_MELTING_SNOW (LDCOMPUTE,                           & ! depends on PF
                                      PRHODREF, PPABST, ZT, ZKA, ZDV, ZCJ, &
                                      PRVT, ZRST/ZPF1D, PCST/ZPF1D, ZLBDS, &
                                      P_RS_CMEL, P_CS_CMEL                 )
   P_RS_CMEL(:) = P_RS_CMEL(:) * ZPF1D(:)
   P_CS_CMEL(:) = P_CS_CMEL(:) * ZPF1D(:)
   !
   PA_RS(:) = PA_RS(:) + P_RS_CMEL(:)
   IF (NMOM_S.GE.2) PA_CS(:) = PA_CS(:) + P_CS_CMEL(:)
   PA_RG(:) = PA_RG(:) - P_RS_CMEL(:)
   IF (NMOM_G.GE.2) PA_CG(:) = PA_CG(:) - P_CS_CMEL(:)

END IF
!
IF (NMOM_R.GE.1) THEN
   CALL LIMA_RAIN_FREEZING (LDCOMPUTE,                                             & ! depends on PF, IF
                            PRHODREF, ZT, ZLVFACT, ZLSFACT,                        &
                            ZRRT/ZPF1D, PCRT/ZPF1D, ZRIT/ZIF1D, PCIT/ZIF1D, ZLBDR, &
                            P_TH_CFRZ, P_RR_CFRZ, P_CR_CFRZ, P_RI_CFRZ, P_CI_CFRZ  )
   P_RR_CFRZ(:) = P_RR_CFRZ(:) * ZIF1D(:)
   P_CR_CFRZ(:) = P_CR_CFRZ(:) * ZIF1D(:)
   P_RI_CFRZ(:) = P_RI_CFRZ(:) * ZIF1D(:)
   P_CI_CFRZ(:) = P_CI_CFRZ(:) * ZIF1D(:)
   P_TH_CFRZ(:) = - P_RR_CFRZ(:) * (ZLSFACT(:)-ZLVFACT(:))
!
   PA_TH(:) = PA_TH(:) + P_TH_CFRZ(:)
   PA_RR(:) = PA_RR(:) + P_RR_CFRZ(:)
   IF (NMOM_R.GE.2) PA_CR(:) = PA_CR(:) + P_CR_CFRZ(:)
   PA_RI(:) = PA_RI(:) + P_RI_CFRZ(:)
   IF (NMOM_I.GE.2) PA_CI(:) = PA_CI(:) + P_CI_CFRZ(:)
   PA_RG(:) = PA_RG(:) - P_RR_CFRZ(:) - P_RI_CFRZ(:)
   IF (NMOM_G.GE.2) PA_CG(:) = PA_CG(:) - P_CR_CFRZ(:)

END IF
!
IF (NMOM_S.GE.1 .AND. NMOM_G.GE.1 .AND. LCIBU) THEN
   !
   ! Conversion melting of snow should account for collected droplets and drops where T>0C, but does not !
   ! Some thermodynamical computations inside, to externalize ?
   !
   CALL LIMA_COLLISIONAL_ICE_BREAKUP (LDCOMPUTE,                                      & ! depends on PF (IF for fragments size)
                                      PRHODREF,                                       &
                                      ZRIT/ZIF1D, ZRST/ZPF1D, ZRGT/ZPF1D, PCIT/ZIF1D, PCST/ZPF1D, PCST/ZPF1D, &
                                      ZLBDS, ZLBDG,                                   &
                                      P_RI_CIBU, P_CI_CIBU                            )
   P_RI_CIBU(:) = P_RI_CIBU(:) * ZPF1D(:)
   P_CI_CIBU(:) = P_CI_CIBU(:) * ZPF1D(:)
   !
   PA_RI(:) = PA_RI(:) + P_RI_CIBU(:)
   IF (NMOM_I.GE.2) PA_CI(:) = PA_CI(:) + P_CI_CIBU(:)
   PA_RS(:) = PA_RS(:) - P_RI_CIBU(:)

END IF
!
IF (NMOM_R.GE.1 .AND. NMOM_I.GE.1 .AND. LRDSF) THEN
   !
   ! Conversion melting of snow should account for collected droplets and drops where T>0C, but does not !
   ! Some thermodynamical computations inside, to externalize ?
   !
   CALL LIMA_RAINDROP_SHATTERING_FREEZING (LDCOMPUTE,                                      & ! depends on PF, IF
                                           PRHODREF,                                       &
                                           ZRRT/ZPF1D, PCRT/ZPF1D, ZRIT/ZIF1D, PCIT/ZIF1D, ZRGT/ZPF1D, &
                                           ZLBDR,                                          &
                                           P_RI_RDSF, P_CI_RDSF                            )
   P_RI_RDSF(:) = P_RI_RDSF(:) * ZIF1D(:)
   P_CI_RDSF(:) = P_CI_RDSF(:) * ZIF1D(:)
   !
   PA_RI(:) = PA_RI(:) + P_RI_RDSF(:)
   IF (NMOM_I.GE.2) PA_CI(:) = PA_CI(:) + P_CI_RDSF(:)
   PA_RG(:) = PA_RG(:) - P_RI_RDSF(:)

END IF
!
IF (NMOM_G.GE.1) THEN
     !
     ! Melting of graupel should account for collected droplets and drops where T>0C, but does not !
     ! Collection and water shedding should also happen where T>0C, but do not !
     ! Hail production as tendency (should be instant to stick to the previous version ?)
     ! Includes Hallett-Mossop  process for riming of droplets by graupel (HMG)
     ! Some thermodynamical computations inside, to externalize ?
     !
   CALL LIMA_GRAUPEL (PTSTEP, LDCOMPUTE,                                     & ! depends on PF, CF, IF
                      PRHODREF, PPABST, ZT, ZKA, ZDV, ZCJ,                   &
                      PRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT,                    &
                      PCCT, PCRT, PCIT, PCST, PCGT,                          &
                      ZLBDC, ZLBDR, ZLBDS, ZLBDG,                            &
                      ZLVFACT, ZLSFACT,                                      &
                      P_TH_WETG, P_RC_WETG, P_CC_WETG, P_RR_WETG, P_CR_WETG, &
                      P_RI_WETG, P_CI_WETG, P_RS_WETG, P_CS_WETG, P_RG_WETG, P_CG_WETG, P_RH_WETG, &
                      P_TH_DRYG, P_RC_DRYG, P_CC_DRYG, P_RR_DRYG, P_CR_DRYG, &
                      P_RI_DRYG, P_CI_DRYG, P_RS_DRYG, P_CS_DRYG, P_RG_DRYG, &
                      P_RI_HMG, P_CI_HMG, P_RG_HMG,                          &
                      P_TH_GMLT, P_RR_GMLT, P_CR_GMLT, P_CG_GMLT,            &
                      PA_TH, PA_RC, PA_CC, PA_RR, PA_CR,                     &
                      PA_RI, PA_CI, PA_RS, PA_CS, PA_RG, PA_CG, PA_RH, PA_CH )
END IF
!
IF (NMOM_H.GE.1) THEN
   CALL LIMA_HAIL_DEPOSITION (LDCOMPUTE, PRHODREF,                                    & ! depends on PF ?
                              ZRHT/ZPF1D, PCHT/ZPF1D, ZSSI, ZLBDH, ZAI, ZCJ, ZLSFACT, &
                              P_TH_DEPH, P_RH_DEPH                                    )
   P_RH_DEPH(:) = P_RH_DEPH(:) * ZPF1D(:)
   P_TH_DEPH(:) = P_RH_DEPH(:) * ZLSFACT(:)
   !
   PA_RV(:) = PA_RV(:) - P_RH_DEPH(:)
   PA_RH(:) = PA_RH(:) + P_RH_DEPH(:)
   PA_TH(:) = PA_TH(:) + P_TH_DEPH(:)
!     CALL LIMA_HAIL_GROWTH   LIMA_HAIL_CONVERSION   LIMA_HAIL_MELTING
   CALL LIMA_HAIL (PTSTEP, LDCOMPUTE,                                     & ! depends on PF, CF, IF
                   PRHODREF, PPABST, ZT, ZKA, ZDV, ZCJ,                   &
                   PRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT, ZRHT,              &
                   PCCT, PCRT, PCIT, PCST, PCGT, PCHT,                    &
                   ZLBDC, ZLBDR, ZLBDS, ZLBDG, ZLBDH,                     &
                   ZLVFACT, ZLSFACT,                                      &
                   P_TH_WETH, P_RC_WETH, P_CC_WETH, P_RR_WETH, P_CR_WETH, &
                   P_RI_WETH, P_CI_WETH, P_RS_WETH, P_CS_WETH, P_RG_WETH, P_CG_WETH, P_RH_WETH, &
                   P_RG_COHG, P_CG_COHG,                                  &
                   P_TH_HMLT, P_RR_HMLT, P_CR_HMLT, P_CH_HMLT,            &
                   PA_TH, PA_RC, PA_CC, PA_RR, PA_CR,                     &
                   PA_RI, PA_CI, PA_RS, PA_CS, PA_RG, PA_CG, PA_RH, PA_CH )
END IF
   !  
END SUBROUTINE LIMA_TENDENCIES
END MODULE MODE_LIMA_TENDENCIES
