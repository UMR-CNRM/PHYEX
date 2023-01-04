!###############################
MODULE MODI_LIMA_TENDENCIES
!###############################
  INTERFACE
     SUBROUTINE LIMA_TENDENCIES (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,                &
                                 PEXNREF, PRHODREF, PPABST, PTHT,                       &
                                 PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT,              &
                                 PCCT, PCRT, PCIT,                                      &
                                 P_TH_HONC, P_RC_HONC, P_CC_HONC,                       & 
                                 P_CC_SELF,                                             & 
                                 P_RC_AUTO, P_CC_AUTO, P_CR_AUTO,                       & 
                                 P_RC_ACCR, P_CC_ACCR,                                  & 
                                 P_CR_SCBU,                                             & 
                                 P_TH_EVAP, P_RR_EVAP,                                  & 
                                 P_RI_CNVI, P_CI_CNVI,                                  & 
                                 P_TH_DEPS, P_RS_DEPS,                                  & 
                                 P_RI_CNVS, P_CI_CNVS,                                  & 
                                 P_RI_AGGS, P_CI_AGGS,                                  & 
                                 P_TH_DEPG, P_RG_DEPG,                                  & 
                                 P_TH_BERFI, P_RC_BERFI,                                & 
                                 P_TH_RIM, P_RC_RIM, P_CC_RIM, P_RS_RIM, P_RG_RIM,      & 
                                 P_RI_HMS, P_CI_HMS, P_RS_HMS,                          & 
                                 P_TH_ACC, P_RR_ACC, P_CR_ACC, P_RS_ACC, P_RG_ACC,      & 
                                 P_RS_CMEL,                                             & 
                                 P_TH_CFRZ, P_RR_CFRZ, P_CR_CFRZ, P_RI_CFRZ, P_CI_CFRZ, & 
                                 P_TH_WETG, P_RC_WETG, P_CC_WETG, P_RR_WETG, P_CR_WETG, & 
                                 P_RI_WETG, P_CI_WETG, P_RS_WETG, P_RG_WETG, P_RH_WETG, & 
                                 P_TH_DRYG, P_RC_DRYG, P_CC_DRYG, P_RR_DRYG, P_CR_DRYG, & 
                                 P_RI_DRYG, P_CI_DRYG, P_RS_DRYG, P_RG_DRYG,            & 
                                 P_RI_HMG, P_CI_HMG, P_RG_HMG,                          & 
                                 P_TH_GMLT, P_RR_GMLT, P_CR_GMLT,                       & 
!!!     Z_RC_WETH, Z_CC_WETH, Z_RR_WETH, Z_CR_WETH,  &           ! wet growth of hail (WETH) : rc, Nc, rr, Nr, ri, Ni, rs, rg, rh, th
!!!     Z_RI_WETH, Z_CI_WETH, Z_RS_WETH, Z_RG_WETH, Z_RH_WETH, & ! wet growth of hail (WETH) : rc, Nc, rr, Nr, ri, Ni, rs, rg, rh, th
!!!     Z_RG_COHG, &                                             ! conversion of hail into graupel (COHG) : rg, rh
!!!     Z_RR_HMLT, Z_CR_HMLT                                     ! hail melting (HMLT) : rr, Nr, rh=-rr, th
                                 PA_TH, PA_RV, PA_RC, PA_CC, PA_RR, PA_CR,              &
                                 PA_RI, PA_CI, PA_RS, PA_RG, PA_RH,                     &
                                 PEVAP3D                                                )
!
REAL,                 INTENT(IN)    :: PTSTEP 
CHARACTER(LEN=*),     INTENT(IN)    :: HFMFILE 
LOGICAL,              INTENT(IN)    :: OCLOSE_OUT 
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
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT      !
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT      ! Number concentrations (/kg)
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
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_EVAP ! evaporation of rain drops (EVAP) : rr, rv=-rr
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CNVI
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CNVI  ! conversion snow -> ice (CNVI) : ri, Ni, rs=-ri
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_DEPS  ! deposition of vapor on snow (DEPS) : rv=-rs, rs, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CNVS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CNVS  ! conversion ice -> snow (CNVS) : ri, Ni, rs=-ri
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
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_RIM
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_RIM
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_RIM
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_RIM
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_RIM   ! cloud droplet riming (RIM) : rc, Nc, rs, rg, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_HMS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_HMS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_HMS   ! hallett mossop snow (HMS) : ri, Ni, rs
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_ACC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_ACC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_ACC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_ACC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_ACC   ! rain accretion on aggregates (ACC) : rr, Nr, rs, rg, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_CMEL  ! conversion-melting (CMEL) : rs, rg=-rs
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CFRZ  ! rain freezing (CFRZ) : rr, Nr, ri, Ni, rg=-rr-ri, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RH_WETG ! wet growth of graupel (WETG) : rc, NC, rr, Nr, ri, Ni, rs, rg, rh, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_DRYG ! dry growth of graupel (DRYG) : rc, Nc, rr, Nr, ri, Ni, rs, rg, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_HMG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_HMG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_HMG  ! hallett mossop graupel (HMG) : ri, Ni, rg
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_GMLT ! graupel melting (GMLT) : rr, Nr, rg=-rr, th
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
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RG
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RH
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PEVAP3D
!
     END SUBROUTINE LIMA_TENDENCIES
  END INTERFACE
END MODULE MODI_LIMA_TENDENCIES
!#####################################################################
!
!#####################################################################
SUBROUTINE LIMA_TENDENCIES (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,                &
                            PEXNREF, PRHODREF, PPABST, PTHT,                       &
                            PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT,              &
                            PCCT, PCRT, PCIT,                                      &
                            P_TH_HONC, P_RC_HONC, P_CC_HONC,                       & 
                            P_CC_SELF,                                             & 
                            P_RC_AUTO, P_CC_AUTO, P_CR_AUTO,                       & 
                            P_RC_ACCR, P_CC_ACCR,                                  & 
                            P_CR_SCBU,                                             & 
                            P_TH_EVAP, P_RR_EVAP,                                  & 
                            P_RI_CNVI, P_CI_CNVI,                                  & 
                            P_TH_DEPS, P_RS_DEPS,                                  & 
                            P_RI_CNVS, P_CI_CNVS,                                  & 
                            P_RI_AGGS, P_CI_AGGS,                                  & 
                            P_TH_DEPG, P_RG_DEPG,                                  & 
                            P_TH_BERFI, P_RC_BERFI,                                & 
                            P_TH_RIM, P_RC_RIM, P_CC_RIM, P_RS_RIM, P_RG_RIM,      & 
                            P_RI_HMS, P_CI_HMS, P_RS_HMS,                          & 
                            P_TH_ACC, P_RR_ACC, P_CR_ACC, P_RS_ACC, P_RG_ACC,      & 
                            P_RS_CMEL,                                             & 
                            P_TH_CFRZ, P_RR_CFRZ, P_CR_CFRZ, P_RI_CFRZ, P_CI_CFRZ, & 
                            P_TH_WETG, P_RC_WETG, P_CC_WETG, P_RR_WETG, P_CR_WETG, & 
                            P_RI_WETG, P_CI_WETG, P_RS_WETG, P_RG_WETG, P_RH_WETG, & 
                            P_TH_DRYG, P_RC_DRYG, P_CC_DRYG, P_RR_DRYG, P_CR_DRYG, & 
                            P_RI_DRYG, P_CI_DRYG, P_RS_DRYG, P_RG_DRYG,            & 
                            P_RI_HMG, P_CI_HMG, P_RG_HMG,                          & 
                            P_TH_GMLT, P_RR_GMLT, P_CR_GMLT,                       & 
!!!     Z_RC_WETH, Z_CC_WETH, Z_RR_WETH, Z_CR_WETH,  &           ! wet growth of hail (WETH) : rc, Nc, rr, Nr, ri, Ni, rs, rg, rh, th
!!!     Z_RI_WETH, Z_CI_WETH, Z_RS_WETH, Z_RG_WETH, Z_RH_WETH, & ! wet growth of hail (WETH) : rc, Nc, rr, Nr, ri, Ni, rs, rg, rh, th
!!!     Z_RG_COHG, &                                             ! conversion of hail into graupel (COHG) : rg, rh
!!!     Z_RR_HMLT, Z_CR_HMLT                                     ! hail melting (HMLT) : rr, Nr, rh=-rr, th
                            PA_TH, PA_RV, PA_RC, PA_CC, PA_RR, PA_CR,              &
                            PA_RI, PA_CI, PA_RS, PA_RG, PA_RH,                     &
                            PEVAP3D                                                )
!     ######################################################################
!!
!!    PURPOSE
!!    -------
!!    AUTHOR
!!    ------
!!    MODIFICATIONS
!!    -------------
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY : XP00, XRD, XRV, XMD, XMV, XCPD, XCPV, XCL, XCI, XLVTT, XLSTT, XTT, &
                                  XALPW, XBETAW, XGAMW, XALPI, XBETAI, XGAMI
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCTMIN,                                                    &
                                  LCOLD_LIMA, LNUCL_LIMA, LSNOW_LIMA, LHAIL_LIMA, LWARM_LIMA, LACTI_LIMA, LRAIN_LIMA
USE MODD_PARAM_LIMA_WARM,  ONLY : XLBC, XLBEXC, XLBR, XLBEXR
USE MODD_PARAM_LIMA_MIXED, ONLY : XLBG, XLBEXG, XLBH, XLBEXH
USE MODD_PARAM_LIMA_COLD,  ONLY : XSCFAC, XLBI, XLBEXI, XLBS, XLBEXS
!
USE MODI_LIMA_DROPLETS_HOM_FREEZING
USE MODI_LIMA_DROPLETS_SELF_COLLECTION
USE MODI_LIMA_DROPLETS_AUTOCONVERSION
USE MODI_LIMA_DROPLETS_ACCRETION
USE MODI_LIMA_DROPS_SELF_COLLECTION
USE MODI_LIMA_RAIN_EVAPORATION
USE MODI_LIMA_ICE_SNOW_DEPOSITION
USE MODI_LIMA_ICE_AGGREGATION_SNOW
USE MODI_LIMA_GRAUPEL_DEPOSITION
USE MODI_LIMA_BERGERON
USE MODI_LIMA_DROPLETS_RIMING_SNOW
USE MODI_LIMA_RAIN_ACCR_SNOW
USE MODI_LIMA_CONVERSION_MELTING_SNOW
USE MODI_LIMA_RAIN_FREEZING
USE MODI_LIMA_GRAUPEL
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                 INTENT(IN)    :: PTSTEP 
CHARACTER(LEN=*),     INTENT(IN)    :: HFMFILE 
LOGICAL,              INTENT(IN)    :: OCLOSE_OUT 
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
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT      !
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT      ! Number concentrations (/kg)
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
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_EVAP ! evaporation of rain drops (EVAP) : rr, rv=-rr
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CNVI
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CNVI  ! conversion snow -> ice (CNVI) : ri, Ni, rs=-ri
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_DEPS  ! deposition of vapor on snow (DEPS) : rv=-rs, rs, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CNVS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CNVS  ! conversion ice -> snow (CNVS) : ri, Ni, rs=-ri
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
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_RIM
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_RIM
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_RIM
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_RIM
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_RIM   ! cloud droplet riming (RIM) : rc, Nc, rs, rg, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_HMS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_HMS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_HMS   ! hallett mossop snow (HMS) : ri, Ni, rs
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_ACC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_ACC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_ACC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_ACC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_ACC   ! rain accretion on aggregates (ACC) : rr, Nr, rs, rg, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_CMEL  ! conversion-melting (CMEL) : rs, rg=-rs
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CFRZ
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CFRZ  ! rain freezing (CFRZ) : rr, Nr, ri, Ni, rg=-rr-ri, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_WETG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RH_WETG ! wet growth of graupel (WETG) : rc, NC, rr, Nr, ri, Ni, rs, rg, rh, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_DRYG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_DRYG ! dry growth of graupel (DRYG) : rc, Nc, rr, Nr, ri, Ni, rs, rg, th
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_HMG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_HMG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_HMG  ! hallett mossop graupel (HMG) : ri, Ni, rg
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RR_GMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_GMLT ! graupel melting (GMLT) : rr, Nr, rg=-rr, th
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
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RG
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RH
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PEVAP3D
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
!-------------------------------------------------------------------------------
! Pre-compute quantities
!
WHERE (LDCOMPUTE(:))
   ZT(:) = PTHT(:) * ( PPABST(:)/XP00 ) ** (XRD/XCPD)
!
   ZW(:) = PEXNREF(:)*( XCPD &
                               +XCPV*PRVT(:) &
                               +XCL*(PRCT(:)+PRRT(:)) &
                               +XCI*(PRIT(:)+PRST(:)+PRGT(:)+PRHT(:)) )
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
!
ZLBDC(:)  = 1.E10
ZLBDC3(:) = 1.E30
WHERE (PRCT(:)>XRTMIN(2) .AND. PCCT(:)>XCTMIN(2) .AND. LDCOMPUTE(:))
   ZLBDC3(:) = XLBC*PCCT(:) / PRCT(:)
   ZLBDC(:)  = ZLBDC3(:)**XLBEXC
END WHERE
ZLBDR(:)  = 1.E10
ZLBDR3(:) = 1.E30
WHERE (PRRT(:)>XRTMIN(3) .AND. PCRT(:)>XCTMIN(3) .AND. LDCOMPUTE(:))
   ZLBDR3(:) = XLBR*PCRT(:) / PRRT(:)
   ZLBDR(:)  = ZLBDR3(:)**XLBEXR
END WHERE
ZLBDI(:)  = 1.E10
WHERE (PRIT(:)>XRTMIN(4) .AND. PCIT(:)>XCTMIN(4) .AND. LDCOMPUTE(:))
   ZLBDI(:) = ( XLBI*PCIT(:) / PRIT(:) )**XLBEXI
END WHERE
ZLBDS(:)  = 1.E10
WHERE (PRST(:)>XRTMIN(5) .AND. LDCOMPUTE(:) )
   ZLBDS(:) = XLBS*( PRHODREF(:)*PRST(:) )**XLBEXS
END WHERE
ZLBDG(:)  = 1.E10
WHERE (PRGT(:)>XRTMIN(6) .AND. LDCOMPUTE(:) )
   ZLBDG(:) = XLBG*( PRHODREF(:)*PRGT(:) )**XLBEXG
END WHERE
ZLBDH(:)  = 1.E10
WHERE (PRHT(:)>XRTMIN(7) .AND. LDCOMPUTE(:) )
   ZLBDH(:) = XLBH*( PRHODREF(:)*PRHT(:) )**XLBEXH
END WHERE
!
!-------------------------------------------------------------------------------
! Call microphysical processes   
!
IF (LCOLD_LIMA .AND. LWARM_LIMA) THEN
   CALL LIMA_DROPLETS_HOM_FREEZING (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,    &
                                    ZT, ZLVFACT, ZLSFACT,                      &
                                    PRCT, PCCT, ZLBDC,                         &
                                    P_TH_HONC, P_RC_HONC, P_CC_HONC,           &
                                    PA_TH, PA_RC, PA_CC, PA_RI, PA_CI          )
END IF
!
IF (LWARM_LIMA) THEN
   CALL LIMA_DROPLETS_SELF_COLLECTION (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                       PRHODREF,                       &
                                       PCCT, ZLBDC3,                   &
                                       P_CC_SELF,                      &
                                       PA_CC                           )
END IF
!
IF (LWARM_LIMA .AND. LRAIN_LIMA) THEN
   CALL LIMA_DROPLETS_AUTOCONVERSION (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                      PRHODREF,                       &
                                      PRCT, ZLBDC, ZLBDR,             &
                                      P_RC_AUTO, P_CC_AUTO, P_CR_AUTO,&
                                      PA_RC, PA_CC, PA_RR, PA_CR      )
END IF
!
IF (LWARM_LIMA .AND. LRAIN_LIMA) THEN
   CALL LIMA_DROPLETS_ACCRETION (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                 PRHODREF,                       &
                                 PRCT, PRRT, PCCT, PCRT,         &
                                 ZLBDC, ZLBDC3, ZLBDR, ZLBDR3,   &
                                 P_RC_ACCR, P_CC_ACCR,           &
                                 PA_RC, PA_CC, PA_RR             )
END IF
!
IF (LWARM_LIMA .AND. LRAIN_LIMA) THEN 
   CALL LIMA_DROPS_SELF_COLLECTION (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                    PRHODREF,                       &
                                    PCRT, ZLBDR, ZLBDR3,            &
                                    P_CR_SCBU,                      &
                                    PA_CR                           )
END IF
!
IF (LWARM_LIMA .AND. LRAIN_LIMA) THEN
   CALL LIMA_RAIN_EVAPORATION (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,     &
                               PRHODREF, ZT, ZLV, ZLVFACT, ZEVSAT, ZRVSAT, &
                               PRVT, PRCT, PRRT, ZLBDR,                    &
                               P_TH_EVAP, P_RR_EVAP,                       &
                               PA_RV, PA_RR, PA_TH,                        &
                               PEVAP3D                                     )
END IF
!
IF (LCOLD_LIMA .AND. LSNOW_LIMA) THEN
   !
   ! Includes vapour deposition on snow, ice -> snow and snow -> ice exchanges
   !
   CALL LIMA_ICE_SNOW_DEPOSITION (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,   &
                                  PRHODREF, ZSSI, ZAI, ZCJ, ZLSFACT,        &
                                  PRIT, PRST, PCIT, ZLBDI, ZLBDS,           &
                                  P_RI_CNVI, P_CI_CNVI,                     &
                                  P_TH_DEPS, P_RS_DEPS,                     &
                                  P_RI_CNVS, P_CI_CNVS,                     &
                                  PA_TH, PA_RV, PA_RI, PA_CI, PA_RS         )
END IF
!
IF (LCOLD_LIMA .AND. LSNOW_LIMA) THEN
   CALL LIMA_ICE_AGGREGATION_SNOW (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                   ZT,                             &
                                   PRIT, PRST, PCIT, ZLBDI, ZLBDS, &
                                   P_RI_AGGS, P_CI_AGGS,           &
                                   PA_RI, PA_CI, PA_RS             )
END IF
!
IF (LWARM_LIMA .AND. LCOLD_LIMA) THEN
   CALL LIMA_GRAUPEL_DEPOSITION (HFMFILE, OCLOSE_OUT, LDCOMPUTE,       &
                                 PRGT, ZSSI, ZLBDG, ZAI, ZCJ, ZLSFACT, &
                                 P_TH_DEPG, P_RG_DEPG,                 &
                                 PA_TH, PA_RV, PA_RG                   )
END IF
!
IF (LWARM_LIMA .AND. LCOLD_LIMA) THEN
   CALL LIMA_BERGERON (HFMFILE, OCLOSE_OUT, LDCOMPUTE,    &
                       PRCT, PRIT, PCIT, ZLBDI,           &
                       ZSSIW, ZAI, ZCJ, ZLVFACT, ZLSFACT, &
                       P_TH_BERFI, P_RC_BERFI,            &
                       PA_TH, PA_RC, PA_RI                )
END IF
!
IF (LWARM_LIMA .AND. LCOLD_LIMA .AND. LSNOW_LIMA) THEN
     !
     ! Graupel production as tendency (or should be tendency + instant to stick to the previous version ?)
     ! Includes the Hallett Mossop process for riming of droplets by snow (HMS)
     !
   CALL LIMA_DROPLETS_RIMING_SNOW (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,                  &
                                   PRHODREF, ZT,                                            &
                                   PRCT, PCCT, PRST, ZLBDC, ZLBDS, ZLVFACT, ZLSFACT,        &
                                   P_TH_RIM, P_RC_RIM, P_CC_RIM, P_RS_RIM, P_RG_RIM,        &
                                   P_RI_HMS, P_CI_HMS, P_RS_HMS,                            &
                                   PA_TH, PA_RC, PA_CC, PA_RI, PA_CI, PA_RS, PA_RG          )
END IF
!
IF (LWARM_LIMA .AND. LRAIN_LIMA .AND. LCOLD_LIMA .AND. LSNOW_LIMA) THEN
   CALL LIMA_RAIN_ACCR_SNOW (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,                  &
                             PRHODREF, ZT,                                            &
                             PRRT, PCRT, PRST, ZLBDR, ZLBDS, ZLVFACT, ZLSFACT,        &
                             P_TH_ACC, P_RR_ACC, P_CR_ACC, P_RS_ACC, P_RG_ACC,        &
                             PA_TH, PA_RR, PA_CR, PA_RS, PA_RG                        )
END IF
!
IF (LWARM_LIMA .AND. LCOLD_LIMA .AND. LSNOW_LIMA) THEN
   !
   ! Conversion melting of snow should account for collected droplets and drops where T>0C, but does not !
   ! Some thermodynamical computations inside, to externalize ?
   !
   CALL LIMA_CONVERSION_MELTING_SNOW (HFMFILE, OCLOSE_OUT, LDCOMPUTE,      &
                                      PRHODREF, PPABST, ZT, ZKA, ZDV, ZCJ, &
                                      PRVT, PRST, ZLBDS,                   &
                                      P_RS_CMEL,                           &
                                      PA_RS, PA_RG                         )
END IF
!
IF (LWARM_LIMA .AND. LRAIN_LIMA .AND. LCOLD_LIMA ) THEN
   CALL LIMA_RAIN_FREEZING (HFMFILE, OCLOSE_OUT, LDCOMPUTE,                        &
                            PRHODREF, ZT, ZLVFACT, ZLSFACT,                        &
                            PRRT, PCRT, PRIT, PCIT, ZLBDR,                         &
                            P_TH_CFRZ, P_RR_CFRZ, P_CR_CFRZ, P_RI_CFRZ, P_CI_CFRZ, &
                            PA_TH, PA_RR, PA_CR, PA_RI, PA_CI, PA_RG               )
END IF
!
IF (LWARM_LIMA .AND. LCOLD_LIMA) THEN
     !
     ! Melting of graupel should account for collected droplets and drops where T>0C, but does not !
     ! Collection and water shedding should also happen where T>0C, but do not !
     ! Hail production as tendency (should be instant to stick to the previous version ?)
     ! Includes Hallett-Mossop  process for riming of droplets by graupel (HMG)
     ! Some thermodynamical computations inside, to externalize ?
     !
   CALL LIMA_GRAUPEL (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,                &
                      PRHODREF, PPABST, ZT, ZKA, ZDV, ZCJ,                   &
                      PRVT, PRCT, PRRT, PRIT, PRST, PRGT,                    &
                      PCCT, PCRT, PCIT,                                      &
                      ZLBDC, ZLBDR, ZLBDS, ZLBDG,                            &
                      ZLVFACT, ZLSFACT,                                      &
                      P_TH_WETG, P_RC_WETG, P_CC_WETG, P_RR_WETG, P_CR_WETG, &
                      P_RI_WETG, P_CI_WETG, P_RS_WETG, P_RG_WETG, P_RH_WETG, &
                      P_TH_DRYG, P_RC_DRYG, P_CC_DRYG, P_RR_DRYG, P_CR_DRYG, &
                      P_RI_DRYG, P_CI_DRYG, P_RS_DRYG, P_RG_DRYG,            &
                      P_RI_HMG, P_CI_HMG, P_RG_HMG,                          &
                      P_TH_GMLT, P_RR_GMLT, P_CR_GMLT,                       &
                      PA_TH, PA_RC, PA_CC, PA_RR, PA_CR,                     &
                      PA_RI, PA_CI, PA_RS, PA_RG, PA_RH                      )
END IF
!
IF (LWARM_LIMA .AND. LCOLD_LIMA .AND. LHAIL_LIMA) THEN
!     CALL LIMA_HAIL_GROWTH

!     CALL LIMA_HAIL_CONVERSION

!     CALL LIMA_HAIL_MELTING
END IF
   !  
END SUBROUTINE LIMA_TENDENCIES
