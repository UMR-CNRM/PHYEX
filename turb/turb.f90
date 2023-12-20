!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
      SUBROUTINE TURB(CST,CSTURB,BUCONF,TURBN,NEBN,D,TLES,            &
              & KRR,KRRL,KRRI,HLBCX,HLBCY,KGRADIENTS,KHALO,           &
              & KSPLIT, OCLOUDMODIFLM, KSV,KSV_LGBEG,KSV_LGEND,       &
              & KSV_LIMA_NR, KSV_LIMA_NS, KSV_LIMA_NG, KSV_LIMA_NH,   &
              & O2D,ONOMIXLG,OFLAT,OCOUPLES,OBLOWSNOW,OIBM,OFLYER,    &
              & OCOMPUTE_SRC, PRSNOW,                                 &
              & OOCEAN,ODEEPOC,ODIAG_IN_RUN,                          &
              & HTURBLEN_CL,HCLOUD,HELEC,                             &
              & PTSTEP,TPFILE,                                        &
              & PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,                         &
              & PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,PCOSSLOPE,PSINSLOPE,    &
              & PRHODJ,PTHVREF,PHGRAD,PZS,                            &
              & PSFTH,PSFRV,PSFSV,PSFU,PSFV,                          &
              & PPABST,PUT,PVT,PWT,PTKET,PSVT,PSRCT,                  &
              & PLENGTHM,PLENGTHH,MFMOIST,                            &
              & PBL_DEPTH,PSBL_DEPTH,                                 &
              & PCEI,PCEI_MIN,PCEI_MAX,PCOEF_AMPL_SAT,                &
              & PTHLT,PRT,                                            &
              & PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,PRTKES,              &
              & PSIGS,                                                &
              & PFLXZTHVMF,PWTH,PWRC,PWSV,PDP,PTP,PTDIFF,PTDISS,      &
              & TBUDGETS, KBUDGETS,                                   &
              & PEDR,PLEM,PRTKEMS,PTPMF,                              &
              & PDRUS_TURB,PDRVS_TURB,                                &
              & PDRTHLS_TURB,PDRRTS_TURB,PDRSVS_TURB,PTR,PDISS,       &
              & PIBM_LS, PIBM_XMUT,                                   &
              & PCURRENT_TKE_DISS, PSSTFL, PSSTFL_C, PSSRFL_C,        &
              & PSSUFL_C, PSSVFL_C,PSSUFL,PSSVFL                      )
!     #################################################################
!
!
!!****  *TURB* - computes the turbulent source terms for the prognostic
!!               variables.
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to compute the source terms in
!!    the evolution equations due to the turbulent mixing.
!!      The source term is computed as the divergence of the turbulent fluxes.
!!    The cartesian fluxes are obtained by a one and a half order closure, based
!!    on a prognostic equation for the Turbulence Kinetic Energy( TKE ). The
!!    system is closed by prescribing a turbulent mixing length. Different
!!    choices are available for this length.
!
!!**  METHOD
!!    ------
!!
!!      The dimensionality of the turbulence parameterization can be chosen by
!!    means of the parameter TURBN%CTURBDIM:
!!           * TURBN%CTURBDIM='1DIM' the parameterization is 1D but can be used in
!!    3D , 2D or 1D simulations. Only the sources associated to the vertical
!!    turbulent fluxes are taken into account.
!!           *  TURBN%CTURBDIM='3DIM' the parameterization is fully 2D or 3D depending
!!    on the model  dimensionality. Of course, it does not make any sense to
!!    activate this option with a 1D model.
!!
!!      The following steps are made:
!!      1- Preliminary computations.
!!      2- The metric coefficients are recovered from the grid knowledge.
!!      3- The mixing length is computed according to its choice:
!!           * TURBN%CTURBLEN='BL89' the Bougeault and Lacarrere algorithm is used.
!!             The mixing length is given by the vertical displacement from its
!!             original level of an air particule having an initial internal
!!             energy equal to its TKE and stopped by the buoyancy forces.
!!             The discrete formulation is second order accurate.
!!           * TURBN%CTURBLEN='DELT' the mixing length is given by the mesh size
!!             depending on the model dimensionality, this length is limited
!!             with the ground distance.
!!           * TURBN%CTURBLEN='DEAR' the mixing length is given by the mesh size
!!             depending on the model dimensionality, this length is limited
!!             with the ground distance and also by the Deardorff mixing length
!!             pertinent in the stable cases.
!!           * TURBN%CTURBLEN='KEPS' the mixing length is deduced from the TKE
!!             dissipation, which becomes a prognostic variable of the model (
!!             Duynkerke formulation).
!!      3'- The cloud mixing length is computed according to HTURBLEN_CLOUD
!!             and emphasized following the CEI index
!!      4- The conservative variables are computed along with Lv/Cp.
!!      5- The turbulent Prandtl numbers are computed from the resolved fields
!!         and TKE
!!      6- The sources associated to the vertical turbulent fluxes are computed
!!      with a temporal scheme allowing a degree of implicitness given by
!!      TURBN%XIMPL, varying from TURBN%XIMPL=0. ( purely explicit scheme) to TURBN%XIMPL=1.
!!      ( purely implicit scheme)
!!      The sources associated to the horizontal fluxes are computed with a
!!      purely explicit temporal scheme. These sources are only computed when
!!      the turbulence parameterization is 2D or 3D( TURBN%CTURBDIM='3DIM' ).
!!      7- The sources for TKE are computed, along with the dissipation of TKE
!!      if TURBN%CTURBLEN='KEPS'.
!!      8- Some turbulence-related quantities are stored in the synchronous
!!      FM-file.
!!      9- The non-conservative variables are retrieved.
!!
!!
!!      The saving of the fields in the synchronous FM-file is controlled by:
!!        * TURBN%LTURB_FLX => saves all the turbulent fluxes and correlations
!!        * TURBN%LTURB_DIAG=> saves the turbulent Prandtl and Schmidt numbers, the
!!                       source terms of TKE and dissipation of TKE
!!
!!    EXTERNAL
!!    --------
!!      SUBROUTINE PRANDTL   : computes the turbulent Prandtl number
!!      SUBROUTINE TURB_VER  : computes the sources from the vertical fluxes
!!      SUBROUTINE TURB_HOR  : computes the sources from the horizontal fluxes
!!      SUBROUTINE TKE_EPS_SOURCES : computes the sources for  TKE and its
!!                                   dissipation
!!      SUBROUTINE BUDGET    : computes and stores the budgets
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!       MODD_PARAMETERS : JPVEXT_TURB  number of marginal vertical points
!!
!!       MODD_CONF      : CCONF model configuration (start/restart)
!!                        L1D   switch for 1D model version
!!                        L2D   switch for 2D model version
!!
!!       MODD_CST  : contains physical constants
!!                    CST%XG   gravity constant
!!                    CST%XRD  Gas constant for dry air
!!                    CST%XRV  Gas constant for vapor
!!
!!       MODD_CTURB : contains turbulence scheme constants
!!                    XCMFS,XCED       to compute the dissipation mixing length
!!                    XTKEMIN  minimum values for the TKE
!!                    CST%XLINI,CST%XLINF      to compute Bougeault-Lacarrere mixing
!!                                     length
!!      Module MODD_BUDGET:
!!         NBUMOD
!!         CBUTYPE
!!         LBU_RU
!!         LBU_RV
!!         LBU_RW
!!         LBU_RTH
!!         LBU_RSV1
!!         LBU_RRV
!!         LBU_RRC
!!         LBU_RRR
!!         LBU_RRI
!!         LBU_RRS
!!         LBU_RRG
!!         LBU_RRH
!!
!!    REFERENCE
!!    ---------
!!      Book 2 of documentation (routine TURB)
!!      Book 1 of documentation (Chapter: Turbulence)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart             * INM and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         05/10/94
!!      Modifications: Feb 14, 1995 (J.Cuxart and J.Stein)
!!                                  Doctorization and Optimization
!!      Modifications: March 21, 1995 (J.M. Carriere)
!!                                  Introduction of cloud water
!!      Modifications: June   1, 1995 (J.Cuxart     )
!!                                  take min(Kz,delta)
!!      Modifications: June   1, 1995 (J.Stein J.Cuxart)
!!                                  remove unnecessary arrays and change Prandtl
!!                                  and Schmidt numbers localizations
!!      Modifications: July  20, 1995 (J.Stein) remove MODI_ground_ocean +
!!                                TZDTCUR + MODD_TIME because they are not used
!!                                change RW in RNP for the outputs
!!      Modifications: August 21, 1995 (Ph. Bougeault)
!!                                  take min(K(z-zsol),delta)
!!      Modifications: Sept 14, 1995 (Ph Bougeault, J. Cuxart)
!!         second order BL89 mixing length computations + add Deardorff length
!!         in the Delta case for stable cases
!!      Modifications: Sept 19, 1995 (J. Stein, J. Cuxart)
!!         define a DEAR case for the mixing length, add MODI_BUDGET and change
!!         some BUDGET calls, add LES tools
!!      Modifications: Oct  16, 1995 (J. Stein) change the budget calls
!!      Modifications: Feb  28, 1996 (J. Stein) optimization +
!!                                              remove min(K(z-zsol),delta)+
!!                                              bug in the tangential fluxes
!!      Modifications: Oct  16, 1996 (J. Stein) change the subgrid condensation
!!                                              scheme + temporal discretization
!!      Modifications: Dec  19, 1996 (J.-P. Pinty) update the budget calls
!!                     Jun  22, 1997 (J. Stein) use the absolute pressure and
!!                                  change the Deardorf length at the surface
!!      Modifications: Apr  27, 1997 (V. Masson) BL89 mix. length computed in
!!                                               a separate routine
!!                     Oct  13, 1999 (J. Stein)  switch for the tgt fluxes
!!                     Jun  24, 1999 (P Jabouille)  Add routine UPDATE_ROTATE_WIND
!!                     Feb  15, 2001 (J. Stein)  remove tgt fluxes
!!                     Mar 8,  2001 (V. Masson) forces the same behaviour near the surface
!!                                              for all mixing lengths
!!                     Nov 06, 2002 (V. Masson) LES budgets
!!                     Nov,    2002 (V. Masson) implement modifications of
!!                                              mixing and dissipative lengths
!!                                              near the surface (according
!!                                              Redelsperger et al 2001)
!!                     Apr,    2003 (V. Masson) bug in Blackadar length
!!                                              bug in LES in 1DIM case
!!                     Feb 20, 2003 (J.-P. Pinty) Add reversible ice processes
!!                     May,26  2004 (P Jabouille) coef for computing dissipative heating
!!                     Sept 2004 (M.Tomasini) Cloud Mixing length modification
!!                                            following the instability
!!                                            criterium CEI calculated in modeln
!!                     May   2006    Remove KEPS
!!                     Sept.2006 (I.Sandu): Modification of the stability criterion for
!!                                 DEAR (theta_v -> theta_l)
!!                     Oct 2007 (J.Pergaud) Add MF contribution for vert. turb. transport
!!                     Oct.2009  (C.Lac) Introduction of different PTSTEP according to the
!!                              advection schemes
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!                     06/2011 (J.escobar ) Bypass Bug with ifort11/12 on  HLBCX,HLBC
!!                     2012-02 Y. Seity,  add possibility to run with reversed
!!                                          vertical levels
!!                     10/2012 (J. Colin) Correct bug in DearDoff for dry simulations
!!                     10/2012 J.Escobar Bypass PGI bug , redefine some allocatable array inplace of automatic
!!                     2014-11 Y. Seity,  add output terms for TKE DDHs budgets
!!                     July 2015 (Wim de Rooy)  modifications to run with RACMO
!!                                              turbulence (TURBN%LHARAT=TRUE)
!!                     04/2016  (C.Lac) correction of negativity for KHKO
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  Q. Rodier      01/2018: introduction of RM17
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  Wim de Rooy    06/2019: update statistical cloud scheme
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!  B. Vie         03/2020: LIMA negativity checks after turbulence, advection and microphysics budgets
!  P. Wautelet 11/06/2020: bugfix: correct PRSVS array indices
!  P. Wautelet + Benoit Vié 06/2020: improve removal of negative scalar variables + adapt the corresponding budgets
!  P. Wautelet 30/06/2020: move removal of negative scalar variables to Sources_neg_correct
!  R. Honnert/V. Masson 02/2021: new mixing length in the grey zone
!  J.L. Redelsperger 03/2021: add Ocean LES case
!  C. Barthe   08/02/2022: add helec in arguments of Sources_neg_correct
! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODE_SHUMAN_PHY, ONLY: MZF_PHY,MXF_PHY,MYF_PHY
USE YOMHOOK ,   ONLY: LHOOK, DR_HOOK, JPHOOK
!
USE MODD_BUDGET,     ONLY:  NBUDGET_U,  NBUDGET_V,  NBUDGET_W,  NBUDGET_TH, NBUDGET_RV, NBUDGET_RC,  &
                             NBUDGET_RI, NBUDGET_SV1, &
                            TBUDGETDATA, TBUDGETCONF_t
USE MODD_CST,        ONLY: CST_t
USE MODD_CTURB,      ONLY: CSTURB_t
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_FIELD,      ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_IO,         ONLY: TFILEDATA
USE MODD_LES,        ONLY: TLES_t
USE MODD_PARAMETERS, ONLY: JPVEXT_TURB, XUNDEF
USE MODD_TURB_n,     ONLY: TURB_t
USE MODD_NEB_n,      ONLY: NEB_t
!
USE MODE_BL89,                ONLY: BL89
USE MODE_BUDGET_PHY,              ONLY: BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY
USE MODE_EMOIST,              ONLY: EMOIST
USE MODE_ETHETA,              ONLY: ETHETA
USE MODE_GRADIENT_U_PHY,      ONLY: GZ_U_UW_PHY
USE MODE_GRADIENT_V_PHY,      ONLY: GZ_V_VW_PHY
USE MODE_GRADIENT_W_PHY,      ONLY: GZ_W_M_PHY
USE MODE_GRADIENT_M_PHY,      ONLY: GZ_M_W_PHY
USE MODE_IBM_MIXINGLENGTH,    ONLY: IBM_MIXINGLENGTH
USE MODE_IO_FIELD_WRITE_PHY,      ONLY: IO_FIELD_WRITE_PHY
USE MODE_RMC01,               ONLY: RMC01
USE MODE_ROTATE_WIND,         ONLY: ROTATE_WIND, UPDATE_ROTATE_WIND
USE MODE_SBL_PHY,             ONLY: LMO
USE MODE_SOURCES_NEG_CORRECT, ONLY: SOURCES_NEG_CORRECT_PHY
USE MODE_TM06,                ONLY: TM06
USE MODE_TKE_EPS_SOURCES,     ONLY: TKE_EPS_SOURCES
USE MODE_TURB_HOR_SPLT,       ONLY: TURB_HOR_SPLT
USE MODE_TURB_VER,            ONLY: TURB_VER
USE MODE_UPDATE_LM,           ONLY: UPDATE_LM
!
USE MODI_LES_MEAN_SUBGRID_PHY
!
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D             ! PHYEX variables dimensions structure
TYPE(CST_t),            INTENT(IN)   :: CST           ! modd_cst general constant structure
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB        ! modd_csturb turb constant structure
TYPE(TBUDGETCONF_t),    INTENT(IN)   :: BUCONF        ! budget structure
TYPE(TURB_t),           INTENT(IN)   :: TURBN         ! modn_turbn (turb namelist) structure
TYPE(NEB_t),            INTENT(IN)   :: NEBN          ! modd_nebn structure
TYPE(TLES_t),           INTENT(INOUT)   :: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KGRADIENTS    ! Number of stored horizontal gradients
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV, KSV_LGBEG, KSV_LGEND ! number of scalar variables
INTEGER,                INTENT(IN)   :: KSV_LIMA_NR,KSV_LIMA_NS,KSV_LIMA_NG,KSV_LIMA_NH
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
INTEGER,                INTENT(IN)   :: KSPLIT        ! number of time-splitting
LOGICAL,                INTENT(IN)   :: OCLOUDMODIFLM ! cloud mixing length modifications
INTEGER,                INTENT(IN)   ::  KHALO        ! Size of the halo for parallel distribution
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  ODEEPOC      ! activates sfc forcing for ideal ocean deep conv
LOGICAL,                INTENT(IN)   ::  OFLYER       ! MesoNH flyer diagnostic
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES version 
LOGICAL,                INTENT(IN)   ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
LOGICAL,                INTENT(IN)   ::  ODIAG_IN_RUN ! switch to activate online diagnostics (mesonh)
LOGICAL,                INTENT(IN)   ::  OIBM         ! switch to modity mixing length near building with IBM
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBLEN_CL  ! kind of cloud mixing length
CHARACTER (LEN=4),      INTENT(IN)   ::  HCLOUD       ! Kind of microphysical scheme
CHARACTER (LEN=4),      INTENT(IN)   ::  HELEC        ! Kind of cloud electricity scheme

REAL,                   INTENT(IN)   ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
REAL,                   INTENT(IN)   ::  PTSTEP       ! timestep
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PDXX,PDYY,PDZZ,PDZX,PDZY
                                        ! metric coefficients
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PZZ       !  physical distance
! between 2 succesive grid points along the K direction
REAL, DIMENSION(D%NIJT),   INTENT(IN)      ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PZS ! orography (for LEONARD terms)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  PRHODJ    ! dry density * Grid size
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  MFMOIST ! moist mass flux dual scheme
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  PTHVREF   ! Virtual Potential
                                        ! Temperature of the reference state
REAL, DIMENSION(D%NIJT,D%NKT,KGRADIENTS),   INTENT(IN) ::  PHGRAD      ! horizontal gradients
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)      ::  PSFTH,PSFRV,   &
! normal surface fluxes of theta and Rv
                                            PSFU,PSFV
! normal surface fluxes of (u,v) parallel to the orography
REAL, DIMENSION(D%NIJT,KSV), INTENT(IN)      ::  PSFSV
! normal surface fluxes of Scalar var.
!
!    prognostic variables at t- deltat
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PPABST      ! Pressure at time t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PUT,PVT,PWT ! wind components
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PTKET       ! TKE
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVT        ! passive scal. var.
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)),   INTENT(IN) ::  PSRCT       ! Second-order flux
                      ! s'rc'/2Sigma_s2 at time t-1 multiplied by Lambda_3
REAL, DIMENSION(MERGE(D%NIJT,0,TURBN%CTOM=='TM06')),INTENT(INOUT) :: PBL_DEPTH  ! BL height for TOMS
REAL, DIMENSION(MERGE(D%NIJT,0,TURBN%LRMC01)),INTENT(INOUT) :: PSBL_DEPTH ! SBL depth for RMC01
!
!    variables for cloud mixing length
REAL, DIMENSION(MERGE(D%NIJT,0,OCLOUDMODIFLM),&
                MERGE(D%NKT,0,OCLOUDMODIFLM)),INTENT(IN)      ::  PCEI
                                                 ! Cloud Entrainment instability
                                                 ! index to emphasize localy
                                                 ! turbulent fluxes
REAL, INTENT(IN)      ::  PCEI_MIN ! minimum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCEI_MAX ! maximum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCOEF_AMPL_SAT ! saturation of the amplification coefficient
!
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PTHLT       ! conservative pot. temp.
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT) ::  PRT         ! water var.  where
                             ! PRT(:,:,:,1) is the conservative mixing ratio
!
! sources of momentum, conservative potential temperature, Turb. Kin. Energy,
! TKE dissipation
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PRUS,PRVS,PRWS,PRTHLS,PRTKES
! Source terms for all water kinds, PRRS(:,:,:,1) is used for the conservative
! mixing ratio
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN),OPTIONAL    ::  PRTKEMS
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT) ::  PRRS
! Source terms for all passive scalar variables
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(INOUT) ::  PRSVS
! Sigma_s at time t+1 : square root of the variance of the deviation to the
! saturation
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(OUT)     ::  PSIGS
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRUS_TURB   ! evolution of rhoJ*U   by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRVS_TURB   ! evolution of rhoJ*V   by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRTHLS_TURB ! evolution of rhoJ*thl by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRRTS_TURB  ! evolution of rhoJ*rt  by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT),OPTIONAL ::  PDRSVS_TURB  ! evolution of rhoJ*Sv  by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  PFLXZTHVMF
!                                           MF contribution for vert. turb. transport
!                                           used in the buoy. prod. of TKE
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PWTH       ! heat flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PWRC       ! cloud water flux
REAL, DIMENSION(D%NIJT,D%NKT,KSV),INTENT(OUT) :: PWSV       ! scalar flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTP        ! Thermal TKE production
                                                   ! MassFlux + turb
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL  :: PTPMF      ! Thermal TKE production
                                                   ! MassFlux Only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PDP        ! Dynamic TKE production
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTDIFF     ! Diffusion TKE term
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTDISS     ! Dissipation TKE term
!
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
LOGICAL, INTENT(IN) :: ONOMIXLG          ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL, INTENT(IN) :: O2D               ! Logical for 2D model version (modd_conf)
!
! length scale from vdfexcu
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PLENGTHM, PLENGTHH
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT), OPTIONAL  :: PEDR  ! EDR
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT), OPTIONAL  :: PLEM  ! Mixing length
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT), OPTIONAL  ::  PTR          ! Transport prod. of TKE
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT), OPTIONAL  ::  PDISS        ! Dissipation of TKE
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT), OPTIONAL  ::  PCURRENT_TKE_DISS ! if ODIAG_IN_RUN in mesonh
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSTFL        ! Time evol Flux of T at sea surface (LOCEAN)
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSTFL_C  ! O-A interface flux for theta(LOCEAN and LCOUPLES)
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSRFL_C  ! O-A interface flux for vapor (LOCEAN and LCOUPLES) 
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSUFL_C        ! Time evol Flux of U at sea surface (LOCEAN)
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSVFL_C  !
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSUFL   
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSVFL  !
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT), OPTIONAL :: PIBM_XMUT ! IBM turbulent viscosity
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN), OPTIONAL  :: PIBM_LS ! IBM Level-set function
!
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT) ::     &
          ZCP,                        &  ! Cp at t-1
          ZEXN,                       &  ! EXN at t-1
          ZT,                         &  ! T at t-1
          ZLOCPEXNM,                  &  ! Lv/Cp/EXNREF at t-1
          ZLM,ZLMW,                   &  ! Turbulent mixing length (+ work array)
          ZLEPS,                      &  ! Dissipative length
          ZTRH,                       &  !
          ZATHETA,ZAMOIST,            &  ! coefficients for s = f (Thetal,Rnp)
          ZCOEF_DISS,                 &  ! 1/(Cph*Exner) for dissipative heating
          ZFRAC_ICE,                  &  ! ri fraction of rc+ri
          ZMWTH,ZMWR,ZMTH2,ZMR2,ZMTHR,&  ! 3rd order moments
          ZFWTH,ZFWR,ZFTH2,ZFR2,ZFTHR,&  ! opposite of verticale derivate of 3rd order moments
          ZTHLM,ZRTKEMS,              &  ! initial potential temp; TKE advective source
          ZSHEAR, ZDUDZ, ZDVDZ,       &  ! horizontal-wind vertical gradient
          ZLVOCPEXNM,ZLSOCPEXNM,      &  ! Lv/Cp/EXNREF and Ls/Cp/EXNREF at t-1
          ZATHETA_ICE,ZAMOIST_ICE,    &  ! coefficients for s = f (Thetal,Rnp)
          ZRVSAT, ZDRVSATDT,          &  ! local array for routine compute_function_thermo
          ZWORK1,ZWORK2,              &  ! working array syntax
          ZETHETA,ZEMOIST,            &  ! coef ETHETA and EMOIST (for DEAR routine)
          ZDTHLDZ,ZDRTDZ,             &  ! dtheta_l/dz, drt_dz used for computing the stablity criterion
          ZCOEF_AMPL,                 &  ! Amplification coefficient of the mixing length
                                         ! when the instability criterium is verified (routine CLOUD_MODIF_LM)
          ZLM_CLOUD                      ! Turbulent mixing length in the clouds (routine CLOUD_MODIF_LM)
!
!
REAL, DIMENSION(D%NIJT,D%NKT,KRR) :: ZRM ! initial mixing ratio
REAL, DIMENSION(D%NIJT) ::  ZTAU11M,ZTAU12M,  &
                                                 ZTAU22M,ZTAU33M,  &
            ! tangential surface fluxes in the axes following the orography
                                                 ZUSLOPE,ZVSLOPE,  &
            ! wind components at the first mass level parallel
            ! to the orography
                                                 ZCDUEFF,          &
            ! - Cd*||u|| where ||u|| is the module of the wind tangential to
            ! orography (ZUSLOPE,ZVSLOPE) at the surface.
                                                 ZUSTAR, ZLMO,     &
                                                 ZRVM, ZSFRV,ZWORK2D
            ! friction velocity, Monin Obuhkov length, work arrays for vapor
!
            ! Virtual Potential Temp. used
            ! in the Deardorff mixing length computation
!
!with LIMA, do not change rain, snow, graupel and hail concentrations (mixing ratio is not changed)
REAL, DIMENSION(D%NIJT,D%NKT,KSV) :: ZRSVS
!
REAL                :: ZEXPL        ! 1-TURBN%XIMPL deg of expl.
REAL                :: ZRVORD       ! RV/RD
REAL                :: ZEPS         ! XMV / XMD
REAL                :: ZD           ! distance to the surface (for routine DELT)
REAL                :: ZVAR         ! Intermediary variable (for routine DEAR)
REAL                :: ZPENTE       ! Slope of the amplification straight line (for routine CLOUD_MODIF_LM)
REAL                :: ZCOEF_AMPL_CEI_NUL! Ordonnate at the origin of the
                                         ! amplification straight line (for routine CLOUD_MODIF_LM)
!
INTEGER             :: IIJB,IIJE,IKB,IKE      ! index value for the
! Beginning and the End of the physical domain for the mass points
INTEGER             :: IKT,IKA,IKU  ! array size in k direction
INTEGER             :: IKL
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain
INTEGER             :: JRR,JK,JSV   ! loop counters
INTEGER             :: JIJ          ! loop counters
REAL                :: ZL0          ! Max. Mixing Length in Blakadar formula
REAL                :: ZALPHA       ! work coefficient :
                                    ! - proportionnality constant between Dz/2 and
!                                   !   BL89 mixing length near the surface
!
REAL :: ZTIME1, ZTIME2
TYPE(TFIELDMETADATA) :: TZFIELD
!
REAL, DIMENSION(D%NIJT,D%NKT,KSV+KRR) :: ZWORKT
REAL, DIMENSION(D%NIJT,D%NKT,KSV+KRR) :: ZWORKS
REAL, DIMENSION(D%NIJT,      KSV+KRR) :: ZWORKSFSV
REAL, DIMENSION(D%NIJT,D%NKT,KSV+KRR) :: ZWORKWSV
!
!*      1.PRELIMINARIES
!         -------------
!
!*      1.1 Set the internal domains, ZEXPL
!
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE,ZHOOK_HANDLE2
IF (LHOOK) CALL DR_HOOK('TURB',0,ZHOOK_HANDLE)
!
IF (TURBN%LHARAT .AND. TURBN%CTURBDIM /= '1DIM') THEN
  CALL ABOR1('TURBN%LHARATU only implemented for option TURBN%CTURBDIM=1DIM!')
ENDIF
IF (TURBN%LHARAT .AND. TLES%LLES_CALL) THEN
  CALL ABOR1('TURBN%LHARATU not implemented for option LLES_CALL')
ENDIF
!
IKT=D%NKT
IKTB=D%NKTB
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IKA=D%NKA
IKU=D%NKU
IKL=D%NKL
IIJE=D%NIJE
IIJB=D%NIJB
!
ZEXPL = 1.- TURBN%XIMPL
ZRVORD= CST%XRV / CST%XRD
!
!Copy data into ZTHLM and ZRM only if needed
IF (TURBN%CTURBLEN=='BL89' .OR. TURBN%CTURBLEN=='RM17' .OR. TURBN%CTURBLEN=='HM21' .OR. TURBN%LRMC01) THEN
  ZTHLM(:,:) = PTHLT(:,:)
  ZRM(:,:,:) = PRT(:,:,:)
END IF
!
!Save LIMA scalar variables sources
ZRSVS(:,:,1:KSV)=PRSVS(:,:,1:KSV)
!
ZWORKT(:,:,1:KSV)=PSVT(:,:,:)
ZWORKS(:,:,1:KSV)=PRSVS(:,:,:)
ZWORKT(:,:,KSV+1:KSV+KRR)=PRT(:,:,:)
ZWORKS(:,:,KSV+1:KSV+KRR)=PRRS(:,:,:)
ZWORKSFSV(:,:)=0.
ZWORKWSV(:,:,:)=0.
ZWORKSFSV(:,1:KSV)=PSFSV(:,:)
!
!----------------------------------------------------------------------------
!
!*      2. COMPUTE CONSERVATIVE VARIABLES AND RELATED QUANTITIES
!          -----------------------------------------------------
!
!*      2.1 Cph at t
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZCP(:,:)=CST%XCPD
!
IF (KRR > 0) ZCP(:,:) = ZCP(:,:) + CST%XCPV * PRT(:,:,1)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
DO JRR = 2,1+KRRL                          ! loop on the liquid components
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)  
  ZCP(:,:)  = ZCP(:,:) + CST%XCL * PRT(:,:,JRR)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END DO
!
DO JRR = 2+KRRL,1+KRRL+KRRI                ! loop on the solid components   
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZCP(:,:)  = ZCP(:,:)  + CST%XCI * PRT(:,:,JRR)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END DO
!
!*      2.2 Exner function at t
!
IF (OOCEAN) THEN
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZEXN(:,:) = 1.
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZEXN(:,:) = (PPABST(:,:)/CST%XP00) ** (CST%XRD/CST%XCPD)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
!*      2.3 dissipative heating coeff a t
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZCOEF_DISS(:,:) = 1/(ZCP(:,:) * ZEXN(:,:))
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
!
ZFRAC_ICE(:,:) = 0.0
ZATHETA(:,:) = 0.0
ZAMOIST(:,:) = 0.0
!
IF (KRRL >=1) THEN
!
!*      2.4 Temperature at t
!
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZT(:,:) =  PTHLT(:,:) * ZEXN(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
!*       2.5 Lv/Cph/Exn
!
  IF ( KRRI >= 1 ) THEN
    IF (NEBN%LSTATNW) THEN
    !wc call new functions depending on statnew
       CALL COMPUTE_FUNCTION_THERMO_NEW_STAT(CST%XALPW,CST%XBETAW,CST%XGAMW,CST%XLVTT,CST%XCL,ZT,ZEXN,ZCP, &
                                 ZLVOCPEXNM,ZAMOIST,ZATHETA)
       CALL COMPUTE_FUNCTION_THERMO_NEW_STAT(CST%XALPI,CST%XBETAI,CST%XGAMI,CST%XLSTT,CST%XCI,ZT,ZEXN,ZCP, &
                                 ZLSOCPEXNM,ZAMOIST_ICE,ZATHETA_ICE)
    ELSE
      CALL COMPUTE_FUNCTION_THERMO(CST%XALPW,CST%XBETAW,CST%XGAMW,CST%XLVTT,CST%XCL,ZT,ZEXN,ZCP, &
                                 ZLVOCPEXNM,ZAMOIST,ZATHETA)
      CALL COMPUTE_FUNCTION_THERMO(CST%XALPI,CST%XBETAI,CST%XGAMI,CST%XLSTT,CST%XCI,ZT,ZEXN,ZCP, &
                                 ZLSOCPEXNM,ZAMOIST_ICE,ZATHETA_ICE)
    ENDIF
!
    !$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
    WHERE(PRT(:,:,2)+PRT(:,:,4)>0.0)
      ZFRAC_ICE(:,:) = PRT(:,:,4) / ( PRT(:,:,2) &
                                          +PRT(:,:,4) )
    END WHERE
    !$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
!
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZLOCPEXNM(:,:) = (1.0-ZFRAC_ICE(:,:))*ZLVOCPEXNM(:,:) &
                           +ZFRAC_ICE(:,:) *ZLSOCPEXNM(:,:)
    ZAMOIST(:,:) = (1.0-ZFRAC_ICE(:,:))*ZAMOIST(:,:) &
                         +ZFRAC_ICE(:,:) *ZAMOIST_ICE(:,:)
    ZATHETA(:,:) = (1.0-ZFRAC_ICE(:,:))*ZATHETA(:,:) &
                         +ZFRAC_ICE(:,:) *ZATHETA_ICE(:,:)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ELSE
    !wc call new stat functions or not
    IF (NEBN%LSTATNW) THEN
      CALL COMPUTE_FUNCTION_THERMO_NEW_STAT(CST%XALPW,CST%XBETAW,CST%XGAMW,CST%XLVTT,CST%XCL,ZT,ZEXN,ZCP, &
                                 ZLOCPEXNM,ZAMOIST,ZATHETA)
    ELSE
      CALL COMPUTE_FUNCTION_THERMO(CST%XALPW,CST%XBETAW,CST%XGAMW,CST%XLVTT,CST%XCL,ZT,ZEXN,ZCP, &
                                   ZLOCPEXNM,ZAMOIST,ZATHETA)
    ENDIF
  END IF
!
!
  IF ( TPFILE%LOPENED .AND. TURBN%LTURB_DIAG ) THEN
    TZFIELD = TFIELDMETADATA(      &
      CMNHNAME   = 'ATHETA',       &
      CSTDNAME   = '',             &
      CLONGNAME  = 'ATHETA',       &
      CUNITS     = 'm',            &
      CDIR       = 'XY',           &
      CCOMMENT   = 'X_Y_Z_ATHETA', &
      NGRID      = 1,              &
      NTYPE      = TYPEREAL,       &
      NDIMS      = 3,              &
      LTIMEDEP   = .TRUE.          )
    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZATHETA)
!
    TZFIELD = TFIELDMETADATA(      &
      CMNHNAME   = 'AMOIST',       &
      CSTDNAME   = '',             &
      CLONGNAME  = 'AMOIST',       &
      CUNITS     = 'm',            &
      CDIR       = 'XY',           &
      CCOMMENT   = 'X_Y_Z_AMOIST', &
      NGRID      = 1,              &
      NTYPE      = TYPEREAL,       &
      NDIMS      = 3,              &
      LTIMEDEP   = .TRUE.          )
    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZAMOIST)
  END IF
!
ELSE
  ZLOCPEXNM(:,:)=0.
END IF              ! loop end on KRRL >= 1
!
! computes conservative variables
!
IF ( KRRL >= 1 ) THEN
  IF ( KRRI >= 1 ) THEN
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ! Rnp at t
    PRT(:,:,1)  = PRT(:,:,1)  + PRT(:,:,2)  & 
                                    + PRT(:,:,4)
    PRRS(:,:,1) = PRRS(:,:,1) + PRRS(:,:,2) & 
                                    + PRRS(:,:,4)
    ! Theta_l at t
    PTHLT(:,:)  = PTHLT(:,:)  - ZLVOCPEXNM(:,:) &
                                    * PRT(:,:,2) &
                                  - ZLSOCPEXNM(:,:) * PRT(:,:,4)
    PRTHLS(:,:) = PRTHLS(:,:) - ZLVOCPEXNM(:,:) &
                                    * PRRS(:,:,2) &
                                  - ZLSOCPEXNM(:,:) * PRRS(:,:,4)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ELSE
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ! Rnp at t
    PRT(:,:,1)  = PRT(:,:,1)  + PRT(:,:,2)
    PRRS(:,:,1) = PRRS(:,:,1) + PRRS(:,:,2)
    ! Theta_l at t
    PTHLT(:,:)  = PTHLT(:,:)  - ZLOCPEXNM(:,:) &
                                    * PRT(:,:,2)
    PRTHLS(:,:) = PRTHLS(:,:) - ZLOCPEXNM(:,:) &
                                    * PRRS(:,:,2)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
END IF
!
!* stores value of conservative variables & wind before turbulence tendency (AROME diag)
IF(PRESENT(PDRUS_TURB)) THEN
  PDRUS_TURB(:,:) = PRUS(:,:)
  PDRVS_TURB(:,:) = PRVS(:,:)
  PDRTHLS_TURB(:,:) = PRTHLS(:,:)
  PDRRTS_TURB(:,:)  = PRRS(:,:,1)
  PDRSVS_TURB(:,:,:)  = PRSVS(:,:,:)
END IF
!----------------------------------------------------------------------------
!
!*      3. MIXING LENGTH : SELECTION AND COMPUTATION
!          -----------------------------------------
!
!
IF (.NOT. TURBN%LHARAT) THEN

SELECT CASE (TURBN%CTURBLEN)
!
!*      3.1 BL89 mixing length
!           ------------------

  CASE ('BL89')
    ZSHEAR(:,:)=0.
    CALL BL89(D,CST,CSTURB,TURBN,PZZ,PDZZ,PTHVREF,ZTHLM,KRR,ZRM,PTKET,ZSHEAR,ZLM,OOCEAN)
!
!*      3.2 RM17 mixing length
!           ------------------

  CASE ('RM17')
    CALL GZ_U_UW_PHY(D,PUT,PDZZ,ZWORK1)
    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    CALL MXF_PHY(D,ZWORK2,ZDUDZ)
    !
    CALL GZ_V_VW_PHY(D,PVT,PDZZ,ZWORK1)
    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    CALL MYF_PHY(D,ZWORK2,ZDVDZ)
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZSHEAR(:,:) = SQRT(ZDUDZ(:,:)*ZDUDZ(:,:) &
                                    + ZDVDZ(:,:)*ZDVDZ(:,:))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL BL89(D,CST,CSTURB,TURBN,PZZ,PDZZ,PTHVREF,ZTHLM,KRR,ZRM,PTKET,ZSHEAR,ZLM,OOCEAN)
!
!*      3.3 Grey-zone combined RM17 & Deardorff mixing lengths
!           --------------------------------------------------

  CASE ('HM21')
    CALL GZ_U_UW_PHY(D,PUT,PDZZ,ZWORK1)
    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    CALL MXF_PHY(D,ZWORK2,ZDUDZ)
    !
    CALL GZ_V_VW_PHY(D,PVT,PDZZ,ZWORK1)
    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    CALL MYF_PHY(D,ZWORK2,ZDVDZ)
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZSHEAR(:,:) = SQRT(ZDUDZ(:,:)*ZDUDZ(:,:) &
                                    + ZDVDZ(:,:)*ZDVDZ(:,:))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL BL89(D,CST,CSTURB,TURBN,PZZ,PDZZ,PTHVREF,ZTHLM,KRR,ZRM,PTKET,ZSHEAR,ZLM,OOCEAN)

    CALL DELT(ZLMW,ODZ=.FALSE.)
    ! The minimum mixing length is chosen between Horizontal grid mesh (not taking into account the vertical grid mesh) and RM17.
    ! For large horizontal grid meshes, this is equal to RM17
    ! For LES grid meshes, this is equivalent to Deardorff : the base mixing lentgh is the horizontal grid mesh,
    !and it is limited by a stability-based length (RM17), as was done in Deardorff length (but taking into account shear as well)
    ! For grid meshes in the grey zone, then this is the smaller of the two.
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZLM(:,:) = MIN(ZLM(:,:),TURBN%XCADAP*ZLMW(:,:))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
!*      3.4 Delta mixing length
!           -------------------
!
  CASE ('DELT')
    CALL DELT(ZLM,ODZ=.TRUE.)
!
!*      3.5 Deardorff mixing length
!           -----------------------
!
  CASE ('DEAR')
    CALL DEAR(ZLM)
!
!*      3.6 Blackadar mixing length
!           -----------------------
!
  CASE ('BLKR')
   ZL0 = 100.
   !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
   ZLM(:,:) = ZL0
   !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

   ZALPHA=0.5**(-1.5)
   !
   DO JK=IKTB,IKTE
     !$mnh_expand_array(JIJ=IIJB:IIJE)
     ZLM(:,JK) = ( 0.5*(PZZ(:,JK)+PZZ(:,JK+IKL)) - &
     & PZZ(:,IKA+JPVEXT_TURB*IKL) ) * PDIRCOSZW(:)
     ZLM(:,JK) = ZALPHA  * ZLM(:,JK) * ZL0 / ( ZL0 + ZALPHA*ZLM(:,JK) )
     !$mnh_end_expand_array(JIJ=IIJB:IIJE)
   END DO
!
   !$mnh_expand_array(JIJ=IIJB:IIJE)
   ZLM(:,IKTB-1) = ZLM(:,IKTB)
   ZLM(:,IKTE+1) = ZLM(:,IKTE)
   !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
!
!
END SELECT
!
!*      3.5 Mixing length modification for cloud
!           -----------------------
IF (OCLOUDMODIFLM) CALL CLOUD_MODIF_LM
ENDIF  ! end LHARRAT

!
!*      3.6 Dissipative length
!           ------------------

IF (TURBN%LHARAT) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZLEPS(:,:)=PLENGTHM(:,:)*(3.75**2.)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  ZLEPS(:,:)=ZLM(:,:)
ENDIF
!
!*      3.7 Correction in the Surface Boundary Layer (Redelsperger 2001)
!           ----------------------------------------
!
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZLMO(:)=XUNDEF
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
IF (TURBN%LRMC01) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZUSTAR(:)=(PSFU(:)**2+PSFV(:)**2)**(0.25)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  IF (KRR>0) THEN
    CALL LMO(D,CST,ZUSTAR,ZTHLM(:,IKB),ZRM(:,IKB,1),PSFTH,PSFRV,ZLMO)
  ELSE
    ZRVM(:)=0.
    ZSFRV(:)=0.
    CALL LMO(D,CST,ZUSTAR,ZTHLM(:,IKB),ZRVM,PSFTH,ZSFRV,ZLMO)
  END IF
  CALL RMC01(D,CST,CSTURB,TURBN,PZZ,PDXX,PDYY,PDZZ,PDIRCOSZW,PSBL_DEPTH,ZLMO,ZLM,ZLEPS)
END IF
!
!RMC01 is only applied on RM17 in HM21
IF (TURBN%CTURBLEN=='HM21') THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZLEPS(:,:) = MIN(ZLEPS(:,:),ZLMW(:,:)*TURBN%XCADAP)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
!*      3.8 Mixing length in external points (used if TURBN%CTURBDIM="3DIM")
!           ----------------------------------------------------------
!
IF (TURBN%CTURBDIM=="3DIM") THEN
  CALL UPDATE_LM(D,HLBCX,HLBCY,ZLM,ZLEPS)
END IF
!
!*      3.9 Mixing length correction if immersed walls
!           ------------------------------------------
!
IF (OIBM) THEN
   CALL IBM_MIXINGLENGTH(D,ZLM,ZLEPS,PIBM_XMUT,PIBM_LS,PTKET)
ENDIF
!----------------------------------------------------------------------------
!
!*      4. GO INTO THE AXES FOLLOWING THE SURFACE
!          --------------------------------------
!
!
!*      4.1 rotate the wind at time t
!
!
!
IF (TURBN%LROTATE_WIND) THEN
  CALL ROTATE_WIND(D,PUT,PVT,PWT,                       &
                     PDIRCOSXW, PDIRCOSYW, PDIRCOSZW,   &
                     PCOSSLOPE,PSINSLOPE,               &
                     PDXX,PDYY,PDZZ,                    &
                     ZUSLOPE,ZVSLOPE                    )
!
  CALL UPDATE_ROTATE_WIND(D,ZUSLOPE,ZVSLOPE,HLBCX,HLBCY)
ELSE
  ZUSLOPE(:)=PUT(:,IKA)
  ZVSLOPE(:)=PVT(:,IKA)
END IF
IF (OOCEAN) THEN
  ZUSLOPE(:)=PUT(:,IKU-1)
  ZVSLOPE(:)=PVT(:,IKU-1)
END IF
!
!
!*      4.2 compute the proportionality coefficient between wind and stress
!
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZCDUEFF(:) =-SQRT ( (PSFU(:)**2 + PSFV(:)**2) /               &
                    (CST%XMNH_TINY + ZUSLOPE(:)**2 + ZVSLOPE(:)**2 ) )
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
!*       4.6 compute the surface tangential fluxes
!
IF (OOCEAN) THEN
  ZTAU11M(:)=0.
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE)        
  ZTAU11M(:) =2./3.*(  (1.+ (PZZ(:,IKB+IKL)-PZZ(:,IKB))  &
                           /(PDZZ(:,IKB+IKL)+PDZZ(:,IKB))  &
                       )   *PTKET(:,IKB)                   &
                     -0.5  *PTKET(:,IKB+IKL)                 &
                    )
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END IF
ZTAU12M(:) =0.0
ZTAU22M(:) =ZTAU11M(:)
ZTAU33M(:) =ZTAU11M(:)
!
!*       4.7 third order terms in temperature and water fluxes and correlations
!            ------------------------------------------------------------------
!
!
ZMWTH(:,:) = 0.     ! w'2th'
ZMWR(:,:)  = 0.     ! w'2r'
ZMTH2(:,:) = 0.     ! w'th'2
ZMR2(:,:)  = 0.     ! w'r'2
ZMTHR(:,:) = 0.     ! w'th'r'
!
IF (TURBN%CTOM=='TM06') THEN
  CALL TM06(D,CST,PTHVREF,PBL_DEPTH,PZZ,PSFTH,ZMWTH,ZMTH2)
!
   CALL GZ_M_W_PHY(D,ZMWTH,PDZZ,ZWORK1)    ! -d(w'2th' )/dz
   CALL GZ_W_M_PHY(D,ZMTH2,PDZZ,ZWORK2)    ! -d(w'th'2 )/dz
   !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
   ZFWTH(:,:) = -ZWORK1(:,:)
   ZFTH2(:,:) = -ZWORK2(:,:)
   !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
  ZFWTH(:,IKTE:) = 0.
  ZFWTH(:,:IKTB) = 0.
  ZFWR(:,:)  = 0.
  ZFTH2(:,IKTE:) = 0.
  ZFTH2(:,:IKTB) = 0.
  ZFR2(:,:)  = 0.
  ZFTHR(:,:) = 0.
ELSE
  ZFWTH(:,:) = 0.
  ZFWR(:,:)  = 0.
  ZFTH2(:,:) = 0.
  ZFR2(:,:)  = 0.
  ZFTHR(:,:) = 0.
ENDIF
!
!----------------------------------------------------------------------------
!
!*      5. TURBULENT SOURCES
!          -----------------
!
IF( BUCONF%LBUDGET_U )  CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_U ), 'VTURB', PRUS(:,:)    )
IF( BUCONF%LBUDGET_V )  CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_V ), 'VTURB', PRVS(:,:)    )
IF( BUCONF%LBUDGET_W )  CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_W ), 'VTURB', PRWS(:,:)    )

IF( BUCONF%LBUDGET_TH ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:,:) + ZLVOCPEXNM(:,:) * PRRS(:,:, 2) &
                                                                          + ZLSOCPEXNM(:,:) * PRRS(:,:, 4) )
  ELSE IF( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:,:) + ZLOCPEXNM(:,:) * PRRS(:,:, 2) )
  ELSE
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:,:) )
  END IF
END IF

IF( BUCONF%LBUDGET_RV ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:,:, 1) - PRRS(:,:, 2) - PRRS(:,:, 4) )
  ELSE IF( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:,:, 1) - PRRS(:,:, 2) )
  ELSE
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:,:, 1) )
  END IF
END IF

IF( BUCONF%LBUDGET_RC ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'VTURB', PRRS  (:,:, 2) )
IF( BUCONF%LBUDGET_RR ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'VTURB', ZWORKS  (:,:, KSV + 3) )
IF( BUCONF%LBUDGET_RI ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'VTURB', PRRS  (:,:, 4) )
IF( BUCONF%LBUDGET_RS ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'VTURB', ZWORKS  (:,:, KSV + 5) )
IF( BUCONF%LBUDGET_RG ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'VTURB', ZWORKS  (:,:, KSV + 6) )
IF( BUCONF%LBUDGET_RH ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RH), 'VTURB', ZWORKS  (:,:, KSV + 7) )

IF( BUCONF%LBUDGET_SV ) THEN
  DO JSV = 1, KSV
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + JSV), 'VTURB', ZWORKS(:,:, JSV) )
  END DO
END IF

CALL TURB_VER(D,CST,CSTURB,TURBN,NEBN,TLES,              &
          KRR,KRRL,KRRI,KGRADIENTS,                      &
          OOCEAN, ODEEPOC, OCOMPUTE_SRC,                 &
          KSV+KRR,KSV_LGBEG,KSV_LGEND,                   &
          ZEXPL, O2D, ONOMIXLG, OFLAT,                   &
          OCOUPLES,OBLOWSNOW,OFLYER, PRSNOW,             &
          PTSTEP,TPFILE,                                 &
          PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,        &
          PCOSSLOPE,PSINSLOPE,                           &
          PRHODJ,PTHVREF,PSFU,PSFV,                      &
          PSFTH,PSFRV,ZWORKSFSV,PSFTH,PSFRV,ZWORKSFSV,   &
          ZCDUEFF,ZTAU11M,ZTAU12M,ZTAU33M,               &
          PUT,PVT,PWT,ZUSLOPE,ZVSLOPE,PTHLT,PRT,ZWORKT,  &
          PTKET,ZLM,PLENGTHM,PLENGTHH,ZLEPS,MFMOIST,     &
          ZLOCPEXNM,ZATHETA,ZAMOIST,PSRCT,ZFRAC_ICE,     &
          ZFWTH,ZFWR,ZFTH2,ZFR2,ZFTHR,PBL_DEPTH,         &
          PSBL_DEPTH,ZLMO,PHGRAD,PZS,                    &
          PRUS,PRVS,PRWS,PRTHLS,PRRS,ZWORKS,             &
          PDP,PTP,PSIGS,PWTH,PWRC,ZWORKWSV,                  &
          PSSTFL, PSSTFL_C, PSSRFL_C,PSSUFL_C,PSSVFL_C,  &
          PSSUFL,PSSVFL                                  )

!IF (HCLOUD == 'LIMA') THEN
!   IF (KSV_LIMA_NR.GT.0) PRSVS(:,:,KSV_LIMA_NR) = ZRSVS(:,:,KSV_LIMA_NR) 
!   IF (KSV_LIMA_NS.GT.0) PRSVS(:,:,KSV_LIMA_NS) = ZRSVS(:,:,KSV_LIMA_NS)
!   IF (KSV_LIMA_NG.GT.0) PRSVS(:,:,KSV_LIMA_NG) = ZRSVS(:,:,KSV_LIMA_NG) 
!   IF (KSV_LIMA_NH.GT.0) PRSVS(:,:,KSV_LIMA_NH) = ZRSVS(:,:,KSV_LIMA_NH)
!END IF

IF( BUCONF%LBUDGET_U ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_U), 'VTURB', PRUS(:,:) )
IF( BUCONF%LBUDGET_V ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_V), 'VTURB', PRVS(:,:) )
IF( BUCONF%LBUDGET_W ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_W), 'VTURB', PRWS(:,:) )

IF( BUCONF%LBUDGET_TH ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:,:) + ZLVOCPEXNM(:,:) * PRRS(:,:, 2) &
                                                                          + ZLSOCPEXNM(:,:) * PRRS(:,:, 4) )
  ELSE IF( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:,:) + ZLOCPEXNM(:,:) * PRRS(:,:, 2) )
  ELSE
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:,:) )
  END IF
END IF

IF( BUCONF%LBUDGET_RV ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:,:, 1) - PRRS(:,:, 2) - PRRS(:,:, 4) )
  ELSE IF( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:,:, 1) - PRRS(:,:, 2) )
  ELSE
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:,:, 1) )
  END IF
END IF

IF( BUCONF%LBUDGET_RC ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'VTURB', PRRS(:,:, 2) )
IF( BUCONF%LBUDGET_RR ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'VTURB', ZWORKS(:,:, KSV + 3) )
IF( BUCONF%LBUDGET_RI ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'VTURB', PRRS(:,:, 4) )
IF( BUCONF%LBUDGET_RS ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'VTURB', ZWORKS(:,:, KSV + 5) )
IF( BUCONF%LBUDGET_RG ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'VTURB', ZWORKS(:,:, KSV + 6) )
IF( BUCONF%LBUDGET_RH ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RH), 'VTURB', ZWORKS(:,:, KSV + 7) )

IF( BUCONF%LBUDGET_SV )  THEN
  DO JSV = 1, KSV
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + JSV), 'VTURB', ZWORKS(:,:, JSV) )
  END DO
END IF
!
IF( TURBN%CTURBDIM == '3DIM' ) THEN
  IF( BUCONF%LBUDGET_U  ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_U ), 'HTURB', PRUS  (:,:) )
  IF( BUCONF%LBUDGET_V  ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_V ), 'HTURB', PRVS  (:,:) )
  IF( BUCONF%LBUDGET_W  ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_W ), 'HTURB', PRWS  (:,:) )

  IF(BUCONF%LBUDGET_TH)  THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:,:) + ZLVOCPEXNM(:,:) * PRRS(:,:, 2) &
                                                                             + ZLSOCPEXNM(:,:) * PRRS(:,:, 4) )
    ELSE IF( KRRL >= 1 ) THEN
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:,:) + ZLOCPEXNM(:,:) * PRRS(:,:, 2) )
    ELSE
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:,:) )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RV ) THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:,:, 1) - PRRS(:,:, 2) - PRRS(:,:, 4) )
    ELSE IF( KRRL >= 1 ) THEN
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:,:, 1) - PRRS(:,:, 2) )
    ELSE
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:,:, 1) )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RC ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'HTURB', PRRS(:,:, 2) )
  IF( BUCONF%LBUDGET_RR ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'HTURB', ZWORKS(:,:, KSV+3) )
  IF( BUCONF%LBUDGET_RI ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'HTURB', PRRS(:,:, 4) )
  IF( BUCONF%LBUDGET_RS ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'HTURB', ZWORKS(:,:, KSV+5) )
  IF( BUCONF%LBUDGET_RG ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'HTURB', ZWORKS(:,:, KSV+6) )
  IF( BUCONF%LBUDGET_RH ) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RH), 'HTURB', ZWORKS(:,:, KSV+7) )

  IF( BUCONF%LBUDGET_SV )  THEN
    DO JSV = 1, KSV
      CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + JSV), 'HTURB', ZWORKS(:,:, JSV) )
    END DO
  END IF
    CALL TURB_HOR_SPLT(D,CST,CSTURB, TURBN, NEBN, TLES,        &
          KSPLIT, KRR, KRRL, KRRI, KSV,KSV_LGBEG,KSV_LGEND,    & 
          KSPLIT, KRR, KRRL, KRRI, KSV+KRR,KSV_LGBEG,KSV_LGEND,& 
          PTSTEP,HLBCX,HLBCY, OFLAT,O2D, ONOMIXLG,             & 
          OOCEAN,OCOMPUTE_SRC,OBLOWSNOW,PRSNOW,                &
          TPFILE, KHALO,                                       &
          PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,                        &
          PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,                       &
          PCOSSLOPE,PSINSLOPE,                                 &
          PRHODJ,PTHVREF,                                      &
          PSFTH,PSFRV,ZWORKSFSV,                               &
          ZCDUEFF,ZTAU11M,ZTAU12M,ZTAU22M,ZTAU33M,             &
          PUT,PVT,PWT,ZUSLOPE,ZVSLOPE,PTHLT,PRT,ZWORKT,        &
          PTKET,ZLM,ZLEPS,                                     &
          ZLOCPEXNM,ZATHETA,ZAMOIST,PSRCT,ZFRAC_ICE,           &
          PDP,PTP,PSIGS,                                       &
          ZTRH,                                                &
          PRUS,PRVS,PRWS,PRTHLS,PRRS,ZWORKS                    )
  !
!  IF (HCLOUD == 'LIMA') THEN
!     IF (KSV_LIMA_NR.GT.0) PRSVS(:,:,KSV_LIMA_NR) = ZRSVS(:,:,KSV_LIMA_NR) 
!     IF (KSV_LIMA_NS.GT.0) PRSVS(:,:,KSV_LIMA_NS) = ZRSVS(:,:,KSV_LIMA_NS)
!     IF (KSV_LIMA_NG.GT.0) PRSVS(:,:,KSV_LIMA_NG) = ZRSVS(:,:,KSV_LIMA_NG) 
!     IF (KSV_LIMA_NH.GT.0) PRSVS(:,:,KSV_LIMA_NH) = ZRSVS(:,:,KSV_LIMA_NH)
!  END IF
  !
  IF( BUCONF%LBUDGET_U ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_U), 'HTURB', PRUS(:,:) )
  IF( BUCONF%LBUDGET_V ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_V), 'HTURB', PRVS(:,:) )
  IF( BUCONF%LBUDGET_W ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_W), 'HTURB', PRWS(:,:) )

  IF( BUCONF%LBUDGET_TH ) THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:,:) + ZLVOCPEXNM(:,:) * PRRS(:,:, 2) &
                                                                            + ZLSOCPEXNM(:,:) * PRRS(:,:, 4) )
    ELSE IF( KRRL >= 1 ) THEN
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:,:) + ZLOCPEXNM(:,:) * PRRS(:,:, 2) )
    ELSE
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:,:) )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RV ) THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:,:, 1) - PRRS(:,:, 2) - PRRS(:,:, 4) )
    ELSE IF( KRRL >= 1 ) THEN
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:,:, 1) - PRRS(:,:, 2) )
    ELSE
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:,:, 1) )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RC ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'HTURB', PRRS(:,:, 2) )
  IF( BUCONF%LBUDGET_RR ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'HTURB', ZWORKS(:,:, KSV+3) )
  IF( BUCONF%LBUDGET_RI ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'HTURB', PRRS(:,:, 4) )
  IF( BUCONF%LBUDGET_RS ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'HTURB', ZWORKS(:,:, KSV+5) )
  IF( BUCONF%LBUDGET_RG ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'HTURB', ZWORKS(:,:, KSV+6) )
  IF( BUCONF%LBUDGET_RH ) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RH), 'HTURB', ZWORKS(:,:, KSV+7) )

  IF( BUCONF%LBUDGET_SV )  THEN
    DO JSV = 1, KSV
      CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + JSV), 'HTURB', ZWORKS(:,:, JSV) )
    END DO
  END IF
END IF
!----------------------------------------------------------------------------
!
!*      6. EVOLUTION OF THE TKE AND ITS DISSIPATION
!          ----------------------------------------
!
!  6.1 Contribution of mass-flux in the TKE buoyancy production if
!      cloud computation is not statistical
CALL MZF_PHY(D,PFLXZTHVMF,ZWORK1)
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PTP(:,:) = PTP(:,:) &
                             + CST%XG / PTHVREF(:,:) * ZWORK1(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

IF(PRESENT(PTPMF))  THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PTPMF(:,:)=CST%XG / PTHVREF(:,:) * ZWORK1(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!  6.2 TKE evolution equation

IF (.NOT. TURBN%LHARAT) THEN
!
IF (BUCONF%LBUDGET_TH)  THEN
  IF ( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS(:,:)+ ZLVOCPEXNM(:,:) * PRRS(:,:,2) &
                                                          & + ZLSOCPEXNM(:,:) * PRRS(:,:,4) )
  ELSE IF ( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS(:,:) + ZLOCPEXNM(:,:) * PRRS(:,:,2) )
  ELSE
    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS(:,:) )
  END IF
END IF
!
IF(PRESENT(PRTKEMS)) THEN
  ZRTKEMS(:,:)=PRTKEMS(:,:)
ELSE
  ZRTKEMS(:,:)=0.
END IF
!
CALL TKE_EPS_SOURCES(D,CST,CSTURB,BUCONF,TURBN,TLES,                    &
                   & PTKET,ZLM,ZLEPS,PDP,ZTRH,                          &
                   & PRHODJ,PDZZ,PDXX,PDYY,PDZX,PDZY,PZZ,               &
                   & PTSTEP,ZEXPL,                                      &
                   & TPFILE,ODIAG_IN_RUN,OOCEAN,                        &
                   & PSFU,PSFV,                                         &
                   & PTP,PRTKES,PRTHLS,ZCOEF_DISS,PTDIFF,PTDISS,ZRTKEMS,&
                   & TBUDGETS,KBUDGETS, PEDR=PEDR, PTR=PTR,PDISS=PDISS, &
                   & PCURRENT_TKE_DISS=PCURRENT_TKE_DISS                )
IF (BUCONF%LBUDGET_TH)  THEN
  IF ( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS(:,:)+ ZLVOCPEXNM(:,:) * PRRS(:,:,2) &
                                                          & + ZLSOCPEXNM(:,:) * PRRS(:,:,4) )
  ELSE IF ( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS(:,:)+ ZLOCPEXNM(:,:) * PRRS(:,:,2) )
  ELSE
    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS(:,:) )
  END IF
END IF
!
ENDIF
!
!----------------------------------------------------------------------------
!
!*      7. STORES SOME INFORMATIONS RELATED TO THE TURBULENCE SCHEME
!          ---------------------------------------------------------
!
IF ( TURBN%LTURB_DIAG .AND. TPFILE%LOPENED ) THEN
!
! stores the mixing length
!
  TZFIELD = TFIELDMETADATA(       &
    CMNHNAME   = 'LM',            &
    CSTDNAME   = '',              &
    CLONGNAME  = 'LM',            &
    CUNITS     = 'm',             &
    CDIR       = 'XY',            &
    CCOMMENT   = 'Mixing length', &
    NGRID      = 1,               &
    NTYPE      = TYPEREAL,        &
    NDIMS      = 3,               &
    LTIMEDEP   = .TRUE.           )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZLM)
!
  IF (KRR /= 0) THEN
!
! stores the conservative potential temperature
!
    TZFIELD = TFIELDMETADATA(                          &
    CMNHNAME   = 'THLM',                               &
    CSTDNAME   = '',                                   &
    CLONGNAME  = 'THLM',                               &
    CUNITS     = 'K',                                  &
    CDIR       = 'XY',                                 &
    CCOMMENT   = 'Conservative potential temperature', &
    NGRID      = 1,                                    &
    NTYPE      = TYPEREAL,                             &
    NDIMS      = 3,                                    &
    LTIMEDEP   = .TRUE.                                )
    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PTHLT)
!
! stores the conservative mixing ratio
!
    TZFIELD = TFIELDMETADATA(                &
    CMNHNAME   = 'RNPM',                     &
    CSTDNAME   = '',                         &
    CLONGNAME  = 'RNPM',                     &
    CUNITS     = 'kg kg-1',                  &
    CDIR       = 'XY',                       &
    CCOMMENT   = 'Conservative mixing ratio',&
    NGRID      = 1,                          &
    NTYPE      = TYPEREAL,                   &
    NDIMS      = 3,                          &
    LTIMEDEP   = .TRUE.                      )
    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PRT(:,:,1))
   END IF
END IF
!
!* stores value of conservative variables & wind before turbulence tendency (AROME only)
IF(PRESENT(PDRUS_TURB)) THEN
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PDRUS_TURB(:,:)   = PRUS(:,:) - PDRUS_TURB(:,:)
  PDRVS_TURB(:,:)   = PRVS(:,:) - PDRVS_TURB(:,:)
  PDRTHLS_TURB(:,:) = PRTHLS(:,:) - PDRTHLS_TURB(:,:)
  PDRRTS_TURB(:,:)  = PRRS(:,:,1) - PDRRTS_TURB(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT,JSV=1:KSV)  
  PDRSVS_TURB(:,:,:)  = PRSVS(:,:,:) - PDRSVS_TURB(:,:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT,JSV=1:KSV)
END IF
!----------------------------------------------------------------------------
!
!*      8. RETRIEVE NON-CONSERVATIVE VARIABLES
!          -----------------------------------
!
IF ( KRRL >= 1 ) THEN
  IF ( KRRI >= 1 ) THEN
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PRT(:,:,1)  = PRT(:,:,1)  - PRT(:,:,2)  &
                                    - PRT(:,:,4)
    PRRS(:,:,1) = PRRS(:,:,1) - PRRS(:,:,2) &
                                    - PRRS(:,:,4)
    PTHLT(:,:)  = PTHLT(:,:)  + ZLVOCPEXNM(:,:) &
                                    * PRT(:,:,2) &
                                    + ZLSOCPEXNM(:,:) * PRT(:,:,4)
    PRTHLS(:,:) = PRTHLS(:,:) + ZLVOCPEXNM(:,:) &
                                    * PRRS(:,:,2) &
                                    + ZLSOCPEXNM(:,:) * PRRS(:,:,4)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
  ELSE
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PRT(:,:,1)  = PRT(:,:,1)  - PRT(:,:,2)
    PRRS(:,:,1) = PRRS(:,:,1) - PRRS(:,:,2)
    PTHLT(:,:)  = PTHLT(:,:)  + ZLOCPEXNM(:,:) &
                                    * PRT(:,:,2)
    PRTHLS(:,:) = PRTHLS(:,:) + ZLOCPEXNM(:,:) &
                                    * PRRS(:,:,2)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
END IF!
!
PRSVS(:,:,:)=ZWORKS(:,:,1:KSV)
IF (TURBN%LTURB_PRECIP) THEN
   IF (KRR.GE.3) PRRS(:,:,3)=ZWORKS(:,:,KSV+3)
   IF (KRR.GE.5) PRRS(:,:,5)=ZWORKS(:,:,KSV+5)
   IF (KRR.GE.6) PRRS(:,:,6)=ZWORKS(:,:,KSV+6)
   IF (KRR.GE.7) PRRS(:,:,7)=ZWORKS(:,:,KSV+7)
END IF
   IF (OFLYER)   PWSV(:,:,:)=ZWORKWSV(:,:,1:KSV)
!
! Remove non-physical negative values (unnecessary in a perfect world) + corresponding budgets
CALL SOURCES_NEG_CORRECT_PHY(D,KSV,HCLOUD,HELEC,'NETUR',KRR,PTSTEP,PPABST,PTHLT,PRT,PRTHLS,PRRS,PRSVS)
!----------------------------------------------------------------------------
!
!*      9. LES averaged surface fluxes
!          ---------------------------
!
IF (TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,PSFTH,TLES%X_LES_Q0)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,PSFRV,TLES%X_LES_E0)
  DO JSV=1,KSV
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,PSFSV(:,JSV),TLES%X_LES_SV0(:,JSV))
  END DO
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,PSFU,TLES%X_LES_UW0)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,PSFV,TLES%X_LES_VW0)
  !
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZWORK2D(:) = (PSFU(:)*PSFU(:)+PSFV(:)*PSFV(:))**0.25
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2D,TLES%X_LES_USTAR)
!----------------------------------------------------------------------------
!
!*     10. LES for 3rd order moments
!          -------------------------
!
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZMWTH,TLES%X_LES_SUBGRID_W2Thl)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZMTH2,TLES%X_LES_SUBGRID_WThl2)
  IF (KRR>0) THEN
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZMWR,TLES%X_LES_SUBGRID_W2Rt)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZMTHR,TLES%X_LES_SUBGRID_WThlRt)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZMR2,TLES%X_LES_SUBGRID_WRt2)
  END IF
!
!----------------------------------------------------------------------------
!
!*     11. LES quantities depending on <w'2> in "1DIM" mode
!          ------------------------------------------------
!
  IF (TURBN%CTURBDIM=="1DIM") THEN
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK1(:,:) = 2./3.*PTKET(:,:)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK1,TLES%X_LES_SUBGRID_U2)
    TLES%X_LES_SUBGRID_V2(:,:,:) = TLES%X_LES_SUBGRID_U2(:,:,:)
    TLES%X_LES_SUBGRID_W2(:,:,:) = TLES%X_LES_SUBGRID_U2(:,:,:)
    !
    CALL GZ_M_W_PHY(D,PTHLT,PDZZ,ZWORK1)
    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK2(:,:)  = 2./3.*PTKET(:,:) *ZWORK2(:,:)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2,TLES%X_LES_RES_ddz_Thl_SBG_W2)
    !
    IF (KRR>=1) THEN
      CALL GZ_M_W_PHY(D,PRT(:,:,1),PDZZ,ZWORK1)
      CALL MZF_PHY(D,ZWORK1,ZWORK2)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZWORK2(:,:)  = 2./3.*PTKET(:,:) *ZWORK2(:,:)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2,TLES%X_LES_RES_ddz_Rt_SBG_W2)
    END IF
    DO JSV=1,KSV
      CALL GZ_M_W_PHY(D,PSVT(:,:,JSV),PDZZ,ZWORK1)
      CALL MZF_PHY(D,ZWORK1,ZWORK2)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZWORK2(:,:)  = 2./3.*PTKET(:,:) *ZWORK2(:,:)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2, TLES%X_LES_RES_ddz_Sv_SBG_W2(:,:,:,JSV))
    END DO
  END IF

!----------------------------------------------------------------------------
!
!*     12. LES mixing end dissipative lengths, presso-correlations
!          -------------------------------------------------------
!
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZLM,TLES%X_LES_SUBGRID_LMix)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZLEPS,TLES%X_LES_SUBGRID_LDiss)
!
!* presso-correlations for subgrid Tke are equal to zero.
!
  ZLEPS(:,:) = 0. !ZLEPS is used as a work array (not used anymore)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZLEPS,TLES%X_LES_SUBGRID_WP)
!
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
IF(PRESENT(PLEM)) PLEM(:,:) = ZLM(:,:)
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TURB',1,ZHOOK_HANDLE)
CONTAINS
!
!     ########################################################################
      SUBROUTINE COMPUTE_FUNCTION_THERMO(PALP,PBETA,PGAM,PLTT,PC,PT,PEXN,PCP,&
                                         PLOCPEXN,PAMOIST,PATHETA            )
!     ########################################################################
!!
!!****  *COMPUTE_FUNCTION_THERMO* routine to compute several thermo functions
!
!!    AUTHOR
!!    ------
!!
!!     JP Pinty      *LA*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   24/02/03
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments
!
REAL,                   INTENT(IN)    :: PALP,PBETA,PGAM,PLTT,PC
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PT,PEXN,PCP
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PLOCPEXN
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PAMOIST,PATHETA
!
!-------------------------------------------------------------------------------
!
  IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO',0,ZHOOK_HANDLE2)
  ZEPS = CST%XMV / CST%XMD
!
!*       1.1 Lv/Cph at  t
!
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PLOCPEXN(:,:) = ( PLTT + (CST%XCPV-PC) *  (PT(:,:)-CST%XTT) ) &
                                     / PCP(:,:)
!
!*      1.2 Saturation vapor pressure at t
!
  ZRVSAT(:,:) =  EXP( PALP - PBETA/PT(:,:) - PGAM*ALOG( PT(:,:) ) )
!
!*      1.3 saturation  mixing ratio at t
!
  ZRVSAT(:,:) =  ZRVSAT(:,:) &
                                    * ZEPS / ( PPABST(:,:) - ZRVSAT(:,:) )
!
!*      1.4 compute the saturation mixing ratio derivative (rvs')
!
  ZDRVSATDT(:,:) = ( PBETA / PT(:,:)  - PGAM ) / PT(:,:)   &
                 * ZRVSAT(:,:) * ( 1. + ZRVSAT(:,:) / ZEPS )
!
!*      1.5 compute Amoist
!
  PAMOIST(:,:)=  0.5 / ( 1.0 + ZDRVSATDT(:,:) * PLOCPEXN(:,:) )
!
!*      1.6 compute Atheta
!
  PATHETA(:,:)= PAMOIST(:,:) * PEXN(:,:) *               &
        ( ( ZRVSAT(:,:) - PRT(:,:,1) ) * PLOCPEXN(:,:) / &
          ( 1. + ZDRVSATDT(:,:) * PLOCPEXN(:,:) )        *               &
          (                                                                  &
           ZRVSAT(:,:) * (1. + ZRVSAT(:,:)/ZEPS)                         &
                        * ( -2.*PBETA/PT(:,:) + PGAM ) / PT(:,:)**2      &
          +ZDRVSATDT(:,:) * (1. + 2. * ZRVSAT(:,:)/ZEPS)                 &
                        * ( PBETA/PT(:,:) - PGAM ) / PT(:,:)             &
          )                                                                  &
         - ZDRVSATDT(:,:)                                                  &
        )
!
!*      1.7 Lv/Cph/Exner at t-1
!
  PLOCPEXN(:,:) = PLOCPEXN(:,:) / PEXN(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO',1,ZHOOK_HANDLE2)
END SUBROUTINE COMPUTE_FUNCTION_THERMO

!     ########################################################################
      SUBROUTINE COMPUTE_FUNCTION_THERMO_NEW_STAT(PALP,PBETA,PGAM,PLTT,PC,PT,PEXN,PCP,&
                                         PLOCPEXN,PAMOIST,PATHETA            )
!     ########################################################################
!!
!!****  *COMPUTE_FUNCTION_THERMO* routine to compute several thermo functions
!
!!    AUTHOR
!!    ------
!!
!!     JP Pinty      *LA*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   24/02/03
!!     Modified: Wim de Rooy 06-02-2019
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments
!
REAL, INTENT(IN)                      :: PALP,PBETA,PGAM,PLTT,PC
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PT,PEXN,PCP
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PLOCPEXN
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PAMOIST,PATHETA
!
!-------------------------------------------------------------------------------
!
  IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO_NEW_STAT',0,ZHOOK_HANDLE2)
  ZEPS = CST%XMV / CST%XMD
!
!*       1.1 Lv/Cph at  t
!
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PLOCPEXN(:,:) = ( PLTT + (CST%XCPV-PC) *  (PT(:,:)-CST%XTT) ) / PCP(:,:)
!
!*      1.2 Saturation vapor pressure at t
!
  ZRVSAT(:,:) =  EXP( PALP - PBETA/PT(:,:) - PGAM*ALOG( PT(:,:) ) )
!
!*      1.3 saturation  mixing ratio at t
!
  ZRVSAT(:,:) =  ZRVSAT(:,:) * ZEPS / ( PPABST(:,:) - ZRVSAT(:,:) )
!
!*      1.4 compute the saturation mixing ratio derivative (rvs')
!
  ZDRVSATDT(:,:) = ( PBETA / PT(:,:)  - PGAM ) / PT(:,:)   &
                 * ZRVSAT(:,:) * ( 1. + ZRVSAT(:,:) / ZEPS )
!
!*      1.5 compute Amoist
!
  PAMOIST(:,:)=  1.0 / ( 1.0 + ZDRVSATDT(:,:) * PLOCPEXN(:,:) )
!
!*      1.6 compute Atheta
!
  PATHETA(:,:)= PAMOIST(:,:) * PEXN(:,:) * ZDRVSATDT(:,:)
!
!*      1.7 Lv/Cph/Exner at t-1
!
  PLOCPEXN(:,:) = PLOCPEXN(:,:) / PEXN(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO_NEW_STAT',1,ZHOOK_HANDLE2)
END SUBROUTINE COMPUTE_FUNCTION_THERMO_NEW_STAT

!
!     ####################
      SUBROUTINE DELT(PLM,ODZ)
!     ####################
!!
!!****  *DELT* routine to compute mixing length for DELT case
!
!!    AUTHOR
!!    ------
!!
!!     M Tomasini      *Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   01/05
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!*       0.1   Declarations of dummy arguments
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PLM
LOGICAL,                INTENT(IN)    :: ODZ
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TURB:DELT',0,ZHOOK_HANDLE2)
!
CALL MXF_PHY(D,PDXX,ZWORK1)
IF (.NOT. O2D) THEN
  CALL MYF_PHY(D,PDYY,ZWORK2)
END IF
!
IF (ODZ) THEN
  ! Dz is take into account in the computation
  DO JK = IKTB,IKTE ! 1D turbulence scheme
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    PLM(:,JK) = PZZ(:,JK+IKL) - PZZ(:,JK)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  END DO
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  PLM(:,IKU) = PLM(:,IKE)
  PLM(:,IKA) = PZZ(:,IKB) - PZZ(:,IKA)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  IF ( TURBN%CTURBDIM /= '1DIM' ) THEN  ! 3D turbulence scheme
    IF ( O2D) THEN
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      PLM(:,:) = SQRT( PLM(:,:)*ZWORK1(:,:) )
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ELSE
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      PLM(:,:) = (PLM(:,:)*ZWORK1(:,:) &
                                   * ZWORK2(:,:) ) ** (1./3.)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    END IF
  END IF
ELSE
  ! Dz not taken into account in computation to assure invariability with vertical grid mesh
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PLM(:,:)=1.E10
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  IF ( TURBN%CTURBDIM /= '1DIM' ) THEN  ! 3D turbulence scheme
    IF ( O2D) THEN
      PLM(:,:) = ZWORK1(:,:)
    ELSE
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      PLM(:,:) = (ZWORK1(:,:)*ZWORK2(:,:) ) ** (1./2.)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    END IF
  END IF
END IF
!
!  mixing length limited by the distance normal to the surface
!  (with the same factor as for BL89)
!
IF (.NOT. TURBN%LRMC01) THEN
  ZALPHA=0.5**(-1.5)
  !
  DO JIJ=IIJB,IIJE
    IF (OOCEAN) THEN
      DO JK=IKTE,IKTB,-1
        ZD=ZALPHA*(PZZ(JIJ,IKTE+1)-PZZ(JIJ,JK))
        IF ( PLM(JIJ,JK)>ZD) THEN
          PLM(JIJ,JK)=ZD
        ELSE
          EXIT
        ENDIF
      END DO
    ELSE
      DO JK=IKTB,IKTE
        ZD=ZALPHA*(0.5*(PZZ(JIJ,JK)+PZZ(JIJ,JK+IKL))&
        -PZZ(JIJ,IKB)) *PDIRCOSZW(JIJ)
        IF ( PLM(JIJ,JK)>ZD) THEN
          PLM(JIJ,JK)=ZD
        ELSE
          EXIT
        ENDIF
      END DO
    ENDIF
  END DO
END IF
!
!$mnh_expand_array(JIJ=IIJB:IIJE)
PLM(:,IKA) = PLM(:,IKB)
PLM(:,IKU) = PLM(:,IKE)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
IF (LHOOK) CALL DR_HOOK('TURB:DELT',1,ZHOOK_HANDLE2)
END SUBROUTINE DELT
!
!     ####################
      SUBROUTINE DEAR(PLM)
!     ####################
!!
!!****  *DEAR* routine to compute mixing length for DEARdorff case
!
!!    AUTHOR
!!    ------
!!
!!     M Tomasini      *Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   01/05
!!      I.Sandu (Sept.2006) : Modification of the stability criterion
!!                            (theta_v -> theta_l)
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!*       0.1   Declarations of dummy arguments
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PLM
!
!-------------------------------------------------------------------------------
!
!   initialize the mixing length with the mesh grid
IF (LHOOK) CALL DR_HOOK('TURB:DEAR',0,ZHOOK_HANDLE2)
IF ( TURBN%CTURBDIM /= '1DIM' ) THEN
  CALL MXF_PHY(D,PDXX,ZWORK1)
  IF (.NOT. O2D) THEN
    CALL MYF_PHY(D,PDYY,ZWORK2)
  END IF
END IF
! 1D turbulence scheme
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
PLM(:,IKTB:IKTE) = PZZ(:,IKTB+IKL:IKTE+IKL) - PZZ(:,IKTB:IKTE)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
!$mnh_expand_array(JIJ=IIJB:IIJE)
PLM(:,IKU) = PLM(:,IKE)
PLM(:,IKA) = PZZ(:,IKB) - PZZ(:,IKA)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
IF ( TURBN%CTURBDIM /= '1DIM' ) THEN  ! 3D turbulence scheme
  IF ( O2D) THEN
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PLM(:,:) = SQRT( PLM(:,:)*ZWORK1(:,:) )
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ELSE
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PLM(:,:) = (PLM(:,:)*ZWORK1(:,:) &
                                 * ZWORK2(:,:) ) ** (1./3.)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
END IF
!   compute a mixing length limited by the stability
!
CALL ETHETA(D,CST,KRR,KRRI,PTHLT,PRT,ZLOCPEXNM,ZATHETA,PSRCT,OOCEAN,OCOMPUTE_SRC,ZETHETA)
CALL EMOIST(D,CST,KRR,KRRI,PTHLT,PRT,ZLOCPEXNM,ZAMOIST,PSRCT,OOCEAN,ZEMOIST)
!
IF (KRR>0) THEN
  DO JK = IKTB+1,IKTE-1
    DO JIJ=IIJB,IIJE
      ZDTHLDZ(JIJ,JK)= 0.5*((PTHLT(JIJ,JK+IKL)-PTHLT(JIJ,JK    ))/PDZZ(JIJ,JK+IKL)+ &
                              (PTHLT(JIJ,JK    )-PTHLT(JIJ,JK-IKL))/PDZZ(JIJ,JK    ))
      ZDRTDZ(JIJ,JK) = 0.5*((PRT(JIJ,JK+IKL,1)-PRT(JIJ,JK    ,1))/PDZZ(JIJ,JK+IKL)+ &
                              (PRT(JIJ,JK    ,1)-PRT(JIJ,JK-IKL,1))/PDZZ(JIJ,JK    ))
      IF (OOCEAN) THEN
        ZVAR=CST%XG*(CST%XALPHAOC*ZDTHLDZ(JIJ,JK)-CST%XBETAOC*ZDRTDZ(JIJ,JK))
      ELSE
        ZVAR=CST%XG/PTHVREF(JIJ,JK)*                                                  &
           (ZETHETA(JIJ,JK)*ZDTHLDZ(JIJ,JK)+ZEMOIST(JIJ,JK)*ZDRTDZ(JIJ,JK))
      END IF
      !
      IF (ZVAR>0.) THEN
        PLM(JIJ,JK)=MAX(CST%XMNH_EPSILON,MIN(PLM(JIJ,JK), &
                      0.76* SQRT(PTKET(JIJ,JK)/ZVAR)))
      END IF
    END DO
  END DO
ELSE! For dry atmos or unsalted ocean runs
  DO JK = IKTB+1,IKTE-1
    DO JIJ=IIJB,IIJE
      ZDTHLDZ(JIJ,JK)= 0.5*((PTHLT(JIJ,JK+IKL)-PTHLT(JIJ,JK    ))/PDZZ(JIJ,JK+IKL)+ &
                              (PTHLT(JIJ,JK    )-PTHLT(JIJ,JK-IKL))/PDZZ(JIJ,JK    ))
      IF (OOCEAN) THEN
        ZVAR= CST%XG*CST%XALPHAOC*ZDTHLDZ(JIJ,JK)
      ELSE
        ZVAR= CST%XG/PTHVREF(JIJ,JK)*ZETHETA(JIJ,JK)*ZDTHLDZ(JIJ,JK)
      END IF
!
      IF (ZVAR>0.) THEN
        PLM(JIJ,JK)=MAX(CST%XMNH_EPSILON,MIN(PLM(JIJ,JK), &
                      0.76* SQRT(PTKET(JIJ,JK)/ZVAR)))
      END IF
    END DO
  END DO
END IF
!  special case near the surface
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZDTHLDZ(:,IKB)=(PTHLT(:,IKB+IKL)-PTHLT(:,IKB))/PDZZ(:,IKB+IKL)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
! For dry simulations
IF (KRR>0) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZDRTDZ(:,IKB)=(PRT(:,IKB+IKL,1)-PRT(:,IKB,1))/PDZZ(:,IKB+IKL)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
ELSE
  ZDRTDZ(:,IKB)=0
ENDIF
!
IF (OOCEAN) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZWORK2D(:)=CST%XG*(CST%XALPHAOC*ZDTHLDZ(:,IKB)-CST%XBETAOC*ZDRTDZ(:,IKB))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZWORK2D(:)=CST%XG/PTHVREF(:,IKB)*                                           &
              (ZETHETA(:,IKB)*ZDTHLDZ(:,IKB)+ZEMOIST(:,IKB)*ZDRTDZ(:,IKB))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END IF
!$mnh_expand_where(JIJ=IIJB:IIJE)
WHERE(ZWORK2D(:)>0.)
  PLM(:,IKB)=MAX(CST%XMNH_EPSILON,MIN( PLM(:,IKB),                 &
                    0.76* SQRT(PTKET(:,IKB)/ZWORK2D(:))))
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE)
!
!  mixing length limited by the distance normal to the surface (with the same factor as for BL89)
!
IF (.NOT. TURBN%LRMC01) THEN
  ZALPHA=0.5**(-1.5)
  !
  DO JIJ=IIJB,IIJE
    IF (OOCEAN) THEN
      DO JK=IKTE,IKTB,-1
        ZD=ZALPHA*(PZZ(JIJ,IKTE+1)-PZZ(JIJ,JK))
        IF ( PLM(JIJ,JK)>ZD) THEN
          PLM(JIJ,JK)=ZD
        ELSE
          EXIT
        ENDIF
      END DO
    ELSE
      DO JK=IKTB,IKTE
        ZD=ZALPHA*(0.5*(PZZ(JIJ,JK)+PZZ(JIJ,JK+IKL))-PZZ(JIJ,IKB)) &
          *PDIRCOSZW(JIJ)
        IF ( PLM(JIJ,JK)>ZD) THEN
          PLM(JIJ,JK)=ZD
        ELSE
          EXIT
        ENDIF
      END DO
    ENDIF
  END DO
END IF
!
!$mnh_expand_array(JIJ=IIJB:IIJE)
PLM(:,IKA) = PLM(:,IKB)
PLM(:,IKE) = PLM(:,IKE-IKL)
PLM(:,IKU) = PLM(:,IKU-IKL)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
IF (LHOOK) CALL DR_HOOK('TURB:DEAR',1,ZHOOK_HANDLE2)
END SUBROUTINE DEAR
!
!     #########################
      SUBROUTINE CLOUD_MODIF_LM
!     #########################
!!
!!*****CLOUD_MODIF_LM routine to:
!!       1/ change the mixing length in the clouds
!!       2/ emphasize the mixing length in the cloud
!!           by the coefficient ZCOEF_AMPL calculated here
!!             when the CEI index is above ZCEI_MIN.
!!
!!
!!      ZCOEF_AMPL ^
!!                 |
!!                 |
!!  ZCOEF_AMPL_SAT -                       ---------- Saturation
!!    (XDUMMY1)    |                      -
!!                 |                     -
!!                 |                    -
!!                 |                   -
!!                 |                  - Amplification
!!                 |                 - straight
!!                 |                - line
!!                 |               -
!!                 |              -
!!                 |             -
!!                 |            -
!!                 |           -
!!               1 ------------
!!                 |
!!                 |
!!               0 -----------|------------|----------> PCEI
!!                 0      ZCEI_MIN     ZCEI_MAX
!!                        (XDUMMY2)    (XDUMMY3)
!!
!!
!!
!!    AUTHOR
!!    ------
!!     M. Tomasini   *CNRM METEO-FRANCE
!!
!!    MODIFICATIONS
!!    -------------
!!     Original   09/07/04
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATION
!              --------------
!
IF (LHOOK) CALL DR_HOOK('TURB:CLOUD_MODIF_LM',0,ZHOOK_HANDLE2)
ZPENTE = ( PCOEF_AMPL_SAT - 1. ) / ( PCEI_MAX - PCEI_MIN )
ZCOEF_AMPL_CEI_NUL = 1. - ZPENTE * PCEI_MIN
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZCOEF_AMPL(:,:) = 1.
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
!*       2.    CALCULATION OF THE AMPLIFICATION COEFFICIENT
!              --------------------------------------------
!
! Saturation
!
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE ( PCEI(:,:)>=PCEI_MAX ) 
  ZCOEF_AMPL(:,:)=PCOEF_AMPL_SAT
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
!
! Between the min and max limits of CEI index, linear variation of the
! amplification coefficient ZCOEF_AMPL as a function of CEI
!
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE ( PCEI(:,:) <  PCEI_MAX .AND. PCEI(:,:) >  PCEI_MIN)
  ZCOEF_AMPL(:,:) = ZPENTE * PCEI(:,:) + ZCOEF_AMPL_CEI_NUL
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
!
!
!*       3.    CALCULATION OF THE MIXING LENGTH IN CLOUDS
!              ------------------------------------------
!
IF (HTURBLEN_CL == TURBN%CTURBLEN) THEN
  ZLM_CLOUD(:,:) = ZLM(:,:)
ELSE
  SELECT CASE (HTURBLEN_CL)
!
!*         3.1 BL89 mixing length
!           ------------------
  CASE ('BL89','RM17','HM21')
    ZSHEAR(:,:)=0.
    CALL BL89(D,CST,CSTURB,TURBN,PZZ,PDZZ,PTHVREF,ZTHLM,KRR,ZRM,PTKET,ZSHEAR,ZLM_CLOUD,OOCEAN)
!
!*         3.2 Delta mixing length
!           -------------------
  CASE ('DELT')
    CALL DELT(ZLM_CLOUD,ODZ=.TRUE.)
!
!*         3.3 Deardorff mixing length
!           -----------------------
  CASE ('DEAR')
    CALL DEAR(ZLM_CLOUD)
!
  END SELECT
ENDIF
!
!*       4.    MODIFICATION OF THE MIXING LENGTH IN THE CLOUDS
!              -----------------------------------------------
!
! Impression before modification of the mixing length
IF ( TURBN%LTURB_DIAG .AND. TPFILE%LOPENED ) THEN
  TZFIELD = TFIELDMETADATA(            &
    CMNHNAME   = 'LM_CLEAR_SKY',       &
    CSTDNAME   = '',                   &
    CLONGNAME  = 'LM_CLEAR_SKY',       &
    CUNITS     = 'm',                  &
    CDIR       = 'XY',                 &
    CCOMMENT   = 'X_Y_Z_LM CLEAR SKY', &
    NGRID      = 1,                    &
    NTYPE      = TYPEREAL,             &
    NDIMS      = 3,                    &
    LTIMEDEP   = .TRUE.                )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZLM)
ENDIF
!
! Amplification of the mixing length when the criteria are verified
!
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE (ZCOEF_AMPL(:,:) /= 1.) 
  ZLM(:,:) = ZCOEF_AMPL(:,:)*ZLM_CLOUD(:,:)
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
!
! Cloud mixing length in the clouds at the points which do not verified the CEI
!
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE (PCEI(:,:) == -1.)
  ZLM(:,:) = ZLM_CLOUD(:,:)
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
!
!
!*       5.    IMPRESSION
!              ----------
!
IF ( TURBN%LTURB_DIAG .AND. TPFILE%LOPENED ) THEN
  TZFIELD = TFIELDMETADATA(         &
    CMNHNAME   = 'COEF_AMPL',       &
    CSTDNAME   = '',                &
    CLONGNAME  = 'COEF_AMPL',       &
    CUNITS     = '1',               &
    CDIR       = 'XY',              &
    CCOMMENT   = 'X_Y_Z_COEF AMPL', &
    NGRID      = 1,                 &
    NTYPE      = TYPEREAL,          &
    NDIMS      = 3,                 &
    LTIMEDEP   = .TRUE.             )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZCOEF_AMPL)
  !
  TZFIELD = TFIELDMETADATA(        &
    CMNHNAME   = 'LM_CLOUD',       &
    CSTDNAME   = '',               &
    CLONGNAME  = 'LM_CLOUD',       &
    CUNITS     = 'm',              &
    CDIR       = 'XY',             &
    CCOMMENT   = 'X_Y_Z_LM CLOUD', &
    NGRID      = 1,                &
    NTYPE      = TYPEREAL,         &
    NDIMS      = 3,                &
    LTIMEDEP   = .TRUE.            )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZLM_CLOUD)
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TURB:CLOUD_MODIF_LM',1,ZHOOK_HANDLE2)
END SUBROUTINE CLOUD_MODIF_LM
!
END SUBROUTINE TURB
