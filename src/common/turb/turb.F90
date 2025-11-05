!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
      SUBROUTINE TURB(CST,CSTURB,BUCONF,TURBN,NEBN,D,TLES,            &
              & KRR,KRRL,KRRI,HLBCX,HLBCY,KGRADIENTSLEO,              &
              & KGRADIENTSGOG,KHALO,                                  &
              & KSPLIT, OCLOUDMODIFLM, KSV,KSV_LGBEG,KSV_LGEND,       &
              & KSV_LIMA_NR, KSV_LIMA_NS, KSV_LIMA_NG, KSV_LIMA_NH,   &
              & O2D,ONOMIXLG,OFLAT,OCOUPLES,OBLOWSNOW,OIBM,OFLYER,    &
              & OCOMPUTE_SRC, PRSNOW,                                 &
              & OOCEAN,ODEEPOC,ODIAG_IN_RUN,                          &
              & HTURBLEN_CL,HCLOUD,HELEC,                             &
              & PTSTEP,TPFILE,                                        &
              & PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,                         &
              & PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,PCOSSLOPE,PSINSLOPE,    &
              & PRHODJ,PTHVREF,PHGRADLEO,PHGRADGOG,PZS,               &
              & PSFTH,PSFRV,PSFSV,PSFU,PSFV,                          &
              & PPABST,PUT,PVT,PWT,PTKET,PSVT,PSRCT,                  &
              & PLENGTHM,PLENGTHH,MFMOIST,                            &
              & PBL_DEPTH,PSBL_DEPTH,                                 &
              & PCEI,PCEI_MIN,PCEI_MAX,PCOEF_AMPL_SAT,                &
              & PTHLT,PRT,                                            &
              & PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,PRTKES,              &
              & PSIGS,                                                &
              & PFLXZTHVMF, PFLXZUMF, PFLXZVMF,                       &
              & PWTH,PWRC,PWSV,PDP,PTP,PTDIFF,PTDISS,      &
              & TBUDGETS, KBUDGETS,                                   &
              & PEDR,PLEM,PRTKEMS,PDPMF,PTPMF,                        &
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
!  A. Marcel Jan 2025: EDMF contribution to dynamic TKE production
! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODE_SHUMAN_PHY, ONLY: MZF_PHY,MXF_PHY,MYF_PHY
USE YOMHOOK ,   ONLY: LHOOK, DR_HOOK, JPHOOK
!
USE MODD_BUDGET,     ONLY:  NBUDGET_U,  NBUDGET_V,  NBUDGET_W,  NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, &
                            NBUDGET_RI, NBUDGET_SV1, NBUDGET_RG, NBUDGET_RH, NBUDGET_RR, NBUDGET_RS, &
                            TBUDGETDATA_PTR, TBUDGETCONF_T
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
USE MODE_CLOUD_MODIF_LM,      ONLY: CLOUD_MODIF_LM
USE MODE_COMPUTE_FUNCTION_THERMO_NEW_STAT, ONLY: COMPUTE_FUNCTION_THERMO_NEW_STAT
USE MODE_COMPUTE_FUNCTION_THERMO, ONLY: COMPUTE_FUNCTION_THERMO
USE MODE_DEAR,                ONLY: DEAR
USE MODE_DELT,                ONLY: DELT
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
USE MODE_MSG,                 ONLY: PRINT_MSG, NVERB_FATAL
!
USE MODI_LES_MEAN_SUBGRID_PHY
USE MODI_SECOND_MNH,          ONLY: SECOND_MNH
!
! These macro are handled by pft_tool.py --craybyPassDOCONCURRENT applied on Cray Rules
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif
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
INTEGER,                INTENT(IN)   :: KGRADIENTSLEO ! Number of stored horizontal gradients for Moeng scheme
INTEGER,                INTENT(IN)   :: KGRADIENTSGOG ! Number of stored horizontal gradients for Goger scheme
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
TYPE(TFILEDATA),        INTENT(INOUT)   ::  TPFILE       ! Output file
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
REAL, DIMENSION(D%NIJT,D%NKT,KGRADIENTSLEO),   INTENT(IN) ::  PHGRADLEO  ! horizontal gradients in Moeng
REAL, DIMENSION(D%NIJT,D%NKT,KGRADIENTSGOG),   INTENT(IN) ::  PHGRADGOG  ! horizontal gradients in Goger
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
REAL, DIMENSION(D%NIJT),INTENT(INOUT) :: PBL_DEPTH  ! BL height for TOMS
REAL, DIMENSION(D%NIJT),INTENT(INOUT) :: PSBL_DEPTH ! SBL depth for RMC01
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
REAL, DIMENSION(MERGE(D%NIJT,0,OFLYER),MERGE(D%NKT,0,OFLYER)), INTENT(OUT)  :: PWTH       ! heat flux
REAL, DIMENSION(MERGE(D%NIJT,0,OFLYER),MERGE(D%NKT,0,OFLYER)), INTENT(OUT)  :: PWRC       ! cloud water flux
REAL, DIMENSION(MERGE(D%NIJT,0,OFLYER),MERGE(D%NKT,0,OFLYER),MERGE(KSV,0,OFLYER)),INTENT(OUT) :: PWSV       ! scalar flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PFLXZUMF   ! MF contribution for vert. turb. transport (dyn prod)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PFLXZVMF   ! MF contribution for vert. turb. transport (dyn prod)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTP        ! Thermal TKE production MassFlux + turb
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL  :: PTPMF      ! Thermal TKE production MassFlux Only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PDP        ! Dynamic TKE production MassFlux + turb
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL  :: PDPMF      ! Dynamic TKE production MassFlux Only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTDIFF     ! Diffusion TKE term
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTDISS     ! Dissipation TKE term
!
TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
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
REAL, DIMENSION(MERGE(D%NIJT,0,ODIAG_IN_RUN),MERGE(D%NKT,0,ODIAG_IN_RUN)), INTENT(INOUT), OPTIONAL  ::  PCURRENT_TKE_DISS ! if ODIAG_IN_RUN in mesonh
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSTFL        ! Time evol Flux of T at sea surface (LOCEAN)
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSTFL_C  ! O-A interface flux for theta(LOCEAN and LCOUPLES)
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSRFL_C  ! O-A interface flux for vapor (LOCEAN and LCOUPLES) 
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSUFL_C        ! Time evol Flux of U at sea surface (LOCEAN)
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSVFL_C  !
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSUFL   
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSVFL  !
!
REAL, DIMENSION(MERGE(D%NIJT,0,OIBM),MERGE(D%NKT,0,OIBM)), INTENT(OUT), OPTIONAL :: PIBM_XMUT ! IBM turbulent viscosity
REAL, DIMENSION(MERGE(D%NIJT,0,OIBM),MERGE(D%NKT,0,OIBM)), INTENT(IN), OPTIONAL  :: PIBM_LS ! IBM Level-set function
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
          ZWORK1,ZWORK2,              &  ! working array syntax
          ZTEMP_BUD
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
LOGICAL :: GOCEAN !Intermediate variable used to work around a Cray compiler bug (CCE 13.0.0)
LOGICAL :: GTURBLEN_BL89_TURBLEN_RM17_TURBLEN_HM21_ORMC01
TYPE(TFIELDMETADATA) :: TZFIELD
!
REAL, DIMENSION(D%NIJT,D%NKT,MERGE(KSV+KRR,KSV,TURBN%LTURB_PRECIP)) :: ZWORKT, ZWORKS
REAL, DIMENSION(D%NIJT,      MERGE(KSV+KRR,KSV,TURBN%LTURB_PRECIP)) :: ZWORKSFSV
REAL, DIMENSION(D%NIJT,D%NKT,MERGE(KSV+KRR,KSV,TURBN%LTURB_PRECIP)) :: ZWORKWSV
INTEGER :: ISV
!
!*      1.PRELIMINARIES
!         -------------
!
!*      1.1 Set the internal domains, ZEXPL
!
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TURB',0,ZHOOK_HANDLE)
!
IF (TURBN%LHARAT .AND. TURBN%CTURBDIM /= '1DIM') THEN
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'TURB', 'TURBN%LHARATU only implemented for option TURBN%CTURBDIM=1DIM!')
ENDIF
IF (TURBN%LHARAT .AND. TLES%LLES_CALL) THEN
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'TURB', 'TURBN%LHARATU not implemented for option LLES_CALL')
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
GOCEAN = OOCEAN
GTURBLEN_BL89_TURBLEN_RM17_TURBLEN_HM21_ORMC01 = & 
TURBN%CTURBLEN=='BL89' .OR. TURBN%CTURBLEN=='RM17' .OR. TURBN%CTURBLEN=='HM21' .OR. TURBN%LRMC01
!Copy data into ZTHLM and ZRM only if needed
IF (GTURBLEN_BL89_TURBLEN_RM17_TURBLEN_HM21_ORMC01) THEN
!$acc kernels present_cr(ZTHLM,ZRM)
  ZTHLM(IIJB:IIJE,1:IKT) = PTHLT(IIJB:IIJE,1:IKT)
  ZRM(IIJB:IIJE,1:IKT,:) = PRT(IIJB:IIJE,1:IKT,:)
!$acc end kernels
END IF
!
!Save LIMA scalar variables sources
ZRSVS(IIJB:IIJE,1:IKT,1:KSV)=PRSVS(IIJB:IIJE,1:IKT,1:KSV)
!
ISV=KSV
IF (TURBN%LTURB_PRECIP) ISV=KSV+KRR
ZWORKT(:,:,1:KSV)=PSVT(:,:,:)
ZWORKS(:,:,1:KSV)=PRSVS(:,:,:)
IF (TURBN%LTURB_PRECIP) ZWORKT(:,:,KSV+1:KSV+KRR)=PRT(:,:,:)
IF (TURBN%LTURB_PRECIP) ZWORKS(:,:,KSV+1:KSV+KRR)=PRRS(:,:,:)
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
!$acc kernels present_cr(ZCOEF_DISS,ZTHLM,ZRM,zcp)
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZCP(IIJB:IIJE,1:IKT)=CST%XCPD
!
IF (KRR > 0) ZCP(IIJB:IIJE,1:IKT) = ZCP(IIJB:IIJE,1:IKT) + CST%XCPV * PRT(IIJB:IIJE,1:IKT,1)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
! PGI20.5 BUG or reproductibility problem , with pointer this loop on JRR parallelize whitout reduction 
!$acc loop seq
DO JRR = 2,1+KRRL                          ! loop on the liquid components
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)  
  ZCP(IIJB:IIJE,1:IKT)  = ZCP(IIJB:IIJE,1:IKT) + CST%XCL * PRT(IIJB:IIJE,1:IKT,JRR)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END DO
!
!$acc loop seq
DO JRR = 2+KRRL,1+KRRL+KRRI                ! loop on the solid components   
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZCP(IIJB:IIJE,1:IKT)  = ZCP(IIJB:IIJE,1:IKT)  + CST%XCI * PRT(IIJB:IIJE,1:IKT,JRR)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END DO
!
!*      2.2 Exner function at t
!
IF (GOCEAN) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZEXN(IIJB:IIJE,1:IKT) = 1.
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZEXN(IIJB:IIJE,1:IKT) = (PPABST(IIJB:IIJE,1:IKT)/CST%XP00) ** (CST%XRD/CST%XCPD)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
!*      2.3 dissipative heating coeff a t
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZCOEF_DISS(IIJB:IIJE,1:IKT) = 1/(ZCP(IIJB:IIJE,1:IKT) * ZEXN(IIJB:IIJE,1:IKT))
!
!
ZFRAC_ICE(IIJB:IIJE,1:IKT) = 0.0
ZATHETA(IIJB:IIJE,1:IKT) = 0.0
ZAMOIST(IIJB:IIJE,1:IKT) = 0.0
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
!
IF (KRRL >=1) THEN
  !
  !*      2.4 Temperature at t
  !
  !$acc kernels present_cr(ZT)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZT(IIJB:IIJE,1:IKT) =  PTHLT(IIJB:IIJE,1:IKT) * ZEXN(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !$acc end kernels
  !
  !*       2.5 Lv/Cph/Exn
  !
  IF ( KRRI >= 1 ) THEN
    IF (NEBN%LSTATNW) THEN
       !wc call new functions depending on statnew
       CALL COMPUTE_FUNCTION_THERMO_NEW_STAT(D, CST, CST%XALPW,CST%XBETAW,CST%XGAMW,CST%XLVTT,CST%XCL,ZT,ZEXN,ZCP, &
                                 ZLVOCPEXNM,ZAMOIST,ZATHETA, PPABST)
       CALL COMPUTE_FUNCTION_THERMO_NEW_STAT(D, CST, CST%XALPI,CST%XBETAI,CST%XGAMI,CST%XLSTT,CST%XCI,ZT,ZEXN,ZCP, &
                                 ZLSOCPEXNM,ZAMOIST_ICE,ZATHETA_ICE, PPABST)
    ELSE
      CALL COMPUTE_FUNCTION_THERMO(D,CST,CST%XALPW,CST%XBETAW,CST%XGAMW,CST%XLVTT,CST%XCL,ZT,ZEXN,ZCP, &
                                 ZLVOCPEXNM,ZAMOIST,ZATHETA,PRT,PPABST,KRR)
      CALL COMPUTE_FUNCTION_THERMO(D,CST,CST%XALPI,CST%XBETAI,CST%XGAMI,CST%XLSTT,CST%XCI,ZT,ZEXN,ZCP, &
                                 ZLSOCPEXNM,ZAMOIST_ICE,ZATHETA_ICE,PRT,PPABST,KRR)
    ENDIF
    !
!$acc kernels present_cr( zamoist, zatheta, zlocpexnm, zlvocpexnm, zlsocpexnm, zamoist_ice, zatheta_ice )
    !$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
    WHERE(PRT(IIJB:IIJE,1:IKT,2)+PRT(IIJB:IIJE,1:IKT,4)>0.0)
      ZFRAC_ICE(IIJB:IIJE,1:IKT) = PRT(IIJB:IIJE,1:IKT,4) / ( PRT(IIJB:IIJE,1:IKT,2) &
                                          +PRT(IIJB:IIJE,1:IKT,4) )
    END WHERE
    !$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZLOCPEXNM(IIJB:IIJE,1:IKT) = (1.0-ZFRAC_ICE(IIJB:IIJE,1:IKT))*ZLVOCPEXNM(IIJB:IIJE,1:IKT) &
                           +ZFRAC_ICE(IIJB:IIJE,1:IKT) *ZLSOCPEXNM(IIJB:IIJE,1:IKT)
    ZAMOIST(IIJB:IIJE,1:IKT) = (1.0-ZFRAC_ICE(IIJB:IIJE,1:IKT))*ZAMOIST(IIJB:IIJE,1:IKT) &
                         +ZFRAC_ICE(IIJB:IIJE,1:IKT) *ZAMOIST_ICE(IIJB:IIJE,1:IKT)
    ZATHETA(IIJB:IIJE,1:IKT) = (1.0-ZFRAC_ICE(IIJB:IIJE,1:IKT))*ZATHETA(IIJB:IIJE,1:IKT) &
                         +ZFRAC_ICE(IIJB:IIJE,1:IKT) *ZATHETA_ICE(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  ELSE
    !wc call new stat functions or not
    IF (NEBN%LSTATNW) THEN
      CALL COMPUTE_FUNCTION_THERMO_NEW_STAT(D,CST,CST%XALPW,CST%XBETAW,CST%XGAMW,CST%XLVTT,CST%XCL,ZT,ZEXN,ZCP, &
                                 ZLOCPEXNM,ZAMOIST,ZATHETA,PPABST)
    ELSE
      CALL COMPUTE_FUNCTION_THERMO(D,CST,CST%XALPW,CST%XBETAW,CST%XGAMW,CST%XLVTT,CST%XCL,ZT,ZEXN,ZCP, &
                                   ZLOCPEXNM,ZAMOIST,ZATHETA,PRT,PPABST,KRR)
    ENDIF
  END IF
  !
  !
  IF ( TPFILE%LOPENED .AND. TURBN%LTURB_DIAG ) THEN
!$acc update self(ZAMOIST,ZATHETA)
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
!$acc kernels present_cr( zlocpexnm )
  ZLOCPEXNM(IIJB:IIJE,1:IKT)=0.
!$acc end kernels
END IF              ! loop end on KRRL >= 1
!
! computes conservative variables
!
!$acc update device(PRRS,PRTHLS)
IF ( KRRL >= 1 ) THEN
!$acc kernels present_cr( zlocpexnm )
  IF ( KRRI >= 1 ) THEN
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ! Rnp at t
    PRT(IIJB:IIJE,1:IKT,1)  = PRT(IIJB:IIJE,1:IKT,1)  + PRT(IIJB:IIJE,1:IKT,2)  & 
                                    + PRT(IIJB:IIJE,1:IKT,4)
    PRRS(IIJB:IIJE,1:IKT,1) = PRRS(IIJB:IIJE,1:IKT,1) + PRRS(IIJB:IIJE,1:IKT,2) & 
                                    + PRRS(IIJB:IIJE,1:IKT,4)
    ! Theta_l at t
    PTHLT(IIJB:IIJE,1:IKT)  = PTHLT(IIJB:IIJE,1:IKT)  - ZLVOCPEXNM(IIJB:IIJE,1:IKT) &
                                    * PRT(IIJB:IIJE,1:IKT,2) &
                                  - ZLSOCPEXNM(IIJB:IIJE,1:IKT) * PRT(IIJB:IIJE,1:IKT,4)
    PRTHLS(IIJB:IIJE,1:IKT) = PRTHLS(IIJB:IIJE,1:IKT) - ZLVOCPEXNM(IIJB:IIJE,1:IKT) &
                                    * PRRS(IIJB:IIJE,1:IKT,2) &
                                  - ZLSOCPEXNM(IIJB:IIJE,1:IKT) * PRRS(IIJB:IIJE,1:IKT,4)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ELSE
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ! Rnp at t
    PRT(IIJB:IIJE,1:IKT,1)  = PRT(IIJB:IIJE,1:IKT,1)  + PRT(IIJB:IIJE,1:IKT,2)
    PRRS(IIJB:IIJE,1:IKT,1) = PRRS(IIJB:IIJE,1:IKT,1) + PRRS(IIJB:IIJE,1:IKT,2)
    ! Theta_l at t
    PTHLT(IIJB:IIJE,1:IKT)  = PTHLT(IIJB:IIJE,1:IKT)  - ZLOCPEXNM(IIJB:IIJE,1:IKT) &
                                    * PRT(IIJB:IIJE,1:IKT,2)
    PRTHLS(IIJB:IIJE,1:IKT) = PRTHLS(IIJB:IIJE,1:IKT) - ZLOCPEXNM(IIJB:IIJE,1:IKT) &
                                    * PRRS(IIJB:IIJE,1:IKT,2)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
!$acc end kernels
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
!$acc kernels present_cr(ZSHEAR)
    ZSHEAR(:,:)=0.
!$acc end kernels
    CALL BL89(D,CST,CSTURB,TURBN,PZZ,PDZZ,PTHVREF,ZTHLM,KRR,ZRM,PTKET,ZSHEAR,ZLM,GOCEAN)
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
    !$acc kernels
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZSHEAR(IIJB:IIJE,1:IKT) = SQRT(ZDUDZ(IIJB:IIJE,1:IKT)*ZDUDZ(IIJB:IIJE,1:IKT) &
                                    + ZDVDZ(IIJB:IIJE,1:IKT)*ZDVDZ(IIJB:IIJE,1:IKT))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL BL89(D,CST,CSTURB,TURBN,PZZ,PDZZ,PTHVREF,ZTHLM,KRR,ZRM,PTKET,ZSHEAR,ZLM,GOCEAN)
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
    !$acc kernels present_cr( ZSHEAR, ZDUDZ, ZDVDZ)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZSHEAR(IIJB:IIJE,1:IKT) = SQRT(ZDUDZ(IIJB:IIJE,1:IKT)*ZDUDZ(IIJB:IIJE,1:IKT) &
                                    + ZDVDZ(IIJB:IIJE,1:IKT)*ZDVDZ(IIJB:IIJE,1:IKT))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL BL89(D,CST,CSTURB,TURBN,PZZ,PDZZ,PTHVREF,ZTHLM,KRR,ZRM,PTKET,ZSHEAR,ZLM,GOCEAN)

    CALL DELT(D, TURBN, O2D, .FALSE., GOCEAN, PZZ, PDYY, PDXX, PDIRCOSZW, ZLMW)
    ! The minimum mixing length is chosen between Horizontal grid mesh (not taking into account the vertical grid mesh) and RM17.
    ! For large horizontal grid meshes, this is equal to RM17
    ! For LES grid meshes, this is equivalent to Deardorff : the base mixing lentgh is the horizontal grid mesh,
    !and it is limited by a stability-based length (RM17), as was done in Deardorff length (but taking into account shear as well)
    ! For grid meshes in the grey zone, then this is the smaller of the two.
    !
    !$acc kernels present_cr(ZLM, ZLMW)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZLM(IIJB:IIJE,1:IKT) = MIN(ZLM(IIJB:IIJE,1:IKT),TURBN%XCADAP*ZLMW(IIJB:IIJE,1:IKT))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
  !
  !*      3.4 Delta mixing length
  !           -------------------
  !
  CASE ('DELT')
    CALL DELT(D, TURBN, O2D, .TRUE., GOCEAN, PZZ, PDYY, PDXX, PDIRCOSZW, ZLM)
  !
  !*      3.5 Deardorff mixing length
  !           -----------------------
  !
  CASE ('DEAR')
    CALL DEAR(D, CST, TURBN, KRR, KRRI, O2D, OCOMPUTE_SRC, GOCEAN, &
    & ZLM, PRT, PDZZ, PZZ, PTKET, PTHVREF, &
    & PTHLT, ZLOCPEXNM, PSRCT, ZAMOIST, PDIRCOSZW, &
    & PDXX, PDYY, ZATHETA)
  !
  !*      3.6 Blackadar mixing length
  !           -----------------------
  !
  CASE ('BLKR')
   !$acc kernels
   ZL0 = 100.
   !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
   ZLM(IIJB:IIJE,1:IKT) = ZL0
   !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

   ZALPHA=0.5**(-1.5)
   !$acc end kernels
   !
   !$acc kernels
   DO JK=IKTB,IKTE
     !$mnh_expand_array(JIJ=IIJB:IIJE)
     ZLM(IIJB:IIJE,JK) = ( 0.5*(PZZ(IIJB:IIJE,JK)+PZZ(IIJB:IIJE,JK+IKL)) - &
     & PZZ(IIJB:IIJE,IKA+JPVEXT_TURB*IKL) ) * PDIRCOSZW(IIJB:IIJE)
     ZLM(IIJB:IIJE,JK) = ZALPHA  * ZLM(IIJB:IIJE,JK) * ZL0 / ( ZL0 + ZALPHA*ZLM(IIJB:IIJE,JK) )
     !$mnh_end_expand_array(JIJ=IIJB:IIJE)
   END DO
   !$acc end kernels
   !
   !$acc kernels
   !$mnh_expand_array(JIJ=IIJB:IIJE)
   ZLM(IIJB:IIJE,IKTB-1) = ZLM(IIJB:IIJE,IKTB)
   ZLM(IIJB:IIJE,IKTE+1) = ZLM(IIJB:IIJE,IKTE)
   !$mnh_end_expand_array(JIJ=IIJB:IIJE)
   !$acc end kernels
   !
   !
   !
END SELECT
!
!*      3.5 Mixing length modification for cloud
!           -----------------------
IF (OCLOUDMODIFLM) CALL CLOUD_MODIF_LM(D, CST, CSTURB, TURBN, TPFILE, TZFIELD, KRR, KRRI, &
  & OCLOUDMODIFLM, GOCEAN, OCOMPUTE_SRC, O2D, HTURBLEN_CL, &
  & PDZZ, PDXX, PDYY, PZZ, &
  & PRT, PTKET, PTHLT, ZTHLM, ZRM, PTHVREF, &
  & ZLOCPEXNM, PSRCT, PCOEF_AMPL_SAT, ZAMOIST, ZATHETA, PDIRCOSZW,  &
  & PCEI, PCEI_MIN, PCEI_MAX, ZLM)
ENDIF  ! end LHARRAT

!
!*      3.6 Dissipative length
!           ------------------

IF (TURBN%LHARAT) THEN
  !$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZLEPS(IIJB:IIJE,1:IKT)=PLENGTHM(IIJB:IIJE,1:IKT)*(3.75**2.)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !$acc end kernels
ELSE
  !$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZLEPS(IIJB:IIJE,1:IKT)=ZLM(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !$acc end kernels
ENDIF
!
!*      3.7 Correction in the Surface Boundary Layer (Redelsperger 2001)
!           ----------------------------------------
!
!$acc kernels
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZLMO(IIJB:IIJE)=XUNDEF
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
IF (TURBN%LRMC01) THEN
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZUSTAR(IIJB:IIJE)=(PSFU(IIJB:IIJE)**2+PSFV(IIJB:IIJE)**2)**(0.25)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
  IF (KRR>0) THEN
    CALL LMO(D,CST,ZUSTAR,ZTHLM(:,IKB),ZRM(:,IKB,1),PSFTH,PSFRV,ZLMO)
  ELSE
  !$acc kernels
    ZRVM(:)=0.
    ZSFRV(:)=0.
    !$acc end kernels
    CALL LMO(D,CST,ZUSTAR,ZTHLM(:,IKB),ZRVM,PSFTH,ZSFRV,ZLMO)
  END IF
  CALL RMC01(D,CST,CSTURB,TURBN,PZZ,PDXX,PDYY,PDZZ,PDIRCOSZW,PSBL_DEPTH,ZLMO,ZLM,ZLEPS)
END IF
!
!RMC01 is only applied on RM17 in HM21
IF (TURBN%CTURBLEN=='HM21') THEN
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZLEPS(IIJB:IIJE,1:IKT) = MIN(ZLEPS(IIJB:IIJE,1:IKT),ZLMW(IIJB:IIJE,1:IKT)*TURBN%XCADAP)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
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
!$acc kernels present_cr(ZUSLOPE,ZVSLOPE)
  ZUSLOPE(IIJB:IIJE)=PUT(IIJB:IIJE,IKA)
  ZVSLOPE(IIJB:IIJE)=PVT(IIJB:IIJE,IKA)
!$acc end kernels
END IF
IF (GOCEAN) THEN
!$acc kernels present_cr(ZUSLOPE,ZVSLOPE)
  ZUSLOPE(IIJB:IIJE)=PUT(IIJB:IIJE,IKU-1)
  ZVSLOPE(IIJB:IIJE)=PVT(IIJB:IIJE,IKU-1)
!$acc end kernels
END IF
!
!
!*      4.2 compute the proportionality coefficient between wind and stress
!
!$acc kernels
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZCDUEFF(IIJB:IIJE) =-SQRT ( (PSFU(IIJB:IIJE)**2 + PSFV(IIJB:IIJE)**2) /               &
                    (CST%XMNH_TINY + ZUSLOPE(IIJB:IIJE)**2 + ZVSLOPE(IIJB:IIJE)**2 ) )
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
!
!*       4.6 compute the surface tangential fluxes
!
!$acc kernels present_cr(ZTAU22M,ZTAU33M)
IF (GOCEAN) THEN
  ZTAU11M(IIJB:IIJE)=0.
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE)        
  ZTAU11M(IIJB:IIJE) =2./3.*(  (1.+ (PZZ(IIJB:IIJE,IKB+IKL)-PZZ(IIJB:IIJE,IKB))  &
                           /(PDZZ(IIJB:IIJE,IKB+IKL)+PDZZ(IIJB:IIJE,IKB))  &
                       )   *PTKET(IIJB:IIJE,IKB)                   &
                     -0.5  *PTKET(IIJB:IIJE,IKB+IKL)                 &
                    )
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END IF
ZTAU12M(IIJB:IIJE) =0.0
ZTAU22M(IIJB:IIJE) =ZTAU11M(IIJB:IIJE)
ZTAU33M(IIJB:IIJE) =ZTAU11M(IIJB:IIJE)
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
!$acc end kernels
!
IF (TURBN%CTOM=='TM06') THEN
  CALL TM06(D,CST,PTHVREF,PBL_DEPTH,PZZ,PSFTH,ZMWTH,ZMTH2)
  !
  CALL GZ_M_W_PHY(D,ZMWTH,PDZZ,ZWORK1)    ! -d(w'2th' )/dz
  CALL GZ_W_M_PHY(D,ZMTH2,PDZZ,ZWORK2)    ! -d(w'th'2 )/dz
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZFWTH(IIJB:IIJE,1:IKT) = -ZWORK1(IIJB:IIJE,1:IKT)
  ZFTH2(IIJB:IIJE,1:IKT) = -ZWORK2(IIJB:IIJE,1:IKT)
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
!$acc kernels present_cr(ZFWTH,ZFWR,ZFTH2,ZFR2,ZFTHR)
  ZFWTH(:,:) = 0.
  ZFWR(:,:)  = 0.
  ZFTH2(:,:) = 0.
  ZFR2(:,:)  = 0.
  ZFTHR(:,:) = 0.
!$acc end kernels
ENDIF
!
!----------------------------------------------------------------------------
!
!*      5. TURBULENT SOURCES
!          -----------------
!
IF( BUCONF%LBUDGET_U )  CALL TBUDGETS(NBUDGET_U )%PTR%INIT_PHY(D, 'VTURB', PRUS    )
IF( BUCONF%LBUDGET_V )  CALL TBUDGETS(NBUDGET_V )%PTR%INIT_PHY(D, 'VTURB', PRVS    )
IF( BUCONF%LBUDGET_W )  CALL TBUDGETS(NBUDGET_W )%PTR%INIT_PHY(D, 'VTURB', PRWS    )

IF( BUCONF%LBUDGET_TH ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRTHLS(:,:) + ZLVOCPEXNM(:,:) * PRRS(:,:, 2) + ZLSOCPEXNM(:,:) * PRRS(:,:, 4) 
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'VTURB', ZTEMP_BUD )
  ELSE IF( KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD, ZLOCPEXNM)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRTHLS(:,:) + ZLOCPEXNM(:,:) * PRRS(:,:, 2)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'VTURB', ZTEMP_BUD )
  ELSE
    CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'VTURB', PRTHLS )
  END IF
END IF

IF( BUCONF%LBUDGET_RV ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRRS(:,:, 1) - PRRS(:,:, 2) - PRRS(:,:, 4) 
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'VTURB', ZTEMP_BUD )
  ELSE IF( KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRRS(:,:, 1) - PRRS(:,:, 2)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'VTURB', ZTEMP_BUD )
  ELSE
    CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'VTURB', PRRS(:,:, 1) )
  END IF
END IF

IF( BUCONF%LBUDGET_RC ) CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'VTURB', PRRS  (:,:, 2) )
IF( BUCONF%LBUDGET_RR ) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'VTURB', PRRS  (:,:, 3) )
IF( BUCONF%LBUDGET_RI ) CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'VTURB', PRRS  (:,:, 4) )
IF( BUCONF%LBUDGET_RS ) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'VTURB', PRRS  (:,:, 5) )
IF( BUCONF%LBUDGET_RG ) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'VTURB', PRRS  (:,:, 6) )
IF( BUCONF%LBUDGET_RH .AND. KRR==7) CALL TBUDGETS(NBUDGET_RH)%PTR%INIT_PHY(D, 'VTURB', PRRS  (:,:, 7) )

IF( BUCONF%LBUDGET_SV ) THEN
  DO JSV = 1, KSV
    CALL TBUDGETS(NBUDGET_SV1 - 1 + JSV)%PTR%INIT_PHY(D, 'VTURB', ZWORKS(:,:, JSV) )
  END DO
END IF
!$acc update device(PRHODJ)
!$acc update_crm device(PRUS,PRVS,PRWS,PRSVS)

CALL TURB_VER(D,CST,CSTURB,TURBN,NEBN,TLES,              &
          KRR,KRRL,KRRI,KGRADIENTSLEO,                   &
          GOCEAN, ODEEPOC, OCOMPUTE_SRC,                 &
          ISV,KSV_LGBEG,KSV_LGEND,                       &
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
          PSBL_DEPTH,ZLMO,PHGRADLEO,PZS,                 &
          PRUS,PRVS,PRWS,PRTHLS,PRRS,ZWORKS,             &
          PDP,PTP,PSIGS,PWTH,PWRC,ZWORKWSV,              &
          PSSTFL, PSSTFL_C, PSSRFL_C,PSSUFL_C,PSSVFL_C,  &
          PSSUFL,PSSVFL                                  )
!$acc update_crm self(PWTH,PWRC,PWSV)
!IF (HCLOUD == 'LIMA') THEN
!   IF (KSV_LIMA_NR.GT.0) PRSVS(:,:,KSV_LIMA_NR) = ZRSVS(:,:,KSV_LIMA_NR) 
!   IF (KSV_LIMA_NS.GT.0) PRSVS(:,:,KSV_LIMA_NS) = ZRSVS(:,:,KSV_LIMA_NS)
!   IF (KSV_LIMA_NG.GT.0) PRSVS(:,:,KSV_LIMA_NG) = ZRSVS(:,:,KSV_LIMA_NG) 
!   IF (KSV_LIMA_NH.GT.0) PRSVS(:,:,KSV_LIMA_NH) = ZRSVS(:,:,KSV_LIMA_NH)
!END IF
IF (TURBN%LTURB_PRECIP) THEN
   IF (KRR.GE.3) PRRS(:,:,3)=ZWORKS(:,:,KSV+3)
   IF (KRR.GE.5) PRRS(:,:,5)=ZWORKS(:,:,KSV+5)
   IF (KRR.GE.6) PRRS(:,:,6)=ZWORKS(:,:,KSV+6)
   IF (KRR.GE.7) PRRS(:,:,7)=ZWORKS(:,:,KSV+7)
END IF

IF (TURBN%LTURB_PRECIP) THEN
  IF( BUCONF%LBUDGET_RR ) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'VTURB', PRRS(:,:, 3) )
  IF( BUCONF%LBUDGET_RS ) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'VTURB', PRRS(:,:, 5) )
  IF( BUCONF%LBUDGET_RG ) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'VTURB', PRRS(:,:, 6) )
  IF( BUCONF%LBUDGET_RH .AND. KRR ==7) CALL TBUDGETS(NBUDGET_RH)%PTR%INIT_PHY(D, 'VTURB', PRRS(:,:, 7) )
  IF (KRR.GE.3) PRRS(:,:,3)=ZWORKS(:,:,KSV+3)
  IF (KRR.GE.5) PRRS(:,:,5)=ZWORKS(:,:,KSV+5)
  IF (KRR.GE.6) PRRS(:,:,6)=ZWORKS(:,:,KSV+6)
  IF (KRR.GE.7) PRRS(:,:,7)=ZWORKS(:,:,KSV+7)
  IF( BUCONF%LBUDGET_RR ) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'VTURB', PRRS(:,:, 3) )
  IF( BUCONF%LBUDGET_RS ) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'VTURB', PRRS(:,:, 5) )
  IF( BUCONF%LBUDGET_RG ) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'VTURB', PRRS(:,:, 6) )
  IF( BUCONF%LBUDGET_RH .AND. KRR ==7) CALL TBUDGETS(NBUDGET_RH)%PTR%END_PHY(D, 'VTURB', PRRS(:,:, 7) )
END IF

IF( BUCONF%LBUDGET_U ) CALL TBUDGETS(NBUDGET_U)%PTR%END_PHY(D, 'VTURB', PRUS )
IF( BUCONF%LBUDGET_V ) CALL TBUDGETS(NBUDGET_V)%PTR%END_PHY(D, 'VTURB', PRVS )
IF( BUCONF%LBUDGET_W ) CALL TBUDGETS(NBUDGET_W)%PTR%END_PHY(D, 'VTURB', PRWS )

IF( BUCONF%LBUDGET_TH ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRTHLS(:,:) + ZLVOCPEXNM(:,:) * PRRS(:,:, 2) + ZLSOCPEXNM(:,:) * PRRS(:,:, 4) 
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'VTURB', ZTEMP_BUD )
  ELSE IF( KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD, ZLOCPEXNM)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRTHLS(:,:) + ZLOCPEXNM(:,:) * PRRS(:,:, 2)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'VTURB', ZTEMP_BUD )
  ELSE
    CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'VTURB', PRTHLS )
  END IF
END IF

IF( BUCONF%LBUDGET_RV ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRRS(:,:, 1) - PRRS(:,:, 2) - PRRS(:,:, 4) 
     !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'VTURB', ZTEMP_BUD )
  ELSE IF( KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRRS(:,:, 1) - PRRS(:,:, 2) 
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'VTURB', ZTEMP_BUD)
  ELSE
    CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'VTURB', PRRS(:,:, 1) )
  END IF
END IF

IF( BUCONF%LBUDGET_RC ) CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'VTURB', PRRS(:,:, 2) )
IF( BUCONF%LBUDGET_RR ) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'VTURB', PRRS(:,:, 3) )
IF( BUCONF%LBUDGET_RI ) CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'VTURB', PRRS(:,:, 4) )
IF( BUCONF%LBUDGET_RS ) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'VTURB', PRRS(:,:, 5) )
IF( BUCONF%LBUDGET_RG ) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'VTURB', PRRS(:,:, 6) )
IF( BUCONF%LBUDGET_RH .AND. KRR==7) CALL TBUDGETS(NBUDGET_RH)%PTR%END_PHY(D, 'VTURB', PRRS(:,:, 7) )

IF( BUCONF%LBUDGET_SV )  THEN
  DO JSV = 1, KSV
    CALL TBUDGETS(NBUDGET_SV1 - 1 + JSV)%PTR%END_PHY(D, 'VTURB', ZWORKS(:,:, JSV) )
  END DO
END IF
!
IF( TURBN%CTURBDIM == '3DIM' ) THEN
  IF( BUCONF%LBUDGET_U  ) CALL TBUDGETS(NBUDGET_U)%PTR%INIT_PHY(D, 'HTURB', PRUS )
  IF( BUCONF%LBUDGET_V  ) CALL TBUDGETS(NBUDGET_V)%PTR%INIT_PHY(D, 'HTURB', PRVS )
  IF( BUCONF%LBUDGET_W  ) CALL TBUDGETS(NBUDGET_W)%PTR%INIT_PHY(D, 'HTURB', PRWS )

  IF(BUCONF%LBUDGET_TH)  THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRTHLS(:,:) + ZLVOCPEXNM(:,:) * PRRS(:,:, 2) + ZLSOCPEXNM(:,:) * PRRS(:,:, 4)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
      CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'HTURB', ZTEMP_BUD )
    ELSE IF( KRRL >= 1 ) THEN
      !$acc kernels present_cr(ZTEMP_BUD)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZTEMP_BUD(:,:) =  PRTHLS(:,:) + ZLOCPEXNM(:,:) * PRRS(:,:, 2)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      !$acc end kernels
      CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'HTURB', ZTEMP_BUD  )
    ELSE
      CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'HTURB', PRTHLS )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RV ) THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
      !$acc kernels present_cr(ZTEMP_BUD)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZTEMP_BUD(:,:) =  PRRS(:,:, 1) - PRRS(:,:, 2) - PRRS(:,:, 4)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      !$acc end kernels
      CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'HTURB', ZTEMP_BUD )
    ELSE IF( KRRL >= 1 ) THEN
      !$acc kernels present_cr(ZTEMP_BUD)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZTEMP_BUD(:,:) =  PRRS(:,:, 1) - PRRS(:,:, 2)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      !$acc end kernels
      CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'HTURB', ZTEMP_BUD )
    ELSE
      CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'HTURB', PRRS(:,:, 1) )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RC ) CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'HTURB', PRRS(:,:, 2) )
  IF( BUCONF%LBUDGET_RR ) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'HTURB', PRRS(:,:, 3) )
  IF( BUCONF%LBUDGET_RI ) CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'HTURB', PRRS(:,:, 4) )
  IF( BUCONF%LBUDGET_RS ) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'HTURB', PRRS(:,:, 5) )
  IF( BUCONF%LBUDGET_RG ) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'HTURB', PRRS(:,:, 6) )
  IF( BUCONF%LBUDGET_RH .AND. KRR==7) CALL TBUDGETS(NBUDGET_RH)%PTR%INIT_PHY(D, 'HTURB', PRRS(:,:, 7) )

  IF( BUCONF%LBUDGET_SV )  THEN
    DO JSV = 1, KSV
      CALL TBUDGETS(NBUDGET_SV1 - 1 + JSV)%PTR%INIT_PHY(D, 'HTURB', ZWORKS(:,:, JSV) )
    END DO
  END IF
    CALL TURB_HOR_SPLT(D,CST,CSTURB, TURBN, NEBN, TLES,        &
          KSPLIT, KRR, KRRL, KRRI, ISV,KSV_LGBEG,KSV_LGEND,    & 
          PTSTEP,HLBCX,HLBCY, OFLAT,O2D, ONOMIXLG,             & 
          GOCEAN,OCOMPUTE_SRC,OBLOWSNOW,PRSNOW,                &
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
  IF (TURBN%LTURB_PRECIP) THEN
    IF( BUCONF%LBUDGET_RR ) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'HTURB', PRRS(:,:, 3) )
    IF( BUCONF%LBUDGET_RS ) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'HTURB', PRRS(:,:, 5) )
    IF( BUCONF%LBUDGET_RG ) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'HTURB', PRRS(:,:, 6) )
    IF( BUCONF%LBUDGET_RH .AND. KRR==7) CALL TBUDGETS(NBUDGET_RH)%PTR%INIT_PHY(D, 'HTURB', PRRS(:,:, 7) )
    IF (KRR.GE.3) PRRS(:,:,3)=ZWORKS(:,:,KSV+3)
    IF (KRR.GE.5) PRRS(:,:,5)=ZWORKS(:,:,KSV+5)
    IF (KRR.GE.6) PRRS(:,:,6)=ZWORKS(:,:,KSV+6)
    IF (KRR.GE.7) PRRS(:,:,7)=ZWORKS(:,:,KSV+7)
    IF( BUCONF%LBUDGET_RR ) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'HTURB', PRRS(:,:, 3) )
    IF( BUCONF%LBUDGET_RS ) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'HTURB', PRRS(:,:, 5) )
    IF( BUCONF%LBUDGET_RG ) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'HTURB', PRRS(:,:, 6) )
    IF( BUCONF%LBUDGET_RH .AND. KRR==7) CALL TBUDGETS(NBUDGET_RH)%PTR%END_PHY(D, 'HTURB', PRRS(:,:, 7) )
  END IF
  !
  IF( BUCONF%LBUDGET_U ) CALL TBUDGETS(NBUDGET_U)%PTR%END_PHY(D, 'HTURB', PRUS )
  IF( BUCONF%LBUDGET_V ) CALL TBUDGETS(NBUDGET_V)%PTR%END_PHY(D, 'HTURB', PRVS )
  IF( BUCONF%LBUDGET_W ) CALL TBUDGETS(NBUDGET_W)%PTR%END_PHY(D, 'HTURB', PRWS )

  IF( BUCONF%LBUDGET_TH ) THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
      !$acc kernels present_cr(ZTEMP_BUD)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZTEMP_BUD(:,:) =  PRTHLS(:,:) + ZLVOCPEXNM(:,:) * PRRS(:,:, 2) + ZLSOCPEXNM(:,:) * PRRS(:,:, 4)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      !$acc end kernels
      CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'HTURB', ZTEMP_BUD )
    ELSE IF( KRRL >= 1 ) THEN
      !$acc kernels present_cr(ZTEMP_BUD)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZTEMP_BUD(:,:) =  PRTHLS(:,:) + ZLOCPEXNM(:,:) * PRRS(:,:, 2)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      !$acc end kernels
      CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'HTURB', ZTEMP_BUD )
    ELSE
      CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'HTURB', PRTHLS )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RV ) THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
      !$acc kernels present_cr(ZTEMP_BUD)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZTEMP_BUD(:,:) =  PRRS(:,:, 1) - PRRS(:,:, 2) - PRRS(:,:, 4)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      !$acc end kernels
      CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'HTURB', ZTEMP_BUD )
    ELSE IF( KRRL >= 1 ) THEN
      !$acc kernels present_cr(ZTEMP_BUD)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZTEMP_BUD(:,:) =  PRRS(:,:, 1) - PRRS(:,:, 2)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      !$acc end kernels
      CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'HTURB', ZTEMP_BUD )
    ELSE
      CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'HTURB', PRRS(:,:, 1) )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RC ) CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'HTURB', PRRS(:,:, 2) )
  IF( BUCONF%LBUDGET_RR ) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'HTURB', PRRS(:,:, 3) )
  IF( BUCONF%LBUDGET_RI ) CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'HTURB', PRRS(:,:, 4) )
  IF( BUCONF%LBUDGET_RS ) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'HTURB', PRRS(:,:, 5) )
  IF( BUCONF%LBUDGET_RG ) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'HTURB', PRRS(:,:, 6) )
  IF( BUCONF%LBUDGET_RH .AND. KRR==7) CALL TBUDGETS(NBUDGET_RH)%PTR%END_PHY(D, 'HTURB', PRRS(:,:, 7) )

  IF( BUCONF%LBUDGET_SV )  THEN
    DO JSV = 1, KSV
      CALL TBUDGETS(NBUDGET_SV1 - 1 + JSV)%PTR%END_PHY(D, 'HTURB', ZWORKS(:,:, JSV) )
    END DO
  END IF
END IF
!$acc update_crm self(PSIGS,PRUS,PRVS,PRWS,PRSVS)
!----------------------------------------------------------------------------
!
!*      6. EVOLUTION OF THE TKE AND ITS DISSIPATION
!          ----------------------------------------
!
!  6.1 Contribution of mass-flux in the TKE buoyancy and dynamical production if
!      cloud computation is not statistical
IF(TURBN%LTHERMMF) THEN
  CALL MZF_PHY(D,PFLXZTHVMF,ZWORK1)
!$acc kernels present_crm(PTP)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PTP(IIJB:IIJE,1:IKT) = PTP(IIJB:IIJE,1:IKT) &
                               + CST%XG / PTHVREF(IIJB:IIJE,1:IKT) * ZWORK1(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  IF(PRESENT(PTPMF)) THEN
!$acc kernels present_crm(PTPMF)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PTPMF(IIJB:IIJE,1:IKT)=CST%XG / PTHVREF(IIJB:IIJE,1:IKT) * ZWORK1(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  ENDIF
ELSEIF(PRESENT(PTPMF))  THEN
!$acc kernels
  PTPMF(IIJB:IIJE,1:IKT)=0.
!$acc end kernels
ENDIF

IF(TURBN%LDYNMF) THEN
  CALL GZ_U_UW_PHY(D,PUT,PDZZ,ZWORK1)
  CALL MXF_PHY(D,ZWORK1,ZWORK2)
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK1(IIJB:IIJE,1:IKT)=ZWORK2(IIJB:IIJE,1:IKT)*PFLXZUMF(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  CALL MZF_PHY(D,ZWORK1,ZWORK2)
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PDP(IIJB:IIJE,1:IKT)=PDP(IIJB:IIJE,1:IKT)-ZWORK2(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  IF(PRESENT(PDPMF)) THEN
!$acc kernels
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PDPMF(IIJB:IIJE,1:IKT)=-ZWORK2(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  ENDIF
!
  CALL GZ_V_VW_PHY(D,PVT,PDZZ,ZWORK1)
  CALL MXF_PHY(D,ZWORK1,ZWORK2)
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK1(IIJB:IIJE,1:IKT)=ZWORK2(IIJB:IIJE,1:IKT)*PFLXZVMF(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  CALL MZF_PHY(D,ZWORK1,ZWORK2)
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PDP(IIJB:IIJE,1:IKT)=PDP(IIJB:IIJE,1:IKT)-ZWORK2(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  IF(PRESENT(PDPMF)) THEN
!$acc kernels
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PDPMF(IIJB:IIJE,1:IKT)=PDPMF(IIJB:IIJE,1:IKT)-ZWORK2(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  ENDIF
ELSEIF(PRESENT(PDPMF)) THEN
!$acc kernels
  PDPMF(IIJB:IIJE,1:IKT)=0.
!$acc end kernels
ENDIF

!  6.2 Horizontal gradients as in Göger et al. (2016)

IF (TURBN%LGOGER) THEN
  ! Add horizontal terms from Göger  et al. (2018)
  ! Increase the Dyn. Prod.
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !* Computation of the horizontal mixing length
  !* Add horizontal terms
  ! DUDX=PHGRADGOG(IIJB:IIJE,1:IKT,1)
  ! DUDY=PHGRADGOG(IIJB:IIJE,1:IKT,2)
  ! DVDX=PHGRADGOG(IIJB:IIJE,1:IKT,3)
  ! DVDY=PHGRADGOG(IIJB:IIJE,1:IKT,4)
  PDP(IIJB:IIJE,1:IKT)=PDP(IIJB:IIJE,1:IKT)+TURBN%XSMAG**2*PDXX(IIJB:IIJE,1:IKT)*PDYY(IIJB:IIJE,1:IKT)* &
                   &(PHGRADGOG(IIJB:IIJE,1:IKT,1)*PHGRADGOG(IIJB:IIJE,1:IKT,1)           &
                   &+PHGRADGOG(IIJB:IIJE,1:IKT,4)*PHGRADGOG(IIJB:IIJE,1:IKT,4)           &
                   &+0.5*(PHGRADGOG(IIJB:IIJE,1:IKT,2)+PHGRADGOG(IIJB:IIJE,1:IKT,3))     &
                   &*(PHGRADGOG(IIJB:IIJE,1:IKT,2)+PHGRADGOG(IIJB:IIJE,1:IKT,3)))**(3./2.)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
ENDIF

!  6.3 TKE evolution equation

IF (.NOT. TURBN%LHARAT) THEN
!
IF (BUCONF%LBUDGET_TH)  THEN
  IF ( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRTHLS(:,:)+ ZLVOCPEXNM(:,:) * PRRS(:,:,2) + ZLSOCPEXNM(:,:) * PRRS(:,:,4)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'DISSH', ZTEMP_BUD )
  ELSE IF ( KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRTHLS(:,:) + ZLOCPEXNM(:,:) * PRRS(:,:,2)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'DISSH', ZTEMP_BUD )
  ELSE
    CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'DISSH', PRTHLS )
  END IF
END IF
!
IF(PRESENT(PRTKEMS)) THEN
!$acc kernels
  ZRTKEMS(:,:)=PRTKEMS(:,:)
!$acc end kernels
ELSE
!$acc kernels
  ZRTKEMS(:,:)=0.
!$acc end kernels
END IF
!
!$acc update_crm device(PRTKES)
CALL TKE_EPS_SOURCES(D,CST,CSTURB,BUCONF,TURBN,TLES,                    &
                   & PTKET,ZLM,ZLEPS,PDP,ZTRH,                          &
                   & PRHODJ,PDZZ,PDXX,PDYY,PDZX,PDZY,PZZ,               &
                   & PTSTEP,ZEXPL,                                      &
                   & TPFILE,ODIAG_IN_RUN,GOCEAN,                        &
                   & PSFU,PSFV,                                         &
                   & PTP,PRTKES,PRTHLS,ZCOEF_DISS,PTDIFF,PTDISS,ZRTKEMS,&
                   & TBUDGETS,KBUDGETS, PEDR=PEDR, PTR=PTR,PDISS=PDISS, &
                   & PCURRENT_TKE_DISS=PCURRENT_TKE_DISS                )
                   !
!$acc update_crm self(PTR,PDISS)
!
!$acc update_crm self(PRTKES)
IF (BUCONF%LBUDGET_TH)  THEN
  IF ( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRTHLS(:,:)+ ZLVOCPEXNM(:,:) * PRRS(:,:,2) + ZLSOCPEXNM(:,:) * PRRS(:,:,4)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'DISSH', ZTEMP_BUD )
  ELSE IF ( KRRL >= 1 ) THEN
    !$acc kernels present_cr(ZTEMP_BUD)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZTEMP_BUD(:,:) =  PRTHLS(:,:) + ZLOCPEXNM(:,:) * PRRS(:,:,2)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$acc end kernels
    CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'DISSH', ZTEMP_BUD )
  ELSE
    CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'DISSH', PRTHLS )
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
!$acc update_crm self(PLEM)
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
!$acc update self(PTHLT)
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
!$acc update self(PRT)
    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PRT(:,:,1))
   END IF
END IF
!
PRSVS(:,:,:)        = ZWORKS(:,:,1:KSV)
IF (OFLYER)   PWSV(:,:,:)=ZWORKWSV(:,:,1:KSV)
!* stores value of conservative variables & wind before turbulence tendency (AROME only)
IF(PRESENT(PDRUS_TURB)) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PDRUS_TURB(IIJB:IIJE,1:IKT)   = PRUS(IIJB:IIJE,1:IKT) - PDRUS_TURB(IIJB:IIJE,1:IKT)
  PDRVS_TURB(IIJB:IIJE,1:IKT)   = PRVS(IIJB:IIJE,1:IKT) - PDRVS_TURB(IIJB:IIJE,1:IKT)
  PDRTHLS_TURB(IIJB:IIJE,1:IKT) = PRTHLS(IIJB:IIJE,1:IKT) - PDRTHLS_TURB(IIJB:IIJE,1:IKT)
  PDRRTS_TURB(IIJB:IIJE,1:IKT)  = PRRS(IIJB:IIJE,1:IKT,1) - PDRRTS_TURB(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT,JSV=1:KSV)  
  PDRSVS_TURB(IIJB:IIJE,1:IKT,:)  = PRSVS(IIJB:IIJE,1:IKT,:) - PDRSVS_TURB(IIJB:IIJE,1:IKT,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT,JSV=1:KSV)
END IF
!----------------------------------------------------------------------------
!
!*      8. RETRIEVE NON-CONSERVATIVE VARIABLES
!          -----------------------------------
!
IF ( KRRL >= 1 ) THEN
  IF ( KRRI >= 1 ) THEN
!$acc kernels present_cr(PRT,PRRS,PTHLT,PRTHLS)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PRT(IIJB:IIJE,1:IKT,1)  = PRT(IIJB:IIJE,1:IKT,1)  - PRT(IIJB:IIJE,1:IKT,2)  &
                                    - PRT(IIJB:IIJE,1:IKT,4)
    PRRS(IIJB:IIJE,1:IKT,1) = PRRS(IIJB:IIJE,1:IKT,1) - PRRS(IIJB:IIJE,1:IKT,2) &
                                    - PRRS(IIJB:IIJE,1:IKT,4)
    PTHLT(IIJB:IIJE,1:IKT)  = PTHLT(IIJB:IIJE,1:IKT)  + ZLVOCPEXNM(IIJB:IIJE,1:IKT) &
                                    * PRT(IIJB:IIJE,1:IKT,2) &
                                    + ZLSOCPEXNM(IIJB:IIJE,1:IKT) * PRT(IIJB:IIJE,1:IKT,4)
    PRTHLS(IIJB:IIJE,1:IKT) = PRTHLS(IIJB:IIJE,1:IKT) + ZLVOCPEXNM(IIJB:IIJE,1:IKT) &
                                    * PRRS(IIJB:IIJE,1:IKT,2) &
                                    + ZLSOCPEXNM(IIJB:IIJE,1:IKT) * PRRS(IIJB:IIJE,1:IKT,4)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
!$acc update self(PRT(:,:,1))
    !
  ELSE
!$acc kernels present_cr(PRT,PRRS,PTHLT,PRTHLS, zlocpexnm )
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PRT(IIJB:IIJE,1:IKT,1)  = PRT(IIJB:IIJE,1:IKT,1)  - PRT(IIJB:IIJE,1:IKT,2)
    PRRS(IIJB:IIJE,1:IKT,1) = PRRS(IIJB:IIJE,1:IKT,1) - PRRS(IIJB:IIJE,1:IKT,2)
    PTHLT(IIJB:IIJE,1:IKT)  = PTHLT(IIJB:IIJE,1:IKT)  + ZLOCPEXNM(IIJB:IIJE,1:IKT) &
                                    * PRT(IIJB:IIJE,1:IKT,2)
    PRTHLS(IIJB:IIJE,1:IKT) = PRTHLS(IIJB:IIJE,1:IKT) + ZLOCPEXNM(IIJB:IIJE,1:IKT) &
                                    * PRRS(IIJB:IIJE,1:IKT,2)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
!$acc update self(PRT(:,:,1))
  END IF
END IF!
!
!
! Remove non-physical negative values (unnecessary in a perfect world) + corresponding budgets
CALL SOURCES_NEG_CORRECT_PHY(D,KSV,HCLOUD,HELEC,'NETUR',KRR,PTSTEP,PPABST,PTHLT,PRT,PRTHLS,PRRS,PRSVS)
!$acc update self( PTHLT ) !PTHLT not modified in Sources_neg_correct
!$acc update self( PRTHLS, PRRS )
!----------------------------------------------------------------------------
!
!*      9. LES averaged surface fluxes
!          ---------------------------
!
IF (TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  !$acc data copy(TLES%X_LES_Q0,TLES%X_LES_E0,TLES%X_LES_SV0,TLES%X_LES_UW0,TLES%X_LES_VW0,TLES%X_LES_USTAR)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,PSFTH,TLES%X_LES_Q0)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,PSFRV,TLES%X_LES_E0)
  DO JSV=1,KSV
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,PSFSV(:,JSV),TLES%X_LES_SV0(:,JSV))
  END DO
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,PSFU,TLES%X_LES_UW0)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,PSFV,TLES%X_LES_VW0)
  !
  !$acc kernels present_cr(ZWORK2D)
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZWORK2D(IIJB:IIJE) = (PSFU(IIJB:IIJE)*PSFU(IIJB:IIJE)+PSFV(IIJB:IIJE)*PSFV(IIJB:IIJE))**0.25
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  !$acc end kernels
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2D,TLES%X_LES_USTAR)
!$acc end data
  !----------------------------------------------------------------------------
  !
  !*     10. LES for 3rd order moments
  !          -------------------------
  !
!$acc data copy(TLES%X_LES_SUBGRID_W2Thl,TLES%X_LES_SUBGRID_WThl2)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZMWTH,TLES%X_LES_SUBGRID_W2Thl)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZMTH2,TLES%X_LES_SUBGRID_WThl2)
!$acc end data
  IF (KRR>0) THEN
!$acc data copy(TLES%X_LES_SUBGRID_W2Rt,TLES%X_LES_SUBGRID_WThlRt,TLES%X_LES_SUBGRID_WRt2)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZMWR,TLES%X_LES_SUBGRID_W2Rt)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZMTHR,TLES%X_LES_SUBGRID_WThlRt)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZMR2,TLES%X_LES_SUBGRID_WRt2)
!$acc end data
  END IF
  !
  !----------------------------------------------------------------------------
  !
  !*     11. LES quantities depending on <w'2> in "1DIM" mode
  !          ------------------------------------------------
  !
  IF (TURBN%CTURBDIM=="1DIM") THEN
!$acc data copy(TLES%X_LES_SUBGRID_U2,TLES%X_LES_SUBGRID_V2,TLES%X_LES_SUBGRID_W2,TLES%X_LES_RES_ddz_Thl_SBG_W2)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK1(IIJB:IIJE,1:IKT) = 2./3.*PTKET(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK1,TLES%X_LES_SUBGRID_U2)
    TLES%X_LES_SUBGRID_V2(:,:,:) = TLES%X_LES_SUBGRID_U2(:,:,:)
    TLES%X_LES_SUBGRID_W2(:,:,:) = TLES%X_LES_SUBGRID_U2(:,:,:)
    !
    CALL GZ_M_W_PHY(D,PTHLT,PDZZ,ZWORK1)
    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK2(IIJB:IIJE,1:IKT)  = 2./3.*PTKET(IIJB:IIJE,1:IKT) *ZWORK2(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2,TLES%X_LES_RES_ddz_Thl_SBG_W2)
!$acc end data
    !
    IF (KRR>=1) THEN
!$acc data copy(TLES%X_LES_RES_ddz_Rt_SBG_W2)
      CALL GZ_M_W_PHY(D,PRT(:,:,1),PDZZ,ZWORK1)
      CALL MZF_PHY(D,ZWORK1,ZWORK2)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZWORK2(IIJB:IIJE,1:IKT)  = 2./3.*PTKET(IIJB:IIJE,1:IKT) *ZWORK2(IIJB:IIJE,1:IKT)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2,TLES%X_LES_RES_ddz_Rt_SBG_W2)
!$acc end data
    END IF
!$acc data copy(TLES%X_LES_RES_ddz_Sv_SBG_W2(:,:,:,1:KSV))    
    DO JSV=1,KSV
      CALL GZ_M_W_PHY(D,PSVT(:,:,JSV),PDZZ,ZWORK1)
      CALL MZF_PHY(D,ZWORK1,ZWORK2)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZWORK2(IIJB:IIJE,1:IKT)  = 2./3.*PTKET(IIJB:IIJE,1:IKT) *ZWORK2(IIJB:IIJE,1:IKT)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2, TLES%X_LES_RES_ddz_Sv_SBG_W2(:,:,:,JSV))
    END DO
!$acc end data
  END IF
  !
  !----------------------------------------------------------------------------
  !
  !*     12. LES mixing end dissipative lengths, presso-correlations
  !          -------------------------------------------------------
  !
!$acc data copy(TLES%X_LES_SUBGRID_LMix,TLES%X_LES_SUBGRID_LDiss,TLES%X_LES_SUBGRID_WP)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZLM,TLES%X_LES_SUBGRID_LMix)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZLEPS,TLES%X_LES_SUBGRID_LDiss)
  !
  !* presso-correlations for subgrid Tke are equal to zero.
  !
!$acc kernels present_cr(ZLEPS)
  ZLEPS(:,:) = 0. !ZLEPS is used as a work array (not used anymore)
!$acc end kernels
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZLEPS,TLES%X_LES_SUBGRID_WP)
!$acc end data
  !
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
IF(PRESENT(PLEM)) PLEM(IIJB:IIJE,IKTB:IKTE) = ZLM(IIJB:IIJE,IKTB:IKTE)
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TURB',1,ZHOOK_HANDLE)
END SUBROUTINE TURB
