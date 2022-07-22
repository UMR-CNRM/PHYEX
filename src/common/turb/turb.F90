!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
      SUBROUTINE TURB(CST,CSTURB,BUCONF,TURBN,D,              &
              & KMI,KRR,KRRL,KRRI,HLBCX,HLBCY,            &
              & KSPLIT,KMODEL_CL,KSV,KSV_LGBEG,KSV_LGEND,HPROGRAM,    &
              & O2D,ONOMIXLG,OFLAT,OLES_CALL,OCOUPLES,OBLOWSNOW,      &
              & OTURB_FLX,OTURB_DIAG,OSUBG_COND,OCOMPUTE_SRC,         &
              & ORMC01,OOCEAN,ODEEPOC,OHARAT, ODIAG_IN_RUN,           &
              & HTURBDIM,HTURBLEN,HTOM,HTURBLEN_CL,HCLOUD,PIMPL,      &
              & PTSTEP,TPFILE,PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,           &
              & PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,PCOSSLOPE,PSINSLOPE,    &
              & PRHODJ,PTHVREF,                                       &
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
              & PCURRENT_TKE_DISS                                     ) 
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
!!    means of the parameter HTURBDIM:
!!           * HTURBDIM='1DIM' the parameterization is 1D but can be used in
!!    3D , 2D or 1D simulations. Only the sources associated to the vertical
!!    turbulent fluxes are taken into account.
!!           *  HTURBDIM='3DIM' the parameterization is fully 2D or 3D depending
!!    on the model  dimensionality. Of course, it does not make any sense to
!!    activate this option with a 1D model. 
!!
!!      The following steps are made:
!!      1- Preliminary computations.
!!      2- The metric coefficients are recovered from the grid knowledge.
!!      3- The mixing length is computed according to its choice:
!!           * HTURBLEN='BL89' the Bougeault and Lacarrere algorithm is used.
!!             The mixing length is given by the vertical displacement from its
!!             original level of an air particule having an initial internal
!!             energy equal to its TKE and stopped by the buoyancy forces.
!!             The discrete formulation is second order accurate.
!!           * HTURBLEN='DELT' the mixing length is given by the mesh size 
!!             depending on the model dimensionality, this length is limited 
!!             with the ground distance.
!!           * HTURBLEN='DEAR' the mixing length is given by the mesh size 
!!             depending on the model dimensionality, this length is limited 
!!             with the ground distance and also by the Deardorff mixing length
!!             pertinent in the stable cases.
!!           * HTURBLEN='KEPS' the mixing length is deduced from the TKE 
!!             dissipation, which becomes a prognostic variable of the model (
!!             Duynkerke formulation).   
!!      3'- The cloud mixing length is computed according to HTURBLEN_CLOUD
!!             and emphasized following the CEI index
!!      4- The conservative variables are computed along with Lv/Cp.
!!      5- The turbulent Prandtl numbers are computed from the resolved fields
!!         and TKE 
!!      6- The sources associated to the vertical turbulent fluxes are computed
!!      with a temporal scheme allowing a degree of implicitness given by 
!!      PIMPL, varying from PIMPL=0. ( purely explicit scheme) to PIMPL=1.
!!      ( purely implicit scheme)
!!      The sources associated to the horizontal fluxes are computed with a
!!      purely explicit temporal scheme. These sources are only computed when
!!      the turbulence parameterization is 2D or 3D( HTURBDIM='3DIM' ).
!!      7- The sources for TKE are computed, along with the dissipation of TKE 
!!      if HTURBLEN='KEPS'.
!!      8- Some turbulence-related quantities are stored in the synchronous 
!!      FM-file.
!!      9- The non-conservative variables are retrieved.  
!!    
!!      
!!      The saving of the fields in the synchronous FM-file is controlled by:
!!        * OTURB_FLX => saves all the turbulent fluxes and correlations
!!        * OTURB_DIAG=> saves the turbulent Prandtl and Schmidt numbers, the
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
!!                                              turbulence (OHARAT=TRUE)
!!                     04/2016  (C.Lac) correction of negativity for KHKO
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  Q. Rodier      01/2018: introduction of RM17
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!  B. Vie         03/2020: LIMA negativity checks after turbulence, advection and microphysics budgets
!  P. Wautelet 11/06/2020: bugfix: correct PRSVS array indices
!  P. Wautelet + Benoit Vié 06/2020: improve removal of negative scalar variables + adapt the corresponding budgets
!  P. Wautelet 30/06/2020: move removal of negative scalar variables to Sources_neg_correct
!  R. Honnert/V. Masson 02/2021: new mixing length in the grey zone
!  J.L. Redelsperger 03/2021: add Ocean LES case
! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
USE MODD_PARAMETERS, ONLY: JPVEXT_TURB, XUNDEF
USE MODD_CST, ONLY: CST_t
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_BUDGET, ONLY:      NBUDGET_U,  NBUDGET_V,  NBUDGET_W,  NBUDGET_TH, NBUDGET_RV, NBUDGET_RC,  &
                            NBUDGET_RR, NBUDGET_RI, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, NBUDGET_SV1, &
                            TBUDGETDATA, TBUDGETCONF_t
USE MODD_FIELD, ONLY: TFIELDDATA,TYPEREAL
USE MODD_IO, ONLY: TFILEDATA
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
USE MODD_LES
USE MODD_IBM_PARAM_n,    ONLY: LIBM, XIBM_LS, XIBM_XMUT
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_TURB_n, ONLY: TURB_t
!
USE MODE_BL89, ONLY: BL89
USE MODE_TURB_VER, ONLY : TURB_VER
USE MODE_ROTATE_WIND, ONLY: ROTATE_WIND
USE MODE_TURB_HOR_SPLT, ONLY: TURB_HOR_SPLT
USE MODE_TKE_EPS_SOURCES, ONLY: TKE_EPS_SOURCES
USE MODE_RMC01, ONLY: RMC01
USE MODE_TM06, ONLY: TM06
USE MODE_UPDATE_LM, ONLY: UPDATE_LM
USE MODE_BUDGET,         ONLY: BUDGET_STORE_INIT, BUDGET_STORE_END
USE MODE_IO_FIELD_WRITE, ONLY: IO_FIELD_WRITE
USE MODE_ll, ONLY: ADD2DFIELD_ll, UPDATE_HALO_ll, CLEANLIST_ll, &
                   LWEST_ll, LEAST_ll, LSOUTH_ll, LNORTH_ll
USE MODE_SBL, ONLY: LMO
USE MODE_SOURCES_NEG_CORRECT, ONLY: SOURCES_NEG_CORRECT
USE MODE_EMOIST, ONLY: EMOIST
USE MODE_ETHETA, ONLY: ETHETA
!
USE MODI_GRADIENT_W
USE MODI_GRADIENT_M
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_IBM_MIXINGLENGTH
USE MODI_LES_MEAN_SUBGRID
USE MODI_SHUMAN, ONLY : MZF, MXF, MYF
USE SHUMAN_PHY, ONLY : MZF_PHY
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
TYPE(TBUDGETCONF_t),    INTENT(IN)   :: BUCONF
TYPE(TURB_t),           INTENT(IN)   :: TURBN
INTEGER,                INTENT(IN)   :: KMI           ! model index number  
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV, KSV_LGBEG, KSV_LGEND ! number of scalar variables
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
INTEGER,                INTENT(IN)   :: KSPLIT        ! number of time-splitting
INTEGER,                INTENT(IN)   :: KMODEL_CL     ! model number for cloud mixing length
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                 ! diagnostic fields in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OSUBG_COND   ! switch for SUBGrid CONDensation
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables
LOGICAL,                INTENT(IN)   ::  ORMC01       ! switch for RMC01 lengths in SBL
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  ODEEPOC      ! activates sfc forcing for ideal ocean deep conv
LOGICAL,                INTENT(IN)   ::  OHARAT       ! switch for LHARATU from AROME
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
LOGICAL,                INTENT(IN)   ::  OLES_CALL    ! compute the LES diagnostics at current time-step
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES version 
LOGICAL,                INTENT(IN)   ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
LOGICAL,                INTENT(IN)   ::  ODIAG_IN_RUN ! switch to activate online diagnostics (mesonh)
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBLEN     ! kind of mixing length
CHARACTER(LEN=4),       INTENT(IN)   ::  HTOM         ! kind of Third Order Moment
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBLEN_CL  ! kind of cloud mixing length
REAL,                   INTENT(IN)   ::  PIMPL        ! degree of implicitness
CHARACTER (LEN=4),      INTENT(IN)   ::  HCLOUD       ! Kind of microphysical scheme
REAL,                   INTENT(IN)   ::  PTSTEP       ! timestep 
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   :: PDXX,PDYY,PDZZ,PDZX,PDZY
                                        ! metric coefficients
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   :: PZZ       !  physical distance 
! between 2 succesive grid points along the K direction
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)      ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)      ::  PRHODJ    ! dry density * Grid size
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)      ::  MFMOIST ! moist mass flux dual scheme
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)      ::  PTHVREF   ! Virtual Potential
                                        ! Temperature of the reference state
!
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)      ::  PSFTH,PSFRV,   &
! normal surface fluxes of theta and Rv 
                                            PSFU,PSFV
! normal surface fluxes of (u,v) parallel to the orography
REAL, DIMENSION(D%NIT,D%NJT,KSV), INTENT(IN)      ::  PSFSV
! normal surface fluxes of Scalar var. 
!
!    prognostic variables at t- deltat
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) ::  PPABST      ! Pressure at time t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) ::  PUT,PVT,PWT ! wind components
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) ::  PTKET       ! TKE
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN) ::  PSVT        ! passive scal. var.
REAL, DIMENSION(MERGE(D%NIT,0,OCOMPUTE_SRC),&
                MERGE(D%NJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)),   INTENT(IN) ::  PSRCT       ! Second-order flux
                      ! s'rc'/2Sigma_s2 at time t-1 multiplied by Lambda_3
REAL, DIMENSION(MERGE(D%NIT,0,HTOM=='TM06'),&
                MERGE(D%NJT,0,HTOM=='TM06')),INTENT(INOUT) :: PBL_DEPTH  ! BL height for TOMS
REAL, DIMENSION(MERGE(D%NIT,0,ORMC01),&
                MERGE(D%NJT,0,ORMC01)),INTENT(INOUT) :: PSBL_DEPTH ! SBL depth for RMC01
!
!    variables for cloud mixing length
REAL, DIMENSION(MERGE(D%NIT,0,KMODEL_CL==KMI .AND. HTURBLEN_CL/='NONE'),&
                MERGE(D%NJT,0,KMODEL_CL==KMI .AND. HTURBLEN_CL/='NONE'),&
                MERGE(D%NKT,0,KMODEL_CL==KMI .AND. HTURBLEN_CL/='NONE')),INTENT(IN)      ::  PCEI 
                                                 ! Cloud Entrainment instability
                                                 ! index to emphasize localy 
                                                 ! turbulent fluxes
REAL, INTENT(IN)      ::  PCEI_MIN ! minimum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCEI_MAX ! maximum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCOEF_AMPL_SAT ! saturation of the amplification coefficient
!
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PTHLT       ! conservative pot. temp.
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(INOUT) ::  PRT         ! water var.  where 
                             ! PRT(:,:,:,1) is the conservative mixing ratio        
!
! sources of momentum, conservative potential temperature, Turb. Kin. Energy, 
! TKE dissipation
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PRUS,PRVS,PRWS,PRTHLS,PRTKES
! Source terms for all water kinds, PRRS(:,:,:,1) is used for the conservative
! mixing ratio
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN),OPTIONAL    ::  PRTKEMS
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(INOUT) ::  PRRS
! Source terms for all passive scalar variables
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(INOUT) ::  PRSVS
! Sigma_s at time t+1 : square root of the variance of the deviation to the 
! saturation
REAL, DIMENSION(MERGE(D%NIT,0,OCOMPUTE_SRC),&
                MERGE(D%NJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(OUT)     ::  PSIGS
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRUS_TURB   ! evolution of rhoJ*U   by turbulence only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRVS_TURB   ! evolution of rhoJ*V   by turbulence only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRTHLS_TURB ! evolution of rhoJ*thl by turbulence only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRRTS_TURB  ! evolution of rhoJ*rt  by turbulence only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(OUT),OPTIONAL ::  PDRSVS_TURB  ! evolution of rhoJ*Sv  by turbulence only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)      ::  PFLXZTHVMF 
!                                           MF contribution for vert. turb. transport
!                                           used in the buoy. prod. of TKE
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PWTH       ! heat flux
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PWRC       ! cloud water flux
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV),INTENT(OUT) :: PWSV       ! scalar flux
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PTP        ! Thermal TKE production
                                                   ! MassFlux + turb
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT),OPTIONAL  :: PTPMF      ! Thermal TKE production
                                                   ! MassFlux Only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PDP        ! Dynamic TKE production
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PTDIFF     ! Diffusion TKE term
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PTDISS     ! Dissipation TKE term
!
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! CPROGRAM is the program currently running (modd_conf)
LOGICAL, INTENT(IN) :: ONOMIXLG          ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL, INTENT(IN) :: O2D               ! Logical for 2D model version (modd_conf)
!
! length scale from vdfexcu
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PLENGTHM, PLENGTHH
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT), OPTIONAL  :: PEDR  ! EDR
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT), OPTIONAL  :: PLEM  ! Mixing length
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(OUT), OPTIONAL  ::  PTR          ! Transport prod. of TKE
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(OUT), OPTIONAL  ::  PDISS        ! Dissipation of TKE
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT), OPTIONAL  ::  PCURRENT_TKE_DISS ! if ODIAG_IN_RUN in mesonh
!
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) ::     &
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
          ZWORK1,                     &  ! working array syntax
          ZETHETA,ZEMOIST,            &  ! coef ETHETA and EMOIST (for DEAR routine)
          ZDTHLDZ,ZDRTDZ,             &  ! dtheta_l/dz, drt_dz used for computing the stablity criterion
          ZCOEF_AMPL,                 &  ! Amplification coefficient of the mixing length
                                         ! when the instability criterium is verified (routine CLOUD_MODIF_LM)
          ZLM_CLOUD                      ! Turbulent mixing length in the clouds (routine CLOUD_MODIF_LM)
!
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR) :: ZRM ! initial mixing ratio 
REAL, DIMENSION(D%NIT,D%NJT) ::  ZTAU11M,ZTAU12M,  &
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
REAL                :: ZEXPL        ! 1-PIMPL deg of expl.
REAL                :: ZRVORD       ! RV/RD
REAL                :: ZEPS         ! XMV / XMD
REAL                :: ZD           ! distance to the surface (for routine DELT)
REAL                :: ZVAR         ! Intermediary variable (for routine DEAR)
REAL                :: ZPENTE       ! Slope of the amplification straight line (for routine CLOUD_MODIF_LM)
REAL                :: ZCOEF_AMPL_CEI_NUL! Ordonnate at the origin of the
                                         ! amplification straight line (for routine CLOUD_MODIF_LM)
!
INTEGER             :: IIE,IIB,IJE,IJB,IKB,IKE      ! index value for the
INTEGER             :: IINFO_ll     ! return code of parallel routine
! Beginning and the End of the physical domain for the mass points
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain 
INTEGER             :: JRR,JK,JSV   ! loop counters
INTEGER             :: JI,JJ        ! loop counters
REAL                :: ZL0          ! Max. Mixing Length in Blakadar formula
REAL                :: ZALPHA       ! work coefficient : 
                                    ! - proportionnality constant between Dz/2 and 
!                                   !   BL89 mixing length near the surface
!
REAL :: ZTIME1, ZTIME2
TYPE(TFIELDDATA) :: TZFIELD
TYPE(LIST_ll), POINTER :: TZFIELDS_ll  ! list of fields to exchange (for UPDATE_ROTATE_WIND)
!
!*      1.PRELIMINARIES
!         -------------
!
!*      1.1 Set the internal domains, ZEXPL 
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE,ZHOOK_HANDLE2
IF (LHOOK) CALL DR_HOOK('TURB',0,ZHOOK_HANDLE)
!
IF (OHARAT .AND. HTURBDIM /= '1DIM') THEN
  CALL ABOR1('OHARATU only implemented for option HTURBDIM=1DIM!')
ENDIF
IF (OHARAT .AND. OLES_CALL) THEN
  CALL ABOR1('OHARATU not implemented for option LLES_CALL')
ENDIF
!
IKT=D%NKT  
IKTB=D%NKTB          
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
!
ZEXPL = 1.- PIMPL
ZRVORD= CST%XRV / CST%XRD
!
!Copy data into ZTHLM and ZRM only if needed
IF (HTURBLEN=='BL89' .OR. HTURBLEN=='RM17' .OR. HTURBLEN=='ADAP' .OR. ORMC01) THEN
  ZTHLM(IIB:IIE,IJB:IJE,1:D%NKT) = PTHLT(IIB:IIE,IJB:IJE,1:D%NKT)
  ZRM(IIB:IIE,IJB:IJE,1:D%NKT,:) = PRT(IIB:IIE,IJB:IJE,1:D%NKT,:)
END IF
!
!----------------------------------------------------------------------------
!
!*      2. COMPUTE CONSERVATIVE VARIABLES AND RELATED QUANTITIES
!          -----------------------------------------------------
!
!*      2.1 Cph at t
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
ZCP(IIB:IIE,IJB:IJE,1:D%NKT)=CST%XCPD
!
IF (KRR > 0) ZCP(IIB:IIE,IJB:IJE,1:D%NKT) = ZCP(IIB:IIE,IJB:IJE,1:D%NKT) + CST%XCPV * PRT(IIB:IIE,IJB:IJE,1:D%NKT,1)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
DO JRR = 2,1+KRRL                          ! loop on the liquid components
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)  
  ZCP(IIB:IIE,IJB:IJE,1:D%NKT)  = ZCP(IIB:IIE,IJB:IJE,1:D%NKT) + CST%XCL * PRT(IIB:IIE,IJB:IJE,1:D%NKT,JRR)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
END DO
!
DO JRR = 2+KRRL,1+KRRL+KRRI                ! loop on the solid components   
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  ZCP(IIB:IIE,IJB:IJE,1:D%NKT)  = ZCP(IIB:IIE,IJB:IJE,1:D%NKT)  + CST%XCI * PRT(IIB:IIE,IJB:IJE,1:D%NKT,JRR)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
END DO
!
!*      2.2 Exner function at t
!
IF (OOCEAN) THEN
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  ZEXN(IIB:IIE,IJB:IJE,1:D%NKT) = 1.
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
ELSE
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  ZEXN(IIB:IIE,IJB:IJE,1:D%NKT) = (PPABST(IIB:IIE,IJB:IJE,1:D%NKT)/CST%XP00) ** (CST%XRD/CST%XCPD)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
END IF
!
!*      2.3 dissipative heating coeff a t
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
ZCOEF_DISS(IIB:IIE,IJB:IJE,1:D%NKT) = 1/(ZCP(IIB:IIE,IJB:IJE,1:D%NKT) * ZEXN(IIB:IIE,IJB:IJE,1:D%NKT)) 
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
!
!
ZFRAC_ICE(IIB:IIE,IJB:IJE,1:D%NKT) = 0.0
ZATHETA(IIB:IIE,IJB:IJE,1:D%NKT) = 0.0
ZAMOIST(IIB:IIE,IJB:IJE,1:D%NKT) = 0.0
!
IF (KRRL >=1) THEN
!
!*      2.4 Temperature at t
!
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  ZT(IIB:IIE,IJB:IJE,1:D%NKT) =  PTHLT(IIB:IIE,IJB:IJE,1:D%NKT) * ZEXN(IIB:IIE,IJB:IJE,1:D%NKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
!
!*       2.5 Lv/Cph/Exn
!
  IF ( KRRI >= 1 ) THEN 
    CALL COMPUTE_FUNCTION_THERMO(CST%XALPW,CST%XBETAW,CST%XGAMW,CST%XLVTT,CST%XCL,ZT,ZEXN,ZCP, &
                                 ZLVOCPEXNM,ZAMOIST,ZATHETA)
    CALL COMPUTE_FUNCTION_THERMO(CST%XALPI,CST%XBETAI,CST%XGAMI,CST%XLSTT,CST%XCI,ZT,ZEXN,ZCP, &
                                 ZLSOCPEXNM,ZAMOIST_ICE,ZATHETA_ICE)
!
    !$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
    WHERE(PRT(IIB:IIE,IJB:IJE,1:D%NKT,2)+PRT(IIB:IIE,IJB:IJE,1:D%NKT,4)>0.0)
      ZFRAC_ICE(IIB:IIE,IJB:IJE,1:D%NKT) = PRT(IIB:IIE,IJB:IJE,1:D%NKT,4) / ( PRT(IIB:IIE,IJB:IJE,1:D%NKT,2) & 
                                          +PRT(IIB:IIE,IJB:IJE,1:D%NKT,4) )
    END WHERE
    !$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
!
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
    ZLOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) = (1.0-ZFRAC_ICE(IIB:IIE,IJB:IJE,1:D%NKT))*ZLVOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) &
                           +ZFRAC_ICE(IIB:IIE,IJB:IJE,1:D%NKT) *ZLSOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT)
    ZAMOIST(IIB:IIE,IJB:IJE,1:D%NKT) = (1.0-ZFRAC_ICE(IIB:IIE,IJB:IJE,1:D%NKT))*ZAMOIST(IIB:IIE,IJB:IJE,1:D%NKT) &
                         +ZFRAC_ICE(IIB:IIE,IJB:IJE,1:D%NKT) *ZAMOIST_ICE(IIB:IIE,IJB:IJE,1:D%NKT)
    ZATHETA(IIB:IIE,IJB:IJE,1:D%NKT) = (1.0-ZFRAC_ICE(IIB:IIE,IJB:IJE,1:D%NKT))*ZATHETA(IIB:IIE,IJB:IJE,1:D%NKT) &
                         +ZFRAC_ICE(IIB:IIE,IJB:IJE,1:D%NKT) *ZATHETA_ICE(IIB:IIE,IJB:IJE,1:D%NKT)
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  ELSE
    CALL COMPUTE_FUNCTION_THERMO(CST%XALPW,CST%XBETAW,CST%XGAMW,CST%XLVTT,CST%XCL,ZT,ZEXN,ZCP, &
                                 ZLOCPEXNM,ZAMOIST,ZATHETA)
  END IF
!
!
  IF ( TPFILE%LOPENED .AND. OTURB_DIAG ) THEN
    TZFIELD%CMNHNAME   = 'ATHETA'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'ATHETA'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_ATHETA'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZATHETA)
! 
    TZFIELD%CMNHNAME   = 'AMOIST'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'AMOIST'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_AMOIST'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZAMOIST)
  END IF
!
ELSE
  ZLOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT)=0.
END IF              ! loop end on KRRL >= 1
!
! computes conservative variables
!
IF ( KRRL >= 1 ) THEN
  IF ( KRRI >= 1 ) THEN
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
    ! Rnp at t
    PRT(IIB:IIE,IJB:IJE,1:D%NKT,1)  = PRT(IIB:IIE,IJB:IJE,1:D%NKT,1)  + PRT(IIB:IIE,IJB:IJE,1:D%NKT,2)  & 
                                    + PRT(IIB:IIE,IJB:IJE,1:D%NKT,4)
    PRRS(IIB:IIE,IJB:IJE,1:D%NKT,1) = PRRS(IIB:IIE,IJB:IJE,1:D%NKT,1) + PRRS(IIB:IIE,IJB:IJE,1:D%NKT,2) & 
                                    + PRRS(IIB:IIE,IJB:IJE,1:D%NKT,4)
    ! Theta_l at t
    PTHLT(IIB:IIE,IJB:IJE,1:D%NKT)  = PTHLT(IIB:IIE,IJB:IJE,1:D%NKT)  - ZLVOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) &
                                    * PRT(IIB:IIE,IJB:IJE,1:D%NKT,2) &
                                  - ZLSOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) * PRT(IIB:IIE,IJB:IJE,1:D%NKT,4)
    PRTHLS(IIB:IIE,IJB:IJE,1:D%NKT) = PRTHLS(IIB:IIE,IJB:IJE,1:D%NKT) - ZLVOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) &
                                    * PRRS(IIB:IIE,IJB:IJE,1:D%NKT,2) &
                                  - ZLSOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) * PRRS(IIB:IIE,IJB:IJE,1:D%NKT,4)
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  ELSE
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
    ! Rnp at t
    PRT(IIB:IIE,IJB:IJE,1:D%NKT,1)  = PRT(IIB:IIE,IJB:IJE,1:D%NKT,1)  + PRT(IIB:IIE,IJB:IJE,1:D%NKT,2) 
    PRRS(IIB:IIE,IJB:IJE,1:D%NKT,1) = PRRS(IIB:IIE,IJB:IJE,1:D%NKT,1) + PRRS(IIB:IIE,IJB:IJE,1:D%NKT,2)
    ! Theta_l at t
    PTHLT(IIB:IIE,IJB:IJE,1:D%NKT)  = PTHLT(IIB:IIE,IJB:IJE,1:D%NKT)  - ZLOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) &
                                    * PRT(IIB:IIE,IJB:IJE,1:D%NKT,2)
    PRTHLS(IIB:IIE,IJB:IJE,1:D%NKT) = PRTHLS(IIB:IIE,IJB:IJE,1:D%NKT) - ZLOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) &
                                    * PRRS(IIB:IIE,IJB:IJE,1:D%NKT,2)
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  END IF
END IF
!
!* stores value of conservative variables & wind before turbulence tendency (AROME diag)
IF(PRESENT(PDRUS_TURB)) THEN
  PDRUS_TURB = PRUS
  PDRVS_TURB = PRVS
  PDRTHLS_TURB = PRTHLS
  PDRRTS_TURB  = PRRS(:,:,:,1)
  PDRSVS_TURB  = PRSVS
END IF
!----------------------------------------------------------------------------
!
!*      3. MIXING LENGTH : SELECTION AND COMPUTATION
!          -----------------------------------------
!
!
IF (.NOT. OHARAT) THEN

SELECT CASE (HTURBLEN)
!
!*      3.1 BL89 mixing length
!           ------------------

  CASE ('BL89')
    ZSHEAR(:,:,:)=0.
    CALL BL89(D,CST,CSTURB,PZZ,PDZZ,PTHVREF,ZTHLM,KRR,ZRM,PTKET,ZSHEAR,ZLM,OOCEAN,HPROGRAM)
!
!*      3.2 RM17 mixing length
!           ------------------

  CASE ('RM17')
    ZDUDZ = MXF(MZF(GZ_U_UW(PUT,PDZZ,D%NKA,D%NKU,D%NKL),D%NKA,D%NKU,D%NKL))
    ZDVDZ = MYF(MZF(GZ_V_VW(PVT,PDZZ,D%NKA,D%NKU,D%NKL),D%NKA,D%NKU,D%NKL))
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    ZSHEAR(:,:,:) = SQRT(ZDUDZ(:,:,:)*ZDUDZ(:,:,:) + ZDVDZ(:,:,:)*ZDVDZ(:,:,:))
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    CALL BL89(D,CST,CSTURB,PZZ,PDZZ,PTHVREF,ZTHLM,KRR,ZRM,PTKET,ZSHEAR,ZLM,OOCEAN,HPROGRAM)
!
!*      3.3 Grey-zone combined RM17 & Deardorff mixing lengths 
!           --------------------------------------------------

  CASE ('ADAP')
    ZDUDZ = MXF(MZF(GZ_U_UW(PUT,PDZZ,D%NKA,D%NKU,D%NKL),D%NKA,D%NKU,D%NKL))
    ZDVDZ = MYF(MZF(GZ_V_VW(PVT,PDZZ,D%NKA,D%NKU,D%NKL),D%NKA,D%NKU,D%NKL))
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    ZSHEAR(:,:,:) = SQRT(ZDUDZ(:,:,:)*ZDUDZ(:,:,:) + ZDVDZ(:,:,:)*ZDVDZ(:,:,:))
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    CALL BL89(D,CST,CSTURB,PZZ,PDZZ,PTHVREF,ZTHLM,KRR,ZRM,PTKET,ZSHEAR,ZLM,OOCEAN,HPROGRAM)

    CALL DELT(ZLMW,ODZ=.FALSE.)
    ! The minimum mixing length is chosen between Horizontal grid mesh (not taking into account the vertical grid mesh) and RM17.
    ! For large horizontal grid meshes, this is equal to RM17
    ! For LES grid meshes, this is equivalent to Deardorff : the base mixing lentgh is the horizontal grid mesh, 
    !                      and it is limited by a stability-based length (RM17), as was done in Deardorff length (but taking into account shear as well)
    ! For grid meshes in the grey zone, then this is the smaller of the two.
    !
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    ZLM(:,:,:) = MIN(ZLM(:,:,:),TURBN%XCADAP*ZLMW(:,:,:))
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
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
   ZLM(:,:,:) = ZL0

   ZALPHA=0.5**(-1.5)
   !
   DO JK=IKTB,IKTE
     !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
     ZLM(:,:,JK) = ( 0.5*(PZZ(:,:,JK)+PZZ(:,:,JK+D%NKL)) - &
     & PZZ(:,:,D%NKA+JPVEXT_TURB*D%NKL) ) * PDIRCOSZW(:,:)
     ZLM(:,:,JK) = ZALPHA  * ZLM(:,:,JK) * ZL0 / ( ZL0 + ZALPHA*ZLM(:,:,JK) )
     !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
   END DO
!
   ZLM(:,:,IKTB-1) = ZLM(:,:,IKTB)
   ZLM(:,:,IKTE+1) = ZLM(:,:,IKTE)
!
!
!
END SELECT
!
!*      3.5 Mixing length modification for cloud
!           -----------------------
IF (KMODEL_CL==KMI .AND. HTURBLEN_CL/='NONE') CALL CLOUD_MODIF_LM
ENDIF  ! end LHARRAT

!
!*      3.6 Dissipative length
!           ------------------

IF (OHARAT) THEN
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  ZLEPS(IIB:IIE,IJB:IJE,1:D%NKT)=PLENGTHM(IIB:IIE,IJB:IJE,1:D%NKT)*(3.75**2.)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
ELSE
  ZLEPS(IIB:IIE,IJB:IJE,1:D%NKT)=ZLM(IIB:IIE,IJB:IJE,1:D%NKT)
ENDIF
!
!*      3.7 Correction in the Surface Boundary Layer (Redelsperger 2001)
!           ----------------------------------------
!
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
ZLMO(:,:)=XUNDEF
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
IF (ORMC01) THEN
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
  ZUSTAR(:,:)=(PSFU(:,:)**2+PSFV(:,:)**2)**(0.25)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
  IF (KRR>0) THEN
    CALL LMO(ZUSTAR,ZTHLM(:,:,IKB),ZRM(:,:,IKB,1),PSFTH,PSFRV,ZLMO)
  ELSE
    ZRVM(:,:)=0.
    ZSFRV(:,:)=0.
    CALL LMO(ZUSTAR,ZTHLM(:,:,IKB),ZRVM,PSFTH,ZSFRV,ZLMO)
  END IF
  CALL RMC01(D,CST,CSTURB,HTURBLEN,PZZ,PDXX,PDYY,PDZZ,PDIRCOSZW,PSBL_DEPTH,ZLMO,ZLM,ZLEPS)
END IF
!
!RMC01 is only applied on RM17 in ADAP
IF (HTURBLEN=='ADAP') THEN
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  ZLEPS(IIB:IIE,IJB:IJE,1:D%NKT) = MIN(ZLEPS(IIB:IIE,IJB:IJE,1:D%NKT),ZLMW(IIB:IIE,IJB:IJE,1:D%NKT)*TURBN%XCADAP)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
END IF
!
!*      3.8 Mixing length in external points (used if HTURBDIM="3DIM")
!           ----------------------------------------------------------
!
IF (HTURBDIM=="3DIM") THEN
  CALL UPDATE_LM(HLBCX,HLBCY,ZLM,ZLEPS)
END IF
!
!*      3.9 Mixing length correction if immersed walls 
!           ------------------------------------------
!
IF (LIBM) THEN
   CALL IBM_MIXINGLENGTH(ZLM,ZLEPS,XIBM_XMUT,XIBM_LS(:,:,:,1),PTKET)
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
IF (HPROGRAM/='AROME ') THEN
  CALL ROTATE_WIND(PUT,PVT,PWT,                       &
                     PDIRCOSXW, PDIRCOSYW, PDIRCOSZW,   &
                     PCOSSLOPE,PSINSLOPE,               &
                     PDXX,PDYY,PDZZ,                    &
                     ZUSLOPE,ZVSLOPE                    )
!
  CALL UPDATE_ROTATE_WIND(ZUSLOPE,ZVSLOPE)
ELSE
  ZUSLOPE=PUT(IIB:IIE,IJB:IJE,D%NKA)
  ZVSLOPE=PVT(IIB:IIE,IJB:IJE,D%NKA)
END IF
!
!
!*      4.2 compute the proportionality coefficient between wind and stress
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
ZCDUEFF(IIB:IIE,IJB:IJE) =-SQRT ( (PSFU(IIB:IIE,IJB:IJE)**2 + PSFV(IIB:IIE,IJB:IJE)**2) /               &
#ifdef REPRO48
                    (1.E-60 + ZUSLOPE(IIB:IIE,IJB:IJE)**2 + ZVSLOPE(IIB:IIE,IJB:IJE)**2 ) )
#else
                    (CST%XMNH_TINY + ZUSLOPE(IIB:IIE,IJB:IJE)**2 + ZVSLOPE(IIB:IIE,IJB:IJE)**2 ) )
#endif                      
!
!*       4.6 compute the surface tangential fluxes
!
ZTAU11M(IIB:IIE,IJB:IJE) =2./3.*(  (1.+ (PZZ(IIB:IIE,IJB:IJE,IKB+D%NKL)-PZZ(IIB:IIE,IJB:IJE,IKB))  &
                           /(PDZZ(IIB:IIE,IJB:IJE,IKB+D%NKL)+PDZZ(IIB:IIE,IJB:IJE,IKB))  &
                       )   *PTKET(IIB:IIE,IJB:IJE,IKB)                   &
                     -0.5  *PTKET(IIB:IIE,IJB:IJE,IKB+D%NKL)                 &
                    )
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
ZTAU12M(IIB:IIE,IJB:IJE) =0.0
ZTAU22M(IIB:IIE,IJB:IJE) =ZTAU11M(IIB:IIE,IJB:IJE)
ZTAU33M(IIB:IIE,IJB:IJE) =ZTAU11M(IIB:IIE,IJB:IJE)
!
!*       4.7 third order terms in temperature and water fluxes and correlations
!            ------------------------------------------------------------------
!
!
ZMWTH(:,:,:) = 0.     ! w'2th'
ZMWR(:,:,:)  = 0.     ! w'2r'
ZMTH2(:,:,:) = 0.     ! w'th'2
ZMR2(:,:,:)  = 0.     ! w'r'2
ZMTHR(:,:,:) = 0.     ! w'th'r'
!
IF (HTOM=='TM06') THEN
  CALL TM06(D%NKA,D%NKU,D%NKL,PTHVREF,PBL_DEPTH,PZZ,PSFTH,ZMWTH,ZMTH2)
!
  ZFWTH = -GZ_M_W(D%NKA,D%NKU,D%NKL,ZMWTH,PDZZ)    ! -d(w'2th' )/dz
  !ZFWR  = -GZ_M_W(D%NKA,D%NKU,D%NKL,ZMWR, PDZZ)    ! -d(w'2r'  )/dz
  ZFTH2 = -GZ_W_M(ZMTH2,PDZZ)    ! -d(w'th'2 )/dz
  !ZFR2  = -GZ_W_M(ZMR2, PDZZ)    ! -d(w'r'2  )/dz
  !ZFTHR = -GZ_W_M(ZMTHR,PDZZ)    ! -d(w'th'r')/dz
!
  ZFWTH(:,:,IKTE:) = 0.
  ZFWTH(:,:,:IKTB) = 0.
  !ZFWR (:,:,IKTE:) = 0.
  !ZFWR (:,:,:IKTB) = 0.
  ZFWR(:,:,:)  = 0.
  ZFTH2(:,:,IKTE:) = 0.
  ZFTH2(:,:,:IKTB) = 0.
  !ZFR2 (:,:,IKTE:) = 0.
  !ZFR2 (:,:,:IKTB) = 0.
  ZFR2(:,:,:)  = 0.
  !ZFTHR(:,:,IKTE:) = 0.
  !ZFTHR(:,:,:IKTB) = 0.
  ZFTHR(:,:,:) = 0.
ELSE
  ZFWTH(:,:,:) = 0.
  ZFWR(:,:,:)  = 0.
  ZFTH2(:,:,:) = 0.
  ZFR2(:,:,:)  = 0.
  ZFTHR(:,:,:) = 0.
ENDIF
!
!----------------------------------------------------------------------------
!
!*      5. TURBULENT SOURCES
!          -----------------
!
IF( BUCONF%LBUDGET_U )  CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_U ), 'VTURB', PRUS(:, :, :)    )
IF( BUCONF%LBUDGET_V )  CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_V ), 'VTURB', PRVS(:, :, :)    )
IF( BUCONF%LBUDGET_W )  CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_W ), 'VTURB', PRWS(:, :, :)    )

IF( BUCONF%LBUDGET_TH ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:, :, :) + ZLVOCPEXNM(:, :, :) * PRRS(:, :, :, 2) &
                                                                          + ZLSOCPEXNM(:, :, :) * PRRS(:, :, :, 4) )
  ELSE IF( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:, :, :) + ZLOCPEXNM(:, :, :) * PRRS(:, :, :, 2) )
  ELSE
    CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:, :, :) )
  END IF
END IF

IF( BUCONF%LBUDGET_RV ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:, :, :, 1) - PRRS(:, :, :, 2) - PRRS(:, :, :, 4) )
  ELSE IF( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:, :, :, 1) - PRRS(:, :, :, 2) )
  ELSE
    CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:, :, :, 1) )
  END IF
END IF

IF( BUCONF%LBUDGET_RC ) CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_RC), 'VTURB', PRRS  (:, :, :, 2) )
IF( BUCONF%LBUDGET_RI ) CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_RI), 'VTURB', PRRS  (:, :, :, 4) )

IF( BUCONF%LBUDGET_SV ) THEN
  DO JSV = 1, KSV
    CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_SV1 - 1 + JSV), 'VTURB', PRSVS(:, :, :, JSV) )
  END DO
END IF

CALL TURB_VER(D, CST,CSTURB,TURBN,KRR, KRRL, KRRI,       &
          OTURB_FLX, OOCEAN, ODEEPOC, OHARAT,OCOMPUTE_SRC,&
          KSV,KSV_LGBEG,KSV_LGEND,                       &
          HTURBDIM,HTOM,PIMPL,ZEXPL,                     &
          HPROGRAM, O2D, ONOMIXLG, OFLAT,                &
          OLES_CALL,OCOUPLES,OBLOWSNOW, ORMC01,          &
          PTSTEP,TPFILE,                                 &
          PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,        &
          PCOSSLOPE,PSINSLOPE,                           &
          PRHODJ,PTHVREF,                                &
          PSFTH,PSFRV,PSFSV,PSFTH,PSFRV,PSFSV,           &
          ZCDUEFF,ZTAU11M,ZTAU12M,ZTAU33M,               &
          PUT,PVT,PWT,ZUSLOPE,ZVSLOPE,PTHLT,PRT,PSVT,    &
          PTKET,ZLM,PLENGTHM,PLENGTHH,ZLEPS,MFMOIST,     &
          ZLOCPEXNM,ZATHETA,ZAMOIST,PSRCT,ZFRAC_ICE,     &
          ZFWTH,ZFWR,ZFTH2,ZFR2,ZFTHR,PBL_DEPTH,         &
          PSBL_DEPTH,ZLMO,                               &
          PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,              &
          PDP,PTP,PSIGS,PWTH,PWRC,PWSV                   )

IF( BUCONF%LBUDGET_U ) CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_U), 'VTURB', PRUS(:, :, :) )
IF( BUCONF%LBUDGET_V ) CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_V), 'VTURB', PRVS(:, :, :) )
IF( BUCONF%LBUDGET_W ) CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_W), 'VTURB', PRWS(:, :, :) )

IF( BUCONF%LBUDGET_TH ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:, :, :) + ZLVOCPEXNM(:, :, :) * PRRS(:, :, :, 2) &
                                                                          + ZLSOCPEXNM(:, :, :) * PRRS(:, :, :, 4) )
  ELSE IF( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:, :, :) + ZLOCPEXNM(:, :, :) * PRRS(:, :, :, 2) )
  ELSE
    CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_TH), 'VTURB', PRTHLS(:, :, :) )
  END IF
END IF

IF( BUCONF%LBUDGET_RV ) THEN
  IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:, :, :, 1) - PRRS(:, :, :, 2) - PRRS(:, :, :, 4) )
  ELSE IF( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:, :, :, 1) - PRRS(:, :, :, 2) )
  ELSE
    CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_RV), 'VTURB', PRRS(:, :, :, 1) )
  END IF
END IF

IF( BUCONF%LBUDGET_RC ) CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_RC), 'VTURB', PRRS(:, :, :, 2) )
IF( BUCONF%LBUDGET_RI ) CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_RI), 'VTURB', PRRS(:, :, :, 4) )

IF( BUCONF%LBUDGET_SV )  THEN
  DO JSV = 1, KSV
    CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_SV1 - 1 + JSV), 'VTURB', PRSVS(:, :, :, JSV) )
  END DO
END IF
!
!Les budgets des termes horizontaux de la turb sont présents dans AROME
! alors que ces termes ne sont pas calculés
#ifdef REPRO48 
#else          
IF( HTURBDIM == '3DIM' ) THEN
#endif
  IF( BUCONF%LBUDGET_U  ) CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_U ), 'HTURB', PRUS  (:, :, :) )
  IF( BUCONF%LBUDGET_V  ) CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_V ), 'HTURB', PRVS  (:, :, :) )
  IF( BUCONF%LBUDGET_W  ) CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_W ), 'HTURB', PRWS  (:, :, :) )

  IF(BUCONF%LBUDGET_TH)  THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
      CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:, :, :) + ZLVOCPEXNM(:, :, :) * PRRS(:, :, :, 2) &
                                                                             + ZLSOCPEXNM(:, :, :) * PRRS(:, :, :, 4) )
    ELSE IF( KRRL >= 1 ) THEN
      CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:, :, :) + ZLOCPEXNM(:, :, :) * PRRS(:, :, :, 2) )
    ELSE
      CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:, :, :) )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RV ) THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
      CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:, :, :, 1) - PRRS(:, :, :, 2) - PRRS(:, :, :, 4) )
    ELSE IF( KRRL >= 1 ) THEN
      CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:, :, :, 1) - PRRS(:, :, :, 2) )
    ELSE
      CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:, :, :, 1) )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RC ) CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_RC), 'HTURB', PRRS(:, :, :, 2) )
  IF( BUCONF%LBUDGET_RI ) CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_RI), 'HTURB', PRRS(:, :, :, 4) )

  IF( BUCONF%LBUDGET_SV )  THEN
    DO JSV = 1, KSV
      CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_SV1 - 1 + JSV), 'HTURB', PRSVS(:, :, :, JSV) )
    END DO
  END IF
!à supprimer une fois le précédent ifdef REPRO48 validé
#ifdef REPRO48
#else
    CALL TURB_HOR_SPLT(D,CST,CSTURB,                           &
          KSPLIT, KRR, KRRL, KRRI, PTSTEP,HLBCX,HLBCY,         &
          OTURB_FLX,OSUBG_COND,OOCEAN,OCOMPUTE_SRC,            &
          TPFILE,                                              &
          PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,                        &
          PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,                       &
          PCOSSLOPE,PSINSLOPE,                                 &
          PRHODJ,PTHVREF,                                      &
          PSFTH,PSFRV,PSFSV,                                   &
          ZCDUEFF,ZTAU11M,ZTAU12M,ZTAU22M,ZTAU33M,             &
          PUT,PVT,PWT,ZUSLOPE,ZVSLOPE,PTHLT,PRT,PSVT,          &
          PTKET,ZLM,ZLEPS,                                     &
          ZLOCPEXNM,ZATHETA,ZAMOIST,PSRCT,ZFRAC_ICE,           &
          PDP,PTP,PSIGS,                                       &
          ZTRH,                                                &
          PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS                     )
#endif
  IF( BUCONF%LBUDGET_U ) CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_U), 'HTURB', PRUS(:, :, :) )
  IF( BUCONF%LBUDGET_V ) CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_V), 'HTURB', PRVS(:, :, :) )
  IF( BUCONF%LBUDGET_W ) CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_W), 'HTURB', PRWS(:, :, :) )

  IF( BUCONF%LBUDGET_TH ) THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
      CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:, :, :) + ZLVOCPEXNM(:, :, :) * PRRS(:, :, :, 2) &
                                                                            + ZLSOCPEXNM(:, :, :) * PRRS(:, :, :, 4) )
    ELSE IF( KRRL >= 1 ) THEN
      CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:, :, :) + ZLOCPEXNM(:, :, :) * PRRS(:, :, :, 2) )
    ELSE
      CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_TH), 'HTURB', PRTHLS(:, :, :) )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RV ) THEN
    IF( KRRI >= 1 .AND. KRRL >= 1 ) THEN
      CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:, :, :, 1) - PRRS(:, :, :, 2) - PRRS(:, :, :, 4) )
    ELSE IF( KRRL >= 1 ) THEN
      CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:, :, :, 1) - PRRS(:, :, :, 2) )
    ELSE
      CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_RV), 'HTURB', PRRS(:, :, :, 1) )
    END IF
  END IF

  IF( BUCONF%LBUDGET_RC ) CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_RC), 'HTURB', PRRS(:, :, :, 2) )
  IF( BUCONF%LBUDGET_RI ) CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_RI), 'HTURB', PRRS(:, :, :, 4) )

  IF( BUCONF%LBUDGET_SV )  THEN
    DO JSV = 1, KSV
      CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_SV1 - 1 + JSV), 'HTURB', PRSVS(:, :, :, JSV) )
    END DO
  END IF
#ifdef REPRO48
#else
END IF
#endif
!----------------------------------------------------------------------------
!
!*      6. EVOLUTION OF THE TKE AND ITS DISSIPATION 
!          ----------------------------------------
!
!  6.1 Contribution of mass-flux in the TKE buoyancy production if 
!      cloud computation is not statistical 
CALL MZF_PHY(D,PFLXZTHVMF,ZWORK1)
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
PTP(IIB:IIE,IJB:IJE,1:D%NKT) = PTP(IIB:IIE,IJB:IJE,1:D%NKT) &
                             + CST%XG / PTHVREF(IIB:IIE,IJB:IJE,1:D%NKT) * ZWORK1(IIB:IIE,IJB:IJE,1:D%NKT)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)

IF(PRESENT(PTPMF))  THEN
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  PTPMF(IIB:IIE,IJB:IJE,1:D%NKT)=CST%XG / PTHVREF(IIB:IIE,IJB:IJE,1:D%NKT) * ZWORK1(IIB:IIE,IJB:IJE,1:D%NKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
END IF
!  6.2 TKE evolution equation

IF (.NOT. OHARAT) THEN
!
IF (BUCONF%LBUDGET_TH)  THEN
  IF ( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS+ ZLVOCPEXNM * PRRS(:,:,:,2) &
                                                          & + ZLSOCPEXNM * PRRS(:,:,:,4) )
  ELSE IF ( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS+ ZLOCPEXNM * PRRS(:,:,:,2) )
  ELSE
    CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS(:, :, :) )
  END IF
END IF
!
IF(PRESENT(PRTKEMS)) THEN
  ZRTKEMS(:,:,:)=PRTKEMS(:,:,:)
ELSE
  ZRTKEMS(:,:,:)=0.
END IF
!
CALL TKE_EPS_SOURCES(D,CST,CSTURB,BUCONF,HPROGRAM,                      &
                   & KMI,PTKET,ZLM,ZLEPS,PDP,ZTRH,                      &
                   & PRHODJ,PDZZ,PDXX,PDYY,PDZX,PDZY,PZZ,               &
                   & PTSTEP,PIMPL,ZEXPL,                                &
                   & HTURBLEN,HTURBDIM,                                 &
                   & TPFILE,OTURB_DIAG,OLES_CALL,ODIAG_IN_RUN,          &
                   & PTP,PRTKES,PRTHLS,ZCOEF_DISS,PTDIFF,PTDISS,ZRTKEMS,&
                   & TBUDGETS,KBUDGETS, PEDR=PEDR, PTR=PTR,PDISS=PDISS, &
                   & PCURRENT_TKE_DISS=PCURRENT_TKE_DISS                )
IF (BUCONF%LBUDGET_TH)  THEN
  IF ( KRRI >= 1 .AND. KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS+ ZLVOCPEXNM * PRRS(:,:,:,2) &
                                                          & + ZLSOCPEXNM * PRRS(:,:,:,4) )
  ELSE IF ( KRRL >= 1 ) THEN
    CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS+ ZLOCPEXNM * PRRS(:,:,:,2) )
  ELSE
    CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_TH), 'DISSH', PRTHLS(:, :, :) )
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
IF ( OTURB_DIAG .AND. TPFILE%LOPENED ) THEN
! 
! stores the mixing length
! 
  TZFIELD%CMNHNAME   = 'LM'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'LM'
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'Mixing length'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZLM)
!
  IF (KRR /= 0) THEN
!
! stores the conservative potential temperature
!
    TZFIELD%CMNHNAME   = 'THLM'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THLM'
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'Conservative potential temperature'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,PTHLT)
!
! stores the conservative mixing ratio
!
    TZFIELD%CMNHNAME   = 'RNPM'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'RNPM'
    TZFIELD%CUNITS     = 'kg kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'Conservative mixing ratio'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_FIELD_WRITE(TPFILE,TZFIELD,PRT(:,:,:,1))
   END IF
END IF
!
!* stores value of conservative variables & wind before turbulence tendency (AROME only)
IF(PRESENT(PDRUS_TURB)) THEN
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  PDRUS_TURB(IIB:IIE,IJB:IJE,1:D%NKT)   = PRUS(IIB:IIE,IJB:IJE,1:D%NKT) - PDRUS_TURB(IIB:IIE,IJB:IJE,1:D%NKT)
  PDRVS_TURB(IIB:IIE,IJB:IJE,1:D%NKT)   = PRVS(IIB:IIE,IJB:IJE,1:D%NKT) - PDRVS_TURB(IIB:IIE,IJB:IJE,1:D%NKT)
  PDRTHLS_TURB(IIB:IIE,IJB:IJE,1:D%NKT) = PRTHLS(IIB:IIE,IJB:IJE,1:D%NKT) - PDRTHLS_TURB(IIB:IIE,IJB:IJE,1:D%NKT)
  PDRRTS_TURB(IIB:IIE,IJB:IJE,1:D%NKT)  = PRRS(IIB:IIE,IJB:IJE,1:D%NKT,1) - PDRRTS_TURB(IIB:IIE,IJB:IJE,1:D%NKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT,JSV=1:KSV)  
  PDRSVS_TURB(IIB:IIE,IJB:IJE,1:D%NKT,:)  = PRSVS(IIB:IIE,IJB:IJE,1:D%NKT,:) - PDRSVS_TURB(IIB:IIE,IJB:IJE,1:D%NKT,:)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT,JSV=1:KSV)
END IF
!----------------------------------------------------------------------------
!
!*      8. RETRIEVE NON-CONSERVATIVE VARIABLES
!          -----------------------------------
!
IF ( KRRL >= 1 ) THEN
  IF ( KRRI >= 1 ) THEN
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
    PRT(IIB:IIE,IJB:IJE,1:D%NKT,1)  = PRT(IIB:IIE,IJB:IJE,1:D%NKT,1)  - PRT(IIB:IIE,IJB:IJE,1:D%NKT,2)  &
                                    - PRT(IIB:IIE,IJB:IJE,1:D%NKT,4)
    PRRS(IIB:IIE,IJB:IJE,1:D%NKT,1) = PRRS(IIB:IIE,IJB:IJE,1:D%NKT,1) - PRRS(IIB:IIE,IJB:IJE,1:D%NKT,2) &
                                    - PRRS(IIB:IIE,IJB:IJE,1:D%NKT,4)
    PTHLT(IIB:IIE,IJB:IJE,1:D%NKT)  = PTHLT(IIB:IIE,IJB:IJE,1:D%NKT)  + ZLVOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) &
                                    * PRT(IIB:IIE,IJB:IJE,1:D%NKT,2) &
                                    + ZLSOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) * PRT(IIB:IIE,IJB:IJE,1:D%NKT,4)
    PRTHLS(IIB:IIE,IJB:IJE,1:D%NKT) = PRTHLS(IIB:IIE,IJB:IJE,1:D%NKT) + ZLVOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) &
                                    * PRRS(IIB:IIE,IJB:IJE,1:D%NKT,2) &
                                    + ZLSOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) * PRRS(IIB:IIE,IJB:IJE,1:D%NKT,4)
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
!
  ELSE
    !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
    PRT(IIB:IIE,IJB:IJE,1:D%NKT,1)  = PRT(IIB:IIE,IJB:IJE,1:D%NKT,1)  - PRT(IIB:IIE,IJB:IJE,1:D%NKT,2) 
    PRRS(IIB:IIE,IJB:IJE,1:D%NKT,1) = PRRS(IIB:IIE,IJB:IJE,1:D%NKT,1) - PRRS(IIB:IIE,IJB:IJE,1:D%NKT,2)
    PTHLT(IIB:IIE,IJB:IJE,1:D%NKT)  = PTHLT(IIB:IIE,IJB:IJE,1:D%NKT)  + ZLOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) &
                                    * PRT(IIB:IIE,IJB:IJE,1:D%NKT,2)
    PRTHLS(IIB:IIE,IJB:IJE,1:D%NKT) = PRTHLS(IIB:IIE,IJB:IJE,1:D%NKT) + ZLOCPEXNM(IIB:IIE,IJB:IJE,1:D%NKT) &
                                    * PRRS(IIB:IIE,IJB:IJE,1:D%NKT,2)
    !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  END IF
END IF

! Remove non-physical negative values (unnecessary in a perfect world) + corresponding budgets
CALL SOURCES_NEG_CORRECT(HCLOUD, 'NETUR',KRR,PTSTEP,PPABST,PTHLT,PRT,PRTHLS,PRRS,PRSVS)
!----------------------------------------------------------------------------
!
!*      9. LES averaged surface fluxes
!          ---------------------------
!
IF (OLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID(PSFTH,X_LES_Q0)
  CALL LES_MEAN_SUBGRID(PSFRV,X_LES_E0)
  DO JSV=1,KSV
    CALL LES_MEAN_SUBGRID(PSFSV(:,:,JSV),X_LES_SV0(:,JSV))
  END DO
  CALL LES_MEAN_SUBGRID(PSFU,X_LES_UW0)
  CALL LES_MEAN_SUBGRID(PSFV,X_LES_VW0)
  CALL LES_MEAN_SUBGRID((PSFU*PSFU+PSFV*PSFV)**0.25,X_LES_USTAR)
!----------------------------------------------------------------------------
!
!*     10. LES for 3rd order moments
!          -------------------------
!
  CALL LES_MEAN_SUBGRID(ZMWTH,X_LES_SUBGRID_W2Thl)
  CALL LES_MEAN_SUBGRID(ZMTH2,X_LES_SUBGRID_WThl2)
  IF (KRR>0) THEN
    CALL LES_MEAN_SUBGRID(ZMWR,X_LES_SUBGRID_W2Rt)
    CALL LES_MEAN_SUBGRID(ZMTHR,X_LES_SUBGRID_WThlRt)
    CALL LES_MEAN_SUBGRID(ZMR2,X_LES_SUBGRID_WRt2)
  END IF
!
!----------------------------------------------------------------------------
!
!*     11. LES quantities depending on <w'2> in "1DIM" mode
!          ------------------------------------------------
!
  IF (HTURBDIM=="1DIM") THEN
    CALL LES_MEAN_SUBGRID(2./3.*PTKET,X_LES_SUBGRID_U2)
    X_LES_SUBGRID_V2(:,:,:) = X_LES_SUBGRID_U2(:,:,:)
    X_LES_SUBGRID_W2(:,:,:) = X_LES_SUBGRID_U2(:,:,:)
    CALL LES_MEAN_SUBGRID(2./3.*PTKET*MZF(GZ_M_W(D%NKA,D%NKU,D%NKL,PTHLT,PDZZ),&
                          D%NKA, D%NKU, D%NKL),X_LES_RES_ddz_Thl_SBG_W2)
    IF (KRR>=1) &
    CALL LES_MEAN_SUBGRID(2./3.*PTKET*MZF(GZ_M_W(D%NKA,D%NKU,D%NKL,PRT(:,:,:,1),PDZZ),&
                         &D%NKA, D%NKU, D%NKL),X_LES_RES_ddz_Rt_SBG_W2)
    DO JSV=1,KSV
      CALL LES_MEAN_SUBGRID(2./3.*PTKET*MZF(GZ_M_W(D%NKA,D%NKU,D%NKL,PSVT(:,:,:,JSV),PDZZ), &
                           &D%NKA, D%NKU, D%NKL), X_LES_RES_ddz_Sv_SBG_W2(:,:,:,JSV))
    END DO
  END IF

!----------------------------------------------------------------------------
!
!*     12. LES mixing end dissipative lengths, presso-correlations
!          -------------------------------------------------------
!
  CALL LES_MEAN_SUBGRID(ZLM,X_LES_SUBGRID_LMix)
  CALL LES_MEAN_SUBGRID(ZLEPS,X_LES_SUBGRID_LDiss)
!
!* presso-correlations for subgrid Tke are equal to zero.
!
  ZLEPS(:,:,:) = 0. !ZLEPS is used as a work array (not used anymore)
  CALL LES_MEAN_SUBGRID(ZLEPS,X_LES_SUBGRID_WP)
!
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
IF(PRESENT(PLEM)) PLEM(IIB:IIE,IJB:IJE,IKTB:IKTE) = ZLM(IIB:IIE,IJB:IJE,IKTB:IKTE)
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TURB',1,ZHOOK_HANDLE)
CONTAINS
!
!
!     ##############################################
      SUBROUTINE UPDATE_ROTATE_WIND(PUSLOPE,PVSLOPE)
!     ##############################################
!!
!!****  *UPDATE_ROTATE_WIND* routine to set rotate wind values at the border
!
!!    AUTHOR
!!    ------
!!
!!     P Jabouille   *CNRM METEO-FRANCE
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   24/06/99
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PUSLOPE,PVSLOPE
! tangential surface fluxes in the axes following the orography
!
IF (LHOOK) CALL DR_HOOK('TURB:UPDATE_ROTATE_WIND',0,ZHOOK_HANDLE2)
!
!*        1  PROLOGUE
!
NULLIFY(TZFIELDS_ll)
!
!         2 Update halo if necessary
!
!!$IF (NHALO == 1) THEN
  CALL ADD2DFIELD_ll( TZFIELDS_ll, PUSLOPE, 'UPDATE_ROTATE_WIND::PUSLOPE' )
  CALL ADD2DFIELD_ll( TZFIELDS_ll, PVSLOPE, 'UPDATE_ROTATE_WIND::PVSLOPE' )
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)
!!$ENDIF
!
!        3 Boundary conditions for non cyclic case
!
IF ( HLBCX(1) /= "CYCL" .AND. LWEST_ll()) THEN
  PUSLOPE(D%NIB-1,:)=PUSLOPE(D%NIB,:)
  PVSLOPE(D%NIB-1,:)=PVSLOPE(D%NIB,:)
END IF
IF ( HLBCX(2) /= "CYCL" .AND. LEAST_ll()) THEN
  PUSLOPE(D%NIE+1,:)=PUSLOPE(D%NIE,:)
  PVSLOPE(D%NIE+1,:)=PVSLOPE(D%NIE,:)
END IF
IF ( HLBCY(1) /= "CYCL" .AND. LSOUTH_ll()) THEN
  PUSLOPE(:,D%NJB-1)=PUSLOPE(:,D%NJB)
  PVSLOPE(:,D%NJB-1)=PVSLOPE(:,D%NJB)
END IF
IF(  HLBCY(2) /= "CYCL" .AND. LNORTH_ll()) THEN
  PUSLOPE(:,D%NJE+1)=PUSLOPE(:,D%NJE)
  PVSLOPE(:,D%NJE+1)=PVSLOPE(:,D%NJE)
END IF
!
IF (LHOOK) CALL DR_HOOK('TURB:UPDATE_ROTATE_WIND',1,ZHOOK_HANDLE2)
!
END SUBROUTINE UPDATE_ROTATE_WIND
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PT,PEXN,PCP
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PLOCPEXN
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PAMOIST,PATHETA
!
!-------------------------------------------------------------------------------
!
  IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO',0,ZHOOK_HANDLE2)
  ZEPS = CST%XMV / CST%XMD
!
!*       1.1 Lv/Cph at  t
!
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
  PLOCPEXN(IIB:IIE,IJB:IJE,1:D%NKT) = ( PLTT + (CST%XCPV-PC) *  (PT(IIB:IIE,IJB:IJE,1:D%NKT)-CST%XTT) ) &
                                     / PCP(IIB:IIE,IJB:IJE,1:D%NKT)
!
!*      1.2 Saturation vapor pressure at t
!
  ZRVSAT(IIB:IIE,IJB:IJE,1:D%NKT) =  EXP( PALP - PBETA/PT(IIB:IIE,IJB:IJE,1:D%NKT) - PGAM*ALOG( PT(IIB:IIE,IJB:IJE,1:D%NKT) ) )
!
!*      1.3 saturation  mixing ratio at t
!
  ZRVSAT(IIB:IIE,IJB:IJE,1:D%NKT) =  ZRVSAT(IIB:IIE,IJB:IJE,1:D%NKT) &
                                    * ZEPS / ( PPABST(IIB:IIE,IJB:IJE,1:D%NKT) - ZRVSAT(IIB:IIE,IJB:IJE,1:D%NKT) )
!
!*      1.4 compute the saturation mixing ratio derivative (rvs')
!
  ZDRVSATDT(IIB:IIE,IJB:IJE,1:D%NKT) = ( PBETA / PT(IIB:IIE,IJB:IJE,1:D%NKT)  - PGAM ) / PT(IIB:IIE,IJB:IJE,1:D%NKT)   &
                 * ZRVSAT(IIB:IIE,IJB:IJE,1:D%NKT) * ( 1. + ZRVSAT(IIB:IIE,IJB:IJE,1:D%NKT) / ZEPS )
!
!*      1.5 compute Amoist
!
  PAMOIST(IIB:IIE,IJB:IJE,1:D%NKT)=  0.5 / ( 1.0 + ZDRVSATDT(IIB:IIE,IJB:IJE,1:D%NKT) * PLOCPEXN(IIB:IIE,IJB:IJE,1:D%NKT) )
!
!*      1.6 compute Atheta
!
  PATHETA(IIB:IIE,IJB:IJE,1:D%NKT)= PAMOIST(IIB:IIE,IJB:IJE,1:D%NKT) * PEXN(IIB:IIE,IJB:IJE,1:D%NKT) *               &
        ( ( ZRVSAT(IIB:IIE,IJB:IJE,1:D%NKT) - PRT(IIB:IIE,IJB:IJE,1:D%NKT,1) ) * PLOCPEXN(IIB:IIE,IJB:IJE,1:D%NKT) / &
          ( 1. + ZDRVSATDT(IIB:IIE,IJB:IJE,1:D%NKT) * PLOCPEXN(IIB:IIE,IJB:IJE,1:D%NKT) )        *               &
          (                                                                  &
           ZRVSAT(IIB:IIE,IJB:IJE,1:D%NKT) * (1. + ZRVSAT(IIB:IIE,IJB:IJE,1:D%NKT)/ZEPS)                         &
                        * ( -2.*PBETA/PT(IIB:IIE,IJB:IJE,1:D%NKT) + PGAM ) / PT(IIB:IIE,IJB:IJE,1:D%NKT)**2      &
          +ZDRVSATDT(IIB:IIE,IJB:IJE,1:D%NKT) * (1. + 2. * ZRVSAT(IIB:IIE,IJB:IJE,1:D%NKT)/ZEPS)                 &
                        * ( PBETA/PT(IIB:IIE,IJB:IJE,1:D%NKT) - PGAM ) / PT(IIB:IIE,IJB:IJE,1:D%NKT)             &
          )                                                                  &
         - ZDRVSATDT(IIB:IIE,IJB:IJE,1:D%NKT)                                                  &
        )
!
!*      1.7 Lv/Cph/Exner at t-1
!
  PLOCPEXN(IIB:IIE,IJB:IJE,1:D%NKT) = PLOCPEXN(IIB:IIE,IJB:IJE,1:D%NKT) / PEXN(IIB:IIE,IJB:IJE,1:D%NKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
!
IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO',1,ZHOOK_HANDLE2)
END SUBROUTINE COMPUTE_FUNCTION_THERMO
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PLM
LOGICAL,                INTENT(IN)    :: ODZ
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TURB:DELT',0,ZHOOK_HANDLE2)
IF (ODZ) THEN
  ! Dz is take into account in the computation
  DO JK = IKTB,IKTE ! 1D turbulence scheme
    PLM(:,:,JK) = PZZ(:,:,JK+D%NKL) - PZZ(:,:,JK)
  END DO
  PLM(:,:,D%NKU) = PLM(:,:,IKE)
  PLM(:,:,D%NKA) = PZZ(:,:,IKB) - PZZ(:,:,D%NKA)
  IF ( HTURBDIM /= '1DIM' ) THEN  ! 3D turbulence scheme
    IF ( O2D) THEN
      PLM(:,:,:) = SQRT( PLM(:,:,:)*MXF(PDXX(:,:,:)) ) 
    ELSE
      PLM(:,:,:) = (PLM(:,:,:)*MXF(PDXX(:,:,:))*MYF(PDYY(:,:,:)) ) ** (1./3.)
    END IF
  END IF
ELSE
  ! Dz not taken into account in computation to assure invariability with vertical grid mesh
  PLM=1.E10
  IF ( HTURBDIM /= '1DIM' ) THEN  ! 3D turbulence scheme
    IF ( O2D) THEN
      PLM(:,:,:) = MXF(PDXX(:,:,:))
    ELSE
      PLM(:,:,:) = (MXF(PDXX(:,:,:))*MYF(PDYY(:,:,:)) ) ** (1./2.)
    END IF
  END IF
END IF
!
!  mixing length limited by the distance normal to the surface 
!  (with the same factor as for BL89)
!
IF (.NOT. ORMC01) THEN
  ZALPHA=0.5**(-1.5)
  !
  DO JJ=1,SIZE(PUT,2)
    DO JI=1,SIZE(PUT,1)
      IF (OOCEAN) THEN
        DO JK=IKTE,IKTB,-1
          ZD=ZALPHA*(PZZ(JI,JJ,IKTE+1)-PZZ(JI,JJ,JK))
          IF ( PLM(JI,JJ,JK)>ZD) THEN
            PLM(JI,JJ,JK)=ZD
          ELSE
            EXIT
          ENDIF
       END DO
      ELSE
        DO JK=IKTB,IKTE
          ZD=ZALPHA*(0.5*(PZZ(JI,JJ,JK)+PZZ(JI,JJ,JK+D%NKL))&
          -PZZ(JI,JJ,IKB)) *PDIRCOSZW(JI,JJ)
          IF ( PLM(JI,JJ,JK)>ZD) THEN
            PLM(JI,JJ,JK)=ZD
          ELSE
            EXIT
          ENDIF
        END DO
      ENDIF   
    END DO
  END DO
END IF
!
PLM(:,:,D%NKA) = PLM(:,:,IKB  )
PLM(:,:,D%NKU  ) = PLM(:,:,IKE)
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PLM
!
!-------------------------------------------------------------------------------
!
!   initialize the mixing length with the mesh grid
IF (LHOOK) CALL DR_HOOK('TURB:DEAR',0,ZHOOK_HANDLE2)
! 1D turbulence scheme
PLM(:,:,IKTB:IKTE) = PZZ(:,:,IKTB+D%NKL:IKTE+D%NKL) - PZZ(:,:,IKTB:IKTE)
PLM(:,:,D%NKU) = PLM(:,:,IKE)
PLM(:,:,D%NKA) = PZZ(:,:,IKB) - PZZ(:,:,D%NKA)
IF ( HTURBDIM /= '1DIM' ) THEN  ! 3D turbulence scheme
  IF ( O2D) THEN
    PLM(:,:,:) = SQRT( PLM(:,:,:)*MXF(PDXX(:,:,:)) )
  ELSE
    PLM(:,:,:) = (PLM(:,:,:)*MXF(PDXX(:,:,:))*MYF(PDYY(:,:,:)) ) ** (1./3.)
  END IF
END IF
!   compute a mixing length limited by the stability
!
CALL ETHETA(D,CST,KRR,KRRI,PTHLT,PRT,ZLOCPEXNM,ZATHETA,PSRCT,OOCEAN,OCOMPUTE_SRC,ZETHETA)
CALL EMOIST(D,CST,KRR,KRRI,PTHLT,PRT,ZLOCPEXNM,ZAMOIST,PSRCT,OOCEAN,ZEMOIST)
!
IF (KRR>0) THEN
  DO JK = IKTB+1,IKTE-1
    DO JJ=1,SIZE(PUT,2)
      DO JI=1,SIZE(PUT,1)
        ZDTHLDZ(JI,JJ,JK)= 0.5*((PTHLT(JI,JJ,JK+D%NKL)-PTHLT(JI,JJ,JK    ))/PDZZ(JI,JJ,JK+D%NKL)+ &
                                (PTHLT(JI,JJ,JK    )-PTHLT(JI,JJ,JK-D%NKL))/PDZZ(JI,JJ,JK    ))
        ZDRTDZ(JI,JJ,JK) = 0.5*((PRT(JI,JJ,JK+D%NKL,1)-PRT(JI,JJ,JK    ,1))/PDZZ(JI,JJ,JK+D%NKL)+ &
                                (PRT(JI,JJ,JK    ,1)-PRT(JI,JJ,JK-D%NKL,1))/PDZZ(JI,JJ,JK    ))
        IF (OOCEAN) THEN
          ZVAR=CST%XG*(CST%XALPHAOC*ZDTHLDZ(JI,JJ,JK)-CST%XBETAOC*ZDRTDZ(JI,JJ,JK))
        ELSE
          ZVAR=CST%XG/PTHVREF(JI,JJ,JK)*                                                  &
             (ZETHETA(JI,JJ,JK)*ZDTHLDZ(JI,JJ,JK)+ZEMOIST(JI,JJ,JK)*ZDRTDZ(JI,JJ,JK))
        END IF
        !
        IF (ZVAR>0.) THEN
          PLM(JI,JJ,JK)=MAX(CST%XMNH_EPSILON,MIN(PLM(JI,JJ,JK), &
                        0.76* SQRT(PTKET(JI,JJ,JK)/ZVAR)))
        END IF
      END DO
    END DO
  END DO
ELSE! For dry atmos or unsalted ocean runs
  DO JK = IKTB+1,IKTE-1
    DO JJ=1,SIZE(PUT,2)
      DO JI=1,SIZE(PUT,1)
        ZDTHLDZ(JI,JJ,JK)= 0.5*((PTHLT(JI,JJ,JK+D%NKL)-PTHLT(JI,JJ,JK    ))/PDZZ(JI,JJ,JK+D%NKL)+ &
                                (PTHLT(JI,JJ,JK    )-PTHLT(JI,JJ,JK-D%NKL))/PDZZ(JI,JJ,JK    ))
        IF (OOCEAN) THEN
          ZVAR= CST%XG*CST%XALPHAOC*ZDTHLDZ(JI,JJ,JK)
        ELSE
          ZVAR= CST%XG/PTHVREF(JI,JJ,JK)*ZETHETA(JI,JJ,JK)*ZDTHLDZ(JI,JJ,JK)
        END IF
!
        IF (ZVAR>0.) THEN
          PLM(JI,JJ,JK)=MAX(CST%XMNH_EPSILON,MIN(PLM(JI,JJ,JK), &
                        0.76* SQRT(PTKET(JI,JJ,JK)/ZVAR)))
        END IF
      END DO
    END DO
  END DO
END IF
!  special case near the surface 
ZDTHLDZ(:,:,IKB)=(PTHLT(:,:,IKB+D%NKL)-PTHLT(:,:,IKB))/PDZZ(:,:,IKB+D%NKL)
! For dry simulations
IF (KRR>0) THEN
  ZDRTDZ(:,:,IKB)=(PRT(:,:,IKB+D%NKL,1)-PRT(:,:,IKB,1))/PDZZ(:,:,IKB+D%NKL)
ELSE
  ZDRTDZ(:,:,IKB)=0
ENDIF
!
IF (OOCEAN) THEN
  ZWORK2D(:,:)=CST%XG*(CST%XALPHAOC*ZDTHLDZ(:,:,IKB)-CST%XBETAOC*ZDRTDZ(:,:,IKB))
ELSE
  ZWORK2D(:,:)=CST%XG/PTHVREF(:,:,IKB)*                                           &
              (ZETHETA(:,:,IKB)*ZDTHLDZ(:,:,IKB)+ZEMOIST(:,:,IKB)*ZDRTDZ(:,:,IKB))
END IF
WHERE(ZWORK2D(:,:)>0.)
  PLM(:,:,IKB)=MAX(CST%XMNH_EPSILON,MIN( PLM(:,:,IKB),                 &
                    0.76* SQRT(PTKET(:,:,IKB)/ZWORK2D(:,:))))
END WHERE
!
!  mixing length limited by the distance normal to the surface (with the same factor as for BL89)
!
IF (.NOT. ORMC01) THEN
  ZALPHA=0.5**(-1.5)
  !
  DO JJ=1,SIZE(PUT,2)
    DO JI=1,SIZE(PUT,1)
      IF (OOCEAN) THEN
        DO JK=IKTE,IKTB,-1
          ZD=ZALPHA*(PZZ(JI,JJ,IKTE+1)-PZZ(JI,JJ,JK))
          IF ( PLM(JI,JJ,JK)>ZD) THEN
            PLM(JI,JJ,JK)=ZD
          ELSE
            EXIT
          ENDIF
        END DO
      ELSE
        DO JK=IKTB,IKTE
          ZD=ZALPHA*(0.5*(PZZ(JI,JJ,JK)+PZZ(JI,JJ,JK+D%NKL))-PZZ(JI,JJ,IKB)) &
            *PDIRCOSZW(JI,JJ)
          IF ( PLM(JI,JJ,JK)>ZD) THEN
            PLM(JI,JJ,JK)=ZD
          ELSE
            EXIT
          ENDIF
        END DO
      ENDIF 
    END DO
  END DO
END IF
!
PLM(:,:,D%NKA) = PLM(:,:,IKB  )
PLM(:,:,IKE  ) = PLM(:,:,IKE-D%NKL)
PLM(:,:,D%NKU  ) = PLM(:,:,D%NKU-D%NKL)
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
ZCOEF_AMPL(:,:,:) = 1.
!
!*       2.    CALCULATION OF THE AMPLIFICATION COEFFICIENT
!              --------------------------------------------
!
! Saturation
!
WHERE ( PCEI(:,:,:)>=PCEI_MAX ) ZCOEF_AMPL(:,:,:)=PCOEF_AMPL_SAT
!
! Between the min and max limits of CEI index, linear variation of the
! amplification coefficient ZCOEF_AMPL as a function of CEI
!
WHERE ( PCEI(:,:,:) <  PCEI_MAX .AND.                                        &
        PCEI(:,:,:) >  PCEI_MIN      )                                       &
        ZCOEF_AMPL(:,:,:) = ZPENTE * PCEI(:,:,:) + ZCOEF_AMPL_CEI_NUL  
!
!
!*       3.    CALCULATION OF THE MIXING LENGTH IN CLOUDS
!              ------------------------------------------
!
IF (HTURBLEN_CL == HTURBLEN) THEN
  ZLM_CLOUD(:,:,:) = ZLM(:,:,:)
ELSE
  SELECT CASE (HTURBLEN_CL)
!
!*         3.1 BL89 mixing length
!           ------------------
  CASE ('BL89','RM17','ADAP')
    ZSHEAR=0.
    CALL BL89(D,CST,CSTURB,PZZ,PDZZ,PTHVREF,ZTHLM,KRR,ZRM,PTKET,ZSHEAR,ZLM_CLOUD,OOCEAN,HPROGRAM)
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
IF ( OTURB_DIAG .AND. TPFILE%LOPENED ) THEN
  TZFIELD%CMNHNAME   = 'LM_CLEAR_SKY'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'LM_CLEAR_SKY'
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_LM CLEAR SKY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZLM)
ENDIF
!
! Amplification of the mixing length when the criteria are verified
!
WHERE (ZCOEF_AMPL(:,:,:) /= 1.) ZLM(:,:,:) = ZCOEF_AMPL(:,:,:)*ZLM_CLOUD(:,:,:)
!
! Cloud mixing length in the clouds at the points which do not verified the CEI
!
WHERE (PCEI(:,:,:) == -1.) ZLM(:,:,:) = ZLM_CLOUD(:,:,:)
!
!
!*       5.    IMPRESSION
!              ----------
!
IF ( OTURB_DIAG .AND. TPFILE%LOPENED ) THEN
  TZFIELD%CMNHNAME   = 'COEF_AMPL'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'COEF_AMPL'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_COEF AMPL'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZCOEF_AMPL)
  !
  TZFIELD%CMNHNAME   = 'LM_CLOUD'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'LM_CLOUD'
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_LM CLOUD'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZLM_CLOUD)
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TURB:CLOUD_MODIF_LM',1,ZHOOK_HANDLE2)
END SUBROUTINE CLOUD_MODIF_LM
!
END SUBROUTINE TURB    
