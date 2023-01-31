!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_VER
IMPLICIT NONE
CONTAINS
SUBROUTINE TURB_VER(D,CST,CSTURB,TURBN,TLES,KRR,KRRL,KRRI,KGRADIENTS,&
                      OOCEAN,ODEEPOC,OCOMPUTE_SRC,                  &
                      KSV,KSV_LGBEG,KSV_LGEND,                      &
                      PEXPL, HPROGRAM, O2D, ONOMIXLG, OFLAT,        &
                      OCOUPLES,OBLOWSNOW,OFLYER,PRSNOW,             & 
                      PTSTEP, TPFILE,                               &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PCOSSLOPE,PSINSLOPE,                          &
                      PRHODJ,PTHVREF,PSFUM,PSFVM,                   &
                      PSFTHM,PSFRM,PSFSVM,PSFTHP,PSFRP,PSFSVP,      &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU33M,              &
                      PUM,PVM,PWM,PUSLOPEM,PVSLOPEM,PTHLM,PRM,PSVM, &
                      PTKEM,PLM,PLENGTHM,PLENGTHH,PLEPS,MFMOIST,    &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,    &
                      PFWTH,PFWR,PFTH2,PFR2,PFTHR,PBL_DEPTH,        &
                      PSBL_DEPTH,PLMO,PHGRAD,PZS,                   &
                      PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,             &
                      PDP,PTP,PSIGS,PWTH,PWRC,PWSV,                 &
                      PSSTFL,PSSTFL_C,PSSRFL_C,PSSUFL_C,PSSVFL_C,   &
                      PSSUFL,PSSVFL                                 )
!     ###############################################################
!
!
!!****  *TURB_VER* -compute the source terms due to the vertical turbulent
!!       fluxes.
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute the vertical turbulent
!     fluxes of the evolutive variables and give back the source 
!     terms to the main program.	In the case of large horizontal meshes,
!     the divergence of these vertical turbulent fluxes represent the whole
!     effect of the turbulence but when the three-dimensionnal version of
!     the turbulence scheme is activated (CTURBDIM="3DIM"), these divergences
!     are completed in the next routine TURB_HOR. 
!		  An arbitrary degree of implicitness has been implemented for the 
!     temporal treatment of these diffusion terms.
!       The vertical boundary conditions are as follows:
!           *  at the bottom, the surface fluxes are prescribed at the same
!              as the other turbulent fluxes
!           *  at the top, the turbulent fluxes are set to 0.
!       It should be noted that the condensation has been implicitely included
!     in this turbulence scheme by using conservative variables and computing
!     the subgrid variance of a statistical variable s indicating the presence 
!     or not of condensation in a given mesh. 
!
!!**  METHOD
!!    ------
!!      1D type calculations are made;
!!      The vertical turbulent fluxes are computed in an off-centered
!!      implicit scheme (a Crank-Nicholson type with coefficients different
!!      than 0.5), which allows to vary the degree of implicitness of the
!!      formulation.
!!      	 The different prognostic variables are treated one by one. 
!!      The contributions of each turbulent fluxes are cumulated into the 
!!      tendency  PRvarS, and into the dynamic and thermal production of 
!!      TKE if necessary.
!!        
!!			 In section 2 and 3, the thermodynamical fields are considered.
!!      Only the turbulent fluxes of the conservative variables
!!      (Thetal and Rnp stored in PRx(:,:,:,1))  are computed. 
!!       Note that the turbulent fluxes at the vertical 
!!      boundaries are given either by the soil scheme for the surface one
!!      ( at the same instant as the others fluxes) and equal to 0 at the 
!!      top of the model. The thermal production is computed by vertically 
!!      averaging the turbulent flux and multiply this flux at the mass point by
!!      a function ETHETA or EMOIST, which preform the transformation from the
!!      conservative variables to the virtual potential temperature. 
!!     
!! 	    In section 4, the variance of the statistical variable
!!      s indicating presence or not of condensation, is determined in function 
!!      of the turbulent moments of the conservative variables and its
!!      squarred root is stored in PSIGS. This information will be completed in 
!!      the horizontal turbulence if the turbulence dimensionality is not 
!!      equal to "1DIM".
!!
!!			 In section 5, the x component of the stress tensor is computed.
!!      The surface flux <u'w'> is computed from the value of the surface
!!      fluxes computed in axes linked to the orography ( i", j" , k"):
!!        i" is parallel to the surface and in the direction of the maximum
!!           slope
!!        j" is also parallel to the surface and in the normal direction of
!!           the maximum slope
!!        k" is the normal to the surface
!!      In order to prevent numerical instability, the implicit scheme has 
!!      been extended to the surface flux regarding to its dependence in 
!!      function of U. The dependence in function of the other components 
!!      introduced by the different rotations is only explicit.
!!      The turbulent fluxes are used to compute the dynamic production of 
!!      TKE. For the last TKE level ( located at PDZZ(:,:,IKB)/2 from the
!!      ground), an harmonic extrapolation from the dynamic production at 
!!      PDZZ(:,:,IKB) is used to avoid an evaluation of the gradient of U
!!      in the surface layer.
!!
!!         In section 6, the same steps are repeated but for the y direction
!!	and in section 7, a diagnostic computation of the W variance is 
!!      performed.
!!
!!         In section 8, the turbulent fluxes for the scalar variables are 
!!      computed by the same way as the conservative thermodynamical variables
!!
!!            
!!    EXTERNAL
!!    --------
!!      GX_U_M, GY_V_M, GZ_W_M :  cartesian gradient operators 
!!      GX_U_UW,GY_V_VW	         (X,Y,Z) represent the direction of the gradient
!!                               _(M,U,...)_ represent the localization of the 
!!                               field to be derivated
!!                               _(M,UW,...) represent the localization of the 
!!                               field	derivated
!!
!!      SUBROUTINE TRIDIAG     : to compute the split implicit evolution
!!                               of a variable located at a mass point
!!
!!      SUBROUTINE TRIDIAG_WIND: to compute the split implicit evolution
!!                               of a variable located at a wind point
!!
!!      FUNCTIONs ETHETA and EMOIST  :  
!!            allows to compute:
!!            - the coefficients for the turbulent correlation between
!!            any variable and the virtual potential temperature, of its 
!!            correlations with the conservative potential temperature and 
!!            the humidity conservative variable:
!!            -------              -------              -------
!!            A' Thv'  =  ETHETA   A' Thl'  +  EMOIST   A' Rnp'  
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST : contains physical constants
!!
!!           CST%XG         : gravity constant
!!
!!      Module MODD_CTURB: contains the set of constants for
!!                        the turbulence scheme
!!
!!           CSTURB%XCMFS,XCMFB : cts for the momentum flux
!!           CSTURB%XCSHF       : ct for the sensible heat flux
!!           CSTURB%XCHF        : ct for the moisture flux
!!           CSTURB%XCTV,CSTURB%XCHV   : cts for the T and moisture variances
!!
!!      Module MODD_PARAMETERS
!!
!!           JPVEXT_TURB     : number of vertical external points
!!
!!
!!    REFERENCE
!!    ---------
!!      Book 1 of documentation (Chapter: Turbulence)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart             * INM and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       August   19, 1994
!!      Modifications: February 14, 1995 (J.Cuxart and J.Stein) 
!!                                  Doctorization and Optimization
!!      Modifications: March 21, 1995 (J.M. Carriere) 
!!                                  Introduction of cloud water
!!      Modifications: June  14, 1995 (J.Cuxart and J. Stein) 
!!                                 Phi3 and Psi3 at w-point + bug in the all
!!                                 or nothing condens. 
!!      Modifications: Sept  15, 1995 (J.Cuxart and J. Stein) 
!!                                 Change the DP computation at the ground
!!      Modifications: October 10, 1995 (J.Cuxart and J. Stein) 
!!                                 Psi for scal var and LES tools
!!      Modifications: November 10, 1995 (J. Stein)
!!                                 change the surface	relations 
!!      Modifications: February 20, 1995 (J. Stein) optimization
!!      Modifications: May 21, 1996 (J. Stein) 
!!                                  bug in the vertical flux of the V wind 
!!                                  component for explicit computation
!!      Modifications: May 21, 1996 (N. wood) 
!!                                  modify the computation of the vertical
!!                                   part or the surface tangential flux
!!      Modifications: May 21, 1996 (P. Jabouille)
!!                                  same modification in the Y direction
!!      
!!      Modifications: Sept 17, 1996 (J. Stein) change the moist case by using
!!                                  Pi instead of Piref + use Atheta and Amoist
!!
!!      Modifications: Nov  24, 1997 (V. Masson) removes the DO loops 
!!      Modifications: Mar  31, 1998 (V. Masson) splits the routine TURB_VER 
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     Feb  20, 2003 (JP Pinty)  Add PFRAC_ICE
!!                     July     2005 (S. Tomas, V. Masson)
!!                                               Add 3rd order moments and
!!                                               implicitation of PHI3, PSI3
!!                     Oct.2009  (C.Lac) Introduction of different PTSTEP according to the
!!                              advection schemes
!!                     Feb. 2012  (Y. Seity) add possibility to run with
!!                                 reversed vertical levels
!!                     10/2012 (J.Escobar) Bypass PGI bug , redefine some allocatable array inplace of automatic
!!                     08/2014 (J.Escobar) Bypass PGI memory leak bug , replace IF statement with IF THEN ENDIF
!!      Modifications: July,    2015  (Wim de Rooy) switch for HARATU (Racmo turbulence scheme)
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! JL Redelsperger 03/2021 : add Ocean LES case
!!--------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY: JPRB
USE YOMHOOK,  ONLY: LHOOK, DR_HOOK
!
USE MODD_CST,            ONLY: CST_t
USE MODD_CTURB,          ONLY: CSTURB_t
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_FIELD,          ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS,     ONLY: JPVEXT_TURB
USE MODD_LES,            ONLY: TLES_t
USE MODD_TURB_n,         ONLY: TURB_t
!
USE MODE_EMOIST,               ONLY: EMOIST
USE MODE_ETHETA,               ONLY: ETHETA
USE MODE_GRADIENT_M_PHY,       ONLY: GZ_M_W_PHY
USE MODE_IO_FIELD_WRITE,       ONLY: IO_FIELD_WRITE_PHY
USE MODE_PRANDTL,              ONLY: PSI_SV, PSI3, PHI3, PRANDTL
USE MODE_SBL_DEPTH,            ONLY: SBL_DEPTH
USE MODE_TURB_VER_DYN_FLUX,    ONLY: TURB_VER_DYN_FLUX
USE MODE_TURB_VER_SV_FLUX,     ONLY: TURB_VER_SV_FLUX
USE MODE_TURB_VER_SV_CORR,     ONLY: TURB_VER_SV_CORR
USE MODE_TURB_VER_THERMO_FLUX, ONLY: TURB_VER_THERMO_FLUX
USE MODE_TURB_VER_THERMO_CORR, ONLY: TURB_VER_THERMO_CORR
!
USE MODI_LES_MEAN_SUBGRID_PHY
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
TYPE(TURB_t),           INTENT(IN)   :: TURBN
TYPE(TLES_t),           INTENT(INOUT):: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KGRADIENTS    ! Number of stored horizontal gradients
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV, KSV_LGBEG, KSV_LGEND ! number of scalar variables
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  ODEEPOC      ! activates sfc forcing for ideal ocean deep conv
LOGICAL,                INTENT(IN)   ::  OFLYER       ! MesoNH flyer diagnostic
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES version 
LOGICAL,                INTENT(IN)   ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
REAL,                   INTENT(IN)   ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
CHARACTER(LEN=6),       INTENT(IN)   ::  HPROGRAM     ! HPROGRAM is the program currently running
LOGICAL,                INTENT(IN)   ::  ONOMIXLG     ! to use turbulence for lagrangian variables
LOGICAL,                INTENT(IN)   ::  O2D          ! Logical for 2D model version
REAL,                   INTENT(IN)   ::  PEXPL        ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! timestep 
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                      ! Metric coefficients
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PZZ          ! altitudes at flux points
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCOSSLOPE    ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSINSLOPE    ! sinus of the angle 
                                      ! between i and the slope vector
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PZS ! orography (for LEONARD terms)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
! MFMOIST used in case of TURBN%LHARATU
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  MFMOIST       ! moist mass flux dual scheme

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTHVREF      ! ref. state Virtual 
                                                      ! Potential Temperature 
REAL, DIMENSION(D%NIJT,D%NKT,KGRADIENTS),INTENT(IN) :: PHGRAD  ! horizontal gradients
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSFUM,PSFVM ! surface fluxes
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSFTHM,PSFRM ! surface fluxes at time
REAL, DIMENSION(D%NIJT,KSV), INTENT(IN)   ::  PSFSVM       ! t - deltat 
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSFTHP,PSFRP ! surface fluxes at time
REAL, DIMENSION(D%NIJT,KSV), INTENT(IN)   ::  PSFSVP       ! t + deltat 
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCDUEFF     ! Cd * || u || at time t
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PUM,PVM,PWM,PTHLM 
  ! Wind and potential temperature at t-Delta t
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN) ::  PRM          ! Mixing ratios 
                                                      ! at t-Delta t
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLM          ! Turb. mixing length   
! PLENGTHM PLENGTHH used in case of TURBN%LHARATU
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLENGTHM     ! Turb. mixing length momentum
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLENGTHH     ! Turb. mixing length heat/moisture 
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLEPS        ! dissipative length
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(IN)   ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFRAC_ICE    ! ri fraction of rc+ri
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFWTH        ! d(w'2th' )/dz
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFWR         ! d(w'2r'  )/dz
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFTH2        ! d(w'th'2 )/dz
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFR2         ! d(w'r'2  )/dz
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFTHR        ! d(w'th'r')/dz
REAL, DIMENSION(MERGE(D%NIJT,0,TURBN%CTOM=='TM06')),INTENT(INOUT) :: PBL_DEPTH  ! BL height for TOMS
REAL, DIMENSION(MERGE(D%NIJT,0,TURBN%LRMC01)),INTENT(INOUT) :: PSBL_DEPTH ! SBL depth for RMC01
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PLMO         ! Monin-Obukhov length
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT)   ::  PRUS, PRVS, PRWS, PRTHLS
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT) ::  PRRS
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(INOUT) ::  PRSVS
                            ! cumulated sources for the prognostic variables
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PDP,PTP   ! Dynamic and thermal
                                                   ! TKE production terms
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(OUT)     ::  PSIGS
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PWTH      ! heat flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PWRC      ! cloud water flux
REAL, DIMENSION(D%NIJT,D%NKT,KSV),INTENT(OUT) :: PWSV       ! scalar flux
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSTFL        ! Time evol Flux of T at sea surface (LOCEAN)!
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSTFL_C  ! O-A interface flux for theta(LOCEAN and LCOUPLES)
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSRFL_C  ! O-A interface flux for vapor (LOCEAN and LCOUPLES) 
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSUFL_C        ! Time evol Flux of U at sea surface (LOCEAN)
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSVFL_C  !
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSUFL   
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSVFL  !
!
!*       0.2  declaration of local variables
!
REAL,  DIMENSION(D%NIJT,D%NKT)  ::  &
       ZBETA,    & ! buoyancy coefficient
       ZSQRT_TKE,& ! sqrt(e)
       ZDTH_DZ,  & ! d(th)/dz
       ZDR_DZ,   & ! d(rt)/dz
       ZRED2TH3, & ! 3D Redeslperger number R*2_th
       ZRED2R3,  & ! 3D Redeslperger number R*2_r
       ZRED2THR3,& ! 3D Redeslperger number R*2_thr
       ZBLL_O_E, & ! beta * Lk * Leps / tke
       ZETHETA,  & ! Coefficient for theta in theta_v computation
       ZEMOIST,  & ! Coefficient for r in theta_v computation
       ZREDTH1,  & ! 1D Redelsperger number for Th
       ZREDR1,   & ! 1D Redelsperger number for r
       ZPHI3,    & ! phi3 Prandtl number
       ZPSI3,    & ! psi3 Prandtl number for vapor
       ZD,       & ! denominator in phi3 terms
       ZWTHV,    & ! buoyancy flux
       ZWU,      & ! (u'w')
       ZWV,      & ! (v'w')
       ZTHLP,    & ! guess of potential temperature due to vert. turbulent flux
       ZRP         ! guess of total water due to vert. turbulent flux
!
REAL, DIMENSION(D%NIJT,D%NKT,KSV)  ::  &
       ZPSI_SV,  & ! Prandtl number for scalars
       ZREDS1,   & ! 1D Redelsperger number R_sv
       ZRED2THS, & ! 3D Redelsperger number R*2_thsv
       ZRED2RS     ! 3D Redelsperger number R*2_rsv
REAL, DIMENSION(D%NIJT,D%NKT)  ::  ZLM
!
LOGICAL :: GUSERV    ! flag to use water vapor
INTEGER :: IKB,IKE,IIJE,IIJB,IKT   ! index value for the Beginning
                                     ! and the End of the physical domain for the mass points
INTEGER :: JSV,JIJ,JK ! loop counter
REAL    :: ZTIME1
REAL    :: ZTIME2
REAL(KIND=JPRB)      :: ZHOOK_HANDLE
TYPE(TFIELDMETADATA) :: TZFIELD
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
!
IF (LHOOK) CALL DR_HOOK('TURB_VER',0,ZHOOK_HANDLE)
!
IKB=D%NKTB
IKE=D%NKTE
IKT=D%NKT
IIJE=D%NIJE
IIJB=D%NIJB
!
!
! 3D Redelsperger numbers
!
!
CALL PRANDTL(D,CST,CSTURB,KRR,KSV,KRRI,TURBN%LTURB_FLX,  &
             TURBN%CTURBDIM,OOCEAN,TURBN%LHARAT,O2D,OCOMPUTE_SRC,&
             TPFILE, OFLAT,                        &
             PDXX,PDYY,PDZZ,PDZX,PDZY,             &
             PTHVREF,PLOCPEXNM,PATHETA,PAMOIST,    &
             PLM,PLEPS,PTKEM,PTHLM,PRM,PSVM,PSRCM, &
             ZREDTH1, ZREDR1,                      &
             ZRED2TH3, ZRED2R3, ZRED2THR3,         &
             ZREDS1,ZRED2THS, ZRED2RS,             &
             ZBLL_O_E,                             &
             ZETHETA, ZEMOIST                      )
!
! Buoyancy coefficient
!
IF (OOCEAN) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZBETA(IIJB:IIJE,1:IKT) = CST%XG*CST%XALPHAOC
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZBETA(IIJB:IIJE,1:IKT) = CST%XG/PTHVREF(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
! Square root of Tke
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZSQRT_TKE(IIJB:IIJE,1:IKT) = SQRT(PTKEM(IIJB:IIJE,1:IKT))
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
! gradients of mean quantities at previous time-step
!
CALL GZ_M_W_PHY(D,PTHLM,PDZZ,ZDTH_DZ)
ZDR_DZ(:,:)  = 0.
IF (KRR>0) CALL GZ_M_W_PHY(D,PRM(:,:,1),PDZZ,ZDR_DZ)
!
!
! Denominator factor in 3rd order terms
!
IF (.NOT. TURBN%LHARAT) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZD(IIJB:IIJE,1:IKT) = (1.+ZREDTH1(IIJB:IIJE,1:IKT)+ZREDR1(IIJB:IIJE,1:IKT)) * &
                               &(1.+0.5*(ZREDTH1(IIJB:IIJE,1:IKT)+ZREDR1(IIJB:IIJE,1:IKT)))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZD(IIJB:IIJE,1:IKT) = 1.
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ENDIF
!
! Phi3 and Psi3 Prandtl numbers
!
GUSERV = KRR/=0
!
CALL PHI3(D,CSTURB,ZREDTH1,ZREDR1,ZRED2TH3,ZRED2R3,ZRED2THR3,TURBN%CTURBDIM,GUSERV,ZPHI3)
IF(KRR/=0) &
CALL PSI3(D,CSTURB,ZREDR1,ZREDTH1,ZRED2R3,ZRED2TH3,ZRED2THR3,TURBN%CTURBDIM,GUSERV,ZPSI3)
!
! Prandtl numbers for scalars
!
CALL PSI_SV(D,CSTURB,KSV,ZREDTH1,ZREDR1,ZREDS1,ZRED2THS,ZRED2RS,ZPHI3,ZPSI3,ZPSI_SV)
!
! LES diagnostics
!
IF (TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZPHI3,TLES%X_LES_SUBGRID_PHI3)
  IF(KRR/=0) THEN
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZPSI3,TLES%X_LES_SUBGRID_PSI3)
  END IF
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!----------------------------------------------------------------------------
!
!
!*       2.   SOURCES OF CONSERVATIVE POTENTIAL TEMPERATURE AND 
!                                                  PARTIAL THERMAL PRODUCTION 
!             ---------------------------------------------------------------
!
!*       3.   SOURCES OF CONSERVATIVE AND CLOUD MIXING RATIO AND 
!                                        COMPLETE THERMAL PRODUCTION 
!             ------------------------------------------------------
!
!*       4.   TURBULENT CORRELATIONS : <w Rc>, <THl THl>, <THl Rnp>, <Rnp Rnp>
!             ----------------------------------------------------------------
!

IF (TURBN%LHARAT) THEN
  ZLM(:,:)=PLENGTHH(:,:)
ELSE
  ZLM(:,:)=PLM(:,:)
ENDIF
!
  CALL  TURB_VER_THERMO_FLUX(D,CST,CSTURB,TURBN,TLES,                 &
                        KRR,KRRL,KRRI,KSV,KGRADIENTS,                 &
                        OOCEAN,ODEEPOC,OFLYER,                        &
                        OCOUPLES,OCOMPUTE_SRC,                        &
                        PEXPL,PTSTEP,HPROGRAM,TPFILE,                 &
                        PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                        PRHODJ,PTHVREF,PHGRAD,PZS,                    &
                        PSFTHM,PSFRM,PSFTHP,PSFRP,                    &
                        PWM,PTHLM,PRM,PSVM,                           &
                        PTKEM,ZLM,PLEPS,                              &
                        PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,    &
                        ZBETA, ZSQRT_TKE, ZDTH_DZ, ZDR_DZ, ZRED2TH3,  &
                        ZRED2R3, ZRED2THR3, ZBLL_O_E, ZETHETA,        &
                        ZEMOIST, ZREDTH1, ZREDR1, ZPHI3, ZPSI3, ZD,   &
                        PFWTH,PFWR,PFTH2,PFR2,PFTHR,                  &
                        MFMOIST,PBL_DEPTH,ZWTHV,                      &
                        PRTHLS,PRRS,ZTHLP,ZRP,PTP,PWTH,PWRC,          &
                        PSSTFL, PSSTFL_C, PSSRFL_C                    )
!
  CALL  TURB_VER_THERMO_CORR(D,CST,CSTURB,TURBN,TLES,                 &
                        KRR,KRRL,KRRI,KSV,                            &
                        OCOMPUTE_SRC,                                 &
                        OCOUPLES,                                     &
                        PEXPL,TPFILE,                                 &
                        PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,           &
                        PRHODJ,PTHVREF,                               &
                        PSFTHM,PSFRM,PSFTHP,PSFRP,                    &
                        PWM,PTHLM,PRM,PSVM,                           &
                        PTKEM,ZLM,PLEPS,                              &
                        PLOCPEXNM,PATHETA,PAMOIST,                    &
                        ZBETA, ZSQRT_TKE, ZDTH_DZ, ZDR_DZ, ZRED2TH3,  &
                        ZRED2R3, ZRED2THR3, ZBLL_O_E, ZETHETA,        &
                        ZEMOIST, ZREDTH1, ZREDR1, ZPHI3, ZPSI3, ZD,   &
                        PFWTH,PFWR,PFTH2,PFR2,PFTHR,                  &
                        ZTHLP,ZRP,MFMOIST,PSIGS                  )
!
!----------------------------------------------------------------------------
!
!
!
!*       5.   SOURCES OF U,W WIND COMPONENTS AND PARTIAL DYNAMIC PRODUCTION 
!             -------------------------------------------------------------
!
!*       6.   SOURCES OF V,W WIND COMPONENTS AND COMPLETE 1D DYNAMIC PRODUCTION 
!             -----------------------------------------------------------------
!
!*       7.   DIAGNOSTIC COMPUTATION OF THE 1D <W W> VARIANCE
!             -----------------------------------------------
!
!
IF (TURBN%LHARAT) ZLM(:,:)=PLENGTHM(:,:)
!
CALL  TURB_VER_DYN_FLUX(D,CST,CSTURB,TURBN,TLES,KSV,O2D,OFLAT,      &
                      KRR,OOCEAN,OCOUPLES,                          &
                      PEXPL,PTSTEP,TPFILE,                          &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PCOSSLOPE,PSINSLOPE,                          &
                      PRHODJ,PSFUM,PSFVM,                           &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU33M,              &
                      PTHLM,PRM,PSVM,PUM,PVM,PWM,PUSLOPEM,PVSLOPEM, &
                      PTKEM,ZLM,MFMOIST,ZWU,ZWV,                    &
                      PRUS,PRVS,PRWS,                               &
                      PDP,PTP,PSSUFL_C,PSSVFL_C,PSSUFL,PSSVFL       )
!
!----------------------------------------------------------------------------
!
!
!*       8.   SOURCES OF PASSIVE SCALAR VARIABLES
!             -----------------------------------
!
IF (TURBN%LHARAT) ZLM(:,:)=PLENGTHH(:,:)
!
IF (KSV>0)                                                          &
CALL  TURB_VER_SV_FLUX(D,CST,CSTURB,TURBN,TLES,ONOMIXLG,            &
                      KSV,KSV_LGBEG,KSV_LGEND,                      &
                      OBLOWSNOW,OFLYER,                             &
                      PEXPL,PTSTEP,TPFILE,PRSNOW,                   &
                      PDZZ,PDIRCOSZW,                               &
                      PRHODJ,PWM,                                   &
                      PSFSVM,PSFSVP,                                &
                      PSVM,                                         &
                      PTKEM,ZLM,MFMOIST,ZPSI_SV,                    &
                      PRSVS,PWSV                                    )
!
!
IF (KSV>0 .AND. TLES%LLES_CALL)                                          &
CALL  TURB_VER_SV_CORR(D,CST,CSTURB,TLES,KRR,KRRL,KRRI,OOCEAN,      &
                      PDZZ,KSV,KSV_LGBEG,KSV_LGEND,ONOMIXLG,        &
                      OBLOWSNOW,OCOMPUTE_SRC,PRSNOW,                &
                      PTHLM,PRM,PTHVREF,                            &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,ZPHI3,ZPSI3,  &
                      PWM,PSVM,                                     &
                      PTKEM,ZLM,PLEPS,ZPSI_SV                       )
!
!
!----------------------------------------------------------------------------
!
!*       9.   DIAGNOSTIC OF Surface Boundary Layer Depth
!             ------------------------------------------
!
IF (TURBN%LRMC01) CALL SBL_DEPTH(D,CSTURB,PZZ,ZWU,ZWV,ZWTHV,PLMO,PSBL_DEPTH)
!
!----------------------------------------------------------------------------
!
!
!*      10.   PRINTS
!             ------
!
!
IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED .AND. .NOT. TURBN%LHARAT) THEN
!
! stores the Turbulent Prandtl number
! 
  TZFIELD = TFIELDMETADATA(                  &
    CMNHNAME   = 'PHI3',                     &
    CSTDNAME   = '',                         &
    CLONGNAME  = 'PHI3',                     &
    CUNITS     = '1',                        &
    CDIR       = 'XY',                       &
    CCOMMENT   = 'Turbulent Prandtl number', &
    NGRID      = 4,                          &
    NTYPE      = TYPEREAL,                   &
    NDIMS      = 3,                          &
    LTIMEDEP   = .TRUE.                      )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZPHI3)
!
! stores the Turbulent Schmidt number
! 
  TZFIELD = TFIELDMETADATA(                  &
    CMNHNAME   = 'PSI3',                     &
    CSTDNAME   = '',                         &
    CLONGNAME  = 'PSI3',                     &
    CUNITS     = '1',                        &
    CDIR       = 'XY',                       &
    CCOMMENT   = 'Turbulent Schmidt number', &
    NGRID      = 4,                          &
    NTYPE      = TYPEREAL,                   &
    NDIMS      = 3,                          &
    LTIMEDEP   = .TRUE.                      )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZPSI3)
!
!
! stores the Turbulent Schmidt number for the scalar variables
! 
  TZFIELD = TFIELDMETADATA(                    &
    CMNHNAME   = 'generic for SV in turb_ver', & !Temporary name to ease identification
    CSTDNAME   = '',                           &
    CUNITS     = '1',                          &
    CDIR       = 'XY',                         &
    NGRID      = 4,                            &
    NTYPE      = TYPEREAL,                     &
    NDIMS      = 3,                            &
    LTIMEDEP   = .TRUE.                        )
  DO JSV=1,KSV
    WRITE(TZFIELD%CMNHNAME, '("PSI_SV_",I3.3)') JSV
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZPSI_SV(:,:,JSV))
  END DO
!
END IF
!
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TURB_VER',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER
END MODULE MODE_TURB_VER 
