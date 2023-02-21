!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_VER_THERMO_FLUX
IMPLICIT NONE
CONTAINS
SUBROUTINE TURB_VER_THERMO_FLUX(D,CST,CSTURB,TURBN,TLES,            &
                      KRR,KRRL,KRRI,KSV,KGRADIENTS,                 &
                      OOCEAN,ODEEPOC,OFLYER,                        &
                      OCOUPLES, OCOMPUTE_SRC,                       &
                      PEXPL,PTSTEP,HPROGRAM,                        &
                      TPFILE,                                       &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PRHODJ,PTHVREF,PHGRAD,PZS,                    &
                      PSFTHM,PSFRM,PSFTHP,PSFRP,                    &
                      PWM,PTHLM,PRM,PSVM,                           &
                      PTKEM,PLM,PLEPS,                              &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,    &
                      PBETA, PSQRT_TKE, PDTH_DZ, PDR_DZ, PRED2TH3,  &
                      PRED2R3, PRED2THR3, PBLL_O_E, PETHETA,        &
                      PEMOIST, PREDTH1, PREDR1, PPHI3, PPSI3, PD,   &
                      PFWTH,PFWR,PFTH2,PFR2,PFTHR,MFMOIST,PBL_DEPTH,&
                      PWTHV,PRTHLS,PRRS,PTHLP,PRP,PTP,PWTH,PWRC,    &
                      PSSTFL, PSSTFL_C, PSSRFL_C                    )
!     ###############################################################
!
!
!!****  *TURB_VER_THERMO_FLUX* -compute the source terms due to the vertical turbulent
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
!!		  and in section 7, a diagnostic computation of the W variance is
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
!!
!!      MXM,MXF,MYM,MYF,MZM,MZF
!!                             :  Shuman functions (mean operators)
!!      DXF,DYF,DZF,DZM
!!                             :  Shuman functions (difference operators)
!!
!!      SUBROUTINE TRIDIAG     : to compute the splitted implicit evolution
!!                               of a variable located at a mass point
!!
!!      SUBROUTINE TRIDIAG_WIND: to compute the splitted implicit evolution
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
!!           JPHEXT     : number of horizontal external points
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
!!      Modifications: Mar  31, 1998 (V. Masson) splits the routine TURB_VER_THERMO_FLUX
!!      Modifications: Oct  18, 2000 (V. Masson) LES computations
!!      Modifications: Dec  01, 2000 (V. Masson) conservation of energy from
!!                                               surface flux in 1DIM case
!!                                               when slopes are present
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     Feb  20, 2003 (JP Pinty)  Add PFRAC_ICE
!!                     May  20, 2003 (JP Pinty)  Correction of ETHETA
!!                                                         and EMOIST calls
!!                     July     2005 (S. Tomas, V. Masson)
!!                                               Add 3rd order moments
!!                                               and implicitation of PHI3 and PSI3
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!                     2012-02 (Y. Seity) add possibility to run with reversed
!!                                             vertical levels
!!      Modifications  July 2015 (Wim de Rooy) TURBN%LHARAT switch
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                     2021 (D. Ricard) last version of HGRAD turbulence scheme
!!                                 Leronard terms instead of Reynolds terms
!!                                 applied to vertical fluxes of r_np and Thl
!!                                 for implicit version of turbulence scheme
!!                                 corrections and cleaning
!!      Modifications: June 2019 (Wim de Rooy) with energycascade, 50MF nog
!!                                             longer necessary
!!                     June 2020 (B. Vie) Patch preventing negative rc and ri in 2.3 and 3.3
!! JL Redelsperger  : 03/2021: Ocean and Autocoupling O-A LES Cases
!!                             Sfc flux shape for LDEEPOC Case
!  P. Wautelet 30/11/2022: compute PWTH and PWRC only when needed
!!--------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1,   ONLY: JPRB
USE MODE_SHUMAN_PHY, ONLY: DZF_PHY, DZM_PHY, MXF_PHY, MYF_PHY, MZF_PHY, MZM_PHY
USE YOMHOOK,    ONLY: LHOOK, DR_HOOK
!
USE MODD_CST,              ONLY: CST_t
USE MODD_CTURB,            ONLY: CSTURB_t
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_t
USE MODD_FIELD,            ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_IO,               ONLY: TFILEDATA
USE MODD_LES,              ONLY: TLES_t
USE MODD_PARAMETERS,       ONLY: JPVEXT_TURB, JPHEXT, XUNDEF
USE MODD_TURB_n,           ONLY: TURB_t
!
USE MODE_GRADIENT_W_PHY, ONLY: GZ_W_M_PHY
USE MODE_IO_FIELD_WRITE_PHY, ONLY: IO_FIELD_WRITE_PHY
USE MODE_PRANDTL
USE MODE_TM06_H,         ONLY: TM06_H
USE MODE_TRIDIAG_THERMO, ONLY: TRIDIAG_THERMO
!
USE MODI_LES_MEAN_SUBGRID_PHY
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
TYPE(TURB_t),           INTENT(IN)   :: TURBN
TYPE(TLES_t),           INTENT(INOUT):: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV           ! number of scalar var.
INTEGER,                INTENT(IN)   :: KGRADIENTS    ! Number of stored horizontal gradients
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  ODEEPOC      ! activates sfc forcing for ideal ocean deep conv
LOGICAL,                INTENT(IN)   ::  OFLYER       ! MesoNH flyer diagnostic
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES version 
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and
REAL,                   INTENT(IN)   ::  PEXPL        ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
CHARACTER(LEN=6),       INTENT(IN)   :: HPROGRAM      ! CPROGRAM is the program currently running (modd_conf)
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PDZZ, PDXX, PDYY, PDZX, PDZY ! Metric coefficients
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the normal to the ground surface
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PZZ          ! altitudes
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTHVREF      ! ref. state Virtual Potential Temperature
!
REAL, DIMENSION(D%NIJT,D%NKT,KGRADIENTS),INTENT(IN) :: PHGRAD  ! horizontal gradients
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PZS ! orography (for LEONARD terms)
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSFTHM,PSFRM ! surface fluxes at time t - deltat
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSFTHP,PSFRP ! surface fluxes at time t + deltat
!
REAL, DIMENSION(D%NIJT,D%NKT),     INTENT(IN) ::  PWM          ! Vertical wind
REAL, DIMENSION(D%NIJT,D%NKT),     INTENT(IN) ::  PTHLM        ! potential temperature at t-Delta t
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN) ::  PRM          ! Mixing ratios at t-Delta t
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVM         ! Mixing ratios
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTKEM        ! TKE at time t
!
! In case TURBN%LHARAT=TRUE, PLM already includes all stability corrections
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLM          ! Turb. mixing length
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLEPS        ! dissipative length
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PATHETA      ! coefficients between
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(MERGE(D%NIT,0,OCOMPUTE_SRC)*MERGE(D%NJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(IN)   ::  PSRCM        ! normalized
! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFRAC_ICE    ! ri fraction of rc+ri
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PBETA        ! buoyancy coefficient
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PSQRT_TKE    ! sqrt(e)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PDTH_DZ      ! d(th)/dz
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PDR_DZ       ! d(rt)/dz
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRED2TH3     ! 3D Redeslperger number R*2_th
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRED2R3      ! 3D Redeslperger number R*2_r
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRED2THR3    ! 3D Redeslperger number R*2_thr
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PBLL_O_E     ! beta * Lk * Leps / tke
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PETHETA      ! Coefficient for theta in theta_v computation
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PEMOIST      ! Coefficient for r in theta_v computation
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PREDTH1      ! 1D Redelsperger number for Th
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PREDR1       ! 1D Redelsperger number for r
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PPHI3        ! Prandtl number for temperature
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PPSI3        ! Prandtl number for vapor
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PD           ! Denominator in Prandtl numbers
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFWTH        ! d(w'2th' )/dz (at flux point)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFWR         ! d(w'2r'  )/dz (at flux point)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFTH2        ! d(w'th'2 )/dz (at mass point)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFR2         ! d(w'r'2  )/dz (at mass point)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFTHR        ! d(w'th'r')/dz (at mass point)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  MFMOIST      ! moist mass flux dual scheme
REAL, DIMENSION(MERGE(D%NIT,0,TURBN%CTOM=='TM06'),&
                MERGE(D%NJT,0,TURBN%CTOM=='TM06')),   INTENT(INOUT)::  PBL_DEPTH    ! BL depth
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PWTHV         ! buoyancy flux
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRTHLS     ! cumulated source for theta
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT) :: PRRS     ! cumulated source for rt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: PTHLP      ! guess of thl at t+ deltat
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: PRP        ! guess of r at t+ deltat
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: PTP       ! Dynamic and thermal
                                                     ! TKE production terms
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: PWTH       ! heat flux
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: PWRC       ! cloud water flux
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSTFL    ! Time evol Flux of T at sea surface (LOCEAN and LCOUPLES)
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSTFL_C  ! O-A interface flux for theta(LOCEAN and LCOUPLES)
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSRFL_C  ! O-A interface flux for vapor (LOCEAN and LCOUPLES) 
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(D%NIJT,D%NKT)  ::  &
       ZA,       & ! work variable for wrc or LES computation
       ZFLXZ,    & ! vertical flux of the treated variable
       ZSOURCE,  & ! source of evolution for the treated variable
       ZKEFF,    & ! effectif diffusion coeff = LT * SQRT( TKE )
       ZF,       & ! Flux in dTh/dt =-dF/dz (evaluated at t-1)(or rt instead of Th)
       ZDFDDTDZ, & ! dF/d(dTh/dz)
       ZDFDDRDZ, & ! dF/d(dr/dz)
       Z3RDMOMENT,&  ! 3 order term in flux or variance equation
       ZF_LEONARD,&  ! Leonard terms
       ZRWTHL,    &
       ZRWRNP,    &
       ZCLD_THOLD,&
       ZALT,      &
       ZWORK1,ZWORK2, &
       ZWORK3,ZWORK4 ! working var. for shuman operators (array syntax)
!
INTEGER             :: IKB,IKE      ! I index values for the Beginning and End
                                    ! mass points of the domain in the 3 direct.
INTEGER             :: IKT,IKA,IKU  ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain
INTEGER             :: JIJ, JK ! loop indexes
INTEGER             :: IIJB, IIJE
INTEGER             :: IKL
!
REAL :: ZTIME1, ZTIME2
REAL :: ZFLPROV
INTEGER           :: JKM          ! vertical index loop
INTEGER           :: JSW
REAL :: ZSWA     ! index for time flux interpolation
!
INTEGER :: IIU, IJU
LOGICAL :: GUSERV   ! flag to use water
LOGICAL :: GFTH2    ! flag to use w'th'2
LOGICAL :: GFWTH    ! flag to use w'2th'
LOGICAL :: GFR2     ! flag to use w'r'2
LOGICAL :: GFWR     ! flag to use w'2r'
LOGICAL :: GFTHR    ! flag to use w'th'r'
TYPE(TFIELDMETADATA) :: TZFIELD
!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_FLUX',0,ZHOOK_HANDLE)
!
! Size for a given proc & a given model
IIU=D%NIT
IJU=D%NJT
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
IKA=D%NKA
IKU=D%NKU
IKTB=D%NKTB
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IKL=D%NKL
!
GUSERV = (KRR/=0)
!
!  compute the coefficients for the uncentred gradient computation near the ground
!
IF (TURBN%LHARAT) THEN
 ! LHARAT so TKE and length scales at half levels!
  !wc 50MF can be omitted with energy cascade included
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZKEFF(IIJB:IIJE,1:IKT) =  PLM(IIJB:IIJE,1:IKT) * SQRT(PTKEM(IIJB:IIJE,1:IKT))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK1(IIJB:IIJE,1:IKT) = PLM(IIJB:IIJE,1:IKT) * SQRT(PTKEM(IIJB:IIJE,1:IKT))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL MZM_PHY(D,ZWORK1,ZKEFF)
ENDIF
!
! Define a cloud mask with ri and rc (used after with a threshold) for Leonard terms
!
IF(TURBN%LLEONARD) THEN
  IF ( KRRL >= 1 ) THEN
    IF ( KRRI >= 1 ) THEN
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZCLD_THOLD(IIJB:IIJE,1:IKT) = PRM(IIJB:IIJE,1:IKT,2) + PRM(IIJB:IIJE,1:IKT,4)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ELSE
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZCLD_THOLD(IIJB:IIJE,1:IKT) = PRM(IIJB:IIJE,1:IKT,2)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    END IF
  END IF
END IF
!
! Flags for 3rd order quantities
!
GFTH2 = .FALSE.
GFR2  = .FALSE.
GFTHR = .FALSE.
GFWTH = .FALSE.
GFWR  = .FALSE.
!
IF (TURBN%CTOM/='NONE') THEN
  GFTH2 = ANY(PFTH2/=0.)
  GFR2  = ANY(PFR2 /=0.) .AND. GUSERV
  GFTHR = ANY(PFTHR/=0.) .AND. GUSERV
  GFWTH = ANY(PFWTH/=0.)
  GFWR  = ANY(PFWR /=0.) .AND. GUSERV
END IF
!----------------------------------------------------------------------------
!
!*       2.   SOURCES OF CONSERVATIVE POTENTIAL TEMPERATURE AND
!                                                  PARTIAL THERMAL PRODUCTION
!             ---------------------------------------------------------------
!
!*       2.1  Splitted value for cons. potential temperature at t+deltat
!
! Compute the turbulent flux F and F' at time t-dt.
!
CALL DZM_PHY(D,PTHLM,ZWORK1)
CALL D_PHI3DTDZ_O_DDTDZ(D,CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,TURBN%CTURBDIM,GUSERV,ZWORK2)
IF (TURBN%LHARAT) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZF(IIJB:IIJE,1:IKT) = -ZKEFF(IIJB:IIJE,1:IKT)*ZWORK1(IIJB:IIJE,1:IKT)/PDZZ(IIJB:IIJE,1:IKT)
  ZDFDDTDZ(IIJB:IIJE,1:IKT) = -ZKEFF(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZF(IIJB:IIJE,1:IKT) = -CSTURB%XCSHF*PPHI3(IIJB:IIJE,1:IKT)*ZKEFF(IIJB:IIJE,1:IKT)& 
                                *ZWORK1(IIJB:IIJE,1:IKT)/PDZZ(IIJB:IIJE,1:IKT)
  ZDFDDTDZ(IIJB:IIJE,1:IKT) = -CSTURB%XCSHF*ZKEFF(IIJB:IIJE,1:IKT)*ZWORK2(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
IF (TURBN%LLEONARD) THEN
! ! Compute the Leonard terms for thl
  CALL MXF_PHY(D,PHGRAD(:,:,1),ZWORK1) ! GX_W_UW(XWT)
  CALL MZM_PHY(D,PHGRAD(:,:,3),ZWORK2) ! GX_M_M(PTHLM
  CALL MYF_PHY(D,PHGRAD(:,:,2),ZWORK3) ! GY_W_VW(PWM)
  CALL MZM_PHY(D,PHGRAD(:,:,4),ZWORK4) ! GY_M_M(PTHLM)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)  
  ZF_LEONARD(IIJB:IIJE,1:IKT)= TURBN%XCOEFHGRADTHL*PDXX(IIJB:IIJE,1:IKT)*PDYY(IIJB:IIJE,1:IKT)/12.0*( &
                     ZWORK1(IIJB:IIJE,1:IKT)*ZWORK2(IIJB:IIJE,1:IKT) &
                   + ZWORK3(IIJB:IIJE,1:IKT)*ZWORK4(IIJB:IIJE,1:IKT))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
! Effect of 3rd order terms in temperature flux (at flux point)
!
! d(w'2th')/dz
IF (GFWTH) THEN
  CALL M3_WTH_W2TH(D,CSTURB,PREDTH1,PREDR1,PD,ZKEFF,PTKEM,Z3RDMOMENT)
  CALL D_M3_WTH_W2TH_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,&
   & PD,PBLL_O_E,PETHETA,ZKEFF,PTKEM,ZWORK1)
!
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZF(IIJB:IIJE,1:IKT)= ZF(IIJB:IIJE,1:IKT) + Z3RDMOMENT(IIJB:IIJE,1:IKT) * PFWTH(IIJB:IIJE,1:IKT)
  ZDFDDTDZ(IIJB:IIJE,1:IKT) = ZDFDDTDZ(IIJB:IIJE,1:IKT) + ZWORK1(IIJB:IIJE,1:IKT) &
                                      * PFWTH(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
! d(w'th'2)/dz
IF (GFTH2) THEN
  CALL M3_WTH_WTH2(D,CSTURB,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA,Z3RDMOMENT)
  CALL D_M3_WTH_WTH2_O_DDTDZ(D,CSTURB,Z3RDMOMENT,PREDTH1,PREDR1,&
    & PD,PBLL_O_E,PETHETA,ZWORK1)
  CALL MZM_PHY(D,PFTH2,ZWORK2)
!
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZF(IIJB:IIJE,1:IKT) = ZF(IIJB:IIJE,1:IKT) + Z3RDMOMENT(IIJB:IIJE,1:IKT) &
                                * ZWORK2(IIJB:IIJE,1:IKT)
  ZDFDDTDZ(IIJB:IIJE,1:IKT) = ZDFDDTDZ(IIJB:IIJE,1:IKT) + ZWORK1(IIJB:IIJE,1:IKT) &
                                      * ZWORK2(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
! d(w'2r')/dz
IF (GFWR) THEN
  CALL M3_WTH_W2R(D,CSTURB,PD,ZKEFF,PTKEM,PBLL_O_E,PEMOIST,PDTH_DZ,ZWORK1)
  CALL D_M3_WTH_W2R_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,ZKEFF,PTKEM,PBLL_O_E,PEMOIST,ZWORK2)
!
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZF(IIJB:IIJE,1:IKT) = ZF(IIJB:IIJE,1:IKT) + ZWORK1(IIJB:IIJE,1:IKT) * PFWR(IIJB:IIJE,1:IKT)
  ZDFDDTDZ(IIJB:IIJE,1:IKT) = ZDFDDTDZ(IIJB:IIJE,1:IKT) + ZWORK2(IIJB:IIJE,1:IKT) &
                                      * PFWR(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
! d(w'r'2)/dz
IF (GFR2) THEN
  CALL M3_WTH_WR2(D,CSTURB,PD,ZKEFF,PTKEM,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PEMOIST,PDTH_DZ,ZWORK1)
  CALL MZM_PHY(D,PFR2,ZWORK2)
  CALL D_M3_WTH_WR2_O_DDTDZ(D,CSTURB,PREDTH1,PREDR1,PD,&
    & ZKEFF,PTKEM,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PEMOIST,ZWORK3)
!
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)    
  ZF(IIJB:IIJE,1:IKT) = ZF(IIJB:IIJE,1:IKT) + ZWORK1(IIJB:IIJE,1:IKT) * ZWORK2(IIJB:IIJE,1:IKT)
  ZDFDDTDZ(IIJB:IIJE,1:IKT) = ZDFDDTDZ(IIJB:IIJE,1:IKT) + ZWORK3(IIJB:IIJE,1:IKT) &
                                      * ZWORK2(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
! d(w'th'r')/dz
IF (GFTHR) THEN
  CALL M3_WTH_WTHR(D,CSTURB,PREDR1,PD,ZKEFF,PTKEM,PSQRT_TKE,PBETA,&
    & PLEPS,PEMOIST,Z3RDMOMENT)
  CALL D_M3_WTH_WTHR_O_DDTDZ(D,CSTURB,Z3RDMOMENT,PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA,ZWORK1)
  CALL MZM_PHY(D,PFTHR, ZWORK2)
!
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZF(IIJB:IIJE,1:IKT) = ZF(IIJB:IIJE,1:IKT) + Z3RDMOMENT(IIJB:IIJE,1:IKT) &
                                * ZWORK2(IIJB:IIJE,1:IKT)
  ZDFDDTDZ(IIJB:IIJE,1:IKT) = ZDFDDTDZ(IIJB:IIJE,1:IKT) + ZWORK1(IIJB:IIJE,1:IKT) &
                                      * ZWORK2(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
! specialcase for surface
IF (OOCEAN) THEN    ! ocean model in coupled case
  !$mnh_expand_array(JIJ=IIJB:IIJE) 
  ZF(IIJB:IIJE,IKE+1) =  PSFTHM(IIJB:IIJE) &
                *0.5* ( 1. + PRHODJ(IIJB:IIJE,IKU)/PRHODJ(IIJB:IIJE,IKE) )
  !$mnh_end_expand_array(JIJ=IIJB:IIJE) 
ELSE ! atmosp bottom
  !*In 3D, a part of the flux goes vertically,
  ! and another goes horizontally (in presence of slopes)
  !*In 1D, part of energy released in horizontal flux is taken into account in the vertical part
  IF (TURBN%CTURBDIM=='3DIM') THEN
    !$mnh_expand_array(JIJ=IIJB:IIJE) 
    ZF(IIJB:IIJE,IKB) = ( TURBN%XIMPL*PSFTHP(IIJB:IIJE) + PEXPL*PSFTHM(IIJB:IIJE) )   &
                       * PDIRCOSZW(IIJB:IIJE)                       &
                       * 0.5 * (1. + PRHODJ(IIJB:IIJE,IKA) / PRHODJ(IIJB:IIJE,IKB))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE) 
  ELSE
    !$mnh_expand_array(JIJ=IIJB:IIJE) 
    ZF(IIJB:IIJE,IKB) = ( TURBN%XIMPL*PSFTHP(IIJB:IIJE) + PEXPL*PSFTHM(IIJB:IIJE) )   &
                       / PDIRCOSZW(IIJB:IIJE)                       &
                       * 0.5 * (1. + PRHODJ(IIJB:IIJE,IKA) / PRHODJ(IIJB:IIJE,IKB))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE) 
  END IF
!
    ! atmos top
#ifdef REPRO48
#else
      ZF(IIJB:IIJE,IKE+1)=0.
#endif
END IF
!
! Compute the split conservative potential temperature at t+deltat
CALL TRIDIAG_THERMO(D,PTHLM,ZF,ZDFDDTDZ,PTSTEP,TURBN%XIMPL,PDZZ,&
                    PRHODJ,PTHLP)
!
! Compute the equivalent tendency for the conservative potential temperature
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)    
ZRWTHL(IIJB:IIJE,1:IKT)= PRHODJ(IIJB:IIJE,1:IKT)*(PTHLP(IIJB:IIJE,1:IKT)-PTHLM(IIJB:IIJE,1:IKT))& 
                                 /PTSTEP
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)    
! replace the flux by the Leonard terms above ZALT and ZCLD_THOLD
IF (TURBN%LLEONARD) THEN
 DO JK=1,IKT
   !$mnh_expand_array(JIJ=IIJB:IIJE)
   ZALT(IIJB:IIJE,JK) = PZZ(IIJB:IIJE,JK)-PZS(IIJB:IIJE) 
   !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
 END DO
 CALL MZM_PHY(D,PRHODJ,ZWORK1)
 !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
 ZWORK2(IIJB:IIJE,1:IKT) = ZWORK1(IIJB:IIJE,1:IKT)*ZF_LEONARD(IIJB:IIJE,1:IKT)
 !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
 CALL GZ_W_M_PHY(D,ZWORK2,PDZZ,ZWORK3)
 !$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
 WHERE ( (ZCLD_THOLD(IIJB:IIJE,1:IKT) >= TURBN%XCLDTHOLD) .AND. ( ZALT(IIJB:IIJE,1:IKT) >= TURBN%XALTHGRAD) )
  ZRWTHL(IIJB:IIJE,1:IKT) = -ZWORK3(IIJB:IIJE,1:IKT)
  PTHLP(IIJB:IIJE,1:IKT)=PTHLM(IIJB:IIJE,1:IKT)+PTSTEP*ZRWTHL(IIJB:IIJE,1:IKT)/PRHODJ(IIJB:IIJE,1:IKT)
 END WHERE
 !$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZWORK1(IIJB:IIJE,1:IKT) = PTHLP(IIJB:IIJE,1:IKT) - PTHLM(IIJB:IIJE,1:IKT)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
CALL DZM_PHY(D,ZWORK1,ZWORK2)
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PRTHLS(IIJB:IIJE,1:IKT)= PRTHLS(IIJB:IIJE,1:IKT)  + ZRWTHL(IIJB:IIJE,1:IKT)
!
!*       2.2  Partial Thermal Production
!
!  Conservative potential temperature flux :
!
!
ZFLXZ(IIJB:IIJE,1:IKT)   = ZF(IIJB:IIJE,1:IKT) + TURBN%XIMPL * ZDFDDTDZ(IIJB:IIJE,1:IKT) * & 
                                   ZWORK2(IIJB:IIJE,1:IKT)/ PDZZ(IIJB:IIJE,1:IKT)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
! replace the flux by the Leonard terms
IF (TURBN%LLEONARD) THEN
 !$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
 WHERE ( (ZCLD_THOLD(IIJB:IIJE,1:IKT) >= TURBN%XCLDTHOLD) .AND. ( ZALT(IIJB:IIJE,1:IKT) >= TURBN%XALTHGRAD) )
  ZFLXZ(IIJB:IIJE,1:IKT) = ZF_LEONARD(IIJB:IIJE,1:IKT)
 END WHERE
 !$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
IF (OOCEAN) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZFLXZ(IIJB:IIJE,IKE+1) = ZFLXZ(IIJB:IIJE,IKE)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZFLXZ(IIJB:IIJE,IKA) = ZFLXZ(IIJB:IIJE,IKB)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END IF
!
IF ( OFLYER ) THEN
  PWTH(:,:IKTB) = XUNDEF
  PWTH(:,IKTE:) = XUNDEF
  !
  DO JK = IKTB + 1, IKTE - 1
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    PWTH(IIJB:IIJE,JK)=0.5*(ZFLXZ(IIJB:IIJE,JK)+ZFLXZ(IIJB:IIJE,JK+IKL))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  END DO
  !
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  PWTH(IIJB:IIJE,IKB)=0.5*(ZFLXZ(IIJB:IIJE,IKB)+ZFLXZ(IIJB:IIJE,IKB+IKL)) 
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  !
  IF (OOCEAN) THEN
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    PWTH(IIJB:IIJE,IKE)=0.5*(ZFLXZ(IIJB:IIJE,IKE)+ZFLXZ(IIJB:IIJE,IKE+IKL))
    PWTH(IIJB:IIJE,IKA)=0.
    PWTH(IIJB:IIJE,IKU)=PWTH(IIJB:IIJE,IKE)! not used
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  ELSE
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    PWTH(IIJB:IIJE,IKA)=0.5*(ZFLXZ(IIJB:IIJE,IKA)+ZFLXZ(IIJB:IIJE,IKA+IKL))
    PWTH(IIJB:IIJE,IKE)=PWTH(IIJB:IIJE,IKE-IKL)
    PWTH(IIJB:IIJE,IKU)=0.
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  END IF
END IF
!
IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
  ! stores the conservative potential temperature vertical flux
  TZFIELD = TFIELDMETADATA(                                         &
   CMNHNAME   = 'THW_FLX',                                          &
   CSTDNAME   = '',                                                 &
   CLONGNAME  = 'THW_FLX',                                          &
   CUNITS     = 'K m s-1',                                          &
   CDIR       = 'XY',                                               &
   CCOMMENT   = 'Conservative potential temperature vertical flux', &
   NGRID      = 4,                                                  &
   NTYPE      = TYPEREAL,                                           &
   NDIMS      = 3,                                                  &
   LTIMEDEP   = .TRUE.                                              )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXZ)
END IF
!
! Contribution of the conservative temperature flux to the buoyancy flux
IF (OOCEAN) THEN
  CALL MZF_PHY(D,ZFLXZ,ZWORK1)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PTP(IIJB:IIJE,1:IKT)= CST%XG*CST%XALPHAOC * ZWORK1(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  IF (KRR /= 0) THEN
    CALL MZM_PHY(D,PETHETA,ZWORK1)
    ZWORK1(IIJB:IIJE,1:IKT) = ZWORK1(IIJB:IIJE,1:IKT) * ZFLXZ(IIJB:IIJE,1:IKT)
    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    !ZWORK1 = MZF( MZM(PETHETA,IKA, IKU, IKL) * ZFLXZ,IKA, IKU, IKL )
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PTP(IIJB:IIJE,1:IKT)  =  PBETA(IIJB:IIJE,1:IKT) * ZWORK2(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    PTP(IIJB:IIJE,IKB)=  PBETA(IIJB:IIJE,IKB) * PETHETA(IIJB:IIJE,IKB) *   &
                   0.5 * ( ZFLXZ(IIJB:IIJE,IKB) + ZFLXZ(IIJB:IIJE,IKB+IKL) )
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  ELSE
    CALL MZF_PHY(D,ZFLXZ,ZWORK1)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PTP(IIJB:IIJE,1:IKT)=  PBETA(IIJB:IIJE,1:IKT) * ZWORK1(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
END IF
!
! Buoyancy flux at flux points
!
CALL MZM_PHY(D,PETHETA,ZWORK1)
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PWTHV(IIJB:IIJE,1:IKT) = ZWORK1(IIJB:IIJE,1:IKT) * ZFLXZ(IIJB:IIJE,1:IKT)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$mnh_expand_array(JIJ=IIJB:IIJE)
PWTHV(IIJB:IIJE,IKB) = PETHETA(IIJB:IIJE,IKB) * ZFLXZ(IIJB:IIJE,IKB)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
IF (OOCEAN) THEN
  ! temperature contribution to Buy flux
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  PWTHV(IIJB:IIJE,IKE) = PETHETA(IIJB:IIJE,IKE) * ZFLXZ(IIJB:IIJE,IKE)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END IF
!*       2.3  Partial vertical divergence of the < Rc w > flux
! Correction for qc and qi negative in AROME
IF(HPROGRAM/='AROME  ') THEN
 IF ( KRRL >= 1 ) THEN
   !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
   ZWORK1(IIJB:IIJE,1:IKT) = ZFLXZ(IIJB:IIJE,1:IKT)/PDZZ(IIJB:IIJE,1:IKT)
   !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
   CALL DZF_PHY(D,ZWORK1,ZWORK2)
   IF ( KRRI >= 1 ) THEN
     !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
     PRRS(IIJB:IIJE,1:IKT,2) = PRRS(IIJB:IIJE,1:IKT,2) -                                        &
                     PRHODJ(IIJB:IIJE,1:IKT)*PATHETA(IIJB:IIJE,1:IKT)*2.*PSRCM(IIJB:IIJE,1:IKT)& 
                     *ZWORK2(IIJB:IIJE,1:IKT) *(1.0-PFRAC_ICE(IIJB:IIJE,1:IKT))
     PRRS(IIJB:IIJE,1:IKT,4) = PRRS(IIJB:IIJE,1:IKT,4) -                                        &
                     PRHODJ(IIJB:IIJE,1:IKT)*PATHETA(IIJB:IIJE,1:IKT)*2.*PSRCM(IIJB:IIJE,1:IKT)&
                     * ZWORK2(IIJB:IIJE,1:IKT)*PFRAC_ICE(IIJB:IIJE,1:IKT)
     !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
   ELSE
     !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
     PRRS(IIJB:IIJE,1:IKT,2) = PRRS(IIJB:IIJE,1:IKT,2) -                                        &
                     PRHODJ(IIJB:IIJE,1:IKT)*PATHETA(IIJB:IIJE,1:IKT)*2.*PSRCM(IIJB:IIJE,1:IKT)&
                     *ZWORK2(IIJB:IIJE,1:IKT)
     !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
   END IF
 END IF
END IF
!
!*       2.4  Storage in LES configuration
!
IF (TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  !
  CALL MZF_PHY(D,ZFLXZ,ZWORK1)
  !
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK1, TLES%X_LES_SUBGRID_WThl )
  !
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK2(IIJB:IIJE,1:IKT) = PWM(IIJB:IIJE,1:IKT)*ZFLXZ(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL MZF_PHY(D,ZWORK2,ZWORK3)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK3, TLES%X_LES_RES_W_SBG_WThl )
  !
  CALL GZ_W_M_PHY(D,PWM,PDZZ,ZWORK2)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK3(IIJB:IIJE,1:IKT) = ZWORK2(IIJB:IIJE,1:IKT) * ZWORK1(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK3, TLES%X_LES_RES_ddxa_W_SBG_UaThl )
  !
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK2(IIJB:IIJE,1:IKT) = PDTH_DZ(IIJB:IIJE,1:IKT)*ZFLXZ(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL MZF_PHY(D,ZWORK2,ZWORK3)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK3, TLES%X_LES_RES_ddxa_Thl_SBG_UaThl )
  !
  CALL MZM_PHY(D,PETHETA,ZWORK2)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK3(IIJB:IIJE,1:IKT) = ZWORK2(IIJB:IIJE,1:IKT) * ZFLXZ(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL MZF_PHY(D,ZWORK3,ZWORK4)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK4, TLES%X_LES_SUBGRID_WThv , .TRUE. ) 
  !
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK2(IIJB:IIJE,1:IKT) = -CSTURB%XCTP*PSQRT_TKE(IIJB:IIJE,1:IKT)/PLM(IIJB:IIJE,1:IKT) &
                                    *ZWORK1(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2, TLES%X_LES_SUBGRID_ThlPz ) 
  !
  IF (KRR>=1) THEN
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK2(IIJB:IIJE,1:IKT) = PDR_DZ(IIJB:IIJE,1:IKT)*ZFLXZ(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL MZF_PHY(D,ZWORK2,ZWORK3)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK3, TLES%X_LES_RES_ddxa_Rt_SBG_UaThl )
  END IF
  !
  !* diagnostic of mixing coefficient for heat
  CALL DZM_PHY(D,PTHLP,ZA)
  !$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
  WHERE (ZA(IIJB:IIJE,1:IKT)==0.) 
    ZA(IIJB:IIJE,1:IKT)=1.E-6
  END WHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZA(IIJB:IIJE,1:IKT) = - ZFLXZ(IIJB:IIJE,1:IKT) / ZA(IIJB:IIJE,1:IKT) * PDZZ(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZA(IIJB:IIJE,IKB) = CSTURB%XCSHF*PPHI3(IIJB:IIJE,IKB)*ZKEFF(IIJB:IIJE,IKB)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  CALL MZF_PHY(D,ZA,ZA)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZA(IIJB:IIJE,1:IKT) = MIN(MAX(ZA(IIJB:IIJE,1:IKT),-1000.),1000.)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZA, TLES%X_LES_SUBGRID_Kh )
  !
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!*       2.5  New boundary layer depth for TOMs
!
IF (TURBN%CTOM=='TM06') CALL TM06_H(D,PTSTEP,PZZ,ZFLXZ,PBL_DEPTH)
!
!----------------------------------------------------------------------------
!
!
!*       3.   SOURCES OF CONSERVATIVE AND CLOUD MIXING RATIO AND
!                                        COMPLETE THERMAL PRODUCTION
!             ------------------------------------------------------
!
!*       3.1  Splitted value for cons. mixing ratio at t+deltat
!
!
IF (KRR /= 0) THEN
  ! Compute the turbulent flux F and F' at time t-dt.
  !
  CALL DZM_PHY(D,PRM(:,:,1),ZWORK1)
 IF (TURBN%LHARAT) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)    
  ZF(IIJB:IIJE,1:IKT) = -ZKEFF(IIJB:IIJE,1:IKT)*ZWORK1(IIJB:IIJE,1:IKT)/PDZZ(IIJB:IIJE,1:IKT)
  ZDFDDRDZ(IIJB:IIJE,1:IKT) = -ZKEFF(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)    
 ELSE
  CALL D_PSI3DRDZ_O_DDRDZ(D,CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,TURBN%CTURBDIM,GUSERV,ZWORK2)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)    
  ZF(IIJB:IIJE,1:IKT) = -CSTURB%XCSHF*PPSI3(IIJB:IIJE,1:IKT)*ZKEFF(IIJB:IIJE,1:IKT)& 
                                *ZWORK1(IIJB:IIJE,1:IKT)/PDZZ(IIJB:IIJE,1:IKT)
  ZDFDDRDZ(IIJB:IIJE,1:IKT) = -CSTURB%XCSHF*ZKEFF(IIJB:IIJE,1:IKT)*ZWORK2(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)    
 ENDIF
  !
  ! Compute Leonard Terms for Cloud mixing ratio
 IF (TURBN%LLEONARD) THEN
  CALL MXF_PHY(D,PHGRAD(:,:,1),ZWORK1) ! GX_W_UW(PWM)
  CALL MZM_PHY(D,PHGRAD(:,:,5),ZWORK2) ! GX_M_M(PRM)
  CALL MYF_PHY(D,PHGRAD(:,:,2),ZWORK3) ! GY_W_VW(PWM)
  CALL MZM_PHY(D,PHGRAD(:,:,6),ZWORK4) ! GY_M_M(PRM)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)  
  ZF_LEONARD(IIJB:IIJE,1:IKT)= TURBN%XCOEFHGRADTHL*PDXX(IIJB:IIJE,1:IKT)*PDYY(IIJB:IIJE,1:IKT)/12.0*( &
                     ZWORK1(IIJB:IIJE,1:IKT)*ZWORK2(IIJB:IIJE,1:IKT) &
                   + ZWORK3(IIJB:IIJE,1:IKT)*ZWORK4(IIJB:IIJE,1:IKT))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
 END IF
  !
  ! Effect of 3rd order terms in temperature flux (at flux point)
  !
  ! d(w'2r')/dz
  IF (GFWR) THEN
    CALL M3_WR_W2R(D,CSTURB,PREDR1,PREDTH1,PD,ZKEFF,PTKEM,Z3RDMOMENT)
    CALL D_M3_WR_W2R_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,&
     & PBLL_O_E,PEMOIST,ZKEFF,PTKEM,ZWORK1)
  !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZF(IIJB:IIJE,1:IKT)= ZF(IIJB:IIJE,1:IKT) + Z3RDMOMENT(IIJB:IIJE,1:IKT) * PFWR(IIJB:IIJE,1:IKT)
    ZDFDDRDZ(IIJB:IIJE,1:IKT) = ZDFDDRDZ(IIJB:IIJE,1:IKT) + ZWORK1(IIJB:IIJE,1:IKT) &
                                        * PFWR(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
  !
  ! d(w'r'2)/dz
  IF (GFR2) THEN
    CALL M3_WR_WR2(D,CSTURB,PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST,Z3RDMOMENT)
    CALL MZM_PHY(D,PFR2,ZWORK1)
    CALL D_M3_WR_WR2_O_DDRDZ(D,CSTURB,Z3RDMOMENT,PREDR1,&
     & PREDTH1,PD,PBLL_O_E,PEMOIST,ZWORK2)
  !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZF(IIJB:IIJE,1:IKT) = ZF(IIJB:IIJE,1:IKT) + Z3RDMOMENT(IIJB:IIJE,1:IKT) &
                                  * ZWORK1(IIJB:IIJE,1:IKT)
    ZDFDDRDZ(IIJB:IIJE,1:IKT) = ZDFDDRDZ(IIJB:IIJE,1:IKT) + ZWORK2(IIJB:IIJE,1:IKT) &
                                        * ZWORK1(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
  !
  ! d(w'2th')/dz
  IF (GFWTH) THEN
    CALL M3_WR_W2TH(D,CSTURB,PD,ZKEFF,&
     & PTKEM,PBLL_O_E,PETHETA,PDR_DZ,ZWORK1)
    CALL D_M3_WR_W2TH_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,&
     & PD,ZKEFF,PTKEM,PBLL_O_E,PETHETA,ZWORK2)
  !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZF(IIJB:IIJE,1:IKT) = ZF(IIJB:IIJE,1:IKT) + ZWORK1(IIJB:IIJE,1:IKT) * PFWTH(IIJB:IIJE,1:IKT)
    ZDFDDRDZ(IIJB:IIJE,1:IKT) = ZDFDDRDZ(IIJB:IIJE,1:IKT) + ZWORK2(IIJB:IIJE,1:IKT) &
                                        * PFWTH(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
  !
  ! d(w'th'2)/dz
  IF (GFTH2) THEN
    CALL MZM_PHY(D,PFTH2,ZWORK1)
    CALL M3_WR_WTH2(D,CSTURB,PD,ZKEFF,PTKEM,&
    & PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA,PDR_DZ,ZWORK2)
    CALL D_M3_WR_WTH2_O_DDRDZ(D,CSTURB,PREDR1,PREDTH1,PD,&
     &ZKEFF,PTKEM,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA,ZWORK3)
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZF(IIJB:IIJE,1:IKT) = ZF(IIJB:IIJE,1:IKT) + ZWORK2(IIJB:IIJE,1:IKT) * ZWORK1(IIJB:IIJE,1:IKT)
    ZDFDDRDZ(IIJB:IIJE,1:IKT) = ZDFDDRDZ(IIJB:IIJE,1:IKT) + ZWORK3(IIJB:IIJE,1:IKT) &
                                        * ZWORK1(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
  !
  ! d(w'th'r')/dz
  IF (GFTHR) THEN
    CALL M3_WR_WTHR(D,CSTURB,PREDTH1,PD,ZKEFF,PTKEM,PSQRT_TKE,PBETA,&
     & PLEPS,PETHETA,Z3RDMOMENT)
    CALL MZM_PHY(D,PFTHR,ZWORK1)
    CALL D_M3_WR_WTHR_O_DDRDZ(D,CSTURB,Z3RDMOMENT,PREDR1, &
     & PREDTH1,PD,PBLL_O_E,PEMOIST,ZWORK2)
  !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZF(IIJB:IIJE,1:IKT) = ZF(IIJB:IIJE,1:IKT) + Z3RDMOMENT(IIJB:IIJE,1:IKT) &
                                  * ZWORK1(IIJB:IIJE,1:IKT)
    ZDFDDRDZ(IIJB:IIJE,1:IKT) = ZDFDDRDZ(IIJB:IIJE,1:IKT) + ZWORK2(IIJB:IIJE,1:IKT) &
                                        * ZWORK1(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
  !
   !special case at sfc
    IF (OOCEAN) THEN
      ! General ocean case
      ! salinity/evap effect to be added later !!!!!
      ZF(IIJB:IIJE,IKE) = 0.
    ELSE ! atmosp case 
    ! atmosp bottom
    !* in 3DIM case, a part of the flux goes vertically, and another goes horizontally
    ! (in presence of slopes)
    !* in 1DIM case, the part of energy released in horizontal flux
    ! is taken into account in the vertical part
    !
    IF (TURBN%CTURBDIM=='3DIM') THEN
      !$mnh_expand_array(JIJ=IIJB:IIJE) 
      ZF(IIJB:IIJE,IKB) = ( TURBN%XIMPL*PSFRP(IIJB:IIJE) + PEXPL*PSFRM(IIJB:IIJE) )       &
                           * PDIRCOSZW(IIJB:IIJE)                       &
                         * 0.5 * (1. + PRHODJ(IIJB:IIJE,IKA) / PRHODJ(IIJB:IIJE,IKB))
      !$mnh_end_expand_array(JIJ=IIJB:IIJE) 
    ELSE
      !$mnh_expand_array(JIJ=IIJB:IIJE) 
      ZF(IIJB:IIJE,IKB) = ( TURBN%XIMPL*PSFRP(IIJB:IIJE) + PEXPL*PSFRM(IIJB:IIJE) )     &
                         / PDIRCOSZW(IIJB:IIJE)                       &
                         * 0.5 * (1. + PRHODJ(IIJB:IIJE,IKA) / PRHODJ(IIJB:IIJE,IKB))
      !$mnh_end_expand_array(JIJ=IIJB:IIJE) 
    END IF
      ! atmos top
#ifdef REPRO48
#else
      ZF(IIJB:IIJE,IKE+1)=0.
#endif
    END IF
  ! Compute the split conservative potential temperature at t+deltat
  CALL TRIDIAG_THERMO(D,PRM(:,:,1),ZF,ZDFDDRDZ,PTSTEP,TURBN%XIMPL,&
                      PDZZ,PRHODJ,PRP)
  !
  ! Compute the equivalent tendency for the conservative mixing ratio
  !
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)    
  ZRWRNP(IIJB:IIJE,1:IKT) = PRHODJ(IIJB:IIJE,1:IKT)*(PRP(IIJB:IIJE,1:IKT)-PRM(IIJB:IIJE,1:IKT,1))& 
                                   /PTSTEP
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)    
  !
  ! replace the flux by the Leonard terms above ZALT and ZCLD_THOLD
  IF (TURBN%LLEONARD) THEN
   CALL MZM_PHY(D,PRHODJ,ZWORK1)
   !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
   ZWORK2(IIJB:IIJE,1:IKT) = ZWORK1(IIJB:IIJE,1:IKT)*ZF_LEONARD(IIJB:IIJE,1:IKT)
   !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
   CALL GZ_W_M_PHY(D,ZWORK2,PDZZ,ZWORK3)
   !$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
   WHERE ( (ZCLD_THOLD(IIJB:IIJE,1:IKT) >= TURBN%XCLDTHOLD ) .AND. ( ZALT(IIJB:IIJE,1:IKT) >= TURBN%XALTHGRAD ) )
    ZRWRNP(IIJB:IIJE,1:IKT) =  -ZWORK3(IIJB:IIJE,1:IKT)
    PRP(IIJB:IIJE,1:IKT)=PRM(IIJB:IIJE,1:IKT,1)+PTSTEP*ZRWTHL(IIJB:IIJE,1:IKT)/PRHODJ(IIJB:IIJE,1:IKT)
   END WHERE
   !$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
  !
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK1(IIJB:IIJE,1:IKT) = PRP(IIJB:IIJE,1:IKT) - PRM(IIJB:IIJE,1:IKT,1)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL DZM_PHY(D,ZWORK1,ZWORK2)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PRRS(IIJB:IIJE,1:IKT,1) = PRRS(IIJB:IIJE,1:IKT,1) + ZRWRNP(IIJB:IIJE,1:IKT)
  !
  !*       3.2  Complete thermal production
  !
  ! cons. mixing ratio flux :
  !
  ZFLXZ(IIJB:IIJE,1:IKT)   = ZF(IIJB:IIJE,1:IKT)                                                &
                 + TURBN%XIMPL * ZDFDDRDZ(IIJB:IIJE,1:IKT) * ZWORK2(IIJB:IIJE,1:IKT) / PDZZ(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !
  ! replace the flux by the Leonard terms above ZALT and ZCLD_THOLD
  IF (TURBN%LLEONARD) THEN
   !$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
   WHERE ( (ZCLD_THOLD(IIJB:IIJE,1:IKT) >= TURBN%XCLDTHOLD ) .AND. ( ZALT(IIJB:IIJE,1:IKT) >= TURBN%XALTHGRAD ) )
    ZFLXZ(IIJB:IIJE,1:IKT) = ZF_LEONARD(IIJB:IIJE,1:IKT)
   END WHERE
   !$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
  !
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZFLXZ(IIJB:IIJE,IKA) = ZFLXZ(IIJB:IIJE,IKB)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  !
  IF (OOCEAN) THEN
    ZFLXZ(IIJB:IIJE,IKU) = ZFLXZ(IIJB:IIJE,IKE)
  END IF
  !
  IF ( OFLYER ) THEN
    DO JK=IKTB+1,IKTE-1
      !$mnh_expand_array(JIJ=IIJB:IIJE)
      PWRC(IIJB:IIJE,JK)=0.5*(ZFLXZ(IIJB:IIJE,JK)+ZFLXZ(IIJB:IIJE,JK+IKL))
      !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    END DO
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    PWRC(IIJB:IIJE,IKB)=0.5*(ZFLXZ(IIJB:IIJE,IKB)+ZFLXZ(IIJB:IIJE,IKB+IKL))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    !
    IF (OOCEAN) THEN
      !$mnh_expand_array(JIJ=IIJB:IIJE)          
      PWRC(IIJB:IIJE,IKE)=0.5*(ZFLXZ(IIJB:IIJE,IKE)+ZFLXZ(IIJB:IIJE,IKE+IKL))
      PWRC(IIJB:IIJE,IKA)=0.
      PWRC(IIJB:IIJE,IKE+1)=ZFLXZ(IIJB:IIJE,IKE+1)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE)    
    ELSE
      !$mnh_expand_array(JIJ=IIJB:IIJE)
      PWRC(IIJB:IIJE,IKA)=0.5*(ZFLXZ(IIJB:IIJE,IKA)+ZFLXZ(IIJB:IIJE,IKA+IKL))
      PWRC(IIJB:IIJE,IKE)=PWRC(IIJB:IIJE,IKE-IKL)
      PWRC(IIJB:IIJE,IKU)=0.
      !$mnh_end_expand_array(JIJ=IIJB:IIJE)    
    END IF
  END IF
  !
  IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
    ! stores the conservative mixing ratio vertical flux
    TZFIELD = TFIELDMETADATA(                                 &
      CMNHNAME   = 'RCONSW_FLX',                              &
      CSTDNAME   = '',                                        &
      CLONGNAME  = 'RCONSW_FLX',                              &
      CUNITS     = 'kg m s-1 kg-1',                           &
      CDIR       = 'XY',                                      &
      CCOMMENT   = 'Conservative mixing ratio vertical flux', &
      NGRID      = 4,                                         &
      NTYPE      = TYPEREAL,                                  &
      NDIMS      = 3,                                         &
      LTIMEDEP   = .TRUE.                                     )
    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXZ)
  END IF
  !
  ! Contribution of the conservative water flux to the Buoyancy flux
  IF (OOCEAN) THEN
     CALL MZF_PHY(D,ZFLXZ,ZWORK1)
     !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
     ZA(IIJB:IIJE,1:IKT)=  -CST%XG*CST%XBETAOC  * ZWORK1(IIJB:IIJE,1:IKT)
     !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ELSE
    CALL MZM_PHY(D,PEMOIST,ZWORK1)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK1(IIJB:IIJE,1:IKT) = ZWORK1(IIJB:IIJE,1:IKT) * ZFLXZ(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZA(IIJB:IIJE,1:IKT)   =  PBETA(IIJB:IIJE,1:IKT) * ZWORK2(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZA(IIJB:IIJE,IKB) =  PBETA(IIJB:IIJE,IKB) * PEMOIST(IIJB:IIJE,IKB) *   &
                   0.5 * ( ZFLXZ(IIJB:IIJE,IKB) + ZFLXZ(IIJB:IIJE,IKB+IKL) )
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PTP(IIJB:IIJE,1:IKT) = PTP(IIJB:IIJE,1:IKT) + ZA(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  END IF
  !
  ! Buoyancy flux at flux points
  !
  CALL MZM_PHY(D,PEMOIST,ZWORK1)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PWTHV(IIJB:IIJE,1:IKT)=PWTHV(IIJB:IIJE,1:IKT) + ZWORK1(IIJB:IIJE,1:IKT) * ZFLXZ(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  PWTHV(IIJB:IIJE,IKB) = PWTHV(IIJB:IIJE,IKB) + PEMOIST(IIJB:IIJE,IKB) * ZFLXZ(IIJB:IIJE,IKB)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  IF (OOCEAN) THEN
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    PWTHV(IIJB:IIJE,IKE) = PWTHV(IIJB:IIJE,IKE) + PEMOIST(IIJB:IIJE,IKE)* ZFLXZ(IIJB:IIJE,IKE)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  END IF
!
!*       3.3  Complete vertical divergence of the < Rc w > flux
! Correction of qc and qi negative for AROME
IF(HPROGRAM/='AROME  ') THEN
   IF ( KRRL >= 1 ) THEN
       !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
       ZWORK2(IIJB:IIJE,1:IKT) = ZFLXZ(IIJB:IIJE,1:IKT) / &
       PDZZ(IIJB:IIJE,1:IKT)
       !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
       CALL DZF_PHY(D,ZWORK2,ZWORK1)
       !
     IF ( KRRI >= 1 ) THEN
       !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
       PRRS(IIJB:IIJE,1:IKT,2) = PRRS(IIJB:IIJE,1:IKT,2) -                                        &
                       PRHODJ(IIJB:IIJE,1:IKT)*PAMOIST(IIJB:IIJE,1:IKT)*2.*PSRCM(IIJB:IIJE,1:IKT)& 
                       *ZWORK1(IIJB:IIJE,1:IKT) *(1.0-PFRAC_ICE(IIJB:IIJE,1:IKT))
       PRRS(IIJB:IIJE,1:IKT,4) = PRRS(IIJB:IIJE,1:IKT,4) -                                        &
                       PRHODJ(IIJB:IIJE,1:IKT)*PAMOIST(IIJB:IIJE,1:IKT)*2.*PSRCM(IIJB:IIJE,1:IKT)&
                       *ZWORK1(IIJB:IIJE,1:IKT) *PFRAC_ICE(IIJB:IIJE,1:IKT)
       !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
     ELSE
       !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
       PRRS(IIJB:IIJE,1:IKT,2) = PRRS(IIJB:IIJE,1:IKT,2) -                                        &
                       PRHODJ(IIJB:IIJE,1:IKT)*PAMOIST(IIJB:IIJE,1:IKT)*2.*PSRCM(IIJB:IIJE,1:IKT)&
                       *ZWORK1(IIJB:IIJE,1:IKT)
       !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
     END IF
   END IF
END IF
!
!*       3.4  Storage in LES configuration
!
  IF (TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    !
    CALL MZF_PHY(D,ZFLXZ,ZWORK1)
    !
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK1, TLES%X_LES_SUBGRID_WRt )
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK2(IIJB:IIJE,1:IKT) = PWM(IIJB:IIJE,1:IKT)*ZFLXZ(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL MZF_PHY(D,ZWORK2,ZWORK3)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK3, TLES%X_LES_RES_W_SBG_WRt )
    !
    CALL GZ_W_M_PHY(D,PWM,PDZZ,ZWORK2)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK3(IIJB:IIJE,1:IKT) = ZWORK2(IIJB:IIJE,1:IKT) * ZWORK1(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK3, TLES%X_LES_RES_ddxa_W_SBG_UaRt )
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK2(IIJB:IIJE,1:IKT) = PDTH_DZ(IIJB:IIJE,1:IKT)*ZFLXZ(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL MZF_PHY(D,ZWORK2,ZWORK3)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK3, TLES%X_LES_RES_ddxa_Thl_SBG_UaRt )
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK2(IIJB:IIJE,1:IKT) = PDR_DZ(IIJB:IIJE,1:IKT)*ZFLXZ(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL MZF_PHY(D,ZWORK2,ZWORK3)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK3, TLES%X_LES_RES_ddxa_Rt_SBG_UaRt )
    !
    CALL MZM_PHY(D,PEMOIST,ZWORK2)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK3(IIJB:IIJE,1:IKT) = ZWORK2(IIJB:IIJE,1:IKT) * ZFLXZ(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL MZF_PHY(D,ZWORK3,ZWORK4)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK4, TLES%X_LES_SUBGRID_WThv , .TRUE. )
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK2(IIJB:IIJE,1:IKT) = -CSTURB%XCTP*PSQRT_TKE(IIJB:IIJE,1:IKT)/PLM(IIJB:IIJE,1:IKT) &
                                      *ZWORK1(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2, TLES%X_LES_SUBGRID_RtPz )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
END IF
!
!----------------------------------------------------------------------------
!
!
!*       4.   TURBULENT CORRELATIONS : <w Rc>
!             -------------------------------
!
!
!*       4.1  <w Rc>
!
IF ( ((TURBN%LTURB_FLX .AND. TPFILE%LOPENED) .OR. TLES%LLES_CALL) .AND. (KRRL > 0) ) THEN
!
! recover the Conservative potential temperature flux :
! With TURBN%LHARAT is true tke and length scales at half levels
! yet modify to use length scale and tke at half levels from vdfexcuhl
 !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
 ZWORK1(IIJB:IIJE,1:IKT) = TURBN%XIMPL * PTHLP(IIJB:IIJE,1:IKT) + PEXPL * PTHLM(IIJB:IIJE,1:IKT)
 !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
 CALL DZM_PHY(D,ZWORK1,ZWORK2)
 IF (TURBN%LHARAT) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZA(IIJB:IIJE,1:IKT)   = ZWORK2(IIJB:IIJE,1:IKT)/ PDZZ(IIJB:IIJE,1:IKT) * &
                                 (-PLM(IIJB:IIJE,1:IKT)*PSQRT_TKE(IIJB:IIJE,1:IKT))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
 ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK1(IIJB:IIJE,1:IKT) = PLM(IIJB:IIJE,1:IKT)*PSQRT_TKE(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL MZM_PHY(D,ZWORK1,ZWORK3)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZA(IIJB:IIJE,1:IKT)   = ZWORK2(IIJB:IIJE,1:IKT)/ PDZZ(IIJB:IIJE,1:IKT) * &
                                  (-PPHI3(IIJB:IIJE,1:IKT)*ZWORK3(IIJB:IIJE,1:IKT)) * CSTURB%XCSHF 
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
 ENDIF
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZA(IIJB:IIJE,IKB) = (TURBN%XIMPL*PSFTHP(IIJB:IIJE) + PEXPL*PSFTHM(IIJB:IIJE)) * PDIRCOSZW(IIJB:IIJE)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  !  
  ! compute <w Rc>
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK1(IIJB:IIJE,1:IKT) = PAMOIST(IIJB:IIJE,1:IKT) * 2.* PSRCM(IIJB:IIJE,1:IKT)
  ZWORK2(IIJB:IIJE,1:IKT) = PATHETA(IIJB:IIJE,1:IKT) * 2.* PSRCM(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL MZM_PHY(D,ZWORK1,ZWORK3)
  CALL MZM_PHY(D,ZWORK2,ZWORK4)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZFLXZ(IIJB:IIJE,1:IKT) = ZWORK3(IIJB:IIJE,1:IKT)* ZFLXZ(IIJB:IIJE,1:IKT) &
                                 + ZWORK4(IIJB:IIJE,1:IKT)* ZA(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZFLXZ(IIJB:IIJE,IKA) = ZFLXZ(IIJB:IIJE,IKB)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  !
  ! store the liquid water mixing ratio vertical flux
  IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
    TZFIELD = TFIELDMETADATA(                                 &
      CMNHNAME   = 'RCW_FLX',                                 &
      CSTDNAME   = '',                                        &
      CLONGNAME  = 'RCW_FLX',                                 &
      CUNITS     = 'kg m s-1 kg-1',                           &
      CDIR       = 'XY',                                      &
      CCOMMENT   = 'Liquid water mixing ratio vertical flux', &
      NGRID      = 4,                                         &
      NTYPE      = TYPEREAL,                                  &
      NDIMS      = 3,                                         &
      LTIMEDEP   = .TRUE.                                     )
    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXZ)
  END IF
  !
! and we store in LES configuration this subgrid flux <w'rc'>
!
  IF (TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL MZF_PHY(D,ZFLXZ,ZWORK1)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK1, TLES%X_LES_SUBGRID_WRc )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
END IF !end of <w Rc>
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_FLUX',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER_THERMO_FLUX
END MODULE MODE_TURB_VER_THERMO_FLUX
