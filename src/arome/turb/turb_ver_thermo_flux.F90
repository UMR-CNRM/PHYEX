!     ######spl
      SUBROUTINE TURB_VER_THERMO_FLUX(KKA,KKU,KKL,KRR,KRRL,KRRI,    &
                      OCLOSE_OUT,OTURB_FLX,HTURBDIM,HTOM,           &
                      PIMPL,PEXPL,                                  &
                      PTSTEP,                                       &
                      HFMFILE,HLUOUT,                               &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PRHODJ,PTHVREF,                               &
                      PSFTHM,PSFRM,PSFTHP,PSFRP,                    &
                      PWM,PTHLM,PRM,PSVM,                           &
                      PTKEM,PLM,PLEPS,                              &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,    &
                      PBETA, PSQRT_TKE, PDTH_DZ, PDR_DZ, PRED2TH3,  &
                      PRED2R3, PRED2THR3, PBLL_O_E, PETHETA,        &
                      PEMOIST, PREDTH1, PREDR1, PPHI3, PPSI3, PD,   &
                      PFWTH,PFWR,PFTH2,PFR2,PFTHR,MFMOIST,PBL_DEPTH,&
                      PWTHV,PRTHLS,PRRS,PTHLP,PRP,PTP,PWTH,PWRC)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      USE MODD_CTURB, ONLY : LHARAT
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
!!           XG         : gravity constant
!!
!!      Module MODD_CTURB: contains the set of constants for
!!                        the turbulence scheme
!!
!!           XCMFS,XCMFB : cts for the momentum flux
!!           XCSHF       : ct for the sensible heat flux
!!           XCHF        : ct for the moisture flux
!!           XCTV,XCHV   : cts for the T and moisture variances
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
!!      Modifications  July 2015 (Wim de Rooy) LHARAT switch
!!--------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_CTURB
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_LES
!
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_GRADIENT_M
USE MODI_SHUMAN , ONLY : DZF, DZM, MZF, MZM
USE MODI_TRIDIAG 
USE MODE_FMWRIT
USE MODI_LES_MEAN_SUBGRID
USE MODI_PRANDTL
USE MODI_TRIDIAG_THERMO
USE MODI_TM06_H
!
USE MODE_PRANDTL
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
!
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=ARO
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
LOGICAL,                INTENT(IN)   ::  OCLOSE_OUT   ! switch for syncronous
                                                      ! file opening       
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
CHARACTER*4,            INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
CHARACTER*4,            INTENT(IN)   ::  HTOM         ! type of Third Order Moment
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
CHARACTER(LEN=*),       INTENT(IN)   ::  HFMFILE      ! Name of the output
                                                      ! FM-file 
CHARACTER(LEN=*),       INTENT(IN)   ::  HLUOUT       ! Output-listing name for
                                                      ! model n
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDZZ, PDXX, PDYY, PDZX, PDZY
                                                      ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PZZ          ! altitudes
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  MFMOIST      ! moist mass flux dual scheme
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHVREF      ! ref. state Virtual 
                                                      ! Potential Temperature 
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSFTHM,PSFRM ! surface fluxes at time
!                                                     ! t - deltat 
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSFTHP,PSFRP ! surface fluxes at time
!                                                     ! t + deltat 
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PWM 
! Vertical wind
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHLM 
! potential temperature at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRM          ! Mixing ratios 
                                                      ! at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM         ! Mixing ratios 
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
!
! In case LHARAT=TRUE, PLM already includes all stability corrections
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLEPS        ! dissipative length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM        ! normalized 
! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFRAC_ICE    ! ri fraction of rc+ri
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PBETA        ! buoyancy coefficient
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSQRT_TKE    ! sqrt(e)
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDTH_DZ      ! d(th)/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDR_DZ       ! d(rt)/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRED2TH3     ! 3D Redeslperger number R*2_th
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRED2R3      ! 3D Redeslperger number R*2_r
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRED2THR3    ! 3D Redeslperger number R*2_thr
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PBLL_O_E     ! beta * Lk * Leps / tke
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PETHETA      ! Coefficient for theta in theta_v computation
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PEMOIST      ! Coefficient for r in theta_v computation
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PREDTH1      ! 1D Redelsperger number for Th
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PREDR1       ! 1D Redelsperger number for r
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PPHI3        ! Prandtl number for temperature
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PPSI3        ! Prandtl number for vapor
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PD           ! Denominator in Prandtl numbers
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFWTH        ! d(w'2th' )/dz (at flux point)
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFWR         ! d(w'2r'  )/dz (at flux point)
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFTH2        ! d(w'th'2 )/dz (at mass point)
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFR2         ! d(w'r'2  )/dz (at mass point)
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFTHR        ! d(w'th'r')/dz (at mass point)
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PBL_DEPTH    ! BL depth
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWTHV         ! buoyancy flux
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHLS     ! cumulated source for theta
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS       ! cumulated source for rt
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PTHLP      ! guess of thl at t+ deltat
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRP        ! guess of r at t+ deltat
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT)::  PTP       ! Dynamic and thermal
                                                     ! TKE production terms
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PWTH       ! heat flux
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PWRC       ! cloud water flux
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))  ::  &
       ZA,       & ! work variable for wrc or LES computation
       ZFLXZ,    & ! vertical flux of the treated variable
       ZSOURCE,  & ! source of evolution for the treated variable
       ZKEFF,    & ! effectif diffusion coeff = LT * SQRT( TKE )
       ZF,       & ! Flux in dTh/dt =-dF/dz (evaluated at t-1)(or rt instead of Th)
       ZDFDDTDZ, & ! dF/d(dTh/dz)
       ZDFDDRDZ, & ! dF/d(dr/dz)
       Z3RDMOMENT  ! 3 order term in flux or variance equation
INTEGER             :: IRESP        ! Return code of FM routines 
INTEGER             :: IGRID        ! C-grid indicator in LFIFM file 
INTEGER             :: ILENCH       ! Length of comment string in LFIFM file
INTEGER             :: IKB,IKE      ! I index values for the Beginning and End
                                    ! mass points of the domain in the 3 direct.
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain 
CHARACTER (LEN=100) :: YCOMMENT     ! comment string in LFIFM file
CHARACTER (LEN=16)  :: YRECFM       ! Name of the desired field in LFIFM file
!
REAL :: ZTIME1, ZTIME2
!
INTEGER :: JK
LOGICAL :: GUSERV   ! flag to use water
LOGICAL :: GFTH2    ! flag to use w'th'2
LOGICAL :: GFWTH    ! flag to use w'2th'
LOGICAL :: GFR2     ! flag to use w'r'2
LOGICAL :: GFWR     ! flag to use w'2r'
LOGICAL :: GFTHR    ! flag to use w'th'r'
!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_FLUX',0,ZHOOK_HANDLE)
IKT  =SIZE(PTHLM,3)  
IKTE =IKT-JPVEXT_TURB  
IKTB =1+JPVEXT_TURB               
IKB=KKA+JPVEXT_TURB*KKL
IKE=KKU-JPVEXT_TURB*KKL
!
GUSERV = (KRR/=0)
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!
IF (LHARAT) THEN
! LHARAT so TKE and length scales at half levels!
ZKEFF(:,:,:) =  PLM(:,:,:) * SQRT(PTKEM(:,:,:)) +50.*MFMOIST(:,:,:)
ELSE
ZKEFF(:,:,:) = MZM(PLM(:,:,:) * SQRT(PTKEM(:,:,:)), KKA, KKU, KKL)
ENDIF
!
!
! Flags for 3rd order quantities
!
GFTH2 = .FALSE.
GFR2  = .FALSE.
GFTHR = .FALSE.
GFWTH = .FALSE.
GFWR  = .FALSE.
!
IF (HTOM/='NONE') THEN
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
IF (LHARAT) THEN
ZF      (:,:,:) = -ZKEFF*DZM(PTHLM, KKA, KKU, KKL)/PDZZ
ZDFDDTDZ(:,:,:) = -ZKEFF
ELSE
ZF      (:,:,:) = -XCSHF*PPHI3*ZKEFF*DZM(PTHLM, KKA, KKU, KKL)/PDZZ
ZDFDDTDZ(:,:,:) = -XCSHF*ZKEFF*D_PHI3DTDZ_O_DDTDZ(PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,GUSERV)
ENDIF
!
! Effect of 3rd order terms in temperature flux (at flux point)
!
! d(w'2th')/dz
IF (GFWTH) THEN
  Z3RDMOMENT= M3_WTH_W2TH(KKA,KKU,KKL,PREDTH1,PREDR1,PD,ZKEFF,PTKEM)
!
  ZF       = ZF       + Z3RDMOMENT * PFWTH
  ZDFDDTDZ = ZDFDDTDZ + D_M3_WTH_W2TH_O_DDTDZ(KKA,KKU,KKL,PREDTH1,PREDR1,&
   & PD,PBLL_O_E,PETHETA,ZKEFF,PTKEM) * PFWTH
END IF
!
! d(w'th'2)/dz
IF (GFTH2) THEN
  Z3RDMOMENT= M3_WTH_WTH2(PREDTH1,PREDR1,PD,PBLL_O_E,PETHETA)
!
  ZF       = ZF       + Z3RDMOMENT * MZM(PFTH2, KKA, KKU, KKL)
  ZDFDDTDZ = ZDFDDTDZ + D_M3_WTH_WTH2_O_DDTDZ(Z3RDMOMENT,PREDTH1,PREDR1,&
    & PD,PBLL_O_E,PETHETA) * MZM(PFTH2, KKA, KKU, KKL)
END IF
!
! d(w'2r')/dz
IF (GFWR) THEN
  ZF       = ZF       + M3_WTH_W2R(KKA,KKU,KKL,PREDTH1,PREDR1,PD,ZKEFF,&
    & PTKEM,PBLL_O_E,PEMOIST,PDTH_DZ) * PFWR
  ZDFDDTDZ = ZDFDDTDZ + D_M3_WTH_W2R_O_DDTDZ(KKA,KKU,KKL,PREDTH1,PREDR1,&
    & PD,ZKEFF,PTKEM,PBLL_O_E,PEMOIST) * PFWR
END IF
!
! d(w'r'2)/dz
IF (GFR2) THEN
  ZF       = ZF       + M3_WTH_WR2(KKA,KKU,KKL,PREDTH1,PREDR1,PD,ZKEFF,PTKEM,&
    & PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PEMOIST,PDTH_DZ) * MZM(PFR2, KKA, KKU, KKL)
  ZDFDDTDZ = ZDFDDTDZ + D_M3_WTH_WR2_O_DDTDZ(KKA,KKU,KKL,PREDTH1,PREDR1,PD,&
    & ZKEFF,PTKEM,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PEMOIST) * MZM(PFR2, KKA, KKU, KKL)
END IF
!
! d(w'th'r')/dz
IF (GFTHR) THEN
  Z3RDMOMENT= M3_WTH_WTHR(KKA,KKU,KKL,PREDR1,PD,ZKEFF,PTKEM,PSQRT_TKE,PBETA,&
    & PLEPS,PEMOIST)
!
  ZF       = ZF       + Z3RDMOMENT * MZM(PFTHR, KKA, KKU, KKL)
  ZDFDDTDZ = ZDFDDTDZ + D_M3_WTH_WTHR_O_DDTDZ(Z3RDMOMENT,PREDTH1,&
    & PREDR1,PD,PBLL_O_E,PETHETA) * MZM(PFTHR, KKA, KKU, KKL)
END IF
!
!* in 3DIM case, a part of the flux goes vertically, and another goes horizontally
! (in presence of slopes)
!* in 1DIM case, the part of energy released in horizontal flux
! is taken into account in the vertical part
!
IF (HTURBDIM=='3DIM') THEN
  ZF(:,:,IKB) = ( PIMPL*PSFTHP(:,:) + PEXPL*PSFTHM(:,:) )   &
                     * PDIRCOSZW(:,:)                       &
                     * 0.5 * (1. + PRHODJ(:,:,KKA) / PRHODJ(:,:,IKB))
ELSE
  ZF(:,:,IKB) = ( PIMPL*PSFTHP(:,:) + PEXPL*PSFTHM(:,:) )   &
                     / PDIRCOSZW(:,:)                       &
                     * 0.5 * (1. + PRHODJ(:,:,KKA) / PRHODJ(:,:,IKB))
END IF
!
! Compute the splitted conservative potential temperature at t+deltat
CALL TRIDIAG_THERMO(KKA,KKU,KKL,PTHLM,ZF,ZDFDDTDZ,PTSTEP,PIMPL,PDZZ,&
                    PRHODJ,PTHLP)
!
! Compute the equivalent tendency for the conservative potential temperature
PRTHLS(:,:,:)= PRTHLS(:,:,:)  +    &
               PRHODJ(:,:,:)*(PTHLP(:,:,:)-PTHLM(:,:,:))/PTSTEP
!
!*       2.2  Partial Thermal Production
!
!  Conservative potential temperature flux : 
!
ZFLXZ(:,:,:)   = ZF                                                &
               + PIMPL * ZDFDDTDZ * DZM(PTHLP - PTHLM, KKA, KKU, KKL) / PDZZ 
!
ZFLXZ(:,:,KKA) = ZFLXZ(:,:,IKB) 
!  
  DO JK=IKTB+1,IKTE-1
   PWTH(:,:,JK)=0.5*(ZFLXZ(:,:,JK)+ZFLXZ(:,:,JK+KKL))
  END DO
  PWTH(:,:,IKB)=0.5*(ZFLXZ(:,:,IKB)+ZFLXZ(:,:,IKB+KKL))
  PWTH(:,:,KKA)=0.5*(ZFLXZ(:,:,KKA)+ZFLXZ(:,:,KKA+KKL))
  PWTH(:,:,IKE)=PWTH(:,:,IKE-KKL)

IF ( OTURB_FLX .AND. OCLOSE_OUT ) THEN
  ! stores the conservative potential temperature vertical flux
  YRECFM  ='THW_FLX'
  YCOMMENT='X_Y_Z_THW_FLX (K*M/S)'
  IGRID   = 4  
  ILENCH=LEN(YCOMMENT) 
  CALL FMWRIT(HFMFILE,YRECFM,HLUOUT,'XY',ZFLXZ,IGRID,ILENCH,YCOMMENT,IRESP)
END IF
!
! Contribution of the conservative temperature flux to the buoyancy flux
IF (KRR /= 0) THEN
  PTP(:,:,:)  =  PBETA * MZF(MZM(PETHETA, KKA, KKU, KKL) * ZFLXZ, KKA, KKU, KKL)
  PTP(:,:,IKB)=  PBETA(:,:,IKB) * PETHETA(:,:,IKB) *   &
                 0.5 * ( ZFLXZ (:,:,IKB) + ZFLXZ (:,:,IKB+KKL) )  
ELSE
  PTP(:,:,:)=  PBETA * MZF(ZFLXZ, KKA, KKU, KKL)
END IF
!
! Buoyancy flux at flux points
! 
PWTHV = MZM(PETHETA, KKA, KKU, KKL) * ZFLXZ
PWTHV(:,:,IKB) = PETHETA(:,:,IKB) * ZFLXZ(:,:,IKB)
!
!*       2.3  Partial vertical divergence of the < Rc w > flux
! Correction for qc and qi negative in AROME 
!IF ( KRRL >= 1 ) THEN
!  IF ( KRRI >= 1 ) THEN
!    PRRS(:,:,:,2) = PRRS(:,:,:,2) -                                        &
!                    DZF(MZM(PRHODJ*PATHETA*2.*PSRCM, KKA, KKU, KKL)*ZFLXZ/PDZZ, KKA, KKU, KKL)       &
!                    *(1.0-PFRAC_ICE(:,:,:))
!    PRRS(:,:,:,4) = PRRS(:,:,:,4) -                                        &
!                    DZF(MZM(PRHODJ*PATHETA*2.*PSRCM, KKA, KKU, KKL)*ZFLXZ/PDZZ, KKA, KKU, KKL)       &
!                    *PFRAC_ICE(:,:,:)
!  ELSE
!    PRRS(:,:,:,2) = PRRS(:,:,:,2) -                                        &
!                    DZF(MZM(PRHODJ*PATHETA*2.*PSRCM, KKA, KKU, KKL)*ZFLXZ/PDZZ, KKA, KKU, KKL) 
!  END IF
!END IF
!
!*       2.4  Storage in LES configuration
! 
IF (LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID(MZF(ZFLXZ, KKA, KKU, KKL), X_LES_SUBGRID_WThl ) 
  CALL LES_MEAN_SUBGRID(MZF(PWM*ZFLXZ, KKA, KKU, KKL), X_LES_RES_W_SBG_WThl )
  CALL LES_MEAN_SUBGRID(GZ_W_M(PWM,PDZZ, KKA, KKU, KKL)*MZF(ZFLXZ, KKA, KKU, KKL),&
      & X_LES_RES_ddxa_W_SBG_UaThl )
  CALL LES_MEAN_SUBGRID(MZF(PDTH_DZ*ZFLXZ, KKA, KKU, KKL), X_LES_RES_ddxa_Thl_SBG_UaThl )
  CALL LES_MEAN_SUBGRID(-XCTP*PSQRT_TKE/PLM*MZF(ZFLXZ, KKA, KKU, KKL), X_LES_SUBGRID_ThlPz ) 
  CALL LES_MEAN_SUBGRID(MZF(MZM(PETHETA, KKA, KKU, KKL)*ZFLXZ, KKA, KKU, KKL), X_LES_SUBGRID_WThv ) 
  IF (KRR>=1) THEN
    CALL LES_MEAN_SUBGRID(MZF(PDR_DZ*ZFLXZ, KKA, KKU, KKL), X_LES_RES_ddxa_Rt_SBG_UaThl )
  END IF
  !* diagnostic of mixing coefficient for heat
  ZA = DZM(PTHLP, KKA, KKU, KKL)
  WHERE (ZA==0.) ZA=1.E-6
  ZA = - ZFLXZ / ZA * PDZZ
  ZA(:,:,IKB) = XCSHF*PPHI3(:,:,IKB)*ZKEFF(:,:,IKB)
  ZA = MZF(ZA, KKA, KKU, KKL)
  ZA = MIN(MAX(ZA,-1000.),1000.)
  CALL LES_MEAN_SUBGRID( ZA, X_LES_SUBGRID_Kh   ) 
  !
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!*       2.5  New boundary layer depth for TOMs
! 
IF (HTOM=='TM06') CALL TM06_H(IKB,IKTB,IKTE,PTSTEP,PZZ,ZFLXZ,PBL_DEPTH)
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
 IF (LHARAT) THEN
  ZF      (:,:,:) = -ZKEFF*DZM(PRM(:,:,:,1), KKA, KKU, KKL)/PDZZ
  ZDFDDRDZ(:,:,:) = -ZKEFF
 ELSE
  ZF      (:,:,:) = -XCSHF*PPSI3*ZKEFF*DZM(PRM(:,:,:,1), KKA, KKU, KKL)/PDZZ
  ZDFDDRDZ(:,:,:) = -XCSHF*ZKEFF*D_PSI3DRDZ_O_DDRDZ(PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,GUSERV)
 ENDIF
  !
  ! Effect of 3rd order terms in temperature flux (at flux point)
  !
  ! d(w'2r')/dz
  IF (GFWR) THEN
    Z3RDMOMENT= M3_WR_W2R(KKA,KKU,KKL,PREDR1,PREDTH1,PD,ZKEFF,PTKEM)
  !
    ZF       = ZF       + Z3RDMOMENT * PFWR
    ZDFDDRDZ = ZDFDDRDZ + D_M3_WR_W2R_O_DDRDZ(KKA,KKU,KKL,PREDR1,PREDTH1,PD,&
     & PBLL_O_E,PEMOIST,ZKEFF,PTKEM) * PFWR
  END IF
  !
  ! d(w'r'2)/dz
  IF (GFR2) THEN
    Z3RDMOMENT= M3_WR_WR2(PREDR1,PREDTH1,PD,PBLL_O_E,PEMOIST)
  !
    ZF       = ZF       + Z3RDMOMENT * MZM(PFR2, KKA, KKU, KKL)
    ZDFDDRDZ = ZDFDDRDZ + D_M3_WR_WR2_O_DDRDZ(Z3RDMOMENT,PREDR1,&
     & PREDTH1,PD,PBLL_O_E,PEMOIST) * MZM(PFR2, KKA, KKU, KKL)
  END IF
  !
  ! d(w'2th')/dz
  IF (GFWTH) THEN
    ZF       = ZF       + M3_WR_W2TH(KKA,KKU,KKL,PREDR1,PREDTH1,PD,ZKEFF,&
     & PTKEM,PBLL_O_E,PETHETA,PDR_DZ) * PFWTH
    ZDFDDRDZ = ZDFDDRDZ + D_M3_WR_W2TH_O_DDRDZ(KKA,KKU,KKL,PREDR1,PREDTH1,& 
     & PD,ZKEFF,PTKEM,PBLL_O_E,PETHETA) * PFWTH
  END IF
  !
  ! d(w'th'2)/dz
  IF (GFTH2) THEN
    ZF       = ZF       + M3_WR_WTH2(KKA,KKU,KKL,PREDR1,PREDTH1,PD,ZKEFF,PTKEM,&
    & PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA,PDR_DZ) * MZM(PFTH2, KKA, KKU, KKL)
    ZDFDDRDZ = ZDFDDRDZ + D_M3_WR_WTH2_O_DDRDZ(KKA,KKU,KKL,PREDR1,PREDTH1,PD,&
     &ZKEFF,PTKEM,PSQRT_TKE,PBLL_O_E,PBETA,PLEPS,PETHETA) * MZM(PFTH2, KKA, KKU, KKL)
  END IF
  !
  ! d(w'th'r')/dz
  IF (GFTHR) THEN
    Z3RDMOMENT= M3_WR_WTHR(KKA,KKU,KKL,PREDTH1,PD,ZKEFF,PTKEM,PSQRT_TKE,PBETA,&
     & PLEPS,PETHETA)
  !
    ZF       = ZF       + Z3RDMOMENT * MZM(PFTHR, KKA, KKU, KKL)
    ZDFDDRDZ = ZDFDDRDZ + D_M3_WR_WTHR_O_DDRDZ(KKA,KKU,KKL,Z3RDMOMENT,PREDR1, &
     & PREDTH1,PD,PBLL_O_E,PEMOIST) * MZM(PFTHR, KKA, KKU, KKL)
  END IF
  !
  !* in 3DIM case, a part of the flux goes vertically, and another goes horizontally
  ! (in presence of slopes)
  !* in 1DIM case, the part of energy released in horizontal flux
  ! is taken into account in the vertical part
  !
  IF (HTURBDIM=='3DIM') THEN
    ZF(:,:,IKB) = ( PIMPL*PSFRP(:,:) + PEXPL*PSFRM(:,:) )       &
                         * PDIRCOSZW(:,:)                       &
                       * 0.5 * (1. + PRHODJ(:,:,KKA) / PRHODJ(:,:,IKB))
  ELSE
    ZF(:,:,IKB) = ( PIMPL*PSFRP(:,:) + PEXPL*PSFRM(:,:) )     &
                       / PDIRCOSZW(:,:)                       &
                       * 0.5 * (1. + PRHODJ(:,:,KKA) / PRHODJ(:,:,IKB))
  END IF
  !
  ! Compute the splitted conservative potential temperature at t+deltat
  CALL TRIDIAG_THERMO(KKA,KKU,KKL,PRM(:,:,:,1),ZF,ZDFDDRDZ,PTSTEP,PIMPL,&
                      PDZZ,PRHODJ,PRP)
  !
  ! Compute the equivalent tendency for the conservative mixing ratio
  PRRS(:,:,:,1) = PRRS(:,:,:,1) + PRHODJ(:,:,:) *     &
                  (PRP(:,:,:)-PRM(:,:,:,1))/PTSTEP
  !
  !*       3.2  Complete thermal production
  !
  ! cons. mixing ratio flux :
  !
  ZFLXZ(:,:,:)   = ZF                                                &
                 + PIMPL * ZDFDDRDZ * DZM(PRP - PRM(:,:,:,1), KKA, KKU, KKL) / PDZZ 
  !
  ZFLXZ(:,:,KKA) = ZFLXZ(:,:,IKB) 
  !
  DO JK=IKTB+1,IKTE-1
   PWRC(:,:,JK)=0.5*(ZFLXZ(:,:,JK)+ZFLXZ(:,:,JK+KKL))
  END DO
  PWRC(:,:,IKB)=0.5*(ZFLXZ(:,:,IKB)+ZFLXZ(:,:,IKB+KKL))
  PWRC(:,:,KKA)=0.5*(ZFLXZ(:,:,KKA)+ZFLXZ(:,:,KKA+KKL))
  PWRC(:,:,IKE)=PWRC(:,:,IKE-KKL)
  !
  !
  IF ( OTURB_FLX .AND. OCLOSE_OUT ) THEN
    ! stores the conservative mixing ratio vertical flux
    YRECFM  ='RCONSW_FLX'
    YCOMMENT='X_Y_Z_RCONSW_FLX (KG*M/S/KG)'
    IGRID   = 4  
    ILENCH=LEN(YCOMMENT) 
    CALL FMWRIT(HFMFILE,YRECFM,HLUOUT,'XY',ZFLXZ,IGRID,ILENCH,YCOMMENT,IRESP)
  END IF
  !
  ! Contribution of the conservative water flux to the Buoyancy flux
  ZA(:,:,:)   =  PBETA * MZF(MZM(PEMOIST, KKA, KKU, KKL) * ZFLXZ, KKA, KKU, KKL)
  ZA(:,:,IKB) =  PBETA(:,:,IKB) * PEMOIST(:,:,IKB) *   &
                 0.5 * ( ZFLXZ (:,:,IKB) + ZFLXZ (:,:,IKB+KKL) )
  PTP(:,:,:) = PTP(:,:,:) + ZA(:,:,:)
  !
  ! Buoyancy flux at flux points
  ! 
  PWTHV          = PWTHV          + MZM(PEMOIST, KKA, KKU, KKL) * ZFLXZ
  PWTHV(:,:,IKB) = PWTHV(:,:,IKB) + PEMOIST(:,:,IKB) * ZFLXZ(:,:,IKB)
!
!*       3.3  Complete vertical divergence of the < Rc w > flux
! Correction of qc and qi negative for AROME
!  IF ( KRRL >= 1 ) THEN
!    IF ( KRRI >= 1 ) THEN
!      PRRS(:,:,:,2) = PRRS(:,:,:,2) -                                        &
!                      DZF(MZM(PRHODJ*PAMOIST*2.*PSRCM, KKA, KKU, KKL)*ZFLXZ/PDZZ, KKA, KKU, KKL)       &
!                      *(1.0-PFRAC_ICE(:,:,:))
!      PRRS(:,:,:,4) = PRRS(:,:,:,4) -                                        &
!                      DZF(MZM(PRHODJ*PAMOIST*2.*PSRCM, KKA, KKU, KKL)*ZFLXZ/PDZZ, KKA, KKU, KKL)       &
!                      *PFRAC_ICE(:,:,:)
!    ELSE
!      PRRS(:,:,:,2) = PRRS(:,:,:,2) -                                        &
!                      DZF(MZM(PRHODJ*PAMOIST*2.*PSRCM, KKA, KKU, KKL)*ZFLXZ/PDZZ, KKA, KKU, KKL) 
!    END IF
!  END IF
!
!*       3.4  Storage in LES configuration
! 
  IF (LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID(MZF(ZFLXZ, KKA, KKU, KKL), X_LES_SUBGRID_WRt ) 
    CALL LES_MEAN_SUBGRID(MZF(PWM*ZFLXZ, KKA, KKU, KKL), X_LES_RES_W_SBG_WRt )
    CALL LES_MEAN_SUBGRID(GZ_W_M(PWM,PDZZ, KKA, KKU, KKL)*MZF(ZFLXZ, KKA, KKU, KKL),&
    & X_LES_RES_ddxa_W_SBG_UaRt )
    CALL LES_MEAN_SUBGRID(MZF(PDTH_DZ*ZFLXZ, KKA, KKU, KKL), X_LES_RES_ddxa_Thl_SBG_UaRt )
    CALL LES_MEAN_SUBGRID(MZF(PDR_DZ*ZFLXZ, KKA, KKU, KKL), X_LES_RES_ddxa_Rt_SBG_UaRt )
    CALL LES_MEAN_SUBGRID(MZF(MZM(PEMOIST, KKA, KKU, KKL)*ZFLXZ, KKA, KKU, KKL), X_LES_SUBGRID_WThv , .TRUE. ) 
    CALL LES_MEAN_SUBGRID(-XCTP*PSQRT_TKE/PLM*MZF(ZFLXZ, KKA, KKU, KKL), X_LES_SUBGRID_RtPz ) 
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
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
IF ( ((OTURB_FLX .AND. OCLOSE_OUT) .OR. LLES_CALL) .AND. (KRRL > 0) ) THEN
!  
! recover the Conservative potential temperature flux : 
! With LHARAT is true tke and length scales at half levels
! yet modify to use length scale and tke at half levels from vdfexcuhl
 IF (LHARAT) THEN
  ZA(:,:,:)   = DZM(PIMPL * PTHLP + PEXPL * PTHLM, KKA, KKU, KKL) / PDZZ *       &
                  (-PLM*PSQRT_TKE) 
 ELSE
  ZA(:,:,:)   = DZM(PIMPL * PTHLP + PEXPL * PTHLM, KKA, KKU, KKL) / PDZZ *       &
                  (-PPHI3*MZM(PLM*PSQRT_TKE, KKA, KKU, KKL)) * XCSHF 
 ENDIF
  ZA(:,:,IKB) = ( PIMPL*PSFTHP(:,:) + PEXPL*PSFTHM(:,:) ) &
               * PDIRCOSZW(:,:)
  !  
  ! compute <w Rc>
  ZFLXZ(:,:,:) = MZM(PAMOIST * 2.* PSRCM, KKA, KKU, KKL) * ZFLXZ(:,:,:) + &
                 MZM(PATHETA * 2.* PSRCM, KKA, KKU, KKL) * ZA(:,:,:)
  ZFLXZ(:,:,KKA) = ZFLXZ(:,:,IKB) 
  !                 
  ! store the liquid water mixing ratio vertical flux
  IF ( OTURB_FLX .AND. OCLOSE_OUT ) THEN
    YRECFM  ='RCW_FLX'
    YCOMMENT='X_Y_Z_RCW_FLX (KG*M/S/KG)'
    IGRID   = 4  
    ILENCH=LEN(YCOMMENT) 
    CALL FMWRIT(HFMFILE,YRECFM,HLUOUT,'XY',ZFLXZ,IGRID,ILENCH,YCOMMENT,IRESP)
  END IF
  !  
! and we store in LES configuration this subgrid flux <w'rc'>
!
  IF (LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MZF(ZFLXZ, KKA, KKU, KKL), X_LES_SUBGRID_WRc ) 
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
END IF !end of <w Rc>
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_FLUX',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER_THERMO_FLUX
