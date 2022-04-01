!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TURB_VER_THERMO_CORR
IMPLICIT NONE
CONTAINS      
SUBROUTINE TURB_VER_THERMO_CORR(D,CST,CSTURB,                       &
                      KRR,KRRL,KRRI,KSV,                            &
                      OTURB_FLX,HTURBDIM,HTOM,OHARAT,OCOMPUTE_SRC,  &
                      OCOUPLES,OLES_CALL,                           &
                      PIMPL,PEXPL,TPFILE,                           &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,           &
                      PRHODJ,PTHVREF,                               &
                      PSFTHM,PSFRM,PSFTHP,PSFRP,                    &
                      PWM,PTHLM,PRM,PSVM,                           &
                      PTKEM,PLM,PLEPS,                              &
                      PLOCPEXNM,PATHETA,PAMOIST,                    &
                      PBETA, PSQRT_TKE, PDTH_DZ, PDR_DZ, PRED2TH3,  &
                      PRED2R3, PRED2THR3, PBLL_O_E, PETHETA,        &
                      PEMOIST, PREDTH1, PREDR1, PPHI3, PPSI3, PD,   &
                      PFWTH,PFWR,PFTH2,PFR2,PFTHR,                  &
                      PTHLP,PRP,MFMOIST,PSIGS                       )
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
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!                     2012-02 (Y. Seity) add possibility to run with reversed 
!!                                              vertical levels
!!      Modifications  July 2015 (Wim de Rooy) OHARAT switch
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!--------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_CST, ONLY: CST_t
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_FIELD,          ONLY: TFIELDDATA, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS, ONLY: JPVEXT_TURB
!USE MODD_CONF
USE MODD_LES
!
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_GRADIENT_M
USE MODI_SHUMAN , ONLY : DZM, MZM, MZF
USE MODI_LES_MEAN_SUBGRID
!
USE MODE_IO_FIELD_WRITE, only: IO_Field_write
USE MODE_PRANDTL
!
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
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KSV           ! number of scalar var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
CHARACTER(len=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
CHARACTER(len=4),       INTENT(IN)   ::  HTOM         ! type of Third Order Moment
LOGICAL,                INTENT(IN)   ::  OHARAT
LOGICAL,                INTENT(IN)   ::  OLES_CALL    ! compute the LES diagnostics at current time-step
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and version 
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PDZZ, PDXX, PDYY, PDZX, PDZY
                                                      ! Metric coefficients
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  MFMOIST      ! moist mass flux dual scheme

REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PTHVREF      ! ref. state Virtual 
                                                      ! Potential Temperature 
!
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PSFTHM,PSFRM ! surface fluxes at time
!                                                     ! t - deltat 
!
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PSFTHP,PSFRP ! surface fluxes at time
!                                                     ! t + deltat 
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PWM 
! Vertical wind
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PTHLM 
! potential temperature at t-Delta t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN) ::  PRM          ! Mixing ratios 
                                                      ! at t-Delta t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN) ::  PSVM         ! Mixing ratios 
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PTKEM        ! TKE at time t
! In case OHARATU=TRUE, PLM already includes all stability corrections
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PLEPS        ! dissipative length   
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PBETA        ! buoyancy coefficient
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PSQRT_TKE    ! sqrt(e)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PDTH_DZ      ! d(th)/dz
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PDR_DZ       ! d(rt)/dz
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRED2TH3     ! 3D Redeslperger number R*2_th
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRED2R3      ! 3D Redeslperger number R*2_r
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRED2THR3    ! 3D Redeslperger number R*2_thr
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PBLL_O_E     ! beta * Lk * Leps / tke
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PETHETA      ! Coefficient for theta in theta_v computation
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PEMOIST      ! Coefficient for r in theta_v computation
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PREDTH1      ! 1D Redelsperger number for Th
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PREDR1       ! 1D Redelsperger number for r
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PPHI3        ! Prandtl number for temperature
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PPSI3        ! Prandtl number for vapor
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PD           ! Denominator in Prandtl numbers
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PFWTH        ! d(w'2th' )/dz (at flux point)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PFWR         ! d(w'2r'  )/dz (at flux point)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PFTH2        ! d(w'th'2 )/dz (at mass point)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PFR2         ! d(w'r'2  )/dz (at mass point)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PFTHR        ! d(w'th'r')/dz (at mass point)
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PTHLP      ! guess of thl at t+ deltat
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN)    :: PRP        ! guess of r at t+ deltat
!
REAL, DIMENSION(MERGE(D%NIT,0,OCOMPUTE_SRC),&
                MERGE(D%NJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)),   INTENT(OUT)  ::  PSIGS     ! Vert. part of Sigma_s at t
!
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  ::  &
       ZA,       & ! work variable for wrc
       ZFLXZ,    & ! vertical flux of the treated variable
       ZSOURCE,  & ! source of evolution for the treated variable
       ZKEFF,    & ! effectif diffusion coeff = LT * SQRT( TKE )
       ZF,       & ! Flux in dTh/dt =-dF/dz (evaluated at t-1)(or rt instead of Th)
       ZDFDDTDZ, & ! dF/d(dTh/dz)
       ZDFDDRDZ, & ! dF/d(dr/dz)
       Z3RDMOMENT, & ! 3 order term in flux or variance equation
! Estimate of full level length and dissipation length scale in case OHARATU
       PLMF,     & ! estimate full level length scale from half levels (sub optimal)
       PLEPSF,   & ! estimate full level diss length scale from half levels (sub optimal)
       ZWORK1,ZWORK2, &
       ZWORK3,ZWORK4 ! working var. for shuman operators (array syntax)

INTEGER             :: IKB,IKE      ! I index values for the Beginning and End
INTEGER             :: IKU  ! array sizes
INTEGER             :: JI, JJ, JK ! loop indexes 

REAL, DIMENSION(D%NIT,D%NJT,MIN(D%NKA+JPVEXT_TURB*D%NKL,D%NKA+JPVEXT_TURB*D%NKL+2*D%NKL):&
                            MAX(D%NKA+JPVEXT_TURB*D%NKL,D%NKA+JPVEXT_TURB*D%NKL+2*D%NKL))&
                    :: ZCOEFF
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground, defined in
                                    ! mass points of the domain in the 3 direct.
!
REAL :: ZTIME1, ZTIME2
!
LOGICAL :: GUSERV   ! flag to use water
LOGICAL :: GFTH2    ! flag to use w'th'2
LOGICAL :: GFWTH    ! flag to use w'2th'
LOGICAL :: GFR2     ! flag to use w'r'2
LOGICAL :: GFWR     ! flag to use w'2r'
LOGICAL :: GFTHR    ! flag to use w'th'r'
TYPE(TFIELDDATA) :: TZFIELD
!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_CORR',0,ZHOOK_HANDLE)
!
IKB=D%NKB
IKE=D%NKE
!
GUSERV = (KRR/=0)
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
ZCOEFF(:,:,IKB+2*D%NKL)= - PDZZ(:,:,IKB+D%NKL) /      &
       ( (PDZZ(:,:,IKB+2*D%NKL)+PDZZ(:,:,IKB+D%NKL)) * PDZZ(:,:,IKB+2*D%NKL) )
ZCOEFF(:,:,IKB+D%NKL)=   (PDZZ(:,:,IKB+2*D%NKL)+PDZZ(:,:,IKB+D%NKL)) /      &
       ( PDZZ(:,:,IKB+D%NKL) * PDZZ(:,:,IKB+2*D%NKL) )
ZCOEFF(:,:,IKB)= - (PDZZ(:,:,IKB+2*D%NKL)+2.*PDZZ(:,:,IKB+D%NKL)) /      &
       ( (PDZZ(:,:,IKB+2*D%NKL)+PDZZ(:,:,IKB+D%NKL)) * PDZZ(:,:,IKB+D%NKL) )
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
!
!
IF (OHARAT) THEN
PLMF=MZF(PLM, D%NKA, D%NKU, D%NKL)
PLEPSF=PLMF
!  function MZF produces -999 for level IKU (82 for 80 levels)
!  so put these to normal value as this level (82) is indeed calculated
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
PLMF(:,:,D%NKT)=0.001
PLEPSF(:,:,D%NKT)=0.001
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
!$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
ZKEFF(:,:,:) = PLM(:,:,:) * SQRT(PTKEM(:,:,:)) + 50*MFMOIST(:,:,:)
!$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
ELSE
ZKEFF(:,:,:) = MZM(PLM(:,:,:) * SQRT(PTKEM(:,:,:)), D%NKA, D%NKU, D%NKL)
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
!
!*       4.   TURBULENT CORRELATIONS : <THl THl>, <THl Rnp>, <Rnp Rnp>
!             --------------------------------------------------------
!
!
!*       4.2  <THl THl> 
!
! Compute the turbulent variance F and F' at time t-dt.
!
IF (OHARAT) THEN
  ZWORK1 = MZF(PDTH_DZ**2, D%NKA, D%NKU, D%NKL)
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  ZF(:,:,:) = PLMF(:,:,:)*PLEPSF(:,:,:)*ZWORK1(:,:,:)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
ELSE
  ZWORK1 = MZF(PPHI3*PDTH_DZ**2, D%NKA, D%NKU, D%NKL)
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  ZF(:,:,:) = CSTURB%XCTV*PLM(:,:,:)*PLEPS(:,:,:)*ZWORK1(:,:,:)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
ENDIF
  ZDFDDTDZ(:,:,:) = 0.     ! this term, because of discretization, is treated separately
  !
  ! Effect of 3rd order terms in temperature flux (at mass point)
  !
  ! d(w'th'2)/dz
  IF (GFTH2) THEN
    ZWORK1 = M3_TH2_WTH2(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE)
    ZWORK2 = D_M3_TH2_WTH2_O_DDTDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,&
     & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA)
    !
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    ZF(:,:,:)       = ZF(:,:,:) + ZWORK1(:,:,:) * PFTH2(:,:,:)
    ZDFDDTDZ(:,:,:) = ZDFDDTDZ(:,:,:) + ZWORK2(:,:,:) * PFTH2(:,:,:)
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  END IF
  !
  ! d(w'2th')/dz
  IF (GFWTH) THEN
    ZWORK1 = M3_TH2_W2TH(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,PD,PDTH_DZ,&
     & PLM,PLEPS,PTKEM)
    ZWORK2 = MZF(PFWTH, D%NKA, D%NKU, D%NKL)
    ZWORK3 = D_M3_TH2_W2TH_O_DDTDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,PD,&
     & PLM,PLEPS,PTKEM,GUSERV)
    !
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    ZF(:,:,:)  = ZF(:,:,:) + ZWORK1(:,:,:) * ZWORK2(:,:,:,)
    ZDFDDTDZ(:,:,:) = ZDFDDTDZ(:,:,:) + ZWORK3(:,:,:) * ZWORK2(:,:,:)
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  END IF
  !
  IF (KRR/=0) THEN
    ! d(w'r'2)/dz
    IF (GFR2) THEN
      ZWORK1 = M3_TH2_WR2(D,CSTURB,D%NKA,D%NKU,D%NKL,PD,PLEPS,PSQRT_TKE,PBLL_O_E,&
       & PEMOIST,PDTH_DZ)
      ZWORK2 = D_M3_TH2_WR2_O_DDTDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,PD,&
       & PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTH_DZ)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZF(:,:,:) = ZF(:,:,:) + ZWORK1(:,:,:) * PFR2(:,:,:)
      ZDFDDTDZ = ZDFDDTDZ + ZWORK2(:,:,:) * PFR2(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF
    !
    ! d(w'2r')/dz
    IF (GFWR) THEN
      ZWORK1 = M3_TH2_W2R(D,CSTURB,D%NKA,D%NKU,D%NKL,PD,PLM,PLEPS,PTKEM,PBLL_O_E,&
       & PEMOIST,PDTH_DZ)
      ZWORK2 = MZF(PFWR, D%NKA, D%NKU, D%NKL)
      ZWORK3 = D_M3_TH2_W2R_O_DDTDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,PD,&
       & PLM,PLEPS,PTKEM,PBLL_O_E,PEMOIST,PDTH_DZ)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZF(:,:,:) = ZF(:,:,:) + ZWORK1(:,:,:) * ZWORK2(:,:,:)
      ZDFDDTDZ(:,:,:) = ZDFDDTDZ(:,:,:) + ZWORK3(:,:,:) * ZWORK1(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF
    !
    ! d(w'th'r')/dz
    IF (GFTHR) THEN
      ZWORK1 = M3_TH2_WTHR(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PD,PLEPS,PSQRT_TKE,&
       & PBLL_O_E,PEMOIST,PDTH_DZ)
      ZWORK2 = D_M3_TH2_WTHR_O_DDTDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,&
       & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTH_DZ)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZF(:,:,:) = ZF(:,:,:) + ZWORK1(:,:,:) * PFTHR(:,:,:)
      ZDFDDTDZ(:,:,:) = ZDFDDTDZ(:,:,:) + ZWORK2(:,:,:) * PFTHR(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF

  END IF
  !
  ZWORK1 = MZF(DZM(PTHLP - PTHLM, D%NKA, D%NKU, D%NKL) / PDZZ, D%NKA, D%NKU, D%NKL)
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  ZFLXZ(:,:,:)   = ZF(:,:,:) + PIMPL * ZDFDDTDZ(:,:,:) * ZWORK1(:,:,:)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  !
  ! special case near the ground ( uncentred gradient )
  IF (OHARAT) THEN
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
    ZFLXZ(:,:,IKB) =  PLMF(:,:,IKB)   &
     * PLEPSF(:,:,IKB)                                         &
     *( PEXPL *                                                  &
     ( ZCOEFF(:,:,IKB+2*D%NKL)*PTHLM(:,:,IKB+2*D%NKL)             &
      +ZCOEFF(:,:,IKB+D%NKL  )*PTHLM(:,:,IKB+D%NKL  )             & 
      +ZCOEFF(:,:,IKB      )*PTHLM(:,:,IKB  )   )**2          &
     +PIMPL *                                                  &
     ( ZCOEFF(:,:,IKB+2*D%NKL)*PTHLP(:,:,IKB+2*D%NKL)             &
      +ZCOEFF(:,:,IKB+D%NKL  )*PTHLP(:,:,IKB+D%NKL  )             &
      +ZCOEFF(:,:,IKB      )*PTHLP(:,:,IKB  )   )**2          &
    ) 
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
   ELSE
     !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
     ZFLXZ(:,:,IKB) = CSTURB%XCTV * PPHI3(:,:,IKB+D%NKL) * PLM(:,:,IKB)   &
     * PLEPS(:,:,IKB)                                         &
     *( PEXPL *                                                  &
     ( ZCOEFF(:,:,IKB+2*D%NKL)*PTHLM(:,:,IKB+2*D%NKL)             &
      +ZCOEFF(:,:,IKB+D%NKL  )*PTHLM(:,:,IKB+D%NKL  )             & 
      +ZCOEFF(:,:,IKB      )*PTHLM(:,:,IKB  )   )**2          &
     +PIMPL *                                                  &
     ( ZCOEFF(:,:,IKB+2*D%NKL)*PTHLP(:,:,IKB+2*D%NKL)             &
      +ZCOEFF(:,:,IKB+D%NKL  )*PTHLP(:,:,IKB+D%NKL  )             &
      +ZCOEFF(:,:,IKB      )*PTHLP(:,:,IKB  )   )**2          &
     )
     !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT) 
   ENDIF
  !
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT) 
  ZFLXZ(:,:,D%NKA) = ZFLXZ(:,:,IKB)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT) 
  !
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  ZFLXZ(:,:,:) = MAX(0., ZFLXZ(:,:,:))
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  !
  IF (KRRL > 0)  THEN
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    PSIGS(:,:,:) = ZFLXZ(:,:,:) * PATHETA(:,:,:)**2
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  END IF
  !
  !
  ! stores <THl THl>  
  IF ( OTURB_FLX .AND. TPFILE%LOPENED ) THEN
    TZFIELD%CMNHNAME   = 'THL_VVAR'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THL_VVAR'
    TZFIELD%CUNITS     = 'K2'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_THL_VVAR'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
  END IF
!
! and we store in LES configuration
!
  IF (OLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID(ZFLXZ, X_LES_SUBGRID_Thl2 ) 
    CALL LES_MEAN_SUBGRID(MZF(PWM, D%NKA, D%NKU, D%NKL)*ZFLXZ, X_LES_RES_W_SBG_Thl2 )
    CALL LES_MEAN_SUBGRID(-2.*CSTURB%XCTD*PSQRT_TKE*ZFLXZ/PLEPS, X_LES_SUBGRID_DISS_Thl2 ) 
    CALL LES_MEAN_SUBGRID(PETHETA*ZFLXZ, X_LES_SUBGRID_ThlThv ) 
    CALL LES_MEAN_SUBGRID(-CSTURB%XA3*PBETA*PETHETA*ZFLXZ, X_LES_SUBGRID_ThlPz, .TRUE. ) 
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
  IF ( KRR /= 0 ) THEN
!
!*       4.3  <THl Rnp>    
!
!
    ! Compute the turbulent variance F and F' at time t-dt.
  IF (OHARAT) THEN
    ZWORK1 = MZF(PDTH_DZ*PDR_DZ, D%NKA, D%NKU, D%NKL)
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    ZF(:,:,:) = PLMF(:,:,:)*PLEPSF(:,:,:)*ZWORK1(:,:,:)
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  ELSE
    ZWORK1 = MZF(0.5*(PPHI3+PPSI3)*PDTH_DZ*PDR_DZ, D%NKA, D%NKU, D%NKL)
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    ZF(:,:,:) = CSTURB%XCTV*PLM(:,:,:)*PLEPS(:,:,:)*ZWORK1(:,:,:)
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  ENDIF
    ZDFDDTDZ(:,:,:) = 0.     ! this term, because of discretization, is treated separately
    ZDFDDRDZ(:,:,:) = 0.     ! this term, because of discretization, is treated separately
    !
    ! Effect of 3rd order terms in temperature flux (at mass point)
    !
    ! d(w'th'2)/dz
    IF (GFTH2) THEN
      ZWORK1 = M3_THR_WTH2(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PD,PLEPS,PSQRT_TKE,&
       & PBLL_O_E,PETHETA,PDR_DZ)
      ZWORK2 = D_M3_THR_WTH2_O_DDTDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,&
       & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ)
      ZWORK3 = D_M3_THR_WTH2_O_DDRDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,&
       & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZF(:,:,:) = ZF(:,:,:) + ZWORK1(:,:,:) * PFTH2(:,:,:)
      ZDFDDTDZ(:,:,:) = ZDFDDTDZ(:,:,:) + ZWORK2(:,:,:) * PFTH2(:,:,:)
      ZDFDDRDZ(:,:,:) = ZDFDDRDZ(:,:,:) + ZWORK3(:,:,:) * PFTH2(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF
    !
    ! d(w'2th')/dz
    IF (GFWTH) THEN
      ZWORK1 = MZF(PFWTH, D%NKA, D%NKU, D%NKL)
      ZWORK2 = M3_THR_W2TH(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PD,PLM,PLEPS,PTKEM,&
       & PDR_DZ)
      ZWORK3 = D_M3_THR_W2TH_O_DDTDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,&
       & PD,PLM,PLEPS,PTKEM,PBLL_O_E,PDR_DZ,PETHETA)
      ZWORK4 = D_M3_THR_W2TH_O_DDRDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,&
       & PD,PLM,PLEPS,PTKEM)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZF(:,:,:) = ZF(:,:,:) + ZWORK2(:,:,:) * ZWORK1(:,:,:)
      ZDFDDTDZ(:,:,:) = ZDFDDTDZ(:,:,:) + ZWORK3(:,:,:) * ZWORK1(:,:,:)
      ZDFDDRDZ(:,:,:) = ZDFDDRDZ(:,:,:) + ZWORK4(:,:,:) * ZWORK1(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF
    !
    ! d(w'r'2)/dz
    IF (GFR2) THEN
      ZWORK1 = M3_THR_WR2(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PD,PLEPS,PSQRT_TKE,&
       & PBLL_O_E,PEMOIST,PDTH_DZ)
      ZWORK2 = D_M3_THR_WR2_O_DDTDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PREDTH1,PD,&
       & PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
      ZWORK3 = D_M3_THR_WR2_O_DDRDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PREDTH1,PD,&
       & PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTH_DZ)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZF(:,:,:) = ZF(:,:,:) + ZWORK1(:,:,:) * PFR2(:,:,:)
      ZDFDDTDZ(:,:,:) = ZDFDDTDZ(:,:,:) + ZWORK2(:,:,:) * PFR2(:,:,:)
      ZDFDDRDZ(:,:,:) = ZDFDDRDZ(:,:,:) + ZWORK3(:,:,:) * PFR2(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF
    !
      ! d(w'2r')/dz
    IF (GFWR) THEN
      ZWORK1 = MZF(PFWR, D%NKA, D%NKU, D%NKL)
      ZWORK2 = M3_THR_W2R(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PD,PLM,PLEPS,PTKEM,&
      & PDTH_DZ)
      ZWORK3 = D_M3_THR_W2R_O_DDTDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PREDTH1,PD,&
      & PLM,PLEPS,PTKEM)
      ZWORK4 = D_M3_THR_W2R_O_DDRDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PREDTH1,PD,&
      & PLM,PLEPS,PTKEM,PBLL_O_E,PDTH_DZ,PEMOIST)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZF(:,:,:) = ZF(:,:,:) + ZWORK2(:,:,:) * ZWORK1(:,:,:)
      ZDFDDTDZ(:,:,:) = ZDFDDTDZ(:,:,:) + ZWORK3(:,:,:) * ZWORK1(:,:,:)
      ZDFDDRDZ(:,:,:) = ZDFDDRDZ(:,:,:) + ZWORK4(:,:,:) * ZWORK1(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF
    !
    ! d(w'th'r')/dz
    IF (GFTHR) THEN
      ZWORK1 = M3_THR_WTHR(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,PD,PLEPS,&
      & PSQRT_TKE)
      ZWORK2 = D_M3_THR_WTHR_O_DDTDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PREDR1,&
      & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA)
      ZWORK3 = D_M3_THR_WTHR_O_DDRDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PREDTH1,&
      & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZF(:,:,:) = ZF(:,:,:) + ZWORK1(:,:,:) * PFTHR(:,:,:)
      ZDFDDTDZ(:,:,:) = ZDFDDTDZ(:,:,:) + ZWORK2(:,:,:) * PFTHR(:,:,:)
      ZDFDDRDZ(:,:,:) = ZDFDDRDZ(:,:,:) + ZWORK3(:,:,:) * PFTHR(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF
    !
    ZWORK1 = MZF(DZM(PTHLP - PTHLM(:,:,:), D%NKA, D%NKU, D%NKL) / PDZZ, D%NKA, D%NKU, D%NKL)
    ZWORK2 = MZF(DZM(PRP   - PRM(:,:,:,1), D%NKA, D%NKU, D%NKL) / PDZZ, D%NKA, D%NKU, D%NKL)
    IF (OHARAT) THEN
      ZWORK3 = MZF(( 2.  & 
                 ) *PDR_DZ  *DZM(PTHLP - PTHLM, D%NKA, D%NKU, D%NKL    ) / PDZZ               &
                +( 2.                                                    &
                 ) *PDTH_DZ *DZM(PRP   - PRM(:,:,:,1), D%NKA, D%NKU, D%NKL) / PDZZ               &
               , D%NKA, D%NKU, D%NKL)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZFLXZ(:,:,:)   = ZF(:,:,:)                                           &
        + PIMPL * PLMF(:,:,:)*PLEPSF(:,:,:)*0.5                          &
          * ZWORK3(:,:,:)                                                &
        + PIMPL * ZDFDDTDZ(:,:,:) * ZWORK1(:,:,:)         &
        + PIMPL * ZDFDDRDZ(:,:,:) * ZWORK2(:,:,:)
       !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    ELSE
      ZWORK3 = MZF(( D_PHI3DTDZ_O_DDTDZ(D,CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,GUSERV) & ! d(phi3*dthdz)/ddthdz term
                  +D_PSI3DTDZ_O_DDTDZ(D,CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,GUSERV) & ! d(psi3*dthdz)/ddthdz term
                 ) *PDR_DZ  *DZM(PTHLP - PTHLM, D%NKA, D%NKU, D%NKL) / PDZZ               &
                +( D_PHI3DRDZ_O_DDRDZ(D,CSTURB,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,GUSERV) & ! d(phi3*drdz )/ddrdz term
                  +D_PSI3DRDZ_O_DDRDZ(D,CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,GUSERV) & ! d(psi3*drdz )/ddrdz term
                 ) *PDTH_DZ *DZM(PRP - PRM(:,:,:,1), D%NKA, D%NKU, D%NKL) / PDZZ               &
               , D%NKA, D%NKU, D%NKL)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZFLXZ(:,:,:)   = ZF(:,:,:)                                                     &
        + PIMPL * CSTURB%XCTV*PLM(:,:,:)*PLEPS(:,:,:)*0.5                                        &
          * ZWORK3(:,:,:)                                                           &
        + PIMPL * ZDFDDTDZ(:,:,:) * ZWORK1(:,:,:)         &
        + PIMPL * ZDFDDRDZ(:,:,:) * ZWORK2(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    ENDIF
    !
    ! special case near the ground ( uncentred gradient )
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
    IF (OHARAT) THEN
    ZFLXZ(:,:,IKB) =                                            & 
    (1. )   &
    *( PEXPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*D%NKL)*PTHLM(:,:,IKB+2*D%NKL)             &
        +ZCOEFF(:,:,IKB+D%NKL  )*PTHLM(:,:,IKB+D%NKL  )             & 
        +ZCOEFF(:,:,IKB      )*PTHLM(:,:,IKB      ))            &
      *( ZCOEFF(:,:,IKB+2*D%NKL)*PRM(:,:,IKB+2*D%NKL,1)             &
        +ZCOEFF(:,:,IKB+D%NKL  )*PRM(:,:,IKB+D%NKL,1  )             & 
        +ZCOEFF(:,:,IKB      )*PRM(:,:,IKB  ,1    ))            &
      +PIMPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*D%NKL)*PTHLP(:,:,IKB+2*D%NKL)             &
        +ZCOEFF(:,:,IKB+D%NKL  )*PTHLP(:,:,IKB+D%NKL  )             &
        +ZCOEFF(:,:,IKB      )*PTHLP(:,:,IKB      ))            &
      *( ZCOEFF(:,:,IKB+2*D%NKL)*PRP(:,:,IKB+2*D%NKL  )             &
        +ZCOEFF(:,:,IKB+D%NKL  )*PRP(:,:,IKB+D%NKL    )             & 
        +ZCOEFF(:,:,IKB      )*PRP(:,:,IKB        ))            &
     ) 
    ELSE 
    ZFLXZ(:,:,IKB) =                                            & 
    (CSTURB%XCHT1 * PPHI3(:,:,IKB+D%NKL) + CSTURB%XCHT2 * PPSI3(:,:,IKB+D%NKL))   &
    *( PEXPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*D%NKL)*PTHLM(:,:,IKB+2*D%NKL)             &
        +ZCOEFF(:,:,IKB+D%NKL  )*PTHLM(:,:,IKB+D%NKL  )             & 
        +ZCOEFF(:,:,IKB      )*PTHLM(:,:,IKB      ))            &
      *( ZCOEFF(:,:,IKB+2*D%NKL)*PRM(:,:,IKB+2*D%NKL,1)             &
        +ZCOEFF(:,:,IKB+D%NKL  )*PRM(:,:,IKB+D%NKL,1  )             & 
        +ZCOEFF(:,:,IKB      )*PRM(:,:,IKB  ,1    ))            &
      +PIMPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*D%NKL)*PTHLP(:,:,IKB+2*D%NKL)             &
        +ZCOEFF(:,:,IKB+D%NKL  )*PTHLP(:,:,IKB+D%NKL  )             &
        +ZCOEFF(:,:,IKB      )*PTHLP(:,:,IKB      ))            &
      *( ZCOEFF(:,:,IKB+2*D%NKL)*PRP(:,:,IKB+2*D%NKL  )             &
        +ZCOEFF(:,:,IKB+D%NKL  )*PRP(:,:,IKB+D%NKL    )             & 
        +ZCOEFF(:,:,IKB      )*PRP(:,:,IKB        ))            &
     ) 
    ENDIF
    !    
    ZFLXZ(:,:,D%NKA) = ZFLXZ(:,:,IKB)
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
    !
      IF ( KRRL > 0 ) THEN
!
!   
!  NB PATHETA is -b in Chaboureau Bechtold 2002 which explains the + sign here
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      PSIGS(:,:,:) = PSIGS(:,:,:) +     &
                     2. * PATHETA(:,:,:) * PAMOIST(:,:,:) * ZFLXZ(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF
    ! stores <THl Rnp>   
    IF ( OTURB_FLX .AND. TPFILE%LOPENED ) THEN
      TZFIELD%CMNHNAME   = 'THLRCONS_VCOR'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = 'THLRCONS_VCOR'
      TZFIELD%CUNITS     = 'K kg kg-1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CCOMMENT   = 'X_Y_Z_THLRCONS_VCOR'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
    END IF
!
! and we store in LES configuration
!
IF (OLES_CALL) THEN
      CALL SECOND_MNH(ZTIME1)
      CALL LES_MEAN_SUBGRID(ZFLXZ, X_LES_SUBGRID_THlRt ) 
      CALL LES_MEAN_SUBGRID(MZF(PWM, D%NKA, D%NKU, D%NKL)*ZFLXZ, X_LES_RES_W_SBG_ThlRt )
      CALL LES_MEAN_SUBGRID(-2.*CSTURB%XCTD*PSQRT_TKE*ZFLXZ/PLEPS, X_LES_SUBGRID_DISS_ThlRt ) 
      CALL LES_MEAN_SUBGRID(PETHETA*ZFLXZ, X_LES_SUBGRID_RtThv ) 
      CALL LES_MEAN_SUBGRID(-CSTURB%XA3*PBETA*PETHETA*ZFLXZ, X_LES_SUBGRID_RtPz, .TRUE. ) 
      CALL LES_MEAN_SUBGRID(PEMOIST*ZFLXZ, X_LES_SUBGRID_ThlThv , .TRUE. ) 
      CALL LES_MEAN_SUBGRID(-CSTURB%XA3*PBETA*PEMOIST*ZFLXZ, X_LES_SUBGRID_ThlPz, .TRUE. ) 
      CALL SECOND_MNH(ZTIME2)
      XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
! 
!
!*       4.4  <Rnp Rnp>
!
!
    ! Compute the turbulent variance F and F' at time t-dt.
IF (OHARAT) THEN
  ZWORK1 = MZF(PDR_DZ**2, D%NKA, D%NKU, D%NKL)
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  ZF(:,:,:) = PLMF(:,:,:)*PLEPSF(:,:,:)*ZWORK1(:,:,:)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
ELSE
  ZWORK1 = MZF(PPSI3*PDR_DZ**2, D%NKA, D%NKU, D%NKL)
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  ZF(:,:,:) = CSTURB%XCTV*PLM(:,:,:)*PLEPS(:,:,:)*ZWORK1(:,:,:)
  !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
ENDIF
    ZDFDDRDZ(:,:,:) = 0.     ! this term, because of discretization, is treated separately
    !
    ! Effect of 3rd order terms in temperature flux (at mass point)
    !
    ! d(w'r'2)/dz
    IF (GFR2) THEN
      ZWORK1 = M3_R2_WR2(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PREDTH1,PD,PLEPS,&
      & PSQRT_TKE)
      ZWORK2 = D_M3_R2_WR2_O_DDRDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PREDTH1,&
      & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZF(:,:,:) = ZF(:,:,:) + ZWORK1(:,:,:) * PFR2(:,:,:)
      ZDFDDRDZ(:,:,:) = ZDFDDRDZ(:,:,:) + ZWORK2(:,:,:) * PFR2(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF
    !
    ! d(w'2r')/dz
    IF (GFWR) THEN
      ZWORK1 = MZF(PFWR, D%NKA, D%NKU, D%NKL)
      ZWORK2 = M3_R2_W2R(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PREDTH1,PD,PDR_DZ,&
      & PLM,PLEPS,PTKEM)
      ZWORK3 = D_M3_R2_W2R_O_DDRDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PREDTH1,&
      & PD,PLM,PLEPS,PTKEM,GUSERV)
      !
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      ZF(:,:,:) = ZF(:,:,:) + ZWORK2(:,:,:) * ZWORK1(:,:,:)
      ZDFDDRDZ(:,:,:) = ZDFDDRDZ(:,:,:) + ZWORK3(:,:,:) * ZWORK1(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF
    !
    IF (KRR/=0) THEN
      ! d(w'r'2)/dz
      IF (GFTH2) THEN
        ZWORK1 = M3_R2_WTH2(D,CSTURB,D%NKA,D%NKU,D%NKL,PD,PLEPS,PSQRT_TKE,&
          & PBLL_O_E,PETHETA,PDR_DZ)
        ZWORK2 = D_M3_R2_WTH2_O_DDRDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,&
          & PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ)
        !
        !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
        ZF(:,:,:) = ZF(:,:,:) + ZWORK1(:,:,:) * PFTH2(:,:,:) 
        ZDFDDRDZ(:,:,:) = ZDFDDRDZ(:,:,:) + ZWORK2(:,:,:) * PFTH2(:,:,:)
        !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      END IF
      !
      ! d(w'2r')/dz
      IF (GFWTH) THEN
        ZWORK1 = MZF(PFWTH, D%NKA, D%NKU, D%NKL)
        ZWORK2 = M3_R2_W2TH(D,CSTURB,D%NKA,D%NKU,D%NKL,PD,PLM,PLEPS,PTKEM,&
          & PBLL_O_E,PETHETA,PDR_DZ)
        ZWORK3 = D_M3_R2_W2TH_O_DDRDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PREDTH1,&
          & PD,PLM,PLEPS,PTKEM,PBLL_O_E,PETHETA,PDR_DZ)
        !
        !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
        ZF(:,:,:) = ZF(:,:,:) + ZWORK2(:,:,:) * ZWORK1(:,:,:)
        ZDFDDRDZ(:,:,:) = ZDFDDRDZ(:,:,:) + ZWORK3(:,:,:) * ZWORK1(:,:,:)
        !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      END IF
      !
      ! d(w'th'r')/dz
      IF (GFTHR) THEN
        ZWORK1 = M3_R2_WTHR(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDTH1,PD,PLEPS,&
          & PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ)
        ZWORK2 = D_M3_R2_WTHR_O_DDRDZ(D,CSTURB,D%NKA,D%NKU,D%NKL,PREDR1,PREDTH1,&
          & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ)
        !
        !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
        ZF(:,:,:) = ZF(:,:,:) + ZWORK1(:,:,:) * PFTHR(:,:,:) 
        ZDFDDRDZ(:,:,:) = ZDFDDRDZ(:,:,:) + ZWORK2(:,:,:) * PFTHR(:,:,:)
        !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      END IF
  
    END IF
    !
  IF (OHARAT) THEN
    ZWORK1 = MZF(2.*PDR_DZ*    &
                 DZM(PRP - PRM(:,:,:,1), D%NKA, D%NKU, D%NKL) / PDZZ, D%NKA, D%NKU, D%NKL)
    ZWORK2 = MZF(DZM(PRP - PRM(:,:,:,1), D%NKA, D%NKU, D%NKL) / PDZZ, D%NKA, D%NKU, D%NKL)
    !
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    ZFLXZ(:,:,:)   = ZF(:,:,:)                                                             &
          + PIMPL * PLMF(:,:,:) *PLEPSF(:,:,:)                                                   &
            * ZWORK1(:,:,:) &
          + PIMPL * ZDFDDRDZ(:,:,:) * ZWORK2(:,:,:)
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
   ELSE
    ZWORK1 = MZF(D_PSI3DRDZ2_O_DDRDZ(D,CSTURB,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,PDR_DZ,HTURBDIM,GUSERV)    &
                 *DZM(PRP - PRM(:,:,:,1), D%NKA, D%NKU, D%NKL) / PDZZ, D%NKA, D%NKU, D%NKL)
    ZWORK2 = MZF(DZM(PRP - PRM(:,:,:,1), D%NKA, D%NKU, D%NKL) / PDZZ, D%NKA, D%NKU, D%NKL)
    !
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    ZFLXZ(:,:,:)   = ZF(:,:,:)                                                             &
          + PIMPL * CSTURB%XCTV*PLM(:,:,:) *PLEPS(:,:,:)                                                 &
            * ZWORK1(:,:,:) &
          + PIMPL * ZDFDDRDZ(:,:,:) * ZWORK2(:,:,:)
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  ENDIF
    !
    ! special case near the ground ( uncentred gradient )
  !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
  IF (OHARAT) THEN
    ZFLXZ(:,:,IKB) =  PLMF(:,:,IKB)   &
        * PLEPSF(:,:,IKB)                                        &
    *( PEXPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*D%NKL)*PRM(:,:,IKB+2*D%NKL,1)             &
        +ZCOEFF(:,:,IKB+D%NKL  )*PRM(:,:,IKB+D%NKL,1  )             & 
        +ZCOEFF(:,:,IKB      )*PRM(:,:,IKB  ,1    ))**2         &
      +PIMPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*D%NKL)*PRP(:,:,IKB+2*D%NKL)               &
        +ZCOEFF(:,:,IKB+D%NKL  )*PRP(:,:,IKB+D%NKL  )               &
        +ZCOEFF(:,:,IKB      )*PRP(:,:,IKB      ))**2           &
     ) 
   ELSE
    ZFLXZ(:,:,IKB) = CSTURB%XCHV * PPSI3(:,:,IKB+D%NKL) * PLM(:,:,IKB)   &
        * PLEPS(:,:,IKB)                                        &
    *( PEXPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*D%NKL)*PRM(:,:,IKB+2*D%NKL,1)             &
        +ZCOEFF(:,:,IKB+D%NKL  )*PRM(:,:,IKB+D%NKL,1  )             & 
        +ZCOEFF(:,:,IKB      )*PRM(:,:,IKB  ,1    ))**2         &
      +PIMPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*D%NKL)*PRP(:,:,IKB+2*D%NKL)               &
        +ZCOEFF(:,:,IKB+D%NKL  )*PRP(:,:,IKB+D%NKL  )               &
        +ZCOEFF(:,:,IKB      )*PRP(:,:,IKB      ))**2           &
     ) 
  ENDIF
    !
    ZFLXZ(:,:,D%NKA) = ZFLXZ(:,:,IKB)
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
    !
    IF ( KRRL > 0 ) THEN
      !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
      PSIGS(:,:,:) = PSIGS(:,:,:) + PAMOIST(:,:,:) **2 * ZFLXZ(:,:,:)
      !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
    END IF
    ! stores <Rnp Rnp>    
    IF ( OTURB_FLX .AND. TPFILE%LOPENED ) THEN
      TZFIELD%CMNHNAME   = 'RTOT_VVAR'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = 'RTOT_VVAR'
      TZFIELD%CUNITS     = 'kg2 kg-2'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CCOMMENT   = 'X_Y_Z_RTOT_VVAR'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
    END IF
    !
    ! and we store in LES configuration
    !
    IF (OLES_CALL) THEN
      CALL SECOND_MNH(ZTIME1)
      CALL LES_MEAN_SUBGRID(ZFLXZ, X_LES_SUBGRID_Rt2 ) 
      CALL LES_MEAN_SUBGRID(MZF(PWM, D%NKA, D%NKU, D%NKL)*ZFLXZ, X_LES_RES_W_SBG_Rt2 )
      CALL LES_MEAN_SUBGRID(PEMOIST*ZFLXZ, X_LES_SUBGRID_RtThv , .TRUE. ) 
      CALL LES_MEAN_SUBGRID(-CSTURB%XA3*PBETA*PEMOIST*ZFLXZ, X_LES_SUBGRID_RtPz, .TRUE. )
      CALL LES_MEAN_SUBGRID(-2.*CSTURB%XCTD*PSQRT_TKE*ZFLXZ/PLEPS, X_LES_SUBGRID_DISS_Rt2 ) 
      CALL SECOND_MNH(ZTIME2)
      XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
    END IF
    !
  END IF  ! end if KRR ne 0
!
!
!        4.5  Vertical part of Sigma_s
!
  IF ( KRRL > 0 ) THEN
    ! Extrapolate PSIGS at the ground and at the top
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
    PSIGS(:,:,D%NKA) = PSIGS(:,:,IKB)
    PSIGS(:,:,D%NKU) = PSIGS(:,:,IKE)
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT)
    !$mnh_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
#ifdef REPRO48
    PSIGS(:,:,:) =  MAX (PSIGS(:,:,:) , 0.)
    PSIGS(:,:,:) =  SQRT(PSIGS(:,:,:))
#else
    PSIGS(:,:,:) =  SQRT( MAX (PSIGS(:,:,:) , 1.E-12) )
#endif
    !$mnh_end_expand_array(JI=1:D%NIT,JJ=1:D%NJT,JK=1:D%NKT)
  END IF

!
!        4.6  Deallocate
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_CORR',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER_THERMO_CORR
END MODULE MODE_TURB_VER_THERMO_CORR
