!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TURB_VER_THERMO_FLUX
IMPLICIT NONE
CONTAINS
      
SUBROUTINE TURB_VER_THERMO_FLUX(KKA,KKU,KKL,KRR,KRRL,KRRI,    &
                      OTURB_FLX,HTURBDIM,HTOM,                      &
                      PIMPL,PEXPL,                                  &
                      PTSTEP,                                       &
                      TPFILE,                                       &
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
                      PWTHV,PRTHLS,PRRS,PTHLP,PRP,PTP,PWTH,PWRC     )
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
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                     2021 (D. Ricard) last version of HGRAD turbulence scheme
!!                                 Leronard terms instead of Reynolds terms
!!                                 applied to vertical fluxes of r_np and Thl
!!                                 for implicit version of turbulence scheme
!!                                 corrections and cleaning
!!                     June 2020 (B. Vie) Patch preventing negative rc and ri in 2.3 and 3.3
!! JL Redelsperger  : 03/2021: Ocean and Autocoupling O-A LES Cases
!!                             Sfc flux shape for LDEEPOC Case
!!--------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
USE MODD_CST
USE MODD_CONF, ONLY: CPROGRAM
USE MODD_CTURB
USE MODD_FIELD,          ONLY: TFIELDDATA, TYPEREAL
USE MODD_GRID_n,         ONLY: XZS, XXHAT, XYHAT
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_METRICS_n,      ONLY: XDXX, XDYY, XDZX, XDZY, XDZZ
USE MODD_PARAMETERS
USE MODD_TURB_n,         ONLY: LHGRAD, XCOEFHGRADTHL, XCOEFHGRADRM, XALTHGRAD, XCLDTHOLD
USE MODD_CONF
USE MODD_LES
USE MODD_DIM_n
USE MODD_DYN_n,          ONLY: LOCEAN
USE MODD_OCEANH
USE MODD_REF,            ONLY: LCOUPLES
USE MODD_TURB_n
USE MODD_FRC
!
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_GRADIENT_M
USE MODI_SHUMAN , ONLY : DZF, DZM, MZF, MZM, MYF, MXF
USE MODI_TRIDIAG 
USE MODI_LES_MEAN_SUBGRID
USE MODI_TRIDIAG_THERMO
USE MODI_TM06_H
!
USE MODE_IO_FIELD_WRITE, ONLY: IO_FIELD_WRITE
USE MODE_PRANDTL
!
USE MODI_SECOND_MNH
USE MODE_ll
USE MODE_GATHER_ll
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
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
CHARACTER(len=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
CHARACTER(len=4),       INTENT(IN)   ::  HTOM         ! type of Third Order Moment
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
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
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PTP       ! Dynamic and thermal
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
       Z3RDMOMENT,&  ! 3 order term in flux or variance equation
       ZF_LEONARD,&  ! Leonard terms
       ZRWTHL,    &
       ZRWRNP,    &
       ZCLD_THOLD
!
REAL,DIMENSION(SIZE(XZS,1),SIZE(XZS,2),KKU)  :: ZALT
!
INTEGER             :: IKB,IKE      ! I index values for the Beginning and End
                                    ! mass points of the domain in the 3 direct.
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain 
INTEGER             :: JI, JJ ! loop indexes 
!
INTEGER                    :: IIB,IJB       ! Lower bounds of the physical
                                            ! sub-domain in x and y directions
INTEGER                    :: IIE,IJE       ! Upper bounds of the physical
                                            ! sub-domain in x and y directions
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZXHAT_ll    !  Position x in the conformal
                                                 ! plane (array on the complete domain)
REAL, DIMENSION(:), ALLOCATABLE   :: ZYHAT_ll    !   Position y in the conformal
                                                 ! plane (array on the complete domain)
!
!
REAL :: ZTIME1, ZTIME2
REAL :: ZDELTAX
REAL :: ZXBEG,ZXEND,ZYBEG,ZYEND ! Forcing size for ocean deep convection
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT)) :: ZDIST ! distance
                                   ! from the center of the cooling               
REAL :: ZFLPROV
INTEGER           :: JKM          ! vertical index loop
INTEGER           :: JSW
REAL :: ZSWA     ! index for time flux interpolation
!
INTEGER :: IIU, IJU
INTEGER :: IRESP
INTEGER :: JK
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
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_FLUX',0,ZHOOK_HANDLE)
!
! Size for a given proc & a given model      
IIU=SIZE(PTHLM,1) 
IJU=SIZE(PTHLM,2)
!
!! Compute Shape of sfc flux for Oceanic Deep Conv Case
! 
IF (LOCEAN .AND. LDEEPOC) THEN
  !*       COMPUTES THE PHYSICAL SUBDOMAIN BOUNDS
  ALLOCATE(ZXHAT_ll(NIMAX_ll+2*JPHEXT),ZYHAT_ll(NJMAX_ll+2*JPHEXT))
  !compute ZXHAT_ll = position in the (0:Lx) domain 1 (Lx=Size of domain1 )
  !compute XXHAT_ll = position in the (L0_subproc,Lx_subproc) domain for the current subproc
  !                                     L0_subproc as referenced in the full domain 1
  CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,IRESP)
  CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,IRESP)
  CALL GET_DIM_EXT_ll('B',IIU,IJU)
  CALL GET_INDICE_ll(IIB,IJB,IIE,IJE,IIU,IJU)
  DO JJ = IJB,IJE
    DO JI = IIB,IIE
      ZDIST(JI,JJ) = SQRT(                         &
      (( (XXHAT(JI)+XXHAT(JI+1))*0.5 - XCENTX_OC ) / XRADX_OC)**2 + &
      (( (XYHAT(JJ)+XYHAT(JJ+1))*0.5 - XCENTY_OC ) / XRADY_OC)**2   &
                                )
    END DO
  END DO
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IF ( ZDIST(JI,JJ) > 1.) XSSTFL(JI,JJ)=0.
    END DO
  END DO
END IF !END DEEP OCEAN CONV CASE
!
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
! Define a cloud mask with ri and rc (used after with a threshold) for Leonard terms
!
IF(LHGRAD) THEN
  IF ( KRRL >= 1 ) THEN
    IF ( KRRI >= 1 ) THEN
      ZCLD_THOLD(:,:,:) = PRM(:,:,:,2) + PRM(:,:,:,4)
    ELSE
      ZCLD_THOLD(:,:,:) = PRM(:,:,:,2)
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
END IF
!
IF (LHGRAD) THEN
 ! Compute the Leonard terms for thl
 ZDELTAX= XXHAT(3) - XXHAT(2)
 ZF_LEONARD (:,:,:)= XCOEFHGRADTHL*ZDELTAX*ZDELTAX/12.0*(      &
                 MXF(GX_W_UW(PWM(:,:,:), XDXX,XDZZ,XDZX,KKA,KKU,KKL))&
                *MZM(GX_M_M(PTHLM(:,:,:),XDXX,XDZZ,XDZX,KKA, KKU, KKL), KKA, KKU, KKL)  &
              +  MYF(GY_W_VW(PWM(:,:,:), XDYY,XDZZ,XDZY,KKA,KKU,KKL))  &
                *MZM(GY_M_M(PTHLM(:,:,:),XDYY,XDZZ,XDZY,KKA, KKU, KKL), KKA, KKU, KKL) )
END IF
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
  ZF       = ZF       + M3_WTH_W2R(KKA,KKU,KKL,PD,ZKEFF,&
    & PTKEM,PBLL_O_E,PEMOIST,PDTH_DZ) * PFWR
  ZDFDDTDZ = ZDFDDTDZ + D_M3_WTH_W2R_O_DDTDZ(KKA,KKU,KKL,PREDTH1,PREDR1,&
    & PD,ZKEFF,PTKEM,PBLL_O_E,PEMOIST) * PFWR
END IF
!
! d(w'r'2)/dz
IF (GFR2) THEN
  ZF       = ZF       + M3_WTH_WR2(KKA,KKU,KKL,PD,ZKEFF,PTKEM,&
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
! compute interface flux
IF (LCOUPLES) THEN   ! Autocoupling O-A LES
  IF (LOCEAN) THEN    ! ocean model in coupled case 
    ZF(:,:,IKE) =  (XSSTFL_C(:,:,1)+XSSRFL_C(:,:,1)) &
                  *0.5* ( 1. + PRHODJ(:,:,KKU)/PRHODJ(:,:,IKE) )
  ELSE                ! atmosph model in coupled case
    ZF(:,:,IKB) =  XSSTFL_C(:,:,1) &
                  *0.5* ( 1. + PRHODJ(:,:,KKA)/PRHODJ(:,:,IKB) )
  ENDIF 
!
ELSE  ! No coupling O and A cases
  ! atmosp bottom
  !*In 3D, a part of the flux goes vertically,
  ! and another goes horizontally (in presence of slopes)
  !*In 1D, part of energy released in horizontal flux is taken into account in the vertical part
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
  IF (LOCEAN) THEN
    ZF(:,:,IKE) = XSSTFL(:,:) *0.5*(1. + PRHODJ(:,:,KKU) / PRHODJ(:,:,IKE))
  ELSE !end ocean case (in nocoupled case)
    ! atmos top
#ifdef REPRO48
#else
      ZF(:,:,IKE)=0.
#endif
  END IF
END IF !end no coupled cases
!
! Compute the split conservative potential temperature at t+deltat
CALL TRIDIAG_THERMO(KKA,KKU,KKL,PTHLM,ZF,ZDFDDTDZ,PTSTEP,PIMPL,PDZZ,&
                    PRHODJ,PTHLP)
!
! Compute the equivalent tendency for the conservative potential temperature
!
ZRWTHL(:,:,:)= PRHODJ(:,:,:)*(PTHLP(:,:,:)-PTHLM(:,:,:))/PTSTEP
! replace the flux by the Leonard terms above ZALT and ZCLD_THOLD
IF (LHGRAD) THEN
 DO JK=1,KKU
  ZALT(:,:,JK) = PZZ(:,:,JK)-XZS(:,:)
 END DO
 WHERE ( (ZCLD_THOLD(:,:,:) >= XCLDTHOLD) .AND. ( ZALT(:,:,:) >= XALTHGRAD) )
  ZRWTHL(:,:,:) = -GZ_W_M(MZM(PRHODJ(:,:,:), KKA, KKU, KKL)*ZF_LEONARD(:,:,:),XDZZ,&
                   KKA, KKU, KKL)
 END WHERE
END IF
!
PRTHLS(:,:,:)= PRTHLS(:,:,:)  + ZRWTHL(:,:,:)
!
!*       2.2  Partial Thermal Production
!
!  Conservative potential temperature flux : 
!
ZFLXZ(:,:,:)   = ZF                                                &
               + PIMPL * ZDFDDTDZ * DZM(PTHLP - PTHLM, KKA, KKU, KKL) / PDZZ
! replace the flux by the Leonard terms
IF (LHGRAD) THEN
 WHERE ( (ZCLD_THOLD(:,:,:) >= XCLDTHOLD) .AND. ( ZALT(:,:,:) >= XALTHGRAD) )
  ZFLXZ(:,:,:) = ZF_LEONARD(:,:,:)
 END WHERE
END IF
!
ZFLXZ(:,:,KKA) = ZFLXZ(:,:,IKB) 
IF (LOCEAN) THEN
  ZFLXZ(:,:,KKU) = ZFLXZ(:,:,IKE)
END IF
!  
DO JK=IKTB+1,IKTE-1
  PWTH(:,:,JK)=0.5*(ZFLXZ(:,:,JK)+ZFLXZ(:,:,JK+KKL))
END DO
!
PWTH(:,:,IKB)=0.5*(ZFLXZ(:,:,IKB)+ZFLXZ(:,:,IKB+KKL)) 
!
IF (LOCEAN) THEN
  PWTH(:,:,IKE)=0.5*(ZFLXZ(:,:,IKE)+ZFLXZ(:,:,IKE+KKL))
  PWTH(:,:,KKA)=0. 
  PWTH(:,:,KKU)=ZFLXZ(:,:,KKU)
ELSE
  PWTH(:,:,KKA)=0.5*(ZFLXZ(:,:,KKA)+ZFLXZ(:,:,KKA+KKL))
  PWTH(:,:,IKE)=PWTH(:,:,IKE-KKL)
END IF
!
IF ( OTURB_FLX .AND. TPFILE%LOPENED ) THEN
  ! stores the conservative potential temperature vertical flux
  TZFIELD%CMNHNAME   = 'THW_FLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'THW_FLX'
  TZFIELD%CUNITS     = 'K m s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'Conservative potential temperature vertical flux'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
END IF
!
! Contribution of the conservative temperature flux to the buoyancy flux
IF (LOCEAN) THEN
  PTP(:,:,:)= XG*XALPHAOC * MZF(ZFLXZ,KKA, KKU, KKL )
ELSE
  IF (KRR /= 0) THEN
    PTP(:,:,:)  =  PBETA * MZF( MZM(PETHETA,KKA, KKU, KKL) * ZFLXZ,KKA, KKU, KKL )
    PTP(:,:,IKB)=  PBETA(:,:,IKB) * PETHETA(:,:,IKB) *   &
                   0.5 * ( ZFLXZ (:,:,IKB) + ZFLXZ (:,:,IKB+KKL) )
  ELSE
    PTP(:,:,:)=  PBETA * MZF( ZFLXZ,KKA, KKU, KKL )
  END IF
END IF 
!
! Buoyancy flux at flux points
! 
PWTHV = MZM(PETHETA, KKA, KKU, KKL) * ZFLXZ
PWTHV(:,:,IKB) = PETHETA(:,:,IKB) * ZFLXZ(:,:,IKB)
!
IF (LOCEAN) THEN
  ! temperature contribution to Buy flux     
  PWTHV(:,:,IKE) = PETHETA(:,:,IKE) * ZFLXZ(:,:,IKE)
END IF
!*       2.3  Partial vertical divergence of the < Rc w > flux
! Correction for qc and qi negative in AROME 
IF(CPROGRAM/='AROME  ') THEN
 IF ( KRRL >= 1 ) THEN
   IF ( KRRI >= 1 ) THEN
     PRRS(:,:,:,2) = PRRS(:,:,:,2) -                                        &
                     PRHODJ*PATHETA*2.*PSRCM*DZF(ZFLXZ/PDZZ, KKA, KKU, KKL) &
                     *(1.0-PFRAC_ICE(:,:,:))
     PRRS(:,:,:,4) = PRRS(:,:,:,4) -                                        &
                     PRHODJ*PATHETA*2.*PSRCM*DZF(ZFLXZ/PDZZ, KKA, KKU, KKL) &
                     *PFRAC_ICE(:,:,:)
   ELSE
     PRRS(:,:,:,2) = PRRS(:,:,:,2) -                                        &
                     PRHODJ*PATHETA*2.*PSRCM*DZF(ZFLXZ/PDZZ, KKA, KKU, KKL)
   END IF
 END IF
END IF
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
  ! Compute Leonard Terms for Cloud mixing ratio
  IF (LHGRAD) THEN
    ZDELTAX= XXHAT(3) - XXHAT(2)
    ZF_LEONARD (:,:,:)= XCOEFHGRADRM*ZDELTAX*ZDELTAX/12.0*(        &
                MXF(GX_W_UW(PWM(:,:,:),  XDXX,XDZZ,XDZX,KKA,KKU,KKL))       &
                *MZM(GX_M_M(PRM(:,:,:,1),XDXX,XDZZ,XDZX,KKA,KKU,KKL),KKA,KKU,KKL) &
                +MYF(GY_W_VW(PWM(:,:,:), XDYY,XDZZ,XDZY,KKA,KKU,KKL))        &
                *MZM(GY_M_M(PRM(:,:,:,1),XDYY,XDZZ,XDZY,KKA,KKU,KKL),KKA,KKU,KKL) )
   END IF
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
    ZF       = ZF       + M3_WR_W2TH(KKA,KKU,KKL,PD,ZKEFF,&
     & PTKEM,PBLL_O_E,PETHETA,PDR_DZ) * PFWTH
    ZDFDDRDZ = ZDFDDRDZ + D_M3_WR_W2TH_O_DDRDZ(KKA,KKU,KKL,PREDR1,PREDTH1,& 
     & PD,ZKEFF,PTKEM,PBLL_O_E,PETHETA) * PFWTH
  END IF
  !
  ! d(w'th'2)/dz
  IF (GFTH2) THEN
    ZF       = ZF       + M3_WR_WTH2(KKA,KKU,KKL,PD,ZKEFF,PTKEM,&
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
  ! compute interface flux
  IF (LCOUPLES) THEN   ! coupling NH O-A
    IF (LOCEAN) THEN    ! ocean model in coupled case
      ! evap effect on salinity to be added later !!!
      ZF(:,:,IKE) =  0.
    ELSE                ! atmosph model in coupled case
      ZF(:,:,IKB) =  0.
      ! AJOUTER FLUX EVAP SUR MODELE ATMOS
    ENDIF
  !
  ELSE  ! No coupling NH OA case
    ! atmosp bottom
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
    IF (LOCEAN) THEN
      ! General ocean case
      ! salinity/evap effect to be added later !!!!!
      ZF(:,:,IKE) = 0.
    ELSE !end ocean case (in nocoupled case)
      ! atmos top
#ifdef REPRO48
#else
      ZF(:,:,IKE)=0.
#endif
    END IF
  END IF!end no coupled cases
  ! Compute the split conservative potential temperature at t+deltat
  CALL TRIDIAG_THERMO(KKA,KKU,KKL,PRM(:,:,:,1),ZF,ZDFDDRDZ,PTSTEP,PIMPL,&
                      PDZZ,PRHODJ,PRP)
  !
  ! Compute the equivalent tendency for the conservative mixing ratio
  !
  ZRWRNP (:,:,:) = PRHODJ(:,:,:)*(PRP(:,:,:)-PRM(:,:,:,1))/PTSTEP
  !
  ! replace the flux by the Leonard terms above ZALT and ZCLD_THOLD
  IF (LHGRAD) THEN
   DO JK=1,KKU
    ZALT(:,:,JK) = PZZ(:,:,JK)-XZS(:,:)
   END DO
   WHERE ( (ZCLD_THOLD(:,:,:) >= XCLDTHOLD ) .AND. ( ZALT(:,:,:) >= XALTHGRAD ) )
    ZRWRNP (:,:,:) =  -GZ_W_M(MZM(PRHODJ(:,:,:),KKA,KKU,KKL)*ZF_LEONARD(:,:,:),XDZZ,KKA,KKU,KKL)
   END WHERE
  END IF
  !
  PRRS(:,:,:,1) = PRRS(:,:,:,1) + ZRWRNP (:,:,:)
  !
  !*       3.2  Complete thermal production
  !
  ! cons. mixing ratio flux :
  !
  ZFLXZ(:,:,:)   = ZF                                                &
                 + PIMPL * ZDFDDRDZ * DZM(PRP - PRM(:,:,:,1), KKA, KKU, KKL) / PDZZ 
  !
  ! replace the flux by the Leonard terms above ZALT and ZCLD_THOLD
  IF (LHGRAD) THEN
   WHERE ( (ZCLD_THOLD(:,:,:) >= XCLDTHOLD ) .AND. ( ZALT(:,:,:) >= XALTHGRAD ) )
    ZFLXZ(:,:,:) = ZF_LEONARD(:,:,:)
   END WHERE
  END IF
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
  IF ( OTURB_FLX .AND. TPFILE%LOPENED ) THEN
    ! stores the conservative mixing ratio vertical flux
    TZFIELD%CMNHNAME   = 'RCONSW_FLX'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'RCONSW_FLX'
    TZFIELD%CUNITS     = 'kg m s-1 kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'Conservative mixing ratio vertical flux'
    TZFIELD%NGRID      = 4
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
  END IF
  !
  ! Contribution of the conservative water flux to the Buoyancy flux
  IF (LOCEAN) THEN
     ZA(:,:,:)=  -XG*XBETAOC  * MZF(ZFLXZ, KKA, KKU, KKL )
  ELSE
    ZA(:,:,:)   =  PBETA * MZF( MZM(PEMOIST, KKA, KKU, KKL) * ZFLXZ,KKA,KKU,KKL )
    ZA(:,:,IKB) =  PBETA(:,:,IKB) * PEMOIST(:,:,IKB) *   &
                   0.5 * ( ZFLXZ (:,:,IKB) + ZFLXZ (:,:,IKB+KKL) )
    PTP(:,:,:) = PTP(:,:,:) + ZA(:,:,:)
  END IF
  !
  ! Buoyancy flux at flux points
  ! 
  PWTHV          = PWTHV          + MZM(PEMOIST, KKA, KKU, KKL) * ZFLXZ
  PWTHV(:,:,IKB) = PWTHV(:,:,IKB) + PEMOIST(:,:,IKB) * ZFLXZ(:,:,IKB)
  IF (LOCEAN) THEN
    PWTHV(:,:,IKE) = PWTHV(:,:,IKE) + PEMOIST(:,:,IKE)* ZFLXZ(:,:,IKE)
  END IF   
!
!*       3.3  Complete vertical divergence of the < Rc w > flux
! Correction of qc and qi negative for AROME
IF(CPROGRAM/='AROME  ') THEN
   IF ( KRRL >= 1 ) THEN
     IF ( KRRI >= 1 ) THEN
       PRRS(:,:,:,2) = PRRS(:,:,:,2) -                                        &
                       PRHODJ*PAMOIST*2.*PSRCM*DZF(ZFLXZ/PDZZ,KKA,KKU,KKL )       &
                       *(1.0-PFRAC_ICE(:,:,:))
       PRRS(:,:,:,4) = PRRS(:,:,:,4) -                                        &
                       PRHODJ*PAMOIST*2.*PSRCM*DZF(ZFLXZ/PDZZ,KKA,KKU,KKL )       &
                       *PFRAC_ICE(:,:,:)
     ELSE
       PRRS(:,:,:,2) = PRRS(:,:,:,2) -                                        &
                       PRHODJ*PAMOIST*2.*PSRCM*DZF(ZFLXZ/PDZZ,KKA,KKU,KKL )
     END IF
   END IF
END IF
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
IF ( ((OTURB_FLX .AND. TPFILE%LOPENED) .OR. LLES_CALL) .AND. (KRRL > 0) ) THEN
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
  IF ( OTURB_FLX .AND. TPFILE%LOPENED ) THEN
    TZFIELD%CMNHNAME   = 'RCW_FLX'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'RCW_FLX'
    TZFIELD%CUNITS     = 'kg m s-1 kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'Liquid water mixing ratio vertical flux'
    TZFIELD%NGRID      = 4
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
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
IF (LOCEAN .AND. LDEEPOC) THEN
  DEALLOCATE(ZXHAT_ll,ZYHAT_ll)
END IF
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_FLUX',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER_THERMO_FLUX
END MODULE MODE_TURB_VER_THERMO_FLUX
