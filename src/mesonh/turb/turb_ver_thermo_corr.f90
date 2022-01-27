!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    #################### 
     MODULE MODI_TURB_VER_THERMO_CORR
!    ####################
!
INTERFACE 
!
      SUBROUTINE TURB_VER_THERMO_CORR(KKA,KKU,KKL,KRR,KRRL,KRRI,    &
                      OTURB_FLX,HTURBDIM,HTOM,                      &
                      PIMPL,PEXPL,                                  &
                      TPFILE,                                       &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,           &
                      PRHODJ,PTHVREF,                               &
                      PSFTHM,PSFRM,PSFTHP,PSFRP,                    &
                      PWM,PTHLM,PRM,PSVM,                           &
                      PTKEM,PLM,PLEPS,                              &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,              &
                      PBETA, PSQRT_TKE, PDTH_DZ, PDR_DZ, PRED2TH3,  &
                      PRED2R3, PRED2THR3, PBLL_O_E, PETHETA,        &
                      PEMOIST, PREDTH1, PREDR1, PPHI3, PPSI3, PD,   &
                      PFWTH,PFWR,PFTH2,PFR2,PFTHR,                  &
                      PTHLP,PRP,PSIGS                          )
!
USE MODD_IO, ONLY: TFILEDATA
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=AR 
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
CHARACTER(len=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
CHARACTER(len=4),       INTENT(IN)   ::  HTOM         ! type of Third Order Moment
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDZZ, PDXX, PDYY, PDZX, PDZY
                                                      ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
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
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLEPS        ! dissipative length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
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
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHLP      ! guess of thl at t+ deltat
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRP        ! guess of r at t+ deltat
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)  ::  PSIGS     ! Vert. part of Sigma_s at t
!
!
!
END SUBROUTINE TURB_VER_THERMO_CORR
!
END INTERFACE
!
END MODULE MODI_TURB_VER_THERMO_CORR
!
!
!     ###############################################################
      SUBROUTINE TURB_VER_THERMO_CORR(KKA,KKU,KKL,KRR, KRRL, KRRI,  &
                      OTURB_FLX,HTURBDIM,HTOM,                      &
                      PIMPL,PEXPL,                                  &
                      TPFILE,                                       &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,           &
                      PRHODJ,PTHVREF,                               &
                      PSFTHM,PSFRM,PSFTHP,PSFRP,                    &
                      PWM,PTHLM,PRM,PSVM,                           &
                      PTKEM,PLM,PLEPS,                              &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,              &
                      PBETA, PSQRT_TKE, PDTH_DZ, PDR_DZ, PRED2TH3,  &
                      PRED2R3, PRED2THR3, PBLL_O_E, PETHETA,        &
                      PEMOIST, PREDTH1, PREDR1, PPHI3, PPSI3, PD,   &
                      PFWTH,PFWR,PFTH2,PFR2,PFTHR,                  &
                      PTHLP,PRP,PSIGS                          )
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
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!                     2012-02 (Y. Seity) add possibility to run with reversed 
!!                                              vertical levels
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!--------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_CTURB
use modd_field,          only: tfielddata, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_LES
!
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_GRADIENT_M
USE MODI_SHUMAN 
USE MODI_TRIDIAG 
USE MODI_LES_MEAN_SUBGRID
USE MODI_PRANDTL
USE MODI_TRIDIAG_THERMO
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
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDZZ, PDXX, PDYY, PDZX, PDZY
                                                      ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
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
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLEPS        ! dissipative length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM        ! normalized 
! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
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
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHLP      ! guess of thl at t+ deltat
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRP        ! guess of r at t+ deltat
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)  ::  PSIGS     ! Vert. part of Sigma_s at t
!
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))  ::  &
       ZA,       & ! work variable for wrc
       ZFLXZ,    & ! vertical flux of the treated variable
       ZSOURCE,  & ! source of evolution for the treated variable
       ZKEFF,    & ! effectif diffusion coeff = LT * SQRT( TKE )
       ZF,       & ! Flux in dTh/dt =-dF/dz (evaluated at t-1)(or rt instead of Th)
       ZDFDDTDZ, & ! dF/d(dTh/dz)
       ZDFDDRDZ, & ! dF/d(dr/dz)
       Z3RDMOMENT  ! 3 order term in flux or variance equation
INTEGER             :: IKB,IKE      ! I index values for the Beginning and End
                                    ! mass points of the domain in the 3 direct.
INTEGER             :: I1,I2        ! For ZCOEFF allocation
REAL, DIMENSION(:,:,:),ALLOCATABLE  :: ZCOEFF
                                    ! coefficients for the uncentred gradient 
                                    ! computation near the ground
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
IKB=KKA+JPVEXT_TURB*KKL
IKE=KKU-JPVEXT_TURB*KKL
I1=MIN(KKA+JPVEXT_TURB*KKL,KKA+JPVEXT_TURB*KKL+2*KKL)
I2=MAX(KKA+JPVEXT_TURB*KKL,KKA+JPVEXT_TURB*KKL+2*KKL)

ALLOCATE(ZCOEFF(SIZE(PDZZ,1),SIZE(PDZZ,2),I1:I2))
!
GUSERV = (KRR/=0)
!
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
ZCOEFF(:,:,IKB+2*KKL)= - PDZZ(:,:,IKB+KKL) /      &
       ( (PDZZ(:,:,IKB+2*KKL)+PDZZ(:,:,IKB+KKL)) * PDZZ(:,:,IKB+2*KKL) )
ZCOEFF(:,:,IKB+KKL)=   (PDZZ(:,:,IKB+2*KKL)+PDZZ(:,:,IKB+KKL)) /      &
       ( PDZZ(:,:,IKB+KKL) * PDZZ(:,:,IKB+2*KKL) )
ZCOEFF(:,:,IKB)= - (PDZZ(:,:,IKB+2*KKL)+2.*PDZZ(:,:,IKB+KKL)) /      &
       ( (PDZZ(:,:,IKB+2*KKL)+PDZZ(:,:,IKB+KKL)) * PDZZ(:,:,IKB+KKL) )
!
ZKEFF(:,:,:) = MZM( PLM(:,:,:) * SQRT(PTKEM(:,:,:)) )
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
  ZF      (:,:,:) = XCTV*PLM*PLEPS*MZF(PPHI3*PDTH_DZ**2)
  ZDFDDTDZ(:,:,:) = 0.     ! this term, because of discretization, is treated separately
  !
  ! Effect of 3rd order terms in temperature flux (at mass point)
  !
  ! d(w'th'2)/dz
  IF (GFTH2) THEN
    ZF       = ZF       + M3_TH2_WTH2(PREDTH1,PREDR1,PD,PLEPS,&
     & PSQRT_TKE) * PFTH2
    ZDFDDTDZ = ZDFDDTDZ + D_M3_TH2_WTH2_O_DDTDZ(PREDTH1,PREDR1,&
     & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA) * PFTH2
  END IF
  !
  ! d(w'2th')/dz
  IF (GFWTH) THEN
    ZF       = ZF       + M3_TH2_W2TH(PREDTH1,PREDR1,PD,PDTH_DZ,&
     & PLM,PLEPS,PTKEM) * MZF(PFWTH)
    ZDFDDTDZ = ZDFDDTDZ + D_M3_TH2_W2TH_O_DDTDZ(PREDTH1,PREDR1,PD,&
     & PLM,PLEPS,PTKEM,GUSERV) * MZF(PFWTH)
  END IF
  !
  IF (KRR/=0) THEN
    ! d(w'r'2)/dz
    IF (GFR2) THEN
      ZF       = ZF       + M3_TH2_WR2(PD,PLEPS,PSQRT_TKE,PBLL_O_E,&
       & PEMOIST,PDTH_DZ) * PFR2
      ZDFDDTDZ = ZDFDDTDZ + D_M3_TH2_WR2_O_DDTDZ(PREDTH1,PREDR1,PD,&
       & PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTH_DZ) * PFR2
    END IF
    !
    ! d(w'2r')/dz
    IF (GFWR) THEN
      ZF       = ZF       + M3_TH2_W2R(PD,PLM,PLEPS,PTKEM,PBLL_O_E,&
       & PEMOIST,PDTH_DZ) * MZF(PFWR)
      ZDFDDTDZ = ZDFDDTDZ + D_M3_TH2_W2R_O_DDTDZ(PREDTH1,PREDR1,PD,&
       & PLM,PLEPS,PTKEM,PBLL_O_E,PEMOIST,PDTH_DZ) * MZF(PFWR)
    END IF
    !
    ! d(w'th'r')/dz
    IF (GFTHR) THEN
      ZF       = ZF       + M3_TH2_WTHR(PREDR1,PD,PLEPS,PSQRT_TKE,&
       & PBLL_O_E,PEMOIST,PDTH_DZ) * PFTHR
      ZDFDDTDZ = ZDFDDTDZ + D_M3_TH2_WTHR_O_DDTDZ(PREDTH1,PREDR1,&
       & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTH_DZ) * PFTHR
    END IF

  END IF
  !
  ZFLXZ(:,:,:)   = ZF                                                              &
  !     + PIMPL * XCTV*PLM*PLEPS                                                   &
  !        *MZF(D_PHI3DTDZ2_O_DDTDZ(PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,PDTH_DZ,HTURBDIM,GUSERV)   &
  !             *DZM(PTHLP - PTHLM) / PDZZ                                        ) &
        + PIMPL * ZDFDDTDZ * MZF(DZM(PTHLP - PTHLM) / PDZZ )
  !
  ! special case near the ground ( uncentred gradient )
  ZFLXZ(:,:,IKB) = XCTV * PPHI3(:,:,IKB+KKL) * PLM(:,:,IKB)   &
     * PLEPS(:,:,IKB)                                         &
  *( PEXPL *                                                  &
     ( ZCOEFF(:,:,IKB+2*KKL)*PTHLM(:,:,IKB+2*KKL)             &
      +ZCOEFF(:,:,IKB+KKL  )*PTHLM(:,:,IKB+KKL  )             & 
      +ZCOEFF(:,:,IKB      )*PTHLM(:,:,IKB  )   )**2          &
    +PIMPL *                                                  &
     ( ZCOEFF(:,:,IKB+2*KKL)*PTHLP(:,:,IKB+2*KKL)             &
      +ZCOEFF(:,:,IKB+KKL  )*PTHLP(:,:,IKB+KKL  )             &
      +ZCOEFF(:,:,IKB      )*PTHLP(:,:,IKB  )   )**2          &
   ) 
  !
  ZFLXZ(:,:,KKA) = ZFLXZ(:,:,IKB) 
  !
  ZFLXZ = MAX(0., ZFLXZ)
  !
  IF (KRRL > 0)  THEN
    PSIGS(:,:,:) = ZFLXZ(:,:,:) * PATHETA(:,:,:)**2 
  END IF
  !
  !
  ! stores <THl THl>  
  IF ( OTURB_FLX .AND. tpfile%lopened ) THEN
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
  IF (LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( ZFLXZ, X_LES_SUBGRID_Thl2 ) 
    CALL LES_MEAN_SUBGRID( MZF(PWM)*ZFLXZ, X_LES_RES_W_SBG_Thl2 )
    CALL LES_MEAN_SUBGRID( -2.*XCTD*PSQRT_TKE*ZFLXZ/PLEPS, X_LES_SUBGRID_DISS_Thl2 ) 
    CALL LES_MEAN_SUBGRID( PETHETA*ZFLXZ, X_LES_SUBGRID_ThlThv ) 
    CALL LES_MEAN_SUBGRID( -XA3*PBETA*PETHETA*ZFLXZ, X_LES_SUBGRID_ThlPz, .TRUE. ) 
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
    ZF      (:,:,:) = XCTV*PLM*PLEPS*MZF(0.5*(PPHI3+PPSI3)*PDTH_DZ*PDR_DZ)
    ZDFDDTDZ(:,:,:) = 0.     ! this term, because of discretization, is treated separately
    ZDFDDRDZ(:,:,:) = 0.     ! this term, because of discretization, is treated separately
    !
    ! Effect of 3rd order terms in temperature flux (at mass point)
    !
    ! d(w'th'2)/dz
    IF (GFTH2) THEN
      ZF       = ZF       + M3_THR_WTH2(PREDR1,PD,PLEPS,PSQRT_TKE,&
       & PBLL_O_E,PETHETA,PDR_DZ) * PFTH2
      ZDFDDTDZ = ZDFDDTDZ + D_M3_THR_WTH2_O_DDTDZ(PREDTH1,PREDR1,&
       & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ) * PFTH2
      ZDFDDRDZ = ZDFDDRDZ + D_M3_THR_WTH2_O_DDRDZ(PREDTH1,PREDR1,&
       & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA) * PFTH2
    END IF
    !
    ! d(w'2th')/dz
    IF (GFWTH) THEN
      ZF       = ZF       + M3_THR_W2TH(PREDR1,PD,PLM,PLEPS,PTKEM,&
       & PDR_DZ) * MZF(PFWTH)
      ZDFDDTDZ = ZDFDDTDZ + D_M3_THR_W2TH_O_DDTDZ(PREDTH1,PREDR1,&
       & PD,PLM,PLEPS,PTKEM,PBLL_O_E,PDR_DZ,PETHETA) * MZF(PFWTH)
      ZDFDDRDZ = ZDFDDRDZ + D_M3_THR_W2TH_O_DDRDZ(PREDTH1,PREDR1,&
       & PD,PLM,PLEPS,PTKEM) * MZF(PFWTH)
    END IF
    !
    ! d(w'r'2)/dz
    IF (GFR2) THEN
      ZF       = ZF       + M3_THR_WR2(PREDTH1,PD,PLEPS,PSQRT_TKE,&
       & PBLL_O_E,PEMOIST,PDTH_DZ) * PFR2
      ZDFDDTDZ = ZDFDDTDZ + D_M3_THR_WR2_O_DDTDZ(PREDR1,PREDTH1,PD,&
       & PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST) * PFR2
      ZDFDDRDZ = ZDFDDRDZ + D_M3_THR_WR2_O_DDRDZ(PREDR1,PREDTH1,PD,&
       & PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTH_DZ) * PFR2
    END IF
    !
      ! d(w'2r')/dz
    IF (GFWR) THEN
      ZF       = ZF       + M3_THR_W2R(PREDTH1,PD,PLM,PLEPS,PTKEM,&
      & PDTH_DZ) * MZF(PFWR)
      ZDFDDTDZ = ZDFDDTDZ + D_M3_THR_W2R_O_DDTDZ(PREDR1,PREDTH1,PD,&
      & PLM,PLEPS,PTKEM) * MZF(PFWR)
      ZDFDDRDZ = ZDFDDRDZ + D_M3_THR_W2R_O_DDRDZ(PREDR1,PREDTH1,PD,&
      & PLM,PLEPS,PTKEM,PBLL_O_E,PDTH_DZ,PEMOIST) * MZF(PFWR)
    END IF
    !
    ! d(w'th'r')/dz
    IF (GFTHR) THEN
      ZF       = ZF       + M3_THR_WTHR(PREDTH1,PREDR1,PD,PLEPS,&
      & PSQRT_TKE) * PFTHR
      ZDFDDTDZ = ZDFDDTDZ + D_M3_THR_WTHR_O_DDTDZ(PREDTH1,PREDR1,&
      & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA) * PFTHR
      ZDFDDRDZ = ZDFDDRDZ + D_M3_THR_WTHR_O_DDRDZ(PREDR1,PREDTH1,&
      & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST) * PFTHR
    END IF
    !
    ZFLXZ(:,:,:)   = ZF                                                     &
        + PIMPL * XCTV*PLM*PLEPS*0.5                                        &
          * MZF( ( D_PHI3DTDZ_O_DDTDZ(PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,GUSERV) & ! d(phi3*dthdz)/ddthdz term
                  +D_PSI3DTDZ_O_DDTDZ(PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,GUSERV) & ! d(psi3*dthdz)/ddthdz term
                 ) *PDR_DZ  *DZM(PTHLP - PTHLM       ) / PDZZ               &
                +( D_PHI3DRDZ_O_DDRDZ(PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,HTURBDIM,GUSERV) & ! d(phi3*drdz )/ddrdz term
                  +D_PSI3DRDZ_O_DDRDZ(PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,HTURBDIM,GUSERV) & ! d(psi3*drdz )/ddrdz term
                 ) *PDTH_DZ *DZM(PRP   - PRM(:,:,:,1)) / PDZZ               &
               )                                                            &
        + PIMPL * ZDFDDTDZ * MZF(DZM(PTHLP - PTHLM(:,:,:)) / PDZZ )         &
        + PIMPL * ZDFDDRDZ * MZF(DZM(PRP   - PRM(:,:,:,1)) / PDZZ )
    !
    ! special case near the ground ( uncentred gradient )
    ZFLXZ(:,:,IKB) =                                            & 
    (XCHT1 * PPHI3(:,:,IKB+KKL) + XCHT2 * PPSI3(:,:,IKB+KKL))   &
    *( PEXPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*KKL)*PTHLM(:,:,IKB+2*KKL)             &
        +ZCOEFF(:,:,IKB+KKL  )*PTHLM(:,:,IKB+KKL  )             & 
        +ZCOEFF(:,:,IKB      )*PTHLM(:,:,IKB      ))            &
      *( ZCOEFF(:,:,IKB+2*KKL)*PRM(:,:,IKB+2*KKL,1)             &
        +ZCOEFF(:,:,IKB+KKL  )*PRM(:,:,IKB+KKL,1  )             & 
        +ZCOEFF(:,:,IKB      )*PRM(:,:,IKB  ,1    ))            &
      +PIMPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*KKL)*PTHLP(:,:,IKB+2*KKL)             &
        +ZCOEFF(:,:,IKB+KKL  )*PTHLP(:,:,IKB+KKL  )             &
        +ZCOEFF(:,:,IKB      )*PTHLP(:,:,IKB      ))            &
      *( ZCOEFF(:,:,IKB+2*KKL)*PRP(:,:,IKB+2*KKL  )             &
        +ZCOEFF(:,:,IKB+KKL  )*PRP(:,:,IKB+KKL    )             & 
        +ZCOEFF(:,:,IKB      )*PRP(:,:,IKB        ))            &
     ) 
    !    
    ZFLXZ(:,:,KKA) = ZFLXZ(:,:,IKB) 
    !
      IF ( KRRL > 0 ) THEN
      PSIGS(:,:,:) = PSIGS(:,:,:) +     &
                     2. * PATHETA(:,:,:) * PAMOIST(:,:,:) * ZFLXZ(:,:,:)
    END IF
    ! stores <THl Rnp>   
    IF ( OTURB_FLX .AND. tpfile%lopened ) THEN
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
    IF (LLES_CALL) THEN
      CALL SECOND_MNH(ZTIME1)
      CALL LES_MEAN_SUBGRID( ZFLXZ, X_LES_SUBGRID_THlRt ) 
      CALL LES_MEAN_SUBGRID( MZF(PWM)*ZFLXZ, X_LES_RES_W_SBG_ThlRt )
      CALL LES_MEAN_SUBGRID( -2.*XCTD*PSQRT_TKE*ZFLXZ/PLEPS, X_LES_SUBGRID_DISS_ThlRt ) 
      CALL LES_MEAN_SUBGRID( PETHETA*ZFLXZ, X_LES_SUBGRID_RtThv ) 
      CALL LES_MEAN_SUBGRID( -XA3*PBETA*PETHETA*ZFLXZ, X_LES_SUBGRID_RtPz, .TRUE. ) 
      CALL LES_MEAN_SUBGRID( PEMOIST*ZFLXZ, X_LES_SUBGRID_ThlThv , .TRUE. ) 
      CALL LES_MEAN_SUBGRID( -XA3*PBETA*PEMOIST*ZFLXZ, X_LES_SUBGRID_ThlPz, .TRUE. ) 
      CALL SECOND_MNH(ZTIME2)
      XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
    END IF
! 
!
!*       4.4  <Rnp Rnp>
!
!
    ! Compute the turbulent variance F and F' at time t-dt.
    ZF      (:,:,:) = XCTV*PLM*PLEPS*MZF(PPSI3*PDR_DZ**2)
    ZDFDDRDZ(:,:,:) = 0.     ! this term, because of discretization, is treated separately
    !
    ! Effect of 3rd order terms in temperature flux (at mass point)
    !
    ! d(w'r'2)/dz
    IF (GFR2) THEN
      ZF       = ZF       + M3_R2_WR2(PREDR1,PREDTH1,PD,PLEPS,&
      & PSQRT_TKE) * PFR2
      ZDFDDRDZ = ZDFDDRDZ + D_M3_R2_WR2_O_DDRDZ(PREDR1,PREDTH1,&
      & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST) * PFR2
    END IF
    !
    ! d(w'2r')/dz
    IF (GFWR) THEN
      ZF       = ZF       + M3_R2_W2R(PREDR1,PREDTH1,PD,PDR_DZ,&
      & PLM,PLEPS,PTKEM) * MZF(PFWR)
      ZDFDDRDZ = ZDFDDRDZ + D_M3_R2_W2R_O_DDRDZ(PREDR1,PREDTH1,&
      & PD,PLM,PLEPS,PTKEM,GUSERV) * MZF(PFWR)
    END IF
    !
    IF (KRR/=0) THEN
      ! d(w'r'2)/dz
      IF (GFTH2) THEN
        ZF       = ZF       + M3_R2_WTH2(PD,PLEPS,PSQRT_TKE,&
        & PBLL_O_E,PETHETA,PDR_DZ) * PFTH2
        ZDFDDRDZ = ZDFDDRDZ + D_M3_R2_WTH2_O_DDRDZ(PREDR1,&
        & PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ) * PFTH2
      END IF
      !
      ! d(w'2r')/dz
      IF (GFWTH) THEN
        ZF       = ZF       + M3_R2_W2TH(PD,PLM,PLEPS,PTKEM,&
        & PBLL_O_E,PETHETA,PDR_DZ) * MZF(PFWTH)
        ZDFDDRDZ = ZDFDDRDZ + D_M3_R2_W2TH_O_DDRDZ(PREDR1,PREDTH1,&
        & PD,PLM,PLEPS,PTKEM,PBLL_O_E,PETHETA,PDR_DZ) * MZF(PFWTH)
      END IF
      !
      ! d(w'th'r')/dz
      IF (GFTHR) THEN
        ZF       = ZF       + M3_R2_WTHR(PREDTH1,PD,PLEPS,&
        & PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ) * PFTHR
        ZDFDDRDZ = ZDFDDRDZ + D_M3_R2_WTHR_O_DDRDZ(PREDR1,PREDTH1,&
        & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ) * PFTHR
      END IF
  
    END IF
    !
    ZFLXZ(:,:,:)   = ZF                                                              &
          + PIMPL * XCTV*PLM*PLEPS                                                   &
            *MZF(D_PSI3DRDZ2_O_DDRDZ(PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,PDR_DZ,HTURBDIM,GUSERV)    &
                 *DZM(PRP - PRM(:,:,:,1)) / PDZZ                                   ) &
          + PIMPL * ZDFDDRDZ * MZF(DZM(PRP - PRM(:,:,:,1)) / PDZZ )
    !
    ! special case near the ground ( uncentred gradient )
    ZFLXZ(:,:,IKB) = XCHV * PPSI3(:,:,IKB+KKL) * PLM(:,:,IKB)   &
        * PLEPS(:,:,IKB)                                        &
    *( PEXPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*KKL)*PRM(:,:,IKB+2*KKL,1)             &
        +ZCOEFF(:,:,IKB+KKL  )*PRM(:,:,IKB+KKL,1  )             & 
        +ZCOEFF(:,:,IKB      )*PRM(:,:,IKB  ,1    ))**2         &
      +PIMPL *                                                  &
       ( ZCOEFF(:,:,IKB+2*KKL)*PRP(:,:,IKB+2*KKL)               &
        +ZCOEFF(:,:,IKB+KKL  )*PRP(:,:,IKB+KKL  )               &
        +ZCOEFF(:,:,IKB      )*PRP(:,:,IKB      ))**2           &
     ) 
    !
    ZFLXZ(:,:,KKA) = ZFLXZ(:,:,IKB) 
    !
    IF ( KRRL > 0 ) THEN
      PSIGS(:,:,:) = PSIGS(:,:,:) + PAMOIST(:,:,:) **2 * ZFLXZ(:,:,:)
    END IF
    ! stores <Rnp Rnp>    
    IF ( OTURB_FLX .AND. tpfile%lopened ) THEN
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
    IF (LLES_CALL) THEN
      CALL SECOND_MNH(ZTIME1)
      CALL LES_MEAN_SUBGRID( ZFLXZ, X_LES_SUBGRID_Rt2 ) 
      CALL LES_MEAN_SUBGRID( MZF(PWM)*ZFLXZ, X_LES_RES_W_SBG_Rt2 )
      CALL LES_MEAN_SUBGRID( PEMOIST*ZFLXZ, X_LES_SUBGRID_RtThv , .TRUE. ) 
      CALL LES_MEAN_SUBGRID( -XA3*PBETA*PEMOIST*ZFLXZ, X_LES_SUBGRID_RtPz, .TRUE. )
      CALL LES_MEAN_SUBGRID( -2.*XCTD*PSQRT_TKE*ZFLXZ/PLEPS, X_LES_SUBGRID_DISS_Rt2 ) 
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
    PSIGS(:,:,KKA) = PSIGS(:,:,IKB)
    PSIGS(:,:,KKU) = PSIGS(:,:,IKE)
    PSIGS(:,:,:) =  SQRT( MAX (PSIGS(:,:,:) , 1.E-12) )
  END IF

!
!        4.6  Deallocate
!
  DEALLOCATE(ZCOEFF)
!----------------------------------------------------------------------------
END SUBROUTINE TURB_VER_THERMO_CORR
