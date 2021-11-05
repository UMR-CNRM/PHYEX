!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    #################### 
     MODULE MODI_TURB_VER 
!    ####################
!
INTERFACE 
!
      SUBROUTINE TURB_VER(KKA,KKU,KKL,KRR,KRRL,KRRI,                &
                      OTURB_FLX,                                    &
                      HTURBDIM,HTOM,PIMPL,PEXPL,                    &
                      PTSTEP, TPFILE,                               &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PCOSSLOPE,PSINSLOPE,                          &
                      PRHODJ,PTHVREF,                               &
                      PSFTHM,PSFRM,PSFSVM,PSFTHP,PSFRP,PSFSVP,      &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU33M,              &
                      PUM,PVM,PWM,PUSLOPEM,PVSLOPEM,PTHLM,PRM,PSVM, &
                      PTKEM,PLM,PLEPS,                              &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,    &
                      PFWTH,PFWR,PFTH2,PFR2,PFTHR,PBL_DEPTH,        &
                      PSBL_DEPTH,PLMO,                              &
                      PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,             &
                      PDP,PTP,PSIGS,PWTH,PWRC,PWSV                  )

!
USE MODD_IO, ONLY: TFILEDATA
!
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
CHARACTER(len=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
CHARACTER(len=4),       INTENT(IN)   ::  HTOM         ! type of Third Order Moment
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! timestep 
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                      ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PZZ          ! altitudes at flux points
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCOSSLOPE    ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSINSLOPE    ! sinus of the angle 
                                      ! between i and the slope vector
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHVREF      ! ref. state Virtual 
                                                      ! Potential Temperature 
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSFTHM,PSFRM ! surface fluxes at time
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVM       ! t - deltat 
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSFTHP,PSFRP ! surface fluxes at time
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVP       ! t + deltat 
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCDUEFF     ! Cd * || u || at time t
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PUM,PVM,PWM,PTHLM 
  ! Wind and potential temperature at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRM          ! Mixing ratios 
                                                      ! at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLEPS        ! dissipative length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exner at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFRAC_ICE    ! ri fraction of rc+ri
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFWTH        ! d(w'2th' )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFWR         ! d(w'2r'  )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFTH2        ! d(w'th'2 )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFR2         ! d(w'r'2  )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFTHR        ! d(w'th'r')/dz
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PBL_DEPTH    ! BL depth
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PSBL_DEPTH   ! SBL depth
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PLMO         ! Monin-Obukhov length
!
REAL, DIMENSION(:,:,:), INTENT(INOUT)   ::  PRUS, PRVS, PRWS, PRTHLS
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS,PRRS
                            ! cumulated sources for the prognostic variables
!
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PDP,PTP   ! Dynamic and thermal
                                                   ! TKE production terms
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PSIGS     ! Vert. part of Sigma_s at t
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWTH       ! heat flux
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWRC       ! cloud water flux
REAL, DIMENSION(:,:,:,:),INTENT(OUT) :: PWSV       ! scalar flux

!
!
END SUBROUTINE TURB_VER
!
END INTERFACE
!
END MODULE MODI_TURB_VER
!
!
!     ###############################################################
      SUBROUTINE TURB_VER(KKA,KKU,KKL,KRR, KRRL, KRRI,              &
                      OTURB_FLX,                                    &
                      HTURBDIM,HTOM,PIMPL,PEXPL,                    &
                      PTSTEP, TPFILE,                               &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PCOSSLOPE,PSINSLOPE,                          &
                      PRHODJ,PTHVREF,                               &
                      PSFTHM,PSFRM,PSFSVM,PSFTHP,PSFRP,PSFSVP,      &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU33M,              &
                      PUM,PVM,PWM,PUSLOPEM,PVSLOPEM,PTHLM,PRM,PSVM, &
                      PTKEM,PLM,PLEPS,                              &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,    &
                      PFWTH,PFWR,PFTH2,PFR2,PFTHR,PBL_DEPTH,        &
                      PSBL_DEPTH,PLMO,                              &
                      PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,             &
                      PDP,PTP,PSIGS,PWTH,PWRC,PWSV                  )
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
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! JL Redelsperger 03/2021 : add Ocean LES case
!!--------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_CTURB
USE MODD_DYN_n,          ONLY: LOCEAN
use modd_field,          only: tfielddata, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_PARAMETERS
USE MODD_LES
USE MODD_NSV,            ONLY: NSV
!
USE MODI_PRANDTL
USE MODI_EMOIST
USE MODI_ETHETA
USE MODI_GRADIENT_M
USE MODI_GRADIENT_W
USE MODI_TURB
USE MODI_TURB_VER_THERMO_FLUX
USE MODI_TURB_VER_THERMO_CORR
USE MODI_TURB_VER_DYN_FLUX
USE MODI_TURB_VER_SV_FLUX
USE MODI_TURB_VER_SV_CORR
USE MODI_LES_MEAN_SUBGRID
USE MODI_SBL_DEPTH
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
REAL,                   INTENT(IN)   ::  PTSTEP       ! timestep 
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                      ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PZZ          ! altitudes at flux points
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCOSSLOPE    ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSINSLOPE    ! sinus of the angle 
                                      ! between i and the slope vector
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHVREF      ! ref. state Virtual 
                                                      ! Potential Temperature 
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSFTHM,PSFRM ! surface fluxes at time
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVM       ! t - deltat 
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSFTHP,PSFRP ! surface fluxes at time
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVP       ! t + deltat 
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCDUEFF     ! Cd * || u || at time t
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PUM,PVM,PWM,PTHLM 
  ! Wind and potential temperature at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRM          ! Mixing ratios 
                                                      ! at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLEPS        ! dissipative length
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFRAC_ICE    ! ri fraction of rc+ri
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFWTH        ! d(w'2th' )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFWR         ! d(w'2r'  )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFTH2        ! d(w'th'2 )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFR2         ! d(w'r'2  )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFTHR        ! d(w'th'r')/dz
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PBL_DEPTH    ! BL depth
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PSBL_DEPTH   ! SBL depth
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PLMO         ! Monin-Obukhov length
!
REAL, DIMENSION(:,:,:), INTENT(INOUT)   ::  PRUS, PRVS, PRWS, PRTHLS
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS,PRRS
                            ! cumulated sources for the prognostic variables
!
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PDP,PTP   ! Dynamic and thermal
                                                   ! TKE production terms
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PSIGS     ! Vert. part of Sigma_s at t
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWTH      ! heat flux
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWRC      ! cloud water flux
REAL, DIMENSION(:,:,:,:),INTENT(OUT) :: PWSV       ! scalar flux

!
!
!
!
!*       0.2  declaration of local variables
!
!JUAN BUG PGI
!!$REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))  ::  &
REAL, ALLOCATABLE, DIMENSION(:,:,:)  ::  &
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

!!$REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3),NSV)  ::  &
REAL, ALLOCATABLE, DIMENSION(:,:,:,:)  ::  &
       ZPSI_SV,  & ! Prandtl number for scalars
       ZREDS1,   & ! 1D Redelsperger number R_sv
       ZRED2THS, & ! 3D Redelsperger number R*2_thsv
       ZRED2RS     ! 3D Redelsperger number R*2_rsv
!
LOGICAL :: GUSERV    ! flag to use water vapor
INTEGER :: IKB,IKE   ! index value for the Beginning
                     ! and the End of the physical domain for the mass points
INTEGER :: JSV       ! loop counter on scalar variables
REAL    :: ZTIME1
REAL    :: ZTIME2
TYPE(TFIELDDATA) :: TZFIELD
!----------------------------------------------------------------------------
ALLOCATE (      ZBETA(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))    ,&
       ZSQRT_TKE(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3)),& 
       ZDTH_DZ(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))  ,&
       ZDR_DZ(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))   ,&
       ZRED2TH3(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3)) ,& 
       ZRED2R3(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))  ,&
       ZRED2THR3(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3)),&
       ZBLL_O_E(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3)) ,&
       ZETHETA(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))  ,&
       ZEMOIST(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))  ,&
       ZREDTH1(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))  ,&
       ZREDR1(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))   ,&
       ZPHI3(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))    ,&
       ZPSI3(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))    ,&
       ZD(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))       ,&
       ZWTHV(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))    ,&
       ZWU(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))      ,&
       ZWV(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))      ,&
       ZTHLP(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))    ,&
       ZRP(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3))     )   

ALLOCATE ( &
 ZPSI_SV(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3),NSV),  &
 ZREDS1(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3),NSV),   &
 ZRED2THS(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3),NSV), &
 ZRED2RS(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3),NSV)   )

!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
!
IKB=KKA+JPVEXT_TURB*KKL
IKE=KKU-JPVEXT_TURB*KKL
!
!
! 3D Redelsperger numbers
!
!
CALL PRANDTL(KKA,KKU,KKL,KRR,KRRI,OTURB_FLX,       &
             HTURBDIM,                             &
             TPFILE,                               &
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
IF (LOCEAN) THEN
  ZBETA = XG*XALPHAOC
ELSE
  ZBETA = XG/PTHVREF
END IF
!
! Square root of Tke
!
ZSQRT_TKE = SQRT(PTKEM)
!
! gradients of mean quantities at previous time-step
!
ZDTH_DZ = GZ_M_W(KKA,KKU,KKL,PTHLM(:,:,:),PDZZ)
ZDR_DZ  = 0.
IF (KRR>0) THEN
ZDR_DZ  = GZ_M_W(KKA,KKU,KKL,PRM(:,:,:,1),PDZZ)
ENDIF
!
!
! Denominator factor in 3rd order terms
!
ZD(:,:,:) = (1.+ZREDTH1+ZREDR1) * (1.+0.5*(ZREDTH1+ZREDR1))
!
! Phi3 and Psi3 Prandtl numbers
!
GUSERV = KRR/=0
!
ZPHI3 = PHI3(ZREDTH1,ZREDR1,ZRED2TH3,ZRED2R3,ZRED2THR3,HTURBDIM,GUSERV)
IF(KRR/=0) THEN
ZPSI3 = PSI3(ZREDR1,ZREDTH1,ZRED2R3,ZRED2TH3,ZRED2THR3,HTURBDIM,GUSERV)
ENDIF
!
! Prandtl numbers for scalars
!
ZPSI_SV = PSI_SV(ZREDTH1,ZREDR1,ZREDS1,ZRED2THS,ZRED2RS,ZPHI3,ZPSI3)
!
! LES diagnostics
!
IF (LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID(ZPHI3,X_LES_SUBGRID_PHI3)
  IF(KRR/=0) THEN
    CALL LES_MEAN_SUBGRID(ZPSI3,X_LES_SUBGRID_PSI3)
  END IF
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
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
!
  CALL  TURB_VER_THERMO_FLUX(KKA,KKU,KKL,KRR,KRRL,KRRI,               &
                        OTURB_FLX,HTURBDIM,HTOM,                      &
                        PIMPL,PEXPL,PTSTEP,                           &
                        TPFILE,                                       &
                        PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                        PRHODJ,PTHVREF,                               &
                        PSFTHM,PSFRM,PSFTHP,PSFRP,                    &
                        PWM,PTHLM,PRM,PSVM,                           &
                        PTKEM,PLM,PLEPS,                              &
                        PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,    &
                        ZBETA, ZSQRT_TKE, ZDTH_DZ, ZDR_DZ, ZRED2TH3,  &
                        ZRED2R3, ZRED2THR3, ZBLL_O_E, ZETHETA,        &
                        ZEMOIST, ZREDTH1, ZREDR1, ZPHI3, ZPSI3, ZD,   &
                        PFWTH,PFWR,PFTH2,PFR2,PFTHR,PBL_DEPTH,ZWTHV,  &
                        PRTHLS,PRRS,ZTHLP,ZRP,PTP,PWTH,PWRC           )
!
  CALL  TURB_VER_THERMO_CORR(KKA,KKU,KKL,KRR,KRRL,KRRI,               &
                        OTURB_FLX,HTURBDIM,HTOM,                      &
                        PIMPL,PEXPL,                                  &
                        TPFILE,                                       &
                        PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,           &
                        PRHODJ,PTHVREF,                               &
                        PSFTHM,PSFRM,PSFTHP,PSFRP,                    &
                        PWM,PTHLM,PRM,PSVM,                           &
                        PTKEM,PLM,PLEPS,                              &
                        PLOCPEXNM,PATHETA,PAMOIST,PSRCM,              &
                        ZBETA, ZSQRT_TKE, ZDTH_DZ, ZDR_DZ, ZRED2TH3,  &
                        ZRED2R3, ZRED2THR3, ZBLL_O_E, ZETHETA,        &
                        ZEMOIST, ZREDTH1, ZREDR1, ZPHI3, ZPSI3, ZD,   &
                        PFWTH,PFWR,PFTH2,PFR2,PFTHR,                  &
                        ZTHLP,ZRP,PSIGS                          )
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
CALL  TURB_VER_DYN_FLUX(KKA,KKU,KKL,                                &
                      OTURB_FLX,KRR,                                &
                      HTURBDIM,PIMPL,PEXPL,PTSTEP,                  &
                      TPFILE,                                       &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PCOSSLOPE,PSINSLOPE,                          &
                      PRHODJ,                                       &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU33M,              &
                      PTHLM,PRM,PSVM,PUM,PVM,PWM,PUSLOPEM,PVSLOPEM, &
                      PTKEM,PLM,ZWU,ZWV,                            &
                      PRUS,PRVS,PRWS,                               &
                      PDP,PTP                                       )
!
!----------------------------------------------------------------------------
!
!
!*       8.   SOURCES OF PASSIVE SCALAR VARIABLES
!             -----------------------------------
!
IF (SIZE(PSVM,4)>0)                                                 &
CALL  TURB_VER_SV_FLUX(KKA,KKU,KKL,                                 &
                      OTURB_FLX,HTURBDIM,                           &
                      PIMPL,PEXPL,PTSTEP,                           &
                      TPFILE,                                       &
                      PDZZ,PDIRCOSZW,                               &
                      PRHODJ,PWM,                                   &
                      PSFSVM,PSFSVP,                                &
                      PSVM,                                         &
                      PTKEM,PLM,ZPSI_SV,                            &
                      PRSVS,PWSV                                    )
!
!
IF (SIZE(PSVM,4)>0 .AND. LLES_CALL)                                 &
CALL  TURB_VER_SV_CORR(KKA,KKU,KKL,KRR,KRRL,KRRI,                   &
                      PDZZ,                                         &
                      PTHLM,PRM,PTHVREF,                            &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,ZPHI3,ZPSI3,  &
                      PWM,PSVM,                                     &
                      PTKEM,PLM,PLEPS,ZPSI_SV                       )
!
!
!----------------------------------------------------------------------------
!
!*       9.   DIAGNOSTIC OF Surface Boundary Layer Depth
!             ------------------------------------------
!
IF (SIZE(PSBL_DEPTH)>0) CALL SBL_DEPTH(IKB,IKE,PZZ,ZWU,ZWV,ZWTHV,PLMO,PSBL_DEPTH)
!
!----------------------------------------------------------------------------
!
!
!*      10.   PRINTS
!             ------
!
!
IF ( OTURB_FLX .AND. tpfile%lopened ) THEN
!
! stores the Turbulent Prandtl number
! 
  TZFIELD%CMNHNAME   = 'PHI3'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'PHI3'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'Turbulent Prandtl number'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZPHI3)
!
! stores the Turbulent Schmidt number
! 
  TZFIELD%CMNHNAME   = 'PSI3'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'PSI3'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'Turbulent Schmidt number'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZPSI3)
!
!
! stores the Turbulent Schmidt number for the scalar variables
! 
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  DO JSV=1,NSV
    WRITE(TZFIELD%CMNHNAME, '("PSI_SV_",I3.3)') JSV
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
    CALL IO_Field_write(TPFILE,TZFIELD,ZPSI_SV(:,:,:,JSV))
  END DO
!
END IF
!
!
!----------------------------------------------------------------------------
END SUBROUTINE TURB_VER                                                                                               
