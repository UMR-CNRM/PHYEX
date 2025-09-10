!MNH_LIC Copyright 1994-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_VER_THERMO_CORR
IMPLICIT NONE
CONTAINS      
SUBROUTINE TURB_VER_THERMO_CORR(D,CST,CSTURB,TURBN,NEBN,TLES,       &
                      KRR,KRRL,KRRI,KSV,                            &
                      OCOMPUTE_SRC,OCOUPLES,                        &
                      PEXPL,TPFILE,                                 &
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
!     terms to the main program.  In the case of large horizontal meshes,
!     the divergence of these vertical turbulent fluxes represent the whole
!     effect of the turbulence but when the three-dimensionnal version of
!     the turbulence scheme is activated (CTURBDIM="3DIM"), these divergences
!     are completed in the next routine TURB_HOR.
!      An arbitrary degree of implicitness has been implemented for the
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
!!         The different prognostic variables are treated one by one.
!!      The contributions of each turbulent fluxes are cumulated into the
!!      tendency  PRvarS, and into the dynamic and thermal production of
!!      TKE if necessary.
!!
!!       In section 2 and 3, the thermodynamical fields are considered.
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
!!       In section 4, the variance of the statistical variable
!!      s indicating presence or not of condensation, is determined in function
!!      of the turbulent moments of the conservative variables and its
!!      squarred root is stored in PSIGS. This information will be completed in
!!      the horizontal turbulence if the turbulence dimensionality is not
!!      equal to "1DIM".
!!
!!       In section 5, the x component of the stress tensor is computed.
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
!!      and in section 7, a diagnostic computation of the W variance is
!!      performed.
!!
!!         In section 8, the turbulent fluxes for the scalar variables are
!!      computed by the same way as the conservative thermodynamical variables
!!
!!
!!    EXTERNAL
!!    --------
!!      GX_U_M, GY_V_M, GZ_W_M :  cartesian gradient operators
!!      GX_U_UW,GY_V_VW           (X,Y,Z) represent the direction of the gradient
!!                               _(M,U,...)_ represent the localization of the
!!                               field to be derivated
!!                               _(M,UW,...) represent the localization of the
!!                               field  derivated
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
!!           TURBN%XCTV,TURBN%XCHV   : cts for the T and moisture variances
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
!!                                 change the surface  relations
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
!!      Modifications  July 2015 (Wim de Rooy) TURBN%LHARAT switch
!!                     04/2016 (M.Moge) Use openACC directives to port the TURB part of Meso-NH on GPU
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Modifications  June 2019 (Wim de Rooy) New set up cloud scheme
!!      Modifications: June 2023 (S. Riette) tunable value for SIGS minimum value
!!--------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODE_SHUMAN_PHY, ONLY: MZM_PHY, MZF_PHY, DZM_PHY
USE YOMHOOK,    ONLY: LHOOK, DR_HOOK, JPHOOK
!
USE MODD_CST,            ONLY: CST_t
USE MODD_CTURB,          ONLY: CSTURB_t
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_FIELD,          ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_LES,            ONLY: TLES_t
USE MODD_PARAMETERS,     ONLY: JPVEXT_TURB
USE MODD_TURB_n,         ONLY: TURB_t
USE MODD_NEB_n,          ONLY: NEB_t
!
USE MODE_IO_FIELD_WRITE_PHY, ONLY: IO_FIELD_WRITE_PHY
USE MODE_PRANDTL
!
USE MODI_LES_MEAN_SUBGRID_PHY
USE MODI_SECOND_MNH
!
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif
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
TYPE(NEB_t),            INTENT(IN)   :: NEBN
TYPE(TLES_t),           INTENT(INOUT):: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KSV           ! number of scalar var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and version 
REAL,                   INTENT(IN)   ::  PEXPL        ! Coef. for temporal disc.
TYPE(TFILEDATA),        INTENT(INOUT)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PDZZ, PDXX, PDYY, PDZX, PDZY
                                                      ! Metric coefficients
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTHVREF      ! ref. state Virtual
                                                      ! Potential Temperature
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSFTHM,PSFRM ! surface fluxes at time
!                                                     ! t - deltat
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSFTHP,PSFRP ! surface fluxes at time
!                                                     ! t + deltat
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PWM
! Vertical wind
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTHLM
! potential temperature at t-Delta t
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN) ::  PRM          ! Mixing ratios
                                                      ! at t-Delta t
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVM         ! Mixing ratios
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTKEM        ! TKE at time t
! In case TURBN%LHARATU=TRUE, PLM already includes all stability corrections
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLM          ! Turb. mixing length
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLEPS        ! dissipative length
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PATHETA      ! coefficients between
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
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
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PTHLP      ! guess of thl at t+ deltat
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRP        ! guess of r at t+ deltat
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  MFMOIST      ! moist mass flux dual scheme
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)),   INTENT(OUT)  ::  PSIGS     ! Vert. part of Sigma_s at t
!
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(D%NIJT,D%NKT)  ::  &
       ZFLXZ,    & ! vertical flux of the treated variable
       ZKEFF,    & ! effectif diffusion coeff = LT * SQRT( TKE )
       ZF,       & ! Flux in dTh/dt =-dF/dz (evaluated at t-1)(or rt instead of Th)
       ZDFDDTDZ, & ! dF/d(dTh/dz)
       ZDFDDRDZ, & ! dF/d(dr/dz)
! Estimate of full level length and dissipation length scale in case TURBN%LHARATU
       PLMF,     & ! estimate full level length scale from half levels (sub optimal)
       PLEPSF,   & ! estimate full level diss length scale from half levels (sub optimal)
       ZWORK1,ZWORK2,&
       ZWORK3,ZWORK4,&
       ZWORK5,ZWORK6,&
       ZWORK7,ZWORK8,&
       ZWKPHIPSI1,ZWKPHIPSI2,&
       ZWKPHIPSI3,ZWKPHIPSI4       ! working var. for shuman operators (array syntax)

INTEGER             :: IIJB, IIJE, IKB,IKE,IKT,IKA ! index value for the mass points of the domain 
INTEGER             :: IKU  ! array sizes
INTEGER             :: IKL
INTEGER             :: JIJ, JK ! loop indexes 

REAL, DIMENSION(D%NIJT,MIN(D%NKA+JPVEXT_TURB*D%NKL,D%NKA+JPVEXT_TURB*D%NKL+2*D%NKL):&
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
TYPE(TFIELDMETADATA) :: TZFIELD
!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
!
REAL, DIMENSION(D%NIJT,D%NKT) ::ZSHUGRADWK1_2D
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMZF2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZSHUGRADWK2_2D
REAL, DIMENSION(D%NIJT,D%NKT) ::ZDZM2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZDZM2D_WORK2
REAL, DIMENSION(D%NIJT,D%NKT) ::ZDZM2D_WORK3
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMZF2D_WORK2
REAL, DIMENSION(D%NIJT,D%NKT) ::ZDZM2D_WORK4
REAL, DIMENSION(D%NIJT,D%NKT) ::ZSHUGRADWK3_2D
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMZF2D_WORK3
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_CORR',0,ZHOOK_HANDLE)
!
IKB=D%NKB
IKE=D%NKE
IKT=D%NKT
IKA=D%NKA
IKU=D%NKU
IKL=D%NKL
IIJE=D%NIJE
IIJB=D%NIJB
!
GUSERV = (KRR/=0)
!
!  compute the coefficients for the uncentred gradient computation near the
!  ground

DO JIJ=IIJB, IIJE
  ZCOEFF(JIJ, IKB+2*IKL)= - PDZZ(JIJ, IKB+IKL) /      &
         ( (PDZZ(JIJ, IKB+2*IKL)+PDZZ(JIJ, IKB+IKL)) * PDZZ(JIJ, IKB+2*IKL) )
  ZCOEFF(JIJ, IKB+IKL)=   (PDZZ(JIJ, IKB+2*IKL)+PDZZ(JIJ, IKB+IKL)) /      &
         ( PDZZ(JIJ, IKB+IKL) * PDZZ(JIJ, IKB+2*IKL) )
  ZCOEFF(JIJ, IKB)= - (PDZZ(JIJ, IKB+2*IKL)+2.*PDZZ(JIJ, IKB+IKL)) /      &
         ( (PDZZ(JIJ, IKB+2*IKL)+PDZZ(JIJ, IKB+IKL)) * PDZZ(JIJ, IKB+IKL) )
END DO

!
!
IF (TURBN%LHARAT) THEN
  CALL MZF_PHY(D,PLM,PLMF)
  !wc Part of the new statistical cloud scheme set up
  IF (NEBN%LSTATNW) THEN
    CALL MZF_PHY(D,PLEPS,PLEPSF)
  ELSE

  DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PLEPSF(JIJ, JK)=PLMF(JIJ, JK)
    END DO
    END DO

  END IF
  !  function MZF produces -999 for level IKU (82 for 80 levels)
  !  so put these to normal value as this level (82) is indeed calculated

  DO JIJ=IIJB, IIJE
    PLMF(JIJ, IKT)=0.001
    PLEPSF(JIJ, IKT)=0.001
  END DO
  ! with energy cascade contribution 50MF term can be omitted
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZKEFF(JIJ, JK) = PLM(JIJ, JK) * SQRT(PTKEM(JIJ, JK))
    END DO
  END DO

ELSE

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK1(JIJ, JK) = PLM(JIJ, JK) * SQRT(PTKEM(JIJ, JK))
    END DO
  END DO

  CALL MZM_PHY(D,ZWORK1,ZKEFF)
ENDIF
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
  GFTH2 = ANY(PFTH2(:,:)/=0.)
  GFR2  = ANY(PFR2(:,:) /=0.) .AND. GUSERV
  GFTHR = ANY(PFTHR(:,:)/=0.) .AND. GUSERV
  GFWTH = ANY(PFWTH(:,:)/=0.)
  GFWR  = ANY(PFWR(:,:) /=0.) .AND. GUSERV
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
  IF (TURBN%LHARAT) THEN

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZWORK1(JIJ, JK)=PDTH_DZ(JIJ, JK)**2
      END DO
    END DO

    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    IF (NEBN%LSTATNW) THEN

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = TURBN%XCTV * & 
                                  PLMF(JIJ, JK)*PLEPSF(JIJ, JK)*ZWORK2(JIJ, JK)
        END DO
      END DO

    ELSE

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = PLMF(JIJ, JK)*PLEPSF(JIJ, JK)*ZWORK2(JIJ, JK)
        END DO
      END DO

    END IF
  ELSE

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZWORK1(JIJ, JK)=PPHI3(JIJ, JK)*PDTH_DZ(JIJ, JK)**2
      END DO
    END DO

    CALL MZF_PHY(D,ZWORK1,ZWORK2)

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZF(JIJ, JK) = TURBN%XCTV*PLM(JIJ, JK)*PLEPS(JIJ, JK)& 
                                     * ZWORK2(JIJ, JK)
      END DO
    END DO

  ENDIF

  ZDFDDTDZ(:,:) = 0.     ! this term, because of discretization, is treated separately

  !
  ! Effect of 3rd order terms in temperature flux (at mass point)
  !
  ! d(w'th'2)/dz
  IF (GFTH2) THEN
    CALL M3_TH2_WTH2(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLEPS,PSQRT_TKE,ZWORK1)
    CALL D_M3_TH2_WTH2_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,&
     & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,ZWORK2)
    !

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZF(JIJ, JK)       = ZF(JIJ, JK) + ZWORK1(JIJ, JK) &
                                          * PFTH2(JIJ, JK)
        ZDFDDTDZ(JIJ, JK) = ZDFDDTDZ(JIJ, JK) + ZWORK2(JIJ, JK) &
                                          * PFTH2(JIJ, JK)
      END DO
    END DO

  END IF
  !
  ! d(w'2th')/dz
  IF (GFWTH) THEN
    CALL M3_TH2_W2TH(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PDTH_DZ,&
     & PLM,PLEPS,PTKEM,ZWORK1)
    CALL MZF_PHY(D,PFWTH,ZWORK2)
    CALL D_M3_TH2_W2TH_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,&
     & PLM,PLEPS,PTKEM,GUSERV,ZWORK3)
    !

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZF(JIJ, JK)  = ZF(JIJ, JK) + ZWORK1(JIJ, JK) &
                                     * ZWORK2(JIJ, JK)
        ZDFDDTDZ(JIJ, JK) = ZDFDDTDZ(JIJ, JK) + ZWORK3(JIJ, JK) &
                                          * ZWORK2(JIJ, JK)
      END DO
    END DO

  END IF
  !
  
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = PTHLP(JIJ, JK) - PTHLM(JIJ, JK)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK2_2D, ZDZM2D_WORK1)
IF (KRR/=0) THEN
    ! d(w'r'2)/dz
    IF (GFR2) THEN
      CALL M3_TH2_WR2(D,CSTURB,TURBN,PD,PLEPS,PSQRT_TKE,PBLL_O_E,&
       & PEMOIST,PDTH_DZ,ZWORK1)
      CALL D_M3_TH2_WR2_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,&
       & PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTH_DZ,ZWORK2)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK1(JIJ, JK) &
                                      * PFR2(JIJ, JK)
          ZDFDDTDZ(JIJ, JK) = ZDFDDTDZ(JIJ, JK) + ZWORK2(JIJ, JK) &
                                            * PFR2(JIJ, JK)
        END DO
      END DO

    END IF
    !
    ! d(w'2r')/dz
    IF (GFWR) THEN
      CALL M3_TH2_W2R(D,CSTURB,TURBN,PD,PLM,PLEPS,PTKEM,PBLL_O_E,&
       & PEMOIST,PDTH_DZ,ZWORK1)
    CALL MZF_PHY(D,PFWR,ZWORK2)
      CALL D_M3_TH2_W2R_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,&
       & PLM,PLEPS,PTKEM,PBLL_O_E,PEMOIST,PDTH_DZ,ZWORK3)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK1(JIJ, JK) &
                                      * ZWORK2(JIJ, JK)
          ZDFDDTDZ(JIJ, JK) = ZDFDDTDZ(JIJ, JK) + ZWORK3(JIJ, JK) &
                                            * ZWORK1(JIJ, JK)
        END DO
      END DO

    END IF
    !
    ! d(w'th'r')/dz
    IF (GFTHR) THEN
      CALL M3_TH2_WTHR(D,CSTURB,TURBN,PREDR1,PD,PLEPS,PSQRT_TKE,&
       & PBLL_O_E,PEMOIST,PDTH_DZ,ZWORK1)
       CALL D_M3_TH2_WTHR_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,&
       & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTH_DZ,ZWORK2)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK1(JIJ, JK) &
                                      * PFTHR(JIJ, JK)
          ZDFDDTDZ(JIJ, JK) = ZDFDDTDZ(JIJ, JK) + ZWORK2(JIJ, JK) &
                                            * PFTHR(JIJ, JK)
        END DO
      END DO

    END IF

  END IF
  
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZDZM2D_WORK1(JIJ, JK) / PDZZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)
!
    
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, JK)   = ZF(JIJ, JK) + TURBN%XIMPL * ZDFDDTDZ(JIJ, JK) &
        * ZMZF2D_WORK1(JIJ, JK)  
  END DO
END DO

!
!
  ! special case near the ground ( uncentred gradient )
  IF (TURBN%LHARAT) THEN

    DO JIJ=IIJB, IIJE
      ZFLXZ(JIJ, IKB) =  PLMF(JIJ, IKB)   &
       * PLEPSF(JIJ, IKB)                                         &
       *( PEXPL *                                                  &
       ( ZCOEFF(JIJ, IKB+2*IKL)*PTHLM(JIJ, IKB+2*IKL)             &
        +ZCOEFF(JIJ, IKB+IKL)*PTHLM(JIJ, IKB+IKL)             & 
        +ZCOEFF(JIJ, IKB)*PTHLM(JIJ, IKB)   )**2          &
       +TURBN%XIMPL *                                                  &
       ( ZCOEFF(JIJ, IKB+2*IKL)*PTHLP(JIJ, IKB+2*IKL)             &
        +ZCOEFF(JIJ, IKB+IKL)*PTHLP(JIJ, IKB+IKL)             &
        +ZCOEFF(JIJ, IKB)*PTHLP(JIJ, IKB)   )**2          &
      ) 
    END DO

    IF (NEBN%LSTATNW) THEN

      DO JIJ=IIJB, IIJE
        ZFLXZ(JIJ, IKB) = TURBN%XCTV * ZFLXZ(JIJ, IKB)
      END DO

    END IF
  ELSE

     DO JIJ=IIJB, IIJE
       ZFLXZ(JIJ, IKB) = TURBN%XCTV * PPHI3(JIJ, IKB+IKL) * PLM(JIJ, IKB)   &
       * PLEPS(JIJ, IKB)                                         &
       *( PEXPL *                                                  &
       ( ZCOEFF(JIJ, IKB+2*IKL)*PTHLM(JIJ, IKB+2*IKL)             &
        +ZCOEFF(JIJ, IKB+IKL)*PTHLM(JIJ, IKB+IKL)             & 
        +ZCOEFF(JIJ, IKB)*PTHLM(JIJ, IKB)   )**2          &
       +TURBN%XIMPL *                                                  &
       ( ZCOEFF(JIJ, IKB+2*IKL)*PTHLP(JIJ, IKB+2*IKL)             &
        +ZCOEFF(JIJ, IKB+IKL)*PTHLP(JIJ, IKB+IKL)             &
        +ZCOEFF(JIJ, IKB)*PTHLP(JIJ, IKB)   )**2          &
       )
     END DO

   ENDIF
  !

  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, IKA) = ZFLXZ(JIJ, IKB)
  END DO

  !

  IF (NEBN%LSTATNW) THEN
    !wc  The variance from the budget eq should be multiplied by 2 here
    !    thl'2=2*L*LEPS*(dthl/dz**2)
    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZFLXZ(JIJ, JK) = MAX(0., 2.*ZFLXZ(JIJ, JK))
      END DO
    END DO
  ELSE
    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZFLXZ(JIJ, JK) = MAX(0., ZFLXZ(JIJ, JK))
      END DO
    END DO
  END IF

  !

  IF (KRRL > 0)  THEN
    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PSIGS(JIJ, JK) = ZFLXZ(JIJ, JK) * PATHETA(JIJ, JK)**2
      END DO
    END DO
  ELSE
    PSIGS(:,:) = 0.
  END IF

  !
  !
  ! stores <THl THl>
  IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
    TZFIELD = TFIELDMETADATA(        &
      CMNHNAME   = 'THL_VVAR',       &
      CSTDNAME   = '',               &
      CLONGNAME  = 'THL_VVAR',       &
      CUNITS     = 'K2',             &
      CDIR       = 'XY',             &
      CCOMMENT   = 'X_Y_Z_THL_VVAR', &
      NGRID      = 1,                &
      NTYPE      = TYPEREAL,         &
      NDIMS      = 3,                &
      LTIMEDEP   = .TRUE.            )

    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXZ)
  END IF
!
! and we store in LES configuration
!
  IF (TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    !
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZFLXZ(:,:), TLES%X_LES_SUBGRID_Thl2 )
    CALL MZF_PHY(D, PWM(:,:), ZMZF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMZF2D_WORK1(:, :)*ZFLXZ(:,:), TLES%X_LES_RES_W_SBG_Thl2 )    
CALL LES_MEAN_SUBGRID_PHY(D, TLES, -2.*CSTURB%XCTD*PSQRT_TKE(:,:)*ZFLXZ(:,:)/PLEPS(:,:), TLES%X_LES_SUBGRID_DISS_Thl2 )
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, PETHETA(:,:)*ZFLXZ(:,:), TLES%X_LES_SUBGRID_ThlThv )
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, -CSTURB%XA3*PBETA(:,:)*PETHETA(:,:)*ZFLXZ(:,:), TLES%X_LES_SUBGRID_ThlPz, .TRUE. )
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
!
  IF ( KRR /= 0 ) THEN
!
!*       4.3  <THl Rnp>
!
!
    ! Compute the turbulent variance F and F' at time t-dt.
  IF (TURBN%LHARAT) THEN

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZWORK1(JIJ, JK) = PDTH_DZ(JIJ, JK)*PDR_DZ(JIJ, JK)
      END DO
    END DO

    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    IF (NEBN%LSTATNW) THEN

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = TURBN%XCTV * &
                                  PLMF(JIJ, JK)*PLEPSF(JIJ, JK)*ZWORK2(JIJ, JK)
        END DO
      END DO

    ELSE

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = PLMF(JIJ, JK)*PLEPSF(JIJ, JK)*ZWORK2(JIJ, JK)
        END DO
      END DO

    END IF
  ELSE
    
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = 0.5*(PPHI3(JIJ, JK)+PPSI3(JIJ, JK))*PDTH_DZ(JIJ, JK)*PDR_DZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZF(JIJ, JK) = TURBN%XCTV*PLM(JIJ, JK)*PLEPS(JIJ, JK)* &
        ZMZF2D_WORK1(JIJ, JK)  
  END DO
END DO

!
ENDIF

    ZDFDDTDZ(:,:) = 0.     ! this term, because of discretization, is treated separately
    ZDFDDRDZ(:,:) = 0.     ! this term, because of discretization, is treated separately

    !
    ! Effect of 3rd order terms in temperature flux (at mass point)
    !
    ! d(w'th'2)/dz
    IF (GFTH2) THEN
      CALL M3_THR_WTH2(D,CSTURB,TURBN,PREDR1,PD,PLEPS,PSQRT_TKE,&
       & PBLL_O_E,PETHETA,PDR_DZ,ZWORK1)
      CALL D_M3_THR_WTH2_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,&
       & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ,ZWORK2)
      CALL D_M3_THR_WTH2_O_DDRDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,&
       & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,ZWORK3)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK1(JIJ, JK) * PFTH2(JIJ, JK)
          ZDFDDTDZ(JIJ, JK) = ZDFDDTDZ(JIJ, JK) + ZWORK2(JIJ, JK) &
                                            * PFTH2(JIJ, JK)
          ZDFDDRDZ(JIJ, JK) = ZDFDDRDZ(JIJ, JK) + ZWORK3(JIJ, JK) &
                                            * PFTH2(JIJ, JK)
        END DO
      END DO

    END IF
    !
    ! d(w'2th')/dz
    IF (GFWTH) THEN
      CALL MZF_PHY(D,PFWTH,ZWORK1)
      CALL M3_THR_W2TH(D,CSTURB,TURBN,PREDR1,PD,PLM,PLEPS,PTKEM,&
       & PDR_DZ,ZWORK2)
      CALL D_M3_THR_W2TH_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,&
       & PD,PLM,PLEPS,PTKEM,PBLL_O_E,PDR_DZ,PETHETA,ZWORK3)
      CALL D_M3_THR_W2TH_O_DDRDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,&
       & PD,PLM,PLEPS,PTKEM,ZWORK4)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK2(JIJ, JK) &
                                      * ZWORK1(JIJ, JK)
          ZDFDDTDZ(JIJ, JK) = ZDFDDTDZ(JIJ, JK) + ZWORK3(JIJ, JK) &
                                            * ZWORK1(JIJ, JK)
          ZDFDDRDZ(JIJ, JK) = ZDFDDRDZ(JIJ, JK) + ZWORK4(JIJ, JK) &
                                            * ZWORK1(JIJ, JK)
        END DO
      END DO

    END IF
    !
    ! d(w'r'2)/dz
    IF (GFR2) THEN
      CALL M3_THR_WR2(D,CSTURB,TURBN,PREDTH1,PD,PLEPS,PSQRT_TKE,&
       & PBLL_O_E,PEMOIST,PDTH_DZ,ZWORK1)
      CALL D_M3_THR_WR2_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,&
       & PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,ZWORK2)
      CALL D_M3_THR_WR2_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,&
       & PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,PDTH_DZ,ZWORK3)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK1(JIJ, JK) * PFR2(JIJ, JK)
          ZDFDDTDZ(JIJ, JK) = ZDFDDTDZ(JIJ, JK) + ZWORK2(JIJ, JK) &
                                            * PFR2(JIJ, JK)
          ZDFDDRDZ(JIJ, JK) = ZDFDDRDZ(JIJ, JK) + ZWORK3(JIJ, JK) &
                                            * PFR2(JIJ, JK)
        END DO
      END DO

    END IF
    !
      ! d(w'2r')/dz
    IF (GFWR) THEN
      CALL MZF_PHY(D,PFWR,ZWORK1)
      CALL M3_THR_W2R(D,CSTURB,TURBN,PREDTH1,PD,PLM,PLEPS,PTKEM,&
      & PDTH_DZ,ZWORK2)
      CALL D_M3_THR_W2R_O_DDTDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,&
      & PLM,PLEPS,PTKEM,ZWORK3)
      CALL D_M3_THR_W2R_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,&
      & PLM,PLEPS,PTKEM,PBLL_O_E,PDTH_DZ,PEMOIST,ZWORK4)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK2(JIJ, JK)*ZWORK1(JIJ, JK)
          ZDFDDTDZ(JIJ, JK) = ZDFDDTDZ(JIJ, JK) + ZWORK3(JIJ, JK) &
                                            * ZWORK1(JIJ, JK)
          ZDFDDRDZ(JIJ, JK) = ZDFDDRDZ(JIJ, JK) + ZWORK4(JIJ, JK) &
                                            * ZWORK1(JIJ, JK)
        END DO
      END DO

    END IF
    !
    ! d(w'th'r')/dz
    IF (GFTHR) THEN
      CALL M3_THR_WTHR(D,CSTURB,TURBN,PREDTH1,PREDR1,PD,PLEPS,&
      & PSQRT_TKE,ZWORK1)
      CALL D_M3_THR_WTHR_O_DDTDZ(D,CSTURB,TURBN,PREDTH1,PREDR1,&
      & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,ZWORK2)
      CALL D_M3_THR_WTHR_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,&
      & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,ZWORK3)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK1(JIJ, JK) * PFTHR(JIJ, JK)
          ZDFDDTDZ(JIJ, JK) = ZDFDDTDZ(JIJ, JK) + ZWORK2(JIJ, JK) &
                                            * PFTHR(JIJ, JK)
          ZDFDDRDZ(JIJ, JK) = ZDFDDRDZ(JIJ, JK) + ZWORK3(JIJ, JK) &
                                            * PFTHR(JIJ, JK)
        END DO
      END DO

    END IF
    !
    IF (TURBN%LHARAT) THEN
      
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = PTHLP(JIJ, JK) - PTHLM(JIJ, JK)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK1_2D, ZDZM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = PRP(JIJ, JK) - PRM(JIJ, JK, 1)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK1_2D, ZDZM2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = PTHLP(JIJ, JK) - PTHLM(JIJ, JK)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK2_2D, ZDZM2D_WORK3)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZDZM2D_WORK3(JIJ, JK) / PDZZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = PRP(JIJ, JK) - PRM(JIJ, JK, 1)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK2_2D, ZDZM2D_WORK4)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZDZM2D_WORK4(JIJ, JK) / PDZZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, JK)   = ZF(JIJ, JK)                  &
            + TURBN%XIMPL * PLMF(JIJ, JK)*PLEPSF(JIJ, JK)*0.5   &
            * (2. *PDR_DZ(JIJ, JK)  *ZDZM2D_WORK1(JIJ, JK) / PDZZ(JIJ, JK)                 &
                   + 2. *PDTH_DZ(JIJ, JK) *ZDZM2D_WORK2(JIJ, JK) / PDZZ(JIJ, JK))            &
            + TURBN%XIMPL * ZDFDDTDZ(JIJ, JK) * ZMZF2D_WORK1(JIJ, JK) &
            + TURBN%XIMPL * ZDFDDRDZ(JIJ, JK) * ZMZF2D_WORK2(JIJ, JK)      
  END DO
END DO

!
IF (NEBN%LSTATNW) THEN    

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            ZFLXZ(JIJ, JK)   = TURBN%XCTV * ZFLXZ(JIJ, JK)
          END DO
        END DO

      END IF
    ELSE
      CALL D_PHI3DTDZ_O_DDTDZ(D,CSTURB,TURBN,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,TURBN%CTURBDIM,GUSERV,ZWKPHIPSI1) 
      ! d(phi3*dthdz)/ddthdz term
      CALL D_PSI3DTDZ_O_DDTDZ(D,CSTURB,TURBN,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,TURBN%CTURBDIM,GUSERV,ZWKPHIPSI2) 
      ! d(psi3*dthdz)/ddthdz term
      CALL D_PHI3DRDZ_O_DDRDZ(D,CSTURB,TURBN,PPHI3,PREDTH1,PREDR1,PRED2TH3,PRED2THR3,TURBN%CTURBDIM,GUSERV,ZWKPHIPSI3)
      ! d(phi3*drdz )/ddrdz term
      CALL D_PSI3DRDZ_O_DDRDZ(D,CSTURB,TURBN,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,TURBN%CTURBDIM,GUSERV,ZWKPHIPSI4)
      ! d(psi3*drdz )/ddrdz term  
      !
      
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = PTHLP(JIJ, JK) - PTHLM(JIJ, JK)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK2_2D, ZDZM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = PRP(JIJ, JK) - PRM(JIJ, JK, 1)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK2_2D, ZDZM2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = (ZWKPHIPSI1(JIJ, JK)+ZWKPHIPSI2(JIJ, JK))& 
            * PDR_DZ(JIJ, JK)*ZDZM2D_WORK1(JIJ, JK) /PDZZ(JIJ, JK) &
                        + (ZWKPHIPSI3(JIJ, JK) + ZWKPHIPSI4(JIJ, JK)) & 
                        *PDTH_DZ(JIJ, JK)*ZDZM2D_WORK2(JIJ, JK)/PDZZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK3_2D(JIJ, JK) = PTHLP(JIJ, JK) - PTHLM(JIJ, JK)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK3_2D, ZDZM2D_WORK3)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZDZM2D_WORK3(JIJ, JK) / PDZZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = PRP(JIJ, JK) - PRM(JIJ, JK, 1)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK2_2D, ZDZM2D_WORK4)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZDZM2D_WORK4(JIJ, JK) / PDZZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK3)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, JK)   = ZF(JIJ, JK)                          &
            + TURBN%XIMPL * TURBN%XCTV*PLM(JIJ, JK)*PLEPS(JIJ, JK)*0.5 &
            * ZMZF2D_WORK1(JIJ, JK)    &
            + TURBN%XIMPL * ZDFDDTDZ(JIJ, JK) * ZMZF2D_WORK2(JIJ, JK)         &
            + TURBN%XIMPL * ZDFDDRDZ(JIJ, JK) * ZMZF2D_WORK3(JIJ, JK)    
  END DO
END DO

!
ENDIF
    !
    ! special case near the ground ( uncentred gradient )
    IF (TURBN%LHARAT) THEN

      DO JIJ=IIJB, IIJE
        ZFLXZ(JIJ, IKB) =                                                            & 
        (1. )                                                                                   &
        *( PEXPL *                                                                              &
           ( ZCOEFF(JIJ, IKB+2*IKL)*PTHLM(JIJ, IKB+2*IKL)             &
            +ZCOEFF(JIJ, IKB+IKL)*PTHLM(JIJ, IKB+IKL)             & 
            +ZCOEFF(JIJ, IKB)*PTHLM(JIJ, IKB))                &
          *( ZCOEFF(JIJ, IKB+2*IKL)*PRM(JIJ, IKB+2*IKL, 1)             &
            +ZCOEFF(JIJ, IKB+IKL)*PRM(JIJ, IKB+IKL, 1)             & 
            +ZCOEFF(JIJ, IKB)*PRM(JIJ, IKB, 1))                &
          +TURBN%XIMPL *                                                                              &
           ( ZCOEFF(JIJ, IKB+2*IKL)*PTHLP(JIJ, IKB+2*IKL)             &
            +ZCOEFF(JIJ, IKB+IKL)*PTHLP(JIJ, IKB+IKL)             &
            +ZCOEFF(JIJ, IKB)*PTHLP(JIJ, IKB))                &
          *( ZCOEFF(JIJ, IKB+2*IKL)*PRP(JIJ, IKB+2*IKL)             &
            +ZCOEFF(JIJ, IKB+IKL)*PRP(JIJ, IKB+IKL)             & 
            +ZCOEFF(JIJ, IKB)*PRP(JIJ, IKB))                &
         )
      END DO

    IF (NEBN%LSTATNW) THEN

      DO JIJ=IIJB, IIJE
        ZFLXZ(JIJ, IKB) = (TURBN%XCHT1 + TURBN%XCHT2) * ZFLXZ(JIJ, IKB)
      END DO

    END IF
    ELSE

      DO JIJ=IIJB, IIJE
        ZFLXZ(JIJ, IKB) =                                                            & 
        (TURBN%XCHT1 * PPHI3(JIJ, IKB+IKL) + TURBN%XCHT2 * PPSI3(JIJ, IKB+IKL))   &
        *( PEXPL *                                                                              &
           ( ZCOEFF(JIJ, IKB+2*IKL)*PTHLM(JIJ, IKB+2*IKL)             &
            +ZCOEFF(JIJ, IKB+IKL)*PTHLM(JIJ, IKB+IKL)             & 
            +ZCOEFF(JIJ, IKB)*PTHLM(JIJ, IKB))                &
          *( ZCOEFF(JIJ, IKB+2*IKL)*PRM(JIJ, IKB+2*IKL, 1)             &
            +ZCOEFF(JIJ, IKB+IKL)*PRM(JIJ, IKB+IKL, 1)             & 
            +ZCOEFF(JIJ, IKB)*PRM(JIJ, IKB, 1))                &
          +TURBN%XIMPL *                                                                              &
           ( ZCOEFF(JIJ, IKB+2*IKL)*PTHLP(JIJ, IKB+2*IKL)             &
            +ZCOEFF(JIJ, IKB+IKL)*PTHLP(JIJ, IKB+IKL)             &
            +ZCOEFF(JIJ, IKB)*PTHLP(JIJ, IKB))                &
          *( ZCOEFF(JIJ, IKB+2*IKL)*PRP(JIJ, IKB+2*IKL)             &
            +ZCOEFF(JIJ, IKB+IKL)*PRP(JIJ, IKB+IKL)             & 
            +ZCOEFF(JIJ, IKB)*PRP(JIJ, IKB))                &
         )
       END DO

    ENDIF

    !    
    DO JIJ=IIJB, IIJE
      ZFLXZ(JIJ, IKA) = ZFLXZ(JIJ, IKB)
    END DO
    !
    IF (NEBN%LSTATNW) THEN
      !wc  The variance from the budget eq should be multiplied by 2 here
      !    e.g. thl'2=2*L*LEPS*(cab)^-1 *(dthl/dz**2)
      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZFLXZ(JIJ, JK) = MIN(0., 2.*ZFLXZ(JIJ, JK))
        END DO
      END DO
    ENDIF

    IF ( KRRL > 0 ) THEN

      IF (NEBN%LSTATNW) THEN
        !wc Part of the new statistical cloud scheme set up. Normal notation so - sign
        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            PSIGS(JIJ, JK) = PSIGS(JIJ, JK) -     &
                           2. * PATHETA(JIJ, JK) * PAMOIST(JIJ, JK) * ZFLXZ(JIJ, JK)
          END DO
        END DO
      ELSE
        !  NB PATHETA is -b in Chaboureau Bechtold 2002 which explains the + sign here
        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            PSIGS(JIJ, JK) = PSIGS(JIJ, JK) +     &
                           2. * PATHETA(JIJ, JK) * PAMOIST(JIJ, JK) * ZFLXZ(JIJ, JK)
          END DO
        END DO
      ENDIF

    END IF
    ! stores <THl Rnp>
    IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
      TZFIELD = TFIELDMETADATA(             &
        CMNHNAME   = 'THLRCONS_VCOR',       &
        CSTDNAME   = '',                    &
        CLONGNAME  = 'THLRCONS_VCOR',       &
        CUNITS     = 'K kg kg-1',           &
        CDIR       = 'XY',                  &
        CCOMMENT   = 'X_Y_Z_THLRCONS_VCOR', &
        NGRID      = 1,                     &
        NTYPE      = TYPEREAL,              &
        NDIMS      = 3,                     &
        LTIMEDEP   = .TRUE.                 )

      CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXZ)
    END IF
!
! and we store in LES configuration
!
IF (TLES%LLES_CALL) THEN
      CALL SECOND_MNH(ZTIME1)
      !
      CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZFLXZ(:,:), TLES%X_LES_SUBGRID_THlRt )
      CALL MZF_PHY(D, PWM(:,:), ZMZF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMZF2D_WORK1(:, :)*ZFLXZ(:,:), TLES%X_LES_RES_W_SBG_ThlRt )      
CALL LES_MEAN_SUBGRID_PHY(D, TLES, -2.*CSTURB%XCTD*PSQRT_TKE(:,:)*ZFLXZ(:,:)/PLEPS(:,:), TLES%X_LES_SUBGRID_DISS_ThlRt )
      CALL LES_MEAN_SUBGRID_PHY(D, TLES, PETHETA(:,:)*ZFLXZ(:,:), TLES%X_LES_SUBGRID_RtThv )
      CALL LES_MEAN_SUBGRID_PHY(D, TLES, -CSTURB%XA3*PBETA(:,:)*PETHETA(:,:)*ZFLXZ(:,:), TLES%X_LES_SUBGRID_RtPz, .TRUE. )
      CALL LES_MEAN_SUBGRID_PHY(D, TLES, PEMOIST(:,:)*ZFLXZ(:,:), TLES%X_LES_SUBGRID_ThlThv , .TRUE. )
      CALL LES_MEAN_SUBGRID_PHY(D, TLES, -CSTURB%XA3*PBETA(:,:)*PEMOIST(:,:)*ZFLXZ(:,:), TLES%X_LES_SUBGRID_ThlPz, .TRUE. )
      !
      CALL SECOND_MNH(ZTIME2)
      TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
!*       4.4  <Rnp Rnp>
!
!
    ! Compute the turbulent variance F and F' at time t-dt.
IF (TURBN%LHARAT) THEN
  
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = PDR_DZ(JIJ, JK)**2
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZF(JIJ, JK) = PLMF(JIJ, JK)*PLEPSF(JIJ, JK)*ZMZF2D_WORK1(JIJ, JK)  
  END DO
END DO

!
IF (NEBN%LSTATNW) THEN

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        ZF(JIJ, JK) = TURBN%XCTV * ZF(JIJ, JK)
      END DO
    END DO

  END IF
ELSE
  
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = PPSI3(JIJ, JK)*PDR_DZ(JIJ, JK)**2
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZF(JIJ, JK) = TURBN%XCTV*PLM(JIJ, JK)*PLEPS(JIJ, JK)& 
                                    *ZMZF2D_WORK1(JIJ, JK)
  END DO
END DO

!
ENDIF

    ZDFDDRDZ(:,:) = 0.     ! this term, because of discretization, is treated separately

    !
    ! Effect of 3rd order terms in temperature flux (at mass point)
    !
    ! d(w'r'2)/dz
    IF (GFR2) THEN
      CALL M3_R2_WR2(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PLEPS,&
      & PSQRT_TKE,ZWORK1)
      CALL D_M3_R2_WR2_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,&
      & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PEMOIST,ZWORK2)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK1(JIJ, JK) * PFR2(JIJ, JK)
          ZDFDDRDZ(JIJ, JK) = ZDFDDRDZ(JIJ, JK) + ZWORK2(JIJ, JK) &
                                            * PFR2(JIJ, JK)
        END DO
      END DO

    END IF
    !
    ! d(w'2r')/dz
    IF (GFWR) THEN
      CALL MZF_PHY(D,PFWR,ZWORK1)
      CALL M3_R2_W2R(D,CSTURB,TURBN,PREDR1,PREDTH1,PD,PDR_DZ,&
      & PLM,PLEPS,PTKEM,ZWORK2)
      CALL D_M3_R2_W2R_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,&
      & PD,PLM,PLEPS,PTKEM,GUSERV,ZWORK3)
      !

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK2(JIJ, JK)*ZWORK1(JIJ, JK)
          ZDFDDRDZ(JIJ, JK) = ZDFDDRDZ(JIJ, JK) + ZWORK3(JIJ, JK) &
                                            * ZWORK1(JIJ, JK)
        END DO
      END DO

    END IF
    !
    IF (KRR/=0) THEN
      ! d(w'r'2)/dz
      IF (GFTH2) THEN
        CALL M3_R2_WTH2(D,CSTURB,TURBN,PD,PLEPS,PSQRT_TKE,&
          & PBLL_O_E,PETHETA,PDR_DZ,ZWORK1)
        CALL D_M3_R2_WTH2_O_DDRDZ(D,CSTURB,TURBN,PREDR1,&
          & PREDTH1,PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ,ZWORK2)
        !

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK1(JIJ, JK)*PFTH2(JIJ, JK) 
            ZDFDDRDZ(JIJ, JK) = ZDFDDRDZ(JIJ, JK) + ZWORK2(JIJ, JK) &
                                              * PFTH2(JIJ, JK)
          END DO
        END DO

      END IF
      !
      ! d(w'2r')/dz
      IF (GFWTH) THEN
      CALL MZF_PHY(D,PFWTH,ZWORK1)
        CALL M3_R2_W2TH(D,CSTURB,TURBN,PD,PLM,PLEPS,PTKEM,&
          & PBLL_O_E,PETHETA,PDR_DZ,ZWORK2)
        CALL D_M3_R2_W2TH_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,&
          & PD,PLM,PLEPS,PTKEM,PBLL_O_E,PETHETA,PDR_DZ,ZWORK3)
        !

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            ZF(JIJ, JK) = ZF(JIJ, JK)+ZWORK2(JIJ, JK)*ZWORK1(JIJ, JK)
            ZDFDDRDZ(JIJ, JK) = ZDFDDRDZ(JIJ, JK) + ZWORK3(JIJ, JK) &
                                              * ZWORK1(JIJ, JK)
          END DO
        END DO

      END IF
      !
      ! d(w'th'r')/dz
      IF (GFTHR) THEN
        CALL M3_R2_WTHR(D,CSTURB,TURBN,PREDTH1,PD,PLEPS,&
          & PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ,ZWORK1)
        CALL D_M3_R2_WTHR_O_DDRDZ(D,CSTURB,TURBN,PREDR1,PREDTH1,&
          & PD,PLEPS,PSQRT_TKE,PBLL_O_E,PETHETA,PDR_DZ,ZWORK2)
        !

        DO JK=1, IKT
          DO JIJ=IIJB, IIJE
            ZF(JIJ, JK) = ZF(JIJ, JK) + ZWORK1(JIJ, JK) &
                                        * PFTHR(JIJ, JK) 
            ZDFDDRDZ(JIJ, JK) = ZDFDDRDZ(JIJ, JK) + ZWORK2(JIJ, JK) &
                                              * PFTHR(JIJ, JK)
          END DO
        END DO

      END IF
    END IF
    !
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK1(JIJ, JK) = PRP(JIJ, JK) - PRM(JIJ, JK, 1)
    END DO
  END DO
  CALL DZM_PHY(D,ZWORK1,ZWORK2)
  IF (TURBN%LHARAT) THEN
    
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = PRP(JIJ, JK) - PRM(JIJ, JK, 1)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK2_2D, ZDZM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = 2.*PDR_DZ(JIJ, JK)* ZDZM2D_WORK1(JIJ, JK) / PDZZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = PRP(JIJ, JK) - PRM(JIJ, JK, 1)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK2_2D, ZDZM2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZDZM2D_WORK2(JIJ, JK) / PDZZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, JK) = ZF(JIJ, JK)                   &
              + TURBN%XIMPL * PLMF(JIJ, JK) *PLEPSF(JIJ, JK) &
                * ZMZF2D_WORK1(JIJ, JK) &
              + TURBN%XIMPL * ZDFDDRDZ(JIJ, JK) * ZMZF2D_WORK2(JIJ, JK)    
  END DO
END DO

!
IF (NEBN%LSTATNW) THEN

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZFLXZ(JIJ, JK) = TURBN%XCTV * ZFLXZ(JIJ, JK)
        END DO
      END DO

     END IF
  ELSE
    CALL D_PSI3DRDZ2_O_DDRDZ(D,CSTURB,TURBN,PPSI3,PREDR1,PREDTH1,PRED2R3,PRED2THR3,PDR_DZ,TURBN%CTURBDIM,GUSERV,ZWKPHIPSI1)
    !
    
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = PRP(JIJ, JK) - PRM(JIJ, JK, 1)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK2_2D, ZDZM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZWKPHIPSI1(JIJ, JK)*ZDZM2D_WORK1(JIJ, JK) / PDZZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = PRP(JIJ, JK) - PRM(JIJ, JK, 1)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK2_2D, ZDZM2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZDZM2D_WORK2(JIJ, JK) / PDZZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, JK) = ZF(JIJ, JK)                             &
              + TURBN%XIMPL * TURBN%XCTV*PLM(JIJ, JK) *PLEPS(JIJ, JK) &
                * ZMZF2D_WORK1(JIJ, JK) &
              + TURBN%XIMPL * ZDFDDRDZ(JIJ, JK) * ZMZF2D_WORK2(JIJ, JK)  
  END DO
END DO

!
ENDIF
    !
    ! special case near the ground ( uncentred gradient )
  IF (TURBN%LHARAT) THEN

    DO JIJ=IIJB, IIJE
      ZFLXZ(JIJ, IKB) =  PLMF(JIJ, IKB)   &
          * PLEPSF(JIJ, IKB)                                        &
      *( PEXPL *                                                  &
         ( ZCOEFF(JIJ, IKB+2*IKL)*PRM(JIJ, IKB+2*IKL, 1)             &
          +ZCOEFF(JIJ, IKB+IKL)*PRM(JIJ, IKB+IKL, 1)             & 
          +ZCOEFF(JIJ, IKB)*PRM(JIJ, IKB, 1))**2         &
        +TURBN%XIMPL *                                                  &
         ( ZCOEFF(JIJ, IKB+2*IKL)*PRP(JIJ, IKB+2*IKL)               &
          +ZCOEFF(JIJ, IKB+IKL)*PRP(JIJ, IKB+IKL)               &
          +ZCOEFF(JIJ, IKB)*PRP(JIJ, IKB))**2           &
      )
    END DO

    IF (NEBN%LSTATNW) THEN

      DO JIJ=IIJB, IIJE
        ZFLXZ(JIJ, IKB) = TURBN%XCHV * ZFLXZ(JIJ, IKB)
      END DO

    END IF 
  ELSE

    DO JIJ=IIJB, IIJE
      ZFLXZ(JIJ, IKB) = TURBN%XCHV * PPSI3(JIJ, IKB+IKL) * PLM(JIJ, IKB)   &
          * PLEPS(JIJ, IKB)                                        &
      *( PEXPL *                                                  &
         ( ZCOEFF(JIJ, IKB+2*IKL)*PRM(JIJ, IKB+2*IKL, 1)             &
          +ZCOEFF(JIJ, IKB+IKL)*PRM(JIJ, IKB+IKL, 1)             & 
          +ZCOEFF(JIJ, IKB)*PRM(JIJ, IKB, 1))**2         &
        +TURBN%XIMPL *                                                  &
         ( ZCOEFF(JIJ, IKB+2*IKL)*PRP(JIJ, IKB+2*IKL)               &
          +ZCOEFF(JIJ, IKB+IKL)*PRP(JIJ, IKB+IKL)               &
          +ZCOEFF(JIJ, IKB)*PRP(JIJ, IKB))**2           &
       ) 
    END DO

  ENDIF
    !

    DO JIJ=IIJB, IIJE
      ZFLXZ(JIJ, IKA) = ZFLXZ(JIJ, IKB)
    END DO
    IF (NEBN%LSTATNW) THEN
      !wc  The variance from the budget eq should be multiplied by 2 here
      !    thl'2=2*L*LEPS*(dthl/dz**2)
      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          ZFLXZ(JIJ, JK) = MAX(0., 2.*ZFLXZ(JIJ, JK))
        END DO
      END DO
    ENDIF

    !
    IF ( KRRL > 0 ) THEN

      DO JK=1, IKT
        DO JIJ=IIJB, IIJE
          PSIGS(JIJ, JK) = PSIGS(JIJ, JK) + PAMOIST(JIJ, JK) **2 &
                                         * ZFLXZ(JIJ, JK)
        END DO
      END DO

    END IF
    ! stores <Rnp Rnp>
    IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
      TZFIELD = TFIELDMETADATA(         &
        CMNHNAME   = 'RTOT_VVAR',       &
        CSTDNAME   = '',                &
        CLONGNAME  = 'RTOT_VVAR',       &
        CUNITS     = 'kg2 kg-2',        &
        CDIR       = 'XY',              &
        CCOMMENT   = 'X_Y_Z_RTOT_VVAR', &
        NGRID      = 1,                 &
        NTYPE      = TYPEREAL,          &
        NDIMS      = 3,                 &
        LTIMEDEP   = .TRUE.             )

      CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXZ)
    END IF
    !
    ! and we store in LES configuration
    !
    IF (TLES%LLES_CALL) THEN
      CALL SECOND_MNH(ZTIME1)
      !
      CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZFLXZ(:,:), TLES%X_LES_SUBGRID_Rt2 )
      CALL MZF_PHY(D, PWM(:,:), ZMZF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMZF2D_WORK1(:, :)*ZFLXZ(:,:), TLES%X_LES_RES_W_SBG_Rt2 )      
CALL LES_MEAN_SUBGRID_PHY(D, TLES, PEMOIST(:,:)*ZFLXZ(:,:), TLES%X_LES_SUBGRID_RtThv , .TRUE. )
      CALL LES_MEAN_SUBGRID_PHY(D, TLES, -CSTURB%XA3*PBETA(:,:)*PEMOIST(:,:)*ZFLXZ(:,:), TLES%X_LES_SUBGRID_RtPz, .TRUE. )
      CALL LES_MEAN_SUBGRID_PHY(D, TLES, -2.*CSTURB%XCTD*PSQRT_TKE(:,:)*ZFLXZ(:,:)/PLEPS(:,:), TLES%X_LES_SUBGRID_DISS_Rt2 )
      !
      CALL SECOND_MNH(ZTIME2)
      TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
    END IF
    !
  END IF  ! end if KRR ne 0
!
!
!        4.5  Vertical part of Sigma_s
!
  IF ( KRRL > 0 ) THEN
    ! Extrapolate PSIGS at the ground and at the top

    DO JIJ=IIJB, IIJE
      PSIGS(JIJ, IKA) = PSIGS(JIJ, IKB)
      PSIGS(JIJ, IKU) = PSIGS(JIJ, IKE)
    END DO
    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PSIGS(JIJ, JK) =  SQRT( MAX (PSIGS(JIJ, JK) , TURBN%XMINSIGS) )
      END DO
    END DO

  END IF

!
!        4.6  Deallocate
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TURB_VER_THERMO_CORR',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER_THERMO_CORR
END MODULE MODE_TURB_VER_THERMO_CORR
