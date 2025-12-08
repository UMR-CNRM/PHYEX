!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TURB_VER_DYN_FLUX
IMPLICIT NONE
CONTAINS
SUBROUTINE TURB_VER_DYN_FLUX(D,CST,CSTURB,TURBN,TLES,KSV,O2D,OFLAT, &
                      KRR,OOCEAN,OCOUPLES,                          &
                      PEXPL,PTSTEP,TPFILE,                          &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PCOSSLOPE,PSINSLOPE,                          &
                      PSEA_UCU,PSEA_VCU,                            &
                      PRHODJ,                                       &
                      PSFUM,PSFVM,                                  &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU33M,              &
                      PTHLM,PRM,PSVM,PUM,PVM,PWM,PUSLOPEM,PVSLOPEM, &
                      PTKEM,PLM,MFMOIST,PWU,PWV,                    &
                      PRUS,PRVS,PRWS,                               &
                      PDP,PTP,PSSUFL_C,PSSVFL_C,PSSUFL,PSSVFL       )
!     ###############################################################
!
!
!!****  *TURB_VER_DYN_FLUX* -compute the source terms due to the vertical turbulent
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
!!      Modifications: Mar  31, 1998 (V. Masson) splits the routine TURB_VER_DYN_FLUX
!!      Modifications: Oct  18, 2000 (J. Stein)  Bug in some computations for IKB level
!!      Modifications: Oct  18, 2000 (V. Masson) LES computations + OFLAT switch
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!      Modifications  July 2015 (Wim de Rooy) TURBN%LHARATU switch
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!      04/2016 (M.Moge) Use openACC directives to port the TURB part of Meso-NH on GPU
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Q. Rodier      17/01/2019 : cleaning : remove cyclic conditions on DP and ZA
!!      Modification   June 2019 (Wim de Rooy) 50*MF term can be removed with
!!                                             inclusion of energy cascade
!! JL Redelsperger 03/2021 : Add Ocean  & O-A Autocoupling LES Cases
!!--------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODE_SHUMAN_PHY
USE YOMHOOK,    ONLY: LHOOK, DR_HOOK, JPHOOK
!
USE MODD_CST,            ONLY: CST_t
USE MODD_CTURB,          ONLY: CSTURB_t
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_FIELD,          ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_LES,            ONLY: TLES_t
USE MODD_PARAMETERS,     ONLY: XUNDEF
USE MODD_TURB_n,         ONLY: TURB_t
!
USE MODE_GRADIENT_U_PHY, ONLY : GZ_U_UW_PHY, GX_U_M_PHY
USE MODE_GRADIENT_V_PHY, ONLY : GZ_V_VW_PHY, GY_V_M_PHY
USE MODE_GRADIENT_W_PHY, ONLY : GX_W_UW_PHY, GY_W_VW_PHY, GZ_W_M_PHY
USE MODE_GRADIENT_M_PHY, ONLY : GX_M_U_PHY, GY_M_V_PHY
USE MODE_IO_FIELD_WRITE_PHY, only: IO_FIELD_WRITE_PHY
USE MODE_ll
USE MODE_TRIDIAG_WIND,   ONLY: TRIDIAG_WIND
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
TYPE(TLES_t),           INTENT(INOUT):: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KSV           ! number of scalar variables
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  O2D          ! Logical for 2D model version (modd_conf)
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES version 
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
REAL,                   INTENT(IN)   ::  PEXPL        ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
TYPE(TFILEDATA),        INTENT(INOUT)::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PDXX, PDYY, PDZZ, PDZX, PDZY
                                                      ! Metric coefficients
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PZZ          ! altitude of flux points
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCOSSLOPE    ! cosinus of the angle
                                      ! between i and the slope vector
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSINSLOPE    ! sinus of the angle
                                      ! between i and the slope vector
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  MFMOIST      ! moist mass flux dual scheme

!
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSFUM,PSFVM !  normal momentum sfc flux
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCDUEFF     ! Cd * || u || at time t
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked
       ! to the maximum slope direction and the surface normal and the binormal
       ! at time t - dt
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSEA_UCU      ! ocean current along X
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSEA_VCU      ! ocean current along Y
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PUM,PVM,PWM, PTHLM
  ! Wind at t-Delta t
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN) ::  PRM
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVM
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PUSLOPEM     ! wind component along the
                                     ! maximum slope direction
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PVSLOPEM     ! wind component along the
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLM          ! Turb. mixing length
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PWU          ! momentum flux u'w'
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PWV          ! momentum flux v'w'
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT)   ::  PRUS, PRVS, PRWS
                            ! cumulated sources for the prognostic variables
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PDP          ! Dynamic TKE production term
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTP          ! Thermal TKE production term
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSUFL_C  ! Time evol Flux of U at sea surface (LOCEAN)
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSVFL_C  !
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSUFL   
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSVFL  !
!
!
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(D%NIJT)  :: ZDIRSINZW ! sinus of the angle
                   ! between the normal and the vertical at the surface
REAL, DIMENSION(D%NIJT):: ZCOEFS, &    ! coeff. for the implicit scheme for the wind at the surface
                          ZWORK11D,ZWORK21D,ZWORK31D,ZWORK41D,ZWORK51D,ZWORK61D
REAL, DIMENSION(D%NIJT,D%NKT)  ::  &
       ZA, &       ! under diagonal elements of the tri-diagonal matrix involved
                   ! in the temporal implicit scheme (also used to store coefficient
                   ! J in Section 5)
       ZRES, &     ! guess of the treated variable at t+ deltat when the turbu-
                   ! lence is the only source of evolution added to the ones
                   ! considered in ZSOURCE
       ZFLXZ,  &   ! vertical flux of the treated variable
       ZSOURCE,  & ! source of evolution for the treated variable
       ZKEFF,    & ! effectif diffusion coeff = LT * SQRT( TKE )
       ZWORK1,ZWORK2,&
       ZWORK3,ZWORK4,&
       ZWORK5        ! working var. for shuman operators (array syntax)
!
INTEGER             :: IIJE,IIJB,IKB,IKE,IKA,IKU ! index value for the mass points of the domain 
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain
INTEGER             :: JSV,JIJ,JK,JI,JJ          ! scalar loop counter
INTEGER             :: IKL
REAL, DIMENSION(D%NIJT)   :: ZCOEFFLXU, &
                             ZCOEFFLXV, ZUSLOPEM, ZVSLOPEM, &
                             ZFLUXSFCU,ZFLUXSFCV
                                    ! coefficients for the surface flux
                                    ! evaluation and copy of PUSLOPEM and
                                    ! PVSLOPEM in local 3D arrays
!
REAL :: ZTIME1, ZTIME2, ZCMFS
TYPE(TFIELDMETADATA) :: TZFIELD
LOGICAL :: GOCEAN !Intermediate variable used to work around a Cray compiler bug (CCE 13.0.0)
!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TURB_VER_DYN_FLUX',0,ZHOOK_HANDLE)
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
GOCEAN= OOCEAN
!
!$acc kernels  present_cr(ZA,ZSOURCE)
ZA(:,:)=XUNDEF
PDP(:,:)=XUNDEF
!
ZSOURCE(:,:) = 0.
ZFLXZ(:,:) = 0.
ZCMFS = CSTURB%XCMFS
IF (TURBN%LHARAT) ZCMFS=1.
!
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZDIRSINZW(IIJB:IIJE) = SQRT(1.-PDIRCOSZW(IIJB:IIJE)**2)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
!  compute the coefficients for the uncentred gradient computation near the
!  ground
!
! With TURBN%LHARATU length scale and TKE are at half levels so remove MZM
!
IF (TURBN%LHARAT) THEN
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZKEFF(IIJB:IIJE,1:IKT) =  PLM(IIJB:IIJE,1:IKT) * SQRT(PTKEM(IIJB:IIJE,1:IKT))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
ELSE
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK1(IIJB:IIJE,1:IKT) = PLM(IIJB:IIJE,1:IKT) * SQRT(PTKEM(IIJB:IIJE,1:IKT))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  CALL MZM_PHY(D,ZWORK1,ZKEFF)
ENDIF
!
!$acc kernels
ZUSLOPEM(IIJB:IIJE)=PUSLOPEM(IIJB:IIJE)
ZVSLOPEM(IIJB:IIJE)=PVSLOPEM(IIJB:IIJE)
ZFLUXSFCU(IIJB:IIJE)=PSFUM(IIJB:IIJE)
ZFLUXSFCV(IIJB:IIJE)=PSFVM(IIJB:IIJE)
!$acc end kernels
!
!----------------------------------------------------------------------------
!
!
!*       5.   SOURCES OF U,W WIND COMPONENTS AND PARTIAL DYNAMIC PRODUCTION
!             -------------------------------------------------------------
!
!*       5.1  Source of U wind component
!
! Preparation of the arguments for TRIDIAG_WIND
!
ZA(:,:)    = -PTSTEP * ZCMFS * MXM( ZKEFF(:,:) ) * &
               MXM(MZM( PRHODJ(:,:) )) / MXM (PDZZ(:,:))**2
!
!
! Compute the source of U wind component
!
! compute the coefficient between the vertical flux and the 2 components of the
! wind following the slope
!$acc kernels
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZCOEFFLXU(IIJB:IIJE) = PCDUEFF(IIJB:IIJE) * (PDIRCOSZW(IIJB:IIJE)**2 - ZDIRSINZW(IIJB:IIJE)**2) &
                                   * PCOSSLOPE(IIJB:IIJE)
ZCOEFFLXV(IIJB:IIJE) = PCDUEFF(IIJB:IIJE) * PDIRCOSZW(IIJB:IIJE) * PSINSLOPE(IIJB:IIJE)
!
! prepare the implicit scheme coefficients for the surface flux
ZCOEFS(IIJB:IIJE)=  ZCOEFFLXU(IIJB:IIJE) * PCOSSLOPE(IIJB:IIJE) * PDIRCOSZW(IIJB:IIJE)  &
                 +ZCOEFFLXV(IIJB:IIJE) * PSINSLOPE(IIJB:IIJE)
!
! average this flux to be located at the U,W vorticity point
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
!
ZCOEFS(:)=MXM(ZCOEFS(:) / PDZZ(:,IKB) )
!
!
!$acc kernels
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=IKTB+1:IKTE-1)
ZSOURCE(IIJB:IIJE,IKTB+1:IKTE-1) = 0.
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=IKTB+1:IKTE-1)
!$acc end kernels
! ZSOURCE= sfc FLUX /DZ
! Sfx flux assumed to be in SI & at vorticity point
!
IF (GOCEAN) THEN  ! Ocean model
  ZSOURCE(:,IKE) = MXM(ZFLUXSFCU(:)/PDZZ(:,IKE)) &
       *0.5 * ( 1. + MXM(PRHODJ(:,IKU)) / MXM(PRHODJ(:,IKE)))
  !
  ! Zero flux at the ocean domain bottom
!$acc kernels
!$mnh_expand_array(JIJ=IIJB:IIJE)
  ZSOURCE(IIJB:IIJE,IKB) = 0.
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
  !
ELSE ! Atmosphere
  ! Compute the explicit tangential flux at the W point    
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE)               
  ZSOURCE(IIJB:IIJE,IKB)     =                                              &
   PTAU11M(IIJB:IIJE) * PCOSSLOPE(IIJB:IIJE) * PDIRCOSZW(IIJB:IIJE) * ZDIRSINZW(IIJB:IIJE) &
   -PTAU12M(IIJB:IIJE) * PSINSLOPE(IIJB:IIJE) * ZDIRSINZW(IIJB:IIJE)                  &
   -PTAU33M(IIJB:IIJE) * PCOSSLOPE(IIJB:IIJE) * ZDIRSINZW(IIJB:IIJE) * PDIRCOSZW(IIJB:IIJE)  
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
!
  ! add the vertical part or the surface flux at the U,W vorticity point
!
    ZSOURCE(:,IKB) =                                  &
      (   MXM( ZSOURCE(:,IKB)   / PDZZ(:,IKB) ) &
       +  MXM( ZCOEFFLXU(:) / PDZZ(:,IKB)       &
               *ZUSLOPEM(:)                           &
              -ZCOEFFLXV(:) / PDZZ(:,IKB)       &
               *ZVSLOPEM(:)                      )    &
       -  ZCOEFS(:) * PUM(:,IKB) * TURBN%XIMPL  &
      ) * 0.5 * ( 1. + MXM(PRHODJ(:,IKA)) / MXM(PRHODJ(:,IKB)) )
!
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZSOURCE(IIJB:IIJE,IKE) = 0.
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
ENDIF
!
! Obtention of the split U at t+ deltat
!
CALL TRIDIAG_WIND(D,PUM,PSEA_UCU,ZA,ZCOEFS,PTSTEP,PEXPL,TURBN%XIMPL,   &
                  MXM(PRHODJ),ZSOURCE,ZRES)
!
!  Compute the equivalent tendency for the U wind component
!
PRUS(:,:)=PRUS(:,:)+MXM(PRHODJ(:,:))*(ZRES(:,:)-PUM(:,:))/PTSTEP
!
!*       5.2  Partial TKE Dynamic Production
!
! vertical flux of the U wind component
!
ZFLXZ(:,:)     = -ZCMFS * MXM(ZKEFF(:,:)) * &
                  DZM(TURBN%XIMPL*ZRES(:,:) + PEXPL*PUM(:,:)) / MXM(PDZZ(:,:))
!
IF (GOCEAN) THEN
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZFLXZ(IIJB:IIJE,IKE+1) = ZFLUXSFCU(IIJB:IIJE)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
ELSE
  ! surface flux
  ZFLXZ(:,IKB)   =   MXM(PDZZ(:,IKB))  *                   &
    ( ZSOURCE(:,IKB)                                       &
     +ZCOEFS(:) * ZRES(:,IKB) * TURBN%XIMPL                &
    ) / 0.5 / ( 1. + MXM(PRHODJ(:,IKA)) / MXM(PRHODJ(:,IKB)) )
  !
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZFLXZ(IIJB:IIJE,IKA) = ZFLXZ(IIJB:IIJE,IKB)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
END IF
!
IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
!$acc update self(ZFLXZ)
  ! stores the U wind component vertical flux
  TZFIELD = TFIELDMETADATA(                        &
    CMNHNAME   = 'UW_VFLX',                        &
    CSTDNAME   = '',                               &
    CLONGNAME  = 'UW_VFLX',                        &
    CUNITS     = 'm2 s-2',                         &
    CDIR       = 'XY',                             &
    CCOMMENT   = 'U wind component vertical flux', &
    NGRID      = 4,                                &
    NTYPE      = TYPEREAL,                         &
    NDIMS      = 3,                                &
    LTIMEDEP   = .TRUE.                            )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXZ)
END IF
!
! first part of total momentum flux
!
!$acc kernels present_cr(PWU)
PWU(IIJB:IIJE,1:IKT) = ZFLXZ(IIJB:IIJE,1:IKT)
!$acc end kernels
!
! Contribution to the TKE dynamic production of TKE
! (computed at mass point)
!
PDP(:,:) = - MZF( MXF ( ZFLXZ * GZ_U_UW(PUM,PDZZ) )  )
!
! Special cases near surface
IF (GOCEAN) THEN
  ! evaluate the dynamic production at w(IKE) and store in PDP(IKE)
  ! before to be extrapolated in tke_eps routine
  PDP(:,IKE) = - MXF (                            &
    ZFLXZ(:,IKE-IKL) * (PUM(:,IKE)-PUM(:,IKE-IKL))  &
                           / MXM(PDZZ(:,IKE-IKL))   &
                           ) 
ELSE ! Atmosphere
  ! evaluate the dynamic production at w(IKB+KKL) in PDP(IKB)
  PDP(:,IKB) = - MXF (                              &
    ZFLXZ(:,IKB+IKL) * (PUM(:,IKB+IKL)-PUM(:,IKB))  &
                         / MXM(PDZZ(:,IKB+IKL))     &
                         )
!
END IF
!
! Storage in the LES configuration
!
IF (TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, MZF(MXF(ZFLXZ)), TLES%X_LES_SUBGRID_WU )
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, MZF(MXF(GZ_U_UW(PUM,PDZZ) &
                          & *ZFLXZ)), TLES%X_LES_RES_ddxa_U_SBG_UaU )
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZCMFS * ZKEFF, TLES%X_LES_SUBGRID_Km )
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!*       5.3  Source of W wind component
!
!
IF(TURBN%CTURBDIM=='3DIM') THEN
  ! Compute the source for the W wind component
                ! used to compute the W source at the ground
!$acc kernels present_cr(ZFLXZ)
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZFLXZ(IIJB:IIJE,IKA) = 2 * ZFLXZ(IIJB:IIJE,IKB) - ZFLXZ(IIJB:IIJE,IKB+IKL) ! extrapolation 
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
 IF (GOCEAN) THEN
!$acc kernels present_cr(ZFLXZ)
   !$mnh_expand_array(JIJ=IIJB:IIJE)
   ZFLXZ(IIJB:IIJE,IKU) = 2 * ZFLXZ(IIJB:IIJE,IKE) - ZFLXZ(IIJB:IIJE,IKE-IKL) ! extrapolation
   !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
 END IF
  !
  IF (.NOT. OFLAT) THEN
    !
    PRWS(:,:)= PRWS                                      &
                -DXF( MZM( MXM(PRHODJ) /PDXX(:,:) )  * ZFLXZ(:,:) )  &
                +DZM( PRHODJ / MZF(PDZZ ) *                &
                      MXF( MZF( ZFLXZ(:,:)*PDZX(:,:) ) / PDXX(:,:) )      &
                    )
  ELSE
    PRWS(:,:)= PRWS(:,:) -DXF( MZM( MXM(PRHODJ) /PDXX(:,:) )  * ZFLXZ(:,:) )
  END IF
  !
  ! Complete the TKE dynamical production with the W wind contribution 
  !
  ZA(:,:)=-MZF( MXF ( ZFLXZ * GX_W_UW(OFLAT,PWM,PDXX,PDZZ,PDZX) )  )
  !
  ! Special cases near surface
  IF (GOCEAN) THEN
    ! evaluate the dynamic production at w(IKE) in PDP(IKE)
  !
    ZA(:,IKE) = - MXF(ZFLXZ(:,IKE) *  DXM(PWM(:,IKE)) &
                            / (0.5*(PDXX(:,IKE-IKL)+PDXX(:,IKE))) )
  !
  ELSE !Atmosphere
    ! evaluate the dynamic production at w(IKB+IKL) in PDP(IKB)
  ZA(:,IKB) = - MXF (                               &
   ZFLXZ(:,IKB+IKL) *                               &
     ( DXM( PWM(:,IKB+IKL) )                        &
      -MXM(  (PWM(:,IKB+2*IKL   )-PWM(:,IKB+IKL))   &
              /(PDZZ(:,IKB+2*IKL)+PDZZ(:,IKB+IKL))  &
            +(PWM(:,IKB+IKL)-PWM(:,IKB  ))          &
              /(PDZZ(:,IKB+IKL)+PDZZ(:,IKB  ))      &
          )                                         &
        * PDZX(:,IKB+IKL)                           &
     ) / (0.5*(PDXX(:,IKB+IKL)+PDXX(:,IKB)))        &
                          )
  END IF
  !
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PDP(IIJB:IIJE,1:IKT)=PDP(IIJB:IIJE,1:IKT)+ZA(IIJB:IIJE,1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  !
  ! Storage in the LES configuration
  !
  IF (TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, MZF(MXF(GX_W_UW(OFLAT, PWM,PDXX,&
      PDZZ,PDZX)*ZFLXZ)), TLES%X_LES_RES_ddxa_W_SBG_UaW )
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, MXF(GX_M_U(OFLAT, PTHLM,PDXX,PDZZ,PDZX)&
      * MZF(ZFLXZ)), TLES%X_LES_RES_ddxa_Thl_SBG_UaW )
    IF (KRR>=1) THEN
      CALL LES_MEAN_SUBGRID_PHY(D, TLES, MXF(GX_U_M(OFLAT, PRM(:,:,1),PDXX,PDZZ,PDZX)&
      *MZF(ZFLXZ)),TLES%X_LES_RES_ddxa_Rt_SBG_UaW )
    END IF
    DO JSV=1,KSV
      CALL LES_MEAN_SUBGRID_PHY(D, TLES, MXF(GX_U_M(OFLAT, PSVM(:,:,JSV),PDXX,PDZZ,&
      PDZX)*MZF(ZFLXZ)),TLES%X_LES_RES_ddxa_Sv_SBG_UaW(:,:,:,JSV) )
    END DO
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
END IF
!
!----------------------------------------------------------------------------
!
!
!*       6.   SOURCES OF V,W WIND COMPONENTS AND COMPLETE 1D TKE DYNAMIC PRODUCTION 
!             -----------------------------------------------------------------
!
!*       6.1  Source of V wind component
!
! Preparation of the arguments for TRIDIAG_WIND
!!
ZA(:,:) = -PTSTEP * ZCMFS * MYM( ZKEFF ) * MYM(MZM( PRHODJ )) / MYM( PDZZ )**2
!
! Compute the source of V wind component
! compute the coefficient between the vertical flux and the 2 components of the
! wind following the slope
!$acc kernels
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZCOEFFLXU(IIJB:IIJE) = PCDUEFF(IIJB:IIJE) * (PDIRCOSZW(IIJB:IIJE)**2 - ZDIRSINZW(IIJB:IIJE)**2) &
                                   * PSINSLOPE(IIJB:IIJE)
ZCOEFFLXV(IIJB:IIJE) = PCDUEFF(IIJB:IIJE) * PDIRCOSZW(IIJB:IIJE) * PCOSSLOPE(IIJB:IIJE)

! prepare the implicit scheme coefficients for the surface flux
ZCOEFS(IIJB:IIJE)=  ZCOEFFLXU(IIJB:IIJE) * PSINSLOPE(IIJB:IIJE) * PDIRCOSZW(IIJB:IIJE)  &
               +ZCOEFFLXV(IIJB:IIJE) * PCOSSLOPE(IIJB:IIJE)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
!
! average this flux to be located at the V,W vorticity point
ZCOEFS(:)=MYM(ZCOEFS(:) / PDZZ(:,IKB) )
!
! No flux in SOURCE TERM NULL OUTSIDE BC
!$acc kernels
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=IKB+1:IKE-1)
ZSOURCE(IIJB:IIJE,IKB+1:IKE-1) = 0.
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=IKB+1:IKE-1)
!$acc end kernels
!
! Surface case
IF (GOCEAN) THEN ! Ocean case
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZCOEFFLXU(IIJB:IIJE) = PCDUEFF(IIJB:IIJE)
  ZCOEFFLXV(IIJB:IIJE) = PCDUEFF(IIJB:IIJE)
  ZCOEFS(IIJB:IIJE)=ZCOEFFLXU(IIJB:IIJE)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
  ! average this flux to be located at the U,W vorticity point
  ZSOURCE(:,IKE) = MYM(ZFLUXSFCV(:) / PDZZ(:,IKE)) &
        *0.5 * ( 1. + MYM(PRHODJ(:,IKU)) / MYM(PRHODJ(:,IKE)))
  !No flux at the ocean domain bottom
!$acc kernels present_cr(ZSOURCE)
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZSOURCE(IIJB:IIJE,IKB) = 0.
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
!
ELSE ! Atmos case
!
 ! compute the explicit tangential flux at the W point
!$acc kernels present_cr(ZSOURCE)
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZSOURCE(IIJB:IIJE,IKB) =                                                                    &
    PTAU11M(IIJB:IIJE) * PSINSLOPE(IIJB:IIJE) * PDIRCOSZW(IIJB:IIJE) * ZDIRSINZW(IIJB:IIJE)   &
   +PTAU12M(IIJB:IIJE) * PCOSSLOPE(IIJB:IIJE) * ZDIRSINZW(IIJB:IIJE)                          &
   -PTAU33M(IIJB:IIJE) * PSINSLOPE(IIJB:IIJE) * ZDIRSINZW(IIJB:IIJE) * PDIRCOSZW(IIJB:IIJE)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
!
    ZSOURCE(:,IKB) =                                      &
      (   MYM( ZSOURCE(:,IKB)   / PDZZ(:,IKB) )     &
       +  MYM( ZCOEFFLXU(:) / PDZZ(:,IKB)           &
              *ZUSLOPEM(:)                                &
              +ZCOEFFLXV(:) / PDZZ(:,IKB)           &
              *ZVSLOPEM(:)                      )         &
       - ZCOEFS(:) * PVM(:,IKB) * TURBN%XIMPL       &
      ) * 0.5 * ( 1. + MYM(PRHODJ(:,IKA)) / MYM(PRHODJ(:,IKB)) )
!
  !No flux at the atmosphere top
!$acc kernels present_cr(ZSOURCE)
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZSOURCE(IIJB:IIJE,IKE) = 0.
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
ENDIF ! End of Ocean or Atmospher Cases
! 
!  Obtention of the split V at t+ deltat 
CALL TRIDIAG_WIND(D,PVM,PSEA_VCU,ZA,ZCOEFS,PTSTEP,PEXPL,TURBN%XIMPL,  &
                  MYM(PRHODJ),ZSOURCE,ZRES)
!
! Compute the equivalent tendency for the V wind component
!
PRVS(:,:)=PRVS(:,:)+MYM(PRHODJ(:,:))*(ZRES(:,:)-PVM(:,:))/PTSTEP
!
!
!*       6.2  Complete 1D dynamic Production
!
!  vertical flux of the V wind component
!
ZFLXZ(:,:)   = -ZCMFS * MYM(ZKEFF) * &
              DZM( TURBN%XIMPL*ZRES + PEXPL*PVM ) / MYM(PDZZ)
!
IF (GOCEAN) THEN
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZFLXZ(IIJB:IIJE,IKE+1)  = ZFLUXSFCV(IIJB:IIJE)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
ELSE
ZFLXZ(:,IKB)   =   MYM(PDZZ(:,IKB))  *                       &
  ( ZSOURCE(:,IKB)                                           &
   +ZCOEFS(:) * ZRES(:,IKB) * TURBN%XIMPL                    &
  ) / 0.5 / ( 1. + MYM(PRHODJ(:,IKA)) / MYM(PRHODJ(:,IKB)) )
  !
!$acc kernels
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZFLXZ(IIJB:IIJE,IKA) = ZFLXZ(IIJB:IIJE,IKB)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels 
END IF
!
IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
  ! stores the V wind component vertical flux
!$acc update self(ZFLXZ)
  TZFIELD = TFIELDMETADATA(                        &
    CMNHNAME   = 'VW_VFLX',                        &
    CSTDNAME   = '',                               &
    CLONGNAME  = 'VW_VFLX',                        &
    CUNITS     = 'm2 s-2',                         &
    CDIR       = 'XY',                             &
    CCOMMENT   = 'V wind component vertical flux', &
    NGRID      = 4,                                &
    NTYPE      = TYPEREAL,                         &
    NDIMS      = 3,                                &
    LTIMEDEP   = .TRUE.                            )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXZ)
END IF
!
! second part of total momentum flux
!
!$acc kernels present_cr(PWV)
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PWV(IIJB:IIJE,1:IKT) = ZFLXZ(IIJB:IIJE,1:IKT)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
!
!  Contribution to the TKE dynamical production 
!    computed at mass point
!
ZA(:,:) = - MZF( MYF ( ZFLXZ * GZ_V_VW(PVM, PDZZ) ) )
!
! Special cases at surface
IF (GOCEAN) THEN
  ! evaluate the dynamic production at w(IKE) in PDP(IKE)
  ! before extrapolation done in routine tke_eps_source
  ZA(:,IKE) = - MYF (                                                  &
   ZFLXZ(:,IKE-IKL) * (PVM(:,IKE)-PVM(:,IKE-IKL))  &
                          / MYM(PDZZ(:,IKE-IKL))                   &
                          )
!
ELSE ! Atmosphere
  ! evaluate the dynamic production at w(IKB+IKL) in PDP(IKB)
ZA(:,IKB)  =                                                 &
                 - MYF (                                          &
ZFLXZ(:,IKB+IKL) * (PVM(:,IKB+IKL)-PVM(:,IKB))  &
                       / MYM(PDZZ(:,IKB+IKL))               &
                       )
END IF
!
!$acc kernels
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PDP(IIJB:IIJE,1:IKT)=PDP(IIJB:IIJE,1:IKT)+ZA(IIJB:IIJE,1:IKT)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
!
! Storage in the LES configuration
!
IF (TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, MZF(MYF(ZFLXZ)), TLES%X_LES_SUBGRID_WV )
  CALL LES_MEAN_SUBGRID_PHY(D, TLES, MZF(MYF(GZ_V_VW(PVM,PDZZ)*&
                    & ZFLXZ)), TLES%X_LES_RES_ddxa_V_SBG_UaV )
  CALL SECOND_MNH(ZTIME2)
  TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
!*       6.3  Source of W wind component
!
IF(TURBN%CTURBDIM=='3DIM') THEN
  ! Compute the source for the W wind component
  IF (GOCEAN) THEN
!$acc kernels present_cr(ZFLXZ)
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZFLXZ(IIJB:IIJE,IKE+IKL) = 2 * ZFLXZ(IIJB:IIJE,IKE) - ZFLXZ(IIJB:IIJE,IKE-IKL) ! extrapolation 
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
  ELSE
!$acc kernels present_cr(ZFLXZ)
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZFLXZ(IIJB:IIJE,IKA) = 2 * ZFLXZ(IIJB:IIJE,IKB) - ZFLXZ(IIJB:IIJE,IKB+IKL) ! extrapolation
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
!$acc end kernels
  END IF
  !
  IF (.NOT. O2D) THEN
    IF (.NOT. OFLAT) THEN
      PRWS(:,:)= PRWS(:,:)                                   &
                  -DYF( MZM( MYM(PRHODJ) /PDYY ) * ZFLXZ )   &
                  +DZM( PRHODJ / MZF(PDZZ ) *                &
                        MYF( MZF( ZFLXZ*PDZY ) / PDYY )      &
                      )
    ELSE
      PRWS(:,:)= PRWS(:,:) -DYF( MZM( MYM(PRHODJ) /PDYY ) * ZFLXZ )
    END IF
  END IF
  ! 
  ! Complete the Dynamical production with the W wind component 
  IF (.NOT. O2D) THEN
    ZA(:,:) = - MZF( MYF ( ZFLXZ(:,:) * GY_W_VW(OFLAT, PWM,PDYY,PDZZ,PDZY) )  )
    ! Special case near surface 
    IF (GOCEAN) THEN
      ! evaluate the dynamic production at w(IKE) and stored in PDP(IKE)
      ZA(:,IKE) = -MYF(ZFLXZ(:,IKE) *  DYM(PWM(:,IKE)) &
                            / (0.5*(PDYY(:,IKE-IKL)+PDYY(:,IKE))))
    ELSE ! Atmosphere
      ! evaluate the dynamic production at w(IKB+IKL) and stored in PDP(IKB)
    ZA(:,IKB) = - MYF (                                &
     ZFLXZ(:,IKB+IKL) *                                &
       ( DYM( PWM(:,IKB+IKL) )                         &
        -MYM(  (PWM(:,IKB+2*IKL)-PWM(:,IKB+IKL))       &
                /(PDZZ(:,IKB+2*IKL)+PDZZ(:,IKB+IKL))   &
              +(PWM(:,IKB+IKL)-PWM(:,IKB  ))           &
                /(PDZZ(:,IKB+IKL)+PDZZ(:,IKB  ))       &
            )                                          &
          * PDZY(:,IKB+IKL)                            &
       ) / (0.5*(PDYY(:,IKB+IKL)+PDYY(:,IKB)))         &
                            )
    !
    END IF
!
!$acc kernels
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PDP(IIJB:IIJE,1:IKT)=PDP(IIJB:IIJE,1:IKT)+ZA(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$acc end kernels
  !
  END IF
  !
  ! Storage in the LES configuration
  !
  IF (TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, MZF(MYF(GY_W_VW(OFLAT, PWM,PDYY,&
    PDZZ,PDZY)*ZFLXZ(:,:))), TLES%X_LES_RES_ddxa_W_SBG_UaW , .TRUE. )
    CALL LES_MEAN_SUBGRID_PHY(D, TLES, MYF(GY_M_V(OFLAT, PTHLM,PDYY,PDZZ,PDZY)&
    *MZF(ZFLXZ)), TLES%X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE. )
    IF (KRR>=1) THEN
      CALL LES_MEAN_SUBGRID_PHY(D, TLES, MYF(GY_V_M(OFLAT, PRM(:,:,1),PDYY,PDZZ,&
      PDZY)*MZF(ZFLXZ(:,:))),TLES%X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE. )
    END IF
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
  !
END IF
!
!
!----------------------------------------------------------------------------
!
!*       7.   DIAGNOSTIC COMPUTATION OF THE 1D <W W> VARIANCE
!             -----------------------------------------------
!
IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED .AND. TURBN%CTURBDIM == '1DIM') THEN
  ZFLXZ(:,:)= (2./3.) * PTKEM(:,:)                     &
     -ZCMFS*PLM(:,:)*SQRT(PTKEM(:,:))*GZ_W_M(PWM,PDZZ)
  ! to be tested &
  !   +XCMFB*(4./3.)*PLM(:,:,:)/SQRT(PTKEM(:,:,:))*PTP(:,:,:)
  ! stores the W variance
  TZFIELD = TFIELDMETADATA(      &
    CMNHNAME   = 'W_VVAR',       &
    CSTDNAME   = '',             &
    CLONGNAME  = 'W_VVAR',       &
    CUNITS     = 'm2 s-2',       &
    CDIR       = 'XY',           &
    CCOMMENT   = 'X_Y_Z_W_VVAR', &
    NGRID      = 1,              &
    NTYPE      = TYPEREAL,       &
    NDIMS      = 3,              &
    LTIMEDEP   = .TRUE.          )
!$acc update self(ZFLXZ)
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXZ(:,:))
END IF
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TURB_VER_DYN_FLUX',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER_DYN_FLUX
END MODULE MODE_TURB_VER_DYN_FLUX
