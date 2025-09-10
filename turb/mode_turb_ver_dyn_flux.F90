!MNH_LIC Copyright 1994-2024 CNRS, Meteo-France and Universite Paul Sabatier
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
TYPE(TFILEDATA),        INTENT(INOUT)   ::  TPFILE       ! Output file
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
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMXM2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMXM2D_WORK2
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMZM2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMXM2D_WORK3
REAL, DIMENSION(D%NIJT) ::ZSHUGRADWK1_1D
REAL, DIMENSION(D%NIJT) ::ZMXM1D_WORK1
REAL, DIMENSION(D%NIJT) ::ZMXM1D_WORK2
REAL, DIMENSION(D%NIJT) ::ZMXM1D_WORK3
REAL, DIMENSION(D%NIJT) ::ZMXM1D_WORK4
REAL, DIMENSION(D%NIJT,D%NKT) ::ZSHUGRADWK1_2D
REAL, DIMENSION(D%NIJT,D%NKT) ::ZDZM2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMZF2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZSHUGRADWK2_2D
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMXF2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZGZ_U_UW2D_WORK1
REAL, DIMENSION(D%NIJT) ::ZMXF1D_WORK1
REAL, DIMENSION(D%NIJT) ::ZSHUGRADWK2_1D
REAL, DIMENSION(D%NIJT,D%NKT) ::ZDXF2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMZF2D_WORK2
REAL, DIMENSION(D%NIJT,D%NKT) ::ZGX_W_UW2D_WORK1
REAL, DIMENSION(D%NIJT) ::ZDXM1D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMYM2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMYM2D_WORK2
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMYM2D_WORK3
REAL, DIMENSION(D%NIJT) ::ZMYM1D_WORK1
REAL, DIMENSION(D%NIJT) ::ZMYM1D_WORK2
REAL, DIMENSION(D%NIJT) ::ZMYM1D_WORK3
REAL, DIMENSION(D%NIJT) ::ZMYM1D_WORK4
REAL, DIMENSION(D%NIJT,D%NKT) ::ZMYF2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZGZ_V_VW2D_WORK1
REAL, DIMENSION(D%NIJT) ::ZMYF1D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZDYF2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZGY_W_VW2D_WORK1
REAL, DIMENSION(D%NIJT) ::ZDYM1D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZGZ_W_M2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZGX_M_U2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZGX_U_M2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZGY_M_V2D_WORK1
REAL, DIMENSION(D%NIJT,D%NKT) ::ZGY_V_M2D_WORK1
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

ZA(:,:)=XUNDEF
PDP(:,:)=XUNDEF
!
ZSOURCE(:,:) = 0.
ZFLXZ(:,:) = 0.
ZCMFS = CSTURB%XCMFS
IF (TURBN%LHARAT) ZCMFS=1.
!
DO JIJ=IIJB, IIJE
  ZDIRSINZW(JIJ) = SQRT(1.-PDIRCOSZW(JIJ)**2)
END DO

!  compute the coefficients for the uncentred gradient computation near the
!  ground
!
! With TURBN%LHARATU length scale and TKE are at half levels so remove MZM
!
IF (TURBN%LHARAT) THEN

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZKEFF(JIJ, JK) =  PLM(JIJ, JK) * SQRT(PTKEM(JIJ, JK))
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

ZUSLOPEM(IIJB:IIJE)=PUSLOPEM(IIJB:IIJE)
ZVSLOPEM(IIJB:IIJE)=PVSLOPEM(IIJB:IIJE)
ZFLUXSFCU(IIJB:IIJE)=PSFUM(IIJB:IIJE)
ZFLUXSFCV(IIJB:IIJE)=PSFVM(IIJB:IIJE)

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
CALL MXM_PHY(D, ZKEFF(:,:), ZMXM2D_WORK1)
CALL MZM_PHY(D, PRHODJ(:,:), ZMZM2D_WORK1)
CALL MXM_PHY(D, ZMZM2D_WORK1, ZMXM2D_WORK2)
CALL MXM_PHY(D, PDZZ(:,:), ZMXM2D_WORK3)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZA(JIJ, JK)    = -PTSTEP * ZCMFS * ZMXM2D_WORK1(JIJ, JK) * &
                   ZMXM2D_WORK2(JIJ, JK) / ZMXM2D_WORK3(JIJ, JK)**2
  END DO
END DO

!
!
!
! Compute the source of U wind component
!
! compute the coefficient between the vertical flux and the 2 components of the
! wind following the slope

DO JIJ=IIJB, IIJE
  ZCOEFFLXU(JIJ) = PCDUEFF(JIJ) * (PDIRCOSZW(JIJ)**2 - ZDIRSINZW(JIJ)**2) &
                                     * PCOSSLOPE(JIJ)
  ZCOEFFLXV(JIJ) = PCDUEFF(JIJ) * PDIRCOSZW(JIJ) * PSINSLOPE(JIJ)
  !
  ! prepare the implicit scheme coefficients for the surface flux
  ZCOEFS(JIJ)=  ZCOEFFLXU(JIJ) * PCOSSLOPE(JIJ) * PDIRCOSZW(JIJ)  &
                   +ZCOEFFLXV(JIJ) * PSINSLOPE(JIJ)
  !
  ! average this flux to be located at the U,W vorticity point
END DO


DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZCOEFS(JIJ) / PDZZ(JIJ, IKB)
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_1D, ZMXM1D_WORK1)
!

DO JIJ=IIJB, IIJE
  ZCOEFS(JIJ)=ZMXM1D_WORK1(JIJ)
END DO

!
!
!

DO JK=IKTB+1, IKTE-1
  DO JIJ=IIJB, IIJE
    ZSOURCE(JIJ, JK) = 0.
  END DO
END DO

! ZSOURCE= sfc FLUX /DZ
! Sfx flux assumed to be in SI & at vorticity point
!
IF (GOCEAN) THEN  ! Ocean model
  
DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZFLUXSFCU(JIJ)/PDZZ(JIJ, IKE)
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_1D, ZMXM1D_WORK1)
CALL MXM2D_PHY(D, PRHODJ(:,IKU), ZMXM1D_WORK2)
CALL MXM2D_PHY(D, PRHODJ(:,IKE), ZMXM1D_WORK3)

DO JIJ=IIJB, IIJE
  ZSOURCE(JIJ, IKE) = ZMXM1D_WORK1(JIJ) &
         *0.5 * ( 1. + ZMXM1D_WORK2(JIJ) / ZMXM1D_WORK3(JIJ))  
END DO

!
!
  ! Zero flux at the ocean domain bottom

DO JIJ=IIJB, IIJE
    ZSOURCE(JIJ, IKB) = 0.
END DO

  !
ELSE ! Atmosphere
  ! Compute the explicit tangential flux at the W point    

  DO JIJ=IIJB, IIJE
    ZSOURCE(JIJ, IKB)     =                                              &
     PTAU11M(JIJ) * PCOSSLOPE(JIJ) * PDIRCOSZW(JIJ) * ZDIRSINZW(JIJ) &
     -PTAU12M(JIJ) * PSINSLOPE(JIJ) * ZDIRSINZW(JIJ)                  &
     -PTAU33M(JIJ) * PCOSSLOPE(JIJ) * ZDIRSINZW(JIJ) * PDIRCOSZW(JIJ)  
  END DO

!
  ! add the vertical part or the surface flux at the U,W vorticity point
!
    
DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZSOURCE(JIJ, IKB)   / PDZZ(JIJ, IKB)
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_1D, ZMXM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZCOEFFLXU(JIJ) / PDZZ(JIJ, IKB)       &
                 *ZUSLOPEM(JIJ)                           &
                -ZCOEFFLXV(JIJ) / PDZZ(JIJ, IKB)       &
                 *ZVSLOPEM(JIJ)
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_1D, ZMXM1D_WORK2)
CALL MXM2D_PHY(D, PRHODJ(:,IKA), ZMXM1D_WORK3)
CALL MXM2D_PHY(D, PRHODJ(:,IKB), ZMXM1D_WORK4)

DO JIJ=IIJB, IIJE
  ZSOURCE(JIJ, IKB) =                                  &
        (   ZMXM1D_WORK1(JIJ) &
         +  ZMXM1D_WORK2(JIJ)    &
         -  ZCOEFS(JIJ) * PUM(JIJ, IKB) * TURBN%XIMPL  &
        ) * 0.5 * ( 1. + ZMXM1D_WORK3(JIJ) / ZMXM1D_WORK4(JIJ) )
END DO

!
!

  DO JIJ=IIJB, IIJE
    ZSOURCE(JIJ, IKE) = 0.
  END DO

ENDIF
!
! Obtention of the split U at t+ deltat
!
CALL MXM_PHY(D, PRHODJ(:, :), ZMXM2D_WORK1)
CALL TRIDIAG_WIND(D,PUM(:, :),ZA(:, :),ZCOEFS(:),PTSTEP,PEXPL,TURBN%XIMPL,   &
                  ZMXM2D_WORK1(:, :),ZSOURCE(:, :),ZRES(:, :))
!
!  Compute the equivalent tendency for the U wind component
!
CALL MXM_PHY(D, PRHODJ(:,:), ZMXM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PRUS(JIJ, JK)=PRUS(JIJ, JK)+ZMXM2D_WORK1(JIJ, JK)*(ZRES(JIJ, JK)-PUM(JIJ, JK))/PTSTEP
  END DO
END DO

!
!
!*       5.2  Partial TKE Dynamic Production
!
! vertical flux of the U wind component
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = TURBN%XIMPL*ZRES(JIJ, JK) + PEXPL*PUM(JIJ, JK)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK1_2D, ZDZM2D_WORK1)
CALL MXM_PHY(D, ZKEFF(:,:), ZMXM2D_WORK1)
CALL MXM_PHY(D, PDZZ(:,:), ZMXM2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, JK)     = -ZCMFS * ZMXM2D_WORK1(JIJ, JK) * &
                      ZDZM2D_WORK1(JIJ, JK) / ZMXM2D_WORK2(JIJ, JK)
  END DO
END DO

!
!
IF (GOCEAN) THEN

  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, IKE+1) = ZFLUXSFCU(JIJ)
  END DO

ELSE
  ! surface flux
  CALL MXM2D_PHY(D, PDZZ(:,IKB), ZMXM1D_WORK1)
CALL MXM2D_PHY(D, PRHODJ(:,IKA), ZMXM1D_WORK2)
CALL MXM2D_PHY(D, PRHODJ(:,IKB), ZMXM1D_WORK3)

DO JIJ=IIJB, IIJE
  ZFLXZ(JIJ, IKB)   =   ZMXM1D_WORK1(JIJ)  *                   &
      ( ZSOURCE(JIJ, IKB)                                       &
       +ZCOEFS(JIJ) * ZRES(JIJ, IKB) * TURBN%XIMPL                &
      ) / 0.5 / ( 1. + ZMXM1D_WORK2(JIJ) / ZMXM1D_WORK3(JIJ) )  
END DO

!
!

  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, IKA) = ZFLXZ(JIJ, IKB)
  END DO

END IF
!
IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN

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

PWU(IIJB:IIJE,1:IKT) = ZFLXZ(IIJB:IIJE,1:IKT)

!
! Contribution to the TKE dynamic production of TKE
! (computed at mass point)
CALL GZ_U_UW_PHY(D, PUM(:, :),PDZZ(:, :), ZGZ_U_UW2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZFLXZ(JIJ, JK) * ZGZ_U_UW2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK2_2D, ZMXF2D_WORK1)
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMXF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PDP(JIJ, JK) = - ZMZF2D_WORK1(JIJ, JK)
  END DO
END DO

!
!
! Special cases near surface
IF (GOCEAN) THEN
  ! evaluate the dynamic production at w(IKE) and store in PDP(IKE)
  ! before to be extrapolated in tke_eps routine
  
DO JIJ=IIJB, IIJE
  ZSHUGRADWK2_1D(JIJ) = PDZZ(JIJ, IKE-IKL)
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK2_1D, ZMXM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZFLXZ(JIJ, IKE-IKL) * (PUM(JIJ, IKE)-PUM(JIJ, IKE-IKL))   / ZMXM1D_WORK1(JIJ)
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_1D, ZMXF1D_WORK1)

DO JIJ=IIJB, IIJE
  PDP(JIJ, IKE) = - ZMXF1D_WORK1(JIJ) 
END DO

!
ELSE ! Atmosphere
  ! evaluate the dynamic production at w(IKB+KKL) in PDP(IKB)
  
DO JIJ=IIJB, IIJE
  ZSHUGRADWK2_1D(JIJ) = PDZZ(JIJ, IKB+IKL)
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK2_1D, ZMXM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZFLXZ(JIJ, IKB+IKL) * (PUM(JIJ, IKB+IKL)-PUM(JIJ, IKB))   / ZMXM1D_WORK1(JIJ)
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_1D, ZMXF1D_WORK1)

DO JIJ=IIJB, IIJE
  PDP(JIJ, IKB) = - ZMXF1D_WORK1(JIJ)
END DO

!
!
END IF
!
! Storage in the LES configuration
!
IF (TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MXF_PHY(D, ZFLXZ(:, :), ZMXF2D_WORK1)
CALL MZF_PHY(D, ZMXF2D_WORK1, ZMZF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMZF2D_WORK1(:, :), TLES%X_LES_SUBGRID_WU )  
CALL GZ_U_UW_PHY(D, PUM(:, :),PDZZ(:, :), ZGZ_U_UW2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZGZ_U_UW2D_WORK1(JIJ, JK) &
                              & *ZFLXZ(JIJ, JK)
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK2_2D, ZMXF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMXF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMZF2D_WORK1(:, :), TLES%X_LES_RES_ddxa_U_SBG_UaU )  
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

  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, IKA) = 2 * ZFLXZ(JIJ, IKB) - ZFLXZ(JIJ, IKB+IKL) ! extrapolation 
  END DO

 IF (GOCEAN) THEN

   DO JIJ=IIJB, IIJE
     ZFLXZ(JIJ, IKU) = 2 * ZFLXZ(JIJ, IKE) - ZFLXZ(JIJ, IKE-IKL) ! extrapolation
   END DO

 END IF
  !
  IF (.NOT. OFLAT) THEN
    !
    CALL MXM_PHY(D, PRHODJ(:, :), ZMXM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZMXM2D_WORK1(JIJ, JK) /PDXX(JIJ, JK)
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_2D, ZMZM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMZM2D_WORK1(JIJ, JK)  * ZFLXZ(JIJ, JK)
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_2D, ZDXF2D_WORK1)
CALL MZF_PHY(D, PDZZ(:, :), ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZFLXZ(JIJ, JK)*PDZX(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK2_2D, ZMZF2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMZF2D_WORK2(JIJ, JK) / PDXX(JIJ, JK)
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = PRHODJ(JIJ, JK) / ZMZF2D_WORK1(JIJ, JK) *                 ZMXF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK1_2D, ZDZM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PRWS(JIJ, JK)= PRWS(JIJ, JK)                                      &
                    -ZDXF2D_WORK1(JIJ, JK)  &
                    +ZDZM2D_WORK1(JIJ, JK)  
  END DO
END DO

!
ELSE
    CALL MXM_PHY(D, PRHODJ(:, :), ZMXM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZMXM2D_WORK1(JIJ, JK) /PDXX(JIJ, JK)
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_2D, ZMZM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMZM2D_WORK1(JIJ, JK)  * ZFLXZ(JIJ, JK)
  END DO
END DO

!
CALL DXF_PHY(D, ZSHUGRADWK1_2D, ZDXF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PRWS(JIJ, JK)= PRWS(JIJ, JK) -ZDXF2D_WORK1(JIJ, JK)  
  END DO
END DO

!
END IF
  !
  ! Complete the TKE dynamical production with the W wind contribution 
  !
  CALL GX_W_UW_PHY(D, OFLAT,PWM(:, :),PDXX(:, :),PDZZ(:, :),PDZX(:, :), ZGX_W_UW2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZFLXZ(JIJ, JK) * ZGX_W_UW2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK2_2D, ZMXF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMXF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZA(JIJ, JK)=-ZMZF2D_WORK1(JIJ, JK)  
  END DO
END DO

!
!
  ! Special cases near surface
  IF (GOCEAN) THEN
    ! evaluate the dynamic production at w(IKE) in PDP(IKE)
  !
    CALL DXM2D_PHY(D, PWM(:,IKE), ZDXM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZFLXZ(JIJ, IKE) *  ZDXM1D_WORK1(JIJ) &
                              / (0.5*(PDXX(JIJ, IKE-IKL)+PDXX(JIJ, IKE)))
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_1D, ZMXF1D_WORK1)

DO JIJ=IIJB, IIJE
  ZA(JIJ, IKE) = - ZMXF1D_WORK1(JIJ)  
END DO

!
!
  ELSE !Atmosphere
    ! evaluate the dynamic production at w(IKB+IKL) in PDP(IKB)
  
DO JIJ=IIJB, IIJE
  ZSHUGRADWK2_1D(JIJ) = PWM(JIJ, IKB+IKL)
END DO

!
CALL DXM2D_PHY(D, ZSHUGRADWK2_1D, ZDXM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = (PWM(JIJ, IKB+2*IKL)-PWM(JIJ, IKB+IKL))    /(PDZZ(JIJ, IKB+2*IKL)+PDZZ(JIJ, IKB+IKL))   +(PWM(JIJ, IKB+IKL)-PWM(JIJ, IKB))           /(PDZZ(JIJ, IKB+IKL)+PDZZ(JIJ, IKB))
END DO

!
CALL MXM2D_PHY(D, ZSHUGRADWK1_1D, ZMXM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZFLXZ(JIJ, IKB+IKL) *                                ( ZDXM1D_WORK1(JIJ)                         -ZMXM1D_WORK1(JIJ)                                           * PDZX(JIJ, IKB+IKL)                           )  / (0.5*(PDXX(JIJ, IKB+IKL)+PDXX(JIJ, IKB)))
END DO

!
CALL MXF2D_PHY(D, ZSHUGRADWK1_1D, ZMXF1D_WORK1)

DO JIJ=IIJB, IIJE
  ZA(JIJ, IKB) = - ZMXF1D_WORK1(JIJ)  
END DO

!
END IF
  !

  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PDP(JIJ, JK)=PDP(JIJ, JK)+ZA(JIJ, JK)
    END DO
  END DO

  !
  ! Storage in the LES configuration
  !
  IF (TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL GX_W_UW_PHY(D, OFLAT, PWM(:, :),PDXX(:, :),&
      PDZZ(:, :),PDZX(:, :), ZGX_W_UW2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZGX_W_UW2D_WORK1(JIJ, JK)*ZFLXZ(JIJ, JK)
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK2_2D, ZMXF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMXF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMZF2D_WORK1(:, :), TLES%X_LES_RES_ddxa_W_SBG_UaW )    
CALL GX_M_U_PHY(D, OFLAT, PTHLM(:, :),PDXX(:, :),PDZZ(:, :),PDZX(:, :), ZGX_M_U2D_WORK1)
CALL MZF_PHY(D, ZFLXZ(:, :), ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZGX_M_U2D_WORK1(JIJ, JK)&
          * ZMZF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMXF2D_WORK1(:, :), TLES%X_LES_RES_ddxa_Thl_SBG_UaW )    
IF (KRR>=1) THEN
      CALL GX_U_M_PHY(D, OFLAT, PRM(:,:,1),PDXX(:, :),PDZZ(:, :),PDZX(:, :), ZGX_U_M2D_WORK1)
CALL MZF_PHY(D, ZFLXZ(:, :), ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZGX_U_M2D_WORK1(JIJ, JK)&
          *ZMZF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMXF2D_WORK1(:, :),TLES%X_LES_RES_ddxa_Rt_SBG_UaW )    
END IF
    DO JSV=1,KSV
      CALL GX_U_M_PHY(D, OFLAT, PSVM(:,:,JSV),PDXX(:, :),PDZZ(:, :),&
      PDZX(:, :), ZGX_U_M2D_WORK1)
CALL MZF_PHY(D, ZFLXZ(:, :), ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZGX_U_M2D_WORK1(JIJ, JK)*ZMZF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MXF_PHY(D, ZSHUGRADWK1_2D, ZMXF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMXF2D_WORK1(:, :),TLES%X_LES_RES_ddxa_Sv_SBG_UaW(:,:,:,JSV) )    
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
CALL MYM_PHY(D, ZKEFF(:, :), ZMYM2D_WORK1)
CALL MZM_PHY(D, PRHODJ(:, :), ZMZM2D_WORK1)
CALL MYM_PHY(D, ZMZM2D_WORK1, ZMYM2D_WORK2)
CALL MYM_PHY(D, PDZZ(:, :), ZMYM2D_WORK3)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZA(JIJ, JK) = -PTSTEP * ZCMFS * ZMYM2D_WORK1(JIJ, JK) * ZMYM2D_WORK2(JIJ, JK) / ZMYM2D_WORK3(JIJ, JK)**2
  END DO
END DO

!
!
! Compute the source of V wind component
! compute the coefficient between the vertical flux and the 2 components of the
! wind following the slope

DO JIJ=IIJB, IIJE
  ZCOEFFLXU(JIJ) = PCDUEFF(JIJ) * (PDIRCOSZW(JIJ)**2 - ZDIRSINZW(JIJ)**2) &
                                     * PSINSLOPE(JIJ)
  ZCOEFFLXV(JIJ) = PCDUEFF(JIJ) * PDIRCOSZW(JIJ) * PCOSSLOPE(JIJ)

  ! prepare the implicit scheme coefficients for the surface flux
  ZCOEFS(JIJ)=  ZCOEFFLXU(JIJ) * PSINSLOPE(JIJ) * PDIRCOSZW(JIJ)  &
                 +ZCOEFFLXV(JIJ) * PCOSSLOPE(JIJ)
END DO

!
! average this flux to be located at the V,W vorticity point

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZCOEFS(JIJ) / PDZZ(JIJ, IKB)
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_1D, ZMYM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZCOEFS(JIJ)=ZMYM1D_WORK1(JIJ)
END DO

!
!
! No flux in SOURCE TERM NULL OUTSIDE BC

DO JK=IKB+1, IKE-1
  DO JIJ=IIJB, IIJE
    ZSOURCE(JIJ, JK) = 0.
  END DO
END DO

!
! Surface case
IF (GOCEAN) THEN ! Ocean case

  DO JIJ=IIJB, IIJE
    ZCOEFFLXU(JIJ) = PCDUEFF(JIJ)
    ZCOEFFLXV(JIJ) = PCDUEFF(JIJ)
    ZCOEFS(JIJ)=ZCOEFFLXU(JIJ)
  END DO

  ! average this flux to be located at the U,W vorticity point
  
DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZFLUXSFCV(JIJ) / PDZZ(JIJ, IKE)
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_1D, ZMYM1D_WORK1)
CALL MYM2D_PHY(D, PRHODJ(:,IKU), ZMYM1D_WORK2)
CALL MYM2D_PHY(D, PRHODJ(:,IKE), ZMYM1D_WORK3)

DO JIJ=IIJB, IIJE
  ZSOURCE(JIJ, IKE) = ZMYM1D_WORK1(JIJ) &
          *0.5 * ( 1. + ZMYM1D_WORK2(JIJ) / ZMYM1D_WORK3(JIJ))  
END DO

!
!No flux at the ocean domain bottom

  DO JIJ=IIJB, IIJE
    ZSOURCE(JIJ, IKB) = 0.
  END DO

!
ELSE ! Atmos case
!
 ! compute the explicit tangential flux at the W point

  DO JIJ=IIJB, IIJE
    ZSOURCE(JIJ, IKB) =                                                                    &
      PTAU11M(JIJ) * PSINSLOPE(JIJ) * PDIRCOSZW(JIJ) * ZDIRSINZW(JIJ)   &
     +PTAU12M(JIJ) * PCOSSLOPE(JIJ) * ZDIRSINZW(JIJ)                          &
     -PTAU33M(JIJ) * PSINSLOPE(JIJ) * ZDIRSINZW(JIJ) * PDIRCOSZW(JIJ)
  END DO

!
    
DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZSOURCE(JIJ, IKB)   / PDZZ(JIJ, IKB)
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_1D, ZMYM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZCOEFFLXU(JIJ) / PDZZ(JIJ, IKB)           &
                *ZUSLOPEM(JIJ)                                &
                +ZCOEFFLXV(JIJ) / PDZZ(JIJ, IKB)           &
                *ZVSLOPEM(JIJ)
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_1D, ZMYM1D_WORK2)
CALL MYM2D_PHY(D, PRHODJ(:,IKA), ZMYM1D_WORK3)
CALL MYM2D_PHY(D, PRHODJ(:,IKB), ZMYM1D_WORK4)

DO JIJ=IIJB, IIJE
  ZSOURCE(JIJ, IKB) =                                      &
        (   ZMYM1D_WORK1(JIJ)     &
         +  ZMYM1D_WORK2(JIJ)         &
         - ZCOEFS(JIJ) * PVM(JIJ, IKB) * TURBN%XIMPL       &
        ) * 0.5 * ( 1. + ZMYM1D_WORK3(JIJ) / ZMYM1D_WORK4(JIJ) )
END DO

!
!
  !No flux at the atmosphere top

  DO JIJ=IIJB, IIJE
    ZSOURCE(JIJ, IKE) = 0.
  END DO

ENDIF ! End of Ocean or Atmospher Cases
! 
!  Obtention of the split V at t+ deltat 
CALL MYM_PHY(D, PRHODJ(:, :), ZMYM2D_WORK1)
CALL TRIDIAG_WIND(D,PVM(:, :),ZA(:, :),ZCOEFS(:),PTSTEP,PEXPL,TURBN%XIMPL,  &
                  ZMYM2D_WORK1(:, :),ZSOURCE(:, :),ZRES(:, :))
!
! Compute the equivalent tendency for the V wind component
!
CALL MYM_PHY(D, PRHODJ(:,:), ZMYM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PRVS(JIJ, JK)=PRVS(JIJ, JK)+ZMYM2D_WORK1(JIJ, JK)*(ZRES(JIJ, JK)-PVM(JIJ, JK))/PTSTEP
  END DO
END DO

!
!
!
!*       6.2  Complete 1D dynamic Production
!
!  vertical flux of the V wind component
!
CALL MYM_PHY(D, ZKEFF(:, :), ZMYM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = TURBN%XIMPL*ZRES(JIJ, JK) + PEXPL*PVM(JIJ, JK)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK1_2D, ZDZM2D_WORK1)
CALL MYM_PHY(D, PDZZ(:, :), ZMYM2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, JK)   = -ZCMFS * ZMYM2D_WORK1(JIJ, JK) * &
                  ZDZM2D_WORK1(JIJ, JK) / ZMYM2D_WORK2(JIJ, JK)
  END DO
END DO

!
!
IF (GOCEAN) THEN

  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, IKE+1)  = ZFLUXSFCV(JIJ)
  END DO

ELSE
CALL MYM2D_PHY(D, PDZZ(:,IKB), ZMYM1D_WORK1)
CALL MYM2D_PHY(D, PRHODJ(:,IKA), ZMYM1D_WORK2)
CALL MYM2D_PHY(D, PRHODJ(:,IKB), ZMYM1D_WORK3)

DO JIJ=IIJB, IIJE
  ZFLXZ(JIJ, IKB)   =   ZMYM1D_WORK1(JIJ)  *                       &
    ( ZSOURCE(JIJ, IKB)                                           &
     +ZCOEFS(JIJ) * ZRES(JIJ, IKB) * TURBN%XIMPL                    &
    ) / 0.5 / ( 1. + ZMYM1D_WORK2(JIJ) / ZMYM1D_WORK3(JIJ) )  
END DO

!
!

  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, IKA) = ZFLXZ(JIJ, IKB)
  END DO

END IF
!
IF ( TURBN%LTURB_FLX .AND. TPFILE%LOPENED ) THEN
  ! stores the V wind component vertical flux

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

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PWV(JIJ, JK) = ZFLXZ(JIJ, JK)
  END DO
END DO

!
!  Contribution to the TKE dynamical production 
!    computed at mass point
!
CALL GZ_V_VW_PHY(D, PVM(:, :), PDZZ(:, :), ZGZ_V_VW2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZFLXZ(JIJ, JK) * ZGZ_V_VW2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK2_2D, ZMYF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMYF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZA(JIJ, JK) = - ZMZF2D_WORK1(JIJ, JK)
  END DO
END DO

!
!
! Special cases at surface
IF (GOCEAN) THEN
  ! evaluate the dynamic production at w(IKE) in PDP(IKE)
  ! before extrapolation done in routine tke_eps_source
  
DO JIJ=IIJB, IIJE
  ZSHUGRADWK2_1D(JIJ) = PDZZ(JIJ, IKE-IKL)
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK2_1D, ZMYM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZFLXZ(JIJ, IKE-IKL) * (PVM(JIJ, IKE)-PVM(JIJ, IKE-IKL))   / ZMYM1D_WORK1(JIJ)
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_1D, ZMYF1D_WORK1)

DO JIJ=IIJB, IIJE
  ZA(JIJ, IKE) = - ZMYF1D_WORK1(JIJ)
END DO

!
!
ELSE ! Atmosphere
  ! evaluate the dynamic production at w(IKB+IKL) in PDP(IKB)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK2_1D(JIJ) = PDZZ(JIJ, IKB+IKL)
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK2_1D, ZMYM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZFLXZ(JIJ, IKB+IKL) * (PVM(JIJ, IKB+IKL)-PVM(JIJ, IKB))   / ZMYM1D_WORK1(JIJ)
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_1D, ZMYF1D_WORK1)

DO JIJ=IIJB, IIJE
  ZA(JIJ, IKB)  =                                                 &
                   - ZMYF1D_WORK1(JIJ)
END DO

!
END IF
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PDP(JIJ, JK)=PDP(JIJ, JK)+ZA(JIJ, JK)
  END DO
END DO

!
! Storage in the LES configuration
!
IF (TLES%LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL MYF_PHY(D, ZFLXZ(:, :), ZMYF2D_WORK1)
CALL MZF_PHY(D, ZMYF2D_WORK1, ZMZF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMZF2D_WORK1(:, :), TLES%X_LES_SUBGRID_WV )  
CALL GZ_V_VW_PHY(D, PVM(:, :),PDZZ(:, :), ZGZ_V_VW2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZGZ_V_VW2D_WORK1(JIJ, JK)*&
                        & ZFLXZ(JIJ, JK)
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK2_2D, ZMYF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMYF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMZF2D_WORK1(:, :), TLES%X_LES_RES_ddxa_V_SBG_UaV )  
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

    DO JIJ=IIJB, IIJE
      ZFLXZ(JIJ, IKE+IKL) = 2 * ZFLXZ(JIJ, IKE) - ZFLXZ(JIJ, IKE-IKL) ! extrapolation 
    END DO

  ELSE

    DO JIJ=IIJB, IIJE
      ZFLXZ(JIJ, IKA) = 2 * ZFLXZ(JIJ, IKB) - ZFLXZ(JIJ, IKB+IKL) ! extrapolation
    END DO

  END IF
  !
  IF (.NOT. O2D) THEN
    IF (.NOT. OFLAT) THEN
      CALL MYM_PHY(D, PRHODJ(:, :), ZMYM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZMYM2D_WORK1(JIJ, JK) /PDYY(JIJ, JK)
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_2D, ZMZM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMZM2D_WORK1(JIJ, JK) * ZFLXZ(JIJ, JK)
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_2D, ZDYF2D_WORK1)
CALL MZF_PHY(D, PDZZ(:, :), ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZFLXZ(JIJ, JK)*PDZY(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK2_2D, ZMZF2D_WORK2)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMZF2D_WORK2(JIJ, JK) / PDYY(JIJ, JK)
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = PRHODJ(JIJ, JK) / ZMZF2D_WORK1(JIJ, JK) *                 ZMYF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL DZM_PHY(D, ZSHUGRADWK1_2D, ZDZM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PRWS(JIJ, JK)= PRWS(JIJ, JK)                                   &
                      -ZDYF2D_WORK1(JIJ, JK)   &
                      +ZDZM2D_WORK1(JIJ, JK)    
  END DO
END DO

!
ELSE
      CALL MYM_PHY(D, PRHODJ(:, :), ZMYM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZMYM2D_WORK1(JIJ, JK) /PDYY(JIJ, JK)
  END DO
END DO

!
CALL MZM_PHY(D, ZSHUGRADWK2_2D, ZMZM2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMZM2D_WORK1(JIJ, JK) * ZFLXZ(JIJ, JK)
  END DO
END DO

!
CALL DYF_PHY(D, ZSHUGRADWK1_2D, ZDYF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    PRWS(JIJ, JK)= PRWS(JIJ, JK) -ZDYF2D_WORK1(JIJ, JK)    
  END DO
END DO

!
END IF
  END IF
  ! 
  ! Complete the Dynamical production with the W wind component 
  IF (.NOT. O2D) THEN
    CALL GY_W_VW_PHY(D, OFLAT, PWM(:, :),PDYY(:, :),PDZZ(:, :),PDZY(:, :), ZGY_W_VW2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZFLXZ(JIJ, JK) * ZGY_W_VW2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK2_2D, ZMYF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMYF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZA(JIJ, JK) = - ZMZF2D_WORK1(JIJ, JK)    
  END DO
END DO

!
! Special case near surface 
    IF (GOCEAN) THEN
      ! evaluate the dynamic production at w(IKE) and stored in PDP(IKE)
      CALL DYM2D_PHY(D, PWM(:,IKE), ZDYM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZFLXZ(JIJ, IKE) *  ZDYM1D_WORK1(JIJ) &
                              / (0.5*(PDYY(JIJ, IKE-IKL)+PDYY(JIJ, IKE)))
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_1D, ZMYF1D_WORK1)

DO JIJ=IIJB, IIJE
  ZA(JIJ, IKE) = -ZMYF1D_WORK1(JIJ)    
END DO

!
ELSE ! Atmosphere
      ! evaluate the dynamic production at w(IKB+IKL) and stored in PDP(IKB)
    
DO JIJ=IIJB, IIJE
  ZSHUGRADWK2_1D(JIJ) = PWM(JIJ, IKB+IKL)
END DO

!
CALL DYM2D_PHY(D, ZSHUGRADWK2_1D, ZDYM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = (PWM(JIJ, IKB+2*IKL)-PWM(JIJ, IKB+IKL))        /(PDZZ(JIJ, IKB+2*IKL)+PDZZ(JIJ, IKB+IKL))    +(PWM(JIJ, IKB+IKL)-PWM(JIJ, IKB))            /(PDZZ(JIJ, IKB+IKL)+PDZZ(JIJ, IKB))
END DO

!
CALL MYM2D_PHY(D, ZSHUGRADWK1_1D, ZMYM1D_WORK1)

DO JIJ=IIJB, IIJE
  ZSHUGRADWK1_1D(JIJ) = ZFLXZ(JIJ, IKB+IKL) *                                 ( ZDYM1D_WORK1(JIJ)                          -ZMYM1D_WORK1(JIJ)                                            * PDZY(JIJ, IKB+IKL)                            )  / (0.5*(PDYY(JIJ, IKB+IKL)+PDYY(JIJ, IKB)))
END DO

!
CALL MYF2D_PHY(D, ZSHUGRADWK1_1D, ZMYF1D_WORK1)

DO JIJ=IIJB, IIJE
  ZA(JIJ, IKB) = - ZMYF1D_WORK1(JIJ)    
END DO

!
!
    END IF
!

    DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PDP(JIJ, JK)=PDP(JIJ, JK)+ZA(JIJ, JK)
      END DO
    END DO

  !
  END IF
  !
  ! Storage in the LES configuration
  !
  IF (TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL GY_W_VW_PHY(D, OFLAT, PWM(:, :),PDYY(:, :),&
    PDZZ(:, :),PDZY(:, :), ZGY_W_VW2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK2_2D(JIJ, JK) = ZGY_W_VW2D_WORK1(JIJ, JK)*ZFLXZ(JIJ, JK)
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK2_2D, ZMYF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZMYF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZSHUGRADWK1_2D, ZMZF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMZF2D_WORK1(:, :), TLES%X_LES_RES_ddxa_W_SBG_UaW , .TRUE. )    
CALL GY_M_V_PHY(D, OFLAT, PTHLM(:, :),PDYY(:, :),PDZZ(:, :),PDZY(:, :), ZGY_M_V2D_WORK1)
CALL MZF_PHY(D, ZFLXZ(:, :), ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZGY_M_V2D_WORK1(JIJ, JK)&
        *ZMZF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMYF2D_WORK1(:, :), TLES%X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE. )    
IF (KRR>=1) THEN
      CALL GY_V_M_PHY(D, OFLAT, PRM(:,:,1),PDYY(:, :),PDZZ(:, :),&
      PDZY(:, :), ZGY_V_M2D_WORK1)
CALL MZF_PHY(D, ZFLXZ(:,:), ZMZF2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZSHUGRADWK1_2D(JIJ, JK) = ZGY_V_M2D_WORK1(JIJ, JK)*ZMZF2D_WORK1(JIJ, JK)
  END DO
END DO

!
CALL MYF_PHY(D, ZSHUGRADWK1_2D, ZMYF2D_WORK1)
CALL LES_MEAN_SUBGRID_PHY(D, TLES, ZMYF2D_WORK1(:, :),TLES%X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE. )    
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
  CALL GZ_W_M_PHY(D, PWM(:, :),PDZZ(:, :), ZGZ_W_M2D_WORK1)

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZFLXZ(JIJ, JK)= (2./3.) * PTKEM(JIJ, JK)                     &
         -ZCMFS*PLM(JIJ, JK)*SQRT(PTKEM(JIJ, JK))*ZGZ_W_M2D_WORK1(JIJ, JK)  
  END DO
END DO

!
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

  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXZ(:,:))
END IF
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TURB_VER_DYN_FLUX',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER_DYN_FLUX
END MODULE MODE_TURB_VER_DYN_FLUX
