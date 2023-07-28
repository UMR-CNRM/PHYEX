!MNH_LIC Copyright 1994-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TURB_VER_SV_FLUX
IMPLICIT NONE
CONTAINS
SUBROUTINE TURB_VER_SV_FLUX(D,CST,CSTURB,TURBN,TLES,ONOMIXLG,       &
                      KSV,KSV_LGBEG,KSV_LGEND,                      &
                      OBLOWSNOW,OFLYER,                             &
                      PEXPL,PTSTEP,TPFILE,PRSNOW,                   &
                      PDZZ,PDIRCOSZW,                               &
                      PRHODJ,PWM,                                   &
                      PSFSVM,PSFSVP,                                &
                      PSVM,                                         &
                      PTKEM,PLM,MFMOIST,PPSI_SV,                    &
                      PRSVS,PWSV                                    )
!
!
!
!!****  *TURB_VER_SV_FLUX* -compute the source terms due to the vertical turbulent
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
!!      Modifications: Mar  31, 1998 (V. Masson) splits the routine TURB_VER_SV_FLUX
!!      Modifications: Dec  01, 2000 (V. Masson) conservation of scalar emission
!!                                               from surface in 1DIM case
!!                                               when slopes are present
!!                     Jun  20, 2001 (J Stein) case of lagragian variables
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!                     Feb 2012(Y. Seity) add possibility to run with reversed
!!                                              vertical levels
!!      Modifications: July 2015 (Wim de Rooy) TURBN%LHARAT switch
!!                     Feb 2017(M. Leriche) add initialisation of ZSOURCE
!!                                   to avoid unknwon values outside physical domain
!!                                   and avoid negative values in sv tendencies
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  Wim de Rooy    06/2019: with energycascade, 50MF nog longer necessary
!  P. Wautelet 30/11/2022: compute PWSV only when needed
!  P. Wautelet 01/06/2023: fix for PWSV
!!--------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1,   ONLY: JPRB
USE MODE_SHUMAN_PHY, ONLY: DZM_PHY, MZM_PHY, MZF_PHY
USE YOMHOOK,    ONLY: LHOOK, DR_HOOK
!
USE MODD_CST,              ONLY: CST_t
USE MODD_CTURB,            ONLY: CSTURB_t
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_t
USE MODD_FIELD,            ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_IO,               ONLY: TFILEDATA
USE MODD_LES,              ONLY: TLES_t
USE MODD_PARAMETERS,       ONLY: JPVEXT_TURB, NMNHNAMELGTMAX
USE MODD_TURB_n,           ONLY: TURB_t
!
USE MODE_EMOIST,         ONLY: EMOIST
USE MODE_ETHETA,         ONLY: ETHETA
USE MODE_GRADIENT_W_PHY, ONLY: GZ_W_M_PHY
USE MODE_GRADIENT_M_PHY, ONLY: GZ_M_W_PHY
USE MODE_IO_FIELD_WRITE_PHY, ONLY: IO_FIELD_WRITE_PHY
USE MODE_TRIDIAG,        ONLY: TRIDIAG
!
USE MODI_LES_MEAN_SUBGRID_PHY
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
TYPE(TURB_t),           INTENT(IN)   :: TURBN
TYPE(TLES_t),           INTENT(INOUT):: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KSV, &
                                       KSV_LGBEG, KSV_LGEND ! number of scalar variables
LOGICAL,                INTENT(IN)   ::  ONOMIXLG     ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL,                INTENT(IN)   ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
LOGICAL,                INTENT(IN)   ::  OFLYER       ! MesoNH flyer diagnostic
REAL,                   INTENT(IN)   ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
REAL,                   INTENT(IN)   ::  PEXPL        ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PDZZ
                                                      ! Metric coefficients
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  MFMOIST       ! moist mf dual scheme

!
REAL, DIMENSION(D%NIJT,KSV), INTENT(IN)   ::  PSFSVM       ! t - deltat
!
REAL, DIMENSION(D%NIJT,KSV), INTENT(IN)   ::  PSFSVP       ! t + deltat
!
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PWM          ! vertical wind
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLM          ! Turb. mixing length
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PPSI_SV      ! Inv.Turb.Sch.for scalars
!
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(INOUT) ::  PRSVS
                            ! cumulated sources for the prognostic variables
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT)  :: PWSV        ! scalar flux
!
!*       0.2  declaration of local variables
!
!
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
       ZWORK3,ZWORK4,&! working var. for shuman operators (array syntax)
       ZMZMRHODJ
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IIJB,IIJE,IKB,IKE,IKA ! index value for the mass points of the domain 
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain
INTEGER             :: IKL
INTEGER             :: JSV          ! loop counters
INTEGER             :: JIJ,JK           ! loop
!
REAL :: ZTIME1, ZTIME2

REAL :: ZCSVP = 4.0  ! constant for scalar flux presso-correlation (RS81)
REAL :: ZCSV          !constant for the scalar flux
!
CHARACTER(LEN=NMNHNAMELGTMAX) :: YMNHNAME
REAL(KIND=JPRB)               :: ZHOOK_HANDLE
TYPE(TFIELDMETADATA)          :: TZFIELD
!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
!
IF (LHOOK) CALL DR_HOOK('TURB_VER_SV_FLUX',0,ZHOOK_HANDLE)
!
IKT=D%NKT  
IKTB=D%NKTB          
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IKA=D%NKA
IKL=D%NKL
IIJE=D%NIJE
IIJB=D%NIJB          
!
IF (TURBN%LHARAT) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZKEFF(:,:) =  PLM(:,:) * SQRT(PTKEM(:,:))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZWORK1(:,:) = PLM(:,:)*SQRT(PTKEM(:,:))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL MZM_PHY(D,ZWORK1,ZKEFF)
ENDIF
!
IF(OBLOWSNOW) THEN
! See Vionnet (PhD, 2012) for a complete discussion around the value of the Schmidt number for blowing snow variables
   ZCSV=CSTURB%XCHF/PRSNOW
ELSE
   ZCSV=CSTURB%XCHF
ENDIF
!----------------------------------------------------------------------------
!
!*       8.   SOURCES OF PASSIVE SCALAR VARIABLES
!             -----------------------------------
!
CALL MZM_PHY(D,PRHODJ,ZMZMRHODJ)
DO JSV=1,KSV
!
  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
!
! Preparation of the arguments for TRIDIAG
    IF (TURBN%LHARAT) THEN
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)  
      ZA(:,:) = -PTSTEP * ZKEFF(:,:) * ZMZMRHODJ(:,:) &
                                   / PDZZ(:,:)**2
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ELSE
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZA(:,:) = -PTSTEP*ZCSV*PPSI_SV(:,:,JSV) *   &
           ZKEFF(:,:) * ZMZMRHODJ(:,:) / PDZZ(:,:)**2
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ENDIF
  ZSOURCE(:,:) = 0.
!
! Compute the sources for the JSVth scalar variable

!* in 3DIM case, a part of the flux goes vertically, and another goes horizontally
! (in presence of slopes)
!* in 1DIM case, the part of energy released in horizontal flux
! is taken into account in the vertical part
  IF (TURBN%CTURBDIM=='3DIM') THEN
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZSOURCE(:,IKB) = (TURBN%XIMPL*PSFSVP(:,JSV) + PEXPL*PSFSVM(:,JSV)) / &
                       PDZZ(:,IKB) * PDIRCOSZW(:)                    &
                     * 0.5 * (1. + PRHODJ(:,IKA) / PRHODJ(:,IKB))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  ELSE
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZSOURCE(:,IKB) = (TURBN%XIMPL*PSFSVP(:,JSV) + PEXPL*PSFSVM(:,JSV)) / &
                       PDZZ(:,IKB) / PDIRCOSZW(:)                    &
                     * 0.5 * (1. + PRHODJ(:,IKA) / PRHODJ(:,IKB))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  END IF
  ZSOURCE(:,IKTB+1:IKTE-1) = 0.
  ZSOURCE(:,IKE) = 0.
!
! Obtention of the split JSV scalar variable at t+ deltat
  CALL TRIDIAG(D,PSVM(:,:,JSV),ZA,PTSTEP,PEXPL,TURBN%XIMPL,PRHODJ,ZSOURCE,ZRES)
!
!  Compute the equivalent tendency for the JSV scalar variable
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PRSVS(:,:,JSV)= PRSVS(:,:,JSV)+    &
                    PRHODJ(:,:)*(ZRES(:,:)-PSVM(:,:,JSV))/PTSTEP
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
  IF ( (TURBN%LTURB_FLX .AND. TPFILE%LOPENED) .OR. TLES%LLES_CALL .OR. OFLYER ) THEN
    ! Diagnostic of the cartesian vertical flux
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK1(:,:) = PLM(:,:)*SQRT(PTKEM(:,:))
    ZWORK2(:,:) = TURBN%XIMPL*ZRES(:,:) + PEXPL*PSVM(:,:,JSV)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL MZM_PHY(D,ZWORK1,ZWORK3)
    CALL DZM_PHY(D,ZWORK2,ZWORK4)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZFLXZ(:,:) = -ZCSV * PPSI_SV(:,:,JSV) * ZWORK3(:,:) & 
                                    / PDZZ(:,:) * &
                  ZWORK4(:,:)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ! surface flux
    !* in 3DIM case, a part of the flux goes vertically, and another goes horizontally
    ! (in presence of slopes)
    !* in 1DIM case, the part of energy released in horizontal flux
    ! is taken into account in the vertical part
    IF (TURBN%CTURBDIM=='3DIM') THEN
      !$mnh_expand_array(JIJ=IIJB:IIJE)
      ZFLXZ(:,IKB) = (TURBN%XIMPL*PSFSVP(:,JSV) + PEXPL*PSFSVM(:,JSV))  &
                       * PDIRCOSZW(:)  
      !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    ELSE
      !$mnh_expand_array(JIJ=IIJB:IIJE)
      ZFLXZ(:,IKB) = (TURBN%XIMPL*PSFSVP(:,JSV) + PEXPL*PSFSVM(:,JSV))  &
                       / PDIRCOSZW(:)
     !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    END IF
    ! extrapolates the flux under the ground so that the vertical average with
    ! the IKB flux gives the ground value
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZFLXZ(:,IKA) = ZFLXZ(:,IKB)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)

    IF ( OFLYER ) THEN
      DO JK=IKTB+1,IKTE-1
        !$mnh_expand_array(JIJ=IIJB:IIJE)
        PWSV(:,JK,JSV)=0.5*(ZFLXZ(:,JK)+ZFLXZ(:,JK+IKL))
        !$mnh_end_expand_array(JIJ=IIJB:IIJE)
      END DO
      !$mnh_expand_array(JIJ=IIJB:IIJE)
      PWSV(:,IKB,JSV)=0.5*(ZFLXZ(:,IKB)+ZFLXZ(:,IKB+IKL))
      PWSV(:,IKE,JSV)=PWSV(:,IKE-IKL,JSV)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    END IF
  END IF
  !
  IF (TURBN%LTURB_FLX .AND. TPFILE%LOPENED) THEN
    ! stores the JSVth vertical flux
    WRITE(YMNHNAME,'("WSV_FLX_",I3.3)') JSV
    TZFIELD = TFIELDMETADATA(                    &
      CMNHNAME   = TRIM( YMNHNAME ),             &
      CSTDNAME   = '',                           &
      CLONGNAME  = TRIM( YMNHNAME ),             &
      CUNITS     = 'SVUNIT m s-1',               &
      CDIR       = 'XY',                         &
      CCOMMENT   = 'X_Y_Z_' // TRIM( YMNHNAME ), &
      NGRID      = 4,                            &
      NTYPE      = TYPEREAL,                     &
      NDIMS      = 3,                            &
      LTIMEDEP   = .TRUE.                        )
    !
    CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZFLXZ)
  END IF
  !
  ! Storage in the LES configuration
  !
  IF (TLES%LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    !
    CALL MZF_PHY(D,ZFLXZ,ZWORK1)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK1, TLES%X_LES_SUBGRID_WSv(:,:,:,JSV) )
    !
    CALL GZ_W_M_PHY(D,PWM,PDZZ,ZWORK2)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK3(:,:) = ZWORK2(:,:) * ZWORK1(:,:)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK3, TLES%X_LES_RES_ddxa_W_SBG_UaSv(:,:,:,JSV) )
    !
    CALL GZ_M_W_PHY(D,PSVM(:,:,JSV),PDZZ,ZWORK1)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK2(:,:) = ZWORK1(:,:) * ZFLXZ(:,:)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL MZF_PHY(D,ZWORK2,ZWORK3)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK3, TLES%X_LES_RES_ddxa_Sv_SBG_UaSv(:,:,:,JSV) )
    !
    CALL MZF_PHY(D,ZFLXZ,ZWORK1)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK2(:,:) = -ZCSVP*SQRT(PTKEM(:,:))/PLM(:,:)*ZWORK1(:,:)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2, TLES%X_LES_SUBGRID_SvPz(:,:,:,JSV) )
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZWORK1(:,:) = PWM(:,:)*ZFLXZ(:,:)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL MZF_PHY(D,ZWORK1,ZWORK2)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2, TLES%X_LES_RES_W_SBG_WSv(:,:,:,JSV) )
    !
    CALL SECOND_MNH(ZTIME2)
    TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
  END IF
  !
END DO   ! end of scalar loop
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TURB_VER_SV_FLUX',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER_SV_FLUX
END MODULE MODE_TURB_VER_SV_FLUX
