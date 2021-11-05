!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    #################### 
     MODULE MODI_TURB_VER_SV_FLUX 
!    ####################
!
INTERFACE 
!
      SUBROUTINE TURB_VER_SV_FLUX(KKA,KKU,KKL,                      &
                      OTURB_FLX,HTURBDIM,                           &
                      PIMPL,PEXPL,                                  &
                      PTSTEP,                                       &
                      TPFILE,                                       &
                      PDZZ,PDIRCOSZW,                               &
                      PRHODJ,PWM,                                   &
                      PSFSVM,PSFSVP,                                &
                      PSVM,                                         &
                      PTKEM,PLM,PPSI_SV,                            &
                      PRSVS,PWSV                                    )
!
USE MODD_IO, ONLY: TFILEDATA
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=AR
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
CHARACTER(len=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDZZ 
                                                      ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVM       ! t - deltat 
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVP       ! t + deltat 
!
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PWM          ! vertical wind
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PPSI_SV      ! Inv.Turb.Sch.for scalars
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS
                            ! cumulated sources for the prognostic variables
REAL, DIMENSION(:,:,:,:), INTENT(OUT)  :: PWSV        ! scalar flux

!
!
END SUBROUTINE TURB_VER_SV_FLUX
!
END INTERFACE
!
END MODULE MODI_TURB_VER_SV_FLUX
!
!
!     ###############################################################
      SUBROUTINE TURB_VER_SV_FLUX(KKA,KKU,KKL,                      &
                      OTURB_FLX,HTURBDIM,                           &
                      PIMPL,PEXPL,                                  &
                      PTSTEP,                                       &
                      TPFILE,                                       &
                      PDZZ,PDIRCOSZW,                               &
                      PRHODJ,PWM,                                   &
                      PSFSVM,PSFSVP,                                &
                      PSVM,                                         &
                      PTKEM,PLM,PPSI_SV,                            &
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
!!                     Feb 2017(M. Leriche) add initialisation of ZSOURCE
!!                                   to avoid unknwon values outside physical domain
!!                                   and avoid negative values in sv tendencies
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
USE MODD_LES
USE MODD_CONF
USE MODD_NSV,            ONLY: XSVMIN, NSV_LGBEG, NSV_LGEND
USE MODD_BLOWSNOW
USE MODE_IO_FIELD_WRITE, only: IO_Field_write
!
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_GRADIENT_M
USE MODI_SHUMAN 
USE MODI_TRIDIAG 
USE MODI_TRIDIAG_WIND 
USE MODI_EMOIST
USE MODI_ETHETA
USE MODI_LES_MEAN_SUBGRID
!
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=ARO
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
CHARACTER(len=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDZZ
                                                      ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVM       ! t - deltat 
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVP       ! t + deltat 
!
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PWM          ! vertical wind
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PPSI_SV      ! Inv.Turb.Sch.for scalars
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS
                            ! cumulated sources for the prognostic variables
REAL, DIMENSION(:,:,:,:), INTENT(OUT)  :: PWSV        ! scalar flux
!
!
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(SIZE(PSVM,1),SIZE(PSVM,2),SIZE(PSVM,3))  ::  &
       ZA, &       ! under diagonal elements of the tri-diagonal matrix involved
                   ! in the temporal implicit scheme (also used to store coefficient
                   ! J in Section 5)
       ZRES, &     ! guess of the treated variable at t+ deltat when the turbu-
                   ! lence is the only source of evolution added to the ones
                   ! considered in ZSOURCE  
       ZFLXZ,  &   ! vertical flux of the treated variable
       ZSOURCE,  & ! source of evolution for the treated variable
       ZKEFF       ! effectif diffusion coeff = LT * SQRT( TKE )
INTEGER             :: IKB,IKE      ! I index values for the Beginning and End
                                    ! mass points of the domain in the 3 direct.
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain 
INTEGER             :: JSV          ! loop counters
INTEGER             :: JK           ! loop
INTEGER             :: ISV          ! number of scalar var.
!
REAL :: ZTIME1, ZTIME2

REAL :: ZCSVP = 4.0  ! constant for scalar flux presso-correlation (RS81)
REAL :: ZCSV          !constant for the scalar flux
!
TYPE(TFIELDDATA)  :: TZFIELD
!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
!
IKB=KKA+JPVEXT_TURB*KKL
IKE=KKU-JPVEXT_TURB*KKL
IKT=SIZE(PSVM,3)
IKTE =IKT-JPVEXT_TURB  
IKTB =1+JPVEXT_TURB               
!
ISV=SIZE(PSVM,4)
!
ZKEFF(:,:,:) = MZM( PLM(:,:,:) * SQRT(PTKEM(:,:,:)) )
!
IF(LBLOWSNOW) THEN
! See Vionnet (PhD, 2012) for a complete discussion around the value of the Schmidt number for blowing snow variables           
   ZCSV= XCHF/XRSNOW
ELSE
   ZCSV= XCHF
ENDIF
!----------------------------------------------------------------------------
!
!*       8.   SOURCES OF PASSIVE SCALAR VARIABLES
!             -----------------------------------
!
DO JSV=1,ISV
!
  IF (LNOMIXLG .AND. JSV >= NSV_LGBEG .AND. JSV<= NSV_LGEND) CYCLE
!
! Preparation of the arguments for TRIDIAG 
  ZA(:,:,:)    = -PTSTEP*ZCSV*PPSI_SV(:,:,:,JSV) *   &
                 ZKEFF * MZM(PRHODJ) /   &
                 PDZZ**2
  ZSOURCE(:,:,:) = 0.
!
! Compute the sources for the JSVth scalar variable

!* in 3DIM case, a part of the flux goes vertically, and another goes horizontally
! (in presence of slopes)
!* in 1DIM case, the part of energy released in horizontal flux
! is taken into account in the vertical part
  IF (HTURBDIM=='3DIM') THEN
    ZSOURCE(:,:,IKB) = (PIMPL*PSFSVP(:,:,JSV) + PEXPL*PSFSVM(:,:,JSV)) / &
                       PDZZ(:,:,IKB) * PDIRCOSZW(:,:)                    &
                     * 0.5 * (1. + PRHODJ(:,:,KKA) / PRHODJ(:,:,IKB))   
  ELSE

    ZSOURCE(:,:,IKB) = (PIMPL*PSFSVP(:,:,JSV) + PEXPL*PSFSVM(:,:,JSV)) / &
                       PDZZ(:,:,IKB) / PDIRCOSZW(:,:)                    &
                     * 0.5 * (1. + PRHODJ(:,:,KKA) / PRHODJ(:,:,IKB))
  END IF
  ZSOURCE(:,:,IKTB+1:IKTE-1) = 0.
  ZSOURCE(:,:,IKE) = 0.
!
! Obtention of the split JSV scalar variable at t+ deltat  
  CALL TRIDIAG(KKA,KKU,KKL,PSVM(:,:,:,JSV),ZA,PTSTEP,PEXPL,PIMPL,PRHODJ,ZSOURCE,ZRES)
!
!  Compute the equivalent tendency for the JSV scalar variable
  PRSVS(:,:,:,JSV)= PRSVS(:,:,:,JSV)+    &
                    PRHODJ(:,:,:)*(ZRES(:,:,:)-PSVM(:,:,:,JSV))/PTSTEP
! PRSVS(:,:,:,JSV)= MAX((PRSVS(:,:,:,JSV)+    &
!                   PRHODJ(:,:,:)*(ZRES(:,:,:)-PSVM(:,:,:,JSV))/PTSTEP),XSVMIN(JSV))
!
  IF ( (OTURB_FLX .AND. tpfile%lopened) .OR. LLES_CALL ) THEN
    ! Diagnostic of the cartesian vertical flux
    !
    ZFLXZ(:,:,:) = -ZCSV * PPSI_SV(:,:,:,JSV) * MZM(PLM*SQRT(PTKEM)) / PDZZ * &
                  DZM( PIMPL*ZRES(:,:,:) + PEXPL*PSVM(:,:,:,JSV) )
    ! surface flux
    !* in 3DIM case, a part of the flux goes vertically, and another goes horizontally
    ! (in presence of slopes)
    !* in 1DIM case, the part of energy released in horizontal flux
    ! is taken into account in the vertical part
    IF (HTURBDIM=='3DIM') THEN
      ZFLXZ(:,:,IKB) = (PIMPL*PSFSVP(:,:,JSV) + PEXPL*PSFSVM(:,:,JSV))  &
                       * PDIRCOSZW(:,:)  
    ELSE
      ZFLXZ(:,:,IKB) = (PIMPL*PSFSVP(:,:,JSV) + PEXPL*PSFSVM(:,:,JSV))  &
                       / PDIRCOSZW(:,:)
    END IF
    ! extrapolates the flux under the ground so that the vertical average with
    ! the IKB flux gives the ground value
    ZFLXZ(:,:,KKA) = ZFLXZ(:,:,IKB)
    DO JK=IKTB+1,IKTE-1
      PWSV(:,:,JK,JSV)=0.5*(ZFLXZ(:,:,JK)+ZFLXZ(:,:,JK+KKL))
    END DO
    PWSV(:,:,IKB,JSV)=0.5*(ZFLXZ(:,:,IKB)+ZFLXZ(:,:,IKB+KKL))
    PWSV(:,:,IKE,JSV)=PWSV(:,:,IKE-KKL,JSV)
 END IF
  !
  IF (OTURB_FLX .AND. tpfile%lopened) THEN
    ! stores the JSVth vertical flux
    WRITE(TZFIELD%CMNHNAME,'("WSV_FLX_",I3.3)') JSV
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    !PW: TODO: use the correct units of the JSV variable (and multiply it by m s-1)
    TZFIELD%CUNITS     = 'SVUNIT m s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
    TZFIELD%NGRID      = 4
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
  END IF
  !
  ! Storage in the LES configuration
  !
  IF (LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MZF(ZFLXZ), X_LES_SUBGRID_WSv(:,:,:,JSV) )
    CALL LES_MEAN_SUBGRID( GZ_W_M(PWM,PDZZ)*MZF(ZFLXZ), &
                           X_LES_RES_ddxa_W_SBG_UaSv(:,:,:,JSV) )
    CALL LES_MEAN_SUBGRID( MZF(GZ_M_W(KKA,KKU,KKL,PSVM(:,:,:,JSV),PDZZ)*ZFLXZ), &
                           X_LES_RES_ddxa_Sv_SBG_UaSv(:,:,:,JSV) )
    CALL LES_MEAN_SUBGRID( -ZCSVP*SQRT(PTKEM)/PLM*MZF(ZFLXZ), X_LES_SUBGRID_SvPz(:,:,:,JSV) )
    CALL LES_MEAN_SUBGRID( MZF(PWM*ZFLXZ), X_LES_RES_W_SBG_WSv(:,:,:,JSV) )
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
  !
END DO   ! end of scalar loop 
!
!----------------------------------------------------------------------------
!
END SUBROUTINE TURB_VER_SV_FLUX
