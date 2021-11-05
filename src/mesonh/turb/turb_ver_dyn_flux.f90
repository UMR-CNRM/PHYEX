!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    #################### 
     MODULE MODI_TURB_VER_DYN_FLUX 
!    ####################
!
INTERFACE 
!
      SUBROUTINE TURB_VER_DYN_FLUX(KKA,KKU,KKL,                     &
                      OTURB_FLX,KRR,                                &
                      HTURBDIM,PIMPL,PEXPL,                         &
                      PTSTEP,                                       &
                      TPFILE,                                       &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PCOSSLOPE,PSINSLOPE,                          &
                      PRHODJ,                                       &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU33M,              &
                      PTHLM,PRM,PSVM,PUM,PVM,PWM,PUSLOPEM,PVSLOPEM, &
                      PTKEM,PLM,PWU,PWV,                            &
                      PRUS,PRVS,PRWS,                               &
                      PDP,PTP                                       )
!
USE MODD_IO, ONLY: TFILEDATA
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=AR 
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
CHARACTER(len=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                      ! Metric coefficients
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PZZ          ! altitude of flux points
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCOSSLOPE    ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSINSLOPE    ! sinus of the angle 
                                      ! between i and the slope vector
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCDUEFF     ! Cd * || u || at time t
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PUM,PVM,PWM,PTHLM
  ! Wind at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRM
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PWU          ! momentum flux u'w'
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PWV          ! momentum flux v'w'
!
REAL, DIMENSION(:,:,:), INTENT(INOUT)   ::  PRUS, PRVS, PRWS
                            ! cumulated sources for the prognostic variables
!
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PDP          ! Dynamic TKE production term
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTP          ! Thermal TKE production term
!
!
!
END SUBROUTINE TURB_VER_DYN_FLUX
!
END INTERFACE
!
END MODULE MODI_TURB_VER_DYN_FLUX
!
!
!     ###############################################################
      SUBROUTINE TURB_VER_DYN_FLUX(KKA,KKU,KKL,                     &
                      OTURB_FLX,KRR,                                &
                      HTURBDIM,PIMPL,PEXPL,                         &
                      PTSTEP,                                       &
                      TPFILE,                                       &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PCOSSLOPE,PSINSLOPE,                          &
                      PRHODJ,                                       &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU33M,              &
                      PTHLM,PRM,PSVM,PUM,PVM,PWM,PUSLOPEM,PVSLOPEM, &
                      PTKEM,PLM,PWU,PWV,                            &
                      PRUS,PRVS,PRWS,                               &
                      PDP,PTP                                       )
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
!!      Modifications: Mar  31, 1998 (V. Masson) splits the routine TURB_VER_DYN_FLUX 
!!      Modifications: Oct  18, 2000 (J. Stein)  Bug in some computations for IKB level
!!      Modifications: Oct  18, 2000 (V. Masson) LES computations + LFLAT switch
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Q. Rodier      17/01/2019 : cleaning : remove cyclic conditions on DP and ZA
!! JL Redelsperger 03/2021 : Add Ocean  & O-A Autocoupling LES Cases
!!--------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_CTURB
USE MODD_DYN_n,          ONLY: LOCEAN
use modd_field,          only: tfielddata, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_LES
USE MODD_NSV
USE MODD_OCEANH
USE MODD_PARAMETERS
USE MODD_REF, ONLY : LCOUPLES
USE MODD_TURB_n
!
!
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_GRADIENT_M
USE MODI_SECOND_MNH
USE MODI_SHUMAN
USE MODI_TRIDIAG 
USE MODI_TRIDIAG_WIND 
USE MODI_LES_MEAN_SUBGRID
!
USE MODE_IO_FIELD_WRITE, only: IO_Field_write
USE MODE_ll
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
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
CHARACTER(len=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                      ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PZZ          ! altitude of flux points
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCOSSLOPE    ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSINSLOPE    ! sinus of the angle 
                                      ! between i and the slope vector
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCDUEFF     ! Cd * || u || at time t
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PUM,PVM,PWM, PTHLM
  ! Wind at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRM
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PWU          ! momentum flux u'w'
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PWV          ! momentum flux v'w'
!
REAL, DIMENSION(:,:,:), INTENT(INOUT)   ::  PRUS, PRVS, PRWS
                            ! cumulated sources for the prognostic variables
!
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PDP          ! Dynamic TKE production term
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTP          ! Thermal TKE production term
!
!
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(SIZE(PUM,1),SIZE(PUM,2))  :: ZDIRSINZW ! sinus of the angle
                   ! between the normal and the vertical at the surface
REAL, DIMENSION(SIZE(PUM,1),SIZE(PUM,2),1):: ZCOEFS    ! coeff. for the 
                   ! implicit scheme for the wind at the surface
REAL, DIMENSION(SIZE(PUM,1),SIZE(PUM,2),SIZE(PUM,3))  ::  &
       ZA, &       ! under diagonal elements of the tri-diagonal matrix involved
                   ! in the temporal implicit scheme (also used to store coefficient
                   ! J in Section 5)
       ZRES, &     ! guess of the treated variable at t+ deltat when the turbu-
                   ! lence is the only source of evolution added to the ones
                   ! considered in ZSOURCE  
       ZFLXZ,  &   ! vertical flux of the treated variable
       ZSOURCE,  & ! source of evolution for the treated variable
       ZKEFF       ! effectif diffusion coeff = LT * SQRT( TKE )
INTEGER             :: IIB,IIE, &   ! I index values for the Beginning and End
                       IJB,IJE, &   ! mass points of the domain in the 3 direct.
                       IKB,IKE      !
INTEGER             :: IKT          ! array size in k direction
INTEGER             :: IKTB,IKTE    ! start, end of k loops in physical domain
INTEGER             :: JSV          ! scalar loop counter
REAL, DIMENSION(SIZE(PDZZ,1),SIZE(PDZZ,2),1) :: ZCOEFFLXU, &
                                    ZCOEFFLXV, ZUSLOPEM, ZVSLOPEM
                                    ! coefficients for the surface flux
                                    ! evaluation and copy of PUSLOPEM and
                                    ! PVSLOPEM in local 3D arrays 
INTEGER             :: IIU,IJU      ! size of array in x,y,z directions
!
REAL :: ZTIME1, ZTIME2
TYPE(TFIELDDATA) :: TZFIELD
!----------------------------------------------------------------------------
!
!*       1.   PRELIMINARIES
!             -------------
ZA=XUNDEF
PDP=XUNDEF
!
IIU=SIZE(PUM,1)
IJU=SIZE(PUM,2)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=KKA+JPVEXT_TURB*KKL
IKE=KKU-JPVEXT_TURB*KKL
IKT=SIZE(PUM,3)          
IKTB=1+JPVEXT_TURB              
IKTE=IKT-JPVEXT_TURB


!
ZSOURCE(:,:,:) = 0.
!
ZDIRSINZW(:,:) = SQRT(1.-PDIRCOSZW(:,:)**2)
!  compute the coefficients for the uncentred gradient computation near the 
!  ground
!
ZKEFF(:,:,:) = MZM( PLM(:,:,:) * SQRT(PTKEM(:,:,:)) )
!
ZUSLOPEM(:,:,1)=PUSLOPEM(:,:)
ZVSLOPEM(:,:,1)=PVSLOPEM(:,:)
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
ZA(:,:,:)    = -PTSTEP * XCMFS *                              &
              MXM( ZKEFF ) * MXM(MZM( PRHODJ )) / &
              MXM( PDZZ )**2
!
!
! Compute the source of U wind component 
!
! compute the coefficient between the vertical flux and the 2 components of the 
! wind following the slope
ZCOEFFLXU(:,:,1) = PCDUEFF(:,:) * (PDIRCOSZW(:,:)**2 - ZDIRSINZW(:,:)**2) &
                                   * PCOSSLOPE(:,:)
ZCOEFFLXV(:,:,1) = PCDUEFF(:,:) * PDIRCOSZW(:,:) * PSINSLOPE(:,:)

! prepare the implicit scheme coefficients for the surface flux
ZCOEFS(:,:,1)=  ZCOEFFLXU(:,:,1) * PCOSSLOPE(:,:) * PDIRCOSZW(:,:)  &
                 +ZCOEFFLXV(:,:,1) * PSINSLOPE(:,:)
!
! average this flux to be located at the U,W vorticity point
ZCOEFS(:,:,1:1)=MXM(ZCOEFS(:,:,1:1) / PDZZ(:,:,IKB:IKB) )
!
!
! ZSOURCE= FLUX /DZ
IF (LOCEAN) THEN  ! OCEAN MODEL ONLY
  ! Sfx flux assumed to be in SI & at vorticity point
  IF (LCOUPLES) THEN  
    ZSOURCE(:,:,IKE:IKE) = XSSUFL_C(:,:,1:1)/PDZZ(:,:,IKE:IKE) &
         *0.5 * ( 1. + MXM(PRHODJ(:,:,KKU:KKU)) / MXM(PRHODJ(:,:,IKE:IKE))) 
  ELSE
    ZSOURCE(:,:,IKE)     = XSSUFL(:,:)
    ZSOURCE(:,:,IKE:IKE) = ZSOURCE (:,:,IKE:IKE) /PDZZ(:,:,IKE:IKE) &
        *0.5 * ( 1. + MXM(PRHODJ(:,:,KKU:KKU)) / MXM(PRHODJ(:,:,IKE:IKE)) )
  ENDIF
  !No flux at the ocean domain bottom
   ZSOURCE(:,:,IKB)           = 0.
   ZSOURCE(:,:,IKTB+1:IKTE-1) = 0
!
ELSE             !ATMOS MODEL ONLY
  IF (LCOUPLES) THEN 
   ZSOURCE(:,:,IKB:IKB) = XSSUFL_C(:,:,1:1)/PDZZ(:,:,IKB:IKB) &
      * 0.5 * ( 1. + MXM(PRHODJ(:,:,KKA:KKA)) / MXM(PRHODJ(:,:,IKB:IKB)) )
  ELSE               
    ! compute the explicit tangential flux at the W point
    ZSOURCE(:,:,IKB)     =                                              &
     PTAU11M(:,:) * PCOSSLOPE(:,:) * PDIRCOSZW(:,:) * ZDIRSINZW(:,:) &
     -PTAU12M(:,:) * PSINSLOPE(:,:) * ZDIRSINZW(:,:)                  &
     -PTAU33M(:,:) * PCOSSLOPE(:,:) * ZDIRSINZW(:,:) * PDIRCOSZW(:,:)  
!
    ! add the vertical part or the surface flux at the U,W vorticity point
!
    ZSOURCE(:,:,IKB:IKB) =                                  &
    (   MXM( ZSOURCE(:,:,IKB:IKB)   / PDZZ(:,:,IKB:IKB) ) &
    +  MXM( ZCOEFFLXU(:,:,1:1) / PDZZ(:,:,IKB:IKB)       &
           *ZUSLOPEM(:,:,1:1)                           &
          -ZCOEFFLXV(:,:,1:1) / PDZZ(:,:,IKB:IKB)       &
           *ZVSLOPEM(:,:,1:1)                      )    &
    -  ZCOEFS(:,:,1:1) * PUM(:,:,IKB:IKB) * PIMPL        &
    ) * 0.5 * ( 1. + MXM(PRHODJ(:,:,KKA:KKA)) / MXM(PRHODJ(:,:,IKB:IKB)) )
  ENDIF 
!
  ZSOURCE(:,:,IKTB+1:IKTE-1) = 0.
  ZSOURCE(:,:,IKE) = 0.
ENDIF !end ocean or atmosphere cases
!
! Obtention of the split U at t+ deltat 
!
CALL TRIDIAG_WIND(KKA,KKU,KKL,PUM,ZA,ZCOEFS(:,:,1),PTSTEP,PEXPL,PIMPL,   &
                  MXM(PRHODJ),ZSOURCE,ZRES)
! 
!  Compute the equivalent tendency for the U wind component
!
PRUS(:,:,:)=PRUS(:,:,:)+MXM(PRHODJ(:,:,:))*(ZRES(:,:,:)-PUM(:,:,:))/PTSTEP
!
!
!*       5.2  Partial Dynamic Production
!
! vertical flux of the U wind component
!
ZFLXZ(:,:,:)     = -XCMFS * MXM(ZKEFF) * &
                  DZM (PIMPL*ZRES + PEXPL*PUM) / MXM(PDZZ)
!
! surface flux 
ZFLXZ(:,:,IKB:IKB)   =   MXM(PDZZ(:,:,IKB:IKB))  *                &
  ( ZSOURCE(:,:,IKB:IKB)                                          &
   +ZCOEFS(:,:,1:1) * ZRES(:,:,IKB:IKB) * PIMPL                   &                
  ) / 0.5 / ( 1. + MXM(PRHODJ(:,:,KKA:KKA)) / MXM(PRHODJ(:,:,IKB:IKB)) )
!
ZFLXZ(:,:,KKA) = ZFLXZ(:,:,IKB) 

IF (LOCEAN) THEN !ocean model at phys sfc (ocean domain top)
  ZFLXZ(:,:,IKE:IKE)   =   MXM(PDZZ(:,:,IKE:IKE))  *                &
                           ZSOURCE(:,:,IKE:IKE)                     &
                           / 0.5 / ( 1. + MXM(PRHODJ(:,:,KKU:KKU)) / MXM(PRHODJ(:,:,IKE:IKE)) )
  ZFLXZ(:,:,KKU) = ZFLXZ(:,:,IKE) 
END IF
!
IF ( OTURB_FLX .AND. tpfile%lopened ) THEN
  ! stores the U wind component vertical flux
  TZFIELD%CMNHNAME   = 'UW_VFLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'UW_VFLX'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'U wind component vertical flux'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
END IF
!
! first part of total momentum flux
!
PWU(:,:,:) = ZFLXZ(:,:,:)
!
! Contribution to the dynamic production of TKE
! compute the dynamic production at the mass point
!
PDP(:,:,:) = - MZF( MXF ( ZFLXZ * GZ_U_UW(PUM,PDZZ) )  )
!
! evaluate the dynamic production at w(IKB+KKL) in PDP(IKB)
PDP(:,:,IKB:IKB) = - MXF (                                                      &
  ZFLXZ(:,:,IKB+KKL:IKB+KKL) * (PUM(:,:,IKB+KKL:IKB+KKL)-PUM(:,:,IKB:IKB))  &
                         / MXM(PDZZ(:,:,IKB+KKL:IKB+KKL))                   &
                         )
!
IF (LOCEAN) THEN
  ! evaluate the dynamic production at w(IKE-KKL) in PDP(IKE)
  PDP(:,:,IKE:IKE) = - MXF (                                                      &
    ZFLXZ(:,:,IKE-KKL:IKE-KKL) * (PUM(:,:,IKE:IKE)-PUM(:,:,IKE-KKL:IKE-KKL))  &
                           / MXM(PDZZ(:,:,IKE-KKL:IKE-KKL))                   &
                           ) 
END IF
!
! Storage in the LES configuration
! 
IF (LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( MZF(MXF(ZFLXZ)), X_LES_SUBGRID_WU )
  CALL LES_MEAN_SUBGRID( MZF(MXF(GZ_U_UW(PUM,PDZZ) &
                          & *ZFLXZ)), X_LES_RES_ddxa_U_SBG_UaU )
  CALL LES_MEAN_SUBGRID( XCMFS * ZKEFF, X_LES_SUBGRID_Km )
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!*       5.3  Source of W wind component
!
!
IF(HTURBDIM=='3DIM') THEN
  ! Compute the source for the W wind component
                ! used to compute the W source at the ground
  ZFLXZ(:,:,KKA) = 2 * ZFLXZ(:,:,IKB) - ZFLXZ(:,:,IKB+KKL) ! extrapolation 
 IF (LOCEAN) THEN
   ZFLXZ(:,:,KKU) = 2 * ZFLXZ(:,:,IKE) - ZFLXZ(:,:,IKE-KKL) ! extrapolation 
 END IF     
     
  !
  IF (.NOT. LFLAT) THEN
    PRWS(:,:,:)= PRWS                                      &
                -DXF( MZM( MXM(PRHODJ) /PDXX )  * ZFLXZ )  &
                +DZM( PRHODJ / MZF(PDZZ ) *                &
                      MXF( MZF( ZFLXZ*PDZX ) / PDXX )      &
                    )
  ELSE
    PRWS(:,:,:)= PRWS -DXF( MZM( MXM(PRHODJ) /PDXX )  * ZFLXZ )
  END IF
  !
  ! Complete the Dynamical production with the W wind component 
  !
  ZA(:,:,:)=-MZF( MXF ( ZFLXZ * GX_W_UW( PWM,PDXX,PDZZ,PDZX) )  )
  !
  !
  ! evaluate the dynamic production at w(IKB+KKL) in PDP(IKB)
  ZA(:,:,IKB:IKB) = - MXF (                                                  &
   ZFLXZ(:,:,IKB+KKL:IKB+KKL) *                                              &
     ( DXM( PWM(:,:,IKB+KKL:IKB+KKL) )                                       &
      -MXM(  (PWM(:,:,IKB+2*KKL:IKB+2*KKL   )-PWM(:,:,IKB+KKL:IKB+KKL))      &
              /(PDZZ(:,:,IKB+2*KKL:IKB+2*KKL)+PDZZ(:,:,IKB+KKL:IKB+KKL))     &
            +(PWM(:,:,IKB+KKL:IKB+KKL)-PWM(:,:,IKB:IKB  ))                   &
              /(PDZZ(:,:,IKB+KKL:IKB+KKL)+PDZZ(:,:,IKB:IKB  ))               &
          )                                                                  &
        * PDZX(:,:,IKB+KKL:IKB+KKL)                                          &
     ) / (0.5*(PDXX(:,:,IKB+KKL:IKB+KKL)+PDXX(:,:,IKB:IKB)))                 &
                          )
  !
IF (LOCEAN) THEN
  ! evaluate the dynamic production at w(IKE-KKL) in PDP(IKE)
  ZA(:,:,IKE:IKE) = - MXF (                                                  &
   ZFLXZ(:,:,IKE-KKL:IKE-KKL) *                                              &
     ( DXM( PWM(:,:,IKE-KKL:IKE-KKL) )                                       &
      -MXM(  (PWM(:,:,IKE-2*KKL:IKE-2*KKL   )-PWM(:,:,IKE-KKL:IKE-KKL))      &
              /(PDZZ(:,:,IKE-2*KKL:IKE-2*KKL)+PDZZ(:,:,IKE-KKL:IKE-KKL))     &
            +(PWM(:,:,IKE-KKL:IKE-KKL)-PWM(:,:,IKE:IKE  ))                   &
              /(PDZZ(:,:,IKE-KKL:IKE-KKL)+PDZZ(:,:,IKE:IKE  ))               &
          )                                                                  &
         * PDZX(:,:,IKE-KKL:IKE-KKL)                                          &
     ) / (0.5*(PDXX(:,:,IKE-KKL:IKE-KKL)+PDXX(:,:,IKE:IKE)))                 &
                          )
END IF
  !
  PDP(:,:,:)=PDP(:,:,:)+ZA(:,:,:)
  !
  ! Storage in the LES configuration
  ! 
  IF (LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MZF(MXF(GX_W_UW(PWM,PDXX,&
      PDZZ,PDZX)*ZFLXZ)), X_LES_RES_ddxa_W_SBG_UaW )
    CALL LES_MEAN_SUBGRID( MXF(GX_M_U(KKA,KKU,KKL,PTHLM,PDXX,PDZZ,PDZX)&
      * MZF(ZFLXZ)), X_LES_RES_ddxa_Thl_SBG_UaW )
    IF (KRR>=1) THEN
      CALL LES_MEAN_SUBGRID(MXF(GX_U_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX)&
      *MZF(ZFLXZ)),X_LES_RES_ddxa_Rt_SBG_UaW )
    END IF
    DO JSV=1,NSV
      CALL LES_MEAN_SUBGRID( MXF(GX_U_M(PSVM(:,:,:,JSV),PDXX,PDZZ,&
      PDZX)*MZF(ZFLXZ)),X_LES_RES_ddxa_Sv_SBG_UaW(:,:,:,JSV) )
    END DO
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
  END IF
END IF
!
!----------------------------------------------------------------------------
!
!
!*       6.   SOURCES OF V,W WIND COMPONENTS AND COMPLETE 1D DYNAMIC PRODUCTION 
!             -----------------------------------------------------------------
!
!*       6.1  Source of V wind component
!
! Preparation of the arguments for TRIDIAG_WIND
!!
ZA(:,:,:)    = - PTSTEP * XCMFS *                              &
              MYM( ZKEFF ) * MYM(MZM( PRHODJ )) /  &
              MYM( PDZZ )**2
!
!
!
! Compute the source of V wind component
! compute the coefficient between the vertical flux and the 2 components of the 
! wind following the slope
ZCOEFFLXU(:,:,1) = PCDUEFF(:,:) * (PDIRCOSZW(:,:)**2 - ZDIRSINZW(:,:)**2) &
                                   * PSINSLOPE(:,:)
ZCOEFFLXV(:,:,1) = PCDUEFF(:,:) * PDIRCOSZW(:,:) * PCOSSLOPE(:,:)

! prepare the implicit scheme coefficients for the surface flux
ZCOEFS(:,:,1)=  ZCOEFFLXU(:,:,1) * PSINSLOPE(:,:) * PDIRCOSZW(:,:)  &
               +ZCOEFFLXV(:,:,1) * PCOSSLOPE(:,:)
!
! average this flux to be located at the V,W vorticity point
ZCOEFS(:,:,1:1)=MYM(ZCOEFS(:,:,1:1) / PDZZ(:,:,IKB:IKB) )
!
IF (LOCEAN) THEN ! Ocean case
  IF (LCOUPLES) THEN
    ZSOURCE(:,:,IKE:IKE) =  XSSVFL_C(:,:,1:1)/PDZZ(:,:,IKE:IKE) &
        *0.5 * ( 1. + MYM(PRHODJ(:,:,KKU:KKU)) / MYM(PRHODJ(:,:,IKE:IKE)) ) 
  ELSE 
    ZSOURCE(:,:,IKE) = XSSVFL(:,:)
    ZSOURCE(:,:,IKE:IKE) = ZSOURCE(:,:,IKE:IKE)/PDZZ(:,:,IKE:IKE) &
        *0.5 * ( 1. + MYM(PRHODJ(:,:,KKU:KKU)) / MYM(PRHODJ(:,:,IKE:IKE)) )
  END IF
  !No flux at the ocean domain bottom
  ZSOURCE(:,:,IKB) = 0.
ELSE ! Atmos case
  IF (.NOT.LCOUPLES) THEN !  only atmosp without coupling
  ! compute the explicit tangential flux at the W point
    ZSOURCE(:,:,IKB)       =                                                  &
      PTAU11M(:,:) * PSINSLOPE(:,:) * PDIRCOSZW(:,:) * ZDIRSINZW(:,:)         &
     +PTAU12M(:,:) * PCOSSLOPE(:,:) * ZDIRSINZW(:,:)                          &
     -PTAU33M(:,:) * PSINSLOPE(:,:) * ZDIRSINZW(:,:) * PDIRCOSZW(:,:) 
!
  ! add the vertical part or the surface flux at the V,W vorticity point
  ZSOURCE(:,:,IKB:IKB) =                                      &
    (   MYM( ZSOURCE(:,:,IKB:IKB)   / PDZZ(:,:,IKB:IKB) )     &
     +  MYM( ZCOEFFLXU(:,:,1:1) / PDZZ(:,:,IKB:IKB)           &
            *ZUSLOPEM(:,:,1:1)                                &
            +ZCOEFFLXV(:,:,1:1) / PDZZ(:,:,IKB:IKB)           &
            *ZVSLOPEM(:,:,1:1)                      )         &
     - ZCOEFS(:,:,1:1) * PVM(:,:,IKB:IKB) * PIMPL             &
    ) * 0.5 * ( 1. + MYM(PRHODJ(:,:,KKA:KKA)) / MYM(PRHODJ(:,:,IKB:IKB)) )
!
  ELSE   !atmosphere when coupling
    ! input flux assumed to be in SI and at vorticity point
    ZSOURCE(:,:,IKB:IKB) =     -XSSVFL_C(:,:,1:1)/(1.*PDZZ(:,:,IKB:IKB)) &
      * 0.5 * ( 1. + MYM(PRHODJ(:,:,KKA:KKA)) / MYM(PRHODJ(:,:,IKB:IKB)) )
  ENDIF
  !No flux at the atmosphere top
  ZSOURCE(:,:,IKE) = 0.
ENDIF ! End of Ocean or Atmospher Cases
ZSOURCE(:,:,IKTB+1:IKTE-1) = 0.
! 
!  Obtention of the split V at t+ deltat 
CALL TRIDIAG_WIND(KKA,KKU,KKL,PVM,ZA,ZCOEFS(:,:,1),PTSTEP,PEXPL,PIMPL,  &
                  MYM(PRHODJ),ZSOURCE,ZRES)
!
! Compute the equivalent tendency for the V wind component
!
PRVS(:,:,:)=PRVS(:,:,:)+MYM(PRHODJ(:,:,:))*(ZRES(:,:,:)-PVM(:,:,:))/PTSTEP
!
!
!*       6.2  Complete 1D dynamic Production
!
!  vertical flux of the V wind component
!
ZFLXZ(:,:,:)   = -XCMFS * MYM(ZKEFF) * &
                DZM( PIMPL*ZRES + PEXPL*PVM ) / MYM(PDZZ)
!
ZFLXZ(:,:,IKB:IKB)   =   MYM(PDZZ(:,:,IKB:IKB))  *                       &
  ( ZSOURCE(:,:,IKB:IKB)                                                 &
   +ZCOEFS(:,:,1:1) * ZRES(:,:,IKB:IKB) * PIMPL                      &      
  ) / 0.5 / ( 1. + MYM(PRHODJ(:,:,KKA:KKA)) / MYM(PRHODJ(:,:,IKB:IKB)) )
!  
!
ZFLXZ(:,:,KKA) = ZFLXZ(:,:,IKB)
!
IF (LOCEAN) THEN
  ZFLXZ(:,:,IKE:IKE)   =   MYM(PDZZ(:,:,IKE:IKE))  *                &
      ZSOURCE(:,:,IKE:IKE)                                          &
      / 0.5 / ( 1. + MYM(PRHODJ(:,:,KKU:KKU)) / MYM(PRHODJ(:,:,IKE:IKE)) )
  ZFLXZ(:,:,KKU) = ZFLXZ(:,:,IKE) 
END IF
!
IF ( OTURB_FLX .AND. tpfile%lopened ) THEN
  ! stores the V wind component vertical flux
  TZFIELD%CMNHNAME   = 'VW_VFLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'VW_VFLX'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'V wind component vertical flux'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
END IF
!
! second part of total momentum flux
!
PWV(:,:,:) = ZFLXZ(:,:,:)
!
!  Contribution to the dynamic production of TKE
! compute the dynamic production contribution at the mass point
!
ZA(:,:,:) = - MZF( MYF ( ZFLXZ * GZ_V_VW(PVM,PDZZ) ) )
!
! evaluate the dynamic production at w(IKB+KKL) in PDP(IKB)
ZA(:,:,IKB:IKB)  =                                                 &
                 - MYF (                                          &
ZFLXZ(:,:,IKB+KKL:IKB+KKL) * (PVM(:,:,IKB+KKL:IKB+KKL)-PVM(:,:,IKB:IKB))  &
                       / MYM(PDZZ(:,:,IKB+KKL:IKB+KKL))               &
                       )
!
IF (LOCEAN) THEN
  ! evaluate the dynamic production at w(IKE-KKL) in PDP(IKE)
  ZA(:,:,IKE:IKE) = - MYF (                                                  &
   ZFLXZ(:,:,IKE-KKL:IKE-KKL) * (PVM(:,:,IKE:IKE)-PVM(:,:,IKE-KKL:IKE-KKL))  &
                          / MYM(PDZZ(:,:,IKE-KKL:IKE-KKL))                   &
                          )
END IF
!
PDP(:,:,:)=PDP(:,:,:)+ZA(:,:,:)
!
! Storage in the LES configuration
!
IF (LLES_CALL) THEN
  CALL SECOND_MNH(ZTIME1)
  CALL LES_MEAN_SUBGRID( MZF(MYF(ZFLXZ)), X_LES_SUBGRID_WV )
  CALL LES_MEAN_SUBGRID( MZF(MYF(GZ_V_VW(PVM,PDZZ)*&
                    & ZFLXZ)), X_LES_RES_ddxa_V_SBG_UaV )
  CALL SECOND_MNH(ZTIME2)
  XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
END IF
!
!
!*       6.3  Source of W wind component 
!
IF(HTURBDIM=='3DIM') THEN
  ! Compute the source for the W wind component
  ZFLXZ(:,:,KKA) = 2 * ZFLXZ(:,:,IKB) - ZFLXZ(:,:,IKB+KKL) ! extrapolation 
  IF (LOCEAN) THEN
    ZFLXZ(:,:,KKU) = 2 * ZFLXZ(:,:,IKE) - ZFLXZ(:,:,IKE-KKL) ! extrapolation 
  END IF
  !
  IF (.NOT. L2D) THEN 
    IF (.NOT. LFLAT) THEN
      PRWS(:,:,:)= PRWS(:,:,:)                               &
                  -DYF( MZM( MYM(PRHODJ) /PDYY ) * ZFLXZ )   &
                  +DZM( PRHODJ / MZF(PDZZ ) *                &
                        MYF( MZF( ZFLXZ*PDZY ) / PDYY )      &
                      )
    ELSE
      PRWS(:,:,:)= PRWS(:,:,:) -DYF( MZM( MYM(PRHODJ) /PDYY ) * ZFLXZ )
    END IF
  END IF
  ! 
  ! Complete the Dynamical production with the W wind component 
  IF (.NOT. L2D) THEN
    ZA(:,:,:) = - MZF( MYF ( ZFLXZ * GY_W_VW( PWM,PDYY,PDZZ,PDZY) )  )
  !
  ! evaluate the dynamic production at w(IKB+KKL) in PDP(IKB)
    ZA(:,:,IKB:IKB) = - MYF (                                              &
     ZFLXZ(:,:,IKB+KKL:IKB+KKL) *                                          &
       ( DYM( PWM(:,:,IKB+KKL:IKB+KKL) )                                   &
        -MYM(  (PWM(:,:,IKB+2*KKL:IKB+2*KKL)-PWM(:,:,IKB+KKL:IKB+KKL))     &
                /(PDZZ(:,:,IKB+2*KKL:IKB+2*KKL)+PDZZ(:,:,IKB+KKL:IKB+KKL)) &
              +(PWM(:,:,IKB+KKL:IKB+KKL)-PWM(:,:,IKB:IKB  ))               &
                /(PDZZ(:,:,IKB+KKL:IKB+KKL)+PDZZ(:,:,IKB:IKB  ))           &
            )                                                              &
          * PDZY(:,:,IKB+KKL:IKB+KKL)                                      &
       ) / (0.5*(PDYY(:,:,IKB+KKL:IKB+KKL)+PDYY(:,:,IKB:IKB)))             &
                            )
  !
    IF (LOCEAN) THEN
     ZA(:,:,IKE:IKE) = - MYF (                                              &
      ZFLXZ(:,:,IKE-KKL:IKE-KKL) *                                          &
        ( DYM( PWM(:,:,IKE-KKL:IKE-KKL) )                                   &
         -MYM(  (PWM(:,:,IKE-2*KKL:IKE-2*KKL)-PWM(:,:,IKE-KKL:IKE-KKL))     &
                 /(PDZZ(:,:,IKE-2*KKL:IKE-2*KKL)+PDZZ(:,:,IKE-KKL:IKE-KKL)) &
               +(PWM(:,:,IKE-KKL:IKE-KKL)-PWM(:,:,IKE:IKE  ))               &
                 /(PDZZ(:,:,IKE-KKL:IKE-KKL)+PDZZ(:,:,IKE:IKE  ))           &
             )                                                              &
           * PDZY(:,:,IKE-KKL:IKE-KKL)                                      &
        ) / (0.5*(PDYY(:,:,IKE-KKL:IKE-KKL)+PDYY(:,:,IKE:IKE)))             &
                            )
    END IF
!    
    PDP(:,:,:)=PDP(:,:,:)+ZA(:,:,:)
  !
  END IF
  !
  ! Storage in the LES configuration
  !
  IF (LLES_CALL) THEN
    CALL SECOND_MNH(ZTIME1)
    CALL LES_MEAN_SUBGRID( MZF(MYF(GY_W_VW(PWM,PDYY,&
    PDZZ,PDZY)*ZFLXZ)), X_LES_RES_ddxa_W_SBG_UaW , .TRUE. )
    CALL LES_MEAN_SUBGRID( MYF(GY_M_V(KKA,KKU,KKL,PTHLM,PDYY,PDZZ,PDZY)&
    *MZF(ZFLXZ)), X_LES_RES_ddxa_Thl_SBG_UaW , .TRUE. )
    IF (KRR>=1) THEN
      CALL LES_MEAN_SUBGRID( MYF(GY_V_M(PRM(:,:,:,1),PDYY,PDZZ,&
      PDZY)*MZF(ZFLXZ)),X_LES_RES_ddxa_Rt_SBG_UaW , .TRUE. )
    END IF
    CALL SECOND_MNH(ZTIME2)
    XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
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
IF ( OTURB_FLX .AND. tpfile%lopened .AND. HTURBDIM == '1DIM') THEN
  ZFLXZ(:,:,:)= (2./3.) * PTKEM(:,:,:)                     &
     -XCMFS*PLM(:,:,:)*SQRT(PTKEM(:,:,:))*GZ_W_M(PWM,PDZZ)
  ! to be tested &
  !   +XCMFB*(4./3.)*PLM(:,:,:)/SQRT(PTKEM(:,:,:))*PTP(:,:,:) 
  ! stores the W variance
  TZFIELD%CMNHNAME   = 'W_VVAR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'W_VVAR'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_W_VVAR'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZ)
END IF
!
!----------------------------------------------------------------------------
!
END SUBROUTINE TURB_VER_DYN_FLUX
