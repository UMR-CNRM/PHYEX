!     ######spl
      SUBROUTINE TKE_EPS_SOURCES(KKA,KKU,KKL,KMI,PTKEM,PLM,PLEPS,PDP,  &
                    & PTRH,PRHODJ,PDZZ,PDXX,PDYY,PDZX,PDZY,PZZ,        &
                    & PTSTEP,PIMPL,PEXPL,                              &
                    & HTURBLEN,HTURBDIM,                               &
                    & HFMFILE,HLUOUT,OCLOSE_OUT,OTURB_DIAG,            &
                 & PTP,PRTKES,PRTHLS,PCOEF_DISS,PTDIFF,PTDISS, &
                 & PEDR,YDDDH, YDLDDH, YDMDDH)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##################################################################
!
!
!!****  *TKE_EPS_SOURCES* - routine to compute the sources of the turbulent 
!!      evolutive variables: TKE and its dissipation when it is taken into 
!!      account. The contribution to the heating of tke dissipation is computed.
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute the sources necessary for
!     the evolution of the turbulent kinetic energy and its dissipation 
!     if necessary.
!
!!**  METHOD
!!    ------
!!      The vertical turbulent flux is computed in an off-centered 
!!    implicit scheme (a Crank-Nicholson type with coefficients different 
!!    than 0.5), which allows to vary the degree of implicitness of the 
!!    formulation.
!!      In high resolution, the horizontal transport terms are also
!!    calculated, but explicitly. 
!!      The evolution of the dissipation as a variable is made if 
!!    the parameter HTURBLEN is set equal to KEPS. The same reasoning 
!!    made for TKE applies.
!!
!!    EXTERNAL
!!    --------
!!      GX_U_M,GY_V_M,GZ_W_M
!!      GX_M_U,GY_M_V          :  Cartesian vertical gradient operators
!!
!!      MXF,MXM.MYF,MYM,MZF,MZM:  Shuman functions (mean operators)
!!      DZF                    :  Shuman functions (difference operators)     
!!
!!      SUBROUTINE TRIDIAG     :  to solve an implicit temporal scheme
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
!!           XCET,XCED  : transport and dissipation cts. for the TKE
!!           XCDP,XCDD,XCDT: constants from the parameterization of
!!                        the K-epsilon equation
!!           XTKEMIN,XEPSMIN : minimum values for the TKE and its
!!                        dissipation
!!
!!      Module MODD_PARAMETERS: 
!!
!!           JPVEXT_TURB
!!      Module MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask 
!!                          'NONE'  ' for no budget
!!         NBUPROCCTR   : process counter used for each budget variable
!!         LBU_RTKE     : logical for budget of RTKE (turbulent kinetic energy)
!!                        .TRUE. = budget of RTKE       
!!                        .FALSE. = no budget of RTKE
!!
!!
!!    REFERENCE
!!    ---------
!!      Book 2 of documentation (routine TKE_EPS_SOURCES)
!!      Book 1 of documentation (Chapter: Turbulence)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart             * INM and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       August 23, 1994
!!      Modifications: Feb 14, 1995 (J.Cuxart and J.Stein)
!!                                  Doctorization and Optimization
!!                     June 29, 1995 (J.Stein) TKE budget
!!                     June 28, 1995 (J.Cuxart) Add LES tools
!!      Modifications: February 29, 1996 (J. Stein) optimization
!!      Modifications: May 6, 1996 (N. Wood) Extend some loops over
!!                                              the outer points
!!      Modifications: August 30, 1996 (P. Jabouille)  calcul ZFLX at the
!!                                                      IKU level
!!                     October 10, 1996 (J.Stein)  set Keff at t-deltat
!!                     Oct 8, 1996 (Cuxart,Sanchez) Var.LES: XETR_TF,XDISS_TF
!!                     December 20, 1996 (J.-P. Pinty) update the CALL BUDGET
!!                     November 24, 1997 (V. Masson) bug in <v'e>
!!                                                   removes the DO loops
!!                     Augu. 9, 1999 (J.Stein) TKE budget correction
!!                     Mar 07  2001 (V. Masson and J. Stein) remove the horizontal 
!!                                         turbulent transports of Tke computation
!!                     Nov 06, 2002 (V. Masson) LES budgets
!!                     July 20, 2003 (J.-P. Pinty P Jabouille) add the dissipative heating
!!                     May   2006    Remove KEPS
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!                     2012-02 Y. Seity,  add possibility to run with reversed 
!!                                    vertical levels
!!                     2014-11 Y. Seity,  add output terms for TKE DDHs budgets
!! --------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_CST
USE MODD_CONF
USE MODD_CTURB
USE MODD_PARAMETERS
USE MODD_BUDGET
USE MODD_LES
USE MODD_DIAG_IN_RUN, ONLY : LDIAG_IN_RUN, XCURRENT_TKE_DISS
!
USE MODI_GRADIENT_M
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_SHUMAN , ONLY : DZM, DZF, MZM, MZF
USE MODI_TRIDIAG 
USE MODI_TRIDIAG_TKE
USE MODI_BUDGET
USE MODE_FMWRIT
USE MODI_LES_MEAN_SUBGRID
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
!!!!!AROME!!USE MODE_ll
!!!!!AROME!!USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
IMPLICIT NONE
!
!
!*       0.1  declarations of arguments
!
!
INTEGER,                 INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                 INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                 INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=ARO

INTEGER,                 INTENT(IN)   ::  KMI          ! model index number  
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PTKEM        ! TKE at t-deltat
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PLM          ! mixing length         
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PLEPS        ! dissipative length
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PDXX,PDYY,PDZZ,PDZX,PDZY
                                                       ! metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PZZ          ! physical height w-pt
REAL,                    INTENT(IN)   ::  PTSTEP       ! Double Time step ( *.5 for
                                                       ! the first time step ) 
REAL,                    INTENT(IN)   ::  PEXPL, PIMPL ! Coef. temporal. disc.
CHARACTER*4,             INTENT(IN)   ::  HTURBDIM     ! dimensionality of the 
                                                       ! turbulence scheme
CHARACTER*4,             INTENT(IN)   ::  HTURBLEN     ! kind of mixing length 
CHARACTER(LEN=*),        INTENT(IN)   ::  HFMFILE      ! Name of the output
                                                       ! FM-file
CHARACTER(LEN=*),        INTENT(IN)   ::  HLUOUT       ! Output-listing name for
                                                       ! model n
LOGICAL,                 INTENT(IN)   ::  OCLOSE_OUT   ! switch for syncronous
                                                       ! file opening
LOGICAL,                 INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                  ! diagnostic fields in the syncronous FM-file
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PDP, PTRH          ! Dyn. prod. of TKE
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PTP          ! Ther. prod. of TKE
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PRTKES       ! RHOD * Jacobian *
                                                       ! TKE at t+deltat
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PRTHLS       ! Source of Theta_l
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PCOEF_DISS   ! 1/(Cph*Exner)
REAL, DIMENSION(:,:,:),  INTENT(OUT)  ::  PTDIFF     ! Diffusion TKE term
REAL, DIMENSION(:,:,:),  INTENT(OUT)  ::  PTDISS     ! Dissipation TKE term
REAL, DIMENSION(:,:,:),  INTENT(OUT)   ::  PEDR         ! EDR 
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PTKEM,1),SIZE(PTKEM,2),SIZE(PTKEM,3))::         &
       ZA,       & ! under diagonal elements of the tri-diagonal matrix involved
                   ! in the temporal implicit scheme
       ZRES,     & ! treated variable at t+ deltat when the turbu-
                   ! lence is the only source of evolution added to the ones
                   ! considered in ZSOURCE. This variable is also used to
                   ! temporarily store some diagnostics stored in FM file
       ZFLX,     & ! horizontal or vertical flux of the treated variable
       ZSOURCE,  & ! source of evolution for the treated variable
       ZTR,      & ! turbulent transport of TKE 
       ZKEFF       ! effectif diffusion coeff = LT * SQRT( TKE )
LOGICAL,DIMENSION(SIZE(PTKEM,1),SIZE(PTKEM,2),SIZE(PTKEM,3)) :: GTKENEG
                   ! 3D mask .T. if TKE < XTKEMIN
INTEGER             :: IIB,IIE,IJB,IJE,IKB,IKE
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain 
INTEGER             :: IIU,IJU,IKU  ! array size in the 3 dimensions 
INTEGER             :: IRESP        ! Return code of FM routines
INTEGER             :: IGRID        ! C-grid indicator in LFIFM file
INTEGER             :: ILENCH       ! Length of comment string in LFIFM file
CHARACTER (LEN=100) :: YCOMMENT     ! comment string in LFIFM file
CHARACTER (LEN=16)  :: YRECFM       ! Name of the desired field in LFIFM file
!
!!!!!AROME!!TYPE(LIST_ll), POINTER :: TZFIELDDISS_ll ! list of fields to exchange
!!!!!AROME!!INTEGER                :: IINFO_ll       ! return code of parallel routine
!

!----------------------------------------------------------------------------
!!!!!AROME!!NULLIFY(TZFIELDDISS_ll)
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TKE_EPS_SOURCES',0,ZHOOK_HANDLE)
IIB=1+JPHEXT
IIU=SIZE(PTKEM,1)
IIE=IIU-JPHEXT
IJB=1+JPHEXT
IJU=SIZE(PTKEM,2)
IJE=IJU-JPHEXT
IKB=KKA+JPVEXT_TURB*KKL
IKE=KKU-JPVEXT_TURB*KKL
!
! compute the effective diffusion coefficient at the mass point
ZKEFF(:,:,:) = PLM(:,:,:) * SQRT(PTKEM(:,:,:)) 
!
!----------------------------------------------------------------------------
!
!*       2.   TKE EQUATION  
!             ------------
!
!*       2.1  Horizontal turbulent explicit transport
!
!
! Complete the sources of TKE with the horizontal turbulent explicit transport
!
IF (HTURBDIM=='3DIM') THEN
  ZTR=PTRH
ELSE
  ZTR=0.
END IF
!
!
!
!*       2.2  Explicit TKE sources except horizontal turbulent transport 
!
!
! extrapolate the dynamic production with a 1/Z law from its value at the 
! W(IKB+1) value stored in PDP(IKB) to the mass localization tke(IKB)
PDP(:,:,IKB) = PDP(:,:,IKB) * (1. + PDZZ(:,:,IKB+KKL)/PDZZ(:,:,IKB))
!
! Compute the source terms for TKE: ( ADVECtion + NUMerical DIFFusion + ..)
! + (Dynamical Production) + (Thermal Production) - (dissipation) 
ZFLX(:,:,:) = XCED * SQRT(PTKEM(:,:,:)) / PLEPS(:,:,:)
ZSOURCE(:,:,:) = PRTKES(:,:,:) / PRHODJ(:,:,:) - PTKEM(:,:,:) / PTSTEP &
   + PDP(:,:,:) + PTP(:,:,:) + ZTR(:,:,:) - PEXPL * ZFLX(:,:,:) * PTKEM(:,:,:)
!
!*       2.2  implicit vertical TKE transport
!
!
! Compute the vector giving the elements just under the diagonal for the 
! matrix inverted in TRIDIAG 
!
ZA(:,:,:)     = - PTSTEP * XCET * &
                MZM(ZKEFF, KKA, KKU, KKL) * MZM(PRHODJ, KKA, KKU, KKL) / PDZZ**2
!
! Compute TKE at time t+deltat: ( stored in ZRES )
!
CALL TRIDIAG_TKE(KKA,KKU,KKL,PTKEM,ZA,PTSTEP,PEXPL,PIMPL,PRHODJ,&
            & ZSOURCE,PTSTEP*ZFLX,ZRES)
!
!* diagnose the dissipation
!
IF (LDIAG_IN_RUN) THEN
  XCURRENT_TKE_DISS = ZFLX(:,:,:) * PTKEM(:,:,:) &
                                  *(PEXPL*PTKEM(:,:,:) + PIMPL*ZRES(:,:,:))
!!!!!AROME!!  CALL ADD3DFIELD_ll(TZFIELDDISS_ll,XCURRENT_TKE_DISS)
!!!!!AROME!!  CALL UPDATE_HALO_ll(TZFIELDDISS_ll,IINFO_ll)
!!!!!AROME!!  CALL CLEANLIST_ll(TZFIELDDISS_ll)
ENDIF
!
! TKE must be greater than its minimum value
!
GTKENEG =  ZRES <= XTKEMIN 
WHERE ( GTKENEG ) 
  ZRES = XTKEMIN
END WHERE

PTDISS(:,:,:) = - ZFLX(:,:,:)*(PEXPL*PTKEM(:,:,:) + PIMPL*ZRES(:,:,:))
!
IF ( LLES_CALL .OR.                         &
     (OTURB_DIAG .AND. OCLOSE_OUT)  ) THEN
!
! Compute the cartesian vertical flux of TKE in ZFLX
!

  ZFLX(:,:,:)   = - XCET * MZM(ZKEFF, KKA, KKU, KKL) *   &
                  DZM(PIMPL * ZRES + PEXPL * PTKEM, KKA, KKU, KKL) / PDZZ
!
  ZFLX(:,:,IKB) = 0.
  ZFLX(:,:,KKA) = 0.
!
! Compute the whole turbulent TRansport of TKE:
!
  ZTR(:,:,:)= ZTR - DZF(MZM(PRHODJ, KKA, KKU, KKL) * ZFLX / PDZZ, KKA, KKU, KKL) /PRHODJ
!
! Storage in the LES configuration
!
  IF (LLES_CALL) THEN
    CALL LES_MEAN_SUBGRID(MZF(ZFLX, KKA, KKU, KKL), X_LES_SUBGRID_WTke )
    CALL LES_MEAN_SUBGRID(-ZTR, X_LES_SUBGRID_ddz_WTke )
  END IF
!
END IF
!
!*       2.4  stores the explicit sources for budget purposes
!
IF (LBUDGET_TKE) THEN
!
! add the dynamical production
!
  PRTKES(:,:,:) = PRTKES(:,:,:) + PDP(:,:,:) * PRHODJ(:,:,:)
  CALL BUDGET (PRTKES(:,:,:),5,'DP_BU_RTKE',YDDDH, YDLDDH, YDMDDH)
!
! add the thermal production
!
  PRTKES(:,:,:) = PRTKES(:,:,:) + PTP(:,:,:) * PRHODJ(:,:,:)
  CALL BUDGET (PRTKES(:,:,:),5,'TP_BU_RTKE',YDDDH, YDLDDH, YDMDDH)
!
! add the dissipation
!
PRTKES(:,:,:) = PRTKES(:,:,:) - XCED * SQRT(PTKEM(:,:,:)) / PLEPS(:,:,:) * &
                (PEXPL*PTKEM(:,:,:) + PIMPL*ZRES(:,:,:)) * PRHODJ(:,:,:)
CALL BUDGET (PRTKES(:,:,:),5,'DISS_BU_RTKE',YDDDH, YDLDDH, YDMDDH)
END IF 
!
!*       2.5  computes the final RTKE and stores the whole turbulent transport
!
PTDIFF(:,:,:) =  ZRES(:,:,:) / PTSTEP - PRTKES(:,:,:)/PRHODJ(:,:,:) &
 & - PDP(:,:,:)- PTP(:,:,:) - PTDISS(:,:,:)

PRTKES(:,:,:) = ZRES(:,:,:) * PRHODJ(:,:,:) / PTSTEP

!
! stores the whole turbulent transport
!
IF (LBUDGET_TKE) CALL BUDGET (PRTKES(:,:,:),5,'TR_BU_RTKE',YDDDH, YDLDDH, YDMDDH)
!
!
!----------------------------------------------------------------------------
!
!*       3.   COMPUTE THE DISSIPATIVE HEATING
!             -------------------------------
!
PRTHLS(:,:,:) = PRTHLS(:,:,:) + XCED * SQRT(PTKEM(:,:,:)) / PLEPS(:,:,:) * &
                (PEXPL*PTKEM(:,:,:) + PIMPL*ZRES(:,:,:)) * PRHODJ(:,:,:) * PCOEF_DISS(:,:,:)
!
!----------------------------------------------------------------------------
!
!*       4.   STORES SOME DIAGNOSTICS
!             -----------------------
!
PEDR(:,:,:)=XCED * (PTKEM(:,:,:)**1.5) / PLEPS(:,:,:)



IF ( OTURB_DIAG .AND. OCLOSE_OUT ) THEN
!
! stores the dynamic production 
!
  YRECFM  ='DP'
  YCOMMENT='X_Y_Z_DP (M**2/S**3)'
  IGRID   = 1
  ILENCH=LEN(YCOMMENT) 
  CALL FMWRIT(HFMFILE,YRECFM,HLUOUT,'XY',PDP,IGRID,ILENCH,YCOMMENT,IRESP)
!
! stores the thermal production 
!
  YRECFM  ='TP'
  YCOMMENT='X_Y_Z_TP (M**2/S**3)'
  IGRID   = 1
  ILENCH=LEN(YCOMMENT)
  CALL FMWRIT(HFMFILE,YRECFM,HLUOUT,'XY',PTP,IGRID,ILENCH,YCOMMENT,IRESP)
!
! stores the whole turbulent transport
!
  YRECFM  ='TR'
  YCOMMENT='X_Y_Z_TR (M**2/S**3)'
  IGRID   = 1
  ILENCH=LEN(YCOMMENT)
  CALL FMWRIT(HFMFILE,YRECFM,HLUOUT,'XY',ZTR,IGRID,ILENCH,YCOMMENT,IRESP)
!
! stores the dissipation of TKE 
!
  YRECFM  ='DISS'
  YCOMMENT='X_Y_Z_DISS (M**2/S**3)'
  IGRID   = 1
  ILENCH=LEN(YCOMMENT)
  ZFLX(:,:,:) =-XCED * (PTKEM(:,:,:)**1.5) / PLEPS(:,:,:) 
  CALL FMWRIT(HFMFILE,YRECFM,HLUOUT,'XY',ZFLX,IGRID,ILENCH,YCOMMENT,IRESP)
END IF
!
! Storage in the LES configuration of the Dynamic Production of TKE and
! the dissipation of TKE 
! 
IF (LLES_CALL ) THEN
  ZFLX(:,:,:) =-XCED * (PTKEM(:,:,:)**1.5) / PLEPS(:,:,:) 
  CALL LES_MEAN_SUBGRID( ZFLX, X_LES_SUBGRID_DISS_Tke )
END IF
! 
!----------------------------------------------------------------------------
! 
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TKE_EPS_SOURCES',1,ZHOOK_HANDLE)
END SUBROUTINE TKE_EPS_SOURCES
