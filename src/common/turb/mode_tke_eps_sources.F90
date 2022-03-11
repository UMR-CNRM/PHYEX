!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TKE_EPS_SOURCES
IMPLICIT NONE
CONTAINS
      SUBROUTINE TKE_EPS_SOURCES(CST,CSTURB,BUCONF,HPROGRAM, &
                    & KKA,KKU,KKL,KMI,PTKEM,PLM,PLEPS,PDP,  &
                    & PTRH,PRHODJ,PDZZ,PDXX,PDYY,PDZX,PDZY,PZZ,        &
                    & PTSTEP,PIMPL,PEXPL,                              &
                    & HTURBLEN,HTURBDIM,                               &
                    & TPFILE,OTURB_DIAG,OLES_CALL,                     &
                    & PTP,PRTKES,PRTHLS,PCOEF_DISS,PTDIFF,PTDISS,PRTKEMS,&
                    & TBUDGETS, KBUDGETS, &
                    & PEDR, PTR,PDISS                                  )
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
!!           CSTURB%XCET,CSTURB%XCED  : transport and dissipation cts. for the TKE
!!           XCDP,XCDD,XCDT: constants from the parameterization of
!!                        the K-epsilon equation
!!           CSTURB%XTKEMIN,XEPSMIN : minimum values for the TKE and its
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
!!                     2015-01 (J. Escobar) missing get_halo(ZRES) for JPHEXT<> 1 
!!     J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
! --------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
USE MODD_ARGSLIST_ll,    ONLY: LIST_ll
USE MODD_BUDGET, ONLY: TBUDGETCONF_t, NBUDGET_TKE, NBUDGET_TH, TBUDGETDATA
USE MODD_CST, ONLY: CST_t
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_DIAG_IN_RUN, ONLY : LDIAG_IN_RUN, XCURRENT_TKE_DISS
USE MODD_FIELD, ONLY: TFIELDDATA, TYPEREAL
USE MODD_IO, ONLY: TFILEDATA
USE MODD_LES
USE MODD_PARAMETERS
!
USE MODE_BUDGET, ONLY: BUDGET_STORE_ADD, BUDGET_STORE_END, BUDGET_STORE_INIT
USE MODE_IO_FIELD_WRITE, ONLY: IO_FIELD_WRITE
USE MODE_ll
!
USE MODI_GET_HALO
USE MODI_GRADIENT_M
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_LES_MEAN_SUBGRID
USE MODE_TRIDIAG_TKE, ONLY: TRIDIAG_TKE
USE MODI_SHUMAN , ONLY : DZM, DZF, MZM, MZF
!
!
IMPLICIT NONE
!
!
!*       0.1  declarations of arguments
!
!
TYPE(CST_t),             INTENT(IN)    :: CST
TYPE(CSTURB_t),          INTENT(IN)    :: CSTURB
TYPE(TBUDGETCONF_t),     INTENT(IN)    :: BUCONF
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
REAL,                    INTENT(IN)   ::  PTSTEP       ! Time step 
REAL,                    INTENT(IN)   ::  PEXPL, PIMPL ! Coef. temporal. disc.
CHARACTER(LEN=4),        INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                       ! turbulence scheme
CHARACTER(LEN=4),        INTENT(IN)   ::  HTURBLEN     ! kind of mixing length
CHARACTER(LEN=6),        INTENT(IN)   ::  HPROGRAM     ! CPROGRAM is the program currently running (modd_conf)
TYPE(TFILEDATA),         INTENT(IN)   ::  TPFILE       ! Output file
LOGICAL,                 INTENT(IN)   ::  OLES_CALL    !
LOGICAL,                 INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                  ! diagnostic fields in the syncronous FM-file
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PDP          ! Dyn. prod. of TKE
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PTRH
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PTP          ! Ther. prod. of TKE
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PRTKES       ! RHOD * Jacobian *
                                                       ! TKE at t+deltat
REAL, DIMENSION(:,:,:),  INTENT(INOUT)::  PRTHLS       ! Source of Theta_l
REAL, DIMENSION(:,:,:),  INTENT(IN)   ::  PCOEF_DISS   ! 1/(Cph*Exner)
REAL, DIMENSION(:,:,:),  INTENT(OUT)  ::  PTDIFF       ! Diffusion TKE term
REAL, DIMENSION(:,:,:),  INTENT(OUT)  ::  PTDISS       ! Dissipation TKE term
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::  PRTKEMS      ! Advection source
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
REAL, DIMENSION(:,:,:),  INTENT(OUT), OPTIONAL  ::  PTR          ! Transport prod. of TKE
REAL, DIMENSION(:,:,:),  INTENT(OUT), OPTIONAL  ::  PDISS        ! Dissipation of TKE
REAL, DIMENSION(:,:,:),  INTENT(OUT), OPTIONAL  ::  PEDR         ! EDR 
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
       ZKEFF,    & ! effectif diffusion coeff = LT * SQRT( TKE )
       ZTR         ! Transport term
LOGICAL,DIMENSION(SIZE(PTKEM,1),SIZE(PTKEM,2),SIZE(PTKEM,3)) :: GTKENEG
                   ! 3D mask .T. if TKE < CSTURB%XTKEMIN
INTEGER             :: IIB,IIE,IJB,IJE,IKB,IKE
                                    ! Index values for the Beginning and End
                                    ! mass points of the domain 
INTEGER             :: IIU,IJU,IKU  ! array size in the 3 dimensions 
!
TYPE(LIST_ll), POINTER :: TZFIELDDISS_ll ! list of fields to exchange
INTEGER                :: IINFO_ll       ! return code of parallel routine
TYPE(TFIELDDATA) :: TZFIELD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
NULLIFY(TZFIELDDISS_ll)
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!

IF (LHOOK) CALL DR_HOOK('TKE_EPS_SOURCES',0,ZHOOK_HANDLE)
!
IIU=SIZE(PTKEM,1)
IJU=SIZE(PTKEM,2)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE,IIU,IJU)
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
! extrapolate the dynamic production with a 1/Z law from its value at the 
! W(IKB+1) value stored in PDP(IKB) to the mass localization tke(IKB)
PDP(:,:,IKB) = PDP(:,:,IKB) * (1. + PDZZ(:,:,IKB+KKL)/PDZZ(:,:,IKB))
!
! Compute the source terms for TKE: ( ADVECtion + NUMerical DIFFusion + ..)
! + (Dynamical Production) + (Thermal Production) - (dissipation) 
ZFLX(:,:,:) = CSTURB%XCED * SQRT(PTKEM(:,:,:)) / PLEPS(:,:,:)
ZSOURCE(:,:,:) = ( PRTKES(:,:,:) +  PRTKEMS(:,:,:) ) / PRHODJ(:,:,:) &
   - PTKEM(:,:,:) / PTSTEP &
   + PDP(:,:,:) + PTP(:,:,:) + ZTR(:,:,:) - PEXPL * ZFLX(:,:,:) * PTKEM(:,:,:)
!
!*       2.2  implicit vertical TKE transport
!
!
! Compute the vector giving the elements just under the diagonal for the 
! matrix inverted in TRIDIAG 
!
ZA(:,:,:)     = - PTSTEP * CSTURB%XCET * &
                MZM(ZKEFF, KKA, KKU, KKL) * MZM(PRHODJ, KKA, KKU, KKL) / PDZZ**2
!
! Compute TKE at time t+deltat: ( stored in ZRES )
!
CALL TRIDIAG_TKE(KKA,KKU,KKL,PTKEM,ZA,PTSTEP,PEXPL,PIMPL,PRHODJ,&
            & ZSOURCE,PTSTEP*ZFLX,ZRES)
CALL GET_HALO(ZRES)
!
!* diagnose the dissipation
!
IF (LDIAG_IN_RUN) THEN
  XCURRENT_TKE_DISS = ZFLX(:,:,:) * PTKEM(:,:,:) &
                                  *(PEXPL*PTKEM(:,:,:) + PIMPL*ZRES(:,:,:))
  CALL ADD3DFIELD_ll( TZFIELDDISS_ll, XCURRENT_TKE_DISS, 'TKE_EPS_SOURCES::XCURRENT_TKE_DISS' )
  CALL UPDATE_HALO_ll(TZFIELDDISS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDDISS_ll)
ENDIF
!
! TKE must be greater than its minimum value
! CL : Now done at the end of the time step in ADVECTION_METSV for MesoNH
IF(HPROGRAM/='MESONH') THEN
 GTKENEG =  ZRES <= CSTURB%XTKEMIN 
 WHERE ( GTKENEG ) 
   ZRES = CSTURB%XTKEMIN
 END WHERE
END IF
PTDISS(:,:,:) = - ZFLX(:,:,:)*(PEXPL*PTKEM(:,:,:) + PIMPL*ZRES(:,:,:))
!
IF ( OLES_CALL .OR.                         &
     (OTURB_DIAG .AND. TPFILE%LOPENED)  ) THEN
!
! Compute the cartesian vertical flux of TKE in ZFLX
!

  ZFLX(:,:,:)   = - CSTURB%XCET * MZM(ZKEFF, KKA, KKU, KKL) *   &
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
  IF (OLES_CALL) THEN
    CALL LES_MEAN_SUBGRID(MZF(ZFLX, KKA, KKU, KKL), X_LES_SUBGRID_WTke )
    CALL LES_MEAN_SUBGRID(-ZTR, X_LES_SUBGRID_ddz_WTke )
  END IF
!
END IF
!
!*       2.4  stores the explicit sources for budget purposes
!
IF (BUCONF%LBUDGET_TKE) THEN
  ! Dynamical production
  CALL BUDGET_STORE_ADD( TBUDGETS(NBUDGET_TKE), 'DP', PDP(:, :, :) * PRHODJ(:, :, :) )
  ! Thermal production
  CALL BUDGET_STORE_ADD( TBUDGETS(NBUDGET_TKE), 'TP', PTP(:, :, :) * PRHODJ(:, :, :) )
  ! Dissipation
  CALL BUDGET_STORE_ADD( TBUDGETS(NBUDGET_TKE), 'DISS',- CSTURB%XCED * SQRT(PTKEM(:,:,:)) / PLEPS(:,:,:) * &
                (PEXPL*PTKEM(:,:,:) + PIMPL*ZRES(:,:,:)) * PRHODJ(:,:,:))
END IF 
!
!*       2.5  computes the final RTKE and stores the whole turbulent transport
!              with the removal of the advection part for MesoNH

!Store the previous source terms in prtkes before initializing the next one
!Should be in IF LBUDGET_TKE only. Was removed out for a correct comput. of PTDIFF in case of LBUDGET_TKE=F in AROME
PRTKES(:,:,:) = PRTKES(:,:,:) + PRHODJ(:,:,:) *                                                           &
                ( PDP(:,:,:) + PTP(:,:,:)                                                                 &
                  - CSTURB%XCED * SQRT(PTKEM(:,:,:)) / PLEPS(:,:,:) * ( PEXPL*PTKEM(:,:,:) + PIMPL*ZRES(:,:,:) ) )
!
PTDIFF(:,:,:) =  ZRES(:,:,:) / PTSTEP - PRTKES(:,:,:)/PRHODJ(:,:,:) &
 & - PDP(:,:,:)- PTP(:,:,:) - PTDISS(:,:,:)

IF (BUCONF%LBUDGET_TKE) CALL BUDGET_STORE_INIT( TBUDGETS(NBUDGET_TKE), 'TR', PRTKES(:, :, :) )
!
PRTKES(:,:,:) = ZRES(:,:,:) * PRHODJ(:,:,:) / PTSTEP -  PRTKEMS(:,:,:)
!
! stores the whole turbulent transport
!
IF (BUCONF%LBUDGET_TKE) CALL BUDGET_STORE_END( TBUDGETS(NBUDGET_TKE), 'TR', PRTKES(:, :, :) )

!----------------------------------------------------------------------------
!
!*       3.   COMPUTE THE DISSIPATIVE HEATING
!             -------------------------------
!
PRTHLS(:,:,:) = PRTHLS(:,:,:) + CSTURB%XCED * SQRT(PTKEM(:,:,:)) / PLEPS(:,:,:) * &
                (PEXPL*PTKEM(:,:,:) + PIMPL*ZRES(:,:,:)) * PRHODJ(:,:,:) * PCOEF_DISS(:,:,:)
!----------------------------------------------------------------------------
!
!*       4.   STORES SOME DIAGNOSTICS
!             -----------------------
!
IF(PRESENT(PDISS)) PDISS(:,:,:) =  -CSTURB%XCED * (PTKEM(:,:,:)**1.5) / PLEPS(:,:,:)
IF(PRESENT(PTR)) PTR=ZTR
IF(PRESENT(PEDR)) PEDR = CSTURB%XCED * (PTKEM(:,:,:)**1.5) / PLEPS(:,:,:)
!
IF ( OTURB_DIAG .AND. TPFILE%LOPENED ) THEN
!
! stores the dynamic production 
!
  TZFIELD%CMNHNAME   = 'DP'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'DP'
  TZFIELD%CUNITS     = 'm2 s-3'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_DP'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,PDP)
!
! stores the thermal production 
!
  TZFIELD%CMNHNAME   = 'TP'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'TP'
  TZFIELD%CUNITS     = 'm2 s-3'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_TP'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,PTP)
!
! stores the whole turbulent transport
!
  TZFIELD%CMNHNAME   = 'TR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'TR'
  TZFIELD%CUNITS     = 'm2 s-3'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_TR'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,ZTR)
!
! stores the dissipation of TKE 
!
  TZFIELD%CMNHNAME   = 'DISS'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'DISS'
  TZFIELD%CUNITS     = 'm2 s-3'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_DISS'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_FIELD_WRITE(TPFILE,TZFIELD,PDISS)
END IF
!
! Storage in the LES configuration of the Dynamic Production of TKE and
! the dissipation of TKE 
! 
IF (OLES_CALL ) THEN
  CALL LES_MEAN_SUBGRID( PDISS, X_LES_SUBGRID_DISS_Tke )
END IF
!
!----------------------------------------------------------------------------
! 
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TKE_EPS_SOURCES',1,ZHOOK_HANDLE)
END SUBROUTINE TKE_EPS_SOURCES
END MODULE MODE_TKE_EPS_SOURCES
