!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_TKE_EPS_SOURCES
IMPLICIT NONE
CONTAINS
      SUBROUTINE TKE_EPS_SOURCES(D,CST,CSTURB,BUCONF,TURBN,TLES,       &
                    & HPROGRAM, KMI,PTKEM,PLM,PLEPS,PDP,               &
                    & PTRH,PRHODJ,PDZZ,PDXX,PDYY,PDZX,PDZY,PZZ,        &
                    & PTSTEP,PEXPL,                                    &
                    & TPFILE,ODIAG_IN_RUN,OOCEAN,                      &
                    & PSFUM,PSFVM,                                     &
                    & PTP,PRTKES,PRTHLS,PCOEF_DISS,PTDIFF,PTDISS,PRTKEMS,&
                    & TBUDGETS, KBUDGETS,                              &
                    & PEDR, PTR,PDISS, PCURRENT_TKE_DISS               )
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
!!    the parameter TURBN%CTURBLEN is set equal to KEPS. The same reasoning 
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
USE PARKIND1,   ONLY: JPRB
USE MODE_SHUMAN_PHY, ONLY: MZM_PHY, MZF_PHY, DZF_PHY, DZM_PHY
USE YOMHOOK,    ONLY: LHOOK, DR_HOOK
!
USE MODD_ARGSLIST_ll,    ONLY: LIST_ll
USE MODD_BUDGET,         ONLY: TBUDGETCONF_t, NBUDGET_TKE, NBUDGET_TH, TBUDGETDATA
USE MODD_CST,            ONLY: CST_t
USE MODD_CTURB,          ONLY: CSTURB_t
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_FIELD,          ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_LES,            ONLY: TLES_t
USE MODD_PARAMETERS,     ONLY: JPVEXT_TURB
USE MODD_TURB_n,         ONLY: TURB_t
!
USE MODE_BUDGET_PHY,         ONLY: BUDGET_STORE_ADD_PHY, BUDGET_STORE_END_PHY, BUDGET_STORE_INIT_PHY
USE MODE_IO_FIELD_WRITE_PHY, ONLY: IO_FIELD_WRITE_PHY
USE MODE_ll
!
USE MODI_GET_HALO
USE MODI_LES_MEAN_SUBGRID_PHY
USE MODE_TRIDIAG_TKE,          ONLY: TRIDIAG_TKE
!
!
IMPLICIT NONE
!
!
!*       0.1  declarations of arguments
!
!
TYPE(DIMPHYEX_t),        INTENT(IN)   :: D
TYPE(CST_t),             INTENT(IN)   :: CST
TYPE(CSTURB_t),          INTENT(IN)   :: CSTURB
TYPE(TBUDGETCONF_t),     INTENT(IN)   :: BUCONF
TYPE(TURB_t),            INTENT(IN)   :: TURBN
TYPE(TLES_t),            INTENT(INOUT):: TLES
CHARACTER(LEN=6),        INTENT(IN)   ::  HPROGRAM     ! CPROGRAM is the program currently running (modd_conf)
INTEGER,                 INTENT(IN)   ::  KMI          ! model index number  
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)   ::  PTKEM        ! TKE at t-deltat
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)   ::  PLM          ! mixing length         
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)   ::  PLEPS        ! dissipative length
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(INOUT)::  PDP          ! Dyn. prod. of TKE
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)   ::  PTRH
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)   ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)   ::  PDXX,PDYY,PDZZ,PDZX,PDZY ! metric coefficients
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)   ::  PZZ          ! physical height w-pt
REAL,                    INTENT(IN)   ::  PTSTEP       ! Time step 
REAL,                    INTENT(IN)   ::  PEXPL        ! Coef. temporal. disc.
TYPE(TFILEDATA),         INTENT(IN)   ::  TPFILE       ! Output file
LOGICAL,                 INTENT(IN)   ::  ODIAG_IN_RUN ! switch to activate online diagnostics (mesonh)
LOGICAL,                 INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
REAL, DIMENSION(D%NIJT),        INTENT(IN)    :: PSFUM,PSFVM ! momentum sfc flux
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)   ::  PTP          ! Ther. prod. of TKE
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(INOUT)::  PRTKES       ! RHOD * Jacobian * TKE at t+deltat
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(INOUT)::  PRTHLS       ! Source of Theta_l
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)   ::  PCOEF_DISS   ! 1/(Cph*Exner)
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT)  ::  PTDIFF       ! Diffusion TKE term
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT)  ::  PTDISS       ! Dissipation TKE term
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)   ::  PRTKEMS      ! Advection source
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER,                        INTENT(IN)   :: KBUDGETS
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT),   OPTIONAL ::  PEDR              ! EDR 
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT),   OPTIONAL ::  PTR               ! Transport prod. of TKE
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT),   OPTIONAL ::  PDISS             ! Dissipation of TKE
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(INOUT), OPTIONAL ::  PCURRENT_TKE_DISS ! if ODIAG_IN_RUN in mesonh
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT) ::         &
       ZA,       & ! under diagonal elements of the tri-diagonal matrix involved
                   ! in the temporal implicit scheme
       ZRES,     & ! treated variable at t+ deltat when the turbu-
                   ! lence is the only source of evolution added to the ones
                   ! considered in ZSOURCE. This variable is also used to
                   ! temporarily store some diagnostics stored in FM file
       ZFLX,     & ! horizontal or vertical flux of the treated variable
       ZSOURCE,  & ! source of evolution for the treated variable
       ZKEFF,    & ! effectif diffusion coeff = LT * SQRT( TKE )
       ZTR,      & ! Transport term
       ZMWORK1,ZMWORK2,& ! working var. for MZM/MZF operators (array syntax)
       ZDWORK1,ZDWORK2   ! working var. for DZM/DZF operators (array syntax)

LOGICAL,DIMENSION(D%NIJT,D%NKT) :: GTKENEG
                   ! 3D mask .T. if TKE < CSTURB%XTKEMIN
INTEGER             :: IIJB,IIJE,IKB,IKE,IKT,IKA,IKL  ! Index value for the mass points of the domain 
!
TYPE(LIST_ll), POINTER :: TZFIELDDISS_ll ! list of fields to exchange
INTEGER                :: IINFO_ll       ! return code of parallel routine
TYPE(TFIELDMETADATA)   :: TZFIELD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JIJ,JK
!
!----------------------------------------------------------------------------
NULLIFY(TZFIELDDISS_ll)
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
IF (LHOOK) CALL DR_HOOK('TKE_EPS_SOURCES',0,ZHOOK_HANDLE)
!
IKB=D%NKB
IKE=D%NKE
IIJB=D%NIJB
IIJE=D%NIJE
IKT=D%NKT
IKA=D%NKA
IKL=D%NKL
!
! compute the effective diffusion coefficient at the mass point
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZKEFF(:,:) = PLM(:,:) * SQRT(PTKEM(:,:))
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
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
IF (TURBN%CTURBDIM=='3DIM') THEN
  ZTR(:,:)=PTRH(:,:)
ELSE
  ZTR(:,:)=0.
END IF
!
!
!*       2.2  Explicit TKE sources except horizontal turbulent transport 
!
! extrapolate the dynamic production with a 1/Z law from its value at the 
IF (OOCEAN) THEN
  ! W(IKE) value stored in PDP(IKE) to the mass localization of tke(IKE)
  !$mnh_expand_array(JIJ=IIJB:IIJE)  
  PDP(:,IKE) = PDP(:,IKE) * (1. + PDZZ(:,IKE)/PDZZ(:,IKE+1))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)  
ELSE
  ! W(IKB+1) value stored in PDP(IKB) to the mass localization tke(IKB)
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  PDP(:,IKB) = PDP(:,IKB) * (1. + PDZZ(:,IKB+IKL)/PDZZ(:,IKB))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END IF
!
! Compute the source terms for TKE: ( ADVECtion + NUMerical DIFFusion + ..)
! + (Dynamical Production) + (Thermal Production) - (dissipation) 
!
CALL MZM_PHY(D,ZKEFF,ZMWORK1)
CALL MZM_PHY(D,PRHODJ,ZMWORK2)
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZFLX(:,:) = CSTURB%XCED * SQRT(PTKEM(:,:)) / PLEPS(:,:)
ZSOURCE(:,:) = ( PRTKES(:,:) +  PRTKEMS(:,:) ) &
                                     / PRHODJ(:,:) - PTKEM(:,:) / PTSTEP &
   + PDP(:,:) + PTP(:,:) + ZTR(:,:) & 
   - PEXPL * ZFLX(:,:) * PTKEM(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
!*       2.2  implicit vertical TKE transport
!
!
! To add here in ZSOURCE surface flux of TKE 
!(assumed to be 0 for ATM, 
IF (OOCEAN) THEN
  !for ocean:wave breaking  simple/very rough param wE = 100 Ustar**3 where ustar is the Tau_atmi/rhocea  
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZSOURCE(:,IKE)=ZSOURCE(:,IKE)-1.E2*((PSFUM(:)**2 + PSFVM(:)**2)**1.5) /PDZZ(:,IKE)  
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)  
END IF
! Compute the vector giving the elements just under the diagonal for the 
! matrix inverted in TRIDIAG 
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZA(:,:)     = - PTSTEP * CSTURB%XCET * ZMWORK1(:,:) & 
                                    * ZMWORK2(:,:) / PDZZ(:,:)**2
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
! Compute TKE at time t+deltat: ( stored in ZRES )
!
CALL TRIDIAG_TKE(D,PTKEM,ZA,PTSTEP,PEXPL,TURBN%XIMPL,PRHODJ,ZSOURCE,PTSTEP*ZFLX,ZRES)
CALL GET_HALO_PHY(D,ZRES)
!
!* diagnose the dissipation
!
IF (ODIAG_IN_RUN) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PCURRENT_TKE_DISS(:,:) = ZFLX(:,:) * PTKEM(:,:) &
                                  *(PEXPL*PTKEM(:,:) + TURBN%XIMPL*ZRES(:,:))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
  CALL ADD2DFIELD_ll(TZFIELDDISS_ll, PCURRENT_TKE_DISS, 'TKE_EPS_SOURCES::PCURRENT_TKE_DISS' )
  CALL UPDATE_HALO_ll(TZFIELDDISS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDDISS_ll)
ENDIF
!
! TKE must be greater than its minimum value
! CL : Now done at the end of the time step in ADVECTION_METSV for MesoNH
IF(HPROGRAM/='MESONH') THEN
 !$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
 GTKENEG(:,:) =  ZRES(:,:) <= CSTURB%XTKEMIN
 WHERE ( GTKENEG(:,:) ) 
   ZRES(:,:) = CSTURB%XTKEMIN
 END WHERE
 !$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PTDISS(:,:) = - ZFLX(:,:)*(PEXPL*PTKEM(:,:) &
                                      + TURBN%XIMPL*ZRES(:,:))
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
IF ( TLES%LLES_CALL .OR.                         &
     (TURBN%LTURB_DIAG .AND. TPFILE%LOPENED)  ) THEN
!
! Compute the cartesian vertical flux of TKE in ZFLX
!
  CALL MZM_PHY(D,ZKEFF,ZMWORK1)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZDWORK1(:,:) = TURBN%XIMPL * ZRES(:,:) + PEXPL * PTKEM(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL DZM_PHY(D,ZDWORK1,ZDWORK2)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZFLX(:,:)  = - CSTURB%XCET * ZMWORK1(:,:) &
                                     * ZDWORK2(:,:) / PDZZ(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
  ZFLX(:,IKB) = 0.
  ZFLX(:,IKA) = 0.
!
! Compute the whole turbulent TRansport of TKE:
!
  CALL MZM_PHY(D,PRHODJ,ZMWORK1)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZMWORK2(:,:) = ZMWORK1(:,:) * ZFLX(:,:) &
                                       / PDZZ(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL DZF_PHY(D,ZMWORK2,ZDWORK1)
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZTR(:,:)= ZTR(:,:) - ZDWORK1(:,:) &
                                  /PRHODJ(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
! Storage in the LES configuration
!
  IF (TLES%LLES_CALL) THEN
    CALL MZF_PHY(D,ZFLX,ZMWORK1)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZMWORK1, TLES%X_LES_SUBGRID_WTke )
    CALL LES_MEAN_SUBGRID_PHY(D,TLES, -ZTR, TLES%X_LES_SUBGRID_ddz_WTke )
  END IF
!
END IF
!
!*       2.4  stores the explicit sources for budget purposes
!
IF (BUCONF%LBUDGET_TKE) THEN
  ! Dynamical production
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZMWORK1(:,:) = PDP(:,:) * PRHODJ(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TKE), 'DP', ZMWORK1)
  !
  ! Thermal production
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZMWORK1(:,:) = PTP(:,:) * PRHODJ(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TKE), 'TP', ZMWORK1)
  !
  ! Dissipation
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZMWORK1(:,:) = -CSTURB%XCED * SQRT(PTKEM(:,:))/PLEPS(:,:) * &
                (PEXPL*PTKEM(:,:) + TURBN%XIMPL*ZRES(:,:))*PRHODJ(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT) 
  CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TKE), 'DISS',ZMWORK1)
END IF 
!
!*       2.5  computes the final RTKE and stores the whole turbulent transport
!              with the removal of the advection part for MesoNH

!Should be in IF LBUDGET_TKE only. Was removed out for a correct comput. of PTDIFF in case of LBUDGET_TKE=F in AROME
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
#ifdef REPRO48
IF (BUCONF%LBUDGET_TKE) THEN
PRTKES(:,:) = PRTKES(:,:) + PDP(:,:) * PRHODJ(:,:)
PRTKES(:,:) = PRTKES(:,:) + PTP(:,:) * PRHODJ(:,:)
PRTKES(:,:) = PRTKES(:,:) - CSTURB%XCED * SQRT(PTKEM(:,:)) / PLEPS(:,:) * &
                (PEXPL*PTKEM(:,:) + TURBN%XIMPL*ZRES(:,:)) * PRHODJ(:,:)
END IF
#else
PRTKES(:,:) = PRTKES(:,:) + PRHODJ(:,:) * &
                ( PDP(:,:) + PTP(:,:)                           &
                  - CSTURB%XCED * SQRT(PTKEM(:,:)) / PLEPS(:,:) &
                  * ( PEXPL*PTKEM(:,:) + TURBN%XIMPL*ZRES(:,:) ) )
#endif
!
PTDIFF(:,:) =  ZRES(:,:) / PTSTEP - PRTKES(:,:)&
                                     /PRHODJ(:,:) &
                           & - PDP(:,:)- PTP(:,:) - PTDISS(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
IF (BUCONF%LBUDGET_TKE) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TKE), 'TR', PRTKES)
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PRTKES(:,:) = ZRES(:,:) * PRHODJ(:,:) / PTSTEP &
                                    -  PRTKEMS(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!
! stores the whole turbulent transport
!
IF (BUCONF%LBUDGET_TKE) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TKE), 'TR', PRTKES)
!
!----------------------------------------------------------------------------
!
!*       3.   COMPUTE THE DISSIPATIVE HEATING
!             -------------------------------
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PRTHLS(:,:) = PRTHLS(:,:) + &
                                    CSTURB%XCED * SQRT(PTKEM(:,:)) / PLEPS(:,:) * &
                (PEXPL*PTKEM(:,:) + TURBN%XIMPL*ZRES(:,:)) &
                * PRHODJ(:,:) * PCOEF_DISS(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!----------------------------------------------------------------------------
!
!*       4.   STORES SOME DIAGNOSTICS
!             -----------------------
!
IF(PRESENT(PTR)) PTR=ZTR
IF(PRESENT(PDISS)) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PDISS(:,:) =  -CSTURB%XCED * (PTKEM(:,:)**1.5) / PLEPS(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
IF(PRESENT(PEDR)) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PEDR(:,:) = CSTURB%XCED * (PTKEM(:,:)**1.5) / PLEPS(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
END IF
!
IF ( TURBN%LTURB_DIAG .AND. TPFILE%LOPENED ) THEN
!
! stores the dynamic production 
!
  TZFIELD = TFIELDMETADATA(  &
    CMNHNAME   = 'DP',       &
    CSTDNAME   = '',         &
    CLONGNAME  = 'DP',       &
    CUNITS     = 'm2 s-3',   &
    CDIR       = 'XY',       &
    CCOMMENT   = 'X_Y_Z_DP', &
    NGRID      = 1,          &
    NTYPE      = TYPEREAL,   &
    NDIMS      = 3,          &
    LTIMEDEP   = .TRUE.      )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PDP)
!
! stores the thermal production 
!
  TZFIELD = TFIELDMETADATA(  &
    CMNHNAME   = 'TP',       &
    CSTDNAME   = '',         &
    CLONGNAME  = 'TP',       &
    CUNITS     = 'm2 s-3',   &
    CDIR       = 'XY',       &
    CCOMMENT   = 'X_Y_Z_TP', &
    NGRID      = 1,          &
    NTYPE      = TYPEREAL,   &
    NDIMS      = 3,          &
    LTIMEDEP   = .TRUE.      )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PTP)
!
! stores the whole turbulent transport
!
  TZFIELD = TFIELDMETADATA(  &
    CMNHNAME   = 'TR',       &
    CSTDNAME   = '',         &
    CLONGNAME  = 'TR',       &
    CUNITS     = 'm2 s-3',   &
    CDIR       = 'XY',       &
    CCOMMENT   = 'X_Y_Z_TR', &
    NGRID      = 1,          &
    NTYPE      = TYPEREAL,   &
    NDIMS      = 3,          &
    LTIMEDEP   = .TRUE.      )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,ZTR)
!
! stores the dissipation of TKE 
!
  TZFIELD = TFIELDMETADATA(    &
    CMNHNAME   = 'DISS',       &
    CSTDNAME   = '',           &
    CLONGNAME  = 'DISS',       &
    CUNITS     = 'm2 s-3',     &
    CDIR       = 'XY',         &
    CCOMMENT   = 'X_Y_Z_DISS', &
    NGRID      = 1,            &
    NTYPE      = TYPEREAL,     &
    NDIMS      = 3,            &
    LTIMEDEP   = .TRUE.        )
  CALL IO_FIELD_WRITE_PHY(D,TPFILE,TZFIELD,PDISS)
END IF
!
! Storage in the LES configuration of the Dynamic Production of TKE and
! the dissipation of TKE 
! 
IF (TLES%LLES_CALL ) THEN
  CALL LES_MEAN_SUBGRID_PHY(D,TLES, PDISS, TLES%X_LES_SUBGRID_DISS_Tke )
END IF
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TKE_EPS_SOURCES',1,ZHOOK_HANDLE)
END SUBROUTINE TKE_EPS_SOURCES
END MODULE MODE_TKE_EPS_SOURCES
