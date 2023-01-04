!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TURB_HOR  
IMPLICIT NONE
CONTAINS
             SUBROUTINE TURB_HOR(D,CST,CSTURB,TURBN,TLES,            &
                      KSPLT, KRR, KRRL, KRRI, PTSTEP,                &
                      KSV, KSV_LGBEG, KSV_LGEND, OFLAT,O2D,ONOMIXLG, &                      
                      OOCEAN,OCOMPUTE_SRC,OBLOWSNOW,                 &
                      TPFILE,                                        &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,PRSNOW,           &
                      PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,                 &
                      PCOSSLOPE,PSINSLOPE,                           &
                      PINV_PDXX, PINV_PDYY, PINV_PDZZ, PMZM_PRHODJ,  &
                      PK,                                            &
                      PRHODJ,PTHVREF,                                & 
                      PSFTHM,PSFRM,PSFSVM,                           &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU22M,PTAU33M,       &
                      PUM,PVM,PWM,PUSLOPEM,PVSLOPEM,PTHLM,PRM,PSVM,  &
                      PTKEM,PLM,PLEPS,                               &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,     &
                      PDP,PTP,PSIGS,                                 &
                      PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS               )
!     ################################################################
!
!
!!****  *TURB_HOR* -routine to compute the source terms in the meso-NH
!!               model equations due to the non-vertical turbulent fluxes.
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute the non-vertical
!     turbulent fluxes of the evolutive variables and give back the 
!     source terms to the main program.
!
!!**  METHOD
!!    ------
!!     Complementary 3D calculations when running at high resolution;
!!    The non-vertical turbulent fluxes are computed explicitly. The 
!!    contributions are cumulated in PRvarS and in DP and TP of TKE
!
! d(rho*T) = -d(rho*u'T'/dxx) -d(-rho*u'T'*dzx/dxx/dzz)
! / dt        / dx             /dz
!!    
!!
!!      Near the bottom of the model, uncentred evaluation of vertical 
!!    gradients are required because no field values are available under 
!!    the level where the gradient must be evaluated. In this case, the 
!!    gradient is computed with a second order accurate uncentred scheme 
!!    according to:
!!
!!        D FF           dzz3                       (dzz3+dzz4)   
!!        ----  = -  ----------------- FF(4)  +  ----------------- FF(3)   
!!        D z         (dzz3+dzz4) dzz4              dzz3 dzz4 
!!  
!!                    dzz4 + 2 dzz3          
!!                -  ----------------- FF(2)
!!                    (dzz3+dzz4) dzz3
!!
!!      where the values are taken from:
!!
!!                  -----    FF(5)
!!                    | 
!!                    |   dzz5
!!                    |    
!!                  -----    FF(4)
!!                    | 
!!                    |   dzz4
!!                    |    
!!                  -----    FF(3)
!!                    | 
!!                    |   dzz3
!!                    |    
!!                  -----    FF(2)    , (D FF / DZ)
!!                    |   dzz2 * 0.5
!!                  -----    ground
!!
!!
!!
!!    EXTERNAL
!!    --------
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
!!           JPVEXT     : number of vertical external points
!!
!!
!!           
!!
!!    REFERENCE
!!    ---------
!!      Book 2 of documentation (routine TURB_HOR)
!!      Book 1 of documentation (Chapter: Turbulence)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart             * INM and Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       Aug 29, 1994
!!      Modifications: Feb 14, 1995 (J.Cuxart and J.Stein) 
!!                                  Doctorization and Optimization
!!                     March 21, 1995 (J.M. Carriere) 
!!                                  Introduction of cloud water
!!                     June  14, 1995 (J. Stein) 
!!                                  rm the ZVTPV computation + bug in the all 
!!                                  or nothing condens. case
!!                     June 28, 1995 (J.Cuxart)  Add the LES tools 
!!                     Sept 19, 1995 (J. Stein) change the surface flux
!!               computations
!!                     Nov  13, 1995 (J. Stein) include the tangential fluxes
!!               bug in <u'w'> at the surface
!!                     Nov  27, 1997 (V. Saravane) spliting of the routine
!!                     Nov  27, 1997 (V. Masson) clearing of the routine
!!                     Nov  06, 2002 (V. Masson) LES budgets
!!                     Feb  20, 2003 (JP Pinty)  Add PFRAC_ICE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST, ONLY : CST_t
USE MODD_CTURB, ONLY : CSTURB_t
USE MODD_TURB_n, ONLY: TURB_t
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_IO, ONLY: TFILEDATA
USE MODD_PARAMETERS
USE MODD_LES, ONLY: TLES_t
!
USE MODE_TURB_HOR_THERMO_FLUX, ONLY: TURB_HOR_THERMO_FLUX
USE MODE_TURB_HOR_THERMO_CORR, ONLY: TURB_HOR_THERMO_CORR
USE MODE_TURB_HOR_DYN_CORR, ONLY: TURB_HOR_DYN_CORR
USE MODE_TURB_HOR_UV, ONLY: TURB_HOR_UV
USE MODE_TURB_HOR_UW, ONLY: TURB_HOR_UW
USE MODE_TURB_HOR_VW, ONLY: TURB_HOR_VW
USE MODE_TURB_HOR_SV_FLUX, ONLY: TURB_HOR_SV_FLUX
USE MODE_TURB_HOR_SV_CORR, ONLY: TURB_HOR_SV_CORR
!
IMPLICIT NONE
!
!
!*       0.1  declaration of arguments
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
TYPE(TURB_t),           INTENT(IN)   :: TURBN
TYPE(TLES_t),           INTENT(INOUT):: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KSPLT         ! current split index
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV,KSV_LGBEG,KSV_LGEND ! number of sv var.
REAL,                   INTENT(IN)   ::  PTSTEP       !
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
LOGICAL,                INTENT(IN)   ::  ONOMIXLG     ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL,                INTENT(IN)   ::  O2D          ! Logical for 2D model version (modd_conf)
LOGICAL,                INTENT(IN)   ::  OOCEAN ! switch for ocean version
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables
LOGICAL,                INTENT(IN)   ::  OBLOWSNOW    ! switch to activate pronostic blowing snow

TYPE(TFILEDATA),          INTENT(IN)    ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PZZ          ! vertical grid
REAL,                     INTENT(IN)    ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
REAL, DIMENSION(:,:),     INTENT(IN)    ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle 
                                      ! between i and the slope vector

REAL, DIMENSION(:,:,:),   INTENT(IN) :: PK           ! Turbulent diffusion doef.
                                                     ! PK = PLM * SQRT(PTKEM)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PINV_PDXX    ! 1./PDXX
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PINV_PDYY    ! 1./PDYY
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PINV_PDZZ    ! 1./PDZZ
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PMZM_PRHODJ  ! MZM(PRHODJ)

REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRHODJ       ! density * grid volume
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHVREF      ! ref. state VPT       
!
REAL, DIMENSION(:,:),     INTENT(IN)    ::  PSFTHM,PSFRM
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PSFSVM       ! surface fluxes
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCDUEFF      ! Cd * || u || at time t
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU22M      ! <vv> in the same axes
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
! Variables at t-1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PUM,PVM,PWM,PTHLM 
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PRM          ! mixing ratios at t-1,
                              !  where PRM(:,:,:,1) = conservative mixing ratio
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PSVM         ! scalar var. at t-1
REAL, DIMENSION(:,:),      INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(:,:),      INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTKEM        ! TKE at time t- dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLM          ! Turb. mixing length
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLEPS        ! dissipative length
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLOCPEXNM    ! Lv(T)/Cp/Exner at time t-1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PAMOIST      ! s and Thetal and Rnp

REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PSRCM
                                  ! normalized 2nd-order flux
                                  ! s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PFRAC_ICE    ! ri fraction of rc+ri
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRUS, PRVS, PRWS, PRTHLS
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS,PRRS   ! var. at t+1 -split-
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PDP,PTP      ! TKE production terms
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PSIGS
                                  ! IN: Vertical part of Sigma_s at t
                                  ! OUT: Total Sigma_s at t
!
!
!
!*       0.2  declaration of local variables
!
! ---------------------------------------------------------------------------
!
!*       1.   PRELIMINARY COMPUTATIONS
!             ------------------------
!
!* Exchange coefficient is limited in order to insure numerical stability
!
!!
!*       2.   < U' THETA'l >
!*       3.   < U' R'np >
!*       4.   < U' TPV' >
!*       5.   < V' THETA'l >
!*       6.   < V' R'np >
!*       7.   < V' TPV' >
!
      CALL      TURB_HOR_THERMO_FLUX(TURBN,TLES,KSPLT, KRR, KRRL, KRRI,&
                      TPFILE,OFLAT, O2D,                             &
                      PK,PINV_PDXX,PINV_PDYY,PINV_PDZZ,PMZM_PRHODJ,  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      PDIRCOSXW,PDIRCOSYW,                           &
                      PRHODJ,                                        &
                      PSFTHM,PSFRM,                                  &
                      PWM,PTHLM,PRM,                                 &
                      PATHETA,PAMOIST,PSRCM,PFRAC_ICE,               &
                      PRTHLS,PRRS                                    )
!
!
!*       8.   TURBULENT CORRELATIONS : <THl THl>, <THl Rnp>, <Rnp Rnp>, Sigma_s
!
      IF (KSPLT==1)                                                  &
      CALL      TURB_HOR_THERMO_CORR(D,CST,TURBN,TLES,               &
                      KRR, KRRL, KRRI,                               &
                      OOCEAN,OCOMPUTE_SRC,O2D,                       &
                      TPFILE,                                        &
                      PINV_PDXX,PINV_PDYY,                           &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      PTHVREF,                                       &
                      PWM,PTHLM,PRM,                                 &
                      PTKEM,PLM,PLEPS,                               &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,               & 
                      PSIGS                                          )
!
!
!*       9.   < U'U'>
!*      10.   < V'V'>
!*      11.   < W'W'>
! 
      CALL       TURB_HOR_DYN_CORR(TURBN,TLES,KSPLT, PTSTEP,         &
                      KRR,KSV,OFLAT, O2D,                            &
                      TPFILE,                                        &
                      PK,PINV_PDZZ,                                  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,                  &
                      PDIRCOSZW,                                     &
                      PCOSSLOPE,PSINSLOPE,                           &
                      PRHODJ,                                        &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU22M,PTAU33M,       &
                      PUM,PVM,PWM, PUSLOPEM,PVSLOPEM,                &
                      PTHLM,PRM,PSVM,                                &
                      PTKEM,PLM,                                     &
                      PDP,PTP,                                       &
                      PRUS,PRVS,PRWS                                 )
!
!
!*      12.   < U'V'>
!
      CALL      TURB_HOR_UV(TURBN,TLES,KSPLT,OFLAT,O2D,              &
                      TPFILE,                                        &
                      PK,PINV_PDXX,PINV_PDYY,PINV_PDZZ,PMZM_PRHODJ,  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,                      &
                      PDIRCOSZW,                                     &
                      PCOSSLOPE,PSINSLOPE,                           &
                      PRHODJ,                                        &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU22M,PTAU33M,       &
                      PUM,PVM,PUSLOPEM,PVSLOPEM,                     &
                      PDP,                                           &
                      PRUS,PRVS                                      )
!
!
!*      13.   < U'W'>
!
      CALL      TURB_HOR_UW(TURBN,TLES,KSPLT,                        &
                      KRR,KSV,OFLAT,                                 &
                      TPFILE,                                        &
                      PK,PINV_PDXX,PINV_PDZZ,PMZM_PRHODJ,            &
                      PDXX,PDZZ,PDZX,                                &
                      PRHODJ,PTHVREF,                                &
                      PUM,PWM,PTHLM,PRM,PSVM,                        &
                      PTKEM,PLM,                                     &
                      PDP,                                           &
                      PRUS,PRWS                                      )
!
!
!*      14.   < V'W'>
!
      CALL      TURB_HOR_VW(TURBN,TLES,KSPLT,                        &
                      KRR,KSV,OFLAT,O2D,                             &
                      TPFILE,                                        &
                      PK,PINV_PDYY,PINV_PDZZ,PMZM_PRHODJ,            &
                      PDYY,PDZZ,PDZY,                                &
                      PRHODJ,PTHVREF,                                &
                      PVM,PWM,PTHLM,PRM,PSVM,                        &
                      PTKEM,PLM,                                     &
                      PDP,                                           &
                      PRVS,PRWS                                      )

!
!
!*      15.   HORIZONTAL FLUXES OF PASSIVE SCALARS
!
      CALL      TURB_HOR_SV_FLUX(TURBN,TLES,KSPLT,OBLOWSNOW,OFLAT,   &
                      TPFILE,KSV_LGBEG,KSV_LGEND,O2D,ONOMIXLG,       &
                      PK,PINV_PDXX,PINV_PDYY,PINV_PDZZ,PMZM_PRHODJ,  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PRSNOW,               &
                      PDIRCOSXW,PDIRCOSYW,                           &
                      PRHODJ,PWM,                                    &
                      PSFSVM,                                        &
                      PSVM,                                          &
                      PRSVS                                          )
!
      IF (KSPLT==1 .AND. TLES%LLES_CALL)                             &
      CALL      TURB_HOR_SV_CORR(D,CST,CSTURB,TLES,                  &
                      KSV,KSV_LGBEG,KSV_LGEND,                       &
                      KRR,KRRL,KRRI,OOCEAN,OCOMPUTE_SRC,OBLOWSNOW,   &
                      ONOMIXLG,O2D,                                  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PRSNOW,               &
                      PLM,PLEPS,PTKEM,PTHVREF,                       &
                      PTHLM,PRM,                                     &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,               &
                      PWM,PSVM                                       )
!
!
END SUBROUTINE TURB_HOR
END MODULE MODE_TURB_HOR
