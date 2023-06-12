!MNH_LIC Copyright 2002-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################
      MODULE MODI_ICE_ADJUST_ELEC
!     ###########################
!
IMPLICIT NONE
INTERFACE
!
      SUBROUTINE ICE_ADJUST_ELEC (KRR, KMI, HRAD, HTURBDIM, HSCONV, HMF_CLOUD,    &
                                  OSUBG_COND, OSIGMAS, PTSTEP,PSIGQSAT,           &
                                  PRHODJ, PEXNREF, PSIGS, PPABST, PZZ,            &
                                  PMFCONV, PCF_MF, PRC_MF, PRI_MF,                &
                                  PRVT, PRCT, PRVS, PRCS, PTHS, PSRCS, PCLDFR ,   &
                                  PRRT, PRRS, PRIT, PRIS, PRST, PRSS, PRGT, PRGS, &
                                  PQPIT, PQPIS, PQCT, PQCS,                       &
                                  PQRT, PQRS, PQIT, PQIS, PQST, PQSS, PQGT, PQGS, &
                                  PQNIT, PQNIS, PRHT, PRHS, PQHT, PQHS            )
IMPLICIT NONE
!
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KMI      ! Model index 
CHARACTER(len=4),         INTENT(IN)    :: HTURBDIM ! Dimensionality of the
                                                    ! turbulence scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSCONV   ! Shallow convection scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HMF_CLOUD! Type of statistical cloud
CHARACTER(len=4),         INTENT(IN)    :: HRAD     ! Radiation scheme name
LOGICAL,                  INTENT(IN)    :: OSUBG_COND ! Switch for Subgrid 
                                                    ! Condensation
LOGICAL,                  INTENT(IN)    :: OSIGMAS  ! Switch for Sigma_s: 
                                                    ! use values computed in CONDENSATION
                                                    ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PTSTEP    ! Double Time step
                                                    ! (single if cold start)
REAL,                     INTENT(IN)   :: PSIGQSAT  ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF ! Reference Exner function
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PSIGS   ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST  ! Absolute Pressure at t        
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PZZ     ! height of model layer
!
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRI_MF! Convective Mass Flux solid mixing ratio
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                   ! s'rc'/2Sigma_s2 at time t+1
                                                   ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PCLDFR  ! Cloud fraction          
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRRS ! Rain water m.r. at t+1
REAL, DIMENSION(:,:,:),   INTENT(INOUT)::  PRIS ! Cloud ice  m.r. at t+1
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRSS ! Aggregate  m.r. at t+1
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRGS ! Graupel    m.r. at t+1
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRRT ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRIT ! Cloud ice  m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRST ! Aggregate  m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRGT ! Graupel    m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) ::  PRHT ! Hail m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) ::  PRHS ! Hail m.r. at t+1
!
REAL, DIMENSION(:,:,:),   INTENT(IN)      ::  PQPIT  ! positive ion m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   ::  PQPIS  !           source
REAL, DIMENSION(:,:,:),   INTENT(IN)      ::  PQNIT  ! negative ion m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   ::  PQNIS  ! source
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PQCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PQCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PQRS ! Rain water m.r. at t+1
REAL, DIMENSION(:,:,:),   INTENT(INOUT)::  PQIS ! Cloud ice  m.r. at t+1
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PQSS ! Aggregate  m.r. at t+1
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PQGS ! Graupel    m.r. at t+1
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PQRT ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PQIT ! Cloud ice  m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PQST ! Aggregate  m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PQGT ! Graupel    m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: PQHT ! Hail  m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: PQHS ! Hail  m.r. at t+1
!
END SUBROUTINE ICE_ADJUST_ELEC
END INTERFACE
END MODULE MODI_ICE_ADJUST_ELEC
!
!     ########################################################################
      SUBROUTINE ICE_ADJUST_ELEC (KRR, KMI, HRAD, HTURBDIM, HSCONV,          &
                             HMF_CLOUD, OSUBG_COND, OSIGMAS, PTSTEP,PSIGQSAT,&
                             PRHODJ, PEXNREF, PSIGS, PPABST, PZZ,            &
                             PMFCONV, PCF_MF, PRC_MF, PRI_MF,                &  
                             PRVT, PRCT, PRVS, PRCS, PTHS, PSRCS, PCLDFR ,   &
                             PRRT, PRRS, PRIT, PRIS, PRST, PRSS, PRGT, PRGS, &
                             PQPIT, PQPIS, PQCT, PQCS,                       &
                             PQRT, PQRS, PQIT, PQIS, PQST, PQSS, PQGT, PQGS, &
                             PQNIT, PQNIS, PRHT, PRHS, PQHT, PQHS            )
!     ########################################################################
!
!!****  *ICE_ADJUST_ELEC* -  compute the ajustment of water vapor in mixed-phase 
!!                           clouds
!!
!!    PURPOSE
!!    -------
!!    The purpose of this routine is to compute the fast microphysical sources
!!    through a saturation ajustement procedure in case of mixed-phase clouds.
!!
!!
!!**  METHOD
!!    ------
!!    Langlois, Tellus, 1973 for the cloudless version.
!!    When cloud water is taken into account, refer to book 1 of the
!!    documentation.
!!
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!         XP00               ! Reference pressure
!!         XMD,XMV            ! Molar mass of dry air and molar mass of vapor
!!         XRD,XRV            ! Gaz constant for dry air, gaz constant for vapor
!!         XCPD,XCPV          ! Cpd (dry air), Cpv (vapor)
!!         XCL                ! Cl (liquid)
!!         XCI                ! Ci (ice)
!!         XTT                ! Triple point temperature
!!         XLVTT              ! Vaporization heat constant
!!         XLSTT              ! Sublimation  heat constant
!!         XALPW,XBETAW,XGAMW ! Constants for saturation vapor over liquid
!!                            !  pressure  function 
!!         XALPI,XBETAI,XGAMI ! Constants for saturation vapor over ice
!!                            !  pressure  function 
!!      Module  MODD_CONF 
!!         CCONF
!!
!!
!!    REFERENCE
!!    ---------
!!      Book 1 and Book2 of documentation ( routine ICE_ADJUST )
!!      Langlois, Tellus, 1973
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty    * Laboratoire d'Aerologie*
!!   
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    2002 
!!      C. Barthe   19/11/09   update to version 4.8.1
!!      M. Chong    Mar. 2010  Add small ions
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet    03/2020: use the new data structures and subroutines for budgets
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,         only: lbudget_th, lbudget_rv, lbudget_rc, lbudget_ri, lbudget_sv,  &
                               NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RI, NBUDGET_SV1, &
                               tbudgets
USE MODD_CONF
USE MODD_CST
USE MODD_ELEC_DESCR, ONLY : XRTMIN_ELEC, XQTMIN, XFC, XFI, XECHARGE
USE MODD_NSV, ONLY : NSV_ELECBEG, NSV_ELECEND
USE MODD_PARAMETERS
USE MODD_RAIN_ICE_DESCR_n, ONLY : XRTMIN, XBI
USE MODD_RAIN_ICE_PARAM_n,   ONLY: RAIN_ICE_PARAMN
USE MODD_NEB_n,            ONLY: NEBN
USE MODD_TURB_n,           ONLY: TURBN
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_t

use mode_budget,          only: Budget_store_init, Budget_store_end
use mode_tools_ll,        only: GET_INDICE_ll
USE MODE_FILL_DIMPHYEX,   ONLY: FILL_DIMPHYEX

USE MODI_CONDENSATION
USE MODI_GET_HALO
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments 
!
INTEGER,                INTENT(IN)    :: KRR      ! Number of moist variables
INTEGER,                INTENT(IN)    :: KMI      ! Model index 
CHARACTER(len=4),       INTENT(IN)    :: HTURBDIM ! Dimensionality of the
                                                  ! turbulence scheme
CHARACTER(LEN=4),       INTENT(IN)    :: HSCONV   ! Shallow convection scheme
CHARACTER(LEN=4),       INTENT(IN)    :: HMF_CLOUD! Type of statistical cloud
CHARACTER(len=4),       INTENT(IN)    :: HRAD     ! Radiation scheme name
LOGICAL,                INTENT(IN)    :: OSUBG_COND ! Switch for Subgrid 
                                                    ! Condensation
LOGICAL,                INTENT(IN)    :: OSIGMAS  ! Switch for Sigma_s: 
                                                  ! use values computed in CONDENSATION
                                                  ! or that from turbulence scheme
REAL,                   INTENT(IN)    :: PTSTEP   ! Double Time step
                                                  ! (single if cold start)
REAL,                     INTENT(IN)   :: PSIGQSAT  ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PSIGS   ! Sigma_s at time t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PPABST  ! Absolute Pressure at t        
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PZZ     ! height of model layer
!
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRI_MF! Convective Mass Flux solid mixing ratio
!
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                 ! s'rc'/2Sigma_s2 at time t+1
                                                 ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PCLDFR  ! Cloud fraction          
!
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRRS ! Rain water m.r. at t+1
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRIS ! Cloud ice  m.r. at t+1
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRSS ! Aggregate  m.r. at t+1
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRGS ! Graupel    m.r. at t+1
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRRT ! Rain water m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRIT ! Cloud ice  m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRST ! Aggregate  m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRGT ! Graupel    m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: PRHT ! Hail m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: PRHS ! Hail m.r. at t+1
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   ::  PQPIT  ! positive ion m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   ::  PQPIS  !           source
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   ::  PQNIT  ! negative ion m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   ::  PQNIS  ! source
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PQCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQCS ! Cloud water m.r. source
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PQRS ! Rain water m.r. at t+1
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQIS ! Cloud ice  m.r. at t+1
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PQSS ! Aggregate  m.r. at t+1
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PQGS ! Graupel    m.r. at t+1
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PQRT ! Rain water m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PQIT ! Cloud ice  m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PQST ! Aggregate  m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PQGT ! Graupel    m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: PQHT ! Hail m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: PQHS ! Hail m.r. at t+1
!
!
!*       0.2   Declarations of local variables :
!
REAL  :: ZEPS  ! Mv/Md
REAL  :: ZT00,ZT0   ! Min and max temperature for the mixed phase liquid and solid water
                    ! for the coeff CND of the barycentric mixing ratio
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) &
                         :: ZEXNS,&  ! guess of the Exner functon at t+1
                            ZT,   &  ! guess of the temperature at t+1
                            ZCPH, &  ! guess of the CPh for the mixing
                            ZLV,  &  ! guess of the Lv at t+1
                            ZLS,  &  ! guess of the Ls at t+1
      ZW1,ZW2,ZW3,ZW4,ZW5,ZW6,ZW7,&  ! Work arrays for intermediate fields
    ZW1_IN, ZW2_IN, ZW3_IN, ZDUM, &
                            ZCND     ! CND=(T-T00)/(T0-T00) cf sc doc and TAO etal (89)
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) &
                         :: ZWE1, &
                            ZWE2
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) &
                         :: ZION_NUMBER,  &  !nearly Nb of elementary charge
                                             ! in hydrometeor charge  
                            ZADD             ! ratio (0 or 1) of ZION_NUMBER 
                                             ! to add to positive
                                             ! or negative ion number
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2)) :: ZSIGQSAT2D
                                             !
INTEGER             :: IIU,IJU,IKU! dimensions of dummy arrays
INTEGER             :: IIB,IJB    ! Horz index values of the first inner mass points
INTEGER             :: IIE,IJE    ! Horz index values of the last inner mass points
INTEGER             :: IKB        ! K index value of the first inner mass point
INTEGER             :: IKE        ! K index value of the last inner mass point
INTEGER             :: JITER,ITERMAX ! iterative loop for first order adjustment
!
LOGICAL             :: LPRETREATMENT, LNEW_ADJUST
!
TYPE(DIMPHYEX_t)    :: D
!
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
IIU = SIZE(PEXNREF,1)
IJU = SIZE(PEXNREF,2)
IKU = SIZE(PEXNREF,3)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
CALL FILL_DIMPHYEX(D, IIU, IJU, IKU)
!
ZEPS = XMV / XMD
!
ITERMAX = 1
!
LPRETREATMENT=.TRUE.     ! FALSE to retreive the previous MASDEV4_1 version
LNEW_ADJUST  =.TRUE.     ! FALSE to retreive the previous MASDEV4_1 version
ZT0  = XTT               ! Usefull if LPRETREATMENT=T or LNEW_ADJUST=T
ZT00 = XTT-40.           ! Usefull if LPRETREATMENT=T or LNEW_ADJUST=T
!
!-------------------------------------------------------------------------------
if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'DEPI', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'DEPI', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'DEPI', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'DEPI', pris(:, :, :) * prhodj(:, :, :) )
if ( lbudget_sv ) then
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg     ), 'DEPI', pqpis(:, :, :) * prhodj(:, :, :) )
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecend     ), 'DEPI', pqnis(:, :, :) * prhodj(:, :, :) )
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'DEPI', pqcs (:, :, :) * prhodj(:, :, :) )
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'DEPI', pqis (:, :, :) * prhodj(:, :, :) )
end if
!
!*       2.     COMPUTE QUANTITIES WITH THE GUESS OF THE FUTURE INSTANT
!               -------------------------------------------------------
!
!*       2.1    estimate the pressure at t+1
!
ZEXNS(:,:,:) = ( PPABST(:,:,:)  / XP00)**(XRD/XCPD)
!
!    beginning of the iterative loop
!
DO JITER = 1, ITERMAX
!
!*       2.2    compute the intermediate temperature at t+1, T*
!  
  ZT(:,:,:) = (PTHS(:,:,:) * PTSTEP) * ZEXNS(:,:,:)
!
!*       2.3    compute the latent heat of vaporization Lv(T*) at t+1
!                   and the latent heat of sublimation  Ls(T*) at t+1
!
  ZLV(:,:,:) = XLVTT + (XCPV - XCL) * (ZT(:,:,:) - XTT)
  ZLS(:,:,:) = XLSTT + (XCPV - XCI) * (ZT(:,:,:) - XTT)
!
!*       2.4    compute the specific heat for moist air (Cph) at t+1
!
  IF     ( KRR == 7 ) THEN
    ZCPH(:,:,:) = XCPD + XCPV *PTSTEP*  PRVS(:,:,:)                             &
                       + XCL  *PTSTEP* (PRCS(:,:,:) + PRRS(:,:,:))              &
                       + XCI  *PTSTEP* (PRIS(:,:,:) + PRSS(:,:,:) + PRGS(:,:,:) &
                                                                  + PRHS(:,:,:))
  ELSE IF( KRR == 6 ) THEN
    ZCPH(:,:,:) = XCPD + XCPV *PTSTEP*  PRVS(:,:,:)                            &
                       + XCL  *PTSTEP* (PRCS(:,:,:) + PRRS(:,:,:))             &
                       + XCI  *PTSTEP* (PRIS(:,:,:) + PRSS(:,:,:) + PRGS(:,:,:))
  ELSE IF( KRR == 5 ) THEN
    ZCPH(:,:,:) = XCPD + XCPV *PTSTEP*  PRVS(:,:,:)                            &
                       + XCL  *PTSTEP* (PRCS(:,:,:) + PRRS(:,:,:))             &
                       + XCI  *PTSTEP* (PRIS(:,:,:) + PRSS(:,:,:))
  ELSE IF( KRR == 3 ) THEN
    ZCPH(:,:,:) = XCPD + XCPV *PTSTEP*  PRVS(:,:,:)              &
                       + XCL  *PTSTEP* (PRCS(:,:,:) + PRRS(:,:,:))
  ELSE IF( KRR == 2 ) THEN
    ZCPH(:,:,:) = XCPD + XCPV *PTSTEP*  PRVS(:,:,:) &
                       + XCL  *PTSTEP*  PRCS(:,:,:)
  END IF
!
!
!*       3.     FIRST ORDER SUBGRID CONDENSATION SCHEME
!               ---------------------------------------
!
  IF (OSUBG_COND) THEN
!
!*       3.1    compute condensate, cloud fraction
!
  !   ZW3=water vapor    ZW1=rc (OUT)  ZW2=ri (OUT)   PSRC= s'rci'/Sigma_s^2
  ! ZW3_IN/ZW2_IN/ZW1_IN (IN)
    ZW3_IN = PRVS * PTSTEP;     ZW1_IN = PRCS * PTSTEP;  ZW2_IN = PRIS * PTSTEP
    ZW3=ZW3_IN; ZW2=ZW2_IN; ZW1=ZW1_IN 
    ZSIGQSAT2D(:,:)=PSIGQSAT
    ZW4 = 1. ! PRODREF is not used if HL variables are not present
!
    CALL CONDENSATION(D, CST, RAIN_ICE_PARAMN, NEBN, TURBN, &
                     &'T', 'CB02', 'CB',                                                  &
                     &PPABST, PZZ, ZW4, ZT, ZW3_IN, ZW3, ZW1_IN, ZW1, ZW2_IN, ZW2,    &
                     &PRRS*PTSTEP, PRSS*PTSTEP, PRGS*PTSTEP, PSIGS, .FALSE., PMFCONV, PCLDFR, PSRCS, .FALSE.,                 &
                     &OSIGMAS, .FALSE.,                                                                 &
                     &ZDUM, ZDUM, ZDUM, ZDUM, ZDUM, ZSIGQSAT2D, &
                     &ZLV, ZLS, ZCPH)
!
!*       3.2    compute the variation of mixing ratio
!
                                                        !         Rc - Rc*
    ZW1(:,:,:) = (ZW1(:,:,:) / PTSTEP) - PRCS(:,:,:)    ! Pcon = ----------
                                                        !         2 Delta t

    ZW2(:,:,:) = (ZW2(:,:,:) / PTSTEP) - PRIS(:,:,:)    ! idem ZW1 but for Ri

  ELSE
!
!
!*       4.     SECOND ORDER ALL OR NOTHING CONDENSATION SCHEME
!                            FOR MIXED-PHASE CLOUD
!               -----------------------------------------------
!
!
!*       4.1    Eventually pretreatment
!
    IF (LPRETREATMENT) THEN
!
!     compute the saturation vapor pressures at t+1
!
      CALL GET_HALO(ZT)
      ZW1(:,:,:) = EXP(XALPW - XBETAW/ZT(:,:,:) - XGAMW*ALOG(ZT(:,:,:))) ! e_sw
      ZW2(:,:,:) = EXP(XALPI - XBETAI/ZT(:,:,:) - XGAMI*ALOG(ZT(:,:,:))) ! e_si
      ZW1(:,:,:) = MIN(PPABST(:,:,:)/2.,ZW1(:,:,:))   ! safety limitation
      ZW2(:,:,:) = MIN(PPABST(:,:,:)/2.,ZW2(:,:,:))   ! safety limitation
!
!     compute the saturation mixing ratios at t+1
!
      ZW3(:,:,:) = ZW1(:,:,:) * ZEPS / &
                  ( PPABST(:,:,:) - ZW1(:,:,:))  ! r_sw
      ZW4(:,:,:) = ZW2(:,:,:) * ZEPS / &
                  ( PPABST(:,:,:) - ZW2(:,:,:))  ! r_si
!
      WHERE(PRVS(:,:,:)*PTSTEP .LT. ZW4(:,:,:) .AND. &
            PRCS(:,:,:) .GT. 0. .AND. ZT(:,:,:) .LT. XTT)
!
!       Subsaturation case with respect to rsi(T,P) (and case rv<0):
!       Evaporation of rc>0 (while enough) to decrease the lack of vapor 
!
        ZW5 (:,:,:)= MIN( PRCS , ZW4(:,:,:)/PTSTEP - PRVS(:,:,:) )                ! RVCNDC
        PRVS(:,:,:)= PRVS(:,:,:) + ZW5(:,:,:) 
        PRCS(:,:,:)= PRCS(:,:,:) - ZW5(:,:,:) 
        PTHS(:,:,:)= PTHS(:,:,:) - ZW5(:,:,:) * ZLV(:,:,:) /(ZCPH(:,:,:)*PEXNREF(:,:,:))
!
      END WHERE
!  
      WHERE (PRVS(:,:,:)*PTSTEP .GT. ZW3(:,:,:))
!
!       Supersaturation case with respect to rsw(T,P):
!       Condensation of the vapor that is left 
!
        ZW5 (:,:,:)= PRVS(:,:,:) - ZW3(:,:,:)/PTSTEP 
        PRVS(:,:,:)= PRVS(:,:,:) - ZW5(:,:,:)                                  ! RVCNDC 
        PRCS(:,:,:)= PRCS(:,:,:) + ZW5(:,:,:) 
        PTHS(:,:,:)= PTHS(:,:,:) + ZW5(:,:,:) * ZLV(:,:,:) /(ZCPH(:,:,:)*PEXNREF(:,:,:))
!
      END WHERE
!
      WHERE (PRCS(:,:,:) .GT. 0. .AND. ZT(:,:,:) .LT. ZT00)
!
!       Treatment of rc>0 if T<T00:
!
        PRIS(:,:,:)= PRIS(:,:,:) + PRCS(:,:,:) 
        PTHS(:,:,:)= PTHS(:,:,:) + PRCS(:,:,:) * &
                    (ZLS(:,:,:) - ZLV(:,:,:)) / (ZCPH(:,:,:) * PEXNREF(:,:,:))
        PRCS(:,:,:)= 0. 
!
      END WHERE
!
!*       4.2    compute the intermediate temperature at t+1, T*
!  
      ZT(:,:,:) = (PTHS(:,:,:) * PTSTEP) * ZEXNS(:,:,:)
!
    END IF  !end PRETREATMENT 
!
!*       4.3    compute the saturation vapor pressures at t+1
!
    ZW1(:,:,:) = EXP(XALPW - XBETAW / ZT(:,:,:) - XGAMW * ALOG(ZT(:,:,:))) ! e_sw
    ZW2(:,:,:) = EXP(XALPI - XBETAI / ZT(:,:,:) - XGAMI * ALOG(ZT(:,:,:))) ! e_si
    ZW1(:,:,:) = MIN(PPABST(:,:,:)/2.,ZW1(:,:,:))   ! safety limitation
    ZW2(:,:,:) = MIN(PPABST(:,:,:)/2.,ZW2(:,:,:))   ! safety limitation
!
!*       4.4    compute the saturation mixing ratios at t+1
!
    ZW3(:,:,:) =  ZW1(:,:,:) * ZEPS     &
               / ( PPABST(:,:,:)  - ZW1(:,:,:) ) ! r_sw
    ZW4(:,:,:) =  ZW2(:,:,:) * ZEPS     &
               / ( PPABST(:,:,:)  - ZW2(:,:,:) ) ! r_si
!
!*       4.5    compute the saturation mixing ratio derivatives (r'_vs)
!
    ZW1(:,:,:) = (( XBETAW/ZT(:,:,:) - XGAMW ) / ZT(:,:,:)) & ! r'_sw
                    * ZW3(:,:,:) * ( 1. + ZW3(:,:,:)/ZEPS )
    ZW2(:,:,:) = (( XBETAI/ZT(:,:,:) - XGAMI ) / ZT(:,:,:)) & ! r'_si
                    * ZW4(:,:,:) * ( 1. + ZW4(:,:,:)/ZEPS )
!
    IF (LNEW_ADJUST) THEN
      ZCND(:,:,:)= (ZT(:,:,:) - ZT00) / (ZT0 - ZT00)       ! Like Tao et al 89
      ZCND(:,:,:)= MAX ( MIN(ZCND(:,:,:),1.) , 0. )
    ELSE
      WHERE ((PRCS(:,:,:)+PRIS(:,:,:)) .GT. 1.0E-20)      &
      ZCND(:,:,:)= PRCS(:,:,:) / (PRCS(:,:,:) + PRIS(:,:,:)) ! Like the original version
    END IF
!
!*       4.5    compute L_v CND + L_s DEP and F'(T)
!
    WHERE ((PRCS(:,:,:)+PRIS(:,:,:)) .GT. 1.0E-20)
!
      ZW5(:,:,:) = ZLS(:,:,:) + (ZLV(:,:,:) - ZLS(:,:,:)) * ZCND(:,:,:)
      ZW6(:,:,:) = ZCPH(:,:,:) * (PRCS(:,:,:) + PRIS(:,:,:)) +          &
                   ZW5(:,:,:)  * (PRCS(:,:,:) * ZW1(:,:,:)              &
                                              + PRIS(:,:,:) * ZW2(:,:,:))
!
!*       4.6    compute Delta 2
!
      ZW7(:,:,:) = (ZW5(:,:,:) / (ZW6(:,:,:) * ZT(:,:,:))) *                                 &
                   (PRCS(:,:,:) * ZW1(:,:,:) *                                               &
        ((-2. * XBETAW + XGAMW * ZT(:,:,:)) / (XBETAW - XGAMW * ZT(:,:,:)) +                 &
               (XBETAW - XGAMW * ZT(:,:,:)) * (1.0 + 2.0 * ZW3(:,:,:) / ZEPS) / ZT(:,:,:)) + &
                    PRIS(:,:,:) * ZW2(:,:,:) *                                               &
        ((-2. * XBETAI + XGAMI * ZT(:,:,:)) / (XBETAI - XGAMI * ZT(:,:,:)) +                 &
               (XBETAI - XGAMI * ZT(:,:,:)) * (1.0 + 2.0 * ZW4(:,:,:) / ZEPS) / ZT(:,:,:)))
!
!*       4.7    compute Delta 1
!
      ZW6(:,:,:) = ZW5(:,:,:) * (PRCS(:,:,:) * ZW3(:,:,:) + PRIS(:,:,:) * ZW4(:,:,:) - &
                                 PRVS(:,:,:) * PTSTEP * (PRCS(:,:,:) + PRIS(:,:,:))) / &
                                 ZW6(:,:,:)
!
!*       4.8    compute the sources
!
      ZW3(:,:,:) = (ZCPH(:,:,:) / ZW5(:,:,:)) *  &
                   (-ZW6(:,:,:) * (1.0 + 0.5 * ZW6(:,:,:) * ZW7(:,:,:))) / PTSTEP
      ZW1(:,:,:) = ZW3(:,:,:) * ZCND(:,:,:)                               ! RVCNDC
      ZW2(:,:,:) = ZW3(:,:,:) * (1.0 - ZCND(:,:,:))                       ! RVDEPI
!
    ELSEWHERE
!
!*       4.9    special case when both r_c and r_i are zero
!
      ZW6(:,:,:) = ZCPH(:,:,:) + ZLV(:,:,:) * ZW1(:,:,:)               ! F'(T)
      ZW7(:,:,:) = (ZLV(:,:,:) / (ZW6(:,:,:) * ZT(:,:,:))) *                & ! Delta 2
                   (ZW1(:,:,:) *                                            &
       ((-2. * XBETAW + XGAMW * ZT(:,:,:)) / (XBETAW - XGAMW * ZT(:,:,:)) + &
              (XBETAW - XGAMW * ZT(:,:,:)) * (1.0 + 2.0 * ZW3(:,:,:) / ZEPS) / ZT(:,:,:)))
      ZW6(:,:,:) = ZLV(:,:,:) * (ZW3(:,:,:) - PRVS(:,:,:) * PTSTEP) / ZW6(:,:,:) ! Delta 1
      ZW1(:,:,:) = (ZCPH(:,:,:) / ZLV(:,:,:)) *                          & ! RVCNDC
                   (-ZW6(:,:,:) * ( 1.0 + 0.5 * ZW6(:,:,:) * ZW7(:,:,:))) / PTSTEP
      ZW2(:,:,:) = 0.0                                                     ! RVDEPI
!
    END WHERE
  END IF 
!
!*       5.     COMPUTE THE SOURCES AND STORES THE CLOUD FRACTION
!               -------------------------------------------------
!
!*       5.1    compute the sources 
!
!*       5.1.1  microphysics
!
  WHERE (ZW1(:,:,:) < 0.0)
    ZW1(:,:,:) = MAX (ZW1(:,:,:), -PRCS(:,:,:))        ! Evaporation rate
  ELSEWHERE
    ZW1(:,:,:) = MIN (ZW1(:,:,:),  PRVS(:,:,:))        ! Condensation rate
  END WHERE
!
  PRVS(:,:,:) = PRVS(:,:,:) - ZW1(:,:,:)
  PRCS(:,:,:) = PRCS(:,:,:) + ZW1(:,:,:)  
  PTHS(:,:,:) = PTHS(:,:,:) +        &
                 ZW1(:,:,:) * ZLV(:,:,:) / (ZCPH(:,:,:) * PEXNREF(:,:,:))
!
  WHERE (ZW2(:,:,:) < 0.0)
    ZW2(:,:,:) = MAX (ZW2(:,:,:), -PRIS(:,:,:))        ! Sublimation rate
  ELSEWHERE
    ZW2(:,:,:) = MIN (ZW2(:,:,:),  PRVS(:,:,:))        ! Deposition rate
  END WHERE
!
  PRVS(:,:,:) = PRVS(:,:,:) - ZW2(:,:,:)
  PRIS(:,:,:) = PRIS(:,:,:) + ZW2(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) +        &
                 ZW2(:,:,:) * ZLS(:,:,:) / (ZCPH(:,:,:) * PEXNREF(:,:,:))
!
!*       5.1.2  electricity
!
  ZWE1(:,:,:) = 0.
  ZWE2(:,:,:) = 0.
!
! the electrical process due to condensation is removed and
! capture of ions by cloud droplets is done in ion_attach_elec routine
!
! evaporation
  WHERE (ABS(PRCT(:,:,:)) > XRTMIN_ELEC(2) .AND. &
         ABS(PQCT(:,:,:)) > XQTMIN(2)      .AND. &
              ZW1(:,:,:)  < -XRTMIN(1))
    ZWE1(:,:,:) = (XFC / 3.) * (PQCT(:,:,:) / PRCT(:,:,:)) * (-ZW1(:,:,:))
    ZION_NUMBER(:,:,:) = ABS(ZWE1(:,:,:)) / XECHARGE
    ZADD(:,:,:) = 0.5 + SIGN(0.5, ZWE1(:,:,:))
    PQPIS(:,:,:) = PQPIS(:,:,:) +  &
                            ZADD(:,:,:) *ZION_NUMBER(:,:,:)
    PQNIS(:,:,:) = PQNIS(:,:,:) +  &
                            (1.-ZADD(:,:,:)) *ZION_NUMBER(:,:,:)
    PQCS(:,:,:) = PQCS(:,:,:) - ZWE1(:,:,:)
  END WHERE
!
! the electrical process due to deposition is removed and
! capture of ions by raindropsp is done in ion_attach_elec routine
!
! sublimation
  WHERE (ABS(PRIT(:,:,:)) > XRTMIN_ELEC(4) .AND. &
         ABS(PQIT(:,:,:)) > XQTMIN(4)      .AND. &
              ZW2(:,:,:)  < -XRTMIN(1))
    ZWE2(:,:,:) = (XFI / XBI) * (PQIT(:,:,:) / PRIT(:,:,:)) * (-ZW2(:,:,:))
    ZION_NUMBER(:,:,:) = ABS(ZWE2(:,:,:)) / XECHARGE
    ZADD(:,:,:) = 0.5 + SIGN(0.5, ZWE2(:,:,:))
    PQPIS(:,:,:) = PQPIS(:,:,:) +  &
                            ZADD(:,:,:) *ZION_NUMBER(:,:,:)
    PQNIS(:,:,:) = PQNIS(:,:,:) +  &
                            (1.-ZADD(:,:,:)) *ZION_NUMBER(:,:,:)
    PQIS(:,:,:) = PQIS(:,:,:) - ZWE2(:,:,:)
  END WHERE
END DO          ! end of the iterative loop
!
!
!*       5.2    compute the cloud fraction PCLDFR
!
IF (.NOT. OSUBG_COND) THEN
  WHERE (PRCS(:,:,:) + PRIS(:,:,:) > 1.E-12 / PTSTEP)
    PCLDFR(:,:,:) = 1.
  ELSEWHERE
    PCLDFR(:,:,:) = 0. 
  ENDWHERE 
  IF (SIZE(PSRCS,3) /= 0) THEN
    PSRCS(:,:,:) = PCLDFR(:,:,:) 
  END IF
ELSE
  IF (HSCONV == 'EDKF' .AND. HMF_CLOUD == 'DIRE') THEN
    PCLDFR(:,:,:) = MIN(1.,PCLDFR(:,:,:)+PCF_MF(:,:,:))
    PRCS(:,:,:)   = PRCS(:,:,:) + PRC_MF(:,:,:) / PTSTEP
    PRIS(:,:,:)   = PRIS(:,:,:)+PRI_MF(:,:,:)/PTSTEP
    PRVS(:,:,:)   = PRVS(:,:,:)- ( PRC_MF(:,:,:) + PRI_MF(:,:,:)) /PTSTEP
    PTHS(:,:,:)   = PTHS(:,:,:) +  ( PRC_MF(:,:,:) * ZLV(:,:,:) +    &
                    PRI_MF(:,:,:) * ZLS(:,:,:) ) / ZCPH(:,:,:) /     &
                    PEXNREF(:,:,:) / PTSTEP
  END IF
ENDIF
!
!
!
!*       6.  STORE THE BUDGET TERMS
!            ----------------------
!
if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'DEPI', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'DEPI', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'DEPI', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'DEPI', pris(:, :, :) * prhodj(:, :, :) )
if ( lbudget_sv ) then
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg     ), 'DEPI', pqpis(:, :, :) * prhodj(:, :, :) )
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecend     ), 'DEPI', pqnis(:, :, :) * prhodj(:, :, :) )
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1 ), 'DEPI', pqcs (:, :, :) * prhodj(:, :, :) )
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3 ), 'DEPI', pqis (:, :, :) * prhodj(:, :, :) )
end if
!------------------------------------------------------------------------------
!
END SUBROUTINE ICE_ADJUST_ELEC 
