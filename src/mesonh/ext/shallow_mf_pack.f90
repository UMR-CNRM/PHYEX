!MNH_LIC Copyright 2010-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    ######################
     MODULE MODI_SHALLOW_MF_PACK
!    ######################
!
INTERFACE
!     #################################################################
      SUBROUTINE SHALLOW_MF_PACK(KRR,KRRL,KRRI,                       &
                TPFILE,PTIME_LES,                                     &
                PTSTEP,                                               &
                PDZZ, PZZ, PDX,PDY,                                   &
                PRHODJ, PRHODREF,                                     &
                PPABSM, PEXN,                                         &
                PSFTH,PSFRV,                                          &
                PTHM,PRM,PUM,PVM,PTKEM,PSVM,                          &
                PRTHS,PRRS,PRUS,PRVS,PRSVS,                           &
                PSIGMF,PRC_MF,PRI_MF,PCF_MF,PFLXZTHVMF  )
!     #################################################################
!!
use MODD_IO,        only: TFILEDATA
use modd_precision, only: MNHTIME
!
!*               1.1  Declaration of Arguments
!
!
INTEGER,                INTENT(IN)   :: KRR        ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL       ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI       ! number of ice water var.
TYPE(TFILEDATA),        INTENT(IN)   :: TPFILE     ! Output file
REAL(kind=MNHTIME),DIMENSION(2), INTENT(OUT)  :: PTIME_LES  ! time spent in LES computations
REAL,                   INTENT(IN)   :: PTSTEP     ! Dynamical timestep 

REAL, DIMENSION(:,:,:), INTENT(IN) ::  PZZ         ! Height of flux point
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PDZZ        ! Metric coefficients
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PRHODJ      ! dry density * Grid size
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PRHODREF    ! dry density of the
                                                     ! reference state
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PPABSM      ! Pressure at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PEXN        ! Exner function at t-dt

REAL, DIMENSION(:,:),   INTENT(IN) ::  PSFTH,PSFRV ! normal surface fluxes of theta and Rv 
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PTHM        ! Theta at t-dt
REAL, DIMENSION(:,:,:,:),INTENT(IN)::  PRM         ! water var. at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PUM,PVM     ! wind components at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PTKEM       ! tke at t-dt

REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM      ! scalar variable a t-dt

REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRUS,PRVS,PRTHS ! Meso-NH sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRRS 
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS            ! Scalar sources 
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PSIGMF,PRC_MF,PRI_MF,PCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PFLXZTHVMF           ! Thermal production for TKE scheme
!
REAL, INTENT(IN) :: PDX,PDY ! Size of mesh in X/Y directions
END SUBROUTINE SHALLOW_MF_PACK

END INTERFACE
!
END MODULE MODI_SHALLOW_MF_PACK

!     #################################################################
      SUBROUTINE SHALLOW_MF_PACK(KRR,KRRL,KRRI,                       &
                TPFILE,PTIME_LES,                                     &
                PTSTEP,                                               &
                PDZZ, PZZ, PDX,PDY,                                   &
                PRHODJ, PRHODREF,                                     &
                PPABSM, PEXN,                                         &
                PSFTH,PSFRV,                                          &
                PTHM,PRM,PUM,PVM,PTKEM,PSVM,                          &
                PRTHS,PRRS,PRUS,PRVS,PRSVS,                           &
                PSIGMF,PRC_MF,PRI_MF,PCF_MF,PFLXZTHVMF  )
!     #################################################################
!!
!!****  *SHALLOW_MF_PACK* - 
!!       
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is
!!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     REFERENCE
!!     ---------
!!
!!
!!     AUTHOR
!!     ------
!!      V.Masson 09/2010
! --------------------------------------------------------------------------
! Modifications:
!  R. Honnert     07/2012: introduction of vertical wind for the height of the thermal
!  M. Leriche     02/2017: avoid negative values for sv tendencies
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  S. Riette      11/2016: support for CFRAC_ICE_SHALLOW_MF
!  P. Wautelet 28/03/2019: use MNHTIME for time measurement variables
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST, ONLY: CST
USE MODD_NEB_n, ONLY: NEBN
USE MODD_TURB_n, ONLY: TURBN
USE MODD_CTURB,  ONLY: CSTURB
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALLN, LMF_FLX
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
!
USE MODD_BUDGET,          ONLY: TBUDGETS,TBUCONF,lbudget_th,nbudget_th
USE MODD_CONF
USE MODD_IO,              ONLY: TFILEDATA
use modd_field,           ONLY: tfieldmetadata, TYPEREAL
USE MODD_NSV,             ONLY: XSVMIN, NSV_LGBEG, NSV_LGEND
USE MODD_PARAMETERS
USE MODD_PARAM_MFSHALL_n
USE modd_precision,       ONLY: MNHTIME

USE mode_budget,          ONLY: Budget_store_init, Budget_store_end, Budget_store_add
USE MODE_IO_FIELD_WRITE,  ONLY: IO_Field_write

USE MODI_DIAGNOS_LES_MF
USE MODI_SHALLOW_MF
USE MODI_SHUMAN
!
IMPLICIT NONE

!*                    0.1  Declaration of Arguments
!
!
!
INTEGER,                INTENT(IN)   :: KRR        ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL       ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI       ! number of ice water var.
TYPE(TFILEDATA),        INTENT(IN)   :: TPFILE     ! Output file
REAL(kind=MNHTIME),DIMENSION(2), INTENT(OUT)  :: PTIME_LES  ! time spent in LES computations
REAL,                   INTENT(IN)   :: PTSTEP     ! Dynamical timestep 

REAL, DIMENSION(:,:,:), INTENT(IN) ::  PZZ         ! Height of flux point
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PDZZ        ! Metric coefficients
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PRHODJ      ! dry density * Grid size
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PRHODREF    ! dry density of the
                                                     ! reference state
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PPABSM      ! Pressure at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PEXN        ! Exner function at t-dt

REAL, DIMENSION(:,:),   INTENT(IN) ::  PSFTH,PSFRV ! normal surface fluxes of theta and Rv 
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PTHM        ! Theta at t-dt
REAL, DIMENSION(:,:,:,:),INTENT(IN)::  PRM         ! water var. at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PUM,PVM     ! wind components at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PTKEM       ! tke at t-dt

REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM        ! scalar variable a t-dt

REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRUS,PRVS,PRTHS ! Meso-NH sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRRS 
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS            ! Scalar sources 
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PSIGMF,PRC_MF,PRI_MF,PCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PFLXZTHVMF           ! Thermal production for TKE scheme
!
REAL, INTENT(IN) :: PDX,PDY ! Size of mesh in X/Y directions
!
!                     0.2  Declaration of local variables
!
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDUDT_TURB   ! tendency of U   by turbulence only
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDVDT_TURB   ! tendency of V   by turbulence only
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDTHLDT_TURB ! tendency of thl by turbulence only
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDRTDT_TURB  ! tendency of rt  by turbulence only
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3),SIZE(PSVM,4)) ::  ZDSVDT_TURB  ! tendency of Sv  by turbulence only
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDUDT_MF   ! tendency of U   by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDVDT_MF   ! tendency of V   by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDTHLDT_MF ! tendency of thl by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDRTDT_MF  ! tendency of Rt by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3),SIZE(PSVM,4)) ::  ZDSVDT_MF  ! tendency of Sv by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZSIGMF,ZRC_MF,ZRI_MF,ZCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFLXZTHVMF           ! Thermal production for TKE scheme
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFLXZTHMF
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFLXZRMF
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFLXZUMF
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFLXZVMF
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZTHL_UP   ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZRT_UP    ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZRV_UP    ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZU_UP     ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZV_UP     ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZRC_UP    ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZRI_UP    ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZTHV_UP   ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZW_UP     ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFRAC_UP  ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZEMF      ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDETR     ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZENTR     ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZUMM     ! wind on mass point
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZVMM     ! wind on mass point
!
INTEGER,DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2))     :: IKLCL,IKETL,IKCTL ! level of LCL,ETL and CTL
INTEGER :: IIU, IJU, IKU, IKB, IKE, IRR, ISV  
INTEGER :: JK,JRR,JSV                          ! Loop counters

LOGICAL :: LSTATNW  !  switch for HARMONIE-AROME turb physics option
                    ! TODO: linked with modd_turbn + init at default_desfmn 

TYPE(TFIELDMETADATA) :: TZFIELD
TYPE(DIMPHYEX_t)     :: YLDIMPHYEXPACK
!------------------------------------------------------------------------
!
!!! 1. Initialisation
CALL FILL_DIMPHYEX(YLDIMPHYEXPACK, SIZE(PZZ,1), SIZE(PZZ,2), SIZE(PZZ,3))
!
! Internal Domain
IIU=SIZE(PTHM,1)
IJU=SIZE(PTHM,2)
IKU=SIZE(PTHM,3)
IKB=1+JPVEXT
IKE=IKU-JPVEXT
!
! number of moist  var
IRR=SIZE(PRM,4)
! number of scalar var
ISV=SIZE(PSVM,4)
!
! wind on mass points
ZUMM=MXF(PUM)
ZVMM=MYF(PVM)
!
!!! 2. Call of the physical parameterization of massflux vertical transport
!
LSTATNW = .FALSE.
!
CALL SHALLOW_MF(YLDIMPHYEXPACK, CST, NEBN, PARAM_MFSHALLN, TURBN, CSTURB,&
                KRR,KRRL,KRRI,ISV,                                    &
                LNOMIXLG,NSV_LGBEG,NSV_LGEND,                         &
                PTSTEP,                                               &
                PDZZ, PZZ,                                            &
                PRHODJ,PRHODREF,                                      &
                PPABSM, PEXN,                                         &
                PSFTH,PSFRV,                                          &
                PTHM,PRM,ZUMM,ZVMM,PTKEM,PSVM,                        &
                ZDUDT_MF,ZDVDT_MF,                                    &
                ZDTHLDT_MF,ZDRTDT_MF,ZDSVDT_MF,                       &
                ZSIGMF,ZRC_MF,ZRI_MF,ZCF_MF,ZFLXZTHVMF,               &
                ZFLXZTHMF,ZFLXZRMF,ZFLXZUMF,ZFLXZVMF,                 &
                ZTHL_UP,ZRT_UP,ZRV_UP,ZRC_UP,ZRI_UP,                  &
                ZU_UP, ZV_UP, ZTHV_UP, ZW_UP,                         &
                ZFRAC_UP,ZEMF,ZDETR,ZENTR,                            &
                IKLCL,IKETL,IKCTL,PDX,PDY,PRSVS,XSVMIN,               &
                TBUCONF, TBUDGETS,SIZE(TBUDGETS)                      )
!
! Fill non-declared-explicit-dimensions output variables
PSIGMF(:,:,:) = ZSIGMF(:,:,:)
PRC_MF(:,:,:) = ZRC_MF(:,:,:)
PRI_MF(:,:,:) = ZRI_MF(:,:,:)
PCF_MF(:,:,:) = ZCF_MF(:,:,:)
PFLXZTHVMF(:,:,:) = ZFLXZTHVMF(:,:,:)
!
!!! 3. Compute source terms for Meso-NH pronostic variables
!!!    ----------------------------------------------------
!
! As the pronostic variable of Meso-Nh are not (yet) the conservative variables
! the thl tendency is put in th and the rt tendency in rv
! the adjustment will do later the repartition between vapor and cloud
PRTHS(:,:,:)  = PRTHS(:,:,:)  +   &
                  PRHODJ(:,:,:)*ZDTHLDT_MF(:,:,:)
PRRS(:,:,:,1) = PRRS(:,:,:,1) +   &
                  PRHODJ(:,:,:)*ZDRTDT_MF(:,:,:)
PRUS(:,:,:)   = PRUS(:,:,:)  +MXM(  &
                  PRHODJ(:,:,:)*ZDUDT_MF(:,:,:))
PRVS(:,:,:)   = PRVS(:,:,:)  +MYM(  &
                  PRHODJ(:,:,:)*ZDVDT_MF(:,:,:))
!
DO JSV=1,ISV 
  IF (LNOMIXLG .AND. JSV >= NSV_LGBEG .AND. JSV<= NSV_LGEND) CYCLE
  PRSVS(:,:,:,JSV)   = MAX((PRSVS(:,:,:,JSV)  +    &
                  PRHODJ(:,:,:)*ZDSVDT_MF(:,:,:,JSV)),XSVMIN(JSV))
END DO     
!
!!! 4. Prints the fluxes in output file
!
IF ( LMF_FLX .AND. tpfile%lopened ) THEN
  ! stores the conservative potential temperature vertical flux
  TZFIELD = TFIELDMETADATA(          &
    CMNHNAME   = 'MF_THW_FLX',       &
    CSTDNAME   = '',                 &
    CLONGNAME  = 'MF_THW_FLX',       &
    CUNITS     = 'K m s-1',          &
    CDIR       = 'XY',               &
    CCOMMENT   = 'X_Y_Z_MF_THW_FLX', &
    NGRID      = 4,                  &
    NTYPE      = TYPEREAL,           &
    NDIMS      = 3,                  &
    LTIMEDEP   = .TRUE.              )
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZTHMF)
  !
  ! stores the conservative mixing ratio vertical flux
    TZFIELD = TFIELDMETADATA(             &
      CMNHNAME   = 'MF_RCONSW_FLX',       &
      CSTDNAME   = '',                    &
      CLONGNAME  = 'MF_RCONSW_FLX',       &
      CUNITS     = 'K m s-1',             &
      CDIR       = 'XY',                  &
      CCOMMENT   = 'X_Y_Z_MF_RCONSW_FLX', &
      NGRID      = 4,                     &
      NTYPE      = TYPEREAL,              &
      NDIMS      = 3,                     &
      LTIMEDEP   = .TRUE.                 )
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZRMF)
  !
  ! stores the theta_v vertical flux
  TZFIELD = TFIELDMETADATA(           &
    CMNHNAME   = 'MF_THVW_FLX',       &
    CSTDNAME   = '',                  &
    CLONGNAME  = 'MF_THVW_FLX',       &
    CUNITS     = 'K m s-1',           &
    CDIR       = 'XY',                &
    CCOMMENT   = 'X_Y_Z_MF_THVW_FLX', &
    NGRID      = 4,                   &
    NTYPE      = TYPEREAL,            &
    NDIMS      = 3,                   &
    LTIMEDEP   = .TRUE.               )
  CALL IO_Field_write(TPFILE,TZFIELD,PFLXZTHVMF)
  !
 IF (PARAM_MFSHALLN%LMIXUV) THEN
  ! stores the U momentum vertical flux
  TZFIELD = TFIELDMETADATA(         &
    CMNHNAME   = 'MF_UW_FLX',       &
    CSTDNAME   = '',                &
    CLONGNAME  = 'MF_UW_FLX',       &
    CUNITS     = 'm2 s-2',          &
    CDIR       = 'XY',              &
    CCOMMENT   = 'X_Y_Z_MF_UW_FLX', &
    NGRID      = 4,                 &
    NTYPE      = TYPEREAL,          &
    NDIMS      = 3,                 &
    LTIMEDEP   = .TRUE.             )
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZUMF)
  !
  ! stores the V momentum vertical flux
  TZFIELD = TFIELDMETADATA(         &
    CMNHNAME   = 'MF_VW_FLX',       &
    CSTDNAME   = '',                &
    CLONGNAME  = 'MF_VW_FLX',       &
    CUNITS     = 'm2 s-2',          &
    CDIR       = 'XY',              &
    CCOMMENT   = 'X_Y_Z_MF_VW_FLX', &
    NGRID      = 4,                 &
    NTYPE      = TYPEREAL,          &
    NDIMS      = 3,                 &
    LTIMEDEP   = .TRUE.             )
  CALL IO_Field_write(TPFILE,TZFIELD,ZFLXZVMF)
  !
 END IF
END IF
!
!!! 5. Externalised LES Diagnostic for Mass Flux Scheme
!!!    ------------------------------------------------
!
      CALL DIAGNOS_LES_MF(IIU,IJU,IKU,PTIME_LES,               &
                          ZTHL_UP,ZRT_UP,ZRV_UP,ZRC_UP,ZRI_UP, &
                          ZU_UP,ZV_UP,ZTHV_UP,ZW_UP,           &
                          ZFRAC_UP,ZEMF,ZDETR,ZENTR,           &
                          ZFLXZTHMF,ZFLXZTHVMF,ZFLXZRMF,       &
                          ZFLXZUMF,ZFLXZVMF,                   &
                          IKLCL,IKETL,IKCTL )
!               
END SUBROUTINE SHALLOW_MF_PACK
