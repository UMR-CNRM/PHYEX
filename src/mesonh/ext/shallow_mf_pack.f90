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
                HMF_UPDRAFT, HMF_CLOUD, OMIXUV,                       &
                OMF_FLX,TPFILE,PTIME_LES,                             &
                PIMPL_MF, PTSTEP,                                     &
                PDZZ, PZZ,                                            &
                PRHODJ, PRHODREF,                                     &
                PPABSM, PEXN,                                         &
                PSFTH,PSFRV,                                          &
                PTHM,PRM,PUM,PVM,PWM,PTKEM,PSVM,                      &
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
CHARACTER (LEN=4),      INTENT(IN)   :: HMF_UPDRAFT! Type of Mass Flux Scheme
                                     ! 'NONE' if no parameterization 
CHARACTER (LEN=4),      INTENT(IN)   :: HMF_CLOUD  ! Type of statistical cloud
                                                   ! scheme
LOGICAL,                INTENT(IN)   :: OMIXUV     ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: OMF_FLX    ! switch to write the
                                                   ! MF fluxes in the synchronous FM-file
TYPE(TFILEDATA),        INTENT(IN)   :: TPFILE     ! Output file
REAL(kind=MNHTIME),DIMENSION(2), INTENT(OUT)  :: PTIME_LES  ! time spent in LES computations
REAL,                   INTENT(IN)   :: PIMPL_MF   ! degre of implicitness
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
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PUM,PVM,PWM ! wind components at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PTKEM       ! tke at t-dt

REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM      ! scalar variable a t-dt

REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRUS,PRVS,PRTHS ! Meso-NH sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRRS 
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS            ! Scalar sources 
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PSIGMF,PRC_MF,PRI_MF,PCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PFLXZTHVMF           ! Thermal production for TKE scheme


END SUBROUTINE SHALLOW_MF_PACK

END INTERFACE
!
END MODULE MODI_SHALLOW_MF_PACK

!     #################################################################
      SUBROUTINE SHALLOW_MF_PACK(KRR,KRRL,KRRI,                       &
                HMF_UPDRAFT, HMF_CLOUD, OMIXUV,                       &
                OMF_FLX,TPFILE,PTIME_LES,                             &
                PIMPL_MF, PTSTEP,                                     &
                PDZZ, PZZ,                                            &
                PRHODJ, PRHODREF,                                     &
                PPABSM, PEXN,                                         &
                PSFTH,PSFRV,                                          &
                PTHM,PRM,PUM,PVM,PWM,PTKEM,PSVM,                      &
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
use modd_budget,          only: lbudget_u, lbudget_v, lbudget_th, lbudget_rv, lbudget_sv,  &
                                NBUDGET_U, NBUDGET_V, NBUDGET_TH, NBUDGET_RV, NBUDGET_SV1, &
                                tbudgets
USE MODD_CONF
USE MODD_CST
USE MODD_IO,              ONLY: TFILEDATA
use modd_field,           only: tfielddata, TYPEREAL
USE MODD_NSV
USE MODD_PARAMETERS
USE MODD_PARAM_ICE,       ONLY: CFRAC_ICE_SHALLOW_MF
USE MODD_PARAM_MFSHALL_n
use modd_precision,       only: MNHTIME

use mode_budget,          only: Budget_store_init, Budget_store_end
USE MODE_IO_FIELD_WRITE,  only: IO_Field_write

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
CHARACTER (LEN=4),      INTENT(IN)   :: HMF_UPDRAFT! Type of Mass Flux Scheme
                                     ! 'NONE' if no parameterization 
CHARACTER (LEN=4),      INTENT(IN)   :: HMF_CLOUD  ! Type of statistical cloud
                                                   ! scheme
LOGICAL,                INTENT(IN)   :: OMIXUV     ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: OMF_FLX    ! switch to write the
                                                   ! MF fluxes in the synchronous FM-file
TYPE(TFILEDATA),        INTENT(IN)   :: TPFILE     ! Output file
REAL(kind=MNHTIME),DIMENSION(2), INTENT(OUT)  :: PTIME_LES  ! time spent in LES computations
REAL,                   INTENT(IN)   :: PIMPL_MF   ! degre of implicitness
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
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PUM,PVM,PWM ! wind components at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PTKEM       ! tke at t-dt

REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM        ! scalar variable a t-dt

REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRUS,PRVS,PRTHS ! Meso-NH sources
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRRS 
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS            ! Scalar sources 
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PSIGMF,PRC_MF,PRI_MF,PCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PFLXZTHVMF           ! Thermal production for TKE scheme
!
!                     0.2  Declaration of local variables
!
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZZZ         ! Height of flux point
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDZZ        ! Metric coefficients
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZRHODJ      ! dry density * Grid size
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZRHODREF    ! dry density of the
                                                     ! reference state
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZPABSM      ! Pressure at time t-1
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZEXN        ! Exner function at t-dt

REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZTHM        ! Theta at t-dt
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3),SIZE(PRM,4)) ::  ZRM         ! water var. at t-dt
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZUM,ZVM,ZWM ! wind components at t-dt
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZTKEM       ! tke at t-dt

REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3),SIZE(PSVM,4)) ::  ZSVM        ! scalar variable a t-dt

REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDUDT_TURB   ! tendency of U   by turbulence only
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDVDT_TURB   ! tendency of V   by turbulence only
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDTHLDT_TURB ! tendency of thl by turbulence only
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDRTDT_TURB  ! tendency of rt  by turbulence only
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3),SIZE(PSVM,4)) ::  ZDSVDT_TURB  ! tendency of Sv  by turbulence only
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDUDT_MF   ! tendency of U   by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDVDT_MF   ! tendency of V   by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDTHLDT_MF ! tendency of thl by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDRTDT_MF  ! tendency of Rt by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3),SIZE(PSVM,4)) ::  ZDSVDT_MF  ! tendency of Sv by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZSIGMF,ZRC_MF,ZRI_MF,ZCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFLXZTHVMF           ! Thermal production for TKE scheme
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFLXZTHMF
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFLXZRMF
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFLXZUMF
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFLXZVMF
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZTHL_UP   ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZRT_UP    ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZRV_UP    ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZU_UP     ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZV_UP     ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZRC_UP    ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZRI_UP    ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZTHV_UP   ! updraft characteristics

REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZTHL_DO   ! downdraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZTHV_DO   ! downdraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZRT_DO    ! downdraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZU_DO     ! downdraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZV_DO     ! downdraft characteristics

REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZW_UP     ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZFRAC_UP  ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZEMF      ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDETR     ! updraft characteristics
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZENTR     ! updraft characteristics
INTEGER,DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2))     :: IKLCL,IKETL,IKCTL ! level of LCL,ETL and CTL
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2)) ::  ZSFTH    ! Surface sensible heat flux
REAL, DIMENSION(SIZE(PTHM,1)*SIZE(PTHM,2)) ::  ZSFRV    ! Surface latent   heat flux
!
!
!* 3D arrays
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZWORK    ! work array
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZUMM     ! wind on mass point
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZVMM     ! wind on mass point
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZWMM     ! wind on mass point
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDUDT    ! tendency of U   by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDVDT    ! tendency of V   by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDTHLDT  ! tendency of thl by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) ::  ZDRTDT   ! tendency of Rt by massflux scheme
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3),SIZE(PSVM,4)) ::  ZDSVDT  ! tendency of Sv by massflux scheme

INTEGER :: IIU, IJU, IKU, IKB, IKE, IRR, ISV  
INTEGER :: JK,JRR,JSV                          ! Loop counters

TYPE(TFIELDDATA) :: TZFIELD
!------------------------------------------------------------------------

!!! 1. Initialisation

! Internal Domain
IIU=SIZE(PTHM,1)
IJU=SIZE(PTHM,2)
IKU=SIZE(PTHM,3)
IKB=1+JPVEXT
IKE=IKU-JPVEXT

! number of moist  var
IRR=SIZE(PRM,4)
! number of scalar var
ISV=SIZE(PSVM,4)

if ( lbudget_u  ) call Budget_store_init( tbudgets(NBUDGET_U ), 'MAFL', prus (:, :, :)    )
if ( lbudget_v  ) call Budget_store_init( tbudgets(NBUDGET_V ), 'MAFL', prvs (:, :, :)    )
if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'MAFL', prths(:, :, :)    )
if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'MAFL', prrs (:, :, :, 1) )
if ( lbudget_sv ) then
  do jsv = 1, isv
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + jsv), 'MAFL', prsvs(:, :, :, jsv) )
  end do
end if

ZSVM(:,:,:) = 0.
!
!
! wind on mass points
ZUMM=MXF(PUM)
ZVMM=MYF(PVM)
ZWMM=MZF(PWM)
!
!!! 2. Pack input variables
!
DO JK=1,IKU
  ZZZ    (:,JK) = RESHAPE(PZZ    (:,:,JK),(/ IIU*IJU /) )
  ZDZZ   (:,JK) = RESHAPE(PDZZ    (:,:,JK),(/ IIU*IJU /) )  
  ZRHODJ (:,JK) = RESHAPE(PRHODJ  (:,:,JK),(/ IIU*IJU /) )  
  ZTHM   (:,JK) = RESHAPE(PTHM   (:,:,JK),(/ IIU*IJU /) )
  ZTKEM  (:,JK) = RESHAPE(PTKEM  (:,:,JK),(/ IIU*IJU /) )
  ZPABSM (:,JK) = RESHAPE(PPABSM (:,:,JK),(/ IIU*IJU /) )
  ZEXN   (:,JK) = RESHAPE(PEXN   (:,:,JK),(/ IIU*IJU /) )
  ZRHODJ (:,JK) = RESHAPE(PRHODJ (:,:,JK),(/ IIU*IJU /) )  
  ZRHODREF(:,JK) = RESHAPE(PRHODREF(:,:,JK),(/ IIU*IJU /) )  
  ZUM    (:,JK) = RESHAPE(ZUMM   (:,:,JK),(/ IIU*IJU /) )
  ZVM    (:,JK) = RESHAPE(ZVMM   (:,:,JK),(/ IIU*IJU /) )
  ZWM    (:,JK) = RESHAPE(ZWMM   (:,:,JK),(/ IIU*IJU /) )
  DO JRR=1,IRR
    ZRM   (:,JK,JRR) = RESHAPE(PRM    (:,:,JK,JRR),(/ IIU*IJU /) ) 
  END DO
  DO JSV=1,ISV 
    IF (LNOMIXLG .AND. JSV >= NSV_LGBEG .AND. JSV<= NSV_LGEND) CYCLE
    ZSVM(:,JK,JSV)   = RESHAPE(PSVM  (:,:,JK,JSV),(/ IIU*IJU /) ) 
  END DO  
END DO

ZSFTH(:)=RESHAPE(PSFTH(:,:),(/ IIU*IJU /) )
ZSFRV(:)=RESHAPE(PSFRV(:,:),(/ IIU*IJU /) )

!!! 3. Call of the physical parameterization of massflux vertical transport

CALL SHALLOW_MF(1,IKU,1,KRR,KRRL,KRRI,                              &
                HMF_UPDRAFT, HMF_CLOUD, CFRAC_ICE_SHALLOW_MF, OMIXUV,                  &
                LNOMIXLG,NSV_LGBEG,NSV_LGEND,                         &
                PIMPL_MF, PTSTEP,                                     &
                ZDZZ, ZZZ,                                            &
                ZRHODJ,ZRHODREF,                                      &
                ZPABSM, ZEXN,                                         &
                ZSFTH,ZSFRV,                                          &
                ZTHM,ZRM,ZUM,ZVM,ZWM,ZTKEM,ZSVM,                      &
                ZDUDT_MF,ZDVDT_MF,                                    &
                ZDTHLDT_MF,ZDRTDT_MF,ZDSVDT_MF,                       &
                ZSIGMF,ZRC_MF,ZRI_MF,ZCF_MF,ZFLXZTHVMF,               &
                ZFLXZTHMF,ZFLXZRMF,ZFLXZUMF,ZFLXZVMF,                 &
                ZTHL_UP,ZRT_UP,ZRV_UP,ZRC_UP,ZRI_UP,                  &
                ZU_UP, ZV_UP, ZTHV_UP, ZW_UP,                         &
                ZTHL_DO,ZTHV_DO,ZRT_DO,ZU_DO, ZV_DO,                  &
                ZFRAC_UP,ZEMF,ZDETR,ZENTR,                            &
                IKLCL,IKETL,IKCTL                                     )

!!! 4. Unpack output variables

ZDTHLDT(:,:,:)=RESHAPE(ZDTHLDT_MF(:,:),(/ IIU,IJU,IKU /) )
ZDRTDT(:,:,:)=RESHAPE(ZDRTDT_MF(:,:),(/ IIU,IJU,IKU /) )
ZDUDT(:,:,:)=RESHAPE(ZDUDT_MF(:,:),(/ IIU,IJU,IKU /) )
ZDVDT(:,:,:)=RESHAPE(ZDVDT_MF(:,:),(/ IIU,IJU,IKU /) )
PSIGMF(:,:,:)=RESHAPE(ZSIGMF(:,:),(/ IIU,IJU,IKU /) )
PRC_MF(:,:,:)=RESHAPE(ZRC_MF(:,:),(/ IIU,IJU,IKU /) )
PRI_MF(:,:,:)=RESHAPE(ZRI_MF(:,:),(/ IIU,IJU,IKU /) )
PCF_MF(:,:,:)=RESHAPE(ZCF_MF(:,:),(/ IIU,IJU,IKU /) )
PFLXZTHVMF(:,:,:)=RESHAPE(ZFLXZTHVMF(:,:),(/ IIU,IJU,IKU /) )
DO JSV=1,ISV 
  IF (LNOMIXLG .AND. JSV >= NSV_LGBEG .AND. JSV<= NSV_LGEND) CYCLE
    ZDSVDT(:,:,:,JSV)   = RESHAPE(ZDSVDT_MF(:,:,JSV),(/ IIU,IJU,IKU /) ) 
END DO  
!
!!! 5. Compute source terms for Meso-NH pronostic variables
!!!    ----------------------------------------------------


! As the pronostic variable of Meso-Nh are not (yet) the conservative variables
! the thl tendency is put in th and the rt tendency in rv
! the adjustment will do later the repartition between vapor and cloud
PRTHS(:,:,:)  = PRTHS(:,:,:)  +   &
                  PRHODJ(:,:,:)*ZDTHLDT(:,:,:)
PRRS(:,:,:,1) = PRRS(:,:,:,1) +   &
                  PRHODJ(:,:,:)*ZDRTDT(:,:,:)
PRUS(:,:,:)   = PRUS(:,:,:)  +MXM(  &
                  PRHODJ(:,:,:)*ZDUDT(:,:,:))
PRVS(:,:,:)   = PRVS(:,:,:)  +MYM(  &
                  PRHODJ(:,:,:)*ZDVDT(:,:,:))

DO JSV=1,ISV 
  IF (LNOMIXLG .AND. JSV >= NSV_LGBEG .AND. JSV<= NSV_LGEND) CYCLE
  PRSVS(:,:,:,JSV)   = MAX((PRSVS(:,:,:,JSV)  +    &
                  PRHODJ(:,:,:)*ZDSVDT(:,:,:,JSV)),XSVMIN(JSV))
END DO     

!!! 7. call to MesoNH budgets
if ( lbudget_u  ) call Budget_store_end( tbudgets(NBUDGET_U ), 'MAFL', prus (:, :, :)    )
if ( lbudget_v  ) call Budget_store_end( tbudgets(NBUDGET_V ), 'MAFL', prvs (:, :, :)    )
if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'MAFL', prths(:, :, :)    )
if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'MAFL', prrs (:, :, :, 1) )
if ( lbudget_sv ) then
  do jsv = 1, isv
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + jsv), 'MAFL', prsvs(:, :, :, jsv) )
  end do
end if

!!! 8. Prints the fluxes in output file
!
IF ( OMF_FLX .AND. tpfile%lopened ) THEN
  ! stores the conservative potential temperature vertical flux
  ZWORK(:,:,:)=RESHAPE(ZFLXZTHMF (:,:),(/ IIU,IJU,IKU /) )
  TZFIELD%CMNHNAME   = 'MF_THW_FLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'MF_THW_FLX'
  TZFIELD%CUNITS     = 'K m s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_MF_THW_FLX'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK)
  !
  ! stores the conservative mixing ratio vertical flux
  ZWORK(:,:,:)=RESHAPE(ZFLXZRMF(:,:),(/ IIU,IJU,IKU /) )
  TZFIELD%CMNHNAME   = 'MF_RCONSW_FLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'MF_RCONSW_FLX'
  TZFIELD%CUNITS     = 'K m s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_MF_RCONSW_FLX'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK)
  !
  ! stores the theta_v vertical flux
  TZFIELD%CMNHNAME   = 'MF_THVW_FLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'MF_THVW_FLX'
  TZFIELD%CUNITS     = 'K m s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_MF_THVW_FLX'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,PFLXZTHVMF)
  !
 IF (OMIXUV) THEN
  ! stores the U momentum vertical flux
  ZWORK(:,:,:)=RESHAPE(ZFLXZUMF(:,:),(/ IIU,IJU,IKU /) )
  TZFIELD%CMNHNAME   = 'MF_UW_FLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'MF_UW_FLX'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_MF_UW_FLX'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK)
  !
  ! stores the V momentum vertical flux
  ZWORK(:,:,:)=RESHAPE(ZFLXZVMF(:,:),(/ IIU,IJU,IKU /) )
  TZFIELD%CMNHNAME   = 'MF_VW_FLX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'MF_VW_FLX'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_MF_VW_FLX'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK)
  !
 END IF
END IF

!!! 9. Externalised LES Diagnostic for Mass Flux Scheme
!!!    ------------------------------------------------

      CALL DIAGNOS_LES_MF(IIU,IJU,IKU,PTIME_LES,               &
                          ZTHL_UP,ZRT_UP,ZRV_UP,ZRC_UP,ZRI_UP, &
                          ZU_UP,ZV_UP,ZTHV_UP,ZW_UP,           &
                          ZFRAC_UP,ZEMF,ZDETR,ZENTR,           &
                          ZFLXZTHMF,ZFLXZTHVMF,ZFLXZRMF,       &
                          ZFLXZUMF,ZFLXZVMF,                   &
                          IKLCL,IKETL,IKCTL )
               

END SUBROUTINE SHALLOW_MF_PACK
