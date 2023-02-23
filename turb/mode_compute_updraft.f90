!MNH_LIC Copyright 2004-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
     MODULE MODE_COMPUTE_UPDRAFT
!    ###########################
!
IMPLICIT NONE
CONTAINS
      SUBROUTINE COMPUTE_UPDRAFT(D,CST,NEB,PARAMMF,TURBN,CSTURB,  &
                                 KSV, HFRAC_ICE,                  &
                                 OENTR_DETR,                      &
                                 ONOMIXLG,KSV_LGBEG,KSV_LGEND,    &
                                 PZZ,PDZZ,                        &
                                 PSFTH,PSFRV,                     &
                                 PPABSM,PRHODREF,PUM,PVM, PTKEM,  &
                                 PTHM,PRVM,PTHLM,PRTM,            &
                                 PSVM,PTHL_UP,PRT_UP,             &
                                 PRV_UP,PRC_UP,PRI_UP,PTHV_UP,    &
                                 PW_UP,PU_UP, PV_UP, PSV_UP,      &
                                 PFRAC_UP,PFRAC_ICE_UP,PRSAT_UP,  &
                                 PEMF,PDETR,PENTR,                &
                                 PBUO_INTEG,KKLCL,KKETL,KKCTL,    &
                                 PDEPTH, PDX, PDY     )

!     #################################################################
!!
!!****  *COMPUTE_UPDRAFT* - calculates caracteristics of the updraft 
!!                         
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to build the updraft model 
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
!!      !!     REFERENCE
!!     ---------
!!       Book 1 of Meso-NH documentation (chapter Turbulence)
!!       Soares et al. 2004 QJ
!!
!!     AUTHOR
!!     ------
!!     J.Pergaud
!!     V.Masson : Optimization 07/2010
!!     S. Riette : 07/2010 : modification for reproducibility  
!!     S. Riette may 2011: ice added, interface modified
!!     S. Riette Jan 2012: support for both order of vertical levels
!!     V.Masson, C.Lac : 02/2011 : SV_UP initialized by a non-zero value
!!     S. Riette Apr 2013: improvement of continuity at the condensation level
!!     R.Honnert Oct 2016 : Add ZSURF and Update with AROME
!!     Q.Rodier  01/2019 : support RM17 mixing length
!!     R.Honnert 01/2019 : add LGZ (reduction of the mass-flux surface closure with the resolution)
!!     S. Riette 06/2022: compute_entr_detr is inlined
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST,             ONLY: CST_t
USE MODD_NEB,             ONLY: NEB_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
USE MODD_TURB_n,          ONLY: TURB_t
USE MODD_CTURB,           ONLY: CSTURB_t
!
USE MODI_SHUMAN_MF, ONLY: MZM_MF, MZF_MF, GZ_M_W_MF

USE MODE_COMPUTE_BL89_ML, ONLY: COMPUTE_BL89_ML
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

IMPLICIT NONE

!*                    1.1  Declaration of Arguments
!
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(NEB_t),            INTENT(IN)   :: NEB
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF
TYPE(TURB_t),           INTENT(IN)   :: TURBN
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
INTEGER,                INTENT(IN)   :: KSV
CHARACTER(LEN=1),       INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
LOGICAL,                INTENT(IN) :: OENTR_DETR! flag to recompute entrainment, detrainment and mass flux
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PZZ       !  Height at the flux point
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PDZZ      !  Metrics coefficient
 
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSFTH,PSFRV
! normal surface fluxes of theta,rv,(u,v) parallel to the orography
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PPABSM     ! Pressure at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PRHODREF   ! dry density of the
                                                  ! reference state
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PUM        ! u mean wind
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PVM        ! v mean wind
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PTKEM      ! TKE at t-dt
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PTHM           ! liquid pot. temp. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PRVM           ! vapor mixing ratio at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PTHLM,PRTM     ! cons. var. at t-dt

REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN)   ::  PSVM           ! scalar var. at t-dt

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  ::  PTHL_UP,PRT_UP   ! updraft properties
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  ::  PU_UP, PV_UP     ! updraft wind components
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT)::  PRV_UP,PRC_UP, & ! updraft rv, rc
                                         PRI_UP,PTHV_UP,& ! updraft ri, THv
                                         PW_UP,PFRAC_UP,& ! updraft w, fraction
                                         PFRAC_ICE_UP,&   ! liquid/solid fraction in updraft
                                         PRSAT_UP         ! Rsat

REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT)  ::  PSV_UP           ! updraft scalar var. 
                                         
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT)::  PEMF,PDETR,PENTR ! Mass_flux,
                                                          ! detrainment,entrainment
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PBUO_INTEG       ! Integrated Buoyancy 
INTEGER, DIMENSION(D%NIJT),  INTENT(INOUT) :: KKLCL,KKETL,KKCTL! LCL, ETL, CTL
REAL, DIMENSION(D%NIJT),     INTENT(OUT)   :: PDEPTH           ! Deepness of cloud
REAL,                   INTENT(IN)    :: PDX, PDY
!                       1.2  Declaration of local variables
!
!
! Mean environment variables at t-dt at flux point
REAL, DIMENSION(D%NIJT,D%NKT) ::    &
                        ZTHM_F,ZRVM_F                 ! Theta,rv of
                                                      ! updraft environnement
REAL, DIMENSION(D%NIJT,D%NKT) ::    &
                        ZRTM_F, ZTHLM_F, ZTKEM_F,&    ! rt, thetal,TKE,pressure,
                        ZUM_F,ZVM_F,ZRHO_F,      &    ! density,momentum
                        ZPRES_F,ZTHVM_F,ZTHVM,   &    ! interpolated at the flux point
                        ZG_O_THVREF,             &    ! g*ThetaV ref
                        ZW_UP2,                  &    ! w**2  of the updraft
                        ZBUO_INTEG_DRY, ZBUO_INTEG_CLD,&! Integrated Buoyancy
                        ZENTR_CLD,ZDETR_CLD           ! wet entrainment and detrainment

REAL, DIMENSION(D%NIJT,D%NKT,KSV) :: &
                        ZSVM_F ! scalar variables 

                        
REAL, DIMENSION(D%NIJT,D%NKT) ::  &
                        ZTH_UP,                  &    ! updraft THETA 
                        ZRC_MIX, ZRI_MIX              ! guess of Rc and Ri for KF mixture

REAL, DIMENSION(D%NIJT,D%NKT) ::  ZCOEF  ! diminution coefficient for too high clouds 
                        
REAL, DIMENSION(D%NIJT)            ::  ZWTHVSURF  ! Surface w'thetav'

REAL  :: ZRDORV       ! RD/RV
REAL  :: ZRVORD       ! RV/RD


REAL, DIMENSION(D%NIJT) :: ZMIX1,ZMIX2,ZMIX3_CLD,ZMIX2_CLD

REAL, DIMENSION(D%NIJT) :: ZLUP         ! Upward Mixing length from the ground

INTEGER  :: JK,JIJ,JSV          ! loop counters

LOGICAL, DIMENSION(D%NIJT) ::  GTEST,GTESTLCL,GTESTETL
                               ! Test if the ascent continue, if LCL or ETL is reached
LOGICAL                          ::  GLMIX 
                               ! To choose upward or downward mixing length
LOGICAL, DIMENSION(D%NIJT)              :: GWORK1
LOGICAL, DIMENSION(D%NIJT,D%NKT) :: GWORK2

INTEGER  :: ITEST

REAL, DIMENSION(D%NIJT) :: ZRC_UP, ZRI_UP, ZRV_UP,&
                                 ZRSATW, ZRSATI,&
                                 ZPART_DRY

REAL  :: ZDEPTH_MAX1, ZDEPTH_MAX2 ! control auto-extinction process

REAL  :: ZTMAX,ZRMAX  ! control value

REAL, DIMENSION(D%NIJT) :: ZSURF
REAL, DIMENSION(D%NIJT,D%NKT) :: ZSHEAR,ZDUDZ,ZDVDZ ! vertical wind shear
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWK
REAL, DIMENSION(D%NIJT,16) :: ZBUF
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!                       1.3  Declaration of additional local variables for compute_entr_detr
!
! Variables for cloudy part
REAL, DIMENSION(D%NIJT) :: ZKIC, ZKIC_F2  ! fraction of env. mass in the muxtures
REAL, DIMENSION(D%NIJT) :: ZEPSI,ZDELTA   ! factor entrainment detrainment
REAL                   :: ZEPSI_CLOUD    ! factor entrainment detrainment
REAL                   :: ZCOEFFMF_CLOUD ! factor for compputing entr. detr.
REAL, DIMENSION(D%NIJT) :: ZMIXTHL,ZMIXRT ! Thetal and rt in the mixtures
REAL, DIMENSION(D%NIJT) :: ZTHMIX         ! Theta and Thetav  of mixtures
REAL, DIMENSION(D%NIJT) :: ZRVMIX,ZRCMIX,ZRIMIX ! mixing ratios in mixtures
REAL, DIMENSION(D%NIJT) :: ZTHVMIX, ZTHVMIX_F2 ! Theta and Thetav  of mixtures
REAL, DIMENSION(D%NIJT) :: ZTHV_UP_F2     ! thv_up at flux point kk+kkl
REAL, DIMENSION(D%NIJT) :: ZRSATW_ED, ZRSATI_ED ! working arrays (mixing ratio at saturation)
REAL, DIMENSION(D%NIJT) :: ZTHV           ! theta V of environment at the bottom of cloudy part  
REAL                   :: ZKIC_INIT      !Initial value of ZKIC
REAL                   :: ZCOTHVU              ! Variation of Thvup between bottom and top of cloudy part

! Variables for dry part
REAL                   :: ZFOESW, ZFOESI       ! saturating vapor pressure
REAL                   :: ZDRSATODP            ! d.Rsat/dP
REAL                   :: ZT                   ! Temperature
REAL                   :: ZWK0D                ! Work array

! Variables for dry and cloudy parts
REAL, DIMENSION(D%NIJT) :: ZCOEFF_MINUS_HALF,&  ! Variation of Thv between mass points kk-kkl and kk
                                  ZCOEFF_PLUS_HALF     ! Variation of Thv between mass points kk and kk+kkl
REAL, DIMENSION(D%NIJT) :: ZPRE                 ! pressure at the bottom of the cloudy part
REAL, DIMENSION(D%NIJT) :: ZG_O_THVREF_ED
REAL, DIMENSION(D%NIJT) :: ZFRAC_ICE            ! fraction of ice
REAL, DIMENSION(D%NIJT) :: ZDZ_STOP,&           ! Exact Height of the LCL above flux level KK
                          ZTHV_MINUS_HALF,&    ! Thv at flux point(kk)  
                          ZTHV_PLUS_HALF       ! Thv at flux point(kk+kkl)
REAL                   :: ZDZ                  ! Delta Z used in computations
INTEGER :: JKLIM
INTEGER :: IIJB,IIJE ! physical horizontal domain indices
INTEGER :: IKT,IKB,IKE,IKL
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAFT',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
IKB=D%NKB
IKE=D%NKE
IKL=D%NKL
!
! Thresholds for the  perturbation of
! theta_l and r_t at the first level of the updraft
ZTMAX=2.0
ZRMAX=1.E-3
!------------------------------------------------------------------------

!                     INITIALISATION

! Initialisation of the constants   
ZRDORV   = CST%XRD / CST%XRV   !=0.622
ZRVORD   = (CST%XRV / CST%XRD) 

ZDEPTH_MAX1=3000. ! clouds with depth inferior to this value are keeped untouched
ZDEPTH_MAX2=4000. ! clouds with depth superior to this value are suppressed

!                 Local variables, internal domain

IF (OENTR_DETR) THEN
  ! Initialisation of intersesting Level :LCL,ETL,CTL
  KKLCL(:)=IKE
  KKETL(:)=IKE
  KKCTL(:)=IKE

  !
  ! Initialisation
  !* udraft governing variables
  PEMF(:,:)=0.
  PDETR(:,:)=0.
  PENTR(:,:)=0.

  ! Initialisation
  !* updraft core variables
  PRV_UP(:,:)=0.
  PRC_UP(:,:)=0.
  PRI_UP(:,:)=0.
  PW_UP(:,:)=0.
  ZTH_UP(:,:)=0.
  PFRAC_UP(:,:)=0.
  PTHV_UP(:,:)=0.

  PBUO_INTEG=0.

  PFRAC_ICE_UP(:,:)=0.
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PRSAT_UP(:,:)=PRVM(:,:) ! should be initialised correctly but is (normaly) not used
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

  !cloud/dry air mixture cloud content
  ZRC_MIX = 0.
  ZRI_MIX = 0.

END IF

! Initialisation of environment variables at t-dt
! variables at flux level
CALL MZM_MF(D, PTHLM(:,:), ZTHLM_F(:,:))
CALL MZM_MF(D, PRTM(:,:), ZRTM_F (:,:))
CALL MZM_MF(D, PUM(:,:), ZUM_F  (:,:))
CALL MZM_MF(D, PVM(:,:), ZVM_F  (:,:))
CALL MZM_MF(D, PTKEM(:,:), ZTKEM_F(:,:))

DO JSV=1,KSV
  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  CALL MZM_MF(D, PSVM(:,:,JSV), ZSVM_F(:,:,JSV))
END DO
!                     
!          Initialisation of updraft characteristics 
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PTHL_UP(:,:)=ZTHLM_F(:,:)
PRT_UP(:,:)=ZRTM_F(:,:)
PU_UP(:,:)=ZUM_F(:,:)
PV_UP(:,:)=ZVM_F(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT,JSV=1:KSV)
PSV_UP(:,:,:)=ZSVM_F(:,:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT,JSV=1:KSV)

! Computation or initialisation of updraft characteristics at the KKB level
! thetal_up,rt_up,thetaV_up, w2,Buoyancy term and mass flux (PEMF)
!$mnh_expand_array(JIJ=IIJB:IIJE)
PTHL_UP(:,IKB)= ZTHLM_F(:,IKB)+ &
                            & MAX(0.,MIN(ZTMAX,(PSFTH(:)/SQRT(ZTKEM_F(:,IKB)))* PARAMMF%XALP_PERT))
PRT_UP(:,IKB) = ZRTM_F(:,IKB)+ &
                            & MAX(0.,MIN(ZRMAX,(PSFRV(:)/SQRT(ZTKEM_F(:,IKB)))* PARAMMF%XALP_PERT)) 
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

IF (OENTR_DETR) THEN
  CALL MZM_MF(D, PTHM (:,:), ZTHM_F (:,:))
  CALL MZM_MF(D, PPABSM(:,:), ZPRES_F(:,:))
  CALL MZM_MF(D, PRHODREF(:,:), ZRHO_F (:,:))
  CALL MZM_MF(D, PRVM(:,:), ZRVM_F (:,:))

  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ! thetav at mass and flux levels
  ZTHVM_F(:,:)=ZTHM_F(:,:)* &
                                    &((1.+ZRVORD*ZRVM_F(:,:))/(1.+ZRTM_F(:,:)))
  ZTHVM(:,:)=PTHM(:,:)* &
                                    &((1.+ZRVORD*PRVM(:,:))/(1.+PRTM(:,:)))

  PTHV_UP(:,:)=ZTHVM_F(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

  ZW_UP2(:,:)=0.
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZW_UP2(:,IKB) = MAX(0.0001,(2./3.)*ZTKEM_F(:,IKB))

  ! Computation of non conservative variable for the KKB level of the updraft
  ! (all or nothing ajustement)
  PRC_UP(:,IKB)=0.
  PRI_UP(:,IKB)=0.
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  CALL TH_R_FROM_THL_RT(CST, NEB, D%NIJT, HFRAC_ICE,PFRAC_ICE_UP(:,IKB),ZPRES_F(:,IKB), &
             PTHL_UP(:,IKB),PRT_UP(:,IKB),ZTH_UP(:,IKB), &
             PRV_UP(:,IKB),PRC_UP(:,IKB),PRI_UP(:,IKB),ZRSATW(:),ZRSATI(:), OOCEAN=.FALSE., &
             PBUF=ZBUF(:,:), KB=D%NIJB, KE=D%NIJE)

  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ! compute updraft thevav and buoyancy term at KKB level
  PTHV_UP(:,IKB) = ZTH_UP(:,IKB)*&
                               & ((1+ZRVORD*PRV_UP(:,IKB))/(1+PRT_UP(:,IKB)))
  ! compute mean rsat in updraft
  PRSAT_UP(:,IKB) = ZRSATW(:)*(1-PFRAC_ICE_UP(:,IKB)) + &
                              & ZRSATI(:)*PFRAC_ICE_UP(:,IKB)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  ! Closure assumption for mass flux at KKB level
  !

  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZG_O_THVREF(:,:)=CST%XG/ZTHVM_F(:,:)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

  ! compute L_up
  GLMIX=.TRUE.
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZTKEM_F(:,IKB)=0.
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  !
  IF(TURBN%CTURBLEN=='RM17') THEN
    CALL GZ_M_W_MF(D, PUM, PDZZ, ZWK)
    CALL MZF_MF(D, ZWK, ZDUDZ)
    CALL GZ_M_W_MF(D, PVM, PDZZ, ZWK)
    CALL MZF_MF(D, ZWK, ZDVDZ)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZSHEAR(:,:) = SQRT(ZDUDZ(:,:)**2 + ZDVDZ(:,:)**2)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ELSE
    ZSHEAR = 0. !no shear in bl89 mixing length
  END IF
  !
#ifdef REPRO48
  CALL COMPUTE_BL89_ML(D, CST, CSTURB, PDZZ,ZTKEM_F(:,IKB),&
                      &ZG_O_THVREF(:,IKB),ZTHVM,IKB,GLMIX,.TRUE.,ZSHEAR,ZLUP)
#else
  CALL COMPUTE_BL89_ML(D, CST, CSTURB, PDZZ,ZTKEM_F(:,IKB),&
                      &ZG_O_THVREF(:,IKB),ZTHVM,IKB,GLMIX,.FALSE.,ZSHEAR,ZLUP)
#endif
  !$mnh_expand_where(JIJ=IIJB:IIJE)
  ZLUP(:)=MAX(ZLUP(:),1.E-10)

  ! Compute Buoyancy flux at the ground
  ZWTHVSURF(:) = (ZTHVM_F(:,IKB)/ZTHM_F(:,IKB))*PSFTH(:)+     &
                (0.61*ZTHM_F(:,IKB))*PSFRV(:)

  ! Mass flux at KKB level (updraft triggered if PSFTH>0.)
  IF (PARAMMF%LGZ) THEN
    IF(PDX==0. .OR. PDY==0.) THEN                                                                                                   
      CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'COMPUTE_UPDRAFT', 'PDX or PDY is NULL with option LGZ!')                                  
    ENDIF
    ZSURF(:)=TANH(PARAMMF%XGZ*SQRT(PDX*PDY)/ZLUP(:))
  ELSE
    ZSURF(:)=1.
  END IF
  WHERE (ZWTHVSURF(:)>0.)
    PEMF(:,IKB) = PARAMMF%XCMF * ZSURF(:) * ZRHO_F(:,IKB) *  &
            ((ZG_O_THVREF(:,IKB))*ZWTHVSURF(:)*ZLUP(:))**(1./3.)
    PFRAC_UP(:,IKB)=MIN(PEMF(:,IKB)/(SQRT(ZW_UP2(:,IKB))*ZRHO_F(:,IKB)), &
                                   &PARAMMF%XFRAC_UP_MAX)
    ZW_UP2(:,IKB)=(PEMF(:,IKB)/(PFRAC_UP(:,IKB)*ZRHO_F(:,IKB)))**2
    GTEST(:)=.TRUE.
  ELSEWHERE
    PEMF(:,IKB) =0.
    GTEST(:)=.FALSE.
  ENDWHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)
ELSE
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  GTEST(:)=PEMF(:,IKB+IKL)>0.
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END IF

!--------------------------------------------------------------------------

!                        3. Vertical ascending loop
!                           -----------------------
!
! If GTEST = T the updraft starts from the KKB level and stops when GTEST becomes F
!
!
GTESTLCL(:)=.FALSE.
GTESTETL(:)=.FALSE.

!       Loop on vertical level

DO JK=IKB,IKE-IKL,IKL

  ! IF the updraft top is reached for all column, stop the loop on levels
  ITEST=COUNT(GTEST(:))
  IF (ITEST==0) CYCLE

  !       Computation of entrainment and detrainment with KF90
  !       parameterization in clouds and LR01 in subcloud layer


  ! to find the LCL (check if JK is LCL or not)
  !$mnh_expand_where(JIJ=IIJB:IIJE)
  WHERE ((PRC_UP(:,JK)+PRI_UP(:,JK)>0.).AND.(.NOT.(GTESTLCL(:))))
      KKLCL(:) = JK           
      GTESTLCL(:)=.TRUE.
  ENDWHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)

  ! COMPUTE PENTR and PDETR at mass level JK
  IF (OENTR_DETR) THEN
    IF(JK/=IKB) THEN
      !$mnh_expand_array(JIJ=IIJB:IIJE)
      ZRC_MIX(:,JK) = ZRC_MIX(:,JK-IKL) ! guess of Rc of mixture
      ZRI_MIX(:,JK) = ZRI_MIX(:,JK-IKL) ! guess of Ri of mixture
      !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    ENDIF
    CALL COMPUTE_ENTR_DETR(D, CST, NEB, PARAMMF, JK,IKB,IKE,IKL,GTEST,GTESTLCL,HFRAC_ICE,PFRAC_ICE_UP(:,JK),&
                           PRHODREF(:,JK),ZPRES_F(:,JK),ZPRES_F(:,JK+IKL),&
                           PZZ(:,:),PDZZ(:,:),ZTHVM(:,:),  &
                           PTHLM(:,:),PRTM(:,:),ZW_UP2(:,:),ZTH_UP(:,JK),   &
                           PTHL_UP(:,JK),PRT_UP(:,JK),ZLUP(:),         &
                           PRC_UP(:,JK),PRI_UP(:,JK),PTHV_UP(:,JK),&
                           PRSAT_UP(:,JK),ZRC_MIX(:,JK),ZRI_MIX(:,JK),                 &
                           PENTR(:,JK),PDETR(:,JK),ZENTR_CLD(:,JK),ZDETR_CLD(:,JK),&
                           ZBUO_INTEG_DRY(:,JK), ZBUO_INTEG_CLD(:,JK), &
                           ZPART_DRY(:)   )
    !$mnh_expand_where(JIJ=IIJB:IIJE)
    PBUO_INTEG(:,JK)=ZBUO_INTEG_DRY(:,JK)+ZBUO_INTEG_CLD(:,JK)

    IF (JK==IKB) THEN
       PDETR(:,JK)=0.
       ZDETR_CLD(:,JK)=0.
    ENDIF   
 
    !       Computation of updraft characteristics at level JK+KKL
    WHERE(GTEST(:))
      ZMIX1(:)=0.5*(PZZ(:,JK+IKL)-PZZ(:,JK))*&
                          &(PENTR(:,JK)-PDETR(:,JK))
      PEMF(:,JK+IKL)=PEMF(:,JK)*EXP(2*ZMIX1(:))
    ENDWHERE
    !$mnh_end_expand_where(JIJ=IIJB:IIJE)
  ELSE !OENTR_DETR
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    GTEST(:) = (PEMF(:,JK+IKL)>0.)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  END IF !OENTR_DETR
  
  ! stop the updraft if MF becomes negative
  !$mnh_expand_where(JIJ=IIJB:IIJE)
  WHERE (GTEST(:).AND.(PEMF(:,JK+IKL)<=0.))
    PEMF(:,JK+IKL)=0.
    KKCTL(:) = JK+IKL
    GTEST(:)=.FALSE.
    PFRAC_ICE_UP(:,JK+IKL)=PFRAC_ICE_UP(:,JK)
    PRSAT_UP(:,JK+IKL)=PRSAT_UP(:,JK)
  ENDWHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)

  ! If the updraft did not stop, compute cons updraft characteritics at jk+KKL
  DO JIJ=IIJB,IIJE
    IF(GTEST(JIJ)) THEN
      ZMIX2(JIJ) = (PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*PENTR(JIJ,JK) !&
      ZMIX3_CLD(JIJ) = (PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*(1.-ZPART_DRY(JIJ))*ZDETR_CLD(JIJ,JK) !&                   
      ZMIX2_CLD(JIJ) = (PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*(1.-ZPART_DRY(JIJ))*ZENTR_CLD(JIJ,JK)
#ifdef REPRO48                  
      PTHL_UP(JIJ,JK+IKL)=(PTHL_UP(JIJ,JK)*(1.-0.5*ZMIX2(JIJ)) + PTHLM(JIJ,JK)*ZMIX2(JIJ)) &
                            /(1.+0.5*ZMIX2(JIJ))   
      PRT_UP(JIJ,JK+IKL) =(PRT_UP (JIJ,JK)*(1.-0.5*ZMIX2(JIJ)) + PRTM(JIJ,JK)*ZMIX2(JIJ))  &
                            /(1.+0.5*ZMIX2(JIJ))
#else
      PTHL_UP(JIJ,JK+IKL)=PTHL_UP(JIJ,JK)*EXP(-ZMIX2(JIJ)) + PTHLM(JIJ,JK)*(1-EXP(-ZMIX2(JIJ)))
      PRT_UP(JIJ,JK+IKL) =PRT_UP (JIJ,JK)*EXP(-ZMIX2(JIJ)) +  PRTM(JIJ,JK)*(1-EXP(-ZMIX2(JIJ)))
#endif
    ENDIF
  ENDDO
  
  IF(PARAMMF%LMIXUV) THEN
    IF(JK/=IKB) THEN
      !$mnh_expand_where(JIJ=IIJB:IIJE)
      WHERE(GTEST(:))
        PU_UP(:,JK+IKL) = (PU_UP(:,JK)*(1-0.5*ZMIX2(:)) + &
                                        &PUM(:,JK)*ZMIX2(:)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(:,JK+IKL)-PZZ(:,JK))*&
                          ((PUM(:,JK+IKL)-PUM(:,JK))/PDZZ(:,JK+IKL)+&
                           (PUM(:,JK)-PUM(:,JK-IKL))/PDZZ(:,JK))        )   &
                          /(1+0.5*ZMIX2(:))
        PV_UP(:,JK+IKL) = (PV_UP(:,JK)*(1-0.5*ZMIX2(:)) + &
                                        &PVM(:,JK)*ZMIX2(:)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(:,JK+IKL)-PZZ(:,JK))*&
                          ((PVM(:,JK+IKL)-PVM(:,JK))/PDZZ(:,JK+IKL)+&
                           (PVM(:,JK)-PVM(:,JK-IKL))/PDZZ(:,JK))    )   &
                          /(1+0.5*ZMIX2(:))
      ENDWHERE
      !$mnh_end_expand_where(JIJ=IIJB:IIJE)
    ELSE
      !$mnh_expand_where(JIJ=IIJB:IIJE)
      WHERE(GTEST(:))
        PU_UP(:,JK+IKL) = (PU_UP(:,JK)*(1-0.5*ZMIX2(:)) + &
                                        &PUM(:,JK)*ZMIX2(:)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(:,JK+IKL)-PZZ(:,JK))*&
                          ((PUM(:,JK+IKL)-PUM(:,JK))/PDZZ(:,JK+IKL))        )   &
                          /(1+0.5*ZMIX2(:))
        PV_UP(:,JK+IKL) = (PV_UP(:,JK)*(1-0.5*ZMIX2(:)) + &
                                        &PVM(:,JK)*ZMIX2(:)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(:,JK+IKL)-PZZ(:,JK))*&
                          ((PVM(:,JK+IKL)-PVM(:,JK))/PDZZ(:,JK+IKL))    )   &
                          /(1+0.5*ZMIX2(:))
      ENDWHERE
      !$mnh_end_expand_where(JIJ=IIJB:IIJE)
    ENDIF
  ENDIF !PARAMMF%LMIXUV
  DO JSV=1,KSV 
    IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
    !$mnh_expand_where(JIJ=IIJB:IIJE)
    WHERE(GTEST(:)) 
      PSV_UP(:,JK+IKL,JSV) = (PSV_UP(:,JK,JSV)*(1-0.5*ZMIX2(:)) + &
                   PSVM(:,JK,JSV)*ZMIX2(:))  /(1+0.5*ZMIX2(:))
    ENDWHERE
    !$mnh_end_expand_where(JIJ=IIJB:IIJE)
  END DO  
  
  IF (OENTR_DETR) THEN

    ! Compute non cons. var. at level JK+KKL
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZRC_UP(:)=PRC_UP(:,JK) ! guess = level just below
    ZRI_UP(:)=PRI_UP(:,JK) ! guess = level just below
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    CALL TH_R_FROM_THL_RT(CST, NEB, D%NIJT, HFRAC_ICE,PFRAC_ICE_UP(:,JK+IKL),ZPRES_F(:,JK+IKL), &
            PTHL_UP(:,JK+IKL),PRT_UP(:,JK+IKL),ZTH_UP(:,JK+IKL),              &
            ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:), OOCEAN=.FALSE., &
            PBUF=ZBUF(:,:), KB=D%NIJB, KE=D%NIJE)
    !$mnh_expand_where(JIJ=IIJB:IIJE)
    WHERE(GTEST(:))
      PRC_UP(:,JK+IKL)=ZRC_UP(:)
      PRV_UP(:,JK+IKL)=ZRV_UP(:)
      PRI_UP(:,JK+IKL)=ZRI_UP(:)
      PRSAT_UP(:,JK+IKL) = ZRSATW(:)*(1-PFRAC_ICE_UP(:,JK+IKL)) + &
                                     & ZRSATI(:)*PFRAC_ICE_UP(:,JK+IKL)
    ENDWHERE

    ! Compute the updraft theta_v, buoyancy and w**2 for level JK+KKL
    WHERE(GTEST(:))
      PTHV_UP(:,JK+IKL) = ZTH_UP(:,JK+IKL)* &
                                    & ((1+ZRVORD*PRV_UP(:,JK+IKL))/(1+PRT_UP(:,JK+IKL)))
      WHERE (ZBUO_INTEG_DRY(:,JK)>0.)
        ZW_UP2(:,JK+IKL)  = ZW_UP2(:,JK) + 2.*(PARAMMF%XABUO-PARAMMF%XBENTR*PARAMMF%XENTR_DRY)* &
                                                                &ZBUO_INTEG_DRY(:,JK)
      ELSEWHERE
        ZW_UP2(:,JK+IKL)  = ZW_UP2(:,JK) + 2.*PARAMMF%XABUO* ZBUO_INTEG_DRY(:,JK)
      ENDWHERE
      ZW_UP2(:,JK+IKL)  = ZW_UP2(:,JK+IKL)*(1.-(PARAMMF%XBDETR*ZMIX3_CLD(:)+ &
                                                                       &PARAMMF%XBENTR*ZMIX2_CLD(:)))&
              /(1.+(PARAMMF%XBDETR*ZMIX3_CLD(:)+PARAMMF%XBENTR*ZMIX2_CLD(:))) &
              +2.*(PARAMMF%XABUO)*ZBUO_INTEG_CLD(:,JK)/ &
              &(1.+(PARAMMF%XBDETR*ZMIX3_CLD(:)+PARAMMF%XBENTR*ZMIX2_CLD(:)))
    ENDWHERE

    ! Test if the updraft has reach the ETL
    WHERE (GTEST(:).AND.(PBUO_INTEG(:,JK)<=0.))
      KKETL(:) = JK+IKL
      GTESTETL(:)=.TRUE.
    ELSEWHERE
      GTESTETL(:)=.FALSE.
    ENDWHERE

    ! Test is we have reached the top of the updraft
    WHERE (GTEST(:).AND.((ZW_UP2(:,JK+IKL)<=0.).OR.(PEMF(:,JK+IKL)<=0.)))
        ZW_UP2(:,JK+IKL)=0.
        PEMF(:,JK+IKL)=0.
        GTEST(:)=.FALSE.
        PTHL_UP(:,JK+IKL)=ZTHLM_F(:,JK+IKL)
        PRT_UP(:,JK+IKL)=ZRTM_F(:,JK+IKL)
        PRC_UP(:,JK+IKL)=0.
        PRI_UP(:,JK+IKL)=0.
        PRV_UP(:,JK+IKL)=0.
        PTHV_UP(:,JK+IKL)=ZTHVM_F(:,JK+IKL)
        PFRAC_UP(:,JK+IKL)=0.
        KKCTL(:)=JK+IKL
    ENDWHERE
 
    ! compute frac_up at JK+KKL
    WHERE (GTEST(:))
      PFRAC_UP(:,JK+IKL)=PEMF(:,JK+IKL)/&
                                      &(SQRT(ZW_UP2(:,JK+IKL))*ZRHO_F(:,JK+IKL))
    ENDWHERE

    ! Updraft fraction must be smaller than XFRAC_UP_MAX
    WHERE (GTEST(:))
      PFRAC_UP(:,JK+IKL)=MIN(PARAMMF%XFRAC_UP_MAX,PFRAC_UP(:,JK+IKL))
    ENDWHERE

    ! When cloudy and non-buoyant, updraft fraction must decrease
    WHERE ((GTEST(:).AND.GTESTETL(:)).AND.GTESTLCL(:))
      PFRAC_UP(:,JK+IKL)=MIN(PFRAC_UP(:,JK+IKL),PFRAC_UP(:,JK))
    ENDWHERE

    ! Mass flux is updated with the new updraft fraction
    IF (OENTR_DETR) PEMF(:,JK+IKL)=PFRAC_UP(:,JK+IKL)*SQRT(ZW_UP2(:,JK+IKL))* &
                                              &ZRHO_F(:,JK+IKL)
    !$mnh_end_expand_where(JIJ=IIJB:IIJE)

  END IF !OENTR_DETR
ENDDO

IF(OENTR_DETR) THEN

  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PW_UP(:,:)=SQRT(ZW_UP2(:,:))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

  !$mnh_expand_array(JIJ=IIJB:IIJE)
  PEMF(:,IKB) =0.
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)

  ! Limits the shallow convection scheme when cloud heigth is higher than 3000m.
  ! To do this, mass flux is multiplied by a coefficient decreasing linearly
  ! from 1 (for clouds of ZDEPTH_MAX1 m of depth) to 0 (for clouds of ZDEPTH_MAX2 m of depth).
  ! This way, all MF fluxes are diminished by this amount.
  ! Diagnosed cloud fraction is also multiplied by the same coefficient.
  !
  DO JIJ=IIJB,IIJE
     PDEPTH(JIJ) = MAX(0., PZZ(JIJ,KKCTL(JIJ)) -  PZZ(JIJ,KKLCL(JIJ)) )
  END DO

  !$mnh_expand_array(JIJ=IIJB:IIJE)
  GWORK1(:)= (GTESTLCL(:) .AND. (PDEPTH(:) > ZDEPTH_MAX1) )
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  DO JK=1,IKT
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    GWORK2(:,JK) = GWORK1(:)
    ZCOEF(:,JK) = (1.-(PDEPTH(:)-ZDEPTH_MAX1)/(ZDEPTH_MAX2-ZDEPTH_MAX1))
    ZCOEF(:,JK)=MIN(MAX(ZCOEF(:,JK),0.),1.)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  ENDDO
  !$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
  WHERE (GWORK2(:,:)) 
    PEMF(:,:)     = PEMF(:,:)     * ZCOEF(:,:)
    PFRAC_UP(:,:) = PFRAC_UP(:,:) * ZCOEF(:,:)
  ENDWHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
ENDIF

IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAFT',1,ZHOOK_HANDLE)
CONTAINS
INCLUDE "th_r_from_thl_rt.func.h"
INCLUDE "compute_frac_ice.func.h"
          SUBROUTINE COMPUTE_ENTR_DETR(D, CST, NEB, PARAMMF,&
                            KK,KKB,KKE,KKL,OTEST,OTESTLCL,&
                            HFRAC_ICE,PFRAC_ICE,PRHODREF,&
                            PPRE_MINUS_HALF,&
                            PPRE_PLUS_HALF,PZZ,PDZZ,&
                            PTHVM,PTHLM,PRTM,PW_UP2,PTH_UP,&
                            PTHL_UP,PRT_UP,PLUP,&
                            PRC_UP,PRI_UP,PTHV_UP,&
                            PRSAT_UP,PRC_MIX,PRI_MIX,      &
                            PENTR,PDETR,PENTR_CLD,PDETR_CLD,&
                            PBUO_INTEG_DRY,PBUO_INTEG_CLD,&
                            PPART_DRY)
!         #############################################################

!!
!!***COMPUTE_ENTR_DETR* - calculates caracteristics of the updraft or downdraft
!!                       using model of the EDMF scheme 
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to compute entrainement and
!!      detrainement at one level of the updraft
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
!!       Book 1 of Meso-NH documentation (chapter Convection)
!!       
!!
!!     AUTHOR
!!     ------
!!    J.Pergaud : 2009
!!
!!    MODIFICATIONS
!!    -------------
!!      Y.Seity (06/2010) Bug correction
!!      V.Masson (09/2010) Optimization
!!      S. Riette april 2011 : ice added, protection against zero divide by Yves Bouteloup
!!                             protection against too big ZPART_DRY, interface modified
!!      S. Riette Jan 2012: support for both order of vertical levels
!!      S. Riette & J. Escobar (11/2013) : remove div by 0 on real*4 case
!!      P.Marguinaud Jun 2012: fix uninitialized variable
!!      P.Marguinaud Nov 2012: fix gfortran bug
!!      S. Riette Apr 2013: bugs correction, rewriting (for optimisation) and
!!                          improvement of continuity at the condensation level
!!      S. Riette Nov 2013: protection against zero divide for min value of dry PDETR
!!      R.Honnert Oct 2016 : Update with AROME
!  P. Wautelet 08/02/2019: bugfix: compute ZEPSI_CLOUD only once and only when it is needed
!!      R. El Khatib 29-Apr-2019 portability fix : compiler may get confused by embricked WHERE statements
!!                          eventually breaking tests with NaN initializations at compile time.
!!                          Replace by IF conditions and traditional DO loops can only improve the performance.
!  P. Wautelet 10/02/2021: bugfix: initialized PPART_DRY everywhere
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST,             ONLY: CST_t
USE MODD_NEB,             ONLY: NEB_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
!
IMPLICIT NONE
!
!                         
!*                    1.1  Declaration of Arguments
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(NEB_t),            INTENT(IN)   :: NEB
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF
!
INTEGER,                INTENT(IN)   :: KK
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
LOGICAL,DIMENSION(D%NIJT),   INTENT(IN)   :: OTEST ! test to see if updraft is running
LOGICAL,DIMENSION(D%NIJT),   INTENT(IN)   :: OTESTLCL !test of condensation 
CHARACTER(LEN=1),       INTENT(IN)   :: HFRAC_ICE ! frac_ice can be compute using
                                              ! Temperature (T) or prescribed
                                              ! (Y)
REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PFRAC_ICE ! fraction of ice
!
!    prognostic variables at t- deltat
!
REAL, DIMENSION(D%NIJT),     INTENT(IN) ::  PRHODREF  !rhodref
REAL, DIMENSION(D%NIJT),     INTENT(IN) ::  PPRE_MINUS_HALF ! Pressure at flux level KK
REAL, DIMENSION(D%NIJT),     INTENT(IN) ::  PPRE_PLUS_HALF ! Pressure at flux level KK+KKL
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PZZ       !  Height at the flux point
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PDZZ       !  metrics coefficient
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PTHVM      ! ThetaV environment 

!
!   thermodynamical variables which are transformed in conservative var.
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PTHLM     ! Thetal
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PRTM      ! total mixing ratio 
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PW_UP2    ! Vertical velocity^2
REAL, DIMENSION(D%NIJT),   INTENT(IN)     ::  PTH_UP,PTHL_UP,PRT_UP  ! updraft properties
REAL, DIMENSION(D%NIJT),   INTENT(IN)     ::  PLUP      ! LUP compute from the ground
REAL, DIMENSION(D%NIJT),   INTENT(IN)     ::  PRC_UP,PRI_UP   ! Updraft cloud content
REAL, DIMENSION(D%NIJT),   INTENT(IN)     ::  PTHV_UP ! Thetav of updraft
REAL, DIMENSION(D%NIJT),   INTENT(IN)     ::  PRSAT_UP ! Mixing ratio at saturation in updraft
REAL, DIMENSION(D%NIJT),   INTENT(INOUT)  ::  PRC_MIX, PRI_MIX      ! Mixture cloud content
REAL, DIMENSION(D%NIJT),   INTENT(OUT)    ::  PENTR     ! Mass flux entrainment of the updraft
REAL, DIMENSION(D%NIJT),   INTENT(OUT)    ::  PDETR     ! Mass flux detrainment of the updraft
REAL, DIMENSION(D%NIJT),   INTENT(OUT)    ::  PENTR_CLD ! Mass flux entrainment of the updraft in cloudy part
REAL, DIMENSION(D%NIJT),   INTENT(OUT)    ::  PDETR_CLD ! Mass flux detrainment of the updraft in cloudy part
REAL, DIMENSION(D%NIJT),   INTENT(OUT)    ::  PBUO_INTEG_DRY, PBUO_INTEG_CLD! Integral Buoyancy
REAL, DIMENSION(D%NIJT),   INTENT(OUT)    ::  PPART_DRY ! ratio of dry part at the transition level
!
!
!                       1.2  Declaration of local variables
!
!     Local array declaration must be put in the compute_updraft subroutine
!     For simplicity all local variables (including scalars) are moved in the compute_updraft subroutine
!

!----------------------------------------------------------------------------------
                        
!                1.3 Initialisation
!                ------------------
  
ZCOEFFMF_CLOUD=PARAMMF%XENTR_MF * CST%XG / PARAMMF%XCRAD_MF
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZG_O_THVREF_ED(:)=CST%XG/PTHVM(:,KK)

ZFRAC_ICE(:)=PFRAC_ICE(:) ! to not modify fraction of ice

ZPRE(:)=PPRE_MINUS_HALF(:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

!                1.4 Estimation of PPART_DRY
DO JIJ=IIJB,IIJE
  IF(OTEST(JIJ) .AND. OTESTLCL(JIJ)) THEN
    !No dry part when condensation level is reached
    PPART_DRY(JIJ)=0.
    ZDZ_STOP(JIJ)=0.
    ZPRE(JIJ)=PPRE_MINUS_HALF(JIJ)
  ELSE IF (OTEST(JIJ) .AND. .NOT. OTESTLCL(JIJ)) THEN
    !Temperature at flux level KK
    ZT=PTH_UP(JIJ)*(PPRE_MINUS_HALF(JIJ)/CST%XP00) ** (CST%XRD/CST%XCPD)
    !Saturating vapor pressure at flux level KK
    ZFOESW = MIN(EXP( CST%XALPW - CST%XBETAW/ZT - CST%XGAMW*LOG(ZT)  ), 0.99*PPRE_MINUS_HALF(JIJ))
    ZFOESI = MIN(EXP( CST%XALPI - CST%XBETAI/ZT - CST%XGAMI*LOG(ZT)  ), 0.99*PPRE_MINUS_HALF(JIJ))
    !Computation of d.Rsat / dP (partial derivations with respect to P and T
    !and use of T=Theta*(P/P0)**(R/Cp) to transform dT into dP with theta_up
    !constant at the vertical)
    ZDRSATODP=(CST%XBETAW/ZT-CST%XGAMW)*(1-ZFRAC_ICE(JIJ))+(CST%XBETAI/ZT-CST%XGAMI)*ZFRAC_ICE(JIJ)
    ZDRSATODP=((CST%XRD/CST%XCPD)*ZDRSATODP-1.)*PRSAT_UP(JIJ)/ &
                &(PPRE_MINUS_HALF(JIJ)-(ZFOESW*(1-ZFRAC_ICE(JIJ)) + ZFOESI*ZFRAC_ICE(JIJ)))
    !Use of d.Rsat / dP and pressure at flux level KK to find pressure (ZPRE)
    !where Rsat is equal to PRT_UP
    ZPRE(JIJ)=PPRE_MINUS_HALF(JIJ)+(PRT_UP(JIJ)-PRSAT_UP(JIJ))/ZDRSATODP
    !Fraction of dry part (computed with pressure and used with heights, no
    !impact found when using log function here and for pressure on flux levels
    !computation)
    PPART_DRY(JIJ)=MAX(0., MIN(1., (PPRE_MINUS_HALF(JIJ)-ZPRE(JIJ))/(PPRE_MINUS_HALF(JIJ)-PPRE_PLUS_HALF(JIJ))))
    !Height above flux level KK of the cloudy part
    ZDZ_STOP(JIJ) = (PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK))*PPART_DRY(JIJ)
  ELSE
    PPART_DRY(JIJ)=0. ! value does not matter, here
  END IF
END DO

!               1.5 Gradient and flux values of thetav
!$mnh_expand_array(JIJ=IIJB:IIJE)
IF(KK/=KKB)THEN
  ZCOEFF_MINUS_HALF(:)=((PTHVM(:,KK)-PTHVM(:,KK-KKL))/PDZZ(:,KK))
  ZTHV_MINUS_HALF(:) = PTHVM(:,KK) - &
                               & ZCOEFF_MINUS_HALF(:)*0.5*(PZZ(:,KK+KKL)-PZZ(:,KK))
ELSE
  ZCOEFF_MINUS_HALF(:)=0.
  ZTHV_MINUS_HALF(:) = PTHVM(:,KK)
ENDIF
ZCOEFF_PLUS_HALF(:)  = ((PTHVM(:,KK+KKL)-PTHVM(:,KK))/PDZZ(:,KK+KKL))
ZTHV_PLUS_HALF(:)  = PTHVM(:,KK) + &
                             & ZCOEFF_PLUS_HALF(:)*0.5*(PZZ(:,KK+KKL)-PZZ(:,KK))
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

!               2  Dry part computation:
!                  Integral buoyancy and computation of PENTR and PDETR for dry part
!               --------------------------------------------------------------------

DO JIJ=IIJB,IIJE
  IF (OTEST(JIJ) .AND. PPART_DRY(JIJ)>0.) THEN
    !Buoyancy computation in two parts to use change of gradient of theta v of environment
    !Between flux level KK and min(mass level, bottom of cloudy part)
    ZDZ=MIN(ZDZ_STOP(JIJ),(PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK))*0.5)
    PBUO_INTEG_DRY(JIJ) = ZG_O_THVREF_ED(JIJ)*ZDZ*&
                (0.5 * (  - ZCOEFF_MINUS_HALF(JIJ))*ZDZ  &
                  - ZTHV_MINUS_HALF(JIJ) + PTHV_UP(JIJ) )

    !Between mass flux KK and bottom of cloudy part (if above mass flux)
    ZDZ=MAX(0., ZDZ_STOP(JIJ)-(PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK))*0.5)
    PBUO_INTEG_DRY(JIJ) = PBUO_INTEG_DRY(JIJ) + ZG_O_THVREF_ED(JIJ)*ZDZ*&
                (0.5 * (  - ZCOEFF_PLUS_HALF(JIJ))*ZDZ &
                  - PTHVM(JIJ,KK) + PTHV_UP(JIJ) )

    !Entr//Detr. computation
    IF (PBUO_INTEG_DRY(JIJ)>=0.) THEN
      PENTR(JIJ) = 0.5/(PARAMMF%XABUO-PARAMMF%XBENTR*PARAMMF%XENTR_DRY)*&
                 LOG(1.+ (2.*(PARAMMF%XABUO-PARAMMF%XBENTR*PARAMMF%XENTR_DRY)/PW_UP2(JIJ,KK))* &
                 PBUO_INTEG_DRY(JIJ))
      PDETR(JIJ) = 0.
    ELSE
      PENTR(JIJ) = 0.
      PDETR(JIJ) = 0.5/(PARAMMF%XABUO)*&
                 LOG(1.+ (2.*(PARAMMF%XABUO)/PW_UP2(JIJ,KK))* &
                 (-PBUO_INTEG_DRY(JIJ)))
    ENDIF
    PENTR(JIJ) = PARAMMF%XENTR_DRY*PENTR(JIJ)/(PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK))    
    PDETR(JIJ) = PARAMMF%XDETR_DRY*PDETR(JIJ)/(PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK))
    !Minimum value of detrainment
    ZWK0D=PLUP(JIJ)-0.5*(PZZ(JIJ,KK)+PZZ(JIJ,KK+KKL))
    ZWK0D=SIGN(MAX(1., ABS(ZWK0D)), ZWK0D) ! ZWK0D must not be zero
    PDETR(JIJ) = MAX(PPART_DRY(JIJ)*PARAMMF%XDETR_LUP/ZWK0D, PDETR(JIJ))
  ELSE
    !No dry part, condensation reached (OTESTLCL)
    PBUO_INTEG_DRY(JIJ) = 0.
    PENTR(JIJ)=0.
    PDETR(JIJ)=0.
  ENDIF
ENDDO

!               3  Wet part computation
!               -----------------------

!               3.1 Integral buoyancy for cloudy part

! Compute theta_v of updraft at flux level KK+KKL                   
!MIX variables are used to avoid declaring new variables
!but we are dealing with updraft and not mixture
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZRCMIX(:)=PRC_UP(:)
ZRIMIX(:)=PRI_UP(:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
CALL TH_R_FROM_THL_RT(CST,NEB,D%NIJT,HFRAC_ICE,ZFRAC_ICE,&
             PPRE_PLUS_HALF,PTHL_UP,PRT_UP,&
             ZTHMIX,ZRVMIX,ZRCMIX,ZRIMIX,&
             ZRSATW_ED, ZRSATI_ED,OOCEAN=.FALSE.,&
             PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZTHV_UP_F2(:) = ZTHMIX(:)*(1.+ZRVORD*ZRVMIX(:))/(1.+PRT_UP(:))
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

! Integral buoyancy for cloudy part
DO JIJ=IIJB,IIJE
  IF(OTEST(JIJ) .AND. PPART_DRY(JIJ)<1.) THEN
    !Gradient of Theta V updraft over the cloudy part, assuming that thetaV updraft don't change
    !between flux level KK and bottom of cloudy part
    ZCOTHVU=(ZTHV_UP_F2(JIJ)-PTHV_UP(JIJ))/((PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK))*(1-PPART_DRY(JIJ)))

    !Computation in two parts to use change of gradient of theta v of environment
    !Between bottom of cloudy part (if under mass level) and mass level KK
    ZDZ=MAX(0., 0.5*(PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK))-ZDZ_STOP(JIJ))
    PBUO_INTEG_CLD(JIJ) = ZG_O_THVREF_ED(JIJ)*ZDZ*&
            (0.5*( ZCOTHVU - ZCOEFF_MINUS_HALF(JIJ))*ZDZ &
              - (PTHVM(JIJ,KK)-ZDZ*ZCOEFF_MINUS_HALF(JIJ)) + PTHV_UP(JIJ) )

    !Between max(mass level, bottom of cloudy part) and flux level KK+KKL
    ZDZ=(PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK))-MAX(ZDZ_STOP(JIJ),0.5*(PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK)))
    PBUO_INTEG_CLD(JIJ) = PBUO_INTEG_CLD(JIJ)+ZG_O_THVREF_ED(JIJ)*ZDZ*&
                      (0.5*( ZCOTHVU - ZCOEFF_PLUS_HALF(JIJ))*ZDZ&
              - (PTHVM(JIJ,KK)+(0.5*((PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK)))-ZDZ)*ZCOEFF_PLUS_HALF(JIJ)) +&
              PTHV_UP(JIJ) )

  ELSE
    !No cloudy part
    PBUO_INTEG_CLD(JIJ)=0.
  END IF
END DO

!               3.2 Critical mixed fraction for KK+KKL flux level (ZKIC_F2) and
!                   for bottom of cloudy part (ZKIC), then a mean for the cloudy part
!                   (put also in ZKIC)
!
!                   computation by estimating unknown  
!                   T^mix r_c^mix and r_i^mix from enthalpy^mix and r_w^mix
!                   We determine the zero crossing of the linear curve
!                   evaluating the derivative using ZMIXF=0.1
                
ZKIC_INIT=0.1  ! starting value for critical mixed fraction for CLoudy Part

!  Compute thetaV of environment at the bottom of cloudy part
!    and cons then non cons. var. of mixture at the bottom of cloudy part

!   JKLIM computed to avoid KKL(KK-KKL) being < KKL*KKB
JKLIM=KKL*MAX(KKL*(KK-KKL),KKL*KKB)
DO JIJ=IIJB,IIJE
  IF(OTEST(JIJ) .AND. PPART_DRY(JIJ)>0.5) THEN
    ZDZ=ZDZ_STOP(JIJ)-0.5*(PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK))
    ZTHV(JIJ)= PTHVM(JIJ,KK)+ZCOEFF_PLUS_HALF(JIJ)*ZDZ
    ZMIXTHL(JIJ) = ZKIC_INIT * &
               (PTHLM(JIJ,KK)+ZDZ*(PTHLM(JIJ,KK+KKL)-PTHLM(JIJ,KK))/PDZZ(JIJ,KK+KKL)) + &
               (1. - ZKIC_INIT)*PTHL_UP(JIJ)
    ZMIXRT(JIJ)  = ZKIC_INIT * &
               (PRTM(JIJ,KK)+ZDZ*(PRTM(JIJ,KK+KKL)-PRTM(JIJ,KK))/PDZZ(JIJ,KK+KKL)) +   &
               (1. - ZKIC_INIT)*PRT_UP(JIJ)
  ELSEIF(OTEST(JIJ)) THEN
    ZDZ=0.5*(PZZ(JIJ,KK+KKL)-PZZ(JIJ,KK))-ZDZ_STOP(JIJ)
    ZTHV(JIJ)= PTHVM(JIJ,KK)-ZCOEFF_MINUS_HALF(JIJ)*ZDZ
    ZMIXTHL(JIJ) = ZKIC_INIT * &
               (PTHLM(JIJ,KK)-ZDZ*(PTHLM(JIJ,KK)-PTHLM(JIJ,JKLIM))/PDZZ(JIJ,KK)) + &
               (1. - ZKIC_INIT)*PTHL_UP(JIJ)
    ZMIXRT(JIJ)  = ZKIC_INIT * &
               (PRTM(JIJ,KK)-ZDZ*(PRTM(JIJ,KK)-PRTM(JIJ,JKLIM))/PDZZ(JIJ,KK)) + &
               (1. - ZKIC_INIT)*PRT_UP(JIJ)
  ELSE
    ZMIXTHL(JIJ) = 300.
    ZMIXRT(JIJ) = 0.1
  ENDIF
ENDDO
CALL TH_R_FROM_THL_RT(CST,NEB,D%NIJT,HFRAC_ICE,ZFRAC_ICE,&
             ZPRE,ZMIXTHL,ZMIXRT,&
             ZTHMIX,ZRVMIX,PRC_MIX,PRI_MIX,&
             ZRSATW_ED, ZRSATI_ED,OOCEAN=.FALSE.,&
             PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZTHVMIX(:) = ZTHMIX(:)*(1.+ZRVORD*ZRVMIX(:))/(1.+ZMIXRT(:))

!  Compute cons then non cons. var. of mixture at the flux level KK+KKL  with initial ZKIC
ZMIXTHL(:) = ZKIC_INIT * 0.5*(PTHLM(:,KK)+PTHLM(:,KK+KKL))+&
                       & (1. - ZKIC_INIT)*PTHL_UP(:)
ZMIXRT(:)  = ZKIC_INIT * 0.5*(PRTM(:,KK)+PRTM(:,KK+KKL))+&
                       & (1. - ZKIC_INIT)*PRT_UP(:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
CALL TH_R_FROM_THL_RT(CST,NEB,D%NIJT,HFRAC_ICE,ZFRAC_ICE,&
             PPRE_PLUS_HALF,ZMIXTHL,ZMIXRT,&
             ZTHMIX,ZRVMIX,PRC_MIX,PRI_MIX,&
             ZRSATW_ED, ZRSATI_ED,OOCEAN=.FALSE.,&
             PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZTHVMIX_F2(:) = ZTHMIX(:)*(1.+ZRVORD*ZRVMIX(:))/(1.+ZMIXRT(:))
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

!Computation of mean ZKIC over the cloudy part
DO JIJ=IIJB,IIJE
  IF (OTEST(JIJ)) THEN
    ! Compute ZKIC at the bottom of cloudy part
    ! Thetav_up at bottom is equal to Thetav_up at flux level KK
    IF (ABS(PTHV_UP(JIJ)-ZTHVMIX(JIJ))<1.E-10) THEN
      ZKIC(JIJ)=1.
    ELSE
      ZKIC(JIJ) = MAX(0.,PTHV_UP(JIJ)-ZTHV(JIJ))*ZKIC_INIT /  &  
                 (PTHV_UP(JIJ)-ZTHVMIX(JIJ))
    END IF
    ! Compute ZKIC_F2 at flux level KK+KKL
    IF (ABS(ZTHV_UP_F2(JIJ)-ZTHVMIX_F2(JIJ))<1.E-10) THEN
      ZKIC_F2(JIJ)=1.
    ELSE
      ZKIC_F2(JIJ) = MAX(0.,ZTHV_UP_F2(JIJ)-ZTHV_PLUS_HALF(JIJ))*ZKIC_INIT /  &  
                 (ZTHV_UP_F2(JIJ)-ZTHVMIX_F2(JIJ))
    END IF
    !Mean ZKIC over the cloudy part
    ZKIC(JIJ)=MAX(MIN(0.5*(ZKIC(JIJ)+ZKIC_F2(JIJ)),1.),0.)
  END IF
END DO

!               3.3 Integration of PDF
!                   According to Kain and Fritsch (1990), we replace delta Mt
!                   in eq. (7) and (8) using eq. (5). Here we compute the ratio
!                   of integrals without computing delta Me

!Constant PDF
!For this PDF, eq. (5) is delta Me=0.5*delta Mt
DO JIJ=IIJB,IIJE
  IF(OTEST(JIJ)) THEN
    ZEPSI(JIJ) = ZKIC(JIJ)**2. !integration multiplied by 2
    ZDELTA(JIJ) = (1.-ZKIC(JIJ))**2. !idem
  ENDIF
ENDDO

!Triangular PDF
!Calculus must be verified before activating this part, but in this state,
!results on ARM case are almost identical
!For this PDF, eq. (5) is also delta Me=0.5*delta Mt
!WHERE(OTEST(:))
!  !Integration multiplied by 2
!  WHERE(ZKIC<0.5)
!    ZEPSI(:)=8.*ZKIC(:)**3/3.
!    ZDELTA(:)=1.-4.*ZKIC(:)**2+8.*ZKIC(:)**3/3.
!  ELSEWHERE
!    ZEPSI(:)=5./3.-4*ZKIC(:)**2+8.*ZKIC(:)**3/3.
!    ZDELTA(:)=8.*(1.-ZKIC(:))**3/3.
!  ENDWHERE
!ENDWHERE

!               3.4 Computation of PENTR and PDETR
DO JIJ=IIJB,IIJE
  IF(OTEST(JIJ)) THEN
    ZEPSI_CLOUD=MIN(ZDELTA(JIJ), ZEPSI(JIJ))
    PENTR_CLD(JIJ) = (1.-PPART_DRY(JIJ))*ZCOEFFMF_CLOUD*PRHODREF(JIJ)*ZEPSI_CLOUD
    PDETR_CLD(JIJ) = (1.-PPART_DRY(JIJ))*ZCOEFFMF_CLOUD*PRHODREF(JIJ)*ZDELTA(JIJ)
    PENTR(JIJ) = PENTR(JIJ)+PENTR_CLD(JIJ)
    PDETR(JIJ) = PDETR(JIJ)+PDETR_CLD(JIJ)
  ELSE
    PENTR_CLD(JIJ) = 0.
    PDETR_CLD(JIJ) = 0.
  ENDIF
ENDDO

END SUBROUTINE COMPUTE_ENTR_DETR
END SUBROUTINE COMPUTE_UPDRAFT
END MODULE MODE_COMPUTE_UPDRAFT

