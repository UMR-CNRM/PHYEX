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
      SUBROUTINE COMPUTE_UPDRAFT(D, CST, NEB, PARAMMF, TURB, CSTURB, &
                                 KSV, HFRAC_ICE,                  &
                                 OENTR_DETR,OMIXUV,               &
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
TYPE(TURB_t),           INTENT(IN)   :: TURB
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
INTEGER,                INTENT(IN)   :: KSV
CHARACTER(LEN=1),       INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
LOGICAL,                INTENT(IN) :: OENTR_DETR! flag to recompute entrainment, detrainment and mass flux
LOGICAL,                INTENT(IN) :: OMIXUV    ! True if mixing of momentum
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

INTEGER  :: JK,JI,JSV          ! loop counters

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

!
IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAFT',0,ZHOOK_HANDLE)

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
  KKLCL(:)=D%NKE
  KKETL(:)=D%NKE
  KKCTL(:)=D%NKE

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
  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
      PRSAT_UP(JI,JK)=PRVM(JI,JK) ! should be initialised correctly but is (normaly) not used
    ENDDO
  ENDDO

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
DO JK=1,D%NKT 
  DO JI=D%NIJB,D%NIJE 
    PTHL_UP(JI,JK)=ZTHLM_F(JI,JK)
    PRT_UP(JI,JK)=ZRTM_F(JI,JK)
    PU_UP(JI,JK)=ZUM_F(JI,JK)
    PV_UP(JI,JK)=ZVM_F(JI,JK)
  ENDDO
ENDDO
DO JSV=1,KSV 
  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
      PSV_UP(JI,JK,JSV)=ZSVM_F(JI,JK,JSV)
    ENDDO
  ENDDO
ENDDO

! Computation or initialisation of updraft characteristics at the KKB level
! thetal_up,rt_up,thetaV_up, w2,Buoyancy term and mass flux (PEMF)
DO JI=D%NIJB,D%NIJE 
  PTHL_UP(JI,D%NKB)= ZTHLM_F(JI,D%NKB)+ &
  & MAX(0.,MIN(ZTMAX,(PSFTH(JI)/SQRT(ZTKEM_F(JI,D%NKB)))* PARAMMF%XALP_PERT))
  PRT_UP(JI,D%NKB) = ZRTM_F(JI,D%NKB)+ &
  & MAX(0.,MIN(ZRMAX,(PSFRV(JI)/SQRT(ZTKEM_F(JI,D%NKB)))* PARAMMF%XALP_PERT)) 
ENDDO

IF (OENTR_DETR) THEN
  CALL MZM_MF(D, PTHM (:,:), ZTHM_F (:,:))
  CALL MZM_MF(D, PPABSM(:,:), ZPRES_F(:,:))
  CALL MZM_MF(D, PRHODREF(:,:), ZRHO_F (:,:))
  CALL MZM_MF(D, PRVM(:,:), ZRVM_F (:,:))

  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
  ! thetav at mass and flux levels
      ZTHVM_F(JI,JK)=ZTHM_F(JI,JK)* &
      &((1.+ZRVORD*ZRVM_F(JI,JK))/(1.+ZRTM_F(JI,JK)))
      ZTHVM(JI,JK)=PTHM(JI,JK)* &
      &((1.+ZRVORD*PRVM(JI,JK))/(1.+PRTM(JI,JK)))
      
      PTHV_UP(JI,JK)=ZTHVM_F(JI,JK)
    ENDDO
  ENDDO

  ZW_UP2(:,:)=0.
  DO JI=D%NIJB,D%NIJE 
    ZW_UP2(JI,D%NKB) = MAX(0.0001,(2./3.)*ZTKEM_F(JI,D%NKB))
    
  ! Computation of non conservative variable for the KKB level of the updraft
  ! (all or nothing ajustement)
    PRC_UP(JI,D%NKB)=0.
    PRI_UP(JI,D%NKB)=0.
  ENDDO
  CALL TH_R_FROM_THL_RT(CST, NEB, D%NIJT, HFRAC_ICE,PFRAC_ICE_UP(:,D%NKB),ZPRES_F(:,D%NKB), &
             PTHL_UP(:,D%NKB),PRT_UP(:,D%NKB),ZTH_UP(:,D%NKB), &
             PRV_UP(:,D%NKB),PRC_UP(:,D%NKB),PRI_UP(:,D%NKB),ZRSATW(:),ZRSATI(:), OOCEAN=.FALSE., &
             PBUF=ZBUF(:,:), KB=D%NIJB, KE=D%NIJE)

  DO JI=D%NIJB,D%NIJE 
  ! compute updraft thevav and buoyancy term at KKB level
    PTHV_UP(JI,D%NKB) = ZTH_UP(JI,D%NKB)*&
    & ((1+ZRVORD*PRV_UP(JI,D%NKB))/(1+PRT_UP(JI,D%NKB)))
  ! compute mean rsat in updraft
    PRSAT_UP(JI,D%NKB) = ZRSATW(JI)*(1-PFRAC_ICE_UP(JI,D%NKB)) + &
    & ZRSATI(JI)*PFRAC_ICE_UP(JI,D%NKB)
  ENDDO
  ! Closure assumption for mass flux at KKB level
  !

  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
      ZG_O_THVREF(JI,JK)=CST%XG/ZTHVM_F(JI,JK)
    ENDDO
  ENDDO

  ! compute L_up
  GLMIX=.TRUE.
  DO JI=D%NIJB,D%NIJE 
    ZTKEM_F(JI,D%NKB)=0.
  ENDDO
  !
  IF(TURB%CTURBLEN=='RM17') THEN
    CALL GZ_M_W_MF(D, PUM, PDZZ, ZWK)
    CALL MZF_MF(D, ZWK, ZDUDZ)
    CALL GZ_M_W_MF(D, PVM, PDZZ, ZWK)
    CALL MZF_MF(D, ZWK, ZDVDZ)
    DO JK=1,D%NKT 
      DO JI=D%NIJB,D%NIJE 
        ZSHEAR(JI,JK) = SQRT(ZDUDZ(JI,JK)**2 + ZDVDZ(JI,JK)**2)
      ENDDO
    ENDDO
  ELSE
    ZSHEAR = 0. !no shear in bl89 mixing length
  END IF
  !
#ifdef REPRO48
  CALL COMPUTE_BL89_ML(D, CST, CSTURB, PDZZ,ZTKEM_F(:,D%NKB),&
                      &ZG_O_THVREF(:,D%NKB),ZTHVM,D%NKB,GLMIX,.TRUE.,ZSHEAR,ZLUP)
#else
  CALL COMPUTE_BL89_ML(D, CST, CSTURB, PDZZ,ZTKEM_F(:,D%NKB),&
                      &ZG_O_THVREF(:,D%NKB),ZTHVM,D%NKB,GLMIX,.FALSE.,ZSHEAR,ZLUP)
#endif
  DO JI=D%NIJB,D%NIJE 
    ZLUP(JI)=MAX(ZLUP(JI),1.E-10)
    
  ! Compute Buoyancy flux at the ground
    ZWTHVSURF(JI) = (ZTHVM_F(JI,D%NKB)/ZTHM_F(JI,D%NKB))*PSFTH(JI)+     &
    (0.61*ZTHM_F(JI,D%NKB))*PSFRV(JI)
    
  ! Mass flux at KKB level (updraft triggered if PSFTH>0.)
    IF (PARAMMF%LGZ) THEN
      IF(PDX==0. .OR. PDY==0.) THEN                                                                                                   
        CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'COMPUTE_UPDRAFT', 'PDX or PDY is NULL with option LGZ!')                                  
      ENDIF
      ZSURF(JI)=TANH(PARAMMF%XGZ*SQRT(PDX*PDY)/ZLUP(JI))
    ELSE
      ZSURF(JI)=1.
    END IF
    IF (ZWTHVSURF(JI)>0.)THEN
      PEMF(JI,D%NKB) = PARAMMF%XCMF * ZSURF(JI) * ZRHO_F(JI,D%NKB) *  &
      ((ZG_O_THVREF(JI,D%NKB))*ZWTHVSURF(JI)*ZLUP(JI))**(1./3.)
      PFRAC_UP(JI,D%NKB)=MIN(PEMF(JI,D%NKB)/(SQRT(ZW_UP2(JI,D%NKB))*ZRHO_F(JI,D%NKB)), &
      &PARAMMF%XFRAC_UP_MAX)
      ZW_UP2(JI,D%NKB)=(PEMF(JI,D%NKB)/(PFRAC_UP(JI,D%NKB)*ZRHO_F(JI,D%NKB)))**2
      GTEST(JI)=.TRUE.
    ELSE
      PEMF(JI,D%NKB) =0.
      GTEST(JI)=.FALSE.
    ENDIF
  ENDDO
ELSE
  DO JI=D%NIJB,D%NIJE 
    GTEST(JI)=PEMF(JI,D%NKB+D%NKL)>0.
  ENDDO
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

DO JK=D%NKB,D%NKE-D%NKL,D%NKL

  ! IF the updraft top is reached for all column, stop the loop on levels
  ITEST=COUNT(GTEST(D%NIJB:D%NIJE))
  IF (ITEST==0) CYCLE

  !       Computation of entrainment and detrainment with KF90
  !       parameterization in clouds and LR01 in subcloud layer


  ! to find the LCL (check if JK is LCL or not)
  DO JI=D%NIJB,D%NIJE 
    IF ((PRC_UP(JI,JK)+PRI_UP(JI,JK)>0.).AND.(.NOT.(GTESTLCL(JI))))THEN
      KKLCL(JI) = JK           
      GTESTLCL(JI)=.TRUE.
    ENDIF
  ENDDO

  ! COMPUTE PENTR and PDETR at mass level JK
  IF (OENTR_DETR) THEN
    IF(JK/=D%NKB) THEN
      DO JI=D%NIJB,D%NIJE 
        ZRC_MIX(JI,JK) = ZRC_MIX(JI,JK-D%NKL) ! guess of Rc of mixture
        ZRI_MIX(JI,JK) = ZRI_MIX(JI,JK-D%NKL) ! guess of Ri of mixture
      ENDDO
    ENDIF
    CALL COMPUTE_ENTR_DETR(D, CST, NEB, PARAMMF, JK,D%NKB,D%NKE,D%NKL,GTEST,GTESTLCL,HFRAC_ICE,PFRAC_ICE_UP(:,JK),&
                           PRHODREF(:,JK),ZPRES_F(:,JK),ZPRES_F(:,JK+D%NKL),&
                           PZZ(:,:),PDZZ(:,:),ZTHVM(:,:),  &
                           PTHLM(:,:),PRTM(:,:),ZW_UP2(:,:),ZTH_UP(:,JK),   &
                           PTHL_UP(:,JK),PRT_UP(:,JK),ZLUP(:),         &
                           PRC_UP(:,JK),PRI_UP(:,JK),PTHV_UP(:,JK),&
                           PRSAT_UP(:,JK),ZRC_MIX(:,JK),ZRI_MIX(:,JK),                 &
                           PENTR(:,JK),PDETR(:,JK),ZENTR_CLD(:,JK),ZDETR_CLD(:,JK),&
                           ZBUO_INTEG_DRY(:,JK), ZBUO_INTEG_CLD(:,JK), &
                           ZPART_DRY(:)   )
    DO JI=D%NIJB,D%NIJE 
      PBUO_INTEG(JI,JK)=ZBUO_INTEG_DRY(JI,JK)+ZBUO_INTEG_CLD(JI,JK)
      
      IF (JK==D%NKB) THEN
        PDETR(JI,JK)=0.
        ZDETR_CLD(JI,JK)=0.
      ENDIF   
      
    !       Computation of updraft characteristics at level JK+KKL
      IF(GTEST(JI))THEN
        ZMIX1(JI)=0.5*(PZZ(JI,JK+D%NKL)-PZZ(JI,JK))*&
        &(PENTR(JI,JK)-PDETR(JI,JK))
        PEMF(JI,JK+D%NKL)=PEMF(JI,JK)*EXP(2*ZMIX1(JI))
      ENDIF
    ENDDO
  ELSE !OENTR_DETR
    DO JI=D%NIJB,D%NIJE 
      GTEST(JI) = (PEMF(JI,JK+D%NKL)>0.)
    ENDDO
  END IF !OENTR_DETR
  
  ! stop the updraft if MF becomes negative
  DO JI=D%NIJB,D%NIJE 
    IF (GTEST(JI).AND.(PEMF(JI,JK+D%NKL)<=0.))THEN
      PEMF(JI,JK+D%NKL)=0.
      KKCTL(JI) = JK+D%NKL
      GTEST(JI)=.FALSE.
      PFRAC_ICE_UP(JI,JK+D%NKL)=PFRAC_ICE_UP(JI,JK)
      PRSAT_UP(JI,JK+D%NKL)=PRSAT_UP(JI,JK)
    ENDIF
  ENDDO

  ! If the updraft did not stop, compute cons updraft characteritics at jk+KKL
  DO JI=D%NIJB,D%NIJE
    IF(GTEST(JI)) THEN
      ZMIX2(JI) = (PZZ(JI,JK+D%NKL)-PZZ(JI,JK))*PENTR(JI,JK) !&
      ZMIX3_CLD(JI) = (PZZ(JI,JK+D%NKL)-PZZ(JI,JK))*(1.-ZPART_DRY(JI))*ZDETR_CLD(JI,JK) !&                   
      ZMIX2_CLD(JI) = (PZZ(JI,JK+D%NKL)-PZZ(JI,JK))*(1.-ZPART_DRY(JI))*ZENTR_CLD(JI,JK)
#ifdef REPRO48                  
      PTHL_UP(JI,JK+D%NKL)=(PTHL_UP(JI,JK)*(1.-0.5*ZMIX2(JI)) + PTHLM(JI,JK)*ZMIX2(JI)) &
                            /(1.+0.5*ZMIX2(JI))   
      PRT_UP(JI,JK+D%NKL) =(PRT_UP (JI,JK)*(1.-0.5*ZMIX2(JI)) + PRTM(JI,JK)*ZMIX2(JI))  &
                            /(1.+0.5*ZMIX2(JI))
#else
      PTHL_UP(JI,JK+D%NKL)=PTHL_UP(JI,JK)*EXP(-ZMIX2(JI)) + PTHLM(JI,JK)*(1-EXP(-ZMIX2(JI)))
      PRT_UP(JI,JK+D%NKL) =PRT_UP (JI,JK)*EXP(-ZMIX2(JI)) +  PRTM(JI,JK)*(1-EXP(-ZMIX2(JI)))
#endif
    ENDIF
  ENDDO
  
  IF(OMIXUV) THEN
    IF(JK/=D%NKB) THEN
      DO JI=D%NIJB,D%NIJE 
        IF(GTEST(JI))THEN
          PU_UP(JI,JK+D%NKL) = (PU_UP(JI,JK)*(1-0.5*ZMIX2(JI)) + &
          &PUM(JI,JK)*ZMIX2(JI)+ &
          0.5*PARAMMF%XPRES_UV*(PZZ(JI,JK+D%NKL)-PZZ(JI,JK))*&
          ((PUM(JI,JK+D%NKL)-PUM(JI,JK))/PDZZ(JI,JK+D%NKL)+&
          (PUM(JI,JK)-PUM(JI,JK-D%NKL))/PDZZ(JI,JK))        )   &
          /(1+0.5*ZMIX2(JI))
          PV_UP(JI,JK+D%NKL) = (PV_UP(JI,JK)*(1-0.5*ZMIX2(JI)) + &
          &PVM(JI,JK)*ZMIX2(JI)+ &
          0.5*PARAMMF%XPRES_UV*(PZZ(JI,JK+D%NKL)-PZZ(JI,JK))*&
          ((PVM(JI,JK+D%NKL)-PVM(JI,JK))/PDZZ(JI,JK+D%NKL)+&
          (PVM(JI,JK)-PVM(JI,JK-D%NKL))/PDZZ(JI,JK))    )   &
          /(1+0.5*ZMIX2(JI))
        ENDIF
      ENDDO
    ELSE
      DO JI=D%NIJB,D%NIJE 
        IF(GTEST(JI))THEN
          PU_UP(JI,JK+D%NKL) = (PU_UP(JI,JK)*(1-0.5*ZMIX2(JI)) + &
          &PUM(JI,JK)*ZMIX2(JI)+ &
          0.5*PARAMMF%XPRES_UV*(PZZ(JI,JK+D%NKL)-PZZ(JI,JK))*&
          ((PUM(JI,JK+D%NKL)-PUM(JI,JK))/PDZZ(JI,JK+D%NKL))        )   &
          /(1+0.5*ZMIX2(JI))
          PV_UP(JI,JK+D%NKL) = (PV_UP(JI,JK)*(1-0.5*ZMIX2(JI)) + &
          &PVM(JI,JK)*ZMIX2(JI)+ &
          0.5*PARAMMF%XPRES_UV*(PZZ(JI,JK+D%NKL)-PZZ(JI,JK))*&
          ((PVM(JI,JK+D%NKL)-PVM(JI,JK))/PDZZ(JI,JK+D%NKL))    )   &
          /(1+0.5*ZMIX2(JI))
        ENDIF
      ENDDO
    ENDIF
  ENDIF !OMIXUV
  DO JSV=1,KSV 
    IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
    DO JI=D%NIJB,D%NIJE 
      IF(GTEST(JI)) THEN
        PSV_UP(JI,JK+D%NKL,JSV) = (PSV_UP(JI,JK,JSV)*(1-0.5*ZMIX2(JI)) + &
        PSVM(JI,JK,JSV)*ZMIX2(JI))  /(1+0.5*ZMIX2(JI))
      ENDIF
    ENDDO
  END DO  
  
  IF (OENTR_DETR) THEN

    ! Compute non cons. var. at level JK+KKL
    DO JI=D%NIJB,D%NIJE 
      ZRC_UP(JI)=PRC_UP(JI,JK) ! guess = level just below
      ZRI_UP(JI)=PRI_UP(JI,JK) ! guess = level just below
    ENDDO
    CALL TH_R_FROM_THL_RT(CST, NEB, D%NIJT, HFRAC_ICE,PFRAC_ICE_UP(:,JK+D%NKL),ZPRES_F(:,JK+D%NKL), &
            PTHL_UP(:,JK+D%NKL),PRT_UP(:,JK+D%NKL),ZTH_UP(:,JK+D%NKL),              &
            ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:), OOCEAN=.FALSE., &
            PBUF=ZBUF(:,:), KB=D%NIJB, KE=D%NIJE)
    DO JI=D%NIJB,D%NIJE 
      IF(GTEST(JI))THEN
        PRC_UP(JI,JK+D%NKL)=ZRC_UP(JI)
        PRV_UP(JI,JK+D%NKL)=ZRV_UP(JI)
        PRI_UP(JI,JK+D%NKL)=ZRI_UP(JI)
        PRSAT_UP(JI,JK+D%NKL) = ZRSATW(JI)*(1-PFRAC_ICE_UP(JI,JK+D%NKL)) + &
        & ZRSATI(JI)*PFRAC_ICE_UP(JI,JK+D%NKL)
      ENDIF
      
    ! Compute the updraft theta_v, buoyancy and w**2 for level JK+KKL
      IF(GTEST(JI))THEN
        PTHV_UP(JI,JK+D%NKL) = ZTH_UP(JI,JK+D%NKL)* &
        & ((1+ZRVORD*PRV_UP(JI,JK+D%NKL))/(1+PRT_UP(JI,JK+D%NKL)))
        IF (ZBUO_INTEG_DRY(JI,JK)>0.)THEN
          ZW_UP2(JI,JK+D%NKL)  = ZW_UP2(JI,JK) + 2.*(PARAMMF%XABUO-PARAMMF%XBENTR*PARAMMF%XENTR_DRY)* &
          &ZBUO_INTEG_DRY(JI,JK)
        ELSE
          ZW_UP2(JI,JK+D%NKL)  = ZW_UP2(JI,JK) + 2.*PARAMMF%XABUO* ZBUO_INTEG_DRY(JI,JK)
        ENDIF
        ZW_UP2(JI,JK+D%NKL)  = ZW_UP2(JI,JK+D%NKL)*(1.-(PARAMMF%XBDETR*ZMIX3_CLD(JI)+ &
        &PARAMMF%XBENTR*ZMIX2_CLD(JI)))&
        /(1.+(PARAMMF%XBDETR*ZMIX3_CLD(JI)+PARAMMF%XBENTR*ZMIX2_CLD(JI))) &
        +2.*(PARAMMF%XABUO)*ZBUO_INTEG_CLD(JI,JK)/ &
        &(1.+(PARAMMF%XBDETR*ZMIX3_CLD(JI)+PARAMMF%XBENTR*ZMIX2_CLD(JI)))
      ENDIF
      
    ! Test if the updraft has reach the ETL
      IF (GTEST(JI).AND.(PBUO_INTEG(JI,JK)<=0.))THEN
        KKETL(JI) = JK+D%NKL
        GTESTETL(JI)=.TRUE.
      ELSE
        GTESTETL(JI)=.FALSE.
      ENDIF
      
    ! Test is we have reached the top of the updraft
      IF (GTEST(JI).AND.((ZW_UP2(JI,JK+D%NKL)<=0.).OR.(PEMF(JI,JK+D%NKL)<=0.)))THEN
        ZW_UP2(JI,JK+D%NKL)=0.
        PEMF(JI,JK+D%NKL)=0.
        GTEST(JI)=.FALSE.
        PTHL_UP(JI,JK+D%NKL)=ZTHLM_F(JI,JK+D%NKL)
        PRT_UP(JI,JK+D%NKL)=ZRTM_F(JI,JK+D%NKL)
        PRC_UP(JI,JK+D%NKL)=0.
        PRI_UP(JI,JK+D%NKL)=0.
        PRV_UP(JI,JK+D%NKL)=0.
        PTHV_UP(JI,JK+D%NKL)=ZTHVM_F(JI,JK+D%NKL)
        PFRAC_UP(JI,JK+D%NKL)=0.
        KKCTL(JI)=JK+D%NKL
      ENDIF
      
    ! compute frac_up at JK+KKL
      IF (GTEST(JI))THEN
        PFRAC_UP(JI,JK+D%NKL)=PEMF(JI,JK+D%NKL)/&
        &(SQRT(ZW_UP2(JI,JK+D%NKL))*ZRHO_F(JI,JK+D%NKL))
      ENDIF
      
    ! Updraft fraction must be smaller than XFRAC_UP_MAX
      IF (GTEST(JI))THEN
        PFRAC_UP(JI,JK+D%NKL)=MIN(PARAMMF%XFRAC_UP_MAX,PFRAC_UP(JI,JK+D%NKL))
      ENDIF
      
    ! When cloudy and non-buoyant, updraft fraction must decrease
      IF ((GTEST(JI).AND.GTESTETL(JI)).AND.GTESTLCL(JI))THEN
        PFRAC_UP(JI,JK+D%NKL)=MIN(PFRAC_UP(JI,JK+D%NKL),PFRAC_UP(JI,JK))
      ENDIF
      
    ! Mass flux is updated with the new updraft fraction
      IF (OENTR_DETR) PEMF(JI,JK+D%NKL)=PFRAC_UP(JI,JK+D%NKL)*SQRT(ZW_UP2(JI,JK+D%NKL))* &
      &ZRHO_F(JI,JK+D%NKL)
    ENDDO

  END IF !OENTR_DETR
ENDDO

IF(OENTR_DETR) THEN

  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
      PW_UP(JI,JK)=SQRT(ZW_UP2(JI,JK))
    ENDDO
  ENDDO

  DO JI=D%NIJB,D%NIJE 
    PEMF(JI,D%NKB) =0.
  ENDDO

  ! Limits the shallow convection scheme when cloud heigth is higher than 3000m.
  ! To do this, mass flux is multiplied by a coefficient decreasing linearly
  ! from 1 (for clouds of ZDEPTH_MAX1 m of depth) to 0 (for clouds of ZDEPTH_MAX2 m of depth).
  ! This way, all MF fluxes are diminished by this amount.
  ! Diagnosed cloud fraction is also multiplied by the same coefficient.
  !
  DO JI=D%NIJB,D%NIJE
     PDEPTH(JI) = MAX(0., PZZ(JI,KKCTL(JI)) -  PZZ(JI,KKLCL(JI)) )
  END DO

  DO JI=D%NIJB,D%NIJE 
    GWORK1(JI)= (GTESTLCL(JI) .AND. (PDEPTH(JI) > ZDEPTH_MAX1) )
  ENDDO
  DO JK=1, D%NKT
    DO JI=D%NIJB,D%NIJE 
      GWORK2(JI,JK) = GWORK1(JI)
      ZCOEF(JI,JK) = (1.-(PDEPTH(JI)-ZDEPTH_MAX1)/(ZDEPTH_MAX2-ZDEPTH_MAX1))
      ZCOEF(JI,JK)=MIN(MAX(ZCOEF(JI,JK),0.),1.)
    ENDDO
  ENDDO
  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
      IF (GWORK2(JI,JK)) THEN
        PEMF(JI,JK)     = PEMF(JI,JK)     * ZCOEF(JI,JK)
        PFRAC_UP(JI,JK) = PFRAC_UP(JI,JK) * ZCOEF(JI,JK)
      ENDIF
    ENDDO
  ENDDO
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
DO JI=D%NIJB,D%NIJE 
  ZG_O_THVREF_ED(JI)=CST%XG/PTHVM(JI,KK)
  
  ZFRAC_ICE(JI)=PFRAC_ICE(JI) ! to not modify fraction of ice
  
  ZPRE(JI)=PPRE_MINUS_HALF(JI)
ENDDO

!                1.4 Estimation of PPART_DRY
DO JI=D%NIJB,D%NIJE
  IF(OTEST(JI) .AND. OTESTLCL(JI)) THEN
    !No dry part when condensation level is reached
    PPART_DRY(JI)=0.
    ZDZ_STOP(JI)=0.
    ZPRE(JI)=PPRE_MINUS_HALF(JI)
  ELSE IF (OTEST(JI) .AND. .NOT. OTESTLCL(JI)) THEN
    !Temperature at flux level KK
    ZT=PTH_UP(JI)*(PPRE_MINUS_HALF(JI)/CST%XP00) ** (CST%XRD/CST%XCPD)
    !Saturating vapor pressure at flux level KK
    ZFOESW = MIN(EXP( CST%XALPW - CST%XBETAW/ZT - CST%XGAMW*LOG(ZT)  ), 0.99*PPRE_MINUS_HALF(JI))
    ZFOESI = MIN(EXP( CST%XALPI - CST%XBETAI/ZT - CST%XGAMI*LOG(ZT)  ), 0.99*PPRE_MINUS_HALF(JI))
    !Computation of d.Rsat / dP (partial derivations with respect to P and T
    !and use of T=Theta*(P/P0)**(R/Cp) to transform dT into dP with theta_up
    !constant at the vertical)
    ZDRSATODP=(CST%XBETAW/ZT-CST%XGAMW)*(1-ZFRAC_ICE(JI))+(CST%XBETAI/ZT-CST%XGAMI)*ZFRAC_ICE(JI)
    ZDRSATODP=((CST%XRD/CST%XCPD)*ZDRSATODP-1.)*PRSAT_UP(JI)/ &
                &(PPRE_MINUS_HALF(JI)-(ZFOESW*(1-ZFRAC_ICE(JI)) + ZFOESI*ZFRAC_ICE(JI)))
    !Use of d.Rsat / dP and pressure at flux level KK to find pressure (ZPRE)
    !where Rsat is equal to PRT_UP
    ZPRE(JI)=PPRE_MINUS_HALF(JI)+(PRT_UP(JI)-PRSAT_UP(JI))/ZDRSATODP
    !Fraction of dry part (computed with pressure and used with heights, no
    !impact found when using log function here and for pressure on flux levels
    !computation)
    PPART_DRY(JI)=MAX(0., MIN(1., (PPRE_MINUS_HALF(JI)-ZPRE(JI))/(PPRE_MINUS_HALF(JI)-PPRE_PLUS_HALF(JI))))
    !Height above flux level KK of the cloudy part
    ZDZ_STOP(JI) = (PZZ(JI,KK+KKL)-PZZ(JI,KK))*PPART_DRY(JI)
  ELSE
    PPART_DRY(JI)=0. ! value does not matter, here
  END IF
END DO

!               1.5 Gradient and flux values of thetav
DO JI=D%NIJB,D%NIJE 
  IF(KK/=KKB)THEN
    ZCOEFF_MINUS_HALF(JI)=((PTHVM(JI,KK)-PTHVM(JI,KK-KKL))/PDZZ(JI,KK))
    ZTHV_MINUS_HALF(JI) = PTHVM(JI,KK) - &
    & ZCOEFF_MINUS_HALF(JI)*0.5*(PZZ(JI,KK+KKL)-PZZ(JI,KK))
  ELSE
    ZCOEFF_MINUS_HALF(JI)=0.
    ZTHV_MINUS_HALF(JI) = PTHVM(JI,KK)
  ENDIF
  ZCOEFF_PLUS_HALF(JI)  = ((PTHVM(JI,KK+KKL)-PTHVM(JI,KK))/PDZZ(JI,KK+KKL))
  ZTHV_PLUS_HALF(JI)  = PTHVM(JI,KK) + &
  & ZCOEFF_PLUS_HALF(JI)*0.5*(PZZ(JI,KK+KKL)-PZZ(JI,KK))
ENDDO

!               2  Dry part computation:
!                  Integral buoyancy and computation of PENTR and PDETR for dry part
!               --------------------------------------------------------------------

DO JI=D%NIJB,D%NIJE
  IF (OTEST(JI) .AND. PPART_DRY(JI)>0.) THEN
    !Buoyancy computation in two parts to use change of gradient of theta v of environment
    !Between flux level KK and min(mass level, bottom of cloudy part)
    ZDZ=MIN(ZDZ_STOP(JI),(PZZ(JI,KK+KKL)-PZZ(JI,KK))*0.5)
    PBUO_INTEG_DRY(JI) = ZG_O_THVREF_ED(JI)*ZDZ*&
                (0.5 * (  - ZCOEFF_MINUS_HALF(JI))*ZDZ  &
                  - ZTHV_MINUS_HALF(JI) + PTHV_UP(JI) )

    !Between mass flux KK and bottom of cloudy part (if above mass flux)
    ZDZ=MAX(0., ZDZ_STOP(JI)-(PZZ(JI,KK+KKL)-PZZ(JI,KK))*0.5)
    PBUO_INTEG_DRY(JI) = PBUO_INTEG_DRY(JI) + ZG_O_THVREF_ED(JI)*ZDZ*&
                (0.5 * (  - ZCOEFF_PLUS_HALF(JI))*ZDZ &
                  - PTHVM(JI,KK) + PTHV_UP(JI) )

    !Entr//Detr. computation
    IF (PBUO_INTEG_DRY(JI)>=0.) THEN
      PENTR(JI) = 0.5/(PARAMMF%XABUO-PARAMMF%XBENTR*PARAMMF%XENTR_DRY)*&
                 LOG(1.+ (2.*(PARAMMF%XABUO-PARAMMF%XBENTR*PARAMMF%XENTR_DRY)/PW_UP2(JI,KK))* &
                 PBUO_INTEG_DRY(JI))
      PDETR(JI) = 0.
    ELSE
      PENTR(JI) = 0.
      PDETR(JI) = 0.5/(PARAMMF%XABUO)*&
                 LOG(1.+ (2.*(PARAMMF%XABUO)/PW_UP2(JI,KK))* &
                 (-PBUO_INTEG_DRY(JI)))
    ENDIF
    PENTR(JI) = PARAMMF%XENTR_DRY*PENTR(JI)/(PZZ(JI,KK+KKL)-PZZ(JI,KK))    
    PDETR(JI) = PARAMMF%XDETR_DRY*PDETR(JI)/(PZZ(JI,KK+KKL)-PZZ(JI,KK))
    !Minimum value of detrainment
    ZWK0D=PLUP(JI)-0.5*(PZZ(JI,KK)+PZZ(JI,KK+KKL))
    ZWK0D=SIGN(MAX(1., ABS(ZWK0D)), ZWK0D) ! ZWK0D must not be zero
    PDETR(JI) = MAX(PPART_DRY(JI)*PARAMMF%XDETR_LUP/ZWK0D, PDETR(JI))
  ELSE
    !No dry part, condensation reached (OTESTLCL)
    PBUO_INTEG_DRY(JI) = 0.
    PENTR(JI)=0.
    PDETR(JI)=0.
  ENDIF
ENDDO

!               3  Wet part computation
!               -----------------------

!               3.1 Integral buoyancy for cloudy part

! Compute theta_v of updraft at flux level KK+KKL                   
!MIX variables are used to avoid declaring new variables
!but we are dealing with updraft and not mixture
DO JI=D%NIJB,D%NIJE 
  ZRCMIX(JI)=PRC_UP(JI)
  ZRIMIX(JI)=PRI_UP(JI)
ENDDO
CALL TH_R_FROM_THL_RT(CST,NEB,D%NIJT,HFRAC_ICE,ZFRAC_ICE,&
             PPRE_PLUS_HALF,PTHL_UP,PRT_UP,&
             ZTHMIX,ZRVMIX,ZRCMIX,ZRIMIX,&
             ZRSATW_ED, ZRSATI_ED,OOCEAN=.FALSE.,&
             PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)
DO JI=D%NIJB,D%NIJE 
  ZTHV_UP_F2(JI) = ZTHMIX(JI)*(1.+ZRVORD*ZRVMIX(JI))/(1.+PRT_UP(JI))
ENDDO

! Integral buoyancy for cloudy part
DO JI=D%NIJB,D%NIJE
  IF(OTEST(JI) .AND. PPART_DRY(JI)<1.) THEN
    !Gradient of Theta V updraft over the cloudy part, assuming that thetaV updraft don't change
    !between flux level KK and bottom of cloudy part
    ZCOTHVU=(ZTHV_UP_F2(JI)-PTHV_UP(JI))/((PZZ(JI,KK+KKL)-PZZ(JI,KK))*(1-PPART_DRY(JI)))

    !Computation in two parts to use change of gradient of theta v of environment
    !Between bottom of cloudy part (if under mass level) and mass level KK
    ZDZ=MAX(0., 0.5*(PZZ(JI,KK+KKL)-PZZ(JI,KK))-ZDZ_STOP(JI))
    PBUO_INTEG_CLD(JI) = ZG_O_THVREF_ED(JI)*ZDZ*&
            (0.5*( ZCOTHVU - ZCOEFF_MINUS_HALF(JI))*ZDZ &
              - (PTHVM(JI,KK)-ZDZ*ZCOEFF_MINUS_HALF(JI)) + PTHV_UP(JI) )

    !Between max(mass level, bottom of cloudy part) and flux level KK+KKL
    ZDZ=(PZZ(JI,KK+KKL)-PZZ(JI,KK))-MAX(ZDZ_STOP(JI),0.5*(PZZ(JI,KK+KKL)-PZZ(JI,KK)))
    PBUO_INTEG_CLD(JI) = PBUO_INTEG_CLD(JI)+ZG_O_THVREF_ED(JI)*ZDZ*&
                      (0.5*( ZCOTHVU - ZCOEFF_PLUS_HALF(JI))*ZDZ&
              - (PTHVM(JI,KK)+(0.5*((PZZ(JI,KK+KKL)-PZZ(JI,KK)))-ZDZ)*ZCOEFF_PLUS_HALF(JI)) +&
              PTHV_UP(JI) )

  ELSE
    !No cloudy part
    PBUO_INTEG_CLD(JI)=0.
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
DO JI=D%NIJB,D%NIJE
  IF(OTEST(JI) .AND. PPART_DRY(JI)>0.5) THEN
    ZDZ=ZDZ_STOP(JI)-0.5*(PZZ(JI,KK+KKL)-PZZ(JI,KK))
    ZTHV(JI)= PTHVM(JI,KK)+ZCOEFF_PLUS_HALF(JI)*ZDZ
    ZMIXTHL(JI) = ZKIC_INIT * &
               (PTHLM(JI,KK)+ZDZ*(PTHLM(JI,KK+KKL)-PTHLM(JI,KK))/PDZZ(JI,KK+KKL)) + &
               (1. - ZKIC_INIT)*PTHL_UP(JI)
    ZMIXRT(JI)  = ZKIC_INIT * &
               (PRTM(JI,KK)+ZDZ*(PRTM(JI,KK+KKL)-PRTM(JI,KK))/PDZZ(JI,KK+KKL)) +   &
               (1. - ZKIC_INIT)*PRT_UP(JI)
  ELSEIF(OTEST(JI)) THEN
    ZDZ=0.5*(PZZ(JI,KK+KKL)-PZZ(JI,KK))-ZDZ_STOP(JI)
    ZTHV(JI)= PTHVM(JI,KK)-ZCOEFF_MINUS_HALF(JI)*ZDZ
    ZMIXTHL(JI) = ZKIC_INIT * &
               (PTHLM(JI,KK)-ZDZ*(PTHLM(JI,KK)-PTHLM(JI,JKLIM))/PDZZ(JI,KK)) + &
               (1. - ZKIC_INIT)*PTHL_UP(JI)
    ZMIXRT(JI)  = ZKIC_INIT * &
               (PRTM(JI,KK)-ZDZ*(PRTM(JI,KK)-PRTM(JI,JKLIM))/PDZZ(JI,KK)) + &
               (1. - ZKIC_INIT)*PRT_UP(JI)
  ELSE
#ifdef REPRO55
    ZMIXTHL(JI) = 0.1
#else
    ZMIXTHL(JI) = 300.
#endif
    ZMIXRT(JI) = 0.1
  ENDIF
ENDDO
CALL TH_R_FROM_THL_RT(CST,NEB,D%NIJT,HFRAC_ICE,ZFRAC_ICE,&
             ZPRE,ZMIXTHL,ZMIXRT,&
             ZTHMIX,ZRVMIX,PRC_MIX,PRI_MIX,&
             ZRSATW_ED, ZRSATI_ED,OOCEAN=.FALSE.,&
             PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)
DO JI=D%NIJB,D%NIJE 
  ZTHVMIX(JI) = ZTHMIX(JI)*(1.+ZRVORD*ZRVMIX(JI))/(1.+ZMIXRT(JI))
  
!  Compute cons then non cons. var. of mixture at the flux level KK+KKL  with initial ZKIC
  ZMIXTHL(JI) = ZKIC_INIT * 0.5*(PTHLM(JI,KK)+PTHLM(JI,KK+KKL))+&
  & (1. - ZKIC_INIT)*PTHL_UP(JI)
  ZMIXRT(JI)  = ZKIC_INIT * 0.5*(PRTM(JI,KK)+PRTM(JI,KK+KKL))+&
  & (1. - ZKIC_INIT)*PRT_UP(JI)
ENDDO
CALL TH_R_FROM_THL_RT(CST,NEB,D%NIJT,HFRAC_ICE,ZFRAC_ICE,&
             PPRE_PLUS_HALF,ZMIXTHL,ZMIXRT,&
             ZTHMIX,ZRVMIX,PRC_MIX,PRI_MIX,&
             ZRSATW_ED, ZRSATI_ED,OOCEAN=.FALSE.,&
             PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)
DO JI=D%NIJB,D%NIJE 
  ZTHVMIX_F2(JI) = ZTHMIX(JI)*(1.+ZRVORD*ZRVMIX(JI))/(1.+ZMIXRT(JI))
ENDDO

!Computation of mean ZKIC over the cloudy part
DO JI=D%NIJB,D%NIJE
  IF (OTEST(JI)) THEN
    ! Compute ZKIC at the bottom of cloudy part
    ! Thetav_up at bottom is equal to Thetav_up at flux level KK
    IF (ABS(PTHV_UP(JI)-ZTHVMIX(JI))<1.E-10) THEN
      ZKIC(JI)=1.
    ELSE
      ZKIC(JI) = MAX(0.,PTHV_UP(JI)-ZTHV(JI))*ZKIC_INIT /  &  
                 (PTHV_UP(JI)-ZTHVMIX(JI))
    END IF
    ! Compute ZKIC_F2 at flux level KK+KKL
    IF (ABS(ZTHV_UP_F2(JI)-ZTHVMIX_F2(JI))<1.E-10) THEN
      ZKIC_F2(JI)=1.
    ELSE
      ZKIC_F2(JI) = MAX(0.,ZTHV_UP_F2(JI)-ZTHV_PLUS_HALF(JI))*ZKIC_INIT /  &  
                 (ZTHV_UP_F2(JI)-ZTHVMIX_F2(JI))
    END IF
    !Mean ZKIC over the cloudy part
    ZKIC(JI)=MAX(MIN(0.5*(ZKIC(JI)+ZKIC_F2(JI)),1.),0.)
  END IF
END DO

!               3.3 Integration of PDF
!                   According to Kain and Fritsch (1990), we replace delta Mt
!                   in eq. (7) and (8) using eq. (5). Here we compute the ratio
!                   of integrals without computing delta Me

!Constant PDF
!For this PDF, eq. (5) is delta Me=0.5*delta Mt
DO JI=D%NIJB,D%NIJE
  IF(OTEST(JI)) THEN
    ZEPSI(JI) = ZKIC(JI)**2. !integration multiplied by 2
    ZDELTA(JI) = (1.-ZKIC(JI))**2. !idem
  ENDIF
ENDDO

!Triangular PDF
!Calculus must be verified before activating this part, but in this state,
!results on ARM case are almost identical
!For this PDF, eq. (5) is also delta Me=0.5*delta Mt
!WHERE(OTEST(D%NIJB:D%NIJE))
!  !Integration multiplied by 2
!  WHERE(ZKIC<0.5)
!    ZEPSI(D%NIJB:D%NIJE)=8.*ZKIC(D%NIJB:D%NIJE)**3/3.
!    ZDELTA(D%NIJB:D%NIJE)=1.-4.*ZKIC(D%NIJB:D%NIJE)**2+8.*ZKIC(D%NIJB:D%NIJE)**3/3.
!  ELSEWHERE
!    ZEPSI(D%NIJB:D%NIJE)=5./3.-4*ZKIC(D%NIJB:D%NIJE)**2+8.*ZKIC(D%NIJB:D%NIJE)**3/3.
!    ZDELTA(D%NIJB:D%NIJE)=8.*(1.-ZKIC(D%NIJB:D%NIJE))**3/3.
!  ENDWHERE
!ENDWHERE

!               3.4 Computation of PENTR and PDETR
DO JI=D%NIJB,D%NIJE
  IF(OTEST(JI)) THEN
    ZEPSI_CLOUD=MIN(ZDELTA(JI), ZEPSI(JI))
    PENTR_CLD(JI) = (1.-PPART_DRY(JI))*ZCOEFFMF_CLOUD*PRHODREF(JI)*ZEPSI_CLOUD
    PDETR_CLD(JI) = (1.-PPART_DRY(JI))*ZCOEFFMF_CLOUD*PRHODREF(JI)*ZDELTA(JI)
    PENTR(JI) = PENTR(JI)+PENTR_CLD(JI)
    PDETR(JI) = PDETR(JI)+PDETR_CLD(JI)
  ELSE
    PENTR_CLD(JI) = 0.
    PDETR_CLD(JI) = 0.
  ENDIF
ENDDO

END SUBROUTINE COMPUTE_ENTR_DETR
END SUBROUTINE COMPUTE_UPDRAFT
END MODULE MODE_COMPUTE_UPDRAFT

