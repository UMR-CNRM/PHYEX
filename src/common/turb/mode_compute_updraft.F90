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
      SUBROUTINE COMPUTE_UPDRAFT(D, CST, NEB, PARAMMF, TURB,      &
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
!
USE MODE_COMPUTE_ENTR_DETR, ONLY: COMPUTE_ENTR_DETR
USE MODE_TH_R_FROM_THL_RT_1D, ONLY: TH_R_FROM_THL_RT_1D
USE MODI_SHUMAN_MF, ONLY: MZM_MF, MZF_MF, GZ_M_W_MF

USE MODE_COMPUTE_BL89_ML, ONLY: COMPUTE_BL89_ML
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
INTEGER,                INTENT(IN)   :: KSV
CHARACTER(LEN=1),       INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
LOGICAL,                INTENT(IN) :: OENTR_DETR! flag to recompute entrainment, detrainment and mass flux
LOGICAL,                INTENT(IN) :: OMIXUV    ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)   :: PZZ       !  Height at the flux point
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)   :: PDZZ      !  Metrics coefficient
 
REAL, DIMENSION(D%NIT),   INTENT(IN)   ::  PSFTH,PSFRV
! normal surface fluxes of theta,rv,(u,v) parallel to the orography
!
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN) ::  PPABSM     ! Pressure at t-dt
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN) ::  PRHODREF   ! dry density of the
                                                  ! reference state
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN) ::  PUM        ! u mean wind
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN) ::  PVM        ! v mean wind
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN) ::  PTKEM      ! TKE at t-dt
!
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   ::  PTHM           ! liquid pot. temp. at t-dt
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   ::  PRVM           ! vapor mixing ratio at t-dt
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   ::  PTHLM,PRTM     ! cons. var. at t-dt

REAL, DIMENSION(D%NIT,D%NKT,KSV), INTENT(IN)   ::  PSVM           ! scalar var. at t-dt

REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT)  ::  PTHL_UP,PRT_UP   ! updraft properties
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT)  ::  PU_UP, PV_UP     ! updraft wind components
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(INOUT)::  PRV_UP,PRC_UP, & ! updraft rv, rc
                                         PRI_UP,PTHV_UP,& ! updraft ri, THv
                                         PW_UP,PFRAC_UP,& ! updraft w, fraction
                                         PFRAC_ICE_UP,&   ! liquid/solid fraction in updraft
                                         PRSAT_UP         ! Rsat

REAL, DIMENSION(D%NIT,D%NKT,KSV), INTENT(OUT)  ::  PSV_UP           ! updraft scalar var. 
                                         
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(INOUT)::  PEMF,PDETR,PENTR ! Mass_flux,
                                                          ! detrainment,entrainment
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(INOUT) :: PBUO_INTEG       ! Integrated Buoyancy 
INTEGER, DIMENSION(D%NIT),  INTENT(INOUT) :: KKLCL,KKETL,KKCTL! LCL, ETL, CTL
REAL, DIMENSION(D%NIT),     INTENT(OUT)   :: PDEPTH           ! Deepness of cloud
REAL,                   INTENT(IN)    :: PDX, PDY
!                       1.2  Declaration of local variables
!
!
! Mean environment variables at t-dt at flux point
REAL, DIMENSION(D%NIT,D%NKT) ::    &
                        ZTHM_F,ZRVM_F                 ! Theta,rv of
                                                      ! updraft environnement
REAL, DIMENSION(D%NIT,D%NKT) ::    &
                        ZRTM_F, ZTHLM_F, ZTKEM_F,&    ! rt, thetal,TKE,pressure,
                        ZUM_F,ZVM_F,ZRHO_F,      &    ! density,momentum
                        ZPRES_F,ZTHVM_F,ZTHVM,   &    ! interpolated at the flux point
                        ZG_O_THVREF,             &    ! g*ThetaV ref
                        ZW_UP2,                  &    ! w**2  of the updraft
                        ZBUO_INTEG_DRY, ZBUO_INTEG_CLD,&! Integrated Buoyancy
                        ZENTR_CLD,ZDETR_CLD           ! wet entrainment and detrainment

REAL, DIMENSION(D%NIT,D%NKT,KSV) :: &
                        ZSVM_F ! scalar variables 

                        
REAL, DIMENSION(D%NIT,D%NKT) ::  &
                        ZTH_UP,                  &    ! updraft THETA 
                        ZRC_MIX, ZRI_MIX              ! guess of Rc and Ri for KF mixture

REAL, DIMENSION(D%NIT,D%NKT) ::  ZCOEF  ! diminution coefficient for too high clouds 
                        
REAL, DIMENSION(D%NIT)            ::  ZWTHVSURF  ! Surface w'thetav'

REAL  :: ZRDORV       ! RD/RV
REAL  :: ZRVORD       ! RV/RD


REAL, DIMENSION(D%NIT) :: ZMIX1,ZMIX2,ZMIX3_CLD,ZMIX2_CLD

REAL, DIMENSION(D%NIT) :: ZLUP         ! Upward Mixing length from the ground

INTEGER  :: JK,JI,JSV          ! loop counters

LOGICAL, DIMENSION(D%NIT) ::  GTEST,GTESTLCL,GTESTETL
                               ! Test if the ascent continue, if LCL or ETL is reached
LOGICAL                          ::  GLMIX 
                               ! To choose upward or downward mixing length
LOGICAL, DIMENSION(D%NIT)              :: GWORK1
LOGICAL, DIMENSION(D%NIT,D%NKT) :: GWORK2

INTEGER  :: ITEST, JLOOP

REAL, DIMENSION(D%NIT) :: ZRC_UP, ZRI_UP, ZRV_UP,&
                                 ZRSATW, ZRSATI,&
                                 ZPART_DRY

REAL  :: ZDEPTH_MAX1, ZDEPTH_MAX2 ! control auto-extinction process

REAL  :: ZTMAX,ZRMAX  ! control value

REAL, DIMENSION(D%NIT) :: ZSURF
REAL, DIMENSION(D%NIT,D%NKT) :: ZSHEAR,ZDUDZ,ZDVDZ ! vertical wind shear
REAL(KIND=JPRB) :: ZHOOK_HANDLE
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
  PRSAT_UP(:,:)=PRVM(:,:) ! should be initialised correctly but is (normaly) not used

  !cloud/dry air mixture cloud content
  ZRC_MIX = 0.
  ZRI_MIX = 0.

END IF

! Initialisation of environment variables at t-dt
! variables at flux level
ZTHLM_F(:,:) = MZM_MF(PTHLM(:,:), D%NKA, D%NKU, D%NKL)
ZRTM_F (:,:) = MZM_MF(PRTM(:,:), D%NKA, D%NKU, D%NKL)
ZUM_F  (:,:) = MZM_MF(PUM(:,:), D%NKA, D%NKU, D%NKL)
ZVM_F  (:,:) = MZM_MF(PVM(:,:), D%NKA, D%NKU, D%NKL)
ZTKEM_F(:,:) = MZM_MF(PTKEM(:,:), D%NKA, D%NKU, D%NKL)

DO JSV=1,KSV
  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  ZSVM_F(:,:,JSV) = MZM_MF(PSVM(:,:,JSV), D%NKA, D%NKU, D%NKL)
END DO
!                     
!          Initialisation of updraft characteristics 
PTHL_UP(:,:)=ZTHLM_F(:,:)
PRT_UP(:,:)=ZRTM_F(:,:)
PU_UP(:,:)=ZUM_F(:,:)
PV_UP(:,:)=ZVM_F(:,:)
PSV_UP(:,:,:)=ZSVM_F(:,:,:)


! Computation or initialisation of updraft characteristics at the KKB level
! thetal_up,rt_up,thetaV_up, w2,Buoyancy term and mass flux (PEMF)

PTHL_UP(:,D%NKB)= ZTHLM_F(:,D%NKB)+MAX(0.,MIN(ZTMAX,(PSFTH(:)/SQRT(ZTKEM_F(:,D%NKB)))*PARAMMF%XALP_PERT))
PRT_UP(:,D%NKB) = ZRTM_F(:,D%NKB)+MAX(0.,MIN(ZRMAX,(PSFRV(:)/SQRT(ZTKEM_F(:,D%NKB)))*PARAMMF%XALP_PERT)) 


IF (OENTR_DETR) THEN
  ZTHM_F (:,:) = MZM_MF(PTHM (:,:), D%NKA, D%NKU, D%NKL)
  ZPRES_F(:,:) = MZM_MF(PPABSM(:,:), D%NKA, D%NKU, D%NKL)
  ZRHO_F (:,:) = MZM_MF(PRHODREF(:,:), D%NKA, D%NKU, D%NKL)
  ZRVM_F (:,:) = MZM_MF(PRVM(:,:), D%NKA, D%NKU, D%NKL)

  ! thetav at mass and flux levels
  ZTHVM_F(:,:)=ZTHM_F(:,:)*((1.+ZRVORD*ZRVM_F(:,:))/(1.+ZRTM_F(:,:)))
  ZTHVM(:,:)=PTHM(:,:)*((1.+ZRVORD*PRVM(:,:))/(1.+PRTM(:,:)))

  PTHV_UP(:,:)=ZTHVM_F(:,:)

  ZW_UP2(:,:)=0.
  ZW_UP2(:,D%NKB) = MAX(0.0001,(2./3.)*ZTKEM_F(:,D%NKB))


  ! Computation of non conservative variable for the KKB level of the updraft
  ! (all or nothing ajustement)
  PRC_UP(:,D%NKB)=0.
  PRI_UP(:,D%NKB)=0.
  CALL TH_R_FROM_THL_RT_1D(CST, NEB, D%NIT, HFRAC_ICE,PFRAC_ICE_UP(:,D%NKB),ZPRES_F(:,D%NKB), &
             PTHL_UP(:,D%NKB),PRT_UP(:,D%NKB),ZTH_UP(:,D%NKB), &
             PRV_UP(:,D%NKB),PRC_UP(:,D%NKB),PRI_UP(:,D%NKB),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.)

  ! compute updraft thevav and buoyancy term at KKB level
  PTHV_UP(:,D%NKB) = ZTH_UP(:,D%NKB)*((1+ZRVORD*PRV_UP(:,D%NKB))/(1+PRT_UP(:,D%NKB)))
  ! compute mean rsat in updraft
  PRSAT_UP(:,D%NKB) = ZRSATW(:)*(1-PFRAC_ICE_UP(:,D%NKB)) + ZRSATI(:)*PFRAC_ICE_UP(:,D%NKB)
                                                            
  ! Closure assumption for mass flux at KKB level
  !

  ZG_O_THVREF(:,:)=CST%XG/ZTHVM_F(:,:)

  ! compute L_up
  GLMIX=.TRUE.
  ZTKEM_F(:,D%NKB)=0.
  !
  IF(TURB%CTURBLEN=='RM17') THEN
    ZDUDZ = MZF_MF(GZ_M_W_MF(PUM,PDZZ, D%NKA, D%NKU, D%NKL), D%NKA, D%NKU, D%NKL)
    ZDVDZ = MZF_MF(GZ_M_W_MF(PVM,PDZZ, D%NKA, D%NKU, D%NKL), D%NKA, D%NKU, D%NKL)
    ZSHEAR = SQRT(ZDUDZ*ZDUDZ + ZDVDZ*ZDVDZ)
  ELSE
    ZSHEAR = 0. !no shear in bl89 mixing length
  END IF
  !
#ifdef REPRO48
  CALL COMPUTE_BL89_ML(D%NKA,D%NKB,D%NKE,D%NKU,D%NKL,PDZZ,ZTKEM_F(:,D%NKB),&
                      &ZG_O_THVREF(:,D%NKB),ZTHVM,D%NKB,GLMIX,.TRUE.,ZSHEAR,ZLUP)
#else
  CALL COMPUTE_BL89_ML(D%NKA,D%NKB,D%NKE,D%NKU,D%NKL,PDZZ,ZTKEM_F(:,D%NKB),&
                      &ZG_O_THVREF(:,D%NKB),ZTHVM,D%NKB,GLMIX,.FALSE.,ZSHEAR,ZLUP)
#endif
  ZLUP(:)=MAX(ZLUP(:),1.E-10)

  ! Compute Buoyancy flux at the ground
  ZWTHVSURF(:) = (ZTHVM_F(:,D%NKB)/ZTHM_F(:,D%NKB))*PSFTH(:)+     &
                (0.61*ZTHM_F(:,D%NKB))*PSFRV(:)

  ! Mass flux at KKB level (updraft triggered if PSFTH>0.)
  IF (PARAMMF%LGZ) THEN
    ZSURF(:)=TANH(PARAMMF%XGZ*SQRT(PDX*PDY)/ZLUP)
  ELSE
    ZSURF(:)=1.
  END IF
  WHERE (ZWTHVSURF(:)>0.)
    PEMF(:,D%NKB) = PARAMMF%XCMF * ZSURF(:) * ZRHO_F(:,D%NKB) *  &
            ((ZG_O_THVREF(:,D%NKB))*ZWTHVSURF*ZLUP)**(1./3.)
    PFRAC_UP(:,D%NKB)=MIN(PEMF(:,D%NKB)/(SQRT(ZW_UP2(:,D%NKB))*ZRHO_F(:,D%NKB)),PARAMMF%XFRAC_UP_MAX)
    ZW_UP2(:,D%NKB)=(PEMF(:,D%NKB)/(PFRAC_UP(:,D%NKB)*ZRHO_F(:,D%NKB)))**2
    GTEST(:)=.TRUE.
  ELSEWHERE
    PEMF(:,D%NKB) =0.
    GTEST(:)=.FALSE.
  ENDWHERE
ELSE
  GTEST(:)=PEMF(:,D%NKB+D%NKL)>0.
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
  ITEST=COUNT(GTEST)
  IF (ITEST==0) CYCLE

!       Computation of entrainment and detrainment with KF90
!       parameterization in clouds and LR01 in subcloud layer


! to find the LCL (check if JK is LCL or not)

  WHERE ((PRC_UP(:,JK)+PRI_UP(:,JK)>0.).AND.(.NOT.(GTESTLCL)))
      KKLCL(:) = JK           
      GTESTLCL(:)=.TRUE.
  ENDWHERE

! COMPUTE PENTR and PDETR at mass level JK
  IF (OENTR_DETR) THEN
    IF(JK/=D%NKB) THEN
      ZRC_MIX(:,JK) = ZRC_MIX(:,JK-D%NKL) ! guess of Rc of mixture
      ZRI_MIX(:,JK) = ZRI_MIX(:,JK-D%NKL) ! guess of Ri of mixture
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
    PBUO_INTEG(:,JK)=ZBUO_INTEG_DRY(:,JK)+ZBUO_INTEG_CLD(:,JK)

    IF (JK==D%NKB) THEN
       PDETR(:,JK)=0.
       ZDETR_CLD(:,JK)=0.
    ENDIF   
 
!       Computation of updraft characteristics at level JK+KKL
    WHERE(GTEST)
      ZMIX1(:)=0.5*(PZZ(:,JK+D%NKL)-PZZ(:,JK))*(PENTR(:,JK)-PDETR(:,JK))
      PEMF(:,JK+D%NKL)=PEMF(:,JK)*EXP(2*ZMIX1(:))
    ENDWHERE
  ELSE
    GTEST(:) = (PEMF(:,JK+D%NKL)>0.)
  END IF 
  

! stop the updraft if MF becomes negative
  WHERE (GTEST.AND.(PEMF(:,JK+D%NKL)<=0.))
    PEMF(:,JK+D%NKL)=0.
    KKCTL(:) = JK+D%NKL
    GTEST(:)=.FALSE.
    PFRAC_ICE_UP(:,JK+D%NKL)=PFRAC_ICE_UP(:,JK)
    PRSAT_UP(:,JK+D%NKL)=PRSAT_UP(:,JK)
  ENDWHERE


! If the updraft did not stop, compute cons updraft characteritics at jk+KKL
  DO JLOOP=1,D%NIT
    IF(GTEST(JLOOP)) THEN
      ZMIX2(JLOOP) = (PZZ(JLOOP,JK+D%NKL)-PZZ(JLOOP,JK))*PENTR(JLOOP,JK) !&
      ZMIX3_CLD(JLOOP) = (PZZ(JLOOP,JK+D%NKL)-PZZ(JLOOP,JK))*(1.-ZPART_DRY(JLOOP))*ZDETR_CLD(JLOOP,JK) !&                   
      ZMIX2_CLD(JLOOP) = (PZZ(JLOOP,JK+D%NKL)-PZZ(JLOOP,JK))*(1.-ZPART_DRY(JLOOP))*ZENTR_CLD(JLOOP,JK)
#ifdef REPRO48                  
      PTHL_UP(JLOOP,JK+D%NKL)=(PTHL_UP(JLOOP,JK)*(1.-0.5*ZMIX2(JLOOP)) + PTHLM(JLOOP,JK)*ZMIX2(JLOOP)) &
                            /(1.+0.5*ZMIX2(JLOOP))   
      PRT_UP(JLOOP,JK+D%NKL) =(PRT_UP (JLOOP,JK)*(1.-0.5*ZMIX2(JLOOP)) + PRTM(JLOOP,JK)*ZMIX2(JLOOP))  &
                            /(1.+0.5*ZMIX2(JLOOP))
#else
      PTHL_UP(JLOOP,JK+D%NKL)=PTHL_UP(JLOOP,JK)*EXP(-ZMIX2(JLOOP)) + PTHLM(JLOOP,JK)*(1-EXP(-ZMIX2(JLOOP)))
      PRT_UP(JLOOP,JK+D%NKL) =PRT_UP (JLOOP,JK)*EXP(-ZMIX2(JLOOP)) +  PRTM(JLOOP,JK)*(1-EXP(-ZMIX2(JLOOP)))
#endif
    ENDIF
  ENDDO
  

  IF(OMIXUV) THEN
    IF(JK/=D%NKB) THEN
      WHERE(GTEST)
        PU_UP(:,JK+D%NKL) = (PU_UP (:,JK)*(1-0.5*ZMIX2(:)) + PUM(:,JK)*ZMIX2(:)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(:,JK+D%NKL)-PZZ(:,JK))*&
                          ((PUM(:,JK+D%NKL)-PUM(:,JK))/PDZZ(:,JK+D%NKL)+&
                           (PUM(:,JK)-PUM(:,JK-D%NKL))/PDZZ(:,JK))        )   &
                          /(1+0.5*ZMIX2(:))
        PV_UP(:,JK+D%NKL) = (PV_UP (:,JK)*(1-0.5*ZMIX2(:)) + PVM(:,JK)*ZMIX2(:)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(:,JK+D%NKL)-PZZ(:,JK))*&
                          ((PVM(:,JK+D%NKL)-PVM(:,JK))/PDZZ(:,JK+D%NKL)+&
                           (PVM(:,JK)-PVM(:,JK-D%NKL))/PDZZ(:,JK))    )   &
                          /(1+0.5*ZMIX2(:))
      ENDWHERE
    ELSE
      WHERE(GTEST)
        PU_UP(:,JK+D%NKL) = (PU_UP (:,JK)*(1-0.5*ZMIX2(:)) + PUM(:,JK)*ZMIX2(:)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(:,JK+D%NKL)-PZZ(:,JK))*&
                          ((PUM(:,JK+D%NKL)-PUM(:,JK))/PDZZ(:,JK+D%NKL))        )   &
                          /(1+0.5*ZMIX2(:))
        PV_UP(:,JK+D%NKL) = (PV_UP (:,JK)*(1-0.5*ZMIX2(:)) + PVM(:,JK)*ZMIX2(:)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(:,JK+D%NKL)-PZZ(:,JK))*&
                          ((PVM(:,JK+D%NKL)-PVM(:,JK))/PDZZ(:,JK+D%NKL))    )   &
                          /(1+0.5*ZMIX2(:))
      ENDWHERE

    ENDIF
  ENDIF
  DO JSV=1,KSV 
     IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
      WHERE(GTEST) 
           PSV_UP(:,JK+D%NKL,JSV) = (PSV_UP (:,JK,JSV)*(1-0.5*ZMIX2(:)) + &
                        PSVM(:,JK,JSV)*ZMIX2(:))  /(1+0.5*ZMIX2(:))
      ENDWHERE                        
  END DO  
  
 IF (OENTR_DETR) THEN

! Compute non cons. var. at level JK+KKL
  ZRC_UP(:)=PRC_UP(:,JK) ! guess = level just below
  ZRI_UP(:)=PRI_UP(:,JK) ! guess = level just below
  CALL TH_R_FROM_THL_RT_1D(CST, NEB, D%NIT, HFRAC_ICE,PFRAC_ICE_UP(:,JK+D%NKL),ZPRES_F(:,JK+D%NKL), &
          PTHL_UP(:,JK+D%NKL),PRT_UP(:,JK+D%NKL),ZTH_UP(:,JK+D%NKL),              &
          ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:), OOCEAN=.FALSE.)
  WHERE(GTEST)
    PRC_UP(:,JK+D%NKL)=ZRC_UP(:)
    PRV_UP(:,JK+D%NKL)=ZRV_UP(:)
    PRI_UP(:,JK+D%NKL)=ZRI_UP(:)
    PRSAT_UP(:,JK+D%NKL) = ZRSATW(:)*(1-PFRAC_ICE_UP(:,JK+D%NKL)) + ZRSATI(:)*PFRAC_ICE_UP(:,JK+D%NKL)
  ENDWHERE
  

! Compute the updraft theta_v, buoyancy and w**2 for level JK+KKL
  WHERE(GTEST)
    PTHV_UP(:,JK+D%NKL) = ZTH_UP(:,JK+D%NKL)*((1+ZRVORD*PRV_UP(:,JK+D%NKL))/(1+PRT_UP(:,JK+D%NKL)))
    WHERE (ZBUO_INTEG_DRY(:,JK)>0.)
      ZW_UP2(:,JK+D%NKL)  = ZW_UP2(:,JK) + 2.*(PARAMMF%XABUO-PARAMMF%XBENTR*PARAMMF%XENTR_DRY)* ZBUO_INTEG_DRY(:,JK)
    ELSEWHERE
      ZW_UP2(:,JK+D%NKL)  = ZW_UP2(:,JK) + 2.*PARAMMF%XABUO* ZBUO_INTEG_DRY(:,JK)
    ENDWHERE
    ZW_UP2(:,JK+D%NKL)  = ZW_UP2(:,JK+D%NKL)*(1.-(PARAMMF%XBDETR*ZMIX3_CLD(:)+PARAMMF%XBENTR*ZMIX2_CLD(:)))&
            /(1.+(PARAMMF%XBDETR*ZMIX3_CLD(:)+PARAMMF%XBENTR*ZMIX2_CLD(:))) &
            +2.*(PARAMMF%XABUO)*ZBUO_INTEG_CLD(:,JK)/(1.+(PARAMMF%XBDETR*ZMIX3_CLD(:)+PARAMMF%XBENTR*ZMIX2_CLD(:)))
 ENDWHERE


  ! Test if the updraft has reach the ETL
  GTESTETL(:)=.FALSE.
  WHERE (GTEST.AND.(PBUO_INTEG(:,JK)<=0.))
      KKETL(:) = JK+D%NKL
      GTESTETL(:)=.TRUE.
  ENDWHERE

  ! Test is we have reached the top of the updraft
  WHERE (GTEST.AND.((ZW_UP2(:,JK+D%NKL)<=0.).OR.(PEMF(:,JK+D%NKL)<=0.)))
      ZW_UP2(:,JK+D%NKL)=0.
      PEMF(:,JK+D%NKL)=0.
      GTEST(:)=.FALSE.
      PTHL_UP(:,JK+D%NKL)=ZTHLM_F(:,JK+D%NKL)
      PRT_UP(:,JK+D%NKL)=ZRTM_F(:,JK+D%NKL)
      PRC_UP(:,JK+D%NKL)=0.
      PRI_UP(:,JK+D%NKL)=0.
      PRV_UP(:,JK+D%NKL)=0.
      PTHV_UP(:,JK+D%NKL)=ZTHVM_F(:,JK+D%NKL)
      PFRAC_UP(:,JK+D%NKL)=0.
      KKCTL(:)=JK+D%NKL
  ENDWHERE
 
  ! compute frac_up at JK+KKL
  WHERE (GTEST)
    PFRAC_UP(:,JK+D%NKL)=PEMF(:,JK+D%NKL)/(SQRT(ZW_UP2(:,JK+D%NKL))*ZRHO_F(:,JK+D%NKL))
  ENDWHERE

  ! Updraft fraction must be smaller than XFRAC_UP_MAX
  WHERE (GTEST)
    PFRAC_UP(:,JK+D%NKL)=MIN(PARAMMF%XFRAC_UP_MAX,PFRAC_UP(:,JK+D%NKL))
  ENDWHERE

  ! When cloudy and non-buoyant, updraft fraction must decrease
  WHERE ((GTEST.AND.GTESTETL).AND.GTESTLCL)
    PFRAC_UP(:,JK+D%NKL)=MIN(PFRAC_UP(:,JK+D%NKL),PFRAC_UP(:,JK))
  ENDWHERE

  ! Mass flux is updated with the new updraft fraction
  IF (OENTR_DETR) PEMF(:,JK+D%NKL)=PFRAC_UP(:,JK+D%NKL)*SQRT(ZW_UP2(:,JK+D%NKL))*ZRHO_F(:,JK+D%NKL)

 END IF

ENDDO

IF(OENTR_DETR) THEN

  PW_UP(:,:)=SQRT(ZW_UP2(:,:))

  PEMF(:,D%NKB) =0.

! Limits the shallow convection scheme when cloud heigth is higher than 3000m.
! To do this, mass flux is multiplied by a coefficient decreasing linearly
! from 1 (for clouds of ZDEPTH_MAX1 m of depth) to 0 (for clouds of ZDEPTH_MAX2 m of depth).
! This way, all MF fluxes are diminished by this amount.
! Diagnosed cloud fraction is also multiplied by the same coefficient.
!
  DO JI=1,D%NIT
     PDEPTH(JI) = MAX(0., PZZ(JI,KKCTL(JI)) -  PZZ(JI,KKLCL(JI)) )
  END DO

  GWORK1(:)= (GTESTLCL(:) .AND. (PDEPTH(:) > ZDEPTH_MAX1) )
  GWORK2(:,:) = SPREAD( GWORK1(:), DIM=2, NCOPIES=MAX(D%NKU,D%NKA) )
  ZCOEF(:,:) = SPREAD( (1.-(PDEPTH(:)-ZDEPTH_MAX1)/(ZDEPTH_MAX2-ZDEPTH_MAX1)), DIM=2, NCOPIES=D%NKT)
  ZCOEF=MIN(MAX(ZCOEF,0.),1.)
  WHERE (GWORK2) 
    PEMF(:,:)     = PEMF(:,:)     * ZCOEF(:,:)
    PFRAC_UP(:,:) = PFRAC_UP(:,:) * ZCOEF(:,:)
  ENDWHERE
ENDIF

IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAFT',1,ZHOOK_HANDLE)
END SUBROUTINE COMPUTE_UPDRAFT
END MODULE MODE_COMPUTE_UPDRAFT
