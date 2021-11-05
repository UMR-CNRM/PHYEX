!MNH_LIC Copyright 2004-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
     MODULE MODI_COMPUTE_UPDRAFT
!    ###########################
!
INTERFACE
!
!     #################################################################
      SUBROUTINE COMPUTE_UPDRAFT(KKA,KKB,KKE,KKU,KKL, HFRAC_ICE,  &
                                 OENTR_DETR,OMIXUV,               &
                                 ONOMIXLG,KSV_LGBEG,KSV_LGEND,    &
                                 PZZ,PDZZ,                        &
                                 PSFTH,PSFRV,                     &
                                 PPABSM,PRHODREF,PUM,PVM,PTKEM,   &
                                 PTHM,PRVM,PTHLM,PRTM,            &
                                 PSVM,PTHL_UP,PRT_UP,             &
                                 PRV_UP,PRC_UP,PRI_UP,PTHV_UP,    &
                                 PW_UP,PU_UP, PV_UP, PSV_UP,      &
                                 PFRAC_UP,PFRAC_ICE_UP,PRSAT_UP,  &
                                 PEMF,PDETR,PENTR,                &
                                 PBUO_INTEG,KKLCL,KKETL,KKCTL,    &
                                 PDEPTH)
!     #################################################################
!
!*                    1.1  Declaration of Arguments
!
!
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
CHARACTER(len=1),       INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
LOGICAL,                INTENT(IN) :: OENTR_DETR! flag to recompute entrainment, detrainment and mass flux
LOGICAL,                INTENT(IN) :: OMIXUV    ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL, DIMENSION(:,:), INTENT(IN)   :: PZZ       !  Height at the flux point
REAL, DIMENSION(:,:), INTENT(IN)   :: PDZZ      !  Metrics coefficient
 
REAL, DIMENSION(:),   INTENT(IN)   ::  PSFTH,PSFRV
! normal surface fluxes of theta,rv,(u,v) parallel to the orography
!
REAL, DIMENSION(:,:),   INTENT(IN) ::  PPABSM     ! Pressure at t-dt
REAL, DIMENSION(:,:),   INTENT(IN) ::  PRHODREF   ! dry density of the
                                                  ! reference state
REAL, DIMENSION(:,:),   INTENT(IN) ::  PUM        ! u mean wind
REAL, DIMENSION(:,:),   INTENT(IN) ::  PVM        ! v mean wind
REAL, DIMENSION(:,:),   INTENT(IN) ::  PTKEM      ! TKE at t-dt
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHM           ! liquid pot. temp. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PRVM           ! vapor mixing ratio at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHLM,PRTM     ! cons. var. at t-dt

REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSVM           ! scalar var. at t-dt

REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PTHL_UP,PRT_UP   ! updraft properties
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PU_UP, PV_UP     ! updraft wind components
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PRV_UP,PRC_UP, & ! updraft rv, rc
                                         PRI_UP,PTHV_UP,& ! updraft ri, THv
                                         PW_UP,PFRAC_UP,& ! updraft w, fraction
                                         PFRAC_ICE_UP,&   ! liquid/solid fraction in updraft
                                         PRSAT_UP         ! Rsat

REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PSV_UP           ! updraft scalar var. 
                                         
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PEMF,PDETR,PENTR ! Mass_flux,
                                                          ! entrainment, detrainment
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PBUO_INTEG       ! Integrated Buoyancy 
INTEGER, DIMENSION(:),  INTENT(INOUT)::  KKLCL,KKETL,KKCTL! LCL, ETL, CTL                                           
REAL, DIMENSION(:),     INTENT(OUT)   :: PDEPTH           ! Deepness of cloud


END SUBROUTINE COMPUTE_UPDRAFT

END INTERFACE
!
END MODULE MODI_COMPUTE_UPDRAFT
!     ######spl
      SUBROUTINE COMPUTE_UPDRAFT(KKA,KKB,KKE,KKU,KKL,HFRAC_ICE,   &
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
                                 PDEPTH     )

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
USE MODD_CST
USE MODD_PARAM_MFSHALL_n
USE MODD_TURB_n, ONLY : CTURBLEN

USE MODI_COMPUTE_ENTR_DETR
USE MODI_TH_R_FROM_THL_RT_1D
USE MODI_SHUMAN_MF

USE MODI_COMPUTE_BL89_ML
USE MODD_GRID_n, ONLY : XDXHAT, XDYHAT


IMPLICIT NONE

!*                    1.1  Declaration of Arguments
!
!
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
CHARACTER(len=1),       INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
LOGICAL,                INTENT(IN) :: OENTR_DETR! flag to recompute entrainment, detrainment and mass flux
LOGICAL,                INTENT(IN) :: OMIXUV    ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL, DIMENSION(:,:), INTENT(IN)   :: PZZ       !  Height at the flux point
REAL, DIMENSION(:,:), INTENT(IN)   :: PDZZ      !  Metrics coefficient
 
REAL, DIMENSION(:),   INTENT(IN)   ::  PSFTH,PSFRV
! normal surface fluxes of theta,rv,(u,v) parallel to the orography
!
REAL, DIMENSION(:,:),   INTENT(IN) ::  PPABSM     ! Pressure at t-dt
REAL, DIMENSION(:,:),   INTENT(IN) ::  PRHODREF   ! dry density of the
                                                  ! reference state
REAL, DIMENSION(:,:),   INTENT(IN) ::  PUM        ! u mean wind
REAL, DIMENSION(:,:),   INTENT(IN) ::  PVM        ! v mean wind
REAL, DIMENSION(:,:),   INTENT(IN) ::  PTKEM      ! TKE at t-dt
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHM           ! liquid pot. temp. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PRVM           ! vapor mixing ratio at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHLM,PRTM     ! cons. var. at t-dt

REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSVM           ! scalar var. at t-dt

REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PTHL_UP,PRT_UP   ! updraft properties
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PU_UP, PV_UP     ! updraft wind components
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PRV_UP,PRC_UP, & ! updraft rv, rc
                                         PRI_UP,PTHV_UP,& ! updraft ri, THv
                                         PW_UP,PFRAC_UP,& ! updraft w, fraction
                                         PFRAC_ICE_UP,&   ! liquid/solid fraction in updraft
                                         PRSAT_UP         ! Rsat

REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PSV_UP           ! updraft scalar var. 
                                         
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PEMF,PDETR,PENTR ! Mass_flux,
                                                          ! detrainment,entrainment
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PBUO_INTEG       ! Integrated Buoyancy 
INTEGER, DIMENSION(:),  INTENT(INOUT) :: KKLCL,KKETL,KKCTL! LCL, ETL, CTL
REAL, DIMENSION(:),     INTENT(OUT)   :: PDEPTH           ! Deepness of cloud
!                       1.2  Declaration of local variables
!
!
! Mean environment variables at t-dt at flux point
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::    &
                        ZTHM_F,ZRVM_F                 ! Theta,rv of
                                                      ! updraft environnement
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::    &
                        ZRTM_F, ZTHLM_F, ZTKEM_F,&    ! rt, thetal,TKE,pressure,
                        ZUM_F,ZVM_F,ZRHO_F,      &    ! density,momentum
                        ZPRES_F,ZTHVM_F,ZTHVM,   &    ! interpolated at the flux point
                        ZG_O_THVREF,             &    ! g*ThetaV ref
                        ZW_UP2,                  &    ! w**2  of the updraft
                        ZBUO_INTEG_DRY, ZBUO_INTEG_CLD,&! Integrated Buoyancy
                        ZENTR_CLD,ZDETR_CLD           ! wet entrainment and detrainment

REAL, DIMENSION(SIZE(PSVM,1),SIZE(PTHM,2),SIZE(PSVM,3)) :: &
                        ZSVM_F ! scalar variables 

                        
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::  &
                        ZTH_UP,                  &    ! updraft THETA 
                        ZRC_MIX, ZRI_MIX              ! guess of Rc and Ri for KF mixture

REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::  ZCOEF  ! diminution coefficient for too high clouds 
                        
REAL, DIMENSION(SIZE(PSFTH,1) )            ::  ZWTHVSURF  ! Surface w'thetav'

REAL  :: ZRDORV       ! RD/RV
REAL  :: ZRVORD       ! RV/RD


REAL, DIMENSION(SIZE(PTHM,1)) :: ZMIX1,ZMIX2,ZMIX3_CLD,ZMIX2_CLD

REAL, DIMENSION(SIZE(PTHM,1)) :: ZLUP         ! Upward Mixing length from the ground

INTEGER  :: ISV                ! Number of scalar variables                               
INTEGER  :: JK,JI,JSV          ! loop counters

LOGICAL, DIMENSION(SIZE(PTHM,1)) ::  GTEST,GTESTLCL,GTESTETL
                               ! Test if the ascent continue, if LCL or ETL is reached
LOGICAL                          ::  GLMIX 
                               ! To choose upward or downward mixing length
LOGICAL, DIMENSION(SIZE(PTHM,1))              :: GWORK1
LOGICAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: GWORK2

INTEGER  :: ITEST,JLOOP

REAL, DIMENSION(SIZE(PTHM,1)) :: ZRC_UP, ZRI_UP, ZRV_UP, ZRSATW, ZRSATI,&
                                 ZPART_DRY

REAL  :: ZDEPTH_MAX1, ZDEPTH_MAX2 ! control auto-extinction process

REAL  :: ZTMAX,ZRMAX  ! control value

REAL, DIMENSION(SIZE(PTHM,1)) :: ZSURF
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZSHEAR,ZDUDZ,ZDVDZ ! vertical wind shear
! Thresholds for the  perturbation of
! theta_l and r_t at the first level of the updraft
ZTMAX=2.0
ZRMAX=1.E-3
!------------------------------------------------------------------------

!                     INITIALISATION

! Initialisation of the constants   
ZRDORV   = XRD / XRV   !=0.622
ZRVORD   = (XRV / XRD) 

ZDEPTH_MAX1=3000. ! clouds with depth inferior to this value are keeped untouched
ZDEPTH_MAX2=4000. ! clouds with depth superior to this value are suppressed

!                 Local variables, internal domain
!number of scalar variables
ISV=SIZE(PSVM,3)

IF (OENTR_DETR) THEN
  ! Initialisation of intersesting Level :LCL,ETL,CTL
  KKLCL(:)=KKE
  KKETL(:)=KKE
  KKCTL(:)=KKE

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
ZTHLM_F(:,:) = MZM_MF(KKA,KKU,KKL,PTHLM(:,:))
ZRTM_F (:,:) = MZM_MF(KKA,KKU,KKL,PRTM(:,:))
ZUM_F  (:,:) = MZM_MF(KKA,KKU,KKL,PUM(:,:))
ZVM_F  (:,:) = MZM_MF(KKA,KKU,KKL,PVM(:,:))
ZTKEM_F(:,:) = MZM_MF(KKA,KKU,KKL,PTKEM(:,:))

DO JSV=1,ISV
  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  ZSVM_F(:,:,JSV) = MZM_MF(KKA,KKU,KKL,PSVM(:,:,JSV))
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

PTHL_UP(:,KKB)= ZTHLM_F(:,KKB)+MAX(0.,MIN(ZTMAX,(PSFTH(:)/SQRT(ZTKEM_F(:,KKB)))*XALP_PERT))
PRT_UP(:,KKB) = ZRTM_F(:,KKB)+MAX(0.,MIN(ZRMAX,(PSFRV(:)/SQRT(ZTKEM_F(:,KKB)))*XALP_PERT)) 


IF (OENTR_DETR) THEN
  ZTHM_F (:,:) = MZM_MF(KKA,KKU,KKL,PTHM (:,:))
  ZPRES_F(:,:) = MZM_MF(KKA,KKU,KKL,PPABSM(:,:))
  ZRHO_F (:,:) = MZM_MF(KKA,KKU,KKL,PRHODREF(:,:))
  ZRVM_F (:,:) = MZM_MF(KKA,KKU,KKL,PRVM(:,:))

  ! thetav at mass and flux levels
  ZTHVM_F(:,:)=ZTHM_F(:,:)*((1.+ZRVORD*ZRVM_F(:,:))/(1.+ZRTM_F(:,:)))
  ZTHVM(:,:)=PTHM(:,:)*((1.+ZRVORD*PRVM(:,:))/(1.+PRTM(:,:)))

  PTHV_UP(:,:)=ZTHVM_F(:,:)

  ZW_UP2(:,:)=0.
  ZW_UP2(:,KKB) = MAX(0.0001,(2./3.)*ZTKEM_F(:,KKB))


  ! Computation of non conservative variable for the KKB level of the updraft
  ! (all or nothing ajustement)
  PRC_UP(:,KKB)=0.
  PRI_UP(:,KKB)=0.
  CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE_UP(:,KKB),ZPRES_F(:,KKB), &
             PTHL_UP(:,KKB),PRT_UP(:,KKB),ZTH_UP(:,KKB), &
             PRV_UP(:,KKB),PRC_UP(:,KKB),PRI_UP(:,KKB),ZRSATW(:),ZRSATI(:))

  ! compute updraft thevav and buoyancy term at KKB level
  PTHV_UP(:,KKB) = ZTH_UP(:,KKB)*((1+ZRVORD*PRV_UP(:,KKB))/(1+PRT_UP(:,KKB)))
  ! compute mean rsat in updraft
  PRSAT_UP(:,KKB) = ZRSATW(:)*(1-PFRAC_ICE_UP(:,KKB)) + ZRSATI(:)*PFRAC_ICE_UP(:,KKB)
                                                            
  ! Closure assumption for mass flux at KKB level
  !

  ZG_O_THVREF=XG/ZTHVM_F

  ! compute L_up
  GLMIX=.TRUE.
  ZTKEM_F(:,KKB)=0.
  !
  IF(CTURBLEN=='RM17') THEN
    ZDUDZ = MZF_MF(KKA,KKU,KKL,GZ_M_W_MF(KKA,KKU,KKL,PUM,PDZZ))
    ZDVDZ = MZF_MF(KKA,KKU,KKL,GZ_M_W_MF(KKA,KKU,KKL,PVM,PDZZ))
    ZSHEAR = SQRT(ZDUDZ*ZDUDZ + ZDVDZ*ZDVDZ)
  ELSE
    ZSHEAR = 0. !no shear in bl89 mixing length
  END IF
 ! 
  CALL COMPUTE_BL89_ML(KKA,KKB,KKE,KKU,KKL,PDZZ,ZTKEM_F(:,KKB),ZG_O_THVREF(:,KKB),ZTHVM,KKB,GLMIX,.FALSE.,ZSHEAR,ZLUP)
  ZLUP(:)=MAX(ZLUP(:),1.E-10)

  ! Compute Buoyancy flux at the ground
  ZWTHVSURF(:) = (ZTHVM_F(:,KKB)/ZTHM_F(:,KKB))*PSFTH(:)+     &
                (0.61*ZTHM_F(:,KKB))*PSFRV(:)

  ! Mass flux at KKB level (updraft triggered if PSFTH>0.)
  IF (LGZ) THEN
    ZSURF(:)=TANH(XGZ*SQRT(XDXHAT(1)*XDYHAT(1))/ZLUP)
  ELSE
    ZSURF(:)=1.
  END IF
  WHERE (ZWTHVSURF(:)>0.)
    PEMF(:,KKB) = XCMF * ZSURF(:) * ZRHO_F(:,KKB) *  &
            ((ZG_O_THVREF(:,KKB))*ZWTHVSURF*ZLUP)**(1./3.)
    PFRAC_UP(:,KKB)=MIN(PEMF(:,KKB)/(SQRT(ZW_UP2(:,KKB))*ZRHO_F(:,KKB)),XFRAC_UP_MAX)
    ZW_UP2(:,KKB)=(PEMF(:,KKB)/(PFRAC_UP(:,KKB)*ZRHO_F(:,KKB)))**2
    GTEST(:)=.TRUE.
  ELSEWHERE
    PEMF(:,KKB) =0.
    GTEST(:)=.FALSE.
  ENDWHERE
ELSE
  GTEST(:)=PEMF(:,KKB+KKL)>0.
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
DO JK=KKB,KKE-KKL,KKL

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
    IF(JK/=KKB) THEN
      ZRC_MIX(:,JK) = ZRC_MIX(:,JK-KKL) ! guess of Rc of mixture
      ZRI_MIX(:,JK) = ZRI_MIX(:,JK-KKL) ! guess of Ri of mixture
    ENDIF
    CALL COMPUTE_ENTR_DETR(JK,KKB,KKE,KKL,GTEST,GTESTLCL,HFRAC_ICE,PFRAC_ICE_UP(:,JK),&
                           PRHODREF(:,JK),ZPRES_F(:,JK),ZPRES_F(:,JK+KKL),&
                           PZZ(:,:),PDZZ(:,:),ZTHVM(:,:),  &
                           PTHLM(:,:),PRTM(:,:),ZW_UP2(:,:),ZTH_UP(:,JK),   &
                           PTHL_UP(:,JK),PRT_UP(:,JK),ZLUP(:),         &
                           PRC_UP(:,JK),PRI_UP(:,JK),PTHV_UP(:,JK),&
                           PRSAT_UP(:,JK),ZRC_MIX(:,JK),ZRI_MIX(:,JK),                 &
                           PENTR(:,JK),PDETR(:,JK),ZENTR_CLD(:,JK),ZDETR_CLD(:,JK),&
                           ZBUO_INTEG_DRY(:,JK), ZBUO_INTEG_CLD(:,JK), &
                           ZPART_DRY(:)   )
    PBUO_INTEG(:,JK)=ZBUO_INTEG_DRY(:,JK)+ZBUO_INTEG_CLD(:,JK)

    IF (JK==KKB) THEN
       PDETR(:,JK)=0.
       ZDETR_CLD(:,JK)=0.
    ENDIF   
 
!       Computation of updraft characteristics at level JK+KKL
    WHERE(GTEST)
      ZMIX1(:)=0.5*(PZZ(:,JK+KKL)-PZZ(:,JK))*(PENTR(:,JK)-PDETR(:,JK))
      PEMF(:,JK+KKL)=PEMF(:,JK)*EXP(2*ZMIX1(:))
    ENDWHERE
  ELSE
    GTEST(:) = (PEMF(:,JK+KKL)>0.)
  END IF 
  

! stop the updraft if MF becomes negative
  WHERE (GTEST.AND.(PEMF(:,JK+KKL)<=0.))
    PEMF(:,JK+KKL)=0.
    KKCTL(:) = JK+KKL
    GTEST(:)=.FALSE.
    PFRAC_ICE_UP(:,JK+KKL)=PFRAC_ICE_UP(:,JK)
    PRSAT_UP(:,JK+KKL)=PRSAT_UP(:,JK)
  ENDWHERE


! If the updraft did not stop, compute cons updraft characteritics at jk+KKL
! WHERE(GTEST)     
  DO JLOOP=1,SIZE(GTEST) 
   IF (GTEST(JLOOP) ) THEN
    ZMIX2(JLOOP) = (PZZ(JLOOP,JK+KKL)-PZZ(JLOOP,JK))*PENTR(JLOOP,JK) !&
    ZMIX3_CLD(JLOOP) = (PZZ(JLOOP,JK+KKL)-PZZ(JLOOP,JK))*(1.-ZPART_DRY(JLOOP))*ZDETR_CLD(JLOOP,JK) !&                   
    ZMIX2_CLD(JLOOP) = (PZZ(JLOOP,JK+KKL)-PZZ(JLOOP,JK))*(1.-ZPART_DRY(JLOOP))*ZENTR_CLD(JLOOP,JK)
                
    !PTHL_UP(JLOOP,JK+KKL)=(PTHL_UP(JLOOP,JK)*(1.-0.5*ZMIX2(JLOOP)) + PTHLM(JLOOP,JK)*ZMIX2(JLOOP)) &
    !                      /(1.+0.5*ZMIX2(JLOOP))   
    !PRT_UP(JLOOP,JK+KKL) =(PRT_UP (JLOOP,JK)*(1.-0.5*ZMIX2(JLOOP)) + PRTM(JLOOP,JK)*ZMIX2(JLOOP))  &
    !                      /(1.+0.5*ZMIX2(JLOOP))

    PTHL_UP(JLOOP,JK+KKL)=PTHL_UP(JLOOP,JK)*EXP(-ZMIX2(JLOOP)) + PTHLM(JLOOP,JK)*(1-EXP(-ZMIX2(JLOOP))) 
    PRT_UP(JLOOP,JK+KKL) =PRT_UP (JLOOP,JK)*EXP(-ZMIX2(JLOOP)) +  PRTM(JLOOP,JK)*(1-EXP(-ZMIX2(JLOOP)))  

   END IF
  END DO
! ENDWHERE
  

  IF(OMIXUV) THEN
    IF(JK/=KKB) THEN
      WHERE(GTEST)
        PU_UP(:,JK+KKL) = (PU_UP (:,JK)*(1-0.5*ZMIX2(:)) + PUM(:,JK)*ZMIX2(:)+ &
                          0.5*XPRES_UV*(PZZ(:,JK+KKL)-PZZ(:,JK))*&
                          ((PUM(:,JK+KKL)-PUM(:,JK))/PDZZ(:,JK+KKL)+&
                           (PUM(:,JK)-PUM(:,JK-KKL))/PDZZ(:,JK))        )   &
                          /(1+0.5*ZMIX2(:))
        PV_UP(:,JK+KKL) = (PV_UP (:,JK)*(1-0.5*ZMIX2(:)) + PVM(:,JK)*ZMIX2(:)+ &
                          0.5*XPRES_UV*(PZZ(:,JK+KKL)-PZZ(:,JK))*&
                          ((PVM(:,JK+KKL)-PVM(:,JK))/PDZZ(:,JK+KKL)+&
                           (PVM(:,JK)-PVM(:,JK-KKL))/PDZZ(:,JK))    )   &
                          /(1+0.5*ZMIX2(:))
      ENDWHERE
    ELSE
      WHERE(GTEST)
        PU_UP(:,JK+KKL) = (PU_UP (:,JK)*(1-0.5*ZMIX2(:)) + PUM(:,JK)*ZMIX2(:)+ &
                          0.5*XPRES_UV*(PZZ(:,JK+KKL)-PZZ(:,JK))*&
                          ((PUM(:,JK+KKL)-PUM(:,JK))/PDZZ(:,JK+KKL))        )   &
                          /(1+0.5*ZMIX2(:))
        PV_UP(:,JK+KKL) = (PV_UP (:,JK)*(1-0.5*ZMIX2(:)) + PVM(:,JK)*ZMIX2(:)+ &
                          0.5*XPRES_UV*(PZZ(:,JK+KKL)-PZZ(:,JK))*&
                          ((PVM(:,JK+KKL)-PVM(:,JK))/PDZZ(:,JK+KKL))    )   &
                          /(1+0.5*ZMIX2(:))
      ENDWHERE

    ENDIF
  ENDIF
  DO JSV=1,ISV 
     IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
      WHERE(GTEST) 
           PSV_UP(:,JK+KKL,JSV) = (PSV_UP (:,JK,JSV)*(1-0.5*ZMIX2(:)) + &
                        PSVM(:,JK,JSV)*ZMIX2(:))  /(1+0.5*ZMIX2(:))
      ENDWHERE                        
  END DO  
  
 IF (OENTR_DETR) THEN

! Compute non cons. var. at level JK+KKL
  ZRC_UP(:)=PRC_UP(:,JK) ! guess = level just below
  ZRI_UP(:)=PRI_UP(:,JK) ! guess = level just below
  CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE_UP(:,JK+KKL),ZPRES_F(:,JK+KKL), &
          PTHL_UP(:,JK+KKL),PRT_UP(:,JK+KKL),ZTH_UP(:,JK+KKL),              &
          ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:))
  WHERE(GTEST)
    PRC_UP(:,JK+KKL)=ZRC_UP(:)
    PRV_UP(:,JK+KKL)=ZRV_UP(:)
    PRI_UP(:,JK+KKL)=ZRI_UP(:)
    PRSAT_UP(:,JK+KKL) = ZRSATW(:)*(1-PFRAC_ICE_UP(:,JK+KKL)) + ZRSATI(:)*PFRAC_ICE_UP(:,JK+KKL)
  ENDWHERE
  

! Compute the updraft theta_v, buoyancy and w**2 for level JK+KKL
  WHERE(GTEST)
    PTHV_UP(:,JK+KKL) = ZTH_UP(:,JK+KKL)*((1+ZRVORD*PRV_UP(:,JK+KKL))/(1+PRT_UP(:,JK+KKL)))
    WHERE (ZBUO_INTEG_DRY(:,JK)>0.)
      ZW_UP2(:,JK+KKL)  = ZW_UP2(:,JK) + 2.*(XABUO-XBENTR*XENTR_DRY)* ZBUO_INTEG_DRY(:,JK)
    ELSEWHERE
      ZW_UP2(:,JK+KKL)  = ZW_UP2(:,JK) + 2.*XABUO* ZBUO_INTEG_DRY(:,JK)
    ENDWHERE
    ZW_UP2(:,JK+KKL)  = ZW_UP2(:,JK+KKL)*(1.-(XBDETR*ZMIX3_CLD(:)+XBENTR*ZMIX2_CLD(:)))&
            /(1.+(XBDETR*ZMIX3_CLD(:)+XBENTR*ZMIX2_CLD(:))) &
            +2.*(XABUO)*ZBUO_INTEG_CLD(:,JK)/(1.+(XBDETR*ZMIX3_CLD(:)+XBENTR*ZMIX2_CLD(:)))
 ENDWHERE


! Test if the updraft has reach the ETL
  GTESTETL(:)=.FALSE.
  WHERE (GTEST.AND.(PBUO_INTEG(:,JK)<=0.))
      KKETL(:) = JK+KKL
      GTESTETL(:)=.TRUE.
  ENDWHERE

! Test is we have reached the top of the updraft
  WHERE (GTEST.AND.((ZW_UP2(:,JK+KKL)<=0.).OR.(PEMF(:,JK+KKL)<=0.)))
      ZW_UP2(:,JK+KKL)=0.
      PEMF(:,JK+KKL)=0.
      GTEST(:)=.FALSE.
      PTHL_UP(:,JK+KKL)=ZTHLM_F(:,JK+KKL)
      PRT_UP(:,JK+KKL)=ZRTM_F(:,JK+KKL)
      PRC_UP(:,JK+KKL)=0.
      PRI_UP(:,JK+KKL)=0.
      PRV_UP(:,JK+KKL)=0.
      PTHV_UP(:,JK+KKL)=ZTHVM_F(:,JK+KKL)
      PFRAC_UP(:,JK+KKL)=0.
      KKCTL(:)=JK+KKL
  ENDWHERE
 
! compute frac_up at JK+KKL
  WHERE (GTEST)
      PFRAC_UP(:,JK+KKL)=PEMF(:,JK+KKL)/(SQRT(ZW_UP2(:,JK+KKL))*ZRHO_F(:,JK+KKL))
  ENDWHERE

! Updraft fraction must be smaller than XFRAC_UP_MAX
  WHERE (GTEST)
      PFRAC_UP(:,JK+KKL)=MIN(XFRAC_UP_MAX,PFRAC_UP(:,JK+KKL))
  ENDWHERE


! When cloudy and non-buoyant, updraft fraction must decrease
  
  WHERE ((GTEST.AND.GTESTETL).AND.GTESTLCL)
    PFRAC_UP(:,JK+KKL)=MIN(PFRAC_UP(:,JK+KKL),PFRAC_UP(:,JK))
  ENDWHERE

! Mass flux is updated with the new updraft fraction
  
  IF (OENTR_DETR) PEMF(:,JK+KKL)=PFRAC_UP(:,JK+KKL)*SQRT(ZW_UP2(:,JK+KKL))*ZRHO_F(:,JK+KKL)

 END IF

ENDDO

IF(OENTR_DETR) THEN

  PW_UP(:,:)=SQRT(ZW_UP2(:,:))

  PEMF(:,KKB) =0.

! Limits the shallow convection scheme when cloud heigth is higher than 3000m.
! To do this, mass flux is multiplied by a coefficient decreasing linearly
! from 1 (for clouds of ZDEPTH_MAX1 m of depth) to 0 (for clouds of ZDEPTH_MAX2 m of depth).
! This way, all MF fluxes are diminished by this amount.
! Diagnosed cloud fraction is also multiplied by the same coefficient.
!
  DO JI=1,SIZE(PTHM,1) 
     PDEPTH(JI) = MAX(0., PZZ(JI,KKCTL(JI)) -  PZZ(JI,KKLCL(JI)) )
  END DO

  GWORK1(:)= (GTESTLCL(:) .AND. (PDEPTH(:) > ZDEPTH_MAX1) )
  GWORK2(:,:) = SPREAD( GWORK1(:), DIM=2, NCOPIES=MAX(KKU,KKA) )
  ZCOEF(:,:) = SPREAD( (1.-(PDEPTH(:)-ZDEPTH_MAX1)/(ZDEPTH_MAX2-ZDEPTH_MAX1)), DIM=2, NCOPIES=SIZE(ZCOEF,2))
  ZCOEF=MIN(MAX(ZCOEF,0.),1.)
  WHERE (GWORK2) 
    PEMF(:,:)     = PEMF(:,:)     * ZCOEF(:,:)
    PFRAC_UP(:,:) = PFRAC_UP(:,:) * ZCOEF(:,:)
  ENDWHERE
ENDIF   
END SUBROUTINE COMPUTE_UPDRAFT
