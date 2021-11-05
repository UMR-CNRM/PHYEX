!MNH_LIC Copyright 2012-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
     MODULE MODI_COMPUTE_UPDRAFT_RHCJ10
!    ###########################
!
INTERFACE
!
!     #################################################################
      SUBROUTINE COMPUTE_UPDRAFT_RHCJ10(KKA,KKB,KKE,KKU,KKL, HFRAC_ICE,  &
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


END SUBROUTINE COMPUTE_UPDRAFT_RHCJ10

END INTERFACE
!
END MODULE MODI_COMPUTE_UPDRAFT_RHCJ10
!
SUBROUTINE COMPUTE_UPDRAFT_RHCJ10(KKA,KKB,KKE,KKU,KKL,HFRAC_ICE,       &
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
!!****  *COMPUTE_UPDRAF_RHCJ10* - calculates caracteristics of the updraft 
!!                         
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to build the updraft following Rio et al (2010)
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
!!       Rio et al (2010) (Boundary Layer Meteorol 135:469-483)
!!
!!     AUTHOR
!!     ------
!!     Y. Bouteloup (2012)
!!     R. Honert Janv 2013 ==> corection of some bugs
!!     Q.Rodier  01/2019 : support RM17 mixing length 
!! --------------------------------------------------------------------------

! WARNING ==>  This updraft is not yet ready to use scalar variables 

!*      0. DECLARATIONS
!          ------------

USE MODD_CST
USE MODD_PARAM_MFSHALL_n
USE MODD_TURB_n, ONLY : CTURBLEN
USE MODI_COMPUTE_ENTR_DETR
USE MODI_TH_R_FROM_THL_RT_1D
USE MODI_SHUMAN_MF

USE MODI_COMPUTE_BL89_ML


IMPLICIT NONE

!*                    1.1  Declaration of Arguments


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

REAL, DIMENSION(:,:),   INTENT(IN) ::  PPABSM     ! Pressure at t-dt
REAL, DIMENSION(:,:),   INTENT(IN) ::  PRHODREF   ! dry density of the
                                                  ! reference state
REAL, DIMENSION(:,:),   INTENT(IN) ::  PUM        ! u mean wind
REAL, DIMENSION(:,:),   INTENT(IN) ::  PVM        ! v mean wind
REAL, DIMENSION(:,:),   INTENT(IN) ::  PTKEM      ! TKE at t-dt
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHM           ! pot. temp. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PRVM           ! vapor mixing ratio at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHLM,PRTM     ! cons. var. at t-dt

REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSVM           ! scalar var. at t-dt

REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PTHL_UP,PRT_UP   ! updraft properties
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PU_UP, PV_UP     ! updraft wind components
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PRV_UP,PRC_UP    ! updraft rv, rc
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PRI_UP           ! updraft ri
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PTHV_UP          ! updraft THv
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PW_UP,PFRAC_UP   ! updraft w, fraction
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PFRAC_ICE_UP     ! liquid/solid fraction in updraft
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PRSAT_UP         ! Rsat

REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PSV_UP           ! updraft scalar var. 
                                         
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PEMF,PDETR,PENTR ! Mass_flux,
                                                          ! detrainment,entrainment
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PBUO_INTEG       ! Integrated Buoyancy 
INTEGER, DIMENSION(:),  INTENT(INOUT)::  KKLCL,KKETL,KKCTL! LCL, ETL, CTL                                     
REAL, DIMENSION(:),     INTENT(OUT)   :: PDEPTH           ! Deepness of cloud
!                       1.2  Declaration of local variables

! Mean environment variables at t-dt at flux point
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::    ZTHM_F,ZRVM_F    ! Theta,rv of
                                                                  ! updraft environnement
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZRTM_F, ZTHLM_F, ZTKEM_F      ! rt, thetal,TKE,pressure,
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZUM_F,ZVM_F,ZRHO_F            ! density,momentum
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZPRES_F,ZTHVM_F,ZTHVM         ! interpolated at the flux point
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZG_O_THVREF                   ! g*ThetaV ref
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZW_UP2                        ! w**2  of the updraft

REAL, DIMENSION(SIZE(PSVM,1),SIZE(PTHM,2),SIZE(PSVM,3)) :: ZSVM_F ! scalar variables 
                        

                        
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZTH_UP                        ! updraft THETA 
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZBUO                          ! Buoyancy 

REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::  ZCOEF  ! diminution coefficient for too high clouds 
                        
REAL, DIMENSION(SIZE(PSFTH,1) )            ::  ZWTHVSURF  ! Surface w'thetav'

REAL  :: ZRDORV       ! RD/RV
REAL  :: ZRVORD       ! RV/RD


REAL, DIMENSION(SIZE(PTHM,1)) :: ZMIX1,ZMIX2,ZMIX3

REAL, DIMENSION(SIZE(PTHM,1)) :: ZLUP         ! Upward Mixing length from the ground


INTEGER  :: ISV                ! Number of scalar variables                               
INTEGER  :: IKU,IIJU           ! array size in k
INTEGER  :: JK,JI,JJ,JSV          ! loop counters

LOGICAL, DIMENSION(SIZE(PTHM,1)) ::  GTEST,GTESTLCL,GTESTETL
                               ! Test if the ascent continue, if LCL or ETL is reached
LOGICAL                          ::  GLMIX 
                               ! To choose upward or downward mixing length
LOGICAL, DIMENSION(SIZE(PTHM,1))              :: GWORK1
LOGICAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: GWORK2

INTEGER  :: ITEST

REAL, DIMENSION(SIZE(PTHM,1)) :: ZRC_UP, ZRI_UP, ZRV_UP, ZRSATW, ZRSATI

REAL,    DIMENSION(SIZE(PTHM,1)) :: ZPHI
REAL,    DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZZDZ,ZZZ

REAL, DIMENSION(SIZE(PTHM,1))              ::  ZTEST,ZDZ,ZWUP_MEAN    ! 
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZCOE,ZWCOE,ZBUCOE
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZDETR_BUO, ZDETR_RT
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZW_MAX               ! w**2  max of the updraft
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZZTOP                ! Top of the updraft
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZBETA1

REAL  :: ZDEPTH_MAX1, ZDEPTH_MAX2 ! control auto-extinction process

REAL  :: ZTMAX,ZRMAX, ZEPS  ! control value

REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZSHEAR,ZDUDZ,ZDVDZ ! vertical wind shear

! Thresholds for the  perturbation of
! theta_l and r_t at the first level of the updraft

ZTMAX=2.0
ZRMAX=1.E-3
ZEPS=1.E-15
!------------------------------------------------------------------------
!                     INITIALISATION

! Initialisation of the constants   
ZRDORV   = XRD / XRV   !=0.622
ZRVORD   = (XRV / XRD) 

! depth are different in compute_updraft (3000. and 4000.) ==> impact is small
ZDEPTH_MAX1=4500. ! clouds with depth infeRIOr to this value are keeped untouched
ZDEPTH_MAX2=5000. ! clouds with depth superior to this value are suppressed


!  Initialisation of ZBETA1 ==> I do not remember why I introduced a KLON array for beta1 !

ZBETA1(:) = XBETA1
   
!                 Local variables, internal domain
! Internal Domain

IKU=SIZE(PTHM,2)
IIJU =SIZE(PTHM,1)
!number of scalar variables
ISV=SIZE(PSVM,3)

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
PRC_UP(:,:)=0.

PW_UP(:,:)=0.
ZTH_UP(:,:)=0.
PFRAC_UP(:,:)=0.
PTHV_UP(:,:)=0.

PBUO_INTEG=0.
ZBUO      =0.

PRI_UP(:,:)=0.
PFRAC_ICE_UP(:,:)=0.
PRSAT_UP(:,:)=PRVM(:,:) ! should be initialised correctly but is (normaly) not used

! Initialisation of environment variables at t-dt
! variables at flux level
ZTHLM_F(:,:) = MZM_MF(KKA,KKU,KKL,PTHLM(:,:))
ZRTM_F (:,:) = MZM_MF(KKA,KKU,KKL,PRTM(:,:))
ZUM_F  (:,:) = MZM_MF(KKA,KKU,KKL,PUM(:,:))
ZVM_F  (:,:) = MZM_MF(KKA,KKU,KKL,PVM(:,:))
ZTKEM_F(:,:) = MZM_MF(KKA,KKU,KKL,PTKEM(:,:))


! This updraft is not yet ready to use scalar variables 
!DO JSV=1,ISV
!  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
!  ZSVM_F(:,:,JSV) = MZM_MF(KKA,KKU,KKL,PSVM(:,:,JSV))
!END DO

!          Initialisation of updraft characteristics 
PTHL_UP(:,:)=ZTHLM_F(:,:)
PRT_UP(:,:)=ZRTM_F(:,:)
PU_UP(:,:)=ZUM_F(:,:)
PV_UP(:,:)=ZVM_F(:,:)
PSV_UP(:,:,:)=0.
! This updraft is not yet ready to use scalar variables 
!IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) then
!    PSV_UP(:,:,:)=ZSVM_F(:,:,:)
!ENDIF

! Computation or initialisation of updraft characteristics at the KKB level
! thetal_up,rt_up,thetaV_up, wï¿½,Buoyancy term and mass flux (PEMF)

!PTHL_UP(:,KKB)= ZTHLM_F(:,KKB)+MAX(0.,MIN(ZTMAX,(PSFTH(:)/SQRT(ZTKEM_F(:,KKB)))*XALP_PERT))
!PRT_UP(:,KKB) = ZRTM_F(:,KKB)+MAX(0.,MIN(ZRMAX,(PSFRV(:)/SQRT(ZTKEM_F(:,KKB)))*XALP_PERT)) 
PTHL_UP(:,KKB)= ZTHLM_F(:,KKB)
PRT_UP(:,KKB) = ZRTM_F(:,KKB)

ZTHM_F (:,:) = MZM_MF(KKA,KKU,KKL,PTHM (:,:))
ZPRES_F(:,:) = MZM_MF(KKA,KKU,KKL,PPABSM(:,:))
ZRHO_F (:,:) = MZM_MF(KKA,KKU,KKL,PRHODREF(:,:))
ZRVM_F (:,:) = MZM_MF(KKA,KKU,KKL,PRVM(:,:))

! thetav at mass and flux levels 
ZTHVM_F(:,:)=ZTHM_F(:,:)*((1.+ZRVORD*ZRVM_F(:,:))/(1.+ZRTM_F(:,:)))
ZTHVM(:,:)=PTHM(:,:)*((1.+ZRVORD*PRVM(:,:))/(1.+PRTM(:,:)))

PTHV_UP(:,:)= ZTHVM_F(:,:)
PRV_UP (:,:)= ZRVM_F (:,:)

ZW_UP2(:,:)=ZEPS
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

!Tout est commente pour tester dans un premier temps la séparation en deux de la 
!  boucle verticale, une pour w et une pour PEMF                                                            

ZG_O_THVREF=XG/ZTHVM_F

! Calcul de la fermeture de Julien Pergaut comme limite max de PHY

DO JK=KKB,KKE-KKL,KKL   !  Vertical loop
  ZZDZ(:,JK)   = MAX(ZEPS,PZZ(:,JK+KKL)-PZZ(:,JK))  ! <== Delta Z between two flux level
  ZZZ(:,JK)    = 0.5*(PZZ(:,JK+KKL)+PZZ(:,JK))     ! <== Hight of mass levels
ENDDO

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
CALL COMPUTE_BL89_ML(KKA,KKB,KKE,KKU,KKL,PDZZ,ZTKEM_F(:,KKB),ZG_O_THVREF(:,KKB), &
                       ZTHVM_F,KKB,GLMIX,.TRUE.,ZSHEAR,ZLUP)
ZLUP(:)=MAX(ZLUP(:),1.E-10)

! Compute Buoyancy flux at the ground
ZWTHVSURF(:) = (ZTHVM_F(:,KKB)/ZTHM_F(:,KKB))*PSFTH(:)+     &
                (0.61*ZTHM_F(:,KKB))*PSFRV(:)

! Mass flux at KKB level (updraft triggered if PSFTH>0.)
WHERE (ZWTHVSURF(:)>0.010)  ! <==  Not 0 Important to have stratocumulus !!!!!
  PEMF(:,KKB) = XCMF * ZRHO_F(:,KKB) * ((ZG_O_THVREF(:,KKB))*ZWTHVSURF*ZLUP)**(1./3.)
  PFRAC_UP(:,KKB)=MIN(PEMF(:,KKB)/(SQRT(ZW_UP2(:,KKB))*ZRHO_F(:,KKB)),XFRAC_UP_MAX)
  ZW_UP2(:,KKB)=(PEMF(:,KKB)/(PFRAC_UP(:,KKB)*ZRHO_F(:,KKB)))**2
  GTEST(:)=.TRUE.
ELSEWHERE
  PEMF(:,KKB) =0.
  GTEST(:)=.FALSE.
ENDWHERE

  
!--------------------------------------------------------------------------

!                        3. Vertical ascending loop
!                           -----------------------
!
! If GTEST = T the updraft starts from the KKB level and stops when GTEST becomes F
!
!
GTESTLCL(:)=.FALSE.
GTESTETL(:)=.FALSE.


!       Loop on vertical level to compute W

ZW_MAX(:)      = 0.
ZZTOP(:)       = 0.
ZPHI(:) = 0.

DO JK=KKB,KKE-KKL,KKL

! IF the updraft top is reached for all column, stop the loop on levels

  ITEST=COUNT(GTEST)
!  IF (ITEST==0) CYCLE  !  <== I do not remember why I removed this ... 

!       Computation of entrainment and detrainment with KF90
!       parameterization in clouds and LR01 in subcloud layer
  
 
! to find the LCL (check if JK is LCL or not)

  WHERE ((PRC_UP(:,JK)+PRI_UP(:,JK)>0.).AND.(.NOT.(GTESTLCL)))
      KKLCL(:) = JK           
      GTESTLCL(:)=.TRUE.
  ENDWHERE


! COMPUTE PENTR and PDETR at mass level JK

    
!  Buoyancy is computed on "flux" levels where updraft variables are known

  ! Compute theta_v of updraft at flux level JK    
    
    ZRC_UP(:)   =PRC_UP(:,JK) ! guess
    ZRI_UP(:)   =PRI_UP(:,JK) ! guess 
    ZRV_UP(:)   =PRV_UP(:,JK)
    CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE_UP(:,JK),&
               PPABSM(:,JK),PTHL_UP(:,JK),PRT_UP(:,JK),&
               ZTH_UP(:,JK),ZRV_UP,ZRC_UP,ZRI_UP,ZRSATW(:),ZRSATI(:))            
    
    WHERE (GTEST)
      PTHV_UP   (:,JK) = ZTH_UP(:,JK)*(1.+ZRVORD*ZRV_UP(:))/(1.+PRT_UP(:,JK))
      ZBUO      (:,JK) = ZG_O_THVREF(:,JK)*(PTHV_UP(:,JK) - ZTHVM_F(:,JK))    
      PBUO_INTEG(:,JK) = ZBUO(:,JK)*(PZZ(:,JK+KKL)-PZZ(:,JK))
      
      ZDZ(:)   = MAX(ZEPS,PZZ(:,JK+KKL)-PZZ(:,JK))
      ZTEST(:) = XA1*ZBUO(:,JK) -  XB*ZW_UP2(:,JK)

      ZCOE(:)      = ZDZ(:)
      WHERE (ZTEST(:)>0.)
        ZCOE(:)    = ZDZ(:)/(1.+ ZBETA1(:))
      ENDWHERE   

!  Convective Vertical speed computation

      ZWCOE(:)         = (1.-XB*ZCOE(:))/(1.+XB*ZCOE(:))
      ZBUCOE(:)        =  2.*ZCOE(:)/(1.+XB*ZCOE(:))

! Second Rachel bug correction (XA1 has been forgotten ... not yet tested ...)
!      ZW_UP2(:,JK+KKL) = MAX(ZEPS,ZW_UP2(:,JK)*ZWCOE(:) + ZBUO(:,JK)*ZBUCOE(:) ) 
      ZW_UP2(:,JK+KKL) = MAX(ZEPS,ZW_UP2(:,JK)*ZWCOE(:) + XA1*ZBUO(:,JK)*ZBUCOE(:) )  
      ZW_MAX(:) = MAX(ZW_MAX(:), SQRT(ZW_UP2(:,JK+KKL)))
      ZWUP_MEAN(:)     = MAX(ZEPS,0.5*(ZW_UP2(:,JK+KKL)+ZW_UP2(:,JK)))
 
!  Entrainement and detrainement

! First Rachel bug correction (Parenthesis around 1+beta1 ==> impact is small)   
     PENTR(:,JK)  = MAX(0.,(ZBETA1(:)/(1.+ZBETA1(:)))*(XA1*ZBUO(:,JK)/ZWUP_MEAN(:)-XB))
     ZDETR_BUO(:) = MAX(0., -(ZBETA1(:)/(1.+ZBETA1(:)))*XA1*ZBUO(:,JK)/ZWUP_MEAN(:))
     ZDETR_RT(:)  = XC*SQRT(MAX(0.,(PRT_UP(:,JK) - ZRTM_F(:,JK))) / MAX(ZEPS,ZRTM_F(:,JK)) / ZWUP_MEAN(:))
     PDETR(:,JK)  = ZDETR_RT(:)+ZDETR_BUO(:)
   
! If the updraft did not stop, compute cons updraft characteritics at jk+1

     ZZTOP(:) = MAX(ZZTOP(:),PZZ(:,JK+KKL))
     ZMIX2(:) = (PZZ(:,JK+KKL)-PZZ(:,JK))*PENTR(:,JK) !&
     ZMIX3(:) = (PZZ(:,JK+KKL)-PZZ(:,JK))*PDETR(:,JK) !&                   
                
     PTHL_UP(:,JK+KKL)=(PTHL_UP(:,JK)*(1.-0.5*ZMIX2(:)) + PTHLM(:,JK)*ZMIX2(:)) &
                          /(1.+0.5*ZMIX2(:))   
     PRT_UP(:,JK+KKL) =(PRT_UP (:,JK)*(1.-0.5*ZMIX2(:)) + PRTM(:,JK)*ZMIX2(:))  &
                          /(1.+0.5*ZMIX2(:))
  ENDWHERE  ! GTEST
  

  IF(OMIXUV) THEN
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
  ENDIF
  
! This updraft is not yet ready to use scalar variables  
!  DO JSV=1,ISV 
!     IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
!      WHERE(GTEST) 
!           PSV_UP(:,JK+KKL,JSV) = (PSV_UP (:,JK,JSV)*(1-0.5*ZMIX2(:)) + &
!                        PSVM(:,JK,JSV)*ZMIX2(:))  /(1+0.5*ZMIX2(:))
!      ENDWHERE                        
!  ENDDO  
  
  
! Compute non cons. var. at level JK+KKL
  ZRC_UP(:)=PRC_UP(:,JK) ! guess = level just below
  ZRI_UP(:)=PRI_UP(:,JK) ! guess = level just below
  ZRV_UP(:)=PRV_UP(:,JK)
  CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE_UP(:,JK+KKL),ZPRES_F(:,JK+KKL), &
          PTHL_UP(:,JK+KKL),PRT_UP(:,JK+KKL),ZTH_UP(:,JK+KKL),              &
          ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:))
  WHERE(GTEST)
    PRC_UP(:,JK+KKL)=ZRC_UP(:)
    PRV_UP(:,JK+KKL)=ZRV_UP(:)
    PRI_UP(:,JK+KKL)=ZRI_UP(:)
    PRSAT_UP(:,JK+KKL) = ZRSATW(:)*(1-PFRAC_ICE_UP(:,JK+KKL)) + ZRSATI(:)*PFRAC_ICE_UP(:,JK+KKL)
  ENDWHERE
  

! Compute the updraft theta_v, buoyancy and w**2 for level JK+1   
 WHERE(GTEST)
      PTHV_UP(:,JK+KKL) = ZTH_UP(:,JK+KKL)*((1+ZRVORD*PRV_UP(:,JK+KKL))/(1+PRT_UP(:,JK+KKL)))
 ENDWHERE

 WHERE(GTEST) 
    ZMIX1(:)=ZZDZ(:,JK)*(PENTR(:,JK)-PDETR(:,JK))
    PEMF(:,JK+KKL)=PEMF(:,JK)*EXP(ZMIX1(:))

! Updraft fraction must be smaller than XFRAC_UP_MAX
    PFRAC_UP(:,JK+KKL)=PEMF(:,JK+KKL)/(SQRT(ZW_UP2(:,JK+KKL))*ZRHO_F(:,JK+KKL))
    PFRAC_UP(:,JK+KKL)=MIN(XFRAC_UP_MAX,PFRAC_UP(:,JK+KKL))
  ENDWHERE  

! Test if the updraft has reach the ETL
  GTESTETL(:)=.FALSE.
  WHERE (GTEST.AND.(PBUO_INTEG(:,JK)<=0.))
      KKETL(:) = JK+KKL
      GTESTETL(:)=.TRUE.
  ENDWHERE


! Test is we have reached the top of the updraft

  WHERE (GTEST.AND.((ZW_UP2(:,JK+KKL)<=ZEPS).OR.(PEMF(:,JK+KKL)<=ZEPS)))
      ZW_UP2   (:,JK+KKL)=ZEPS
      PEMF     (:,JK+KKL)=0.
      GTEST    (:)       =.FALSE.
      PTHL_UP  (:,JK+KKL)=ZTHLM_F(:,JK+KKL)
      PRT_UP   (:,JK+KKL)=ZRTM_F(:,JK+KKL)
      PRC_UP   (:,JK+KKL)=0.
      PRI_UP   (:,JK+KKL)=0.
      PRV_UP   (:,JK+KKL)=ZRVM_F (:,JK+KKL)
      PTHV_UP  (:,JK+KKL)=ZTHVM_F(:,JK+KKL)
      PFRAC_UP (:,JK+KKL)=0.
      KKCTL    (:)       =JK+KKL
      
  ENDWHERE


ENDDO   ! Fin de la boucle verticale 

PW_UP(:,:)=SQRT(ZW_UP2(:,:))
PEMF(:,KKB) =0.

! Limits the shallow convection scheme when cloud heigth is higher than 3000m.
! To do this, mass flux is multiplied by a coefficient decreasing linearly
! from 1 (for clouds of 3000m of depth) to 0 (for clouds of 4000m of depth).
! This way, all MF fluxes are diminished by this amount.
! Diagnosed cloud fraction is also multiplied by the same coefficient.
!
DO JI=1,SIZE(PTHM,1) 
   PDEPTH(JI) = MAX(0., PZZ(JI,KKCTL(JI)) -  PZZ(JI,KKLCL(JI)) )
END DO

GWORK1(:)= (GTESTLCL(:) .AND. (PDEPTH(:) > ZDEPTH_MAX1) )
GWORK2(:,:) = SPREAD( GWORK1(:), DIM=2, NCOPIES=IKU )
ZCOEF(:,:) = SPREAD( (1.-(PDEPTH(:)-ZDEPTH_MAX1)/(ZDEPTH_MAX2-ZDEPTH_MAX1)), DIM=2, NCOPIES=IKU)
ZCOEF=MIN(MAX(ZCOEF,0.),1.)
WHERE (GWORK2) 
   PEMF(:,:)     = PEMF(:,:)     * ZCOEF(:,:) 
   PFRAC_UP(:,:) = PFRAC_UP(:,:) * ZCOEF(:,:) 
ENDWHERE

END SUBROUTINE COMPUTE_UPDRAFT_RHCJ10


