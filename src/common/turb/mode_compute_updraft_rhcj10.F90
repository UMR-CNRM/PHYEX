!MNH_LIC Copyright 2012-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
     MODULE MODE_COMPUTE_UPDRAFT_RHCJ10
!    ###########################
!
IMPLICIT NONE
CONTAINS
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
!!****  *COMPUTE_UPDRAFT_RHCJ10* - calculates caracteristics of the updraft 
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
!!     R. El Khatib 15-Oct-2014 Optimization
!!     Q.Rodier  01/2019 : support RM17 mixing length
!! --------------------------------------------------------------------------

! WARNING ==>  This updraft is not yet ready to use scalar variables 

!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_PARAM_MFSHALL_n
USE MODD_TURB_n, ONLY : CTURBLEN
USE MODE_TH_R_FROM_THL_RT_1D, ONLY: TH_R_FROM_THL_RT_1D
USE MODI_SHUMAN_MF, ONLY: MZF_MF, MZM_MF, GZ_M_W_MF

USE MODE_COMPUTE_BL89_ML, ONLY: COMPUTE_BL89_ML
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK


IMPLICIT NONE

!*                    1.1  Declaration of Arguments
!
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
CHARACTER(LEN=1),       INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
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
!REAL, DIMENSION(:,:),   INTENT(IN)   ::  PEXNM       ! Exner function at t-dt
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
!
! Mean environment variables at t-dt at flux point
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::    ZTHM_F,ZRVM_F    ! Theta,rv of
                                                                  ! updraft environnement
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZRTM_F, ZTHLM_F, ZTKEM_F      ! rt, thetal,TKE,pressure,
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZUM_F,ZVM_F,ZRHO_F            ! density,momentum
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZPRES_F,ZTHVM_F               ! interpolated at the flux point
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZG_O_THVREF                   ! g*ThetaV ref
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZW_UP2                        ! w**2  of the updraft

REAL, DIMENSION(SIZE(PSVM,1),SIZE(PTHM,2),SIZE(PSVM,3)) :: ZSVM_F ! scalar variables 
                        

                        
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZTH_UP                        ! updraft THETA 
!REAL, DIMENSION(SIZE(PTHM,1))              :: ZT_UP                         ! updraft T
!REAL, DIMENSION(SIZE(PTHM,1))              :: ZLVOCPEXN                     ! updraft L
!REAL, DIMENSION(SIZE(PTHM,1))              :: ZCP                           ! updraft cp
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZBUO                          ! Buoyancy 
!REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZTHS_UP,ZTHSM

REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::  ZCOEF  ! diminution coefficient for too high clouds 
                        
REAL                                       ::  ZWTHVSURF  ! Surface w'thetav'

REAL  :: ZRVORD       ! RV/RD


REAL, DIMENSION(SIZE(PTHM,1)) :: ZMIX1,ZMIX2

REAL, DIMENSION(SIZE(PTHM,1)) :: ZLUP         ! Upward Mixing length from the ground

INTEGER  :: ISV                ! Number of scalar variables                               
INTEGER  :: IKU,IIJU           ! array size in k
INTEGER  :: JK,JI,JSV          ! loop counters

LOGICAL, DIMENSION(SIZE(PTHM,1)) ::  GTEST,GTESTLCL
                               ! Test if the ascent continue, if LCL or ETL is reached
LOGICAL                          ::  GLMIX 
                               ! To choose upward or downward mixing length
LOGICAL, DIMENSION(SIZE(PTHM,1))              :: GWORK1
LOGICAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: GWORK2

INTEGER  :: ITEST

REAL, DIMENSION(SIZE(PTHM,1)) :: ZRC_UP, ZRI_UP, ZRV_UP, ZRSATW, ZRSATI

REAL,    DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZZDZ

REAL, DIMENSION(SIZE(PTHM,1))              ::  ZTEST,ZDZ,ZWUP_MEAN    ! 
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZCOE,ZWCOE,ZBUCOE
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZDETR_BUO, ZDETR_RT
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZW_MAX               ! w**2  max of the updraft
REAL, DIMENSION(SIZE(PTHM,1))              ::  ZZTOP                ! Top of the updraft
!REAL, DIMENSION(SIZE(PTHM,1))              ::  ZQTM,ZQT_UP

REAL  :: ZDEPTH_MAX1, ZDEPTH_MAX2 ! control auto-extinction process

REAL  :: ZTMAX,ZRMAX, ZEPS  ! control value

REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZSHEAR,ZDUDZ,ZDVDZ ! vertical wind shear
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAFT_RHCJ10',0,ZHOOK_HANDLE)

! Thresholds for the  perturbation of
! theta_l and r_t at the first level of the updraft

ZTMAX=2.0
ZRMAX=1.E-3
ZEPS=1.E-15
!------------------------------------------------------------------------
!                     INITIALISATION

! Initialisation of the constants   
ZRVORD   = (XRV / XRD) 

! depth are different in compute_updraft (3000. and 4000.) ==> impact is small
ZDEPTH_MAX1=4500. ! clouds with depth infeRIOr to this value are keeped untouched
ZDEPTH_MAX2=5000. ! clouds with depth superior to this value are suppressed


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

!no ice cloud coded yet
PRI_UP(:,:)=0.
PFRAC_ICE_UP(:,:)=0.
PRSAT_UP(:,:)=PRVM(:,:) ! should be initialised correctly but is (normaly) not used

! Initialisation of environment variables at t-dt

! variables at flux level
ZTHLM_F(:,:) = MZM_MF(PTHLM(:,:), KKA, KKU, KKL)
ZRTM_F (:,:) = MZM_MF(PRTM(:,:), KKA, KKU, KKL)
ZUM_F  (:,:) = MZM_MF(PUM(:,:), KKA, KKU, KKL)
ZVM_F  (:,:) = MZM_MF(PVM(:,:), KKA, KKU, KKL)
ZTKEM_F(:,:) = MZM_MF(PTKEM(:,:), KKA, KKU, KKL)

! This updraft is not yet ready to use scalar variables
!DO JSV=1,ISV
!  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
! *** SR merge AROME/Meso-nh: following two lines come from the AROME version
!   ZSVM_F(:,KKB:IKU,JSV) = 0.5*(PSVM(:,KKB:IKU,JSV)+PSVM(:,1:IKU-1,JSV))
!   ZSVM_F(:,1,JSV)       = ZSVM_F(:,KKB,JSV) 
! *** the following single line comes from the Meso-NH version
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
! thetal_up,rt_up,thetaV_up, w,Buoyancy term and mass flux (PEMF)

DO JI=1,IIJU
  !PTHL_UP(JI,KKB)= ZTHLM_F(JI,KKB)+MAX(0.,MIN(ZTMAX,(PSFTH(JI)/SQRT(ZTKEM_F(JI,KKB)))*XALP_PERT))
  !PRT_UP(JI,KKB) = ZRTM_F(JI,KKB)+MAX(0.,MIN(ZRMAX,(PSFRV(JI)/SQRT(ZTKEM_F(JI,KKB)))*XALP_PERT)) 
  PTHL_UP(JI,KKB)= ZTHLM_F(JI,KKB)
  PRT_UP(JI,KKB) = ZRTM_F(JI,KKB)
  !ZQT_UP(JI) = PRT_UP(JI,KKB)/(1.+PRT_UP(JI,KKB))
  !ZTHS_UP(JI,KKB)=PTHL_UP(JI,KKB)*(1.+XLAMBDA_MF*ZQT_UP(JI))
ENDDO

ZTHM_F (:,:) = MZM_MF(PTHM (:,:), KKA, KKU, KKL)
ZPRES_F(:,:) = MZM_MF(PPABSM(:,:), KKA, KKU, KKL)
ZRHO_F (:,:) = MZM_MF(PRHODREF(:,:), KKA, KKU, KKL)
ZRVM_F (:,:) = MZM_MF(PRVM(:,:), KKA, KKU, KKL)

! thetav at mass and flux levels 
DO JK=1,IKU
  DO JI=1,IIJU
    ZTHVM_F(JI,JK)=ZTHM_F(JI,JK)*((1.+ZRVORD*ZRVM_F(JI,JK))/(1.+ZRTM_F(JI,JK)))
  ENDDO
ENDDO

PTHV_UP(:,:)= ZTHVM_F(:,:)
PRV_UP (:,:)= ZRVM_F (:,:)

ZW_UP2(:,:)=ZEPS
!ZW_UP2(:,KKB) = MAX(0.0001,(3./6.)*ZTKEM_F(:,KKB))
ZW_UP2(:,KKB) = MAX(0.0001,(2./3.)*ZTKEM_F(:,KKB))

! Computation of non conservative variable for the KKB level of the updraft
! (all or nothing ajustement)

PRC_UP(:,KKB)=0.
PRI_UP(:,KKB)=0.
CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE_UP(:,KKB),ZPRES_F(:,KKB), &
             PTHL_UP(:,KKB),PRT_UP(:,KKB),ZTH_UP(:,KKB), &
             PRV_UP(:,KKB),PRC_UP(:,KKB),PRI_UP(:,KKB),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.)

DO JI=1,IIJU
  ! compute updraft thevav and buoyancy term at KKB level             
  PTHV_UP(JI,KKB) = ZTH_UP(JI,KKB)*((1+ZRVORD*PRV_UP(JI,KKB))/(1+PRT_UP(JI,KKB))) 
  ! compute mean rsat in updraft
  PRSAT_UP(JI,KKB) = ZRSATW(JI)*(1-PFRAC_ICE_UP(JI,KKB)) + ZRSATI(JI)*PFRAC_ICE_UP(JI,KKB)
ENDDO

!Tout est commente pour tester dans un premier temps la separation en deux de la 
!  boucle verticale, une pour w et une pour PEMF                                                            

ZG_O_THVREF=XG/ZTHVM_F

! Calcul de la fermeture de Julien Pergaut comme limite max de PHY

DO JK=KKB,KKE-KKL,KKL   !  Vertical loop
  DO JI=1,IIJU
    ZZDZ(JI,JK)   = MAX(ZEPS,PZZ(JI,JK+KKL)-PZZ(JI,JK))  ! <== Delta Z between two flux level
  ENDDO
ENDDO

! compute L_up
GLMIX=.TRUE.
ZTKEM_F(:,KKB)=0.
!
IF(CTURBLEN=='RM17') THEN
  ZDUDZ = MZF_MF(GZ_M_W_MF(PUM,PDZZ, KKA, KKU, KKL), KKA, KKU, KKL)
  ZDVDZ = MZF_MF(GZ_M_W_MF(PVM,PDZZ, KKA, KKU, KKL), KKA, KKU, KKL)
  ZSHEAR = SQRT(ZDUDZ*ZDUDZ + ZDVDZ*ZDVDZ)
ELSE
  ZSHEAR = 0. !no shear in bl89 mixing length
END IF
!
CALL COMPUTE_BL89_ML(KKA,KKB,KKE,KKU,KKL,PDZZ,ZTKEM_F(:,KKB),ZG_O_THVREF(:,KKB), &
                       ZTHVM_F,KKB,GLMIX,.TRUE.,ZSHEAR,ZLUP)
ZLUP(:)=MAX(ZLUP(:),1.E-10)

DO JI=1,IIJU
  ! Compute Buoyancy flux at the ground
  ZWTHVSURF = (ZTHVM_F(JI,KKB)/ZTHM_F(JI,KKB))*PSFTH(JI)+     &
              (0.61*ZTHM_F(JI,KKB))*PSFRV(JI)

  ! Mass flux at KKB level (updraft triggered if PSFTH>0.)
  IF (ZWTHVSURF>0.010) THEN ! <==  Not 0 Important to have stratocumulus !!!!!
    PEMF(JI,KKB) = XCMF * ZRHO_F(JI,KKB) * ((ZG_O_THVREF(JI,KKB))*ZWTHVSURF*ZLUP(JI))**(1./3.)
    PFRAC_UP(JI,KKB)=MIN(PEMF(JI,KKB)/(SQRT(ZW_UP2(JI,KKB))*ZRHO_F(JI,KKB)),XFRAC_UP_MAX)
    !PEMF(JI,KKB) = ZRHO_F(JI,KKB)*PFRAC_UP(JI,KKB)*SQRT(ZW_UP2(JI,KKB))
    ZW_UP2(JI,KKB)=(PEMF(JI,KKB)/(PFRAC_UP(JI,KKB)*ZRHO_F(JI,KKB)))**2
    GTEST(JI)=.TRUE.
  ELSE
    PEMF(JI,KKB) =0.
    GTEST(JI)=.FALSE.
  ENDIF
ENDDO

  
!--------------------------------------------------------------------------

!                        3. Vertical ascending loop
!                           -----------------------
!
! If GTEST = T the updraft starts from the KKB level and stops when GTEST becomes F
!
!
GTESTLCL(:)=.FALSE.


!       Loop on vertical level to compute W

ZW_MAX(:)      = 0.
ZZTOP(:)       = 0.

DO JK=KKB,KKE-KKL,KKL

! IF the updraft top is reached for all column, stop the loop on levels

  !ITEST=COUNT(GTEST)
  !IF (ITEST==0) CYCLE

!       Computation of entrainment and detrainment with KF90
!       parameterization in clouds and LR01 in subcloud layer
  
 
! to find the LCL (check if JK is LCL or not)

  DO JI=1,IIJU
    IF ((PRC_UP(JI,JK)+PRI_UP(JI,JK)>0.).AND.(.NOT.(GTESTLCL(JI)))) THEN
      KKLCL(JI) = JK           
      GTESTLCL(JI)=.TRUE.
    ENDIF
  ENDDO


! COMPUTE PENTR and PDETR at mass level JK

    
!  Buoyancy is computed on "flux" levels where updraft variables are known

  ! Compute theta_v of updraft at flux level JK    
    
    ZRC_UP(:)   =PRC_UP(:,JK) ! guess
    ZRI_UP(:)   =PRI_UP(:,JK) ! guess 
    ZRV_UP(:)   =PRV_UP(:,JK)
    CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE_UP(:,JK),&
               PPABSM(:,JK),PTHL_UP(:,JK),PRT_UP(:,JK),&
               ZTH_UP(:,JK),ZRV_UP,ZRC_UP,ZRI_UP,ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.)            
    
  DO JI=1,IIJU
    IF (GTEST(JI)) THEN
      PTHV_UP   (JI,JK) = ZTH_UP(JI,JK)*(1.+ZRVORD*ZRV_UP(JI))/(1.+PRT_UP(JI,JK))
      ZBUO      (JI,JK) = ZG_O_THVREF(JI,JK)*(PTHV_UP(JI,JK) - ZTHVM_F(JI,JK))    
      PBUO_INTEG(JI,JK) = ZBUO(JI,JK)*(PZZ(JI,JK+KKL)-PZZ(JI,JK))
      
      ZDZ(JI)   = MAX(ZEPS,PZZ(JI,JK+KKL)-PZZ(JI,JK))
      ZTEST(JI) = XA1*ZBUO(JI,JK) -  XB*ZW_UP2(JI,JK)

      !  Ancien calcul de la vitesse
      ZCOE(JI)      = ZDZ(JI)
      IF (ZTEST(JI)>0.) THEN
        ZCOE(JI)    = ZDZ(JI)/(1.+ XBETA1)
      ENDIF

      !  Convective Vertical speed computation
      ZWCOE(JI)         = (1.-XB*ZCOE(JI))/(1.+XB*ZCOE(JI))
      ZBUCOE(JI)        =  2.*ZCOE(JI)/(1.+XB*ZCOE(JI))

      ! Second Rachel bug correction (XA1 has been forgotten)
      ZW_UP2(JI,JK+KKL) = MAX(ZEPS,ZW_UP2(JI,JK)*ZWCOE(JI) + XA1*ZBUO(JI,JK)*ZBUCOE(JI) )  
      ZW_MAX(JI) = MAX(ZW_MAX(JI), SQRT(ZW_UP2(JI,JK+KKL)))
      ZWUP_MEAN(JI)     = MAX(ZEPS,0.5*(ZW_UP2(JI,JK+KKL)+ZW_UP2(JI,JK)))
 
      !  Entrainement and detrainement

      ! First Rachel bug correction (Parenthesis around 1+beta1 ==> impact is small)   
      PENTR(JI,JK)  = MAX(0.,(XBETA1/(1.+XBETA1))*(XA1*ZBUO(JI,JK)/ZWUP_MEAN(JI)-XB))
      ZDETR_BUO(JI) = MAX(0., -(XBETA1/(1.+XBETA1))*XA1*ZBUO(JI,JK)/ZWUP_MEAN(JI))
      ZDETR_RT(JI)  = XC*SQRT(MAX(0.,(PRT_UP(JI,JK) - ZRTM_F(JI,JK))) / MAX(ZEPS,ZRTM_F(JI,JK)) / ZWUP_MEAN(JI))
      PDETR(JI,JK)  = ZDETR_RT(JI)+ZDETR_BUO(JI)
   
      ! If the updraft did not stop, compute cons updraft characteritics at jk+1
      ZZTOP(JI) = MAX(ZZTOP(JI),PZZ(JI,JK+KKL))
      ZMIX2(JI) = (PZZ(JI,JK+KKL)-PZZ(JI,JK))*PENTR(JI,JK) !&

      !ZQTM(JI) = PRTM(JI,JK)/(1.+PRTM(JI,JK))            
      !ZTHSM(JI,JK) = PTHLM(JI,JK)*(1.+XLAMBDA_MF*ZQTM(JI))
      !ZTHS_UP(JI,JK+KKL)=(ZTHS_UP(JI,JK)*(1.-0.5*ZMIX2(JI)) + ZTHSM(JI,JK)*ZMIX2(JI)) &
      !                     /(1.+0.5*ZMIX2(JI))
      PRT_UP(JI,JK+KKL) =(PRT_UP (JI,JK)*(1.-0.5*ZMIX2(JI)) + PRTM(JI,JK)*ZMIX2(JI))  &
                           /(1.+0.5*ZMIX2(JI))
      !ZQT_UP(JI) = PRT_UP(JI,JK+KKL)/(1.+PRT_UP(JI,JK+KKL))
      !PTHL_UP(JI,JK+KKL)=ZTHS_UP(JI,JK+KKL)/(1.+XLAMBDA_MF*ZQT_UP(JI))
      PTHL_UP(JI,JK+KKL)=(PTHL_UP(JI,JK)*(1.-0.5*ZMIX2(JI)) + PTHLM(JI,JK)*ZMIX2(JI)) &
                           /(1.+0.5*ZMIX2(JI))
    ENDIF  ! GTEST
  ENDDO
  

  IF(OMIXUV) THEN
    IF(JK/=KKB) THEN
      DO JI=1,IIJU
        IF(GTEST(JI)) THEN
          PU_UP(JI,JK+KKL) = (PU_UP (JI,JK)*(1-0.5*ZMIX2(JI)) + PUM(JI,JK)*ZMIX2(JI)+ &
                            0.5*XPRES_UV*(PZZ(JI,JK+KKL)-PZZ(JI,JK))*&
                            ((PUM(JI,JK+KKL)-PUM(JI,JK))/PDZZ(JI,JK+KKL)+&
                             (PUM(JI,JK)-PUM(JI,JK-KKL))/PDZZ(JI,JK))        )   &
                            /(1+0.5*ZMIX2(JI))
          PV_UP(JI,JK+KKL) = (PV_UP (JI,JK)*(1-0.5*ZMIX2(JI)) + PVM(JI,JK)*ZMIX2(JI)+ &
                            0.5*XPRES_UV*(PZZ(JI,JK+KKL)-PZZ(JI,JK))*&
                            ((PVM(JI,JK+KKL)-PVM(JI,JK))/PDZZ(JI,JK+KKL)+&
                             (PVM(JI,JK)-PVM(JI,JK-KKL))/PDZZ(JI,JK))    )   &
                            /(1+0.5*ZMIX2(JI))
        ENDIF
      ENDDO
    ELSE
      DO JI=1,IIJU
        IF(GTEST(JI)) THEN
          PU_UP(JI,JK+KKL) = (PU_UP (JI,JK)*(1-0.5*ZMIX2(JI)) + PUM(JI,JK)*ZMIX2(JI)+ &
                            0.5*XPRES_UV*(PZZ(JI,JK+KKL)-PZZ(JI,JK))*&
                            ((PUM(JI,JK+KKL)-PUM(JI,JK))/PDZZ(JI,JK+KKL))        ) &
                            /(1+0.5*ZMIX2(JI))
          PV_UP(JI,JK+KKL) = (PV_UP (JI,JK)*(1-0.5*ZMIX2(JI)) + PVM(JI,JK)*ZMIX2(JI)+ &
                            0.5*XPRES_UV*(PZZ(JI,JK+KKL)-PZZ(JI,JK))*&
                            ((PVM(JI,JK+KKL)-PVM(JI,JK))/PDZZ(JI,JK+KKL))    )   &
                            /(1+0.5*ZMIX2(JI))
        ENDIF
      ENDDO
    ENDIF
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
          ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.)

  DO JI=1,IIJU
    IF(GTEST(JI)) THEN
      !ZT_UP(JI) = ZTH_UP(JI,JK+KKL)*PEXNM(JI,JK+KKL)
      !ZCP(JI) = XCPD + XCL * ZRC_UP(JI)
      !ZLVOCPEXN(JI)=(XLVTT + (XCPV-XCL) *  (ZT_UP(JI)-XTT) ) / ZCP(JI) / PEXNM(JI,JK+KKL)
      !PRC_UP(JI,JK+KKL)=MIN(0.5E-3,ZRC_UP(JI))  ! On ne peut depasser 0.5 g/kg (autoconversion donc elimination !)
      !PTHL_UP(JI,JK+KKL) = PTHL_UP(JI,JK+KKL)+ZLVOCPEXN(JI)*(ZRC_UP(JI)-PRC_UP(JI,JK+KKL))
      PRC_UP(JI,JK+KKL)=ZRC_UP(JI)
      PRV_UP(JI,JK+KKL)=ZRV_UP(JI)
      PRI_UP(JI,JK+KKL)=ZRI_UP(JI)
      !PRT_UP(JI,JK+KKL)  = PRC_UP(JI,JK+KKL) + PRV_UP(JI,JK+KKL)
      PRSAT_UP(JI,JK+KKL) = ZRSATW(JI)*(1-PFRAC_ICE_UP(JI,JK+KKL)) + ZRSATI(JI)*PFRAC_ICE_UP(JI,JK+KKL)

      ! Compute the updraft theta_v, buoyancy and w**2 for level JK+1   
      !PTHV_UP(:,JK+KKL) = PTH_UP(:,JK+KKL)*((1+ZRVORD*PRV_UP(:,JK+KKL))/(1+PRT_UP(:,JK+KKL)))
      !PTHV_UP(JI,JK+KKL) = ZTH_UP(JI,JK+KKL)*(1.+0.608*PRV_UP(JI,JK+KKL) - PRC_UP(JI,JK+KKL))
      !! A corriger pour utiliser q et non r !!!!      
      !ZMIX1(JI)=ZZDZ(JI,JK)*(PENTR(JI,JK)-PDETR(JI,JK))
      PTHV_UP(JI,JK+KKL) = ZTH_UP(JI,JK+KKL)*((1+ZRVORD*PRV_UP(JI,JK+KKL))/(1+PRT_UP(JI,JK+KKL)))
      ZMIX1(JI)=ZZDZ(JI,JK)*(PENTR(JI,JK)-PDETR(JI,JK))
    ENDIF
  ENDDO

  DO JI=1,IIJU
    IF(GTEST(JI)) THEN
      PEMF(JI,JK+KKL)=PEMF(JI,JK)*EXP(ZMIX1(JI))
    ENDIF
  ENDDO

  DO JI=1,IIJU
    IF(GTEST(JI)) THEN
      ! Updraft fraction must be smaller than XFRAC_UP_MAX
      PFRAC_UP(JI,JK+KKL)=MIN(XFRAC_UP_MAX, &
                             &PEMF(JI,JK+KKL)/(SQRT(ZW_UP2(JI,JK+KKL))*ZRHO_F(JI,JK+KKL)))
      !PEMF(JI,JK+KKL) = ZRHO_F(JI,JK+KKL)*PFRAC_UP(JI,JK+KKL)*SQRT(ZW_UP2(JI,JK+KKL))
    ENDIF
  ENDDO

! Test if the updraft has reach the ETL
  DO JI=1,IIJU
    IF (GTEST(JI) .AND. (PBUO_INTEG(JI,JK)<=0.)) THEN
      KKETL(JI) = JK+KKL
    ENDIF
  ENDDO


! Test is we have reached the top of the updraft
  DO JI=1,IIJU
    IF (GTEST(JI) .AND. ((ZW_UP2(JI,JK+KKL)<=ZEPS).OR.(PEMF(JI,JK+KKL)<=ZEPS))) THEN
      ZW_UP2   (JI,JK+KKL)=ZEPS
      PEMF     (JI,JK+KKL)=0.
      GTEST    (JI)       =.FALSE.
      PTHL_UP  (JI,JK+KKL)=ZTHLM_F(JI,JK+KKL)
      PRT_UP   (JI,JK+KKL)=ZRTM_F(JI,JK+KKL)
      PRC_UP   (JI,JK+KKL)=0.
      PRI_UP   (JI,JK+KKL)=0.
      PRV_UP   (JI,JK+KKL)=ZRVM_F (JI,JK+KKL)
      PTHV_UP  (JI,JK+KKL)=ZTHVM_F(JI,JK+KKL)
      PFRAC_UP (JI,JK+KKL)=0.
      KKCTL    (JI)       =JK+KKL
    ENDIF
  ENDDO

ENDDO   ! Fin de la boucle verticale 

PW_UP(:,:)=SQRT(ZW_UP2(:,:))
PEMF(:,KKB) =0.

! Limits the shallow convection scheme when cloud heigth is higher than 3000m.
! To do this, mass flux is multiplied by a coefficient decreasing linearly
! from 1 (for clouds of 3000m of depth) to 0 (for clouds of 4000m of depth).
! This way, all MF fluxes are diminished by this amount.
! Diagnosed cloud fraction is also multiplied by the same coefficient.
!
DO JI=1,IIJU
   PDEPTH(JI) = MAX(0., PZZ(JI,KKCTL(JI)) -  PZZ(JI,KKLCL(JI)) )
ENDDO

GWORK1(:)= (GTESTLCL(:) .AND. (PDEPTH(:) > ZDEPTH_MAX1) )
GWORK2(:,:) = SPREAD( GWORK1(:), DIM=2, NCOPIES=IKU )
ZCOEF(:,:) = SPREAD( (1.-(PDEPTH(:)-ZDEPTH_MAX1)/(ZDEPTH_MAX2-ZDEPTH_MAX1)), DIM=2, NCOPIES=IKU)
ZCOEF(:,:)=MIN(MAX(ZCOEF(:,:),0.),1.)
DO JK=1, IKU
  DO JI=1,IIJU
    IF (GWORK2(JI,JK)) THEN
      PEMF(JI,JK)     = PEMF(JI,JK)     * ZCOEF(JI,JK) 
      PFRAC_UP(JI,JK) = PFRAC_UP(JI,JK) * ZCOEF(JI,JK) 
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAFT_RHCJ10',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_UPDRAFT_RHCJ10
END MODULE MODE_COMPUTE_UPDRAFT_RHCJ10
