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
SUBROUTINE COMPUTE_UPDRAFT_RHCJ10(D,CST,NEB,PARAMMF,TURBN,CSTURB, &
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
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST,             ONLY: CST_t
USE MODD_NEB,             ONLY: NEB_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
USE MODD_TURB_n,          ONLY: TURB_t
USE MODD_CTURB,           ONLY: CSTURB_t
!
USE MODI_SHUMAN_MF, ONLY: MZF_MF, MZM_MF, GZ_M_W_MF

USE MODE_COMPUTE_BL89_ML, ONLY: COMPUTE_BL89_ML
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK


IMPLICIT NONE

!*                    1.1  Declaration of Arguments
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
LOGICAL,                INTENT(IN)   :: OENTR_DETR! flag to recompute entrainment, detrainment and mass flux
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PZZ       !  Height at the flux point
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PDZZ      !  Metrics coefficient
 
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSFTH,PSFRV
! normal surface fluxes of theta,rv,(u,v) parallel to the orography

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PPABSM     ! Pressure at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PRHODREF   ! dry density of the
                                                  ! reference state
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PUM        ! u mean wind
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PVM        ! v mean wind
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PTKEM      ! TKE at t-dt
!
!REAL, DIMENSION(:,:),   INTENT(IN)   ::  PEXNM       ! Exner function at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PTHM           ! pot. temp. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PRVM           ! vapor mixing ratio at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PTHLM,PRTM     ! cons. var. at t-dt

REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN)   ::  PSVM           ! scalar var. at t-dt

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  ::  PTHL_UP,PRT_UP   ! updraft properties
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  ::  PU_UP, PV_UP     ! updraft wind components
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT)::  PRV_UP,PRC_UP    ! updraft rv, rc
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT)::  PRI_UP           ! updraft ri
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT)::  PTHV_UP          ! updraft THv
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT)::  PW_UP,PFRAC_UP   ! updraft w, fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT)::  PFRAC_ICE_UP     ! liquid/solid fraction in updraft
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT)::  PRSAT_UP         ! Rsat

REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT)  ::  PSV_UP           ! updraft scalar var. 
                                         
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT)::  PEMF,PDETR,PENTR ! Mass_flux,
                                                          ! detrainment,entrainment
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PBUO_INTEG       ! Integrated Buoyancy 
INTEGER, DIMENSION(D%NIJT),  INTENT(INOUT)::  KKLCL,KKETL,KKCTL! LCL, ETL, CTL                                     
REAL, DIMENSION(D%NIJT),     INTENT(OUT)   :: PDEPTH           ! Deepness of cloud
!                       1.2  Declaration of local variables
!
! Mean environment variables at t-dt at flux point
REAL, DIMENSION(D%NIJT,D%NKT) ::    ZTHM_F,ZRVM_F    ! Theta,rv of
                                                                  ! updraft environnement
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRTM_F, ZTHLM_F, ZTKEM_F      ! rt, thetal,TKE,pressure,
REAL, DIMENSION(D%NIJT,D%NKT) :: ZUM_F,ZVM_F,ZRHO_F            ! density,momentum
REAL, DIMENSION(D%NIJT,D%NKT) :: ZPRES_F,ZTHVM_F               ! interpolated at the flux point
REAL, DIMENSION(D%NIJT,D%NKT) :: ZG_O_THVREF                   ! g*ThetaV ref
REAL, DIMENSION(D%NIJT,D%NKT) :: ZW_UP2                        ! w**2  of the updraft

REAL, DIMENSION(D%NIJT,D%NKT,KSV) :: ZSVM_F ! scalar variables 
                        

                        
REAL, DIMENSION(D%NIJT,D%NKT) :: ZTH_UP                        ! updraft THETA 
!REAL, DIMENSION(SIZE(PTHM,1))              :: ZT_UP                         ! updraft T
!REAL, DIMENSION(SIZE(PTHM,1))              :: ZLVOCPEXN                     ! updraft L
!REAL, DIMENSION(SIZE(PTHM,1))              :: ZCP                           ! updraft cp
REAL, DIMENSION(D%NIJT,D%NKT) :: ZBUO                          ! Buoyancy 
!REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZTHS_UP,ZTHSM

REAL, DIMENSION(D%NIJT,D%NKT) ::  ZCOEF  ! diminution coefficient for too high clouds 
                        
REAL                                       ::  ZWTHVSURF  ! Surface w'thetav'

REAL  :: ZRVORD       ! RV/RD


REAL, DIMENSION(D%NIJT) :: ZMIX1,ZMIX2

REAL, DIMENSION(D%NIJT) :: ZLUP         ! Upward Mixing length from the ground

INTEGER  :: JK,JIJ,JSV          ! loop counters
INTEGER :: IIJB,IIJE ! physical horizontal domain indices
INTEGER :: IKT,IKB,IKE,IKL
LOGICAL, DIMENSION(D%NIJT) ::  GTEST,GTESTLCL
                               ! Test if the ascent continue, if LCL or ETL is reached
LOGICAL                          ::  GLMIX 
                               ! To choose upward or downward mixing length
LOGICAL, DIMENSION(D%NIJT)              :: GWORK1
LOGICAL, DIMENSION(D%NIJT,D%NKT) :: GWORK2

INTEGER  :: ITEST

REAL, DIMENSION(D%NIJT) :: ZRC_UP, ZRI_UP, ZRV_UP, ZRSATW, ZRSATI

REAL,    DIMENSION(D%NIJT,D%NKT) :: ZZDZ

REAL, DIMENSION(D%NIJT)              ::  ZTEST,ZDZ,ZWUP_MEAN    ! 
REAL, DIMENSION(D%NIJT)              ::  ZCOE,ZWCOE,ZBUCOE
REAL, DIMENSION(D%NIJT)              ::  ZDETR_BUO, ZDETR_RT
REAL, DIMENSION(D%NIJT)              ::  ZW_MAX               ! w**2  max of the updraft
REAL, DIMENSION(D%NIJT)              ::  ZZTOP                ! Top of the updraft
!REAL, DIMENSION(SIZE(PTHM,1))              ::  ZQTM,ZQT_UP

REAL  :: ZDEPTH_MAX1, ZDEPTH_MAX2 ! control auto-extinction process

REAL  :: ZTMAX,ZRMAX, ZEPS  ! control value

REAL, DIMENSION(D%NIJT,D%NKT) :: ZSHEAR,ZDUDZ,ZDVDZ ! vertical wind shear
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWK
REAL, DIMENSION(D%NIJT,16)    :: ZBUF
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAFT_RHCJ10',0,ZHOOK_HANDLE)
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
ZEPS=1.E-15
!------------------------------------------------------------------------
!                     INITIALISATION

! Initialisation of the constants   
ZRVORD   = (CST%XRV / CST%XRD) 

! depth are different in compute_updraft (3000. and 4000.) ==> impact is small
ZDEPTH_MAX1=4500. ! clouds with depth infeRIOr to this value are keeped untouched
ZDEPTH_MAX2=5000. ! clouds with depth superior to this value are suppressed


!                 Local variables, internal domain

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
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PRSAT_UP(:,:)=PRVM(:,:) ! should be initialised correctly but is (normaly) not used
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

! Initialisation of environment variables at t-dt

! variables at flux level
CALL MZM_MF(D, PTHLM(:,:), ZTHLM_F(:,:))
CALL MZM_MF(D, PRTM(:,:), ZRTM_F(:,:))
CALL MZM_MF(D, PUM(:,:), ZUM_F(:,:))
CALL MZM_MF(D, PVM(:,:), ZVM_F(:,:))
CALL MZM_MF(D, PTKEM(:,:), ZTKEM_F(:,:))

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
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PTHL_UP(:,:)=ZTHLM_F(:,:)
PRT_UP(:,:)=ZRTM_F(:,:)
PU_UP(:,:)=ZUM_F(:,:)
PV_UP(:,:)=ZVM_F(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PSV_UP(:,:,:)=0.
! This updraft is not yet ready to use scalar variables
!IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) then
!    PSV_UP(:,:,:)=ZSVM_F(:,:,:)
!ENDIF

! Computation or initialisation of updraft characteristics at the KKB level
! thetal_up,rt_up,thetaV_up, w,Buoyancy term and mass flux (PEMF)

DO JIJ=IIJB,IIJE
  !PTHL_UP(JIJ,KKB)= ZTHLM_F(JIJ,KKB)+MAX(0.,MIN(ZTMAX,(PSFTH(JIJ)/SQRT(ZTKEM_F(JIJ,KKB)))*XALP_PERT))
  !PRT_UP(JIJ,KKB) = ZRTM_F(JIJ,KKB)+MAX(0.,MIN(ZRMAX,(PSFRV(JIJ)/SQRT(ZTKEM_F(JIJ,KKB)))*XALP_PERT)) 
  PTHL_UP(JIJ,IKB)= ZTHLM_F(JIJ,IKB)
  PRT_UP(JIJ,IKB) = ZRTM_F(JIJ,IKB)
  !ZQT_UP(JIJ) = PRT_UP(JIJ,KKB)/(1.+PRT_UP(JIJ,KKB))
  !ZTHS_UP(JIJ,KKB)=PTHL_UP(JIJ,KKB)*(1.+XLAMBDA_MF*ZQT_UP(JIJ))
ENDDO

CALL MZM_MF(D, PTHM (:,:), ZTHM_F(:,:))
CALL MZM_MF(D, PPABSM(:,:), ZPRES_F(:,:))
CALL MZM_MF(D, PRHODREF(:,:), ZRHO_F(:,:))
CALL MZM_MF(D, PRVM(:,:), ZRVM_F(:,:))

! thetav at mass and flux levels 
DO JK=1,IKT
  DO JIJ=D%NIB,D%NIJE
    ZTHVM_F(JIJ,JK)=ZTHM_F(JIJ,JK)*((1.+ZRVORD*ZRVM_F(JIJ,JK))/(1.+ZRTM_F(JIJ,JK)))
  ENDDO
ENDDO

!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PTHV_UP(:,:)= ZTHVM_F(:,:)
PRV_UP(:,:)= ZRVM_F(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

ZW_UP2(:,:)=ZEPS
!$mnh_expand_array(JIJ=IIJB:IIJE)
!ZW_UP2(:,KKB) = MAX(0.0001,(3./6.)*ZTKEM_F(:,KKB))
ZW_UP2(:,IKB) = MAX(0.0001,(2./3.)*ZTKEM_F(:,IKB))
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

! Computation of non conservative variable for the KKB level of the updraft
! (all or nothing ajustement)

!$mnh_expand_array(JIJ=IIJB:IIJE)
PRC_UP(:,IKB)=0.
PRI_UP(:,IKB)=0.
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
CALL TH_R_FROM_THL_RT(CST,NEB,D%NIJT,HFRAC_ICE,PFRAC_ICE_UP(:,IKB),ZPRES_F(:,IKB), &
             PTHL_UP(:,IKB),PRT_UP(:,IKB),ZTH_UP(:,IKB), &
             PRV_UP(:,IKB),PRC_UP(:,IKB),PRI_UP(:,IKB),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.,&
             PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)

DO JIJ=IIJB,IIJE
  ! compute updraft thevav and buoyancy term at KKB level             
  PTHV_UP(JIJ,IKB) = ZTH_UP(JIJ,IKB)*((1+ZRVORD*PRV_UP(JIJ,IKB))/(1+PRT_UP(JIJ,IKB))) 
  ! compute mean rsat in updraft
  PRSAT_UP(JIJ,IKB) = ZRSATW(JIJ)*(1-PFRAC_ICE_UP(JIJ,IKB)) + ZRSATI(JIJ)*PFRAC_ICE_UP(JIJ,IKB)
ENDDO

!Tout est commente pour tester dans un premier temps la separation en deux de la 
!  boucle verticale, une pour w et une pour PEMF                                                            

!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZG_O_THVREF(:,:)=CST%XG/ZTHVM_F(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

! Calcul de la fermeture de Julien Pergaut comme limite max de PHY

DO JK=IKB,IKE-IKL,IKL   !  Vertical loop
  DO JIJ=IIJB,IIJE
    ZZDZ(JIJ,JK)   = MAX(ZEPS,PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))  ! <== Delta Z between two flux level
  ENDDO
ENDDO

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
  ZSHEAR(:,:) = 0. !no shear in bl89 mixing length
END IF
!
CALL COMPUTE_BL89_ML(D, CST, CSTURB, PDZZ,ZTKEM_F(:,IKB),ZG_O_THVREF(:,IKB), &
                       ZTHVM_F,IKB,GLMIX,.TRUE.,ZSHEAR,ZLUP)
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZLUP(:)=MAX(ZLUP(:),1.E-10)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

DO JIJ=IIJB,IIJE
  ! Compute Buoyancy flux at the ground
  ZWTHVSURF = (ZTHVM_F(JIJ,IKB)/ZTHM_F(JIJ,IKB))*PSFTH(JIJ)+     &
              (0.61*ZTHM_F(JIJ,IKB))*PSFRV(JIJ)

  ! Mass flux at KKB level (updraft triggered if PSFTH>0.)
  IF (ZWTHVSURF>0.010) THEN ! <==  Not 0 Important to have stratocumulus !!!!!
    PEMF(JIJ,IKB) = PARAMMF%XCMF * ZRHO_F(JIJ,IKB) * ((ZG_O_THVREF(JIJ,IKB))*ZWTHVSURF*ZLUP(JIJ))**(1./3.)
    PFRAC_UP(JIJ,IKB)=MIN(PEMF(JIJ,IKB)/(SQRT(ZW_UP2(JIJ,IKB))*ZRHO_F(JIJ,IKB)),PARAMMF%XFRAC_UP_MAX)
    !PEMF(JIJ,KKB) = ZRHO_F(JIJ,KKB)*PFRAC_UP(JIJ,KKB)*SQRT(ZW_UP2(JIJ,KKB))
    ZW_UP2(JIJ,IKB)=(PEMF(JIJ,IKB)/(PFRAC_UP(JIJ,IKB)*ZRHO_F(JIJ,IKB)))**2
    GTEST(JIJ)=.TRUE.
  ELSE
    PEMF(JIJ,IKB) =0.
    GTEST(JIJ)=.FALSE.
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

DO JK=IKB,IKE-IKL,IKL

! IF the updraft top is reached for all column, stop the loop on levels

  !ITEST=COUNT(GTEST)
  !IF (ITEST==0) CYCLE

!       Computation of entrainment and detrainment with KF90
!       parameterization in clouds and LR01 in subcloud layer
  
 
! to find the LCL (check if JK is LCL or not)

  DO JIJ=IIJB,IIJE
    IF ((PRC_UP(JIJ,JK)+PRI_UP(JIJ,JK)>0.).AND.(.NOT.(GTESTLCL(JIJ)))) THEN
      KKLCL(JIJ) = JK           
      GTESTLCL(JIJ)=.TRUE.
    ENDIF
  ENDDO


! COMPUTE PENTR and PDETR at mass level JK

    
!  Buoyancy is computed on "flux" levels where updraft variables are known

  ! Compute theta_v of updraft at flux level JK    
    
    !$mnh_expand_array(JIJ=IIJB:IIJE)
    ZRC_UP(:)   =PRC_UP(:,JK) ! guess
    ZRI_UP(:)   =PRI_UP(:,JK) ! guess 
    ZRV_UP(:)   =PRV_UP(:,JK)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE)
    CALL TH_R_FROM_THL_RT(CST,NEB, D%NIJT, HFRAC_ICE,PFRAC_ICE_UP(:,JK),&
               PPABSM(:,JK),PTHL_UP(:,JK),PRT_UP(:,JK),&
               ZTH_UP(:,JK),ZRV_UP,ZRC_UP,ZRI_UP,ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.,&
               PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)
    
  DO JIJ=IIJB,IIJE
    IF (GTEST(JIJ)) THEN
      PTHV_UP(JIJ,JK)    = ZTH_UP(JIJ,JK)*(1.+ZRVORD*ZRV_UP(JIJ))/(1.+PRT_UP(JIJ,JK))
      ZBUO(JIJ,JK)       = ZG_O_THVREF(JIJ,JK)*(PTHV_UP(JIJ,JK) - ZTHVM_F(JIJ,JK))    
      PBUO_INTEG(JIJ,JK) = ZBUO(JIJ,JK)*(PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))
      
      ZDZ(JIJ)   = MAX(ZEPS,PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))
      ZTEST(JIJ) = PARAMMF%XA1*ZBUO(JIJ,JK) -  PARAMMF%XB*ZW_UP2(JIJ,JK)

      !  Ancien calcul de la vitesse
      ZCOE(JIJ)      = ZDZ(JIJ)
      IF (ZTEST(JIJ)>0.) THEN
        ZCOE(JIJ)    = ZDZ(JIJ)/(1.+ PARAMMF%XBETA1)
      ENDIF

      !  Convective Vertical speed computation
      ZWCOE(JIJ)         = (1.-PARAMMF%XB*ZCOE(JIJ))/(1.+PARAMMF%XB*ZCOE(JIJ))
      ZBUCOE(JIJ)        =  2.*ZCOE(JIJ)/(1.+PARAMMF%XB*ZCOE(JIJ))

      ! Second Rachel bug correction (XA1 has been forgotten)
      ZW_UP2(JIJ,JK+IKL) = MAX(ZEPS,ZW_UP2(JIJ,JK)*ZWCOE(JIJ) + PARAMMF%XA1*ZBUO(JIJ,JK)*ZBUCOE(JIJ) )  
      ZW_MAX(JIJ) = MAX(ZW_MAX(JIJ), SQRT(ZW_UP2(JIJ,JK+IKL)))
      ZWUP_MEAN(JIJ)     = MAX(ZEPS,0.5*(ZW_UP2(JIJ,JK+IKL)+ZW_UP2(JIJ,JK)))
 
      !  Entrainement and detrainement

      ! First Rachel bug correction (Parenthesis around 1+beta1 ==> impact is small)   
      PENTR(JIJ,JK)  = MAX(0.,(PARAMMF%XBETA1/(1.+PARAMMF%XBETA1))*(PARAMMF%XA1*ZBUO(JIJ,JK)/ZWUP_MEAN(JIJ)-PARAMMF%XB))
      ZDETR_BUO(JIJ) = MAX(0., -(PARAMMF%XBETA1/(1.+PARAMMF%XBETA1))*PARAMMF%XA1*ZBUO(JIJ,JK)/ZWUP_MEAN(JIJ))
      ZDETR_RT(JIJ)  = PARAMMF%XC*SQRT(MAX(0.,(PRT_UP(JIJ,JK) - ZRTM_F(JIJ,JK))) / MAX(ZEPS,ZRTM_F(JIJ,JK)) / ZWUP_MEAN(JIJ))
      PDETR(JIJ,JK)  = ZDETR_RT(JIJ)+ZDETR_BUO(JIJ)
   
      ! If the updraft did not stop, compute cons updraft characteritics at jk+1
      ZZTOP(JIJ) = MAX(ZZTOP(JIJ),PZZ(JIJ,JK+IKL))
      ZMIX2(JIJ) = (PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*PENTR(JIJ,JK) !&

      !ZQTM(JIJ) = PRTM(JIJ,JK)/(1.+PRTM(JIJ,JK))            
      !ZTHSM(JIJ,JK) = PTHLM(JIJ,JK)*(1.+XLAMBDA_MF*ZQTM(JIJ))
      !ZTHS_UP(JIJ,JK+KKL)=(ZTHS_UP(JIJ,JK)*(1.-0.5*ZMIX2(JIJ)) + ZTHSM(JIJ,JK)*ZMIX2(JIJ)) &
      !                     /(1.+0.5*ZMIX2(JIJ))
      PRT_UP(JIJ,JK+IKL) =(PRT_UP (JIJ,JK)*(1.-0.5*ZMIX2(JIJ)) + PRTM(JIJ,JK)*ZMIX2(JIJ))  &
                           /(1.+0.5*ZMIX2(JIJ))
      !ZQT_UP(JIJ) = PRT_UP(JIJ,JK+KKL)/(1.+PRT_UP(JIJ,JK+KKL))
      !PTHL_UP(JIJ,JK+KKL)=ZTHS_UP(JIJ,JK+KKL)/(1.+XLAMBDA_MF*ZQT_UP(JIJ))
      PTHL_UP(JIJ,JK+IKL)=(PTHL_UP(JIJ,JK)*(1.-0.5*ZMIX2(JIJ)) + PTHLM(JIJ,JK)*ZMIX2(JIJ)) &
                           /(1.+0.5*ZMIX2(JIJ))
    ENDIF  ! GTEST
  ENDDO
  

  IF(PARAMMF%LMIXUV) THEN
    IF(JK/=IKB) THEN
      DO JIJ=IIJB,IIJE
        IF(GTEST(JIJ)) THEN
          PU_UP(JIJ,JK+IKL) = (PU_UP (JIJ,JK)*(1-0.5*ZMIX2(JIJ)) + PUM(JIJ,JK)*ZMIX2(JIJ)+ &
                            0.5*PARAMMF%XPRES_UV*(PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*&
                            ((PUM(JIJ,JK+IKL)-PUM(JIJ,JK))/PDZZ(JIJ,JK+IKL)+&
                             (PUM(JIJ,JK)-PUM(JIJ,JK-IKL))/PDZZ(JIJ,JK))        )   &
                            /(1+0.5*ZMIX2(JIJ))
          PV_UP(JIJ,JK+IKL) = (PV_UP (JIJ,JK)*(1-0.5*ZMIX2(JIJ)) + PVM(JIJ,JK)*ZMIX2(JIJ)+ &
                            0.5*PARAMMF%XPRES_UV*(PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*&
                            ((PVM(JIJ,JK+IKL)-PVM(JIJ,JK))/PDZZ(JIJ,JK+IKL)+&
                             (PVM(JIJ,JK)-PVM(JIJ,JK-IKL))/PDZZ(JIJ,JK))    )   &
                            /(1+0.5*ZMIX2(JIJ))
        ENDIF
      ENDDO
    ELSE
      DO JIJ=IIJB,IIJE
        IF(GTEST(JIJ)) THEN
          PU_UP(JIJ,JK+IKL) = (PU_UP (JIJ,JK)*(1-0.5*ZMIX2(JIJ)) + PUM(JIJ,JK)*ZMIX2(JIJ)+ &
                            0.5*PARAMMF%XPRES_UV*(PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*&
                            ((PUM(JIJ,JK+IKL)-PUM(JIJ,JK))/PDZZ(JIJ,JK+IKL))        ) &
                            /(1+0.5*ZMIX2(JIJ))
          PV_UP(JIJ,JK+IKL) = (PV_UP (JIJ,JK)*(1-0.5*ZMIX2(JIJ)) + PVM(JIJ,JK)*ZMIX2(JIJ)+ &
                            0.5*PARAMMF%XPRES_UV*(PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*&
                            ((PVM(JIJ,JK+IKL)-PVM(JIJ,JK))/PDZZ(JIJ,JK+IKL))    )   &
                            /(1+0.5*ZMIX2(JIJ))
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
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZRC_UP(:)=PRC_UP(:,JK) ! guess = level just below
  ZRI_UP(:)=PRI_UP(:,JK) ! guess = level just below
  ZRV_UP(:)=PRV_UP(:,JK)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  CALL TH_R_FROM_THL_RT(CST,NEB, D%NIJT, HFRAC_ICE,PFRAC_ICE_UP(:,JK+IKL),ZPRES_F(:,JK+IKL), &
          PTHL_UP(:,JK+IKL),PRT_UP(:,JK+IKL),ZTH_UP(:,JK+IKL),              &
          ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.,&
          PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)

  DO JIJ=IIJB,IIJE
    IF(GTEST(JIJ)) THEN
      !ZT_UP(JIJ) = ZTH_UP(JIJ,JK+KKL)*PEXNM(JIJ,JK+KKL)
      !ZCP(JIJ) = XCPD + XCL * ZRC_UP(JIJ)
      !ZLVOCPEXN(JIJ)=(XLVTT + (XCPV-XCL) *  (ZT_UP(JIJ)-XTT) ) / ZCP(JIJ) / PEXNM(JIJ,JK+KKL)
      !PRC_UP(JIJ,JK+KKL)=MIN(0.5E-3,ZRC_UP(JIJ))  ! On ne peut depasser 0.5 g/kg (autoconversion donc elimination !)
      !PTHL_UP(JIJ,JK+KKL) = PTHL_UP(JIJ,JK+KKL)+ZLVOCPEXN(JIJ)*(ZRC_UP(JIJ)-PRC_UP(JIJ,JK+KKL))
      PRC_UP(JIJ,JK+IKL)=ZRC_UP(JIJ)
      PRV_UP(JIJ,JK+IKL)=ZRV_UP(JIJ)
      PRI_UP(JIJ,JK+IKL)=ZRI_UP(JIJ)
      !PRT_UP(JIJ,JK+KKL)  = PRC_UP(JIJ,JK+KKL) + PRV_UP(JIJ,JK+KKL)
      PRSAT_UP(JIJ,JK+IKL) = ZRSATW(JIJ)*(1-PFRAC_ICE_UP(JIJ,JK+IKL)) + ZRSATI(JIJ)*PFRAC_ICE_UP(JIJ,JK+IKL)

      ! Compute the updraft theta_v, buoyancy and w**2 for level JK+1   
      !PTHV_UP(:,JK+KKL) = PTH_UP(:,JK+KKL)*((1+ZRVORD*PRV_UP(:,JK+KKL))/(1+PRT_UP(:,JK+KKL)))
      !PTHV_UP(JIJ,JK+KKL) = ZTH_UP(JIJ,JK+KKL)*(1.+0.608*PRV_UP(JIJ,JK+KKL) - PRC_UP(JIJ,JK+KKL))
      !! A corriger pour utiliser q et non r !!!!      
      !ZMIX1(JIJ)=ZZDZ(JIJ,JK)*(PENTR(JIJ,JK)-PDETR(JIJ,JK))
      PTHV_UP(JIJ,JK+IKL) = ZTH_UP(JIJ,JK+IKL)*((1+ZRVORD*PRV_UP(JIJ,JK+IKL))/(1+PRT_UP(JIJ,JK+IKL)))
      ZMIX1(JIJ)=ZZDZ(JIJ,JK)*(PENTR(JIJ,JK)-PDETR(JIJ,JK))
    ENDIF
  ENDDO

  DO JIJ=IIJB,IIJE
    IF(GTEST(JIJ)) THEN
      PEMF(JIJ,JK+IKL)=PEMF(JIJ,JK)*EXP(ZMIX1(JIJ))
    ENDIF
  ENDDO

  DO JIJ=IIJB,IIJE
    IF(GTEST(JIJ)) THEN
      ! Updraft fraction must be smaller than XFRAC_UP_MAX
      PFRAC_UP(JIJ,JK+IKL)=MIN(PARAMMF%XFRAC_UP_MAX, &
                             &PEMF(JIJ,JK+IKL)/(SQRT(ZW_UP2(JIJ,JK+IKL))*ZRHO_F(JIJ,JK+IKL)))
      !PEMF(JIJ,JK+KKL) = ZRHO_F(JIJ,JK+KKL)*PFRAC_UP(JIJ,JK+KKL)*SQRT(ZW_UP2(JIJ,JK+KKL))
    ENDIF
  ENDDO

! Test if the updraft has reach the ETL
  DO JIJ=IIJB,IIJE
    IF (GTEST(JIJ) .AND. (PBUO_INTEG(JIJ,JK)<=0.)) THEN
      KKETL(JIJ) = JK+IKL
    ENDIF
  ENDDO


! Test is we have reached the top of the updraft
  DO JIJ=IIJB,IIJE
    IF (GTEST(JIJ) .AND. ((ZW_UP2(JIJ,JK+IKL)<=ZEPS).OR.(PEMF(JIJ,JK+IKL)<=ZEPS))) THEN
      ZW_UP2   (JIJ,JK+IKL)=ZEPS
      PEMF     (JIJ,JK+IKL)=0.
      GTEST    (JIJ)       =.FALSE.
      PTHL_UP  (JIJ,JK+IKL)=ZTHLM_F(JIJ,JK+IKL)
      PRT_UP   (JIJ,JK+IKL)=ZRTM_F(JIJ,JK+IKL)
      PRC_UP   (JIJ,JK+IKL)=0.
      PRI_UP   (JIJ,JK+IKL)=0.
      PRV_UP   (JIJ,JK+IKL)=ZRVM_F (JIJ,JK+IKL)
      PTHV_UP  (JIJ,JK+IKL)=ZTHVM_F(JIJ,JK+IKL)
      PFRAC_UP (JIJ,JK+IKL)=0.
      KKCTL    (JIJ)       =JK+IKL
    ENDIF
  ENDDO

ENDDO   ! Fin de la boucle verticale 

!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PW_UP(:,:)=SQRT(ZW_UP2(:,:))
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!$mnh_expand_array(JIJ=IIJB:IIJE)
PEMF(:,IKB) =0.
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

! Limits the shallow convection scheme when cloud heigth is higher than 3000m.
! To do this, mass flux is multiplied by a coefficient decreasing linearly
! from 1 (for clouds of 3000m of depth) to 0 (for clouds of 4000m of depth).
! This way, all MF fluxes are diminished by this amount.
! Diagnosed cloud fraction is also multiplied by the same coefficient.
!
DO JIJ=IIJB,IIJE
   PDEPTH(JIJ) = MAX(0., PZZ(JIJ,KKCTL(JIJ)) -  PZZ(JIJ,KKLCL(JIJ)) )
ENDDO

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
DO JK=1,IKT
  DO JIJ=IIJB,IIJE
    IF (GWORK2(JIJ,JK)) THEN
      PEMF(JIJ,JK)     = PEMF(JIJ,JK)     * ZCOEF(JIJ,JK) 
      PFRAC_UP(JIJ,JK) = PFRAC_UP(JIJ,JK) * ZCOEF(JIJ,JK) 
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAFT_RHCJ10',1,ZHOOK_HANDLE)
!
CONTAINS
INCLUDE "th_r_from_thl_rt.func.h"
INCLUDE "compute_frac_ice.func.h"
!
END SUBROUTINE COMPUTE_UPDRAFT_RHCJ10
END MODULE MODE_COMPUTE_UPDRAFT_RHCJ10
