!MNH_LIC Copyright 2012-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
     MODULE MODE_COMPUTE_UPDRAFT_RAHA
!    ###########################
!
IMPLICIT NONE
CONTAINS
      SUBROUTINE COMPUTE_UPDRAFT_RAHA(D, CST, NEB, PARAMMF,          &
                                 KSV, HFRAC_ICE, OENTR_DETR,         &
                                 ONOMIXLG,KSV_LGBEG,KSV_LGEND,       &
                                 PZZ,PDZZ,                           &
                                 PSFTH,PSFRV,                        &
                                 PPABSM,PRHODREF,PUM,PVM, PTKEM,     &
                                 PEXNM,PTHM,PRVM,PTHLM,PRTM,         &
                                 PSVM,PTHL_UP,PRT_UP,                &
                                 PRV_UP,PRC_UP,PRI_UP,PTHV_UP,       &
                                 PW_UP,PU_UP, PV_UP, PSV_UP,         &
                                 PFRAC_UP,PFRAC_ICE_UP,PRSAT_UP,     &
                                 PEMF,PDETR,PENTR,                   &
                                 PBUO_INTEG,KKLCL,KKETL,KKCTL,       &
                                 PDEPTH     )

!     #################################################################
!!
!!****  *COMPUTE_UPDRAF_RAHA* - calculates caracteristics of the updraft 
!!                         
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to build the updraft following Rio et al (2010)
!!      Same as compute_updraft_rhcj10 exept the use of Hourdin et al closure
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
!!       Hourdin et al (xxxx)
!!
!!     AUTHOR
!!     ------
!!     Y. Bouteloup (2012)
!!     R. Honnert Janv 2013 ==> corection of some coding bugs
!!     Y. Bouteloup Janv 2014 ==> Allow the use of loops in the both direction
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST,             ONLY: CST_t
USE MODD_NEB,             ONLY: NEB_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
!
USE MODI_SHUMAN_MF, ONLY: MZM_MF
!
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

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PEXNM       ! Exner function at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PTHM           ! liquid pot. temp. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PRVM           ! vapor mixing ratio at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PTHLM,PRTM     ! cons. var. at t-dt

REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN)   ::  PSVM           ! scalar var. at t-dt

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  ::  PTHL_UP,PRT_UP   ! updraft properties
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  ::  PU_UP, PV_UP     ! updraft wind components
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT)::  PRV_UP,PRC_UP    ! updraft rv, rc
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT)::  PRI_UP,PTHV_UP   ! updraft ri, THv
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
!
! Mean environment variables at t-dt at flux point
REAL, DIMENSION(D%NIJT,D%NKT) ::    ZTHM_F,ZRVM_F,ZRCM_F    ! Theta,rv of
                                                                  ! updraft environnement
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRTM_F, ZTHLM_F, ZTKEM_F      ! rt, thetal,TKE,pressure,
REAL, DIMENSION(D%NIJT,D%NKT) :: ZUM_F,ZVM_F,ZRHO_F            ! density,momentum
REAL, DIMENSION(D%NIJT,D%NKT) :: ZPRES_F,ZTHVM_F,ZTHVM         ! interpolated at the flux point
REAL, DIMENSION(D%NIJT,D%NKT) :: ZG_O_THVREF                   ! g*ThetaV ref
REAL, DIMENSION(D%NIJT,D%NKT) :: ZW_UP2                        ! w**2  of the updraft

REAL, DIMENSION(D%NIJT,D%NKT,KSV) :: ZSVM_F ! scalar variables 
                        

                        
REAL, DIMENSION(D%NIJT,D%NKT) :: ZTH_UP                        ! updraft THETA 
REAL, DIMENSION(D%NIJT)              :: ZT_UP                         ! updraft T
REAL, DIMENSION(D%NIJT)              :: ZLVOCPEXN                     ! updraft L
REAL, DIMENSION(D%NIJT)              :: ZCP                           ! updraft cp
REAL, DIMENSION(D%NIJT,D%NKT) :: ZBUO                          ! Buoyancy 
REAL, DIMENSION(D%NIJT,D%NKT) :: ZTHS_UP,ZTHSM

REAL, DIMENSION(D%NIJT,D%NKT) ::  ZCOEF  ! diminution coefficient for too high clouds 
                        
REAL, DIMENSION(D%NIJT)            ::  ZWTHVSURF  ! Surface w'thetav'

REAL  :: ZRDORV       ! RD/RV
REAL  :: ZRVORD       ! RV/RD


REAL, DIMENSION(D%NIJT) :: ZMIX1,ZMIX2,ZMIX3

REAL, DIMENSION(D%NIJT) :: ZLUP         ! Upward Mixing length from the ground

REAL, DIMENSION(D%NIJT) :: ZDEPTH       ! Deepness limit for cloud

INTEGER  :: JK,JIJ,JSV          ! loop counters
INTEGER :: IIJB,IIJE ! physical horizontal domain indices
INTEGEr :: IKT,IKB,IKE,IKL

LOGICAL, DIMENSION(D%NIJT) ::  GTEST,GTESTLCL,GTESTETL
                               ! Test if the ascent continue, if LCL or ETL is reached
LOGICAL                          ::  GLMIX 
                               ! To choose upward or downward mixing length
LOGICAL, DIMENSION(D%NIJT)              :: GWORK1
LOGICAL, DIMENSION(D%NIJT,D%NKT) :: GWORK2


INTEGER  :: ITEST

REAL, DIMENSION(D%NIJT) :: ZRC_UP, ZRI_UP, ZRV_UP, ZWP2, ZRSATW, ZRSATI

LOGICAL, DIMENSION(D%NIJT) :: GTEST_FER
REAL,    DIMENSION(D%NIJT) :: ZPHI,ZALIM_STAR_TOT
REAL,    DIMENSION(D%NIJT,D%NKT) :: ZDTHETASDZ,ZALIM_STAR,ZZDZ,ZZZ
INTEGER, DIMENSION(D%NIJT) :: IALIM

REAL, DIMENSION(D%NIJT)              ::  ZTEST,ZDZ,ZWUP_MEAN    ! 
REAL, DIMENSION(D%NIJT)              ::  ZCOE,ZWCOE,ZBUCOE
REAL, DIMENSION(D%NIJT)              ::  ZDETR_BUO, ZDETR_RT
REAL, DIMENSION(D%NIJT)              ::  ZW_MAX               ! w**2  max of the updraft
REAL, DIMENSION(D%NIJT)              ::  ZZTOP                ! Top of the updraft
REAL, DIMENSION(D%NIJT)              ::  ZA,ZB,ZQTM,ZQT_UP

REAL  :: ZDEPTH_MAX1, ZDEPTH_MAX2 ! control auto-extinction process

REAL  :: ZTMAX,ZRMAX, ZEPS  ! control value
REAL, DIMENSION(D%NIJT,16)           ::  ZBUF

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAF_RAHA',0,ZHOOK_HANDLE)
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
ZRDORV   = CST%XRD / CST%XRV   !=0.622
ZRVORD   = (CST%XRV / CST%XRD) 

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
PRV_UP(:,:)=0.
PRC_UP(:,:)=0.

PW_UP(:,:)=0.
ZTH_UP(:,:)=0.
PFRAC_UP(:,:)=0.
PTHV_UP(:,:)=0.

PBUO_INTEG(:,:)=0.
ZBUO(:,:)      =0.

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

!DO JSV=1,ISV 
! IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
!   ZSVM_F(:,KKB:IKU,JSV) = 0.5*(PSVM(:,KKB:IKU,JSV)+PSVM(:,1:IKU-1,JSV))
!   ZSVM_F(:,1,JSV)       = ZSVM_F(:,KKB,JSV) 
!END DO  

!          Initialisation of updraft characteristics 
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PTHL_UP(:,:)=ZTHLM_F(:,:)
PRT_UP(:,:)=ZRTM_F(:,:)
PU_UP(:,:)=ZUM_F(:,:)
PV_UP(:,:)=ZVM_F(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PSV_UP(:,:,:)=0.
!IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) then
!    PSV_UP(:,:,:)=ZSVM_F(:,:,:)
!ENDIF

! Computation or initialisation of updraft characteristics at the KKB level
! thetal_up,rt_up,thetaV_up, w�,Buoyancy term and mass flux (PEMF)

!$mnh_expand_array(JIJ=IIJB:IIJE)
PTHL_UP(:,IKB)= ZTHLM_F(:,IKB)+ &
                          & MAX(0.,MIN(ZTMAX,(PSFTH(:)/SQRT(ZTKEM_F(:,IKB)))*PARAMMF%XALP_PERT))
PRT_UP(:,IKB) = ZRTM_F(:,IKB)+ &
                          & MAX(0.,MIN(ZRMAX,(PSFRV(:)/SQRT(ZTKEM_F(:,IKB)))*PARAMMF%XALP_PERT)) 

ZQT_UP(:) = PRT_UP(:,IKB)/(1.+PRT_UP(:,IKB))
ZTHS_UP(:,IKB)=PTHL_UP(:,IKB)*(1.+PARAMMF%XLAMBDA_MF*ZQT_UP(:))
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

CALL MZM_MF(D, PTHM (:,:), ZTHM_F(:,:))
CALL MZM_MF(D, PPABSM(:,:), ZPRES_F(:,:))
CALL MZM_MF(D, PRHODREF(:,:), ZRHO_F(:,:))
CALL MZM_MF(D, PRVM(:,:), ZRVM_F(:,:))

!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
! thetav at mass and flux levels 
ZTHVM_F(:,:)=ZTHM_F(:,:)*((1.+ZRVORD*ZRVM_F(:,:))/&
                                                             &(1.+ZRTM_F(:,:)))
ZTHVM(:,:)=PTHM(:,:)*((1.+ZRVORD*PRVM(:,:))/(1.+PRTM(:,:)))

PTHV_UP(:,:)= ZTHVM_F(:,:)
PRV_UP(:,:) = ZRVM_F(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

ZW_UP2(:,:)=ZEPS
!$mnh_expand_array(JIJ=IIJB:IIJE)
ZW_UP2(:,IKB) = MAX(0.0001,(1./6.)*ZTKEM_F(:,IKB))
GTEST(:) = (ZW_UP2(:,IKB) > ZEPS)  
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

! Computation of non conservative variable for the KKB level of the updraft
! (all or nothing ajustement)
!$mnh_expand_array(JIJ=IIJB:IIJE)
PRC_UP(:,IKB)=0.
PRI_UP(:,IKB)=0.
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

CALL TH_R_FROM_THL_RT(CST, NEB, D%NIJT, HFRAC_ICE,PFRAC_ICE_UP(:,IKB),ZPRES_F(:,IKB), &
             PTHL_UP(:,IKB),PRT_UP(:,IKB),ZTH_UP(:,IKB), &
             PRV_UP(:,IKB),PRC_UP(:,IKB),PRI_UP(:,IKB),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.,&
             PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)

!$mnh_expand_array(JIJ=IIJB:IIJE)
! compute updraft thevav and buoyancy term at KKB level             
PTHV_UP(:,IKB) = ZTH_UP(:,IKB)*((1+ZRVORD*PRV_UP(:,IKB))/(1+PRT_UP(:,IKB))) 
! compute mean rsat in updraft
PRSAT_UP(:,IKB) = ZRSATW(:)*(1-PFRAC_ICE_UP(:,IKB)) + &
                            & ZRSATI(:)*PFRAC_ICE_UP(:,IKB)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

!Tout est commente pour tester dans un premier temps la s�paration en deux de la 
!  boucle verticale, une pour w et une pour PEMF                                                            
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZG_O_THVREF(:,:)=CST%XG/ZTHVM_F(:,:)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

!  Definition de l'alimentation au sens de la fermeture de Hourdin et al

ZALIM_STAR(:,:)   = 0.
ZALIM_STAR_TOT(:) = 0.    ! <== Normalization of ZALIM_STAR
IALIM(:)          = IKB   ! <== Top level of the alimentation layer

DO JK=IKB,IKE-IKL,IKL   !  Vertical loop
  !$mnh_expand_where(JIJ=IIJB:IIJE)
  ZZDZ(:,JK)   = MAX(ZEPS,PZZ(:,JK+IKL)-PZZ(:,JK))       ! <== Delta Z between two flux level
  ZZZ(:,JK)    = MAX(0.,0.5*(PZZ(:,JK+IKL)+PZZ(:,JK)) )  ! <== Hight of mass levels
  ZDTHETASDZ(:,JK) = (ZTHVM_F(:,JK)-ZTHVM_F(:,JK+IKL))   ! <== Delta theta_v
  
  WHERE ((ZTHVM_F(:,JK+IKL)<ZTHVM_F(:,JK)) .AND. &
        &(ZTHVM_F(:,IKB)>=ZTHVM_F(:,JK)))
     ZALIM_STAR(:,JK)  = SQRT(ZZZ(:,JK))*ZDTHETASDZ(:,JK)/ZZDZ(:,JK)
     ZALIM_STAR_TOT(:) = ZALIM_STAR_TOT(:)+ZALIM_STAR(:,JK)*ZZDZ(:,JK)
     IALIM(:)          = JK
  ENDWHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)
ENDDO

! Normalization of ZALIM_STAR
DO JK=IKB,IKE-IKL,IKL   !  Vertical loop   
   !$mnh_expand_where(JIJ=IIJB:IIJE)
   WHERE (ZALIM_STAR_TOT(:) > ZEPS)
      ZALIM_STAR(:,JK)  = ZALIM_STAR(:,JK)/ZALIM_STAR_TOT(:)
   ENDWHERE
   !$mnh_end_expand_where(JIJ=IIJB:IIJE)
ENDDO
ZALIM_STAR_TOT(:) = 0.


! --------- END of alimentation calculation  ---------------------------------------      
  

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


DO JK=IKB,IKE-IKL,IKL
  !$mnh_expand_where(JIJ=IIJB:IIJE)
  ! IF the updraft top is reached for all column, stop the loop on levels

  !ITEST=COUNT(GTEST(:))
  !IF (ITEST==0) CYCLE

  !       Computation of entrainment and detrainment with KF90
  !       parameterization in clouds and LR01 in subcloud layer

  ! to find the LCL (check if JK is LCL or not)
  WHERE ((PRC_UP(:,JK)+PRI_UP(:,JK)>0.).AND.(.NOT.(GTESTLCL(:))))
    KKLCL(:) = JK           
    GTESTLCL(:)=.TRUE.
  ENDWHERE

  ! COMPUTE PENTR and PDETR at mass level JK

    
  !  Buoyancy is computed on "flux" levels where updraft variables are known

  ! Compute theta_v of updraft at flux level JK    
    
  ZRC_UP(:)  = PRC_UP(:,JK)
  ZRI_UP(:)  = PRI_UP(:,JK) ! guess 
  ZRV_UP(:)  = PRV_UP(:,JK)
  ZBUO(:,JK) = ZG_O_THVREF(:,JK)*(PTHV_UP(:,JK) - ZTHVM_F(:,JK))
  PBUO_INTEG(:,JK) = ZBUO(:,JK)*(PZZ(:,JK+IKL)-PZZ(:,JK))
  
  ZDZ(:)   = MAX(ZEPS,PZZ(:,JK+IKL)-PZZ(:,JK))
  ZTEST(:) = PARAMMF%XA1*ZBUO(:,JK) -  PARAMMF%XB*ZW_UP2(:,JK)

  ZCOE(:)      = ZDZ(:)
  WHERE (ZTEST(:)>0.)
    ZCOE(:)    = ZDZ(:)/(1.+ PARAMMF%XBETA1)
  ENDWHERE   

  !  Calcul de la vitesse

  ZWCOE(:)         = (1.-PARAMMF%XB*ZCOE(:))/(1.+PARAMMF%XB*ZCOE(:))
  ZBUCOE(:)        =  2.*ZCOE(:)/(1.+PARAMMF%XB*ZCOE(:))

  ZW_UP2(:,JK+IKL) = MAX(ZEPS,ZW_UP2(:,JK)*ZWCOE(:) + &
                                         &PARAMMF%XA1*ZBUO(:,JK)*ZBUCOE(:))
  ZW_MAX(:) = MAX(ZW_MAX(:), SQRT(ZW_UP2(:,JK+IKL)))
  ZWUP_MEAN(:)     = MAX(ZEPS,0.5*(ZW_UP2(:,JK+IKL)+ZW_UP2(:,JK)))

  !  Entrainement et detrainement

  PENTR(:,JK)  = MAX(0.,(PARAMMF%XBETA1/(1.+PARAMMF%XBETA1))* &
                                 &(PARAMMF%XA1*ZBUO(:,JK)/ZWUP_MEAN(:)-PARAMMF%XB))
  
  ZDETR_BUO(:) = MAX(0., -(PARAMMF%XBETA1/(1.+PARAMMF%XBETA1))*PARAMMF%XA1*ZBUO(:,JK)/ &
                                   &ZWUP_MEAN(:))
  ZDETR_RT(:)  = PARAMMF%XC*SQRT(MAX(0.,(PRT_UP(:,JK) - ZRTM_F(:,JK))) / &
                                          &MAX(ZEPS,ZRTM_F(:,JK)) / ZWUP_MEAN(:))
  PDETR(:,JK)  = ZDETR_RT(:)+ZDETR_BUO(:)

   
  ! If the updraft did not stop, compute cons updraft characteritics at jk+1
  WHERE(GTEST(:))     
    ZZTOP(:) = MAX(ZZTOP(:),PZZ(:,JK+IKL))
    ZMIX2(:) = (PZZ(:,JK+IKL)-PZZ(:,JK))*PENTR(:,JK) !&
    ZMIX3(:) = (PZZ(:,JK+IKL)-PZZ(:,JK))*PDETR(:,JK) !&                   
           
    ZQTM(:) = PRTM(:,JK)/(1.+PRTM(:,JK))            
    ZTHSM(:,JK) = PTHLM(:,JK)*(1.+PARAMMF%XLAMBDA_MF*ZQTM(:))
    ZTHS_UP(:,JK+IKL)=(ZTHS_UP(:,JK)*(1.-0.5*ZMIX2(:)) + &
                                    &ZTHSM(:,JK)*ZMIX2(:))&
                          /(1.+0.5*ZMIX2(:))   
    PRT_UP(:,JK+IKL)=(PRT_UP(:,JK)*(1.-0.5*ZMIX2(:)) + &
                                   &PRTM(:,JK)*ZMIX2(:))  &
                          /(1.+0.5*ZMIX2(:))
    ZQT_UP(:) = PRT_UP(:,JK+IKL)/(1.+PRT_UP(:,JK+IKL))
    PTHL_UP(:,JK+IKL)=ZTHS_UP(:,JK+IKL)/(1.+PARAMMF%XLAMBDA_MF*ZQT_UP(:))                      
  ENDWHERE
  

  IF(PARAMMF%LMIXUV) THEN
    IF(JK/=IKB) THEN
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
    ELSE
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

    ENDIF
  ENDIF
  !DO JSV=1,ISV 
  !   IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  !    WHERE(GTEST(:)) 
  !         PSV_UP(:,JK+KKL,JSV) = (PSV_UP(:,JK,JSV)*(1-0.5*ZMIX2(:)) + &
  !                      PSVM(:,JK,JSV)*ZMIX2(:))  /(1+0.5*ZMIX2(:))
  !    ENDWHERE                        
  !ENDDO  
  
  
  ! Compute non cons. var. at level JK+KKL
  ZRC_UP(:)=PRC_UP(:,JK) ! guess = level just below
  ZRI_UP(:)=PRI_UP(:,JK) ! guess = level just below
  ZRV_UP(:)=PRV_UP(:,JK)
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)
  CALL TH_R_FROM_THL_RT(CST,NEB, D%NIJT, HFRAC_ICE,PFRAC_ICE_UP(:,JK+IKL),ZPRES_F(:,JK+IKL), &
          PTHL_UP(:,JK+IKL),PRT_UP(:,JK+IKL),ZTH_UP(:,JK+IKL),              &
          ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.,&
          PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)
  !$mnh_expand_where(JIJ=IIJB:IIJE)
  WHERE(GTEST(:))
    ZT_UP(:) = ZTH_UP(:,JK+IKL)*PEXNM(:,JK+IKL)
    ZCP(:) = CST%XCPD + CST%XCL * ZRC_UP(:)
    ZLVOCPEXN(:)=(CST%XLVTT + (CST%XCPV-CST%XCL) *  (ZT_UP(:)-CST%XTT) ) / &
                          &ZCP(:) / PEXNM(:,JK+IKL)
    PRC_UP(:,JK+IKL)=MIN(0.5E-3,ZRC_UP(:))  ! On ne peut depasser 0.5 g/kg (autoconversion donc elimination !)
    PTHL_UP(:,JK+IKL) = PTHL_UP(:,JK+IKL)+ &
                                  & ZLVOCPEXN(:)*(ZRC_UP(:)-PRC_UP(:,JK+IKL))
    PRV_UP(:,JK+IKL)=ZRV_UP(:)
    PRI_UP(:,JK+IKL)=ZRI_UP(:)
    PRT_UP(:,JK+IKL)  = PRC_UP(:,JK+IKL) + PRV_UP(:,JK+IKL)
    PRSAT_UP(:,JK+IKL) = ZRSATW(:)*(1-PFRAC_ICE_UP(:,JK+IKL)) + &
                                   & ZRSATI(:)*PFRAC_ICE_UP(:,JK+IKL)
  ENDWHERE
  

  ! Compute the updraft theta_v, buoyancy and w**2 for level JK+1   
  WHERE(GTEST(:))
    !PTHV_UP(:,JK+KKL) = ZTH_UP(:,JK+KKL)*((1+ZRVORD*PRV_UP(:,JK+KKL))/(1+PRT_UP(:,JK+KKL)))
    PTHV_UP(:,JK+IKL) = ZTH_UP(:,JK+IKL)* &
                                  & (1.+0.608*PRV_UP(:,JK+IKL) - PRC_UP(:,JK+IKL))
  ENDWHERE


  ! Test if the updraft has reach the ETL
  GTESTETL(:)=.FALSE.
  WHERE (GTEST(:).AND.(PBUO_INTEG(:,JK)<=0.))
    KKETL(:) = JK+IKL
    GTESTETL(:)=.TRUE.
  ENDWHERE

  ! Test is we have reached the top of the updraft
  WHERE (GTEST(:).AND.((ZW_UP2(:,JK+IKL)<=ZEPS)))
    ZW_UP2(:,JK+IKL)=ZEPS
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
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)
ENDDO 

! Closure assumption for mass flux at KKB+1 level (Mass flux is supposed to be 0 at KKB level !)
!   Hourdin et al 2002 formulation


!$mnh_expand_array(JIJ=IIJB:IIJE)
ZZTOP(:) = MAX(ZZTOP(:),ZEPS)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)

DO JK=IKB+IKL,IKE-IKL,IKL !  Vertical loop
  !$mnh_expand_where(JIJ=IIJB:IIJE)
   WHERE(JK<=IALIM(:))
     ZALIM_STAR_TOT(:) = ZALIM_STAR_TOT(:) + ZALIM_STAR(:,JK)**2* &
                                                               & ZZDZ(:,JK)/PRHODREF(:,JK)
   ENDWHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)
ENDDO   

!$mnh_expand_where(JIJ=IIJB:IIJE)
WHERE (ZALIM_STAR_TOT(:)*ZZTOP(:) > ZEPS)
 ZPHI(:) =  ZW_MAX(:)/(PARAMMF%XR*ZZTOP(:)*ZALIM_STAR_TOT(:))
ENDWHERE

GTEST(:) = .TRUE.
PEMF(:,IKB+IKL) = ZPHI(:)*ZZDZ(:,IKB)*ZALIM_STAR(:,IKB)
! Updraft fraction must be smaller than XFRAC_UP_MAX
PFRAC_UP(:,IKB+IKL)=PEMF(:,IKB+IKL)/ &
                                 &(SQRT(ZW_UP2(:,IKB+IKL))*ZRHO_F(:,IKB+IKL))
PFRAC_UP(:,IKB+IKL)=MIN(PARAMMF%XFRAC_UP_MAX,PFRAC_UP(:,IKB+IKL))
PEMF(:,IKB+IKL) = ZRHO_F(:,IKB+IKL)*PFRAC_UP(:,IKB+IKL)* &
                              & SQRT(ZW_UP2(:,IKB+IKL))
!$mnh_end_expand_where(JIJ=IIJB:IIJE)

DO JK=IKB+IKL,IKE-IKL,IKL !  Vertical loop
  !$mnh_expand_where(JIJ=IIJB:IIJE)
     
  GTEST(:) = (ZW_UP2(:,JK) > ZEPS)  

  WHERE (GTEST(:))
    WHERE(JK<IALIM(:))
      PEMF(:,JK+IKL) = MAX(0.,PEMF(:,JK) + ZPHI(:)*ZZDZ(:,JK)* &
                                                                & (PENTR(:,JK) - PDETR(:,JK)))
    ELSEWHERE
      ZMIX1(:)=ZZDZ(:,JK)*(PENTR(:,JK)-PDETR(:,JK))
      PEMF(:,JK+IKL)=PEMF(:,JK)*EXP(ZMIX1(:))
    ENDWHERE

! Updraft fraction must be smaller than XFRAC_UP_MAX
    PFRAC_UP(:,JK+IKL)=PEMF(:,JK+IKL)/&
                                    &(SQRT(ZW_UP2(:,JK+IKL))*ZRHO_F(:,JK+IKL))
    PFRAC_UP(:,JK+IKL)=MIN(PARAMMF%XFRAC_UP_MAX,PFRAC_UP(:,JK+IKL))
    PEMF(:,JK+IKL) = ZRHO_F(:,JK+IKL)*PFRAC_UP(:,JK+IKL)*&
                                 & SQRT(ZW_UP2(:,JK+IKL))
  ENDWHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)
ENDDO

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
END DO

!$mnh_expand_array(JIJ=IIJB:IIJE)
GWORK1(:)= (GTESTLCL(:) .AND. (PDEPTH(:) > ZDEPTH_MAX1) )
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
DO JK=1,D%NKT
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

IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAF_RAHA',1,ZHOOK_HANDLE)
!
CONTAINS
INCLUDE "th_r_from_thl_rt.func.h"
INCLUDE "compute_frac_ice.func.h"
!
END SUBROUTINE COMPUTE_UPDRAFT_RAHA
END MODULE MODE_COMPUTE_UPDRAFT_RAHA
