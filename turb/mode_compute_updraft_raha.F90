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
                                 KSV, HFRAC_ICE, OENTR_DETR, OMIXUV, &
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

INTEGER  :: JK,JI,JJ,JSV          ! loop counters

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

PW_UP(:,:)=0.
ZTH_UP(:,:)=0.
PFRAC_UP(:,:)=0.
PTHV_UP(:,:)=0.

PBUO_INTEG(:,:)=0.
ZBUO(:,:)      =0.

!no ice cloud coded yet 
PRI_UP(:,:)=0.
PFRAC_ICE_UP(:,:)=0.
DO JK=1,D%NKT 
  DO JI=D%NIJB,D%NIJE 
    PRSAT_UP(JI,JK)=PRVM(JI,JK) ! should be initialised correctly but is (normaly) not used
  ENDDO
ENDDO

! Initialisation of environment variables at t-dt

! variables at flux level
CALL MZM_MF(D, PTHLM(:,:), ZTHLM_F(:,:))
CALL MZM_MF(D, PRTM(:,:), ZRTM_F(:,:))
CALL MZM_MF(D, PUM(:,:), ZUM_F(:,:))
CALL MZM_MF(D, PVM(:,:), ZVM_F(:,:))
CALL MZM_MF(D, PTKEM(:,:), ZTKEM_F(:,:))

!DO JSV=1,ISV 
! IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
!   ZSVM_F(D%NIJB:D%NIJE,KKB:IKU,JSV) = 0.5*(PSVM(D%NIJB:D%NIJE,KKB:IKU,JSV)+PSVM(D%NIJB:D%NIJE,1:IKU-1,JSV))
!   ZSVM_F(D%NIJB:D%NIJE,1,JSV)       = ZSVM_F(D%NIJB:D%NIJE,KKB,JSV) 
!END DO  

!          Initialisation of updraft characteristics 
DO JK=1,D%NKT 
  DO JI=D%NIJB,D%NIJE 
    PTHL_UP(JI,JK)=ZTHLM_F(JI,JK)
    PRT_UP(JI,JK)=ZRTM_F(JI,JK)
    PU_UP(JI,JK)=ZUM_F(JI,JK)
    PV_UP(JI,JK)=ZVM_F(JI,JK)
  ENDDO
ENDDO
PSV_UP(:,:,:)=0.
!IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) then
!    PSV_UP(D%NIJB:D%NIJE,:,:)=ZSVM_F(D%NIJB:D%NIJE,:,:)
!ENDIF

! Computation or initialisation of updraft characteristics at the KKB level
! thetal_up,rt_up,thetaV_up, w�,Buoyancy term and mass flux (PEMF)

DO JI=D%NIJB,D%NIJE 
  PTHL_UP(JI,D%NKB)= ZTHLM_F(JI,D%NKB)+ &
  & MAX(0.,MIN(ZTMAX,(PSFTH(JI)/SQRT(ZTKEM_F(JI,D%NKB)))*PARAMMF%XALP_PERT))
  PRT_UP(JI,D%NKB) = ZRTM_F(JI,D%NKB)+ &
  & MAX(0.,MIN(ZRMAX,(PSFRV(JI)/SQRT(ZTKEM_F(JI,D%NKB)))*PARAMMF%XALP_PERT)) 
  
  ZQT_UP(JI) = PRT_UP(JI,D%NKB)/(1.+PRT_UP(JI,D%NKB))
  ZTHS_UP(JI,D%NKB)=PTHL_UP(JI,D%NKB)*(1.+PARAMMF%XLAMBDA_MF*ZQT_UP(JI))
ENDDO

CALL MZM_MF(D, PTHM (:,:), ZTHM_F(:,:))
CALL MZM_MF(D, PPABSM(:,:), ZPRES_F(:,:))
CALL MZM_MF(D, PRHODREF(:,:), ZRHO_F(:,:))
CALL MZM_MF(D, PRVM(:,:), ZRVM_F(:,:))

DO JK=1,D%NKT 
  DO JI=D%NIJB,D%NIJE 
! thetav at mass and flux levels 
    ZTHVM_F(JI,JK)=ZTHM_F(JI,JK)*((1.+ZRVORD*ZRVM_F(JI,JK))/&
    &(1.+ZRTM_F(JI,JK)))
    ZTHVM(JI,JK)=PTHM(JI,JK)*((1.+ZRVORD*PRVM(JI,JK))/(1.+PRTM(JI,JK)))
    
    PTHV_UP(JI,JK)= ZTHVM_F(JI,JK)
    PRV_UP(JI,JK) = ZRVM_F(JI,JK)
  ENDDO
ENDDO

ZW_UP2(:,:)=ZEPS
DO JI=D%NIJB,D%NIJE 
  ZW_UP2(JI,D%NKB) = MAX(0.0001,(1./6.)*ZTKEM_F(JI,D%NKB))
  GTEST(JI) = (ZW_UP2(JI,D%NKB) > ZEPS)  
ENDDO

! Computation of non conservative variable for the KKB level of the updraft
! (all or nothing ajustement)
DO JI=D%NIJB,D%NIJE 
  PRC_UP(JI,D%NKB)=0.
  PRI_UP(JI,D%NKB)=0.
ENDDO

CALL TH_R_FROM_THL_RT(CST, NEB, D%NIJT, HFRAC_ICE,PFRAC_ICE_UP(:,D%NKB),ZPRES_F(:,D%NKB), &
             PTHL_UP(:,D%NKB),PRT_UP(:,D%NKB),ZTH_UP(:,D%NKB), &
             PRV_UP(:,D%NKB),PRC_UP(:,D%NKB),PRI_UP(:,D%NKB),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.,&
             PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)

DO JI=D%NIJB,D%NIJE 
! compute updraft thevav and buoyancy term at KKB level             
  PTHV_UP(JI,D%NKB) = ZTH_UP(JI,D%NKB)*((1+ZRVORD*PRV_UP(JI,D%NKB))/(1+PRT_UP(JI,D%NKB))) 
! compute mean rsat in updraft
  PRSAT_UP(JI,D%NKB) = ZRSATW(JI)*(1-PFRAC_ICE_UP(JI,D%NKB)) + &
  & ZRSATI(JI)*PFRAC_ICE_UP(JI,D%NKB)
ENDDO

!Tout est commente pour tester dans un premier temps la s�paration en deux de la 
!  boucle verticale, une pour w et une pour PEMF                                                            
DO JK=1,D%NKT 
  DO JI=D%NIJB,D%NIJE 
    ZG_O_THVREF(JI,JK)=CST%XG/ZTHVM_F(JI,JK)
  ENDDO
ENDDO

!  Definition de l'alimentation au sens de la fermeture de Hourdin et al

ZALIM_STAR(:,:)   = 0.
ZALIM_STAR_TOT(:) = 0.    ! <== Normalization of ZALIM_STAR
IALIM(:)          = D%NKB   ! <== Top level of the alimentation layer

DO JK=D%NKB,D%NKE-D%NKL,D%NKL   !  Vertical loop
  DO JI=D%NIJB,D%NIJE 
    ZZDZ(JI,JK)   = MAX(ZEPS,PZZ(JI,JK+D%NKL)-PZZ(JI,JK))       ! <== Delta Z between two flux level
    ZZZ(JI,JK)    = MAX(0.,0.5*(PZZ(JI,JK+D%NKL)+PZZ(JI,JK)) )  ! <== Hight of mass levels
    ZDTHETASDZ(JI,JK) = (ZTHVM_F(JI,JK)-ZTHVM_F(JI,JK+D%NKL))   ! <== Delta theta_v
    
    IF ((ZTHVM_F(JI,JK+D%NKL)<ZTHVM_F(JI,JK)) .AND. &
        &(ZTHVM_F(JI,D%NKB)>=ZTHVM_F(JI,JK)))THEN
    ZALIM_STAR(JI,JK)  = SQRT(ZZZ(JI,JK))*ZDTHETASDZ(JI,JK)/ZZDZ(JI,JK)
    ZALIM_STAR_TOT(JI) = ZALIM_STAR_TOT(JI)+ZALIM_STAR(JI,JK)*ZZDZ(JI,JK)
    IALIM(JI)          = JK
  ENDIF
ENDDO
ENDDO

! Normalization of ZALIM_STAR
DO JK=D%NKB,D%NKE-D%NKL,D%NKL   !  Vertical loop   
DO JI=D%NIJB,D%NIJE 
  IF (ZALIM_STAR_TOT(JI) > ZEPS)THEN
    ZALIM_STAR(JI,JK)  = ZALIM_STAR(JI,JK)/ZALIM_STAR_TOT(JI)
  ENDIF
ENDDO
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


DO JK=D%NKB,D%NKE-D%NKL,D%NKL
DO JI=D%NIJB,D%NIJE 
  ! IF the updraft top is reached for all column, stop the loop on levels
  
  !ITEST=COUNT(GTEST(JI))
  !IF (ITEST==0) CYCLE
  
  !       Computation of entrainment and detrainment with KF90
  !       parameterization in clouds and LR01 in subcloud layer
  
  ! to find the LCL (check if JK is LCL or not)
  IF ((PRC_UP(JI,JK)+PRI_UP(JI,JK)>0.).AND.(.NOT.(GTESTLCL(JI))))THEN
    KKLCL(JI) = JK           
    GTESTLCL(JI)=.TRUE.
  ENDIF
  
  ! COMPUTE PENTR and PDETR at mass level JK
  
  
  !  Buoyancy is computed on "flux" levels where updraft variables are known
  
  ! Compute theta_v of updraft at flux level JK    
  
  ZRC_UP(JI)  = PRC_UP(JI,JK)
  ZRI_UP(JI)  = PRI_UP(JI,JK) ! guess 
  ZRV_UP(JI)  = PRV_UP(JI,JK)
  ZBUO(JI,JK) = ZG_O_THVREF(JI,JK)*(PTHV_UP(JI,JK) - ZTHVM_F(JI,JK))
  PBUO_INTEG(JI,JK) = ZBUO(JI,JK)*(PZZ(JI,JK+D%NKL)-PZZ(JI,JK))
  
  ZDZ(JI)   = MAX(ZEPS,PZZ(JI,JK+D%NKL)-PZZ(JI,JK))
  ZTEST(JI) = PARAMMF%XA1*ZBUO(JI,JK) -  PARAMMF%XB*ZW_UP2(JI,JK)
  
  ZCOE(JI)      = ZDZ(JI)
  IF (ZTEST(JI)>0.)THEN
    ZCOE(JI)    = ZDZ(JI)/(1.+ PARAMMF%XBETA1)
  ENDIF   
  
  !  Calcul de la vitesse
  
  ZWCOE(JI)         = (1.-PARAMMF%XB*ZCOE(JI))/(1.+PARAMMF%XB*ZCOE(JI))
  ZBUCOE(JI)        =  2.*ZCOE(JI)/(1.+PARAMMF%XB*ZCOE(JI))
  
  ZW_UP2(JI,JK+D%NKL) = MAX(ZEPS,ZW_UP2(JI,JK)*ZWCOE(JI) + &
  &PARAMMF%XA1*ZBUO(JI,JK)*ZBUCOE(JI))
  ZW_MAX(JI) = MAX(ZW_MAX(JI), SQRT(ZW_UP2(JI,JK+D%NKL)))
  ZWUP_MEAN(JI)     = MAX(ZEPS,0.5*(ZW_UP2(JI,JK+D%NKL)+ZW_UP2(JI,JK)))
  
  !  Entrainement et detrainement
  
  PENTR(JI,JK)  = MAX(0.,(PARAMMF%XBETA1/(1.+PARAMMF%XBETA1))* &
  &(PARAMMF%XA1*ZBUO(JI,JK)/ZWUP_MEAN(JI)-PARAMMF%XB))
  
  ZDETR_BUO(JI) = MAX(0., -(PARAMMF%XBETA1/(1.+PARAMMF%XBETA1))*PARAMMF%XA1*ZBUO(JI,JK)/ &
  &ZWUP_MEAN(JI))
  ZDETR_RT(JI)  = PARAMMF%XC*SQRT(MAX(0.,(PRT_UP(JI,JK) - ZRTM_F(JI,JK))) / &
  &MAX(ZEPS,ZRTM_F(JI,JK)) / ZWUP_MEAN(JI))
  PDETR(JI,JK)  = ZDETR_RT(JI)+ZDETR_BUO(JI)
  
  
  ! If the updraft did not stop, compute cons updraft characteritics at jk+1
  IF(GTEST(JI))     THEN
    ZZTOP(JI) = MAX(ZZTOP(JI),PZZ(JI,JK+D%NKL))
    ZMIX2(JI) = (PZZ(JI,JK+D%NKL)-PZZ(JI,JK))*PENTR(JI,JK) !&
    ZMIX3(JI) = (PZZ(JI,JK+D%NKL)-PZZ(JI,JK))*PDETR(JI,JK) !&                   
    
    ZQTM(JI) = PRTM(JI,JK)/(1.+PRTM(JI,JK))            
    ZTHSM(JI,JK) = PTHLM(JI,JK)*(1.+PARAMMF%XLAMBDA_MF*ZQTM(JI))
    ZTHS_UP(JI,JK+D%NKL)=(ZTHS_UP(JI,JK)*(1.-0.5*ZMIX2(JI)) + &
    &ZTHSM(JI,JK)*ZMIX2(JI))&
    /(1.+0.5*ZMIX2(JI))   
    PRT_UP(JI,JK+D%NKL)=(PRT_UP(JI,JK)*(1.-0.5*ZMIX2(JI)) + &
    &PRTM(JI,JK)*ZMIX2(JI))  &
    /(1.+0.5*ZMIX2(JI))
    ZQT_UP(JI) = PRT_UP(JI,JK+D%NKL)/(1.+PRT_UP(JI,JK+D%NKL))
    PTHL_UP(JI,JK+D%NKL)=ZTHS_UP(JI,JK+D%NKL)/(1.+PARAMMF%XLAMBDA_MF*ZQT_UP(JI))                      
  ENDIF
  
  
  IF(OMIXUV) THEN
    IF(JK/=D%NKB) THEN
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
    ELSE
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
      
    ENDIF
  ENDIF
  !DO JSV=1,ISV 
  !   IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  !    IF(GTEST(JI)) THEN
  !         PSV_UP(JI,JK+KKL,JSV) = (PSV_UP(JI,JK,JSV)*(1-0.5*ZMIX2(JI)) + &
  !                      PSVM(JI,JK,JSV)*ZMIX2(JI))  /(1+0.5*ZMIX2(JI))
  !    ENDIF                        
  !ENDDO  
  
  
  ! Compute non cons. var. at level JK+KKL
  ZRC_UP(JI)=PRC_UP(JI,JK) ! guess = level just below
  ZRI_UP(JI)=PRI_UP(JI,JK) ! guess = level just below
  ZRV_UP(JI)=PRV_UP(JI,JK)
ENDDO
  CALL TH_R_FROM_THL_RT(CST,NEB, D%NIJT, HFRAC_ICE,PFRAC_ICE_UP(:,JK+D%NKL),ZPRES_F(:,JK+D%NKL), &
          PTHL_UP(:,JK+D%NKL),PRT_UP(:,JK+D%NKL),ZTH_UP(:,JK+D%NKL),              &
          ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.,&
          PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)
DO JI=D%NIJB,D%NIJE 
  IF(GTEST(JI))THEN
    ZT_UP(JI) = ZTH_UP(JI,JK+D%NKL)*PEXNM(JI,JK+D%NKL)
    ZCP(JI) = CST%XCPD + CST%XCL * ZRC_UP(JI)
    ZLVOCPEXN(JI)=(CST%XLVTT + (CST%XCPV-CST%XCL) *  (ZT_UP(JI)-CST%XTT) ) / &
    &ZCP(JI) / PEXNM(JI,JK+D%NKL)
    PRC_UP(JI,JK+D%NKL)=MIN(0.5E-3,ZRC_UP(JI))  ! On ne peut depasser 0.5 g/kg (autoconversion donc elimination !)
    PTHL_UP(JI,JK+D%NKL) = PTHL_UP(JI,JK+D%NKL)+ &
    & ZLVOCPEXN(JI)*(ZRC_UP(JI)-PRC_UP(JI,JK+D%NKL))
    PRV_UP(JI,JK+D%NKL)=ZRV_UP(JI)
    PRI_UP(JI,JK+D%NKL)=ZRI_UP(JI)
    PRT_UP(JI,JK+D%NKL)  = PRC_UP(JI,JK+D%NKL) + PRV_UP(JI,JK+D%NKL)
    PRSAT_UP(JI,JK+D%NKL) = ZRSATW(JI)*(1-PFRAC_ICE_UP(JI,JK+D%NKL)) + &
    & ZRSATI(JI)*PFRAC_ICE_UP(JI,JK+D%NKL)
  ENDIF
  
  
  ! Compute the updraft theta_v, buoyancy and w**2 for level JK+1   
  IF(GTEST(JI))THEN
    !PTHV_UP(JI,JK+KKL) = ZTH_UP(JI,JK+KKL)*((1+ZRVORD*PRV_UP(JI,JK+KKL))/(1+PRT_UP(JI,JK+KKL)))
    PTHV_UP(JI,JK+D%NKL) = ZTH_UP(JI,JK+D%NKL)* &
    & (1.+0.608*PRV_UP(JI,JK+D%NKL) - PRC_UP(JI,JK+D%NKL))
  ENDIF
  
  
  ! Test if the updraft has reach the ETL
  GTESTETL(JI)=.FALSE.
  IF (GTEST(JI).AND.(PBUO_INTEG(JI,JK)<=0.))THEN
    KKETL(JI) = JK+D%NKL
    GTESTETL(JI)=.TRUE.
  ENDIF
  
  ! Test is we have reached the top of the updraft
  IF (GTEST(JI).AND.((ZW_UP2(JI,JK+D%NKL)<=ZEPS)))THEN
    ZW_UP2(JI,JK+D%NKL)=ZEPS
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
ENDDO
ENDDO 

! Closure assumption for mass flux at KKB+1 level (Mass flux is supposed to be 0 at KKB level !)
!   Hourdin et al 2002 formulation


DO JI=D%NIJB,D%NIJE 
ZZTOP(JI) = MAX(ZZTOP(JI),ZEPS)
ENDDO

DO JK=D%NKB+D%NKL,D%NKE-D%NKL,D%NKL !  Vertical loop
DO JI=D%NIJB,D%NIJE 
  IF(JK<=IALIM(JI))THEN
    ZALIM_STAR_TOT(JI) = ZALIM_STAR_TOT(JI) + ZALIM_STAR(JI,JK)**2* &
    & ZZDZ(JI,JK)/PRHODREF(JI,JK)
  ENDIF
ENDDO
ENDDO   

DO JI=D%NIJB,D%NIJE 
IF (ZALIM_STAR_TOT(JI)*ZZTOP(JI) > ZEPS)THEN
  ZPHI(JI) =  ZW_MAX(JI)/(PARAMMF%XR*ZZTOP(JI)*ZALIM_STAR_TOT(JI))
ENDIF

GTEST(JI) = .TRUE.
PEMF(JI,D%NKB+D%NKL) = ZPHI(JI)*ZZDZ(JI,D%NKB)*ZALIM_STAR(JI,D%NKB)
! Updraft fraction must be smaller than XFRAC_UP_MAX
PFRAC_UP(JI,D%NKB+D%NKL)=PEMF(JI,D%NKB+D%NKL)/ &
                                 &(SQRT(ZW_UP2(JI,D%NKB+D%NKL))*ZRHO_F(JI,D%NKB+D%NKL))
PFRAC_UP(JI,D%NKB+D%NKL)=MIN(PARAMMF%XFRAC_UP_MAX,PFRAC_UP(JI,D%NKB+D%NKL))
PEMF(JI,D%NKB+D%NKL) = ZRHO_F(JI,D%NKB+D%NKL)*PFRAC_UP(JI,D%NKB+D%NKL)* &
                              & SQRT(ZW_UP2(JI,D%NKB+D%NKL))
ENDDO

DO JK=D%NKB+D%NKL,D%NKE-D%NKL,D%NKL !  Vertical loop
DO JI=D%NIJB,D%NIJE 
  
  GTEST(JI) = (ZW_UP2(JI,JK) > ZEPS)  
  
  IF (GTEST(JI))THEN
    IF(JK<IALIM(JI))THEN
      PEMF(JI,JK+D%NKL) = MAX(0.,PEMF(JI,JK) + ZPHI(JI)*ZZDZ(JI,JK)* &
      & (PENTR(JI,JK) - PDETR(JI,JK)))
    ELSE
      ZMIX1(JI)=ZZDZ(JI,JK)*(PENTR(JI,JK)-PDETR(JI,JK))
      PEMF(JI,JK+D%NKL)=PEMF(JI,JK)*EXP(ZMIX1(JI))
    ENDIF
    
! Updraft fraction must be smaller than XFRAC_UP_MAX
    PFRAC_UP(JI,JK+D%NKL)=PEMF(JI,JK+D%NKL)/&
    &(SQRT(ZW_UP2(JI,JK+D%NKL))*ZRHO_F(JI,JK+D%NKL))
    PFRAC_UP(JI,JK+D%NKL)=MIN(PARAMMF%XFRAC_UP_MAX,PFRAC_UP(JI,JK+D%NKL))
    PEMF(JI,JK+D%NKL) = ZRHO_F(JI,JK+D%NKL)*PFRAC_UP(JI,JK+D%NKL)*&
    & SQRT(ZW_UP2(JI,JK+D%NKL))
  ENDIF
ENDDO
ENDDO

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
! from 1 (for clouds of 3000m of depth) to 0 (for clouds of 4000m of depth).
! This way, all MF fluxes are diminished by this amount.
! Diagnosed cloud fraction is also multiplied by the same coefficient.
!
DO JI=D%NIJB,D%NIJE 
  PDEPTH(JI) = MAX(0., PZZ(JI,KKCTL(JI)) -  PZZ(JI,KKLCL(JI)) )
END DO

DO JI=D%NIJB,D%NIJE 
GWORK1(JI)= (GTESTLCL(JI) .AND. (PDEPTH(JI) > ZDEPTH_MAX1) )
ENDDO
DO JK=1,D%NKT
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

IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAF_RAHA',1,ZHOOK_HANDLE)
!
CONTAINS
INCLUDE "th_r_from_thl_rt.func.h"
INCLUDE "compute_frac_ice.func.h"
!
END SUBROUTINE COMPUTE_UPDRAFT_RAHA
END MODULE MODE_COMPUTE_UPDRAFT_RAHA
