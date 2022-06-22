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

REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   ::  PEXNM       ! Exner function at t-dt
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   ::  PTHM           ! liquid pot. temp. at t-dt
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   ::  PRVM           ! vapor mixing ratio at t-dt
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   ::  PTHLM,PRTM     ! cons. var. at t-dt

REAL, DIMENSION(D%NIT,D%NKT,KSV), INTENT(IN)   ::  PSVM           ! scalar var. at t-dt

REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT)  ::  PTHL_UP,PRT_UP   ! updraft properties
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT)  ::  PU_UP, PV_UP     ! updraft wind components
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(INOUT)::  PRV_UP,PRC_UP    ! updraft rv, rc
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(INOUT)::  PRI_UP,PTHV_UP   ! updraft ri, THv
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(INOUT)::  PW_UP,PFRAC_UP   ! updraft w, fraction
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(INOUT)::  PFRAC_ICE_UP     ! liquid/solid fraction in updraft
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(INOUT)::  PRSAT_UP         ! Rsat

REAL, DIMENSION(D%NIT,D%NKT,KSV), INTENT(OUT)  ::  PSV_UP           ! updraft scalar var. 
                                         
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(INOUT)::  PEMF,PDETR,PENTR ! Mass_flux,
                                                          ! detrainment,entrainment
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(INOUT) :: PBUO_INTEG       ! Integrated Buoyancy 
INTEGER, DIMENSION(D%NIT),  INTENT(INOUT)::  KKLCL,KKETL,KKCTL! LCL, ETL, CTL                                     
REAL, DIMENSION(D%NIT),     INTENT(OUT)   :: PDEPTH           ! Deepness of cloud
!                       1.2  Declaration of local variables
!
!
! Mean environment variables at t-dt at flux point
REAL, DIMENSION(D%NIT,D%NKT) ::    ZTHM_F,ZRVM_F,ZRCM_F    ! Theta,rv of
                                                                  ! updraft environnement
REAL, DIMENSION(D%NIT,D%NKT) :: ZRTM_F, ZTHLM_F, ZTKEM_F      ! rt, thetal,TKE,pressure,
REAL, DIMENSION(D%NIT,D%NKT) :: ZUM_F,ZVM_F,ZRHO_F            ! density,momentum
REAL, DIMENSION(D%NIT,D%NKT) :: ZPRES_F,ZTHVM_F,ZTHVM         ! interpolated at the flux point
REAL, DIMENSION(D%NIT,D%NKT) :: ZG_O_THVREF                   ! g*ThetaV ref
REAL, DIMENSION(D%NIT,D%NKT) :: ZW_UP2                        ! w**2  of the updraft

REAL, DIMENSION(D%NIT,D%NKT,KSV) :: ZSVM_F ! scalar variables 
                        

                        
REAL, DIMENSION(D%NIT,D%NKT) :: ZTH_UP                        ! updraft THETA 
REAL, DIMENSION(D%NIT)              :: ZT_UP                         ! updraft T
REAL, DIMENSION(D%NIT)              :: ZLVOCPEXN                     ! updraft L
REAL, DIMENSION(D%NIT)              :: ZCP                           ! updraft cp
REAL, DIMENSION(D%NIT,D%NKT) :: ZBUO                          ! Buoyancy 
REAL, DIMENSION(D%NIT,D%NKT) :: ZTHS_UP,ZTHSM

REAL, DIMENSION(D%NIT,D%NKT) ::  ZCOEF  ! diminution coefficient for too high clouds 
                        
REAL, DIMENSION(D%NIT)            ::  ZWTHVSURF  ! Surface w'thetav'

REAL  :: ZRDORV       ! RD/RV
REAL  :: ZRVORD       ! RV/RD


REAL, DIMENSION(D%NIT) :: ZMIX1,ZMIX2,ZMIX3

REAL, DIMENSION(D%NIT) :: ZLUP         ! Upward Mixing length from the ground

REAL, DIMENSION(D%NIT) :: ZDEPTH       ! Deepness limit for cloud

INTEGER  :: JK,JI,JJ,JSV          ! loop counters

LOGICAL, DIMENSION(D%NIT) ::  GTEST,GTESTLCL,GTESTETL
                               ! Test if the ascent continue, if LCL or ETL is reached
LOGICAL                          ::  GLMIX 
                               ! To choose upward or downward mixing length
LOGICAL, DIMENSION(D%NIT)              :: GWORK1
LOGICAL, DIMENSION(D%NIT,D%NKT) :: GWORK2


INTEGER  :: ITEST

REAL, DIMENSION(D%NIT) :: ZRC_UP, ZRI_UP, ZRV_UP, ZWP2, ZRSATW, ZRSATI

LOGICAL, DIMENSION(D%NIT) :: GTEST_FER
REAL,    DIMENSION(D%NIT) :: ZPHI,ZALIM_STAR_TOT
REAL,    DIMENSION(D%NIT,D%NKT) :: ZDTHETASDZ,ZALIM_STAR,ZZDZ,ZZZ
INTEGER, DIMENSION(D%NIT) :: IALIM

REAL, DIMENSION(D%NIT)              ::  ZTEST,ZDZ,ZWUP_MEAN    ! 
REAL, DIMENSION(D%NIT)              ::  ZCOE,ZWCOE,ZBUCOE
REAL, DIMENSION(D%NIT)              ::  ZDETR_BUO, ZDETR_RT
REAL, DIMENSION(D%NIT)              ::  ZW_MAX               ! w**2  max of the updraft
REAL, DIMENSION(D%NIT)              ::  ZZTOP                ! Top of the updraft
REAL, DIMENSION(D%NIT)              ::  ZA,ZB,ZQTM,ZQT_UP

REAL  :: ZDEPTH_MAX1, ZDEPTH_MAX2 ! control auto-extinction process

REAL  :: ZTMAX,ZRMAX, ZEPS  ! control value
REAL, DIMENSION(D%NIT,16)           ::  ZBUF

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
!$mnh_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
PRSAT_UP(D%NIB:D%NIE,:)=PRVM(D%NIB:D%NIE,:) ! should be initialised correctly but is (normaly) not used
!$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)

! Initialisation of environment variables at t-dt

! variables at flux level
CALL MZM_MF(D, PTHLM(:,:), ZTHLM_F(:,:))
CALL MZM_MF(D, PRTM(:,:), ZRTM_F(:,:))
CALL MZM_MF(D, PUM(:,:), ZUM_F(:,:))
CALL MZM_MF(D, PVM(:,:), ZVM_F(:,:))
CALL MZM_MF(D, PTKEM(:,:), ZTKEM_F(:,:))

!DO JSV=1,ISV 
! IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
!   ZSVM_F(D%NIB:D%NIE,KKB:IKU,JSV) = 0.5*(PSVM(D%NIB:D%NIE,KKB:IKU,JSV)+PSVM(D%NIB:D%NIE,1:IKU-1,JSV))
!   ZSVM_F(D%NIB:D%NIE,1,JSV)       = ZSVM_F(D%NIB:D%NIE,KKB,JSV) 
!END DO  

!          Initialisation of updraft characteristics 
!$mnh_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
PTHL_UP(D%NIB:D%NIE,:)=ZTHLM_F(D%NIB:D%NIE,:)
PRT_UP(D%NIB:D%NIE,:)=ZRTM_F(D%NIB:D%NIE,:)
PU_UP(D%NIB:D%NIE,:)=ZUM_F(D%NIB:D%NIE,:)
PV_UP(D%NIB:D%NIE,:)=ZVM_F(D%NIB:D%NIE,:)
!$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
PSV_UP(:,:,:)=0.
!IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) then
!    PSV_UP(D%NIB:D%NIE,:,:)=ZSVM_F(D%NIB:D%NIE,:,:)
!ENDIF

! Computation or initialisation of updraft characteristics at the KKB level
! thetal_up,rt_up,thetaV_up, w�,Buoyancy term and mass flux (PEMF)

!$mnh_expand_array(JI=D%NIB:D%NIE)
PTHL_UP(D%NIB:D%NIE,D%NKB)= ZTHLM_F(D%NIB:D%NIE,D%NKB)+ &
                          & MAX(0.,MIN(ZTMAX,(PSFTH(D%NIB:D%NIE)/SQRT(ZTKEM_F(D%NIB:D%NIE,D%NKB)))*PARAMMF%XALP_PERT))
PRT_UP(D%NIB:D%NIE,D%NKB) = ZRTM_F(D%NIB:D%NIE,D%NKB)+ &
                          & MAX(0.,MIN(ZRMAX,(PSFRV(D%NIB:D%NIE)/SQRT(ZTKEM_F(D%NIB:D%NIE,D%NKB)))*PARAMMF%XALP_PERT)) 

ZQT_UP(D%NIB:D%NIE) = PRT_UP(D%NIB:D%NIE,D%NKB)/(1.+PRT_UP(D%NIB:D%NIE,D%NKB))
ZTHS_UP(D%NIB:D%NIE,D%NKB)=PTHL_UP(D%NIB:D%NIE,D%NKB)*(1.+PARAMMF%XLAMBDA_MF*ZQT_UP(D%NIB:D%NIE))
!$mnh_end_expand_array(JI=D%NIB:D%NIE)

CALL MZM_MF(D, PTHM (:,:), ZTHM_F(:,:))
CALL MZM_MF(D, PPABSM(:,:), ZPRES_F(:,:))
CALL MZM_MF(D, PRHODREF(:,:), ZRHO_F(:,:))
CALL MZM_MF(D, PRVM(:,:), ZRVM_F(:,:))

!$mnh_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
! thetav at mass and flux levels 
ZTHVM_F(D%NIB:D%NIE,:)=ZTHM_F(D%NIB:D%NIE,:)*((1.+ZRVORD*ZRVM_F(D%NIB:D%NIE,:))/(1.+ZRTM_F(D%NIB:D%NIE,:)))
ZTHVM(D%NIB:D%NIE,:)=PTHM(D%NIB:D%NIE,:)*((1.+ZRVORD*PRVM(D%NIB:D%NIE,:))/(1.+PRTM(D%NIB:D%NIE,:)))

PTHV_UP(D%NIB:D%NIE,:)= ZTHVM_F(D%NIB:D%NIE,:)
PRV_UP(D%NIB:D%NIE,:) = ZRVM_F(D%NIB:D%NIE,:)
!$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)

ZW_UP2(:,:)=ZEPS
!$mnh_expand_array(JI=D%NIB:D%NIE)
ZW_UP2(D%NIB:D%NIE,D%NKB) = MAX(0.0001,(1./6.)*ZTKEM_F(D%NIB:D%NIE,D%NKB))
GTEST(D%NIB:D%NIE) = (ZW_UP2(D%NIB:D%NIE,D%NKB) > ZEPS)  
!$mnh_end_expand_array(JI=D%NIB:D%NIE)

! Computation of non conservative variable for the KKB level of the updraft
! (all or nothing ajustement)
!$mnh_expand_array(JI=D%NIB:D%NIE)
PRC_UP(D%NIB:D%NIE,D%NKB)=0.
PRI_UP(D%NIB:D%NIE,D%NKB)=0.
!$mnh_end_expand_array(JI=D%NIB:D%NIE)

CALL TH_R_FROM_THL_RT(CST, NEB, D%NIT, HFRAC_ICE,PFRAC_ICE_UP(:,D%NKB),ZPRES_F(:,D%NKB), &
             PTHL_UP(:,D%NKB),PRT_UP(:,D%NKB),ZTH_UP(:,D%NKB), &
             PRV_UP(:,D%NKB),PRC_UP(:,D%NKB),PRI_UP(:,D%NKB),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.,&
             PBUF=ZBUF, KB=D%NIB, KE=D%NIE)

!$mnh_expand_array(JI=D%NIB:D%NIE)
! compute updraft thevav and buoyancy term at KKB level             
PTHV_UP(D%NIB:D%NIE,D%NKB) = ZTH_UP(D%NIB:D%NIE,D%NKB)*((1+ZRVORD*PRV_UP(D%NIB:D%NIE,D%NKB))/(1+PRT_UP(D%NIB:D%NIE,D%NKB))) 
! compute mean rsat in updraft
PRSAT_UP(D%NIB:D%NIE,D%NKB) = ZRSATW(D%NIB:D%NIE)*(1-PFRAC_ICE_UP(D%NIB:D%NIE,D%NKB)) + &
                            & ZRSATI(D%NIB:D%NIE)*PFRAC_ICE_UP(D%NIB:D%NIE,D%NKB)
!$mnh_end_expand_array(JI=D%NIB:D%NIE)

!Tout est commente pour tester dans un premier temps la s�paration en deux de la 
!  boucle verticale, une pour w et une pour PEMF                                                            
!$mnh_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
ZG_O_THVREF(D%NIB:D%NIE,:)=CST%XG/ZTHVM_F(D%NIB:D%NIE,:)
!$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)

!  Definition de l'alimentation au sens de la fermeture de Hourdin et al

ZALIM_STAR(:,:)   = 0.
ZALIM_STAR_TOT(:) = 0.    ! <== Normalization of ZALIM_STAR
IALIM(:)          = D%NKB   ! <== Top level of the alimentation layer

DO JK=D%NKB,D%NKE-D%NKL,D%NKL   !  Vertical loop
  !$mnh_expand_where(JI=D%NIB:D%NIE)
  ZZDZ(D%NIB:D%NIE,JK)   = MAX(ZEPS,PZZ(D%NIB:D%NIE,JK+D%NKL)-PZZ(D%NIB:D%NIE,JK))       ! <== Delta Z between two flux level
  ZZZ(D%NIB:D%NIE,JK)    = MAX(0.,0.5*(PZZ(D%NIB:D%NIE,JK+D%NKL)+PZZ(D%NIB:D%NIE,JK)) )  ! <== Hight of mass levels
  ZDTHETASDZ(D%NIB:D%NIE,JK) = (ZTHVM_F(D%NIB:D%NIE,JK)-ZTHVM_F(D%NIB:D%NIE,JK+D%NKL))   ! <== Delta theta_v
  
  WHERE ((ZTHVM_F(D%NIB:D%NIE,JK+D%NKL)<ZTHVM_F(D%NIB:D%NIE,JK)) .AND. (ZTHVM_F(D%NIB:D%NIE,D%NKB)>=ZTHVM_F(D%NIB:D%NIE,JK)))
     ZALIM_STAR(D%NIB:D%NIE,JK)  = SQRT(ZZZ(D%NIB:D%NIE,JK))*ZDTHETASDZ(D%NIB:D%NIE,JK)/ZZDZ(D%NIB:D%NIE,JK)
     ZALIM_STAR_TOT(D%NIB:D%NIE) = ZALIM_STAR_TOT(D%NIB:D%NIE)+ZALIM_STAR(D%NIB:D%NIE,JK)*ZZDZ(D%NIB:D%NIE,JK)
     IALIM(D%NIB:D%NIE)          = JK
  ENDWHERE
  !$mnh_end_expand_where(JI=D%NIB:D%NIE)
ENDDO

! Normalization of ZALIM_STAR
DO JK=D%NKB,D%NKE-D%NKL,D%NKL   !  Vertical loop   
   !$mnh_expand_where(JI=D%NIB:D%NIE)
   WHERE (ZALIM_STAR_TOT(D%NIB:D%NIE) > ZEPS)
      ZALIM_STAR(D%NIB:D%NIE,JK)  = ZALIM_STAR(D%NIB:D%NIE,JK)/ZALIM_STAR_TOT(D%NIB:D%NIE)
   ENDWHERE
   !$mnh_end_expand_where(JI=D%NIB:D%NIE)
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
  !$mnh_expand_where(JI=D%NIB:D%NIE)
  ! IF the updraft top is reached for all column, stop the loop on levels

  !ITEST=COUNT(GTEST(D%NIB:D%NIE))
  !IF (ITEST==0) CYCLE

  !       Computation of entrainment and detrainment with KF90
  !       parameterization in clouds and LR01 in subcloud layer

  ! to find the LCL (check if JK is LCL or not)
  WHERE ((PRC_UP(D%NIB:D%NIE,JK)+PRI_UP(D%NIB:D%NIE,JK)>0.).AND.(.NOT.(GTESTLCL(D%NIB:D%NIE))))
    KKLCL(D%NIB:D%NIE) = JK           
    GTESTLCL(D%NIB:D%NIE)=.TRUE.
  ENDWHERE

  ! COMPUTE PENTR and PDETR at mass level JK

    
  !  Buoyancy is computed on "flux" levels where updraft variables are known

  ! Compute theta_v of updraft at flux level JK    
    
  ZRC_UP(D%NIB:D%NIE)  = PRC_UP(D%NIB:D%NIE,JK)
  ZRI_UP(D%NIB:D%NIE)  = PRI_UP(D%NIB:D%NIE,JK) ! guess 
  ZRV_UP(D%NIB:D%NIE)  = PRV_UP(D%NIB:D%NIE,JK)
  ZBUO(D%NIB:D%NIE,JK) = ZG_O_THVREF(D%NIB:D%NIE,JK)*(PTHV_UP(D%NIB:D%NIE,JK) - ZTHVM_F(D%NIB:D%NIE,JK))
  PBUO_INTEG(D%NIB:D%NIE,JK) = ZBUO(D%NIB:D%NIE,JK)*(PZZ(D%NIB:D%NIE,JK+D%NKL)-PZZ(D%NIB:D%NIE,JK))
  
  ZDZ(D%NIB:D%NIE)   = MAX(ZEPS,PZZ(D%NIB:D%NIE,JK+D%NKL)-PZZ(D%NIB:D%NIE,JK))
  ZTEST(D%NIB:D%NIE) = PARAMMF%XA1*ZBUO(D%NIB:D%NIE,JK) -  PARAMMF%XB*ZW_UP2(D%NIB:D%NIE,JK)

  ZCOE(D%NIB:D%NIE)      = ZDZ(D%NIB:D%NIE)
  WHERE (ZTEST(D%NIB:D%NIE)>0.)
    ZCOE(D%NIB:D%NIE)    = ZDZ(D%NIB:D%NIE)/(1.+ PARAMMF%XBETA1)
  ENDWHERE   

  !  Calcul de la vitesse

  ZWCOE(D%NIB:D%NIE)         = (1.-PARAMMF%XB*ZCOE(D%NIB:D%NIE))/(1.+PARAMMF%XB*ZCOE(D%NIB:D%NIE))
  ZBUCOE(D%NIB:D%NIE)        =  2.*ZCOE(D%NIB:D%NIE)/(1.+PARAMMF%XB*ZCOE(D%NIB:D%NIE))

  ZW_UP2(D%NIB:D%NIE,JK+D%NKL) = MAX(ZEPS,ZW_UP2(D%NIB:D%NIE,JK)*ZWCOE(D%NIB:D%NIE) + &
                                         &PARAMMF%XA1*ZBUO(D%NIB:D%NIE,JK)*ZBUCOE(D%NIB:D%NIE))
  ZW_MAX(D%NIB:D%NIE) = MAX(ZW_MAX(D%NIB:D%NIE), SQRT(ZW_UP2(D%NIB:D%NIE,JK+D%NKL)))
  ZWUP_MEAN(D%NIB:D%NIE)     = MAX(ZEPS,0.5*(ZW_UP2(D%NIB:D%NIE,JK+D%NKL)+ZW_UP2(D%NIB:D%NIE,JK)))

  !  Entrainement et detrainement

  PENTR(D%NIB:D%NIE,JK)  = MAX(0.,(PARAMMF%XBETA1/(1.+PARAMMF%XBETA1))* &
                                 &(PARAMMF%XA1*ZBUO(D%NIB:D%NIE,JK)/ZWUP_MEAN(D%NIB:D%NIE)-PARAMMF%XB))
  
  ZDETR_BUO(D%NIB:D%NIE) = MAX(0., -(PARAMMF%XBETA1/(1.+PARAMMF%XBETA1))*PARAMMF%XA1*ZBUO(D%NIB:D%NIE,JK)/ &
                                   &ZWUP_MEAN(D%NIB:D%NIE))
  ZDETR_RT(D%NIB:D%NIE)  = PARAMMF%XC*SQRT(MAX(0.,(PRT_UP(D%NIB:D%NIE,JK) - ZRTM_F(D%NIB:D%NIE,JK))) / &
                                          &MAX(ZEPS,ZRTM_F(D%NIB:D%NIE,JK)) / ZWUP_MEAN(D%NIB:D%NIE))
  PDETR(D%NIB:D%NIE,JK)  = ZDETR_RT(D%NIB:D%NIE)+ZDETR_BUO(D%NIB:D%NIE)

   
  ! If the updraft did not stop, compute cons updraft characteritics at jk+1
  WHERE(GTEST(D%NIB:D%NIE))     
    ZZTOP(D%NIB:D%NIE) = MAX(ZZTOP(D%NIB:D%NIE),PZZ(D%NIB:D%NIE,JK+D%NKL))
    ZMIX2(D%NIB:D%NIE) = (PZZ(D%NIB:D%NIE,JK+D%NKL)-PZZ(D%NIB:D%NIE,JK))*PENTR(D%NIB:D%NIE,JK) !&
    ZMIX3(D%NIB:D%NIE) = (PZZ(D%NIB:D%NIE,JK+D%NKL)-PZZ(D%NIB:D%NIE,JK))*PDETR(D%NIB:D%NIE,JK) !&                   
           
    ZQTM(D%NIB:D%NIE) = PRTM(D%NIB:D%NIE,JK)/(1.+PRTM(D%NIB:D%NIE,JK))            
    ZTHSM(D%NIB:D%NIE,JK) = PTHLM(D%NIB:D%NIE,JK)*(1.+PARAMMF%XLAMBDA_MF*ZQTM(D%NIB:D%NIE))
    ZTHS_UP(D%NIB:D%NIE,JK+D%NKL)=(ZTHS_UP(D%NIB:D%NIE,JK)*(1.-0.5*ZMIX2(D%NIB:D%NIE)) + ZTHSM(D%NIB:D%NIE,JK)*ZMIX2(D%NIB:D%NIE))&
                          /(1.+0.5*ZMIX2(D%NIB:D%NIE))   
    PRT_UP(D%NIB:D%NIE,JK+D%NKL)=(PRT_UP(D%NIB:D%NIE,JK)*(1.-0.5*ZMIX2(D%NIB:D%NIE)) + PRTM(D%NIB:D%NIE,JK)*ZMIX2(D%NIB:D%NIE))  &
                          /(1.+0.5*ZMIX2(D%NIB:D%NIE))
    ZQT_UP(D%NIB:D%NIE) = PRT_UP(D%NIB:D%NIE,JK+D%NKL)/(1.+PRT_UP(D%NIB:D%NIE,JK+D%NKL))
    PTHL_UP(D%NIB:D%NIE,JK+D%NKL)=ZTHS_UP(D%NIB:D%NIE,JK+D%NKL)/(1.+PARAMMF%XLAMBDA_MF*ZQT_UP(D%NIB:D%NIE))                      
  ENDWHERE
  

  IF(OMIXUV) THEN
    IF(JK/=D%NKB) THEN
      WHERE(GTEST(D%NIB:D%NIE))
        PU_UP(D%NIB:D%NIE,JK+D%NKL) = (PU_UP(D%NIB:D%NIE,JK)*(1-0.5*ZMIX2(D%NIB:D%NIE)) + PUM(D%NIB:D%NIE,JK)*ZMIX2(D%NIB:D%NIE)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(D%NIB:D%NIE,JK+D%NKL)-PZZ(D%NIB:D%NIE,JK))*&
                          ((PUM(D%NIB:D%NIE,JK+D%NKL)-PUM(D%NIB:D%NIE,JK))/PDZZ(D%NIB:D%NIE,JK+D%NKL)+&
                           (PUM(D%NIB:D%NIE,JK)-PUM(D%NIB:D%NIE,JK-D%NKL))/PDZZ(D%NIB:D%NIE,JK))        )   &
                          /(1+0.5*ZMIX2(D%NIB:D%NIE))
        PV_UP(D%NIB:D%NIE,JK+D%NKL) = (PV_UP(D%NIB:D%NIE,JK)*(1-0.5*ZMIX2(D%NIB:D%NIE)) + PVM(D%NIB:D%NIE,JK)*ZMIX2(D%NIB:D%NIE)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(D%NIB:D%NIE,JK+D%NKL)-PZZ(D%NIB:D%NIE,JK))*&
                          ((PVM(D%NIB:D%NIE,JK+D%NKL)-PVM(D%NIB:D%NIE,JK))/PDZZ(D%NIB:D%NIE,JK+D%NKL)+&
                           (PVM(D%NIB:D%NIE,JK)-PVM(D%NIB:D%NIE,JK-D%NKL))/PDZZ(D%NIB:D%NIE,JK))    )   &
                          /(1+0.5*ZMIX2(D%NIB:D%NIE))
      ENDWHERE
    ELSE
      WHERE(GTEST(D%NIB:D%NIE))
        PU_UP(D%NIB:D%NIE,JK+D%NKL) = (PU_UP(D%NIB:D%NIE,JK)*(1-0.5*ZMIX2(D%NIB:D%NIE)) + PUM(D%NIB:D%NIE,JK)*ZMIX2(D%NIB:D%NIE)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(D%NIB:D%NIE,JK+D%NKL)-PZZ(D%NIB:D%NIE,JK))*&
                          ((PUM(D%NIB:D%NIE,JK+D%NKL)-PUM(D%NIB:D%NIE,JK))/PDZZ(D%NIB:D%NIE,JK+D%NKL))        )   &
                          /(1+0.5*ZMIX2(D%NIB:D%NIE))
        PV_UP(D%NIB:D%NIE,JK+D%NKL) = (PV_UP(D%NIB:D%NIE,JK)*(1-0.5*ZMIX2(D%NIB:D%NIE)) + PVM(D%NIB:D%NIE,JK)*ZMIX2(D%NIB:D%NIE)+ &
                          0.5*PARAMMF%XPRES_UV*(PZZ(D%NIB:D%NIE,JK+D%NKL)-PZZ(D%NIB:D%NIE,JK))*&
                          ((PVM(D%NIB:D%NIE,JK+D%NKL)-PVM(D%NIB:D%NIE,JK))/PDZZ(D%NIB:D%NIE,JK+D%NKL))    )   &
                          /(1+0.5*ZMIX2(D%NIB:D%NIE))
      ENDWHERE

    ENDIF
  ENDIF
  !DO JSV=1,ISV 
  !   IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  !    WHERE(GTEST(D%NIB:D%NIE)) 
  !         PSV_UP(D%NIB:D%NIE,JK+KKL,JSV) = (PSV_UP(D%NIB:D%NIE,JK,JSV)*(1-0.5*ZMIX2(D%NIB:D%NIE)) + &
  !                      PSVM(D%NIB:D%NIE,JK,JSV)*ZMIX2(D%NIB:D%NIE))  /(1+0.5*ZMIX2(D%NIB:D%NIE))
  !    ENDWHERE                        
  !ENDDO  
  
  
  ! Compute non cons. var. at level JK+KKL
  ZRC_UP(D%NIB:D%NIE)=PRC_UP(D%NIB:D%NIE,JK) ! guess = level just below
  ZRI_UP(D%NIB:D%NIE)=PRI_UP(D%NIB:D%NIE,JK) ! guess = level just below
  ZRV_UP(D%NIB:D%NIE)=PRV_UP(D%NIB:D%NIE,JK)
  !$mnh_end_expand_where(JI=D%NIB:D%NIE)
  CALL TH_R_FROM_THL_RT(CST,NEB, D%NIT, HFRAC_ICE,PFRAC_ICE_UP(:,JK+D%NKL),ZPRES_F(:,JK+D%NKL), &
          PTHL_UP(:,JK+D%NKL),PRT_UP(:,JK+D%NKL),ZTH_UP(:,JK+D%NKL),              &
          ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.,&
          PBUF=ZBUF, KB=D%NIB, KE=D%NIE)
  !$mnh_expand_where(JI=D%NIB:D%NIE)
  WHERE(GTEST(D%NIB:D%NIE))
    ZT_UP(D%NIB:D%NIE) = ZTH_UP(D%NIB:D%NIE,JK+D%NKL)*PEXNM(D%NIB:D%NIE,JK+D%NKL)
    ZCP(D%NIB:D%NIE) = CST%XCPD + CST%XCL * ZRC_UP(D%NIB:D%NIE)
    ZLVOCPEXN(D%NIB:D%NIE)=(CST%XLVTT + (CST%XCPV-CST%XCL) *  (ZT_UP(D%NIB:D%NIE)-CST%XTT) ) / &
                          &ZCP(D%NIB:D%NIE) / PEXNM(D%NIB:D%NIE,JK+D%NKL)
    PRC_UP(D%NIB:D%NIE,JK+D%NKL)=MIN(0.5E-3,ZRC_UP(D%NIB:D%NIE))  ! On ne peut depasser 0.5 g/kg (autoconversion donc elimination !)
    PTHL_UP(D%NIB:D%NIE,JK+D%NKL) = PTHL_UP(D%NIB:D%NIE,JK+D%NKL)+ &
                                  & ZLVOCPEXN(D%NIB:D%NIE)*(ZRC_UP(D%NIB:D%NIE)-PRC_UP(D%NIB:D%NIE,JK+D%NKL))
    PRV_UP(D%NIB:D%NIE,JK+D%NKL)=ZRV_UP(D%NIB:D%NIE)
    PRI_UP(D%NIB:D%NIE,JK+D%NKL)=ZRI_UP(D%NIB:D%NIE)
    PRT_UP(D%NIB:D%NIE,JK+D%NKL)  = PRC_UP(D%NIB:D%NIE,JK+D%NKL) + PRV_UP(D%NIB:D%NIE,JK+D%NKL)
    PRSAT_UP(D%NIB:D%NIE,JK+D%NKL) = ZRSATW(D%NIB:D%NIE)*(1-PFRAC_ICE_UP(D%NIB:D%NIE,JK+D%NKL)) + &
                                   & ZRSATI(D%NIB:D%NIE)*PFRAC_ICE_UP(D%NIB:D%NIE,JK+D%NKL)
  ENDWHERE
  

  ! Compute the updraft theta_v, buoyancy and w**2 for level JK+1   
  WHERE(GTEST(D%NIB:D%NIE))
    !PTHV_UP(D%NIB:D%NIE,JK+KKL) = ZTH_UP(D%NIB:D%NIE,JK+KKL)*((1+ZRVORD*PRV_UP(D%NIB:D%NIE,JK+KKL))/(1+PRT_UP(D%NIB:D%NIE,JK+KKL)))
    PTHV_UP(D%NIB:D%NIE,JK+D%NKL) = ZTH_UP(D%NIB:D%NIE,JK+D%NKL)* &
                                  & (1.+0.608*PRV_UP(D%NIB:D%NIE,JK+D%NKL) - PRC_UP(D%NIB:D%NIE,JK+D%NKL))
  ENDWHERE


  ! Test if the updraft has reach the ETL
  GTESTETL(D%NIB:D%NIE)=.FALSE.
  WHERE (GTEST(D%NIB:D%NIE).AND.(PBUO_INTEG(D%NIB:D%NIE,JK)<=0.))
    KKETL(D%NIB:D%NIE) = JK+D%NKL
    GTESTETL(D%NIB:D%NIE)=.TRUE.
  ENDWHERE

  ! Test is we have reached the top of the updraft
  WHERE (GTEST(D%NIB:D%NIE).AND.((ZW_UP2(D%NIB:D%NIE,JK+D%NKL)<=ZEPS)))
    ZW_UP2(D%NIB:D%NIE,JK+D%NKL)=ZEPS
    GTEST(D%NIB:D%NIE)=.FALSE.
    PTHL_UP(D%NIB:D%NIE,JK+D%NKL)=ZTHLM_F(D%NIB:D%NIE,JK+D%NKL)
    PRT_UP(D%NIB:D%NIE,JK+D%NKL)=ZRTM_F(D%NIB:D%NIE,JK+D%NKL)
    PRC_UP(D%NIB:D%NIE,JK+D%NKL)=0.
    PRI_UP(D%NIB:D%NIE,JK+D%NKL)=0.
    PRV_UP(D%NIB:D%NIE,JK+D%NKL)=0.
    PTHV_UP(D%NIB:D%NIE,JK+D%NKL)=ZTHVM_F(D%NIB:D%NIE,JK+D%NKL)
    PFRAC_UP(D%NIB:D%NIE,JK+D%NKL)=0.
    KKCTL(D%NIB:D%NIE)=JK+D%NKL
  ENDWHERE
  !$mnh_end_expand_where(JI=D%NIB:D%NIE)
ENDDO 

! Closure assumption for mass flux at KKB+1 level (Mass flux is supposed to be 0 at KKB level !)                                                 
!   Hourdin et al 2002 formulation


!$mnh_expand_array(JI=D%NIB:D%NIE)
ZZTOP(D%NIB:D%NIE) = MAX(ZZTOP(D%NIB:D%NIE),ZEPS)
!$mnh_end_expand_array(JI=D%NIB:D%NIE)

DO JK=D%NKB+D%NKL,D%NKE-D%NKL,D%NKL !  Vertical loop
  !$mnh_expand_where(JI=D%NIB:D%NIE)
   WHERE(JK<=IALIM(D%NIB:D%NIE))
     ZALIM_STAR_TOT(D%NIB:D%NIE) = ZALIM_STAR_TOT(D%NIB:D%NIE) + ZALIM_STAR(D%NIB:D%NIE,JK)**2* &
                                                               & ZZDZ(D%NIB:D%NIE,JK)/PRHODREF(D%NIB:D%NIE,JK)
   ENDWHERE
  !$mnh_end_expand_where(JI=D%NIB:D%NIE)
ENDDO   

!$mnh_expand_where(JI=D%NIB:D%NIE)
WHERE (ZALIM_STAR_TOT(D%NIB:D%NIE)*ZZTOP(D%NIB:D%NIE) > ZEPS)
 ZPHI(D%NIB:D%NIE) =  ZW_MAX(D%NIB:D%NIE)/(PARAMMF%XR*ZZTOP(D%NIB:D%NIE)*ZALIM_STAR_TOT(D%NIB:D%NIE))
ENDWHERE

GTEST(D%NIB:D%NIE) = .TRUE.
PEMF(D%NIB:D%NIE,D%NKB+D%NKL) = ZPHI(D%NIB:D%NIE)*ZZDZ(D%NIB:D%NIE,D%NKB)*ZALIM_STAR(D%NIB:D%NIE,D%NKB)
! Updraft fraction must be smaller than XFRAC_UP_MAX
PFRAC_UP(D%NIB:D%NIE,D%NKB+D%NKL)=PEMF(D%NIB:D%NIE,D%NKB+D%NKL)/ &
                                 &(SQRT(ZW_UP2(D%NIB:D%NIE,D%NKB+D%NKL))*ZRHO_F(D%NIB:D%NIE,D%NKB+D%NKL))
PFRAC_UP(D%NIB:D%NIE,D%NKB+D%NKL)=MIN(PARAMMF%XFRAC_UP_MAX,PFRAC_UP(D%NIB:D%NIE,D%NKB+D%NKL))
PEMF(D%NIB:D%NIE,D%NKB+D%NKL) = ZRHO_F(D%NIB:D%NIE,D%NKB+D%NKL)*PFRAC_UP(D%NIB:D%NIE,D%NKB+D%NKL)* &
                              & SQRT(ZW_UP2(D%NIB:D%NIE,D%NKB+D%NKL))
!$mnh_end_expand_where(JI=D%NIB:D%NIE)

DO JK=D%NKB+D%NKL,D%NKE-D%NKL,D%NKL !  Vertical loop
  !$mnh_expand_where(JI=D%NIB:D%NIE)
     
  GTEST(D%NIB:D%NIE) = (ZW_UP2(D%NIB:D%NIE,JK) > ZEPS)  

  WHERE (GTEST(D%NIB:D%NIE))
    WHERE(JK<IALIM(D%NIB:D%NIE))
      PEMF(D%NIB:D%NIE,JK+D%NKL) = MAX(0.,PEMF(D%NIB:D%NIE,JK) + ZPHI(D%NIB:D%NIE)*ZZDZ(D%NIB:D%NIE,JK)* &
                                                                & (PENTR(D%NIB:D%NIE,JK) - PDETR(D%NIB:D%NIE,JK)))
    ELSEWHERE
      ZMIX1(D%NIB:D%NIE)=ZZDZ(D%NIB:D%NIE,JK)*(PENTR(D%NIB:D%NIE,JK)-PDETR(D%NIB:D%NIE,JK))
      PEMF(D%NIB:D%NIE,JK+D%NKL)=PEMF(D%NIB:D%NIE,JK)*EXP(ZMIX1(D%NIB:D%NIE))
    ENDWHERE

! Updraft fraction must be smaller than XFRAC_UP_MAX
    PFRAC_UP(D%NIB:D%NIE,JK+D%NKL)=PEMF(D%NIB:D%NIE,JK+D%NKL)/(SQRT(ZW_UP2(D%NIB:D%NIE,JK+D%NKL))*ZRHO_F(D%NIB:D%NIE,JK+D%NKL))
    PFRAC_UP(D%NIB:D%NIE,JK+D%NKL)=MIN(PARAMMF%XFRAC_UP_MAX,PFRAC_UP(D%NIB:D%NIE,JK+D%NKL))
    PEMF(D%NIB:D%NIE,JK+D%NKL) = ZRHO_F(D%NIB:D%NIE,JK+D%NKL)*PFRAC_UP(D%NIB:D%NIE,JK+D%NKL)*SQRT(ZW_UP2(D%NIB:D%NIE,JK+D%NKL))
  ENDWHERE
  !$mnh_end_expand_where(JI=D%NIB:D%NIE)
ENDDO

!$mnh_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
PW_UP(D%NIB:D%NIE,:)=SQRT(ZW_UP2(D%NIB:D%NIE,:))
!$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
!$mnh_expand_array(JI=D%NIB:D%NIE)
PEMF(D%NIB:D%NIE,D%NKB) =0.
!$mnh_end_expand_array(JI=D%NIB:D%NIE)

! Limits the shallow convection scheme when cloud heigth is higher than 3000m.
! To do this, mass flux is multiplied by a coefficient decreasing linearly
! from 1 (for clouds of 3000m of depth) to 0 (for clouds of 4000m of depth).
! This way, all MF fluxes are diminished by this amount.
! Diagnosed cloud fraction is also multiplied by the same coefficient.
!
DO JI=D%NIB,D%NIE 
  PDEPTH(JI) = MAX(0., PZZ(JI,KKCTL(JI)) -  PZZ(JI,KKLCL(JI)) )
END DO

!$mnh_expand_array(JI=D%NIB:D%NIE)
GWORK1(D%NIB:D%NIE)= (GTESTLCL(D%NIB:D%NIE) .AND. (PDEPTH(D%NIB:D%NIE) > ZDEPTH_MAX1) )
!$mnh_end_expand_array(JI=D%NIB:D%NIE)
DO JK=1,D%NKT
  !$mnh_expand_array(JI=D%NIB:D%NIE)
  GWORK2(D%NIB:D%NIE,JK) = GWORK1(D%NIB:D%NIE)
  ZCOEF(D%NIB:D%NIE,JK) = (1.-(PDEPTH(D%NIB:D%NIE)-ZDEPTH_MAX1)/(ZDEPTH_MAX2-ZDEPTH_MAX1))
  ZCOEF(D%NIB:D%NIE,JK)=MIN(MAX(ZCOEF(D%NIB:D%NIE,JK),0.),1.)
  !$mnh_end_expand_array(JI=D%NIB:D%NIE)
ENDDO
!$mnh_expand_where(JI=D%NIB:D%NIE,JK=1:D%NKT)
WHERE (GWORK2(D%NIB:D%NIE,:)) 
  PEMF(D%NIB:D%NIE,:)     = PEMF(D%NIB:D%NIE,:)     * ZCOEF(D%NIB:D%NIE,:)
  PFRAC_UP(D%NIB:D%NIE,:) = PFRAC_UP(D%NIB:D%NIE,:) * ZCOEF(D%NIB:D%NIE,:)
ENDWHERE
!$mnh_end_expand_where(JI=D%NIB:D%NIE,JK=1:D%NKT)

IF (LHOOK) CALL DR_HOOK('COMPUTE_UPDRAF_RAHA',1,ZHOOK_HANDLE)
!
CONTAINS
INCLUDE "th_r_from_thl_rt.func.h"
INCLUDE "compute_frac_ice.func.h"
!
END SUBROUTINE COMPUTE_UPDRAFT_RAHA
END MODULE MODE_COMPUTE_UPDRAFT_RAHA
