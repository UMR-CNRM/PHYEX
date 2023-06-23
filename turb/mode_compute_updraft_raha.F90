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
      SUBROUTINE COMPUTE_UPDRAFT_RAHA(D, CST, NEBN, PARAMMF,         &
                                 KSV, OENTR_DETR,                    &
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
USE MODD_NEB_n,           ONLY: NEB_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
!
USE MODI_SHUMAN_MF, ONLY: MZM_MF
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

!*                    1.1  Declaration of Arguments
!
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(NEB_t),            INTENT(IN)   :: NEBN
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF
INTEGER,                INTENT(IN)   :: KSV
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
REAL, DIMENSION(D%NIJT,D%NKT) :: ZTHM_F,ZRVM_F                 ! Theta,rv of updraft environnement
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRTM_F, ZTHLM_F, ZTKEM_F      ! rt, thetal,TKE,pressure,
REAL, DIMENSION(D%NIJT,D%NKT) :: ZUM_F,ZVM_F,ZRHO_F            ! density,momentum
REAL, DIMENSION(D%NIJT,D%NKT) :: ZPRES_F,ZTHVM_F,ZTHVM         ! interpolated at the flux point
REAL, DIMENSION(D%NIJT,D%NKT) :: ZG_O_THVREF                   ! g*ThetaV ref
REAL, DIMENSION(D%NIJT,D%NKT) :: ZW_UP2                        ! w**2  of the updraft

REAL, DIMENSION(D%NIJT,D%NKT) :: ZTH_UP                        ! updraft THETA 
REAL, DIMENSION(D%NIJT)              :: ZT_UP                         ! updraft T
REAL, DIMENSION(D%NIJT)              :: ZLVOCPEXN                     ! updraft L
REAL, DIMENSION(D%NIJT)              :: ZCP                           ! updraft cp
REAL, DIMENSION(D%NIJT,D%NKT) :: ZBUO                          ! Buoyancy 
REAL, DIMENSION(D%NIJT,D%NKT) :: ZTHS_UP,ZTHSM

REAL, DIMENSION(D%NIJT,D%NKT) ::  ZCOEF  ! diminution coefficient for too high clouds 
                        
REAL  :: ZRDORV       ! RD/RV
REAL  :: ZRVORD       ! RV/RD


REAL, DIMENSION(D%NIJT) :: ZMIX1,ZMIX2,ZMIX3


INTEGER  :: JK,JIJ          ! loop counters
INTEGER :: IIJB,IIJE ! physical horizontal domain indices
INTEGEr :: IKT,IKB,IKE,IKL

LOGICAL, DIMENSION(D%NIJT) ::  GTEST,GTESTLCL,GTESTETL
                               ! Test if the ascent continue, if LCL or ETL is reached
LOGICAL, DIMENSION(D%NIJT)              :: GWORK1
LOGICAL, DIMENSION(D%NIJT,D%NKT) :: GWORK2



REAL, DIMENSION(D%NIJT) :: ZRC_UP, ZRI_UP, ZRV_UP, ZRSATW, ZRSATI

REAL,    DIMENSION(D%NIJT) :: ZPHI,ZALIM_STAR_TOT
REAL,    DIMENSION(D%NIJT,D%NKT) :: ZDTHETASDZ,ZALIM_STAR,ZZDZ,ZZZ
INTEGER, DIMENSION(D%NIJT) :: IALIM

REAL, DIMENSION(D%NIJT)              ::  ZTEST,ZDZ,ZWUP_MEAN    ! 
REAL, DIMENSION(D%NIJT)              ::  ZCOE,ZWCOE,ZBUCOE
REAL, DIMENSION(D%NIJT)              ::  ZDETR_BUO, ZDETR_RT
REAL, DIMENSION(D%NIJT)              ::  ZW_MAX               ! w**2  max of the updraft
REAL, DIMENSION(D%NIJT)              ::  ZZTOP                ! Top of the updraft
REAL, DIMENSION(D%NIJT)              ::  ZQTM,ZQT_UP

REAL  :: ZDEPTH_MAX1, ZDEPTH_MAX2 ! control auto-extinction process

REAL  :: ZTMAX,ZRMAX, ZEPS  ! control value
REAL, DIMENSION(D%NIJT,16)           ::  ZBUF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
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
DO JK=1,IKT 
  DO JIJ=IIJB,IIJE 
    PRSAT_UP(JIJ,JK)=PRVM(JIJ,JK) ! should be initialised correctly but is (normaly) not used
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
!   ZSVM_F(IIJB:IIJE,KKB:IKU,JSV) = 0.5*(PSVM(IIJB:IIJE,KKB:IKU,JSV)+PSVM(IIJB:IIJE,1:IKU-1,JSV))
!   ZSVM_F(IIJB:IIJE,1,JSV)       = ZSVM_F(IIJB:IIJE,KKB,JSV) 
!END DO  

!          Initialisation of updraft characteristics 
DO JK=1,IKT 
  DO JIJ=IIJB,IIJE 
    PTHL_UP(JIJ,JK)=ZTHLM_F(JIJ,JK)
    PRT_UP(JIJ,JK)=ZRTM_F(JIJ,JK)
    PU_UP(JIJ,JK)=ZUM_F(JIJ,JK)
    PV_UP(JIJ,JK)=ZVM_F(JIJ,JK)
  ENDDO
ENDDO
PSV_UP(:,:,:)=0.
!IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) then
!    PSV_UP(IIJB:IIJE,:,:)=ZSVM_F(IIJB:IIJE,:,:)
!ENDIF

! Computation or initialisation of updraft characteristics at the KKB level
! thetal_up,rt_up,thetaV_up, w�,Buoyancy term and mass flux (PEMF)

DO JIJ=IIJB,IIJE 
  PTHL_UP(JIJ,IKB)= ZTHLM_F(JIJ,IKB)+ &
  & MAX(0.,MIN(ZTMAX,(PSFTH(JIJ)/SQRT(ZTKEM_F(JIJ,IKB)))*PARAMMF%XALP_PERT))
  PRT_UP(JIJ,IKB) = ZRTM_F(JIJ,IKB)+ &
  & MAX(0.,MIN(ZRMAX,(PSFRV(JIJ)/SQRT(ZTKEM_F(JIJ,IKB)))*PARAMMF%XALP_PERT)) 
  
  ZQT_UP(JIJ) = PRT_UP(JIJ,IKB)/(1.+PRT_UP(JIJ,IKB))
  ZTHS_UP(JIJ,IKB)=PTHL_UP(JIJ,IKB)*(1.+PARAMMF%XLAMBDA_MF*ZQT_UP(JIJ))
ENDDO

CALL MZM_MF(D, PTHM (:,:), ZTHM_F(:,:))
CALL MZM_MF(D, PPABSM(:,:), ZPRES_F(:,:))
CALL MZM_MF(D, PRHODREF(:,:), ZRHO_F(:,:))
CALL MZM_MF(D, PRVM(:,:), ZRVM_F(:,:))

DO JK=1,IKT 
  DO JIJ=IIJB,IIJE 
! thetav at mass and flux levels 
    ZTHVM_F(JIJ,JK)=ZTHM_F(JIJ,JK)*((1.+ZRVORD*ZRVM_F(JIJ,JK))/&
    &(1.+ZRTM_F(JIJ,JK)))
    ZTHVM(JIJ,JK)=PTHM(JIJ,JK)*((1.+ZRVORD*PRVM(JIJ,JK))/(1.+PRTM(JIJ,JK)))
    
    PTHV_UP(JIJ,JK)= ZTHVM_F(JIJ,JK)
    PRV_UP(JIJ,JK) = ZRVM_F(JIJ,JK)
  ENDDO
ENDDO

ZW_UP2(:,:)=ZEPS
DO JIJ=IIJB,IIJE 
  ZW_UP2(JIJ,IKB) = MAX(0.0001,(1./6.)*ZTKEM_F(JIJ,IKB))
  GTEST(JIJ) = (ZW_UP2(JIJ,IKB) > ZEPS)  
ENDDO

! Computation of non conservative variable for the KKB level of the updraft
! (all or nothing ajustement)
DO JIJ=IIJB,IIJE 
  PRC_UP(JIJ,IKB)=0.
  PRI_UP(JIJ,IKB)=0.
ENDDO

CALL TH_R_FROM_THL_RT(CST, NEBN, D%NIJT, NEBN%CFRAC_ICE_SHALLOW_MF,PFRAC_ICE_UP(:,IKB),ZPRES_F(:,IKB), &
             PTHL_UP(:,IKB),PRT_UP(:,IKB),ZTH_UP(:,IKB), &
             PRV_UP(:,IKB),PRC_UP(:,IKB),PRI_UP(:,IKB),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.,&
             PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)

DO JIJ=IIJB,IIJE 
! compute updraft thevav and buoyancy term at KKB level             
  PTHV_UP(JIJ,IKB) = ZTH_UP(JIJ,IKB)*((1+ZRVORD*PRV_UP(JIJ,IKB))/(1+PRT_UP(JIJ,IKB))) 
! compute mean rsat in updraft
  PRSAT_UP(JIJ,IKB) = ZRSATW(JIJ)*(1-PFRAC_ICE_UP(JIJ,IKB)) + &
  & ZRSATI(JIJ)*PFRAC_ICE_UP(JIJ,IKB)
ENDDO

!Tout est commente pour tester dans un premier temps la s�paration en deux de la 
!  boucle verticale, une pour w et une pour PEMF                                                            
DO JK=1,IKT 
  DO JIJ=IIJB,IIJE 
    ZG_O_THVREF(JIJ,JK)=CST%XG/ZTHVM_F(JIJ,JK)
  ENDDO
ENDDO

!  Definition de l'alimentation au sens de la fermeture de Hourdin et al

ZALIM_STAR(:,:)   = 0.
ZALIM_STAR_TOT(:) = 0.    ! <== Normalization of ZALIM_STAR
IALIM(:)          = IKB   ! <== Top level of the alimentation layer

DO JK=IKB,IKE-IKL,IKL   !  Vertical loop
  DO JIJ=IIJB,IIJE 
    ZZDZ(JIJ,JK)   = MAX(ZEPS,PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))       ! <== Delta Z between two flux level
    ZZZ(JIJ,JK)    = MAX(0.,0.5*(PZZ(JIJ,JK+IKL)+PZZ(JIJ,JK)) )  ! <== Hight of mass levels
    ZDTHETASDZ(JIJ,JK) = (ZTHVM_F(JIJ,JK)-ZTHVM_F(JIJ,JK+IKL))   ! <== Delta theta_v
    
    IF ((ZTHVM_F(JIJ,JK+IKL)<ZTHVM_F(JIJ,JK)) .AND. &
        &(ZTHVM_F(JIJ,IKB)>=ZTHVM_F(JIJ,JK)))THEN
    ZALIM_STAR(JIJ,JK)  = SQRT(ZZZ(JIJ,JK))*ZDTHETASDZ(JIJ,JK)/ZZDZ(JIJ,JK)
    ZALIM_STAR_TOT(JIJ) = ZALIM_STAR_TOT(JIJ)+ZALIM_STAR(JIJ,JK)*ZZDZ(JIJ,JK)
    IALIM(JIJ)          = JK
  ENDIF
ENDDO
ENDDO

! Normalization of ZALIM_STAR
DO JK=IKB,IKE-IKL,IKL   !  Vertical loop   
DO JIJ=IIJB,IIJE 
  IF (ZALIM_STAR_TOT(JIJ) > ZEPS)THEN
    ZALIM_STAR(JIJ,JK)  = ZALIM_STAR(JIJ,JK)/ZALIM_STAR_TOT(JIJ)
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


DO JK=IKB,IKE-IKL,IKL
DO JIJ=IIJB,IIJE 
  ! IF the updraft top is reached for all column, stop the loop on levels
  
  !ITEST=COUNT(GTEST(JIJ))
  !IF (ITEST==0) CYCLE
  
  !       Computation of entrainment and detrainment with KF90
  !       parameterization in clouds and LR01 in subcloud layer
  
  ! to find the LCL (check if JK is LCL or not)
  IF ((PRC_UP(JIJ,JK)+PRI_UP(JIJ,JK)>0.).AND.(.NOT.(GTESTLCL(JIJ))))THEN
    KKLCL(JIJ) = JK           
    GTESTLCL(JIJ)=.TRUE.
  ENDIF
  
  ! COMPUTE PENTR and PDETR at mass level JK
  
  
  !  Buoyancy is computed on "flux" levels where updraft variables are known
  
  ! Compute theta_v of updraft at flux level JK    
  
  ZRC_UP(JIJ)  = PRC_UP(JIJ,JK)
  ZRI_UP(JIJ)  = PRI_UP(JIJ,JK) ! guess 
  ZRV_UP(JIJ)  = PRV_UP(JIJ,JK)
  ZBUO(JIJ,JK) = ZG_O_THVREF(JIJ,JK)*(PTHV_UP(JIJ,JK) - ZTHVM_F(JIJ,JK))
  PBUO_INTEG(JIJ,JK) = ZBUO(JIJ,JK)*(PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))
  
  ZDZ(JIJ)   = MAX(ZEPS,PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))
  ZTEST(JIJ) = PARAMMF%XA1*ZBUO(JIJ,JK) -  PARAMMF%XB*ZW_UP2(JIJ,JK)
  
  ZCOE(JIJ)      = ZDZ(JIJ)
  IF (ZTEST(JIJ)>0.)THEN
    ZCOE(JIJ)    = ZDZ(JIJ)/(1.+ PARAMMF%XBETA1)
  ENDIF   
  
  !  Calcul de la vitesse
  
  ZWCOE(JIJ)         = (1.-PARAMMF%XB*ZCOE(JIJ))/(1.+PARAMMF%XB*ZCOE(JIJ))
  ZBUCOE(JIJ)        =  2.*ZCOE(JIJ)/(1.+PARAMMF%XB*ZCOE(JIJ))
  
  ZW_UP2(JIJ,JK+IKL) = MAX(ZEPS,ZW_UP2(JIJ,JK)*ZWCOE(JIJ) + &
  &PARAMMF%XA1*ZBUO(JIJ,JK)*ZBUCOE(JIJ))
  ZW_MAX(JIJ) = MAX(ZW_MAX(JIJ), SQRT(ZW_UP2(JIJ,JK+IKL)))
  ZWUP_MEAN(JIJ)     = MAX(ZEPS,0.5*(ZW_UP2(JIJ,JK+IKL)+ZW_UP2(JIJ,JK)))
  
  !  Entrainement et detrainement
  
  PENTR(JIJ,JK)  = MAX(0.,(PARAMMF%XBETA1/(1.+PARAMMF%XBETA1))* &
  &(PARAMMF%XA1*ZBUO(JIJ,JK)/ZWUP_MEAN(JIJ)-PARAMMF%XB))
  
  ZDETR_BUO(JIJ) = MAX(0., -(PARAMMF%XBETA1/(1.+PARAMMF%XBETA1))*PARAMMF%XA1*ZBUO(JIJ,JK)/ &
  &ZWUP_MEAN(JIJ))
  ZDETR_RT(JIJ)  = PARAMMF%XC*SQRT(MAX(0.,(PRT_UP(JIJ,JK) - ZRTM_F(JIJ,JK))) / &
  &MAX(ZEPS,ZRTM_F(JIJ,JK)) / ZWUP_MEAN(JIJ))
  PDETR(JIJ,JK)  = ZDETR_RT(JIJ)+ZDETR_BUO(JIJ)
  
  
  ! If the updraft did not stop, compute cons updraft characteritics at jk+1
  IF(GTEST(JIJ))     THEN
    ZZTOP(JIJ) = MAX(ZZTOP(JIJ),PZZ(JIJ,JK+IKL))
    ZMIX2(JIJ) = (PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*PENTR(JIJ,JK) !&
    ZMIX3(JIJ) = (PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*PDETR(JIJ,JK) !&                   
    
    ZQTM(JIJ) = PRTM(JIJ,JK)/(1.+PRTM(JIJ,JK))            
    ZTHSM(JIJ,JK) = PTHLM(JIJ,JK)*(1.+PARAMMF%XLAMBDA_MF*ZQTM(JIJ))
    ZTHS_UP(JIJ,JK+IKL)=(ZTHS_UP(JIJ,JK)*(1.-0.5*ZMIX2(JIJ)) + &
    &ZTHSM(JIJ,JK)*ZMIX2(JIJ))&
    /(1.+0.5*ZMIX2(JIJ))   
    PRT_UP(JIJ,JK+IKL)=(PRT_UP(JIJ,JK)*(1.-0.5*ZMIX2(JIJ)) + &
    &PRTM(JIJ,JK)*ZMIX2(JIJ))  &
    /(1.+0.5*ZMIX2(JIJ))
    ZQT_UP(JIJ) = PRT_UP(JIJ,JK+IKL)/(1.+PRT_UP(JIJ,JK+IKL))
    PTHL_UP(JIJ,JK+IKL)=ZTHS_UP(JIJ,JK+IKL)/(1.+PARAMMF%XLAMBDA_MF*ZQT_UP(JIJ))                      
  ENDIF
  
  
  IF(PARAMMF%LMIXUV) THEN
    IF(JK/=IKB) THEN
      IF(GTEST(JIJ))THEN
        PU_UP(JIJ,JK+IKL) = (PU_UP(JIJ,JK)*(1-0.5*ZMIX2(JIJ)) + &
        &PUM(JIJ,JK)*ZMIX2(JIJ)+ &
        0.5*PARAMMF%XPRES_UV*(PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*&
        ((PUM(JIJ,JK+IKL)-PUM(JIJ,JK))/PDZZ(JIJ,JK+IKL)+&
        (PUM(JIJ,JK)-PUM(JIJ,JK-IKL))/PDZZ(JIJ,JK))        )   &
        /(1+0.5*ZMIX2(JIJ))
        PV_UP(JIJ,JK+IKL) = (PV_UP(JIJ,JK)*(1-0.5*ZMIX2(JIJ)) + &
        &PVM(JIJ,JK)*ZMIX2(JIJ)+ &
        0.5*PARAMMF%XPRES_UV*(PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*&
        ((PVM(JIJ,JK+IKL)-PVM(JIJ,JK))/PDZZ(JIJ,JK+IKL)+&
        (PVM(JIJ,JK)-PVM(JIJ,JK-IKL))/PDZZ(JIJ,JK))    )   &
        /(1+0.5*ZMIX2(JIJ))
      ENDIF
    ELSE
      IF(GTEST(JIJ))THEN
        PU_UP(JIJ,JK+IKL) = (PU_UP(JIJ,JK)*(1-0.5*ZMIX2(JIJ)) + &
        &PUM(JIJ,JK)*ZMIX2(JIJ)+ &
        0.5*PARAMMF%XPRES_UV*(PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*&
        ((PUM(JIJ,JK+IKL)-PUM(JIJ,JK))/PDZZ(JIJ,JK+IKL))        )   &
        /(1+0.5*ZMIX2(JIJ))
        PV_UP(JIJ,JK+IKL) = (PV_UP(JIJ,JK)*(1-0.5*ZMIX2(JIJ)) + &
        &PVM(JIJ,JK)*ZMIX2(JIJ)+ &
        0.5*PARAMMF%XPRES_UV*(PZZ(JIJ,JK+IKL)-PZZ(JIJ,JK))*&
        ((PVM(JIJ,JK+IKL)-PVM(JIJ,JK))/PDZZ(JIJ,JK+IKL))    )   &
        /(1+0.5*ZMIX2(JIJ))
      ENDIF
      
    ENDIF
  ENDIF
  !DO JSV=1,ISV 
  !   IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  !    IF(GTEST(JIJ)) THEN
  !         PSV_UP(JIJ,JK+KKL,JSV) = (PSV_UP(JIJ,JK,JSV)*(1-0.5*ZMIX2(JIJ)) + &
  !                      PSVM(JIJ,JK,JSV)*ZMIX2(JIJ))  /(1+0.5*ZMIX2(JIJ))
  !    ENDIF                        
  !ENDDO  
  
  
  ! Compute non cons. var. at level JK+KKL
  ZRC_UP(JIJ)=PRC_UP(JIJ,JK) ! guess = level just below
  ZRI_UP(JIJ)=PRI_UP(JIJ,JK) ! guess = level just below
  ZRV_UP(JIJ)=PRV_UP(JIJ,JK)
ENDDO
  CALL TH_R_FROM_THL_RT(CST,NEBN, D%NIJT, NEBN%CFRAC_ICE_SHALLOW_MF,PFRAC_ICE_UP(:,JK+IKL),ZPRES_F(:,JK+IKL), &
          PTHL_UP(:,JK+IKL),PRT_UP(:,JK+IKL),ZTH_UP(:,JK+IKL),              &
          ZRV_UP(:),ZRC_UP(:),ZRI_UP(:),ZRSATW(:),ZRSATI(:),OOCEAN=.FALSE.,&
          PBUF=ZBUF, KB=D%NIJB, KE=D%NIJE)
DO JIJ=IIJB,IIJE 
  IF(GTEST(JIJ))THEN
    ZT_UP(JIJ) = ZTH_UP(JIJ,JK+IKL)*PEXNM(JIJ,JK+IKL)
    ZCP(JIJ) = CST%XCPD + CST%XCL * ZRC_UP(JIJ)
    ZLVOCPEXN(JIJ)=(CST%XLVTT + (CST%XCPV-CST%XCL) *  (ZT_UP(JIJ)-CST%XTT) ) / &
    &ZCP(JIJ) / PEXNM(JIJ,JK+IKL)
    PRC_UP(JIJ,JK+IKL)=MIN(0.5E-3,ZRC_UP(JIJ))  ! On ne peut depasser 0.5 g/kg (autoconversion donc elimination !)
    PTHL_UP(JIJ,JK+IKL) = PTHL_UP(JIJ,JK+IKL)+ &
    & ZLVOCPEXN(JIJ)*(ZRC_UP(JIJ)-PRC_UP(JIJ,JK+IKL))
    PRV_UP(JIJ,JK+IKL)=ZRV_UP(JIJ)
    PRI_UP(JIJ,JK+IKL)=ZRI_UP(JIJ)
    PRT_UP(JIJ,JK+IKL)  = PRC_UP(JIJ,JK+IKL) + PRV_UP(JIJ,JK+IKL)
    PRSAT_UP(JIJ,JK+IKL) = ZRSATW(JIJ)*(1-PFRAC_ICE_UP(JIJ,JK+IKL)) + &
    & ZRSATI(JIJ)*PFRAC_ICE_UP(JIJ,JK+IKL)
  ENDIF
  
  
  ! Compute the updraft theta_v, buoyancy and w**2 for level JK+1   
  IF(GTEST(JIJ))THEN
    !PTHV_UP(JIJ,JK+KKL) = ZTH_UP(JIJ,JK+KKL)*((1+ZRVORD*PRV_UP(JIJ,JK+KKL))/(1+PRT_UP(JIJ,JK+KKL)))
    PTHV_UP(JIJ,JK+IKL) = ZTH_UP(JIJ,JK+IKL)* &
    & (1.+0.608*PRV_UP(JIJ,JK+IKL) - PRC_UP(JIJ,JK+IKL))
  ENDIF
  
  
  ! Test if the updraft has reach the ETL
  GTESTETL(JIJ)=.FALSE.
  IF (GTEST(JIJ).AND.(PBUO_INTEG(JIJ,JK)<=0.))THEN
    KKETL(JIJ) = JK+IKL
    GTESTETL(JIJ)=.TRUE.
  ENDIF
  
  ! Test is we have reached the top of the updraft
  IF (GTEST(JIJ).AND.((ZW_UP2(JIJ,JK+IKL)<=ZEPS)))THEN
    ZW_UP2(JIJ,JK+IKL)=ZEPS
    GTEST(JIJ)=.FALSE.
    PTHL_UP(JIJ,JK+IKL)=ZTHLM_F(JIJ,JK+IKL)
    PRT_UP(JIJ,JK+IKL)=ZRTM_F(JIJ,JK+IKL)
    PRC_UP(JIJ,JK+IKL)=0.
    PRI_UP(JIJ,JK+IKL)=0.
    PRV_UP(JIJ,JK+IKL)=0.
    PTHV_UP(JIJ,JK+IKL)=ZTHVM_F(JIJ,JK+IKL)
    PFRAC_UP(JIJ,JK+IKL)=0.
    KKCTL(JIJ)=JK+IKL
  ENDIF
ENDDO
ENDDO 

! Closure assumption for mass flux at KKB+1 level (Mass flux is supposed to be 0 at KKB level !)
!   Hourdin et al 2002 formulation


DO JIJ=IIJB,IIJE 
ZZTOP(JIJ) = MAX(ZZTOP(JIJ),ZEPS)
ENDDO

DO JK=IKB+IKL,IKE-IKL,IKL !  Vertical loop
DO JIJ=IIJB,IIJE 
  IF(JK<=IALIM(JIJ))THEN
    ZALIM_STAR_TOT(JIJ) = ZALIM_STAR_TOT(JIJ) + ZALIM_STAR(JIJ,JK)**2* &
    & ZZDZ(JIJ,JK)/PRHODREF(JIJ,JK)
  ENDIF
ENDDO
ENDDO   

DO JIJ=IIJB,IIJE 
IF (ZALIM_STAR_TOT(JIJ)*ZZTOP(JIJ) > ZEPS)THEN
  ZPHI(JIJ) =  ZW_MAX(JIJ)/(PARAMMF%XR*ZZTOP(JIJ)*ZALIM_STAR_TOT(JIJ))
ENDIF

GTEST(JIJ) = .TRUE.
PEMF(JIJ,IKB+IKL) = ZPHI(JIJ)*ZZDZ(JIJ,IKB)*ZALIM_STAR(JIJ,IKB)
! Updraft fraction must be smaller than XFRAC_UP_MAX
PFRAC_UP(JIJ,IKB+IKL)=PEMF(JIJ,IKB+IKL)/ &
                                 &(SQRT(ZW_UP2(JIJ,IKB+IKL))*ZRHO_F(JIJ,IKB+IKL))
PFRAC_UP(JIJ,IKB+IKL)=MIN(PARAMMF%XFRAC_UP_MAX,PFRAC_UP(JIJ,IKB+IKL))
PEMF(JIJ,IKB+IKL) = ZRHO_F(JIJ,IKB+IKL)*PFRAC_UP(JIJ,IKB+IKL)* &
                              & SQRT(ZW_UP2(JIJ,IKB+IKL))
ENDDO

DO JK=IKB+IKL,IKE-IKL,IKL !  Vertical loop
DO JIJ=IIJB,IIJE 
  
  GTEST(JIJ) = (ZW_UP2(JIJ,JK) > ZEPS)  
  
  IF (GTEST(JIJ))THEN
    IF(JK<IALIM(JIJ))THEN
      PEMF(JIJ,JK+IKL) = MAX(0.,PEMF(JIJ,JK) + ZPHI(JIJ)*ZZDZ(JIJ,JK)* &
      & (PENTR(JIJ,JK) - PDETR(JIJ,JK)))
    ELSE
      ZMIX1(JIJ)=ZZDZ(JIJ,JK)*(PENTR(JIJ,JK)-PDETR(JIJ,JK))
      PEMF(JIJ,JK+IKL)=PEMF(JIJ,JK)*EXP(ZMIX1(JIJ))
    ENDIF
    
! Updraft fraction must be smaller than XFRAC_UP_MAX
    PFRAC_UP(JIJ,JK+IKL)=PEMF(JIJ,JK+IKL)/&
    &(SQRT(ZW_UP2(JIJ,JK+IKL))*ZRHO_F(JIJ,JK+IKL))
    PFRAC_UP(JIJ,JK+IKL)=MIN(PARAMMF%XFRAC_UP_MAX,PFRAC_UP(JIJ,JK+IKL))
    PEMF(JIJ,JK+IKL) = ZRHO_F(JIJ,JK+IKL)*PFRAC_UP(JIJ,JK+IKL)*&
    & SQRT(ZW_UP2(JIJ,JK+IKL))
  ENDIF
ENDDO
ENDDO

DO JK=1,IKT 
DO JIJ=IIJB,IIJE 
  PW_UP(JIJ,JK)=SQRT(ZW_UP2(JIJ,JK))
ENDDO
ENDDO
DO JIJ=IIJB,IIJE 
PEMF(JIJ,IKB) =0.
ENDDO

! Limits the shallow convection scheme when cloud heigth is higher than 3000m.
! To do this, mass flux is multiplied by a coefficient decreasing linearly
! from 1 (for clouds of 3000m of depth) to 0 (for clouds of 4000m of depth).
! This way, all MF fluxes are diminished by this amount.
! Diagnosed cloud fraction is also multiplied by the same coefficient.
!
DO JIJ=IIJB,IIJE 
  PDEPTH(JIJ) = MAX(0., PZZ(JIJ,KKCTL(JIJ)) -  PZZ(JIJ,KKLCL(JIJ)) )
END DO

DO JIJ=IIJB,IIJE 
GWORK1(JIJ)= (GTESTLCL(JIJ) .AND. (PDEPTH(JIJ) > ZDEPTH_MAX1) )
ENDDO
DO JK=1,D%NKT
DO JIJ=IIJB,IIJE 
  GWORK2(JIJ,JK) = GWORK1(JIJ)
  ZCOEF(JIJ,JK) = (1.-(PDEPTH(JIJ)-ZDEPTH_MAX1)/(ZDEPTH_MAX2-ZDEPTH_MAX1))
  ZCOEF(JIJ,JK)=MIN(MAX(ZCOEF(JIJ,JK),0.),1.)
ENDDO
ENDDO
DO JK=1,IKT 
DO JIJ=IIJB,IIJE 
  IF (GWORK2(JIJ,JK)) THEN
    PEMF(JIJ,JK)     = PEMF(JIJ,JK)     * ZCOEF(JIJ,JK)
    PFRAC_UP(JIJ,JK) = PFRAC_UP(JIJ,JK) * ZCOEF(JIJ,JK)
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
