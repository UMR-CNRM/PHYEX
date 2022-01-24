!     ######spl
      SUBROUTINE SHALLOW_MF(KKA,KKU,KKL,KRR,KRRL,KRRI,                &
                HMF_UPDRAFT, HMF_CLOUD, HFRAC_ICE, OMIXUV,            &
                ONOMIXLG,KSV_LGBEG,KSV_LGEND,                         &
                PIMPL_MF, PTSTEP, PTSTEP_MET, PTSTEP_SV,              &
                PDZZ, PZZ,                                            &
                PRHODJ, PRHODREF,                                     &
                PPABSM, PEXNM,                                        &
                PSFTH,PSFRV,                                          &
                PTHM,PRM,PUM,PVM,PTKEM,PSVM,                          &
                PDUDT_MF,PDVDT_MF,                                    &
                PDTHLDT_MF,PDRTDT_MF,PDSVDT_MF,                       &
                PSIGMF,PRC_MF,PRI_MF,PCF_MF,PFLXZTHVMF,               &
                PFLXZTHMF,PFLXZRMF,PFLXZUMF,PFLXZVMF,                 &
                PTHL_UP,PRT_UP,PRV_UP,PRC_UP,PRI_UP,                  &
                PU_UP, PV_UP, PTHV_UP, PW_UP,                         &
                PFRAC_UP,PEMF,PDETR,PENTR,                            &
                KKLCL,KKETL,KKCTL                                     )

      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #################################################################
!!
!!****  *SHALLOW_MF* - 
!!       
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is
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
!!     REFERENCE
!!     ---------
!!
!!
!!     AUTHOR
!!     ------
!!     J.Pergaud
!!
!!    MODIFICATIONS
!!    -------------
!!      Original
!!      V.Masson 09/2010 : optimization
!!      S. Riette 18 May 2010 interface changed due to ice correction
!!      S.Riette DUAL case
!!      S. Riette Jan 2012: support for both order of vertical levels
!!      S. Riette Nov 2016: HFRAC_ICE support
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_PARAMETERS, ONLY: JPVEXT
USE MODD_CMFSHALL

USE MODE_THL_RT_FROM_TH_R_MF, ONLY: THL_RT_FROM_TH_R_MF
USE MODI_COMPUTE_UPDRAFT
USE MODI_COMPUTE_UPDRAFT_RHCJ10
USE MODI_COMPUTE_UPDRAFT_RAHA
USE MODI_MF_TURB
USE MODI_MF_TURB_EXPL
USE MODI_COMPUTE_MF_CLOUD
USE MODI_COMPUTE_FRAC_ICE
!
IMPLICIT NONE

!*                    0.1  Declaration of Arguments
!
!
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER,                INTENT(IN)   :: KRR          ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL         ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI         ! number of ice water var.
CHARACTER (LEN=4),      INTENT(IN)   :: HMF_UPDRAFT  ! Type of Mass Flux Scheme
                                     ! 'NONE' if no parameterization 
CHARACTER (LEN=4),      INTENT(IN)   :: HMF_CLOUD    ! Type of statistical cloud
                                                     ! scheme
CHARACTER*1,            INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
LOGICAL,                INTENT(IN)   :: OMIXUV    ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL,                   INTENT(IN)   :: PIMPL_MF     ! degre of implicitness
REAL,              INTENT(IN)     ::  PTSTEP   ! Dynamical timestep 
REAL,              INTENT(IN)     ::  PTSTEP_MET! Timestep for meteorological variables                        
REAL,              INTENT(IN)     ::  PTSTEP_SV! Timestep for tracer variables

REAL, DIMENSION(:,:),   INTENT(IN) ::  PZZ         ! Height of flux point
REAL, DIMENSION(:,:),   INTENT(IN) ::  PDZZ        ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN) ::  PRHODJ      ! dry density * Grid size
REAL, DIMENSION(:,:),   INTENT(IN) ::  PRHODREF    ! dry density of the
                                                     ! reference state
REAL, DIMENSION(:,:),   INTENT(IN) ::  PPABSM      ! Pressure at time t-1
REAL, DIMENSION(:,:),   INTENT(IN) ::  PEXNM       ! Exner function at t-dt

REAL, DIMENSION(:),   INTENT(IN)   ::  PSFTH,PSFRV ! normal surface fluxes of theta and Rv 
REAL, DIMENSION(:,:), INTENT(IN)   ::  PTHM        ! Theta at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PRM         ! water var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN) ::  PUM,PVM     ! wind components at t-dt
REAL, DIMENSION(:,:),   INTENT(IN) ::  PTKEM       ! tke at t-dt

REAL, DIMENSION(:,:,:), INTENT(IN) ::  PSVM        ! scalar variable a t-dt

REAL, DIMENSION(:,:),   INTENT(OUT)::  PDUDT_MF     ! tendency of U   by massflux scheme
REAL, DIMENSION(:,:),   INTENT(OUT)::  PDVDT_MF     ! tendency of V   by massflux scheme
REAL, DIMENSION(:,:),   INTENT(OUT)::  PDTHLDT_MF   ! tendency of thl by massflux scheme
REAL, DIMENSION(:,:),   INTENT(OUT)::  PDRTDT_MF    ! tendency of rt  by massflux scheme
REAL, DIMENSION(:,:,:), INTENT(OUT)::  PDSVDT_MF    ! tendency of Sv  by massflux scheme

REAL, DIMENSION(:,:), INTENT(OUT)     ::  PSIGMF,PRC_MF,PRI_MF,PCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(:,:), INTENT(OUT)     ::  PFLXZTHVMF           ! Thermal production for TKE scheme
REAL, DIMENSION(:,:), INTENT(OUT)     ::  PFLXZTHMF
REAL, DIMENSION(:,:), INTENT(OUT)     ::  PFLXZRMF
REAL, DIMENSION(:,:), INTENT(OUT)     ::  PFLXZUMF
REAL, DIMENSION(:,:), INTENT(OUT)     ::  PFLXZVMF
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PTHL_UP   ! Thl updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PRT_UP    ! Rt  updraft characteristics
REAL, DIMENSION(:,:), INTENT(OUT) ::  PRV_UP    ! Vapor updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PU_UP     ! U wind updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PV_UP     ! V wind updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PRC_UP    ! cloud content updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PRI_UP    ! ice content   updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PTHV_UP   ! Thv   updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PW_UP     ! vertical speed updraft characteristics
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PFRAC_UP  ! updraft fraction
REAL, DIMENSION(:,:), INTENT(INOUT) ::  PEMF      ! updraft mass flux
REAL, DIMENSION(:,:), INTENT(OUT) ::  PDETR     ! updraft detrainment
REAL, DIMENSION(:,:), INTENT(OUT) ::  PENTR     ! updraft entrainment
INTEGER,DIMENSION(:), INTENT(OUT) :: KKLCL,KKETL,KKCTL ! level of LCL,ETL and CTL
!
!                     0.2  Declaration of local variables
!
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) ::     &
          ZTHLM,                                  & !
          ZRTM,                                   & !
          ZTHVM,                                  & !
          ZEMF_O_RHODREF,                         & ! entrainment/detrainment
          ZTHVDT,ZTHDT,ZRVDT,                     & ! tendencies
          ZBUO_INTEG                                ! integrated buoyancy
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZFRAC_ICE

REAL, DIMENSION(SIZE(PSVM,1),SIZE(PSVM,2),SIZE(PSVM,3)) ::  &
                                          ZSV_UP,&  ! updraft scalar var.
                                          ZFLXZSVMF ! Flux     
REAL, DIMENSION(SIZE(PTHM,1)) :: ZDEPTH             ! Deepness of cloud
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZFRAC_ICE_UP ! liquid/solid fraction in updraft
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZRSAT_UP ! Rsat in updraft

LOGICAL :: GENTR_DETR  ! flag to recompute entrainment, detrainment and mass flux
INTEGER :: IKB         ! near ground physical index
INTEGER :: IKE         ! uppest atmosphere physical index
!------------------------------------------------------------------------

!!! 1. Initialisation

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SHALLOW_MF',0,ZHOOK_HANDLE)

! vertical boundaries
IKB=KKA+KKL*JPVEXT
IKE=KKU-KKL*JPVEXT

! updraft governing variables
IF (HMF_UPDRAFT == 'EDKF') THEN
  PENTR      = 1.E20
  PDETR      = 1.E20
  PEMF       = 1.E20
  ZBUO_INTEG = 1.E20
ENDIF

! Thermodynamics functions
ZFRAC_ICE(:,:) = 0.
WHERE(PRM(:,:,2)+PRM(:,:,4) > 1.E-20)
  ZFRAC_ICE(:,:) = PRM(:,:,4) / (PRM(:,:,2)+PRM(:,:,4))
ENDWHERE
CALL COMPUTE_FRAC_ICE(HFRAC_ICE,ZFRAC_ICE(:,:),PTHM(:,:)*PEXNM(:,:))

! Conservative variables at t-dt
CALL THL_RT_FROM_TH_R_MF(KRR,KRRL,KRRI,    &
                         PTHM, PRM, PEXNM, &
                         ZTHLM, ZRTM       )

! Virtual potential temperature at t-dt
ZTHVM(:,:) = PTHM(:,:)*((1.+XRV / XRD *PRM(:,:,1))/(1.+ZRTM(:,:))) 

! 
!!! 2. Compute updraft
!!!    ---------------
!
IF (HMF_UPDRAFT == 'EDKF') THEN
  GENTR_DETR = .TRUE.
  CALL COMPUTE_UPDRAFT(KKA,IKB,IKE,KKU,KKL,HFRAC_ICE,GENTR_DETR,OMIXUV,&
                       ONOMIXLG,KSV_LGBEG,KSV_LGEND,             &
                       PZZ,PDZZ,                                 &
                       PSFTH,PSFRV,PPABSM,PRHODREF,              &
                       PUM,PVM,PTKEM,                            &
                       PTHM,PRM(:,:,1),ZTHLM,ZRTM,PSVM,          &
                       PTHL_UP,PRT_UP,PRV_UP,PRC_UP,PRI_UP,      &
                       PTHV_UP, PW_UP, PU_UP, PV_UP, ZSV_UP,     &
                       PFRAC_UP,ZFRAC_ICE_UP,ZRSAT_UP,PEMF,PDETR,&
                       PENTR,ZBUO_INTEG,KKLCL,KKETL,KKCTL,ZDEPTH )
ELSEIF (HMF_UPDRAFT == 'RHCJ') THEN
   CALL COMPUTE_UPDRAFT_RHCJ10(KKA,IKB,IKE,KKU,KKL,HFRAC_ICE,GENTR_DETR,OMIXUV,&
                       ONOMIXLG,KSV_LGBEG,KSV_LGEND,             &
                       PZZ,PDZZ,                                 &
                       PSFTH,PSFRV,PPABSM,PRHODREF,              &
                       PUM,PVM,PTKEM,                            &
                       PEXNM,PTHM,PRM(:,:,1),ZTHLM,ZRTM,PSVM,    &
                       PTHL_UP,PRT_UP,PRV_UP,PRC_UP,PRI_UP,      &
                       PTHV_UP, PW_UP, PU_UP, PV_UP, ZSV_UP,     &
                       PFRAC_UP,ZFRAC_ICE_UP,ZRSAT_UP,PEMF,PDETR,&
                       PENTR,ZBUO_INTEG,KKLCL,KKETL,KKCTL,ZDEPTH )
ELSEIF (HMF_UPDRAFT == 'RAHA') THEN
   CALL COMPUTE_UPDRAFT_RAHA(KKA,IKB,IKE,KKU,KKL,HFRAC_ICE,GENTR_DETR,OMIXUV,&
                       ONOMIXLG,KSV_LGBEG,KSV_LGEND,             &
                       PZZ,PDZZ,                                 &
                       PSFTH,PSFRV,PPABSM,PRHODREF,              &
                       PUM,PVM,PTKEM,                            &
                       PEXNM,PTHM,PRM(:,:,1),ZTHLM,ZRTM,PSVM,    &
                       PTHL_UP,PRT_UP,PRV_UP,PRC_UP,PRI_UP,      &
                       PTHV_UP, PW_UP, PU_UP, PV_UP, ZSV_UP,     &
                       PFRAC_UP,ZFRAC_ICE_UP,ZRSAT_UP,PEMF,PDETR,&
                       PENTR,ZBUO_INTEG,KKLCL,KKETL,KKCTL,ZDEPTH )
ELSEIF (HMF_UPDRAFT == 'DUAL') THEN
  !Updraft characteristics are already computed and received by interface
ELSE
  WRITE(*,*) ' STOP'                                                     
  WRITE(*,*) ' NO UPDRAFT MODEL FOR EDKF : CMF_UPDRAFT =',HMF_UPDRAFT 
  CALL ABORT
  STOP
ENDIF

!!! 5. Compute diagnostic convective cloud fraction and content
!!!    --------------------------------------------------------
!
CALL COMPUTE_MF_CLOUD(KKA,IKB,IKE,KKU,KKL,KRR,KRRL,KRRI,&
                      HMF_CLOUD,ZFRAC_ICE,              &
                      PRC_UP,PRI_UP,PEMF,               &
                      PTHL_UP,PRT_UP,PFRAC_UP,          &
                      PTHV_UP,ZFRAC_ICE_UP,             &
                      ZRSAT_UP,PEXNM,ZTHLM,ZRTM,        &
                      PTHM, ZTHVM, PRM,                 &
                      PDZZ,PZZ,KKLCL,                   &
                      PPABSM,PRHODREF,                  &
                      PRC_MF,PRI_MF,PCF_MF,PSIGMF,ZDEPTH)


!!! 3. Compute fluxes of conservative variables and their divergence = tendency
!!!    ------------------------------------------------------------------------
!
ZEMF_O_RHODREF=PEMF/PRHODREF

IF ( PIMPL_MF > 1.E-10 ) THEN  
CALL MF_TURB(KKA, IKB, IKE, KKU, KKL, OMIXUV,                         &
             ONOMIXLG,KSV_LGBEG,KSV_LGEND,                            &
             PIMPL_MF, PTSTEP, PTSTEP_MET, PTSTEP_SV,                 &
             PDZZ,                                                    &
             PRHODJ,                                                  &
             ZTHLM,ZTHVM,ZRTM,PUM,PVM,PSVM,                           &
             PDTHLDT_MF,PDRTDT_MF,PDUDT_MF,PDVDT_MF,PDSVDT_MF,        &
             ZEMF_O_RHODREF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,ZSV_UP,&
             PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF,         &
             ZFLXZSVMF                                                )
ELSE
  CALL MF_TURB_EXPL(KKA, IKB, IKE, KKU, KKL, OMIXUV,                     &
           PRHODJ,                                                       &
           ZTHLM,ZTHVM,ZRTM,PUM,PVM,                                     &
           PDTHLDT_MF,PDRTDT_MF,PDUDT_MF,PDVDT_MF,                       &
           ZEMF_O_RHODREF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,            &
           PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF)
ENDIF

! security in the case HMF_UPDRAFT = 'DUAL'
! to be modified if 'DUAL' is evolving (momentum mixing for example)
IF( HMF_UPDRAFT == 'DUAL') THEN
  ! Now thetav_up from vdfhghtnn is used!
  PFLXZTHVMF=0.
  ! Yes/No UV mixing!
!  PDUDT_MF=0.
!  PDVDT_MF=0.
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SHALLOW_MF',1,ZHOOK_HANDLE)
END SUBROUTINE SHALLOW_MF
