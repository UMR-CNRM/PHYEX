!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ################################################################
      SUBROUTINE SHALLOW_MF(D, CST, NEB, PARAMMF, TURB, CSTURB,       &
                KRR, KRRL, KRRI, KSV,                                 &
                HMF_UPDRAFT, HMF_CLOUD, HFRAC_ICE, OMIXUV, OSTATNW,   &
                ONOMIXLG,KSV_LGBEG,KSV_LGEND,                         &
                PIMPL_MF, PTSTEP,                                     &
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
                KKLCL,KKETL,KKCTL,PDX,PDY                             )

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
!!      R.Honnert 07/2012 : elemnts of Rio according to Bouteloup
!!      R.Honnert 07/2012 : MF gray zone 
!!      R.Honnert 10/2016 : SURF=gray zone initilisation + EDKF  
!!      R.Honnert 10/2016 : Update with Arome
!!      S. Riette Nov 2016: HFRAC_ICE support
!!      Philippe Wautelet 28/05/2018: corrected truncated integer division (2/3 -> 2./3.)
!!      Q.Rodier  01/2019 : support RM17 mixing length
!!      R.Honnert 1/2019  : remove SURF 
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  R. Honnert     04/2021: remove HRIO and BOUT schemes
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
USE MODD_CTURB,           ONLY: CSTURB_t
!
USE MODE_THL_RT_FROM_TH_R_MF, ONLY: THL_RT_FROM_TH_R_MF
USE MODE_COMPUTE_UPDRAFT, ONLY: COMPUTE_UPDRAFT
USE MODE_COMPUTE_UPDRAFT_RHCJ10, ONLY: COMPUTE_UPDRAFT_RHCJ10
USE MODE_COMPUTE_UPDRAFT_RAHA, ONLY: COMPUTE_UPDRAFT_RAHA
USE MODE_MF_TURB, ONLY: MF_TURB
USE MODE_MF_TURB_EXPL, ONLY: MF_TURB_EXPL
USE MODE_COMPUTE_MF_CLOUD, ONLY: COMPUTE_MF_CLOUD
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE

!*                    0.1  Declaration of Arguments
!
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(NEB_t),            INTENT(IN)   :: NEB
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF
TYPE(TURB_t),           INTENT(IN)   :: TURB
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
INTEGER,                INTENT(IN)   :: KRR          ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL         ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI         ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV
CHARACTER(LEN=4),       INTENT(IN)   :: HMF_UPDRAFT  ! Type of Mass Flux Scheme
                                     ! 'NONE' if no parameterization 
CHARACTER(LEN=4),       INTENT(IN)   :: HMF_CLOUD    ! Type of statistical cloud
                                                     ! scheme
CHARACTER(LEN=1),       INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
LOGICAL,                INTENT(IN)   :: OSTATNW      ! cloud scheme inclues convect. covar. contrib
LOGICAL,                INTENT(IN)   :: OMIXUV    ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL,                   INTENT(IN)   :: PIMPL_MF     ! degre of implicitness
REAL,                   INTENT(IN)   :: PTSTEP    ! Dynamical timestep 

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PZZ         ! Height of flux point
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PDZZ        ! Metric coefficients
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PRHODJ      ! dry density * Grid size
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PRHODREF    ! dry density of the
                                                     ! reference state
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PPABSM      ! Pressure at time t-1
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PEXNM       ! Exner function at t-dt

REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSFTH,PSFRV ! normal surface fluxes of theta and Rv 
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTHM        ! Theta at t-dt
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN) ::  PRM         ! water var. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PUM,PVM     ! wind components at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PTKEM       ! tke at t-dt

REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVM        ! scalar variable a t-dt

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)::  PDUDT_MF     ! tendency of U   by massflux scheme
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)::  PDVDT_MF     ! tendency of V   by massflux scheme
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)::  PDTHLDT_MF   ! tendency of thl by massflux scheme
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)::  PDRTDT_MF    ! tendency of rt  by massflux scheme
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT)::  PDSVDT_MF    ! tendency of Sv  by massflux scheme

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PSIGMF,PRC_MF,PRI_MF,PCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PFLXZTHVMF           ! Thermal production for TKE scheme
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PFLXZTHMF
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PFLXZRMF
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PFLXZUMF
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     ::  PFLXZVMF
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PTHL_UP   ! Thl updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PRT_UP    ! Rt  updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PRV_UP    ! Vapor updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PU_UP     ! U wind updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PV_UP     ! V wind updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PRC_UP    ! cloud content updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PRI_UP    ! ice content   updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PTHV_UP   ! Thv   updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PW_UP     ! vertical speed updraft characteristics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PFRAC_UP  ! updraft fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) ::  PEMF      ! updraft mass flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) ::  PDETR     ! updraft detrainment
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) ::  PENTR     ! updraft entrainment
INTEGER,DIMENSION(D%NIJT), INTENT(OUT) :: KKLCL,KKETL,KKCTL ! level of LCL,ETL and CTL
REAL,                 INTENT(IN)  :: PDX, PDY
!
!                     0.2  Declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT) ::     &
          ZTHLM,                                  & !
          ZRTM,                                   & !
          ZTHVM,                                  & !
          ZEMF_O_RHODREF,                         & ! entrainment/detrainment
          ZBUO_INTEG                                ! integrated buoyancy
REAL, DIMENSION(D%NIJT,D%NKT) :: ZFRAC_ICE
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWK

REAL, DIMENSION(D%NIJT,D%NKT,KSV) ::  &
                                          ZSV_UP,&  ! updraft scalar var.
                                          ZFLXZSVMF ! Flux     
REAL, DIMENSION(D%NIJT) :: ZDEPTH             ! Deepness of cloud
REAL, DIMENSION(D%NIJT,D%NKT) :: ZFRAC_ICE_UP ! liquid/solid fraction in updraft
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRSAT_UP ! Rsat in updraft

LOGICAL :: GENTR_DETR  ! flag to recompute entrainment, detrainment and mass flux
INTEGER, DIMENSION(D%NIJT,D%NKT) :: IERR
INTEGER :: JI, JK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------

!!! 1. Initialisation
IF (LHOOK) CALL DR_HOOK('SHALLOW_MF',0,ZHOOK_HANDLE)

! updraft governing variables
IF (HMF_UPDRAFT == 'EDKF'  .OR. HMF_UPDRAFT == 'RHCJ') THEN
  PENTR      = 1.E20
  PDETR      = 1.E20
  PEMF       = 1.E20
  ZBUO_INTEG = 1.E20
ENDIF

! Thermodynamics functions
ZFRAC_ICE(:,:) = 0.
IF (KRR.GE.4) THEN
  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
      IF(PRM(JI,JK,2)+PRM(JI,JK,4) > 1.E-20)THEN
        ZFRAC_ICE(JI,JK) = PRM(JI,JK,4) / (PRM(JI,JK,2)+PRM(JI,JK,4))
      ENDIF
    ENDDO
  ENDDO
ENDIF
DO JK=1,D%NKT 
  DO JI=D%NIJB,D%NIJE 
    ZWK(JI,JK)=PTHM(JI,JK)*PEXNM(JI,JK)
  ENDDO
ENDDO
CALL COMPUTE_FRAC_ICE(HFRAC_ICE,NEB,ZFRAC_ICE(:,:),ZWK(:,:), IERR(:,:))

! Conservative variables at t-dt
CALL THL_RT_FROM_TH_R_MF(D, CST, KRR,KRRL,KRRI,    &
                         PTHM, PRM, PEXNM, &
                         ZTHLM, ZRTM       )

! Virtual potential temperature at t-dt
DO JK=1,D%NKT 
  DO JI=D%NIJB,D%NIJE 
    ZTHVM(JI,JK) = PTHM(JI,JK)*&
    & ((1.+CST%XRV / CST%XRD *PRM(JI,JK,1))/(1.+ZRTM(JI,JK))) 
  ENDDO
ENDDO
! 
!!! 2. Compute updraft
!!!    ---------------
!
IF (HMF_UPDRAFT == 'EDKF') THEN
  GENTR_DETR = .TRUE.
  CALL COMPUTE_UPDRAFT(D, CST, NEB, PARAMMF, TURB, CSTURB,       &
                       KSV, HFRAC_ICE, GENTR_DETR, OMIXUV,       &
                       ONOMIXLG,KSV_LGBEG,KSV_LGEND,             &
                       PZZ,PDZZ,                                 &
                       PSFTH,PSFRV,PPABSM,PRHODREF,              &
                       PUM,PVM,PTKEM,                            &
                       PTHM,PRM(:,:,1),ZTHLM,ZRTM,PSVM,          &
                       PTHL_UP,PRT_UP,PRV_UP,PRC_UP,PRI_UP,      &
                       PTHV_UP, PW_UP, PU_UP, PV_UP, ZSV_UP,     &
                       PFRAC_UP,ZFRAC_ICE_UP,ZRSAT_UP,PEMF,PDETR,&
                       PENTR,ZBUO_INTEG,KKLCL,KKETL,KKCTL,ZDEPTH,&
                       PDX,PDY)
ELSEIF (HMF_UPDRAFT == 'RHCJ') THEN
  GENTR_DETR = .TRUE.
  CALL COMPUTE_UPDRAFT_RHCJ10(D, CST, NEB, PARAMMF, TURB, CSTURB,&
                       KSV, HFRAC_ICE, GENTR_DETR, OMIXUV,       &
                       ONOMIXLG,KSV_LGBEG,KSV_LGEND,             &
                       PZZ,PDZZ,                                 &
                       PSFTH,PSFRV,PPABSM,PRHODREF,              &
                       PUM,PVM,PTKEM,                            &
                       PTHM,PRM(:,:,1),ZTHLM,ZRTM,PSVM,          &
                       PTHL_UP,PRT_UP,PRV_UP,PRC_UP,PRI_UP,      &
                       PTHV_UP, PW_UP, PU_UP, PV_UP, ZSV_UP,     &
                       PFRAC_UP,ZFRAC_ICE_UP,ZRSAT_UP,PEMF,PDETR,&
                       PENTR,ZBUO_INTEG,KKLCL,KKETL,KKCTL,ZDEPTH )
ELSEIF (HMF_UPDRAFT == 'RAHA') THEN
   CALL COMPUTE_UPDRAFT_RAHA(D, CST, NEB, PARAMMF,               &
                       KSV, HFRAC_ICE, GENTR_DETR, OMIXUV,       &
                       ONOMIXLG,KSV_LGBEG,KSV_LGEND,             &
                       PZZ,PDZZ,                                 &
                       PSFTH,PSFRV,                              &
                       PPABSM,PRHODREF,PUM,PVM,PTKEM,            &
                       PEXNM,PTHM,PRM(:,:,1),ZTHLM,ZRTM,         &
                       PSVM,PTHL_UP,PRT_UP,                      &
                       PRV_UP,PRC_UP,PRI_UP, PTHV_UP,            &
                       PW_UP, PU_UP, PV_UP, ZSV_UP,              &
                       PFRAC_UP,ZFRAC_ICE_UP,ZRSAT_UP,           &
                       PEMF,PDETR,PENTR,                         &
                       ZBUO_INTEG,KKLCL,KKETL,KKCTL,             &
                       ZDEPTH )
ELSEIF (HMF_UPDRAFT == 'DUAL') THEN
  !Updraft characteristics are already computed and received by interface
ELSE
  CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'SHALLOW_MF', 'no updraft model for EDKF: CMF_UPDRAFT='//TRIM(HMF_UPDRAFT) )
ENDIF

!!! 5. Compute diagnostic convective cloud fraction and content
!!!    --------------------------------------------------------
!
CALL COMPUTE_MF_CLOUD(D, CST, CSTURB, PARAMMF, OSTATNW, &
                      KRR, KRRL, KRRI,                  &
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
DO JK=1,D%NKT 
  DO JI=D%NIJB,D%NIJE 
    ZEMF_O_RHODREF(JI,JK)=PEMF(JI,JK)/PRHODREF(JI,JK)
  ENDDO
ENDDO

IF ( PIMPL_MF > 1.E-10 ) THEN  
  CALL MF_TURB(D, KSV, OMIXUV,                     &
             ONOMIXLG,KSV_LGBEG,KSV_LGEND,                            &
             PIMPL_MF, PTSTEP,                                        &
             PDZZ,                                                    &
             PRHODJ,                                                  &
             ZTHLM,ZTHVM,ZRTM,PUM,PVM,PSVM,                           &
             PDTHLDT_MF,PDRTDT_MF,PDUDT_MF,PDVDT_MF,PDSVDT_MF,        &
             ZEMF_O_RHODREF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,ZSV_UP,&
             PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF,         &
             ZFLXZSVMF                                                )
ELSE
  CALL MF_TURB_EXPL(D, PARAMMF, OMIXUV,                 &
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
!
CONTAINS
INCLUDE "compute_frac_ice.func.h"
!
END SUBROUTINE SHALLOW_MF
