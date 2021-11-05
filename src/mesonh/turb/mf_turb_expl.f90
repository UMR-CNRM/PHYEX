!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!    ######################
     MODULE MODI_MF_TURB_EXPL
!    ######################
!
INTERFACE
!
!     #################################################################
      SUBROUTINE MF_TURB_EXPL(KKA,KKB,KKE,KKU,KKL,OMIXUV,             &
                PRHODJ,                                               &
                PTHLM,PTHVM,PRTM,PUM,PVM,                             &
                PTHLDT,PRTDT,PUDT,PVDT,                               &
                PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,              &
                PFLXZTHLMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF)
!     #################################################################
!
!*                    1.1  Declaration of Arguments
!
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
LOGICAL,                INTENT(IN)   :: OMIXUV      ! True if mixing of momentum

REAL, DIMENSION(:,:), INTENT(IN)   :: PRHODJ      ! dry density * Grid size

!   Conservative var. at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHLM        ! conservative pot. temp.
REAL, DIMENSION(:,:), INTENT(IN) ::  PRTM         ! water var.  where 

!  Virtual potential temperature at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHVM 
!  Momentum at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PUM
REAL, DIMENSION(:,:), INTENT(IN) ::  PVM
!
! Tendencies of conservative variables
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PTHLDT

REAL, DIMENSION(:,:),   INTENT(OUT) ::  PRTDT 

! Tendencies of momentum
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PUDT
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PVDT

! Updraft characteritics
REAL, DIMENSION(:,:), INTENT(IN)   ::  PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP

! Fluxes
REAL, DIMENSION(:,:), INTENT(OUT)  ::  PFLXZTHLMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF

END SUBROUTINE MF_TURB_EXPL

END INTERFACE
!
END MODULE MODI_MF_TURB_EXPL
!

!     ######spl
      SUBROUTINE MF_TURB_EXPL(KKA,KKB,KKE,KKU,KKL,OMIXUV,             &
                PRHODJ,                                               &
                PTHLM,PTHVM,PRTM,PUM,PVM,                             &
                PTHLDT,PRTDT,PUDT,PVDT,                               &
                PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,              &
                PFLXZTHLMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF)

!     #################################################################
!
!
!!****  *MF_TURB_EXPL* - computes the MF_turbulent source terms for the prognostic
!!                       variables (when PIMPL=0)
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to compute the source terms in 
!!    the evolution equations due to the MF turbulent mixing. 
!!      The source term is computed as the divergence of the turbulent fluxes.
!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!     
!!
!!    MODIFICATIONS
!!    -------------
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------

USE MODD_PARAM_MFSHALL_n
USE MODI_SHUMAN_MF

IMPLICIT NONE


!*      0.1  declarations of arguments


INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
LOGICAL,                INTENT(IN)   :: OMIXUV      ! True if mixing of momentum

REAL, DIMENSION(:,:), INTENT(IN)   :: PRHODJ      ! dry density * Grid size

!   Conservative var. at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHLM        ! conservative pot. temp.
REAL, DIMENSION(:,:), INTENT(IN) ::  PRTM         ! water var.  where 

!  Virtual potential temperature at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHVM 
!  Momentum at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PUM
REAL, DIMENSION(:,:), INTENT(IN) ::  PVM
!
! Tendencies of conservative variables
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PTHLDT

REAL, DIMENSION(:,:),   INTENT(OUT) ::  PRTDT 

! Tendencies of momentum
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PUDT
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PVDT

! Updraft characteritics
REAL, DIMENSION(:,:), INTENT(IN)   ::  PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP

! Fluxes
REAL, DIMENSION(:,:), INTENT(OUT)  ::  PFLXZTHLMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF

REAL, DIMENSION(SIZE(PFLXZTHLMF,1),SIZE(PFLXZTHLMF,2)) :: ZFLXZTHSMF,ZTHS_UP,ZTHSM  ! Theta S flux
REAL, DIMENSION(SIZE(PFLXZTHLMF,1),SIZE(PFLXZTHLMF,2)) :: ZQT_UP,ZQTM,ZTHSDT,ZQTDT
REAL, DIMENSION(SIZE(PFLXZTHLMF,1),SIZE(PFLXZTHLMF,2)) :: ZTHLM_F,ZRTM_F

INTEGER                              :: JK            ! loop counter

!----------------------------------------------------------------------------
!
!*      1.PRELIMINARIES
!         -------------

PFLXZRMF   = 0.
PFLXZTHVMF = 0.
PFLXZTHLMF = 0.
PFLXZUMF   = 0.
PFLXZVMF   = 0.
PTHLDT = 0.
PRTDT  = 0.
PUDT   = 0.
PVDT   = 0.

!
!----------------------------------------------------------------------------
!
!*      2. COMPUTE THE MEAN FLUX OF CONSERVATIVE VARIABLES at time t-dt
!          (equation (3) of Soares et al)
!          + THE MEAN FLUX OF THETA_V (buoyancy flux)
!          -----------------------------------------------
!   ( Resulting fluxes are in flux level (w-point) as PEMF and PTHL_UP )

ZRTM_F (:,:) = MZM_MF(KKA,KKU,KKL,PRTM (:,:))
ZTHLM_F(:,:) = MZM_MF(KKA,KKU,KKL,PTHLM(:,:))
ZQTM   (:,:) = ZRTM_F (:,:)/(1.+ZRTM_F (:,:))
ZQT_UP (:,:) = PRT_UP (:,:)/(1.+PRT_UP (:,:))
ZTHS_UP(:,:) = PTHL_UP(:,:)*(1.+XLAMBDA_MF*ZQT_UP(:,:))
ZTHSM  (:,:) = ZTHLM_F(:,:)*(1.+XLAMBDA_MF*ZQTM(:,:))

PFLXZTHLMF(:,:)  = PEMF(:,:)*(PTHL_UP(:,:)-MZM_MF(KKA,KKU,KKL,PTHLM(:,:)))  ! ThetaL
PFLXZRMF(:,:)    = PEMF(:,:)*(PRT_UP (:,:)-MZM_MF(KKA,KKU,KKL,PRTM (:,:)))  ! Rt
PFLXZTHVMF(:,:)  = PEMF(:,:)*(PTHV_UP(:,:)-MZM_MF(KKA,KKU,KKL,PTHVM(:,:)))  ! ThetaV

ZFLXZTHSMF(:,:)  = PEMF(:,:)*(ZTHS_UP(:,:)-ZTHSM(:,:))    ! Theta S flux

IF (OMIXUV) THEN
  PFLXZUMF(:,:) =  PEMF(:,:)*(PU_UP(:,:)-MZM_MF(KKA,KKU,KKL,PUM(:,:)))  ! U
  PFLXZVMF(:,:) =  PEMF(:,:)*(PV_UP(:,:)-MZM_MF(KKA,KKU,KKL,PVM(:,:)))  ! V
ELSE
  PFLXZUMF(:,:) = 0.
  PFLXZVMF(:,:) = 0.
ENDIF


!----------------------------------------------------------------------------
!
!*      3. COMPUTE TENDENCIES OF CONSERVATIVE VARIABLES (or treated as such...)
!          (explicit formulation)
!          --------------------------------------------

DO JK=KKB,KKE-KKL,KKL
!  PTHLDT(:,JK) = (PFLXZTHLMF(:,JK  ) - PFLXZTHLMF(:,JK+KKL)) / PRHODJ(:,JK)
  PRTDT (:,JK) = (PFLXZRMF  (:,JK  ) - PFLXZRMF  (:,JK+KKL)) / PRHODJ(:,JK)
  ZQTDT (:,JK) = PRTDT (:,JK)/(1.+ ZRTM_F (:,JK)*ZRTM_F (:,JK))
  ZTHSDT(:,JK) = (ZFLXZTHSMF(:,JK  ) - ZFLXZTHSMF(:,JK+KKL)) / PRHODJ(:,JK)
  PTHLDT(:,JK) = ZTHSDT(:,JK)/(1.+XLAMBDA_MF*ZQTM(:,JK)) - ZTHLM_F(:,JK)*XLAMBDA_MF*ZQTDT(:,JK)
END DO

IF (OMIXUV) THEN
  DO JK=KKB,KKE-KKL,KKL
    PUDT(:,JK) = (PFLXZUMF(:,JK  ) - PFLXZUMF(:,JK+KKL)) / PRHODJ(:,JK)
    PVDT(:,JK) = (PFLXZVMF(:,JK  ) - PFLXZVMF(:,JK+KKL)) / PRHODJ(:,JK)
  END DO
ENDIF  


END SUBROUTINE MF_TURB_EXPL
