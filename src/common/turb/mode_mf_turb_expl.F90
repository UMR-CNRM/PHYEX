!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!    ######################
     MODULE MODE_MF_TURB_EXPL
!    ######################
IMPLICIT NONE
CONTAINS
      SUBROUTINE MF_TURB_EXPL(D, PARAMMF, OMIXUV,                     &
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
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODI_SHUMAN_MF, ONLY: MZM_MF

IMPLICIT NONE


!*      0.1  declarations of arguments


TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF
LOGICAL,                INTENT(IN)   :: OMIXUV      ! True if mixing of momentum

REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)   :: PRHODJ      ! dry density * Grid size

!   Conservative var. at t-dt
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) ::  PTHLM        ! conservative pot. temp.
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) ::  PRTM         ! water var.  where 

!  Virtual potential temperature at t-dt
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) ::  PTHVM 
!  Momentum at t-dt
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) ::  PUM
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) ::  PVM
!
! Tendencies of conservative variables
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT) ::  PTHLDT

REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT) ::  PRTDT 

! Tendencies of momentum
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT) ::  PUDT
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT) ::  PVDT

! Updraft characteritics
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)   ::  PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP

! Fluxes
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  ::  PFLXZTHLMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF

REAL, DIMENSION(D%NIT,D%NKT) :: ZFLXZTHSMF,ZTHS_UP,ZTHSM  ! Theta S flux
REAL, DIMENSION(D%NIT,D%NKT) :: ZQT_UP,ZQTM,ZTHSDT,ZQTDT
REAL, DIMENSION(D%NIT,D%NKT) :: ZTHLM_F,ZRTM_F

INTEGER                              :: JK, JI            ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!----------------------------------------------------------------------------
!
!*      1.PRELIMINARIES
!         -------------

IF (LHOOK) CALL DR_HOOK('MF_TURB_EXPL',0,ZHOOK_HANDLE)

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

CALL MZM_MF(D, PRTM (:,:), ZRTM_F(:,:))
CALL MZM_MF(D, PTHLM(:,:), ZTHLM_F(:,:))
!$mnh_expand_array(JI=D%NIB:D%NIE,JK=D%NKTB:D%NKTE)
ZQTM(D%NIB:D%NIE,:)   = ZRTM_F(D%NIB:D%NIE,:)/(1.+ZRTM_F(D%NIB:D%NIE,:))
ZQT_UP(D%NIB:D%NIE,:) = PRT_UP(D%NIB:D%NIE,:)/(1.+PRT_UP(D%NIB:D%NIE,:))
ZTHS_UP(D%NIB:D%NIE,:)= PTHL_UP(D%NIB:D%NIE,:)*(1.+PARAMMF%XLAMBDA_MF*ZQT_UP(D%NIB:D%NIE,:))
ZTHSM(D%NIB:D%NIE,:)  = ZTHLM_F(D%NIB:D%NIE,:)*(1.+PARAMMF%XLAMBDA_MF*ZQTM(D%NIB:D%NIE,:))
!$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=D%NKTB:D%NKTE)

CALL MZM_MF(D, PTHLM(:,:), PFLXZTHLMF(:,:))
CALL MZM_MF(D, PRTM(:,:), PFLXZRMF(:,:))
CALL MZM_MF(D, PTHVM(:,:), PFLXZTHVMF(:,:))
!$mnh_expand_array(JI=D%NIB:D%NIE,JK=D%NKTB:D%NKTE)
PFLXZTHLMF(D%NIB:D%NIE,:)  = PEMF(D%NIB:D%NIE,:)*(PTHL_UP(D%NIB:D%NIE,:)-PFLXZTHLMF(D%NIB:D%NIE,:))  ! ThetaL
PFLXZRMF(D%NIB:D%NIE,:)    = PEMF(D%NIB:D%NIE,:)*(PRT_UP(D%NIB:D%NIE,:)-PFLXZRMF(D%NIB:D%NIE,:))  ! Rt
PFLXZTHVMF(D%NIB:D%NIE,:)  = PEMF(D%NIB:D%NIE,:)*(PTHV_UP(D%NIB:D%NIE,:)-PFLXZTHVMF(D%NIB:D%NIE,:))  ! ThetaV

ZFLXZTHSMF(D%NIB:D%NIE,:)  = PEMF(D%NIB:D%NIE,:)*(ZTHS_UP(D%NIB:D%NIE,:)-ZTHSM(D%NIB:D%NIE,:))    ! Theta S flux
!$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=D%NKTB:D%NKTE)

IF (OMIXUV) THEN
  CALL MZM_MF(D, PUM(:,:), PFLXZUMF(:,:))
  CALL MZM_MF(D, PVM(:,:), PFLXZVMF(:,:))
  !$mnh_expand_array(JI=D%NIB:D%NIE,JK=D%NKTB:D%NKTE)
  PFLXZUMF(D%NIB:D%NIE,:) =  PEMF(D%NIB:D%NIE,:)*(PU_UP(D%NIB:D%NIE,:)-PFLXZUMF(D%NIB:D%NIE,:))  ! U
  PFLXZVMF(D%NIB:D%NIE,:) =  PEMF(D%NIB:D%NIE,:)*(PV_UP(D%NIB:D%NIE,:)-PFLXZVMF(D%NIB:D%NIE,:))  ! V
  !$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=D%NKTB:D%NKTE)
ELSE
  PFLXZUMF(:,:) = 0.
  PFLXZVMF(:,:) = 0.
ENDIF


!----------------------------------------------------------------------------
!
!*      3. COMPUTE TENDENCIES OF CONSERVATIVE VARIABLES (or treated as such...)
!          (explicit formulation)
!          --------------------------------------------

DO JK=D%NKB,D%NKE-D%NKL,D%NKL
  DO JI=D%NIB,D%NIE
    !PTHLDT(JI,JK) = (PFLXZTHLMF(JI,JK  ) - PFLXZTHLMF(JI,JK+D%NKL)) / PRHODJ(JI,JK)
    PRTDT(JI,JK) = (PFLXZRMF(JI,JK) - PFLXZRMF(JI,JK+D%NKL)) / PRHODJ(JI,JK)
    ZQTDT(JI,JK) = PRTDT(JI,JK)/(1.+ ZRTM_F(JI,JK)*ZRTM_F(JI,JK))
    ZTHSDT(JI,JK)= (ZFLXZTHSMF(JI,JK) - ZFLXZTHSMF(JI,JK+D%NKL)) / PRHODJ(JI,JK)
    PTHLDT(JI,JK) = ZTHSDT(JI,JK)/(1.+PARAMMF%XLAMBDA_MF*ZQTM(JI,JK)) - ZTHLM_F(JI,JK)*PARAMMF%XLAMBDA_MF*ZQTDT(JI,JK)
  ENDDO
END DO

IF (OMIXUV) THEN
  DO JK=D%NKB,D%NKE-D%NKL,D%NKL
    DO JI=D%NIB,D%NIE
      PUDT(JI,JK) = (PFLXZUMF(JI,JK) - PFLXZUMF(JI,JK+D%NKL)) / PRHODJ(JI,JK)
      PVDT(JI,JK) = (PFLXZVMF(JI,JK) - PFLXZVMF(JI,JK+D%NKL)) / PRHODJ(JI,JK)
    ENDDO
  END DO
ENDIF  


IF (LHOOK) CALL DR_HOOK('MF_TURB_EXPL',1,ZHOOK_HANDLE)
END SUBROUTINE MF_TURB_EXPL
END MODULE MODE_MF_TURB_EXPL
