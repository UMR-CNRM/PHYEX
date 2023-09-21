!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!    ######################
     MODULE MODE_MF_TURB_EXPL
!    ######################
IMPLICIT NONE
CONTAINS
      SUBROUTINE MF_TURB_EXPL(D, PARAMMF,                             &
                PRHODJ,PTHLM,PTHVM,PRTM,PUM,PVM,                      &
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
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODI_SHUMAN_MF, ONLY: MZM_MF

IMPLICIT NONE


!*      0.1  declarations of arguments


TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PRHODJ      ! dry density * Grid size

!   Conservative var. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PTHLM        ! conservative pot. temp.
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PRTM         ! water var.  where 

!  Virtual potential temperature at t-dt
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PTHVM 
!  Momentum at t-dt
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PUM
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PVM
!
! Tendencies of conservative variables
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) ::  PTHLDT

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) ::  PRTDT 

! Tendencies of momentum
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) ::  PUDT
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) ::  PVDT

! Updraft characteritics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP

! Fluxes
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PFLXZTHLMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF

REAL, DIMENSION(D%NIJT,D%NKT) :: ZFLXZTHSMF,ZTHS_UP,ZTHSM  ! Theta S flux
REAL, DIMENSION(D%NIJT,D%NKT) :: ZQT_UP,ZQTM,ZTHSDT,ZQTDT
REAL, DIMENSION(D%NIJT,D%NKT) :: ZTHLM_F,ZRTM_F

INTEGER                              :: JK, JIJ            ! loop counter
INTEGER :: IIJB,IIJE ! physical horizontal domain indices
INTEGER :: IKT,IKB,IKE,IKL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!----------------------------------------------------------------------------
!
!*      1.PRELIMINARIES
!         -------------

IF (LHOOK) CALL DR_HOOK('MF_TURB_EXPL',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
IKB=D%NKB
IKE=D%NKE
IKL=D%NKL
!
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
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZQTM(:,:)   = ZRTM_F(:,:)/(1.+ZRTM_F(:,:))
ZQT_UP(:,:) = PRT_UP(:,:)/(1.+PRT_UP(:,:))
ZTHS_UP(:,:)= PTHL_UP(:,:)*(1.+PARAMMF%XLAMBDA_MF*ZQT_UP(:,:))
ZTHSM(:,:)  = ZTHLM_F(:,:)*(1.+PARAMMF%XLAMBDA_MF*ZQTM(:,:))
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

CALL MZM_MF(D, PTHLM(:,:), PFLXZTHLMF(:,:))
CALL MZM_MF(D, PRTM(:,:), PFLXZRMF(:,:))
CALL MZM_MF(D, PTHVM(:,:), PFLXZTHVMF(:,:))
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PFLXZTHLMF(:,:)  = PEMF(:,:)*(PTHL_UP(:,:)-PFLXZTHLMF(:,:))  ! ThetaL
PFLXZRMF(:,:)    = PEMF(:,:)*(PRT_UP(:,:)-PFLXZRMF(:,:))  ! Rt
PFLXZTHVMF(:,:)  = PEMF(:,:)*(PTHV_UP(:,:)-PFLXZTHVMF(:,:))  ! ThetaV

ZFLXZTHSMF(:,:)  = PEMF(:,:)*(ZTHS_UP(:,:)-ZTHSM(:,:))    ! Theta S flux
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

IF (PARAMMF%LMIXUV) THEN
  CALL MZM_MF(D, PUM(:,:), PFLXZUMF(:,:))
  CALL MZM_MF(D, PVM(:,:), PFLXZVMF(:,:))
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PFLXZUMF(:,:) =  PEMF(:,:)*(PU_UP(:,:)-PFLXZUMF(:,:))  ! U
  PFLXZVMF(:,:) =  PEMF(:,:)*(PV_UP(:,:)-PFLXZVMF(:,:))  ! V
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  PFLXZUMF(:,:) = 0.
  PFLXZVMF(:,:) = 0.
ENDIF


!----------------------------------------------------------------------------
!
!*      3. COMPUTE TENDENCIES OF CONSERVATIVE VARIABLES (or treated as such...)
!          (explicit formulation)
!          --------------------------------------------

DO JK=IKB,IKE-IKL,IKL
  DO JIJ=IIJB,IIJE
    !PTHLDT(JIJ,JK) = (PFLXZTHLMF(JIJ,JK  ) - PFLXZTHLMF(JIJ,JK+IKL)) / PRHODJ(JIJ,JK)
    PRTDT(JIJ,JK) = (PFLXZRMF(JIJ,JK) - PFLXZRMF(JIJ,JK+IKL)) / PRHODJ(JIJ,JK)
    ZQTDT(JIJ,JK) = PRTDT(JIJ,JK)/(1.+ ZRTM_F(JIJ,JK)*ZRTM_F(JIJ,JK))
    ZTHSDT(JIJ,JK)= (ZFLXZTHSMF(JIJ,JK) - ZFLXZTHSMF(JIJ,JK+IKL)) / PRHODJ(JIJ,JK)
    PTHLDT(JIJ,JK) = ZTHSDT(JIJ,JK)/(1.+PARAMMF%XLAMBDA_MF*ZQTM(JIJ,JK)) - ZTHLM_F(JIJ,JK)*PARAMMF%XLAMBDA_MF*ZQTDT(JIJ,JK)
  ENDDO
END DO

IF (PARAMMF%LMIXUV) THEN
  DO JK=IKB,IKE-IKL,IKL
    DO JIJ=IIJB,IIJE
      PUDT(JIJ,JK) = (PFLXZUMF(JIJ,JK) - PFLXZUMF(JIJ,JK+IKL)) / PRHODJ(JIJ,JK)
      PVDT(JIJ,JK) = (PFLXZVMF(JIJ,JK) - PFLXZVMF(JIJ,JK+IKL)) / PRHODJ(JIJ,JK)
    ENDDO
  END DO
ENDIF  


IF (LHOOK) CALL DR_HOOK('MF_TURB_EXPL',1,ZHOOK_HANDLE)
END SUBROUTINE MF_TURB_EXPL
END MODULE MODE_MF_TURB_EXPL
