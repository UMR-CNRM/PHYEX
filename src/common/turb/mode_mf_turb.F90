!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!    ######################
     MODULE MODE_MF_TURB
!    ######################
!
IMPLICIT NONE
CONTAINS
      SUBROUTINE MF_TURB(D, KSV, PARAMMF, &
                ONOMIXLG,KSV_LGBEG,KSV_LGEND,                         &
                PTSTEP, PDZZ, PRHODJ, &
                PTHLM,PTHVM,PRTM,PUM,PVM,PSVM,PTKEM,                  &
                PTHLDT,PRTDT,PUDT,PVDT,PSVDT,PTKEDT,                  &
                PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,PSV_UP,PTKE_UP, &
                PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF,      &
                PFLXZSVMF,PFLXZTKEMF,PFRAC_UP                         )

!     #################################################################
!
!
!!****  *MF_TURB* - computes the MF_turbulent source terms for the prognostic
!!                  variables. 
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
!!  10/2009     (C.Lac)        Introduction of different PTSTEP according to the
!!                              advection schemes
!!  09/2010     (V.Masson)     Optimization
!!     S. Riette Jan 2012: support for both order of vertical levels
!!                         suppression of useless initialisations
!!     A. Marcel Jan 2025: TKE mixing
!!      A. Marcel Jan 2025: relaxation of the small fraction assumption
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
!
USE MODI_SHUMAN_MF, ONLY: MZM_MF
USE MODE_TRIDIAG_MASSFLUX, ONLY: TRIDIAG_MASSFLUX
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
INTEGER,                INTENT(IN)   :: KSV
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL,                 INTENT(IN)     ::  PTSTEP   ! Dynamical timestep 
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PDZZ        ! metric coefficients

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PRHODJ      ! dry density * Grid size

!   Conservative var. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PTHLM        ! conservative pot. temp.
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PRTM         ! water var.  where 
!  Virtual potential temperature at t-dt
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PTHVM 
!  Momentum at t-dt
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PUM
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PVM
!  scalar variables at t-dt
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVM
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) ::  PTKEM ! TKE at t-dt
!
! Tendencies of conservative variables
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) ::  PTHLDT

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) ::  PRTDT 
! Tendencies of momentum
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) ::  PUDT
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) ::  PVDT
! Tendencies of scalar variables
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT) ::  PSVDT
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) ::  PTKEDT ! Tendencies of TKE


! Updraft characteritics
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,PTKE_UP
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PFRAC_UP
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSV_UP
! Fluxes
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF,PFLXZTKEMF

REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT)::  PFLXZSVMF
!
!
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!

REAL, DIMENSION(D%NIJT,D%NKT) :: ZVARS
REAL, DIMENSION(D%NIJT,D%NKT) :: ZMEMF !-PEMF
INTEGER :: JSV          !number of scalar variables and Loop counter
INTEGER :: JIJ, JK
INTEGER :: IIJB,IIJE ! physical horizontal domain indices
INTEGER :: IKT
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRELAX_ALPHA_COEFF ! 1/(1-fraction)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
!
!*      1.PRELIMINARIES
!         -------------
!
IF (LHOOK) CALL DR_HOOK('MF_TURB',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
PFLXZSVMF(:,:,:) = 0.
PSVDT(:,:,:) = 0.
!
IF(PARAMMF%LRELAX_ALPHA_MF) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  ZRELAX_ALPHA_COEFF(IIJB:IIJE, 1:IKT)=1./(1.-PFRAC_UP(IIJB:IIJE, 1:IKT))
  ZMEMF(IIJB:IIJE, 1:IKT)=-PEMF(IIJB:IIJE, 1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE, 1:IKT)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  ZRELAX_ALPHA_COEFF(IIJB:IIJE, 1:IKT)=1.
  ZMEMF(:,:) = -PEMF(:,:)
ENDIF
!
!----------------------------------------------------------------------------
!
!*      2. COMPUTE THE MEAN FLUX OF CONSERVATIVE VARIABLES at time t-dt
!          (equation (3) of Soares et al)
!          + THE MEAN FLUX OF THETA_V (buoyancy flux)
!          -----------------------------------------------
!   ( Resulting fluxes are in flux level (w-point) as PEMF and PTHL_UP )
!

CALL MZM_MF(D, PTHLM(:,:), PFLXZTHMF(:,:))
CALL MZM_MF(D, PRTM(:,:), PFLXZRMF(:,:))
CALL MZM_MF(D, PTHVM(:,:), PFLXZTHVMF(:,:))

!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PFLXZTHMF(IIJB:IIJE,1:IKT) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                           & (PTHL_UP(IIJB:IIJE,1:IKT)-PFLXZTHMF(IIJB:IIJE,1:IKT))
PFLXZRMF(IIJB:IIJE,1:IKT) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                          & (PRT_UP(IIJB:IIJE,1:IKT)-PFLXZRMF(IIJB:IIJE,1:IKT))
PFLXZTHVMF(IIJB:IIJE,1:IKT) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                            & (PTHV_UP(IIJB:IIJE,1:IKT)-PFLXZTHVMF(IIJB:IIJE,1:IKT))
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

IF (PARAMMF%LMIXUV) THEN
  CALL MZM_MF(D, PUM(:,:), PFLXZUMF(:,:))
  CALL MZM_MF(D, PVM(:,:), PFLXZVMF(:,:))
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PFLXZUMF(IIJB:IIJE,1:IKT) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                            & (PU_UP(IIJB:IIJE,1:IKT)-PFLXZUMF(IIJB:IIJE,1:IKT))
  PFLXZVMF(IIJB:IIJE,1:IKT) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                            & (PV_UP(IIJB:IIJE,1:IKT)-PFLXZVMF(IIJB:IIJE,1:IKT))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  PFLXZUMF(:,:) = 0.
  PFLXZVMF(:,:) = 0.
ENDIF

IF (PARAMMF%LMIXTKE) THEN
  CALL MZM_MF(D, PTKEM(:,:), PFLXZTKEMF(:,:))
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PFLXZTKEMF(IIJB:IIJE,1:IKT) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                              & (PTKE_UP(IIJB:IIJE,1:IKT)-PFLXZTKEMF(IIJB:IIJE,1:IKT))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  PFLXZTKEMF(:,:) = 0.
ENDIF
!
!
!----------------------------------------------------------------------------
!
!*      3. COMPUTE TENDENCIES OF CONSERVATIVE VARIABLES (or treated as such...)
!          (implicit formulation)
!          --------------------------------------------
!

!
!
! 3.1 Compute the tendency for the conservative potential temperature
!     (PDZZ and flux in w-point and PRHODJ is mass point, result in mass point)
!
CALL TRIDIAG_MASSFLUX(D,PTHLM,PFLXZTHMF,ZMEMF,PTSTEP,PARAMMF%XIMPL_MF,  &
                      PDZZ,PRHODJ,ZVARS )
! compute new flux and THL tendency
CALL MZM_MF(D, ZVARS(:,:), PFLXZTHMF(:,:))
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PFLXZTHMF(IIJB:IIJE,1:IKT) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                           & (PTHL_UP(IIJB:IIJE,1:IKT)-PFLXZTHMF(IIJB:IIJE,1:IKT))
PTHLDT(IIJB:IIJE,1:IKT)= (ZVARS(IIJB:IIJE,1:IKT)-PTHLM(IIJB:IIJE,1:IKT))/PTSTEP
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

!
! 3.2 Compute the tendency for the conservative mixing ratio
!
CALL TRIDIAG_MASSFLUX(D,PRTM(:,:),PFLXZRMF,ZMEMF,PTSTEP,PARAMMF%XIMPL_MF,  &
                                 PDZZ,PRHODJ,ZVARS )
! compute new flux and RT tendency
CALL MZM_MF(D, ZVARS(:,:), PFLXZRMF(:,:))
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
PFLXZRMF(IIJB:IIJE,1:IKT) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                          & (PRT_UP(IIJB:IIJE,1:IKT)-PFLXZRMF(IIJB:IIJE,1:IKT))
PRTDT(IIJB:IIJE,1:IKT) = (ZVARS(IIJB:IIJE,1:IKT)-PRTM(IIJB:IIJE,1:IKT))/PTSTEP
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
!

IF (PARAMMF%LMIXUV) THEN
  !
  ! 3.3 Compute the tendency for the (non conservative but treated as it) zonal momentum
  !     (PDZZ and flux in w-point and PRHODJ is mass point, result in mass point)
  !

  CALL TRIDIAG_MASSFLUX(D,PUM,PFLXZUMF,ZMEMF,PTSTEP,PARAMMF%XIMPL_MF,  &
                                 PDZZ,PRHODJ,ZVARS )
  ! compute new flux and U tendency
  CALL MZM_MF(D, ZVARS(:,:), PFLXZUMF(:,:))
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PFLXZUMF(IIJB:IIJE,1:IKT) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                            & (PU_UP(IIJB:IIJE,1:IKT)-PFLXZUMF(IIJB:IIJE,1:IKT))
  PUDT(IIJB:IIJE,1:IKT)= (ZVARS(IIJB:IIJE,1:IKT)-PUM(IIJB:IIJE,1:IKT))/PTSTEP
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !
  !
  ! 3.4 Compute the tendency for the (non conservative but treated as it for the time beiing)
  !                                  meridian momentum
  !     (PDZZ and flux in w-point and PRHODJ is mass point, result in mass point)
  !
  CALL TRIDIAG_MASSFLUX(D,PVM,PFLXZVMF,ZMEMF,PTSTEP,PARAMMF%XIMPL_MF,  &
                                 PDZZ,PRHODJ,ZVARS )
  ! compute new flux and V tendency
  CALL MZM_MF(D, ZVARS(:,:), PFLXZVMF(:,:))
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PFLXZVMF(IIJB:IIJE,1:IKT) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                            & (PV_UP(IIJB:IIJE,1:IKT)-PFLXZVMF(IIJB:IIJE,1:IKT))
  PVDT(IIJB:IIJE,1:IKT)= (ZVARS(IIJB:IIJE,1:IKT)-PVM(IIJB:IIJE,1:IKT))/PTSTEP
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  PUDT(:,:)=0.
  PVDT(:,:)=0.
ENDIF

DO JSV=1,KSV 

  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  
  !*     compute mean flux of scalar variables at time t-dt
  !   ( Resulting fluxes are in flux level (w-point) as PEMF and PTHL_UP )

  CALL MZM_MF(D, PSVM(:,:,JSV), PFLXZSVMF(:,:,JSV))
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PFLXZSVMF(IIJB:IIJE,1:IKT,JSV) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                                 & (PSV_UP(IIJB:IIJE,1:IKT,JSV)-PFLXZSVMF(IIJB:IIJE,1:IKT,JSV))
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  !
  ! 3.5 Compute the tendency for scalar variables
  !     (PDZZ and flux in w-point and PRHODJ is mass point, result in mass point)
  !
  CALL TRIDIAG_MASSFLUX(D,PSVM(:,:,JSV),PFLXZSVMF(:,:,JSV),&
                        ZMEMF,PTSTEP,PARAMMF%XIMPL_MF,PDZZ,PRHODJ,ZVARS )
  ! compute new flux and Sv tendency
  CALL MZM_MF(D, ZVARS, PFLXZSVMF(:,:,JSV))
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
  PFLXZSVMF(IIJB:IIJE,1:IKT,JSV) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                                 & (PSV_UP(IIJB:IIJE,1:IKT,JSV)-PFLXZSVMF(IIJB:IIJE,1:IKT,JSV))
  PSVDT(IIJB:IIJE,1:IKT,JSV)= (ZVARS(IIJB:IIJE,1:IKT)-PSVM(IIJB:IIJE,1:IKT,JSV))/PTSTEP
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

ENDDO
!
!
! 3.6 Compute the tendency for the TKE
!     (PDZZ and flux in w-point and PRHODJ is mass point, result in mass point)
!
IF (PARAMMF%LMIXTKE) THEN                                                                                                          
  CALL TRIDIAG_MASSFLUX(D,PTKEM,PFLXZTKEMF,ZMEMF,PTSTEP,PARAMMF%XIMPL_MF,PDZZ,PRHODJ,ZVARS)

  ! compute new flux and TKE tendency
  CALL MZM_MF(D, ZVARS(:,:), PFLXZTKEMF(:,:))
  PFLXZTKEMF(IIJB:IIJE,1:IKT) = PEMF(IIJB:IIJE,1:IKT)*ZRELAX_ALPHA_COEFF(IIJB:IIJE,1:IKT)*&
                              & (PTKE_UP(IIJB:IIJE,1:IKT)-PFLXZTKEMF(IIJB:IIJE,1:IKT))
  PTKEDT(IIJB:IIJE,1:IKT)= (ZVARS(IIJB:IIJE,1:IKT)-PTKEM(IIJB:IIJE,1:IKT))/PTSTEP
ELSE
  PTKEDT(:,:)=0.
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MF_TURB',1,ZHOOK_HANDLE)
END SUBROUTINE MF_TURB    
END MODULE MODE_MF_TURB
