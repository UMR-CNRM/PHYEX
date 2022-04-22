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
      SUBROUTINE MF_TURB(D, KSV, OMIXUV,                  &
                ONOMIXLG,KSV_LGBEG,KSV_LGEND,                         &
                PIMPL, PTSTEP,                                        &
                PDZZ,                                                 &
                PRHODJ,                                               &
                PTHLM,PTHVM,PRTM,PUM,PVM,PSVM,                        &
                PTHLDT,PRTDT,PUDT,PVDT,PSVDT,                         &
                PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,PSV_UP,       &
                PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF,      &
                PFLXZSVMF                                             )

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
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
!
USE MODI_SHUMAN_MF, ONLY: MZM_MF
USE MODE_TRIDIAG_MASSFLUX, ONLY: TRIDIAG_MASSFLUX
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
INTEGER,                INTENT(IN)   :: KSV
LOGICAL,                INTENT(IN)   :: OMIXUV      ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL,                   INTENT(IN)   :: PIMPL       ! degree of implicitness
REAL,                 INTENT(IN)     ::  PTSTEP   ! Dynamical timestep 
!
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)   :: PDZZ        ! metric coefficients

REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)   :: PRHODJ      ! dry density * Grid size

!   Conservative var. at t-dt
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) ::  PTHLM        ! conservative pot. temp.
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) ::  PRTM         ! water var.  where 
!  Virtual potential temperature at t-dt
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) ::  PTHVM 
!  Momentum at t-dt
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) ::  PUM
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) ::  PVM
!  scalar variables at t-dt
REAL, DIMENSION(D%NIT,D%NKT,KSV), INTENT(IN) ::  PSVM
!
! Tendencies of conservative variables
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT) ::  PTHLDT

REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT) ::  PRTDT 
! Tendencies of momentum
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT) ::  PUDT
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT) ::  PVDT
! Tendencies of scalar variables
REAL, DIMENSION(D%NIT,D%NKT,KSV), INTENT(OUT) ::  PSVDT


! Updraft characteritics
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)   ::  PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP
REAL, DIMENSION(D%NIT,D%NKT,KSV), INTENT(IN) ::  PSV_UP
! Fluxes
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  ::  PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF

REAL, DIMENSION(D%NIT,D%NKT,KSV), INTENT(OUT)::  PFLXZSVMF
!
!
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!

REAL, DIMENSION(D%NIT,D%NKT) :: ZVARS
INTEGER :: JSV          !number of scalar variables and Loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
!
!*      1.PRELIMINARIES
!         -------------
!
IF (LHOOK) CALL DR_HOOK('MF_TURB',0,ZHOOK_HANDLE)
!
PFLXZSVMF = 0.
PSVDT = 0.

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
PFLXZTHMF(:,:) = PEMF(:,:)*(PTHL_UP(:,:)-PFLXZTHMF(:,:))

CALL MZM_MF(D, PRTM(:,:), PFLXZRMF(:,:))
PFLXZRMF(:,:) =  PEMF(:,:)*(PRT_UP(:,:)-PFLXZRMF(:,:))

CALL MZM_MF(D, PTHVM(:,:), PFLXZTHVMF(:,:))
PFLXZTHVMF(:,:) = PEMF(:,:)*(PTHV_UP(:,:)-PFLXZTHVMF(:,:))

IF (OMIXUV) THEN
  CALL MZM_MF(D, PUM(:,:), PFLXZUMF(:,:))
  PFLXZUMF(:,:) =  PEMF(:,:)*(PU_UP(:,:)-PFLXZUMF(:,:))
  CALL MZM_MF(D, PVM(:,:), PFLXZVMF(:,:))
  PFLXZVMF(:,:) =  PEMF(:,:)*(PV_UP(:,:)-PFLXZVMF(:,:))
ELSE
  PFLXZUMF(:,:) = 0.
  PFLXZVMF(:,:) = 0.
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
CALL TRIDIAG_MASSFLUX(D,PTHLM,PFLXZTHMF,-PEMF,PTSTEP,PIMPL,  &
                      PDZZ,PRHODJ,ZVARS )
! compute new flux
CALL MZM_MF(D, ZVARS(:,:), PFLXZTHMF(:,:))
PFLXZTHMF(:,:) = PEMF(:,:)*(PTHL_UP(:,:)-PFLXZTHMF(:,:))

!!! compute THL tendency
!
PTHLDT(:,:)= (ZVARS(:,:)-PTHLM(:,:))/PTSTEP

!
! 3.2 Compute the tendency for the conservative mixing ratio
!
CALL TRIDIAG_MASSFLUX(D,PRTM(:,:),PFLXZRMF,-PEMF,PTSTEP,PIMPL,  &
                                 PDZZ,PRHODJ,ZVARS )
! compute new flux
CALL MZM_MF(D, ZVARS(:,:), PFLXZRMF(:,:))
PFLXZRMF(:,:) =  PEMF(:,:)*(PRT_UP(:,:)-PFLXZRMF(:,:))

!!! compute RT tendency
PRTDT(:,:) = (ZVARS(:,:)-PRTM(:,:))/PTSTEP
!

IF (OMIXUV) THEN
  !
  ! 3.3 Compute the tendency for the (non conservative but treated as it) zonal momentum
  !     (PDZZ and flux in w-point and PRHODJ is mass point, result in mass point)
  !

  CALL TRIDIAG_MASSFLUX(D,PUM,PFLXZUMF,-PEMF,PTSTEP,PIMPL,  &
                                 PDZZ,PRHODJ,ZVARS )
  ! compute new flux
  CALL MZM_MF(D, ZVARS(:,:), PFLXZUMF(:,:))
  PFLXZUMF(:,:) = PEMF(:,:)*(PU_UP(:,:)-PFLXZUMF(:,:))

  ! compute U tendency
  PUDT(:,:)= (ZVARS(:,:)-PUM(:,:))/PTSTEP

  !
  !
  ! 3.4 Compute the tendency for the (non conservative but treated as it for the time beiing)
  !                                  meridian momentum
  !     (PDZZ and flux in w-point and PRHODJ is mass point, result in mass point)
  !
  CALL TRIDIAG_MASSFLUX(D,PVM,PFLXZVMF,-PEMF,PTSTEP,PIMPL,  &
                                 PDZZ,PRHODJ,ZVARS )
  ! compute new flux
  CALL MZM_MF(D, ZVARS(:,:), PFLXZVMF(:,:))
  PFLXZVMF(:,:) = PEMF(:,:)*(PV_UP(:,:)-PFLXZVMF(:,:))

  ! compute V tendency
  PVDT(:,:)= (ZVARS(:,:)-PVM(:,:))/PTSTEP
ELSE
  PUDT(:,:)=0.
  PVDT(:,:)=0.
ENDIF

DO JSV=1,KSV 

  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  
  !*     compute mean flux of scalar variables at time t-dt
  !   ( Resulting fluxes are in flux level (w-point) as PEMF and PTHL_UP )

  CALL MZM_MF(D, PSVM(:,:,JSV), PFLXZSVMF(:,:,JSV))
  PFLXZSVMF(:,:,JSV) = PEMF(:,:)*(PSV_UP(:,:,JSV)-PFLXZSVMF(:,:,JSV))
  
  !
  ! 3.5 Compute the tendency for scalar variables
  !     (PDZZ and flux in w-point and PRHODJ is mass point, result in mass point)
  !
  CALL TRIDIAG_MASSFLUX(D,PSVM(:,:,JSV),PFLXZSVMF(:,:,JSV),&
                        -PEMF,PTSTEP,PIMPL,PDZZ,PRHODJ,ZVARS )
  ! compute new flux
  CALL MZM_MF(D, ZVARS, PFLXZSVMF(:,:,JSV))
  PFLXZSVMF(:,:,JSV) = PEMF(:,:)*(PSV_UP(:,:,JSV)-PFLXZSVMF(:,:,JSV))

  ! compute Sv tendency
  PSVDT(:,:,JSV)= (ZVARS(:,:)-PSVM(:,:,JSV))/PTSTEP

ENDDO
!
IF (LHOOK) CALL DR_HOOK('MF_TURB',1,ZHOOK_HANDLE)
END SUBROUTINE MF_TURB    
END MODULE MODE_MF_TURB
