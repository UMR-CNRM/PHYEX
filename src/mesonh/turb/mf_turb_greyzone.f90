!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
     MODULE MODI_MF_TURB_GREYZONE
!    ######################
!
INTERFACE
!     #################################################################
      SUBROUTINE MF_TURB_GREYZONE(KKA,KKB,KKE,KKU,KKL,OMIXUV,                  &
                ONOMIXLG,KSV_LGBEG,KSV_LGEND,                         &
                PIMPL, PTSTEP,                              &
                PDZZ,                                                 &
                PRHODJ,                                               &
                PTHLM,PTHVM,PRTM,PUM,PVM,PSVM,                        &
                PTHLDT,PRTDT,PUDT,PVDT,PSVDT,                         &
                PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,PSV_UP,       &
                     PTHL_DO,PTHV_DO,PRT_DO,PU_DO,PV_DO,PSV_DO,       &
                PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF,      &
                PFLXZSVMF                                             )

!     #################################################################
!
!               
!*               1.1  Declaration of Arguments
!
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise

LOGICAL,                INTENT(IN)   :: OMIXUV      ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL,                   INTENT(IN)   :: PIMPL       ! degree of implicitness
REAL,                 INTENT(IN)     ::  PTSTEP   ! Dynamical timestep 
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PDZZ        ! metric coefficients

REAL, DIMENSION(:,:), INTENT(IN)   ::  PRHODJ    ! dry density * Grid size

!   Conservative var. at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHLM       ! conservative pot. temp.
REAL, DIMENSION(:,:), INTENT(IN) ::  PRTM         ! water var.  where 
!  Virtual potential temperature at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHVM 
!  Momentum at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PUM
REAL, DIMENSION(:,:), INTENT(IN) ::  PVM
!  scalar variables at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PSVM
!
! Tendencies of conservative variables
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PTHLDT

REAL, DIMENSION(:,:),   INTENT(OUT) ::  PRTDT 
! Tendencies of momentum
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PUDT
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PVDT
! Tendencies of scalar variables
REAL, DIMENSION(:,:,:), INTENT(OUT) ::  PSVDT


! Updraft characteritics
REAL, DIMENSION(:,:), INTENT(IN)   ::  PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP
REAL, DIMENSION(:,:), INTENT(IN)   ::       PTHL_DO,PTHV_DO,PRT_DO,PU_DO,PV_DO
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PSV_UP, PSV_DO
! Fluxes
REAL, DIMENSION(:,:), INTENT(OUT)  ::  PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF

REAL, DIMENSION(:,:,:), INTENT(OUT)::  PFLXZSVMF

END SUBROUTINE MF_TURB_GREYZONE

END INTERFACE
!
END MODULE MODI_MF_TURB_GREYZONE
!     #################################################################
      SUBROUTINE MF_TURB_GREYZONE(KKA, KKB, KKE, KKU, KKL,OMIXUV,                  &
                ONOMIXLG,KSV_LGBEG,KSV_LGEND,                         &
                PIMPL, PTSTEP,                                        &
                PDZZ,                                                 &
                PRHODJ,                                               &
                PTHLM,PTHVM,PRTM,PUM,PVM,PSVM,                        &
                PTHLDT,PRTDT,PUDT,PVDT,PSVDT,                         &
                PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP,PSV_UP,       &
                     PTHL_DO,PTHV_DO,PRT_DO,PU_DO,PV_DO,PSV_DO,       &
                PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF,      &
                PFLXZSVMF                                             )

!     #################################################################
!
!
!!****  *MF_TURB_GREYZONE* - computes the MF_turbulent source terms for the prognostic
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
USE MODD_PARAM_MFSHALL_n
!
USE MODI_SHUMAN_MF
USE MODI_TRIDIAG_MASSFLUX
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
LOGICAL,                INTENT(IN)   :: OMIXUV      ! True if mixing of momentum
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer
REAL,                   INTENT(IN)   :: PIMPL       ! degree of implicitness
REAL,                 INTENT(IN)     ::  PTSTEP   ! Dynamical timestep 
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PDZZ        ! metric coefficients

REAL, DIMENSION(:,:), INTENT(IN)   :: PRHODJ      ! dry density * Grid size

!   Conservative var. at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHLM        ! conservative pot. temp.
REAL, DIMENSION(:,:), INTENT(IN) ::  PRTM         ! water var.  where 
!  Virtual potential temperature at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHVM 
!  Momentum at t-dt
REAL, DIMENSION(:,:), INTENT(IN) ::  PUM
REAL, DIMENSION(:,:), INTENT(IN) ::  PVM
!  scalar variables at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PSVM
!
! Tendencies of conservative variables
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PTHLDT

REAL, DIMENSION(:,:),   INTENT(OUT) ::  PRTDT 
! Tendencies of momentum
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PUDT
REAL, DIMENSION(:,:),   INTENT(OUT) ::  PVDT
! Tendencies of scalar variables
REAL, DIMENSION(:,:,:), INTENT(OUT) ::  PSVDT


! Updraft/environment characteritics
REAL, DIMENSION(:,:), INTENT(IN)   ::  PEMF,PTHL_UP,PTHV_UP,PRT_UP,PU_UP,PV_UP
REAL, DIMENSION(:,:), INTENT(IN)   ::  PTHL_DO,PTHV_DO,PRT_DO,PU_DO,PV_DO
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PSV_UP,PSV_DO
! Fluxes
REAL, DIMENSION(:,:), INTENT(OUT)  ::  PFLXZTHMF,PFLXZTHVMF,PFLXZRMF,PFLXZUMF,PFLXZVMF

REAL, DIMENSION(:,:,:), INTENT(OUT)::  PFLXZSVMF
!
!
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!

REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2)) :: ZVARS

!
INTEGER :: ISV,JSV          !number of scalar variables and Loop counter
!
!----------------------------------------------------------------------------
!
!*      1.PRELIMINARIES
!         -------------
!
!
! number of scalar var
ISV=SIZE(PSVM,3)

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
! downdraft data are on the flux points
PFLXZTHMF(:,:) = PEMF(:,:)*(PTHL_UP(:,:)-PTHL_DO(:,:))

PFLXZRMF(:,:) =  PEMF(:,:)*(PRT_UP(:,:)-PRT_DO(:,:))

PFLXZTHVMF(:,:) = PEMF(:,:)*(PTHV_UP(:,:)-PTHV_DO(:,:))

IF (OMIXUV) THEN
  PFLXZUMF(:,:) =  PEMF(:,:)*(PU_UP(:,:)-PU_DO(:,:))
  PFLXZVMF(:,:) =  PEMF(:,:)*(PV_UP(:,:)-PV_DO(:,:))
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
CALL TRIDIAG_MASSFLUX(KKA,KKB,KKE,KKU,KKL,PTHLM,PFLXZTHMF,-PEMF,PTSTEP,PIMPL,  &
                      PDZZ,PRHODJ,ZVARS )
! compute new flux
!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pourquoi on le recalcule ici alors qu'il n'est pas utilisé ailleurs
! sauf pour l'écriture ?
! Est ce que ZVARS est au point de masse pour qu'il doivent être remis au point
! de flux ?
!!!!!!!!!!!!!!!!!!!!!!!!!!
PFLXZTHMF(:,:) = PEMF(:,:)*(PTHL_UP(:,:)-MZM_MF(KKA,KKU,KKL,ZVARS(:,:)))

!!! compute THL tendency
!
PTHLDT(:,:)= (ZVARS(:,:)-PTHLM(:,:))/PTSTEP

!
! 3.2 Compute the tendency for the conservative mixing ratio
!
CALL TRIDIAG_MASSFLUX(KKA,KKB,KKE,KKU,KKL,PRTM(:,:),PFLXZRMF,-PEMF,PTSTEP,PIMPL,  &
                                 PDZZ,PRHODJ,ZVARS )
! compute new flux
PFLXZRMF(:,:) =  PEMF(:,:)*(PRT_UP(:,:)-MZM_MF(KKA,KKU,KKL,ZVARS(:,:)))

!!! compute RT tendency
PRTDT(:,:) = (ZVARS(:,:)-PRTM(:,:))/PTSTEP
!

IF (OMIXUV) THEN
  !
  ! 3.3 Compute the tendency for the (non conservative but treated as it) zonal momentum
  !     (PDZZ and flux in w-point and PRHODJ is mass point, result in mass point)
  !

  CALL TRIDIAG_MASSFLUX(KKA,KKB,KKE,KKU,KKL,PUM,PFLXZUMF,-PEMF,PTSTEP,PIMPL,  &
                                 PDZZ,PRHODJ,ZVARS )
  ! compute new flux
  PFLXZUMF(:,:) = PEMF(:,:)*(PU_UP(:,:)-MZM_MF(KKA,KKU,KKL,ZVARS(:,:)))

  ! compute U tendency
  PUDT(:,:)= (ZVARS(:,:)-PUM(:,:))/PTSTEP

  !
  !
  ! 3.4 Compute the tendency for the (non conservative but treated as it for the time beiing)
  !                                  meridian momentum
  !     (PDZZ and flux in w-point and PRHODJ is mass point, result in mass point)
  !
  CALL TRIDIAG_MASSFLUX(KKA,KKB,KKE,KKU,KKL,PVM,PFLXZVMF,-PEMF,PTSTEP,PIMPL,  &
                                 PDZZ,PRHODJ,ZVARS )
  ! compute new flux
  PFLXZVMF(:,:) = PEMF(:,:)*(PV_UP(:,:)-MZM_MF(KKA,KKU,KKL,ZVARS(:,:)))

  ! compute V tendency
  PVDT(:,:)= (ZVARS(:,:)-PVM(:,:))/PTSTEP
ELSE
  PUDT(:,:)=0.
  PVDT(:,:)=0.
ENDIF

DO JSV=1,ISV 

  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  
  !*     compute mean flux of scalar variables at time t-dt
  !   ( Resulting fluxes are in flux level (w-point) as PEMF and PTHL_UP )

  PFLXZSVMF(:,:,JSV) = PEMF(:,:)*(PSV_UP(:,:,JSV)-MZM_MF(KKA,KKU,KKL,PSVM(:,:,JSV)))
  
  !
  ! 3.5 Compute the tendency for scalar variables
  !     (PDZZ and flux in w-point and PRHODJ is mass point, result in mass point)
  !
  CALL TRIDIAG_MASSFLUX(KKA,KKB,KKE,KKU,KKL,PSVM(:,:,JSV),PFLXZSVMF(:,:,JSV),&
                        -PEMF,PTSTEP,PIMPL,PDZZ,PRHODJ,ZVARS )
  ! compute new flux
  PFLXZSVMF(:,:,JSV) = PEMF(:,:)*(PSV_UP(:,:,JSV)-MZM_MF(KKA,KKU,KKL,ZVARS))

  ! compute Sv tendency
  PSVDT(:,:,JSV)= (ZVARS(:,:)-PSVM(:,:,JSV))/PTSTEP

ENDDO
!
END SUBROUTINE MF_TURB_GREYZONE    
