!MNH_LIC Copyright 2022-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################
      MODULE MODI_ION_SOURCE_ELEC
!     ###########################
!
INTERFACE
      SUBROUTINE ION_SOURCE_ELEC (KTCOUNT, KRR, HLBCX, HLBCY, &
                                  PRHODREF, PRHODJ, PRT,      &
                                  PSVT, PSVS,                 &
                                  PEFIELDU, PEFIELDV, PEFIELDW)
!
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! Temporal loop counter
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ   ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF ! Reference dry air density
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT      ! Moist variables at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT     ! Scalar variable at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS     ! Scalar variable sources
!
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PEFIELDU  !  Electric field
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PEFIELDV  !    components
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PEFIELDW  ! along x, y and z
!
END SUBROUTINE ION_SOURCE_ELEC
END INTERFACE
END MODULE MODI_ION_SOURCE_ELEC
!
!     #########################################################
      SUBROUTINE ION_SOURCE_ELEC (KTCOUNT, KRR, HLBCX, HLBCY, &
                                  PRHODREF, PRHODJ, PRT,      &
                                  PSVT, PSVS,                 &
                                  PEFIELDU, PEFIELDV, PEFIELDW)
!     #########################################################
!!
!!****  * -  compute the ion source from drift motion and cosmic rays
!!
!!    AUTHOR
!!    ------
!!      Christelle Barthe  * LAERO *
!!      extracted from resolved_elecn in MNH versions < 5-5-0
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/02/2022
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use mode_budget,     only: Budget_store_add
USE MODE_ELEC_ll
!
use modd_budget,     only: lbudget_sv, NBUDGET_SV1, tbudgets
USE MODD_NSV,        ONLY: NSV_ELECBEG, NSV_ELECEND ! Scalar variables for budgets
USE MODD_ELEC_DESCR, ONLY: XECHARGE
USE MODD_ELEC_n,     ONLY: XIONSOURCEFW
!
USE MODI_ION_DRIFT
USE MODI_TO_ELEC_FIELD_n
!
IMPLICIT NONE
!
!*       0.1   Declaration of dummy arguments
!
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! Temporal loop counter
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX,HLBCY   ! X and Y-direc. LBC type
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ   ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF ! Reference dry air density
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT      ! Moist variables at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT     ! Scalar variable at time t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS     ! Scalar variable sources
!
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PEFIELDU  !  Electric field
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PEFIELDV  !    components
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PEFIELDW  ! along x, y and z
!
!
!*       0.2   Declaration of local variables
!
REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) :: ZDRIFTP ! positive ion drift
REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) :: ZDRIFTN ! negative ion drift
!
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE ELECTRIC FIELD AT MASS POINTS
!               -----------------------------------------
!
PSVT(:,:,:,1)     =  XECHARGE * PSVT(:,:,:,1)    ! 1/kg --> C/kg
PSVT(:,:,:,KRR+1) = -XECHARGE * PSVT(:,:,:,KRR+1)
!
CALL TO_ELEC_FIELD_n (PRT, PSVT(:,:,:,1:KRR+1), PRHODJ, &
                      KTCOUNT, KRR,                     &
                      PEFIELDU, PEFIELDV, PEFIELDW      )
!
PSVT(:,:,:,1)     =  PSVT(:,:,:,1)     / XECHARGE    ! back to 1/kg 
PSVT(:,:,:,KRR+1) = -PSVT(:,:,:,KRR+1) / XECHARGE
!
!-------------------------------------------------------------------------------
!
!*       2.     ION SOURCE FROM DRIFT MOTION AND COSMIC RAYS
!               --------------------------------------------
!
!*       2.1    Compute source term from -/+(Div (N.mu E)) at mass points, 
!               N positive or negative ion number per kg of air (= PSVT)
!               This is a contribution of drift motion to Source PSVS for ions
!               in 1/(kg.s)
!
!CALL MYPROC_ELEC_ll (IPROC)  ! CB : utile ?
!
CALL ION_DRIFT(KRR, ZDRIFTP, ZDRIFTN, PSVT, HLBCX, HLBCY)
!
PSVS(:,:,:,1)     = PSVS(:,:,:,1)     + ZDRIFTP(:,:,:)
PSVS(:,:,:,KRR+1) = PSVS(:,:,:,KRR+1) + ZDRIFTN(:,:,:)
!
if ( lbudget_sv ) then
  call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg), 'DRIFT', zdriftp(:, :, :) * prhodj(:, :, :) )
  call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecend), 'DRIFT', zdriftn(:, :, :) * prhodj(:, :, :) )
end if
!
!
!*       2.2    Add Cosmic Ray source
!
PSVS(:,:,:,1)     = PSVS(:,:,:,1)     + XIONSOURCEFW(:,:,:) / PRHODREF(:,:,:)
PSVS(:,:,:,KRR+1) = PSVS(:,:,:,KRR+1) + XIONSOURCEFW(:,:,:) / PRHODREF(:,:,:)
!
if ( lbudget_sv ) then
  call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg), 'CORAY', xionsourcefw(:,:,:)/prhodref(:,:,:) * prhodj(:, :, :) )
  call Budget_store_add( tbudgets(NBUDGET_SV1 - 1 + nsv_elecend), 'CORAY', xionsourcefw(:,:,:)/prhodref(:,:,:) * prhodj(:, :, :) )
end if
!
END SUBROUTINE ION_SOURCE_ELEC
