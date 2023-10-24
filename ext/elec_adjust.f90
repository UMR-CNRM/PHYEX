!MNH_LIC Copyright 2022-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_ELEC_ADJUST
!     #######################
!
INTERFACE
!
      SUBROUTINE ELEC_ADJUST (KRR, PRHODJ, HCLOUD, HBUNAME, &
                              PRC, PRI, PQC, PQI,           &
                              PQCS, PQIS, PQPIS, PQNIS, PCND, PDEP)
!
INTEGER,                INTENT(IN)    :: KRR      ! Number of moist variables
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRHODJ
CHARACTER(LEN=4),       INTENT(IN)    :: HCLOUD   ! Kind of microphysical scheme
CHARACTER(len=*),       INTENT(IN)    :: HBUNAME  ! Name of the budget
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRC      ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRI      ! Cloud ice  m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PQC      ! Cloud water charge density to adjust
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PQI      ! Cloud ice  charge density to adjust
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQCS     ! Cloud water charge density source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQIS     ! Cloud ice  charge density source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQPIS    ! Positive ion charge density source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQNIS    ! Negative ion charge density source
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PCND     ! Rate of condensation
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PDEP     ! Rate of sublimation
!
END SUBROUTINE ELEC_ADJUST
!
END INTERFACE
!
END MODULE MODI_ELEC_ADJUST

!     #############################################################
      SUBROUTINE ELEC_ADJUST (KRR, PRHODJ, HCLOUD, HBUNAME, &
                              PRC, PRI, PQC, PQI,           &
                              PQCS, PQIS, PQPIS, PQNIS, PCND, PDEP)
!     #############################################################
!
!!****  *ELEC_ADJUST* -  compute the exchange of electric charges associated with
!!                       condensation and sublimation of ice crystals.
!!                       The capture of ions by cloud droplets and ice crystals is done 
!!                       in the ion_attach_elec routine.
!!
!!    PURPOSE
!!    -------
!!    The purpose of this routine is to compute the fast microphysical sources
!!    through a saturation ajustement procedure in case of mixed-phase clouds.
!!
!!
!!**  METHOD
!!    ------
!!    The rate of charge exchanged is computed proportionnaly to the rate of mass
!!    exchanged during sublimation and condensation. The sublimation and condensation rates
!!    are computed in ICE3/4 or LIMA.
!!     
!!
!!    AUTHOR
!!    ------
!!      C. Barthe    * Laboratoire d'Aerologie*
!!   
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/09/2022 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,          only: lbudget_sv, NBUDGET_SV1, tbudgets
!
USE MODD_ELEC_DESCR,      ONLY: XRTMIN_ELEC, XQTMIN, XFC, XFI, XECHARGE
USE MODD_PARAM_LIMA_COLD, ONLY: XBI_L=>XBI
USE MODD_RAIN_ICE_DESCR_n,ONLY: XRTMIN_I=>XRTMIN, XBI_I=>XBI
USE MODD_PARAM_LIMA,      ONLY: XRTMIN_L=>XRTMIN
USE MODD_NSV,             ONLY: NSV_ELECBEG
!
use mode_budget,       only: Budget_store_init, Budget_store_end
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments
!
INTEGER,                INTENT(IN)    :: KRR      ! Number of moist variables
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRHODJ
CHARACTER(LEN=4),       INTENT(IN)    :: HCLOUD   ! Kind of microphysical scheme
CHARACTER(len=*),       INTENT(IN)    :: HBUNAME  ! Name of the budget
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRC      ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PRI      ! Cloud ice  m.r. at t
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PQC      ! Cloud water charge density to adjust
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PQI      ! Cloud ice  charge density to adjust
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQCS     ! Cloud water charge density source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQIS     ! Cloud ice  charge density source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQPIS    ! Positive ion charge density source
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PQNIS    ! Negative ion charge density source
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PCND     ! Rate of condensation
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PDEP     ! Rate of sublimation
!
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) :: ZWELEC      ! Work array
REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) :: ZION_NUMBER ! Nb of elementary charge in hydrometeor charge
REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) :: ZADD        ! Ratio (0/1) of ZION_NUMBER to add to positive 
                                                                             ! or negative ion nb
REAL :: ZBI
REAL, DIMENSION(KRR) :: ZRTMIN
!
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
! choose the right parameters between ice3 and lima
IF (HCLOUD(1:3) == 'ICE') THEN
  ZRTMIN(1:KRR) = XRTMIN_I(1:KRR)
  ZBI = XBI_I
ELSE IF (HCLOUD == 'LIMA') THEN
  ZRTMIN(1:KRR) = XRTMIN_L(1:KRR)
  ZBI = XBI_L
END IF
!
if ( lbudget_sv ) then
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg ),      trim( hbuname ), pqpis(:, :, :) * prhodj(:, :, :) )
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1),   trim( hbuname ), pqcs(:, :, :)  * prhodj(:, :, :) )
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3),   trim( hbuname ), pqis(:, :, :)  * prhodj(:, :, :) )
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + krr), trim( hbuname ), pqnis(:, :, :) * prhodj(:, :, :) )
end if
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE SOURCES FOR ELECTRIC CHARGES
!               ----------------------------------------
!
ZWELEC(:,:,:) = 0.
!
!*       2.1    Evaporation of cloud droplets
!
WHERE (ABS(PRC(:,:,:)) > XRTMIN_ELEC(2) .AND. &
       ABS(PQC(:,:,:)) > XQTMIN(2)      .AND. &
           PCND(:,:,:) < -ZRTMIN(1))
  ZWELEC(:,:,:) = (XFC / 3.) * (PQC(:,:,:) / PRC(:,:,:)) * (-PCND(:,:,:))
  ! nb of elementary charges in hydrometeor charge
  ZION_NUMBER(:,:,:) = ABS(ZWELEC(:,:,:)) / XECHARGE
  ! ratio (0 or 1) of the number of ions to add to positive or negative ion number
  ZADD(:,:,:) = 0.5 + SIGN(0.5, ZWELEC(:,:,:))
  !
  PQPIS(:,:,:) = PQPIS(:,:,:) +       ZADD(:,:,:)  * ZION_NUMBER(:,:,:)
  PQNIS(:,:,:) = PQNIS(:,:,:) + (1. - ZADD(:,:,:)) * ZION_NUMBER(:,:,:)
  PQCS(:,:,:)  = PQCS(:,:,:)  - ZWELEC(:,:,:)
END WHERE
!
!
!*       2.2    Sublimation of ice crystals
!
WHERE (ABS(PRI(:,:,:)) > XRTMIN_ELEC(4) .AND. &
       ABS(PQI(:,:,:)) > XQTMIN(4)      .AND. &
           PDEP(:,:,:) < -ZRTMIN(1))
  ZWELEC(:,:,:) = (XFI / ZBI) * (PQI(:,:,:) / PRI(:,:,:)) * (-PDEP(:,:,:))
  ZION_NUMBER(:,:,:) = ABS(ZWELEC(:,:,:)) / XECHARGE
  ZADD(:,:,:) = 0.5 + SIGN(0.5, ZWELEC(:,:,:))
  !
  PQPIS(:,:,:) = PQPIS(:,:,:) +       ZADD(:,:,:)  * ZION_NUMBER(:,:,:)
  PQNIS(:,:,:) = PQNIS(:,:,:) + (1. - ZADD(:,:,:)) * ZION_NUMBER(:,:,:)
  PQIS(:,:,:)  = PQIS(:,:,:)  - ZWELEC(:,:,:)
END WHERE
!
!-------------------------------------------------------------------------------
!
!*       3.     STORE THE BUDGET TERMS
!               ----------------------

if ( lbudget_sv ) then
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg ),      trim( hbuname ), pqpis(:, :, :) * prhodj(:, :, :) )
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 1),   trim( hbuname ), pqcs(:, :, :)  * prhodj(:, :, :) )
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + 3),   trim( hbuname ), pqis(:, :, :)  * prhodj(:, :, :) )
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_elecbeg + krr), trim( hbuname ), pqnis(:, :, :) * prhodj(:, :, :) )
end if
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ELEC_ADJUST
