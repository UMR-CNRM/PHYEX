!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 les 2006/08/30 18:41:10
!-----------------------------------------------------------------
!      #####################
MODULE MODI_LES_MEAN_SUBGRID_PHY
!      #####################
!
INTERFACE LES_MEAN_SUBGRID_PHY
!

SUBROUTINE LES_MEAN_SUBGRID_3D_PHY(D,PA, PA_MEAN, OSUM)
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_LES, ONLY: NLES_CURRENT_TCOUNT
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),          INTENT(IN)    :: D
REAL,    DIMENSION(:,:), INTENT(IN)    :: PA
REAL,    DIMENSION(:,:,:), INTENT(INOUT) :: PA_MEAN
LOGICAL, OPTIONAL,         INTENT(IN)    :: OSUM
!
END SUBROUTINE LES_MEAN_SUBGRID_3D_PHY
!
SUBROUTINE LES_MEAN_SUBGRID_SURF_PHY(D,PA, PA_MEAN, OSUM)
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),          INTENT(IN)    :: D
REAL,    DIMENSION(:), INTENT(IN)    :: PA
REAL,    DIMENSION(:), INTENT(INOUT) :: PA_MEAN
LOGICAL, OPTIONAL,         INTENT(IN)    :: OSUM
!
END SUBROUTINE LES_MEAN_SUBGRID_SURF_PHY
!
END INTERFACE LES_MEAN_SUBGRID_PHY
!
END MODULE MODI_LES_MEAN_SUBGRID_PHY
!
!     ##############################################
      SUBROUTINE LES_MEAN_SUBGRID_3D_PHY(D, PA, PA_MEAN, OSUM)
!     ##############################################
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),          INTENT(IN)    :: D
REAL,    DIMENSION(:,:), INTENT(IN)    :: PA
REAL,    DIMENSION(:,:,:), INTENT(INOUT) :: PA_MEAN
LOGICAL, OPTIONAL,         INTENT(IN)    :: OSUM
!
END SUBROUTINE LES_MEAN_SUBGRID_3D_PHY
!
!     ##############################################
      SUBROUTINE LES_MEAN_SUBGRID_SURF_PHY(D, PA, PA_MEAN, OSUM)
!     ##############################################
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),          INTENT(IN)    :: D
REAL,    DIMENSION(:), INTENT(IN)    :: PA
REAL,    DIMENSION(:), INTENT(INOUT) :: PA_MEAN
LOGICAL, OPTIONAL,         INTENT(IN)    :: OSUM
!
END SUBROUTINE LES_MEAN_SUBGRID_SURF_PHY
