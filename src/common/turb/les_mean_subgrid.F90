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
!
MODULE MODI_LES_MEAN_SUBGRID
!      #####################
!
INTERFACE LES_MEAN_SUBGRID
!
      SUBROUTINE LES_MEAN_SUBGRID_3D(PA, PA_MEAN, OSUM)

REAL,    DIMENSION(:,:,:), INTENT(IN)    :: PA
!
REAL,    DIMENSION(:,:,:), INTENT(INOUT) :: PA_MEAN
!
LOGICAL, OPTIONAL,         INTENT(IN)    :: OSUM
!
END SUBROUTINE LES_MEAN_SUBGRID_3D
!

      SUBROUTINE LES_MEAN_SUBGRID_SURF(PA, PA_MEAN, OSUM)

REAL,    DIMENSION(:,:), INTENT(IN)    :: PA
!
REAL,    DIMENSION(:),   INTENT(INOUT) :: PA_MEAN
!
LOGICAL, OPTIONAL,       INTENT(IN)    :: OSUM
!
END SUBROUTINE LES_MEAN_SUBGRID_SURF
!
END INTERFACE
!
END MODULE MODI_LES_MEAN_SUBGRID
!
!     ##############################################
      SUBROUTINE LES_MEAN_SUBGRID_3D(PA, PA_MEAN, OSUM)
!     ##############################################
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:,:), INTENT(IN)    :: PA
REAL,    DIMENSION(:,:,:), INTENT(INOUT) :: PA_MEAN
LOGICAL, OPTIONAL,         INTENT(IN)    :: OSUM
!
CALL ABORT ! AROME should not used this mesonh routine, if yes, check LLES_CALL
!
END SUBROUTINE LES_MEAN_SUBGRID_3D
!
!     ##############################################
      SUBROUTINE LES_MEAN_SUBGRID_SURF(PA, PA_MEAN, OSUM)
!     ##############################################
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:), INTENT(IN)    :: PA
REAL,    DIMENSION(:),   INTENT(INOUT) :: PA_MEAN
LOGICAL, OPTIONAL,       INTENT(IN)    :: OSUM
!
!
CALL ABORT ! AROME should not used this mesonh routine, if yes, check LLES_CALL
!
END SUBROUTINE LES_MEAN_SUBGRID_SURF
