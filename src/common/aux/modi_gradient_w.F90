!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_GRADIENT_W
!     ######################
!
! Dummy interfaces for horizontal turbulence

IMPLICIT NONE
CONTAINS
!
FUNCTION GZ_W_M(PA,PDZZ)      RESULT(PGZ_W_M)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_W_M ! result mass point
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_W', 'Prohibited call')
END FUNCTION GZ_W_M
!            
FUNCTION GX_W_UW(OFLAT,PA,PDXX,PDZZ,PDZX)      RESULT(PGX_W_UW)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_W_UW ! result UW point
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_W', 'Prohibited call')
END FUNCTION GX_W_UW
!            
FUNCTION GY_W_VW(OFLAT,PA,PDYY,PDZZ,PDZY)      RESULT(PGY_W_VW)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_W_VW ! result VW point
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_W', 'Prohibited call')
END FUNCTION GY_W_VW
!
END MODULE MODI_GRADIENT_W
