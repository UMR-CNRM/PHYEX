!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_GRADIENT_V
!     ######################
!
! Dummy interfaces for horizontal turbulence

IMPLICIT NONE
CONTAINS
!
FUNCTION GY_V_M(OFLAT,PA,PDYY,PDZZ,PDZY)      RESULT(PGY_V_M)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_V_M ! result mass point
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_V', 'Prohibited call')
END FUNCTION GY_V_M
!           
FUNCTION GX_V_UV(OFLAT,PA,PDXX,PDZZ,PDZX)      RESULT(PGX_V_UV)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_V_UV ! result UV point
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_V', 'Prohibited call')
END FUNCTION GX_V_UV
!
FUNCTION GZ_V_VW(PA,PDZZ)      RESULT(PGZ_V_VW)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_V_VW ! result VW point
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_V', 'Prohibited call')
END FUNCTION GZ_V_VW
!
END MODULE MODI_GRADIENT_V
