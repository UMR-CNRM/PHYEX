!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_GRADIENT_M
!     ###################### 

! Dummy interfaces for horizontal turbulence

!
IMPLICIT NONE
CONTAINS
!
FUNCTION GX_M_M(OFLAT,PA,PDXX,PDZZ,PDZX)      RESULT(PGX_M_M)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_M_M ! result mass point
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_M', 'Prohibited call')
END FUNCTION GX_M_M
!
FUNCTION GY_M_M(OFLAT,PA,PDYY,PDZZ,PDZY)      RESULT(PGY_M_M)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_M_M ! result mass point
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_M', 'Prohibited call')
END FUNCTION GY_M_M
!
FUNCTION GZ_M_M(PA,PDZZ)      RESULT(PGZ_M_M)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_M_M ! result mass point
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_M', 'Prohibited call')
END FUNCTION GZ_M_M
!
FUNCTION GX_M_U(OFLAT,PY,PDXX,PDZZ,PDZX) RESULT(PGX_M_U)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX                   ! d*xx
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZX                   ! d*zx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   ! d*zz
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGX_M_U  ! result at flux side
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_M', 'Prohibited call')
END FUNCTION GX_M_U
!
FUNCTION GY_M_V(OFLAT,PY,PDYY,PDZZ,PDZY) RESULT(PGY_M_V)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY                   !d*yy
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZY                   !d*zy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   !d*zz
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGY_M_V  ! result at flux side
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_M', 'Prohibited call')
END FUNCTION GY_M_V
!
FUNCTION GZ_M_W(PY,PDZZ) RESULT(PGZ_M_W)
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY: NVERB_FATAL
IMPLICIT NONE
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   !d*zz
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGZ_M_W  ! result at flux side
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'MODI_GRADIENT_M', 'Prohibited call')
END FUNCTION GZ_M_W
!
END MODULE MODI_GRADIENT_M
