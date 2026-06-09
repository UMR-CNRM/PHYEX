!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_GRADIENT_M
!     ###################### 
!
IMPLICIT NONE
INTERFACE
!
!
FUNCTION GX_M_M(PA,PDXX,PDZZ,PDZX,KKA,KKU,KL)      RESULT(PGX_M_M)
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
INTEGER, INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes (AROME)
INTEGER, INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise (AROME)
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_M_M ! result mass point
!
END FUNCTION GX_M_M
!
!
FUNCTION GY_M_M(PA,PDYY,PDZZ,PDZY,KKA,KKU,KL)      RESULT(PGY_M_M)
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
INTEGER, INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes (AROME)
INTEGER, INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise (AROME)
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_M_M ! result mass point
!
END FUNCTION GY_M_M
!
!
FUNCTION GZ_M_M(PA,PDZZ,KKA,KKU,KL)      RESULT(PGZ_M_M)
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
INTEGER, INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes (AROME)
INTEGER, INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise (AROME)
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_M_M ! result mass point
!
END FUNCTION GZ_M_M
!
      FUNCTION GX_M_U(KKA,KKU,KL,PY,PDXX,PDZZ,PDZX) RESULT(PGX_M_U)
!  
IMPLICIT NONE
!
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX                   ! d*xx
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZX                   ! d*zx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   ! d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGX_M_U  ! result at flux
                                                              ! side
END FUNCTION GX_M_U
!
!
      FUNCTION GY_M_V(KKA,KKU,KL,PY,PDYY,PDZZ,PDZY) RESULT(PGY_M_V)
!
IMPLICIT NONE
!
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY                   !d*yy
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZY                   !d*zy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGY_M_V  ! result at flux
                                                              ! side
END FUNCTION GY_M_V
!
      FUNCTION GZ_M_W(KKA, KKU, KL,PY,PDZZ) RESULT(PGZ_M_W)
!  
IMPLICIT NONE
!
                                                          ! Metric coefficient:
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGZ_M_W  ! result at flux
                                                              ! side
!
END FUNCTION GZ_M_W
!
END INTERFACE
!
END MODULE MODI_GRADIENT_M
