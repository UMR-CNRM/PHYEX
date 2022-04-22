!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##################
      MODULE MODI_SHUMAN_MF
!     ##################
!
INTERFACE
!
SUBROUTINE DZF_MF(D, PA, PDZF)
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
IMPLICIT NONE
TYPE(DIMPHYEX_t),             INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PA     ! variable at flux side
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PDZF   ! result at mass
                                                 ! localization
END SUBROUTINE DZF_MF
!
SUBROUTINE DZM_MF(D, PA, PDZM)
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
IMPLICIT NONE
TYPE(DIMPHYEX_t),             INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PDZM   ! result at flux
                                                 ! side
END SUBROUTINE DZM_MF
!
SUBROUTINE MZF_MF(D, PA, PMZF)
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
IMPLICIT NONE
TYPE(DIMPHYEX_t),             INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PA     ! variable at flux side
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PMZF   ! result at mass
                                                 ! localization
END SUBROUTINE MZF_MF
!
SUBROUTINE MZM_MF(D, PA, PMZM)
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
IMPLICIT NONE
TYPE(DIMPHYEX_t),             INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PMZM   ! result at flux localization
END SUBROUTINE MZM_MF
!
SUBROUTINE GZ_M_W_MF(D, PY, PDZZ, PGZ_M_W)
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
IMPLICIT NONE
TYPE(DIMPHYEX_t),             INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PDZZ ! Metric coefficient d*zz
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PY   ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PGZ_M_W  ! result at flux side
END SUBROUTINE GZ_M_W_MF
!
END INTERFACE
!
END MODULE MODI_SHUMAN_MF
!
!     ###############################
      SUBROUTINE MZF_MF(D, PA, PMZF)
!     ###############################
!
!!****  *MZF* -  SHUMAN_MF operator : mean operator in z direction for a
!!                                 variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean
!     along the z direction (K index) for a field PA localized at a z-flux
!     point (w point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------
!!        The result PMZF(:,:,k) is defined by 0.5*(PA(:,:,k)+PA(:,:,k+1))
!!        At k=size(PA,3), PMZF(:,:,k) is defined by PA(:,:,k).
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN_MF operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94
!!                   optimisation                 20/08/00 J. Escobar
!!      S. Riette, Jan 2012: Simplification and suppression of array overflow
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),             INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PA     ! variable at flux side
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PMZF   ! result at mass
                                                 ! localization
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK             ! Loop index in z direction
!
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZF
!              ------------------
!
DO JK=2,D%NKT-1
  PMZF(:,JK) = 0.5*( PA(:,JK)+PA(:,JK+D%NKL) )
END DO
PMZF(:,D%NKA) = 0.5*( PA(:,D%NKA)+PA(:,D%NKA+D%NKL) )
PMZF(:,D%NKU) = PA(:,D%NKU)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MZF_MF
!     ###############################
      SUBROUTINE MZM_MF(D, PA, PMZM)
!     ###############################
!
!!****  *MZM* -  SHUMAN_MF operator : mean operator in z direction for a
!!                                 mass variable
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean
!     along the z direction (K index) for a field PA localized at a mass
!     point. The result is localized at a z-flux point (w point).
!
!!**  METHOD
!!    ------
!!        The result PMZM(:,:,k) is defined by 0.5*(PA(:,:,k)+PA(:,:,k-1))
!!        At k=1, PMZM(:,:,1) is defined by PA(:,:,1).
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN_MF operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94
!!                   optimisation                 20/08/00 J. Escobar
!!      S. Riette, Jan 2012: Simplification and suppression of array overflow
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),             INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PMZM   ! result at flux localization
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK             ! Loop index in z direction
!
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZM
!              ------------------
!
DO JK=2,D%NKT-1
  PMZM(:,JK) = 0.5*( PA(:,JK)+PA(:,JK-D%NKL) )
END DO
PMZM(:,D%NKA) = PA(:,D%NKA)
PMZM(:,D%NKU) = 0.5*( PA(:,D%NKU)+PA(:,D%NKU-D%NKL) )
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MZM_MF
!     ###############################
      SUBROUTINE DZF_MF(D, PA, PDZF)
!     ###############################
!
!!****  *DZF* -  SHUMAN_MF operator : finite difference operator in z direction
!!                                  for a variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference
!     along the z direction (K index) for a field PA localized at a z-flux
!     point (w point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------
!!        The result PDZF(:,:,k) is defined by (PA(:,:,k+1)-PA(:,:,k))
!!        At k=size(PA,3), PDZF(:,:,k) is defined by 0.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN_MF operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94
!!                   optimisation                 20/08/00 J. Escobar
!!      S. Riette, Jan 2012: Simplification and suppression of array overflow
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),             INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PA     ! variable at flux side
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PDZF   ! result at mass
                                                 ! localization
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK           ! Loop index in z direction
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DZF
!              ------------------
!
DO JK=2,D%NKT-1
  PDZF(:,JK) = PA(:,JK+D%NKL) - PA(:,JK)
END DO
PDZF(:,D%NKA) = PA(:,D%NKA+D%NKL) - PA(:,D%NKA)
PDZF(:,D%NKU) = 0.
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DZF_MF
!     ###############################
      SUBROUTINE DZM_MF(D, PA, PDZM)
!     ###############################
!
!!****  *DZM* -  SHUMAN_MF operator : finite difference operator in z direction
!!                                  for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference
!     along the z direction (K index) for a field PA localized at a mass
!     point. The result is localized at a z-flux point (w point).
!
!!**  METHOD
!!    ------
!!        The result PDZM(:,j,:) is defined by (PA(:,:,k)-PA(:,:,k-1))
!!        At k=1, PDZM(:,:,k) is defined by 0.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN_MF operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94
!!                   optimisation                 20/08/00 J. Escobar
!!      S. Riette, Jan 2012: Simplification and suppression of array overflow
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),             INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PDZM   ! result at flux
                                                 ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK            ! Loop index in z direction
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DZM
!              ------------------
!
DO JK=2,D%NKT-1
  PDZM(:,JK) = PA(:,JK) - PA(:,JK-D%NKL)
END DO
PDZM(:,D%NKA) = 0.
PDZM(:,D%NKU) = PA(:,D%NKU) - PA(:,D%NKU-D%NKL)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DZM_MF

!     ###############################
      SUBROUTINE GZ_M_W_MF(D, PY, PDZZ, PGZ_M_W)
!     ###############################
!
!!****  *GZ_M_W * - Compute the gradient along z direction for a
!!       variable localized at a mass point
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!                    dzm(PY)
!       PGZ_M_W =    -------
!                     d*zz
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!    S.Riette moving of code previously in compute_mf_cloud code
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25 Aug 2011
!!      S. Riette, Jan 2012: Simplification and suppression of array overflow
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!!
!
!-------------------------------------------------------------------------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),             INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PDZZ ! Metric coefficient d*zz
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)  :: PY   ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PGZ_M_W  ! result at flux side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER  JK
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG Z
!              -----------------------------
!
DO JK=2,D%NKT-1
  PGZ_M_W(:,JK) = (PY(:,JK) - PY(:,JK-D%NKL)) / PDZZ(:,JK)
END DO
PGZ_M_W(:,D%NKA) = 0.
PGZ_M_W(:,D%NKU) = (PY(:,D%NKU) - PY(:,D%NKU-D%NKL)) / PDZZ(:,D%NKU)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GZ_M_W_MF
