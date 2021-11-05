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
FUNCTION DZF_MF(KKA,KKU,KKL,PA)  RESULT(PDZF)
INTEGER,              INTENT(IN)       :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)       :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN)       :: PA     ! variable at flux
                                                 !  side
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: PDZF   ! result at mass
                                                 ! localization 
END FUNCTION DZF_MF
!
FUNCTION DZM_MF(KKA,KKU,KKL,PA)  RESULT(PDZM)
INTEGER,              INTENT(IN)       :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)       :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN)       :: PA     ! variable at mass
                                                 ! localization
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: PDZM   ! result at flux
                                                 ! side
END FUNCTION DZM_MF
!
FUNCTION MZF_MF(KKA,KKU,KKL,PA)  RESULT(PMZF)
INTEGER,              INTENT(IN)       :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)       :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN)       :: PA     ! variable at flux
                                                 !  side
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: PMZF   ! result at mass
                                                 ! localization 
END FUNCTION MZF_MF
!
FUNCTION MZM_MF(KKA,KKU,KKL,PA)  RESULT(PMZM)
INTEGER,              INTENT(IN)       :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)       :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN)       :: PA     ! variable at mass localization
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: PMZM   ! result at flux localization 
END FUNCTION MZM_MF
!
FUNCTION GZ_M_W_MF(KKA,KKU,KKL,PY,PDZZ) RESULT(PGZ_M_W)
INTEGER,              INTENT(IN)       :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)  :: KKL  ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN)  :: PDZZ ! Metric coefficient d*zz
REAL, DIMENSION(:,:), INTENT(IN)  :: PY   ! variable at mass localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2)) :: PGZ_M_W  ! result at flux side
END FUNCTION GZ_M_W_MF
!
END INTERFACE
!
END MODULE MODI_SHUMAN_MF
!
!     ###############################
      FUNCTION MZF_MF(KKA,KKU,KKL,PA)  RESULT(PMZF)
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
!!	V. Ducrocq       * Meteo France *
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
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
INTEGER,              INTENT(IN)       :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)       :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN)       :: PA     ! variable at flux
                                                 !  side
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: PMZF   ! result at mass
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
DO JK=2,SIZE(PA,2)-1
  PMZF(:,JK) = 0.5*( PA(:,JK)+PA(:,JK+KKL) )
END DO
PMZF(:,KKA) = 0.5*( PA(:,KKA)+PA(:,KKA+KKL) )
PMZF(:,KKU) = PA(:,KKU)
!
!-------------------------------------------------------------------------------
!
END FUNCTION MZF_MF
!     ###############################
      FUNCTION MZM_MF(KKA,KKU,KKL,PA)  RESULT(PMZM)
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
!!	V. Ducrocq       * Meteo France *
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
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
INTEGER,              INTENT(IN)       :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)       :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN)       :: PA     ! variable at mass localization
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: PMZM   ! result at flux localization 
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
DO JK=2,SIZE(PA,2)-1
  PMZM(:,JK) = 0.5*( PA(:,JK)+PA(:,JK-KKL) )
END DO
PMZM(:,KKA) = PA(:,KKA)
PMZM(:,KKU) = 0.5*( PA(:,KKU)+PA(:,KKU-KKL) )
!
!-------------------------------------------------------------------------------
!
END FUNCTION MZM_MF
!     ###############################
      FUNCTION DZF_MF(KKA,KKU,KKL,PA)  RESULT(PDZF)
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
!!	V. Ducrocq       * Meteo France *
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
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
INTEGER,              INTENT(IN)       :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)       :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN)       :: PA     ! variable at flux
                                                 !  side
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: PDZF   ! result at mass
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
DO JK=2,SIZE(PA,2)-1
  PDZF(:,JK) = PA(:,JK+KKL) - PA(:,JK)
END DO
PDZF(:,KKA) = PA(:,KKA+KKL) - PA(:,KKA)
PDZF(:,KKU) = 0.
!
!-------------------------------------------------------------------------------
!
END FUNCTION DZF_MF
!     ###############################
      FUNCTION DZM_MF(KKA,KKU,KKL,PA)  RESULT(PDZM)
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
!!	V. Ducrocq       * Meteo France *
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
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
INTEGER,              INTENT(IN)       :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)       :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN)       :: PA     ! variable at mass
                                                 ! localization
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: PDZM   ! result at flux
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
DO JK=2,SIZE(PA,2)-1
  PDZM(:,JK) = PA(:,JK) - PA(:,JK-KKL)
END DO
PDZM(:,KKA) = 0.
PDZM(:,KKU) = PA(:,KKU) - PA(:,KKU-KKL)
!
!-------------------------------------------------------------------------------
!
END FUNCTION DZM_MF

!     ###############################
      FUNCTION GZ_M_W_MF(KKA,KKU,KKL,PY,PDZZ) RESULT(PGZ_M_W)
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
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
INTEGER,              INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)  :: KKL  ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:), INTENT(IN)  :: PDZZ ! Metric coefficient d*zz
REAL, DIMENSION(:,:), INTENT(IN)  :: PY   ! variable at mass localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2)) :: PGZ_M_W  ! result at flux side
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
DO JK=2,SIZE(PY,2)-1
  PGZ_M_W(:,JK) = (PY(:,JK) - PY(:,JK-KKL)) / PDZZ(:,JK)
END DO
PGZ_M_W(:,KKA) = 0.
PGZ_M_W(:,KKU) = (PY(:,KKU) - PY(:,KKU-KKL)) / PDZZ(:,KKU)
!
!-------------------------------------------------------------------------------
!
END FUNCTION GZ_M_W_MF
