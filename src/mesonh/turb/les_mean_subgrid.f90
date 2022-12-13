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
!
!!****  *LES_MEAN_SUBGRID* computes the average of one subgrid
!!                         field on one processor
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/02/00
!!      V. Masson        06/11/02 use of 2D masks
!!       C.Lac           10/2014 : Correction on user masks
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_LES
!
USE MODI_LES_VER_INT
USE MODI_LES_MEAN_1PROC
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:,:), INTENT(IN)    :: PA
!
REAL,    DIMENSION(:,:,:), INTENT(INOUT) :: PA_MEAN
!
LOGICAL, OPTIONAL,         INTENT(IN)    :: OSUM
!
!
!       0.2  declaration of local variables
!
REAL,    DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA_MEAN,1)) :: ZA_LES
REAL,    DIMENSION(SIZE(PA_MEAN,1))                       :: ZA_MEAN
REAL,    DIMENSION(SIZE(PA_MEAN,1))                       :: ZA_MEAN_OLD
LOGICAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA_MEAN,1)) :: GMASK
!
INTEGER, DIMENSION(SIZE(PA_MEAN,1)) :: IAVG_PTS
INTEGER, DIMENSION(SIZE(PA_MEAN,1)) :: IUND_PTS
!
INTEGER                             :: IMASK     ! mask counter
INTEGER                             :: JI        ! loop control
!-------------------------------------------------------------------------------
!
IF (.NOT. LLES_CALL) RETURN
!
ZA_MEAN_OLD(:) = 0.
!-------------------------------------------------------------------------------
!
!* interpolation on LES vertical levels.
!
CALL LES_VER_INT(PA,ZA_LES)
!
!* subgrid computations on cartesian mask
!  --------------------------------------
!
IMASK = 1
!
!* averaging on the current processor domain of the subgrid variable
!
CALL LES_MEAN_1PROC(ZA_LES, LLES_CURRENT_CART_MASK(:,:,:), ZA_MEAN, IAVG_PTS, IUND_PTS)
!
IF (PRESENT(OSUM)) THEN
  IF (OSUM) ZA_MEAN_OLD(:) = PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK)
END IF
PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK) = ZA_MEAN_OLD(:) + ZA_MEAN(:)
!
!-------------------------------------------------------------------------------
!
!* subgrid computations on nebulosity/clear-sky masks
!  --------------------------------------------------
!
IF (LLES_NEB_MASK) THEN
!
!* on nebulosity mask
!  ------------------
!
  IMASK = IMASK + 1
!
  GMASK(:,:,:) = LLES_CURRENT_NEB_MASK (:,:,:) .AND.  LLES_CURRENT_CART_MASK(:,:,:)
!
!* averaging on the current processor domain of the subgrid variable
!
  CALL LES_MEAN_1PROC(ZA_LES, GMASK, ZA_MEAN, IAVG_PTS, IUND_PTS)
!
  IF (PRESENT(OSUM)) THEN
    IF (OSUM) ZA_MEAN_OLD(:) = PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK)
  END IF
  PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK) = ZA_MEAN_OLD(:) + ZA_MEAN(:)
!
!* on clear-sky mask
!  -----------------
!
  IMASK = IMASK + 1
!
  GMASK(:,:,:) = (.NOT. LLES_CURRENT_NEB_MASK (:,:,:)) .AND. LLES_CURRENT_CART_MASK(:,:,:)
!
!* averaging on the current processor domain of the subgrid variable
!
  CALL LES_MEAN_1PROC(ZA_LES, GMASK, ZA_MEAN, IAVG_PTS, IUND_PTS)
!
  IF (PRESENT(OSUM)) THEN
    IF (OSUM) ZA_MEAN_OLD(:) = PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK)
  END IF
  PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK) = ZA_MEAN_OLD(:) + ZA_MEAN(:)
!
END IF
!
!-------------------------------------------------------------------------------
!
!* subgrid computations on core/no core masks
!  --------------------------------------------------------------
!
IF (LLES_CORE_MASK) THEN
!
!* on core mask
!  ------------
!
  IMASK = IMASK + 1
!
  GMASK(:,:,:) = LLES_CURRENT_CORE_MASK(:,:,:) .AND. LLES_CURRENT_CART_MASK(:,:,:)
!
!* averaging on the current processor domain of the subgrid variable
!
  CALL LES_MEAN_1PROC(ZA_LES, GMASK, ZA_MEAN, IAVG_PTS, IUND_PTS)
!
  IF (PRESENT(OSUM)) THEN
    IF (OSUM) ZA_MEAN_OLD(:) = PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK)
  END IF
  PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK) = ZA_MEAN_OLD(:) + ZA_MEAN(:)
!
!* on NO core mask
!  ------------------------
!
  IMASK = IMASK + 1
!
  GMASK(:,:,:) = (.NOT. LLES_CURRENT_CORE_MASK(:,:,:)) .AND. LLES_CURRENT_CART_MASK(:,:,:)
!
!* averaging on the current processor domain of the subgrid variable
!
  CALL LES_MEAN_1PROC(ZA_LES, GMASK, ZA_MEAN, IAVG_PTS, IUND_PTS)
!
  IF (PRESENT(OSUM)) THEN
    IF (OSUM) ZA_MEAN_OLD(:) = PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK)
  END IF
  PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK) = ZA_MEAN_OLD(:) + ZA_MEAN(:)
END IF
!
!-------------------------------------------------------------------------------
!
!* subgrid computations on conditional sampling mask
!  -------------------------------------------------
!
IF (LLES_CS_MASK) THEN
  IMASK = IMASK + 1
!
  GMASK(:,:,:) = LLES_CURRENT_CS1_MASK(:,:,:)
!
!* averaging on the current processor domain of the subgrid variable
!
  CALL LES_MEAN_1PROC(ZA_LES, GMASK, ZA_MEAN, IAVG_PTS, IUND_PTS)
!
  IF (PRESENT(OSUM)) THEN
    IF (OSUM) ZA_MEAN_OLD(:) = PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK)
  END IF
  PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK) = ZA_MEAN_OLD(:) + ZA_MEAN(:)
!
  IMASK = IMASK + 1
!
  GMASK(:,:,:) = LLES_CURRENT_CS2_MASK(:,:,:)
!
!* averaging on the current processor domain of the subgrid variable
!
  CALL LES_MEAN_1PROC(ZA_LES, GMASK, ZA_MEAN, IAVG_PTS, IUND_PTS)
!
  IF (PRESENT(OSUM)) THEN
    IF (OSUM) ZA_MEAN_OLD(:) = PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK)
  END IF
  PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK) = ZA_MEAN_OLD(:) + ZA_MEAN(:)
!
  IMASK = IMASK + 1
!
  GMASK(:,:,:) = LLES_CURRENT_CS3_MASK(:,:,:)
!
!* averaging on the current processor domain of the subgrid variable
!
  CALL LES_MEAN_1PROC(ZA_LES, GMASK, ZA_MEAN, IAVG_PTS, IUND_PTS)
!
  IF (PRESENT(OSUM)) THEN
    IF (OSUM) ZA_MEAN_OLD(:) = PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK)
  END IF
  PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK) = ZA_MEAN_OLD(:) + ZA_MEAN(:)
END IF
!
!-------------------------------------------------------------------------------
!
!* subgrid computations on user mask
!  ---------------------------------
!
IF (LLES_MY_MASK) THEN
 DO JI=1,NLES_MASKS_USER
  IMASK = IMASK + 1
!
  GMASK(:,:,:) = LLES_CURRENT_MY_MASKS(:,:,:,JI)
!
!* averaging on the current processor domain of the subgrid variable
!
  CALL LES_MEAN_1PROC(ZA_LES, GMASK, ZA_MEAN, IAVG_PTS, IUND_PTS)
!
  IF (PRESENT(OSUM)) THEN
    IF (OSUM) ZA_MEAN_OLD(:) = PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK)
  END IF
  PA_MEAN(:,NLES_CURRENT_TCOUNT,IMASK) = ZA_MEAN_OLD(:) + ZA_MEAN(:)
 END DO
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_MEAN_SUBGRID_3D
!
!     ##############################################
      SUBROUTINE LES_MEAN_SUBGRID_SURF(PA, PA_MEAN, OSUM)
!     ##############################################
!
!
!!****  *LES_MEAN_SUBGRID* computes the average of one subgrid
!!                         field on one processor
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/02/00
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_LES
!
USE MODI_LES_MEAN_1PROC
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:), INTENT(IN)    :: PA
!
REAL,    DIMENSION(:),   INTENT(INOUT) :: PA_MEAN
!
LOGICAL, OPTIONAL,       INTENT(IN)    :: OSUM
!
!
!       0.2  declaration of local variables
!
REAL    :: ZA_MEAN
REAL    :: ZA_MEAN_OLD
!
INTEGER :: IAVG_PTS
INTEGER :: IUND_PTS
!
!-------------------------------------------------------------------------------
!
IF (.NOT. LLES_CALL) RETURN
!
ZA_MEAN_OLD = 0.
IF (PRESENT(OSUM)) THEN
  IF (OSUM) ZA_MEAN_OLD = PA_MEAN(NLES_CURRENT_TCOUNT)
END IF
!-------------------------------------------------------------------------------
!
!* subgrid computations on cartesian mask
!  --------------------------------------
!
!* averaging on the current processor domain of the subgrid variable
!
CALL LES_MEAN_1PROC(PA, LLES_CURRENT_CART_MASK(:,:,1), ZA_MEAN, IAVG_PTS, IUND_PTS)
!
PA_MEAN(NLES_CURRENT_TCOUNT) = ZA_MEAN_OLD + ZA_MEAN
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_MEAN_SUBGRID_SURF
