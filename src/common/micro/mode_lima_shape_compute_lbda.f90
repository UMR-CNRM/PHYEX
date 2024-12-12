! 
MODULE MODE_LIMA_SHAPE_COMPUTE_LBDA
!
IMPLICIT NONE
CONTAINS
!
! ################################################################
  SUBROUTINE LIMA_SHAPE_COMPUTE_LBDA (KLOOP,                     &
                                      PRI, PCI_SHAPE, PRI_SHAPE, &
                                      PLBDAI_SHAPE)
! ################################################################
!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!
!!    AUTHOR
!!    ------
!!      C. Barthe  * LACy *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/04/2018
!!      C. Barthe      26/11/2019 : supprime PRIS --> inutile
!!      C. Barthe      24/01/2024 : integration into 5-7-0
!!
!-----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_PARAM_LIMA, ONLY : XRTMIN, XCTMIN, NNB_CRYSTAL_SHAPE
USE MODD_PARAM_LIMA_COLD, ONLY : XLBI_SHAPE, XLBEXI_SHAPE, XBI_SHAPE, &
                                 XFSEDRI_TOT_SHAPE
!
IMPLICIT NONE
!
!*      0.1   Declaration of dummy arguments
!
INTEGER,                 INTENT(IN)    :: KLOOP        ! dimension of the arrays
REAL,    DIMENSION(:),   INTENT(IN)    :: PRI          ! Total cloud ice m.r. (source or at t)
REAL,    DIMENSION(:,:), INTENT(IN)    :: PCI_SHAPE    ! Pristine ice conc. per shape (source or at t)
REAL,    DIMENSION(:,:), INTENT(INOUT) :: PRI_SHAPE    ! pristine ice m.r. per shape (source or at t)
REAL,    DIMENSION(:,:), INTENT(INOUT) :: PLBDAI_SHAPE ! lambda per habit
!
!
!*      0.2   Declaration of local variables
!
INTEGER :: II, JSH, &    ! loop counter
           ICOUNT_SHAPE  ! count the number of points where a shape dominates
!
INTEGER, DIMENSION(KLOOP)                   :: ISHAPE_MAX  ! Nb of the dominant habit per mesh
REAL,    DIMENSION(KLOOP)                   :: ZCI         ! Total nb conc. of ice crystals
REAL,    DIMENSION(NNB_CRYSTAL_SHAPE)       :: ZLBDAI_AVG  ! mean value of lambda per habit
REAL,    DIMENSION(KLOOP)                   :: ZLBDAI
REAL,    DIMENSION(KLOOP)                   :: ZAUX, ZSUM
REAL,    DIMENSION(KLOOP)                   :: ZONEOVER_VAR
REAL,    DIMENSION(KLOOP,NNB_CRYSTAL_SHAPE) :: ZCSHAPE  ! nb conc. ratio for each habit
REAL,    DIMENSION(KLOOP,NNB_CRYSTAL_SHAPE) :: ZRSHAPE  ! mass ratio for each habit
!
!-----------------------------------------------------------------------------
!
!*     1.    FIND THE DOMINANT SHAPE IN EACH GRID MESH
!            -----------------------------------------
!
! total ice crystal number concentration
ZCI = SUM(PCI_SHAPE,DIM=2)

ZONEOVER_VAR(:) = 0.
WHERE (ZCI(:) .GT. 0.0) ZONEOVER_VAR(:) = 1.0 / ZCI(:)

ZCSHAPE(:,:) = 0.0
DO JSH = 1, NNB_CRYSTAL_SHAPE
! compute the ratio (in terms of nb conc.): shape/tot
  WHERE ((ZCI(:) .GT. 0.0) .AND. (PCI_SHAPE(:,JSH) .GT. 0.0))
    ZCSHAPE(:,JSH) = MIN(PCI_SHAPE(:,JSH)*ZONEOVER_VAR(:), 1.0)
  END WHERE
END DO
!
!
!*      2.    COMPUTE THE MIXING RATIO PER HABIT
!             ----------------------------------
!
PLBDAI_SHAPE(:,:) = 1.E10
PRI_SHAPE(:,:) = 0
DO JSH = 1, NNB_CRYSTAL_SHAPE
  WHERE (PRI(:) > XRTMIN(4))
    PRI_SHAPE(:,JSH) = ZCSHAPE(:,JSH) * PRI(:)
  END WHERE
!
!
!*      3.    COMPUTE A MEAN LAMBDA PER HABIT
!             -------------------------------
!
!
  WHERE ((PRI_SHAPE(:,JSH) > XRTMIN(4)) .AND. (PCI_SHAPE(:,JSH)> XCTMIN(4)))
    PLBDAI_SHAPE(:,JSH) =  (XLBI_SHAPE(JSH) * PCI_SHAPE(:,JSH) / &
                            PRI_SHAPE(:,JSH))**XLBEXI_SHAPE(JSH)
  END WHERE
END DO
!
!
END SUBROUTINE LIMA_SHAPE_COMPUTE_LBDA
!
!-------------------------------------------------------------------------------
!
! ###################################################################
  SUBROUTINE LIMA_SHAPE_COMPUTE_LBDA_3D (PRI, PCI_SHAPE, PRI_SHAPE, &
                                         PLBDAI_SHAPE)
! ###################################################################
!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!
!!    AUTHOR
!!    ------
!!      C. Barthe  * LACy *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/04/2018
!!      C. Barthe      26/11/2019 : supprime PRIS --> inutile
!!      C. Barthe      24/01/2024 : integration into 5-7-0
!!
!-----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT
USE MODD_PARAM_LIMA, ONLY : XRTMIN, XCTMIN, NNB_CRYSTAL_SHAPE
USE MODD_PARAM_LIMA_COLD, ONLY : XLBI_SHAPE, XLBEXI_SHAPE, XBI_SHAPE, &
                                 XFSEDRI_TOT_SHAPE
!
IMPLICIT NONE
!
!*      0.1   Declaration of dummy arguments
!
REAL,    DIMENSION(:,:,:),   INTENT(IN)    :: PRI          ! Total cloud ice m.r. (source or at t)
REAL,    DIMENSION(:,:,:,:), INTENT(IN)    :: PCI_SHAPE    ! Pristine ice conc. per shape (source or at t)
REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PRI_SHAPE    ! pristine ice m.r. per shape (source or at t)
REAL,    DIMENSION(:,:,:,:), INTENT(INOUT) :: PLBDAI_SHAPE ! lambda per habit
!
!
!*      0.2   Declaration of local variables
!
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE   ! Physical domain
INTEGER :: II, IJ, IK, JSH, &    ! loop counter
           ICOUNT_SHAPE, &  ! count the number of points where a shape dominates
           ISHMAX
!
INTEGER, DIMENSION(SIZE(PRI,1),SIZE(PRI,2),SIZE(PRI,3)) :: ISHAPE_MAX  ! Nb of the dominant habit per mesh
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2),SIZE(PRI,3)) :: ZCI         ! Total nb conc. of ice crystals
REAL,    DIMENSION(NNB_CRYSTAL_SHAPE)                   :: ZLBDAI_AVG  ! mean value of lambda per habit
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2),SIZE(PRI,3)) :: ZLBDAI
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2),SIZE(PRI,3)) :: ZAUX, ZSUM
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2),SIZE(PRI,3)) :: ZONEOVER_VAR
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2),SIZE(PRI,3),NNB_CRYSTAL_SHAPE) :: ZCSHAPE  ! nb conc. ratio for each habit
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2),SIZE(PRI,3),NNB_CRYSTAL_SHAPE) :: ZRSHAPE  ! mass ratio for each habit
!
!-----------------------------------------------------------------------------
!
! Physical domain
IIB = 1 + JPHEXT
IIE = SIZE(PRI,1) - JPHEXT
IJB = 1 + JPHEXT
IJE = SIZE(PRI,2) - JPHEXT
IKB = 1 + JPVEXT
IKE = SIZE(PRI,3) - JPVEXT
!
!*     1.    FIND THE DOMINANT SHAPE IN EACH GRID MESH
!            -----------------------------------------
!
! total ice crystal number concentration
ZCI(:,:,:) = SUM(PCI_SHAPE,DIM=4)

ZONEOVER_VAR(:,:,:) = 0.
WHERE (ZCI(:,:,:) .GT. 0.0) ZONEOVER_VAR(:,:,:) = 1.0 / ZCI(:,:,:)

ZCSHAPE(:,:,:,:) = 0.0
DO JSH = 1, NNB_CRYSTAL_SHAPE
! compute the ratio (in terms of nb conc.): shape/tot
  WHERE ((ZCI(:,:,:) .GT. 0.0) .AND. (PCI_SHAPE(:,:,:,JSH) .GT. 0.0))
    ZCSHAPE(:,:,:,JSH) = MIN(PCI_SHAPE(:,:,:,JSH)*ZONEOVER_VAR(:,:,:), 1.0)
  END WHERE
END DO
!
PLBDAI_SHAPE(:,:,:,:) = 1.E10
PRI_SHAPE(:,:,:,:) = 0.
!
!
!*      2.    COMPUTE THE MIXING RATIO PER HABIT
!             ----------------------------------
!
PLBDAI_SHAPE(:,:,:,:) = 1.E10
PRI_SHAPE(:,:,:,:) = 0
DO JSH = 1, NNB_CRYSTAL_SHAPE

  WHERE (PRI(:,:,:) > XRTMIN(4))
    PRI_SHAPE(:,:,:,JSH) = ZCSHAPE(:,:,:,JSH) * PRI(:,:,:)
  END WHERE
!
!
!*      3.    COMPUTE A MEAN LAMBDA PER HABIT
!             -------------------------------
!
!
  WHERE (PRI_SHAPE(:,:,:,JSH) > XRTMIN(4) .AND. PCI_SHAPE(:,:,:,JSH) > XCTMIN(4))
    PLBDAI_SHAPE(:,:,:,JSH) =  (XLBI_SHAPE(JSH) * PCI_SHAPE(:,:,:,JSH) / &
                            PRI_SHAPE(:,:,:,JSH))**XLBEXI_SHAPE(JSH)
  END WHERE
END DO
!
END SUBROUTINE LIMA_SHAPE_COMPUTE_LBDA_3D
!
!-------------------------------------------------------------------------------
!
END MODULE MODE_LIMA_SHAPE_COMPUTE_LBDA
