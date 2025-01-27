! 
MODULE MODE_LIMA_SHAPE_COMPUTE_LBDA
!
IMPLICIT NONE
CONTAINS
!
! ################################################################
  SUBROUTINE LIMA_SHAPE_COMPUTE_LBDA (LIMAP, LIMAC, KLOOP,       &
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
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*      0.1   Declaration of dummy arguments
!
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
INTEGER,                 INTENT(IN)    :: KLOOP        ! dimension of the arrays
REAL,    DIMENSION(KLOOP),   INTENT(IN)    :: PRI          ! Total cloud ice m.r. (source or at t)
REAL,    DIMENSION(KLOOP,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN)    :: PCI_SHAPE    ! Pristine ice conc. per shape (source or at t)
REAL,    DIMENSION(KLOOP,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(INOUT) :: PRI_SHAPE    ! pristine ice m.r. per shape (source or at t)
REAL,    DIMENSION(KLOOP,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(INOUT) :: PLBDAI_SHAPE ! lambda per habit
!
!
!*      0.2   Declaration of local variables
!
INTEGER :: II, ISH, &    ! loop counter
           ICOUNT_SHAPE  ! count the number of points where a shape dominates
!
INTEGER, DIMENSION(KLOOP)                   :: ISHAPE_MAX  ! Nb of the dominant habit per mesh
REAL,    DIMENSION(KLOOP)                   :: ZCI         ! Total nb conc. of ice crystals
REAL,    DIMENSION(LIMAP%NNB_CRYSTAL_SHAPE)       :: ZLBDAI_AVG  ! mean value of lambda per habit
REAL,    DIMENSION(KLOOP)                   :: ZLBDAI
REAL,    DIMENSION(KLOOP)                   :: ZAUX, ZSUM
REAL,    DIMENSION(KLOOP)                   :: ZONEOVER_VAR
REAL,    DIMENSION(KLOOP,LIMAP%NNB_CRYSTAL_SHAPE) :: ZCSHAPE  ! nb conc. ratio for each habit
REAL,    DIMENSION(KLOOP,LIMAP%NNB_CRYSTAL_SHAPE) :: ZRSHAPE  ! mass ratio for each habit
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------------------
!
!*     1.    FIND THE DOMINANT SHAPE IN EACH GRID MESH
!            -----------------------------------------
!
! total ice crystal number concentration
IF (LHOOK) CALL DR_HOOK('LIMA_SHAPE_COMPUTE_LBDA', 0, ZHOOK_HANDLE)
ZCI = SUM(PCI_SHAPE,DIM=2)

ZONEOVER_VAR(:) = 0.
WHERE (ZCI(:) .GT. 0.0) ZONEOVER_VAR(:) = 1.0 / ZCI(:)

ZCSHAPE(:,:) = 0.0
DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
! compute the ratio (in terms of nb conc.): shape/tot
  WHERE ((ZCI(:) .GT. 0.0) .AND. (PCI_SHAPE(:,ISH) .GT. 0.0))
    ZCSHAPE(:,ISH) = MIN(PCI_SHAPE(:,ISH)*ZONEOVER_VAR(:), 1.0)
  END WHERE
END DO
!
!
!*      2.    COMPUTE THE MIXING RATIO PER HABIT
!             ----------------------------------
!
PLBDAI_SHAPE(:,:) = 1.E10
PRI_SHAPE(:,:) = 0
DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
  WHERE (PRI(:) > LIMAP%XRTMIN(4))
    PRI_SHAPE(:,ISH) = ZCSHAPE(:,ISH) * PRI(:)
  END WHERE
!
!
!*      3.    COMPUTE A MEAN LAMBDA PER HABIT
!             -------------------------------
!
!
  WHERE ((PRI_SHAPE(:,ISH) > LIMAP%XRTMIN(4)) .AND. (PCI_SHAPE(:,ISH)> LIMAP%XCTMIN(4)))
    PLBDAI_SHAPE(:,ISH) =  (LIMAC%XLBI_SHAPE(ISH) * PCI_SHAPE(:,ISH) / &
                            PRI_SHAPE(:,ISH))**LIMAC%XLBEXI_SHAPE(ISH)
  END WHERE
END DO
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_SHAPE_COMPUTE_LBDA', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_SHAPE_COMPUTE_LBDA
!
!-------------------------------------------------------------------------------
!
! ###################################################################
  SUBROUTINE LIMA_SHAPE_COMPUTE_LBDA_3D (D, LIMAP, LIMAC, PRI, PCI_SHAPE, PRI_SHAPE, &
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
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*      0.1   Declaration of dummy arguments
!
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
!
REAL,    DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRI          ! Total cloud ice m.r. (source or at t)
REAL,    DIMENSION(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN)    :: PCI_SHAPE    ! Pristine ice conc. per shape (source or at t)
REAL,    DIMENSION(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(INOUT) :: PRI_SHAPE    ! pristine ice m.r. per shape (source or at t)
REAL,    DIMENSION(D%NIJT,D%NKT,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(INOUT) :: PLBDAI_SHAPE ! lambda per habit
!
!
!*      0.2   Declaration of local variables
!
INTEGER :: II, IJ, IK, ISH, &    ! loop counter
           ICOUNT_SHAPE, &  ! count the number of points where a shape dominates
           ISHMAX
!
INTEGER, DIMENSION(SIZE(PRI,1),SIZE(PRI,2)) :: ISHAPE_MAX  ! Nb of the dominant habit per mesh
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2)) :: ZCI         ! Total nb conc. of ice crystals
REAL,    DIMENSION(LIMAP%NNB_CRYSTAL_SHAPE) :: ZLBDAI_AVG  ! mean value of lambda per habit
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2)) :: ZLBDAI
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2)) :: ZAUX, ZSUM
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2)) :: ZONEOVER_VAR
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2),LIMAP%NNB_CRYSTAL_SHAPE) :: ZCSHAPE  ! nb conc. ratio for each habit
REAL,    DIMENSION(SIZE(PRI,1),SIZE(PRI,2),LIMAP%NNB_CRYSTAL_SHAPE) :: ZRSHAPE  ! mass ratio for each habit
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------------------
!
!
!*     1.    FIND THE DOMINANT SHAPE IN EACH GRID MESH
!            -----------------------------------------
!
! total ice crystal number concentration
IF (LHOOK) CALL DR_HOOK('LIMA_SHAPE_COMPUTE_LBDA_3D', 0, ZHOOK_HANDLE)
ZCI(:,:) = SUM(PCI_SHAPE,DIM=3)

ZONEOVER_VAR(:,:) = 0.
WHERE (ZCI(:,:) .GT. 0.0) ZONEOVER_VAR(:,:) = 1.0 / ZCI(:,:)

ZCSHAPE(:,:,:) = 0.0
DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
! compute the ratio (in terms of nb conc.): shape/tot
  WHERE ((ZCI(:,:) .GT. 0.0) .AND. (PCI_SHAPE(:,:,ISH) .GT. 0.0))
    ZCSHAPE(:,:,ISH) = MIN(PCI_SHAPE(:,:,ISH)*ZONEOVER_VAR(:,:), 1.0)
  END WHERE
END DO
!
PLBDAI_SHAPE(:,:,:) = 1.E10
PRI_SHAPE(:,:,:) = 0.
!
!
!*      2.    COMPUTE THE MIXING RATIO PER HABIT
!             ----------------------------------
!
PLBDAI_SHAPE(:,:,:) = 1.E10
PRI_SHAPE(:,:,:) = 0
DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE

  WHERE (PRI(:,:) > LIMAP%XRTMIN(4))
    PRI_SHAPE(:,:,ISH) = ZCSHAPE(:,:,ISH) * PRI(:,:)
  END WHERE
!
!
!*      3.    COMPUTE A MEAN LAMBDA PER HABIT
!             -------------------------------
!
!
  WHERE (PRI_SHAPE(:,:,ISH) > LIMAP%XRTMIN(4) .AND. PCI_SHAPE(:,:,ISH) > LIMAP%XCTMIN(4))
    PLBDAI_SHAPE(:,:,ISH) =  (LIMAC%XLBI_SHAPE(ISH) * PCI_SHAPE(:,:,ISH) / &
                            PRI_SHAPE(:,:,ISH))**LIMAC%XLBEXI_SHAPE(ISH)
  END WHERE
END DO
!
IF (LHOOK) CALL DR_HOOK('LIMA_SHAPE_COMPUTE_LBDA_3D', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_SHAPE_COMPUTE_LBDA_3D
!
!-------------------------------------------------------------------------------
!
END MODULE MODE_LIMA_SHAPE_COMPUTE_LBDA
