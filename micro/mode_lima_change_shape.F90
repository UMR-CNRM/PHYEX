MODULE MODE_LIMA_CHANGE_SHAPE
!
  IMPLICIT NONE
!
CONTAINS
!
! ######################################################################
  SUBROUTINE LIMA_CHANGE_SHAPE(CST, LIMAP, KSIZE, &
                               PT, PCTMIN, PDEP, PCI_SHAPE, P_SHCI_HACH)
! ######################################################################
!
!!    PURPOSE
!!    -------
!!    The purpose of this routine is to compute the change of shape of small ice
!!    crystals, when plates are in a temperature range of column growth and vice
!!    versa. Plates and columns turn into irregulars definitely.
!! 
!!    QUESTIONS : -  Rajouter une valeur seuil sur la quantité de vapeur d'eau
!!                déposée sur les cristaux (PDEP en entrée de la routine pour
!!                l'instant mais qui ne sert pas pour l'instant)?
!!                - 
!!
!!    METHOD
!!    ------
!!
!!    AUTHOR
!!    ------
!!      M. Claeys  * LACy *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       07/11/2018
!!      C. Barthe  02/2024: adapt to MNH-V5-7-0 + clean the code
!!      C. Barthe  04/2024: irregular no more used, bullet-rosettes instead (shape n°3)
!
!-----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK

USE MODD_CST, ONLY : XTT
!
IMPLICIT NONE
!
!
!*      0.1   Declaration of dummy arguments
!
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(CST_T),INTENT(IN)::CST
INTEGER,              INTENT(IN)    :: KSIZE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT  ! Temp at t (K)
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PDEP  ! water vapor deposited on ice crystals
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCTMIN
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN)    :: PCI_SHAPE  ! Nb conc. of ice crystal shapes
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(INOUT) :: P_SHCI_HACH ! for budget
!
!
!*      0.2   Declaration of local variables
!
REAL, DIMENSION(SIZE(PT))    :: ZTC ! temperature (C)
REAL, DIMENSION(SIZE(PT))    :: ZFACTOR_PLA_IRR !  
REAL, DIMENSION(SIZE(PT))    :: ZFACTOR_COL_IRR !
!
INTEGER, DIMENSION(SIZE(PT)) :: IDX_CHEN
INTEGER                      :: II  ! index
!
CHARACTER(LEN=3), DIMENSION(SIZE(PT)) :: YSHAPE_NAME ! Name of the shape corresponding to 
                                                     ! the temperature growth regime, at t
!
! Chen 90
REAL, DIMENSION(121)      :: ZTEMP_CHEN ! Lookup table T (Chen and Lamb, 1994)
REAL, DIMENSION(SIZE(PT)) :: ZGAM_CHEN_T
REAL, DIMENSION(121)      :: ZGAM_CHEN ! Lookup table gamma(T) (Chen and Lamb, 1994)
                                       ! Index_GAM_CHEN = MAX( 1,MIN( 121,INT((XTT-T)*4.0)+1 ) )
                                       ! XTT = 273.16 and T is the temperature
!
REAL, DIMENSION(SIZE(PT)) :: ZDGAMMA_COL_IRR !
REAL, DIMENSION(SIZE(PT)) :: ZDGAMMA_PLA_IRR !
!
REAL, DIMENSION(SIZE(PT)) :: ZZW1, ZZW2
!
REAL                      :: ZGAMMA_COL, ZGAMMA_PLA, ZGAMMA_IRR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
! 
!-----------------------------------------------------------------------------
!
!*       0.   PRELIMINARIES
!             -------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_CHANGE_SHAPE', 0, ZHOOK_HANDLE)
ZZW1(:) = 0.
ZZW2(:) = 0.
!
! Change temperature unit: K --> C
ZTC(:) = PT(:) - CST%XTT
!
! Find the ice crystal shape associated with temperature in each grid mesh
YSHAPE_NAME(:) = 'IRR'
WHERE ( (ZTC(:) .GT. -3.0) .OR.    &
       ((ZTC(:) .LE. -9.0) .AND. (ZTC(:) .GT. -20.0)))
  YSHAPE_NAME(:) = 'PLA'
END WHERE
!
WHERE (((ZTC(:) .LE. -3.0) .AND. (ZTC(:) .GT. -9.0)) .OR. &
        (ZTC(:) .LE. -22))
  YSHAPE_NAME(:) = 'COL'
END WHERE
!
!
!*       1.   COMPUTE GAMMA(T) 
!             ----------------
!
! (marine : à mettre dans modd_lima_cold_mixed ?)
!
!*       1.1  Table for Gamma and temperature (Chen and Lamb, 1994)
!
ZTEMP_CHEN(:) = (/ &
& 0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, &
& 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.25, &  
& 5.50, 5.75, 6.00, 6.25, 6.50, 6.75, 7.00, 7.25, 7.50, 7.75, 8.00, &
& 8.25, 8.50, 8.75, 9.00, 9.25, 9.50, 9.75, 10.00, 10.25, 10.50, 10.75, & 
& 11.00, 11.25, 11.50, 11.75, 12.00, 12.25, 12.50, 12.75, 13.00, 13.25, 13.50, &
& 13.75, 14.00, 14.25, 14.50, 14.75, 15.00, 15.25, 15.50, 15.75, 16.00, 16.25, &  
& 16.50, 16.75, 17.00, 17.25, 17.50, 17.75, 18.00, 18.25, 18.50, 18.75, 19.00, &
& 19.25, 19.50, 19.75, 20.00, 20.25, 20.50, 20.75, 21.00, 21.25, 21.50, 21.75, & 
& 22.00, 22.25, 22.50, 22.75, 23.00, 23.25, 23.50, 23.75, 24.00, 24.25, 24.50, &
& 24.75, 25.00, 25.25, 25.50, 25.75, 26.00, 26.25, 26.50, 26.75, 27.00, 27.25, &  
& 27.50, 27.75, 28.00, 28.25, 28.50, 28.75, 29.00, 29.25, 29.50, 29.75, 30.00 /)
!
ZTEMP_CHEN(:) = -ZTEMP_CHEN(:)
!
ZGAM_CHEN(1:121) = (/ &
& 1.00, 0.98, 0.96, 0.94, 0.92, 0.90, 0.88, 0.86, 0.83, 0.81, 0.78, &
& 0.76, 0.70, 0.54, 0.47, 0.52, 0.63, 0.81, 1.10, 1.48, 1.91, 2.09, &
& 2.29, 2.40, 2.45, 2.43, 2.37, 2.29, 2.14, 2.00, 1.86, 1.74, 1.62, &
& 1.51, 1.40, 1.29, 1.19, 1.10, 1.00, 0.92, 0.85, 0.79, 0.72, 0.67, &
& 0.62, 0.58, 0.54, 0.50, 0.47, 0.44, 0.41, 0.38, 0.35, 0.33, 0.32, &
& 0.30, 0.29, 0.29, 0.28, 0.28, 0.28, 0.28, 0.28, 0.29, 0.29, 0.30, &
& 0.31, 0.32, 0.33, 0.35, 0.37, 0.39, 0.43, 0.46, 0.49, 0.52, 0.56, &
& 0.61, 0.66, 0.72, 0.79, 0.86, 0.95, 1.05, 1.15, 1.26, 1.38, 1.50, &
& 1.60, 1.70, 1.78, 1.84, 1.88, 1.91, 1.91, 1.88, 1.86, 1.84, 1.80, &
& 1.74, 1.70, 1.64, 1.58, 1.55, 1.51, 1.48, 1.45, 1.43, 1.41, 1.39, &
& 1.38, 1.36, 1.35, 1.34, 1.33, 1.32, 1.31, 1.30, 1.29, 1.29, 1.28 /) 
! Chen data from 0 to -30C and given each -0.25C
!
!
!*       1.2  Find indexes for Gamma
!
DO II = 1, SIZE(PT,1)
  IDX_CHEN(II) = MINLOC(ABS(ZTEMP_CHEN(:)-ZTC(II)),1)
  IF (ZTC(II) .GT. ZTEMP_CHEN(IDX_CHEN(II)) .AND. IDX_CHEN(II) .GT. 1) THEN
    IDX_CHEN(II) = IDX_CHEN(II) - 1
  END IF
END DO
!
! linear interpolation to find Gamma given the temperature
!
ZGAM_CHEN_T(:) = ZGAM_CHEN(IDX_CHEN(:))
!
WHERE (IDX_CHEN(:) .LT. 121 .AND. IDX_CHEN(:) .GT. 0)
  ZGAM_CHEN_T(:) = ZGAM_CHEN(IDX_CHEN(:)) +                              &
                                    (ZTC(:) - ZTEMP_CHEN(IDX_CHEN(:))) * &
                  (ZGAM_CHEN(IDX_CHEN(:)+1) - ZGAM_CHEN(IDX_CHEN(:))) /  &
                 (ZTEMP_CHEN(IDX_CHEN(:)+1) - ZTEMP_CHEN(IDX_CHEN(:)) )
END WHERE
!
WHERE (IDX_CHEN(:) .EQ. SIZE(ZTEMP_CHEN,1))
  ZGAM_CHEN_T(:) = ZGAM_CHEN(SIZE(ZTEMP_CHEN,1))
END WHERE
!
! Gamma values for plates, columns and irregulars (from Figure 3, Chen and Lamb, 1984)
ZGAMMA_PLA = 0.5
ZGAMMA_COL = 1.5
ZGAMMA_IRR = 1.0
!
!
!*       2.   RATE OF TRANSFERT BETWEEN SHAPES
!             --------------------------------
!!++cb-- 2/11/20 : rajout d'une condition sur la temperature (transfert
!uniquement pour les temperatures > -20C, au dela, polycristaux)
!
! proportion of colums to transfer to irregulars
ZFACTOR_COL_IRR(:) = 0.
WHERE (YSHAPE_NAME(:) == 'PLA' .AND. PCI_SHAPE(:,2) .GT. PCTMIN(4))
  ZFACTOR_COL_IRR(:) = (ZGAM_CHEN_T(:) - ZGAMMA_IRR) / (ZGAMMA_PLA - ZGAMMA_IRR)
END WHERE
!
! the factor must be between 0 and 1
ZFACTOR_COL_IRR(:) = MAX(ZFACTOR_COL_IRR(:), 0.0) 
ZFACTOR_COL_IRR(:) = MIN(ZFACTOR_COL_IRR(:), 1.0)
!
! proportion of plates to transfer to irregulars
ZFACTOR_PLA_IRR(:) = 0.
WHERE (YSHAPE_NAME(:) == 'COL' .AND. PCI_SHAPE(:,1) .GT. PCTMIN(4))
  ZFACTOR_PLA_IRR(:) = (ZGAM_CHEN_T(:) - ZGAMMA_IRR) / (ZGAMMA_COL - ZGAMMA_IRR)
END WHERE
!
! the factor must be between 0 and 1
ZFACTOR_PLA_IRR(:) = MAX(ZFACTOR_PLA_IRR(:), 0.0)
ZFACTOR_PLA_IRR(:) = MIN(ZFACTOR_PLA_IRR(:), 1.0)
!
!
!
!*       4.   CHANGE IN NUMBER CONCENTRATION
!             ------------------------------
! 
ZZW1(:) = PCI_SHAPE(:,1) * ZFACTOR_PLA_IRR(:)
ZZW2(:) = PCI_SHAPE(:,2) * ZFACTOR_COL_IRR(:)
!
! for budgets
P_SHCI_HACH(:,1) = -ZZW1(:)
P_SHCI_HACH(:,2) = -ZZW2(:)
P_SHCI_HACH(:,3) = ZZW1(:) + ZZW2(:)  ! add to bullet-rosettes
!
IF (LHOOK) CALL DR_HOOK('LIMA_CHANGE_SHAPE', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_CHANGE_SHAPE
!
END MODULE MODE_LIMA_CHANGE_SHAPE
