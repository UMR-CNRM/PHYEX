!      #################################
       MODULE MODI_LIMA_DROPLETS_SELF_COLLECTION
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                             PRHODREF,                       &
                                             PCCT, PLBDC3,                   &
                                             P_CC_SELF,                      &
                                             PA_CC                           )
!
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE 
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC3  ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_SELF
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CC
!
END SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION
END INTERFACE
END MODULE MODI_LIMA_DROPLETS_SELF_COLLECTION
!
!     ######################################################################
      SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                                PRHODREF,                       &
                                                PCCT, PLBDC3,                   &
                                                P_CC_SELF,                      &
                                                PA_CC                           )
!     ######################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the cold-phase homogeneous
!!    freezing of CCN, droplets and drops (T<-35°C)
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      C. Barthe  * LACy*   jan. 2014  add budgets
!!      B.Vie 10/2016 Bug zero division
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA,      ONLY : XCTMIN
USE MODD_PARAM_LIMA_WARM, ONLY : XSELFC
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE 
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC3  ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_SELF
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CC
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PCCT)) :: ZW ! work arrays
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PRELIMINARY COMPUTATIONS
!	        ------------------------
!
!
P_CC_SELF(:)=0.
!
WHERE( PCCT(:)>XCTMIN(2) .AND. LDCOMPUTE(:) )
   ZW(:) = XSELFC*(PCCT(:)/PLBDC3(:))**2 * PRHODREF(:) ! analytical integration
   P_CC_SELF(:) = - ZW(:)
   PA_CC(:) = PA_CC(:) + P_CC_SELF(:)
END WHERE
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION
