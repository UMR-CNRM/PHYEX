!      #################################
       MODULE MODI_LIMA_BERGERON
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_BERGERON (HFMFILE, OCLOSE_OUT, LDCOMPUTE,    &
                             PRCT, PRIT, PCIT, PLBDI,           &
                             PSSIW, PAI, PCJ, PLVFACT, PLSFACT, &
                             P_TH_BERFI, P_RC_BERFI,            &
                             PA_TH, PA_RC, PA_RI                )
!
CHARACTER(LEN=*),     INTENT(IN)    :: HFMFILE 
LOGICAL,              INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDI   ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PSSIW   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PAI     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_BERFI
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_BERFI
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RI
!!
END SUBROUTINE LIMA_BERGERON
END INTERFACE
END MODULE MODI_LIMA_BERGERON
!
!     ######################################################################
      SUBROUTINE LIMA_BERGERON(HFMFILE, OCLOSE_OUT, LDCOMPUTE,    &
                               PRCT, PRIT, PCIT, PLBDI,           &
                               PSSIW, PAI, PCJ, PLVFACT, PLSFACT, &
                               P_TH_BERFI, P_RC_BERFI,            &
                               PA_TH, PA_RC, PA_RI                )
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
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCTMIN
USE MODD_PARAM_LIMA_COLD,  ONLY : XDI, X0DEPI, X2DEPI
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=*),     INTENT(IN)    :: HFMFILE 
LOGICAL,              INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDI   ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PSSIW   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PAI     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_BERFI
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_BERFI
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RI
!
!*       0.2   Declarations of local variables :
!
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PRELIMINARY COMPUTATIONS
!	        ------------------------
!
P_TH_BERFI(:) = 0.0
P_RC_BERFI(:) = 0.0
!
WHERE( (PRCT(:)>XRTMIN(2)) .AND. (PRIT(:)>XRTMIN(4)) .AND. (PCIT(:)>XCTMIN(4)) .AND. LDCOMPUTE(:))
!   ZZW(:) = EXP( (XALPW-XALPI) - (XBETAW-XBETAI)/ZZT(:)          &
!        - (XGAMW-XGAMI)*ALOG(ZZT(:)) ) -1.0 
! supersaturation over ice at water saturation
   P_RC_BERFI(:) = - ( PSSIW(:) / PAI(:) ) * PCIT(:) *        &
        ( X0DEPI/PLBDI(:)+X2DEPI*PCJ(:)*PCJ(:)/PLBDI(:)**(XDI+2.0) )
   P_TH_BERFI(:) = - P_RC_BERFI(:)*(PLSFACT(:)-PLVFACT(:))
END WHERE
!
PA_RC(:) = PA_RC(:) + P_RC_BERFI(:)
PA_RI(:) = PA_RI(:) - P_RC_BERFI(:)
PA_TH(:) = PA_TH(:) + P_TH_BERFI(:)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_BERGERON
