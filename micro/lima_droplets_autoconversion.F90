!      #################################
       MODULE MODI_LIMA_DROPLETS_AUTOCONVERSION
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_DROPLETS_AUTOCONVERSION (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                            PRHODREF,                       &
                                            PRCT, PLBDC, PLBDR,             &
                                            P_RC_AUTO, P_CC_AUTO, P_CR_AUTO,&
                                            PA_RC, PA_CC, PA_RR, PA_CR      )
!
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE 
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR   ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_AUTO
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_AUTO
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_AUTO
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RR
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CR
!
END SUBROUTINE LIMA_DROPLETS_AUTOCONVERSION
END INTERFACE
END MODULE MODI_LIMA_DROPLETS_AUTOCONVERSION
!
!     ######################################################################
      SUBROUTINE LIMA_DROPLETS_AUTOCONVERSION (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                               PRHODREF,                       &
                                               PRCT, PLBDC, PLBDR,             &
                                               P_RC_AUTO, P_CC_AUTO, P_CR_AUTO,&
                                               PA_RC, PA_CC, PA_RR, PA_CR      )
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
USE MODD_PARAM_LIMA,      ONLY : XRTMIN
USE MODD_PARAM_LIMA_WARM, ONLY : XLAUTR, XAUTO1, XLAUTR_THRESHOLD, &
                                 XITAUTR, XAUTO2, XITAUTR_THRESHOLD, &
                                 XACCR4, XACCR5, XACCR3, XACCR1, XAC
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
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR   ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_AUTO
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_AUTO
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CR_AUTO
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RR
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CR
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRCT)) :: ZW1, ZW2, ZW3 ! work arrays
!
!-------------------------------------------------------------------------------
!
!
!
!*       3. Autoconversion of cloud droplets (Berry-Reinhardt parameterization)
!   	 ----------------------------------------------------------------------
!
!
!
P_RC_AUTO(:) = 0.0
P_CC_AUTO(:) = 0.0
P_CR_AUTO(:) = 0.0
!
ZW3(:) = 0.0
ZW2(:) = 0.0
ZW1(:) = 0.0
WHERE( PRCT(:)>XRTMIN(2) .AND. LDCOMPUTE(:) )
   ZW2(:) = MAX( 0.0, &
                     XLAUTR*PRHODREF(:)*PRCT(:)*(XAUTO1/PLBDC(:)**4-XLAUTR_THRESHOLD) ) ! L 
!
   ZW3(:) = MAX( 0.0, &
                     XITAUTR*ZW2(:)*PRCT(:)*(XAUTO2/PLBDC(:)-XITAUTR_THRESHOLD) ) ! L/tau
!
   P_RC_AUTO(:) = - ZW3(:)
!
   ZW1(:) = MIN( MIN( 1.2E4, &
                          (XACCR4/PLBDC(:)-XACCR5)/XACCR3 ), &
                          PLBDR(:)/XACCR1 ) ! D**-1 threshold diameter for 
                                                ! switching the autoconversion regimes
                                                ! min (80 microns, D_h, D_r)
   ZW3(:) = ZW3(:) * MAX( 0.0,ZW1(:) )**3 / XAC 
!
   P_CC_AUTO(:) = 0.
   P_CR_AUTO(:) = ZW3(:)
!
   PA_RC(:) = PA_RC(:) + P_RC_AUTO(:)
   PA_CC(:) = PA_CC(:) 
   PA_RR(:) = PA_RR(:) - P_RC_AUTO(:)
   PA_CR(:) = PA_CR(:) + P_CR_AUTO(:)
END WHERE
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPLETS_AUTOCONVERSION
