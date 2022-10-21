!      #################################
       MODULE MODI_LIMA_DROPLETS_HOM_FREEZING
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_DROPLETS_HOM_FREEZING (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                          PT, PLVFACT, PLSFACT,                   &
                                          PRCT, PCCT, PLBDC,                      &
                                          P_TH_HONC, P_RC_HONC, P_CC_HONC,        &
                                          PA_TH, PA_RC, PA_CC, PA_RI, PA_CI       )
!
REAL,                 INTENT(IN)    :: PTSTEP 
CHARACTER(LEN=*),     INTENT(IN)    :: HFMFILE 
LOGICAL,              INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PT       ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT  ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT  ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! Cloud water lambda
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_HONC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_HONC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_HONC
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CI
!
END SUBROUTINE LIMA_DROPLETS_HOM_FREEZING
END INTERFACE
END MODULE MODI_LIMA_DROPLETS_HOM_FREEZING
!
!     ######################################################################
      SUBROUTINE LIMA_DROPLETS_HOM_FREEZING (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                             PT, PLVFACT, PLSFACT,                   &
                                             PRCT, PCCT, PLBDC,                      &
                                             P_TH_HONC, P_RC_HONC, P_CC_HONC,        &
                                             PA_TH, PA_RC, PA_CC, PA_RI, PA_CI       )
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
USE MODD_CST,             ONLY : XTT
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, XNUC
USE MODD_PARAM_LIMA_COLD, ONLY : XC_HONC, XTEXP1_HONC, XTEXP2_HONC, XTEXP3_HONC,   &
                                 XTEXP4_HONC, XTEXP5_HONC 
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                 INTENT(IN)    :: PTSTEP 
CHARACTER(LEN=*),     INTENT(IN)    :: HFMFILE 
LOGICAL,              INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PT       ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT  ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT  ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! Cloud water lambda
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_HONC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_HONC
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_HONC
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CC
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CI
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PT)) ::  ZZW, ZZX, ZZY, ZTCELSIUS
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PRELIMINARY COMPUTATIONS
!	        ------------------------
!
!
P_TH_HONC(:) = 0.
P_RC_HONC(:) = 0.
P_CC_HONC(:) = 0.
!
WHERE ( (PT(:)<XTT-35.0) .AND. (PCCT(:)>XCTMIN(2)) .AND. (PRCT(:)>XRTMIN(2)) )
   ZTCELSIUS(:) = PT(:)-XTT                                    ! T [°C]
!
!
!
!-------------------------------------------------------------------------------
!
!
!*       3.     Cloud droplets homogeneous freezing
!	        -----------------------------------
!
!
!  Compute the droplet homogeneous nucleation source: RCHONI
!                 -> Pruppacher(1995)
!
   ZZW(:) = 0.0
   ZZX(:) = 0.0
   ZZY(:) = 0.0

   ZZX(:) = 1.0 / ( 1.0 + (XC_HONC/PLBDC(:))*PTSTEP*                      &
           EXP( XTEXP1_HONC + ZTCELSIUS(:)*(                        &
           XTEXP2_HONC + ZTCELSIUS(:)*(                        &
           XTEXP3_HONC + ZTCELSIUS(:)*(                        &
           XTEXP4_HONC + ZTCELSIUS(:)*XTEXP5_HONC))) ) )**XNUC
!
   ZZW(:) = PCCT(:) * (1.0 - ZZX(:))                                  ! CCHONI
   ZZY(:) = PRCT(:) * (1.0 - ZZX(:))                                  ! RCHONI
!
   P_RC_HONC(:) = - ZZY(:)/PTSTEP
   P_CC_HONC(:) = - ZZW(:)/PTSTEP
   P_TH_HONC(:) = P_RC_HONC(:) * (PLSFACT(:)-PLVFACT(:))
!
   PA_TH(:) = PA_TH(:) + P_TH_HONC(:)
   PA_RC(:) = PA_RC(:) + P_RC_HONC(:)
   PA_CC(:) = PA_CC(:) + P_CC_HONC(:)
   PA_RI(:) = PA_RI(:) - P_RC_HONC(:)
   PA_CI(:) = PA_CI(:) - P_CC_HONC(:)
!
END WHERE
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPLETS_HOM_FREEZING
