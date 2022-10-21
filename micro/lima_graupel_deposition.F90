!      #################################
       MODULE MODI_LIMA_GRAUPEL_DEPOSITION
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_GRAUPEL_DEPOSITION (HFMFILE, OCLOSE_OUT, LDCOMPUTE,       &
                                       PRGT, PSSI, PLBDG, PAI, PCJ, PLSFACT, &
                                       P_TH_DEPG, P_RG_DEPG,                 &
                                       PA_TH, PA_RV, PA_RG                   )
!
CHARACTER(LEN=*),     INTENT(IN)    :: HFMFILE 
LOGICAL,              INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PSSI    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDG   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PAI     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_DEPG
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RV
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RG
!!
END SUBROUTINE LIMA_GRAUPEL_DEPOSITION
END INTERFACE
END MODULE MODI_LIMA_GRAUPEL_DEPOSITION
!
!     ######################################################################
      SUBROUTINE LIMA_GRAUPEL_DEPOSITION (HFMFILE, OCLOSE_OUT, LDCOMPUTE,       &
                                          PRGT, PSSI, PLBDG, PAI, PCJ, PLSFACT, &
                                          P_TH_DEPG, P_RG_DEPG,                 &
                                          PA_TH, PA_RV, PA_RG                   )
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
USE MODD_PARAM_LIMA,       ONLY : XRTMIN
USE MODD_PARAM_LIMA_MIXED, ONLY : X0DEPG, XEX0DEPG, X1DEPG, XEX1DEPG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE 
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PSSI    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDG   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PAI     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPG
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RG_DEPG
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RV
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RG
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
P_TH_DEPG(:) = 0.0
P_RG_DEPG(:) = 0.0
WHERE ( (PRGT(:)>XRTMIN(6)) .AND. LDCOMPUTE(:) )
   !Correction BVIE RHODREF
   !      ZZW(:) = ( ZSSI(:)/(ZRHODREF(:)*ZAI(:)) ) *                               &
   P_RG_DEPG(:) = ( PSSI(:)/(PAI(:)) ) *                               &
        ( X0DEPG*PLBDG(:)**XEX0DEPG + X1DEPG*PCJ(:)*PLBDG(:)**XEX1DEPG )
   P_TH_DEPG(:) = P_RG_DEPG(:)*PLSFACT(:)
END WHERE
!
PA_RV(:) = PA_RV(:) - P_RG_DEPG(:)
PA_RG(:) = PA_RG(:) + P_RG_DEPG(:)
PA_TH(:) = PA_TH(:) + P_TH_DEPG(:)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_GRAUPEL_DEPOSITION
