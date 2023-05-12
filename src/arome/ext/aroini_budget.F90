!     ######spl
SUBROUTINE AROINI_BUDGET(LDBU_ENABLE)
!
!**** *AROINI_BUDGET*   - Initialize common meso_NH MODD_ used in BUDGET for AROME

!     Purpose.
!     --------
!            Set implicit default values for MODD_BUDGET for the use in AROME 

!**   Interface.
!     ----------
!        *CALL* *AROINI_BUDGET

!        Explicit arguments :
!        --------------------
!       None

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        See documentation
!        To use budgets in DDH for AROME, budget must have type CART. First dimension is NPROMA and
!        second dimension is 1. Budgets are reset after each tipe step. Processes not used in AROME are
!        marked with 3

!     Externals.
!     ----------

!     Reference.
!     ----------
!
!     Author.
!     -------
!        Y. Seity 

!     Modifications.
!     --------------
!        Original :    03-12-12
!        T. Kovacic    05-04-27   Initialization for DDH
!        S. Riette     July-22    Simplification
!     ------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
USE MODD_BUDGET, ONLY : LBU_ENABLE, LBUDGET_U, LBUDGET_V, LBUDGET_W, LBUDGET_SV, &
                      & LBUDGET_TH, LBUDGET_TKE, &
                      & LBUDGET_RV, LBUDGET_RC, LBUDGET_RR, LBUDGET_RI, LBUDGET_RS, &
                      & LBUDGET_RG, LBUDGET_RH, TBUCONF_ASSOCIATE
!
IMPLICIT NONE
!
!*       0.1   declarations of argument
!
!
LOGICAL, INTENT(IN) :: LDBU_ENABLE
!
!*       0.2   Declarations of local variables :
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('AROINI_BUDGET',0,ZHOOK_HANDLE)
!
CALL TBUCONF_ASSOCIATE()
!
LBU_ENABLE = LDBU_ENABLE
!
LBUDGET_U  = LBU_ENABLE                                                                                                                
LBUDGET_V  = LBU_ENABLE                                                                                                                
LBUDGET_W  = LBU_ENABLE                                                                                                                
LBUDGET_TH = LBU_ENABLE                                                                                                                
LBUDGET_TKE= LBU_ENABLE                                                                                                               
LBUDGET_RV = LBU_ENABLE                                                                                                                
LBUDGET_RC = LBU_ENABLE                                                                                                                
LBUDGET_RR = LBU_ENABLE                                                                                                               
LBUDGET_RI = LBU_ENABLE                                                                                                               
LBUDGET_RS = LBU_ENABLE                                                                                                               
LBUDGET_RG = LBU_ENABLE                                                                                                               
LBUDGET_RH = LBU_ENABLE                                                                                                               
LBUDGET_SV = .FALSE. 
!
IF (LHOOK) CALL DR_HOOK('AROINI_BUDGET',1,ZHOOK_HANDLE)
!
END SUBROUTINE AROINI_BUDGET
