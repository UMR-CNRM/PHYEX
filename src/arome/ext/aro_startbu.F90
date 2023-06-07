!     ######spl
SUBROUTINE ARO_STARTBU( KIDIA, KFDIA, KLEV, KRR,KSV,PRHODJ,&
                        & PRUS,PRVS,PRWS,PRTHS,PRRS,PRTKES,YDDDH, YDLDDH, YDMDDH)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

!     Purpose.
!     --------
!            Prepare budget arrays at the start of budget calculations.

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

!     Externals.
!     ----------

!     Reference.
!     ----------
!
!     Author.
!     -------
!        T. Kovacic

!     Modifications.
!     --------------
!        Original :    05-05-06
!        19-Sept-08: O.Riviere Removal of unecessary part for new diagnostic data flow
!        30-Janv-19: F.Voitus  new DDH superstructure + RR budget correction 
!     ------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
USE MODE_BUDGET_PHY, ONLY: BUDGET_DDH
USE DDH_MIX     , ONLY  : TYP_DDH
USE YOMLDDH     , ONLY  : TLDDH
USE YOMMDDH     , ONLY  : TMDDH

!
IMPLICIT NONE
!
!*       0.1   declarations of argument
!
INTEGER,                  INTENT(IN)  :: KIDIA
INTEGER,                  INTENT(IN)  :: KFDIA
INTEGER,                  INTENT(IN)  :: KLEV
INTEGER,                  INTENT(IN)  :: KRR     ! Number of moist variables
INTEGER,                  INTENT(IN)  :: KSV     ! Number of Scalar Variables
!
REAL, DIMENSION(KFDIA,1,KLEV),   INTENT(IN)  :: PRHODJ         ! (Rho) dry * Jacobian
!
REAL, DIMENSION(KFDIA,1,KLEV),   INTENT(IN) :: PRUS, PRVS, PRWS         ! Source
REAL, DIMENSION(KFDIA,1,KLEV),   INTENT(IN) :: PRTHS, PRTKES            !   -
REAL, DIMENSION(KFDIA,1,KLEV,KRR), INTENT(IN) :: PRRS              !  terms

TYPE(TYP_DDH)        , INTENT(INOUT) :: YDDDH
TYPE(TLDDH)        , INTENT(IN) :: YDLDDH
TYPE(TMDDH)        , INTENT(IN) :: YDMDDH

!
!
!*       0.2   Declarations of local variables :
!

LOGICAL   ::  LL_BUDGET_RR
INTEGER   ::  JR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE




    LL_BUDGET_RR = (LBUDGET_RV).OR.(LBUDGET_RC).OR.(LBUDGET_RR) &
                 & .OR.(LBUDGET_RI).OR.(LBUDGET_RS) & 
                 & .OR.(LBUDGET_RG).OR.(LBUDGET_RH) 

!
    IF (LHOOK) CALL DR_HOOK('ARO_STARTBU',0,ZHOOK_HANDLE)
    
    IF (LBUDGET_U)   CALL BUDGET_DDH (PRUS(:,:,:)*PRHODJ(:,:,:),1,'INIF_BU_RU',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_V)   CALL BUDGET_DDH (PRVS(:,:,:)*PRHODJ(:,:,:),2,'INIF_BU_RV',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_W)   CALL BUDGET_DDH (PRWS(:,:,:)*PRHODJ(:,:,:),3,'INIF_BU_RW',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_TH)  CALL BUDGET_DDH (PRTHS(:,:,:)*PRHODJ(:,:,:),4,'INIF_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_TKE) CALL BUDGET_DDH (PRTKES(:,:,:)*PRHODJ(:,:,:),5,'INIF_BU_RTKE',YDDDH, YDLDDH, YDMDDH)
    
    IF (LL_BUDGET_RR) THEN
       DO JR = 1,KRR
          CALL BUDGET_DDH (PRRS(:,:,:,JR)*PRHODJ(:,:,:),5+JR,'INIF_BU_RR',YDDDH, YDLDDH, YDMDDH)
       END DO
    END IF        
    

!
IF (LHOOK) CALL DR_HOOK('ARO_STARTBU',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_STARTBU
