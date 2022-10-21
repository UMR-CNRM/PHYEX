!      #################################
       MODULE MODI_LIMA_ICE_AGGREGATION_SNOW
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_ICE_AGGREGATION_SNOW (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                         PT,                             &
                                         PRIT, PRST, PCIT, PLBDI, PLBDS, &
                                         P_RI_AGGS, P_CI_AGGS,           &
                                         PA_RI, PA_CI, PA_RS             )
!
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE 
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PT
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT
REAL, DIMENSION(:),   INTENT(IN)    :: PRST
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDI 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_AGGS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_AGGS
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RS
!
END SUBROUTINE LIMA_ICE_AGGREGATION_SNOW
END INTERFACE
END MODULE MODI_LIMA_ICE_AGGREGATION_SNOW
!
!     ######################################################################
      SUBROUTINE LIMA_ICE_AGGREGATION_SNOW (HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                            PT,                             &
                                            PRIT, PRST, PCIT, PLBDI, PLBDS, &
                                            P_RI_AGGS, P_CI_AGGS,           &
                                            PA_RI, PA_CI, PA_RS             )
!     ######################################################################
!
!!    PURPOSE
!!    -------
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
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN
USE MODD_PARAM_LIMA_COLD, ONLY : XBI, XCCS, XCXS, XCOLEXIS, XAGGS_CLARGE1, XAGGS_CLARGE2, &
                                 XAGGS_RLARGE1, XAGGS_RLARGE2
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE 
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PT
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT
REAL, DIMENSION(:),   INTENT(IN)    :: PRST
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDI 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_AGGS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_AGGS
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RS
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRIT)) :: ZZW1, ZZW2, ZZW3 ! work arrays
!
!-------------------------------------------------------------------------------
!
!
!*       2.4    Aggregation of r_i on r_s: CIAGGS and RIAGGS
!        ---------------------------------------------------
!
ZZW1(:) = 0.
ZZW2(:) = 0.
ZZW3(:) = 0.
!
P_RI_AGGS(:) = 0.
P_CI_AGGS(:) = 0.
!
!!$print *, "aggregation of i on s"
!!$print *, "ri", PRIT(156,2,22)
!!$print *, "Ni", PCIT(156,2,22)
!!$print *, "rs", PRST(156,2,22)
!!$print *, "lambda i", PLBDI(156,2,22)
!!$print *, "lambda s", PLBDS(156,2,22)
!!$print *, "T", PT(156,2,22)
!!$print *, "C1", XAGGS_CLARGE1
!!$print *, "C2", XAGGS_CLARGE2
!!$print *, "R1", XAGGS_RLARGE1
!!$print *, "R2", XAGGS_RLARGE2
!
WHERE ( (PRIT(:)>XRTMIN(4)) .AND. (PRST(:)>XRTMIN(5)) .AND. LDCOMPUTE(:) )
   ZZW1(:) = (PLBDI(:) / PLBDS(:))**3
   ZZW2(:) = (PCIT(:)*(XCCS*PLBDS(:)**XCXS)*EXP( XCOLEXIS*(PT(:)-XTT) )) &
        / (PLBDI(:)**3)
   ZZW3(:) = ZZW2(:)*(XAGGS_CLARGE1+XAGGS_CLARGE2*ZZW1(:))
!
   P_CI_AGGS(:) = - ZZW3(:)
!
   ZZW2(:) = ZZW2(:) / PLBDI(:)**XBI
   ZZW2(:) = ZZW2(:)*(XAGGS_RLARGE1+XAGGS_RLARGE2*ZZW1(:))
!
   P_RI_AGGS(:) = - ZZW2(:)
END WHERE
!
!!$print *, "tendance ci", P_CI_AGGS(156,2,22)
!!$print *, "tendance ri", P_RI_AGGS(156,2,22)
!
PA_RI(:) = PA_RI(:) + P_RI_AGGS(:)
PA_CI(:) = PA_CI(:) + P_CI_AGGS(:)
PA_RS(:) = PA_RS(:) - P_RI_AGGS(:)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_ICE_AGGREGATION_SNOW
