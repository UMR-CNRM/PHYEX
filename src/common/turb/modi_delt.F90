!     ######spl
     MODULE MODI_DELT
!    ################ 
!
IMPLICIT NONE
INTERFACE
!
  SUBROUTINE DELT (PLM, PWORK2, D, PWORK1, O2D, PZZ, PDYY, PDIRCOSZW,  &
  & GOCEAN, TURBN, PDXX, ODZ)
    !     ####################
    !!
    !!****  *DELT* routine to compute mixing length for DELT case
    !
    !!    AUTHOR
    !!    ------
    !!
    !!     M Tomasini      *Meteo-France
    !!
    !!    MODIFICATIONS
    !!    -------------
    !!      Original   01/05
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !              ------------
    !
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_CST, ONLY: CST_t
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
    USE MODD_TURB_n, ONLY: TURB_t
    USE MODE_SHUMAN_PHY, ONLY: MXF_PHY,MYF_PHY

    IMPLICIT NONE
    !*       0.1   Declarations of dummy arguments
    !
    TYPE(DIMPHYEX_t), INTENT(IN) :: D
    REAL, INTENT(OUT), DIMENSION(D%NIJT, D%NKT) :: PLM
    LOGICAL, INTENT(IN) :: ODZ
    REAL, INTENT(INOUT) :: PWORK2(D%NIJT, D%NKT)
    REAL, INTENT(INOUT) :: PWORK1(D%NIJT, D%NKT)
    LOGICAL, INTENT(IN) :: O2D
    REAL, INTENT(IN) :: PZZ(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PDYY(D%NIJT, D%NKT)
    REAL, INTENT(IN) :: PDIRCOSZW(D%NIJT)
    LOGICAL, INTENT(INOUT) :: GOCEAN
    TYPE(TURB_t), INTENT(IN) :: TURBN
    REAL, INTENT(IN) :: PDXX(D%NIJT, D%NKT)

  END SUBROUTINE DELT

  END INTERFACE

  END MODULE MODI_DELT
