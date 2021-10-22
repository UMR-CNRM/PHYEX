!     ######spl
     MODULE MODI_COMPUTE_BL89_ML
!    ###########################

INTERFACE

!     ###################################################################
      SUBROUTINE COMPUTE_BL89_ML(KKA,KKB,KKE,KKU,KKL,PDZZ2D, &
             PTKEM_DEP,PG_O_THVREF,PVPT,KK,OUPORDN,OFLUX,PLWORK)
!     ###################################################################

!*               1.1  Declaration of Arguments

INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:),   INTENT(IN)  :: PDZZ2D
REAL, DIMENSION(:),     INTENT(IN)  :: PTKEM_DEP
REAL, DIMENSION(:),     INTENT(IN)  :: PG_O_THVREF
REAL, DIMENSION(:,:),   INTENT(IN)  :: PVPT
INTEGER,                INTENT(IN)  :: KK
LOGICAL,                INTENT(IN)  :: OUPORDN
LOGICAL,                INTENT(IN)  :: OFLUX         ! Computation must be done from flux level
REAL, DIMENSION(:),     INTENT(OUT) :: PLWORK

END SUBROUTINE COMPUTE_BL89_ML

END INTERFACE
!
END MODULE MODI_COMPUTE_BL89_ML
