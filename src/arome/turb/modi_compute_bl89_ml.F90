!     ######spl
     MODULE MODI_COMPUTE_BL89_ML
!    ###########################

INTERFACE

!     ###################################################################
      SUBROUTINE COMPUTE_BL89_ML(KKA,KKB,KKE,KKU,KKL,PDZZ2D, &
             PTKEM_DEP,PG_O_THVREF,PVPT,KK,OUPORDN,OFLUX,PSHEAR,PLWORK)
!     ###################################################################

!*               1.1  Declaration of Arguments

INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:),   INTENT(IN)  :: PDZZ2D        ! height difference between two mass levels
REAL, DIMENSION(:),     INTENT(IN)  :: PTKEM_DEP     ! TKE to consume
REAL, DIMENSION(:),     INTENT(IN)  :: PG_O_THVREF   ! g/ThetaVRef at the departure point
REAL, DIMENSION(:,:),   INTENT(IN)  :: PVPT          ! ThetaV on mass levels
INTEGER,                INTENT(IN)  :: KK            ! index of departure level
LOGICAL,                INTENT(IN)  :: OUPORDN       ! switch to compute upward (true) or
                                                     !   downward (false) mixing length
LOGICAL,                INTENT(IN)  :: OFLUX         ! Computation must be done from flux level
REAL, DIMENSION(:),     INTENT(OUT) :: PLWORK        ! Resulting mixing length
REAL, DIMENSION(:,:),   INTENT(IN)  :: PSHEAR        ! vertical wind shear for RM17 mixing length

END SUBROUTINE COMPUTE_BL89_ML

END INTERFACE
!
END MODULE MODI_COMPUTE_BL89_ML
