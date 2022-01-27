!     ######spl
     MODULE MODI_SBL_DEPTH  
!    ################ 
!
INTERFACE
!
      SUBROUTINE SBL_DEPTH(KKB,KKE,PZZ,PFLXU,PFLXV,PWTHV,PLMO,PSBL_DEPTH)
!
INTEGER,                INTENT(IN)    :: KKB       ! first physical level
INTEGER,                INTENT(IN)    :: KKE       ! upper physical level
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PZZ       ! altitude of flux levels
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFLXU     ! u'w'
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFLXV     ! v'w'
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PWTHV     ! buoyancy flux
REAL, DIMENSION(:,:),   INTENT(IN)    :: PLMO      ! Monin-Obukhov length
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSBL_DEPTH! boundary layer height
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SBL_DEPTH
!
END INTERFACE
!
END MODULE MODI_SBL_DEPTH
