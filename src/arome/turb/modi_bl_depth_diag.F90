!     ######spl
     MODULE MODI_BL_DEPTH_DIAG  
!    ################ 
!
INTERFACE BL_DEPTH_DIAG  
!
!
      FUNCTION BL_DEPTH_DIAG_3D(KKB,KKE,PSURF,PZS,PFLUX,PZZ,PFTOP_O_FSURF)

INTEGER,                INTENT(IN)           :: KKB          ! bottom point
INTEGER,                INTENT(IN)           :: KKE          ! top point
REAL, DIMENSION(:,:),   INTENT(IN)           :: PSURF        ! surface flux
REAL, DIMENSION(:,:),   INTENT(IN)           :: PZS          ! orography
REAL, DIMENSION(:,:,:), INTENT(IN)           :: PFLUX        ! flux
REAL, DIMENSION(:,:,:), INTENT(IN)           :: PZZ          ! altitude of flux points
REAL,                   INTENT(IN)           :: PFTOP_O_FSURF! Flux at BL top / Surface flux
REAL, DIMENSION(SIZE(PSURF,1),SIZE(PSURF,2)) :: BL_DEPTH_DIAG_3D
!
END FUNCTION BL_DEPTH_DIAG_3D
!
!
      FUNCTION BL_DEPTH_DIAG_1D(KKB,KKE,PSURF,PZS,PFLUX,PZZ,PFTOP_O_FSURF)
INTEGER,                INTENT(IN)           :: KKB          ! bottom point
INTEGER,                INTENT(IN)           :: KKE          ! top point
REAL,                   INTENT(IN)           :: PSURF        ! surface flux
REAL,                   INTENT(IN)           :: PZS          ! orography
REAL, DIMENSION(:),     INTENT(IN)           :: PFLUX        ! flux
REAL, DIMENSION(:),     INTENT(IN)           :: PZZ          ! altitude of flux points
REAL,                   INTENT(IN)           :: PFTOP_O_FSURF! Flux at BL top / Surface flux
REAL                                         :: BL_DEPTH_DIAG_1D
!
END FUNCTION BL_DEPTH_DIAG_1D
!
END INTERFACE
!
END MODULE MODI_BL_DEPTH_DIAG
