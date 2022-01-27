!     ######spl
FUNCTION BL_DEPTH_DIAG_1D(KKB,KKE,PSURF,PZS,PFLUX,PZZ,PFTOP_O_FSURF)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
USE MODI_BL_DEPTH_DIAG_3D
IMPLICIT NONE
!
INTEGER,                INTENT(IN)           :: KKB          ! bottom point
INTEGER,                INTENT(IN)           :: KKE          ! top point
REAL,                   INTENT(IN)           :: PSURF        ! surface flux
REAL,                   INTENT(IN)           :: PZS          ! orography
REAL, DIMENSION(:),     INTENT(IN)           :: PFLUX        ! flux
REAL, DIMENSION(:),     INTENT(IN)           :: PZZ          ! altitude of flux points
REAL,                   INTENT(IN)           :: PFTOP_O_FSURF! Flux at BL top / Surface flux
REAL                                         :: BL_DEPTH_DIAG_1D
!
REAL, DIMENSION(1,1)             :: ZSURF
REAL, DIMENSION(1,1)             :: ZZS
REAL, DIMENSION(1,1,SIZE(PFLUX)) :: ZFLUX
REAL, DIMENSION(1,1,SIZE(PZZ))   :: ZZZ
REAL, DIMENSION(1,1)             :: ZBL_DEPTH_DIAG
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_1D',0,ZHOOK_HANDLE)
ZSURF        = PSURF
ZZS          = PZS
ZFLUX(1,1,:) = PFLUX(:)
ZZZ  (1,1,:) = PZZ  (:)
!
ZBL_DEPTH_DIAG = BL_DEPTH_DIAG_3D(KKB,KKE,ZSURF,ZZS,ZFLUX,ZZZ,PFTOP_O_FSURF)
!
BL_DEPTH_DIAG_1D = ZBL_DEPTH_DIAG(1,1)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_1D',1,ZHOOK_HANDLE)
END FUNCTION BL_DEPTH_DIAG_1D
