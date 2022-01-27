!     ######spl
      MODULE MODI_GRADIENT_W
!     ######################
!
INTERFACE
!
!            
FUNCTION GZ_W_M(PA,PDZZ, KKA, KKU, KL)      RESULT(PGZ_W_M)
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_W_M ! result mass point
!
END FUNCTION GZ_W_M
!            
FUNCTION GX_W_UW(PA,PDXX,PDZZ,PDZX, KKA, KKU, KL)      RESULT(PGX_W_UW)
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_W_UW ! result UW point
!
END FUNCTION GX_W_UW
!
!            
FUNCTION GY_W_VW(PA,PDXX,PDZZ,PDZX, KKA, KKU, KL)      RESULT(PGY_W_VW)
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_W_VW ! result VW point
!
END FUNCTION GY_W_VW
!
!
END INTERFACE
!
END MODULE MODI_GRADIENT_W
