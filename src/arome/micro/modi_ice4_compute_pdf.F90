MODULE MODI_ICE4_COMPUTE_PDF
INTERFACE
SUBROUTINE ICE4_COMPUTE_PDF(KSIZE, HSUBG_AUCV, HSUBG_PR_PDF, &
                            PRHODREF, PRCT, PCF, PSIGMA_RC,&
                            PHLC_HCF, PHLC_LCF, PHLC_HRC, PHLC_LRC, PRF)
IMPLICIT NONE
INTEGER,                INTENT(IN)  :: KSIZE
CHARACTER(LEN=4),       INTENT(IN)  :: HSUBG_AUCV     ! Kind of Subgrid autoconversion method
CHARACTER*80,           INTENT(IN)  :: HSUBG_PR_PDF   ! pdf for subgrid precipitation
REAL, DIMENSION(KSIZE), INTENT(IN)  :: PRHODREF   ! Reference density
REAL, DIMENSION(KSIZE), INTENT(IN)  :: PRCT       ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE), INTENT(IN)  :: PCF        ! Cloud fraction
REAL, DIMENSION(KSIZE), INTENT(IN)  :: PSIGMA_RC  ! Standard deviation of rc at time t
REAL, DIMENSION(KSIZE), INTENT(OUT) :: PHLC_HCF   ! HLCLOUDS : fraction of High Cloud Fraction in grid
REAL, DIMENSION(KSIZE), INTENT(OUT) :: PHLC_LCF   ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                                  !    note that PCF = PHLC_HCF + PHLC_LCF
REAL, DIMENSION(KSIZE), INTENT(OUT) :: PHLC_HRC   ! HLCLOUDS : LWC that is High LWC in grid
REAL, DIMENSION(KSIZE), INTENT(OUT) :: PHLC_LRC   ! HLCLOUDS : LWC that is Low  LWC in grid
                                                  !    note that PRC = PHLC_HRC + PHLC_LRC
REAL, DIMENSION(KSIZE), INTENT(OUT) :: PRF        ! Rain fraction
END SUBROUTINE ICE4_COMPUTE_PDF
END INTERFACE
END MODULE MODI_ICE4_COMPUTE_PDF
