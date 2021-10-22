SUBROUTINE ICE4_COMPUTE_PDF(KSIZE, HSUBG_AUCV, HSUBG_PR_PDF, &
                            PRHODREF, PRCT, PCF, PSIGMA_RC,&
                            PHLC_HCF, PHLC_LCF, PHLC_HRC, PHLC_LRC, PRF)
!!
!!**  PURPOSE
!!    -------
!!      Computes the pdf used to split cloud into high and low content parts
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the plitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!
!
!*      0. DECLARATIONS
!          ------------
!
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
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
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KSIZE) :: ZRCRAUTC,      & !RC value to begin rain formation =XCRIAUTC/RHODREF
                          ZHLC_RCMAX,    & !HLCLOUDS : maximum value for RC in distribution
                          ZHLC_LRCLOCAL, & !HLCLOUDS : LWC that is Low  LWC local in LCF
                          ZHLC_HRCLOCAL    !HLCLOUDS : LWC that is High LWC local in HCF
                                                    !    note that ZRC/CF = ZHLC_HRCLOCAL+ ZHLC_LRCLOCAL
                                                    !                     = PHLC_HRC/HCF+ PHLC_LRC/LCF
REAL :: ZCOEFFRCM
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_COMPUTE_PDF', 0, ZHOOK_HANDLE)!

!Cloud water split between high and low content part is done according to autoconversion option
ZRCRAUTC(:)=XCRIAUTC/PRHODREF(:) ! Autoconversion rc threshold
IF(HSUBG_AUCV=='NONE') THEN
  !Cloud water is entirely in low or high part
  WHERE(PRCT(:)>ZRCRAUTC(:))
    PHLC_HCF(:)=1.
    PHLC_LCF(:)=0.
    PHLC_HRC(:)=PRCT(:)
    PHLC_LRC(:)=0.
    PRF(:)     =1.
  ELSEWHERE(PRCT(:)>XRTMIN(2))
    PHLC_HCF(:)=0.
    PHLC_LCF(:)=1.
    PHLC_HRC(:)=0.
    PHLC_LRC(:)=PRCT(:)
    PRF(:)     =0.
  ELSEWHERE
    PHLC_HCF(:)=0.
    PHLC_LCF(:)=0.
    PHLC_HRC(:)=0.
    PHLC_LRC(:)=0.
    PRF(:)     =0.
  END WHERE

ELSEIF(HSUBG_AUCV=='CLFR') THEN
  !Cloud water is only in the cloudy part and entirely in low or high part
  WHERE(PCF(:)>0.)
    WHERE(PRCT(:)/PCF(:)>ZRCRAUTC(:))
      PHLC_HCF(:)=PCF(:)
      PHLC_LCF(:)=0.
      PHLC_HRC(:)=PRCT(:)
      PHLC_LRC(:)=0.
      PRF(:)     =PCF(:)
    ELSEWHERE(PRCT(:)>XRTMIN(2))
      PHLC_HCF(:)=0.
      PHLC_LCF(:)=PCF(:)
      PHLC_HRC(:)=0.0
      PHLC_LRC(:)=PRCT(:)
      PRF(:)     =0.
    ELSEWHERE
      PHLC_HCF(:)=0.
      PHLC_LCF(:)=0.
      PHLC_HRC(:)=0.
      PHLC_LRC(:)=0.
      PRF(:)     =0.
    END WHERE
  ELSEWHERE
    PHLC_HCF(:)=0.
    PHLC_LCF(:)=0.
    PHLC_HRC(:)=0.
    PHLC_LRC(:)=0.
    PRF(:)     =0.
  END WHERE

ELSEIF(HSUBG_AUCV=='PDF ') THEN
  !Cloud water is split between high and low part according to a PDF
  !    'HLCRECTPDF'    : rectangular PDF form
  !    'HLCTRIANGPDF'  : triangular PDF form
  !    'HLCQUADRAPDF'  : second order quadratic PDF form
  !    'HLCISOTRIPDF'  : isocele triangular PDF
  !    'SIGM'          : Redelsperger and Sommeria (1986)
  IF(HSUBG_PR_PDF=='SIGM') THEN
    ! Redelsperger and Sommeria (1986) but organised according to Turner (2011, 2012)
    WHERE (PRCT(:)>ZRCRAUTC(:)+PSIGMA_RC(:))
      PHLC_HCF(:)=1.
      PHLC_LCF(:)=0.
      PHLC_HRC(:)=PRCT(:)
      PHLC_LRC(:)=0.
      PRF(:)     =1.
    ELSEWHERE(PRCT(:)> (ZRCRAUTC(:)-PSIGMA_RC(:)) .AND. &
            & PRCT(:)<=(ZRCRAUTC(:)+PSIGMA_RC(:))       )
      PHLC_HCF(:)=(PRCT(:)+PSIGMA_RC(:)-ZRCRAUTC(:))/ &
                  &(2.*PSIGMA_RC(:))
      PHLC_LCF(:)=MAX(0., PCF(:)-PHLC_HCF(:))
      PHLC_HRC(:)=(PRCT(:)+PSIGMA_RC(:)-ZRCRAUTC(:))* &
                  &(PRCT(:)+PSIGMA_RC(:)+ZRCRAUTC(:))/ &
                  &(4.*PSIGMA_RC(:))
      PHLC_LRC(:)=MAX(0., PRCT(:)-PHLC_HRC(:))
      PRF(:)     =PHLC_HCF(:)
    ELSEWHERE(PRCT(:)>XRTMIN(2) .AND. PCF(:)>0.)
      PHLC_HCF(:)=0.
      PHLC_LCF(:)=PCF(:)
      PHLC_HRC(:)=0.
      PHLC_LRC(:)=PRCT(:)
      PRF(:)     =0.
    ELSEWHERE
      PHLC_HCF(:)=0.
      PHLC_LCF(:)=0.
      PHLC_HRC(:)=0.
      PHLC_LRC(:)=0.
      PRF(:)     =0.
    END WHERE
   ! Turner (2011, 2012)
  ELSEIF(HSUBG_PR_PDF=='HLCRECTPDF' .OR. HSUBG_PR_PDF=='HLCISOTRIPDF' .OR. &
         &HSUBG_PR_PDF=='HLCTRIANGPDF' .OR. HSUBG_PR_PDF=='HLCQUADRAPDF') THEN
    ! Calculate maximum value r_cM from PDF forms
    IF(HSUBG_PR_PDF=='HLCRECTPDF' .OR. HSUBG_PR_PDF=='HLCISOTRIPDF') THEN
      ZCOEFFRCM=2.
    ELSE IF(HSUBG_PR_PDF=='HLCTRIANGPDF') THEN
      ZCOEFFRCM=3.
    ELSE IF(HSUBG_PR_PDF=='HLCQUADRAPDF') THEN
      ZCOEFFRCM=4.
    END IF
    WHERE(PRCT(:).GT.0. .AND. PCF(:).GT.0.)
      ZHLC_RCMAX(:)=ZCOEFFRCM*PRCT(:)/PCF(:)
    END WHERE
    ! Split available water and cloud fraction in two parts
    ! Calculate local mean values int he low and high parts for the 3 PDF forms:
    IF(HSUBG_PR_PDF=='HLCRECTPDF') THEN
      WHERE(PRCT(:).GT.0. .AND. PCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
        ZHLC_LRCLOCAL(:)=0.5*ZRCRAUTC(:)
        ZHLC_HRCLOCAL(:)=( ZHLC_RCMAX(:) + ZRCRAUTC(:))/2.0
      END WHERE
    ELSE IF(HSUBG_PR_PDF=='HLCTRIANGPDF') THEN
      WHERE(PRCT(:).GT.0. .AND. PCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
        ZHLC_LRCLOCAL(:)=( ZRCRAUTC(:) *(3.0 * ZHLC_RCMAX(:) - 2.0 * ZRCRAUTC(:) ) ) &
                        / (3.0 * (2.0 * ZHLC_RCMAX(:) - ZRCRAUTC(:)  ) )
        ZHLC_HRCLOCAL(:)=(ZHLC_RCMAX(:) + 2.0*ZRCRAUTC(:)) / 3.0
      END WHERE
    ELSE IF(HSUBG_PR_PDF=='HLCQUADRAPDF') THEN
      WHERE(PRCT(:).GT.0. .AND. PCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
        ZHLC_LRCLOCAL(:)=(3.0 *ZRCRAUTC(:)**3 - 8.0 *ZRCRAUTC(:)**2 * ZHLC_RCMAX(:) &
                        + 6.0*ZRCRAUTC(:) *ZHLC_RCMAX(:)**2 ) &
                        / &
                        (4.0* ZRCRAUTC(:)**2 -12.0*ZRCRAUTC(:) *ZHLC_RCMAX(:) &
                        + 12.0 * ZHLC_RCMAX(:)**2 )
        ZHLC_HRCLOCAL(:)=(ZHLC_RCMAX(:) + 3.0*ZRCRAUTC(:))/4.0
      END WHERE
    ELSE IF(HSUBG_PR_PDF=='HLCISOTRIPDF') THEN
      WHERE(PRCT(:).GT.0. .AND. PCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
        WHERE((PRCT(:) / PCF(:)).LE.ZRCRAUTC(:))
          ZHLC_LRCLOCAL(:)=( (ZHLC_RCMAX(:))**3 &
                          -(12.0 * (ZHLC_RCMAX(:))*(ZRCRAUTC(:))**2) &
                          +(8.0 * ZRCRAUTC(:)**3) ) &
                          /( (6.0 * (ZHLC_RCMAX(:))**2) &
                          -(24.0 * (ZHLC_RCMAX(:)) * ZRCRAUTC(:)) &
                          +(12.0 * ZRCRAUTC(:)**2) )
          ZHLC_HRCLOCAL(:)=( ZHLC_RCMAX(:) + 2.0 * ZRCRAUTC(:) )/3.0
        ELSEWHERE
          ZHLC_LRCLOCAL(:)=(2.0/3.0) * ZRCRAUTC(:)
          ZHLC_HRCLOCAL(:)=(3.0*ZHLC_RCMAX(:)**3 - 8.0*ZRCRAUTC(:)**3) &
                          / (6.0 * ZHLC_RCMAX(:)**2 - 12.0*ZRCRAUTC(:)**2)
        END WHERE
      END WHERE
    END IF
    ! Compare r_cM  to r_cR to know if cloud water content is high enough to split in two parts or not
    WHERE (PRCT(:).GT.0. .AND. PCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
      ! Calculate final values for LCF and HCF:
      PHLC_LCF(:)=PCF(:) &
                    *(ZHLC_HRCLOCAL- &
                    (PRCT(:) / PCF(:))) &
                    / (ZHLC_HRCLOCAL-ZHLC_LRCLOCAL)
      PHLC_HCF(:)=MAX(0., PCF(:)-PHLC_LCF(:))
      !
      ! Calculate final values for LRC and HRC:
      PHLC_LRC(:)=ZHLC_LRCLOCAL*PHLC_LCF(:)
      PHLC_HRC(:)=MAX(0., PRCT(:)-PHLC_LRC(:))
    ELSEWHERE (PRCT(:).GT.0. .AND. PCF(:).GT.0. .AND. ZHLC_RCMAX(:).LE.ZRCRAUTC(:))
      ! Put all available cloud water and his fraction in the low part
      PHLC_LCF(:)=PCF(:)
      PHLC_HCF(:)=0.
      PHLC_LRC(:)=PRCT(:)
      PHLC_HRC(:)=0.
    ELSEWHERE
      PHLC_LCF(:)=0.
      PHLC_HCF(:)=0.
      PHLC_LRC(:)=0.
      PHLC_HRC(:)=0.
    END WHERE
     PRF(:)=PHLC_HCF(:) !Precipitation fraction
   ELSE
    !wrong HSUBG_PR_PDF case
    CALL ABORT
    STOP 'wrong HSUBG_PR_PDF case'
   ENDIF
ELSE
  !wrong HSUBG_AUCV case
  CALL ABORT
  STOP 'wrong HSUBG_AUCV case'
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ICE4_COMPUTE_PDF', 1, ZHOOK_HANDLE)
END SUBROUTINE ICE4_COMPUTE_PDF
