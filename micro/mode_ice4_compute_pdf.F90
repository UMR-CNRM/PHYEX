!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_COMPUTE_PDF
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_COMPUTE_PDF(CST, ICEP, ICED, D, HSUBG_AUCV_RC, HSUBG_AUCV_RI, HSUBG_PR_PDF, &
                            LDMICRO, PRHODREF, PRCT, PRIT, PCF, PT, PSIGMA_RC,&
                            PHLC_HCF, PHLC_LCF, PHLC_HRC, PHLC_LRC, &
                            PHLI_HCF, PHLI_LCF, PHLI_HRI, PHLI_LRI, PRF)

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
!!      S. Riette Sept 23: LDMICRO mask
!
!
!*      0. DECLARATIONS
!          ------------
!
!
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_CST,            ONLY: CST_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
USE MODE_MSG, ONLY: NVERB_FATAL, PRINT_MSG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
CHARACTER(LEN=4),       INTENT(IN)  :: HSUBG_AUCV_RC     ! Kind of Subgrid autoconversion method for cloud water
CHARACTER(LEN=80),      INTENT(IN)  :: HSUBG_AUCV_RI     ! Kind of Subgrid autoconversion method for cloud ice
CHARACTER(LEN=80),      INTENT(IN)  :: HSUBG_PR_PDF   ! pdf for subgrid precipitation
LOGICAL, DIMENSION(D%NIJT), INTENT(IN)  :: LDMICRO    ! Computation mask
REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PRHODREF   ! Reference density
REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PRCT       ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PRIT       ! Ice Crystal m.r. at t
REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PCF        ! Cloud fraction
REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PT         ! Temperature
REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PSIGMA_RC  ! Standard deviation of rc at time t
!Note for INTENT STATUS: in 'ADJU' case the PHL?_??? variables must be able to "cross" the subroutine untouched
REAL, DIMENSION(D%NIJT), INTENT(INOUT) :: PHLC_HCF   ! HLCLOUDS : fraction of High Cloud Fraction in grid
REAL, DIMENSION(D%NIJT), INTENT(INOUT) :: PHLC_LCF   ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                                  !    note that PCF = PHLC_HCF + PHLC_LCF
REAL, DIMENSION(D%NIJT), INTENT(INOUT) :: PHLC_HRC   ! HLCLOUDS : LWC that is High LWC in grid
REAL, DIMENSION(D%NIJT), INTENT(INOUT) :: PHLC_LRC   ! HLCLOUDS : LWC that is Low  LWC in grid
                                                  !    note that PRC = PHLC_HRC + PHLC_LRC
REAL, DIMENSION(D%NIJT), INTENT(INOUT) :: PHLI_HCF
REAL, DIMENSION(D%NIJT), INTENT(INOUT) :: PHLI_LCF
REAL, DIMENSION(D%NIJT), INTENT(INOUT) :: PHLI_HRI
REAL, DIMENSION(D%NIJT), INTENT(INOUT) :: PHLI_LRI
REAL, DIMENSION(D%NIJT), INTENT(OUT) :: PRF        ! Rain fraction
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT) :: ZRCRAUTC,      & !RC value to begin rain formation =XCRIAUTC/RHODREF
                          ZCRIAUTI,      & !RI value to begin snow formation
                          ZHLC_RCMAX,    & !HLCLOUDS : maximum value for RC in distribution
                          ZHLC_LRCLOCAL, & !HLCLOUDS : LWC that is Low  LWC local in LCF
                          ZHLC_HRCLOCAL, & !HLCLOUDS : LWC that is High LWC local in HCF
                                                    !    note that ZRC/CF = ZHLC_HRCLOCAL+ ZHLC_LRCLOCAL
                                                    !                     = PHLC_HRC/HCF+ PHLC_LRC/LCF
                          ZSUMRC, ZSUMRI
REAL :: ZCOEFFRCM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JIJ
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_COMPUTE_PDF', 0, ZHOOK_HANDLE)
!
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
#endif
!
!Cloud water split between high and low content part is done according to autoconversion option

DO JIJ=D%NIJB, D%NIJE
  IF (LDMICRO(JIJ)) THEN
    ZRCRAUTC(JIJ)=ICEP%XCRIAUTC/PRHODREF(JIJ) ! Autoconversion rc threshold
  ELSE
    ZRCRAUTC(JIJ)=0.
  END IF
END DO

IF(HSUBG_AUCV_RC=='NONE') THEN
  !Cloud water is entirely in low or high part

 DO JIJ=D%NIJB, D%NIJE
    IF (.NOT. LDMICRO(JIJ)) THEN
      ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
      PHLC_HCF(JIJ)=0.
      PHLC_LCF(JIJ)=0.
      PHLC_HRC(JIJ)=0.
      PHLC_LRC(JIJ)=0.
    ELSE IF (PRCT(JIJ)>ZRCRAUTC(JIJ) .AND. LDMICRO(JIJ)) THEN
      PHLC_HCF(JIJ)=1.
      PHLC_LCF(JIJ)=0.
      PHLC_HRC(JIJ)=PRCT(JIJ)
      PHLC_LRC(JIJ)=0.
    ELSE IF (PRCT(JIJ)>ICED%XRTMIN(2) .AND. LDMICRO(JIJ)) THEN
      PHLC_HCF(JIJ)=0.
      PHLC_LCF(JIJ)=1.
      PHLC_HRC(JIJ)=0.
      PHLC_LRC(JIJ)=PRCT(JIJ)
    ELSE
      PHLC_HCF(JIJ)=0.
      PHLC_LCF(JIJ)=0.
      PHLC_HRC(JIJ)=0.
      PHLC_LRC(JIJ)=0.
    END IF
  END DO

ELSEIF(HSUBG_AUCV_RC=='CLFR') THEN
  !Cloud water is only in the cloudy part and entirely in low or high part

 DO JIJ=D%NIJB, D%NIJE
    IF (.NOT. LDMICRO(JIJ)) THEN
      ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
      PHLC_HCF(JIJ)=0.
      PHLC_LCF(JIJ)=0.
      PHLC_HRC(JIJ)=0.
      PHLC_LRC(JIJ)=0.
    ELSE IF (PCF(JIJ)>0. .AND. PRCT(JIJ)>ZRCRAUTC(JIJ)*PCF(JIJ) .AND. LDMICRO(JIJ)) THEN
      PHLC_HCF(JIJ)=PCF(JIJ)
      PHLC_LCF(JIJ)=0.
      PHLC_HRC(JIJ)=PRCT(JIJ)
      PHLC_LRC(JIJ)=0.
    ELSE IF (PCF(JIJ)>0. .AND. PRCT(JIJ)>ICED%XRTMIN(2) .AND. LDMICRO(JIJ)) THEN
      PHLC_HCF(JIJ)=0.
      PHLC_LCF(JIJ)=PCF(JIJ)
      PHLC_HRC(JIJ)=0.0
      PHLC_LRC(JIJ)=PRCT(JIJ)
    ELSE
      PHLC_HCF(JIJ)=0.
      PHLC_LCF(JIJ)=0.
      PHLC_HRC(JIJ)=0.
      PHLC_LRC(JIJ)=0.
    END IF
  END DO

ELSEIF(HSUBG_AUCV_RC=='ADJU') THEN

  DO JIJ=D%NIJB, D%NIJE
    IF (LDMICRO(JIJ)) THEN
      ZSUMRC(JIJ)=PHLC_LRC(JIJ)+PHLC_HRC(JIJ)
    ELSE
      ZSUMRC(JIJ)=0.
    END IF
    IF (ZSUMRC(JIJ) .GT. 1.E-20 .AND. LDMICRO(JIJ)) THEN
      PHLC_LRC(JIJ)=PHLC_LRC(JIJ)*PRCT(JIJ)/ZSUMRC(JIJ)
      PHLC_HRC(JIJ)=PHLC_HRC(JIJ)*PRCT(JIJ)/ZSUMRC(JIJ)
    ELSE
      PHLC_LRC(JIJ)=0.
      PHLC_HRC(JIJ)=0.
    END IF
  END DO

ELSEIF(HSUBG_AUCV_RC=='PDF') THEN
  !Cloud water is split between high and low part according to a PDF
  !    'HLCRECTPDF'    : rectangular PDF form
  !    'HLCTRIANGPDF'  : triangular PDF form
  !    'HLCQUADRAPDF'  : second order quadratic PDF form
  !    'HLCISOTRIPDF'  : isocele triangular PDF
  !    'SIGM'          : Redelsperger and Sommeria (1986)
  IF(HSUBG_PR_PDF=='SIGM') THEN
    ! Redelsperger and Sommeria (1986) but organised according to Turner (2011, 2012)

    DO JIJ=D%NIJB, D%NIJE
      IF (.NOT. LDMICRO(JIJ)) THEN
        ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
        PHLC_HCF(JIJ)=0.
        PHLC_LCF(JIJ)=0.
        PHLC_HRC(JIJ)=0.
        PHLC_LRC(JIJ)=0.
      ELSE IF (PRCT(JIJ)>ZRCRAUTC(JIJ)+PSIGMA_RC(JIJ) .AND. LDMICRO(JIJ)) THEN
        PHLC_HCF(JIJ)=1.
        PHLC_LCF(JIJ)=0.
        PHLC_HRC(JIJ)=PRCT(JIJ)
        PHLC_LRC(JIJ)=0.
      ELSE IF (PRCT(JIJ)> (ZRCRAUTC(JIJ)-PSIGMA_RC(JIJ)) .AND. PRCT(JIJ)<=(ZRCRAUTC(JIJ)+PSIGMA_RC(JIJ)) .AND. LDMICRO(JIJ)) THEN
        PHLC_HCF(JIJ)=(PRCT(JIJ)+PSIGMA_RC(JIJ)-ZRCRAUTC(JIJ))/ &
                    &(2.*PSIGMA_RC(JIJ))
        PHLC_LCF(JIJ)=MAX(0., PCF(JIJ)-PHLC_HCF(JIJ))
        PHLC_HRC(JIJ)=(PRCT(JIJ)+PSIGMA_RC(JIJ)-ZRCRAUTC(JIJ))* &
                    &(PRCT(JIJ)+PSIGMA_RC(JIJ)+ZRCRAUTC(JIJ))/ &
                    &(4.*PSIGMA_RC(JIJ))
        PHLC_LRC(JIJ)=MAX(0., PRCT(JIJ)-PHLC_HRC(JIJ))
      ELSE IF (PRCT(JIJ)>ICED%XRTMIN(2) .AND. PCF(JIJ)>0. .AND. LDMICRO(JIJ)) THEN
        PHLC_HCF(JIJ)=0.
        PHLC_LCF(JIJ)=PCF(JIJ)
        PHLC_HRC(JIJ)=0.
        PHLC_LRC(JIJ)=PRCT(JIJ)
      ELSE
        PHLC_HCF(JIJ)=0.
        PHLC_LCF(JIJ)=0.
        PHLC_HRC(JIJ)=0.
        PHLC_LRC(JIJ)=0.
      END IF
    END DO

  ELSEIF(HSUBG_PR_PDF=='HLCRECTPDF' .OR. HSUBG_PR_PDF=='HLCISOTRIPDF' .OR. &
         &HSUBG_PR_PDF=='HLCTRIANGPDF' .OR. HSUBG_PR_PDF=='HLCQUADRAPDF') THEN
    ! Turner (2011, 2012)
    ! Calculate maximum value r_cM from PDF forms
    IF(HSUBG_PR_PDF=='HLCRECTPDF' .OR. HSUBG_PR_PDF=='HLCISOTRIPDF') THEN
      ZCOEFFRCM=2.
    ELSE IF(HSUBG_PR_PDF=='HLCTRIANGPDF') THEN
      ZCOEFFRCM=3.
    ELSE IF(HSUBG_PR_PDF=='HLCQUADRAPDF') THEN
      ZCOEFFRCM=4.
    END IF

    DO JIJ=D%NIJB, D%NIJE
      IF (.NOT. LDMICRO(JIJ)) THEN
        ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
        ZHLC_RCMAX(JIJ)=0.
      ELSE IF (PRCT(JIJ).GT.0. .AND. PCF(JIJ).GT.0. .AND. LDMICRO(JIJ)) THEN
        ZHLC_RCMAX(JIJ)=ZCOEFFRCM*PRCT(JIJ)/PCF(JIJ)
      ELSE
        ZHLC_RCMAX(JIJ)=0.
      END IF
    END DO

    ! Split available water and cloud fraction in two parts
    ! Calculate local mean values int he low and high parts for the 3 PDF forms:
    IF(HSUBG_PR_PDF=='HLCRECTPDF') THEN
      DO JIJ=D%NIJB, D%NIJE
        IF (.NOT. LDMICRO(JIJ)) THEN
          ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
          ZHLC_LRCLOCAL(JIJ)=0.
          ZHLC_HRCLOCAL(JIJ)=0
        ELSE IF (PRCT(JIJ).GT.0. .AND. PCF(JIJ).GT.0. .AND. ZHLC_RCMAX(JIJ).GT.ZRCRAUTC(JIJ) .AND. LDMICRO(JIJ)) THEN
          ZHLC_LRCLOCAL(JIJ)=0.5*ZRCRAUTC(JIJ)
          ZHLC_HRCLOCAL(JIJ)=( ZHLC_RCMAX(JIJ) + ZRCRAUTC(JIJ))/2.0
        ELSE
          ZHLC_LRCLOCAL(JIJ)=0.
          ZHLC_HRCLOCAL(JIJ)=0.
        END IF
      END DO
    ELSE IF(HSUBG_PR_PDF=='HLCTRIANGPDF') THEN
      DO JIJ=D%NIJB, D%NIJE
        IF (.NOT. LDMICRO(JIJ)) THEN
          ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
          ZHLC_LRCLOCAL(JIJ)=0.
          ZHLC_HRCLOCAL(JIJ)=0
        ELSE IF (PRCT(JIJ).GT.0. .AND. PCF(JIJ).GT.0. .AND. ZHLC_RCMAX(JIJ).GT.ZRCRAUTC(JIJ) .AND. LDMICRO(JIJ)) THEN
          ZHLC_LRCLOCAL(JIJ)=( ZRCRAUTC(JIJ) *(3.0 * ZHLC_RCMAX(JIJ) - 2.0 * ZRCRAUTC(JIJ) ) ) &
                          / (3.0 * (2.0 * ZHLC_RCMAX(JIJ) - ZRCRAUTC(JIJ)  ) )
          ZHLC_HRCLOCAL(JIJ)=(ZHLC_RCMAX(JIJ) + 2.0*ZRCRAUTC(JIJ)) / 3.0
        ELSE
          ZHLC_LRCLOCAL(JIJ)=0.
          ZHLC_HRCLOCAL(JIJ)=0.
        END IF
      END DO
    ELSE IF(HSUBG_PR_PDF=='HLCQUADRAPDF') THEN
      DO JIJ=D%NIJB, D%NIJE
        IF (.NOT. LDMICRO(JIJ)) THEN
          ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
          ZHLC_LRCLOCAL(JIJ)=0.
          ZHLC_HRCLOCAL(JIJ)=0
        ELSE IF (PRCT(JIJ).GT.0. .AND. PCF(JIJ).GT.0. .AND. ZHLC_RCMAX(JIJ).GT.ZRCRAUTC(JIJ) .AND. LDMICRO(JIJ)) THEN
          ZHLC_LRCLOCAL(JIJ)=(3.0 *ZRCRAUTC(JIJ)**3 - 8.0 *ZRCRAUTC(JIJ)**2 * ZHLC_RCMAX(JIJ) &
                          + 6.0*ZRCRAUTC(JIJ) *ZHLC_RCMAX(JIJ)**2 ) &
                          / &
                          (4.0* ZRCRAUTC(JIJ)**2 -12.0*ZRCRAUTC(JIJ) *ZHLC_RCMAX(JIJ) &
                          + 12.0 * ZHLC_RCMAX(JIJ)**2 )
          ZHLC_HRCLOCAL(JIJ)=(ZHLC_RCMAX(JIJ) + 3.0*ZRCRAUTC(JIJ))/4.0
        ELSE
          ZHLC_LRCLOCAL(JIJ)=0.
          ZHLC_HRCLOCAL(JIJ)=0.
        END IF
      END DO
    ELSE IF(HSUBG_PR_PDF=='HLCISOTRIPDF') THEN
      DO JIJ=D%NIJB, D%NIJE
        IF (.NOT. LDMICRO(JIJ)) THEN
          ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
          ZHLC_LRCLOCAL(JIJ)=0.
          ZHLC_HRCLOCAL(JIJ)=0
        ELSE IF (PRCT(JIJ).LE.ZRCRAUTC(JIJ)*PCF(JIJ) .AND. &
              &PRCT(JIJ).GT.0. .AND. PCF(JIJ).GT.0. .AND. &
              &ZHLC_RCMAX(JIJ).GT.ZRCRAUTC(JIJ) .AND. LDMICRO(JIJ)) THEN
          ZHLC_LRCLOCAL(JIJ)=( (ZHLC_RCMAX(JIJ))**3 &
                          -(12.0 * (ZHLC_RCMAX(JIJ))*(ZRCRAUTC(JIJ))**2) &
                          +(8.0 * ZRCRAUTC(JIJ)**3) ) &
                          /( (6.0 * (ZHLC_RCMAX(JIJ))**2) &
                          -(24.0 * (ZHLC_RCMAX(JIJ)) * ZRCRAUTC(JIJ)) &
                          +(12.0 * ZRCRAUTC(JIJ)**2) )
          ZHLC_HRCLOCAL(JIJ)=( ZHLC_RCMAX(JIJ) + 2.0 * ZRCRAUTC(JIJ) )/3.0
        ELSE IF (PRCT(JIJ).GT.0. .AND. PCF(JIJ).GT.0. .AND. ZHLC_RCMAX(JIJ).GT.ZRCRAUTC(JIJ) .AND. LDMICRO(JIJ)) THEN
          ZHLC_LRCLOCAL(JIJ)=(2.0/3.0) * ZRCRAUTC(JIJ)
          ZHLC_HRCLOCAL(JIJ)=(3.0*ZHLC_RCMAX(JIJ)**3 - 8.0*ZRCRAUTC(JIJ)**3) &
                          / (6.0 * ZHLC_RCMAX(JIJ)**2 - 12.0*ZRCRAUTC(JIJ)**2)
        ELSE
          ZHLC_LRCLOCAL(JIJ)=0.
          ZHLC_HRCLOCAL(JIJ)=0.
        END IF
      END DO
    END IF
    ! Compare r_cM  to r_cR to know if cloud water content is high enough to split in two parts or not
    DO JIJ=D%NIJB, D%NIJE
      IF (.NOT. LDMICRO(JIJ)) THEN
        ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
        PHLC_HCF(JIJ)=0.
        PHLC_LCF(JIJ)=0.
        PHLC_HRC(JIJ)=0.
        PHLC_LRC(JIJ)=0.
      ELSE IF (PRCT(JIJ).GT.0. .AND. PCF(JIJ).GT.0. .AND. ZHLC_RCMAX(JIJ).GT.ZRCRAUTC(JIJ) .AND. LDMICRO(JIJ)) THEN
        ! Calculate final values for LCF and HCF:
        PHLC_LCF(JIJ)=PCF(JIJ) &
                      *(ZHLC_HRCLOCAL(JIJ)- &
                      (PRCT(JIJ) / PCF(JIJ))) &
                      / (ZHLC_HRCLOCAL(JIJ)-ZHLC_LRCLOCAL(JIJ))
        PHLC_HCF(JIJ)=MAX(0., PCF(JIJ)-PHLC_LCF(JIJ))
        !
        ! Calculate final values for LRC and HRC:
        PHLC_LRC(JIJ)=ZHLC_LRCLOCAL(JIJ)*PHLC_LCF(JIJ)
        PHLC_HRC(JIJ)=MAX(0., PRCT(JIJ)-PHLC_LRC(JIJ))
      ELSE IF (PRCT(JIJ).GT.0. .AND. PCF(JIJ).GT.0. .AND. ZHLC_RCMAX(JIJ).LE.ZRCRAUTC(JIJ) .AND. LDMICRO(JIJ)) THEN
        ! Put all available cloud water and his fraction in the low part
        PHLC_LCF(JIJ)=PCF(JIJ)
        PHLC_HCF(JIJ)=0.
        PHLC_LRC(JIJ)=PRCT(JIJ)
        PHLC_HRC(JIJ)=0.
      ELSE
        PHLC_LCF(JIJ)=0.
        PHLC_HCF(JIJ)=0.
        PHLC_LRC(JIJ)=0.
        PHLC_HRC(JIJ)=0.
      END IF
    END DO

  ELSE
    CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_COMPUTE_PDF','wrong HSUBG_PR_PDF case')
  ENDIF
ELSE
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_COMPUTE_PDF','wrong HSUBG_AUCV case')
ENDIF
!
!Ice water split between high and low content part is done according to autoconversion option

DO JIJ=D%NIJB, D%NIJE
  IF (LDMICRO(JIJ)) THEN
    ZCRIAUTI(JIJ)=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(PT(JIJ)-CST%XTT)+ICEP%XBCRIAUTI)) ! Autoconversion ri threshold
  ELSE
    ZCRIAUTI(JIJ)=0.
  END IF
END DO

IF(HSUBG_AUCV_RI=='NONE') THEN
  !Cloud water is entirely in low or high part

  DO JIJ=D%NIJB, D%NIJE
    IF (.NOT. LDMICRO(JIJ)) THEN
      ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
      PHLI_HCF(JIJ)=0.
      PHLI_LCF(JIJ)=0.
      PHLI_HRI(JIJ)=0.
      PHLI_LRI(JIJ)=0.
    ELSE IF (PRIT(JIJ)>ZCRIAUTI(JIJ) .AND. LDMICRO(JIJ)) THEN
      PHLI_HCF(JIJ)=1.
      PHLI_LCF(JIJ)=0.
      PHLI_HRI(JIJ)=PRIT(JIJ)
      PHLI_LRI(JIJ)=0.
    ELSE IF (PRIT(JIJ)>ICED%XRTMIN(4) .AND. LDMICRO(JIJ)) THEN
      PHLI_HCF(JIJ)=0.
      PHLI_LCF(JIJ)=1.
      PHLI_HRI(JIJ)=0.
      PHLI_LRI(JIJ)=PRIT(JIJ)
    ELSE
      PHLI_HCF(JIJ)=0.
      PHLI_LCF(JIJ)=0.
      PHLI_HRI(JIJ)=0.
      PHLI_LRI(JIJ)=0.
    END IF
  END DO

ELSEIF(HSUBG_AUCV_RI=='CLFR') THEN
  !Cloud water is only in the cloudy part and entirely in low or high part

  DO JIJ=D%NIJB, D%NIJE
    IF (.NOT. LDMICRO(JIJ)) THEN
      ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
      PHLI_HCF(JIJ)=0.
      PHLI_LCF(JIJ)=0.
      PHLI_HRI(JIJ)=0.
      PHLI_LRI(JIJ)=0.
    ELSE IF (PCF(JIJ)>0. .AND. PRIT(JIJ)>ZCRIAUTI(JIJ)*PCF(JIJ) .AND. LDMICRO(JIJ)) THEN
      PHLI_HCF(JIJ)=PCF(JIJ)
      PHLI_LCF(JIJ)=0.
      PHLI_HRI(JIJ)=PRIT(JIJ)
      PHLI_LRI(JIJ)=0.
    ELSE IF (PCF(JIJ)>0. .AND. PRIT(JIJ)>ICED%XRTMIN(4) .AND. LDMICRO(JIJ)) THEN
      PHLI_HCF(JIJ)=0.
      PHLI_LCF(JIJ)=PCF(JIJ)
      PHLI_HRI(JIJ)=0.0
      PHLI_LRI(JIJ)=PRIT(JIJ)
    ELSE
      PHLI_HCF(JIJ)=0.
      PHLI_LCF(JIJ)=0.
      PHLI_HRI(JIJ)=0.
      PHLI_LRI(JIJ)=0.
    END IF
  END DO

ELSEIF(HSUBG_AUCV_RI=='ADJU') THEN

  DO JIJ=D%NIJB, D%NIJE
    IF (LDMICRO(JIJ)) THEN
      ZSUMRI(JIJ)=PHLI_LRI(JIJ)+PHLI_HRI(JIJ)
    ELSE
      ZSUMRI(JIJ)=0.
    END IF
    IF (ZSUMRI(JIJ) .GT. 1.E-20 .AND. LDMICRO(JIJ)) THEN
      PHLI_LRI(JIJ)=PHLI_LRI(JIJ)*PRIT(JIJ)/ZSUMRI(JIJ)
      PHLI_HRI(JIJ)=PHLI_HRI(JIJ)*PRIT(JIJ)/ZSUMRI(JIJ)
    ELSE
      PHLI_LRI(JIJ)=0.
      PHLI_HRI(JIJ)=0.
    END IF
  END DO

ELSE
  !wrong HSUBG_AUCV_RI case
  CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'ICE4_COMPUTE_PDF', 'wrong HSUBG_AUCV_RI case' )
ENDIF
!

DO JIJ=D%NIJB, D%NIJE
  IF (LDMICRO(JIJ)) THEN
    PRF(JIJ)=MAX(PHLC_HCF(JIJ),PHLI_HCF(JIJ))
  ELSE
    PRF(JIJ)=0.
  END IF
END DO

!
IF (LHOOK) CALL DR_HOOK('ICE4_COMPUTE_PDF', 1, ZHOOK_HANDLE)
END SUBROUTINE ICE4_COMPUTE_PDF

END MODULE MODE_ICE4_COMPUTE_PDF
