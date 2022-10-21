!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_COMPUTE_PDF
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_COMPUTE_PDF(CST, ICEP, ICED, KSIZE, HSUBG_AUCV_RC, HSUBG_AUCV_RI, HSUBG_PR_PDF, &
                            PRHODREF, PRCT, PRIT, PCF, PT, PSIGMA_RC,&
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
!!
!
!
!*      0. DECLARATIONS
!          ------------
!
!
USE MODD_CST,            ONLY: CST_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
USE MODE_MSG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
INTEGER,                INTENT(IN)  :: KSIZE
CHARACTER(LEN=4),       INTENT(IN)  :: HSUBG_AUCV_RC     ! Kind of Subgrid autoconversion method for cloud water
CHARACTER(LEN=80),      INTENT(IN)  :: HSUBG_AUCV_RI     ! Kind of Subgrid autoconversion method for cloud ice
CHARACTER(LEN=80),      INTENT(IN)  :: HSUBG_PR_PDF   ! pdf for subgrid precipitation
REAL, DIMENSION(KSIZE), INTENT(IN)  :: PRHODREF   ! Reference density
REAL, DIMENSION(KSIZE), INTENT(IN)  :: PRCT       ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE), INTENT(IN)  :: PRIT       ! Ice Crystal m.r. at t
REAL, DIMENSION(KSIZE), INTENT(IN)  :: PCF        ! Cloud fraction
REAL, DIMENSION(KSIZE), INTENT(IN)  :: PT         ! Temperature
REAL, DIMENSION(KSIZE), INTENT(IN)  :: PSIGMA_RC  ! Standard deviation of rc at time t
!Note for INTENT STATUS: in 'ADJU' case the PHL?_??? variables must be able to "cross" the subroutine untouched
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PHLC_HCF   ! HLCLOUDS : fraction of High Cloud Fraction in grid
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PHLC_LCF   ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                                  !    note that PCF = PHLC_HCF + PHLC_LCF
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PHLC_HRC   ! HLCLOUDS : LWC that is High LWC in grid
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PHLC_LRC   ! HLCLOUDS : LWC that is Low  LWC in grid
                                                  !    note that PRC = PHLC_HRC + PHLC_LRC
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PHLI_HCF
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PHLI_LCF
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PHLI_HRI
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PHLI_LRI
REAL, DIMENSION(KSIZE), INTENT(OUT) :: PRF        ! Rain fraction
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KSIZE) :: ZRCRAUTC,      & !RC value to begin rain formation =XCRIAUTC/RHODREF
                          ZCRIAUTI,      & !RI value to begin snow formation
                          ZHLC_RCMAX,    & !HLCLOUDS : maximum value for RC in distribution
                          ZHLC_LRCLOCAL, & !HLCLOUDS : LWC that is Low  LWC local in LCF
                          ZHLC_HRCLOCAL, & !HLCLOUDS : LWC that is High LWC local in HCF
                                                    !    note that ZRC/CF = ZHLC_HRCLOCAL+ ZHLC_LRCLOCAL
                                                    !                     = PHLC_HRC/HCF+ PHLC_LRC/LCF
                          ZSUMRC, ZSUMRI
REAL :: ZCOEFFRCM
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JI
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_COMPUTE_PDF', 0, ZHOOK_HANDLE)!

!Cloud water split between high and low content part is done according to autoconversion option
ZRCRAUTC(:)=ICEP%XCRIAUTC/PRHODREF(:) ! Autoconversion rc threshold
IF(HSUBG_AUCV_RC=='NONE') THEN
  !Cloud water is entirely in low or high part
  DO JI=1,KSIZE 
    IF(PRCT(JI)>ZRCRAUTC(JI))THEN
      PHLC_HCF(JI)=1.
      PHLC_LCF(JI)=0.
      PHLC_HRC(JI)=PRCT(JI)
      PHLC_LRC(JI)=0.
    ELSEIF(PRCT(JI)>ICED%XRTMIN(2))THEN
      PHLC_HCF(JI)=0.
      PHLC_LCF(JI)=1.
      PHLC_HRC(JI)=0.
      PHLC_LRC(JI)=PRCT(JI)
    ELSE
      PHLC_HCF(JI)=0.
      PHLC_LCF(JI)=0.
      PHLC_HRC(JI)=0.
      PHLC_LRC(JI)=0.
    ENDIF
  ENDDO

ELSEIF(HSUBG_AUCV_RC=='CLFR') THEN
  !Cloud water is only in the cloudy part and entirely in low or high part
  DO JI=1,KSIZE 
    IF(PCF(JI)>0. .AND. PRCT(JI)>ZRCRAUTC(JI)*PCF(JI))THEN
      PHLC_HCF(JI)=PCF(JI)
      PHLC_LCF(JI)=0.
      PHLC_HRC(JI)=PRCT(JI)
      PHLC_LRC(JI)=0.
    ELSEIF(PCF(JI)>0. .AND. PRCT(JI)>ICED%XRTMIN(2))THEN
      PHLC_HCF(JI)=0.
      PHLC_LCF(JI)=PCF(JI)
      PHLC_HRC(JI)=0.0
      PHLC_LRC(JI)=PRCT(JI)
    ELSE
      PHLC_HCF(JI)=0.
      PHLC_LCF(JI)=0.
      PHLC_HRC(JI)=0.
      PHLC_LRC(JI)=0.
    ENDIF
  ENDDO
ELSEIF(HSUBG_AUCV_RC=='ADJU') THEN
  DO JI=1,KSIZE 
    ZSUMRC(JI)=PHLC_LRC(JI)+PHLC_HRC(JI)
    IF(ZSUMRC(JI) .GT. 0.)THEN
      PHLC_LRC(JI)=PHLC_LRC(JI)*PRCT(JI)/ZSUMRC(JI)
      PHLC_HRC(JI)=PHLC_HRC(JI)*PRCT(JI)/ZSUMRC(JI)
    ELSE
      PHLC_LRC(JI)=0.
      PHLC_HRC(JI)=0.
    ENDIF
  ENDDO
ELSEIF(HSUBG_AUCV_RC=='PDF ') THEN
  !Cloud water is split between high and low part according to a PDF
  !    'HLCRECTPDF'    : rectangular PDF form
  !    'HLCTRIANGPDF'  : triangular PDF form
  !    'HLCQUADRAPDF'  : second order quadratic PDF form
  !    'HLCISOTRIPDF'  : isocele triangular PDF
  !    'SIGM'          : Redelsperger and Sommeria (1986)
  IF(HSUBG_PR_PDF=='SIGM') THEN
    ! Redelsperger and Sommeria (1986) but organised according to Turner (2011, 2012)
    DO JI=1,KSIZE 
      IF (PRCT(JI)>ZRCRAUTC(JI)+PSIGMA_RC(JI))THEN
        PHLC_HCF(JI)=1.
        PHLC_LCF(JI)=0.
        PHLC_HRC(JI)=PRCT(JI)
        PHLC_LRC(JI)=0.
      ELSEIF(PRCT(JI)> (ZRCRAUTC(JI)-PSIGMA_RC(JI)) .AND. PRCT(JI)<=(ZRCRAUTC(JI)+PSIGMA_RC(JI))       )THEN
        PHLC_HCF(JI)=(PRCT(JI)+PSIGMA_RC(JI)-ZRCRAUTC(JI))/ &
        &(2.*PSIGMA_RC(JI))
        PHLC_LCF(JI)=MAX(0., PCF(JI)-PHLC_HCF(JI))
        PHLC_HRC(JI)=(PRCT(JI)+PSIGMA_RC(JI)-ZRCRAUTC(JI))* &
        &(PRCT(JI)+PSIGMA_RC(JI)+ZRCRAUTC(JI))/ &
        &(4.*PSIGMA_RC(JI))
        PHLC_LRC(JI)=MAX(0., PRCT(JI)-PHLC_HRC(JI))
      ELSEIF(PRCT(JI)>ICED%XRTMIN(2) .AND. PCF(JI)>0.)THEN
        PHLC_HCF(JI)=0.
        PHLC_LCF(JI)=PCF(JI)
        PHLC_HRC(JI)=0.
        PHLC_LRC(JI)=PRCT(JI)
      ELSE
        PHLC_HCF(JI)=0.
        PHLC_LCF(JI)=0.
        PHLC_HRC(JI)=0.
        PHLC_LRC(JI)=0.
      ENDIF
    ENDDO
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
    DO JI=1,KSIZE 
      IF(PRCT(JI).GT.0. .AND. PCF(JI).GT.0.)THEN
        ZHLC_RCMAX(JI)=ZCOEFFRCM*PRCT(JI)/PCF(JI)
      ENDIF
    ! Split available water and cloud fraction in two parts
    ! Calculate local mean values int he low and high parts for the 3 PDF forms:
      IF(HSUBG_PR_PDF=='HLCRECTPDF') THEN
        IF(PRCT(JI).GT.0. .AND. PCF(JI).GT.0. .AND. ZHLC_RCMAX(JI).GT.ZRCRAUTC(JI))THEN
          ZHLC_LRCLOCAL(JI)=0.5*ZRCRAUTC(JI)
          ZHLC_HRCLOCAL(JI)=( ZHLC_RCMAX(JI) + ZRCRAUTC(JI))/2.0
        ENDIF
      ELSE IF(HSUBG_PR_PDF=='HLCTRIANGPDF') THEN
        IF(PRCT(JI).GT.0. .AND. PCF(JI).GT.0. .AND. ZHLC_RCMAX(JI).GT.ZRCRAUTC(JI))THEN
          ZHLC_LRCLOCAL(JI)=( ZRCRAUTC(JI) *(3.0 * ZHLC_RCMAX(JI) - 2.0 * ZRCRAUTC(JI) ) ) &
          / (3.0 * (2.0 * ZHLC_RCMAX(JI) - ZRCRAUTC(JI)  ) )
          ZHLC_HRCLOCAL(JI)=(ZHLC_RCMAX(JI) + 2.0*ZRCRAUTC(JI)) / 3.0
        ENDIF
      ELSE IF(HSUBG_PR_PDF=='HLCQUADRAPDF') THEN
        IF(PRCT(JI).GT.0. .AND. PCF(JI).GT.0. .AND. ZHLC_RCMAX(JI).GT.ZRCRAUTC(JI))THEN
          ZHLC_LRCLOCAL(JI)=(3.0 *ZRCRAUTC(JI)**3 - 8.0 *ZRCRAUTC(JI)**2 * ZHLC_RCMAX(JI) &
          + 6.0*ZRCRAUTC(JI) *ZHLC_RCMAX(JI)**2 ) &
          / &
          (4.0* ZRCRAUTC(JI)**2 -12.0*ZRCRAUTC(JI) *ZHLC_RCMAX(JI) &
          + 12.0 * ZHLC_RCMAX(JI)**2 )
          ZHLC_HRCLOCAL(JI)=(ZHLC_RCMAX(JI) + 3.0*ZRCRAUTC(JI))/4.0
        ENDIF
      ELSE IF(HSUBG_PR_PDF=='HLCISOTRIPDF') THEN
        IF(PRCT(JI).GT.0. .AND. PCF(JI).GT.0. .AND. ZHLC_RCMAX(JI).GT.ZRCRAUTC(JI))THEN
          IF((PRCT(JI) / PCF(JI)).LE.ZRCRAUTC(JI))THEN
            ZHLC_LRCLOCAL(JI)=( (ZHLC_RCMAX(JI))**3 &
            -(12.0 * (ZHLC_RCMAX(JI))*(ZRCRAUTC(JI))**2) &
            +(8.0 * ZRCRAUTC(JI)**3) ) &
            /( (6.0 * (ZHLC_RCMAX(JI))**2) &
            -(24.0 * (ZHLC_RCMAX(JI)) * ZRCRAUTC(JI)) &
            +(12.0 * ZRCRAUTC(JI)**2) )
            ZHLC_HRCLOCAL(JI)=( ZHLC_RCMAX(JI) + 2.0 * ZRCRAUTC(JI) )/3.0
          ELSE
            ZHLC_LRCLOCAL(JI)=(2.0/3.0) * ZRCRAUTC(JI)
            ZHLC_HRCLOCAL(JI)=(3.0*ZHLC_RCMAX(JI)**3 - 8.0*ZRCRAUTC(JI)**3) &
            / (6.0 * ZHLC_RCMAX(JI)**2 - 12.0*ZRCRAUTC(JI)**2)
          ENDIF
        ENDIF
      END IF
    ! Compare r_cM  to r_cR to know if cloud water content is high enough to split in two parts or not
      IF (PRCT(JI).GT.0. .AND. PCF(JI).GT.0. .AND. ZHLC_RCMAX(JI).GT.ZRCRAUTC(JI))THEN
      ! Calculate final values for LCF and HCF:
        PHLC_LCF(JI)=PCF(JI) &
        *(ZHLC_HRCLOCAL(JI)- &
        (PRCT(JI) / PCF(JI))) &
        / (ZHLC_HRCLOCAL(JI)-ZHLC_LRCLOCAL(JI))
        PHLC_HCF(JI)=MAX(0., PCF(JI)-PHLC_LCF(JI))
      !
      ! Calculate final values for LRC and HRC:
        PHLC_LRC(JI)=ZHLC_LRCLOCAL(JI)*PHLC_LCF(JI)
        PHLC_HRC(JI)=MAX(0., PRCT(JI)-PHLC_LRC(JI))
      ELSEIF (PRCT(JI).GT.0. .AND. PCF(JI).GT.0. .AND. ZHLC_RCMAX(JI).LE.ZRCRAUTC(JI))THEN
      ! Put all available cloud water and his fraction in the low part
        PHLC_LCF(JI)=PCF(JI)
        PHLC_HCF(JI)=0.
        PHLC_LRC(JI)=PRCT(JI)
        PHLC_HRC(JI)=0.
      ELSE
        PHLC_LCF(JI)=0.
        PHLC_HCF(JI)=0.
        PHLC_LRC(JI)=0.
        PHLC_HRC(JI)=0.
      ENDIF
    ENDDO
  ELSE
    CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_COMPUTE_PDF','wrong HSUBG_PR_PDF case')
  ENDIF
ELSE
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_COMPUTE_PDF','wrong HSUBG_AUCV case')
ENDIF
!
!Ice water split between high and low content part is done according to autoconversion option
DO JI=1,KSIZE 
  ZCRIAUTI(JI)=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(PT(JI)-CST%XTT)+ICEP%XBCRIAUTI)) ! Autoconversion ri threshold
ENDDO
IF(HSUBG_AUCV_RI=='NONE') THEN
  DO JI=1,KSIZE 
  !Cloud water is entirely in low or high part
    IF(PRIT(JI)>ZCRIAUTI(JI))THEN
      PHLI_HCF(JI)=1.
      PHLI_LCF(JI)=0.
      PHLI_HRI(JI)=PRIT(JI)
      PHLI_LRI(JI)=0.
    ELSEIF(PRIT(JI)>ICED%XRTMIN(4))THEN
      PHLI_HCF(JI)=0.
      PHLI_LCF(JI)=1.
      PHLI_HRI(JI)=0.
      PHLI_LRI(JI)=PRIT(JI)
    ELSE
      PHLI_HCF(JI)=0.
      PHLI_LCF(JI)=0.
      PHLI_HRI(JI)=0.
      PHLI_LRI(JI)=0.
    ENDIF
  ENDDO
ELSEIF(HSUBG_AUCV_RI=='CLFR') THEN
  !Cloud water is only in the cloudy part and entirely in low or high part
  DO JI=1,KSIZE 
    IF(PCF(JI)>0. .AND. PRIT(JI)>ZCRIAUTI(JI)*PCF(JI))THEN
      PHLI_HCF(JI)=PCF(JI)
      PHLI_LCF(JI)=0.
      PHLI_HRI(JI)=PRIT(JI)
      PHLI_LRI(JI)=0.
    ELSEIF(PCF(JI)>0. .AND. PRIT(JI)>ICED%XRTMIN(4))THEN
      PHLI_HCF(JI)=0.
      PHLI_LCF(JI)=PCF(JI)
      PHLI_HRI(JI)=0.0
      PHLI_LRI(JI)=PRIT(JI)
    ELSE
      PHLI_HCF(JI)=0.
      PHLI_LCF(JI)=0.
      PHLI_HRI(JI)=0.
      PHLI_LRI(JI)=0.
    ENDIF
  ENDDO
ELSEIF(HSUBG_AUCV_RI=='ADJU') THEN
  DO JI=1,KSIZE 
    ZSUMRI(JI)=PHLI_LRI(JI)+PHLI_HRI(JI)
    IF(ZSUMRI(JI) .GT. 0.)THEN
      PHLI_LRI(JI)=PHLI_LRI(JI)*PRIT(JI)/ZSUMRI(JI)
      PHLI_HRI(JI)=PHLI_HRI(JI)*PRIT(JI)/ZSUMRI(JI)
    ELSE
      PHLI_LRI(JI)=0.
      PHLI_HRI(JI)=0.
    ENDIF
  ENDDO
ELSE
  !wrong HSUBG_AUCV_RI case
  CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'ICE4_COMPUTE_PDF', 'wrong HSUBG_AUCV_RI case' )
ENDIF
!
DO JI=1,KSIZE 
#ifdef REPRO48
  PRF(JI)=PHLC_HCF(JI)
#else
  PRF(JI)=MAX(PHLC_HCF(JI),PHLI_HCF(JI))
#endif
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ICE4_COMPUTE_PDF', 1, ZHOOK_HANDLE)
END SUBROUTINE ICE4_COMPUTE_PDF
END MODULE MODE_ICE4_COMPUTE_PDF
