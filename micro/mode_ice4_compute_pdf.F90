!MNH_LIC Copyright 1994-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_COMPUTE_PDF
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_COMPUTE_PDF(CST, ICEP, ICED, KSIZE, HSUBG_AUCV_RC, HSUBG_AUCV_RI, HSUBG_PR_PDF, &
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
USE MODD_CST,            ONLY: CST_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
LOGICAL, DIMENSION(KSIZE), INTENT(IN)  :: LDMICRO    ! Computation mask
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JL
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_COMPUTE_PDF', 0, ZHOOK_HANDLE)!

!Cloud water split between high and low content part is done according to autoconversion option

DO JL=1, KSIZE
  IF (LDMICRO(JL)) THEN
    ZRCRAUTC(JL)=ICEP%XCRIAUTC/PRHODREF(JL) ! Autoconversion rc threshold
  ELSE
    ZRCRAUTC(JL)=0.
  END IF
END DO

IF(HSUBG_AUCV_RC=='NONE') THEN
  !Cloud water is entirely in low or high part

 DO JL=1, KSIZE
    IF (PRCT(JL)>ZRCRAUTC(JL) .AND. LDMICRO(JL)) THEN
      PHLC_HCF(JL)=1.
      PHLC_LCF(JL)=0.
      PHLC_HRC(JL)=PRCT(JL)
      PHLC_LRC(JL)=0.
    ELSE IF (PRCT(JL)>ICED%XRTMIN(2) .AND. LDMICRO(JL)) THEN
      PHLC_HCF(JL)=0.
      PHLC_LCF(JL)=1.
      PHLC_HRC(JL)=0.
      PHLC_LRC(JL)=PRCT(JL)
    ELSE
      PHLC_HCF(JL)=0.
      PHLC_LCF(JL)=0.
      PHLC_HRC(JL)=0.
      PHLC_LRC(JL)=0.
    END IF
  END DO

ELSEIF(HSUBG_AUCV_RC=='CLFR') THEN
  !Cloud water is only in the cloudy part and entirely in low or high part

 DO JL=1, KSIZE
    IF (PCF(JL)>0. .AND. PRCT(JL)>ZRCRAUTC(JL)*PCF(JL) .AND. LDMICRO(JL)) THEN
      PHLC_HCF(JL)=PCF(JL)
      PHLC_LCF(JL)=0.
      PHLC_HRC(JL)=PRCT(JL)
      PHLC_LRC(JL)=0.
    ELSE IF (PCF(JL)>0. .AND. PRCT(JL)>ICED%XRTMIN(2) .AND. LDMICRO(JL)) THEN
      PHLC_HCF(JL)=0.
      PHLC_LCF(JL)=PCF(JL)
      PHLC_HRC(JL)=0.0
      PHLC_LRC(JL)=PRCT(JL)
    ELSE
      PHLC_HCF(JL)=0.
      PHLC_LCF(JL)=0.
      PHLC_HRC(JL)=0.
      PHLC_LRC(JL)=0.
    END IF
  END DO

ELSEIF(HSUBG_AUCV_RC=='ADJU') THEN

  DO JL=1, KSIZE
    IF (LDMICRO(JL)) THEN
      ZSUMRC(JL)=PHLC_LRC(JL)+PHLC_HRC(JL)
    ELSE
      ZSUMRC(JL)=0.
    END IF
    IF (ZSUMRC(JL) .GT. 1.E-20 .AND. LDMICRO(JL)) THEN
      PHLC_LRC(JL)=PHLC_LRC(JL)*PRCT(JL)/ZSUMRC(JL)
      PHLC_HRC(JL)=PHLC_HRC(JL)*PRCT(JL)/ZSUMRC(JL)
    ELSE
      PHLC_LRC(JL)=0.
      PHLC_HRC(JL)=0.
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

    DO JL=1, KSIZE
      IF (PRCT(JL)>ZRCRAUTC(JL)+PSIGMA_RC(JL) .AND. LDMICRO(JL)) THEN
        PHLC_HCF(JL)=1.
        PHLC_LCF(JL)=0.
        PHLC_HRC(JL)=PRCT(JL)
        PHLC_LRC(JL)=0.
      ELSE IF (PRCT(JL)> (ZRCRAUTC(JL)-PSIGMA_RC(JL)) .AND. PRCT(JL)<=(ZRCRAUTC(JL)+PSIGMA_RC(JL)) .AND. LDMICRO(JL)) THEN
        PHLC_HCF(JL)=(PRCT(JL)+PSIGMA_RC(JL)-ZRCRAUTC(JL))/ &
                    &(2.*PSIGMA_RC(JL))
        PHLC_LCF(JL)=MAX(0., PCF(JL)-PHLC_HCF(JL))
        PHLC_HRC(JL)=(PRCT(JL)+PSIGMA_RC(JL)-ZRCRAUTC(JL))* &
                    &(PRCT(JL)+PSIGMA_RC(JL)+ZRCRAUTC(JL))/ &
                    &(4.*PSIGMA_RC(JL))
        PHLC_LRC(JL)=MAX(0., PRCT(JL)-PHLC_HRC(JL))
      ELSE IF (PRCT(JL)>ICED%XRTMIN(2) .AND. PCF(JL)>0. .AND. LDMICRO(JL)) THEN
        PHLC_HCF(JL)=0.
        PHLC_LCF(JL)=PCF(JL)
        PHLC_HRC(JL)=0.
        PHLC_LRC(JL)=PRCT(JL)
      ELSE
        PHLC_HCF(JL)=0.
        PHLC_LCF(JL)=0.
        PHLC_HRC(JL)=0.
        PHLC_LRC(JL)=0.
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

    DO JL=1, KSIZE
      IF (PRCT(JL).GT.0. .AND. PCF(JL).GT.0. .AND. LDMICRO(JL)) THEN
        ZHLC_RCMAX(JL)=ZCOEFFRCM*PRCT(JL)/PCF(JL)
      ELSE
        ZHLC_RCMAX(JL)=0.
      END IF
    END DO

    ! Split available water and cloud fraction in two parts
    ! Calculate local mean values int he low and high parts for the 3 PDF forms:
    IF(HSUBG_PR_PDF=='HLCRECTPDF') THEN
      DO JL=1, KSIZE
        IF (PRCT(JL).GT.0. .AND. PCF(JL).GT.0. .AND. ZHLC_RCMAX(JL).GT.ZRCRAUTC(JL) .AND. LDMICRO(JL)) THEN
          ZHLC_LRCLOCAL(JL)=0.5*ZRCRAUTC(JL)
          ZHLC_HRCLOCAL(JL)=( ZHLC_RCMAX(JL) + ZRCRAUTC(JL))/2.0
        ELSE
          ZHLC_LRCLOCAL(JL)=0.
          ZHLC_HRCLOCAL(JL)=0.
        END IF
      END DO
    ELSE IF(HSUBG_PR_PDF=='HLCTRIANGPDF') THEN
      DO JL=1, KSIZE
        IF (PRCT(JL).GT.0. .AND. PCF(JL).GT.0. .AND. ZHLC_RCMAX(JL).GT.ZRCRAUTC(JL) .AND. LDMICRO(JL)) THEN
          ZHLC_LRCLOCAL(JL)=( ZRCRAUTC(JL) *(3.0 * ZHLC_RCMAX(JL) - 2.0 * ZRCRAUTC(JL) ) ) &
                          / (3.0 * (2.0 * ZHLC_RCMAX(JL) - ZRCRAUTC(JL)  ) )
          ZHLC_HRCLOCAL(JL)=(ZHLC_RCMAX(JL) + 2.0*ZRCRAUTC(JL)) / 3.0
        ELSE
          ZHLC_LRCLOCAL(JL)=0.
          ZHLC_HRCLOCAL(JL)=0.
        END IF
      END DO
    ELSE IF(HSUBG_PR_PDF=='HLCQUADRAPDF') THEN
      DO JL=1, KSIZE
        IF (PRCT(JL).GT.0. .AND. PCF(JL).GT.0. .AND. ZHLC_RCMAX(JL).GT.ZRCRAUTC(JL) .AND. LDMICRO(JL)) THEN
          ZHLC_LRCLOCAL(JL)=(3.0 *ZRCRAUTC(JL)**3 - 8.0 *ZRCRAUTC(JL)**2 * ZHLC_RCMAX(JL) &
                          + 6.0*ZRCRAUTC(JL) *ZHLC_RCMAX(JL)**2 ) &
                          / &
                          (4.0* ZRCRAUTC(JL)**2 -12.0*ZRCRAUTC(JL) *ZHLC_RCMAX(JL) &
                          + 12.0 * ZHLC_RCMAX(JL)**2 )
          ZHLC_HRCLOCAL(JL)=(ZHLC_RCMAX(JL) + 3.0*ZRCRAUTC(JL))/4.0
        ELSE
          ZHLC_LRCLOCAL(JL)=0.
          ZHLC_HRCLOCAL(JL)=0.
        END IF
      END DO
    ELSE IF(HSUBG_PR_PDF=='HLCISOTRIPDF') THEN
      DO JL=1, KSIZE
        IF (PRCT(JL).LE.ZRCRAUTC(JL)*PCF(JL) .AND. &
              &PRCT(JL).GT.0. .AND. PCF(JL).GT.0. .AND. &
              &ZHLC_RCMAX(JL).GT.ZRCRAUTC(JL) .AND. LDMICRO(JL)) THEN
          ZHLC_LRCLOCAL(JL)=( (ZHLC_RCMAX(JL))**3 &
                          -(12.0 * (ZHLC_RCMAX(JL))*(ZRCRAUTC(JL))**2) &
                          +(8.0 * ZRCRAUTC(JL)**3) ) &
                          /( (6.0 * (ZHLC_RCMAX(JL))**2) &
                          -(24.0 * (ZHLC_RCMAX(JL)) * ZRCRAUTC(JL)) &
                          +(12.0 * ZRCRAUTC(JL)**2) )
          ZHLC_HRCLOCAL(JL)=( ZHLC_RCMAX(JL) + 2.0 * ZRCRAUTC(JL) )/3.0
        ELSE IF (PRCT(JL).GT.0. .AND. PCF(JL).GT.0. .AND. ZHLC_RCMAX(JL).GT.ZRCRAUTC(JL) .AND. LDMICRO(JL)) THEN
          ZHLC_LRCLOCAL(JL)=(2.0/3.0) * ZRCRAUTC(JL)
          ZHLC_HRCLOCAL(JL)=(3.0*ZHLC_RCMAX(JL)**3 - 8.0*ZRCRAUTC(JL)**3) &
                          / (6.0 * ZHLC_RCMAX(JL)**2 - 12.0*ZRCRAUTC(JL)**2)
        ELSE
          ZHLC_LRCLOCAL(JL)=0.
          ZHLC_HRCLOCAL(JL)=0.
        END IF
      END DO
    END IF
    ! Compare r_cM  to r_cR to know if cloud water content is high enough to split in two parts or not
    DO JL=1, KSIZE
      IF (PRCT(JL).GT.0. .AND. PCF(JL).GT.0. .AND. ZHLC_RCMAX(JL).GT.ZRCRAUTC(JL) .AND. LDMICRO(JL)) THEN
        ! Calculate final values for LCF and HCF:
        PHLC_LCF(JL)=PCF(JL) &
                      *(ZHLC_HRCLOCAL(JL)- &
                      (PRCT(JL) / PCF(JL))) &
                      / (ZHLC_HRCLOCAL(JL)-ZHLC_LRCLOCAL(JL))
        PHLC_HCF(JL)=MAX(0., PCF(JL)-PHLC_LCF(JL))
        !
        ! Calculate final values for LRC and HRC:
        PHLC_LRC(JL)=ZHLC_LRCLOCAL(JL)*PHLC_LCF(JL)
        PHLC_HRC(JL)=MAX(0., PRCT(JL)-PHLC_LRC(JL))
      ELSE IF (PRCT(JL).GT.0. .AND. PCF(JL).GT.0. .AND. ZHLC_RCMAX(JL).LE.ZRCRAUTC(JL) .AND. LDMICRO(JL)) THEN
        ! Put all available cloud water and his fraction in the low part
        PHLC_LCF(JL)=PCF(JL)
        PHLC_HCF(JL)=0.
        PHLC_LRC(JL)=PRCT(JL)
        PHLC_HRC(JL)=0.
      ELSE
        PHLC_LCF(JL)=0.
        PHLC_HCF(JL)=0.
        PHLC_LRC(JL)=0.
        PHLC_HRC(JL)=0.
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

DO JL=1, KSIZE
  IF (LDMICRO(JL)) THEN
    ZCRIAUTI(JL)=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(PT(JL)-CST%XTT)+ICEP%XBCRIAUTI)) ! Autoconversion ri threshold
  ELSE
    ZCRIAUTI(JL)=0.
  END IF
END DO

IF(HSUBG_AUCV_RI=='NONE') THEN
  !Cloud water is entirely in low or high part

  DO JL=1, KSIZE
    IF (PRIT(JL)>ZCRIAUTI(JL) .AND. LDMICRO(JL)) THEN
      PHLI_HCF(JL)=1.
      PHLI_LCF(JL)=0.
      PHLI_HRI(JL)=PRIT(JL)
      PHLI_LRI(JL)=0.
    ELSE IF (PRIT(JL)>ICED%XRTMIN(4) .AND. LDMICRO(JL)) THEN
      PHLI_HCF(JL)=0.
      PHLI_LCF(JL)=1.
      PHLI_HRI(JL)=0.
      PHLI_LRI(JL)=PRIT(JL)
    ELSE
      PHLI_HCF(JL)=0.
      PHLI_LCF(JL)=0.
      PHLI_HRI(JL)=0.
      PHLI_LRI(JL)=0.
    END IF
  END DO

ELSEIF(HSUBG_AUCV_RI=='CLFR') THEN
  !Cloud water is only in the cloudy part and entirely in low or high part

  DO JL=1, KSIZE
    IF (PCF(JL)>0. .AND. PRIT(JL)>ZCRIAUTI(JL)*PCF(JL) .AND. LDMICRO(JL)) THEN
      PHLI_HCF(JL)=PCF(JL)
      PHLI_LCF(JL)=0.
      PHLI_HRI(JL)=PRIT(JL)
      PHLI_LRI(JL)=0.
    ELSE IF (PCF(JL)>0. .AND. PRIT(JL)>ICED%XRTMIN(4) .AND. LDMICRO(JL)) THEN
      PHLI_HCF(JL)=0.
      PHLI_LCF(JL)=PCF(JL)
      PHLI_HRI(JL)=0.0
      PHLI_LRI(JL)=PRIT(JL)
    ELSE
      PHLI_HCF(JL)=0.
      PHLI_LCF(JL)=0.
      PHLI_HRI(JL)=0.
      PHLI_LRI(JL)=0.
    END IF
  END DO

ELSEIF(HSUBG_AUCV_RI=='ADJU') THEN

  DO JL=1, KSIZE
    IF (LDMICRO(JL)) THEN
      ZSUMRI(JL)=PHLI_LRI(JL)+PHLI_HRI(JL)
    ELSE
      ZSUMRI(JL)=0.
    END IF
    IF (ZSUMRI(JL) .GT. 1.E-20 .AND. LDMICRO(JL)) THEN
      PHLI_LRI(JL)=PHLI_LRI(JL)*PRIT(JL)/ZSUMRI(JL)
      PHLI_HRI(JL)=PHLI_HRI(JL)*PRIT(JL)/ZSUMRI(JL)
    ELSE
      PHLI_LRI(JL)=0.
      PHLI_HRI(JL)=0.
    END IF
  END DO

ELSE
  !wrong HSUBG_AUCV_RI case
  CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'ICE4_COMPUTE_PDF', 'wrong HSUBG_AUCV_RI case' )
ENDIF
!

DO JL=1, KSIZE
  IF (LDMICRO(JL)) THEN
    PRF(JL)=MAX(PHLC_HCF(JL),PHLI_HCF(JL))
  ELSE
    PRF(JL)=0.
  END IF
END DO

!
IF (LHOOK) CALL DR_HOOK('ICE4_COMPUTE_PDF', 1, ZHOOK_HANDLE)
END SUBROUTINE ICE4_COMPUTE_PDF

END MODULE MODE_ICE4_COMPUTE_PDF
