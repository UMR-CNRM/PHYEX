!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! ******* TO BE INCLUDED IN THE *CONTAINS* OF A SUBROUTINE, IN ORDER TO EASE AUTOMATIC INLINING ******
! => Don't use drHook !!!
!
ELEMENTAL SUBROUTINE ICE4_COMPUTE_PDF(CST, ICEP, ICED, HSUBG_AUCV_RC, HSUBG_AUCV_RI, HSUBG_PR_PDF, &
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
!!      S. Riette Nov 26: elemental version
!
!
!*      0. DECLARATIONS
!          ------------
!
!
USE MODD_CST,            ONLY: CST_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
!
IMPLICIT NONE
!

!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
CHARACTER(LEN=4),       INTENT(IN)  :: HSUBG_AUCV_RC     ! Kind of Subgrid autoconversion method for cloud water
CHARACTER(LEN=80),      INTENT(IN)  :: HSUBG_AUCV_RI     ! Kind of Subgrid autoconversion method for cloud ice
CHARACTER(LEN=80),      INTENT(IN)  :: HSUBG_PR_PDF   ! pdf for subgrid precipitation
LOGICAL, INTENT(IN)  :: LDMICRO    ! Computation mask
REAL,    INTENT(IN)  :: PRHODREF   ! Reference density
REAL,    INTENT(IN)  :: PRCT       ! Cloud water m.r. at t
REAL,    INTENT(IN)  :: PRIT       ! Ice Crystal m.r. at t
REAL,    INTENT(IN)  :: PCF        ! Cloud fraction
REAL,    INTENT(IN)  :: PT         ! Temperature
REAL,    INTENT(IN)  :: PSIGMA_RC  ! Standard deviation of rc at time t
!Note for INTENT STATUS: in 'ADJU' case the PHL?_??? variables must be able to "cross" the subroutine untouched
REAL,    INTENT(INOUT) :: PHLC_HCF   ! HLCLOUDS : fraction of High Cloud Fraction in grid
REAL,    INTENT(INOUT) :: PHLC_LCF   ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                     !    note that PCF = PHLC_HCF + PHLC_LCF
REAL,    INTENT(INOUT) :: PHLC_HRC   ! HLCLOUDS : LWC that is High LWC in grid
REAL,    INTENT(INOUT) :: PHLC_LRC   ! HLCLOUDS : LWC that is Low  LWC in grid
                                     !    note that PRC = PHLC_HRC + PHLC_LRC
REAL,    INTENT(INOUT) :: PHLI_HCF
REAL,    INTENT(INOUT) :: PHLI_LCF
REAL,    INTENT(INOUT) :: PHLI_HRI
REAL,    INTENT(INOUT) :: PHLI_LRI
REAL,    INTENT(OUT)   :: PRF        ! Rain fraction
!
!*       0.2  declaration of local variables
!
REAL :: ZRCRAUTC,      & !RC value to begin rain formation =XCRIAUTC/RHODREF
        ZCRIAUTI,      & !RI value to begin snow formation
        ZHLC_RCMAX,    & !HLCLOUDS : maximum value for RC in distribution
        ZHLC_LRCLOCAL, & !HLCLOUDS : LWC that is Low  LWC local in LCF
        ZHLC_HRCLOCAL, & !HLCLOUDS : LWC that is High LWC local in HCF
                         !    note that ZRC/CF = ZHLC_HRCLOCAL+ ZHLC_LRCLOCAL
                         !                     = PHLC_HRC/HCF+ PHLC_LRC/LCF
        ZSUMRC, ZSUMRI
REAL :: ZCOEFFRCM
!-------------------------------------------------------------------------------
!
!Cloud water split between high and low content part is done according to autoconversion option
IF (LDMICRO) THEN
  ZRCRAUTC=ICEP%XCRIAUTC/PRHODREF ! Autoconversion rc threshold
ELSE
  ZRCRAUTC=0.
END IF
IF(HSUBG_AUCV_RC=='NONE') THEN
  !Cloud water is entirely in low or high part
  IF(.NOT. LDMICRO) THEN
    PHLC_HCF=0.
    PHLC_LCF=0.
    PHLC_HRC=0.
    PHLC_LRC=0.
  ELSE IF(PRCT>ZRCRAUTC) THEN
    PHLC_HCF=1.
    PHLC_LCF=0.
    PHLC_HRC=PRCT
    PHLC_LRC=0.
  ELSE IF(PRCT>ICED%XRTMIN(2)) THEN
    PHLC_HCF=0.
    PHLC_LCF=1.
    PHLC_HRC=0.
    PHLC_LRC=PRCT
  ELSE
    PHLC_HCF=0.
    PHLC_LCF=0.
    PHLC_HRC=0.
    PHLC_LRC=0.
  END IF
ELSEIF(HSUBG_AUCV_RC=='CLFR') THEN
  !Cloud water is only in the cloudy part and entirely in low or high part
  IF(.NOT. LDMICRO) THEN
    PHLC_HCF=0.
    PHLC_LCF=0.
    PHLC_HRC=0.
    PHLC_LRC=0.
  ELSE IF(PCF>0. .AND. PRCT>ZRCRAUTC*PCF) THEN
    PHLC_HCF=PCF
    PHLC_LCF=0.
    PHLC_HRC=PRCT
    PHLC_LRC=0.
  ELSE IF(PCF>0. .AND. PRCT>ICED%XRTMIN(2)) THEN
    PHLC_HCF=0.
    PHLC_LCF=PCF
    PHLC_HRC=0.0
    PHLC_LRC=PRCT
  ELSE
    PHLC_HCF=0.
    PHLC_LCF=0.
    PHLC_HRC=0.
    PHLC_LRC=0.
  END IF
ELSEIF(HSUBG_AUCV_RC=='ADJU') THEN
  IF(LDMICRO) THEN
    ZSUMRC=PHLC_LRC+PHLC_HRC
  ELSE
    ZSUMRC=0.
  ENDIF
  IF(ZSUMRC .GT. 1.E-20 .AND. LDMICRO) THEN
    PHLC_LRC=PHLC_LRC*PRCT/ZSUMRC
    PHLC_HRC=PHLC_HRC*PRCT/ZSUMRC
  ELSE
    PHLC_LRC=0.
    PHLC_HRC=0.
  ENDIF
ELSEIF(HSUBG_AUCV_RC=='PDF') THEN
  !Cloud water is split between high and low part according to a PDF
  !    'HLCRECTPDF'    : rectangular PDF form
  !    'HLCTRIANGPDF'  : triangular PDF form
  !    'HLCQUADRAPDF'  : second order quadratic PDF form
  !    'HLCISOTRIPDF'  : isocele triangular PDF
  !    'SIGM'          : Redelsperger and Sommeria (1986)
  IF(HSUBG_PR_PDF=='SIGM') THEN
    ! Redelsperger and Sommeria (1986) but organised according to Turner (2011, 2012)
    IF(.NOT. LDMICRO) THEN
      PHLC_HCF=0.
      PHLC_LCF=0.
      PHLC_HRC=0.
      PHLC_LRC=0.
    ELSE IF(PRCT>ZRCRAUTC+PSIGMA_RC) THEN
      PHLC_HCF=1.
      PHLC_LCF=0.
      PHLC_HRC=PRCT
      PHLC_LRC=0.
    ELSE IF(PRCT> (ZRCRAUTC-PSIGMA_RC) .AND. PRCT<=(ZRCRAUTC+PSIGMA_RC)) THEN
      PHLC_HCF=(PRCT+PSIGMA_RC-ZRCRAUTC)/ &
               &(2.*PSIGMA_RC)
      PHLC_LCF=MAX(0., PCF-PHLC_HCF)
      PHLC_HRC=(PRCT+PSIGMA_RC-ZRCRAUTC)* &
               &(PRCT+PSIGMA_RC+ZRCRAUTC)/ &
               &(4.*PSIGMA_RC)
      PHLC_LRC=MAX(0., PRCT-PHLC_HRC)
    ELSE IF(PRCT>ICED%XRTMIN(2) .AND. PCF>0.) THEN
      PHLC_HCF=0.
      PHLC_LCF=PCF
      PHLC_HRC=0.
      PHLC_LRC=PRCT
    ELSE
      PHLC_HCF=0.
      PHLC_LCF=0.
      PHLC_HRC=0.
      PHLC_LRC=0.
    END IF
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
    IF(.NOT. LDMICRO) THEN
      ZHLC_RCMAX=0.
    ELSE IF(PRCT.GT.0. .AND. PCF.GT.0.) THEN
      ZHLC_RCMAX=ZCOEFFRCM*PRCT/PCF
    ELSE
      ZHLC_RCMAX=0.
    END IF

    ! Split available water and cloud fraction in two parts
    ! Calculate local mean values int he low and high parts for the 3 PDF forms:
    IF(HSUBG_PR_PDF=='HLCRECTPDF') THEN
      IF(.NOT. LDMICRO) THEN
        ZHLC_LRCLOCAL=0.
        ZHLC_HRCLOCAL=0
      ELSE IF(PRCT.GT.0. .AND. PCF.GT.0. .AND. ZHLC_RCMAX.GT.ZRCRAUTC) THEN
        ZHLC_LRCLOCAL=0.5*ZRCRAUTC
        ZHLC_HRCLOCAL=( ZHLC_RCMAX + ZRCRAUTC)/2.0
      ELSE
        ZHLC_LRCLOCAL=0.
        ZHLC_HRCLOCAL=0.
      END IF
    ELSE IF(HSUBG_PR_PDF=='HLCTRIANGPDF') THEN
      IF(.NOT. LDMICRO) THEN
        ZHLC_LRCLOCAL=0.
        ZHLC_HRCLOCAL=0
      ELSE IF(PRCT.GT.0. .AND. PCF.GT.0. .AND. ZHLC_RCMAX.GT.ZRCRAUTC) THEN
        ZHLC_LRCLOCAL=(ZRCRAUTC *(3.0 * ZHLC_RCMAX - 2.0 * ZRCRAUTC)) &
                      / (3.0 * (2.0 * ZHLC_RCMAX - ZRCRAUTC  ) )
        ZHLC_HRCLOCAL=(ZHLC_RCMAX + 2.0*ZRCRAUTC) / 3.0
      ELSE
        ZHLC_LRCLOCAL=0.
        ZHLC_HRCLOCAL=0.
      END IF
    ELSE IF(HSUBG_PR_PDF=='HLCQUADRAPDF') THEN
      IF(.NOT. LDMICRO) THEN
        ZHLC_LRCLOCAL=0.
        ZHLC_HRCLOCAL=0
      ELSE IF(PRCT.GT.0. .AND. PCF.GT.0. .AND. ZHLC_RCMAX.GT.ZRCRAUTC) THEN
        ZHLC_LRCLOCAL=(3.0 *ZRCRAUTC**3 &
                      - 8.0 *ZRCRAUTC**2 * ZHLC_RCMAX &
                      + 6.0*ZRCRAUTC *ZHLC_RCMAX**2 ) &
                      / &
                      (4.0* ZRCRAUTC**2 -12.0*ZRCRAUTC *ZHLC_RCMAX &
                      + 12.0 * ZHLC_RCMAX**2 )
        ZHLC_HRCLOCAL=(ZHLC_RCMAX + 3.0*ZRCRAUTC)/4.0
      ELSE
        ZHLC_LRCLOCAL=0.
        ZHLC_HRCLOCAL=0.
      END IF
    ELSE IF(HSUBG_PR_PDF=='HLCISOTRIPDF') THEN
      IF(.NOT. LDMICRO) THEN
        ZHLC_LRCLOCAL=0.
        ZHLC_HRCLOCAL=0
      ELSE IF(PRCT.LE.ZRCRAUTC*PCF .AND. &
             &PRCT.GT.0. .AND. PCF.GT.0. .AND. &
             &ZHLC_RCMAX.GT.ZRCRAUTC) THEN
        ZHLC_LRCLOCAL=( (ZHLC_RCMAX)**3 &
                      -(12.0 * (ZHLC_RCMAX)*(ZRCRAUTC)**2) &
                      +(8.0 * ZRCRAUTC**3) ) &
                      /( (6.0 * (ZHLC_RCMAX)**2) &
                      -(24.0 * (ZHLC_RCMAX) * ZRCRAUTC) &
                      +(12.0 * ZRCRAUTC**2) )
        ZHLC_HRCLOCAL=( ZHLC_RCMAX + 2.0 * ZRCRAUTC )/3.0
      ELSE IF(PRCT.GT.0. .AND. PCF.GT.0. .AND. ZHLC_RCMAX.GT.ZRCRAUTC) THEN
        ZHLC_LRCLOCAL=(2.0/3.0) * ZRCRAUTC
        ZHLC_HRCLOCAL=(3.0*ZHLC_RCMAX**3 - 8.0*ZRCRAUTC**3) &
                      / (6.0 * ZHLC_RCMAX**2 - 12.0*ZRCRAUTC**2)
      ELSE
        ZHLC_LRCLOCAL=0.
        ZHLC_HRCLOCAL=0.
      END IF
    END IF
    ! Compare r_cM  to r_cR to know if cloud water content is high enough to split in two parts or not
    IF(.NOT. LDMICRO) THEN
      PHLC_HCF=0.
      PHLC_LCF=0.
      PHLC_HRC=0.
      PHLC_LRC=0.
    ELSE IF(PRCT.GT.0. .AND. PCF.GT.0. .AND. ZHLC_RCMAX.GT.ZRCRAUTC) THEN
      ! Calculate final values for LCF and HCF:
      PHLC_LCF=PCF &
               *(ZHLC_HRCLOCAL- &
               (PRCT / PCF)) &
               / (ZHLC_HRCLOCAL-ZHLC_LRCLOCAL)
      PHLC_HCF=MAX(0., PCF-PHLC_LCF)
      !
      ! Calculate final values for LRC and HRC:
      PHLC_LRC=ZHLC_LRCLOCAL*PHLC_LCF
      PHLC_HRC=MAX(0., PRCT-PHLC_LRC)
    ELSE IF(PRCT.GT.0. .AND. PCF.GT.0. .AND. ZHLC_RCMAX.LE.ZRCRAUTC) THEN
      ! Put all available cloud water and his fraction in the low part
      PHLC_LCF=PCF
      PHLC_HCF=0.
      PHLC_LRC=PRCT
      PHLC_HRC=0.
    ELSE
      PHLC_LCF=0.
      PHLC_HCF=0.
      PHLC_LRC=0.
      PHLC_HRC=0.
    END IF
  ENDIF
ENDIF
!
!Ice water split between high and low content part is done according to autoconversion option
IF(LDMICRO) THEN
  ! Autoconversion ri threshold
  ZCRIAUTI=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(PT-CST%XTT)+ICEP%XBCRIAUTI))
ELSE
  ZCRIAUTI=0.
ENDIF
IF(HSUBG_AUCV_RI=='NONE') THEN
  !Cloud water is entirely in low or high part
  IF(.NOT. LDMICRO) THEN
    PHLI_HCF=0.
    PHLI_LCF=0.
    PHLI_HRI=0.
    PHLI_LRI=0.
  ELSE IF(PRIT>ZCRIAUTI) THEN
    PHLI_HCF=1.
    PHLI_LCF=0.
    PHLI_HRI=PRIT
    PHLI_LRI=0.
  ELSE IF(PRIT>ICED%XRTMIN(4)) THEN
    PHLI_HCF=0.
    PHLI_LCF=1.
    PHLI_HRI=0.
    PHLI_LRI=PRIT
  ELSE
    PHLI_HCF=0.
    PHLI_LCF=0.
    PHLI_HRI=0.
    PHLI_LRI=0.
  END IF
ELSEIF(HSUBG_AUCV_RI=='CLFR') THEN
  !Cloud water is only in the cloudy part and entirely in low or high part
  IF(.NOT. LDMICRO) THEN
    PHLI_HCF=0.
    PHLI_LCF=0.
    PHLI_HRI=0.
    PHLI_LRI=0.
  ELSE IF(PCF>0. .AND. PRIT>ZCRIAUTI*PCF) THEN
    PHLI_HCF=PCF
    PHLI_LCF=0.
    PHLI_HRI=PRIT
    PHLI_LRI=0.
  ELSE IF(PCF>0. .AND. PRIT>ICED%XRTMIN(4)) THEN
    PHLI_HCF=0.
    PHLI_LCF=PCF
    PHLI_HRI=0.0
    PHLI_LRI=PRIT
  ELSE
    PHLI_HCF=0.
    PHLI_LCF=0.
    PHLI_HRI=0.
    PHLI_LRI=0.
  END IF
ELSEIF(HSUBG_AUCV_RI=='ADJU') THEN
  IF(LDMICRO) THEN
    ZSUMRI=PHLI_LRI+PHLI_HRI
  ELSE
    ZSUMRI=0.
  ENDIF
  IF(ZSUMRI .GT. 1.E-20 .AND. LDMICRO) THEN
    PHLI_LRI=PHLI_LRI*PRIT/ZSUMRI
    PHLI_HRI=PHLI_HRI*PRIT/ZSUMRI
  ELSE
    PHLI_LRI=0.
    PHLI_HRI=0.
  ENDIF
ENDIF
!
IF(LDMICRO) THEN
  PRF=MAX(PHLC_HCF,PHLI_HCF)
ELSE
  PRF=0.
ENDIF
!
END SUBROUTINE ICE4_COMPUTE_PDF
