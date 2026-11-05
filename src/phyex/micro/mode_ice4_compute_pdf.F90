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

!$ACDC singlecolumn

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
USE MODE_MSG, ONLY: PRINT_MSG
USE MODD_IO, ONLY:NVERB_FATAL
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
!$acc kernels
!$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
WHERE (LDMICRO(D%NIJB:D%NIJE))
  ZRCRAUTC(D%NIJB:D%NIJE)=ICEP%XCRIAUTC/PRHODREF(D%NIJB:D%NIJE) ! Autoconversion rc threshold
ELSEWHERE
  ZRCRAUTC(D%NIJB:D%NIJE)=0.
END WHERE
!$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
IF(HSUBG_AUCV_RC=='NONE') THEN
  !Cloud water is entirely in low or high part
!$acc kernels
 !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
  WHERE(.NOT. LDMICRO(D%NIJB:D%NIJE))
    ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
    PHLC_HCF(D%NIJB:D%NIJE)=0.
    PHLC_LCF(D%NIJB:D%NIJE)=0.
    PHLC_HRC(D%NIJB:D%NIJE)=0.
    PHLC_LRC(D%NIJB:D%NIJE)=0.
  ELSEWHERE(PRCT(D%NIJB:D%NIJE)>ZRCRAUTC(D%NIJB:D%NIJE) .AND. LDMICRO(D%NIJB:D%NIJE))
    PHLC_HCF(D%NIJB:D%NIJE)=1.
    PHLC_LCF(D%NIJB:D%NIJE)=0.
    PHLC_HRC(D%NIJB:D%NIJE)=PRCT(D%NIJB:D%NIJE)
    PHLC_LRC(D%NIJB:D%NIJE)=0.
  ELSEWHERE(PRCT(D%NIJB:D%NIJE)>ICED%XRTMIN(2) .AND. LDMICRO(D%NIJB:D%NIJE))
    PHLC_HCF(D%NIJB:D%NIJE)=0.
    PHLC_LCF(D%NIJB:D%NIJE)=1.
    PHLC_HRC(D%NIJB:D%NIJE)=0.
    PHLC_LRC(D%NIJB:D%NIJE)=PRCT(D%NIJB:D%NIJE)
  ELSEWHERE
    PHLC_HCF(D%NIJB:D%NIJE)=0.
    PHLC_LCF(D%NIJB:D%NIJE)=0.
    PHLC_HRC(D%NIJB:D%NIJE)=0.
    PHLC_LRC(D%NIJB:D%NIJE)=0.
  END WHERE
  !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
ELSEIF(HSUBG_AUCV_RC=='CLFR') THEN
  !Cloud water is only in the cloudy part and entirely in low or high part
!$acc kernels
 !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
  WHERE(.NOT. LDMICRO(D%NIJB:D%NIJE))
    ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
    PHLC_HCF(D%NIJB:D%NIJE)=0.
    PHLC_LCF(D%NIJB:D%NIJE)=0.
    PHLC_HRC(D%NIJB:D%NIJE)=0.
    PHLC_LRC(D%NIJB:D%NIJE)=0.
  ELSEWHERE(PCF(D%NIJB:D%NIJE)>0. .AND. PRCT(D%NIJB:D%NIJE)>ZRCRAUTC(D%NIJB:D%NIJE)*PCF(D%NIJB:D%NIJE) &
           &.AND. LDMICRO(D%NIJB:D%NIJE))
    PHLC_HCF(D%NIJB:D%NIJE)=PCF(D%NIJB:D%NIJE)
    PHLC_LCF(D%NIJB:D%NIJE)=0.
    PHLC_HRC(D%NIJB:D%NIJE)=PRCT(D%NIJB:D%NIJE)
    PHLC_LRC(D%NIJB:D%NIJE)=0.
  ELSEWHERE(PCF(D%NIJB:D%NIJE)>0. .AND. PRCT(D%NIJB:D%NIJE)>ICED%XRTMIN(2) .AND. LDMICRO(D%NIJB:D%NIJE))
    PHLC_HCF(D%NIJB:D%NIJE)=0.
    PHLC_LCF(D%NIJB:D%NIJE)=PCF(D%NIJB:D%NIJE)
    PHLC_HRC(D%NIJB:D%NIJE)=0.0
    PHLC_LRC(D%NIJB:D%NIJE)=PRCT(D%NIJB:D%NIJE)
  ELSEWHERE
    PHLC_HCF(D%NIJB:D%NIJE)=0.
    PHLC_LCF(D%NIJB:D%NIJE)=0.
    PHLC_HRC(D%NIJB:D%NIJE)=0.
    PHLC_LRC(D%NIJB:D%NIJE)=0.
  END WHERE
  !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
ELSEIF(HSUBG_AUCV_RC=='ADJU') THEN
!$acc kernels
  !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
  WHERE(LDMICRO(D%NIJB:D%NIJE))
    ZSUMRC(D%NIJB:D%NIJE)=PHLC_LRC(D%NIJB:D%NIJE)+PHLC_HRC(D%NIJB:D%NIJE)
  ELSEWHERE
    ZSUMRC(D%NIJB:D%NIJE)=0.
  ENDWHERE
  WHERE(ZSUMRC(D%NIJB:D%NIJE) .GT. 1.E-20 .AND. LDMICRO(D%NIJB:D%NIJE))
    PHLC_LRC(D%NIJB:D%NIJE)=PHLC_LRC(D%NIJB:D%NIJE)*PRCT(D%NIJB:D%NIJE)/ZSUMRC(D%NIJB:D%NIJE)
    PHLC_HRC(D%NIJB:D%NIJE)=PHLC_HRC(D%NIJB:D%NIJE)*PRCT(D%NIJB:D%NIJE)/ZSUMRC(D%NIJB:D%NIJE)
  ELSEWHERE
    PHLC_LRC(D%NIJB:D%NIJE)=0.
    PHLC_HRC(D%NIJB:D%NIJE)=0.
  ENDWHERE
  !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
ELSEIF(HSUBG_AUCV_RC=='PDF') THEN
  !Cloud water is split between high and low part according to a PDF
  !    'HLCRECTPDF'    : rectangular PDF form
  !    'HLCTRIANGPDF'  : triangular PDF form
  !    'HLCQUADRAPDF'  : second order quadratic PDF form
  !    'HLCISOTRIPDF'  : isocele triangular PDF
  !    'SIGM'          : Redelsperger and Sommeria (1986)
  IF(HSUBG_PR_PDF=='SIGM') THEN
    ! Redelsperger and Sommeria (1986) but organised according to Turner (2011, 2012)
!$acc kernels
    !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
    WHERE(.NOT. LDMICRO(D%NIJB:D%NIJE))
      ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
      PHLC_HCF(D%NIJB:D%NIJE)=0.
      PHLC_LCF(D%NIJB:D%NIJE)=0.
      PHLC_HRC(D%NIJB:D%NIJE)=0.
      PHLC_LRC(D%NIJB:D%NIJE)=0.
    ELSEWHERE(PRCT(D%NIJB:D%NIJE)>ZRCRAUTC(D%NIJB:D%NIJE)+PSIGMA_RC(D%NIJB:D%NIJE) .AND. LDMICRO(D%NIJB:D%NIJE))
      PHLC_HCF(D%NIJB:D%NIJE)=1.
      PHLC_LCF(D%NIJB:D%NIJE)=0.
      PHLC_HRC(D%NIJB:D%NIJE)=PRCT(D%NIJB:D%NIJE)
      PHLC_LRC(D%NIJB:D%NIJE)=0.
    ELSEWHERE(PRCT(D%NIJB:D%NIJE)> (ZRCRAUTC(D%NIJB:D%NIJE)-PSIGMA_RC(D%NIJB:D%NIJE)) .AND. &
             &PRCT(D%NIJB:D%NIJE)<=(ZRCRAUTC(D%NIJB:D%NIJE)+PSIGMA_RC(D%NIJB:D%NIJE)) .AND. LDMICRO(D%NIJB:D%NIJE))
      PHLC_HCF(D%NIJB:D%NIJE)=(PRCT(D%NIJB:D%NIJE)+PSIGMA_RC(D%NIJB:D%NIJE)-ZRCRAUTC(D%NIJB:D%NIJE))/ &
                  &(2.*PSIGMA_RC(D%NIJB:D%NIJE))
      PHLC_LCF(D%NIJB:D%NIJE)=MAX(0., PCF(D%NIJB:D%NIJE)-PHLC_HCF(D%NIJB:D%NIJE))
      PHLC_HRC(D%NIJB:D%NIJE)=(PRCT(D%NIJB:D%NIJE)+PSIGMA_RC(D%NIJB:D%NIJE)-ZRCRAUTC(D%NIJB:D%NIJE))* &
                  &(PRCT(D%NIJB:D%NIJE)+PSIGMA_RC(D%NIJB:D%NIJE)+ZRCRAUTC(D%NIJB:D%NIJE))/ &
                  &(4.*PSIGMA_RC(D%NIJB:D%NIJE))
      PHLC_LRC(D%NIJB:D%NIJE)=MAX(0., PRCT(D%NIJB:D%NIJE)-PHLC_HRC(D%NIJB:D%NIJE))
    ELSEWHERE(PRCT(D%NIJB:D%NIJE)>ICED%XRTMIN(2) .AND. PCF(D%NIJB:D%NIJE)>0. .AND. LDMICRO(D%NIJB:D%NIJE))
      PHLC_HCF(D%NIJB:D%NIJE)=0.
      PHLC_LCF(D%NIJB:D%NIJE)=PCF(D%NIJB:D%NIJE)
      PHLC_HRC(D%NIJB:D%NIJE)=0.
      PHLC_LRC(D%NIJB:D%NIJE)=PRCT(D%NIJB:D%NIJE)
    ELSEWHERE
      PHLC_HCF(D%NIJB:D%NIJE)=0.
      PHLC_LCF(D%NIJB:D%NIJE)=0.
      PHLC_HRC(D%NIJB:D%NIJE)=0.
      PHLC_LRC(D%NIJB:D%NIJE)=0.
    END WHERE
    !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
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
!$acc kernels
    !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
    WHERE(.NOT. LDMICRO(D%NIJB:D%NIJE))
      ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
      ZHLC_RCMAX(D%NIJB:D%NIJE)=0.
    ELSEWHERE(PRCT(D%NIJB:D%NIJE).GT.0. .AND. PCF(D%NIJB:D%NIJE).GT.0. .AND. LDMICRO(D%NIJB:D%NIJE))
      ZHLC_RCMAX(D%NIJB:D%NIJE)=ZCOEFFRCM*PRCT(D%NIJB:D%NIJE)/PCF(D%NIJB:D%NIJE)
    ELSEWHERE
      ZHLC_RCMAX(D%NIJB:D%NIJE)=0.
    END WHERE
    !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)

    ! Split available water and cloud fraction in two parts
    ! Calculate local mean values int he low and high parts for the 3 PDF forms:
    IF(HSUBG_PR_PDF=='HLCRECTPDF') THEN
      !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
      WHERE(.NOT. LDMICRO(D%NIJB:D%NIJE))
        ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=0.
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=0
      ELSEWHERE(PRCT(D%NIJB:D%NIJE).GT.0. .AND. PCF(D%NIJB:D%NIJE).GT.0. .AND. &
               &ZHLC_RCMAX(D%NIJB:D%NIJE).GT.ZRCRAUTC(D%NIJB:D%NIJE) .AND. LDMICRO(D%NIJB:D%NIJE))
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=0.5*ZRCRAUTC(D%NIJB:D%NIJE)
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=( ZHLC_RCMAX(D%NIJB:D%NIJE) + ZRCRAUTC(D%NIJB:D%NIJE))/2.0
      ELSEWHERE
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=0.
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=0.
      END WHERE
      !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
    ELSE IF(HSUBG_PR_PDF=='HLCTRIANGPDF') THEN
      !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
      WHERE(.NOT. LDMICRO(D%NIJB:D%NIJE))
        ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=0.
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=0
      ELSEWHERE(PRCT(D%NIJB:D%NIJE).GT.0. .AND. PCF(D%NIJB:D%NIJE).GT.0. .AND. &
               &ZHLC_RCMAX(D%NIJB:D%NIJE).GT.ZRCRAUTC(D%NIJB:D%NIJE) .AND. LDMICRO(D%NIJB:D%NIJE))
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=(ZRCRAUTC(D%NIJB:D%NIJE) *(3.0 * ZHLC_RCMAX(D%NIJB:D%NIJE) - 2.0 * ZRCRAUTC(D%NIJB:D%NIJE))) &
                        / (3.0 * (2.0 * ZHLC_RCMAX(D%NIJB:D%NIJE) - ZRCRAUTC(D%NIJB:D%NIJE)  ) )
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=(ZHLC_RCMAX(D%NIJB:D%NIJE) + 2.0*ZRCRAUTC(D%NIJB:D%NIJE)) / 3.0
      ELSEWHERE
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=0.
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=0.
      END WHERE
      !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
    ELSE IF(HSUBG_PR_PDF=='HLCQUADRAPDF') THEN
      !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
      WHERE(.NOT. LDMICRO(D%NIJB:D%NIJE))
        ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=0.
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=0
      ELSEWHERE(PRCT(D%NIJB:D%NIJE).GT.0. .AND. PCF(D%NIJB:D%NIJE).GT.0. .AND. &
               &ZHLC_RCMAX(D%NIJB:D%NIJE).GT.ZRCRAUTC(D%NIJB:D%NIJE) .AND. LDMICRO(D%NIJB:D%NIJE))
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=(3.0 *ZRCRAUTC(D%NIJB:D%NIJE)**3 &
                        - 8.0 *ZRCRAUTC(D%NIJB:D%NIJE)**2 * ZHLC_RCMAX(D%NIJB:D%NIJE) &
                        + 6.0*ZRCRAUTC(D%NIJB:D%NIJE) *ZHLC_RCMAX(D%NIJB:D%NIJE)**2 ) &
                        / &
                        (4.0* ZRCRAUTC(D%NIJB:D%NIJE)**2 -12.0*ZRCRAUTC(D%NIJB:D%NIJE) *ZHLC_RCMAX(D%NIJB:D%NIJE) &
                        + 12.0 * ZHLC_RCMAX(D%NIJB:D%NIJE)**2 )
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=(ZHLC_RCMAX(D%NIJB:D%NIJE) + 3.0*ZRCRAUTC(D%NIJB:D%NIJE))/4.0
      ELSEWHERE
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=0.
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=0.
      END WHERE
      !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
    ELSE IF(HSUBG_PR_PDF=='HLCISOTRIPDF') THEN
      !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
      WHERE(.NOT. LDMICRO(D%NIJB:D%NIJE))
        ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=0.
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=0
      ELSEWHERE (PRCT(D%NIJB:D%NIJE).LE.ZRCRAUTC(D%NIJB:D%NIJE)*PCF(D%NIJB:D%NIJE) .AND. &
            &PRCT(D%NIJB:D%NIJE).GT.0. .AND. PCF(D%NIJB:D%NIJE).GT.0. .AND. &
            &ZHLC_RCMAX(D%NIJB:D%NIJE).GT.ZRCRAUTC(D%NIJB:D%NIJE) .AND. LDMICRO(D%NIJB:D%NIJE))
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=( (ZHLC_RCMAX(D%NIJB:D%NIJE))**3 &
                        -(12.0 * (ZHLC_RCMAX(D%NIJB:D%NIJE))*(ZRCRAUTC(D%NIJB:D%NIJE))**2) &
                        +(8.0 * ZRCRAUTC(D%NIJB:D%NIJE)**3) ) &
                        /( (6.0 * (ZHLC_RCMAX(D%NIJB:D%NIJE))**2) &
                        -(24.0 * (ZHLC_RCMAX(D%NIJB:D%NIJE)) * ZRCRAUTC(D%NIJB:D%NIJE)) &
                        +(12.0 * ZRCRAUTC(D%NIJB:D%NIJE)**2) )
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=( ZHLC_RCMAX(D%NIJB:D%NIJE) + 2.0 * ZRCRAUTC(D%NIJB:D%NIJE) )/3.0
      ELSEWHERE(PRCT(D%NIJB:D%NIJE).GT.0. .AND. PCF(D%NIJB:D%NIJE).GT.0. .AND. &
               &ZHLC_RCMAX(D%NIJB:D%NIJE).GT.ZRCRAUTC(D%NIJB:D%NIJE) .AND. LDMICRO(D%NIJB:D%NIJE))
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=(2.0/3.0) * ZRCRAUTC(D%NIJB:D%NIJE)
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=(3.0*ZHLC_RCMAX(D%NIJB:D%NIJE)**3 - 8.0*ZRCRAUTC(D%NIJB:D%NIJE)**3) &
                        / (6.0 * ZHLC_RCMAX(D%NIJB:D%NIJE)**2 - 12.0*ZRCRAUTC(D%NIJB:D%NIJE)**2)
      ELSEWHERE
        ZHLC_LRCLOCAL(D%NIJB:D%NIJE)=0.
        ZHLC_HRCLOCAL(D%NIJB:D%NIJE)=0.
      END WHERE
      !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
    END IF
    ! Compare r_cM  to r_cR to know if cloud water content is high enough to split in two parts or not
    !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
    WHERE(.NOT. LDMICRO(D%NIJB:D%NIJE))
      ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
      PHLC_HCF(D%NIJB:D%NIJE)=0.
      PHLC_LCF(D%NIJB:D%NIJE)=0.
      PHLC_HRC(D%NIJB:D%NIJE)=0.
      PHLC_LRC(D%NIJB:D%NIJE)=0.
    ELSEWHERE (PRCT(D%NIJB:D%NIJE).GT.0. .AND. PCF(D%NIJB:D%NIJE).GT.0. .AND. &
              &ZHLC_RCMAX(D%NIJB:D%NIJE).GT.ZRCRAUTC(D%NIJB:D%NIJE) .AND. LDMICRO(D%NIJB:D%NIJE))
      ! Calculate final values for LCF and HCF:
      PHLC_LCF(D%NIJB:D%NIJE)=PCF(D%NIJB:D%NIJE) &
                    *(ZHLC_HRCLOCAL(D%NIJB:D%NIJE)- &
                    (PRCT(D%NIJB:D%NIJE) / PCF(D%NIJB:D%NIJE))) &
                    / (ZHLC_HRCLOCAL(D%NIJB:D%NIJE)-ZHLC_LRCLOCAL(D%NIJB:D%NIJE))
      PHLC_HCF(D%NIJB:D%NIJE)=MAX(0., PCF(D%NIJB:D%NIJE)-PHLC_LCF(D%NIJB:D%NIJE))
      !
      ! Calculate final values for LRC and HRC:
      PHLC_LRC(D%NIJB:D%NIJE)=ZHLC_LRCLOCAL(D%NIJB:D%NIJE)*PHLC_LCF(D%NIJB:D%NIJE)
      PHLC_HRC(D%NIJB:D%NIJE)=MAX(0., PRCT(D%NIJB:D%NIJE)-PHLC_LRC(D%NIJB:D%NIJE))
    ELSEWHERE (PRCT(D%NIJB:D%NIJE).GT.0. .AND. PCF(D%NIJB:D%NIJE).GT.0. .AND. &
              &ZHLC_RCMAX(D%NIJB:D%NIJE).LE.ZRCRAUTC(D%NIJB:D%NIJE) .AND. LDMICRO(D%NIJB:D%NIJE))
      ! Put all available cloud water and his fraction in the low part
      PHLC_LCF(D%NIJB:D%NIJE)=PCF(D%NIJB:D%NIJE)
      PHLC_HCF(D%NIJB:D%NIJE)=0.
      PHLC_LRC(D%NIJB:D%NIJE)=PRCT(D%NIJB:D%NIJE)
      PHLC_HRC(D%NIJB:D%NIJE)=0.
    ELSEWHERE
      PHLC_LCF(D%NIJB:D%NIJE)=0.
      PHLC_HCF(D%NIJB:D%NIJE)=0.
      PHLC_LRC(D%NIJB:D%NIJE)=0.
      PHLC_HRC(D%NIJB:D%NIJE)=0.
    END WHERE
    !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
  ELSE
    CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_COMPUTE_PDF','wrong HSUBG_PR_PDF case')
  ENDIF
ELSE
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_COMPUTE_PDF','wrong HSUBG_AUCV case')
ENDIF
!
!Ice water split between high and low content part is done according to autoconversion option
!$acc kernels
!$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
WHERE(LDMICRO(D%NIJB:D%NIJE))
  ! Autoconversion ri threshold
  ZCRIAUTI(D%NIJB:D%NIJE)=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(PT(D%NIJB:D%NIJE)-CST%XTT)+ICEP%XBCRIAUTI))
ELSEWHERE
  ZCRIAUTI(D%NIJB:D%NIJE)=0.
ENDWHERE
!$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
IF(HSUBG_AUCV_RI=='NONE') THEN
  !Cloud water is entirely in low or high part
!$acc kernels
  !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
  WHERE(.NOT. LDMICRO(D%NIJB:D%NIJE))
    ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
    PHLI_HCF(D%NIJB:D%NIJE)=0.
    PHLI_LCF(D%NIJB:D%NIJE)=0.
    PHLI_HRI(D%NIJB:D%NIJE)=0.
    PHLI_LRI(D%NIJB:D%NIJE)=0.
  ELSEWHERE(PRIT(D%NIJB:D%NIJE)>ZCRIAUTI(D%NIJB:D%NIJE) .AND. LDMICRO(D%NIJB:D%NIJE))
    PHLI_HCF(D%NIJB:D%NIJE)=1.
    PHLI_LCF(D%NIJB:D%NIJE)=0.
    PHLI_HRI(D%NIJB:D%NIJE)=PRIT(D%NIJB:D%NIJE)
    PHLI_LRI(D%NIJB:D%NIJE)=0.
  ELSEWHERE(PRIT(D%NIJB:D%NIJE)>ICED%XRTMIN(4) .AND. LDMICRO(D%NIJB:D%NIJE))
    PHLI_HCF(D%NIJB:D%NIJE)=0.
    PHLI_LCF(D%NIJB:D%NIJE)=1.
    PHLI_HRI(D%NIJB:D%NIJE)=0.
    PHLI_LRI(D%NIJB:D%NIJE)=PRIT(D%NIJB:D%NIJE)
  ELSEWHERE
    PHLI_HCF(D%NIJB:D%NIJE)=0.
    PHLI_LCF(D%NIJB:D%NIJE)=0.
    PHLI_HRI(D%NIJB:D%NIJE)=0.
    PHLI_LRI(D%NIJB:D%NIJE)=0.
  END WHERE
  !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
ELSEIF(HSUBG_AUCV_RI=='CLFR') THEN
  !Cloud water is only in the cloudy part and entirely in low or high part
!$acc kernels
  !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
  WHERE(.NOT. LDMICRO(D%NIJB:D%NIJE))
    ! Needed to prevent evaluation, in AROME, of the next elseif (after mnh_expand transformation) condition
    PHLI_HCF(D%NIJB:D%NIJE)=0.
    PHLI_LCF(D%NIJB:D%NIJE)=0.
    PHLI_HRI(D%NIJB:D%NIJE)=0.
    PHLI_LRI(D%NIJB:D%NIJE)=0.
  ELSEWHERE(PCF(D%NIJB:D%NIJE)>0. .AND. PRIT(D%NIJB:D%NIJE)>ZCRIAUTI(D%NIJB:D%NIJE)*PCF(D%NIJB:D%NIJE) &
           &.AND. LDMICRO(D%NIJB:D%NIJE))
    PHLI_HCF(D%NIJB:D%NIJE)=PCF(D%NIJB:D%NIJE)
    PHLI_LCF(D%NIJB:D%NIJE)=0.
    PHLI_HRI(D%NIJB:D%NIJE)=PRIT(D%NIJB:D%NIJE)
    PHLI_LRI(D%NIJB:D%NIJE)=0.
  ELSEWHERE(PCF(D%NIJB:D%NIJE)>0. .AND. PRIT(D%NIJB:D%NIJE)>ICED%XRTMIN(4) .AND. LDMICRO(D%NIJB:D%NIJE))
    PHLI_HCF(D%NIJB:D%NIJE)=0.
    PHLI_LCF(D%NIJB:D%NIJE)=PCF(D%NIJB:D%NIJE)
    PHLI_HRI(D%NIJB:D%NIJE)=0.0
    PHLI_LRI(D%NIJB:D%NIJE)=PRIT(D%NIJB:D%NIJE)
  ELSEWHERE
    PHLI_HCF(D%NIJB:D%NIJE)=0.
    PHLI_LCF(D%NIJB:D%NIJE)=0.
    PHLI_HRI(D%NIJB:D%NIJE)=0.
    PHLI_LRI(D%NIJB:D%NIJE)=0.
  END WHERE
  !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
ELSEIF(HSUBG_AUCV_RI=='ADJU') THEN
!$acc kernels
  !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
  WHERE(LDMICRO(D%NIJB:D%NIJE))
    ZSUMRI(D%NIJB:D%NIJE)=PHLI_LRI(D%NIJB:D%NIJE)+PHLI_HRI(D%NIJB:D%NIJE)
  ELSEWHERE
    ZSUMRI(D%NIJB:D%NIJE)=0.
  ENDWHERE
  WHERE(ZSUMRI(D%NIJB:D%NIJE) .GT. 1.E-20 .AND. LDMICRO(D%NIJB:D%NIJE))
    PHLI_LRI(D%NIJB:D%NIJE)=PHLI_LRI(D%NIJB:D%NIJE)*PRIT(D%NIJB:D%NIJE)/ZSUMRI(D%NIJB:D%NIJE)
    PHLI_HRI(D%NIJB:D%NIJE)=PHLI_HRI(D%NIJB:D%NIJE)*PRIT(D%NIJB:D%NIJE)/ZSUMRI(D%NIJB:D%NIJE)
  ELSEWHERE
    PHLI_LRI(D%NIJB:D%NIJE)=0.
    PHLI_HRI(D%NIJB:D%NIJE)=0.
  ENDWHERE
  !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
ELSE
  !wrong HSUBG_AUCV_RI case
  CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'ICE4_COMPUTE_PDF', 'wrong HSUBG_AUCV_RI case' )
ENDIF
!
!$acc kernels
!$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
WHERE(LDMICRO(D%NIJB:D%NIJE))
  PRF(D%NIJB:D%NIJE)=MAX(PHLC_HCF(D%NIJB:D%NIJE),PHLI_HCF(D%NIJB:D%NIJE))
ELSEWHERE
  PRF(D%NIJB:D%NIJE)=0.
ENDWHERE
!$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
!
IF (LHOOK) CALL DR_HOOK('ICE4_COMPUTE_PDF', 1, ZHOOK_HANDLE)
END SUBROUTINE ICE4_COMPUTE_PDF

END MODULE MODE_ICE4_COMPUTE_PDF
