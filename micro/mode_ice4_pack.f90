!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_PACK
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_PACK(D, CST, PARAMI, ICEP, ICED, BUCONF,                   &
                     KPROMA, KSIZE, KSIZE2,                                &
                     PTSTEP, KRR, OSAVE_MICRO, ODMICRO, OELEC, PEXN,       &
                     PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,      &
                     PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,               &
                     PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,             &
                     PEVAP3D,                                              &
                     PRAINFR, PSIGS,                                       &
                     PRVHENI, PLVFACT, PLSFACT,                            &
                     PWR,                                                  &
                     TBUDGETS, KBUDGETS,                                   &
                     PMICRO_TEND, PLATHAM_IAGGS, PRHS                      )
!     ######################################################################
!
!!****  * -  compute the explicit microphysical sources
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to pack arrays to compute
!!      the microphysics tendencies
!!
!!
!!    METHOD
!!    ------
!!      Pack arrays by chuncks
!!
!!
!!    MODIFICATIONS
!!    -------------
!!     R. El Khatib 28-Apr-2023 Fix (and re-enable) the cache-blocking mechanism on top of phyex
!  -----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_BUDGET,         ONLY: TBUDGETDATA, TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_FIELDS_ADDRESS, ONLY : & ! common fields adress
      & ITH,     & ! Potential temperature
      & IRV,     & ! Water vapor
      & IRC,     & ! Cloud water
      & IRR,     & ! Rain water
      & IRI,     & ! Pristine ice
      & IRS,     & ! Snow/aggregate
      & IRG,     & ! Graupel
      & IRH,     & ! Hail
      & IBUNUM,       & ! Number of tendency terms
      & IBUNUM_EXTRA    ! Number of extra tendency terms

USE MODE_MSG,            ONLY: PRINT_MSG, NVERB_FATAL

USE MODE_ICE4_STEPPING, ONLY: ICE4_STEPPING
USE MODE_ICE4_BUDGETS, ONLY: ICE4_BUDGETS
!
IMPLICIT NONE

!NOTES ON SIZES
!If we pack:
! - KSIZE is the number of relevant point (with mixing ratio different from 0)
! - KPROMA is the size of bloc of points
! - ZSIZE2 has the same value as KPROMA
!If we do not pack:
! - KSIZE is the total number of points
! - KPROMA is null for memory saving
! - KSIZE2 has the same value as KSIZE
!
!When we do not pack, we can transmit directly the 3D arrays to the ice4_stepping subroutine, we do not need
!to copy the values. It is why KPROMA is null because we do not need these arrays.
!But some arrays must me manipulated before being transmitted and we need temporary arrays for this.
!KSIZE2 is used for those arrays that must be dimensioned KPROMA if we pack or with the total size if not.


!
!*       0.1   Declarations of dummy arguments :
!
!
!
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),        INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
INTEGER,                  INTENT(IN)    :: KPROMA ! cache-blocking factor for microphysic loop
INTEGER,                  INTENT(IN)    :: KSIZE
INTEGER,                  INTENT(IN)    :: KSIZE2
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL,                  INTENT(IN)    :: OSAVE_MICRO  ! If true, microphysical tendencies are saved
LOGICAL,                  INTENT(IN)    :: OELEC        ! if true, cloud electricity is activated
LOGICAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: ODMICRO ! mask to limit computation
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLC_HRC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLC_HCF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLI_HRI
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLI_HCF
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PEVAP3D! Rain evap profile
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRAINFR !Precipitation fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRVHENI ! heterogeneous nucleation
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PLVFACT
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PLSFACT
REAL, DIMENSION(D%NIJT,D%NKT,0:7), INTENT(INOUT) :: PWR
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
REAL, DIMENSION(MERGE(D%NIJT,0,OSAVE_MICRO),MERGE(D%NKT,0,OSAVE_MICRO),MERGE(IBUNUM-IBUNUM_EXTRA,0,OSAVE_MICRO)), &
                                          INTENT(INOUT) :: PMICRO_TEND  ! Microphysical tendencies
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), &
                                          INTENT(IN)    :: PLATHAM_IAGGS  ! E Function to simulate
                                                                          ! enhancement of IAGGS
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL,  INTENT(INOUT) :: PRHS    ! Hail m.r. source
!
!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
INTEGER :: JIJ, JK
INTEGER :: IKT, IKTB, IKTE, IIJT, IIJB, IIJE
INTEGER :: ISTIJ, ISTK
!
LOGICAL :: GEXT_TEND
!
!Output packed total mixing ratio change (for budgets only)
REAL, DIMENSION(KSIZE, IBUNUM-IBUNUM_EXTRA) :: ZBU_PACK
!
!For packing
INTEGER :: IMICRO ! Case r_x>0 locations
INTEGER :: JL, JV
REAL, DIMENSION(KPROMA) :: &
                        & ZCIT,     & ! Pristine ice conc. at t
                        & ZRHODREF, & ! RHO Dry REFerence
                        & ZPRES,    & ! Pressure
                        & ZEXN,     & ! EXNer Pressure
                        & ZCF,      & ! Cloud fraction
                        & ZHLC_HCF, & ! HLCLOUDS : fraction of High Cloud Fraction in grid
                        & ZHLC_HRC, & ! HLCLOUDS : LWC that is High LWC in grid
                        & ZHLI_HCF, &
                        & ZHLI_HRI, &
                        & ZRREVAV
REAL, DIMENSION(KSIZE2) :: ZSIGMA_RC ! Standard deviation of rc at time t
LOGICAL, DIMENSION(KPROMA) :: LLMICRO
!
!Output packed tendencies (for budgets only)
REAL, DIMENSION(KPROMA, IBUNUM-IBUNUM_EXTRA) :: ZBU_SUM
!
!For mixing-ratio-splitting
REAL, DIMENSION(KPROMA,0:7) :: ZVART !Packed variables
REAL, DIMENSION(KSIZE2,0:7) :: ZEXTPK   !To take into acount external tendencies inside the splitting
!
!For retroaction of E on IAGGS
REAL, DIMENSION(MERGE(KPROMA,0,OELEC)) :: ZLATHAM_IAGGS
!
INTEGER, DIMENSION(KPROMA) :: I1,I2 ! Used to replace the COUNT and PACK intrinsics on variables
INTEGER, DIMENSION(KSIZE) :: I1TOT, I2TOT ! Used to replace the COUNT and PACK intrinsics
!
INTEGER :: IC, JMICRO
LOGICAL :: LLSIGMA_RC, LL_AUCV_ADJU
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_PACK', 0, ZHOOK_HANDLE)
!
!*       1.     GENERALITIES
!               ------------
!
IKT=D%NKT
IKTB=D%NKTB
IKTE=D%NKTE
IIJT=D%NIJT
IIJB=D%NIJB
IIJE=D%NIJE
GEXT_TEND=.TRUE.
LLSIGMA_RC=(PARAMI%CSUBG_AUCV_RC=='PDF ' .AND. PARAMI%CSUBG_PR_PDF=='SIGM')
LL_AUCV_ADJU=(PARAMI%CSUBG_AUCV_RC=='ADJU' .OR. PARAMI%CSUBG_AUCV_RI=='ADJU')
!
IF(PARAMI%LPACK_MICRO) THEN
  IF(KPROMA /= KSIZE .AND. (PARAMI%CSUBG_RR_EVAP=='PRFR' .OR. PARAMI%CSUBG_RC_RR_ACCR=='PRFR')) THEN
    CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'RAIN_ICE', 'For now, KPROMA must be equal to KSIZE, see comments in code for explanation')
    ! Microphyscs was optimized by introducing chunks of KPROMA size
    ! Thus, in ice4_tendencies, the 1D array represent only a fraction of the points where microphysical species are present
    ! We cannot rebuild the entire 3D arrays in the subroutine, so we cannot call ice4_rainfr_vert in it
    ! A solution would be to suppress optimisation in this case by setting KPROMA=KSIZE in rain_ice
    ! Another solution would be to compute column by column?
    ! Another one would be to cut tendencies in 3 parts: before rainfr_vert, rainfr_vert, after rainfr_vert
  ENDIF
  !
  IF(BUCONF%LBU_ENABLE .OR. OSAVE_MICRO) THEN
    DO JV=1, IBUNUM-IBUNUM_EXTRA
      ZBU_PACK(:, JV)=0.
    ENDDO
  ENDIF
  !
  !*       2.     POINT SELECTION
  !               ---------------
  !
  !  optimization by looking for locations where
  !  the microphysical fields are larger than a minimal value only !!!
  !
  IF (KSIZE /= COUNT(ODMICRO(:,:))) THEN
      CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'ICE4_PACK', 'ICE4_PACK : KSIZE /= COUNT(ODMICRO)')
  ENDIF
  
  IF (KSIZE > 0) THEN
    !
    !*       3.     CACHE-BLOCKING LOOP
    !               -------------------
    !
  
    ! starting indexes :
    IC=0
    ISTK=IKTB
    ISTIJ=IIJB
  
    DO JMICRO=1,KSIZE,KPROMA
  
      IMICRO=MIN(KPROMA,KSIZE-JMICRO+1)
      !
      !*       4.     PACKING
      !               -------
      !
  
      ! Setup packing parameters
      OUTER_LOOP: DO JK = ISTK, IKTE
        IF (ANY(ODMICRO(:,JK))) THEN
          DO JIJ = ISTIJ, IIJE
            IF (ODMICRO(JIJ,JK)) THEN
              IC=IC+1
              LLMICRO(IC)=.TRUE.
              ! Initialization of variables in packed format :
              ZVART(IC, ITH)=PWR(JIJ, JK, ITH)
              ZVART(IC, IRV)=PWR(JIJ, JK, IRV)
              ZVART(IC, IRC)=PWR(JIJ, JK, IRC)
              ZVART(IC, IRR)=PWR(JIJ, JK, IRR)
              ZVART(IC, IRI)=PWR(JIJ, JK, IRI)
              ZVART(IC, IRS)=PWR(JIJ, JK, IRS)
              ZVART(IC, IRG)=PWR(JIJ, JK, IRG)
              IF (KRR==7) THEN
                ZVART(IC, IRH)=PWR(JIJ, JK, IRH)
              ENDIF
              IF (GEXT_TEND) THEN
                !The th tendency is not related to a mixing ratio change, there is no exn/exnref issue here
                ZEXTPK(IC, ITH)=PTHS(JIJ, JK)
                ZEXTPK(IC, IRV)=PRVS(JIJ, JK)
                ZEXTPK(IC, IRC)=PRCS(JIJ, JK)
                ZEXTPK(IC, IRR)=PRRS(JIJ, JK)
                ZEXTPK(IC, IRI)=PRIS(JIJ, JK)
                ZEXTPK(IC, IRS)=PRSS(JIJ, JK)
                ZEXTPK(IC, IRG)=PRGS(JIJ, JK)
                IF (KRR==7) THEN
                  ZEXTPK(IC, IRH)=PRHS(JIJ, JK)
                ENDIF
              ENDIF
              ZCIT       (IC)=PCIT    (JIJ, JK)
              ZCF        (IC)=PCLDFR  (JIJ, JK)
              ZRHODREF   (IC)=PRHODREF(JIJ, JK)
              ZPRES      (IC)=PPABST  (JIJ, JK)
              ZEXN       (IC)=PEXN    (JIJ, JK)
              IF(LLSIGMA_RC) THEN
                ZSIGMA_RC(IC)=PSIGS   (JIJ, JK)
              ENDIF
              IF (LL_AUCV_ADJU) THEN
                ZHLC_HCF(IC) = PHLC_HCF(JIJ, JK)
                ZHLC_HRC(IC) = PHLC_HRC(JIJ, JK)
                ZHLI_HCF(IC) = PHLI_HCF(JIJ, JK)
                ZHLI_HRI(IC) = PHLI_HRI(JIJ, JK)
              ENDIF
              IF (OELEC) ZLATHAM_IAGGS(IC) = PLATHAM_IAGGS(JIJ, JK)
              ! Save indices for later usages:
              I1(IC) = JIJ
              I2(IC) = JK
              IF(BUCONF%LBU_ENABLE .OR. OSAVE_MICRO) THEN
                I1TOT(JMICRO+IC-1)=JIJ
                I2TOT(JMICRO+IC-1)=JK
              ENDIF
              IF (IC==IMICRO) THEN
                ! the end of the chunk has been reached, then reset the starting index :
                ISTIJ=JIJ+1
                IF (ISTIJ <= IIJE) THEN
                  ISTK=JK
                ELSE
                  ! end of line, restart from 1 and increment upper loop
                  ISTIJ=D%NIJB
                  ISTK=JK+1
                  IF (ISTK > IKTE) THEN
                    ! end of line, restart from 1
                    ISTK=IKTB
                  ENDIF
                ENDIF
                IC=0
                EXIT OUTER_LOOP
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        ! restart inner loop on JIJ :
        ISTIJ=IIJB
      ENDDO OUTER_LOOP
  
      !
      !*       5.     TENDENCIES COMPUTATION
      !               ----------------------
      !
      CALL ICE4_STEPPING(D, CST, PARAMI, ICEP, ICED, BUCONF, &
                        &LLSIGMA_RC, LL_AUCV_ADJU, GEXT_TEND, &
                        &KPROMA, IMICRO, LLMICRO, PTSTEP, &
                        &KRR, OSAVE_MICRO, OELEC, &
                        &ZEXN, ZRHODREF, I1, I2, &
                        &ZPRES, ZCF, ZSIGMA_RC, &
                        &ZCIT, &
                        &ZVART, &
                        &ZHLC_HCF, ZHLC_HRC, &
                        &ZHLI_HCF, ZHLI_HRI, PRAINFR, &
                        &ZEXTPK, ZBU_SUM, ZRREVAV, &
                        &ZLATHAM_IAGGS)
      !
      !*       6.     UNPACKING
      !               ---------
      !
      DO JL=1, IMICRO
        PCIT  (I1(JL),I2(JL))=ZCIT   (JL)
        IF(PARAMI%LWARM) THEN
          PEVAP3D(I1(JL),I2(JL))=ZRREVAV(JL)
        ENDIF
        PWR(I1(JL),I2(JL),IRV)=ZVART(JL, IRV)
        PWR(I1(JL),I2(JL),IRC)=ZVART(JL, IRC)
        PWR(I1(JL),I2(JL),IRR)=ZVART(JL, IRR)
        PWR(I1(JL),I2(JL),IRI)=ZVART(JL, IRI)
        PWR(I1(JL),I2(JL),IRS)=ZVART(JL, IRS)
        PWR(I1(JL),I2(JL),IRG)=ZVART(JL, IRG)
        IF (KRR==7) THEN
          PWR(I1(JL),I2(JL),IRH)=ZVART(JL, IRH)
        ENDIF
      ENDDO
      IF(BUCONF%LBU_ENABLE .OR. OSAVE_MICRO) THEN
        DO JV=1, IBUNUM-IBUNUM_EXTRA
          DO JL=1, IMICRO
            ZBU_PACK(JMICRO+JL-1, JV) = ZBU_SUM(JL, JV)
          ENDDO
        ENDDO
      ENDIF
  
  
    ENDDO ! JMICRO
  ENDIF ! KSIZE > 0

ELSE ! PARAMI%LPACK_MICRO
  !We assume, here, that points outside the physical domain of the model (extral levels,
  !horizontal points in the halo) contain valid values, sufficiently valid to be used in tests
  !such as "PTHT(JL)>ZTHRESHOLD .AND. LLMICRO(JL)". In these tests, LLMICRO(JL) will be evaluated
  !to .FALSE. on these kind of points but valid values for PTHT are needed to prevent crash.
  !
  IF (KSIZE /= D%NIJT*D%NKT) THEN
      CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'ICE4_PACK', 'ICE4_PACK : KSIZE /= NIJT*NKT')
  ENDIF

  !Some arrays must be copied. In order not to waste memory, we re-use temporary arrays
  !declared for the pack case.
  IC=0
  DO JK = 1, IKT
    DO JIJ = 1, IIJT
      IC=IC+1
      I1TOT(IC)=JIJ
      I2TOT(IC)=JK
      IF (GEXT_TEND) THEN
        !The th tendency is not related to a mixing ratio change, there is no exn/exnref issue here
        ZEXTPK(IC, ITH)=PTHS(JIJ, JK)
        ZEXTPK(IC, IRV)=PRVS(JIJ, JK)
        ZEXTPK(IC, IRC)=PRCS(JIJ, JK)
        ZEXTPK(IC, IRR)=PRRS(JIJ, JK)
        ZEXTPK(IC, IRI)=PRIS(JIJ, JK)
        ZEXTPK(IC, IRS)=PRSS(JIJ, JK)
        ZEXTPK(IC, IRG)=PRGS(JIJ, JK)
        IF (KRR==7) THEN
          ZEXTPK(IC, IRH)=PRHS(JIJ, JK)
        ENDIF
      ENDIF
      IF(LLSIGMA_RC) THEN
        !Copy needed because sigma is modified in ice4_stepping
        ZSIGMA_RC(IC)=PSIGS(JIJ, JK)
      ENDIF
    ENDDO
  ENDDO
  !
  !When PARAMI%LPACK_MICRO=T, values on the extra levels are not given to ice4_stepping,
  !so there was not filled in rain_ice.
  !When PARAMI%LPACK_MICRO=F, we need to complement the work done in rain_ice to provide
  !valid values on these levels.
  !The same applies for the first points and last points on the horizontal dimension.
  IF (IKTB /= 1) THEN
    DO JK=1, IKTB-1
      PWR(:, JK, :)=PWR(:, IKTB, :)
    ENDDO
  ENDIF
  IF (IKTE /= IKT) THEN
    DO JK=IKTE+1, IKT
      PWR(:, JK, :)=PWR(:, IKTE, :)
    ENDDO
  ENDIF
  IF (IIJB /= 1) THEN
    DO JIJ=1, IIJB-1
      PWR(JIJ, :, :)=PWR(IIJB, :, :)
    ENDDO
  ENDIF
  IF (IIJE /= IIJT) THEN
    DO JIJ=IIJE+1, IIJT
      PWR(JIJ, :, :)=PWR(IIJE, :, :) 
    ENDDO
  ENDIF
  !
  !*       5bis.  TENDENCIES COMPUTATION
  !               ----------------------
  !
  CALL ICE4_STEPPING(D, CST, PARAMI, ICEP, ICED, BUCONF, &
                    &LLSIGMA_RC, LL_AUCV_ADJU, GEXT_TEND, &
                    &KSIZE, KSIZE, ODMICRO, PTSTEP, &
                    &KRR, OSAVE_MICRO, OELEC, &
                    &PEXN, PRHODREF, I1TOT, I2TOT, &
                    &PPABST, PCLDFR, ZSIGMA_RC, &
                    &PCIT, &
                    &PWR, &
                    &PHLC_HCF, PHLC_HRC, &
                    &PHLI_HCF, PHLI_HRI, PRAINFR, &
                    &ZEXTPK, ZBU_PACK, PEVAP3D, &
                    &ZLATHAM_IAGGS)

ENDIF ! PARAMI%LPACK_MICRO
!
!
!*       6.     SAVE MICROPHYSICAL TENDENCIES USED BY OTHER PHYSICAL PARAMETERIZATIONS
!               ----------------------------------------------------------------------
!
IF (OSAVE_MICRO) THEN
  DO JV = 1, IBUNUM-IBUNUM_EXTRA
    DO JL = 1, KSIZE 
      PMICRO_TEND(I1TOT(JL),I2TOT(JL),JV) = ZBU_PACK(JL,JV)
    ENDDO
  ENDDO        
END IF
!
!
!*       7.     BUDGETS
!               -------
!
IF(BUCONF%LBU_ENABLE) THEN
  !Budgets for the different processes
  CALL ICE4_BUDGETS(D, PARAMI, BUCONF, KSIZE, PTSTEP, KRR, I1TOT, I2TOT, &
                    PLVFACT, PLSFACT, PRHODJ, PEXNREF, &
                    PRVHENI, ZBU_PACK, &
                    TBUDGETS, KBUDGETS)
ENDIF
IF (LHOOK) CALL DR_HOOK('ICE4_PACK', 1, ZHOOK_HANDLE)
END SUBROUTINE ICE4_PACK
END MODULE MODE_ICE4_PACK
