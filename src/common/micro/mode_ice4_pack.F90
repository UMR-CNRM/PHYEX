!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_PACK
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_PACK(D, CST, PARAMI, ICEP, ICED, BUCONF, &
                    &KPROMA, KSIZE, PTSTEP, &
                    &KRR, OSAVE_MICRO, LDMICRO, OELEC, &
                    &PEXN, PRHODREF, PPABST, PCIT, PCLDFR, &
                    &PHLC_HCF, PHLC_HRC, PHLI_HCF, PHLI_HRI, &
                    &PTHS, PRS, PRREVAV, PRAINFR, PSIGS, PTHT, PRT, &
                    &PICLDFR, PZZZ, PCONC3D, PSSIO, PSSIU, PIFR, &
                    &PBUDGETS, PLATHAM_IAGGS)
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
!!     S. Riette Sept 23: all 3D arrays are suppressed from ice4_stepping
!  -----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_BUDGET,         ONLY: TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_FIELDS_ADDRESS, ONLY : & ! common fields adress
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
!
IMPLICIT NONE

!NOTES ON SIZES
!If we pack:
! - KSIZE is the number of relevant point (with mixing ratio different from 0)
! - KPROMA is the size of bloc of points
!If we do not pack:
! - KSIZE is the total number of points
! - KPROMA is null for memory saving
!
!When we do not pack, we can transmit directly the 3D arrays to the ice4_stepping subroutine, we do not need
!to copy the values. It is why KPROMA is null because we do not need these arrays.


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
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL,                  INTENT(IN)    :: OSAVE_MICRO  ! If true, microphysical tendencies are saved
LOGICAL,                  INTENT(IN)    :: OELEC        ! if true, cloud electricity is activated
LOGICAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: LDMICRO ! mask to limit computation
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLC_HCF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLC_HRC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLI_HCF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLI_HRI
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(D%NIJT,D%NKT,KRR),   INTENT(INOUT) :: PRS    ! m.r. source
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRREVAV! Rain evap profile
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRAINFR !Precipitation fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PICLDFR
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PZZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCONC3D
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSSIO   ! Super-saturation with respect to ice in the supersaturated fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSSIU   ! Sub-saturation with respect to ice in the  subsaturated fraction 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PIFR    ! Ratio cloud ice moist part to dry part 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHT
REAL, DIMENSION(D%NIJT,D%NKT,7), INTENT(INOUT) :: PRT
REAL, DIMENSION(MERGE(D%NIJT,0,OSAVE_MICRO .OR. BUCONF%LBU_ENABLE), &
                MERGE(D%NKT,0,OSAVE_MICRO .OR. BUCONF%LBU_ENABLE), &
                MERGE(IBUNUM-IBUNUM_EXTRA,0,OSAVE_MICRO .OR. BUCONF%LBU_ENABLE)), &
                                          INTENT(OUT) :: PBUDGETS  ! Microphysical tendencies
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), &
                                          INTENT(IN)    :: PLATHAM_IAGGS  ! E Function to simulate
                                                                          ! enhancement of IAGGS
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
!For packing
INTEGER :: IMICRO ! Case r_x>0 locations
INTEGER :: JL, JV
REAL, DIMENSION(KPROMA) :: &
                        & ZCIT,     & ! Pristine ice conc. at t
                        & ZICLDFR,  & ! Ice cloud fraction
                        & ZZZZ,     & ! model level height
                        & ZCONC3D,  & ! Cloud droplet number concentration at t
                        & ZSSIO,    & ! Super-saturation with respect to ice in the supersaturated fraction
                        & ZSSIU,    & ! Sub-saturation with respect to ice in the  subsaturated fraction
                        & ZIFR,     & ! Ratio cloud ice moist part to dry part
                        & ZRHODREF, & ! RHO Dry REFerence
                        & ZPABST,   & ! Pressure
                        & ZEXN,     & ! EXNer Pressure
                        & ZCLDFR,   & ! Cloud fraction
                        & ZHLC_HCF, & ! HLCLOUDS : fraction of High Cloud Fraction in grid
                        & ZHLC_HRC, & ! HLCLOUDS : LWC that is High LWC in grid
                        & ZHLI_HCF, &
                        & ZHLI_HRI, &
                        & ZRAINFR,  &
                        & ZRREVAV,  &
                        & ZSIGS,    & ! Standard deviation of rc at time t
                        & ZTHT,     &
                        & ZTHS
LOGICAL, DIMENSION(KPROMA) :: LLMICRO
!
!Output packed tendencies (for budgets only)
REAL, DIMENSION(KPROMA, IBUNUM-IBUNUM_EXTRA) :: ZBUDGETS
!
REAL, DIMENSION(KPROMA,7) :: ZRT !Packed variables
REAL, DIMENSION(KPROMA,7) :: ZRS !To take into acount external tendencies inside the splitting
!
!For retroaction of E on IAGGS
REAL, DIMENSION(MERGE(KPROMA,0,OELEC)) :: ZLATHAM_IAGGS
!
INTEGER, DIMENSION(KPROMA) :: I1,I2 ! Used to replace the COUNT and PACK intrinsics on variables
!
INTEGER :: IC, JMICRO, IDX
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
LLSIGMA_RC=(PARAMI%CSUBG_AUCV_RC=='PDF ' .AND. PARAMI%CSUBG_PR_PDF=='SIGM')
LL_AUCV_ADJU=(PARAMI%CSUBG_AUCV_RC=='ADJU' .OR. PARAMI%CSUBG_AUCV_RI=='ADJU')
!
IF(PARAMI%LPACK_MICRO) THEN
  !
  !*       2.     POINT SELECTION
  !               ---------------
  !
  !  optimization by looking for locations where
  !  the microphysical fields are larger than a minimal value only !!!
  !
  IF (KSIZE /= COUNT(LDMICRO(IIJB:IIJE,IKTB:IKTE))) THEN
      CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'ICE4_PACK', 'ICE4_PACK : KSIZE /= COUNT(LDMICRO)')
  ENDIF

  PBUDGETS(:,:,:)=0.
  
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
!$acc kernels
!$acc loop seq
      OUTER_LOOP: DO JK = ISTK, IKTE
        IF (ANY(LDMICRO(IIJB:IIJE,JK))) THEN
          !$acc loop gang vector independent
          DO JIJ = ISTIJ, IIJE
            IF (LDMICRO(JIJ,JK)) THEN
              !$acc atomic capture
              IC=IC+1
              IDX=IC !change of variable to use acc atomic capture (IC is shared)
              !$acc end atomic
              LLMICRO(IDX)=.TRUE.
              ! Initialization of variables in packed format :
              ZTHT(IDX)=PTHT(JIJ, JK)
              ZRT(IDX, IRV)=PRT(JIJ, JK, IRV)
              ZRT(IDX, IRC)=PRT(JIJ, JK, IRC)
              ZRT(IDX, IRR)=PRT(JIJ, JK, IRR)
              ZRT(IDX, IRI)=PRT(JIJ, JK, IRI)
              ZRT(IDX, IRS)=PRT(JIJ, JK, IRS)
              ZRT(IDX, IRG)=PRT(JIJ, JK, IRG)
              IF (KRR==7) THEN
                ZRT(IDX, IRH)=PRT(JIJ, JK, IRH)
              ENDIF
              IF (PARAMI%LEXT_TEND) THEN
                !The th tendency is not related to a mixing ratio change, there is no exn/exnref issue here
                ZTHS(IDX)=PTHS(JIJ, JK)
                ZRS(IDX, IRV)=PRS(JIJ, JK, IRV)
                ZRS(IDX, IRC)=PRS(JIJ, JK, IRC)
                ZRS(IDX, IRR)=PRS(JIJ, JK, IRR)
                ZRS(IDX, IRI)=PRS(JIJ, JK, IRI)
                ZRS(IDX, IRS)=PRS(JIJ, JK, IRS)
                ZRS(IDX, IRG)=PRS(JIJ, JK, IRG)
                IF (KRR==7) THEN
                  ZRS(IDX, IRH)=PRS(JIJ, JK, IRH)
                ENDIF
              ENDIF
              ZCIT       (IDX)=PCIT    (JIJ, JK)
              ZCLDFR     (IDX)=PCLDFR  (JIJ, JK)
              ZRHODREF   (IDX)=PRHODREF(JIJ, JK)
              ZPABST     (IDX)=PPABST  (JIJ, JK)
              ZEXN       (IDX)=PEXN    (JIJ, JK)
              ZICLDFR    (IC)=PICLDFR (JIJ, JK)
              ZZZZ       (IC)=PZZZ    (JIJ, JK)
              ZCONC3D    (IC)=PCONC3D (JIJ, JK)
              ZSSIO      (IC)=PSSIO   (JIJ, JK)
              ZSSIU      (IC)=PSSIU   (JIJ, JK)
              ZIFR       (IC)=PIFR    (JIJ, JK)
              IF(LLSIGMA_RC) THEN
                ZSIGS(IDX)    =PSIGS   (JIJ, JK)
              ENDIF
              IF (LL_AUCV_ADJU) THEN
                ZHLC_HCF(IDX) = PHLC_HCF(JIJ, JK)
                ZHLC_HRC(IDX) = PHLC_HRC(JIJ, JK)
                ZHLI_HCF(IDX) = PHLI_HCF(JIJ, JK)
                ZHLI_HRI(IDX) = PHLI_HRI(JIJ, JK)
              ENDIF
              ZRAINFR(IC)=PRAINFR(JIJ, JK)
              IF (OELEC) ZLATHAM_IAGGS(IC) = PLATHAM_IAGGS(JIJ, JK)
              ! Save indices for later usages:
              I1(IDX) = JIJ
              I2(IDX) = JK
              IF (IC==IMICRO) THEN
                ! the end of the chunk has been reached, then reset the starting index :
                ISTIJ=JIJ+1
                IF (ISTIJ <= IIJE) THEN
                  !$acc atomic write
                  ISTK=JK
                  !$acc end atomic
                ELSE
                  ! end of line, restart from 1 and increment upper loop
                  ISTIJ=D%NIJB
                  !$acc atomic write
                  ISTK=JK+1
                  !$acc end atomic
                  IF (ISTK > IKTE) THEN
                    ! end of line, restart from 1
                    !$acc atomic write
                    ISTK=IKTB
                    !$acc end atomic
                  ENDIF
                ENDIF
#ifndef MNH_OPENACC
                IC=0
                EXIT OUTER_LOOP
#endif
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        ! restart inner loop on JIJ :
        ISTIJ=IIJB
      ENDDO OUTER_LOOP
  !$acc end kernels
      !
      !*       5.     TENDENCIES COMPUTATION
      !               ----------------------
      !
      CALL ICE4_STEPPING(CST, PARAMI, ICEP, ICED, BUCONF, &
                        &KPROMA, IMICRO, PTSTEP, &
                        &KRR, OSAVE_MICRO, LLMICRO, OELEC, &
                        &ZEXN, ZRHODREF, &
                        &ZPABST, ZCIT, ZCLDFR, &
                        &ZHLC_HCF, ZHLC_HRC, &
                        &ZHLI_HCF, ZHLI_HRI, &
                        &ZTHS, ZRS, ZRREVAV, &
                        &ZRAINFR, ZSIGS, &
                        &ZTHT, ZRT, &
                        &ZICLDFR, ZZZZ, ZCONC3D, &
                        &ZSSIO, ZSSIU, ZIFR, &
                        &ZBUDGETS, &
                        &ZLATHAM_IAGGS)
      !
      !*       6.     UNPACKING
      !               ---------
      !
!$acc kernels
!$acc loop independent
      DO JL=1, IMICRO
        PCIT  (I1(JL),I2(JL))=ZCIT   (JL)
        PRREVAV(I1(JL),I2(JL))=ZRREVAV(JL)
        PRT(I1(JL),I2(JL),IRV)=ZRT(JL, IRV)
        PRT(I1(JL),I2(JL),IRC)=ZRT(JL, IRC)
        PRT(I1(JL),I2(JL),IRR)=ZRT(JL, IRR)
        PRT(I1(JL),I2(JL),IRI)=ZRT(JL, IRI)
        PRT(I1(JL),I2(JL),IRS)=ZRT(JL, IRS)
        PRT(I1(JL),I2(JL),IRG)=ZRT(JL, IRG)
        IF (KRR==7) THEN
          PRT(I1(JL),I2(JL),IRH)=ZRT(JL, IRH)
        ENDIF
        PRAINFR(I1(JL),I2(JL))=ZRAINFR(JL)
      ENDDO
!$acc end kernels 
      IF(BUCONF%LBU_ENABLE .OR. OSAVE_MICRO) THEN
!$acc kernels
!$acc loop independent collapse(2)
        DO JV=1, IBUNUM-IBUNUM_EXTRA
          DO JL=1, IMICRO
            PBUDGETS(I1(JL),I2(JL),JV)=ZBUDGETS(JL, JV)
          ENDDO
        ENDDO
!$acc end kernels
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
  !
  !When PARAMI%LPACK_MICRO=T, values on the extra levels are not given to ice4_stepping,
  !so there was not filled in rain_ice.
  !When PARAMI%LPACK_MICRO=F, we need to complement the work done in rain_ice to provide
  !valid values on these levels.
  !The same applies for the first points and last points on the horizontal dimension.
  IF (IKTB /= 1) THEN
    DO JK=1, IKTB-1
      PRT(:, JK, :)=PRT(:, IKTB, :)
    ENDDO
  ENDIF
  IF (IKTE /= IKT) THEN
    DO JK=IKTE+1, IKT
      PRT(:, JK, :)=PRT(:, IKTE, :)
    ENDDO
  ENDIF
  IF (IIJB /= 1) THEN
    DO JIJ=1, IIJB-1
      PRT(JIJ, :, :)=PRT(IIJB, :, :)
    ENDDO
  ENDIF
  IF (IIJE /= IIJT) THEN
    DO JIJ=IIJE+1, IIJT
      PRT(JIJ, :, :)=PRT(IIJE, :, :) 
    ENDDO
  ENDIF
  !
  !*       5bis.  TENDENCIES COMPUTATION
  !               ----------------------
  !
  CALL ICE4_STEPPING(CST, PARAMI, ICEP, ICED, BUCONF, &
                    &KSIZE, D%NIJB, PTSTEP, &
                    &KRR, OSAVE_MICRO, LDMICRO, OELEC, &
                    &PEXN, PRHODREF, &
                    &PPABST, PCIT, PCLDFR, &
                    &PHLC_HCF, PHLC_HRC, &
                    &PHLI_HCF, PHLI_HRI,  &
                    &PTHS, PRS, PRREVAV, &
                    &PRAINFR, PSIGS, &
                    &PTHT, PRT, &
                    &PICLDFR, PZZZ, PCONC3D, &
                    &PSSIO, PSSIU, PIFR, &
                    &PBUDGETS, &
                    &PLATHAM_IAGGS)

ENDIF ! PARAMI%LPACK_MICRO
!
IF (LHOOK) CALL DR_HOOK('ICE4_PACK', 1, ZHOOK_HANDLE)
END SUBROUTINE ICE4_PACK
END MODULE MODE_ICE4_PACK
