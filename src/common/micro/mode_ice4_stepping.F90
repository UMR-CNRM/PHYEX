!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_STEPPING
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_STEPPING(D, CST, PARAMI, ICEP, ICED, BUCONF, &
                        &LDSIGMA_RC, LDAUCV_ADJU, LDEXT_TEND, &
                        &KPROMA, KMICRO, LDMICRO, PTSTEP, &
                        &KRR, OSAVE_MICRO, OELEC, &
                        &PEXN, PRHODREF, K1, K2, &
                        &PPRES, PCF, PSIGMA_RC, &
                        &PCIT, &
                        &PVART, &
                        &PHLC_HCF, PHLC_HRC, &
                        &PHLI_HCF, PHLI_HRI, PRAINFR, &
                        &PEXTPK, PBU_SUM, PRREVAV, &
                        &PLATHAM_IAGGS)
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
!!     R. El Khatib 03-May-2023 Replace OMP SIMD loops by explicit loops : more portable and even slightly faster
!!     S. Riette Sept 23: 3D arrays suppressed from ice4_tendencies
!  -----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t

USE MODD_BUDGET,         ONLY: TBUDGETCONF_t
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
      & IBUNUM,       & ! Number of tendency terms
      & IBUNUM_MR,    & ! Number of tendency terms expressed as mixing ratio changes
      & IBUNUM_EXTRA, & ! Number of extra tendency terms
      & IRREVAV,      & ! Index for the evaporation tendency
      & IBUEXTRAIND     ! Index indirection

USE MODE_ICE4_TENDENCIES, ONLY: ICE4_TENDENCIES
!
IMPLICIT NONE
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
LOGICAL,                  INTENT(IN)    :: LDSIGMA_RC
LOGICAL,                  INTENT(IN)    :: LDAUCV_ADJU
LOGICAL,                  INTENT(IN)    :: LDEXT_TEND
INTEGER,                  INTENT(IN)    :: KPROMA ! cache-blocking factor for microphysic loop
INTEGER,                  INTENT(IN)    :: KMICRO ! Case r_x>0 locations
LOGICAL, DIMENSION(KPROMA), INTENT(IN)  :: LDMICRO
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL,                  INTENT(IN)    :: OSAVE_MICRO   ! if true, save the microphysical tendencies
LOGICAL,                  INTENT(IN)    :: OELEC         ! if true, cloud electricity is activated
!
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PEXN    ! Exner function
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PRHODREF! Reference density
INTEGER, DIMENSION(KPROMA),                     INTENT(IN)    :: K1,K2 ! Used to replace the COUNT and PACK intrinsics on variables
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PPRES
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PCF ! Cloud fraction
REAL,    DIMENSION(KPROMA),                     INTENT(INOUT) :: PSIGMA_RC
REAL,    DIMENSION(KPROMA),                     INTENT(INOUT) :: PCIT
REAL,    DIMENSION(KPROMA,0:7),                 INTENT(INOUT) :: PVART !Packed variables
REAL,    DIMENSION(KPROMA),                     INTENT(INOUT) :: PHLC_HRC
REAL,    DIMENSION(KPROMA),                     INTENT(INOUT) :: PHLC_HCF
REAL,    DIMENSION(KPROMA),                     INTENT(INOUT) :: PHLI_HRI
REAL,    DIMENSION(KPROMA),                     INTENT(INOUT) :: PHLI_HCF
REAL,    DIMENSION(D%NIJT,D%NKT),               INTENT(INOUT) :: PRAINFR
REAL,    DIMENSION(KPROMA,0:7),                 INTENT(INOUT) :: PEXTPK !To take into acount external tendencies inside the splitting
REAL,    DIMENSION(KPROMA, IBUNUM-IBUNUM_EXTRA),INTENT(OUT)   :: PBU_SUM
REAL,    DIMENSION(KPROMA),                     INTENT(OUT)   :: PRREVAV
REAL,    DIMENSION(MERGE(KPROMA,0,OELEC)),      INTENT(IN)    :: PLATHAM_IAGGS ! E Function to simulate
                                                                               ! enhancement of IAGGS
!
!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
LOGICAL :: LSOFT ! Must we really compute tendencies or only adjust them to new T variables
INTEGER :: INB_ITER_MAX ! Maximum number of iterations (with real tendencies computation)
REAL :: ZTSTEP ! length of sub-timestep in case of time splitting
REAL :: ZINV_TSTEP ! Inverse ov PTSTEP
REAL :: ZTIME_THRESHOLD ! Time to reach threshold
!
INTEGER :: JL, JV, JJV
REAL, DIMENSION(KPROMA) :: &
                        & ZTIME, & ! Current integration time (starts with 0 and ends with PTSTEP)
                        & ZMAXTIME, & ! Time on which we can apply the current tendencies
                        & ZTIME_LASTCALL, &     ! Integration time when last tendecies call has been done
                        & ZSSI,     &
                        & ZZT,      & ! Temperature
                        & ZLSFACT,  & ! L_s/(Pi*C_ph)
                        & ZLVFACT,  & ! L_v/(Pi*C_ph)
                        & ZHLC_LCF, & ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                      !    note that PCF = PHLC_HCF + ZHLC_LCF
                        & ZHLC_LRC, & ! HLCLOUDS : LWC that is Low  LWC in grid
                                      !    note that ZRC = PHLC_HRC + ZHLC_LRC
                        & ZHLI_LCF, &
                        & ZHLI_LRI
LOGICAL, DIMENSION(KPROMA) :: LLCOMPUTE ! .TRUE. or points where we must compute tendencies,
!
!Output packed tendencies (for budgets only)
REAL, DIMENSION(KPROMA, IBUNUM) :: ZBU_INST
!
!For mixing-ratio-splitting
LOGICAL :: LLCPZ0RT
REAL :: ZTIME_THRESHOLD1D(KPROMA) ! Time to reach threshold
REAL, DIMENSION(KPROMA, KRR) :: Z0RT ! Mixing-ratios at the beginig of the current loop
!
REAL, DIMENSION(KPROMA,0:7) :: ZA, ZB
!
REAL, DIMENSION(KPROMA, 8) :: ZRS_TEND, ZRG_TEND
REAL, DIMENSION(KPROMA,10) :: ZRH_TEND

INTEGER, DIMENSION(KPROMA) :: IITER    ! Number of iterations done (with real tendencies computation)
!
REAL, DIMENSION(KPROMA) :: ZSUM2, ZMAXB
REAL :: ZDEVIDE, ZX
!
LOGICAL :: LL_ANY_ITER
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_STEPPING', 0, ZHOOK_HANDLE)
!
!*       1.     GENERALITIES
!               ------------
!
ZINV_TSTEP=1./PTSTEP
!
IF(BUCONF%LBU_ENABLE .OR. OSAVE_MICRO) THEN
  DO JV=1, IBUNUM-IBUNUM_EXTRA
    PBU_SUM(:, JV)=0.
  ENDDO
ENDIF

!Maximum number of iterations
!We only count real iterations (those for which we *compute* tendencies)
INB_ITER_MAX=PARAMI%NMAXITER_MICRO
IF(PARAMI%XTSTEP_TS/=0.)THEN
  INB_ITER_MAX=MAX(1, INT(PTSTEP/PARAMI%XTSTEP_TS)) !At least the number of iterations needed for the time-splitting
  ZTSTEP=PTSTEP/INB_ITER_MAX
  INB_ITER_MAX=MAX(PARAMI%NMAXITER_MICRO, INB_ITER_MAX) !For the case XMRSTEP/=0. at the same time
ENDIF

IF (LDEXT_TEND) THEN
  DO JV=0, KRR
    DO JL=1, KMICRO
      PEXTPK(JL, JV)=PEXTPK(JL, JV)-PVART(JL, JV)*ZINV_TSTEP
    ENDDO
  ENDDO
ENDIF
IF (LDSIGMA_RC) THEN
  DO JL=1, KMICRO
    PSIGMA_RC(JL)=PSIGMA_RC(JL)*2.
  ENDDO 
ENDIF
IF (LDAUCV_ADJU) THEN
  DO JL=1, KMICRO
    ZHLC_LRC(JL) = PVART(JL, IRC) - PHLC_HRC(JL)
    ZHLI_LRI(JL) = PVART(JL, IRI) - PHLI_HRI(JL)
    IF(PVART(JL, IRC)>0.) THEN
      ZHLC_LCF(JL) = PCF(JL)- PHLC_HCF(JL)
    ELSE
      ZHLC_LCF(JL)=0.
    ENDIF
    IF(PVART(JL, IRI)>0.) THEN
      ZHLI_LCF(JL) = PCF(JL)- PHLI_HCF(JL)
    ELSE
      ZHLI_LCF(JL)=0.
    ENDIF
  ENDDO
ENDIF

!-------------------------------------------------------------------------------
!
!***       4.4 temporal loop
!
!
IITER(1:KMICRO)=0
DO JL=1, KMICRO
  IF(LDMICRO(JL)) THEN
    ZTIME(JL)=0. ! Current integration time (all points may have a different integration time)
  ELSE
    ZTIME(JL)=PTSTEP ! Nothing to do on this point, it has already reached the end of the timestep
  ENDIF
ENDDO

DO WHILE(ANY(ZTIME(1:KMICRO)<PTSTEP)) ! Loop to *really* compute tendencies

  IF(PARAMI%XTSTEP_TS/=0.) THEN
    ! In this case we need to remember the time when tendencies were computed
    ! because when time has evolved more than a limit, we must re-compute tendencies
    ZTIME_LASTCALL(1:KMICRO)=ZTIME(1:KMICRO)
  ENDIF
  DO JL=1, KMICRO
    IF (ZTIME(JL) < PTSTEP) THEN
      LLCOMPUTE(JL)=.TRUE. ! Computation (.TRUE.) only for points for which integration time has not reached the timestep
      IITER(JL)=IITER(JL)+1
    ELSE
      LLCOMPUTE(JL)=.FALSE.
    ENDIF
  ENDDO
  LL_ANY_ITER=ANY(IITER(1:KMICRO) < INB_ITER_MAX)
  LLCPZ0RT=.TRUE.
  LSOFT=.FALSE. ! We *really* compute the tendencies

  DO WHILE(ANY(LLCOMPUTE(1:KMICRO))) ! Loop to adjust tendencies when we cross the 0°C or when a species disappears
    ZSUM2(1:KMICRO)=PVART(1:KMICRO, IRI)
    DO JV=IRI+1,KRR
      DO JL=1, KMICRO
        ZSUM2(JL)=ZSUM2(JL)+PVART(JL, JV)
      ENDDO
    ENDDO
    DO JL=1, KMICRO
      ZDEVIDE=(CST%XCPD + CST%XCPV*PVART(JL, IRV) + CST%XCL*(PVART(JL, IRC)+PVART(JL, IRR)) + CST%XCI*ZSUM2(JL)) * PEXN(JL)
      ZZT(JL) = PVART(JL, ITH) * PEXN(JL)
      ZLSFACT(JL)=(CST%XLSTT+(CST%XCPV-CST%XCI)*(ZZT(JL)-CST%XTT)) / ZDEVIDE
      ZLVFACT(JL)=(CST%XLVTT+(CST%XCPV-CST%XCL)*(ZZT(JL)-CST%XTT)) / ZDEVIDE
    ENDDO
    !-------------------------------------------------------------------------------
    !
    !***       4.5 Effective tendencies computation
    !
    !
    ! Tendencies are *really* computed when LSOFT==.FALSE. and only adjusted otherwise
    CALL ICE4_TENDENCIES(CST, PARAMI, ICEP, ICED, BUCONF, &
                        &KPROMA, KMICRO, &
                        &KRR, LSOFT, LLCOMPUTE, &
                        &OSAVE_MICRO, OELEC, &
                        &PEXN, PRHODREF, ZLVFACT, ZLSFACT, K1, K2, &
                        &PPRES, PCF, PSIGMA_RC, &
                        &PCIT, &
                        &ZZT, PVART, &
                        &PLATHAM_IAGGS, &
                        &ZBU_INST, &
                        &ZRS_TEND, ZRG_TEND, ZRH_TEND, ZSSI, &
                        &ZA, ZB, &
                        &PHLC_HCF, ZHLC_LCF, PHLC_HRC, ZHLC_LRC, &
                        &PHLI_HCF, ZHLI_LCF, PHLI_HRI, ZHLI_LRI, PRAINFR)

    ! External tendencies
    IF(LDEXT_TEND) THEN
      DO JV=0, KRR
        DO JL=1, KMICRO
          ZA(JL, JV) = ZA(JL, JV) + PEXTPK(JL, JV)
        ENDDO
      ENDDO
    ENDIF
    !-------------------------------------------------------------------------------
    !
    !***       4.6 Time integration
    !
    !
    ! If we can, we shall use these tendencies until the end of the timestep
    DO JL=1, KMICRO
      IF(LLCOMPUTE(JL)) THEN
        ZMAXTIME(JL)=(PTSTEP-ZTIME(JL)) ! Remaining time until the end of the timestep
      ELSE
        ZMAXTIME(JL)=0.
      ENDIF
    ENDDO

    !We need to adjust tendencies when temperature reaches 0
    IF(PARAMI%LFEEDBACKT) THEN
      DO JL=1, KMICRO
        !Is ZB(:, ITH) enough to change temperature sign?
        ZX=CST%XTT/PEXN(JL)
        IF ((PVART(JL, ITH) - ZX) * (PVART(JL, ITH) + ZB(JL, ITH) - ZX) < 0.) THEN
          ZMAXTIME(JL)=0.
        ENDIF
        !Can ZA(:, ITH) make temperature change of sign?
        IF (ABS(ZA(JL,ITH)) > 1.E-20 ) THEN
          ZTIME_THRESHOLD=(ZX - ZB(JL, ITH) - PVART(JL, ITH))/ZA(JL, ITH)
          IF (ZTIME_THRESHOLD > 0.) THEN
            ZMAXTIME(JL)=MIN(ZMAXTIME(JL), ZTIME_THRESHOLD)
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    !We need to adjust tendencies when a species disappears
    !When a species is missing, only the external tendencies can be negative (and we must keep track of it)
    DO JV=1, KRR
      DO JL=1, KMICRO
        IF (ZA(JL, JV) < -1.E-20 .AND. PVART(JL, JV) > ICED%XRTMIN(JV)) THEN
          ZMAXTIME(JL)=MIN(ZMAXTIME(JL), -(ZB(JL, JV)+PVART(JL, JV))/ZA(JL, JV))
          ZMAXTIME(JL)=MAX(ZMAXTIME(JL), CST%XMNH_TINY) !to prevent rounding errors
        ENDIF
      ENDDO
    ENDDO

    !We stop when the end of the timestep is reached
    DO JL=1, KMICRO
      IF (ZTIME(JL)+ZMAXTIME(JL) >= PTSTEP) THEN
        LLCOMPUTE(JL)=.FALSE.
      ENDIF
    ENDDO
    !We must recompute tendencies when the end of the sub-timestep is reached
    IF (PARAMI%XTSTEP_TS/=0.) THEN
      DO JL=1, KMICRO
        IF ((IITER(JL) < INB_ITER_MAX) .AND. (ZTIME(JL)+ZMAXTIME(JL) > ZTIME_LASTCALL(JL)+ZTSTEP)) THEN
          ZMAXTIME(JL)=ZTIME_LASTCALL(JL)-ZTIME(JL)+ZTSTEP
          LLCOMPUTE(JL)=.FALSE.
        ENDIF
      ENDDO
    ENDIF

    !We must recompute tendencies when the maximum allowed change is reached
    !When a species is missing, only the external tendencies can be active and we do not want to recompute
    !the microphysical tendencies when external tendencies are negative (results won't change because species was already missing)
    IF (PARAMI%XMRSTEP/=0.) THEN
      IF (LL_ANY_ITER) THEN
        ! In this case we need to remember the initial mixing ratios used to compute the tendencies
        ! because when mixing ratio has evolved more than a threshold, we must re-compute tendencies
        ! Thus, at first iteration (ie when LLCPZ0RT=.TRUE.) we copy PVART into Z0RT
        DO JV=1,KRR
          IF (LLCPZ0RT) THEN
            Z0RT(1:KMICRO, JV)=PVART(1:KMICRO, JV)
          ENDIF
          DO JL=1, KMICRO
            IF (IITER(JL)<INB_ITER_MAX .AND. ABS(ZA(JL,JV))>1.E-20) THEN
              ZTIME_THRESHOLD1D(JL)=(SIGN(1., ZA(JL, JV))*PARAMI%XMRSTEP+ &
                                    &Z0RT(JL, JV)-PVART(JL, JV)-ZB(JL, JV))/ZA(JL, JV)
            ELSE
              ZTIME_THRESHOLD1D(JL)=-1.
            ENDIF
          ENDDO
          DO JL=1, KMICRO
            IF (ZTIME_THRESHOLD1D(JL)>=0 .AND. ZTIME_THRESHOLD1D(JL)<ZMAXTIME(JL) .AND. &
               &(PVART(JL, JV)>ICED%XRTMIN(JV) .OR. ZA(JL, JV)>0.)) THEN
              ZMAXTIME(JL)=MIN(ZMAXTIME(JL), ZTIME_THRESHOLD1D(JL))
              LLCOMPUTE(JL)=.FALSE.
            ENDIF
          ENDDO
          IF (JV == 1) THEN
            DO JL=1, KMICRO
              ZMAXB(JL)=ABS(ZB(JL, JV))
            ENDDO
          ELSE
            DO JL=1, KMICRO
              ZMAXB(JL)=MAX(ZMAXB(JL), ABS(ZB(JL, JV)))
            ENDDO
          ENDIF
        ENDDO
        LLCPZ0RT=.FALSE.
        DO JL=1, KMICRO
          IF (IITER(JL)<INB_ITER_MAX .AND. ZMAXB(JL)>PARAMI%XMRSTEP) THEN
            ZMAXTIME(JL)=0.
            LLCOMPUTE(JL)=.FALSE.
          ENDIF
        ENDDO
      ENDIF ! LL_ANY_ITER
    ENDIF ! XMRSTEP/=0.
    !-------------------------------------------------------------------------------
    !
    !***       4.7 New values of variables for next iteration
    !
    !
    DO JV=0, KRR
      DO JL=1, KMICRO
        IF(LDMICRO(JL)) THEN
          PVART(JL, JV)=PVART(JL, JV)+ZA(JL, JV)*ZMAXTIME(JL)+ZB(JL, JV)
        ENDIF
      ENDDO
    ENDDO
    DO JL=1, KMICRO
      IF (PVART(JL,IRI)<=0. .AND. LDMICRO(JL)) PCIT(JL) = 0.
      ZTIME(JL)=ZTIME(JL)+ZMAXTIME(JL)
    ENDDO
    !-------------------------------------------------------------------------------
    !
    !***       4.8 Mixing ratio change due to each process
    !
    IF(BUCONF%LBU_ENABLE .OR. OSAVE_MICRO) THEN
      !Mixing ratio change due to a tendency
      DO JV=1, IBUNUM-IBUNUM_MR-IBUNUM_EXTRA
        DO JL=1, KMICRO
          PBU_SUM(JL, JV) = PBU_SUM(JL, JV) + ZBU_INST(JL, JV)*ZMAXTIME(JL)
        ENDDO
      ENDDO

      !Mixing ratio change due to a mixing ratio change
      DO JV=IBUNUM-IBUNUM_MR-IBUNUM_EXTRA+1, IBUNUM-IBUNUM_EXTRA
        DO JL=1, KMICRO
          PBU_SUM(JL, JV) = PBU_SUM(JL, JV) + ZBU_INST(JL, JV)
        ENDDO
      ENDDO

      !Extra contribution as a mixing ratio change
      DO JV=IBUNUM-IBUNUM_EXTRA+1, IBUNUM
        JJV=IBUEXTRAIND(JV)
        DO JL=1, KMICRO
          PBU_SUM(JL, JJV) = PBU_SUM(JL, JJV) + ZBU_INST(JL, JV)
        ENDDO
      ENDDO
    ENDIF
    !-------------------------------------------------------------------------------
    !
    !***       4.9 Next loop
    !
    LSOFT=.TRUE. ! We try to adjust tendencies (inner while loop)
  ENDDO !Iterations on tendency computations (WHILE ANY(LLCOMPUTE))
ENDDO !Temporal loop

IF(LDEXT_TEND) THEN
  !Z..T variables contain the external tendency, we substract it
  DO JV=0, KRR
    DO JL=1, KMICRO
      IF(LDMICRO(JL)) THEN
        PVART(JL, JV) = PVART(JL, JV) - PEXTPK(JL, JV) * PTSTEP
      ENDIF
    ENDDO
  ENDDO
ENDIF
DO JL=1, KMICRO
  PRREVAV(JL)=ZBU_INST(JL, IRREVAV)
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ICE4_STEPPING', 1, ZHOOK_HANDLE)
END SUBROUTINE ICE4_STEPPING
END MODULE MODE_ICE4_STEPPING
