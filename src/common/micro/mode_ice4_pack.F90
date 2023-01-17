!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_PACK
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_PACK(D, CST, PARAMI, ICEP, ICED, BUCONF,                   &
                     KPROMA, KSIZE,                                        &
                     HSUBG_AUCV_RC, HSUBG_AUCV_RI,                         &
                     PTSTEP, KRR, ODMICRO, PEXN,                           &
                     PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,      &
                     PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,               &
                     PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                   &
                     PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                     PEVAP3D,                                              &
                     PRAINFR, PSIGS,                                       &
                     PRVHENI, PLVFACT, PLSFACT,                            &
                     PWR,                                                  &
                     TBUDGETS, KBUDGETS,                                   &
                     PRHT, PRHS                                            )
!  -----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_BUDGET,         ONLY: TBUDGETDATA, TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
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
      & IBUNUM_MR,    & ! Number of tendency terms expressed as mixing ratio changes
      & IBUNUM_EXTRA, & ! Number of extra tendency terms
      & IRREVAV,      & ! Index for the evaporation tendency
      & IBUEXTRAIND     ! Index indirection

USE MODE_MSG,            ONLY: PRINT_MSG, NVERB_FATAL

USE MODE_ICE4_TENDENCIES, ONLY: ICE4_TENDENCIES
USE MODE_ICE4_BUDGETS, ONLY: ICE4_BUDGETS
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
INTEGER,                  INTENT(IN)    :: KPROMA ! cache-blocking factor for microphysic loop
INTEGER,                  INTENT(IN)    :: KSIZE
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Kind of Subgrid autoconversion method
CHARACTER(LEN=80),        INTENT(IN)    :: HSUBG_AUCV_RI ! Kind of Subgrid autoconversion method
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
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
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PHLC_HRC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PHLC_HCF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PHLI_HRI
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PHLI_HCF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: PEVAP3D! Rain evap profile
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)   :: PRAINFR !Precipitation fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRVHENI ! heterogeneous nucleation
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PLVFACT
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PLSFACT
REAL, DIMENSION(D%NIJT,D%NKT,0:7), INTENT(OUT) :: PWR
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL,  INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL,  INTENT(INOUT) :: PRHS    ! Hail m.r. source
!
!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
INTEGER :: JIJ, JK
INTEGER :: IKTB, IKTE, IIJB, IIJE
INTEGER :: ISTIJ, ISTK
!
LOGICAL :: GEXT_TEND
LOGICAL :: LSOFT ! Must we really compute tendencies or only adjust them to new T variables
INTEGER :: INB_ITER_MAX ! Maximum number of iterations (with real tendencies computation)
REAL :: ZTSTEP ! length of sub-timestep in case of time splitting
REAL :: ZINV_TSTEP ! Inverse ov PTSTEP
REAL :: ZTIME_THRESHOLD ! Time to reach threshold
!
!Output packed total mixing ratio change (for budgets only)
REAL, DIMENSION(KSIZE, IBUNUM-IBUNUM_EXTRA) :: ZBU_PACK
!
!For packing
INTEGER :: IMICRO ! Case r_x>0 locations
INTEGER :: JL, JV, JJV
REAL, DIMENSION(KPROMA) :: ZTIME ! Current integration time (starts with 0 and ends with PTSTEP)
REAL, DIMENSION(KPROMA) :: &
                        & ZMAXTIME, & ! Time on which we can apply the current tendencies
                        & ZTIME_LASTCALL, &     ! Integration time when last tendecies call has been done
                        & ZSSI,     &
                        & ZCIT,     & ! Pristine ice conc. at t
                        & ZRHODREF, & ! RHO Dry REFerence
                        & ZZT,      & ! Temperature
                        & ZPRES,    & ! Pressure
                        & ZEXN,     & ! EXNer Pressure
                        & ZLSFACT,  & ! L_s/(Pi*C_ph)
                        & ZLVFACT,  & ! L_v/(Pi*C_ph)
                        & ZSIGMA_RC,& ! Standard deviation of rc at time t
                        & ZCF,      & ! Cloud fraction
                        & ZHLC_HCF, & ! HLCLOUDS : fraction of High Cloud Fraction in grid
                        & ZHLC_LCF, & ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                      !    note that ZCF = ZHLC_HCF + ZHLC_LCF
                        & ZHLC_HRC, & ! HLCLOUDS : LWC that is High LWC in grid
                        & ZHLC_LRC, & ! HLCLOUDS : LWC that is Low  LWC in grid
                                      !    note that ZRC = ZHLC_HRC + ZHLC_LRC
                        & ZHLI_HCF, &
                        & ZHLI_LCF, &
                        & ZHLI_HRI, &
                        & ZHLI_LRI
LOGICAL, DIMENSION(KPROMA) :: LLCOMPUTE ! .TRUE. or points where we must compute tendencies,
!
!Output packed tendencies (for budgets only)
REAL, DIMENSION(KPROMA, IBUNUM-IBUNUM_EXTRA) :: ZBU_SUM
REAL, DIMENSION(KPROMA, IBUNUM) :: ZBU_INST
!
!For mixing-ratio-splitting
LOGICAL :: LLCPZ0RT
REAL :: ZTIME_THRESHOLD1D(KPROMA) ! Time to reach threshold
REAL, DIMENSION(KPROMA, KRR) :: Z0RT ! Mixing-ratios at the beginig of the current loop
!
REAL, DIMENSION(KPROMA,0:7) :: &
                        & ZVART, & !Packed variables
                        & ZEXTPK, & !To take into acount external tendencies inside the splitting
                        & ZA, ZB
!
REAL, DIMENSION(KPROMA, 8) :: ZRS_TEND, ZRG_TEND
REAL, DIMENSION(KPROMA,10) :: ZRH_TEND

INTEGER, DIMENSION(KPROMA) :: &
                       & I1,I2, & ! Used to replace the COUNT and PACK intrinsics on variables
                       & IITER    ! Number of iterations done (with real tendencies computation)
INTEGER, DIMENSION(KSIZE) :: I1TOT, I2TOT ! Used to replace the COUNT and PACK intrinsics
!
REAL, DIMENSION(KPROMA) :: ZSUM2, ZMAXB
REAL :: ZDEVIDE, ZX
!
INTEGER :: IC, JMICRO
LOGICAL :: LLSIGMA_RC, LL_ANY_ITER, LL_AUCV_ADJU
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_PACK', 0, ZHOOK_HANDLE)
!
!*       1.     GENERALITIES
!               ------------
!
IKTB=D%NKTB
IKTE=D%NKTE
IIJB=D%NIJB
IIJE=D%NIJE
GEXT_TEND=.TRUE.
ZINV_TSTEP=1./PTSTEP
!
IF(KPROMA /= KSIZE) THEN
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'RAIN_ICE', 'For now, KPROMA must be equal to KSIZE, see code for explanation')
  ! Microphyscs was optimized by introducing chunks of KPROMA size
  ! Thus, in ice4_tendencies, the 1D array represent only a fraction of the points where microphisical species are present
  ! We cannot rebuild the entire 3D arrays in the subroutine, so we cannot call ice4_rainfr_vert in it
  ! A solution would be to suppress optimisation in this case by setting KPROMA=KSIZE in rain_ice
  ! Another solution would be to compute column by column?
  ! Another one would be to cut tendencies in 3 parts: before rainfr_vert, rainfr_vert, after rainfr_vert
ENDIF
!
IF(BUCONF%LBU_ENABLE) THEN
  DO JV=1, IBUNUM-IBUNUM_EXTRA
    ZBU_PACK(:, JV)=0.
    ZBU_SUM(:, JV)=0.
  ENDDO
ENDIF
!-------------------------------------------------------------------------------
!
!***       4.1 Point selection
!
!  optimization by looking for locations where
!  the microphysical fields are larger than a minimal value only !!!
!
IF (KSIZE /= COUNT(ODMICRO(IIJB:IIJE,IKTB:IKTE))) THEN
    CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'RAIN_ICE', 'RAIN_ICE : KSIZE /= COUNT(ODMICRO)')
ENDIF

IF (KSIZE > 0) THEN

  !Maximum number of iterations
  !We only count real iterations (those for which we *compute* tendencies)
  INB_ITER_MAX=PARAMI%NMAXITER
  IF(PARAMI%XTSTEP_TS/=0.)THEN
    INB_ITER_MAX=MAX(1, INT(PTSTEP/PARAMI%XTSTEP_TS)) !At least the number of iterations needed for the time-splitting
    ZTSTEP=PTSTEP/INB_ITER_MAX
    INB_ITER_MAX=MAX(PARAMI%NMAXITER, INB_ITER_MAX) !For the case XMRSTEP/=0. at the same time
  ENDIF
  LLSIGMA_RC=(HSUBG_AUCV_RC=='PDF ' .AND. PARAMI%CSUBG_PR_PDF=='SIGM')
  LL_AUCV_ADJU=(HSUBG_AUCV_RC=='ADJU' .OR. HSUBG_AUCV_RI=='ADJU')

  !-------------------------------------------------------------------------------
  !
  !***       4.2 Cache-blocking loop
  !

  ! starting indexes :
  IC=0
  ISTK=IKTB
  ISTIJ=IIJB

  DO JMICRO=1,KSIZE,KPROMA

    IMICRO=MIN(KPROMA,KSIZE-JMICRO+1)
    !-------------------------------------------------------------------------------
    !
    !***       4.3 Packing
    !

    ! Setup packing parameters
    OUTER_LOOP: DO JK = ISTK, IKTE
      IF (ANY(ODMICRO(:,JK))) THEN
        DO JIJ = ISTIJ, IIJE
          IF (ODMICRO(JIJ,JK)) THEN
            IC=IC+1
            ! Initialization of variables in packed format :
            ZVART(IC, ITH)=PTHT(JIJ, JK)
            ZVART(IC, IRV)=PRVT(JIJ, JK)
            ZVART(IC, IRC)=PRCT(JIJ, JK)
            ZVART(IC, IRR)=PRRT(JIJ, JK)
            ZVART(IC, IRI)=PRIT(JIJ, JK)
            ZVART(IC, IRS)=PRST(JIJ, JK)
            ZVART(IC, IRG)=PRGT(JIJ, JK)
            IF (KRR==7) THEN
              ZVART(IC, IRH)=PRHT(JIJ, JK)
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
            ! Save indices for later usages:
            I1(IC) = JIJ
            I2(IC) = JK
            I1TOT(JMICRO+IC-1)=JIJ
            I2TOT(JMICRO+IC-1)=JK
            IF (IC==IMICRO) THEN
              ! the end of the chunk has been reached, then reset the starting index :
              ISTIJ=JIJ+1
              IF (ISTIJ <= IIJE) THEN
                ISTK=JK
              ELSE
                ! end of line, restart from 1 and increment upper loop
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

    IF (GEXT_TEND) THEN
      DO JV=0, KRR
        DO JL=1, IMICRO
          ZEXTPK(JL, JV)=ZEXTPK(JL, JV)-ZVART(JL, JV)*ZINV_TSTEP
        ENDDO
      ENDDO
    ENDIF
    IF (LLSIGMA_RC) THEN
      DO JL=1, IMICRO
        ZSIGMA_RC(JL)=ZSIGMA_RC(JL)*2.
      ENDDO 
    ENDIF
    IF (LL_AUCV_ADJU) THEN
      DO JL=1, IMICRO
        ZHLC_LRC(JL) = ZVART(JL, IRC) - ZHLC_HRC(JL)
        ZHLI_LRI(JL) = ZVART(JL, IRI) - ZHLI_HRI(JL)
        IF(ZVART(JL, IRC)>0.) THEN
          ZHLC_LCF(JL) = ZCF(JL)- ZHLC_HCF(JL)
        ELSE
          ZHLC_LCF(JL)=0.
        ENDIF
        IF(ZVART(JL, IRI)>0.) THEN
          ZHLI_LCF(JL) = ZCF(JL)- ZHLI_HCF(JL)
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
    IITER(1:IMICRO)=0
    ZTIME(1:IMICRO)=0. ! Current integration time (all points may have a different integration time)

    DO WHILE(ANY(ZTIME(1:IMICRO)<PTSTEP)) ! Loop to *really* compute tendencies

      IF(PARAMI%XTSTEP_TS/=0.) THEN
        ! In this case we need to remember the time when tendencies were computed
        ! because when time has evolved more than a limit, we must re-compute tendencies
        ZTIME_LASTCALL(1:IMICRO)=ZTIME(1:IMICRO)
      ENDIF
      DO JL=1, IMICRO
        IF (ZTIME(JL) < PTSTEP) THEN
          LLCOMPUTE(JL)=.TRUE. ! Computation (.TRUE.) only for points for which integration time has not reached the timestep
          IITER(JL)=IITER(JL)+1
        ELSE
          LLCOMPUTE(JL)=.FALSE.
        ENDIF
      ENDDO
      LL_ANY_ITER=ANY(IITER(1:IMICRO) < INB_ITER_MAX)
      LLCPZ0RT=.TRUE.
      LSOFT=.FALSE. ! We *really* compute the tendencies

      DO WHILE(ANY(LLCOMPUTE(1:IMICRO))) ! Loop to adjust tendencies when we cross the 0°C or when a species disappears
!$OMP SIMD
        DO JL=1, IMICRO
          ZSUM2(JL)=SUM(ZVART(JL,IRI:KRR))
        ENDDO
        DO JL=1, IMICRO
          ZDEVIDE=(CST%XCPD + CST%XCPV*ZVART(JL, IRV) + CST%XCL*(ZVART(JL, IRC)+ZVART(JL, IRR)) + CST%XCI*ZSUM2(JL)) * ZEXN(JL)
          ZZT(JL) = ZVART(JL, ITH) * ZEXN(JL)
          ZLSFACT(JL)=(CST%XLSTT+(CST%XCPV-CST%XCI)*(ZZT(JL)-CST%XTT)) / ZDEVIDE
          ZLVFACT(JL)=(CST%XLVTT+(CST%XCPV-CST%XCL)*(ZZT(JL)-CST%XTT)) / ZDEVIDE
        ENDDO
        !-------------------------------------------------------------------------------
        !
        !***       4.5 Effective tendencies computation
        !
        !
        ! Tendencies are *really* computed when LSOFT==.FALSE. and only adjusted otherwise
        CALL ICE4_TENDENCIES(D, CST, PARAMI, ICEP, ICED, BUCONF, &
                            &KPROMA, IMICRO, &
                            &KRR, LSOFT, LLCOMPUTE, &
                            &HSUBG_AUCV_RC, HSUBG_AUCV_RI, &
                            &ZEXN, ZRHODREF, ZLVFACT, ZLSFACT, I1, I2, &
                            &ZPRES, ZCF, ZSIGMA_RC, &
                            &ZCIT, &
                            &ZZT, ZVART, &
                            &ZBU_INST, &
                            &ZRS_TEND, ZRG_TEND, ZRH_TEND, ZSSI, &
                            &ZA, ZB, &
                            &ZHLC_HCF, ZHLC_LCF, ZHLC_HRC, ZHLC_LRC, &
                            &ZHLI_HCF, ZHLI_LCF, ZHLI_HRI, ZHLI_LRI, PRAINFR)

        ! External tendencies
        IF(GEXT_TEND) THEN
          DO JV=0, KRR
            DO JL=1, IMICRO
              ZA(JL, JV) = ZA(JL, JV) + ZEXTPK(JL, JV)
            ENDDO
          ENDDO
        ENDIF
        !-------------------------------------------------------------------------------
        !
        !***       4.6 Time integration
        !
        !
        ! If we can, we shall use these tendencies until the end of the timestep
        DO JL=1, IMICRO
          IF(LLCOMPUTE(JL)) THEN
            ZMAXTIME(JL)=(PTSTEP-ZTIME(JL)) ! Remaining time until the end of the timestep
          ELSE
            ZMAXTIME(JL)=0.
          ENDIF
        ENDDO

        !We need to adjust tendencies when temperature reaches 0
        IF(PARAMI%LFEEDBACKT) THEN
          DO JL=1, IMICRO
            !Is ZB(:, ITH) enough to change temperature sign?
            ZX=CST%XTT/ZEXN(JL)
            IF ((ZVART(JL, ITH) - ZX) * (ZVART(JL, ITH) + ZB(JL, ITH) - ZX) < 0.) THEN
              ZMAXTIME(JL)=0.
            ENDIF
            !Can ZA(:, ITH) make temperature change of sign?
            IF (ABS(ZA(JL,ITH)) > 1.E-20 ) THEN
              ZTIME_THRESHOLD=(ZX - ZB(JL, ITH) - ZVART(JL, ITH))/ZA(JL, ITH)
              IF (ZTIME_THRESHOLD > 0.) THEN
                ZMAXTIME(JL)=MIN(ZMAXTIME(JL), ZTIME_THRESHOLD)
              ENDIF
            ENDIF
          ENDDO
        ENDIF

        !We need to adjust tendencies when a species disappears
        !When a species is missing, only the external tendencies can be negative (and we must keep track of it)
        DO JV=1, KRR
          DO JL=1, IMICRO
            IF (ZA(JL, JV) < -1.E-20 .AND. ZVART(JL, JV) > ICED%XRTMIN(JV)) THEN
              ZMAXTIME(JL)=MIN(ZMAXTIME(JL), -(ZB(JL, JV)+ZVART(JL, JV))/ZA(JL, JV))
            ENDIF
          ENDDO
        ENDDO

        !We stop when the end of the timestep is reached
        DO JL=1, IMICRO
          IF (ZTIME(JL)+ZMAXTIME(JL) >= PTSTEP) THEN
            LLCOMPUTE(JL)=.FALSE.
          ENDIF
        ENDDO
        !We must recompute tendencies when the end of the sub-timestep is reached
        IF (PARAMI%XTSTEP_TS/=0.) THEN
          DO JL=1, IMICRO
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
            ! Thus, at first iteration (ie when LLCPZ0RT=.TRUE.) we copy ZVART into Z0RT
            DO JV=1,KRR
              IF (LLCPZ0RT) Z0RT(1:IMICRO, JV)=ZVART(1:IMICRO, JV)
              DO JL=1, IMICRO
                IF (IITER(JL)<INB_ITER_MAX .AND. ABS(ZA(JL,JV))>1.E-20) THEN
                  ZTIME_THRESHOLD1D(JL)=(SIGN(1., ZA(JL, JV))*PARAMI%XMRSTEP+ &
                                        &Z0RT(JL, JV)-ZVART(JL, JV)-ZB(JL, JV))/ZA(JL, JV)
                ELSE
                  ZTIME_THRESHOLD1D(JL)=-1.
                ENDIF
              ENDDO
              DO JL=1, IMICRO
                IF (ZTIME_THRESHOLD1D(JL)>=0 .AND. ZTIME_THRESHOLD1D(JL)<ZMAXTIME(JL) .AND. &
                   &(ZVART(JL, JV)>ICED%XRTMIN(JV) .OR. ZA(JL, JV)>0.)) THEN
                  ZMAXTIME(JL)=MIN(ZMAXTIME(JL), ZTIME_THRESHOLD1D(JL))
                  LLCOMPUTE(JL)=.FALSE.
                ENDIF
              ENDDO
            ENDDO
            LLCPZ0RT=.FALSE.
!$OMP SIMD
            DO JL=1,IMICRO
              ZMAXB(JL)=MAXVAL(ABS(ZB(JL,1:KRR)))
            ENDDO
            DO JL=1, IMICRO
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
          DO JL=1, IMICRO
            ZVART(JL, JV)=ZVART(JL, JV)+ZA(JL, JV)*ZMAXTIME(JL)+ZB(JL, JV)
          ENDDO
        ENDDO
        DO JL=1, IMICRO
#ifdef REPRO55
          ZCIT(JL)=ZCIT(JL) * MAX(0., -SIGN(1., -ZVART(JL,IRI)))
#else
          IF (ZVART(JL,IRI)<=0.) ZCIT(JL) = 0.
#endif
          ZTIME(JL)=ZTIME(JL)+ZMAXTIME(JL)
        ENDDO
        !-------------------------------------------------------------------------------
        !
        !***       4.8 Mixing ratio change due to each process
        !
        IF(BUCONF%LBU_ENABLE) THEN
          !Mixing ratio change due to a tendency
          DO JV=1, IBUNUM-IBUNUM_MR-IBUNUM_EXTRA
            DO JL=1, IMICRO
              ZBU_SUM(JL, JV) = ZBU_SUM(JL, JV) + ZBU_INST(JL, JV)*ZMAXTIME(JL)
            ENDDO
          ENDDO

          !Mixing ratio change due to a mixing ratio change
          DO JV=IBUNUM-IBUNUM_MR-IBUNUM_EXTRA+1, IBUNUM-IBUNUM_EXTRA
            DO JL=1, IMICRO
              ZBU_SUM(JL, JV) = ZBU_SUM(JL, JV) + ZBU_INST(JL, JV)
            ENDDO
          ENDDO

          !Extra contribution as a mixing ratio change
          DO JV=IBUNUM-IBUNUM_EXTRA+1, IBUNUM
            JJV=IBUEXTRAIND(JV)
            DO JL=1, IMICRO
              ZBU_SUM(JL, JJV) = ZBU_SUM(JL, JJV) + ZBU_INST(JL, JV)
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

    IF(GEXT_TEND) THEN
      !Z..T variables contain the external tendency, we substract it
      DO JV=0, KRR
        DO JL=1, IMICRO
          ZVART(JL, JV) = ZVART(JL, JV) - ZEXTPK(JL, JV) * PTSTEP
        ENDDO
      ENDDO
    ENDIF

    !-------------------------------------------------------------------------------
    !
    !***       4.10 Unpacking
    !
    !
    DO JL=1, IMICRO
      PCIT  (I1(JL),I2(JL))=ZCIT   (JL)
      IF(PARAMI%LWARM) THEN
        PEVAP3D(I1(JL),I2(JL))=ZBU_INST(JL, IRREVAV)
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
    IF(BUCONF%LBU_ENABLE) THEN
      DO JV=1, IBUNUM-IBUNUM_EXTRA
        DO JL=1, IMICRO
          ZBU_PACK(JMICRO+JL-1, JV) = ZBU_SUM(JL, JV)
        ENDDO
      ENDDO
    ENDIF


  ENDDO ! JMICRO
ENDIF ! KSIZE > 0
!-------------------------------------------------------------------------------
!
!***       4.11 Budgets
!
!
IF(BUCONF%LBU_ENABLE) THEN
  !Budgets for the different processes
  CALL ICE4_BUDGETS(D, PARAMI, BUCONF, KSIZE, KPROMA, PTSTEP, KRR, I1TOT, I2TOT, &
                    PLVFACT, PLSFACT, PRHODJ, PEXNREF, &
                    PRVHENI, ZBU_PACK, &
                    TBUDGETS, KBUDGETS)
ENDIF
IF (LHOOK) CALL DR_HOOK('ICE4_PACK', 1, ZHOOK_HANDLE)
END SUBROUTINE ICE4_PACK
END MODULE MODE_ICE4_PACK
