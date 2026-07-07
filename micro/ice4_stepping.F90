!MNH_LIC Copyright 1995-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
SUBROUTINE ICE4_STEPPING(D, CST, PARAMI, ICEP, ICED, BUCONF, PTSTEP, &
                        &KRR, OSAVE_MICRO, LDMICRO, OELEC, &
                        &PEXN, PRHODREF, PPABST, PCIT, PCLDFR, &
                        &PHLC_HCF, PHLC_HRC, PHLI_HCF, PHLI_HRI, &
                        &PTHS, PRS, PRREVAV, PRAINFR, PSIGS, PTHT, PRT, &
                        &PICLDFR, PZZZ, PCONC3D, PSSIO, PSSIU, PIFR, &
                        &PBUDGETS, PLATHAM_IAGGS)

!$ACDC singlecolumn --nocreate-interface

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
      & IBUNUM,       & ! Number of tendency terms
      & IBUNUM_MR,    & ! Number of tendency terms expressed as mixing ratio changes
      & IBUNUM_EXTRA, & ! Number of extra tendency terms
      & IRREVAV,      & ! Index for the evaporation tendency
      & IBUEXTRAIND
! Index indirection
!
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
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL,                  INTENT(IN)    :: OSAVE_MICRO   ! if true, save the microphysical tendencies
LOGICAL, DIMENSION(D%NIJT), INTENT(IN)  :: LDMICRO
LOGICAL,                  INTENT(IN)    :: OELEC         ! if true, cloud electricity is activated
!
REAL,    DIMENSION(D%NIJT),                     INTENT(IN)    :: PEXN    ! Exner function
REAL,    DIMENSION(D%NIJT),                     INTENT(IN)    :: PRHODREF! Reference density
REAL,    DIMENSION(D%NIJT),                     INTENT(IN)    :: PPABST

REAL,    DIMENSION(D%NIJT),                     INTENT(INOUT) :: PCIT

REAL,    DIMENSION(D%NIJT),                     INTENT(IN)    :: PCLDFR ! Cloud fraction
REAL,    DIMENSION(D%NIJT),                     INTENT(IN)    :: PICLDFR ! Ice cloud fraction
REAL,    DIMENSION(D%NIJT),                     INTENT(IN)    :: PZZZ    ! Model level height
REAL,    DIMENSION(D%NIJT),                     INTENT(IN)    :: PCONC3D ! Cloud croplet number concentration
REAL,    DIMENSION(D%NIJT),                     INTENT(IN)    :: PSSIO   ! Super-saturation with respect to ice in the supersaturated fraction
REAL,    DIMENSION(D%NIJT),                     INTENT(IN)    :: PSSIU   ! Sub-saturation with respect to ice in the  subsaturated fraction
REAL,    DIMENSION(D%NIJT),                     INTENT(IN)    :: PIFR    ! Ratio cloud ice moist part to dry part
REAL,    DIMENSION(MERGE(D%NIJT,0,PARAMI%CSUBG_AUCV_RC=='ADJU' .OR. PARAMI%CSUBG_AUCV_RI=='ADJU')), INTENT(INOUT) :: PHLC_HRC, &
                                                                                                                  & PHLC_HCF, &
                                                                                                                  & PHLI_HRI, &
                                                                                                                  & PHLI_HCF
REAL,    DIMENSION(MERGE(D%NIJT,0,PARAMI%LEXT_TEND)),   INTENT(IN) :: PTHS !To take into acount external tendencies inside the splitting
REAL,    DIMENSION(MERGE(D%NIJT,0,PARAMI%LEXT_TEND),7), INTENT(IN) :: PRS !To take into acount external tendencies inside the splitting
REAL,    DIMENSION(D%NIJT),                     INTENT(OUT)   :: PRREVAV
REAL,    DIMENSION(D%NIJT),                     INTENT(INOUT) :: PRAINFR
REAL,    DIMENSION(MERGE(D%NIJT,0,PARAMI%CSUBG_AUCV_RC=='PDF ' .AND. PARAMI%CSUBG_PR_PDF=='SIGM')), INTENT(IN) :: PSIGS
REAL,    DIMENSION(D%NIJT),                     INTENT(INOUT) :: PTHT
REAL,    DIMENSION(D%NIJT,7),                   INTENT(INOUT) :: PRT !Packed variables
REAL,    DIMENSION(MERGE(D%NIJT,0,BUCONF%LBU_ENABLE .OR. OSAVE_MICRO), &
                   MERGE(IBUNUM-IBUNUM_EXTRA,0,BUCONF%LBU_ENABLE .OR. OSAVE_MICRO)),INTENT(OUT)   :: PBUDGETS
REAL,    DIMENSION(MERGE(D%NIJT,0,OELEC)),      INTENT(IN)    :: PLATHAM_IAGGS ! E Function to simulate
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
INTEGER :: JIJ, JV, JJV
REAL, DIMENSION(D%NIJT) :: &
                        & ZTIME, & ! Current integration time (starts with 0 and ends with PTSTEP)
                        & ZMAXTIME, & ! Time on which we can apply the current tendencies
                        & ZTIME_LASTCALL, &     ! Integration time when last tendecies call has been done
                        & ZSSI,     &
                        & ZZT,      & ! Temperature
                        & ZLSFACT,  & ! L_s/(Pi*C_ph)
                        & ZLVFACT,  & ! L_v/(Pi*C_ph)
                        & ZHLC_LCF, & ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                      !    note that PCLDFR = PHLC_HCF + ZHLC_LCF
                        & ZHLC_LRC, & ! HLCLOUDS : LWC that is Low  LWC in grid
                                      !    note that ZRC = PHLC_HRC + ZHLC_LRC
                        & ZHLI_LCF, &
                        & ZHLI_LRI, &
                        & ZSIGMA_RC, &
                        & ZEXTTH
REAL, DIMENSION(D%NIJT,7) :: ZEXTPK
LOGICAL, DIMENSION(D%NIJT) :: LLCOMPUTE ! .TRUE. or points where we must compute tendencies,
REAL, DIMENSION(SIZE(ICED%XRTMIN))   :: ZRSMIN
!
!Output packed tendencies (for budgets only)
REAL, DIMENSION(D%NIJT, IBUNUM) :: ZBU_INST
!
!For mixing-ratio-splitting
LOGICAL :: LLCPZ0RT
REAL :: ZTIME_THRESHOLD1D(D%NIJT) ! Time to reach threshold
REAL, DIMENSION(D%NIJT, KRR) :: Z0RT ! Mixing-ratios at the beginig of the current loop
!
REAL, DIMENSION(D%NIJT,KRR) :: ZA, ZB
REAL, DIMENSION(D%NIJT)   :: ZATH, ZBTH
!
REAL, DIMENSION(D%NIJT, 8) :: ZRS_TEND, ZRG_TEND
REAL, DIMENSION(D%NIJT,10) :: ZRH_TEND

INTEGER, DIMENSION(D%NIJT) :: IITER    ! Number of iterations done (with real tendencies computation)
!
REAL, DIMENSION(D%NIJT) :: ZSUM2, ZMAXB
REAL :: ZDEVIDE, ZX
!
LOGICAL :: LL_ANY_ITER
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_STEPPING', 0, ZHOOK_HANDLE)
!
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
#endif
!
!*       1.     GENERALITIES
!               ------------
!
ZINV_TSTEP=1./PTSTEP
ZRSMIN = ICED%XRTMIN
!
IF(BUCONF%LBU_ENABLE .OR. OSAVE_MICRO) THEN
  DO JV=1, IBUNUM-IBUNUM_EXTRA
    PBUDGETS(:, JV)=0.
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

IF (PARAMI%LEXT_TEND) THEN

  DO JIJ=D%NIJB, D%NIJE 
    ZEXTTH(JIJ)=PTHS(JIJ)-PTHT(JIJ)*ZINV_TSTEP
  ENDDO



  DO JV=1, KRR
    DO JIJ=D%NIJB, D%NIJE 
      ZEXTPK(JIJ, JV)=PRS(JIJ, JV)-PRT(JIJ, JV)*ZINV_TSTEP
    ENDDO
  ENDDO

ENDIF
IF (PARAMI%CSUBG_AUCV_RC=='PDF ' .AND. PARAMI%CSUBG_PR_PDF=='SIGM') THEN


  DO JIJ=D%NIJB, D%NIJE 
    ZSIGMA_RC(JIJ)=PSIGS(JIJ)*2.
  ENDDO

ENDIF
IF (PARAMI%CSUBG_AUCV_RC=='ADJU' .OR. PARAMI%CSUBG_AUCV_RI=='ADJU') THEN


  DO JIJ=D%NIJB, D%NIJE 
    ZHLC_LRC(JIJ) = PRT(JIJ, IRC) - PHLC_HRC(JIJ)
    ZHLI_LRI(JIJ) = PRT(JIJ, IRI) - PHLI_HRI(JIJ)
    IF(PRT(JIJ, IRC)>0.) THEN
      ZHLC_LCF(JIJ) = PCLDFR(JIJ)- PHLC_HCF(JIJ)
    ELSE
      ZHLC_LCF(JIJ)=0.
    ENDIF
    IF(PRT(JIJ, IRI)>0.) THEN
      ZHLI_LCF(JIJ) = PCLDFR(JIJ)- PHLI_HCF(JIJ)
    ELSE
      ZHLI_LCF(JIJ)=0.
    ENDIF
  ENDDO

ENDIF

!-------------------------------------------------------------------------------
!
!***       4.4 temporal loop
!
!

IITER(D%NIJB:D%NIJE )=0

DO JIJ=D%NIJB, D%NIJE 
  IF(LDMICRO(JIJ)) THEN
    ZTIME(JIJ)=0. ! Current integration time (all points may have a different integration time)
  ELSE
    ZTIME(JIJ)=PTSTEP ! Nothing to do on this point, it has already reached the end of the timestep
  ENDIF
ENDDO


DO WHILE(ANY(ZTIME(D%NIJB:D%NIJE )<PTSTEP)) ! Loop to *really* compute tendencies

  IF(PARAMI%XTSTEP_TS/=0.) THEN
    ! In this case we need to remember the time when tendencies were computed
    ! because when time has evolved more than a limit, we must re-compute tendencies
    ZTIME_LASTCALL(D%NIJB:D%NIJE )=ZTIME(D%NIJB:D%NIJE )
  ENDIF


  DO JIJ=D%NIJB, D%NIJE 
    IF (ZTIME(JIJ) < PTSTEP) THEN
      LLCOMPUTE(JIJ)=.TRUE. ! Computation (.TRUE.) only for points for which integration time has not reached the timestep
      IITER(JIJ)=IITER(JIJ)+1
    ELSE
      LLCOMPUTE(JIJ)=.FALSE.
    ENDIF
  ENDDO

  LL_ANY_ITER=ANY(IITER(D%NIJB:D%NIJE ) < INB_ITER_MAX)
  LLCPZ0RT=.TRUE.
  LSOFT=.FALSE. ! We *really* compute the tendencies

  DO WHILE(ANY(LLCOMPUTE(D%NIJB:D%NIJE ))) ! Loop to adjust tendencies when we cross the 0°C or when a species disappears

    ZSUM2(D%NIJB:D%NIJE )=PRT(D%NIJB:D%NIJE , IRI)

    DO JV=IRI+1,KRR
      DO JIJ=D%NIJB, D%NIJE 
        ZSUM2(JIJ)=ZSUM2(JIJ)+PRT(JIJ, JV)
      ENDDO
    ENDDO



    DO JIJ=D%NIJB, D%NIJE 
      ZDEVIDE=(CST%XCPD + CST%XCPV*PRT(JIJ, IRV) + CST%XCL*(PRT(JIJ, IRC)+PRT(JIJ, IRR)) + CST%XCI*ZSUM2(JIJ)) * PEXN(JIJ)
      ZZT(JIJ) = PTHT(JIJ) * PEXN(JIJ)
      ZLSFACT(JIJ)=(CST%XLSTT+(CST%XCPV-CST%XCI)*(ZZT(JIJ)-CST%XTT)) / ZDEVIDE
      ZLVFACT(JIJ)=(CST%XLVTT+(CST%XCPV-CST%XCL)*(ZZT(JIJ)-CST%XTT)) / ZDEVIDE
    ENDDO

    !-------------------------------------------------------------------------------
    !
    !***       4.5 Effective tendencies computation
    !
    !
    ! Tendencies are *really* computed when LSOFT==.FALSE. and only adjusted otherwise
    CALL ICE4_TENDENCIES(CST, PARAMI, ICEP, ICED, BUCONF, &
                        &D, &
                        &KRR, LSOFT, LLCOMPUTE, &
                        &OSAVE_MICRO, OELEC, &
                        &PEXN, PRHODREF, ZLVFACT, ZLSFACT, &
                        &PPABST, PCLDFR, ZSIGMA_RC, &
                        &PCIT, &
                        &ZZT, PICLDFR, PZZZ, PCONC3D, &
                        &PSSIO, PSSIU, PIFR, PTHT, PRT, &
                        &PLATHAM_IAGGS, &
                        &ZBU_INST, &
                        &ZRS_TEND, ZRG_TEND, ZRH_TEND, ZSSI, &
                        &ZA, ZB, ZATH, ZBTH, &
                        &PHLC_HCF, ZHLC_LCF, PHLC_HRC, ZHLC_LRC, &
                        &PHLI_HCF, ZHLI_LCF, PHLI_HRI, ZHLI_LRI, PRAINFR)

    ! External tendencies
    IF(PARAMI%LEXT_TEND) THEN

      DO JIJ=D%NIJB, D%NIJE 
        ZATH(JIJ) = ZATH(JIJ) + ZEXTTH(JIJ)
      ENDDO



      DO JV=1, KRR
        DO JIJ=D%NIJB, D%NIJE 
          ZA(JIJ, JV) = ZA(JIJ, JV) + ZEXTPK(JIJ, JV)
        ENDDO
      ENDDO

    ENDIF
    !-------------------------------------------------------------------------------
    !
    !***       4.6 Time integration
    !
    !
    ! If we can, we shall use these tendencies until the end of the timestep


    DO JIJ=D%NIJB, D%NIJE 
      IF(LLCOMPUTE(JIJ)) THEN
        ZMAXTIME(JIJ)=(PTSTEP-ZTIME(JIJ)) ! Remaining time until the end of the timestep
      ELSE
        ZMAXTIME(JIJ)=0.
      ENDIF
    ENDDO

    !We need to adjust tendencies when temperature reaches 0
    IF(PARAMI%LFEEDBACKT) THEN

      !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE  )
      DO JIJ=D%NIJB, D%NIJE 
        !Is ZB(:, ITH) enough to change temperature sign?
        ZX=CST%XTT/PEXN(JIJ)
        IF ((PTHT(JIJ) - ZX) * (PTHT(JIJ) + ZBTH(JIJ) - ZX) < 0.) THEN
          ZMAXTIME(JIJ)=0.
        ENDIF
        !Can ZATH(:) make temperature change of sign?
        IF (ABS(ZATH(JIJ)) > 1.E-20 ) THEN
          ZTIME_THRESHOLD=(ZX - ZBTH(JIJ) - PTHT(JIJ))/ZATH(JIJ)
          IF (ZTIME_THRESHOLD > 0.) THEN
            ZMAXTIME(JIJ)=MIN(ZMAXTIME(JIJ), ZTIME_THRESHOLD)
          ENDIF
        ENDIF
      ENDDO
      !$mnh_end_do()  

    ENDIF

    !We need to adjust tendencies when a species disappears
    !When a species is missing, only the external tendencies can be negative (and we must keep track of it)

    DO JV=1, KRR

      DO JIJ=D%NIJB, D%NIJE 
        IF (ZA(JIJ, JV) < -1.E-20 .AND. PRT(JIJ, JV) > ZRSMIN(JV)) THEN
          ZMAXTIME(JIJ)=MIN(ZMAXTIME(JIJ), -(ZB(JIJ, JV)+PRT(JIJ, JV))/ZA(JIJ, JV))
          ZMAXTIME(JIJ)=MAX(ZMAXTIME(JIJ), CST%XMNH_TINY) !to prevent rounding errors
        ENDIF
      ENDDO
    ENDDO

    !We stop when the end of the timestep is reached


    DO JIJ=D%NIJB, D%NIJE 
      IF (ZTIME(JIJ)+ZMAXTIME(JIJ) >= PTSTEP) THEN
        LLCOMPUTE(JIJ)=.FALSE.
      ENDIF
    ENDDO

    !We must recompute tendencies when the end of the sub-timestep is reached
    IF (PARAMI%XTSTEP_TS/=0.) THEN


      DO JIJ=D%NIJB, D%NIJE 
        IF ((IITER(JIJ) < INB_ITER_MAX) .AND. (ZTIME(JIJ)+ZMAXTIME(JIJ) > ZTIME_LASTCALL(JIJ)+ZTSTEP)) THEN
          ZMAXTIME(JIJ)=ZTIME_LASTCALL(JIJ)-ZTIME(JIJ)+ZTSTEP
          LLCOMPUTE(JIJ)=.FALSE.
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
        ! Thus, at first iteration (ie when LLCPZ0RT=.TRUE.) we copy PRT into Z0RT
        DO JV=1,KRR

          IF (LLCPZ0RT) THEN
            Z0RT(D%NIJB:D%NIJE , JV)=PRT(D%NIJB:D%NIJE , JV)
          ENDIF

          DO JIJ=D%NIJB, D%NIJE 
            IF (IITER(JIJ)<INB_ITER_MAX .AND. ABS(ZA(JIJ,JV))>1.E-20) THEN
              ZTIME_THRESHOLD1D(JIJ)=(SIGN(1., ZA(JIJ, JV))*PARAMI%XMRSTEP+ &
                                    &Z0RT(JIJ, JV)-PRT(JIJ, JV)-ZB(JIJ, JV))/ZA(JIJ, JV)
            ELSE
              ZTIME_THRESHOLD1D(JIJ)=-1.
            ENDIF
          ENDDO

          DO JIJ=D%NIJB, D%NIJE 
            IF (ZTIME_THRESHOLD1D(JIJ)>=0 .AND. ZTIME_THRESHOLD1D(JIJ)<ZMAXTIME(JIJ) .AND. &
               &(PRT(JIJ, JV)>ZRSMIN(JV) .OR. ZA(JIJ, JV)>0.)) THEN
              ZMAXTIME(JIJ)=MIN(ZMAXTIME(JIJ), ZTIME_THRESHOLD1D(JIJ))
              LLCOMPUTE(JIJ)=.FALSE.
            ENDIF
          ENDDO
          IF (JV == 1) THEN

            DO JIJ=D%NIJB, D%NIJE 
              ZMAXB(JIJ)=ABS(ZB(JIJ, JV))
            ENDDO
          ELSE

            DO JIJ=D%NIJB, D%NIJE 
              ZMAXB(JIJ)=MAX(ZMAXB(JIJ), ABS(ZB(JIJ, JV)))
            ENDDO
          ENDIF

        ENDDO
        LLCPZ0RT=.FALSE.


        DO JIJ=D%NIJB, D%NIJE 
          IF (IITER(JIJ)<INB_ITER_MAX .AND. ZMAXB(JIJ)>PARAMI%XMRSTEP) THEN
            ZMAXTIME(JIJ)=0.
            LLCOMPUTE(JIJ)=.FALSE.
          ENDIF
        ENDDO

      ENDIF ! LL_ANY_ITER
    ENDIF ! XMRSTEP/=0.
    !-------------------------------------------------------------------------------
    !
    !***       4.7 New values of variables for next iteration
    !
    !

    DO JIJ=D%NIJB, D%NIJE 
      IF(LDMICRO(JIJ)) THEN
        PTHT(JIJ)=PTHT(JIJ)+ZATH(JIJ)*ZMAXTIME(JIJ)+ZBTH(JIJ)
      ENDIF
    ENDDO



    DO JV=1, KRR
      DO JIJ=D%NIJB, D%NIJE 
        IF(LDMICRO(JIJ)) THEN
          PRT(JIJ, JV)=PRT(JIJ, JV)+ZA(JIJ, JV)*ZMAXTIME(JIJ)+ZB(JIJ, JV)
        ENDIF
      ENDDO
    ENDDO



    DO JIJ=D%NIJB, D%NIJE 
      IF (PRT(JIJ,IRI)<=0. .AND. LDMICRO(JIJ)) PCIT(JIJ) = 0.
      ZTIME(JIJ)=ZTIME(JIJ)+ZMAXTIME(JIJ)
    ENDDO

    !-------------------------------------------------------------------------------
    !
    !***       4.8 Mixing ratio change due to each process
    !
    IF(BUCONF%LBU_ENABLE .OR. OSAVE_MICRO) THEN
      !Mixing ratio change due to a tendency


      DO JV=1, IBUNUM-IBUNUM_MR-IBUNUM_EXTRA
        DO JIJ=D%NIJB, D%NIJE 
          PBUDGETS(JIJ, JV) = PBUDGETS(JIJ, JV) + ZBU_INST(JIJ, JV)*ZMAXTIME(JIJ)
        ENDDO
      ENDDO


      !Mixing ratio change due to a mixing ratio change


      DO JV=IBUNUM-IBUNUM_MR-IBUNUM_EXTRA+1, IBUNUM-IBUNUM_EXTRA
        DO JIJ=D%NIJB, D%NIJE 
          PBUDGETS(JIJ, JV) = PBUDGETS(JIJ, JV) + ZBU_INST(JIJ, JV)
        ENDDO
      ENDDO


      !Extra contribution as a mixing ratio change
      DO JV=IBUNUM-IBUNUM_EXTRA+1, IBUNUM
        JJV=IBUEXTRAIND(JV)


        DO JIJ=D%NIJB, D%NIJE 
          PBUDGETS(JIJ, JJV) = PBUDGETS(JIJ, JJV) + ZBU_INST(JIJ, JV)
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

IF(PARAMI%LEXT_TEND) THEN
  !Z..T variables contain the external tendency, we substract it

  DO JIJ=D%NIJB, D%NIJE 
    IF(LDMICRO(JIJ)) THEN
      PTHT(JIJ)=PTHT(JIJ) - ZEXTTH(JIJ) * PTSTEP
    ENDIF
  ENDDO



  DO JV=1, KRR
    DO JIJ=D%NIJB, D%NIJE 
      IF(LDMICRO(JIJ)) THEN
        PRT(JIJ, JV) = PRT(JIJ, JV) - ZEXTPK(JIJ, JV) * PTSTEP
      ENDIF
    ENDDO
  ENDDO

ENDIF

DO JIJ=D%NIJB, D%NIJE 
  PRREVAV(JIJ)=ZBU_INST(JIJ, IRREVAV)
ENDDO

!
IF (LHOOK) CALL DR_HOOK('ICE4_STEPPING', 1, ZHOOK_HANDLE)
END SUBROUTINE ICE4_STEPPING
