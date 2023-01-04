!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_ICE4_SLOW
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_SLOW(CST, ICEP, ICED, KPROMA, KSIZE, LDSOFT, LDCOMPUTE, PRHODREF, PT, &
                     &PSSI, PLVFACT, PLSFACT, &
                     &PRVT, PRCT, PRIT, PRST, PRGT, &
                     &PLBDAS, PLBDAG, &
                     &PAI, PCJ, PHLI_HCF, PHLI_HRI,&
                     &PRCHONI, PRVDEPS, PRIAGGS, PRIAUTS, PRVDEPG)
!!
!!**  PURPOSE
!!    -------
!!      Computes the slow process
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the splitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     R. El Khatib 24-Aug-2021 Optimizations
!  J. Wurtz       03/2022: New snow characteristics with LSNOW_T
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: CST_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),                  INTENT(IN)    :: CST
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),       INTENT(IN)    :: ICED
INTEGER,                      INTENT(IN)    :: KPROMA, KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KPROMA),   INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PSSI     ! Supersaturation over ice
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRVT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRGT     ! Graupel/hail m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PAI      ! Thermodynamical function
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PHLI_HCF !
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PHLI_HRI !
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRCHONI  ! Homogeneous nucleation
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRVDEPS  ! Deposition on r_s
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRIAGGS  ! Aggregation on r_s
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRIAUTS  ! Autoconversion of r_i for r_s production
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRVDEPG  ! Deposition on r_g
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KPROMA) :: ZCRIAUTI
REAL                    :: ZTIMAUTIC
INTEGER                 :: JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_SLOW', 0, ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
!*       3.2     compute the homogeneous nucleation source: RCHONI
!
DO JL=1, KSIZE
  IF(PT(JL)<CST%XTT-35.0 .AND. PRCT(JL)>ICED%XRTMIN(2) .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRCHONI(JL) = MIN(1000.,ICEP%XHON*PRHODREF(JL)*PRCT(JL)       &
                                 *EXP( ICEP%XALPHA3*(PT(JL)-CST%XTT)-ICEP%XBETA3 ))
    ENDIF
  ELSE
    PRCHONI(JL) = 0.
  ENDIF
ENDDO
!
!*       3.4    compute the deposition, aggregation and autoconversion sources
!
!
!*       3.4.2  compute the riming-conversion of r_c for r_i production: RCAUTI
!
!  ZZW(:) = 0.0
!  ZTIMAUTIC = SQRT( ICEP%XTIMAUTI*ICEP%XTIMAUTC )
!  WHERE ( (PRCT(:)>0.0) .AND. (PRIT(:)>0.0) .AND. (PRCS(:)>0.0) )
!    ZZW(:) = MIN( PRCS(:),ZTIMAUTIC * MAX( SQRT( PRIT(:)*PRCT(:) ),0.0 ) )
!    PRIS(:) = PRIS(:) + ZZW(:)
!    PRCS(:) = PRCS(:) - ZZW(:)
!    PTHS(:) = PTHS(:) + ZZW(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(RCAUTI))
!  END WHERE
!
!*       3.4.3  compute the deposition on r_s: RVDEPS
!
DO JL=1, KSIZE
  IF(PRVT(JL)>ICED%XRTMIN(1) .AND. PRST(JL)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
#if defined(REPRO48) || defined(REPRO55)
      PRVDEPS(JL) = ( PSSI(JL)/(PRHODREF(JL)*PAI(JL)) ) *                               &
                 ( ICEP%X0DEPS*PLBDAS(JL)**ICEP%XEX0DEPS + ICEP%X1DEPS*PCJ(JL)*PLBDAS(JL)**ICEP%XEX1DEPS )
#else
  PRVDEPS(JL) = ( PRST(JL)*(PSSI(JL)/PAI(JL)) ) *                               &
                 ( ICEP%X0DEPS*PLBDAS(JL)**(ICED%XBS+ICEP%XEX0DEPS) + ICEP%X1DEPS*PCJ(JL) * &
                 (1+0.5*(ICED%XFVELOS/PLBDAS(JL))**ICED%XALPHAS)**(-ICED%XNUS+ICEP%XEX1DEPS/ICED%XALPHAS) &
                  *(PLBDAS(JL))**(ICED%XBS+ICEP%XEX1DEPS) )
#endif
    ENDIF
  ELSE
    PRVDEPS(JL) = 0.
  ENDIF
ENDDO
!
!*       3.4.4  compute the aggregation on r_s: RIAGGS
!
DO JL=1, KSIZE
  IF(PRIT(JL)>ICED%XRTMIN(4) .AND. PRST(JL)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRIAGGS(JL) = ICEP%XFIAGGS * EXP( ICEP%XCOLEXIS*(PT(JL)-CST%XTT) ) &
                         * PRIT(JL)                      &
#if defined(REPRO48) || defined(REPRO55)
                         * PLBDAS(JL)**ICEP%XEXIAGGS          &
                         * PRHODREF(JL)**(-ICED%XCEXVT)
#else
                         * PRST(JL) * (1+(ICED%XFVELOS/PLBDAS(JL))**ICED%XALPHAS)**&
                         (-ICED%XNUS+ICEP%XEXIAGGS/ICED%XALPHAS) &
                         * PRHODREF(JL)**(-ICED%XCEXVT+1.) &
                         * ((PLBDAS(JL))**(ICED%XBS+ICEP%XEXIAGGS))
#endif
    ENDIF
  ELSE
    PRIAGGS(JL) = 0.
  ENDIF
ENDDO
!
!*       3.4.5  compute the autoconversion of r_i for r_s production: RIAUTS
!
DO JL=1, KSIZE
#ifdef REPRO48
  !This was wrong because, with this formulation and in the LDSOFT case, PRIAUTS
  !was not set to 0 when ri is inferior to the autoconversion threshold
  IF(PRIT(JL)>ICED%XRTMIN(4) .AND. LDCOMPUTE(JL)) THEN
#else
  IF(PHLI_HRI(JL)>ICED%XRTMIN(4) .AND. LDCOMPUTE(JL)) THEN
#endif
    IF(.NOT. LDSOFT) THEN
      !ZCRIAUTI(:)=MIN(ICEP%XCRIAUTI,10**(0.06*(PT(:)-CST%XTT)-3.5))
      ZCRIAUTI(JL)=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(PT(JL)-CST%XTT)+ICEP%XBCRIAUTI))
      PRIAUTS(JL) = ICEP%XTIMAUTI * EXP( ICEP%XTEXAUTI*(PT(JL)-CST%XTT) ) &
                                  * MAX(PHLI_HRI(JL)-ZCRIAUTI(JL)*PHLI_HCF(JL), 0.)
    ENDIF
  ELSE
    PRIAUTS(JL) = 0.
  ENDIF
ENDDO
!
!*       3.4.6  compute the deposition on r_g: RVDEPG
!
!
DO JL=1, KSIZE
  IF(PRVT(JL)>ICED%XRTMIN(1) .AND. PRGT(JL)>ICED%XRTMIN(6) .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRVDEPG(JL) = ( PSSI(JL)/(PRHODREF(JL)*PAI(JL)) ) *                               &
                 ( ICEP%X0DEPG*PLBDAG(JL)**ICEP%XEX0DEPG + ICEP%X1DEPG*PCJ(JL)*PLBDAG(JL)**ICEP%XEX1DEPG )
    ENDIF
  ELSE
    PRVDEPG(JL) = 0.
  ENDIF
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ICE4_SLOW', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_SLOW
END MODULE MODE_ICE4_SLOW
