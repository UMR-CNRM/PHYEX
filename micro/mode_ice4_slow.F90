!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_ICE4_SLOW
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_SLOW(CST, PARAMI, ICEP, ICED, D, LDSOFT, OELEC, LDCOMPUTE, PRHODREF, PT, &
                     &PSSI, PLVFACT, PLSFACT, &
                     &PRVT, PRCT, PRIT, PRST, PRGT, &
                     &PLBDAS, PLBDAG, &
                     &PAI, PCJ, PHLI_HCF, PHLI_HRI,&
                     &PLATHAM_IAGGS, &
                     &PRCHONI, PRVDEPS, PRIAGGS, PRIAUTS, PRVDEPG)

!$ACDC singlecolumn

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
!  K.I.Ivarsson,  02/2023: Introduce possibility to reduce graupel in non-convective conditions and
!        high ice super saturation (set by namelist)
!  C. Barthe      06/2023: Add retroaction of electric field on IAGGS
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),                  INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),            INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),       INTENT(IN)    :: ICED
TYPE(DIMPHYEX_t),             INTENT(IN)    :: D
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL,                      INTENT(IN)    :: OELEC
LOGICAL, DIMENSION(D%NIJT),   INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PSSI     ! Supersaturation over ice
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLSFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRVT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRGT     ! Graupel/hail m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PAI      ! Thermodynamical function
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PHLI_HCF !
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PHLI_HRI !
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC)), INTENT(IN) :: PLATHAM_IAGGS ! enhancement factor of IAGGS due to Efield
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRCHONI  ! Homogeneous nucleation
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRVDEPS  ! Deposition on r_s
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRIAGGS  ! Aggregation on r_s
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRIAUTS  ! Autoconversion of r_i for r_s production
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRVDEPG  ! Deposition on r_g
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT) :: ZCRIAUTI
INTEGER                 :: JIJ
REAL            :: ZREDGR,ZREDSN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_SLOW', 0, ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*       3.1  Tuning possibilties (only for OCND2)

ZREDGR  = 1.      ! Tuning of the deposition of graupel, 1. is ref. value 
ZREDSN  = 1.      ! Tuning of the deposition of snow, 1. is ref. value

IF(PARAMI%LOCND2) THEN
  IF(.NOT. PARAMI%LMODICEDEP) THEN
    ZREDGR  = ICEP%XFRMIN(39)  ! Tuning factor, may be /= 1. 
    ZREDSN  = ICEP%XFRMIN(40)  ! Tuning factor, may be /= 1.
  ENDIF
ENDIF
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
#endif
!
!*       3.2     compute the homogeneous nucleation source: RCHONI
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PT(JIJ)<CST%XTT-35.0 .AND. PRCT(JIJ)>ICED%XRTMIN(2) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRCHONI(JIJ) = MIN(1000.,ICEP%XHON*PRHODREF(JIJ)*PRCT(JIJ)       &
                                 *EXP( ICEP%XALPHA3*(PT(JIJ)-CST%XTT)-ICEP%XBETA3 ))
    ENDIF
  ELSE
    PRCHONI(JIJ) = 0.
  ENDIF
ENDDO
!$mnh_end_do()

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

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRVT(JIJ)>ICED%XRTMIN(1) .AND. PRST(JIJ)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      IF(.NOT. ICEP%LNEWCOEFF) THEN
        PRVDEPS(JIJ) = ( PSSI(JIJ)/(PRHODREF(JIJ)*PAI(JIJ)) ) *                               &
                   ( ICEP%X0DEPS*PLBDAS(JIJ)**ICEP%XEX0DEPS + ICEP%X1DEPS*PCJ(JIJ)*PLBDAS(JIJ)**ICEP%XEX1DEPS )
        PRVDEPS(JIJ) = PRVDEPS(JIJ)*ZREDSN
      ELSE
        PRVDEPS(JIJ) = ( PRST(JIJ)*(PSSI(JIJ)/PAI(JIJ)) ) *                               &
                      ( ICEP%X0DEPS*PLBDAS(JIJ)**(ICED%XBS+ICEP%XEX0DEPS) + ICEP%X1DEPS*PCJ(JIJ) * &
                      (1+0.5*(ICED%XFVELOS/PLBDAS(JIJ))**ICED%XALPHAS)**(-ICED%XNUS+ICEP%XEX1DEPS/ICED%XALPHAS) &
                       *(PLBDAS(JIJ))**(ICED%XBS+ICEP%XEX1DEPS) )
      ENDIF
    ENDIF
  ELSE
    PRVDEPS(JIJ) = 0.
  ENDIF
ENDDO
!$mnh_end_do()

!
!*       3.4.4  compute the aggregation on r_s: RIAGGS
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRIT(JIJ)>ICED%XRTMIN(4) .AND. PRST(JIJ)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      IF(.NOT. ICEP%LNEWCOEFF) THEN
        PRIAGGS(JIJ) = ICEP%XFIAGGS * EXP( ICEP%XCOLEXIS*(PT(JIJ)-CST%XTT) ) &
                           * PRIT(JIJ)                      &
                           * PLBDAS(JIJ)**ICEP%XEXIAGGS          &
                           * PRHODREF(JIJ)**(-ICED%XCEXVT)
      ELSE
        PRIAGGS(JIJ) = ICEP%XFIAGGS * EXP( ICEP%XCOLEXIS*(PT(JIJ)-CST%XTT) ) &
                           * PRIT(JIJ)                      &
                           * PRST(JIJ) * (1+(ICED%XFVELOS/PLBDAS(JIJ))**ICED%XALPHAS)**&
                           (-ICED%XNUS+ICEP%XEXIAGGS/ICED%XALPHAS) &
                           * PRHODREF(JIJ)**(-ICED%XCEXVT+1.) &
                           * ((PLBDAS(JIJ))**(ICED%XBS+ICEP%XEXIAGGS))
      ENDIF
      IF (OELEC) PRIAGGS(JIJ) = PRIAGGS(JIJ) * PLATHAM_IAGGS(JIJ)
    ENDIF
  ELSE
    PRIAGGS(JIJ) = 0.
  ENDIF
ENDDO
!$mnh_end_do()

!
!*       3.4.5  compute the autoconversion of r_i for r_s production: RIAUTS
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PHLI_HRI(JIJ)>ICED%XRTMIN(4) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      !ZCRIAUTI(:)=MIN(ICEP%XCRIAUTI,10**(0.06*(PT(:)-CST%XTT)-3.5))
      ZCRIAUTI(JIJ)=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(PT(JIJ)-CST%XTT)+ICEP%XBCRIAUTI))
      PRIAUTS(JIJ) = ICEP%XTIMAUTI * EXP( ICEP%XTEXAUTI*(PT(JIJ)-CST%XTT) ) &
                                  * MAX(PHLI_HRI(JIJ)-ZCRIAUTI(JIJ)*PHLI_HCF(JIJ), 0.)
    ENDIF
  ELSE
    PRIAUTS(JIJ) = 0.
  ENDIF
ENDDO
!$mnh_end_do()

!
!*       3.4.6  compute the deposition on r_g: RVDEPG
!
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRVT(JIJ)>ICED%XRTMIN(1) .AND. PRGT(JIJ)>ICED%XRTMIN(6) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRVDEPG(JIJ) = ( PSSI(JIJ)/(PRHODREF(JIJ)*PAI(JIJ)) ) *                               &
                 ( ICEP%X0DEPG*PLBDAG(JIJ)**ICEP%XEX0DEPG + ICEP%X1DEPG*PCJ(JIJ)*PLBDAG(JIJ)**ICEP%XEX1DEPG )
      PRVDEPG(JIJ) = PRVDEPG(JIJ)*ZREDGR
    ENDIF
  ELSE
    PRVDEPG(JIJ) = 0.
  ENDIF
ENDDO
!$mnh_end_do()

!
IF (LHOOK) CALL DR_HOOK('ICE4_SLOW', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_SLOW
END MODULE MODE_ICE4_SLOW
