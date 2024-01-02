!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
ELEMENTAL SUBROUTINE ICE4_NUCLEATION(CST, PARAMI, ICEP, ICED, ODCOMPUTE, &
                           PTHT, PPABST, PRHODREF, PEXN, PLSFACT, PT, &
                           PRVT, &
                           PCIT, PRVHENI_MR)
!!
!!**  PURPOSE
!!    -------
!!      Computes the nucleation
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the splitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     R. El Khatib 24-Aug-2021 Optimizations
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),        INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
LOGICAL, INTENT(IN)    :: ODCOMPUTE
REAL,    INTENT(IN)    :: PTHT    ! Theta at t
REAL,    INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL,    INTENT(IN)    :: PRHODREF! Reference density
REAL,    INTENT(IN)    :: PEXN    ! Exner function
REAL,    INTENT(IN)    :: PLSFACT
REAL,    INTENT(IN)    :: PT      ! Temperature at time t
REAL,    INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL,    INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL,    INTENT(OUT)   :: PRVHENI_MR ! Mixing ratio change due to the heterogeneous nucleation
!
!*       0.2  declaration of local variables
!
REAL  :: ZW ! work array
LOGICAL  :: GNEGT  ! Test where to compute the HEN process
REAL   :: ZZW,      & ! Work array
                           ZUSW,     & ! Undersaturation over water
                           ZSSI        ! Supersaturation over ice
!-------------------------------------------------------------------------------
!
!
  IF (ODCOMPUTE) THEN
GNEGT=PT<CST%XTT .AND. PRVT>ICED%XRTMIN(1)
ELSE
GNEGT=.FALSE.
END IF

ZUSW=0.
ZZW=0.

IF (GNEGT) THEN
 ZZW=ALOG(PT)
 ZUSW=EXP(CST%XALPW - CST%XBETAW/PT - CST%XGAMW*ZZW)          ! es_w
 ZZW=EXP(CST%XALPI - CST%XBETAI/PT - CST%XGAMI*ZZW)           ! es_i
END IF
  
ZSSI=0.
IF (GNEGT) THEN
  ZZW=MIN(PPABST/2., ZZW)             ! safety limitation
  ZSSI=PRVT*(PPABST-ZZW) / (CST%XEPSILO*ZZW) - 1.0
                                               ! Supersaturation over ice
  ZUSW=MIN(PPABST/2., ZUSW)            ! safety limitation
  ZUSW=(ZUSW/ZZW)*((PPABST-ZZW)/(PPABST-ZUSW)) - 1.0
                             ! Supersaturation of saturated water vapor over ice
  ZSSI=MIN(ZSSI, ZUSW) ! limitation of SSi according to SSw=0
END IF
  
ZZW=0.

IF(GNEGT) THEN
  IF(PT<CST%XTT-5.0 .AND. ZSSI>0.0) THEN
    ZZW=ICEP%XNU20*EXP(ICEP%XALPHA2*ZSSI-ICEP%XBETA2)
  ELSEIF(PT<=CST%XTT-2.0 .AND. PT>=CST%XTT-5.0 .AND. ZSSI>0.0) THEN
    ZZW=MAX(ICEP%XNU20*EXP(-ICEP%XBETA2 ), &                                                                                       
                ICEP%XNU10*EXP(-ICEP%XBETA1*(PT-CST%XTT))*(ZSSI/ZUSW)**ICEP%XALPHA1)
  ENDIF
ENDIF
IF (GNEGT) THEN
  ZZW=ZZW-PCIT
  ZZW=MIN(ZZW, 50.E3) ! limitation provisoire a 50 l^-1
END IF

PRVHENI_MR=0.

IF (GNEGT) THEN
  PRVHENI_MR=MAX(ZZW, 0.0)*ICEP%XMNU0/PRHODREF
  PRVHENI_MR=MIN(PRVT, PRVHENI_MR)
END IF
!Limitation due to 0 crossing of temperature
IF(PARAMI%LFEEDBACKT) THEN
  ZW=0.
  IF (GNEGT) THEN
    ZW=MIN(PRVHENI_MR, &
              MAX(0., (CST%XTT/PEXN-PTHT)/PLSFACT)) / &
              MAX(PRVHENI_MR, 1.E-20)
  END IF
  PRVHENI_MR=PRVHENI_MR*ZW
  ZZW=ZZW*ZW
ENDIF
IF (GNEGT) THEN
  PCIT=MAX(ZZW+PCIT, PCIT)
END IF
END SUBROUTINE ICE4_NUCLEATION
