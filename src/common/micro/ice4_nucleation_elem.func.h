!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
ELEMENTAL SUBROUTINE ICE4_NUCLEATION_ELEM(CST, PARAMI, ICEP, ICED, ODCOMPUTE, &
                           PTHT, PPABST, PRHODREF, PEXN, PLSFACT, PT, &
                           PRVT, &
                           PCIT, PRVHENI_MR)
! ******* TO BE INCLUDED IN THE *CONTAINS* OF A SUBROUTINE, IN ORDER TO EASE AUTOMATIC INLINING ******
! => Don't use drHook !!!
!

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
!!     S. Riette Feb 2022: as an include file
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
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
REAL :: ZW ! work array
LOGICAL :: GNEGT  ! Test where to compute the HEN process
REAL  :: ZZW,      & ! Work scalar
         ZUSW,     & ! Undersaturation over water
         ZSSI        ! Supersaturation over ice
!-------------------------------------------------------------------------------
!
GNEGT=PT<CST%XTT .AND. PRVT>ICED%XRTMIN(1) .AND. ODCOMPUTE

PRVHENI_MR=0.
IF(GNEGT) THEN
  ZZW=ALOG(PT)
  ZUSW=EXP(CST%XALPW - CST%XBETAW/PT - CST%XGAMW*ZZW)          ! es_w
  ZZW=EXP(CST%XALPI - CST%XBETAI/PT - CST%XGAMI*ZZW)           ! es_i

  ZZW=MIN(PPABST/2., ZZW)             ! safety limitation
  ZSSI=PRVT*(PPABST-ZZW) / (CST%XEPSILO*ZZW) - 1.0 ! Supersaturation over ice
  ZUSW=MIN(PPABST/2., ZUSW)            ! safety limitation
  ZUSW=(ZUSW/ZZW)*((PPABST-ZZW)/(PPABST-ZUSW)) - 1.0
                             ! Supersaturation of saturated water vapor over ice
  !
  !*       3.1     compute the heterogeneous nucleation source RVHENI
  !
  !*       3.1.1   compute the cloud ice concentration
  !
  ZSSI=MIN(ZSSI, ZUSW) ! limitation of SSi according to SSw=0

  IF(PT<CST%XTT-5. .AND. ZSSI>0.) THEN
    ZZW=ICEP%XNU20*EXP(ICEP%XALPHA2*ZSSI-ICEP%XBETA2)
  ELSEIF(PT<=CST%XTT-2. .AND. PT>=CST%XTT-5. .AND. ZSSI>0.) THEN
    ZZW=MAX(ICEP%XNU20*EXP(-ICEP%XBETA2 ), &
            ICEP%XNU10*EXP(-ICEP%XBETA1*(PT-CST%XTT))*(ZSSI/ZUSW)**ICEP%XALPHA1)
  ELSE
    ZZW=0.
  ENDIF

  ZZW=ZZW-PCIT
  ZZW=MIN(ZZW, 50.E3) ! limitation provisoire a 50 l^-1
  !
  !*       3.1.2   update the r_i and r_v mixing ratios
  !
  PRVHENI_MR=MAX(ZZW, 0.0)*ICEP%XMNU0/PRHODREF
  PRVHENI_MR=MIN(PRVT, PRVHENI_MR)
  !
  !Limitation due to 0 crossing of temperature
  !
  IF(PARAMI%LFEEDBACKT) THEN
    ZW=MIN(PRVHENI_MR, MAX(0., (CST%XTT/PEXN-PTHT)/PLSFACT)) / &
       MAX(PRVHENI_MR, 1.E-20)
    PRVHENI_MR=PRVHENI_MR*ZW
    ZZW=ZZW*ZW
  ENDIF
  !
  PCIT=MAX(ZZW+PCIT, PCIT)
ENDIF
!
END SUBROUTINE ICE4_NUCLEATION_ELEM
