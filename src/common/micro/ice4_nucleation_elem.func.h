!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
ELEMENTAL SUBROUTINE ICE4_NUCLEATION_ELEM(ODCOMPUTE, &
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
USE MODD_CST,            ONLY: XALPI, XALPW, XBETAI, XBETAW, XGAMI, XGAMW, XMD, XMV, XTT, XEPSILO
USE MODD_PARAM_ICE,      ONLY: LFEEDBACKT
USE MODD_RAIN_ICE_PARAM, ONLY: XALPHA1, XALPHA2, XBETA1, XBETA2, XMNU0, XNU10, XNU20
USE MODD_RAIN_ICE_DESCR, ONLY: XRTMIN
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
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
GNEGT=PT<XTT .AND. PRVT>XRTMIN(1) .AND. ODCOMPUTE

PRVHENI_MR=0.
IF(GNEGT) THEN
  ZZW=ALOG(PT)
  ZUSW=EXP(XALPW - XBETAW/PT - XGAMW*ZZW)          ! es_w
  ZZW=EXP(XALPI - XBETAI/PT - XGAMI*ZZW)           ! es_i

  ZZW=MIN(PPABST/2., ZZW)             ! safety limitation
  ZSSI=PRVT*(PPABST-ZZW) / (XEPSILO*ZZW) - 1.0 ! Supersaturation over ice
  ZUSW=MIN(PPABST/2., ZUSW)            ! safety limitation
  ZUSW=(ZUSW/ZZW)*((PPABST-ZZW)/(PPABST-ZUSW)) - 1.0
                             ! Supersaturation of saturated water vapor over ice
  !
  !*       3.1     compute the heterogeneous nucleation source RVHENI
  !
  !*       3.1.1   compute the cloud ice concentration
  !
  ZSSI=MIN(ZSSI, ZUSW) ! limitation of SSi according to SSw=0

  IF(PT<XTT-5. .AND. ZSSI>0.) THEN
    ZZW=XNU20*EXP(XALPHA2*ZSSI-XBETA2)
  ELSEIF(PT<=XTT-2. .AND. PT>=XTT-5. .AND. ZSSI>0.) THEN
    ZZW=MAX(XNU20*EXP(-XBETA2 ), &
            XNU10*EXP(-XBETA1*(PT-XTT))*(ZSSI/ZUSW)**XALPHA1)
  ELSE
    ZZW=0.
  ENDIF

  ZZW=ZZW-PCIT
  ZZW=MIN(ZZW, 50.E3) ! limitation provisoire a 50 l^-1
  !
  !*       3.1.2   update the r_i and r_v mixing ratios
  !
  PRVHENI_MR=MAX(ZZW, 0.0)*XMNU0/PRHODREF
  PRVHENI_MR=MIN(PRVT, PRVHENI_MR)
  !
  !Limitation due to 0 crossing of temperature
  !
  IF(LFEEDBACKT) THEN
    ZW=MIN(PRVHENI_MR, MAX(0., (XTT/PEXN-PTHT)/PLSFACT)) / &
       MAX(PRVHENI_MR, 1.E-20)
    PRVHENI_MR=PRVHENI_MR*ZW
    ZZW=ZZW*ZW
  ENDIF
  !
  PCIT=MAX(ZZW+PCIT, PCIT)
ENDIF
!
END SUBROUTINE ICE4_NUCLEATION_ELEM
