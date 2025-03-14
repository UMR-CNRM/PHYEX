!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_ICE4_NUCLEATION
IMPLICIT NONE
CONTAINS
SUBROUTINE LIMA_ICE4_NUCLEATION(LIMAP, LIMAC, CST, KSIZE, &
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
USE MODD_CST,            ONLY: CST_T
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_T),              INTENT(IN)    :: CST
INTEGER,                  INTENT(IN)    :: KSIZE
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PTHT    ! Theta at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT      ! Temperature at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: PRVHENI_MR ! Mixing ratio change due to the heterogeneous nucleation
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KSIZE) :: ZW ! work array
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
LOGICAL, DIMENSION(KSIZE) :: GNEGT  ! Test where to compute the HEN process
REAL, DIMENSION(KSIZE)  :: ZZW,      & ! Work array
                           ZUSW,     & ! Undersaturation over water
                           ZSSI        ! Supersaturation over ice
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
INTEGER :: II
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE4_NUCLEATION', 0, ZHOOK_HANDLE)!
!
DO II=1, KSIZE
  GNEGT(II)=PT(II)<CST%XTT .AND. PRVT(II)>LIMAP%XRTMIN(1)
END DO

ZUSW(:)=0.
ZZW(:)=0.
DO II=1, KSIZE
  IF (GNEGT(II)) THEN
    ZZW(II)=ALOG(PT(II))
    ZUSW(II)=EXP(CST%XALPW - CST%XBETAW/PT(II) - CST%XGAMW*ZZW(II))          ! es_w
    ZZW(II)=EXP(CST%XALPI - CST%XBETAI/PT(II) - CST%XGAMI*ZZW(II))           ! es_i
  END IF
END DO

ZSSI(:)=0.
DO II=1, KSIZE
  IF (GNEGT(II)) THEN
    ZZW(II)=MIN(PPABST(II)/2., ZZW(II))             ! safety limitation
    ZSSI(II)=PRVT(II)*(PPABST(II)-ZZW(II)) / (CST%XEPSILO*ZZW(II)) - 1.0
                                                 ! Supersaturation over ice
    ZUSW(II)=MIN(PPABST(II)/2., ZUSW(II))            ! safety limitation
    ZUSW(II)=(ZUSW(II)/ZZW(II))*((PPABST(II)-ZZW(II))/(PPABST(II)-ZUSW(II))) - 1.0
                               ! Supersaturation of saturated water vapor over ice
    !
    !*       3.1     compute the heterogeneous nucleation source RVHENI
    !
    !*       3.1.1   compute the cloud ice concentration
    !
    ZSSI(II)=MIN(ZSSI(II), ZUSW(II)) ! limitation of SSi according to SSw=0
  END IF
END DO

ZZW(:)=0.
DO II=1,KSIZE
  IF(GNEGT(II)) THEN
    IF(PT(II)<CST%XTT-5.0 .AND. ZSSI(II)>0.0) THEN
      ZZW(II)=LIMAC%XNU20*EXP(LIMAC%XALPHA2*ZSSI(II)-LIMAC%XBETA2)
    ELSEIF(PT(II)<=CST%XTT-2.0 .AND. PT(II)>=CST%XTT-5.0 .AND. ZSSI(II)>0.0) THEN
      ZZW(II)=MAX(LIMAC%XNU20*EXP(-LIMAC%XBETA2 ), &                                                                                       
                  LIMAC%XNU10*EXP(-LIMAC%XBETA1*(PT(II)-CST%XTT))*(ZSSI(II)/ZUSW(II))**LIMAC%XALPHA1)
    ENDIF
  ENDIF
ENDDO
DO II=1, KSIZE
  IF (GNEGT(II)) THEN
    ZZW(II)=ZZW(II)-PCIT(II)
    ZZW(II)=MIN(ZZW(II), 50.E3) ! limitation provisoire a 50 l^-1
  END IF
END DO

PRVHENI_MR(:)=0.
DO II=1, KSIZE
  IF (GNEGT(II)) THEN
    !
    !*       3.1.2   update the r_i and r_v mixing ratios
    !
    PRVHENI_MR(II)=MAX(ZZW(II), 0.0)*LIMAC%XMNU0/PRHODREF(II)
    PRVHENI_MR(II)=MIN(PRVT(II), PRVHENI_MR(II))
  END IF
END DO
!Limitation due to 0 crossing of temperature
IF(LIMAP%LFEEDBACKT) THEN
  ZW(:)=0.
  DO II=1, KSIZE
    IF (GNEGT(II)) THEN
      ZW(II)=MIN(PRVHENI_MR(II), &
                MAX(0., (CST%XTT/PEXN(II)-PTHT(II))/PLSFACT(II))) / &
                MAX(PRVHENI_MR(II), 1.E-20)
    END IF
    PRVHENI_MR(II)=PRVHENI_MR(II)*ZW(II)
    ZZW(II)=ZZW(II)*ZW(II)
  END DO
ENDIF
DO II=1, KSIZE
  IF (GNEGT(II)) THEN
    PCIT(II)=MAX(ZZW(II)+PCIT(II), PCIT(II))
  END IF
END DO
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE4_NUCLEATION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_ICE4_NUCLEATION
END MODULE MODE_LIMA_ICE4_NUCLEATION
