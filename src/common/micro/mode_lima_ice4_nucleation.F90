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
!$mnh_expand_where(II=1:KSIZE)
GNEGT(1:KSIZE)=PT(1:KSIZE)<CST%XTT .AND. PRVT(1:KSIZE)>LIMAP%XRTMIN(1)
!$mnh_end_expand_where(II=1:KSIZE)

ZUSW(:)=0.
ZZW(:)=0.
!$mnh_expand_where(II=1:KSIZE)
WHERE(GNEGT(1:KSIZE))
  ZZW(1:KSIZE)=LOG(PT(1:KSIZE))
  ZUSW(1:KSIZE)=EXP(CST%XALPW - CST%XBETAW/PT(1:KSIZE) - CST%XGAMW*ZZW(1:KSIZE))          ! es_w
  ZZW(1:KSIZE)=EXP(CST%XALPI - CST%XBETAI/PT(1:KSIZE) - CST%XGAMI*ZZW(1:KSIZE))           ! es_i
END WHERE
!$mnh_end_expand_where(II=1:KSIZE)

ZSSI(:)=0.
!$mnh_expand_where(II=1:KSIZE)
WHERE(GNEGT(1:KSIZE))
  ZZW(1:KSIZE)=MIN(PPABST(1:KSIZE)/2., ZZW(1:KSIZE))             ! safety limitation
  ZSSI(1:KSIZE)=PRVT(1:KSIZE)*(PPABST(1:KSIZE)-ZZW(1:KSIZE)) / (CST%XEPSILO*ZZW(1:KSIZE)) - 1.0
                                               ! Supersaturation over ice
  ZUSW(1:KSIZE)=MIN(PPABST(1:KSIZE)/2., ZUSW(1:KSIZE))            ! safety limitation
  ZUSW(1:KSIZE)=(ZUSW(1:KSIZE)/ZZW(1:KSIZE))*((PPABST(1:KSIZE)-ZZW(1:KSIZE))/(PPABST(1:KSIZE)-ZUSW(1:KSIZE))) - 1.0
                             ! Supersaturation of saturated water vapor over ice
  !
  !*       3.1     compute the heterogeneous nucleation source RVHENI
  !
  !*       3.1.1   compute the cloud ice concentration
  !
  ZSSI(1:KSIZE)=MIN(ZSSI(1:KSIZE), ZUSW(1:KSIZE)) ! limitation of SSi according to SSw=0
END WHERE
!$mnh_end_expand_where(II=1:KSIZE)

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
!$mnh_expand_where(II=1:KSIZE)
WHERE(GNEGT(1:KSIZE))
  ! convert between m-3 (ICE3) and kg-1 (LIMA)
  ZZW(1:KSIZE)=ZZW(1:KSIZE)-PCIT(1:KSIZE)*PRHODREF(1:KSIZE)
  ZZW(1:KSIZE)=MIN(ZZW(1:KSIZE), 50.E3) ! limitation provisoire a 50 l^-1
END WHERE
!$mnh_end_expand_where(II=1:KSIZE)

PRVHENI_MR(:)=0.
!$mnh_expand_where(II=1:KSIZE)
WHERE(GNEGT(1:KSIZE))
  !
  !*       3.1.2   update the r_i and r_v mixing ratios
  !
  PRVHENI_MR(1:KSIZE)=MAX(ZZW(1:KSIZE), 0.0)*LIMAC%XMNU0/PRHODREF(1:KSIZE)
  PRVHENI_MR(1:KSIZE)=MIN(PRVT(1:KSIZE), PRVHENI_MR(1:KSIZE))
END WHERE
!$mnh_end_expand_where(II=1:KSIZE)
!Limitation due to 0 crossing of temperature
IF(LIMAP%LFEEDBACKT) THEN
  ZW(:)=0.
  !$mnh_expand_where(II=1:KSIZE)
  WHERE(GNEGT(1:KSIZE))
    ZW(1:KSIZE)=MIN(PRVHENI_MR(1:KSIZE), &
              MAX(0., (CST%XTT/PEXN(1:KSIZE)-PTHT(1:KSIZE))/PLSFACT(1:KSIZE))) / &
              MAX(PRVHENI_MR(1:KSIZE), 1.E-20)
  END WHERE
  PRVHENI_MR(1:KSIZE)=PRVHENI_MR(1:KSIZE)*ZW(1:KSIZE)
  ZZW(1:KSIZE)=ZZW(1:KSIZE)*ZW(1:KSIZE)
  !$mnh_end_expand_where(II=1:KSIZE)
ENDIF
!$mnh_expand_where(II=1:KSIZE)
WHERE(GNEGT(1:KSIZE))
  ! convert from m-3 (ICE3) to kg-1 (LIMA)
  PCIT(1:KSIZE)=MAX(ZZW(1:KSIZE)/PRHODREF(1:KSIZE)+PCIT(1:KSIZE), PCIT(1:KSIZE))
END WHERE
!$mnh_end_expand_where(II=1:KSIZE)
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE4_NUCLEATION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_ICE4_NUCLEATION
END MODULE MODE_LIMA_ICE4_NUCLEATION
