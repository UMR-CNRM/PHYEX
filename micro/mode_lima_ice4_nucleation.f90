!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_ICE4_NUCLEATION
IMPLICIT NONE
CONTAINS
SUBROUTINE LIMA_ICE4_NUCLEATION(CST, KSIZE, &
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
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_PARAM_LIMA_COLD, ONLY : XALPHA1, XBETA1, XALPHA2, XBETA2, XNU10, XNU20, XMNU0
USE MODD_PARAM_LIMA, ONLY: LFEEDBACKT, XRTMIN
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN)    :: CST
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
LOGICAL, DIMENSION(KSIZE) :: GNEGT  ! Test where to compute the HEN process
REAL, DIMENSION(KSIZE)  :: ZZW,      & ! Work array
                           ZUSW,     & ! Undersaturation over water
                           ZSSI        ! Supersaturation over ice
INTEGER :: JI
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE4_NUCLEATION', 0, ZHOOK_HANDLE)!
!
!$mnh_expand_where(JI=1:KSIZE)
GNEGT(:)=PT(:)<CST%XTT .AND. PRVT(:)>XRTMIN(1)
!$mnh_end_expand_where(JI=1:KSIZE)

ZUSW(:)=0.
ZZW(:)=0.
!$mnh_expand_where(JI=1:KSIZE)
WHERE(GNEGT(:))
  ZZW(:)=ALOG(PT(:))
  ZUSW(:)=EXP(CST%XALPW - CST%XBETAW/PT(:) - CST%XGAMW*ZZW(:))          ! es_w
  ZZW(:)=EXP(CST%XALPI - CST%XBETAI/PT(:) - CST%XGAMI*ZZW(:))           ! es_i
END WHERE
!$mnh_end_expand_where(JI=1:KSIZE)

ZSSI(:)=0.
!$mnh_expand_where(JI=1:KSIZE)
WHERE(GNEGT(:))
  ZZW(:)=MIN(PPABST(:)/2., ZZW(:))             ! safety limitation
  ZSSI(:)=PRVT(:)*(PPABST(:)-ZZW(:)) / (CST%XEPSILO*ZZW(:)) - 1.0
                                               ! Supersaturation over ice
  ZUSW(:)=MIN(PPABST(:)/2., ZUSW(:))            ! safety limitation
  ZUSW(:)=(ZUSW(:)/ZZW(:))*((PPABST(:)-ZZW(:))/(PPABST(:)-ZUSW(:))) - 1.0
                             ! Supersaturation of saturated water vapor over ice
  !
  !*       3.1     compute the heterogeneous nucleation source RVHENI
  !
  !*       3.1.1   compute the cloud ice concentration
  !
  ZSSI(:)=MIN(ZSSI(:), ZUSW(:)) ! limitation of SSi according to SSw=0
END WHERE
!$mnh_end_expand_where(JI=1:KSIZE)

ZZW(:)=0.
DO JI=1,KSIZE
  IF(GNEGT(JI)) THEN
    IF(PT(JI)<CST%XTT-5.0 .AND. ZSSI(JI)>0.0) THEN
      ZZW(JI)=XNU20*EXP(XALPHA2*ZSSI(JI)-XBETA2)
    ELSEIF(PT(JI)<=CST%XTT-2.0 .AND. PT(JI)>=CST%XTT-5.0 .AND. ZSSI(JI)>0.0) THEN
      ZZW(JI)=MAX(XNU20*EXP(-XBETA2 ), &                                                                                       
                  XNU10*EXP(-XBETA1*(PT(JI)-CST%XTT))*(ZSSI(JI)/ZUSW(JI))**XALPHA1)
    ENDIF
  ENDIF
ENDDO
!$mnh_expand_where(JI=1:KSIZE)
WHERE(GNEGT(:))
  ZZW(:)=ZZW(:)-PCIT(:)
  ZZW(:)=MIN(ZZW(:), 50.E3) ! limitation provisoire a 50 l^-1
END WHERE
!$mnh_end_expand_where(JI=1:KSIZE)

PRVHENI_MR(:)=0.
!$mnh_expand_where(JI=1:KSIZE)
WHERE(GNEGT(:))
  !
  !*       3.1.2   update the r_i and r_v mixing ratios
  !
  PRVHENI_MR(:)=MAX(ZZW(:), 0.0)*XMNU0/PRHODREF(:)
  PRVHENI_MR(:)=MIN(PRVT(:), PRVHENI_MR(:))
END WHERE
!$mnh_end_expand_where(JI=1:KSIZE)
!Limitation due to 0 crossing of temperature
IF(LFEEDBACKT) THEN
  ZW(:)=0.
  !$mnh_expand_where(JI=1:KSIZE)
  WHERE(GNEGT(:))
    ZW(:)=MIN(PRVHENI_MR(:), &
              MAX(0., (CST%XTT/PEXN(:)-PTHT(:))/PLSFACT(:))) / &
              MAX(PRVHENI_MR(:), 1.E-20)
  END WHERE
  PRVHENI_MR(:)=PRVHENI_MR(:)*ZW(:)
  ZZW(:)=ZZW(:)*ZW(:)
  !$mnh_end_expand_where(JI=1:KSIZE)
ENDIF
!$mnh_expand_where(JI=1:KSIZE)
WHERE(GNEGT(:))
  PCIT(:)=MAX(ZZW(:)+PCIT(:), PCIT(:))
END WHERE
!$mnh_end_expand_where(JI=1:KSIZE)
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE4_NUCLEATION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_ICE4_NUCLEATION
END MODULE MODE_LIMA_ICE4_NUCLEATION
