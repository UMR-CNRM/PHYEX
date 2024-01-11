!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
SUBROUTINE ICE4_NUCLEATION(CST, PARAMI, ICEP, ICED, KSIZE, ODCOMPUTE, &
                           PTHT, PPABST, PRHODREF, PEXN, PLSFACT, PT, &
                           PRVT, &
                           PCIT, PRVHENI_MR, PBUF, LDBUF)
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
INTEGER,                  INTENT(IN)    :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PTHT    ! Theta at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT      ! Temperature at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: PRVHENI_MR ! Mixing ratio change due to the heterogeneous nucleation
REAL, DIMENSION(KSIZE, 4),INTENT(OUT)   :: PBUF    ! Buffer
LOGICAL, DIMENSION(KSIZE),INTENT(OUT)   :: LDBUF   ! Buffer
!
!*       0.2  declaration of local variables
!
!This routine will be included in a CONTAIN part of other routine
!For GPU source-to-source transformation, this routine cannot declare
!local arrays. We must use buffer (declared in the calling routine)
INTEGER :: JL
INTEGER, PARAMETER :: IW1=1, IW2=2, IUSW=3, ISSI=4
!-------------------------------------------------------------------------------
!
!
!$mnh_expand_where(JL=1:KSIZE)
WHERE(ODCOMPUTE(:))
  LDBUF(:)=PT(:)<CST%XTT .AND. PRVT(:)>ICED%XRTMIN(1)
ELSEWHERE
  LDBUF(:)=.FALSE.
ENDWHERE
!$mnh_end_expand_where(JL=1:KSIZE)

PBUF(:, IUSW)=0.
PBUF(:, IW2)=0.
!$mnh_expand_where(JL=1:KSIZE)
WHERE(LDBUF(:))
  PBUF(:, IW2)=ALOG(PT(:))
  PBUF(:, IUSW)=EXP(CST%XALPW - CST%XBETAW/PT(:) - CST%XGAMW*PBUF(:, IW2))          ! es_w
  PBUF(:, IW2)=EXP(CST%XALPI - CST%XBETAI/PT(:) - CST%XGAMI*PBUF(:, IW2))           ! es_i
END WHERE
!$mnh_end_expand_where(JL=1:KSIZE)

PBUF(:, ISSI)=0.
!$mnh_expand_where(JL=1:KSIZE)
WHERE(LDBUF(:))
  PBUF(:, IW2)=MIN(PPABST(:)/2., PBUF(:, IW2))             ! safety limitation
  PBUF(:, ISSI)=PRVT(:)*(PPABST(:)-PBUF(:, IW2)) / (CST%XEPSILO*PBUF(:, IW2)) - 1.0
                                               ! Supersaturation over ice
  PBUF(:, IUSW)=MIN(PPABST(:)/2., PBUF(:, IUSW))            ! safety limitation
  PBUF(:, IUSW)=(PBUF(:, IUSW)/PBUF(:, IW2))*((PPABST(:)-PBUF(:, IW2))/(PPABST(:)-PBUF(:, IUSW))) - 1.0
                             ! Supersaturation of saturated water vapor over ice
  PBUF(:, ISSI)=MIN(PBUF(:, ISSI), PBUF(:, IUSW)) ! limitation of SSi according to SSw=0
END WHERE
!$mnh_end_expand_where(JL=1:KSIZE)

PBUF(:, IW2)=0.
DO JL=1,KSIZE
  IF(LDBUF(JL)) THEN
    IF(PT(JL)<CST%XTT-5.0 .AND. PBUF(JL, ISSI)>0.0) THEN
      PBUF(JL, IW2)=ICEP%XNU20*EXP(ICEP%XALPHA2*PBUF(JL, ISSI)-ICEP%XBETA2)
    ELSEIF(PT(JL)<=CST%XTT-2.0 .AND. PT(JL)>=CST%XTT-5.0 .AND. PBUF(JL, ISSI)>0.0) THEN
      PBUF(JL, IW2)=MAX(ICEP%XNU20*EXP(-ICEP%XBETA2 ), &                                                                                       
                  ICEP%XNU10*EXP(-ICEP%XBETA1*(PT(JL)-CST%XTT))*(PBUF(JL, ISSI)/PBUF(JL, IUSW))**ICEP%XALPHA1)
    ENDIF
  ENDIF
ENDDO
!$mnh_expand_where(JL=1:KSIZE)
WHERE(LDBUF(:))
  PBUF(:, IW2)=PBUF(:, IW2)-PCIT(:)
  PBUF(:, IW2)=MIN(PBUF(:, IW2), 50.E3) ! limitation provisoire a 50 l^-1
END WHERE
!$mnh_end_expand_where(JL=1:KSIZE)

PRVHENI_MR(:)=0.
!$mnh_expand_where(JL=1:KSIZE)
WHERE(LDBUF(:))
  PRVHENI_MR(:)=MAX(PBUF(:, IW2), 0.0)*ICEP%XMNU0/PRHODREF(:)
  PRVHENI_MR(:)=MIN(PRVT(:), PRVHENI_MR(:))
END WHERE
!$mnh_end_expand_where(JL=1:KSIZE)
!Limitation due to 0 crossing of temperature
IF(PARAMI%LFEEDBACKT) THEN
  PBUF(:, IW1)=0.
  !$mnh_expand_where(JL=1:KSIZE)
  WHERE(LDBUF(:))
    PBUF(:, IW1)=MIN(PRVHENI_MR(:), &
              MAX(0., (CST%XTT/PEXN(:)-PTHT(:))/PLSFACT(:))) / &
              MAX(PRVHENI_MR(:), 1.E-20)
  END WHERE
  PRVHENI_MR(:)=PRVHENI_MR(:)*PBUF(:, IW1)
  PBUF(:, IW2)=PBUF(:, IW2)*PBUF(:, IW1)
  !$mnh_end_expand_where(JL=1:KSIZE)
ENDIF
!$mnh_expand_where(JL=1:KSIZE)
WHERE(LDBUF(:))
  PCIT(:)=MAX(PBUF(:, IW2)+PCIT(:), PCIT(:))
END WHERE
!$mnh_end_expand_where(JL=1:KSIZE)
!
END SUBROUTINE ICE4_NUCLEATION
