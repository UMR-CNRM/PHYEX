!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_RRHONG
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_RRHONG(CST, PARAMI, ICED, KPROMA, KSIZE, LDCOMPUTE, &
                       &PEXN, PLVFACT, PLSFACT, &
                       &PT,   PRRT, &
                       &PTHT, &
                       &PRRHONG_MR)
!!
!!**  PURPOSE
!!    -------
!!      Computes the RRHONG process
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the splitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),        INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
INTEGER, INTENT(IN) :: KPROMA, KSIZE
LOGICAL, DIMENSION(KPROMA),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PEXN     ! Exner function
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PTHT     ! Theta at t
REAL, DIMENSION(KPROMA),       INTENT(OUT)   :: PRRHONG_MR ! Mixing ratio change due to spontaneous freezing
!
!*       0.2  declaration of local variables
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JL
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_RRHONG',0,ZHOOK_HANDLE)
!
!*       3.3     compute the spontaneous freezing source: RRHONG
!
DO JL=1, KSIZE
  IF(PT(JL)<CST%XTT-35.0 .AND. PRRT(JL)>ICED%XRTMIN(3) .AND. LDCOMPUTE(JL)) THEN
    PRRHONG_MR(JL)=PRRT(JL)
    IF(PARAMI%LFEEDBACKT) THEN
      !Limitation due to -35 crossing of temperature
      PRRHONG_MR(JL)=MIN(PRRHONG_MR(JL), MAX(0., ((CST%XTT-35.)/PEXN(JL)-PTHT(JL))/(PLSFACT(JL)-PLVFACT(JL))))
    ENDIF
  ELSE
    PRRHONG_MR(JL)=0.
  ENDIF
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ICE4_RRHONG', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_RRHONG
END MODULE MODE_ICE4_RRHONG
