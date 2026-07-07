!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_RRHONG
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_RRHONG(CST, PARAMI, ICED, D, LDCOMPUTE, &
                       &PEXN, PLVFACT, PLSFACT, &
                       &PT,   PRRT, &
                       &PTHT, &
                       &PRRHONG_MR)

!$ACDC singlecolumn --dummy

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
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
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
TYPE(DIMPHYEX_t),INTENT(IN) :: D
LOGICAL, DIMENSION(D%NIJT),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PEXN     ! Exner function
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PTHT     ! Theta at t
REAL, DIMENSION(D%NIJT),       INTENT(OUT)   :: PRRHONG_MR ! Mixing ratio change due to spontaneous freezing
!
!*       0.2  declaration of local variables
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JIJ
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_RRHONG',0,ZHOOK_HANDLE)
!
!*       3.3     compute the spontaneous freezing source: RRHONG
!


DO JIJ=D%NIJB, D%NIJE
  IF(PT(JIJ)<CST%XTT-35.0 .AND. PRRT(JIJ)>ICED%XRTMIN(3) .AND. LDCOMPUTE(JIJ)) THEN
    PRRHONG_MR(JIJ)=PRRT(JIJ)
    IF(PARAMI%LFEEDBACKT) THEN
      !Limitation due to -35 crossing of temperature
      PRRHONG_MR(JIJ)=MIN(PRRHONG_MR(JIJ), MAX(0., ((CST%XTT-35.)/PEXN(JIJ)-PTHT(JIJ))/(PLSFACT(JIJ)-PLVFACT(JIJ))))
    ENDIF
  ELSE
    PRRHONG_MR(JIJ)=0.
  ENDIF
ENDDO

!
IF (LHOOK) CALL DR_HOOK('ICE4_RRHONG', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_RRHONG
END MODULE MODE_ICE4_RRHONG
