!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_RIMLTC
IMPLICIT NONE
CONTAINS

SUBROUTINE ICE4_RIMLTC(CST, PARAMI, D, LDCOMPUTE, &
                       &PEXN, PLVFACT, PLSFACT, &
                       &PT, &
                       &PTHT, PRIT, &
                       &PRIMLTC_MR)

!$ACDC singlecolumn --dummy

!!
!!**  PURPOSE
!!    -------
!!      Computes the RIMLTC process
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
USE MODD_DIMPHYEX,  ONLY: DIMPHYEX_t
USE MODD_CST,       ONLY: CST_t
USE MODD_PARAM_ICE_n, ONLY: PARAM_ICE_t
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),                  INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),            INTENT(IN)    :: PARAMI
TYPE(DIMPHYEX_t),             INTENT(IN)    :: D
LOGICAL, DIMENSION(D%NIJT),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PEXN     ! Exner function
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PTHT     ! Theta at t
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PRIT     ! Cloud ice at t
REAL, DIMENSION(D%NIJT),       INTENT(OUT)   :: PRIMLTC_MR ! Mixing ratio change due to cloud ice melting
!
!*       0.2  declaration of local variables
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JIJ
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_RIMLTC',0,ZHOOK_HANDLE)
!
!*       7.1    cloud ice melting
!


DO JIJ=D%NIJB, D%NIJE
  IF(PRIT(JIJ)>0. .AND. PT(JIJ)>CST%XTT .AND. LDCOMPUTE(JIJ)) THEN
    PRIMLTC_MR(JIJ)=PRIT(JIJ)
    IF(PARAMI%LFEEDBACKT) THEN
      !Limitation due to 0 crossing of temperature
      PRIMLTC_MR(JIJ)=MIN(PRIMLTC_MR(JIJ), MAX(0., (PTHT(JIJ)-CST%XTT/PEXN(JIJ)) / (PLSFACT(JIJ)-PLVFACT(JIJ))))
    ENDIF
  ELSE
    PRIMLTC_MR(JIJ)=0.
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('ICE4_RIMLTC', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_RIMLTC
END MODULE MODE_ICE4_RIMLTC
