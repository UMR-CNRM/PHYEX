!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_RIMLTC
IMPLICIT NONE
CONTAINS

SUBROUTINE ICE4_RIMLTC(CST, PARAMI, KSIZE, PCOMPUTE, &
                       &PEXN, PLVFACT, PLSFACT, &
                       &PT, &
                       &PTHT, PRIT, &
                       &PRIMLTC_MR)
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
USE MODD_CST,       ONLY: CST_t
USE MODD_PARAM_ICE, ONLY: PARAM_ICE_t
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),                  INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),            INTENT(IN)    :: PARAMI
INTEGER,                      INTENT(IN)    :: KSIZE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PEXN     ! Exner function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PTHT     ! Theta at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Cloud ice at t
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIMLTC_MR ! Mixing ratio change due to cloud ice melting
!
!*       0.2  declaration of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL, DIMENSION(KSIZE) :: ZMASK
INTEGER :: JL
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_RIMLTC',0,ZHOOK_HANDLE)
!
!*       7.1    cloud ice melting
!
PRIMLTC_MR(:)=0.
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., -PRIT(JL))) * & ! PRIT(:)>0.
           &MAX(0., -SIGN(1., CST%XTT-PT(JL))) * & ! PT(:)>XTT
           &PCOMPUTE(JL)
  PRIMLTC_MR(JL)=PRIT(JL) * ZMASK(JL)
ENDDO

IF(PARAMI%LFEEDBACKT) THEN
  !Limitation due to 0 crossing of temperature
  DO JL=1, KSIZE
    PRIMLTC_MR(JL)=MIN(PRIMLTC_MR(JL), MAX(0., (PTHT(JL)-CST%XTT/PEXN(JL)) / (PLSFACT(JL)-PLVFACT(JL))))
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('ICE4_RIMLTC', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_RIMLTC
END MODULE MODE_ICE4_RIMLTC
