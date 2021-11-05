!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_ICE4_RIMLTC
INTERFACE
SUBROUTINE ICE4_RIMLTC(KSIZE, LDSOFT, PCOMPUTE, &
                       &PEXN, PLVFACT, PLSFACT, &
                       &PT, &
                       &PTHT, PRIT, &
                       &PRIMLTC_MR, PB_TH, PB_RC, PB_RI)
IMPLICIT NONE
INTEGER, INTENT(IN) :: KSIZE
LOGICAL,                  INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PEXN     ! Exner function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PTHT     ! Theta at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Cloud ice at t
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIMLTC_MR ! Mixing ratio change due to cloud ice melting
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RI
END SUBROUTINE ICE4_RIMLTC
END INTERFACE
END MODULE MODI_ICE4_RIMLTC
SUBROUTINE ICE4_RIMLTC(KSIZE, LDSOFT, PCOMPUTE, &
                       &PEXN, PLVFACT, PLSFACT, &
                       &PT, &
                       &PTHT, PRIT, &
                       &PRIMLTC_MR, PB_TH, PB_RC, PB_RI)
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
USE MODD_CST,       ONLY: XTT
USE MODD_PARAM_ICE, ONLY: LFEEDBACKT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PEXN     ! Exner function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PTHT     ! Theta at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Cloud ice at t
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIMLTC_MR ! Mixing ratio change due to cloud ice melting
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RI
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KSIZE) :: ZMASK
INTEGER :: JL
!
!-------------------------------------------------------------------------------
!
!*       7.1    cloud ice melting
!
PRIMLTC_MR(:)=0.
IF(.NOT. LDSOFT) THEN
  DO JL=1, KSIZE
    ZMASK(JL)=MAX(0., -SIGN(1., -PRIT(JL))) * & ! PRIT(:)>0.
             &MAX(0., -SIGN(1., XTT-PT(JL))) * & ! PT(:)>XTT
             &PCOMPUTE(JL)
    PRIMLTC_MR(JL)=PRIT(JL) * ZMASK(JL)
  ENDDO

  IF(LFEEDBACKT) THEN
    !Limitation due to 0 crossing of temperature
    DO JL=1, KSIZE
      PRIMLTC_MR(JL)=MIN(PRIMLTC_MR(JL), MAX(0., (PTHT(JL)-XTT/PEXN(JL)) / (PLSFACT(JL)-PLVFACT(JL))))
    ENDDO
  ENDIF
ENDIF
DO JL=1, KSIZE
  PB_RC(JL) = PB_RC(JL) + PRIMLTC_MR(JL)
  PB_RI(JL) = PB_RI(JL) - PRIMLTC_MR(JL)
  PB_TH(JL) = PB_TH(JL) - PRIMLTC_MR(JL)*(PLSFACT(JL)-PLVFACT(JL))
ENDDO
!
!
END SUBROUTINE ICE4_RIMLTC
