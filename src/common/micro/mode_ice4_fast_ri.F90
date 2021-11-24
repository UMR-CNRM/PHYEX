!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_ICE4_FAST_RI
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RI(KSIZE, LDSOFT, PCOMPUTE, &
                       &PRHODREF, PLVFACT, PLSFACT, &
                       &PAI, PCJ, PCIT, &
                       &PSSI, &
                       &PRCT, PRIT, &
                       &PRCBERI, PA_TH, PA_RC, PA_RI)
!!
!!**  PURPOSE
!!    -------
!!      Computes the fast ri process
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the splitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!     S. Riette, 11/2021: loop instead of array syntax
!!
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_RAIN_ICE_PARAM, ONLY: X0DEPI,X2DEPI
USE MODD_RAIN_ICE_DESCR, ONLY: XDI,XLBEXI,XLBI,XRTMIN
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PAI      ! Thermodynamical function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PSSI     ! Supersaturation over ice
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCBERI  ! Bergeron-Findeisen effect
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
!
!*       0.2  declaration of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JL
LOGICAL :: LMASK
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RI',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*       7.2    Bergeron-Findeisen effect: RCBERI
!
DO JL=1, KSIZE
  LMASK = PSSI(JL)>0. .AND. PRCT(JL)>XRTMIN(2) .AND. &
        & PRIT(JL)>XRTMIN(4) .AND. PCIT(JL)>0. .AND. &
        & PCOMPUTE(JL)==1
  IF(LMASK) THEN
    IF(.NOT. LDSOFT) THEN
      PRCBERI(JL) = MIN(1.E8, XLBI*(PRHODREF(JL)*PRIT(JL)/PCIT(JL))**XLBEXI) ! Lbda_i
      PRCBERI(JL) = ( PSSI(JL) / (PRHODREF(JL)*PAI(JL)) ) * PCIT(JL) * &
                    ( X0DEPI/PRCBERI(JL) + X2DEPI*PCJ(JL)*PCJ(JL)/PRCBERI(JL)**(XDI+2.0) )
    ENDIF
    PA_RC(JL) = PA_RC(JL) - PRCBERI(JL)
    PA_RI(JL) = PA_RI(JL) + PRCBERI(JL)
    PA_TH(JL) = PA_TH(JL) + PRCBERI(JL)*(PLSFACT(JL)-PLVFACT(JL))
  ELSE
    PRCBERI(JL) = 0.
  ENDIF
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RI', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_FAST_RI
END MODULE MODE_ICE4_FAST_RI
