!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_ICE4_FAST_RI
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RI(ICEP, ICED, KPROMA, KSIZE, LDSOFT, LDCOMPUTE, &
                       &PRHODREF, PLVFACT, PLSFACT, &
                       &PAI, PCJ, PCIT, &
                       &PSSI, &
                       &PRCT, PRIT, &
                       &PRCBERI)
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
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),       INTENT(IN)    :: ICED
INTEGER,                      INTENT(IN)    :: KPROMA, KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KPROMA),   INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PAI      ! Thermodynamical function
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PSSI     ! Supersaturation over ice
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRCBERI  ! Bergeron-Findeisen effect
!
!*       0.2  declaration of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JL
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RI',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*       7.2    Bergeron-Findeisen effect: RCBERI
!
DO JL=1, KSIZE
  IF(PSSI(JL)>0. .AND. PRCT(JL)>ICED%XRTMIN(2) .AND. PRIT(JL)>ICED%XRTMIN(4) &
#ifdef REPRO48
     .AND. PCIT(JL)>0. .AND. LDCOMPUTE(JL)) THEN
#else
     .AND. PCIT(JL)>1.E-20 .AND. LDCOMPUTE(JL)) THEN
#endif
    IF(.NOT. LDSOFT) THEN
      PRCBERI(JL) = MIN(1.E8, ICED%XLBI*(PRHODREF(JL)*PRIT(JL)/PCIT(JL))**ICED%XLBEXI) ! Lbda_i
      PRCBERI(JL) = ( PSSI(JL) / (PRHODREF(JL)*PAI(JL)) ) * PCIT(JL) * &
                    ( ICEP%X0DEPI/PRCBERI(JL) + ICEP%X2DEPI*PCJ(JL)*PCJ(JL)/PRCBERI(JL)**(ICED%XDI+2.0) )
    ENDIF
  ELSE
    PRCBERI(JL) = 0.
  ENDIF
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RI', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_FAST_RI
END MODULE MODE_ICE4_FAST_RI
