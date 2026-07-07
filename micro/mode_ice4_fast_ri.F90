!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_ICE4_FAST_RI
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RI(ICEP, ICED, D, LDSOFT, LDCOMPUTE, &
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
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),       INTENT(IN)    :: ICED
TYPE(DIMPHYEX_t),             INTENT(IN)    :: D
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(D%NIJT),   INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLSFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PAI      ! Thermodynamical function
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PSSI     ! Supersaturation over ice
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRCBERI  ! Bergeron-Findeisen effect
!
!*       0.2  declaration of local variables
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JIJ
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RI',0,ZHOOK_HANDLE)
!
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
#endif
!-------------------------------------------------------------------------------
!
!*       7.2    Bergeron-Findeisen effect: RCBERI
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PSSI(JIJ)>0. .AND. PRCT(JIJ)>ICED%XRTMIN(2) .AND. PRIT(JIJ)>ICED%XRTMIN(4) &
     .AND. PCIT(JIJ)>1.E-20 .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRCBERI(JIJ) = MIN(1.E8, ICED%XLBI*(PRHODREF(JIJ)*PRIT(JIJ)/PCIT(JIJ))**ICED%XLBEXI) ! Lbda_i
      PRCBERI(JIJ) = ( PSSI(JIJ) / (PRHODREF(JIJ)*PAI(JIJ)) ) * PCIT(JIJ) * &
                    ( ICEP%X0DEPI/PRCBERI(JIJ) + ICEP%X2DEPI*PCJ(JIJ)*PCJ(JIJ)/PRCBERI(JIJ)**(ICED%XDI+2.0) )
    ENDIF
  ELSE
    PRCBERI(JIJ) = 0.
  ENDIF
ENDDO
!$mnh_end_do()

!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RI', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_FAST_RI
END MODULE MODE_ICE4_FAST_RI
