!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODI_ICE4_FAST_RI
INTERFACE
SUBROUTINE ICE4_FAST_RI(KSIZE, LDSOFT, PCOMPUTE, &
                       &PRHODREF, PLVFACT, PLSFACT, &
                       &PAI, PCJ, PCIT, &
                       &PSSI, &
                       &PRCT, PRIT, &
                       &PRCBERI, PA_TH, PA_RC, PA_RI)
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_RAIN_ICE_PARAM
USE MODD_RAIN_ICE_DESCR
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
END SUBROUTINE ICE4_FAST_RI
END INTERFACE
END MODULE MODI_ICE4_FAST_RI
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
!!
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_RAIN_ICE_DESCR, ONLY: XDI,XLBEXI,XLBI,XRTMIN
USE MODD_RAIN_ICE_PARAM, ONLY: X0DEPI,X2DEPI
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
REAL, DIMENSION(KSIZE) :: ZMASK
INTEGER :: JL
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!*       7.2    Bergeron-Findeisen effect: RCBERI
!
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., -PSSI(JL))) * &          ! PSSI(:)>0.
           &MAX(0., -SIGN(1., XRTMIN(2)-PRCT(JL))) * & ! PRCT(:)>XRTMIN(2)
           &MAX(0., -SIGN(1., XRTMIN(4)-PRIT(JL))) * & ! PRIT(:)>XRTMIN(4)
           &MAX(0., -SIGN(1., 1.E-20-PCIT(JL))) * &          ! PCIT(:)>0.
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRCBERI(JL) = PRCBERI(JL) * ZMASK(JL)
  ENDDO
ELSE
  PRCBERI(:) = 0.
  WHERE(ZMASK(:)==1.)
    PRCBERI(:) = MIN(1.E8, XLBI*(PRHODREF(:)*PRIT(:)/PCIT(:))**XLBEXI) ! Lbda_i
    PRCBERI(:) = ( PSSI(:) / (PRHODREF(:)*PAI(:)) ) * PCIT(:) * &
                 ( X0DEPI/PRCBERI(:) + X2DEPI*PCJ(:)*PCJ(:)/PRCBERI(:)**(XDI+2.0) )
  END WHERE
ENDIF
DO JL=1, KSIZE
  PA_RC(JL) = PA_RC(JL) - PRCBERI(JL)
  PA_RI(JL) = PA_RI(JL) + PRCBERI(JL)
  PA_TH(JL) = PA_TH(JL) + PRCBERI(JL)*(PLSFACT(JL)-PLVFACT(JL))
ENDDO
!
!
END SUBROUTINE ICE4_FAST_RI
