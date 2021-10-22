!     ######spl
      SUBROUTINE CH_FCN(PTIME,PCONC,PDCDT,KMI,KVECNPT,KEQ)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ########################################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
!!*** *CH_FCN*
!!
!!    PURPOSE
!!    -------
!       calculation of first derivative for the chemical reaction mechanism
!!
!!**  METHOD
!!    ------
!!      For each prognostical chemical species the first derivative is
!!    calculated as defined by the chemical reaction mechanism.
!!    The reaction rates and other user-defined auxiliary variables are
!!    transfered in the TYPE(CCSTYPE) variable TPK%.
!!    The subroutine PRODLOSS is called in order to calculate P and L
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Karsten Suhre (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 26/07/96
!!    Modified 05/05/98: Vectorization (Vincent Crassier & KS)
!!    Modified 31/10/03: New interface for better MesoNH compilation (D. Gazen)
!!
!!----------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_PRODLOSS
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
REAL,    INTENT(IN)                          :: PTIME
INTEGER, INTENT(IN)                          :: KVECNPT
INTEGER, INTENT(IN)                          :: KEQ
REAL,    INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC
REAL,    INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PDCDT
INTEGER, INTENT(IN)                          :: KMI
!!
!!    LOCAL VARIABLES
!!    ---------------
REAL, DIMENSION(KVECNPT,KEQ)     :: ZPROD, ZLOSS
!!----------------------------------------------------------------------
!!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_FCN',0,ZHOOK_HANDLE)
CALL CH_PRODLOSS(PTIME,PCONC,ZPROD,ZLOSS,KMI,KVECNPT,KEQ)
PDCDT(:,:) = ZPROD(:,:) - PCONC(:,:) * ZLOSS(:,:)
IF (LHOOK) CALL DR_HOOK('CH_FCN',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE CH_FCN
