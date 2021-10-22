!     ######spl
      MODULE MODD_CH_JVALUES_n
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ########################
!!
!!*** *MODD_CH_JVALUES$n*
!!
!!    PURPOSE
!!    -------
!       This module contains the photolysis rates (JVALUES) calculated
!     by the radiative transfer code TUV39
!!
!!**  AUTHOR
!!    ------
!!    K. Suhre      *Laboratoire d'Aerologie*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 05/03/97
!!    14/02/99 (K. Suhre) add surface albedo and column dobson to parameter list
!!    16/02/99 (K. Suhre) add offline option
!!    20/01/01 (C. Mari) 3D interpolation of J values 
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CH_JVALUES_t
!
  REAL, DIMENSION(:,:,:,:), POINTER :: XJVALUES=>NULL()  
                              ! photolysis rates after time interpolation
  LOGICAL :: GSFIRSTCALL = .TRUE. ! flag for initialization on first call
!
!-----------------------------------------------------------------------------
END TYPE CH_JVALUES_t

TYPE(CH_JVALUES_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_JVALUES_MODEL

REAL, POINTER, DIMENSION(:,:,:,:) :: XJVALUES=>NULL()
LOGICAL, POINTER :: GSFIRSTCALL=>NULL()

CONTAINS

SUBROUTINE CH_JVALUES_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODD_CH_JVALUES_N:CH_JVALUES_GOTO_MODEL',0,ZHOOK_HANDLE)
CH_JVALUES_MODEL(KFROM)%XJVALUES=>XJVALUES
!
! Current model is set to model KTO
XJVALUES=>CH_JVALUES_MODEL(KTO)%XJVALUES
GSFIRSTCALL=>CH_JVALUES_MODEL(KTO)%GSFIRSTCALL

IF (LHOOK) CALL DR_HOOK('MODD_CH_JVALUES_N:CH_JVALUES_GOTO_MODEL',1,ZHOOK_HANDLE)
END SUBROUTINE CH_JVALUES_GOTO_MODEL

END MODULE MODD_CH_JVALUES_n
