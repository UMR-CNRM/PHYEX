!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
      FUNCTION GENERAL_GAMMA(PALPHA,PNU,PLBDA,PX)  RESULT(PGENERAL_GAMMA)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###################################################################
!
!
!!****  *GENERAL_GAMMA * -  Generalized gamma  function
!!
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the Generalized gamma
!    function of its argument.
!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODULE MODI_GAMMA : computation of the Gamma function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine CONDENS)
!!
!!
!!    AUTHOR
!!    ------
!!      Jean-Pierre Pinty *LA/OMP*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     7/11/95
!
!*       0. DECLARATIONS
!           ------------
!
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
REAL, INTENT(IN)                     :: PALPHA
REAL, INTENT(IN)                     :: PNU
REAL, INTENT(IN)                     :: PLBDA
REAL, INTENT(IN)                     :: PX
REAL                                 :: PGENERAL_GAMMA
!
!*       0.2 declarations of local variables
!
REAL                                 :: ZARG,ZPOWER
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GENERAL_GAMMA',0,ZHOOK_HANDLE)
ZARG   = PLBDA*PX
ZPOWER = PALPHA*PNU - 1.0
!
PGENERAL_GAMMA = (PALPHA/GAMMA(PNU))*(ZARG**ZPOWER)*PLBDA*EXP(-(ZARG**PALPHA))
IF (LHOOK) CALL DR_HOOK('GENERAL_GAMMA',1,ZHOOK_HANDLE)
RETURN
!
END FUNCTION GENERAL_GAMMA
