!MNH_LIC Copyright 1995-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!--------------------------------------------------------------------------
! ******* TO BE INCLUDED IN THE *CONTAINS* OF A SUBROUTINE OR FUNCTION ******
! => Don't use drHook !!!
! The GAMMA function defined here hides the intrinsic function of the same name.
! GAMMA_DBLE is the same function for double precision arguments (whatever the
! precision of the default REAL).
!
!     ######################################
      ELEMENTAL FUNCTION GAMMA(PX)  RESULT(PGAMMA)
!     ######################################
!
!
!!****  *GAMMA * -  Gamma  function
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
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Press, Teukolsky, Vetterling and Flannery: Numerical Recipes, 206-207
!!
!!    AUTHOR
!!    ------
!!      Jean-Pierre Pinty *LA/OMP*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     7/11/95
!!      C. Barthe    9/11/09  add a function for 1D arguments
!  P. Wautelet 22/06/2022: GAMMA_X0D is now declared PURE
!  S. Riette   Nov 2026: elemental function in gamma.h (replaces GAMMA_X0D and GAMMA_X1D)
!
!*       0. DECLARATIONS
!           ------------
!
IMPLICIT NONE



!
!*       0.1 declarations of arguments and result
!
REAL, INTENT(IN)                     :: PX
REAL                                 :: PGAMMA
!
!*       0.2 declarations of local variables
!
INTEGER                              :: JJ ! Loop index
REAL                                 :: ZSER,ZSTP,ZTMP,ZX,ZY
REAL, DIMENSION(6), PARAMETER        :: ZCOEF=(/76.18009172947146, &
                                               &-86.50532032941677, &
                                               & 24.01409824083091, &
                                               & -1.231739572450155, &
                                               &  0.1208650973866179E-2, &
                                               & -0.5395239384953E-5/)
REAL                                 :: ZPI
!
!-------------------------------------------------------------------------------
!
!*       1. SOME CONSTANTS
!           --------------
!
ZSTP     =  2.5066282746310005
!
ZPI = 3.141592654
!
!-------------------------------------------------------------------------------
!
!*       2. COMPUTE GAMMA
!           -------------
!
IF (PX .LT. 0.) THEN
  ZX = 1. - PX
ELSE
  ZX = PX
END IF
ZY = ZX
ZTMP =  ZX + 5.5
ZTMP = (ZX + 0.5) * LOG(ZTMP) - ZTMP
ZSER = 1.000000000190015
!
DO JJ = 1, 6
  ZY = ZY + 1.0
  ZSER = ZSER + ZCOEF(JJ) / ZY
END DO
!
IF (PX .LT. 0.) THEN
  PGAMMA = ZPI / SIN(ZPI*PX) / EXP(ZTMP + LOG(ZSTP*ZSER/ZX))
ELSE
  PGAMMA = EXP(ZTMP + LOG(ZSTP*ZSER/ZX))
END IF
!
END FUNCTION GAMMA
!
!-------------------------------------------------------------------------------
!
!     ######################################
      ELEMENTAL FUNCTION GAMMA_DBLE(PX)  RESULT(PGAMMA)
!     ######################################
!
!!****  *GAMMA_DBLE * -  Gamma  function for double precision arguments
!!
!!    PURPOSE
!!    -------
!       Same computation as the GAMMA function, done in double precision whatever
!     the precision of the default REAL (the results are identical to those of
!     GAMMA when the default REAL is in double precision)
!
!*       0. DECLARATIONS
!           ------------
!
USE MODD_PRECISION, ONLY: MNHREAL64
!
IMPLICIT NONE



!
!*       0.1 declarations of arguments and result
!
REAL(KIND=MNHREAL64), INTENT(IN)          :: PX
REAL(KIND=MNHREAL64)                      :: PGAMMA
!
!*       0.2 declarations of local variables
!
INTEGER                              :: JJ ! Loop index
REAL(KIND=MNHREAL64)                      :: ZSER,ZSTP,ZTMP,ZX,ZY
REAL(KIND=MNHREAL64), DIMENSION(6), PARAMETER :: ZCOEF=(/76.18009172947146_MNHREAL64, &
                                                   &-86.50532032941677_MNHREAL64, &
                                                   & 24.01409824083091_MNHREAL64, &
                                                   & -1.231739572450155_MNHREAL64, &
                                                   &  0.1208650973866179E-2_MNHREAL64, &
                                                   & -0.5395239384953E-5_MNHREAL64/)
REAL(KIND=MNHREAL64)                      :: ZPI
!
!-------------------------------------------------------------------------------
!
!*       1. SOME CONSTANTS
!           --------------
!
ZSTP     =  2.5066282746310005_MNHREAL64
!
ZPI = 3.141592654_MNHREAL64
!
!-------------------------------------------------------------------------------
!
!*       2. COMPUTE GAMMA
!           -------------
!
IF (PX .LT. 0._MNHREAL64) THEN
  ZX = 1._MNHREAL64 - PX
ELSE
  ZX = PX
END IF
ZY = ZX
ZTMP =  ZX + 5.5_MNHREAL64
ZTMP = (ZX + 0.5_MNHREAL64) * LOG(ZTMP) - ZTMP
ZSER = 1.000000000190015_MNHREAL64
!
DO JJ = 1, 6
  ZY = ZY + 1.0_MNHREAL64
  ZSER = ZSER + ZCOEF(JJ) / ZY
END DO
!
IF (PX .LT. 0._MNHREAL64) THEN
  PGAMMA = ZPI / SIN(ZPI*PX) / EXP(ZTMP + LOG(ZSTP*ZSER/ZX))
ELSE
  PGAMMA = EXP(ZTMP + LOG(ZSTP*ZSER/ZX))
END IF
!
END FUNCTION GAMMA_DBLE
