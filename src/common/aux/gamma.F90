!MNH_LIC Copyright 1995-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!--------------------------------------------------------------------------
!
!
!*       1.   FUNCTION GAMMA FOR SCALAR VARIABLE
! 
!
!     ######################################
      PURE FUNCTION GAMMA_X0D(PX)  RESULT(PGAMMA)
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
!
!*       0. DECLARATIONS
!           ------------
!
IMPLICIT NONE

!$acc routine seq

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
RETURN
!
END FUNCTION GAMMA_X0D
!
!-------------------------------------------------------------------------------
!
!
!*       1.   FUNCTION GAMMA FOR 1D ARRAY
! 
!
!     ######################################
      PURE FUNCTION GAMMA_X1D(PX)  RESULT(PGAMMA)
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
!  P. Wautelet 22/06/2022: GAMMA_X1D is now declared PURE
!!
!-------------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!           ------------
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
REAL, DIMENSION(:), INTENT(IN)       :: PX
REAL, DIMENSION(SIZE(PX))            :: PGAMMA
!
!*       0.2 declarations of local variables
!
INTEGER                              :: JJ ! Loop index
INTEGER                              :: JI ! Loop index
REAL                                 :: ZSER, ZSTP, ZTMP, ZX, ZY
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
DO JI = 1, SIZE(PX)
  IF (PX(JI) .LT. 0.) THEN
    ZX = 1. - PX(JI)
  ELSE 
    ZX = PX(JI)
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
  IF (PX(JI) .LT. 0.) THEN
    PGAMMA = ZPI / SIN(ZPI*PX(JI)) / EXP(ZTMP + LOG(ZSTP*ZSER/ZX))
  ELSE
    PGAMMA = EXP(ZTMP + LOG(ZSTP*ZSER/ZX))
  END IF
END DO
RETURN
!
END FUNCTION GAMMA_X1D
