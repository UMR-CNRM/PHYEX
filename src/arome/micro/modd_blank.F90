!     ######spl
      MODULE MODD_BLANK
!     #################
!
!!****  *MODD_BLANK* -  Declarative module for MesoNH developpers namelist
!!
!!    PURPOSE
!!    -------
!!
!!      Offer dummy real, integer, logical and character variables for
!!    test and debugging purposes.
!!
!!**  METHOD
!!    ------
!!
!!      Eight dummy real, integer, logical and character*80 variables are
!!    defined and passed through the namelist read operations. None of the
!!    MesoNH routines uses any of those variables. When a developper choses
!!    to introduce temporarily a parameter to some subroutine, he has to
!!    introduce a USE MODD_BLANK statement into that subroutine. Then he
!!    can use any of the variables defined here and change them easily via
!!    the namelist input.
!!
!!    REFERENCE
!!    ---------
!!      None
!!
!!    AUTHOR
!!    ------
!!      K. Suhre   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    Original 25/04/96
!!      updated     17/11/00  (P Jabouille) Use dummy array
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY : JPDUMMY
!
IMPLICIT NONE
!
REAL, SAVE         :: XDUMMY1, XDUMMY2, XDUMMY3, XDUMMY4, &
                      XDUMMY5, XDUMMY6, XDUMMY7, XDUMMY8
INTEGER, SAVE      :: NDUMMY1, NDUMMY2, NDUMMY3, NDUMMY4, &
                      NDUMMY5, NDUMMY6, NDUMMY7, NDUMMY8
LOGICAL, SAVE      :: LDUMMY1, LDUMMY2, LDUMMY3, LDUMMY4, &
                      LDUMMY5, LDUMMY6, LDUMMY7, LDUMMY8
CHARACTER*80, SAVE :: CDUMMY1, CDUMMY2, CDUMMY3, CDUMMY4, &
                      CDUMMY5, CDUMMY6, CDUMMY7, CDUMMY8
!
REAL,    SAVE, DIMENSION(JPDUMMY) :: XDUMMY
INTEGER, SAVE, DIMENSION(JPDUMMY) :: NDUMMY
LOGICAL, SAVE, DIMENSION(JPDUMMY) :: LDUMMY
CHARACTER*80, SAVE, DIMENSION(JPDUMMY) :: CDUMMY
!
END MODULE MODD_BLANK
