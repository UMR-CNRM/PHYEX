!     ######spl
      MODULE MODD_PARAMETERS
!     ######################
!
!!****  *MODD_PARAMETERS* - declaration of parameter variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the variables
!     which have the PARAMETER attribute
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_PARAMETER)
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    4/07/94
!!      Modification 10/03/95 (I.Mallet)   add the coupling files maximum number
!!      Modification 10/04/95 (Ph. Hereil) add the budget related informations
!!      Modification 15/03/99 (V. Masson)  add default value
!!      Modification 17/11/00 (P.Jabouille) add the dummy array size
!!      Modification 22/01/01 (D.Gazen) change JPSVMAX from 100 to 200
!!                                         and JPBUMAX from 120 to 250
!!      Modification 17/05/04 (P.Jabouille) add JPOUTMAX
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: JPHEXT = 0      ! Horizontal External points number
INTEGER, PARAMETER :: JPVEXT = 0      ! Vertical External points number
INTEGER, PARAMETER :: JPVEXT_TURB = 1      ! Vertical External points number
INTEGER, PARAMETER :: JPMODELMAX = 8  ! Maximum allowed number of nested models
INTEGER, PARAMETER :: JPCPLFILEMAX = 24 ! Maximum allowed number of CouPLing FILEs
INTEGER, PARAMETER :: JPBUMAX= 250     ! Maximum of allowed budgets 
INTEGER, PARAMETER :: JPBUPROMAX = 60 ! Maximum of allowed processes for all
                                      ! budgets
INTEGER, PARAMETER :: JPRIMMAX = 6    ! Maximum number of points for the
                       ! horizontal relaxation for the outermost verticals
INTEGER, PARAMETER :: JPSVMAX  = 200  ! Maximum number of scalar variables
!
!
REAL,    PARAMETER :: XUNDEF = 1.E+20   ! default value for undefined or unused
!                                     ! field.
INTEGER, PARAMETER :: NUNDEF = 1E+9    ! default value for undefined or unused
!                                     ! field.
INTEGER, PARAMETER :: JPDUMMY  = 20   ! Size of dummy array
!
INTEGER, PARAMETER :: JPOUTMAX = 192 ! Maximum allowed number of OUTput files
!
INTEGER, PARAMETER :: JPLIMACCNMAX = 10 ! Maximum allowed number of CCN modes in LIMA
INTEGER, PARAMETER :: JPLIMAIFNMAX = 10 ! Maximum allowed number of IFN modes in LIMA
!
END MODULE MODD_PARAMETERS
