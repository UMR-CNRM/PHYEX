!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     ###########################
      MODULE MODD_PARAMETERS_DEP
!     ###########################
!
!!****  *MODD_PARAMETERS_DEP* - declaration of parameter variables
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
!!      01-02-2011 M. Mokhtari Adaptation modd_parameter under modd_parameters_dep for Aladin
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: JPHEXT = 0      ! Horizontal External points number
INTEGER, PARAMETER :: JPVEXT = 1      ! Vertical External points number
INTEGER, PARAMETER :: JPMODELMAX = 8  ! Maximum allowed number of nested models
INTEGER, PARAMETER :: JPCPLFILEMAX = 8 ! Maximum allowed number of CouPLing FILEs
INTEGER, PARAMETER :: JPBUMAX= 250     ! Maximum of allowed budgets
INTEGER, PARAMETER :: JPBUPROMAX = 40 ! Maximum of allowed processes for all
                                      ! budgets
INTEGER, PARAMETER :: JPRIMMAX = 6    ! Maximum number of points for the
                       ! horizontal relaxation for the outermost verticals
INTEGER, PARAMETER :: JPSVMAX  = 200  ! Maximum number of scalar variables


REAL,    PARAMETER :: XUNDEF = 1.E+20   ! default value for undefined or unused
                                      ! field.
INTEGER, PARAMETER :: NUNDEF = 1E+9    ! default value for undefined or unused
                                      ! field.
INTEGER, PARAMETER :: JPDUMMY  = 20   ! Size of dummy array

INTEGER, PARAMETER :: JPOUTMAX = 48  ! Maximum allowed number of OUTput files

END MODULE MODD_PARAMETERS_DEP
