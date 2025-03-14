!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
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
INTEGER, PARAMETER :: JPSVNAMELGTMAX = 10 ! Maximum length of a scalar variable name (do not set to less than 10)
!
!
REAL,    PARAMETER :: XUNDEF = 1.E+20   ! default value for undefined or unused
!                                       ! field.
REAL,    PARAMETER :: XNEGUNDEF = -999. ! default value for undefined or unused
!                                       ! field (negative value guaranteed)
INTEGER, PARAMETER :: NUNDEF = 1E+9     ! default value for undefined or unused
!                                       ! field.
INTEGER, PARAMETER :: NNEGUNDEF = -999  ! default value for undefined or unused
!                                       ! field (negative value guaranteed)
INTEGER, PARAMETER :: JPDUMMY  = 20   ! Size of dummy array
!
INTEGER, PARAMETER :: JPOUTMAX = 192 ! Maximum allowed number of OUTput files
INTEGER, PARAMETER :: JPOUTVARMAX = 192 ! Maximum allowed number of variables in an output file
!
INTEGER, PARAMETER :: NBUNAMELGTMAX  = 32  ! Maximum length of a budget name
INTEGER, PARAMETER :: NCOMMENTLGTMAX = 100 ! Maximum length of a comment
INTEGER, PARAMETER :: NMNHNAMELGTMAX = 32  ! Maximum length of a MNH variable name
INTEGER, PARAMETER :: NSTDNAMELGTMAX = 64  ! Maximum length of the standard name of a variable (CF convention)
!
INTEGER, PARAMETER :: NDIRNAMELGTMAX = 512 ! Maximum length of a directory name
INTEGER, PARAMETER :: NFILENAMELGTMAX = 32 ! Maximum length of a file name (must be at least NFILENAMELGTMAXLFI)
INTEGER, PARAMETER :: NFILENAMELGTMAXLFI = 28 ! Maximum length of a file name in LFI file (this is necessary
                                              ! to keep backward compatibility), MUST BE 28
!
INTEGER, PARAMETER :: NLFIMAXCOMMENTLENGTH = 100 ! Length of comments in LFI files
!
INTEGER, PARAMETER :: JPLIMACCNMAX = 10 ! Maximum allowed number of CCN modes in LIMA
INTEGER, PARAMETER :: JPLIMAIFNMAX = 10 ! Maximum allowed number of IFN modes in LIMA
INTEGER, PARAMETER :: NNBCRYSTALMAX = 4 ! Maximum allowed number of IFN modes in LIMA
!
INTEGER, PARAMETER :: NGRIDUNKNOWN = -1 ! Unknown Arakawa grid number
INTEGER, PARAMETER :: NEXPNAMELGTMAX    = 32 ! should be at least 5
INTEGER, PARAMETER :: NSEGNAMELGTMAX    = 32 ! should be at least 5
!
END MODULE MODD_PARAMETERS
