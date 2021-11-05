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
!!	V. Ducrocq   *Meteo France*
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
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!      B.VIE 2016 LIMA
! P. Wautelet: 05/2016-04/2018: new data structures and calls for I/O
! Q. Rodier   29/03/2019: increase maximum number of outputs to 999
! P. Wautelet 17/01/2020: add NBUNAMELGTMAX and NCOMMENTLGTMAX parameters
! P. Wautelet 13/03/2020: remove JPBUMAX and JPBUPROMAX
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!JUAN CYCLK
!INTEGER, PARAMETER :: JPHEXT = 3     ! Horizontal External points number
INTEGER,SAVE      :: JPHEXT = 1     ! Horizontal External points number
!
!JUAN CYCLK
INTEGER, PARAMETER :: JPVEXT = 1      ! Vertical External points number
INTEGER, PARAMETER :: JPVEXT_TURB = 1      ! Vertical External points number
INTEGER, PARAMETER :: JPMODELMAX = 8  ! Maximum allowed number of nested models 
INTEGER, PARAMETER :: JPCPLFILEMAX = 24 ! Maximum allowed number of CouPLing FILEs
INTEGER, PARAMETER :: JPRIMMAX = 6    ! Maximum number of points for the
                       ! horizontal relaxation for the outermost verticals
INTEGER, PARAMETER :: JPSVMAX  = 200  ! Maximum number of scalar variables
INTEGER, PARAMETER :: JPSVNAMELGTMAX = 10 ! Maximum length of a scalar variable name (do not set to less than 10)
!
!
REAL,    PARAMETER :: XUNDEF = 999.     ! default value for undefined or unused
!                                       ! field.
REAL,    PARAMETER :: XNEGUNDEF = -999. ! default value for undefined or unused
!                                       ! field (negative value guaranteed)
INTEGER, PARAMETER :: NUNDEF = 999      ! default value for undefined or unused
!                                       ! field.
INTEGER, PARAMETER :: NNEGUNDEF = -999  ! default value for undefined or unused
!                                       ! field (negative value guaranteed)
INTEGER, PARAMETER :: JPDUMMY  = 20   ! Size of dummy array
!
INTEGER, PARAMETER :: JPOUTMAX = 999    ! Maximum allowed number of OUTput files
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
!
INTEGER, PARAMETER :: NGRIDUNKNOWN = -1 ! Unknown Arakawa grid number
!
END MODULE MODD_PARAMETERS
