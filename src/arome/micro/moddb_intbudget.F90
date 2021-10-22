!     ######spl
      MODULE MODDB_INTBUDGET
!     ##################
!
!!****  *MODDB_INTBUDGET* - module for interfacing MNH's budgets with DDH
!!
!!    PURPOSE
!!    -------
!!       Passing some arrays from apl_arome.f90 to mpa subroutines: RHODJ,QDM and EXNREFM
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      O. Riviere
!!
!!    MODIFICATIONS
!!    -------------
!!      21/06/08
!!      18/09/17  F.Voitus 
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
IMPLICIT NONE

TYPE PTR_VAR
REAL,DIMENSION(:,:),POINTER:: VARMULT
END TYPE PTR_VAR

INTEGER,SAVE::NBUDGET
INTEGER,SAVE::NLON
INTEGER,SAVE::NLEV
REAL,DIMENSION(:,:),SAVE,ALLOCATABLE,TARGET:: TCON2
REAL,DIMENSION(:,:),SAVE,ALLOCATABLE,TARGET:: TCON3
REAL,DIMENSION(:,:),SAVE,ALLOCATABLE,TARGET:: TCON1
CHARACTER(LEN=2),DIMENSION(13):: CVARNAME=       &
(/"UU","VV","WW","CT","KK","QV","QL","QR","QI","QS","QG","QH","SV"/)
TYPE(PTR_VAR),DIMENSION(13)::TAB_VARMULT
REAL,DIMENSION(:,:,:,:),SAVE,ALLOCATABLE:: TVARSM

END MODULE MODDB_INTBUDGET
