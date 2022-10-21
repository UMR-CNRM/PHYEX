!     ######spl
      MODULE MODD_LUNIT
!     #################
!
!!****  *MODD_LUNIT* - declaration of names and logical unit numbers of files
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the
!     logical unit numbers  of  output file for all models.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_LUNIT)
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/05/94
!!      V. Masson   01/2004 add file names for use in externalized surface!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
INTEGER :: ILUOUT
CHARACTER(LEN=16),SAVE :: CLUOUT0    ! Name of output_listing file
CHARACTER(LEN=28),SAVE :: COUTFMFILE ! name of the output FM-file being written
CHARACTER(LEN=28),SAVE :: CPGDFILE   ! name of the PGD file for PREP_REAL_CASE
!
END MODULE MODD_LUNIT
