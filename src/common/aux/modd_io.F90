MODULE MODD_IO
USE MODD_PARAMETERS, ONLY: NFILENAMELGTMAX
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: NVERB_NO=0, NVERB_FATAL=1, NVERB_ERROR=2, NVERB_WARNING=3, NVERB_INFO=4, NVERB_DEBUG=5
INTEGER, SAVE :: N_ABORT_LEVEL = NVERB_ERROR
!
!Structure describing the characteristics of a file
TYPE TFILEDATA
  CHARACTER(LEN=NFILENAMELGTMAX) :: CNAME = '' !Filename
  CHARACTER(LEN=:),ALLOCATABLE   :: CDIRNAME   !Directory name
  CHARACTER(LEN=13) :: CTYPE   = "UNKNOWN" !Filetype (PGD, MNH, DES, NML...)
  CHARACTER(LEN=7)  :: CFORMAT = "UNKNOWN" !Fileformat (NETCDF4, LFI, LFICDF4...)
  CHARACTER(LEN=7)  :: CMODE   = "UNKNOWN" !Opening mode (read, write...)
  LOGICAL           :: LOPENED = .FALSE.   !Is the file opened
  INTEGER           :: NOPEN_CURRENT = 0   !Number of times the file is currently opened (several opens without close are allowed)
  INTEGER           :: NOPEN   = 0         !Number of times the file has been opened (during the current execution)
  INTEGER           :: NCLOSE  = 0         !Number of times the file has been closed (during the current execution)
  !
  INTEGER           :: NMASTER_RANK  = -1      !Rank of the master process (no meaning if LMULTIMASTERS=.T.)
  INTEGER           :: NMPICOMM      = -1      !MPI communicator used for IO on this file
  LOGICAL           :: LMASTER       = .FALSE. !True if process is master of the file (process that open/read/write/close)
  LOGICAL           :: LMULTIMASTERS = .FALSE. !True if several processes may access the file
  INTEGER           :: NSUBFILES_IOZ = 0       !Number of sub-files (Z-split files based on this file)
                                               !For example if 2 sub-files and this file is abcd,
                                               !the 2 sub-files are abcd.Z001 and abcd.Z002
!  TYPE(TFILE_ELT),DIMENSION(:),ALLOCATABLE :: TFILES_IOZ !Corresponding Z-split files
  !
  INTEGER              :: NMODEL = 0              !Model number corresponding to the file (field not always set)
  INTEGER,DIMENSION(3) :: NMNHVERSION = (/0,0,0/) !MesoNH version used to create the file
  !
  INTEGER :: NLU = -1 ! logical unit number
END TYPE TFILEDATA
ENDMODULE MODD_IO
