!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 10/01/2019: use NEWUNIT argument of OPEN (removed ISTDOUT, ISTDERR, added NNULLUNIT, CNULLFILE)
!  P. Wautelet 21/01/2019: add LIO_ALLOW_NO_BACKUP and LIO_NO_WRITE to modd_io_ll to allow to disable writes (for bench purposes)
!  P. Wautelet 07/02/2019: force TYPE to a known value for IO_File_add2list
!  P. Wautelet 12/03/2019: add TMAINFILE field in TFILEDATA
!  P. Wautelet 17/01/2020: add 'BUD' category for Print_msg + corresponding namelist variables
!  P. Wautelet 22/09/2020: add ldimreduced in tfiledata
!  P. Wautelet 10/11/2020: new data structures for netCDF dimensions
!-----------------------------------------------------------------

#define MNH_REDUCE_DIMENSIONS_IN_FILES 1

MODULE MODD_IO
!
use modd_netcdf,     only: tdimsnc
USE MODD_PARAMETERS, ONLY: NDIRNAMELGTMAX, NFILENAMELGTMAX
use modd_precision,  only: CDFINT, LFIINT
!
IMPLICIT NONE
!
!
INTEGER, PARAMETER :: NVERB_NO=0, NVERB_FATAL=1, NVERB_ERROR=2, NVERB_WARNING=3, NVERB_INFO=4, NVERB_DEBUG=5

INTEGER                     :: NNULLUNIT = -1  ! /dev/null fortran unit, value set in IO_Init
CHARACTER(LEN=*), PARAMETER :: CNULLFILE = "/dev/null"

INTEGER, SAVE :: NIO_RANK ! Rank of IO process
INTEGER, SAVE :: ISP     !! Actual proc number
INTEGER, SAVE :: ISNPROC !! Total number of allocated processes
LOGICAL, SAVE :: GSMONOPROC = .FALSE. !! True if sequential execution (ISNPROC = 1) 

LOGICAL, SAVE :: L1D   = .FALSE. ! TRUE if 1D model version
LOGICAL, SAVE :: L2D   = .FALSE. ! TRUE if 2D model version
LOGICAL, SAVE :: LPACK = .FALSE. ! TRUE if FM compression occurs in 1D or 2D model version

LOGICAL, SAVE :: LIOCDF4    = .FALSE. ! TRUE will enable full NetCDF4 (HDF5) I/O support
LOGICAL, SAVE :: LLFIOUT    = .FALSE. ! TRUE will also force LFI output when LIOCDF4 is on (debug only)  
LOGICAL, SAVE :: LLFIREAD   = .FALSE. ! TRUE will force LFI read (instead of NetCDF) when LIOCDF4 is on (debug only)  

LOGICAL, SAVE :: LVERB_OUTLST = .TRUE.  ! TRUE will PRINT_MSG in OUTPUT_LISTINGn files
LOGICAL, SAVE :: LVERB_STDOUT = .FALSE. ! TRUE will also PRINT_MSG on standard output
LOGICAL, SAVE :: LVERB_ALLPRC = .FALSE. ! FALSE: only process 0 do PRINT_MSG, TRUE: all processes
INTEGER, SAVE :: NBUD_VERB        = NVERB_INFO    ! Verbosity level for budgets
INTEGER, SAVE :: NBUD_ABORT_LEVEL = NVERB_ERROR   ! Level of budget error necessary to force stop of application
INTEGER, SAVE :: NIO_VERB        = NVERB_INFO    ! Verbosity level for IO
INTEGER, SAVE :: NIO_ABORT_LEVEL = NVERB_ERROR   ! Level of IO error necessary to force stop of application

INTEGER, SAVE :: NGEN_VERB        = NVERB_INFO    ! Verbosity level for 'GEN' (generic) messages
INTEGER, SAVE :: NGEN_ABORT_LEVEL = NVERB_ERROR   ! Level of 'GEN' error necessary to force stop of application

CHARACTER(LEN=NDIRNAMELGTMAX) :: CIO_DIR = '' ! Directory for IO

logical, save :: LIO_ALLOW_NO_BACKUP = .false. ! Allow to have no valid backup time (useful for some tests)
logical, save :: LIO_NO_WRITE        = .false. ! Disable file writes (useful for benchs)

!Structure containing one pointer to a file
!Useful to create arrays of pointers to files
TYPE TFILE_ELT
  TYPE(TFILEDATA),POINTER :: TFILE => NULL()
END TYPE TFILE_ELT

!Structure describing the characteristics of an output or a backup
TYPE TOUTBAK
  INTEGER           :: NID = -1     !Backup number
  INTEGER           :: NSTEP        !Timestep number
  REAL              :: XTIME        !Time from start of the segment (in seconds and rounded to a timestep)
  INTEGER           :: NOUTDAD = -1 !Index of the corresponding dad file (file with same time)
  TYPE(TFILEDATA),POINTER :: TFILE => NULL() !Corresponding file
  TYPE(TFILE_ELT),DIMENSION(:),ALLOCATABLE :: TFILE_IOZ !Corresponding Z-split files
  INTEGER,DIMENSION(:),POINTER :: NFIELDLIST => NULL() !List of the fields to read or write
END TYPE TOUTBAK

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
#if ( MNH_REDUCE_DIMENSIONS_IN_FILES == 1 )
  logical           :: ldimreduced   = .true.  !True if number of dimensions of fields can be reduced (for 2D simulations)
#else
  logical           :: ldimreduced   = .false. !True if number of dimensions of fields can be reduced (for 2D simulations)
#endif
  !
  INTEGER           :: NSUBFILES_IOZ = 0       !Number of sub-files (Z-split files based on this file)
                                               !For example if 2 sub-files and this file is abcd,
                                               !the 2 sub-files are abcd.Z001 and abcd.Z002
  TYPE(TFILE_ELT),DIMENSION(:),ALLOCATABLE :: TFILES_IOZ !Corresponding Z-split files
  !
  INTEGER              :: NMODEL = 0              !Model number corresponding to the file (field not always set)
  INTEGER,DIMENSION(3) :: NMNHVERSION = (/0,0,0/) !MesoNH version used to create the file
  !
#ifdef MNH_IOLFI
  ! Fields for LFI files
  INTEGER(KIND=LFIINT) :: NLFININAR = 0  !Number of articles of the LFI file (only accurate if file opened in read mode)
  INTEGER(KIND=LFIINT) :: NLFINPRAR = 0  !Number of predicted articles of the LFI file (non crucial)
  INTEGER              :: NLFITYPE  = -1 !Type of the file (used to generate list of files to transfers)
  INTEGER              :: NLFIVERB  = 1  !LFI verbosity level
  INTEGER(KIND=LFIINT) :: NLFIFLU   = -1 !File identifier
#endif
  !
#ifdef MNH_IOCDF4
  ! Fields for netCDF files
  INTEGER(KIND=CDFINT)   :: NNCID = -1 !File identifier (corresponding to the actual group)
  INTEGER(KIND=CDFINT)   :: NNCNAR = 0 !Number of articles of the netCDF file (only accurate if file opened in read mode)
  LOGICAL                :: LNCREDUCE_FLOAT_PRECISION = .FALSE. ! Reduce the precision of floats to single precision
                                                                ! instead of double precision
  LOGICAL                :: LNCCOMPRESS = .FALSE. ! Do compression on fields
  INTEGER(KIND=CDFINT)   :: NNCCOMPRESS_LEVEL = 0 ! Compression level
  type(tdimsnc), pointer :: tncdims => Null()     ! Dimensions of netCDF file
#endif
  !
  !Fields for other files
  INTEGER :: NLU = -1                      !Logical unit number
  INTEGER :: NRECL = -1                    !Fortran RECL (record length)
  CHARACTER(LEN=11) :: CFORM   = "UNKNOWN" !Fortran FORM (FORMATTED/UNFORMATTED)
  CHARACTER(LEN=10) :: CACCESS = "UNKNOWN" !Fortran ACCESS (DIRECT/SEQUENTIAL/STREAM)
  !
  TYPE(TFILEDATA),POINTER :: TDADFILE   => NULL() !Corresponding dad file
  TYPE(TFILEDATA),POINTER :: TDESFILE   => NULL() !Corresponding .des file
  TYPE(TFILEDATA),POINTER :: TDATAFILE  => NULL() !Corresponding data file (if .des file)
  TYPE(TFILEDATA),POINTER :: TMAINFILE  => NULL() !Corresponding main file if the file is an sub-file
  !
  TYPE(TFILEDATA),POINTER :: TFILE_PREV => NULL()
  TYPE(TFILEDATA),POINTER :: TFILE_NEXT => NULL()
END TYPE TFILEDATA

!Structure containing a pointer to a file (useful to create arrays of pointers to files)
TYPE TPTR2FILE
  TYPE(TFILEDATA),POINTER :: TZFILE => NULL()
END TYPE

TYPE(TFILEDATA),POINTER,SAVE :: TFILE_FIRST => NULL()
TYPE(TFILEDATA),POINTER,SAVE :: TFILE_LAST  => NULL()

TYPE(TFILEDATA),POINTER,SAVE :: TFILE_SURFEX  => NULL() !Pointer used to find the file used when writing SURFEX fields in write_surf_mnh.f90

TYPE(TFILEDATA),POINTER,SAVE :: TFILE_OUTPUTLISTING  => NULL() !Pointer used to point to the file used when writing to OUTPUT_LISTINGn file

!Non existing file which can be used as a dummy target
TYPE(TFILEDATA),TARGET, SAVE :: TFILE_DUMMY = TFILEDATA(CNAME="dummy",CDIRNAME=NULL(),TFILES_IOZ=NULL())

END MODULE MODD_IO
