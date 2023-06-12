!MNH_LIC Copyright 2004-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################
      PROGRAM PREP_SURFEX
!     #############################
!
!!****  *PREP_SURFEX* - program to write an initial FM file from real case
!!                                situation containing only surface fields.
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/2004 (P. Le Moigne)
!!    10/10/2011  J.Escobar call INI_PARAZ_ll
!!  06/2016     (G.Delautier) phasage surfex 8
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 07/02/2019: force TYPE to a known value for IO_File_add2list
!!          2021  B.Vie     LIMA - CAMS coupling
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF,        ONLY : CPROGRAM,&
                             L1D, L2D, LPACK
USE MODD_CONF_n,      ONLY : CSTORAGE_TYPE
USE MODD_IO,          ONLY : TFILEDATA, TFILE_SURFEX
USE MODD_LUNIT,       ONLY : TPGDFILE, TLUOUT0
USE MODD_LUNIT_n,     ONLY : CINIFILE, TINIFILE
USE MODD_MNH_SURFEX_n
USE MODD_PARAMETERS,  ONLY : JPMODELMAX,JPHEXT,JPVEXT, NUNDEF, XUNDEF
USE MODD_TIME_n,      ONLY : TDTCUR
!
use mode_field,            only: Ini_field_list, Ini_field_scalars
USE MODE_FINALIZE_MNH,     only: FINALIZE_MNH
USE MODE_IO,               only: IO_Init
USE MODE_IO_FIELD_READ,    only: IO_Field_read
USE MODE_IO_FIELD_WRITE,   only: IO_Field_write, IO_Header_write
USE MODE_IO_FILE,          only: IO_File_close, IO_File_open
USE MODE_IO_MANAGE_STRUCT, only: IO_File_add2list
USE MODE_ll
USE MODE_MSG
USE MODE_MODELN_HANDLER
USE MODE_SPLITTINGZ_ll
!
USE MODI_OPEN_PRC_FILES
USE MODI_PREP_SURF_MNH
USE MODI_READ_ALL_NAMELISTS
USE MODI_VERSION
USE MODE_INI_CST, ONLY: INI_CST
!
IMPLICIT NONE
!
!*       0.1   Declaration of local variables
!              ------------------------------
!
CHARACTER(LEN=28)     :: YATMFILE        ! name of the Atmospheric file
CHARACTER(LEN=6)      :: YATMFILETYPE    ! type of the Atmospheric file
CHARACTER(LEN=28)     :: YCHEMFILE       ! name of the Chemical file (not used)
CHARACTER(LEN=6)      :: YCHEMFILETYPE   ! type of the Chemical file (not used)
CHARACTER(LEN=28)     :: YCAMSFILE       ! name of the input CAMS file
CHARACTER(LEN=6)      :: YCAMSFILETYPE   ! type of the input CAMS file
CHARACTER(LEN=28)     :: YSURFFILE       ! name of the Surface file (not used)
CHARACTER(LEN=6)      :: YSURFFILETYPE   ! type of the Surface file (not used)
CHARACTER(LEN=28)     :: YPGDFILE        ! name of the physiographic data
!                                        ! file
!
!* file management variables and counters
!
INTEGER               :: ILUOUT0         ! logical unit for listing file
INTEGER               :: IRESP           ! return code in FM routines
!
INTEGER               :: IINFO_ll        ! return code of // routines
CHARACTER (LEN=100)   :: HCOMMENT
INTEGER               :: II, IJ, IGRID, ILENGTH
!
TYPE(TFILEDATA),POINTER :: TZATMFILE => NULL()
TYPE(TFILEDATA),POINTER :: TZPRE_REAL1FILE => NULL()
!
!-------------------------------------------------------------------------------
!
!
!*       1.    SET DEFAULT VALUES
!              ------------------
!
CALL GOTO_MODEL(1)
!
CALL VERSION
CPROGRAM='REAL  '
CSTORAGE_TYPE='SU'
!
!-------------------------------------------------------------------------------
!
!*       2.    OPENNING OF THE FILES
!              ---------------------
CALL IO_Init()
!
CALL OPEN_PRC_FILES(TZPRE_REAL1FILE,YATMFILE, YATMFILETYPE,TZATMFILE &
                                   ,YCHEMFILE,YCHEMFILETYPE &
                                   ,YSURFFILE,YSURFFILETYPE &
                                   ,YPGDFILE,TPGDFILE       &
                                   ,YCAMSFILE,YCAMSFILETYPE)
ILUOUT0 = TLUOUT0%NLU
!
!-------------------------------------------------------------------------------
!
!*       3.    INITIALIZATION OF PHYSICAL CONSTANTS
!              ------------------------------------
!
CALL INI_CST
!
!-------------------------------------------------------------------------------
!
!*       4.    READING OF NAMELIST
!              -------------------
!
!*       4.1   reading of configuration variables
!
CALL IO_File_close(TZPRE_REAL1FILE)
!
!*       4.2   reading of values of some configuration variables in namelist
!
CALL INI_FIELD_LIST()
!
CALL INI_FIELD_SCALARS()
!
CALL IO_Field_read(TPGDFILE,'IMAX',II)
CALL IO_Field_read(TPGDFILE,'JMAX',IJ)
CALL SET_JP_ll(JPMODELMAX,JPHEXT,JPVEXT,JPHEXT)
CALL SET_DAD0_ll()
CALL SET_DIM_ll(II, IJ, 1)
CALL SET_LBX_ll('OPEN',1)
CALL SET_LBY_ll('OPEN', 1)
CALL SET_XRATIO_ll(1, 1)
CALL SET_YRATIO_ll(1, 1)
CALL SET_XOR_ll(1, 1)
CALL SET_XEND_ll(II+2*JPHEXT, 1)
CALL SET_YOR_ll(1, 1)
CALL SET_YEND_ll(IJ+2*JPHEXT, 1)
CALL SET_DAD_ll(0, 1)
!JUANZ CALL INI_PARA_ll(IINFO_ll)
CALL INI_PARAZ_ll(IINFO_ll)
!
!-------------------------------------------------------------------------------
!
!
!*       5.    PREPARATION OF SURFACE FIELDS
!              -----------------------------
!
!* reading of date
!
IF (YATMFILETYPE=='MESONH') THEN
  CALL IO_File_add2list(TZATMFILE,TRIM(YATMFILE),'MNH','READ',KLFITYPE=1,KLFIVERB=1)
  CALL IO_File_open(TZATMFILE)
  CALL IO_Field_read(TZATMFILE,'DTCUR',TDTCUR)
  CALL IO_File_close(TZATMFILE)
ELSE
  TDTCUR%nyear  = NUNDEF
  TDTCUR%nmonth = NUNDEF
  TDTCUR%nday   = NUNDEF
  TDTCUR%xtime  = XUNDEF
END IF
!
CALL SURFEX_ALLOC_LIST(1)
YSURF_CUR => YSURF_LIST(1)
CALL READ_ALL_NAMELISTS(YSURF_CUR,'MESONH','PRE',.FALSE.)
CALL GOTO_SURFEX(1)
!
CALL IO_File_add2list(TINIFILE,TRIM(CINIFILE),'PGD','WRITE',KLFITYPE=1,KLFIVERB=1)
!The open is done later in PREP_SURF_MNH when domain dimensions are known
!
TFILE_SURFEX => TINIFILE
CALL PREP_SURF_MNH(YATMFILE,YATMFILETYPE,OINIFILEOPEN=.TRUE.)
NULLIFY(TFILE_SURFEX)
!
!-------------------------------------------------------------------------------
!
CALL IO_Header_write(TINIFILE)
CALL IO_Field_write(TINIFILE,'SURF','EXTE')
CALL IO_Field_write(TINIFILE,'L1D', L1D)
CALL IO_Field_write(TINIFILE,'L2D', L2D)
CALL IO_Field_write(TINIFILE,'PACK',LPACK)
!
!-------------------------------------------------------------------------------
WRITE(ILUOUT0,*) ' '
WRITE(ILUOUT0,*) '----------------------------------'
WRITE(ILUOUT0,*) '|                                |'
WRITE(ILUOUT0,*) '|   PREP_SURFEX ends correctly   |'
WRITE(ILUOUT0,*) '|                                |'
WRITE(ILUOUT0,*) '----------------------------------'
CALL IO_File_close(TINIFILE)
!
CALL FINALIZE_MNH()
!-------------------------------------------------------------------------------
!
END PROGRAM PREP_SURFEX
