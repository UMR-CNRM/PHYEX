!MNH_LIC Copyright 1996-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      PROGRAM XY_TO_LATLON
!     ####################
!
!!**** *XY_TO_LATLON* program to compute latitude and longiude from x and y
!!                    for a MESONH file
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!    module MODE_GRIDPROJ : contains projection routines
!!                         SM_LATLON and SM_XYHAT
!!                           
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    module MODD_GRID     : variables for projection:
!!                         XLAT0,XLON0,XRPK,XBETA
!!
!!    module MODD_PGDDIM  : specify the dimentions of the data arrays:
!!                         NPGDIMAX and NPGDJMAX
!!
!!    module MODD_PGDGRID : grid variables:
!!                         XPGDLONOR,XPGDLATOR: longitude and latitude of the
!!                           origine point for the conformal projection.
!!                         XPGDXHAT,XPGDYHAT: position x,y in the conformal plane   
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson          Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     26/01/96
!!
!!    no transfer of the file when closing   Dec. 09, 1996 (V.Masson)
!!    + changes call to READ_HGRID
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 07/02/2019: force TYPE to a known value for IO_File_add2list
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 14/04/2020: add missing initializations (XY_TO_LATLON was not working)
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
use MODD_CONF,             only: CPROGRAM
USE MODD_DIM_n
USE MODD_GRID
USE MODD_IO,               ONLY: TFILEDATA
USE MODD_PGDDIM
USE MODD_PGDGRID
USE MODD_PARAMETERS
USE MODD_LUNIT
!
USE MODE_FIELD,            ONLY: INI_FIELD_LIST
USE MODE_GRIDPROJ
USE MODE_INIT_ll,          only: SET_DIM_ll, SET_JP_ll
USE MODE_IO,               only: IO_Config_set, IO_Init
use MODE_IO_FIELD_READ,    only: IO_Field_read
USE MODE_IO_FILE,          only: IO_File_close, IO_File_open
USE MODE_IO_MANAGE_STRUCT, only: IO_File_add2list
USE MODE_MODELN_HANDLER,   ONLY: GOTO_MODEL
use MODE_SPLITTINGZ_ll
!
USE MODE_INI_CST, ONLY: INI_CST
USE MODI_READ_HGRID
!
USE MODN_CONFIO,           ONLY: NAM_CONFIO
!
IMPLICIT NONE
!
!*    0.2    Declaration of variables
!            ------------------------
!
CHARACTER(LEN=28) :: YINIFILE  ! name of input FM file
CHARACTER(LEN=28) :: YNAME     ! true name of input FM file
CHARACTER(LEN=28) :: YDAD      ! name of dad of input FM file
CHARACTER(LEN=2)  :: YSTORAGE_TYPE
INTEGER :: INAM                ! Logical unit for namelist file
INTEGER :: ILUOUT0             ! Logical unit for output file.
INTEGER :: IRESP               ! Return-code if problem eraised.
REAL    :: ZI,ZJ               ! input positions of the point
INTEGER :: II,IJ               ! integer positions of the point
REAL    :: ZXHAT               ! output conformal coodinate x
REAL    :: ZYHAT               ! output conformal coodinate y
REAL    :: ZLAT                ! output latitude
REAL    :: ZLON                ! output longitude
TYPE(TFILEDATA),POINTER :: TZINIFILE => NULL()
TYPE(TFILEDATA),POINTER :: TZNMLFILE => NULL()
!
!*    0.3    Declaration of namelists
!            ------------------------
!
NAMELIST/NAM_INIFILE/ YINIFILE
!----------------------------------------------------------------------------
!
    WRITE(*,*) '+---------------------------------+'
    WRITE(*,*) '|      program xy_to_latlon       |'
    WRITE(*,*) '+---------------------------------+'
    WRITE(*,*) ''
    WRITE(*,*) 'Warning: I and J are integer for flux points'
!
!*    1.     Initializations
!            ---------------
!
CALL GOTO_MODEL(1)
!
CALL VERSION()
!
CPROGRAM='LAT2XY'
!
CALL IO_Init()
!
CALL INI_CST()
!
CALL INI_FIELD_LIST()
!
!*    2.     Reading of namelist file
!            ------------------------
!
CALL IO_File_add2list(TZNMLFILE,'XY2LATLON1.nam','NML','READ')
CALL IO_File_open(TZNMLFILE)
INAM=TZNMLFILE%NLU
READ(INAM,NAM_INIFILE)
!
READ(INAM,NAM_CONFIO)
CALL IO_Config_set()
CALL IO_File_close(TZNMLFILE)
!
!*    1.     Opening of MESONH file
!            ----------------------
!
CALL IO_File_add2list(TZINIFILE,TRIM(YINIFILE),'MNH','READ',KLFITYPE=2,KLFIVERB=2)
CALL IO_File_open(TZINIFILE)
!
CALL IO_Field_read(TZINIFILE,'IMAX',  NIMAX)
CALL IO_Field_read(TZINIFILE,'JMAX',  NJMAX)
NKMAX = 1
CALL IO_Field_read(TZINIFILE,'JPHEXT',JPHEXT)
!
CALL SET_JP_ll(1,JPHEXT,JPVEXT,JPHEXT)
CALL SET_DIM_ll(NIMAX, NJMAX, NKMAX)
CALL INI_PARAZ_ll(IRESP)
!
!*    2.     Reading of MESONH file
!            ----------------------
!
CALL READ_HGRID(0,TZINIFILE,YNAME,YDAD,YSTORAGE_TYPE)
!
!*    3.     Closing of MESONH file
!            ----------------------
!
CALL IO_File_close(TZINIFILE)
!
!-------------------------------------------------------------------------------
!
!*    4.     Reading of I and J
!            ------------------
!
DO
  WRITE(*,*) '-------------------------------------------------------------------'
  WRITE(*,*) 'please enter index I  (real, quit or q to stop):'
  READ(*,*,ERR=1) ZI
  WRITE(*,*) 'please enter index J  (real, quit or q to stop):'
  READ(*,*,ERR=1) ZJ
!
  II=MAX(MIN(INT(ZI),NPGDIMAX+2*JPHEXT-1),1)
  IJ=MAX(MIN(INT(ZJ),NPGDJMAX+2*JPHEXT-1),1)
  ZXHAT=XPGDXHAT(II) + (ZI-REAL(II)) * ( XPGDXHAT(II+1) - XPGDXHAT(II) )
  ZYHAT=XPGDYHAT(IJ) + (ZJ-REAL(IJ)) * ( XPGDYHAT(IJ+1) - XPGDYHAT(IJ) )
!
  WRITE(*,*) 'x=', ZXHAT
  WRITE(*,*) 'y=', ZYHAT
!
  CALL SM_LATLON(XPGDLATOR,XPGDLONOR,   &
                 ZXHAT,ZYHAT,ZLAT,ZLON)
!
  WRITE(*,*) 'lat=', ZLAT
  WRITE(*,*) 'lon=', ZLON
END DO
1   WRITE(*,*) 'good bye'
!
!-------------------------------------------------------------------------------
!
END PROGRAM XY_TO_LATLON
