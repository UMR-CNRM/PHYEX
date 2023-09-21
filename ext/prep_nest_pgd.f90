!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #####################
      PROGRAM PREP_NEST_PGD
!     #####################
!
!!****  *PREP_NEST_PGD* - to make coherent pgd files for nesting
!!
!!    PURPOSE
!!    -------
!!
!!       The purpose of this program is to prepare pgd files with which
!!       nesting can be performed. A pgd file must be coherent with its
!!       father:
!!         The average of orography of fine model on each of its father grid
!!       mesh must be the same as its father orography.
!!
!!       All the pgd files are read at the begining of the program,
!!       then they are checked, and recursively, the orography of a father
!!       is replaced by the averaged orography from ist son.
!!
!!       The control data are given in the namelist file PRE_NEST.nam
!!
!! &NAM_NEST_PGD1 CPGD='coarser model' /
!! &NAM_NEST_PGD2 CPGD='medium model' , IDAD=1 /
!! &NAM_NEST_PGD3 CPGD='medium model' , IDAD=1 /
!! &NAM_NEST_PGD4 CPGD='fine model' , IDAD=2 /
!! &NAM_NEST_PGD5 CPGD='fine model' , IDAD=2 /
!! &NAM_NEST_PGD6 CPGD='fine model' , IDAD=3 /
!! &NAM_NEST_PGD7 CPGD='very fine model' , IDAD=6 /
!! &NAM_NEST_PGD8 CPGD='very very fine model' , IDAD=7 /
!!
!!        In each namelist is given the name of the pgd file, and the number
!!      of its father. This one MUST be smaller.
!!        There is one output file for each input file, with the suffix
!!      '.nest' added at the end of the file name (even if the file has not
!!      been changed).
!!
!!        In the case of the namelist above, one obtain something like:
!!
!!   +----------------------------------------------------------+
!!   |                                                 1        |
!!   |   +-----------------------+                              |
!!   |   |                    2  |                              |
!!   |   |                       |                              |
!!   |   |              +-+      |                              |
!!   |   | +-------+    |5|      |   +-----------------------+  |
!!   |   | |  4    |    +-+      |   |   +----------+     3  |  |
!!   |   | +-------+             |   |   |+------+ 6|        |  |
!!   |   +-----------------------+   |   || +-+ 7|  |        |  |
!!   |                               |   || |8|  |  |        |  |
!!   |                               |   || +-+  |  |        |  |
!!   |                               |   |+------+  |        |  |
!!   |                               |   +----------+        |  |
!!   |                               +-----------------------+  |
!!   +----------------------------------------------------------+
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
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
!!      Original    26/09/95
!!                  30/07/97 (Masson) split of mode_lfifm_pgd
!!                  2014 (M.Faivre)
!!                  06/2015 (M.Moge) parallelization 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      J.Escobar : 19/04/2016 : Pb IOZ/NETCDF , missing OPARALLELIO=.FALSE. for PGD files
!!  06/2016     (G.Delautier) phasage surfex 8
!!      P.Wautelet : 08/07/2016 : removed MNH_NCWRIT define
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 07/02/2019: force TYPE to a known value for IO_File_add2list
!  P. Wautelet 06/07/2021: use FINALIZE_MNH
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST
USE MODD_DIM_n
USE MODD_IO,               ONLY: TFILE_SURFEX, TPTR2FILE
USE MODD_GRID_n,           ONLY: XZSMT
USE MODD_LUNIT,            ONLY: TPGDFILE,TLUOUT0,TOUTDATAFILE
USE MODD_MNH_SURFEX_n
USE MODD_NESTING
USE MODD_PARAMETERS
USE MODD_VAR_ll,           ONLY: NPROC, IP, NMNH_COMM_WORLD
!
use mode_field,            only: Ini_field_list
USE MODE_FINALIZE_MNH,     only: FINALIZE_MNH
USE MODE_IO,               only: IO_Init, IO_Pack_set
USE MODE_IO_FIELD_READ,    only: IO_Field_read
USE MODE_IO_FIELD_WRITE,   only: IO_Field_write, IO_Header_write
USE MODE_IO_FILE,          only: IO_File_close, IO_File_open
USE MODE_IO_MANAGE_STRUCT, only: IO_File_add2list
USE MODE_ll
USE MODE_MNH_WORLD,        ONLY: INIT_NMNH_COMM_WORLD
USE MODE_MODELN_HANDLER
USE MODE_MPPDB
USE MODE_SPLITTINGZ_ll,    ONLY: INI_PARAZ_ll
!
USE MODI_DEFINE_MASK_n
USE MODI_INIT_HORGRID_ll_n
USE MODI_INIT_PGD_SURF_ATM
USE MODI_NEST_FIELD_n
USE MODI_NEST_ZSMT_n
USE MODI_OPEN_NESTPGD_FILES
USE MODI_READ_ALL_NAMELISTS
USE MODI_READ_HGRID
USE MODI_RETRIEVE1_NEST_INFO_n
USE MODI_VERSION
USE MODI_WRITE_PGD_SURF_ATM_N
USE MODE_INI_CST, ONLY: INI_CST
!
IMPLICIT NONE
!
!*       0.1   Declaration of local variables
!              ------------------------------
!
INTEGER, DIMENSION(JPMODELMAX) :: NXSIZE   ! number of grid points for each model
INTEGER, DIMENSION(JPMODELMAX) :: NYSIZE   ! in x and y-directions
                                           ! relatively to its father grid
!
INTEGER                        :: ILUOUT0
INTEGER                        :: IINFO_ll ! return code of // routines
INTEGER                        :: JPGD     ! loop control
CHARACTER(LEN=28)              :: YMY_NAME,YDAD_NAME
CHARACTER(LEN=2)               :: YSTORAGE_TYPE
LOGICAL, DIMENSION(JPMODELMAX) :: L1D_ALL  ! Flag for      1D conf. for each PGD
LOGICAL, DIMENSION(JPMODELMAX) :: L2D_ALL  ! Flag for      2D conf. for each PGD
LOGICAL, DIMENSION(JPMODELMAX) :: LPACK_ALL! Flag for packing conf. for each PGD
!
INTEGER                        :: JTIME,ITIME
INTEGER                        :: IIMAX,IJMAX,IKMAX
INTEGER                        :: IDXRATIO,IDYRATIO
INTEGER                        :: IDAD
INTEGER                        :: II
LOGICAL     :: GISINIT
!
TYPE(TPTR2FILE),DIMENSION(:),ALLOCATABLE        :: TZFILEPGD     ! Input  PGD files
TYPE(TPTR2FILE),DIMENSION(:),ALLOCATABLE,TARGET :: TZFILENESTPGD ! Output PGD files
!
!-------------------------------------------------------------------------------
!
CALL MPPDB_INIT()
!
CALL VERSION
CPROGRAM='NESPGD'
!
CALL IO_Init()
!!$CALL SET_JP_ll(JPMODELMAX,JPHEXT,JPVEXT,JPHEXT)
!
!*       1.    INITIALIZATION OF PHYSICAL CONSTANTS
!              ------------------------------------
!
CALL INI_CST
!
!-------------------------------------------------------------------------------
!
!*       2.    OPENING OF THE FILES
!              ---------------------
!
NVERB=1
!
CALL OPEN_NESTPGD_FILES(TZFILEPGD,TZFILENESTPGD)
CALL SET_JP_ll(JPMODELMAX,JPHEXT,JPVEXT,JPHEXT)
!
ILUOUT0 = TLUOUT0%NLU
!
CALL SURFEX_ALLOC_LIST(NMODEL)
YSURF_CUR => YSURF_LIST(1)
CALL READ_ALL_NAMELISTS(YSURF_CUR,'MESONH','PRE',.FALSE.)
!
!-------------------------------------------------------------------------------
!
!*       3.    READING OF THE GRIDS
!              --------------------
!
CALL INI_FIELD_LIST()
!
CALL SET_DAD0_ll()
DO JPGD=1,NMODEL
  ! read and set dimensions and ratios of model JPGD
  CALL IO_Field_read(TZFILEPGD(JPGD)%TZFILE,'IMAX',   IIMAX)
  CALL IO_Field_read(TZFILEPGD(JPGD)%TZFILE,'JMAX',   IJMAX)
  CALL IO_Field_read(TZFILEPGD(JPGD)%TZFILE,'DXRATIO',NDXRATIO_ALL(JPGD))
  CALL IO_Field_read(TZFILEPGD(JPGD)%TZFILE,'DYRATIO',NDYRATIO_ALL(JPGD))
  CALL IO_Field_read(TZFILEPGD(JPGD)%TZFILE,'XSIZE',  NXSIZE(JPGD))
  CALL IO_Field_read(TZFILEPGD(JPGD)%TZFILE,'YSIZE',  NYSIZE(JPGD))
  CALL IO_Field_read(TZFILEPGD(JPGD)%TZFILE,'XOR',    NXOR_ALL(JPGD))
  CALL IO_Field_read(TZFILEPGD(JPGD)%TZFILE,'YOR',    NYOR_ALL(JPGD))
  CALL SET_DIM_ll(IIMAX, IJMAX, 1)
  ! compute origin and end of local subdomain of model JPGD
  ! initialize variables from MODD_NESTING, origin and end of global model JPGD in coordinates of its father
  IF ( NDAD(JPGD) > 0 ) THEN
    NXEND_ALL(JPGD) = NXOR_ALL(JPGD) + NXSIZE(JPGD) - 1 + 2*JPHEXT
    NYEND_ALL(JPGD) = NYOR_ALL(JPGD) + NYSIZE(JPGD) - 1 + 2*JPHEXT
  ELSE  ! this is not a son model
    NXOR_ALL(JPGD) = 1
    NXEND_ALL(JPGD) = IIMAX+2*JPHEXT
    NYOR_ALL(JPGD) = 1
    NYEND_ALL(JPGD) = IJMAX+2*JPHEXT
    NDXRATIO_ALL(JPGD) = 1
    NDYRATIO_ALL(JPGD) = 1
  ENDIF
  ! initialize variables from MODD_DIM_ll, origin and end of global model JPGD in coordinates of its father
  CALL SET_XOR_ll(NXOR_ALL(JPGD), JPGD)
  CALL SET_XEND_ll(NXEND_ALL(JPGD), JPGD)
  CALL SET_YOR_ll(NYOR_ALL(JPGD), JPGD)
  CALL SET_YEND_ll(NYEND_ALL(JPGD), JPGD)
  ! set the father model of model JPGD
! set MODD_NESTING::NDAD using MODD_DIM_ll::NDAD
! MODD_DIM_ll::NDAD was filled in OPEN_NESTPGD_FILES
  CALL SET_DAD_ll(NDAD(JPGD), JPGD)
  ! set the ratio of model JPGD in MODD_DIM_ll
  CALL SET_XRATIO_ll(NDXRATIO_ALL(JPGD), JPGD)
  CALL SET_YRATIO_ll(NDYRATIO_ALL(JPGD), JPGD)
END DO
!
! reading of the grids
!
  CALL SET_DIM_ll(NXEND_ALL(1)-NXOR_ALL(1)+1-2*JPHEXT, NYEND_ALL(1)-NYOR_ALL(1)+1-2*JPHEXT, 1) 
  CALL INI_PARAZ_ll(IINFO_ll)
DO JPGD=1,NMODEL
  CALL GOTO_MODEL(JPGD)
  CALL GO_TOMODEL_ll(JPGD,IINFO_ll)
  CALL GOTO_SURFEX(JPGD)
  CALL IO_Field_read(TZFILEPGD(JPGD)%TZFILE,'L1D', L1D_ALL(JPGD))
  CALL IO_Field_read(TZFILEPGD(JPGD)%TZFILE,'L2D', L2D_ALL(JPGD))
  CALL IO_Field_read(TZFILEPGD(JPGD)%TZFILE,'PACK',LPACK_ALL(JPGD))
  CALL IO_Pack_set(L1D_ALL(JPGD),L2D_ALL(JPGD),LPACK_ALL(JPGD))
  CALL READ_HGRID(JPGD,TZFILEPGD(JPGD)%TZFILE,YMY_NAME,YDAD_NAME,YSTORAGE_TYPE)
  CSTORAGE_TYPE='PG'
END DO
  CALL INI_PARAZ_ll(IINFO_ll)
!
!-------------------------------------------------------------------------------
!
!*       5.    MASKS DEFINITIONS
!              -----------------
!

DO JPGD=1,NMODEL
  CALL GOTO_SURFEX(JPGD)
  CALL GOTO_MODEL(JPGD)
  CALL GO_TOMODEL_ll(JPGD,IINFO_ll)
!!$  CALL INIT_HORGRID_ll_n()
  CALL DEFINE_MASK_n()
END DO
!
!-------------------------------------------------------------------------------
!
!*       6.    MODIFICATION OF OROGRAPHY
!              -------------------------
!
WRITE(ILUOUT0,FMT=*)
WRITE(ILUOUT0,FMT=*) 'field ZS   of all models'
DO JPGD=NMODEL,1,-1
  CALL GOTO_MODEL(JPGD)
  CALL GO_TOMODEL_ll(JPGD,IINFO_ll)
  CALL GOTO_SURFEX(JPGD)
  CALL NEST_FIELD_n('ZS    ')
END DO
!
! *** Adaptation of smooth topography for SLEVE coordinate
!
WRITE(ILUOUT0,FMT=*)
WRITE(ILUOUT0,FMT=*) 'field ZSMT of all models'
DO JPGD=1,NMODEL
  CALL GOTO_MODEL(JPGD)
  CALL GO_TOMODEL_ll(JPGD,IINFO_ll)
  CALL GOTO_SURFEX(JPGD)
  CALL NEST_ZSMT_n('ZSMT  ')
END DO

!
!-------------------------------------------------------------------------------
!
!*       7.    SURFACE FIELDS READING
!              ----------------------
!
DO JPGD=1,NMODEL
  IF (LEN_TRIM(TZFILEPGD(JPGD)%TZFILE%CNAME)>0) THEN
    CALL GO_TOMODEL_ll(JPGD,IINFO_ll)
    TPGDFILE => TZFILEPGD(JPGD)%TZFILE
    CALL GOTO_MODEL(JPGD)
    CALL GOTO_SURFEX(JPGD)
    CALL INIT_PGD_SURF_ATM(YSURF_CUR,'MESONH','PGD',                         &
         '                            ','      ',&
         NUNDEF,NUNDEF,NUNDEF,XUNDEF             )
  END IF
END DO
!
!-------------------------------------------------------------------------------
!
!*       8.    MODIFICATION OF OROGRAPHY
!              -------------------------
!
DO JPGD=1,NMODEL
  CALL GOTO_MODEL(JPGD)
  CALL GO_TOMODEL_ll(JPGD,IINFO_ll)
  CALL GOTO_SURFEX(JPGD)
  CALL MNHPUT_ZS_n
END DO
!
!-------------------------------------------------------------------------------
!
!*      10.    SURFACE FIELDS WRITING
!              ----------------------
!
DO JPGD=1,NMODEL
  CALL GO_TOMODEL_ll(JPGD,IINFO_ll)
  TPGDFILE     => TZFILEPGD(JPGD)%TZFILE
  TOUTDATAFILE => TZFILENESTPGD(JPGD)%TZFILE
  CALL GOTO_MODEL(JPGD)
  !Open done here because grid dimensions have to be known
  CALL IO_File_open(TZFILENESTPGD(JPGD)%TZFILE)
  CALL GOTO_SURFEX(JPGD)
  TFILE_SURFEX => TZFILENESTPGD(JPGD)%TZFILE
  CALL WRITE_PGD_SURF_ATM_n(YSURF_CUR,'MESONH')
  NULLIFY(TFILE_SURFEX)
  CALL IO_Field_write(TZFILENESTPGD(JPGD)%TZFILE,'ZSMT',XZSMT)
END DO
!
!-------------------------------------------------------------------------------
!
!*      12.    Write configuration variables in the output file
!              ------------------------------------------------
!
!
DO JPGD=1,NMODEL
  CALL IO_Header_write(TZFILENESTPGD(JPGD)%TZFILE)
  IF ( ASSOCIATED(TZFILENESTPGD(JPGD)%TZFILE%TDADFILE) ) THEN
    CALL IO_Field_write(TZFILENESTPGD(JPGD)%TZFILE,'DXRATIO',NDXRATIO_ALL(JPGD))
    CALL IO_Field_write(TZFILENESTPGD(JPGD)%TZFILE,'DYRATIO',NDYRATIO_ALL(JPGD))
    CALL IO_Field_write(TZFILENESTPGD(JPGD)%TZFILE,'XOR',    NXOR_ALL(JPGD))
    CALL IO_Field_write(TZFILENESTPGD(JPGD)%TZFILE,'YOR',    NYOR_ALL(JPGD))
  END IF
  CALL IO_Field_write(TZFILENESTPGD(JPGD)%TZFILE,'SURF',  'EXTE')
  CALL IO_Field_write(TZFILENESTPGD(JPGD)%TZFILE,'L1D',   L1D_ALL(JPGD))
  CALL IO_Field_write(TZFILENESTPGD(JPGD)%TZFILE,'L2D',   L2D_ALL(JPGD))
  CALL IO_Field_write(TZFILENESTPGD(JPGD)%TZFILE,'PACK',  LPACK_ALL(JPGD))
  CALL IO_Field_write(TZFILENESTPGD(JPGD)%TZFILE,'JPHEXT',JPHEXT)
END DO
!
!-------------------------------------------------------------------------------
!
!*      13.    CLOSING OF THE FILES
!              --------------------
!
DO JPGD=1,NMODEL
  CALL IO_File_close(TZFILEPGD(JPGD)%TZFILE)
  CALL IO_File_close(TZFILENESTPGD(JPGD)%TZFILE)
END DO
!
!* loop to spare enough time to transfer commands before end of program
ITIME=0
DO JTIME=1,1000000
  ITIME=ITIME+1
END DO
!-------------------------------------------------------------------------------
!
!*      12.    EPILOGUE
!              --------
!
WRITE(ILUOUT0,FMT=*)
WRITE(ILUOUT0,FMT=*) '************************************************'
WRITE(ILUOUT0,FMT=*) '* PREP_NEST_PGD: PREP_NEST_PGD ends correctly. *'
WRITE(ILUOUT0,FMT=*) '************************************************'
!
!-------------------------------------------------------------------------------
!
!*      10.    FINALIZE THE PARALLEL SESSION
!              -----------------------------
!
CALL FINALIZE_MNH()

! CALL END_PARA_ll(IINFO_ll)
!
! CALL SURFEX_DEALLO_LIST
!
!-------------------------------------------------------------------------------

END PROGRAM PREP_NEST_PGD
