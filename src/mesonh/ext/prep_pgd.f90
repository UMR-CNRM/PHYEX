!MNH_LIC Copyright 1995-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ################
      PROGRAM PREP_PGD
!     ################
!!
!!    PURPOSE
!!    -------
!!   This program prepares the physiographic data fields.
!!
!!    METHOD
!!    ------
!!
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
!!    AUTHOR
!!    ------
!!
!!    F. Mereyde                  Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     21/07/95
!!    Modification 26/07/95       Treatment of orography and subgrid-scale
!!                                orography roughness length (V. Masson)
!!    Modification 22/05/96       Variable CSTORAGE_TYPE (V. Masson)
!!    Modification 25/05/96       Modification of splines, correction on z0rel
!!                                and set limits for some surface varaibles
!!    Modification 12/06/96       Treatment of a rare case for ZPGDZ0EFF (Masson)
!!    Modification 22/11/96       removes the filtering. It will have to be
!!                                performed in ADVANCED_PREP_PGD (Masson)
!!    Modification 15/03/99       **** MAJOR MODIFICATION **** (Masson)
!!                                PGD fields are now defined from the cover
!!                                type fractions in the grid meshes
!!                                User can still include its own data, and
!!                                even additional (dummy) fields
!!    Modificatio 06/00           patch approach, for vegetation related variable (Solmon/Masson)
!                                  averaging is performed on subclass(=patch) of nature
!!                08/03/01        add chemical emission treatment (D.Gazen)
!!    Modification 15/10/01       allow namelists in different orders (I.Mallet)
!!
!!                                ################################
!!    MODIFICATION 13/10/03       EXTERNALIZED VERSION (V. Masson)
!!                                ################################
!!    J.Escobar    4/04/2008      Improve checking --> add STATUS=OLD in open_ll(PRE_PGD1.nam,...
!!
!!    Modification 30/03/2012     Add NAM_NCOUT for netcdf output (S.Bielli)
!!    S.Bielli     23/04/2014     supress writing of LAt and LON in NETCDF case
!!    S.Bielli     20/11/2014     add writing of LAt and LON in NETCDF case
!!    M.Moge       01/03/2015     use MPPDB + SPLIT_GRID is now called in PGD_GRID. Here we extend 
!!                                the new grid on the halo with EXTEND_GRID_ON_HALO (M.Moge)
!!    M.Moge          06/2015     write NDXRATIO,NDYRATIO,NXSIZE,NYSIZE,NXOR,NYOR in .lfi output file
!!    J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!    J.Escobar : 05/10/2015 : missing JPHEXT for LAT/LON/ZS/ZSMT writing
!!    M.Moge          11/2015     disable the creation of files on multiple 
!!                                Z-levels when using parallel IO for PREP_PGD
!!    06/2016     (G.Delautier) phasage surfex 8
!!    P.Wautelet : 08/07/2016 : removed MNH_NCWRIT define
!!    10/2016    (S.Faroux S.Bielli) correction for NHALO=0
!!  01/2018      (G.Delautier) SURFEX 8.1
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!  Q. Rodier 01/2019 : add a new filtering for very high slopes in NAM_ZSFILTER
!  P. Wautelet 07/02/2019: force TYPE to a known value for IO_File_add2list
!  P. Wautelet 07/02/2019: remove OPARALLELIO argument from open and close files subroutines
!                          (nsubfiles_ioz is now determined in IO_File_add2list)
!  P. Wautelet 14/02/2019: remove CLUOUT/CLUOUT0 and associated variables
!  P. Wautelet 06/07/2021: use FINALIZE_MNH
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_CONF,   ONLY : CPROGRAM, L1D, L2D, LPACK, LCARTESIAN
USE MODD_CONF_n,ONLY : CSTORAGE_TYPE
USE MODD_LUNIT,  ONLY : TLUOUT0
USE MODD_LUNIT_n,ONLY : LUNIT_MODEL
USE MODD_PARAMETERS, ONLY : XUNDEF
USE MODD_IO,               only: TFILEDATA, TFILE_OUTPUTLISTING, TFILE_SURFEX
use modd_precision,   only: LFIINT
USE MODD_IO_SURF_MNH, ONLY : NHALO
USE MODD_SPAWN, ONLY : NDXRATIO,NDYRATIO,NXSIZE,NYSIZE,NXOR,NYOR
!
use mode_field,            only: Ini_field_list
USE MODE_FINALIZE_MNH,     only: FINALIZE_MNH
USE MODE_INI_CST,          ONLY: INI_CST
USE MODE_IO,               only: IO_Config_set, IO_Init
USE MODE_IO_FIELD_WRITE,   only: IO_Field_write, IO_Header_write
USE MODE_IO_FILE,          only: IO_File_close, IO_File_open
USE MODE_IO_MANAGE_STRUCT, only: IO_File_add2list
use mode_ll
USE MODE_MODELN_HANDLER
USE MODE_MSG
USE MODE_POS
!
USE MODI_ZSMT_PGD
!
!JUAN
USE MODN_CONFZ
USE MODD_PARAMETERS, ONLY : JPHEXT  
USE MODD_CONF, ONLY       : NHALO_CONF_MNH => NHALO
!JUAN
!
USE MODI_READ_ALL_NAMELISTS
USE MODI_VERSION
USE MODI_PGD_GRID_SURF_ATM
USE MODI_SPLIT_GRID
USE MODI_PGD_SURF_ATM
USE MODI_WRITE_PGD_SURF_ATM_N
USE MODD_MNH_SURFEX_n
!
USE MODE_MPPDB
USE MODI_EXTEND_GRID_ON_HALO
!
USE MODN_CONFIO, ONLY : NAM_CONFIO
!
IMPLICIT NONE
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: IRESP    ! return code for I/O
INTEGER :: ILUOUT0
INTEGER :: ILUNAM
LOGICAL :: GFOUND
CHARACTER(LEN=28) :: YDAD     =' '        ! name of dad of input FM file
CHARACTER(LEN=28) :: CPGDFILE ='PGDFILE'  ! name of the output file
CHARACTER(LEN=100) :: YMSG
INTEGER           :: NZSFILTER=1          ! number of iteration for filter for fine   orography
INTEGER           :: NLOCZSFILTER=3       ! number of iteration for filter of local fine orography
LOGICAL           :: LHSLOP=.FALSE.       ! filtering of local slopes higher than XHSLOP   
REAL              :: XHSLOP=1.0           ! slopes where the local fine filtering is applied   
INTEGER           :: NSLEVE   =12         ! number of iteration for filter for smooth orography
REAL              :: XSMOOTH_ZS = XUNDEF  ! optional uniform smooth orography for SLEVE coordinate
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZWORK ! work array for lat and lon reshape
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZWORK_LAT ! work array for lat and lon reshape
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZWORK_LON ! work array for lat and lon reshape
INTEGER           :: IIMAX, IJMAX
INTEGER           :: NHALO_MNH 
TYPE(TFILEDATA),POINTER :: TZFILE     => NULL()
TYPE(TFILEDATA),POINTER :: TZNMLFILE  => NULL() ! Namelist file
!
NAMELIST/NAM_PGDFILE/CPGDFILE, NHALO
NAMELIST/NAM_ZSFILTER/NZSFILTER,NLOCZSFILTER,LHSLOP,XHSLOP
NAMELIST/NAM_SLEVE/NSLEVE, XSMOOTH_ZS
NAMELIST/NAM_CONF_PGD/JPHEXT, NHALO_MNH
!------------------------------------------------------------------------------
!
CALL MPPDB_INIT()
!
CPROGRAM='PGD   '
!
!*    1.      Set default names and parallelized I/O
!             --------------------------------------
!
CALL IO_Init()
!
NHALO=15
!
CALL IO_File_add2list(TLUOUT0,'OUTPUT_LISTING0','OUTPUTLISTING','WRITE')
CALL IO_File_open(TLUOUT0)
!
!Set output file for PRINT_MSG
TFILE_OUTPUTLISTING => TLUOUT0
!
LUNIT_MODEL(1)%TLUOUT => TLUOUT0
ILUOUT0=TLUOUT0%NLU
!
!JUAN
CALL IO_File_add2list(TZNMLFILE,'PRE_PGD1.nam','NML','READ')
CALL IO_File_open(TZNMLFILE,KRESP=IRESP)
ILUNAM = TZNMLFILE%NLU
IF (IRESP.NE.0 ) THEN
  WRITE(YMSG,*) 'file PRE_PGD1.nam not found, IRESP=', IRESP
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_PGD',YMSG)
ENDIF
!JUAN

CALL POSNAM( TZNMLFILE, 'NAM_PGDFILE', GFOUND )
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_PGDFILE)
CALL POSNAM( TZNMLFILE, 'NAM_ZSFILTER', GFOUND )
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_ZSFILTER)
CALL POSNAM( TZNMLFILE, 'NAM_SLEVE', GFOUND )
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_SLEVE)
!JUANZ
CALL POSNAM( TZNMLFILE, 'NAM_CONFZ', GFOUND )
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_CONFZ)
CALL POSNAM( TZNMLFILE, 'NAM_CONF_PGD', GFOUND )
IF (GFOUND) THEN
   NHALO_MNH = NHALO_CONF_MNH
   READ(UNIT=ILUNAM,NML=NAM_CONF_PGD)
   NHALO_CONF_MNH = NHALO_MNH
ENDIF
!JUANZ
CALL POSNAM( TZNMLFILE, 'NAM_CONFIO', GFOUND )
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_CONFIO)
CALL IO_Config_set()
!
CALL IO_File_close(TZNMLFILE)
!
!
CALL SURFEX_ALLOC_LIST(1)
YSURF_CUR => YSURF_LIST(1)
CALL READ_ALL_NAMELISTS(YSURF_CUR,'MESONH','PRE',.FALSE.)
!
CALL INI_FIELD_LIST()
!
CALL GOTO_MODEL(1)
CALL GOTO_SURFEX(1)
!
CALL VERSION
CSTORAGE_TYPE = 'PG'
!
CALL INI_CST
!
!
!*    2.      Preparation of surface physiographic fields
!             -------------------------------------------
!
!*            Initializes the grid
!             --------------------
! 
CALL PGD_GRID_SURF_ATM(YSURF_CUR%UG, YSURF_CUR%U,YSURF_CUR%GCP,'MESONH',&
                       '                            ','      ',.FALSE.,HDIR='-')
!
CALL EXTEND_GRID_ON_HALO('MESONH',YSURF_CUR%UG, YSURF_CUR%U,&
        YSURF_CUR%UG%G%NGRID_PAR, YSURF_CUR%UG%G%XGRID_PAR)
!
!
!*            Initializes all physiographic fields
!             ------------------------------------
!
CALL PGD_SURF_ATM(YSURF_CUR,'MESONH','                            ','      ',.FALSE.)
!
!
!*    3.      Writes the physiographic fields
!             -------------------------------
!
CALL IO_File_add2list(TZFILE,CPGDFILE,'PGD','WRITE',KLFINPRAR=INT(1,KIND=LFIINT),KLFITYPE=1,KLFIVERB=5)
!
CALL IO_File_open(TZFILE)
!
CALL IO_Header_write(TZFILE)
!
CALL IO_Field_write(TZFILE,'SURF','EXTE')
CALL IO_Field_write(TZFILE,'L1D', L1D)
CALL IO_Field_write(TZFILE,'L2D', L2D)
CALL IO_Field_write(TZFILE,'PACK',LPACK)
IF ( NDXRATIO <= 0 .AND. NDYRATIO <= 0 ) THEN
  NDXRATIO = 1
  NDYRATIO = 1
ENDIF
IF ( NXSIZE < 0 .AND. NYSIZE < 0 ) THEN
  NXSIZE = 0
  NYSIZE = 0
ENDIF
IF ( NXOR <= 0 .AND. NYOR <= 0 ) THEN
  NXOR = 1
  NYOR = 1
ENDIF
CALL IO_Field_write(TZFILE,'DXRATIO',NDXRATIO)
CALL IO_Field_write(TZFILE,'DYRATIO',NDYRATIO)
CALL IO_Field_write(TZFILE,'XSIZE',  NXSIZE)
CALL IO_Field_write(TZFILE,'YSIZE',  NYSIZE)
CALL IO_Field_write(TZFILE,'XOR',    NXOR)
CALL IO_Field_write(TZFILE,'YOR',    NYOR)
CALL IO_Field_write(TZFILE,'JPHEXT', JPHEXT)
!
TFILE_SURFEX => TZFILE
ALLOCATE(YSURF_CUR%DUO%CSELECT(0))
CALL WRITE_PGD_SURF_ATM_n(YSURF_CUR,'MESONH')
NULLIFY(TFILE_SURFEX) !Probably not necessary
!
!*    4.      Computes and writes smooth orography for SLEVE coordinate
!             ---------------------------------------------------------
CALL ZSMT_PGD(TZFILE,NZSFILTER,NSLEVE,NLOCZSFILTER,LHSLOP,XHSLOP,XSMOOTH_ZS)
!
IF (.NOT.LCARTESIAN) THEN
!!!! WRITE LAT and LON
   CALL GET_DIM_PHYS_ll('B',IIMAX,IJMAX)
   ALLOCATE(ZWORK(IIMAX+NHALO*2,IJMAX+NHALO*2))
   ALLOCATE(ZWORK_LAT(IIMAX+2*JPHEXT,IJMAX+2*JPHEXT))
   ALLOCATE(ZWORK_LON(IIMAX+2*JPHEXT,IJMAX+2*JPHEXT))
   ZWORK=RESHAPE(YSURF_CUR%UG%G%XLAT, (/ (IIMAX+NHALO*2),(IJMAX+NHALO*2) /) )
   IF (NHALO/=0) THEN   
     ZWORK_LAT=ZWORK(NHALO:(IIMAX+NHALO+1),NHALO:(IJMAX+NHALO+1))
   ELSE
     ZWORK_LAT(2:IIMAX+1,2:IJMAX+1)=ZWORK
     ZWORK_LAT(1,:) = ZWORK_LAT(2,:)
     ZWORK_LAT(IIMAX+2,:) = ZWORK_LAT(IIMAX+1,:)
     ZWORK_LAT(:,1) = ZWORK_LAT(:,2)
     ZWORK_LAT(:,IJMAX+2) = ZWORK_LAT(:,IJMAX+1)
   ENDIF
   ZWORK=RESHAPE(YSURF_CUR%UG%G%XLON, (/ IIMAX+NHALO*2,IJMAX+NHALO*2 /) )
   IF (NHALO/=0) THEN   
     ZWORK_LON=ZWORK(NHALO:(IIMAX+NHALO+1),NHALO:(IJMAX+NHALO+1))
   ELSE
     ZWORK_LON(2:IIMAX+1,2:IJMAX+1)=ZWORK
     ZWORK_LON(1,:) = ZWORK_LON(2,:)
     ZWORK_LON(IIMAX+2,:) = ZWORK_LON(IIMAX+1,:)
     ZWORK_LON(:,1) = ZWORK_LON(:,2)
     ZWORK_LON(:,IJMAX+2) = ZWORK_LON(:,IJMAX+1)           
   ENDIF   
   CALL IO_Field_write(TZFILE,'LAT',ZWORK_LAT)
   CALL IO_Field_write(TZFILE,'LON',ZWORK_LON)
   !
   DEALLOCATE(ZWORK,ZWORK_LAT,ZWORK_LON)
END IF
!
!
WRITE(ILUOUT0,*)
WRITE(ILUOUT0,*) '***************************'
WRITE(ILUOUT0,*) '* PREP_PGD ends correctly *'
WRITE(ILUOUT0,*) '***************************'
!
!*    6.      Close parallelized I/O
!             ----------------------
!
CALL IO_File_close(TZFILE)
!
CALL FINALIZE_MNH()
!
!-------------------------------------------------------------------------------
!
END PROGRAM PREP_PGD
