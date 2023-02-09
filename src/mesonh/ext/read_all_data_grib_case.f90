!MNH_LIC Copyright 1998-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################################
      MODULE MODI_READ_ALL_DATA_GRIB_CASE
!     #################################
INTERFACE
SUBROUTINE READ_ALL_DATA_GRIB_CASE(HFILE,TPPRE_REAL1,HGRIB,TPPGDFILE,    &
                    PTIME_HORI,KVERB,ODUMMY_REAL                         ) 
!
USE MODD_IO, ONLY: TFILEDATA
!
CHARACTER(LEN=4),  INTENT(IN)    :: HFILE       ! which file ('ATM0','ATM1' or 'CHEM')
TYPE(TFILEDATA),POINTER,INTENT(INOUT) :: TPPRE_REAL1 ! PRE_REAL1 file
CHARACTER(LEN=28), INTENT(IN)    :: HGRIB       ! name of the GRIB file
TYPE(TFILEDATA),   INTENT(IN)    :: TPPGDFILE   ! physiographic data file
INTEGER,           INTENT(IN)    :: KVERB       ! verbosity level
LOGICAL,           INTENT(IN)    :: ODUMMY_REAL ! flag to interpolate dummy fields
REAL,              INTENT(INOUT) :: PTIME_HORI  ! time spent in hor. interpolations
!
END SUBROUTINE READ_ALL_DATA_GRIB_CASE
!
END INTERFACE
END MODULE MODI_READ_ALL_DATA_GRIB_CASE
!     ##########################################################################
      SUBROUTINE READ_ALL_DATA_GRIB_CASE(HFILE,TPPRE_REAL1,HGRIB,TPPGDFILE,    &
                       PTIME_HORI,KVERB,ODUMMY_REAL                            )
!     ##########################################################################
!
!!****  *READ_ALL_DATA_GRIB_CASE* - reads data for the initialization of real cases.
!!
!!    PURPOSE
!!    -------
!     This routine reads the two input files :
!       The PGD which is closed after reading
!       The GRIB file
!     Projection is read in READ_LFIFM_PGD (MODD_GRID).
!     Grid and definition of large domain are read in PGD file and Grib files.
!     The PGD files are also read in READ_LFIFM_PGD.
!     The PGD file is closed.
!     The MESO-NH domain is defined from PRE_REAL1.nam inputs in SET_SUBDOMAIN_CEP.
!     Vertical grid is defined in READ_VER_GRID.
!     PGD fields are stored on MESO-NH domain (in TRUNC_PGD).
!!
!!**  METHOD
!!    ------
!!  0. Declarations
!!    1. Declaration of arguments
!!    2. Declaration of local variables
!!  1. Read PGD file
!!    1. Domain restriction
!!    2. Coordinate conversion to lat,lon system
!!  2. Read Grib fields
!!  3. Vertical grid
!!  4. Free all temporary allocations
!!
!!    EXTERNAL
!!    --------
!!    subroutine READ_LFIFM_PGD    : to read PGD file
!!    subroutine SET_SUBDOMAIN     : to define the horizontal MESO-NH domain.
!!    subroutine READ_VER_GRID     : to read the vertical grid in namelist file.
!!    subroutine HORIBL            : horizontal bilinear interpolation
!!    subroutine XYTOLATLON        : projection from conformal to lat,lon
!!
!!    Module     MODI_SET_SUBDOMAIN     : interface for subroutine SET_SUBDOMAIN
!!    Module     MODI_READ_VER_GRID     : interface for subroutine READ_VER_GRID
!!    Module     MODI_HORIBL            : interface for subroutine HORIBL
!!    Module     MODI_XYTOLATLON        : interface for subroutine XYTOLATLON
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     : contains logical unit names for all models
!!         TLUOUT0 : name of output-listing
!!      Module MODD_PGDDIM    : contains dimension of PGD fields
!!         NPGDIMAX: dimension along x (no external point)
!!         NPGDJMAX: dimension along y (no external point)
!!      Module MODD_PARAMETERS
!!         JPHEXT
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 1 : Informations on ISBA model (soil moisture)
!!      "Encoding and decoding Grib data", John D.Chambers, ECMWF, October 95
!!      "A guide to Grib", John D.Stackpole, National weather service, March 94
!!
!!    AUTHOR
!!    ------
!!
!!      J. Pettre and V. Bousquet
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/11/98
!!                  15/03/99 (V. Masson) phasing with new PGD fields
!!                  21/04/99 (V. Masson) bug in mask definitions for max Y index
!!                  22/04/99 (V. Masson) optimizer bug in u,v loop
!!                                       --> splitting of the loop
!!                                       and splitting of the routine in more
!!                                       contains
!!                  28/05/99 (V. Bousquet) bug in wind interpolated variable for
!!                                         Arpege
!!                  31/05/99 (V. Masson) set pressure points (given on a regular grid at ECMWF)
!!                             on orography points (assuming the last are included in the former)
!!                                       pressure computation from parameters A and B
!!                                        (instead of interpolation from grib grid)
!!                  20/07/00 (V. Masson) increase the threshold for land_sea index
!!                  22/11/00 (P. Tulet) add INTERPOL_SV to initialize SV fields
!!                           (I. Mallet) from MOCAGE model (IMODE=3)
!!                  01/02/01 (D. Gazen) add INI_NSV
!!                  18/05/01 (P. Jabouille) problem with 129 grib code
!!                  05/12/01 (I. Mallet) add Aladin reunion model
!!                  02/10/02 (I. Mallet) 2 orography fields for CEP (SFC, ML=1)
!!                  01/12/03 (D. Gazen)   change Chemical scheme interface
!!                  01/2004  (V. Masson) removes surface (externalization)
!!                  01/06/02 (O.Nuissier) filtering of tropical cyclone
!!                  01/05/04 (P. Tulet) add INTERPOL_SV to initialize SV dust
!!                                                         and aerosol fields
!!                  08/06/2010 (G. Tanguy) replace GRIBEX by GRIB_API : change
!!                                         of all the subroutine
!!                  05/12/2016 (G.Delautier) length of HGRID for grib_api > 1.14
!!                  08/03/2018 (P.Wautelet)  replace ADD_FORECAST_TO_DATE by DATETIME_CORRECTDATE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!         Pergaud  : 2018 add GFS
!!                   01/2019 (G.Delautier via Q.Rodier) for GRIB2 ARPEGE and AROME from EPYGRAM
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 14/03/2019: correct ZWS when variable not present in file
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  Q. Rodier   16/09/2019: switch of GRIB number ID for orography in ARPEGE/AROME in EPyGrAM
!  Q. Rodier   27/01/2020: switch of GRIB number ID for orography and hydrometeors in ARPEGE/AROME in EPyGrAM v1.3.7
!  Q. Rodier   21/04/2020: correction GFS u and v wind component written in the right vertical order
!  Q. Rodier   02/09/2020: Read and interpol geopotential height for interpolation on isobaric surface Grid of NCEP
!  P. Wautelet 09/03/2021: move some chemistry initializations to ini_nsv
!JP Chaboureau 02/08/2021: add ERA5 reanalysis in pressure levels
!JP Chaboureau 18/10/2022: correction on vertical level for GFS and ERA5 reanalyses in pressure levels
!-------------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!------------
!
USE MODE_DATETIME
USE MODE_IO_FILE,    ONLY: IO_File_close
USE MODE_MSG
USE MODE_TIME
USE MODE_THERMO
USE MODE_TOOLS, ONLY: UPCASE
use mode_tools_ll,   only: GET_DIM_EXT_ll
!
USE MODI_READ_HGRID_n
USE MODI_READ_VER_GRID
USE MODI_XYTOLATLON
USE MODI_HORIBL
USE MODI_INI_NSV
USE MODI_REMOVAL_VORTEX
USE MODI_INI_CTURB
USE MODI_CH_OPEN_INPUT
!
USE MODD_IO, ONLY: TFILEDATA
USE MODD_FIELD_n, ONLY: XZWS, XZWS_DEFAULT
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST
USE MODD_LUNIT
USE MODD_PARAMETERS
USE MODD_GRID
USE MODD_GRID_n
USE MODD_DIM_n
USE MODD_PARAM_n, ONLY : CTURB
USE MODD_TIME
USE MODD_TIME_n
USE MODD_CH_MNHC_n, ONLY : LUSECHEM,LUSECHAQ,LUSECHIC,LCH_PH
USE MODD_CH_M9_n,   ONLY : NEQ ,  CNAMES
USE MODD_CH_AEROSOL, ONLY: CORGANIC, NCARB, NSOA, NSP, LORILAM,&
                           JPMODE, LVARSIGI, LVARSIGJ
USE MODD_NSV      , ONLY : NSV
USE MODD_HURR_CONF, ONLY : LFILTERING,CFILTERING
USE MODD_PREP_REAL
USE MODE_MODELN_HANDLER
!JUAN REALZ
USE MODE_MPPDB
!JUAN REALZ
!
USE GRIB_API
!
IMPLICIT NONE
!
!* 0.1. Declaration of arguments
!       ------------------------
!
CHARACTER(LEN=4),  INTENT(IN)    :: HFILE       ! which file ('ATM0','ATM1' or 'CHEM')
TYPE(TFILEDATA),POINTER,INTENT(INOUT) :: TPPRE_REAL1! PRE_REAL1 file
CHARACTER(LEN=28), INTENT(IN)    :: HGRIB       ! name of the GRIB file
TYPE(TFILEDATA),   INTENT(IN)    :: TPPGDFILE   ! physiographic data file
INTEGER,           INTENT(IN)    :: KVERB       ! verbosity level
LOGICAL,           INTENT(IN)    :: ODUMMY_REAL ! flag to interpolate dummy fields
REAL,              INTENT(INOUT) :: PTIME_HORI  ! time spent in hor. interpolations
!
!* 0.2 Declaration of local variables
!      ------------------------------
! General purpose variables
INTEGER                            :: ILUOUT0       ! Unit used for output msg.
INTEGER                            :: IRESP         ! Return code of FM-routines
INTEGER                            :: IRET          ! Return code from subroutines
INTEGER(KIND=kindOfInt)            :: IRET_GRIB     ! Return code from subroutines
INTEGER, PARAMETER                 :: JP_GFS=31     ! number of pressure levels for GFS model
INTEGER, PARAMETER                 :: JP_ERA=37     ! number of pressure levels for ERA5 reanalysis
REAL                               :: ZA,ZB,ZC      ! Dummy variables
REAL                               :: ZD,ZE,ZF      !  |
REAL                               :: ZTEMP         !  |
INTEGER                            :: JI,JJ         ! Dummy counters
INTEGER                            :: JLOOP1,JLOOP2 !  |
INTEGER                            :: JLOOP3,JLOOP4 !  |
INTEGER                            :: JLOOP         !  |
! Variables used by the PGD reader
CHARACTER(LEN=28)                  :: YPGD_NAME     ! not used - dummy argument
CHARACTER(LEN=28)                  :: YPGD_DAD_NAME ! not used - dummy argument
CHARACTER(LEN=2)                   :: YPGD_TYPE     ! not used - dummy argument
! PGD Grib definition variables
INTEGER                            :: INO           ! Number of points of the grid
INTEGER                            :: IIU           ! Number of points along X
INTEGER                            :: IJU           ! Number of points along Y
REAL, DIMENSION(:), ALLOCATABLE    :: ZXOUT         ! mapping PGD -> Grib (lon.)
REAL, DIMENSION(:), ALLOCATABLE    :: ZYOUT         ! mapping PGD -> Grib (lat.)
REAL, DIMENSION(:), ALLOCATABLE    :: ZLONOUT       ! mapping PGD -> Grib (lon.)
REAL, DIMENSION(:), ALLOCATABLE    :: ZLATOUT       ! mapping PGD -> Grib (lat.)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZXM           ! X of PGD mass points
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZYM           ! Y of PGD mass points
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZLATM         ! Lat of PGD mass points
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZLONM         ! Lon of PGD mass points
! Variable involved in the task of reading the grib file
INTEGER(KIND=kindOfInt)            :: IUNIT         ! unit of the grib file
CHARACTER(LEN=50)                  :: HGRID         ! type of grid
INTEGER                            :: IPAR          ! Parameter identifier
INTEGER                            :: ITYP          ! type of level (Grib code table 3)
INTEGER                            :: ILEV1         ! level definition
INTEGER                            :: ILEV2         ! level definition
INTEGER                            :: IMODEL        ! Type of Grib file :
                                                    ! 0 -> ECMWF
                                                    ! 1 -> METEO FRANCE - ALADIN/AROME
                                                    ! 2 -> METEO FRANCE - ALADIN-REUNION
                                                    ! 3 -> METEO FRANCE - ARPEGE
                                                    ! 4 -> METEO FRANCE - ARPEGE
                                                    ! 5 -> METEO FRANCE - MOCAGE
                                                    ! 10 -> NCEP - GFS
INTEGER                            :: ICENTER       ! number of center
INTEGER                            :: ISIZE         ! size of grib message
INTEGER(KIND=kindOfInt)            :: ICOUNT        ! number of messages in the file
INTEGER(KIND=kindOfInt),DIMENSION(:),ALLOCATABLE   :: IGRIB         ! number of the grib in memory
INTEGER                            :: INUM ,INUM_ZS ! number of a grib message 
REAL,DIMENSION(:),ALLOCATABLE      :: ZPARAM        ! parameter of grib grid
INTEGER,DIMENSION(:),ALLOCATABLE   :: IINLO         ! longitude of grib grid
INTEGER(KIND=kindOfInt),DIMENSION(:),ALLOCATABLE   :: IINLO_GRIB         ! longitude of grib grid
REAL,DIMENSION(:),ALLOCATABLE      :: ZPARAM_ZS     ! parameter of grib grid for ZS
INTEGER,DIMENSION(:),ALLOCATABLE   :: IINLO_ZS      ! longitude of grib grid for ZS
REAL,DIMENSION(:),ALLOCATABLE      :: ZVALUE        ! Intermediate array
REAL,DIMENSION(:),ALLOCATABLE      :: ZOUT          ! Intermediate arrays
! Grib grid definition variables
INTEGER                            :: INI           ! Number of points
INTEGER                            :: INLEVEL       ! Number of levels
INTEGER                            :: ISTARTLEVEL   ! First level (0 or 1)
TYPE(DATE_TIME)                    :: TPTCUR        ! Date & time of the grib file data
INTEGER                            :: ITWOZS
! surface pressure
REAL, DIMENSION(:), ALLOCATABLE    :: ZPS_G         ! Grib data : Ps
REAL, DIMENSION(:), ALLOCATABLE    :: ZLNPS_G       ! Grib data : ln(Ps)
REAL, DIMENSION(:), ALLOCATABLE    :: ZWORK_LNPS    ! Grib data on zs grid: ln(Ps)
INTEGER                            :: INJ,INJ_ZS
! orography
CHARACTER(LEN=50)                  :: HGRID_ZS      ! type of grid
!
! Reading and projection of the wind vectors u, v
REAL                               :: ZALPHA        ! Angle of rotation
REAL, DIMENSION(:), ALLOCATABLE    :: ZTU_LS        ! Intermediate array for U
REAL, DIMENSION(:), ALLOCATABLE    :: ZTV_LS        !  |                     V
REAL                               :: ZLATPOLE      ! Arpege stretching pole latitude
REAL                               :: ZLONPOLE      ! Arpege stretching pole longitude
REAL                               :: ZLAT,ZLON     ! Lat,lon of current point
REAL                               :: ZCOS,ZSIN     ! cos,sin of rotation matrix
REAL, DIMENSION(:), ALLOCATABLE    :: ZTU0_LS       ! Arpege temp array for U
REAL, DIMENSION(:), ALLOCATABLE    :: ZTV0_LS       !  |                    V
!
! variables for hurricane filtering
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTVF_LS,ZMSLP_LS
REAL                               :: ZGAMREF       ! Standard atmosphere lapse rate (K/m)
! date
INTEGER  :: ITIME
INTEGER  :: IDATE
INTEGER  :: ITIMESTEP
CHARACTER(LEN=10) :: CSTEPUNIT
CHARACTER(LEN=15)  :: YVAL
!chemistery field
CHARACTER(LEN=16)                  :: YPRE_MOC="PRE_MOC1.nam"
INTEGER, DIMENSION(:), ALLOCATABLE :: INUMGRIB, INUMLEV  ! grib
INTEGER, DIMENSION(:), ALLOCATABLE :: INUMLEV1, INUMLEV2 !numbers
INTEGER                            :: IMOC
INTEGER                            :: IVAR
INTEGER                            :: ICHANNEL
INTEGER                            :: INDX
INTEGER                            :: INACT
CHARACTER(LEN=40)                  :: YINPLINE        ! input line
CHARACTER(LEN=16)                  :: YFIELD
CHARACTER, PARAMETER               :: YPTAB = CHAR(9) ! TAB character is ASCII : 9
CHARACTER, PARAMETER               :: YPCOM = CHAR(44)! COMma character is ASCII : 44
CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE :: YMNHNAME ! species names
INTEGER                            :: JN, JNREAL ! loop control variables
CHARACTER(LEN=40)                  :: YFORMAT
CHARACTER(LEN=100)                 :: YMSG
! temperature and humidity
INTEGER                             :: IT,IQ
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZPF_G    ! Pressure (flux point)
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZPM_G    ! Pressure (mass point)
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZEXNF_G  ! Exner fct. (flux point)
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZEXNM_G  ! Exner fct. (mass point)
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZGH_G     ! Geopotential Height
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZT_G     ! Temperature
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZQ_G     ! Specific humidity
REAL, DIMENSION(:), ALLOCATABLE     :: ZH_G     ! Relative humidity
REAL, DIMENSION(:), ALLOCATABLE     :: ZTHV_G   ! Theta V
REAL, DIMENSION(:), ALLOCATABLE     :: ZRV_G    ! Vapor mixing ratio
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZPF_LS   ! Pressure (flux point)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZPM_LS   ! Pressure (mass point)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZEXNF_LS ! Exner fct. (flux point)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZEXNM_LS ! Exner fct. (mass point)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZH_LS    ! Relative humidity
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZRV_LS   ! Vapor mixing ratio
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTHV_LS  ! Theta V
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTEV_LS  ! T V
REAL, DIMENSION(:), ALLOCATABLE     :: ZPV      ! vertical level in grib file
INTEGER                             :: IPVPRESENT ,IPV
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZR_DUM
INTEGER                             :: IMI
TYPE(TFILEDATA),POINTER             :: TZFILE
INTEGER, DIMENSION(JP_GFS)    :: IP_GFS   ! list of pressure levels for GFS model
INTEGER, DIMENSION(JP_ERA)    :: IP_ERA   ! list of pressure levels for ERA5 reanalysis
INTEGER :: IVERSION,ILEVTYPE
LOGICAL                       :: GFIND    ! to test if sea wave height is found
!---------------------------------------------------------------------------------------
IP_GFS=(/1000,975,950,925,900,850,800,750,700,650,600,550,500,450,400,350,300,&
           250,200,150,100,70,50,30,20,10,7,5,3,2,1/)
IP_ERA=(/1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,&
           400,350,300,250,225,200,175,150,125,100,70,50,30,20,10,7,5,3,2,1/)
!
TZFILE => NULL()
!
IMI = GET_CURRENT_MODEL_INDEX()
!
!* 1. READ PGD FILE
!     -------------
!
ILUOUT0 = TLUOUT0%NLU
CALL READ_HGRID_n(TPPGDFILE,YPGD_NAME,YPGD_DAD_NAME,YPGD_TYPE)
!
! 1.1 Domain restriction
!
!JUAN REALZ
CALL GET_DIM_EXT_ll('B',IIU,IJU)
!!$IIU=NIMAX + 2*JPHEXT
!!$IJU=NJMAX + 2*JPHEXT
!JUAN REALZ
INO = IIU * IJU
!
!
! 1.2 Coordinate conversion to lat,lon system
!
ALLOCATE (ZXM(IIU,IJU))
ALLOCATE (ZYM(IIU,IJU))
ALLOCATE (ZLONM(IIU,IJU))
ALLOCATE (ZLATM(IIU,IJU))
ZXM(1:IIU-1,1) = (XXHAT(1:IIU-1) + XXHAT(2:IIU) ) / 2.
ZXM(IIU,1)     = XXHAT(IIU) - XXHAT(IIU-1) + ZXM(IIU-1,1)
ZXM(:,2:IJU)   = SPREAD(ZXM(:,1),2,IJU-1)
ZYM(1,1:IJU-1) = (XYHAT(1:IJU-1) + XYHAT(2:IJU)) / 2.
ZYM(1,IJU)     = XYHAT(IJU) - XYHAT(IJU-1) + ZYM(1,IJU-1)
ZYM(2:IIU,:)   = SPREAD(ZYM(1,:),1,IIU-1)
CALL SM_XYTOLATLON_A (XLAT0,XLON0,XRPK,XLATORI,XLONORI,ZXM,ZYM,ZLATM,ZLONM, &
                      IIU,IJU)
ALLOCATE (ZLONOUT(INO))
ALLOCATE (ZLATOUT(INO))
JLOOP1 = 0
DO JJ = 1, IJU
  ZLONOUT(JLOOP1+1:JLOOP1+IIU) = ZLONM(1:IIU,JJ)
  ZLATOUT(JLOOP1+1:JLOOP1+IIU) = ZLATM(1:IIU,JJ)
  JLOOP1 = JLOOP1 + IIU
ENDDO
DEALLOCATE (ZLATM)
DEALLOCATE (ZLONM)
DEALLOCATE (ZYM)
DEALLOCATE (ZXM)
!
ALLOCATE (ZXOUT(INO))
ALLOCATE (ZYOUT(INO))
!
!---------------------------------------------------------------------------------------
!
!* 2. READ GRIB FIELDS
!     ----------------
!
IF (HFILE(1:3)=='ATM' .OR. HFILE=='CHEM') THEN
  WRITE (ILUOUT0,'(A,A4)') ' -- Grib reader started for ',HFILE
ELSE
  CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE','bad input argument')
END IF
!
!* 2.1 Charge in memory the grib messages
!
! open grib file
CALL GRIB_OPEN_FILE(IUNIT,HGRIB,'R',IRET_GRIB)
IF (IRET_GRIB /= 0) THEN
  WRITE(YMSG,*) 'Error opening the grib file ',TRIM(HGRIB),', error code ', IRET_GRIB
  CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
END IF
! count the messages in the file
CALL GRIB_COUNT_IN_FILE(IUNIT,ICOUNT,IRET_GRIB)
IF (IRET_GRIB /= 0) THEN
  WRITE(YMSG,*) 'Error in reading the grib file ',TRIM(HGRIB),', error code ', IRET_GRIB
  CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
END IF
ALLOCATE(IGRIB(ICOUNT))
! initialize the tabular with a negativ number 
! ( all the IGRIB will be  different )
IGRIB(:)=-12
!charge all the message in memory
DO JLOOP=1,ICOUNT
CALL GRIB_NEW_FROM_FILE(IUNIT,IGRIB(JLOOP),IRET_GRIB)
IF (IRET_GRIB /= 0) THEN
  WRITE(YMSG,*) 'Error in reading the grib file - ILOOP=',JLOOP,' - error code ', IRET_GRIB
  CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
END IF
END DO
! close the grib file
CALL GRIB_CLOSE_FILE(IUNIT)
!
!
!---------------------------------------------------------------------------------------
!* 2.2 Research center with the first message
!---------------------------------------------------------------------------------------
!
CALL GRIB_GET(IGRIB(1),'centre',ICENTER,IRET_GRIB)
IF (IRET_GRIB /= 0) THEN
  WRITE(YMSG,*) 'Error in reading center - error code ', IRET_GRIB
  CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
END IF
CALL GRIB_GET(IGRIB(1),'typeOfGrid',HGRID,IRET_GRIB)
IF (IRET_GRIB /= 0) THEN
  WRITE(YMSG,*) 'Error in reading type of grid - error code ', IRET_GRIB
  CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
END IF
!
IMODEL = -1
SELECT CASE (ICENTER)
  CASE (98)
    WRITE (ILUOUT0,'(A)') &
        ' | Grib file from European Center for Medium-range Weather Forecast'
    IMODEL = 0
    ALLOCATE(ZPARAM(6))
  CASE (85)
    SELECT CASE (HGRID)      
      CASE('lambert')
        WRITE (ILUOUT0,'(A)') ' | Grib file from French Weather Service - Arome france model'
        CALL GRIB_GET(IGRIB(1),'editionNumber',IVERSION,IRET_GRIB)
        IF (IVERSION==2) THEN
          IMODEL = 6 ! GRIB2 since summer 2018 (epygram)
        ELSE
          IMODEL = 1 ! GRIB1 befor summer 2018 (lfi2mv) 
        ENDIF
        ALLOCATE(ZPARAM(10))
      CASE('mercator')
        WRITE (ILUOUT0,'(A)') ' | Grib file from French Weather Service - Aladin reunion model'
      IMODEL = 2
      ALLOCATE(ZPARAM(9))

    CASE('unknown_PLPresent','reduced_stretched_rotated_gg')
      WRITE (ILUOUT0,'(A)') ' | Grib file from French Weather Service - Arpege model'
      ALLOCATE(ZPARAM(10))
      CALL GRIB_GET(IGRIB(1),'editionNumber',IVERSION,IRET_GRIB)
        IF (IVERSION==2) THEN
          IMODEL = 7 ! GRIB2 since summer 2018 (epygram)
        ELSE
          IMODEL = 3 ! GRIB1 befor summer 2018 (lfi2mv)
        ENDIF

    CASE('regular_gg')
      WRITE (ILUOUT0,'(A)') ' | Grib file from French Weather Service - Arpege model'
      WRITE (ILUOUT0,'(A)') 'but same grid as ECMWF model (unstretched)'
      IMODEL = 4
      ALLOCATE(ZPARAM(10))
    CASE('regular_ll')  
      WRITE (ILUOUT0,'(A)') ' | Grib file from French Weather Service - Mocage model'
      IMODEL = 5
       ALLOCATE(ZPARAM(6))
    END SELECT    
  CASE (7)
    WRITE (ILUOUT0,'(A)') ' | Grib file from National Center for Environmental Prediction'
    IMODEL = 10
    ALLOCATE(ZPARAM(6))
END SELECT
IF (IMODEL==-1) THEN
  CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE','unsupported Grib file format')
END IF
!
!---------------------------------------------------------------------------------------
!* 2.3 Read and interpol orography 
!---------------------------------------------------------------------------------------
!
WRITE (ILUOUT0,'(A)') ' | Searching orography'
SELECT CASE (IMODEL)
  CASE(0) ! ECMWF
     CALL SEARCH_FIELD(IGRIB,INUM_ZS,KPARAM=129)
     IF(INUM_ZS < 0) THEN
       WRITE (ILUOUT0,'(A)')'Orography is missing - abort'
     END IF 
  CASE(3,4,5) ! arpege et mocage
      CALL SEARCH_FIELD(IGRIB,INUM_ZS,KPARAM=8)
       IF(INUM_ZS < 0) THEN
         WRITE (ILUOUT0,'(A)')'Orography is missing - abort'
       ENDIF 
  CASE(1,2) ! aladin et aladin reunion
      CALL SEARCH_FIELD(IGRIB,INUM_ZS,KPARAM=6)
       IF(INUM_ZS < 0) THEN
         WRITE (ILUOUT0,'(A)')'Orography is missing - abort'
       ENDIF 
  CASE(6,7) !  arpege and arome GRIB2
      CALL SEARCH_FIELD(IGRIB,INUM_ZS,KDIS=0,KCAT=3,KNUMBER=4)
       IF(INUM_ZS < 0) THEN
         CALL SEARCH_FIELD(IGRIB,INUM_ZS,KDIS=0,KCAT=193,KNUMBER=5)
         IF(INUM_ZS < 0) THEN
           WRITE (ILUOUT0,'(A)')'Orography is missing - abort'
         END IF
       ENDIF 
  CASE(10) ! NCEP
       CALL SEARCH_FIELD(IGRIB,INUM_ZS,KDIS=0,KCAT=3,KNUMBER=5,KTFFS=1)
        IF(INUM_ZS < 0) THEN
          WRITE (ILUOUT0,'(A)')'Orography is missing - abort'
        ENDIF    
END SELECT
ZPARAM(:)=-999.
CALL GRIB_GET(IGRIB(INUM_ZS),'Nj',INJ,IRET_GRIB)
ALLOCATE(IINLO(INJ))
CALL COORDINATE_CONVERSION(IMODEL,IGRIB(INUM_ZS),IIU,IJU,ZLONOUT,ZLATOUT,&
        ZXOUT,ZYOUT,INI,ZPARAM,IINLO)
ALLOCATE(ZPARAM_ZS(SIZE(ZPARAM)))
ZPARAM_ZS=ZPARAM
ALLOCATE(IINLO_ZS(SIZE(IINLO)))
IINLO_ZS=IINLO
CALL GRIB_GET_SIZE(IGRIB(INUM_ZS),'values',ISIZE)
ALLOCATE(ZVALUE(ISIZE))
CALL GRIB_GET(IGRIB(INUM_ZS),'values',ZVALUE)
ALLOCATE(ZOUT(INO))
CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
            ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
DEALLOCATE(IINLO)
DEALLOCATE(ZVALUE)
!
IF (IMODEL/=10) THEN ! others than NCEP
  ! Data given in archives are multiplied by the gravity acceleration
  ZOUT = ZOUT / XG
END IF
!
! Stores the field in a 2 dimension array
IF (HFILE(1:3)=='ATM') THEN
  ALLOCATE (XZS_LS(IIU,IJU))
  ALLOCATE (XZSMT_LS(IIU,IJU))
  CALL ARRAY_1D_TO_2D (INO,ZOUT,IIU,IJU,XZS_LS)
  XZSMT_LS = XZS_LS   ! no smooth orography in this case
ELSE IF (HFILE=='CHEM') THEN
  ALLOCATE (XZS_SV_LS(IIU,IJU))
  CALL ARRAY_1D_TO_2D (INO,ZOUT,IIU,IJU,XZS_SV_LS)
END IF
DEALLOCATE (ZOUT)
!
!---------------------------------------------------------------------------------------
!* 2.3 bis Read and interpol Sea Wave significant height
!---------------------------------------------------------------------------------------
WRITE (ILUOUT0,'(A)') ' | Searching sea wave significant height'
SELECT CASE (IMODEL)
  CASE(0) ! ECMWF
    ALLOCATE (XZWS(IIU,IJU))
    GFIND=.FALSE.
    !    
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=140229)
    IF(INUM < 0) THEN
      CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=229)
      !
      IF(INUM < 0) THEN
        WRITE (YVAL,'( E15.8 )') XZWS_DEFAULT
        WRITE (ILUOUT0,'(A)')' | !!! WARNING !!! Sea wave height is missing in '// &
               'the GRIB file - the default value of '//TRIM(YVAL)//' meters is used'
        XZWS = XZWS_DEFAULT
      ELSE
        GFIND=.TRUE.
      END IF
    ELSE
      GFIND=.TRUE. 
    END IF
  !
  IF (GFIND) THEN
    !!!!!!!!!!! Faire en sorte de le faire que pour le CASE(0)
    ! Sea wave significant height disponible uniquement pour ECMWF
    ! recuperation du tableau de valeurs
    CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
    ALLOCATE(IINLO(INJ))
    CALL COORDINATE_CONVERSION(IMODEL,IGRIB(INUM),IIU,IJU,ZLONOUT,ZLATOUT,&
          ZXOUT,ZYOUT,INI,ZPARAM,IINLO)
    ALLOCATE(ZVALUE(ISIZE))
    CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
    ! Change 9999 value to -1
    WHERE(ZVALUE.EQ.9999.) ZVALUE=0.
    ALLOCATE(ZOUT(INO))
    CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
              ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
    DEALLOCATE(IINLO)
    DEALLOCATE(ZVALUE)
    ! Stores the field in a 2 dimension array
    CALL ARRAY_1D_TO_2D (INO,ZOUT,IIU,IJU,XZWS)
    DEALLOCATE (ZOUT)
  END IF
END SELECT
!
!---------------------------------------------------------------------------------------
!* 2.4 Interpolation surface pressure
!---------------------------------------------------------------------------------------
!
!*  2.4.1  Read pressure
!   
WRITE (ILUOUT0,'(A)') ' | Searching pressure'

SELECT CASE (IMODEL)
  CASE(0) ! ECMWF
     CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=152)
     IF( INUM < 0 ) THEN
         WRITE (ILUOUT0,'(A)') ' | Logarithm of surface pressure is missing. It is then supposed that'
         WRITE (ILUOUT0,'(A)') ' | this ECMWF file has atmospheric fields on pressure levels (e.g. ERA5)'
         CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=134)
         IMODEL = 11
     END IF
  CASE(1,2,3,4,5) ! arpege mocage aladin et aladin reunion
      CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=1)
  CASE(6,7) ! NEW AROME,ARPEGE
     CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=3,KNUMBER=0)
  CASE(10) ! NCEP
      CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=134)
END SELECT
IF( INUM < 0 ) call Print_msg( NVERB_FATAL, 'IO', 'READ_ALL_DATA_GRIB_CASE', 'surface pressure is missing' )
! recuperation du tableau de valeurs
CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
ALLOCATE(ZVALUE(ISIZE))
CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
! determination des tableaux ZPS_G et ZLNPS_G
SELECT CASE (IMODEL)
  CASE(0,6,7) ! ECMWF
    ALLOCATE (ZPS_G  (ISIZE))
    ALLOCATE (ZLNPS_G(ISIZE))
    ZLNPS_G(:) =     ZVALUE(1:ISIZE)
    ZPS_G  (:) = EXP(ZVALUE(1:ISIZE))
  CASE(1,2,3,4,5,10,11) ! arpege mocage aladin aladin-reunion NCEP ERA5
    ALLOCATE (ZPS_G  (ISIZE))
    ALLOCATE (ZLNPS_G(ISIZE))
    ZPS_G  (:) =     ZVALUE(1:ISIZE)
    ZLNPS_G(:) = LOG(ZVALUE(1:ISIZE)) 
END SELECT
DEALLOCATE (ZVALUE)
!
!*  2.4.2  Removes pressure points not on orography points
!       (if pressure is on a regular grid)
CALL GRIB_GET(IGRIB(INUM),'typeOfGrid',HGRID)
CALL GRIB_GET(IGRIB(INUM_ZS),'typeOfGrid',HGRID_ZS)
CALL GRIB_GET(IGRIB(INUM),'Nj',INJ)
CALL GRIB_GET(IGRIB(INUM_ZS),'Nj',INJ_ZS)
!
IF ( HGRID(1:7)=='regular' .AND. HGRID_ZS(1:7)=='reduced' .AND.&
     INJ == INJ_ZS) THEN
  call Print_msg( NVERB_FATAL, 'IO', 'READ_ALL_DATA_GRIB_CASE', &
                  'HGRID(1:7)==regular .AND. HGRID_ZS(1:7)==reduced .AND. INJ == INJ_ZS' )
ELSE
   ALLOCATE(ZWORK_LNPS(SIZE(ZLNPS_G)))
   ZWORK_LNPS(:) = ZLNPS_G(:)
ENDIF
!
IF (HFILE(1:3)=='ATM') THEN
  ALLOCATE (XPS_LS(IIU,IJU))
ELSE IF (HFILE=='CHEM') THEN
  ALLOCATE (XPS_SV_LS(IIU,IJU))
END IF
!
ALLOCATE(IINLO(INJ))
CALL COORDINATE_CONVERSION(IMODEL,IGRIB(INUM),IIU,IJU,ZLONOUT,ZLATOUT,&
        ZXOUT,ZYOUT,INI,ZPARAM,IINLO)
ALLOCATE(ZOUT(INO))
CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI,&
         ZWORK_LNPS,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE. )
DEALLOCATE(ZWORK_LNPS)
DEALLOCATE(IINLO)
!
!* 2.4.3 conversion to surface pressure
!
ZOUT=EXP(ZOUT)
IF (HFILE(1:3)=='ATM') THEN
  CALL ARRAY_1D_TO_2D (INO,ZOUT,IIU,IJU,XPS_LS(:,:))
ELSE IF (HFILE=='CHEM') THEN
  CALL ARRAY_1D_TO_2D (INO,ZOUT,IIU,IJU,XPS_SV_LS(:,:))
END IF
!JUAN REALZ
CALL MPPDB_CHECK2D(XZS_LS,"XZS_LS",PRECISION)
!JUAN REALZ
DEALLOCATE (ZOUT)
DEALLOCATE (ZLNPS_G)
!
!---------------------------------------------------------------------------------------
!* 2.5 Interpolation temperature and humidity
!---------------------------------------------------------------------------------------
!
!* 2.5.1 Read T and Q on each level
!
WRITE (ILUOUT0,'(A)') ' | Reading T and Q fields'
!
IF (IMODEL/=10.AND.IMODEL/=11) THEN
  SELECT CASE (IMODEL)
    CASE(0) ! ECMWF
          ISTARTLEVEL=1
          IT=130
          IQ=133
    CASE(1,2,3,4,5) ! arpege mocage aladin et aladin reunion
          IT=11
          IQ=51
          ISTARTLEVEL=1
    CASE(6,7) !GRIB2 AROME AND ARPEGE
     IT=130
     IQ=133
     ISTARTLEVEL=1
  END SELECT
        
  CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IT,KLEV1=ISTARTLEVEL)
  IF(INUM < 0) THEN 
    ISTARTLEVEL=0
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IT,KLEV1=ISTARTLEVEL)
  ENDIF
  IF(INUM < 0) call Print_msg( NVERB_FATAL, 'IO', 'READ_ALL_DATA_GRIB_CASE', 'air temperature is missing' )
  CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IQ,KLEV1=ISTARTLEVEL)
  IF(INUM < 0) call Print_msg( NVERB_FATAL, 'IO', 'READ_ALL_DATA_GRIB_CASE', 'atmospheric specific humidity is missing' )
ELSEIF (IMODEL==10) THEN ! NCEP
  ISTARTLEVEL=1000
  IT=130
  IQ=157
  CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IT,KLEV1=ISTARTLEVEL)
  IF(INUM < 0) call Print_msg( NVERB_FATAL, 'IO', 'READ_ALL_DATA_GRIB_CASE', 'air temperature is missing' )
  CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IQ,KLEV1=ISTARTLEVEL)
  IF(INUM < 0) call Print_msg( NVERB_FATAL, 'IO', 'READ_ALL_DATA_GRIB_CASE', 'atmospheric relative humidity is missing' )
ELSE ! ERA5
  ISTARTLEVEL=1000
  IT=130
  IQ=133
  CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IT,KLEV1=ISTARTLEVEL)
  IF(INUM < 0) call Print_msg( NVERB_FATAL, 'IO', 'READ_ALL_DATA_GRIB_CASE', 'air temperature is missing' )
  CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IQ,KLEV1=ISTARTLEVEL)
  IF(INUM < 0) call Print_msg( NVERB_FATAL, 'IO', 'READ_ALL_DATA_GRIB_CASE', 'atmospheric specific humidity is missing' )
ENDIF
!
IF (IMODEL/=10.AND.IMODEL/=11) THEN ! others than NCEP AND ERA5
  CALL GRIB_GET(IGRIB(INUM),'NV',INLEVEL)
  INLEVEL = NINT(INLEVEL / 2.) - 1
  CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
ELSE
  IF (IMODEL==10) THEN
    INLEVEL=JP_GFS
  ELSE
    INLEVEL=JP_ERA
  END IF
END IF
!
ALLOCATE (ZT_G(ISIZE,INLEVEL))
ALLOCATE (ZQ_G(ISIZE,INLEVEL))
!
IF (IMODEL/=10.AND.IMODEL/=11) THEN ! others than NCEP AND ERA5
  DO JLOOP1=1, INLEVEL
    ILEV1 = JLOOP1-1+ISTARTLEVEL
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IQ,KLEV1=ILEV1)
    IF (INUM< 0) THEN
      WRITE(YMSG,*) 'atmospheric humidity level ',JLOOP1,' is missing'
      CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
    END IF
    CALL GRIB_GET(IGRIB(INUM),'values',ZQ_G(:,INLEVEL-JLOOP1+1))
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IT,KLEV1=ILEV1)
    IF (INUM< 0) THEN
      WRITE(YMSG,*) 'atmospheric temperature level ',JLOOP1,' is missing'
      CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
    END IF
    CALL GRIB_GET(IGRIB(INUM),'values',ZT_G(:,INLEVEL-JLOOP1+1))
    CALL GRIB_GET(IGRIB(INUM),'Nj',INJ,IRET_GRIB)
  END DO
ELSEIF (IMODEL==10) THEN ! NCEP
  DO JLOOP1=1, INLEVEL
    ILEV1 = IP_GFS(JLOOP1)
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IQ,KLEV1=ILEV1)
    IF (INUM< 0) THEN
      WRITE(YMSG,*) 'atmospheric humidity level ',JLOOP1,' is missing'
      CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
    END IF
    CALL GRIB_GET(IGRIB(INUM),'values',ZQ_G(:,JLOOP1),IRET_GRIB)
    CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=0,KNUMBER=0,KLEV1=ILEV1,KTFFS=100)
    IF (INUM< 0) THEN
      WRITE(YMSG,*) 'atmospheric temperature level ',JLOOP1,' is missing'
      CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
    END IF
    CALL GRIB_GET(IGRIB(INUM),'values',ZT_G(:,JLOOP1),IRET_GRIB)
    CALL GRIB_GET(IGRIB(INUM),'Nj',INJ,IRET_GRIB)
  END DO
ELSE ! ERA5
  DO JLOOP1=1, INLEVEL
    ILEV1 = IP_ERA(JLOOP1)
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IQ,KLEV1=ILEV1)
    IF (INUM< 0) THEN
      WRITE(YMSG,*) 'atmospheric humidity level ',JLOOP1,' is missing'
      CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
    END IF
    CALL GRIB_GET(IGRIB(INUM),'values',ZQ_G(:,JLOOP1),IRET_GRIB)
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IT,KLEV1=ILEV1)
    IF (INUM< 0) THEN
      WRITE(YMSG,*) 'atmospheric temperature level ',JLOOP1,' is missing'
      CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
    END IF
    CALL GRIB_GET(IGRIB(INUM),'values',ZT_G(:,JLOOP1),IRET_GRIB)
    CALL GRIB_GET(IGRIB(INUM),'Nj',INJ,IRET_GRIB)
  END DO
END IF

ALLOCATE(IINLO(INJ))
CALL COORDINATE_CONVERSION(IMODEL,IGRIB(INUM),IIU,IJU,ZLONOUT,ZLATOUT,&
        ZXOUT,ZYOUT,INI,ZPARAM,IINLO)
!
!*  2.5.2  Load level definition parameters A and B
!
IF (IMODEL/=10.AND.IMODEL/=11) THEN ! others than NCEP AND ERA5

  IF (HFILE(1:3)=='ATM') THEN
    XP00_LS = 101325.
  ELSE IF (HFILE=='CHEM') THEN
    XP00_SV_LS = 101325.
  END IF
!
  IF (INLEVEL > 0) THEN
    IF (HFILE(1:3)=='ATM') THEN
      ALLOCATE (XA_LS(INLEVEL))
      ALLOCATE (XB_LS(INLEVEL))
    ELSE IF (HFILE=='CHEM') THEN
      ALLOCATE (XA_SV_LS(INLEVEL))
      ALLOCATE (XB_SV_LS(INLEVEL))
    END IF
!
    CALL GRIB_GET(IGRIB(INUM),'PVPresent',IPVPRESENT)
    IF (IPVPRESENT==1) THEN
       CALL GRIB_GET_SIZE(IGRIB(INUM),'pv',IPV)
       ALLOCATE(ZPV(IPV))
       CALL GRIB_GET(IGRIB(INUM),'pv',ZPV)
    ELSE
      CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE','there is no PV value in this message')
    ENDIF
    SELECT CASE (IMODEL)
      CASE (0,3,4,6,7)
        DO JLOOP1 = 1, INLEVEL
          XA_LS(1 + INLEVEL - JLOOP1) = ZPV(1+JLOOP1) / XP00_LS
          XB_LS(1 + INLEVEL - JLOOP1) = ZPV(2+INLEVEL+JLOOP1)
        END DO
      CASE (1,2)
        JLOOP2 = 2
        DO JLOOP1 = 1, INLEVEL
          JLOOP2 = JLOOP2 + 1
          XA_LS(1 + INLEVEL - JLOOP1) = ZPV(JLOOP2)
          JLOOP2 = JLOOP2 + 1
          XB_LS(1 + INLEVEL - JLOOP1) = ZPV(JLOOP2)
        END DO
      CASE (5)
        DO JLOOP1 = 1, INLEVEL
          IF (HFILE(1:3)=='ATM') THEN
            XA_LS(1 + INLEVEL - JLOOP1) = ZPV(1+        JLOOP1) / XP00_LS**2 
            XB_LS(1 + INLEVEL - JLOOP1) = ZPV(2+INLEVEL+JLOOP1)
          ELSE IF (HFILE=='CHEM') THEN
            XA_SV_LS(1 + INLEVEL - JLOOP1) = ZPV(1+        JLOOP1) / XP00_LS**2 
            XB_SV_LS(1 + INLEVEL - JLOOP1) = ZPV(2+INLEVEL+JLOOP1)
          END IF
        END DO
    END SELECT
  ELSE
    CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE','level definition section is missing')
  END IF
ELSE
  ALLOCATE (XA_LS(INLEVEL))
  ALLOCATE (XB_LS(0))
  IF (IMODEL==10) THEN
    XA_LS = 100.*IP_GFS
  ELSE
    XA_LS = 100.*IP_ERA
  END IF
END IF
!
!*  2.5.3  Compute atmospheric pressure on grib grid
!
WRITE (ILUOUT0,'(A)') ' | Atmospheric pressure on Grib grid is being computed'

ALLOCATE (ZPF_G(INI,INLEVEL))
IF (IMODEL/=10.AND.IMODEL/=11) THEN ! others than NCEP and ERA5
  IF (HFILE(1:3)=='ATM') THEN
    ZPF_G(:,:) = SPREAD(XA_LS,1,INI)*XP00_LS + &
    SPREAD(XB_LS,1,INI)*SPREAD(ZPS_G(1:INI),2,INLEVEL)
  ELSE IF (HFILE=='CHEM') THEN
    ZPF_G(:,:) = SPREAD(XA_SV_LS,1,INI)*XP00_SV_LS + &
    SPREAD(XB_SV_LS,1,INI)*SPREAD(ZPS_G(1:INI),2,INLEVEL)
  END IF
ELSE
  IF (IMODEL==10) THEN
    ZPF_G(:,:) = 100.*SPREAD(IP_GFS,1,INI)
  ELSE
    ZPF_G(:,:) = 100.*SPREAD(IP_ERA,1,INI)
  END IF
END IF
DEALLOCATE (ZPS_G)
!
ALLOCATE (ZEXNF_G(INI,INLEVEL))
ZEXNF_G(:,:) = (ZPF_G(:,:)/XP00)**(XRD/XCPD)
!
ALLOCATE (ZEXNM_G(INI,INLEVEL))
ZEXNM_G(:,1:INLEVEL-1) = (ZEXNF_G(:,1:INLEVEL-1)-ZEXNF_G(:,2:INLEVEL)) / &
                     (LOG(ZEXNF_G(:,1:INLEVEL-1))-LOG(ZEXNF_G(:,2:INLEVEL)))
ZEXNM_G(:,INLEVEL) = (ZPF_G(:,INLEVEL)/2./XP00)**(XRD/XCPD)
!
IF (IMODEL==10.OR.IMODEL==11) ZEXNM_G(:,:)=ZEXNF_G(:,:) ! for GFS and ERA5 on pressure levels
!
DEALLOCATE (ZEXNF_G)
DEALLOCATE (ZPF_G)
!
ALLOCATE (ZPM_G(INI,INLEVEL))
ZPM_G(:,:) = XP00*(ZEXNM_G(:,:))**(XCPD/XRD)
!
!*  2.5.4  Compute the relative humidity and the interpolate Thetav, P, Q and H
!
IF (IMODEL==1) THEN
  ! search cloud_water in Arome case (forecast)
  ISTARTLEVEL = 1
  IPAR=246
  CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR,KLEV1=ISTARTLEVEL)
  IF (INUM < 0) THEN
    ISTARTLEVEL = 0
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR,KLEV1=ISTARTLEVEL)
  END IF
  IF (INUM > 0) THEN
    WRITE (ILUOUT0,'(A)') ' | Grib file from French Weather Service - Arome model (forecast)'
    LCPL_AROME=.TRUE.
    NRR=6
  END IF
  ! search also turbulent kinetic energy 
  ISTARTLEVEL = 1
  IPAR=251
  CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR,KLEV1=ISTARTLEVEL)
  IF (INUM < 0) THEN
    ISTARTLEVEL = 0
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR,KLEV1=ISTARTLEVEL)
  END IF
  IF (INUM > 0) CTURB='TKEL'
END IF

IF (IMODEL==6) THEN ! GRIB2 AROME
! search cloud_water in Arome case (forecast)
  ISTARTLEVEL = 1
  CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=6,KNUMBER=6,KLEV1=ISTARTLEVEL)
  IF (INUM < 0) THEN
    ISTARTLEVEL = 0
    CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=6,KNUMBER=6,KLEV1=ISTARTLEVEL)
  END IF
  IF (INUM < 0) THEN
    ISTARTLEVEL = 1
    CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=1,KNUMBER=0,KLEV1=ISTARTLEVEL)
  END IF
  IF (INUM > 0) THEN
    WRITE (ILUOUT0,'(A)') ' | Grib file from French Weather Service - Arome model (forecast)'
    LCPL_AROME=.TRUE.
    NRR=6
  END IF
  ! search also turbulent kinetic energy 
  ISTARTLEVEL = 1
  CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=19,KNUMBER=11,KLEV1=ISTARTLEVEL)
  IF (INUM < 0) THEN
    ISTARTLEVEL = 0
    CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=19,KNUMBER=11,KLEV1=ISTARTLEVEL)
  END IF
  IF (INUM > 0) CTURB='TKEL'
END IF
!
!
WRITE (ILUOUT0,'(A)') ' | Computing relative humidity on each level'
ALLOCATE (ZH_G(INI))
ALLOCATE (ZH_LS(IIU,IJU,INLEVEL))
IF (HFILE(1:3)=='ATM') THEN
  ALLOCATE (XT_LS(IIU,IJU,INLEVEL))
  ALLOCATE (XQ_LS(IIU,IJU,INLEVEL,NRR)) ; XQ_LS=0.
ELSE IF (HFILE=='CHEM') THEN
  ALLOCATE (XT_SV_LS(IIU,IJU,INLEVEL))
  ALLOCATE (XQ_SV_LS(IIU,IJU,INLEVEL,1))
END IF
IF (CTURB=='TKEL') THEN
  IF (ALLOCATED(XTKE_LS)) DEALLOCATE(XTKE_LS)
  ALLOCATE(XTKE_LS(IIU,IJU,INLEVEL)) ; XTKE_LS=0.
  CALL INI_CTURB
ELSE
  IF (ALLOCATED(XTKE_LS)) DEALLOCATE(XTKE_LS)
  ALLOCATE(XTKE_LS(0,0,0))
END IF
ALLOCATE (ZTHV_LS(IIU,IJU,INLEVEL))
ALLOCATE (ZTHV_G(INI))
ALLOCATE (ZRV_G(INI))
ALLOCATE (ZOUT(INO))
IF (IMODEL/=10) THEN ! others than NCEP
  DO JLOOP1=1, INLEVEL
    !
    ! Compute Theta V and relative humidity on grib grid
    !
    !   (1/rv) = (1/q)  -  1
    !   Thetav = T . (P0/P)^(Rd/Cpd) . ( (1 + (Rv/Rd).rv) / (1 + rv) )
    !   Hu = P / ( ( (Rd/Rv) . ((1/rv) - 1) + 1) . Es(T) )
    !
    ZRV_G(:) = 1. / (1./MAX(ZQ_G(:,JLOOP1),1.E-12) - 1.)
    !
    ZTHV_G(:)=ZT_G(:,JLOOP1) * ((XP00/ZPM_G(:,JLOOP1))**(XRD/XCPD)) * &
                             ((1. + XRV*ZRV_G(:)/XRD) / (1. + ZRV_G(:)))
    !
    ZH_G(1:INI) = 100. * ZPM_G(:,JLOOP1) / ( (XRD/XRV)*(1./ZRV_G(:)+1.)*SM_FOES(ZT_G(:,JLOOP1)) )
    ZH_G(:) = MAX(MIN(ZH_G(:),100.),0.)
    !
    ! Interpolation : H           
    CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
        ZH_G,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
    CALL ARRAY_1D_TO_2D (INO,ZOUT,IIU,IJU,ZH_LS(:,:,JLOOP1))
    ZH_LS(:,:,JLOOP1) = MAX(MIN(ZH_LS(:,:,JLOOP1),100.),0.)
    !
    ! interpolation : Theta V
    CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
        ZTHV_G,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
    CALL ARRAY_1D_TO_2D (INO,ZOUT,IIU,IJU,ZTHV_LS(:,:,JLOOP1))
    !
  END DO
ELSE !GFS and ERA5 on pressure levels
  DO JLOOP1=1, INLEVEL
    ZH_G(:)  =ZQ_G(:,JLOOP1)
    ZRV_G(:) = (XRD/XRV)*SM_FOES(ZT_G(:,JLOOP1))*0.01*ZH_G(:) &
                        /(ZPM_G(:,JLOOP1) -SM_FOES(ZT_G(:,JLOOP1))*0.01*ZH_G(:))
    ZTHV_G(:)=ZT_G(:,JLOOP1) * ((XP00/ZPM_G(:,JLOOP1))**(XRD/XCPD)) * &
                               ((1. + XRV*ZRV_G(:)/XRD) / (1. + ZRV_G(:)))
    !
    ! Interpolation : H           
    CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
          ZH_G,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
    CALL ARRAY_1D_TO_2D (INO,ZOUT,IIU,IJU,ZH_LS(:,:,JLOOP1))
    ZH_LS(:,:,JLOOP1) = MAX(MIN(ZH_LS(:,:,JLOOP1),100.),0.)
    !
    ! interpolation : Theta V
    CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
          ZTHV_G,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
    CALL ARRAY_1D_TO_2D (INO,ZOUT,IIU,IJU,ZTHV_LS(:,:,JLOOP1))
    !
  END DO
END IF
  
DEALLOCATE (ZOUT)


!---------------------------------------------------------------------------------------
!* 2.5.4.2 Read and interpol geopotential height for interpolation on isobaric surface Grid of NCEP 
!---------------------------------------------------------------------------------------
!
ALLOCATE (ZGH_G(ISIZE,INLEVEL))
!
IF (IMODEL==10.OR.IMODEL==11) THEN !NCEP or ERA5 with pressure grid only
 DO JLOOP1=1, INLEVEL
  IF (IMODEL==10) THEN
    ILEV1 = IP_GFS(JLOOP1)
    CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=3,KNUMBER=5,KLEV1=ILEV1)
  ELSE
    ILEV1 = IP_ERA(JLOOP1)
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=129,KLEV1=ILEV1)
  END IF
    IF (INUM< 0) THEN
    !callabortstop
      WRITE(YMSG,*) 'Geopotential height level ',JLOOP1,' is missing'
      CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_ALL_DATA_GRIB_CASE',YMSG)
    END IF
  !
  CALL GRIB_GET(IGRIB(INUM),'values',ZGH_G(:,JLOOP1),IRET_GRIB)
  CALL GRIB_GET(IGRIB(INUM),'Nj',INJ,IRET_GRIB)
  !
  IF (IMODEL/=10) THEN
    ! Data given in archives are multiplied by the gravity acceleration
    ZGH_G(:,JLOOP1) = ZGH_G(:,JLOOP1) / XG
  END IF
  !
  END DO
 !
 CALL COORDINATE_CONVERSION(IMODEL,IGRIB(INUM_ZS),IIU,IJU,ZLONOUT,ZLATOUT,&
          ZXOUT,ZYOUT,INI,ZPARAM,IINLO)
 !
 ALLOCATE(ZOUT(INO))
 ALLOCATE(XGH_LS(IIU,IJU,INLEVEL))
 !
 DO JLOOP1=1, INLEVEL
  CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
              ZGH_G(:,JLOOP1),INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
  CALL ARRAY_1D_TO_2D (INO,ZOUT,IIU,IJU,XGH_LS(:,:,JLOOP1))
 END DO
 DEALLOCATE(ZOUT)
END IF
!
!*  2.5.5  Compute atmospheric pressure on MESO-NH grid
!
WRITE (ILUOUT0,'(A)') ' | Atmospheric pressure on the Meso-NH grid is being computed'
ALLOCATE (ZPF_LS(IIU,IJU,INLEVEL))
IF (IMODEL/=10.AND.IMODEL/=11) THEN ! others than NCEP and ERA5
  IF (HFILE(1:3)=='ATM') THEN
    ZPF_LS(:,:,:) = SPREAD(SPREAD(XA_LS,1,IIU),2,IJU)*XP00_LS + &
    SPREAD(SPREAD(XB_LS,1,IIU),2,IJU)*SPREAD(XPS_LS,3,INLEVEL)
  ELSE IF (HFILE=='CHEM') THEN
    ZPF_LS(:,:,:) = SPREAD(SPREAD(XA_SV_LS,1,IIU),2,IJU)*XP00_LS + &
    SPREAD(SPREAD(XB_SV_LS,1,IIU),2,IJU)*SPREAD(XPS_SV_LS,3,INLEVEL)
  END IF
ELSE
  IF(IMODEL==10) THEN
    ZPF_LS(:,:,:) = 100.*SPREAD(SPREAD(IP_GFS,1,IIU),2,IJU)
  ELSE
    ZPF_LS(:,:,:) = 100.*SPREAD(SPREAD(IP_ERA,1,IIU),2,IJU)
  END IF
END IF  
!
ALLOCATE (ZEXNF_LS(IIU,IJU,INLEVEL))
ZEXNF_LS(:,:,:) = (ZPF_LS(:,:,:)/XP00)**(XRD/XCPD)
!
ALLOCATE (ZEXNM_LS(IIU,IJU,INLEVEL))
ZEXNM_LS(:,:,1:INLEVEL-1) = (ZEXNF_LS(:,:,1:INLEVEL-1)-ZEXNF_LS(:,:,2:INLEVEL)) / &
                     (LOG(ZEXNF_LS(:,:,1:INLEVEL-1))-LOG(ZEXNF_LS(:,:,2:INLEVEL)))
ZEXNM_LS(:,:,INLEVEL) = (ZPF_LS(:,:,INLEVEL)/2./XP00)**(XRD/XCPD)
!
IF (IMODEL==10.OR.IMODEL==11) ZEXNM_LS(:,:,:)=ZEXNF_LS(:,:,:) ! for GFS and ERA5 on pressure levels
!
DEALLOCATE (ZEXNF_LS)
DEALLOCATE (ZPF_LS)
!
ALLOCATE (ZPM_LS(IIU,IJU,INLEVEL))
ZPM_LS(:,:,:) = XP00*(ZEXNM_LS(:,:,:))**(XCPD/XRD)
!
!*  2.5.6 Compute the vapor mixing ratio and the final specific humdity
!
!  The vapor mixing ratio is calculated by an interating process on rv and
!  Thetav. Have a look to MODE_THERMO for further informations.
ALLOCATE (ZR_DUM(IIU,IJU,INLEVEL,1))
ALLOCATE (ZRV_LS(IIU,IJU,INLEVEL))
ALLOCATE (ZTEV_LS(IIU,IJU,INLEVEL))
ZTEV_LS(:,:,:) = ZTHV_LS(:,:,:) * ZEXNM_LS(:,:,:)
ZRV_LS(:,:,:) = SM_PMR_HU(ZPM_LS(:,:,:),                &
ZTEV_LS(:,:,:),ZH_LS(:,:,:),ZR_DUM(:,:,:,:),KITERMAX=100)
IF (HFILE(1:3)=='ATM') THEN
  XQ_LS(:,:,:,1) = ZRV_LS(:,:,:) / (1. + ZRV_LS(:,:,:))
ELSE IF (HFILE=='CHEM') THEN
  XQ_SV_LS(:,:,:,1) = ZRV_LS(:,:,:) / (1. + ZRV_LS(:,:,:))
ENDIF
!JUAN
CALL  MPPDB_CHECK3D(XQ_LS(:,:,:,1),"XQ_LS",PRECISION)
!JUAN
DEALLOCATE (ZTEV_LS)
DEALLOCATE (ZH_LS)
DEALLOCATE (ZR_DUM)
!
!*  2.5.7 Compute T from the interpolated Theta V
!
!   T = Thetav . (P/P0)^(Rd/Cpd) . ((1 + rv) / (1 + (Rv/Rd).rv))
!!
IF (HFILE(1:3)=='ATM') THEN
  XT_LS(:,:,:) = ZTHV_LS(:,:,:) * ZEXNM_LS(:,:,:) &
                              * ((1.+ZRV_LS(:,:,:))/(1.+(XRV/XRD)*ZRV_LS(:,:,:)))
 CALL MPPDB_CHECK3D(XT_LS,"XT_LS",PRECISION)
ELSE IF (HFILE=='CHEM') THEN
  XT_SV_LS(:,:,:) = ZTHV_LS(:,:,:) * ZEXNM_LS(:,:,:) &
                              * ((1.+ZRV_LS(:,:,:))/(1.+(XRV/XRD)*ZRV_LS(:,:,:)))
 CALL MPPDB_CHECK3D(XT_SV_LS,"XT_SV_LS",PRECISION)
ENDIF
!
DEALLOCATE (ZRV_LS)
DEALLOCATE (ZTHV_LS)
DEALLOCATE (ZPM_LS)
DEALLOCATE (ZEXNM_LS)
!
!*  2.5.8 Read the other specific ratios (if Arome model)
!
IF (NRR >1) THEN
  IF (IMODEL==1) THEN
    WRITE (ILUOUT0,'(A)') ' | Reading Q fields (except humidity)'
    DO JLOOP2=1,NRR-1
      IPAR=246+JLOOP2-1
      DO JLOOP1=1, INLEVEL
        ILEV1 = JLOOP1-1+ISTARTLEVEL
        CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR,KLEV1=ILEV1)

        IF (INUM < 0) THEN
          WRITE(YMSG,*) 'Specific ratio ',IPAR,' at level ',JLOOP1,' is missing'
          CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
        END IF
        CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
        ALLOCATE(ZVALUE(ISIZE))
        CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
        ALLOCATE(ZOUT(INO))
        CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
                    ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
        CALL ARRAY_1D_TO_2D(INO,ZOUT,IIU,IJU,XQ_LS(:,:,INLEVEL-JLOOP1+1,1+JLOOP2))
        DEALLOCATE (ZVALUE)
        DEALLOCATE (ZOUT)
      END DO
    END DO 
  ELSE ! GRIB2 AROME IMODEL =6
    WRITE (ILUOUT0,'(A)') ' | Reading Q fields (except humidity)'
    DO JLOOP1=1, INLEVEL
      ILEV1 = JLOOP1-1+ISTARTLEVEL
      CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=1,KNUMBER=83,KLEV1=ILEV1)

      IF (INUM < 0) THEN
        WRITE(YMSG,*) 'Specific ratio ',IPAR,' at level ',JLOOP1,' is missing'
        CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
      END IF
      CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
      ALLOCATE(ZVALUE(ISIZE))
      CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
      ALLOCATE(ZOUT(INO))
      CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
        ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
      CALL ARRAY_1D_TO_2D(INO,ZOUT,IIU,IJU,XQ_LS(:,:,INLEVEL-JLOOP1+1,2))
      DEALLOCATE (ZVALUE)
      DEALLOCATE (ZOUT)
    END DO

    DO JLOOP1=1, INLEVEL
      ILEV1 = JLOOP1-1+ISTARTLEVEL
      CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=1,KNUMBER=85,KLEV1=ILEV1)

      IF (INUM < 0) THEN
        WRITE(YMSG,*) 'Specific ratio for rain at level ',JLOOP1,' is missing'
        CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
      END IF
      CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
      ALLOCATE(ZVALUE(ISIZE))
      CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
      ALLOCATE(ZOUT(INO))
      CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
        ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
      CALL ARRAY_1D_TO_2D(INO,ZOUT,IIU,IJU,XQ_LS(:,:,INLEVEL-JLOOP1+1,3))
      DEALLOCATE (ZVALUE)
      DEALLOCATE (ZOUT)
    END DO  


    DO JLOOP1=1, INLEVEL
      ILEV1 = JLOOP1-1+ISTARTLEVEL
      CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=1,KNUMBER=84,KLEV1=ILEV1)
      IF (INUM < 0) THEN
        WRITE(YMSG,*) 'Specific ratio for ICE at level ',JLOOP1,' is missing'
        CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
      END IF
      CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
      ALLOCATE(ZVALUE(ISIZE))
      CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
      ALLOCATE(ZOUT(INO))
      CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
        ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
      CALL ARRAY_1D_TO_2D(INO,ZOUT,IIU,IJU,XQ_LS(:,:,INLEVEL-JLOOP1+1,4))
      DEALLOCATE (ZVALUE)
      DEALLOCATE (ZOUT)
    END DO


    DO JLOOP1=1, INLEVEL
      ILEV1 = JLOOP1-1+ISTARTLEVEL
      CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=1,KNUMBER=86,KLEV1=ILEV1)
      IF (INUM < 0) THEN
        WRITE(YMSG,*) 'Specific ratio ',IPAR,' at level ',JLOOP1,' is missing'
        CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
      END IF
      CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
      ALLOCATE(ZVALUE(ISIZE))
      CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
      ALLOCATE(ZOUT(INO))
      CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
        ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
      CALL ARRAY_1D_TO_2D(INO,ZOUT,IIU,IJU,XQ_LS(:,:,INLEVEL-JLOOP1+1,5))
      DEALLOCATE (ZVALUE)
      DEALLOCATE (ZOUT)
    END DO


    DO JLOOP1=1, INLEVEL
      ILEV1 = JLOOP1-1+ISTARTLEVEL
      CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=1,KNUMBER=201,KLEV1=ILEV1)
      IF (INUM < 0) THEN
        WRITE(YMSG,*) 'Specific ratio ',IPAR,' at level ',JLOOP1,' is missing'
        CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
      END IF
      CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
      ALLOCATE(ZVALUE(ISIZE))
      CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
      ALLOCATE(ZOUT(INO))
      CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
        ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
      CALL ARRAY_1D_TO_2D(INO,ZOUT,IIU,IJU,XQ_LS(:,:,INLEVEL-JLOOP1+1,6))
      DEALLOCATE (ZVALUE)
      DEALLOCATE (ZOUT)
    END DO
  END IF
END IF
!
IF (CTURB=='TKEL') THEN
  WRITE (ILUOUT0,'(A)') ' | Reading TKE field'
  DO JLOOP1=1, INLEVEL
    ILEV1 = JLOOP1-1+ISTARTLEVEL
    IF (IMODEL==1) THEN
      CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=251,KLEV1=ILEV1)
    ELSE ! case 6 new arome
      CALL SEARCH_FIELD(IGRIB,INUM,KDIS=0,KCAT=19,KNUMBER=11,KLEV1=ILEV1)
    END IF
    IF (INUM <  0) THEN
      WRITE(YMSG,*) 'TKE at level ',JLOOP1,' is missing'
      CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
    END IF
    CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
    ALLOCATE(ZVALUE(ISIZE))
    CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
    ALLOCATE(ZOUT(INO))
    CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
        ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE.)
    CALL ARRAY_1D_TO_2D(INO,ZOUT,IIU,IJU,XTKE_LS(:,:,INLEVEL-JLOOP1+1))
    DEALLOCATE (ZVALUE)
    DEALLOCATE (ZOUT)
  END DO
END IF
DEALLOCATE(IINLO)
!
!---------------------------------------------------------------------------------------
!* 2.6 Interpolation of MOCAGE variable
!---------------------------------------------------------------------------------------

IF (IMODEL==5) THEN 
  LUSECHEM = .TRUE.
  IF (LORILAM) THEN
    CORGANIC = "MPMPO"
    LVARSIGI = .TRUE.
    LVARSIGJ = .TRUE.
  END IF
  ! initialise NSV_* variables
  CALL INI_NSV(IMI)
  IF( HFILE=='ATM0' ) THEN
    ALLOCATE (XSV_LS(IIU,IJU,INLEVEL,NSV))
  ELSE IF (HFILE=='CHEM' ) THEN
    DEALLOCATE(XSV_LS)
    ALLOCATE (XSV_LS(IIU,IJU,INLEVEL,NSV))
  ELSE
    CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE','Mocage model: Bad input argument in read_all_data_grib_case')
  END IF
  XSV_LS(:,:,:,:) = 0.
  ILEV1=-1
!
  WRITE (ILUOUT0,'(A,A4,A)') ' | Reading Mocage species (ppp) from ',HFILE,' file'
!
!*       2.6.1  read mocage species
!
! open input file
  CALL CH_OPEN_INPUT(YPRE_MOC, "MOC2MESONH", TZFILE, ILUOUT0, KVERB)
  ICHANNEL = TZFILE%NLU
!
! read number of mocage species to transfer into mesonh
  READ(ICHANNEL, *) IMOC
  IF (KVERB >= 5) WRITE(ILUOUT0,*) "number of mocage species to transfer into mesonh : ", IMOC
!
! read data input format
  READ(ICHANNEL,"(A)") YFORMAT
  YFORMAT=UPCASE(YFORMAT)
  IF (KVERB >= 5) WRITE(ILUOUT0,*) "input format is: ", YFORMAT
!
! allocate fields
  ALLOCATE(YMNHNAME(IMOC))
  ALLOCATE(INUMGRIB(IMOC))
!
! read variables names and Grib code
  IF (INDEX(YFORMAT,'A') < INDEX(YFORMAT,'I')) THEN
    DO JI = 1, IMOC
      READ(ICHANNEL,YFORMAT) YMNHNAME(JI), INUMGRIB(JI)
      WRITE(ILUOUT0,YFORMAT) YMNHNAME(JI), INUMGRIB(JI)
    END DO
  ELSE
    DO JI = 1, IMOC
      READ(ICHANNEL,YFORMAT) INUMGRIB(JI), YMNHNAME(JI)
      WRITE(ILUOUT0,YFORMAT) INUMGRIB(JI), YMNHNAME(JI)
    END DO
  ENDIF
  !
  ! close file
  CALL IO_File_close(TZFILE)
  TZFILE => NULL()
  !
  !*  2.6.2   exchange mocage values onto prognostic variables XSV_LS
  !
  IF (KVERB >= 10) WRITE(ILUOUT0,'(A,I4)') ' NEQ=',NEQ
  !
  DO JNREAL = 1, NEQ
    INACT = 0
    search_loop2 : DO JN = 1, IMOC
      IF (CNAMES(JNREAL) .EQ. YMNHNAME(JN)) THEN
        INACT = JN
        EXIT search_loop2
      END IF
    END DO search_loop2
    IF (INACT .NE. 0) THEN
      DO JLOOP1=1, INLEVEL
        ILEV1 = JLOOP1
        CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=INUMGRIB(JN),KLEV1=ILEV1)
        IF (INUM <  0) THEN
          WRITE(YMSG,*) 'Atmospheric ',INUMGRIB(JN),' grib chemical species level ',JLOOP1,' is missing'
          CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
        END IF
        CALL GRIB_GET(IGRIB(INUM),'Nj',INJ,IRET_GRIB)
        ALLOCATE(IINLO(INJ))
        CALL COORDINATE_CONVERSION(IMODEL,IGRIB(INUM),IIU,IJU,ZLONOUT,ZLATOUT,&
          ZXOUT,ZYOUT,INI,ZPARAM,IINLO)
        CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
        ALLOCATE(ZVALUE(ISIZE))
        CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
        ALLOCATE(ZOUT(INO))
        CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
           ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.TRUE. )
        CALL ARRAY_1D_TO_2D(INO,ZOUT,IIU,IJU, &
                            XSV_LS(:,:,INLEVEL-JLOOP1+1,JNREAL))
        DEALLOCATE (ZVALUE)
        DEALLOCATE (ZOUT)
        DEALLOCATE(IINLO)
      END DO
    END IF
  END DO
  XSV_LS(:,:,:,:) = MAX(XSV_LS(:,:,:,:),0.)
ELSE
  LORILAM  = .FALSE. 
  LUSECHEM = .FALSE. 
  ! initialise NSV_* variables
  CALL INI_NSV(1)
  IF (NSV > 0) THEN 
    ALLOCATE (XSV_LS(IIU,IJU,INLEVEL,NSV))
    XSV_LS(:,:,:,:) = 0.
  ELSE 
    ALLOCATE(XSV_LS(0,0,0,0))
  END IF
END IF
!
!---------------------------------------------------------------------------------------
!* 2.7 Search, read, interpolate and project winds
!---------------------------------------------------------------------------------------
!
! The way winds are processed depends upon the type of archive :
!
! -> ECMWF, NCEP
!   Winds are projected from a standard lat,lon grid to MesoNH grid. This correcponds to
!   a rotation of an angle :
!    Alpha = k.(L-L0) - Beta
!   k,L0 and Beta definiiton parameter of MesoNH grid
!   L longitude of the point of the tangent coordinate system
!
! -> Aladin
!   The grid used by Aladin files is the same than the one of MesoNH. !
! So we have 2 sets  of parameters :
!     k   L0   Beta      for MesoNH
!     k'  L0'  Beta'     for Aladin (Beta'=0 for Aladin)
!   We applied twice the formula seen for standard lat,lon grid and we get :
!    Alpha = k.(L-L0) - Beta - k'.(L-L0')
!
! -> Arpege
! Arpege winds are given on the tangent coordinate system of the stretched grid.
! Therefore they have first to be projected on a standard lat,lon grid. This is done
! before the interpolation. The projection formulas were given by Meteo France.
! After this projection, the file is simil
!
IF (HFILE(1:3)=='ATM') THEN
ISTARTLEVEL = 1
ALLOCATE (XU_LS(IIU,IJU,INLEVEL))
ALLOCATE (XV_LS(IIU,IJU,INLEVEL))
ALLOCATE (ZTU_LS(INO))
ALLOCATE (ZTV_LS(INO))
!
SELECT CASE (IMODEL)
  CASE (0,6,7)
    IPAR = 131
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR,KLEV1=ISTARTLEVEL)
    IF (INUM< 0) THEN
      ISTARTLEVEL = 0
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR,KLEV1=ISTARTLEVEL)
    END IF
  CASE (1,2,3)
    IPAR  = 33
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR,KLEV1=ISTARTLEVEL)
    IF (INUM < 0) THEN
      ISTARTLEVEL = 0
      CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR,KLEV1=ISTARTLEVEL)
    END IF
  CASE (10,11)
    IPAR = 131
    ISTARTLEVEL = 1
END SELECT

DO JLOOP1 = ISTARTLEVEL, ISTARTLEVEL+INLEVEL-1
  IF (IMODEL/=10.AND.IMODEL/=11) THEN ! others than NCEP and ERA5
    ILEV1 = JLOOP1
  ELSE
    IF(IMODEL==10) THEN
      ILEV1 = IP_GFS(INLEVEL+ISTARTLEVEL-JLOOP1)
    ELSE
      ILEV1 = IP_ERA(INLEVEL+ISTARTLEVEL-JLOOP1)
    END IF
  END IF
  ! read component u 
  CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR,KLEV1=ILEV1)
  IF (INUM < 0) THEN
    WRITE(YMSG,*) 'wind vector component "u" at level ',JLOOP1,' is missing'
    CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
  END IF
  CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
  ALLOCATE(ZVALUE(ISIZE))
  CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
  IF (IMODEL==3.OR.(IMODEL==7)) THEN
    ALLOCATE(ZTU0_LS(INI))
    ZTU0_LS(:) = ZVALUE(:)
  ELSE
    CALL GRIB_GET(IGRIB(INUM),'Nj',INJ,IRET_GRIB)
    IF(ALLOCATED(IINLO)) DEALLOCATE (IINLO)
    ALLOCATE(IINLO(INJ))      
    CALL COORDINATE_CONVERSION(IMODEL,IGRIB(INUM),IIU,IJU,ZLONOUT,ZLATOUT,&
          ZXOUT,ZYOUT,INI,ZPARAM,IINLO)
    ALLOCATE(ZOUT(INO))
    CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
           ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.TRUE.,PTIME_HORI,.FALSE. )
    ZTU_LS(:) = ZOUT(:)
    DEALLOCATE(IINLO)
    DEALLOCATE(ZOUT)
  END IF
  DEALLOCATE (ZVALUE)
  ! read component v and perform interpolation if not Arpege grid
  IF (IMODEL/=10.AND.IMODEL/=11) THEN ! others than NCEP and ERA5
    ILEV1 = JLOOP1
  ELSE
    IF(IMODEL==10) THEN
      ILEV1 = IP_GFS(INLEVEL+ISTARTLEVEL-JLOOP1)
    ELSE
      ILEV1 = IP_ERA(INLEVEL+ISTARTLEVEL-JLOOP1)
    END IF
  END IF
  CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR+1,KLEV1=ILEV1)
  IF (INUM < 0) THEN
    WRITE(YMSG,*) 'wind vector component "v" at level ',JLOOP1,' is missing'
    CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
  END IF
  CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
  ALLOCATE(ZVALUE(ISIZE))
  CALL GRIB_GET(IGRIB(INUM),'values',ZVALUE)
  IF ((IMODEL==3).OR.(IMODEL==7)) THEN 
    CALL GRIB_GET(IGRIB(INUM),'Nj',INJ,IRET_GRIB)
    ALLOCATE(IINLO(INJ))         
    CALL COORDINATE_CONVERSION(IMODEL,IGRIB(INUM),IIU,IJU,ZLONOUT,ZLATOUT,&
          ZXOUT,ZYOUT,INI,ZPARAM,IINLO)     
    ALLOCATE(ZTV0_LS(INI))
    ZTV0_LS(:) = ZVALUE(:)
  ELSE
     CALL GRIB_GET(IGRIB(INUM),'Nj',INJ,IRET_GRIB)
    ALLOCATE(IINLO(INJ))          
    CALL COORDINATE_CONVERSION(IMODEL,IGRIB(INUM),IIU,IJU,ZLONOUT,ZLATOUT,&
          ZXOUT,ZYOUT,INI,ZPARAM,IINLO)
    ALLOCATE(ZOUT(INO))
    CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
           ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.TRUE.,PTIME_HORI,.FALSE.)
    ZTV_LS(:) = ZOUT(:)
    DEALLOCATE(ZOUT)
  END IF
  DEALLOCATE (ZVALUE)
   ! interpolations for arpege grid
  IF ((IMODEL==3).OR.(IMODEL==7)) THEN
    ! Comes back to real winds instead of stretched winds
    ! (but still with components according to Arpege grid axes)
    ZLATPOLE = ZPARAM(7) * XPI/180.
    ZLONPOLE = ZPARAM(8) * XPI/180.
    ZC       = ZPARAM(9)
    ZD       = ZC * ZC
    JLOOP3 = 0
    JLOOP4 = 1
    ZLAT = ZPARAM(3) * XPI / 180.
    DO JLOOP2=1, INI
      ZLON = JLOOP3 * 2. * XPI / IINLO(JLOOP4)
      ! Compute the scale factor
      ZA = ((1.+ZD) - (1.-ZD)*SIN(ZLAT)) / (2. * ZC)
      ZTU0_LS(JLOOP2) = ZTU0_LS(JLOOP2) * ZA
      ZTV0_LS(JLOOP2) = ZTV0_LS(JLOOP2) * ZA
      ! next parallel
      JLOOP3 = JLOOP3 + 1
      IF (JLOOP3 == IINLO(JLOOP4)) THEN
        JLOOP3 = 0
        ZLAT = ZLAT + (((ZPARAM(5)-ZPARAM(3))/(ZPARAM(2)-1)) * XPI / 180.)
        JLOOP4 = JLOOP4 + 1
      END IF
    END DO
    !
    ! interpolation
    CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,&
                INI,ZTU0_LS,INO,ZXOUT,ZYOUT,ZTU_LS,.TRUE.,PTIME_HORI,.FALSE.)
    CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,&
                INI,ZTV0_LS,INO,ZXOUT,ZYOUT,ZTV_LS,.TRUE.,PTIME_HORI,.FALSE.)
    DEALLOCATE(IINLO)
    !
    ! Rotation of the components from Arpege grid axes to real sphere axes
    !
    DO JLOOP2=1, INO
      ZLAT = ZYOUT(JLOOP2) * XPI / 180.
      ZLON = ZXOUT(JLOOP2) * XPI / 180.
      ! Compute the rotation matrix
      ZA = (ZD+1.) + (ZD-1.)*SIN(ZLAT)
      ZB = (ZD-1.) + (ZD+1.)*SIN(ZLAT)
      ZE = 2.*ZC*COS(ZLATPOLE)*COS(ZLAT)*COS(ZLON) + ZB*SIN(ZLATPOLE)
      IF (ABS(ZE) .GE. ABS(ZA)) THEN
        ZF = -2.*ZC*COS(ZLATPOLE)/ ( COS(ZLAT)* ((ZD+1.)+(ZD-1.)*SIN(ZLATPOLE)) )
        ZSIN = -ZF*SIN(ZLONPOLE-ZLON)
        ZCOS =  ZF*COS(ZLONPOLE-ZLON)
      ELSE
        ZF = 1. / SQRT(ZA*ZA - ZE*ZE)
        ZSIN = -COS(ZLATPOLE)*SIN(ZLON)*ZA*ZF
        ZCOS = (2.*ZC*SIN(ZLATPOLE)*COS(ZLAT)-ZB*COS(ZLATPOLE)*COS(ZLON))*ZF
      ENDIF
      ZTEMP = ZTU_LS(JLOOP2)
      ZTU_LS(JLOOP2) = ZCOS*ZTEMP - ZSIN*ZTV_LS(JLOOP2)
      ZTV_LS(JLOOP2) = ZSIN*ZTEMP + ZCOS*ZTV_LS(JLOOP2)
    END DO
  END IF
  !
  ! Rotation of the components from the real sphere axes (Arpege, CEP)
  ! or model axes (Aladin) to MESO-NH axes
  !
  JLOOP4=0
  DO JJ=1,IJU
    DO JI=1,IIU
      JLOOP4=JLOOP4+1
      IF (IMODEL==2 .OR. IMODEL==1 ) THEN 
          IF (IMODEL==2) THEN          ! ALADIN  REUNION
              ZALPHA=0
          ELSE  !ALADIN
              ZALPHA = (XRPK*(ZLONOUT(JLOOP4)-XLON0)-XBETA) - &
              (SIN(ZPARAM(9)*XPI/180.)*(ZLONOUT(JLOOP4)-ZPARAM(10)))
          ENDIF
      ELSE  ! CEP, ARPEGE (after projection)
          ZALPHA = XRPK*(ZLONOUT(JLOOP4)-XLON0)-XBETA
      ENDIF
      ZALPHA = ZALPHA * XPI / 180.
      XU_LS(JI,JJ,INLEVEL+ISTARTLEVEL-JLOOP1)= &
          ZTU_LS(JLOOP4)*COS(ZALPHA) - ZTV_LS(JLOOP4)*SIN(ZALPHA)
      XV_LS(JI,JJ,INLEVEL+ISTARTLEVEL-JLOOP1)= &
          ZTU_LS(JLOOP4)*SIN(ZALPHA) + ZTV_LS(JLOOP4)*COS(ZALPHA)
    ENDDO
  ENDDO
  IF ((IMODEL==3).OR.(IMODEL==7)) THEN ! deallocation of Arpege arrays
    DEALLOCATE (ZTU0_LS)
    DEALLOCATE (ZTV0_LS)
  END IF
END DO
DEALLOCATE (ZTU_LS)
DEALLOCATE (ZTV_LS)
IF(ALLOCATED(IINLO)) DEALLOCATE (IINLO)
END IF
!
!-------------------------------------------------------------------------------
!* 2.8 Filter the characteristics of the large-scale vortex
!-------------------------------------------------------------------------------
IF (HFILE(1:3)=='ATM' .AND. LFILTERING) THEN
  WRITE (ILUOUT0,'(A)') ' | Starting the filtering of the fields to remove large-scale vortex'
  IF (INDEX(CFILTERING,'Q')/=0) THEN
    WRITE (ILUOUT0,'(A)') ' -> Filtering of Q is now available!'       
    WRITE (ILUOUT0,'(A,A5)') ' CFILTERING= ',CFILTERING
  ENDIF
  !
  IF (INDEX(CFILTERING,'P')/=0) THEN
    ! compute reduced surface pressure
    ALLOCATE(ZTVF_LS(IIU,IJU),ZMSLP_LS(IIU,IJU))
    ! compute pressure reduced to first level above mean sea level
    !                                        (rather than above ground level)
    ZGAMREF=-6.5E-3
    !virtual temperature at the first level above ground
    ZTVF_LS(:,:) = XT_LS(:,:,1)*(1.+XQ_LS(:,:,1,1)*(XRV/XRD-1))
    !virtual temperature averaged between first level above ground
    !                                 and first level above sea level
    ZTVF_LS(:,:) = ZTVF_LS(:,:)-0.5*ZGAMREF*XZS_LS(:,:)
    ZMSLP_LS(:,:)=XPS_LS(:,:)*EXP(XG*XZS_LS(:,:)/(XRD*ZTVF_LS(:,:)))
  ENDIF
  !
  IF (INDEX(CFILTERING,'P')==0) THEN
    IF (INDEX(CFILTERING,'Q')==0) THEN
      CALL REMOVAL_VORTEX(XZS_LS,XU_LS,XV_LS,XT_LS)
    ELSE
      CALL REMOVAL_VORTEX(XZS_LS,XU_LS,XV_LS,XT_LS,PQ_LS=XQ_LS(:,:,:,1))
    ENDIF
  ELSE
    IF (INDEX(CFILTERING,'Q')==0) THEN
      CALL REMOVAL_VORTEX(XZS_LS,XU_LS,XV_LS,XT_LS,PPS_LS=ZMSLP_LS)
    ELSE
      CALL REMOVAL_VORTEX(XZS_LS,XU_LS,XV_LS,XT_LS,PQ_LS=XQ_LS(:,:,:,1),PPS_LS=ZMSLP_LS)
    ENDIF
    XPS_LS(:,:)  = ZMSLP_LS(:,:)*EXP(-XG*XZS_LS(:,:)/(XRD*ZTVF_LS(:,:)))
    DEALLOCATE(ZTVF_LS,ZMSLP_LS)
  ENDIF
  !
END IF
!
!---------------------------------------------------------------------------------------
!* 2.9 Read date
!---------------------------------------------------------------------------------------
!
WRITE (ILUOUT0,'(A)') ' | Reading date'
CALL GRIB_GET(IGRIB(INUM),'dataDate',IDATE,IRET_GRIB)
CALL GRIB_GET(IGRIB(INUM),'dataTime',ITIME,IRET_GRIB)
TPTCUR%xtime=ITIME/100*3600+(ITIME-(ITIME/100)*100)*60
TPTCUR%nyear=IDATE/10000
TPTCUR%nmonth=INT((IDATE-TPTCUR%nyear*10000)/100)
TPTCUR%nday=IDATE-TPTCUR%nyear*10000-TPTCUR%nmonth*100
CALL GRIB_GET(IGRIB(INUM),'startStep',ITIMESTEP,IRET_GRIB)
CALL GRIB_GET(IGRIB(INUM),'stepUnits',CSTEPUNIT,IRET_GRIB)
IF (IMODEL==0.OR.IMODEL==11) THEN
  ITWOZS=0
  IF ((TPTCUR%nyear ==2000).AND.(TPTCUR%nmonth  >11)) ITWOZS=1
  IF ((TPTCUR%nyear ==2000).AND.(TPTCUR%nmonth ==11)) THEN
    IF ( (TPTCUR%nday   >20 )   .OR.   &
        ((TPTCUR%nday  ==20 ).AND.(TPTCUR%xtime >=64800 ))) ITWOZS=1
  END IF
  IF ( TPTCUR%nyear ==2001 )  ITWOZS=1
  IF ((TPTCUR%nyear ==2002).AND.(TPTCUR%nmonth  <11)) ITWOZS=1
  IF ((TPTCUR%nyear ==2002).AND.(TPTCUR%nmonth ==11)) THEN
    IF ( (TPTCUR%nday   <24 )   .OR.   &
        ((TPTCUR%nday  ==25 ).AND.(TPTCUR%xtime  <64800 ))) ITWOZS=1
  END IF
  IF (ITWOZS==1) &
    WRITE(ILUOUT0,*) ' Check that both orography fields on 1st model level and on surface are used.'
END IF

CALL MPPDB_CHECK3D(XU_LS,"XU_LS",PRECISION)
CALL MPPDB_CHECK3D(XV_LS,"XV_LS",PRECISION)

SELECT CASE (CSTEPUNIT)       ! Time unit indicator
  CASE ('h')                    !hour
    TPTCUR%xtime = TPTCUR%xtime + ITIMESTEP*3600.
  CASE ('m')                    !minute
    TPTCUR%xtime = TPTCUR%xtime + ITIMESTEP*60.
  CASE ('s')                    !minute
    TPTCUR%xtime = TPTCUR%xtime + ITIMESTEP
   CASE DEFAULT
    WRITE (ILUOUT0,'(A,A,A)') ' | error CSTEPUNIT=',CSTEPUNIT, ' is different of s,m or h'
END SELECT
CALL DATETIME_CORRECTDATE(TPTCUR)
IF (HFILE(1:3)=='ATM') THEN
  CALL SM_PRINT_TIME(TPTCUR,TLUOUT0,'MESONH current date')
  TDTCUR = TPTCUR
  TDTMOD = TPTCUR
  TDTSEG = TPTCUR
  TDTEXP = TPTCUR
ELSE IF (HFILE=='CHEM') THEN
  CALL SM_PRINT_TIME(TPTCUR,TLUOUT0,'current date in MesoNH format')
ENDIF
!
!-------------------------------------------------------------------------------
!* 2.10 Read and interpolate dummy fields listed in free-format part of nml file
!-------------------------------------------------------------------------------
IF (ODUMMY_REAL) THEN
  !
  WRITE (ILUOUT0,'(A)') ' | Try to read 2D dummy fields'
  !
  !*       2.10.1   read 2D dummy fields
  !
  ! close file
  CALL IO_File_close(TPPRE_REAL1)
  ! open input file
  CALL CH_OPEN_INPUT(TPPRE_REAL1%CNAME, "DUMMY_2D", TZFILE, ILUOUT0, KVERB)
  ICHANNEL = TZFILE%NLU
  !
  ! read number of dummy 2D fields to transfer into mesonh
  READ(ICHANNEL, *) IMOC
  IF (KVERB >= 5) WRITE(ILUOUT0,*) "number of dummy fields to transfer into Mesonh : ", IMOC
  ALLOCATE(XDUMMY_2D(IIU,IJU,IMOC),CDUMMY_2D(IMOC))
  ALLOCATE(INUMGRIB(IMOC),INUMLEV(IMOC),INUMLEV1(IMOC),INUMLEV2(IMOC))
  INUMLEV(:)=-1 ; INUMLEV1(:)=-1 ; INUMLEV2(:)=-1
  !
  IVAR=0
  ! read variables names and Grib codes
  DO JI = 1, IMOC
    READ(ICHANNEL,'(A)') YINPLINE
    YINPLINE= TRIM(ADJUSTL(YINPLINE))
    IF (LEN_TRIM(YINPLINE) == 0) CYCLE ! skip blank line
    ! transform tab and comma character into blank
    DO JJ=1,LEN_TRIM(YINPLINE)
      IF (YINPLINE(JJ:JJ)==YPTAB .OR. YINPLINE(JJ:JJ)==YPCOM) YINPLINE(JJ:JJ)= ' '
    END DO
    IF (KVERB >= 10) WRITE(ILUOUT0,*) 'Line read : ', YINPLINE
    ! extract field name
    INDX= INDEX(YINPLINE,' ')
    YFIELD= YINPLINE(1:INDX-1)
    IF (KVERB >= 5) WRITE(ILUOUT0,*) 'Field being treated : ', YFIELD
    ITYP=105
    ILEV1=-1
    ILEV2=-1
    ! extract the parameter indicator
    YINPLINE= ADJUSTL(YINPLINE(INDX:))
    INDX= INDEX(YINPLINE,' ')
    IF (INDX == 1) THEN
      WRITE(ILUOUT0,*) ' Parameter indicator is missing. ',YFIELD,' not treated.'
      CYCLE
    END IF
    IVAR=IVAR+1
    READ(YINPLINE(1:INDX-1),*) IPAR
    IF (NVERB>=5) WRITE(ILUOUT0,*) ' Parameter indicator: ',IPAR
    ! extract the level indicator (optional)
    YINPLINE= ADJUSTL(YINPLINE(INDX:))
    INDX= INDEX(YINPLINE,' ')
    IF (INDX /= 1) THEN
      READ(YINPLINE(1:INDX-1),*) ITYP
      IF (NVERB>=5) WRITE(ILUOUT0,*) ' Level indicator is indicated: ',ITYP 
    END IF
    ! extract the first level value (optional)
    YINPLINE= ADJUSTL(YINPLINE(INDX:))
    INDX= INDEX(YINPLINE,' ')
    IF (INDX /= 1) THEN
      READ(YINPLINE(1:INDX-1),*) ILEV1
      IF (NVERB>=5) WRITE(ILUOUT0,*) ' Level1 value is indicated: ',ILEV1
    END IF
    ! extract the second level value (optional)
    YINPLINE= ADJUSTL(YINPLINE(INDX:))
    INDX= INDEX(YINPLINE,' ')
    IF (INDX /= 1) THEN
      READ(YINPLINE(1:INDX-1),*) ILEV2
      IF (NVERB>=5) WRITE(ILUOUT0,*) ' Level2 value is indicated: ',ILEV2
    END IF
    !
    CDUMMY_2D(IVAR)=YFIELD ; INUMGRIB(IVAR)=IPAR
    INUMLEV(IVAR)=ITYP     ; INUMLEV1(IVAR)=ILEV1 ; INUMLEV2(IVAR)=ILEV2
    !
  END DO
  !
  CALL IO_File_close(TZFILE)
  TZFILE => NULL()
  !
  IF (NVERB>=10) THEN
    WRITE(ILUOUT0,*) CDUMMY_2D(1:IVAR)
    WRITE(ILUOUT0,*) INUMGRIB(1:IVAR)
    WRITE(ILUOUT0,*) INUMLEV(1:IVAR)
    WRITE(ILUOUT0,*) INUMLEV1(1:IVAR)
    WRITE(ILUOUT0,*) INUMLEV2(1:IVAR)
  END IF
  !
  IF (IVAR /= IMOC) THEN
    WRITE (ILUOUT0,'(A,I3,A,I3,A)') ' -> Number of correct lines (',IVAR,') is different of ',IMOC,' - abort'
    WRITE(YMSG,*) 'number of correct lines (',IVAR,') is different of ',IMOC
    CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
  END IF
  !
  !* 2.10.2 read and interpolate variables onto dummy variables XDUMMY_2D
  !
  DO JI = 1, IMOC
    WRITE(ILUOUT0,'(A,4(I3,1X))') CDUMMY_2D(JI),INUMGRIB(JI),INUMLEV(JI),INUMLEV1(JI),INUMLEV2(JI)
    CALL SEARCH_FIELD(IGRIB,INUM,KPARAM=IPAR,KLEV1=ILEV1)
    IF (INUM < 0) THEN
      WRITE (ILUOUT0,'(A,I3,A,I2,A)') ' -> 2D field ',INUMGRIB(JI),' is missing - abort'
      WRITE(YMSG,*) '2D field ',INUMGRIB(JI),' is missing'
      CALL PRINT_MSG(NVERB_FATAL,'IO','READ_ALL_DATA_GRIB_CASE',YMSG)
    END IF
    CALL GRIB_GET(IGRIB(INUM),'Nj',INJ,IRET_GRIB)
    ALLOCATE(IINLO(INJ))
    CALL COORDINATE_CONVERSION(IMODEL,IGRIB(INUM),IIU,IJU,ZLONOUT,ZLATOUT,&
        ZXOUT,ZYOUT,INI,ZPARAM,IINLO)
    CALL GRIB_GET_SIZE(IGRIB(INUM),'values',ISIZE)
    ALLOCATE(ZVALUE(ISIZE))
    CALL GRIB_GET(IGRIB(INUM_ZS),'values',ZVALUE)
    ALLOCATE(ZOUT(INO))
    CALL HORIBL(ZPARAM(3),ZPARAM(4),ZPARAM(5),ZPARAM(6),INT(ZPARAM(2)),IINLO,INI, &
         ZVALUE,INO,ZXOUT,ZYOUT,ZOUT,.FALSE.,PTIME_HORI,.FALSE. )
    DEALLOCATE(IINLO)
    DEALLOCATE(ZVALUE)
    CALL ARRAY_1D_TO_2D(INO,ZOUT,IIU,IJU,XDUMMY_2D(:,:,JI))
    DEALLOCATE (ZOUT)
  END DO
!
ENDIF
!
!---------------------------------------------------------------------------------------
!
!* 3. VERTICAL GRID
!
IF (HFILE(1:3)=='ATM') THEN
  WRITE (ILUOUT0,'(A)') ' | Reading of vertical grid in progress'
  CALL READ_VER_GRID(TPPRE_REAL1)
END IF

!
!---------------------------------------------------------------------------------------
!
!* 4. Free all temporary allocations
!
DEALLOCATE (ZLATOUT)
DEALLOCATE (ZLONOUT)
DEALLOCATE (ZXOUT)
DEALLOCATE (ZYOUT)
DEALLOCATE(ZPARAM)
DEALLOCATE(ZPARAM_ZS)
DEALLOCATE(IINLO_ZS)
DO JLOOP=1,ICOUNT
  CALL GRIB_RELEASE(IGRIB(JLOOP))
ENDDO
DEALLOCATE(IGRIB)

WRITE (ILUOUT0,'(A,A4,A)') ' -- Grib decoder for ',HFILE,' file ended successfully'
!
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!
!

!
CONTAINS
!
!
!     ##########################################################################
      SUBROUTINE ARRAY_1D_TO_2D (KN1,P1,KL1,KL2,P2)
!     ##########################################################################
!
!       Small routine used to store a linear array into a 2 dimension array
!
IMPLICIT NONE
INTEGER,                INTENT(IN)  :: KN1
REAL,DIMENSION(KN1),    INTENT(IN)  :: P1
INTEGER,                INTENT(IN)  :: KL1
INTEGER,                INTENT(IN)  :: KL2
REAL,DIMENSION(KL1,KL2),INTENT(OUT) :: P2
INTEGER                 :: JLOOP1_A1T2
INTEGER                 :: JLOOP2_A1T2
INTEGER                 :: JPOS_A1T2
!
IF (KN1 < KL1*KL2) THEN
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ARRAY_1D_TO_2D','sizes do not match')
END IF
JPOS_A1T2 = 1
DO JLOOP2_A1T2 = 1, KL2
  DO JLOOP1_A1T2 = 1, KL1
    P2(JLOOP1_A1T2,JLOOP2_A1T2) = P1(JPOS_A1T2)
    JPOS_A1T2 = JPOS_A1T2 + 1
  END DO
END DO
END SUBROUTINE ARRAY_1D_TO_2D
!
!
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!#################################################################################
SUBROUTINE SEARCH_FIELD(KGRIB,KNUM,KPARAM,KDIS,KCAT,KNUMBER,KLEV1,KTFFS)
!#################################################################################
! search the grib message corresponding to KPARAM,KLTYPE,KLEV1,KLEV2 in all 
! the KGIRB messages
!
USE MODD_LUNIT
USE GRIB_API
!
IMPLICIT NONE
!
!
INTEGER(KIND=kindOfInt),DIMENSION(:),INTENT(IN) :: KGRIB ! number of grib messages
INTEGER,INTENT(OUT)             :: KNUM  ! number of the message researched
INTEGER,INTENT(IN),OPTIONAL     :: KPARAM ! INdicator of parameter/paramId
INTEGER,INTENT(IN),OPTIONAL     :: KDIS ! Discipline (GRIB2)
INTEGER,INTENT(IN),OPTIONAL     :: KCAT ! Catégorie (GRIB2)
INTEGER,INTENT(IN),OPTIONAL     :: KNUMBER ! parameterNumber (GRIB2)
INTEGER,INTENT(IN),OPTIONAL     :: KLEV1  ! Level 
INTEGER,INTENT(IN),OPTIONAL     :: KTFFS  ! TypeOfFirstFixedSurface 
!
! Declaration of local variables
!
INTEGER :: IFOUND  ! Number of correct parameters
INTEGER :: ISEARCH  ! Number of correct parameters to find
INTEGER :: IRET    ! error code 
INTEGER :: IPARAM,IDIS,ICAT,INUMBER,ITFFS
INTEGER :: ILEV1   ! Level parameter 1
INTEGER :: JLOOP   ! Dummy counter
INTEGER :: IVERSION
! Variables used to display messages
INTEGER :: ILUOUT0   ! Logical unit number of the listing
!
ILUOUT0 = TLUOUT0%NLU
!
ISEARCH=0
! Initialize as not found
KNUM = -1
!
IF (PRESENT(KPARAM)) ISEARCH=ISEARCH+1
IF (PRESENT(KDIS)) ISEARCH=ISEARCH+1
IF (PRESENT(KCAT)) ISEARCH=ISEARCH+1
IF (PRESENT(KNUMBER)) ISEARCH=ISEARCH+1
IF (PRESENT(KLEV1)) ISEARCH=ISEARCH+1
IF(PRESENT(KTFFS)) ISEARCH=ISEARCH+1
!
DO JLOOP=1,SIZE(KGRIB)
      IFOUND = 0
      ! 
      CALL GRIB_GET(KGRIB(JLOOP),'editionNumber',IVERSION,IRET_GRIB)
      IF (IRET_GRIB >   0) THEN
        WRITE (ILUOUT0,'(A)')' | Error encountered in the Grib file, skipping field'
        CYCLE
      ELSE IF (IRET_GRIB == -6) THEN
        WRITE (ILUOUT0,'(A)')' | ECMWF pseudo-Grib data encountered, skipping field'
        CYCLE
      ENDIF
      !
     IF (PRESENT(KTFFS)) THEN
        CALL GRIB_GET(KGRIB(JLOOP),'typeOfFirstFixedSurface',ITFFS,IRET_GRIB)
        IF (IRET_GRIB >   0) THEN
          WRITE (ILUOUT0,'(A)')' | Error encountered in the Grib file, skipping field'
          CYCLE
        ELSE IF (IRET_GRIB == -6) THEN
          WRITE (ILUOUT0,'(A)')' | ECMWF pseudo-Grib data encountered, skipping field'
          CYCLE
        ENDIF
        IF (ITFFS==KTFFS) THEN
          IFOUND = IFOUND + 1
        ELSE
          CYCLE
        ENDIF
      ENDIF

      IF (PRESENT(KPARAM)) THEN
        IF (IVERSION == 2) THEN
          CALL GRIB_GET(KGRIB(JLOOP),'paramId',IPARAM,IRET_GRIB)
        ELSE
          CALL GRIB_GET(KGRIB(JLOOP),'indicatorOfParameter',IPARAM,IRET_GRIB)
        ENDIF
        IF (IRET_GRIB >   0) THEN
          WRITE (ILUOUT0,'(A)')' | Error encountered in the Grib file, skipping field'
          CYCLE
        ELSE IF (IRET_GRIB == -6) THEN
          WRITE (ILUOUT0,'(A)')' | ECMWF pseudo-Grib data encountered, skipping field'
          CYCLE
        ENDIF
        IF (IPARAM==KPARAM) THEN
                IFOUND = IFOUND + 1
        ELSE
          CYCLE
        ENDIF
      ENDIF
      !
      IF (PRESENT(KDIS)) THEN
        CALL GRIB_GET(KGRIB(JLOOP),'discipline',IDIS,IRET_GRIB)
        IF (IRET_GRIB >   0) THEN
          WRITE (ILUOUT0,'(A)')' | Error encountered in the Grib file, skipping field'
          CYCLE
        ELSE IF (IRET_GRIB == -6) THEN
          WRITE (ILUOUT0,'(A)')' | ECMWF pseudo-Grib data encountered, skipping field'
          CYCLE
        ENDIF
        IF (IDIS==KDIS) THEN
          IFOUND = IFOUND + 1
        ELSE
          CYCLE
        ENDIF
      ENDIF
      IF (PRESENT(KCAT)) THEN
        CALL GRIB_GET(KGRIB(JLOOP),'parameterCategory',ICAT,IRET_GRIB)
        IF (IRET_GRIB >   0) THEN
          WRITE (ILUOUT0,'(A)')' | Error encountered in the Grib file, skipping field'
          CYCLE
        ELSE IF (IRET_GRIB == -6) THEN
          WRITE (ILUOUT0,'(A)')' | ECMWF pseudo-Grib data encountered, skipping field'
          CYCLE
        ENDIF
        IF (ICAT==KCAT) THEN
          IFOUND = IFOUND + 1
        ELSE
          CYCLE
        ENDIF
      ENDIF
      IF (PRESENT(KNUMBER)) THEN      
        CALL GRIB_GET(KGRIB(JLOOP),'parameterNumber',INUMBER,IRET_GRIB)
        IF (IRET_GRIB >   0) THEN
          WRITE (ILUOUT0,'(A)')' | Error encountered in the Grib file, skipping field'
          CYCLE
        ELSE IF (IRET_GRIB == -6) THEN
          WRITE (ILUOUT0,'(A)')' | ECMWF pseudo-Grib data encountered, skipping field'
          CYCLE
        ENDIF
        IF (INUMBER==KNUMBER) THEN
          IFOUND = IFOUND + 1
        ELSE
          CYCLE
        ENDIF
      ENDIF
      !
      IF(PRESENT(KLEV1)) THEN
        CALL GRIB_GET(KGRIB(JLOOP),'topLevel',ILEV1,IRET_GRIB)
        IF (IRET_GRIB >   0) THEN
          WRITE (ILUOUT0,'(A)')' | Error encountered in the Grib file, skipping field'
          CYCLE
        ELSE IF (IRET_GRIB == -6) THEN
          WRITE (ILUOUT0,'(A)')' | ECMWF pseudo-Grib data encountered, skipping field'
          CYCLE
        ENDIF
        IF (ILEV1==KLEV1) THEN
          IFOUND = IFOUND + 1
        ELSE
          CYCLE
        ENDIF
      ENDIF
      !
      IF (IFOUND == ISEARCH) THEN
          KNUM=JLOOP
          EXIT
      ELSE  ! field not found
          KNUM=-1
      END IF
END DO
!
END SUBROUTINE SEARCH_FIELD
!#################################################################################
SUBROUTINE COORDINATE_CONVERSION(KMODEL,KGRIB,KNOLON,KNOLARG,&
                           PLONOUT,PLATOUT,PLXOUT,PLYOUT,KNI,PPARAM,KINLO)
!#################################################################################
!perform coordinate conversion from lat/lon system to x,y (depends on the grib
! type)
!!   AUTHOR
!!   ------
!!
!!   G. Tanguy
!!
!!   MODIFICATIONS
!!   -------------
!!
!!   Original       08/06/2010

USE MODD_CST
USE MODI_LATLONTOXY
USE GRIB_API
!
IMPLICIT NONE
!
!
INTEGER(KIND=kindOfInt),INTENT(IN)                            :: KGRIB  ! number of the grib message
INTEGER,INTENT(IN)                            :: KMODEL ! number of the model
INTEGER,INTENT(OUT)                           :: KNI    ! number of points
INTEGER,INTENT(IN)                            :: KNOLON,KNOLARG  ! Number of output points
REAL,DIMENSION( KNOLON*KNOLARG),INTENT(IN)    :: PLONOUT ! Output coordinate,
REAL,DIMENSION( KNOLON*KNOLARG),INTENT(IN)    :: PLATOUT ! lat/lon system
REAL,DIMENSION( KNOLON*KNOLARG),INTENT(INOUT) :: PLXOUT  ! Converted output coordinates
REAL,DIMENSION( KNOLON*KNOLARG),INTENT(INOUT) :: PLYOUT  ! (depends on Grib Grid type)
REAL,DIMENSION(:),INTENT(INOUT)               :: PPARAM  ! output parameters of
!                                                the grid to avoid many calculations
INTEGER,DIMENSION(:),INTENT(INOUT)            :: KINLO   ! Number of points along a parallel
!===============================
INTEGER                                       :: IINLA ! Number of points along a meridian
INTEGER                                       :: JLOOP1,JLOOP2 ! Dummy counter
INTEGER                                       :: INO ! Number of output points
REAL                                          :: ZILA1  ! Grib first point latitude
REAL                                          :: ZILO1  ! Grib first point longitude
REAL                                          :: ZILA2  ! Grib last point latitude
REAL                                          :: ZILO2  ! Grib last point longitude
REAL                                          :: ZILASP ! Grib streching pole lat
REAL                                          :: ZILOSP ! Grib streching pole lon
LOGICAL                                       :: GREADY ! Used to test if projection is needed
INTEGER                                       :: ILENX ! nb points in X
INTEGER                                       :: ILENY ! nb points in Y
INTEGER                                       :: IEARTH  ! 
REAL                                          :: ZSTRECH ! streching of arpege grid
INTEGER(KIND=kindOfInt)                                       :: IMISSING ! dummy variable
! Aladin projection
REAL                                          :: ZALALAT0  ! Grid definition parameters
REAL                                          :: ZALALON0  !  |
REAL                                          :: ZALALATOR !  |
REAL                                          :: ZALALONOR !  |
REAL                                          :: ZALARPK   !  |
REAL, DIMENSION(:,:), ALLOCATABLE             :: ZXM       ! Intermediate arrays
REAL, DIMENSION(:,:), ALLOCATABLE             :: ZYM       !  |
REAL, DIMENSION(:,:), ALLOCATABLE             :: ZLONM     !  |
REAL, DIMENSION(:,:), ALLOCATABLE             :: ZLATM     !  |
! CEP projection
REAL, DIMENSION(:), ALLOCATABLE               :: ZLATGRIB
REAL, DIMENSION(:), ALLOCATABLE               :: ZLONGRIB
INTEGER                                       :: INBLATGRIB,INBLONGRIB
!JUAN
INTEGER(KIND=kindOfInt),DIMENSION(:),ALLOCATABLE            :: INLO_GRIB   ! Number of points along a parallel
!JUAN
!
!--------------------------------------------------------------------------------
!
!JUAN
ALLOCATE(INLO_GRIB(SIZE(KINLO)))
!JUAN
INO= KNOLON*KNOLARG
SELECT CASE (KMODEL)
!
CASE(0,5,11) ! CEP/MOCAGE/ERA5
! en theorie il faut ces 4 lignes
!  CALL GRIB_GET(KGRIB,'latitudeOfFirstGridPointInDegrees',ZILA1)
!  CALL GRIB_GET(KGRIB,'longitudeOfFirstGridPointInDegrees',ZILO1)
!  CALL GRIB_GET(KGRIB,'latitudeOfLastGridPointInDegrees',ZILA2)
!  CALL GRIB_GET(KGRIB,'longitudeOfLastGridPointInDegrees',ZILO2)
! pourtant au passage de GRIB1 a GRIB2 les arrondi etait fait differement
! et on n'obtenais pas les meme resultat entre un fichier grib1 et le meme 
! convertit en GRIB2
! Du coup en faisant ce qui suit on prend une valeur recalculee par grib_api
! suivant l'ordre N de la gausienne donc plus precise et donc la meme entre le
! GRIB1 et le GRIB2
  CALL GRIB_GET(KGRIB,'Nj',IINLA,IRET_GRIB)
  CALL GRIB_GET_SIZE(KGRIB,'latitudes',INBLATGRIB)
  CALL GRIB_GET_SIZE(KGRIB,'longitudes',INBLONGRIB)
  ALLOCATE(ZLATGRIB(INBLATGRIB))
  ALLOCATE(ZLONGRIB(INBLONGRIB))
  CALL GRIB_GET(KGRIB,'latitudes',ZLATGRIB)
  CALL GRIB_GET(KGRIB,'longitudes',ZLONGRIB)
  ZILA1=MAXVAL(ZLATGRIB)
  ZILO1=MINVAL(ZLONGRIB)
  ZILA2=MINVAL(ZLATGRIB)
  ZILO2=MAXVAL(ZLONGRIB)
  KNI=0
  CALL GRIB_IS_MISSING(KGRIB,'pl',IMISSING,IRET_GRIB)
  IF (IRET_GRIB /= 0 .OR. IMISSING==1)  THEN   ! pl not present
    CALL GRIB_GET(KGRIB,'Ni',INLO_GRIB(1),IRET_GRIB)
    INLO_GRIB(2:)=INLO_GRIB(1)
    KNI=IINLA*INLO_GRIB(1)
    GREADY= (PPARAM(1)==INLO_GRIB(1) .AND. PPARAM(2)==IINLA .AND.&
      PPARAM(3)==ZILA1 .AND. PPARAM(4)==ZILO1 .AND.&
      PPARAM(5)==ZILA2 .AND. PPARAM(6)==ZILO2)
   PPARAM(1)=INLO_GRIB(1)
   PPARAM(2)=IINLA
   PPARAM(3)=ZILA1
   PPARAM(4)=ZILO1 
   PPARAM(5)=ZILA2
   PPARAM(6)=ZILO2
  ELSE ! pl present in the grib
    CALL GRIB_GET(KGRIB,'pl',INLO_GRIB)
    DO JLOOP1=1 ,IINLA
      KNI = KNI + INLO_GRIB(JLOOP1)
    ENDDO
    GREADY= (PPARAM(1)==INLO_GRIB(1) .AND.&
      PPARAM(3)==ZILA1 .AND. PPARAM(4)==ZILO1 .AND.&
      PPARAM(5)==ZILA2 .AND. PPARAM(6)==ZILO2)
    PPARAM(1)=INLO_GRIB(1)
    PPARAM(2)=IINLA
    PPARAM(3)=ZILA1
    PPARAM(4)=ZILO1 
    PPARAM(5)=ZILA2
    PPARAM(6)=ZILO2    
  END IF
  IF (.NOT. GREADY) THEN
    PLXOUT=PLONOUT
    PLYOUT=PLATOUT
  ENDIF
!
CASE(1,6) ! ALADIN
!
  CALL GRIB_GET(KGRIB,'Nj',IINLA,IRET_GRIB)
  CALL GRIB_GET(KGRIB,'Ni',INLO_GRIB(1),IRET_GRIB)
  INLO_GRIB(2:)=INLO_GRIB(1)
  CALL GRIB_GET(KGRIB,'DxInMetres',ILENX)
  CALL GRIB_GET(KGRIB,'DyInMetres',ILENY)
  CALL GRIB_GET(KGRIB,'LoVInDegrees',ZALALON0)
  CALL GRIB_GET(KGRIB,'Latin1InDegrees',ZALALAT0)
  KNI = IINLA*INLO_GRIB(1)
  ZILA1 = 0.
  ZILO1 = 0.
  ZILA2 = ZILA1 + (IINLA   -1)*ILENY
  ZILO2 = ZILO1 + (INLO_GRIB(1)-1)*ILENX
  GREADY= (PPARAM(1)==INLO_GRIB(1) .AND. PPARAM(2)==IINLA .AND.&
      PPARAM(3)==ZILA1 .AND. PPARAM(4)==ZILO1 .AND.&
      PPARAM(5)==ZILA2 .AND. PPARAM(6)==ZILO2.AND.&
      PPARAM(7)==ILENX .AND. PPARAM(8)==ILENY.AND.&
      PPARAM(9)==ZALALAT0 .AND. PPARAM(10)==ZALALON0)
  IF(.NOT. GREADY) THEN
     PPARAM(1)=INLO_GRIB(1)
     PPARAM(2)=IINLA
     PPARAM(3)=ZILA1
     PPARAM(4)=ZILO1 
     PPARAM(5)=ZILA2
     PPARAM(6)=ZILO2    
     PPARAM(7)=ILENX
     PPARAM(8)=ILENY 
     PPARAM(9)=ZALALAT0
     PPARAM(10)=ZALALON0
!        
     IF (ZALALON0 > 180.) ZALALON0 = ZALALON0 - 360.
     CALL GRIB_GET(KGRIB,'latitudeOfFirstGridPointInDegrees',ZALALATOR)
    CALL GRIB_GET(KGRIB,'longitudeOfFirstGridPointInDegrees',ZALALONOR)
     IF (ZALALONOR > 180.) ZALALONOR = ZALALONOR - 360.
     ZALARPK = SIN(ZALALAT0/180.*XPI)
     ALLOCATE (ZXM(KNOLON,KNOLARG))
     ALLOCATE (ZYM(KNOLON,KNOLARG))
     ALLOCATE (ZLONM(KNOLON,KNOLARG))
     ALLOCATE (ZLATM(KNOLON,KNOLARG))
     JLOOP1=0
     DO JLOOP2=1, KNOLARG
       ZLONM(1:KNOLON,JLOOP2) = PLONOUT(1+JLOOP1:KNOLON+JLOOP1)
       ZLATM(1:KNOLON,JLOOP2) = PLATOUT(1+JLOOP1:KNOLON+JLOOP1)
       JLOOP1 = JLOOP1+KNOLON
     END DO
     CALL SM_LATLONTOXY_A (ZALALAT0,ZALALON0,ZALARPK,ZALALATOR,ZALALONOR, &
                           ZXM,ZYM,ZLATM,ZLONM,KNOLON,KNOLARG,6367470.)
     JLOOP1=0
     DO JLOOP2=1, KNOLARG
        PLXOUT(1+JLOOP1:KNOLON+JLOOP1)=ZXM(1:KNOLON,JLOOP2)
        PLYOUT(1+JLOOP1:KNOLON+JLOOP1)=ZYM(1:KNOLON,JLOOP2)
        JLOOP1 = JLOOP1+KNOLON
     ENDDO
     DEALLOCATE (ZLATM)
     DEALLOCATE (ZLONM)
     DEALLOCATE (ZYM)
     DEALLOCATE (ZXM)
  END IF
!
CASE(2) ! ALADIN REUNION
!
  CALL GRIB_GET(KGRIB,'Nj',IINLA,IRET_GRIB)
  CALL GRIB_GET(KGRIB,'Ni',INLO_GRIB(1),IRET_GRIB)
  INLO_GRIB(2:)=INLO_GRIB(1)
  CALL GRIB_GET(KGRIB,'DiInMetres',ILENX)
  CALL GRIB_GET(KGRIB,'DjInMetres',ILENY)
  CALL GRIB_GET(KGRIB,'LaDInDegrees',ZALALAT0)
  KNI = IINLA*INLO_GRIB(1)
  ZILA1 = 0.
  ZILO1 = 0.
  ZILA2 = ZILA1 + (IINLA   -1)*ILENY
  ZILO2 = ZILO1 + (INLO_GRIB(1)-1)*ILENX
  GREADY= (PPARAM(1)==INLO_GRIB(1) .AND. PPARAM(2)==IINLA .AND.&
      PPARAM(3)==ZILA1 .AND. PPARAM(4)==ZILO1 .AND.&
      PPARAM(5)==ZILA2 .AND. PPARAM(6)==ZILO2.AND.&
      PPARAM(7)==ILENX .AND. PPARAM(8)==ILENY.AND.&
      PPARAM(9)==ZALALAT0)
  IF(.NOT. GREADY) THEN
    PPARAM(1)=INLO_GRIB(1)
    PPARAM(2)=IINLA
    PPARAM(3)=ZILA1
    PPARAM(4)=ZILO1 
    PPARAM(5)=ZILA2
    PPARAM(6)=ZILO2    
    PPARAM(7)=ILENX
    PPARAM(8)=ILENY 
    PPARAM(9)=ZALALAT0
    ZALALON0 = 0.
    CALL GRIB_GET(KGRIB,'latitudeOfFirstGridPointInDegrees',ZALALATOR)
    CALL GRIB_GET(KGRIB,'longitudeOfFirstGridPointInDegrees',ZALALONOR)
    IF (ZALALONOR > 180.) ZALALONOR = ZALALONOR - 360.
    ZALARPK = 0
    ALLOCATE (ZXM(KNOLON,KNOLARG))
    ALLOCATE (ZYM(KNOLON,KNOLARG))
    ALLOCATE (ZLONM(KNOLON,KNOLARG))
    ALLOCATE (ZLATM(KNOLON,KNOLARG))
    JLOOP1=0
    DO JLOOP2=1, KNOLARG
       ZLONM(1:KNOLON,JLOOP2) = PLONOUT(1+JLOOP1:KNOLON+JLOOP1)
       ZLATM(1:KNOLON,JLOOP2) = PLATOUT(1+JLOOP1:KNOLON+JLOOP1)
       JLOOP1 = JLOOP1+KNOLON
    END DO
    CALL GRIB_GET(KGRIB,'earthIsOblate',IEARTH)
    IF (IEARTH==0) THEN
        CALL SM_LATLONTOXY_A (ZALALAT0,ZALALON0,ZALARPK,ZALALATOR,ZALALONOR, &
                              ZXM,ZYM,ZLATM,ZLONM,KNOLON,KNOLARG,6367470.)
    ELSE
        CALL SM_LATLONTOXY_A (ZALALAT0,ZALALON0,ZALARPK,ZALALATOR,ZALALONOR, &
                              ZXM,ZYM,ZLATM,ZLONM,KNOLON,KNOLARG)
    END IF
    JLOOP1=0
    DO JLOOP2=1, KNOLARG
       PLXOUT(1+JLOOP1:KNOLON+JLOOP1)=ZXM(1:KNOLON,JLOOP2)
       PLYOUT(1+JLOOP1:KNOLON+JLOOP1)=ZYM(1:KNOLON,JLOOP2)
       JLOOP1 = JLOOP1+KNOLON
    ENDDO
    DEALLOCATE (ZLATM)
    DEALLOCATE (ZLONM)
    DEALLOCATE (ZYM)
    DEALLOCATE (ZXM)       
  END IF
!
CASE(3,4,7) ! ARPEGE
!
!print*,"=========COORDINATE CONVERSION CASE ARPEGE ============="
! PROBLEME AVEC LES GRIB d'EPYGRAM
! dans longitudeOfLastGridPointInDegrees on la la longitude du dernier point du
! tableau (donc au pole sud)  
! dans les GRIB1 ont avait la valeur max du tableau des longitude (donc à
! l'equateur) 
  CALL GRIB_GET(KGRIB,'latitudeOfFirstGridPointInDegrees',ZILA1)
  CALL GRIB_GET(KGRIB,'longitudeOfFirstGridPointInDegrees',ZILO1)
  CALL GRIB_GET(KGRIB,'latitudeOfLastGridPointInDegrees',ZILA2)
  CALL GRIB_GET(KGRIB,'longitudeOfLastGridPointInDegrees',ZILO2)
  CALL GRIB_GET(KGRIB,'latitudeOfStretchingPoleInDegrees',ZILASP)
  CALL GRIB_GET(KGRIB,'longitudeOfStretchingPoleInDegrees',ZILOSP)
  CALL GRIB_GET(KGRIB,'stretchingFactor',ZSTRECH)
  CALL GRIB_GET(KGRIB,'Nj',IINLA,IRET_GRIB)
!
  KNI=0
  CALL GRIB_IS_MISSING(KGRIB,'pl',IRET_GRIB)
  IF (IRET_GRIB == 1)  THEN !  regular
     CALL GRIB_GET(KGRIB,'Ni',INLO_GRIB(1),IRET_GRIB)
     INLO_GRIB(2:)=INLO_GRIB(1)
     KNI=IINLA*INLO_GRIB(1)
     GREADY= (PPARAM(1)==INLO_GRIB(1) .AND. PPARAM(2)==IINLA .AND.&
       PPARAM(3)==ZILA1 .AND. PPARAM(4)==ZILO1 .AND.&
       PPARAM(5)==ZILA2 .AND. PPARAM(6)==ZILO2 .AND.&
       PPARAM(7)==ZILASP .AND. PPARAM(8)==ZILOSP .AND.&
       PPARAM(9)==ZSTRECH)
  ELSE !  quasi-regular
     CALL GRIB_GET(KGRIB,'pl',INLO_GRIB)
     DO JLOOP1=1 ,IINLA
        KNI = KNI + INLO_GRIB(JLOOP1)
     ENDDO
     ZILO2=360.-360./(MAXVAL(INLO_GRIB))
     GREADY= (PPARAM(1)==INLO_GRIB(1) .AND.&
       PPARAM(3)==ZILA1 .AND. PPARAM(4)==ZILO1 .AND.&
       PPARAM(5)==ZILA2 .AND. PPARAM(6)==ZILO2 .AND.&
       PPARAM(7)==ZILASP .AND. PPARAM(8)==ZILOSP .AND.&
       PPARAM(9)==ZSTRECH)
  END IF
!
  IF (.NOT. GREADY) THEN
    CALL ARPEGE_STRETCH_A(INO,ZILASP,ZILOSP, &
        ZSTRECH,PLATOUT,PLONOUT,PLYOUT,PLXOUT)
     PPARAM(1)=INLO_GRIB(1)
     PPARAM(2)=IINLA
     PPARAM(3)=ZILA1
     PPARAM(4)=ZILO1 
     PPARAM(5)=ZILA2
     PPARAM(6)=ZILO2
     PPARAM(7)=ZILASP
     PPARAM(8)=ZILOSP
     PPARAM(9)=ZSTRECH    
  END IF
!
CASE(10) ! NCEP
!
  CALL GRIB_GET(KGRIB,'latitudeOfFirstGridPointInDegrees',ZILA1)
  CALL GRIB_GET(KGRIB,'longitudeOfFirstGridPointInDegrees',ZILO1)
  CALL GRIB_GET(KGRIB,'latitudeOfLastGridPointInDegrees',ZILA2)
  CALL GRIB_GET(KGRIB,'longitudeOfLastGridPointInDegrees',ZILO2)
  CALL GRIB_GET(KGRIB,'Nj',IINLA,IRET_GRIB)
  CALL GRIB_GET(KGRIB,'Ni',INLO_GRIB(1),IRET_GRIB)
  INLO_GRIB(2:)=INLO_GRIB(1)
  KNI=IINLA*INLO_GRIB(1)
  GREADY= (PPARAM(1)==INLO_GRIB(1) .AND. PPARAM(2)==IINLA .AND.&
           PPARAM(3)==ZILA1 .AND. PPARAM(4)==ZILO1 .AND.&
           PPARAM(5)==ZILA2 .AND. PPARAM(6)==ZILO2)
  PPARAM(1)=INLO_GRIB(1)
  PPARAM(2)=IINLA
  PPARAM(3)=ZILA1
  PPARAM(4)=ZILO1 
  PPARAM(5)=ZILA2
  PPARAM(6)=ZILO2
  IF (.NOT. GREADY) THEN
    PLXOUT=PLONOUT
    PLYOUT=PLATOUT
  ENDIF
END SELECT
!JUAN
KINLO=INLO_GRIB
!JUAN
END SUBROUTINE COORDINATE_CONVERSION
!
!     ###################################################################
      SUBROUTINE ARPEGE_STRETCH_A(KN,PLAP,PLOP,PCOEF,PLAR,PLOR,PLAC,PLOC)
!     ###################################################################
!!****  *ARPEGE_STRETCH_A* - Projection to Arpege stretched grid
!!
!!   PURPOSE
!!   -------
!!
!!   Projection from standard Lat,Lon grid to Arpege stretched grid
!!
!!   METHOD
!!   ------
!!
!!   The projection is defined in two steps :
!!    1. A rotation to place the stretching pole at the north pole
!!    2. The stretching
!!   This routine is a basic implementation of the informations founded in
!!     'Note de travail Arpege nr.3'
!!     'Transformation de coordonnees'
!!     J.F.Geleyn 1988
!!   This document describes a slightly different transformation in 3 steps. Only the
!!   two first steps are to be taken in account (at the time of writing this paper has
!!   not been updated).
!!
!!   EXTERNAL
!!   --------
!!
!!   Module MODD_CST
!!     XPI
!!
!!   IMPLICIT ARGUMENTS
!!   ------------------
!!
!!   REFERENCE
!!   ---------
!!
!!   This routine is based on :
!!     'Note de travail ARPEGE' number 3
!!     by J.F. GELEYN (may 1988)
!!
!!   AUTHOR
!!   ------
!!
!!   V.Bousquet
!!
!!   MODIFICATIONS
!!   -------------
!!
!!   Original       07/01/1999
!!
!----------------------------------------------------------------------------
!
!*   0. DECLARATIONS
!    ---------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*   0.1. Declaration of arguments
!    -----------------------------
!
INTEGER,             INTENT(IN)  :: KN            ! Number of points to convert
REAL,                INTENT(IN)  :: PLAP          ! Latitude of stretching pole
REAL,                INTENT(IN)  :: PLOP          ! Longitude of stretching pole
REAL,                INTENT(IN)  :: PCOEF         ! Stretching coefficient
REAL, DIMENSION(KN), INTENT(IN)  :: PLAR          ! Lat. of points
REAL, DIMENSION(KN), INTENT(IN)  :: PLOR          ! Lon. of points
REAL, DIMENSION(KN), INTENT(OUT) :: PLAC          ! Computed pseudo-lat. of points
REAL, DIMENSION(KN), INTENT(OUT) :: PLOC          ! Computed pseudo-lon. of points
!
!*   0.2. Declaration of local variables
!    -----------------------------------
!
REAL                       :: ZSINSTRETCHLA ! Sine of stretching point lat.
REAL                       :: ZSINSTRETCHLO ! Sine of stretching point lon.
REAL                       :: ZCOSSTRETCHLA ! Cosine of stretching point lat.
REAL                       :: ZCOSSTRETCHLO ! Cosine of stretching point lon.
REAL                       :: ZSINLA        ! Sine of computed point latitude
REAL                       :: ZSINLO        ! Sine of computed point longitude
REAL                       :: ZCOSLA        ! Cosine of computed point latitude
REAL                       :: ZCOSLO        ! Cosine of computed point longitude
REAL                       :: ZSINLAS       ! Sine of point's pseudo-latitude
REAL                       :: ZSINLOS       ! Sine of point's pseudo-longitude
REAL                       :: ZCOSLOS       ! Cosine of point's pseudo-lon.
REAL                       :: ZA,ZB,ZD      ! Dummy variables used for
REAL                       :: ZX,ZY         ! computations
!
INTEGER                    :: JLOOP1        ! Dummy loop counter
!
!----------------------------------------------------------------------------
!
ZSINSTRETCHLA = SIN(PLAP*XPI/180.)
ZCOSSTRETCHLA = COS(PLAP*XPI/180.)
ZSINSTRETCHLO = SIN(PLOP*XPI/180.)
ZCOSSTRETCHLO = COS(PLOP*XPI/180.)
! L = longitude (0 = Greenwich, + toward east)
! l = latitude (90 = N.P., -90 = S.P.)
! p stands for stretching pole
PLAC(:) = PLAR(:) * XPI / 180.
PLOC(:) = PLOR(:) * XPI / 180.
! A = 1 + c.c
ZA = 1. + PCOEF*PCOEF
! B = 1 - c.c
ZB = 1. - PCOEF*PCOEF
DO JLOOP1=1, KN
  ZSINLA = SIN(PLAC(JLOOP1))
  ZCOSLA = COS(PLAC(JLOOP1))
  ZSINLO = SIN(PLOC(JLOOP1))
  ZCOSLO = COS(PLOC(JLOOP1))
  ! X = cos(Lp-L)
  ZX = ZCOSLO*ZCOSSTRETCHLO + ZSINLO*ZSINSTRETCHLO
  ! Y = sin(Lp-L)
  ZY = ZSINSTRETCHLO*ZCOSLO - ZSINLO*ZCOSSTRETCHLO
  ! D = (1+c.c) + (1-c.c)(sin lp.sin l + cos lp.cos l.cos(Lp-L))
  ZD = ZA + ZB*(ZSINSTRETCHLA*ZSINLA+ZCOSSTRETCHLA*ZCOSLA*ZX)
  !          (1-c.c)+(1+c.c)((sin lp.sin l + cos lp.cos l.cos(Lp-L))
  ! sin lr = -------------------------------------------------------
  !                                  D
  ZSINLAS = (ZB + ZA*(ZSINSTRETCHLA*ZSINLA+ZCOSSTRETCHLA*ZCOSLA*ZX)) / ZD
  ! D' = D * cos lr
  ZD = ZD * (AMAX1(1e-6,SQRT(1.-ZSINLAS*ZSINLAS)))
  !          2.c.(cos lp.sin l - sin lp.cos l.cos(Lp-L))
  ! cos Lr = -------------------------------------------
  !                              D'
  ZCOSLOS = 2.*PCOEF*(ZCOSSTRETCHLA*ZSINLA-ZSINSTRETCHLA*ZCOSLA*ZX) / ZD
  !          2.c.cos l.cos(Lp-L)
  ! sin Lr = -------------------
  !                  D'
  ZSINLOS = 2.*PCOEF*(ZCOSLA*ZY) / ZD
  ! saturations (corrects calculation errors)
  ZSINLAS = MAX(ZSINLAS,-1.)
  ZSINLAS = MIN(ZSINLAS, 1.)
  ZCOSLOS = MAX(ZCOSLOS,-1.)
  ZCOSLOS = MIN(ZCOSLOS, 1.)
  ! back from sine & cosine
  PLAC(JLOOP1) = ASIN(ZSINLAS)
  IF (ZSINLOS>0) THEN
    PLOC(JLOOP1) =  ACOS(ZCOSLOS)
  ELSE
    PLOC(JLOOP1) = -ACOS(ZCOSLOS)
  ENDIF
ENDDO
PLOC(:) = PLOC(:) * 180. / XPI
PLAC(:) = PLAC(:) * 180. / XPI
RETURN
END SUBROUTINE ARPEGE_STRETCH_A
!
!
END SUBROUTINE READ_ALL_DATA_GRIB_CASE
