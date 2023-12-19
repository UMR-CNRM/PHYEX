!MNH_LIC Copyright 1995-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      PROGRAM PREP_REAL_CASE
!     ######################
!
!!****  *PREP_REAL_CASE* - program to write an initial FM file from real case
!!                         situation.
!!
!!    PURPOSE
!!    -------
!!
!!       The purpose of this program is to prepare an initial meso-NH file
!!    (LFIFM and DESFM files) filled by some fields of a real situation. 
!!    General data are given by the MESO-NH user in the namelist file 
!!    'PRE_REAL1.nam'. The fields are obtained from three sources:
!!        - an atmospheric input file, which can be:
!!          * an Aladin file, itself obtained from an Arpege file with
!!            the Aladin routine "FULLPOS".
!!          * a grib file (ECMWF, Grib Arpege or Grib Aladin)
!!          * a MESONH file
!!        - an physiographic data file.
!!
!!       1) Fields obtained from the Atmospheric file:
!!          -----------------------------------------
!!
!!      - the projection parameters (checked with PGD file): 
!!                 reference latitude and longitude
!!                 parameter of projection
!!                 angle of rotation of the domain
!!
!!      - the horizontal grid definition (checked with PGD file):
!!                 grid mesh
!!                 latitude and longitude of the reference point 
!!                                   (with data from PRE_REAL1.nam)
!!
!!      - thermodynamical 3D and 2D fields:
!!                 potential temperature
!!                 vapor mixing ratio
!!
!!      - dynamical fields:
!!                 three components of the wind
!!
!!      - reference anelastic state variables:
!!                 profile of virtual potential temperature
!!                 profile of dry density
!!                 Exner function at model top
!!
!!      - total dry air mass
!!
!!
!!       2) Fields obtained from the physiographic data file:
!!          ------------------------------------------------
!!
!!      - the projection parameters: 
!!                 reference latitude and longitude
!!                 parameter of projection
!!                 angle of rotation of the domain
!!
!!      - the horizontal grid definition:
!!                 grid mesh
!!                 latitude and longitude of the reference point 
!!                                   (with data from PRE_REAL1.nam)
!!      - physiografic fields: (orographic, vegetation, soil and radiation fields)
!!
!!
!!       3) Data obtained from the namelist file PRE_REAL1.nam:
!!          --------------------------------------------------
!! 
!!      - type of equations system
!!      - vertical grid definition
!!      - number of points in x and y directions
!!      - level of verbosity
!!      - name of the different files
!!
!!
!!**  METHOD
!!    ------
!!      In this program, once the MESO-NH domain is calculated, all the
!!      2D or 3D fields are computed on the MESO-NH horizontal domain WITH 
!!      the external points. This is particularly important for the large
!!      scale fields during the MESO-NH run.
!!
!!   1) The following PREP_REAL_CASE program:
!!
!!      - set default values for global variables which will be written in 
!!        DESFM file (by calling DEFAULT_DESFM1); lateral boundary conditions
!!        are open.
!!
!!      - opens the different files (by calling OPEN_PRC_FILES).
!!
!!      - initializes physical constants (by calling INI_CST).
!!
!!      - initializes the horizontal domain from the data read in the 
!!        descriptive part of the Aladin file and the directives read in the
!!        namelist file (routines READ_GENERAL and SET_SUBDOMAIN in 
!!        READ_ALL_DATA). This MESO-NH domain is a part of the Aladin domain.
!!
!!      - initializes global variables from namelists and the MESO-NH
!!        vertical grid definition variables in the namelist file
!!        (routine READ_VER_GRID).
!!
!!      - initializes the physiographic 2D fields from the physiographic data
!!        file, in particular the MESO-NH orography.
!!
!!      - reads the 3D and 2D variable fields in the Grib file
!!        (routine READ_ALL_DATA_GRIB_CASE),
!!        if  HATMFILETYPE='GRIBEX':
!!           absolute temperature
!!           specific humidity
!!           horizontal contravariant wind
!!           surface pressure
!!           large scale orography
!!
!!      - reads the 3D and 2D variable fields in the input MESONH file
!!        (routine READ_ALL_DATA_MESONH_CASE), if HATMFILETYPE='MESONH':
!!           potential temperature
!!           vapor mixing ratio
!!           horizontal wind
!!           other mixing ratios
!!           turbulence prognostic and semi-prognostic variables
!!           large scale orography
!!
!!      - computes some geometric variables (routines SM_GRIDPROJ and METRICS),
!!        in particular:
!!          * altitude 3D array
!!          * metric coefficients
!!          * jacobian
!!
!!      - initializes MESO-NH thermodynamical fields:
!!          * changes of variables (routine VER_PREP_mmmmmm_CASE):
!!              absolute temperature --> virtual potential temperature
!!              specific humidity    --> vapor mixing ratio
!!          * interpolates/extrapolates the fields from the large scale 
!!            orography to the MESO-NH one (routine VER_INT_THERMO in
!!            VER_THERMO, by using a shifting function method).
!!            in water vapor case, the interpolations are always performed
!!            on relative humidity.
!!          * the pressure is computed on each grid by integration of the
!!            hydrostatic equation from bottom or top. When input atmospheric
!!            file is a MESO-NH one, information about the difference between
!!            hydrostatic pressure and total pressure is kept and interpolated
!!            during the entire PREP_REAL_CASE process.
!!          * interpolates the fields to the MESO-NH vertical grid
!!            (also by routine VER_INT_THERMO in VER_THERMO).
!!          * computes the potential temperature (routine VER_THERMO).
!!          * sets to zero the mixing ratios, except the vapor mixing ratio
!!            (VER_THERMO).
!!
!!      - initializes the reference anelastic state variables (routine SET_REFZ
!!        in VER_THERMO).
!!
!!      - computes the total dry air mass (routine DRY_MASS in VER_THERMO).
!!
!!      - initializes MESO-NH dynamical variables:
!!          * changes Aladin contravariant wind into true horizontal wind
!!            (in subroutine VER_PREP).
!!          * interpolates/extrapolates the momentum from the large scale 
!!            orography to the MESO-NH one (routine VER_INT_DYN in
!!            VER_DYN, by using a shifting function method).
!!          * interpolates the fields to the MESO-NH vertical grid
!!            (also by routine VER_INT_DYN in VER_DYN). The fields
!!            are located on a horizontal Arakawa A-grid, as the Aladin fields.
!!          * The momentum is interpolated to the Arakawa C-grid 
!!            (routine VER_DYN).
!!          * A first guess of the vertical momentum, verifying the 
!!            uncompressible continuity equation and the material lower boundary
!!            condition against the ground, is computed (routine WGUESS).
!!          * computes the final non-divergent wind field (routine 
!!            ANEL_BALANCE).
!!
!!      - copies the interpolated fields also at t-dt and in the large scale
!!        fields (routine INI_PROG_VAR).
!!  
!!      - writes the DESFM and LFIFM files (routines WRITE_DESFM1 and
!!        WRITE_LFIFM1).
!!
!!
!!   2) Some conventions are used in this program and its subroutines because
!!      of the number of different grids and fields:
!!
!!      - subscripts:
!!          * the subscripts I and J are used for all the horizontal grid.
!!          * the subcript K is used for the MESO-NH vertical grid (increasing
!!            from bottom to top).
!!          * the subscript L is used for the Aladin or input Mesonh grids
!!            (increasing from bottom to top).
!!
!!      - suffixes:
!!          * _LS: 
!!              If used for a geographic or horizontal grid definition variable,
!!               this variable is connected to the large horizontal domain.
!!              If used for a surface variable, this variable corresponds to
!!               the large scale orography, and therefore will be modified.
!!              If used for another variable, this variable is discretized
!!               on the Aladin or input MESONH file vertical grid 
!!               (large-scale orography with input vertical discretization,
!!                either coming from eta levels or input Gal-Chen grid).
!!          * _MX:
!!              Such a variable is discretized on the mixed grid.
!!              (large-scale orography with output Gal-Chen vertical grid
!!               discretization)
!!          * _SH:
!!              Such a variable is discretized on the shifted grid.
!!              (fine orography with a shifted vertical grid, NOT Gal-Chen)
!!          * no suffix:
!!              The variable is discretized on the MESO-NH grid.
!!              (fine orography with output Gal-Chen vertical grid discretization)
!!
!!      - additional pre-suffixes: (for pressure, Exner and altitude fields)
!!          * MASS:
!!              The variable is discretized on a mass point
!!          * FLUX:
!!              The variable is discretized on a flux point
!!
!!
!!      - names of variables: for a physical variable VAR:
!!          * pVARs        is the variable itself.
!!          * pRHODVARs    is the variable multiplied by the dry density rhod.
!!          * pRHODJVARs   is the variable multiplied by the dry density rhod
!!                         and the Jacobian.            
!!          * pRVARs       is the variable multiplied by rhod_ref, the anelastic
!!                         reference state dry density and the Jacobian.
!!        where p and s are the appropriate prefix and suffix.
!!
!!      - allocation of arrays: the arrays are allocated
!!          * just before their initialization for the general arrays stored in 
!!            modules.
!!          * in the subroutine in which they are declared for the local arrays 
!!            in a subroutine.
!!          * in the routine in which they are initialized for the arrays 
!!            defined in the monitor PREP_REAL_CASE. In this case they are in 
!!            fact passed as pointer to the subroutines to allow their 
!!            dynamical allocation (exception which confirms the rule: ZJ).
!!
!!
!!    EXTERNAL
!!    --------
!!
!!      Routine DEFAULT_DESFM1 : to set default values for variables which can be
!!                              contained in DESFM file.
!!      Routine OPEN_PRC_FILES: to open all files.
!!      Routine INI_CST       : to initialize physical constants.
!!      Routine READ_ALL_DATA_GRIB_CASE : to read all input data.
!!      Routine READ_ALL_DATA_MESONH_CASE : to read all input data.
!!      Routine SM_GRIDPROJ   : to compute some grid variables, in case of 
!!                              conformal projection.
!!      Routine METRICS       : to compute metric coefficients.
!!      Routine VER_PREP_GRIBEX_CASE      : to prepare the interpolations.
!!      Routine VER_PREP_MESONH_CASE      : to prepare the interpolations.
!!      Routine VER_THERMO    : to perform the interpolation of thermodynamical
!!                              variables.
!!      Routine VER_DYN       : to perform the interpolation of dynamical
!!                              variables.
!!      Routine INI_PROG_VAR  : to initialize the prognostic varaibles not yet
!!                              initialized
!!      Routine WRITE_DESFM1  : to write a DESFM file.
!!      Routine WRITE_LFIFM1  : to write a LFIFM file.
!!      Routine IO_File_close : to close a FM-file (DESFM + LFIFM).
!!
!!      Module MODE_GRIDPROJ  : contains conformal projection routines
!!    
!!      Module MODI_DEFAULT_DESFM1  : interface module for routine DEFAULT_DESFM1
!!      Module MODI_OPEN_PRC_FILES : interface module for routine OPEN_PRC_FILES
!!      Module MODI_READ_ALL_DATA_MESONH_CASE : interface module for routine
!!                                              READ_ALL_DATA_MESONH_CASE
!!      Module MODI_METRICS        : interface module for routine METRICS
!!      Module MODI_VER_PREP_GRIBEX_CASE      : interface module for routine
!!                                              VER_PREP_GRIBEX_CASE
!!      Module MODI_VER_PREP_MESONH_CASE      : interface module for routine
!!                                              VER_PREP_MESONH_CASE
!!      Module MODI_VER_THERMO     : interface module for routine VER_THERMO
!!      Module MODI_VER_DYN        : interface module for routine VER_DYN
!!      Module MODI_INI_PROG_VAR   : interface module for routine INI_PROG_VAR
!!      Module MODI_WRITE_DESFM1   : interface module for routine WRITE_DESFM1
!!      Module MODI_WRITE_LFIFM1   : interface module for routine WRITE_LFIFM1
!!      
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF  : contains configuration variables for all models.
!!         NVERB   : verbosity level for output-listing
!!      Module MODD_CONF1 : contains configuration variables for model 1.
!!         NRR     : number of moist variables
!!      Module MODD_LUNIT : contains logical unit and names of files.
!!      Module MODD_LUNIT : contains logical unit and names of files (model1).
!!         CINIFILE: name of the FM file which will be used for the MESO-NH run.
!!      Module MODD_GRID1 : contains grid variables.
!!         XLAT    : latitude of the grid points
!!         XLON    : longitudeof the grid points
!!         XXHAT   : position xhat in the conformal plane
!!         XYHAT   : position yhat in the conformal plane
!!         XDXHAT  : horizontal local meshlength on the conformal plane
!!         XDYHAT  : horizontal local meshlength on the conformal plane
!!         XZS     : MESO-NH orography
!!         XZZ     : altitude
!!         XZHAT   : height zhat
!!         XMAP    : map factor
!!      Module MODD_LBC1 : contains declaration of lateral boundary conditions
!!        CLBCX   : X-direction LBC type at left(1) and right(2) boundaries 
!!        CLBCY   : Y-direction LBC type at left(1) and right(2) boundaries   
!!      Module MODD_PARAM1 : contains declaration of the parameterizations' names
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
!!      Original    01/01/95
!!                  Sept. 21, 1995  (J.Stein and V.Masson) surface pressure
!!                  Jan.  09, 1996  (V. Masson) pressure function deduced from
!!                                  hydrostatic pressure
!!                  Jan.  31, 1996  (V. Masson) possibility to initialize
!!                                  atmospheric fields from MESONH file
!!                  Mar.  18, 1996  (V. Masson) new vertical extrapolation of Ts
!!                                  in case of initialization with MESONH file
!!                  Apr   17, 1996  (J. Stein ) change the DEFAULT_DESFM CALL
!!                  May   25, 1996  (V. Masson) Variable CSTORAGE_TYPE
!!                  Aug   26, 1996  (V. Masson) Only thinshell approximation is
!!                                  currently available.
!!                  Sept  24, 1996  (V. Masson) add writing of varaibles for
!!                                  nesting ('DAD_NAME', 'DXRATIO', 'DYRATIO')
!!                  Oct   11, 1996  (V. Masson) L1D and L2D configurations
!!                  Oct   28, 1996  (V. Masson) add deallocations and NVERB
!!                                  default set to 1
!!                  Dec   02, 1996  (V. Masson) vertical interpolation of
!!                                  surface fields in aladin case
!!                  Dec   12, 1996  (V. Masson) add LS vertical velocity
!!                  Jan   16, 1997  (J. Stein)  Durran's anelastic system
!!                  May   07, 1997  (V. Masson) add LS tke
!!                  Jun   27, 1997  (V. Masson) add absolute pressure
!!                  Jul   09, 1997  (V. Masson) add namelist NAM_REAL_CONF
!!                  Jul   10, 1997  (V. Masson) add LS epsilon
!!                  Aug   25, 1997  (V. Masson) add computing time analysis
!!                  Jan   20, 1998  (J. Stein)  add LB and LS fields
!!                  Apr,  30, 1998  (V. Masson) Large scale VEG and LAI
!!                  Jun,  04, 1998  (V. Masson) Large scale D2 and Aladin ISBA
!!                                              files
!!                  Jun,  04, 1998  (V. Masson) Add new soil interface var.
!!                  Jan   20, 1999  (J. Stein)  add a Boundaries call
!!                  March 15  1999  (J. Pettre, V. Bousquet and V. Masson)
!!                                              initialization from GRIB files
!!                  Jul       2000 (F.solmon/V.Masson) Adaptation for patch 
!!                                            according to GRIB or MESONH case
!!                  Nov   22, 2000  (P.Tulet, I. Mallet) initialization
!!                                              from GRIB MOCAGE file
!!                  Fev   01, 2001 (D.Gazen) add module MODD_NSV for NSV variable
!!                  Jul   02, 2001 (J.Stein) add LCARTESIAN case
!!                  Oct   15, 2001 (I.Mallet) allow namelists in different orders
!!                  Dec       2003 (V.Masson) removes surface calls
!!                  Jun   01, 2002 (O.Nuissier) filtering of tropical cyclone
!!                  Aou   09, 2005 (D.Barbary) add CDADATMFILE CDADBOGFILE
!!                   May   2006    Remove KEPS
!!                  Feb   02, 2012 (C. Mari) interpolation from MOZART
!!                                  add call to READ_CHEM_NETCDF_CASE &
!!                                  VER_PREP_NETCDF_CASE 
!!                  Mar   2012    Add NAM_NCOUT for netcdf output
!!                  July  2013     (Bosseur & Filippi) Adds Forefire
!!                  Mars  2014     (J.Escobar) Missing 'full' UPDATE_METRICS for arp2lfi // run
!!                   April 2014     (G.TANGUY) Add LCOUPLING
!!                        2014     (M.Faivre)
!!                  Fevr  2015     (M.Moge) Cleaning up
!!                  Aug   2015     (M.Moge) removing EXTRAPOL on XDXX and XDYY in part 8
!!    J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!    M.Leriche        2015 : add LUSECHEM  dans NAM_CH_CONF 
!!                  Feb   02, 2012 (C. Mari & BV) interpolation from CAMS
!!                                  add call to READ_CAMS_NETCDF_CASE &
!!                                  VER_PREP_NETCDF_CASE 
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!      Modification    02/2016  (JP Pinty) Convert CAMS mix ratio to nbr conc
!
!!  06/2016     (G.Delautier) phasage surfex 8
!!    P.Wautelet : 08/07/2016 : removed MNH_NCWRIT define
!!     B.VIE 2016 : LIMA
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 07/02/2019: force TYPE to a known value for IO_File_add2list
!  P. Wautelet 14/02/2019: remove CLUOUT/CLUOUT0 and associated variables
!  S. Bielli      02/2019: sea salt: significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 20/03/2019: missing use MODI_INIT_SALT
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  T.Nagel        02/2021: add IBM
!  P. Wautelet 06/07/2021: use FINALIZE_MNH
!!  M. Leriche 26/01/2022: add reading of CAMS reanalysis for chemistry
!!                         and/or for LIMA
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET,           ONLY: TBUCONF_ASSOCIATE
USE MODD_CH_M9_n
USE MODD_CH_MNHC_n,        ONLY: LUSECHAQ_n=>LUSECHAQ,LUSECHIC_n=>LUSECHIC, LUSECHEM_n=>LUSECHEM
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST
USE MODD_DIM_n
!UPG*PT
USE MODD_CH_AEROSOL
USE MODD_DUST,             ONLY:  LDUST, NMODE_DST, CRGUNITD, XINISIG, XINIRADIUS, XN0MIN,&
                                  LDSTCAMS
!UPG*PT

USE MODD_DYN_n,            CPRESOPT_n=>CPRESOPT, LRES_n=>LRES, XRES_n=>XRES , NITR_n=>NITR
USE MODD_FIELD_n
USE MODD_GR_FIELD_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_HURR_CONF
USE MODD_IBM_LSF,          ONLY: CIBM_TYPE, LIBM_LSF, NIBM_SMOOTH, XIBM_SMOOTH
USE MODD_IBM_PARAM_n,      ONLY: XIBM_LS
USE MODD_IO,               ONLY: TFILEDATA, TFILE_SURFEX
USE MODD_LBC_n
USE MODD_LES,              ONLY: LES_ASSOCIATE
USE MODD_LSFIELD_n
USE MODD_LUNIT,            ONLY: TPGDFILE,TLUOUT0,TOUTDATAFILE
USE MODD_LUNIT_n,          ONLY: CINIFILE,TINIFILE,TLUOUT
USE MODD_METRICS_n
USE MODD_MNH_SURFEX_n
USE MODD_NESTING
USE MODD_NSV
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA,       ONLY: PARAM_LIMA_INIT, NMOD_CCN, NMOD_IFN
USE MODD_PARAM_n
USE MODD_PREP_REAL
USE MODD_REF_n
!UPG*PT
USE MODD_SALT,             ONLY:  LSALT, NMODE_SLT, CRGUNITS, XINISIG_SLT, XINIRADIUS_SLT, XN0MIN_SLT,&
                                  LSLTCAMS
USE MODD_CH_AERO_n,        ONLY:  XM3D, XRHOP3D, XSIG3D, XRG3D, XN3D, XCTOTA3D
!UPG*PT
USE MODD_TURB_n
!
USE MODE_EXTRAPOL
use mode_field,            only: Alloc_field_scalars, Ini_field_list, Ini_field_scalars
USE MODE_FINALIZE_MNH,     only: FINALIZE_MNH
USE MODE_GRIDCART
USE MODE_GRIDPROJ
USE MODE_INI_CST,          ONLY: INI_CST
USE MODE_IO,               only: IO_Init
USE MODE_IO_FIELD_READ,    only: IO_Field_read
USE MODE_IO_FIELD_WRITE,   only: IO_Header_write
USE MODE_IO_FILE,          only: IO_File_close, IO_File_open
USE MODE_IO_MANAGE_STRUCT, only: IO_File_add2list, IO_File_find_byname
USE MODE_ll
USE MODE_MODELN_HANDLER
USE MODE_MPPDB
USE MODE_MSG
USE MODE_POS
USE MODE_SPLITTINGZ_ll
!
USE MODI_BOUNDARIES
USE MODI_COMPARE_DAD
USE MODI_DEALLOCATE_MODEL1
USE MODI_DEALLOC_PARA_LL
USE MODI_DEFAULT_DESFM_n
USE MODI_ERROR_ON_TEMPERATURE
USE MODI_IBM_INIT_LS
USE MODI_INI_PROG_VAR
USE MODI_INIT_SALT
USE MODI_LIMA_MIXRAT_TO_NCONC
USE MODI_METRICS
USE MODI_MNHREAD_ZS_DUMMY_n
USE MODI_MNHWRITE_ZS_DUMMY_n
USE MODI_OPEN_PRC_FILES
USE MODI_PREP_SURF_MNH
USE MODI_PRESSURE_IN_PREP
USE MODI_READ_ALL_DATA_GRIB_CASE
USE MODI_READ_ALL_DATA_MESONH_CASE
USE MODI_READ_ALL_NAMELISTS
!UPG*PT
!USE MODI_READ_CAMS_DATA_NETCDF_CASE
!USE MODI_READ_CHEM_DATA_NETCDF_CASE
USE MODI_READ_CHEM_DATA_MOZART_CASE
USE MODI_READ_CHEM_DATA_CAMS_CASE
USE MODI_READ_LIMA_DATA_NETCDF_CASE
USE MODI_AER2LIMA
USE MODI_CH_AER_EQM_INIT_n
!UPG*PT
USE MODI_READ_VER_GRID
USE MODI_SECOND_MNH
USE MODI_SET_REF
USE MODI_UPDATE_METRICS
USE MODI_VER_DYN
USE MODI_VER_PREP_GRIBEX_CASE
USE MODI_VER_PREP_MESONH_CASE
USE MODI_VER_PREP_NETCDF_CASE
USE MODI_VERSION
USE MODI_VER_THERMO
USE MODI_WRITE_DESFM_n
USE MODI_WRITE_LFIFM_n
!
USE MODN_CONF,             ONLY: JPHEXT , NHALO
USE MODN_CONFZ
!
IMPLICIT NONE
!
!*       0.1   Declaration of local variables
!              ------------------------------
!
CHARACTER(LEN=28)              :: YATMFILE    ! name of the Atmospheric file
CHARACTER(LEN=6)               :: YATMFILETYPE! type of the Atmospheric file
CHARACTER(LEN=28)              :: YCHEMFILE    ! name of the Chemical file
CHARACTER(LEN=6)               :: YCHEMFILETYPE! type of the Chemical file
!UP*PT
!CHARACTER(LEN=28)              :: YCAMSFILE    ! name of the input CAMS file
!CHARACTER(LEN=6)               :: YCAMSFILETYPE! type of the input CAMS file
CHARACTER(LEN=28)              :: YLIMAFILE    ! name of the input MACC file
CHARACTER(LEN=6)               :: YLIMAFILETYPE! type of the input MACC file
!UP*PT
CHARACTER(LEN=28)              :: YSURFFILE    ! name of the Surface file
CHARACTER(LEN=6)               :: YSURFFILETYPE! type of the Surface file
CHARACTER(LEN=28)              :: YPGDFILE    ! name of the physiographic data
!                                             ! file
!
CHARACTER(LEN=28)              :: YDAD_NAME   ! true name of the atmospheric file
!
!* other variables
!
REAL,DIMENSION(:,:,:), ALLOCATABLE:: ZJ       ! Jacobian
!
!* file management variables and counters
!
INTEGER                           :: ILUOUT0  ! logical unit for listing file
INTEGER                           :: IPRE_REAL1 ! logical unit for namelist file
INTEGER                           :: IRESP    ! return code in FM routines
LOGICAL                           :: GFOUND   ! Return code when searching namelist
INTEGER                           :: NIU,NJU,NKU   ! Upper bounds in x,y,z directions
!
REAL :: ZSTART, ZEND, ZTIME1, ZTIME2, ZTOT, ZALL ! for computing time analysis
REAL :: ZMISC, ZREAD, ZHORI, ZPREP, ZSURF, ZTHERMO, ZDYN, ZDIAG, ZWRITE
REAL :: ZDG                                      ! diagnostics time in routines
INTEGER                           :: IINFO_ll    ! return code of // routines
! Namelist model variables
CHARACTER(LEN=5)                  :: CPRESOPT
INTEGER                           :: NITR
LOGICAL                           :: LRES
REAL                              :: XRES
LOGICAL                           :: LSHIFT      ! flag to perform vertical shift or not.
LOGICAL                           :: LDUMMY_REAL ! flag to read and interpolate
                                                 !dummy fields from GRIBex file
INTEGER                           :: JRR      ! loop counter for moist var.
LOGICAL  :: LUSECHAQ
LOGICAL  :: LUSECHIC
LOGICAL  :: LUSECHEM
INTEGER :: JN
!
TYPE(TFILEDATA),POINTER :: TZATMFILE => NULL()
TYPE(TFILEDATA),POINTER :: TZPRE_REAL1FILE => NULL()
!
!
!*       0.3   Declaration of namelists
!              ------------------------
!
NAMELIST/NAM_REAL_CONF/ NVERB, CEQNSYS, CPRESOPT, LSHIFT, LDUMMY_REAL, &
                        LRES, XRES, NITR,LCOUPLING, NHALO , JPHEXT
! Filtering and balancing of the large-scale and radar tropical cyclone
NAMELIST/NAM_HURR_CONF/ LFILTERING, CFILTERING,   &
XLAMBDA, NK, XLATGUESS, XLONGUESS, XBOXWIND, XRADGUESS, NPHIL, NDIAG_FILT,   &
NLEVELR0,LBOGUSSING,               & 
XLATBOG, XLONBOG, XVTMAXSURF, XRADWINDSURF,       &
XMAX, XC, XRHO_Z, XRHO_ZZ, XB_0, XBETA_Z, XBETA_ZZ,&
XANGCONV0, XANGCONV1000, XANGCONV2000,            &
                        CDADATMFILE, CDADBOGFILE
 NAMELIST/NAM_AERO_CONF/ LORILAM, LINITPM, LDUST, XINIRADIUSI, XINIRADIUSJ,&
                         XINISIGI, XINISIGJ, XN0IMIN, XN0JMIN, CRGUNIT, CRGUNITD,&
                         LSALT, CRGUNITS, NMODE_DST, XINISIG, XINIRADIUS, XN0MIN,&
                         XINISIG_SLT, XINIRADIUS_SLT, XN0MIN_SLT, NMODE_SLT, &
                         LDSTCAMS, LSLTCAMS,CACTCCN,CCLOUD, NMOD_IFN, NMOD_CCN, LAERINIT

NAMELIST/NAM_CH_CONF/ LUSECHAQ,LUSECHIC,LUSECHEM
!
NAMELIST/NAM_IBM_LSF/ LIBM_LSF, CIBM_TYPE, NIBM_SMOOTH, XIBM_SMOOTH
!
! name of dad of input FM file
INTEGER                :: II, IJ, IGRID, ILENGTH
CHARACTER (LEN=100)    :: HCOMMENT
TYPE(LIST_ll), POINTER :: TZFIELDS_ll=>NULL()   ! list of fields to exchange
!
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZLBXRHO, ZLBYRHO
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZLBXZZ, ZLBYZZ
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZLBXPABST, ZLBYPABST
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZLBXRM, ZLBYRM
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZLBXTHM, ZLBYTHM
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZLBXSVM, ZLBYSVM
!
INTEGER :: ILBX,ILBY,IIB,IJB,IIE,IJE
LOGICAL :: GAERINIT
!-------------------------------------------------------------------------------
!
CALL MPPDB_INIT()
!
CALL GOTO_MODEL(1,ONOFIELDLIST=.TRUE.)
!
ZDIAG = 0.
CALL SECOND_MNH (ZSTART)
!
ZHORI  = 0.
ZSURF  = 0.
ZTIME1 = ZSTART
!
!*       1.    SET DEFAULT VALUES
!              ------------------
!
CALL VERSION
CPROGRAM='REAL  '
!
CALL ALLOC_FIELD_SCALARS()
CALL TBUCONF_ASSOCIATE()
CALL LES_ASSOCIATE()
CALL DEFAULT_DESFM_n(1)
NRR=1
IDX_RVT = 1
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
!UPG*PT
!                                   ,YCAMSFILE,YCAMSFILETYPE)
                                   ,YLIMAFILE,YLIMAFILETYPE)
!UPG*PT
ILUOUT0 = TLUOUT0%NLU
TLUOUT => TLUOUT0
!
IF (YATMFILETYPE=='MESONH') THEN
  LSHIFT = .FALSE.
ELSE IF (YATMFILETYPE=='GRIBEX') THEN
  LSHIFT = .TRUE.
ELSE
  LSHIFT = .TRUE.
  WRITE(ILUOUT0,FMT=*) 'HATMFILETYPE WAS SET TO: '//TRIM(YATMFILETYPE)
  WRITE(ILUOUT0,FMT=*) 'ONLY TWO VALUES POSSIBLE FOR HATMFILETYPE:'
  WRITE(ILUOUT0,FMT=*) 'EITHER MESONH OR GRIBEX'
  WRITE(ILUOUT0,FMT=*) '-> JOB ABORTED'
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_REAL_CASE','')
END IF
!                             
LCPL_AROME=.FALSE.
LCOUPLING=.FALSE.
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
IPRE_REAL1 = TZPRE_REAL1FILE%NLU
!
CALL INIT_NMLVAR
CALL POSNAM( TZPRE_REAL1FILE, 'NAM_REAL_CONF', GFOUND )
IF (GFOUND) READ(IPRE_REAL1,NAM_REAL_CONF)
CALL PARAM_LIMA_INIT(CPROGRAM, TZPRE_REAL1FILE, .FALSE., ILUOUT0, .FALSE., .TRUE., .FALSE., 0)
!
CALL INI_FIELD_LIST()
!
CALL INI_FIELD_SCALARS()
!
!*       4.2   reading of values of some configuration variables in namelist
!
!
!JUAN REALZ from prep_surfex
!
IF (YATMFILETYPE == 'GRIBEX') THEN
!
!*       4.1  Vertical Spatial grid 
!
CALL INIT_NMLVAR()
CALL READ_VER_GRID(TZPRE_REAL1FILE)
!
CALL IO_Field_read(TPGDFILE,'IMAX',NIMAX)
CALL IO_Field_read(TPGDFILE,'JMAX',NJMAX)
!
NIMAX_ll=NIMAX   !! _ll variables are global variables
NJMAX_ll=NJMAX   !! but the old names are kept in PRE_IDEA1.nam file
!
CALL SET_JP_ll(JPMODELMAX,JPHEXT,JPVEXT,JPHEXT)
CALL SET_DAD0_ll()
!JUAN 4/04/2014 correction for PREP_REAL_CASE on Gribex files 
!CALL SET_DIM_ll(NIMAX_ll, NJMAX_ll, 128)
CALL SET_DIM_ll(NIMAX_ll, NJMAX_ll, NKMAX)
CALL SET_LBX_ll('OPEN',1)
CALL SET_LBY_ll('OPEN', 1)
CALL SET_XRATIO_ll(1, 1)
CALL SET_YRATIO_ll(1, 1)
CALL SET_XOR_ll(1, 1)
CALL SET_XEND_ll(NIMAX_ll+2*JPHEXT, 1)
CALL SET_YOR_ll(1, 1)
CALL SET_YEND_ll(NJMAX_ll+2*JPHEXT, 1)
CALL SET_DAD_ll(0, 1)
!JUANZ
!CALL INI_PARA_ll(IINFO_ll)
CALL INI_PARAZ_ll(IINFO_ll)
!JUANZ

!
! sizes of arrays of the extended sub-domain
!
CALL GET_DIM_PHYS_ll('B',NIMAX,NJMAX)
!!$CALL GET_DIM_EXT_ll('B',NIU,NJU)
!!$CALL GET_INDICE_ll(NIB,NJB,NIE,NJE)
!!$CALL GET_OR_ll('B',IXOR,IYOR)
ENDIF
!JUAN REALZ
!
LDUMMY_REAL= .FALSE.
LFILTERING= .FALSE.
CFILTERING= 'UVT  '
XLATGUESS= XUNDEF ; XLONGUESS= XUNDEF ; XBOXWIND=XUNDEF; XRADGUESS= XUNDEF
NK=50 ; XLAMBDA=0.2 ; NPHIL=24
NLEVELR0=15
NDIAG_FILT=-1
LBOGUSSING= .FALSE.
XLATBOG= XUNDEF ; XLONBOG= XUNDEF
XVTMAXSURF= XUNDEF ; XRADWINDSURF= XUNDEF
XMAX=16000. ; XC=0.7 ; XRHO_Z=-0.3 ; XRHO_ZZ=0.9
XB_0=1.65 ; XBETA_Z=-0.5 ; XBETA_ZZ=0.35
XANGCONV0=0. ; XANGCONV1000=0. ; XANGCONV2000=0.
CDADATMFILE=' ' ; CDADBOGFILE=' '
!
CALL INIT_NMLVAR
CALL POSNAM( TZPRE_REAL1FILE, 'NAM_REAL_CONF', GFOUND )
IF (GFOUND) READ(IPRE_REAL1,NAM_REAL_CONF)
CALL POSNAM( TZPRE_REAL1FILE, 'NAM_HURR_CONF', GFOUND )
IF (GFOUND) READ(IPRE_REAL1,NAM_HURR_CONF)
CALL POSNAM( TZPRE_REAL1FILE, 'NAM_CH_CONF', GFOUND )
IF (GFOUND) READ(UNIT=IPRE_REAL1,NML=NAM_CH_CONF)
CALL UPDATE_MODD_FROM_NMLVAR
CALL POSNAM( TZPRE_REAL1FILE, 'NAM_AERO_CONF', GFOUND )
IF (GFOUND) READ(IPRE_REAL1,NAM_AERO_CONF)
CALL POSNAM( TZPRE_REAL1FILE, 'NAM_CONFZ', GFOUND )
IF (GFOUND) READ(UNIT=IPRE_REAL1,NML=NAM_CONFZ)
CALL POSNAM( TZPRE_REAL1FILE, 'NAM_IBM_LSF' , GFOUND )
IF (GFOUND) READ(UNIT=IPRE_REAL1,NML=NAM_IBM_LSF)
!
GAERINIT = LAERINIT

! Sea salt
CALL INIT_SALT
!
!*       4.3   set soil scheme to ISBA for initialization from GRIB
!
IF (YATMFILETYPE=='GRIBEX') THEN
  CLBCX(:) ='OPEN'
  CLBCY(:) ='OPEN'
END IF
!
CALL SECOND_MNH(ZTIME2)
ZMISC = ZTIME2 - ZTIME1
!-------------------------------------------------------------------------------
!
!*       5.    READING OF THE INPUT DATA
!              -------------------------
!
ZTIME1 = ZTIME2
!
IF (YATMFILETYPE=='MESONH') THEN
  CALL READ_ALL_DATA_MESONH_CASE(TZPRE_REAL1FILE,YATMFILE,TPGDFILE,YDAD_NAME)
ELSE IF (YATMFILETYPE=='GRIBEX') THEN
  IF(LEN_TRIM(YCHEMFILE)>0 .AND. YCHEMFILETYPE=='GRIBEX')THEN
    CALL READ_ALL_DATA_GRIB_CASE('ATM1',TZPRE_REAL1FILE,YATMFILE,TPGDFILE,ZHORI,NVERB,LDUMMY_REAL)
  ELSE
    CALL READ_ALL_DATA_GRIB_CASE('ATM0',TZPRE_REAL1FILE,YATMFILE,TPGDFILE,ZHORI,NVERB,LDUMMY_REAL)
  END IF
!
  YDAD_NAME=' '
END IF
LAERINIT = GAERINIT 
!
IF (NIMAX==1 .AND. NJMAX==1) THEN
  L1D=.TRUE.
  L2D=.FALSE.
ELSE IF (NJMAX==1) THEN
  L1D=.FALSE.
  L2D=.TRUE.
ELSE
  L1D=.FALSE.
  L2D=.FALSE.
END IF
!
! UPG*PT
!*       5.1   reading of the input chemical data
!
!IF(LEN_TRIM(YCHEMFILE)>0)THEN
!  ! read again Nam_aero_conf
!  CALL POSNAM( TZPRE_REAL1FILE, 'NAM_AERO_CONF', GFOUND )
!  IF (GFOUND) READ(IPRE_REAL1,NAM_AERO_CONF)
!  IF(YCHEMFILETYPE=='GRIBEX') &
!  CALL READ_ALL_DATA_GRIB_CASE('CHEM',TZPRE_REAL1FILE,YCHEMFILE,TPGDFILE,ZHORI,NVERB,LDUMMY_REAL)
!  IF (YCHEMFILETYPE=='NETCDF') &
!  CALL READ_CHEM_DATA_NETCDF_CASE(TZPRE_REAL1FILE,YCHEMFILE,TPGDFILE,ZHORI,NVERB,LDUMMY_REAL)
!END IF
!
!*       5.2   reading the input CAMS data
!
!IF(LEN_TRIM(YCAMSFILE)>0)THEN
!   IF(YCAMSFILETYPE=='NETCDF') THEN
!      CALL READ_CAMS_DATA_NETCDF_CASE(TZPRE_REAL1FILE,YCAMSFILE,TPGDFILE,ZHORI,NVERB,LDUMMY_REAL)
!   ELSE
!      CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_REAL_CASE','CANNOT READ CAMS GRIB FILES YET')
!   END IF
!END IF
!*       5.1   reading CAMS or MACC files for init LIMA
!
IF(LEN_TRIM(YLIMAFILE)>0)THEN
  IF(YLIMAFILETYPE=='NETCDF') THEN
    CALL READ_LIMA_DATA_NETCDF_CASE(TZPRE_REAL1FILE,YLIMAFILE,TPGDFILE,ZHORI,NVERB,LDUMMY_REAL)
  ELSE
    WRITE(ILUOUT0,FMT=*)
    !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_REAL_CASE','Pb in MACC/CAMS file')
    STOP
  END IF
END IF
!
!*       5.2   reading of the input chemical data + dusts + salts if needed
!
IF(LEN_TRIM(YCHEMFILE)>0)THEN
  ! read again Nam_aero_conf
  CALL POSNAM( TZPRE_REAL1FILE, 'NAM_AERO_CONF', GFOUND )
  IF (GFOUND) READ(IPRE_REAL1,NAM_AERO_CONF)
  IF(YCHEMFILETYPE=='GRIBEX') &
  CALL READ_ALL_DATA_GRIB_CASE('CHEM',TZPRE_REAL1FILE,YCHEMFILE,TPGDFILE,ZHORI,NVERB,LDUMMY_REAL)
  IF (YCHEMFILETYPE=='MOZART') &
  CALL READ_CHEM_DATA_MOZART_CASE(TZPRE_REAL1FILE,YCHEMFILE,TPGDFILE,ZHORI,NVERB,LDUMMY_REAL)
  IF (YCHEMFILETYPE=='CAMSEU') &
  CALL READ_CHEM_DATA_CAMS_CASE(TZPRE_REAL1FILE,YCHEMFILE,TPGDFILE,ZHORI,NVERB, &
                                LDUMMY_REAL,LUSECHEM)
END IF

!UPG*PT
!
CALL IO_File_close(TZPRE_REAL1FILE)
!
CALL SECOND_MNH(ZTIME2)
ZREAD = ZTIME2 - ZTIME1 - ZHORI
!-------------------------------------------------------------------------------
!
CALL IO_File_add2list(TINIFILE,CINIFILE,'MNH','WRITE',KLFITYPE=1,KLFIVERB=NVERB)
CALL IO_File_open(TINIFILE)
!
ZTIME1=ZTIME2
!
!*       6.    CONFIGURATION VARIABLES
!              -----------------------
!
!*       6.1   imposed values of some other configuration variables
!
CDCONV='NONE'
CSCONV='NONE'
CRAD='NONE'
CCONF='START'
NRIMX=6
NRIMY=6
LHORELAX_UVWTH=.TRUE.
LHORELAX_RV=LUSERV
LHORELAX_RC=LUSERC
LHORELAX_RR=LUSERR
LHORELAX_RI=LUSERI
LHORELAX_RS=LUSERS
LHORELAX_RG=LUSERG
LHORELAX_RH=LUSERH
LHORELAX_SV(:)=.FALSE.
LHORELAX_SVC2R2 = (NSV_C2R2 > 0)
LHORELAX_SVC1R3 = (NSV_C1R3 > 0)
LHORELAX_SVLIMA = (NSV_LIMA > 0)
LHORELAX_SVELEC = (NSV_ELEC > 0)
LHORELAX_SVCHEM = (NSV_CHEM > 0)
LHORELAX_SVCHIC = (NSV_CHIC > 0)
LHORELAX_SVDST  = (NSV_DST > 0)
LHORELAX_SVSLT  = (NSV_SLT > 0)
LHORELAX_SVAER  = (NSV_AER > 0)
LHORELAX_SVPP   = (NSV_PP > 0)
#ifdef MNH_FOREFIRE
LHORELAX_SVFF   = (NSV_FF > 0)
#endif
LHORELAX_SVCS   = (NSV_CS > 0)

LHORELAX_SVLG   = .FALSE.
LHORELAX_SV(1:NSV)=.TRUE.
IF ( CTURB /= 'NONE') THEN
  LHORELAX_TKE = .TRUE.
ELSE
  LHORELAX_TKE = .FALSE.
END IF
!
!
CSTORAGE_TYPE='TT'
!-------------------------------------------------------------------------------
!
!*       8.    COMPUTATION OF GEOMETRIC VARIABLES
!              ----------------------------------
!
ZTIME1 = ZTIME2
!
ALLOCATE(XMAP(SIZE(XXHAT),SIZE(XYHAT)))
ALLOCATE(XLAT(SIZE(XXHAT),SIZE(XYHAT)))
ALLOCATE(XLON(SIZE(XXHAT),SIZE(XYHAT)))
ALLOCATE(XDXHAT(SIZE(XXHAT)))
ALLOCATE(XDYHAT(SIZE(XYHAT)))
ALLOCATE(XZZ(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)))
ALLOCATE(ZJ(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)))
!
IF (LCARTESIAN) THEN
  CALL SM_GRIDCART(XXHAT,XYHAT,XZHAT,XZS,LSLEVE,XLEN1,XLEN2,XZSMT,XDXHAT,XDYHAT,XZZ,ZJ)
  XMAP=1.
ELSE
  CALL SM_GRIDPROJ( XXHAT, XYHAT, XZHAT, XXHATM, XYHATM, XZS,      &
                    LSLEVE, XLEN1, XLEN2, XZSMT, XLATORI, XLONORI, &
                    XMAP, XLAT, XLON, XDXHAT, XDYHAT, XZZ, ZJ      )
END IF
!
CALL MPPDB_CHECK2D(XZS,"prep_real_case8:XZS",PRECISION)
CALL MPPDB_CHECK2D(XMAP,"prep_real_case8:XMAP",PRECISION)
CALL MPPDB_CHECK2D(XLAT,"prep_real_case8:XLAT",PRECISION)
CALL MPPDB_CHECK2D(XLON,"prep_real_case8:XLON",PRECISION)
CALL MPPDB_CHECK3D(XZZ,"prep_real_case8:XZZ",PRECISION)
CALL MPPDB_CHECK3D(ZJ,"prep_real_case8:ZJ",PRECISION)
!
ALLOCATE(XDXX(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)))
ALLOCATE(XDYY(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)))
ALLOCATE(XDZX(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)))
ALLOCATE(XDZY(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)))
ALLOCATE(XDZZ(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)))
!
!20131024 add update halo
!=> corrects on PDXX calculation in metrics and XDXX !!
CALL ADD3DFIELD_ll( TZFIELDS_ll, XZZ, 'PREP_REAL_CASE::XZZ' )
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
CALL METRICS(XMAP,XDXHAT,XDYHAT,XZZ,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
CALL MPPDB_CHECK3D(XDXX,"prc8-beforeupdate_metrics:PDXX",PRECISION)
CALL MPPDB_CHECK3D(XDYY,"prc8-beforeupdate_metrics:PDYY",PRECISION)
CALL MPPDB_CHECK3D(XDZX,"prc8-beforeupdate_metrics:PDZX",PRECISION)
CALL MPPDB_CHECK3D(XDZY,"prc8-beforeupdate_metrics:PDZY",PRECISION)
!
CALL UPDATE_METRICS(CLBCX,CLBCY,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
!20131112 add update_halo for XDYY and XDZY!!
CALL ADD3DFIELD_ll( TZFIELDS_ll, XDXX, 'PREP_REAL_CASE::XDXX' )
CALL ADD3DFIELD_ll( TZFIELDS_ll, XDZX, 'PREP_REAL_CASE::XDZX' )
CALL ADD3DFIELD_ll( TZFIELDS_ll, XDYY, 'PREP_REAL_CASE::XDYY' )
CALL ADD3DFIELD_ll( TZFIELDS_ll, XDZY, 'PREP_REAL_CASE::XDZY' )
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)

!CALL EXTRAPOL('W',XDXX,XDZX)
!CALL EXTRAPOL('S',XDYY,XDZY)

CALL SECOND_MNH(ZTIME2)

ZMISC = ZMISC + ZTIME2 - ZTIME1
!-------------------------------------------------------------------------------
!
!*       9.    PREPARATION OF THE VERTICAL SHIFT AND INTERPOLATION
!              ---------------------------------------------------
!
ZTIME1 = ZTIME2
!
IF (YATMFILETYPE=='GRIBEX') THEN
  CALL VER_PREP_GRIBEX_CASE('ATM ',ZDG)
ELSE IF (YATMFILETYPE=='MESONH') THEN
  CALL VER_PREP_MESONH_CASE(ZDG)
END IF
!
IF (LEN_TRIM(YCHEMFILE)>0 .AND. YCHEMFILETYPE=='GRIBEX') THEN
  CALL VER_PREP_GRIBEX_CASE('CHEM',ZDG)
END IF
!UPG*PT
!IF ((LEN_TRIM(YCHEMFILE)>0 .AND. YCHEMFILETYPE=='NETCDF') .OR. &
!   (LEN_TRIM(YCAMSFILE)>0 .AND. YCAMSFILETYPE=='NETCDF')) THEN
!  CALL VER_PREP_NETCDF_CASE(ZDG)
!END IF
IF (LEN_TRIM(YCHEMFILE)>0 .AND. ((YCHEMFILETYPE=='MOZART').OR. &
                                 (YCHEMFILETYPE=='CAMSEU'))) THEN
  CALL VER_PREP_NETCDF_CASE(ZDG,XSV_LS)

  DEALLOCATE(XSV_LS)
END IF
!
IF (LEN_TRIM(YLIMAFILE)>0 .AND. YLIMAFILETYPE=='NETCDF') THEN
  CALL VER_PREP_NETCDF_CASE(ZDG,XSV_LS_LIMA)
  DEALLOCATE(XSV_LS_LIMA)
END IF
!UPG*PT
!
CALL SECOND_MNH(ZTIME2)
ZPREP = ZTIME2 - ZTIME1 - ZDG
ZDIAG = ZDIAG + ZDG
!-------------------------------------------------------------------------------
!
!*      10.    VERTICAL INTERPOLATION OF ALL THERMODYNAMICAL VARIABLES
!              -------------------------------------------------------
!
ZTIME1 = ZTIME2
!
ALLOCATE(XPSURF(SIZE(XXHAT),SIZE(XYHAT)))
!
CALL EXTRAPOL('E',XEXNTOP2D)
IF (YATMFILETYPE=='GRIBEX') THEN
  CALL VER_THERMO(TINIFILE,LSHIFT,XTHV_MX,XR_MX,XZS_LS,XZSMT_LS,XZMASS_MX,XZFLUX_MX,XPMHP_MX,ZJ, &
                  XDXX,XDYY,XEXNTOP2D,XPSURF,ZDG                               )
ELSE IF (YATMFILETYPE=='MESONH') THEN
  CALL VER_THERMO(TINIFILE,LSHIFT,XTHV_MX,XR_MX,XZS_LS,XZSMT_LS,XZMASS_MX,XZFLUX_MX,XPMHP_MX,ZJ, &
                  XDXX,XDYY,XEXNTOP2D,XPSURF,ZDG,                              &
                  XLSTH_MX,XLSRV_MX                                            )
END IF
!
CALL SECOND_MNH(ZTIME2)
ZTHERMO = ZTIME2 - ZTIME1 - ZDG
ZDIAG = ZDIAG + ZDG
!-------------------------------------------------------------------------------
!
!*      12.    VERTICAL INTERPOLATION OF DYNAMICAL VARIABLES
!              ---------------------------------------------
!
ZTIME1 = ZTIME2
IF (YATMFILETYPE=='GRIBEX') THEN
  CALL VER_DYN(LSHIFT,XU_MX,XV_MX,XW_MX,XRHOD_MX,XZFLUX_MX,XZMASS_MX,XZS_LS,     &
               XDXX,XDYY,XDZZ,XDZX,XDZY,ZJ,YATMFILETYPE                          )
ELSE IF (YATMFILETYPE=='MESONH') THEN
  CALL VER_DYN(LSHIFT,XU_MX,XV_MX,XW_MX,XRHOD_MX,XZFLUX_MX,XZMASS_MX,XZS_LS,     &
               XDXX,XDYY,XDZZ,XDZX,XDZY,ZJ,YATMFILETYPE,                         &
               XLSU_MX,XLSV_MX,XLSW_MX                                           )
END IF
!
!
IF (ALLOCATED(XTHV_MX)) DEALLOCATE(XTHV_MX)
IF (ALLOCATED(XR_MX)) DEALLOCATE(XR_MX)
IF (ALLOCATED(XPMHP_MX)) DEALLOCATE(XPMHP_MX)
IF (ALLOCATED(XU_MX)) DEALLOCATE(XU_MX)
IF (ALLOCATED(XV_MX)) DEALLOCATE(XV_MX)
IF (ALLOCATED(XW_MX)) DEALLOCATE(XW_MX)
IF (ALLOCATED(XLSTH_MX)) DEALLOCATE(XLSTH_MX)
IF (ALLOCATED(XLSRV_MX)) DEALLOCATE(XLSRV_MX)
IF (ALLOCATED(XLSU_MX)) DEALLOCATE(XLSU_MX)
IF (ALLOCATED(XLSV_MX)) DEALLOCATE(XLSV_MX)
IF (ALLOCATED(XLSW_MX)) DEALLOCATE(XLSW_MX)
IF (ALLOCATED(XZFLUX_MX)) DEALLOCATE(XZFLUX_MX)
IF (ALLOCATED(XZMASS_MX)) DEALLOCATE(XZMASS_MX)
IF (ALLOCATED(XRHOD_MX)) DEALLOCATE(XRHOD_MX)
IF (ALLOCATED(XEXNTOP2D)) DEALLOCATE(XEXNTOP2D)
IF (ALLOCATED(XZS_LS)) DEALLOCATE(XZS_LS)
IF (ALLOCATED(XZSMT_LS)) DEALLOCATE(XZSMT_LS)
!
!-------------------------------------------------------------------------------
!
!*      13.    ANELASTIC CORRECTION
!              --------------------
!
CALL PRESSURE_IN_PREP(XDXX,XDYY,XDZX,XDZY,XDZZ)
!
CALL SECOND_MNH(ZTIME2)
ZDYN = ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*      14.    INITIALIZATION OF THE REMAINING PROGNOSTIC VARIABLES (COPIES)
!              -------------------------------------------------------------
!
ZTIME1 = ZTIME2
!
IF(LEN_TRIM(YCHEMFILE)>0 .AND. YCHEMFILETYPE=='MESONH')THEN
  CALL INI_PROG_VAR(XTKE_MX,XSV_MX,YCHEMFILE)
  LHORELAX_SVCHEM = (NSV_CHEM > 0)
  LHORELAX_SVCHIC = (NSV_CHIC > 0)
  LHORELAX_SVDST  = (NSV_DST > 0)
  LHORELAX_SVSLT  = (NSV_SLT > 0)
  LHORELAX_SVAER  = (NSV_AER > 0)
ELSE
!
!UPG*PT
!IF (LEN_TRIM(YCAMSFILE)>0 .AND. YCAMSFILETYPE=='NETCDF') THEN
IF (LEN_TRIM(YLIMAFILE)>0 .AND. YLIMAFILETYPE=='NETCDF') THEN
!UPG*PT
  CALL LIMA_MIXRAT_TO_NCONC(XPABST, XTHT, XRT(:,:,:,1), XSV_MX)
END IF
!
  CALL INI_PROG_VAR(XTKE_MX,XSV_MX)
END IF
!

! Initialization of ORILAM variables
IF (LORILAM) THEN
  IF (.NOT.(ASSOCIATED(XN3D)))      ALLOCATE(XN3D(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),JPMODE))
  IF (.NOT.(ASSOCIATED(XRG3D)))     ALLOCATE(XRG3D(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),JPMODE))
  IF (.NOT.(ASSOCIATED(XSIG3D)))    ALLOCATE(XSIG3D(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),JPMODE))
  IF (.NOT.(ASSOCIATED(XRHOP3D)))   ALLOCATE(XRHOP3D(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),JPMODE))
  IF (.NOT.(ASSOCIATED(XM3D)))      ALLOCATE(XM3D(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),JPMODE*3))
  IF (.NOT.(ASSOCIATED(XCTOTA3D))) &
    ALLOCATE(XCTOTA3D(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NSP+NCARB+NSOA,JPMODE))

  CALL CH_AER_EQM_INIT_n(XSVT(:,:,:,NSV_CHEMBEG:NSV_CHEMEND),&
                         XSVT(:,:,:,NSV_AERBEG:NSV_AEREND),&
                         XM3D,XRHOP3D,XSIG3D,&
                         XRG3D,XN3D, XRHODREF, XCTOTA3D)
END IF
!
! Initialization LIMA variables by ORILAM
IF (CCLOUD == 'LIMA' .AND. ((LORILAM).OR.(LDUST).OR.(LSALT))) THEN

    ! Init LIMA by ORILAM
    CALL AER2LIMA(XSVT, XRHODREF, XRT(:,:,:,1), XPABST, XTHT,XZZ)

    ! Init LB LIMA by ORILAM
    ALLOCATE(ZLBXRHO(SIZE(XLBXSVM,1), SIZE(XLBXSVM,2), SIZE(XLBXSVM,3)))
    ALLOCATE(ZLBYRHO(SIZE(XLBYSVM,1), SIZE(XLBYSVM,2), SIZE(XLBYSVM,3)))
    ALLOCATE(ZLBXPABST(SIZE(XLBXSVM,1), SIZE(XLBXSVM,2), SIZE(XLBXSVM,3)))
    ALLOCATE(ZLBYPABST(SIZE(XLBYSVM,1), SIZE(XLBYSVM,2), SIZE(XLBYSVM,3)))
    ALLOCATE(ZLBXTHM(SIZE(XLBXSVM,1), SIZE(XLBXSVM,2), SIZE(XLBXSVM,3)))
    ALLOCATE(ZLBYTHM(SIZE(XLBYSVM,1), SIZE(XLBYSVM,2), SIZE(XLBYSVM,3)))
    ALLOCATE(ZLBXZZ(SIZE(XLBXSVM,1), SIZE(XLBXSVM,2), SIZE(XLBXSVM,3)))
    ALLOCATE(ZLBYZZ(SIZE(XLBYSVM,1), SIZE(XLBYSVM,2), SIZE(XLBYSVM,3)))
    ALLOCATE(ZLBXRM(SIZE(XLBXSVM,1), SIZE(XLBXSVM,2), SIZE(XLBXSVM,3)))
    ALLOCATE(ZLBYRM(SIZE(XLBYSVM,1), SIZE(XLBYSVM,2), SIZE(XLBYSVM,3)))
    ALLOCATE(ZLBXSVM(SIZE(XLBXSVM,1), SIZE(XLBXSVM,2), SIZE(XLBXSVM,3), SIZE(XLBXSVM,4)))
    ALLOCATE(ZLBYSVM(SIZE(XLBYSVM,1), SIZE(XLBYSVM,2), SIZE(XLBYSVM,3), SIZE(XLBXSVM,4)))

    ILBX=SIZE(XLBXSVM,1)/2-JPHEXT
    ILBY=SIZE(XLBYSVM,2)/2-JPHEXT

    CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)

    ZLBXRHO(1:ILBX+1,:,:)        = XRHODREF(IIB-1:IIB-1+ILBX,:,:)
    ZLBXRHO(ILBX+2:2*ILBX+2,:,:) = XRHODREF(IIE+1-ILBX:IIE+1,:,:)
    ZLBYRHO(:,1:ILBY+1,:)        = XRHODREF(:,IJB-1:IJB-1+ILBY,:)
    ZLBYRHO(:,ILBY+2:2*ILBY+2,:) = XRHODREF(:,IJE+1-ILBY:IJE+1,:)
    ZLBXPABST(1:ILBX+1,:,:)        = XPABST(IIB-1:IIB-1+ILBX,:,:)
    ZLBXPABST(ILBX+2:2*ILBX+2,:,:) = XPABST(IIE+1-ILBX:IIE+1,:,:)
    ZLBYPABST(:,1:ILBY+1,:)        = XPABST(:,IJB-1:IJB-1+ILBY,:)
    ZLBYPABST(:,ILBY+2:2*ILBY+2,:) = XPABST(:,IJE+1-ILBY:IJE+1,:)
    ZLBXTHM(1:ILBX+1,:,:)          = XTHT(IIB-1:IIB-1+ILBX,:,:)
    ZLBXTHM(ILBX+2:2*ILBX+2,:,:)   = XTHT(IIE+1-ILBX:IIE+1,:,:)
    ZLBYTHM(:,1:ILBY+1,:)          = XTHT(:,IJB-1:IJB-1+ILBY,:)
    ZLBYTHM(:,ILBY+2:2*ILBY+2,:)   = XTHT(:,IJE+1-ILBY:IJE+1,:)
    ZLBXZZ(1:ILBX+1,:,:)        = XZZ(IIB-1:IIB-1+ILBX,:,:)
    ZLBXZZ(ILBX+2:2*ILBX+2,:,:) = XZZ(IIE+1-ILBX:IIE+1,:,:)
    ZLBYZZ(:,1:ILBY+1,:)        = XZZ(:,IJB-1:IJB-1+ILBY,:)
    ZLBYZZ(:,ILBY+2:2*ILBY+2,:) = XZZ(:,IJE+1-ILBY:IJE+1,:)
    ZLBXSVM(1:ILBX+1,:,:,:)        = XSVT(IIB-1:IIB-1+ILBX,:,:,:)
    ZLBXSVM(ILBX+2:2*ILBX+2,:,:,:) = XSVT(IIE+1-ILBX:IIE+1,:,:,:)
    ZLBYSVM(:,1:ILBY+1,:,:)        = XSVT(:,IJB-1:IJB-1+ILBY,:,:)
    ZLBYSVM(:,ILBY+2:2*ILBY+2,:,:) = XSVT(:,IJE+1-ILBY:IJE+1,:,:)
    ZLBXRM(1:ILBX+1,:,:)          = XRT(IIB-1:IIB-1+ILBX,:,:,1)
    ZLBXRM(ILBX+2:2*ILBX+2,:,:)   = XRT(IIE+1-ILBX:IIE+1,:,:,1)
    ZLBYRM(:,1:ILBY+1,:)          = XRT(:,IJB-1:IJB-1+ILBY,:,1)
    ZLBYRM(:,ILBY+2:2*ILBY+2,:)   = XRT(:,IJE+1-ILBY:IJE+1,:,1)


    CALL AER2LIMA(ZLBXSVM, ZLBXRHO, ZLBXRM(:,:,:), ZLBXPABST, ZLBXTHM, ZLBXZZ)
    CALL AER2LIMA(ZLBYSVM, ZLBYRHO, ZLBYRM(:,:,:), ZLBYPABST, ZLBYTHM, ZLBYZZ)

    DEALLOCATE(ZLBXRHO)
    DEALLOCATE(ZLBYRHO)
    DEALLOCATE(ZLBXPABST)
    DEALLOCATE(ZLBYPABST)
    DEALLOCATE(ZLBXTHM)
    DEALLOCATE(ZLBYTHM)
    DEALLOCATE(ZLBXZZ)
    DEALLOCATE(ZLBYZZ)
    DEALLOCATE(ZLBXRM)
    DEALLOCATE(ZLBYRM)
    DEALLOCATE(ZLBXSVM)
    DEALLOCATE(ZLBYSVM)
END IF
!
IF (ALLOCATED(XSV_MX)) DEALLOCATE(XSV_MX)
IF (ALLOCATED(XTKE_MX)) DEALLOCATE(XTKE_MX)
!
CALL BOUNDARIES (                                                 &
          0.,CLBCX,CLBCY,NRR,NSV,1,                               &
          XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,   &
          XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,   &
          XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,   &
          XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,   &
          XRHODJ,XRHODREF,                                        &
          XUT, XVT, XWT, XTHT, XTKET, XRT, XSVT, XSRCT            )
!
CALL SECOND_MNH(ZTIME2)
ZMISC = ZMISC + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*      15.    Error on temperature during interpolations
!              ------------------------------------------
!
ZTIME1 = ZTIME2
!
IF (YATMFILETYPE=='GRIBEX' .AND. NVERB>1) THEN
  CALL ERROR_ON_TEMPERATURE(XT_LS,XPMASS_LS,XPABST,XPS_LS,XPSURF)
END IF
!
IF (YATMFILETYPE=='GRIBEX') THEN
  DEALLOCATE(XT_LS)
  DEALLOCATE(XPMASS_LS)
  DEALLOCATE(XPS_LS)
END IF
!
IF (ALLOCATED(XPSURF)) DEALLOCATE(XPSURF)
!
CALL SECOND_MNH(ZTIME2)
ZDIAG =  ZDIAG + ZTIME2 - ZTIME1
!-------------------------------------------------------------------------------
!
!*       16.    INITIALIZE LEVELSET FOR IBM
!              ---------------------------
!
IF (LIBM_LSF) THEN
  !
  IF (.NOT.LCARTESIAN) THEN
    CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_IDEAL_CASE','IBM can only be used with cartesian coordinates')
  ENDIF
  !
  CALL GET_DIM_EXT_ll('B',NIU,NJU)
  NKU=NKMAX+2*JPVEXT
  !
  ALLOCATE(XIBM_LS(NIU,NJU,NKU,4))
  !
  CALL IBM_INIT_LS(XIBM_LS)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      17.    WRITING OF THE MESO-NH FM-FILE
!              ------------------------------
!
ZTIME1 = ZTIME2
!
CSTORAGE_TYPE='TT'
IF (YATMFILETYPE=='GRIBEX') THEN
  CSURF = "EXTE"
  DO JRR=1,NRR
    IF (JRR==1) THEN
      LUSERV=.TRUE.
      IDX_RVT = JRR
    END IF
    IF (JRR==2) THEN
      LUSERC=.TRUE.
      IDX_RCT = JRR
    END IF
    IF (JRR==3) THEN
      LUSERR=.TRUE.
      IDX_RRT = JRR
    END IF
    IF (JRR==4) THEN
      LUSERI=.TRUE.
      IDX_RIT = JRR
    END IF
    IF (JRR==5) THEN
      LUSERS=.TRUE.
      IDX_RST = JRR
    END IF
    IF (JRR==6) THEN
      LUSERG=.TRUE.
      IDX_RGT = JRR
    END IF
    IF (JRR==7) THEN
      LUSERH=.TRUE.
      IDX_RHT = JRR
    END IF
  END DO
END IF
!
CALL WRITE_DESFM_n(1,TINIFILE)
CALL IO_Header_write(TINIFILE,HDAD_NAME=YDAD_NAME)
CALL WRITE_LFIFM_n(TINIFILE,YDAD_NAME)
! 
CALL SECOND_MNH(ZTIME2)
ZWRITE = ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*      18.    OROGRAPHIC and DUMMY PHYSIOGRAPHIC FIELDS
!              -----------------------------------------
!
!* reading in the PGD file
!
CALL MNHREAD_ZS_DUMMY_n(TPGDFILE)
!
!* writing in the output file
!
TOUTDATAFILE => TINIFILE
CALL MNHWRITE_ZS_DUMMY_n(TINIFILE)
!
CALL DEALLOCATE_MODEL1(3)
!
IF (YATMFILETYPE=='MESONH'.AND. YATMFILE/=YPGDFILE) THEN
  CALL IO_File_find_byname(TRIM(YATMFILE),TZATMFILE,IRESP)
  CALL IO_File_close(TZATMFILE)
END IF
!-------------------------------------------------------------------------------
!
!*      19.    INTERPOLATION OF SURFACE VARIABLES
!              ----------------------------------
!
IF (.NOT. LCOUPLING ) THEN
  ZTIME1 = ZTIME2
!
  IF (CSURF=="EXTE") THEN
    IF (YATMFILETYPE/='MESONH') THEN
      CALL SURFEX_ALLOC_LIST(1)
      YSURF_CUR => YSURF_LIST(1)
      CALL READ_ALL_NAMELISTS(YSURF_CUR,'MESONH','PRE',.FALSE.)
    ENDIF
    CALL GOTO_SURFEX(1)
    TFILE_SURFEX => TINIFILE
    CALL PREP_SURF_MNH(YSURFFILE,YSURFFILETYPE)
    NULLIFY(TFILE_SURFEX)
  ENDIF
!
  CALL SECOND_MNH(ZTIME2)
  ZSURF = ZSURF + ZTIME2 - ZTIME1
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      20.    EPILOGUE
!              --------
!
WRITE(ILUOUT0,*)
WRITE(ILUOUT0,*)
WRITE(ILUOUT0,*)
WRITE(ILUOUT0,*) '**************************************************'
WRITE(ILUOUT0,*) '* PREP_REAL_CASE: PREP_REAL_CASE ends correctly. *'
WRITE(ILUOUT0,*) '**************************************************'
WRITE(ILUOUT0,*)
!
!-------------------------------------------------------------------------------
!
CALL SECOND_MNH (ZEND)
!
ZTOT        = ZEND - ZSTART          ! for computing time analysis
!
ZALL = ZMISC + ZREAD + ZHORI + ZPREP + ZTHERMO + ZSURF + ZDYN + ZDIAG + ZWRITE 
!
WRITE(ILUOUT0,*)
WRITE(ILUOUT0,*) ' ------------------------------------------------------------ '
WRITE(ILUOUT0,*) '|                                                            |'
WRITE(ILUOUT0,*) '|          COMPUTING TIME ANALYSIS in PREP_REAL_CASE         |'
WRITE(ILUOUT0,*) '|                                                            |'
WRITE(ILUOUT0,*) '|------------------------------------------------------------|'
WRITE(ILUOUT0,*) '|                     |                   |                  |'
WRITE(ILUOUT0,*) '|    ROUTINE NAME     |     CPU-TIME      |   PERCENTAGE %   |'
WRITE(ILUOUT0,*) '|                     |                   |                  |'
WRITE(ILUOUT0,*) '|---------------------|-------------------|------------------|'
WRITE(ILUOUT0,*) '|                     |                   |                  |'
WRITE(UNIT=ILUOUT0,FMT=2) ZREAD, 100.*ZREAD/ZTOT
WRITE(UNIT=ILUOUT0,FMT=9) ZHORI, 100.*ZHORI/ZTOT
WRITE(UNIT=ILUOUT0,FMT=3) ZPREP, 100.*ZPREP/ZTOT
WRITE(UNIT=ILUOUT0,FMT=4) ZTHERMO, 100.*ZTHERMO/ZTOT
WRITE(UNIT=ILUOUT0,FMT=6) ZDYN, 100.*ZDYN/ZTOT
WRITE(UNIT=ILUOUT0,FMT=7) ZDIAG, 100.*ZDIAG/ZTOT
WRITE(UNIT=ILUOUT0,FMT=8) ZWRITE, 100.*ZWRITE/ZTOT
WRITE(UNIT=ILUOUT0,FMT=1) ZMISC, 100.*ZMISC/ZTOT
WRITE(UNIT=ILUOUT0,FMT=5) ZSURF, 100.*ZSURF/ZTOT
!
WRITE(UNIT=ILUOUT0,FMT=10) ZTOT   , 100.*ZALL/ZTOT
WRITE(ILUOUT0,*) ' ------------------------------------------------------------ '
!
!                  FORMATS
!                  -------
!
2  FORMAT(' |   READING OF DATA   |     ',F8.3,'      |     ',F8.3,'     |')
9  FORMAT(' | HOR. INTERPOLATIONS |     ',F8.3,'      |     ',F8.3,'     |')
3  FORMAT(' |      VER_PREP       |     ',F8.3,'      |     ',F8.3,'     |')
4  FORMAT(' |     VER_THERMO      |     ',F8.3,'      |     ',F8.3,'     |')
6  FORMAT(' |      VER_DYN        |     ',F8.3,'      |     ',F8.3,'     |')
7  FORMAT(' |    DIAGNOSTICS      |     ',F8.3,'      |     ',F8.3,'     |')
8  FORMAT(' |       WRITE         |     ',F8.3,'      |     ',F8.3,'     |')
1  FORMAT(' |   MISCELLANEOUS     |     ',F8.3,'      |     ',F8.3,'     |')
5  FORMAT(' |      SURFACE        |     ',F8.3,'      |     ',F8.3,'     |')
10 FORMAT(' |   PREP_REAL_CASE    |     ',F8.3,'      |     ',F8.3,'     |')
!
!-------------------------------------------------------------------------------
!
IF (LEN_TRIM(YDAD_NAME)>0) THEN
  WRITE(ILUOUT0,*) ' '
  WRITE(ILUOUT0,*) ' ------------------------------------------------------------'
  WRITE(ILUOUT0,*) '|  Nesting allowed                                           |'
  WRITE(ILUOUT0,*) '|  DAD_NAME="',YDAD_NAME,'"                   |'
  WRITE(ILUOUT0,*) ' ------------------------------------------------------------'
  WRITE(ILUOUT0,*) ' '
ELSE
  WRITE(ILUOUT0,*) ' '
  WRITE(ILUOUT0,*) ' ------------------------------------------------------------'
  WRITE(ILUOUT0,*) '|  Nesting not allowed with a larger-scale model.            |'
  WRITE(ILUOUT0,*) '|  The new file can only be used as model number 1           |'
  WRITE(ILUOUT0,*) ' ------------------------------------------------------------'
  WRITE(ILUOUT0,*) ' '
END IF
!
!-------------------------------------------------------------------------------
!
CALL IO_File_close(TINIFILE)
CALL IO_File_close(TPGDFILE)
!
CALL FINALIZE_MNH()
!
!-------------------------------------------------------------------------------
!
CONTAINS 

SUBROUTINE INIT_NMLVAR
CPRESOPT=CPRESOPT_n
LRES=LRES_n
XRES=XRES_n
NITR=NITR_n
LUSECHAQ=LUSECHAQ_n
LUSECHIC=LUSECHIC_n
LUSECHEM=LUSECHEM_n
END SUBROUTINE INIT_NMLVAR

SUBROUTINE UPDATE_MODD_FROM_NMLVAR
CPRESOPT_n=CPRESOPT
LRES_n=LRES
XRES_n=XRES
NITR_n=NITR
LUSECHAQ_n=LUSECHAQ
LUSECHIC_n=LUSECHIC
LUSECHEM_n=LUSECHEM
END SUBROUTINE UPDATE_MODD_FROM_NMLVAR

END PROGRAM PREP_REAL_CASE
