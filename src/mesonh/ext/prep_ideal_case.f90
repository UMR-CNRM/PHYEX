!MNH_LIC Copyright 1994-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      PROGRAM PREP_IDEAL_CASE
!     #######################
!
!!****  *PREP_IDEAL_CASE* - program to write an initial FM-file 
!!
!!    PURPOSE
!!    -------
!       The purpose of this program is to prepare an initial meso-NH file
!     (LFIFM and DESFM files) filled with some idealized fields.    
!
!      ---- The present version can provide two types of fields:
!
!      1) CIDEAL = 'CSTN' : 3D fields derived  from a vertical profile with
!         ---------------   n levels of constant moist Brunt Vaisala frequency
!             The vertical profile is read in EXPRE file.                 
!             These fields can be used for model runs 
!
!      2) CIDEAL = 'RSOU' : 3D fields derived from a radiosounding.
!          --------------- 
!             The radiosounding is read in EXPRE file. 
!             The following kind of data  is permitted :
!                  YKIND = 'STANDARD'  :   Zsol, Psol, Tsol, TDsol
!                                         (Pressure, dd, ff) , 
!                                         (Pressure, T, Td)
!                  YKIND = 'PUVTHVMR'  : zsol, Psol, Thvsol, Rsol
!                                        (Pressure, U, V) , 
!                                        (Pressure, THv, R)
!                  YKIND = 'PUVTHVHU'  :  zsol, Psol, Thvsol, Husol
!                                         (Pressure, U, V) , 
!                                         (Pressure, THv, Hu)
!                  YKIND = 'ZUVTHVHU'  :  zsol, Psol, Thvsol, Husol
!                                         (height, U, V) , 
!                                         (height, THv, Hu)
!                  YKIND = 'ZUVTHVMR'  :  zsol, Psol, Thvsol, Rsol
!                                         (height, U, V) , 
!                                         (height, THv, R)
!                  YKIND = 'PUVTHDMR'  : zsol, Psol, Thdsol, Rsol
!                                         (Pressure, U, V) , 
!                                         (Pressure, THd, R)
!                  YKIND = 'PUVTHDHU'  : zsol, Psol, Thdsol, Husol
!                                         (Pressure, U, V) , 
!                                         (Pressure, THd, Hu)
!                  YKIND = 'ZUVTHDMR'  :  zsol, Psol, Thdsol, Rsol
!                                         (height, U, V) , 
!                                         (height, THd, R)
!                  YKIND = 'ZUVTHLMR'  :  zsol, Psol, Thdsol, Rsol
!                                         (height, U, V) , 
!                                         (height, THl, Rt)
!
!             These fields can be used for model runs 
!
!      Cases (1) and (2) can be balanced
!      (geostrophic, hydrostatic  and anelastic balances) if desired.
!
!      ---- The orography can be flat (YZS='FLAT'), but also 
!      sine-shaped (YZS='SINE') or  bell-shaped (YZS='BELL')
!
!      ---- The U(z)  profile given in the RSOU and CSTN cases can
!      be multiplied (CUFUN="Y*Z") by a function of y (function FUNUY)  
!      The V(z) profile  given in the RSOU and CSTN cases can
!      be multiplied (CVFUN="X*Z") by a function of x (function FUNVX). 
!      If it is not the case, i.e. U(y,z)=U(z) then CUFUN="ZZZ" and 
!      CVFUN="ZZZ" for V(y,z)=V(z). Instead of these separable forms,
!      non-separables functions FUNUYZ (CUFUN="Y,Z")  and FUNVXZ (CVFUN="X,Z") 
!      can be used to specify the wind components.
!
!!**  METHOD
!!    ------
!!      The directives and data to perform the preparation of the initial FM
!!    file are stored in EXPRE file. This file is composed  of two parts : 
!!          - a namelists-format  part which is present in all cases
!!          - a free-format  part which contains data in cases 
!!       of discretised orography (CZS='DATA')
!!       of radiosounding (CIDEAL='RSOU') or Nv=cste  profile (CIDEAL='CSTN')
!!       of forced version (LFORCING=.TRUE.)
!!    
!!
!!      The following  PREP_IDEAL_CASE program  :
!!
!!             - initializes physical constants by calling INI_CST 
!!
!!             - sets default values for global variables which will be 
!!     written  in DESFM file and for variables in EXPRE file (namelists part)
!!     which will be written in LFIFM file.    
!!
!!             - reads the namelists part of EXPRE file which gives 
!!     informations about the preinitialization to perform,
!!
!!             - allocates memory for arrays, 
!!
!!             - initializes fields depending on the 
!!              directives  (CIDEAL in namelist NAM_CONF_PRE) :
!!  
!!                * grid variables : 
!!                  The gridpoints are regularly spaced by XDELTAX, XDELTAY.
!!               The grid is stretched along the z direction, the mesh varies 
!!               from XDZGRD near the ground to XDZTOP near the top and the 
!!               weigthing function is a TANH function characterized by its 
!!               center and width above and under this center
!!                  The orography is initialized following the kind of orography
!!               (YZS in namelist NAM_CONF_PRE) and the degrees of freedom :
!!                     sine-shape ---> ZHMAX, IEXPX,IEXPY
!!                     bell-shape ---> ZHMAX, ZAX,ZAY,IIZS,IJZS
!!                  The horizontal grid variables are initialized following
!!                the kind of geometry (LCARTESIAN in namelist NAM_CONF_PRE) 
!!                and the grid parameters XLAT0,XLON0,XBETA in both geometries
!!                and XRPK,XLONORI,XLATORI  in conformal projection.
!!                  In the  case of initialization from a radiosounding, the
!!                date and time is read in free-part of the EXPRE file. In other
!!                cases year, month and day are set to NUNDEF and time to 0.
!!
!!               * prognostic fields : 
!!
!!                     U,V,W, Theta and r. are first determined. They are
!!                multiplied by rhoj after the anelastic reference state 
!!                computation.
!!                     For the CSTN and RSOU cases, the determination of 
!!                Theta and rv is performed  respectively by SET_RSOU
!!                and by SET_CSTN which call the common routine SET_MASS. 
!!                These three routines have  the following actions :
!!          ---   The input vertical profile   is converted in 
!!                variables (U,V,thetav,r) and  interpolated
!!                on a mixed grid (with VERT_COORD) as in PREP_REAL_CASE 
!!          ---   A variation of the u-wind component( x-model axis component) 
!!                 is possible in y direction, a variation of the v-wind component 
!!                (y-model axis component) is possible in x direction.
!!          ---   Thetav could be computed with thermal wind balance
!!                (LGEOSBAL=.TRUE. with call of SET_GEOSBAL)                 
!!          ---   The mass fields (theta and r ) and the wind components are 
!!                then interpolated on the model grid with orography  as in
!!                PREP_REAL_CASE with the option LSHIFT                
!!          ---   An  anelastic correction is  applied in PRESSURE_IN_PREP in
!!                the case of non-vanishing orography.    
!!            
!!               * anelastic reference state variables :
!!
!!                   1D reference state : 
!!                     RSOU and CSTN cases : rhorefz and thvrefz are computed 
!!                         by   SET_REFZ (called by SET_MASS).
!!                         They are deduced from thetav and r on the model grid
!!                         without orography.
!!                   The 3D reference state is  computed by SET_REF   
!!            
!!               * The total mass of dry air is computed by TOTAL_DMASS              
!!
!!             - writes the DESFM file, 
!!
!!             - writes the LFIFM file . 
!!
!!    EXTERNAL
!!    --------
!!      DEFAULT_DESFM : to set default values for variables which can be 
!!                      contained in DESFM file
!!      DEFAULT_EXPRE : to  set default values for other global variables 
!!                      which can be contained in namelist-part of EXPRE file
!!      Module MODE_GRIDPROJ : contains conformal projection routines
!!           SM_GRIDPROJ   : to compute some grid variables, in
!!                           case of conformal projection.
!!      Module MODE_GRIDCART : contains cartesian geometry routines
!!           SM_GRIDCART   : to compute some grid variables, in
!!                           case of cartesian geometry.
!!      SET_RSOU      : to initialize mass fields from a radiosounding
!!      SET_CSTN      : to initialize mass fields from a vertical profile of 
!!                      n layers of Nv=cste 
!!      SET_REF       : to compute  rhoJ 
!!      RESSURE_IN_PREP : to apply an anelastic correction in the case of
!!                        non-vanishing orography 
!!      IO_File_open : to open a FM-file (DESFM + LFIFM)
!!      WRITE_DESFM   : to write the  DESFM file
!!      WRI_LFIFM     : to write the   LFIFM file  
!!      IO_File_close : to close a FM-file (DESFM + LFIFM)
!!
!!      MXM,MYM,MZM   : Shuman operators
!!      WGUESS        : to compute W with the continuity equation from 
!!                      the U,V values 
!!
!!
!!      
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains parameters
!!      Module MODD_DIM1      : contains dimensions 
!!      Module MODD_CONF       : contains  configuration variables for 
!!                                all models
!!      Module MODD_CST        : contains physical constants
!!      Module MODD_GRID       : contains grid variables  for all models
!!      Module MODD_GRID1     : contains grid variables
!!      Module MODD_TIME      : contains time variables for all models  
!!      Module MODD_TIME1     : contains time variables  
!!      Module MODD_REF        : contains reference state variables for
!!                               all models
!!      Module MODD_REF1      : contains reference state variables 
!!      Module MODD_LUNIT      : contains variables which concern names
!!                            and logical unit numbers of files  for all models
!!      Module MODD_FIELD1    : contains prognostics  variables
!!      Module MODD_GR_FIELD1 : contains the surface prognostic variables 
!!      Module MODD_LSFIELD1    : contains Larger Scale fields
!!      Module MODD_DYN1        : contains dynamic control variables for model 1
!!      Module MODD_LBC1        : contains lbc control variables for model 1
!!
!!
!!      Module MODN_CONF1    : contains  configuration variables for model 1
!!                               and the NAMELIST list
!!      Module MODN_LUNIT1    : contains variables which concern names
!!                               and logical unit numbers of files and 
!!                               the NAMELIST list
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (program PREP_IDEAL_CASE)
!!    
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                            05/05/94
!!      updated                V. Ducrocq   27/06/94   
!!      updated                P.M.         27/07/94
!!      updated                V. Ducrocq   23/08/94 
!!      updated                V. Ducrocq   01/09/94 
!!      namelist changes       J. Stein     26/10/94 
!!      namelist changes       J. Stein     04/11/94 
!!      remove the second step of the geostrophic balance 14/11/94 (J.Stein)
!!      add grid stretching in the z direction + Larger scale fields +
!!      cleaning                                           6/12/94 (J.Stein) 
!!      periodize the orography and the grid sizes in the periodic case
!!                                                        19/12/94 (J.Stein) 
!!      correct a bug in the Larger Scale Fields initialization
!!                                                        19/12/94 (J.Stein) 
!!      add the vertical grid stretching                  02/01/95 (J. Stein)
!!      Total mass of dry air computation                 02/01/95 (J.P.Lafore) 
!!      add the 1D switch                                 13/01/95 (J. Stein)
!!      enforce a regular vertical grid if desired        18/01/95 (J. Stein)
!!      add the tdtcur initialization                     26/01/95 (J. Stein)
!!      bug in the test of the type of RS localization    25/02/95 (J. Stein)
!!      remove R from the historical variables            16/03/95 (J. Stein)
!!      error on the grid stretching                      30/06/95 (J. Stein)
!!      add the soil fields                               01/09/95 (S.Belair)
!!      change the streching function  and the wind guess
!!        (J. Stein and V.Masson)                         21/09/95 
!!      reset to FALSE LUSERC,..,LUSERH                   12/12/95 (J. Stein)
!!      enforce the RS localization in 1D and 2D config.
!!      + add the 'TSZ0' option for the soil variables    28/01/96 (J. Stein)
!!      initialization of domain from center point        31/01/96 (V. Masson)
!!      add the constant file reading                     05/02/96 (J. Stein)
!!      enter vertical model levels values                20/10/95 (T.Montmerle)
!!      add LFORCING option                               19/02/96 (K. Suhre)
!!      modify structure of NAM_CONF_PRE                  20/02/96 (J.-P. Pinty)
!!      default of the domain center when use of pgd file 12/03/96 (V. Masson)
!!      change the surface initialization                 20/03/96 ( Stein,
!!                                                    Bougeault, Kastendeutsch )
!!      change the DEFAULT_DESFMN CALL                    17/04/96 ( Lafore )
!!      set the STORAGE_TYPE to 'TT' (a single instant)   30/04/96 (Stein, 
!!                                                    Jabouille)
!!      new wguess to spread  the divergence              15/05/96 (Stein)
!!      set LTHINSHELL to TRUE + return to the old wguess 29/08/96 (Stein)
!!      MY_NAME and DAD_NAME writing for nesting          30/07/96 (Lafore)
!!      MY_NAME and DAD_NAME reading in pgd file          26/09/96 (Masson)
!!       and reading of pgd grid in a new routine
!!      XXHAT and XYHAT are set to 0. at origine point    02/10/96 (Masson)
!!      add LTHINSHELL in namelist NAM_CONF_PRE           08/10/96 (Masson)
!!      restores use of TS and T2                         26/11/96 (Masson)
!!      value  XUNDEF for soil and vegetation fields on sea 27/11/96 (Masson)
!!      use of HUG and HU2 in both ISBA and TSZ0 cases    04/12/96 (Masson)
!!      add initialization of chemical variables          06/08/96 (K. Suhre)
!!      add MANUAL option for the terrain elevation       12/12/96 (J.-P. Pinty)
!!      set DATA instead of MANUAL for the terrain
!!      elevation option
!!      add new anelastic equations' systems              29/06/97 (Stein)
!!      split mode_lfifm_pgd                              29/07/97 (Masson)
!!      add directional z0 and subgrid scale orography    31/07/97 (Masson)
!!      separates surface treatment in PREP_IDEAL_SURF    15/03/99 (Masson)
!!      new PGD fields allocations                        15/03/99 (Masson)
!!      iterative call to pressure solver                 15/03/99 (Masson)
!!      removes TSZ0 case                                 04/01/00 (Masson)
!!      parallelization                                   18/06/00 (Pinty)
!!      adaptation for patch approach                     02/07/00 (Solmon/Masson)
!!      bug in W LB field on Y direction                  05/03/01 (Stein)
!!      add module MODD_NSV for NSV variable              01/02/01 (D. Gazen) 
!!      allow namelists in different orders               15/10/01 (I. Mallet)
!!      allow LUSERC and LUSERI in 1D configuration       05/06/02 (P. Jabouille)
!!      add  ZUVTHLMR case (move in set_rsou latter)      05/12/02 Jabouille/Masson
!!      move LHORELAX_SV (after INI_NSV)                  30/04/04 (Pinty)
!!      Correction Parallel bug IBEG & IDEND  evalution   13/11/08 J.Escobar
!!      add the option LSHIFT for interpolation of        26/10/10 (G.Tanguy)
!!      correction for XHAT & parallelizarion of ZSDATA   23/09/11 J.Escobar
!!      the vertical profile (as in PREP_REAL_CASE)
!!      add use MODI of SURFEX routines                   10/10/111 J.Escobar
!!
!!      For 2D modeling: 
!!      Initialization of ADVFRC profiles (SET_ADVFRC)    06/2010 (P.Peyrille)
!!      when LDUMMY(2)=T in PRE_IDEA1.nam 
!!      USE MODDB_ADVFRC_n for grid-nesting               02*2012 (M. Tomasini)
!!      LBOUSS in MODD_REF                                07/2013 (C.Lac)
!!      Correction for ZS in PGD file                     04/2014 (G. TANGUY)
!!      Bug : remove NC WRITE_HGRID                       05/2014 (S. Bielli via J.Escobar )
!!      BUG if ZFRC and ZFRC_ADV or ZFRC_REL are used together  11/2014 (G. Delautier)
!!      Bug : detected with cray compiler ,
!!                  missing '&' in continuation string  3/12/2014 J.Escobar
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  06/2016     (G.Delautier) phasage surfex 8
!!      P.Wautelet : 08/07/2016 : removed MNH_NCWRIT define
!!  01/2018      (G.Delautier) SURFEX 8.1
!  P. Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 07/02/2019: force TYPE to a known value for IO_File_add2list
!  P. Wautelet 14/02/2019: remove CLUOUT/CLUOUT0 and associated variables
!  P. Wautelet 28/03/2019: use MNHTIME for time measurement variables
!  P. Wautelet 28/03/2019: use TFILE instead of unit number for set_iluout_timing
!  P. Wautelet 19/04/2019: removed unused dummy arguments and variables
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  F. Auguste     02/2021: add IBM
!  P. Wautelet 09/03/2021: move some chemistry initializations to ini_nsv
!  Jean-Luc Redelsperger 03/2021: ocean LES case
!  P. Wautelet 06/07/2021: use FINALIZE_MNH
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS       ! Declarative modules
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_BUDGET,           ONLY: TBUCONF_ASSOCIATE
USE MODD_DIM_n
USE MODD_CONF
USE MODD_CST
USE MODD_GRID
USE MODD_GRID_n
USE MODD_IBM_LSF,     ONLY: CIBM_TYPE, LIBM_LSF, NIBM_SMOOTH, XIBM_SMOOTH
USE MODD_IBM_PARAM_n, ONLY: XIBM_LS
USE MODD_METRICS_n
USE MODD_LES, ONLY : LES_ASSOCIATE
USE MODD_PGDDIM
USE MODD_PGDGRID
USE MODD_TIME
USE MODD_TIME_n
USE MODD_REF
USE MODD_REF_n
USE MODD_LUNIT
USE MODD_FIELD_n
USE MODD_DYN_n
USE MODD_LBC_n
USE MODD_LSFIELD_n
USE MODD_PARAM_n
USE MODD_CH_MNHC_n,        ONLY:  LUSECHEM, LUSECHAQ, LUSECHIC, LCH_PH, LCH_INIT_FIELD
USE MODD_CH_AEROSOL,ONLY:  LORILAM, CORGANIC, LVARSIGI, LVARSIGJ, LINITPM, XINIRADIUSI, &
                           XINIRADIUSJ, XINISIGI, XINISIGJ, XN0IMIN, XN0JMIN, CRGUNIT
USE MODD_DUST,      ONLY:  LDUST, NMODE_DST, CRGUNITD, XINISIG, XINIRADIUS, XN0MIN 
USE MODD_SALT,      ONLY:  LSALT, NMODE_SLT, CRGUNITS, XINISIG_SLT, XINIRADIUS_SLT, XN0MIN_SLT
USE MODD_VAR_ll,    ONLY:  NPROC
USE MODD_LUNIT,     ONLY:  TLUOUT0, TOUTDATAFILE
USE MODD_LUNIT_n
USE MODD_IO,        ONLY: TFILE_DUMMY, TFILE_OUTPUTLISTING
USE MODD_CONF_n
USE MODD_NSV,       ONLY: NSV, NSV_ASSOCIATE
use modd_precision, only: LFIINT, MNHREAL_MPI, MNHTIME
!
USE MODN_BLANK_n
!
USE MODE_FINALIZE_MNH,     only: FINALIZE_MNH
USE MODE_THERMO
USE MODE_POS
USE MODE_GRIDCART         ! Executive modules
USE MODE_GRIDPROJ
USE MODE_GATHER_ll
USE MODE_IO,               only: IO_Config_set, IO_Init, IO_Pack_set
USE MODE_IO_FIELD_READ,    only: IO_Field_read
USE MODE_IO_FIELD_WRITE,   only: IO_Field_write, IO_Header_write
USE MODE_IO_FILE,          only: IO_File_close, IO_File_open
USE MODE_IO_MANAGE_STRUCT, only: IO_File_add2list
USE MODE_ll
USE MODE_MODELN_HANDLER
use mode_field,            only: Alloc_field_scalars, Ini_field_list, Ini_field_scalars
USE MODE_MSG
USE MODE_SET_GRID,         only: INTERP_HORGRID_TO_MASSPOINTS, STORE_GLOB_HORGRID
!
USE MODI_DEFAULT_DESFM_n    ! Interface modules
USE MODI_DEFAULT_EXPRE
USE MODI_IBM_INIT_LS
USE MODI_READ_HGRID
USE MODI_SHUMAN
USE MODI_SET_RSOU
USE MODI_SET_CSTN
USE MODI_SET_FRC
USE MODI_PRESSURE_IN_PREP
USE MODI_WRITE_DESFM_n
USE MODI_WRITE_LFIFM_n
USE MODI_METRICS
USE MODI_UPDATE_METRICS
USE MODI_SET_REF
USE MODI_SET_PERTURB
USE MODI_TOTAL_DMASS
USE MODI_CH_INIT_FIELD_n
USE MODI_INI_NSV
USE MODI_READ_PRE_IDEA_NAM_n
USE MODI_ZSMT_PIC
USE MODI_ZSMT_PGD
USE MODI_READ_VER_GRID
USE MODI_READ_ALL_NAMELISTS
USE MODI_PGD_GRID_SURF_ATM
USE MODI_SPLIT_GRID
USE MODI_PGD_SURF_ATM
USE MODI_ICE_ADJUST_BIS
USE MODI_WRITE_PGD_SURF_ATM_n
USE MODI_PREP_SURF_MNH
USE MODI_INIT_SALT
USE MODI_AER2LIMA
USE MODD_PARAM_LIMA
!
!JUAN
USE MODE_SPLITTINGZ_ll
USE MODD_SUB_MODEL_n
USE MODE_MNH_TIMING
USE MODN_CONFZ
!JUAN
!
USE MODI_VERSION
USE MODI_INIT_PGD_SURF_ATM
USE MODI_WRITE_SURF_ATM_N
USE MODD_MNH_SURFEX_n
! Modif ADVFRC
USE MODD_2D_FRC
USE MODD_ADVFRC_n     ! Modif for grid-nesting
USE MODI_SETADVFRC
USE MODD_RELFRC_n     ! Modif for grid-nesting
USE MODI_SET_RELFRC
!
USE MODE_INI_CST, ONLY: INI_CST
USE MODD_NEB_n, ONLY: NEBN
USE MODI_WRITE_HGRID
USE MODD_MPIF
USE MODD_VAR_ll
USE MODD_IO, ONLY: TFILEDATA,TFILE_SURFEX
!
USE MODE_MPPDB
!
USE MODD_GET_n
!
USE MODN_CONFIO, ONLY : NAM_CONFIO
!
IMPLICIT NONE
!
!*       0.1  Declarations of global variables not declared in the modules
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XJ ! Jacobian
REAL :: XLATCEN=XUNDEF, XLONCEN=XUNDEF ! latitude and longitude of the center of
                                     ! the domain for initialization. This 
                                     ! point is vertical vorticity point
                                     !          ------------------------
REAL :: XDELTAX=0.5E4, XDELTAY=0.5E4 ! horizontal mesh lengths  
                                     !  used to determine  XXHAT,XYHAT
!
INTEGER :: NLUPRE,NLUOUT           ! Logical unit numbers for EXPRE file
                                   ! and for output_listing file
INTEGER :: NRESP                   ! return code in FM routines
INTEGER :: NTYPE                   ! type of file (cpio or not)
INTEGER(KIND=LFIINT) :: NNPRAR     ! number of articles predicted in the LFIFM file
LOGICAL :: GFOUND                  ! Return code when searching namelist
!
INTEGER :: JLOOP,JILOOP,JJLOOP     ! Loop indexes
!
INTEGER :: NIB,NJB,NKB             ! Begining useful area  in x,y,z directions
INTEGER :: NIE,NJE                 ! Ending useful area  in x,y directions
INTEGER :: NIU,NJU,NKU             ! Upper bounds in x,y,z directions
CHARACTER(LEN=4)   :: CIDEAL ='CSTN'     ! kind of idealized fields
                                         ! 'CSTN' : Nv=cste case 
                                         ! 'RSOU' : radiosounding case
CHARACTER(LEN=4)   :: CZS    ='FLAT'     ! orography selector
                                         ! 'FLAT' : zero orography
                                         ! 'SINE' : sine-shaped orography 
                                         ! 'BELL' : bell-shaped orography 
REAL    :: XHMAX=XUNDEF            ! Maximum height for orography
REAL    :: NEXPX=3,NEXPY=1         ! Exponents for  orography in case of CZS='SINE'
REAL    :: XAX= 1.E4, XAY=1.E4     ! Widths for orography in case CZS='BELL'
                                   ! along x and y 
INTEGER :: NIZS = 5, NJZS = 5      ! Localization of the center in 
                                   ! case CZS ='BELL' 
!
!*       0.1.1 Declarations of local variables for N=cste and 
!              radiosounding cases :
!
INTEGER            :: NYEAR,NMONTH,NDAY ! year, month and day in EXPRE file
REAL               :: XTIME             ! time in EXPRE file
LOGICAL            :: LPERTURB =.FALSE. ! Logical to add a perturbation to 
                                        ! a basic state 
LOGICAL            :: LGEOSBAL =.FALSE. ! Logical to satisfy the geostrophic
                                        ! balance
                                        ! .TRUE. for geostrophic balance
                                        ! .FALSE. to ignore this balance
LOGICAL            :: LSHIFT   =.FALSE.  ! flag to perform vertical shift or not.        
CHARACTER(LEN=3)   :: CFUNU ='ZZZ'      ! CHARACTER STRING for variation of
                                        ! U in y direction
                                        ! 'ZZZ'  : U = U(Z)
                                        ! 'Y*Z'  : U = F(Y) * U(Z)
                                        ! 'Y,Z'  : U = G(Y,Z)
CHARACTER(LEN=3)   :: CFUNV ='ZZZ'      ! CHARACTER STRING for variation of
                                        ! V in x direction
                                        ! 'ZZZ'  : V = V(Z)
                                        ! 'Y*Z'  : V = F(X) * V(Z)
                                        ! 'Y,Z'  : V = G(X,Z)
CHARACTER(LEN=6)   :: CTYPELOC='IJGRID' ! Type of informations  used to give the
                                        ! localization of vertical profile
                                        ! 'IJGRID'  for (i,j) point  on index space
                                        ! 'XYHATM' for (x,y) coordinates on
                                        !  conformal or cartesian plane
                                        ! 'LATLON' for (latitude,longitude) on
                                        !   spherical earth  
REAL               :: XLATLOC= 45., XLONLOC=0.
                                        ! Latitude and longitude of the vertical
                                        ! profile localization  (used in case 
                                        ! CTYPELOC='LATLON') 
REAL               :: XXHATLOC=2.E4, XYHATLOC=2.E4 
                                        ! (x,y) of the vertical profile
                                        ! localization  (used in cases 
                                        ! CTYPELOC='LATLON' and 'XYHATM') 
INTEGER, DIMENSION(1) :: NILOC=4, NJLOC=4 
                                        ! (i,j) of the vertical profile
                                        ! localization 
!
!
REAL,DIMENSION(:,:,:),ALLOCATABLE   :: XCORIOZ ! Coriolis parameter (this
                                                 ! is exceptionnaly a 3D array
                                                 ! for computing needs)
!
!
!*       0.1.2 Declarations of local variables used when a PhysioGraphic Data
!              file is used :
!
INTEGER             :: JSV                      ! loop index on scalar var.
CHARACTER(LEN=28)   :: CPGD_FILE=' '            ! Physio-Graphic Data file name
LOGICAL  :: LREAD_ZS = .TRUE.,                & ! switch to use orography 
                                                ! coming from the PGD file
            LREAD_GROUND_PARAM = .TRUE.         ! switch to use soil parameters
                                                ! useful for the soil scheme
                                                ! coming from the PGD file

INTEGER           :: NSLEVE   =12         ! number of iteration for smooth orography
REAL              :: XSMOOTH_ZS = XUNDEF  ! optional uniform smooth orography for SLEVE coordinate
CHARACTER(LEN=28) :: YPGD_NAME, YPGD_DAD_NAME   ! general information
CHARACTER(LEN=2)  :: YPGD_TYPE
!
INTEGER           :: IINFO_ll                   ! return code of // routines
TYPE(LIST_ll), POINTER :: TZ_FIELDS_ll           ! list of metric coefficient fields
!
INTEGER :: IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU     ! dimensions of the
INTEGER :: IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2       ! West-east LB arrays
INTEGER :: IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV     ! dimensions of the
INTEGER :: IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2       ! North-south LB arrays
INTEGER :: IBEG,IEND,IXOR,IXDIM,IYOR,IYDIM,ILBX,ILBY
!
REAL, DIMENSION(:,:,:), ALLOCATABLE ::ZTHL,ZT,ZRT,ZFRAC_ICE,&
                                      ZEXN,ZLVOCPEXN,ZLSOCPEXN,ZCPH, &
                                      ZRSATW, ZRSATI
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZBUF
                                 ! variables for adjustement
REAL                :: ZDIST
!
!JUAN TIMING
REAL(kind=MNHTIME), DIMENSION(2) :: ZTIME1, ZTIME2, ZEND, ZTOT
CHARACTER                 :: YMI
INTEGER                   :: IMI
!JUAN TIMING
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZZS_ll
INTEGER                           :: IJ
!
REAL              :: ZZS_MAX, ZZS_MAX_ll
INTEGER           :: IJPHEXT
!
TYPE(TFILEDATA),POINTER :: TZEXPREFILE  => NULL()
!
!
!*       0.2  Namelist declarations
!
NAMELIST/NAM_CONF_PRE/ LTHINSHELL,LCARTESIAN,    &! Declarations in MODD_CONF
                       LPACK,                    &!
                       NVERB,CIDEAL,CZS,         &!+global variables initialized
                       LBOUSS,LOCEAN,LPERTURB,   &! at their declarations
                       LFORCING,CEQNSYS,         &! at their declarations
                       LSHIFT,L2D_ADV_FRC,L2D_REL_FRC, &
                       NHALO , JPHEXT
NAMELIST/NAM_GRID_PRE/ XLON0,XLAT0,            & ! Declarations in MODD_GRID
                       XBETA,XRPK,             & 
                       XLONORI,XLATORI
NAMELIST/NAM_GRIDH_PRE/ XLATCEN,XLONCEN,       & ! local variables  initialized
                 XDELTAX,XDELTAY,              & ! at their declarations
                 XHMAX,NEXPX,NEXPY,            &
                 XAX,XAY,NIZS,NJZS
NAMELIST/NAM_VPROF_PRE/LGEOSBAL, CFUNU,CFUNV,   &! global variables initialized
                     CTYPELOC,XLATLOC,XLONLOC,  &!  at their declarations
                     XXHATLOC,XYHATLOC,NILOC,NJLOC
NAMELIST/NAM_REAL_PGD/CPGD_FILE,                 & ! Physio-Graphic Data file
                                                   !  name
                      LREAD_ZS,                  & ! switch to use orography 
                                                   ! coming from the PGD file
                      LREAD_GROUND_PARAM
NAMELIST/NAM_SLEVE/NSLEVE, XSMOOTH_ZS
!
!*       0.3  Auxillary Namelist declarations
!
NAMELIST/NAM_AERO_PRE/ LORILAM, LINITPM, XINIRADIUSI, XINIRADIUSJ, &
                       XINISIGI, XINISIGJ, XN0IMIN, XN0JMIN, CRGUNIT, &
                       LDUST, LSALT, CRGUNITD, CRGUNITS,&
                       NMODE_DST, XINISIG, XINIRADIUS, XN0MIN,&
                       XINISIG_SLT, XINIRADIUS_SLT, XN0MIN_SLT, &
                       NMODE_SLT
!
NAMELIST/NAM_IBM_LSF/ LIBM_LSF, CIBM_TYPE, NIBM_SMOOTH, XIBM_SMOOTH
!
!-------------------------------------------------------------------------------
!
!*       0.    PROLOGUE
!              --------
CALL MPPDB_INIT()
!
CALL GOTO_MODEL(1)
!
CALL IO_Init()
NULLIFY(TZ_FIELDS_ll)
CALL VERSION
CPROGRAM='IDEAL '
!
!JUAN TIMING
  XT_START     = 0.0_MNHTIME
  XT_STORE     = 0.0_MNHTIME
!
  CALL SECOND_MNH2(ZEND)
!
!JUAN TIMING
!
!*       1.    INITIALIZE PHYSICAL CONSTANTS :         
!              ------------------------------
!
NVERB = 5
CALL INI_CST
!
!-------------------------------------------------------------------------------
!
!
!*  	 2.    SET DEFAULT VALUES  :  
!              --------------------
!
!
!*       2.1  For variables in DESFM file
!
CALL ALLOC_FIELD_SCALARS()
CALL TBUCONF_ASSOCIATE()
CALL LES_ASSOCIATE()
CALL DEFAULT_DESFM_n(1)
CALL NSV_ASSOCIATE()
!
CSURF = "NONE"
!
!
!*       2.2  For other global variables in EXPRE file
!
CALL DEFAULT_EXPRE
!-------------------------------------------------------------------------------
!
!*  	 3.    READ THE EXPRE FILE :  
!              --------------------
!
!*       3.1   initialize logical unit numbers (EXPRE and output-listing files)
!              and open these files :
! 
! 
CALL IO_File_add2list(TLUOUT0,'OUTPUT_LISTING1','OUTPUTLISTING','WRITE')
CALL IO_File_open(TLUOUT0)
NLUOUT = TLUOUT0%NLU
!Set output files for PRINT_MSG
TLUOUT              => TLUOUT0
TFILE_OUTPUTLISTING => TLUOUT0
!
CALL IO_File_add2list(TZEXPREFILE,'PRE_IDEA1.nam','NML','READ')
CALL IO_File_open(TZEXPREFILE)
NLUPRE=TZEXPREFILE%NLU
!
!*       3.2   read in NLUPRE the namelist informations
!
WRITE(NLUOUT,FMT=*) 'attempt to read ',TRIM(TZEXPREFILE%CNAME),' file'
CALL POSNAM( TZEXPREFILE, 'NAM_REAL_PGD', GFOUND )
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_REAL_PGD)
!
!
CALL POSNAM( TZEXPREFILE, 'NAM_CONF_PRE', GFOUND )
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_CONF_PRE)
!JUANZ
CALL POSNAM( TZEXPREFILE, 'NAM_CONFZ', GFOUND )
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_CONFZ)
!JUANZ
CALL POSNAM( TZEXPREFILE, 'NAM_CONFIO', GFOUND )
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_CONFIO)
CALL IO_Config_set()
CALL POSNAM( TZEXPREFILE, 'NAM_GRID_PRE', GFOUND )
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_GRID_PRE)
CALL POSNAM( TZEXPREFILE, 'NAM_GRIDH_PRE', GFOUND )
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_GRIDH_PRE)
CALL POSNAM( TZEXPREFILE, 'NAM_VPROF_PRE', GFOUND )
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_VPROF_PRE)
CALL POSNAM( TZEXPREFILE, 'NAM_BLANKN', GFOUND )
CALL INIT_NAM_BLANKn
IF (GFOUND) THEN
  READ(UNIT=NLUPRE,NML=NAM_BLANKn)
  CALL UPDATE_NAM_BLANKn
END IF
CALL READ_PRE_IDEA_NAM_n( TZEXPREFILE )
CALL POSNAM( TZEXPREFILE, 'NAM_AERO_PRE', GFOUND )
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_AERO_PRE)
CALL POSNAM( TZEXPREFILE, 'NAM_IBM_LSF', GFOUND )
IF (GFOUND) READ(UNIT=NLUPRE,NML=NAM_IBM_LSF )
!
CALL INI_FIELD_LIST()
!
CALL INI_FIELD_SCALARS()
! Sea salt
CALL INIT_SALT
!
IF( LEN_TRIM(CPGD_FILE) /= 0 ) THEN 
  ! open the PGD_FILE
  CALL IO_File_add2list(TPGDFILE,TRIM(CPGD_FILE),'PGD','READ',KLFINPRAR=NNPRAR,KLFITYPE=2,KLFIVERB=NVERB)
  CALL IO_File_open(TPGDFILE)

  ! read the grid in the PGD file
  CALL IO_Field_read(TPGDFILE,'IMAX',  NIMAX)
  CALL IO_Field_read(TPGDFILE,'JMAX',  NJMAX)
  CALL IO_Field_read(TPGDFILE,'JPHEXT',IJPHEXT)

  IF ( CPGD_FILE /= CINIFILEPGD) THEN
     WRITE(NLUOUT,FMT=*) ' WARNING : in PRE_IDEA1.nam, in NAM_LUNITn you&
          & have CINIFILEPGD= ',CINIFILEPGD
     WRITE(NLUOUT,FMT=*) ' whereas in NAM_REAL_PGD you have CPGD_FILE = '&
          ,CPGD_FILE
     WRITE(NLUOUT,FMT=*) ' '
     WRITE(NLUOUT,FMT=*) ' CINIFILEPGD HAS BEEN SET TO  ',CPGD_FILE        
     CINIFILEPGD=CPGD_FILE
  END IF
  IF ( IJPHEXT .NE. JPHEXT ) THEN
     WRITE(NLUOUT,FMT=*) ' PREP_IDEAL_CASE : JPHEXT in PRE_IDEA1.nam/NAM_CONF_PRE ( or default value )&
        & JPHEXT=',JPHEXT
     WRITE(NLUOUT,FMT=*) ' different from PGD files=', CINIFILEPGD,' value JPHEXT=',IJPHEXT
     WRITE(NLUOUT,FMT=*) '-> JOB ABORTED'
     CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_IDEAL_CASE','')
     !WRITE(NLUOUT,FMT=*) ' JPHEXT HAS BEEN SET TO ', IJPHEXT
     !IJPHEXT = JPHEXT
  END IF
END IF
!
NIMAX_ll=NIMAX   !! _ll variables are global variables
NJMAX_ll=NJMAX   !! but the old names are kept in PRE_IDEA1.nam file
!
!*       3.3   check some parameters:
!
L1D=.FALSE. ; L2D=.FALSE.
!
IF ((NIMAX == 1).OR.(NJMAX == 1)) THEN 
  L2D=.TRUE.
  NJMAX_ll=1
  NIMAX_ll=MAX(NIMAX,NJMAX)
  WRITE(NLUOUT,FMT=*) ' NJMAX HAS BEEN SET TO 1 SINCE 2D INITIAL FILE IS REQUIRED &
                   & (L2D=TRUE) )' 
END IF
!
IF ((NIMAX == 1).AND.(NJMAX == 1)) THEN 
  L1D=.TRUE.
  NIMAX_ll = 1
  NJMAX_ll = 1
  WRITE(NLUOUT,FMT=*) ' 1D INITIAL FILE IS REQUIRED (L1D=TRUE) ' 
END IF
!
IF(.NOT. L1D) THEN
  LHORELAX_UVWTH=.TRUE.
  LHORELAX_RV=.TRUE.
ENDIF
!
NRIMX= MIN(JPRIMMAX,NIMAX_ll/2)
!
IF (L2D) THEN
  NRIMY=0
ELSE
  NRIMY= MIN(JPRIMMAX,NJMAX_ll/2)
END IF
!
IF (L1D) THEN
  NRIMX=0
  NRIMY=0
END IF
!
IF (L1D .AND. ( LPERTURB .OR. LGEOSBAL .OR.                &
               (.NOT. LCARTESIAN ) .OR. (.NOT. LTHINSHELL) ))THEN 
  LGEOSBAL   = .FALSE.
  LPERTURB   = .FALSE.
  LCARTESIAN = .TRUE. 
  LTHINSHELL = .TRUE. 
  WRITE(NLUOUT,FMT=*) ' LGEOSBAL AND LPERTURB HAVE BEEN SET TO FALSE &
                      & AND LCARTESIAN AND LTHINSHELL TO TRUE        &
                      & SINCE 1D INITIAL FILE IS REQUIRED (L1D=TRUE)' 
END IF
!
IF (LGEOSBAL .AND. LSHIFT ) THEN
  LSHIFT=.FALSE.
  WRITE(NLUOUT,FMT=*) ' LSHIFT HAS BEEN SET TO FALSE SINCE &
                        & LGEOSBAL=.TRUE. IS REQUIRED '
END IF
!                      
!*       3.4   compute the number of moist variables :
!
IF (.NOT.LUSERV) THEN
  LUSERV = .TRUE.
  WRITE(NLUOUT,FMT=*) ' LUSERV HAS BEEN RESET TO TRUE, SINCE A MOIST VARIABLE &
                   & IS PRESENT IN EXPRE FILE (CIDEAL = RSOU OR CSTN)' 
END IF
!
IF((LUSERI .OR. LUSERC).AND. (CIDEAL /= 'RSOU')) THEN
  !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_IDEAL_CASE','use of hydrometeors is only allowed in RSOU case')
ENDIF
IF (LUSERI) THEN
  LUSERC =.TRUE.
  LUSERR =.TRUE.
  LUSERI =.TRUE.
  LUSERS =.TRUE.
  LUSERG =.TRUE.
  LUSERH =.FALSE.
  CCLOUD='ICE3'
ELSEIF(LUSERC) THEN
  LUSERR =.FALSE.
  LUSERI =.FALSE.
  LUSERS =.FALSE.
  LUSERG =.FALSE.
  LUSERH =.FALSE.
  CCLOUD='REVE'
ELSE
  LUSERC =.FALSE.
  LUSERR =.FALSE.
  LUSERI =.FALSE.
  LUSERS =.FALSE.
  LUSERG =.FALSE.
  LUSERH =.FALSE.
  LHORELAX_RC=.FALSE.
  LHORELAX_RR=.FALSE.
  LHORELAX_RI=.FALSE.
  LHORELAX_RS=.FALSE.
  LHORELAX_RG=.FALSE.
  LHORELAX_RH=.FALSE.
  CCLOUD='NONE'
!
END IF
!
NRR=0
IF (LUSERV) THEN
  NRR=NRR+1
  IDX_RVT = NRR
END IF
IF (LUSERC) THEN
  NRR=NRR+1
  IDX_RCT = NRR
END IF
IF (LUSERR) THEN
  NRR=NRR+1
  IDX_RRT = NRR
END IF
IF (LUSERI) THEN
  NRR=NRR+1
  IDX_RIT = NRR
END IF
IF (LUSERS) THEN
  NRR=NRR+1
  IDX_RST = NRR
END IF
IF (LUSERG) THEN
  NRR=NRR+1
  IDX_RGT = NRR
END IF
IF (LUSERH) THEN
  NRR=NRR+1
  IDX_RHT = NRR
END IF
!
! NRR=4 for RSOU case because RI and Rc always computed
IF (CIDEAL == 'RSOU' .AND. NRR < 4 ) NRR=4
!                   
!
!*       3.5   Chemistry
!
IF (LORILAM .OR. LCH_INIT_FIELD) THEN
  LUSECHEM = .TRUE.
  IF (LORILAM) THEN
    CORGANIC = "MPMPO"
    LVARSIGI = .TRUE.
    LVARSIGJ = .TRUE.
  END IF
END IF
! initialise NSV_* variables
CALL INI_NSV(1)
LHORELAX_SV(:)=.FALSE.
IF(.NOT. L1D) LHORELAX_SV(1:NSV)=.TRUE.
!
!-------------------------------------------------------------------------------
!
!*       4.    ALLOCATE MEMORY FOR ARRAYS :  
!   	       ----------------------------
!
!*       4.1  Vertical Spatial grid 
!
CALL READ_VER_GRID(TZEXPREFILE)
!
!*       4.2  Initialize parallel variables and compute array's dimensions
!
!
IF(LGEOSBAL) THEN
  CALL SET_SPLITTING_ll('XSPLITTING')  ! required for integration of thermal wind balance
ELSE
  CALL SET_SPLITTING_ll('BSPLITTING')
ENDIF
CALL SET_JP_ll(1,JPHEXT,JPVEXT,JPHEXT)
CALL SET_DAD0_ll()
CALL SET_DIM_ll(NIMAX_ll, NJMAX_ll, NKMAX)
CALL IO_Pack_set(L1D,L2D,LPACK)
CALL SET_LBX_ll(CLBCX(1), 1)
CALL SET_LBY_ll(CLBCY(1), 1)
CALL SET_XRATIO_ll(1, 1)
CALL SET_YRATIO_ll(1, 1)
CALL SET_XOR_ll(1, 1)
CALL SET_XEND_ll(NIMAX_ll+2*JPHEXT, 1)
CALL SET_YOR_ll(1, 1)
CALL SET_YEND_ll(NJMAX_ll+2*JPHEXT, 1)
CALL SET_DAD_ll(0, 1)
CALL INI_PARAZ_ll(IINFO_ll)
!
! sizes of arrays of the extended sub-domain
!
CALL GET_DIM_EXT_ll('B',NIU,NJU)
CALL GET_DIM_PHYS_ll('B',NIMAX,NJMAX)
CALL GET_INDICE_ll(NIB,NJB,NIE,NJE)
CALL GET_OR_ll('B',IXOR,IYOR)
NKB=1+JPVEXT
NKU=NKMAX+2*JPVEXT
!
!*       4.3  Global variables absent from the modules :
!
ALLOCATE(XJ(NIU,NJU,NKU))
SELECT CASE(CIDEAL)
  CASE('RSOU','CSTN')
    IF (LGEOSBAL) ALLOCATE(XCORIOZ(NIU,NJU,NKU))  ! exceptionally a 3D array  
  CASE DEFAULT                      ! undefined preinitialization
   !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_IDEAL_CASE','CIDEAL is not correctly defined')
END SELECT 
!
!*       4.4   Prognostic variables at M instant (module MODD_FIELD1):
!
ALLOCATE(XUT(NIU,NJU,NKU))
ALLOCATE(XVT(NIU,NJU,NKU))
ALLOCATE(XWT(NIU,NJU,NKU))
ALLOCATE(XTHT(NIU,NJU,NKU))
ALLOCATE(XPABST(NIU,NJU,NKU))
ALLOCATE(XRT(NIU,NJU,NKU,NRR))
ALLOCATE(XSVT(NIU,NJU,NKU,NSV))
!
!*       4.5   Grid variables (module MODD_GRID1 and MODD_METRICS1):
!
ALLOCATE(XMAP(NIU,NJU))
ALLOCATE(XLAT(NIU,NJU))
ALLOCATE(XLON(NIU,NJU))
ALLOCATE(XDXHAT(NIU),XDYHAT(NJU))
IF (LEN_TRIM(CPGD_FILE)==0) ALLOCATE(XZS(NIU,NJU))
IF (LEN_TRIM(CPGD_FILE)==0) ALLOCATE(ZZS_ll(NIMAX_ll))
IF (LEN_TRIM(CPGD_FILE)==0) ALLOCATE(XZSMT(NIU,NJU))
ALLOCATE(XZZ(NIU,NJU,NKU))
!
ALLOCATE(XDXX(NIU,NJU,NKU))
ALLOCATE(XDYY(NIU,NJU,NKU))
ALLOCATE(XDZX(NIU,NJU,NKU))
ALLOCATE(XDZY(NIU,NJU,NKU))
ALLOCATE(XDZZ(NIU,NJU,NKU))
!
!*       4.6   Reference state variables (modules MODD_REF and MODD_REF1):
!
ALLOCATE(XRHODREFZ(NKU),XTHVREFZ(NKU))
XTHVREFZ(:)=0.0
IF (LCOUPLES) THEN
  ! Arrays for reference state different in ocean and atmosphere
  ALLOCATE(XRHODREFZO(NKU),XTHVREFZO(NKU))
  XTHVREFZO(:)=0.0
END IF
IF(CEQNSYS == 'DUR') THEN
  ALLOCATE(XRVREF(NIU,NJU,NKU))
ELSE
  ALLOCATE(XRVREF(0,0,0))
END IF
ALLOCATE(XRHODREF(NIU,NJU,NKU),XTHVREF(NIU,NJU,NKU),XEXNREF(NIU,NJU,NKU))
ALLOCATE(XRHODJ(NIU,NJU,NKU))
!
!*       4.7   Larger Scale fields (modules MODD_LSFIELD1):
!
ALLOCATE(XLSUM(NIU,NJU,NKU))
ALLOCATE(XLSVM(NIU,NJU,NKU))
ALLOCATE(XLSWM(NIU,NJU,NKU))
ALLOCATE(XLSTHM(NIU,NJU,NKU))
IF ( NRR >= 1) THEN
  ALLOCATE(XLSRVM(NIU,NJU,NKU))
ELSE
  ALLOCATE(XLSRVM(0,0,0))
ENDIF
!
!  allocate lateral boundary field used for coupling
!
IF ( L1D) THEN                         ! 1D case
!
  NSIZELBX_ll=0
  NSIZELBXU_ll=0
  NSIZELBY_ll=0
  NSIZELBYV_ll=0
  NSIZELBXTKE_ll=0
  NSIZELBXR_ll=0
  NSIZELBXSV_ll=0
  NSIZELBYTKE_ll=0
  NSIZELBYR_ll=0
  NSIZELBYSV_ll=0
  ALLOCATE(XLBXUM(0,0,0))
  ALLOCATE(XLBYUM(0,0,0))
  ALLOCATE(XLBXVM(0,0,0))
  ALLOCATE(XLBYVM(0,0,0))
  ALLOCATE(XLBXWM(0,0,0))
  ALLOCATE(XLBYWM(0,0,0))
  ALLOCATE(XLBXTHM(0,0,0))
  ALLOCATE(XLBYTHM(0,0,0))
  ALLOCATE(XLBXTKEM(0,0,0))
  ALLOCATE(XLBYTKEM(0,0,0))
  ALLOCATE(XLBXRM(0,0,0,0))
  ALLOCATE(XLBYRM(0,0,0,0))
  ALLOCATE(XLBXSVM(0,0,0,0))
  ALLOCATE(XLBYSVM(0,0,0,0))
!
ELSEIF( L2D ) THEN             ! 2D case (not yet parallelized)
!                                          
  CALL GET_SIZEX_LB(NIMAX_ll,NJMAX_ll,NRIMX,               &
                    IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU, &
                    IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2)
  NSIZELBY_ll=0
  NSIZELBYV_ll=0
  NSIZELBYTKE_ll=0
  NSIZELBYR_ll=0
  NSIZELBYSV_ll=0
  ALLOCATE(XLBYUM(0,0,0))
  ALLOCATE(XLBYVM(0,0,0))
  ALLOCATE(XLBYWM(0,0,0))
  ALLOCATE(XLBYTHM(0,0,0))
  ALLOCATE(XLBYTKEM(0,0,0))
  ALLOCATE(XLBYRM(0,0,0,0))
  ALLOCATE(XLBYSVM(0,0,0,0))
  !
  IF ( LHORELAX_UVWTH ) THEN
!JUAN A REVOIR TODO_JPHEXT
! <<<<<<< prep_ideal_case.f90
    ! NSIZELBX_ll=2*NRIMX+2
    ! NSIZELBXU_ll=2*NRIMX+2
    ALLOCATE(XLBXUM(IISIZEXFU,NJU,NKU))
    ALLOCATE(XLBXVM(IISIZEXF,NJU,NKU))
    ALLOCATE(XLBXWM(IISIZEXF,NJU,NKU))
    ALLOCATE(XLBXTHM(IISIZEXF,NJU,NKU))
! =======
    NSIZELBX_ll=2*NRIMX+2*JPHEXT
    NSIZELBXU_ll=2*NRIMX+2*JPHEXT
    ! ALLOCATE(XLBXUM(2*NRIMX+2*JPHEXT,NJU,NKU))
    ! ALLOCATE(XLBXVM(2*NRIMX+2*JPHEXT,NJU,NKU))
    ! ALLOCATE(XLBXWM(2*NRIMX+2*JPHEXT,NJU,NKU))
    ! ALLOCATE(XLBXTHM(2*NRIMX+2*JPHEXT,NJU,NKU))
! >>>>>>> 1.3.2.4.2.3.2.14.2.8.2.11.2.2
  ELSE
    NSIZELBX_ll= 2*JPHEXT     ! 2
    NSIZELBXU_ll=2*(JPHEXT+1) ! 4 
    ALLOCATE(XLBXUM(NSIZELBXU_ll,NJU,NKU))
    ALLOCATE(XLBXVM(NSIZELBX_ll,NJU,NKU))
    ALLOCATE(XLBXWM(NSIZELBX_ll,NJU,NKU))
    ALLOCATE(XLBXTHM(NSIZELBX_ll,NJU,NKU))
  END IF  
  !
  IF ( NRR > 0 ) THEN
    IF (       LHORELAX_RV .OR. LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI    &
          .OR. LHORELAX_RS .OR. LHORELAX_RG .OR. LHORELAX_RH                     &
       ) THEN 
!JUAN A REVOIR TODO_JPHEXT
! <<<<<<< prep_ideal_case.f90
      ! NSIZELBXR_ll=2* NRIMX+2
      ALLOCATE(XLBXRM(IISIZEXF,NJU,NKU,NRR))
! =======
      NSIZELBXR_ll=2*NRIMX+2*JPHEXT
      ! ALLOCATE(XLBXRM(2*NRIMX+2*JPHEXT,NJU,NKU,NRR))
! >>>>>>> 1.3.2.4.2.3.2.14.2.8.2.11.2.2
    ELSE
      NSIZELBXR_ll=2*JPHEXT ! 2
      ALLOCATE(XLBXRM(NSIZELBXR_ll,NJU,NKU,NRR))
    ENDIF
  ELSE
    NSIZELBXR_ll=0
    ALLOCATE(XLBXRM(0,0,0,0))
  END IF
  !
  IF ( NSV > 0 ) THEN 
    IF ( ANY( LHORELAX_SV(:)) ) THEN
!JUAN A REVOIR TODO_JPHEXT
! <<<<<<< prep_ideal_case.f90
      ! NSIZELBXSV_ll=2* NRIMX+2
      ALLOCATE(XLBXSVM(IISIZEXF,NJU,NKU,NSV))
! =======
      NSIZELBXSV_ll=2*NRIMX+2*JPHEXT
      ! ALLOCATE(XLBXSVM(2*NRIMX+2*JPHEXT,NJU,NKU,NSV))
! >>>>>>> 1.3.2.4.2.3.2.14.2.8.2.11.2.2
    ELSE
      NSIZELBXSV_ll=2*JPHEXT ! 2
      ALLOCATE(XLBXSVM(NSIZELBXSV_ll,NJU,NKU,NSV))
    END IF
  ELSE
    NSIZELBXSV_ll=0
    ALLOCATE(XLBXSVM(0,0,0,0))
  END IF
!
ELSE                                   ! 3D case
!
  CALL GET_SIZEX_LB(NIMAX_ll,NJMAX_ll,NRIMX,               &
                    IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU, &
                    IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2)
  CALL GET_SIZEY_LB(NIMAX_ll,NJMAX_ll,NRIMY,               &
                    IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV, &
                    IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2)
!
  IF ( LHORELAX_UVWTH ) THEN
    NSIZELBX_ll=2*NRIMX+2*JPHEXT
    NSIZELBXU_ll=2*NRIMX+2*JPHEXT
    NSIZELBY_ll=2*NRIMY+2*JPHEXT
    NSIZELBYV_ll=2*NRIMY+2*JPHEXT
    ALLOCATE(XLBXUM(IISIZEXFU,IJSIZEXFU,NKU))
    ALLOCATE(XLBYUM(IISIZEYF,IJSIZEYF,NKU))
    ALLOCATE(XLBXVM(IISIZEXF,IJSIZEXF,NKU))
    ALLOCATE(XLBYVM(IISIZEYFV,IJSIZEYFV,NKU))
    ALLOCATE(XLBXWM(IISIZEXF,IJSIZEXF,NKU))
    ALLOCATE(XLBYWM(IISIZEYF,IJSIZEYF,NKU))
    ALLOCATE(XLBXTHM(IISIZEXF,IJSIZEXF,NKU))
    ALLOCATE(XLBYTHM(IISIZEYF,IJSIZEYF,NKU))
  ELSE
    NSIZELBX_ll=2*JPHEXT      ! 2
    NSIZELBXU_ll=2*(JPHEXT+1) ! 4
    NSIZELBY_ll=2*JPHEXT      ! 2
    NSIZELBYV_ll=2*(JPHEXT+1) ! 4
    ALLOCATE(XLBXUM(IISIZEX4,IJSIZEX4,NKU))
    ALLOCATE(XLBYUM(IISIZEY2,IJSIZEY2,NKU))
    ALLOCATE(XLBXVM(IISIZEX2,IJSIZEX2,NKU))
    ALLOCATE(XLBYVM(IISIZEY4,IJSIZEY4,NKU))
    ALLOCATE(XLBXWM(IISIZEX2,IJSIZEX2,NKU))
    ALLOCATE(XLBYWM(IISIZEY2,IJSIZEY2,NKU))
    ALLOCATE(XLBXTHM(IISIZEX2,IJSIZEX2,NKU))
    ALLOCATE(XLBYTHM(IISIZEY2,IJSIZEY2,NKU))
  END IF  
  !
  IF ( NRR > 0 ) THEN
    IF (       LHORELAX_RV .OR. LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI    &
          .OR. LHORELAX_RS .OR. LHORELAX_RG .OR. LHORELAX_RH                     &
       ) THEN 
      NSIZELBXR_ll=2*NRIMX+2*JPHEXT
      NSIZELBYR_ll=2*NRIMY+2*JPHEXT
      ALLOCATE(XLBXRM(IISIZEXF,IJSIZEXF,NKU,NRR))
      ALLOCATE(XLBYRM(IISIZEYF,IJSIZEYF,NKU,NRR))
    ELSE
      NSIZELBXR_ll=2*JPHEXT    ! 2
      NSIZELBYR_ll=2*JPHEXT    ! 2
      ALLOCATE(XLBXRM(IISIZEX2,IJSIZEX2,NKU,NRR))
      ALLOCATE(XLBYRM(IISIZEY2,IJSIZEY2,NKU,NRR))
    ENDIF
  ELSE
    NSIZELBXR_ll=0
    NSIZELBYR_ll=0
    ALLOCATE(XLBXRM(0,0,0,0))
    ALLOCATE(XLBYRM(0,0,0,0))
  END IF
  !
  IF ( NSV > 0 ) THEN 
    IF ( ANY( LHORELAX_SV(:)) ) THEN
      NSIZELBXSV_ll=2*NRIMX+2*JPHEXT
      NSIZELBYSV_ll=2*NRIMY+2*JPHEXT
      ALLOCATE(XLBXSVM(IISIZEXF,IJSIZEXF,NKU,NSV))
      ALLOCATE(XLBYSVM(IISIZEYF,IJSIZEYF,NKU,NSV))
    ELSE
      NSIZELBXSV_ll=2*JPHEXT    ! 2
      NSIZELBYSV_ll=2*JPHEXT    ! 2
      ALLOCATE(XLBXSVM(IISIZEX2,IJSIZEX2,NKU,NSV))
      ALLOCATE(XLBYSVM(IISIZEY2,IJSIZEY2,NKU,NSV))
    END IF
  ELSE
    NSIZELBXSV_ll=0
    NSIZELBYSV_ll=0
    ALLOCATE(XLBXSVM(0,0,0,0))
    ALLOCATE(XLBYSVM(0,0,0,0))
  END IF
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       5.     INITIALIZE ALL THE MODEL VARIABLES
!   	        ----------------------------------
!
!
!*       5.1    Grid variables and RS localization:
!
!*       5.1.1  Horizontal Spatial grid :
!
IF( LEN_TRIM(CPGD_FILE) /= 0 ) THEN 
!--------------------------------------------------------
! the MESONH horizontal grid will be read in the PGD_FILE 
!--------------------------------------------------------
  CALL READ_HGRID(1,TPGDFILE,YPGD_NAME,YPGD_DAD_NAME,YPGD_TYPE)
! control the cartesian option
  IF( LCARTESIAN ) THEN
     WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE : IN GENERAL, THE USE OF A PGD_FILE &
                & IMPLIES THAT YOU MUST TAKE INTO ACCOUNT THE EARTH SPHERICITY'
     WRITE(NLUOUT,FMT=*) 'NEVERTHELESS, LCARTESIAN HAS BEEN KEPT TO TRUE'
  END IF   
!
!* use of the externalized surface
!
  CSURF = "EXTE"
!
! determine whether the model is flat or no
!
  ZZS_MAX = ABS( MAXVAL(XZS(NIB:NIU-JPHEXT,NJB:NJU-JPHEXT)))
  CALL MPI_ALLREDUCE(ZZS_MAX, ZZS_MAX_ll, 1, MNHREAL_MPI, MPI_MAX,  &
                     NMNH_COMM_WORLD,IINFO_ll)
  IF( ABS(ZZS_MAX_ll)  < 1.E-10 ) THEN
    LFLAT=.TRUE.
  ELSE
    LFLAT=.FALSE.
  END IF
!

ELSE
!------------------------------------------------------------------------
! the MESONH horizontal grid is built from the PRE_IDEA1.nam informations
!------------------------------------------------------------------------
!
  ALLOCATE( XXHAT(NIU),  XYHAT(NJU)  )
  ALLOCATE( XXHATM(NIU), XYHATM(NJU) )
!
! define the grid localization at the earth surface by the central point
! coordinates
!
  IF (XLONCEN/=XUNDEF .OR. XLATCEN/=XUNDEF) THEN
    IF (XLONCEN/=XUNDEF .AND. XLATCEN/=XUNDEF) THEN 
!
! it should be noted that XLATCEN and XLONCEN refer to a vertical
! vorticity point and (XLATORI, XLONORI) refer to the mass point of
! conformal coordinates (0,0). This is to allow the centering of the model in
! a non-cyclic  configuration regarding to XLATCEN or XLONCEN.
!
      CALL SM_LATLON(XLATCEN,XLONCEN,                     &
                       -XDELTAX*(NIMAX_ll/2-0.5+JPHEXT),  &
                       -XDELTAY*(NJMAX_ll/2-0.5+JPHEXT),  &
                       XLATORI,XLONORI)
!
      WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE : XLATORI=' , XLATORI, &
                          ' XLONORI= ', XLONORI
    ELSE
   !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_IDEAL_CASE',&
                     'latitude and longitude of the center point must be initialized alltogether or not')
    END IF
  END IF
!
  IF (NPROC > 1) THEN
    CALL GET_DIM_EXT_ll('B',IXDIM,IYDIM)
    IBEG = IXOR-JPHEXT-1
    IEND = IBEG+IXDIM-1
    XXHAT(:) = (/ (REAL(JLOOP)*XDELTAX, JLOOP=IBEG,IEND) /)
    IBEG = IYOR-JPHEXT-1
    IEND = IBEG+IYDIM-1
    XYHAT(:) = (/ (REAL(JLOOP)*XDELTAY, JLOOP=IBEG,IEND) /)
!
  ELSE
    XXHAT(:) = (/ (REAL(JLOOP-NIB)*XDELTAX, JLOOP=1,NIU) /)
    XYHAT(:) = (/ (REAL(JLOOP-NJB)*XDELTAY, JLOOP=1,NJU) /)
  END IF

  ! Interpolations of positions to mass points
  CALL INTERP_HORGRID_TO_MASSPOINTS( XXHAT, XYHAT, XXHATM, XYHATM )

  ! Collect global domain boundaries
  CALL STORE_GLOB_HORGRID( XXHAT, XYHAT, XXHATM, XYHATM, XXHAT_ll, XYHAT_ll, XXHATM_ll, XYHATM_ll, XHAT_BOUND, XHATM_BOUND )

END IF
!
!*       5.1.2  Orography and Gal-Chen Sommerville transformation :
!
IF (    LEN_TRIM(CPGD_FILE) == 0  .OR. .NOT. LREAD_ZS) THEN
  SELECT CASE(CZS)     ! 'FLAT' or 'SINE' or 'BELL'
  CASE('FLAT')
    LFLAT = .TRUE.
    IF (XHMAX==XUNDEF) THEN
      XZS(:,:) = 0.
    ELSE
      XZS(:,:) = XHMAX
    END IF
  CASE('SINE')       ! sinus-shaped orography 
    IF (XHMAX==XUNDEF) XHMAX=300.
    LFLAT    =.FALSE.
    XZS(:,:) = XHMAX          &      ! three-dimensional case   
    *SPREAD((/((SIN((XPI/(NIMAX_ll+2*JPHEXT-1))*JLOOP)**2)**NEXPX,JLOOP=IXOR-1,IXOR+NIU-2)/),2,NJU) &
    *SPREAD((/((SIN((XPI/(NJMAX_ll+2*JPHEXT-1))*JLOOP)**2)**NEXPY,JLOOP=IYOR-1,IYOR+NJU-2)/),1,NIU)
    IF(L1D) THEN                     ! one-dimensional case
      XZS(:,:) = XHMAX 
    END IF        
  CASE('BELL')       ! bell-shaped orography 
    IF (XHMAX==XUNDEF) XHMAX=300.
    LFLAT = .FALSE.
    IF(.NOT.L2D) THEN                ! three-dimensional case
      XZS(:,:) = XHMAX  / ( 1.                                           &
        + ( (SPREAD(XXHAT(1:NIU),2,NJU) - REAL(NIZS) * XDELTAX) /XAX ) **2  &
        + ( (SPREAD(XYHAT(1:NJU),1,NIU) - REAL(NJZS) * XDELTAY) /XAY ) **2  ) **1.5
    ELSE                             ! two-dimensional case
      XZS(:,:) = XHMAX  / ( 1.                                          &
        + ( (SPREAD(XXHAT(1:NIU),2,NJU) - REAL(NIZS) * XDELTAX) /XAX ) **2 )
    ENDIF
    IF(L1D) THEN                     ! one-dimensional case
      XZS(:,:) = XHMAX 
    END IF        
  CASE('COSI')       ! (1+cosine)**4 shape
    IF (XHMAX==XUNDEF) XHMAX=800.
    LFLAT = .FALSE.
    IF(L2D) THEN                     ! two-dimensional case
      DO JILOOP = 1, NIU
        ZDIST = XXHAT(JILOOP)-REAL(NIZS)*XDELTAX
        IF( ABS(ZDIST)<(4.0*XAX) ) THEN
          XZS(JILOOP,:) = (XHMAX/16.0)*( 1.0 + COS((XPI*ZDIST)/(4.0*XAX)) )**4
        ELSE
          XZS(JILOOP,:) = 0.0
        ENDIF
      END DO
    ENDIF
  CASE('SCHA')       ! exp(-(x/a)**2)*cosine(pi*x/lambda)**2 shape
    IF (XHMAX==XUNDEF) XHMAX=800.
    LFLAT = .FALSE.
    IF(L2D) THEN                     ! two-dimensional case
      DO JILOOP = 1, NIU
        ZDIST = XXHAT(JILOOP)-REAL(NIZS)*XDELTAX
        IF( ABS(ZDIST)<(4.0*XAX) ) THEN
          XZS(JILOOP,:) = XHMAX*EXP(-(ZDIST/XAY)**2)*COS((XPI*ZDIST)/XAX)**2
        ELSE
          XZS(JILOOP,:) = 0.0
        ENDIF
      END DO
    ENDIF
  CASE('AGNE')       ! h*a**2/(x**2+a**2) shape
    LFLAT = .FALSE.
    IF(L2D) THEN                     ! two-dimensional case
      DO JILOOP = 1, NIU
        ZDIST = XXHAT(JILOOP)-REAL(NIZS)*XDELTAX
          XZS(JILOOP,:) = XHMAX*(XAX**2)/(XAX**2+ZDIST**2)
      END DO
		ELSE		! three dimensionnal case - infinite profile in y direction
			DO JILOOP = 1, NIU
        ZDIST = XXHAT(JILOOP)-REAL(NIZS)*XDELTAX
          XZS(JILOOP,:) = XHMAX*(XAX**2)/(XAX**2+ZDIST**2)
      END DO
    ENDIF

  CASE('DATA')       ! discretized orography
    LFLAT    =.FALSE.
    WRITE(NLUOUT,FMT=*) 'CZS="DATA",   ATTEMPT TO READ ARRAY     &
                    &XZS(NIB:NIU-JPHEXT:1,NJU-JPHEXT:NJB:-1) &
                    &starting from the first index'
    CALL POSKEY(NLUPRE,NLUOUT,'ZSDATA')
    DO JJLOOP = NJMAX_ll+2*JPHEXT-1,JPHEXT+1,-1    ! input like a map prior the sounding
      READ(NLUPRE,FMT=*) ZZS_ll
      IF ( ( JJLOOP <= ( NJU-JPHEXT + IYOR-1 ) ) .AND. ( JJLOOP >= ( NJB + IYOR-1 ) ) ) THEN
         IJ    = JJLOOP - ( IYOR-1 )
         XZS(NIB:NIU-JPHEXT,IJ) = ZZS_ll(IXOR:IXOR + NIU-JPHEXT - NIB )
      END IF
    END DO
!
  CASE DEFAULT   ! undefined  shape of orography
   !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_IDEAL_CASE','erroneous ground type')
  END SELECT
!
  CALL ADD2DFIELD_ll( TZ_FIELDS_ll, XZS, 'PREP_IDEAL_CASE::XZS' )
  CALL UPDATE_HALO_ll(TZ_FIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZ_FIELDS_ll)
!
END IF
!
!IF( ( LEN_TRIM(CPGD_FILE) /= 0 ) .AND. .NOT.LFLAT .AND. &
! ((CLBCX(1) /= "OPEN" ) .OR. &
! (CLBCX(2) /= "OPEN" ) .OR. (CLBCY(1) /= "OPEN" ) .OR. &
! (CLBCY(2) /= "OPEN" )) )  THEN 
!   !callabortstop
!  CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_IDEAL_CASE','with a PGD file, you cannot be in a cyclic LBC')
!END IF
!
IF (LWEST_ll())  THEN
  DO JILOOP = 1,JPHEXT
    XZS(JILOOP,:) = XZS(NIB,:)
  END DO
END IF
IF (LEAST_ll()) THEN
  DO JILOOP = NIU-JPHEXT+1,NIU
    XZS(JILOOP,:)=XZS(NIU-JPHEXT,:)
  END DO
END IF
IF (LSOUTH_ll()) THEN
  DO JJLOOP = 1,JPHEXT
    XZS(:,JJLOOP)=XZS(:,NJB)
  END DO
END IF
IF (LNORTH_ll()) THEN
  DO JJLOOP =NJU-JPHEXT+1,NJU
    XZS(:,JJLOOP)=XZS(:,NJU-JPHEXT)
  END DO
END IF
!
IF ( LEN_TRIM(CPGD_FILE) == 0  .OR. .NOT. LREAD_ZS) THEN
  IF (LSLEVE) THEN
    CALL ZSMT_PIC(NSLEVE,XSMOOTH_ZS)
  ELSE
    XZSMT(:,:) = 0.
  END IF
END IF
!
IF (LCARTESIAN) THEN
  CALL SM_GRIDCART(XXHAT,XYHAT,XZHAT,XZS,LSLEVE,XLEN1,XLEN2,XZSMT,XDXHAT,XDYHAT,XZZ,XJ)
  XMAP=1.
ELSE
  CALL SM_GRIDPROJ( XXHAT, XYHAT, XZHAT, XXHATM, XYHATM, XZS,      &
                    LSLEVE, XLEN1, XLEN2, XZSMT, XLATORI, XLONORI, &
                    XMAP, XLAT, XLON, XDXHAT, XDYHAT, XZZ, XJ      )
END IF
!*       5.4.1  metrics coefficients and update halos:
!
CALL METRICS(XMAP,XDXHAT,XDYHAT,XZZ,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
CALL UPDATE_METRICS(CLBCX,CLBCY,XDXX,XDYY,XDZX,XDZY,XDZZ)                         
!
!*       5.1.3  Compute the localization in index space of the vertical profile
!               in CSTN and RSOU cases  :
!
IF (CTYPELOC =='LATLON' ) THEN  
  IF (.NOT.LCARTESIAN) THEN                            ! compute (x,y) if 
    CALL SM_XYHAT(XLATORI,XLONORI,                 &   ! the localization 
                  XLATLOC,XLONLOC,XXHATLOC,XYHATLOC)   ! is given in latitude 
  ELSE                                                 ! and longitude
    WRITE(NLUOUT,FMT=*) 'CTYPELOC CANNOT BE LATLON IN CARTESIAN GEOMETRY'
    WRITE(NLUOUT,FMT=*) '-> JOB ABORTED'
   !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_IDEAL_CASE','CTYPELOC cannot be LATLON in cartesian geometry')
  END IF 
END IF  
!
IF (CTYPELOC /= 'IJGRID') THEN                                               
  NILOC = MINLOC(ABS(XXHATLOC-XXHAT_ll(:)))
  NJLOC = MINLOC(ABS(XYHATLOC-XYHAT_ll(:)))
END IF
!
IF ( L1D .AND. ( NILOC(1) /= 1 .OR. NJLOC(1) /= 1 ) ) THEN
  NILOC = 1
  NJLOC = 1
  WRITE(NLUOUT,FMT=*) 'FOR 1D CONFIGURATION, THE RS INFORMATIONS ARE TAKEN AT &
                      & I=1 AND J=1 (CENTRAL VERTICAL WITHOUT HALO)'
END IF
!
IF ( L2D .AND. ( NJLOC(1) /= 1 ) ) THEN
  NJLOC = 1
  WRITE(NLUOUT,FMT=*) 'FOR 2D CONFIGURATION, THE RS INFORMATIONS ARE TAKEN AT &
                      & J=1 (CENTRAL PLANE WITHOUT HALO)'
END IF
!
!*       5.2    Prognostic variables (not multiplied by  rhoJ) : u,v,w,theta,r
!               and 1D anelastic reference state
!
!
!*       5.2.1  Use a Radiosounding : CIDEAL='RSOU''
!
IF (CIDEAL == 'RSOU') THEN
  WRITE(NLUOUT,FMT=*) 'CIDEAL="RSOU", attempt to read DATE'
  CALL POSKEY(NLUPRE,NLUOUT,'RSOU')
  READ(NLUPRE,FMT=*)  NYEAR,NMONTH,NDAY,XTIME
  TDTCUR = DATE_TIME(NYEAR,NMONTH,NDAY,XTIME)
  TDTEXP = TDTCUR
  TDTSEG = TDTCUR
  TDTMOD = TDTCUR
  WRITE(NLUOUT,FMT=*) 'CIDEAL="RSOU", ATTEMPT TO PROCESS THE SOUNDING DATA'
  IF (LGEOSBAL) THEN
    CALL SET_RSOU(TFILE_DUMMY,TZEXPREFILE,CFUNU,CFUNV,NILOC(1),NJLOC(1),LBOUSS, &
                  XJ,LSHIFT,XCORIOZ)
  ELSE
    CALL SET_RSOU(TFILE_DUMMY,TZEXPREFILE,CFUNU,CFUNV,NILOC(1),NJLOC(1),LBOUSS, &
                  XJ,LSHIFT)
  END IF
!
!*       5.2.2  N=cste  and U(z) : CIDEAL='CSTN'
!
ELSE IF (CIDEAL == 'CSTN') THEN
  WRITE(NLUOUT,FMT=*) 'CIDEAL="CSTN", attempt to read DATE'
  CALL POSKEY(NLUPRE,NLUOUT,'CSTN')
  READ(NLUPRE,FMT=*)  NYEAR,NMONTH,NDAY,XTIME
  TDTCUR = DATE_TIME(NYEAR,NMONTH,NDAY,XTIME)
  TDTEXP = TDTCUR
  TDTSEG = TDTCUR
  TDTMOD = TDTCUR
  WRITE(NLUOUT,FMT=*) 'CIDEAL="CSTN", ATTEMPT TO PROCESS THE SOUNDING DATA'
  IF (LGEOSBAL) THEN
    CALL SET_CSTN(TFILE_DUMMY,TZEXPREFILE,CFUNU,CFUNV,NILOC(1),NJLOC(1),LBOUSS, &
                  XJ,LSHIFT,XCORIOZ)
  ELSE
    CALL SET_CSTN(TFILE_DUMMY,TZEXPREFILE,CFUNU,CFUNV,NILOC(1),NJLOC(1),LBOUSS, &
                  XJ,LSHIFT)
  END IF
!
END IF 
!
!*       5.3    Forcing variables
!
IF (LFORCING) THEN
  WRITE(NLUOUT,FMT=*) 'FORCING IS ENABLED, ATTEMPT TO SET FORCING FIELDS'
  CALL POSKEY(NLUPRE,NLUOUT,'ZFRC ','PFRC')
  CALL SET_FRC(TZEXPREFILE)
END IF
!
!! ---------------------------------------------------------------------
! Modif PP ADV FRC
! 5.4.2 initialize profiles for adv forcings
IF (L2D_ADV_FRC) THEN
    WRITE(NLUOUT,FMT=*) 'L2D_ADV_FRC IS SET TO  TRUE'
    WRITE(NLUOUT,FMT=*) 'ADVECTING FORCING USED IS USER MADE, NOT STANDARD ONE ' 
    WRITE(NLUOUT,FMT=*) 'IT IS FOR 2D IDEALIZED WAM STUDY ONLY ' 
   CALL POSKEY(NLUPRE,NLUOUT,'ZFRC_ADV')
   CALL SET_ADVFRC(TZEXPREFILE)
ENDIF
IF (L2D_REL_FRC) THEN
    WRITE(NLUOUT,FMT=*) 'L2D_REL_FRC IS SET TO  TRUE'
    WRITE(NLUOUT,FMT=*) 'RELAXATION FORCING USED IS USER MADE, NOT STANDARD ONE ' 
    WRITE(NLUOUT,FMT=*) 'IT IS FOR 2D IDEALIZED WAM STUDY ONLY ' 
   CALL POSKEY(NLUPRE,NLUOUT,'ZFRC_REL')
   CALL SET_RELFRC(TZEXPREFILE)
ENDIF
!*       5.4    3D Reference state variables :
!
!
!*       5.4.1  metrics coefficients and update halos:
!
CALL METRICS(XMAP,XDXHAT,XDYHAT,XZZ,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
CALL UPDATE_METRICS(CLBCX,CLBCY,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
!*       5.4.2  3D reference state :
!
CALL SET_REF( 0, TFILE_DUMMY,                            &
              XZZ, XZHATM, XJ, XDXX, XDYY, CLBCX, CLBCY, &
              XREFMASS, XMASS_O_PHI0, XLINMASS,          &
              XRHODREF, XTHVREF, XRVREF, XEXNREF, XRHODJ )
!
!
!*       5.5.1  Absolute pressure :
!
!
!*       5.5.2  Total mass of dry air Md computation :
!
CALL TOTAL_DMASS(XJ,XRHODREF,XDRYMASST)
!
!
!*       5.6    Complete prognostic variables (multipliy by  rhoJ) at time t :
!
! U grid   : gridpoint 2
IF (LWEST_ll())  XUT(1,:,:)    = 2.*XUT(2,:,:) - XUT(3,:,:)
! V grid   : gridpoint 3
IF (LSOUTH_ll())  XVT(:,1,:)    = 2.*XVT(:,2,:) - XVT(:,3,:)
! SV : gridpoint 1
XSVT(:,:,:,:) = 0.
!
!
!*       5.7   Larger scale fields initialization :
!
XLSUM(:,:,:) = XUT(:,:,:)        ! these fields do not satisfy the 
XLSVM(:,:,:) = XVT(:,:,:)        ! lower boundary condition but are 
XLSWM(:,:,:) = XWT(:,:,:)        ! in equilibrium
XLSTHM(:,:,:)= XTHT(:,:,:)
XLSRVM(:,:,:)= XRT(:,:,:,1)
!
! enforce the vertical homogeneity under the ground and above the top of
! the model for the LS fields
!
XLSUM(:,:,NKB-1)=XLSUM(:,:,NKB)
XLSUM(:,:,NKU)=XLSUM(:,:,NKU-1)
XLSVM(:,:,NKB-1)=XLSVM(:,:,NKB)
XLSVM(:,:,NKU)=XLSVM(:,:,NKU-1)
XLSWM(:,:,NKB-1)=XLSWM(:,:,NKB)
XLSWM(:,:,NKU)=XLSWM(:,:,NKU-1)
XLSTHM(:,:,NKB-1)=XLSTHM(:,:,NKB)
XLSTHM(:,:,NKU)=XLSTHM(:,:,NKU-1)
IF ( NRR > 0 ) THEN
  XLSRVM(:,:,NKB-1)=XLSRVM(:,:,NKB)
  XLSRVM(:,:,NKU)=XLSRVM(:,:,NKU-1)
END IF
!
ILBX=SIZE(XLBXUM,1)
ILBY=SIZE(XLBYUM,2)
IF(LWEST_ll() .AND. .NOT. L1D) THEN
  XLBXUM(1:NRIMX+JPHEXT,        :,:)     = XUT(2:NRIMX+JPHEXT+1,        :,:)
  XLBXVM(1:NRIMX+JPHEXT,        :,:)     = XVT(1:NRIMX+JPHEXT,        :,:)
  XLBXWM(1:NRIMX+JPHEXT,        :,:)     = XWT(1:NRIMX+JPHEXT,        :,:)
  XLBXTHM(1:NRIMX+JPHEXT,        :,:)   = XTHT(1:NRIMX+JPHEXT,        :,:)
  XLBXRM(1:NRIMX+JPHEXT,        :,:,:)   = XRT(1:NRIMX+JPHEXT,        :,:,:)
ENDIF
IF(LEAST_ll() .AND. .NOT. L1D) THEN
  XLBXUM(ILBX-NRIMX-JPHEXT+1:ILBX,:,:)     = XUT(NIU-NRIMX-JPHEXT+1:NIU,    :,:)
  XLBXVM(ILBX-NRIMX-JPHEXT+1:ILBX,:,:)     = XVT(NIU-NRIMX-JPHEXT+1:NIU,    :,:)
  XLBXWM(ILBX-NRIMX-JPHEXT+1:ILBX,:,:)     = XWT(NIU-NRIMX-JPHEXT+1:NIU,    :,:)
  XLBXTHM(ILBX-NRIMX-JPHEXT+1:ILBX,:,:)   = XTHT(NIU-NRIMX-JPHEXT+1:NIU,    :,:)
  XLBXRM(ILBX-NRIMX-JPHEXT+1:ILBX,:,:,:)   = XRT(NIU-NRIMX-JPHEXT+1:NIU,    :,:,:)
ENDIF
IF(LSOUTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) THEN
  XLBYUM(:,1:NRIMY+JPHEXT,        :)     = XUT(:,1:NRIMY+JPHEXT,      :)
  XLBYVM(:,1:NRIMY+JPHEXT,        :)     = XVT(:,2:NRIMY+JPHEXT+1,      :)
  XLBYWM(:,1:NRIMY+JPHEXT,        :)     = XWT(:,1:NRIMY+JPHEXT,  :)
  XLBYTHM(:,1:NRIMY+JPHEXT,        :)    = XTHT(:,1:NRIMY+JPHEXT,      :)
  XLBYRM(:,1:NRIMY+JPHEXT,        :,:)   = XRT(:,1:NRIMY+JPHEXT,      :,:)
ENDIF
IF(LNORTH_ll().AND. .NOT. L1D .AND. .NOT. L2D) THEN
  XLBYUM(:,ILBY-NRIMY-JPHEXT+1:ILBY,:)     = XUT(:,NJU-NRIMY-JPHEXT+1:NJU,  :)
  XLBYVM(:,ILBY-NRIMY-JPHEXT+1:ILBY,:)     = XVT(:,NJU-NRIMY-JPHEXT+1:NJU,  :)
  XLBYWM(:,ILBY-NRIMY-JPHEXT+1:ILBY,:)     = XWT(:,NJU-NRIMY-JPHEXT+1:NJU,  :)
  XLBYTHM(:,ILBY-NRIMY-JPHEXT+1:ILBY,:)    = XTHT(:,NJU-NRIMY-JPHEXT+1:NJU,  :)
  XLBYRM(:,ILBY-NRIMY-JPHEXT+1:ILBY,:,:)   = XRT(:,NJU-NRIMY-JPHEXT+1:NJU,  :,:)
ENDIF
DO JSV = 1, NSV
  IF(LWEST_ll() .AND. .NOT. L1D) &
  XLBXSVM(1:NRIMX+JPHEXT,        :,:,JSV)   = XSVT(1:NRIMX+JPHEXT,        :,:,JSV)
  IF(LEAST_ll() .AND. .NOT. L1D) &
  XLBXSVM(ILBX-NRIMX-JPHEXT+1:ILBX,:,:,JSV)   = XSVT(NIU-NRIMX-JPHEXT+1:NIU,    :,:,JSV)
  IF(LSOUTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) &
  XLBYSVM(:,1:NRIMY+JPHEXT,        :,JSV)   = XSVT(:,1:NRIMY+JPHEXT,      :,JSV)
  IF(LNORTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) &
  XLBYSVM(:,ILBY-NRIMY-JPHEXT+1:ILBY,:,JSV)   = XSVT(:,NJU-NRIMY-JPHEXT+1:NJU,  :,JSV)
END DO
!
!
!*       5.8   Add a perturbation to a basic state :
!
IF(LPERTURB) CALL SET_PERTURB(TZEXPREFILE)
!
!
!*       5.9   Anelastic correction and pressure:
!
IF (.NOT.LOCEAN) THEN
  CALL ICE_ADJUST_BIS(XPABST,XTHT,XRT)
  IF ( .NOT. L1D ) CALL PRESSURE_IN_PREP(XDXX,XDYY,XDZX,XDZY,XDZZ)
  CALL ICE_ADJUST_BIS(XPABST,XTHT,XRT)
END IF
!
!
!*       5.10  Compute THETA, vapor and cloud mixing ratio
!
IF (CIDEAL == 'RSOU') THEN
    ALLOCATE(ZEXN(NIU,NJU,NKU))         
  ALLOCATE(ZT(NIU,NJU,NKU))  
  ALLOCATE(ZTHL(NIU,NJU,NKU))              
  ALLOCATE(ZRT(NIU,NJU,NKU))              
  ALLOCATE(ZCPH(NIU,NJU,NKU))        
  ALLOCATE(ZLVOCPEXN(NIU,NJU,NKU))        
  ALLOCATE(ZLSOCPEXN(NIU,NJU,NKU))  
  ALLOCATE(ZFRAC_ICE(NIU,NJU,NKU))
  ALLOCATE(ZRSATW(NIU,NJU,NKU))
  ALLOCATE(ZRSATI(NIU,NJU,NKU))             
  ALLOCATE(ZBUF(NIU,NJU,NKU,16))
  ZRT=XRT(:,:,:,1)+XRT(:,:,:,2)+XRT(:,:,:,4)
IF (LOCEAN) THEN
  ZEXN(:,:,:)= 1.
  ZT=XTHT
  ZTHL=XTHT
  ZCPH=XCPD+ XCPV * XRT(:,:,:,1)
  ZLVOCPEXN = XLVTT
  ZLSOCPEXN = XLSTT
ELSE
  ZEXN=(XPABST/XP00) ** (XRD/XCPD)
  ZT=XTHT*(XPABST/XP00)**(XRD/XCPD)
  ZCPH=XCPD+ XCPV * XRT(:,:,:,1)+ XCL *XRT(:,:,:,2)  + XCI * XRT(:,:,:,4)
  ZLVOCPEXN = (XLVTT + (XCPV-XCL) * (ZT-XTT))/(ZCPH*ZEXN)
  ZLSOCPEXN = (XLSTT + (XCPV-XCI) * (ZT-XTT))/(ZCPH*ZEXN)
  ZTHL=XTHT-ZLVOCPEXN*XRT(:,:,:,2)-ZLSOCPEXN*XRT(:,:,:,4)
  CALL TH_R_FROM_THL_RT(CST, NEBN, SIZE(ZFRAC_ICE), 'T',ZFRAC_ICE,XPABST,ZTHL,ZRT,XTHT,XRT(:,:,:,1), &
                        XRT(:,:,:,2),XRT(:,:,:,4),ZRSATW, ZRSATI,OOCEAN=.FALSE.,&
                        PBUF=ZBUF)
END IF
  DEALLOCATE(ZEXN)         
  DEALLOCATE(ZT)       
  DEALLOCATE(ZCPH)        
  DEALLOCATE(ZLVOCPEXN)        
  DEALLOCATE(ZLSOCPEXN)
  DEALLOCATE(ZTHL) 
  DEALLOCATE(ZRT)
  DEALLOCATE(ZBUF)
! Coherence test
  IF ((.NOT. LUSERI) ) THEN
    IF (MAXVAL(XRT(:,:,:,4))/= 0) THEN
       WRITE(NLUOUT,FMT=*) "*********************************"             
       WRITE(NLUOUT,FMT=*) 'WARNING'      
       WRITE(NLUOUT,FMT=*) 'YOU HAVE LUSERI=FALSE '
       WRITE(NLUOUT,FMT=*) ' BUT WITH YOUR RADIOSOUNDING Ri/=0'
       WRITE(NLUOUT,FMT=*) MINVAL(XRT(:,:,:,4)),MAXVAL(XRT(:,:,:,4))
       WRITE(NLUOUT,FMT=*) "*********************************"       
    ENDIF  
  ENDIF
  IF ((.NOT. LUSERC)) THEN
    IF (MAXVAL(XRT(:,:,:,2))/= 0) THEN          
      WRITE(NLUOUT,FMT=*) "*********************************"
      WRITE(NLUOUT,FMT=*) 'WARNING'    
      WRITE(NLUOUT,FMT=*) 'YOU HAVE LUSERC=FALSE '
      WRITE(NLUOUT,FMT=*) 'BUT WITH YOUR RADIOSOUNDING RC/=0'
      WRITE(NLUOUT,FMT=*) MINVAL(XRT(:,:,:,2)),MAXVAL(XRT(:,:,:,2))      
      WRITE(NLUOUT,FMT=*) "*********************************"
    ENDIF  
  ENDIF
      ! on remet les bonnes valeurs pour NRR
  IF(CCLOUD=='NONE') NRR=1
  IF(CCLOUD=='REVE') NRR=2
END IF
!
!-------------------------------------------------------------------------------
!
!*  	 6.    INITIALIZE SCALAR VARIABLES FOR CHEMISTRY
!   	       -----------------------------------------
!
!  before calling chemistry
CCONF = 'START'
CSTORAGE_TYPE='TT'                  
CALL IO_File_close(TZEXPREFILE)  ! Close the EXPRE file
!
IF ( LCH_INIT_FIELD ) CALL CH_INIT_FIELD_n(1, NLUOUT, NVERB)
!
! Initialization LIMA variables by ORILAM 
IF (CCLOUD == 'LIMA' .AND. ((LORILAM).OR.(LDUST).OR.(LSALT))) &
    CALL AER2LIMA(XSVT, XRHODREF, XRT(:,:,:,1), XPABST, XTHT, XZZ)
!-------------------------------------------------------------------------------
!
!*  	 7.    INITIALIZE LEVELSET FOR IBM
!   	       ---------------------------
!
IF (LIBM_LSF) THEN
  !
  ! In their current state, the IBM can only be used in
  ! combination with cartesian coordinates and flat orography.
  !
  IF ((CZS.NE."FLAT").OR.(.NOT.LCARTESIAN)) THEN
    CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_IDEAL_CASE','IBM can only be used with flat ground')
  ENDIF
  !
  ALLOCATE(XIBM_LS(NIU,NJU,NKU,4))
  !
  CALL IBM_INIT_LS(XIBM_LS)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!*   	 8.    WRITE THE FMFILE 
!   	       ----------------
!
CALL SECOND_MNH2(ZTIME1)
!
NNPRAR = 22 + 2*(NRR+NSV)   &    ! 22 = number of grid variables + reference 
       + 8 + 17                  ! state variables + dimension variables
                                 ! 2*(8+NRR+NSV) + 1 = number of prognostic
                                 ! variables at time t and t-dt
NTYPE=1
!
CALL IO_File_add2list(TINIFILE,TRIM(CINIFILE),'MNH','WRITE',KLFINPRAR=NNPRAR,KLFITYPE=NTYPE,KLFIVERB=NVERB)
!
CALL IO_File_open(TINIFILE)
!
CALL IO_Header_write(TINIFILE)
!
CALL WRITE_DESFM_n(1,TINIFILE)
!
CALL WRITE_LFIFM_n(TINIFILE,'')  ! There is no DAD model for PREP_IDEAL_CASE
!
CALL SECOND_MNH2(ZTIME2)
!
XT_STORE = XT_STORE + ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*     9.     EXTERNALIZED SURFACE
!             --------------------
!
!
IF (CSURF =='EXTE') THEN
  IF (LEN_TRIM(CINIFILEPGD)==0) THEN
    IF (LEN_TRIM(CPGD_FILE)/=0) THEN
      CINIFILEPGD=CPGD_FILE
    ELSE
      !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','PREP_IDEAL_CASE','CINIFILEPGD needed in NAM_LUNITn')
    ENDIF      
  ENDIF
  CALL SURFEX_ALLOC_LIST(1)
  YSURF_CUR => YSURF_LIST(1)
  CALL READ_ALL_NAMELISTS(YSURF_CUR,'MESONH','PRE',.FALSE.)              
  ! Switch to model 1 surface variables
  CALL GOTO_SURFEX(1)
  !* definition of physiographic fields
  ! computed ...
  IF (LEN_TRIM(CPGD_FILE)==0 .OR. .NOT. LREAD_GROUND_PARAM) THEN
    TPGDFILE => TINIFILE
    CALL PGD_GRID_SURF_ATM(YSURF_CUR%UG, YSURF_CUR%U,YSURF_CUR%GCP,'MESONH',TINIFILE%CNAME,'MESONH',.TRUE.,HDIR='-')
    CALL PGD_SURF_ATM     (YSURF_CUR,'MESONH',TINIFILE%CNAME,'MESONH',.TRUE.)
    CALL IO_File_add2list(TINIFILEPGD,TRIM(CINIFILEPGD),'PGD','WRITE',KLFINPRAR=NNPRAR,KLFITYPE=NTYPE,KLFIVERB=NVERB)
    CALL IO_File_open (TINIFILEPGD)
    TPGDFILE => TINIFILEPGD
  ELSE
  ! ... or read from file.
    CALL INIT_PGD_SURF_ATM( YSURF_CUR, 'MESONH', 'PGD',               &
                            '                            ', '      ', &
                            TDTCUR%nyear, TDTCUR%nmonth,              &
                            TDTCUR%nday, TDTCUR%xtime                 )
!
  END IF
  !
  !* forces orography from atmospheric file
  IF (.NOT. LREAD_ZS) CALL MNHPUT_ZS_n
  !
  ! on ecrit un nouveau fichier PGD que s'il n'existe pas
  IF (LEN_TRIM(CPGD_FILE)==0 .OR. .NOT. LREAD_GROUND_PARAM) THEN
    !* writing of physiographic fields in the file
    CSTORAGE_TYPE='PG'
    !
    CALL IO_Header_write(TINIFILEPGD)
    CALL IO_Field_write(TINIFILEPGD,'JPHEXT', JPHEXT)    
    CALL IO_Field_write(TINIFILEPGD,'SURF','EXTE')
    CALL IO_Field_write(TINIFILEPGD,'L1D', L1D)
    CALL IO_Field_write(TINIFILEPGD,'L2D', L2D)
    CALL IO_Field_write(TINIFILEPGD,'PACK',LPACK)
    CALL WRITE_HGRID(1,TINIFILEPGD)
    !
    TOUTDATAFILE => TINIFILEPGD
    !
    TFILE_SURFEX => TINIFILEPGD
    ALLOCATE(YSURF_CUR%DUO%CSELECT(0))    
    CALL WRITE_PGD_SURF_ATM_n(YSURF_CUR,'MESONH')
    NULLIFY(TFILE_SURFEX)
    CSTORAGE_TYPE='TT'
  ENDIF
  !
  !
  !* rereading of physiographic fields and definition of prognostic fields
  !* writing of all surface fields
  TOUTDATAFILE => TINIFILE
  TFILE_SURFEX => TINIFILE
  CALL PREP_SURF_MNH('                            ','      ')
  NULLIFY(TFILE_SURFEX)
ELSE
  CSURF = "NONE"
END IF
!
!-------------------------------------------------------------------------------
!
!*     10.     CLOSES THE FILE
!             ---------------
!
IF (CSURF =='EXTE' .AND. (LEN_TRIM(CPGD_FILE)==0 .OR. .NOT. LREAD_GROUND_PARAM)) THEN
  CALL IO_File_close(TINIFILEPGD)
ENDIF
CALL IO_File_close(TINIFILE)
IF( LEN_TRIM(CPGD_FILE) /= 0 ) THEN
  CALL IO_File_close(TPGDFILE)
ENDIF
!
!
!-------------------------------------------------------------------------------
!
!*      11.    PRINTS ON OUTPUT-LISTING
!              ------------------------
!
IF (NVERB >= 5) THEN
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: LCARTESIAN,CIDEAL,CZS=', &
                                    LCARTESIAN,CIDEAL,CZS 
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: LUSERV=',LUSERV
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: XLON0,XLAT0,XBETA,XRPK,XLONORI,XLATORI=', &
                                    XLON0,XLAT0,XBETA,XRPK,XLONORI,XLATORI
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: XDELTAX,XDELTAY=',XDELTAX,XDELTAY
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: NVERB=',NVERB
  IF(LCARTESIAN) THEN
    WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: No map projection used.'
  ELSE
    IF (XRPK == 1.) THEN
      WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: Polar stereo used.'
    ELSE IF (XRPK == 0.) THEN
      WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: Mercator used.'
    ELSE
      WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: Lambert used, cone factor=',XRPK 
    END IF
  END IF
END IF
!
IF (NVERB >= 5) THEN
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: IIB, IJB, IKB=',NIB,NJB,NKB
  WRITE(NLUOUT,FMT=*) 'PREP_IDEAL_CASE: IIU, IJU, IKU=',NIU,NJU,NKU
END IF
!
!
!*       28.1   print statistics!
!
  !
  CALL SECOND_MNH2(ZTIME2)
  XT_START=XT_START+ZTIME2-ZEND
  !
  ! Set File Timing OUTPUT
  !
  CALL SET_ILUOUT_TIMING(TLUOUT0)
  !
  ! Compute global time
  !
  CALL TIME_STAT_ll(XT_START,ZTOT)
  !
  !
  IMI = 1
  CALL TIME_HEADER_ll(IMI)
  !
  CALL TIME_STAT_ll(XT_STORE,ZTOT,      ' STORE-FIELDS','=')
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')  
  WRITE(YMI,FMT="(I0)") IMI
  CALL TIME_STAT_ll(XT_START,ZTOT,      ' MODEL'//YMI,'+')
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')
  CALL  TIMING_SEPARATOR('+')
WRITE(NLUOUT,FMT=*) ' '
WRITE(NLUOUT,FMT=*) '****************************************************'
WRITE(NLUOUT,FMT=*) '* PREP_IDEAL_CASE: PREP_IDEAL_CASE ENDS CORRECTLY. *'
WRITE(NLUOUT,FMT=*) '****************************************************'
!
CALL FINALIZE_MNH()
!
!
CONTAINS
INCLUDE "th_r_from_thl_rt.func.h"
INCLUDE "compute_frac_ice.func.h"
END PROGRAM PREP_IDEAL_CASE
