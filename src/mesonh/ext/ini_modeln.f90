!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_INI_MODEL_n
!     #######################
!
INTERFACE
!
       SUBROUTINE INI_MODEL_n(KMI,TPINIFILE)
!
USE MODD_IO, ONLY: TFILEDATA
!
INTEGER,          INTENT(IN)   :: KMI       ! Model Index
TYPE(TFILEDATA),  INTENT(IN)   :: TPINIFILE ! Initial file
!
END SUBROUTINE INI_MODEL_n
!
END INTERFACE
!
END MODULE MODI_INI_MODEL_n
!     ############################################
      SUBROUTINE INI_MODEL_n(KMI,TPINIFILE)
!     ############################################
!
!!****  *INI_MODEL_n* - routine to initialize the nested model _n
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize the variables
!     of the nested model _n.
!
!!**  METHOD
!!    ------
!!      The initialization of the model _n is performed as follows :
!!       - Memory for arrays are then allocated :
!!            * If  turbulence kinetic energy variable is not needed
!!    (CTURB='NONE'),  XTKET, XTKEM and XTKES are zero-size arrays.
!!            * If  dissipation of TKE variable is not needed
!!    (CTURBLEN /='KEPS'),  XEPST, XEPSM and XREPSS are zero-size arrays.
!!            * Memory for mixing ratio arrays is allocated according to the
!!     value of logicals LUSERn (the number NRR of moist variables is deduced).
!!            * The latitude (XLAT), longitude (XLON) and map factor (XMAP)
!!     arrays are zero-size arrays if Cartesian geometry (LCARTESIAN=.TRUE.)
!!            * Memory for reference state without orography ( XRHODREFZ and
!!     XTHVREFZ) is only allocated in INI_MODEL1
!!            * The horizontal Coriolis parameters (XCORIOX and XCORIOY) arrays
!!     are  zero-size arrays if thinshell approximation (LTHINSHELL=.TRUE.)
!!            * The Curvature coefficients (XCURVX and XCURVY) arrays
!!     are zero-size arrays if Cartesian geometry (LCARTESIAN=.TRUE.)
!!            * Memory for the Jacobian (ZJ) local array is allocated
!!     (This variable is computed in SET_GRID and used in SET_REF).
!!       - The spatial and temporal grid variables are initialized by SET_GRID.
!!       - The metric coefficients are computed by METRICS (they are using in
!!     the SET-REF call).
!!       - The prognostic variables and are read in initial
!!    LFIFM file (in READ_FIELD)
!!       - The reference state variables  are initialized by SET_REF.
!!       - The temporal indexes of the outputs are computed by SET_OUTPUT_TIMES
!!       - The large scale sources are computed in case of coupling case by
!!    INI_CPL.
!!       - The initialization of the parameters needed for the dynamics
!!         of the model n is realized in INI_DYNAMICS.
!!       - Then the initial file (DESFM+LFIFM files) is closed by IO_File_close.
!!       - The initialization of the parameters needed for the ECMWF radiation
!!         code is realized in INI_RADIATIONS.
!!       - The contents of the scalar variables are overwritten by
!!         the chemistry initialization subroutine CH_INIT_FIELDn when
!!         the flags LUSECHEM and LCH_INIT_FIELD are set to TRUE.
!!         This allows easy initialization of the chemical fields at a
!!         restart of the model.
!!
!!    EXTERNAL
!!    --------
!!      SET_DIM     : to initialize dimensions
!!      SET_GRID    : to initialize grid
!!      METRICS     : to compute metric coefficients
!!      READ_FIELD  : to initialize field
!!      FMCLOS      : to close a FM-file
!!      SET_REF     : to initialize reference state for anelastic approximation
!!      INI_DYNAMICS: to initialize parameters for the dynamics
!!      INI_TKE_EPS : to initialize the TKE
!!      SET_DIRCOS  : to compute the director cosinus of the orography
!!      INI_RADIATIONS : to initialize radiation computations
!!      CH_INIT_CCS: to initialize the chemical core system
!!      CH_INIT_FIELDn: to (re)initialize the scalar variables
!!      INI_DEEP_CONVECTION : to initialize the deep convection scheme
!!      CLEANLIST_ll : deaalocate a list
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_PARAMETERS : contains declaration of parameter variables
!!         JPHEXT : Horizontal external points number
!!         JPVEXT : Vertical external points number
!!
!!      Module MODD_MODD_DYN   : contains declaration of parameters
!!                               for the dynamics
!!      Module MODD_CONF       : contains declaration of configuration variables
!!                               for all models
!!         NMODEL     : Number of nested models
!!         NVERB      : Level of informations on output-listing
!!                          0 for minimum  prints
!!                          5 for intermediate level of prints
!!                         10 for maximum  prints
!!
!!      Module MODD_REF        : contains declaration of reference state
!!                               variables for all models
!!      Module MODD_FIELD_n    : contains declaration of prognostic fields
!!      Module MODD_LSFIELD_n  : contains declaration of Larger Scale fields
!!      Module MODD_GRID_n     : contains declaration of spatial grid variables
!!      Module MODD_TIME_n     : contains declaration of temporal grid variables
!!      Module MODD_REF_n      : contains declaration of reference state
!!                               variables
!!      Module MODD_CURVCOR_n  : contains declaration of curvature and Coriolis
!!                               variables
!!      Module MODD_BUDGET     : contains declarations of the budget parameters
!!      Module MODD_RADIATIONS_n:contains declaration of the variables of the
!!                               radiation interface scheme
!!      Module MODD_STAND_ATM  : contains declaration of the 5 standard
!!                               atmospheres used for the ECMWF-radiation code
!!      Module MODD_FRC        : contains declaration of the control variables
!!                               and of the forcing fields
!!      Module MODD_CH_MNHC_n   : contains the control parameters for chemistry
!!      Module MODD_DEEP_CONVECTION_n: contains declaration of the variables of
!!                                     the deep convection scheme
!!
!!
!!
!!
!!      Module MODN_CONF_n     : contains declaration of namelist NAM_CONFn and
!!                             uses module MODD_CONF_n (configuration variables)
!!      Module MODN_LUNIT_n    : contains declaration of namelist NAM_LUNITn and
!!                             uses module MODD_LUNIT_n (Logical units)
!!      Module MODN_DYN_n      : contains declaration of namelist NAM_DYNn and
!!                             uses module MODD_DYN_n (control of dynamics)
!!      Module MODN_PARAM_n    : contains declaration of namelist NAM_PARAMn and
!!                             uses module MODD_PARAM_n (control of physical
!!                             parameterization)
!!      Module MODN_LBC_n      : contains declaration of namelist NAM_LBCn and
!!                             uses module MODD_LBC_n (lateral boundaries)
!!      Module MODN_TURB_n     : contains declaration of namelist NAM_TURBn and
!!                             uses module MODD_TURB_n (turbulence scheme)
!!      Module MODN_PARAM_RAD_n: contains declaration of namelist NAM_PARAM_RADn
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine INI_MODEL_n)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     10/06/94
!!      Modification 17/10/94  (Stein)  For LCORIO
!!      Modification 20/10/94  (Stein)  For SET_GRID and NAMOUTN
!!      Modification 26/10/94  (Stein)  Modifications of the namelist names
!!      Modification 10/11/94  (Lafore) allocatation of tke fields
!!      Modification 22/11/94  (Stein)  change the READ_FIELDS call ( add
!!                                      pressure function
!!      Modification 06/12/94  (Stein)  add the LS fields
!!                   12/12/94  (Stein)  rename END_INI in INI_DYNAMICS
!!      Modification 09/01/95  (Stein)  add the turbulence scheme
!!      Modification Jan 19, 1995 (J. Cuxart) add the TKE initialization
!!                   Jan 23, 1995 (J. Stein ) remove the condition
!!                             LTHINSHELL=T LCARTESIAN=T => stop
!!      Modification Feb 16, 1995 (I.Mallet) add the METRICS call and
!!                                      change the SET_REF call (add
!!                                      the lineic mass)
!!      Modification Mar 10, 1995 (I. Mallet) add the COUPLING initialization
!!                   June 29,1995 (Ph. Hereil, J. Stein) add the budget init.
!!      Modification Sept. 1, 1995 (S. Belair) Reading of the surface variables
!!                                 and parameters for ISBA (i.e., add a
!!                                 CALL READ_GR_FIELD)
!!      Modification 18/08/95     (J.P.Lafore)   time step change case
!!                   25/09/95     (J. Cuxart and J.Stein)   add LES variables
!!                                and the diachronic file initialization
!!      Modification Sept 20,1995 (Lafore) coupling for the dry mass Md
!!      Modification Sept. 12, 1995 (J.-P. Pinty) add the initialization of
!!                                      the ECMWF radiation code
!!      Modification Sept. 13, 1995 (J.-P. Pinty) control the allocation of the
!!                                      arrays of MODD_GR_FIELD_n
!!      Modification Nove. 17, 1995 (J.Stein) control of the control !!
!!                   March 01, 1996 (J. Stein) add the cloud fraction
!!                   April 03, 1996 (J. Stein) unify the ISBA and TSZ0 cases
!!      Modification 13/12/95 (M. Georgelin) add the forcing variables in
!!                                           the call read_field, and their
!!                                           allocation.
!!                   Mai   23, 1996 (J. Stein) allocate XSEA in the TSZ0 case
!!                   June  11, 1996 (V. Masson) add XSILT and XLAKE of
!!                                              MODD_GR_FIELD_n
!!                   August 7, 1996 (K. Suhre)  add (re)initialization of
!!                                              chemistry
!!                   Octo. 11, 1996 (J. Stein ) add XSRCT and XSRCM
!!                   October 8, 1996 (J. Cuxart, E. Sanchez) Moist LES diagnostics
!!                                     and control on TKE initialization.
!!      Modification 19/12/96 (J.-P. Pinty) add the ice parameterization and
!!                                          the precipitation fields
!!      Modification 11/01/97 (J.-P. Pinty) add the deep convection
!!                   Nov.   1, 1996 (V. Masson) Read the vertical grid kind
!!                   Nov.  20, 1996 (V. Masson) control of convection calling time
!!                   July  16, 1996 (J.P.Lafore) update of EXSEG file reading
!!                   Oct.  08, 1996 (J.P.Lafore, V.Masson)
!!                                       MY_NAME and DAD_NAME reading and check
!!                   Oct.  30, 1996 (J.P.Lafore) resolution ratio reading for nesting
!!                                       and Bikhardt interpolation coef. initialization
!!                   Nov.  22, 1996 (J.P.Lafore) allocation of LS sources for nesting
!!                   Feb.  26, 1997 (J.P.Lafore) allocation of "surfacic" LS fields
!!                   March 10, 1997 (J.P.Lafore) forcing only for model 1
!!                   June  22, 1997 (J. Stein)   add the absolute pressure
!!                   July  09, 1997 (V. Masson)  add directional z0 and SSO
!!                   Aug.  18, 1997 (V. Masson)  consistency between storage
!!                                               type and CCONF
!!                   Dec.  22, 1997 (J. Stein)   add the LS field spawning
!!                   Jan.  24, 1998 (P.Bechtold) change MODD_FRC and MODD_DEEP_CONVECTION
!!                   Dec.  24, 1997 (V.Masson)   directional z0 parameters
!!                   Aug.  13, 1998 (V. Ducrocq P Jabouille)   //
!!                   Mai.  26, 1998 (J. Stein)   remove NXEND,NYEND
!!                   Feb.   1, 1999 (J. Stein)   compute the Bikhardt
!!                                       interpolation coeff. before the call to set_grid
!!                   April  5, 1999 (V. Ducrocq) change the DXRATIO_ALL init.
!!                   April  12, 1999 (J. Stein)  cleaning + INI_SPAWN_LS
!!                   Apr.   7, 1999 (P Jabouille) store the metric coefficients
!!                                                in modd_metrics_n
!!                   Jui.   15,1999 (P Jabouille) split the routines in two parts
!!                   Jan.   04,2000 (V. Masson)   removes the TSZ0 case
!!                   Apr.   15,2000 (P Jabouille) parallelization of grid nesting
!!                   Aug.   20,2000 (J Stein    ) tranpose XBFY
!!                   Jui    01,2000 (F.solmon )   adapatation for patch approach
!!                   Jun.   15,2000 (J.-P. Pinty) add C2R2 initialization
!!                   Nov.  15,2000 (V.Masson) use of ini_modeln in prep_real_case
!!                   Nov.  15,2000 (V.Masson) call of LES routines
!!                   Nov.  15,2000 (V.Masson) aircraft and balloon initialization routines
!!                   Jan.  22,2001 (D.Gazen) update_nsv set NSV_* var. for current model
!!                   Mar.  04,2002 (V.Ducrocq) initialization to temporal series
!!                   Mar.  15,2002 (F.Solmon) modification of ini_radiation interface
!!                   Nov.  29,2002 (JP Pinty) add C3R5, ICE2, ICE4, ELEC
!!                   Jan.  2004    (V.Masson) externalization of surface
!!                   May   2006    Remove KEPS
!!                   Apr.  2010    (M. Leriche) add pH for aqueous phase chemistry
!!                   Jul.  2010    (M. Leriche) add Ice phase chemistry
!!                   Oct.  2010  (J.Escobar) check if local domain not to small for NRIMX NRIMY
!!                   Nov.  2010  (J.Escobar) PGI BUG , add SIZE(CSV) to init_ground routine
!!                   Nov.  2009    (C. Barthe) add call to INI_ELEC_n
!!                   Mar.  2010    (M. Chong) add small ions
!!                   Apr.  2011    (M. Chong) correction of RESTART (ELEC)
!!                   June  2011  (B.Aouizerats) Prognostic aerosols
!!                   June  2011  (P.Aumond) Drag of the vegetation
!!                                         + Mean fields
!!                   July  2013  (Bosseur & Filippi) Adds Forefire
!!                   P. Tulet      Nov 2014 accumulated moles of aqueous species that fall at the surface
!!                   JAn.  2015  (F. Brosse) bug in allocate XACPRAQ
!!                   Dec 2014 (C.Lac) : For reproducibility START/RESTA
!!                   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!                   V. Masson     Feb 2015 replaces, for aerosols, cover fractions by sea, town, bare soil fractions
!!                   J.Escobar : 19/04/2016 : Pb IOZ/NETCDF , missing OPARALLELIO=.FALSE. for PGD files
!!                   J.Escobar : 01/06/2016 : correct check limit of NRIM versus local subdomain size IDIM
!!                   06/2016     (G.Delautier) phasage surfex 8
!!                   Modification    01/2016  (JP Pinty) Add LIMA
!!                   Aug.  2016 (J.Pianezze) Add SFX_OASIS_READ_NAM function from SurfEx
!!                   M.Leriche 2016 Chemistry
!!                   10/2016 M.Mazoyer New KHKO output fields
!!                      10/2016 (C.Lac) Add max values
!!       F. Brosse   Oct.  2016 add prod/loss terms computation for chemistry
!!                   M.Leriche 2016 Chemistry
!!                   M.Leriche 10/02/17 prevent negative values in LBX(Y)SVS
!!                   M.Leriche 01/07/2017 Add DIAG chimical surface fluxes
!!                   09/2017 Q.Rodier add LTEND_UV_FRC
!!                   02/2018 Q.Libois ECRAD
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                   V. Vionnet : 18/07/2017 : add blowing snow scheme
!!                   01/18 J.Colin Add DRAG
!  P. Wautelet 29/01/2019: bug: add missing zero-size allocations
!  P. Wautelet 07/02/2019: force TYPE to a known value for IO_File_add2list
!  P. Wautelet 13/02/2019: initialize XALBUV even if no radiation (needed in CH_INTERP_JVALUES)
!  P. Wautelet 13/02/2019: removed PPABSM and PTSTEP dummy arguments of READ_FIELD
!  P. Wautelet 14/02/2019: remove CLUOUT/CLUOUT0 and associated variables
!  P. Wautelet 14/02/2019: remove HINIFILE dummy argument from INI_RADIATIONS_ECMWF/ECRAD
!!                   02/2019 C.Lac add rain fraction as an output field
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 14/03/2019: correct ZWS when variable not present in file (set to XZWS_DEFAULT)
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 19/04/2019: removed unused dummy arguments and variables
!  P. Wautelet 07/06/2019: allocate lookup tables for optical properties only when needed
!  P. Wautelet 13/09/2019: budget: simplify and modernize date/time management
!  C. Lac         11/2019: correction in the drag formula and application to building in addition to tree
!  S. Riette      04/2020: XHL* fields
!  F. Auguste     02/2021: add IBM
!  T.Nigel        02/2021: add turbulence recycling
! J.L.Redelsperger 06/2011: OCEAN case
! A. Costes       12/2021: Blaze fire model
!---------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
#ifdef MNH_ECRAD
USE YOERDI,                 only: RCCO2
#endif

USE MODD_2D_FRC
USE MODD_ADVFRC_n
USE MODD_ADV_n
use MODD_AEROSET,           only: POLYTAU, POLYSSA, POLYG
USE MODD_ARGSLIST_ll,       only: LIST_ll
USE MODD_BIKHARDT_n
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
USE MODD_BUDGET
USE MODD_CH_AERO_n,         only: XSOLORG,XMI
USE MODD_CH_AEROSOL,        only: LORILAM
USE MODD_CH_BUDGET_n
USE MODD_CH_FLX_n,          only: XCHFLX
USE MODD_CH_M9_n,           only:NNONZEROTERMS
USE MODD_CH_MNHC_n,         only: LUSECHEM, LUSECHAQ, LUSECHIC, LCH_INIT_FIELD, &
                                  LCH_CONV_LINOX, XCH_TUV_DOBNEW, LCH_PH
USE MODD_CH_PH_n
USE MODD_CH_PRODLOSSTOT_n
USE MODD_CLOUD_MF_n
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST
USE MODD_CTURB
USE MODD_CURVCOR_n
USE MODD_DEEP_CONVECTION_n
USE MODD_DEF_EDDY_FLUX_n   ! for VT and WT fluxes
USE MODD_DEF_EDDYUV_FLUX_n ! FOR UV
USE MODD_DIAG_FLAG,         only: LCHEMDIAG, CSPEC_BU_DIAG
USE MODD_DIM_n
USE MODD_DRAG_n
USE MODD_DRAGTREE_n
USE MODD_DUST
use MODD_DUST_OPT_LKT,      only: NMAX_RADIUS_LKT_DUST=>NMAX_RADIUS_LKT, NMAX_SIGMA_LKT_DUST=>NMAX_SIGMA_LKT,               &
                                  NMAX_WVL_SW_DUST=>NMAX_WVL_SW,                                                            &
                                  XEXT_COEFF_WVL_LKT_DUST=>XEXT_COEFF_WVL_LKT, XEXT_COEFF_550_LKT_DUST=>XEXT_COEFF_550_LKT, &
                                  XPIZA_LKT_DUST=>XPIZA_LKT, XCGA_LKT_DUST=>XCGA_LKT
USE MODD_DYN
USE MODD_DYN_n
USE MODD_DYNZD
USE MODD_DYNZD_n
USE MODD_ELEC_n,            only: XCION_POS_FW, XCION_NEG_FW
USE MODD_EOL_MAIN
USE MODD_FIELD_n
#ifdef MNH_FOREFIRE
USE MODD_FOREFIRE
USE MODD_FOREFIRE_n
#endif
USE MODD_FRC
USE MODD_FRC_n
USE MODD_GET_n
USE MODD_GRID_n
USE MODD_GRID,              only: XLONORI,XLATORI
USE MODD_IBM_PARAM_n,       only: LIBM, XIBM_IEPS, XIBM_LS, XIBM_XMUT
USE MODD_IO,                only: CIO_DIR, TFILEDATA, TFILE_DUMMY
USE MODD_IO_SURF_MNH,       only: IO_SURF_MNH_MODEL
USE MODD_LATZ_EDFLX
USE MODD_LBC_n,             only: CLBCX, CLBCY
use modd_les
USE MODD_LSFIELD_n
USE MODD_LUNIT_n
USE MODD_MEAN_FIELD
USE MODD_MEAN_FIELD_n
USE MODD_METRICS_n
USE MODD_MNH_SURFEX_n
USE MODD_NESTING,           only: CDAD_NAME, NDAD, NDT_2_WAY, NDTRATIO, NDXRATIO_ALL, NDYRATIO_ALL
USE MODD_NSV
USE MODD_NSV
USE MODD_NUDGING_n,         only: LNUDGING
USE MODD_OCEANH
USE MODD_OUT_n
USE MODD_PARAMETERS
USE MODD_PARAM_KAFR_n
USE MODD_PARAM_MFSHALL_n
USE MODD_PARAM_n
USE MODD_PARAM_RAD_n,       only: CAER, CAOP, CLW
USE MODD_PASPOL
USE MODD_PASPOL_n
USE MODD_PAST_FIELD_n
use modd_precision,         only: LFIINT
USE MODD_RADIATIONS_n
USE MODD_RECYCL_PARAM_n
USE MODD_REF
USE MODD_REF_n
USE MODD_RELFRC_n
use MODD_SALT,              only: LSALT
use MODD_SALT_OPT_LKT,      only: NMAX_RADIUS_LKT_SALT=>NMAX_RADIUS_LKT, NMAX_SIGMA_LKT_SALT=>NMAX_SIGMA_LKT,               &
                                  NMAX_WVL_SW_SALT=>NMAX_WVL_SW,                                                            &
                                  XEXT_COEFF_WVL_LKT_SALT=>XEXT_COEFF_WVL_LKT, XEXT_COEFF_550_LKT_SALT=>XEXT_COEFF_550_LKT, &
                                  XPIZA_LKT_SALT=>XPIZA_LKT, XCGA_LKT_SALT=>XCGA_LKT
USE MODD_SERIES,            only: LSERIES
USE MODD_SHADOWS_n
USE MODD_STAND_ATM,         only: XSTROATM, XSMLSATM, XSMLWATM, XSPOSATM, XSPOWATM
USE MODD_TIME
USE MODD_TIME_n
USE MODD_TURB_CLOUD,        only: NMODEL_CLOUD, CTURBLEN_CLOUD,XCEI
USE MODD_TURB_n
USE MODD_VAR_ll,            only: IP

USE MODE_GATHER_ll
use mode_ini_budget,        only: Budget_preallocate, Ini_budget
USE MODE_INI_ONE_WAY_n
USE MODE_IO
USE MODE_IO_FIELD_READ,     only: IO_Field_read
USE MODE_IO_FILE,           only: IO_File_open
USE MODE_IO_MANAGE_STRUCT,  only: IO_File_add2list
USE MODE_ll
USE MODE_MODELN_HANDLER
USE MODE_MPPDB
USE MODE_MSG
USE MODE_SPLITTINGZ_ll,     only: GET_DIM_EXTZ_ll
USE MODE_TYPE_ZDIFFU
USE MODE_FIELD,             ONLY: INI_FIELD_LIST

USE MODI_CH_AER_MOD_INIT
USE MODI_CH_INIT_BUDGET_n
USE MODI_CH_INIT_FIELD_n
USE MODI_CH_INIT_JVALUES
USE MODI_CH_INIT_PRODLOSSTOT_n
USE MODI_GET_SIZEX_LB
USE MODI_GET_SIZEY_LB
USE MODI_INI_AEROSET1
USE MODI_INI_AEROSET2
USE MODI_INI_AEROSET3
USE MODI_INI_AEROSET4
USE MODI_INI_AEROSET5
USE MODI_INI_AEROSET6
USE MODI_INI_AIRCRAFT_BALLOON
USE MODI_INI_AIRCRAFT_BALLOON
USE MODI_INI_BIKHARDT_n
USE MODI_INI_CPL
USE MODI_INI_DEEP_CONVECTION
USE MODI_INI_DRAG
USE MODI_INI_DYNAMICS
USE MODI_INI_ELEC_n
USE MODI_INI_EOL_ADNR
USE MODI_INI_EOL_ALM
USE MODI_INI_LES_N
USE MODI_INI_LG
USE MODI_INI_LW_SETUP
USE MODI_INI_MICRO_n
USE MODI_INI_POSPROFILER_n
USE MODI_INI_RADIATIONS
USE MODI_INI_RADIATIONS_ECMWF
USE MODI_INI_RADIATIONS_ECRAD
USE MODI_INI_SERIES_N
USE MODI_INI_SPAWN_LS_n
USE MODI_INI_SURF_RAD
USE MODI_INI_SURFSTATION_n
USE MODI_INI_SW_SETUP
USE MODI_INIT_AEROSOL_PROPERTIES
#ifdef MNH_FOREFIRE
USE MODI_INIT_FOREFIRE_n
#endif
USE MODI_INIT_GROUND_PARAM_n
USE MODI_INI_TKE_EPS
USE MODI_METRICS
USE MODI_MNHGET_SURF_PARAM_n
USE MODI_MNHREAD_ZS_DUMMY_n
USE MODI_READ_FIELD
USE MODI_SET_DIRCOS
USE MODI_SET_GRID
USE MODI_SET_REF
#ifdef CPLOASIS
USE MODI_SFX_OASIS_READ_NAM
#endif
USE MODI_SUNPOS_n
USE MODI_SURF_SOLAR_GEOM
USE MODI_UPDATE_METRICS
USE MODI_UPDATE_NSV
#ifdef MNH_ECRAD
#if ( VER_ECRAD == 140 )
USE YOERDI   , ONLY :RCCO2
#endif
#endif
!
USE MODD_FIRE
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
INTEGER,          INTENT(IN)   :: KMI       ! Model Index
TYPE(TFILEDATA),  INTENT(IN)   :: TPINIFILE ! Initial file
!
!*       0.2   declarations of local variables
!
REAL, PARAMETER :: NALBUV_DEFAULT = 0.01 ! Arbitrary low value for XALBUV
!
INTEGER             :: JSV     ! Loop index
INTEGER             :: IRESP   ! Return code of FM routines
INTEGER             :: ILUOUT  ! Logical unit number of output-listing
CHARACTER(LEN=28)   :: YNAME
INTEGER             :: IIU     ! Upper dimension in x direction (local)
INTEGER             :: IJU     ! Upper dimension in y direction (local)
INTEGER             :: IIU_ll  ! Upper dimension in x direction (global)
INTEGER             :: IJU_ll  ! Upper dimension in y direction (global)
INTEGER             :: IKU     ! Upper dimension in z direction
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZJ ! Jacobian
LOGICAL             :: GINIDCONV ! logical switch for the deep convection
                               ! initialization
LOGICAL             :: GINIRAD ! logical switch for the radiation
                               ! initialization
logical             :: gles    ! Logical to determine if LES diagnostics are enabled
!
!
TYPE(LIST_ll), POINTER :: TZINITHALO2D_ll ! pointer for the list of 2D fields
                                      !  which must be communicated in INIT
TYPE(LIST_ll), POINTER :: TZINITHALO3D_ll ! pointer for the list of 3D fields
                                      !  which must be communicated in INIT
!
INTEGER :: IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU     ! dimensions of the
INTEGER :: IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2       ! West-east LB arrays
INTEGER :: IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV     ! dimensions of the
INTEGER :: IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2       ! North-south LB arrays
INTEGER :: IINFO_ll  ! Return code of //routines
INTEGER :: IIY,IJY
INTEGER :: IIU_B,IJU_B
INTEGER :: IIU_SXP2_YP1_Z_ll,IJU_SXP2_YP1_Z_ll,IKU_SXP2_YP1_Z_ll
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZCO2   ! CO2 concentration near the surface
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZSEA   ! sea fraction
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZTOWN  ! town fraction
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZBARE  ! bare soil fraction
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDIR_ALB ! direct albedo
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZSCA_ALB ! diffuse albedo
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZEMIS    ! emissivity
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZTSRAD   ! surface temperature
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZIBM_LS ! LevelSet IBM
!
!
INTEGER, DIMENSION(:,:),ALLOCATABLE :: IINDEX   ! indices of non-zero terms
INTEGER, DIMENSION(:),ALLOCATABLE   :: IIND
INTEGER                             :: JM, JT
!
!------------------------------------------
! Dummy pointers needed to correct an ifort Bug
REAL, DIMENSION(:), POINTER :: DPTR_XZHAT
REAL, DIMENSION(:), POINTER :: DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4
REAL, DIMENSION(:), POINTER :: DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4
REAL, DIMENSION(:), POINTER :: DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4
REAL, DIMENSION(:), POINTER :: DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4
CHARACTER(LEN=4), DIMENSION(:), POINTER :: DPTR_CLBCX,DPTR_CLBCY
INTEGER, DIMENSION(:,:,:), POINTER :: DPTR_NKLIN_LBXU,DPTR_NKLIN_LBYU,DPTR_NKLIN_LBXV,DPTR_NKLIN_LBYV
INTEGER, DIMENSION(:,:,:), POINTER :: DPTR_NKLIN_LBXW,DPTR_NKLIN_LBYW,DPTR_NKLIN_LBXM,DPTR_NKLIN_LBYM
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXU,DPTR_XCOEFLIN_LBYU
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXV,DPTR_XCOEFLIN_LBYV
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXW,DPTR_XCOEFLIN_LBYW
REAL, DIMENSION(:,:,:), POINTER :: DPTR_XCOEFLIN_LBXM,DPTR_XCOEFLIN_LBYM
REAL, DIMENSION(:,:,:),   POINTER :: DPTR_XLBXUM,DPTR_XLBYUM,DPTR_XLBXVM,DPTR_XLBYVM
REAL, DIMENSION(:,:,:),   POINTER :: DPTR_XLBXWM,DPTR_XLBYWM,DPTR_XLBXTHM,DPTR_XLBYTHM
REAL, DIMENSION(:,:,:),   POINTER :: DPTR_XLBXTKEM,DPTR_XLBYTKEM
REAL, DIMENSION(:,:,:,:),   POINTER :: DPTR_XLBXSVM,DPTR_XLBYSVM
REAL, DIMENSION(:,:,:,:), POINTER :: DPTR_XLBXRM,DPTR_XLBYRM
REAL, DIMENSION(:,:,:),   POINTER ::  DPTR_XZZ
REAL, DIMENSION(:,:,:), POINTER ::   DPTR_XLSUM,DPTR_XLSVM,DPTR_XLSWM,DPTR_XLSTHM,DPTR_XLSRVM
REAL, DIMENSION(:,:,:), POINTER ::   DPTR_XLSUS,DPTR_XLSVS,DPTR_XLSWS,DPTR_XLSTHS,DPTR_XLSRVS
REAL, DIMENSION(:,:),   POINTER ::   DPTR_XLSZWSM,DPTR_XLSZWSS
!
INTEGER                         ::  IIB,IJB,IIE,IJE,IDIMX,IDIMY,IMI
! Fire model
INTEGER             :: INBPARAMSENSIBLE, INBPARAMLATENT
!-------------------------------------------------------------------------------
!
!*       0.    PROLOGUE
!              --------
! Compute relaxation coefficients without changing INI_DYNAMICS nor RELAXDEF
!
IF (CCLOUD == 'LIMA') THEN
  LHORELAX_SVC1R3=LHORELAX_SVLIMA
END IF
!
! UPDATE CONSTANTS FOR OCEAN MODEL
IF (LOCEAN) THEN
  XP00=XP00OCEAN
  XTH00=XTH00OCEAN
END IF
!
!
NULLIFY(TZINITHALO2D_ll)
NULLIFY(TZINITHALO3D_ll)
!
!*       1.    RETRIEVE LOGICAL UNIT NUMBER
!              ----------------------------
!
ILUOUT = TLUOUT%NLU
!
!-------------------------------------------------------------------------------
!
!*       2.   END OF READING
!             --------------
!*       2.1  Read number of forcing fields
!
IF (LFORCING) THEN ! Retrieve the number of time-dependent forcings.
  CALL IO_Field_read(TPINIFILE,'FRC',NFRC,IRESP)
  IF ( (IRESP /= 0) .OR. (NFRC <=0) ) THEN
    WRITE(ILUOUT,'(A/A)') &
     "INI_MODEL_n ERROR: you want to read forcing variables from FMfile", &
     "                   but no fields have been found by IO_Field_read"
!callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_MODEL_n','')
  END IF
END IF
!
! Modif PP for time evolving adv forcing
  IF ( L2D_ADV_FRC ) THEN ! Retrieve the number of time-dependent forcings.
    WRITE(ILUOUT,FMT=*) "INI_MODEL_n ENTER ADV_FORCING"
    CALL IO_Field_read(TPINIFILE,'NADVFRC1',NADVFRC,IRESP)
    IF ( (IRESP /= 0) .OR. (NADVFRC <=0) ) THEN
      WRITE(ILUOUT,'(A/A)') &
      "INI_MODELn ERROR: you want to read forcing ADV variables from FMfile", &
      "                   but no fields have been found by IO_Field_read"
    !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_MODEL_n','')
    END IF
    WRITE(ILUOUT,*) 'NADVFRC = ', NADVFRC
END IF
!
IF ( L2D_REL_FRC ) THEN ! Retrieve the number of time-dependent forcings.
    WRITE(ILUOUT,FMT=*) "INI_MODEL_n ENTER REL_FORCING"
    CALL IO_Field_read(TPINIFILE,'NRELFRC1',NRELFRC,IRESP)
    IF ( (IRESP /= 0) .OR. (NRELFRC <=0) ) THEN
      WRITE(ILUOUT,'(A/A)') &
      "INI_MODELn ERROR: you want to read forcing REL variables from FMfile", &
      "                   but no fields have been found by IO_Field_read"
    !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_MODEL_n','')
    END IF
    WRITE(ILUOUT,*) 'NRELFRC = ', NRELFRC
END IF
!*       2.2  Checks the position of vertical absorbing layer
!
IKU=NKMAX+2*JPVEXT
!
ALLOCATE(XZHAT(IKU))
CALL IO_Field_read(TPINIFILE,'ZHAT',XZHAT)
CALL IO_Field_read(TPINIFILE,'ZTOP',XZTOP)
IF (XALZBOT>=XZHAT(IKU) .AND. LVE_RELAX) THEN
  WRITE(ILUOUT,FMT=*) "INI_MODEL_n ERROR: you want to use vertical relaxation"
  WRITE(ILUOUT,FMT=*) "                  but bottom of layer XALZBOT(",XALZBOT,")"
  WRITE(ILUOUT,FMT=*) "                  is upper than model top    (",XZHAT(IKU),")"
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_MODEL_n','')
END IF
IF (LVE_RELAX) THEN
 IF (XALZBOT>=XZHAT(IKU-4) ) THEN
  WRITE(ILUOUT,FMT=*) "INI_MODEL_n WARNING: you want to use vertical relaxation"
  WRITE(ILUOUT,FMT=*) "                    but the layer defined by XALZBOT(",XALZBOT,")"
  WRITE(ILUOUT,FMT=*) "                    contains less than 5 model levels"
 END IF
END IF
DEALLOCATE(XZHAT)
!
!*       2.3  Compute sizes of arrays of the extended sub-domain
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IIU_ll=NIMAX_ll + 2 * JPHEXT
IJU_ll=NJMAX_ll + 2 * JPHEXT
! initialize NIMAX and NJMAX for not updated versions regarding the parallelism
! spawning,...
CALL GET_DIM_PHYS_ll('B',NIMAX,NJMAX)
!
CALL GET_INDICE_ll( IIB,IJB,IIE,IJE)
IDIMX = IIE - IIB + 1
IDIMY = IJE - IJB + 1
!
NRR=0
NRRL=0
NRRI=0
IF (CGETRVT /= 'SKIP' ) THEN
  NRR = NRR+1
  IDX_RVT = NRR
END IF
IF (CGETRCT /= 'SKIP' ) THEN
  NRR = NRR+1
  NRRL = NRRL+1
  IDX_RCT = NRR
END IF
IF (CGETRRT /= 'SKIP' ) THEN
  NRR = NRR+1
  NRRL = NRRL+1
  IDX_RRT = NRR
END IF
IF (CGETRIT /= 'SKIP' ) THEN
  NRR = NRR+1
  NRRI = NRRI+1
  IDX_RIT = NRR
END IF
IF (CGETRST /= 'SKIP' ) THEN
  NRR = NRR+1
  NRRI = NRRI+1
  IDX_RST = NRR
END IF
IF (CGETRGT /= 'SKIP' ) THEN
  NRR = NRR+1
  NRRI = NRRI+1
  IDX_RGT = NRR
END IF
IF (CGETRHT /= 'SKIP' ) THEN
  NRR = NRR+1
  NRRI = NRRI+1
  IDX_RHT = NRR
END IF
IF (NVERB >= 5) THEN
  WRITE (UNIT=ILUOUT,FMT='("THERE ARE ",I2," WATER VARIABLES")') NRR
  WRITE (UNIT=ILUOUT,FMT='("THERE ARE ",I2," LIQUID VARIABLES")') NRRL
  WRITE (UNIT=ILUOUT,FMT='("THERE ARE ",I2," SOLID VARIABLES")') NRRI
END IF
!
!*       2.4  Update NSV and floating indices for the current model
!
!
CALL UPDATE_NSV(KMI)
!
!-------------------------------------------------------------------------------
!
!*       3.    ALLOCATE  MEMORY
!              -----------------
! * Module RECYCL
!
IF (LRECYCL) THEN
!
  NR_COUNT = 0
!
  ALLOCATE(XUMEANW(IJU,IKU,INT(XNUMBELT)))      ; XUMEANW  = 0.0
  ALLOCATE(XVMEANW(IJU,IKU,INT(XNUMBELT)))      ; XVMEANW  = 0.0
  ALLOCATE(XWMEANW(IJU,IKU,INT(XNUMBELT)))      ; XWMEANW  = 0.0
  ALLOCATE(XUMEANN(IIU,IKU,INT(XNUMBELT)))      ; XUMEANN  = 0.0
  ALLOCATE(XVMEANN(IIU,IKU,INT(XNUMBELT)))      ; XVMEANN  = 0.0
  ALLOCATE(XWMEANN(IIU,IKU,INT(XNUMBELT)))      ; XWMEANN  = 0.0
  ALLOCATE(XUMEANE(IJU,IKU,INT(XNUMBELT)))      ; XUMEANE  = 0.0
  ALLOCATE(XVMEANE(IJU,IKU,INT(XNUMBELT)))      ; XVMEANE  = 0.0
  ALLOCATE(XWMEANE(IJU,IKU,INT(XNUMBELT)))      ; XWMEANE  = 0.0
  ALLOCATE(XUMEANS(IIU,IKU,INT(XNUMBELT)))      ; XUMEANS  = 0.0
  ALLOCATE(XVMEANS(IIU,IKU,INT(XNUMBELT)))      ; XVMEANS  = 0.0
  ALLOCATE(XWMEANS(IIU,IKU,INT(XNUMBELT)))      ; XWMEANS  = 0.0
  ALLOCATE(XTBV(IIU,IJU,IKU))                   ; XTBV  = 0.0
ELSE
  ALLOCATE(XUMEANW(0,0,0))
  ALLOCATE(XVMEANW(0,0,0))
  ALLOCATE(XWMEANW(0,0,0))
  ALLOCATE(XUMEANN(0,0,0))
  ALLOCATE(XVMEANN(0,0,0))
  ALLOCATE(XWMEANN(0,0,0))
  ALLOCATE(XUMEANE(0,0,0))
  ALLOCATE(XVMEANE(0,0,0))
  ALLOCATE(XWMEANE(0,0,0))
  ALLOCATE(XUMEANS(0,0,0))
  ALLOCATE(XVMEANS(0,0,0))
  ALLOCATE(XWMEANS(0,0,0))
  ALLOCATE(XTBV   (0,0,0))
END IF
!
!
!*       3.1   Module MODD_FIELD_n
!
IF (LMEAN_FIELD) THEN
!
  MEAN_COUNT = 0
!
  ALLOCATE(XUM_MEAN(IIU,IJU,IKU))      ; XUM_MEAN  = 0.0
  ALLOCATE(XVM_MEAN(IIU,IJU,IKU))      ; XVM_MEAN  = 0.0
  ALLOCATE(XWM_MEAN(IIU,IJU,IKU))      ; XWM_MEAN  = 0.0
  ALLOCATE(XTHM_MEAN(IIU,IJU,IKU))     ; XTHM_MEAN = 0.0
  ALLOCATE(XTEMPM_MEAN(IIU,IJU,IKU))   ; XTEMPM_MEAN = 0.0
  ALLOCATE(XSVT_MEAN(IIU,IJU,IKU))     ; XSVT_MEAN  = 0.0
  IF (CTURB/='NONE') THEN
    ALLOCATE(XTKEM_MEAN(IIU,IJU,IKU))
    XTKEM_MEAN = 0.0
  ELSE
    ALLOCATE(XTKEM_MEAN(0,0,0))
  END IF
  ALLOCATE(XPABSM_MEAN(IIU,IJU,IKU))   ; XPABSM_MEAN = 0.0
!
  ALLOCATE(XU2_M2(IIU,IJU,IKU))      ; XU2_M2  = 0.0
!
  ALLOCATE(XU2_M2(IIU,IJU,IKU))      ; XU2_M2  = 0.0
  ALLOCATE(XV2_M2(IIU,IJU,IKU))      ; XV2_M2  = 0.0
  ALLOCATE(XW2_M2(IIU,IJU,IKU))      ; XW2_M2  = 0.0
  ALLOCATE(XTH2_M2(IIU,IJU,IKU))     ; XTH2_M2 = 0.0
  ALLOCATE(XTEMP2_M2(IIU,IJU,IKU))   ; XTEMP2_M2 = 0.0
  ALLOCATE(XPABS2_M2(IIU,IJU,IKU))   ; XPABS2_M2 = 0.0
!
  IF (LCOV_FIELD) THEN
    ALLOCATE(XUV_MEAN(IIU,IJU,IKU))    ; XUV_MEAN  = 0.0
    ALLOCATE(XUW_MEAN(IIU,IJU,IKU))    ; XUW_MEAN  = 0.0
    ALLOCATE(XVW_MEAN(IIU,IJU,IKU))    ; XVW_MEAN  = 0.0
    ALLOCATE(XWTH_MEAN(IIU,IJU,IKU))   ; XWTH_MEAN  = 0.0
  END IF
!
  ALLOCATE(XUM_MAX(IIU,IJU,IKU))      ; XUM_MAX  = -1.E20
  ALLOCATE(XVM_MAX(IIU,IJU,IKU))      ; XVM_MAX  = -1.E20
  ALLOCATE(XWM_MAX(IIU,IJU,IKU))      ; XWM_MAX  = -1.E20
  ALLOCATE(XTHM_MAX(IIU,IJU,IKU))     ; XTHM_MAX = 0.0
  ALLOCATE(XTEMPM_MAX(IIU,IJU,IKU))   ; XTEMPM_MAX = 0.0
  IF (CTURB/='NONE') THEN
    ALLOCATE(XTKEM_MAX(IIU,IJU,IKU))
    XTKEM_MAX = 0.0
  ELSE
    ALLOCATE(XTKEM_MAX(0,0,0))
  END IF
  ALLOCATE(XPABSM_MAX(IIU,IJU,IKU))   ; XPABSM_MAX = 0.0
ELSE
!
  ALLOCATE(XUM_MEAN(0,0,0))
  ALLOCATE(XVM_MEAN(0,0,0))
  ALLOCATE(XWM_MEAN(0,0,0))
  ALLOCATE(XTHM_MEAN(0,0,0))
  ALLOCATE(XTEMPM_MEAN(0,0,0))
  ALLOCATE(XSVT_MEAN(0,0,0))
  ALLOCATE(XTKEM_MEAN(0,0,0))
  ALLOCATE(XPABSM_MEAN(0,0,0))
!
  ALLOCATE(XU2_M2(0,0,0))
  ALLOCATE(XV2_M2(0,0,0))
  ALLOCATE(XW2_M2(0,0,0))
  ALLOCATE(XTH2_M2(0,0,0))
  ALLOCATE(XTEMP2_M2(0,0,0))
  ALLOCATE(XPABS2_M2(0,0,0))
!
  IF (LCOV_FIELD) THEN
    ALLOCATE(XUV_MEAN(0,0,0))
    ALLOCATE(XUW_MEAN(0,0,0))
    ALLOCATE(XVW_MEAN(0,0,0))
    ALLOCATE(XWTH_MEAN(0,0,0))
  END IF
!
  ALLOCATE(XUM_MAX(0,0,0))
  ALLOCATE(XVM_MAX(0,0,0))
  ALLOCATE(XWM_MAX(0,0,0))
  ALLOCATE(XTHM_MAX(0,0,0))
  ALLOCATE(XTEMPM_MAX(0,0,0))
  ALLOCATE(XTKEM_MAX(0,0,0))
  ALLOCATE(XPABSM_MAX(0,0,0))
END IF
!
IF ((CUVW_ADV_SCHEME(1:3)=='CEN') .AND. (CTEMP_SCHEME == 'LEFR') ) THEN
  ALLOCATE(XUM(IIU,IJU,IKU))
  ALLOCATE(XVM(IIU,IJU,IKU))
  ALLOCATE(XWM(IIU,IJU,IKU))
  ALLOCATE(XDUM(IIU,IJU,IKU))
  ALLOCATE(XDVM(IIU,IJU,IKU))
  ALLOCATE(XDWM(IIU,IJU,IKU))
  IF (CCONF == 'START') THEN
    XUM  = 0.0
    XVM  = 0.0
    XWM  = 0.0
    XDUM  = 0.0
    XDVM  = 0.0
    XDWM  = 0.0
  END IF
ELSE
  ALLOCATE(XUM(0,0,0))
  ALLOCATE(XVM(0,0,0))
  ALLOCATE(XWM(0,0,0))
  ALLOCATE(XDUM(0,0,0))
  ALLOCATE(XDVM(0,0,0))
  ALLOCATE(XDWM(0,0,0))
END IF
!
ALLOCATE(XUT(IIU,IJU,IKU))      ; XUT  = 0.0
ALLOCATE(XVT(IIU,IJU,IKU))      ; XVT  = 0.0
ALLOCATE(XWT(IIU,IJU,IKU))      ; XWT  = 0.0
ALLOCATE(XTHT(IIU,IJU,IKU))     ; XTHT = 0.0
ALLOCATE(XRUS(IIU,IJU,IKU))     ; XRUS = 0.0
ALLOCATE(XRVS(IIU,IJU,IKU))     ; XRVS = 0.0
ALLOCATE(XRWS(IIU,IJU,IKU))     ; XRWS = 0.0
ALLOCATE(XRUS_PRES(IIU,IJU,IKU)); XRUS_PRES = 0.0
ALLOCATE(XRVS_PRES(IIU,IJU,IKU)); XRVS_PRES = 0.0
ALLOCATE(XRWS_PRES(IIU,IJU,IKU)); XRWS_PRES = 0.0
ALLOCATE(XRTHS(IIU,IJU,IKU))    ; XRTHS = 0.0
ALLOCATE(XRTHS_CLD(IIU,IJU,IKU)); XRTHS_CLD = 0.0

IF ( LIBM ) THEN
  ALLOCATE(ZIBM_LS(IIU,IJU,IKU))  ; ZIBM_LS = 0.0
  ALLOCATE(XIBM_XMUT(IIU,IJU,IKU)); XIBM_XMUT = 0.0
ELSE
  ALLOCATE(ZIBM_LS  (0,0,0))
  ALLOCATE(XIBM_XMUT(0,0,0))
END IF

IF ( LRECYCL ) THEN
  ALLOCATE(XFLUCTUNW(IJU,IKU))    ; XFLUCTUNW  = 0.0
  ALLOCATE(XFLUCTVNN(IIU,IKU))    ; XFLUCTVNN  = 0.0
  ALLOCATE(XFLUCTUTN(IIU,IKU))    ; XFLUCTUTN  = 0.0
  ALLOCATE(XFLUCTVTW(IJU,IKU))    ; XFLUCTVTW  = 0.0
  ALLOCATE(XFLUCTUNE(IJU,IKU))    ; XFLUCTUNE  = 0.0
  ALLOCATE(XFLUCTVNS(IIU,IKU))    ; XFLUCTVNS  = 0.0
  ALLOCATE(XFLUCTUTS(IIU,IKU))    ; XFLUCTUTS  = 0.0
  ALLOCATE(XFLUCTVTE(IJU,IKU))    ; XFLUCTVTE  = 0.0
  ALLOCATE(XFLUCTWTW(IJU,IKU))    ; XFLUCTWTW  = 0.0
  ALLOCATE(XFLUCTWTN(IIU,IKU))    ; XFLUCTWTN  = 0.0
  ALLOCATE(XFLUCTWTE(IJU,IKU))    ; XFLUCTWTE  = 0.0
  ALLOCATE(XFLUCTWTS(IIU,IKU))    ; XFLUCTWTS  = 0.0
ELSE
  ALLOCATE(XFLUCTUNW(0,0))
  ALLOCATE(XFLUCTVNN(0,0))
  ALLOCATE(XFLUCTUTN(0,0))
  ALLOCATE(XFLUCTVTW(0,0))
  ALLOCATE(XFLUCTUNE(0,0))
  ALLOCATE(XFLUCTVNS(0,0))
  ALLOCATE(XFLUCTUTS(0,0))
  ALLOCATE(XFLUCTVTE(0,0))
  ALLOCATE(XFLUCTWTW(0,0))
  ALLOCATE(XFLUCTWTN(0,0))
  ALLOCATE(XFLUCTWTE(0,0))
  ALLOCATE(XFLUCTWTS(0,0))
END IF
!
IF (CTURB /= 'NONE') THEN
  ALLOCATE(XTKET(IIU,IJU,IKU))
  ALLOCATE(XRTKES(IIU,IJU,IKU))
  ALLOCATE(XRTKEMS(IIU,IJU,IKU)); XRTKEMS = 0.0
  ALLOCATE(XWTHVMF(IIU,IJU,IKU))
  ALLOCATE(XDYP(IIU,IJU,IKU))
  ALLOCATE(XTHP(IIU,IJU,IKU))
  ALLOCATE(XTR(IIU,IJU,IKU))
  ALLOCATE(XDISS(IIU,IJU,IKU))
  ALLOCATE(XLEM(IIU,IJU,IKU))
  XTKEMIN=XKEMIN
  XCED   =XCEDIS
ELSE
  ALLOCATE(XTKET(0,0,0))
  ALLOCATE(XRTKES(0,0,0))
  ALLOCATE(XRTKEMS(0,0,0))
  ALLOCATE(XWTHVMF(0,0,0))
  ALLOCATE(XDYP(0,0,0))
  ALLOCATE(XTHP(0,0,0))
  ALLOCATE(XTR(0,0,0))
  ALLOCATE(XDISS(0,0,0))
  ALLOCATE(XLEM(0,0,0))
END IF
IF (CTOM == 'TM06') THEN
  ALLOCATE(XBL_DEPTH(IIU,IJU))
ELSE
  ALLOCATE(XBL_DEPTH(0,0))
END IF
IF (LRMC01) THEN
  ALLOCATE(XSBL_DEPTH(IIU,IJU))
ELSE
  ALLOCATE(XSBL_DEPTH(0,0))
END IF
!
ALLOCATE(XPABSM(IIU,IJU,IKU)) ; XPABSM = 0.0
ALLOCATE(XPABST(IIU,IJU,IKU)) ; XPABST = 0.0
!
ALLOCATE(XRT(IIU,IJU,IKU,NRR)) ;     XRT = 0.0
ALLOCATE(XRRS(IIU,IJU,IKU,NRR)) ;    XRRS = 0.0
ALLOCATE(XRRS_CLD(IIU,IJU,IKU,NRR)); XRRS_CLD = 0.0
!
IF (CTURB /= 'NONE' .AND. NRR>1) THEN
  ALLOCATE(XSRCT(IIU,IJU,IKU))
  ALLOCATE(XSIGS(IIU,IJU,IKU))
ELSE
  ALLOCATE(XSRCT(0,0,0))
  ALLOCATE(XSIGS(0,0,0))
END IF
IF (CCLOUD == 'ICE3'.OR.CCLOUD == 'ICE4') THEN
  ALLOCATE(XHLC_HRC(IIU,IJU,IKU))
  ALLOCATE(XHLC_HCF(IIU,IJU,IKU))
  ALLOCATE(XHLI_HRI(IIU,IJU,IKU))
  ALLOCATE(XHLI_HCF(IIU,IJU,IKU))
  XHLC_HRC(:,:,:)=0.
  XHLC_HCF(:,:,:)=0.
  XHLI_HRI(:,:,:)=0.
  XHLI_HCF(:,:,:)=0.
ELSE
  ALLOCATE(XHLC_HRC(0,0,0))
  ALLOCATE(XHLC_HCF(0,0,0))
  ALLOCATE(XHLI_HRI(0,0,0))
  ALLOCATE(XHLI_HCF(0,0,0))
END IF
!
IF (NRR>1) THEN
  ALLOCATE(XCLDFR(IIU,IJU,IKU));  XCLDFR (:, :, :) = 0.
  ALLOCATE(XICEFR(IIU,IJU,IKU));  XICEFR (:, :, :) = 0.
  ALLOCATE(XRAINFR(IIU,IJU,IKU)); XRAINFR(:, :, :) = 0.
ELSE
  ALLOCATE(XCLDFR(0,0,0))
  ALLOCATE(XICEFR(0,0,0))
  ALLOCATE(XRAINFR(0,0,0))
END IF
!
ALLOCATE(XSVT(IIU,IJU,IKU,NSV)) ;     XSVT  = 0.
ALLOCATE(XRSVS(IIU,IJU,IKU,NSV));     XRSVS = 0.
ALLOCATE(XRSVS_CLD(IIU,IJU,IKU,NSV)); XRSVS_CLD = 0.0
ALLOCATE(XZWS(IIU,IJU)) ;             XZWS(:,:) = XZWS_DEFAULT
!
IF (LPASPOL) THEN
  ALLOCATE( XATC(IIU,IJU,IKU,NSV_PP) )
  XATC = 0.
ELSE
  ALLOCATE( XATC(0,0,0,0))
END IF
!
IF(LBLOWSNOW) THEN
  ALLOCATE(XSNWCANO(IIU,IJU,NBLOWSNOW_2D))
  ALLOCATE(XRSNWCANOS(IIU,IJU,NBLOWSNOW_2D))
  XSNWCANO(:,:,:) = 0.0
  XRSNWCANOS(:,:,:) = 0.0
ELSE
  ALLOCATE(XSNWCANO(0,0,0))
  ALLOCATE(XRSNWCANOS(0,0,0))
END IF
!
!*       3.2   Module MODD_GRID_n and MODD_METRICS_n
!
IF (LCARTESIAN) THEN
  ALLOCATE(XLON(0,0))
  ALLOCATE(XLAT(0,0))
  ALLOCATE(XMAP(0,0))
ELSE
  ALLOCATE(XLON(IIU,IJU))
  ALLOCATE(XLAT(IIU,IJU))
  ALLOCATE(XMAP(IIU,IJU))
END IF
ALLOCATE(XXHAT(IIU))
ALLOCATE(XDXHAT(IIU))
ALLOCATE(XYHAT(IJU))
ALLOCATE(XDYHAT(IJU))
ALLOCATE(XZS(IIU,IJU))
ALLOCATE(XZSMT(IIU,IJU))
ALLOCATE(XZZ(IIU,IJU,IKU))
ALLOCATE(XZHAT(IKU))
ALLOCATE(XDIRCOSZW(IIU,IJU))
ALLOCATE(XDIRCOSXW(IIU,IJU))
ALLOCATE(XDIRCOSYW(IIU,IJU))
ALLOCATE(XCOSSLOPE(IIU,IJU))
ALLOCATE(XSINSLOPE(IIU,IJU))
!
ALLOCATE(XDXX(IIU,IJU,IKU))
ALLOCATE(XDYY(IIU,IJU,IKU))
ALLOCATE(XDZX(IIU,IJU,IKU))
ALLOCATE(XDZY(IIU,IJU,IKU))
ALLOCATE(XDZZ(IIU,IJU,IKU))
!
!*       3.3   Modules MODD_REF and  MODD_REF_n
!
! Different reference states for Ocean and Atmosphere models
!  For the moment, same reference states for O and A
!IF ((KMI == 1).OR.LCOUPLES) THEN
IF (KMI==1) THEN
  ALLOCATE(XRHODREFZ(IKU),XTHVREFZ(IKU))
ELSE IF (LCOUPLES) THEN
! in coupled O-A case, need different variables for ocean
  ALLOCATE(XRHODREFZO(IKU),XTHVREFZO(IKU))
ELSE
  !Do not allocate XRHODREFZ and XTHVREFZ because they are the same on all grids (not 'n' variables)
END IF
!
ALLOCATE(XPHIT(IIU,IJU,IKU))
ALLOCATE(XRHODREF(IIU,IJU,IKU))
ALLOCATE(XTHVREF(IIU,IJU,IKU))
ALLOCATE(XEXNREF(IIU,IJU,IKU))
ALLOCATE(XRHODJ(IIU,IJU,IKU))
IF (CEQNSYS=='DUR' .AND. LUSERV) THEN
  ALLOCATE(XRVREF(IIU,IJU,IKU))
ELSE
  ALLOCATE(XRVREF(0,0,0))
END IF
!
!*       3.4   Module MODD_CURVCOR_n
!
IF (LTHINSHELL) THEN
  ALLOCATE(XCORIOX(0,0))
  ALLOCATE(XCORIOY(0,0))
ELSE
  ALLOCATE(XCORIOX(IIU,IJU))
  ALLOCATE(XCORIOY(IIU,IJU))
END IF
  ALLOCATE(XCORIOZ(IIU,IJU))
IF (LCARTESIAN) THEN
  ALLOCATE(XCURVX(0,0))
  ALLOCATE(XCURVY(0,0))
ELSE
  ALLOCATE(XCURVX(IIU,IJU))
  ALLOCATE(XCURVY(IIU,IJU))
END IF
!
!*       3.5   Module MODD_DYN_n
!
CALL GET_DIM_EXT_ll('Y',IIY,IJY)
IF (L2D) THEN
  ALLOCATE(XBFY(IIY,IJY,IKU))
ELSE
  ALLOCATE(XBFY(IJY,IIY,IKU)) ! transposition needed by the optimisation of the
                              ! FFT solver
END IF
CALL GET_DIM_EXT_ll('B',IIU_B,IJU_B)
ALLOCATE(XBFB(IIU_B,IJU_B,IKU))
CALL GET_DIM_EXTZ_ll('SXP2_YP1_Z',IIU_SXP2_YP1_Z_ll,IJU_SXP2_YP1_Z_ll,IKU_SXP2_YP1_Z_ll)
ALLOCATE(XBF_SXP2_YP1_Z(IIU_SXP2_YP1_Z_ll,IJU_SXP2_YP1_Z_ll,IKU_SXP2_YP1_Z_ll))
ALLOCATE(XAF(IKU),XCF(IKU))
ALLOCATE(XTRIGSX(3*IIU_ll))
ALLOCATE(XTRIGSY(3*IJU_ll))
ALLOCATE(XRHOM(IKU))
ALLOCATE(XALK(IKU))
ALLOCATE(XALKW(IKU))
ALLOCATE(XALKBAS(IKU))
ALLOCATE(XALKWBAS(IKU))
!
IF ( LHORELAX_UVWTH .OR. LHORELAX_RV .OR.                                  &
     LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI  .OR. LHORELAX_RS  .OR. &
     LHORELAX_RG .OR. LHORELAX_RH .OR. LHORELAX_TKE .OR.                   &
     ANY(LHORELAX_SV) ) THEN
  ALLOCATE(XKURELAX(IIU,IJU))
  ALLOCATE(XKVRELAX(IIU,IJU))
  ALLOCATE(XKWRELAX(IIU,IJU))
  ALLOCATE(LMASK_RELAX(IIU,IJU))
ELSE
  ALLOCATE(XKURELAX(0,0))
  ALLOCATE(XKVRELAX(0,0))
  ALLOCATE(XKWRELAX(0,0))
  ALLOCATE(LMASK_RELAX(0,0))
END IF
!
! Additional fields for truly horizontal diffusion (Module MODD_DYNZD$n)
IF (LZDIFFU) THEN
  CALL INIT_TYPE_ZDIFFU_HALO2(XZDIFFU_HALO2)
ELSE
  CALL INIT_TYPE_ZDIFFU_HALO2(XZDIFFU_HALO2,0)
ENDIF
!
!*       3.6   Larger Scale variables (Module MODD_LSFIELD$n)
!
!
! upper relaxation part
!
ALLOCATE(XLSUM(IIU,IJU,IKU))    ; XLSUM  = 0.0
ALLOCATE(XLSVM(IIU,IJU,IKU))    ; XLSVM  = 0.0
ALLOCATE(XLSWM(IIU,IJU,IKU))    ; XLSWM  = 0.0
ALLOCATE(XLSTHM(IIU,IJU,IKU))   ; XLSTHM = 0.0
IF ( NRR > 0 ) THEN
  ALLOCATE(XLSRVM(IIU,IJU,IKU)) ; XLSRVM = 0.0
ELSE
  ALLOCATE(XLSRVM(0,0,0))
END IF
ALLOCATE(XLSZWSM(IIU,IJU)) ; XLSZWSM = -1.
!
!  lbc part
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
ELSEIF( L2D ) THEN                         ! 2D case
!
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
  CALL GET_SIZEX_LB(NIMAX_ll,NJMAX_ll,NRIMX,  &
       IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU, &
       IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2)
!
  IF ( LHORELAX_UVWTH ) THEN
    NSIZELBX_ll=2*NRIMX+2*JPHEXT
    NSIZELBXU_ll=2*NRIMX+2*JPHEXT
    ALLOCATE(XLBXUM(IISIZEXFU,IJSIZEXFU,IKU))
    ALLOCATE(XLBXVM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBXWM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBXTHM(IISIZEXF,IJSIZEXF,IKU))
  ELSE
    NSIZELBX_ll=2*JPHEXT      ! 2
    NSIZELBXU_ll=2*(JPHEXT+1) ! 4
    ALLOCATE(XLBXUM(IISIZEX4,IJSIZEX4,IKU))
    ALLOCATE(XLBXVM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBXWM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBXTHM(IISIZEX2,IJSIZEX2,IKU))
  END IF
!
  IF (CTURB /= 'NONE') THEN
    IF ( LHORELAX_TKE) THEN
      NSIZELBXTKE_ll=2* NRIMX+2*JPHEXT
      ALLOCATE(XLBXTKEM(IISIZEXF,IJSIZEXF,IKU))
    ELSE
      NSIZELBXTKE_ll=2*JPHEXT  ! 2
      ALLOCATE(XLBXTKEM(IISIZEX2,IJSIZEX2,IKU))
    END IF
  ELSE
    NSIZELBXTKE_ll=0
    ALLOCATE(XLBXTKEM(0,0,0))
  END IF
  !
  IF ( NRR > 0 ) THEN
    IF (LHORELAX_RV .OR. LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI    &
         .OR. LHORELAX_RS .OR. LHORELAX_RG .OR. LHORELAX_RH               &
       ) THEN
      NSIZELBXR_ll=2* NRIMX+2*JPHEXT
      ALLOCATE(XLBXRM(IISIZEXF,IJSIZEXF,IKU,NRR))
    ELSE
      NSIZELBXR_ll=2*JPHEXT  ! 2
      ALLOCATE(XLBXRM(IISIZEX2,IJSIZEX2,IKU,NRR))
    ENDIF
  ELSE
    NSIZELBXR_ll=0
    ALLOCATE(XLBXRM(0,0,0,0))
  END IF
  !
  IF ( NSV > 0 ) THEN
    IF ( ANY( LHORELAX_SV(:)) ) THEN
      NSIZELBXSV_ll=2* NRIMX+2*JPHEXT
      ALLOCATE(XLBXSVM(IISIZEXF,IJSIZEXF,IKU,NSV))
    ELSE
      NSIZELBXSV_ll=2*JPHEXT  ! 2
      ALLOCATE(XLBXSVM(IISIZEX2,IJSIZEX2,IKU,NSV))
    END IF
  ELSE
    NSIZELBXSV_ll=0
    ALLOCATE(XLBXSVM(0,0,0,0))
  END IF
!
ELSE                                   ! 3D case
!
!
  CALL GET_SIZEX_LB(NIMAX_ll,NJMAX_ll,NRIMX,               &
                    IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU, &
                    IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2)
  CALL GET_SIZEY_LB(NIMAX_ll,NJMAX_ll,NRIMY,               &
                    IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV, &
                    IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2)
!
! check if local domain not to small for NRIMX NRIMY
!
  IF ( CLBCX(1) /= 'CYCL' )  THEN
     IF ( NRIMX .GT. IDIMX )   THEN
        WRITE(*,'(A,I8,A/A,2I8,/A)') "Processor=", IP-1, &
             " :: INI_MODEL_n ERROR:  ( NRIMX  > IDIMX )  ", &
             " Local domain to small for relaxation NRIMX,IDIMX ", &
             NRIMX,IDIMX ,&
             " change relaxation parameters or number of processors "
        call Print_msg(NVERB_FATAL,'GEN','INI_MODEL_n','')
     END IF
  END IF
  IF ( CLBCY(1) /= 'CYCL' ) THEN
     IF ( NRIMY .GT. IDIMY )  THEN
        WRITE(*,'(A,I8,A/A,2I8,/A)') "Processor=", IP-1, &
             " :: INI_MODEL_n ERROR:  ( NRIMY > IDIMY )  ", &
             " Local domain to small for relaxation NRIMY,IDIMY ", &
             NRIMY,IDIMY ,&
             " change relaxation parameters or number of processors "
        call Print_msg(NVERB_FATAL,'GEN','INI_MODEL_n','')
     END IF
  END IF
IF ( LHORELAX_UVWTH ) THEN
    NSIZELBX_ll=2*NRIMX+2*JPHEXT
    NSIZELBXU_ll=2*NRIMX+2*JPHEXT
    NSIZELBY_ll=2*NRIMY+2*JPHEXT
    NSIZELBYV_ll=2*NRIMY+2*JPHEXT
    ALLOCATE(XLBXUM(IISIZEXFU,IJSIZEXFU,IKU))
    ALLOCATE(XLBYUM(IISIZEYF,IJSIZEYF,IKU))
    ALLOCATE(XLBXVM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBYVM(IISIZEYFV,IJSIZEYFV,IKU))
    ALLOCATE(XLBXWM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBYWM(IISIZEYF,IJSIZEYF,IKU))
    ALLOCATE(XLBXTHM(IISIZEXF,IJSIZEXF,IKU))
    ALLOCATE(XLBYTHM(IISIZEYF,IJSIZEYF,IKU))
  ELSE
    NSIZELBX_ll=2*JPHEXT  ! 2
    NSIZELBXU_ll=2*(JPHEXT+1) ! 4
    NSIZELBY_ll=2*JPHEXT  ! 2
    NSIZELBYV_ll=2*(JPHEXT+1) ! 4
    ALLOCATE(XLBXUM(IISIZEX4,IJSIZEX4,IKU))
    ALLOCATE(XLBYUM(IISIZEY2,IJSIZEY2,IKU))
    ALLOCATE(XLBXVM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBYVM(IISIZEY4,IJSIZEY4,IKU))
    ALLOCATE(XLBXWM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBYWM(IISIZEY2,IJSIZEY2,IKU))
    ALLOCATE(XLBXTHM(IISIZEX2,IJSIZEX2,IKU))
    ALLOCATE(XLBYTHM(IISIZEY2,IJSIZEY2,IKU))
  END IF
  !
  IF (CTURB /= 'NONE') THEN
    IF ( LHORELAX_TKE) THEN
      NSIZELBXTKE_ll=2*NRIMX+2*JPHEXT
      NSIZELBYTKE_ll=2*NRIMY+2*JPHEXT
      ALLOCATE(XLBXTKEM(IISIZEXF,IJSIZEXF,IKU))
      ALLOCATE(XLBYTKEM(IISIZEYF,IJSIZEYF,IKU))
    ELSE
      NSIZELBXTKE_ll=2*JPHEXT  ! 2
      NSIZELBYTKE_ll=2*JPHEXT  ! 2
      ALLOCATE(XLBXTKEM(IISIZEX2,IJSIZEX2,IKU))
      ALLOCATE(XLBYTKEM(IISIZEY2,IJSIZEY2,IKU))
    END IF
  ELSE
    NSIZELBXTKE_ll=0
    NSIZELBYTKE_ll=0
    ALLOCATE(XLBXTKEM(0,0,0))
    ALLOCATE(XLBYTKEM(0,0,0))
  END IF
  !
  IF ( NRR > 0 ) THEN
    IF (LHORELAX_RV .OR. LHORELAX_RC .OR. LHORELAX_RR .OR. LHORELAX_RI    &
          .OR. LHORELAX_RS .OR. LHORELAX_RG .OR. LHORELAX_RH              &
       ) THEN
      NSIZELBXR_ll=2*NRIMX+2*JPHEXT
      NSIZELBYR_ll=2*NRIMY+2*JPHEXT
      ALLOCATE(XLBXRM(IISIZEXF,IJSIZEXF,IKU,NRR))
      ALLOCATE(XLBYRM(IISIZEYF,IJSIZEYF,IKU,NRR))
    ELSE
      NSIZELBXR_ll=2*JPHEXT  ! 2
      NSIZELBYR_ll=2*JPHEXT  ! 2
      ALLOCATE(XLBXRM(IISIZEX2,IJSIZEX2,IKU,NRR))
      ALLOCATE(XLBYRM(IISIZEY2,IJSIZEY2,IKU,NRR))
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
      ALLOCATE(XLBXSVM(IISIZEXF,IJSIZEXF,IKU,NSV))
      ALLOCATE(XLBYSVM(IISIZEYF,IJSIZEYF,IKU,NSV))
    ELSE
      NSIZELBXSV_ll=2*JPHEXT  ! 2
      NSIZELBYSV_ll=2*JPHEXT  ! 2
      ALLOCATE(XLBXSVM(IISIZEX2,IJSIZEX2,IKU,NSV))
      ALLOCATE(XLBYSVM(IISIZEY2,IJSIZEY2,IKU,NSV))
    END IF
  ELSE
    NSIZELBXSV_ll=0
    NSIZELBYSV_ll=0
    ALLOCATE(XLBXSVM(0,0,0,0))
    ALLOCATE(XLBYSVM(0,0,0,0))
  END IF
END IF      ! END OF THE IF STRUCTURE ON THE MODEL DIMENSION
!
!
IF ( KMI > 1 ) THEN
  ! it has been assumed that the THeta field used the largest rim area compared
  ! to the others prognostic variables, if it is not the case, you must change
  ! these lines
  ALLOCATE(XCOEFLIN_LBXM(SIZE(XLBXTHM,1),SIZE(XLBXTHM,2),SIZE(XLBXTHM,3)))
  ALLOCATE(   NKLIN_LBXM(SIZE(XLBXTHM,1),SIZE(XLBXTHM,2),SIZE(XLBXTHM,3)))
  ALLOCATE(XCOEFLIN_LBYM(SIZE(XLBYTHM,1),SIZE(XLBYTHM,2),SIZE(XLBYTHM,3)))
  ALLOCATE(   NKLIN_LBYM(SIZE(XLBYTHM,1),SIZE(XLBYTHM,2),SIZE(XLBYTHM,3)))
  ALLOCATE(XCOEFLIN_LBXU(SIZE(XLBXUM,1),SIZE(XLBXUM,2),SIZE(XLBXUM,3)))
  ALLOCATE(   NKLIN_LBXU(SIZE(XLBXUM,1),SIZE(XLBXUM,2),SIZE(XLBXUM,3)))
  ALLOCATE(XCOEFLIN_LBYU(SIZE(XLBYUM,1),SIZE(XLBYUM,2),SIZE(XLBYUM,3)))
  ALLOCATE(   NKLIN_LBYU(SIZE(XLBYUM,1),SIZE(XLBYUM,2),SIZE(XLBYUM,3)))
  ALLOCATE(XCOEFLIN_LBXV(SIZE(XLBXVM,1),SIZE(XLBXVM,2),SIZE(XLBXVM,3)))
  ALLOCATE(   NKLIN_LBXV(SIZE(XLBXVM,1),SIZE(XLBXVM,2),SIZE(XLBXVM,3)))
  ALLOCATE(XCOEFLIN_LBYV(SIZE(XLBYVM,1),SIZE(XLBYVM,2),SIZE(XLBYVM,3)))
  ALLOCATE(   NKLIN_LBYV(SIZE(XLBYVM,1),SIZE(XLBYVM,2),SIZE(XLBYVM,3)))
  ALLOCATE(XCOEFLIN_LBXW(SIZE(XLBXWM,1),SIZE(XLBXWM,2),SIZE(XLBXWM,3)))
  ALLOCATE(   NKLIN_LBXW(SIZE(XLBXWM,1),SIZE(XLBXWM,2),SIZE(XLBXWM,3)))
  ALLOCATE(XCOEFLIN_LBYW(SIZE(XLBYWM,1),SIZE(XLBYWM,2),SIZE(XLBYWM,3)))
  ALLOCATE(   NKLIN_LBYW(SIZE(XLBYWM,1),SIZE(XLBYWM,2),SIZE(XLBYWM,3)))
ELSE
  ALLOCATE(XCOEFLIN_LBXM(0,0,0))
  ALLOCATE(   NKLIN_LBXM(0,0,0))
  ALLOCATE(XCOEFLIN_LBYM(0,0,0))
  ALLOCATE(   NKLIN_LBYM(0,0,0))
  ALLOCATE(XCOEFLIN_LBXU(0,0,0))
  ALLOCATE(   NKLIN_LBXU(0,0,0))
  ALLOCATE(XCOEFLIN_LBYU(0,0,0))
  ALLOCATE(   NKLIN_LBYU(0,0,0))
  ALLOCATE(XCOEFLIN_LBXV(0,0,0))
  ALLOCATE(   NKLIN_LBXV(0,0,0))
  ALLOCATE(XCOEFLIN_LBYV(0,0,0))
  ALLOCATE(   NKLIN_LBYV(0,0,0))
  ALLOCATE(XCOEFLIN_LBXW(0,0,0))
  ALLOCATE(   NKLIN_LBXW(0,0,0))
  ALLOCATE(XCOEFLIN_LBYW(0,0,0))
  ALLOCATE(   NKLIN_LBYW(0,0,0))
END IF
!
!  allocation of the LS fields for vertical relaxation and numerical diffusion
IF( .NOT. LSTEADYLS )  THEN
!
  ALLOCATE(XLSUS(SIZE(XLSUM,1),SIZE(XLSUM,2),SIZE(XLSUM,3)))
  ALLOCATE(XLSVS(SIZE(XLSVM,1),SIZE(XLSVM,2),SIZE(XLSVM,3)))
  ALLOCATE(XLSWS(SIZE(XLSWM,1),SIZE(XLSWM,2),SIZE(XLSWM,3)))
  ALLOCATE(XLSTHS(SIZE(XLSTHM,1),SIZE(XLSTHM,2),SIZE(XLSTHM,3)))
  ALLOCATE(XLSRVS(SIZE(XLSRVM,1),SIZE(XLSRVM,2),SIZE(XLSRVM,3)))
  ALLOCATE(XLSZWSS(SIZE(XLSZWSM,1),SIZE(XLSZWSM,2)))
!
ELSE
!
  ALLOCATE(XLSUS(0,0,0))
  ALLOCATE(XLSVS(0,0,0))
  ALLOCATE(XLSWS(0,0,0))
  ALLOCATE(XLSTHS(0,0,0))
  ALLOCATE(XLSRVS(0,0,0))
  ALLOCATE(XLSZWSS(0,0))
!
END IF
!  allocation of the LB fields for horizontal relaxation and Lateral Boundaries
IF( .NOT. ( LSTEADYLS .AND. KMI==1 ) )  THEN
!
  ALLOCATE(XLBXTKES(SIZE(XLBXTKEM,1),SIZE(XLBXTKEM,2),SIZE(XLBXTKEM,3)))
  ALLOCATE(XLBYTKES(SIZE(XLBYTKEM,1),SIZE(XLBYTKEM,2),SIZE(XLBYTKEM,3)))
  ALLOCATE(XLBXUS(SIZE(XLBXUM,1),SIZE(XLBXUM,2),SIZE(XLBXUM,3)))
  ALLOCATE(XLBYUS(SIZE(XLBYUM,1),SIZE(XLBYUM,2),SIZE(XLBYUM,3)))
  ALLOCATE(XLBXVS(SIZE(XLBXVM,1),SIZE(XLBXVM,2),SIZE(XLBXVM,3)))
  ALLOCATE(XLBYVS(SIZE(XLBYVM,1),SIZE(XLBYVM,2),SIZE(XLBYVM,3)))
  ALLOCATE(XLBXWS(SIZE(XLBXWM,1),SIZE(XLBXWM,2),SIZE(XLBXWM,3)))
  ALLOCATE(XLBYWS(SIZE(XLBYWM,1),SIZE(XLBYWM,2),SIZE(XLBYWM,3)))
  ALLOCATE(XLBXTHS(SIZE(XLBXTHM,1),SIZE(XLBXTHM,2),SIZE(XLBXTHM,3)))
  ALLOCATE(XLBYTHS(SIZE(XLBYTHM,1),SIZE(XLBYTHM,2),SIZE(XLBYTHM,3)))
  ALLOCATE(XLBXRS(SIZE(XLBXRM,1),SIZE(XLBXRM,2),SIZE(XLBXRM,3),SIZE(XLBXRM,4)))
  ALLOCATE(XLBYRS(SIZE(XLBYRM,1),SIZE(XLBYRM,2),SIZE(XLBYRM,3),SIZE(XLBYRM,4)))
  ALLOCATE(XLBXSVS(SIZE(XLBXSVM,1),SIZE(XLBXSVM,2),SIZE(XLBXSVM,3),SIZE(XLBXSVM,4)))
  ALLOCATE(XLBYSVS(SIZE(XLBYSVM,1),SIZE(XLBYSVM,2),SIZE(XLBYSVM,3),SIZE(XLBYSVM,4)))
!
ELSE
!
  ALLOCATE(XLBXTKES(0,0,0))
  ALLOCATE(XLBYTKES(0,0,0))
  ALLOCATE(XLBXUS(0,0,0))
  ALLOCATE(XLBYUS(0,0,0))
  ALLOCATE(XLBXVS(0,0,0))
  ALLOCATE(XLBYVS(0,0,0))
  ALLOCATE(XLBXWS(0,0,0))
  ALLOCATE(XLBYWS(0,0,0))
  ALLOCATE(XLBXTHS(0,0,0))
  ALLOCATE(XLBYTHS(0,0,0))
  ALLOCATE(XLBXRS(0,0,0,0))
  ALLOCATE(XLBYRS(0,0,0,0))
  ALLOCATE(XLBXSVS(0,0,0,0))
  ALLOCATE(XLBYSVS(0,0,0,0))
!
END IF
!
!
!*       3.7   Module MODD_RADIATIONS_n (except XOZON and XAER)
!
! Initialization of SW bands
NSWB_OLD = 6 ! Number of bands in ECMWF original scheme (from Fouquart et Bonnel (1980))
             ! then modified through INI_RADIATIONS_ECMWF but remains equal to 6 practically

#ifdef MNH_ECRAD
#if ( VER_ECRAD == 140 )
NLWB_OLD = 16 ! For XEMIS initialization (should be spectral in the future)
#endif
#endif

NLWB_MNH = 16 ! For XEMIS initialization (should be spectral in the future)

IF (CRAD == 'ECRA') THEN
    NSWB_MNH = 14
#ifdef MNH_ECRAD
#if ( VER_ECRAD == 140 )
    NLWB_MNH = 16
#endif
#endif
ELSE
    NSWB_MNH = NSWB_OLD
#ifdef MNH_ECRAD
#if ( VER_ECRAD == 140 )
    NLWB_MNH = NLWB_OLD
#endif
#endif
END IF

ALLOCATE(XSW_BANDS (NSWB_MNH))
ALLOCATE(XLW_BANDS (NLWB_MNH))
ALLOCATE(XZENITH   (IIU,IJU))
ALLOCATE(XAZIM     (IIU,IJU))
ALLOCATE(XALBUV    (IIU,IJU))
XALBUV(:,:) = NALBUV_DEFAULT !Set to an arbitrary low value (XALBUV is needed in CH_INTERP_JVALUES even if no radiation)
ALLOCATE(XDIRSRFSWD(IIU,IJU,NSWB_MNH))
ALLOCATE(XSCAFLASWD(IIU,IJU,NSWB_MNH))
ALLOCATE(XFLALWD   (IIU,IJU))
!
IF (CRAD /= 'NONE') THEN
  ALLOCATE(XSLOPANG(IIU,IJU))
  ALLOCATE(XSLOPAZI(IIU,IJU))
  ALLOCATE(XDTHRAD(IIU,IJU,IKU))
  ALLOCATE(XDIRFLASWD(IIU,IJU,NSWB_MNH))
  ALLOCATE(XDIR_ALB(IIU,IJU,NSWB_MNH))
  ALLOCATE(XSCA_ALB(IIU,IJU,NSWB_MNH))
  ALLOCATE(XEMIS  (IIU,IJU,NLWB_MNH))
  ALLOCATE(XTSRAD (IIU,IJU))    ; XTSRAD = 0.0
  ALLOCATE(XSEA (IIU,IJU))
  ALLOCATE(XZS_XY (IIU,IJU))
  ALLOCATE(NCLEARCOL_TM1(IIU,IJU))
  ALLOCATE(XSWU(IIU,IJU,IKU))
  ALLOCATE(XSWD(IIU,IJU,IKU))
  ALLOCATE(XLWU(IIU,IJU,IKU))
  ALLOCATE(XLWD(IIU,IJU,IKU))
  ALLOCATE(XDTHRADSW(IIU,IJU,IKU))
  ALLOCATE(XDTHRADLW(IIU,IJU,IKU))
  ALLOCATE(XRADEFF(IIU,IJU,IKU))
ELSE
  ALLOCATE(XSLOPANG(0,0))
  ALLOCATE(XSLOPAZI(0,0))
  ALLOCATE(XDTHRAD(0,0,0))
  ALLOCATE(XDIRFLASWD(0,0,0))
  ALLOCATE(XDIR_ALB(0,0,0))
  ALLOCATE(XSCA_ALB(0,0,0))
  ALLOCATE(XEMIS  (0,0,0))
  ALLOCATE(XTSRAD (0,0))
  ALLOCATE(XSEA (0,0))
  ALLOCATE(XZS_XY (0,0))
  ALLOCATE(NCLEARCOL_TM1(0,0))
  ALLOCATE(XSWU(0,0,0))
  ALLOCATE(XSWD(0,0,0))
  ALLOCATE(XLWU(0,0,0))
  ALLOCATE(XLWD(0,0,0))
  ALLOCATE(XDTHRADSW(0,0,0))
  ALLOCATE(XDTHRADLW(0,0,0))
  ALLOCATE(XRADEFF(0,0,0))
END IF

IF (CRAD == 'ECMW' .OR. CRAD == 'ECRA') THEN
  ALLOCATE(XSTROATM(31,6))
  ALLOCATE(XSMLSATM(31,6))
  ALLOCATE(XSMLWATM(31,6))
  ALLOCATE(XSPOSATM(31,6))
  ALLOCATE(XSPOWATM(31,6))
  ALLOCATE(XSTATM(31,6))
ELSE
  ALLOCATE(XSTROATM(0,0))
  ALLOCATE(XSMLSATM(0,0))
  ALLOCATE(XSMLWATM(0,0))
  ALLOCATE(XSPOSATM(0,0))
  ALLOCATE(XSPOWATM(0,0))
  ALLOCATE(XSTATM(0,0))
END IF
!
!*       3.8   Module MODD_DEEP_CONVECTION_n
!
IF (CDCONV /= 'NONE' .OR. CSCONV == 'KAFR') THEN
  ALLOCATE(NCOUNTCONV(IIU,IJU))
  ALLOCATE(XDTHCONV(IIU,IJU,IKU))
  ALLOCATE(XDRVCONV(IIU,IJU,IKU))
  ALLOCATE(XDRCCONV(IIU,IJU,IKU))
  ALLOCATE(XDRICONV(IIU,IJU,IKU))
  ALLOCATE(XPRCONV(IIU,IJU))
  ALLOCATE(XPACCONV(IIU,IJU))
  ALLOCATE(XPRSCONV(IIU,IJU))
  ! diagnostics
  IF (LCH_CONV_LINOX) THEN
    ALLOCATE(XIC_RATE(IIU,IJU))
    ALLOCATE(XCG_RATE(IIU,IJU))
    ALLOCATE(XIC_TOTAL_NUMBER(IIU,IJU))
    ALLOCATE(XCG_TOTAL_NUMBER(IIU,IJU))
  ELSE
    ALLOCATE(XIC_RATE(0,0))
    ALLOCATE(XCG_RATE(0,0))
    ALLOCATE(XIC_TOTAL_NUMBER(0,0))
    ALLOCATE(XCG_TOTAL_NUMBER(0,0))
  END IF
  IF ( LDIAGCONV )  THEN
    ALLOCATE(XUMFCONV(IIU,IJU,IKU))
    ALLOCATE(XDMFCONV(IIU,IJU,IKU))
    ALLOCATE(XPRLFLXCONV(IIU,IJU,IKU))
    ALLOCATE(XPRSFLXCONV(IIU,IJU,IKU))
    ALLOCATE(XCAPE(IIU,IJU))
    ALLOCATE(NCLTOPCONV(IIU,IJU))
    ALLOCATE(NCLBASCONV(IIU,IJU))
  ELSE
    ALLOCATE(XUMFCONV(0,0,0))
    ALLOCATE(XDMFCONV(0,0,0))
    ALLOCATE(XPRLFLXCONV(0,0,0))
    ALLOCATE(XPRSFLXCONV(0,0,0))
    ALLOCATE(XCAPE(0,0))
    ALLOCATE(NCLTOPCONV(0,0))
    ALLOCATE(NCLBASCONV(0,0))
  END IF
ELSE
  ALLOCATE(NCOUNTCONV(0,0))
  ALLOCATE(XDTHCONV(0,0,0))
  ALLOCATE(XDRVCONV(0,0,0))
  ALLOCATE(XDRCCONV(0,0,0))
  ALLOCATE(XDRICONV(0,0,0))
  ALLOCATE(XPRCONV(0,0))
  ALLOCATE(XPACCONV(0,0))
  ALLOCATE(XPRSCONV(0,0))
  ALLOCATE(XIC_RATE(0,0))
  ALLOCATE(XCG_RATE(0,0))
  ALLOCATE(XIC_TOTAL_NUMBER(0,0))
  ALLOCATE(XCG_TOTAL_NUMBER(0,0))
  ALLOCATE(XUMFCONV(0,0,0))
  ALLOCATE(XDMFCONV(0,0,0))
  ALLOCATE(XPRLFLXCONV(0,0,0))
  ALLOCATE(XPRSFLXCONV(0,0,0))
  ALLOCATE(XCAPE(0,0))
  ALLOCATE(NCLTOPCONV(0,0))
  ALLOCATE(NCLBASCONV(0,0))
END IF
!
IF ((CDCONV == 'KAFR' .OR. CSCONV == 'KAFR') &
    .AND. LSUBG_COND .AND. LSIG_CONV) THEN
  ALLOCATE(XMFCONV(IIU,IJU,IKU))
ELSE
  ALLOCATE(XMFCONV(0,0,0))
ENDIF
!
IF ((CDCONV == 'KAFR' .OR. CSCONV == 'KAFR') &
    .AND. LCHTRANS .AND. NSV > 0 ) THEN
  ALLOCATE(XDSVCONV(IIU,IJU,IKU,NSV))
ELSE
  ALLOCATE(XDSVCONV(0,0,0,0))
END IF
!
ALLOCATE(XCF_MF(IIU,IJU,IKU)) ; XCF_MF=0.0
ALLOCATE(XRC_MF(IIU,IJU,IKU)) ; XRC_MF=0.0
ALLOCATE(XRI_MF(IIU,IJU,IKU)) ; XRI_MF=0.0
!
!*       3.9   Local variables
!
ALLOCATE(ZJ(IIU,IJU,IKU))
!
!*      3.10 Forcing variables (Module MODD_FRC and MODD_FRCn)
!
IF ( LFORCING ) THEN
  ALLOCATE(XWTFRC(IIU,IJU,IKU)) ; XWTFRC = XUNDEF
  ALLOCATE(XUFRC_PAST(IIU,IJU,IKU)) ; XUFRC_PAST = XUNDEF
  ALLOCATE(XVFRC_PAST(IIU,IJU,IKU)) ; XVFRC_PAST = XUNDEF
ELSE
  ALLOCATE(XWTFRC(0,0,0))
  ALLOCATE(XUFRC_PAST(0,0,0))
  ALLOCATE(XVFRC_PAST(0,0,0))
END IF
!
IF (KMI == 1) THEN
  IF ( LFORCING ) THEN
    ALLOCATE(TDTFRC(NFRC))
    ALLOCATE(XUFRC(IKU,NFRC))
    ALLOCATE(XVFRC(IKU,NFRC))
    ALLOCATE(XWFRC(IKU,NFRC))
    ALLOCATE(XTHFRC(IKU,NFRC))
    ALLOCATE(XRVFRC(IKU,NFRC))
    ALLOCATE(XTENDTHFRC(IKU,NFRC))
    ALLOCATE(XTENDRVFRC(IKU,NFRC))
    ALLOCATE(XGXTHFRC(IKU,NFRC))
    ALLOCATE(XGYTHFRC(IKU,NFRC))
    ALLOCATE(XPGROUNDFRC(NFRC))
    ALLOCATE(XTENDUFRC(IKU,NFRC))
    ALLOCATE(XTENDVFRC(IKU,NFRC))
  ELSE
    ALLOCATE(TDTFRC(0))
    ALLOCATE(XUFRC(0,0))
    ALLOCATE(XVFRC(0,0))
    ALLOCATE(XWFRC(0,0))
    ALLOCATE(XTHFRC(0,0))
    ALLOCATE(XRVFRC(0,0))
    ALLOCATE(XTENDTHFRC(0,0))
    ALLOCATE(XTENDRVFRC(0,0))
    ALLOCATE(XGXTHFRC(0,0))
    ALLOCATE(XGYTHFRC(0,0))
    ALLOCATE(XPGROUNDFRC(0))
    ALLOCATE(XTENDUFRC(0,0))
    ALLOCATE(XTENDVFRC(0,0))
  END IF
ELSE
  !Do not allocate because they are the same on all grids (not 'n' variables)
END IF
! ----------------------------------------------------------------------
!
IF (L2D_ADV_FRC) THEN
  WRITE(ILUOUT,*) 'L2D_ADV_FRC IS SET TO', L2D_ADV_FRC
  WRITE(ILUOUT,*) 'ADV FRC WILL BE SET'
  ALLOCATE(TDTADVFRC(NADVFRC))
  ALLOCATE(XDTHFRC(IIU,IJU,IKU,NADVFRC))  ; XDTHFRC=0.
  ALLOCATE(XDRVFRC(IIU,IJU,IKU,NADVFRC))  ; XDRVFRC=0.
ELSE
  ALLOCATE(TDTADVFRC(0))
  ALLOCATE(XDTHFRC(0,0,0,0))
  ALLOCATE(XDRVFRC(0,0,0,0))
ENDIF

IF (L2D_REL_FRC) THEN
  WRITE(ILUOUT,*) 'L2D_REL_FRC IS SET TO', L2D_REL_FRC
  WRITE(ILUOUT,*) 'REL FRC WILL BE SET'
  ALLOCATE(TDTRELFRC(NRELFRC))
  ALLOCATE(XTHREL(IIU,IJU,IKU,NRELFRC))  ; XTHREL=0.
  ALLOCATE(XRVREL(IIU,IJU,IKU,NRELFRC))  ; XRVREL=0.
ELSE
  ALLOCATE(TDTRELFRC(0))
  ALLOCATE(XTHREL(0,0,0,0))
  ALLOCATE(XRVREL(0,0,0,0))
ENDIF
!
!*      4.11 BIS: Eddy fluxes allocation
!
IF ( LTH_FLX ) THEN
  ALLOCATE(XVTH_FLUX_M(IIU,IJU,IKU)) ; XVTH_FLUX_M = 0.
  ALLOCATE(XWTH_FLUX_M(IIU,IJU,IKU)) ; XWTH_FLUX_M = 0.
  IF (KMI /= 1) THEN
    ALLOCATE(XRTHS_EDDY_FLUX(IIU,IJU,IKU))
    XRTHS_EDDY_FLUX = 0.
  ELSE
    ALLOCATE(XRTHS_EDDY_FLUX(0,0,0))
  ENDIF
ELSE
  ALLOCATE(XVTH_FLUX_M(0,0,0))
  ALLOCATE(XWTH_FLUX_M(0,0,0))
  ALLOCATE(XRTHS_EDDY_FLUX(0,0,0))
END IF
!
IF ( LUV_FLX) THEN
  ALLOCATE(XVU_FLUX_M(IIU,IJU,IKU)) ; XVU_FLUX_M  = 0.
  IF (KMI /= 1) THEN
    ALLOCATE(XRVS_EDDY_FLUX(IIU,IJU,IKU))
    XRVS_EDDY_FLUX = 0.
  ELSE
    ALLOCATE(XRVS_EDDY_FLUX(0,0,0))
  ENDIF
ELSE
  ALLOCATE(XVU_FLUX_M(0,0,0))
  ALLOCATE(XRVS_EDDY_FLUX(0,0,0))
END IF
!
!*      3.11   Module MODD_ICE_CONC_n
!
IF (     (CCLOUD == 'ICE3'.OR.CCLOUD == 'ICE4') .AND.   &
     (CPROGRAM == 'DIAG  '.OR.CPROGRAM == 'MESONH')) THEN
  ALLOCATE(XCIT(IIU,IJU,IKU))
ELSE
  ALLOCATE(XCIT(0,0,0))
END IF
!
IF ( CCLOUD == 'KHKO' .OR. CCLOUD == 'C2R2') THEN
   ALLOCATE(XSUPSAT(IIU,IJU,IKU))
   ALLOCATE(XNACT(IIU,IJU,IKU))
   ALLOCATE(XNPRO(IIU,IJU,IKU))
   ALLOCATE(XSSPRO(IIU,IJU,IKU))
ELSE
   ALLOCATE(XSUPSAT(0,0,0))
   ALLOCATE(XNACT(0,0,0))
   ALLOCATE(XNPRO(0,0,0))
   ALLOCATE(XSSPRO(0,0,0))
END IF
!
!*      3.12   Module MODD_TURB_CLOUD
!
IF (.NOT.(ALLOCATED(XCEI))) ALLOCATE(XCEI(0,0,0))
IF (KMI == NMODEL_CLOUD .AND. CTURBLEN_CLOUD/='NONE' ) THEN
  DEALLOCATE(XCEI)
  ALLOCATE(XCEI(IIU,IJU,IKU))
ENDIF
!
!*      3.13  Module MODD_CH_PH_n
!
IF (LUSECHAQ.AND.(CPROGRAM == 'DIAG  '.OR.CPROGRAM == 'MESONH')) THEN
  IF (LCH_PH) THEN
    ALLOCATE(XPHC(IIU,IJU,IKU))
    IF (NRRL==2) THEN
      ALLOCATE(XPHR(IIU,IJU,IKU))
      ALLOCATE(XACPHR(IIU,IJU))
      XACPHR(:,:) =  0.
    ENDIF
  ENDIF
  IF (NRRL==2) THEN
    ALLOCATE(XACPRAQ(IIU,IJU,NSV_CHAC/2))
    XACPRAQ(:,:,:) = 0.
  ENDIF
ENDIF
IF (.NOT.(ASSOCIATED(XPHC))) ALLOCATE(XPHC(0,0,0))
IF (.NOT.(ASSOCIATED(XPHR))) ALLOCATE(XPHR(0,0,0))
IF (.NOT.(ASSOCIATED(XACPHR))) ALLOCATE(XACPHR(0,0))
IF (.NOT.(ASSOCIATED(XACPRAQ))) ALLOCATE(XACPRAQ(0,0,0))
IF ((LUSECHEM).AND.(CPROGRAM == 'DIAG  ')) THEN
  ALLOCATE(XCHFLX(IIU,IJU,NSV_CHEM))
  XCHFLX(:,:,:) = 0.
ELSE
  ALLOCATE(XCHFLX(0,0,0))
END IF
!
!*          3.14 Module MODD_DRAG
!
IF (LDRAG) THEN
  ALLOCATE(XDRAG(IIU,IJU))
ELSE
  ALLOCATE(XDRAG(0,0))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.    INITIALIZE BUDGET VARIABLES
!              ---------------------------
!
gles = lles_mean .or. lles_resolved  .or. lles_subgrid .or. lles_updraft &
                 .or. lles_downdraft .or. lles_spectra
!Called if budgets are enabled via NAM_BUDGET
!or if LES budgets are enabled via NAM_LES (condition on kmi==1 to call it max once)
if ( ( cbutype /= "NONE" .and. nbumod == kmi ) .or. ( ( gles .or. lcheck ) .and. kmi == 1 ) ) THEN
  call Budget_preallocate()
end if

IF ( CBUTYPE /= "NONE" .AND. NBUMOD == KMI ) THEN
  CALL Ini_budget(ILUOUT,XTSTEP,NSV,NRR,                                      &
             LNUMDIFU,LNUMDIFTH,LNUMDIFSV,                                    &
             LHORELAX_UVWTH,LHORELAX_RV, LHORELAX_RC,LHORELAX_RR,             &
             LHORELAX_RI,LHORELAX_RS,LHORELAX_RG, LHORELAX_RH,LHORELAX_TKE,   &
             LHORELAX_SV, LVE_RELAX, LVE_RELAX_GRD,                           &
             LCHTRANS,LNUDGING,LDRAGTREE,LDEPOTREE,LMAIN_EOL,                 &
             CRAD,CDCONV,CSCONV,CTURB,CTURBDIM,CCLOUD                         )
END IF
!
!-------------------------------------------------------------------------------
!
!
!*       5.    INITIALIZE INTERPOLATION COEFFICIENTS
!
CALL INI_BIKHARDT_n (NDXRATIO_ALL(KMI),NDYRATIO_ALL(KMI),KMI)
!
!-------------------------------------------------------------------------------
!
!*       6.     BUILT THE GENERIC OUTPUT NAME
!               ----------------------------
!
IF (KMI == 1) THEN
  DO IMI = 1 , NMODEL
    WRITE(IO_SURF_MNH_MODEL(IMI)%COUTFILE,'(A,".",I1,".",A)') CEXP,IMI,TRIM(ADJUSTL(CSEG))
    WRITE(YNAME, '(A,".",I1,".",A)') CEXP,IMI,TRIM(ADJUSTL(CSEG))//'.000'
    CALL IO_File_add2list(LUNIT_MODEL(IMI)%TDIAFILE,YNAME,'MNHDIACHRONIC','WRITE', &
                          HDIRNAME=CIO_DIR,                                        &
                          KLFINPRAR=INT(50,KIND=LFIINT),KLFITYPE=1,KLFIVERB=NVERB, &
                          TPDADFILE=LUNIT_MODEL(NDAD(IMI))%TDIAFILE )
  END DO
  !
  TDIAFILE => LUNIT_MODEL(KMI)%TDIAFILE !Necessary because no call to GOTO_MODEL before needing it
  !
  IF (CPROGRAM=='MESONH') THEN
    IF ( NDAD(KMI) == 1)  CDAD_NAME(KMI) = CEXP//'.1.'//CSEG
    IF ( NDAD(KMI) == 2)  CDAD_NAME(KMI) = CEXP//'.2.'//CSEG
    IF ( NDAD(KMI) == 3)  CDAD_NAME(KMI) = CEXP//'.3.'//CSEG
    IF ( NDAD(KMI) == 4)  CDAD_NAME(KMI) = CEXP//'.4.'//CSEG
    IF ( NDAD(KMI) == 5)  CDAD_NAME(KMI) = CEXP//'.5.'//CSEG
    IF ( NDAD(KMI) == 6)  CDAD_NAME(KMI) = CEXP//'.6.'//CSEG
    IF ( NDAD(KMI) == 7)  CDAD_NAME(KMI) = CEXP//'.7.'//CSEG
    IF ( NDAD(KMI) == 8)  CDAD_NAME(KMI) = CEXP//'.8.'//CSEG
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!*       7.    INITIALIZE GRIDS AND METRIC COEFFICIENTS
!              ----------------------------------------
!
CALL SET_GRID(KMI,TPINIFILE,IKU,NIMAX_ll,NJMAX_ll,                       &
              XTSTEP,XSEGLEN,                                            &
              XLONORI,XLATORI,XLON,XLAT,                                 &
              XXHAT,XYHAT,XDXHAT,XDYHAT, XMAP,                           &
              XZS,XZZ,XZHAT,XZTOP,LSLEVE,XLEN1,XLEN2,XZSMT,              &
              ZJ,                                                        &
              TDTMOD,TDTCUR,NSTOP,NBAK_NUMB,NOUT_NUMB,TBACKUPN,TOUTPUTN)
!
CALL METRICS(XMAP,XDXHAT,XDYHAT,XZZ,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
!* update halos of metric coefficients
!
!
CALL UPDATE_METRICS(CLBCX,CLBCY,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
!
CALL SET_DIRCOS(CLBCX,CLBCY,XDXX,XDYY,XDZX,XDZY,TZINITHALO2D_ll,   &
                XDIRCOSXW,XDIRCOSYW,XDIRCOSZW,XCOSSLOPE,XSINSLOPE  )
!
! grid nesting initializations
IF ( KMI == 1 ) THEN
  XTSTEP_MODEL1=XTSTEP
END IF
!
NDT_2_WAY(KMI)=4
!
!-------------------------------------------------------------------------------
!
!*      8.    INITIALIZE DATA FOR JVALUES AND AEROSOLS
!
IF ( LUSECHEM .OR. LCHEMDIAG ) THEN
  IF ((KMI==1).AND.(CPROGRAM == "MESONH".OR.CPROGRAM == "DIAG  "))  &
    CALL CH_INIT_JVALUES(TDTCUR%nday, TDTCUR%nmonth,      &
                         TDTCUR%nyear, ILUOUT, XCH_TUV_DOBNEW)
!
  IF (LORILAM) THEN
    CALL CH_AER_MOD_INIT
  ENDIF
END IF
IF (.NOT.(ASSOCIATED(XMI))) ALLOCATE(XMI(0,0,0,0))
IF (.NOT.(ASSOCIATED(XSOLORG))) ALLOCATE(XSOLORG(0,0,0,0))
!
IF (CCLOUD=='LIMA') CALL INIT_AEROSOL_PROPERTIES
!
!
!
!
!-------------------------------------------------------------------------------
!
!*      9.    FIRE initializations
!              --------------------
!
IF(LBLAZE) THEN
  !
  !     9.1 Array allocation
  !         ----------------
  !
  ! Level Set function
  ALLOCATE(XLSPHI(IIU,IJU,NREFINX*NREFINY));      XLSPHI(:,:,:) = 0.

  ! BMap array
  ! BMap default value
  ! -1 = The fire is not here yet
  ALLOCATE(XBMAP(IIU,IJU,NREFINX*NREFINY));       XBMAP(:,:,:) = -1.

  ! A array
  ALLOCATE(XFMRFA(IIU,IJU,NREFINX*NREFINY));      XFMRFA(:,:,:) = 0.

  ! Wf0 array
  ALLOCATE(XFMWF0(IIU,IJU,NREFINX*NREFINY));      XFMWF0(:,:,:) = 0.

  ! R0 array
  ALLOCATE(XFMR0(IIU,IJU,NREFINX*NREFINY));       XFMR0(:,:,:) = 0.

  ! r00 array
  ALLOCATE(XFMR00(IIU,IJU,NREFINX*NREFINY));      XFMR00(:,:,:) = 0.

  ! Ignition
  ! Default value as 1E6 : Ignition long after simulation end time
  ! 1E6 should be enough as it is more than 11 days
  ALLOCATE(XFMIGNITION(IIU,IJU,NREFINX*NREFINY)); XFMIGNITION(:,:,:) = 1.E6

  ! Fuel type
  ALLOCATE(XFMFUELTYPE(IIU,IJU,NREFINX*NREFINY)); XFMFUELTYPE(:,:,:) = 0.

  ! Residence time function
  ALLOCATE(XFIRETAU(IIU,IJU,NREFINX*NREFINY));    XFIRETAU(:,:,:) = 0.

  ! Rate of spread with wind
  ALLOCATE(XFIRERW(IIU,IJU,NREFINX*NREFINY));    XFIRERW(:,:,:) = 0.

  ! Sensible heat flux parameters
  ! get number of parameters
  SELECT CASE(CHEAT_FLUX_MODEL)
  CASE('CST')
    ! 1 parameter for model : nominal injection value
    INBPARAMSENSIBLE = 1

  CASE('EXP')
    ! 2 parameters for model : Max value and characteristic time
    INBPARAMSENSIBLE = 2

  CASE('EXS')
    ! 3 parameters for model : Max value and characteristic time, smoldering injection value
    INBPARAMSENSIBLE = 3
  END SELECT

  ALLOCATE(XFLUXPARAMH(IIU,IJU,NREFINX*NREFINY,INBPARAMSENSIBLE));
  XFLUXPARAMH(:,:,:,:) = 0.

  ! Latent heat flux parameters
  ! get number of parameters
  SELECT CASE(CLATENT_FLUX_MODEL)
  CASE('CST')
    ! 1 parameter for model : nominal injection value
    INBPARAMLATENT = 1
    
  CASE('EXP')
    ! 2 parameters for model : Max value and characteristic time
    INBPARAMLATENT = 2
  END SELECT

  ALLOCATE(XFLUXPARAMW(IIU,IJU,NREFINX*NREFINY,INBPARAMLATENT));
  XFLUXPARAMW(:,:,:,:) = 0.

  ! Available Sensible energy
  ALLOCATE(XFMASE(IIU,IJU,NREFINX*NREFINY)); XFMASE(:,:,:) = 0.

  ! Available Latent energy
  ALLOCATE(XFMAWC(IIU,IJU,NREFINX*NREFINY)); XFMAWC(:,:,:) = 0.

  ! Walking Ignition map (Arrival time matrix for ignition)
  ALLOCATE(XFMWALKIG(IIU,IJU,NREFINX*NREFINY)); XFMWALKIG(:,:,:) = -1.

  ! Sensible heat flux (W/m2)
  ALLOCATE(XFMFLUXHDH(IIU,IJU,NREFINX*NREFINY)); XFMFLUXHDH(:,:,:) = 0.

  ! Latent heat flux (kg/s/m2)
  ALLOCATE(XFMFLUXHDW(IIU,IJU,NREFINX*NREFINY)); XFMFLUXHDW(:,:,:) = 0.

  ! filtered wind on front normal (m/s)
  ALLOCATE(XFMHWS(IIU,IJU,NREFINX*NREFINY)); XFMHWS(:,:,:) = 0.

  ! filtered wind U (m/s)
  ALLOCATE(XFMWINDU(IIU,IJU,NREFINX*NREFINY)); XFMWINDU(:,:,:) = 0.

  ! filtered wind V (m/s)
  ALLOCATE(XFMWINDV(IIU,IJU,NREFINX*NREFINY)); XFMWINDV(:,:,:) = 0.

  ! filtered wind W (m/s)
  ALLOCATE(XFMWINDW(IIU,IJU,NREFINX*NREFINY)); XFMWINDW(:,:,:) = 0.

  ! Gradient of Level Set on x
  ALLOCATE(XGRADLSPHIX(IIU,IJU,NREFINX*NREFINY)); XGRADLSPHIX(:,:,:) = 0.

  ! Gradient of Level Set on y
  ALLOCATE(XGRADLSPHIY(IIU,IJU,NREFINX*NREFINY)); XGRADLSPHIY(:,:,:) = 0.

  ! Wind for fire
  ALLOCATE(XFIREWIND(IIU,IJU,NREFINX*NREFINY)); XFIREWIND(:,:,:) = 0.

  ! Orographic gradient on fire mesh
  ALLOCATE(XFMGRADOROX(IIU,IJU,NREFINX*NREFINY));      XFMGRADOROX(:,:,:) = 0.
  ALLOCATE(XFMGRADOROY(IIU,IJU,NREFINX*NREFINY));      XFMGRADOROY(:,:,:) = 0.
  !
  !     9.2 Array 2d fire mesh allocation
  !         -----------------------------
  !
  ! Level Set 2d
  ALLOCATE(XLSPHI2D(IIU*NREFINX,IJU*NREFINY)); XLSPHI2D(:,:) = 0.
  ! Gradient of Level Set on x 2d
  ALLOCATE(XGRADLSPHIX2D(IIU*NREFINX,IJU*NREFINY)); XGRADLSPHIX2D(:,:) = 0.

  ! Gradient of Level Set on y 2d
  ALLOCATE(XGRADLSPHIY2D(IIU*NREFINX,IJU*NREFINY)); XGRADLSPHIY2D(:,:) = 0.

  ! Level Set mask on x 2d
  ALLOCATE(XGRADMASKX(IIU*NREFINX,IJU*NREFINY)); XGRADMASKX(:,:) = 0.

  ! Level Set mask on y 2d
  ALLOCATE(XGRADMASKY(IIU*NREFINX,IJU*NREFINY)); XGRADMASKY(:,:) = 0.

  ! burnt surface ratio 2d
  ALLOCATE(XSURFRATIO2D(IIU*NREFINX,IJU*NREFINY)); XSURFRATIO2D(:,:) = 0.

  ! Level Set diffusuon x 2d
  ALLOCATE(XLSDIFFUX2D(IIU*NREFINX,IJU*NREFINY)); XLSDIFFUX2D(:,:) = 0.

  ! Level Set diffusion y 2d
  ALLOCATE(XLSDIFFUY2D(IIU*NREFINX,IJU*NREFINY)); XLSDIFFUY2D(:,:) = 0.

  ! ROS diffusion 2d
  ALLOCATE(XFIRERW2D(IIU*NREFINX,IJU*NREFINY)); XFIRERW2D(:,:) = 0.
  !
  !     9.3 Compute fire mesh size
  !         ----------------------
  !
  XFIREMESHSIZE(1) = (XXHAT(2) - XXHAT(1)) / REAL(NREFINX)
  XFIREMESHSIZE(2) = (XYHAT(2) - XYHAT(1)) / REAL(NREFINY)
  !
ELSE
  !
  !     9.4 Default allocation
  !          ------------------
  !
  ! 3d array
  ALLOCATE(XLSPHI(0,0,0))
  ALLOCATE(XBMAP(0,0,0))
  ALLOCATE(XFMRFA(0,0,0))
  ALLOCATE(XFMR0(0,0,0))
  ALLOCATE(XFMWF0(0,0,0))
  ALLOCATE(XFMR00(0,0,0))
  ALLOCATE(XFMIGNITION(0,0,0))
  ALLOCATE(XFMFUELTYPE(0,0,0))
  ALLOCATE(XFIRETAU(0,0,0))
  ALLOCATE(XFIRERW(0,0,0))
  ALLOCATE(XFLUXPARAMH(0,0,0,0))
  ALLOCATE(XFLUXPARAMW(0,0,0,0))
  ALLOCATE(XFMASE(0,0,0))
  ALLOCATE(XFMAWC(0,0,0))
  ALLOCATE(XFMWALKIG(0,0,0))
  ALLOCATE(XFMFLUXHDH(0,0,0))
  ALLOCATE(XFMFLUXHDW(0,0,0))
  ALLOCATE(XGRADLSPHIX(0,0,0))
  ALLOCATE(XGRADLSPHIY(0,0,0))
  ALLOCATE(XFIREWIND(0,0,0))
  ALLOCATE(XFMGRADOROX(0,0,0))
  ALLOCATE(XFMGRADOROY(0,0,0))
  ! 2d array
  ALLOCATE(XLSPHI2D(0,0))
  ALLOCATE(XGRADLSPHIX2D(0,0))
  ALLOCATE(XGRADLSPHIY2D(0,0))
  ALLOCATE(XGRADMASKX(0,0))
  ALLOCATE(XGRADMASKY(0,0))
  ALLOCATE(XSURFRATIO2D(0,0))
  ALLOCATE(XLSDIFFUX2D(0,0))
  ALLOCATE(XLSDIFFUY2D(0,0))
  ALLOCATE(XFIRERW2D(0,0))
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       9.    INITIALIZE THE PROGNOSTIC FIELDS
!              --------------------------------
!
CALL MPPDB_CHECK3D(XUT,"INI_MODEL_N-before read_field::XUT",PRECISION)
CALL READ_FIELD(KMI,TPINIFILE,IIU,IJU,IKU,                                    &
                CGETTKET,CGETRVT,CGETRCT,CGETRRT,CGETRIT,CGETCIT,CGETZWS,     &
                CGETRST,CGETRGT,CGETRHT,CGETSVT,CGETSRCT,CGETSIGS,CGETCLDFR,  &
                CGETICEFR, CGETBL_DEPTH,CGETSBL_DEPTH,CGETPHC,CGETPHR,        &
                CUVW_ADV_SCHEME, CTEMP_SCHEME,                                &
                NSIZELBX_ll, NSIZELBXU_ll, NSIZELBY_ll, NSIZELBYV_ll,         &
                NSIZELBXTKE_ll,NSIZELBYTKE_ll,                                &
                NSIZELBXR_ll,NSIZELBYR_ll,NSIZELBXSV_ll,NSIZELBYSV_ll,        &
                XUM,XVM,XWM,XDUM,XDVM,XDWM,                                   &
                XUT,XVT,XWT,XTHT,XPABST,XTKET,XRTKEMS,                        &
                XRT,XSVT,XZWS,XCIT,XDRYMASST,XDRYMASSS,                       &
                XSIGS,XSRCT,XCLDFR,XICEFR, XBL_DEPTH,XSBL_DEPTH,XWTHVMF,      &
                XPHC,XPHR, XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,XLSZWSM,           &
                XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,                        &
                XLBXRM,XLBXSVM,                                               &
                XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,                        &
                XLBYRM,XLBYSVM,                                               &
                NFRC,TDTFRC,XUFRC,XVFRC,XWFRC,XTHFRC,XRVFRC,                  &
                XTENDTHFRC,XTENDRVFRC,XGXTHFRC,XGYTHFRC,                      &
                XPGROUNDFRC, XATC,                                            &
                XTENDUFRC, XTENDVFRC,                                         &
                NADVFRC,TDTADVFRC,XDTHFRC,XDRVFRC,                            &
                NRELFRC,TDTRELFRC,XTHREL,XRVREL,                              &
                XVTH_FLUX_M,XWTH_FLUX_M,XVU_FLUX_M,                           &
                XRUS_PRES,XRVS_PRES,XRWS_PRES,XRTHS_CLD,XRRS_CLD,XRSVS_CLD,   &
                ZIBM_LS,XIBM_XMUT,XUMEANW,XVMEANW,XWMEANW,XUMEANN,XVMEANN,    &
                XWMEANN,XUMEANE,XVMEANE,XWMEANE,XUMEANS,XVMEANS,XWMEANS,      &
                XLSPHI, XBMAP, XFMASE, XFMAWC, XFMWINDU, XFMWINDV, XFMWINDW, XFMHWS )

!
!-------------------------------------------------------------------------------
!
!
!*        10.  INITIALIZE REFERENCE STATE
!              ---------------------------
!
!
CALL SET_REF(KMI,TPINIFILE,                         &
             XZZ,XZHAT,ZJ,XDXX,XDYY,CLBCX,CLBCY,    &
             XREFMASS,XMASS_O_PHI0,XLINMASS,        &
             XRHODREF,XTHVREF,XRVREF,XEXNREF,XRHODJ )
!
!-------------------------------------------------------------------------------
!
!*       10.1    INITIALIZE THE TURBULENCE VARIABLES
!               -----------------------------------
!
IF ((CTURB == 'TKEL').AND.(CCONF=='START')) THEN
  CALL MPPDB_CHECK3D(XUT,"INI_MODEL_N-before ini_tke_eps::XUT",PRECISION)
  CALL INI_TKE_EPS(CGETTKET,XTHVREF,XZZ, &
                   XUT,XVT,XTHT,                  &
                   XTKET,TZINITHALO3D_ll    )
  CALL MPPDB_CHECK3D(XUT,"INI_MODEL_N-after ini_tke_eps::XUT",PRECISION)
END IF
!
!
!*       10.2   INITIALIZE THE LES VARIABLES
!               ----------------------------
!
CALL INI_LES_n
!
!-------------------------------------------------------------------------------
!
!*       11.    INITIALIZE THE SOURCE OF TOTAL DRY MASS Md
!               ------------------------------------------
!
IF((KMI==1).AND.LSTEADYLS .AND. (CCONF=='START') ) THEN
   XDRYMASSS = 0.
END IF
!
!-------------------------------------------------------------------------------
!
!*       12.    INITIALIZE THE MICROPHYSICS
!               ----------------------------
!
IF (CELEC == 'NONE') THEN
  CALL INI_MICRO_n(TPINIFILE,ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*       13.    INITIALIZE THE ATMOSPHERIC ELECTRICITY
!               --------------------------------------
!
ELSE
  CALL INI_ELEC_n(ILUOUT, CELEC, CCLOUD, TPINIFILE, &
                  XTSTEP, XZZ,                      &
                  XDXX, XDYY, XDZZ, XDZX, XDZY      )
!
  WRITE (UNIT=ILUOUT,&
  FMT='(/,"ELECTRIC VARIABLES ARE BETWEEN INDEX",I2," AND ",I2)')&
  NSV_ELECBEG, NSV_ELECEND
!
    IF( CGETSVT(NSV_ELECBEG)=='INIT' ) THEN
      XSVT(:,:,:,NSV_ELECBEG) = XCION_POS_FW(:,:,:)                  ! Nb/kg
      XSVT(:,:,:,NSV_ELECEND) = XCION_NEG_FW(:,:,:)
!
      XSVT(:,:,:,NSV_ELECBEG+1:NSV_ELECEND-1) = 0.0
    ELSE  ! Convert elec_variables per m3 into elec_variables per kg of air
      DO JSV = NSV_ELECBEG, NSV_ELECEND
         XSVT(:,:,:,JSV) = XSVT(:,:,:,JSV) / XRHODREF(:,:,:)
      ENDDO
    END IF
END IF
!
!-------------------------------------------------------------------------------
!
!*       14.   INITIALIZE THE LARGE SCALE SOURCES
!              ----------------------------------
!
IF ((KMI==1).AND.(.NOT. LSTEADYLS)) THEN
  CALL MPPDB_CHECK3D(XUT,"INI_MODEL_N-before ini_cpl::XUT",PRECISION)
  CALL INI_CPL(NSTOP,XTSTEP,LSTEADYLS,CCONF,                                  &
               CGETTKET,                                                      &
               CGETRVT,CGETRCT,CGETRRT,CGETRIT,                               &
               CGETRST,CGETRGT,CGETRHT,CGETSVT,LCH_INIT_FIELD,                &
               NSV,NIMAX_ll,NJMAX_ll,                                         &
               NSIZELBX_ll,NSIZELBXU_ll,NSIZELBY_ll,NSIZELBYV_ll,             &
               NSIZELBXTKE_ll,NSIZELBYTKE_ll,                                 &
               NSIZELBXR_ll,NSIZELBYR_ll,NSIZELBXSV_ll,NSIZELBYSV_ll,         &
               XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,XLSZWSM,XDRYMASST,             &
               XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,          &
               XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,          &
               XLSUS,XLSVS,XLSWS,XLSTHS,XLSRVS,XLSZWSS,XDRYMASSS,             &
               XLBXUS,XLBXVS,XLBXWS,XLBXTHS,XLBXTKES,XLBXRS,XLBXSVS,          &
               XLBYUS,XLBYVS,XLBYWS,XLBYTHS,XLBYTKES,XLBYRS,XLBYSVS           )
  CALL MPPDB_CHECK3D(XUT,"INI_MODEL_N-after ini_cpl::XUT",PRECISION)
!
      DO JSV=NSV_CHEMBEG,NSV_CHEMEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_LNOXBEG,NSV_LNOXEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_AERBEG,NSV_AEREND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_DSTBEG,NSV_DSTEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_DSTDEPBEG,NSV_DSTDEPEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_SLTBEG,NSV_SLTEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_SLTDEPBEG,NSV_SLTDEPEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
      DO JSV=NSV_PPBEG,NSV_PPEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
#ifdef MNH_FOREFIRE
      DO JSV=NSV_FFBEG,NSV_FFEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
      !
#endif
! Blaze smoke
DO JSV=NSV_FIREBEG,NSV_FIREEND
  XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
  XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
ENDDO
!
      DO JSV=NSV_CSBEG,NSV_CSEND
        XLBXSVS(:,:,:,JSV)=MAX(XLBXSVS(:,:,:,JSV),0.)
        XLBYSVS(:,:,:,JSV)=MAX(XLBYSVS(:,:,:,JSV),0.)
      ENDDO
!
END IF
!
IF ( KMI > 1) THEN
  ! Use dummy pointers to correct an ifort BUG
  DPTR_XBMX1=>XBMX1
  DPTR_XBMX2=>XBMX2
  DPTR_XBMX3=>XBMX3
  DPTR_XBMX4=>XBMX4
  DPTR_XBMY1=>XBMY1
  DPTR_XBMY2=>XBMY2
  DPTR_XBMY3=>XBMY3
  DPTR_XBMY4=>XBMY4
  DPTR_XBFX1=>XBFX1
  DPTR_XBFX2=>XBFX2
  DPTR_XBFX3=>XBFX3
  DPTR_XBFX4=>XBFX4
  DPTR_XBFY1=>XBFY1
  DPTR_XBFY2=>XBFY2
  DPTR_XBFY3=>XBFY3
  DPTR_XBFY4=>XBFY4
  DPTR_CLBCX=>CLBCX
  DPTR_CLBCY=>CLBCY
  !
  DPTR_XZZ=>XZZ
  DPTR_XZHAT=>XZHAT
  DPTR_XLSUM=>XLSUM
  DPTR_XLSVM=>XLSVM
  DPTR_XLSWM=>XLSWM
  DPTR_XLSTHM=>XLSTHM
  DPTR_XLSRVM=>XLSRVM
  DPTR_XLSZWSM=>XLSZWSM
  DPTR_XLSUS=>XLSUS
  DPTR_XLSVS=>XLSVS
  DPTR_XLSWS=>XLSWS
  DPTR_XLSTHS=>XLSTHS
  DPTR_XLSRVS=>XLSRVS
  DPTR_XLSZWSS=>XLSZWSS
  !
  DPTR_NKLIN_LBXU=>NKLIN_LBXU
  DPTR_XCOEFLIN_LBXU=>XCOEFLIN_LBXU
  DPTR_NKLIN_LBYU=>NKLIN_LBYU
  DPTR_XCOEFLIN_LBYU=>XCOEFLIN_LBYU
  DPTR_NKLIN_LBXV=>NKLIN_LBXV
  DPTR_XCOEFLIN_LBXV=>XCOEFLIN_LBXV
  DPTR_NKLIN_LBYV=>NKLIN_LBYV
  DPTR_XCOEFLIN_LBYV=>XCOEFLIN_LBYV
  DPTR_NKLIN_LBXW=>NKLIN_LBXW
  DPTR_XCOEFLIN_LBXW=>XCOEFLIN_LBXW
  DPTR_NKLIN_LBYW=>NKLIN_LBYW
  DPTR_XCOEFLIN_LBYW=>XCOEFLIN_LBYW
  DPTR_NKLIN_LBXM=>NKLIN_LBXM
  DPTR_XCOEFLIN_LBXM=>XCOEFLIN_LBXM
  DPTR_NKLIN_LBYM=>NKLIN_LBYM
  DPTR_XCOEFLIN_LBYM=>XCOEFLIN_LBYM
  !
  CALL INI_SPAWN_LS_n(NDAD(KMI),XTSTEP,KMI,                                 &
       DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4,DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4,      &
       DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4,DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4,      &
       NDXRATIO_ALL(KMI),NDYRATIO_ALL(KMI),                  &
       DPTR_CLBCX,DPTR_CLBCY,DPTR_XZZ,DPTR_XZHAT,                                &
       LSLEVE,XLEN1,XLEN2,                                   &
       DPTR_XLSUM,DPTR_XLSVM,DPTR_XLSWM,DPTR_XLSTHM,DPTR_XLSRVM,DPTR_XLSZWSM,         &
       DPTR_XLSUS,DPTR_XLSVS,DPTR_XLSWS,DPTR_XLSTHS,DPTR_XLSRVS,DPTR_XLSZWSS,                      &
       DPTR_NKLIN_LBXU,DPTR_XCOEFLIN_LBXU,DPTR_NKLIN_LBYU,DPTR_XCOEFLIN_LBYU,    &
       DPTR_NKLIN_LBXV,DPTR_XCOEFLIN_LBXV,DPTR_NKLIN_LBYV,DPTR_XCOEFLIN_LBYV,    &
       DPTR_NKLIN_LBXW,DPTR_XCOEFLIN_LBXW,DPTR_NKLIN_LBYW,DPTR_XCOEFLIN_LBYW,    &
       DPTR_NKLIN_LBXM,DPTR_XCOEFLIN_LBXM,DPTR_NKLIN_LBYM,DPTR_XCOEFLIN_LBYM     )
  !
  DPTR_XLBXUM=>XLBXUM
  DPTR_XLBYUM=>XLBYUM
  DPTR_XLBXVM=>XLBXVM
  DPTR_XLBYVM=>XLBYVM
  DPTR_XLBXWM=>XLBXWM
  DPTR_XLBYWM=>XLBYWM
  DPTR_XLBXTHM=>XLBXTHM
  DPTR_XLBYTHM=>XLBYTHM
  DPTR_XLBXTKEM=>XLBXTKEM
  DPTR_XLBYTKEM=>XLBYTKEM
  DPTR_XLBXRM=>XLBXRM
  DPTR_XLBYRM=>XLBYRM
  DPTR_XLBXSVM=>XLBXSVM
  DPTR_XLBYSVM=>XLBYSVM
  IF (CCONF=='START')  THEN
  CALL INI_ONE_WAY_n(NDAD(KMI),KMI,                        &
       DPTR_XBMX1,DPTR_XBMX2,DPTR_XBMX3,DPTR_XBMX4,DPTR_XBMY1,DPTR_XBMY2,DPTR_XBMY3,DPTR_XBMY4,        &
       DPTR_XBFX1,DPTR_XBFX2,DPTR_XBFX3,DPTR_XBFX4,DPTR_XBFY1,DPTR_XBFY2,DPTR_XBFY3,DPTR_XBFY4,        &
       NDXRATIO_ALL(KMI),NDYRATIO_ALL(KMI),      &
       DPTR_CLBCX,DPTR_CLBCY,NRIMX,NRIMY,                                &
       DPTR_NKLIN_LBXU,DPTR_XCOEFLIN_LBXU,DPTR_NKLIN_LBYU,DPTR_XCOEFLIN_LBYU,      &
       DPTR_NKLIN_LBXV,DPTR_XCOEFLIN_LBXV,DPTR_NKLIN_LBYV,DPTR_XCOEFLIN_LBYV,      &
       DPTR_NKLIN_LBXW,DPTR_XCOEFLIN_LBXW,DPTR_NKLIN_LBYW,DPTR_XCOEFLIN_LBYW,      &
       DPTR_NKLIN_LBXM,DPTR_XCOEFLIN_LBXM,DPTR_NKLIN_LBYM,DPTR_XCOEFLIN_LBYM,      &
       CCLOUD, LUSECHAQ, LUSECHIC,                                                 &
       DPTR_XLBXUM,DPTR_XLBYUM,DPTR_XLBXVM,DPTR_XLBYVM,DPTR_XLBXWM,DPTR_XLBYWM,    &
       DPTR_XLBXTHM,DPTR_XLBYTHM,                                                  &
       DPTR_XLBXTKEM,DPTR_XLBYTKEM,                                                &
       DPTR_XLBXRM,DPTR_XLBYRM,DPTR_XLBXSVM,DPTR_XLBYSVM                           )
   ENDIF
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       15.    INITIALIZE THE SCALAR VARIABLES
!               -------------------------------
!
IF (LLG .AND. LINIT_LG .AND. CPROGRAM=='MESONH') &
  CALL INI_LG(XXHAT,XYHAT,XZZ,XSVT,XLBXSVM,XLBYSVM)

!
!-------------------------------------------------------------------------------
!
!*       16.    INITIALIZE THE PARAMETERS FOR THE DYNAMICS
!               ------------------------------------------
!
CALL INI_DYNAMICS(XLON,XLAT,XRHODJ,XTHVREF,XMAP,XZZ,XDXHAT,XDYHAT,            &
             XZHAT,CLBCX,CLBCY,XTSTEP,                                        &
             LVE_RELAX,LVE_RELAX_GRD,LHORELAX_UVWTH,LHORELAX_RV,              &
             LHORELAX_RC,LHORELAX_RR,LHORELAX_RI,LHORELAX_RS,LHORELAX_RG,     &
             LHORELAX_RH,LHORELAX_TKE,LHORELAX_SV,                            &
             LHORELAX_SVC2R2,LHORELAX_SVC1R3,LHORELAX_SVELEC,LHORELAX_SVLG,   &
             LHORELAX_SVCHEM,LHORELAX_SVAER,LHORELAX_SVDST,LHORELAX_SVSLT,    &
             LHORELAX_SVPP,LHORELAX_SVCS,LHORELAX_SVCHIC,LHORELAX_SVSNW,      &
#ifdef MNH_FOREFIRE
             LHORELAX_SVFF,                                                   &
#endif
             XRIMKMAX,NRIMX,NRIMY,                                            &
             XALKTOP,XALKGRD,XALZBOT,XALZBAS,                                 &
             XT4DIFU,XT4DIFTH,XT4DIFSV,                                       &
             XCORIOX,XCORIOY,XCORIOZ,XCURVX,XCURVY,                           &
             XDXHATM,XDYHATM,XRHOM,XAF,XBFY,XCF,XTRIGSX,XTRIGSY,NIFAXX,NIFAXY,&
             XALK,XALKW,NALBOT,XALKBAS,XALKWBAS,NALBAS,                       &
             LMASK_RELAX,XKURELAX,XKVRELAX,XKWRELAX,                          &
             XDK2U,XDK4U,XDK2TH,XDK4TH,XDK2SV,XDK4SV,                         &
             LZDIFFU,XZDIFFU_HALO2,                                           &
             XBFB,XBF_SXP2_YP1_Z                                              )
!
!
!*      16.1 Initialize the XDRAG array
!              -------------
IF (LDRAG) THEN
   CALL INI_DRAG(LMOUNT,XZS,XHSTART,NSTART,XDRAG)
ENDIF
!*      16.2 Initialize the LevelSet function
!              -------------
IF (LIBM) THEN
  ALLOCATE(XIBM_LS(IIU,IJU,IKU,4)) ; XIBM_LS  = -XIBM_IEPS
  XIBM_LS(:,:,:,1)=ZIBM_LS(:,:,:)
  DEALLOCATE(ZIBM_LS)
ENDIF
!-------------------------------------------------------------------------------
!
!*      17.    SURFACE FIELDS
!              --------------
!
!*      17.1   Radiative setup
!              ---------------
!
IF (CRAD   /= 'NONE') THEN
  IF (CGETRAD =='INIT') THEN
    GINIRAD  =.TRUE.
  ELSE
    GINIRAD  =.FALSE.
  END IF
  CALL INI_RADIATIONS(TPINIFILE,GINIRAD,TDTCUR,TDTEXP,XZZ, &
                      XDXX, XDYY,                         &
                      XSINDEL,XCOSDEL,XTSIDER,XCORSOL,    &
                      XSLOPANG,XSLOPAZI,                  &
                      XDTHRAD,XDIRFLASWD,XSCAFLASWD,      &
                      XFLALWD,XDIRSRFSWD,NCLEARCOL_TM1,   &
                      XZENITH,XAZIM,                      &
                      TDTRAD_FULL,TDTRAD_CLONLY,          &
                      TZINITHALO2D_ll,                    &
                      XRADEFF,XSWU,XSWD,XLWU,             &
                      XLWD,XDTHRADSW,XDTHRADLW           )
  !
  IF (GINIRAD) CALL SUNPOS_n(XZENITH,PAZIMSOL=XAZIM)
  CALL SURF_SOLAR_GEOM    (XZS, XZS_XY)
  !
  ALLOCATE(XXHAT_ll                 (IIU_ll))
  ALLOCATE(XYHAT_ll                 (IJU_ll))
  ALLOCATE(XZS_ll                   (IIU_ll,IJU_ll))
  ALLOCATE(XZS_XY_ll                (IIU_ll,IJU_ll))
  !
  CALL GATHERALL_FIELD_ll('XY',XZS,XZS_ll,IRESP)
  CALL GATHERALL_FIELD_ll('XY',XZS_XY,XZS_XY_ll,IRESP)
  CALL GATHERALL_FIELD_ll('XX',XXHAT,XXHAT_ll,IRESP)
  CALL GATHERALL_FIELD_ll('YY',XYHAT,XYHAT_ll,IRESP)
  XZS_MAX_ll=MAXVAL(XZS_ll)
ELSE
  XAZIM       = XPI
  XZENITH     = XPI/2.
  XDIRSRFSWD  = 0.
  XSCAFLASWD  = 0.
  XFLALWD     = 300.  ! W/m2
  XTSIDER     = 0.
END IF
!
!
CALL INI_SW_SETUP (CRAD,NSWB_MNH,XSW_BANDS)
CALL INI_LW_SETUP (CRAD,NLWB_MNH,XLW_BANDS)
!
!
!       17.1.1 Special initialisation for CO2 content
!              CO2 (molar mass=44) horizontally and vertically homogeneous at 360 ppm
!
XCCO2 = 360.0E-06 * 44.0E-03 / XMD
#ifdef MNH_ECRAD
RCCO2 = 360.0E-06 * 44.0E-03 / XMD
#endif
!
!
!*      17.2   Externalized surface fields
!              ---------------------------
!
ALLOCATE(ZCO2(IIU,IJU))
ZCO2(:,:) = XCCO2
!

ALLOCATE(ZDIR_ALB(IIU,IJU,NSWB_MNH))
ALLOCATE(ZSCA_ALB(IIU,IJU,NSWB_MNH))
ALLOCATE(ZEMIS  (IIU,IJU,NLWB_MNH))
ALLOCATE(ZTSRAD (IIU,IJU))
!
IF (LCOUPLES.AND.(KMI>1))THEN
  CSURF ="NONE"
ELSE
  IF ((TPINIFILE%NMNHVERSION(1)==4 .AND. TPINIFILE%NMNHVERSION(2)>=6) .OR. TPINIFILE%NMNHVERSION(1)>4) THEN
    CALL IO_Field_read(TPINIFILE,'SURF',CSURF)
  ELSE
    CSURF = "EXTE"
  END IF
END IF
!
!
IF (CSURF=='EXTE' .AND. (CPROGRAM=='MESONH' .OR. CPROGRAM=='DIAG  ')) THEN
  ! ouverture du fichier PGD
  IF  ( LEN_TRIM(CINIFILEPGD) > 0 ) THEN
    CALL IO_File_add2list(TINIFILEPGD,TRIM(CINIFILEPGD),'PGD','READ',KLFITYPE=2,KLFIVERB=NVERB)
    CALL IO_File_open(TINIFILEPGD,KRESP=IRESP)
    LUNIT_MODEL(KMI)%TINIFILEPGD => TINIFILEPGD
    IF (IRESP/=0) THEN
      WRITE(ILUOUT,FMT=*) "INI_MODEL_n ERROR TO OPEN THE FILE CINIFILEPGD=",CINIFILEPGD
      WRITE(ILUOUT,FMT=*) "CHECK YOUR NAMELIST NAM_LUNITn"
    !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_MODEL_n','')
    ENDIF
  ELSE
  ! case after a spawning
    CINIFILEPGD = TPINIFILE%CNAME
  END IF
  !
  CALL GOTO_SURFEX(KMI)
#ifdef CPLOASIS
  CALL SFX_OASIS_READ_NAM(CPROGRAM,XTSTEP)
  WRITE(*,*) 'SFX-OASIS: READ NAM_SFX_SEA_CPL OK'
#endif
  !* initialization of surface
  CALL INIT_GROUND_PARAM_n ('ALL',SIZE(CSV),CSV,ZCO2,                             &
                            XZENITH,XAZIM,XSW_BANDS,XLW_BANDS,ZDIR_ALB,ZSCA_ALB,  &
                            ZEMIS,ZTSRAD                                )
  !
  IF (SIZE(XEMIS)>0) THEN
    XDIR_ALB = ZDIR_ALB
    XSCA_ALB = ZSCA_ALB
    XEMIS    = ZEMIS
    XTSRAD   = ZTSRAD
    CALL MNHGET_SURF_PARAM_n (PSEA=XSEA)
  END IF
ELSE
  !* fields not physically necessary, but must be initialized
  IF (SIZE(XEMIS)>0) THEN
    XDIR_ALB = 0.
    XSCA_ALB = 0.
    XEMIS    = 1.
    XTSRAD   = XTT
    XSEA     = 1.
  END IF
END IF
IF (CSURF=='EXTE' .AND. (CPROGRAM=='SPAWN ')) THEN
  ! ouverture du fichier PGD
  CALL IO_File_add2list(TINIFILEPGD,TRIM(CINIFILEPGD),'PGD','READ',KLFITYPE=2,KLFIVERB=NVERB)
  CALL IO_File_open(TINIFILEPGD,KRESP=IRESP)
  LUNIT_MODEL(KMI)%TINIFILEPGD => TINIFILEPGD
  IF (IRESP/=0) THEN
    WRITE(ILUOUT,FMT=*) "INI_MODEL_n ERROR TO OPEN THE FILE CINIFILEPGD=",CINIFILEPGD
    WRITE(ILUOUT,FMT=*) "CHECK YOUR NAMELIST NAM_LUNIT2_SPA"
    !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_MODEL_n','')
  ENDIF
ENDIF
!
IF (.NOT.ASSOCIATED(TINIFILEPGD)) TINIFILEPGD => TFILE_DUMMY
!
  !* special case after spawning in prep_real_case
IF (CSURF=='EXRM' .AND. CPROGRAM=='REAL  ') CSURF = 'EXTE'
!
DEALLOCATE(ZDIR_ALB)
DEALLOCATE(ZSCA_ALB)
DEALLOCATE(ZEMIS   )
DEALLOCATE(ZTSRAD  )
!
DEALLOCATE(ZCO2)
!
!
!* in a RESTART case, reads surface radiative quantities in the MESONH file
!
IF ((CRAD  == 'ECMW' .OR. CRAD  == 'ECRA') .AND. CGETRAD=='READ') THEN
  CALL INI_SURF_RAD(TPINIFILE, XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD)
END IF
!
!
!*      17.3   Mesonh fields
!              -------------
!
IF (CPROGRAM/='REAL  ') CALL MNHREAD_ZS_DUMMY_n(TINIFILEPGD)
!
!-------------------------------------------------------------------------------
!
!*       18.    INITIALIZE THE PARAMETERS FOR THE PHYSICS
!               -----------------------------------------
!
IF (CRAD   == 'ECMW') THEN
!
!* get cover mask for aerosols
!
  IF (CPROGRAM=='MESONH' .OR. CPROGRAM=='DIAG  ') THEN
    ALLOCATE(ZSEA(IIU,IJU))
    ALLOCATE(ZTOWN(IIU,IJU))
    ALLOCATE(ZBARE(IIU,IJU))
    IF (CSURF=='EXTE') THEN
      CALL GOTO_SURFEX(KMI)
      CALL MNHGET_SURF_PARAM_n(PSEA=ZSEA,PTOWN=ZTOWN,PBARE=ZBARE)
    ELSE
      ZSEA (:,:) = 1.
      ZTOWN(:,:) = 0.
      ZBARE(:,:) = 0.
    END IF
!
    IF ( CAOP=='EXPL' .AND. LDUST .AND. KMI==1) THEN
      ALLOCATE( XEXT_COEFF_WVL_LKT_DUST( NMAX_RADIUS_LKT_DUST, NMAX_SIGMA_LKT_DUST, NMAX_WVL_SW_DUST ) )
      ALLOCATE( XEXT_COEFF_550_LKT_DUST( NMAX_RADIUS_LKT_DUST, NMAX_SIGMA_LKT_DUST )                   )
      ALLOCATE( XPIZA_LKT_DUST         ( NMAX_RADIUS_LKT_DUST, NMAX_SIGMA_LKT_DUST, NMAX_WVL_SW_DUST ) )
      ALLOCATE( XCGA_LKT_DUST          ( NMAX_RADIUS_LKT_DUST, NMAX_SIGMA_LKT_DUST, NMAX_WVL_SW_DUST ) )
    END IF
!
    IF ( CAOP=='EXPL' .AND. LSALT .AND. KMI==1) THEN
      ALLOCATE( XEXT_COEFF_WVL_LKT_SALT( NMAX_RADIUS_LKT_SALT, NMAX_SIGMA_LKT_SALT, NMAX_WVL_SW_SALT ) )
      ALLOCATE( XEXT_COEFF_550_LKT_SALT( NMAX_RADIUS_LKT_SALT, NMAX_SIGMA_LKT_SALT )                   )
      ALLOCATE( XPIZA_LKT_SALT         ( NMAX_RADIUS_LKT_SALT, NMAX_SIGMA_LKT_SALT, NMAX_WVL_SW_SALT ) )
      ALLOCATE( XCGA_LKT_SALT          ( NMAX_RADIUS_LKT_SALT, NMAX_SIGMA_LKT_SALT, NMAX_WVL_SW_SALT ) )
    END IF
!
    CALL INI_RADIATIONS_ECMWF (XZHAT,XPABST,XTHT,XTSRAD,XLAT,XLON,TDTCUR,TDTEXP,       &
                               CLW,NDLON,NFLEV,NFLUX,NRAD,NSWB_OLD,CAER,NAER,NSTATM,   &
                               XSTATM,ZSEA,ZTOWN,ZBARE,XOZON, XAER,XDST_WL, LSUBG_COND )
!
    DEALLOCATE(ZSEA,ZTOWN,ZBARE)
    ALLOCATE (XAER_CLIM(SIZE(XAER,1),SIZE(XAER,2),SIZE(XAER,3),SIZE(XAER,4)))
    XAER_CLIM(:,:,:,:) =XAER(:,:,:,:)
!
  END IF

ELSE IF (CRAD   == 'ECRA') THEN
#ifdef MNH_ECRAD
!* get cover mask for aerosols
!
  IF (CPROGRAM=='MESONH' .OR. CPROGRAM=='DIAG  ') THEN
    ALLOCATE(ZSEA(IIU,IJU))
    ALLOCATE(ZTOWN(IIU,IJU))
    ALLOCATE(ZBARE(IIU,IJU))
    IF (CSURF=='EXTE') THEN
      CALL GOTO_SURFEX(KMI)
      CALL MNHGET_SURF_PARAM_n(PSEA=ZSEA,PTOWN=ZTOWN,PBARE=ZBARE)
    ELSE
      ZSEA (:,:) = 1.
      ZTOWN(:,:) = 0.
      ZBARE(:,:) = 0.
    END IF
!
    CALL INI_RADIATIONS_ECRAD (XZHAT,XPABST,XTHT,XTSRAD,XLAT,XLON,TDTCUR,TDTEXP,       &
                               CLW,NDLON,NFLEV,NFLUX,NRAD,NSWB_OLD,CAER,NAER,NSTATM,   &
                               XSTATM,ZSEA,ZTOWN,ZBARE,XOZON, XAER,XDST_WL, LSUBG_COND )

    DEALLOCATE(ZSEA,ZTOWN,ZBARE)
    ALLOCATE (XAER_CLIM(SIZE(XAER,1),SIZE(XAER,2),SIZE(XAER,3),SIZE(XAER,4)))
    XAER_CLIM(:,:,:,:) = XAER(:,:,:,:)
!
  END IF
#endif
ELSE
  ALLOCATE (XOZON(0,0,0))
  ALLOCATE (XAER(0,0,0,0))
  ALLOCATE (XDST_WL(0,0,0,0))
  ALLOCATE (XAER_CLIM(0,0,0,0))
END IF
!
!
!
IF (CDCONV /= 'NONE' .OR. CSCONV == 'KAFR') THEN
  IF (CGETCONV=='INIT') THEN
    GINIDCONV=.TRUE.
  ELSE
    GINIDCONV=.FALSE.
  END IF
!
!  commensurability between convection calling time and time step
!
  XDTCONV=XTSTEP*REAL( INT( (MIN(XDTCONV,1800.)+1.E-10)/XTSTEP ) )
  XDTCONV=MAX( XDTCONV, XTSTEP )
  IF (NVERB>=10) THEN
    WRITE(ILUOUT,*) 'XDTCONV has been set to : ',XDTCONV
  END IF
  CALL INI_DEEP_CONVECTION (TPINIFILE,GINIDCONV,TDTCUR,                      &
                           NCOUNTCONV,XDTHCONV,XDRVCONV,XDRCCONV,            &
                           XDRICONV,XPRCONV,XPRSCONV,XPACCONV,               &
                           XUMFCONV,XDMFCONV,XMFCONV,XPRLFLXCONV,XPRSFLXCONV,&
                           XCAPE,NCLTOPCONV,NCLBASCONV,                      &
                           TDTDCONV, CGETSVCONV, XDSVCONV,                   &
                           LCH_CONV_LINOX, XIC_RATE, XCG_RATE,               &
                           XIC_TOTAL_NUMBER, XCG_TOTAL_NUMBER                )

END IF
!
!-------------------------------------------------------------------------------
!
!
!*      19.    ALLOCATION OF THE TEMPORAL SERIES
!              ---------------------------------
!
IF (LSERIES .AND. CPROGRAM/='DIAG  ') CALL INI_SERIES_n
!
!-------------------------------------------------------------------------------
!
!
!*      20.   (re)initialize scalar variables
!             -------------------------------
!
!
IF ( LUSECHEM .OR. LCHEMDIAG ) THEN
  IF (CPROGRAM=='MESONH'.AND.CCONF=='RESTA') LCH_INIT_FIELD =.FALSE.
  IF (CPROGRAM=='MESONH'.OR. CPROGRAM=='DIAG  ' .OR. CPROGRAM=='IDEAL ') &
        CALL CH_INIT_FIELD_n(KMI, ILUOUT, NVERB)
END IF
!
!-------------------------------------------------------------------------------
!
!*      21.    UPDATE HALO
!              -----------
!
!
CALL UPDATE_HALO_ll(TZINITHALO3D_ll,IINFO_ll)
CALL UPDATE_HALO_ll(TZINITHALO2D_ll,IINFO_ll)
CALL CLEANLIST_ll(TZINITHALO3D_ll)
CALL CLEANLIST_ll(TZINITHALO2D_ll)
!
!
!-------------------------------------------------------------------------------
!
!*      22.    DEALLOCATION
!              -------------
!
DEALLOCATE(ZJ)
!
DEALLOCATE(XSTROATM)
DEALLOCATE(XSMLSATM)
DEALLOCATE(XSMLWATM)
DEALLOCATE(XSPOSATM)
DEALLOCATE(XSPOWATM)
!
!-------------------------------------------------------------------------------
!
!*      23.     BALLOON and AIRCRAFT initializations
!              ------------------------------------
!
CALL INI_AIRCRAFT_BALLOON(TPINIFILE,XTSTEP, TDTSEG, XSEGLEN, NRR, NSV, &
                          IKU,CTURB=="TKEL" ,                          &
                          XLATORI, XLONORI                             )
!
!-------------------------------------------------------------------------------
!
!*      24.     STATION initializations
!              -----------------------
!
CALL INI_SURFSTATION_n(XTSTEP, XSEGLEN, NRR, NSV, &
                       CTURB=="TKEL" , KMI,       &
                       XLATORI, XLONORI           )
!
!-------------------------------------------------------------------------------
!
!*      25.     PROFILER initializations
!              ------------------------
!
CALL INI_POSPROFILER_n(XTSTEP, XSEGLEN, NRR, NSV,  &
                       CTURB=="TKEL",              &
                       XLATORI, XLONORI            )
!
!-------------------------------------------------------------------------------
!
!*      26.     Prognostic aerosols
!              ------------------------
!
IF ( ( CRAD=='ECMW' .OR. CRAD=='ECRA' ) .AND. CAOP=='EXPL' .AND. LORILAM ) THEN
  ALLOCATE(POLYTAU(6,10,8,6,13))
  ALLOCATE(POLYSSA(6,10,8,6,13))
  ALLOCATE(POLYG  (6,10,8,6,13))
  CALL INI_AEROSET1
  CALL INI_AEROSET2
  CALL INI_AEROSET3
  CALL INI_AEROSET4
  CALL INI_AEROSET5
  CALL INI_AEROSET6
END IF
#ifdef MNH_FOREFIRE
!
!-------------------------------------------------------------------------------
!
!*      27.    FOREFIRE initializations
!              ------------------------
!

! Coupling with ForeFire if resolution is low enough
!---------------------------------------------------
IF ( LFOREFIRE .AND. 0.5*(XXHAT(2)-XXHAT(1)+XYHAT(2)-XYHAT(1)) < COUPLINGRES ) THEN
  FFCOUPLING = .TRUE.
ELSE
  FFCOUPLING = .FALSE.
ENDIF

! Initializing the ForeFire variables
!------------------------------------
IF ( LFOREFIRE ) THEN
  CALL INIT_FOREFIRE_n(KMI, ILUOUT, IP &
        , TDTCUR%nyear, TDTCUR%nmonth, TDTCUR%nday, TDTCUR%xtime, XTSTEP)
END IF
#endif

!-------------------------------------------------------------------------------
!
!*      30.   Total production/Loss for chemical species
!
IF (LCHEMDIAG)  THEN
  CALL CH_INIT_PRODLOSSTOT_n(ILUOUT)
  IF (NEQ_PLT>0) THEN
    ALLOCATE(XPROD(IIU,IJU,IKU,NEQ_PLT))
    ALLOCATE(XLOSS(IIU,IJU,IKU,NEQ_PLT))
    XPROD=0.0
    XLOSS=0.0
  ELSE
    ALLOCATE(XPROD(0,0,0,0))
    ALLOCATE(XLOSS(0,0,0,0))
  END IF
ELSE
  ALLOCATE(XPROD(0,0,0,0))
  ALLOCATE(XLOSS(0,0,0,0))
END IF
!
!-------------------------------------------------------------------------------
!
!*     31. Extended production/loss terms for chemical species
!
IF (LCHEMDIAG) THEN
  CALL CH_INIT_BUDGET_n(ILUOUT)
  IF (NEQ_BUDGET>0) THEN
    ALLOCATE(IINDEX(2,NNONZEROTERMS))
    ALLOCATE(IIND(NEQ_BUDGET))
    CALL CH_NONZEROTERMS(KMI,IINDEX,NNONZEROTERMS)
    ALLOCATE(XTCHEM(NEQ_BUDGET))
    DO JM=1,NEQ_BUDGET
      IIND(JM)=COUNT((IINDEX(1,:))==NSPEC_BUDGET(JM))
      ALLOCATE(XTCHEM(JM)%NB_REAC(IIND(JM)))
      ALLOCATE(XTCHEM(JM)%XB_REAC(IIU,IJU,IKU,IIND(JM)))
    END DO
    DEALLOCATE(IIND)
    DEALLOCATE(IINDEX)
  ELSE
    ALLOCATE(XTCHEM(0))
  END IF
ELSE
  ALLOCATE(XTCHEM(0))
END IF
!-------------------------------------------------------------------------------
!
!*     32. Wind turbine
!
IF (LMAIN_EOL .AND. KMI == NMODEL_EOL) THEN
 ALLOCATE(XFX_RG(IIU,IJU,IKU))
 ALLOCATE(XFY_RG(IIU,IJU,IKU))
 ALLOCATE(XFZ_RG(IIU,IJU,IKU))
 ALLOCATE(XFX_SMR_RG(IIU,IJU,IKU))
 ALLOCATE(XFY_SMR_RG(IIU,IJU,IKU))
 ALLOCATE(XFZ_SMR_RG(IIU,IJU,IKU))
 SELECT CASE(CMETH_EOL)
  CASE('ADNR')
   CALL INI_EOL_ADNR
  CASE('ALM')
   CALL INI_EOL_ALM(XDXX,XDYY)
 END SELECT
END IF
!
!*     33.  Auto-coupling Atmos-Ocean LES NH
!
IF (LCOUPLES) THEN
 ALLOCATE(XSSUFL_C(IIU,IJU,1)); XSSUFL_C=0.0
 ALLOCATE(XSSVFL_C(IIU,IJU,1)); XSSVFL_C=0.0
 ALLOCATE(XSSTFL_C(IIU,IJU,1)); XSSTFL_C=0.0
 ALLOCATE(XSSRFL_C(IIU,IJU,1)); XSSRFL_C=0.
ELSE
 ALLOCATE(XSSUFL_C(0,0,0))
 ALLOCATE(XSSVFL_C(0,0,0))
 ALLOCATE(XSSTFL_C(0,0,0))
 ALLOCATE(XSSRFL_C(0,0,0))
END IF
!
END SUBROUTINE INI_MODEL_n
