!MNH_LIC Copyright 1994-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################### 
      MODULE MODI_READ_EXSEG_n
!     ######################
!
INTERFACE
!
      SUBROUTINE READ_EXSEG_n(KMI,TPEXSEGFILE,HCONF,OFLAT,OUSERV,                  &
                   OUSERC,OUSERR,OUSERI,OUSECI,OUSERS,OUSERG,OUSERH,               &
                   OUSECHEM,OUSECHAQ,OUSECHIC,OCH_PH,OCH_CONV_LINOX,OSALT,         &
                   ODEPOS_SLT, ODUST,ODEPOS_DST, OCHTRANS,                         &
                   OORILAM,ODEPOS_AER, OLG,OPASPOL, OFIRE,                         &
#ifdef MNH_FOREFIRE
                   OFOREFIRE,                                                      &
#endif
                   OLNOX_EXPLICIT,                                                 &
                   OCONDSAMP,OBLOWSNOW,                                            &
                   KRIMX,KRIMY, KSV_USER,                                          &
                   HTURB,HTOM,ORMC01,HRAD,HDCONV,HSCONV,HCLOUD,HELEC,              &
                   HEQNSYS,PTSTEP_ALL,HINIFILEPGD                                  )
!
USE MODD_IO,   ONLY: TFILEDATA
!
INTEGER,            INTENT(IN) :: KMI    ! Model index
TYPE(TFILEDATA),    INTENT(IN) :: TPEXSEGFILE ! EXSEG file
!     The following variables are read by READ_DESFM in DESFM descriptor : 
CHARACTER (LEN=*),  INTENT(IN) :: HCONF  ! configuration var. linked to FMfile
LOGICAL,            INTENT(IN) :: OFLAT  ! Logical for zero orography
LOGICAL,            INTENT(IN) :: OUSERV,OUSERC,OUSERR,OUSERI,OUSERS, &
                                  OUSERG,OUSERH  ! kind of moist variables in 
                                                 ! FMfile
LOGICAL,            INTENT(IN) :: OUSECI         ! ice concentration in
                                                 ! FMfile
LOGICAL,            INTENT(IN) :: OUSECHEM       ! Chemical FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OUSECHAQ       ! Aqueous chemical FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OUSECHIC       ! Ice chemical FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OCH_PH         ! pH FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OCH_CONV_LINOX ! LiNOx FLAG in FMFILE
LOGICAL,            INTENT(IN) :: ODUST          ! Dust FLAG in FMFILE
LOGICAL,DIMENSION(:), INTENT(IN)  :: ODEPOS_DST  ! Dust wet deposition FLAG in FMFILE     
LOGICAL,DIMENSION(:), INTENT(IN)  :: ODEPOS_SLT  ! Sea Salt wet deposition FLAG in FMFILE     
LOGICAL,DIMENSION(:), INTENT(IN)  :: ODEPOS_AER  ! Orilam wet deposition FLAG in FMFILE     
LOGICAL,            INTENT(IN) :: OSALT          ! Sea Salt FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OORILAM        ! Orilam FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OPASPOL        ! Passive pollutant FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OFIRE          ! Blaze FLAG in FMFILE
#ifdef MNH_FOREFIRE
LOGICAL,            INTENT(IN) :: OFOREFIRE      ! ForeFire FLAG in FMFILE
#endif
LOGICAL,            INTENT(IN) :: OLNOX_EXPLICIT ! explicit LNOx FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OCONDSAMP      ! Conditional sampling FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OBLOWSNOW     ! Blowing snow FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OCHTRANS       ! LCHTRANS FLAG in FMFILE

LOGICAL,            INTENT(IN) :: OLG            ! lagrangian FLAG in FMFILE
INTEGER,            INTENT(IN) :: KRIMX, KRIMY   ! number of points for the
                       ! horizontal relaxation for the outermost verticals
INTEGER,            INTENT(IN) :: KSV_USER    ! number of additional scalar
                                         ! variables in FMfile 
CHARACTER (LEN=*),  INTENT(IN) :: HTURB  ! Kind of turbulence parameterization
                                         ! used to  produce FMFILE
CHARACTER (LEN=*),  INTENT(IN) :: HTOM   ! Kind of third order moment
LOGICAL,            INTENT(IN) :: ORMC01 ! flag for RMC01 SBL computations
CHARACTER (LEN=*),  INTENT(IN) :: HRAD   ! Kind of radiation scheme
CHARACTER (LEN=4),  INTENT(IN) :: HDCONV ! Kind of deep convection scheme
CHARACTER (LEN=4),  INTENT(IN) :: HSCONV ! Kind of shallow convection scheme
CHARACTER (LEN=4),  INTENT(IN) :: HCLOUD ! Kind of microphysical scheme
CHARACTER (LEN=4),  INTENT(IN) :: HELEC  ! Kind of electrical scheme
CHARACTER (LEN=*),  INTENT(IN) :: HEQNSYS! type of equations' system
REAL,DIMENSION(:),  INTENT(INOUT):: PTSTEP_ALL ! Time STEP of ALL models
CHARACTER (LEN=*),  INTENT(IN) :: HINIFILEPGD ! name of PGD file
!
END SUBROUTINE READ_EXSEG_n
!
END INTERFACE
!
END MODULE MODI_READ_EXSEG_n
!
!
!     #########################################################################
      SUBROUTINE READ_EXSEG_n(KMI,TPEXSEGFILE,HCONF,OFLAT,OUSERV,                  &
                   OUSERC,OUSERR,OUSERI,OUSECI,OUSERS,OUSERG,OUSERH,               &
                   OUSECHEM,OUSECHAQ,OUSECHIC,OCH_PH,OCH_CONV_LINOX,OSALT,         &
                   ODEPOS_SLT, ODUST,ODEPOS_DST, OCHTRANS,                         &
                   OORILAM,ODEPOS_AER, OLG,OPASPOL, OFIRE,                         &
#ifdef MNH_FOREFIRE
                   OFOREFIRE,                                                      &
#endif
                   OLNOX_EXPLICIT,                                                 &
                   OCONDSAMP, OBLOWSNOW,                                           &
                   KRIMX,KRIMY, KSV_USER,                                          &
                   HTURB,HTOM,ORMC01,HRAD,HDCONV,HSCONV,HCLOUD,HELEC,              &
                   HEQNSYS,PTSTEP_ALL,HINIFILEPGD                                  )
!     #########################################################################
!
!!****  *READ_EXSEG_n * - routine to read  the descriptor file EXSEG
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to read the descriptor file called 
!     EXSEG and to control the coherence with FMfile data . 
!       
!!
!!**  METHOD
!!    ------
!!      The descriptor file is read. Namelists (NAMXXXn) which contain
!!    variables linked to one nested model are at the beginning of the file.
!!    Namelists (NAMXXX) which contain variables common to all models
!!    are at the end of the file. When the  model index is different from 1, 
!!    the end of the file (namelists NAMXXX) is not read. 
!!
!!      Coherence between the initial file (description read in DESFM file)
!!    and the segment to perform (description read in EXSEG file)
!!    is checked for segment achievement  configurations 
!!    or  postprocessing configuration. The get indicators are set according 
!!      to the following check :
!!
!!      - segment achievement and preinit configurations :
!!
!!       * if there is no turbulence kinetic energy in  initial
!!      file (HTURB='NONE'), and the segment to perform requires a turbulence
!!      parameterization (CTURB /= 'NONE'), the get indicators for turbulence 
!!      kinetic energy variables are set to 'INIT'; i.e. these variables will be
!!      set equal to  zero by READ_FIELD according to the get indicators.
!!       * The same procedure is applied to the dissipation of TKE. 
!!       * if there is no moist variables RRn  in  initial file (OUSERn=.FALSE.)
!!      and the segment to perform requires  moist variables RRn
!!      (LUSERn=.TRUE.),  the get indicators for moist variables RRn are set
!!      equal to 'INIT'; i.e. these variables will be set equal to  zero by 
!!      READ_FIELD according to the get indicators.
!!       * if there are KSV_USER additional scalar variables in initial file and the
!!      segment to perform needs more than KSV_USER additional variables,  the get
!!      indicators for these (NSV_USER-KSV_USER) additional scalar variables are set 
!!      equal to 'INIT'; i.e. these variables will be set equal to  zero by 
!!      READ_FIELD according to the get indicators. If the segment to perform
!!      needs less additional scalar variables than there are in initial file,
!!      the get indicators for these (KSV_USER - NSV_USER) additional scalar variables are
!!      set equal to 'SKIP'.      
!!       * warning messages are printed if the fields in initial file are the
!!      same at time t and t-dt (HCONF='START') and  a leap-frog advance 
!!      at first time step will be used  for the segment to perform
!!      (CCONF='RESTA'); It is likewise when HCONF='RESTA' and CCONF='START'.
!!       * A warning message is printed if the orography in initial file is zero
!!      (OFLAT=.TRUE.) and the segment to perform considers no-zero orography
!!      (LFLAT=.FALSE.). It is likewise for LFLAT=.TRUE. and OFLAT=.FALSE..
!!       If the segment to perform requires zero orography (LFLAT=.TRUE.), the
!!       orography (XZS) will not read in initial file but  set equal to zero 
!!       by SET_GRID.  
!!       * check of the depths of the Lateral Damping Layer in x and y 
!!       direction is performed
!!       * If some coupling files are specified, LSTEADYLS is set to T
!!       * If no coupling files are specified, LSTEADYLS is set to F
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!     Module MODN_CONF : CCONF,LTHINSHELL,LFLAT,NMODEL,NVERB     
!!
!!      Module MODN_DYN : LCORIO, LZDIFFU
!!
!!      Module MODN_NESTING : NDAD(m),NDTRATIO(m),XWAY(m)
!!
!!      Module MODN_BUDGET : CBUTYPE,XBULEN
!!
!!      Module MODN_CONF1 : LUSERV,LUSERC,LUSERR,LUSERI,LUSERS,LUSERG,LUSERH,CSEG
!!
!!      Module MODN_DYN1 : XTSTEP,CPRESOPT,NITR,XRELAX
!!
!!      Module MODD_ADV1 : CMET_ADV_SCHEME,CSV_ADV_SCHEME,CUVW_ADV_SCHEME,NLITER
!!
!!      Module MODN_PARAM1 : CTURB,CRAD,CDCONV,CSCONV
!!
!!      Module MODN_LUNIT1 : 
!!      Module MODN_LBC1 : CLBCX,CLBCY,NLBLX,NLBLY,XCPHASE,XPOND
!!
!!      Module MODN_TURB_n : CTURBLEN,CTURBDIM
!!
!!      Module MODD_GET1: 
!!        CGETTKEM,CGETTKET,
!!        CGETRVM,CGETRCM,CGETRRM,CGETRIM,CGETRSM,CGETRGM,CGETRHM
!!        CGETRVT,CGETRCT,CGETRRT,CGETRIT,CGETRST,CGETRGT,CGETRHT,CGETSVM
!!        CGETSVT,CGETSIGS,CGETSRCM,CGETSRCT 
!!        NCPL_NBR,NCPL_TIMES,NCPL_CUR
!!      Module MODN_LES  : contains declaration of the control parameters
!!                                for Large Eddy Simulations' storages
!!                                for the forcing
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the  documentation (routine READ_EXSEG_n)
!!      
!!
!!    AUTHOR
!!    ------
!!  	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/06/94 
!!      Modification   26/10/94  (Stein)  remove NAM_GET from the Namelists
!!      present in DESFM + change the namelist names
!!      Modification   22/11/94  (Stein)  add GET indicator for phi
!!      Modification   21/12/94  (Stein)  add GET indicator for LS fields  
!!      Modification   06/01/95  (Stein)  bug in the test for Scalar Var.
!!      Modifications  09/01/95  (Stein)  add the turbulence scheme
!!      Modifications  09/01/95  (Stein)  add the 1D switch
!!      Modifications  10/03/95  (Mallet) add coherence in coupling case
!!      Modifications  16/03/95  (Stein)  remove R from the historical variables
!!      Modifications  01/03/95  (Hereil) add the budget namelists
!!      Modifications  16/06/95  (Stein)  coherence control for the
!!                microphysical scheme + remove the wrong messge for RESTA conf
!!      Modifications  30/06/95  (Stein)  conditionnal reading of the fields
!!                used by the moist turbulence scheme
!!      Modifications  12/09/95  (Pinty)  add the radiation scheme
!!      Modification   06/02/96  (J.Vila) implement scalar advection schemes
!!      Modifications  24/02/96  (Stein)  change the default value for CCPLFILE
!!      Modifications  02/05/96  (Stein Jabouille) change the Z0SEA activation
!!      Modifications  24/05/96  (Stein)  change the SRC SIGS control
!!      Modifications  08/09/96  (Masson) the coupling file names are reset to
!!                               default value " " before reading in EXSEG1.nam
!!                               to avoid extra non-existant coupling files
!!
!!      Modifications  25/04/95  (K.Suhre)add namelist NAM_BLANK
!!                                        add read for LFORCING
!!                     25/04/95  (K.Suhre)add namelist NAM_FRC
!!                                        and switch checking
!!                     06/08/96  (K.Suhre)add namelist NAM_CH_MNHCn
!!                                        and NAM_CH_SOLVER
!!      Modifications  10/10/96  (Stein)  change SRC into SRCM and SRCT
!!      Modifications  11/04/96  (Pinty)  add the rain-ice microphysical scheme
!!      Modifications  11/01/97  (Pinty)  add the deep convection scheme
!!      Modifications  22/05/97  (Lafore) gridnesting implementation
!!      Modifications  22/06/97  (Stein)  add the absolute pressure + cleaning
!!      Modifications  25/08/97  (Masson) add tests on surface schemes
!!                     22/10/97  (Stein) remove the RIMX /= 0 control
!!                                       + new namelist + cleaning
!!      Modifications  17/04/98  (Masson) add tests on character variables
!!      Modification   15/03/99  (Masson) add tests on PROGRAM
!!      Modification   04/01/00  (Masson) removes TSZ0 case
!!      Modification   04/06/00  (Pinty)  add C2R2 scheme
!!                     11/12/00  (Tomasini) add CSEA_FLUX to MODD_PARAMn
!!                                          delete the test on SST_FRC only in 1D
!!      Modification   22/01/01  (Gazen)  change NSV,KSV to NSV_USER,KSV_USER and add
!!                                        NSV_* variables initialization
!!      Modification   15/10/01  (Mallet) allow namelists in different orders
!!      Modification   18/03/02  (Solmon)  new radiation scheme test
!!      Modification   29/11/02 (JP Pinty) add C3R5, ICE2, ICE4, ELEC
!!      Modification   06/11/02  (Masson)  new LES BL height diagnostic
!!      Modification   06/11/02  (Jabouille)  remove LTHINSHELL LFORCING test
!!      Modification   01/12/03  (Gazen)   change Chemical scheme interface
!!      Modification   01/2004   (Masson) removes surface (externalization)
!!      Modification   01/2005   (Masson) removes 1D and 2D switches
!!      Modification   04/2005   (Tulet)  add dust, orilam
!!      Modification   03/2006   (O.Geoffroy) Add KHKO scheme
!!      Modification   04/2006   (Maric)  include 4th order advection scheme
!!      Modification   05/2006   (Masson) add nudging
!!      Modification   05/2006   Remove KEPS
!!      Modification   04/2006   (Maric)  include PPM advection scheme
!!      Modification   04/2006   (J.Escobar) Bug dollarn add CALL UPDATE_NAM_CONFN  
!!      Modifications  01/2007   (Malardel,Pergaud)  add the MF shallow
!!                               convection scheme MODN_PARAM_MFSHALL_n
!!      Modification   09/2009   (J.Escobar) add more info on relaxation problems
!!      Modification   09/2011   (J.Escobar) re-add 'ZRESI' choose
!!      Modification   12/2011   (C.Lac) Adaptation to FIT temporal scheme 
!!      Modification   12/2012   (S.Bielli) add NAM_NCOUT for netcdf output (removed 08/07/2016)
!!      Modification   02/2012   (Pialat/Tulet) add ForeFire
!!      Modification   02/2012   (T.Lunet) add of new Runge-Kutta methods
!!      Modification   01/2015   (C. Barthe) add explicit LNOx
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      M.Leriche 18/12/2015 : bug chimie glace dans prep_real_case
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!      Modification   02/2016   (M.Leriche) treat gas and aq. chemicals separately
!!      P.Wautelet 08/07/2016 : removed MNH_NCWRIT define
!!      Modification   10/2016    (C.LAC) Add OSPLIT_WENO + Add droplet
!!                                deposition + Add max values
!!      Modification   11/2016   (Ph. Wautelet) Allocate/initialise some output/backup structures
!!      Modification   03/2017   (JP Chaboureau) Fix the initialization of
!!                                               LUSERx-type variables for LIMA
!!      M.Leriche      06/2017 for spawn and prep_real avoid abort if wet dep for 
!!                             aerosol and no cloud scheme defined
!!      Q.Libois       02/2018  ECRAD
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Modification   07/2017   (V. Vionnet) add blowing snow scheme
!!      Modification   01/2019   (Q. Rodier) define XCEDIS depending on BL89 or RM17 mixing length
!!      Modification   01/2019   (P. Wautelet) bugs correction: incorrect writes
!!      Modification   01/2019   (R. Honnert) remove SURF in CMF_UPDRAFT
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  C. Lac         11/2019: correction in the drag formula and application to building in addition to tree
!  Q. Rodier      03/2020: add abort if use of any LHORELAX and cyclic conditions
!  F.Auguste      02/2021: add IBM
!  T.Nagel        02/2021: add turbulence recycling
!  E.Jezequel     02/2021: add stations read from CSV file
!  P. Wautelet 09/03/2021: simplify allocation of scalar variable names
!  P. Wautelet 09/03/2021: move some chemistry initializations to ini_nsv
!  P. Wautelet 10/03/2021: move scalar variable name initializations to ini_nsv
!  R. Honnert  23/04/2021: add HM21 mixing length and delete HRIO and BOUT from CMF_UPDRAFT
!  S. Riette   11/05/2021  HighLow cloud
!  A. Costes      12/2021: add Blaze fire model
!  R. Schoetter   12/2021:  multi-level coupling between MesoNH and SURFEX
!  P. Wautelet 27/04/2022: add namelist for profilers
!  P. Wautelet 24/06/2022: remove check on CSTORAGE_TYPE for restart of ForeFire variables
!  P. Wautelet 13/07/2022: add namelist for flyers and balloons
!  P. Wautelet 19/08/2022: add namelist for aircrafts
!------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_AIRCRAFT_BALLOON, ONLY: NAIRCRAFTS, NBALLOONS
USE MODD_BLOWSNOW
USE MODD_BUDGET
USE MODD_CH_AEROSOL
USE MODD_CH_M9_n, ONLY : NEQ
USE MODD_CONDSAMP
USE MODD_CONF
USE MODD_CONF_n,  ONLY: CSTORAGE_TYPE
USE MODD_CONFZ
! USE MODD_DRAG_n
USE MODD_DUST
USE MODD_DYN
USE MODD_DYN_n, ONLY : LHORELAX_SVLIMA, LHORELAX_SVFIRE
#ifdef MNH_FOREFIRE
USE MODD_FOREFIRE
#endif
USE MODD_GET_n
USE MODD_GR_FIELD_n
USE MODD_IO,   ONLY: TFILEDATA
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_NSV,NSV_USER_n=>NSV_USER
USE MODD_PARAMETERS
USE MODD_PASPOL
USE MODD_SALT
USE MODD_VAR_ll,  ONLY: NPROC
USE MODD_VISCOSITY

USE MODE_MSG
USE MODE_POS

USE MODI_INI_NSV
USE MODI_TEST_NAM_VAR

USE MODN_2D_FRC
USE MODN_ADV_n      ! The final filling of these modules for the model n is
USE MODN_AIRCRAFTS, ONLY: AIRCRAFTS_NML_ALLOCATE, NAM_AIRCRAFTS
USE MODN_BACKUP
USE MODN_BALLOONS,  ONLY: BALLOONS_NML_ALLOCATE, NAM_BALLOONS
USE MODN_BLANK_n
USE MODN_BLOWSNOW
USE MODN_BLOWSNOW_n
USE MODN_BUDGET
USE MODN_CH_MNHC_n
USE MODN_CH_ORILAM
USE MODN_CH_SOLVER_n
USE MODN_CONDSAMP
USE MODN_CONF
USE MODN_CONF_n
USE MODN_CONFZ
USE MODN_DRAGBLDG_n
USE MODN_COUPLING_LEVELS_n
USE MODN_DRAG_n
USE MODN_DRAGTREE_n
USE MODN_DUST
USE MODN_DYN
USE MODN_DYN_n      ! to avoid the duplication of this routine for each model.
USE MODN_ELEC
USE MODN_EOL
USE MODN_EOL_ADNR
USE MODN_EOL_ALM
USE MODN_FIRE_n
USE MODN_FLYERS
#ifdef MNH_FOREFIRE
USE MODN_FOREFIRE
#endif
USE MODN_FRC
USE MODN_IBM_PARAM_n
USE MODN_LATZ_EDFLX
USE MODN_LBC_n      ! routine is used for each nested model. This has been done
USE MODN_LES
USE MODN_LUNIT_n
USE MODN_MEAN
USE MODN_NESTING
USE MODN_NUDGING_n
USE MODN_OUTPUT
USE MODN_PARAM_C1R3, ONLY : NAM_PARAM_C1R3, CPRISTINE_ICE_C1R3,    &
                            CHEVRIMED_ICE_C1R3
USE MODN_PARAM_C2R2, ONLY : EPARAM_CCN=>HPARAM_CCN, EINI_CCN=>HINI_CCN, &
                            WNUC=>XNUC, WALPHAC=>XALPHAC, NAM_PARAM_C2R2
USE MODN_PARAM_ECRAD_n
USE MODD_PARAM_ICE_n, ONLY : PARAM_ICEN_INIT, PARAM_ICEN, CSUBG_AUCV_RC, CSUBG_AUCV_RI
USE MODN_PARAM_KAFR_n
USE MODD_PARAM_LIMA, ONLY : FINI_CCN=>HINI_CCN,PARAM_LIMA_INIT,NMOD_CCN,LSCAV, &
                            CPRISTINE_ICE_LIMA, CHEVRIMED_ICE_LIMA, NMOD_IFN, NMOD_IMM, &
                            LACTI, LNUCL, XALPHAC, XNUC, LMEYERS, &
                            LPTSPLIT, LSPRO, LADJ, LKHKO, &
                            NMOM_C, NMOM_R, NMOM_I, NMOM_S, NMOM_G, NMOM_H
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALLN_INIT
USE MODN_PARAM_n    ! realized in subroutine ini_model n
USE MODN_PARAM_RAD_n
USE MODN_PASPOL
USE MODN_PROFILER_n, LDIAG_SURFRAD_PROF => LDIAG_SURFRAD
USE MODN_RECYCL_PARAM_n
USE MODN_SALT
USE MODN_SERIES
USE MODN_SERIES_n
USE MODN_STATION_n,  LDIAG_SURFRAD_STAT => LDIAG_SURFRAD
USE MODD_TURB_n, ONLY: TURBN_INIT, CTOM, CTURBDIM, LRMC01, LHARAT, &
                       LCLOUDMODIFLM, CTURBLEN_CLOUD, XCEI_MIN, XCEI_MAX
USE MODD_NEB_n, ONLY: NEBN_INIT, LSIGMAS, LSUBG_COND, CCONDENS, LSTATNW
USE MODN_VISCOSITY
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
!
INTEGER,            INTENT(IN) :: KMI    ! Model index
TYPE(TFILEDATA),    INTENT(IN) :: TPEXSEGFILE ! EXSEG file
!     The following variables are read by READ_DESFM in DESFM descriptor : 
CHARACTER (LEN=*),  INTENT(IN) :: HCONF  ! configuration var. linked to FMfile
LOGICAL,            INTENT(IN) :: OFLAT  ! Logical for zero orography
LOGICAL,            INTENT(IN) :: OUSERV,OUSERC,OUSERR,OUSERI,OUSERS, &
                                  OUSERG,OUSERH  ! kind of moist variables in 
                                                 ! FMfile
LOGICAL,            INTENT(IN) :: OUSECI         ! ice concentration in
                                                 ! FMfile
LOGICAL,            INTENT(IN) :: OUSECHEM       ! Chemical FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OUSECHAQ       ! Aqueous chemical FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OUSECHIC       ! Ice chemical FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OCH_PH         ! pH FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OCH_CONV_LINOX ! LiNOx FLAG in FMFILE
LOGICAL,            INTENT(IN) :: ODUST          ! Dust FLAG in FMFILE
LOGICAL,DIMENSION(:), INTENT(IN) :: ODEPOS_DST   ! Dust Deposition FLAG in FMFILE
LOGICAL,DIMENSION(:), INTENT(IN) :: ODEPOS_SLT   ! Sea Salt wet deposition FLAG in FMFILE     
LOGICAL,DIMENSION(:), INTENT(IN) :: ODEPOS_AER   ! Orilam wet deposition FLAG in FMFILE     
LOGICAL,            INTENT(IN) :: OSALT          ! Sea Salt FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OORILAM        ! Orilam FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OPASPOL        ! Passive pollutant FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OFIRE          ! Blaze FLAG in FMFILE
#ifdef MNH_FOREFIRE
LOGICAL,            INTENT(IN) :: OFOREFIRE      ! ForeFire FLAG in FMFILE
#endif
LOGICAL,            INTENT(IN) :: OLNOX_EXPLICIT ! explicit LNOx FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OCONDSAMP      ! Conditional sampling FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OCHTRANS       ! LCHTRANS FLAG in FMFILE
LOGICAL,            INTENT(IN) :: OBLOWSNOW     ! Blowing snow FLAG in FMFILE

LOGICAL,            INTENT(IN) :: OLG            ! lagrangian FLAG in FMFILE
INTEGER,            INTENT(IN) :: KRIMX, KRIMY   ! number of points for the
                       ! horizontal relaxation for the outermost verticals
INTEGER,            INTENT(IN) :: KSV_USER    ! number of additional scalar
                                         ! variables in FMfile 
CHARACTER (LEN=*),  INTENT(IN) :: HTURB  ! Kind of turbulence parameterization
                                         ! used to  produce FMFILE
CHARACTER (LEN=*),  INTENT(IN) :: HTOM   ! Kind of third order moment
LOGICAL,            INTENT(IN) :: ORMC01 ! flag for RMC01 SBL computations
CHARACTER (LEN=*),  INTENT(IN) :: HRAD   ! Kind of radiation scheme
CHARACTER (LEN=4),  INTENT(IN) :: HDCONV ! Kind of deep convection scheme
CHARACTER (LEN=4),  INTENT(IN) :: HSCONV ! Kind of shallow convection scheme
CHARACTER (LEN=4),  INTENT(IN) :: HCLOUD ! Kind of microphysical scheme
CHARACTER (LEN=4),  INTENT(IN) :: HELEC  ! Kind of electrical scheme
CHARACTER (LEN=*),  INTENT(IN) :: HEQNSYS! type of equations' system
REAL,DIMENSION(:),  INTENT(INOUT):: PTSTEP_ALL ! Time STEP of ALL models
CHARACTER (LEN=*),  INTENT(IN) :: HINIFILEPGD ! name of PGD file
!
!*       0.2   declarations of local variables
!
CHARACTER(LEN=3) :: YMODEL
INTEGER :: ILUSEG,ILUOUT ! logical unit numbers of EXSEG file and outputlisting
INTEGER :: JS,JCI,JI,JSV       ! Loop indexes 
LOGICAL :: GRELAX              
LOGICAL :: GFOUND              ! Return code when searching namelist
!
!-------------------------------------------------------------------------------
!
!*       1.    READ EXSEG FILE
!              ---------------
!
CALL PRINT_MSG(NVERB_DEBUG,'IO','READ_EXSEG_n','called for '//TRIM(TPEXSEGFILE%CNAME))
!
ILUSEG = TPEXSEGFILE%NLU
ILUOUT = TLUOUT%NLU
!
CALL INIT_NAM_LUNITN
CCPLFILE(:)="                            "
CALL INIT_NAM_CONFN
CALL INIT_NAM_DYNN
CALL INIT_NAM_ADVN
CALL INIT_NAM_DRAGTREEN
CALL INIT_NAM_DRAGBLDGN
CALL INIT_NAM_COUPLING_LEVELSN
CALL INIT_NAM_PARAMN
CALL INIT_NAM_PARAM_RADN
#ifdef MNH_ECRAD
CALL INIT_NAM_PARAM_ECRADN
#endif
CALL INIT_NAM_PARAM_KAFRN
CALL INIT_NAM_LBCN
CALL INIT_NAM_NUDGINGN
CALL INIT_NAM_BLANKN
CALL INIT_NAM_DRAGN
CALL INIT_NAM_IBM_PARAMN
CALL INIT_NAM_RECYCL_PARAMN
CALL INIT_NAM_CH_MNHCN
CALL INIT_NAM_CH_SOLVERN
CALL INIT_NAM_SERIESN
CALL INIT_NAM_BLOWSNOWN
CALL INIT_NAM_PROFILERn
CALL INIT_NAM_STATIONn
CALL INIT_NAM_FIREn
!
WRITE(UNIT=ILUOUT,FMT="(/,'READING THE EXSEG.NAM FILE')")
CALL POSNAM( TPEXSEGFILE, 'NAM_LUNITN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_LUNITn)
CALL POSNAM( TPEXSEGFILE, 'NAM_CONFN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_CONFn)
CALL POSNAM( TPEXSEGFILE, 'NAM_DYNN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_DYNn)
CALL POSNAM( TPEXSEGFILE, 'NAM_ADVN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_ADVn)
CALL POSNAM( TPEXSEGFILE, 'NAM_PARAMN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_PARAMn)
CALL POSNAM( TPEXSEGFILE, 'NAM_PARAM_RADN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_PARAM_RADn)
#ifdef MNH_ECRAD
CALL POSNAM( TPEXSEGFILE, 'NAM_PARAM_ECRADN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_PARAM_ECRADn)
#endif
CALL POSNAM( TPEXSEGFILE, 'NAM_PARAM_KAFRN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_PARAM_KAFRn)
CALL PARAM_MFSHALLN_INIT(CPROGRAM, TPEXSEGFILE, .FALSE., ILUOUT, .FALSE., .TRUE., .FALSE., 0)
CALL POSNAM( TPEXSEGFILE, 'NAM_LBCN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_LBCn)
CALL POSNAM( TPEXSEGFILE, 'NAM_NUDGINGN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_NUDGINGn)
CALL TURBN_INIT(CPROGRAM, TPEXSEGFILE, .FALSE., ILUOUT, .FALSE., .TRUE., .FALSE., 0)
CALL NEBN_INIT(CPROGRAM, TPEXSEGFILE, .FALSE., ILUOUT, .FALSE., .TRUE., .FALSE., 0)
CALL PARAM_ICEN_INIT(CPROGRAM, TPEXSEGFILE, .FALSE., ILUOUT, .FALSE., .TRUE., .FALSE., 0)
CALL POSNAM( TPEXSEGFILE, 'NAM_DRAGN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_DRAGn)
CALL POSNAM( TPEXSEGFILE, 'NAM_IBM_PARAMN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_IBM_PARAMn)
CALL POSNAM( TPEXSEGFILE, 'NAM_RECYCL_PARAMN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_RECYCL_PARAMn)
CALL POSNAM( TPEXSEGFILE, 'NAM_CH_MNHCN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_CH_MNHCn)
CALL POSNAM( TPEXSEGFILE, 'NAM_CH_SOLVERN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_CH_SOLVERn)
CALL POSNAM( TPEXSEGFILE, 'NAM_SERIESN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_SERIESn)
CALL POSNAM( TPEXSEGFILE, 'NAM_BLANKN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_BLANKn)
CALL POSNAM( TPEXSEGFILE, 'NAM_BLOWSNOWN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_BLOWSNOWn)
CALL POSNAM( TPEXSEGFILE, 'NAM_DRAGTREEN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_DRAGTREEn)
CALL POSNAM( TPEXSEGFILE, 'NAM_DRAGBLDGN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_DRAGBLDGn)
CALL POSNAM( TPEXSEGFILE,'NAM_COUPLING_LEVELSN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_COUPLING_LEVELSn)
CALL POSNAM( TPEXSEGFILE, 'NAM_EOL', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_EOL)
CALL POSNAM( TPEXSEGFILE, 'NAM_EOL_ADNR', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_EOL_ADNR)
CALL POSNAM( TPEXSEGFILE, 'NAM_EOL_ALM', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_EOL_ALM)
CALL POSNAM( TPEXSEGFILE, 'NAM_PROFILERN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_PROFILERn)
CALL POSNAM( TPEXSEGFILE, 'NAM_STATIONN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_STATIONn)
CALL POSNAM( TPEXSEGFILE, 'NAM_FIREN', GFOUND )
IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_FIREn)
!
IF (KMI == 1) THEN                                               
  WRITE(UNIT=ILUOUT,FMT="(' namelists common to all the models ')")
  CALL POSNAM( TPEXSEGFILE, 'NAM_CONF', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_CONF)
  CALL POSNAM( TPEXSEGFILE, 'NAM_CONFZ', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_CONFZ)
  CALL POSNAM( TPEXSEGFILE, 'NAM_DYN', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_DYN)
  CALL POSNAM( TPEXSEGFILE, 'NAM_NESTING', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_NESTING)
  CALL POSNAM( TPEXSEGFILE, 'NAM_BACKUP', GFOUND )
  IF (GFOUND) THEN
    !Should have been allocated before in READ_DESFM_n
    IF (.NOT.ALLOCATED(XBAK_TIME)) THEN
      ALLOCATE(XBAK_TIME(NMODEL,JPOUTMAX))
      XBAK_TIME(:,:) = XNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(XOUT_TIME)) THEN
      ALLOCATE(XOUT_TIME(NMODEL,JPOUTMAX)) !Allocate *OUT* variables to prevent
      XOUT_TIME(:,:) = XNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(NBAK_STEP)) THEN
      ALLOCATE(NBAK_STEP(NMODEL,JPOUTMAX))
      NBAK_STEP(:,:) = NNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(NOUT_STEP)) THEN
      ALLOCATE(NOUT_STEP(NMODEL,JPOUTMAX)) !problems if NAM_OUTPUT does not exist
      NOUT_STEP(:,:) = NNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(COUT_VAR)) THEN
      ALLOCATE(COUT_VAR (NMODEL,JPOUTVARMAX))
      COUT_VAR(:,:)  = ''
    END IF
    READ(UNIT=ILUSEG,NML=NAM_BACKUP)
  ELSE
    CALL POSNAM( TPEXSEGFILE, 'NAM_FMOUT', GFOUND )
    IF (GFOUND) THEN
      CALL PRINT_MSG(NVERB_FATAL,'IO','READ_EXSEG_n','use namelist NAM_BACKUP instead of namelist NAM_FMOUT')
    ELSE
      IF (CPROGRAM=='MESONH') CALL PRINT_MSG(NVERB_ERROR,'IO','READ_EXSEG_n','namelist NAM_BACKUP not found')
    END IF
  END IF
  CALL POSNAM( TPEXSEGFILE, 'NAM_OUTPUT', GFOUND )
  IF (GFOUND) THEN
    !Should have been allocated before in READ_DESFM_n
    IF (.NOT.ALLOCATED(XBAK_TIME)) THEN
      ALLOCATE(XBAK_TIME(NMODEL,JPOUTMAX)) !Allocate *BAK* variables to prevent
      XBAK_TIME(:,:) = XNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(XOUT_TIME)) THEN
      ALLOCATE(XOUT_TIME(NMODEL,JPOUTMAX))
      XOUT_TIME(:,:) = XNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(NBAK_STEP)) THEN
      ALLOCATE(NBAK_STEP(NMODEL,JPOUTMAX)) !problems if NAM_BACKUP does not exist
      NBAK_STEP(:,:) = NNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(NOUT_STEP)) THEN
      ALLOCATE(NOUT_STEP(NMODEL,JPOUTMAX))
      NOUT_STEP(:,:) = NNEGUNDEF
    END IF
    IF (.NOT.ALLOCATED(COUT_VAR)) THEN
      ALLOCATE(COUT_VAR (NMODEL,JPOUTVARMAX))
      COUT_VAR(:,:)  = ''
    END IF
    READ(UNIT=ILUSEG,NML=NAM_OUTPUT)
  END IF
  CALL POSNAM( TPEXSEGFILE, 'NAM_BUDGET', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_BUDGET)

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RU', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RU ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RU was already allocated' )
      DEALLOCATE( CBULIST_RU )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RU(NBULISTMAXLINES) )
    CBULIST_RU(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RU)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RU(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RV', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RV ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RV was already allocated' )
      DEALLOCATE( CBULIST_RV )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RV(NBULISTMAXLINES) )
    CBULIST_RV(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RV)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RV(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RW', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RW ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RW was already allocated' )
      DEALLOCATE( CBULIST_RW )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RW(NBULISTMAXLINES) )
    CBULIST_RW(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RW)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RW(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RTH', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RTH ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RTH was already allocated' )
      DEALLOCATE( CBULIST_RTH )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RTH(NBULISTMAXLINES) )
    CBULIST_RTH(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RTH)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RTH(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RTKE', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RTKE ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RTKE was already allocated' )
      DEALLOCATE( CBULIST_RTKE )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RTKE(NBULISTMAXLINES) )
    CBULIST_RTKE(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RTKE)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RTKE(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RRV', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RRV ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RRV was already allocated' )
      DEALLOCATE( CBULIST_RRV )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRV(NBULISTMAXLINES) )
    CBULIST_RRV(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RRV)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRV(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RRC', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RRC ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RRC was already allocated' )
      DEALLOCATE( CBULIST_RRC )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRC(NBULISTMAXLINES) )
    CBULIST_RRC(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RRC)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRC(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RRR', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RRR ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RRR was already allocated' )
      DEALLOCATE( CBULIST_RRR )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRR(NBULISTMAXLINES) )
    CBULIST_RRR(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RRR)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRR(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RRI', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RRI ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RRI was already allocated' )
      DEALLOCATE( CBULIST_RRI )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRI(NBULISTMAXLINES) )
    CBULIST_RRI(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RRI)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRI(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RRS', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RRS ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RRS was already allocated' )
      DEALLOCATE( CBULIST_RRS )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRS(NBULISTMAXLINES) )
    CBULIST_RRS(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RRS)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRS(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RRG', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RRG ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RRG was already allocated' )
      DEALLOCATE( CBULIST_RRG )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRG(NBULISTMAXLINES) )
    CBULIST_RRG(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RRG)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRG(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RRH', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RRH ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RRH was already allocated' )
      DEALLOCATE( CBULIST_RRH )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRH(NBULISTMAXLINES) )
    CBULIST_RRH(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RRH)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RRH(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_BU_RSV', GFOUND )
  IF (GFOUND) THEN
    IF ( ALLOCATED( CBULIST_RSV ) ) THEN
      CALL Print_msg( NVERB_WARNING, 'IO', 'READ_EXSEG_n', 'unexpected: CBULIST_RSV was already allocated' )
      DEALLOCATE( CBULIST_RSV )
    END IF
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RSV(NBULISTMAXLINES) )
    CBULIST_RSV(:) = ''
    READ(UNIT=ILUSEG,NML=NAM_BU_RSV)
  ELSE
    ALLOCATE( CHARACTER(LEN=NBULISTMAXLEN) :: CBULIST_RSV(0) )
  END IF

  CALL POSNAM( TPEXSEGFILE, 'NAM_LES', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_LES)
  CALL POSNAM( TPEXSEGFILE, 'NAM_MEAN', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_MEAN)
  CALL POSNAM( TPEXSEGFILE, 'NAM_PDF', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_PDF)
  CALL POSNAM( TPEXSEGFILE, 'NAM_FRC', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_FRC)
  CALL POSNAM( TPEXSEGFILE, 'NAM_PARAM_C2R2', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_PARAM_C2R2)
  CALL POSNAM( TPEXSEGFILE, 'NAM_PARAM_C1R3', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_PARAM_C1R3)
  CALL PARAM_LIMA_INIT(CPROGRAM, TPEXSEGFILE, .FALSE., ILUOUT, .FALSE., .TRUE., .FALSE., 0)
  CALL POSNAM( TPEXSEGFILE, 'NAM_ELEC', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_ELEC)
  CALL POSNAM( TPEXSEGFILE, 'NAM_SERIES', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_SERIES)
  CALL POSNAM( TPEXSEGFILE, 'NAM_CH_ORILAM', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_CH_ORILAM)
  CALL POSNAM( TPEXSEGFILE, 'NAM_DUST', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_DUST)
  CALL POSNAM( TPEXSEGFILE, 'NAM_SALT', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_SALT)
  CALL POSNAM( TPEXSEGFILE, 'NAM_PASPOL', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_PASPOL)
#ifdef MNH_FOREFIRE
  CALL POSNAM( TPEXSEGFILE, 'NAM_FOREFIRE', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_FOREFIRE)
#endif
  CALL POSNAM( TPEXSEGFILE, 'NAM_CONDSAMP', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_CONDSAMP)
  CALL POSNAM( TPEXSEGFILE, 'NAM_2D_FRC', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_2D_FRC)
  CALL POSNAM( TPEXSEGFILE, 'NAM_LATZ_EDFLX', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_LATZ_EDFLX)
  CALL POSNAM( TPEXSEGFILE, 'NAM_BLOWSNOW', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_BLOWSNOW)
  CALL POSNAM( TPEXSEGFILE, 'NAM_VISC', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_VISC)

  CALL POSNAM( TPEXSEGFILE, 'NAM_FLYERS', GFOUND )
  IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_FLYERS)

  IF ( NAIRCRAFTS > 0 ) THEN
    CALL AIRCRAFTS_NML_ALLOCATE( NAIRCRAFTS )
    CALL POSNAM( TPEXSEGFILE, 'NAM_AIRCRAFTS', GFOUND )
    IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_AIRCRAFTS)
  END IF

  IF ( NBALLOONS > 0 ) THEN
    CALL BALLOONS_NML_ALLOCATE( NBALLOONS )
    CALL POSNAM( TPEXSEGFILE, 'NAM_BALLOONS', GFOUND )
    IF (GFOUND) READ(UNIT=ILUSEG,NML=NAM_BALLOONS)
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
CALL TEST_NAM_VAR(ILUOUT,'CPRESOPT',CPRESOPT,'RICHA','CGRAD','CRESI','ZRESI')
!
CALL TEST_NAM_VAR(ILUOUT,'CUVW_ADV_SCHEME',CUVW_ADV_SCHEME, &
       'CEN4TH','CEN2ND','WENO_K' ) 
CALL TEST_NAM_VAR(ILUOUT,'CMET_ADV_SCHEME',CMET_ADV_SCHEME, &
      &'PPM_00','PPM_01','PPM_02')
CALL TEST_NAM_VAR(ILUOUT,'CSV_ADV_SCHEME',CSV_ADV_SCHEME,   &
      &'PPM_00','PPM_01','PPM_02')
CALL TEST_NAM_VAR(ILUOUT,'CTEMP_SCHEME',CTEMP_SCHEME,       &
  &'RK11','RK21','RK33','RKC4','RK53','RK4B','RK62','RK65','NP32','SP32','LEFR')
!
CALL TEST_NAM_VAR(ILUOUT,'CTURB',CTURB,'NONE','TKEL')
CALL TEST_NAM_VAR(ILUOUT,'CRAD',CRAD,'NONE','FIXE','ECMW',&
#ifdef MNH_ECRAD
                 'ECRA',&
#endif
                 'TOPA')
CALL TEST_NAM_VAR(ILUOUT,'CCLOUD',CCLOUD,'NONE','REVE','KESS',  &
      & 'ICE3','ICE4','C2R2','C3R5','KHKO','LIMA')
CALL TEST_NAM_VAR(ILUOUT,'CDCONV',CDCONV,'NONE','KAFR')
CALL TEST_NAM_VAR(ILUOUT,'CSCONV',CSCONV,'NONE','KAFR','EDKF')
CALL TEST_NAM_VAR(ILUOUT,'CELEC',CELEC,'NONE','ELE3','ELE4')
!
CALL TEST_NAM_VAR(ILUOUT,'CAER',CAER,'TANR','TEGE','SURF','NONE')
CALL TEST_NAM_VAR(ILUOUT,'CAOP',CAOP,'CLIM','EXPL')
CALL TEST_NAM_VAR(ILUOUT,'CLW',CLW,'RRTM','MORC')
CALL TEST_NAM_VAR(ILUOUT,'CEFRADL',CEFRADL,'PRES','OCLN','MART','C2R2','LIMA')
CALL TEST_NAM_VAR(ILUOUT,'CEFRADI',CEFRADI,'FX40','LIOU','SURI','C3R5','LIMA')
CALL TEST_NAM_VAR(ILUOUT,'COPWLW',COPWLW,'SAVI','SMSH','LILI','MALA')
CALL TEST_NAM_VAR(ILUOUT,'COPILW',COPILW,'FULI','EBCU','SMSH','FU98')
CALL TEST_NAM_VAR(ILUOUT,'COPWSW',COPWSW,'SLIN','FOUQ','MALA')
CALL TEST_NAM_VAR(ILUOUT,'COPISW',COPISW,'FULI','EBCU','FU96')
!
CALL TEST_NAM_VAR(ILUOUT,'CLBCX(1)',CLBCX(1),'CYCL','WALL','OPEN')
CALL TEST_NAM_VAR(ILUOUT,'CLBCX(2)',CLBCX(2),'CYCL','WALL','OPEN')
CALL TEST_NAM_VAR(ILUOUT,'CLBCY(1)',CLBCY(1),'CYCL','WALL','OPEN')
CALL TEST_NAM_VAR(ILUOUT,'CLBCY(2)',CLBCY(2),'CYCL','WALL','OPEN')
!
CALL TURBN_INIT(CPROGRAM, TPEXSEGFILE, .FALSE., ILUOUT, .FALSE., .FALSE., .TRUE., 0)
CALL NEBN_INIT(CPROGRAM, TPEXSEGFILE, .FALSE., ILUOUT, .FALSE., .FALSE., .TRUE., 0)
!
CALL TEST_NAM_VAR(ILUOUT,'CCH_TDISCRETIZATION',CCH_TDISCRETIZATION, &
                 'SPLIT     ','CENTER    ','LAGGED    ')
!
CALL TEST_NAM_VAR(ILUOUT,'CCONF',CCONF,'START','RESTA')
CALL TEST_NAM_VAR(ILUOUT,'CEQNSYS',CEQNSYS,'LHE','DUR','MAE')
CALL TEST_NAM_VAR(ILUOUT,'CSPLIT',CSPLIT,'BSPLITTING','XSPLITTING','YSPLITTING')
!
CALL TEST_NAM_VAR(ILUOUT,'CBUTYPE',CBUTYPE,'NONE','CART','MASK')
!
CALL TEST_NAM_VAR(ILUOUT,'CRELAX_HEIGHT_TYPE',CRELAX_HEIGHT_TYPE,'FIXE','THGR')
!
CALL TEST_NAM_VAR(ILUOUT,'CLES_NORM_TYPE',CLES_NORM_TYPE,'NONE','CONV','EKMA','MOBU')
CALL TEST_NAM_VAR(ILUOUT,'CBL_HEIGHT_DEF',CBL_HEIGHT_DEF,'TKE','KE','WTV','FRI','DTH')
CALL TEST_NAM_VAR(ILUOUT,'CTURBLEN_CLOUD',CTURBLEN_CLOUD,'NONE','DEAR','DELT','BL89')
!
!   The test on the mass flux scheme for shallow convection
!
CALL PARAM_MFSHALLN_INIT(CPROGRAM, TPEXSEGFILE, .FALSE., ILUOUT, .FALSE., .FALSE., .TRUE., 0)
!
!   The test on the CSOLVER name is made elsewhere
!
CALL PARAM_ICEN_INIT(CPROGRAM, TPEXSEGFILE, .FALSE., ILUOUT, .FALSE., .FALSE., .TRUE., 0)
IF( CCLOUD == 'C3R5' ) THEN
  CALL TEST_NAM_VAR(ILUOUT,'CPRISTINE_ICE_C1R3',CPRISTINE_ICE_C1R3, &
                                                'PLAT','COLU','BURO')
  CALL TEST_NAM_VAR(ILUOUT,'CHEVRIMED_ICE_C1R3',CHEVRIMED_ICE_C1R3, &
                                                'GRAU','HAIL')
END IF
!
IF( CCLOUD == 'LIMA' ) THEN
  CALL PARAM_LIMA_INIT(CPROGRAM, TPEXSEGFILE, .FALSE., ILUOUT, .FALSE., .FALSE., .TRUE., 0)
END IF
! Blaze
CALL UPDATE_NAM_FIREn
IF (LBLAZE) THEN
  ! Blaze is only allowed on finer model(s)
  DO JI = 1, NMODEL
    IF ( JI /= KMI .AND. NDAD(JI) == KMI ) THEN
      WRITE( YMODEL, '( I3 )' ) JI
      CMNHMSG(1) = 'Blaze fire model only allowed on finer model'
      CMNHMSG(2) = '=> disabled on model ' // YMODEL
      CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'READ_EXSEG_n' )
      LBLAZE = .FALSE.
    END IF
  END DO
  CALL TEST_NAM_VAR(ILUOUT,'CPROPAG_MODEL',CPROPAG_MODEL,'SANTONI2011')
  CALL TEST_NAM_VAR(ILUOUT,'CHEAT_FLUX_MODEL',CHEAT_FLUX_MODEL,'CST','EXP','EXS')
  CALL TEST_NAM_VAR(ILUOUT,'CLATENT_FLUX_MODEL',CLATENT_FLUX_MODEL,'CST','EXP')
  CALL TEST_NAM_VAR(ILUOUT,'CFIRE_CPL_MODE',CFIRE_CPL_MODE,'2WAYCPL','FIR2ATM','ATM2FIR')
  CALL TEST_NAM_VAR(ILUOUT,'CWINDFILTER',CWINDFILTER,'EWAM','WLIM')
END IF
!
IF(LBLOWSNOW) THEN
       CALL TEST_NAM_VAR(ILUOUT,'CSNOWSEDIM',CSNOWSEDIM,'NONE','MITC','CARR','TABC')
       IF (XALPHA_SNOW .NE. 3 .AND. CSNOWSEDIM=='TABC') THEN
         WRITE(ILUOUT,*) '*****************************************'
         WRITE(ILUOUT,*) '* XALPHA_SNW must be set to 3 when                '
         WRITE(ILUOUT,*) '* CSNOWSEDIM = TABC                                 '
         WRITE(ILUOUT,*) '* Update the look-up table in BLOWSNOW_SEDIM_LKT1D    '
         WRITE(ILUOUT,*) '* to use TABC with a different value of XEMIALPHA_SNW'
         WRITE(ILUOUT,*) '*****************************************'
         !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
       ENDIF
END IF
! Consistency checks between phyex modules
IF ((CSUBG_AUCV_RC == 'ADJU' .OR. CSUBG_AUCV_RI == 'ADJU') .AND. CCONDENS /= 'GAUS') THEN 
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'READ_EXSEGN', & 
                &"CSUBG_AUCV_RC and/or CSUBG_AUCV_RI cannot be 'ADJU' if CCONDENS is not 'GAUS'") 
ENDIF
IF (.NOT. LHARAT .AND. LSTATNW) THEN
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'READ_EXSEGN', &
                &'LSTATNW only tested in combination with HARATU and EDMFm!')
ENDIF
!
!-------------------------------------------------------------------------------!
!*       2.    FIRST INITIALIZATIONS
!              ---------------------
!
!*       2.1   Time step in gridnesting case
!
IF (KMI /= 1 .AND. NDAD(KMI) /= KMI)  THEN
  XTSTEP = PTSTEP_ALL(NDAD(KMI)) / NDTRATIO(KMI)
END IF
PTSTEP_ALL(KMI) = XTSTEP
!
!*       2.2    Fill the global configuration module 
!
! Check coherence between the microphysical scheme and water species and 
!initialize the logicals LUSERn 
!
SELECT CASE ( CCLOUD )
  CASE ( 'NONE' )
    IF (.NOT. ( (.NOT. LUSERC) .AND. (.NOT. LUSERR) .AND. (.NOT. LUSERI) .AND. &
                (.NOT. LUSERS) .AND. (.NOT. LUSERG) .AND. (.NOT. LUSERH)       &
              ) .AND. CPROGRAM=='MESONH' )  THEN
!
      LUSERC=.FALSE.
      LUSERR=.FALSE.; LUSERI=.FALSE.
      LUSERS=.FALSE.; LUSERG=.FALSE. 
      LUSERH=.FALSE.
!
    END IF
!
    IF (CSUBG_AUCV_RC == 'SIGM')  THEN
!
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE THE SUBGRID AUTOCONVERSION SCHEME '                                         
      WRITE(UNIT=ILUOUT,FMT=*) ' WITHOUT MICROPHYSICS'                         
      WRITE(UNIT=ILUOUT,FMT=*) ' CSUBG_AUCV IS PUT TO "NONE"'                            
!
      CSUBG_AUCV_RC = 'NONE'
!
    END IF
!
  CASE ( 'REVE' ) 
    IF (.NOT. ( LUSERV .AND. LUSERC .AND. (.NOT. LUSERR) .AND. (.NOT. LUSERI) &
               .AND. (.NOT. LUSERS) .AND. (.NOT. LUSERG) .AND. (.NOT. LUSERH) &
              ) )  THEN
!
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE A REVERSIBLE MICROPHYSICAL "   ,&
      &" SCHEME. YOU WILL ONLY HAVE VAPOR AND CLOUD WATER ",/,                 &
      &" LUSERV AND LUSERC ARE TO TRUE AND THE OTHERS TO FALSE ")') 
!
      LUSERV=.TRUE. ; LUSERC=.TRUE.
      LUSERR=.FALSE.; LUSERI=.FALSE.
      LUSERS=.FALSE.; LUSERG=.FALSE. 
      LUSERH=.FALSE.
    END IF
!
    IF (CSUBG_AUCV_RC == 'SIGM')  THEN
!
      WRITE(UNIT=ILUOUT,FMT=9003) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE BOTH A REVERSIBLE MICROPHYSICAL SCHEME '
      WRITE(UNIT=ILUOUT,FMT=*) ' AND THE SUBGRID AUTOCONVERSION SCHEME '       
      WRITE(UNIT=ILUOUT,FMT=*) 'BUT YOU DO NOT HAVE RAIN in the "REVE" SCHEME'
      WRITE(UNIT=ILUOUT,FMT=*) ' CSUBG_AUCV_RC IS PUT TO "NONE"'                            
!
      CSUBG_AUCV_RC = 'NONE'
!
    END IF
!
  CASE ( 'KESS' )
    IF (.NOT. ( LUSERV .AND. LUSERC .AND. LUSERR .AND. (.NOT. LUSERI) .AND. &
              (.NOT. LUSERS) .AND. (.NOT. LUSERG) .AND. (.NOT. LUSERH)      &
              ) )  THEN
!
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE A KESSLER MICROPHYSICAL "   , &
      &" SCHEME. YOU WILL ONLY HAVE VAPOR, CLOUD WATER AND RAIN ",/,           &
      &" LUSERV, LUSERC AND LUSERR ARE SET TO TRUE AND THE OTHERS TO FALSE ")') 
!
      LUSERV=.TRUE. ; LUSERC=.TRUE. ; LUSERR=.TRUE.
      LUSERI=.FALSE.; LUSERS=.FALSE.
      LUSERG=.FALSE.; LUSERH=.FALSE.
    END IF
!
    IF (CSUBG_AUCV_RC == 'SIGM')  THEN
!
      WRITE(UNIT=ILUOUT,FMT=9003) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE BOTH A KESSLER MICROPHYSICAL SCHEME '
      WRITE(UNIT=ILUOUT,FMT=*) ' AND THE SUBGRID AUTOCONVERSION SCHEME USING'
      WRITE(UNIT=ILUOUT,FMT=*) 'SIGMA_RC.'
      WRITE(UNIT=ILUOUT,FMT=*) 'THIS IS NOT YET AVAILABLE.'
      WRITE(UNIT=ILUOUT,FMT=*) 'SET CSUBG_AUCV_RC TO "CLFR" or "NONE"  OR CCLOUD TO "ICE3"'
!
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
    END IF
!
  CASE ( 'ICE3' )
    IF (.NOT. ( LUSERV .AND. LUSERC .AND. LUSERR .AND. LUSERI .AND. LUSECI &
                       .AND. LUSERS .AND. LUSERG .AND. (.NOT. LUSERH))     &
                .AND. CPROGRAM=='MESONH' )  THEN
      !
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE THE ice3 SIMPLE MIXED PHASE'
      WRITE(UNIT=ILUOUT,FMT=*) 'MICROPHYSICAL SCHEME. YOU WILL ONLY HAVE VAPOR, CLOUD WATER,'
      WRITE(UNIT=ILUOUT,FMT=*) 'RAIN WATER, CLOUD ICE (MIXING RATIO AND CONCENTRATION)'
      WRITE(UNIT=ILUOUT,FMT=*) 'SNOW-AGGREGATES AND GRAUPELN.'
      WRITE(UNIT=ILUOUT,FMT=*) 'LUSERV,LUSERC,LUSERR,LUSERI,LUSECI,LUSERS,LUSERG ARE SET TO TRUE'
      WRITE(UNIT=ILUOUT,FMT=*) 'AND LUSERH TO FALSE'
!
      LUSERV=.TRUE. ; LUSERC=.TRUE. ; LUSERR=.TRUE.
      LUSERI=.TRUE. ; LUSECI=.TRUE.
      LUSERS=.TRUE. ; LUSERG=.TRUE. 
      LUSERH=.FALSE.
    END IF    
!
    IF (CSUBG_AUCV_RC == 'SIGM' .AND. .NOT. LSUBG_COND)  THEN
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE SUBGRID AUTOCONVERSION SCHEME'
      WRITE(UNIT=ILUOUT,FMT=*) ' WITHOUT THE SUBGRID CONDENSATION SCHEME.'
      WRITE(UNIT=ILUOUT,FMT=*) 'THIS IS NOT ALLOWED: CSUBG_AUCV_RC is SET to NONE' 
      CSUBG_AUCV_RC='NONE' 
    END IF
!
    IF (CSUBG_AUCV_RC == 'CLFR' .AND. CSCONV /= 'EDKF')  THEN
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE SUBGRID AUTOCONVERSION SCHEME'
      WRITE(UNIT=ILUOUT,FMT=*) 'WITH THE CONVECTIVE CLOUD FRACTION WITHOUT EDKF'
      WRITE(UNIT=ILUOUT,FMT=*) 'THIS IS NOT ALLOWED: CSUBG_AUCV_RC is SET to NONE' 
      CSUBG_AUCV_RC='NONE' 
    END IF
!
  CASE ( 'ICE4' )
    IF (.NOT. ( LUSERV .AND. LUSERC .AND. LUSERR .AND. LUSERI .AND. LUSECI &
                       .AND. LUSERS .AND. LUSERG .AND. LUSERH)             &
                .AND. CPROGRAM=='MESONH' )  THEN
      !
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE THE ice4 SIMPLE MIXED PHASE'
      WRITE(UNIT=ILUOUT,FMT=*) 'MICROPHYSICAL SCHEME. YOU WILL ONLY HAVE VAPOR, CLOUD WATER,'
      WRITE(UNIT=ILUOUT,FMT=*) 'RAIN WATER, CLOUD ICE (MIXING RATIO AND CONCENTRATION)'
      WRITE(UNIT=ILUOUT,FMT=*) 'SNOW-AGGREGATES, GRAUPELN AND HAILSTONES.'
      WRITE(UNIT=ILUOUT,FMT=*) 'LUSERV,LUSERC,LUSERR,LUSERI,LUSECI,LUSERS,LUSERG'
      WRITE(UNIT=ILUOUT,FMT=*) 'AND LUSERH ARE SET TO TRUE'
!
      LUSERV=.TRUE. ; LUSERC=.TRUE. ; LUSERR=.TRUE.
      LUSERI=.TRUE. ; LUSECI=.TRUE.
      LUSERS=.TRUE. ; LUSERG=.TRUE. ; LUSERH=.TRUE.
    END IF
!
    IF (CSUBG_AUCV_RC /= 'NONE' .AND. .NOT. LSUBG_COND)  THEN
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE SUBGRID AUTOCONVERSION SCHEME'
      WRITE(UNIT=ILUOUT,FMT=*) ' WITHOUT THE SUBGRID CONDENSATION SCHEME.'
      WRITE(UNIT=ILUOUT,FMT=*) 'THIS IS NOT ALLOWED: CSUBG_AUCV_RC is SET to NONE' 
      CSUBG_AUCV_RC='NONE' 
    END IF
!
  CASE ( 'C2R2','C3R5', 'KHKO' )
    IF (( EPARAM_CCN == 'XXX') .OR. (EINI_CCN  == 'XXX')) THEN
          WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE A 2-MOMENT MICROPHYSICAL ",    &
         &" SCHEME BUT YOU DIDNT FILL CORRECTLY NAM_PARAM_C2R2", &
         &" YOU HAVE TO FILL HPARAM_CCN and HINI_CCN ")')
 !callabortstop
        CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
    END IF
    IF (HCLOUD == 'NONE') THEN
      CGETCLOUD = 'SKIP'
    ELSE IF (HCLOUD == 'REVE' ) THEN
      CGETCLOUD = 'INI1'
    ELSE IF (HCLOUD == 'KESS' ) THEN
      CGETCLOUD = 'INI2'
    ELSE IF (HCLOUD == 'ICE3' ) THEN
      IF (CCLOUD == 'C3R5') THEN
        CGETCLOUD = 'INI2'
      ELSE
        WRITE(UNIT=ILUOUT,FMT=9003) KMI
        WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE WARM MICROPHYSICAL ",    &
        &" SCHEME BUT YOU WERE USING THE ICE3 SCHEME PREVIOUSLY.",/,          &
        &" AS THIS IS A LITTLE BIT STUPID IT IS NOT AUTHORIZED !!!")')
 !callabortstop
        CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
      END IF
    ELSE
      CGETCLOUD = 'READ' ! This is automatically done
    END IF
!
    IF ((CCLOUD == 'C2R2' ).OR. (CCLOUD == 'KHKO' )) THEN
      IF (.NOT. ( LUSERV .AND. LUSERC .AND. LUSERR .AND. (.NOT. LUSERI) .AND. &
                (.NOT. LUSERS) .AND. (.NOT. LUSERG) .AND. (.NOT. LUSERH)      &
                ) )  THEN
        WRITE(UNIT=ILUOUT,FMT=9002) KMI
        WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE C2R2 MICROPHYSICAL ",    &
        &" SCHEME. YOU WILL ONLY HAVE VAPOR, CLOUD WATER AND RAIN ",/,        &
        &"LUSERV, LUSERC AND LUSERR ARE SET TO TRUE AND THE OTHERS TO FALSE ")')
!
        LUSERV=.TRUE. ; LUSERC=.TRUE. ; LUSERR=.TRUE.
        LUSERI=.FALSE.; LUSERS=.FALSE.
        LUSERG=.FALSE.; LUSERH=.FALSE.
      END IF
    ELSE IF (CCLOUD == 'C3R5') THEN
      IF (.NOT. ( LUSERV .AND. LUSERC .AND. LUSERR .AND. LUSERI .AND. &
                  LUSERS .AND. LUSERG .AND. (.NOT. LUSERH)            &
                ) )  THEN
        WRITE(UNIT=ILUOUT,FMT=9002) KMI
        WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE C3R5 MICROPHYS. SCHEME.",&
        &" YOU WILL HAVE VAPOR, CLOUD WATER/ICE, RAIN, SNOW AND GRAUPEL ",/,  &
        &"LUSERV, LUSERC, LUSERR, LUSERI, LUSERS, LUSERG ARE SET TO TRUE")' )
!
        LUSERV=.TRUE. ; LUSERC=.TRUE. ; LUSERR=.TRUE.
        LUSERI=.TRUE. ; LUSECI=.TRUE.
        LUSERS=.TRUE. ; LUSERG=.TRUE.
        LUSERH=.FALSE.
      END IF
    ELSE IF (CCLOUD == 'LIMA') THEN
      IF (.NOT. ( LUSERV .AND. LUSERC .AND. LUSERR .AND. LUSERI .AND. &
                  LUSERS .AND. LUSERG .AND. (.NOT. LUSERH)            &
                ) )  THEN
        WRITE(UNIT=ILUOUT,FMT=9002) KMI
        WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE LIMA MICROPHYS. SCHEME.",&
        &" YOU WILL HAVE VAPOR, CLOUD WATER/ICE, RAIN, SNOW AND GRAUPEL ",/,  &
        &"LUSERV, LUSERC, LUSERR, LUSERI, LUSERS, LUSERG ARE SET TO TRUE")' )
!
        LUSERV=.TRUE. ; LUSERC=.TRUE. ; LUSERR=.TRUE.
        LUSERI=.TRUE. ; LUSECI=.TRUE.        
        LUSERS=.TRUE. ; LUSERG=.TRUE. 
        LUSERH=.FALSE.
      END IF
    END IF
!
    IF (LSUBG_COND)  THEN
      WRITE(UNIT=ILUOUT,FMT=9003) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE BOTH THE SIMPLE MIXED PHASE'
      WRITE(UNIT=ILUOUT,FMT=*) 'MICROPHYS. SCHEME AND THE SUBGRID COND. SCHEME.'
      WRITE(UNIT=ILUOUT,FMT=*) 'THIS IS NOT YET AVAILABLE.'
      WRITE(UNIT=ILUOUT,FMT=*) 'SET LSUBG_COND TO FALSE OR CCLOUD TO "REVE", "KESS"' 
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
    END IF
!
    IF ( CEFRADL /= 'C2R2') THEN
      WRITE(UNIT=ILUOUT,FMT=9001) KMI
      WRITE(UNIT=ILUOUT,FMT=*) ' YOU DID NOT CHOOSE CEFRADL=C2R2 FOR RADIATION'
      WRITE(UNIT=ILUOUT,FMT=*) ' IT IS ADVISED TO USE CEFRADL=C2R2 '                       
      WRITE(UNIT=ILUOUT,FMT=*) ' WITH A 2-MOMENT MICROPHYSICAL SCHEME'                     
    END IF
!
    IF ( CCLOUD == 'C3R5' .AND. CEFRADI /= 'C3R5') THEN
      WRITE(UNIT=ILUOUT,FMT=9001) KMI
      WRITE(UNIT=ILUOUT,FMT=*) ' YOU DID NOT CHOOSE CEFRADI=C3R5 FOR RADIATION'
      WRITE(UNIT=ILUOUT,FMT=*) ' IT IS ADVISED TO USE CEFRADI=C3R5 '                       
      WRITE(UNIT=ILUOUT,FMT=*) ' WITH A 2-MOMENT MICROPHYSICAL SCHEME'                     
    END IF
!
   IF ( WALPHAC /= 3.0 .OR. WNUC /= 2.0) THEN
      WRITE(UNIT=ILUOUT,FMT=9001) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'IT IS ADVISED TO USE XALPHAC=3. and XNUC=2.'
      WRITE(UNIT=ILUOUT,FMT=*) 'FOR STRATOCUMULUS WITH KHKO SCHEME. '
   END IF
!
   IF ( CEFRADL /= 'C2R2') THEN
      WRITE(UNIT=ILUOUT,FMT=9001) KMI
      WRITE(UNIT=ILUOUT,FMT=*) ' YOU DID NOT CHOOSE CEFRADL=C2R2 FOR RADIATION'
      WRITE(UNIT=ILUOUT,FMT=*) ' IT IS ADVISED TO USE CEFRADL=C2R2 '   
      WRITE(UNIT=ILUOUT,FMT=*) ' WITH A 2-MOMENT MICROPHYSICAL SCHEME'
   END IF
!
  CASE ( 'LIMA')
    IF (HCLOUD == 'NONE') THEN
      CGETCLOUD = 'SKIP'
    ELSE IF (HCLOUD == 'REVE' ) THEN
      CGETCLOUD = 'INI1'
    ELSE IF (HCLOUD == 'KESS' ) THEN
      CGETCLOUD = 'INI2'
    ELSE IF (HCLOUD == 'ICE3' ) THEN
      CGETCLOUD = 'INI2'
    ELSE
      CGETCLOUD = 'READ' ! This is automatically done
    END IF
!
    IF (NMOM_C.GE.1) THEN
      LUSERV=.TRUE. ; LUSERC=.TRUE. ; LUSERR=.TRUE.
      LUSERI=.FALSE.; LUSERS=.FALSE. ; LUSERG=.FALSE.; LUSERH=.FALSE.
    END IF
!
    IF (NMOM_I.GE.1) THEN
      LUSERV=.TRUE. ; LUSERC=.TRUE. ; LUSERR=.TRUE.
      LUSERI=.TRUE. ; LUSERS=.TRUE. ; LUSERG=.TRUE.
      LUSERH= NMOM_H.GE.1
    END IF
    !
    IF (LSPRO) LADJ=.FALSE.
    IF (.NOT.LPTSPLIT) THEN
       IF (NMOM_C==1) NMOM_C=2
       IF (NMOM_R==1) NMOM_R=2
       IF (NMOM_I==1) NMOM_I=2
       IF (NMOM_S==2 .OR. NMOM_G==2 .OR. NMOM_H==2) THEN
          NMOM_S=2
          NMOM_G=2
          IF (NMOM_H.GE.1) NMOM_H=2
       END IF
    END IF
!
    IF (LSUBG_COND .AND. (.NOT. LPTSPLIT)) THEN
      WRITE(UNIT=ILUOUT,FMT=9001) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'YOU MUST USE LPTSPLIT=T with CCLOUD=LIMA'
      WRITE(UNIT=ILUOUT,FMT=*) 'AND LSUBG_COND '
      CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','use LPTSPLIT=T with LIMA and LSUBG_COND=T')
    END IF
!
    IF (LSUBG_COND .AND. (.NOT. LADJ)) THEN
      WRITE(UNIT=ILUOUT,FMT=9001) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'YOU MUST USE LADJ=T with CCLOUD=LIMA'
      WRITE(UNIT=ILUOUT,FMT=*) 'AND LSUBG_COND '
      CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','use LADJ=T with LIMA and LSUBG_COND=T')
    END IF
!    
    IF ( LKHKO .AND. (XALPHAC /= 3.0 .OR. XNUC /= 2.0) ) THEN
      WRITE(UNIT=ILUOUT,FMT=9001) KMI
      WRITE(UNIT=ILUOUT,FMT=*) 'IT IS ADVISED TO USE XALPHAC=3. and XNUC=2.'
      WRITE(UNIT=ILUOUT,FMT=*) 'FOR STRATOCUMULUS. '
    END IF
!
    IF ( CEFRADL /= 'LIMA') THEN
      WRITE(UNIT=ILUOUT,FMT=9001) KMI
      WRITE(UNIT=ILUOUT,FMT=*) ' YOU DID NOT CHOOSE CEFRADL=LIMA FOR RADIATION'
      WRITE(UNIT=ILUOUT,FMT=*) ' IT IS ADVISED TO USE CEFRADL=LIMA '
      WRITE(UNIT=ILUOUT,FMT=*) ' WITH A 2-MOMENT MICROPHYSICAL SCHEME "LIMA"'
    END IF
!
END SELECT
!
LUSERV_G(KMI) = LUSERV
LUSERC_G(KMI) = LUSERC
LUSERR_G(KMI) = LUSERR
LUSERI_G(KMI) = LUSERI
LUSERS_G(KMI) = LUSERS
LUSERG_G(KMI) = LUSERG
LUSERH_G(KMI) = LUSERH
LUSETKE(KMI) = (CTURB /= 'NONE')
!
!-------------------------------------------------------------------------------
!
!*       2.3     Chemical and NSV_* variables initializations
!
CALL UPDATE_NAM_IBM_PARAMN
CALL UPDATE_NAM_RECYCL_PARAMN
CALL UPDATE_NAM_PARAMN
CALL UPDATE_NAM_DYNN
CALL UPDATE_NAM_CONFN
!
IF (LORILAM .AND. .NOT. LUSECHEM) THEN
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU CANNOT USE ORILAM AEROSOL SCHEME WITHOUT  '
  WRITE(ILUOUT,FMT=*) 'CHEMICAL GASEOUS CHEMISTRY                    '
  WRITE(ILUOUT,FMT=*) 'THEREFORE LUSECHEM IS SET TO TRUE    '
  LUSECHEM=.TRUE.
END IF
!
IF (LUSECHAQ.AND.(.NOT.LUSECHEM))  THEN
  WRITE(UNIT=ILUOUT,FMT=9003) KMI
  WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE AQUEOUS PHASE CHEMISTRY'
  WRITE(UNIT=ILUOUT,FMT=*) 'BUT THE CHEMISTRY IS NOT ACTIVATED'
  WRITE(UNIT=ILUOUT,FMT=*) 'SET LUSECHEM TO TRUE IF YOU WANT REALLY USE CHEMISTRY' 
  WRITE(UNIT=ILUOUT,FMT=*) 'OR SET LUSECHAQ TO FALSE IF YOU DO NOT WANT USE IT' 
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
END IF
IF (LUSECHAQ.AND.(.NOT.LUSERC).AND.CPROGRAM=='MESONH')  THEN
  WRITE(UNIT=ILUOUT,FMT=9003) KMI
  WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE AQUEOUS PHASE CHEMISTRY'
  WRITE(UNIT=ILUOUT,FMT=*) 'BUT CLOUD MICROPHYSICS IS NOT ACTIVATED'
  WRITE(UNIT=ILUOUT,FMT=*) 'LUSECHAQ IS SET TO FALSE' 
  LUSECHAQ = .FALSE.
END IF
IF (LUSECHAQ.AND.CCLOUD(1:3) == 'ICE'.AND. .NOT. LUSECHIC) THEN
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE AQUEOUS PHASE CHEMISTRY'
  WRITE(UNIT=ILUOUT,FMT=*) 'WITH MIXED PHASE CLOUD MICROPHYSICS'
  WRITE(UNIT=ILUOUT,FMT=*) 'SET LUSECHIC TO TRUE IF YOU WANT TO ACTIVATE'
  WRITE(UNIT=ILUOUT,FMT=*) 'ICE PHASE CHEMICAL SPECIES'
  IF (LCH_RET_ICE) THEN
    WRITE(UNIT=ILUOUT,FMT=*) 'LCH_RET_ICE TRUE MEANS ALL SOLUBLE'
    WRITE(UNIT=ILUOUT,FMT=*) 'GASES ARE RETAINED IN ICE PHASE'
    WRITE(UNIT=ILUOUT,FMT=*) 'WHEN SUPERCOOLED WATER FREEZES'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=*) 'LCH_RET_ICE FALSE MEANS ALL SOLUBLE'
    WRITE(UNIT=ILUOUT,FMT=*) 'GASES GO BACK TO THE GAS PHASE WHEN'
    WRITE(UNIT=ILUOUT,FMT=*) 'SUPERCOOLED WATER FREEZES'
  ENDIF
ENDIF
IF (LUSECHIC.AND. .NOT. CCLOUD(1:3) == 'ICE'.AND.CPROGRAM=='MESONH') THEN
  WRITE(UNIT=ILUOUT,FMT=9003) KMI
  WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE ICE PHASE CHEMISTRY'
  WRITE(UNIT=ILUOUT,FMT=*) 'BUT MIXED PHASE CLOUD MICROPHYSICS IS NOT ACTIVATED'
  WRITE(UNIT=ILUOUT,FMT=*) 'LUSECHIC IS SET TO FALSE'  
  LUSECHIC= .FALSE.
ENDIF
IF (LCH_PH.AND. (.NOT. LUSECHAQ)) THEN
  WRITE(UNIT=ILUOUT,FMT=9003) KMI
  WRITE(UNIT=ILUOUT,FMT=*) 'DIAGNOSTIC PH COMPUTATION IS ACTIVATED'
  WRITE(UNIT=ILUOUT,FMT=*) 'BUT AQUEOUS PHASE CHEMISTRY IS NOT ACTIVATED'
  WRITE(UNIT=ILUOUT,FMT=*) 'SET LUSECHAQ TO TRUE IF YOU WANT TO ACTIVATE IT'
  WRITE(UNIT=ILUOUT,FMT=*) 'LCH_PH IS SET TO FALSE'
  LCH_PH= .FALSE.
ENDIF
IF (LUSECHIC.AND.(.NOT.LUSECHAQ))  THEN
  WRITE(UNIT=ILUOUT,FMT=9003) KMI
  WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE ICE PHASE CHEMISTRY'
  WRITE(UNIT=ILUOUT,FMT=*) 'BUT THE AQUEOUS PHASE CHEMISTRY IS NOT ACTIVATED'
  WRITE(UNIT=ILUOUT,FMT=*) 'SET LUSECHAQ TO TRUE IF YOU WANT REALLY USE CLOUD CHEMISTRY'
  WRITE(UNIT=ILUOUT,FMT=*) 'OR SET LUSECHIC TO FALSE IF YOU DO NOT WANT USE IT'
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
END IF
IF ((LUSECHIC).AND.(LCH_RET_ICE)) THEN
  WRITE(UNIT=ILUOUT,FMT=9003) KMI
  WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE RETENTION OF SOLUBLE GASES IN ICE'
  WRITE(UNIT=ILUOUT,FMT=*) 'BUT THE ICE PHASE CHEMISTRY IS ACTIVATED'
  WRITE(UNIT=ILUOUT,FMT=*) 'FLAG LCH_RET_ICE IS ONLY USES WHEN LUSECHIC IS SET'
  WRITE(UNIT=ILUOUT,FMT=*) 'TO FALSE IE NO CHEMICAL SPECIES IN ICE' 
ENDIF
!
CALL UPDATE_NAM_CH_MNHCN
CALL INI_NSV(KMI)
!
! From this point, all NSV* variables contain valid values for model KMI
! 
DO JSV = 1,NSV
 LUSESV(JSV,KMI) = .TRUE.
END DO
!
IF ( CAOP=='EXPL' .AND. .NOT.LDUST .AND. .NOT.LORILAM          &
                  .AND. .NOT.LSALT .AND. .NOT.(CCLOUD=='LIMA') ) THEN
      WRITE(UNIT=ILUOUT,FMT=9001) KMI
      WRITE(UNIT=ILUOUT,FMT=*) ' YOU WANT TO USE EXPLICIT AEROSOL OPTICAL '       
      WRITE(UNIT=ILUOUT,FMT=*) 'PROPERTIES BUT YOU DONT HAVE DUST OR '            
      WRITE(UNIT=ILUOUT,FMT=*) 'AEROSOL OR SALT THEREFORE CAOP=CLIM'            
      CAOP='CLIM'
END IF      
!-------------------------------------------------------------------------------
!
!*       3.    CHECK COHERENCE BETWEEN EXSEG VARIABLES AND FMFILE ATTRIBUTES
!              -------------------------------------------------------------
!
!
!*       3.1  Turbulence variable 
!
IF ((CTURB /= 'NONE').AND.(HTURB == 'NONE')) THEN
  CGETTKET ='INIT'
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT=*)'YOU WANT TO USE TURBULENCE KINETIC ENERGY TKE'
  WRITE(UNIT=ILUOUT,FMT=*)'WHEREAS IT  IS NOT IN INITIAL FMFILE'
  WRITE(UNIT=ILUOUT,FMT=*)'TKE WILL BE INITIALIZED TO ZERO'
ELSE 
  IF (CTURB /= 'NONE') THEN 
    CGETTKET ='READ'
    IF ((CCONF=='START') .AND. CPROGRAM /= 'DIAG') CGETTKET='INIT' 
  ELSE
   CGETTKET ='SKIP'
  END IF
END IF
!
!
IF ((CTOM == 'TM06').AND.(HTOM /= 'TM06')) THEN
  CGETBL_DEPTH ='INIT'
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT=*)'YOU WANT TO USE BL DEPTH FOR THIRD ORDER MOMENTS'
  WRITE(UNIT=ILUOUT,FMT=*)'WHEREAS IT IS NOT IN INITIAL FMFILE'
  WRITE(UNIT=ILUOUT,FMT=*)'IT WILL BE INITIALIZED TO ZERO'
ELSE 
  IF (CTOM == 'TM06') THEN 
    CGETBL_DEPTH ='READ'
  ELSE 
    CGETBL_DEPTH ='SKIP'
  END IF
END IF
!
IF (LRMC01 .AND. .NOT. ORMC01) THEN
  CGETSBL_DEPTH ='INIT'
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT=*)'YOU WANT TO USE SBL DEPTH FOR RMC01'
  WRITE(UNIT=ILUOUT,FMT=*)'WHEREAS IT IS NOT IN INITIAL FMFILE'
  WRITE(UNIT=ILUOUT,FMT=*)'IT WILL BE INITIALIZED TO ZERO'
ELSE 
  IF (LRMC01) THEN 
    CGETSBL_DEPTH ='READ'
  ELSE 
    CGETSBL_DEPTH ='SKIP'
  END IF
END IF
!
!
!*       3.2  Moist  variables 
!
IF (LUSERV.AND. (.NOT.OUSERV)) THEN
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE VAPOR VARIABLE Rv WHEREAS IT  ", &
                  & "IS NOT IN INITIAL FMFILE",/,                           &
                  & "Rv WILL BE INITIALIZED TO ZERO")')
  CGETRVT='INIT'
ELSE                                                                
  IF (LUSERV) THEN
    CGETRVT='READ'
  ELSE
    CGETRVT='SKIP'
  END IF
END IF
!
IF (LUSERC.AND. (.NOT.OUSERC)) THEN
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE CLOUD VARIABLE Rc WHEREAS IT ",  &
                 &   " IS NOT IN INITIAL FMFILE",/,                        &
                 &   "Rc WILL BE INITIALIZED TO ZERO")')
  CGETRCT='INIT'
ELSE
  IF (LUSERC) THEN
    CGETRCT='READ'
!   IF(CCONF=='START') CGETRCT='INIT' 
  ELSE
    CGETRCT='SKIP'
  END IF
END IF
!
IF (LUSERR.AND. (.NOT.OUSERR)) THEN
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE RAIN VARIABLE Rr WHEREAS IT ", &
                  &  "IS NOT IN INITIAL FMFILE",/,                     &
                  &  " Rr WILL BE INITIALIZED TO ZERO")')
 
  CGETRRT='INIT'
ELSE
  IF (LUSERR) THEN
    CGETRRT='READ'
!   IF( (CCONF=='START').AND. CPROGRAM /= 'DIAG') CGETRRT='INIT' 
  ELSE
    CGETRRT='SKIP'
  END IF
END IF
!
IF (LUSERI.AND. (.NOT.OUSERI)) THEN
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE ICE VARIABLE Ri WHEREAS IT ", &
               &    "IS NOT IN INITIAL FMFILE",/,                      &
               &     " Ri WILL BE INITIALIZED TO ZERO")')
  CGETRIT='INIT'
ELSE
  IF (LUSERI) THEN
    CGETRIT='READ'
!   IF(CCONF=='START') CGETRIT='INIT' 
  ELSE
    CGETRIT='SKIP'
  END IF
END IF
!
IF (LUSECI.AND. (.NOT.OUSECI)) THEN
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE ICE CONC. VARIABLE Ci WHEREAS IT ",&
               &          "IS NOT IN INITIAL FMFILE",/,                       &
               &          "   Ci WILL BE INITIALIZED TO ZERO")')
  CGETCIT='INIT'
ELSE
  IF (LUSECI) THEN
    CGETCIT='READ'
  ELSE
    CGETCIT='SKIP'
  END IF
END IF
!
IF (LUSERS.AND. (.NOT.OUSERS)) THEN
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE SNOW VARIABLE Rs WHEREAS IT ",&
                  &  "IS NOT IN INITIAL FMFILE",/,                       &
                  &  " Rs WILL BE INITIALIZED TO ZERO")')
  CGETRST='INIT'
ELSE
  IF (LUSERS) THEN
    CGETRST='READ'
!   IF ( (CCONF=='START').AND. CPROGRAM /= 'DIAG') CGETRST='INIT' 
  ELSE
    CGETRST='SKIP'
  END IF
END IF
!
IF (LUSERG.AND. (.NOT.OUSERG)) THEN
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE GRAUPEL VARIABLE Rg WHEREAS ",&
                   & " IT IS NOTIN INITIAL FMFILE",/,                    &
                   & "Rg WILL BE INITIALIZED TO ZERO")')
  CGETRGT='INIT'
ELSE
  IF (LUSERG) THEN
    CGETRGT='READ'
!   IF ( (CCONF=='START') .AND. CPROGRAM /= 'DIAG') CGETRGT='INIT' 
  ELSE
    CGETRGT='SKIP'
  END IF
END IF
!
IF (LUSERH.AND. (.NOT.OUSERH)) THEN
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE HAIL VARIABLE Rh WHEREAS",&
                  &  "IT IS NOT IN INITIAL FMFILE",/,                &
                  & " Rh WILL BE INITIALIZED TO ZERO")')
   CGETRHT='INIT'
ELSE
  IF (LUSERH) THEN
    CGETRHT='READ'
!   IF ( (CCONF=='START') .AND. CPROGRAM /= 'DIAG') CGETRHT='INIT' 
  ELSE
    CGETRHT='SKIP'
  END IF
END IF
!
IF (LUSERC.AND. (.NOT.OUSERC)) THEN
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT=*) 'THE CLOUD FRACTION WILL BE INITIALIZED ACCORDING'
  WRITE(UNIT=ILUOUT,FMT=*) 'TO CLOUD MIXING RATIO VALUE OR SET TO 0'
  CGETCLDFR = 'INIT'
ELSE
  IF ( LUSERC ) THEN
    CGETCLDFR = 'READ'
    IF ( (CCONF=='START') .AND. CPROGRAM /= 'DIAG') CGETCLDFR='INIT' 
  ELSE
    CGETCLDFR = 'SKIP'
  END IF
END IF
!
IF (LUSERI.AND. (.NOT.OUSERI)) THEN
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT=*) 'THE ICE CLOUD FRACTION WILL BE INITIALIZED ACCORDING'
  WRITE(UNIT=ILUOUT,FMT=*) 'TO CLOUD MIXING RATIO VALUE OR SET TO 0'
  CGETICEFR = 'INIT'
ELSE
  IF ( LUSERI ) THEN
    CGETICEFR = 'READ'
    IF ( (CCONF=='START') .AND. CPROGRAM /= 'DIAG') CGETICEFR='INIT' 
  ELSE
    CGETICEFR = 'SKIP'
  END IF
END IF
!
!
!*       3.3  Moist turbulence
!
IF ( LUSERC .AND. CTURB /= 'NONE' ) THEN
  IF ( .NOT. (OUSERC .AND. HTURB /= 'NONE') ) THEN
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE MOIST TURBULENCE WHEREAS IT ",/, &
                 &   " WAS NOT THE CASE FOR THE INITIAL FMFILE GENERATION",/, &
                 &    "SRC AND SIGS ARE INITIALIZED TO 0")')
    CGETSRCT ='INIT'
    CGETSIGS ='INIT'
  ELSE
    CGETSRCT ='READ'
    IF ( (CCONF=='START') .AND. CPROGRAM /= 'DIAG') CGETSRCT ='INIT'
    CGETSIGS ='READ'
  END IF
ELSE
  CGETSRCT ='SKIP'
  CGETSIGS ='SKIP'
END IF
!
IF(LCLOUDMODIFLM .AND. CTURBLEN_CLOUD/='NONE') THEN
  IF (CTURB=='NONE' .OR. .NOT.LUSERC) THEN
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO COMPUTE A MIXING LENGTH FOR CLOUD=", &
                 &   ", WHEREAS YOU DO NOT SPECIFY A TURBULENCE SCHEME OR ",  &
                 &   "USE OF RC,",/," CTURBLEN_CLOUD IS SET TO NONE")')       &
                 CTURBLEN_CLOUD
    CTURBLEN_CLOUD='NONE'
  END IF
  IF( XCEI_MIN > XCEI_MAX ) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT='("PROBLEM OF CEI LIMITS FOR CLOUD MIXING  ",/, &
                 &   "LENGTH COMPUTATION: XCEI_MIN=",E9.3,", XCEI_MAX=",E9.3)')&
                 XCEI_MIN,XCEI_MAX
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
END IF
!
IF ( LSIGMAS ) THEN
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE SIGMA_S FROM TURBULENCE SCHEME",/, &
                 &   " IN ICE SUBGRID CONDENSATION, SO YOUR SIGMA_S"/, &
                 &   " MIGHT BE SMALL ABOVE PBL DEPENDING ON LENGTH SCALE")')
END IF
!
IF (LSUBG_COND .AND. CTURB=='NONE' ) THEN
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE SUBGRID CONDENSATION'
  WRITE(UNIT=ILUOUT,FMT=*) ' WITHOUT TURBULENCE '                           
  WRITE(UNIT=ILUOUT,FMT=*) 'THIS IS NOT ALLOWED: LSUBG_COND is SET to FALSE'
  LSUBG_COND=.FALSE.
END IF
!
IF (L1D .AND. CTURB/='NONE' .AND. CTURBDIM == '3DIM') THEN
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE 3D TURBULENCE IN 1D CONFIGURATION '
  WRITE(UNIT=ILUOUT,FMT=*) 'THIS IS NOT POSSIBLE: CTURBDIM IS SET TO 1DIM'
  CTURBDIM = '1DIM'
END IF
!
!*       3.4  Additional scalar variables 
!
IF (NSV_USER == KSV_USER) THEN
  DO  JS = 1,KSV_USER             ! to read all the variables in initial file 
    CGETSVT(JS)='READ'            ! and to initialize them 
!   IF(CCONF=='START')CGETSVT(JS)='INIT'       ! with  these values  
  END DO
ELSEIF (NSV_USER > KSV_USER) THEN
 IF (KSV_USER == 0) THEN
    CGETSVT(1:NSV_USER)='INIT'
 ELSE
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE MORE ADDITIONAL SCALAR  " ,&
  &" VARIABLES THAN THERE ARE IN INITIAL FMFILE",/,                  &
  & "THE SUPPLEMENTARY VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
  DO  JS = 1,KSV_USER             ! to read all the variables in initial file 
    CGETSVT(JS)='READ'            ! and to initialize them
!   IF(CCONF=='START')CGETSVT(JS)='INIT'        ! with  these values
  END DO
  DO JS = KSV_USER+1, NSV_USER    ! to initialize to zero supplementary
    CGETSVT(JS)='INIT'            ! initial file)
  END DO
 END IF
ELSE
  WRITE(UNIT=ILUOUT,FMT=9000) KMI
  WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE LESS ADDITIONAL SCALAR  " ,&
  &" VARIABLES THAN THERE ARE IN INITIAL FMFILE")')
  DO  JS = 1,NSV_USER             ! to read the first NSV_USER variables in initial file 
    CGETSVT(JS)='READ'            ! and to initialize with these values
!   IF(CCONF=='START') CGETSVT(JS)='INIT'  
  END DO 
  DO  JS = NSV_USER + 1, KSV_USER ! to skip the last (KSV_USER-NSV_USER) variables
    CGETSVT(JS)='SKIP' 
  END DO
END IF
!
! C2R2 and KHKO SV case
!
IF (CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'KHKO') THEN 
  IF (HCLOUD == 'C2R2' .OR. HCLOUD == 'C3R5' .OR. HCLOUD == 'KHKO') THEN
    CGETSVT(NSV_C2R2BEG:NSV_C2R2END)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_C2R2BEG:NSV_C2R2END)='INIT'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR C2R2 &
         & (or KHKO) SCHEME IN INITIAL FMFILE",/,&
         & "THE C2R2 (or KHKO) VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_C2R2BEG:NSV_C2R2END)='INIT'
  END IF
END IF
!
! C3R5 SV case
!
IF (CCLOUD == 'C3R5') THEN 
  IF (HCLOUD == 'C3R5') THEN
    CGETSVT(NSV_C1R3BEG:NSV_C1R3END)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_C1R3BEG:NSV_C1R3END)='INIT'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR C3R5 &
         &SCHEME IN INITIAL FMFILE",/,&
         & "THE C1R3 VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_C1R3BEG:NSV_C1R3END)='INIT'
  END IF
END IF
!
! LIMA SV case
!
IF (CCLOUD == 'LIMA') THEN
  IF (HCLOUD == 'LIMA') THEN
    CGETSVT(NSV_LIMA_BEG:NSV_LIMA_END)='READ'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR LIMA &
           & SCHEME IN INITIAL FMFILE",/,&
           & "THE LIMA VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_LIMA_BEG:NSV_LIMA_END)='INIT'
  END IF
END IF
!
! Electrical SV case
!
IF (CELEC /= 'NONE') THEN 
  IF (HELEC /= 'NONE') THEN
    CGETSVT(NSV_ELECBEG:NSV_ELECEND)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_ELECBEG:NSV_ELECEND)='INIT'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR ELECTRICAL &
         &SCHEME IN INITIAL FMFILE",/,&
         & "THE ELECTRICAL VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_ELECBEG:NSV_ELECEND)='INIT'
  END IF
END IF
!
! (explicit) LINOx SV case 
!
IF (CELEC /= 'NONE' .AND. LLNOX_EXPLICIT) THEN
  IF (HELEC /= 'NONE' .AND. OLNOX_EXPLICIT) THEN
    CGETSVT(NSV_LNOXBEG:NSV_LNOXEND)='READ' 
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR LINOX &
         & IN INITIAL FMFILE",/,& 
         & "THE LINOX VARIABLES HAVE BEEN INITIALIZED TO ZERO ")')  
    CGETSVT(NSV_LNOXBEG:NSV_LNOXEND)='INIT' 
  END IF
END IF
!
! Chemical SV case (excluding aqueous chemical species)
!
IF (LUSECHEM) THEN
  IF (OUSECHEM) THEN
    CGETSVT(NSV_CHGSBEG:NSV_CHGSEND)='READ'
    IF(CCONF=='START' .AND. LCH_INIT_FIELD ) CGETSVT(NSV_CHGSBEG:NSV_CHGSEND)='INIT'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR CHEMICAL &
         &SCHEME IN INITIAL FMFILE",/,&
         & "THE CHEMICAL VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_CHGSBEG:NSV_CHGSEND)='INIT'
  END IF
END IF
! add aqueous chemical species
IF (LUSECHAQ) THEN
  IF (OUSECHAQ) THEN
    CGETSVT(NSV_CHACBEG:NSV_CHACEND)='READ'
!    IF(CCONF=='START') CGETSVT(NSV_CHACBEG:NSV_CHACEND)='INIT'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR CHEMICAL &
         &SCHEME IN AQUEOUS PHASE IN INITIAL FMFILE",/,&
         & "THE AQUEOUS PHASE CHEMICAL VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_CHACBEG:NSV_CHACEND)='INIT'
  END IF
END IF
! add ice phase chemical species
IF (LUSECHIC) THEN
  IF (OUSECHIC) THEN
    CGETSVT(NSV_CHICBEG:NSV_CHICEND)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_CHICBEG:NSV_CHICEND)='INIT'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR CHEMICAL &
         &SPECIES IN ICE PHASE IN INITIAL FMFILE",/,&
         & "THE ICE PHASE CHEMICAL VARIABLES HAVE BEEN INITIALIZED TO ZERO ")')
    CGETSVT(NSV_CHICBEG:NSV_CHICEND)='INIT'
  END IF
END IF
! pH values = diagnostics
IF (LCH_PH .AND. .NOT. OCH_PH) THEN
  CGETPHC ='INIT'  !will be initialized to XCH_PHINIT
  IF (LUSERR) THEN
    CGETPHR = 'INIT' !idem
  ELSE 
    CGETPHR = 'SKIP'
  ENDIF
ELSE
  IF (LCH_PH) THEN
    CGETPHC ='READ'
    IF (LUSERR) THEN
      CGETPHR = 'READ'
    ELSE
      CGETPHR = 'SKIP'
    ENDIF
  ELSE
    CGETPHC ='SKIP'
    CGETPHR ='SKIP'
  END IF
END IF
!
! Dust case
!
IF (LDUST) THEN
  IF (ODUST) THEN
    CGETSVT(NSV_DSTBEG:NSV_DSTEND)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_DSTBEG:NSV_DSTEND)='INIT'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR DUST &
         &SCHEME IN INITIAL FMFILE",/,&
         & "THE DUST VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_DSTBEG:NSV_DSTEND)='INIT'
  END IF
  IF (LDEPOS_DST(KMI)) THEN

          !UPG *PT
  IF((CCLOUD /= 'ICE3').AND.(CCLOUD /= 'ICE4').AND.(CCLOUD /= 'KESS')&
  .AND.(CCLOUD /= 'KHKO').AND.(CCLOUD /= 'C2R2').AND.(CCLOUD /= 'LIMA').AND. &
       (CPROGRAM/='SPAWN').AND.(CPROGRAM/='REAL'))  THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT='("ERROR: WET DEPOSITION OF DUST IS ONLY CODED FOR THE",/,&
         & "MICROPHYSICAL SCHEME as ICE3, ICE4, KESS, KHKO, LIMA and C2R2")') 
          !UPG *PT
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF  

   IF (ODEPOS_DST(KMI) ) THEN
    CGETSVT(NSV_DSTDEPBEG:NSV_DSTDEPEND)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_DSTDEPBEG:NSV_DSTDEPEND)='INIT'    
   ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR RAIN and CLOUD DUST &
         &  SCHEME IN INITIAL FMFILE",/,&
         & "THE MOIST DUST VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_DSTDEPBEG:NSV_DSTDEPEND)='INIT'    
   END IF
  END IF  

  IF(NMODE_DST.GT.3 .OR. NMODE_DST.LT.1) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT='("DUST MODES MUST BE BETWEEN 1 and 3 ")') 
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF     
END IF
!
! Sea Salt case
!
IF (LSALT) THEN
  IF (OSALT) THEN
    CGETSVT(NSV_SLTBEG:NSV_SLTEND)='READ'
    CGETZWS='READ'
!   IF(CCONF=='START') CGETSVT(NSV_SLTBEG:NSV_SLTEND)='INIT'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR SALT &
         &SCHEME IN INITIAL FMFILE",/,&
         & "THE SALT VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_SLTBEG:NSV_SLTEND)='INIT'
    CGETZWS='INIT'
  END IF
  IF (LDEPOS_SLT(KMI)) THEN

          !UPG*PT
  IF((CCLOUD /= 'ICE3').AND.(CCLOUD /= 'ICE4').AND.(CCLOUD /= 'KESS')&
  !.AND.(CCLOUD /= 'KHKO').AND.(CCLOUD /= 'C2R2').AND.                &
  .AND.(CCLOUD /= 'KHKO').AND.(CCLOUD /= 'C2R2').AND.(CCLOUD /= 'LIMA').AND. &
       (CPROGRAM/='SPAWN').AND.(CPROGRAM/='REAL'))  THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT='("ERROR: WET DEPOSITION OF SEA SALT AEROSOLS IS ONLY CODED FOR THE",/,&
         & "MICROPHYSICAL SCHEME as ICE3, ICE4, KESS, KHKO, LIMA and C2R2")') 
          !UPG*PT
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF  

   IF (ODEPOS_SLT(KMI) ) THEN
    CGETSVT(NSV_SLTDEPBEG:NSV_SLTDEPEND)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_SLTDEPBEG:NSV_SLTDEPEND)='INIT'    
   ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR RAIN and CLOUD SEA SALT &
         &  SCHEME IN INITIAL FMFILE",/,&
         & "THE MOIST SEA SALT VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_SLTDEPBEG:NSV_SLTDEPEND)='INIT'    
   END IF
  END IF
  IF(NMODE_SLT.GT.8 .OR. NMODE_SLT.LT.1) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT='("SALT MODES MUST BE BETWEEN 1 and 8 ")') 
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF     
END IF 
!
! Orilam SV case
!
IF (LORILAM) THEN
  IF (OORILAM) THEN
    CGETSVT(NSV_AERBEG:NSV_AEREND)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_AERBEG:NSV_AEREND)='INIT'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR AEROSOL &
         &SCHEME IN INITIAL FMFILE",/,&
         & "THE AEROSOLS VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_AERBEG:NSV_AEREND)='INIT'
  END IF
  IF (LDEPOS_AER(KMI)) THEN

          !UPG*PT
  IF((CCLOUD /= 'ICE3').AND.(CCLOUD /= 'ICE4').AND.(CCLOUD /= 'KESS')&
  .AND.(CCLOUD /= 'KHKO').AND.(CCLOUD /= 'C2R2').AND.(CCLOUD /= 'LIMA').AND. &
  !.AND.(CCLOUD /= 'KHKO').AND.(CCLOUD /= 'C2R2').AND.                &
       (CPROGRAM/='SPAWN').AND.(CPROGRAM/='REAL'))  THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT='("ERROR: WET DEPOSITION OF ORILAM AEROSOLS IS ONLY CODED FOR THE",/,&
         & "MICROPHYSICAL SCHEME as ICE3, ICE4, KESS, KHKO, LIMA and C2R2")') 
          !UPG*PT
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF  

   IF (ODEPOS_AER(KMI) ) THEN
    CGETSVT(NSV_AERDEPBEG:NSV_AERDEPEND)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_AERDEPBEG:NSV_AERDEPEND)='INIT'    
   ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR RAIN and IN CLOUD  &
         &  AEROSOL SCHEME IN INITIAL FMFILE",/,&
         & "THE MOIST AEROSOL VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_AERDEPBEG:NSV_AERDEPEND)='INIT'    
   END IF
  END IF
END IF
!
! Lagrangian variables
!
IF (LINIT_LG .AND. .NOT.(LLG)) THEN
      WRITE(UNIT=ILUOUT,FMT=9001) KMI
      WRITE(UNIT=ILUOUT,FMT='("IT IS INCOHERENT TO HAVE LINIT_LG=.T. AND LLG=.F.",/,&
      & "IF YOU WANT LAGRANGIAN TRACERS CHANGE LLG TO .T. ")')
ENDIF
IF (LLG) THEN
  IF (OLG .AND. .NOT.(LINIT_LG .AND. CPROGRAM=='MESONH')) THEN
    CGETSVT(NSV_LGBEG:NSV_LGEND)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_LGBEG:NSV_LGEND)='INIT'
  ELSE
    IF(.NOT.(LINIT_LG) .AND. CPROGRAM=='MESONH') THEN
      WRITE(UNIT=ILUOUT,FMT=9001) KMI
      WRITE(UNIT=ILUOUT,FMT='("THERE IS NO LAGRANGIAN VARIABLES IN INITIAL FMFILE",/,&
                       & "THE LAGRANGIAN VARIABLES HAVE BEEN REINITIALIZED")')
      LINIT_LG=.TRUE.
    ENDIF
    CGETSVT(NSV_LGBEG:NSV_LGEND)='INIT'
  END IF
END IF
!
!
! LINOx SV case
!
IF (.NOT.LUSECHEM .AND. LCH_CONV_LINOX) THEN
  IF (.NOT.OUSECHEM .AND. OCH_CONV_LINOX) THEN
    CGETSVT(NSV_LNOXBEG:NSV_LNOXEND)='READ'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9002) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR LINOX &
         &IN INITIAL FMFILE",/,&
         & "THE LINOX VARIABLES HAVE BEEN INITIALIZED TO ZERO ")') 
    CGETSVT(NSV_LNOXBEG:NSV_LNOXEND)='INIT'
  END IF
END IF
!
! Passive pollutant case
!
IF (LPASPOL) THEN
  IF (OPASPOL) THEN
    CGETSVT(NSV_PPBEG:NSV_PPEND)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_PPBEG:NSV_PPEND)='INIT'    
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO PASSIVE SCALAR VARIABLES IN INITIAL FMFILE",/,&
                       & "THE VARIABLES HAVE BEEN INITIALIZED TO ZERO")')
    CGETSVT(NSV_PPBEG:NSV_PPEND)='INIT'
  END IF
END IF
!
#ifdef MNH_FOREFIRE
! ForeFire
!
IF (LFOREFIRE) THEN
  IF (OFOREFIRE) THEN
    CGETSVT(NSV_FFBEG:NSV_FFEND)='READ'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO FOREFIRE SCALAR VARIABLES IN INITIAL FMFILE",/,&
                       & "THE VARIABLES HAVE BEEN INITIALIZED TO ZERO")')
    CGETSVT(NSV_FFBEG:NSV_FFEND)='INIT'
  END IF
END IF
#endif
! Blaze smoke
!
IF (LBLAZE) THEN
  IF (OFIRE) THEN
    CGETSVT(NSV_FIREBEG:NSV_FIREEND)='READ'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO BLAZE SCALAR VARIABLES IN INITIAL FMFILE",/,&
                       & "THE VARIABLES HAVE BEEN INITIALIZED TO ZERO")')
    CGETSVT(NSV_FIREBEG:NSV_FIREEND)='INIT'
  END IF
END IF
!
! Conditional sampling case
!
IF (LCONDSAMP) THEN
  IF (OCONDSAMP) THEN
    CGETSVT(NSV_CSBEG:NSV_CSEND)='READ'
!   IF(CCONF=='START') CGETSVT(NSV_CSBEG:NSV_CSEND)='INIT'       
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO PASSIVE SCALAR VARIABLES IN INITIAL FMFILE",/,&
                       & "THE VARIABLES HAVE BEEN INITIALIZED TO ZERO")')
    CGETSVT(NSV_CSBEG:NSV_CSEND)='INIT'
  END IF
END IF
!
! Blowing snow scheme
!
IF (LBLOWSNOW) THEN
  IF (OBLOWSNOW) THEN
    CGETSVT(NSV_SNWBEG:NSV_SNWEND)='READ'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='("THERE IS NO SCALAR VARIABLES FOR BLOWING SNOW &
         &SCHEME IN INITIAL FMFILE",/,&
         & "THE BLOWING SNOW VARIABLES HAVE BEEN INITIALIZED TO ZERO ")')
    CGETSVT(NSV_SNWBEG:NSV_SNWEND)='INIT'
  END IF
END IF
!
!
!
!*       3.5  Check coherence between the radiation control parameters
!
IF( CRAD == 'ECMW' .AND. CPROGRAM=='MESONH' ) THEN
  IF(CLW == 'RRTM' .AND. COPILW == 'SMSH') THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT=*) 'the SMSH parametrisation of LW optical properties for cloud ice'
    WRITE(UNIT=ILUOUT,FMT=*) '(COPILW) can not be used with RRTM radiation scheme'
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  ENDIF
  IF(CLW == 'MORC' .AND. COPWLW == 'LILI') THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT=*) 'the LILI parametrisation of LW optical properties for cloud water'
    WRITE(UNIT=ILUOUT,FMT=*) '(COPWLW) can not be used with MORC radiation scheme'
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  ENDIF
  IF( .NOT. LSUBG_COND) THEN
    WRITE(UNIT=ILUOUT,FMT=9000) KMI
    WRITE(UNIT=ILUOUT,FMT=*) 'YOU DO NOT WANT TO USE SUBGRID CONDENSATION'             
    WRITE(UNIT=ILUOUT,FMT=*) 'THE OVERLAP OPTION IS NOVLP=5 IN ini_radconf.f90'
  ELSE IF (CLW == 'MORC') THEN
    WRITE(UNIT=ILUOUT,FMT=9000) KMI
    WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE MORCRETTE LW SCHEME'                   
    WRITE(UNIT=ILUOUT,FMT=*) 'THE OVERLAP OPTION IS NOVLP=5 IN ini_radconf.f90'
  ELSE
    WRITE(UNIT=ILUOUT,FMT=9000) KMI
    WRITE(UNIT=ILUOUT,FMT=*) 'THE OVERLAP OPTION IS NOVLP=6 IN ini_radconf.f90'
  ENDIF
!
  IF( LCLEAR_SKY .AND. XDTRAD_CLONLY /= XDTRAD) THEN
    ! Check the validity of the LCLEAR_SKY approximation
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT=*) 'YOU WANT TO USE BOTH THE CLEAR-SKY APPROXIMATION'
    WRITE(UNIT=ILUOUT,FMT=*) '(i.e. AVERAGE THE WHOLE CLOUDFREE VERTICALS BUT KEEP'
    WRITE(UNIT=ILUOUT,FMT=*) 'ALL THE CLOUDY VERTICALS) AND'
    WRITE(UNIT=ILUOUT,FMT=*) 'THE CLOUD-ONLY APPROXIMATION (i.e. YOU CALL MORE OFTEN THE'
    WRITE(UNIT=ILUOUT,FMT=*) 'RADIATIONS FOR THE CLOUDY VERTICALS THAN FOR CLOUDFREE ONES).'
    WRITE(UNIT=ILUOUT,FMT=*) 'THIS IS NOT POSSIBLE, SO CHOOSE BETWEEN :'
    WRITE(UNIT=ILUOUT,FMT=*) 'XDTRAD_CLONLY = XDTRAD and LCLEAR_SKY = FALSE'
!
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
!
  IF( XDTRAD_CLONLY > XDTRAD ) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT='("BAD USE OF THE CLOUD-ONLY APPROXIMATION   " ,&
    &" XDTRAD SHOULD BE LARGER THAN XDTRAD_CLONLY                      ")')
!
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
!
  IF(( XDTRAD < XTSTEP ).OR. ( XDTRAD_CLONLY < XTSTEP )) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT='("THE RADIATION CALL XDTRAD OR XDTRAD_CLONLY " ,&
    &" IS MORE FREQUENT THAN THE TIME STEP SO ADJUST XDTRAD OR XDTRAD_CLONLY ")')
!
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
END IF
!
IF ( CRAD /= 'NONE' .AND. CPROGRAM=='MESONH' ) THEN
  CGETRAD='READ'
  IF( HRAD == 'NONE' .AND. CCONF=='RESTA') THEN
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT=*) 'YOU ARE PERFORMING A RESTART. FOR THIS SEGMENT, YOU ARE USING A RADIATION'
    WRITE(UNIT=ILUOUT,FMT=*) 'SCHEME AND NO RADIATION SCHEME WAS USED FOR THE PREVIOUS SEGMENT.'
    CGETRAD='INIT'
  END IF
  IF(CCONF=='START') THEN
    CGETRAD='INIT'
  END IF
  IF(CCONF=='RESTA' .AND. (.NOT. LAERO_FT) .AND.  (.NOT. LORILAM) &
                    .AND. (.NOT. LSALT)  .AND. (.NOT. LDUST)) THEN
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT=*) '!!! WARNING !!! FOR REPRODUCTIBILITY BETWEEN START and START+RESTART,'
    WRITE(UNIT=ILUOUT,FMT=*) 'YOU MUST USE LAERO_FT=T WITH CAER=TEGE IF CCONF=RESTA IN ALL SEGMENTS'
    WRITE(UNIT=ILUOUT,FMT=*) 'TO UPDATE THE OZONE AND AEROSOLS CLIMATOLOGY USED BY THE RADIATION CODE;'
  END IF
END IF
!
!        3.6  check the initialization of the deep convection scheme
!
IF ( (CDCONV /= 'KAFR') .AND. &
      (CSCONV /= 'KAFR') .AND. LCHTRANS ) THEN
   WRITE(UNIT=ILUOUT,FMT=9003) KMI
   WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE LCHTRANS OPTION= ",&
    &"CONVECTIVE TRANSPORT OF TRACERS BUT  IT CAN ONLY",&
    &"BE USED FOR THE KAIN FRITSCH SCHEME ")')
 !callabortstop
   CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
END IF
!
SELECT CASE ( CDCONV )
  CASE( 'KAFR' )
    IF (.NOT. ( LUSERV ) ) THEN
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE KAIN-FRITSCH DEEP CONV. ",&
      &" SCHEME. YOU MUST HAVE VAPOR ",/,"LUSERV IS SET TO TRUE ")')
      LUSERV=.TRUE.
    ELSE IF (.NOT. ( LUSERI ) ) THEN
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE KAIN-FRITSCH",&
      &" DEEP CONV. SCHEME. BUT THE DETRAINED CLOUD ICE WILL BE ADDED TO   ",&
      &" THE CLOUD WATER  ")')
    ELSE IF (.NOT. ( LUSERI.AND.LUSERC ) ) THEN
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE KAIN-FRITSCH",&
      &" DEEP CONV. SCHEME. BUT THE DETRAINED CLOUD WATER AND CLOUD ICE    ",&
      &" WILL BE ADDED TO THE WATER VAPOR FIELD  ")')
    END IF
    IF ( LCHTRANS .AND. NSV == 0 ) THEN
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE LCHTRANS OPTION= ",&
      &"CONVECTIVE TRANSPORT OF TRACERS BUT YOUR TRACER ",&
      &"NUMBER NSV IS ZERO ",/,"LCHTRANS IS SET TO FALSE")')
      LCHTRANS=.FALSE.
    END IF
END SELECT
!
IF ( CDCONV == 'KAFR' .AND. LCHTRANS .AND. NSV > 0 ) THEN
  IF( OCHTRANS ) THEN
    CGETSVCONV='READ'
  ELSE  
    CGETSVCONV='INIT'
  END IF
END IF
!
SELECT CASE ( CSCONV )
  CASE( 'KAFR' )
    IF (.NOT. ( LUSERV ) ) THEN
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE KAIN-FRITSCH SHALLOW CONV. ",&
      &" SCHEME. YOU MUST HAVE VAPOR ",/,"LUSERV IS SET TO TRUE ")')
      LUSERV=.TRUE.
    ELSE IF (.NOT. ( LUSERI ) ) THEN
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE KAIN-FRITSCH",&
      &" SHALLOW CONV. SCHEME. BUT THE DETRAINED CLOUD ICE WILL BE ADDED TO   ",&
      &" THE CLOUD WATER  ")')
    ELSE IF (.NOT. ( LUSERI.AND.LUSERC ) ) THEN
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE KAIN-FRITSCH",&
      &" SHALLOW CONV. SCHEME. BUT THE DETRAINED CLOUD WATER AND CLOUD ICE    ",&
      &" WILL BE ADDED TO THE WATER VAPOR FIELD  ")')
    END IF
    IF ( LCHTRANS .AND. NSV == 0 ) THEN
      WRITE(UNIT=ILUOUT,FMT=9002) KMI
      WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE LCHTRANS OPTION= ",&
      &"CONVECTIVE TRANSPORT OF TRACERS BUT YOUR TRACER ",&
      &"NUMBER NSV IS ZERO ",/,"LCHTRANS IS SET TO FALSE")')
      LCHTRANS=.FALSE.
    END IF
 CASE( 'EDKF' )
    IF (CTURB == 'NONE' ) THEN
      WRITE(UNIT=ILUOUT,FMT=9003) KMI
      WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE THE EDKF ", & 
      &"SHALLOW CONVECTION WITHOUT TURBULENCE SCHEME : ", &
      &"IT IS NOT POSSIBLE")')
!
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
    END IF
END SELECT
!
!
CGETCONV = 'SKIP'
!
IF ( (CDCONV /= 'NONE' .OR. CSCONV == 'KAFR' ) .AND. CPROGRAM=='MESONH') THEN
  CGETCONV = 'READ'
  IF( HDCONV == 'NONE' .AND. CCONF=='RESTA') THEN
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(UNIT=ILUOUT,FMT='(" YOU ARE PERFORMING A RESTART. FOR THIS  ",&
     &" SEGMENT, YOU ARE USING A DEEP CONVECTION SCHEME AND NO DEEP    ",&
     &" CONVECTION SCHEME WAS USED FOR THE PREVIOUS SEGMENT. ")')
!
    CGETCONV = 'INIT'
  END IF
  IF(CCONF=='START') THEN
    CGETCONV = 'INIT'
  END IF
END IF
!
!*       3.7  configuration and model version
!
IF (KMI == 1) THEN                                                
! 
  IF (L1D.AND.(CLBCX(1)/='CYCL'.AND.CLBCX(2)/='CYCL' &
          .AND.CLBCY(1)/='CYCL'.AND.CLBCY(2)/='CYCL')) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE A 1D MODEL VERSION WITH NON-CYCL",&
                & "CLBCX OR CLBCY VALUES")')
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
  IF (L2D.AND.(CLBCY(1)/='CYCL'.AND.CLBCY(2)/='CYCL')) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE A 2D MODEL VERSION WITH NON-CYCL",&
                & " CLBCY VALUES")')
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
  !
 IF ( (.NOT. LCARTESIAN) .AND. ( LCORIO) .AND. (.NOT. LGEOST_UV_FRC) ) THEN
   WRITE(UNIT=ILUOUT,FMT=9002) KMI
   WRITE(UNIT=ILUOUT,FMT='("BE CAREFUL YOU COULD HAVE SPURIOUS MOTIONS      " ,&
        & " NEAR THE LBC AS LCORIO=T and  LGEOST_UV_FRC=F")')
 END IF
  !
  IF ((.NOT.LFLAT).AND.OFLAT) THEN                                      
    WRITE(UNIT=ILUOUT,FMT=9002) KMI 
    WRITE(UNIT=ILUOUT,FMT=*) 'ZERO OROGRAPHY IN INITIAL FILE'
    WRITE(UNIT=ILUOUT,FMT=*) '***** ALL TERMS HAVE BEEN NEVERTHELESS COMPUTED WITHOUT SIMPLIFICATION*****'
    WRITE(UNIT=ILUOUT,FMT=*) 'THIS SHOULD LEAD TO ERRORS IN THE PRESSURE COMPUTATION'     
  END IF
  IF (LFLAT.AND.(.NOT.OFLAT)) THEN
    WRITE(UNIT=ILUOUT,FMT=9002) KMI
    WRITE(UNIT=ILUOUT,FMT='(" OROGRAPHY IS NOT EQUAL TO ZERO ",            &
          & "IN INITIAL FILE" ,/,                                          &
          & "******* OROGRAPHY HAS BEEN SET TO ZERO *********",/,          &
          & "ACCORDING TO ZERO OROGRAPHY, SIMPLIFICATIONS  HAVE  ",        &
          & "BEEN MADE IN  COMPUTATIONS")')
  END IF 
END IF
!
!*       3.8  System of equations
!
IF ( HEQNSYS /= CEQNSYS ) THEN
  WRITE(UNIT=ILUOUT,FMT=9003) KMI
  WRITE(ILUOUT,FMT=*) 'YOU HAVE CHANGED THE SYSTEM OF EQUATIONS'
  WRITE(ILUOUT,FMT=*) 'THE ANELASTIC CONSTRAINT IS PERHAPS CHANGED :'
  WRITE(ILUOUT,FMT=*) 'FOR THE INITIAL FILE YOU HAVE USED ',HEQNSYS
  WRITE(ILUOUT,FMT=*) 'FOR THE RUN YOU PLAN TO USE ',CEQNSYS
  WRITE(ILUOUT,FMT=*) 'THIS CAN LEAD TO A NUMERICAL EXPLOSION IN THE FIRST TIME STEPS'
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
END IF
!
!        3.9  Numerical schemes
!
IF ( (CUVW_ADV_SCHEME == 'CEN4TH') .AND. &
      (CTEMP_SCHEME /= 'LEFR') .AND. (CTEMP_SCHEME /= 'RKC4') ) THEN
   WRITE(UNIT=ILUOUT,FMT=9003) KMI
   WRITE(UNIT=ILUOUT,FMT='("CEN4TH SCHEME HAS TO BE USED WITH ",&
    &"CTEMP_SCHEME = LEFR of RKC4 ONLY")')   
 !callabortstop
   CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
END IF
!
IF ( (CUVW_ADV_SCHEME == 'WENO_K') .AND. LNUMDIFU ) THEN
   WRITE(UNIT=ILUOUT,FMT=9002) KMI
   WRITE(UNIT=ILUOUT,FMT='("YOU WANT TO USE NUMERICAL DIFFUSION ",&
    &"WITH WENO SCHEME ALREADY DIFFUSIVE")')   
END IF
!-------------------------------------------------------------------------------
!
!*       4.    CHECK COHERENCE BETWEEN EXSEG VARIABLES
!              ---------------------------------------
!        
!*       4.1  coherence between coupling variables in EXSEG file  
!                      
IF (KMI == 1) THEN
  NCPL_NBR = 0
  DO JCI = 1,JPCPLFILEMAX
    IF (LEN_TRIM(CCPLFILE(JCI)) /= 0) THEN        ! Finds the number 
      NCPL_NBR = NCPL_NBR + 1                     ! of coupling files
    ENDIF
    IF (JCI/=JPCPLFILEMAX) THEN                   ! Deplaces the coupling files
      IF ((LEN_TRIM(CCPLFILE(JCI)) == 0) .AND.   &! names if one missing
          (LEN_TRIM(CCPLFILE(JCI+1)) /= 0)) THEN
        DO JI=JCI,JPCPLFILEMAX-1
          CCPLFILE(JI)=CCPLFILE(JI+1)
        END DO
        CCPLFILE(JPCPLFILEMAX)='    '
      END IF
    END IF
  END DO
!
  IF (NCPL_NBR /= 0) THEN         
    LSTEADYLS = .FALSE.
  ELSE
    LSTEADYLS = .TRUE.
  ENDIF 
END IF
!        
!*       4.3   check consistency in forcing switches
!
IF ( LFORCING ) THEN
  IF ( LRELAX_THRV_FRC .AND. ( LTEND_THRV_FRC .OR. LGEOST_TH_FRC ) ) THEN
    WRITE(UNIT=ILUOUT,FMT=9002) KMI
    WRITE(ILUOUT,FMT=*) 'YOU CHOSE A TEMPERATURE AND HUMIDITY RELAXATION'
    WRITE(ILUOUT,FMT=*) 'TOGETHER WITH TENDENCY OR GEOSTROPHIC FORCING'
    WRITE(ILUOUT,FMT=*) &
  'YOU MIGHT CHECK YOUR SWITCHES: LRELAX_THRV_FRC, LTEND_THRV_FRC, AND'
    WRITE(ILUOUT,FMT=*) 'LGEOST_TH_FRC'
  END IF
!
  IF ( LRELAX_UV_FRC .AND. LRELAX_UVMEAN_FRC) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(ILUOUT,FMT=*) 'YOU MUST CHOOSE BETWEEN A RELAXATION APPLIED TO'
    WRITE(ILUOUT,FMT=*) 'THE 3D FULL WIND FIELD (LRELAX_UV_FRC) OR'
    WRITE(ILUOUT,FMT=*) 'THE HORIZONTAL MEAN WIND (LRELAX_UVMEAN_FRC)'
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
!
  IF ( (LRELAX_UV_FRC .OR. LRELAX_UVMEAN_FRC) .AND. LGEOST_UV_FRC ) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(ILUOUT,FMT=*) 'YOU MUST NOT USE A WIND RELAXATION' 
    WRITE(ILUOUT,FMT=*) 'TOGETHER WITH A GEOSTROPHIC FORCING'
    WRITE(ILUOUT,FMT=*) 'CHECK SWITCHES: LRELAX_UV_FRC, LRELAX_UVMEAN_FRC, LGEOST_UV_FRC'
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
!
  IF ( CRELAX_HEIGHT_TYPE.NE."FIXE" .AND. CRELAX_HEIGHT_TYPE.NE."THGR" ) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(ILUOUT,FMT=*) 'CRELAX_HEIGHT_TYPE MUST BE EITHER "FIXE" OR "THGR"'
    WRITE(ILUOUT,FMT=*) 'BUT IT IS "', CRELAX_HEIGHT_TYPE, '"'
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
!
  IF ( .NOT.LCORIO .AND. LGEOST_UV_FRC ) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(ILUOUT,FMT=*) 'YOU CANNOT HAVE A GEOSTROPHIC FORCING WITHOUT'
    WRITE(ILUOUT,FMT=*) 'ACTIVATING LCORIOLIS OPTION'
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
!
  IF ( LPGROUND_FRC ) THEN
    WRITE(ILUOUT,FMT=*) 'SURFACE PRESSURE FORCING NOT YET IMPLEMENTED'
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
!
END IF
!
IF (LTRANS .AND. .NOT. LFLAT ) THEN
    WRITE(UNIT=ILUOUT,FMT=9003) KMI
    WRITE(ILUOUT,FMT=*) 'YOU ASK FOR A CONSTANT SPEED DOMAIN TRANSLATION '
    WRITE(ILUOUT,FMT=*) 'BUT NOT IN THE FLAT TERRAIN CASE:'
    WRITE(ILUOUT,FMT=*) 'THIS IS NOT ALLOWED ACTUALLY'
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
END IF
!
!*       4.4  Check the coherence between the LUSERn and LHORELAX
!
IF (.NOT. LUSERV .AND. LHORELAX_RV) THEN
  LHORELAX_RV=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX RV FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RV=FALSE'
END IF
!
IF (.NOT. LUSERC .AND. LHORELAX_RC) THEN
  LHORELAX_RC=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX RC FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RC=FALSE'
END IF
!
IF (.NOT. LUSERR .AND. LHORELAX_RR) THEN
  LHORELAX_RR=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX RR FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RR=FALSE'
END IF
!
IF (.NOT. LUSERI .AND. LHORELAX_RI) THEN
  LHORELAX_RI=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX RI FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RI=FALSE'
END IF
!
IF (.NOT. LUSERS .AND. LHORELAX_RS) THEN
  LHORELAX_RS=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX RS FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RS=FALSE'
END IF
!
IF (.NOT. LUSERG .AND. LHORELAX_RG) THEN
  LHORELAX_RG=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX RG FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RG=FALSE'
END IF
!
IF (.NOT. LUSERH .AND. LHORELAX_RH) THEN
  LHORELAX_RH=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX RH FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RH=FALSE'
END IF
!
IF (CTURB=='NONE' .AND. LHORELAX_TKE) THEN
  LHORELAX_TKE=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX TKE FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_TKE=FALSE'
END IF
!
!
IF (CCLOUD/='C2R2'  .AND. CCLOUD/='KHKO'  .AND. LHORELAX_SVC2R2) THEN
  LHORELAX_SVC2R2=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX C2R2 or KHKO FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVC2R2=FALSE'
END IF
!
IF (CCLOUD/='C3R5' .AND. LHORELAX_SVC1R3) THEN
  LHORELAX_SVC1R3=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX C3R5 FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVC1R3=FALSE'
END IF
!
IF (CCLOUD/='LIMA' .AND. LHORELAX_SVLIMA) THEN
  LHORELAX_SVLIMA=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX LIMA FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVLIMA=FALSE'
END IF
!
IF (CELEC(1:3) /= 'ELE' .AND. LHORELAX_SVELEC) THEN
  LHORELAX_SVELEC=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX ELEC FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVELEC=FALSE'
END IF
!
IF (.NOT. LUSECHEM .AND. LHORELAX_SVCHEM) THEN
  LHORELAX_SVCHEM=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX CHEM FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVCHEM=FALSE'
END IF
!
IF (.NOT. LUSECHIC .AND. LHORELAX_SVCHIC) THEN
  LHORELAX_SVCHIC=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX ICE CHEM FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVCHIC=FALSE'
END IF
!
IF (.NOT. LORILAM .AND. LHORELAX_SVAER) THEN
  LHORELAX_SVAER=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX AEROSOL FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVAER=FALSE'
END IF

IF (.NOT. LDUST .AND. LHORELAX_SVDST) THEN
  LHORELAX_SVDST=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX DUST FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVDST=FALSE'
END IF

IF (.NOT. LSALT .AND. LHORELAX_SVSLT) THEN
  LHORELAX_SVSLT=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX SEA SALT FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVSLT=FALSE'
END IF

IF (.NOT. LPASPOL .AND. LHORELAX_SVPP) THEN
  LHORELAX_SVPP=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX PASSIVE POLLUTANT FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVPP=FALSE'
END IF
#ifdef MNH_FOREFIRE
IF (.NOT. LFOREFIRE .AND. LHORELAX_SVFF) THEN
  LHORELAX_SVFF=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX FOREFIRE FLUXES BUT THEY DO NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVFF=FALSE'
END IF
#endif
IF (.NOT. LBLAZE .AND. LHORELAX_SVFIRE) THEN
  LHORELAX_SVFIRE=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX BLAZE FLUXES BUT THEY DO NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVFIRE=FALSE'
END IF
IF (.NOT. LCONDSAMP .AND. LHORELAX_SVCS) THEN
  LHORELAX_SVCS=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX CONDITIONAL SAMPLING FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVCS=FALSE'
END IF

IF (.NOT. LBLOWSNOW .AND. LHORELAX_SVSNW) THEN
  LHORELAX_SVSNW=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX BLOWING SNOW FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SVSNW=FALSE'
END IF

IF (ANY(LHORELAX_SV(NSV+1:))) THEN
  LHORELAX_SV(NSV+1:)=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX SV(NSV+1:) FIELD BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SV(NSV+1:)=FALSE'
END IF
!
!*       4.5   check the number of points for the horizontal relaxation
!
IF ( NRIMX > KRIMX .AND. .NOT.LHORELAX_SVELEC ) THEN
  NRIMX = KRIMX 
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO USE A LARGER NUMBER OF POINTS '
  WRITE(ILUOUT,FMT=*) 'FOR THE HORIZONTAL RELAXATION THAN THE '
  WRITE(ILUOUT,FMT=*) 'CORRESPONDING NUMBER OF LARGE SCALE FIELDS:'
  WRITE(ILUOUT,FMT=*) 'IT IS THEREFORE REDUCED TO NRIMX =',NRIMX
END IF
!
IF ( L2D .AND. KRIMY>0 ) THEN
  NRIMY = 0 
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO USE A 2D MODEL THEREFORE  NRIMY=0 '
END IF
!
IF ( NRIMY > KRIMY .AND. .NOT.LHORELAX_SVELEC ) THEN
  NRIMY = KRIMY 
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO USE A LARGER NUMBER OF POINTS '
  WRITE(ILUOUT,FMT=*) 'FOR THE HORIZONTAL RELAXATION THAN THE '
  WRITE(ILUOUT,FMT=*) 'CORRESPONDING NUMBER OF LARGE SCALE FIELDS:'
  WRITE(ILUOUT,FMT=*) 'IT IS THEREFORE REDUCED TO NRIMY =',NRIMY
END IF
!
IF ( (.NOT. LHORELAX_UVWTH) .AND. (.NOT.(ANY(LHORELAX_SV))) .AND.  &
     (.NOT. LHORELAX_SVC2R2).AND. (.NOT. LHORELAX_SVC1R3)   .AND.  &
     (.NOT. LHORELAX_SVLIMA).AND.                                  &
     (.NOT. LHORELAX_SVELEC).AND. (.NOT. LHORELAX_SVCHEM)   .AND.  &
     (.NOT. LHORELAX_SVLG)  .AND. (.NOT. LHORELAX_SVPP)     .AND.  &
     (.NOT. LHORELAX_SVCS)  .AND. (.NOT. LHORELAX_SVFIRE)   .AND.  &
#ifdef MNH_FOREFIRE
     (.NOT. LHORELAX_SVFF)  .AND.                                  &
#endif
     (.NOT. LHORELAX_RV)    .AND. (.NOT. LHORELAX_RC)       .AND.  &
     (.NOT. LHORELAX_RR)    .AND. (.NOT. LHORELAX_RI)       .AND.  &
     (.NOT. LHORELAX_RS)    .AND. (.NOT. LHORELAX_RG)       .AND.  &
     (.NOT. LHORELAX_RH)    .AND. (.NOT. LHORELAX_TKE)      .AND.  &
     (.NOT. LHORELAX_SVCHIC).AND.                                  &
                                  (NRIMX /= 0 .OR. NRIMY /= 0))  THEN  
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU DO NOT WANT TO USE THE HORIZONTAL RELAXATION '
  WRITE(ILUOUT,FMT=*) 'THEREFORE NRIMX=NRIMY=0     '
  NRIMX=0
  NRIMY=0
END IF
!
IF ((LHORELAX_UVWTH  .OR. LHORELAX_SVPP   .OR.  &
     LHORELAX_SVCS   .OR. LHORELAX_SVFIRE .OR.  &
#ifdef MNH_FOREFIRE
     LHORELAX_SVFF   .OR.                       &
#endif
     LHORELAX_SVC2R2 .OR. LHORELAX_SVC1R3 .OR.  &
     LHORELAX_SVLIMA .OR.                       &
     LHORELAX_SVELEC .OR. LHORELAX_SVCHEM .OR.  &
     LHORELAX_SVLG   .OR. ANY(LHORELAX_SV) .OR. &
     LHORELAX_RV     .OR. LHORELAX_RC .OR.      &
     LHORELAX_RR     .OR. LHORELAX_RI .OR.      &
     LHORELAX_RG     .OR. LHORELAX_RS .OR.      &
     LHORELAX_RH     .OR. LHORELAX_TKE.OR.      &
     LHORELAX_SVCHIC )                          &
     .AND. (NRIMX==0 .OR. (NRIMY==0 .AND. .NOT.(L2D) ))) THEN
  WRITE(UNIT=ILUOUT,FMT=9003) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO USE THE HORIZONTAL RELAXATION '
  WRITE(ILUOUT,FMT=*) 'BUT NRIMX OR NRIMY=0 CHANGE YOUR VALUES   '
  WRITE(ILUOUT,FMT=*) "LHORELAX_UVWTH=",LHORELAX_UVWTH
  WRITE(ILUOUT,FMT=*) "LHORELAX_SVC2R2=",LHORELAX_SVC2R2
  WRITE(ILUOUT,FMT=*) "LHORELAX_SVC1R3=",LHORELAX_SVC1R3
  WRITE(ILUOUT,FMT=*) "LHORELAX_SVLIMA=",LHORELAX_SVLIMA
  WRITE(ILUOUT,FMT=*) "LHORELAX_SVELEC=",LHORELAX_SVELEC
  WRITE(ILUOUT,FMT=*) "LHORELAX_SVCHEM=",LHORELAX_SVCHEM
  WRITE(ILUOUT,FMT=*) "LHORELAX_SVCHIC=",LHORELAX_SVCHIC
  WRITE(ILUOUT,FMT=*) "LHORELAX_SVLG=",LHORELAX_SVLG
  WRITE(ILUOUT,FMT=*) "LHORELAX_SVPP=",LHORELAX_SVPP
  WRITE(ILUOUT,FMT=*) "LHORELAX_SVFIRE=",LHORELAX_SVFIRE
#ifdef MNH_FOREFIRE
  WRITE(ILUOUT,FMT=*) "LHORELAX_SVFF=",LHORELAX_SVFF
#endif
  WRITE(ILUOUT,FMT=*) "LHORELAX_SVCS=",LHORELAX_SVCS
  WRITE(ILUOUT,FMT=*) "LHORELAX_SV=",LHORELAX_SV
  WRITE(ILUOUT,FMT=*) "LHORELAX_RV=",LHORELAX_RV
  WRITE(ILUOUT,FMT=*) "LHORELAX_RC=",LHORELAX_RC
  WRITE(ILUOUT,FMT=*) "LHORELAX_RR=",LHORELAX_RR
  WRITE(ILUOUT,FMT=*) "LHORELAX_RI=",LHORELAX_RI
  WRITE(ILUOUT,FMT=*) "LHORELAX_RG=",LHORELAX_RG
  WRITE(ILUOUT,FMT=*) "LHORELAX_RS=",LHORELAX_RS
  WRITE(ILUOUT,FMT=*) "LHORELAX_RH=",LHORELAX_RH
  WRITE(ILUOUT,FMT=*) "LHORELAX_TKE=", LHORELAX_TKE
  WRITE(ILUOUT,FMT=*) "NRIMX=",NRIMX
  WRITE(ILUOUT,FMT=*) "NRIMY=",NRIMY
  WRITE(ILUOUT,FMT=*) "L2D=",L2D
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
END IF
! 
IF ((LHORELAX_UVWTH  .OR. LHORELAX_SVPP  .OR.   &
     LHORELAX_SVCS   .OR. LHORELAX_SVFIRE .OR.  &
#ifdef MNH_FOREFIRE
     LHORELAX_SVFF   .OR.                       &
#endif
     LHORELAX_SVC2R2 .OR. LHORELAX_SVC1R3 .OR.  &
     LHORELAX_SVLIMA .OR.                       &
     LHORELAX_SVELEC .OR. LHORELAX_SVCHEM .OR.  &
     LHORELAX_SVLG   .OR. ANY(LHORELAX_SV) .OR. &
     LHORELAX_RV     .OR. LHORELAX_RC .OR.      &
     LHORELAX_RR     .OR. LHORELAX_RI .OR.      &
     LHORELAX_RG     .OR. LHORELAX_RS .OR.      &
     LHORELAX_RH     .OR. LHORELAX_TKE.OR.      &
     LHORELAX_SVCHIC )                          &
     .AND. (KMI /=1)) THEN
  WRITE(UNIT=ILUOUT,FMT=9003) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO USE THE HORIZONTAL RELAXATION '
  WRITE(ILUOUT,FMT=*) 'FOR A NESTED MODEL BUT THE COUPLING IS ALREADY DONE' 
  WRITE(ILUOUT,FMT=*) 'BY THE GRID NESTING. CHANGE LHORELAX TO FALSE'
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
END IF
!
IF ((LHORELAX_UVWTH  .OR. LHORELAX_SVPP  .OR.   &
     LHORELAX_SVCS   .OR. LHORELAX_SVFIRE .OR.  &
#ifdef MNH_FOREFIRE
     LHORELAX_SVFF   .OR.                       &
#endif
     LHORELAX_SVC2R2 .OR. LHORELAX_SVC1R3 .OR.  &
     LHORELAX_SVLIMA .OR.                       &
     LHORELAX_SVELEC .OR. LHORELAX_SVCHEM .OR.  &
     LHORELAX_SVLG   .OR. ANY(LHORELAX_SV) .OR. &
     LHORELAX_RV     .OR. LHORELAX_RC .OR.      &
     LHORELAX_RR     .OR. LHORELAX_RI .OR.      &
     LHORELAX_RG     .OR. LHORELAX_RS .OR.      &
     LHORELAX_RH     .OR. LHORELAX_TKE.OR.      &
     LHORELAX_SVCHIC )                          &
     .AND. (CLBCX(1)=='CYCL'.OR.CLBCX(2)=='CYCL' &
          .OR.CLBCY(1)=='CYCL'.OR.CLBCY(2)=='CYCL')) THEN
  WRITE(UNIT=ILUOUT,FMT=9003) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO USE THE HORIZONTAL RELAXATION '
  WRITE(ILUOUT,FMT=*) 'FOR CYCLIC CLBCX OR CLBCY VALUES' 
  WRITE(ILUOUT,FMT=*) 'CHANGE LHORELAX TO FALSE'
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
END IF
!
IF (KMI==1) THEN
  GRELAX = .NOT.(OUSERV)  .AND.  LUSERV  .AND. LHORELAX_RV
ELSE
  GRELAX = .NOT.(LUSERV_G(NDAD(KMI)))  .AND.  LUSERV_G(KMI).AND. LHORELAX_RV 
END IF 
!
IF ( GRELAX )  THEN
  LHORELAX_RV=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX THE RV FIELD'
  WRITE(ILUOUT,FMT=*) 'TOWARDS THE LARGE SCALE FIELD OF MODEL',NDAD(KMI)
  WRITE(ILUOUT,FMT=*) 'BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RV=FALSE'
END IF
!
IF (KMI==1) THEN
  GRELAX = .NOT.(OUSERC)  .AND.  LUSERC  .AND. LHORELAX_RC
ELSE
  GRELAX = .NOT.(LUSERC_G(NDAD(KMI)))  .AND.  LUSERC_G(KMI).AND. LHORELAX_RC 
END IF 
!
IF ( GRELAX )  THEN
  LHORELAX_RC=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX THE RC FIELD'
  WRITE(ILUOUT,FMT=*) 'TOWARDS THE LARGE SCALE FIELD OF MODEL',NDAD(KMI)
  WRITE(ILUOUT,FMT=*) 'BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RC=FALSE'
END IF
!
IF (KMI==1) THEN
  GRELAX = .NOT.(OUSERR)  .AND.  LUSERR  .AND. LHORELAX_RR
ELSE
  GRELAX = .NOT.(LUSERR_G(NDAD(KMI)))  .AND.  LUSERR_G(KMI).AND. LHORELAX_RR 
END IF 
!
IF ( GRELAX )  THEN
  LHORELAX_RR=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX THE RR FIELD'
  WRITE(ILUOUT,FMT=*) 'TOWARDS THE LARGE SCALE FIELD OF MODEL',NDAD(KMI)
  WRITE(ILUOUT,FMT=*) 'BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RR=FALSE'
END IF
!
IF (KMI==1) THEN
  GRELAX = .NOT.(OUSERI)  .AND.  LUSERI  .AND. LHORELAX_RI
ELSE
  GRELAX = .NOT.(LUSERI_G(NDAD(KMI)))  .AND.  LUSERI_G(KMI).AND. LHORELAX_RI 
END IF 
!
IF ( GRELAX )  THEN
  LHORELAX_RI=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX THE RI FIELD'
  WRITE(ILUOUT,FMT=*) 'TOWARDS THE LARGE SCALE FIELD OF MODEL',NDAD(KMI)
  WRITE(ILUOUT,FMT=*) 'BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RI=FALSE'
END IF
!
IF (KMI==1) THEN
  GRELAX = .NOT.(OUSERG)  .AND.  LUSERG  .AND. LHORELAX_RG
ELSE
  GRELAX = .NOT.(LUSERG_G(NDAD(KMI)))  .AND.  LUSERG_G(KMI).AND. LHORELAX_RG 
END IF 
!
IF ( GRELAX )  THEN
  LHORELAX_RG=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX THE RG FIELD'
  WRITE(ILUOUT,FMT=*) 'TOWARDS THE LARGE SCALE FIELD OF MODEL',NDAD(KMI)
  WRITE(ILUOUT,FMT=*) 'BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RG=FALSE'
END IF
!
IF (KMI==1) THEN
  GRELAX = .NOT.(OUSERH)  .AND.  LUSERH  .AND. LHORELAX_RH
ELSE
  GRELAX = .NOT.(LUSERH_G(NDAD(KMI)))  .AND.  LUSERH_G(KMI).AND. LHORELAX_RH 
END IF 
!
IF ( GRELAX )  THEN
  LHORELAX_RH=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX THE RH FIELD'
  WRITE(ILUOUT,FMT=*) 'TOWARDS THE LARGE SCALE FIELD OF MODEL',NDAD(KMI)
  WRITE(ILUOUT,FMT=*) 'BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RH=FALSE'
END IF
!
IF (KMI==1) THEN
  GRELAX = .NOT.(OUSERS)  .AND.  LUSERS  .AND. LHORELAX_RS
ELSE
  GRELAX = .NOT.(LUSERS_G(NDAD(KMI)))  .AND.  LUSERS_G(KMI).AND. LHORELAX_RS 
END IF 
!
IF ( GRELAX )  THEN
  LHORELAX_RS=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX THE RS FIELD'
  WRITE(ILUOUT,FMT=*) 'TOWARDS THE LARGE SCALE FIELD OF MODEL',NDAD(KMI)
  WRITE(ILUOUT,FMT=*) 'BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_RS=FALSE'
END IF
!
IF (KMI==1) THEN
  GRELAX = HTURB=='NONE'  .AND.  LUSETKE(1).AND. LHORELAX_TKE
ELSE
  GRELAX = .NOT.(LUSETKE(NDAD(KMI)))  .AND.  LUSETKE(KMI) .AND. LHORELAX_TKE 
END IF 
!
IF ( GRELAX )  THEN
  LHORELAX_TKE=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX THE TKE FIELD'
  WRITE(ILUOUT,FMT=*) 'TOWARDS THE LARGE SCALE FIELD OF MODEL',NDAD(KMI)
  WRITE(ILUOUT,FMT=*) 'BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_TKE=FALSE'
END IF
!
!
DO JSV = 1,NSV_USER
!
  IF (KMI==1) THEN
    GRELAX = KSV_USER<JSV  .AND.  LUSESV(JSV,1).AND.  LHORELAX_SV(JSV) 
  ELSE
    GRELAX = .NOT.(LUSESV(JSV,NDAD(KMI))) .AND. LUSESV(JSV,KMI) .AND. LHORELAX_SV(JSV)   
  END IF 
  !
  IF ( GRELAX )  THEN
  LHORELAX_SV(JSV)=.FALSE.
  WRITE(UNIT=ILUOUT,FMT=9001) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO RELAX THE ',JSV,'  SV FIELD'
  WRITE(ILUOUT,FMT=*) 'TOWARDS THE LARGE SCALE FIELD OF MODEL',NDAD(KMI)
  WRITE(ILUOUT,FMT=*) 'BUT IT DOES NOT EXIST.'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LHORELAX_SV(',JSV,')=FALSE'
 END IF
END DO
!
!*       4.6   consistency in LES diagnostics choices
!
IF (CLES_NORM_TYPE=='EKMA' .AND. .NOT. LCORIO) THEN
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO USE THE EKMAN NORMALIZATION'
  WRITE(ILUOUT,FMT=*) 'BUT CORIOLIS FORCE IS NOT USED (LCORIO=.FALSE.)'
  WRITE(ILUOUT,FMT=*) 'THEN, NO NORMALIZATION IS PERFORMED'
  CLES_NORM_TYPE='NONE'
END IF
!
!*       4.7  Check the coherence with LNUMDIFF
!
IF (L1D .AND. (LNUMDIFU .OR. LNUMDIFTH .OR. LNUMDIFSV) ) THEN
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU WANT TO USE HORIZONTAL DIFFUSION '
  WRITE(ILUOUT,FMT=*) 'BUT YOU ARE IN A COLUMN MODEL (L1D=.TRUE.).'
  WRITE(ILUOUT,FMT=*) 'THEREFORE LNUMDIFU and LNUMDIFTH and LNUMDIFSV'
  WRITE(ILUOUT,FMT=*) 'ARE SET TO FALSE'
  LNUMDIFU=.FALSE.
  LNUMDIFTH=.FALSE.
  LNUMDIFSV=.FALSE.
END IF
!
IF (.NOT. LNUMDIFTH .AND. LZDIFFU) THEN
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(ILUOUT,FMT=*) 'YOU DO NOT WANT TO USE HORIZONTAL DIFFUSION (LNUMDIFTH=F)'
  WRITE(ILUOUT,FMT=*) 'BUT YOU WANT TO USE Z-NUMERICAL DIFFUSION  '
  WRITE(ILUOUT,FMT=*) 'THEREFORE LNUMDIFTH IS SET TO TRUE'
  LNUMDIFTH=.TRUE.
END IF
!
!*       4.8  Other
!
IF (XTNUDGING < 4.*XTSTEP) THEN
  XTNUDGING = 4.*XTSTEP
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(UNIT=ILUOUT,FMT='("TIME SCALE FOR NUDGING CAN NOT BE SMALLER THAN",  &
                 &   " FOUR TIMES THE TIME STEP")')
  WRITE(ILUOUT,FMT=*) 'XTNUDGING is SET TO ',XTNUDGING
END IF
!
!
IF (XWAY(KMI) == 3. ) THEN
  XWAY(KMI) = 2.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(UNIT=ILUOUT,FMT='("XWAY=3 DOES NOT EXIST ANYMORE; ", &                               
                 &   " IT IS REPLACED BY XWAY=2 ")')
END IF
!
IF ( (KMI == 1) .AND. XWAY(KMI) /= 0. ) THEN
  XWAY(KMI) = 0.
  WRITE(UNIT=ILUOUT,FMT=9002) KMI
  WRITE(UNIT=ILUOUT,FMT='("XWAY MUST BE EQUAL TO 0 FOR DAD MODEL")')   
END IF
!
!JUANZ ZRESI solver need BSPLITTING 
IF ( CPRESOPT == 'ZRESI' .AND. CSPLIT /= 'BSPLITTING' ) THEN
  WRITE(UNIT=ILUOUT,FMT=9003) KMI
  WRITE(UNIT=ILUOUT,FMT='("Paralleliez in Z solver CPRESOPT=ZRESI need also CSPLIT=BSPLITTING ")')
  WRITE(ILUOUT,FMT=*) ' ERROR you have to set also CSPLIT=BSPLITTING '
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
END IF
!
IF ( LEN_TRIM(HINIFILEPGD)>0 ) THEN
  IF ( CINIFILEPGD/=HINIFILEPGD ) THEN 
    WRITE(UNIT=ILUOUT,FMT=9001) KMI
    WRITE(ILUOUT,FMT=*) ' ERROR : in EXSEG1.nam, in NAM_LUNITn you have CINIFILEPGD= ',CINIFILEPGD
    WRITE(ILUOUT,FMT=*) ' whereas in .des you have CINIFILEPGD= ',HINIFILEPGD
    WRITE(ILUOUT,FMT=*) ' Please check your Namelist '
    WRITE(ILUOUT,FMT=*) ' For example, you may have specified the un-nested PGD file instead of the nested PGD file '
    WRITE(ILUOUT,FMT=*) 
    WRITE(ILUOUT,FMT=*) '###############'
    WRITE(ILUOUT,FMT=*) ' MESONH ABORTS'
    WRITE(ILUOUT,FMT=*) '###############'
    WRITE(ILUOUT,FMT=*) 
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_EXSEG_n','')
  END IF
ELSE
  CINIFILEPGD = ''
!* note that after a spawning, there is no value for CINIFILEPGD in the .des file,
!  so the checking cannot be made if the user starts a simulation directly from
!  a spawned file (without the prep_real_case stage)
END IF
!-------------------------------------------------------------------------------
!
!*       5.    WE DO NOT FORGET TO UPDATE ALL DOLLARN NAMELIST VARIABLES
!              ---------------------------------------------------------
!
CALL UPDATE_NAM_LUNITN
CALL UPDATE_NAM_CONFN
CALL UPDATE_NAM_DRAGTREEN
CALL UPDATE_NAM_DRAGBLDGN
CALL UPDATE_NAM_COUPLING_LEVELSN
CALL UPDATE_NAM_DYNN
CALL UPDATE_NAM_ADVN
CALL UPDATE_NAM_PARAMN
CALL UPDATE_NAM_PARAM_RADN
#ifdef MNH_ECRAD
CALL UPDATE_NAM_PARAM_ECRADN
#endif
CALL UPDATE_NAM_PARAM_KAFRN
CALL UPDATE_NAM_LBCN
CALL UPDATE_NAM_NUDGINGN
CALL UPDATE_NAM_BLANKN
CALL UPDATE_NAM_CH_MNHCN
CALL UPDATE_NAM_CH_SOLVERN
CALL UPDATE_NAM_SERIESN
CALL UPDATE_NAM_BLOWSNOWN
CALL UPDATE_NAM_PROFILERn
CALL UPDATE_NAM_STATIONn
CALL UPDATE_NAM_FIREn
!-------------------------------------------------------------------------------
WRITE(UNIT=ILUOUT,FMT='(/)')
!-------------------------------------------------------------------------------
!
!*       6.    FORMATS
!              -------
!
9000  FORMAT(/,'NOTE  IN READ_EXSEG FOR MODEL ', I2, ' : ',/, &
             '--------------------------------')
9001  FORMAT(/,'CAUTION ERROR IN READ_EXSEG FOR MODEL ', I2,' : ',/, &
             '----------------------------------------' )
9002  FORMAT(/,'WARNING IN READ_EXSEG FOR MODEL ', I2,' : ',/, &
             '----------------------------------' )
9003  FORMAT(/,'FATAL ERROR IN READ_EXSEG FOR MODEL ', I2,' : ',/, &
             '--------------------------------------' )
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_EXSEG_n
