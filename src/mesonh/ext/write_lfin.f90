!MNH_LIC Copyright 1994-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_WRITE_LFIFM_n
!     #########################
!
INTERFACE
!
SUBROUTINE WRITE_LFIFM_n(TPFILE,HDADFILE)
!
USE MODD_IO, ONLY: TFILEDATA
!
IMPLICIT NONE
!
TYPE(TFILEDATA), INTENT(INOUT) :: TPFILE   ! File characteristics
CHARACTER(LEN=*),INTENT(IN) :: HDADFILE ! Corresponding FM-file name of its DAD model
END SUBROUTINE WRITE_LFIFM_n
!
END INTERFACE
!
END MODULE MODI_WRITE_LFIFM_n
!
!
!     ##########################################
      SUBROUTINE WRITE_LFIFM_n(TPFILE,HDADFILE)
!     ##########################################
!
!!****  *WRITE_LFIFM_n* - routine to write a LFIFM file for model $n
!!
!!    PURPOSE
!!    -------
!        The purpose of this routine is to write an initial LFIFM File 
!     of name YFMFILE//'.lfi' with the FM routines.  
!
!!**  METHOD
!!    ------
!!      The data are written in the LFIFM file :
!!        - dimensions
!!        - grid variables
!!        - configuration variables
!!        - prognostic variables at time t and t-dt
!!        - 1D anelastic reference state
!!
!!      The localization on the model grid is also indicated :
!!
!!        IGRID = 1 for mass grid point
!!        IGRID = 2 for U grid point
!!        IGRID = 3 for V grid point
!!        IGRID = 4 for w grid point
!!        IGRID = 0 for meaningless case
!!          
!!
!!    EXTERNAL
!!    --------
!!      WRITE_BALLOON_n : routine to write balloon records
!!      WRITE_LB_n : routine to write LB fields
!!      FMWRIT     : FM-routine to write a record
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_DIM_n   : contains dimensions
!!      Module MODD_TIME    : contains time variables for all models
!!      Module MODD_TIME_n   : contains time variables 
!!      Module MODD_GRID    : contains spatial grid variables for all models
!!      Module MODD_GRID_n : contains spatial grid variables
!!      Module MODD_REF     : contains reference state variables
!!      Module MODD_LUNIT_n: contains logical unit variables.
!!      Module MODD_CONF    : contains configuration variables for all models
!!      Module MODD_CONF_n  : contains configuration variables
!!      Module MODD_FIELD_n  : contains prognostic variables
!!      Module MODD_GR_FIELD_n : contains surface prognostic variables
!!      Module MODD_LSFIELD_n  : contains Larger Scale variables
!!      Module MODD_PARAM_n    : contains parameterization options
!!      Module MODD_TURB_n    : contains turbulence options
!!      Module MODD_FRC    : contains forcing variables
!!      Module MODD_DEEP_CONVECTION_n : contains deep convection tendencies
!!      Module MODD_PARAM_KAFR_n : contains configuration
!!      Module MODD_AIRCRAFT_BALLOON : contains balloon and aircraft variables
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!    V. Ducrocq   *Meteo France* 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/05/94 
!!       V. Ducrocq    27/06/94                  
!!       J.Stein       20/10/94 (name of the FMFILE)
!!       J.Stein       06/12/94 add the LS fields   
!!       J.P. Lafore   09/01/95 add the DRYMASST
!!       J.Stein       20/01/95 add TKE and change the ycomment for the water 
!!                              variables       
!!       J.Stein       23/01/95 add a TKE switch and MODD_PARAM_n
!!       J.Stein       16/03/95 remove R from the historical variables
!!       J.Stein       20/03/95 add the EPS var.  
!!       J.Stein       30/06/95 add the variables related to the subgrid condens
!!       S. Belair     01/09/95 add surface variables and ground parameters
!!       J.-P. Pinty   15/09/95 add the radiation parameters
!!       J.Stein       23/01/96 add the TSZ0 option for the surface scheme
!!       M.Georgelin   13/12/95 add the forcing variables 
!!       J.-P. Pinty   15/02/96 add external control for the forcing
!!       J.Stein P.Bougeault  15/03/96 add the cloud fraction and change the
!!                                     surface parameters for TSZ0 option
!!       J.Stein P.Jabouille  30/04/96 add the storage type
!!       J.Stein P.Jabouille  20/05/96 switch for XSIGS and XSRC 
!!       J.Stein              10/10/96 change Xsrc into XSRCM and XRCT
!!       J.P. Lafore          30/07/96 add YFMFILE and HDADFILE writing
!!                                     corresponding to MY_NAME and DAD_NAME (for nesting)
!!       V.Masson      08/10/96 add LTHINSHELL
!!       J.-P. Pinty   15/12/96 add the microphysics (ice)
!!       J.-P. Pinty   11/01/97 add the deep convection
!!       J.-P. Pinty   27/01/97 split the recording of the SV array
!!       J.-P. Pinty   29/01/97 set recording of PRCONV and PACCONV in mm/h and
!!                                                         mm respectively
!!       J. Viviand    04/02/97 convert precipitation rates in mm/h
!!       J.P. Lafore   25/11/96 resolution ratio and position for nesting
!!       J.P. Lafore   26/02/97 adding of "surfacic" LS fields
!!       J.Stein       22/06/97 use the absolute pressure
!!       V.Masson      09/07/97 add directional z0 and Subgrid-Scale Orography
!!       V.Masson      18/08/97 call to fmwrit directly with dates and strings
!!       J.Stein       22/10/97 add the LB fields for U,V,W, THETA, RV....
!!       P.Bechtold    24/01/98 add convective tracer tendencies
!!       P.Jabouille   15/10/98 //
!!       P.Jabouille   25/05/99 replace 'DTRAD_CLONLY' by 'DTRAD_CLLY' (size too long)
!!       J. Stein      20/05/98 remove NXEND and NYEND
!!       V. Masson     04/01/00 remove TSZ0 option
!!       P. Jabouille  03/04/00 write XCIT only for MESONH program
!!       K. Suhre      03/12/99 add chemical variable names                         
!        F.solmon /V.Masson   06/00 adapt for patch surface variables
!!       D.Gazen       22/01/01 use MODD_NSV and add names to scalar variables
!!       G.Jaubert     06/06/01 add Balloon current positions
!!       P.Jabouille   10/04/02 extra radiative surface flux
!!       J.-P. Pinty   29/11/02 add C3R5, ICE2, ICE4, CELEC
!!       V. Masson     01/2004  removes surface (externalization)
!!                     05/2006  Remove KEPS
!!       J. escobar    02/09/2009 missing YDIR for CLDFR variable
!!                     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!       P. Aumond     12/2009 Mean_UM,...
!!       M. Leriche    16/07/10 add ice phase chemical species
!!       C. Barthe     Jan. 2011  add diagnostics for elec
!!       J. Escobar    Feb. 2012  replace MINVAL/MAXVAL by MIN_ll/MAX_ll in OUTPUT_LISTING
!!       P.Peyrille    06/12 2D west african monsoon: ADV forcing and fluxes writing
!!                     AEROSOLS and ozone vertical distribution are also written
!!       M.Tomasini    06/12 2D west african monsoon: nesting for ADV forcing writing
!!       Pialat/Tulet  15/02/2012 add ForeFire variables
!!       J. Escobar    Mars 2014 , missing YDIR="XY" in 1.6 for tendencies fields 
!!       P. Tulet      Nov 2014 accumulated moles of aqueous species that fall at the surface
!!       M.Faivre      2014
!!       C.Lac         Dec.2014 writing past wind fields for centred advection
!!       J.-P. Pinty   Jan 2015 add LNOx and flash map diagnostics
!!       J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!       P. Tulet & M. Leriche    Nov 2015 add mean pH value in the rain at the surface
!!       Modification    01/2016  (JP Pinty) Add LIMA
!!       M.Mazoyer     04/16 : Add supersaturation fields
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 11/07/2016: remove MNH_NCWRIT define
!  V. Vionnet     07/2017: add blowing snow variables
!  JP Chaboureau 27/11/2017: add wind tendency forcing
!  Q. Libois      02/2018: move Diagnostic related to the radiations in radiations.f90
!  P. Wautelet 11/01/2019: bug correction in write XBL_DEPTH->XSBL_DEPTH
!  C. Lac      18/02/2019: add rain fraction as an output field
!  S. Bielli      02/2019: Sea salt: significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Tulet       02/2020: correction for dust and sea salts
!  B. Vie         06/2020: add prognostic supersaturation for LIMA
!  PA. Joulin     12/2020: add wind turbine outputs
!  F. Auguste     02/2021: add IBM
!  T. Nagel       02/2021: add turbulence recycling
!  P. Wautelet 10/03/2021: use scalar variable names for dust and salt
!  P. Wautelet 11/03/2021: bugfix: correct name for NSV_LIMA_IMM_NUCL
!  J.L. Redelsperger 03/2021: add OCEAN and auto-coupled O-A LES cases
!  R. Schoetter   12/2021: adds humidity and other mean diagnostics
!  A. Costes      12/2021: add Blaze fire model
!  P. Wautelet 04/02/2022: use TSVLIST to manage metadata of scalar variables
!  E. Jezequel    11/2022: add covariances from MEAN fields
!  H. Toumi       09/2022: add ADR
!  PA. Joulin     04/2023: update EOL metadata management
!  M. Mandement   01/2024: add max 10 m wind gust speed formulations
!  A. Marcel Jan 2025: EDMF contribution to dynamic TKE production
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_2D_FRC
USE MODD_ADVFRC_n
USE MODD_ADV_n,            ONLY: CUVW_ADV_SCHEME, XRTKEMS, CTEMP_SCHEME, LSPLIT_CFL
USE MODD_AIRCRAFT_BALLOON, ONLY: LFLYER
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
USE MODD_CH_AEROSOL
USE MODD_CH_M9_n
USE MODD_CH_MNHC_n,       ONLY: LUSECHEM,LCH_CONV_LINOX, &
                                LUSECHAQ,LUSECHIC,LCH_PH, XCH_PHINIT
USE MODD_CH_PH_n
USE MODD_CLOUDPAR
USE MODD_CONDSAMP
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST
USE MODD_DEEP_CONVECTION_n
USE MODD_DEF_EDDY_FLUX_n
USE MODD_DEF_EDDYUV_FLUX_n
USE MODD_DIM_n
USE MODD_DUMMY_GR_FIELD_n
USE MODD_DUST
USE MODD_DYN_n
USE MODD_ELEC_DESCR,      ONLY: LLNOX_EXPLICIT
USE MODD_ELEC_FLASH
USE MODD_ELEC_n
USE MODD_EOL_ADNR
USE MODD_EOL_ADR
USE MODD_EOL_ALM
USE MODD_EOL_MAIN
USE MODD_EOL_SHARED_IO
USE MODD_FIELD_n
use modd_field,       only: NMNHDIM_UNUSED, tfieldmetadata, tfieldlist, NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED, &
                            TYPECHAR, TYPEDATE, TYPEINT, TYPELOG, TYPEREAL
USE MODD_FIRE_n
#ifdef MNH_FOREFIRE
USE MODD_FOREFIRE
#endif
USE MODD_FRC
USE MODD_GR_FIELD_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_HURR_CONF, ONLY: LFILTERING,CFILTERING,NDIAG_FILT
USE MODD_HURR_FIELD_n
USE MODD_IBM_LSF,         ONLY: LIBM_LSF
USE MODD_IBM_PARAM_n,     ONLY: LIBM, XIBM_LS, XIBM_XMUT
USE MODD_IO, ONLY: TFILEDATA
USE MODD_LATZ_EDFLX
USE MODD_LIMA_PRECIP_SCAVENGING_n
USE MODD_LSFIELD_n
USE MODD_LUNIT_n
USE MODD_MEAN_FIELD_n
USE MODD_NESTING
USE MODD_NSV
USE MODD_OCEANH
USE MODD_PARAM_C2R2,      ONLY: LSUPSAT
USE MODD_PARAMETERS
USE MODD_PARAM_KAFR_n,      ONLY: LCHTRANS
USE MODD_PARAM_LIMA,     ONLY: LSCAV, LAERO_MASS
USE MODD_PARAM_n
USE MODD_PASPOL
USE MODD_PAST_FIELD_n
USE MODD_PRECIP_n
USE MODD_PREP_REAL, ONLY: CDUMMY_2D, XDUMMY_2D
USE MODD_RADIATIONS_n,   ONLY : XDTHRAD, NCLEARCOL_TM1, XFLALWD, &
                                XZENITH, XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD, &
                                XDIRSRFSWD, XSCAFLASWD, XDIRFLASWD, XAZIM
USE MODD_RECYCL_PARAM_n
USE MODD_REF
USE MODD_REF_n,  ONLY : XRHODREF
USE MODD_RELFRC_n
USE MODD_SALT
USE MODD_TIME
USE MODD_TIME_n
USE MODD_TURB_n

USE MODE_EXTRAPOL
use mode_field, only: Find_field_id_from_mnhname
USE MODE_GRIDPROJ
USE MODE_IO_FIELD_WRITE, only: IO_Field_write
USE MODE_IO_FILE,        only: IO_File_close
USE MODE_ll
USE MODE_MODELN_HANDLER
USE MODE_MPPDB
USE MODE_MSG
USE MODE_TOOLS, ONLY: UPCASE
USE MODE_WRITE_BALLOON_n, ONLY: WRITE_BALLOON_n

USE MODI_CH_AER_REALLFI_n
USE MODI_DUST_FILTER
USE MODI_DUSTLFI_n
USE MODI_SALT_FILTER
USE MODI_SALTLFI_n
USE MODI_SHUMAN
USE MODI_WRITE_LB_n

IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
TYPE(TFILEDATA), INTENT(INOUT) :: TPFILE   ! File characteristics
CHARACTER(LEN=*),INTENT(IN) :: HDADFILE ! Corresponding FM-file name of its DAD model
!
!*       0.2   Declarations of local variables
!
INTEGER           :: ILUOUT         ! logical unit
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears 
                                    !in LFI subroutines at the open of the file              
!
INTEGER           :: JSV            ! loop index for scalar variables
!
CHARACTER(LEN=3)  :: YFRC           ! to mark the time of the forcing
CHARACTER(LEN=1)  :: YFRC1          ! to mark the time of the forcing
INTEGER           :: JT             ! loop index
!
REAL,DIMENSION(:,:), ALLOCATABLE  :: ZWORK2D     ! Working array
REAL,DIMENSION(:,:,:), ALLOCATABLE  :: ZWORK3D     ! Working array
!
REAL                              :: ZLATOR, ZLONOR ! geographical coordinates of 1st mass point
INTEGER :: IMI ! Current model index
!
INTEGER           :: INFO_ll
INTEGER           :: JI,JJ,JK   ! loop index
INTEGER           :: IIU,IJU,IKU,IIB,IJB,IKB,IIE,IJE,IKE ! Arrays bounds
!
INTEGER              :: IDX
INTEGER              :: IID
TYPE(TFIELDMETADATA) :: TZFIELD
!-------------------------------------------------------------------------------
!
!*      0. Initialization
!
IMI = GET_CURRENT_MODEL_INDEX()
!
ILUOUT=TLUOUT%NLU

IDX = 1
!
ALLOCATE(ZWORK2D(SIZE(XTHT,1),SIZE(XTHT,2)))
ALLOCATE(ZWORK3D(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)))
!
!*       0.2     ARRAYS BOUNDS INITIALIZATION
!
IIU=NIMAX+2*JPHEXT
IJU=NJMAX+2*JPHEXT
IKU=NKMAX+2*JPVEXT
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=IKU-JPVEXT
!
!*       1.     WRITES IN THE LFI FILE
! 
!
!*       1.0    File and HDADFILE writing :
!
CALL IO_Field_write(TPFILE,'FILETYPE',TPFILE%CTYPE)
!
IF (LEN_TRIM(HDADFILE)>0) THEN
  CALL IO_Field_write(TPFILE,'DXRATIO',NDXRATIO_ALL(IMI))
  CALL IO_Field_write(TPFILE,'DYRATIO',NDYRATIO_ALL(IMI))
  CALL IO_Field_write(TPFILE,'XOR',    NXOR_ALL(IMI))
  CALL IO_Field_write(TPFILE,'YOR',    NYOR_ALL(IMI))
END IF
!
!*       1.1    Type and Dimensions :
!
CALL IO_Field_write(TPFILE,'IMAX',NIMAX_ll)
CALL IO_Field_write(TPFILE,'JMAX',NJMAX_ll)
CALL IO_Field_write(TPFILE,'KMAX',NKMAX)
!
CALL IO_Field_write(TPFILE,'JPHEXT',JPHEXT)
!
!*       1.2    Grid variables :
!
IF (.NOT.LCARTESIAN) THEN
  CALL IO_Field_write(TPFILE,'RPK',   XRPK)
  CALL IO_Field_write(TPFILE,'LONORI',XLONORI)
  CALL IO_Field_write(TPFILE,'LATORI',XLATORI)
! 
!* diagnostic of 1st mass point
!
  CALL SM_LATLON( XLATORI, XLONORI, XHATM_BOUND(NEXTE_XMIN), XHATM_BOUND(NEXTE_YMIN), ZLATOR, ZLONOR )
!
  CALL IO_Field_write(TPFILE,'LONOR',ZLONOR)
  CALL IO_Field_write(TPFILE,'LATOR',ZLATOR)
END IF 
!
CALL IO_Field_write(TPFILE,'THINSHELL',LTHINSHELL)
CALL IO_Field_write(TPFILE,'LAT0',XLAT0)
CALL IO_Field_write(TPFILE,'LON0',XLON0)
CALL IO_Field_write(TPFILE,'BETA',XBETA)
!
CALL IO_Field_write(TPFILE,'XHAT',XXHAT)
CALL IO_Field_write(TPFILE,'YHAT',XYHAT)
CALL IO_Field_write(TPFILE,'ZHAT',XZHAT)
CALL IO_Field_write(TPFILE,'ZTOP',XZTOP)
!
IF (.NOT.LCARTESIAN) THEN
  CALL IO_Field_write(TPFILE,'LAT',XLAT)
  CALL IO_Field_write(TPFILE,'LON',XLON)
END IF
!
CALL IO_Field_write(TPFILE,'ZS',   XZS)
IF(ASSOCIATED(XZWS)) THEN
  CALL IO_Field_write(TPFILE,'ZWS',  XZWS)
END IF
CALL IO_Field_write(TPFILE,'ZSMT', XZSMT)
CALL IO_Field_write(TPFILE,'SLEVE',LSLEVE)
!
IF (LSLEVE) THEN
  CALL IO_Field_write(TPFILE,'LEN1',XLEN1)
  CALL IO_Field_write(TPFILE,'LEN2',XLEN2)
END IF
!
!
CALL IO_Field_write(TPFILE,'DTMOD',TDTMOD)
CALL IO_Field_write(TPFILE,'DTCUR',TDTCUR)
CALL IO_Field_write(TPFILE,'DTEXP',TDTEXP)
CALL IO_Field_write(TPFILE,'DTSEG',TDTSEG)
!
!*       1.3    Configuration  variables :
!
CALL IO_Field_write(TPFILE,'L1D',      L1D)
CALL IO_Field_write(TPFILE,'L2D',      L2D)
CALL IO_Field_write(TPFILE,'PACK',     LPACK)
CALL IO_Field_write(TPFILE,'CARTESIAN',LCARTESIAN)
CALL IO_Field_write(TPFILE,'LBOUSS',   LBOUSS)
CALL IO_Field_write(TPFILE,'LOCEAN',   LOCEAN)
CALL IO_Field_write(TPFILE,'LCOUPLES', LCOUPLES)
!
CALL IO_Field_write(TPFILE,'SURF',     CSURF)
CALL IO_Field_write(TPFILE,'CPL_AROME',LCPL_AROME)
CALL IO_Field_write(TPFILE,'COUPLING', LCOUPLING)
!
TZFIELD = TFIELDMETADATA(   &
  CMNHNAME   = 'RECYCLING', &
  CLONGNAME  = 'RECYCLING', &
  CSTDNAME   = '',          &
  CUNITS     = '',          &
  CDIR       = '--',        &
  CCOMMENT   = '',          &
  NGRID      = 1,           &
  NTYPE      = TYPELOG,     &
  NDIMS      = 0,           &
  LTIMEDEP   = .FALSE.      )
CALL IO_Field_write(TPFILE,TZFIELD,LRECYCL)
!
!*       1.4    Prognostic variables :
!
!
!*       1.4.1  Time t:
!
!20131128 check XUT-> X_Y_W_U wind component for PRC
!  CALL EXTRAPOL('W',XUT)
!  CALL EXTRAPOL('E',XUT)
!  CALL EXTRAPOL('N',XUT)
!  CALL EXTRAPOL('S',XUT)
CALL MPPDB_CHECK3D(XUT,"write_lfifmn before IO_Field_write::XUT",PRECISION)
CALL IO_Field_write(TPFILE,'UT',XUT)
CALL MPPDB_CHECK3D(XUT,"write_lfifmn after IO_Field_write::XUT",PRECISION)
!
!20131128 check XVT-> X_Y_W_V wind component for PRC
CALL MPPDB_CHECK3D(XVT,"write_lfifmn::XVT",PRECISION)
!
CALL IO_Field_write(TPFILE,'VT',XVT)
CALL IO_Field_write(TPFILE,'WT',XWT)
!
CALL IO_Field_write(TPFILE,'THT',XTHT)
!
IF (LBLAZE) THEN
  CALL IO_Field_write( TPFILE, 'FMREFINRATIOX', NREFINX )
  CALL IO_Field_write( TPFILE, 'FMREFINRATIOY', NREFINY )
  CALL IO_Field_write( TPFILE, 'FMPHI',  XLSPHI )
  CALL IO_Field_write( TPFILE, 'FMBMAP', XBMAP )
  CALL IO_Field_write( TPFILE, 'FMROS0', XFMR0 )
  CALL IO_Field_write( TPFILE, 'FMROS', XFIRERW )
  CALL IO_Field_write( TPFILE, 'FMASE', XFMASE )
  CALL IO_Field_write( TPFILE, 'FMAWC', XFMAWC )
  CALL IO_Field_write( TPFILE, 'FMFLUXHDH', XFMFLUXHDH )
  CALL IO_Field_write( TPFILE, 'FMFLUXHDW', XFMFLUXHDW )
  IF (LWINDFILTER .AND. CWINDFILTER=='WLIM') THEN
    CALL IO_Field_write( TPFILE, 'FMHWS', XFMHWS )
  ELSE
    CALL IO_Field_write( TPFILE, 'FMWINDU', XFMWINDU )
    CALL IO_Field_write( TPFILE, 'FMWINDV', XFMWINDV )
    CALL IO_Field_write( TPFILE, 'FMWINDW', XFMWINDW )
  END IF
  CALL IO_Field_write( TPFILE, 'FMGRADOROX', XFMGRADOROX )
  CALL IO_Field_write( TPFILE, 'FMGRADOROY', XFMGRADOROY )
END IF
!
!*       1.4.2  Time t-dt:
!
IF ( (CUVW_ADV_SCHEME == 'CEN4TH') .AND. (CTEMP_SCHEME == 'LEFR') ) THEN
  CALL IO_Field_write(TPFILE,'UM', XUM)
  CALL IO_Field_write(TPFILE,'VM', XVM)
  CALL IO_Field_write(TPFILE,'WM', XWM)
  CALL IO_Field_write(TPFILE,'DUM',XDUM)
  CALL IO_Field_write(TPFILE,'DVM',XDVM)
  CALL IO_Field_write(TPFILE,'DWM',XDWM)
END IF
!
IF (LIBM .OR. LIBM_LSF) THEN
  !
  TZFIELD = TFIELDMETADATA(                       &
    CMNHNAME  = 'LSFP',                           &
    CLONGNAME = 'LSFP',                           &
    CSTDNAME  = '',                               &
    CUNITS    = 'm',                              &
    CDIR      = 'XY',                             &
    NGRID     = 1,                                &
    NTYPE     = TYPEREAL,                         &
    NDIMS     = 3,                                &
    LTIMEDEP  = .TRUE.,                           &
    CCOMMENT  = 'Level Set Function at mass node' )
  !
  CALL IO_Field_write(TPFILE,TZFIELD,XIBM_LS(:,:,:,1))
  !
  IF ( CPROGRAM == 'MESONH' ) THEN
    TZFIELD = TFIELDMETADATA( &
      CMNHNAME   = 'XMUT',    &
      CLONGNAME  = 'XMUT',    &
      CSTDNAME   = '',        &
      CUNITS     = 'm2 s-1',  &
      CDIR       = 'XY',      &
      NGRID      = 1,         &
      NTYPE      = TYPEREAL,  &
      NDIMS      = 3,         &
      LTIMEDEP   = .TRUE.     )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XIBM_XMUT)
  END IF
   !
ENDIF
!
IF (LRECYCL) THEN
  !
  TZFIELD = TFIELDMETADATA(                                  &
    CMNHNAME   = 'RCOUNT',                                   &
    CLONGNAME  = 'RCOUNT',                                   &
    CSTDNAME   = '',                                         &
    CUNITS     = '',                                         &
    CDIR       = '--',                                       &
    NGRID      = 1,                                          &
    NTYPE      = TYPEINT,                                    &
    NDIMS      = 0,                                          &
    LTIMEDEP   = .TRUE.,                                     &
    CCOMMENT   = 'Incremental counter for averaging purpose' )
  CALL IO_Field_write(TPFILE,TZFIELD,NR_COUNT)
  !
  IF (LRECYCLW) THEN
    TZFIELD = TFIELDMETADATA(                                   &
      CMNHNAME   = 'URECYCLW',                                  &
      CLONGNAME  = 'URECYCLW',                                  &
      CSTDNAME   = '',                                          &
      CUNITS     = 'm s-1',                                     &
      CDIR       = 'XY',                                        &
      NGRID      = 2,                                           &
      NTYPE      = TYPEREAL,                                    &
      NDIMS      = 3,                                           &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                      &
      CCOMMENT   = 'UMEAN-WEST side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XUMEANW(:,:,:))
    !
    TZFIELD = TFIELDMETADATA(                                   &
      CMNHNAME   = 'VRECYCLW',                                  &
      CLONGNAME  = 'VRECYCLW',                                  &
      CSTDNAME   = '',                                          &
      CUNITS     = 'm s-1',                                     &
      CDIR       = 'XY',                                        &
      NGRID      = 3,                                           &
      NTYPE      = TYPEREAL,                                    &
      NDIMS      = 3,                                           &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                      &
      CCOMMENT   = 'VMEAN-WEST side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XVMEANW(:,:,:))
    !
    TZFIELD = TFIELDMETADATA(                                   &
      CMNHNAME   = 'WRECYCLW',                                  &
      CLONGNAME  = 'WRECYCLW',                                  &
      CSTDNAME   = '',                                          &
      CUNITS     = 'm s-1',                                     &
      CDIR       = 'XY',                                        &
      NGRID      = 4,                                           &
      NTYPE      = TYPEREAL,                                    &
      NDIMS      = 3,                                           &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                      &
      CCOMMENT   = 'WMEAN-WEST side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XWMEANW(:,:,:))
    !
  ENDIF  
  IF (LRECYCLN) THEN
    TZFIELD = TFIELDMETADATA(                                    &
      CMNHNAME   = 'URECYCLN',                                   &
      CLONGNAME  = 'URECYCLN',                                   &
      CSTDNAME   = '',                                           &
      CUNITS     = 'm s-1',                                      &
      CDIR       = 'XY',                                         &
      NGRID      = 2,                                            &
      NTYPE      = TYPEREAL,                                     &
      NDIMS      = 3,                                            &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                       &
      CCOMMENT   = 'UMEAN-NORTH side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XUMEANN(:,:,:))
    !
    TZFIELD = TFIELDMETADATA(                                    &
      CMNHNAME   = 'VRECYCLN',                                   &
      CLONGNAME  = 'VRECYCLN',                                   &
      CSTDNAME   = '',                                           &
      CUNITS     = 'm s-1',                                      &
      CDIR       = 'XY',                                         &
      NGRID      = 3,                                            &
      NTYPE      = TYPEREAL,                                     &
      NDIMS      = 3,                                            &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                       &
      CCOMMENT   = 'VMEAN-NORTH side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XVMEANN(:,:,:))
    !
    TZFIELD = TFIELDMETADATA(                                    &
      CMNHNAME   = 'WRECYCLN',                                   &
      CLONGNAME  = 'WRECYCLN',                                   &
      CSTDNAME   = '',                                           &
      CUNITS     = 'm s-1',                                      &
      CDIR       = 'XY',                                         &
      NGRID      = 4,                                            &
      NTYPE      = TYPEREAL,                                     &
      NDIMS      = 3,                                            &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                       &
      CCOMMENT   = 'WMEAN-NORTH side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XWMEANN(:,:,:))
    !
  ENDIF
  IF (LRECYCLE) THEN
    TZFIELD = TFIELDMETADATA(                                   &
      CMNHNAME   = 'URECYCLE',                                  &
      CLONGNAME  = 'URECYCLE',                                  &
      CSTDNAME   = '',                                          &
      CUNITS     = 'm s-1',                                     &
      CDIR       = 'XY',                                        &
      NGRID      = 2,                                           &
      NTYPE      = TYPEREAL,                                    &
      NDIMS      = 3,                                           &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                      &
      CCOMMENT   = 'UMEAN-EAST side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XUMEANE(:,:,:))
    !
    TZFIELD = TFIELDMETADATA(                                   &
      CMNHNAME   = 'VRECYCLE',                                  &
      CLONGNAME  = 'VRECYCLE',                                  &
      CSTDNAME   = '',                                          &
      CUNITS     = 'm s-1',                                     &
      CDIR       = 'XY',                                        &
      NGRID      = 3,                                           &
      NTYPE      = TYPEREAL,                                    &
      NDIMS      = 3,                                           &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                      &
      CCOMMENT   = 'VMEAN-EAST side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XVMEANE(:,:,:))
    !
    TZFIELD = TFIELDMETADATA(                                   &
      CMNHNAME   = 'WRECYCLE',                                  &
      CLONGNAME  = 'WRECYCLE',                                  &
      CSTDNAME   = '',                                          &
      CUNITS     = 'm s-1',                                     &
      CDIR       = 'XY',                                        &
      NGRID      = 4,                                           &
      NTYPE      = TYPEREAL,                                    &
      NDIMS      = 3,                                           &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                      &
      CCOMMENT   = 'WMEAN-EAST side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XWMEANE(:,:,:))
    !
  ENDIF
  IF (LRECYCLS) THEN
    TZFIELD = TFIELDMETADATA(                                    &
      CMNHNAME   = 'URECYCLS',                                   &
      CLONGNAME  = 'URECYCLS',                                   &
      CSTDNAME   = '',                                           &
      CUNITS     = 'm s-1',                                      &
      CDIR       = 'XY',                                         &
      NGRID      = 2,                                            &
      NTYPE      = TYPEREAL,                                     &
      NDIMS      = 3,                                            &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                       &
      CCOMMENT   = 'UMEAN-SOUTH side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XUMEANS(:,:,:))
    !
    TZFIELD = TFIELDMETADATA(                                    &
      CMNHNAME   = 'VRECYCLS',                                   &
      CLONGNAME  = 'VRECYCLS',                                   &
      CSTDNAME   = '',                                           &
      CUNITS     = 'm s-1',                                      &
      CDIR       = 'XY',                                         &
      NGRID      = 3,                                            &
      NTYPE      = TYPEREAL,                                     &
      NDIMS      = 3,                                            &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                       &
      CCOMMENT   = 'VMEAN-SOUTH side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XVMEANS(:,:,:))
    !
    TZFIELD = TFIELDMETADATA(                                    &
      CMNHNAME   = 'WRECYCLS',                                   &
      CLONGNAME  = 'WRECYCLS',                                   &
      CSTDNAME   = '',                                           &
      CUNITS     = 'm s-1',                                      &
      CDIR       = 'XY',                                         &
      NGRID      = 4,                                            &
      NTYPE      = TYPEREAL,                                     &
      NDIMS      = 3,                                            &
      NDIMLIST   = [ NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED ], &
      LTIMEDEP   = .TRUE.,                                       &
      CCOMMENT   = 'WMEAN-SOUTH side plan for recycling purpose' )
    !
    CALL IO_Field_write(TPFILE,TZFIELD,XWMEANS(:,:,:))
    !
  ENDIF  
ENDIF
!
IF (MEAN_COUNT /= 0) THEN
!
  TZFIELD = TFIELDMETADATA(                          &
    CMNHNAME   = 'generic for mean_count variables', & !Temporary name to ease identification
    CSTDNAME   = '',                                 &
    CDIR       = 'XY',                               &
    NTYPE      = TYPEREAL,                           &
    NGRID      = 2,                                  &
    NDIMS      = 3,                                  &
    LTIMEDEP   = .TRUE.                              )
!
  TZFIELD%CMNHNAME   = 'UMME'
  TZFIELD%CLONGNAME  = 'UMME'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_U component of mean wind'
  ZWORK3D = XUM_MEAN
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'U2ME'
  TZFIELD%CLONGNAME  = 'U2ME'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_U component of mean wind variance'
  ZWORK3D = XU2_M2/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'UMMA'
  TZFIELD%CLONGNAME  = 'UMMA'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_U component of max wind'
  CALL IO_Field_write(TPFILE,TZFIELD,XUM_MAX)
!
  IF (LCOV_FIELD) THEN
    !
    TZFIELD%CMNHNAME   = 'UVME'
    TZFIELD%CLONGNAME  = 'UVME'
    TZFIELD%CUNITS     = 'm2 s-2'
    TZFIELD%CCOMMENT   = 'X_Y_Z_UV component of mean wind variance'
    ZWORK3D = XUV_MEAN/MEAN_COUNT
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
    !
    TZFIELD%CMNHNAME   = 'UWME'
    TZFIELD%CLONGNAME  = 'UWME'
    TZFIELD%CUNITS     = 'm2 s-2'
    TZFIELD%CCOMMENT   = 'X_Y_Z_UW component of mean wind variance'
    ZWORK3D = XUW_MEAN/MEAN_COUNT
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
    !
    TZFIELD%CMNHNAME   = 'VWME'
    TZFIELD%CLONGNAME  = 'VWME'
    TZFIELD%CUNITS     = 'm2 s-2'
    TZFIELD%CCOMMENT   = 'X_Y_Z_VW component of mean wind variance'
    ZWORK3D = XVW_MEAN/MEAN_COUNT
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
    !
    TZFIELD%CMNHNAME   = 'WTHME'
    TZFIELD%CLONGNAME  = 'WTHME'
    TZFIELD%CUNITS     = 'm2 s-2'
    TZFIELD%CCOMMENT   = 'X_Y_Z_WTH component of mean wind variance'
    ZWORK3D = XWTH_MEAN/MEAN_COUNT
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
    !
  END IF
  !
  IF (LMINMAX_VORT) THEN
  ! Min and max vorticity
  !
    TZFIELD%CMNHNAME   = 'UM1_MAX'
    TZFIELD%CLONGNAME  = 'UM1_MAX'
    TZFIELD%CUNITS     = 's-1'
    TZFIELD%CCOMMENT   = 'X_Y_Z_x_maximum relative vorticity'
    CALL IO_Field_write(TPFILE,TZFIELD,XUM1_MAX)
  !
  !
    TZFIELD%CMNHNAME   = 'UM1_MIN'
    TZFIELD%CLONGNAME  = 'UM1_MIN'
    TZFIELD%CUNITS     = 's-1'
    TZFIELD%CCOMMENT   = 'X_Y_Z_x_maximum relative vorticity'
    CALL IO_Field_write(TPFILE,TZFIELD,XUM1_MIN)
  END IF
  !
  !
  TZFIELD = TFIELDMETADATA(                          &
    CMNHNAME   = 'generic for mean_count variables', & !Temporary name to ease identification
    CSTDNAME   = '',                                 &
    CDIR       = 'XY',                               &
    NTYPE      = TYPEREAL,                           &
    NGRID      = 3,                                  &
    NDIMS      = 3,                                  &
    LTIMEDEP   = .TRUE.                              )
!
  TZFIELD%CMNHNAME   = 'VMME'
  TZFIELD%CLONGNAME  = 'VMME'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_V component of mean wind'
  ZWORK3D = XVM_MEAN
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'V2ME'
  TZFIELD%CLONGNAME  = 'V2ME'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_V component of mean wind variance'
  ZWORK3D = XV2_M2/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'VMMA'
  TZFIELD%CLONGNAME  = 'VMMA'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_V component of max wind'
  CALL IO_Field_write(TPFILE,TZFIELD,XVM_MAX)
!
  IF (LMINMAX_VORT) THEN
  ! Min and max vorticity
  !
    TZFIELD%CMNHNAME   = 'VM1_MAX'
    TZFIELD%CLONGNAME  = 'VM1_MAX'
    TZFIELD%CUNITS     = 's-1'
    TZFIELD%CCOMMENT   = 'X_Y_Z_y_maximum relative vorticity'
    CALL IO_Field_write(TPFILE,TZFIELD,XVM1_MAX)
!
    TZFIELD%CMNHNAME   = 'VM1_MIN'
    TZFIELD%CLONGNAME  = 'VM1_MIN'
    TZFIELD%CUNITS     = 's-1'
    TZFIELD%CCOMMENT   = 'X_Y_Z_y_minimum relative vorticity'
    CALL IO_Field_write(TPFILE,TZFIELD,XVM1_MIN)
  END IF
!
  TZFIELD = TFIELDMETADATA(                          &
    CMNHNAME   = 'generic for mean_count variables', & !Temporary name to ease identification
    CSTDNAME   = '',                                 &
    CDIR       = 'XY',                               &
    NTYPE      = TYPEREAL,                           &
    NGRID      = 4,                                  &
    NDIMS      = 3,                                  &
    LTIMEDEP   = .TRUE.                              )
!
  TZFIELD%CMNHNAME   = 'WMME'
  TZFIELD%CLONGNAME  = 'WMME'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_vertical mean wind'
  ZWORK3D = XWM_MEAN
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'W2ME'
  TZFIELD%CLONGNAME  = 'W2ME'
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_vertical mean wind variance'
  ZWORK3D = XW2_M2/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'WMMA'
  TZFIELD%CLONGNAME  = 'WMMA'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_vertical max wind'
  CALL IO_Field_write(TPFILE,TZFIELD,XWM_MAX)
!
  TZFIELD%CMNHNAME   = 'WMMI'
  TZFIELD%CLONGNAME  = 'WMMI'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_vertical min wind'
  CALL IO_Field_write(TPFILE,TZFIELD,XWM_MIN)
!
  IF (LMINMAX_VORT) THEN
  ! Min and max vorticity
  !
    TZFIELD%CMNHNAME   = 'WM1_MAX'
    TZFIELD%CLONGNAME  = 'WM1_MAX'
    TZFIELD%CUNITS     = 's-1'
    TZFIELD%CCOMMENT   = 'X_Y_Z_z_maximum relative vorticity'
    CALL IO_Field_write(TPFILE,TZFIELD,XWM1_MAX)
!
    TZFIELD%CMNHNAME   = 'WM1_MIN'
    TZFIELD%CLONGNAME  = 'WM1_MIN'
    TZFIELD%CUNITS     = 's-1'
    TZFIELD%CCOMMENT   = 'X_Y_Z_z_maximum relative vorticity'
    CALL IO_Field_write(TPFILE,TZFIELD,XWM1_MIN)
  END IF
!
  !
  ! Calculation of mean horizontal wind speed and
  ! wind direction based on average components
  !
  XWIFF_MEAN = SQRT((MXF(XUM_MEAN)/MEAN_COUNT)**2 + (MYF(XVM_MEAN)/MEAN_COUNT)**2)
  XWIDD_MEAN = 180.0 + (90.0 - 180.0*ATAN2(MYF(XVM_MEAN)/MEAN_COUNT,MXF(XUM_MEAN)/MEAN_COUNT)/XPI)
  !
  WHERE (XWIDD_MEAN(:,:,:).GT.360.0)
     XWIDD_MEAN(:,:,:) = XWIDD_MEAN(:,:,:) - 360.0
  ENDWHERE
  !
  IF ( MINVAL(XWIDD_MEAN) < 0. .OR. MAXVAL(XWIDD_MEAN) > 360. ) &
    call Print_msg( NVERB_FATAL, 'GEN', 'WRITE_LFIFM_n', 'Wrong wind direction' )
  !
  TZFIELD%CMNHNAME   = 'WIFFME'
  TZFIELD%CLONGNAME  = 'WIFFME'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_horizontal mean wind speed'
  CALL IO_Field_write(TPFILE,TZFIELD,XWIFF_MEAN)
  !
  TZFIELD%CMNHNAME   = 'WIDDME'
  TZFIELD%CLONGNAME  = 'WIDDME'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_horizontal mean wind direction (degrees from north)'
  CALL IO_Field_write(TPFILE,TZFIELD,XWIDD_MEAN)
  !
  TZFIELD%CMNHNAME   = 'WIFFMAX'
  TZFIELD%CLONGNAME  = 'WIFFMAX'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_horizontal maximum wind speed'
  CALL IO_Field_write(TPFILE,TZFIELD,XWIFF_MAX)
  XWIFF_MAX(:,:,:)=XNEGUNDEF
  !
  TZFIELD%CMNHNAME   = 'WIDDMAX'
  TZFIELD%CLONGNAME  = 'WIDDMAX'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_horizontal maximum wind direction'
  CALL IO_Field_write(TPFILE,TZFIELD,XWIDD_MAX)
  XWIDD_MAX(:,:,:)=XNEGUNDEF
!
  TZFIELD%NGRID      = 1
  TZFIELD = TFIELDMETADATA(                          &
    CMNHNAME   = 'generic for mean_count variables', & !Temporary name to ease identification
    CSTDNAME   = '',                                 &
    CDIR       = 'XY',                               &
    NTYPE      = TYPEREAL,                           &
    NGRID      = 1,                                  &
    NDIMS      = 3,                                  &
    LTIMEDEP   = .TRUE.                              )
!
  TZFIELD%CMNHNAME   = 'CMME'
  TZFIELD%CLONGNAME  = 'CMME'
  TZFIELD%CUNITS     = 'kg kg-1'
  TZFIELD%CCOMMENT   = 'mean Passive scalar'
  ZWORK3D = XSVT_MEAN/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'THMME'
  TZFIELD%CLONGNAME  = 'THMME'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean potential temperature'
  ZWORK3D = XTHM_MEAN
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'TH2ME'
  TZFIELD%CLONGNAME  = 'TH2ME'
  TZFIELD%CUNITS     = 'K2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean potential temperature variance'
  ZWORK3D = XTH2_M2/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'THMMA'
  TZFIELD%CLONGNAME  = 'THMMA'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CCOMMENT   = 'X_Y_Z_max potential temperature'
  CALL IO_Field_write(TPFILE,TZFIELD,XTHM_MAX)
!
  TZFIELD%CMNHNAME   = 'TEMPMME'
  TZFIELD%CLONGNAME  = 'TEMPMME'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean temperature'
  ZWORK3D = XTEMPM_MEAN
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'TEMP2ME'
  TZFIELD%CLONGNAME  = 'TEMP2ME'
  TZFIELD%CUNITS     = 'K2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean temperature variance'
  ZWORK3D = XTEMP2_M2/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'TEMPMMA'
  TZFIELD%CLONGNAME  = 'TEMPMMA'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CCOMMENT   = 'X_Y_Z_max temperature'
  CALL IO_Field_write(TPFILE,TZFIELD,XTEMPM_MAX)
!
  TZFIELD%CMNHNAME   = 'QSPECME'
  TZFIELD%CLONGNAME  = 'QSPECME'
  TZFIELD%CUNITS     = 'kg kg-1'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean specific humidity'
  ZWORK3D = XQ_MEAN/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'RHME'
  TZFIELD%CLONGNAME  = 'RHME'
  TZFIELD%CUNITS     = 'percent'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean relative humidity, water'
  ZWORK3D = XRH_W_MEAN/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'RHME_ICE'
  TZFIELD%CLONGNAME  = 'RHME_ICE'
  TZFIELD%CUNITS     = 'percent'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean relative humidity, ice'
  ZWORK3D = XRH_I_MEAN/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'RHME_WEIG'
  TZFIELD%CLONGNAME  = 'RHME_WEIG'
  TZFIELD%CUNITS     = 'percent'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean relative humidity, weighted'
  ZWORK3D = XRH_P_MEAN/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'RHMAXME'
  TZFIELD%CLONGNAME  = 'RHMAXME'
  TZFIELD%CUNITS     = 'percent'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean column maximum relative humidity, water'
  ZWORK2D = XRH_W_MAXCOL_MEAN/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'RHMAXME_ICE'
  TZFIELD%CLONGNAME  = 'RHMAXME_ICE'
  TZFIELD%CUNITS     = 'percent'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean column maximum relative humidity, ice'
  ZWORK2D = XRH_I_MAXCOL_MEAN/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'RHMAXME_WEIG'
  TZFIELD%CLONGNAME  = 'RHMAXME_WEIG'
  TZFIELD%CUNITS     = 'percent'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean column maximum relative humidity, weighted'
  ZWORK2D = XRH_P_MAXCOL_MEAN/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'PABSMME'
  TZFIELD%CLONGNAME  = 'PABSMME'
  TZFIELD%CUNITS     = 'Pa'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean ABSolute Pressure'
  ZWORK3D = XPABSM_MEAN
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
  TZFIELD%CMNHNAME   = 'PABS2ME'
  TZFIELD%CLONGNAME  = 'PABS2ME'
  TZFIELD%CUNITS     = 'Pa2'
  TZFIELD%CCOMMENT   = 'X_Y_Z_mean ABSolute Pressure variance'
  ZWORK3D = XPABS2_M2/MEAN_COUNT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
  !
  TZFIELD%CMNHNAME   = 'PABSMMA'
  TZFIELD%CLONGNAME  = 'PABSMMA'
  TZFIELD%CUNITS     = 'Pa'
  TZFIELD%CCOMMENT   = 'X_Y_Z_max ABSolute Pressure'
  CALL IO_Field_write(TPFILE,TZFIELD,XPABSM_MAX)
!
  IF (CTURB /= 'NONE') THEN
    TZFIELD%CMNHNAME   = 'TKEMME'
    TZFIELD%CLONGNAME  = 'TKEMME'
    TZFIELD%CUNITS     = 'm2 s-2'
    TZFIELD%CCOMMENT   = 'X_Y_Z_mean kinetic energy'
    ZWORK3D= XTKEM_MEAN/MEAN_COUNT
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
!
    TZFIELD%CMNHNAME   = 'TKEMMA'
    TZFIELD%CLONGNAME  = 'TKEMMA'
    TZFIELD%CUNITS     = 'm2 s-2'
    TZFIELD%CCOMMENT   = 'X_Y_Z_max kinetic energy'
    CALL IO_Field_write(TPFILE,TZFIELD,XTKEM_MAX)
  END IF
!
  IF (LMINMAX_WINDFFTKE) THEN
    TZFIELD = TFIELDMETADATA(                          &
      CMNHNAME   = 'generic for mean_count variables', & !Temporary name to ease identification
      CSTDNAME   = '',                                 &
      CDIR       = 'XY',                               &
      NTYPE      = TYPEREAL,                           &
      NGRID      = 1,                                  &
      NDIMS      = 2,                                  &
      LTIMEDEP   = .TRUE.                              )
  !
    TZFIELD%CMNHNAME   = 'WMOD10MAX_MA'
    TZFIELD%CLONGNAME  = 'WMOD10MAX_MA'
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CCOMMENT   = 'X_Y_max WMOD10MAX'
    CALL IO_Field_write(TPFILE,TZFIELD,XWMOD10MAX_MAX)
  !
    IF (CTURB /= 'NONE') THEN
     TZFIELD%CMNHNAME   = 'XTKEMAX_MA'
     TZFIELD%CLONGNAME  = 'XTKEMAX_MA'
     TZFIELD%CUNITS     = 'm s-1'
     TZFIELD%CCOMMENT   = 'X_Y_max XTKEMAX'
     CALL IO_Field_write(TPFILE,TZFIELD,XTKEMAX_MAX)
  !
     TZFIELD%CMNHNAME   = 'XTKE10MAX_MA'
     TZFIELD%CLONGNAME  = 'XTKE10MAX_MA'
     TZFIELD%CUNITS     = 'm s-1'
     TZFIELD%CCOMMENT   = 'X_Y_max XTKE10MAX'
     CALL IO_Field_write(TPFILE,TZFIELD,XTKE10MAX_MAX)
  !
     TZFIELD%CMNHNAME   = 'XTKE20MAX_MA'
     TZFIELD%CLONGNAME  = 'XTKE20MAX_MA'
     TZFIELD%CUNITS     = 'm s-1'
     TZFIELD%CCOMMENT   = 'X_Y_max XTKE20MAX'
     CALL IO_Field_write(TPFILE,TZFIELD,XTKE20MAX_MAX)
    END IF
  !
    TZFIELD%CMNHNAME   = 'FF10MAX_MA'
    TZFIELD%CLONGNAME  = 'FF10MAX_MA'
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CCOMMENT   = 'X_Y_max FF10MAX'
    CALL IO_Field_write(TPFILE,TZFIELD,XFF10MAX_MAX)
  !
    TZFIELD%CMNHNAME   = 'FF10MAX2_MA'
    TZFIELD%CLONGNAME  = 'FF10MAX2_MA'
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CCOMMENT   = 'X_Y_max FF10MAX2'
    CALL IO_Field_write(TPFILE,TZFIELD,XFF10MAX2_MAX)
  !
    TZFIELD%CMNHNAME   = 'FF10MAX_AROME_MA'
    TZFIELD%CLONGNAME  = 'FF10MAX_AROME_MA'
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CCOMMENT   = 'X_Y_max FF10MAX_AROME'
    CALL IO_Field_write(TPFILE,TZFIELD,XFF10MAX_AROME_MAX)
  !
  END IF

  IF (LMINMAX_MSLP) THEN
  ! Min and max sea level pressure
    TZFIELD%CMNHNAME   = 'MSLP_MAX'
    TZFIELD%CLONGNAME  = 'MSLP_MAX'
    TZFIELD%CUNITS     = 'hPa'
    TZFIELD%CCOMMENT   = 'X_Y_max Mean Sea Level Pressure'
    CALL IO_Field_write(TPFILE,TZFIELD,XMSLP_MAX)
  !
    TZFIELD%CMNHNAME   = 'MSLP_MIN'
    TZFIELD%CLONGNAME  = 'MSLP_MIN'
    TZFIELD%CUNITS     = 'hPa'
    TZFIELD%CCOMMENT   = 'X_Y_min Mean Sea Level Pressure'
    CALL IO_Field_write(TPFILE,TZFIELD,XMSLP_MIN)
  !
  END IF
!
  ! Max updraft helicity
  IF (LUH_MAX) THEN
    TZFIELD = TFIELDMETADATA(                          &
      CMNHNAME   = 'generic for mean_count variables', & !Temporary name to ease identification
      CSTDNAME   = '',                                 &
      CDIR       = 'XY',                               &
      NTYPE      = TYPEREAL,                           &
      NGRID      = 4,                                  &
      NDIMS      = 2,                                  &
      LTIMEDEP   = .TRUE.                              )
!
    TZFIELD%CMNHNAME   = 'UH_MAX'
    TZFIELD%CLONGNAME  = 'UH_MAX'
    TZFIELD%CUNITS     = 'm2 s-2'
    TZFIELD%CCOMMENT   = 'X_Y_max Updraft Helicity'
    CALL IO_Field_write(TPFILE,TZFIELD,XUH_MAX)
  END IF
!
END IF
!
!
IF (CTURB /= 'NONE') THEN
  CALL IO_Field_write(TPFILE,'TKET',XTKET)
  IF (CPROGRAM == 'MESONH' .AND. LSPLIT_CFL) CALL IO_Field_write(TPFILE,'TKEMS',XRTKEMS)
END IF
!
!
CALL IO_Field_write(TPFILE,'PABST',XPABST)
!
IF (NRR >=1) THEN
  IF (LUSERV) CALL IO_Field_write(TPFILE,'RVT',XRT(:,:,:,IDX_RVT))
  IF (LUSERC) THEN
    CALL IO_Field_write(TPFILE,'RCT',XRT(:,:,:,IDX_RCT))
    WRITE (ILUOUT,*) IDX_RCT,' RC min-max ',MIN_ll(XRT(:,:,:,IDX_RCT),INFO_ll),MAX_ll(XRT(:,:,:,IDX_RCT),INFO_ll)
  END IF
  IF (LUSERR) THEN
    CALL IO_Field_write(TPFILE,'RRT',XRT(:,:,:,IDX_RRT))
    WRITE (ILUOUT,*) IDX_RRT,' RR min-max ',MIN_ll(XRT(:,:,:,IDX_RRT),INFO_ll),MAX_ll(XRT(:,:,:,IDX_RRT),INFO_ll)
  END IF 
  IF (LUSERI) THEN
    CALL IO_Field_write(TPFILE,'RIT',XRT(:,:,:,IDX_RIT))
    WRITE (ILUOUT,*) IDX_RIT,' RI min-max ',MIN_ll(XRT(:,:,:,IDX_RIT),INFO_ll),MAX_ll(XRT(:,:,:,IDX_RIT),INFO_ll)
    IF ( CPROGRAM == 'MESONH' .AND. CCLOUD(1:3) == 'ICE') THEN
      CALL IO_Field_write(TPFILE,'CIT',XCIT(:,:,:))
    END IF
  END IF 
  IF (LUSERS) THEN
    CALL IO_Field_write(TPFILE,'RST',XRT(:,:,:,IDX_RST))
    WRITE (ILUOUT,*) IDX_RST,' RS min-max ',MINVAL(XRT(:,:,:,IDX_RST)),MAXVAL(XRT(:,:,:,IDX_RST))
  END IF
  IF (LUSERG) THEN
    CALL IO_Field_write(TPFILE,'RGT',XRT(:,:,:,IDX_RGT))
    WRITE (ILUOUT,*) IDX_RGT,' RG min-max ',MINVAL(XRT(:,:,:,IDX_RGT)),MAXVAL(XRT(:,:,:,IDX_RGT))
  END IF 
  IF (LUSERH) CALL IO_Field_write(TPFILE,'RHT',XRT(:,:,:,IDX_RHT))
END IF
!
IF (NSV >= 1 ) THEN
  ! aerosol scalar variables
  IF ( LORILAM ) THEN
    IF ((CPROGRAM == 'REAL  ').AND.(NSV_AER > 1).AND.(IMI==1).AND.(LAERINIT))  &
      CALL CH_AER_REALLFI_n(XSVT(:,:,:,NSV_AERBEG:NSV_AEREND),XSVT(:,:,:,NSV_CHEMBEG-1+JP_CH_CO), XRHODREF)
    IF ((CPROGRAM == 'IDEAL ').AND.(NSV_AER > 1).AND.(IMI==1))  &
      CALL CH_AER_REALLFI_n(XSVT(:,:,:,NSV_AERBEG:NSV_AEREND),XSVT(:,:,:,NSV_CHEMBEG-1+JP_CH_CO),  XRHODREF)
  END IF

  ! dust scalar variables
  IF ( LDUST ) THEN
    IF ((CPROGRAM == 'REAL  ').AND.(NSV_DST > 1).AND.(IMI==1).AND.(LDSTINIT).AND.(.NOT.LDSTCAMS)) &
      CALL DUSTLFI_n(XSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND), XRHODREF)
    IF ((CPROGRAM == 'IDEAL ').AND.(NSV_DST > 1).AND.(IMI==1)) &
      CALL DUSTLFI_n(XSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND), XRHODREF)
    !At this point, we have the tracer array in order of importance, i.e.
    !if mode 2 is most important it will occupy place 1-3 of XSVT
    CALL DUST_FILTER(XSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND), XRHODREF)
  END IF

  ! sea salt scalar variables
  IF ( LSALT ) THEN
    IF ((CPROGRAM == 'REAL  ').AND.(NSV_SLT > 1).AND.(IMI==1).AND.(LSLTINIT).AND.(.NOT.LSLTCAMS)) &
      CALL SALTLFI_n(XSVT(:,:,:,NSV_SLTBEG:NSV_SLTEND), XRHODREF, XZZ)
    IF ((CPROGRAM == 'IDEAL ').AND.(NSV_SLT > 1).AND.(IMI==1)) &
      CALL SALTLFI_n(XSVT(:,:,:,NSV_SLTBEG:NSV_SLTEND), XRHODREF, XZZ)
    !At this point, we have the tracer array in order of importance, i.e.
    !if mode 2 is most important it will occupy place 1-3 of XSVT
    CALL SALT_FILTER(XSVT(:,:,:,NSV_SLTBEG:NSV_SLTEND), XRHODREF)
  END IF

  !Store all scalar variables
  DO JSV = 1, NSV
    CALL IO_Field_write( TPFILE, TSVLIST(JSV), XSVT(:,:,:,JSV) )
  END DO

  IF (LSCAV .AND. LAERO_MASS) THEN
  IF (ASSOCIATED(XINPAP)) THEN
  IF (SIZE(XINPAP) /= 0 ) THEN
    CALL IO_Field_write(TPFILE,'INPAP',XINPAP)
    !
    ZWORK2D(:,:)  = XRHOLW*XINPRR(:,:)*XSVT(:,:,2,NSV_LIMA_SCAVMASS)/ &
                                        max( 1.e-20,XRT(:,:,2,3) ) !~2=at ground level
    TZFIELD = TFIELDMETADATA(                                      &
      CMNHNAME   = 'INPBP',                                        &
      CSTDNAME   = '',                                             &
      CLONGNAME  = 'INPBP',                                        &
      CUNITS     = 'kg m-2 s-1',                                   &
      CDIR       = 'XY',                                           &
      CCOMMENT   = 'X_Y_INstantaneous Precipitating Aerosol Rate', &
      NGRID      = 1,                                              &
      NTYPE      = TYPEREAL,                                       &
      NDIMS      = 2,                                              &
      LTIMEDEP   = .TRUE.                                          )
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK2D)
    !
    CALL IO_Field_write(TPFILE,'ACPAP',XACPAP)
  END IF
  END IF
  END IF

IF ((LORILAM).AND.(CPROGRAM == 'MESONH')) THEN
DO JSV = 1 , NSV_AER
    TZFIELD = TFIELDMETADATA(                             &
      CMNHNAME   = 'FLX_'//TRIM(UPCASE(CAERONAMES(JSV))), &
      CSTDNAME   = '',                                    &
      CLONGNAME  = 'FLX_'//TRIM(UPCASE(CAERONAMES(JSV))), &
      CUNITS     = 'kg m-2 s-1',                          &
      CDIR       = 'XY',                                  &
      CCOMMENT   = 'Aerosol mass flux',                   &
      NGRID      = 1,                                     &
      NTYPE      = TYPEREAL,                              &
      NDIMS      = 2,                                     &
      LTIMEDEP   = .TRUE.                                 )
    CALL IO_Field_write(TPFILE,TZFIELD,XFLX_AER(:,:,JSV))

    TZFIELD = TFIELDMETADATA(                                          &
      CMNHNAME   = 'FLXT_'//TRIM(UPCASE(CAERONAMES(JSV))),             &
      CSTDNAME   = '',                                                 &
      CLONGNAME  = 'FLXT_'//TRIM(UPCASE(CAERONAMES(JSV))),             &
      CUNITS     = 'kg m-2',                                           &
      CDIR       = 'XY',                                               &
      CCOMMENT   = 'Integrated aerosol mass flux since start/restart', &
      NGRID      = 1,                                                  &
      NTYPE      = TYPEREAL,                                           &
      NDIMS      = 2,                                                  &
      LTIMEDEP   = .TRUE.                                              )
    CALL IO_Field_write(TPFILE,TZFIELD,XFLXT_AER(:,:,JSV))
END DO
END IF

IF ((LSALT).AND.(CPROGRAM == 'MESONH')) THEN
DO JSV = 1 , NMODE_SLT
    WRITE (YFRC1,'(I1.1)') JSV

    TZFIELD = TFIELDMETADATA(              &
      CMNHNAME   = 'FLX_SLT'//YFRC1,       &
      CSTDNAME   = '',                     &
      CLONGNAME  = 'FLX_SLT'//YFRC1,       &
      CUNITS     = 'm-2 s-1',              &
      CDIR       = 'XY',                   &
      CCOMMENT   = 'Sea salt number flux', &
      NGRID      = 1,                      &
      NTYPE      = TYPEREAL,               &
      NDIMS      = 2,                      &
      LTIMEDEP   = .TRUE.                  )
    CALL IO_Field_write(TPFILE,TZFIELD,XFLX_SLT(:,:,JSV))

    TZFIELD = TFIELDMETADATA(                                             &
      CMNHNAME   = 'FLXT_SLT'//YFRC1,                                     &
      CSTDNAME   = '',                                                    &
      CLONGNAME  = 'FLXT_SLT'//YFRC1,                                     &
      CUNITS     = 'm-2',                                                 &
      CDIR       = 'XY',                                                  &
      CCOMMENT   = 'Integrated sea salt number flux since start/restart', &
      NGRID      = 1,                                                     &
      NTYPE      = TYPEREAL,                                              &
      NDIMS      = 2,                                                     &
      LTIMEDEP   = .TRUE.                                                 )
    CALL IO_Field_write(TPFILE,TZFIELD,XFLXT_SLT(:,:,JSV))
END DO
END IF

IF ((LUSECHEM).AND.(CPROGRAM == 'MESONH')) THEN
    TZFIELD = TFIELDMETADATA(            &
      CMNHNAME   = 'FLX_DMS',            &
      CSTDNAME   = '',                   &
      CLONGNAME  = 'FLX_DMS',            &
      CUNITS     = 'kg m-2 s-1',         &
      CDIR       = 'XY',                 &
      CCOMMENT   = 'Sea salt mass flux', &
      NGRID      = 1,                    &
      NTYPE      = TYPEREAL,             &
      NDIMS      = 2,                    &
      LTIMEDEP   = .TRUE.                )
    CALL IO_Field_write(TPFILE,TZFIELD,XFLX_DMS(:,:))

    TZFIELD = TFIELDMETADATA(                                           &
      CMNHNAME   = 'FLXT_DMS',                                          &
      CSTDNAME   = '',                                                  &
      CLONGNAME  = 'FLXT_DMS',                                          &
      CUNITS     = 'kg m-2',                                            &
      CDIR       = 'XY',                                                &
      CCOMMENT   = 'Integrated sea salt mass flux since start/restart', &
      NGRID      = 1,                                                   &
      NTYPE      = TYPEREAL,                                            &
      NDIMS      = 2,                                                   &
      LTIMEDEP   = .TRUE.                                               )
    CALL IO_Field_write(TPFILE,TZFIELD,XFLXT_DMS(:,:))
END IF



  ! electrical scalar variables
  IF (CELEC /= 'NONE') THEN
    CALL IO_Field_write(TPFILE,'EFIELDU',XEFIELDU)
    CALL IO_Field_write(TPFILE,'EFIELDV',XEFIELDV)
    CALL IO_Field_write(TPFILE,'EFIELDW',XEFIELDW)
 !
    TZFIELD = TFIELDMETADATA(       &
      CMNHNAME   = 'EMODULE',       &
      CSTDNAME   = '',              &
      CLONGNAME  = 'EMODULE',       &
      CUNITS     = 'V m-1',         &
      CDIR       = 'XY',            &
      CCOMMENT   = 'X_Y_Z_EMODULE', &
      NGRID      = 1,               &
      NTYPE      = TYPEREAL,        &
      NDIMS      = 3,               &
      LTIMEDEP   = .TRUE.           )
    ZWORK3D(:,:,:) = (XEFIELDU**2 + XEFIELDV**2 + XEFIELDW**2)**0.5
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK3D)
 !
    CALL FIND_FIELD_ID_FROM_MNHNAME('NI_IAGGS',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'pC m-3 s-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XNI_IAGGS*1.E12)
 !
    CALL FIND_FIELD_ID_FROM_MNHNAME('NI_IDRYG',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'pC m-3 s-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XNI_IDRYG*1.E12)
 !
    CALL FIND_FIELD_ID_FROM_MNHNAME('NI_SDRYG',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'pC m-3 s-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XNI_SDRYG*1.E12)
 !
    CALL FIND_FIELD_ID_FROM_MNHNAME('INDUC_CG',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'pC m-3 s-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XIND_RATE*1.E12)
 !
    CALL IO_Field_write(TPFILE,'TRIG_IC',   NMAP_TRIG_IC)
    CALL IO_Field_write(TPFILE,'IMPACT_CG', NMAP_IMPACT_CG)
    CALL IO_Field_write(TPFILE,'AREA_CG',   NMAP_2DAREA_CG)
    CALL IO_Field_write(TPFILE,'AREA_IC',   NMAP_2DAREA_IC)
    CALL IO_Field_write(TPFILE,'FLASH_3DCG',NMAP_3DCG)
    CALL IO_Field_write(TPFILE,'FLASH_3DIC',NMAP_3DIC)
  END IF
!
  IF ( ((CCLOUD == 'KHKO') .OR.(CCLOUD == 'C2R2')) .AND. (.NOT. LSUPSAT)) THEN
    CALL IO_Field_write(TPFILE,'SUPSATMAX',XSUPSAT(:,:,:))
    CALL IO_Field_write(TPFILE,'NACT',     XNACT(:,:,:))
  END IF
  IF ( ((CCLOUD == 'KHKO') .OR.(CCLOUD == 'C2R2')) .AND. LSUPSAT) THEN
    CALL IO_Field_write(TPFILE,'SSPRO',XSSPRO(:,:,:))
    CALL IO_Field_write(TPFILE,'NPRO', XNPRO(:,:,:))
  END IF
!
  IF (LUSECHEM) THEN
    IF (LUSECHAQ.AND.NRR>=3) THEN ! accumulated moles of aqueous species that fall at the surface (mol/m2)
      DO JSV = NSV_CHACBEG + NSV_CHAC / 2, NSV_CHACEND
        TZFIELD = TSVLIST(JSV)
        TZFIELD%CMNHNAME   = 'ACPR_' // TRIM( TZFIELD%CMNHNAME )
        TZFIELD%CLONGNAME  = 'ACPR_' // TRIM( TZFIELD%CLONGNAME )
        TZFIELD%CUNITS     = 'mol m-2'
        TZFIELD%CCOMMENT   = 'X_Y_Accumulated moles of aqueous species at the surface'
        TZFIELD%NDIMS      = 2
        TZFIELD%NDIMLIST(3)  = TZFIELD%NDIMLIST(4) ! Necessary if LTIMEDEP=.TRUE.
        TZFIELD%NDIMLIST(4:) = NMNHDIM_UNUSED
        ZWORK2D(:,:)  = XACPRAQ(:,:,JSV-NSV_CHACBEG-NSV_CHAC/2+1)
        CALL IO_Field_write(TPFILE,TZFIELD,ZWORK2D)
      END DO
    END IF
    IF (LUSECHAQ.AND.LCH_PH) THEN  ! pH values in cloud
      CALL IO_Field_write(TPFILE,'PHC',XPHC)
      IF (NRR>=3) THEN
        CALL IO_Field_write(TPFILE,'PHR',XPHR)
        ! compute mean pH in accumulated surface water
        !ZWORK2D(:,:) = 10**(-XCH_PHINIT)
        WHERE (XACPRR > 0.)
          ZWORK2D(:,:) =  XACPHR(:,:) *1E3 / XACPRR(:,:) ! moles of H+ / l of water
        ELSE WHERE
          ZWORK2D(:,:) = XUNDEF
        END WHERE
        WHERE ((ZWORK2D(:,:) < 1E-1).AND.(ZWORK2D(:,:) > 1E-14))
          ZWORK2D(:,:) = -LOG10(ZWORK2D(:,:))           ! mean pH of surface water
        END WHERE
        TZFIELD = TFIELDMETADATA(     &
          CMNHNAME   = 'MEANPHR',     &
          CSTDNAME   = '',            &
          CLONGNAME  = 'MEANPHR',     &
          CUNITS     = '1',           &
          CDIR       = 'XY',          &
          CCOMMENT   = 'X_Y_MEAN_PH', &
          NGRID      = 1,             &
          NTYPE      = TYPEREAL,      &
          NDIMS      = 2,             &
          LTIMEDEP = .TRUE.           )
        CALL IO_Field_write(TPFILE,TZFIELD,ZWORK2D)
      ENDIF
    ENDIF
  ENDIF
  !
  TZFIELD = TFIELDMETADATA(                      &
    CMNHNAME   = 'NSVCHEM',                      &
    CSTDNAME   = '',                             &
    CLONGNAME  = 'NSVCHEM',                      &
    CUNITS     = '',                             &
    CDIR       = '--',                           &
    CCOMMENT   = 'Number of chemical variables', &
    NGRID      = 0,                              &
    NTYPE      = TYPEINT,                        &
    NDIMS      = 0,                              &
    LTIMEDEP   = .FALSE.                         )
  CALL IO_Field_write(TPFILE,TZFIELD,NSV_CHEM_LIST)
  !
  IF ( NSV_CHEM_LIST > 0 ) THEN
    TZFIELD = TFIELDMETADATA(                    &
      CMNHNAME   = 'SV_CHEM_LIST',               &
      CSTDNAME   = '',                           &
      CLONGNAME  = 'SV_CHEM_LIST',               &
      CUNITS     = '',                           &
      CDIR       = '--',                         &
      CCOMMENT   = 'List of chemical variables', &
      NGRID      = 0,                            &
      NTYPE      = TYPECHAR,                     &
      NDIMS      = 1,                            &
      LTIMEDEP   = .FALSE.                       )
    CALL IO_Field_write(TPFILE,TZFIELD,CSV_CHEM_LIST)
  END IF
END IF
!
CALL IO_Field_write(TPFILE,'LSUM', XLSUM)
CALL IO_Field_write(TPFILE,'LSVM', XLSVM)
CALL IO_Field_write(TPFILE,'LSWM', XLSWM)
CALL IO_Field_write(TPFILE,'LSTHM',XLSTHM)
IF (LUSERV) CALL IO_Field_write(TPFILE,'LSRVM',XLSRVM)
!
CALL WRITE_LB_n(TPFILE)
!
!
CALL IO_Field_write(TPFILE,'DRYMASST',XDRYMASST)
IF (CPROGRAM == 'MESONH') THEN
  CALL IO_Field_write(TPFILE,'DRYMASSS',XDRYMASSS)
ELSE
  CALL IO_Field_write(TPFILE,'DRYMASSS',0.)
END IF
!
IF( CTURB /= 'NONE' .AND. CTOM=='TM06') THEN
  CALL IO_Field_write(TPFILE,'BL_DEPTH',XBL_DEPTH)
END IF
!
IF( CTURB /= 'NONE' .AND. LRMC01) THEN
  CALL IO_Field_write(TPFILE,'SBL_DEPTH',XSBL_DEPTH)
END IF
!
IF( CTURB /= 'NONE' .AND. CSCONV == 'EDKF' .AND.(CPROGRAM == 'MESONH' .OR. CPROGRAM == 'DIAG')) THEN
  CALL IO_Field_write(TPFILE,'WTHVMF',XWTHVMF)
  CALL IO_Field_write(TPFILE,'WUMF',XWUMF)
  CALL IO_Field_write(TPFILE,'WVMF',XWVMF)
END IF
!
IF( NRR > 1 .AND. CTURB /= 'NONE' ) THEN
  CALL IO_Field_write(TPFILE,'SRCT',XSRCT)
  CALL IO_Field_write(TPFILE,'SIGS',XSIGS)
END IF
!
!*       1.5    Reference state variables :
!
IF (LCOUPLES.AND.LOCEAN) THEN
  CALL IO_Field_write(TPFILE,'RHOREFZ',XRHODREFZO)
  CALL IO_Field_write(TPFILE,'THVREFZ',XTHVREFZO)
  CALL IO_Field_write(TPFILE,'EXNTOP', XEXNTOPO)
ELSE
  CALL IO_Field_write(TPFILE,'RHOREFZ',XRHODREFZ)
  CALL IO_Field_write(TPFILE,'THVREFZ',XTHVREFZ)
  CALL IO_Field_write(TPFILE,'EXNTOP', XEXNTOP)
END IF
!
!
!*       1.6  Tendencies                                         
!
IF (CPROGRAM == 'MESONH') THEN
  IF (CTEMP_SCHEME/='LEFR') THEN
    CALL IO_Field_write(TPFILE,'US_PRES',XRUS_PRES)
    CALL IO_Field_write(TPFILE,'VS_PRES',XRVS_PRES)
    CALL IO_Field_write(TPFILE,'WS_PRES',XRWS_PRES)
  END IF
  IF (LSPLIT_CFL) THEN
    CALL IO_Field_write(TPFILE,'THS_CLD',XRTHS_CLD)
!
    IF (NRR >=1) THEN
      IF (LUSERV) CALL IO_Field_write(TPFILE,'RVS_CLD',XRRS_CLD(:,:,:,IDX_RVT))
      IF (LUSERC) CALL IO_Field_write(TPFILE,'RCS_CLD',XRRS_CLD(:,:,:,IDX_RCT))
      IF (LUSERR) CALL IO_Field_write(TPFILE,'RRS_CLD',XRRS_CLD(:,:,:,IDX_RRT))
      IF (LUSERI) CALL IO_Field_write(TPFILE,'RIS_CLD',XRRS_CLD(:,:,:,IDX_RIT))
      IF (LUSERS) CALL IO_Field_write(TPFILE,'RSS_CLD',XRRS_CLD(:,:,:,IDX_RST))
      IF (LUSERG) CALL IO_Field_write(TPFILE,'RGS_CLD',XRRS_CLD(:,:,:,IDX_RGT))
      IF (LUSERH) CALL IO_Field_write(TPFILE,'RHS_CLD',XRRS_CLD(:,:,:,IDX_RHT))
    END IF 
  END IF
END IF 
!
!IF (LSPLIT_CFL) THEN
! IF (NSV >=1) THEN
!    DO JSV = NSV_C2R2BEG,NSV_C2R2END
!     IF (JSV == NSV_C2R2BEG ) THEN
!       TZFIELD = TFIELDMETADATA(       &
!         CMNHNAME   = 'RSVS_CLD1',     &
!         CSTDNAME   = '',              &
!         CLONGNAME  = 'RSVS_CLD1',     &
!         CUNITS     = '1',             &
!         CDIR       = 'XY',            &
!         CCOMMENT   = 'X_Y_Z_RHS_CLD', &
!         NGRID      = 1,               &
!         NTYPE      = TYPEREAL,        &
!         NDIMS      = 3,               &
!         LTIMEDEP   = .TRUE.           )
!       CALL IO_Field_write(TPFILE,TZFIELD,XRRS_CLD(:,:,:,IRR))
!     END IF
!     IF (JSV == NSV_C2R2END ) THEN
!       TZFIELD = TFIELDMETADATA(       &
!         CMNHNAME   = 'RSVS_CLD2',     &
!         CSTDNAME   = '',              &
!         CLONGNAME  = 'RSVS_CLD2',     &
!         CUNITS     = '1',             &
!         CDIR       = 'XY',            &
!         CCOMMENT   = 'X_Y_Z_RHS_CLD', &
!         NGRID      = 1,               &
!         NTYPE      = TYPEREAL,        &
!         NDIMS      = 3,               &
!         LTIMEDEP   = .TRUE.           )
!       CALL IO_Field_write(TPFILE,TZFIELD,XRRS_CLD(:,:,:,IRR))
!     END IF
!    END DO
! END IF
!ENDIF
!
!*       1.8    Diagnostic variables related to the radiations
!
!
IF (CRAD /= 'NONE') THEN
  CALL IO_Field_write(TPFILE,'DTRAD_FULL',TDTRAD_FULL)
  CALL IO_Field_write(TPFILE,'DTRAD_CLLY',TDTRAD_CLONLY)
!
  CALL IO_Field_write(TPFILE,'DTHRAD',      XDTHRAD)
  CALL IO_Field_write(TPFILE,'FLALWD',      XFLALWD)
  CALL IO_Field_write(TPFILE,'DIRFLASWD',   XDIRFLASWD)
  CALL IO_Field_write(TPFILE,'SCAFLASWD',   XSCAFLASWD)
  CALL IO_Field_write(TPFILE,'DIRSRFSWD',   XDIRSRFSWD)
  CALL IO_Field_write(TPFILE,'CLEARCOL_TM1',NCLEARCOL_TM1)
  CALL IO_Field_write(TPFILE,'ZENITH',      XZENITH)
  CALL IO_Field_write(TPFILE,'AZIM',        XAZIM)
  CALL IO_Field_write(TPFILE,'DIR_ALB',     XDIR_ALB)
  CALL IO_Field_write(TPFILE,'SCA_ALB',     XSCA_ALB)
  !
  CALL PRINT_MSG(NVERB_INFO,'IO','WRITE_LFIFM_n','EMIS: writing only first band')
  CALL FIND_FIELD_ID_FROM_MNHNAME('EMIS',IID,IRESP)
  TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
  TZFIELD%NDIMS = 2
  TZFIELD%NDIMLIST(3) = TZFIELD%NDIMLIST(4)
  TZFIELD%NDIMLIST(4) = NMNHDIM_UNUSED
  CALL IO_Field_write(TPFILE,TZFIELD,XEMIS(:,:,1))
  !
  CALL IO_Field_write(TPFILE,'TSRAD',       XTSRAD)
ENDIF
!
IF (NRR > 1 .AND. CPROGRAM == 'MESONH') THEN
  CALL IO_Field_write(TPFILE,'CLDFR',XCLDFR)
  CALL IO_Field_write(TPFILE,'ICEFR',XICEFR)
  CALL IO_Field_write(TPFILE,'RAINFR',XRAINFR)
END IF
!
!
!*       1.9     Diagnostic variables related to deep convection
!
!
IF (CDCONV /= 'NONE' .OR. CSCONV == 'KAFR') THEN
!
! 
!
  CALL IO_Field_write(TPFILE,'DTDCONV',  TDTDCONV)
  CALL IO_Field_write(TPFILE,'COUNTCONV',NCOUNTCONV)
  CALL IO_Field_write(TPFILE,'DTHCONV',  XDTHCONV)
  CALL IO_Field_write(TPFILE,'DRVCONV',  XDRVCONV)
  CALL IO_Field_write(TPFILE,'DRCCONV',  XDRCCONV)
  CALL IO_Field_write(TPFILE,'DRICONV',  XDRICONV)
!
  CALL FIND_FIELD_ID_FROM_MNHNAME('PRCONV',IID,IRESP)
  TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
  TZFIELD%CUNITS = 'mm hour-1'
  CALL IO_Field_write(TPFILE,TZFIELD,XPRCONV*3.6E6)
!
  CALL FIND_FIELD_ID_FROM_MNHNAME('PACCONV',IID,IRESP)
  TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
  TZFIELD%CUNITS = 'mm'
  CALL IO_Field_write(TPFILE,TZFIELD,XPACCONV*1.0E3)
!
  CALL FIND_FIELD_ID_FROM_MNHNAME('PRSCONV',IID,IRESP)
  TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
  TZFIELD%CUNITS = 'mm hour-1'
  CALL IO_Field_write(TPFILE,TZFIELD,XPRSCONV*3.6E6)
!
  IF ( LCH_CONV_LINOX ) THEN
    CALL IO_Field_write(TPFILE,'IC_RATE',    XIC_RATE)
    CALL IO_Field_write(TPFILE,'CG_RATE',    XCG_RATE)
    CALL IO_Field_write(TPFILE,'IC_TOTAL_NB',XIC_TOTAL_NUMBER)
    CALL IO_Field_write(TPFILE,'CG_TOTAL_NB',XCG_TOTAL_NUMBER)
  END IF
!
  IF ( LCHTRANS .AND. NSV > 0 ) THEN
   ! scalar variables are recorded
   ! individually in the file
    TZFIELD = TFIELDMETADATA(     &
      CMNHNAME   = 'generic for DSVCONV', & !Temporary name to ease identification
      CUNITS     = 's-1',         &
      CDIR       = 'XY',          &
      NGRID      = 1,             &
      NTYPE      = TYPEREAL,      &
      NDIMS      = 3,             &
      LTIMEDEP   = .TRUE.         )

    DO JSV = 1, NSV
      TZFIELD%CMNHNAME   = 'DSVCONV_' // TRIM( TSVLIST(JSV)%CMNHNAME )
      TZFIELD%CLONGNAME  = 'DSVCONV_' // TRIM( TSVLIST(JSV)%CLONGNAME )
      TZFIELD%CCOMMENT   = 'Convective tendency for ' // TRIM( TSVLIST(JSV)%CMNHNAME )
      CALL IO_Field_write( TPFILE, TZFIELD, XDSVCONV(:,:,:,JSV) )
    END DO
  END IF
!
END IF
!
!
!*       1.10   Diagnostic variables related to the precipitations
!
IF ( (CPROGRAM /= 'IDEAL') .AND. (CPROGRAM /= 'REAL')  )THEN
  IF (ASSOCIATED(XINPRC)) THEN
  IF (SIZE(XINPRC) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRC',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XINPRC*3.6E6)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRC',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_write(TPFILE,TZFIELD,XACPRC*1.0E3)
!
  ENDIF
  ENDIF
!
  IF (ASSOCIATED(XINDEP)) THEN
  IF (SIZE(XINDEP) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INDEP',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XINDEP*3.6E6)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACDEP',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_write(TPFILE,TZFIELD,XACDEP*1.0E3)
!
  ENDIF
  ENDIF
!
  IF (ASSOCIATED(XINPRR)) THEN
  IF (SIZE(XINPRR) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRR',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XINPRR*3.6E6)
!
    CALL IO_Field_write(TPFILE,'INPRR3D',XINPRR3D)
    CALL IO_Field_write(TPFILE,'EVAP3D', XEVAP3D)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRR',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_write(TPFILE,TZFIELD,XACPRR*1.0E3)
!
  ENDIF
  ENDIF
!
  IF (ASSOCIATED(XINPRS)) THEN
  IF (SIZE(XINPRS) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRS',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XINPRS*3.6E6)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRS',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_write(TPFILE,TZFIELD,XACPRS*1.0E3)
  END IF
  END IF
!
  IF (ASSOCIATED(XINPRG)) THEN
  IF (SIZE(XINPRG) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRG',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XINPRG*3.6E6)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRG',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_write(TPFILE,TZFIELD,XACPRG*1.0E3)
  END IF
  END IF
!
  IF (ASSOCIATED(XINPRH)) THEN
  IF (SIZE(XINPRH) /= 0 ) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRH',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XINPRH*3.6E6)
!
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRH',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_write(TPFILE,TZFIELD,XACPRH*1.0E3)
  ENDIF
  ENDIF
!
  IF (ASSOCIATED(XINPRS)) THEN
  IF (SIZE(XINPRS) /= 0 ) THEN
    ZWORK2D = XINPRR + XINPRS
    IF (SIZE(XINPRG) /= 0 ) ZWORK2D = ZWORK2D + XINPRG
    IF (SIZE(XINPRH) /= 0 ) ZWORK2D = ZWORK2D + XINPRH
    IF (SIZE(XINPRC) /= 0 ) ZWORK2D = ZWORK2D + XINPRC
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRT',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK2D*3.6E6)
!
    ZWORK2D = XACPRR + XACPRS
    IF (SIZE(XINPRG) /= 0 ) ZWORK2D = ZWORK2D + XACPRG
    IF (SIZE(XINPRH) /= 0 ) ZWORK2D = ZWORK2D + XACPRH
    IF (SIZE(XINPRC) /= 0 ) ZWORK2D = ZWORK2D + XACPRC
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRT',IID,IRESP)
    TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK2D*1.0E3)
  END IF
  END IF
!
END IF
!
IF(LBLOWSNOW) THEN
  IF (ASSOCIATED(XSNWSUBL3D)) THEN
    IF (SIZE(XSNWSUBL3D) /= 0 ) THEN
      TZFIELD = TFIELDMETADATA(                                             &
        CMNHNAME   = 'SNWSUBL3D',                                           &
        CSTDNAME   = '',                                                    &
        CLONGNAME  = 'SNWSUBL3D',                                           &
        CUNITS     = 'kg m-3 s-1',                                          &
        CDIR       = 'XY',                                                  &
        CCOMMENT   = 'X_Y_INstantaneous 3D Drifting snow sublimation flux', &
        NGRID      = 1,                                                     &
        NTYPE      = TYPEREAL,                                              &
        NDIMS      = 3,                                                     &
        LTIMEDEP   = .TRUE.                                                 )
      CALL IO_Field_write(TPFILE,TZFIELD,XSNWSUBL3D(:,:,:))
      ZWORK2D(:,:) = 0.
      DO JK = IKB,IKE
        ZWORK2D(:,:) = ZWORK2D(:,:)+XSNWSUBL3D(:,:,JK) * &
                    (XZZ(:,:,JK+1)-XZZ(:,:,JK))/XRHOLW*3600*24
      END DO
      ZWORK2D(:,:) = ZWORK2D(:,:)*1000. ! vapor water in mm unit
      !
      TZFIELD = TFIELDMETADATA(                                 &
        CMNHNAME   = 'COL_SNWSUBL',                             &
        CSTDNAME   = '',                                        &
        CLONGNAME  = 'COL_SNWSUBL',                             &
        CUNITS     = 'mm day-1',                                &
        CDIR       = 'XY',                                      &
        CCOMMENT   = 'X_Y_Column Sublimation Rate (mmSWE/day)', &
        NGRID      = 4,                                         &
        NTYPE      = TYPEREAL,                                  &
        NDIMS      = 2,                                         &
        LTIMEDEP   = .TRUE.                                     )
      CALL IO_Field_write(TPFILE,TZFIELD,ZWORK2D(:,:))
    END IF
  END IF
ENDIF
!
!*       1.11   Ocean LES variables
!
IF ((.NOT.LCOUPLES).AND.LOCEAN) THEN
  CALL IO_Field_write(TPFILE,'NFRCLT',NFRCLT)
  CALL IO_Field_write(TPFILE,'NINFRT',NINFRT)
  !
  TZFIELD = TFIELDMETADATA(                               &
    CMNHNAME   = 'SSUFL_T',                               &
    CSTDNAME   = '',                                      &
    CLONGNAME  = 'SSUFL',                                 &
    CUNITS     = 'kg m-1 s-1',                            &
    CDIR       = '--',                                    &
    CCOMMENT   = 'sfc stress along U to force ocean LES', &
    NGRID      = 0,                                       &
    NTYPE      = TYPEREAL,                                &
    NDIMS      = 1,                                       &
    LTIMEDEP   = .FALSE.                                  )
  CALL IO_Field_write(TPFILE,TZFIELD,XSSUFL_T(:))
  !
  TZFIELD = TFIELDMETADATA(                               &
    CMNHNAME   = 'SSVFL_T',                               &
    CSTDNAME   = '',                                      &
    CLONGNAME  = 'SSVFL',                                 &
    CUNITS     = 'kg m-1 s-1',                            &
    CDIR       = '--',                                    &
    CCOMMENT   = 'sfc stress along V to force ocean LES', &
    NGRID      = 0,                                       &
    NTYPE      = TYPEREAL,                                &
    NDIMS      = 1,                                       &
    LTIMEDEP   = .FALSE.                                  )
  CALL IO_Field_write(TPFILE,TZFIELD,XSSVFL_T(:))
  !
  TZFIELD = TFIELDMETADATA(                                &
    CMNHNAME   = 'SSTFL_T',                                &
    CSTDNAME   = '',                                       &
    CLONGNAME  = 'SSTFL',                                  &
    CUNITS     = 'kg m3 K m s-1',                          &
    CDIR       = '--',                                     &
    CCOMMENT   = 'sfc total heat flux to force ocean LES', &
    NGRID      = 0,                                        &
    NTYPE      = TYPEREAL,                                 &
    NDIMS      = 1,                                        &
    LTIMEDEP   = .FALSE.                                   )
  CALL IO_Field_write(TPFILE,TZFIELD,XSSTFL_T(:))
  !
  TZFIELD = TFIELDMETADATA(                           &
    CMNHNAME   = 'SSOLA_T',                           &
    CSTDNAME   = '',                                  &
    CLONGNAME  = 'SSOLA',                             &
    CUNITS     = 'kg m3 K m s-1',                     &
    CDIR       = '--',                                &
    CCOMMENT   = 'sfc solar flux to force ocean LES', &
    NGRID      = 0,                                   &
    NTYPE      = TYPEREAL,                            &
    NDIMS      = 1,                                   &
    LTIMEDEP   = .FALSE.                              )
  CALL IO_Field_write(TPFILE,TZFIELD,XSSOLA_T(:))
  !
END IF ! ocean sfc forcing end    
!
!*       1.12   Forcing variables
!
IF (LFORCING) THEN
!
  CALL IO_Field_write(TPFILE,'FRC',NFRC)
!
  DO JT=1,NFRC
!
    WRITE (YFRC,'(I3.3)') JT
!
    TZFIELD = TFIELDMETADATA(                             &
      CMNHNAME   = 'DTFRC'//YFRC,                         &
      CSTDNAME   = '',                                    &
      CLONGNAME  = 'DTFRC'//YFRC,                         &
      CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S', &
      CDIR       = '--',                                  &
      CCOMMENT   = 'Date of forcing profile '//YFRC,      &
      NGRID      = 0,                                     &
      NTYPE      = TYPEDATE,                              &
      NDIMS      = 0,                                     &
      LTIMEDEP   = .FALSE.                                )
    CALL IO_Field_write(TPFILE,TZFIELD,TDTFRC(JT))
!
    TZFIELD = TFIELDMETADATA(                                    &
      CMNHNAME   = 'UFRC'//YFRC,                                 &
      CSTDNAME   = '',                                           &
      CLONGNAME  = 'UFRC'//YFRC,                                 &
      CUNITS     = 'm s-1',                                      &
      CDIR       = '--',                                         &
      CCOMMENT   = 'Zonal component of horizontal forcing wind', &
      NGRID      = 1,                                            &
      NTYPE      = TYPEREAL,                                     &
      NDIMS      = 1,                                            &
      LTIMEDEP   = .FALSE.                                       )
    CALL IO_Field_write(TPFILE,TZFIELD,XUFRC(:,JT))
!
    TZFIELD = TFIELDMETADATA(                                       &
      CMNHNAME   = 'VFRC'//YFRC,                                    &
      CSTDNAME   = '',                                              &
      CLONGNAME  = 'VFRC'//YFRC,                                    &
      CUNITS     = 'm s-1',                                         &
      CDIR       = '--',                                            &
      CCOMMENT   = 'Meridian component of horizontal forcing wind', &
      NGRID      = 1,                                               &
      NTYPE      = TYPEREAL,                                        &
      NDIMS      = 1,                                               &
      LTIMEDEP   = .FALSE.                                          )
    CALL IO_Field_write(TPFILE,TZFIELD,XVFRC(:,JT))
!
    TZFIELD = TFIELDMETADATA(               &
      CMNHNAME   = 'WFRC'//YFRC,            &
      CSTDNAME   = '',                      &
      CLONGNAME  = 'WFRC'//YFRC,            &
      CUNITS     = 'm s-1',                 &
      CDIR       = '--',                    &
      CCOMMENT   = 'Vertical forcing wind', &
      NGRID      = 4,                       &
      NTYPE      = TYPEREAL,                &
      NDIMS      = 1,                       &
      LTIMEDEP   = .FALSE.                  )
    CALL IO_Field_write(TPFILE,TZFIELD,XWFRC(:,JT))
!
    TZFIELD = TFIELDMETADATA(                       &
      CMNHNAME   = 'THFRC'//YFRC,                   &
      CSTDNAME   = '',                              &
      CLONGNAME  = 'THFRC'//YFRC,                   &
      CUNITS     = 'K',                             &
      CDIR       = '--',                            &
      CCOMMENT   = 'Forcing potential temperature', &
      NGRID      = 1,                               &
      NTYPE      = TYPEREAL,                        &
      NDIMS      = 1,                               &
      LTIMEDEP   = .FALSE.                          )
    CALL IO_Field_write(TPFILE,TZFIELD,XTHFRC(:,JT))
!
    TZFIELD = TFIELDMETADATA(                    &
      CMNHNAME   = 'RVFRC'//YFRC,                &
      CSTDNAME   = '',                           &
      CLONGNAME  = 'RVFRC'//YFRC,                &
      CUNITS     = 'kg kg-1',                    &
      CDIR       = '--',                         &
      CCOMMENT   = 'Forcing vapor mixing ratio', &
      NGRID      = 1,                            &
      NTYPE      = TYPEREAL,                     &
      NDIMS      = 1,                            &
      LTIMEDEP   = .FALSE.                       )
    CALL IO_Field_write(TPFILE,TZFIELD,XRVFRC(:,JT))
!
    TZFIELD = TFIELDMETADATA(                                                &
      CMNHNAME   = 'TENDTHFRC'//YFRC,                                        &
      CSTDNAME   = '',                                                       &
      CLONGNAME  = 'TENDTHFRC'//YFRC,                                        &
      CUNITS     = 'K s-1',                                                  &
      CDIR       = '--',                                                     &
      CCOMMENT   = 'Large-scale potential temperature tendency for forcing', &
      NGRID      = 1,                                                        &
      NTYPE      = TYPEREAL,                                                 &
      NDIMS      = 1,                                                        &
      LTIMEDEP   = .FALSE.                                                   )
    CALL IO_Field_write(TPFILE,TZFIELD,XTENDTHFRC(:,JT))
!
    TZFIELD = TFIELDMETADATA(                                             &
      CMNHNAME   = 'TENDRVFRC'//YFRC,                                     &
      CSTDNAME   = '',                                                    &
      CLONGNAME  = 'TENDRVFRC'//YFRC,                                     &
      CUNITS     = 'kg kg-1 s-1',                                         &
      CDIR       = '--',                                                  &
      CCOMMENT   = 'Large-scale vapor mixing ratio tendency for forcing', &
      NGRID      = 1,                                                     &
      NTYPE      = TYPEREAL,                                              &
      NDIMS      = 1,                                                     &
      LTIMEDEP   = .FALSE.                                                )
    CALL IO_Field_write(TPFILE,TZFIELD,XTENDRVFRC(:,JT))
!
    TZFIELD = TFIELDMETADATA(                                                &
      CMNHNAME   = 'GXTHFRC'//YFRC,                                          &
      CSTDNAME   = '',                                                       &
      CLONGNAME  = 'GXTHFRC'//YFRC,                                          &
      CUNITS     = 'K m-1',                                                  &
      CDIR       = '--',                                                     &
      CCOMMENT   = 'Large-scale potential temperature gradient for forcing', &
      NGRID      = 1,                                                        &
      NTYPE      = TYPEREAL,                                                 &
      NDIMS      = 1,                                                        &
      LTIMEDEP   = .FALSE.                                                   )
    CALL IO_Field_write(TPFILE,TZFIELD,XGXTHFRC(:,JT))
!
    TZFIELD = TFIELDMETADATA(                                                &
      CMNHNAME   = 'GYTHFRC'//YFRC,                                          &
      CSTDNAME   = '',                                                       &
      CLONGNAME  = 'GYTHFRC'//YFRC,                                          &
      CUNITS     = 'K m-1',                                                  &
      CDIR       = '--',                                                     &
      CCOMMENT   = 'Large-scale potential temperature gradient for forcing', &
      NGRID      = 1,                                                        &
      NTYPE      = TYPEREAL,                                                 &
      NDIMS      = 1,                                                        &
      LTIMEDEP   = .FALSE.                                                   )
    CALL IO_Field_write(TPFILE,TZFIELD,XGYTHFRC(:,JT))
!
    TZFIELD = TFIELDMETADATA(                 &
      CMNHNAME   = 'PGROUNDFRC'//YFRC,        &
      CSTDNAME   = '',                        &
      CLONGNAME  = 'PGROUNDFRC'//YFRC,        &
      CUNITS     = 'Pa',                      &
      CDIR       = '--',                      &
      CCOMMENT   = 'Forcing ground pressure', &
      NGRID      = 1,                         &
      NTYPE      = TYPEREAL,                  &
      NDIMS      = 0,                         &
      LTIMEDEP   = .FALSE.                    )
    CALL IO_Field_write(TPFILE,TZFIELD,XPGROUNDFRC(JT))
!
    TZFIELD = TFIELDMETADATA(                            &
      CMNHNAME   = 'TENDUFRC'//YFRC,                     &
      CSTDNAME   = '',                                   &
      CLONGNAME  = 'TENDUFRC'//YFRC,                     &
      CUNITS     = 'm s-1',                              &
      CDIR       = '--',                                 &
      CCOMMENT   = 'Large-scale U tendency for forcing', &
      NGRID      = 1,                                    &
      NTYPE      = TYPEREAL,                             &
      NDIMS      = 1,                                    &
      LTIMEDEP   = .FALSE.                               )
    CALL IO_Field_write(TPFILE,TZFIELD,XTENDUFRC(:,JT))
!
    TZFIELD = TFIELDMETADATA(                            &
      CMNHNAME   = 'TENDVFRC'//YFRC,                     &
      CSTDNAME   = '',                                   &
      CLONGNAME  = 'TENDVFRC'//YFRC,                     &
      CUNITS     = 'm s-1',                              &
      CDIR       = '--',                                 &
      CCOMMENT   = 'Large-scale V tendency for forcing', &
      NGRID      = 1,                                    &
      NTYPE      = TYPEREAL,                             &
      NDIMS      = 1,                                    &
      LTIMEDEP   = .FALSE.                               )
    CALL IO_Field_write(TPFILE,TZFIELD,XTENDVFRC(:,JT))
!
  END DO
!
!
END IF
!
! -------------------------------------------------------------------------
IF ( L2D_ADV_FRC ) THEN
!
  TZFIELD = TFIELDMETADATA(                    &
    CMNHNAME   = 'NADVFRC1',                   &
    CSTDNAME   = '',                           &
    CLONGNAME  = 'NADVFRC1',                   &
    CUNITS     = '1',                          &
    CDIR       = '--',                         &
    CCOMMENT   = 'Number of forcing profiles', &
    NGRID      = 0,                            &
    NTYPE      = TYPEINT,                      &
    NDIMS      = 0,                            &
    LTIMEDEP   = .FALSE.                       )
  CALL IO_Field_write(TPFILE,TZFIELD,NADVFRC)
!
  DO JT=1,NADVFRC
!
    WRITE (YFRC,'(I3.3)') JT
!
    TZFIELD = TFIELDMETADATA(                                       &
      CMNHNAME   = 'DTADV'//YFRC,                                   &
      CSTDNAME   = '',                                              &
      CLONGNAME  = 'DTADV'//YFRC,                                   &
      CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S',           &
      CDIR       = '--',                                            &
      CCOMMENT   = 'Date and time of the advecting forcing '//YFRC, &
      NGRID      = 0,                                               &
      NTYPE      = TYPEDATE,                                        &
      NDIMS      = 0,                                               &
      LTIMEDEP   = .FALSE.                                          )
    CALL IO_Field_write(TPFILE,TZFIELD,TDTADVFRC(JT))
!
    TZFIELD = TFIELDMETADATA(      &
      CMNHNAME   = 'TH_ADV'//YFRC, &
      CSTDNAME   = '',             &
      CLONGNAME  = 'TH_ADV'//YFRC, &
      CUNITS     = 'K s-1',        &
      CDIR       = '--',           &
      CCOMMENT   = '',             &
      NGRID      = 1,              &
      NTYPE      = TYPEREAL,       &
      NDIMS      = 3,              &
      LTIMEDEP   = .FALSE.         )
    CALL IO_Field_write(TPFILE,TZFIELD,XDTHFRC(:,:,:,JT))
!    
    TZFIELD = TFIELDMETADATA(     &
      CMNHNAME   = 'Q_ADV'//YFRC, &
      CSTDNAME   = '',            &
      CLONGNAME  = 'Q_ADV'//YFRC, &
      CUNITS     = 'kg kg-1 s-1', &
      CDIR       = '--',          &
      CCOMMENT   = '',            &
      NGRID      = 1,             &
      NTYPE      = TYPEREAL,      &
      NDIMS      = 3,             &
      LTIMEDEP   = .FALSE.        )
    CALL IO_Field_write(TPFILE,TZFIELD,XDRVFRC(:,:,:,JT))
!
  ENDDO
ENDIF
!
IF ( L2D_REL_FRC ) THEN
!
  TZFIELD = TFIELDMETADATA(                    &
    CMNHNAME   = 'NRELFRC1',                   &
    CSTDNAME   = '',                           &
    CLONGNAME  = 'NRELFRC1',                   &
    CUNITS     = '1',                          &
    CDIR       = '--',                         &
    CCOMMENT   = 'Number of forcing profiles', &
    NGRID      = 0,                            &
    NTYPE      = TYPEINT,                      &
    NDIMS      = 0,                            &
    LTIMEDEP   = .FALSE.                       )
  CALL IO_Field_write(TPFILE,TZFIELD,NRELFRC)
!
  DO JT=1,NRELFRC
!
    WRITE (YFRC,'(I3.3)') JT
!
    TZFIELD = TFIELDMETADATA(                                        &
      CMNHNAME   = 'DTREL'//YFRC,                                    &
      CSTDNAME   = '',                                               &
      CLONGNAME  = 'DTREL'//YFRC,                                    &
      CUNITS     = 'seconds since YYYY-MM-DD HH:MM:SS.S',            &
      CDIR       = '--',                                             &
      CCOMMENT   = 'Date and time of the relaxation forcing '//YFRC, &
      NGRID      = 0,                                                &
      NTYPE      = TYPEDATE,                                         &
      NDIMS      = 0,                                                &
      LTIMEDEP   = .FALSE.                                           )
    CALL IO_Field_write(TPFILE,TZFIELD,TDTRELFRC(JT))
!                                                                
    TZFIELD = TFIELDMETADATA(      &
      CMNHNAME   = 'TH_REL'//YFRC, &
      CSTDNAME   = '',             &
      CLONGNAME  = 'TH_REL'//YFRC, &
      CUNITS     = 'K',            &
      CDIR       = '--',           &
      CCOMMENT   = '',             &
      NGRID      = 1,              &
      NTYPE      = TYPEREAL,       &
      NDIMS      = 3,              &
      LTIMEDEP   = .FALSE.         )
    CALL IO_Field_write(TPFILE,TZFIELD,XTHREL(:,:,:,JT))
!    
    TZFIELD = TFIELDMETADATA(     &
      CMNHNAME   = 'Q_REL'//YFRC, &
      CSTDNAME   = '',            &
      CLONGNAME  = 'Q_REL'//YFRC, &
      CUNITS     = 'kg kg-1',     &
      CDIR       = '--',          &
      CCOMMENT   = '',            &
      NGRID      = 1,             &
      NTYPE      = TYPEREAL,      &
      NDIMS      = 3,             &
      LTIMEDEP   = .FALSE.        )
    CALL IO_Field_write(TPFILE,TZFIELD,XRVREL(:,:,:,JT))
!
  ENDDO
ENDIF
!
!*       1.13   Eddy Fluxes variables    ! Modif PP
!
IF ( LTH_FLX ) THEN
   CALL IO_Field_write(TPFILE,'VT_FLX',XVTH_FLUX_M)
   CALL IO_Field_write(TPFILE,'WT_FLX',XWTH_FLUX_M)
END IF
!
IF ( LUV_FLX) CALL IO_Field_write(TPFILE,'VU_FLX',XVU_FLUX_M)
!
!*       1.14   Balloon variables
!
!
! Write balloon coordinates in backup file to allow restart with current balloon position
IF (LFLYER) CALL WRITE_BALLOON_n(TPFILE)
!
!
!*       1.15    Filtered variables for hurricane initialization
!
!
IF ( CPROGRAM=='REAL  ' ) THEN
  IF (LFILTERING) THEN
  !
    IF (NDIAG_FILT >=0) THEN
!
!             i) Total fields (TOT=BASIC+TOTDIS)
!
      CALL IO_Field_write(TPFILE,'UT15',   XUTOT)
      CALL IO_Field_write(TPFILE,'VT15',   XVTOT)
      CALL IO_Field_write(TPFILE,'TEMPTOT',XTTOT)
      IF (INDEX(CFILTERING,'P')/=0) CALL IO_Field_write(TPFILE,'PRESTOT',XPTOT)
      IF (INDEX(CFILTERING,'Q')/=0) CALL IO_Field_write(TPFILE,'HUMTOT', XQTOT)
!
!             ii) Environmental fields (ENV=TOT-VORDIS)
!
      CALL IO_Field_write(TPFILE,'UT16',   XUENV)
      CALL IO_Field_write(TPFILE,'VT16',   XVENV)
      CALL IO_Field_write(TPFILE,'TEMPENV',XTENV)
      IF (INDEX(CFILTERING,'P')/=0) CALL IO_Field_write(TPFILE,'PRESENV',XPENV)
      IF (INDEX(CFILTERING,'Q')/=0) CALL IO_Field_write(TPFILE,'HUMENV', XQENV)
!
    END IF
    IF (NDIAG_FILT >=1) THEN
!
!             iii) Basic (filtered) fields
!
      CALL IO_Field_write(TPFILE,'UT17',   XUBASIC)
      CALL IO_Field_write(TPFILE,'VT17',   XVBASIC)
      CALL IO_Field_write(TPFILE,'TEMPBAS',XTBASIC)
      IF (INDEX(CFILTERING,'P')/=0) CALL IO_Field_write(TPFILE,'PRESBAS',XPBASIC)
      IF (INDEX(CFILTERING,'Q')/=0) CALL IO_Field_write(TPFILE,'HUMBAS', XQBASIC)
    END IF
    IF (NDIAG_FILT >=2) THEN
!
!             iv) Total disturbance tangential wind
!
      CALL IO_Field_write(TPFILE,'VTDIS',XVTDIS)
!
    END IF
!
  END IF
!
!*       1.16    Dummy variables in PREP_REAL_CASE
!
  IF (ALLOCATED(CDUMMY_2D)) THEN
    TZFIELD = TFIELDMETADATA(                       &
      CMNHNAME   = 'generic for CDUMMY_2D variables', & !Temporary name to ease identification
      CSTDNAME   = '',                                &
      CUNITS     = '',                                &
      CDIR       = 'XY',                              &
      NGRID      = 1,                                 &
      NTYPE      = TYPEREAL,                          &
      NDIMS      = 2,                                 &
      LTIMEDEP   = .TRUE.                             )
    !
    DO JI = 1, SIZE( XDUMMY_2D, 3 )
      TZFIELD%CMNHNAME   = ADJUSTL(CDUMMY_2D(JI))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
      CALL IO_Field_write(TPFILE,TZFIELD,XDUMMY_2D(:,:,JI))
    END DO
  END IF
!
END IF
!
!*       1.17    Wind turbine variables 
!
!             i) Main
!
IF (LMAIN_EOL .AND. IMI == NMODEL_EOL) THEN
  TZFIELD = TFIELDMETADATA(                            &
    CMNHNAME   = 'generic for wind turbine var',       & !Temporary name to ease identification
    CSTDNAME   = '',                                   &
    CUNITS     = 'N',                                  &
    CDIR       = 'XY',                                 &
    NGRID      = 1,                                    &
    NTYPE      = TYPEREAL,                             &
    NDIMS      = 3,                                    &
    LTIMEDEP   = .TRUE.                                )
!
  TZFIELD%CMNHNAME   = 'FX_RG'
  TZFIELD%CLONGNAME  = 'FX_RG'
  TZFIELD%CCOMMENT   = 'X-component field of aerodynamic force (wind->rotor) in global frame (N)'
  CALL IO_Field_write(TPFILE,TZFIELD,XFX_RG)
!
  TZFIELD%CMNHNAME   = 'FY_RG'
  TZFIELD%CLONGNAME  = 'FY_RG'
  TZFIELD%CCOMMENT   = 'Y-component field of aerodynamic force (wind->rotor) in global frame (N)'
  CALL IO_Field_write(TPFILE,TZFIELD,XFY_RG)
!
  TZFIELD%CMNHNAME   = 'FZ_RG'
  TZFIELD%CLONGNAME  = 'FZ_RG'
  TZFIELD%CCOMMENT   = 'Z-component field of aerodynamic force (wind->rotor) in global frame (N)'
  CALL IO_Field_write(TPFILE,TZFIELD,XFZ_RG)
!
  TZFIELD%CMNHNAME   = 'FX_SMR_RG'
  TZFIELD%CLONGNAME  = 'FX_SMR_RG'
  TZFIELD%CCOMMENT   = 'X-component field of smeared aerodynamic force (wind->rotor) in global frame (N)'
  TZFIELD%CCOMMENT   = ''
  CALL IO_Field_write(TPFILE,TZFIELD,XFX_SMR_RG)
!
  TZFIELD%CMNHNAME   = 'FY_SMR_RG'
  TZFIELD%CLONGNAME  = 'FY_SMR_RG'
  TZFIELD%CCOMMENT   = 'Y-component field of smeared aerodynamic force (wind->rotor) in global frame (N)'
  CALL IO_Field_write(TPFILE,TZFIELD,XFY_SMR_RG)
!
  TZFIELD%CMNHNAME   = 'FZ_SMR_RG'
  TZFIELD%CLONGNAME  = 'FZ_SMR_RG'
  TZFIELD%CCOMMENT   = 'Z-component field of smeared aerodynamic force (wind->rotor) in global frame (N)'
  CALL IO_Field_write(TPFILE,TZFIELD,XFZ_SMR_RG)
!
SELECT CASE(CMETH_EOL)
!
!             ii) Actuator Disk without Rotation model
!
  CASE('ADNR') ! Actuator Disc Non-Rotating
!
    TZFIELD = TFIELDMETADATA(                    &
      CMNHNAME   = 'generic for ADNR variables', & !Temporary name to ease identification
      CSTDNAME   = '',                           &
      CUNITS     = '1',                          &
      CDIR       = '--',                         &
      NGRID      = 1,                            &
      NTYPE      = TYPEREAL,                     &
      NDIMS      = 1,                            &
      LTIMEDEP   = .TRUE.                        )
!
    TZFIELD%CMNHNAME   = 'A_INDU'
    TZFIELD%CLONGNAME  = 'INDUCTION_FACTOR'
    TZFIELD%CCOMMENT   = 'Induction factor (1)'
    CALL IO_Field_write(TPFILE,TZFIELD,XA_INDU)
!
    TZFIELD%CMNHNAME   = 'CT_D'
    TZFIELD%CLONGNAME  = 'CTHRUST_D'
    TZFIELD%CCOMMENT   = 'Thrust coefficient at disk (1),    &
                          used with wind speed at disk'
    CALL IO_Field_write(TPFILE,TZFIELD,XCT_D)
!
    TZFIELD%CMNHNAME   = 'THRUT'
    TZFIELD%CLONGNAME  = 'THRUSTT_EOL'
    TZFIELD%CUNITS     = 'N'
    TZFIELD%CCOMMENT   = 'RID instantaneous thrust of the wind turbines (N)'
    CALL IO_Field_write(TPFILE,TZFIELD,XTHRUT)
!
    IF (MEAN_COUNT /= 0) THEN

      TZFIELD%CMNHNAME   = 'THRUMME'
      TZFIELD%CLONGNAME  = 'MEAN_THRUST_EOL'
      TZFIELD%CUNITS     = 'N'
      TZFIELD%CCOMMENT   = 'RID mean thrust of the wind turbines (N)'
      CALL IO_Field_write(TPFILE,TZFIELD,XTHRU_SUM/MEAN_COUNT)
!
    END IF
!
!             iii) Actuator Disc with Rotation Model
!
  CASE('ADR') ! Actuator Disc with Rotation
!
! * 1D Variables (rotor id)
    TZFIELD = TFIELDMETADATA(                   &
      CMNHNAME   = '1D ADR var: (rot)',         & !Temporary name to ease identification
      CSTDNAME   = '',                          &
      CDIR       = '--',                        &
      NGRID      = 1,                           &
      NTYPE      = TYPEREAL,                    &
      NDIMS      = 1,                           &
      LTIMEDEP   = .TRUE.                       ) 
!
    TZFIELD%CMNHNAME   = 'THRUT'
    TZFIELD%CLONGNAME  = 'THRUSTT_EOL'
    TZFIELD%CUNITS     = 'N'
    TZFIELD%CCOMMENT   = 'RID instantaneous thrust (N) of wind turbines'
    CALL IO_Field_write(TPFILE,TZFIELD,XTHRUT)
!
    TZFIELD%CMNHNAME   = 'TORQT'
    TZFIELD%CLONGNAME  = 'TORQUET_EOL'
    TZFIELD%CUNITS     = 'Nm'
    TZFIELD%CCOMMENT   = 'RID instantaneous torque (Nm) of wind turbines'
    CALL IO_Field_write(TPFILE,TZFIELD,XTORQT)
!
    TZFIELD%CMNHNAME   = 'POWT'
    TZFIELD%CLONGNAME  = 'POWERT_EOL'
    TZFIELD%CUNITS     = 'W'
    TZFIELD%CCOMMENT   = 'RID instantaneous power (W) of wind turbines'
    CALL IO_Field_write(TPFILE,TZFIELD,XPOWT)
!
!
! * 3D Variables (rotor id, azimuthal id, radial id)
    TZFIELD = TFIELDMETADATA(                           &
      CMNHNAME   = '3D ADR var: (rot,azi,rad)',         & 
      CSTDNAME   = '',                                  &
      CDIR       = '--',                                &
      NGRID      = 1,                                   &
      NTYPE      = TYPEREAL,                            &
      NDIMS      = 3,                                   &
      LTIMEDEP   = .TRUE.                               )
!
    TZFIELD%CMNHNAME   = 'ELT_RAD'
    TZFIELD%CLONGNAME  = 'ELT_RAD'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CCOMMENT   = 'RID_AZI_RAD radius (m) of wind turbine blade elements'
    CALL IO_Field_write(TPFILE,TZFIELD,XELT_RAD)
!
    TZFIELD%CMNHNAME   = 'ELT_AZI'
    TZFIELD%CLONGNAME  = 'ELT_AZI'
    TZFIELD%CUNITS     = 'rad'
    TZFIELD%CCOMMENT   = 'RID_AZI_RAD Azimutal angle (rad) of wind turbine blade elements'
    CALL IO_Field_write(TPFILE,TZFIELD,XELT_AZI)
!
    TZFIELD%CMNHNAME   = 'AOA'
    TZFIELD%CLONGNAME  = 'ANGLE OF ATTACK'
    TZFIELD%CUNITS     = 'rad'
    TZFIELD%CCOMMENT   = 'RID_AZI_RAD instantaneous angle of attack (rad)'
    CALL IO_Field_write(TPFILE,TZFIELD,XAOA_GLB)
!
    TZFIELD%CMNHNAME   = 'FLIFT'
    TZFIELD%CLONGNAME  = 'LIFT FORCE'
    TZFIELD%CUNITS     = 'N'
    TZFIELD%CCOMMENT   = 'RID_AZI_RAD instantaneous lift (N) in relative frame'
    CALL IO_Field_write(TPFILE,TZFIELD,XFLIFT_GLB)
!
    TZFIELD%CMNHNAME   = 'FDRAG'
    TZFIELD%CLONGNAME  = 'DRAG FORCE'
    TZFIELD%CUNITS     = 'N'
    TZFIELD%CCOMMENT   = 'RID_AZI_RAD instantaneous drag (N) in relative frame'
    CALL IO_Field_write(TPFILE,TZFIELD,XFDRAG_GLB)
!
! * 4D Variables (rotor id, azimuthal id, radial id, xyz)
    TZFIELD = TFIELDMETADATA(                               &
      CMNHNAME   = '4D ADR var: (rot,azi,rad,xyz)',         & 
      CSTDNAME   = '',                                      &
      CDIR       = '--',                                    &
      NGRID      = 1,                                       &
      NTYPE      = TYPEREAL,                                &
      NDIMS      = 4,                                       &
      LTIMEDEP   = .TRUE.                                   )
!
    TZFIELD%CMNHNAME   = 'FAERO_RA'
    TZFIELD%CLONGNAME  = 'AERODYNAMIC FORCE RA'
    TZFIELD%CUNITS     = 'N'
    TZFIELD%CCOMMENT   = 'RID_AZI_RAD_XYZ instantaneous forces (N) in RA'
    CALL IO_Field_write(TPFILE,TZFIELD,XFAERO_RA_GLB)
!
    TZFIELD%CMNHNAME   = 'FAERO_RG'
    TZFIELD%CLONGNAME  = 'AERODYNAMIC FORCE RG'
    TZFIELD%CUNITS     = 'N'
    TZFIELD%CCOMMENT   = 'RID_AZI_RAD_XYZ instantaneous forces (N) in RG'
    CALL IO_Field_write(TPFILE,TZFIELD,XFAERO_RG_GLB)
!
! * Blade Equivalent Variables (rotor id, radial id, xyz)
    TZFIELD = TFIELDMETADATA(                                  &
      CMNHNAME   = 'AOA_BLEQ',                                 &
      CSTDNAME   = '',                                         &
      CLONGNAME  = 'BLADE_EQ_ANGLE_OF_ATTACK',                 &
      CUNITS     = 'rad',                                      &
      CDIR       = '--',                                       &
      CCOMMENT   = 'RID_RAD blade eq. angle of attack (rad)',  &
      NGRID      = 1,                                          &
      NTYPE      = TYPEREAL,                                   &
      NDIMS      = 2,                                          &
      LTIMEDEP   = .TRUE.                                      )
    CALL IO_Field_write(TPFILE,TZFIELD,XAOA_BLEQ_GLB)

    TZFIELD = TFIELDMETADATA(                                  &
      CMNHNAME   = 'FAERO_BLEQ_RA',                            &
      CSTDNAME   = '',                                         &
      CLONGNAME  = 'BLADE_EQ_AERO_FORCE_RA',                   &
      CUNITS     = 'N',                                        &
      CDIR       = '--',                                       &
      CCOMMENT   = 'RID_RAD_XYZ blade eq. forces (N) in RA',   &
      NGRID      = 1,                                          &
      NTYPE      = TYPEREAL,                                   &
      NDIMS      = 3,                                          &
      LTIMEDEP   = .TRUE.                                      )
    CALL IO_Field_write(TPFILE,TZFIELD,XFAERO_BLEQ_RA_GLB)
!
    IF (MEAN_COUNT /= 0) THEN
!
! * 1D Variables (rotor id)
      TZFIELD = TFIELDMETADATA(                        &
        CMNHNAME   = '1D ADR mean var: (rot)',         & !Temporary name to ease identification
        CSTDNAME   = '',                               &
        CDIR       = '--',                             &
        NGRID      = 1,                                &
        NTYPE      = TYPEREAL,                         &
        NDIMS      = 1,                                &
        LTIMEDEP   = .TRUE.                            ) 
!
      TZFIELD%CMNHNAME   = 'THRUMME'
      TZFIELD%CLONGNAME  = 'MEAN_THRUST_EOL'
      TZFIELD%CUNITS     = 'N'
      TZFIELD%CCOMMENT   = 'RID mean thrust of the wind turbines (N)'
      CALL IO_Field_write(TPFILE,TZFIELD,XTHRU_SUM/MEAN_COUNT)
!
      TZFIELD%CMNHNAME   = 'TORQMME'
      TZFIELD%CLONGNAME  = 'MEAN_TORQUE_EOL'
      TZFIELD%CUNITS     = 'Nm'
      TZFIELD%CCOMMENT   = 'RID mean torque of the wind turbines (Nm)'
      CALL IO_Field_write(TPFILE,TZFIELD,XTORQ_SUM/MEAN_COUNT)
!
      TZFIELD%CMNHNAME   = 'POWMME'
      TZFIELD%CLONGNAME  = 'MEAN_POWER_EOL'
      TZFIELD%CUNITS     = 'W'
      TZFIELD%CCOMMENT   = 'RID mean power of the wind turbines (W)'
      CALL IO_Field_write(TPFILE,TZFIELD,XPOW_SUM/MEAN_COUNT)
!
      TZFIELD = TFIELDMETADATA(                                &
        CMNHNAME   = 'AOAMME',                                 &
        CSTDNAME   = '',                                       &
        CLONGNAME  = 'MEAN_ANGLE_OF_ATTACK',                   &
        CUNITS     = 'rad',                                    &
        CDIR       = '--',                                     &
        CCOMMENT   = 'RID_AZI_RAD mean angle of attack (rad)', &
        NGRID      = 1,                                        &
        NTYPE      = TYPEREAL,                                 &
        NDIMS      = 3,                                        &
        LTIMEDEP   = .TRUE.                                    )
      CALL IO_Field_write(TPFILE,TZFIELD,XAOA_SUM/MEAN_COUNT)
!
      TZFIELD = TFIELDMETADATA(                               &
        CMNHNAME   = 'FAEROMME_RA',                           &
        CSTDNAME   = '',                                      &
        CLONGNAME  = 'MEAN_AERODYNAMIC_FORCE_RA',             &
        CUNITS     = 'N',                                     &
        CDIR       = '--',                                    &
        CCOMMENT   = 'RID_AZI_RAD_XYZ mean forces (N) in RA', &
        NGRID      = 1,                                       &
        NTYPE      = TYPEREAL,                                &
        NDIMS      = 4,                                       &
        LTIMEDEP   = .TRUE.                                   )
      CALL IO_Field_write(TPFILE,TZFIELD,XFAERO_RA_SUM/MEAN_COUNT)
!
! * Blade Equivalent Variables (rotor id, radial id, xyz)
      TZFIELD = TFIELDMETADATA(                                      &
        CMNHNAME   = 'AOAMME_BLEQ',                                  &
        CSTDNAME   = '',                                             &
        CLONGNAME  = 'MEAN_BLADE_EQ_AOA',                            &
        CUNITS     = 'rad',                                          &
        CDIR       = '--',                                           &
        CCOMMENT   = 'RID_RAD mean blade eq. angle of attack (rad)', &
        NGRID      = 1,                                              &
        NTYPE      = TYPEREAL,                                       &
        NDIMS      = 2,                                              &
        LTIMEDEP   = .TRUE.                                          )
      CALL IO_Field_write(TPFILE,TZFIELD,XAOA_BLEQ_SUM/MEAN_COUNT)

      TZFIELD = TFIELDMETADATA(                                       &
        CMNHNAME   = 'FAEROMME_BLEQ_RA',                              &
        CSTDNAME   = '',                                              &
        CLONGNAME  = 'MEAN_BLADE_EQ_AERO_F_RA',                       &
        CUNITS     = 'N',                                             &
        CDIR       = '--',                                            &
        CCOMMENT   = 'RID_RAD_XYZ mean blade eq. forces (N) in RA',   &
        NGRID      = 1,                                               &
        NTYPE      = TYPEREAL,                                        &
        NDIMS      = 3,                                               &
        LTIMEDEP   = .TRUE.                                           )
      CALL IO_Field_write(TPFILE,TZFIELD,XFAERO_BLEQ_RA_SUM/MEAN_COUNT)
!
    END IF
!
!             iv) Actuator Line Model
!
  CASE('ALM') ! Actuator Line Method
!
! * 1D Variables (rotor id)
    TZFIELD = TFIELDMETADATA(                   &
      CMNHNAME   = '1D ALM var: (rot)',         & !Temporary name to ease identification
      CSTDNAME   = '',                          &
      CDIR       = '--',                        &
      NGRID      = 1,                           &
      NTYPE      = TYPEREAL,                    &
      NDIMS      = 1,                           &
      LTIMEDEP   = .TRUE.                       )
!
    TZFIELD%CMNHNAME   = 'THRUT'
    TZFIELD%CLONGNAME  = 'THRUSTT_EOL'
    TZFIELD%CUNITS     = 'N'
    TZFIELD%CCOMMENT   = 'RID instantaneous thrust (N) of wind turbines'
    CALL IO_Field_write(TPFILE,TZFIELD,XTHRUT)
!
    TZFIELD%CMNHNAME   = 'TORQT'
    TZFIELD%CLONGNAME  = 'TORQUET_EOL'
    TZFIELD%CUNITS     = 'Nm'
    TZFIELD%CCOMMENT   = 'RID instantaneous torque (Nm) of wind turbines'
    CALL IO_Field_write(TPFILE,TZFIELD,XTORQT)
!
    TZFIELD%CMNHNAME   = 'POWT'
    TZFIELD%CLONGNAME  = 'POWERT_EOL'
    TZFIELD%CUNITS     = 'W'
    TZFIELD%CCOMMENT   = 'RID instantaneous power (W) of wind turbines'
    CALL IO_Field_write(TPFILE,TZFIELD,XPOWT)
!
! * 3D Variables (rotor id, blade id, radial id)
    TZFIELD = TFIELDMETADATA(                             &
      CMNHNAME   = '3D ALM var: (rot,blade,rad)',         & 
      CSTDNAME   = '',                                    &
      CDIR       = '--',                                  &
      NGRID      = 1,                                     &
      NTYPE      = TYPEREAL,                              &
      NDIMS      = 3,                                     &
      LTIMEDEP   = .TRUE.                                 )
!
    TZFIELD%CMNHNAME   = 'ELT_RAD'
    TZFIELD%CLONGNAME  = 'ELT_RAD'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CCOMMENT   = 'RID_BID_EID radius (m) of wind turbine blade elements'
    CALL IO_Field_write(TPFILE,TZFIELD,XELT_RAD)
!
    TZFIELD%CMNHNAME   = 'AOA'
    TZFIELD%CLONGNAME  = 'ANGLE OF ATTACK'
    TZFIELD%CUNITS     = 'rad'
    TZFIELD%CCOMMENT   = 'RID_BID_EID instantaneous angle of attack (rad)'
    CALL IO_Field_write(TPFILE,TZFIELD,XAOA_GLB)
!
    TZFIELD%CMNHNAME   = 'FLIFT'
    TZFIELD%CLONGNAME  = 'LIFT FORCE'
    TZFIELD%CUNITS     = 'N'
    TZFIELD%CCOMMENT   = 'RID_BID_EID instantaneous lift (N) in relative frame'
    CALL IO_Field_write(TPFILE,TZFIELD,XFLIFT_GLB)
!
    TZFIELD%CMNHNAME   = 'FDRAG'
    TZFIELD%CLONGNAME  = 'DRAG FORCE'
    TZFIELD%CUNITS     = 'N'
    TZFIELD%CCOMMENT   = 'RID_BID_EID instantaneous drag (N) in relative frame'
    CALL IO_Field_write(TPFILE,TZFIELD,XFDRAG_GLB)
!
! * 4D Variables (rotor id, azimuthal id, radial id, xyz)
    TZFIELD = TFIELDMETADATA(                                 &
      CMNHNAME   = '4D ALM var: (rot,blade,rad,xyz)',         & 
      CSTDNAME   = '',                                        &
      CDIR       = '--',                                      &
      NGRID      = 1,                                         &
      NTYPE      = TYPEREAL,                                  &
      NDIMS      = 4,                                         &
      LTIMEDEP   = .TRUE.                                     )
!
    TZFIELD%CMNHNAME   = 'FAERO_RE'
    TZFIELD%CLONGNAME  = 'AERODYNAMIC FORCE RE'
    TZFIELD%CUNITS     = 'N'
    TZFIELD%CCOMMENT   = 'RID_BID_EID_XYZ instantaneous forces (N) in RE'
    CALL IO_Field_write(TPFILE,TZFIELD,XFAERO_RE_GLB)
!
    TZFIELD%CMNHNAME   = 'FAERO_RG'
    TZFIELD%CLONGNAME  = 'AERODYNAMIC FORCE RG'
    TZFIELD%CUNITS     = 'N'
    TZFIELD%CCOMMENT   = 'RID_BID_EID_XYZ instantaneous forces (N) in RG'
    CALL IO_Field_write(TPFILE,TZFIELD,XFAERO_RG_GLB)
!
    IF (MEAN_COUNT /= 0) THEN
!
! * 1D Variables (rotor id)
      TZFIELD = TFIELDMETADATA(                        &
        CMNHNAME   = '1D ALM mean var: (rot)',         & !Temporary name to ease identification
        CSTDNAME   = '',                               &
        CDIR       = '--',                             &
        NGRID      = 1,                                &
        NTYPE      = TYPEREAL,                         &
        NDIMS      = 1,                                &
        LTIMEDEP   = .TRUE.                            )
!
      TZFIELD%CMNHNAME   = 'THRUMME'
      TZFIELD%CLONGNAME  = 'MEAN_THRUST_EOL'
      TZFIELD%CUNITS     = 'N'
      TZFIELD%CCOMMENT   = 'RID mean thrust of the wind turbines (N)'
      CALL IO_Field_write(TPFILE,TZFIELD,XTHRU_SUM/MEAN_COUNT)
!
      TZFIELD%CMNHNAME   = 'TORQMME'
      TZFIELD%CLONGNAME  = 'MEAN_TORQUE_EOL'
      TZFIELD%CUNITS     = 'Nm'
      TZFIELD%CCOMMENT   = 'RID mean torque of the wind turbines (Nm)'
      CALL IO_Field_write(TPFILE,TZFIELD,XTORQ_SUM/MEAN_COUNT)
!
      TZFIELD%CMNHNAME   = 'POWMME'
      TZFIELD%CLONGNAME  = 'MEAN_POWER_EOL'
      TZFIELD%CUNITS     = 'W'
      TZFIELD%CCOMMENT   = 'RID mean power of the wind turbines (W)'
      CALL IO_Field_write(TPFILE,TZFIELD,XPOW_SUM/MEAN_COUNT)
!
      TZFIELD = TFIELDMETADATA(                                &
        CMNHNAME   = 'AOAMME',                                 &
        CSTDNAME   = '',                                       &
        CLONGNAME  = 'MEAN_ANGLE_OF_ATTACK',                   &
        CUNITS     = 'rad',                                    &
        CDIR       = '--',                                     &
        CCOMMENT   = 'RID_BID_EID mean angle of attack (rad)', &
        NGRID      = 1,                                        &
        NTYPE      = TYPEREAL,                                 &
        NDIMS      = 3,                                        &
        LTIMEDEP   = .TRUE.                                    )
      CALL IO_Field_write(TPFILE,TZFIELD,XAOA_SUM/MEAN_COUNT)
!
      TZFIELD = TFIELDMETADATA(                               &
        CMNHNAME   = 'FAEROMME_RE',                           &
        CSTDNAME   = '',                                      &
        CLONGNAME  = 'MEAN_AERODYNAMIC_FORCE_RE',             &
        CUNITS     = 'N',                                     &
        CDIR       = '--',                                    &
        CCOMMENT   = 'RID_BID_EID_XYZ mean forces (N) in RE', &
        NGRID      = 1,                                       &
        NTYPE      = TYPEREAL,                                &
        NDIMS      = 4,                                       &
        LTIMEDEP   = .TRUE.                                   )
      CALL IO_Field_write(TPFILE,TZFIELD,XFAERO_RE_SUM/MEAN_COUNT)
!
    END IF
!
  END SELECT
END IF 
!
DEALLOCATE(ZWORK2D,ZWORK3D)
!
!-------------------------------------------------------------------------------!
!
END SUBROUTINE WRITE_LFIFM_n  
