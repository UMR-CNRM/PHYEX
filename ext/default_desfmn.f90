!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################
      MODULE MODI_DEFAULT_DESFM_n
!     ###########################
!
INTERFACE
!
SUBROUTINE DEFAULT_DESFM_n(KMI)
INTEGER,         INTENT(IN)  :: KMI       ! Model index
END SUBROUTINE DEFAULT_DESFM_n
!
END INTERFACE
!
END MODULE MODI_DEFAULT_DESFM_n
!
!
!
!     ###############################
      SUBROUTINE DEFAULT_DESFM_n(KMI)
!     ###############################
!
!!****  *DEFAULT_DESFM_n * - set default values for descriptive variables of
!!                         model KMI
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set default values for the variables
!     in descriptor files by filling the corresponding variables which
!     are stored in modules.
!
!
!!**  METHOD
!!    ------
!!      Each variable in modules, which can be initialized by  reading its
!!    value in the descriptor file is set to a default value.
!!     When this routine is used during INIT, the modules of the first model
!!   are used to temporarily store  the variables associated with a nested
!!   model.
!!     When this routine is used during  SPAWNING, the modules of a second
!!   model must be initialized.
!!     Default values for variables common to all models are set only
!!   at the first call of DEFAULT_DESFM_n (i.e. when KMI=1)
!!
!!
!!    EXTERNAL
!!    --------
!!     NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : JPHEXT,JPVEXT
!!
!!      Module MODD_CONF       : CCONF,L2D,L1D,LFLAT,NMODEL,NVERB
!!
!!      Module MODD_DYN        : XSEGLEN,XASSELIN,LCORIO,LNUMDIFF
!!                               XALKTOP,XALZBOT
!!
!!      Module MODD_BAKOUT
!!
!!      Module MODD_NESTING    : NDAD(m),NDTRATIO(m),XWAY(m)
!!
!!      Module MODD_CONF_n    : LUSERV,LUSERC,LUSERR,LUSERI,LUSERS
!!                              LUSERG,LUSERH,CSEG,CEXP
!!
!!      Module MODD_LUNIT_n   : CINIFILE,CCPLFILE
!!
!!
!!      Module MODD_DYN_n     : XTSTEP,CPRESOPT,NITR,XRELAX,LHO_RELAX
!!         LVE_RELAX,XRIMKMAX,NRIMX,NRIMY
!!
!!      Module MODD_ADV_n : CUVW_ADV_SCHEME,CMET_ADV_SCHEME,CSV_ADV_SCHEME,NLITER
!!
!!      Module MODD_PARAM_n : CTURB,CRAD,CDCONV,CSCONV
!!
!!      Module MODD_LBC_n : CLBCX, CLBCY,NLBLX,NLBLY,XCPHASE,XCPHASE_PBL,XPOND
!!
!!      Module MODD_TURB_n : XIMPL,CTURBLEN,CTURBDIM,LTURB_FLX,LTURB_DIAG,LSUBG_COND
!!                           LTGT_FLX
!!
!!
!!      Module MODD_PARAM_RAD_n:
!!          XDTRAD,XDTRAD_CLONLY,LCLEAR_SKY,NRAD_COLNBR, NRAD_DIAG
!!
!!      Module MODD_BUDGET : CBUTYPE,NBUMOD,XBULEN,NBUKL, NBUKH,LBU_KCP,XBUWRI
!!         NBUIL, NBUIH,NBUJL, NBUJH,LBU_ICP,LBU_JCP,NBUMASK
!!
!!      Module MODD_BLANK_n:
!!
!!          XDUMMYi, NDUMMYi, LDUMMYi, CDUMMYi
!!
!!      Module MODD_FRC :
!!
!!          LGEOST_UV_FRC,LGEOST_TH_FRC,LTEND_THRV_FRC
!!          LVERT_MOTION_FRC,LRELAX_THRV_FRC,LRELAX_UV_FRC,LRELAX_UVMEAN_FRC,
!!          XRELAX_TIME_FRC
!!          XRELAX_HEIGHT_FRC,CRELAX_HEIGHT_TYPE,LTRANS,XUTRANS,XVTRANS,
!!          LPGROUND_FRC
!!
!!      Module MODD_PARAM_ICE :
!!
!!          LWARM,CPRISTINE_ICE
!!
!!      Module MODD_PARAM_KAFR_n :
!!
!!          XDTCONV,LREFRESH_ALL,LDOWN,NICE,LCHTRANS
!!
!!      Module MODD_PARAM_MFSHALL_n :
!!
!!         CMF_UPDRAFT,LMIXUV,CMF_CLOUD,XIMPL_MF,LMF_FLX              
!!
!!
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (routine DEFAULT_DESFM_n)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      02/06/94
!!      Modifications 17/10/94  (Stein)  For LCORIO
!!      Modifications 06/12/94  (Stein)  remove LBOUSS+add LABSLAYER, LNUMDIFF
!!                                       ,LSTEADYLS
!!      Modifications 06/12/94  (Stein)  remove LABSLAYER, add LHO_RELAX,
!!                                       LVE_RELAX, NRIMX, NRIMY, XRIMKMAX
!!      Modifications 09/01/95  (Lafore) add LSTEADY_DMASS
!!      Modifications 09/01/95  (Stein)  add the turbulence scheme namelist
!!      Modifications 09/01/95  (Stein)  add the 1D switch
!!      Modifications 10/03/95  (Mallet) add the coupling files
!!                    29/06/95  ( Stein, Nicolau, Hereil) add the budgets
!!      Modifications 25/09/95  ( Stein )add the LES tools
!!      Modifications 25/10/95  ( Stein )add the radiations
!!      Modifications 23/10/95  (Vila, lafore) new scalar advection scheme
!!      Modifications 24/02/96  (Stein)  change the default value for CCPLFILE
!!      Modifications 12/02/96  (Lafore) transformation to DEFAULT_DESFM_n for
!!                                       spawning
!!      Modifications 25/04/96  (Suhre)  add the blank module
!!      Modifications 29/07/96  (Pinty&Suhre) add module MODD_FRC
!!      Modifications 11/04/96  (Pinty)  add the rain-ice scheme and modify
!!                                       the split arrays in MODD_PARAM_RAD_n
!!      Modifications 11/01/97  (Pinty)  add the deep convection scheme
!!      Modifications 24/11/96  (Masson)  add LREFRESH_ALL in deep convection
!!      Modifications 12/02/96  (Lafore) transformation to DEFAULT_DESFM_n for spawning
!!      Modifications 22/07/96  (Lafore) gridnesting implementation
!!      Modifications 29/07/96  (Lafore) add the module MODD_FMOUT (renamed MODD_BAKOUT)
!!      Modifications 23/06/97  (Stein)  add the equation system name
!!      Modifications 10/07/97  (Masson) add MODD_PARAM_GROUNDn : CROUGH
!!      Modifications 28/07/97  (Masson) remove LREFRESH_ALL and LSTEADY_DMASS
!!      Modifications 08/10/97  (Stein)  switch (_n=1) to initialize the
!!                                       parameters common to all models
!!      Modifications 24/01/98 (Bechtold) add LREFRESH_ALL, LCHTRANS,
!!                                         LTEND_THRV_FR and LSST_FRC
!!      Modifications 18/07/99  (Stein)  add LRAD_DIAG
!!      Modification  15/03/99 (Masson)  use of XUNDEF
!!      Modification  11/12/00 (Tomasini) Add CSEA_FLUX to MODD_PARAMn
!!      Modification  22/01/01 (Gazen) delete NSV and add LHORELAX_SVC2R2
!!                                     LHORELAX_SVCHEM,LHORELAX_SVLG
!!      Modification 15/03/02 (Solmon) radiation scheme: remove NSPOT and add
!!                                   default for aerosol and cloud rad. prop. control
!!      Modification 22/05/02 (Jabouille) put chimical default here
!!      Modification 01/2004  (Masson) removes surface (externalization)
!!                      09/04 (M. Tomasini) New namelist to modify the
!!                                             Cloud mixing length
!!                   07/05 (P.Tulet) New namelists for dust and aerosol
!!      Modification 01/2007  (Malardel, Pergaud) Add MODD_PARAM_MFSHALL_n
!!      Modification 10/2009  (Aumond) Add user multimasks for LES
!!      Modification 10/2009  (Aumond) Add MEAN_FIELD              
!!      Modification 12/04/07 (Leriche) add LUSECHAQ for aqueous chemistry
!!      Modification 30/05/07 (Leriche) add LCH_PH and XCH_PHINIT for pH
!!      Modification 25/04/08 (Leriche) add XRTMIN_AQ LWC threshold for aq. chemistry
!!                   16/07/10           add LHORELAX_SVIC
!!                   16/09/10           add LUSECHIC
!!                   13/01/11           add LCH_RET_ICE
!!                   01/07/11 (F.Couvreux) Add CONDSAMP
!!                   01/07/11 (B.Aouizerats) Add CAOP    
!!                   07/2013  (C.Lac) add WENO, LCHECK              
!!                   07/2013  (Bosseur & Filippi) adds Forefire
!!                   08/2015  (Redelsperger & Pianezze) add XPOND coefficient for LBC
!!      Modification 24/03/16 (Leriche) remove LCH_SURFACE_FLUX 
!!                                      put NCH_VEC_LENGTH = 50 instead of 1000
!!
!!                   04/2016 (C.LAC) negative contribution to the budget split between advection, turbulence and microphysics for KHKO/C2R2
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!      Modification 24/03/16 (Leriche) remove LCH_SURFACE_FLUX 
!!                                      put NCH_VEC_LENGTH = 50 instead of 1000
!!                   10/2016 (C.Lac) VSIGQSAT change from 0 to 0.02 for coherence with AROME
!!                    10/2016 (C.Lac) Add droplet deposition
!!                   10/2016  (R.Honnert and S.Riette) : Improvement of EDKF and adaptation to the grey zone
!!                   10/2016  (F Brosse) add prod/loss terms computation for chemistry
!!                   07/2017  (V. Masson) adds time step for output files writing.
!!                   09/2017 Q.Rodier add LTEND_UV_FRC
!!                   02/2018 Q.Libois ECRAD
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                   01/2018  (S. Riette) new budgets and variables for ICE3/ICE4
!!                   01/2018 (J.Colin) add VISC and DRAG
!!                   07/2017 (V. Vionnet) add blowing snow variables
!!                   01/2019 (R. Honnert) add reduction of the mass-flux surface closure with the resolution
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!!                   05/2019 F.Brient add tracer emission from the top of the boundary-layer
!!                   11/2019 C.Lac correction in the drag formula and application to building in addition to tree
!  P. Wautelet 17/04/2020: move budgets switch values into modd_budget
!  P. Wautelet 30/06/2020: add NNETURSV, NNEADVSV and NNECONSV variables
!  F. Auguste, T. Nagel 02/2021: add IBM defaults parameters
!  T. Nagel       02/2021: add turbulence recycling defaults parameters
!  P-A Joulin  21/05/2021: add Wind turbines
!  S. Riette   21/05/2021: add options to PDF subgrid scheme
!  D. Ricard      05/2021: add the contribution of Leonard terms in the turbulence scheme
!  JL Redelsperger 06/2021: add parameters allowing to active idealized oceanic convection
!  B. Vie         06/2021: add prognostic supersaturation for LIMA
!  Q. Rodier      06/2021: modify default value to LGZ=F (grey-zone corr.), LSEDI and OSEDC=T (LIMA sedimentation)
!  F. Couvreux    06/2021: add LRELAX_UVMEAN_FRC
!  Q. Rodier      07/2021: modify XPOND=1
!  A. Costes      12/2021: Blaze fire model
!  C. Barthe      03/2022: add CIBU and RDSF options in LIMA
!  Delbeke/Vie    03/2022: KHKO option in LIMA
!  P. Wautelet 27/04/2022: add namelist for profilers
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_PARAMETERS
USE MODD_CONF             !        For INIT only DEFAULT_DESFM1
USE MODD_CONFZ
USE MODD_DYN
USE MODD_NESTING
USE MODD_BAKOUT
USE MODD_SERIES
USE MODD_CONF_n           ! modules used to set the default values is only
USE MODD_LUNIT_n          ! the one corresponding to model 1. These memory
USE MODD_DIM_n            ! addresses will then be filled by the values read in
USE MODD_DYN_n            ! the DESFM corresponding to model n which may have
USE MODD_ADV_n            ! missing values. This is why we affect default values.
USE MODD_PARAM_n          !      For SPAWNING DEFAULT_DESFM2 is also used
USE MODD_LBC_n
USE MODD_OUT_n
USE MODD_TURB_n
USE MODD_BUDGET
USE MODD_LES
USE MODD_PARAM_RAD_n
#ifdef MNH_ECRAD
USE MODD_PARAM_ECRAD_n
#if ( VER_ECRAD == 140 ) 
USE MODD_RADIATIONS_n , ONLY : NSWB_MNH, NLWB_MNH
#endif
#endif
USE MODD_BLANK_n
USE MODD_FRC
USE MODD_PARAM_ICE
USE MODD_PARAM_C2R2
USE MODD_TURB_CLOUD
USE MODD_PARAM_KAFR_n
USE MODD_PARAM_MFSHALL_n
USE MODD_CH_MNHC_n
USE MODD_SERIES_n
USE MODD_NUDGING_n
USE MODD_CH_AEROSOL
USE MODD_DUST
USE MODD_SALT
USE MODD_PASPOL
USE MODD_CONDSAMP
USE MODD_MEAN_FIELD
USE MODD_DRAGTREE_n
USE MODD_DRAGBLDG_n
USE MODD_EOL_MAIN
USE MODD_EOL_ADNR
USE MODD_EOL_ALM
USE MODD_EOL_SHARED_IO
USE MODD_ALLPROFILER_n
USE MODD_ALLSTATION_n
!
!
USE MODD_PARAM_LIMA, ONLY : LCOLD, LNUCL, LSEDI, LHHONI, LSNOW, LHAIL, LMEYERS,         &
                            NMOM_I, NMOM_S, NMOM_G, NMOM_H,                             &
                            NMOD_IFN, XIFN_CONC, LIFN_HOM, CIFN_SPECIES,                &
                            CINT_MIXING, NMOD_IMM, NIND_SPECIE, LMURAKAMI,              &
                            YSNOW_T=>LSNOW_T, CPRISTINE_ICE_LIMA, CHEVRIMED_ICE_LIMA,   &
                            XFACTNUC_DEP, XFACTNUC_CON,                                 &
                            OWARM=>LWARM, LACTI, ORAIN=>LRAIN, OSEDC=>LSEDC,            &
                            OACTIT=>LACTIT, LBOUND, LSPRO, LADJ, LKHKO, NMOM_C, NMOM_R, &
                            NMOD_CCN, XCCN_CONC, LKESSLERAC,                            &
                            LCCN_HOM, CCCN_MODES,                                       &
                            YALPHAR=>XALPHAR, YNUR=>XNUR,                               &
                            YALPHAC=>XALPHAC, YNUC=>XNUC, CINI_CCN=>HINI_CCN,           &
                            CTYPE_CCN=>HTYPE_CCN, YFSOLUB_CCN=>XFSOLUB_CCN,             &
                            YACTEMP_CCN=>XACTEMP_CCN, YAERDIFF=>XAERDIFF,               &
                            YAERHEIGHT=>XAERHEIGHT,                                     &
                            LSCAV, LAERO_MASS, NPHILLIPS,                               &
                            LCIBU, XNDEBRIS_CIBU, LRDSF,                                &
                            ODEPOC=>LDEPOC, OVDEPOC=>XVDEPOC, OACTTKE=>LACTTKE,         &
                            LPTSPLIT, L_LFEEDBACKT=>LFEEDBACKT, L_NMAXITER=>NMAXITER,   &
                            L_XMRSTEP=>XMRSTEP, L_XTSTEP_TS=>XTSTEP_TS
!
USE MODD_LATZ_EDFLX
USE MODD_2D_FRC
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
USE MODD_DRAG_n
USE MODD_VISCOSITY
USE MODD_RECYCL_PARAM_n
USE MODD_IBM_PARAM_n
USE MODD_IBM_LSF
#ifdef MNH_FOREFIRE
USE MODD_FOREFIRE
#endif
USE MODD_FIRE
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER,         INTENT(IN)  :: KMI       ! Model index
!
!*       0.2   declaration of local variables
!
INTEGER             :: JM      ! loop index
!
!-------------------------------------------------------------------------------
!
!*      1.    SET DEFAULT VALUES FOR MODD_LUNIT_n :
!             ----------------------------------
!
! CINIFILE='INIFILE'
CINIFILEPGD='' !Necessary to keep this line to prevent problems with spawning
CCPLFILE(:)='    '
!
!-------------------------------------------------------------------------------
!
!*      2.    SET DEFAULT VALUES FOR MODD_CONF AND MODD_CONF_n :
!             ------------------------------------------------
!
IF (KMI == 1) THEN
  CCONF      ='START'
  LTHINSHELL = .FALSE.
  L2D        = .FALSE.
  L1D        = .FALSE.
  LFLAT      = .FALSE.
  NMODEL     = 1
  CEQNSYS    = 'DUR'
  NVERB      = 5
  CEXP       = 'EXP01'
  CSEG       = 'SEG01'
  LFORCING   = .FALSE.
  L2D_ADV_FRC= .FALSE.
  L2D_REL_FRC= .FALSE. 
  XRELAX_HEIGHT_BOT = 0.
  XRELAX_HEIGHT_TOP = 30000.
  XRELAX_TIME = 864000.
  LPACK      = .TRUE.
  NHALO      = 1
#ifdef MNH_SX5
  CSPLIT     ='YSPLITTING'   ! NEC vectoriel architecture , low number of PROC
#else
  CSPLIT     ='BSPLITTING'   ! Scalaire architecture , high number of PROC
#endif
  NZ_PROC    = 0               !JUAN Z_SPLITTING :: number of proc in Z splitting
  NZ_SPLITTING = 10            !JUAN Z_SPLITTING :: for debug NZ=1=flat_inv;  NZ=10=flat_invz; NZ=1+2 the two
  LLG        = .FALSE.
  LINIT_LG   = .FALSE.
  CINIT_LG   = 'FMOUT'
  LNOMIXLG   = .FALSE.
  LCHECK     = .FALSE.
END IF
!
CCLOUD    = 'NONE'
LUSERV    = .TRUE.
LUSERC    = .FALSE.
LUSERR    = .FALSE.
LUSERI    = .FALSE.
LUSERS    = .FALSE.
LUSERG    = .FALSE.
LUSERH    = .FALSE.
LOCEAN    = .FALSE.
!NSV      = 0
!NSV_USER = 0
LUSECI    = .FALSE.
!
!-------------------------------------------------------------------------------
!
!*      3.    SET DEFAULT VALUES FOR MODD_DYN AND MODD_DYN_n :
!             -----------------------------------------------
!
IF (KMI == 1) THEN
  XSEGLEN   = 43200.
  XASSELIN  = 0.2
  XASSELIN_SV = 0.02
  LCORIO    = .TRUE.
  LNUMDIFU  = .TRUE.
  LNUMDIFTH = .FALSE.
  LNUMDIFSV = .FALSE.
  XALZBOT   = 4000.
  XALKTOP   = 0.01
  XALKGRD   = 0.01
  XALZBAS   = 0.01
END IF
!
XTSTEP    = 60.
CPRESOPT  = 'CRESI'
NITR      = 4
LITRADJ   = .TRUE.
LRES      = .FALSE.
XRES      = 1.E-07
XRELAX    = 1.
LVE_RELAX = .FALSE.
LVE_RELAX_GRD = .FALSE.
XRIMKMAX  = 0.01 / XTSTEP
XT4DIFU   = 1800.          
XT4DIFTH  = 1800.           
XT4DIFSV  = 1800.           
!
IF (KMI == 1) THEN   ! for model 1 we have a Large scale information
  NRIMX = JPRIMMAX   ! for U,V,W,TH,Rv used for the hor. relaxation
  NRIMY = JPRIMMAX
ELSE
  NRIMX = 0          ! for inner models we use only surfacic fields to
  NRIMY = 0          ! give the lbc and no hor. relaxation is used
END IF
!
LHORELAX_UVWTH = .FALSE.
LHORELAX_RV = .FALSE.
LHORELAX_RC = .FALSE. ! for all these fields, no large scale is usally available
LHORELAX_RR = .FALSE. ! for model 1 and for inner models, we only use surfacic
LHORELAX_RS = .FALSE. ! fiels ( no hor. relax. )
LHORELAX_RI = .FALSE.
LHORELAX_RG = .FALSE.
LHORELAX_RH = .FALSE.
LHORELAX_TKE = .FALSE.
LHORELAX_SV(:) = .FALSE.
LHORELAX_SVC2R2 = .FALSE.
LHORELAX_SVC1R3 = .FALSE.
LHORELAX_SVELEC = .FALSE.
LHORELAX_SVLG   = .FALSE.
LHORELAX_SVCHEM = .FALSE.
LHORELAX_SVCHIC = .FALSE.
LHORELAX_SVDST  = .FALSE.
LHORELAX_SVSLT  = .FALSE.
LHORELAX_SVPP   = .FALSE.
LHORELAX_SVCS   = .FALSE.
LHORELAX_SVAER  = .FALSE.
!
LHORELAX_SVLIMA = .FALSE.
!
#ifdef MNH_FOREFIRE
LHORELAX_SVFF   = .FALSE.
#endif
LHORELAX_SVSNW  = .FALSE.
LHORELAX_SVFIRE = .FALSE.
!
!
!-------------------------------------------------------------------------------
!
!*      4.    SET DEFAULT VALUES FOR MODD_NESTING :
!             -----------------------------------
!
IF (KMI == 1) THEN
  NDAD(1)=1
  DO JM=2,JPMODELMAX
    NDAD(JM)  = JM - 1
  END DO
  NDTRATIO(:) = 1
  XWAY(:)     = 2.      ! two-way interactive gridnesting
  XWAY(1)     = 0.      ! except for model 1
END IF
!
!-------------------------------------------------------------------------------
!
!*      5.    SET DEFAULT VALUES FOR MODD_ADV_n :
!             ----------------------------------
!
CUVW_ADV_SCHEME =  'CEN4TH'
CMET_ADV_SCHEME =  'PPM_01'
CSV_ADV_SCHEME  =  'PPM_01'
CTEMP_SCHEME    =  'RKC4'        
NWENO_ORDER = 3
NSPLIT      = 1
LSPLIT_CFL  = .TRUE.
LSPLIT_WENO = .TRUE.
XSPLIT_CFL  = 0.8
LCFL_WRIT   = .FALSE.
!
!-------------------------------------------------------------------------------
!
!*      6.    SET DEFAULT VALUES FOR MODD_PARAM_n :
!             -----------------------------------
!
CTURB   = 'NONE'
CRAD    = 'NONE'
CDCONV  = 'NONE'
CSCONV  = 'NONE'
CELEC   = 'NONE'
CACTCCN = 'NONE'
!
!-------------------------------------------------------------------------------
!
!*      7.    SET DEFAULT VALUES FOR MODD_LBC_n :
!             ---------------------------------
!
CLBCX(1) ='CYCL'
CLBCX(2) ='CYCL'
CLBCY(1) ='CYCL'
CLBCY(2) ='CYCL'
NLBLX(:) = 1
NLBLY(:) = 1
XCPHASE = 20.
XCPHASE_PBL = 0.
XCARPKMAX = XUNDEF
XPOND = 1.0
!
!-------------------------------------------------------------------------------
!
!*      8.    SET DEFAULT VALUES FOR MODD_NUDGING_n :
!             ---------------------------------
!
LNUDGING = .FALSE.
XTNUDGING = 21600.
!
!-------------------------------------------------------------------------------
!
!*      9.    SET DEFAULT VALUES FOR MODD_BAKOUT and MODD_OUT_n :
!             ------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
!
!*      10.   SET DEFAULT VALUES FOR MODD_TURB_n :
!             ----------------------------------
!
XIMPL     = 1.
XKEMIN    = 0.01
XCEDIS    = 0.84
XCADAP    = 0.5
CTURBLEN  = 'BL89'
CTURBDIM  = '1DIM'
LTURB_FLX =.FALSE.
LTURB_DIAG=.FALSE.
LSUBG_COND=.FALSE.
CSUBG_AUCV='NONE' 
CSUBG_AUCV_RI='NONE'
LSIGMAS   =.TRUE.
LSIG_CONV =.FALSE.
LRMC01    =.FALSE.
CTOM      ='NONE'
VSIGQSAT  = 0.02
CCONDENS='CB02'
CLAMBDA3='CB'
CSUBG_MF_PDF='TRIANGLE'
LLEONARD =.FALSE.
XCOEFHGRADTHL = 1.0
XCOEFHGRADRM = 1.0
XALTHGRAD = 2000.0
XCLDTHOLD = -1.0

!-------------------------------------------------------------------------------
!
!*      10b.   SET DEFAULT VALUES FOR MODD_DRAGTREE :
!             ----------------------------------
!
LDRAGTREE = .FALSE.
LDEPOTREE = .FALSE.
XVDEPOTREE = 0.02 ! 2 cm/s 
!------------------------------------------------------------------------------
!
!*      10c.   SET DEFAULT VALUES FOR MODD_DRAGB
!             ----------------------------------
!
LDRAGBLDG = .FALSE.
!
!*      10d.   SET DEFAULT VALUES FOR MODD_EOL* :
!             ----------------------------------
!
!       10d.i) MODD_EOL_MAIN
!
LMAIN_EOL  = .FALSE.
CMETH_EOL  = 'ADNR'
CSMEAR     = '3LIN'
NMODEL_EOL = 1
!
!       10d.ii) MODD_EOL_SHARED_IO
!
CFARM_CSVDATA     = 'data_farm.csv'
CTURBINE_CSVDATA  = 'data_turbine.csv'
CBLADE_CSVDATA    = 'data_blade.csv' 
CAIRFOIL_CSVDATA  = 'data_airfoil.csv'
!
CINTERP           = 'CLS'
!
!       10d.iii) MODD_EOL_ALM
!
NNB_BLAELT        =  42
LTIMESPLIT        = .FALSE.
LTIPLOSSG         = .TRUE.
LTECOUTPTS        = .FALSE.
!
!------------------------------------------------------------------------------
!*      10.e   SET DEFAULT VALUES FOR MODD_ALLPROFILER_n :
!             ----------------------------------
!
NNUMB_PROF    = 0
XSTEP_PROF    = 60.0
XX_PROF(:)    = XUNDEF
XY_PROF(:)    = XUNDEF
XZ_PROF(:)    = XUNDEF
XLAT_PROF(:)  = XUNDEF
XLON_PROF(:)  = XUNDEF
CNAME_PROF(:) = ''
CFILE_PROF    = 'NO_INPUT_CSV'
! LDIAG_SURFRAD = .TRUE.
!------------------------------------------------------------------------------
!*      10.f   SET DEFAULT VALUES FOR MODD_ALLSTATION_n :
!             ----------------------------------
!
NNUMB_STAT    = 0
XSTEP_STAT    = 60.0
XX_STAT(:)    = XUNDEF
XY_STAT(:)    = XUNDEF
XZ_STAT(:)    = XUNDEF
XLAT_STAT(:)  = XUNDEF
XLON_STAT(:)  = XUNDEF
CNAME_STAT(:) = ''
CFILE_STAT    = 'NO_INPUT_CSV'
LDIAG_SURFRAD = .TRUE.
!
!-------------------------------------------------------------------------------
!
!*      11.   SET DEFAULT VALUES FOR MODD_BUDGET :
!             ------------------------------------
!
!       11.1 General  budget variables
!
IF (KMI == 1) THEN
  CBUTYPE = 'NONE'
  NBUMOD = 1
  XBULEN = XSEGLEN
  XBUWRI = XSEGLEN
  NBUKL  = 1
  NBUKH  = 0        
  LBU_KCP = .TRUE.
!
!       11.2 Variables for the cartesian box
!
  NBUIL = 1
  NBUIH = 0      
  NBUJL = 1
  NBUJH = 0          
  LBU_ICP = .TRUE.
  LBU_JCP = .TRUE.
!
!       11.3 Variables for the mask
!
  NBUMASK = 1
END IF
!
!-------------------------------------------------------------------------------
!
!*      12.    SET DEFAULT VALUES FOR MODD_LES :
!             ---------------------------------
!
IF (KMI == 1) THEN
  LLES_MEAN               = .FALSE.
  LLES_RESOLVED           = .FALSE.
  LLES_SUBGRID            = .FALSE.
  LLES_UPDRAFT            = .FALSE.
  LLES_DOWNDRAFT          = .FALSE.
  LLES_SPECTRA            = .FALSE.
!
  NLES_LEVELS             = NUNDEF
  XLES_ALTITUDES          = XUNDEF
  NSPECTRA_LEVELS         = NUNDEF
  XSPECTRA_ALTITUDES      = XUNDEF
  NLES_TEMP_SERIE_I       = NUNDEF
  NLES_TEMP_SERIE_J       = NUNDEF
  NLES_TEMP_SERIE_Z       = NUNDEF
  CLES_NORM_TYPE          = 'NONE'
  CBL_HEIGHT_DEF          = 'KE'
  XLES_TEMP_SAMPLING      = XUNDEF
  XLES_TEMP_MEAN_START    = XUNDEF
  XLES_TEMP_MEAN_END      = XUNDEF
  XLES_TEMP_MEAN_STEP     = 3600.
  LLES_CART_MASK          = .FALSE.
  NLES_IINF               = NUNDEF
  NLES_ISUP               = NUNDEF
  NLES_JINF               = NUNDEF
  NLES_JSUP               = NUNDEF
  LLES_NEB_MASK           = .FALSE.
  LLES_CORE_MASK          = .FALSE.
  LLES_MY_MASK            = .FALSE.
  NLES_MASKS_USER         = NUNDEF
  LLES_CS_MASK            = .FALSE.

  LLES_PDF               = .FALSE.
  NPDF                   = 1
  XTH_PDF_MIN            = 270.
  XTH_PDF_MAX            = 350.
  XW_PDF_MIN             = -10.
  XW_PDF_MAX             = 10.
  XTHV_PDF_MIN           = 270.
  XTHV_PDF_MAX           = 350.
  XRV_PDF_MIN            = 0.
  XRV_PDF_MAX            = 20.
  XRC_PDF_MIN            = 0.
  XRC_PDF_MAX            = 1.
  XRR_PDF_MIN            = 0.
  XRR_PDF_MAX            = 1.
  XRI_PDF_MIN            = 0.
  XRI_PDF_MAX            = 1.
  XRS_PDF_MIN            = 0.
  XRS_PDF_MAX            = 1.
  XRG_PDF_MIN            = 0.
  XRG_PDF_MAX            = 1.
  XRT_PDF_MIN            = 0.
  XRT_PDF_MAX            = 20.
  XTHL_PDF_MIN           = 270.
  XTHL_PDF_MAX           = 350.
END IF
!
!-------------------------------------------------------------------------------
!
!*      13.   SET DEFAULT VALUES FOR MODD_PARAM_RAD_n :
!             ---------------------------------------
!
XDTRAD        = XTSTEP
XDTRAD_CLONLY = XTSTEP
LCLEAR_SKY    =.FALSE.
NRAD_COLNBR   = 1000
NRAD_DIAG     = 0
CLW ='RRTM'
CAER='SURF'
CAOP='CLIM'
CEFRADL='MART'
CEFRADI='LIOU'
COPWSW = 'FOUQ'
COPISW = 'EBCU'
COPWLW = 'SMSH'
COPILW = 'EBCU'
XFUDG = 1.
LAERO_FT=.FALSE.
LFIX_DAT=.FALSE.
!
#ifdef MNH_ECRAD
!*      13bis.   SET DEFAULT VALUES FOR MODD_PARAM_ECRAD_n :
!             ---------------------------------------
!
#if ( VER_ECRAD == 101 )
NSWSOLVER = 0           ! 0: 'McICA 1: 'SPARTACUS' 2: 'SPARTACUS' + 3D effect                            
NLWSOLVER = 0           ! 0: 'McICA 1: 'SPARTACUS' 2: 'SPARTACUS' + 3D effect 
#endif
#if ( VER_ECRAD == 140 )
LSPEC_ALB = .FALSE.
LSPEC_EMISS = .FALSE.


!ALLOCATE(USER_ALB_DIFF(NSWB_MNH))
!ALLOCATE(USER_ALB_DIR(NSWB_MNH))
!ALLOCATE(USER_EMISS(NLWB_MNH))
!PRINT*,USER_ALB_DIFF
!USER_ALB_DIFF = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
!USER_ALB_DIR = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
!USER_EMISS = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
SURF_TYPE="SNOW"

NLWSOLVER = 1           ! 0: 'McICA 1: 'SPARTACUS' 2: 'SPARTACUS' + 3D effect 
NSWSOLVER = 1          ! 0: 'McICA 1: 'SPARTACUS' 2: 'SPARTACUS' + 3D effect                            
#endif
! LEFF3D         = .TRUE.
! LSIDEM         = .TRUE.
NREG           = 3            ! Number of cloudy regions (3=TripleClouds)
! LLWCSCA        = .TRUE.       ! LW cloud scattering
! LLWASCA        = .TRUE.       ! LW aerosols scattering
NLWSCATTERING = 2
NAERMACC = 0
! CGAS           = 'RRTMG-IFS'  ! Gas optics model
NOVLP          = 1     ! overlap assumption ; 0= 'Max-Ran' ; 1= 'Exp-Ran'; 2 = 'Exp-Exp'
NLIQOPT        = 3            ! 1: 'Monochromatic', 2: 'HuStamnes', 3: 'SOCRATES', 4: 'Slingo' 
NICEOPT        = 3            ! 1: 'Monochromatic', 2: 'Fu-PSRAD', 3: 'Fu-IFS', 4: 'Baran', 5: 'Baran2016', 6: 'Baran2017' 
! LSW_ML_E       = .FALSE.
! LLW_ML_E       = .FALSE.
! LPSRAD         = .FALSE.
! 
NRADLP = 1                    ! 0: ERA-15, 1: Zhang and Rossow, 2: Martin (1994) et Woods (2000)
NRADIP = 1                    ! 0: 40 mum, 1: Liou and Ou (1994), 2: Liou and Ou (1994) improved, 3: Sun and Rikus (1999)
XCLOUD_FRAC_STD = 1.0_JPRB ! change to 0.75 for more realistic distribution
#endif
!-------------------------------------------------------------------------------
!
!*      14.   SET DEFAULT VALUES FOR MODD_BLANK_n :
!             -----------------------------------
!
XDUMMY1       = 0.
XDUMMY2       = 0.
XDUMMY3       = 0.
XDUMMY4       = 0.
XDUMMY5       = 0.
XDUMMY6       = 0.
XDUMMY7       = 0.
XDUMMY8       = 0.
!
NDUMMY1       = 0
NDUMMY2       = 0
NDUMMY3       = 0
NDUMMY4       = 0
NDUMMY5       = 0
NDUMMY6       = 0
NDUMMY7       = 0
NDUMMY8       = 0
!
LDUMMY1       = .TRUE.
LDUMMY2       = .TRUE.
LDUMMY3       = .TRUE.
LDUMMY4       = .TRUE.
LDUMMY5       = .TRUE.
LDUMMY6       = .TRUE.
LDUMMY7       = .TRUE.
LDUMMY8       = .TRUE.
!
CDUMMY1       = ' '
CDUMMY2       = ' '
CDUMMY3       = ' '
CDUMMY4       = ' '
CDUMMY5       = ' '
CDUMMY6       = ' '
CDUMMY7       = ' '
CDUMMY8       = ' '
!
!------------------------------------------------------------------------------
!
!*      15.   SET DEFAULT VALUES FOR MODD_FRC :
!             ---------------------------------
!
IF (KMI == 1) THEN
  LGEOST_UV_FRC      = .FALSE.
  LGEOST_TH_FRC      = .FALSE.
  LTEND_THRV_FRC      = .FALSE.
  LTEND_UV_FRC       = .FALSE. 
  LVERT_MOTION_FRC   = .FALSE.
  LRELAX_THRV_FRC    = .FALSE.
  LRELAX_UV_FRC      = .FALSE.
  LRELAX_UVMEAN_FRC  = .FALSE.
  XRELAX_TIME_FRC    = 10800.
  XRELAX_HEIGHT_FRC  = 0.
  CRELAX_HEIGHT_TYPE = "FIXE"
  LTRANS             = .FALSE.
  XUTRANS            = 0.0
  XVTRANS            = 0.0
  LPGROUND_FRC       = .FALSE.
  LDEEPOC   = .FALSE.
  XCENTX_OC = 16000.
  XCENTY_OC = 16000.
  XRADX_OC  =  8000.  
  XRADY_OC  =  8000.
END IF
!
!-------------------------------------------------------------------------------
!
!
!*      16.   SET DEFAULT VALUES FOR MODD_PARAM_ICE :
!             ---------------------------------------
!
IF (KMI == 1) THEN
  LRED    = .TRUE.
  LWARM = .TRUE.
  CPRISTINE_ICE = 'PLAT'
  LSEDIC  = .TRUE.
  LCONVHG = .FALSE.
  CSEDIM  = 'SPLI'
  LFEEDBACKT = .TRUE.
  LEVLIMIT = .TRUE.
  LNULLWETG = .TRUE.
  LWETGPOST = .TRUE.
  LNULLWETH = .TRUE.
  LWETHPOST = .TRUE.
  CSNOWRIMING = 'M90 '
  CSUBG_RC_RR_ACCR = 'NONE'
  CSUBG_RR_EVAP = 'NONE'
  CSUBG_PR_PDF = 'SIGM'
  XFRACM90 = 0.1
  LCRFLIMIT = .TRUE.
  NMAXITER = 5
  XMRSTEP = 0.00005
  XTSTEP_TS = 0.
  LADJ_BEFORE = .TRUE.
  LADJ_AFTER = .TRUE.
  CFRAC_ICE_ADJUST = 'S'
  XSPLIT_MAXCFL = 0.8
  CFRAC_ICE_SHALLOW_MF = 'S'
  LSEDIM_AFTER = .FALSE.
  LDEPOSC = .FALSE.
  XVDEPOSC= 0.02 ! 2 cm/s
  LSNOW_T=.FALSE.
  LPACK_INTERP=.TRUE.
  LPACK_MICRO=.TRUE. ! Meso-NH does not work with LPACK_MICRO=.FALSE.
END IF
!
!-------------------------------------------------------------------------------
!
!
!*      17.   SET DEFAULT VALUES FOR MODD_PARAM_KAFR_n :
!             --------------------------------------------
!
XDTCONV       = MAX( 300.0,XTSTEP )
NICE          = 1
LREFRESH_ALL  = .TRUE.
LCHTRANS      = .FALSE.
LDOWN         = .TRUE.
LSETTADJ      = .FALSE.
XTADJD        = 3600.
XTADJS        = 10800.
LDIAGCONV     = .FALSE.
NENSM         = 0
!
!-------------------------------------------------------------------------------
!
!
!*      18.   SET DEFAULT VALUES FOR MODD_PARAM_MFSHALL_n :
!             --------------------------------------------
!
XIMPL_MF    = 1.
CMF_UPDRAFT = 'EDKF'
CMF_CLOUD   = 'DIRE'
LMIXUV      = .TRUE. 
LMF_FLX     = .FALSE.
!
XALP_PERT   = 0.3  
XABUO       = 1.     
XBENTR      = 1.   
XBDETR      = 0.   
XCMF        = 0.065 
XENTR_MF    = 0.035
XCRAD_MF    = 50.  
XENTR_DRY   = 0.55 
XDETR_DRY   = 10.  
XDETR_LUP   = 1.  
XKCF_MF     = 2.75
XKRC_MF     = 1.   
XTAUSIGMF   = 600.  
XPRES_UV    = 0.5  
XFRAC_UP_MAX= 0.33
XALPHA_MF = 2.     
XSIGMA_MF = 20.  
!
XA1    =  2./3.  
XB     =  0.002       
XC     =  0.012     
XBETA1 =  0.9 
XR     =  2.
XLAMBDA_MF=  0.
LGZ    =  .FALSE.
XGZ    =  1.83 ! between 1.83 and 1.33
!
!-------------------------------------------------------------------------------
!
!*      19.   SET DEFAULT VALUES FOR MODD_PARAM_C2R2 :
!             ----------------------------------------
!
IF (KMI == 1) THEN
  XNUC    = 1.0
  XALPHAC = 3.0
  XNUR    = 2.0
  XALPHAR = 1.0
!
  LRAIN   = .TRUE.
  LSEDC   = .TRUE.
  LACTIT  = .FALSE.
  LSUPSAT = .FALSE.
  LDEPOC  = .FALSE.
  XVDEPOC = 0.02 ! 2 cm/s
  LACTTKE = .TRUE.
!
  HPARAM_CCN = 'XXX'
  HINI_CCN   = 'XXX'
  HTYPE_CCN  = 'X'
!
  XCHEN      = 0.0
  XKHEN      = 0.0
  XMUHEN     = 0.0
  XBETAHEN   = 0.0
!
  XCONC_CCN   = 0.0
  XAERDIFF    = 0.0
  XAERHEIGHT  = 2000
  XR_MEAN_CCN = 0.0
  XLOGSIG_CCN = 0.0
  XFSOLUB_CCN = 1.0
  XACTEMP_CCN = 280.
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      19.BIS   SET DEFAULT VALUES FOR MODD_PARAM_LIMA :
!                ----------------------------------------
!
IF (KMI == 1) THEN
   LPTSPLIT     = .TRUE.
   L_LFEEDBACKT = .TRUE.
   L_NMAXITER   = 1
   L_XMRSTEP    = 0.
   L_XTSTEP_TS  = 0.
!
   YNUC    = 1.0
   YALPHAC = 3.0
   YNUR    = 2.0
   YALPHAR = 1.0
!
   OWARM  = .TRUE.
   LACTI  = .TRUE.
   ORAIN  = .TRUE.
   OSEDC  = .TRUE.
   OACTIT = .FALSE.
   LADJ   = .TRUE.
   LSPRO  = .FALSE.
   LKHKO  = .FALSE.
   ODEPOC = .TRUE.
   LBOUND = .FALSE.
   OACTTKE = .TRUE.
   LKESSLERAC = .FALSE.
!
   NMOM_C = 2
   NMOM_R = 2
!
   OVDEPOC = 0.02 ! 2 cm/s
!
   CINI_CCN   = 'AER'
   CTYPE_CCN(:) = 'M'
!
   YAERDIFF    = 0.0
   YAERHEIGHT  = 2000.
!  YR_MEAN_CCN = 0.0   ! In case of 'CCN' initialization
!  YLOGSIG_CCN = 0.0
   YFSOLUB_CCN = 1.0
   YACTEMP_CCN = 280.
!
   NMOD_CCN = 1
!
!* AP Scavenging
!
   LSCAV      = .FALSE.
   LAERO_MASS = .FALSE.
!
   LCCN_HOM = .TRUE.
   CCCN_MODES = 'COPT'
   XCCN_CONC(:)=300.
!
   LHHONI = .FALSE.
   LCOLD  = .TRUE.
   LNUCL  = .TRUE.
   LSEDI  = .TRUE.
   LSNOW  = .TRUE.
   LHAIL  = .FALSE.
   YSNOW_T = .FALSE.
   LMURAKAMI = .TRUE.
   CPRISTINE_ICE_LIMA = 'PLAT'
   CHEVRIMED_ICE_LIMA = 'GRAU'
   XFACTNUC_DEP = 1.0  
   XFACTNUC_CON = 1.0
   NMOM_I = 2
   NMOM_S = 1
   NMOM_G = 1
   NMOM_H = 1
   NMOD_IFN = 1
   NIND_SPECIE = 1
   LMEYERS = .FALSE.
   LIFN_HOM = .TRUE.
   CIFN_SPECIES = 'PHILLIPS'
   CINT_MIXING = 'DM2'
   XIFN_CONC(:) = 100.
   NMOD_IMM = 0
   NPHILLIPS=8
   LCIBU = .FALSE.
   XNDEBRIS_CIBU = 50.0
   LRDSF = .FALSE.
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      20.   SET DEFAULT VALUES FOR MODD_CH_MNHC_n
!             -------------------------------------
!
LUSECHEM            = .FALSE.
LUSECHAQ            = .FALSE.
LUSECHIC            = .FALSE.
LCH_INIT_FIELD      = .FALSE.
LCH_CONV_SCAV       = .FALSE.
LCH_CONV_LINOX      = .FALSE.
LCH_PH              = .FALSE.
LCH_RET_ICE         = .FALSE.
XCH_PHINIT          = 5.2
XRTMIN_AQ           = 5.e-8
CCHEM_INPUT_FILE    = 'EXSEG1.nam'
CCH_TDISCRETIZATION = 'SPLIT'
NCH_SUBSTEPS        = 1
LCH_TUV_ONLINE      = .FALSE.
CCH_TUV_LOOKUP      = 'PHOTO.TUV39'
CCH_TUV_CLOUDS      = 'NONE'
XCH_TUV_ALBNEW      = -1.
XCH_TUV_DOBNEW      = -1.
XCH_TUV_TUPDATE     = 600.
CCH_VEC_METHOD      = 'MAX'
NCH_VEC_LENGTH      = 50
XCH_TS1D_TSTEP      = 600.
CCH_TS1D_COMMENT    = 'no comment'
CCH_TS1D_FILENAME   = 'IO1D'
CSPEC_PRODLOSS      = ''
CSPEC_BUDGET        = ''
!
!-------------------------------------------------------------------------------
!
!*      21.   SET DEFAULT VALUES FOR MODD_SERIES AND MODD_SERIE_n
!             ---------------------------------------------------
!
IF (KMI == 1) THEN
  LSERIES      = .FALSE.
  LMASKLANDSEA = .FALSE.
  LWMINMAX     = .FALSE.
  LSURF        = .FALSE.
ENDIF
!
NIBOXL = 1 !+ JPHEXT
NIBOXH = 1 !+ 2*JPHEXT
NJBOXL = 1 !+ JPHEXT
NJBOXH = 1 !+ 2*JPHEXT
NKCLS  = 1 !+ JPVEXT
NKLOW  = 1 !+ JPVEXT
NKMID  = 1 !+ JPVEXT
NKUP   = 1 !+ JPVEXT
NKCLA  = 1 !+ JPVEXT
NBJSLICE = 1
NJSLICEL(:) = 1 !+ JPHEXT
NJSLICEH(:) = 1 !+ 2*JPHEXT
NFREQSERIES  = INT(XSEGLEN /(100.*XTSTEP) )
NFREQSERIES  = MAX(NFREQSERIES,1)
!
!-------------------------------------------------------------------------------
!
!*      22.   SET DEFAULT VALUES FOR MODD_TURB_CLOUD
!             --------------------------------------
!
IF (KMI == 1) THEN
  NMODEL_CLOUD = NUNDEF
  CTURBLEN_CLOUD = 'DELT'
  XCOEF_AMPL_SAT = 5.
  XCEI_MIN = 0.001E-06
  XCEI_MAX = 0.01E-06
ENDIF
!-------------------------------------------------------------------------------
!
!*      22.   SET DEFAULT VALUES FOR MODD_MEAN_FIELD
!             --------------------------------------
!
IF (KMI == 1) THEN
  LMEAN_FIELD = .FALSE.
  LCOV_FIELD  = .FALSE.
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      22.   SET DEFAULT VALUES FOR MODD_AEROSOL
!             -----------------------------------
IF (KMI == 1) THEN ! other values are defined in modd_ch_aerosol
!
! aerosol lognormal parameterization

LVARSIGI  = .FALSE.   ! switch to active pronostic dispersion for I mode
LVARSIGJ  = .FALSE.   ! switch to active pronostic dispersion for J mode
LHETEROSO4 = .FALSE.  ! switch to active sulfates heteronegeous
                      ! production
LSEDIMAERO = .FALSE.  ! switch to active aerosol sedimentation
LAERINIT   = .FALSE.  ! switch to initialize aerosol in arome
CMINERAL      = "NONE"   ! mineral equilibrium scheme
CORGANIC      = "NONE"   ! mineral equilibrium scheme
CNUCLEATION   = "NONE" ! sulfates nucleation scheme
LDEPOS_AER(:) = .FALSE.

ENDIF

!*      23.   SET DEFAULT VALUES FOR MODD_DUST and MODD_SALT
!             ----------------------------------------------
!
IF (KMI == 1) THEN ! other values initialized in modd_dust
  LDUST      = .FALSE.
  NMODE_DST  = 3
  LVARSIG    = .FALSE.
  LSEDIMDUST = .FALSE.
  LDEPOS_DST(:) = .FALSE.

  LSALT      = .FALSE.
  LVARSIG_SLT= .FALSE.
  LSEDIMSALT = .FALSE.
  LDEPOS_SLT(:)     = .FALSE.
ENDIF
!
!-------------------------------------------------------------------------------
!
!
!*      24.   SET DEFAULT VALUES FOR MODD_PASPOL           
!             ----------------------------------
!
! other values initialized in modd_paspol
!
IF (KMI == 1) THEN 
  LPASPOL    = .FALSE.
  NRELEASE   = 0
  CPPINIT(:) ='1PT'
  XPPLAT(:)  = 0.
  XPPLON (:) = 0.
  XPPMASS(:) = 0.
  XPPBOT(:)  = 0.
  XPPTOP(:)  = 0.
  CPPT1(:)   = "20010921090000"
  CPPT2(:)   = "20010921090000"
  CPPT3(:)   = "20010921091500"
  CPPT4(:)   = "20010921091500"
ENDIF
!
!-------------------------------------------------------------------------------
!
!
!*      25.   SET DEFAULT VALUES FOR MODD_CONDSAMP         
!             ----------------------------------
!
! other values initialized in modd_condsamp
!
IF (KMI == 1) THEN 
  LCONDSAMP  = .FALSE.
  NCONDSAMP  = 3
  XRADIO(:)  = 900.
  XSCAL(:)   = 1.
  XHEIGHT_BASE = 100.
  XDEPTH_BASE  = 100.
  XHEIGHT_TOP  = 100.
  XDEPTH_TOP   = 100.
  NFINDTOP     = 0
  XTHVP        = 0.25
  LTPLUS       = .TRUE.
ENDIF  
!-------------------------------------------------------------------------------
!
!
!*      26.   SET DEFAULT VALUES FOR MODD_LATZ_EDFLX      
!             ----------------------------------
!
IF (KMI == 1) THEN 
  LUV_FLX=.FALSE.
  XUV_FLX1=3.E+14
  XUV_FLX2=0.
  LTH_FLX=.FALSE.
  XTH_FLX=0.75
ENDIF  
#ifdef MNH_FOREFIRE
!-------------------------------------------------------------------------------
!
!*      27.   SET DEFAULT VALUES FOR MODD_FOREFIRE         
!             ----------------------------------
!
! other values initialized in modd_forefire
!
IF (KMI == 1) THEN 
  LFOREFIRE = .FALSE.
  LFFCHEM = .FALSE.
  COUPLINGRES = 100.
  NFFSCALARS = 0
ENDIF  
#endif                 
!-------------------------------------------------------------------------------
!
!*      28.   SET DEFAULT VALUES FOR MODD_BLOWSNOW AND  MODD_BLOWSNOW_n       
!             ----------------------------------------
! 
IF (KMI == 1) THEN
   LBLOWSNOW  = .FALSE.
   XALPHA_SNOW  = 3.
   XRSNOW       = 4.
   CSNOWSEDIM  = 'TABC'
END IF
LSNOWSUBL = .FALSE.
!
!
!-------------------------------------------------------------------------------
!
!*      29.   SET DEFAULT VALUES FOR MODD_VISC           
!             ----------------------------------
!
! other values initialized in modd_VISC
!
IF (KMI == 1) THEN
  LVISC    = .FALSE.
  LVISC_UVW    = .FALSE.
  LVISC_TH    = .FALSE.
  LVISC_SV    = .FALSE.
  LVISC_R    = .FALSE.
  XMU_V   = 0.
  XPRANDTL = 0.
ENDIF
!
!-------------------------------------------------------------------------------
!
!
!*      30.   SET DEFAULT VALUES FOR MODD_DRAG           
!             ----------------------------------
!
! other values initialized in modd_DRAG
!
IF (KMI == 1) THEN
  LDRAG    = .FALSE.
  LMOUNT   = .FALSE.
  NSTART = 1
  XHSTART = 0.
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      31.   SET DEFAULT VALUES FOR MODD_IBM_PARAMn         
!             --------------------------------------
!                                                           
  LIBM           = .FALSE.
  LIBM_TROUBLE   = .FALSE.
  CIBM_ADV       = 'NOTHIN'
  XIBM_EPSI      = 1.E-9
  XIBM_IEPS      = 1.E+9
  NIBM_ITR       = 8  
  XIBM_RUG       = 0.01    ! (m^1.s^-0)
  XIBM_VISC      = 1.56e-5 ! (m^2.s^-1)
  XIBM_CNU       = 0.06    ! (m^0.s^-0)  

  NIBM_LAYER_P = 2
  NIBM_LAYER_Q = 2
  NIBM_LAYER_R = 2
  NIBM_LAYER_S = 2
  NIBM_LAYER_T = 2
  NIBM_LAYER_E = 2
  NIBM_LAYER_V = 2

  XIBM_RADIUS_P = 2.
  XIBM_RADIUS_Q = 2.
  XIBM_RADIUS_R = 2.
  XIBM_RADIUS_S = 2.
  XIBM_RADIUS_T = 2.
  XIBM_RADIUS_E = 2.
  XIBM_RADIUS_V = 2.

  XIBM_POWERS_P = 1.
  XIBM_POWERS_Q = 1.
  XIBM_POWERS_R = 1.
  XIBM_POWERS_S = 1.
  XIBM_POWERS_T = 1.
  XIBM_POWERS_E = 1.
  XIBM_POWERS_V = 1.

  CIBM_MODE_INTE3_P = 'LAI'
  CIBM_MODE_INTE3_Q = 'LAI'
  CIBM_MODE_INTE3_R = 'LAI'
  CIBM_MODE_INTE3_S = 'LAI'
  CIBM_MODE_INTE3_T = 'LAI'
  CIBM_MODE_INTE3_E = 'LAI'
  CIBM_MODE_INTE3_V = 'LAI'

  CIBM_MODE_INTE1_P = 'CL2'
  CIBM_MODE_INTE1_Q = 'CL2'
  CIBM_MODE_INTE1_R = 'CL2'
  CIBM_MODE_INTE1_S = 'CL2'
  CIBM_MODE_INTE1_T = 'CL2'
  CIBM_MODE_INTE1_E = 'CL2'
  CIBM_MODE_INTE1NV = 'CL2'
  CIBM_MODE_INTE1TV = 'CL2'
  CIBM_MODE_INTE1CV = 'CL2'

  CIBM_MODE_BOUND_P = 'SYM'
  CIBM_MODE_BOUND_Q = 'SYM'
  CIBM_MODE_BOUND_R = 'SYM'
  CIBM_MODE_BOUND_S = 'SYM'
  CIBM_MODE_BOUND_T = 'SYM'
  CIBM_MODE_BOUND_E = 'SYM'
  CIBM_MODE_BOUNT_V = 'ASY'
  CIBM_MODE_BOUNN_V = 'ASY'
  CIBM_MODE_BOUNC_V = 'ASY'

  XIBM_FORC_BOUND_P = 0.
  XIBM_FORC_BOUND_Q = 0.
  XIBM_FORC_BOUND_R = 0.
  XIBM_FORC_BOUND_S = 0. 
  XIBM_FORC_BOUND_T = 0.
  XIBM_FORC_BOUND_E = 0.
  XIBM_FORC_BOUNN_V = 0.
  XIBM_FORC_BOUNT_V = 0.
  XIBM_FORC_BOUNC_V = 0.

  CIBM_TYPE_BOUND_P = 'NEU'
  CIBM_TYPE_BOUND_Q = 'NEU'
  CIBM_TYPE_BOUND_R = 'NEU'
  CIBM_TYPE_BOUND_S = 'NEU'
  CIBM_TYPE_BOUND_T = 'NEU'
  CIBM_TYPE_BOUND_E = 'NEU'
  CIBM_TYPE_BOUNT_V = 'DIR'
  CIBM_TYPE_BOUNN_V = 'DIR'
  CIBM_TYPE_BOUNC_V = 'DIR'

  CIBM_FORC_BOUND_P = 'CST'
  CIBM_FORC_BOUND_Q = 'CST'
  CIBM_FORC_BOUND_R = 'CST'
  CIBM_FORC_BOUND_S = 'CST'
  CIBM_FORC_BOUND_T = 'CST'
  CIBM_FORC_BOUND_E = 'CST'
  CIBM_FORC_BOUNN_V = 'CST'
  CIBM_FORC_BOUNT_V = 'CST'
  CIBM_FORC_BOUNC_V = 'CST'
  CIBM_FORC_BOUNR_V = 'CST'

!
!-------------------------------------------------------------------------------
!
!*      32.   SET DEFAULT VALUES FOR MODD_RECYCL_PARAMn         
!             --------------------------------------
!
  LRECYCL  = .FALSE.
  LRECYCLN = .FALSE.
  LRECYCLW = .FALSE.
  LRECYCLE = .FALSE.
  LRECYCLS = .FALSE.
  XDRECYCLN  = 0.
  XARECYCLN  = 0.
  XDRECYCLW  = 0.
  XARECYCLW  = 0.
  XDRECYCLS  = 0.
  XARECYCLS  = 0.
  XDRECYCLE  = 0.
  XARECYCLE  = 0.
  XTMOY      = 0.
  XTMOYCOUNT = 0.
  XNUMBELT   = 28.
  XRCOEFF    = 0.2
  XTBVTOP    = 500.
  XTBVBOT    = 300.
!
!-------------------------------------------------------------------------------
!
!*      33.   SET DEFAULT VALUES FOR MODD_FIRE
!             --------------------------------
!
! Blaze fire model namelist
!
IF (KMI == 1) THEN
  LBLAZE = .FALSE.	            ! Flag for Fire model use, default FALSE
  !
  CPROPAG_MODEL = 'SANTONI2011' ! Fire propagation model (default SANTONI2011)
  !
  CHEAT_FLUX_MODEL   = 'EXS'	  ! Sensible heat flux injection model (default EXS)
  CLATENT_FLUX_MODEL = 'EXP'	  ! latent heat flux injection model (default EXP)
  XFERR = 0.8                   ! Energy released in flamming stage (only for EXP)
  !
  CFIRE_CPL_MODE = '2WAYCPL'    ! Coupling mode (default 2way coupled)
  CBMAPFILE = CINIFILE          ! File name of BMAP for FIR2ATM mode
  LINTERPWIND = .TRUE.          ! Horizontal interpolation of wind
  LSGBAWEIGHT = .FALSE.         ! Flag for use of weighted average method for SubGrid Burning Area computation
  !
  NFIRE_WENO_ORDER = 3	        ! Weno order (1,3,5)
  NFIRE_RK_ORDER = 3	          ! Runge Kutta order (1,2,3,4)
  !
  NREFINX = 1	                  ! Refinement ratio X
  NREFINY	= 1                   ! Refinement ratio Y
  !
  XCFLMAXFIRE = 0.8             ! Max CFL on fire mesh
  XLSDIFFUSION = 0.1            ! Numerical diffusion of LevelSet
  XROSDIFFUSION = 0.05          ! Numerical diffusion of ROS
  !
  XFLUXZEXT = 3.                ! Flux distribution on vertical caracteristic length
  XFLUXZMAX = 4. * XFLUXZEXT    ! Flux distribution on vertical max injetion height
  !
  XFLXCOEFTMP = 1.              ! Flux multiplicator. For testing
  !
  LWINDFILTER = .FALSE.         ! Fire wind filtering flag
  CWINDFILTER = 'EWAM'	        ! Wind filter method (EWAM or WLIM)
  XEWAMTAU    = 20.             ! Time averaging constant for EWAM method (s)
  XWLIMUTH    = 8.              ! Thresehold wind value for WLIM method (m/s)
  XWLIMUTMAX  = 9.              ! Maximum wind value for WLIM method (m/s) (needs to be >= XWLIMUTH )
  !
  NNBSMOKETRACER = 1            ! Nb of smoke tracers
  !
  NWINDSLOPECPLMODE = 0         ! Flag for use of wind/slope in ROS (0 = wind + slope, 1 = wind only, 2 = slope only (U0=0))
  !
  !
  !
  !! DO NOT CHANGE BELOW PARAMETERS
  XFIREMESHSIZE(:) = 0.         ! Fire mesh size (dxf,dyf)
  LRESTA_ASE = .FALSE.          ! Flag for using ASE in RESTA file
  LRESTA_AWC = .FALSE.          ! Flag for using AWC in RESTA file
  LRESTA_EWAM = .FALSE.         ! Flag for using EWAM in RESTA file
  LRESTA_WLIM = .FALSE.         ! Flag for using WLIM in RESTA file
ENDIF

!-------------------------------------------------------------------------------
END SUBROUTINE DEFAULT_DESFM_n
