!MNH_LIC Copyright 1994-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_READ_FIELD
!     ######################
!
INTERFACE 
!
      SUBROUTINE READ_FIELD(KOCEMI,TPINIFILE,KIU,KJU,KKU,                    &
            HGETTKET,HGETRVT,HGETRCT,HGETRRT,HGETRIT,HGETCIT,HGETZWS,        &
            HGETRST,HGETRGT,HGETRHT,HGETSVT,HGETSRCT,HGETSIGS,HGETCLDFR,HGETICEFR, &
            HGETBL_DEPTH,HGETSBL_DEPTH,HGETPHC,HGETPHR,HUVW_ADV_SCHEME,      &
            HTEMP_SCHEME,KSIZELBX_ll,KSIZELBXU_ll,KSIZELBY_ll,KSIZELBYV_ll,  &
            KSIZELBXTKE_ll,KSIZELBYTKE_ll,                                   &
            KSIZELBXR_ll,KSIZELBYR_ll,KSIZELBXSV_ll,KSIZELBYSV_ll,           &
            PUM,PVM,PWM,PDUM,PDVM,PDWM,                                      &
            PUT,PVT,PWT,PTHT,PPABST,PTKET,PRTKEMS,                           &
            PRT,PSVT,PZWS,PCIT,PDRYMASST,PDRYMASSS,                          &            
            PSIGS,PSRCT,PCLDFR,PICEFR,PBL_DEPTH,PSBL_DEPTH,PWTHVMF,PPHC,PPHR, &
            PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM, PLSZWSM,                        &
            PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,            &
            PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,            &
            KFRC,TPDTFRC,PUFRC,PVFRC,PWFRC,PTHFRC,PRVFRC,                    &
            PTENDTHFRC,PTENDRVFRC,PGXTHFRC,PGYTHFRC,PPGROUNDFRC,PATC,        &
            PTENDUFRC,PTENDVFRC,                                             &
            KADVFRC,TPDTADVFRC,PDTHFRC,PDRVFRC,                              &
            KRELFRC,TPDTRELFRC, PTHREL, PRVREL,                              &
            PVTH_FLUX_M,PWTH_FLUX_M,PVU_FLUX_M,                              &
            PRUS_PRES,PRVS_PRES,PRWS_PRES,PRTHS_CLD,PRRS_CLD,PRSVS_CLD,      &
            PIBM_LSF,PIBM_XMUT,PUMEANW,PVMEANW,PWMEANW,PUMEANN,PVMEANN,      &
            PWMEANN,PUMEANE,PVMEANE,PWMEANE,PUMEANS,PVMEANS,PWMEANS,         &
            PLSPHI,PBMAP,PFMASE,PFMAWC,PFMWINDU,PFMWINDV,PFMWINDW,PFMHWS     )
!
USE MODD_IO, ONLY : TFILEDATA
USE MODD_TIME ! for type DATE_TIME
!
!
INTEGER,                   INTENT(IN)  :: KOCEMI !Ocan model index
TYPE(TFILEDATA),           INTENT(IN)  :: TPINIFILE    !Initial file
INTEGER,                   INTENT(IN)  :: KIU, KJU, KKU
                             ! array sizes in x, y and z  directions
! 
CHARACTER (LEN=*),         INTENT(IN)  :: HGETTKET,                          &
                                          HGETRVT,HGETRCT,HGETRRT,           &
                                          HGETRIT,HGETRST,HGETRGT,HGETRHT,   & 
                                          HGETCIT,HGETSRCT, HGETZWS,         &
                                          HGETSIGS, HGETCLDFR, HGETICEFR,    &
                                          HGETBL_DEPTH, HGETSBL_DEPTH,       &
                                          HGETPHC, HGETPHR
CHARACTER (LEN=*), DIMENSION(:),INTENT(IN)  :: HGETSVT
!
! GET indicators to know wether a given  variable should or not be read in the
! FM file at time t-deltat and t
CHARACTER(LEN=6),         INTENT(IN)    :: HUVW_ADV_SCHEME ! advection scheme for wind
CHARACTER(LEN=4),         INTENT(IN)    :: HTEMP_SCHEME ! advection scheme for wind
!
! sizes of the West-east total LB area
INTEGER, INTENT(IN) :: KSIZELBX_ll,KSIZELBXU_ll      ! for T,V,W and u 
INTEGER, INTENT(IN) :: KSIZELBXTKE_ll                ! for TKE 
INTEGER, INTENT(IN) :: KSIZELBXR_ll,KSIZELBXSV_ll    ! for Rx and SV    
! sizes of the North-south total LB area
INTEGER, INTENT(IN) :: KSIZELBY_ll,KSIZELBYV_ll      ! for T,U,W  and v
INTEGER, INTENT(IN) :: KSIZELBYTKE_ll                ! for TKE
INTEGER, INTENT(IN) :: KSIZELBYR_ll,KSIZELBYSV_ll    ! for Rx and SV 
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PUM,PVM,PWM     ! U,V,W at t-dt
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PDUM,PDVM,PDWM  ! Difference on U,V,W 
                                                          ! between t+dt and t-dt
REAL, DIMENSION(:,:),      INTENT(OUT) :: PBL_DEPTH       ! BL depth
REAL, DIMENSION(:,:),      INTENT(OUT) :: PSBL_DEPTH      ! SBL depth
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PWTHVMF         ! MassFlux buoyancy flux
!
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PUT,PVT,PWT     ! U,V,W at t
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PTHT,PTKET      ! theta, tke and
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PRTKEMS         ! tke adv source
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PPABST          ! pressure at t
REAL, DIMENSION(:,:,:,:),  INTENT(OUT) :: PRT,PSVT        ! moist and scalar
                                                          ! variables at t
REAL, DIMENSION(:,:),      INTENT(INOUT) :: PZWS
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PSRCT           ! turbulent flux
                                                          !  <s'Rc'> at t 
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PCIT            ! ice conc. at t
REAL,                      INTENT(OUT) :: PDRYMASST       ! Md(t)
REAL,                      INTENT(OUT) :: PDRYMASSS       ! d Md(t) / dt
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PSIGS           ! =sqrt(<s's'>) for the
                                                          ! Subgrid Condensation
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PCLDFR          ! cloud fraction  
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PICEFR          ! cloud fraction  
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PPHC            ! pH value in cloud water  
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PPHR            ! pH value in rainwater  
! Larger Scale fields
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLSUM,PLSVM,PLSWM    ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLSTHM,  PLSRVM      ! Mass
! LB fields
REAL, DIMENSION(:,:),            INTENT(OUT) :: PLSZWSM              ! significant height of sea waves
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXUM,PLBXVM,PLBXWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYUM,PLBYVM,PLBYWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTKEM          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTKEM
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBXRM  ,PLBXSVM  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBYRM  ,PLBYSVM  ! in x and y-dir.
! Forcing fields
INTEGER,                        INTENT(IN)    :: KFRC              ! number of forcing
TYPE (DATE_TIME), DIMENSION(:), INTENT(OUT)   :: TPDTFRC           ! date of forcing profs.
REAL, DIMENSION(:,:),           INTENT(OUT)   :: PUFRC,PVFRC,PWFRC ! forcing variables
REAL, DIMENSION(:,:),           INTENT(OUT)   :: PTHFRC,PRVFRC
REAL, DIMENSION(:,:),           INTENT(OUT) :: PTENDUFRC,PTENDVFRC
REAL, DIMENSION(:,:),           INTENT(OUT)   :: PTENDTHFRC,PTENDRVFRC,PGXTHFRC,PGYTHFRC
REAL, DIMENSION(:),             INTENT(OUT)   :: PPGROUNDFRC
REAL, DIMENSION(:,:,:,:),       INTENT(OUT)   :: PATC
INTEGER,                        INTENT(IN)    :: KADVFRC           ! number of forcing
TYPE (DATE_TIME), DIMENSION(:), INTENT(OUT)   :: TPDTADVFRC        ! date of forcing profs.
REAL, DIMENSION(:,:,:,:),       INTENT(OUT)   :: PDTHFRC, PDRVFRC
INTEGER,                        INTENT(IN)    :: KRELFRC           ! number of forcing
TYPE (DATE_TIME), DIMENSION(:), INTENT(OUT)   :: TPDTRELFRC        ! date of forcing profs.
REAL, DIMENSION(:,:,:,:),       INTENT(OUT)   :: PTHREL, PRVREL
REAL, DIMENSION(:,:,:),         INTENT(OUT)   :: PVTH_FLUX_M,PWTH_FLUX_M,PVU_FLUX_M ! Eddy fluxes
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PRUS_PRES, PRVS_PRES, PRWS_PRES
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PRTHS_CLD
REAL, DIMENSION(:,:,:,:),       INTENT(INOUT) :: PRRS_CLD, PRSVS_CLD
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PIBM_LSF,PIBM_XMUT
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PUMEANW,PVMEANW,PWMEANW
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PUMEANN,PVMEANN,PWMEANN
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PUMEANE,PVMEANE,PWMEANE
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PUMEANS,PVMEANS,PWMEANS
!
! Fire Model fields
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLSPHI    ! Fire Model Level Set function Phi [-]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PBMAP     ! Fire Model Burning map [s]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMASE    ! Fire Model Available Sensible Energy [J/m2]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMAWC    ! Fire Model Available Water Content [kg/m2]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMWINDU  ! Fire Model filtered u wind [m/s]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMWINDV  ! Fire Model filtered v wind [m/s]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMWINDW  ! Fire Model filtered w wind [m/s]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMHWS    ! Fire Model filtered horizontal wind speed [m/s]
!
END SUBROUTINE READ_FIELD
!
END INTERFACE
!
END MODULE MODI_READ_FIELD
!
!     ########################################################################
      SUBROUTINE READ_FIELD(KOCEMI,TPINIFILE,KIU,KJU,KKU,                    &
            HGETTKET,HGETRVT,HGETRCT,HGETRRT,HGETRIT,HGETCIT,HGETZWS,        &
            HGETRST,HGETRGT,HGETRHT,HGETSVT,HGETSRCT,HGETSIGS,HGETCLDFR,HGETICEFR, &
            HGETBL_DEPTH,HGETSBL_DEPTH,HGETPHC,HGETPHR,HUVW_ADV_SCHEME,      &
            HTEMP_SCHEME,KSIZELBX_ll,KSIZELBXU_ll,KSIZELBY_ll,KSIZELBYV_ll,  &
            KSIZELBXTKE_ll,KSIZELBYTKE_ll,                                   &
            KSIZELBXR_ll,KSIZELBYR_ll,KSIZELBXSV_ll,KSIZELBYSV_ll,           &
            PUM,PVM,PWM,PDUM,PDVM,PDWM,                                      &
            PUT,PVT,PWT,PTHT,PPABST,PTKET,PRTKEMS,                           &
            PRT,PSVT,PZWS,PCIT,PDRYMASST,PDRYMASSS,                          &
            PSIGS,PSRCT,PCLDFR,PICEFR,PBL_DEPTH,PSBL_DEPTH,PWTHVMF,PPHC,PPHR, &
            PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM,                         &
            PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,            &
            PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,            &
            KFRC,TPDTFRC,PUFRC,PVFRC,PWFRC,PTHFRC,PRVFRC,                    &
            PTENDTHFRC,PTENDRVFRC,PGXTHFRC,PGYTHFRC,PPGROUNDFRC,PATC,        &
            PTENDUFRC,PTENDVFRC,                                             &
            KADVFRC,TPDTADVFRC,PDTHFRC,PDRVFRC,                              &
            KRELFRC,TPDTRELFRC, PTHREL, PRVREL,                              &
            PVTH_FLUX_M,PWTH_FLUX_M,PVU_FLUX_M,                              &
            PRUS_PRES,PRVS_PRES,PRWS_PRES,PRTHS_CLD,PRRS_CLD,PRSVS_CLD,      &
            PIBM_LSF,PIBM_XMUT,PUMEANW,PVMEANW,PWMEANW,PUMEANN,PVMEANN,      &
            PWMEANN,PUMEANE,PVMEANE,PWMEANE,PUMEANS,PVMEANS,PWMEANS,         &
            PLSPHI,PBMAP,PFMASE,PFMAWC,PFMWINDU,PFMWINDV,PFMWINDW,PFMHWS     )
!     ########################################################################
!
!!****  *READ_FIELD* - routine to read prognostic and surface fields
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize  prognostic and 
!     surface fields by reading their  value in initial file or by setting 
!     them to a fixed value.
!
!!**  METHOD
!!    ------
!!      According to the get indicators, the prognostics fields are :
!!          - initialized by reading their value in the LFIFM file 
!!    if the corresponding indicators are equal to 'READ'  
!!          - initialized to zero if the corresponding indicators 
!!    are equal to 'INIT'
!!          -  not initialized if their corresponding indicators 
!!    are equal to 'SKIP'
!!
!!      In case of time step change, all fields at t-dt are (linearly)
!!    interpolated to get a consistant initial state before the segment 
!!    integration  
!!
!!    EXTERNAL
!!    --------
!!      FMREAD   : to read data in LFIFM file
!!      INI_LS   : to initialize larger scale fields
!!      INI_LB   : to initialize "2D" surfacic LB fields 
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_CONF   : NVERB,CCONF,CPROGRAM
!!
!!      Module MODD_CTURB :  XTKEMIN
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (routine READ_FIELD)
!!      
!!
!!    AUTHOR
!!    ------
!!  	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        15/06/94 
!!      modification    22/11/94  add the pressure function (J.Stein)
!!      modification    22/11/94  add the LS fields         (J.Stein)
!!      modification    06/01/95  add Md(t)                 (J.P.Lafore)
!!                      26/03/95  add EPS var               (J. Cuxart)
!!                      30/06/95  add var related to the Subgrid condensation
!!                                                                   (J.Stein)
!!                      18/08/95  time step change case     (J.P.Lafore)
!!                      01/03/96  add the cloud fraction    (J. Stein)
!!     modification     13/12/95  add fmread of the forcing variables 
!!                                                          (M.Georgelin)
!!     modification     13/02/96  external control of the forcing (J.-P. Pinty)
!!                      11/04/96  add the ice concentration (J.-P. Pinty)
!!                      27/01/97  read ISVR 3D fields of SV (J.-P. Pinty)
!!                      26/02/97  "surfacic" LS fieds  introduction (J.P.Lafore)
!!          (V MASSON)  03/03/97  positivity control for time step change 
!!                      10/04/97  proper treatment of minima for LS-fields (J.P.Lafore)
!!           J. Stein   22/06/97  use the absolute pressure
!!           J. Stein   22/10/97  cleaning + add the LB fields for u,v,w,theta,Rv
!!          P. Bechtold 22/01/98  add SST and surface pressure forcing
!!          V. Ducrocq  14/08/98  //,  remove KIINF,KJINF,KISUP,KJSUP,
!!                                     and introduce INI_LS and INI_LB
!!          J. Stein    22/01/99  add the reading of STORAGE_TYPE to improve
!!                                the START case when the file contains 2
!!                                instants MT
!!          D. Gazen    22/01/01  use MODD_NSV to handle NSV floating indices
!!                                for the current model
!!          V. Masson   01/2004   removes surface (externalization)
!!       J.-P. Pinty    06/05/04  treat NSV_* for C1R3 and ELEC
!!                      05/06     Remove EPS
!!          M. Leriche  04/10     add pH in cloud water and rainwater
!!          M. Leriche  07/10     treat NSV_* for ice phase chemical species
!!          C.Lac       11/11     Suppress all the t-Dt fields
!!          M.Tomasini, 
!!          P. Peyrille   06/12   2D west african monsoon : add reading of ADV forcing and addy fluxes 
!!          C.Lac       03/13     add prognostic supersaturation for C2R2/KHKO
!!          Bosseur & Filippi 07/13 Adds Forefire
!!          M. Leriche  11/14     correct bug in pH initialization
!!          C.Lac       12/14     correction for reproducibility START/RESTA
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!          M. Leriche  02/16     treat gas and aq. chemicals separately
!!          C.Lac        10/16 CEN4TH with RKC4 + Correction on RK loop
!!                   09/2017 Q.Rodier add LTEND_UV_FRC
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  V. Vionnet       07/17:   add blowing snow scheme
!  P. Wautelet    01/2019:  corrected intent of PDUM,PDVM,PDWM (OUT->INOUT)
!  P. Wautelet 13/02/2019: removed PPABSM and PTSTEP dummy arguments (bugfix: PPABSM was intent(OUT))
!  S. Bielli      02/2019:  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 14/03/2019: correct ZWS when variable not present in file
!  M. Leriche  10/06/2019: in restart case read all immersion modes for LIMA
!  B. Vie         06/2020: Add prognostic supersaturation for LIMA
!  F. Auguste     02/2021: add fields necessary for IBM
!  T. Nagel       02/2021: add fields necessary for turbulence recycling
!  JL. Redelsperger 03/2021:  add necessary variables for Ocean LES case
!  A. Costes      12/2021: add Blaze fire model
!  P. Wautelet 04/02/2022: use TSVLIST to manage metadata of scalar variables
!!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_2D_FRC,          ONLY: L2D_ADV_FRC, L2D_REL_FRC
USE MODD_ADV_n,           ONLY: CTEMP_SCHEME, LSPLIT_CFL
USE MODD_BLOWSNOW_n,      ONLY: XSNWCANO
USE MODD_CONF,            ONLY: CCONF, CPROGRAM, L1D, LFORCING, NVERB
USE MODD_CONF_n,          ONLY: IDX_RVT, IDX_RCT, IDX_RRT, IDX_RIT, IDX_RST, IDX_RGT, IDX_RHT
USE MODD_CST,             ONLY: XALPW, XBETAW, XCPD, XGAMW, XMD, XMV, XP00, XRD
USE MODD_TURB_n,          ONLY: XTKEMIN
USE MODD_DYN_n,           ONLY: LOCEAN
use modd_field,           only: tfieldmetadata, tfieldlist, NMNHDIM_NI, NMNHDIM_NJ, NMNHDIM_NOTLISTED, &
                                TYPEDATE, TYPEREAL, TYPELOG, TYPEINT
USE MODD_FIELD_n,         only: XZWS_DEFAULT
USE MODD_FIRE_n,          ONLY: CWINDFILTER, LBLAZE, LRESTA_ASE, LRESTA_AWC, LRESTA_EWAM, LRESTA_WLIM, LWINDFILTER
USE MODD_IBM_PARAM_n,     ONLY: LIBM
USE MODD_IO,              ONLY: TFILEDATA
USE MODD_LATZ_EDFLX,      ONLY: LTH_FLX, LUV_FLX
USE MODD_LUNIT_N,         ONLY: TLUOUT
USE MODD_NSV,             ONLY: NSV, NSV_C2R2BEG, NSV_C2R2END, NSV_CSBEG, NSV_CSEND, &
#ifdef MNH_FOREFIRE
                                NSV_FFBEG, NSV_FFEND,                                &
#endif
                                NSV_PPBEG, NSV_PPEND, NSV_SNW, NSV_USER, TSVLIST
USE MODD_OCEANH,          ONLY: NFRCLT, NINFRT, XSSOLA_T, XSSUFL_T, XSSTFL_T, XSSVFL_T
USE MODD_PARAM_C2R2,      ONLY: LSUPSAT
USE MODD_PARAMETERS,      ONLY: XUNDEF
USE MODD_PARAM_n,         ONLY: CSCONV
USE MODD_RECYCL_PARAM_n,  ONLY: LRECYCLE, LRECYCLN, LRECYCLS, LRECYCLW, NR_COUNT
USE MODD_REF,             ONLY: LCOUPLES
USE MODD_TIME,            ONLY: DATE_TIME
!
use mode_field,           only: Find_field_id_from_mnhname
USE MODE_IO_FIELD_READ,   only: IO_Field_read
USE MODE_MSG
USE MODE_TOOLS,           ONLY: UPCASE
!
USE MODI_INI_LB
USE MODI_INI_LS
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
!
INTEGER,                   INTENT(IN)  :: KOCEMI !Ocan model index
TYPE(TFILEDATA),           INTENT(IN)  :: TPINIFILE    !Initial file
INTEGER,                   INTENT(IN)  :: KIU, KJU, KKU
                             ! array sizes in x, y and z  directions
! 
CHARACTER (LEN=*),         INTENT(IN)  :: HGETTKET,                          &
                                          HGETRVT,HGETRCT,HGETRRT,           &
                                          HGETRIT,HGETRST,HGETRGT,HGETRHT,   & 
                                          HGETCIT,HGETSRCT, HGETZWS,         &
                                          HGETSIGS, HGETCLDFR, HGETICEFR,    &
                                          HGETBL_DEPTH, HGETSBL_DEPTH,       &
                                          HGETPHC, HGETPHR
CHARACTER (LEN=*), DIMENSION(:),INTENT(IN)  :: HGETSVT
!
! GET indicators to know wether a given  variable should or not be read in the
! FM file at time t-deltat and t
!
CHARACTER(LEN=6),         INTENT(IN)    :: HUVW_ADV_SCHEME ! advection scheme for wind
CHARACTER(LEN=4),         INTENT(IN)    :: HTEMP_SCHEME ! advection scheme for wind
!
! sizes of the West-east total LB area
INTEGER, INTENT(IN) :: KSIZELBX_ll,KSIZELBXU_ll      ! for T,V,W and u 
INTEGER, INTENT(IN) :: KSIZELBXTKE_ll                ! for TKE 
INTEGER, INTENT(IN) :: KSIZELBXR_ll,KSIZELBXSV_ll    ! for Rx and SV    
! sizes of the North-south total LB area
INTEGER, INTENT(IN) :: KSIZELBY_ll,KSIZELBYV_ll      ! for T,U,W  and v
INTEGER, INTENT(IN) :: KSIZELBYTKE_ll                ! for TKE
INTEGER, INTENT(IN) :: KSIZELBYR_ll,KSIZELBYSV_ll    ! for Rx and SV 
!
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PUM,PVM,PWM     ! U,V,W at t-dt
REAL, DIMENSION(:,:,:),    INTENT(INOUT) :: PDUM,PDVM,PDWM  ! Difference on U,V,W
                                                          ! between t+dt and t-dt
REAL, DIMENSION(:,:),      INTENT(OUT) :: PBL_DEPTH       ! BL depth
REAL, DIMENSION(:,:),      INTENT(OUT) :: PSBL_DEPTH      ! SBL depth
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PWTHVMF         ! MassFlux buoyancy flux
!
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PUT,PVT,PWT     ! U,V,W at t
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PTHT,PTKET      ! theta, tke and
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PRTKEMS         ! tke adv source
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PPABST          ! pressure at t
REAL, DIMENSION(:,:,:,:),  INTENT(OUT) :: PRT,PSVT        ! moist and scalar
                                                          ! variables at t
REAL, DIMENSION(:,:),      INTENT(INOUT) :: PZWS
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PSRCT           ! turbulent flux
                                                          !  <s'Rc'> at t 
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PCIT            ! ice conc. at t
REAL,                      INTENT(OUT) :: PDRYMASST       ! Md(t)
REAL,                      INTENT(OUT) :: PDRYMASSS       ! d Md(t) / dt
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PSIGS           ! =sqrt(<s's'>) for the
                                                          ! Subgrid Condensation
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PCLDFR          ! cloud fraction  
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PICEFR          ! cloud fraction  
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PPHC            ! pH value in cloud water  
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PPHR            ! pH value in rainwater  
!
!
! Larger Scale fields
REAL, DIMENSION(:,:),            INTENT(OUT) :: PLSZWSM              ! significant height of sea waves
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLSUM,PLSVM,PLSWM    ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLSTHM,  PLSRVM      ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXUM,PLBXVM,PLBXWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYUM,PLBYVM,PLBYWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTKEM          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTKEM
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBXRM  ,PLBXSVM  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBYRM  ,PLBYSVM  ! in x and y-dir.
!
!       
! Forcing fields
INTEGER,                        INTENT(IN)    :: KFRC              ! number of forcing
TYPE (DATE_TIME), DIMENSION(:), INTENT(OUT)   :: TPDTFRC           ! date of forcing profs.
REAL, DIMENSION(:,:),           INTENT(OUT)   :: PUFRC,PVFRC,PWFRC ! forcing variables
REAL, DIMENSION(:,:),           INTENT(OUT)   :: PTHFRC,PRVFRC
REAL, DIMENSION(:,:),           INTENT(OUT) :: PTENDUFRC,PTENDVFRC
REAL, DIMENSION(:,:),           INTENT(OUT)   :: PTENDTHFRC,PTENDRVFRC,PGXTHFRC,PGYTHFRC
REAL, DIMENSION(:),             INTENT(OUT)   :: PPGROUNDFRC
REAL, DIMENSION(:,:,:,:),       INTENT(OUT)   :: PATC
INTEGER,                        INTENT(IN)    :: KADVFRC           ! number of forcing
TYPE (DATE_TIME), DIMENSION(:), INTENT(OUT)   :: TPDTADVFRC        ! date of forcing profs.
REAL, DIMENSION(:,:,:,:),       INTENT(OUT)   :: PDTHFRC, PDRVFRC
INTEGER,                        INTENT(IN)    :: KRELFRC           ! number of forcing
TYPE (DATE_TIME), DIMENSION(:), INTENT(OUT)   :: TPDTRELFRC        ! date of forcing profs.
REAL, DIMENSION(:,:,:,:),       INTENT(OUT)   :: PTHREL, PRVREL
REAL, DIMENSION(:,:,:),         INTENT(OUT)   :: PVTH_FLUX_M,PWTH_FLUX_M,PVU_FLUX_M ! Eddy fluxes
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PRUS_PRES, PRVS_PRES, PRWS_PRES
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PRTHS_CLD
REAL, DIMENSION(:,:,:,:),       INTENT(INOUT) :: PRRS_CLD, PRSVS_CLD
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PIBM_LSF          ! LSF for IBM
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PIBM_XMUT         ! Turbulent viscosity
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PUMEANW,PVMEANW,PWMEANW ! Velocity average at West boundary
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PUMEANN,PVMEANN,PWMEANN ! Velocity average at North boundary
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PUMEANE,PVMEANE,PWMEANE ! Velocity average at East boundary
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PUMEANS,PVMEANS,PWMEANS ! Velocity average at South boundary
! Fire Model fields
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLSPHI    ! Fire Model Level Set function Phi [-]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PBMAP     ! Fire Model Burning map [s]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMASE    ! Fire Model Available Sensible Energy [J/m2]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMAWC    ! Fire Model Available Water Content [kg/m2]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMWINDU  ! Fire Model filtered u wind [m/s]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMWINDV  ! Fire Model filtered v wind [m/s]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMWINDW  ! Fire Model filtered v wind [m/s]
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PFMHWS    ! Fire Model filtered horizontal wind speed [m/s]
!
!*       0.2   declarations of local variables
!
INTEGER                      :: IID
INTEGER                      :: ILUOUT       ! Unit number for prints
INTEGER                      :: IRESP
INTEGER                      :: ISV          ! total number of  scalar variables
INTEGER                      :: JSV          ! Loop index for additional scalar variables
INTEGER                      :: JKLOOP,JRR   ! Loop indexes
INTEGER                      :: IIUP,IJUP    ! size  of working window arrays
INTEGER                      :: JT           ! loop index
LOGICAL                      :: GLSOURCE     ! switch for the source term (for ini_ls and ini_lb)
LOGICAL                      :: ZLRECYCL     ! switch if turbulence recycling is activated
LOGICAL                      :: GOLDFILEFORMAT
CHARACTER(LEN=3)             :: YFRC         ! To mark the different forcing dates
CHARACTER(LEN=3)             :: YNUM3
CHARACTER(LEN=15)            :: YVAL
REAL, DIMENSION(KIU,KJU,KKU) :: ZWORK        ! to compute supersaturation
TYPE(TFIELDMETADATA)         :: TZFIELD
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATION
!              ---------------
!
GLSOURCE=.FALSE.
ZWORK = 0.0
!
!If TPINIFILE file was written with a MesoNH version < 5.6, some variables had different names or were not available
GOLDFILEFORMAT = (        TPINIFILE%NMNHVERSION(1) < 5                                       &
                   .OR. ( TPINIFILE%NMNHVERSION(1) == 5 .AND. TPINIFILE%NMNHVERSION(2) < 6 ) )
!-------------------------------------------------------------------------------
!
!*       2.    READ PROGNOSTIC VARIABLES
!              -------------------------
!
!*       2.1  Time t:
!
IF (TPINIFILE%NMNHVERSION(1)<5) THEN
  CALL FIND_FIELD_ID_FROM_MNHNAME('UT',IID,IRESP)
  TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
  TZFIELD%CMNHNAME = 'UM'
  CALL IO_Field_read(TPINIFILE,TZFIELD,PUT)
  !
  CALL FIND_FIELD_ID_FROM_MNHNAME('VT',IID,IRESP)
  TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
  TZFIELD%CMNHNAME = 'VM'
  CALL IO_Field_read(TPINIFILE,TZFIELD,PVT)
  !
  CALL FIND_FIELD_ID_FROM_MNHNAME('WT',IID,IRESP)
  TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
  TZFIELD%CMNHNAME = 'WM'
  CALL IO_Field_read(TPINIFILE,TZFIELD,PWT)
  !
  CALL FIND_FIELD_ID_FROM_MNHNAME('THT',IID,IRESP)
  TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
  TZFIELD%CMNHNAME = 'THM'
  CALL IO_Field_read(TPINIFILE,TZFIELD,PTHT)
  !
  CALL FIND_FIELD_ID_FROM_MNHNAME('PABST',IID,IRESP)
  TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
  TZFIELD%CMNHNAME = 'PABSM'
  CALL IO_Field_read(TPINIFILE,TZFIELD,PPABST)
ELSE
  CALL IO_Field_read(TPINIFILE,'UT',PUT)
  CALL IO_Field_read(TPINIFILE,'VT',PVT)
  CALL IO_Field_read(TPINIFILE,'WT',PWT)
  CALL IO_Field_read(TPINIFILE,'THT',PTHT)
  CALL IO_Field_read(TPINIFILE,'PABST',PPABST)
ENDIF
!
SELECT CASE(HGETTKET)                   
  CASE('READ')
    IF (TPINIFILE%NMNHVERSION(1)<5) THEN
      CALL FIND_FIELD_ID_FROM_MNHNAME('TKET',IID,IRESP)
      TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
      TZFIELD%CMNHNAME = 'TKEM'
      CALL IO_Field_read(TPINIFILE,TZFIELD,PTKET)
    ELSE
      CALL IO_Field_read(TPINIFILE,'TKET',PTKET)
    END IF
    IF ( ( (TPINIFILE%NMNHVERSION(1)==5 .AND. TPINIFILE%NMNHVERSION(2)>0) .OR. TPINIFILE%NMNHVERSION(1)>5 ) &
        .AND. (CCONF == 'RESTA') .AND. LSPLIT_CFL) THEN
      CALL IO_Field_read(TPINIFILE,'TKEMS',PRTKEMS)
    END IF
  CASE('INIT')
    PTKET(:,:,:)   = XTKEMIN
    PRTKEMS(:,:,:) = 0.
END SELECT 
!
SELECT CASE(HGETZWS)
  CASE('READ')
    CALL IO_Field_read(TPINIFILE,'ZWS',PZWS,IRESP)
    !If the field ZWS is not in the file, set its value to XZWS_DEFAULT
    !ZWS is present in files since MesoNH 5.4.2
    IF ( IRESP/=0 ) THEN
      WRITE (YVAL,'( E15.8 )') XZWS_DEFAULT
      CALL PRINT_MSG(NVERB_WARNING,'IO','READ_FIELD','ZWS not found in file: using default value: '//TRIM(YVAL)//' m')
      PZWS(:,:) = XZWS_DEFAULT
    END IF

  CASE('INIT')
    PZWS(:,:)=0.
END SELECT 
!
SELECT CASE(HGETRVT)             ! vapor
  CASE('READ')
    IF (TPINIFILE%NMNHVERSION(1)<5) THEN
      CALL FIND_FIELD_ID_FROM_MNHNAME('RVT',IID,IRESP)
      TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
      TZFIELD%CMNHNAME = 'RVM'
      CALL IO_Field_read(TPINIFILE,TZFIELD,PRT(:,:,:,IDX_RVT))
    ELSE
      CALL IO_Field_read(TPINIFILE,'RVT',PRT(:,:,:,IDX_RVT))
    END IF
  CASE('INIT')
    PRT(:,:,:,IDX_RVT) = 0.
END SELECT 
!
SELECT CASE(HGETRCT)             ! cloud 
  CASE('READ') 
    IF (TPINIFILE%NMNHVERSION(1)<5) THEN
      CALL FIND_FIELD_ID_FROM_MNHNAME('RCT',IID,IRESP)
      TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
      TZFIELD%CMNHNAME = 'RCM'
      CALL IO_Field_read(TPINIFILE,TZFIELD,PRT(:,:,:,IDX_RCT))
    ELSE
      CALL IO_Field_read(TPINIFILE,'RCT',PRT(:,:,:,IDX_RCT))
    END IF
  CASE('INIT')
    PRT(:,:,:,IDX_RCT) = 0.
END SELECT
!
SELECT CASE(HGETRRT)             ! rain 
  CASE('READ') 
    IF (TPINIFILE%NMNHVERSION(1)<5) THEN
      CALL FIND_FIELD_ID_FROM_MNHNAME('RRT',IID,IRESP)
      TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
      TZFIELD%CMNHNAME = 'RRM'
      CALL IO_Field_read(TPINIFILE,TZFIELD,PRT(:,:,:,IDX_RRT))
    ELSE
      CALL IO_Field_read(TPINIFILE,'RRT',PRT(:,:,:,IDX_RRT))
    END IF 
  CASE('INIT')
    PRT(:,:,:,IDX_RRT) = 0.
END SELECT
!
SELECT CASE(HGETRIT)             ! cloud ice
  CASE('READ') 
    IF (TPINIFILE%NMNHVERSION(1)<5) THEN
      CALL FIND_FIELD_ID_FROM_MNHNAME('RIT',IID,IRESP)
      TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
      TZFIELD%CMNHNAME = 'RIM'
      CALL IO_Field_read(TPINIFILE,TZFIELD,PRT(:,:,:,IDX_RIT))
    ELSE
      CALL IO_Field_read(TPINIFILE,'RIT',PRT(:,:,:,IDX_RIT))
    END IF 
  CASE('INIT')
    PRT(:,:,:,IDX_RIT) = 0.
END SELECT
!
SELECT CASE(HGETRST)             ! snow
  CASE('READ')
    IF (TPINIFILE%NMNHVERSION(1)<5) THEN
      CALL FIND_FIELD_ID_FROM_MNHNAME('RST',IID,IRESP)
      TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
      TZFIELD%CMNHNAME = 'RSM'
      CALL IO_Field_read(TPINIFILE,TZFIELD,PRT(:,:,:,IDX_RST))
    ELSE
      CALL IO_Field_read(TPINIFILE,'RST',PRT(:,:,:,IDX_RST))
    END IF 
  CASE('INIT')
    PRT(:,:,:,IDX_RST) = 0.
END SELECT
!
SELECT CASE(HGETRGT)             ! graupel
  CASE('READ') 
    IF (TPINIFILE%NMNHVERSION(1)<5) THEN
      CALL FIND_FIELD_ID_FROM_MNHNAME('RGT',IID,IRESP)
      TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
      TZFIELD%CMNHNAME = 'RGM'
      CALL IO_Field_read(TPINIFILE,TZFIELD,PRT(:,:,:,IDX_RGT))
    ELSE
      CALL IO_Field_read(TPINIFILE,'RGT',PRT(:,:,:,IDX_RGT))
    END IF 
  CASE('INIT')
    PRT(:,:,:,IDX_RGT) = 0.
END SELECT
!
SELECT CASE(HGETRHT)             ! hail
  CASE('READ') 
    IF (TPINIFILE%NMNHVERSION(1)<5) THEN
      CALL FIND_FIELD_ID_FROM_MNHNAME('RHT',IID,IRESP)
      TZFIELD = TFIELDMETADATA( TFIELDLIST(IID) )
      TZFIELD%CMNHNAME = 'RHM'
      CALL IO_Field_read(TPINIFILE,TZFIELD,PRT(:,:,:,IDX_RHT))
    ELSE
      CALL IO_Field_read(TPINIFILE,'RHT',PRT(:,:,:,IDX_RHT))
    END IF 
  CASE('INIT')
    PRT(:,:,:,IDX_RHT) = 0.
END SELECT
!
SELECT CASE(HGETCIT)             ! ice concentration
  CASE('READ')
    IF (SIZE(PCIT) /= 0 ) CALL IO_Field_read(TPINIFILE,'CIT',PCIT)
  CASE('INIT')
    PCIT(:,:,:)=0.
END SELECT
!
IF (LIBM .AND. CPROGRAM=='MESONH') THEN
   !
  TZFIELD = TFIELDMETADATA( &
    CMNHNAME   = 'LSFP',    &
    CLONGNAME  = 'LSFP',    &
    CSTDNAME   = '',        &
    CUNITS     = 'm',       &
    CDIR       = 'XY',      &
    NGRID      = 1,         &
    NTYPE      = TYPEREAL,  &
    NDIMS      = 3,         &
    LTIMEDEP   = .TRUE.     )
   !
   CALL IO_Field_read(TPINIFILE,TZFIELD,PIBM_LSF)
   !
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
   CALL IO_Field_read(TPINIFILE,TZFIELD,PIBM_XMUT)
   !
ENDIF
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
CALL IO_Field_read(TPINIFILE,TZFIELD,ZLRECYCL,IRESP)
!If field not found (file from older version of MesoNH) => set ZLRECYCL to false
IF ( IRESP /= 0 ) ZLRECYCL = .FALSE.

IF (ZLRECYCL) THEN
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
  CALL IO_Field_read(TPINIFILE,TZFIELD,NR_COUNT)
  !
  IF (NR_COUNT .NE. 0) THEN
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PUMEANW)
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PVMEANW)
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PWMEANW)
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PUMEANN)
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PVMEANN)
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PWMEANN)
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PUMEANE)
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PVMEANE)
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PWMEANE)
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PUMEANS)
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PVMEANS)
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
      CALL IO_Field_read(TPINIFILE,TZFIELD,PWMEANS)
    ENDIF
  ENDIF  
ENDIF

! Blaze fire model
IF (LBLAZE .AND. CCONF=='RESTA') THEN
  ! Blaze is not compliant with MNHVERSION(1)<5
  ! Blaze begins with MNH 5.3.1
  CALL IO_Field_read(TPINIFILE,'FMPHI',PLSPHI,IRESP)
  IF (IRESP /= 0) PLSPHI(:,:,:) = 0.
  CALL IO_Field_read(TPINIFILE,'FMBMAP',PBMAP,IRESP)
  IF (IRESP /= 0) PBMAP(:,:,:) = -1.
  CALL IO_Field_read(TPINIFILE,'FMASE',PFMASE,IRESP)
  IF(IRESP == 0) THEN
    ! flag for the use of restart value for ASE initialization
    LRESTA_ASE = .TRUE.
  ELSE
    CALL PRINT_MSG( NVERB_WARNING, 'IO', 'READ_FIELD', 'PFMASE set to 0' )
    PFMASE(:,:,:) = 0.
  END IF
  CALL IO_Field_read(TPINIFILE,'FMAWC',PFMAWC,IRESP)
  ! flag for the use of restart value for AWC initialization
  IF(IRESP == 0) THEN
    LRESTA_AWC = .TRUE.
  ELSE
    CALL PRINT_MSG( NVERB_WARNING, 'IO', 'READ_FIELD', 'PFMAWC set to 0' )
    PFMAWC(:,:,:) = 0.
  END IF
  ! read wind on fire grid if present
  IF (LWINDFILTER) THEN
    ! read in file only if wind filtering is required
    SELECT CASE(CWINDFILTER)
    CASE('EWAM')
      ! read u
      CALL IO_Field_read(TPINIFILE,'FMWINDU',PFMWINDU,IRESP)
      ! flag for EWAM filtered u wind
      IF(IRESP == 0) THEN
        LRESTA_EWAM = .TRUE.
      ELSE
        CALL PRINT_MSG( NVERB_WARNING, 'IO', 'READ_FIELD', 'PFMWINDU set to 0' )
        PFMWINDU(:,:,:) = 0.
      END IF
      ! read v
      CALL IO_Field_read(TPINIFILE,'FMWINDV',PFMWINDV,IRESP)
      ! flag for EWAM filtered v wind
      IF(IRESP == 0 .AND. LRESTA_EWAM) THEN
        ! u and v fields found
        LRESTA_EWAM = .TRUE.
      ELSE
        ! u or v fields NOT found
        LRESTA_EWAM = .FALSE.
      END IF
      IF (IRESP /= 0) THEN
        CALL PRINT_MSG( NVERB_WARNING, 'IO', 'READ_FIELD', 'PFMWINDV set to 0' )
        PFMWINDV(:,:,:) = 0.
      END IF
      ! read w
      CALL IO_Field_read(TPINIFILE,'FMWINDW',PFMWINDW,IRESP)
      ! flag for EWAM filtered w wind
      IF(IRESP == 0 .AND. LRESTA_EWAM) THEN
        ! u and v and w fields found
        LRESTA_EWAM = .TRUE.
      ELSE
        ! u or v or w fields NOT found
        LRESTA_EWAM = .FALSE.
      END IF
      IF (IRESP /= 0) THEN
        CALL PRINT_MSG( NVERB_WARNING, 'IO', 'READ_FIELD', 'PFMWINDW set to 0' )
        PFMWINDW(:,:,:) = 0.
      END IF

    CASE('WLIM')
      CALL IO_Field_read(TPINIFILE,'FMHWS',PFMHWS,IRESP)
      ! flag for WLIM filtered horizontal wind speed
      IF(IRESP == 0) THEN
        LRESTA_WLIM = .TRUE.
      ELSE
        CALL PRINT_MSG( NVERB_WARNING, 'IO', 'READ_FIELD', 'PFMHWS set to 0' )
        PFMHWS(:,:,:) = 0.
      END IF
    END SELECT
  END IF
END IF
!
!  Scalar Variables Reading : Users, C2R2, C1R3, LIMA, ELEC, Chemical SV
!
ISV= SIZE(PSVT,4)
!
DO JSV = 1, NSV              ! initialize according to the get indicators
  SELECT CASE( HGETSVT(JSV) )
    CASE ('READ')
      TZFIELD = TSVLIST(JSV)

      IF ( GOLDFILEFORMAT ) THEN
        IF ( ( JSV >= 1         .AND. JSV <= NSV_USER  ) .OR. &
             ( JSV >= NSV_PPBEG .AND. JSV <= NSV_PPEND ) .OR. &
#ifdef MNH_FOREFIRE
             ( JSV >= NSV_FFBEG .AND. JSV <= NSV_FFEND ) .OR. &
#endif
             ( JSV >= NSV_CSBEG .AND. JSV <= NSV_CSEND )      ) THEN
          !Some variables were written with an other name in MesoNH < 5.6
          WRITE(TZFIELD%CMNHNAME,'(A3,I3.3)')'SVT',JSV
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CSTDNAME   = ''
          TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
        ELSE
          !Scalar variables were written with a T suffix in older versions
          TZFIELD%CMNHNAME  = TRIM( TZFIELD%CMNHNAME )  // 'T'
          TZFIELD%CLONGNAME = TRIM( TZFIELD%CLONGNAME ) // 'T'
        END IF
      END IF

      CALL IO_Field_read( TPINIFILE, TZFIELD, PSVT(:,:,:,JSV), IRESP )

      IF ( IRESP /= 0 ) THEN
        CALL PRINT_MSG( NVERB_WARNING, 'IO', 'READ_FIELD', 'PSVT set to 0 for ' // TRIM( TZFIELD%CMNHNAME ) )
        PSVT(:,:,:,JSV) = 0.
      END IF

    CASE ('INIT')
      PSVT(:,:,:,JSV) = 0.

      IF ( JSV == NSV_C2R2END ) THEN
        IF ( LSUPSAT .AND. (HGETRVT == 'READ') ) THEN
          ZWORK(:,:,:) = (PPABST(:,:,:)/XP00 )**(XRD/XCPD)
          ZWORK(:,:,:) = PTHT(:,:,:)*ZWORK(:,:,:)
          ZWORK(:,:,:) = EXP(XALPW-XBETAW/ZWORK(:,:,:)-XGAMW*LOG(ZWORK(:,:,:)))
          !rvsat
          ZWORK(:,:,:) = (XMV / XMD)*ZWORK(:,:,:)/(PPABST(:,:,:)-ZWORK(:,:,:))
          ZWORK(:,:,:) = PRT(:,:,:,IDX_RVT)/ZWORK(:,:,:)
          PSVT(:,:,:,NSV_C2R2END ) = ZWORK(:,:,:)
        END IF
      END IF

  END SELECT
END DO

DO JSV = NSV_PPBEG, NSV_PPEND
  SELECT CASE( HGETSVT(JSV) )
    CASE ('READ')
      WRITE( YNUM3, '( I3.3 )' ) JSV

      TZFIELD = TFIELDMETADATA(            &
        CMNHNAME   = 'ATC' // YNUM3,       &
        CSTDNAME   = '',                   &
        CLONGNAME  = 'ATC' // YNUM3,       &
        CCOMMENT   = 'X_Y_Z_ATC' // YNUM3, &
        CUNITS     = 'm-3',                &
        CDIR       = 'XY',                 &
        NGRID      = 1,                    &
        NTYPE      = TYPEREAL,             &
        NDIMS      = 3,                    &
        LTIMEDEP   = .TRUE.                )

      CALL IO_Field_read( TPINIFILE, TZFIELD, PATC(:,:,:,JSV-NSV_PPBEG+1), IRESP )

      IF ( IRESP /= 0 ) THEN
        PATC(:,:,:,JSV-NSV_PPBEG+1) = 0.
      ENDIF

    CASE ('INIT')
      PATC(:,:,:,JSV-NSV_PPBEG+1) = 0.

  END SELECT
END DO

IF ( NSV_SNW >= 1 ) THEN
  TZFIELD = TFIELDMETADATA(                &
    CMNHNAME   = 'generic for SNOWCANO_M', &
    CUNITS     = 'kg kg-1',                &
    CDIR       = 'XY',                     &
    NGRID      = 1,                        &
    NTYPE      = TYPEREAL,                 &
    NDIMS      = 2,                        &
    LTIMEDEP   = .TRUE.                    )
  DO JSV = 1, NSV_SNW
    SELECT CASE(HGETSVT(JSV))
      CASE ('READ')
        WRITE(TZFIELD%CMNHNAME,'(A10,I3.3)')'SNOWCANO_M',JSV
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A8,I3.3)') 'X_Y_Z_','SNOWCANO',JSV
        CALL IO_Field_read( TPINIFILE, TZFIELD, XSNWCANO(:,:,JSV) )
      CASE ('INIT')
        XSNWCANO(:,:,JSV) = 0.
    END SELECT
  END DO
END IF
!
IF (CCONF == 'RESTA') THEN
  IF (CTEMP_SCHEME/='LEFR') THEN
    CALL IO_Field_read(TPINIFILE,'US_PRES',PRUS_PRES)
    CALL IO_Field_read(TPINIFILE,'VS_PRES',PRVS_PRES)
    CALL IO_Field_read(TPINIFILE,'WS_PRES',PRWS_PRES)
  END IF
  IF (LSPLIT_CFL) THEN
    CALL IO_Field_read(TPINIFILE,'THS_CLD',PRTHS_CLD)
    DO JRR = 1, SIZE(PRT,4)
      SELECT CASE(JRR)
        CASE (1)
          CALL IO_Field_read(TPINIFILE,'RVS_CLD',PRRS_CLD(:,:,:,JRR))
        CASE (2)
          CALL IO_Field_read(TPINIFILE,'RCS_CLD',PRRS_CLD(:,:,:,JRR))
        CASE (3)
          CALL IO_Field_read(TPINIFILE,'RRS_CLD',PRRS_CLD(:,:,:,JRR))
        CASE (4)
          CALL IO_Field_read(TPINIFILE,'RIS_CLD',PRRS_CLD(:,:,:,JRR))
        CASE (5)
          CALL IO_Field_read(TPINIFILE,'RSS_CLD',PRRS_CLD(:,:,:,JRR))
        CASE (6)
          CALL IO_Field_read(TPINIFILE,'RGS_CLD',PRRS_CLD(:,:,:,JRR))
        CASE (7)
          CALL IO_Field_read(TPINIFILE,'RHS_CLD',PRRS_CLD(:,:,:,JRR))
        CASE DEFAULT
          CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_FIELD','PRT is too big')
      END SELECT
    END DO
    DO JSV = NSV_C2R2BEG,NSV_C2R2END
      IF (JSV == NSV_C2R2BEG ) THEN
        TZFIELD = TFIELDMETADATA(       &
          CMNHNAME   = 'RSVS_CLD1',     &
          CSTDNAME   = '',              &
          CLONGNAME  = 'RSVS_CLD1',     &
          CUNITS     = '1',             &
          CDIR       = 'XY',            &
          CCOMMENT   = 'X_Y_Z_RHS_CLD', &
          NGRID      = 1,               &
          NTYPE      = TYPEREAL,        &
          NDIMS      = 3,               &
          LTIMEDEP   = .TRUE.           )
        CALL IO_Field_read(TPINIFILE,TZFIELD,PRSVS_CLD(:,:,:,JSV))
      END IF
      IF (JSV == NSV_C2R2BEG ) THEN
        TZFIELD = TFIELDMETADATA(       &
          CMNHNAME   = 'RSVS_CLD2',     &
          CSTDNAME   = '',              &
          CLONGNAME  = 'RSVS_CLD2',     &
          CUNITS     = '1',             &
          CDIR       = 'XY',            &
          CCOMMENT   = 'X_Y_Z_RHS_CLD', &
          NGRID      = 1,               &
          NTYPE      = TYPEREAL,        &
          NDIMS      = 3,               &
          LTIMEDEP   = .TRUE.           )
        CALL IO_Field_read(TPINIFILE,TZFIELD,PRSVS_CLD(:,:,:,JSV))
      END IF
    END DO
  END IF
END IF
!
!*       2.1  Time t-dt:
!
IF (CPROGRAM=='MESONH' .AND. HUVW_ADV_SCHEME(1:3)=='CEN' .AND. &
        HTEMP_SCHEME == 'LEFR' ) THEN
  IF (CCONF=='RESTA') THEN
    CALL IO_Field_read(TPINIFILE,'UM', PUM)
    CALL IO_Field_read(TPINIFILE,'VM', PVM)
    CALL IO_Field_read(TPINIFILE,'WM', PWM)
    CALL IO_Field_read(TPINIFILE,'DUM',PDUM)
    CALL IO_Field_read(TPINIFILE,'DVM',PDVM)
    CALL IO_Field_read(TPINIFILE,'DWM',PDWM)
  ELSE
    PUM = PUT
    PVM = PVT
    PWM = PWT
  END IF
END IF
!
!*       2.2a  3D LS fields  
!
!
CALL INI_LS(TPINIFILE,HGETRVT,GLSOURCE,PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM)
!
!
!*       2.2b  2D "surfacic" LB fields   
!
!
CALL INI_LB(TPINIFILE,GLSOURCE,ISV,                                   &
     KSIZELBX_ll,KSIZELBXU_ll,KSIZELBY_ll,KSIZELBYV_ll,               &
     KSIZELBXTKE_ll,KSIZELBYTKE_ll,                                   &
     KSIZELBXR_ll,KSIZELBYR_ll,KSIZELBXSV_ll,KSIZELBYSV_ll,           &
     HGETTKET,HGETRVT,HGETRCT,HGETRRT,HGETRIT,HGETRST,                &
     HGETRGT,HGETRHT,HGETSVT,                                         &
     PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,            &
     PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM             )
!
!
!*       2.3  Some special variables:
!
CALL IO_Field_read(TPINIFILE,'DRYMASST',PDRYMASST) ! dry mass
IF (CCONF=='RESTA') THEN
  CALL IO_Field_read(TPINIFILE,'DRYMASSS',PDRYMASSS,IRESP) ! dry mass tendency

  ! DRYMASSS was not written in backup files before MesoNH 5.5.1
  IF ( IRESP /= 0 ) THEN
    CALL PRINT_MSG( NVERB_WARNING, 'IO', 'READ_FIELD', 'PDRYMASSS set to 0 for ' // TRIM( TZFIELD%CMNHNAME ) )
    PDRYMASSS = 0.
  END IF
ELSE
  PDRYMASSS=XUNDEF   ! should not be used
END IF
!
SELECT CASE(HGETSRCT)                ! turbulent flux SRC at time t
  CASE('READ')
    CALL IO_Field_read(TPINIFILE,'SRCT',PSRCT)
  CASE('INIT')
    PSRCT(:,:,:)=0.
END SELECT
!
SELECT CASE(HGETSIGS)                ! subgrid condensation
  CASE('READ')
    CALL IO_Field_read(TPINIFILE,'SIGS',PSIGS)
  CASE('INIT')
    PSIGS(:,:,:)=0.
END SELECT
!
SELECT CASE(HGETPHC)             ! pH in cloud water
  CASE('READ')
    CALL IO_Field_read(TPINIFILE,'PHC',PPHC)
  CASE('INIT')
    PPHC(:,:,:)=0.
END SELECT
!
SELECT CASE(HGETPHR)             ! pH in rainwater
  CASE('READ')
    CALL IO_Field_read(TPINIFILE,'PHR',PPHR)
  CASE('INIT')
    PPHR(:,:,:)=0.
END SELECT
!
IRESP=0
IF(HGETCLDFR=='READ') THEN           ! cloud fraction
  CALL IO_Field_read(TPINIFILE,'CLDFR',PCLDFR,IRESP)
ENDIF
IF(HGETCLDFR=='INIT' .OR. IRESP /= 0) THEN
  IF(SIZE(PRT,4) > 3) THEN
    WHERE(PRT(:,:,:,2)+PRT(:,:,:,4) > 1.E-30)
      PCLDFR(:,:,:) = 1.
    ELSEWHERE
      PCLDFR(:,:,:) = 0.
    ENDWHERE
  ELSE
    WHERE(PRT(:,:,:,2) > 1.E-30)
      PCLDFR(:,:,:) = 1.
    ELSEWHERE
      PCLDFR(:,:,:) = 0.
    ENDWHERE
  ENDIF
ENDIF
!
IRESP=0
IF(HGETICEFR=='READ') THEN           ! cloud fraction
  CALL IO_Field_read(TPINIFILE,'ICEFR',PICEFR,IRESP)
ENDIF
IF(HGETCLDFR=='INIT' .OR. IRESP /= 0) THEN
  IF(SIZE(PRT,4) > 3) THEN
    WHERE(PRT(:,:,:,4) > 1.E-30)
       PICEFR(:,:,:) = 1.
    ELSEWHERE
      PICEFR(:,:,:) = 0.
    ENDWHERE
  ELSE
     PICEFR(:,:,:) = 0.
  ENDIF
ENDIF
!
!* boundary layer depth
!
IF (HGETBL_DEPTH=='READ') THEN
  CALL IO_Field_read(TPINIFILE,'BL_DEPTH',PBL_DEPTH)
ELSE
  PBL_DEPTH(:,:)=XUNDEF
END IF
!
!* surface boundary layer depth
!
IF (HGETSBL_DEPTH=='READ') THEN
  CALL IO_Field_read(TPINIFILE,'SBL_DEPTH',PSBL_DEPTH)
ELSE
  PSBL_DEPTH(:,:)=0.
END IF
!
!* Contribution from MAss Flux parameterizations to vert. flux of buoyancy
!
SELECT CASE(HGETTKET)                   
  CASE('READ') 
    IF (CSCONV=='EDKF') THEN 
      CALL IO_Field_read(TPINIFILE,'WTHVMF',PWTHVMF)
    ELSE
      PWTHVMF(:,:,:)=0
    ENDIF
  CASE('INIT')
    PWTHVMF(:,:,:)=0.
END SELECT 
!-------------------------------------------------------------------------------
!
!*       2.4   READ FORCING VARIABLES
!              ----------------------
!
! READ FIELD ONLY FOR MODEL1 (identical for all model in GN)
IF (LOCEAN .AND. (.NOT.LCOUPLES) .AND. (KOCEMI==1)) THEN
!
  CALL IO_Field_read(TPINIFILE,'NFRCLT',NFRCLT)
  CALL IO_Field_read(TPINIFILE,'NINFRT',NINFRT)
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
  ALLOCATE(XSSUFL_T(NFRCLT))
  CALL IO_Field_read(TPINIFILE,TZFIELD,XSSUFL_T(:))
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
  ALLOCATE(XSSVFL_T(NFRCLT))
  CALL IO_Field_read(TPINIFILE,TZFIELD,XSSVFL_T(:))
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
  ALLOCATE(XSSTFL_T(NFRCLT))
  CALL IO_Field_read(TPINIFILE,TZFIELD,XSSTFL_T(:))
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
  ALLOCATE(XSSOLA_T(NFRCLT))
  CALL IO_Field_read(TPINIFILE,TZFIELD,XSSOLA_T(:))
!
END IF ! ocean sfc forcing end    

!
IF ( LFORCING ) THEN
  DO JT=1,KFRC
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,TPDTFRC(JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PUFRC(:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PVFRC(:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PWFRC(:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PTHFRC(:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PRVFRC(:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PTENDTHFRC(:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PTENDRVFRC(:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PGXTHFRC(:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PGYTHFRC(:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PPGROUNDFRC(JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PTENDUFRC(:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PTENDVFRC(:,JT))
  END DO
END IF
!
!-------------------------------------------------------------------------------
IF (L2D_ADV_FRC) THEN

  DO JT=1,KADVFRC  
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,TPDTADVFRC(JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PDTHFRC(:,:,:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PDRVFRC(:,:,:,JT))
  ENDDO
ENDIF
!
IF (L2D_REL_FRC) THEN

  DO JT=1,KRELFRC  
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,TPDTRELFRC(JT))
    !
    ! Relaxation
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PTHREL(:,:,:,JT))
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
    CALL IO_Field_read(TPINIFILE,TZFIELD,PRVREL(:,:,:,JT))
  ENDDO
ENDIF
!
IF (LUV_FLX) THEN
  IF ( CCONF /= 'START' .OR. CPROGRAM=='SPAWN ' ) THEN
    CALL IO_Field_read(TPINIFILE,'VU_FLX',PVU_FLUX_M)
  ELSE IF (CCONF == 'START') THEN
    PVU_FLUX_M(:,:,:)=0.
  END IF
ENDIF
!
IF (LTH_FLX) THEN
  IF ( CCONF /= 'START' .OR. CPROGRAM=='SPAWN ' ) THEN
    CALL IO_Field_read(TPINIFILE,'VT_FLX',PVTH_FLUX_M)
    CALL IO_Field_read(TPINIFILE,'WT_FLX',PWTH_FLUX_M)
   ELSE IF (CCONF == 'START') THEN
       PWTH_FLUX_M(:,:,:)=0.
       PVTH_FLUX_M(:,:,:)=0.
   END IF
ENDIF
!
!-------------------------------------------------------------------------------
!
!
!*       3.    PRINT ON OUTPUT-LISTING
!              ----------------------
!
IF (NVERB >= 10 .AND. .NOT. L1D) THEN
  IIUP = SIZE(PUT,1)
  IJUP = SIZE(PVT,2) 
  ILUOUT= TLUOUT%NLU
! 
  WRITE(ILUOUT,FMT=*) 'READ_FIELD: Some PUT values:'
  WRITE(ILUOUT,FMT=*) '(1,1,JK)   (IIU/2,IJU/2,JK)   (IIU,IJU,JK)    JK  '
  DO JKLOOP=1,KKU
    WRITE(ILUOUT,FMT=*) PUT(1,1,JKLOOP),PUT(IIUP/2,IJUP/2,JKLOOP), &
    PUT(IIUP,KJU,JKLOOP),JKLOOP    
  END DO
!
  WRITE(ILUOUT,FMT=*) 'READ_FIELD: Some PVT values:'
  WRITE(ILUOUT,FMT=*) '(1,1,JK)   (IIU/2,IJU/2,JK)   (IIU,IJU,JK)    JK  '
  DO JKLOOP=1,KKU
    WRITE(ILUOUT,FMT=*) PVT(1,1,JKLOOP),PVT(IIUP/2,IJUP/2,JKLOOP), &
    PVT(IIUP,IJUP,JKLOOP),JKLOOP    
  END DO
!
  WRITE(ILUOUT,FMT=*) 'READ_FIELD: Some PWT values:'
  WRITE(ILUOUT,FMT=*) '(1,1,JK)   (IIU/2,IJU/2,JK)   (IIU,IJU,JK)    JK  '
  DO JKLOOP=1,KKU
    WRITE(ILUOUT,FMT=*) PWT(1,1,JKLOOP),PWT(IIUP/2,IJUP/2,JKLOOP), &
    PWT(IIUP,IJUP,JKLOOP),JKLOOP    
  END DO
!
  WRITE(ILUOUT,FMT=*) 'READ_FIELD: Some PTHT values:'
  WRITE(ILUOUT,FMT=*) '(1,1,JK)   (IIU/2,IJU/2,JK)   (IIU,IJU,JK)    JK  '
  DO JKLOOP=1,KKU
    WRITE(ILUOUT,FMT=*) PTHT(1,1,JKLOOP),PTHT(IIUP/2,IJUP/2,JKLOOP), &
    PTHT(IIUP,IJUP,JKLOOP),JKLOOP    
  END DO
!
  IF(SIZE(PTKET,1) /=0) THEN
    WRITE(ILUOUT,FMT=*) 'READ_FIELD: Some PTKET values:'
    WRITE(ILUOUT,FMT=*) '(1,1,JK)   (IIU/2,IJU/2,JK)   (IIU,IJU,JK)    JK  '
    DO JKLOOP=1,KKU
      WRITE(ILUOUT,FMT=*) PTKET(1,1,JKLOOP),PTKET(IIUP/2,IJUP/2,JKLOOP), &
      PTKET(IIUP,IJUP,JKLOOP),JKLOOP    
    END DO
  END IF
!
  IF (SIZE(PRT,4) /= 0) THEN
    WRITE(ILUOUT,FMT=*) 'READ_FIELD: Some PRT values:'
    DO JRR = 1, SIZE(PRT,4)
      WRITE(ILUOUT,FMT=*) 'JRR = ',JRR
      WRITE(ILUOUT,FMT=*) '(1,1,JK)   (IIU/2,IJU/2,JK)   (IIU,IJU,JK)    JK  '
      DO JKLOOP=1,KKU
        WRITE(ILUOUT,FMT=*) PRT(1,1,JKLOOP,JRR),PRT(IIUP/2,IJUP/2,JKLOOP,JRR), &
        PRT(IIUP,IJUP,JKLOOP,JRR),JKLOOP    
      END DO
    END DO
!
  END IF   
!
  IF (SIZE(PSVT,4) /= 0) THEN
    WRITE(ILUOUT,FMT=*) 'READ_FIELD: Some PSVT values:'
    DO JRR = 1, SIZE(PSVT,4)
      WRITE(ILUOUT,FMT=*) 'JRR = ',JRR
      WRITE(ILUOUT,FMT=*) '(1,1,JK)   (IIU/2,IJU/2,JK)   (IIU,IJU,JK)    JK  '
      DO JKLOOP=1,KKU
        WRITE(ILUOUT,FMT=*) PSVT(1,1,JKLOOP,JRR),PSVT(IIUP/2,IJUP/2,JKLOOP,JRR), &
        PSVT(IIUP,IJUP,JKLOOP,JRR),JKLOOP    
      END DO
    END DO
!
  END IF   
END IF 
!-------------------------------------------------------------------------------
! 
!
END SUBROUTINE READ_FIELD
