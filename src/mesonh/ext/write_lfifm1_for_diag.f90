!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!################################
MODULE MODI_WRITE_LFIFM1_FOR_DIAG
!################################
INTERFACE
      SUBROUTINE WRITE_LFIFM1_FOR_DIAG(TPFILE,HDADFILE)
!
USE MODD_IO, ONLY: TFILEDATA
!
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE       ! outpput data file
CHARACTER(LEN=28), INTENT(IN) :: HDADFILE     ! corresponding FM-file name of 
                                              ! its DAD model
!
END SUBROUTINE WRITE_LFIFM1_FOR_DIAG
END INTERFACE
END MODULE MODI_WRITE_LFIFM1_FOR_DIAG
!
!     ##################################################
      SUBROUTINE WRITE_LFIFM1_FOR_DIAG(TPFILE,HDADFILE)
!     ##################################################
!
!!****  *WRITE_LFIFM1* - routine to write a LFIFM file for model 1
!!
!!    PURPOSE
!!    -------
!        The purpose of this routine is to write an initial LFIFM File 
!     of name YFMFILE2//'.lfi' with the FM routines.  
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
!!      FMWRIT : FM-routine to write a record
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_DIM1   : contains dimensions
!!      Module MODD_TIME1   : contains time variables and uses MODD_TIME
!!      Module MODD_GRID    : contains spatial grid variables for all models
!!      Module MODD_GRID1 : contains spatial grid variables
!!      Module MODD_REF     : contains reference state variables
!!      Module MODD_LUNIT1: contains logical unit variables.
!!      Module MODD_CONF    : contains configuration variables for all models
!!      Module MODD_CONF1  : contains configuration variables
!!      Module MODD_FIELD1  : contains prognostic variables
!!      Module MODD_GR_FIELD1 : contains surface prognostic variables
!!      Module MODD_LSFIELD1  : contains Larger Scale variables
!!      Module MODD_PARAM1    : contains parameterization options
!!      Module MODD_TURB1    : contains turbulence options
!!      Module MODD_FRC    : contains forcing variables
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!  	V. Ducrocq   *Meteo France* 
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
!!       J.Stein       23/01/95 add a TKE switch and MODD_PARAM1
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
!!       J.P. Lafore          30/07/96 add YFMFILE2 and HDADFILE writing
!!                                     corresponding to MY_NAME and DAD_NAME (for nesting)
!!       V.Masson             08/10/96 add LTHINSHELL
!!       J.-P. Pinty   15/12/96 add the microphysics (ice)
!!       J.-P. Pinty   11/01/97 add the deep convection
!!       J.-P. Pinty   27/01/97 split the recording of the SV array
!!       J.-P. Pinty   29/01/97 set recording of PRCONV and PACCONV in mm/h and
!!                                                         mm respectively
!!       J. Viviand    04/02/97 convert precipitation rates in mm/h
!!       P. Hereil     04/12/97 add the calculation of cloud top and moist PV
!!       P.Hereil N Asencio 3/02/98 add the calculation of  precipitation on large scale grid mesh
!!       N Asencio 2/10/98 suppress flux calculation if start file
!!       V Masson 25/11/98 places dummy arguments in module MODD_DIAG_FLAG
!!       V Masson 04/01/00 removes TSZ0 option
!!       J.-P. Pinty   29/11/02 add C3R5, ICE2, ICE4, CELEC
!!       V Masson 01/2004  removes surface (externalization)
!!       P. Tulet 01/2005   add dust, orilam
!!       M. Leriche 04/2007 add aqueous concentration in M
!!       O. Caumont 03/2008 add simulation of radar observations
!!       O. Caumont 14/09/2009 modifications to allow for polar outputs (radar diagnostics)
!!       October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!       G. Tanguy  10/2009 add possibility to run radar after 
!!                          PREP_REAL_CASE with AROME
!!       O. Caumont 01/2011 [radar diagnostics] add control check for NMAX; revise comments
!!       O. Caumont 05/2011 [radar diagnostics] change output format
!!       G.Tanguy/ JP Pinty/ JP Chabureau 18/05/2011 : add lidar simulator
!!       S.Bielli 12/2012 : add latitude and longitude
!!       F. Duffourg 02/2013 : add new fields
!!      J.Escobar 21/03/2013: for HALOK get correctly local array dim/bound
!!       J. escobar 27/03/2014 : write LAT/LON only in not CARTESIAN case
!!       G.Delautier    2014 : remplace MODD_RAIN_C2R2_PARAM par MODD_RAIN_C2R2_KHKO_PARAM
!!       C. Augros 2014 : new radar simulator (T matrice)
!!       D.Ricard 2015 : add THETAES + POVOES  (LMOIST_ES=T)
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!       C.Lac  04/2016 : add visibility and droplet deposition
!! 10/2017      (G.Delautier) New boundary layer height : replace LBLTOP by CBLTOP 
!!       T.Dauhut      10/2017 : add parallel 3D clustering
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!       D.Ricard and P.Marquet 2016-2017 : THETAL + THETAS1 POVOS1 or THETAS2 POVOS2
!!                                        if  LMOIST_L     LMOIST_S1   or  LMOIST_S2
!  P. Wautelet 08/02/2019: minor bug: compute ZWORK36 only when needed
!  S  Bielli      02/2019: sea salt: significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 18/03/2020: remove ICE2 option
!      B. Vie          06/2020 Add prognostic supersaturation for LIMA
!  P. Wautelet 11/03/2021: bugfix: correct name for NSV_LIMA_IMM_NUCL
!  J.L Redelsperger 03/2021 Adding OCEAN LES Case and Autocoupled O-A LES 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIM_n
USE MODD_CONF
USE MODD_CONF_n
use modd_field,           only: tfielddata, tfieldlist, TYPEINT, TYPEREAL
USE MODD_GRID
USE MODD_GRID_n
USE MODD_IO, ONLY : TFILEDATA
USE MODD_METRICS_n
USE MODD_TIME
USE MODD_TIME_n
USE MODD_DYN_n
USE MODD_FIELD_n
USE MODD_GR_FIELD_n
USE MODD_LSFIELD_n
USE MODD_PARAM_n
USE MODD_CURVCOR_n
USE MODD_REF
USE MODD_REF_n
USE MODD_LUNIT, ONLY : TLUOUT0
USE MODD_LUNIT_n
USE MODD_TURB_n
USE MODD_RADIATIONS_n
USE MODD_FRC
USE MODD_PRECIP_n
USE MODD_CST
USE MODD_CLOUDPAR
USE MODD_DEEP_CONVECTION_n
USE MODD_PARAM_KAFR_n
USE MODD_NESTING
USE MODD_PARAMETERS
USE MODD_DIAG_FLAG
USE MODD_NSV
USE MODD_CH_M9_n,         ONLY : CNAMES, NEQAQ
USE MODD_RAIN_C2R2_DESCR, ONLY : C2R2NAMES
USE MODD_ICE_C1R3_DESCR,  ONLY : C1R3NAMES
USE MODD_ELEC_DESCR, ONLY : CELECNAMES
USE MODD_RAIN_C2R2_KHKO_PARAM
USE MODD_ICE_C1R3_PARAM
USE MODD_PARAM_ICE,       ONLY : LSEDIC
USE MODD_PARAM_LIMA,      ONLY : NMOD_CCN, NMOD_IFN, NMOD_IMM, NINDICE_CCN_IMM,&
                                 LSCAV, LHHONI, LAERO_MASS,                    &
                                 LLIMA_DIAG,                                   &
                                 NSPECIE, XMDIAM_IFN, XSIGMA_IFN, ZFRAC=>XFRAC,&
                                 XR_MEAN_CCN, XLOGSIG_CCN, NMOM_S, NMOM_G, NMOM_H 
USE MODD_PARAM_LIMA_WARM, ONLY : CLIMA_WARM_CONC, CAERO_MASS
USE MODD_PARAM_LIMA_COLD, ONLY : CLIMA_COLD_CONC
USE MODD_LG,              ONLY : CLGNAMES
USE MODD_PASPOL,          ONLY : LPASPOL
USE MODD_CONDSAMP,        ONLY : LCONDSAMP
!
USE MODD_DIAG_FLAG
USE MODD_RADAR, ONLY: XLAT_RAD,XELEV,&
     XSTEP_RAD,NBRAD,NBELEV,NBAZIM,NBSTEPMAX,&
     NCURV_INTERPOL,LATT,LCART_RAD,NPTS_H,NPTS_V,XGRID,&
     LREFR,LDNDZ,NMAX,CNAME_RAD,NDIFF,&
     XLON_RAD,XALT_RAD,XLAM_RAD,XDT_RAD,LWBSCS,LWREFL
use modd_precision,       only: MNHREAL_MPI
!
USE MODI_RADAR_SIMULATOR
!
USE MODD_DUST
USE MODD_CSTS_DUST
USE MODD_SALT
USE MODD_BLOWSNOW
USE MODD_CH_AEROSOL
USE MODD_CH_AERO_n
USE MODD_CH_MNHC_n
USE MODE_DUST_PSD
USE MODE_SALT_PSD
USE MODE_BLOWSNOW_PSD
USE MODE_AERO_PSD
USE MODI_GRADIENT_M
USE MODI_GRADIENT_W
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_SHUMAN
USE MODI_RADAR_RAIN_ICE
USE MODI_INI_RADAR
USE MODI_COMPUTE_MEAN_PRECIP
USE MODI_UV_TO_ZONAL_AND_MERID
USE MODI_CALCSOUND
USE MODI_FREE_ATM_PROFILE
USE MODI_GPS_ZENITH
USE MODI_CONTRAV
!
use mode_field,            only: Find_field_id_from_mnhname
USE MODE_GRIDPROJ
USE MODE_GATHER_ll
USE MODE_IO_FIELD_WRITE,   only: IO_Field_write
USE MODE_IO_FILE,          only: IO_File_close, IO_File_open
USE MODE_IO_MANAGE_STRUCT, only: IO_File_add2list
USE MODE_ll
use mode_msg
USE MODE_THERMO
USE MODE_TOOLS,            ONLY: UPCASE
USE MODE_MODELN_HANDLER
USE MODI_LIDAR
USE MODI_CLUSTERING
!
USE MODD_MPIF
USE MODD_VAR_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE       ! outpput data file
CHARACTER(LEN=28), INTENT(IN) :: HDADFILE     ! corresponding FM-file name of 
                                              ! its DAD model
!
!*       0.2   Declarations of local variables
!
INTEGER           :: IRESP          ! return-code for the file routines 
!
CHARACTER(LEN=3)  :: YFRC           ! to mark the time of the forcing
CHARACTER(LEN=31) :: YFGRI          ! file name for GPS stations
!
INTEGER           :: IIU,IJU,IKU,IIB,IJB,IKB,IIE,IJE,IKE ! Arrays bounds
! 
INTEGER                :: JLOOP,JI,JJ,JK,JSV,JT,JH,JV,JEL    ! loop index
INTEGER :: IMI ! Current model index
! 
REAL :: ZRV_OV_RD !  XRV / XRD
REAL :: ZGAMREF   ! Standard atmosphere lapse rate (K/m)
REAL :: ZX0D      ! work real scalar  
REAL :: ZLATOR, ZLONOR ! geographical coordinates of 1st mass point
REAL :: ZXHATM, ZYHATM ! conformal    coordinates of 1st mass point
REAL, DIMENSION(:), ALLOCATABLE   :: ZXHAT_ll    !  Position x in the conformal
                                                 ! plane (array on the complete domain)
REAL, DIMENSION(:), ALLOCATABLE   :: ZYHAT_ll    !   Position y in the conformal
                                                 ! plane (array on the complete domain)
!
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))  :: ZPOVO
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))  :: ZTEMP
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))  :: ZVOX,ZVOY,ZVOZ 
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))  :: ZCORIOZ 
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))  :: ZWORK31,ZWORK32
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))  :: ZWORK33,ZWORK34
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2))               :: ZWORK21,ZWORK22
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2))               :: ZWORK23,ZWORK24
REAL,DIMENSION(:,:,:,:,:), ALLOCATABLE                  :: ZWORK42 ! reflectivity on a cartesian grid (PREFL_CART)
REAL,DIMENSION(:,:,:,:,:), ALLOCATABLE                  :: ZWORK42_BIS
REAL,DIMENSION(:,:,:), ALLOCATABLE                      :: ZWORK43 ! latlon coordinates of cartesian grid points (PLATLON)
REAL,DIMENSION(:,:,:), ALLOCATABLE                      :: ZPHI,ZTHETAE,ZTHETAV
REAL,DIMENSION(:,:,:), ALLOCATABLE                      :: ZTHETAES,ZTHETAL,ZTHETAS1,ZTHETAS2
REAL,DIMENSION(:,:,:), ALLOCATABLE                      :: ZVISIKUN,ZVISIGUL,ZVISIZHA 
INTEGER, DIMENSION(:,:), ALLOCATABLE                    :: IWORK1
integer :: ICURR,INBOUT,IERR
!
REAL,DIMENSION(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NSP+NCARB+NSOA,JPMODE):: ZPTOTA
REAL,DIMENSION(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_DST*2):: ZSDSTDEP
REAL,DIMENSION(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_SLT*2):: ZSSLTDEP
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZSIG_DST, ZRG_DST, ZN0_DST
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZSIG_SLT, ZRG_SLT, ZN0_SLT
REAL,DIMENSION(:,:,:), ALLOCATABLE    :: ZBET_SNW, ZRG_SNW
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZMA_SNW
REAL,DIMENSION(:,:,:), ALLOCATABLE  :: ZRHOT, ZTMP ! work array
!
! GBOTUP = True does clustering from bottom up to top, False top down to surface
LOGICAL                                                   :: GBOTUP ! clustering propagation
LOGICAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) :: GCLOUD ! mask
INTEGER,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) :: ICLUSTERID, ICLUSTERLV
REAL,   DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) :: ZCLDSIZE

!ECRITURE DANS UN FICHIER ASCII DE RESULTATS 
!INITIALISATION DU NOM DE FICHIER CREE EN PARALLELE AVEC CELUI LFI
TYPE(TFILEDATA),POINTER :: TZRSFILE
INTEGER :: ILURS
CHARACTER(LEN=32) :: YRS
CHARACTER(LEN=3),DIMENSION(:),ALLOCATABLE  :: YRAD
CHARACTER(LEN=2*INT(NBSTEPMAX*XSTEP_RAD/XGRID)*2*9+1), DIMENSION(:), ALLOCATABLE :: CLATLON
CHARACTER(LEN=2*9) :: CBUFFER
CHARACTER(LEN=4)  :: YELEV
CHARACTER(LEN=3)  :: YGRID_SIZE
INTEGER :: IEL,IIELV
CHARACTER(LEN=5)  :: YVIEW   ! Upward or Downward integration
INTEGER           :: IACCMODE
!
!-------------------------------------------------------------------------------
INTEGER :: IAUX ! work variable 
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) :: ZWORK35,ZWORK36, ZW1, ZW2, ZW3
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2))              :: ZWORK25,ZWORK26
REAL    :: ZEAU ! Mean precipitable water
INTEGER, DIMENSION(SIZE(XZZ,1),SIZE(XZZ,2))          ::IKTOP ! level in which is the altitude 3000m
REAL, DIMENSION(SIZE(XZZ,1),SIZE(XZZ,2),SIZE(XZZ,3)) :: ZDELTAZ ! interval (m) between two levels K
INTEGER :: ILUOUT0 ! Logical unit number for output-listing
!
CHARACTER(LEN=2)  :: INDICE
CHARACTER(LEN=100) :: YMSG
INTEGER           :: IID
TYPE(TFIELDDATA)              :: TZFIELD
TYPE(TFIELDDATA),DIMENSION(2) :: TZFIELD2
!
! LIMA LIDAR
REAL,DIMENSION(:,:,:,:), ALLOCATABLE :: ZTMP1, ZTMP2, ZTMP3, ZTMP4
!
! hauteur couche limite
REAL,DIMENSION(:,:,:),ALLOCATABLE :: ZZZ_GRID1
REAL,DIMENSION(:,:),ALLOCATABLE :: ZTHVSOL,ZSHMIX
REAL,DIMENSION(:,:,:),ALLOCATABLE :: ZZONWIND,ZMERWIND,ZFFWIND2,ZRIB
!
!-------------------------------------------------------------------------------
!
!*       0.     ARRAYS BOUNDS INITIALIZATION
!
CALL GET_DIM_EXT_ll ('B',IIU,IJU)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKU=NKMAX+2*JPVEXT
IKB=1+JPVEXT
IKE=IKU-JPVEXT

IMI = GET_CURRENT_MODEL_INDEX()
ILUOUT0 = TLUOUT0%NLU
TZRSFILE => NULL()
!-------------------------------------------------------------------------------
!
!*       1.     WRITES IN THE LFI FILE
!               ---------------------- 
!
!*       1.0    TPFILE%CNAME and HDADFILE :
!
CALL IO_Field_write(TPFILE,'MASDEV',   NMASDEV)
CALL IO_Field_write(TPFILE,'BUGFIX',   NBUGFIX)
CALL IO_Field_write(TPFILE,'BIBUSER',  CBIBUSER)
CALL IO_Field_write(TPFILE,'PROGRAM',  CPROGRAM)
!
CALL IO_Field_write(TPFILE,'L1D',      L1D)
CALL IO_Field_write(TPFILE,'L2D',      L2D)
CALL IO_Field_write(TPFILE,'PACK',     LPACK)
!
CALL IO_Field_write(TPFILE,'MY_NAME',  TPFILE%CNAME)
CALL IO_Field_write(TPFILE,'DAD_NAME', HDADFILE)
!
IF (LEN_TRIM(HDADFILE)>0) THEN
  CALL IO_Field_write(TPFILE,'DXRATIO',NDXRATIO_ALL(1))
  CALL IO_Field_write(TPFILE,'DYRATIO',NDYRATIO_ALL(1))
  CALL IO_Field_write(TPFILE,'XOR',    NXOR_ALL(1))
  CALL IO_Field_write(TPFILE,'YOR',    NYOR_ALL(1))
END IF
!
CALL IO_Field_write(TPFILE,'SURF',     CSURF)
!
!*       1.1    Type and Dimensions :
!
CALL IO_Field_write(TPFILE,'STORAGE_TYPE','DI')
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
  ALLOCATE(ZXHAT_ll(NIMAX_ll+ 2 * JPHEXT),ZYHAT_ll(NJMAX_ll+2 * JPHEXT))
  CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,IRESP) !//
  CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,IRESP) !//
  ZXHATM = 0.5 * (ZXHAT_ll(1)+ZXHAT_ll(2))
  ZYHATM = 0.5 * (ZYHAT_ll(1)+ZYHAT_ll(2))
  CALL SM_LATLON(XLATORI,XLONORI,ZXHATM,ZYHATM,ZLATOR,ZLONOR)
  DEALLOCATE(ZXHAT_ll,ZYHAT_ll)
!
  CALL IO_Field_write(TPFILE,'LONOR',ZLONOR)
  CALL IO_Field_write(TPFILE,'LATOR',ZLATOR)
!
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
CALL IO_Field_write(TPFILE,'ZS',   XZS)
CALL IO_Field_write(TPFILE,'ZWS',  XZWS)
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
CALL IO_Field_write(TPFILE,'CARTESIAN',LCARTESIAN)
CALL IO_Field_write(TPFILE,'LBOUSS',   LBOUSS)
CALL IO_Field_write(TPFILE,'LOCEAN',   LOCEAN)
CALL IO_Field_write(TPFILE,'LCOUPLES', LCOUPLES)
!
IF (LCARTESIAN .AND. LWIND_ZM) THEN
  LWIND_ZM=.FALSE.
  PRINT*,'YOU ARE IN CARTESIAN GEOMETRY SO LWIND_ZM IS FORCED TO FALSE'
END IF
!*       1.4    Reference state variables :
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
CALL IO_Field_write(TPFILE,'RHODREF',XRHODREF)
CALL IO_Field_write(TPFILE,'THVREF', XTHVREF)
!
!
!*       1.5    Variables necessary for plots
!
! PABST,THT,POVOM for cross sections at constant pressure 
! level or constant theta level or constant PV level
!
IF (INDEX(CISO,'PR') /= 0) THEN
  CALL IO_Field_write(TPFILE,'PABST',XPABST)
END IF
!
IF (INDEX(CISO,'TK') /= 0) THEN
  CALL IO_Field_write(TPFILE,'THT',XTHT)
END IF
!
ZCORIOZ(:,:,:)=SPREAD( XCORIOZ(:,:),DIM=3,NCOPIES=IKU )
ZVOX(:,:,:)=GY_W_VW(XWT,XDYY,XDZZ,XDZY)-GZ_V_VW(XVT,XDZZ)
ZVOX(:,:,2)=ZVOX(:,:,3)
ZVOY(:,:,:)=GZ_U_UW(XUT,XDZZ)-GX_W_UW(XWT,XDXX,XDZZ,XDZX)
ZVOY(:,:,2)=ZVOY(:,:,3)
ZVOZ(:,:,:)=GX_V_UV(XVT,XDXX,XDZZ,XDZX)-GY_U_UV(XUT,XDYY,XDZZ,XDZY)
ZVOZ(:,:,2)=ZVOZ(:,:,3)
ZVOZ(:,:,1)=ZVOZ(:,:,3)
ZWORK31(:,:,:)=GX_M_M(XTHT,XDXX,XDZZ,XDZX)
ZWORK32(:,:,:)=GY_M_M(XTHT,XDYY,XDZZ,XDZY)
ZWORK33(:,:,:)=GZ_M_M(XTHT,XDZZ)
ZPOVO(:,:,:)= ZWORK31(:,:,:)*MZF(MYF(ZVOX(:,:,:)))     &
             + ZWORK32(:,:,:)*MZF(MXF(ZVOY(:,:,:)))     &
             + ZWORK33(:,:,:)*(MYF(MXF(ZVOZ(:,:,:))) + ZCORIOZ(:,:,:))
ZPOVO(:,:,:)= ZPOVO(:,:,:)*1E6/XRHODREF(:,:,:)
ZPOVO(:,:,1)  =-1.E+11
ZPOVO(:,:,IKU)=-1.E+11
IF (INDEX(CISO,'EV') /= 0) THEN
  TZFIELD%CMNHNAME   = 'POVOT'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'POVOT'
  TZFIELD%CUNITS     = 'PVU' ! 1 PVU = 1e-6 m^2 s^-1 K kg^-1
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_POtential VOrticity'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZPOVO)
END IF
!
!
IF (LVAR_RS) THEN
  CALL IO_Field_write(TPFILE,'UT',XUT)
  CALL IO_Field_write(TPFILE,'VT',XVT)
  !
  IF (LWIND_ZM) THEN
    TZFIELD2(1)%CMNHNAME   = 'UM_ZM'
    TZFIELD2(1)%CSTDNAME   = ''
    TZFIELD2(1)%CLONGNAME  = 'UM_ZM'
    TZFIELD2(1)%CUNITS     = 'm s-1'
    TZFIELD2(1)%CDIR       = 'XY'
    TZFIELD2(1)%CCOMMENT   = 'Zonal component of horizontal wind'
    TZFIELD2(1)%NGRID      = 2
    TZFIELD2(1)%NTYPE      = TYPEREAL
    TZFIELD2(1)%NDIMS      = 3
    TZFIELD2(1)%LTIMEDEP   = .TRUE.
    !
    TZFIELD2(2)%CMNHNAME   = 'VM_ZM'
    TZFIELD2(2)%CSTDNAME   = ''
    TZFIELD2(2)%CLONGNAME  = 'VM_ZM'
    TZFIELD2(2)%CUNITS     = 'm s-1'
    TZFIELD2(2)%CDIR       = 'XY'
    TZFIELD2(2)%CCOMMENT   = 'Meridian component of horizontal wind'
    TZFIELD2(2)%NGRID      = 3
    TZFIELD2(2)%NTYPE      = TYPEREAL
    TZFIELD2(2)%NDIMS      = 3
    TZFIELD2(2)%LTIMEDEP   = .TRUE.
    !
    CALL UV_TO_ZONAL_AND_MERID(XUT,XVT,23,TPFILE=TPFILE,TZFIELDS=TZFIELD2)
  END IF
  !
  CALL IO_Field_write(TPFILE,'WT',XWT)
  !
  !   write mixing ratio for water vapor required to plot radio-soundings
  !
  IF (LUSERV) THEN
    CALL IO_Field_write(TPFILE,'RVT',XRT(:,:,:,IDX_RVT))
  END IF
END IF
!
!*   Latitude and Longitude arrays
!
IF (.NOT.LCARTESIAN) THEN
  CALL IO_Field_write(TPFILE,'LAT',XLAT)
  CALL IO_Field_write(TPFILE,'LON',XLON)
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       1.6    Other pronostic variables
!
ZTEMP(:,:,:)=XTHT(:,:,:)*(XPABST(:,:,:)/ XP00) **(XRD/XCPD)
!
IF (LVAR_TURB) THEN
  IF (CTURB /= 'NONE') THEN
    CALL IO_Field_write(TPFILE,'TKET',XTKET)
    !
    IF( NRR > 1 ) THEN
      CALL IO_Field_write(TPFILE,'SRCT',XSRCT)
      CALL IO_Field_write(TPFILE,'SIGS',XSIGS)
    END IF
    ! 
    IF(CTOM=='TM06') THEN
      CALL IO_Field_write(TPFILE,'BL_DEPTH',XBL_DEPTH)
    END IF
  END IF
END IF
!
!* Rains
!
IF (LVAR_PR .AND. LUSERR .AND. SIZE(XINPRR)>0 ) THEN
  !
  ! explicit species
  !
  CALL FIND_FIELD_ID_FROM_MNHNAME('INPRR',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CUNITS = 'mm hour-1'
  CALL IO_Field_write(TPFILE,TZFIELD,XINPRR*3.6E6)
  !
  CALL IO_Field_write(TPFILE,'INPRR3D',XINPRR3D)
  CALL IO_Field_write(TPFILE,'EVAP3D', XEVAP3D)
  !
  CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRR',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CUNITS = 'mm'
  CALL IO_Field_write(TPFILE,TZFIELD,XACPRR*1.0E3)
  !
  IF (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR.&
      CCLOUD == 'KHKO' .OR. CCLOUD == 'LIMA') THEN 
    IF (SIZE(XINPRC) /= 0 ) THEN
      CALL FIND_FIELD_ID_FROM_MNHNAME('INPRC',IID,IRESP)
      TZFIELD = TFIELDLIST(IID)
      TZFIELD%CUNITS = 'mm hour-1'
      CALL IO_Field_write(TPFILE,TZFIELD,XINPRC*3.6E6)
      !
      CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRC',IID,IRESP)
      TZFIELD = TFIELDLIST(IID)
      TZFIELD%CUNITS = 'mm'
      CALL IO_Field_write(TPFILE,TZFIELD,XACPRC*1.0E3)
    END IF 
    IF (SIZE(XINDEP) /= 0 ) THEN
      CALL FIND_FIELD_ID_FROM_MNHNAME('INDEP',IID,IRESP)
      TZFIELD = TFIELDLIST(IID)
      TZFIELD%CUNITS = 'mm hour-1'
      CALL IO_Field_write(TPFILE,TZFIELD,XINDEP*3.6E6)
      !
      CALL FIND_FIELD_ID_FROM_MNHNAME('ACDEP',IID,IRESP)
      TZFIELD = TFIELDLIST(IID)
      TZFIELD%CUNITS = 'mm'
      CALL IO_Field_write(TPFILE,TZFIELD,XACDEP*1.0E3)
    END IF 
  END IF 
  IF (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'LIMA') THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRS',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XINPRS*3.6E6)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRS',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_write(TPFILE,TZFIELD,XACPRS*1.0E3)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRG',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XINPRG*3.6E6)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRG',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_write(TPFILE,TZFIELD,XACPRG*1.0E3)
  !
    IF (SIZE(XINPRH) /= 0 ) THEN
      CALL FIND_FIELD_ID_FROM_MNHNAME('INPRH',IID,IRESP)
      TZFIELD = TFIELDLIST(IID)
      TZFIELD%CUNITS = 'mm hour-1'
      CALL IO_Field_write(TPFILE,TZFIELD,XINPRH*3.6E6)
      !
      CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRH',IID,IRESP)
      TZFIELD = TFIELDLIST(IID)
      TZFIELD%CUNITS = 'mm'
      CALL IO_Field_write(TPFILE,TZFIELD,XACPRH*1.0E3)
    ENDIF
  !
    ZWORK21(:,:) = XINPRR(:,:) + XINPRS(:,:) + XINPRG(:,:)
    IF (SIZE(XINPRC) /= 0 ) &     
      ZWORK21(:,:) = ZWORK21(:,:) + XINPRC(:,:)
    IF (SIZE(XINPRH) /= 0 ) &       
      ZWORK21(:,:) = ZWORK21(:,:) + XINPRH(:,:)
    CALL FIND_FIELD_ID_FROM_MNHNAME('INPRT',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21*3.6E6)
  !
    ZWORK21(:,:) = XACPRR(:,:) + XACPRS(:,:) + XACPRG(:,:)
    IF (SIZE(XINPRC) /= 0 ) &      
      ZWORK21(:,:) = ZWORK21(:,:) + XACPRC(:,:)
    IF (SIZE(XINPRH) /= 0 ) &        
      ZWORK21(:,:) = ZWORK21(:,:) + XACPRH(:,:)
  !
    CALL FIND_FIELD_ID_FROM_MNHNAME('ACPRT',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21*1.0E3)
  !
  END IF
  !
  !* Convective rain
  !
  IF (CDCONV /= 'NONE') THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('PRCONV',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XPRCONV*3.6E6)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('PACCONV',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm'
    CALL IO_Field_write(TPFILE,TZFIELD,XPACCONV*1.0E3)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('PRSCONV',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'mm hour-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XPRSCONV*3.6E6)
  END IF
END IF
IF (LVAR_PR ) THEN
  !Precipitable water in kg/m**2 
  ZWORK21(:,:) = 0.
  ZWORK22(:,:) = 0.
  ZWORK23(:,:) = 0.
  ZWORK31(:,:,:) = DZF(XZZ(:,:,:))
  DO JK = IKB,IKE
    !* Calcul de qtot
    IF  (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'LIMA') THEN
      ZWORK23(IIB:IIE,IJB:IJE) = XRT(IIB:IIE,IJB:IJE,JK,1) + &
      XRT(IIB:IIE,IJB:IJE,JK,2) + XRT(IIB:IIE,IJB:IJE,JK,3) + &
      XRT(IIB:IIE,IJB:IJE,JK,4) + XRT(IIB:IIE,IJB:IJE,JK,5) + &
      XRT(IIB:IIE,IJB:IJE,JK,6)
    ELSE
      ZWORK23(IIB:IIE,IJB:IJE) = XRT(IIB:IIE,IJB:IJE,JK,1)
    ENDIF
    !* Calcul de l'eau precipitable
    ZWORK21(IIB:IIE,IJB:IJE)=XRHODREF(IIB:IIE,IJB:IJE,JK)* &
    ZWORK23(IIB:IIE,IJB:IJE)* ZWORK31(IIB:IIE,IJB:IJE,JK)
    !* Sum 
    ZWORK22(IIB:IIE,IJB:IJE) = ZWORK22(IIB:IIE,IJB:IJE)+ZWORK21(IIB:IIE,IJB:IJE)
    ZWORK21(:,:) = 0.
    ZWORK23(:,:) = 0.
  END DO
  !* Precipitable water in kg/m**2
  TZFIELD%CMNHNAME   = 'PRECIP_WAT'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'PRECIP_WAT'
  TZFIELD%CUNITS     = 'kg m-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK22)
ENDIF
!
!
!* Flux d'humidite et d'hydrometeores
IF (LHU_FLX) THEN
  ZWORK35(:,:,:) = XRHODREF(:,:,:) * XRT(:,:,:,1)
  ZWORK31(:,:,:) = MXM(ZWORK35(:,:,:)) * XUT(:,:,:)
  ZWORK32(:,:,:) = MYM(ZWORK35(:,:,:)) * XVT(:,:,:)
  ZWORK35(:,:,:) = GX_U_M(ZWORK31,XDXX,XDZZ,XDZX) + GY_V_M(ZWORK32,XDYY,XDZZ,XDZY)
  IF  (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'LIMA') THEN
    ZWORK36(:,:,:) = ZWORK35(:,:,:) + XRHODREF(:,:,:) * (XRT(:,:,:,2) + &
    XRT(:,:,:,3) + XRT(:,:,:,4) + XRT(:,:,:,5) + XRT(:,:,:,6))
    ZWORK33(:,:,:) = MXM(ZWORK36(:,:,:)) * XUT(:,:,:)
    ZWORK34(:,:,:) = MYM(ZWORK36(:,:,:)) * XVT(:,:,:)
    ZWORK36(:,:,:) = GX_U_M(ZWORK33,XDXX,XDZZ,XDZX) + GY_V_M(ZWORK34,XDYY,XDZZ,XDZY)
  ENDIF
  !
  ! Integration sur 3000 m
  !
  IKTOP(:,:)=0
  DO JK=1,IKU-1
    WHERE (((XZZ(:,:,JK) -XZS(:,:))<= 3000.0) .AND. ((XZZ(:,:,JK+1) -XZS(:,:))> 3000.0))
      IKTOP(:,:)=JK
    END WHERE
  END DO
  ZDELTAZ(:,:,:)=DZF(XZZ)
  ZWORK21(:,:) = 0.
  ZWORK22(:,:) = 0.
  ZWORK25(:,:) = 0.  
  DO JJ=1,IJU
    DO JI=1,IIU
      IAUX=IKTOP(JI,JJ)
      DO JK=IKB,IAUX-1 
        ZWORK21(JI,JJ) = ZWORK21(JI,JJ) + ZWORK31(JI,JJ,JK) * ZDELTAZ(JI,JJ,JK)
        ZWORK22(JI,JJ) = ZWORK22(JI,JJ) + ZWORK32(JI,JJ,JK) * ZDELTAZ(JI,JJ,JK)
        ZWORK25(JI,JJ) = ZWORK25(JI,JJ) + ZWORK35(JI,JJ,JK) * ZDELTAZ(JI,JJ,JK)
      ENDDO
      IF (IAUX >= IKB) THEN
        ZDELTAZ(JI,JJ,IAUX)= 3000. - (XZZ(JI,JJ,IAUX) -XZS(JI,JJ))
        ZWORK21(JI,JJ) = ZWORK21(JI,JJ) + ZWORK31(JI,JJ,IAUX) * ZDELTAZ(JI,JJ,IAUX) 
        ZWORK22(JI,JJ) = ZWORK22(JI,JJ) + ZWORK32(JI,JJ,IAUX) * ZDELTAZ(JI,JJ,IAUX)
        ZWORK25(JI,JJ) = ZWORK25(JI,JJ) + ZWORK35(JI,JJ,IAUX) * ZDELTAZ(JI,JJ,IAUX)
      ENDIF
    ENDDO
  ENDDO
  IF  (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'LIMA') THEN
    ZWORK23(:,:) = 0.
    ZWORK24(:,:) = 0.
    ZWORK26(:,:) = 0.
    DO JJ=1,IJU
      DO JI=1,IIU
        IAUX=IKTOP(JI,JJ)
        DO JK=IKB,IAUX-1 
          ZWORK23(JI,JJ) = ZWORK23(JI,JJ) + ZWORK33(JI,JJ,JK) * ZDELTAZ(JI,JJ,JK)
          ZWORK24(JI,JJ) = ZWORK24(JI,JJ) + ZWORK34(JI,JJ,JK) * ZDELTAZ(JI,JJ,JK)
          ZWORK26(JI,JJ) = ZWORK26(JI,JJ) + ZWORK36(JI,JJ,JK) * ZDELTAZ(JI,JJ,JK)
        ENDDO
        IF (IAUX >= IKB) THEN
          ZDELTAZ(JI,JJ,IAUX)= 3000. - (XZZ(JI,JJ,IAUX) -XZS(JI,JJ))
          ZWORK23(JI,JJ) = ZWORK23(JI,JJ) + ZWORK33(JI,JJ,IAUX) * ZDELTAZ(JI,JJ,IAUX) 
          ZWORK24(JI,JJ) = ZWORK24(JI,JJ) + ZWORK34(JI,JJ,IAUX) * ZDELTAZ(JI,JJ,IAUX)
          ZWORK26(JI,JJ) = ZWORK26(JI,JJ) + ZWORK36(JI,JJ,IAUX) * ZDELTAZ(JI,JJ,IAUX)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  ! Ecriture
  !  composantes U et V du flux surfacique d'humidite 
  TZFIELD%CMNHNAME   = 'UM90'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'UM90'
  TZFIELD%CUNITS     = 'kg s-1 m-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = ''
  TZFIELD%NGRID      = 2
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  !  
  TZFIELD%CMNHNAME   = 'VM90'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'VM90'
  TZFIELD%CUNITS     = 'kg s-1 m-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = ''
  TZFIELD%NGRID      = 3
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK32)
  !  composantes U et V du flux d'humidite integre sur 3000 metres
  TZFIELD%CMNHNAME   = 'UM91'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'UM91'
  TZFIELD%CUNITS     = 'kg s-1 m-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = ''
  TZFIELD%NGRID      = 2
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  !
  TZFIELD%CMNHNAME   = 'VM91'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'VM91'
  TZFIELD%CUNITS     = 'kg s-1 m-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = ''
  TZFIELD%NGRID      = 3
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK22)
  !
  !   Convergence d'humidite
  TZFIELD%CMNHNAME   = 'HMCONV'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'HMCONV'
  TZFIELD%CUNITS     = 'kg s-1 m-3'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Horizontal CONVergence of moisture flux'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,-ZWORK35)
  !
  !   Convergence d'humidite integre sur 3000 metres
  TZFIELD%CMNHNAME   = 'HMCONV3000'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'HMCONV3000'
  TZFIELD%CUNITS     = 'kg s-1 m-3'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Horizontal CONVergence of moisture flux'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,-ZWORK25)
  !
  IF  (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'LIMA') THEN
    !  composantes U et V du flux surfacique d'hydrometeores  
    TZFIELD%CMNHNAME   = 'UM92'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'UM92'
    TZFIELD%CUNITS     = 'kg s-1 m-2'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = ''
    TZFIELD%NGRID      = 2
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK33)
    ! 
    TZFIELD%CMNHNAME   = 'VM92'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'VM92'
    TZFIELD%CUNITS     = 'kg s-1 m-2'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = ''
    TZFIELD%NGRID      = 3
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK34)
    !  composantes U et V du flux d'hydrometeores integre sur 3000 metres
    TZFIELD%CMNHNAME   = 'UM93'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'UM93'
    TZFIELD%CUNITS     = 'kg s-1 m-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = ''
    TZFIELD%NGRID      = 2
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK23)
    ! 
    TZFIELD%CMNHNAME   = 'VM93'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'VM93'
    TZFIELD%CUNITS     = 'kg s-1 m-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = ''
    TZFIELD%NGRID      = 3
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK24)
    !   Convergence d'hydrometeores
    TZFIELD%CMNHNAME   = 'HMCONV_TT'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'HMCONV_TT'
    TZFIELD%CUNITS     = 'kg s-1 m-3'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Horizontal CONVergence of hydrometeor flux'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,-ZWORK36)
    !   Convergence d'hydrometeores integre sur 3000 metres
    TZFIELD%CMNHNAME   = 'HMCONV3000_TT'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'HMCONV3000_TT'
    TZFIELD%CUNITS     = 'kg s-1 m-3'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Horizontal CONVergence of hydrometeor flux'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,-ZWORK26)
  ENDIF
ENDIF
!
!* Moist variables
!
IF (LVAR_MRW .OR. LLIMA_DIAG) THEN
  IF (NRR >=1) THEN
    ! Moist variables are written individually in file
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    IF (LUSERV) THEN
      TZFIELD%CMNHNAME   = 'MRV'
      TZFIELD%CLONGNAME  = 'MRV'
      TZFIELD%CUNITS     = 'g kg-1'
      TZFIELD%CCOMMENT   = 'X_Y_Z_MRV'
      CALL IO_Field_write(TPFILE,TZFIELD,XRT(:,:,:,IDX_RVT)*1.E3)
    END IF
    IF (LUSERC) THEN
      TZFIELD%CMNHNAME   = 'MRC'
      TZFIELD%CLONGNAME  = 'MRC'
      TZFIELD%CUNITS     = 'g kg-1'
      TZFIELD%CCOMMENT   = 'X_Y_Z_MRC'
      CALL IO_Field_write(TPFILE,TZFIELD,XRT(:,:,:,IDX_RCT)*1.E3)
!
      TZFIELD%CMNHNAME   = 'VRC'
      TZFIELD%CLONGNAME  = 'VRC'
      TZFIELD%CUNITS     = '1' !vol/vol
      TZFIELD%CCOMMENT   = 'X_Y_Z_VRC (vol/vol)'
      CALL IO_Field_write(TPFILE,TZFIELD,XRT(:,:,:,IDX_RCT)*XRHODREF(:,:,:)/1.E3)
    END IF
    IF (LUSERR) THEN
      TZFIELD%CMNHNAME   = 'MRR'
      TZFIELD%CLONGNAME  = 'MRR'
      TZFIELD%CUNITS     = 'g kg-1'
      TZFIELD%CCOMMENT   = 'X_Y_Z_MRR'
      CALL IO_Field_write(TPFILE,TZFIELD,XRT(:,:,:,IDX_RRT)*1.E3)
!
      TZFIELD%CMNHNAME   = 'VRR'
      TZFIELD%CLONGNAME  = 'VRR'
      TZFIELD%CUNITS     = '1' !vol/vol
      TZFIELD%CCOMMENT   = 'X_Y_Z_VRR (vol/vol)'
      CALL IO_Field_write(TPFILE,TZFIELD,XRT(:,:,:,IDX_RRT)*XRHODREF(:,:,:)/1.E3)
    END IF
    IF (LUSERI) THEN
      TZFIELD%CMNHNAME   = 'MRI'
      TZFIELD%CLONGNAME  = 'MRI'
      TZFIELD%CUNITS     = 'g kg-1'
      TZFIELD%CCOMMENT   = 'X_Y_Z_MRI'
      CALL IO_Field_write(TPFILE,TZFIELD,XRT(:,:,:,IDX_RIT)*1.E3)
!
      IF (LUSECI) THEN
        CALL IO_Field_write(TPFILE,'CIT',XCIT(:,:,:))
      END IF
    END IF
    IF (LUSERS) THEN
      TZFIELD%CMNHNAME   = 'MRS'
      TZFIELD%CLONGNAME  = 'MRS'
      TZFIELD%CUNITS     = 'g kg-1'
      TZFIELD%CCOMMENT   = 'X_Y_Z_MRS'
      CALL IO_Field_write(TPFILE,TZFIELD,XRT(:,:,:,IDX_RST)*1.E3)
    END IF
    IF (LUSERG) THEN
      TZFIELD%CMNHNAME   = 'MRG'
      TZFIELD%CLONGNAME  = 'MRG'
      TZFIELD%CUNITS     = 'g kg-1'
      TZFIELD%CCOMMENT   = 'X_Y_Z_MRG'
      CALL IO_Field_write(TPFILE,TZFIELD,XRT(:,:,:,IDX_RGT)*1.E3)
    END IF
    IF (LUSERH) THEN
      TZFIELD%CMNHNAME   = 'MRH'
      TZFIELD%CLONGNAME  = 'MRH'
      TZFIELD%CUNITS     = 'g kg-1'
      TZFIELD%CCOMMENT   = 'X_Y_Z_MRH'
      CALL IO_Field_write(TPFILE,TZFIELD,XRT(:,:,:,IDX_RHT)*1.E3)
    END IF
  END IF
END IF
!
!* Scalar Variables
!
! User scalar variables
! individually in the file
IF (LVAR_MRSV) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'g kg-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = 1,NSV_USER
    WRITE(TZFIELD%CMNHNAME,'(A4,I3.3)')'MRSV',JSV
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3)')'X_Y_Z_','MRSV',JSV
    CALL IO_Field_write(TPFILE,TZFIELD,XSVT(:,:,:,JSV)*1.E3)
  END DO
END IF
! microphysical C2R2 scheme scalar variables
IF(LVAR_MRW) THEN
  IF (NSV_C2R2END>=NSV_C2R2BEG) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_C2R2BEG,NSV_C2R2END
      TZFIELD%CMNHNAME   = TRIM(C2R2NAMES(JSV-NSV_C2R2BEG+1))//'T'
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      IF (JSV < NSV_C2R2END) THEN
        TZFIELD%CUNITS     = 'cm-3'
        ZWORK31(:,:,:)=XSVT(:,:,:,JSV)*1.E-6
      ELSE
        TZFIELD%CUNITS     = 'l-1'
        ZWORK31(:,:,:)=XSVT(:,:,:,JSV)*1.E-3
      END IF
      WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3)')'X_Y_Z_','MRSV',JSV
      CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    END DO
  END IF
  ! microphysical C3R5 scheme additional scalar variables
  IF (NSV_C1R3END>=NSV_C1R3BEG) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'l-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_C1R3BEG,NSV_C1R3END
      TZFIELD%CMNHNAME   = TRIM(C1R3NAMES(JSV-NSV_C1R3BEG+1))//'T'
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
      CALL IO_Field_write(TPFILE,TZFIELD,XSVT(:,:,:,JSV)*1.E-3)
    END DO
  END IF
END IF
!
! microphysical LIMA scheme scalar variables
!
IF (LLIMA_DIAG) THEN
  IF (NSV_LIMA_END>=NSV_LIMA_BEG) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
  END IF
  !
  DO JSV = NSV_LIMA_BEG,NSV_LIMA_END
    !
    TZFIELD%CUNITS     = 'cm-3'
    WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
    !
! Nc
    IF (JSV .EQ. NSV_LIMA_NC) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_WARM_CONC(1))//'T'
    END IF
! Nr
    IF (JSV .EQ. NSV_LIMA_NR) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_WARM_CONC(2))//'T'
    END IF
! N CCN free
    IF (JSV .GE. NSV_LIMA_CCN_FREE .AND. JSV .LT. NSV_LIMA_CCN_ACTI) THEN
      WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_CCN_FREE + 1)
      TZFIELD%CMNHNAME   = TRIM(CLIMA_WARM_CONC(3))//INDICE//'T'
    END IF
! N CCN acti
    IF (JSV .GE. NSV_LIMA_CCN_ACTI .AND. JSV .LT. NSV_LIMA_CCN_ACTI + NMOD_CCN) THEN
      WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_CCN_ACTI + 1)
      TZFIELD%CMNHNAME   = TRIM(CLIMA_WARM_CONC(4))//INDICE//'T'
    END IF
! Scavenging
    IF (JSV .EQ. NSV_LIMA_SCAVMASS) THEN
      TZFIELD%CMNHNAME   = TRIM(CAERO_MASS(1))//'T'
      TZFIELD%CUNITS     = 'kg cm-3'
    END IF
! Ni
    IF (JSV .EQ. NSV_LIMA_NI) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_CONC(1))//'T'
    END IF
! Ns
    IF (JSV .EQ. NSV_LIMA_NS) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_CONC(2))//'T'
    END IF
! Ng
    IF (JSV .EQ. NSV_LIMA_NG) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_CONC(3))//'T'
    END IF
! Nh
    IF (JSV .EQ. NSV_LIMA_NH) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_CONC(4))//'T'
    END IF
! N IFN free
    IF (JSV .GE. NSV_LIMA_IFN_FREE .AND. JSV .LT. NSV_LIMA_IFN_NUCL) THEN
      WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_IFN_FREE + 1)
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_CONC(5))//INDICE//'T'
    END IF
! N IFN nucl
    IF (JSV .GE. NSV_LIMA_IFN_NUCL .AND. JSV .LT. NSV_LIMA_IFN_NUCL + NMOD_IFN) THEN
      WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_IFN_NUCL + 1)
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_CONC(6))//INDICE//'T'
    END IF
! N IMM nucl
    IF (JSV .GE. NSV_LIMA_IMM_NUCL .AND. JSV .LT. NSV_LIMA_IMM_NUCL + NMOD_IMM) THEN
      WRITE(INDICE,'(I2.2)')(NINDICE_CCN_IMM(JSV - NSV_LIMA_IMM_NUCL + 1))
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_CONC(7))//INDICE//'T'
    END IF
! Hom. freez. of CCN
    IF (JSV .EQ. NSV_LIMA_HOM_HAZE) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_COLD_CONC(8))//'T'
    END IF
    !
! Supersaturation          
    IF (JSV .EQ. NSV_LIMA_SPRO) THEN
      TZFIELD%CMNHNAME   = TRIM(CLIMA_WARM_CONC(5))//'T'
    END IF
    !
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    ZWORK31(:,:,:)=XSVT(:,:,:,JSV)*1.E-6*XRHODREF(:,:,:)
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END DO
!
  IF (LUSERC) THEN
    TZFIELD%CMNHNAME   = 'LWC'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'LWC'
    TZFIELD%CUNITS     = 'g m-3'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_LWC'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    ZWORK31(:,:,:)=XRT(:,:,:,2)*1.E3*XRHODREF(:,:,:)
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END IF
!
  IF (LUSERI) THEN
    TZFIELD%CMNHNAME   = 'IWC'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'IWC'
    TZFIELD%CUNITS     = 'g m-3'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_MRI'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    ZWORK31(:,:,:)=XRT(:,:,:,4)*1.E3*XRHODREF(:,:,:)
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END IF
!
END IF
!
! chemical scalar variables in gas phase ppbv
IF (LCHEMDIAG) THEN
  DO JSV = NSV_CHGSBEG,NSV_CHGSEND
    TZFIELD%CMNHNAME   = TRIM(UPCASE(CNAMES(JSV-NSV_CHGSBEG+1)))//'T'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'ppbv'
    TZFIELD%CDIR       = 'XY'
    WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3)')'X_Y_Z_','CHIM',JSV
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,XSVT(:,:,:,JSV)*1.E9)
  END DO
END IF
IF (LCHAQDIAG) THEN    !aqueous concentration in M
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'M'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  ZWORK31(:,:,:)=0.
  DO JSV = NSV_CHACBEG, NSV_CHACBEG-1+NEQAQ/2   !cloud water
    TZFIELD%CMNHNAME   = TRIM(CNAMES(JSV-NSV_CHACBEG+NSV_CHGS+1))//'T'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3)')'X_Y_Z_','CHAQ',JSV
    WHERE(((XRT(:,:,:,2)*XRHODREF(:,:,:))/1.e3) .GE. XRTMIN_AQ)
      ZWORK31(:,:,:)=(XSVT(:,:,:,JSV)*1000.)/(XMD*1.E+3*XRT(:,:,:,2))
    ENDWHERE
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END DO
  !
  ZWORK31(:,:,:)=0.
  DO JSV = NSV_CHACBEG+NEQAQ/2, NSV_CHACEND    !rain water
    TZFIELD%CMNHNAME   = TRIM(CNAMES(JSV-NSV_CHACBEG+NSV_CHGS+1))//'T'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3)')'X_Y_Z_','CHAQ',JSV
    WHERE(((XRT(:,:,:,3)*XRHODREF(:,:,:))/1.e3) .GE. XRTMIN_AQ)
      ZWORK31(:,:,:)=(XSVT(:,:,:,JSV)*1000.)/(XMD*1.E+3*XRT(:,:,:,3))
    ENDWHERE
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END DO
!  ZWORK31(:,:,:)=0.
!  DO JSV = NSV_CHICBEG,NSV_CHICEND   ! ice phase
!    TZFIELD%CMNHNAME   = TRIM(CICNAMES(JSV-NSV_CHICBEG+1))//'T'
!    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
!    WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3,A4)')'X_Y_Z_','CHIC',JSV,' (M)'
!    WHERE(((XRT(:,:,:,3)*XRHODREF(:,:,:))/1.e3) .GE. XRTMIN_AQ)
!      ZWORK31(:,:,:)=(XSVT(:,:,:,JSV)*1000.)/(XMD*1.E+3*XRT(:,:,:,3))
!    ENDWHERE
!    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!  END DO
END IF

! Passive polluant scalar variables
IF (LPASPOL) THEN
  ALLOCATE(ZRHOT( SIZE(XTHT,1), SIZE(XTHT,2),SIZE(XTHT,3)))
  ALLOCATE(ZTMP(  SIZE(XTHT,1), SIZE(XTHT,2),SIZE(XTHT,3)))
!
!*	Density                                          
!
  ZRHOT(:,:,:)=XPABST(:,:,:)/(XRD*XTHT(:,:,:)*((XPABST(:,:,:)/XP00)**(XRD/XCPD)))
!
!*	Conversion g/m3.
!
  ZRHOT(:,:,:)=ZRHOT(:,:,:)*1000.0
  !
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'g m-3'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = 1,NSV_PP
    ZTMP(:,:,:)=ABS( XSVT(:,:,:,JSV+NSV_PPBEG-1)*ZRHOT(:,:,:) )
    WRITE(TZFIELD%CMNHNAME,'(A3,I3.3)')'PPT',JSV
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
    CALL IO_Field_write(TPFILE,TZFIELD,ZTMP)
  END DO
  DEALLOCATE(ZTMP)
  DEALLOCATE(ZRHOT)
END IF
! Conditional sampling variables
IF (LCONDSAMP) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_CSBEG,NSV_CSEND
    WRITE(TZFIELD%CMNHNAME,'(A3,I3.3)')'CST',JSV
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)
    CALL IO_Field_write(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
  END DO
END IF
!
!  Blowing snow variables
!
IF(LBLOWSNOW) THEN
  TZFIELD%CMNHNAME   = 'SNWSUBL3D'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
  TZFIELD%CUNITS     = 'kg m-3 s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_INstantaneous 3D Drifting snow sublimation flux'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,XSNWSUBL3D(:,:,:))
  !
  ZWORK21(:,:) = 0.
  DO JK = IKB,IKE
    ZWORK21(:,:) = ZWORK21(:,:)+XSNWSUBL3D(:,:,JK) * &
                  (XZZ(:,:,JK+1)-XZZ(:,:,JK))/XRHOLW*3600*24
  END DO
  ZWORK21(:,:) = ZWORK21(:,:)*1000. ! vapor water in mm unit
  !
  TZFIELD%CMNHNAME   = 'COL_SNWSUBL'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
  TZFIELD%CUNITS     = 'mm day-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Column Sublimation Rate (mmSWE/day)'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21(:,:))
  !
  IF(.NOT.ALLOCATED(ZBET_SNW)) &
        ALLOCATE(ZBET_SNW(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3)))
  IF(.NOT.ALLOCATED(ZRG_SNW))  &
    ALLOCATE(ZRG_SNW(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3)))
  IF(.NOT.ALLOCATED(ZMA_SNW))  &
    ALLOCATE(ZMA_SNW(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3),NBLOWSNOW3D))
  !
  CALL PPP2SNOW(XSVT(:,:,:,NSV_SNWBEG:NSV_SNWEND),XRHODREF,&
               PBET3D=ZBET_SNW, PRG3D=ZRG_SNW, PM3D=ZMA_SNW)
  !
  TZFIELD%CMNHNAME   = 'SNWRGA'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'RG (mean) SNOW'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZRG_SNW(:,:,:))
  !
  TZFIELD%CMNHNAME   = 'SNWBETA'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'BETA SNOW'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZBET_SNW(:,:,:))
  !
  TZFIELD%CMNHNAME   = 'SNWNOA'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
  TZFIELD%CUNITS     = 'm-3'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'NUM CONC SNOW (#/m3)'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZMA_SNW(:,:,:,1))
  !
  TZFIELD%CMNHNAME   = 'SNWMASS'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
  TZFIELD%CUNITS     = 'kg m-3'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'MASS CONC SNOW'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZMA_SNW(:,:,:,2))
  !
  ZWORK21(:,:) = 0.
  DO JK = IKB,IKE
    ZWORK21(:,:) = ZWORK21(:,:)+ZMA_SNW(:,:,JK,2) * &
                   (XZZ(:,:,JK+1)-XZZ(:,:,JK))/XRHOLW
  END DO
  ZWORK21(:,:) = ZWORK21(:,:)*1000. ! vapor water in mm unit
  TZFIELD%CMNHNAME   = 'THDS'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
  TZFIELD%CUNITS     = 'mm'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_THickness of Drifting Snow (mm SWE)'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21(:,:))
END IF
! Lagrangian variables
IF (LTRAJ) THEN
  TZFIELD%CSTDNAME   = ''
  !PW TODO: check units
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_LGBEG,NSV_LGEND
    TZFIELD%CMNHNAME   = TRIM(CLGNAMES(JSV-NSV_LGBEG+1))//'T'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A6,A20,I3.3,A4)')'X_Y_Z_','Lagrangian variable ',JSV,' (M)'
    CALL IO_Field_write(TPFILE,TZFIELD,XSVT(:,:,:,JSV))
  END DO
  ! X coordinate
  DO JK=1,IKU
    DO JJ=1,IJU
      DO JI=1,IIU-1
       ZWORK31(JI,JJ,JK)=0.5*(XXHAT(JI)+XXHAT(JI+1))
      END DO
      ZWORK31(IIU,JJ,JK)=2.*ZWORK31(IIU-1,JJ,JK) - ZWORK31(IIU-2,JJ,JK)
    END DO
  END DO
  TZFIELD%CMNHNAME   = 'X'
  TZFIELD%CLONGNAME  = 'X'
  TZFIELD%CCOMMENT   = 'X_Y_Z_X coordinate'
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  ! Y coordinate
  DO JK=1,IKU
    DO JI=1,IIU
      DO JJ=1,IJU-1
        ZWORK31(JI,JJ,JK)=0.5*(XYHAT(JJ)+XYHAT(JJ+1))
      END DO
      ZWORK31(JI,IJU,JK)=2.*ZWORK31(JI,IJU-1,JK) - ZWORK31(JI,IJU-2,JK)
    END DO
  END DO
  TZFIELD%CMNHNAME   = 'Y'
  TZFIELD%CLONGNAME  = 'Y'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Y coordinate'
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
END IF
! linox scalar variables
IF (.NOT.(LUSECHEM .OR. LCHEMDIAG) .AND. LCH_CONV_LINOX) THEN
  DO JSV = NSV_LNOXBEG,NSV_LNOXEND
    TZFIELD%CMNHNAME   = 'LINOXT'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'ppbv'
    TZFIELD%CDIR       = 'XY'
    WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3)')'X_Y_Z_','LNOX',JSV
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,XSVT(:,:,:,JSV)*1.E9)
  END DO
END IF
IF (LELECDIAG .AND. CELEC .NE. "NONE") THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_ELECBEG,NSV_ELECEND
    TZFIELD%CMNHNAME   = TRIM(CELECNAMES(JSV-NSV_ELECBEG+1))//'T'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    IF (JSV .GT. NSV_ELECBEG .AND. JSV .LT. NSV_ELECEND) THEN 
      TZFIELD%CUNITS     = 'C m-3'
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','SVT',JSV
    ELSE
      TZFIELD%CUNITS     = 'm-3'
      WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3,A8)')'X_Y_Z_','SVT',JSV,' (nb ions/m3)'
    END IF
    ZWORK31(:,:,:)=XSVT(:,:,:,JSV) * XRHODREF(:,:,:)  ! C/kg --> C/m3
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END DO
END IF
! Sea Salt variables
IF (LSALT) THEN
  IF(.NOT.ALLOCATED(ZSIG_SLT)) &
    ALLOCATE(ZSIG_SLT(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3), NMODE_SLT))
  IF(.NOT.ALLOCATED(ZRG_SLT))  &
    ALLOCATE(ZRG_SLT(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3), NMODE_SLT))
  IF(.NOT.ALLOCATED(ZN0_SLT))  &
    ALLOCATE(ZN0_SLT(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3), NMODE_SLT))
  !
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CUNITS     = 'ppbv'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_SLTBEG,NSV_SLTEND
    TZFIELD%CMNHNAME   = TRIM(UPCASE(CSALTNAMES(JSV-NSV_SLTBEG+1)))//'T'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3)')'X_Y_Z_','SALT',JSV
    CALL IO_Field_write(TPFILE,TZFIELD,XSVT(:,:,:,JSV)*1.E9)
  END DO
  !
  CALL PPP2SALT(XSVT(:,:,:,NSV_SLTBEG:NSV_SLTEND),XRHODREF,&
               PSIG3D=ZSIG_SLT, PRG3D=ZRG_SLT, PN3D=ZN0_SLT)
  !
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JJ=1,NMODE_SLT
    WRITE(TZFIELD%CMNHNAME,'(A6,I1)')'SLTRGA',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'um'
    WRITE(TZFIELD%CCOMMENT,'(A18,I1)')'RG (nb) SALT MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZRG_SLT(:,:,:,JJ))
    !
    WRITE(TZFIELD%CMNHNAME,'(A7,I1)')'SLTRGAM',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'um'
    WRITE(TZFIELD%CCOMMENT,'(A17,I1)')'RG (m) SALT MODE ',JJ
    ZWORK31(:,:,:)=ZRG_SLT(:,:,:,JJ) / (EXP(-3.*(LOG(ZSIG_SLT(:,:,:,JJ)))**2))
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    !
    WRITE(TZFIELD%CMNHNAME,'(A6,I1)')'SLTN0A',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm-3'
    WRITE(TZFIELD%CCOMMENT,'(A13,I1)')'N0 SALT MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZN0_SLT(:,:,:,JJ))
    !
    WRITE(TZFIELD%CMNHNAME,'(A7,I1)')'SLTSIGA',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = '1'
    WRITE(TZFIELD%CCOMMENT,'(A16,I1)')'SIGMA SALT MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZSIG_SLT(:,:,:,JJ))
    !SALT MASS CONCENTRATION
    WRITE(TZFIELD%CMNHNAME,'(A4,I1)')'SLTMSS',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'ug m-3'
    WRITE(TZFIELD%CCOMMENT,'(A14,I1)')'MASSCONC MODE ',JJ
    ZWORK31(:,:,:)= ZN0_SLT(:,:,:,JJ)*4./3.*3.14*2500.*1e9 & !kg-->ug
       * (ZRG_SLT(:,:,:,JJ)**3)*1.d-18 &  !um-->m
       * exp(4.5*log(ZSIG_SLT(:,:,:,JJ))*log(ZSIG_SLT(:,:,:,JJ)))
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    !SALT BURDEN (g/m2)
    ZWORK21(:,:)=0.0
    DO JK=IKB,IKE
      ZWORK31(:,:,JK) = ZWORK31(:,:,JK) *(XZZ(:,:,JK+1)-XZZ(:,:,JK))      &
                       *1.d-6 ! Convert to ug/m2-->g/m2 in each layer
    END DO
    DO JK=IKB,IKE
      DO JT=IJB,IJE
        DO JI=IIB,IIE
           ZWORK21(JI,JT)=ZWORK21(JI,JT)+ZWORK31(JI,JT,JK)
        ENDDO
      ENDDO
    ENDDO
    WRITE(TZFIELD%CMNHNAME,'(A7,I1)')'SLTBRDN',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'g m-2'
    WRITE(TZFIELD%CCOMMENT,'(A6,I1)')'BURDEN',JJ
    TZFIELD%NDIMS      = 2
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
    !
    TZFIELD%NDIMS      = 3
  ENDDO
END IF
IF (LSALT.AND.LDEPOS_SLT(IMI)) THEN
  !
  ZSSLTDEP=XSVT(:,:,:,NSV_SLTDEPBEG:NSV_SLTDEPEND)
  !
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppbv'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = 1,NSV_SLTDEP
    TZFIELD%CMNHNAME   = TRIM(UPCASE(CDESLTNAMES(JSV)))//'T'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3)')'X_Y_Z_','SALTDEP',JSV
    CALL IO_Field_write(TPFILE,TZFIELD,ZSSLTDEP(:,:,:,JSV)*1.E9)
  END DO
  !
  DO JJ=1,NMODE_SLT
    ! FOR CLOUDS
    WRITE(TZFIELD%CMNHNAME,'(A9,I1)')'SLTDEPN0A',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A16,I1)')'N0 DUSTDEP MODE ',JJ
    TZFIELD%CUNITS     = 'm-3'
    ! CLOUD: CALCULATE MOMENT 3 FROM TOTAL AEROSOL MASS
    ZWORK31(:,:,:) = ZSSLTDEP(:,:,:,JJ)         &!==>molec_{aer}/molec_{air}
                     *(XMOLARWEIGHT_DUST/XMD)   &!==>kg_{aer}/kg_{air}
                     *XRHODREF(:,:,:)           &!==>kg_{aer}/m3_{air}
                     /XDENSITY_DUST             &!==>m3_{aer}/m3_{air}
                     *XM3TOUM3                  &!==>um3_{aer}/m3_{air}
                     /(XPI*4./3.)                !==>um3_{aer}/m3_{air}
            !==>volume 3rd moment
    !CLOUD: CALCULATE MOMENT 0 FROM DISPERSION AND MEAN RADIUS
    ZWORK31(:,:,:) =  ZWORK31(:,:,:)/                        &
                      ((ZRG_SLT(:,:,:,JJ)**3)*               &
                      EXP(4.5 * LOG(ZSIG_SLT(:,:,:,JJ))**2))
    !CLOUD: RETURN TO CONCENTRATION #/m3
    ZWORK31(:,:,:)= ZWORK31(:,:,:) *   XMD/ &
                    (XAVOGADRO*XRHODREF(:,:,:))
    !CLOUD:  Get number concentration (#/molec_{air}==>#/m3)
    ZWORK31(:,:,:)=                   &
                    ZWORK31(:,:,:)    & !#/molec_{air}
                    * XAVOGADRO       & !==>#/mole
                    / XMD             & !==>#/kg_{air}
                    * XRHODREF(:,:,:)   !==>#/m3
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    ! CLOUD:   DUST MASS CONCENTRATION
    WRITE(TZFIELD%CMNHNAME,'(A9,I1)')'SLTDEPMSS',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A17,I1)')'DEPMASSCONC MODE ',JJ
    TZFIELD%CUNITS     = 'ug m-3'
    ZWORK31(:,:,:)= ZWORK31(:,:,:)*4./3.*3.14*2500.*1e9 & !kg-->ug
                    * (ZRG_SLT(:,:,:,JJ)**3)*1.d-18     &  !um-->m
                    * exp(4.5*log(ZSIG_SLT(:,:,:,JJ))*log(ZSIG_SLT(:,:,:,JJ)))
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    !   FOR RAIN DROPS
    WRITE(TZFIELD%CMNHNAME,'(A9,I1)')'SLTDEPN0A',JJ+NMODE_SLT
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A16,I1)')'N0 DUSTDEP MODE ',JJ+NMODE_SLT
    TZFIELD%CUNITS     = 'm-3'
    ! RAIN: CALCULATE MOMENT 3 FROM TOTAL AEROSOL MASS
    ZWORK31(:,:,:) = ZSSLTDEP(:,:,:,JJ+NMODE_SLT)  &!==>molec_{aer}/molec_{air}
                     *(XMOLARWEIGHT_DUST/XMD)      &!==>kg_{aer}/kg_{air}
                     *XRHODREF(:,:,:)              &!==>kg_{aer}/m3_{air}
                     /XDENSITY_DUST                &!==>m3_{aer}/m3_{air}
                     *XM3TOUM3                     &!==>um3_{aer}/m3_{air}
                     /(XPI*4./3.)                   !==>um3_{aer}/m3_{air}
            !==>volume 3rd moment   
    !RAIN: CALCULATE MOMENT 0 FROM DISPERSION AND MEAN RADIUS
    ZWORK31(:,:,:)= ZWORK31(:,:,:)/                        &
                    ((ZRG_SLT(:,:,:,JJ)**3)*               &
                    EXP(4.5 * LOG(ZSIG_SLT(:,:,:,JJ))**2))
    !RAIN: RETURN TO CONCENTRATION #/m3
    ZWORK31(:,:,:)= ZWORK31(:,:,:) *   XMD/ &
                    (XAVOGADRO*XRHODREF(:,:,:))
    !RAIN: Get number concentration (#/molec_{air}==>#/m3)
    ZWORK31(:,:,:)=                   &
                    ZWORK31(:,:,:)    & !#/molec_{air}
                    * XAVOGADRO       & !==>#/mole
                    / XMD             & !==>#/kg_{air}
                    * XRHODREF(:,:,:)   !==>#/m3
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    ! RAIN:   DUST MASS CONCENTRATION
    WRITE(TZFIELD%CMNHNAME,'(A9,I1)')'SLTDEPMSS',JJ+NMODE_SLT
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A17,I1)')'DEPMASSCONC MODE ',JJ+NMODE_SLT
    TZFIELD%CUNITS     = 'ug m-3'
    ZWORK31(:,:,:)= ZWORK31(:,:,:)*4./3.*3.14*2500.*1e9 & !kg-->ug
                    * (ZRG_SLT(:,:,:,JJ)**3)*1.d-18     &  !um-->m
                    * exp(4.5*log(ZSIG_SLT(:,:,:,JJ))*log(ZSIG_SLT(:,:,:,JJ)))
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END DO
!
END IF
! Dust variables
IF (LDUST) THEN
  IF(.NOT.ALLOCATED(ZSIG_DST)) &
    ALLOCATE(ZSIG_DST(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3), NMODE_DST))
  IF(.NOT.ALLOCATED(ZRG_DST))  &
    ALLOCATE(ZRG_DST(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3), NMODE_DST))
  IF(.NOT.ALLOCATED(ZN0_DST))  &
    ALLOCATE(ZN0_DST(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3), NMODE_DST))
  !
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppbv'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_DSTBEG,NSV_DSTEND
    TZFIELD%CMNHNAME   = TRIM(UPCASE(CDUSTNAMES(JSV-NSV_DSTBEG+1)))//'T'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3)')'X_Y_Z_','DUST',JSV
    CALL IO_Field_write(TPFILE,TZFIELD,XSVT(:,:,:,JSV)*1.E9)
  END DO
  !
  CALL PPP2DUST(XSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND),XRHODREF,&
               PSIG3D=ZSIG_DST, PRG3D=ZRG_DST, PN3D=ZN0_DST)
  DO JJ=1,NMODE_DST
    WRITE(TZFIELD%CMNHNAME,'(A6,I1)')'DSTRGA',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'um'
    WRITE(TZFIELD%CCOMMENT,'(A18,I1)')'RG (nb) DUST MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZRG_DST(:,:,:,JJ))
    !
    WRITE(TZFIELD%CMNHNAME,'(A7,I1)')'DSTRGAM',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'um'
    WRITE(TZFIELD%CCOMMENT,'(A17,I1)')'RG (m) DUST MODE ',JJ
    ZWORK31(:,:,:)=ZRG_DST(:,:,:,JJ) / (EXP(-3.*(LOG(ZSIG_DST(:,:,:,JJ)))**2))
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    !
    WRITE(TZFIELD%CMNHNAME,'(A6,I1)')'DSTN0A',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm-3'
    WRITE(TZFIELD%CCOMMENT,'(A13,I1)')'N0 DUST MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZN0_DST(:,:,:,JJ))
    !
    WRITE(TZFIELD%CMNHNAME,'(A7,I1)')'DSTSIGA',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = '1'
    WRITE(TZFIELD%CCOMMENT,'(A16,I1)')'SIGMA DUST MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZSIG_DST(:,:,:,JJ))
    !DUST MASS CONCENTRATION
    WRITE(TZFIELD%CMNHNAME,'(A4,I1)')'DSTMSS',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'ug m-3'
    WRITE(TZFIELD%CCOMMENT,'(A14,I1)')'MASSCONC MODE ',JJ
    ZWORK31(:,:,:)= ZN0_DST(:,:,:,JJ)*4./3.*3.14*2500.*1e9 & !kg-->ug
       * (ZRG_DST(:,:,:,JJ)**3)*1.d-18 &  !um-->m
       * exp(4.5*log(ZSIG_DST(:,:,:,JJ))*log(ZSIG_DST(:,:,:,JJ)))
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    !DUST BURDEN (g/m2)
    ZWORK21(:,:)=0.0
    DO JK=IKB,IKE
      ZWORK31(:,:,JK) = ZWORK31(:,:,JK) *(XZZ(:,:,JK+1)-XZZ(:,:,JK))      &
                       *1.d-6 ! Convert to ug/m2-->g/m2 in each layer
    END DO
    DO JK=IKB,IKE
      DO JT=IJB,IJE
        DO JI=IIB,IIE
           ZWORK21(JI,JT)=ZWORK21(JI,JT)+ZWORK31(JI,JT,JK)
        ENDDO
      ENDDO
    ENDDO
    WRITE(TZFIELD%CMNHNAME,'(A7,I1)')'DSTBRDN',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'g m-2'
    WRITE(TZFIELD%CCOMMENT,'(A6,I1)')'BURDEN',JJ
    TZFIELD%NDIMS      = 2
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
    !
    TZFIELD%NDIMS      = 3
  ENDDO
END IF
IF (LDUST.AND.LDEPOS_DST(IMI)) THEN
  !
  ZSDSTDEP=XSVT(:,:,:,NSV_DSTDEPBEG:NSV_DSTDEPEND)
  !
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppbv'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = 1,NSV_DSTDEP
    TZFIELD%CMNHNAME   = TRIM(UPCASE(CDEDSTNAMES(JSV)))//'T'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3)')'X_Y_Z_','DUSTDEP',JSV
    CALL IO_Field_write(TPFILE,TZFIELD,ZSDSTDEP(:,:,:,JSV)*1.E9)
  END DO
  !
  DO JJ=1,NMODE_DST
    ! FOR CLOUDS
    WRITE(TZFIELD%CMNHNAME,'(A9,I1)')'DSTDEPN0A',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A16,I1)')'N0 DUSTDEP MODE ',JJ
    TZFIELD%CUNITS     = 'm-3'
    ! CLOUD: CALCULATE MOMENT 3 FROM TOTAL AEROSOL MASS
    ZWORK31(:,:,:) = ZSDSTDEP(:,:,:,JJ)         &!==>molec_{aer}/molec_{air}
                     *(XMOLARWEIGHT_DUST/XMD)   &!==>kg_{aer}/kg_{air}
                     *XRHODREF(:,:,:)           &!==>kg_{aer}/m3_{air}
                     /XDENSITY_DUST             &!==>m3_{aer}/m3_{air}
                     *XM3TOUM3                  &!==>um3_{aer}/m3_{air}
                     /(XPI*4./3.)                !==>um3_{aer}/m3_{air}
            !==>volume 3rd moment
    !CLOUD: CALCULATE MOMENT 0 FROM DISPERSION AND MEAN RADIUS
    ZWORK31(:,:,:)=  ZWORK31(:,:,:)/      &
                    ((ZRG_DST(:,:,:,JJ)**3)*      &
                    EXP(4.5 * LOG(ZSIG_DST(:,:,:,JJ))**2))
    !CLOUD: RETURN TO CONCENTRATION #/m3
    ZWORK31(:,:,:)= ZWORK31(:,:,:) *   XMD/ &
                     (XAVOGADRO*XRHODREF(:,:,:))
    !CLOUD:  Get number concentration (#/molec_{air}==>#/m3)
    ZWORK31(:,:,:)=                         &
                    ZWORK31(:,:,:)                  & !#/molec_{air}
                    * XAVOGADRO                     & !==>#/mole
                    / XMD                           & !==>#/kg_{air}
                    * XRHODREF(:,:,:)                 !==>#/m3  
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    ! CLOUD:   DUST MASS CONCENTRATION
    WRITE(TZFIELD%CMNHNAME,'(A9,I1)')'DSTDEPMSS',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A17,I1)')'DEPMASSCONC MODE ',JJ
    TZFIELD%CUNITS     = 'ug m-3'
    ZWORK31(:,:,:)= ZWORK31(:,:,:)*4./3.*3.14*2500.*1e9 & !kg-->ug
          * (ZRG_DST(:,:,:,JJ)**3)*1.d-18               &  !um-->m
          * exp(4.5*log(ZSIG_DST(:,:,:,JJ))*log(ZSIG_DST(:,:,:,JJ)))
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    !   FOR RAIN DROPS
    WRITE(TZFIELD%CMNHNAME,'(A9,I1)')'DSTDEPN0A',JJ+NMODE_DST
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A16,I1)')'N0 DUSTDEP MODE ',JJ+NMODE_DST
    TZFIELD%CUNITS     = 'm-3'
    ! RAIN: CALCULATE MOMENT 3 FROM TOTAL AEROSOL MASS
    ZWORK31(:,:,:)=ZSDSTDEP(:,:,:,JJ+NMODE_DST)  &!==>molec_{aer}/molec_{air}
            *(XMOLARWEIGHT_DUST/XMD)             &!==>kg_{aer}/kg_{air}
            *XRHODREF(:,:,:)                     &!==>kg_{aer}/m3_{air}
            *(1.d0/XDENSITY_DUST)                &!==>m3_{aer}/m3_{air}
            *XM3TOUM3                            &!==>um3_{aer}/m3_{air}
            /(XPI*4./3.)                          !==>um3_{aer}/m3_{air}
            !==>volume 3rd moment   
    !RAIN: CALCULATE MOMENT 0 FROM DISPERSION AND MEAN RADIUS
    ZWORK31(:,:,:)= ZWORK31(:,:,:)/                 &
             ((ZRG_DST(:,:,:,JJ)**3)*               &
              EXP(4.5 * LOG(ZSIG_DST(:,:,:,JJ))**2))
    !RAIN: RETURN TO CONCENTRATION #/m3
    ZWORK31(:,:,:)= ZWORK31(:,:,:) *   XMD/ &
                    (XAVOGADRO*XRHODREF(:,:,:))
    !RAIN: Get number concentration (#/molec_{air}==>#/m3)
    ZWORK31(:,:,:)=                   &
                    ZWORK31(:,:,:)    & !#/molec_{air}
                    * XAVOGADRO       & !==>#/mole
                    / XMD             & !==>#/kg_{air}
                    * XRHODREF(:,:,:)   !==>#/m3
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    ! RAIN:   DUST MASS CONCENTRATION
    WRITE(TZFIELD%CMNHNAME,'(A9,I1)')'DSTDEPMSS',JJ+NMODE_DST
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A17,I1)')'DEPMASSCONC MODE ',JJ+NMODE_DST
    TZFIELD%CUNITS     = 'ug m-3'
    ZWORK31(:,:,:)= ZWORK31(:,:,:)*4./3.*3.14*2500.*1e9 & !kg-->ug
                    * (ZRG_DST(:,:,:,JJ)**3)*1.d-18     &  !um-->m
                    * exp(4.5*log(ZSIG_DST(:,:,:,JJ))*log(ZSIG_DST(:,:,:,JJ)))
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END DO
!
END IF
! Aerosol
IF ((LCHEMDIAG).AND.(LORILAM).AND.(LUSECHEM)) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppbv'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_AERBEG,NSV_AEREND
    TZFIELD%CMNHNAME   = TRIM(UPCASE(CAERONAMES(JSV-NSV_AERBEG+1)))//'T'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A6,A4,I3.3)')'X_Y_Z_','AERO',JSV
    CALL IO_Field_write(TPFILE,TZFIELD,XSVT(:,:,:,JSV)*1.E9)
  END DO
  !
  IF (.NOT.(ASSOCIATED(XN3D)))   &
    ALLOCATE(XN3D(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),JPMODE))
  IF (.NOT.(ASSOCIATED(XRG3D)))  &
    ALLOCATE(XRG3D(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),JPMODE))
  IF (.NOT.(ASSOCIATED(XSIG3D))) &
    ALLOCATE(XSIG3D(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),JPMODE))
  !
  CALL  PPP2AERO(XSVT(:,:,:,NSV_AERBEG:NSV_AEREND), XRHODREF, &
                 PSIG3D=XSIG3D, PRG3D=XRG3D, PN3D=XN3D, PCTOTA=ZPTOTA) 
  DO JJ=1,JPMODE
    WRITE(TZFIELD%CMNHNAME,'(A3,I1)')'RGA',JJ
    TZFIELD%CLONGNAME  = 'RGA'
    TZFIELD%CUNITS     = 'um'
    WRITE(TZFIELD%CCOMMENT,'(A21,I1)')'RG (nb) AEROSOL MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,XRG3D(:,:,:,JJ))
    !
    WRITE(TZFIELD%CMNHNAME,'(A4,I1)')'RGAM',JJ
    TZFIELD%CLONGNAME  = 'RGAM'
    TZFIELD%CUNITS     = 'um'
    WRITE(TZFIELD%CCOMMENT,'(A20,I1)')'RG (m) AEROSOL MODE ',JJ
    ZWORK31(:,:,:)=XRG3D(:,:,:,JJ) / (EXP(-3.*(LOG(XSIG3D(:,:,:,JJ)))**2))
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    !
    WRITE(TZFIELD%CMNHNAME,'(A3,I1)')'N0A',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'cm-3'
    WRITE(TZFIELD%CCOMMENT,'(A16,I1)')'N0 AEROSOL MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,XN3D(:,:,:,JJ)*1.E-6)
    !
    WRITE(TZFIELD%CMNHNAME,'(A4,I1)')'SIGA',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = '1'
    WRITE(TZFIELD%CCOMMENT,'(A19,I1)')'SIGMA AEROSOL MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,XSIG3D(:,:,:,JJ))
    !
    WRITE(TZFIELD%CMNHNAME,'(A4,I1)')'MSO4',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'ug m-3'
    WRITE(TZFIELD%CCOMMENT,'(A22,I1)')'MASS SO4 AEROSOL MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_SO4,JJ))
    !
    WRITE(TZFIELD%CMNHNAME,'(A4,I1)')'MNO3',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'ug m-3'
    WRITE(TZFIELD%CCOMMENT,'(A22,I1)')'MASS NO3 AEROSOL MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_NO3,JJ))
    !
    WRITE(TZFIELD%CMNHNAME,'(A4,I1)')'MNH3',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'ug m-3'
    WRITE(TZFIELD%CCOMMENT,'(A22,I1)')'MASS NH3 AEROSOL MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_NH3,JJ))
    !
    WRITE(TZFIELD%CMNHNAME,'(A4,I1)')'MH2O',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'ug m-3'
    WRITE(TZFIELD%CCOMMENT,'(A22,I1)')'MASS H2O AEROSOL MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_H2O,JJ))
    !
    IF (NSOA .EQ. 10) THEN
      WRITE(TZFIELD%CMNHNAME,'(A5,I1)')'MSOA1',JJ
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'ug m-3'
      WRITE(TZFIELD%CCOMMENT,'(A23,I1)')'MASS SOA1 AEROSOL MODE ',JJ
      CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_SOA1,JJ))
      !
      WRITE(TZFIELD%CMNHNAME,'(A5,I1)')'MSOA2',JJ
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'ug m-3'
      WRITE(TZFIELD%CCOMMENT,'(A23,I1)')'MASS SOA2 AEROSOL MODE ',JJ
      CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_SOA2,JJ))
      !
      WRITE(TZFIELD%CMNHNAME,'(A5,I1)')'MSOA3',JJ
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'ug m-3'
      WRITE(TZFIELD%CCOMMENT,'(A23,I1)')'MASS SOA3 AEROSOL MODE ',JJ
      CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_SOA3,JJ))
      !
      WRITE(TZFIELD%CMNHNAME,'(A5,I1)')'MSOA4',JJ
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'ug m-3'
      WRITE(TZFIELD%CCOMMENT,'(A23,I1)')'MASS SOA4 AEROSOL MODE ',JJ
      CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_SOA4,JJ))
      !
      WRITE(TZFIELD%CMNHNAME,'(A5,I1)')'MSOA5',JJ
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'ug m-3'
      WRITE(TZFIELD%CCOMMENT,'(A23,I1)')'MASS SOA5 AEROSOL MODE ',JJ
      CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_SOA5,JJ))
      !
      WRITE(TZFIELD%CMNHNAME,'(A5,I1)')'MSOA6',JJ
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'ug m-3'
      WRITE(TZFIELD%CCOMMENT,'(A23,I1)')'MASS SOA6 AEROSOL MODE ',JJ
      CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_SOA6,JJ))
      !
      WRITE(TZFIELD%CMNHNAME,'(A5,I1)')'MSOA7',JJ
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'ug m-3'
      WRITE(TZFIELD%CCOMMENT,'(A23,I1)')'MASS SOA7 AEROSOL MODE ',JJ
      CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_SOA7,JJ))
      !
      WRITE(TZFIELD%CMNHNAME,'(A5,I1)')'MSOA8',JJ
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'ug m-3'
      WRITE(TZFIELD%CCOMMENT,'(A23,I1)')'MASS SOA8 AEROSOL MODE ',JJ
      CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_SOA8,JJ))
      !
      WRITE(TZFIELD%CMNHNAME,'(A5,I1)')'MSOA9',JJ
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'ug m-3'
      WRITE(TZFIELD%CCOMMENT,'(A23,I1)')'MASS SOA9 AEROSOL MODE ',JJ
      CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_SOA9,JJ))
      !
      WRITE(TZFIELD%CMNHNAME,'(A6,I1)')'MSOA10',JJ
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      TZFIELD%CUNITS     = 'ug m-3'
      WRITE(TZFIELD%CCOMMENT,'(A24,I1)')'MASS SOA10 AEROSOL MODE ',JJ
      CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_SOA10,JJ))
    END IF
    !
    WRITE(TZFIELD%CMNHNAME,'(A3,I1)')'MOC',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'ug m-3'
    WRITE(TZFIELD%CCOMMENT,'(A21,I1)')'MASS OC AEROSOL MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_OC,JJ))
    !
    WRITE(TZFIELD%CMNHNAME,'(A3,I1)')'MBC',JJ
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'ug m-3'
    WRITE(TZFIELD%CCOMMENT,'(A21,I1)')'MASS BC AEROSOL MODE ',JJ
    CALL IO_Field_write(TPFILE,TZFIELD,ZPTOTA(:,:,:,JP_AER_BC,JJ))
  ENDDO
END IF
!
!* Large Scale variables
!
IF (LVAR_LS) THEN
  CALL IO_Field_write(TPFILE,'LSUM', XLSUM)
  CALL IO_Field_write(TPFILE,'LSVM', XLSVM)
  !
  IF (LWIND_ZM) THEN
    TZFIELD2(1)%CMNHNAME   = 'LSUM_ZM'
    TZFIELD2(1)%CSTDNAME   = ''
    TZFIELD2(1)%CLONGNAME  = 'LSUM_ZM'
    TZFIELD2(1)%CUNITS     = 'm s-1'
    TZFIELD2(1)%CDIR       = 'XY'
    TZFIELD2(1)%CCOMMENT   = 'Large Scale Zonal component of horizontal wind'
    TZFIELD2(1)%NGRID      = 2
    TZFIELD2(1)%NTYPE      = TYPEREAL
    TZFIELD2(1)%NDIMS      = 3
    TZFIELD2(1)%LTIMEDEP   = .TRUE.
    !
    TZFIELD2(2)%CMNHNAME   = 'LSVM_ZM'
    TZFIELD2(2)%CSTDNAME   = ''
    TZFIELD2(2)%CLONGNAME  = 'LSVM_ZM'
    TZFIELD2(2)%CUNITS     = 'm s-1'
    TZFIELD2(2)%CDIR       = 'XY'
    TZFIELD2(2)%CCOMMENT   = 'Large Scale Meridian component of horizontal wind'
    TZFIELD2(2)%NGRID      = 3
    TZFIELD2(2)%NTYPE      = TYPEREAL
    TZFIELD2(2)%NDIMS      = 3
    TZFIELD2(2)%LTIMEDEP   = .TRUE.
    !
    CALL UV_TO_ZONAL_AND_MERID(XLSUM,XLSVM,23,TPFILE=TPFILE,TZFIELDS=TZFIELD2)
  ENDIF
  !
  CALL IO_Field_write(TPFILE,'LSWM', XLSWM)
  CALL IO_Field_write(TPFILE,'LSTHM',XLSTHM)
!
  IF (LUSERV) THEN
    CALL FIND_FIELD_ID_FROM_MNHNAME('LSRVM',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CUNITS = 'g kg-1'
    CALL IO_Field_write(TPFILE,TZFIELD,XLSRVM(:,:,:)*1.E3)
  END IF
END IF
!
!* Forcing variables
!
IF (LVAR_FRC .AND. LFORCING) THEN
!
  DO JT=1,NFRC
    WRITE (YFRC,'(I3.3)') JT
!
    TZFIELD%CMNHNAME   = 'UFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Zonal component of horizontal forcing wind'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_Field_write(TPFILE,TZFIELD,XUFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'VFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Meridian component of horizontal forcing wind'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_Field_write(TPFILE,TZFIELD,XVFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'WFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Vertical forcing wind'
    TZFIELD%NGRID      = 4
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_Field_write(TPFILE,TZFIELD,XWFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'THFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Forcing potential temperature'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_Field_write(TPFILE,TZFIELD,XTHFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'RVFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'kg kg-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Forcing vapor mixing ratio'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_Field_write(TPFILE,TZFIELD,XRVFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'TENDTHFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Large-scale potential temperature tendency for forcing'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_Field_write(TPFILE,TZFIELD,XTENDTHFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'TENDRVFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'kg kg-1 s-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Large-scale vapor mixing ratio tendency for forcing'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_Field_write(TPFILE,TZFIELD,XTENDRVFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'GXTHFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K m-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Large-scale potential temperature gradient for forcing'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_Field_write(TPFILE,TZFIELD,XGXTHFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'GYTHFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K m-1'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Large-scale potential temperature gradient for forcing'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 1
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_Field_write(TPFILE,TZFIELD,XGYTHFRC(:,JT))
!
    TZFIELD%CMNHNAME   = 'PGROUNDFRC'//YFRC
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'Pa'
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = 'Forcing ground pressure'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 0
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_Field_write(TPFILE,TZFIELD,XPGROUNDFRC(JT))
!
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
!*       1.7    Some diagnostic variables
!
IF (LTPZH .OR. LCOREF) THEN
!
!* Temperature in celsius
  TZFIELD%CMNHNAME   = 'TEMP'
  TZFIELD%CSTDNAME   = 'air_temperature'
  TZFIELD%CLONGNAME  = 'TEMP'
  TZFIELD%CUNITS     = 'celsius'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_TEMPerature'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  ZWORK31(:,:,:)=ZTEMP(:,:,:) - XTT
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!
!* Pressure in hPa        
  CALL FIND_FIELD_ID_FROM_MNHNAME('PABST',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CMNHNAME   = 'PRES'
  TZFIELD%CUNITS = 'hPa'
  CALL IO_Field_write(TPFILE,TZFIELD,XPABST(:,:,:)*1E-2)
!
!* Geopotential in meters
  CALL IO_Field_write(TPFILE,'ALT',XZZ)
!
!* Relative humidity in percent
  IF (LUSERV) THEN
    ZWORK31(:,:,:)=SM_FOES(ZTEMP(:,:,:))
    ZWORK33(:,:,:)=ZWORK31(:,:,:)
    ZWORK31(:,:,:)=(XMV/XMD)*ZWORK31(:,:,:)/(XPABST(:,:,:)-ZWORK31(:,:,:))
    ZWORK32(:,:,:)=100.*XRT(:,:,:,1)/ZWORK31(:,:,:)
    IF (CCLOUD(1:3) =='ICE' .OR. CCLOUD =='C3R5' .OR. CCLOUD == 'LIMA')  THEN
      WHERE ( ZTEMP(:,:,:)< XTT)
        ZWORK31(:,:,:) = EXP( XALPI - XBETAI/ZTEMP(:,:,:) &
                       - XGAMI*ALOG(ZTEMP(:,:,:)) ) !saturation over ice
        ZWORK33(:,:,:)=ZWORK31(:,:,:)
        ZWORK31(:,:,:)=(XMV/XMD)*ZWORK31(:,:,:)/(XPABST(:,:,:)-ZWORK31(:,:,:))
        ZWORK32(:,:,:)=100.*XRT(:,:,:,1)/ZWORK31(:,:,:)
      END WHERE
    END IF
    !
    TZFIELD%CMNHNAME   = 'REHU'
    TZFIELD%CSTDNAME   = 'relative_humidity'
    TZFIELD%CLONGNAME  = 'REHU'
    TZFIELD%CUNITS     = 'percent'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_RElative HUmidity'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK32)
    !
    TZFIELD%CMNHNAME   = 'VPRES'
    TZFIELD%CSTDNAME   = 'water_vapor_partial_pressure_in_air'
    TZFIELD%CLONGNAME  = 'VPRES'
    TZFIELD%CUNITS     = 'hPa'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Vapor PRESsure'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    ZWORK33(:,:,:)=ZWORK33(:,:,:)*ZWORK32(:,:,:)*1E-4
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK33)
    !
    IF (LCOREF) THEN
      ZWORK33(:,:,:)=(77.6*( XPABST(:,:,:)*1E-2                &
                            +ZWORK33(:,:,:)*4810/ZTEMP(:,:,:)) &
                      -6*ZWORK33(:,:,:)                        )/ZTEMP(:,:,:)
      TZFIELD%CMNHNAME   = 'COREF'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = 'COREF'
      TZFIELD%CUNITS     = '1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CCOMMENT   = 'X_Y_Z_REFraction COindex (N-units)'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      CALL IO_Field_write(TPFILE,TZFIELD,ZWORK33)
      !
      ZWORK33(:,:,:)=ZWORK33(:,:,:)+MZF(XZZ(:,:,:))*1E6/XRADIUS
      TZFIELD%CMNHNAME   = 'MCOREF'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = 'MCOREF'
      TZFIELD%CUNITS     = '1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CCOMMENT   = 'X_Y_Z_Modified REFraction COindex (M-units)'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      CALL IO_Field_write(TPFILE,TZFIELD,ZWORK33)
    END IF
  ELSE
    PRINT*, 'NO WATER VAPOR IN ',TPFILE%CNAME,' RELATIVE HUMIDITY IS NOT COMPUTED'
  END IF
!
END IF
!
!-------------------------------------------------------------------------------
!
!* Virtual potential temperature
!
IF ( LMOIST_V .OR. LMSLP .OR. CBLTOP/='NONE' ) THEN
  ALLOCATE(ZTHETAV(IIU,IJU,IKU))
!
  IF(NRR > 0) THEN
!   compute the ratio : 1 + total water mass / dry air mass
    ZRV_OV_RD = XRV / XRD
    ZTHETAV(:,:,:) = 1. + XRT(:,:,:,1)
    DO JLOOP = 2,1+NRRL+NRRI                
      ZTHETAV(:,:,:) = ZTHETAV(:,:,:) + XRT(:,:,:,JLOOP)
    END DO
! compute the virtual potential temperature when water is present in any form
    ZTHETAV(:,:,:) = XTHT(:,:,:) * (1.+XRT(:,:,:,1)*ZRV_OV_RD) / ZTHETAV(:,:,:)
  ELSE
! compute the virtual potential temperature when water is absent
    ZTHETAV(:,:,:) = XTHT(:,:,:)
  END IF
!
  IF (LMOIST_V .AND. NRR > 0) THEN
! Virtual potential temperature
    TZFIELD%CMNHNAME   = 'THETAV'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THETAV'
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Virtual potential temperature'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZTHETAV)
  END IF
!
END IF
!
!-------------------------------------------------------------------------------
!
!* Fog Visibility
!
IF (LVISI) THEN
!
  IF ((CCLOUD /= 'NONE') .AND. (CCLOUD /='REVE')) ALLOCATE(ZVISIKUN(IIU,IJU,IKU))
  IF ((CCLOUD == 'C2R2') .OR. (CCLOUD =='KHKO')) THEN
    ALLOCATE(ZVISIGUL(IIU,IJU,IKU))
    ALLOCATE(ZVISIZHA(IIU,IJU,IKU))
  END IF
!
  IF ((CCLOUD /= 'NONE') .AND. (CCLOUD /='REVE')) THEN
    ZVISIKUN(:,:,:) = 10000.
    WHERE ( XRT(:,:,:,2) >= 1E-08 )
     ZVISIKUN(:,:,:) =0.027/(XRT(:,:,:,2)*XRHODREF(:,:,:))**0.88*1000.
    END WHERE
! Visibity Kunkel                     
    TZFIELD%CMNHNAME   = 'VISIKUN'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'VISIKUN'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Visibility Kunkel'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZVISIKUN)
!
    IF ((CCLOUD == 'C2R2') .OR. (CCLOUD =='KHKO')) THEN
      ZVISIGUL(:,:,:) = 10000.
      ZVISIZHA(:,:,:) = 10000.
      WHERE ( (XRT(:,:,:,2) >= 1E-08 ) .AND. (XSVT(:,:,:,NSV_C2R2BEG+1) >=0.001 ) )
       ZVISIGUL(:,:,:) =1.002/(XRT(:,:,:,2)*XRHODREF(:,:,:)*XSVT(:,:,:,NSV_C2R2BEG+1))**0.6473*1000.
       ZVISIZHA(:,:,:) =0.187/(XRT(:,:,:,2)*XRHODREF(:,:,:)*XSVT(:,:,:,NSV_C2R2BEG+1))**0.34*1000.
      END WHERE
! Visibity Gultepe                    
      TZFIELD%CMNHNAME   = 'VISIGUL'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = 'VISIGUL'
      TZFIELD%CUNITS     = 'm'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CCOMMENT   = 'X_Y_Z_Visibility Gultepe'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
      CALL IO_Field_write(TPFILE,TZFIELD,ZVISIGUL)
! Visibity Zhang                      
      TZFIELD%CMNHNAME   = 'VISIZHA'
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CLONGNAME  = 'VISIZHA'
      TZFIELD%CUNITS     = 'm'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%CCOMMENT   = 'X_Y_Z_Visibility Zhang'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
      CALL IO_Field_write(TPFILE,TZFIELD,ZVISIZHA)
!
      DEALLOCATE(ZVISIGUL,ZVISIZHA)
    END IF
    DEALLOCATE(ZVISIKUN)
  END IF
!
END IF
!
!-------------------------------------------------------------------------------
!
!* Thetae computation according eq.(21), (43) of Bolton 1980 (MWR108,p 1046-1053)
!
IF (( LMOIST_E .OR. LBV_FR ) .AND. (NRR>0)) THEN
  ALLOCATE(ZTHETAE(IIU,IJU,IKU))
  !
  ZWORK31(:,:,:) = MAX(XRT(:,:,:,1),1.E-10)
  ZTHETAE(:,:,:)= (    2840./                                          &
         (3.5*ALOG(XTHT(:,:,:)*( XPABST(:,:,:)/XP00 )**(XRD/XCPD)  )   &
         - ALOG( XPABST(:,:,:)*0.01*ZWORK31(:,:,:) / ( 0.622+ZWORK31(:,:,:) ) ) &
         -4.805   )    ) + 55.
  ZTHETAE(:,:,:)= XTHT(:,:,:) * EXP( (3376. / ZTHETAE(:,:,:) - 2.54)  &
                 *ZWORK31(:,:,:) *(1. +0.81 *ZWORK31(:,:,:)) )
!
  IF (LMOIST_E) THEN
    TZFIELD%CMNHNAME   = 'THETAE'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THETAE'
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Equivalent potential temperature'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZTHETAE)
  END IF
END IF
!-------------------------------------------------------------------------------
!
!* Thetaes computation 
!
IF (LMOIST_ES .AND. (NRR>0)) THEN
  ALLOCATE(ZTHETAES(IIU,IJU,IKU))
  ZWORK31(:,:,:) = MAX(QSAT(ZTEMP(:,:,:),XPABST(:,:,:)),1.E-10)
  ZTHETAES(:,:,:)= (    2840./                                          &
       (3.5*ALOG(XTHT(:,:,:)*( XPABST(:,:,:)/XP00 )**(XRD/XCPD)  )   &
       - ALOG( XPABST(:,:,:)*0.01*ZWORK31(:,:,:) / ( 0.622+ZWORK31(:,:,:) ) ) &
       -4.805   )    ) + 55.
  ZTHETAES(:,:,:)= XTHT(:,:,:) * EXP( (3376. / ZTHETAE(:,:,:) - 2.54)  &
               *ZWORK31(:,:,:) *(1. +0.81 *ZWORK31(:,:,:)) )
  TZFIELD%CMNHNAME   = 'THETAES'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'THETAES'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Equivalent Saturated potential temperature'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZTHETAES)
ENDIF
!
!-------------------------------------------------------------------------------
!* The Liquid-Water potential temperature (Betts, 1973)
!  (also needed for THETAS1 or THETAS2)
!
IF ( LMOIST_L .OR. LMOIST_S1 .OR. LMOIST_S2 ) THEN
!
  ALLOCATE(ZTHETAL(IIU,IJU,IKU))
!
  IF(NRR > 1) THEN
!  The latent heat of Vaporization:
    ZWORK31(:,:,:) = XLVTT + (XCPV-XCL)*(ZTEMP(:,:,:)-XTT)
!  The latent heat of Sublimation:
    ZWORK32(:,:,:) = XLSTT + (XCPV-XCI)*(ZTEMP(:,:,:)-XTT)
!  The numerator in the exponential 
!  and the total water mixing ratio:
    ZTHETAL(:,:,:) = 0.0
    ZWORK33(:,:,:) = XRT(:,:,:,1)
    DO JLOOP = 2,1+NRRL              
       ZTHETAL(:,:,:) = ZTHETAL(:,:,:) + XRT(:,:,:,JLOOP)*ZWORK31(:,:,:)
       ZWORK33(:,:,:) = ZWORK33(:,:,:) + XRT(:,:,:,JLOOP)
    END DO
    DO JLOOP = 1+NRRL+1,1+NRRL+NRRI              
       ZTHETAL(:,:,:) = ZTHETAL(:,:,:) + XRT(:,:,:,JLOOP)*ZWORK32(:,:,:)
       ZWORK33(:,:,:) = ZWORK33(:,:,:) + XRT(:,:,:,JLOOP)
    END DO
!   compute the liquid-water potential temperature 
!   theta_l = theta * exp[ -(L_vap * ql + L_sub * qi) / (c_pd * T) ]
!   when water is present in any form:
    ZTHETAL(:,:,:) = XTHT(:,:,:) &
            * exp(-ZTHETAL(:,:,:)/(1.0+ZWORK33(:,:,:))/XCPD/ZTEMP(:,:,:))
  ELSE
!   compute the liquid-water potential temperature 
!    when water is absent:
    ZTHETAL(:,:,:) = XTHT(:,:,:)
  END IF
!
  IF (LMOIST_L .AND. NRR > 0) THEN
    ! Liquid-Water potential temperature
    TZFIELD%CMNHNAME   = 'THETAL'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THETAL'
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Liquid water potential temperature'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZTHETAL)
  END IF
!
END IF
!
!-------------------------------------------------------------------------------
!
!* The Moist-air Entropy potential temperature (Marquet, QJ2011, HDR2016)
!
IF ( LMOIST_S1 .OR. LMOIST_S2 ) THEN
  IF (LMOIST_S1) THEN
    ALLOCATE(ZTHETAS1(IIU,IJU,IKU))
  END IF
  IF (LMOIST_S2) THEN
    ALLOCATE(ZTHETAS2(IIU,IJU,IKU))
  END IF
!
! The total water (ZWORK31) and condensed water (ZWORK32) mixing ratios:
  ZWORK32(:,:,:) = 0.0
  IF(NRR > 0) THEN
    DO JLOOP = 2,1+NRRL+NRRI  
       ZWORK32(:,:,:) = ZWORK32(:,:,:) + XRT(:,:,:,JLOOP)
    END DO
  END IF
  ZWORK31(:,:,:) =  ZWORK32(:,:,:) + XRT(:,:,:,1)
!
  IF (LMOIST_S1) THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! thetas1 = thetal * exp[ 5.87 * qt ] ; with qt=rt/(1+rt)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ZTHETAS1(:,:,:) = ZTHETAL(:,:,:) * &
     exp( 5.87*ZWORK31(:,:,:)/(1.0+ZWORK31(:,:,:)) )
  END IF
  IF (LMOIST_S2) THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! thetas2 = thetal * exp[ (5.87-0.46*ln(rv/0.0124)-0.46*qc) * qt ]
!           where qt=rt/(1+rt)  and qc=rc/(1+rt)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ZWORK33(:,:,:) =  5.87 - 0.46 * log(MAX(XRT(:,:,:,1),1.E-10)/0.0124)
    ZTHETAS2(:,:,:) = ZTHETAL(:,:,:) * &
     exp( ZWORK33(:,:,:)*ZWORK31(:,:,:)/(1.0+ZWORK31(:,:,:)) &
                  - 0.46*ZWORK32(:,:,:)/(1.0+ZWORK31(:,:,:)) )
  END IF
  IF (LMOIST_S1) THEN
! The Moist-air Entropy potential temperature (1st order)
    TZFIELD%CMNHNAME   = 'THETAS1'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THETAS1'
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Moist air Entropy (1st order) potential temperature'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZTHETAS1)
  END IF
  IF (LMOIST_S2) THEN
! The Moist-air Entropy potential temperature (2nd order)
    TZFIELD%CMNHNAME   = 'THETAS2'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THETAS2'
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Moist air Entropy (2nd order) potential temperature'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZTHETAS2)  
  END IF
!
END IF
!
!-------------------------------------------------------------------------------
!!
!
!* Vorticity quantities
!
IF (LVORT) THEN
! Vorticity x
  ZWORK31(:,:,:)=MYF(MZF(MXM(ZVOX(:,:,:))))
  TZFIELD%CMNHNAME   = 'UM1'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'UM1'
  TZFIELD%CUNITS     = 's-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_x component of vorticity'
  TZFIELD%NGRID      = 2
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!    
! Vorticity y
  ZWORK32(:,:,:)=MZF(MXF(MYM(ZVOY(:,:,:))))
  TZFIELD%CMNHNAME   = 'VM1'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'VM1'
  TZFIELD%CUNITS     = 's-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_y component of vorticity'
  TZFIELD%NGRID      = 3
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK32)
  !
  IF (LWIND_ZM) THEN
    TZFIELD2(1)%CMNHNAME   = 'UM1_ZM'
    TZFIELD2(1)%CSTDNAME   = ''
    TZFIELD2(1)%CLONGNAME  = 'UM1_ZM'
    TZFIELD2(1)%CUNITS     = 'm s-1'
    TZFIELD2(1)%CDIR       = 'XY'
    TZFIELD2(1)%CCOMMENT   = 'Zonal component of horizontal vorticity'
    TZFIELD2(1)%NGRID      = 2
    TZFIELD2(1)%NTYPE      = TYPEREAL
    TZFIELD2(1)%NDIMS      = 3
    TZFIELD2(1)%LTIMEDEP   = .TRUE.
    !
    TZFIELD2(2)%CMNHNAME   = 'VM1_ZM'
    TZFIELD2(2)%CSTDNAME   = ''
    TZFIELD2(2)%CLONGNAME  = 'VM1_ZM'
    TZFIELD2(2)%CUNITS     = 'm s-1'
    TZFIELD2(2)%CDIR       = 'XY'
    TZFIELD2(2)%CCOMMENT   = 'Meridian component of horizontal vorticity'
    TZFIELD2(2)%NGRID      = 3
    TZFIELD2(2)%NTYPE      = TYPEREAL
    TZFIELD2(2)%NDIMS      = 3
    TZFIELD2(2)%LTIMEDEP   = .TRUE.
    !
    CALL UV_TO_ZONAL_AND_MERID(ZWORK31,ZWORK32,23,TPFILE=TPFILE,TZFIELDS=TZFIELD2)
  ENDIF
!    
! Vorticity z
  ZWORK31(:,:,:)=MXF(MYF(MZM(ZVOZ(:,:,:))))
  TZFIELD%CMNHNAME   = 'WM1'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'WM1'
  TZFIELD%CUNITS     = 's-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_relative vorticity'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!
! Absolute Vorticity 
  ZWORK31(:,:,:)=MYF(MXF(ZVOZ(:,:,:))) + ZCORIOZ(:,:,:)
  TZFIELD%CMNHNAME   = 'ABVOR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ABVOR'
  TZFIELD%CUNITS     = 's-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_z ABsolute VORticity'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!
END IF
!    
IF ( LMEAN_POVO ) THEN 
  !
  ALLOCATE(IWORK1(SIZE(XTHT,1),SIZE(XTHT,2)))
  !
  IWORK1(:,:)=0
  ZWORK21(:,:)=0.
  IF (XMEAN_POVO(1)>XMEAN_POVO(2)) THEN
    !Invert values (smallest must be first)
    ZX0D = XMEAN_POVO(1)
    XMEAN_POVO(1) = XMEAN_POVO(2)
    XMEAN_POVO(2) = ZX0D
  END IF
  DO JK=IKB,IKE
    WHERE((XPABST(:,:,JK)>XMEAN_POVO(1)).AND.(XPABST(:,:,JK)<XMEAN_POVO(2)))
      ZWORK21(:,:)=ZWORK21(:,:)+ZPOVO(:,:,JK)
      IWORK1(:,:)=IWORK1(:,:)+1
    END WHERE
  END DO
  WHERE (IWORK1(:,:)>0) ZWORK21(:,:)=ZWORK21(:,:)/REAL( IWORK1(:,:) )
  TZFIELD%CMNHNAME   = 'MEAN_POVO'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'MEAN_POVO'
  TZFIELD%CUNITS     = 'PVU' ! 1 PVU = 1e-6 m^2 s^-1 K kg^-1
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_MEAN of POtential VOrticity'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
END IF
!
! Virtual Potential Vorticity in PV units
IF (LMOIST_V .AND. (NRR>0) ) THEN
  ZWORK31(:,:,:)=GX_M_M(ZTHETAV,XDXX,XDZZ,XDZX)
  ZWORK32(:,:,:)=GY_M_M(ZTHETAV,XDYY,XDZZ,XDZY)
  ZWORK33(:,:,:)=GZ_M_M(ZTHETAV,XDZZ)
  ZWORK34(:,:,:)= ZWORK31(:,:,:)*MZF(MYF(ZVOX(:,:,:)))     &
               + ZWORK32(:,:,:)*MZF(MXF(ZVOY(:,:,:)))     &
               + ZWORK33(:,:,:)*(MYF(MXF(ZVOZ(:,:,:))) + ZCORIOZ(:,:,:))
  ZWORK34(:,:,:)=ZWORK34(:,:,:)*1E6/XRHODREF(:,:,:)
  TZFIELD%CMNHNAME   = 'POVOV'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'POVOV'
  TZFIELD%CUNITS     = 'PVU' ! 1 PVU = 1e-6 m^2 s^-1 K kg^-1
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Virtual POtential VOrticity'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK34)
!
  IF (LMEAN_POVO) THEN
    IWORK1(:,:)=0
    ZWORK21(:,:)=0.
    DO JK=IKB,IKE
      WHERE((XPABST(:,:,JK)>XMEAN_POVO(1)).AND.(XPABST(:,:,JK)<XMEAN_POVO(2)))
          ZWORK21(:,:)=ZWORK21(:,:)+ZWORK34(:,:,JK)
          IWORK1(:,:)=IWORK1(:,:)+1
      END WHERE
    END DO
    WHERE(IWORK1(:,:)>0) ZWORK21(:,:)=ZWORK21(:,:)/REAL( IWORK1(:,:) )
    TZFIELD%CMNHNAME   = 'MEAN_POVOV'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'MEAN_POVOV'
    TZFIELD%CUNITS     = 'PVU' ! 1 PVU = 1e-6 m^2 s^-1 K kg^-1
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_MEAN of Virtual POtential VOrticity'
    TZFIELD%NGRID      = 4
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  END IF
END IF
!
! Equivalent Potential Vorticity in PV units
IF (LMOIST_E .AND. (NRR>0) ) THEN
!
  ZWORK31(:,:,:)=GX_M_M(ZTHETAE,XDXX,XDZZ,XDZX)
  ZWORK32(:,:,:)=GY_M_M(ZTHETAE,XDYY,XDZZ,XDZY)
  ZWORK33(:,:,:)=GZ_M_M(ZTHETAE,XDZZ)
  ZWORK34(:,:,:)= ZWORK31(:,:,:)*MZF(MYF(ZVOX(:,:,:)))     &
                + ZWORK32(:,:,:)*MZF(MXF(ZVOY(:,:,:)))     &
                + ZWORK33(:,:,:)*(MYF(MXF(ZVOZ(:,:,:))) + ZCORIOZ(:,:,:))
  ZWORK34(:,:,:)=ZWORK34(:,:,:)*1E6/XRHODREF(:,:,:)
  TZFIELD%CMNHNAME   = 'POVOE'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'POVOE'
  TZFIELD%CUNITS     = 'PVU' ! 1 PVU = 1e-6 m^2 s^-1 K kg^-1
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Equivalent POtential VOrticity'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK34)
!
  IF (LMEAN_POVO) THEN
    IWORK1(:,:)=0
    ZWORK21(:,:)=0.
    DO JK=IKB,IKE
      WHERE((XPABST(:,:,JK)>XMEAN_POVO(1)).AND.(XPABST(:,:,JK)<XMEAN_POVO(2)))
        ZWORK21(:,:)=ZWORK21(:,:)+ZWORK34(:,:,JK)
        IWORK1(:,:)=IWORK1(:,:)+1
      END WHERE
    END DO
    WHERE(IWORK1(:,:)>0) ZWORK21(:,:)=ZWORK21(:,:)/REAL( IWORK1(:,:) )
    TZFIELD%CMNHNAME   = 'MEAN_POVOE'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'MEAN_POVOE'
    TZFIELD%CUNITS     = 'PVU' ! 1 PVU = 1e-6 m^2 s^-1 K kg^-1
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_MEAN of Equivalent POtential VOrticity'
    TZFIELD%NGRID      = 4
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
    DEALLOCATE(IWORK1)
  END IF 
  !
END IF
!
! Equivalent Saturated Potential Vorticity in PV units
IF (LMOIST_ES .AND. (NRR>0) ) THEN
  ZWORK31(:,:,:)=GX_M_M(ZTHETAES,XDXX,XDZZ,XDZX)
  ZWORK32(:,:,:)=GY_M_M(ZTHETAES,XDYY,XDZZ,XDZY)
  ZWORK33(:,:,:)=GZ_M_M(ZTHETAES,XDZZ)
  ZWORK34(:,:,:)= ZWORK31(:,:,:)*MZF(MYF(ZVOX(:,:,:)))     &
                + ZWORK32(:,:,:)*MZF(MXF(ZVOY(:,:,:)))     &
                + ZWORK33(:,:,:)*(MYF(MXF(ZVOZ(:,:,:))) + ZCORIOZ(:,:,:))
  ZWORK34(:,:,:)=ZWORK34(:,:,:)*1E6/XRHODREF(:,:,:)
  TZFIELD%CMNHNAME   = 'POVOES'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'POVOES'
  TZFIELD%CUNITS     = 'PVU' ! 1 PVU = 1e-6 m^2 s^-1 K kg^-1
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Equivalent Saturated POtential VOrticity'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK34)
ENDIF
!
!
!-------------------------------------------------------------------------------
!
!* Horizontal divergence
!
IF (LDIV) THEN
!
  ZWORK31=GX_U_M(XUT,XDXX,XDZZ,XDZX) + GY_V_M(XVT,XDYY,XDZZ,XDZY)
  TZFIELD%CMNHNAME   = 'HDIV'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'HDIV'
  TZFIELD%CUNITS     = 's-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Horizontal DIVergence'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!
  IF (LUSERV) THEN
    TZFIELD%CMNHNAME   = 'HMDIV'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'HMDIV'
    TZFIELD%CUNITS     = 'kg m-3 s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Horizontal Moisture DIVergence HMDIV'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    ZWORK31=MXM(XRHODREF*XRT(:,:,:,1))*XUT
    ZWORK32=MYM(XRHODREF*XRT(:,:,:,1))*XVT
    ZWORK33=GX_U_M(ZWORK31,XDXX,XDZZ,XDZX) + GY_V_M(ZWORK32,XDYY,XDZZ,XDZY)
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK33)
  END IF
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Clustering
!
IF (LCLSTR) THEN
  GCLOUD(:,:,:)=.FALSE.
  GBOTUP=LBOTUP
  IF (CFIELD=='W') THEN
  WHERE(XWT(:,:,:).GT.XTHRES) GCLOUD(:,:,:)=.TRUE.
  END IF
  IF (CFIELD=='CLOUD') THEN
  WHERE((XRT(:,:,:,2)+XRT(:,:,:,4)+XRT(:,:,:,5)+XRT(:,:,:,6)).GT.XTHRES) GCLOUD(:,:,:)=.TRUE.
  END IF
  PRINT *,'CALL CLUSTERING COUNT(GCLOUD)=',COUNT(GCLOUD)
  CALL CLUSTERING(GBOTUP,GCLOUD,XWT,ICLUSTERID,ICLUSTERLV,ZCLDSIZE)
  PRINT *,'GOT OUT OF CLUSTERING'
  !
  TZFIELD%CMNHNAME   = 'CLUSTERID'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'CLUSTERID'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_CLUSTER (ID NUMBER)'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEINT
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ICLUSTERID)
  !
  TZFIELD%CMNHNAME   = 'CLUSTERLV'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'CLUSTERLV'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_CLUSTER (BASE OR TOP LEVEL)'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEINT
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ICLUSTERLV)
  !
  TZFIELD%CMNHNAME   = 'CLDSIZE'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'CLDSIZE'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_CLDSIZE (HOR. SECTION)'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZCLDSIZE)
END IF
!
!-------------------------------------------------------------------------------
!
!* Geostrophic and Ageostrophic wind (m/s)
!
IF (LGEO .OR. LAGEO) THEN
  ALLOCATE(ZPHI(IIU,IJU,IKU))
  IF(CEQNSYS=='MAE' .OR. CEQNSYS=='DUR') THEN
    ZPHI(:,:,:)=(XPABST(:,:,:)/XP00)**(XRD/XCPD)-XEXNREF(:,:,:)
    !
    ZPHI(1,1,:)=2*ZPHI(1,2,:)-ZPHI(1,3,:)
    ZPHI(1,IJU,:)=2*ZPHI(1,IJU-1,:)-ZPHI(1,IJU-2,:)
    ZPHI(IIU,1,:)=2*ZPHI(IIU,2,:)-ZPHI(IIU,3,:)
    ZPHI(IIU,IJU,:)=2*ZPHI(IIU,IJU-1,:)-ZPHI(IIU,IJU-2,:)
    ZWORK31(:,:,:)=-MXM(GY_M_M(ZPHI,XDYY,XDZZ,XDZY)*XCPD*XTHVREF/ZCORIOZ)
    !
    ZPHI(1,1,:)=2*ZPHI(2,1,:)-ZPHI(3,1,:)
    ZPHI(IIU,1,:)=2*ZPHI(IIU-1,1,:)-ZPHI(IIU-2,1,:)
    ZPHI(1,IJU,:)=2*ZPHI(2,IJU,:)-ZPHI(3,IJU,:)
    ZPHI(IIU,IJU,:)=2*ZPHI(IIU-1,IJU,:)-ZPHI(IIU-2,IJU,:)
    ZWORK32(:,:,:)=MYM(GX_M_M(ZPHI,XDXX,XDZZ,XDZX)*XCPD*XTHVREF/ZCORIOZ)
  !
  ELSE IF(CEQNSYS=='LHE') THEN
    ZPHI(:,:,:)= ((XPABST(:,:,:)/XP00)**(XRD/XCPD)-XEXNREF(:,:,:))   &
               * XCPD * XTHVREF(:,:,:)
    !
    ZPHI(1,1,:)=2*ZPHI(1,2,:)-ZPHI(1,3,:)
    ZPHI(1,IJU,:)=2*ZPHI(1,IJU-1,:)-ZPHI(1,IJU-2,:)
    ZPHI(IIU,1,:)=2*ZPHI(IIU,2,:)-ZPHI(IIU,3,:)
    ZPHI(IIU,IJU,:)=2*ZPHI(IIU,IJU-1,:)-ZPHI(IIU,IJU-2,:)
    ZWORK31(:,:,:)=-MXM(GY_M_M(ZPHI,XDYY,XDZZ,XDZY)/ZCORIOZ)
    !
    ZPHI(1,1,:)=2*ZPHI(2,1,:)-ZPHI(3,1,:)
    ZPHI(IIU,1,:)=2*ZPHI(IIU-1,1,:)-ZPHI(IIU-2,1,:)
    ZPHI(1,IJU,:)=2*ZPHI(2,IJU,:)-ZPHI(3,IJU,:)
    ZPHI(IIU,IJU,:)=2*ZPHI(IIU-1,IJU,:)-ZPHI(IIU-2,IJU,:)
    ZWORK32(:,:,:)=MYM(GX_M_M(ZPHI,XDXX,XDZZ,XDZX)/ZCORIOZ)
  END IF
  DEALLOCATE(ZPHI)
!
  IF (LGEO) THEN 
    TZFIELD%CMNHNAME   = 'UM88'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'UM88'
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_U component of GEOstrophic wind'
    TZFIELD%NGRID      = 2
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
! 
    TZFIELD%CMNHNAME   = 'VM88'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'VM88'
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_V component of GEOstrophic wind'
    TZFIELD%NGRID      = 3
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK32)
    !
    IF (LWIND_ZM) THEN
      TZFIELD2(1)%CMNHNAME   = 'UM88_ZM'
      TZFIELD2(1)%CSTDNAME   = ''
      TZFIELD2(1)%CLONGNAME  = 'UM88_ZM'
      TZFIELD2(1)%CUNITS     = 'm s-1'
      TZFIELD2(1)%CDIR       = 'XY'
      TZFIELD2(1)%CCOMMENT   = 'Zonal component of GEOstrophic wind'
      TZFIELD2(1)%NGRID      = 2
      TZFIELD2(1)%NTYPE      = TYPEREAL
      TZFIELD2(1)%NDIMS      = 3
      TZFIELD2(1)%LTIMEDEP   = .TRUE.
      !
      TZFIELD2(2)%CMNHNAME   = 'VM88_ZM'
      TZFIELD2(2)%CSTDNAME   = ''
      TZFIELD2(2)%CLONGNAME  = 'VM88_ZM'
      TZFIELD2(2)%CUNITS     = 'm s-1'
      TZFIELD2(2)%CDIR       = 'XY'
      TZFIELD2(2)%CCOMMENT   = 'Meridian component of GEOstrophic wind'
      TZFIELD2(2)%NGRID      = 3
      TZFIELD2(2)%NTYPE      = TYPEREAL
      TZFIELD2(2)%NDIMS      = 3
      TZFIELD2(2)%LTIMEDEP   = .TRUE.
      !
      CALL UV_TO_ZONAL_AND_MERID(ZWORK31,ZWORK32,23,TPFILE=TPFILE,TZFIELDS=TZFIELD2)
    ENDIF
!
! wm necessary to plot vertical cross sections of wind vectors
    CALL FIND_FIELD_ID_FROM_MNHNAME('WT',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CMNHNAME  = 'WM88'
    TZFIELD%CLONGNAME = 'WM88'
    CALL IO_Field_write(TPFILE,TZFIELD,XWT)
  END IF
!
  IF (LAGEO) THEN
    ZWORK31(:,:,:)=XUT(:,:,:)-ZWORK31(:,:,:)
    ZWORK32(:,:,:)=XVT(:,:,:)-ZWORK32(:,:,:)
    !
    TZFIELD%CMNHNAME   = 'UM89'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'UM89'
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_U component of AGEOstrophic wind'
    TZFIELD%NGRID      = 2
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
    !
    TZFIELD%CMNHNAME   = 'VM89'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'VM89'
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_V component of AGEOstrophic wind'
    TZFIELD%NGRID      = 3
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK32)
    !
    IF (LWIND_ZM) THEN
      TZFIELD2(1)%CMNHNAME   = 'UM89_ZM'
      TZFIELD2(1)%CSTDNAME   = ''
      TZFIELD2(1)%CLONGNAME  = 'UM89_ZM'
      TZFIELD2(1)%CUNITS     = 'm s-1'
      TZFIELD2(1)%CDIR       = 'XY'
      TZFIELD2(1)%CCOMMENT   = 'Zonal component of AGEOstrophic wind'
      TZFIELD2(1)%NGRID      = 2
      TZFIELD2(1)%NTYPE      = TYPEREAL
      TZFIELD2(1)%NDIMS      = 3
      TZFIELD2(1)%LTIMEDEP   = .TRUE.
      !
      TZFIELD2(2)%CMNHNAME   = 'VM89_ZM'
      TZFIELD2(2)%CSTDNAME   = ''
      TZFIELD2(2)%CLONGNAME  = 'VM89_ZM'
      TZFIELD2(2)%CUNITS     = 'm s-1'
      TZFIELD2(2)%CDIR       = 'XY'
      TZFIELD2(2)%CCOMMENT   = 'Meridian component of AGEOstrophic wind'
      TZFIELD2(2)%NGRID      = 3
      TZFIELD2(2)%NTYPE      = TYPEREAL
      TZFIELD2(2)%NDIMS      = 3
      TZFIELD2(2)%LTIMEDEP   = .TRUE.
      !
      CALL UV_TO_ZONAL_AND_MERID(ZWORK31,ZWORK32,23,TPFILE=TPFILE,TZFIELDS=TZFIELD2)
    ENDIF
!
! wm necessary to plot vertical cross sections of wind vectors
    CALL FIND_FIELD_ID_FROM_MNHNAME('WT',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CMNHNAME  = 'WM89'
    TZFIELD%CLONGNAME = 'WM89'
    CALL IO_Field_write(TPFILE,TZFIELD,XWT)
  END IF
!
END IF
!
!-------------------------------------------------------------------------------
!
!* Contravariant wind field
!
!
IF(LWIND_CONTRAV) THEN!$
  CALL CONTRAV ((/"TEST","TEST"/),(/"TEST","TEST"/),XUT,XVT,XWT,XDXX,XDYY,XDZZ,XDZX,XDZY, &
                ZWORK31,ZWORK32,ZWORK33,2)
  !
  TZFIELD%CMNHNAME   = 'WNORM'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'WNORM'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_W surface normal wind'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK33)
END IF
!-------------------------------------------------------------------------------
!
!* Mean Sea Level Pressure in hPa   
!
IF (LMSLP) THEN
  ZGAMREF=-6.5E-3
!  Exner function at the first mass point
  ZWORK21(:,:) = (XPABST(:,:,IKB) /XP00)**(XRD/XCPD)
!  virtual temperature at the first mass point
  ZWORK21(:,:) = ZWORK21(:,:) * ZTHETAV(:,:,IKB)
!  virtual temperature at ground level
  ZWORK21(:,:) = ZWORK21(:,:) - ZGAMREF*((XZZ(:,:,IKB)+XZZ(:,:,IKB+1))/2.-XZS(:,:))
!  virtual temperature at sea level
  ZWORK22(:,:) = ZWORK21(:,:) - ZGAMREF*XZS(:,:)
!  average underground virtual temperature
  ZWORK22(:,:) = 0.5*(ZWORK21(:,:)+ZWORK22(:,:))
!  surface pressure
  ZWORK21(:,:) = ( XPABST(:,:,IKB) + XPABST(:,:,IKB-1) )*.5
!  sea level pressure (hPa)
  ZWORK22(:,:) = 1.E-2*ZWORK21(:,:)*EXP(XG*XZS(:,:)/(XRD*ZWORK22(:,:)))
!
  TZFIELD%CMNHNAME   = 'MSLP'
  TZFIELD%CSTDNAME   = 'air_pressure_at_sea_level'
  TZFIELD%CLONGNAME  = 'MSLP'
  TZFIELD%CUNITS     = 'hPa'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Mean Sea Level Pressure'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK22)
END IF
!-------------------------------------------------------------------------------
!
!* Vapor, cloud water and ice thickness
!
IF (LTHW) THEN
!
  ZWORK21(:,:) = 0.
  IF(SIZE(XRT,4)>=1)THEN
    DO JK = IKB,IKE
      ZWORK21(:,:) = ZWORK21(:,:)+XRHODREF(:,:,JK)*XRT(:,:,JK,1) * &
                     (XZZ(:,:,JK+1)-XZZ(:,:,JK))/XRHOLW
    END DO
    ZWORK21(:,:) = ZWORK21(:,:)*1000. ! vapor water in mm unit
    TZFIELD%CMNHNAME   = 'THVW'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THVW'
    TZFIELD%CUNITS     = 'mm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_THickness of Vapor Water'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  END IF
  !
  ZWORK21(:,:) = 0.
  IF(SIZE(XRT,4)>=2)THEN
    DO JK = IKB,IKE
      ZWORK21(:,:) = ZWORK21(:,:)+XRHODREF(:,:,JK)*XRT(:,:,JK,2) * &
                     (XZZ(:,:,JK+1)-XZZ(:,:,JK))/XRHOLW
    END DO
    ZWORK21(:,:) = ZWORK21(:,:)*1000. ! cloud water in mm unit
    TZFIELD%CMNHNAME   = 'THCW'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THCW'
    TZFIELD%CUNITS     = 'mm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_THickness of Cloud Water'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  END IF
  !
  ZWORK21(:,:) = 0.
  IF(SIZE(XRT,4)>=3)THEN
    DO JK = IKB,IKE
      ZWORK21(:,:) = ZWORK21(:,:)+XRHODREF(:,:,JK)*XRT(:,:,JK,3) * &
                     (XZZ(:,:,JK+1)-XZZ(:,:,JK))/XRHOLW
    END DO
    ZWORK21(:,:) = ZWORK21(:,:)*1000. ! rain water in mm unit
    TZFIELD%CMNHNAME   = 'THRW'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THRW'
    TZFIELD%CUNITS     = 'mm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_THickness of Rain Water'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  END IF
  !
  ZWORK21(:,:)   = 0.
  IF(SIZE(XRT,4)>=4)THEN
    DO JK = IKB,IKE
      ZWORK21(:,:) = ZWORK21(:,:)+XRHODREF(:,:,JK)*XRT(:,:,JK,4) * &
                   (XZZ(:,:,JK+1)-XZZ(:,:,JK))/XRHOLW
    END DO
    ZWORK21(:,:) = ZWORK21(:,:)*1000.   ! ice thickness in mm unit
    TZFIELD%CMNHNAME   = 'THIC'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THIC'
    TZFIELD%CUNITS     = 'mm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_THickness of ICe'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  END IF
  !
  ZWORK21(:,:)   = 0.
  IF(SIZE(XRT,4)>=5)THEN
    DO JK = IKB,IKE
      ZWORK21(:,:) = ZWORK21(:,:)+XRHODREF(:,:,JK)*XRT(:,:,JK,5) * &
                   (XZZ(:,:,JK+1)-XZZ(:,:,JK))/XRHOLW
    END DO
    ZWORK21(:,:) = ZWORK21(:,:)*1000.   ! snow thickness in mm unit
    TZFIELD%CMNHNAME   = 'THSN'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THSN'
    TZFIELD%CUNITS     = 'mm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_THickness of SNow'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  END IF
  !
  ZWORK21(:,:)   = 0.
  IF(SIZE(XRT,4)>=6)THEN
    DO JK = IKB,IKE
      ZWORK21(:,:) = ZWORK21(:,:)+XRHODREF(:,:,JK)*XRT(:,:,JK,6) * &
                   (XZZ(:,:,JK+1)-XZZ(:,:,JK))/XRHOLW
    END DO
    ZWORK21(:,:) = ZWORK21(:,:)*1000.   ! graupel thickness in mm unit
    TZFIELD%CMNHNAME   = 'THGR'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THGR'
    TZFIELD%CUNITS     = 'mm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_THickness of GRaupel'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  END IF
  !
  ZWORK21(:,:)   = 0.
  IF(SIZE(XRT,4)>=7)THEN
    DO JK = IKB,IKE
      ZWORK21(:,:) = ZWORK21(:,:)+XRHODREF(:,:,JK)*XRT(:,:,JK,7) * &
                   (XZZ(:,:,JK+1)-XZZ(:,:,JK))/XRHOLW
    END DO
    ZWORK21(:,:) = ZWORK21(:,:)*1000.   ! hail thickness in mm unit
    TZFIELD%CMNHNAME   = 'THHA'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'THHA'
    TZFIELD%CUNITS     = 'mm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_THickness of HAil'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!* Accumulated and instantaneous total precip rates  in mm and mm/h
!
IF (LTOTAL_PR .AND. SIZE (XACPRR)>0 ) THEN
  ZWORK21(:,:) = 0.
  !
  IF (LUSERR) THEN
    ZWORK21(:,:) = XACPRR(:,:)*1E3      
  END IF
  IF (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'LIMA') THEN
    ZWORK21(:,:) = ZWORK21(:,:) + (XACPRS(:,:) + XACPRG(:,:))*1E3
    IF (SIZE(XINPRC) /= 0 ) &         
      ZWORK21(:,:) = ZWORK21(:,:) + XACPRC(:,:) *1E3
    IF (SIZE(XINPRH) /= 0 ) &        
      ZWORK21(:,:) = ZWORK21(:,:) + XACPRH(:,:) *1E3
  END IF
  IF (CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'KHKO' &
                                             .OR. CCLOUD == 'LIMA' ) THEN
    IF (SIZE(XINPRC) /= 0 ) &         
      ZWORK21(:,:) = ZWORK21(:,:) + XACPRC(:,:) *1E3
  END IF
  IF (CDCONV /= 'NONE') THEN
    ZWORK21(:,:) = ZWORK21(:,:) + XPACCONV(:,:)*1E3    
  END IF
  IF (LUSERR .OR. CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C3R5' .OR. &
                  CCLOUD == 'LIMA' .OR. CDCONV /= 'NONE') THEN
    TZFIELD%CMNHNAME   = 'ACTOPR'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'ACTOPR'
    TZFIELD%CUNITS     = 'mm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_ACccumulated TOtal Precipitation Rate'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  ELSE
    PRINT * ,'YOU WANT TO COMPUTE THE ACCUMULATED RAIN'
    PRINT * ,'BUT NO RAIN IS PRESENT IN THE MODEL' 
  END IF
  ! 
  ! calculation of the mean accumulated precipitations in the mesh-grid of a 
  !large-scale model
  IF (LMEAN_PR .AND. LUSERR) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'mm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Large Scale ACccumulated TOtal Precipitation Rate'
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JK=1,SIZE(XMEAN_PR),2      
      IF (XMEAN_PR(JK) .NE. XUNDEF .AND. XMEAN_PR(JK+1) .NE. XUNDEF) THEN
        PRINT * ,'MEAN accumulated RAIN: GRID ', XMEAN_PR(JK), XMEAN_PR(JK+1)
        CALL COMPUTE_MEAN_PRECIP(ZWORK21,XMEAN_PR(JK:JK+1),ZWORK22,TZFIELD%NGRID)
        !
        JI=INT(XMEAN_PR(JK))
        JJ=INT(XMEAN_PR(JK+1))
        WRITE(TZFIELD%CMNHNAME,'(A9,2I2.2)')'LS_ACTOPR',JI,JJ
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        CALL IO_Field_write(TPFILE,TZFIELD,ZWORK22)
      END IF
    END DO
    !
  END IF
  !
  !
  ZWORK21(:,:) = 0.
  !    
  IF (LUSERR) THEN
    ZWORK21(:,:) = XINPRR(:,:)*3.6E6      
  END IF
  IF (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'LIMA') THEN
    ZWORK21(:,:) = ZWORK21(:,:) + (XINPRS(:,:) + XINPRG(:,:))*3.6E6
    IF (SIZE(XINPRC) /= 0 ) &      
      ZWORK21(:,:) = ZWORK21(:,:) + XINPRC(:,:) *3.6E6       
    IF (SIZE(XINPRH) /= 0 ) &      
      ZWORK21(:,:) = ZWORK21(:,:) + XINPRH(:,:) *3.6E6       
  END IF
  IF (CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'KHKO' &
                                             .OR. CCLOUD == 'LIMA' ) THEN
    IF (SIZE(XINPRC) /= 0 ) &         
      ZWORK21(:,:) = ZWORK21(:,:) + XINPRC(:,:) *3.6E6        
  END IF
  IF (CDCONV /= 'NONE') THEN
    ZWORK21(:,:) = ZWORK21(:,:) + XPRCONV(:,:)*3.6E6  
  END IF
  IF (LUSERR .OR. CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C3R5' .OR. &
                  CCLOUD == 'LIMA' .OR. CDCONV /= 'NONE') THEN
    TZFIELD%CMNHNAME   = 'INTOPR'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'INTOPR'
    TZFIELD%CUNITS     = 'mm hour-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_INstantaneous TOtal Precipitation Rate'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  ELSE
    PRINT * ,'YOU WANT TO COMPUTE THE RAIN RATE'
    PRINT * ,'BUT NO RAIN IS PRESENT IN THE MODEL' 
  END IF
!
  ! calculation of the mean instantaneous precipitations in the mesh-grid of a 
  ! large-scale model
  IF (LMEAN_PR .AND. LUSERR) THEN
    CALL COMPUTE_MEAN_PRECIP(ZWORK21,XMEAN_PR,ZWORK22,TZFIELD%NGRID)
!
    TZFIELD%CMNHNAME   = 'LS_INTOPR'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'LS_INTOPR'
    TZFIELD%CUNITS     = 'mm hour-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Large Scale INstantaneous TOtal Precipitation Rate'
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK22)
  END IF
!
END IF
!
!-------------------------------------------------------------------------------
!
!* CAPEMAX, CINMAX (corresponding to CAPEMAX), CAPE, CIN, DCAPE, VKE in J/kg
!
IF (NCAPE >=0 .AND. LUSERV) THEN
   ZWORK31(:,:,:) = XRT(:,:,:,1) * 1000.  ! vapour mixing ratio in g/kg
   ZWORK32(:,:,:)=0.0
   ZWORK33(:,:,:)=0.0
   ZWORK34(:,:,:)=0.0
   CALL CALCSOUND( XPABST(:,:,IKB:IKE)* 0.01 ,ZTEMP(:,:,IKB:IKE)- XTT, &
                   ZWORK31(:,:,IKB:IKE),                               &
                   ZWORK32(:,:,IKB:IKE),ZWORK33(:,:,IKB:IKE),          &
                   ZWORK34(:,:,IKB:IKE),ZWORK21,ZWORK22                )
  !
  TZFIELD%CMNHNAME   = 'CAPEMAX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'CAPEMAX'
  TZFIELD%CUNITS     = 'J kg-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_MAX of Convective Available Potential Energy'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  !
  TZFIELD%CMNHNAME   = 'CINMAX'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'CINMAX'
  TZFIELD%CUNITS     = 'J kg-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_MAX of Convective INhibition energy'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK22)
  !
  IF (NCAPE >=1) THEN
    TZFIELD%CMNHNAME   = 'CAPE3D'
    TZFIELD%CSTDNAME   = 'atmosphere_convective_available_potential_energy'
    TZFIELD%CLONGNAME  = 'CAPE3D'
    TZFIELD%CUNITS     = 'J kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Convective Available Potential Energy'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK32)
    !
    TZFIELD%CMNHNAME   = 'CIN3D'
    TZFIELD%CSTDNAME   = 'atmosphere_convective_inhibition'
    TZFIELD%CLONGNAME  = 'CIN3D'
    TZFIELD%CUNITS     = 'J kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Convective INhibition energy'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK33)
    !
    TZFIELD%CMNHNAME   = 'DCAPE3D'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'DCAPE3D'
    TZFIELD%CUNITS     = 'J kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = ''
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK34)
  END IF
  !
  IF (NCAPE >=2) THEN
    ZWORK31(:,:,1:IKU-1)= 0.5*(XWT(:,:,1:IKU-1)+XWT(:,:,2:IKU))
    ZWORK31(:,:,IKU)    = 0.
    ZWORK31=0.5*ZWORK31**2
    !
    TZFIELD%CMNHNAME   = 'VKE'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'VKE'
    TZFIELD%CUNITS     = 'J kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Vertical Kinetic Energy'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END IF
ENDIF
!
!-------------------------------------------------------------------------------
!
!* B-V frequency to assess thermal tropopause
!
IF (LBV_FR) THEN
  ZWORK32(:,:,:)=DZM(XTHT(:,:,:))/ MZM(XTHT(:,:,:))
  DO JK=1,IKU
   DO JJ=1,IJU
    DO JI=1,IIU
      IF(ZWORK32(JI,JJ,JK)<0.) THEN
        ZWORK31(JI,JJ,JK)= -1.*SQRT( ABS( XG*ZWORK32(JI,JJ,JK)/ XDZZ(JI,JJ,JK) ) )
      ELSE
        ZWORK31(JI,JJ,JK)= SQRT( ABS( XG*ZWORK32(JI,JJ,JK)/ XDZZ(JI,JJ,JK) ) )
      END IF
    ENDDO
   ENDDO
  ENDDO
  !
  TZFIELD%CMNHNAME   = 'BV'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'BV'
  TZFIELD%CUNITS     = 's-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Brunt-Vaissala frequency'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!  
  IF (NRR > 0) THEN
    ZWORK32(:,:,:)=DZM(ZTHETAE(:,:,:))/ MZM(ZTHETAE(:,:,:))
    DO JK=1,IKU
     DO JJ=1,IJU
      DO JI=1,IIU
        IF (ZWORK32(JI,JJ,JK)<0.) THEN
          ZWORK31(JI,JJ,JK)= -1.*SQRT( ABS( XG*ZWORK32(JI,JJ,JK)/ XDZZ(JI,JJ,JK) ) )
        ELSE
          ZWORK31(JI,JJ,JK)= SQRT( ABS( XG*ZWORK32(JI,JJ,JK)/ XDZZ(JI,JJ,JK) ) )
        END IF
      ENDDO
     ENDDO
    ENDDO
!
    TZFIELD%CMNHNAME   = 'BVE'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'BVE'
    TZFIELD%CUNITS     = 's-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Equivalent Brunt-Vaissala frequency'
    TZFIELD%NGRID      = 4
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END IF
END IF
!
IF(ALLOCATED(ZTHETAE)) DEALLOCATE(ZTHETAE)
IF(ALLOCATED(ZTHETAES)) DEALLOCATE(ZTHETAES)
!-------------------------------------------------------------------------------
!
!* GPS synthetic ZTD, ZHD, ZWD 
!
IF ( NGPS>=0 ) THEN
  !  surface temperature
  ZGAMREF=-6.5E-3
  ZWORK21(:,:) = ZTEMP(:,:,IKB) - ZGAMREF*((XZZ(:,:,IKB)+XZZ(:,:,IKB+1))/2.-XZS(:,:))
  !
  YFGRI=ADJUSTL(ADJUSTR(TPFILE%CNAME)//'GPS')
  CALL GPS_ZENITH (YFGRI,XRT(:,:,:,1),ZTEMP,XPABST,ZWORK21,ZWORK22,ZWORK23,ZWORK24)    
  !
  TZFIELD%CMNHNAME   = 'ZTD'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ZTD'
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Zenithal Total Delay'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK22)
  !
  IF (NGPS>=1) THEN
    TZFIELD%CMNHNAME   = 'ZHD'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'ZHD'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Zenithal Hydrostatic Delay'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK23)
    !
    TZFIELD%CMNHNAME   = 'ZWD'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'ZWD'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Zenithal Wet Delay'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK24)
    !
  END IF
  !
END IF
!
!-------------------------------------------------------------------------------
!
!* Radar reflectivities
!
IF(LRADAR .AND. LUSERR) THEN
! CASE  PREP_REAL_CASE after arome
  IF (CCLOUD=='NONE' .OR. CCLOUD=='KESS') THEN
    DEALLOCATE(XCIT)
    ALLOCATE(XCIT(IIU,IJU,IKU))    
    XCIT(:,:,:)=800.
    CALL INI_RADAR('PLAT')
  ELSE IF (CCLOUD=='LIMA') THEN
    DEALLOCATE(XCIT)
    ALLOCATE(XCIT(IIU,IJU,IKU))    
    XCIT(:,:,:)=XSVT(:,:,:,NSV_LIMA_NI)
    CALL INI_RADAR('PLAT')
  END IF
! 
  IF (NMOM_S.GE.2) ZW1(:,:,:)=XSVT(:,:,:,NSV_LIMA_NS)
  IF (NMOM_G.GE.2) ZW2(:,:,:)=XSVT(:,:,:,NSV_LIMA_NG)
  IF (NMOM_H.GE.2) ZW3(:,:,:)=XSVT(:,:,:,NSV_LIMA_NH)
  IF (NVERSION_RAD == 1) THEN 
! original version of radar diagnostics 
      WRITE(ILUOUT0,*) 'radar diagnostics from RADAR_RAIN_ICE routine'
    IF (CCLOUD=='LIMA') THEN
        CALL RADAR_RAIN_ICE (XRT, XCIT, XRHODREF, ZTEMP, ZWORK31, ZWORK32, &
                         ZWORK33, ZWORK34,XSVT(:,:,:,NSV_LIMA_NR),         &
                         ZW1(:,:,:),ZW2(:,:,:),ZW3(:,:,:) )
    ELSE          
      CALL RADAR_RAIN_ICE (XRT, XCIT, XRHODREF, ZTEMP, ZWORK31, ZWORK32, &
                                                         ZWORK33, ZWORK34 )
  ENDIF 
!
  TZFIELD%CMNHNAME   = 'RARE'
  TZFIELD%CSTDNAME   = 'equivalent_reflectivity_factor'
  TZFIELD%CLONGNAME  = 'RARE'
  TZFIELD%CUNITS     = 'dBZ'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_RAdar REflectivity'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!
  TZFIELD%CMNHNAME   = 'VDOP'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'VDOP'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_radar DOPpler fall speed'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK32)
!
  TZFIELD%CMNHNAME   = 'ZDR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ZDR'
  TZFIELD%CUNITS     = 'dBZ'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Differential polar Reflectivity'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK33)
!
  TZFIELD%CMNHNAME   = 'KDP'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'KDP'
  TZFIELD%CUNITS     = 'degree km-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Differential Phase Reflectivity'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK34)
!
   ELSE 
    !
    WRITE(ILUOUT0,*) 'radar diagnostics from RADAR_SIMULATOR routine'
    
    NBRAD=COUNT(XLAT_RAD(:) /= XUNDEF)
    NMAX=INT(NBSTEPMAX*XSTEP_RAD/XGRID)
    IF(NBSTEPMAX*XSTEP_RAD/XGRID/=NMAX .AND. (LCART_RAD)) THEN
      CALL PRINT_MSG(NVERB_FATAL,'GEN','WRITE_LFIFM1_FOR_DIAG', &
                    'NBSTEPMAX*XSTEP_RAD/XGRID is not an integer; please choose another combination')
    ENDIF
    DO JI=1,NBRAD
       NBELEV(JI)=COUNT(XELEV(JI,:) /= XUNDEF)
       WRITE(ILUOUT0,*) 'Number of ELEVATIONS : ', NBELEV(JI), 'FOR RADAR:', JI
    END DO
    IIELV=MAXVAL(NBELEV(1:NBRAD))
    WRITE(ILUOUT0,*) 'Maximum number of ELEVATIONS',IIELV
    WRITE(ILUOUT0,*) 'YOU HAVE ASKED FOR ', NBRAD, 'RADARS'
    !
    IF (LCART_RAD) NBAZIM=8*NMAX  ! number of azimuths 
    WRITE(ILUOUT0,*) ' Number of AZIMUTHS : ', NBAZIM
    IF (LCART_RAD) THEN
      ALLOCATE(ZWORK43(NBRAD,4*NMAX,2*NMAX))
    ELSE
      ALLOCATE(ZWORK43(1,NBAZIM,1))
    END IF
!! Some controls...
    IF(NBRAD/=COUNT(XLON_RAD(:) /= XUNDEF).OR.NBRAD/=COUNT(XALT_RAD(:) /= XUNDEF).OR. &
       NBRAD/=COUNT(XLAM_RAD(:) /= XUNDEF).OR.NBRAD/=COUNT(XDT_RAD(:) /= XUNDEF).OR. &
       NBRAD/=COUNT(CNAME_RAD(:) /= "UNDEF")) THEN
      CALL PRINT_MSG(NVERB_FATAL,'GEN','WRITE_LFIFM1_FOR_DIAG','inconsistency in DIAG1.nam')
    END IF
    IF(NCURV_INTERPOL==0.AND.(LREFR.OR.LDNDZ)) THEN
      LREFR=.FALSE.
      LDNDZ=.FALSE.
      WRITE(ILUOUT0,*) "Warning: cannot output refractivity nor its vertical gradient when NCURV_INTERPOL=0"
    END IF
    IF(MOD(NPTS_H,2)==0) THEN
      NPTS_H=NPTS_H+1
      WRITE(ILUOUT0,*) "Warning: NPTS_H has to be ODD. Setting it to ",NPTS_H
    END IF
    IF(MOD(NPTS_V,2)==0) THEN
      NPTS_V=NPTS_V+1
      WRITE(ILUOUT0,*) "Warning: NPTS_V has to be ODD. Setting it to ",NPTS_V
    END IF
    IF(LWBSCS.AND.LWREFL) THEN
      LWREFL=.FALSE.
      WRITE(ILUOUT0,*) "Warning: LWREFL cannot be set to .TRUE. if LWBSCS is also set to .TRUE.. Setting LWREFL to .FALSE.."
    END IF
    IF(CCLOUD=="LIMA" .AND. NDIFF/=7) THEN
      WRITE(YMSG,*) 'NDIFF=',NDIFF,' not available with CCLOUD=LIMA'
      CALL PRINT_MSG(NVERB_FATAL,'GEN','WRITE_LFIFM1_FOR_DIAG',YMSG)
    END IF
    INBOUT=28 !28: Temperature + RHR, RHS, RHG, ZDA, ZDS, ZDG, KDR, KDS, KDG      
    IF (CCLOUD=='LIMA') INBOUT=INBOUT+1 ! rain concentration CRT
    IF(LREFR) INBOUT=INBOUT+1 !+refractivity
    IF(LDNDZ) INBOUT=INBOUT+1 !+refractivity vertical gradient 
    IF(LATT)  INBOUT=INBOUT+12 !+AER-AEG AVR-AVG (vertical specific attenuation) and ATR-ATG 
    IF ( CCLOUD=='ICE4' ) THEN
      INBOUT=INBOUT+5 ! HAIL ZEH RHH ZDH KDH M_H 
      IF (LATT) THEN
        INBOUT=INBOUT+3  ! AEH AVH ATH
      ENDIF
    END IF
    WRITE(ILUOUT0,*) "Nombre de variables dans ZWORK42 en sortie de radar_simulator:",INBOUT

    IF (LCART_RAD) THEN
      ALLOCATE(ZWORK42(NBRAD,IIELV,2*NMAX,2*NMAX,INBOUT))
    ELSE
      ALLOCATE(ZWORK42(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,INBOUT))
      ALLOCATE(ZWORK42_BIS(NBRAD,IIELV,NBAZIM,NBSTEPMAX+1,INBOUT))
    END IF
    !
    IF (CCLOUD=='LIMA') THEN
      CALL RADAR_SIMULATOR(XUT,XVT,XWT,XRT,XSVT(:,:,:,NSV_LIMA_NI),XRHODREF,&
                      ZTEMP,XPABST,ZWORK42,ZWORK43,XSVT(:,:,:,NSV_LIMA_NR))
    ELSE ! ICE3
      CALL RADAR_SIMULATOR(XUT,XVT,XWT,XRT,XCIT,XRHODREF,ZTEMP,XPABSM,ZWORK42,ZWORK43)      
    ENDIF
    ALLOCATE(YRAD(INBOUT))
    YRAD(1:8)=(/"ZHH","ZDR","KDP","CSR","ZER","ZEI","ZES","ZEG"/)
    ICURR=9
    IF (CCLOUD=='ICE4') THEN
      YRAD(ICURR)="ZEH"
      ICURR=ICURR+1      
    END IF
    YRAD(ICURR)="VRU"
    ICURR=ICURR+1
    IF(LATT) THEN
      IF (CCLOUD=='ICE4') THEN
        YRAD(ICURR:ICURR+14)=(/"AER","AEI","AES","AEG","AEH","AVR","AVI","AVS","AVG","AVH","ATR","ATI","ATS","ATG","ATH"/)
        ICURR=ICURR+15
      ELSE
        YRAD(ICURR:ICURR+11)=(/"AER","AEI","AES","AEG","AVR","AVI","AVS","AVG","ATR","ATI","ATS","ATG"/)
        ICURR=ICURR+12
      END IF
    END IF
    YRAD(ICURR:ICURR+2)=(/"RHV","PDP","DHV"/)
    ICURR=ICURR+3
    YRAD(ICURR:ICURR+2)=(/"RHR","RHS","RHG"/)
    ICURR=ICURR+3
    IF (CCLOUD=='ICE4') THEN
      YRAD(ICURR)="RHH"
      ICURR=ICURR+1      
    END IF
    YRAD(ICURR:ICURR+2)=(/"ZDA","ZDS","ZDG"/)
    ICURR=ICURR+3
    IF (CCLOUD=='ICE4') THEN
      YRAD(ICURR)="ZDH"
      ICURR=ICURR+1      
    END IF
    YRAD(ICURR:ICURR+2)=(/"KDR","KDS","KDG"/)
    ICURR=ICURR+3
    IF (CCLOUD=='ICE4') THEN
      YRAD(ICURR)="KDH"
      ICURR=ICURR+1      
    END IF
    YRAD(ICURR:ICURR+4)=(/"HAS","M_R","M_I","M_S","M_G"/)
    ICURR=ICURR+5
    IF (CCLOUD=='ICE4') THEN
      YRAD(ICURR)="M_H"
      ICURR=ICURR+1      
    END IF
    YRAD(ICURR:ICURR+1)=(/"CIT","TEM"/)
    ICURR=ICURR+2
    IF (CCLOUD=='LIMA') THEN
      YRAD(ICURR)="CRT"
      ICURR=ICURR+1
    ENDIF
    IF(LREFR) THEN
      YRAD(ICURR)="RFR"
      ICURR=ICURR+1
    END IF
    IF(LDNDZ) THEN
      YRAD(ICURR)="DNZ"
      ICURR=ICURR+1
    END IF
    IF (LCART_RAD) THEN
      DO JI=1,NBRAD
        IEL=NBELEV(JI)
        ! writing latlon in internal files 
        ALLOCATE(CLATLON(2*NMAX))
        CLATLON=""
        DO JV=2*NMAX,1,-1
          DO JH=1,2*NMAX
            WRITE(CBUFFER,'(2(f8.3,1X))') ZWORK43(JI,2*JH-1,JV),ZWORK43(JI,2*JH,JV)
                  CLATLON(JV)=TRIM(CLATLON(JV)) // " " // TRIM(CBUFFER)
          END DO
          CLATLON(JV)=TRIM(ADJUSTL(CLATLON(JV)))
        END DO            
        DO JEL=1,IEL
          WRITE(YELEV,'(I2.2,A1,I1.1)') FLOOR(XELEV(JI,JEL)),'.',&
                 INT(ANINT(10.*XELEV(JI,JEL))-10*INT(XELEV(JI,JEL)))
          WRITE(YGRID_SIZE,'(I3.3)') 2*NMAX
          DO JJ=1,SIZE(ZWORK42(:,:,:,:,:),5)
            YRS=YRAD(JJ)//CNAME_RAD(JI)(1:3)//YELEV//YGRID_SIZE//TRIM(TPFILE%CNAME)
            CALL IO_File_add2list(TZRSFILE,YRS,'TXT','WRITE',KRECL=8192)
            CALL IO_File_open(TZRSFILE,HSTATUS='NEW')
            ILURS = TZRSFILE%NLU
            WRITE(ILURS,'(A,4F12.6,2I5)') '**domaine LATLON ',ZWORK43(JI,1,1),ZWORK43(JI,4*NMAX-1,2*NMAX), &
                  ZWORK43(JI,2,1),ZWORK43(JI,4*NMAX,2*NMAX),2*NMAX,2*NMAX !! HEADER
            DO JV=2*NMAX,1,-1
              DO JH=1,2*NMAX
                WRITE(ILURS,'(E11.5,1X)',ADVANCE='NO') ZWORK42(JI,JEL,JH,JV,JJ)
              END DO
              WRITE(ILURS,*) ''
            END DO
                  
            DO JV=2*NMAX,1,-1
              WRITE(ILURS,*) CLATLON(JV)
            END DO                  
            CALL IO_File_close(TZRSFILE)
            TZRSFILE => NULL()
          END DO               
        END DO
        DEALLOCATE(CLATLON)
      END DO
    ELSE ! polar output
       CALL MPI_ALLREDUCE(ZWORK42, ZWORK42_BIS, SIZE(ZWORK42), MNHREAL_MPI, MPI_MAX, NMNH_COMM_WORLD, IERR)
      DO JI=1,NBRAD
        IEL=NBELEV(JI)
        DO JEL=1,IEL
          WRITE(YELEV,'(I2.2,A1,I1.1)') FLOOR(XELEV(JI,JEL)),'.',&
                INT(ANINT(10.*XELEV(JI,JEL))-10*INT(XELEV(JI,JEL)))
          DO JJ=1,SIZE(ZWORK42(:,:,:,:,:),5)
            YRS="P"//YRAD(JJ)//CNAME_RAD(JI)(1:3)//YELEV//TRIM(TPFILE%CNAME)
            CALL IO_File_add2list(TZRSFILE,YRS,'TXT','WRITE')
            CALL IO_File_open(TZRSFILE)
            ILURS = TZRSFILE%NLU
            DO JH=1,NBAZIM
              DO JV=1,NBSTEPMAX+1
                WRITE(ILURS,"(F15.7)") ZWORK42_BIS(JI,JEL,JH,JV,JJ)
              END DO
            END DO                                    
            CALL IO_File_close(TZRSFILE)
            TZRSFILE => NULL()
          END DO
        END DO
      END DO
    END IF !polar output
    DEALLOCATE(ZWORK42,ZWORK43)
   END IF
END IF
!
IF (LLIDAR) THEN
  PRINT *,'CALL LIDAR/RADAR with TPFILE%CNAME =',TPFILE%CNAME
  YVIEW='     '
  YVIEW=TRIM(CVIEW_LIDAR)
  PRINT *,'CVIEW_LIDAR REQUESTED ',YVIEW
  IF (YVIEW/='NADIR'.AND.YVIEW/='ZENIT') YVIEW='NADIR'
  PRINT *,'CVIEW_LIDAR USED ',YVIEW
  PRINT *,'XALT_LIDAR REQUESTED (m) ',XALT_LIDAR
  PRINT *,'XWVL_LIDAR REQUESTED (m) ',XWVL_LIDAR
  IF (XWVL_LIDAR==XUNDEF) XWVL_LIDAR=0.532E-6
  IF (XWVL_LIDAR<1.E-7.OR.XWVL_LIDAR>2.E-6) THEN
    PRINT *,'CAUTION: THE XWVL_LIDAR REQUESTED IS OUTSIDE THE USUAL RANGE'
    XWVL_LIDAR=0.532E-6
  ENDIF
  PRINT *,'XWVL_LIDAR USED (m) ',XWVL_LIDAR
!
  IF (LDUST) THEN
    IACCMODE=MIN(2,NMODE_DST)
    ALLOCATE(ZTMP1(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3), 1))
    ALLOCATE(ZTMP2(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3), 1))
    ALLOCATE(ZTMP3(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3), 1))
    ZTMP1(:,:,:,1)=ZN0_DST(:,:,:,IACCMODE)
    ZTMP2(:,:,:,1)=ZRG_DST(:,:,:,IACCMODE)
    ZTMP3(:,:,:,1)=ZSIG_DST(:,:,:,IACCMODE)
    SELECT CASE ( CCLOUD )
    CASE('KESS''ICE3','ICE4')
      CALL LIDAR(CCLOUD, YVIEW, XALT_LIDAR, XWVL_LIDAR, XZZ, XRHODREF, ZTEMP, XCLDFR, &
                 XRT, ZWORK31, ZWORK32,                                        &
                 PDSTC=ZTMP1,                                                  &
                 PDSTD=ZTMP2,                                                  &
                 PDSTS=ZTMP3)
    CASE('C2R2')
      CALL LIDAR(CCLOUD, YVIEW, XALT_LIDAR, XWVL_LIDAR, XZZ, XRHODREF, ZTEMP, XCLDFR, &
                 XRT, ZWORK31, ZWORK32,                                        &
                 PCT=XSVT(:,:,:,NSV_C2R2BEG+1:NSV_C2R2END),                    &
                 PDSTC=ZTMP1,                                                  &
                 PDSTD=ZTMP2,                                                  &
                 PDSTS=ZTMP3)
    CASE('C3R5')
      CALL LIDAR(CCLOUD, YVIEW, XALT_LIDAR, XWVL_LIDAR, XZZ, XRHODREF, ZTEMP, XCLDFR, &
                 XRT, ZWORK31, ZWORK32,                                        &
                 PCT=XSVT(:,:,:,NSV_C2R2BEG+1:NSV_C1R3END-1),                  &
                 PDSTC=ZTMP1,                                                  &
                 PDSTD=ZTMP2,                                                  &
                 PDSTS=ZTMP3)
    CASE('LIMA')
! PCT(2) = droplets (3)=drops (4)=ice crystals
       ALLOCATE(ZTMP4(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3), 4))
       ZTMP4(:,:,:,1)=0.
       ZTMP4(:,:,:,2)=XSVT(:,:,:,NSV_LIMA_NC)
       ZTMP4(:,:,:,3)=XSVT(:,:,:,NSV_LIMA_NR)
       ZTMP4(:,:,:,4)=XSVT(:,:,:,NSV_LIMA_NI)
!
       CALL LIDAR(CCLOUD, YVIEW, XALT_LIDAR, XWVL_LIDAR, XZZ, XRHODREF, ZTEMP, MAX(XCLDFR,XICEFR),&
            XRT, ZWORK31, ZWORK32,                                  &
            PCT=ZTMP4,                            &
            PDSTC=ZTMP1,                          &
            PDSTD=ZTMP2,                          &
            PDSTS=ZTMP3)
!
    END SELECT
  ELSE
    SELECT CASE ( CCLOUD )
    CASE('KESS','ICE3','ICE4')
      CALL LIDAR(CCLOUD, YVIEW, XALT_LIDAR, XWVL_LIDAR, XZZ, XRHODREF, ZTEMP, XCLDFR, &
           XRT, ZWORK31, ZWORK32)
    CASE('C2R2')
      CALL LIDAR(CCLOUD, YVIEW, XALT_LIDAR, XWVL_LIDAR, XZZ, XRHODREF, ZTEMP, XCLDFR, &
           XRT, ZWORK31, ZWORK32,                                  &
           PCT=XSVT(:,:,:,NSV_C2R2BEG+1:NSV_C2R2END))
    CASE('C3R5')
      CALL LIDAR(CCLOUD, YVIEW, XALT_LIDAR, XWVL_LIDAR, XZZ, XRHODREF, ZTEMP, XCLDFR, &
           XRT, ZWORK31, ZWORK32,                                  &
           PCT=XSVT(:,:,:,NSV_C2R2BEG+1:NSV_C1R3END-1))
    CASE('LIMA')
! PCT(2) = droplets (3)=drops (4)=ice crystals
       ALLOCATE(ZTMP4(SIZE(XSVT,1), SIZE(XSVT,2), SIZE(XSVT,3), 4))
       ZTMP4(:,:,:,1)=0.
       ZTMP4(:,:,:,2)=XSVT(:,:,:,NSV_LIMA_NC)
       ZTMP4(:,:,:,3)=XSVT(:,:,:,NSV_LIMA_NR)
       ZTMP4(:,:,:,4)=XSVT(:,:,:,NSV_LIMA_NI)
!
       CALL LIDAR(CCLOUD, YVIEW, XALT_LIDAR, XWVL_LIDAR, XZZ, XRHODREF, ZTEMP, MAX(XCLDFR,XICEFR),&
            XRT, ZWORK31, ZWORK32,                                  &
            PCT=ZTMP4)
    END SELECT
  ENDIF
!
  IF( ALLOCATED(ZTMP1) ) DEALLOCATE(ZTMP1)
  IF( ALLOCATED(ZTMP2) ) DEALLOCATE(ZTMP2)
  IF( ALLOCATED(ZTMP3) ) DEALLOCATE(ZTMP3)
  IF( ALLOCATED(ZTMP4) ) DEALLOCATE(ZTMP4)
!
  TZFIELD%CMNHNAME   = 'LIDAR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'LIDAR'
  TZFIELD%CUNITS     = 'm-1 sr-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Normalized_Lidar_Profile'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!
  TZFIELD%CMNHNAME   = 'LIPAR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'LIPAR'
  TZFIELD%CUNITS     = 'm-1 sr-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_Particle_Lidar_Profile'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK32)
!
END IF
!
!-------------------------------------------------------------------------------
!
!* Height of boundary layer
!
IF (CBLTOP == 'THETA') THEN
  !
  ! methode de la parcelle
  !
  ALLOCATE(ZSHMIX(IIU,IJU))

  ZWORK31(:,:,1:IKU-1)=0.5*(XZZ(:,:,1:IKU-1)+XZZ(:,:,2:IKU))
  ZWORK31(:,:,IKU)=2.*ZWORK31(:,:,IKU-1)-ZWORK31(:,:,IKU-2)
  ZWORK21(:,:) = ZTHETAV(:,:,IKB)+0.5
  ZSHMIX(:,:)  = 0.0
  DO JJ=1,IJU
    DO JI=1,IIU 
      DO JK=IKB,IKE
        IF ( ZTHETAV(JI,JJ,JK).GT.ZWORK21(JI,JJ) ) THEN
          ZSHMIX(JI,JJ) =  ZWORK31(JI,JJ,JK-1)  &
                        +( ZWORK31(JI,JJ,JK) - ZWORK31 (JI,JJ,JK-1) )  &
                        /( ZTHETAV(JI,JJ,JK) - ZTHETAV(JI,JJ,JK-1) )  &
                        *( ZWORK21(JI,JJ)    - ZTHETAV(JI,JJ,JK-1) )
          EXIT
        END IF
      END DO
    END DO
  END DO
  ZSHMIX(:,:)=ZSHMIX(:,:)-XZS(:,:)
  ZSHMIX(:,:)=MAX(ZSHMIX(:,:),50.0)
  !
  TZFIELD%CMNHNAME   = 'HBLTOP'
  TZFIELD%CSTDNAME   = 'atmosphere_boundary_layer_thickness'
  TZFIELD%CLONGNAME  = 'HBLTOP'
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'Height of Boundary Layer TOP'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZSHMIX)
  !
  DEALLOCATE(ZSHMIX)
ELSEIF (CBLTOP == 'RICHA') THEN
  !
  ! methode du "bulk Richardson number"
  !
  ALLOCATE(ZRIB(IIU,IJU,IKU))
  ALLOCATE(ZSHMIX(IIU,IJU))

  ZWORK31(:,:,1:IKU-1)=0.5*(XZZ(:,:,1:IKU-1)+XZZ(:,:,2:IKU))
  ZWORK31(:,:,IKU)=2.*ZWORK31(:,:,IKU-1)-ZWORK31(:,:,IKU-2)
  ZWORK32=MXF(XUT)
  ZWORK33=MYF(XVT)
  ZWORK34=ZWORK32**2+ZWORK33**2
  DO JK=IKB,IKE  
    ZRIB(:,:,JK)=XG*ZWORK31(:,:,JK)*(ZTHETAV(:,:,JK)-ZTHETAV(:,:,IKB))/(ZTHETAV(:,:,IKB)*ZWORK34(:,:,JK))
  ENDDO
  ZSHMIX=0.0
  DO JJ=1,IJU
    DO JI=1,IIU 
      DO JK=IKB,IKE
        IF ( ZRIB(JI,JJ,JK).GT.0.25 ) THEN
          ZSHMIX(JI,JJ) = ZWORK31(JI,JJ,JK-1)  &
                        +( ZWORK31(JI,JJ,JK) - ZWORK31(JI,JJ,JK-1) )  &
                        *( 0.25 - ZRIB(JI,JJ,JK-1) )  &
                        /( ZRIB(JI,JJ,JK)    - ZRIB(JI,JJ,JK-1) )
          EXIT
        END IF
      END DO
    END DO
  END DO
  ZSHMIX(:,:)=ZSHMIX(:,:)-XZS(:,:)
  !
  TZFIELD%CMNHNAME   = 'HBLTOP'
  TZFIELD%CSTDNAME   = 'atmosphere_boundary_layer_thickness'
  TZFIELD%CLONGNAME  = 'HBLTOP'
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'Height of Boundary Layer TOP'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZSHMIX)
  !
  DEALLOCATE(ZRIB,ZSHMIX)
ENDIF
  ! used before 5-3-1 version
  !
  !ZGAMREF=3.5E-3 ! K/m
  !ZWORK31(:,:,1:IKU-1)=0.5*(XZZ(:,:,1:IKU-1)+XZZ(:,:,2:IKU))
  !ZWORK31(:,:,IKU)=2.*ZWORK31(:,:,IKU-1)-ZWORK31(:,:,IKU-2)
  !CALL FREE_ATM_PROFILE(ZTHETAV,ZWORK31,XZS,XZSMT,ZGAMREF,ZWORK32,ZWORK33)
!
IF (ALLOCATED(ZTHETAV)) DEALLOCATE(ZTHETAV)
!
!
!* Ligthning
!
IF ( LCH_CONV_LINOX ) THEN 
  CALL IO_Field_write(TPFILE,'IC_RATE',    XIC_RATE)
  CALL IO_Field_write(TPFILE,'CG_RATE',    XCG_RATE)
  CALL IO_Field_write(TPFILE,'IC_TOTAL_NB',XIC_TOTAL_NUMBER)
  CALL IO_Field_write(TPFILE,'CG_TOTAL_NB',XCG_TOTAL_NUMBER)
END IF
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       1.8    My own  variables :
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_LFIFM1_FOR_DIAG  
