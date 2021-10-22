!     ######spl
      MODULE MODD_CH_MNHC_n
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #####################
!!
!!*** *MODD_CH_MNHC$n*
!!
!!    PURPOSE
!!    -------
!       This module contains all parameters that control the
!     chemical part of MesoNH.
!!
!!**  AUTHOR
!!    ------
!!    K. Suhre      *Laboratoire d'Aerollogie*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/05/95
!!    27/07/96 (K. Suhre) restructured
!!    30/11/99 (K. Suhre) add new parameters
!!    16/11/00 (C. Mari)  add new parameters
!!    28/05/02 (C. Mari)  move default values to default_desfmn
!!    13/07/03 (J.-P. Pinty) add flag for lightning production of NOx
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CH_MNHC_t
!
! switch which indicates whether chemistry is used or not
!
  LOGICAL :: LUSECHEM
!
!* Initialization
!
  LOGICAL :: LCH_INIT_FIELD ! flag indicating whether initialization
                            ! of chemical fields shall be done during MesoNH run using
                            ! CH_INIT_FIELD (overwrites initial values from FM-files)
                            ! or not
!
!* Surface options
!
  LOGICAL :: LCH_SURFACE_FLUX  ! flag indicating whether surface flux
                               ! for chemical species shall be calculated using
                               ! CH_SURFACE_FLUX or not (dry deposition and emission)
!
!* Scavenging
!
  LOGICAL :: LCH_CONV_SCAV
                 ! flag for calculation of scavenging
                 ! by convective precipitations (active only if LCHTRANS=.TRUE.)
!
  CHARACTER(LEN=10) :: CCH_EXPLICIT_SCAV
               ! wet deposition method (not yet implemented)
!
!*
  CHARACTER(LEN=10) :: CCH_SCHEME
                 ! name of chemical scheme
!
  CHARACTER(LEN=80) :: CCHEM_INPUT_FILE
                 ! name of general
                 ! purpose ASCII input file (handeled by CH_OPEN_INPUT)
!
  CHARACTER(LEN=10) :: CCH_TDISCRETIZATION
                 ! temporal discretization:
                 ! "SPLIT"  : use time-splitting, input fields for solver are
                 !            scalar variables at t+dt (derived from XRSVS)
                 ! "CENTER" : input fields for solver are
                 !            scalar variables at t (XSVT)
                 ! "LAGGED" : input fields for solver are
                 !            scalar variables at t-dt (XSVM)
!
  INTEGER           :: NCH_SUBSTEPS
                 ! number of chemical timesteps to be taken during one
                 ! double timestep of MesoNH (MesoNH integrates with timesteps
                 ! of lenght 2*XTSTEP using leapfrog), the timestep of the
                 ! solver will be calculated as
                 ! ZDTSOLVER = 2*XTSTEP / NCH_SUBSTEPS
!* LiNOx
!
LOGICAL :: LCH_CONV_LINOX
                 ! flag for calculation of NOx production by lightnings
!
!* photolysis rates (TUV)
!
  LOGICAL      :: LCH_TUV_ONLINE  ! switch online/lookup table
  CHARACTER*80 :: CCH_TUV_LOOKUP  ! name of lookup table file
  CHARACTER*4  :: CCH_TUV_CLOUDS  ! method for calculating the
                                ! impact of clouds on radiation
                                ! "FOUQ" (model clouds, only 1-D)
  REAL :: XCH_TUV_ALBNEW  ! surface albedo (if negative the albedo
                        ! will be read from DATAX/albedo.dat)
  REAL :: XCH_TUV_DOBNEW  ! scaling factor for ozone column dobson
                        ! (if negative, no scaling will be performed,
                        ! note: the O3 profile will be read from
                        ! DATAX/O3.profile, if this file is empty, the
                        ! US standard O3 profile will be used)
  REAL :: XCH_TUV_TUPDATE ! update frequency for TUV (in seconds)
!
!* vectorization
!
  CHARACTER*3 :: CCH_VEC_METHOD          ! type of vectorization mask
                                       ! 'MAX' take NCH_VEC_LENGTH points
                                       ! 'TOT' take all grid points
                                       ! 'HOR' take horizontal layers
                                       ! 'VER' take vertical columns
  INTEGER     :: NCH_VEC_LENGTH          ! number of points for 'MAX' option
!
!* 1-D time series
!
  REAL         :: XCH_TS1D_TSTEP         ! time between two call to write_ts1d
  CHARACTER*80 :: CCH_TS1D_COMMENT       ! comment for write_ts1d
  CHARACTER*80 :: CCH_TS1D_FILENAME      ! filename for write_ts1d files
!
END TYPE CH_MNHC_t

TYPE(CH_MNHC_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_MNHC_MODEL

LOGICAL, POINTER :: LUSECHEM=>NULL()
LOGICAL, POINTER :: LCH_INIT_FIELD=>NULL()
LOGICAL, POINTER :: LCH_SURFACE_FLUX=>NULL()
LOGICAL, POINTER :: LCH_CONV_SCAV=>NULL()
CHARACTER(LEN=10), POINTER :: CCH_EXPLICIT_SCAV=>NULL()
CHARACTER(LEN=10), POINTER :: CCH_SCHEME=>NULL()
CHARACTER(LEN=80), POINTER :: CCHEM_INPUT_FILE=>NULL()
CHARACTER(LEN=10), POINTER :: CCH_TDISCRETIZATION=>NULL()
INTEGER, POINTER :: NCH_SUBSTEPS=>NULL()
LOGICAL, POINTER :: LCH_CONV_LINOX=>NULL()
LOGICAL, POINTER :: LCH_TUV_ONLINE=>NULL()
CHARACTER*80, POINTER :: CCH_TUV_LOOKUP=>NULL()
CHARACTER*4, POINTER :: CCH_TUV_CLOUDS=>NULL()
REAL, POINTER :: XCH_TUV_ALBNEW=>NULL()
REAL, POINTER :: XCH_TUV_DOBNEW=>NULL()
REAL, POINTER :: XCH_TUV_TUPDATE=>NULL()
CHARACTER*3, POINTER :: CCH_VEC_METHOD=>NULL()
INTEGER, POINTER :: NCH_VEC_LENGTH=>NULL()
REAL, POINTER :: XCH_TS1D_TSTEP=>NULL()
CHARACTER*80, POINTER :: CCH_TS1D_COMMENT=>NULL()
CHARACTER*80, POINTER :: CCH_TS1D_FILENAME=>NULL()

CONTAINS

SUBROUTINE CH_MNHC_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODD_CH_MNHC_N:CH_MNHC_GOTO_MODEL',0,ZHOOK_HANDLE)
LUSECHEM=>CH_MNHC_MODEL(KTO)%LUSECHEM
LCH_INIT_FIELD=>CH_MNHC_MODEL(KTO)%LCH_INIT_FIELD
LCH_SURFACE_FLUX=>CH_MNHC_MODEL(KTO)%LCH_SURFACE_FLUX
LCH_CONV_SCAV=>CH_MNHC_MODEL(KTO)%LCH_CONV_SCAV
CCH_EXPLICIT_SCAV=>CH_MNHC_MODEL(KTO)%CCH_EXPLICIT_SCAV
CCH_SCHEME=>CH_MNHC_MODEL(KTO)%CCH_SCHEME
CCHEM_INPUT_FILE=>CH_MNHC_MODEL(KTO)%CCHEM_INPUT_FILE
CCH_TDISCRETIZATION=>CH_MNHC_MODEL(KTO)%CCH_TDISCRETIZATION
NCH_SUBSTEPS=>CH_MNHC_MODEL(KTO)%NCH_SUBSTEPS
LCH_CONV_LINOX=>CH_MNHC_MODEL(KTO)%LCH_CONV_LINOX
LCH_TUV_ONLINE=>CH_MNHC_MODEL(KTO)%LCH_TUV_ONLINE
CCH_TUV_LOOKUP=>CH_MNHC_MODEL(KTO)%CCH_TUV_LOOKUP
CCH_TUV_CLOUDS=>CH_MNHC_MODEL(KTO)%CCH_TUV_CLOUDS
XCH_TUV_ALBNEW=>CH_MNHC_MODEL(KTO)%XCH_TUV_ALBNEW
XCH_TUV_DOBNEW=>CH_MNHC_MODEL(KTO)%XCH_TUV_DOBNEW
XCH_TUV_TUPDATE=>CH_MNHC_MODEL(KTO)%XCH_TUV_TUPDATE
CCH_VEC_METHOD=>CH_MNHC_MODEL(KTO)%CCH_VEC_METHOD
NCH_VEC_LENGTH=>CH_MNHC_MODEL(KTO)%NCH_VEC_LENGTH
XCH_TS1D_TSTEP=>CH_MNHC_MODEL(KTO)%XCH_TS1D_TSTEP
CCH_TS1D_COMMENT=>CH_MNHC_MODEL(KTO)%CCH_TS1D_COMMENT
CCH_TS1D_FILENAME=>CH_MNHC_MODEL(KTO)%CCH_TS1D_FILENAME

IF (LHOOK) CALL DR_HOOK('MODD_CH_MNHC_N:CH_MNHC_GOTO_MODEL',1,ZHOOK_HANDLE)
END SUBROUTINE CH_MNHC_GOTO_MODEL

END MODULE MODD_CH_MNHC_n
