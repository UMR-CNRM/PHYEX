!     ######spl
      MODULE MODD_CH_MODEL0D
!!    ######################
!!
!!*** *MODD_CH_MODEL0D*
!!
!!    PURPOSE
!!    -------
!       contains all control variables of the box model,
!     all timing variables, file names, formats and I/O channels
!     used by the box-model are declared here
!!
!!**  AUTHOR
!!    ------
!!    K. Suhre     *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 02/03/95
!!    27/07/96 (K. Suhre) restructured
!!    16/02/99 (K. Suhre) add new parameters for TUV
!!              XCH_TUV_DOBNEW  : scaling factor fir column dobson
!!              XCH_TUV_ALBNEW  : surface albedo
!!              XCH_TUV_TUPDATE : update frequency for TUV
!!              LCH_TUV_ONLINE  : switch between online/lookup table
!!              CCH_TUV_LOOKUP  : name of lookup table file
!!              CCH_TUV_CLOUDS  : method for the impact of clouds on radiation
!!    18/02/99  LCH_SURFACE0D   : apply surface emission / deposition fluxes
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
!
IMPLICIT NONE
SAVE
!
!*       0.1  integration and I/O timing control
!
REAL :: XTBEGIN = 0.0   ! begin the simulation in seconds
REAL :: XTEND   = 0.0   ! end the simulation in seconds
REAL :: XDTOUT  = 900.  ! timestep for printout of intermediate results
REAL :: XDTDIAG = 3600. ! timestep for calculation of diagnostics
!
REAL :: XTSIMUL       ! time of simulation
REAL :: XDTACT = 40.  ! timestep of the simulation
REAL :: XTNEXTOUT     ! time of next result output
REAL :: XTNEXTDIAG    ! time of next diagnosis
REAL :: XTNEXTMETEO   ! time of next meteo update
!
!*       0.2  file names, formats and I/O channels
!
CHARACTER*128 :: CINITFILE  = "CHCONTROL1.nam" ! name of initial value file
CHARACTER*128 :: CMETEOFILE = "CHCONTROL1.nam" ! meteo update file
!
CHARACTER*128 :: COUTFILE    = "BOX.OUT"      ! name of final output file
CHARACTER*128 :: CRESULTFILE = "BOX.RESULT"   ! regular output file
CHARACTER*128 :: CDIAGFILE   = "BOX.DIAG"     ! diagnostics output file
!
CHARACTER*80  :: CRUNID        = "no runid specified" ! runid for output file
CHARACTER*40  :: CRESULTFORMAT = "(5E16.8)" ! Format for results
CHARACTER*40  :: CDIAGFORMAT   = "(5E16.8)" ! Format for diagnostics
!
INTEGER :: NFILEIO   = 42 ! channel to be used for all intermediate file I/O
INTEGER :: NMETEOIO  = 43 ! channel to be used for all meteo I/O
INTEGER :: NRESULTIO = 44 ! channel to be used for all regular result file I/O
INTEGER :: NDIAGIO   = 45 ! channel to be used for all diagnostics file I/O
!
!*       0.3  verbosity level
!
INTEGER :: NVERB     = 5  ! verbosity level: 0 (lowest) <= NVERB <= 10 (highest)
!
!*       0.4  parameters for TUV
!
LOGICAL      :: LCH_TUV_ONLINE = .TRUE.        ! switch online/lookup table
CHARACTER*80 :: CCH_TUV_LOOKUP = "PHOTO.TUV39" ! name of lookup table file
CHARACTER*4  :: CCH_TUV_CLOUDS = "NONE"        ! method for calculating the
                                               ! impact of clouds on radiation
                                               ! "FOUQ" (model clouds, only 1-D)
                                               ! "RADM" (parameterized, for 3-D)
REAL :: XCH_TUV_ALBNEW  = -1.  ! surface albedo (if negative the albedo
                                               ! will be read from DATAX/albedo.dat)
REAL :: XCH_TUV_DOBNEW  = -1.  ! scaling factor for ozone column dobson
                               ! (if negative, no scaling will be performed,
                               ! note: the O3 profile will be read from
                               ! DATAX/O3.profile, if this file is empty, the
                               ! US standard O3 profile will be used)
REAL :: XCH_TUV_TUPDATE = 600. ! update frequency for TUV (in seconds)
!
LOGICAL :: LCH_SURFACE0D = .FALSE. ! switch to activate surface fluxes
!
CHARACTER(LEN=80) :: CCHEM_INPUT_FILE
END MODULE MODD_CH_MODEL0D
