!     ######spl
      MODULE MODD_CH_M9_SCHEME
!!    ########################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
!!
!!*** *MODD_CH_M9_SCHEME*
!!
!!    PURPOSE
!!    -------
!     definition of variables and types for the chemical core system
!!
!!**  METHOD
!!    ------
!!      All constants and auxiliary variables are stored in one common
!!    data type (CCSTYPE). This allows to pass them all as one single
!!    variable in the argument lists of the CCS.
!!    The constants NEQ and NREAC are duplicated here in order to avoid
!!    decouple the CCS from the other modules of MNHC.
!!    Variables to be transfered from the meteorological part are stored
!!    in the data type METEOTRANSTYPE (number, value and name).
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Karsten Suhre (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 26/07/96
!!    Modified 05/05/98: Vectorization (Vincent Crassier & KS)
!!    Modified 31/10/03: New interface for better MesoNH compilation (D. Gazen)
!!
!!----------------------------------------------------------------------
!!    DECLARATIONS
!!    ------------
IMPLICIT NONE
INTEGER, PARAMETER :: NEQ           = 40 ! number of prognostic chemical species
INTEGER, PARAMETER :: NREAC         = 128 ! number of chemical reactions
INTEGER, PARAMETER :: NMETEOVARS    =  10 ! number of meteorological variables
INTEGER, PARAMETER :: NNONZEROTERMS = 593 ! number of non-zero terms returned by CH_TERMS
CHARACTER(LEN=32),  DIMENSION(NEQ),   TARGET :: CNAMES ! names of the species
CHARACTER(LEN=32),  DIMENSION(NREAC), TARGET :: CREACS ! the reaction rate names
CHARACTER(LEN=256), DIMENSION(NREAC), TARGET :: CFULLREACS ! the full reactions
!
TYPE CCSTYPE ! reaction rates and auxiliary variables
 REAL,DIMENSION(:),POINTER :: K001,K002,K003,K004,K005,K006,K007,K008,K009,K010&
&,K011,K012,K013,K014,K015,K016,K017,K018,K019,K020,K021,K022,K023,K024,K025,K0&
&26,K027,K028,K029,K030,K031,K032,K033,K034,K035,K036,K037,K038,K039,K040,K041,&
&K042,K043,K044,K045,K046,K047,K048,K049,K050,K051,K052,K053,K054,K055,K056,K05&
&7,K058,K059,K060,K061,K062,K063,K064,K065,K066,K067,K068,K069,K070,K071,K072,K&
&073,K074,K075,K076,K077,K078,K079,K080,K081,K082,K083,K084,K085,K086,K087,K088&
&,K089,K090,K091,K092,K093,K094,K095,K096,K097,K098,K099,K100,K101,K102,K103,K1&
&04,K105,K106,K107,K108,K109,K110,K111,K112,K113,K114,K115,K116,K117,K118,K119,&
&K120,K121,K122,K123,K124,K125,K126,K127,K128
! output channel (NOUT) and verbosity level (NVERB)
INTEGER,DIMENSION(:),POINTER :: NOUT=>NULL(),NVERB=>NULL()
! auxiliary variables defined by the user, if any (e.g. O2, N2, H2O)
! /BEGIN_MODULE/
 !
 ! supplementary variables of the CCS that are to be placed into
 ! the TYPE definition of TPK (to be addressed e.g. as TPK%O2):
 !
 INTEGER,DIMENSION(:),POINTER :: MODELLEVEL       ! index of the model level (1 for box model)
 REAL,DIMENSION(:),POINTER    :: T,              &! temperature (K)
            PRESSURE,       &! pressure (atm)
            M,              &! air density (molec/cm3)
            H2O,            &! conc. of water molecules (molec/cm3)
            CLOUDWATER,     &! cloud water (kg/kg)
            O2, N2, H2,     &! conc. of oxigen nitrogen, hydrogen (molec/cm3)
            OH, O1D, O3P,   &! (molec/cm3) at equilibrium (fast species)
            LON,            &! longitude of curtrent grid point (degree)
            LAT              ! latitude of curtrent grid point (degree)

 INTEGER,DIMENSION(:),POINTER :: YEAR, MONTH, DAY ! starting date of experiment (~DTEXP)
 !
! /END_MODULE/
END TYPE CCSTYPE
!
! Use array of CCSTYPE to handle the 8 possible models :
! TACCS(i) refers to the CCSTYPE variable of the ith model
! You should declare a TYPE(CCSTYPE) pointer variable TZK to point to
! TACCS(i) in each subroutine that deals with CCSTYPE variables :
!
! TYPE(CCSTYPE),POINTER :: TZK
!
! TZK=>TACCS(KMI)
!
TYPE(CCSTYPE), DIMENSION(8), TARGET, SAVE :: TACCS ! 8 models
!
! list of chemical species indices
INTEGER, PARAMETER :: JP_O3 = 1
INTEGER, PARAMETER :: JP_H2O2 = 2
INTEGER, PARAMETER :: JP_NO = 3
INTEGER, PARAMETER :: JP_NO2 = 4
INTEGER, PARAMETER :: JP_NO3 = 5
INTEGER, PARAMETER :: JP_N2O5 = 6
INTEGER, PARAMETER :: JP_HONO = 7
INTEGER, PARAMETER :: JP_HNO3 = 8
INTEGER, PARAMETER :: JP_HNO4 = 9
INTEGER, PARAMETER :: JP_NH3 = 10
INTEGER, PARAMETER :: JP_SO2 = 11
INTEGER, PARAMETER :: JP_SULF = 12
INTEGER, PARAMETER :: JP_CO = 13
INTEGER, PARAMETER :: JP_OH = 14
INTEGER, PARAMETER :: JP_HO2 = 15
INTEGER, PARAMETER :: JP_CH4 = 16
INTEGER, PARAMETER :: JP_ETH = 17
INTEGER, PARAMETER :: JP_ALKA = 18
INTEGER, PARAMETER :: JP_ALKE = 19
INTEGER, PARAMETER :: JP_BIO = 20
INTEGER, PARAMETER :: JP_ARO = 21
INTEGER, PARAMETER :: JP_HCHO = 22
INTEGER, PARAMETER :: JP_ALD = 23
INTEGER, PARAMETER :: JP_KET = 24
INTEGER, PARAMETER :: JP_CARBO = 25
INTEGER, PARAMETER :: JP_ONIT = 26
INTEGER, PARAMETER :: JP_PAN = 27
INTEGER, PARAMETER :: JP_OP1 = 28
INTEGER, PARAMETER :: JP_OP2 = 29
INTEGER, PARAMETER :: JP_ORA2 = 30
INTEGER, PARAMETER :: JP_MO2 = 31
INTEGER, PARAMETER :: JP_ALKAP = 32
INTEGER, PARAMETER :: JP_ALKEP = 33
INTEGER, PARAMETER :: JP_BIOP = 34
INTEGER, PARAMETER :: JP_PHO = 35
INTEGER, PARAMETER :: JP_ADD = 36
INTEGER, PARAMETER :: JP_AROP = 37
INTEGER, PARAMETER :: JP_CARBOP = 38
INTEGER, PARAMETER :: JP_OLN = 39
INTEGER, PARAMETER :: JP_XO2 = 40
!
END MODULE MODD_CH_M9_SCHEME
