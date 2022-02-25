!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_SET_RSOU
!     ####################
!
INTERFACE
!
      SUBROUTINE SET_RSOU(TPFILE,TPEXPREFILE,HFUNU,HFUNV,KILOC,KJLOC,OBOUSS,&
                          PJ,OSHIFT,PCORIOZ)
!
USE MODD_IO, ONLY : TFILEDATA
!
TYPE(TFILEDATA),        INTENT(IN)  :: TPFILE ! outpput data file
TYPE(TFILEDATA),        INTENT(IN)  :: TPEXPREFILE ! input data file
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNU  ! type of variation of U
                                              ! in y direction
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNV  ! type of variation of V
                                              ! in x direction
INTEGER,                INTENT(IN)  :: KILOC  ! I Localisation of vertical profile
INTEGER,                INTENT(IN)  :: KJLOC  ! J Localisation of vertical profile
LOGICAL,                INTENT(IN)  :: OBOUSS ! logical switch for Boussinesq version
REAL, DIMENSION(:,:,:), INTENT(IN) :: PJ ! jacobien 
LOGICAL,                INTENT(IN)  :: OSHIFT ! logical switch for vertical shift
!
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: PCORIOZ ! Coriolis parameter
                                               ! (exceptionnaly 3D array)
!
END SUBROUTINE SET_RSOU
!
END INTERFACE
!
END MODULE MODI_SET_RSOU
!
!     ########################################################################
      SUBROUTINE SET_RSOU(TPFILE,TPEXPREFILE,HFUNU,HFUNV,KILOC,KJLOC,OBOUSS, &
                          PJ,OSHIFT,PCORIOZ)
!     ########################################################################
!
!!****  *SET_RSOU * -  to initialize mass fiels from a radiosounding 
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is  to initialize the mass field (theta,r,
!     thetavrefz,rhorefz) on model grid from a radiosounding located at point 
!     (KILOC,KJLOC).
!
!        The free-formatted part of EXPRE file contains the radiosounding data.The data 
!     are stored in following order :
!        
!         - year,month,day, time (these variables are read in PREINIT program)
!         - kind of data in EXPRE file (see below for more explanations about
!                                       YKIND)
!         - ZGROUND
!         - PGROUND
!         - temperature variable at ground ( depending on the data Kind )
!         - moist variable  at ground ( depending on the data Kind )
!         - number of wind data levels  ( variable ILEVELU)
!         - height  ,  dd  , ff   |
!            or           or      |  ILEVELU times
!           pressure,  U   , V    |
!         - number of mass levels  ( variable ILEVELM), including the ground
!           level
!         - height  ,  T  ,          Td              |
!            or        or            or              |  (ILEVELM-1) times
!           pressure,  THeta_Dry  , Mixing Ratio     |
!                      or            or              |
!                      THeta_V    , relative HUmidity|  
!
!       NB : the first mass level is at ground
! 
!       The following kind of data  is permitted :
!                  YKIND = 'STANDARD'  :   ZGROUND, PGROUND, TGROUND, TDGROUND
!                                         (Pressure, dd, ff) , 
!                                         (Pressure, T, Td)
!                  YKIND = 'PUVTHVMR'  : zGROUND, PGROUND, ThvGROUND, RGROUND
!                                        (Pressure, U, V) , 
!                                        (Pressure, THv, R)
!                  YKIND = 'PUVTHVHU'  :  zGROUND, PGROUND, ThvGROUND, HuGROUND
!                                         (Pressure, U, V) , 
!                                         (Pressure, THv, Hu)
!                  YKIND = 'ZUVTHVHU'  :  zGROUND, PGROUND, ThvGROUND, HuGROUND
!                                         (height, U, V) , 
!                                         (height, THv, Hu)
!                  YKIND = 'ZUVTHVMR'  :  zGROUND, PGROUND, ThvGROUND, RGROUND
!                                         (height, U, V) , 
!                                         (height, THv, R)
!                  YKIND = 'PUVTHDMR'  : zGROUND, PGROUND, ThdGROUND, RGROUND
!                                         (Pressure, U, V) , 
!                                         (Pressure, THd, R)
!                  YKIND = 'PUVTHDHU'  : zGROUND, PGROUND, ThdGROUND, HuGROUND
!                                         (Pressure, U, V) , 
!                                         (Pressure, THd, Hu)
!                  YKIND = 'ZUVTHDMR'  :  zGROUND, PGROUND, ThdGROUND,
!                  RGROUND
!                                         (height, U, V) , 
!                                         (height, THd, R)
!                  YKIND = 'PUVTHU'  :   ZGROUND, PGROUND, TGROUND, HuGROUND
!                                         (Pressure, U, V) , 
!                                         (Pressure, T, Hu)
!
!  For ocean-LES case  the following kind of data is permitted
!         
!                  YKIND = 'IDEALOCE' : ZGROUND (Water depth),PGROUND(Sfc Atmos Press),
!                                        TGROUND (SST), RGROUND (SSS)
!                                         (Depth , U, V) starting from sfc   
!                                         (Depth,  T, S)
!                                     (Time, LE, H, SW_d,SW_u,LW_d,LW_u,Stress_X,Stress_Y)
!
!                  YKIND = 'STANDOCE' :  (Depth , Temp, Salinity, U, V) starting from sfc   
!                                        (Time, LE, H, SW_d,SW_u,LW_d,LW_u,Stress_X,Stress_Y)     
!        
!!**  METHOD
!!    ------
!!      The radiosounding is first read, then data are converted in order to
!!    always obtain  the following variables  (case YKIND = 'ZUVTHVMR') :
!!       (height,U,V) and (height,Thetav,r) which are the model variables.
!!      That is to say : 
!!         - YKIND = 'STANDARD' : 
!!           dd,ff converted in U,V
!!           Td + pressure ----> r
!!           T,r ---> Tv  + pressure ----> thetav
!!           Pressure + thetav + ZGROUND ----> height (for mass levels)
!!           Thetav at mass levels ----> thetav at wind levels
!!           Pressure + thetav + ZGROUND + PGROUND ---->height (for wind levels)
!!         - YKIND = 'PUVTHVMR' : 
!!           Pressure + thetav + ZGROUND ----> height (for mass levels)
!!           Thetav at mass levels ----> thetav at wind levels
!!           Pressure + thetav + ZGROUND + PGROUND ---->height (for wind levels)
!!         - YKIND = 'PUVTHVHU' : 
!!           thetav + pressure ----> Tv +pressure +Hu ----> r 
!!           Pressure + thetav + ZGROUND ----> height (for mass levels)
!!           Thetav at mass levels ----> thetav at wind levels
!!           Pressure + thetav + ZGROUND + PGROUND ---->height (for wind levels)
!!         - YKIND = 'ZUVTHVHU' : 
!!           height +thetav + PGROUND -----> pressure (for mass levels)
!!           thetav + pressure ----> Tv +pressure +Hu ----> r 
!!         - YKIND = 'PUVTHDVMR' : 
!!           thetad + r ----> thetav
!!           pressure + thetav + ZGROUND ----> height (for mass levels)
!!           Thetav at mass levels ----> thetav at wind levels
!!           Pressure + thetav + ZGROUND + PGROUND ---->height (for wind levels)
!!         - YKIND = 'PUVTHDHU' :  
!!           thetad + pressure -----> T
!!           T + pressure + Hu -----> r 
!!           thetad + r -----> thetav
!!           pressure + thetav + ZGROUND ----> height (for mass levels)
!!           Thetav at mass levels ----> thetav at wind levels
!!           Pressure + thetav + ZGROUND + PGROUND ---->height (for wind levels)
!!         - YKIND = 'ZUVTHDHU' :  
!!           thetad + r -----> thetav
!!         - YKIND = 'PUVTHU' :  
!!           T + pressure -----> thetad
!!           T + pressure + Hu  -----> r
!!           thetad + r -----> thetav
!!           pressure + thetav + ZGROUND ----> height (for mass levels)
!!           Thetav at mass levels ----> thetav at wind levels
!!           
!!      The following basic formula are used :
!!                Rd  es(Td)
!!            r = --  ----------
!!                Rv  P - es(Td)
!!
!!                 1 + (Rv/Rd) r
!!            Tv = --------------  T  
!!                    1 + r
!!
!!                          P00    Rd/Cpd             1 + (Rv/Rd) r
!!            Thetav = Tv ( ---- )         = Thetad ( --------------)
!!                           P                           1 + r
!!       The integration of hydrostatic relation is used to compute height from
!!   pressure  and vice-versa. This is done by HEIGHT_PRESS and PRESS_HEIGHT 
!!   routines.  
!!   
!!      Then, these data are interpolated on a vertical grid which is  
!!    a mixed grid calaculated with VERT_COORD from the vertical levels of MNH
!!    grid and with a constant ororgraphy equal to the altitude of the vertical
!!    profile (ZZGROUND) (It permits to keep low levels information with a
!!    shifting function (as in PREP_REAL_CASE))
!!
!!      Then, the 3D mass and wind fields are deduced in SET_MASS
!!     
!!
!!    EXTERNAL
!!    --------
!!      SET_MASS : to compute mass field on 3D-model grid
!!      Module MODE_THERMO : contains thermodynamic routines
!!         SM_FOES : To compute saturation vapor pressure from 
!!                   temperature
!!         SM_PMR_HU : to compute vapor mixing ratio from pressure, virtual
!!                    temperature and relative humidity
!!      HEIGHT_PRESS : to compute height from pressure and thetav 
!!                    by integration of hydrostatic relation
!!      PRESS_HEIGHT : to compute pressure from height and thetav
!!                    by integration of hydrostatic relation
!!      THETAVPU_THETAVPM : to interpolate thetav on wind levels
!!                          from thetav on mass levels
!!
!!      Module MODI_HEIGHT_PRESS      : interface for function HEIGHT_PRESS
!!      Module MODI_PRESS_HEIGHT      : interface for function PRESS_HEIGHT
!!      Module MODI_THETAVPU_THETAVPM : interface for function
!!                                      THETAVPU_THETVPM
!!      Module MODI_SET_MASS          : interface for subroutine SET_MASS
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST     : contains physical constants
!!        XPI  : Pi
!!        XRV  : Gas constant for vapor 
!!        XRD  : Gas constant for dry air
!!        XCPD : Specific heat for dry air at constant pressure
!!
!!      Module MODD_LUNIT1  : contains logical unit names 
!!        TLUOUT : name of output-listing
!!
!!      Module MODD_CONF    : contains configuration variables for all models. 
!!        NVERB : verbosity level for output-listing
!!
!!      Module MODD_GRID1  : contains grid variables
!!        XZHAT : height of w-levels of vertical model grid without orography
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (routine SET_RSOU)
!!    
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/08/94
!!      J.Stein     06/12/94  change the way to prescribe the horizontal wind
!!                            variations + cleaning
!!      J.Stein     18/01/95  bug corrections in the ILEVELM readings
!!      J.Stein     16/04/95  put the same names of the declarative modules
!!                            in the descriptive part
!!      J.Stein     30/01/96  use the RS ground pressure to initialize the
!!                            hydrostatic pressure computation
!!      V.Masson    02/09/96  add allocation of ZTHVU in two cases
!!      P.Jabouille 14/02/96  bug in extrapolation of ZMRM below the first level
!! Jabouille/Masson 05/12/02  add  ZUVTHLMR case and hydrometeor initialization
!!      P.Jabouille 29/10/03  add  hydrometeor initialization for ZUVTHDMR case
!!      G. Tanguy   26/10/10  change the interpolation of the RS : we use now a
!!                            mixed grid (PREP_REAL_CASE method)
!!                            add PUVTHU case
!!      V.Masson    12/08/13  Parallelization of the initilization profile
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 19/04/2019: removed unused dummy arguments and variables
!  JL Redelsperger 01/2021: Ocean LES cases added
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST
USE MODD_DYN_n,      ONLY: LOCEAN
USE MODD_FIELD_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_LUNIT_n,    ONLY: TLUOUT
USE MODD_IO,         ONLY: TFILEDATA
USE MODD_NETCDF
USE MODD_OCEANH
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODD_TYPE_DATE
!
USE MODE_ll
USE MODE_MSG
USE MODE_THERMO
!
USE MODI_COMPUTE_EXNER_FROM_GROUND
USE MODI_HEIGHT_PRESS
USE MODI_PRESS_HEIGHT
USE MODI_SET_MASS
USE MODI_SHUMAN
USE MODI_THETAVPU_THETAVPM
USE MODI_TH_R_FROM_THL_RT_1D
USE MODI_VERT_COORD
!
USE NETCDF          ! for reading the NR files 
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments :
!
TYPE(TFILEDATA),        INTENT(IN)  :: TPFILE ! outpput data file
TYPE(TFILEDATA),        INTENT(IN)  :: TPEXPREFILE ! input data file
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNU  ! type of variation of U
                                              ! in y direction
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNV  ! type of variation of V
                                              ! in x direction
INTEGER,                INTENT(IN)  :: KILOC  ! I Localisation of vertical profile
INTEGER,                INTENT(IN)  :: KJLOC  ! J Localisation of vertical profile
LOGICAL,                INTENT(IN)  :: OBOUSS ! logical switch for Boussinesq version
LOGICAL,                INTENT(IN)  :: OSHIFT ! logical switch for vertical shift
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: PCORIOZ ! Coriolis parameter
                                               ! (exceptionnaly 3D array)
REAL, DIMENSION(:,:,:), INTENT(IN) :: PJ ! jacobien 
!
!
!*       0.2   Declarations of local variables :
!
INTEGER                         :: ILUPRE ! logical unit number of the EXPRE return code
INTEGER                         :: ILUOUT    ! Logical unit number for output-listing
! local variables for reading sea sfc flux forcing for ocean model
INTEGER                   :: IFRCLT 
REAL, DIMENSION(:), ALLOCATABLE :: ZSSUFL_T,ZSSVFL_T,ZSSTFL_T,ZSSOLA_T !
TYPE (DATE_TIME), DIMENSION(:), ALLOCATABLE :: ZFRCLT ! date/time of sea surface forcings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  variables read in EXPRE file at the RS/CTD levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHARACTER(LEN=8)                :: YKIND     ! Kind of  variables in 
                                             ! EXPRE FILE
INTEGER                         :: ILEVELU   ! number of wind levels 
REAL, DIMENSION(:), ALLOCATABLE :: ZHEIGHTU  ! Height at wind levels 
REAL, DIMENSION(:), ALLOCATABLE :: ZPRESSU   ! Pressure at wind levels
REAL, DIMENSION(:), ALLOCATABLE :: ZTHVU     ! Thetav at wind levels
REAL, DIMENSION(:), ALLOCATABLE :: ZU,ZV     ! wind components
REAL, DIMENSION(:), ALLOCATABLE :: ZDD,ZFF   ! dd (direction) and ff(force)
                                             !     for wind
REAL                            :: ZZGROUND,ZPGROUND ! height and Pressure at ground  
REAL                            :: ZTGROUND,ZTHVGROUND,ZTHDGROUND,ZTHLGROUND,    & 
                                   ZTDGROUND,ZMRGROUND,ZHUGROUND        
                                                  ! temperature and moisture
                                                  ! variables at ground
INTEGER                         :: ILEVELM   ! number of mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZHEIGHTM  ! Height at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZPRESSM   ! Pressure at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZTHV      ! Thetav at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZTHD      ! Theta (dry) at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZTHL      ! Thetal at mass levels  
REAL, DIMENSION(:), ALLOCATABLE :: ZTH       ! Theta at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZT        ! Temperature at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZMR       ! Vapor mixing ratio at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZMRC      ! cloud mixing ratio at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZMRI      ! ice mixing ratio or cloud concentration
REAL, DIMENSION(:), ALLOCATABLE :: ZRT       ! total mixing ratio  
REAL, DIMENSION(:), ALLOCATABLE :: ZPRESS    ! pressure at mass level 
REAL, DIMENSION(:), ALLOCATABLE :: ZHU       ! relative humidity at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZTD       ! Td at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZTV       ! Tv at mass levels
REAL, DIMENSION(:), ALLOCATABLE :: ZEXN      
REAL, DIMENSION(:), ALLOCATABLE :: ZCPH     
REAL, DIMENSION(:), ALLOCATABLE :: ZLVOCPEXN     
REAL, DIMENSION(:), ALLOCATABLE :: ZLSOCPEXN    
REAL, DIMENSION(SIZE(XZHAT))    :: ZZFLUX_PROFILE ! altitude of flux points on the initialization columns
REAL, DIMENSION(SIZE(XZHAT))    :: ZZMASS_PROFILE ! altitude of mass points on the initialization columns
!
!  fields on the grid of the model without orography
!
REAL, DIMENSION(SIZE(XZHAT))    :: ZUW,ZVW ! Wind at w model grid levels
REAL, DIMENSION(SIZE(XZHAT))    :: ZMRM   ! vapor mixing ratio at mass model
                                          !grid levels 
REAL, DIMENSION(SIZE(XZHAT))    :: ZMRCM,ZMRIM
REAL, DIMENSION(SIZE(XZHAT))    :: ZTHVM  ! Temperature at mass model grid levels
REAL, DIMENSION(SIZE(XZHAT))    :: ZTHLM  ! Thetal at mass model grid levels
REAL, DIMENSION(SIZE(XZHAT))    :: ZTHM  ! Thetal at mass model grid levels
REAL, DIMENSION(SIZE(XZHAT))    :: ZRHODM   ! density at mass model grid level
REAL, DIMENSION(:), ALLOCATABLE :: ZMRT    ! Total Vapor mixing ratio at mass levels on mixed grid
REAL, DIMENSION(:), ALLOCATABLE :: ZEXNMASS  ! exner fonction at mass level
REAL, DIMENSION(:), ALLOCATABLE :: ZEXNFLUX  ! exner fonction at flux level
REAL                            :: ZEXNSURF  ! exner fonction at surface
REAL, DIMENSION(:), ALLOCATABLE :: ZPREFLUX   ! pressure at flux model grid level
REAL, DIMENSION(:), ALLOCATABLE :: ZFRAC_ICE ! ice fraction
REAL, DIMENSION(:), ALLOCATABLE :: ZRSATW, ZRSATI             
REAL                            :: ZDZSDH,ZDZ1SDH,ZDZ2SDH ! interpolation
                                                          ! working arrays
!
INTEGER         :: JK,JKLEV,JKU,JKM,JKT,JJ,JI,JO,JLOOP  ! Loop indexes
INTEGER         :: IKU                ! Upper bound in z direction
REAL            :: ZRDSCPD,ZRADSDG, & ! Rd/Cpd, Pi/180.,
                   ZRVSRD,ZRDSRV,   & ! Rv/Rd, Rd/Rv
                   ZPTOP              ! Pressure at domain top 
LOGICAL         :: GUSERC             ! use of input data cloud
INTEGER         :: IIB, IIE, IJB, IJE
INTEGER         :: IXOR_ll, IYOR_ll
INTEGER         :: IINFO_ll
LOGICAL         :: GPROFILE_IN_PROC   ! T : initialization profile is in current processor
!
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT))   ::ZZS_LS
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) ::ZZFLUX_MX,ZZMASS_MX ! mixed grid
!-------------------------------------------------------------------------------
! For standard ocean version, reading external files
CHARACTER(LEN=256) :: yinfile, yinfisf ! files to be read
INTEGER :: IDX
INTEGER(KIND=CDFINT) :: INZ, INLATI, INLONGI
INTEGER(KIND=CDFINT) :: incid, ivarid, idimid, idimlen
REAL, DIMENSION(:,:,:),     ALLOCATABLE :: ZOC_TEMPERATURE,ZOC_SALINITY,ZOC_U,ZOC_V
REAL, DIMENSION(:),     ALLOCATABLE :: ZOC_DEPTH  
REAL, DIMENSION(:),     ALLOCATABLE :: ZOC_LE,ZOC_H
REAL, DIMENSION(:),     ALLOCATABLE :: ZOC_SW_DOWN,ZOC_SW_UP,ZOC_LW_DOWN,ZOC_LW_UP
REAL, DIMENSION(:),     ALLOCATABLE :: ZOC_TAUX,ZOC_TAUY

!--------------------------------------------------------------------------------
!
!*	 1.     PROLOGUE : INITIALIZE SOME CONSTANTS, RETRIEVE LOGICAL
!               UNIT NUMBERS AND READ KIND OF DATA IN EXPRE FILE
!	        -------------------------------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
CALL GET_OR_ll('B',IXOR_ll,IYOR_ll)
!
!*       1.1  initialize some constants 
!
ZRDSCPD = XRD / XCPD
ZRADSDG = XPI/180.
ZRVSRD  = XRV/XRD
ZRDSRV = XRD/XRV
!
!*       1.2  Retrieve logical unit numbers 
!                         
ILUPRE = TPEXPREFILE%NLU
ILUOUT = TLUOUT%NLU
!
!*       1.3  Read data kind in EXPRE file 
!
READ(ILUPRE,*) YKIND
WRITE(ILUOUT,*) 'YKIND read in set_rsou: ', YKIND
!
IF(LUSERC .AND. YKIND/='PUVTHDMR' .AND. YKIND/='ZUVTHDMR' .AND.  YKIND/='ZUVTHLMR') THEN 
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_RSOU','hydrometeors are not allowed for YKIND = '//trim(YKIND))
ENDIF
!
IF(YKIND=='ZUVTHLMR' .AND. .NOT. LUSERC) THEN
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_RSOU','LUSERC=T is required for YKIND=ZUVTHLMR')
ENDIF
!
GUSERC=.FALSE.
IF(LUSERC .AND. (YKIND == 'PUVTHDMR' .OR. YKIND == 'ZUVTHDMR')) GUSERC=.TRUE.
!-------------------------------------------------------------------------------
!
!*	 2.     READ DATA AND CONVERT IN (height,U,V), (height,Thetav,r)
!	        --------------------------------------------------------
!
SELECT CASE(YKIND)
!   
!     2.0.1 Ocean case 1
!   
  CASE ('IDEALOCE')
!
    XP00=XP00OCEAN
    ! Read data in PRE_IDEA1.nam
    ! Surface      
    WRITE(ILUOUT,FMT=*) 'Reading data for ideal ocean :IDEALOCE'   
    READ(ILUPRE,*) ZPTOP           ! P_atmosphere at sfc =P top domain
    READ(ILUPRE,*) ZTGROUND        ! SST 
    READ(ILUPRE,*) ZMRGROUND       ! SSS
    WRITE(ILUOUT,FMT=*) 'Patm SST SSS', ZPTOP,ZTGROUND,ZMRGROUND
    READ(ILUPRE,*) ILEVELU         ! Read number of Current levels
    ! Allocate required memory 
    ALLOCATE(ZHEIGHTU(ILEVELU),ZU(ILEVELU),ZV(ILEVELU))
    ALLOCATE(ZOC_U(ILEVELU,1,1),ZOC_V(ILEVELU,1,1))
    WRITE(ILUOUT,FMT=*) 'Level number for Current in data', ILEVELU
    ! Read U and V at each wind level
    DO JKU = 1,ILEVELU
      READ(ILUPRE,*) ZHEIGHTU(JKU),ZOC_U(JKU,1,1),ZOC_V(JKU,1,1)
      ! WRITE(ILUOUT,FMT=*) 'Leveldata D(m) under sfc: U_cur, V_cur', JKU, ZHEIGHTU(JKU),ZU(JKU),ZV(JKU)
    END DO
    DO JKU=1,ILEVELU
      ! Z axis reoriented as in the model
      IDX     = ILEVELU-JKU+1
      ZU(JKU) = ZOC_U(IDX,1,1)
      ZV(JKU) = ZOC_V(IDX,1,1)
      ! ZHEIGHT used only in set_ rsou, defined as such ZHEIGHT(ILEVELM)=H_model
      ! Z oriented in same time to have a model domain axis going
      ! from 0m (ocean bottom/model bottom) towards H (ocean sfc/model top)
    END DO
    ! Read number of mass levels  
    READ(ILUPRE,*) ILEVELM
    ! Allocate required memory 
    ALLOCATE(ZOC_DEPTH(ILEVELM))
    ALLOCATE(ZHEIGHTM(ILEVELM))           
    ALLOCATE(ZTHL(ILEVELM),ZTH(ILEVELM),ZTHV(ILEVELM)) 
    ALLOCATE(ZMR(ILEVELM),ZRT(ILEVELM))
    ALLOCATE(ZOC_TEMPERATURE(ILEVELM,1,1),ZOC_SALINITY(ILEVELM,1,1))
    ! Read T and S at each mass level 
    DO JKM= 2,ILEVELM
      READ(ILUPRE,*) ZOC_DEPTH(JKM),ZOC_TEMPERATURE(JKM,1,1),ZOC_SALINITY(JKM,1,1)
    END DO
    ! Complete the mass arrays with the ground informations read in EXPRE file
    ZOC_DEPTH(1)    = 0.                    
    ZOC_TEMPERATURE(1,1,1)= ZTGROUND
    ZOC_SALINITY(1,1,1)= ZMRGROUND
    !!!!!!!!!!!!!!!!!!!!!!!!Inversing Axis!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Going from the data (axis downward i.e inverse model) grid to the model grid (axis upward)
    ! Uniform bathymetry; depth goes from ocean sfc downwards  (data grid)
    ! ZHEIGHT goes from the model domain bottom up to the sfc ocean (top of model domain)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ZZGROUND   = 0.
    ZTGROUND   = ZOC_TEMPERATURE(ILEVELM,1,1) 
    ZMRGROUND  = ZOC_SALINITY(ILEVELM,1,1)
    DO JKM= 1,ILEVELM
      ! Z upward axis (oriented as in the model), i.e.
      ! going from 0m (ocean bottom/model bottom) upward to H (ocean sfc/model top)
      ! ZHEIGHT used only in set_ rsou, defined as such ZHEIGHT(ILEVELM)=H_model
      IDX          = ILEVELM-JKM+1
      ZTH(JKM)      = ZOC_TEMPERATURE(IDX,1,1)
      ZMR(JKM)     = ZOC_SALINITY(IDX,1,1)
      ZHEIGHTM(JKM)= ZOC_DEPTH(ILEVELM)- ZOC_DEPTH(IDX)
      WRITE(ILUOUT,FMT=*) 'Model oriented initial data: JKM IDX depth T S ZHEIGHTM', &
                      JKM,IDX,ZOC_DEPTH(IDX),ZTH(JKM),ZMR(JKM),ZHEIGHTM(JKM)
    END DO
    ! mass levels of the RS
    ZTHV = ZTH ! TV==THETA=TL
    ZTHL = ZTH
    ZRT  = ZMR  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! READ Sea Surface Forcing !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Reading the forcings from prep_idea1.nam
    READ(ILUPRE,*) IFRCLT  ! Number of time-dependent forcing 
    IF (IFRCLT > 99*8) THEN
      ! CAUTION: number of forcing times is limited by the WRITE format 99(8E10.3)
      !          and also by the name of forcing variables (format I3.3)
      !          You have to modify those if you need more forcing times
     CALL PRINT_MSG(NVERB_FATAL,'IO','SET_RSOU','maximum forcing times NFRCLT is 99*8')
    END IF
!
    WRITE(UNIT=ILUOUT,FMT='(" THERE ARE ",I2," SFC FLUX FORCINGs AT:")') IFRCLT
    ALLOCATE(ZFRCLT(IFRCLT))
    ALLOCATE(ZSSUFL_T(IFRCLT)); ZSSUFL_T = 0.0
    ALLOCATE(ZSSVFL_T(IFRCLT)); ZSSVFL_T = 0.0
    ALLOCATE(ZSSTFL_T(IFRCLT)); ZSSTFL_T = 0.0
    ALLOCATE(ZSSOLA_T(IFRCLT)); ZSSOLA_T = 0.0
    DO JKT = 1,IFRCLT
      WRITE(ILUOUT,FMT='(A, I4)') "SET_RSOU/Reading Sea Surface forcing: Number=", JKT
      READ(ILUPRE,*) ZFRCLT(JKT)%nyear, ZFRCLT(JKT)%nmonth, &
                     ZFRCLT(JKT)%nday,  ZFRCLT(JKT)%xtime
      READ(ILUPRE,*) ZSSUFL_T(JKT)
      READ(ILUPRE,*) ZSSVFL_T(JKT)
      READ(ILUPRE,*) ZSSTFL_T(JKT)
      READ(ILUPRE,*) ZSSOLA_T(JKT)
    END DO
!
    DO JKT = 1 , IFRCLT
      WRITE(UNIT=ILUOUT,FMT='(F9.0, "s, date:", I3, "/", I3, "/", I5)') &
                 ZFRCLT(JKT)%xtime,         ZFRCLT(JKT)%nday,   &
                 ZFRCLT(JKT)%nmonth, ZFRCLT(JKT)%nyear
    END DO
    NINFRT= INT(ZFRCLT(2)%xtime)
    WRITE(ILUOUT,FMT='(A)') &
         "Number U-Stress, V-Stress, Heat turb Flux, Solar Flux Interval(s)",NINFRT
    DO JKT = 1, IFRCLT
      WRITE(ILUOUT,FMT='(I10,99(3F10.2))') JKT, ZSSUFL_T(JKT),ZSSVFL_T(JKT),ZSSTFL_T(JKT) 
    END DO
    NFRCLT = IFRCLT
    ALLOCATE(TFRCLT(NFRCLT))
    ALLOCATE(XSSUFL_T(NFRCLT));XSSUFL_T(:)=0.
    ALLOCATE(XSSVFL_T(NFRCLT));XSSVFL_T(:)=0.
    ALLOCATE(XSSTFL_T(NFRCLT));XSSTFL_T(:)=0.
    ALLOCATE(XSSOLA_T(NFRCLT));XSSOLA_T(:)=0.
!
    DO JKT=1,NFRCLT
      TFRCLT(JKT)= ZFRCLT(JKT)
      XSSUFL_T(JKT)=ZSSUFL_T(JKT)/XRH00OCEAN
      XSSVFL_T(JKT)=ZSSVFL_T(JKT)/XRH00OCEAN
      ! working in SI
      XSSTFL_T(JKT)=ZSSTFL_T(JKT) /(3900.*XRH00OCEAN)
      XSSOLA_T(JKT)=ZSSOLA_T(JKT) /(3900.*XRH00OCEAN)
    END DO   
    DEALLOCATE(ZFRCLT)
    DEALLOCATE(ZSSUFL_T)
    DEALLOCATE(ZSSVFL_T)
    DEALLOCATE(ZSSTFL_T)
    DEALLOCATE(ZSSOLA_T)
!
!--------------------------------------------------------------------------------   
! 2.0.2  Ocean standard initialize from netcdf files
!        U,V,T,S at Z levels + Forcings at model TOP (sea surface) 
!--------------------------------------------------------------------------------   
!
  CASE ('STANDOCE')
!   
    XP00=XP00OCEAN
    READ(ILUPRE,*) ZPTOP           ! P_atmosphere at sfc =P top domain
    READ(ILUPRE,*) YINFILE, YINFISF
    WRITE(ILUOUT,FMT=*) 'Netcdf files to read:', YINFILE, YINFISF
    ! Open file containing initial profiles
    CALL check(nf90_open(yinfile,NF90_NOWRITE,incid), "opening NC file")
    ! Reading dimensions and lengths
    CALL check( nf90_inq_dimid(incid, "depth",idimid), "getting depth  dimension id" )
    CALL check( nf90_inquire_dimension(incid, idimid,             len=INZ),     "getting INZ"   )
    CALL check( nf90_inquire_dimension(incid, INT(2,KIND=CDFINT), len=INLONGI), "getting NLONG" )
    CALL check( nf90_inquire_dimension(incid, INT(1,KIND=CDFINT), len=INLATI),  "getting NLAT"  )
!   
    WRITE(ILUOUT,FMT=*) 'NB LEVLS READ INZ, NLONG NLAT ', INZ, INLONGI,INLATI
    ALLOCATE(ZOC_TEMPERATURE(INLATI,INLONGI,INZ),ZOC_SALINITY(INLATI,INLONGI,INZ))
    ALLOCATE(ZOC_U(INLATI,INLONGI,INZ),ZOC_V(INLATI,INLONGI,INZ))
    ALLOCATE(ZOC_DEPTH(INZ))
    WRITE(ILUOUT,FMT=*) 'NETCDF READING ==> Temp'
    CALL check(nf90_inq_varid(incid,"temperature",ivarid), "getting temp ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_TEMPERATURE), "reading temp")
    WRITE(ILUOUT,FMT=*) 'Netcdf Reading ==> salinity'
    CALL check(nf90_inq_varid(incid,"salinity",ivarid), "getting salinity ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_SALINITY), "reading salinity")
    WRITE(ILUOUT,FMT=*) 'Netcdf ==> Reading depth'
    CALL check(nf90_inq_varid(incid,"depth",ivarid), "getting depth ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_DEPTH), "reading depth")
    WRITE(ILUOUT,FMT=*) 'depth: max min ', MAXVAL(ZOC_DEPTH),MINVAL(ZOC_DEPTH)
    WRITE(ILUOUT,FMT=*) 'depth 1 nz: ', ZOC_DEPTH(1),ZOC_DEPTH(INZ)
    WRITE(ILUOUT,FMT=*) 'Netcdf Reading ==> Currents'
    CALL check(nf90_inq_varid(incid,"u",ivarid), "getting u ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_U), "reading u")
    CALL check(nf90_inq_varid(incid,"v",ivarid), "getting v ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_V), "reading v")
    CALL check(nf90_close(incid), "closing yinfile")
    WRITE(ILUOUT,FMT=*) 'End of initial file reading'
!
    DO JKM=1,INZ
     ZOC_TEMPERATURE(1,1,JKM)=ZOC_TEMPERATURE(1,1,JKM)+273.15
     WRITE(ILUOUT,FMT=*) 'Z T(Kelvin) S(Sverdup) U V K',&
     JKM,ZOC_DEPTH(JKM),ZOC_TEMPERATURE(1,1,JKM),ZOC_SALINITY(1,1,JKM),ZOC_U(1,1,JKM),ZOC_V(1,1,JKM), JKM
    ENDDO
    !  number of data levels
    ILEVELM=INZ
    ! Model bottom
    ZTGROUND  = ZOC_TEMPERATURE(1,1,ILEVELM)
    ZMRGROUND = ZOC_SALINITY(1,1,ILEVELM)
    ZZGROUND=0.
    ! Allocate required memory
    ALLOCATE(ZHEIGHTM(ILEVELM))           
    ALLOCATE(ZT(ILEVELM))
    ALLOCATE(ZTV(ILEVELM))
    ALLOCATE(ZMR(ILEVELM))
    ALLOCATE(ZTHV(ILEVELM)) 
    ALLOCATE(ZTHL(ILEVELM))
    ALLOCATE(ZRT(ILEVELM))
    !  Going from the inverse model grid (data) to the normal one
    DO JKM= 1,ILEVELM
      ! Z axis reoriented as in the model
      IDX     = ILEVELM-JKM+1
      ZT(JKM) = ZOC_TEMPERATURE(1,1,IDX)
      ZMR(JKM)  =  ZOC_SALINITY(1,1,IDX)
      ! ZHEIGHT used only in set_ rsou, defined as such ZHEIGHT(ILEVELM)=H_model
      ! Z oriented in same time to have a model domain axis going
      ! from 0m (ocean bottom/model bottom) towards H (ocean sfc/model top)
      ! translation/inversion
      ZHEIGHTM(JKM) = -ZOC_DEPTH(IDX) + ZOC_DEPTH(ILEVELM)
      WRITE(ILUOUT,FMT=*) 'End gridmodel comput: JKM IDX depth T S ZHEIGHTM', &
      JKM,IDX,ZOC_DEPTH(IDX),ZT(JKM),ZMR(JKM),ZHEIGHTM(JKM)
    END DO
    ! complete ther variables
    ZTV  = ZT
    ZTHV = ZT
    ZRT  = ZMR
    ZTHL = ZT
    ZTH  = ZT
    ! INIT --- U V -----
    ILEVELU = INZ               ! Same nb of levels for u,v,T,S
    !Assume that current and temp are given at same level
    ALLOCATE(ZHEIGHTU(ILEVELU))           
    ALLOCATE(ZU(ILEVELU),ZV(ILEVELU))
    ZHEIGHTU=ZHEIGHTM
    DO JKM= 1,ILEVELU
      ! Z axis reoriented as in the model
      IDX     = ILEVELU-JKM+1
      ZU(JKM) = ZOC_U(1,1,IDX)
      ZV(JKM) = ZOC_V(1,1,IDX)
      ! ZHEIGHT used only in set_ rsou, defined as such ZHEIGHT(ILEVELM)=H_model
      ! Z oriented in same time to have a model domain axis going
      ! from 0m (ocean bottom/model bottom) towards H (ocean sfc/model top)
    END DO
!
    DEALLOCATE(ZOC_TEMPERATURE)
    DEALLOCATE(ZOC_SALINITY)
    DEALLOCATE(ZOC_U)
    DEALLOCATE(ZOC_V)
    DEALLOCATE(ZOC_DEPTH)
!
    ! Reading/initializing  surface forcings
!
    WRITE(ILUOUT,FMT=*) 'netcdf sfc forcings file to be read:',yinfisf
    ! Open of sfc forcing file
    CALL check(nf90_open(yinfisf,NF90_NOWRITE,incid), "opening NC file")
    ! Reading dimension and length
    CALL check( nf90_inq_dimid(incid,"t",idimid), "getting  time dimension id" )
    CALL check( nf90_inquire_dimension(incid, idimid, len=idimlen), "getting idimlen "  )
!
    WRITE(ILUOUT,FMT=*) 'nb sfc-forcing time idimlen=',idimlen
    ALLOCATE(ZOC_LE(idimlen))
    ALLOCATE(ZOC_H(idimlen))
    ALLOCATE(ZOC_SW_DOWN(idimlen))
    ALLOCATE(ZOC_SW_UP(idimlen))
    ALLOCATE(ZOC_LW_DOWN(idimlen))
    ALLOCATE(ZOC_LW_UP(idimlen))
    ALLOCATE(ZOC_TAUX(idimlen))
    ALLOCATE(ZOC_TAUY(idimlen))
!
    WRITE(ILUOUT,FMT=*)'Netcdf Reading ==> LE'  
    CALL check(nf90_inq_varid(incid,"LE",ivarid), "getting LE ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_LE), "reading LE flux")
    WRITE(ILUOUT,FMT=*)'Netcdf Reading ==> H'  
    CALL check(nf90_inq_varid(incid,"H",ivarid), "getting H ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_H), "reading H flux")
    WRITE(ILUOUT,FMT=*) 'Netcdf Reading ==> SW_DOWN'  
    CALL check(nf90_inq_varid(incid,"SW_DOWN",ivarid), "getting SW_DOWN ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_SW_DOWN), "reading SW_DOWN")
    WRITE(ILUOUT,FMT=*) 'Netcdf Reading ==> SW_UP'  
    CALL check(nf90_inq_varid(incid,"SW_UP",ivarid), "getting SW_UP ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_SW_UP), "reading SW_UP")
    WRITE(ILUOUT,FMT=*) 'Netcdf Reading ==> LW_DOWN'  
    CALL check(nf90_inq_varid(incid,"LW_DOWN",ivarid), "getting LW_DOWN ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_LW_DOWN), "reading LW_DOWN")
    WRITE(ILUOUT,FMT=*) 'Netcdf Reading ==> LW_UP'  
    CALL check(nf90_inq_varid(incid,"LW_UP",ivarid), "getting LW_UP ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_LW_UP), "reading LW_UP")
    WRITE(ILUOUT,FMT=*) 'Netcdf Reading ==> TAUX'  
    CALL check(nf90_inq_varid(incid,"TAUX",ivarid), "getting TAUX ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_TAUX), "reading TAUX")
    WRITE(ILUOUT,FMT=*) 'Netcdf Reading ==> TAUY'  
    CALL check(nf90_inq_varid(incid,"TAUY",ivarid), "getting TAUY ivarid")
    CALL check(nf90_get_var(incid,ivarid,ZOC_TAUY), "reading TAUY")
    CALL check(nf90_close(incid), "closing yinfifs")
!
    WRITE(ILUOUT,FMT=*) '  Forcing-Number    LE     H     SW_down     SW_up    LW_down   LW_up TauX TauY' 
    DO JKM = 1, idimlen
      WRITE(ILUOUT,FMT=*) JKM, ZOC_LE(JKM), ZOC_H(JKM),ZOC_SW_DOWN(JKM),ZOC_SW_UP(JKM),&
                          ZOC_LW_DOWN(JKM),ZOC_LW_UP(JKM),ZOC_TAUX(JKM),ZOC_TAUY(JKM)   
    ENDDO
    ! IFRCLT FORCINGS at sea surface
    IFRCLT=idimlen
    ALLOCATE(ZFRCLT(IFRCLT)) 
    ALLOCATE(ZSSUFL_T(IFRCLT)); ZSSUFL_T = 0.0
    ALLOCATE(ZSSVFL_T(IFRCLT)); ZSSVFL_T = 0.0
    ALLOCATE(ZSSTFL_T(IFRCLT)); ZSSTFL_T = 0.0
    ALLOCATE(ZSSOLA_T(IFRCLT)); ZSSOLA_T = 0.0
    DO JKT=1,IFRCLT
      ! Initial file for CINDY-DYNAMO: all fluxes correspond to the absolute value (>0)
      ! modele ocean: axe z dirigé du bas vers la sfc de l'océan
      ! => flux dirigé vers le haut (positif ocean vers l'atmopshere i.e. bas vers le haut)
      ZSSOLA_T(JKT)=ZOC_SW_DOWN(JKT)-ZOC_SW_UP(JKT)
      ZSSTFL_T(JKT)=(ZOC_LW_DOWN(JKT)-ZOC_LW_UP(JKT)-ZOC_LE(JKT)-ZOC_H(JKT))
      ! assume that Tau given on file is along Ox
      ! rho_air UW_air = rho_ocean UW_ocean= N/m2
      ! uw_ocean
      ZSSUFL_T(JKT)=ZOC_TAUX(JKT)
      ZSSVFL_T(JKT)=ZOC_TAUY(JKT)
      WRITE(ILUOUT,FMT=*) 'Forcing Nb Sol NSol UW_oc VW',&
                          JKT,ZSSOLA_T(JKT),ZSSTFL_T(JKT),ZSSUFL_T(JKT),ZSSVFL_T(JKT) 
    ENDDO
    ! Allocate and Writing the corresponding variables in module MODD_OCEAN_FRC
    NFRCLT=IFRCLT
    ! value to read later on file ? 
    NINFRT=600
    ALLOCATE(TFRCLT(NFRCLT))
    ALLOCATE(XSSUFL_T(NFRCLT));XSSUFL_T(:)=0.
    ALLOCATE(XSSVFL_T(NFRCLT));XSSVFL_T(:)=0.
    ALLOCATE(XSSTFL_T(NFRCLT));XSSTFL_T(:)=0.
    ALLOCATE(XSSOLA_T(NFRCLT));XSSOLA_T(:)=0.
    ! on passe en unités SI, signe, etc pour le modele ocean
    !  W/m2 => SI :  /(CP_mer * rho_mer)
    ! a revoir dans tt le code pour mettre de svaleurs plus exactes
    DO JKT=1,NFRCLT
      TFRCLT(JKT)= ZFRCLT(JKT)
      XSSUFL_T(JKT)=ZSSUFL_T(JKT)/XRH00OCEAN
      XSSVFL_T(JKT)=ZSSVFL_T(JKT)/XRH00OCEAN
      XSSTFL_T(JKT)=ZSSTFL_T(JKT) /(3900.*XRH00OCEAN)
      XSSOLA_T(JKT)=ZSSOLA_T(JKT) /(3900.*XRH00OCEAN)
    END DO   
    DEALLOCATE(ZFRCLT)
    DEALLOCATE(ZSSUFL_T)
    DEALLOCATE(ZSSVFL_T)
    DEALLOCATE(ZSSTFL_T)
    DEALLOCATE(ZSSOLA_T)
    DEALLOCATE(ZOC_LE)    
    DEALLOCATE(ZOC_H)    
    DEALLOCATE(ZOC_SW_DOWN)    
    DEALLOCATE(ZOC_SW_UP)    
    DEALLOCATE(ZOC_LW_DOWN)    
    DEALLOCATE(ZOC_LW_UP)    
    DEALLOCATE(ZOC_TAUX)
    DEALLOCATE(ZOC_TAUY)
    ! END OCEAN STANDARD
!
!
!*       2.1  ATMOSPHERIC STANDARD case : ZGROUND, PGROUND, TGROUND, TDGROUND
!                               (Pressure, dd, ff) , 
!                               (Pressure, T, Td)
!
  CASE ('STANDARD')
    
    READ(ILUPRE,*) ZZGROUND                 ! Read data at ground level
    READ(ILUPRE,*) ZPGROUND
    READ(ILUPRE,*) ZTGROUND
    READ(ILUPRE,*) ZTDGROUND
!
    READ(ILUPRE,*) ILEVELU               ! Read number of wind levels
    ALLOCATE(ZPRESSU(ILEVELU))           ! Allocate memory for arrays to be read
    ALLOCATE(ZDD(ILEVELU),ZFF(ILEVELU))  
    ALLOCATE(ZHEIGHTU(ILEVELU))          ! Allocate memory for needed 
    ALLOCATE(ZU(ILEVELU),ZV(ILEVELU))    !  arrays
    ALLOCATE(ZTHVU(ILEVELU))             ! Allocate memory for intermediate
                                         !  arrays 
!
    DO JKU = 1,ILEVELU                   ! Read data at wind levels
      READ(ILUPRE,*) ZPRESSU(JKU),ZDD(JKU),ZFF(JKU)  
    END DO  
!
    READ(ILUPRE,*) ILEVELM           ! Read number of mass levels
                                     ! including the ground level
    ALLOCATE(ZPRESSM(ILEVELM))     ! Allocate memory for arrays to be read
    ALLOCATE(ZT(ILEVELM))                                   
    ALLOCATE(ZTD(ILEVELM))                                  
    ALLOCATE(ZHEIGHTM(ILEVELM))    ! Allocate memory for needed 
    ALLOCATE(ZTHV(ILEVELM))        !  arrays
    ALLOCATE(ZMR(ILEVELM))
    ALLOCATE(ZTV(ILEVELM))         ! Allocate memory for intermediate arrays
    ALLOCATE(ZTHL(ILEVELM))
    ALLOCATE(ZRT(ILEVELM))                                              
!
!   
    DO JKM= 2,ILEVELM                ! Read data at mass levels
      READ(ILUPRE,*) ZPRESSM(JKM),ZT(JKM),ZTD(JKM)
    END DO
    ZPRESSM(1)=ZPGROUND                 ! Mass level 1 is at the ground
    ZT(1)=ZTGROUND
    ZTD(1)=ZTDGROUND 
!   
!     recover the North-South and West-East wind components
    ZU(:) = ZFF(:)*COS(ZRADSDG*(270.-ZDD(:)) )              
    ZV(:) = ZFF(:)*SIN(ZRADSDG*(270.-ZDD(:)) )
!
!     compute vapor mixing ratio 
    ZMR(:)  = SM_FOES(ZTD(:))                            &      
            / ( (ZPRESSM(:) - SM_FOES(ZTD(:))) * ZRVSRD )       
! 
!     compute Tv
    ZTV(:) = ZT(:) * (1. + ZRVSRD * ZMR(:))/(1.+ZMR(:))         
!
!     compute thetav
    ZTHV(:) = ZTV(:) * (XP00/ ZPRESSM(:)) **(ZRDSCPD)         
!
!     compute height at the mass levels of the RS
    ZHEIGHTM(:) = HEIGHT_PRESS(ZPRESSM,ZTHV,ZPGROUND,ZTHV(1),ZZGROUND) 
!
!     compute thetav and height at the wind levels of the RS
    ZTHVU(:) = THETAVPU_THETAVPM(ZPRESSM,ZPRESSU,ZTHV)        
    ZHEIGHTU(:) = HEIGHT_PRESS(ZPRESSU,ZTHVU,ZPGROUND,ZTHV(1),ZZGROUND)
!
!   Compute Thetal and Rt
    ZRT(:)=ZMR(:)
    ZTHL(:)=ZTHV(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))
!
!*       2.2  PUVTHVMR case : zGROUND, PGROUND, ThvGROUND, RGROUND
!                               (Pressure, U, V) , 
!                               (Pressure, THv, R)
!
  CASE ('PUVTHVMR') 
!  
!     Read data at ground level      
    READ(ILUPRE,*) ZZGROUND                       
    READ(ILUPRE,*) ZPGROUND
    READ(ILUPRE,*) ZTHVGROUND
    READ(ILUPRE,*) ZMRGROUND
!
!     Read number of wind levels
    READ(ILUPRE,*) ILEVELU  
! 
!     Allocate the required memory               
    ALLOCATE(ZPRESSU(ILEVELU))              
    ALLOCATE(ZU(ILEVELU),ZV(ILEVELU))    
    ALLOCATE(ZHEIGHTU(ILEVELU))             
    ALLOCATE(ZTHVU(ILEVELU))                
!
!     Read the data at each wind level of the RS   
    DO JKU =1,ILEVELU
      READ(ILUPRE,*) ZPRESSU(JKU),ZU(JKU),ZV(JKU)   
    END DO
! 
!     Read number of mass levels
    READ(ILUPRE,*) ILEVELM
!
!     Allocate the required memory     
    ALLOCATE(ZPRESSM(ILEVELM))          
    ALLOCATE(ZHEIGHTM(ILEVELM))
    ALLOCATE(ZTHV(ILEVELM))
    ALLOCATE(ZMR(ILEVELM))
    ALLOCATE(ZTHL(ILEVELM))
    ALLOCATE(ZRT(ILEVELM))
!   
!     Read the data at each mass level of the RS   
    DO JKM = 2,ILEVELM
      READ(ILUPRE,*) ZPRESSM(JKM),ZTHV(JKM),ZMR(JKM)
    END DO
!
!     Complete the mass arrays with the ground informations read in EXPRE file
    ZPRESSM(1) = ZPGROUND
    ZTHV(1)    = ZTHVGROUND
    ZMR(1)     = ZMRGROUND 
!
!     Compute height of the mass levels of the RS
    ZHEIGHTM(:) = HEIGHT_PRESS(ZPRESSM,ZTHV,ZPGROUND,ZTHV(1),ZZGROUND) 
!                                                                 
!     Compute thetav and heigth at the wind levels of the RS 
    ZTHVU(:) = THETAVPU_THETAVPM(ZPRESSM,ZPRESSU,ZTHV)         
    ZHEIGHTU(:) = HEIGHT_PRESS(ZPRESSU,ZTHVU,ZPGROUND,ZTHV(1),ZZGROUND)
!
! on interpole thetal(=theta quand il n'y a pas d'eau liquide) et r total
    ZRT(:)=ZMR(:)
    ZTHL(:)=ZTHV(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))
!
!*       2.3  PUVTHVHU case :  zGROUND, PGROUND, ThvGROUND, HuGROUND
!                               (Pressure, U, V) , 
!                               (Pressure, THv, Hu)
!
  CASE ('PUVTHVHU')
! 
!     Read data at ground level      
    READ(ILUPRE,*) ZZGROUND                
    READ(ILUPRE,*) ZPGROUND
    READ(ILUPRE,*) ZTHVGROUND
    READ(ILUPRE,*) ZHUGROUND
!
!     Read number of wind levels
    READ(ILUPRE,*) ILEVELU
!
!     Allocate the required memory     
    ALLOCATE(ZPRESSU(ILEVELU))             
    ALLOCATE(ZU(ILEVELU),ZV(ILEVELU))
    ALLOCATE(ZHEIGHTU(ILEVELU))            
    ALLOCATE(ZTHVU(ILEVELU))               
!                                           
!     Read the data at each wind level of the RS   
    DO JKU =1,ILEVELU
      READ(ILUPRE,*) ZPRESSU(JKU),ZU(JKU),ZV(JKU)   
    END DO
! 
!     Read number of mass levels
    READ(ILUPRE,*) ILEVELM
!
!     Allocate the required memory  
    ALLOCATE(ZPRESSM(ILEVELM))           
    ALLOCATE(ZTHV(ILEVELM))
    ALLOCATE(ZHU(ILEVELM))
    ALLOCATE(ZHEIGHTM(ILEVELM))          
    ALLOCATE(ZMR(ILEVELM))
    ALLOCATE(ZTV(ILEVELM))     
    ALLOCATE(ZTHL(ILEVELM))
    ALLOCATE(ZRT(ILEVELM))
!   
!     Read the data at each mass level of the RS   
    DO JKM = 2,ILEVELM
      READ(ILUPRE,*) ZPRESSM(JKM),ZTHV(JKM),ZHU(JKM)
    END DO
!
!     Complete the mass arrays with the ground informations read in EXPRE file
    ZPRESSM(1) = ZPGROUND                    ! Mass level 1 is at the ground
    ZTHV(1)    = ZTHVGROUND
    ZHU(1)     = ZHUGROUND
! 
!     Compute Tv
    ZTV(:)=ZTHV(:) * (ZPRESSM(:) / XP00) ** ZRDSCPD  
!
!     Compte mixing ratio      
    ZMR(:)=SM_PMR_HU(ZPRESSM(:),ZTV(:),ZHU(:),SPREAD(ZMR(:),2,1))     
!                                                          
!     Compute height of the mass levels of the RS
    ZHEIGHTM(:) = HEIGHT_PRESS(ZPRESSM,ZTHV,ZPGROUND,ZTHV(1),ZZGROUND)
!                                                                 
!     Compute thetav and height of the wind levels of the RS                                                                
    ZTHVU(:) = THETAVPU_THETAVPM(ZPRESSM,ZPRESSU,ZTHV)   
    ZHEIGHTU(:) = HEIGHT_PRESS(ZPRESSU,ZTHVU,ZPGROUND,ZTHV(1),ZZGROUND)
! 
! on interpole thetal(=theta quand il n'y a pas d'eau liquide) et r total
    ZRT(:)=ZMR(:)
    ZTHL(:)=ZTHV(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))
!
!*       2.4  ZUVTHVHU case :  zGROUND, PGROUND, ThvGROUND, HuGROUND
!                               (height, U, V) , 
!                               (height, THv, Hu)
!
  CASE ('ZUVTHVHU')
!     Read data at ground level      
    READ(ILUPRE,*) ZZGROUND                
    READ(ILUPRE,*) ZPGROUND
    READ(ILUPRE,*) ZTHVGROUND
    READ(ILUPRE,*) ZHUGROUND
!
!     Read number of wind levels
    READ(ILUPRE,*) ILEVELU
!
!     Allocate the required memory 
    ALLOCATE(ZPRESSU(ILEVELU))
    ALLOCATE(ZHEIGHTU(ILEVELU))           
    ALLOCATE(ZU(ILEVELU),ZV(ILEVELU))
!
!
!     Read the data at each wind level of the RS   
    DO JKU = 1,ILEVELU
      READ(ILUPRE,*) ZHEIGHTU(JKU),ZU(JKU),ZV(JKU)  
    END DO
! 
!     Read number of mass levels
    READ(ILUPRE,*) ILEVELM
!
!     Allocate the required memory 
    ALLOCATE(ZHEIGHTM(ILEVELM))          
    ALLOCATE(ZTHV(ILEVELM))
    ALLOCATE(ZHU(ILEVELM))
    ALLOCATE(ZMR(ILEVELM))               
    ALLOCATE(ZPRESSM(ILEVELM))           
    ALLOCATE(ZTV(ILEVELM))
    ALLOCATE(ZTHL(ILEVELM))
    ALLOCATE(ZRT(ILEVELM))    
!   
!     Read the data at each mass level of the RS
    DO JKM = 2,ILEVELM
      READ(ILUPRE,*) ZHEIGHTM(JKM),ZTHV(JKM),ZHU(JKM)
    END DO 
!
!     Complete the mass arrays with the ground informations read in EXPRE file
    ZHEIGHTM(1) = ZZGROUND                   ! Mass level 1 is at the ground 
    ZTHV(1)    = ZTHVGROUND
    ZHU(1)     = ZHUGROUND
!
!     Compute Pressure at the mass levels of the RS
    ZPRESSM= PRESS_HEIGHT(ZHEIGHTM,ZTHV,ZPGROUND,ZTHV(1),ZHEIGHTM(1)) 
!                                                       
!     Compute Tv and the mixing ratio at the mass levels of the RS
    ZTV(:)=ZTHV(:) * (ZPRESSM(:) / XP00) ** ZRDSCPD   
    ZMR(:)=SM_PMR_HU(ZPRESSM(:),ZTV(:),ZHU(:),SPREAD(ZMR(:),2,1)) 
!
! on interpole thetal(=theta quand il n'y a pas d'eau liquide) et r total
    ZRT(:)=ZMR(:)
    ZTHL(:)=ZTHV(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))
!
!
!*       2.5  ZUVTHVMR case :  zGROUND, PGROUND, ThvGROUND, RGROUND
!                               (height, U, V) , 
!                               (height, THv, R)
!
!
  CASE ('ZUVTHVMR')
!     Read data at ground level      
    READ(ILUPRE,*) ZZGROUND
    READ(ILUPRE,*) ZPGROUND
    READ(ILUPRE,*) ZTHVGROUND
    READ(ILUPRE,*) ZMRGROUND
!
!     Read number of wind levels
    READ(ILUPRE,*) ILEVELU
!
!     Allocate the required memory 
    ALLOCATE(ZHEIGHTU(ILEVELU))         
    ALLOCATE(ZU(ILEVELU),ZV(ILEVELU))
!
!     Read the data at each wind level of the RS
    DO JKU = 1,ILEVELU
      READ(ILUPRE,*) ZHEIGHTU(JKU),ZU(JKU),ZV(JKU)  
    END DO
! 
!     Read number of mass levels
    READ(ILUPRE,*) ILEVELM
!
!     Allocate the required memory 
    ALLOCATE(ZHEIGHTM(ILEVELM))           
    ALLOCATE(ZTHV(ILEVELM))
    ALLOCATE(ZMR(ILEVELM))
    ALLOCATE(ZTHL(ILEVELM))
    ALLOCATE(ZRT(ILEVELM))    
!
!     Read the data at each mass level of the RS
    DO JKM=2,ILEVELM    
      READ(ILUPRE,*) ZHEIGHTM(JKM),ZTHV(JKM),ZMR(JKM) 
    END DO
!
!     Complete the mass arrays with the ground informations read in EXPRE file
    ZHEIGHTM(1)= ZZGROUND                      ! Mass level 1 is at the ground
    ZTHV(1)    = ZTHVGROUND
    ZMR(1)     = ZMRGROUND
! on interpole thetal(=theta quand il n'y a pas d'eau liquide) et r total
    ZRT(:)=ZMR(:)
    ZTHL(:)=ZTHV(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))    
!
!
!*       2.6  PUVTHDMR case : zGROUND, PGROUND, ThdGROUND, RGROUND
!                               (Pressure, U, V) , 
!                               (Pressure, THd, R)
!
  CASE ('PUVTHDMR')
!     Read data at ground level      
    READ(ILUPRE,*) ZZGROUND               
    READ(ILUPRE,*) ZPGROUND
    READ(ILUPRE,*) ZTHDGROUND
    READ(ILUPRE,*) ZMRGROUND
!
!     Read number of wind levels
    READ(ILUPRE,*) ILEVELU
!
!     Allocate the required memory 
    ALLOCATE(ZPRESSU(ILEVELU))          
    ALLOCATE(ZU(ILEVELU),ZV(ILEVELU))
    ALLOCATE(ZTHVU(ILEVELU))
    ALLOCATE(ZHEIGHTU(ILEVELU))         
!
!     Read the data at each wind level of the RS
    DO JKU =1,ILEVELU    
      READ(ILUPRE,*) ZPRESSU(JKU),ZU(JKU),ZV(JKU)
    END DO
!
!     Read number of mass levels 
    READ(ILUPRE,*) ILEVELM
!
!     Allocate the required memory 
    ALLOCATE(ZPRESSM(ILEVELM))         
    ALLOCATE(ZTHD(ILEVELM))
    ALLOCATE(ZMR(ILEVELM))
    ALLOCATE(ZMRC(ILEVELM))
    ZMRC=0
    ALLOCATE(ZMRI(ILEVELM))
    ZMRI=0
    ALLOCATE(ZHEIGHTM(ILEVELM))        
    ALLOCATE(ZTHV(ILEVELM))
    ALLOCATE(ZTHL(ILEVELM))
    ALLOCATE(ZRT(ILEVELM))     
!
!     Read the data at each mass level of the RS
    DO JKM=2,ILEVELM
      IF(LUSERI) THEN
        READ(ILUPRE,*) ZPRESSM(JKM),ZTHD(JKM),ZMR(JKM),ZMRC(JKM),ZMRI(JKM)
      ELSEIF (GUSERC) THEN
        READ(ILUPRE,*) ZPRESSM(JKM),ZTHD(JKM),ZMR(JKM),ZMRC(JKM)
      ELSE
        READ(ILUPRE,*) ZPRESSM(JKM),ZTHD(JKM),ZMR(JKM)
      ENDIF
    END DO
!
!     Complete the mass arrays with the ground informations read in EXPRE file
    ZPRESSM(1) = ZPGROUND                     ! Mass level 1 is at  the ground
    ZTHD(1)    = ZTHDGROUND
    ZMR(1)     = ZMRGROUND 
    IF(GUSERC) ZMRC(1) = ZMRC(2)
    IF(LUSERI) ZMRI(1) = ZMRI(2)
!
!     Compute thetav at the mass levels of the RS
      ZTHV(:) = ZTHD(:) * (1. + ZRVSRD *ZMR(:))/(1.+ZMR(:)+ZMRC(:)+ZMRI(:))
!
!     Compute the heights at the mass levels of the RS
    ZHEIGHTM(:) = HEIGHT_PRESS(ZPRESSM,ZTHV,ZPGROUND,ZTHV(1),ZZGROUND)
!
!     Compute thetav and heights of the wind levels
    ZTHVU(:) = THETAVPU_THETAVPM(ZPRESSM,ZPRESSU,ZTHV)   
    ZHEIGHTU(:) = HEIGHT_PRESS(ZPRESSU,ZTHVU,ZPGROUND,ZTHV(1),ZZGROUND)
!
!    Compute Theta l and Rt  
    IF (.NOT. GUSERC  .AND. .NOT. LUSERI) THEN
      ZRT(:)=ZMR(:)
      ZTHL(:)=ZTHV(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))   
    ELSE
       ALLOCATE(ZEXN(ILEVELM))         
       ALLOCATE(ZT(ILEVELM))       
       ALLOCATE(ZCPH(ILEVELM))        
       ALLOCATE(ZLVOCPEXN(ILEVELM))        
       ALLOCATE(ZLSOCPEXN(ILEVELM))               
       ZRT(:)=ZMR(:)+ZMRI(:)+ZMRC(:)
       ZEXN(:)=(ZPRESSM/XP00) ** (XRD/XCPD)
       ZT(:)=ZTHV*(ZPRESSM(:)/XP00)**(ZRDSCPD)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))
       ZCPH(:)=XCPD+ XCPV * ZMR(:)+ XCL *ZMRC(:)  + XCI * ZMRI(:)
       ZLVOCPEXN(:) = (XLVTT + (XCPV-XCL) * (ZT(:)-XTT))/(ZCPH*ZEXN(:))
       ZLSOCPEXN(:) = (XLSTT + (XCPV-XCI) * (ZT(:)-XTT))/(ZCPH*ZEXN(:))
       ZTHL(:)=ZTHV(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))-ZLVOCPEXN(:)*ZMRC(:)-ZLSOCPEXN(:)*ZMRI(:)
       DEALLOCATE(ZEXN)         
       DEALLOCATE(ZT)       
       DEALLOCATE(ZCPH)        
       DEALLOCATE(ZLVOCPEXN)        
       DEALLOCATE(ZLSOCPEXN)
    ENDIF   
!
!
!*       2.7  PUVTHDHU case :  zGROUND, PGROUND, ThdGROUND, HuGROUND
!                               (Pressure, U, V) , 
!                               (Pressure, THd, Hu)
!
  CASE ('PUVTHDHU')
!     Read data at ground level      
    READ(ILUPRE,*) ZZGROUND              
    READ(ILUPRE,*) ZPGROUND
    READ(ILUPRE,*) ZTHDGROUND
    READ(ILUPRE,*) ZHUGROUND
!
!     Read number of wind levels  
    READ(ILUPRE,*) ILEVELU
!
!     Allocate the required memory
    ALLOCATE(ZPRESSU(ILEVELU))         
    ALLOCATE(ZU(ILEVELU),ZV(ILEVELU))
    ALLOCATE(ZTHVU(ILEVELU))
    ALLOCATE(ZHEIGHTU(ILEVELU))        
!
!     Read the data at each wind level of the RS
    DO JKU = 1,ILEVELU
      READ(ILUPRE,*) ZPRESSU(JKU),ZU(JKU),ZV(JKU)
    END DO
!
!     Read number of mass levels 
    READ(ILUPRE,*) ILEVELM
!
! Allocate the required memory
    ALLOCATE(ZPRESSM(ILEVELM))          
    ALLOCATE(ZTHD(ILEVELM))
    ALLOCATE(ZHU(ILEVELM))
    ALLOCATE(ZHEIGHTM(ILEVELM))          
    ALLOCATE(ZTHV(ILEVELM))
    ALLOCATE(ZMR(ILEVELM))
    ALLOCATE(ZT(ILEVELM)) 
    ALLOCATE(ZTHL(ILEVELM))
    ALLOCATE(ZRT(ILEVELM))                                          
!   
!     Read the data at each mass level of the RS
    DO JKM =2,ILEVELM
      READ(ILUPRE,*) ZPRESSM(JKM),ZTHD(JKM), ZHU(JKM)
    END DO
!     Complete the mass arrays with the ground informations read in EXPRE file
    ZPRESSM(1) = ZPGROUND                     ! Mass level 1 is at the  ground
    ZTHD(1)    = ZTHDGROUND
    ZHU(1)     = ZHUGROUND
! 
    ZT(:) = ZTHD(:) * (ZPRESSM(:)/XP00)**ZRDSCPD  ! compute T and mixing ratio
    ZMR(:) = ZRDSRV*SM_FOES(ZT(:))/((ZPRESSM(:)*100./ZHU(:)) -SM_FOES(ZT(:)))
                         
!     Compute thetav at the mass levels of the RS
    ZTHV(:) = ZTHD(:)  * (1. + ZRVSRD *ZMR(:))/(1.+ZMR(:)) 
!
!     Compute height at mass levels
    ZHEIGHTM(:) = HEIGHT_PRESS(ZPRESSM,ZTHV,ZPGROUND,ZTHV(1),ZZGROUND) 
!
!     Compute thetav and heights of the wind levels 
    ZTHVU(:) = THETAVPU_THETAVPM(ZPRESSM,ZPRESSU,ZTHV) 
    ZHEIGHTU(:) = HEIGHT_PRESS(ZPRESSU,ZTHVU,ZPGROUND,ZTHV(1),ZZGROUND)   
!
!   Compute thetal and Rt
    ZRT(:)=ZMR(:)
    ZTHL(:)=ZTHV(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))      
!
!*       2.8  ZUVTHDMR case :  zGROUND, PGROUND, ThdGROUND, RGROUND
!                               (height, U, V) , 
!                               (height, THd, R)
!
  CASE ('ZUVTHDMR')
! Read data at ground level
    READ(ILUPRE,*) ZZGROUND                  
    READ(ILUPRE,*) ZPGROUND
    READ(ILUPRE,*) ZTHDGROUND
    READ(ILUPRE,*) ZMRGROUND
!
!     Read number of wind levels  
    READ(ILUPRE,*) ILEVELU
!
!   Allocate required memory 
    ALLOCATE(ZHEIGHTU(ILEVELU))           
    ALLOCATE(ZU(ILEVELU),ZV(ILEVELU))
!
!     Read the data at each wind level of the RS
    DO JKU = 1,ILEVELU
      READ(ILUPRE,*) ZHEIGHTU(JKU),ZU(JKU),ZV(JKU)
    END DO
! 
!     Read number of mass levels  
     READ(ILUPRE,*) ILEVELM
!
!   Allocate required memory 
    ALLOCATE(ZHEIGHTM(ILEVELM))           
    ALLOCATE(ZTHD(ILEVELM))
    ALLOCATE(ZMR(ILEVELM))
    ALLOCATE(ZTHV(ILEVELM)) 
    ALLOCATE(ZMRC(ILEVELM))
    ZMRC=0
    ALLOCATE(ZMRI(ILEVELM))
    ZMRI=0
    ALLOCATE(ZTHL(ILEVELM))
    ALLOCATE(ZRT(ILEVELM))
!
!     Read the data at each mass level of the RS
    DO JKM= 2,ILEVELM
      IF(LUSERI) THEN
        READ(ILUPRE,*) ZHEIGHTM(JKM),ZTHD(JKM),ZMR(JKM),ZMRC(JKM),ZMRI(JKM)
      ELSEIF (GUSERC) THEN
        READ(ILUPRE,*) ZHEIGHTM(JKM),ZTHD(JKM),ZMR(JKM),ZMRC(JKM)
      ELSE
        READ(ILUPRE,*) ZHEIGHTM(JKM),ZTHD(JKM),ZMR(JKM)
      ENDIF
    END DO
!     Complete the mass arrays with the ground informations read in EXPRE file
    ZHEIGHTM(1) = ZZGROUND                     ! Mass level 1 is at ground
    ZTHD(1)     = ZTHDGROUND
    ZMR(1)      = ZMRGROUND
    IF(GUSERC) ZMRC(1) = ZMRC(2)
    IF(LUSERI) ZMRI(1) = ZMRI(2)
!     Compute thetav at the mass levels of the RS
    IF(LUSERI) THEN
      ZTHV(:) = ZTHD(:) * (1. + ZRVSRD *ZMR(:))/(1.+ZMR(:)+ZMRC(:)+ZMRI(:))
    ELSEIF (GUSERC) THEN
      ZTHV(:) = ZTHD(:) * (1. + ZRVSRD *ZMR(:))/(1.+ZMR(:)+ZMRC(:))
    ELSE
      ZTHV(:) = ZTHD(:) * (1. + ZRVSRD *ZMR(:))/(1.+ZMR(:))
    ENDIF    
!    
!    Compute Theta l and Rt  
    IF (.NOT. GUSERC  .AND. .NOT. LUSERI) THEN
      ZRT(:)=ZMR(:)
      ZTHL(:)=ZTHV(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))   
    ELSE
       ALLOCATE(ZEXN(ILEVELM))  
       ALLOCATE(ZEXNFLUX(ILEVELM))        
       ALLOCATE(ZT(ILEVELM))       
       ALLOCATE(ZCPH(ILEVELM))        
       ALLOCATE(ZLVOCPEXN(ILEVELM))        
       ALLOCATE(ZLSOCPEXN(ILEVELM))  
       ZRT(:)=ZMR(:)+ZMRI(:)+ZMRC(:)
       ZEXNSURF=(ZPGROUND/XP00) ** (XRD/XCPD)
       CALL COMPUTE_EXNER_FROM_GROUND(ZTHV,ZHEIGHTM,ZEXNSURF,ZEXNFLUX,ZEXN)
       ZT(:)=ZTHV*ZEXN(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))
       ZCPH(:)=XCPD+ XCPV * ZMR(:)+ XCL *ZMRC(:)  + XCI * ZMRI(:)
       ZLVOCPEXN(:) = (XLVTT + (XCPV-XCL) * (ZT(:)-XTT))/(ZCPH*ZEXN(:))
       ZLSOCPEXN(:) = (XLSTT + (XCPV-XCI) * (ZT(:)-XTT))/(ZCPH*ZEXN(:))
       ZTHL(:)=ZTHV(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))-ZLVOCPEXN(:)*ZMRC(:)-ZLSOCPEXN(:)*ZMRI(:)
       DEALLOCATE(ZEXN)         
       DEALLOCATE(ZEXNFLUX)         
       DEALLOCATE(ZT)       
       DEALLOCATE(ZCPH)        
       DEALLOCATE(ZLVOCPEXN)        
       DEALLOCATE(ZLSOCPEXN)
    ENDIF  
!
! 2.9 ZUVTHLMR case : zGROUND, PGROUND, ThdGROUND, RGROUND
!                               (height, U, V)
!                               (height, THL, Rt)

!
  CASE ('ZUVTHLMR')
! Read data at ground level
    READ(ILUPRE,*) ZZGROUND                  
    READ(ILUPRE,*) ZPGROUND
    READ(ILUPRE,*) ZTHLGROUND
    READ(ILUPRE,*) ZMRGROUND
!
!   Read number of wind levels  
    READ(ILUPRE,*) ILEVELU
!
!   Allocate required memory 
    ALLOCATE(ZHEIGHTU(ILEVELU))           
    ALLOCATE(ZU(ILEVELU),ZV(ILEVELU))
!
!   Read the data at each wind level of the RS
    DO JKU = 1,ILEVELU
      READ(ILUPRE,*) ZHEIGHTU(JKU),ZU(JKU),ZV(JKU)
    END DO
! 
!   Read number of mass levels  
    READ(ILUPRE,*) ILEVELM
!
!   Allocate required memory 
    ALLOCATE(ZHEIGHTM(ILEVELM))           
    ALLOCATE(ZTHL(ILEVELM))
    ALLOCATE(ZTH(ILEVELM))
    ALLOCATE(ZMR(ILEVELM))
    ALLOCATE(ZTHV(ILEVELM)) 
    ALLOCATE(ZMRC(ILEVELM))
    ZMRC=0
    ALLOCATE(ZMRI(ILEVELM))
    ZMRI=0    
    ALLOCATE(ZRT(ILEVELM))
!
!   Read the data at each mass level of the RS
    DO JKM= 2,ILEVELM
!     IF(LUSERI) THEN
!       READ(ILUPRE,*) ZHEIGHTM(JKM),ZTHL(JKM),ZMR(JKM),ZMRC(JKM),ZMRI(JKM)
!     ELSEIF (GUSERC) THEN
      IF (GUSERC) THEN
        READ(ILUPRE,*) ZHEIGHTM(JKM),ZTHL(JKM),ZMR(JKM),ZMRC(JKM)
      ELSE
        READ(ILUPRE,*) ZHEIGHTM(JKM),ZTHL(JKM),ZMR(JKM)
      ENDIF
    END DO
!   Complete the mass arrays with the ground informations read in EXPRE file
    ZHEIGHTM(1) = ZZGROUND                     ! Mass level 1 is at ground
    ZTHL(1)     = ZTHLGROUND
    ZMR(1)      = ZMRGROUND
    IF(GUSERC) ZMRC(1) = ZMRC(2)
!   IF(LUSERI) ZMRI(1) = ZMRI(2)
!
!   Compute Rt
    ZRT(:)=ZMR+ZMRC+ZMRI
!
!*       2.10  PUVTHU case :  zGROUND, PGROUND, TempGROUND, HuGROUND
!                               (Pressure, U, V) , 
!                               (Pressure, Temp, Hu)
!
  CASE ('PUVTHU')
!     Read data at ground level      
    READ(ILUPRE,*) ZZGROUND              
    READ(ILUPRE,*) ZPGROUND
    READ(ILUPRE,*) ZTGROUND
    READ(ILUPRE,*) ZHUGROUND
!
!     Read number of wind levels  
    READ(ILUPRE,*) ILEVELU
!
!     Allocate the required memory
    ALLOCATE(ZPRESSU(ILEVELU))         
    ALLOCATE(ZU(ILEVELU),ZV(ILEVELU))
    ALLOCATE(ZTHVU(ILEVELU))
    ALLOCATE(ZHEIGHTU(ILEVELU))        
!
!     Read the data at each wind level of the RS
    DO JKU = 1,ILEVELU
      READ(ILUPRE,*) ZPRESSU(JKU),ZU(JKU),ZV(JKU)
    END DO
!
!     Read number of mass levels 
    READ(ILUPRE,*) ILEVELM
!
! Allocate the required memory
    ALLOCATE(ZPRESSM(ILEVELM))          
    ALLOCATE(ZTHD(ILEVELM))
    ALLOCATE(ZHU(ILEVELM))
    ALLOCATE(ZHEIGHTM(ILEVELM))          
    ALLOCATE(ZTHV(ILEVELM))
    ALLOCATE(ZMR(ILEVELM))
    ALLOCATE(ZT(ILEVELM))  
    ALLOCATE(ZTHL(ILEVELM))
    ALLOCATE(ZRT(ILEVELM))

!   
!     Read the data at each mass level of the RS
    DO JKM =2,ILEVELM
      READ(ILUPRE,*) ZPRESSM(JKM),ZT(JKM), ZHU(JKM)
    END DO
!     Complete the mass arrays with the ground informations read in EXPRE file
    ZPRESSM(1) = ZPGROUND                     ! Mass level 1 is at the  ground
    ZT(1)    = ZTGROUND
    ZHU(1)     = ZHUGROUND
! 
    ZTHD(:) = ZT(:) / (ZPRESSM(:)/XP00)**ZRDSCPD  ! compute THD and mixing ratio
    ZMR(:) = ZRDSRV*SM_FOES(ZT(:))/((ZPRESSM(:)*100./ZHU(:)) -SM_FOES(ZT(:)))
!     Compute thetav at the mass levels of the RS
    ZTHV(:) = ZTHD(:)  * (1. + ZRVSRD *ZMR(:))/(1.+ZMR(:)) 
!
!     Compute height at mass levels
    ZHEIGHTM(:) = HEIGHT_PRESS(ZPRESSM,ZTHV,ZPGROUND,ZTHV(1),ZZGROUND) 
!
!     Compute thetav and heights of the wind levels 
    ZTHVU(:) = THETAVPU_THETAVPM(ZPRESSM,ZPRESSU,ZTHV) 
    ZHEIGHTU(:) = HEIGHT_PRESS(ZPRESSU,ZTHVU,ZPGROUND,ZTHV(1),ZZGROUND)   
!
! on interpole thetal(=theta quand il n'y a pas d'eau liquide) et r total
     ZRT(:)=ZMR(:)
     ZTHL(:)=ZTHV(:)*(1+ZRT(:))/(1+ZRVSRD*ZRT(:))
  CASE DEFAULT
    CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_RSOU','data type YKIND='//TRIM(YKIND)//' in PREFILE unknown')
END SELECT
!
!-------------------------------------------------------------------------------
!
!*	 3.     INTERPOLATE ON THE  VERTICAL MIXED MODEL GRID 
!	        ---------------------------------------------------------
!
!
! 
IKU=SIZE(XZHAT)
!
!*       3.1    Compute mixed grid
! 
IF (PRESENT(PCORIOZ)) THEN 
! LGEOSBAL=T (no shift allowed, MNH grid without ororgraphy)
  ZZS_LS(:,:)=0
ELSE
  IF (OSHIFT) THEN
    ZZS_LS(:,:)=ZZGROUND
  ELSE
    ZZS_LS(:,:)=0
  ENDIF
ENDIF
CALL VERT_COORD(LSLEVE,ZZS_LS,ZZS_LS,XLEN1,XLEN2,XZHAT,ZZFLUX_MX)
ZZMASS_MX(:,:,:)=MZF(ZZFLUX_MX)
ZZMASS_MX(:,:,IKU)=1.5*ZZFLUX_MX(:,:,IKU)-0.5*ZZFLUX_MX(:,:,IKU-1)
!
!*        3.2  Interpolate and extrapolate U and V on w- mixed grid levels
!
!* vertical grid at initialization profile location
GPROFILE_IN_PROC=(KILOC+JPHEXT-IXOR_ll+1>=IIB .AND. KILOC+JPHEXT-IXOR_ll+1<=IIE) &
         & .AND. (KJLOC+JPHEXT-IYOR_ll+1>=IJB .AND. KJLOC+JPHEXT-IYOR_ll+1<=IJE)
!
IF (GPROFILE_IN_PROC) THEN
  ZZMASS_PROFILE(:) = ZZMASS_MX(KILOC+JPHEXT-IXOR_ll+1,KJLOC+JPHEXT-IYOR_ll+1,:)
  ZZFLUX_PROFILE(:) = ZZFLUX_MX(KILOC+JPHEXT-IXOR_ll+1,KJLOC+JPHEXT-IYOR_ll+1,:)
ELSE
  ZZMASS_PROFILE(:) = 0.
  ZZFLUX_PROFILE(:) = 0.
END IF
DO JK = 1,IKU
  CALL REDUCESUM_ll(ZZMASS_PROFILE(JK), IINFO_ll)
  CALL REDUCESUM_ll(ZZFLUX_PROFILE(JK), IINFO_ll)
END DO

! interpolation of U and V
DO JK = 1,IKU
  IF (ZZFLUX_PROFILE(JK) <= ZHEIGHTU(1)) THEN      ! extrapolation below the first level
          ZDZSDH  = (ZZFLUX_PROFILE(JK) - ZHEIGHTU(1)) / (ZHEIGHTU(2) - ZHEIGHTU(1))
          ZUW(JK) = ZU(1) + (ZU(2) - ZU(1)) * ZDZSDH
          ZVW(JK) = ZV(1) + (ZV(2) - ZV(1)) * ZDZSDH 
  ELSE IF (ZZFLUX_PROFILE(JK) > ZHEIGHTU(ILEVELU) ) THEN  ! extrapolation above the last
    ZDZSDH  = (ZZFLUX_PROFILE(JK) - ZHEIGHTU(ILEVELU))                   &   ! level
            / (ZHEIGHTU(ILEVELU) - ZHEIGHTU(ILEVELU-1))
    ZUW(JK) = ZU(ILEVELU) + (ZU(ILEVELU) -ZU(ILEVELU -1)) * ZDZSDH
    ZVW(JK) = ZV(ILEVELU) + (ZV(ILEVELU) -ZV(ILEVELU -1)) * ZDZSDH
  ELSE                       ! interpolation between the first and last levels
    DO JKLEV = 1,ILEVELU-1
      IF ( (ZZFLUX_PROFILE(JK) > ZHEIGHTU(JKLEV)).AND.                  &
           (ZZFLUX_PROFILE(JK) <= ZHEIGHTU(JKLEV+1)) )THEN 
        ZDZ1SDH = (ZZFLUX_PROFILE(JK) - ZHEIGHTU(JKLEV))                &
              / (ZHEIGHTU(JKLEV+1)-ZHEIGHTU(JKLEV)) 
        ZDZ2SDH = (ZHEIGHTU(JKLEV+1) - ZZFLUX_PROFILE(JK) )             &
              / (ZHEIGHTU(JKLEV+1)-ZHEIGHTU(JKLEV)) 
        ZUW(JK) = (ZU(JKLEV) * ZDZ2SDH) + (ZU(JKLEV+1) *ZDZ1SDH)
        ZVW(JK) = (ZV(JKLEV) * ZDZ2SDH) + (ZV(JKLEV+1) *ZDZ1SDH)
      END IF        
    END DO
  END IF
END DO
!
!*       3.3    Interpolate and extrapolate Thetav and r on  mass mixed grid levels
!
ZMRCM=0
ZMRIM=0
DO JK = 1,IKU
  IF (ZZMASS_PROFILE(JK) <= ZHEIGHTM(1)) THEN    ! extrapolation below the first level
    ZDZSDH  = (ZZMASS_PROFILE(JK) - ZHEIGHTM(1)) / (ZHEIGHTM(2) - ZHEIGHTM(1))
    ZTHLM(JK) = ZTHL(1) + (ZTHL(2) - ZTHL(1)) * ZDZSDH
    ZMRM(JK)  = ZRT(1) + (ZRT(2) - ZRT(1)) * ZDZSDH 
    IF (GUSERC) ZMRCM(JK)  = ZMRC(1) + (ZMRC(2) - ZMRC(1)) * ZDZSDH 
    IF (LUSERI) ZMRIM(JK)  = ZMRI(1) + (ZMRI(2) - ZMRI(1)) * ZDZSDH 
  ELSE IF (ZZMASS_PROFILE(JK) > ZHEIGHTM(ILEVELM) ) THEN  ! extrapolation above the last
    ZDZSDH  = (ZZMASS_PROFILE(JK) - ZHEIGHTM(ILEVELM))                   &  ! level
            / (ZHEIGHTM(ILEVELM) - ZHEIGHTM(ILEVELM-1))
    ZTHLM(JK) = ZTHL(ILEVELM) + (ZTHL(ILEVELM) -ZTHL(ILEVELM -1)) * ZDZSDH
    ZMRM(JK) = ZRT(ILEVELM) + (ZRT(ILEVELM) -ZRT(ILEVELM -1)) * ZDZSDH
    IF (GUSERC) ZMRCM(JK) = ZMRC(ILEVELM) + (ZMRC(ILEVELM) -ZMRC(ILEVELM -1)) * ZDZSDH
    IF (LUSERI) ZMRIM(JK) = ZMRI(ILEVELM) + (ZMRI(ILEVELM) -ZMRI(ILEVELM -1)) * ZDZSDH
  ELSE  ! interpolation between the first and last levels
    DO JKLEV = 1,ILEVELM-1
      IF ( (ZZMASS_PROFILE(JK) > ZHEIGHTM(JKLEV)).AND.                  &
           (ZZMASS_PROFILE(JK) <= ZHEIGHTM(JKLEV+1)) )THEN 
        ZDZ1SDH = (ZZMASS_PROFILE(JK) - ZHEIGHTM(JKLEV))                 &
              / (ZHEIGHTM(JKLEV+1)-ZHEIGHTM(JKLEV)) 
        ZDZ2SDH = (ZHEIGHTM(JKLEV+1) - ZZMASS_PROFILE(JK) )              &
              / (ZHEIGHTM(JKLEV+1)-ZHEIGHTM(JKLEV)) 
        ZTHLM(JK) = (ZTHL(JKLEV) * ZDZ2SDH) + (ZTHL(JKLEV+1) *ZDZ1SDH)
        ZMRM(JK)  = (ZRT(JKLEV) * ZDZ2SDH) + (ZRT(JKLEV+1) *ZDZ1SDH)
        IF (GUSERC) ZMRCM(JK)  = (ZMRC(JKLEV) * ZDZ2SDH) + (ZMRC(JKLEV+1) *ZDZ1SDH)
        IF (LUSERI) ZMRIM(JK)  = (ZMRI(JKLEV) * ZDZ2SDH) + (ZMRI(JKLEV+1) *ZDZ1SDH)
      END IF 
    END DO
  END IF
END DO
!
! Compute thetaV rv ri and Rc with adjustement
ALLOCATE(ZEXNFLUX(IKU))
ALLOCATE(ZEXNMASS(IKU))
ALLOCATE(ZPRESS(IKU))
ALLOCATE(ZPREFLUX(IKU))
ALLOCATE(ZFRAC_ICE(IKU))
ALLOCATE(ZRSATW(IKU))
ALLOCATE(ZRSATI(IKU))
ALLOCATE(ZMRT(IKU))
ZMRT=ZMRM+ZMRCM+ZMRIM
ZTHVM=ZTHLM
!
IF (LOCEAN) THEN
  ZRHODM(:)=XRH00OCEAN*(1.-XALPHAOC*(ZTHLM(:) - XTH00OCEAN)&
          +XBETAOC* (ZMRM(:)  - XSA00OCEAN))
  ZPREFLUX(IKU)=ZPTOP
  DO JK=IKU-1,2,-1
    ZPREFLUX(JK) = ZPREFLUX(JK+1) + XG*ZRHODM(JK)*(ZZFLUX_PROFILE(JK+1)-ZZFLUX_PROFILE(JK))
  END DO
  ZPGROUND=ZPREFLUX(2)
  WRITE(ILUOUT,FMT=*)'ZPGROUND i.e. Pressure at ocean domain bottom',ZPGROUND
  ZTHM=ZTHVM
ELSE
! Atmospheric case   
  ZEXNSURF=(ZPGROUND/XP00)**(XRD/XCPD)
  DO JLOOP=1,20 ! loop for pression 
    CALL COMPUTE_EXNER_FROM_GROUND(ZTHVM,ZZMASS_PROFILE(:),ZEXNSURF,ZEXNFLUX,ZEXNMASS)
    ZPRESS(:)=XP00*(ZEXNMASS(:))**(XCPD/XRD)
    CALL TH_R_FROM_THL_RT_1D('T',ZFRAC_ICE,ZPRESS,ZTHLM,ZMRT,ZTHM,ZMRM,ZMRCM,ZMRIM, &
                              ZRSATW, ZRSATI)
     ZTHVM(:)=ZTHM(:)*(1.+XRV/XRD*ZMRM(:))/(1.+(ZMRM(:)+ZMRIM(:)+ZMRCM(:)))
  ENDDO
ENDIF
!
DEALLOCATE(ZEXNFLUX)
DEALLOCATE(ZEXNMASS)
DEALLOCATE(ZPRESS)
DEALLOCATE(ZFRAC_ICE)
DEALLOCATE(ZRSATW)
DEALLOCATE(ZRSATI)       
DEALLOCATE(ZMRT)
!-------------------------------------------------------------------------------
!
!* 4.     COMPUTE FIELDS ON THE MODEL GRID (WITH OROGRAPHY)
!         -------------------------------------------------
CALL SET_MASS(TPFILE,GPROFILE_IN_PROC, ZZFLUX_PROFILE,                      &
              KILOC+JPHEXT,KJLOC+JPHEXT,ZZS_LS,ZZMASS_MX,ZZFLUX_MX,ZPGROUND,&
              ZTHVM,ZMRM,ZUW,ZVW,OSHIFT,OBOUSS,PJ,HFUNU,HFUNV,              &
              PMRCM=ZMRCM,PMRIM=ZMRIM,PCORIOZ=PCORIOZ)
!
DEALLOCATE(ZPREFLUX)
DEALLOCATE(ZHEIGHTM)          
DEALLOCATE(ZTHV)
DEALLOCATE(ZMR)
DEALLOCATE(ZTHL)
!-------------------------------------------------------------------------------
CONTAINS
  SUBROUTINE CHECK( ISTATUS, YLOC )
    INTEGER(KIND=CDFINT), INTENT(IN) :: ISTATUS
    CHARACTER(LEN=*),     INTENT(IN) :: YLOC

    IF( ISTATUS /= NF90_NOERR ) THEN
      CALL PRINT_MSG( NVERB_ERROR, 'IO', 'SET_RSOU', 'error at ' // Trim( yloc) // ': ' // NF90_STRERROR( ISTATUS ) )
    END IF
  END SUBROUTINE check
!
END SUBROUTINE SET_RSOU
