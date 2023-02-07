!MNH_LIC Copyright 1995-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    ########################
     MODULE MODI_RADIATIONS   
!    ########################
!
CONTAINS
!
!   ############################################################################
    SUBROUTINE RADIATIONS (TPFILE,OCLEAR_SKY,OCLOUD_ONLY,                      &
               KCLEARCOL_TM1,HEFRADL,HEFRADI,HOPWSW,HOPISW,HOPWLW,HOPILW,      &
               PFUDG, KDLON, KFLEV, KRAD_DIAG, KFLUX, KRAD, KAER, KSWB_OLD,    &
               KSWB_MNH,KLWB_MNH, KSTATM,KRAD_COLNBR,PCOSZEN,PSEA, PCORSOL,    &
               PDIR_ALB, PSCA_ALB,PEMIS, PCLDFR, PCCO2, PTSRAD, PSTATM,        &
               PTHT, PRT, PPABST, POZON, PAER, PDST_WL, PAER_CLIM, PSVT,       &
               PDTHRAD, PSRFLWD, PSRFSWD_DIR,PSRFSWD_DIF, PRHODREF, PZZ,       &
               PRADEFF, PSWU, PSWD, PLWU,PLWD, PDTHRADSW, PDTHRADLW            )
!   ############################################################################
!
!!****  *RADIATIONS * - routine to call the SW and LW radiation calculations
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to prepare the temperature, water vapor
!!    liquid water, cloud fraction, ozone profiles for the ECMWF radiation
!!    calculations. There is a great number of available radiative fluxes in
!!    the output, but only the potential temperature radiative tendency and the
!!    SW and LW surface fluxes are provided in the output of the routine.
!!    Two simplified computations are available (switches OCLEAR_SKY and
!!    OCLOUD_ONLY). When OCLOUD_ONLY is .TRUE. the computations are performed
!!    for the cloudy columns only. Furthermore with OCLEAR_SKY being .TRUE.
!!    the clear sky columns are averaged and the computations are made for
!!    the cloudy columns plus a single ensemble-mean clear sky column.
!!
!!**  METHOD
!!    ------
!!      First the temperature, water vapor, liquid water, cloud fraction
!!    and  profile arrays are built using the current model fields and
!!    the standard atmosphere for the upper layer filling.
!!    The standard atmosphere is used between the levels IKUP and
!!    KFLEV where KFLEV is the number of vertical levels for the radiation 
!!    computations.    
!!    The aerosols optical thickness and the ozone fields come directly
!!    from ini_radiation step (climatlogies used) and are already defined for KFLEV. 
!!    Surface parameter ( albedo, emiss ) are also defined from current surface fields.
!!    In the case of clear-sky or cloud-only approximations, the cloudy
!!    columns are selected by testing the vertically integrated cloud fraction
!!    and the radiation computations are performed for these columns plus the
!!    mean clear-sky one. In addition, columns where cloud have disapeared are determined
!!    by saving cloud trace between radiation step and they are also recalculated
!!    in cloud only step. In all case, the sun position correponds to  the centered
!!    time between 2 full radiation steps (determined in physparam).
!!      Then the ECMWF radiation package is called and the radiative
!!    heating/cooling tendancies are reformatted in case of partial
!!    computations.  In case of "cloud-only approximation" the only cloudy
!!    column radiative fields are updated.
!!
!!    EXTERNAL
!!    --------
!!      Subroutine ECMWF_RADIATION_VERS2 : ECMWF interface calling radiation routines
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST  : constants
!!        XP00 : reference pressure
!!        XCPD : calorific capacity of dry air at constant pressure
!!        XRD  : gas constant for dry air
!!      Module MODD_PARAMETERS : parameters
!!        JPHEXT : Extra columns on the horizontal boundaries
!!        JPVEXT : Extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine RADIATIONS )
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/02/95 
!!      J.Stein     20/12/95 add the array splitting in order to save memory
!!      J.-P. Pinty 19/11/96 change the split arrays, specific humidity
!!                           and add the ice phase
!!      J.Stein     22/06/97 use of the absolute pressure
!!      P.Jabouille 31/07/97 impose a zero humidity for dry simulation
!!      V.Masson    22/09/97 case of clear-sky approx. with no clear-sky column
!!      V.Masson    07/11/97 half level pressure defined from averaged Exner
!!                           function
!!      V.Masson    07/11/97 modification of junction between standard atm
!!                           and model for half level variables (top model
!!                           pressure and temperatures are used preferentially
!!                           to atm standard profile for the first point).
!!      P.Jabouille 24/08/98 impose positivity for ZQLAVE
!!      J.-P. Pinty 29/01/98 add storage for diagnostics
!!      J. Stein    18/07/99 add the ORAD_DIAG switch and keep inside the
!!                           subroutine the partial tendencies 
!!
!!      F.Solmon    04/03/01  MAJOR MODIFICATIONS, updated version of ECMWF radiation scheme
!!      P.Jabouille 05/05/03 bug in humidity conversion
!!      Y.Seity     25/08/03  KSWB=6 for SW direct and scattered surface 
!!                            downward fluxes used in surface scheme. 
!!      P. Tulet    01/20/05  climatologic SSA
!!      A. Grini    05/20/05  dust direct effect (optical properties)
!!      V.Masson, C.Lac 08/10 Correction of inversion of Diffuse and direct albedo
!!      B.Aouizerats 2010     Explicit aerosol optical properties
!!      C.Lac       11/2015   Correction on aerosols
!!      B.Vie            /13  LIMA
!!      J.Escobar 30/03/2017  : Management of compilation of ECMWF_RAD in REAL*8 with MNH_REAL=R4
!!      J.Escobar 29/06/2017  : Check if Pressure Decreasing with height <-> elsif PB & STOP 
!!      Q.LIBOIS  06/2017     : correction on CLOUD_ONLY
!!      Q.Libois  02/2018     : ECRAD
!!      Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      J.Escobar 28/06/2018 : Reproductible parallelisation of CLOUD_ONLY case
!!      J.Escobar 20/07/2018 : for real*4 compilation, convert with REAL(X) argument to SUM_DD... 
!!      P.Wautelet 22/01/2019: use standard FLUSH statement instead of non standard intrinsics
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 06/09/2022: small fix: GSURF_CLOUD was not set outside of physical domain
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE PARKIND1,         ONLY: JPRB
USE OYOESW    , ONLY : RTAUA    ,RPIZA    ,RCGA
!
USE MODD_CH_AEROSOL,  ONLY: LORILAM
USE MODD_CONF,        ONLY: LCARTESIAN
USE MODD_CST
USE MODD_DUST,        ONLY: LDUST
use modd_field,          only: tfieldmetadata, TYPEREAL
USE MODD_GRID ,       ONLY: XLAT0, XLON0
USE MODD_GRID_n ,     ONLY: XLAT, XLON
USE MODD_IO,          ONLY: TFILEDATA
USE MODD_LUNIT_n,     ONLY: TLUOUT
USE MODD_NSV,         ONLY: NSV_C2R2,NSV_C2R2BEG,NSV_C2R2END,     &
                            NSV_C1R3,NSV_C1R3BEG,NSV_C1R3END,     &
                            NSV_DSTBEG, NSV_DSTEND,               &
                            NSV_AERBEG, NSV_AEREND,               &
                            NSV_SLTBEG, NSV_SLTEND,               &
                            NSV_LIMA,NSV_LIMA_BEG,NSV_LIMA_END,   &
                            NSV_LIMA_NC, NSV_LIMA_NR, NSV_LIMA_NI
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA
USE MODD_PARAM_n,     ONLY: CCLOUD, CRAD
USE MODD_PARAM_RAD_n, ONLY: CAOP
USE MODD_RAIN_ICE_DESCR
USE MODD_SALT,        ONLY: LSALT
USE MODD_TIME
!
USE MODE_DUSTOPT
USE MODE_IO_FIELD_WRITE, only: IO_Field_write
USE MODE_ll
use mode_msg
USE MODE_REPRO_SUM,      ONLY : SUM_DD_R2_R1_ll,SUM_DD_R1_ll
!
#ifdef MNH_PGI
USE MODE_PACK_PGI
#endif
USE MODE_SALTOPT
USE MODE_SUM_ll,          ONLY: MIN_ll
USE MODE_SUM2_ll,         ONLY: GMINLOC_ll
USE MODE_THERMO
!
USE MODI_AEROOPT_GET
USE MODI_ECMWF_RADIATION_VERS2
USE MODI_ECRAD_INTERFACE
USE MODD_VAR_ll,      ONLY: IP
!  
IMPLICIT NONE
!
!*       0.1   DECLARATIONS OF DUMMY ARGUMENTS :
!
TYPE(TFILEDATA),  INTENT(IN)         :: TPFILE    ! Output file
LOGICAL, INTENT(IN)                  :: OCLOUD_ONLY! flag for the cloud column
                                                   !    computations only
LOGICAL, INTENT(IN)                  :: OCLEAR_SKY ! 
INTEGER, INTENT(IN)                  :: KDLON   ! number of columns where the
                                                ! radiation calculations are
                                                !       performed
INTEGER, INTENT(IN)                  :: KFLEV   ! number of vertical levels
                                                !    where the radiation
                                                ! calculations are performed
INTEGER, INTENT(IN)                  :: KRAD_DIAG  ! index for the number of
                                                   !  fields in the output
INTEGER, INTENT(IN)                  :: KFLUX   ! number of top and ground 
                                                ! fluxes for the ZFLUX array
INTEGER, INTENT(IN)                  :: KRAD    ! number of satellite radiances
                                                ! for the ZRAD and ZRADCS arrays
INTEGER, INTENT(IN)                  :: KAER    ! number of AERosol classes

INTEGER, INTENT(IN)                  :: KSWB_OLD    ! number of SW band ECMWF 
INTEGER, INTENT(IN)                  :: KSWB_MNH    ! number of SW band ECRAD
INTEGER, INTENT(IN)                  :: KLWB_MNH    ! number of LW band ECRAD
INTEGER, INTENT(IN)                  :: KSTATM  ! index of the standard 
                                                ! atmosphere level just above
                                                !      the model top
INTEGER, INTENT(IN)                  :: KRAD_COLNBR ! factor by which the memory
                                                    ! is split
                                                    !
                                               !Choice of :             
CHARACTER (LEN=*), INTENT (IN)       :: HEFRADL ! 
CHARACTER (LEN=*), INTENT (IN)       :: HEFRADI ! 
CHARACTER (LEN=*), INTENT (IN)       :: HOPWSW !cloud water SW optical properties   
CHARACTER (LEN=*), INTENT (IN)       :: HOPISW !ice water SW optical properties 
CHARACTER (LEN=*), INTENT (IN)       :: HOPWLW !cloud water LW optical properties
CHARACTER (LEN=*), INTENT (IN)       :: HOPILW !ice water  LW optical properties
REAL,               INTENT(IN)       :: PFUDG  ! subgrid cloud inhomogenity factor
REAL, DIMENSION(:,:),     INTENT(IN) :: PCOSZEN ! COS(zenithal solar angle)
REAL,                     INTENT(IN) :: PCORSOL ! SOLar constant CORrection
REAL, DIMENSION(:,:),     INTENT(IN) :: PSEA    ! Land-sea mask
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PDIR_ALB! Surface direct ALBedo
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PSCA_ALB! Surface diffuse ALBedo
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PEMIS   ! Surface IR EMISsivity
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PCLDFR  ! CLouD FRaction
REAL,                     INTENT(IN) :: PCCO2   ! CO2 content
REAL, DIMENSION(:,:),     INTENT(IN) :: PTSRAD  ! RADiative Surface Temperature
REAL, DIMENSION(:,:),     INTENT(IN) :: PSTATM  ! selected standard atmosphere
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT    ! THeta at t
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRT     ! moist variables at t (humidity, cloud water, rain water, ice water)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PPABST  ! pressure at t
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PSVT    ! scalar variable ( C2R2 and C1R3  particle)
!
REAL, DIMENSION(:,:,:),   POINTER    :: POZON   ! OZONE field from clim.
REAL, DIMENSION(:,:,:,:), POINTER    :: PAER    ! AERosols optical thickness from clim. 
REAL, DIMENSION(:,:,:,:), POINTER    :: PDST_WL    ! AERosols Extinction by wavelength . 
REAL, DIMENSION(:,:,:,:), POINTER    :: PAER_CLIM    ! AERosols optical thickness from clim.
                                                ! note : the vertical dimension of 
                                                ! these fields include the "radiation levels"
                                                ! above domain top
                                                ! 
                                                 
REAL, DIMENSION(:,:,:), INTENT(IN)   :: PRHODREF ![kg/m3] air density
REAL, DIMENSION(:,:,:), INTENT(IN)   :: PZZ      ![m] height of layers

INTEGER, DIMENSION(:,:), INTENT(INOUT)  :: KCLEARCOL_TM1 ! trace of cloud/clear col
                                                         ! at the previous radiation step
!                                                 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PDTHRAD ! THeta RADiative Tendancy
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PSRFLWD ! Downward SuRFace LW Flux
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSRFSWD_DIR ! Downward SuRFace SW Flux DIRect 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSRFSWD_DIF ! Downward SuRFace SW Flux DIFfuse 
REAL, DIMENSION(:,:,:),     INTENT(INOUT) :: PSWU ! upward SW Flux 
REAL, DIMENSION(:,:,:),     INTENT(INOUT) :: PSWD ! downward SW Flux 
REAL, DIMENSION(:,:,:),     INTENT(INOUT) :: PLWU ! upward LW Flux 
REAL, DIMENSION(:,:,:),     INTENT(INOUT) :: PLWD ! downward LW Flux 
REAL, DIMENSION(:,:,:),     INTENT(INOUT) :: PDTHRADSW ! dthrad sw 
REAL, DIMENSION(:,:,:),     INTENT(INOUT) :: PDTHRADLW !  dthradsw
REAL, DIMENSION(:,:,:),     INTENT(INOUT) :: PRADEFF ! effective radius
!
!
!*       0.2   DECLARATIONS OF LOCAL VARIABLES
!
LOGICAL                         :: GNOCL     ! .TRUE. when no cloud is present
                                             !     with OCLEAR_SKY .TRUE.
LOGICAL                         :: GAOP      ! .TRUE. when CAOP='EXPL'
LOGICAL, DIMENSION(KDLON,KFLEV) :: GCLOUD    ! .TRUE. for the cloudy columns
LOGICAL, DIMENSION(KFLEV,KDLON) :: GCLOUDT   ! transpose of the GCLOUD array
LOGICAL, DIMENSION(KDLON)       :: GCLEAR_2D ! .TRUE. for the clear-sky columns
LOGICAL, DIMENSION(KDLON,KFLEV) :: GCLEAR    ! .TRUE. for all the levels of the 
                                             !                clear-sky columns
LOGICAL, DIMENSION(KDLON,KSWB_MNH)  :: GCLEAR_SWB! .TRUE. for all the bands of the  
                                             !                clear-sky columns
INTEGER, DIMENSION(:), ALLOCATABLE :: ICLEAR_2D_TM1 !
!
INTEGER :: JI,JJ,JK,JK1,JK2,JKRAD,JALBS! loop indices
!
INTEGER :: IIB           ! I index value of the first inner mass point
INTEGER :: IJB           ! J index value of the first inner mass point
INTEGER :: IKB           ! K index value of the first inner mass point
INTEGER :: IIE           ! I index value of the last inner mass point
INTEGER :: IJE           ! J index value of the last inner mass point
INTEGER :: IKE           ! K index value of the last inner mass point
INTEGER :: IKU           ! array size for the third  index
INTEGER :: IIJ           ! reformatted array index
INTEGER :: IKSTAE        ! level number of the STAndard atmosphere array
INTEGER :: IKUP          ! vertical level above which STAndard atmosphere data
                         ! are filled in
!
INTEGER :: ICLEAR_COL    ! number of    clear-sky columns
INTEGER :: ICLOUD_COL    ! number of    cloudy    columns
INTEGER :: ICLOUD        ! number of levels corresponding of the cloudy columns
INTEGER :: IDIM          ! effective number of columns for which the radiation
                         ! code is run
INTEGER :: INIR          ! index corresponding to NIR fisrt band (in SW)
!
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZTAVE    ! mean-layer temperature
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZTAVE_RAD    ! mean-layer temperature
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZPAVE    ! mean-layer pressure
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZPAVE_RAD    ! mean-layer pressure
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZQSAVE   ! saturation specific humidity
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZQVAVE   ! mean-layer specific humidity
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZQLAVE   ! Liquid water KG/KG
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZQRAVE   ! Rain water  KG/KG
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZQIAVE   ! Ice water Kg/KG
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZQLWC   ! liquid water content kg/m3
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZQRWC   ! Rain water  content kg/m3
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZQIWC   ! ice water content  kg/m3
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCFAVE   ! mean-layer cloud fraction
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZO3AVE   ! mean-layer ozone content 
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZPRES_HL ! half-level pressure
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZT_HL    ! half-level temperature
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZDPRES   ! layer pressure thickness
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCCT_C2R2! Cloud water Concentarion (C2R2)
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCRT_C2R2! Rain water Concentarion (C2R2)
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCIT_C1R3! Ice water Concentarion (C2R2)
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCCT_LIMA! Cloud water Concentration(LIMA)
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCRT_LIMA! Rain water Concentration(LIMA)
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCIT_LIMA! Ice water Concentration(LIMA)
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE :: ZAER     ! aerosol optical thickness
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZALBP    ! spectral surface albedo for direct radiations
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZALBD    ! spectral surface albedo for diffuse radiations 
REAL(KIND=JPRB), DIMENSION (:,:),  ALLOCATABLE :: ZEMIS    ! surface LW  emissivity 
REAL(KIND=JPRB), DIMENSION (:,:), ALLOCATABLE  :: ZEMIW    ! surface LW  WINDOW emissivity
REAL(KIND=JPRB), DIMENSION(:), ALLOCATABLE     :: ZTS      ! reformatted surface PTSRAD array 
REAL(KIND=JPRB), DIMENSION(:), ALLOCATABLE     :: ZLSM     ! reformatted land sea mask
REAL(KIND=JPRB), DIMENSION(:),   ALLOCATABLE   :: ZRMU0    ! Reformatted ZMU0 array
REAL(KIND=JPRB)                     :: ZRII0    ! corrected solar constant
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDTLW    ! LW temperature tendency
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDTSW    ! SW temperature tendency
REAL, DIMENSION(:,:), ALLOCATABLE :: ZNFLW_CS ! CLEAR-SKY LW NET FLUXES
REAL, DIMENSION(:,:), ALLOCATABLE :: ZNFLW    ! TOTAL LW NET FLUXES
REAL, DIMENSION(:,:), ALLOCATABLE :: ZNFSW_CS ! CLEAR-SKY SW NET FLUXES
REAL, DIMENSION(:,:), ALLOCATABLE :: ZNFSW    ! TOTAL SW NET FLUXES
REAL, DIMENSION(:,:), ALLOCATABLE :: ZFLUX_TOP_GND_IRVISNIR ! Top and 
                                                            ! Ground radiative FLUXes
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLUX_SW_DOWN ! DowNward SW Flux profiles
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLUX_SW_UP   ! UPward   SW Flux profiles
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZFLUX_LW      !          LW Flux profiles
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZDTLW_CS ! LW Clear-Sky temp. tendency
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZDTSW_CS ! SW Clear-Sky temp. tendency
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLUX_TOP_GND_IRVISNIR_CS ! Top and
                                                  !  Ground Clear-Sky radiative FLUXes
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZSFSWDIR !surface SW direct flux
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZSFSWDIF !surface SW diffuse flux

REAL, DIMENSION(:),     ALLOCATABLE :: ZPLAN_ALB_VIS, ZPLAN_ALB_NIR
                        ! PLANetary ALBedo in VISible, Near-InfraRed regions
REAL, DIMENSION(:),     ALLOCATABLE :: ZPLAN_TRA_VIS, ZPLAN_TRA_NIR
                        ! PLANetary TRANsmission in VISible, Near-InfraRed regions
REAL, DIMENSION(:),     ALLOCATABLE :: ZPLAN_ABS_VIS, ZPLAN_ABS_NIR
                        ! PLANetary ABSorption in VISible, Near-InfraRed regions
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZEFCL_LWD, ZEFCL_LWU
                        ! EFective  DOWNward and UPward LW nebulosity (equivalent emissivities)
                        ! undefined if RRTM is used for LW
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLWP, ZFIWP
                        ! Liquid and Ice Water Path
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZRADLP, ZRADIP
                        ! Cloud liquid water and ice effective radius
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZEFCL_RRTM, ZCLSW_TOTAL
                        ! effective LW nebulosity ( RRTM case) 
                        ! and SW CLoud fraction for mixed phase clouds
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTAU_TOTAL, ZOMEGA_TOTAL, ZCG_TOTAL
                        ! effective optical thickness, single scattering albedo
                        ! and asymetry factor for mixed phase clouds
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLUX_SW_DOWN_CS, ZFLUX_SW_UP_CS
                        ! Clear-Sky  DowNward and UPward   SW Flux profiles
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZFLUX_LW_CS
                        ! Thicknes of the mesh
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZDZ
!
REAL, DIMENSION(KDLON,KFLEV) :: ZZDTSW ! SW diabatic heating
REAL, DIMENSION(KDLON,KFLEV) :: ZZDTLW ! LW diabatic heating
REAL, DIMENSION(KDLON)       :: ZZTGVIS! SW surface flux in the VIS band
REAL, DIMENSION(KDLON)       :: ZZTGNIR! SW surface flux in the NIR band
REAL, DIMENSION(KDLON)       :: ZZTGIR ! LW surface flux in the IR bands
REAL, DIMENSION(KDLON,SIZE(PSRFSWD_DIR,3)) :: ZZSFSWDIR
!                                      ! SW direct surface flux   
REAL, DIMENSION(KDLON,SIZE(PSRFSWD_DIR,3)) :: ZZSFSWDIF
!                                      ! SW diffuse surface flux   
!
REAL, DIMENSION(KDLON)       :: ZCLOUD ! vertically summed cloud fraction
!
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZEXNT ! Exner function
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2))    :: ZLWD    ! surface Downward LW flux
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PSRFSWD_DIR,3)) :: ZSWDDIR ! surface
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PSRFSWD_DIR,3)) :: ZSWDDIF ! surface Downward SW diffuse flux
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3),KSWB_OLD) :: ZPIZAZ ! Aerosols SSA
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3),KSWB_OLD) :: ZTAUAZ ! Aerosols Optical Detph
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3),KSWB_OLD) :: ZCGAZ  ! Aerosols Asymetric factor
REAL :: ZZTGVISC    ! downward surface SW flux (VIS band) for clear_sky
REAL :: ZZTGNIRC    ! downward surface SW flux (NIR band) for clear_sky
REAL :: ZZTGIRC     ! downward surface LW flux for clear_sky
REAL, DIMENSION(SIZE(PSRFSWD_DIR,3)) :: ZZSFSWDIRC
!                   ! downward surface SW direct flux for clear sky
REAL, DIMENSION(SIZE(PSRFSWD_DIR,3)) :: ZZSFSWDIFC
!                   ! downward surface SW diffuse flux for clear sky
REAL, DIMENSION(KFLEV) :: ZT_CLEAR  ! ensemble mean clear-sky temperature
REAL, DIMENSION(KFLEV) :: ZP_CLEAR  ! ensemble mean clear-sky temperature
REAL, DIMENSION(KFLEV) :: ZQV_CLEAR ! ensemble mean clear-sky specific humidity
REAL, DIMENSION(KFLEV) :: ZOZ_CLEAR ! ensemble mean clear-sky ozone
REAL, DIMENSION(KFLEV) :: ZHP_CLEAR ! ensemble mean clear-sky half-lev. pression
REAL, DIMENSION(KFLEV) :: ZHT_CLEAR ! ensemble mean clear-sky half-lev. temp.
REAL, DIMENSION(KFLEV) :: ZDP_CLEAR ! ensemble mean clear-sky pressure thickness
REAL, DIMENSION(KFLEV,KAER) :: ZAER_CLEAR  ! ensemble mean clear-sky aerosols optical thickness
REAL, DIMENSION(KSWB_MNH)       :: ZALBP_CLEAR ! ensemble mean clear-sky surface albedo (parallel)
REAL, DIMENSION(KSWB_MNH)       :: ZALBD_CLEAR ! ensemble mean clear-sky surface albedo (diffuse)
REAL                        :: ZEMIS_CLEAR ! ensemble mean clear-sky surface emissivity
REAL                        :: ZEMIW_CLEAR ! ensemble mean clear-sky LW window
REAL                        :: ZRMU0_CLEAR ! ensemble mean clear-sky MU0
REAL                        :: ZTS_CLEAR   ! ensemble mean clear-sky surface temperature.
REAL                        :: ZLSM_CLEAR  !  ensemble mean clear-sky land sea-mask  
REAL                        :: ZLAT_CLEAR,ZLON_CLEAR
!
!work arrays
REAL, DIMENSION(:),   ALLOCATABLE :: ZWORK1, ZWORK2, ZWORK3, ZWORK
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK4, ZWORK1AER, ZWORK2AER, ZWORK_GRID
LOGICAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2)) :: ZWORKL
!
!  split arrays used to split the memory required by the ECMWF_radiation 
!  subroutine, the fields have the same meaning as their complete counterpart
!
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZALBP_SPLIT, ZALBD_SPLIT
REAL(KIND=JPRB), DIMENSION(:),     ALLOCATABLE :: ZEMIS_SPLIT, ZEMIW_SPLIT
REAL(KIND=JPRB), DIMENSION(:),     ALLOCATABLE :: ZRMU0_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZCFAVE_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZO3AVE_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZT_HL_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZPRES_HL_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZTAVE_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZPAVE_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE :: ZAER_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZDPRES_SPLIT
REAL(KIND=JPRB), DIMENSION(:),     ALLOCATABLE :: ZLSM_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZQVAVE_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZQSAVE_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZQLAVE_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZQIAVE_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZQRAVE_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZQRWC_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZQLWC_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZQIWC_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:),   ALLOCATABLE :: ZDZ_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCCT_C2R2_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCRT_C2R2_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCIT_C1R3_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCCT_LIMA_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCRT_LIMA_SPLIT
REAL(KIND=JPRB), DIMENSION(:,:), ALLOCATABLE   :: ZCIT_LIMA_SPLIT
REAL(KIND=JPRB), DIMENSION(:),     ALLOCATABLE :: ZTS_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZSFSWDIR_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZSFSWDIF_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZNFLW_CS_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZNFLW_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZNFSW_CS_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZNFSW_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZDTLW_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZDTSW_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLUX_TOP_GND_IRVISNIR_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLUX_SW_DOWN_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLUX_SW_UP_SPLIT
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZFLUX_LW_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZDTLW_CS_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZDTSW_CS_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLUX_TOP_GND_IRVISNIR_CS_SPLIT
REAL, DIMENSION(:),     ALLOCATABLE :: ZPLAN_ALB_VIS_SPLIT
REAL, DIMENSION(:),     ALLOCATABLE :: ZPLAN_ALB_NIR_SPLIT
REAL, DIMENSION(:),     ALLOCATABLE :: ZPLAN_TRA_VIS_SPLIT
REAL, DIMENSION(:),     ALLOCATABLE :: ZPLAN_TRA_NIR_SPLIT
REAL, DIMENSION(:),     ALLOCATABLE :: ZPLAN_ABS_VIS_SPLIT
REAL, DIMENSION(:),     ALLOCATABLE :: ZPLAN_ABS_NIR_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZEFCL_LWD_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZEFCL_LWU_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLWP_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFIWP_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZRADLP_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZRADIP_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZEFCL_RRTM_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZCLSW_TOTAL_SPLIT
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTAU_TOTAL_SPLIT
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZOMEGA_TOTAL_SPLIT
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCG_TOTAL_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLUX_SW_DOWN_CS_SPLIT
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZFLUX_SW_UP_CS_SPLIT
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZFLUX_LW_CS_SPLIT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZPIZA_EQ_TMP        !Single scattering albedo of aerosols (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZIR        !Real part of the aerosol refractive index(lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZII        !Imaginary part of the aerosol refractive index (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCGA_EQ_TMP         !Assymetry factor aerosols            (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZTAUREL_EQ_TMP      !tau/tau_{550} aerosols               (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZPIZA_DST_TMP        !Single scattering albedo of dust (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCGA_DST_TMP         !Assymetry factor dust            (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZTAUREL_DST_TMP      !tau/tau_{550} dust               (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZPIZA_AER_TMP        !Single scattering albedo of aerosol from ORILAM (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCGA_AER_TMP         !Assymetry factor aerosol from ORILAM            (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZTAUREL_AER_TMP      !tau/tau_{550} aerosol from ORILAM               (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZPIZA_SLT_TMP        !Single scattering albedo of sea salt (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCGA_SLT_TMP         !Assymetry factor of sea salt            (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZTAUREL_SLT_TMP      !tau/tau_{550} of sea salt               (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: PAER_AER      !tau/tau_{550} aerosol from ORILAM               (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: PAER_SLT      !tau/tau_{550} sea salt               (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: PAER_DST     !tau/tau_{550} dust               (lon,lat,lev,wvl)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTAU550_EQ_TMP      !tau/tau_{550} aerosols               (lon,lat,lev,wvl)
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE   :: ZPIZA_EQ            !Single scattering albedo of aerosols (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE   :: ZCGA_EQ             !Assymetry factor aerosols            (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE   :: ZTAUREL_EQ          !tau/tau_{550} aerosols               (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE   :: ZPIZA_EQ_SPLIT      !Single scattering albedo of aerosols (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE   :: ZCGA_EQ_SPLIT       !Assymetry factor aerosols            (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(:,:,:), ALLOCATABLE   :: ZTAUREL_EQ_SPLIT    !tau/tau_{550} aerosols               (points,lev,wvl)
REAL, DIMENSION(KFLEV,KSWB_OLD)           :: ZPIZA_EQ_CLEAR      !Single scattering albedo of aerosols (lev,wvl)
REAL, DIMENSION(KFLEV,KSWB_OLD)           :: ZCGA_EQ_CLEAR       !Assymetry factor aerosols            (lev,wvl)
REAL, DIMENSION(KFLEV,KSWB_OLD)           :: ZTAUREL_EQ_CLEAR    !tau/tau_{550} aerosols               (lev,wvl)
INTEGER                               :: WVL_IDX              !Counter for wavelength

!
INTEGER  :: JI_SPLIT          ! loop on the split array
INTEGER  :: INUM_CALL         ! number of CALL of the radiation scheme
INTEGER  :: IDIM_EFF          ! effective number of air-columns to compute
INTEGER  :: IDIM_RESIDUE      ! number of remaining air-columns to compute
INTEGER  :: IBEG, IEND        ! auxiliary indices
!
!
REAL, DIMENSION(SIZE(PDTHRAD,1),SIZE(PDTHRAD,2),SIZE(PDTHRAD,3)) &
     :: ZDTRAD_LW! LW temperature tendency
REAL, DIMENSION(SIZE(PDTHRAD,1),SIZE(PDTHRAD,2),SIZE(PDTHRAD,3)) &
     :: ZDTRAD_SW! SW temperature tendency
INTEGER             :: ILUOUT       ! Logical unit number for output-listing
INTEGER             :: IRESP        ! Return code of FM routines
REAL, DIMENSION(SIZE(PDTHRAD,1),SIZE(PDTHRAD,2),SIZE(PDTHRAD,3)) &
     :: ZSTORE_3D, ZSTORE_3D2! 3D work array for storage
REAL, DIMENSION(SIZE(PDTHRAD,1),SIZE(PDTHRAD,2)) &
     :: ZSTORE_2D   ! 2D work array for storage!
INTEGER                         :: JBAND       ! Solar band index
CHARACTER (LEN=4), DIMENSION(KSWB_OLD) :: YBAND_NAME  ! Solar band name
CHARACTER (LEN=2)               :: YDIR        ! Type of the data field
!
INTEGER :: ISWB ! number of SW spectral bands (between radiations and surface schemes)
INTEGER :: JSWB ! loop on SW spectral bands
INTEGER :: JAE  ! loop on aerosol class
TYPE(TFIELDMeTaDATA) :: TZFIELD2D, TZFIELD3D
!
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZDZPABST
REAL :: ZMINVAL
INTEGER, DIMENSION(3) :: IMINLOC
INTEGER :: IINFO_ll
LOGICAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2)) :: GCLOUD_SURF
!
REAL(KIND=JPRB), DIMENSION(:),   ALLOCATABLE :: ZLON,ZLAT
REAL(KIND=JPRB), DIMENSION(:),   ALLOCATABLE :: ZLON_SPLIT,ZLAT_SPLIT
!
INTEGER                            :: ICLEAR_COL_ll
INTEGER, DIMENSION(:), ALLOCATABLE :: INDEX_ICLEAR_COL
REAL, DIMENSION(KFLEV)             :: ZT_CLEAR_DD  ! ensemble mean clear-sky temperature
REAL                               :: ZCLEAR_COL_ll , ZDLON_ll
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!*       1.    COMPUTE DIMENSIONS OF ARRAYS AND OTHER INDICES
!              ----------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)  ! this definition must be coherent with
                                      ! the one used in ini_radiations routine
IKU = SIZE(PTHT,3)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
!
IKSTAE = SIZE(PSTATM,1)
IKUP   = IKE-JPVEXT+1
! 
ISWB   = SIZE(PSRFSWD_DIR,3)
!
!-------------------------------------------------------------------------------
!*       1.1   CHECK PRESSURE DECREASING
!              -------------------------
ZDZPABST(:,:,1:IKU-1) = PPABST(:,:,1:IKU-1) - PPABST(:,:,2:IKU)
ZDZPABST(:,:,IKU) = ZDZPABST(:,:,IKU-1)
!
ZMINVAL=MIN_ll(ZDZPABST,IINFO_ll)
!
IF ( ZMINVAL <= 0.0 ) THEN
   ILUOUT = TLUOUT%NLU
   IMINLOC=GMINLOC_ll( ZDZPABST )
   WRITE(ILUOUT,*) ' radiation.f90 STOP :: SOMETHING WRONG WITH PRESSURE , ZDZPABST <= 0.0 '  
   WRITE(ILUOUT,*) ' radiation :: ZDZPABST ', ZMINVAL,' located at ',   IMINLOC
   FLUSH(unit=ILUOUT)
   call Print_msg( NVERB_FATAL, 'GEN', 'RADIATIONS', 'something wrong with pressure: ZDZPABST <= 0.0' )

ENDIF
!------------------------------------------------------------------------------
ALLOCATE(ZLAT(KDLON))
ALLOCATE(ZLON(KDLON))
IF(LCARTESIAN) THEN
  ZLAT(:) = XLAT0*(XPI/180.)
  ZLON(:) = XLON0*(XPI/180.)
ELSE
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZLAT(IIJ) =  XLAT(JI,JJ)*(XPI/180.)
        ZLON(IIJ) =  XLON(JI,JJ)*(XPI/180.)
    END DO
  END DO
END IF
!-------------------------------------------------------------------------------
!
!*       2.    INITIALIZES THE MEAN-LAYER VARIABLES
!              ------------------------------------
!
ZEXNT(:,:,:)= ( PPABST(:,:,:)/XP00 ) ** (XRD/XCPD)
!
! Columns where radiation is computed are put on a single line
ALLOCATE(ZTAVE(KDLON,KFLEV))
ALLOCATE(ZQVAVE(KDLON,KFLEV))
ALLOCATE(ZQLAVE(KDLON,KFLEV))
ALLOCATE(ZQIAVE(KDLON,KFLEV))
ALLOCATE(ZCFAVE(KDLON,KFLEV))
ALLOCATE(ZQRAVE(KDLON,KFLEV))
ALLOCATE(ZQLWC(KDLON,KFLEV))
ALLOCATE(ZQIWC(KDLON,KFLEV))
ALLOCATE(ZQRWC(KDLON,KFLEV))
ALLOCATE(ZDZ(KDLON,KFLEV))
!
ZQVAVE(:,:) = 0.0
ZQLAVE(:,:) = 0.0
ZQIAVE(:,:) = 0.0
ZQRAVE(:,:) = 0.0
ZCFAVE(:,:) = 0.0
ZQLWC(:,:) = 0.0
ZQIWC(:,:) = 0.0
ZQRWC(:,:) = 0.0
ZDZ(:,:)=0.0
!
!COMPUTE THE MESH SIZE
DO JK=IKB,IKE
  JKRAD = JK-JPVEXT
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZDZ(IIJ,JKRAD)  =  PZZ(JI,JJ,JK+1) - PZZ(JI,JJ,JK)
      ZTAVE(IIJ,JKRAD)  = PTHT(JI,JJ,JK)*ZEXNT(JI,JJ,JK) ! Conversion potential temperature -> actual temperature
    END DO
  END DO
END DO
!
!  Check if the humidity mixing ratio is available
!
IF( SIZE(PRT(:,:,:,:),4) >= 1 ) THEN
  DO JK=IKB,IKE
    JKRAD = JK-JPVEXT
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZQVAVE(IIJ,JKRAD) =MAX(0., PRT(JI,JJ,JK,1))
      END DO
    END DO
  END DO
END IF
!
!  Check if the cloudwater mixing ratio is available
!
IF( SIZE(PRT(:,:,:,:),4) >= 2 ) THEN
  DO JK=IKB,IKE
    JKRAD = JK-JPVEXT
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZQLAVE(IIJ,JKRAD) = MAX(0.,PRT(JI,JJ,JK,2))
        ZQLWC(IIJ,JKRAD) = MAX(0.,PRT(JI,JJ,JK,2)*PRHODREF(JI,JJ,JK))
        ZCFAVE(IIJ,JKRAD) = PCLDFR(JI,JJ,JK)
      END DO
    END DO
  END DO
END IF
!
!  Check if the rainwater mixing ratio is available
!
IF( SIZE(PRT(:,:,:,:),4) >= 3 ) THEN
  DO JK=IKB,IKE
    JKRAD = JK-JPVEXT
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZQRWC(IIJ,JKRAD) = MAX(0.,PRT(JI,JJ,JK,3)*PRHODREF(JI,JJ,JK))
        ZQRAVE(IIJ,JKRAD) = MAX(0.,PRT(JI,JJ,JK,3))
      END DO
    END DO
  END DO
END IF
!
!  Check if the cloudice mixing ratio is available
!
IF( SIZE(PRT(:,:,:,:),4) >= 4 ) THEN
  DO JK=IKB,IKE
    JKRAD = JK-JPVEXT
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZQIWC(IIJ,JKRAD) = MAX(0.,PRT(JI,JJ,JK,4)*PRHODREF(JI,JJ,JK))
!        ZQIAVE(IIJ,JKRAD) = MAX( PRT(JI,JJ,JK,4)-XRTMIN(4),0.0 )
        ZQIAVE(IIJ,JKRAD) = MAX( PRT(JI,JJ,JK,4),0.0 )
      END DO
    END DO
  END DO
END IF
!
!  Standard atmosphere extension
!
DO JK=IKUP,KFLEV
  JK1 = (KSTATM-1)+(JK-IKUP)
  JK2 = JK1+1
  ZTAVE(:,JK)  = 0.5*( PSTATM(JK1,3)+PSTATM(JK2,3) )
  ZQVAVE(:,JK) = 0.5*( PSTATM(JK1,5)/PSTATM(JK1,4)+   &
                 PSTATM(JK2,5)/PSTATM(JK2,4)    )
END DO
!
!        2.1 pronostic water concentation fields (C2R2 coupling) 
!
IF( NSV_C2R2 /= 0 ) THEN
  ALLOCATE (ZCCT_C2R2(KDLON, KFLEV))
  ALLOCATE (ZCRT_C2R2(KDLON, KFLEV))
  ZCCT_C2R2(:, :) = 0.
  ZCRT_C2R2 (:,:) = 0.
  DO JK=IKB,IKE
    JKRAD = JK-JPVEXT
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZCCT_C2R2 (IIJ,JKRAD) = MAX(0.,PSVT(JI,JJ,JK,NSV_C2R2BEG+1))
        ZCRT_C2R2 (IIJ,JKRAD) = MAX(0.,PSVT(JI,JJ,JK,NSV_C2R2BEG+2))
      END DO
    END DO
  END DO
ELSE 
  ALLOCATE (ZCCT_C2R2(0,0))
  ALLOCATE (ZCRT_C2R2(0,0))
END IF
!
IF( NSV_C1R3 /= 0 ) THEN
  ALLOCATE (ZCIT_C1R3(KDLON, KFLEV))
  ZCIT_C1R3 (:,:) = 0.
  DO JK=IKB,IKE
    JKRAD = JK-JPVEXT
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZCIT_C1R3 (IIJ,JKRAD) = MAX(0.,PSVT(JI,JJ,JK,NSV_C1R3BEG))
      END DO
    END DO
  END DO
ELSE 
  ALLOCATE (ZCIT_C1R3(0,0))
END IF
!
!
!        2.1*bis pronostic water concentation fields (LIMA coupling) 
!
IF( CCLOUD == 'LIMA' ) THEN
   ALLOCATE (ZCCT_LIMA(KDLON, KFLEV))
   ALLOCATE (ZCRT_LIMA(KDLON, KFLEV))
   ALLOCATE (ZCIT_LIMA(KDLON, KFLEV))
   ZCCT_LIMA(:, :) = 0.
   ZCRT_LIMA (:,:) = 0.
   ZCIT_LIMA (:,:) = 0.
   DO JK=IKB,IKE
      JKRAD = JK-JPVEXT
      DO JJ=IJB,IJE
         DO JI=IIB,IIE
            IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
            IF (NMOM_C.GE.2) ZCCT_LIMA(IIJ,JKRAD) = MAX(0.,PSVT(JI,JJ,JK,NSV_LIMA_NC))
            IF (NMOM_R.GE.2) ZCRT_LIMA(IIJ,JKRAD) = MAX(0.,PSVT(JI,JJ,JK,NSV_LIMA_NR))
            IF (NMOM_I.GE.2) ZCIT_LIMA(IIJ,JKRAD) = MAX(0.,PSVT(JI,JJ,JK,NSV_LIMA_NI))
         END DO
      END DO
   END DO
END IF
!
!-------------------------------------------------------------------------------
!
!*       3.    INITIALIZES THE HALF-LEVEL VARIABLES
!  	           ------------------------------------
!
ALLOCATE(ZPRES_HL(KDLON,KFLEV+1))
ALLOCATE(ZT_HL(KDLON,KFLEV+1))
!
DO JK=IKB,IKE+1
  JKRAD = JK-JPVEXT
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZPRES_HL(IIJ,JKRAD) = XP00 * (0.5*(ZEXNT(JI,JJ,JK)+ZEXNT(JI,JJ,JK-1)))**(XCPD/XRD)
    END DO
  END DO
END DO

!  Standard atmosphere extension - pressure
!* begining at ikup+1 level allows to use a model domain higher than 50km
!
DO JK=IKUP+1,KFLEV+1
  JK1 = (KSTATM-1)+(JK-IKUP)
  ZPRES_HL(:,JK) = PSTATM(JK1,2)*100.0 ! mb -> Pa
END DO
!
!  Surface temperature at the first level
!  and surface radiative temperature
ALLOCATE(ZTS(KDLON))
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZT_HL(IIJ,1) = PTSRAD(JI,JJ)
    ZTS(IIJ) = PTSRAD(JI,JJ)
  END DO
END DO
!
!  Temperature at half levels
!
ZT_HL(:,2:IKE-JPVEXT) = 0.5*(ZTAVE(:,1:IKE-JPVEXT-1)+ZTAVE(:,2:IKE-JPVEXT))
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZT_HL(IIJ,IKE-JPVEXT+1)  =  0.5*PTHT(JI,JJ,IKE  )*ZEXNT(JI,JJ,IKE  ) &
         + 0.5*PTHT(JI,JJ,IKE+1)*ZEXNT(JI,JJ,IKE+1)
  END DO
END DO
!
!  Standard atmosphere extension - temperature
!* begining at ikup+1 level allows to use a model domain higher than 50km
!
DO JK=IKUP+1,KFLEV+1
  JK1 = (KSTATM-1)+(JK-IKUP)
  ZT_HL(:,JK) = PSTATM(JK1,3)
END DO
!
!mean layer pressure and layer differential pressure (from half level variables)
!
ALLOCATE(ZPAVE(KDLON,KFLEV))
ALLOCATE(ZDPRES(KDLON,KFLEV))
DO JKRAD=1,KFLEV
  ZPAVE(:,JKRAD)=0.5*(ZPRES_HL(:,JKRAD)+ZPRES_HL(:,JKRAD+1))
  ZDPRES(:,JKRAD)=ZPRES_HL(:,JKRAD)-ZPRES_HL(:,JKRAD+1)
END DO
!-----------------------------------------------------------------------
!*       4.    INITIALIZES THE AEROSOLS and OZONE PROFILES from climatology
!	           -------------------------------------------
!
!        4.1    AEROSOL optical thickness
! EXPL -> defined online, otherwise climatology
IF (CAOP=='EXPL') THEN
   GAOP = .TRUE.
ELSE
   GAOP = .FALSE.
ENDIF
!
IF (CAOP=='EXPL') THEN
   ALLOCATE(ZPIZA_EQ_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD))
   ALLOCATE(ZCGA_EQ_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD))
   ALLOCATE(ZTAUREL_EQ_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD))

   ALLOCATE(ZPIZA_DST_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD))
   ALLOCATE(ZCGA_DST_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD))
   ALLOCATE(ZTAUREL_DST_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD)) 
   ALLOCATE(PAER_DST(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3)))

   ALLOCATE(ZPIZA_AER_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD))
   ALLOCATE(ZCGA_AER_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD))
   ALLOCATE(ZTAUREL_AER_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD))
   ALLOCATE(PAER_AER(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3)))

   ALLOCATE(ZPIZA_SLT_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD))
   ALLOCATE(ZCGA_SLT_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD))
   ALLOCATE(ZTAUREL_SLT_TMP(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3),KSWB_OLD))
   ALLOCATE(PAER_SLT(SIZE(PAER,1),SIZE(PAER,2),SIZE(PAER,3)))
   

   ALLOCATE(ZII(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),KSWB_OLD))
   ALLOCATE(ZIR(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),KSWB_OLD))

  ZPIZA_EQ_TMP = 0.
  ZCGA_EQ_TMP = 0.
  ZTAUREL_EQ_TMP = 0.

  ZPIZA_DST_TMP = 0.
  ZCGA_DST_TMP = 0.
  ZTAUREL_DST_TMP = 0

  ZPIZA_SLT_TMP = 0.
  ZCGA_SLT_TMP = 0.
  ZTAUREL_SLT_TMP = 0

  ZPIZA_AER_TMP = 0.
  ZCGA_AER_TMP = 0.
  ZTAUREL_AER_TMP = 0

  PAER_DST=0.
  PAER_SLT=0.
  PAER_AER=0.
  
 IF (LORILAM) THEN
   CALL AEROOPT_GET(                             &
        PSVT(IIB:IIE,IJB:IJE,:,NSV_AERBEG:NSV_AEREND)        &  !I [ppv]  aerosols concentration
        ,PZZ(IIB:IIE,IJB:IJE,:)                   &  !I [m] height of layers
        ,PRHODREF(IIB:IIE,IJB:IJE,:)              &  !I [kg/m3] density of air
        ,ZPIZA_AER_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,:)   &  !O [-] single scattering albedo of aerosols
        ,ZCGA_AER_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,:)    &  !O [-] assymetry factor for aerosols
        ,ZTAUREL_AER_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,:) &  !O [-] opt.depth(wvl=lambda)/opt.depth(wvl=550nm)
        ,PAER_AER(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT)            &  !O [-] optical depth of aerosols at wvl=550nm
        ,KSWB_OLD                                    &  !I |nbr] number of shortwave bands
        ,ZIR(IIB:IIE,IJB:IJE,:,:) &  !O [-] opt.depth(wvl=lambda)/opt.depth(wvl=550nm)
        ,ZII(IIB:IIE,IJB:IJE,:,:) &  !O [-] opt.depth(wvl=lambda)/opt.depth(wvl=550nm)
        )
 ENDIF
 IF(LDUST) THEN
   CALL DUSTOPT_GET(                             &
        PSVT(IIB:IIE,IJB:IJE,:,NSV_DSTBEG:NSV_DSTEND)        &  !I [ppv] Dust scalar concentration
        ,PZZ(IIB:IIE,IJB:IJE,:)                   &  !I [m] height of layers
        ,PRHODREF(IIB:IIE,IJB:IJE,:)              &  !I [kg/m3] density of air
        ,ZPIZA_DST_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,:)   &  !O [-] single scattering albedo of dust
        ,ZCGA_DST_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,:)    &  !O [-] assymetry factor for dust
        ,ZTAUREL_DST_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,:) &  !O [-] opt.depth(wvl=lambda)/opt.depth(wvl=550nm)
        ,PAER_DST(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT)            &  !O [-] optical depth of dust at wvl=550nm
        ,KSWB_OLD                                   &  !I |nbr] number of shortwave bands
        )
   DO WVL_IDX=1,KSWB_OLD
     PDST_WL(:,:,:,WVL_IDX) = ZTAUREL_DST_TMP(:,:,:,WVL_IDX)* PAER(:,:,:,3)             
   ENDDO
 ENDIF
 IF(LSALT) THEN
   CALL SALTOPT_GET(                             &
        PSVT(IIB:IIE,IJB:IJE,:,NSV_SLTBEG:NSV_SLTEND)        &  !I [ppv] sea salt scalar concentration
        ,PZZ(IIB:IIE,IJB:IJE,:)                   &  !I [m] height of layers
        ,PRHODREF(IIB:IIE,IJB:IJE,:)              &  !I [kg/m3] density of air
        ,PTHT(IIB:IIE,IJB:IJE,:)                  &  !I [K] potential temperature
        ,PPABST(IIB:IIE,IJB:IJE,:)                &  !I [hPa] pressure
        ,PRT(IIB:IIE,IJB:IJE,:,:)                 &  !I [kg/kg] water mixing ratio
        ,ZPIZA_SLT_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,:)   &  !O [-] single scattering albedo of sea salt
        ,ZCGA_SLT_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,:)    &  !O [-] assymetry factor for sea salt
        ,ZTAUREL_SLT_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,:) &  !O [-] opt.depth(wvl=lambda)/opt.depth(wvl=550nm)
        ,PAER_SLT(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT)            &  !O [-] optical depth of sea salt at wvl=550nm
        ,KSWB_OLD                                    &  !I |nbr] number of shortwave bands
        )
 ENDIF

 ZTAUREL_EQ_TMP(:,:,:,:)=ZTAUREL_DST_TMP(:,:,:,:)+ZTAUREL_AER_TMP(:,:,:,:)+ZTAUREL_SLT_TMP(:,:,:,:)
 
 PAER(:,:,:,2)=PAER_SLT(:,:,:)
 PAER(:,:,:,3)=PAER_DST(:,:,:)
 PAER(:,:,:,4)=PAER_AER(:,:,:)


 WHERE (ZTAUREL_EQ_TMP(:,:,:,:).GT.0.0)
  ZPIZA_EQ_TMP(:,:,:,:)=(ZPIZA_DST_TMP(:,:,:,:)*ZTAUREL_DST_TMP(:,:,:,:)+&
                    ZPIZA_AER_TMP(:,:,:,:)*ZTAUREL_AER_TMP(:,:,:,:)+&
                    ZPIZA_SLT_TMP(:,:,:,:)*ZTAUREL_SLT_TMP(:,:,:,:))/&
                    ZTAUREL_EQ_TMP(:,:,:,:) 
 END WHERE
 WHERE ((ZTAUREL_EQ_TMP(:,:,:,:).GT.0.0).AND.(ZPIZA_EQ_TMP(:,:,:,:).GT.0.0))
  ZCGA_EQ_TMP(:,:,:,:)=(ZPIZA_DST_TMP(:,:,:,:)*ZTAUREL_DST_TMP(:,:,:,:)*ZCGA_DST_TMP(:,:,:,:)+&
                   ZPIZA_AER_TMP(:,:,:,:)*ZTAUREL_AER_TMP(:,:,:,:)*ZCGA_AER_TMP(:,:,:,:)+&
                   ZPIZA_SLT_TMP(:,:,:,:)*ZTAUREL_SLT_TMP(:,:,:,:)*ZCGA_SLT_TMP(:,:,:,:))/&
                   (ZTAUREL_EQ_TMP(:,:,:,:)*ZPIZA_EQ_TMP(:,:,:,:))
 END WHERE

 ZTAUREL_EQ_TMP(:,:,:,:)=max(1.E-8,ZTAUREL_EQ_TMP(:,:,:,:))
 ZCGA_EQ_TMP(:,:,:,:)=max(1.E-8,ZCGA_EQ_TMP(:,:,:,:))
 ZPIZA_EQ_TMP(:,:,:,:)=max(1.E-8,ZPIZA_EQ_TMP(:,:,:,:))
 PAER(:,:,:,3)=max(1.E-8,PAER(:,:,:,3))
 ZPIZA_EQ_TMP(:,:,:,:)=min(0.99,ZPIZA_EQ_TMP(:,:,:,:))


ENDIF      
!
! Computes SSA, optical depth and assymetry factor for clear sky (aerosols)
ZTAUAZ(:,:,:,:) = 0.
ZPIZAZ(:,:,:,:) = 0.
ZCGAZ(:,:,:,:)  = 0.
DO WVL_IDX=1,KSWB_OLD
 DO JAE=1,KAER
      !Special optical properties for dust
      IF (CAOP=='EXPL'.AND.(JAE==3)) THEN
      !Ponderation of aerosol optical in case of explicit optical factor
      !ti
        ZTAUAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX)= ZTAUAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX) + &
                                                 PAER(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,JAE) * &
                                       ZTAUREL_EQ_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,WVL_IDX) 
      !wi*ti
        ZPIZAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX)= ZPIZAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX) + &
                                  PAER(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,JAE)                * &
                                  ZTAUREL_EQ_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,WVL_IDX) * &
                                  ZPIZA_EQ_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,WVL_IDX)
      !wi*ti*gi
        ZCGAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX) = ZCGAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX) + &
                                 PAER(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,JAE)                * &
                                 ZTAUREL_EQ_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,WVL_IDX) * &
                                 ZPIZA_EQ_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,WVL_IDX)   * &
                                 ZCGA_EQ_TMP(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,WVL_IDX)
      ELSE

      !Ponderation of aerosol optical properties 
      !ti
        ZTAUAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX)=ZTAUAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX)+&
             PAER(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,JAE) * RTAUA(WVL_IDX,JAE)
      !wi*ti
        ZPIZAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX)=ZPIZAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX)+&
                                               PAER(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,JAE) *&
                                        RTAUA(WVL_IDX,JAE)*RPIZA(WVL_IDX,JAE)
      !wi*ti*gi
        ZCGAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX) =  ZCGAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX) +&
                                               PAER(IIB:IIE,IJB:IJE,IKB-JPVEXT:IKE-JPVEXT,JAE)   *&
                        RTAUA(WVL_IDX,JAE)*RPIZA(WVL_IDX,JAE)*RCGA(WVL_IDX,JAE)
           ENDIF
 ENDDO
! assymetry factor:

ZCGAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX) = ZCGAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX)  / &
                                   ZPIZAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX)
! SSA:
ZPIZAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX) = ZPIZAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX) / &
                                    ZTAUAZ(IIB:IIE,IJB:IJE,IKB:IKE,WVL_IDX)
ENDDO
!

!
ALLOCATE(ZAER(KDLON,KFLEV,KAER))
! Aerosol classes
! 1=Continental   2=Maritime   3=Desert     4=Urban     5=Volcanic 6=Stratos.Bckgnd
! Loaded from climatology
DO JJ=IJB,IJE
   DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZAER (IIJ,:,:) = PAER_CLIM  (JI,JJ,:,:)
   END DO
END DO
IF ((CAOP=='EXPL') .AND. LDUST ) THEN
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZAER (IIJ,:,3) = PAER  (JI,JJ,:,3)
    END DO
  END DO
END IF
IF ((CAOP=='EXPL') .AND. LSALT ) THEN
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZAER (IIJ,:,2) = PAER  (JI,JJ,:,2)
    END DO
  END DO
END IF
IF ((CAOP=='EXPL') .AND. LORILAM ) THEN
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZAER (IIJ,:,4) = PAER  (JI,JJ,:,4)
    END DO
  END DO
END IF
!
ALLOCATE(ZPIZA_EQ(KDLON,KFLEV,KSWB_OLD))
ALLOCATE(ZCGA_EQ(KDLON,KFLEV,KSWB_OLD))
ALLOCATE(ZTAUREL_EQ(KDLON,KFLEV,KSWB_OLD))
IF(CAOP=='EXPL')THEN
    !Transform from vector of type #lon #lat #lev #wvl
    !to vectors of type #points, #levs, #wavelengths
  DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZPIZA_EQ(IIJ,:,:) = ZPIZA_EQ_TMP(JI,JJ,:,:)
    ZCGA_EQ(IIJ,:,:)= ZCGA_EQ_TMP(JI,JJ,:,:)
    ZTAUREL_EQ(IIJ,:,:)=ZTAUREL_EQ_TMP(JI,JJ,:,:)
  END DO
  END DO
  DEALLOCATE(ZPIZA_EQ_TMP)
  DEALLOCATE(ZCGA_EQ_TMP)
  DEALLOCATE(ZTAUREL_EQ_TMP)
  DEALLOCATE(ZPIZA_DST_TMP)
  DEALLOCATE(ZCGA_DST_TMP)
  DEALLOCATE(ZTAUREL_DST_TMP)  
  DEALLOCATE(ZPIZA_AER_TMP)
  DEALLOCATE(ZCGA_AER_TMP)
  DEALLOCATE(ZTAUREL_AER_TMP)
  DEALLOCATE(ZPIZA_SLT_TMP)
  DEALLOCATE(ZCGA_SLT_TMP)
  DEALLOCATE(ZTAUREL_SLT_TMP)
  DEALLOCATE(PAER_DST)
  DEALLOCATE(PAER_AER)
  DEALLOCATE(PAER_SLT)
  DEALLOCATE(ZIR)
  DEALLOCATE(ZII)
END IF


!
!      4.2   OZONE content 
!
ALLOCATE(ZO3AVE(KDLON,KFLEV))
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZO3AVE(IIJ,:)  = POZON (JI,JJ,:)           
  END DO
END DO
#ifdef MNH_ECRAD
#if ( VER_ECRAD == 140 )
POZON = POZON
#endif
#endif
!
!-------------------------------------------------------------------------------
!
!*       5.    CALLS THE E.C.M.W.F. RADIATION CODE
!	           -----------------------------------
!
!
!*       5.1   INITIALIZES 2D AND SURFACE FIELDS
!
ALLOCATE(ZRMU0(KDLON))
ALLOCATE(ZLSM(KDLON))
! 
ALLOCATE(ZALBP(KDLON,KSWB_MNH))
ALLOCATE(ZALBD(KDLON,KSWB_MNH))
!
ALLOCATE(ZEMIS(KDLON,KLWB_MNH))
ALLOCATE(ZEMIW(KDLON,KLWB_MNH))
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZEMIS(IIJ,:)   = PEMIS(JI,JJ,:)
    ZRMU0(IIJ)    = PCOSZEN(JI,JJ)
    ZLSM(IIJ)     = 1.0 - PSEA(JI,JJ)  
  END DO
END DO  
!
! spectral albedo
!
IF ( SIZE(PDIR_ALB,3)==1 ) THEN
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      !  sw direct and diffuse albedos
      ZALBP(IIJ,:)  = PDIR_ALB(JI,JJ,1)
      ZALBD(IIJ,:)  = PSCA_ALB(JI,JJ,1)
      !
    END DO
  END DO
ELSE  
  DO JK=1, SIZE(PDIR_ALB,3)
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
         IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
         !  sw direct and diffuse albedos
         ZALBP(IIJ,JK)  = PDIR_ALB(JI,JJ,JK)
         ZALBD(IIJ,JK)  = PSCA_ALB(JI,JJ,JK)
       ENDDO
     END DO
   ENDDO  
END IF
!
!
! LW emissivity
ZEMIW(:,:)= ZEMIS(:,:)
!
!solar constant
ZRII0= PCORSOL*XI0  ! solar constant multiplied by seasonal variations due to Earth-Sun distance
!
!
!*       5.2   ACCOUNTS FOR THE CLEAR-SKY APPROXIMATION
!
!  Performs the horizontal average of the fields when no cloud
!
ZCLOUD(:) = SUM( ZCFAVE(:,:),DIM=2 ) ! one where no cloud on the vertical
!
! MODIF option CLLY      
ALLOCATE ( ICLEAR_2D_TM1(KDLON) )
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ICLEAR_2D_TM1(IIJ) = KCLEARCOL_TM1(JI,JJ)
  END DO
END DO
!
IF(OCLOUD_ONLY .OR. OCLEAR_SKY) THEN
  !
  GCLEAR_2D(:) = .TRUE.
  WHERE( (ZCLOUD(:) > 0.0) .OR. (ICLEAR_2D_TM1(:)==0) )  ! FALSE on cloudy columns
    GCLEAR_2D(:) = .FALSE.
  END WHERE
  !
  ICLEAR_COL = COUNT( GCLEAR_2D(:) )  ! number of clear sky columns
  !
  ALLOCATE(INDEX_ICLEAR_COL(ICLEAR_COL))
  IIJ = 0
  DO JI=1,KDLON
     IF ( GCLEAR_2D(JI) ) THEN
         IIJ = IIJ + 1
         INDEX_ICLEAR_COL(IIJ) = JI
     END IF
  END DO

  IF( ICLEAR_COL == KDLON ) THEN ! No cloud case so only the mean clear-sky
!!$    GCLEAR_2D(1) = .FALSE.       !           column is selected
!!$    ICLEAR_COL = KDLON-1
    GNOCL = .TRUE.               ! TRUE if no cloud at all
  ELSE
    GNOCL = .FALSE.
  END IF

  GCLEAR(:,:) = SPREAD( GCLEAR_2D(:),DIM=2,NCOPIES=KFLEV )  ! vertical extension of clear columns 2D map
  ICLOUD_COL = KDLON - ICLEAR_COL                           ! number of  cloudy columns
!
  ZCLEAR_COL_ll = REAL(ICLEAR_COL)
  CALL REDUCESUM_ll(ZCLEAR_COL_ll,IINFO_ll)
  !ZDLON_ll = KDLON
  !CALL REDUCESUM_ll(ZDLON_ll,IINFO_ll)

  !IF (IP == 1 )  
  !print*,",RADIATIOn COULD_ONLY=OCLOUD_ONLY,OCLEAR_SKY,ZCLEAR_COL_ll,ICLEAR_COL,ICLOUD_COL,KDON,ZDLON_ll,GNOCL=", &
  !     OCLOUD_ONLY,OCLEAR_SKY,ZCLEAR_COL_ll,ICLEAR_COL,ICLOUD_COL,KDLON,ZDLON_ll,GNOCL
!
!!$  IF( ICLEAR_COL /=0 ) THEN ! at least one clear-sky column exists -> average profiles on clear columns
  IF( ZCLEAR_COL_ll /= 0.0  ) THEN ! at least one clear-sky column exists -> average profiles on clear columns
    ZT_CLEAR(:)  = SUM_DD_R2_R1_ll(ZTAVE(INDEX_ICLEAR_COL(:),:)) / ZCLEAR_COL_ll
    ZP_CLEAR(:)  = SUM_DD_R2_R1_ll(ZPAVE(INDEX_ICLEAR_COL(:),:)) / ZCLEAR_COL_ll
    ZQV_CLEAR(:) = SUM_DD_R2_R1_ll(REAL(ZQVAVE(INDEX_ICLEAR_COL(:),:))) / ZCLEAR_COL_ll
    ZOZ_CLEAR(:) = SUM_DD_R2_R1_ll(REAL(ZO3AVE(INDEX_ICLEAR_COL(:),:))) / ZCLEAR_COL_ll
    ZDP_CLEAR(:) = SUM_DD_R2_R1_ll(REAL(ZDPRES(INDEX_ICLEAR_COL(:),:))) / ZCLEAR_COL_ll
 
    DO JK1=1,KAER
      ZAER_CLEAR(:,JK1) = SUM_DD_R2_R1_ll(REAL(ZAER(INDEX_ICLEAR_COL(:),:,JK1))) / ZCLEAR_COL_ll
    END DO
    !Get an average value for the clear column
    IF(CAOP=='EXPL')THEN
       DO WVL_IDX=1,KSWB_OLD
          ZPIZA_EQ_CLEAR(:,WVL_IDX)   = SUM_DD_R2_R1_ll(REAL(ZPIZA_EQ(  INDEX_ICLEAR_COL(:),:,WVL_IDX))) / ZCLEAR_COL_ll
          ZCGA_EQ_CLEAR(:,WVL_IDX)    = SUM_DD_R2_R1_ll(REAL(ZCGA_EQ(   INDEX_ICLEAR_COL(:),:,WVL_IDX))) / ZCLEAR_COL_ll
          ZTAUREL_EQ_CLEAR(:,WVL_IDX) = SUM_DD_R2_R1_ll(REAL(ZTAUREL_EQ(INDEX_ICLEAR_COL(:),:,WVL_IDX))) / ZCLEAR_COL_ll
       ENDDO
    ENDIF   
    !
    ZHP_CLEAR(1:KFLEV) = SUM_DD_R2_R1_ll(REAL(ZPRES_HL(INDEX_ICLEAR_COL(:),1:KFLEV))) / ZCLEAR_COL_ll
    ZHT_CLEAR(1:KFLEV) = SUM_DD_R2_R1_ll(REAL(ZT_HL   (INDEX_ICLEAR_COL(:),1:KFLEV))) / ZCLEAR_COL_ll
    ! 
    ZALBP_CLEAR(:) = SUM_DD_R2_R1_ll(REAL(ZALBP(INDEX_ICLEAR_COL(:),:))) / ZCLEAR_COL_ll
    ZALBD_CLEAR(:) = SUM_DD_R2_R1_ll(REAL(ZALBD(INDEX_ICLEAR_COL(:),:))) / ZCLEAR_COL_ll
    ! 
    ZEMIS_CLEAR = SUM_DD_R1_ll(REAL(ZEMIS(INDEX_ICLEAR_COL(:),1))) / ZCLEAR_COL_ll
    ZEMIW_CLEAR = SUM_DD_R1_ll(REAL(ZEMIW(INDEX_ICLEAR_COL(:),1))) / ZCLEAR_COL_ll
    ZRMU0_CLEAR = SUM_DD_R1_ll(REAL(ZRMU0(INDEX_ICLEAR_COL(:))))   / ZCLEAR_COL_ll
    ZTS_CLEAR   = SUM_DD_R1_ll(REAL(ZTS(INDEX_ICLEAR_COL(:))))     / ZCLEAR_COL_ll
    ZLSM_CLEAR  = SUM_DD_R1_ll(REAL(ZLSM(INDEX_ICLEAR_COL(:))))    / ZCLEAR_COL_ll
    ZLAT_CLEAR  = SUM_DD_R1_ll(REAL(ZLAT(INDEX_ICLEAR_COL(:))))    / ZCLEAR_COL_ll
    ZLON_CLEAR  = SUM_DD_R1_ll(REAL(ZLON(INDEX_ICLEAR_COL(:))))    / ZCLEAR_COL_ll 
!
  ELSE ! no clear columns -> the first column is chosen, without physical meaning: it will not be
    ! unpacked after the call to the radiation ecmwf routine
    ZT_CLEAR(:)  = ZTAVE(1,:)
    ZP_CLEAR(:)  = ZPAVE(1,:)
    ZQV_CLEAR(:) = ZQVAVE(1,:)
    ZOZ_CLEAR(:) = ZO3AVE(1,:)
    ZDP_CLEAR(:) = ZDPRES(1,:)
    ZAER_CLEAR(:,:) = ZAER(1,:,:)
    IF(CAOP=='EXPL')THEN
       ZPIZA_EQ_CLEAR(:,:)=ZPIZA_EQ(1,:,:)
       ZCGA_EQ_CLEAR(:,:)=ZCGA_EQ(1,:,:)
       ZTAUREL_EQ_CLEAR(:,:)=ZTAUREL_EQ(1,:,:)
    ENDIF
!
    ZHP_CLEAR(1:KFLEV)  = ZPRES_HL(1,1:KFLEV)
    ZHT_CLEAR(1:KFLEV)  = ZT_HL(1,1:KFLEV)
    ZALBP_CLEAR(:) = ZALBP(1,:)
    ZALBD_CLEAR(:) = ZALBD(1,:)
!
    ZEMIS_CLEAR  = ZEMIS(1,1)
    ZEMIW_CLEAR  = ZEMIW(1,1) 
    ZRMU0_CLEAR  = ZRMU0(1)
    ZTS_CLEAR    = ZTS(1) 
    ZLSM_CLEAR   = ZLSM(1) 
    ZLAT_CLEAR   = ZLAT(1)
    ZLON_CLEAR   = ZLON(1)
  END IF
  !
  GCLOUD(:,:) = .NOT.GCLEAR(:,:) ! .true. where the column is cloudy
  GCLOUDT(:,:)=TRANSPOSE(GCLOUD(:,:))
  ICLOUD = ICLOUD_COL*KFLEV ! total number of voxels in cloudy columns
  ALLOCATE(ZWORK1(ICLOUD))
  ALLOCATE(ZWORK2(ICLOUD+KFLEV)) !  allocation for the KFLEV levels of 
                                 !  the ICLOUD cloudy columns
                                 !  and of the KFLEV levels of the clear sky one
  !
  ! temperature profiles
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZTAVE(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= ZT_CLEAR(1:)      !   and the single clear_sky one
  DEALLOCATE(ZTAVE)
  ALLOCATE(ZTAVE(ICLOUD_COL+1,KFLEV))
  ZTAVE(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  !
  ! vapor mixing ratio profiles
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZQVAVE(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= ZQV_CLEAR(1:)      !   and the single clear_sky one
  DEALLOCATE(ZQVAVE)
  ALLOCATE(ZQVAVE(ICLOUD_COL+1,KFLEV))
  ZQVAVE(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  !
  ! mesh size 
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZDZ(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= 0.0              !   and the single clear_sky one
  DEALLOCATE(ZDZ)
  ALLOCATE(ZDZ(ICLOUD_COL+1,KFLEV))
  ZDZ(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  !
  !
  ! liquid water mixing ratio profiles
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZQLAVE(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= 0.0              !   and the single clear_sky one
  DEALLOCATE(ZQLAVE)
  ALLOCATE(ZQLAVE(ICLOUD_COL+1,KFLEV))
  ZQLAVE(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  !
  !rain 
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZQRAVE(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= 0.0              !   and the single clear_sky one
  DEALLOCATE(ZQRAVE)
  ALLOCATE(ZQRAVE(ICLOUD_COL+1,KFLEV))
  ZQRAVE(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  ! 
  ! ice water mixing ratio profiles
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZQIAVE(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= 0.0      !   and the single clear_sky one
  DEALLOCATE(ZQIAVE)
  ALLOCATE(ZQIAVE(ICLOUD_COL+1,KFLEV))
  ZQIAVE(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  !
  !
  ! liquid water mixing ratio profiles
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZQLWC(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= 0.0              !   and the single clear_sky one
  DEALLOCATE(ZQLWC)
  ALLOCATE(ZQLWC(ICLOUD_COL+1,KFLEV))
  ZQLWC(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  !
  !rain 
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZQRWC(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= 0.0              !   and the single clear_sky one
  DEALLOCATE(ZQRWC)
  ALLOCATE(ZQRWC(ICLOUD_COL+1,KFLEV))
  ZQRWC(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  ! 
  ! ice water mixing ratio profiles
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZQIWC(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= 0.0      !   and the single clear_sky one
  DEALLOCATE(ZQIWC)
  ALLOCATE(ZQIWC(ICLOUD_COL+1,KFLEV))
  ZQIWC(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  !
  !
  ! cloud fraction profiles
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZCFAVE(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= 0.0      !   and the single clear_sky one
  DEALLOCATE(ZCFAVE)
  ALLOCATE(ZCFAVE(ICLOUD_COL+1,KFLEV))
  ZCFAVE(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  !
  ! C2R2 water particle concentration
  !
  IF ( SIZE(ZCCT_C2R2) > 0 )  THEN
    ZWORK1(:) = PACK( TRANSPOSE(ZCCT_C2R2(:,:)),MASK=GCLOUDT(:,:) )
    ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
    ZWORK2(ICLOUD+1:)= 0.0      !   and the single clear_sky one
    DEALLOCATE(ZCCT_C2R2)
    ALLOCATE(ZCCT_C2R2(ICLOUD_COL+1,KFLEV))
    ZCCT_C2R2 (:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  ENDIF
  IF ( SIZE (ZCRT_C2R2) > 0 )  THEN
    ZWORK1(:) = PACK( TRANSPOSE(ZCRT_C2R2(:,:)),MASK=GCLOUDT(:,:) )
    ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
    ZWORK2(ICLOUD+1:)= 0.0      !   and the single clear_sky one
    DEALLOCATE(ZCRT_C2R2)
    ALLOCATE(ZCRT_C2R2(ICLOUD_COL+1,KFLEV))
    ZCRT_C2R2 (:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  ENDIF
  IF ( SIZE (ZCIT_C1R3) > 0)  THEN
    ZWORK1(:) = PACK( TRANSPOSE(ZCIT_C1R3(:,:)),MASK=GCLOUDT(:,:) )
    ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
    ZWORK2(ICLOUD+1:)= 0.0      !   and the single clear_sky one
    DEALLOCATE(ZCIT_C1R3)
    ALLOCATE(ZCIT_C1R3(ICLOUD_COL+1,KFLEV))
    ZCIT_C1R3 (:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  ENDIF
  !
  ! LIMA water particle concentration
  !
  IF( CCLOUD == 'LIMA' ) THEN
     ZWORK1(:) = PACK( TRANSPOSE(ZCCT_LIMA(:,:)),MASK=GCLOUDT(:,:) )
     ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
     ZWORK2(ICLOUD+1:)= 0.0      !   and the single clear_sky one
     DEALLOCATE(ZCCT_LIMA)
     ALLOCATE(ZCCT_LIMA(ICLOUD_COL+1,KFLEV))
     ZCCT_LIMA (:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
!
     ZWORK1(:) = PACK( TRANSPOSE(ZCRT_LIMA(:,:)),MASK=GCLOUDT(:,:) )
     ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
     ZWORK2(ICLOUD+1:)= 0.0      !   and the single clear_sky one
     DEALLOCATE(ZCRT_LIMA)
     ALLOCATE(ZCRT_LIMA(ICLOUD_COL+1,KFLEV))
     ZCRT_LIMA (:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
!
     ZWORK1(:) = PACK( TRANSPOSE(ZCIT_LIMA(:,:)),MASK=GCLOUDT(:,:) )
     ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
     ZWORK2(ICLOUD+1:)= 0.0      !   and the single clear_sky one
     DEALLOCATE(ZCIT_LIMA)
     ALLOCATE(ZCIT_LIMA(ICLOUD_COL+1,KFLEV))
     ZCIT_LIMA (:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  ENDIF
  !
  ! ozone content profiles
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZO3AVE(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= ZOZ_CLEAR(1:)    !   and the single clear_sky one
  DEALLOCATE(ZO3AVE)
  ALLOCATE(ZO3AVE(ICLOUD_COL+1,KFLEV))
  ZO3AVE(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZPAVE(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= ZP_CLEAR(1:)     !   and the single clear_sky one
  DEALLOCATE(ZPAVE)
  ALLOCATE(ZPAVE(ICLOUD_COL+1,KFLEV))
  ZPAVE(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  !
  !pressure thickness
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZDPRES(:,:)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD) ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= ZDP_CLEAR(1:)    !   and the single clear_sky one
  DEALLOCATE(ZDPRES)
  ALLOCATE(ZDPRES(ICLOUD_COL+1,KFLEV))
  ZDPRES(:,:) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  !
  !aerosols
  !
  ALLOCATE(ZWORK1AER(ICLOUD,KAER))
  ALLOCATE(ZWORK2AER(ICLOUD+KFLEV,KAER))
  DO JK=1,KAER
    ZWORK1AER(:,JK) = PACK( TRANSPOSE(ZAER(:,:,JK)),MASK=GCLOUDT(:,:) )
    ZWORK2AER(1:ICLOUD,JK)=ZWORK1AER(:,JK)
    ZWORK2AER(ICLOUD+1:,JK)=ZAER_CLEAR(:,JK)
  END DO
  DEALLOCATE(ZAER)
  ALLOCATE(ZAER(ICLOUD_COL+1,KFLEV,KAER))
  DO JK=1,KAER
    ZAER(:,:,JK) = TRANSPOSE( RESHAPE( ZWORK2AER(:,JK),(/KFLEV,ICLOUD_COL+1/) ) )
  END DO
  DEALLOCATE (ZWORK1AER)
  DEALLOCATE (ZWORK2AER)
  !
  IF(CAOP=='EXPL')THEN
     ALLOCATE(ZWORK1AER(ICLOUD,KSWB_OLD))        !New vector with value for all cld. points
     ALLOCATE(ZWORK2AER(ICLOUD+KFLEV,KSWB_OLD))  !New vector with value for all cld.points + 1 clr column
     !Single scattering albedo
     DO WVL_IDX=1,KSWB_OLD
        ZWORK1AER(:,WVL_IDX) = PACK( TRANSPOSE(ZPIZA_EQ(:,:,WVL_IDX)),MASK=GCLOUDT(:,:) )
        ZWORK2AER(1:ICLOUD,WVL_IDX) = ZWORK1AER(:,WVL_IDX)
        ZWORK2AER(ICLOUD+1:,WVL_IDX) = ZPIZA_EQ_CLEAR(:,WVL_IDX)
     ENDDO
     DEALLOCATE(ZPIZA_EQ)
     ALLOCATE(ZPIZA_EQ(ICLOUD_COL+1,KFLEV,KSWB_OLD))
     DO WVL_IDX=1,KSWB_OLD
        ZPIZA_EQ(:,:,WVL_IDX) = TRANSPOSE( RESHAPE( ZWORK2AER(:,WVL_IDX),(/KFLEV,ICLOUD_COL+1/) ) )
     ENDDO
     !Assymetry factor
     DO WVL_IDX=1,KSWB_OLD
        ZWORK1AER(:,WVL_IDX) = PACK(TRANSPOSE(ZCGA_EQ(:,:,WVL_IDX)), MASK=GCLOUDT(:,:))
        ZWORK2AER(1:ICLOUD,WVL_IDX) = ZWORK1AER(:,WVL_IDX)
        ZWORK2AER(ICLOUD+1:,WVL_IDX) = ZCGA_EQ_CLEAR(:,WVL_IDX)
     ENDDO
     DEALLOCATE(ZCGA_EQ)
     ALLOCATE(ZCGA_EQ(ICLOUD_COL+1,KFLEV,KSWB_OLD))
     DO WVL_IDX=1,KSWB_OLD
        ZCGA_EQ(:,:,WVL_IDX) = TRANSPOSE(RESHAPE(ZWORK2AER(:,WVL_IDX),(/KFLEV,ICLOUD_COL+1/)))
     ENDDO
     !Relative wavelength-distributed optical depth
     DO WVL_IDX=1,KSWB_OLD
        ZWORK1AER(:,WVL_IDX) =  PACK(TRANSPOSE(ZTAUREL_EQ(:,:,WVL_IDX)), MASK=GCLOUDT(:,:))
        ZWORK2AER(1:ICLOUD,WVL_IDX) = ZWORK1AER(:,WVL_IDX)
        ZWORK2AER(ICLOUD+1:,WVL_IDX) = ZTAUREL_EQ_CLEAR(:,WVL_IDX)
     ENDDO
     DEALLOCATE(ZTAUREL_EQ)
     ALLOCATE(ZTAUREL_EQ(ICLOUD_COL+1,KFLEV,KSWB_OLD))
     DO WVL_IDX=1,KSWB_OLD
        ZTAUREL_EQ(:,:,WVL_IDX) = TRANSPOSE(RESHAPE(ZWORK2AER(:,WVL_IDX),(/KFLEV,ICLOUD_COL+1/)))
     ENDDO
     DEALLOCATE(ZWORK1AER)
     DEALLOCATE(ZWORK2AER)
  ELSE
     DEALLOCATE(ZPIZA_EQ)
     ALLOCATE(ZPIZA_EQ(ICLOUD_COL+1,KFLEV,KSWB_OLD))
     DEALLOCATE(ZCGA_EQ)
     ALLOCATE(ZCGA_EQ(ICLOUD_COL+1,KFLEV,KSWB_OLD))
     DEALLOCATE(ZTAUREL_EQ)
     ALLOCATE(ZTAUREL_EQ(ICLOUD_COL+1,KFLEV,KSWB_OLD))
  ENDIF !Check on LDUST
  
  ! half-level variables
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZPRES_HL(:,1:KFLEV)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD)  ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= ZHP_CLEAR(1:)     !   and the single clear_sky one
  DEALLOCATE(ZPRES_HL)
  ALLOCATE(ZPRES_HL(ICLOUD_COL+1,KFLEV+1))
  ZPRES_HL(:,1:KFLEV) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  ZPRES_HL(:,KFLEV+1) = PSTATM(IKSTAE,2)*100.0
  !
  ZWORK1(:) = PACK( TRANSPOSE(ZT_HL(:,1:KFLEV)),MASK=GCLOUDT(:,:) )
  ZWORK2(1:ICLOUD) = ZWORK1(1:ICLOUD)  ! fills the ICLOUD_COL cloudy columns
  ZWORK2(ICLOUD+1:)= ZHT_CLEAR(1:)     !   and the single clear_sky one
  DEALLOCATE(ZT_HL)
  ALLOCATE(ZT_HL(ICLOUD_COL+1,KFLEV+1))
  ZT_HL(:,1:KFLEV) = TRANSPOSE( RESHAPE( ZWORK2(:),(/KFLEV,ICLOUD_COL+1/) ) )
  ZT_HL(:,KFLEV+1) = PSTATM(IKSTAE,3)
  !
  ! surface fields
  !
  ALLOCATE(ZWORK3(ICLOUD_COL))
  ALLOCATE(ZWORK4(ICLOUD_COL,KSWB_MNH))
  ALLOCATE(ZWORK(KDLON))
  DO JALBS=1,KSWB_MNH
    ZWORK(:)  = ZALBP(:,JALBS)
    ZWORK3(:) = PACK( ZWORK(:),MASK=.NOT.GCLEAR_2D(:) )
    ZWORK4(:,JALBS) = ZWORK3(:)
  END DO
  DEALLOCATE(ZALBP)
  ALLOCATE(ZALBP(ICLOUD_COL+1,KSWB_MNH))
  ZALBP(1:ICLOUD_COL,:) = ZWORK4(1:ICLOUD_COL,:)
  ZALBP(ICLOUD_COL+1,:) = ZALBP_CLEAR(:)
  !
  DO JALBS=1,KSWB_MNH
    ZWORK(:)  = ZALBD(:,JALBS)
    ZWORK3(:) = PACK( ZWORK(:),MASK=.NOT.GCLEAR_2D(:) )
    ZWORK4(:,JALBS) = ZWORK3(:)
  END DO
  DEALLOCATE(ZALBD)
  ALLOCATE(ZALBD(ICLOUD_COL+1,KSWB_MNH))
  ZALBD(1:ICLOUD_COL,:) = ZWORK4(1:ICLOUD_COL,:)  
  ZALBD(ICLOUD_COL+1,:) = ZALBD_CLEAR(:)  
  !
  DEALLOCATE(ZWORK4)
  !
  ZWORK3(:) = PACK( ZEMIS(:,1),MASK=.NOT.GCLEAR_2D(:) )
  DEALLOCATE(ZEMIS)
  ALLOCATE(ZEMIS(ICLOUD_COL+1,1))
  ZEMIS(1:ICLOUD_COL,1) = ZWORK3(1:ICLOUD_COL)
  ZEMIS(ICLOUD_COL+1,1) = ZEMIS_CLEAR
  !
  !
  ZWORK3(:) = PACK( ZEMIW(:,1),MASK=.NOT.GCLEAR_2D(:) )
  DEALLOCATE(ZEMIW)
  ALLOCATE(ZEMIW(ICLOUD_COL+1,1))
  ZEMIW(1:ICLOUD_COL,1) = ZWORK3(1:ICLOUD_COL)
  ZEMIW(ICLOUD_COL+1,1) = ZEMIW_CLEAR
  ! 
  !
  ZWORK3(:) = PACK( ZRMU0(:),MASK=.NOT.GCLEAR_2D(:) )
  DEALLOCATE(ZRMU0)
  ALLOCATE(ZRMU0(ICLOUD_COL+1))
  ZRMU0(1:ICLOUD_COL) = ZWORK3(1:ICLOUD_COL)
  ZRMU0(ICLOUD_COL+1) = ZRMU0_CLEAR
  !
  ZWORK3(:) = PACK( ZLSM(:),MASK=.NOT.GCLEAR_2D(:) )
  DEALLOCATE(ZLSM)
  ALLOCATE(ZLSM(ICLOUD_COL+1))
  ZLSM(1:ICLOUD_COL) = ZWORK3(1:ICLOUD_COL)
  ZLSM (ICLOUD_COL+1)= ZLSM_CLEAR
  ! 
  ZWORK3(:) = PACK( ZLAT(:),MASK=.NOT.GCLEAR_2D(:) )
  DEALLOCATE(ZLAT)
  ALLOCATE(ZLAT(ICLOUD_COL+1))
  ZLAT(1:ICLOUD_COL) = ZWORK3(1:ICLOUD_COL)
  ZLAT (ICLOUD_COL+1)= ZLAT_CLEAR
  ! 
  ZWORK3(:) = PACK( ZLON(:),MASK=.NOT.GCLEAR_2D(:) )
  DEALLOCATE(ZLON)
  ALLOCATE(ZLON(ICLOUD_COL+1))
  ZLON(1:ICLOUD_COL) = ZWORK3(1:ICLOUD_COL)
  ZLON (ICLOUD_COL+1)= ZLON_CLEAR
  !   
  ZWORK3(:) = PACK( ZTS(:),MASK=.NOT.GCLEAR_2D(:) )
  DEALLOCATE(ZTS)
  ALLOCATE(ZTS(ICLOUD_COL+1))
  ZTS(1:ICLOUD_COL) = ZWORK3(1:ICLOUD_COL)
  ZTS(ICLOUD_COL+1) = ZTS_CLEAR
  !
  DEALLOCATE(ZWORK1)
  DEALLOCATE(ZWORK2)
  DEALLOCATE(ZWORK3)
  DEALLOCATE(ZWORK)
  !  
  IDIM = ICLOUD_COL +1 ! Number of columns where RT is computed 
!
ELSE 
  !
  !*       5.3   RADIATION COMPUTATIONS FOR THE FULL COLUMN NUMBER (KDLON)
  !
  IDIM = KDLON
END IF
!
! initialisation of cloud trace for the next radiation time step
! (if unchanged columns are not recomputed)
WHERE ( ZCLOUD(:) <= 0.0 )
  ICLEAR_2D_TM1(:) = 1
ELSEWHERE
  ICLEAR_2D_TM1(:) = 0
END WHERE
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    KCLEARCOL_TM1(JI,JJ) = ICLEAR_2D_TM1(IIJ) ! output to be saved for next time step
  END DO
END DO
! 
!
!*       5.4  VERTICAL grid modification(up-down) for compatibility with ECMWF 
!             radiation vertical grid. ALLOCATION of the outputs.  
!            
!             
ALLOCATE (ZWORK_GRID(SIZE(ZPRES_HL,1),KFLEV+1))
!
!half level pressure
ZWORK_GRID(:,:)=ZPRES_HL(:,:)
DO JKRAD=1, KFLEV+1
  JK1=(KFLEV+1)+1-JKRAD
  ZPRES_HL(:,JKRAD) = ZWORK_GRID(:,JK1)
END DO
!
!half level temperature
ZWORK_GRID(:,:)=ZT_HL(:,:)
DO  JKRAD=1, KFLEV+1
  JK1=(KFLEV+1)+1-JKRAD
  ZT_HL(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO
!
DEALLOCATE(ZWORK_GRID)
!
!mean layer variables
!-------------------------------------
ALLOCATE(ZWORK_GRID(SIZE(ZTAVE,1),KFLEV))
!
!mean layer temperature
ZWORK_GRID(:,:)=ZTAVE(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZTAVE(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO
!
!mean layer pressure
ZWORK_GRID(:,:)=ZPAVE(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZPAVE(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO
!
!mean layer pressure thickness
ZWORK_GRID(:,:)=ZDPRES(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZDPRES(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO
!
!mesh size
ZWORK_GRID(:,:)=ZDZ(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZDZ(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO

!mean layer cloud fraction
ZWORK_GRID(:,:)=ZCFAVE(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZCFAVE(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO
!
!mean layer water vapor mixing ratio
ZWORK_GRID(:,:)=ZQVAVE(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZQVAVE(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO
!
!ice
ZWORK_GRID(:,:)=ZQIAVE(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZQIAVE(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO
!
!liquid water
ZWORK_GRID(:,:)=ZQLAVE(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZQLAVE(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO


!rain water
ZWORK_GRID(:,:)=ZQRAVE(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZQRAVE(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO
!
!ice water content
ZWORK_GRID(:,:)=ZQIWC(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZQIWC(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO
!
!liquid water content
ZWORK_GRID(:,:)=ZQLWC(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZQLWC(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO


!rain water content
ZWORK_GRID(:,:)=ZQRWC(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZQRWC(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO


!C2R2 water particle concentration
!
IF (SIZE(ZCCT_C2R2) > 0) THEN
  ZWORK_GRID(:,:)=ZCCT_C2R2(:,:)
  DO JKRAD=1, KFLEV
    JK1=KFLEV+1-JKRAD
    ZCCT_C2R2(:,JKRAD)=ZWORK_GRID(:,JK1)
  END DO
END IF
IF (SIZE(ZCRT_C2R2) > 0) THEN
  ZWORK_GRID(:,:)=ZCRT_C2R2(:,:)
  DO JKRAD=1, KFLEV
    JK1=KFLEV+1-JKRAD
    ZCRT_C2R2(:,JKRAD)=ZWORK_GRID(:,JK1)
  END DO
END IF
IF (SIZE(ZCIT_C1R3) > 0) THEN
  ZWORK_GRID(:,:)=ZCIT_C1R3(:,:)
  DO JKRAD=1, KFLEV
    JK1=KFLEV+1-JKRAD
    ZCIT_C1R3(:,JKRAD)=ZWORK_GRID(:,JK1)
  END DO
END IF
!
!LIMA water particle concentration
!
IF( CCLOUD == 'LIMA' ) THEN
   ZWORK_GRID(:,:)=ZCCT_LIMA(:,:)
   DO JKRAD=1, KFLEV
      JK1=KFLEV+1-JKRAD
      ZCCT_LIMA(:,JKRAD)=ZWORK_GRID(:,JK1)
   END DO
!
   ZWORK_GRID(:,:)=ZCRT_LIMA(:,:)
   DO JKRAD=1, KFLEV
      JK1=KFLEV+1-JKRAD
      ZCRT_LIMA(:,JKRAD)=ZWORK_GRID(:,JK1)
   END DO
!
   ZWORK_GRID(:,:)=ZCIT_LIMA(:,:)
   DO JKRAD=1, KFLEV
      JK1=KFLEV+1-JKRAD
      ZCIT_LIMA(:,JKRAD)=ZWORK_GRID(:,JK1)
   END DO
END IF
!
!ozone content 
ZWORK_GRID(:,:)=ZO3AVE(:,:)
DO JKRAD=1, KFLEV
  JK1=KFLEV+1-JKRAD
  ZO3AVE(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO
!
!aerosol optical depth 
DO JI=1,KAER
  ZWORK_GRID(:,:)=ZAER(:,:,JI)
  DO JKRAD=1, KFLEV
    JK1=KFLEV+1-JKRAD
    ZAER(:,JKRAD,JI)=ZWORK_GRID(:,JK1)
  END DO
END DO
IF (CAOP=='EXPL') THEN
!TURN MORE FIELDS UPSIDE DOWN...
!Dust single scattering albedo
DO JI=1,KSWB_OLD
   ZWORK_GRID(:,:)=ZPIZA_EQ(:,:,JI)
   DO JKRAD=1,KFLEV
      JK1=KFLEV+1-JKRAD
      ZPIZA_EQ(:,JKRAD,JI)=ZWORK_GRID(:,JK1)
   ENDDO
ENDDO
!Dust asymmetry factor
DO JI=1,KSWB_OLD
   ZWORK_GRID(:,:)=ZCGA_EQ(:,:,JI)
   DO JKRAD=1,KFLEV
      JK1=KFLEV+1-JKRAD
      ZCGA_EQ(:,JKRAD,JI)=ZWORK_GRID(:,JK1)
   ENDDO
ENDDO
DO JI=1,KSWB_OLD
   ZWORK_GRID(:,:)=ZTAUREL_EQ(:,:,JI)
   DO JKRAD=1,KFLEV
      JK1=KFLEV+1-JKRAD
      ZTAUREL_EQ(:,JKRAD,JI)=ZWORK_GRID(:,JK1)
   ENDDO
ENDDO 

END IF

!
DEALLOCATE(ZWORK_GRID)
!
!mean layer saturation specific humidity
!
ALLOCATE(ZQSAVE(SIZE(ZTAVE,1),SIZE(ZTAVE,2)))
!
WHERE (ZTAVE(:,:) > XTT)
  ZQSAVE(:,:) = QSAT(ZTAVE, ZPAVE)
ELSEWHERE
  ZQSAVE(:,:) = QSATI(ZTAVE, ZPAVE)
END WHERE
!
! allocations for the radiation code outputs
!
ALLOCATE(ZDTLW(IDIM,KFLEV))
ALLOCATE(ZDTSW(IDIM,KFLEV))
ALLOCATE(ZFLUX_TOP_GND_IRVISNIR(IDIM,KFLUX))
ALLOCATE(ZSFSWDIR(IDIM,ISWB))
ALLOCATE(ZSFSWDIF(IDIM,ISWB))
ALLOCATE(ZDTLW_CS(IDIM,KFLEV))
ALLOCATE(ZDTSW_CS(IDIM,KFLEV))
ALLOCATE(ZFLUX_TOP_GND_IRVISNIR_CS(IDIM,KFLUX))
!
!
ALLOCATE(ZFLUX_LW(IDIM,2,KFLEV+1))
ALLOCATE(ZFLUX_SW_DOWN(IDIM,KFLEV+1))
ALLOCATE(ZFLUX_SW_UP(IDIM,KFLEV+1))
ALLOCATE(ZRADLP(IDIM,KFLEV))
IF( KRAD_DIAG >= 1) THEN
  ALLOCATE(ZNFLW(IDIM,KFLEV+1))
  ALLOCATE(ZNFSW(IDIM,KFLEV+1))
ELSE
  ALLOCATE(ZNFLW(0,0))
  ALLOCATE(ZNFSW(0,0))
END IF
! 
IF( KRAD_DIAG >= 2) THEN
  ALLOCATE(ZFLUX_SW_DOWN_CS(IDIM,KFLEV+1))
  ALLOCATE(ZFLUX_SW_UP_CS(IDIM,KFLEV+1))
  ALLOCATE(ZFLUX_LW_CS(IDIM,2,KFLEV+1))
  ALLOCATE(ZNFLW_CS(IDIM,KFLEV+1))
  ALLOCATE(ZNFSW_CS(IDIM,KFLEV+1))
ELSE
  ALLOCATE(ZFLUX_SW_DOWN_CS(0,0))
  ALLOCATE(ZFLUX_SW_UP_CS(0,0))
  ALLOCATE(ZFLUX_LW_CS(0,0,0))
  ALLOCATE(ZNFSW_CS(0,0))
  ALLOCATE(ZNFLW_CS(0,0))
END IF
!
IF( KRAD_DIAG >= 3) THEN
  ALLOCATE(ZPLAN_ALB_VIS(IDIM))
  ALLOCATE(ZPLAN_ALB_NIR(IDIM))
  ALLOCATE(ZPLAN_TRA_VIS(IDIM))
  ALLOCATE(ZPLAN_TRA_NIR(IDIM))
  ALLOCATE(ZPLAN_ABS_VIS(IDIM))
  ALLOCATE(ZPLAN_ABS_NIR(IDIM))
ELSE
  ALLOCATE(ZPLAN_ALB_VIS(0))
  ALLOCATE(ZPLAN_ALB_NIR(0))
  ALLOCATE(ZPLAN_TRA_VIS(0))
  ALLOCATE(ZPLAN_TRA_NIR(0))
  ALLOCATE(ZPLAN_ABS_VIS(0))
  ALLOCATE(ZPLAN_ABS_NIR(0))
END IF
!
IF( KRAD_DIAG >= 4) THEN
  ALLOCATE(ZEFCL_RRTM(IDIM,KFLEV))
  ALLOCATE(ZCLSW_TOTAL(IDIM,KFLEV))
  ALLOCATE(ZTAU_TOTAL(IDIM,KSWB_OLD,KFLEV))
  ALLOCATE(ZOMEGA_TOTAL(IDIM,KSWB_OLD,KFLEV))
  ALLOCATE(ZCG_TOTAL(IDIM,KSWB_OLD,KFLEV))
  ALLOCATE(ZEFCL_LWD(IDIM,KFLEV))
  ALLOCATE(ZEFCL_LWU(IDIM,KFLEV))
  ALLOCATE(ZFLWP(IDIM,KFLEV))
  ALLOCATE(ZFIWP(IDIM,KFLEV))  
  ALLOCATE(ZRADIP(IDIM,KFLEV)) 
ELSE
  ALLOCATE(ZEFCL_RRTM(0,0))
  ALLOCATE(ZCLSW_TOTAL(0,0))
  ALLOCATE(ZTAU_TOTAL(0,0,0))
  ALLOCATE(ZOMEGA_TOTAL(0,0,0))
  ALLOCATE(ZCG_TOTAL(0,0,0))
  ALLOCATE(ZEFCL_LWD(0,0))
  ALLOCATE(ZEFCL_LWU(0,0))
  ALLOCATE(ZFLWP(0,0))
  ALLOCATE(ZFIWP(0,0))
  ALLOCATE(ZRADIP(0,0))
END IF
!
!*       5.6   CALLS THE ECMWF_RADIATION ROUTINES
!
! mixing ratio -> specific humidity conversion (for ECMWF routine)
! mixing ratio = mv/md ; specific humidity = mv/(mv+md)

ZQVAVE(:,:) = ZQVAVE(:,:) / (1.+ZQVAVE(:,:)) ! Because 
! ZAER = 1e-5*ZAER
! ZO3AVE = 1e-5*ZO3AVE!
IF( IDIM <= KRAD_COLNBR ) THEN 
!
! there is less than KRAD_COLNBR columns to be considered therefore
! no split of the arrays is performed
! Note that radiation scheme only takes scalar emissivities so only fist value of the spectral emissivity is taken
 ALLOCATE(ZTAVE_RAD(SIZE(ZTAVE,1),SIZE(ZTAVE,2)))
 ALLOCATE(ZPAVE_RAD(SIZE(ZPAVE,1),SIZE(ZPAVE,2)))
 ZTAVE_RAD = ZTAVE
 ZPAVE_RAD = ZPAVE
 IF (CCLOUD == 'LIMA') THEN
  IF (CRAD == "ECMW") THEN
    CALL ECMWF_RADIATION_VERS2  ( IDIM ,KFLEV, KRAD_DIAG, KAER,     &      
        ZDZ,HEFRADL,HEFRADI,HOPWSW, HOPISW, HOPWLW, HOPILW,PFUDG,      &
        ZRII0, ZAER , ZALBD, ZALBP, ZPRES_HL, ZPAVE_RAD,               &
        PCCO2, ZCFAVE, ZDPRES, ZEMIS(:,1), ZEMIW(:,1), ZLSM, ZRMU0,          &
        ZO3AVE , ZQVAVE, ZQIAVE ,ZQIWC,ZQLAVE,ZQLWC, ZQSAVE, ZQRAVE,  ZQRWC,  &
        ZT_HL,ZTAVE_RAD, ZTS, ZCCT_LIMA, ZCRT_LIMA, ZCIT_LIMA,         &
        ZNFLW_CS, ZNFLW, ZNFSW_CS,ZNFSW,                           &
        ZDTLW, ZDTSW, ZFLUX_TOP_GND_IRVISNIR,                      &
        ZSFSWDIR, ZSFSWDIF,                                        &
        ZFLUX_SW_DOWN, ZFLUX_SW_UP, ZFLUX_LW ,                     &
        ZDTLW_CS, ZDTSW_CS, ZFLUX_TOP_GND_IRVISNIR_CS,             &
        ZFLUX_SW_DOWN_CS, ZFLUX_SW_UP_CS, ZFLUX_LW_CS,             &           
        ZPLAN_ALB_VIS,ZPLAN_ALB_NIR, ZPLAN_TRA_VIS, ZPLAN_TRA_NIR, &
        ZPLAN_ABS_VIS, ZPLAN_ABS_NIR,ZEFCL_LWD, ZEFCL_LWU,         &
        ZFLWP, ZFIWP,ZRADLP, ZRADIP,ZEFCL_RRTM,  ZCLSW_TOTAL,  ZTAU_TOTAL,  &
        ZOMEGA_TOTAL,ZCG_TOTAL,                                    &
        GAOP, ZPIZA_EQ,ZCGA_EQ,ZTAUREL_EQ                       )

        
  ELSE IF (CRAD == "ECRA") THEN  
    CALL ECRAD_INTERFACE  ( IDIM ,KFLEV, KRAD_DIAG, KAER,     &      
        ZDZ,HEFRADL,HEFRADI,HOPWSW, HOPISW, HOPWLW, HOPILW,PFUDG,      &
        ZRII0, ZAER , ZALBD, ZALBP, ZPRES_HL, ZPAVE_RAD,               &
        PCCO2, ZCFAVE, ZDPRES, ZEMIS(:,1), ZEMIW(:,1), ZLSM, ZRMU0,          &
        ZO3AVE , ZQVAVE, ZQIAVE ,ZQIWC,ZQLAVE,ZQLWC, ZQSAVE, ZQRAVE,  ZQRWC,  &
        ZT_HL,ZTAVE_RAD, ZTS, ZCCT_LIMA, ZCRT_LIMA, ZCIT_LIMA,         &
        ZNFLW, ZNFSW, ZNFLW_CS, ZNFSW_CS,                           &
        ZDTLW, ZDTSW, ZFLUX_TOP_GND_IRVISNIR,                      &
        ZSFSWDIR, ZSFSWDIF,                                        &
        ZFLUX_SW_DOWN, ZFLUX_SW_UP, ZFLUX_LW ,                     &
        ZDTLW_CS, ZDTSW_CS, ZFLUX_TOP_GND_IRVISNIR_CS,             &
        ZFLUX_SW_DOWN_CS, ZFLUX_SW_UP_CS, ZFLUX_LW_CS,             &           
        ZPLAN_ALB_VIS,ZPLAN_ALB_NIR, ZPLAN_TRA_VIS, ZPLAN_TRA_NIR, &
        ZPLAN_ABS_VIS, ZPLAN_ABS_NIR,ZEFCL_LWD, ZEFCL_LWU,         &
        ZFLWP, ZFIWP,ZRADLP, ZRADIP,ZEFCL_RRTM,  ZCLSW_TOTAL,  ZTAU_TOTAL,  &
        ZOMEGA_TOTAL,ZCG_TOTAL,                                    &
        GAOP, ZPIZA_EQ,ZCGA_EQ,ZTAUREL_EQ,ZLAT,ZLON                )      
  ENDIF      
 
 ELSE
   IF (CRAD == "ECMW") THEN 
     CALL ECMWF_RADIATION_VERS2  ( IDIM ,KFLEV, KRAD_DIAG, KAER,     &      
        ZDZ,HEFRADL,HEFRADI,HOPWSW, HOPISW, HOPWLW, HOPILW,PFUDG,      &
        ZRII0, ZAER , ZALBD, ZALBP, ZPRES_HL, ZPAVE_RAD,               &
        PCCO2, ZCFAVE, ZDPRES, ZEMIS(:,1), ZEMIW(:,1), ZLSM, ZRMU0,          &
        ZO3AVE , ZQVAVE, ZQIAVE ,ZQIWC,ZQLAVE,ZQLWC, ZQSAVE, ZQRAVE,  ZQRWC,  &
        ZT_HL,ZTAVE_RAD, ZTS, ZCCT_C2R2, ZCRT_C2R2, ZCIT_C1R3,         &
        ZNFLW_CS, ZNFLW, ZNFSW_CS,ZNFSW,                           &
        ZDTLW, ZDTSW, ZFLUX_TOP_GND_IRVISNIR,                      &
        ZSFSWDIR, ZSFSWDIF,                                        &
        ZFLUX_SW_DOWN, ZFLUX_SW_UP, ZFLUX_LW ,                     &
        ZDTLW_CS, ZDTSW_CS, ZFLUX_TOP_GND_IRVISNIR_CS,             &
        ZFLUX_SW_DOWN_CS, ZFLUX_SW_UP_CS, ZFLUX_LW_CS,             &           
        ZPLAN_ALB_VIS,ZPLAN_ALB_NIR, ZPLAN_TRA_VIS, ZPLAN_TRA_NIR, &
        ZPLAN_ABS_VIS, ZPLAN_ABS_NIR,ZEFCL_LWD, ZEFCL_LWU,         &
        ZFLWP, ZFIWP,ZRADLP, ZRADIP,ZEFCL_RRTM,  ZCLSW_TOTAL,  ZTAU_TOTAL,  &
        ZOMEGA_TOTAL,ZCG_TOTAL,                                    &
        GAOP, ZPIZA_EQ,ZCGA_EQ,ZTAUREL_EQ                       )
   
   ELSE IF (CRAD == "ECRA") THEN 
     CALL ECRAD_INTERFACE  ( IDIM ,KFLEV, KRAD_DIAG, KAER,                   &      
        ZDZ,HEFRADL,HEFRADI,HOPWSW, HOPISW, HOPWLW, HOPILW,PFUDG,      &
        ZRII0, ZAER , ZALBD, ZALBP, ZPRES_HL, ZPAVE_RAD,               &
        PCCO2, ZCFAVE, ZDPRES, ZEMIS(:,1), ZEMIW(:,1), ZLSM, ZRMU0,          &
        ZO3AVE , ZQVAVE, ZQIAVE ,ZQIWC,ZQLAVE,ZQLWC, ZQSAVE, ZQRAVE,  ZQRWC,  &
        ZT_HL,ZTAVE_RAD, ZTS, ZCCT_C2R2, ZCRT_C2R2, ZCIT_C1R3,         &
        ZNFLW, ZNFSW, ZNFLW_CS, ZNFSW_CS,                           &
        ZDTLW, ZDTSW, ZFLUX_TOP_GND_IRVISNIR,                      &
        ZSFSWDIR, ZSFSWDIF,                                        &
        ZFLUX_SW_DOWN, ZFLUX_SW_UP, ZFLUX_LW ,                     &
        ZDTLW_CS, ZDTSW_CS, ZFLUX_TOP_GND_IRVISNIR_CS,             &
        ZFLUX_SW_DOWN_CS, ZFLUX_SW_UP_CS, ZFLUX_LW_CS,             &           
        ZPLAN_ALB_VIS,ZPLAN_ALB_NIR, ZPLAN_TRA_VIS, ZPLAN_TRA_NIR, &
        ZPLAN_ABS_VIS, ZPLAN_ABS_NIR,ZEFCL_LWD, ZEFCL_LWU,         &
        ZFLWP, ZFIWP,ZRADLP, ZRADIP,ZEFCL_RRTM,  ZCLSW_TOTAL,  ZTAU_TOTAL,  &
        ZOMEGA_TOTAL,ZCG_TOTAL,                                    &
        GAOP, ZPIZA_EQ,ZCGA_EQ,ZTAUREL_EQ ,ZLAT,ZLON               )
  END IF      
   
   
 END IF
 DEALLOCATE(ZTAVE_RAD,ZPAVE_RAD)
!
ELSE
!
! the splitting of the arrays will be performed
!
  INUM_CALL = CEILING( REAL( IDIM ) / REAL( KRAD_COLNBR ) )
  IDIM_RESIDUE = IDIM
!
  DO JI_SPLIT = 1 , INUM_CALL
    IDIM_EFF = MIN( IDIM_RESIDUE,KRAD_COLNBR )
    !
    IF( JI_SPLIT == 1 .OR. JI_SPLIT == INUM_CALL ) THEN       
      ALLOCATE(  ZALBP_SPLIT(IDIM_EFF,KSWB_MNH))
      ALLOCATE(  ZALBD_SPLIT(IDIM_EFF,KSWB_MNH))  
      ALLOCATE(  ZEMIS_SPLIT(IDIM_EFF))
      ALLOCATE(  ZEMIW_SPLIT(IDIM_EFF))
      ALLOCATE(  ZRMU0_SPLIT(IDIM_EFF))
      ALLOCATE(  ZLAT_SPLIT(IDIM_EFF))
      ALLOCATE(  ZLON_SPLIT(IDIM_EFF))
      ALLOCATE(  ZCFAVE_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZO3AVE_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZT_HL_SPLIT(IDIM_EFF,KFLEV+1))
      ALLOCATE(  ZPRES_HL_SPLIT(IDIM_EFF,KFLEV+1))
      ALLOCATE(  ZDZ_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZQLAVE_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZQIAVE_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZQRAVE_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZQLWC_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZQIWC_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZQRWC_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZQVAVE_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZTAVE_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZPAVE_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZAER_SPLIT( IDIM_EFF,KFLEV,KAER))
      ALLOCATE(  ZPIZA_EQ_SPLIT(IDIM_EFF,KFLEV,KSWB_OLD))
      ALLOCATE(  ZCGA_EQ_SPLIT(IDIM_EFF,KFLEV,KSWB_OLD))
      ALLOCATE(  ZTAUREL_EQ_SPLIT(IDIM_EFF,KFLEV,KSWB_OLD))
      ALLOCATE(  ZDPRES_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZLSM_SPLIT(IDIM_EFF))
      ALLOCATE(  ZQSAVE_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZTS_SPLIT(IDIM_EFF))
      ! output pronostic       
      ALLOCATE(  ZDTLW_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZDTSW_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZFLUX_TOP_GND_IRVISNIR_SPLIT(IDIM_EFF,KFLUX))
      ALLOCATE(  ZSFSWDIR_SPLIT(IDIM_EFF,ISWB))
      ALLOCATE(  ZSFSWDIF_SPLIT(IDIM_EFF,ISWB))
      ALLOCATE(  ZDTLW_CS_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZDTSW_CS_SPLIT(IDIM_EFF,KFLEV))
      ALLOCATE(  ZFLUX_TOP_GND_IRVISNIR_CS_SPLIT(IDIM_EFF,KFLUX))
!
      ALLOCATE(  ZFLUX_LW_SPLIT(IDIM_EFF,2,KFLEV+1))
      ALLOCATE(  ZFLUX_SW_DOWN_SPLIT(IDIM_EFF,KFLEV+1))
      ALLOCATE(  ZFLUX_SW_UP_SPLIT(IDIM_EFF,KFLEV+1))
      ALLOCATE(  ZRADLP_SPLIT(IDIM_EFF,KFLEV))
      IF(KRAD_DIAG >=1) THEN
        ALLOCATE(  ZNFSW_SPLIT(IDIM_EFF,KFLEV+1))
        ALLOCATE(  ZNFLW_SPLIT(IDIM_EFF,KFLEV+1))
      ELSE
        ALLOCATE(  ZNFSW_SPLIT(0,0))
        ALLOCATE(  ZNFLW_SPLIT(0,0))
      END IF
!
      IF( KRAD_DIAG >= 2) THEN      
        ALLOCATE(  ZFLUX_SW_DOWN_CS_SPLIT(IDIM_EFF,KFLEV+1))
        ALLOCATE(  ZFLUX_SW_UP_CS_SPLIT(IDIM_EFF,KFLEV+1))
        ALLOCATE(  ZFLUX_LW_CS_SPLIT(IDIM_EFF,2,KFLEV+1))
        ALLOCATE(   ZNFSW_CS_SPLIT(IDIM_EFF,KFLEV+1))
        ALLOCATE(   ZNFLW_CS_SPLIT(IDIM_EFF,KFLEV+1))
      ELSE
        ALLOCATE(  ZFLUX_SW_DOWN_CS_SPLIT(0,0))
        ALLOCATE(  ZFLUX_SW_UP_CS_SPLIT(0,0))
        ALLOCATE(  ZFLUX_LW_CS_SPLIT(0,0,0))
        ALLOCATE(  ZNFSW_CS_SPLIT(0,0))
        ALLOCATE(  ZNFLW_CS_SPLIT(0,0))
      END IF
!
      IF( KRAD_DIAG >= 3) THEN
        ALLOCATE(  ZPLAN_ALB_VIS_SPLIT(IDIM_EFF))
        ALLOCATE(  ZPLAN_ALB_NIR_SPLIT(IDIM_EFF))
        ALLOCATE(  ZPLAN_TRA_VIS_SPLIT(IDIM_EFF))
        ALLOCATE(  ZPLAN_TRA_NIR_SPLIT(IDIM_EFF))
        ALLOCATE(  ZPLAN_ABS_VIS_SPLIT(IDIM_EFF))
        ALLOCATE(  ZPLAN_ABS_NIR_SPLIT(IDIM_EFF))
      ELSE
        ALLOCATE(  ZPLAN_ALB_VIS_SPLIT(0))
        ALLOCATE(  ZPLAN_ALB_NIR_SPLIT(0))
        ALLOCATE(  ZPLAN_TRA_VIS_SPLIT(0))
        ALLOCATE(  ZPLAN_TRA_NIR_SPLIT(0))
        ALLOCATE(  ZPLAN_ABS_VIS_SPLIT(0))
        ALLOCATE(  ZPLAN_ABS_NIR_SPLIT(0))
      END IF
!
      IF( KRAD_DIAG >= 4) THEN
        ALLOCATE(  ZEFCL_RRTM_SPLIT(IDIM_EFF,KFLEV))
        ALLOCATE(  ZCLSW_TOTAL_SPLIT(IDIM_EFF,KFLEV))
        ALLOCATE(  ZTAU_TOTAL_SPLIT(IDIM_EFF,KSWB_OLD,KFLEV))
        ALLOCATE(  ZOMEGA_TOTAL_SPLIT(IDIM_EFF,KSWB_OLD,KFLEV))
        ALLOCATE(  ZCG_TOTAL_SPLIT(IDIM_EFF,KSWB_OLD,KFLEV))
        ALLOCATE(  ZEFCL_LWD_SPLIT(IDIM_EFF,KFLEV))
        ALLOCATE(  ZEFCL_LWU_SPLIT(IDIM_EFF,KFLEV))
        ALLOCATE(  ZFLWP_SPLIT(IDIM_EFF,KFLEV))
        ALLOCATE(  ZFIWP_SPLIT(IDIM_EFF,KFLEV))
        ALLOCATE(  ZRADIP_SPLIT(IDIM_EFF,KFLEV))
      ELSE
        ALLOCATE(  ZEFCL_RRTM_SPLIT(0,0))
        ALLOCATE(  ZCLSW_TOTAL_SPLIT(0,0))
        ALLOCATE(  ZTAU_TOTAL_SPLIT(0,0,0))
        ALLOCATE(  ZOMEGA_TOTAL_SPLIT(0,0,0))
        ALLOCATE(  ZCG_TOTAL_SPLIT(0,0,0))
        ALLOCATE(  ZEFCL_LWD_SPLIT(0,0))
        ALLOCATE(  ZEFCL_LWU_SPLIT(0,0))
        ALLOCATE(  ZFLWP_SPLIT(0,0))
        ALLOCATE(  ZFIWP_SPLIT(0,0))
        ALLOCATE(  ZRADIP_SPLIT(0,0))
      END IF
!
! C2R2 coupling
!
      IF (SIZE (ZCCT_C2R2) > 0)  THEN
        ALLOCATE (ZCCT_C2R2_SPLIT(IDIM_EFF,KFLEV))
      ELSE
        ALLOCATE (ZCCT_C2R2_SPLIT(0,0))
      END IF
!
      IF (SIZE (ZCRT_C2R2) > 0)  THEN
        ALLOCATE (ZCRT_C2R2_SPLIT(IDIM_EFF,KFLEV))
      ELSE
        ALLOCATE (ZCRT_C2R2_SPLIT(0,0))
      END IF
!
      IF (SIZE (ZCIT_C1R3) > 0)  THEN
        ALLOCATE (ZCIT_C1R3_SPLIT(IDIM_EFF,KFLEV))
      ELSE
        ALLOCATE (ZCIT_C1R3_SPLIT(0,0))
      END IF
!
! LIMA coupling
!
      IF( CCLOUD == 'LIMA' ) THEN
          ALLOCATE (ZCCT_LIMA_SPLIT(IDIM_EFF,KFLEV))
          ALLOCATE (ZCRT_LIMA_SPLIT(IDIM_EFF,KFLEV))
          ALLOCATE (ZCIT_LIMA_SPLIT(IDIM_EFF,KFLEV))
      END IF
    END IF
! 
! fill the split arrays with their values taken from the full arrays 
!
    IBEG = IDIM-IDIM_RESIDUE+1
    IEND = IBEG+IDIM_EFF-1
!
    ZALBP_SPLIT(:,:) = ZALBP( IBEG:IEND ,:)
    ZALBD_SPLIT(:,:) = ZALBD( IBEG:IEND ,:)
    ZEMIS_SPLIT(:) = ZEMIS ( IBEG:IEND,1 )
    ZEMIW_SPLIT(:) = ZEMIW ( IBEG:IEND,1 )
    ZRMU0_SPLIT(:)    = ZRMU0 ( IBEG:IEND )
    ZLAT_SPLIT(:)    = ZLAT ( IBEG:IEND )
    ZLON_SPLIT(:)    = ZLON ( IBEG:IEND )
    ZCFAVE_SPLIT(:,:) = ZCFAVE( IBEG:IEND ,:)
    ZO3AVE_SPLIT(:,:) = ZO3AVE( IBEG:IEND ,:)
    ZT_HL_SPLIT(:,:)    = ZT_HL( IBEG:IEND ,:)
    ZPRES_HL_SPLIT(:,:) = ZPRES_HL( IBEG:IEND ,:)
    ZQLAVE_SPLIT(:,:) = ZQLAVE( IBEG:IEND , :)
    ZDZ_SPLIT(:,:) = ZDZ( IBEG:IEND , :)
    ZQIAVE_SPLIT(:,:) = ZQIAVE( IBEG:IEND ,:)
    ZQRAVE_SPLIT (:,:) = ZQRAVE (IBEG:IEND ,:)
    ZQLWC_SPLIT(:,:) = ZQLWC( IBEG:IEND , :)
    ZQIWC_SPLIT(:,:) = ZQIWC( IBEG:IEND ,:)
    ZQRWC_SPLIT(:,:) = ZQRWC (IBEG:IEND ,:)
    ZQVAVE_SPLIT(:,:) = ZQVAVE( IBEG:IEND ,:)
    ZTAVE_SPLIT(:,:)  = ZTAVE ( IBEG:IEND ,:)
    ZPAVE_SPLIT(:,:)  = ZPAVE ( IBEG:IEND ,:)
    ZAER_SPLIT (:,:,:)  = ZAER  ( IBEG:IEND ,:,:)
    IF(CAOP=='EXPL')THEN
       ZPIZA_EQ_SPLIT(:,:,:)=ZPIZA_EQ(IBEG:IEND,:,:)
       ZCGA_EQ_SPLIT(:,:,:)=ZCGA_EQ(IBEG:IEND,:,:)
       ZTAUREL_EQ_SPLIT(:,:,:)=ZTAUREL_EQ(IBEG:IEND,:,:)
    ENDIF
    ZDPRES_SPLIT(:,:)  = ZDPRES (IBEG:IEND ,:)
    ZLSM_SPLIT (:)    = ZLSM (IBEG:IEND)
    ZQSAVE_SPLIT (:,:) = ZQSAVE (IBEG:IEND ,:)
    ZTS_SPLIT (:) = ZTS (IBEG:IEND)
!
!  CALL the ECMWF radiation with the split array
!
  IF (CCLOUD == 'LIMA') THEN
! LIMA concentrations
     ZCCT_LIMA_SPLIT(:,:) = ZCCT_LIMA (IBEG:IEND ,:)
     ZCRT_LIMA_SPLIT(:,:) = ZCRT_LIMA (IBEG:IEND ,:)
     ZCIT_LIMA_SPLIT(:,:) = ZCIT_LIMA (IBEG:IEND ,:)
     
   IF (CRAD == "ECMW") THEN  
!
        CALL ECMWF_RADIATION_VERS2  ( IDIM_EFF , KFLEV, KRAD_DIAG, KAER,               &
                ZDZ_SPLIT,HEFRADL,HEFRADI,HOPWSW, HOPISW, HOPWLW, HOPILW,PFUDG,          &
                ZRII0, ZAER_SPLIT , ZALBD_SPLIT, ZALBP_SPLIT, ZPRES_HL_SPLIT,            &
                ZPAVE_SPLIT,PCCO2, ZCFAVE_SPLIT, ZDPRES_SPLIT, ZEMIS_SPLIT, ZEMIW_SPLIT, &
                ZLSM_SPLIT, ZRMU0_SPLIT,ZO3AVE_SPLIT , ZQVAVE_SPLIT, ZQIAVE_SPLIT ,ZQIWC_SPLIT,      &
                ZQLAVE_SPLIT,ZQLWC_SPLIT,ZQSAVE_SPLIT, ZQRAVE_SPLIT,ZQRWC_SPLIT,  ZT_HL_SPLIT,      &
                ZTAVE_SPLIT, ZTS_SPLIT, ZCCT_LIMA_SPLIT,ZCRT_LIMA_SPLIT,ZCIT_LIMA_SPLIT, &
                ZNFLW_CS_SPLIT, ZNFLW_SPLIT, ZNFSW_CS_SPLIT,ZNFSW_SPLIT,                 &
                ZDTLW_SPLIT, ZDTSW_SPLIT, ZFLUX_TOP_GND_IRVISNIR_SPLIT,                  &
                ZSFSWDIR_SPLIT, ZSFSWDIF_SPLIT,                                          &
                ZFLUX_SW_DOWN_SPLIT, ZFLUX_SW_UP_SPLIT, ZFLUX_LW_SPLIT ,                 &
                ZDTLW_CS_SPLIT, ZDTSW_CS_SPLIT, ZFLUX_TOP_GND_IRVISNIR_CS_SPLIT,         &
                ZFLUX_SW_DOWN_CS_SPLIT, ZFLUX_SW_UP_CS_SPLIT, ZFLUX_LW_CS_SPLIT,         &
                ZPLAN_ALB_VIS_SPLIT,ZPLAN_ALB_NIR_SPLIT, ZPLAN_TRA_VIS_SPLIT,            &
                ZPLAN_TRA_NIR_SPLIT, ZPLAN_ABS_VIS_SPLIT, ZPLAN_ABS_NIR_SPLIT,           &
                ZEFCL_LWD_SPLIT, ZEFCL_LWU_SPLIT, ZFLWP_SPLIT,ZFIWP_SPLIT,               &
                ZRADLP_SPLIT,ZRADIP_SPLIT,ZEFCL_RRTM_SPLIT, ZCLSW_TOTAL_SPLIT,           &
                ZTAU_TOTAL_SPLIT,ZOMEGA_TOTAL_SPLIT, ZCG_TOTAL_SPLIT,                    &
                GAOP,ZPIZA_EQ_SPLIT,ZCGA_EQ_SPLIT,ZTAUREL_EQ_SPLIT  )
                
   ELSE IF (CRAD == "ECRA") THEN  
        CALL ECRAD_INTERFACE  ( IDIM_EFF ,KFLEV, KRAD_DIAG, KAER,                   &      
            ZDZ,HEFRADL,HEFRADI,HOPWSW, HOPISW, HOPWLW, HOPILW,PFUDG,      &
            ZRII0, ZAER_SPLIT , ZALBD_SPLIT, ZALBP_SPLIT, ZPRES_HL_SPLIT, ZPAVE_SPLIT,               &
            PCCO2, ZCFAVE_SPLIT, ZDPRES_SPLIT, ZEMIS_SPLIT, ZEMIW_SPLIT, ZLSM_SPLIT, ZRMU0_SPLIT,          &
            ZO3AVE_SPLIT , ZQVAVE_SPLIT, ZQIAVE_SPLIT ,ZQIWC_SPLIT,ZQLAVE_SPLIT,ZQLWC_SPLIT, &
            ZQSAVE_SPLIT, ZQRAVE_SPLIT,  ZQRWC_SPLIT,  &
            ZT_HL_SPLIT,ZTAVE_SPLIT, ZTS_SPLIT, ZCCT_LIMA_SPLIT, &
            ZCRT_LIMA_SPLIT, ZCIT_LIMA_SPLIT,         &
            ZNFLW_SPLIT, ZNFSW_SPLIT, ZNFLW_CS_SPLIT, ZNFSW_CS_SPLIT,                           &
            ZDTLW_SPLIT, ZDTSW_SPLIT, ZFLUX_TOP_GND_IRVISNIR_SPLIT,                      &
            ZSFSWDIR_SPLIT, ZSFSWDIF_SPLIT,                                        &
            ZFLUX_SW_DOWN_SPLIT, ZFLUX_SW_UP_SPLIT, ZFLUX_LW_SPLIT ,                     &
            ZDTLW_CS_SPLIT, ZDTSW_CS_SPLIT, ZFLUX_TOP_GND_IRVISNIR_CS_SPLIT,             &
            ZFLUX_SW_DOWN_CS_SPLIT, ZFLUX_SW_UP_CS_SPLIT, ZFLUX_LW_CS_SPLIT,             &           
            ZPLAN_ALB_VIS_SPLIT,ZPLAN_ALB_NIR_SPLIT, ZPLAN_TRA_VIS_SPLIT, ZPLAN_TRA_NIR_SPLIT, &
            ZPLAN_ABS_VIS_SPLIT, ZPLAN_ABS_NIR_SPLIT,ZEFCL_LWD_SPLIT, ZEFCL_LWU_SPLIT,         &
            ZFLWP_SPLIT, ZFIWP_SPLIT,ZRADLP_SPLIT, ZRADIP_SPLIT, &
            ZEFCL_RRTM_SPLIT,  ZCLSW_TOTAL_SPLIT,  ZTAU_TOTAL_SPLIT,  &
            ZOMEGA_TOTAL_SPLIT,ZCG_TOTAL_SPLIT,                                    &
            GAOP, ZPIZA_EQ_SPLIT,ZCGA_EQ_SPLIT,ZTAUREL_EQ_SPLIT,ZLAT_SPLIT,ZLON_SPLIT    )
  END IF              
  ELSE
! C2R2 concentrations
    IF (SIZE (ZCCT_C2R2) > 0)  ZCCT_C2R2_SPLIT(:,:) = ZCCT_C2R2 (IBEG:IEND ,:)
    IF (SIZE (ZCRT_C2R2) > 0)  ZCRT_C2R2_SPLIT(:,:) = ZCRT_C2R2 (IBEG:IEND ,:)  
    IF (SIZE (ZCIT_C1R3) > 0)  ZCIT_C1R3_SPLIT(:,:) = ZCIT_C1R3 (IBEG:IEND ,:)
    IF (CRAD == "ECMW") THEN   
        CALL ECMWF_RADIATION_VERS2  ( IDIM_EFF , KFLEV, KRAD_DIAG, KAER,              &    
                ZDZ_SPLIT,HEFRADL,HEFRADI,HOPWSW, HOPISW, HOPWLW, HOPILW,PFUDG,                    &
                ZRII0, ZAER_SPLIT , ZALBD_SPLIT, ZALBP_SPLIT, ZPRES_HL_SPLIT,            &
                ZPAVE_SPLIT,PCCO2, ZCFAVE_SPLIT, ZDPRES_SPLIT, ZEMIS_SPLIT, ZEMIW_SPLIT, &
                ZLSM_SPLIT, ZRMU0_SPLIT,ZO3AVE_SPLIT , ZQVAVE_SPLIT, ZQIAVE_SPLIT ,ZQIWC_SPLIT,      & 
                ZQLAVE_SPLIT,ZQLWC_SPLIT,ZQSAVE_SPLIT, ZQRAVE_SPLIT,ZQRWC_SPLIT,  ZT_HL_SPLIT,      &
                ZTAVE_SPLIT, ZTS_SPLIT, ZCCT_C2R2_SPLIT,ZCRT_C2R2_SPLIT,ZCIT_C1R3_SPLIT, & 
                ZNFLW_CS_SPLIT, ZNFLW_SPLIT, ZNFSW_CS_SPLIT,ZNFSW_SPLIT,                 &          
                ZDTLW_SPLIT, ZDTSW_SPLIT, ZFLUX_TOP_GND_IRVISNIR_SPLIT,                  &
                ZSFSWDIR_SPLIT, ZSFSWDIF_SPLIT,                                          &
                ZFLUX_SW_DOWN_SPLIT, ZFLUX_SW_UP_SPLIT, ZFLUX_LW_SPLIT ,                 &
                ZDTLW_CS_SPLIT, ZDTSW_CS_SPLIT, ZFLUX_TOP_GND_IRVISNIR_CS_SPLIT,         &
                ZFLUX_SW_DOWN_CS_SPLIT, ZFLUX_SW_UP_CS_SPLIT, ZFLUX_LW_CS_SPLIT,         & 
                ZPLAN_ALB_VIS_SPLIT,ZPLAN_ALB_NIR_SPLIT, ZPLAN_TRA_VIS_SPLIT,            &
                ZPLAN_TRA_NIR_SPLIT, ZPLAN_ABS_VIS_SPLIT, ZPLAN_ABS_NIR_SPLIT,           &
                ZEFCL_LWD_SPLIT, ZEFCL_LWU_SPLIT, ZFLWP_SPLIT,ZFIWP_SPLIT,               &
                ZRADLP_SPLIT,ZRADIP_SPLIT,ZEFCL_RRTM_SPLIT, ZCLSW_TOTAL_SPLIT,           &
                ZTAU_TOTAL_SPLIT,ZOMEGA_TOTAL_SPLIT, ZCG_TOTAL_SPLIT,                    &
                GAOP,ZPIZA_EQ_SPLIT,ZCGA_EQ_SPLIT,ZTAUREL_EQ_SPLIT  )
         
    ELSE IF (CRAD == "ECRA") THEN  
       CALL ECRAD_INTERFACE  ( IDIM_EFF ,KFLEV, KRAD_DIAG, KAER,                   &      
            ZDZ_SPLIT,HEFRADL,HEFRADI,HOPWSW, HOPISW, HOPWLW, HOPILW,PFUDG,      &
            ZRII0, ZAER_SPLIT , ZALBD_SPLIT, ZALBP_SPLIT, ZPRES_HL_SPLIT, ZPAVE_SPLIT,               &
            PCCO2, ZCFAVE_SPLIT, ZDPRES_SPLIT, ZEMIS_SPLIT, ZEMIW_SPLIT, ZLSM_SPLIT, ZRMU0_SPLIT,          &
            ZO3AVE_SPLIT , ZQVAVE_SPLIT, ZQIAVE_SPLIT ,ZQIWC_SPLIT,ZQLAVE_SPLIT,ZQLWC_SPLIT, &
            ZQSAVE_SPLIT, ZQRAVE_SPLIT,  ZQRWC_SPLIT,  &
            ZT_HL_SPLIT,ZTAVE_SPLIT, ZTS_SPLIT, ZCCT_C2R2_SPLIT, &
            ZCRT_C2R2_SPLIT, ZCIT_C1R3_SPLIT,         &
            ZNFLW_SPLIT, ZNFSW_SPLIT, ZNFLW_CS_SPLIT, ZNFSW_CS_SPLIT,                           &
            ZDTLW_SPLIT, ZDTSW_SPLIT, ZFLUX_TOP_GND_IRVISNIR_SPLIT,                      &
            ZSFSWDIR_SPLIT, ZSFSWDIF_SPLIT,                                        &
            ZFLUX_SW_DOWN_SPLIT, ZFLUX_SW_UP_SPLIT, ZFLUX_LW_SPLIT ,                     &
            ZDTLW_CS_SPLIT, ZDTSW_CS_SPLIT, ZFLUX_TOP_GND_IRVISNIR_CS_SPLIT,             &
            ZFLUX_SW_DOWN_CS_SPLIT, ZFLUX_SW_UP_CS_SPLIT, ZFLUX_LW_CS_SPLIT,             &           
            ZPLAN_ALB_VIS_SPLIT,ZPLAN_ALB_NIR_SPLIT, ZPLAN_TRA_VIS_SPLIT, ZPLAN_TRA_NIR_SPLIT, &
            ZPLAN_ABS_VIS_SPLIT, ZPLAN_ABS_NIR_SPLIT,ZEFCL_LWD_SPLIT, ZEFCL_LWU_SPLIT,         &
            ZFLWP_SPLIT, ZFIWP_SPLIT,ZRADLP_SPLIT, ZRADIP_SPLIT, &
            ZEFCL_RRTM_SPLIT,  ZCLSW_TOTAL_SPLIT,  ZTAU_TOTAL_SPLIT,  &
            ZOMEGA_TOTAL_SPLIT,ZCG_TOTAL_SPLIT,                                    &
            GAOP, ZPIZA_EQ_SPLIT,ZCGA_EQ_SPLIT,ZTAUREL_EQ_SPLIT,ZLAT_SPLIT,ZLON_SPLIT    )
    END IF                   
    END IF 
!
! fill the full output arrays with the split arrays
!
    ZDTLW( IBEG:IEND ,:)  =  ZDTLW_SPLIT(:,:)  
    ZDTSW( IBEG:IEND ,:)  =  ZDTSW_SPLIT(:,:) 
    ZFLUX_TOP_GND_IRVISNIR( IBEG:IEND ,:)=  ZFLUX_TOP_GND_IRVISNIR_SPLIT(:,:) 
    ZSFSWDIR (IBEG:IEND,:)  = ZSFSWDIR_SPLIT(:,:)
    ZSFSWDIF (IBEG:IEND,:)  = ZSFSWDIF_SPLIT(:,:)
!
    ZDTLW_CS( IBEG:IEND ,:) =  ZDTLW_CS_SPLIT(:,:)
    ZDTSW_CS( IBEG:IEND ,:) =  ZDTSW_CS_SPLIT(:,:)
    ZFLUX_TOP_GND_IRVISNIR_CS( IBEG:IEND ,:) =                     &
         ZFLUX_TOP_GND_IRVISNIR_CS_SPLIT(:,:)
    ZFLUX_LW( IBEG:IEND ,:,:)    =  ZFLUX_LW_SPLIT(:,:,:) 
    ZFLUX_SW_DOWN( IBEG:IEND ,:) =  ZFLUX_SW_DOWN_SPLIT(:,:)
    ZFLUX_SW_UP( IBEG:IEND ,:)   =  ZFLUX_SW_UP_SPLIT(:,:)
    ZRADLP( IBEG:IEND ,:) = ZRADLP_SPLIT(:,:)
    IF ( tpfile%lopened ) THEN
      IF( KRAD_DIAG >= 1) THEN
        ZNFLW(IBEG:IEND ,:)= ZNFLW_SPLIT(:,:)
        ZNFSW(IBEG:IEND ,:)= ZNFSW_SPLIT(:,:)
        IF( KRAD_DIAG >= 2) THEN
          ZFLUX_SW_DOWN_CS( IBEG:IEND ,:) = ZFLUX_SW_DOWN_CS_SPLIT(:,:)
          ZFLUX_SW_UP_CS( IBEG:IEND ,:)   = ZFLUX_SW_UP_CS_SPLIT(:,:)
          ZFLUX_LW_CS( IBEG:IEND ,:,:)    = ZFLUX_LW_CS_SPLIT(:,:,:)
          ZNFLW_CS(IBEG:IEND ,:)= ZNFLW_CS_SPLIT(:,:)
          ZNFSW_CS(IBEG:IEND ,:)= ZNFSW_CS_SPLIT(:,:)
          IF( KRAD_DIAG >= 3) THEN
            ZPLAN_ALB_VIS( IBEG:IEND ) = ZPLAN_ALB_VIS_SPLIT(:)
            ZPLAN_ALB_NIR( IBEG:IEND ) = ZPLAN_ALB_NIR_SPLIT(:)
            ZPLAN_TRA_VIS( IBEG:IEND ) = ZPLAN_TRA_VIS_SPLIT(:)
            ZPLAN_TRA_NIR( IBEG:IEND ) = ZPLAN_TRA_NIR_SPLIT(:)
            ZPLAN_ABS_VIS( IBEG:IEND ) = ZPLAN_ABS_VIS_SPLIT(:)
            ZPLAN_ABS_NIR( IBEG:IEND ) = ZPLAN_ABS_NIR_SPLIT(:)          
            IF( KRAD_DIAG >= 4) THEN
              ZEFCL_LWD( IBEG:IEND ,:) = ZEFCL_LWD_SPLIT(:,:)
              ZEFCL_LWU( IBEG:IEND ,:)   = ZEFCL_LWU_SPLIT(:,:)
              ZFLWP( IBEG:IEND ,:) = ZFLWP_SPLIT(:,:)
              ZFIWP( IBEG:IEND ,:) = ZFIWP_SPLIT(:,:)
              ZRADIP( IBEG:IEND ,:) = ZRADIP_SPLIT(:,:)
              ZEFCL_RRTM( IBEG:IEND ,:) = ZEFCL_RRTM_SPLIT(:,:)
              ZCLSW_TOTAL( IBEG:IEND ,:) = ZCLSW_TOTAL_SPLIT(:,:)
              ZTAU_TOTAL( IBEG:IEND ,:,:)  = ZTAU_TOTAL_SPLIT(:,:,:)
              ZOMEGA_TOTAL( IBEG:IEND ,:,:)= ZOMEGA_TOTAL_SPLIT(:,:,:)
              ZCG_TOTAL( IBEG:IEND ,:,:)   = ZCG_TOTAL_SPLIT(:,:,:)                
            END IF
          END IF
        END IF
      END IF
    END IF
!
    IDIM_RESIDUE = IDIM_RESIDUE - IDIM_EFF
!
! desallocation of the split arrays
!
    IF( JI_SPLIT >= INUM_CALL-1 ) THEN
      DEALLOCATE(  ZALBP_SPLIT )
      DEALLOCATE(  ZALBD_SPLIT )  
      DEALLOCATE(  ZEMIS_SPLIT  )
      DEALLOCATE(  ZEMIW_SPLIT  )
      DEALLOCATE(  ZLAT_SPLIT  )
      DEALLOCATE(  ZLON_SPLIT  )
      DEALLOCATE(  ZRMU0_SPLIT      )
      DEALLOCATE(  ZCFAVE_SPLIT     )
      DEALLOCATE(  ZO3AVE_SPLIT     )
      DEALLOCATE(  ZT_HL_SPLIT      )
      DEALLOCATE(  ZPRES_HL_SPLIT   )
      DEALLOCATE(  ZDZ_SPLIT     )
      DEALLOCATE(  ZQLAVE_SPLIT     )
      DEALLOCATE(  ZQIAVE_SPLIT     )
      DEALLOCATE(  ZQVAVE_SPLIT     )
      DEALLOCATE(  ZTAVE_SPLIT      )
      DEALLOCATE(  ZPAVE_SPLIT      )
      DEALLOCATE(  ZAER_SPLIT       )
      DEALLOCATE(  ZDPRES_SPLIT     )
      DEALLOCATE(  ZLSM_SPLIT       )
      DEALLOCATE(  ZQSAVE_SPLIT     )
      DEALLOCATE(  ZQRAVE_SPLIT  )
      DEALLOCATE(  ZQLWC_SPLIT     )
      DEALLOCATE(  ZQRWC_SPLIT     )
      DEALLOCATE(  ZQIWC_SPLIT     )
      IF ( ALLOCATED( ZCCT_C2R2_SPLIT ) ) DEALLOCATE(   ZCCT_C2R2_SPLIT  )
      IF ( ALLOCATED( ZCRT_C2R2_SPLIT ) ) DEALLOCATE(   ZCRT_C2R2_SPLIT  )
      IF ( ALLOCATED( ZCIT_C1R3_SPLIT ) ) DEALLOCATE(   ZCIT_C1R3_SPLIT  )
      IF ( ALLOCATED( ZCCT_LIMA_SPLIT ) ) DEALLOCATE(   ZCCT_LIMA_SPLIT  )
      IF ( ALLOCATED( ZCRT_LIMA_SPLIT ) ) DEALLOCATE(   ZCRT_LIMA_SPLIT  )
      IF ( ALLOCATED( ZCIT_LIMA_SPLIT ) ) DEALLOCATE(   ZCIT_LIMA_SPLIT  )
      DEALLOCATE(  ZTS_SPLIT    )
      DEALLOCATE(   ZNFLW_CS_SPLIT)
      DEALLOCATE(   ZNFLW_SPLIT)
      DEALLOCATE(   ZNFSW_CS_SPLIT)
      DEALLOCATE(   ZNFSW_SPLIT)
      DEALLOCATE(ZDTLW_SPLIT)
      DEALLOCATE(ZDTSW_SPLIT)
      DEALLOCATE(ZFLUX_TOP_GND_IRVISNIR_SPLIT)
      DEALLOCATE(ZSFSWDIR_SPLIT)
      DEALLOCATE(ZSFSWDIF_SPLIT)
      DEALLOCATE(ZFLUX_SW_DOWN_SPLIT)
      DEALLOCATE(ZFLUX_SW_UP_SPLIT)
      DEALLOCATE(ZFLUX_LW_SPLIT)
      DEALLOCATE(ZDTLW_CS_SPLIT)
      DEALLOCATE(ZDTSW_CS_SPLIT)
      DEALLOCATE(ZFLUX_TOP_GND_IRVISNIR_CS_SPLIT)
      DEALLOCATE(ZPLAN_ALB_VIS_SPLIT)
      DEALLOCATE(ZPLAN_ALB_NIR_SPLIT)
      DEALLOCATE(ZPLAN_TRA_VIS_SPLIT)
      DEALLOCATE(ZPLAN_TRA_NIR_SPLIT)
      DEALLOCATE(ZPLAN_ABS_VIS_SPLIT)
      DEALLOCATE(ZPLAN_ABS_NIR_SPLIT)
      DEALLOCATE(ZEFCL_LWD_SPLIT)
      DEALLOCATE(ZEFCL_LWU_SPLIT)
      DEALLOCATE(ZFLWP_SPLIT)
      DEALLOCATE(ZRADLP_SPLIT)
      DEALLOCATE(ZRADIP_SPLIT)
      DEALLOCATE(ZFIWP_SPLIT)
      DEALLOCATE(ZEFCL_RRTM_SPLIT)
      DEALLOCATE(ZCLSW_TOTAL_SPLIT)
      DEALLOCATE(ZTAU_TOTAL_SPLIT)
      DEALLOCATE(ZOMEGA_TOTAL_SPLIT)
      DEALLOCATE(ZCG_TOTAL_SPLIT)
      DEALLOCATE(ZFLUX_SW_DOWN_CS_SPLIT)
      DEALLOCATE(ZFLUX_SW_UP_CS_SPLIT)
      DEALLOCATE(ZFLUX_LW_CS_SPLIT)
      DEALLOCATE(ZPIZA_EQ_SPLIT)
      DEALLOCATE(ZCGA_EQ_SPLIT)
      DEALLOCATE(ZTAUREL_EQ_SPLIT)
    END IF
  END DO
END IF

!
DEALLOCATE(ZTAVE)
DEALLOCATE(ZPAVE)
DEALLOCATE(ZQVAVE)
DEALLOCATE(ZQLAVE)
DEALLOCATE(ZDZ)
DEALLOCATE(ZQIAVE)
DEALLOCATE(ZCFAVE)
DEALLOCATE(ZPRES_HL)
DEALLOCATE(ZT_HL)
DEALLOCATE(ZRMU0) 
DEALLOCATE(ZLSM)
DEALLOCATE(ZQSAVE)
DEALLOCATE(ZAER)
DEALLOCATE(ZPIZA_EQ)
DEALLOCATE(ZCGA_EQ)
DEALLOCATE(ZTAUREL_EQ)
DEALLOCATE(ZDPRES)
DEALLOCATE(ZCCT_C2R2)
DEALLOCATE(ZCRT_C2R2)
DEALLOCATE(ZCIT_C1R3)
DEALLOCATE(ZLAT)
DEALLOCATE(ZLON)
IF (CCLOUD == 'LIMA') THEN 
  DEALLOCATE(ZCCT_LIMA)
  DEALLOCATE(ZCRT_LIMA)
  DEALLOCATE(ZCIT_LIMA)
END IF
!
DEALLOCATE(ZTS)
DEALLOCATE(ZALBP)
DEALLOCATE(ZALBD)
DEALLOCATE(ZEMIS)
DEALLOCATE(ZEMIW)
DEALLOCATE(ZQRAVE)
DEALLOCATE(ZQLWC)
DEALLOCATE(ZQIWC)
DEALLOCATE(ZQRWC)
DEALLOCATE(ICLEAR_2D_TM1)
!
!*       5.6   UNCOMPRESSES THE OUTPUT FIELD IN CASE OF 
!                      CLEAR-SKY APPROXIMATION
!
IF(OCLEAR_SKY .OR. OCLOUD_ONLY) THEN
  ALLOCATE(ZWORK1(ICLOUD))
  ALLOCATE(ZWORK2(ICLOUD+KFLEV)) !       allocation for the KFLEV levels of 
  ALLOCATE(ZWORK4(KFLEV,KDLON))
  ZWORK2(:) = PACK( TRANSPOSE(ZDTLW(:,:)),MASK=.TRUE. )
!
  DO JK=1,KFLEV
    ZWORK4(JK,:) = ZWORK2(ICLOUD+JK)
  END DO
  ZWORK1(1:ICLOUD) = ZWORK2(1:ICLOUD)
  ZZDTLW(:,:) = TRANSPOSE( UNPACK( ZWORK1(:),MASK=GCLOUDT(:,:)   &
       ,FIELD=ZWORK4(:,:) ) )
  !
  ZWORK2(:) = PACK( TRANSPOSE(ZDTSW(:,:)),MASK=.TRUE. )
  DO JK=1,KFLEV
    ZWORK4(JK,:) = ZWORK2(ICLOUD+JK)
  END DO
  ZWORK1(1:ICLOUD) = ZWORK2(1:ICLOUD)
  ZZDTSW(:,:) = TRANSPOSE( UNPACK( ZWORK1(:),MASK=GCLOUDT(:,:)   &
       ,FIELD=ZWORK4(:,:) ) )
  !
  DEALLOCATE(ZWORK1)
  DEALLOCATE(ZWORK2)
  DEALLOCATE(ZWORK4)
  !
  ZZTGVISC   = ZFLUX_TOP_GND_IRVISNIR(ICLOUD_COL+1,5)
  !
  ZZTGVIS(:) = UNPACK( ZFLUX_TOP_GND_IRVISNIR(:,5),MASK=.NOT.GCLEAR_2D(:), &
       FIELD=ZZTGVISC  )
  ZZTGNIRC   = ZFLUX_TOP_GND_IRVISNIR(ICLOUD_COL+1,6)
  !
  ZZTGNIR(:) = UNPACK( ZFLUX_TOP_GND_IRVISNIR(:,6),MASK=.NOT.GCLEAR_2D(:), &
       FIELD=ZZTGNIRC )
  ZZTGIRC    = ZFLUX_TOP_GND_IRVISNIR(ICLOUD_COL+1,4)
  !
  ZZTGIR (:) = UNPACK( ZFLUX_TOP_GND_IRVISNIR(:,4),MASK=.NOT.GCLEAR_2D(:), &
       FIELD=ZZTGIRC  )
  !
  DO JSWB=1,ISWB
    ZZSFSWDIRC(JSWB) = ZSFSWDIR (ICLOUD_COL+1,JSWB)
    !
    ZZSFSWDIR(:,JSWB) =  UNPACK(ZSFSWDIR (:,JSWB),MASK=.NOT.GCLEAR_2D(:), &
         FIELD= ZZSFSWDIRC(JSWB)  ) 
    !
    ZZSFSWDIFC(JSWB) = ZSFSWDIF (ICLOUD_COL+1,JSWB)
    !
    ZZSFSWDIF(:,JSWB) =  UNPACK(ZSFSWDIF (:,JSWB),MASK=.NOT.GCLEAR_2D(:), &
         FIELD= ZZSFSWDIFC(JSWB)  )
  END DO
!
!  No cloud case
!
  IF( GNOCL ) THEN
          IF (SIZE(ZZDTLW,1)>1) THEN
             ZZDTLW(1,:)= ZZDTLW(2,:)
          ENDIF
          IF (SIZE(ZZDTSW,1)>1) THEN
             ZZDTSW(1,:)= ZZDTSW(2,:)
          ENDIF
    ZZTGVIS(1) = ZZTGVISC
    ZZTGNIR(1) = ZZTGNIRC
    ZZTGIR(1)  = ZZTGIRC
    ZZSFSWDIR(1,:) =  ZZSFSWDIRC(:)
    ZZSFSWDIF(1,:) =  ZZSFSWDIFC(:)
  END IF
ELSE
  ZZDTLW(:,:) = ZDTLW(:,:)
  ZZDTSW(:,:) = ZDTSW(:,:)
  ZZTGVIS(:)  = ZFLUX_TOP_GND_IRVISNIR(:,5)
  ZZTGNIR(:)  = ZFLUX_TOP_GND_IRVISNIR(:,6)
  ZZTGIR(:)   = ZFLUX_TOP_GND_IRVISNIR(:,4)
  ZZSFSWDIR(:,:) =  ZSFSWDIR(:,:)
  ZZSFSWDIF(:,:) =  ZSFSWDIF(:,:) 
END IF
!
DEALLOCATE(ZDTLW)
DEALLOCATE(ZDTSW)
DEALLOCATE(ZSFSWDIR)
DEALLOCATE(ZSFSWDIF)
!
!--------------------------------------------------------------------------------------------
!
!*       6.    COMPUTES THE RADIATIVE SOURCES AND THE DOWNWARD SURFACE FLUXES in 2D horizontal 
!              ------------------------------------------------------------------------------
!
!  Computes the SW and LW radiative tendencies
!  note : tendencies in K/s for MNH (from K/day)
!
ZDTRAD_LW(:,:,:)=0.0
ZDTRAD_SW(:,:,:)=0.0
DO JK=IKB,IKE
  JKRAD= JK-JPVEXT
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZDTRAD_LW(JI,JJ,JK) = ZZDTLW(IIJ,JKRAD)/XDAY  ! XDAY from  modd_cst (day duration in s)
      ZDTRAD_SW(JI,JJ,JK) = ZZDTSW(IIJ,JKRAD)/XDAY      
    END DO
  END DO
END DO
!
!  Computes the downward SW and LW surface fluxes + diffuse and direct contribution
!
ZLWD(:,:)=0.
ZSWDDIR(:,:,:)=0.
ZSWDDIF(:,:,:)=0.
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZLWD(JI,JJ) = ZZTGIR(IIJ)
    ZSWDDIR(JI,JJ,:) = ZZSFSWDIR (IIJ,:)
    ZSWDDIF(JI,JJ,:) = ZZSFSWDIF (IIJ,:)
  END DO
END DO
!
!final  THETA_radiative tendency and surface fluxes 
!
IF(OCLOUD_ONLY) THEN

  GCLOUD_SURF(:,:) = .FALSE.
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        GCLOUD_SURF(JI,JJ) = GCLOUD(IIJ,1)
    END DO
  END DO
 
  ZWORKL(:,:) = GCLOUD_SURF(:,:)

  DO JK = IKB,IKE
    WHERE( ZWORKL(:,:) )
      PDTHRAD(:,:,JK) = (ZDTRAD_LW(:,:,JK)+ZDTRAD_SW(:,:,JK))/ZEXNT(:,:,JK)
    ENDWHERE
  END DO
  !
  WHERE( ZWORKL(:,:) )
    PSRFLWD(:,:) = ZLWD(:,:)
  ENDWHERE
  DO JSWB=1,ISWB
    WHERE( ZWORKL(:,:) )
      PSRFSWD_DIR (:,:,JSWB) = ZSWDDIR(:,:,JSWB)
      PSRFSWD_DIF (:,:,JSWB) = ZSWDDIF(:,:,JSWB)
    END WHERE
  END DO
ELSE
  PDTHRAD(:,:,:) = (ZDTRAD_LW(:,:,:)+ZDTRAD_SW(:,:,:))/ZEXNT(:,:,:)  ! tendency in potential temperature
  PDTHRADSW(:,:,:) = ZDTRAD_SW(:,:,:)/ZEXNT(:,:,:)
  PDTHRADLW(:,:,:) = ZDTRAD_LW(:,:,:)/ZEXNT(:,:,:)
  PSRFLWD(:,:) = ZLWD(:,:)
  DO JSWB=1,ISWB
    PSRFSWD_DIR (:,:,JSWB) = ZSWDDIR(:,:,JSWB)
    PSRFSWD_DIF (:,:,JSWB) = ZSWDDIF(:,:,JSWB)
  END DO
!
!sw and lw fluxes 
!
  DO JK=IKB,IKE
   JKRAD = JK - JPVEXT
   DO JJ=IJB,IJE
    DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          PSWU(JI,JJ,JK) = ZFLUX_SW_UP(IIJ,JKRAD)
          PSWD(JI,JJ,JK) = ZFLUX_SW_DOWN(IIJ,JKRAD)
          PLWU(JI,JJ,JK) = ZFLUX_LW(IIJ,1,JKRAD)
          PLWD(JI,JJ,JK) = -ZFLUX_LW(IIJ,2,JKRAD)  ! in ECMWF all fluxes are upward
    END DO
   END DO
  END DO
!!!effective radius
  DO JK=IKB,IKE
   JKRAD = JK - JPVEXT
   DO JJ=IJB,IJE
    DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          PRADEFF(JI,JJ,JK) = ZRADLP(IIJ,JKRAD)
    END DO
   END DO
  END DO
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       7.    STORE SOME ADDITIONNAL RADIATIVE FIELDS
!              ---------------------------------------
!
IF( tpfile%lopened .AND. (KRAD_DIAG >= 1) ) THEN
  ZSTORE_3D(:,:,:) = 0.0
  ZSTORE_3D2(:,:,:) = 0.0
  ZSTORE_2D(:,:)   = 0.0
  !
  TZFIELD2D = TFIELDMETADATA(                 &
    CMNHNAME   = 'generic 2D for radiations', & !Temporary name to ease identification
    CSTDNAME   = '',                          &
    CDIR       = 'XY',                        &
    NGRID      = 1,                           &
    NTYPE      = TYPEREAL,                    &
    NDIMS      = 2,                           &
    LTIMEDEP   = .TRUE.                       )

  TZFIELD3D = TFIELDMETADATA(                 &
    CMNHNAME   = 'generic 3D for radiations', & !Temporary name to ease identification
    CSTDNAME   = '',                          &
    CDIR       = 'XY',                        &
    NGRID      = 1,                           &
    NTYPE      = TYPEREAL,                    &
    NDIMS      = 3,                           &
    LTIMEDEP   = .TRUE.                       )

  IF( KRAD_DIAG >= 1) THEN
    !
    ILUOUT = TLUOUT%NLU
    WRITE(UNIT=ILUOUT,FMT='(/," STORE ADDITIONNAL RADIATIVE FIELDS:", &
         & " KRAD_DIAG=",I1,/)') KRAD_DIAG
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZFLUX_SW_DOWN(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'SWF_DOWN'
    TZFIELD3D%CLONGNAME  = 'SWF_DOWN'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_SWF_DOWN'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
!
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZFLUX_SW_UP(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'SWF_UP'
    TZFIELD3D%CLONGNAME  = 'SWF_UP'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_SWF_UP'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
!
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = -ZFLUX_LW(IIJ,2,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'LWF_DOWN'
    TZFIELD3D%CLONGNAME  = 'LWF_DOWN'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_LWF_DOWN'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
!
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZFLUX_LW(IIJ,1,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'LWF_UP'
    TZFIELD3D%CLONGNAME  = 'LWF_UP'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_LWF_UP'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
!
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZNFLW(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'LWF_NET'
    TZFIELD3D%CLONGNAME  = 'LWF_NET'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_LWF_NET'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
!
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZNFSW(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'SWF_NET'
    TZFIELD3D%CLONGNAME  = 'SWF_NET'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_SWF_NET'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
!
    DO JK=IKB,IKE
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZSTORE_3D(JI,JJ,JK) = ZDTRAD_LW (JI,JJ,JK)*XDAY
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'DTRAD_LW'
    TZFIELD3D%CLONGNAME  = 'DTRAD_LW'
    TZFIELD3D%CUNITS     = 'K day-1'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_DTRAD_LW'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
!
    DO JK=IKB,IKE
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZSTORE_3D(JI,JJ,JK) = ZDTRAD_SW (JI,JJ,JK)*XDAY
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'DTRAD_SW'
    TZFIELD3D%CLONGNAME  = 'DTRAD_SW'
    TZFIELD3D%CUNITS     = 'K day-1'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_DTRAD_SW'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
!
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZFLUX_TOP_GND_IRVISNIR(IIJ,5)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'RADSWD_VIS'
    TZFIELD2D%CLONGNAME  = 'RADSWD_VIS'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_RADSWD_VIS'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
!
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZFLUX_TOP_GND_IRVISNIR(IIJ,6)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'RADSWD_NIR'
    TZFIELD2D%CLONGNAME  = 'RADSWD_NIR'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_RADSWD_NIR'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
    !
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZFLUX_TOP_GND_IRVISNIR(IIJ,4)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'RADLWD'
    TZFIELD2D%CLONGNAME  = 'RADLWD'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_RADLWD'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
  END IF
  !
  !
  IF( KRAD_DIAG >= 2) THEN
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZFLUX_SW_DOWN_CS(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'SWF_DOWN_CS'
    TZFIELD3D%CLONGNAME  = 'SWF_DOWN_CS'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_SWF_DOWN_CS'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZFLUX_SW_UP_CS(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'SWF_UP_CS'
    TZFIELD3D%CLONGNAME  = 'SWF_UP_CS'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_SWF_UP_CS'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = -ZFLUX_LW_CS(IIJ,2,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'LWF_DOWN_CS'
    TZFIELD3D%CLONGNAME  = 'LWF_DOWN_CS'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_LWF_DOWN_CS'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZFLUX_LW_CS(IIJ,1,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'LWF_UP_CS'
    TZFIELD3D%CLONGNAME  = 'LWF_UP_CS'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_LWF_UP_CS'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZNFLW_CS(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'LWF_NET_CS'
    TZFIELD3D%CLONGNAME  = 'LWF_NET_CS'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_LWF_NET_CS'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZNFSW_CS(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'SWF_NET_CS'
    TZFIELD3D%CLONGNAME  = 'SWF_NET_CS'
    TZFIELD3D%CUNITS     = 'W m-2'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_SWF_NET_CS'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK-JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZDTSW_CS(IIJ,JKRAD) 
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'DTRAD_SW_CS'
    TZFIELD3D%CLONGNAME  = 'DTRAD_SW_CS'
    TZFIELD3D%CUNITS     = 'K day-1'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_DTRAD_SW_CS'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK-JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZDTLW_CS(IIJ,JKRAD) 
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'DTRAD_LW_CS'
    TZFIELD3D%CLONGNAME  = 'DTRAD_LW_CS'
    TZFIELD3D%CUNITS     = 'K day-1'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_DTRAD_LW_CS'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZFLUX_TOP_GND_IRVISNIR_CS(IIJ,5)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'RADSWD_VIS_CS'
    TZFIELD2D%CLONGNAME  = 'RADSWD_VIS_CS'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_RADSWD_VIS_CS'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
    !
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZFLUX_TOP_GND_IRVISNIR_CS(IIJ,6)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'RADSWD_NIR_CS'
    TZFIELD2D%CLONGNAME  = 'RADSWD_NIR_CS'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_RADSWD_NIR_CS'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
    !
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZFLUX_TOP_GND_IRVISNIR_CS(IIJ,4)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'RADLWD_CS'
    TZFIELD2D%CLONGNAME  = 'RADLWD_CS'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_RADLWD_CS'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
  END IF
  !
  !
  IF( KRAD_DIAG >= 3) THEN
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZPLAN_ALB_VIS(IIJ)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'PLAN_ALB_VIS'
    TZFIELD2D%CLONGNAME  = 'PLAN_ALB_VIS'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_PLAN_ALB_VIS'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
    !
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZPLAN_ALB_NIR(IIJ)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'PLAN_ALB_NIR'
    TZFIELD2D%CLONGNAME  = 'PLAN_ALB_NIR'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_PLAN_ALB_NIR'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
    !
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZPLAN_TRA_VIS(IIJ)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'PLAN_TRA_VIS'
    TZFIELD2D%CLONGNAME  = 'PLAN_TRA_VIS'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_PLAN_TRA_VIS'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
    !
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZPLAN_TRA_NIR(IIJ)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'PLAN_TRA_NIR'
    TZFIELD2D%CLONGNAME  = 'PLAN_TRA_NIR'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_PLAN_TRA_NIR'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
    !
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZPLAN_ABS_VIS(IIJ)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'PLAN_ABS_VIS'
    TZFIELD2D%CLONGNAME  = 'PLAN_ABS_VIS'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_PLAN_ABS_VIS'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
    !
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZSTORE_2D(JI,JJ) = ZPLAN_ABS_NIR(IIJ)
      END DO
    END DO
    TZFIELD2D%CMNHNAME   = 'PLAN_ABS_NIR'
    TZFIELD2D%CLONGNAME  = 'PLAN_ABS_NIR'
    TZFIELD2D%CUNITS     = ''
    TZFIELD2D%CCOMMENT   = 'X_Y_PLAN_ABS_NIR'
    CALL IO_Field_write(TPFILE,TZFIELD2D,ZSTORE_2D)
    !
    !
  END IF
!
!
  IF( KRAD_DIAG >= 4) THEN
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZEFCL_LWD(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'EFNEB_DOWN'
    TZFIELD3D%CLONGNAME  = 'EFNEB_DOWN'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_EFNEB_DOWN'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZEFCL_LWU(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'EFNEB_UP'
    TZFIELD3D%CLONGNAME  = 'EFNEB_UP'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_EFNEB_UP'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZFLWP(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'FLWP'
    TZFIELD3D%CLONGNAME  = 'FLWP'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_FLWP'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZFIWP(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'FIWP'
    TZFIELD3D%CLONGNAME  = 'FIWP'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_FIWP'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZRADLP(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'EFRADL'
    TZFIELD3D%CLONGNAME  = 'EFRADL'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_RAD_microm'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    ! 
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZRADIP(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'EFRADI'
    TZFIELD3D%CLONGNAME  = 'EFRADI'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_RAD_microm'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZCLSW_TOTAL(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'SW_NEB'
    TZFIELD3D%CLONGNAME  = 'SW_NEB'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_SW_NEB'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZEFCL_RRTM(IIJ,JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'RRTM_LW_NEB'
    TZFIELD3D%CLONGNAME  = 'RRTM_LW_NEB'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_LW_NEB'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    !
    ! spectral bands
    IF (KSWB_OLD==6) THEN
      INIR = 4
    ELSE
      INIR = 2
    END IF

    DO JBAND=1,INIR-1
      WRITE(YBAND_NAME(JBAND),'(A3,I1)') 'VIS', JBAND
    END DO
    DO JBAND= INIR, KSWB_OLD
      WRITE(YBAND_NAME(JBAND),'(A3,I1)') 'NIR', JBAND
    END DO
!
    DO JBAND=1,KSWB_OLD
      TZFIELD3D%CMNHNAME   = 'ODAER_'//YBAND_NAME(JBAND)
      TZFIELD3D%CLONGNAME  = 'ODAER_'//YBAND_NAME(JBAND)
      TZFIELD3D%CUNITS     = ''
      TZFIELD3D%CCOMMENT   = 'X_Y_Z_OD_'//YBAND_NAME(JBAND)
      CALL IO_Field_write(TPFILE,TZFIELD3D,ZTAUAZ(:,:,:,JBAND))
      !
      TZFIELD3D%CMNHNAME   = 'SSAAER_'//YBAND_NAME(JBAND)
      TZFIELD3D%CLONGNAME  = 'SSAAER_'//YBAND_NAME(JBAND)
      TZFIELD3D%CUNITS     = ''
      TZFIELD3D%CCOMMENT   = 'X_Y_Z_SSA_'//YBAND_NAME(JBAND)
      CALL IO_Field_write(TPFILE,TZFIELD3D,ZPIZAZ(:,:,:,JBAND))
      !
      TZFIELD3D%CMNHNAME   = 'GAER_'//YBAND_NAME(JBAND)
      TZFIELD3D%CLONGNAME  = 'GAER_'//YBAND_NAME(JBAND)
      TZFIELD3D%CUNITS     = ''
      TZFIELD3D%CCOMMENT   = 'X_Y_Z_G_'//YBAND_NAME(JBAND)
      CALL IO_Field_write(TPFILE,TZFIELD3D,ZCGAZ(:,:,:,JBAND))
    ENDDO

    DO JBAND=1,KSWB_OLD
      DO JK=IKB,IKE
        JKRAD = JK - JPVEXT
        DO JJ=IJB,IJE
          DO JI=IIB,IIE
            IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
            ZSTORE_3D(JI,JJ,JK) = ZTAU_TOTAL(IIJ,JBAND,JKRAD)
          END DO
        END DO
      END DO
      TZFIELD3D%CMNHNAME   = 'OTH_'//YBAND_NAME(JBAND)
      TZFIELD3D%CLONGNAME  = 'OTH_'//YBAND_NAME(JBAND)
      TZFIELD3D%CUNITS     = ''
      TZFIELD3D%CCOMMENT   = 'X_Y_Z_OTH_'//YBAND_NAME(JBAND)
      CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
      !
      DO JK=IKB,IKE
        JKRAD = JK - JPVEXT
        DO JJ=IJB,IJE
          DO JI=IIB,IIE
            IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
            ZSTORE_3D(JI,JJ,JK) = ZOMEGA_TOTAL(IIJ,JBAND,JKRAD)
          END DO
        END DO
      END DO
      TZFIELD3D%CMNHNAME   = 'SSA_'//YBAND_NAME(JBAND)
      TZFIELD3D%CLONGNAME  = 'SSA_'//YBAND_NAME(JBAND)
      TZFIELD3D%CUNITS     = ''
      TZFIELD3D%CCOMMENT   = 'X_Y_Z_SSA_'//YBAND_NAME(JBAND)
      CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
      !
      DO JK=IKB,IKE
        JKRAD = JK - JPVEXT
        DO JJ=IJB,IJE
          DO JI=IIB,IIE
            IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
            ZSTORE_3D(JI,JJ,JK) = ZCG_TOTAL(IIJ,JBAND,JKRAD)
          END DO
        END DO
      END DO
      TZFIELD3D%CMNHNAME   = 'ASF_'//YBAND_NAME(JBAND)
      TZFIELD3D%CLONGNAME  = 'ASF_'//YBAND_NAME(JBAND)
      TZFIELD3D%CUNITS     = ''
      TZFIELD3D%CCOMMENT   = 'X_Y_Z_ASF_'//YBAND_NAME(JBAND)
      CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
    END DO
  END IF
  !
  !
  IF (KRAD_DIAG >= 5)   THEN
!
! OZONE and AER optical thickness climato  entering the ecmwf_radiation_vers2
! note the vertical grid is re-inversed for graphic !   
    DO JK=IKB,IKE
      JKRAD = KFLEV+1 - JK + JPVEXT               
      DO JJ=IJB,IJE
        DO JI=IIB,IIE 
          IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
          ZSTORE_3D(JI,JJ,JK) = ZO3AVE(IIJ, JKRAD)
        END DO
      END DO
    END DO
    TZFIELD3D%CMNHNAME   = 'O3CLIM'
    TZFIELD3D%CLONGNAME  = 'O3CLIM'
    TZFIELD3D%CUNITS     = 'Pa Pa-1'
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_O3'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D)
! 
!cumulated optical thickness of aerosols
!cumul begin from the top of the domain, not from the TOA !      
!
!land 
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZSTORE_3D(JI,JJ,JK) = PAER(JI,JJ,JKRAD,1)
        END DO
      END DO
    END DO
!
    ZSTORE_2D (:,:) = 0.
    DO JK=IKB,IKE
      JK1=IKE-JK+IKB 
      ZSTORE_2D(:,:) = ZSTORE_2D(:,:) + ZSTORE_3D(:,:,JK1)
      ZSTORE_3D2(:,:,JK1) = ZSTORE_2D(:,:)  
    END DO
    TZFIELD3D%CMNHNAME   = 'CUM_AER_LAND'
    TZFIELD3D%CLONGNAME  = 'CUM_AER_LAND'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_CUM_AER_OPT'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D2)
!
! sea
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZSTORE_3D(JI,JJ,JK) = PAER(JI,JJ,JKRAD,2)
        END DO
      END DO
    END DO
!sum
    ZSTORE_2D (:,:) = 0.
    DO JK=IKB,IKE
      JK1=IKE-JK+IKB 
      ZSTORE_2D(:,:) = ZSTORE_2D(:,:) + ZSTORE_3D(:,:,JK1)
      ZSTORE_3D2(:,:,JK1) = ZSTORE_2D(:,:)  
    END DO
!
    TZFIELD3D%CMNHNAME   = 'CUM_AER_SEA'
    TZFIELD3D%CLONGNAME  = 'CUM_AER_SEA'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_CUM_AER_OPT'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D2)
!
! desert
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZSTORE_3D(JI,JJ,JK) = PAER(JI,JJ,JKRAD,3)
        END DO
      END DO
    END DO
!sum     
    ZSTORE_2D (:,:) = 0.
    DO JK=IKB,IKE
      JK1=IKE-JK+IKB 
      ZSTORE_2D(:,:) = ZSTORE_2D(:,:) + ZSTORE_3D(:,:,JK1)
      ZSTORE_3D2(:,:,JK1) = ZSTORE_2D(:,:)  
    END DO
!    
    TZFIELD3D%CMNHNAME   = 'CUM_AER_DES'
    TZFIELD3D%CLONGNAME  = 'CUM_AER_DES'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_CUM_AER_OPT'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D2)
!
! urban
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZSTORE_3D(JI,JJ,JK) = PAER(JI,JJ,JKRAD,4)
        END DO
      END DO
    END DO
!sum      
    ZSTORE_2D (:,:) = 0.
    DO JK=IKB,IKE
      JK1=IKE-JK+IKB 
      ZSTORE_2D(:,:) = ZSTORE_2D(:,:) + ZSTORE_3D(:,:,JK1)
      ZSTORE_3D2(:,:,JK1) = ZSTORE_2D(:,:)  
    END DO
!
    TZFIELD3D%CMNHNAME   = 'CUM_AER_URB'
    TZFIELD3D%CLONGNAME  = 'CUM_AER_URB'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_CUM_AER_OPT'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D2)
!
! Volcanoes
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZSTORE_3D(JI,JJ,JK) = PAER(JI,JJ,JKRAD,5)
        END DO
      END DO
    END DO
!sum         
    ZSTORE_2D (:,:) = 0.
    DO JK=IKB,IKE
      JK1=IKE-JK+IKB 
      ZSTORE_2D(:,:) = ZSTORE_2D(:,:) + ZSTORE_3D(:,:,JK1)
      ZSTORE_3D2(:,:,JK1) = ZSTORE_2D(:,:)  
    END DO
!
    TZFIELD3D%CMNHNAME   = 'CUM_AER_VOL'
    TZFIELD3D%CLONGNAME  = 'CUM_AER_VOL'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_CUM_AER_OPT'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D2)
!
! stratospheric background
    DO JK=IKB,IKE
      JKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZSTORE_3D(JI,JJ,JK) = PAER(JI,JJ,JKRAD,6)
        END DO
      END DO
    END DO
!sum      
    ZSTORE_2D (:,:) = 0.
    DO JK=IKB,IKE
      JK1=IKE-JK+IKB 
      ZSTORE_2D(:,:) = ZSTORE_2D(:,:) + ZSTORE_3D(:,:,JK1)
      ZSTORE_3D2(:,:,JK1) = ZSTORE_2D(:,:)  
    END DO
!
    TZFIELD3D%CMNHNAME   = 'CUM_AER_STRB'
    TZFIELD3D%CLONGNAME  = 'CUM_AER_STRB'
    TZFIELD3D%CUNITS     = ''
    TZFIELD3D%CCOMMENT   = 'X_Y_Z_CUM_AER_OPT'
    CALL IO_Field_write(TPFILE,TZFIELD3D,ZSTORE_3D2)
  ENDIF
END IF
!
DEALLOCATE(ZNFLW_CS)
DEALLOCATE(ZNFLW)
DEALLOCATE(ZNFSW_CS)
DEALLOCATE(ZNFSW)
DEALLOCATE(ZFLUX_TOP_GND_IRVISNIR)
DEALLOCATE(ZFLUX_SW_DOWN)
DEALLOCATE(ZFLUX_SW_UP)
DEALLOCATE(ZFLUX_LW)
DEALLOCATE(ZDTLW_CS)
DEALLOCATE(ZDTSW_CS)
DEALLOCATE(ZFLUX_TOP_GND_IRVISNIR_CS)
DEALLOCATE(ZPLAN_ALB_VIS)
DEALLOCATE(ZPLAN_ALB_NIR)
DEALLOCATE(ZPLAN_TRA_VIS)
DEALLOCATE(ZPLAN_TRA_NIR)
DEALLOCATE(ZPLAN_ABS_VIS)
DEALLOCATE(ZPLAN_ABS_NIR)
DEALLOCATE(ZEFCL_LWD)
DEALLOCATE(ZEFCL_LWU)
DEALLOCATE(ZFLWP)
DEALLOCATE(ZFIWP)
DEALLOCATE(ZRADLP)
DEALLOCATE(ZRADIP)
DEALLOCATE(ZEFCL_RRTM)
DEALLOCATE(ZCLSW_TOTAL)
DEALLOCATE(ZTAU_TOTAL)
DEALLOCATE(ZOMEGA_TOTAL)
DEALLOCATE(ZCG_TOTAL)
DEALLOCATE(ZFLUX_SW_DOWN_CS)
DEALLOCATE(ZFLUX_SW_UP_CS)
DEALLOCATE(ZFLUX_LW_CS)
DEALLOCATE(ZO3AVE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE RADIATIONS
!
END MODULE MODI_RADIATIONS  
