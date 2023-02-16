!MNH_LIC Copyright 2000-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ##########################
MODULE MODI_AIRCRAFT_BALLOON_EVOL
!      ##########################
!
INTERFACE
!
      SUBROUTINE AIRCRAFT_BALLOON_EVOL(PTSTEP,               &
                       PXHAT, PYHAT, PZ,                     &
                       PMAP, PLONOR, PLATOR,                 &
                       PU, PV, PW, PP, PTH, PR, PSV, PTKE,   &
                       PTS, PRHODREF, PCIT,TPFLYER, PSEA     )
!
USE MODD_AIRCRAFT_BALLOON
!
REAL,                     INTENT(IN)     :: PTSTEP ! time step
REAL, DIMENSION(:),       INTENT(IN)     :: PXHAT  ! x coordinate
REAL, DIMENSION(:),       INTENT(IN)     :: PYHAT  ! y coordinate
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PZ     ! z array
REAL, DIMENSION(:,:),     INTENT(IN)     :: PMAP   ! map factor
REAL,                     INTENT(IN)     :: PLONOR ! origine longitude
REAL,                     INTENT(IN)     :: PLATOR ! origine latitude
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PU     ! horizontal wind X component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PV     ! horizontal wind Y component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PW     ! vertical wind
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PP     ! pressure
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTH    ! potential temperature
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PR     ! water mixing ratios
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PSV    ! Scalar variables
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTKE   ! turbulent kinetic energy
REAL, DIMENSION(:,:),     INTENT(IN)     :: PTS    ! surface temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PRHODREF ! dry air density of the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PCIT     ! pristine ice concentration
!
TYPE(FLYER),              INTENT(INOUT)  :: TPFLYER! balloon/aircraft
REAL, DIMENSION(:,:),OPTIONAL,INTENT(IN)     :: PSEA
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AIRCRAFT_BALLOON_EVOL
!
END INTERFACE
!
END MODULE MODI_AIRCRAFT_BALLOON_EVOL
!
!     ########################################################
      SUBROUTINE AIRCRAFT_BALLOON_EVOL(PTSTEP,               &
                       PXHAT, PYHAT, PZ,                     &
                       PMAP, PLONOR, PLATOR,                 &
                       PU, PV, PW, PP, PTH, PR, PSV, PTKE,   &
                       PTS, PRHODREF, PCIT,TPFLYER, PSEA     )
!     ########################################################
!
!
!!****  *AIRCRAFT_BALLOON_EVOL* - (advects and) stores 
!!                                balloons/aircrafts in the model
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!    
!! 1) All the balloons are tested. If the current balloon is
!!     a) in the current model
!!     b) not crashed
!!   the following computations are done.
!!
!! 2) The balloon position is computed.
!!       Interpolations at balloon positions are performed according to mass
!! points (because density is computed here for iso-density balloons).
!! Therefore, all model variables are used at mass points. Shuman averaging
!! are performed on X, Y, Z, U, V, W.
!!
!! 3) Storage of balloon data
!!       If storage is asked for this time-step, the data are recorded in the
!! balloon time-series.
!!
!! 4) Balloon advection
!!       If the balloon is launched, it is advected according its type
!!    a) iso-density balloons are advected following horizontal wind.
!!          the slope of the iso-density surfaces is neglected.
!!    b) radio-sounding balloons are advected according to all wind velocities.
!!          the vertical ascent speed is added to the vertical wind speed.  
!!    c) Constant Volume balloons are advected according to all wind velocities.
!!          the vertical ascent speed is computed using the balloon equation
!!  
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Valery Masson             * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!     Original 15/05/2000
!!     Apr,19, 2001 (G.Jaubert) add CVBALL type
!!     March, 2008 (P.Lacarrere) Add 3D fluxes
!!     Dec,12, 2008 (M. Leriche) move ZTDIST out from if.not.(tpflyer%fly)
!!     Dec,15, 2008 (V. Masson) correct do while aircraft move
!!     March, 2013 (O.Caumont) add radar reflectivities
!!     April, 2014 (C.Lac) allow RARE calculation only if CCLOUD=ICE3
!!     May, 2014 (O.Caumont) modify RARE for hydrometeors containing ice
!!                           add bright band calculation for RARE
!!     Feb, 2015 (C.Lac) Correction to prevent aircraft crash 
!!     July, 2015 (O.Nuissier/F.Duffourg) Add microphysics diagnostic for
!!                                      aircraft, ballon and profiler
!!      October, 2016 (G.DELAUTIER) LIMA
!!     March,28, 2018 (P. Wautelet) replace TEMPORAL_DIST by DATETIME_DISTANCE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 13/09/2019: budget: simplify and modernize date/time management
!  P. Wautelet 01/10/2020: bugfix: initialize GSTORE
!  P. Wautelet 14/01/2021: bugfixes: -ZXCOEF and ZYCOEF were not computed if CVBALL
!                                    -PCIT was used if CCLOUD/=ICEx (not allocated)
!                                    -PSEA was always used even if not allocated (CSURF/=EXTE)
!                                    -do not use PMAP if cartesian domain
!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_AIRCRAFT_BALLOON
USE MODD_CONF
USE MODD_CST
USE MODD_DIAG_IN_RUN
USE MODD_GRID
USE MODD_LUNIT_n,          ONLY: TLUOUT
USE MODD_NESTING
USE MODD_NSV,              ONLY : NSV_LIMA_NI,NSV_LIMA_NR,NSV_LIMA_NC
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA,       ONLY: XALPHAR_L=>XALPHAR,XNUR_L=>XNUR,XALPHAS_L=>XALPHAS,XNUS_L=>XNUS,&
                                 XALPHAG_L=>XALPHAG,XNUG_L=>XNUG, XALPHAI_L=>XALPHAI,XNUI_L=>XNUI,&
                                 XRTMIN_L=>XRTMIN,XALPHAC_L=>XALPHAC,XNUC_L=>XNUC,LSNOW_T_L=>LSNOW_T
USE MODD_PARAM_LIMA_COLD,  ONLY: XDI_L=>XDI,XLBEXI_L=>XLBEXI,XLBI_L=>XLBI,XAI_L=>XAI,XBI_L=>XBI,XC_I_L=>XC_I,&
                                 XLBEXS_L=>XLBEXS,XLBS_L=>XLBS,XCCS_L=>XCCS,&
                                 XAS_L=>XAS,XBS_L=>XBS,XCXS_L=>XCXS,XNS_L=>XNS,        &
                                 XLBDAS_MAX,XLBDAS_MIN
USE MODD_PARAM_LIMA_MIXED, ONLY: XDG_L=>XDG,XLBEXG_L=>XLBEXG,XLBG_L=>XLBG,XCCG_L=>XCCG,&
                                 XAG_L=>XAG,XBG_L=>XBG,XCXG_L=>XCXG,XCG_L=>XCG
USE MODD_PARAM_LIMA_WARM,  ONLY: XLBEXR_L=>XLBEXR,XLBR_L=>XLBR,XBR_L=>XBR,XAR_L=>XAR,&
                                 XBC_L=>XBC,XAC_L=>XAC
USE MODD_PARAM_n,          ONLY: CCLOUD, CSURF
USE MODD_PARAM_ICE,        ONLY: LSNOW_T_I=>LSNOW_T
USE MODD_RAIN_ICE_DESCR,   ONLY: XALPHAR_I=>XALPHAR,XNUR_I=>XNUR,XLBEXR_I=>XLBEXR,&
                                 XLBR_I=>XLBR,XCCR_I=>XCCR,XBR_I=>XBR,XAR_I=>XAR,&
                                 XALPHAC_I=>XALPHAC,XNUC_I=>XNUC,&
                                 XLBC_I=>XLBC,XBC_I=>XBC,XAC_I=>XAC,&
                                 XALPHAC2_I=>XALPHAC2,XNUC2_I=>XNUC2,&
                                 XALPHAS_I=>XALPHAS,XNUS_I=>XNUS,XLBEXS_I=>XLBEXS,&
                                 XLBS_I=>XLBS,XCCS_I=>XCCS,XAS_I=>XAS,XBS_I=>XBS,XCXS_I=>XCXS,XNS_I=>XNS,&
                                 XALPHAG_I=>XALPHAG,XNUG_I=>XNUG,XDG_I=>XDG,XLBEXG_I=>XLBEXG,&
                                 XLBG_I=>XLBG,XCCG_I=>XCCG,XAG_I=>XAG,XBG_I=>XBG,XCXG_I=>XCXG,XCG_I=>XCG,&
                                 XALPHAI_I=>XALPHAI,XNUI_I=>XNUI,XDI_I=>XDI,XLBEXI_I=>XLBEXI,&
                                 XLBI_I=>XLBI,XAI_I=>XAI,XBI_I=>XBI,XC_I_I=>XC_I,&
                                 XRTMIN_I=>XRTMIN,XCONC_LAND,XCONC_SEA
USE MODD_REF_n,            ONLY: XRHODREF
USE MODD_TIME,             only: tdtexp
USE MODD_TIME_n,           only: tdtcur
USE MODD_TURB_FLUX_AIRCRAFT_BALLOON
!
USE MODE_DATETIME
USE MODE_FGAU,             ONLY: GAULAG
USE MODE_FSCATTER,         ONLY: QEPSW,QEPSI,BHMIE,MOMG,MG
USE MODE_GRIDPROJ
USE MODE_ll
USE MODE_MSG
!
USE MODI_GAMMA,            ONLY: GAMMA
USE MODI_WATER_SUM
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
REAL,                     INTENT(IN)     :: PTSTEP ! time step
REAL, DIMENSION(:),       INTENT(IN)     :: PXHAT  ! x coordinate
REAL, DIMENSION(:),       INTENT(IN)     :: PYHAT  ! y coordinate
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PZ     ! z array
REAL, DIMENSION(:,:),     INTENT(IN)     :: PMAP   ! map factor
REAL,                     INTENT(IN)     :: PLONOR ! origine longitude
REAL,                     INTENT(IN)     :: PLATOR ! origine latitude
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PU     ! horizontal wind X component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PV     ! horizontal wind Y component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PW     ! vertical wind
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PP     ! pressure
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTH    ! potential temperature
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PR     ! water mixing ratios
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PSV    ! Scalar variables
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTKE   ! turbulent kinetic energy
REAL, DIMENSION(:,:),     INTENT(IN)     :: PTS    ! surface temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PRHODREF ! dry air density of the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PCIT     ! pristine ice concentration
!
TYPE(FLYER),              INTENT(INOUT)  :: TPFLYER! balloon/aircraft
REAL, DIMENSION(:,:),     INTENT(IN)     :: PSEA
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!
INTEGER :: IMI        ! model index
REAL    :: ZTHIS_PROC ! 1 if balloon is currently treated by this proc., else 0
!
INTEGER :: IIB        ! current processor domain sizes
INTEGER :: IJB
INTEGER :: IIE
INTEGER :: IJE
INTEGER :: IIU
INTEGER :: IJU
INTEGER :: IKB
INTEGER :: IKE
INTEGER :: IKU
!
INTEGER :: JK         ! loop index
!
REAL, DIMENSION(SIZE(PXHAT))        :: ZXHATM ! mass point coordinates
REAL, DIMENSION(SIZE(PYHAT))        :: ZYHATM ! mass point coordinates
!
REAL, DIMENSION(2,2,SIZE(PZ,3))     :: ZZM    ! mass point coordinates
REAL, DIMENSION(2,2,SIZE(PZ,3))     :: ZZU    ! U points z coordinates
REAL, DIMENSION(2,2,SIZE(PZ,3))     :: ZZV    ! V points z coordinates
REAL, DIMENSION(2,2,SIZE(PZ,3))     :: ZWM    ! mass point wind
!
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZTHV   ! virtual potential temperature
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZTV    ! virtual temperature
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZTEMP  ! temperature
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZEXN   ! Exner function
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZRHO   ! air density
REAL                                :: ZFLYER_EXN ! balloon/aircraft Exner func.
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZTHW_FLUX  !       
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZRCW_FLUX  !
REAL, DIMENSION(2,2,SIZE(PSV,3),SIZE(PSV,4))    :: ZSVW_FLUX
!
REAL    :: ZTDIST   ! time until launch (sec)
LOGICAL :: GLAUNCH  ! launch/takeoff is effective at this time-step (if true)
LOGICAL :: GSTORE   ! storage occurs at this time step
!
INTEGER :: II       ! mass balloon position (x index)
INTEGER :: IJ       ! mass balloon position (y index)
INTEGER :: IU       ! U flux point balloon position (x index)
INTEGER :: IV       ! V flux point balloon position (y index)
INTEGER :: IDU      ! difference between IU and II
INTEGER :: IDV      ! difference between IV and IJ
!
INTEGER :: IK00     ! balloon position for II  , IJ
INTEGER :: IK01     ! balloon position for II  , IJ+1
INTEGER :: IK10     ! balloon position for II+1, IJ
INTEGER :: IK11     ! balloon position for II+1, IJ+1
INTEGER :: IU00     ! balloon position for IU  , IJ
INTEGER :: IU01     ! balloon position for IU  , IJ+1
INTEGER :: IU10     ! balloon position for IU+1, IJ
INTEGER :: IU11     ! balloon position for IU+1, IJ+1
INTEGER :: IV00     ! balloon position for II  , IV
INTEGER :: IV01     ! balloon position for II  , IV+1
INTEGER :: IV10     ! balloon position for II+1, IV
INTEGER :: IV11     ! balloon position for II+1, IV+1
!
REAL :: ZXCOEF      ! X direction interpolation coefficient
REAL :: ZUCOEF      ! X direction interpolation coefficient (for U)
REAL :: ZYCOEF      ! Y direction interpolation coefficient
REAL :: ZVCOEF      ! Y direction interpolation coefficient (for V)
!
REAL :: ZZCOEF00    ! Z direction interpolation coefficient for II  , IJ
REAL :: ZZCOEF01    ! Z direction interpolation coefficient for II  , IJ+1
REAL :: ZZCOEF10    ! Z direction interpolation coefficient for II+1, IJ
REAL :: ZZCOEF11    ! Z direction interpolation coefficient for II+1, IJ+1
REAL :: ZUCOEF00    ! Z direction interpolation coefficient for IU  , IJ
REAL :: ZUCOEF01    ! Z direction interpolation coefficient for IU  , IJ+1
REAL :: ZUCOEF10    ! Z direction interpolation coefficient for IU+1, IJ
REAL :: ZUCOEF11    ! Z direction interpolation coefficient for IU+1, IJ+1
REAL :: ZVCOEF00    ! Z direction interpolation coefficient for II  , IV
REAL :: ZVCOEF01    ! Z direction interpolation coefficient for II  , IV+1
REAL :: ZVCOEF10    ! Z direction interpolation coefficient for II+1, IV
REAL :: ZVCOEF11    ! Z direction interpolation coefficient for II+1, IV+1
!
INTEGER :: IN       ! time index
INTEGER :: JLOOP,JLOOP2    ! loop counter
!
REAL    :: ZU_BAL   ! horizontal wind speed at balloon location (along x)
REAL    :: ZV_BAL   ! horizontal wind speed at balloon location (along y)
REAL    :: ZW_BAL   ! vertical   wind speed at balloon location (along z)
REAL    :: ZMAP     ! map factor at balloon location
REAL    :: ZGAM     ! rotation between meso-nh base and spherical lat-lon base.
INTEGER :: IL       ! flight segment index
REAL    :: ZSEG_FRAC! fraction of flight in the current segment
REAL    :: ZRO_BAL  ! air density at balloon location
!
INTEGER :: IINFO_ll ! return code
INTEGER :: ILUOUT   ! logical unit
INTEGER :: IRESP    ! return code
!
! specific to cloud radar
REAL, DIMENSION(SIZE(PR,3))    :: ZTEMPZ! vertical profile of temperature
REAL, DIMENSION(SIZE(PR,3))    :: ZRHODREFZ ! vertical profile of dry air density of the reference state
REAL, DIMENSION(SIZE(PR,3))    :: ZCIT     ! pristine ice concentration
REAL, DIMENSION(SIZE(PR,3))    :: ZCCI,ZCCR,ZCCC     ! ICE,RAIN CLOUD concentration (LIMA)
REAL, DIMENSION(SIZE(PR,1),SIZE(PR,2),SIZE(PR,3))    :: ZR   
REAL, DIMENSION(SIZE(PR,3),SIZE(PR,4)+1) :: ZRZ  ! vertical profile of hydrometeor mixing ratios
REAL                           :: ZA,ZB,ZCC,ZCX,ZALPHA,ZNU,ZLB,ZLBEX,ZNS,ZRHOHYD   ! generic microphysical parameters
INTEGER                        :: JJ    ! loop counter for quadrature
COMPLEX                        :: QMW,QMI,QM,QB,QEPSIW,QEPSWI   ! dielectric parameter
REAL                           :: ZAETOT,ZAETMP,ZREFLOC,ZQSCA,ZQBACK,ZQEXT ! temporary scattering parameters
REAL,DIMENSION(:),ALLOCATABLE  :: ZAELOC,ZZMZ ! temporary arrays
INTEGER                        :: JPTS_GAULAG=7 ! number of points for Gauss-Laguerre quadrature
REAL                           :: ZLBDA   ! slope distribution parameter
REAL                           :: ZN   ! number concentration
REAL                           :: ZFRAC_ICE  ! ice water fraction
REAL                           :: ZDELTA_EQUIV ! mass-equivalent Gauss-Laguerre point
REAL                           :: ZFW ! liquid fraction
REAL                           :: ZFPW ! weight for mixed-phase reflectivity
REAL,DIMENSION(:),ALLOCATABLE  :: ZX,ZW ! Gauss-Laguerre points and weights
REAL,DIMENSION(:),ALLOCATABLE  :: ZRTMIN ! local values for XRTMIN
LOGICAL                        :: GCALC
!----------------------------------------------------------------------------
!
!*      1.   PRELIMINARIES
!            -------------
!
IF(.NOT. ALLOCATED(XTHW_FLUX)) &
ALLOCATE(XTHW_FLUX(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3)))
IF(.NOT. ALLOCATED(XRCW_FLUX)) &
ALLOCATE(XRCW_FLUX(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3)))
IF(.NOT. ALLOCATED(XSVW_FLUX)) &
ALLOCATE(XSVW_FLUX(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3),SIZE(PSV,4)))
ILUOUT = TLUOUT%NLU
!
ZR = 0.
GSTORE = .FALSE.
!
!*      1.0  initialization of processor test
!            --------------------------------
!
ZTHIS_PROC=0.
!
!
!*      1.1  test on the model
!            -----------------
!
CALL GET_MODEL_NUMBER_ll  (IMI)
!
!
IF (TPFLYER%MODEL  /= 'FIX' .AND. COUNT(NDAD(:) == IMI) /= 0 &
   .AND. ( TPFLYER%NMODEL == IMI .OR. NDAD(TPFLYER%NMODEL) == IMI ) &
   .AND. TPFLYER%X_CUR /= XUNDEF .AND. TPFLYER%Y_CUR /= XUNDEF &
   .AND.  TPFLYER%FLY .AND. .NOT. TPFLYER%CRASH &
   .AND. CPROGRAM == 'MESONH' ) THEN
  CALL FLYER_CHANGE_MODEL(IMI)
ENDIF
!
IF ( TPFLYER%NMODEL /= IMI ) RETURN
!
!----------------------------------------------------------------------------
!
!*      2.   PRELIMINARIES-2
!            -------------
!
!*      2.1  Indices
!            -------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB =   1   + JPVEXT
IKE = SIZE(PZ,3) - JPVEXT
IKU = SIZE(PZ,3)
!
!
!*      2.2  Interpolations of model variables to mass points
!            ------------------------------------------------
!
IIU=SIZE(PXHAT)
IJU=SIZE(PYHAT)
!
ZXHATM(1:IIU-1)=0.5*PXHAT(1:IIU-1)+0.5*PXHAT(2:IIU  )
ZXHATM(  IIU  )=1.5*PXHAT(  IIU  )-0.5*PXHAT(  IIU-1)
!
ZYHATM(1:IJU-1)=0.5*PYHAT(1:IJU-1)+0.5*PYHAT(2:IJU  )
ZYHATM(  IJU  )=1.5*PYHAT(  IJU  )-0.5*PYHAT(  IJU-1)
!----------------------------------------------------------------------------
!
!*      2.3  Compute time until launch by comparison of dates and times
!            ----------------------------------------------------------
!
CALL DATETIME_DISTANCE(TPFLYER%LAUNCH,TDTCUR,ZTDIST)
!
!*      3.   LAUNCH
!            ------
!
GLAUNCH     = .FALSE.
!
!
IF (.NOT. TPFLYER%FLY) THEN
!
!
!*      3.1  comparison of dates and times
!            -----------------------------
!
!  CALL DATETIME_DISTANCE(TPFLYER%LAUNCH,TDTCUR,ZTDIST)
!
!*      3.2  launch/takeoff is effective
!            ---------------------------
!
  IF (ZTDIST >= - PTSTEP ) THEN
    IF (TPFLYER%TYPE=='AIRCRA') THEN
!
!*     3.2.1 Determination of flight segment
!            -------------------------------
!
      TPFLYER%SEGCURN = 1
      IL = TPFLYER%SEGCURN
      !
      TPFLYER%SEGCURT = ZTDIST
      !
      DO WHILE (TPFLYER%SEGCURT>TPFLYER%SEGTIME(IL) .AND. IL <= TPFLYER%SEG)
        TPFLYER%SEGCURN = TPFLYER%SEGCURN + 1
        IL = TPFLYER%SEGCURN
        TPFLYER%SEGCURT = TPFLYER%SEGCURT - TPFLYER%SEGTIME(IL-1)
        IF (IL>TPFLYER%SEG) EXIT
      END DO
      !
      !* end of flight
      !
      IF (IL > TPFLYER%SEG) THEN
        TPFLYER%FLY=.FALSE.
      ELSE
        TPFLYER%FLY = .TRUE.
        GLAUNCH     = .TRUE.
        TPFLYER%CRASH=.FALSE.
        IF (ZTDIST <= PTSTEP ) THEN
          WRITE(ILUOUT,*) '-------------------------------------------------------------------'
          WRITE(ILUOUT,*) 'Aircraft ',TPFLYER%TITLE,' takes off the   ', &
                      TDTCUR%nday,'/',TDTCUR%nmonth,'/',                 &
                      TDTCUR%nyear,' at ',NINT(TDTCUR%xtime),' sec.'
          WRITE(ILUOUT,*) '-------------------------------------------------------------------'
        ENDIF
      ENDIF
    ELSE IF (ZTDIST <= PTSTEP ) THEN
      TPFLYER%FLY = .TRUE.
      GLAUNCH     = .TRUE.
      WRITE(ILUOUT,*) '-------------------------------------------------------------------'
      WRITE(ILUOUT,*) 'Balloon  ',TPFLYER%TITLE,' is launched the ', &
                    TDTCUR%nday,'/',TDTCUR%nmonth,'/',               &
                    TDTCUR%nyear,' at ',NINT(TDTCUR%xtime),' sec.'
      WRITE(ILUOUT,*) '-------------------------------------------------------------------'
    END IF
!
!*      3.3  Initial horizontal positions
!            ----------------------------
!
    IF (TPFLYER%TYPE=='RADIOS' .OR. TPFLYER%TYPE=='ISODEN' .OR. TPFLYER%TYPE=='CVBALL') THEN
      TPFLYER%X_CUR = TPFLYER%XLAUNCH
      TPFLYER%Y_CUR = TPFLYER%YLAUNCH
    END IF
    IF (TPFLYER%TYPE=='AIRCRA') THEN
!
!
!*       3.3.2 Determination of initial position
!              -----------------------------
!
      IF (TPFLYER%FLY) THEN
        ZSEG_FRAC = TPFLYER%SEGCURT / TPFLYER%SEGTIME(IL)
        !
        TPFLYER%X_CUR = (1.-ZSEG_FRAC) * TPFLYER%SEGX(IL  ) &
                      +     ZSEG_FRAC  * TPFLYER%SEGX(IL+1)
        TPFLYER%Y_CUR = (1.-ZSEG_FRAC) * TPFLYER%SEGY(IL  ) &
                      +     ZSEG_FRAC  * TPFLYER%SEGY(IL+1)
      END IF
!
    END IF
  END IF
END IF
!
!*      3.4  instant of storage
!            ------------------
!
IF ( TPFLYER%T_CUR == XUNDEF ) TPFLYER%T_CUR = TPFLYER%STEP - PTSTEP
!
TPFLYER%T_CUR = TPFLYER%T_CUR + PTSTEP
!
IF ( TPFLYER%T_CUR >= TPFLYER%STEP - 1.E-10 ) THEN
  GSTORE = .TRUE.
  TPFLYER%T_CUR = TPFLYER%T_CUR - TPFLYER%STEP
  TPFLYER%N_CUR = TPFLYER%N_CUR + 1
END IF
!
IF (GSTORE) THEN
  IN = TPFLYER%N_CUR
#if 0
  tpflyer%tpdates(in)%nyear  = tdtexp%nyear
  tpflyer%tpdates(in)%nmonth = tdtexp%nmonth
  tpflyer%tpdates(in)%nday   = tdtexp%nday
  tpflyer%tpdates(in)%xtime  = tdtexp%xtime + ( in - 1 ) * tpflyer%step
#else
  tpflyer%tpdates(in) = tdtcur
#endif
END IF
!
IF ( TPFLYER%FLY) THEN
!
!----------------------------------------------------------------------------
!
!*      4.   FLYER POSITION
!            --------------
!
!*      4.1  X position
!            ----------
!
  IU=COUNT( PXHAT (:)<=TPFLYER%X_CUR )
  II=COUNT( ZXHATM(:)<=TPFLYER%X_CUR )
!
  IF (IU<IIB   .AND. LWEST_ll()) THEN
    IF (TPFLYER%MODEL == 'FIX' .OR. TPFLYER%NMODEL == 1 ) THEN 
      TPFLYER%CRASH=.TRUE.
    ELSE
      II=IIB
      IU=IIB
    END IF
  END IF
  IF (IU>IIE     .AND. LEAST_ll())  THEN
    IF (TPFLYER%MODEL == 'FIX'  .OR. TPFLYER%NMODEL == 1) THEN 
      TPFLYER%CRASH=.TRUE.
    ELSE
      II=IIE
      IU=IIE
    END IF
  END IF
!
!
!*      4.2  Y position
!            ----------
!
  IV=COUNT( PYHAT (:)<=TPFLYER%Y_CUR )
  IJ=COUNT( ZYHATM(:)<=TPFLYER%Y_CUR )
!
  IF (IV<IJB   .AND. LSOUTH_ll()) THEN
    IF (TPFLYER%MODEL == 'FIX'  .OR. TPFLYER%NMODEL == 1) THEN 
      TPFLYER%CRASH=.TRUE.
    ELSE
      IJ=IJB
      IV=IJB
    END IF
  END IF
  IF (IV>IJE     .AND. LNORTH_ll()) THEN
    IF (TPFLYER%MODEL == 'FIX'  .OR. TPFLYER%NMODEL == 1) THEN 
      TPFLYER%CRASH=.TRUE.
    ELSE
      IJ=IJE
      IV=IJE
    END IF
  END IF
!
!
!*      4.3  Position of balloon according to processors
!            -------------------------------------------
!
  IF (IU>=IIB .AND. IU<=IIE .AND. IV>=IJB .AND. IV<=IJE) ZTHIS_PROC=1.
!
!
!*      4.4  Computations only on correct processor
!            --------------------------------------
!
!----------------------------------------------------------------------------
  IF (ZTHIS_PROC>0. .AND. .NOT. TPFLYER%CRASH) THEN
!----------------------------------------------------------------------------
!
!*      4.5  Interpolations of model variables to mass points
!            ------------------------------------------------
!

    ZZM(:,:,1:IKU-1)=0.5 *PZ(II  :II+1,IJ  :IJ+1,1:IKU-1)+0.5 *PZ(II  :II+1,IJ  :IJ+1,2:IKU  )
    ZZM(:,:,  IKU  )=1.5 *PZ(II  :II+1,IJ  :IJ+1,  IKU-1)-0.5 *PZ(II  :II+1,IJ  :IJ+1,  IKU-2)
!
    IDU = IU - II
    ZZU(:,:,1:IKU-1)=0.25*PZ(IDU+II-1:IDU+II,  IJ  :IJ+1,1:IKU-1)+0.25*PZ(IDU+II-1:IDU+II  ,IJ  :IJ+1,2:IKU  ) &
                  +0.25*PZ(IDU+II  :IDU+II+1,IJ  :IJ+1,1:IKU-1)+0.25*PZ(IDU+II  :IDU+II+1,IJ  :IJ+1,2:IKU  )
    ZZU(:,:,  IKU  )=0.75*PZ(IDU+II-1:IDU+II  ,IJ  :IJ+1,  IKU-1)-0.25*PZ(IDU+II-1:IDU+II  ,IJ  :IJ+1,  IKU-2) &
                  +0.75*PZ(IDU+II  :IDU+II+1,IJ  :IJ+1,  IKU-1)-0.25*PZ(IDU+II  :IDU+II+1,IJ  :IJ+1,  IKU-2)

    IDV = IV - IJ 
    ZZV(:,:,1:IKU-1)=0.25*PZ(II  :II+1,IDV+IJ-1:IDV+IJ  ,1:IKU-1)+0.25*PZ(II  :II+1,IDV+IJ-1:IDV+IJ  ,2:IKU  ) &
                  +0.25*PZ(II  :II+1,IDV+IJ  :IDV+IJ+1,1:IKU-1)+0.25*PZ(II  :II+1,IDV+IJ  :IDV+IJ+1,2:IKU  )
    ZZV(:,:,  IKU  )=0.75*PZ(II  :II+1,IDV+IJ-1:IDV+IJ  ,  IKU-1)-0.25*PZ(II  :II+1,IDV+IJ-1:IDV+IJ  ,  IKU-2) &
                  +0.75*PZ(II  :II+1,IDV+IJ  :IDV+IJ+1,  IKU-1)-0.25*PZ(II  :II+1,IDV+IJ  :IDV+IJ+1,  IKU-2)
!
!
    ZWM(:,:,1:IKU-1)=0.5*PW(II:II+1,IJ:IJ+1,1:IKU-1)+0.5*PW(II:II+1,IJ:IJ+1,2:IKU  )
    ZWM(:,:,  IKU  )=1.5*PW(II:II+1,IJ:IJ+1,  IKU-1)-0.5*PW(II:II+1,IJ:IJ+1,  IKU-2)
!
!----------------------------------------------------------------------------
!
!*      5.   BALLOON/AIRCRAFT VERTICAL POSITION
!            ----------------------------------
!
!
!*      5.1  Density
!            -------
!
    ZEXN(:,:,:    ) = (PP(II:II+1,IJ:IJ+1,:)/XP00)**(XRD/XCPD)
    DO JK=IKB-1,1,-1
      ZEXN(:,:,JK) = 1.5 * ZEXN(:,:,JK+1) - 0.5 * ZEXN(:,:,JK+2)
    END DO
    DO JK=IKE+1,IKU
      ZEXN(:,:,JK) = 1.5 * ZEXN(:,:,JK-1) - 0.5 * ZEXN(:,:,JK-2)
    END DO
    !
    IF (TPFLYER%TYPE=='ISODEN' .OR. TPFLYER%TYPE=='CVBALL' &
        .OR. TPFLYER%TYPE=='AIRCRA' ) THEN
      ZTHV(:,:,:) = PTH(II:II+1,IJ:IJ+1,:)
      IF (SIZE(PR,4)>0)                                                     &
      ZTHV(:,:,:) = ZTHV(:,:,:) * ( 1. + XRV/XRD*PR(II:II+1,IJ:IJ+1,:,1) )  &
                                / ( 1. + WATER_SUM(PR(II:II+1,IJ:IJ+1,:,:)) )
      !
      ZTV (:,:,:) = ZTHV(:,:,:) * ZEXN(:,:,:)
      ZRHO(:,:,:) = PP(II:II+1,IJ:IJ+1,:) / (XRD*ZTV(:,:,:))
      DO JK=IKB-1,1,-1
        ZRHO(:,:,JK) = 1.5 * ZRHO(:,:,JK+1) - 0.5 * ZRHO(:,:,JK+2)
      END DO
      DO JK=IKE+1,IKU
        ZRHO(:,:,JK) = 1.5 * ZRHO(:,:,JK-1) - 0.5 * ZRHO(:,:,JK-2)
      END DO
     ZTHW_FLUX(:,:,:) = ZRHO(:,:,:)*XCPD *XTHW_FLUX(II:II+1,IJ:IJ+1,:)
     ZRCW_FLUX(:,:,:) = ZRHO(:,:,:)*XLVTT*XRCW_FLUX(II:II+1,IJ:IJ+1,:)
     ZSVW_FLUX(:,:,:,:) = XSVW_FLUX(II:II+1,IJ:IJ+1,:,:)
    END IF

!
!*      5.2  Initial vertical positions
!            --------------------------
!
    IF (GLAUNCH) THEN
!
!*      5.2.1 Iso-density balloon
!
      IF (TPFLYER%TYPE=='ISODEN') THEN
        ZXCOEF = (TPFLYER%X_CUR - ZXHATM(II)) / (ZXHATM(II+1) - ZXHATM(II))
        ZXCOEF = MAX (0.,MIN(ZXCOEF,1.))
        ZYCOEF = (TPFLYER%Y_CUR - ZYHATM(IJ)) / (ZYHATM(IJ+1) - ZYHATM(IJ))
        ZYCOEF = MAX (0.,MIN(ZYCOEF,1.))
        IF ( TPFLYER%ALT /= XUNDEF ) THEN
          IK00 = MAX ( COUNT (TPFLYER%ALT >= ZZM(1,1,:)), 1)
          IK01 = MAX ( COUNT (TPFLYER%ALT >= ZZM(1,2,:)), 1)
          IK10 = MAX ( COUNT (TPFLYER%ALT >= ZZM(2,1,:)), 1)
          IK11 = MAX ( COUNT (TPFLYER%ALT >= ZZM(2,2,:)), 1)
          ZZCOEF00 = (TPFLYER%ALT - ZZM(1,1,IK00)) / ( ZZM(1,1,IK00+1) - ZZM(1,1,IK00))
          ZZCOEF01 = (TPFLYER%ALT - ZZM(1,2,IK01)) / ( ZZM(1,2,IK01+1) - ZZM(1,2,IK01))
          ZZCOEF10 = (TPFLYER%ALT - ZZM(2,1,IK10)) / ( ZZM(2,1,IK10+1) - ZZM(2,1,IK10))
          ZZCOEF11 = (TPFLYER%ALT - ZZM(2,2,IK11)) / ( ZZM(2,2,IK11+1) - ZZM(2,2,IK11))
          TPFLYER%RHO = FLYER_INTERP(ZRHO)
        ELSE IF ( TPFLYER%PRES /= XUNDEF ) THEN
          ZFLYER_EXN = (TPFLYER%PRES/XP00)**(XRD/XCPD)
          IK00 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(1,1,:)), 1)
          IK01 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(1,2,:)), 1)
          IK10 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(2,1,:)), 1)
          IK11 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(2,2,:)), 1)
          ZZCOEF00 = (ZFLYER_EXN - ZEXN(1,1,IK00)) / ( ZEXN(1,1,IK00+1) - ZEXN(1,1,IK00))
          ZZCOEF01 = (ZFLYER_EXN - ZEXN(1,2,IK01)) / ( ZEXN(1,2,IK01+1) - ZEXN(1,2,IK01))
          ZZCOEF10 = (ZFLYER_EXN - ZEXN(2,1,IK10)) / ( ZEXN(2,1,IK10+1) - ZEXN(2,1,IK10))
          ZZCOEF11 = (ZFLYER_EXN - ZEXN(2,2,IK11)) / ( ZEXN(2,2,IK11+1) - ZEXN(2,2,IK11))
          TPFLYER%RHO = FLYER_INTERP(ZRHO)
        ELSE
          WRITE(ILUOUT,*) 'Error in balloon initial position (balloon ',TPFLYER%TITLE,' )'
          WRITE(ILUOUT,*) 'neither initial ALTITUDE or PRESsure is given'
          WRITE(ILUOUT,*) 'Check your INI_BALLOON routine'
!callabortstop
          CALL PRINT_MSG(NVERB_FATAL,'GEN','AIRCRAFT_BALLOON_EVOL','')
        END IF
      END IF
!
!*      5.2.2 Radiosounding balloon
!
      IF (TPFLYER%TYPE=='RADIOS') THEN
        TPFLYER%Z_CUR = TPFLYER%ALT
        TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(1,1,IKB) )
        TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(2,1,IKB) )
        TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(1,2,IKB) )
        TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(2,2,IKB) )
      END IF

!*      5.2.3 Aircraft
!
      IF (TPFLYER%TYPE=='AIRCRA') THEN
       IF (TPFLYER%ALTDEF) THEN
         TPFLYER%P_CUR = (1.-ZSEG_FRAC) * TPFLYER%SEGP(IL  ) &
                      +     ZSEG_FRAC  * TPFLYER%SEGP(IL+1)
       ELSE 
         TPFLYER%Z_CUR = (1.-ZSEG_FRAC) * TPFLYER%SEGZ(IL ) &
                      +     ZSEG_FRAC  * TPFLYER%SEGZ(IL +1 )
       END IF
      END IF
!
!*      5.2.4 Constant Volume Balloon
!
      IF (TPFLYER%TYPE=='CVBALL') THEN
        ZXCOEF = (TPFLYER%X_CUR - ZXHATM(II)) / (ZXHATM(II+1) - ZXHATM(II))
        ZXCOEF = MAX (0.,MIN(ZXCOEF,1.))
        ZYCOEF = (TPFLYER%Y_CUR - ZYHATM(IJ)) / (ZYHATM(IJ+1) - ZYHATM(IJ))
        ZYCOEF = MAX (0.,MIN(ZYCOEF,1.))
        IF ( TPFLYER%ALT /= XUNDEF ) THEN
          IK00 = MAX ( COUNT (TPFLYER%ALT >= ZZM(1,1,:)), 1)
          IK01 = MAX ( COUNT (TPFLYER%ALT >= ZZM(1,2,:)), 1)
          IK10 = MAX ( COUNT (TPFLYER%ALT >= ZZM(2,1,:)), 1)
          IK11 = MAX ( COUNT (TPFLYER%ALT >= ZZM(2,2,:)), 1)
          IF (IK00*IK01*IK10*IK11 .EQ. 0) THEN
            TPFLYER%Z_CUR = TPFLYER%ALT
            TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(1,1,IKB) )
            TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(2,1,IKB) )
            TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(1,2,IKB) )
            TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(2,2,IKB) )
          ELSE
            ZZCOEF00 = (TPFLYER%ALT - ZZM(1,1,IK00)) / ( ZZM(1,1,IK00+1) - ZZM(1,1,IK00))
            ZZCOEF01 = (TPFLYER%ALT - ZZM(1,2,IK01)) / ( ZZM(1,2,IK01+1) - ZZM(1,2,IK01))
            ZZCOEF10 = (TPFLYER%ALT - ZZM(2,1,IK10)) / ( ZZM(2,1,IK10+1) - ZZM(2,1,IK10))
            ZZCOEF11 = (TPFLYER%ALT - ZZM(2,2,IK11)) / ( ZZM(2,2,IK11+1) - ZZM(2,2,IK11))
            TPFLYER%RHO = FLYER_INTERP(ZRHO)
            TPFLYER%Z_CUR = FLYER_INTERP(ZZM)
          END IF
        ELSE IF ( TPFLYER%PRES /= XUNDEF ) THEN
          ZFLYER_EXN = (TPFLYER%PRES/XP00)**(XRD/XCPD)
          IK00 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(1,1,:)), 1)
          IK01 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(1,2,:)), 1)
          IK10 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(2,1,:)), 1)
          IK11 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(2,2,:)), 1)
          IF (IK00*IK01*IK10*IK11 .EQ. 0) THEN
            TPFLYER%Z_CUR = ZZM(1,1,IKB)
            TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(2,1,IKB) )
            TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(1,2,IKB) )
            TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(2,2,IKB) )
          ELSE
            ZZCOEF00 = (ZFLYER_EXN - ZEXN(1,1,IK00)) / ( ZEXN(1,1,IK00+1) - ZEXN(1,1,IK00))
            ZZCOEF01 = (ZFLYER_EXN - ZEXN(1,2,IK01)) / ( ZEXN(1,2,IK01+1) - ZEXN(1,2,IK01))
            ZZCOEF10 = (ZFLYER_EXN - ZEXN(2,1,IK10)) / ( ZEXN(2,1,IK10+1) - ZEXN(2,1,IK10))
            ZZCOEF11 = (ZFLYER_EXN - ZEXN(2,2,IK11)) / ( ZEXN(2,2,IK11+1) - ZEXN(2,2,IK11))
            TPFLYER%RHO = FLYER_INTERP(ZRHO)
            TPFLYER%Z_CUR = FLYER_INTERP(ZZM)
          END IF
        ELSE
          TPFLYER%RHO = TPFLYER%MASS / TPFLYER%VOLUME
          IK00 = MAX ( COUNT (TPFLYER%RHO <= ZRHO(1,1,:)), 1)
          IK01 = MAX ( COUNT (TPFLYER%RHO <= ZRHO(1,2,:)), 1)
          IK10 = MAX ( COUNT (TPFLYER%RHO <= ZRHO(2,1,:)), 1)
          IK11 = MAX ( COUNT (TPFLYER%RHO <= ZRHO(2,2,:)), 1)
          IF (IK00*IK01*IK10*IK11 .EQ. 0) THEN
            TPFLYER%Z_CUR = ZZM(1,1,IKB)
            TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(2,1,IKB) )
            TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(1,2,IKB) )
            TPFLYER%Z_CUR = MAX ( TPFLYER%Z_CUR , ZZM(2,2,IKB) )
          ELSE
            ZZCOEF00 = (TPFLYER%RHO - ZRHO(1,1,IK00)) / ( ZRHO(1,1,IK00+1) - ZRHO(1,1,IK00))
            ZZCOEF01 = (TPFLYER%RHO - ZRHO(1,2,IK01)) / ( ZRHO(1,2,IK01+1) - ZRHO(1,2,IK01))
            ZZCOEF10 = (TPFLYER%RHO - ZRHO(2,1,IK10)) / ( ZRHO(2,1,IK10+1) - ZRHO(2,1,IK10))
            ZZCOEF11 = (TPFLYER%RHO - ZRHO(2,2,IK11)) / ( ZRHO(2,2,IK11+1) - ZRHO(2,2,IK11))
            TPFLYER%Z_CUR = FLYER_INTERP(ZZM)
          END IF
        END IF
      END IF
    END IF
!
!
!
!*      5.3  Vertical position
!            -----------------
!
    IF (TPFLYER%TYPE=='ISODEN') THEN
      IK00 = MAX ( COUNT (TPFLYER%RHO <= ZRHO(1,1,:)), 1)
      IK01 = MAX ( COUNT (TPFLYER%RHO <= ZRHO(1,2,:)), 1)
      IK10 = MAX ( COUNT (TPFLYER%RHO <= ZRHO(2,1,:)), 1)
      IK11 = MAX ( COUNT (TPFLYER%RHO <= ZRHO(2,2,:)), 1)
    ELSE IF (TPFLYER%TYPE=='RADIOS' .OR. TPFLYER%TYPE=='CVBALL') THEN
      IK00 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZM(1,1,:)), 1)
      IK01 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZM(1,2,:)), 1)
      IK10 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZM(2,1,:)), 1)
      IK11 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZM(2,2,:)), 1)
    ELSE IF (TPFLYER%TYPE=='AIRCRA') THEN
            IF (TPFLYER%ALTDEF) THEN
              ZFLYER_EXN = (TPFLYER%P_CUR/XP00)**(XRD/XCPD)
              IK00 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(1,1,:)), 1)
              IK01 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(1,2,:)), 1)
              IK10 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(2,1,:)), 1)
              IK11 = MAX ( COUNT (ZFLYER_EXN <= ZEXN(2,2,:)), 1)
            ELSE
              IK00 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZM(1,1,:)), 1)
              IK01 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZM(1,2,:)), 1)
              IK10 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZM(2,1,:)), 1)
              IK11 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZM(2,2,:)), 1)
            END IF
    END IF
    IK00 = MAX ( IK00, IKB )
    IK01 = MAX ( IK01, IKB )
    IK10 = MAX ( IK10, IKB )
    IK11 = MAX ( IK11, IKB )
!
!
!*      5.4  Crash of the balloon
!            --------------------
!
!
    IF (IK00 <  IKB .OR. IK01 <  IKB .OR. IK10 <  IKB .OR. IK11 <  IKB .OR. &
        IK00 >= IKE .OR. IK01 >= IKE .OR. IK10 >= IKE .OR. IK11 >= IKE  ) THEN
      TPFLYER%CRASH=.TRUE.
    END IF
!
  END IF
!
!
  IF (TPFLYER%CRASH) THEN
    TPFLYER%FLY = .FALSE.
    IF (TPFLYER%TYPE=='AIRCRA' .AND. .NOT. GLAUNCH ) THEN
      WRITE(ILUOUT,*) 'Aircraft ',TPFLYER%TITLE,' flew out of the domain the ', &
                    TDTCUR%nday,'/',TDTCUR%nmonth,'/',                          &
                    TDTCUR%nyear,' at ',TDTCUR%xtime,' sec.'
    ELSE IF (TPFLYER%TYPE /= 'AIRCRA') THEN
      WRITE(ILUOUT,*) 'Balloon ',TPFLYER%TITLE,' crashed the ',                 &
                    TDTCUR%nday,'/',TDTCUR%nmonth,'/',                          &
                    TDTCUR%nyear,' at ',TDTCUR%xtime,' sec.'
    END IF
  ELSE
    IF (TPFLYER%TYPE=='AIRCRA' .AND. .NOT. GLAUNCH .AND. ZTDIST > PTSTEP ) THEN
      WRITE(ILUOUT,*) '-------------------------------------------------------------------'
      WRITE(ILUOUT,*) 'Aircraft ',TPFLYER%TITLE,' flies  in leg',TPFLYER%SEGCURN ,' the ',  &
        TDTCUR%nday,'/',TDTCUR%nmonth,'/',      &
        TDTCUR%nyear,' at ',NINT(TDTCUR%xtime),' sec.'
      WRITE(ILUOUT,*) '-------------------------------------------------------------------'
    ENDIF
!
!----------------------------------------------------------------------------
    IF (ZTHIS_PROC>0.) THEN
!----------------------------------------------------------------------------
!
!*      6.   INITIALIZATIONS FOR INTERPOLATIONS
!            ----------------------------------
!
!*      6.1  Interpolation coefficient for X
!            -------------------------------
!
      ZXCOEF = (TPFLYER%X_CUR - ZXHATM(II)) / (ZXHATM(II+1) - ZXHATM(II))
      ZXCOEF = MAX (0.,MIN(ZXCOEF,1.))
!
!
!*      6.2  Interpolation coefficient for y
!            -------------------------------
!
      ZYCOEF = (TPFLYER%Y_CUR - ZYHATM(IJ)) / (ZYHATM(IJ+1) - ZYHATM(IJ))
      ZYCOEF = MAX (0.,MIN(ZYCOEF,1.))
!
!
!*      6.3  Interpolation coefficients for the 4 suroundings verticals
!            ----------------------------------------------------------
!
      IF (TPFLYER%TYPE=='ISODEN') THEN
        ZZCOEF00 = (TPFLYER%RHO - ZRHO(1,1,IK00)) / ( ZRHO(1,1,IK00+1) - ZRHO(1,1,IK00) )
        ZZCOEF01 = (TPFLYER%RHO - ZRHO(1,2,IK01)) / ( ZRHO(1,2,IK01+1) - ZRHO(1,2,IK01) )
        ZZCOEF10 = (TPFLYER%RHO - ZRHO(2,1,IK10)) / ( ZRHO(2,1,IK10+1) - ZRHO(2,1,IK10) )
        ZZCOEF11 = (TPFLYER%RHO - ZRHO(2,2,IK11)) / ( ZRHO(2,2,IK11+1) - ZRHO(2,2,IK11) )
        TPFLYER%Z_CUR = FLYER_INTERP(ZZM)
      ELSE IF (TPFLYER%TYPE=='RADIOS' .OR. TPFLYER%TYPE=='CVBALL') THEN
        ZZCOEF00 = (TPFLYER%Z_CUR - ZZM(1,1,IK00)) / ( ZZM(1,1,IK00+1) - ZZM(1,1,IK00) )
        ZZCOEF01 = (TPFLYER%Z_CUR - ZZM(1,2,IK01)) / ( ZZM(1,2,IK01+1) - ZZM(1,2,IK01) )
        ZZCOEF10 = (TPFLYER%Z_CUR - ZZM(2,1,IK10)) / ( ZZM(2,1,IK10+1) - ZZM(2,1,IK10) )
        ZZCOEF11 = (TPFLYER%Z_CUR - ZZM(2,2,IK11)) / ( ZZM(2,2,IK11+1) - ZZM(2,2,IK11) )
      ELSE IF (TPFLYER%TYPE=='AIRCRA') THEN
              IF (TPFLYER%ALTDEF) THEN
        ZZCOEF00 = (ZFLYER_EXN - ZEXN(1,1,IK00)) / ( ZEXN(1,1,IK00+1) - ZEXN(1,1,IK00) )
        ZZCOEF01 = (ZFLYER_EXN - ZEXN(1,2,IK01)) / ( ZEXN(1,2,IK01+1) - ZEXN(1,2,IK01) )
        ZZCOEF10 = (ZFLYER_EXN - ZEXN(2,1,IK10)) / ( ZEXN(2,1,IK10+1) - ZEXN(2,1,IK10) )
        ZZCOEF11 = (ZFLYER_EXN - ZEXN(2,2,IK11)) / ( ZEXN(2,2,IK11+1) - ZEXN(2,2,IK11) )
        TPFLYER%Z_CUR = FLYER_INTERP(ZZM)
                      ELSE
        ZZCOEF00 = (TPFLYER%Z_CUR - ZZM(1,1,IK00)) / ( ZZM(1,1,IK00+1) - ZZM(1,1,IK00) )
        ZZCOEF01 = (TPFLYER%Z_CUR - ZZM(1,2,IK01)) / ( ZZM(1,2,IK01+1) - ZZM(1,2,IK01) )
        ZZCOEF10 = (TPFLYER%Z_CUR - ZZM(2,1,IK10)) / ( ZZM(2,1,IK10+1) - ZZM(2,1,IK10) )
        ZZCOEF11 = (TPFLYER%Z_CUR - ZZM(2,2,IK11)) / ( ZZM(2,2,IK11+1) - ZZM(2,2,IK11) )
        TPFLYER%P_CUR = FLYER_INTERP(PP)
              END IF
      END IF
!
!----------------------------------------------------------------------------
!
!*      7.   INITIALIZATIONS FOR INTERPOLATIONS OF U AND V
!            ---------------------------------------------
!
!*      7.1  Interpolation coefficient for X (for U)
!            -------------------------------
!
      ZUCOEF = (TPFLYER%X_CUR - PXHAT(IU)) / (PXHAT(IU+1) - PXHAT(IU))
      ZUCOEF = MAX(0.,MIN(ZUCOEF,1.))
!
!
!*      7.2  Interpolation coefficient for y (for V)
!            -------------------------------
!
      ZVCOEF = (TPFLYER%Y_CUR - PYHAT(IV)) / (PYHAT(IV+1) - PYHAT(IV))
      ZVCOEF = MAX(0.,MIN(ZVCOEF,1.))
!
!
!*      7.3  Interpolation coefficients for the 4 suroundings verticals (for U)
!            ----------------------------------------------------------
!
      IU00 = MAX( COUNT (TPFLYER%Z_CUR >= ZZU(1,1,:)), 1)
      IU01 = MAX( COUNT (TPFLYER%Z_CUR >= ZZU(1,2,:)), 1)
      IU10 = MAX( COUNT (TPFLYER%Z_CUR >= ZZU(2,1,:)), 1)
      IU11 = MAX( COUNT (TPFLYER%Z_CUR >= ZZU(2,2,:)), 1)
      ZUCOEF00 = (TPFLYER%Z_CUR - ZZU(1,1,IU00)) / ( ZZU(1,1,IU00+1) - ZZU(1,1,IU00) )
      ZUCOEF01 = (TPFLYER%Z_CUR - ZZU(1,2,IU01)) / ( ZZU(1,2,IU01+1) - ZZU(1,2,IU01) )
      ZUCOEF10 = (TPFLYER%Z_CUR - ZZU(2,1,IU10)) / ( ZZU(2,1,IU10+1) - ZZU(2,1,IU10) )
      ZUCOEF11 = (TPFLYER%Z_CUR - ZZU(2,2,IU11)) / ( ZZU(2,2,IU11+1) - ZZU(2,2,IU11) )
!
!
!*      7.4  Interpolation coefficients for the 4 suroundings verticals (for V)
!            ----------------------------------------------------------
!

      IV00 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZV(1,1,:)), 1)
      IV01 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZV(1,2,:)), 1)
      IV10 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZV(2,1,:)), 1)
      IV11 = MAX ( COUNT (TPFLYER%Z_CUR >= ZZV(2,2,:)), 1)
      ZVCOEF00 = (TPFLYER%Z_CUR - ZZV(1,1,IV00)) / ( ZZV(1,1,IV00+1) - ZZV(1,1,IV00) )
      ZVCOEF01 = (TPFLYER%Z_CUR - ZZV(1,2,IV01)) / ( ZZV(1,2,IV01+1) - ZZV(1,2,IV01) )
      ZVCOEF10 = (TPFLYER%Z_CUR - ZZV(2,1,IV10)) / ( ZZV(2,1,IV10+1) - ZZV(2,1,IV10) )
      ZVCOEF11 = (TPFLYER%Z_CUR - ZZV(2,2,IV11)) / ( ZZV(2,2,IV11+1) - ZZV(2,2,IV11) )
!
!----------------------------------------------------------------------------
!
!*      8.   DATA RECORDING
!            --------------
!
      IF ( GSTORE ) THEN
        TPFLYER%X   (IN) = TPFLYER%X_CUR
        TPFLYER%Y   (IN) = TPFLYER%Y_CUR
        TPFLYER%Z   (IN) = TPFLYER%Z_CUR
        !
        CALL SM_LATLON(PLATOR,PLONOR,          &
                     TPFLYER%X_CUR, TPFLYER%Y_CUR,       &
                     TPFLYER%YLAT(IN), TPFLYER%XLON(IN)  )
        !
        ZU_BAL = FLYER_INTERP_U(PU)
        ZV_BAL = FLYER_INTERP_V(PV)
        ZGAM   = (XRPK * (TPFLYER%XLON(IN) - XLON0) - XBETA)*(XPI/180.)
        TPFLYER%ZON (IN) = ZU_BAL * COS(ZGAM) + ZV_BAL * SIN(ZGAM)
        TPFLYER%MER (IN) = - ZU_BAL * SIN(ZGAM) + ZV_BAL * COS(ZGAM)
        !
        TPFLYER%W   (IN) = FLYER_INTERP(ZWM)
        TPFLYER%TH  (IN) = FLYER_INTERP(PTH)
        !
        ZFLYER_EXN = FLYER_INTERP(ZEXN)
        TPFLYER%P   (IN) = XP00 * ZFLYER_EXN**(XCPD/XRD)
        !
        DO JLOOP=1,SIZE(PR,4)
          TPFLYER%R   (IN,JLOOP) = FLYER_INTERP(PR(:,:,:,JLOOP))
          IF (JLOOP>=2) ZR(:,:,:) = ZR(:,:,:) + PR(:,:,:,JLOOP)
        END DO
        DO JLOOP=1,SIZE(PSV,4)
          TPFLYER%SV  (IN,JLOOP) = FLYER_INTERP(PSV(:,:,:,JLOOP))
        END DO
        TPFLYER%RTZ  (IN,:) = FLYER_INTERPZ(ZR(:,:,:))
        DO JLOOP=1,SIZE(PR,4)
          TPFLYER%RZ  (IN,:,JLOOP) = FLYER_INTERPZ(PR(:,:,:,JLOOP))
        END DO
        ! Fin Modifs ON
        TPFLYER%FFZ  (IN,:) = FLYER_INTERPZ(SQRT(PU**2+PV**2))
        IF (CCLOUD=="LIMA") THEN                                  
          TPFLYER%CIZ  (IN,:) = FLYER_INTERPZ(PSV(:,:,:,NSV_LIMA_NI))  
          TPFLYER%CCZ  (IN,:) = FLYER_INTERPZ(PSV(:,:,:,NSV_LIMA_NC))  
          TPFLYER%CRZ  (IN,:) = FLYER_INTERPZ(PSV(:,:,:,NSV_LIMA_NR))  
        ELSE IF ( CCLOUD=="ICE3" .OR. CCLOUD=="ICE4" ) THEN
          TPFLYER%CIZ  (IN,:) = FLYER_INTERPZ(PCIT(:,:,:))
        ENDIF             
        ! initialization CRARE and CRARE_ATT + LWC and IWC
        TPFLYER%CRARE(IN,:) = 0.
        TPFLYER%CRARE_ATT(IN,:) = 0.
        TPFLYER%LWCZ  (IN,:) = 0.
        TPFLYER%IWCZ  (IN,:) = 0.
      IF (CCLOUD=="LIMA" .OR. CCLOUD=="ICE3" ) THEN ! only for ICE3 and LIMA
       TPFLYER%LWCZ  (IN,:) = FLYER_INTERPZ((PR(:,:,:,2)+PR(:,:,:,3))*PRHODREF(:,:,:))
       TPFLYER%IWCZ  (IN,:) = FLYER_INTERPZ((PR(:,:,:,4)+PR(:,:,:,5)+PR(:,:,:,6))*PRHODREF(:,:,:))
       ZTEMPZ(:)=FLYER_INTERPZ(PTH(II:II+1,IJ:IJ+1,:) * ZEXN(:,:,:))
        ZRHODREFZ(:)=FLYER_INTERPZ(PRHODREF(:,:,:))
        IF (CCLOUD=="LIMA") THEN
          ZCCI(:)=FLYER_INTERPZ(PSV(:,:,:,NSV_LIMA_NI))
          ZCCR(:)=FLYER_INTERPZ(PSV(:,:,:,NSV_LIMA_NR))
          ZCCC(:)=FLYER_INTERPZ(PSV(:,:,:,NSV_LIMA_NC))
        ELSE
          ZCIT(:)=FLYER_INTERPZ(PCIT(:,:,:))
        ENDIF
        DO JLOOP=3,6
          ZRZ(:,JLOOP)=FLYER_INTERPZ(PR(:,:,:,JLOOP))
        END DO
        if ( csurf == 'EXTE' ) then
          DO JK=1,IKU
            ZRZ(JK,2)=FLYER_INTERP_2D(PR(:,:,JK,2)*PSEA(:,:))       ! becomes cloud mixing ratio over sea
            ZRZ(JK,7)=FLYER_INTERP_2D(PR(:,:,JK,2)*(1.-PSEA(:,:)))  ! becomes cloud mixing ratio over land
          END DO
        else
          !if csurf/='EXTE', psea is not allocated
          DO JK=1,IKU
            ZRZ(JK,2)=FLYER_INTERP_2D(PR(:,:,JK,2))
            ZRZ(JK,7) = 0.
          END DO
        end if
        ALLOCATE(ZAELOC(IKU))
        !
        ZAELOC(:)=0.
        ! initialization of quadrature points and weights
        ALLOCATE(ZX(JPTS_GAULAG),ZW(JPTS_GAULAG))
        CALL GAULAG(JPTS_GAULAG,ZX,ZW) ! for integration over diameters
        ! initialize minimum values
        ALLOCATE(ZRTMIN(SIZE(PR,4)+1))
        IF (CCLOUD == 'LIMA') THEN
          ZRTMIN(2)=XRTMIN_L(2) ! cloud water over sea
          ZRTMIN(3)=XRTMIN_L(3)
          ZRTMIN(4)=XRTMIN_L(4)
          ZRTMIN(5)=1E-10
          ZRTMIN(6)=XRTMIN_L(6)
          ZRTMIN(7)=XRTMIN_L(2) ! cloud water over land
        ELSE
          ZRTMIN(2)=XRTMIN_I(2) ! cloud water over sea
          ZRTMIN(3)=XRTMIN_I(3)
          ZRTMIN(4)=XRTMIN_I(4)
          ZRTMIN(5)=1E-10
          ZRTMIN(6)=XRTMIN_I(6)
          ZRTMIN(7)=XRTMIN_I(2) ! cloud water over land
        ENDIF
        ! compute cloud radar reflectivity from vertical profiles of temperature and mixing ratios
        DO JK=1,IKU
          QMW=SQRT(QEPSW(ZTEMPZ(JK),XLIGHTSPEED/XLAM_CRAD))
          QMI=SQRT(QEPSI(ZTEMPZ(JK),XLIGHTSPEED/XLAM_CRAD))
          DO JLOOP=2,7
            IF (CCLOUD == 'LIMA') THEN
              GCALC=(ZRZ(JK,JLOOP)>ZRTMIN(JLOOP).AND.(JLOOP.NE.4.OR.ZCCI(JK)>0.).AND.&
                    (JLOOP.NE.3.OR.ZCCR(JK)>0.).AND.((JLOOP.NE.2.AND. JLOOP.NE.7).OR.ZCCC(JK)>0.))
            ELSE
              GCALC=(ZRZ(JK,JLOOP)>ZRTMIN(JLOOP).AND.(JLOOP.NE.4.OR.ZCIT(JK)>0.))
            ENDIF
            IF(GCALC) THEN
              SELECT CASE(JLOOP)
                CASE(2) ! cloud water over sea
                  IF (CCLOUD == 'LIMA') THEN
                    ZA=XAC_L
                    ZB=XBC_L
                    ZCC=ZCCC(JK)*ZRHODREFZ(JK)
                    ZCX=0.
                    ZALPHA=XALPHAC_L
                    ZNU=XNUC_L
                    ZLBEX=1.0/(ZCX-ZB)
                    ZLB=( ZA*ZCC*MOMG(ZALPHA,ZNU,ZB) )**(-ZLBEX)
                  ELSE
                    ZA=XAC_I
                    ZB=XBC_I
                    ZCC=XCONC_SEA
                    ZCX=0.
                    ZALPHA=XALPHAC2_I
                    ZNU=XNUC2_I
                    ZLBEX=1.0/(ZCX-ZB)
                    ZLB=( ZA*ZCC*MOMG(ZALPHA,ZNU,ZB) )**(-ZLBEX)
                  ENDIF
                CASE(3) ! rain water
                  IF (CCLOUD == 'LIMA') THEN
                    ZA=XAR_L
                    ZB=XBR_L
                    ZCC=ZCCR(JK)*ZRHODREFZ(JK)
                    ZCX=0.
                    ZALPHA=XALPHAR_L
                    ZNU=XNUR_L
                    ZLBEX=1.0/(ZCX-ZB)
                    ZLB=( ZA*ZCC*MOMG(ZALPHA,ZNU,ZB) )**(-ZLBEX)
                  ELSE
                    ZA=XAR_I
                    ZB=XBR_I
                    ZCC=XCCR_I
                    ZCX=-1.
                    ZALPHA=XALPHAR_I
                    ZNU=XNUR_I
                    ZLB=XLBR_I
                    ZLBEX=XLBEXR_I
                  ENDIF
                CASE(4) ! pristine ice
                  IF (CCLOUD == 'LIMA') THEN
                    ZA=XAI_L
                    ZB=XBI_L
                    ZCC=ZCCI(JK)*ZRHODREFZ(JK)
                    ZCX=0.
                    ZALPHA=XALPHAI_L
                    ZNU=XNUI_L
                    ZLBEX=1.0/(ZCX-ZB)
                    ZLB=( ZA*ZCC*MOMG(ZALPHA,ZNU,ZB) )**(-ZLBEX) ! because ZCC not included in XLBI
                    ZFW=0
                  ELSE
                    ZA=XAI_I
                    ZB=XBI_I
                    ZCC=ZCIT(JK)
                    ZCX=0.
                    ZALPHA=XALPHAI_I
                    ZNU=XNUI_I
                    ZLBEX=XLBEXI_I
                    ZLB=XLBI_I*ZCC**(-ZLBEX) ! because ZCC not included in XLBI
                    ZFW=0
                  ENDIF                          
                CASE(5) ! snow
                  IF (CCLOUD == 'LIMA') THEN
                    ZA=XAS_L
                    ZB=XBS_L
                    ZCC=XCCS_L
                    ZCX=XCXS_L
                    ZALPHA=XALPHAS_L
                    ZNU=XNUS_L
                    ZNS=XNS_L
                    ZLB=XLBS_L
                    ZLBEX=XLBEXS_L
                    ZFW=0
                  ELSE
                    ZA=XAS_I
                    ZB=XBS_I
                    ZCC=XCCS_I
                    ZCX=XCXS_I
                    ZALPHA=XALPHAS_I
                    ZNU=XNUS_I
                    ZNS=XNS_I
                    ZLB=XLBS_I
                    ZLBEX=XLBEXS_I
                    ZFW=0
                  ENDIF
                CASE(6) ! graupel
                  !If temperature between -10 and 10°C and Mr and Mg over min threshold: melting graupel
                  ! with liquid water fraction Fw=Mr/(Mr+Mg) else dry graupel (Fw=0)    
                  IF( ZTEMPZ(JK) > XTT-10 .AND. ZTEMPZ(JK) < XTT+10 &
                    .AND. ZRZ(JK,3) > ZRTMIN(3) ) THEN
                    ZFW=ZRZ(JK,3)/(ZRZ(JK,3)+ZRZ(JK,JLOOP))
                  ELSE
                    ZFW=0
                  ENDIF
                  IF (CCLOUD == 'LIMA') THEN
                    ZA=XAG_L
                    ZB=XBG_L
                    ZCC=XCCG_L
                    ZCX=XCXG_L
                    ZALPHA=XALPHAG_L
                    ZNU=XNUG_L
                    ZLB=XLBG_L
                    ZLBEX=XLBEXG_L
                  ELSE
                    ZA=XAG_I
                    ZB=XBG_I
                    ZCC=XCCG_I
                    ZCX=XCXG_I
                    ZALPHA=XALPHAG_I
                    ZNU=XNUG_I
                    ZLB=XLBG_I
                    ZLBEX=XLBEXG_I
                  ENDIF                          
                CASE(7) ! cloud water over land
                  IF (CCLOUD == 'LIMA') THEN
                    ZA=XAC_L
                    ZB=XBC_L
                    ZCC=ZCCC(JK)*ZRHODREFZ(JK)
                    ZCX=0.
                    ZALPHA=XALPHAC_L
                    ZNU=XNUC_L
                    ZLBEX=1.0/(ZCX-ZB)
                    ZLB=( ZA*ZCC*MOMG(ZALPHA,ZNU,ZB) )**(-ZLBEX)
                  ELSE
                    ZA=XAC_I
                    ZB=XBC_I
                    ZCC=XCONC_LAND
                    ZCX=0.
                    ZALPHA=XALPHAC_I
                    ZNU=XNUC_I
                    ZLBEX=1.0/(ZCX-ZB)
                    ZLB=( ZA*ZCC*MOMG(ZALPHA,ZNU,ZB) )**(-ZLBEX)
                  ENDIF
              END SELECT
              IF (JLOOP.EQ.5 .AND. ( (CCLOUD=='LIMA'.AND.LSNOW_T_L).OR. &
                                     (CCLOUD=='ICE3'.AND.LSNOW_T_I) ) ) THEN
                 IF (ZTEMPZ(JK)>-10.) THEN
                    ZLBDA = MAX(MIN(XLBDAS_MAX, 10**(14.554-0.0423*(ZTEMPZ(JK)+273.15))),XLBDAS_MIN)
                 ELSE
                    ZLBDA = MAX(MIN(XLBDAS_MAX, 10**(6.226-0.0106*(ZTEMPZ(JK)+273.15))),XLBDAS_MIN)
                 END IF
                 ZN=ZNS*ZRHODREFZ(JK)*ZRZ(JK,JLOOP)*ZLBDA**ZB
              ELSE
                 ZLBDA=ZLB*(ZRHODREFZ(JK)*ZRZ(JK,JLOOP))**ZLBEX
                 ZN=ZCC*ZLBDA**ZCX
              END IF
              ZREFLOC=0.
              ZAETMP=0.
              DO JJ=1,JPTS_GAULAG ! Gauss-Laguerre quadrature
                ZDELTA_EQUIV=ZX(JJ)**(1./ZALPHA)/ZLBDA
                SELECT CASE(JLOOP)
                  CASE(2,3,7)
                    QM=QMW
                  CASE(4,5,6)
                    ! pristine ice, snow, dry graupel
                    ZRHOHYD=MIN(6.*ZA*ZDELTA_EQUIV**(ZB-3.)/XPI,.92*XRHOLW)
                    QM=sqrt(MG(QMI**2,CMPLX(1,0),ZRHOHYD/.92/XRHOLW))
                    
                    ! water inclusions in ice in air
                    QEPSWI=MG(QMW**2,QM**2,ZFW)
                    ! ice in air inclusions in water
                    QEPSIW=MG(QM**2,QMW**2,1.-ZFW)
                  
                    !MG weighted rule (Matrosov 2008)
                    IF(ZFW .LT. 0.37) THEN
                      ZFPW=0
                    ELSE IF(ZFW .GT. 0.63) THEN
                      ZFPW=1
                    ELSE
                      ZFPW=(ZFW-0.37)/(0.63-0.37)
                    ENDIF  
                    QM=sqrt(QEPSWI*(1.-ZFPW)+QEPSIW*ZFPW)
                END SELECT
                CALL BHMIE(XPI/XLAM_CRAD*ZDELTA_EQUIV,QM,ZQEXT,ZQSCA,ZQBACK)
                ZREFLOC=ZREFLOC+ZQBACK*ZX(JJ)**(ZNU-1)*ZDELTA_EQUIV**2*ZW(JJ)
                ZAETMP =ZAETMP +ZQEXT *ZX(JJ)**(ZNU-1)*ZDELTA_EQUIV**2*ZW(JJ)
              END DO
              ZREFLOC=ZREFLOC*(XLAM_CRAD/XPI)**4*ZN/(4.*GAMMA(ZNU)*.93)
              ZAETMP=ZAETMP  *           XPI    *ZN/(4.*GAMMA(ZNU))
              TPFLYER%CRARE(IN,JK)=TPFLYER%CRARE(IN,JK)+ZREFLOC
              ZAELOC(JK)=ZAELOC(JK)+ZAETMP
            END IF

          END DO

        END DO

        ! apply attenuation
        ALLOCATE(ZZMZ(IKU))
        ZZMZ(:)=FLYER_INTERPZ(ZZM(:,:,:))
        ! nadir
        ZAETOT=1.
        DO JK=COUNT(TPFLYER%Z_CUR >= ZZMZ(:)),1,-1
          IF(JK.EQ.COUNT(TPFLYER%Z_CUR >= ZZMZ(:))) THEN
            IF(TPFLYER%Z_CUR<=ZZMZ(JK)+.5*(ZZMZ(JK+1)-ZZMZ(JK))) THEN
              ! only attenuation from ZAELOC(JK)
              ZAETOT=ZAETOT*EXP(-2.*(ZAELOC(JK)*(TPFLYER%Z_CUR-ZZMZ(JK))))
            ELSE
              ! attenuation from ZAELOC(JK) and ZAELOC(JK+1)
              ZAETOT=ZAETOT*EXP(-2.*(ZAELOC(JK+1)*(TPFLYER%Z_CUR-.5*(ZZMZ(JK+1)+ZZMZ(JK))) &
                +ZAELOC(JK)*.5*(ZZMZ(JK+1)-ZZMZ(JK))))
            END IF
          ELSE
            ! attenuation from ZAELOC(JK) and ZAELOC(JK+1)
            ZAETOT=ZAETOT*EXP(-(ZAELOC(JK+1)+ZAELOC(JK))*(ZZMZ(JK+1)-ZZMZ(JK)))
          END IF
          TPFLYER%CRARE_ATT(IN,JK)=TPFLYER%CRARE(IN,JK)*ZAETOT
        END DO
        ! zenith
        ZAETOT=1.
        DO JK = MAX(COUNT(TPFLYER%Z_CUR >= ZZMZ(:)),1)+1,IKU
          IF ( JK .EQ. (MAX(COUNT(TPFLYER%Z_CUR >= ZZMZ(:)),1)+1) ) THEN        
            IF(TPFLYER%Z_CUR>=ZZMZ(JK)-.5*(ZZMZ(JK)-ZZMZ(JK-1))) THEN
              ! only attenuation from ZAELOC(JK)
              ZAETOT=ZAETOT*EXP(-2.*(ZAELOC(JK)*(ZZMZ(JK)-TPFLYER%Z_CUR)))
            ELSE
              ! attenuation from ZAELOC(JK) and ZAELOC(JK-1)
              ZAETOT=ZAETOT*EXP(-2.*(ZAELOC(JK-1)*(.5*(ZZMZ(JK)+ZZMZ(JK-1))-TPFLYER%Z_CUR) &
                +ZAELOC(JK)*.5*(ZZMZ(JK)-ZZMZ(JK-1))))
            END IF
          ELSE
            ! attenuation from ZAELOC(JK) and ZAELOC(JK-1)
            ZAETOT=ZAETOT*EXP(-(ZAELOC(JK-1)+ZAELOC(JK))*(ZZMZ(JK)-ZZMZ(JK-1)))
          END IF
          TPFLYER%CRARE_ATT(IN,JK)=TPFLYER%CRARE(IN,JK)*ZAETOT
        END DO
        TPFLYER%ZZ  (IN,:) = ZZMZ(:)
        DEALLOCATE(ZZMZ,ZAELOC)
        ! m^3 → mm^6/m^3 → dBZ
        WHERE(TPFLYER%CRARE(IN,:)>0)
          TPFLYER%CRARE(IN,:)=10.*LOG10(1.E18*TPFLYER%CRARE(IN,:))
        ELSEWHERE
          TPFLYER%CRARE(IN,:)=XUNDEF
        END WHERE
        WHERE(TPFLYER%CRARE_ATT(IN,:)>0)
          TPFLYER%CRARE_ATT(IN,:)=10.*LOG10(1.E18*TPFLYER%CRARE_ATT(IN,:))
        ELSEWHERE
          TPFLYER%CRARE_ATT(IN,:)=XUNDEF
        END WHERE
        DEALLOCATE(ZX,ZW,ZRTMIN)
      END IF ! end LOOP ICE3
        ! vertical wind
        TPFLYER%WZ  (IN,:) = FLYER_INTERPZ(ZWM(:,:,:))
        IF (SIZE(PTKE)>0) TPFLYER%TKE  (IN)    = FLYER_INTERP(PTKE)
        IF (SIZE(PTS) >0) TPFLYER%TSRAD(IN)    = FLYER_INTERP_2D(PTS)
        IF (LDIAG_IN_RUN) TPFLYER%TKE_DISS(IN) = FLYER_INTERP(XCURRENT_TKE_DISS)
        TPFLYER%ZS(IN)  = FLYER_INTERP_2D(PZ(:,:,1+JPVEXT))
        TPFLYER%THW_FLUX(IN) = FLYER_INTERP(ZTHW_FLUX)
        TPFLYER%RCW_FLUX(IN) = FLYER_INTERP(ZRCW_FLUX)
        DO JLOOP=1,SIZE(PSV,4)
         TPFLYER%SVW_FLUX(IN,JLOOP) = FLYER_INTERP(ZSVW_FLUX(:,:,:,JLOOP))
        END DO
      END IF
!
!----------------------------------------------------------------------------
!
!*      9.   BALLOON ADVECTION
!            -----------------
!
      IF (TPFLYER%TYPE=='RADIOS' .OR. TPFLYER%TYPE=='ISODEN' .OR. TPFLYER%TYPE=='CVBALL') THEN
        ZU_BAL = FLYER_INTERP_U(PU)
        ZV_BAL = FLYER_INTERP_V(PV)
        if ( .not. lcartesian ) then
          ZMAP = FLYER_INTERP_2D(PMAP)
        else
          ZMAP = 1.
        end if
        !
        TPFLYER%X_CUR = TPFLYER%X_CUR   +   ZU_BAL * PTSTEP * ZMAP
        TPFLYER%Y_CUR = TPFLYER%Y_CUR   +   ZV_BAL * PTSTEP * ZMAP
      END IF
      !
      IF (TPFLYER%TYPE=='RADIOS') THEN
        ZW_BAL = FLYER_INTERP(ZWM)
        TPFLYER%Z_CUR = TPFLYER%Z_CUR + ( ZW_BAL + TPFLYER%WASCENT ) * PTSTEP
      END IF
      !
      IF (TPFLYER%TYPE=='CVBALL') THEN
        ZW_BAL = FLYER_INTERP(ZWM)
        ZRO_BAL = FLYER_INTERP(ZRHO)
        ! calculation with a time step of 1 second or less
        IF (INT(PTSTEP) .GT. 1 ) THEN
          DO JK=1,INT(PTSTEP)
            TPFLYER%WASCENT = TPFLYER%WASCENT &
              -  ( 1. / (1. + TPFLYER%INDDRAG ) ) * 1. * &
                 ( XG * ( ( TPFLYER%MASS / TPFLYER%VOLUME ) - ZRO_BAL ) / ( TPFLYER%MASS / TPFLYER%VOLUME ) &
                    + TPFLYER%WASCENT * ABS ( TPFLYER%WASCENT ) * &
                      TPFLYER%DIAMETER * TPFLYER%AERODRAG / ( 2. * TPFLYER%VOLUME ) &
                  )
            TPFLYER%Z_CUR = TPFLYER%Z_CUR + ( ZW_BAL + TPFLYER%WASCENT ) * 1.
          END DO
        END IF
        IF (PTSTEP .GT. INT(PTSTEP)) THEN
            TPFLYER%WASCENT = TPFLYER%WASCENT &
              -  ( 1. / (1. + TPFLYER%INDDRAG ) ) * (PTSTEP-INT(PTSTEP)) * &
                 ( XG * ( ( TPFLYER%MASS / TPFLYER%VOLUME ) - ZRO_BAL ) / ( TPFLYER%MASS / TPFLYER%VOLUME ) &
                    + TPFLYER%WASCENT * ABS ( TPFLYER%WASCENT ) * &
                      TPFLYER%DIAMETER * TPFLYER%AERODRAG / ( 2. * TPFLYER%VOLUME ) &
                  )
            TPFLYER%Z_CUR = TPFLYER%Z_CUR + ( ZW_BAL + TPFLYER%WASCENT ) * (PTSTEP-INT(PTSTEP))
        END IF
      END IF
!
!----------------------------------------------------------------------------
  END IF
!----------------------------------------------------------------------------
!
!*     10.   AIRCRAFT MOVE (computations done on all processors, to limit exchanges)
!            -------------
!
    IF (TPFLYER%TYPE=='AIRCRA') THEN
!
!
!*     10.1  Determination of flight segment
!            -------------------------------
!
      IL = TPFLYER%SEGCURN
      !
      TPFLYER%SEGCURT = TPFLYER%SEGCURT + PTSTEP
      !
       DO WHILE (TPFLYER%SEGCURT>TPFLYER%SEGTIME(IL))
         TPFLYER%SEGCURN = TPFLYER%SEGCURN + 1
         IL = TPFLYER%SEGCURN
         TPFLYER%SEGCURT = TPFLYER%SEGCURT - TPFLYER%SEGTIME(IL-1)
         IF (IL>TPFLYER%SEG) EXIT
      END DO 
!      DO WHILE (TPFLYER%SEGCURT>TPFLYER%SEGTIME(IL) .AND. IL <= TPFLYER%SEG)
!        TPFLYER%SEGCURN = TPFLYER%SEGCURN + 1
!        IL = TPFLYER%SEGCURN
!        TPFLYER%SEGCURT = TPFLYER%SEGCURT - TPFLYER%SEGTIME(IL-1)
!      END DO
      !
      !* end of flight
      !
      IF (IL > TPFLYER%SEG) TPFLYER%FLY=.FALSE.
!
!
!*     10.2  Determination of new position
!            -----------------------------
!
      IF (TPFLYER%FLY) THEN
        ZSEG_FRAC = TPFLYER%SEGCURT / TPFLYER%SEGTIME(IL)
        !
        TPFLYER%X_CUR = (1.-ZSEG_FRAC) * TPFLYER%SEGX(IL  ) &
                      +     ZSEG_FRAC  * TPFLYER%SEGX(IL+1)
        TPFLYER%Y_CUR = (1.-ZSEG_FRAC) * TPFLYER%SEGY(IL  ) &
                      +     ZSEG_FRAC  * TPFLYER%SEGY(IL+1)
          IF (TPFLYER%ALTDEF) THEN
             TPFLYER%P_CUR = (1.-ZSEG_FRAC) * TPFLYER%SEGP(IL  ) &
                      +     ZSEG_FRAC  * TPFLYER%SEGP(IL+1)
          ELSE
             TPFLYER%Z_CUR = (1.-ZSEG_FRAC) * TPFLYER%SEGZ(IL  ) &
                      +     ZSEG_FRAC  * TPFLYER%SEGZ(IL+1) 
          END IF
      END IF
    !
    END IF
  !
  END IF
! 
END IF
!
!----------------------------------------------------------------------------
!
!*     11.   EXCHANGE OF INFORMATION BETWEEN PROCESSORS
!            ------------------------------------------
!
!*     11.1  current position
!            ----------------
!
CALL DISTRIBUTE_FLYER_L(TPFLYER%FLY)
CALL DISTRIBUTE_FLYER_L(TPFLYER%CRASH)
CALL DISTRIBUTE_FLYER(TPFLYER%X_CUR)
CALL DISTRIBUTE_FLYER(TPFLYER%Y_CUR)
IF (TPFLYER%TYPE=='CVBALL') THEN
  CALL DISTRIBUTE_FLYER(TPFLYER%Z_CUR)
  CALL DISTRIBUTE_FLYER(TPFLYER%WASCENT)
ELSE
  IF (TPFLYER%TYPE=='RADIOS') CALL DISTRIBUTE_FLYER(TPFLYER%Z_CUR)
  IF (TPFLYER%TYPE=='AIRCRA') THEN
     IF (TPFLYER%ALTDEF) THEN
        CALL DISTRIBUTE_FLYER(TPFLYER%P_CUR)
     ELSE
        CALL DISTRIBUTE_FLYER(TPFLYER%Z_CUR)
     ENDIF
  END IF
  IF (TPFLYER%TYPE=='ISODEN' ) CALL DISTRIBUTE_FLYER(TPFLYER%RHO)
END IF
!
!*     11.2  data stored
!            -----------
!
IF ( GSTORE ) THEN
  CALL DISTRIBUTE_FLYER(TPFLYER%X   (IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%Y   (IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%Z   (IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%XLON(IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%YLAT(IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%ZON (IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%MER (IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%W   (IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%P   (IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%TH  (IN))
  DO JLOOP=1,SIZE(PR,4)
    CALL DISTRIBUTE_FLYER(TPFLYER%R   (IN,JLOOP))
  END DO
  DO JLOOP=1,SIZE(PSV,4)
    CALL DISTRIBUTE_FLYER(TPFLYER%SV  (IN,JLOOP))
  END DO
  DO JLOOP=1,IKU              
    CALL DISTRIBUTE_FLYER(TPFLYER%RTZ (IN,JLOOP))
    DO JLOOP2=1,SIZE(PR,4)
      CALL DISTRIBUTE_FLYER(TPFLYER%RZ (IN,JLOOP,JLOOP2))
    ENDDO
    CALL DISTRIBUTE_FLYER(TPFLYER%FFZ (IN,JLOOP))
    CALL DISTRIBUTE_FLYER(TPFLYER%CIZ (IN,JLOOP))
    IF (CCLOUD== 'LIMA' ) THEN
      CALL DISTRIBUTE_FLYER(TPFLYER%CRZ (IN,JLOOP))
      CALL DISTRIBUTE_FLYER(TPFLYER%CCZ (IN,JLOOP))      
    ENDIF
    CALL DISTRIBUTE_FLYER(TPFLYER%IWCZ (IN,JLOOP))
    CALL DISTRIBUTE_FLYER(TPFLYER%LWCZ (IN,JLOOP))
    CALL DISTRIBUTE_FLYER(TPFLYER%CRARE (IN,JLOOP))
    CALL DISTRIBUTE_FLYER(TPFLYER%CRARE_ATT (IN,JLOOP))
    CALL DISTRIBUTE_FLYER(TPFLYER%WZ (IN,JLOOP))
    CALL DISTRIBUTE_FLYER(TPFLYER%ZZ (IN,JLOOP))
  END DO
  IF (SIZE(PTKE)>0) CALL DISTRIBUTE_FLYER(TPFLYER%TKE  (IN))
  IF (SIZE(PTS) >0) CALL DISTRIBUTE_FLYER(TPFLYER%TSRAD(IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%ZS  (IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%THW_FLUX(IN))
  CALL DISTRIBUTE_FLYER(TPFLYER%RCW_FLUX(IN))
  DO JLOOP=1,SIZE(PSV,4)
    CALL DISTRIBUTE_FLYER(TPFLYER%SVW_FLUX(IN,JLOOP))
  END DO
END IF
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
CONTAINS
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
FUNCTION FLYER_INTERP(PA) RESULT(PB)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA
REAL                               :: PB
!
INTEGER :: JI, JJ
!
IF (SIZE(PA,1)==2) THEN
  JI=1
  JJ=1
ELSE
  JI=II
  JJ=IJ
END IF
!
PB = (1.- ZYCOEF) * (1.-ZXCOEF) * ( (1.-ZZCOEF00) * PA(JI  ,JJ  ,IK00) + ZZCOEF00 * PA(JI  ,JJ  ,IK00+1)) &
   + (1.- ZYCOEF) * (   ZXCOEF) * ( (1.-ZZCOEF10) * PA(JI+1,JJ  ,IK10) + ZZCOEF10 * PA(JI+1,JJ  ,IK10+1)) &
   + (    ZYCOEF) * (1.-ZXCOEF) * ( (1.-ZZCOEF01) * PA(JI  ,JJ+1,IK01) + ZZCOEF01 * PA(JI  ,JJ+1,IK01+1)) &
   + (    ZYCOEF) * (   ZXCOEF) * ( (1.-ZZCOEF11) * PA(JI+1,JJ+1,IK11) + ZZCOEF11 * PA(JI+1,JJ+1,IK11+1))
!
END FUNCTION FLYER_INTERP
!----------------------------------------------------------------------------
FUNCTION FLYER_INTERPZ(PA) RESULT(PB)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA
REAL, DIMENSION(SIZE(PA,3))        :: PB
!
INTEGER :: JI, JJ, JK
!
IF (SIZE(PA,1)==2) THEN
  JI=1
  JJ=1
ELSE
  JI=II
  JJ=IJ
END IF
!
!
DO JK=1,SIZE(PA,3)
 IF ( (PA(JI,JJ,JK) /= XUNDEF) .AND. (PA(JI+1,JJ,JK) /= XUNDEF) .AND. &
      (PA(JI,JJ+1,JK) /= XUNDEF) .AND. (PA(JI+1,JJ+1,JK) /= XUNDEF) ) THEN
   PB(JK) = (1.-ZYCOEF) * (1.-ZXCOEF) *  PA(JI,JJ,JK) + &
        (1.-ZYCOEF) * (ZXCOEF) *  PA(JI+1,JJ,JK)  + &
        (ZYCOEF) * (1.-ZXCOEF) *  PA(JI,JJ+1,JK)  + &
        (ZYCOEF) * (ZXCOEF) *  PA(JI+1,JJ+1,JK) 
 ELSE
   PB(JK) = XUNDEF 
 END IF    
END DO
!
END FUNCTION FLYER_INTERPZ
!----------------------------------------------------------------------------
FUNCTION FLYER_INTERP_U(PA) RESULT(PB)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA
REAL                               :: PB
!
INTEGER :: JI, JJ
!
IF (SIZE(PA,1)==2) THEN
  JI=1
  JJ=1
ELSE
  JI=IU
  JJ=IJ
END IF
!
PB = (1.- ZYCOEF) * (1.-ZUCOEF) * ( (1.-ZUCOEF00) * PA(JI  ,JJ  ,IU00) + ZUCOEF00 * PA(JI  ,JJ  ,IU00+1)) &
   + (1.- ZYCOEF) * (   ZUCOEF) * ( (1.-ZUCOEF10) * PA(JI+1,JJ  ,IU10) + ZUCOEF10 * PA(JI+1,JJ  ,IU10+1)) &
   + (    ZYCOEF) * (1.-ZUCOEF) * ( (1.-ZUCOEF01) * PA(JI  ,JJ+1,IU01) + ZUCOEF01 * PA(JI  ,JJ+1,IU01+1)) &
   + (    ZYCOEF) * (   ZUCOEF) * ( (1.-ZUCOEF11) * PA(JI+1,JJ+1,IU11) + ZUCOEF11 * PA(JI+1,JJ+1,IU11+1))
!
END FUNCTION FLYER_INTERP_U
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
FUNCTION FLYER_INTERP_V(PA) RESULT(PB)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA
REAL                               :: PB
!
INTEGER :: JI, JJ
!
IF (SIZE(PA,1)==2) THEN
  JI=1
  JJ=1
ELSE
  JI=II
  JJ=IV
END IF
!
PB = (1.- ZVCOEF) * (1.-ZXCOEF) * ( (1.-ZVCOEF00) * PA(JI  ,JJ  ,IV00) + ZVCOEF00 * PA(JI  ,JJ  ,IV00+1)) &
   + (1.- ZVCOEF) * (   ZXCOEF) * ( (1.-ZVCOEF10) * PA(JI+1,JJ  ,IV10) + ZVCOEF10 * PA(JI+1,JJ  ,IV10+1)) &
   + (    ZVCOEF) * (1.-ZXCOEF) * ( (1.-ZVCOEF01) * PA(JI  ,JJ+1,IV01) + ZVCOEF01 * PA(JI  ,JJ+1,IV01+1)) &
   + (    ZVCOEF) * (   ZXCOEF) * ( (1.-ZVCOEF11) * PA(JI+1,JJ+1,IV11) + ZVCOEF11 * PA(JI+1,JJ+1,IV11+1))
!
END FUNCTION FLYER_INTERP_V
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
FUNCTION FLYER_INTERP_2D(PA) RESULT(PB)
!
REAL, DIMENSION(:,:), INTENT(IN) :: PA
REAL                             :: PB
!
INTEGER :: JI, JJ
!
IF (SIZE(PA,1)==2) THEN
  JI=1
  JJ=1
ELSE
  JI=II
  JJ=IJ
END IF
!
PB = (1.- ZYCOEF) * (1.-ZXCOEF) * PA(JI  ,JJ  ) &
   + (1.- ZYCOEF) * (   ZXCOEF) * PA(JI+1,JJ  ) &
   + (    ZYCOEF) * (1.-ZXCOEF) * PA(JI  ,JJ+1) &
   + (    ZYCOEF) * (   ZXCOEF) * PA(JI+1,JJ+1)
!
END FUNCTION FLYER_INTERP_2D
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE DISTRIBUTE_FLYER(PA)
!
REAL, INTENT(INOUT) :: PA
!
PA = PA * ZTHIS_PROC
CALL REDUCESUM_ll(PA,IINFO_ll)
!
END SUBROUTINE DISTRIBUTE_FLYER
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE DISTRIBUTE_FLYER_N(KA)
!
INTEGER, INTENT(INOUT) :: KA
REAL                   :: ZA
!
ZA=KA
!
ZA = ZA * ZTHIS_PROC
CALL REDUCESUM_ll(ZA,IINFO_ll)
!
IF (NINT(ZA)/=0) KA=NINT(ZA)
!
END SUBROUTINE DISTRIBUTE_FLYER_N
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE DISTRIBUTE_FLYER_L(OA)
!
LOGICAL, INTENT(INOUT) :: OA
REAL                   :: ZA
!
ZA=0.
IF (OA) ZA=1.
!
CALL REDUCESUM_ll(ZA,IINFO_ll)
!
IF (ZA==0.) THEN
  OA=.FALSE.
ELSE
  OA=.TRUE.
END IF
!
END SUBROUTINE DISTRIBUTE_FLYER_L
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE FLYER_CHANGE_MODEL(IMI)
!
INTEGER, INTENT(IN) :: IMI ! model index
!
INTEGER :: IMK      ! kid model index
INTEGER :: IMODEL   ! TPFLYER model index at the beginning of the subroutine
INTEGER :: IU       ! U flux point balloon position (x index)
INTEGER :: IV       ! V flux point balloon position (y index)
INTEGER :: IU_ABS   ! U flux point balloon  position (in the model)
INTEGER :: IV_ABS   ! V flux point balloon position (in the model)
INTEGER :: IXOR     ! Origin's coordinates of the extended 2 way subdomain
INTEGER :: IYOR     ! Origin's coordinates of the extended 2 way subdomain
INTEGER :: IIB      ! current processor domain sizes
INTEGER :: IJB
INTEGER :: IIE
INTEGER :: IJE
!
!
IMODEL=TPFLYER%NMODEL
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IU=COUNT( PXHAT (:)<=TPFLYER%X_CUR )
IV=COUNT( PYHAT (:)<=TPFLYER%Y_CUR )
ZTHIS_PROC=0.
IF (IU>=IIB .AND. IU<=IIE .AND. IV>=IJB .AND. IV<=IJE) ZTHIS_PROC=1.
IF (ZTHIS_PROC .EQ. 1) THEN
  CALL GET_OR_LL('B',IXOR,IYOR)
  IU_ABS=IU + IXOR - 1
  IV_ABS=IV + IYOR - 1 
  !
  IF (TPFLYER%NMODEL == IMI) THEN
    !
    ! go to the kid model if the flyer location is inside
    ! ------------------
    !
    DO IMK=IMI+1,NMODEL
      IF (NDAD(IMK) == IMI .AND. &
         IU_ABS>=NXOR_ALL(IMK)  .AND. IU_ABS<=NXEND_ALL(IMK)  .AND. &
         IV_ABS>=NYOR_ALL(IMK)  .AND. IV_ABS<=NYEND_ALL(IMK) ) THEN
        TPFLYER%NMODEL = IMK
        !
      END IF
    END DO
    !
  ELSE
    !
    ! come from the kid model if the flyer location is outside
    ! ------------------
    !
    IMK = TPFLYER%NMODEL
    IF (IU_ABS<NXOR_ALL(IMK)  .OR. IU_ABS>NXEND_ALL(IMK)  .OR. &
         IV_ABS<NYOR_ALL(IMK)  .OR. IV_ABS>NYEND_ALL(IMK) ) THEN
        TPFLYER%NMODEL = IMI
        !
    END IF
  END IF
END IF
!
! send the information to all the processors
! ----------------------------------------
!
CALL DISTRIBUTE_FLYER_N(TPFLYER%NMODEL)
ZTHIS_PROC=0.
!
! print if the model changes
!---------------------------------
IF (TPFLYER%NMODEL /= IMODEL) THEN
   IF (NDAD(IMODEL) == TPFLYER%NMODEL) THEN
      WRITE(ILUOUT,*) '-------------------------------------------------------------------'
      WRITE(ILUOUT,*) TPFLYER%TITLE,' comes from model ',IMODEL,' in  model ',&
             TPFLYER%NMODEL,' at ',NINT(TDTCUR%xtime),' sec.'
      WRITE(ILUOUT,*) '-------------------------------------------------------------------'
   ELSE
      WRITE(ILUOUT,*) '-------------------------------------------------------------------'
      WRITE(ILUOUT,*) TPFLYER%TITLE,' goes from model ',IMODEL,' to  model ',&
             TPFLYER%NMODEL,' at ',NINT(TDTCUR%xtime),' sec.'
      WRITE(ILUOUT,*) '-------------------------------------------------------------------'
   ENDIF
ENDIF
!
!
END SUBROUTINE FLYER_CHANGE_MODEL
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
END SUBROUTINE AIRCRAFT_BALLOON_EVOL
