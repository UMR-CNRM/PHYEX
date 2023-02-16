!MNH_LIC Copyright 2002-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ##########################
MODULE MODI_PROFILER_n
!      ##########################
!
INTERFACE
!
      SUBROUTINE PROFILER_n(PTSTEP,                               &
                            PXHAT, PYHAT, PZ,PRHODREF,            &
                            PU, PV, PW, PTH, PR, PSV, PTKE,       &
                            PTS, PP, PAER, PCLDFR, PCIT, PSEA)
!
REAL,                     INTENT(IN)     :: PTSTEP ! time step
REAL, DIMENSION(:),       INTENT(IN)     :: PXHAT  ! x coordinate
REAL, DIMENSION(:),       INTENT(IN)     :: PYHAT  ! y coordinate
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PZ     ! z array
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PU     ! horizontal wind X component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PV     ! horizontal wind Y component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PW     ! vertical wind
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTH    ! potential temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PRHODREF                            
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PR     ! water mixing ratios
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PSV    ! Scalar variables
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTKE   ! turbulent kinetic energy
REAL, DIMENSION(:,:),     INTENT(IN)     :: PTS    ! surface temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PP     ! pressure
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PAER   ! aerosol extinction
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PCLDFR ! cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PCIT   ! ice concentration
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN)     :: PSEA   ! for radar 
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PROFILER_n
!
END INTERFACE
!
END MODULE MODI_PROFILER_n
!
!  ########################################################
      SUBROUTINE PROFILER_n(PTSTEP,                               &
                            PXHAT, PYHAT, PZ,PRHODREF,            &
                            PU, PV, PW, PTH, PR, PSV, PTKE,       &
                            PTS, PP, PAER, PCLDFR, PCIT, PSEA)
!     ########################################################
!
!
!
!!****  *PROFILER_n* - (advects and) stores 
!!                                stations/s in the model
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!    
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
!!      Pierre TULET / Valery Masson             * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!     Original 15/02/2002
!!     March 2013 : C.Lac : Corrections for 1D + new fields (RARE,THV,DD,FF)
!!     April 2014 : C.Lac : Call RADAR only if ICE3   
!!     C.Lac 10/2016  Add visibility diagnostic
!!     March,28, 2018 (P. Wautelet) replace TEMPORAL_DIST by DATETIME_DISTANCE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 13/09/2019: budget: simplify and modernize date/time management
!  M. Taufour  05/07/2021: modify RARE for hydrometeors containing ice and add bright band calculation for RARE
!  P. Wautelet 09/02/2022: add message when some variables not computed
!                          + bugfix: put values in variables in this case
!                          + move some operations outside a do loop
! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_DIAG_IN_RUN
USE MODD_GRID
USE MODD_SUB_PROFILER_n
USE MODD_NSV
USE MODD_PARAMETERS
USE MODD_PARAM_n,        ONLY: CCLOUD, CRAD, CSURF
USE MODD_PROFILER_n
USE MODD_TIME,           only: tdtexp
USE MODD_TIME_n,         only: tdtcur
!
USE MODE_ll
USE MODE_MSG
!
USE MODI_GPS_ZENITH_GRID
USE MODI_LIDAR
USE MODI_RADAR_RAIN_ICE
USE MODI_WATER_SUM
USE MODE_FGAU,             ONLY : GAULAG
USE MODE_FSCATTER,         ONLY: QEPSW,QEPSI,BHMIE,MOMG,MG
USE MODD_PARAM_LIMA,       ONLY: XALPHAR_L=>XALPHAR,XNUR_L=>XNUR,XALPHAS_L=>XALPHAS,XNUS_L=>XNUS,&
                                 XALPHAG_L=>XALPHAG,XNUG_L=>XNUG, XALPHAI_L=>XALPHAI,XNUI_L=>XNUI,&
                                 XRTMIN_L=>XRTMIN,XALPHAC_L=>XALPHAC,XNUC_L=>XNUC, LSNOW_T_L=>LSNOW_T
USE MODD_PARAM_LIMA_COLD,  ONLY: XDI_L=>XDI,XLBEXI_L=>XLBEXI,XLBI_L=>XLBI,XAI_L=>XAI,XBI_L=>XBI,XC_I_L=>XC_I,&
                                 XLBEXS_L=>XLBEXS,XLBS_L=>XLBS,XCCS_L=>XCCS,XNS_L=>XNS,&
                                 XAS_L=>XAS,XBS_L=>XBS,XCXS_L=>XCXS,        &
                                 XLBDAS_MIN,XLBDAS_MAX
USE MODD_PARAM_LIMA_MIXED, ONLY: XDG_L=>XDG,XLBEXG_L=>XLBEXG,XLBG_L=>XLBG,XCCG_L=>XCCG,&
                                 XAG_L=>XAG,XBG_L=>XBG,XCXG_L=>XCXG,XCG_L=>XCG
USE MODD_PARAM_LIMA_WARM,  ONLY: XLBEXR_L=>XLBEXR,XLBR_L=>XLBR,XBR_L=>XBR,XAR_L=>XAR,&
                                 XBC_L=>XBC,XAC_L=>XAC
USE MODD_PARAM_ICE,        ONLY: LSNOW_T_I=>LSNOW_T
USE MODD_RAIN_ICE_DESCR,   ONLY: XALPHAR_I=>XALPHAR,XNUR_I=>XNUR,XLBEXR_I=>XLBEXR,&
                                 XLBR_I=>XLBR,XCCR_I=>XCCR,XBR_I=>XBR,XAR_I=>XAR,&
                                 XALPHAC_I=>XALPHAC,XNUC_I=>XNUC,&
                                 XLBC_I=>XLBC,XBC_I=>XBC,XAC_I=>XAC,&
                                 XALPHAC2_I=>XALPHAC2,XNUC2_I=>XNUC2,&
                                 XALPHAS_I=>XALPHAS,XNUS_I=>XNUS,XLBEXS_I=>XLBEXS,&
                                 XLBS_I=>XLBS,XCCS_I=>XCCS,XNS_I=>XNS,XAS_I=>XAS,XBS_I=>XBS,XCXS_I=>XCXS,&
                                 XALPHAG_I=>XALPHAG,XNUG_I=>XNUG,XDG_I=>XDG,XLBEXG_I=>XLBEXG,&
                                 XLBG_I=>XLBG,XCCG_I=>XCCG,XAG_I=>XAG,XBG_I=>XBG,XCXG_I=>XCXG,XCG_I=>XCG,&
                                 XALPHAI_I=>XALPHAI,XNUI_I=>XNUI,XDI_I=>XDI,XLBEXI_I=>XLBEXI,&
                                 XLBI_I=>XLBI,XAI_I=>XAI,XBI_I=>XBI,XC_I_I=>XC_I,&
                                 XRTMIN_I=>XRTMIN,XCONC_LAND,XCONC_SEA
!
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
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PU     ! horizontal wind X component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PV     ! horizontal wind Y component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PW     ! vertical wind
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTH    ! potential temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PRHODREF                            
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PR     ! water mixing ratios
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PSV    ! Scalar variables
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTKE   ! turbulent kinetic energy
REAL, DIMENSION(:,:),     INTENT(IN)     :: PTS    ! surface temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PP     ! pressure
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PAER   ! aerosol extinction
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PCLDFR ! cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PCIT   ! ice concentration
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN)     :: PSEA   ! for radar 
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!
INTEGER :: IIB        ! current processor domain sizes
INTEGER :: IJB
INTEGER :: IKB
INTEGER :: IIE
INTEGER :: IJE
INTEGER :: IKE
INTEGER :: IIU
INTEGER :: IJU
INTEGER :: IKU
!
!
REAL, DIMENSION(SIZE(PXHAT))        :: ZXHATM ! mass point coordinates
REAL, DIMENSION(SIZE(PYHAT))        :: ZYHATM ! mass point coordinates
!
REAL, DIMENSION(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3),SIZE(PSV,4))  :: ZWORK  
REAL, DIMENSION(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3),SIZE(PAER,4))  :: ZWORK2  
!
LOGICAL :: GSTORE   ! storage occurs at this time step
!
INTEGER :: IN     ! time index
INTEGER :: JSV    ! loop counter
INTEGER :: JK     ! loop
INTEGER :: IKRAD
!
REAL,DIMENSION(SIZE(PZ,3)) :: ZU_PROFILER ! horizontal wind speed profile at station location (along x)
REAL,DIMENSION(SIZE(PZ,3)) :: ZV_PROFILER ! horizontal wind speed profile at station location (along y)
REAL,DIMENSION(SIZE(PZ,3)) :: ZFF         ! horizontal wind speed profile at station location 
REAL,DIMENSION(SIZE(PZ,3)) :: ZDD         ! horizontal wind speed profile at station location 
REAL,DIMENSION(SIZE(PZ,3)) :: ZRHOD       ! dry air density in moist mixing profile at station location
REAL,DIMENSION(SIZE(PZ,3)) :: ZRV         ! water vapour mixing ratio profile at station location
REAL,DIMENSION(SIZE(PZ,3)) :: ZT          ! temperature profile at station location
REAL,DIMENSION(SIZE(PZ,3)) :: ZTV         ! virtual temperature profile at station location
REAL,DIMENSION(SIZE(PZ,3)) :: ZPRES       ! pressure profile at station location
REAL,DIMENSION(SIZE(PZ,3)) :: ZE          ! water vapour partial pressure profile at station location
REAL,DIMENSION(SIZE(PZ,3)) :: ZZ          ! altitude of model levels at station location
REAL,DIMENSION(SIZE(PZ,3)-1) :: ZZHATM      ! altitude of mass point levels at station location
REAL                       :: ZGAM        ! rotation between meso-nh base and spherical lat-lon base.
!
REAL                       :: XZS_GPS       ! GPS station altitude
REAL                       :: ZIWV        ! integrated water vapour at station location
REAL                       :: ZZM_STAT      ! altitude at station location
REAL                       :: ZTM_STAT      ! temperature at station location
REAL                       :: ZTV_STAT      ! virtual temperature at station location
REAL                       :: ZPM_STAT      ! pressure at station location
REAL                       :: ZEM_STAT      ! water vapour partial pressure at station location
REAL                       :: ZZTD_PROFILER ! ZTD at station location
REAL                       :: ZZHD_PROFILER ! ZHD at station location
REAL                       :: ZZWD_PROFILER ! ZWD at station location
REAL                       :: ZZHDR         ! ZHD correction at station location
REAL                       :: ZZWDR         ! ZWD correction at station location
!
INTEGER                    :: IINFO_ll    ! return code
INTEGER                    :: ILUOUT      ! logical unit
INTEGER                    :: IRESP       ! return code
INTEGER                    :: I           ! loop for stations
!
REAL,DIMENSION(SIZE(PTH,1),SIZE(PTH,2))              :: ZZTD,ZZHD,ZZWD
REAL,DIMENSION(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3))  :: ZTEMP,ZRARE,ZTHV,ZTEMPV
REAL,DIMENSION(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3))  :: ZWORK32,ZWORK33,ZWORK34
REAL,DIMENSION(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3))  :: ZVISI,ZVISIKUN
REAL ::  ZK1,ZK2,ZK3            ! k1, k2 and K3 atmospheric refractivity constants
REAL  :: ZRDSRV                 ! XRD/XRV
!
! specific to cloud radar
INTEGER                        :: JLOOP,JLOOP2    ! loop counter
REAL, DIMENSION(SIZE(PR,3))    :: ZTEMPZ! vertical profile of temperature
REAL, DIMENSION(SIZE(PR,3))    :: ZRHODREFZ ! vertical profile of dry air density of the reference state
REAL, DIMENSION(SIZE(PR,3))    :: ZCIT     ! pristine ice concentration
REAL, DIMENSION(SIZE(PR,3))    :: ZCCI,ZCCR,ZCCC     ! ICE,RAIN CLOUD concentration (LIMA)
REAL, DIMENSION(SIZE(PR,1),SIZE(PR,2),SIZE(PR,3))    :: ZR   
REAL, DIMENSION(SIZE(PR,3),SIZE(PR,4)+1) :: ZRZ  ! vertical profile of hydrometeor mixing ratios
REAL                           :: ZA,ZB,ZCC,ZCX,ZALPHA,ZNU,ZLB,ZLBEX,ZRHOHYD,XLAM_CRAD,ZNS   ! generic microphysical parameters
INTEGER                        :: JJ    ! loop counter for quadrature
COMPLEX                        :: QMW,QMI,QM,QB,QEPSIW,QEPSWI   ! dielectric parameter
REAL                           :: ZAETOT,ZAETMP,ZREFLOC,ZQSCA,ZQBACK,ZQEXT ! temporary scattering parameters
REAL,DIMENSION(:),ALLOCATABLE  :: ZAELOC,ZZMZ ! temporary arrays
INTEGER                        :: JPTS_GAULAG=9 ! number of points for Gauss-Laguerre quadrature
REAL                           :: ZLBDA   ! slope distribution parameter
REAL                           :: ZN   ! number cocentration
REAL                           :: ZFRAC_ICE  ! ice water fraction
REAL                           :: ZDELTA_EQUIV ! mass-equivalent Gauss-Laguerre point
REAL                           :: ZFW ! liquid fraction
REAL                           :: ZFPW ! weight for mixed-phase reflectivity
REAL,DIMENSION(:),ALLOCATABLE  :: ZX,ZW ! Gauss-Laguerre points and weights
REAL,DIMENSION(:),ALLOCATABLE  :: ZRTMIN ! local values for XRTMIN
LOGICAL                        :: GCALC
!----------------------------------------------------------------------------
!
!*      2.   PRELIMINARIES
!            -------------
!
!*      2.0   Refractivity coeficients
!             ------------------------
! Bevis et al. (1994)
ZK1 = 0.776       ! K/Pa
ZK2 = 0.704       ! K/Pa
ZK3 = 3739.       ! K2/Pa
ZRDSRV=XRD/XRV
!
XLAM_CRAD        = 3.154E-3 ! (in m) <=> 95.04 GHz = Rasta cloud radar frequency
!*      2.1  Indices
!            -------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKU =   SIZE(PZ,3)     ! nombre de niveaux sur la verticale
IKB = JPVEXT+1
IKE = IKU-JPVEXT
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
!
!----------------------------------------------------------------------------
!
!
!*      3.4  instant of storage
!            ------------------
!
IF ( TPROFILER%T_CUR == XUNDEF ) TPROFILER%T_CUR = TPROFILER%STEP - PTSTEP
!
TPROFILER%T_CUR = TPROFILER%T_CUR + PTSTEP
!
IF ( TPROFILER%T_CUR >= TPROFILER%STEP - 1.E-10 ) THEN
  GSTORE = .TRUE.
  TPROFILER%T_CUR = TPROFILER%T_CUR - TPROFILER%STEP
  TPROFILER%N_CUR = TPROFILER%N_CUR + 1
  IN = TPROFILER%N_CUR
ELSE
  GSTORE = .FALSE.
END IF
!
IF (GSTORE) THEN
#if 0
  tprofiler%tpdates(in)%date%year  = tdtexp%date%year
  tprofiler%tpdates(in)%date%month = tdtexp%date%month
  tprofiler%tpdates(in)%date%day   = tdtexp%date%day
  tprofiler%tpdates(in)%xtime      = tdtexp%xtime + ( in - 1 ) * tprofiler%step
#else
  tprofiler%tpdates(in) = tdtcur
#endif
END IF
!
!
!----------------------------------------------------------------------------
!
!*      4.   PROFILER POSITION
!            --------------
!
!*      4.0  initialization of processor test
!            --------------------------------
IF (GPROFILERFIRSTCALL) THEN
GPROFILERFIRSTCALL=.FALSE.
!
 IF (.NOT.(ASSOCIATED(ZTHIS_PROCS))) ALLOCATE(ZTHIS_PROCS(NUMBPROFILER))
!
IF (.NOT.(ASSOCIATED(II))) ALLOCATE(II(NUMBPROFILER))
IF (.NOT.(ASSOCIATED(IJ))) ALLOCATE(IJ(NUMBPROFILER))
IF (.NOT.(ASSOCIATED(IV))) ALLOCATE(IV(NUMBPROFILER))
IF (.NOT.(ASSOCIATED(IU))) ALLOCATE(IU(NUMBPROFILER))
IF (.NOT.(ASSOCIATED(ZXCOEF))) ALLOCATE(ZXCOEF(NUMBPROFILER))
IF (.NOT.(ASSOCIATED(ZUCOEF))) ALLOCATE(ZUCOEF(NUMBPROFILER))
IF (.NOT.(ASSOCIATED(ZYCOEF))) ALLOCATE(ZYCOEF(NUMBPROFILER))
IF (.NOT.(ASSOCIATED(ZVCOEF))) ALLOCATE(ZVCOEF(NUMBPROFILER))
!
ZXCOEF(:)=XUNDEF
ZUCOEF(:)=XUNDEF
ZYCOEF(:)=XUNDEF
ZVCOEF(:)=XUNDEF
!
DO I=1,NUMBPROFILER

ZTHIS_PROCS(I)=0.
!
!*      4.1  X position
!            ----------
!
IU(I)=COUNT( PXHAT (:)<=TPROFILER%X(I) )
II(I)=COUNT( ZXHATM(:)<=TPROFILER%X(I) )
!
IF (II(I)<=IIB-1   .AND. LWEST_ll() .AND. .NOT. L1D) TPROFILER%ERROR(I)=.TRUE.
IF (II(I)>=IIE     .AND. LEAST_ll() .AND. .NOT. L1D) TPROFILER%ERROR(I)=.TRUE.
!
!
!*      4.2  Y position
!            ----------
!
IV(I)=COUNT( PYHAT (:)<=TPROFILER%Y(I) )
IJ(I)=COUNT( ZYHATM(:)<=TPROFILER%Y(I) )
!
IF (IJ(I)<=IJB-1   .AND. LSOUTH_ll() .AND. .NOT. L1D) TPROFILER%ERROR(I)=.TRUE.
IF (IJ(I)>=IJE     .AND. LNORTH_ll() .AND. .NOT. L1D) TPROFILER%ERROR(I)=.TRUE.
!
!
!*      4.3  Position of station according to processors
!            -------------------------------------------
!
IF (IU(I)>=IIB .AND. IU(I)<=IIE .AND. IV(I)>=IJB .AND. IV(I)<=IJE) ZTHIS_PROCS(I)=1.
IF (L1D) ZTHIS_PROCS(I)=1.
!
!
!*      4.4  Computations only on correct processor
!            --------------------------------------
ZXCOEF(I) = 0.
ZYCOEF(I) = 0.
ZUCOEF(I) = 0.         
ZVCOEF(I) = 0.
IF (ZTHIS_PROCS(I) >0. .AND. .NOT. L1D) THEN
!
!*      6.1  Interpolation coefficient for X
!            -------------------------------
!
  ZXCOEF(I) = (TPROFILER%X(I) - ZXHATM(II(I))) / (ZXHATM(II(I)+1) - ZXHATM(II(I)))
!
!
!
!*      6.2  Interpolation coefficient for y
!            -------------------------------
!
  ZYCOEF(I) = (TPROFILER%Y(I) - ZYHATM(IJ(I))) / (ZYHATM(IJ(I)+1) - ZYHATM(IJ(I)))
!
!----------------------------------------------------------------------------
!
!*      7.   INITIALIZATIONS FOR INTERPOLATIONS OF U AND V
!            ---------------------------------------------
!
!*      7.1  Interpolation coefficient for X (for U)
!            -------------------------------
!
  ZUCOEF(I) = (TPROFILER%X(I) - PXHAT(IU(I))) / (PXHAT(IU(I)+1) - PXHAT(IU(I)))
!
!*      7.2  Interpolation coefficient for y (for V)
!            -------------------------------
!
  ZVCOEF(I) = (TPROFILER%Y(I) - PYHAT(IV(I))) / (PYHAT(IV(I)+1) - PYHAT(IV(I)))
!
END IF
ENDDO
END IF
!----------------------------------------------------------------------------
!
!*      8.   DATA RECORDING
!            --------------
!
ZTEMP(:,:,:)=PTH(:,:,:)*(PP(:,:,:)/ XP00) **(XRD/XCPD)
! Theta_v
ZTHV(:,:,:) = PTH(:,:,:) / (1.+WATER_SUM(PR(:,:,:,:)))*(1.+PR(:,:,:,1)/ZRDSRV)
! virtual temperature
ZTEMPV(:,:,:)=ZTHV(:,:,:)*(PP(:,:,:)/ XP00) **(XRD/XCPD)
CALL GPS_ZENITH_GRID(PR(:,:,:,1),ZTEMP,PP,ZZTD,ZZHD,ZZWD)
! Kunkel formulation
IF (SIZE(PR,4) >= 2) THEN
  ZVISIKUN(:,:,:) = 10E5  !default value
  WHERE ( PR(:,:,:,2) /=0 )
    ZVISIKUN(:,:,:) =0.027/(10**(-8)+(PR(:,:,:,2)/(1+PR(:,:,:,2))*PRHODREF(:,:,:)*1000))**0.88
  END WHERE
END IF
! Gultepe formulation
IF ((SIZE(PR,4) >= 2) .AND. NSV_C2R2END /= 0 ) THEN 
  WHERE ( (PR(:,:,:,2) /=0. ) .AND. (PSV(:,:,:,NSV_C2R2BEG+1) /=0. ) )
    ZVISI(:,:,:) =1.002/(PR(:,:,:,2)*PRHODREF(:,:,:)*PSV(:,:,:,NSV_C2R2BEG+1))**0.6473
  END WHERE
END IF
!
IF (GSTORE) THEN
  DO I=1,NUMBPROFILER
    IF ((ZTHIS_PROCS(I)==1.).AND.(.NOT. TPROFILER%ERROR(I))) THEN
      !
      ZZ(:)                  = PROFILER_INTERP(PZ)
      ZRHOD(:)               = PROFILER_INTERP(PRHODREF)
      ZPRES(:)               = PROFILER_INTERP(PP)
      ZU_PROFILER(:)         = PROFILER_INTERP_U(PU)
      ZV_PROFILER(:)         = PROFILER_INTERP_V(PV)
      ZGAM                   = (XRPK * (TPROFILER%LON(I) - XLON0) - XBETA)*(XPI/180.)
      ZFF(:)                 = SQRT(ZU_PROFILER(:)**2 + ZV_PROFILER(:)**2)
      DO JK=1,IKU
       IF (ZU_PROFILER(JK) >=0. .AND. ZV_PROFILER(JK) > 0.) &
         ZDD(JK) = ATAN(ABS(ZU_PROFILER(JK)/ZV_PROFILER(JK))) * 180./XPI + 180.
       IF (ZU_PROFILER(JK) >0. .AND. ZV_PROFILER(JK) <= 0.) &
         ZDD(JK) = ATAN(ABS(ZV_PROFILER(JK)/ZU_PROFILER(JK))) * 180./XPI + 270.
       IF (ZU_PROFILER(JK) <=0. .AND. ZV_PROFILER(JK) < 0.) &
         ZDD(JK) = ATAN(ABS(ZU_PROFILER(JK)/ZV_PROFILER(JK))) * 180./XPI 
       IF (ZU_PROFILER(JK) <0. .AND. ZV_PROFILER(JK) >= 0.) &
         ZDD(JK) = ATAN(ABS(ZV_PROFILER(JK)/ZU_PROFILER(JK))) * 180./XPI + 90.
       IF (ZU_PROFILER(JK) == 0. .AND. ZV_PROFILER(JK) == 0.) &
         ZDD(JK) = XUNDEF
      END DO
      ! GPS IWV and ZTD
      XZS_GPS=TPROFILER%ALT(I)
      IF ( ABS( ZZ(IKB)-XZS_GPS ) < 150 ) THEN ! distance between real and model orography ok
        ZRV(:)                 = PROFILER_INTERP(PR(:,:,:,1))
        ZT(:)                  = PROFILER_INTERP(ZTEMP)
        ZE(:)                  = ZPRES(:)*ZRV(:)/(ZRDSRV+ZRV(:))
        ZTV(:)                 = PROFILER_INTERP(ZTEMPV)
        ZZTD_PROFILER          = PROFILER_INTERP_2D(ZZTD)
        ZZHD_PROFILER          = PROFILER_INTERP_2D(ZZHD)
        ZZWD_PROFILER          = PROFILER_INTERP_2D(ZZWD)
        ZIWV = 0.
        DO JK=IKB,IKE
         ZIWV=ZIWV+ZRHOD(JK)*ZRV(JK)*(ZZ(JK+1)-ZZ(JK))
        END DO
        IF (ZZ(IKB) < XZS_GPS) THEN ! station above the model orography     
          DO JK=IKB+1,IKE
            IF ( ZZ(JK) < XZS_GPS) THEN ! whole layer to remove
              ZZHDR=( 1.E-6 * ZK1 * ZPRES(JK-1) * ( ZZ(JK) - ZZ(JK-1) ) / ZTV(JK-1))
              ZZWDR=( 1.E-6 *  ( (ZK2-ZRDSRV*ZK1) + ( ZK3/ZT(JK-1) ) ) * &
                 ZE(JK-1)* ( ZZ(JK) - ZZ(JK-1) ) / ZT(JK-1) )
              ZZHD_PROFILER=ZZHD_PROFILER-ZZHDR
              ZZWD_PROFILER=ZZWD_PROFILER-ZZWDR
              ZZTD_PROFILER=ZZTD_PROFILER-ZZHDR-ZZWDR
            ELSE                       ! partial layer to remove
              ZZHDR=( 1.E-6 * ZK1 * ZPRES(JK-1) * ( XZS_GPS - ZZ(JK-1) ) / ZTV(JK-1)) 
              ZZWDR=( 1.E-6 *  ( (ZK2-ZRDSRV*ZK1) + ( ZK3/ZT(JK-1) ) ) * &
                 ZE(JK-1)* ( XZS_GPS - ZZ(JK-1) ) / ZT(JK-1) ) 
              ZZHD_PROFILER=ZZHD_PROFILER-ZZHDR
              ZZWD_PROFILER=ZZWD_PROFILER-ZZWDR
              ZZTD_PROFILER=ZZTD_PROFILER-ZZHDR-ZZWDR
              EXIT
            END IF
          END DO 
        ELSE ! station below the model orography
! Extrapolate variables below the model orography assuming constant T&Tv gradients,
! constant rv and hydrostatic law
          ZZHATM(:)=0.5*(ZZ(1:IKU-1)+ZZ(2:IKU))
          ZZM_STAT=0.5*(XZS_GPS+ZZ(IKB))
          ZTM_STAT=ZT(IKB) + ( (ZZM_STAT-ZZHATM(IKB))*&
             ( ZT(IKB)- ZT(IKB+1) )/(ZZHATM(IKB)-ZZHATM(IKB+1)) )
          ZTV_STAT=ZTV(IKB) + ( (ZZM_STAT-ZZHATM(IKB))*&
             ( ZTV(IKB)- ZTV(IKB+1) )/(ZZHATM(IKB)-ZZHATM(IKB+1)) )
          ZPM_STAT = ZPRES(IKB) * EXP(XG *(ZZM_STAT-ZZHATM(IKB))&
             /(XRD* 0.5 *(ZTV_STAT+ZTV(IKB))))
          ZEM_STAT = ZPM_STAT * ZRV(IKB) / ( ZRDSRV + ZRV(IKB) )
! add contribution below the model orography        
          ZZHDR=( 1.E-6 * ZK1 * ZPM_STAT * ( ZZ(IKB) - XZS_GPS ) / ZTV_STAT )
          ZZWDR=( 1.E-6 * ( (ZK2-ZRDSRV*ZK1) + (ZK3/ZTM_STAT) )&
             * ZEM_STAT* ( ZZ(IKB) - XZS_GPS ) / ZTM_STAT )
          ZZHD_PROFILER=ZZHD_PROFILER+ZZHDR
          ZZWD_PROFILER=ZZWD_PROFILER+ZZWDR
          ZZTD_PROFILER=ZZTD_PROFILER+ZZHDR+ZZWDR
        END IF
        TPROFILER%IWV(IN,I)= ZIWV
        TPROFILER%ZTD(IN,I)= ZZTD_PROFILER
        TPROFILER%ZWD(IN,I)= ZZWD_PROFILER
        TPROFILER%ZHD(IN,I)= ZZHD_PROFILER
      ELSE
        CMNHMSG(1) = 'altitude of profiler ' // TRIM( TPROFILER%NAME(I) ) // ' is too far from orography'
        CMNHMSG(2) = 'some variables are therefore not computed (IWV, ZTD, ZWD, ZHD)'
        CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'PROFILER_n' )
        TPROFILER%IWV(IN,I)= XUNDEF
        TPROFILER%ZTD(IN,I)= XUNDEF
        TPROFILER%ZWD(IN,I)= XUNDEF
        TPROFILER%ZHD(IN,I)= XUNDEF
      END IF
      TPROFILER%ZON (IN,:,I) = ZU_PROFILER(:) * COS(ZGAM) + ZV_PROFILER(:) * SIN(ZGAM)
      TPROFILER%MER (IN,:,I) = - ZU_PROFILER(:) * SIN(ZGAM) + ZV_PROFILER(:) * COS(ZGAM)
      TPROFILER%FF  (IN,:,I) = ZFF(:)                                                     
      TPROFILER%DD  (IN,:,I) = ZDD(:)   
      TPROFILER%W   (IN,:,I) = PROFILER_INTERP(PW)
      TPROFILER%TH  (IN,:,I) = PROFILER_INTERP(PTH)
      TPROFILER%THV (IN,:,I) = PROFILER_INTERP(ZTHV)
      TPROFILER%VISI(IN,:,I) = PROFILER_INTERP(ZVISI)
      TPROFILER%VISIKUN(IN,:,I) = PROFILER_INTERP(ZVISIKUN)
      TPROFILER%ZZ  (IN,:,I) = ZZ(:)
      TPROFILER%RHOD(IN,:,I) = ZRHOD(:)
      TPROFILER%CIZ(IN,:,I) = PROFILER_INTERP(PCIT)
! add RARE
        ! initialization CRARE and CRARE_ATT + LWC and IWC
      TPROFILER%CRARE(IN,:,I) = 0.
      TPROFILER%CRARE_ATT(IN,:,I) = 0.
      TPROFILER%LWCZ  (IN,:,I) = 0.
      TPROFILER%IWCZ  (IN,:,I) = 0.
      IF (CCLOUD=="LIMA" .OR. CCLOUD=="ICE3" ) THEN ! only for ICE3 and LIMA
       TPROFILER%LWCZ  (IN,:,I) = PROFILER_INTERP((PR(:,:,:,2)+PR(:,:,:,3))*PRHODREF(:,:,:))
       TPROFILER%IWCZ  (IN,:,I) = PROFILER_INTERP((PR(:,:,:,4)+PR(:,:,:,5)+PR(:,:,:,6))*PRHODREF(:,:,:))
       ZTEMPZ(:)=PROFILER_INTERP(ZTEMP(:,:,:))
       ZRHODREFZ(:)=PROFILER_INTERP(PRHODREF(:,:,:))
       ZCIT(:)=PROFILER_INTERP(PCIT(:,:,:))
       IF (CCLOUD=="LIMA") THEN
          ZCCI(:)=PROFILER_INTERP(PSV(:,:,:,NSV_LIMA_NI))
          ZCCR(:)=PROFILER_INTERP(PSV(:,:,:,NSV_LIMA_NR))
          ZCCC(:)=PROFILER_INTERP(PSV(:,:,:,NSV_LIMA_NC))
       ENDIF
       DO JLOOP=3,6
          ZRZ(:,JLOOP)=PROFILER_INTERP(PR(:,:,:,JLOOP))
       END DO
       IF (CSURF=="EXTE") THEN
         DO JK=1,IKU
            ZRZ(JK,2)=PROFILER_INTERP_2D(PR(:,:,JK,2)*PSEA(:,:))       ! becomes cloud mixing ratio over sea
            ZRZ(JK,7)=PROFILER_INTERP_2D(PR(:,:,JK,2)*(1.-PSEA(:,:)))  ! becomes cloud mixing ratio over land
         END DO
       ELSE
          ZRZ(:,2)=PROFILER_INTERP(PR(:,:,:,2))
          ZRZ(:,7)=0.
       END IF
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
        ! compute cloud radar reflectivity from vertical profiles of temperature
        ! and mixing ratios
       DO JK=1,IKU
         QMW=SQRT(QEPSW(ZTEMPZ(JK),XLIGHTSPEED/XLAM_CRAD))
         QMI=SQRT(QEPSI(ZTEMPZ(JK),XLIGHTSPEED/XLAM_CRAD))
         DO JLOOP=2,7
           IF (CCLOUD == 'LIMA') THEN
              GCALC=(ZRZ(JK,JLOOP)>ZRTMIN(JLOOP).AND.(JLOOP.NE.4.OR.ZCCI(JK)>0.).AND.&
                    (JLOOP.NE.3.OR.ZCCR(JK)>0.).AND.((JLOOP.NE.2.AND.JLOOP.NE.7).OR.ZCCC(JK)>0.))
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
                  !If temperature between -10 and 10B0C and Mr and Mg over min
                  !threshold: melting graupel
                  ! with liquid water fraction Fw=Mr/(Mr+Mg) else dry graupel
                  ! (Fw=0)    
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
             TPROFILER%CRARE(IN,JK,I)=TPROFILER%CRARE(IN,JK,I)+ZREFLOC
             ZAELOC(JK)=ZAELOC(JK)+ZAETMP
           END IF
         END DO
       END DO
     ! apply attenuation
       ALLOCATE(ZZMZ(IKU))
       ZZMZ = ZZ(:) ! PROFILER_INTERP(ZZM(:,:,:))
!        ZZMZ(1)=ZZM_STAT
        ! zenith
       ZAETOT=1.
       DO JK = 2,IKU
            ! attenuation from ZAELOC(JK) and ZAELOC(JK-1)
        ZAETOT=ZAETOT*EXP(-(ZAELOC(JK-1)+ZAELOC(JK))*(ZZMZ(JK)-ZZMZ(JK-1)))
        TPROFILER%CRARE_ATT(IN,JK,I)=TPROFILER%CRARE(IN,JK,I)*ZAETOT
       END DO
!        TPROFILER%ZZ  (IN,:,I) = ZZMZ(:)
       DEALLOCATE(ZZMZ,ZAELOC)
        ! m^3 bmm^6/m^3 bdBZ
       WHERE(TPROFILER%CRARE(IN,:,I)>0)
          TPROFILER%CRARE(IN,:,I)=10.*LOG10(1.E18*TPROFILER%CRARE(IN,:,I))
       ELSEWHERE
          TPROFILER%CRARE(IN,:,I)=XUNDEF
       END WHERE
       WHERE(TPROFILER%CRARE_ATT(IN,:,I)>0)
          TPROFILER%CRARE_ATT(IN,:,I)=10.*LOG10(1.E18*TPROFILER%CRARE_ATT(IN,:,I))
       ELSEWHERE
          TPROFILER%CRARE_ATT(IN,:,I)=XUNDEF
       END WHERE
       DEALLOCATE(ZX,ZW,ZRTMIN)
     END IF ! end LOOP ICE3
! end add RARE
!!
      IF (.NOT. L1D) THEN
        TPROFILER%P   (IN,:,I) = PROFILER_INTERP(PP(II(I):II(I)+1,IJ(I):IJ(I)+1,:))
      ELSE
        TPROFILER%P   (IN,:,I) = PROFILER_INTERP(PP)
      END IF
      !
      DO JSV=1,SIZE(PR,4)
        TPROFILER%R   (IN,:,I,JSV) = PROFILER_INTERP(PR(:,:,:,JSV))
      END DO
        ZWORK(:,:,:,:)=PSV(:,:,:,:)
        ZWORK(:,:,1,:)=PSV(:,:,2,:)
      DO JSV=1,SIZE(PSV,4)
        TPROFILER%SV  (IN,:,I,JSV) = PROFILER_INTERP(ZWORK(:,:,:,JSV))
      END DO
      ZWORK2(:,:,:,:) = 0.
      DO JK=IKB,IKE
        IKRAD = JK - JPVEXT
        ZWORK2(:,:,JK,:)=PAER(:,:,IKRAD,:)
      ENDDO
      DO JSV=1,SIZE(PAER,4)
        TPROFILER%AER(IN,:,I,JSV) = PROFILER_INTERP(ZWORK2(:,:,:,JSV))
      ENDDO
      IF (SIZE(PTKE)>0) TPROFILER%TKE  (IN,:,I) = PROFILER_INTERP(PTKE)
      !
      IF (LDIAG_IN_RUN) THEN
        TPROFILER%T2M   (IN,I)     = PROFILER_INTERP_2D(XCURRENT_T2M   )
        TPROFILER%Q2M   (IN,I)     = PROFILER_INTERP_2D(XCURRENT_Q2M   )
        TPROFILER%HU2M  (IN,I)     = PROFILER_INTERP_2D(XCURRENT_HU2M  )
        TPROFILER%ZON10M(IN,I)     = PROFILER_INTERP_2D(XCURRENT_ZON10M)
        TPROFILER%MER10M(IN,I)     = PROFILER_INTERP_2D(XCURRENT_MER10M)
        TPROFILER%RN    (IN,I)     = PROFILER_INTERP_2D(XCURRENT_RN    )
        TPROFILER%H     (IN,I)     = PROFILER_INTERP_2D(XCURRENT_H     )
        TPROFILER%LE    (IN,I)     = PROFILER_INTERP_2D(XCURRENT_LE    )
        TPROFILER%LEI   (IN,I)     = PROFILER_INTERP_2D(XCURRENT_LEI   )        
        TPROFILER%GFLUX (IN,I)     = PROFILER_INTERP_2D(XCURRENT_GFLUX )
       IF (CRAD /= 'NONE') THEN
        TPROFILER%SWD   (IN,I)     = PROFILER_INTERP_2D(XCURRENT_SWD   )
        TPROFILER%SWU   (IN,I)     = PROFILER_INTERP_2D(XCURRENT_SWU   )
        TPROFILER%LWD   (IN,I)     = PROFILER_INTERP_2D(XCURRENT_LWD   )
        TPROFILER%LWU   (IN,I)     = PROFILER_INTERP_2D(XCURRENT_LWU   )
       END IF
        TPROFILER%TKE_DISS(IN,:,I) = PROFILER_INTERP(XCURRENT_TKE_DISS)
      ENDIF
    ENDIF
!
!----------------------------------------------------------------------------
!
!*     11.   EXCHANGE OF INFORMATION BETWEEN PROCESSORS
!            ------------------------------------------
!
!*     11.2  data stored
!            -----------
!
  CALL DISTRIBUTE_PROFILER(TPROFILER%X   (I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%Y   (I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%LON (I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%LAT (I))
  !
  IF (LDIAG_IN_RUN) THEN
    CALL DISTRIBUTE_PROFILER(TPROFILER%T2M   (IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%Q2M   (IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%HU2M  (IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%ZON10M(IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%MER10M(IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%RN    (IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%H     (IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%LE    (IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%LEI   (IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%GFLUX (IN,I))
   IF (CRAD /= 'NONE') THEN
    CALL DISTRIBUTE_PROFILER(TPROFILER%LWD   (IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%LWU   (IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%SWD   (IN,I))
    CALL DISTRIBUTE_PROFILER(TPROFILER%SWU   (IN,I))
   ENDIF
  ENDIF
 DO JK=1,IKU
  CALL DISTRIBUTE_PROFILER(TPROFILER%ZON (IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%MER (IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%FF  (IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%DD  (IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%W   (IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%P   (IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%ZZ  (IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%TH  (IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%THV (IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%VISI(IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%VISIKUN(IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%RHOD(IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%CRARE(IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%CRARE_ATT(IN,JK,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%CIZ(IN,JK,I))

  !
  IF (LDIAG_IN_RUN) CALL DISTRIBUTE_PROFILER(TPROFILER%TKE_DISS(IN,JK,I))
  !
  DO JSV=1,SIZE(PR,4)
    CALL DISTRIBUTE_PROFILER(TPROFILER%R   (IN,JK,I,JSV))
  END DO
  DO JSV=1,SIZE(PSV,4)
    CALL DISTRIBUTE_PROFILER(TPROFILER%SV  (IN,JK,I,JSV))
  END DO
  IF (SIZE(PTKE)>0) CALL DISTRIBUTE_PROFILER(TPROFILER%TKE  (IN,JK,I))
 ENDDO

  CALL DISTRIBUTE_PROFILER(TPROFILER%IWV(IN,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%ZTD(IN,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%ZHD(IN,I))
  CALL DISTRIBUTE_PROFILER(TPROFILER%ZWD(IN,I))
ENDDO
!
END IF
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
CONTAINS
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
FUNCTION PROFILER_INTERP_2D(PA) RESULT(PB)
!
REAL, DIMENSION(:,:), INTENT(IN) :: PA
REAL                             :: PB
!
INTEGER :: JI, JJ
!
IF (SIZE(PA,1)==2) THEN
  JI=1
  JJ=1
ELSEIF (L1D) THEN
     JI=2
     JJ=2
ELSE
  JI=II(I)
  JJ=IJ(I)
END IF
!
!
PB = (1.-ZYCOEF(I)) * (1.-ZXCOEF(I)) *  PA(JI,JJ)    + &
     (1.-ZYCOEF(I)) *    (ZXCOEF(I)) *  PA(JI+1,JJ)  + &
     (   ZYCOEF(I)) * (1.-ZXCOEF(I)) *  PA(JI,JJ+1)  + &
     (   ZYCOEF(I)) *    (ZXCOEF(I)) *  PA(JI+1,JJ+1)
!
END FUNCTION PROFILER_INTERP_2D
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
FUNCTION PROFILER_INTERP(PA) RESULT(PB)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA
REAL, DIMENSION(SIZE(PA,3))        :: PB
!
INTEGER :: JI, JJ,JK
!
IF (SIZE(PA,1)==2) THEN
  JI=1
  JJ=1
ELSEIF (L1D) THEN
     JI=2
     JJ=2
ELSE
  JI=II(I)
  JJ=IJ(I)
END IF
!
!
DO JK=1,SIZE(PA,3)
 IF ( (PA(JI,JJ,JK) /= XUNDEF) .AND. (PA(JI+1,JJ,JK) /= XUNDEF) .AND. &
      (PA(JI,JJ+1,JK) /= XUNDEF) .AND. (PA(JI+1,JJ+1,JK) /= XUNDEF) ) THEN
    PB(JK) = (1.-ZYCOEF(I)) * (1.-ZXCOEF(I)) *  PA(JI,JJ,JK) + &
            (1.-ZYCOEF(I)) * (ZXCOEF(I)) *  PA(JI+1,JJ,JK)  + &
            (ZYCOEF(I)) * (1.-ZXCOEF(I)) *  PA(JI,JJ+1,JK)  + &
            (ZYCOEF(I)) * (ZXCOEF(I)) *  PA(JI+1,JJ+1,JK) 
 ELSE
    PB(JK) = XUNDEF 
 END IF
END DO
!
END FUNCTION PROFILER_INTERP
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
FUNCTION PROFILER_INTERP_U(PA) RESULT(PB)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA
REAL, DIMENSION(SIZE(PA,3))        :: PB
!
INTEGER :: JI, JJ
!
IF (SIZE(PA,1)==2) THEN
  JI=1
  JJ=1
ELSEIF (L1D) THEN
     JI=2
     JJ=2
ELSE
  JI=IU(I)
  JJ=IJ(I)
END IF
!
PB(:) = (1.- ZYCOEF(I)) * (1.-ZUCOEF(I)) * PA(JI  ,JJ  ,:) &
   +    (1.- ZYCOEF(I)) * (   ZUCOEF(I)) * PA(JI+1,JJ  ,:) &
   + (    ZYCOEF(I)) * (1.-ZUCOEF(I)) *    PA(JI  ,JJ+1,:) &
   + (    ZYCOEF(I)) * (   ZUCOEF(I)) *    PA(JI+1,JJ+1,:)
!
END FUNCTION PROFILER_INTERP_U
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
FUNCTION PROFILER_INTERP_V(PA) RESULT(PB)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA
REAL, DIMENSION(SIZE(PA,3))        :: PB
!
INTEGER :: JI, JJ
!
IF (SIZE(PA,1)==2) THEN
  JI=1
  JJ=1
ELSEIF (L1D) THEN
     JI=2
     JJ=2
ELSE
  JI=II(I)
  JJ=IV(I)
END IF
!
PB(:) = (1.- ZVCOEF(I)) * (1.-ZXCOEF(I)) * PA(JI  ,JJ  ,:) &
   + (1.- ZVCOEF(I)) * (   ZXCOEF(I)) *  PA(JI+1,JJ  ,:)   &
   + (    ZVCOEF(I)) * (1.-ZXCOEF(I)) *  PA(JI  ,JJ+1,:)   &
   + (    ZVCOEF(I)) * (   ZXCOEF(I)) *  PA(JI+1,JJ+1,:) 
!
END FUNCTION PROFILER_INTERP_V
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE DISTRIBUTE_PROFILER(PAS)
!
REAL, INTENT(INOUT) :: PAS
!
PAS = PAS * ZTHIS_PROCS(I)

CALL REDUCESUM_ll(PAS,IINFO_ll)
!
END SUBROUTINE DISTRIBUTE_PROFILER
!----------------------------------------------------------------------------
!
END SUBROUTINE PROFILER_n
