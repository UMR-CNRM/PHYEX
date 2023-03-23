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
      SUBROUTINE PROFILER_n( PZ, PRHODREF,                   &
                             PU, PV, PW, PTH, PR, PSV, PTKE, &
                             PTS, PP, PAER, PCIT, PSEA )
!
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
!     ########################################################
      SUBROUTINE PROFILER_n( PZ, PRHODREF,                   &
                             PU, PV, PW, PTH, PR, PSV, PTKE, &
                             PTS, PP, PAER, PCIT, PSEA )
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
!  P. Wautelet    04/2022: restructure profilers for better performance, reduce memory usage and correct some problems/bugs
! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,              ONLY: XCPD, XG, XLAM_CRAD, XLIGHTSPEED, XP00, XPI, XRD, XRHOLW, XRV, XTT
USE MODD_DIAG_IN_RUN
USE MODD_GRID,             ONLY: XBETA, XLON0, XRPK
USE MODD_NSV,              ONLY: NSV_C2R2, NSV_C2R2BEG, NSV_LIMA_NC, NSV_LIMA_NI, NSV_LIMA_NR
USE MODD_PARAMETERS,       ONLY: JPVEXT, XUNDEF
USE MODD_PARAM_ICE,        ONLY: LSNOW_T_I => LSNOW_T
USE MODD_PARAM_LIMA,       ONLY: LSNOW_T_L => LSNOW_T,                                                       &
                                 XALPHAR_L => XALPHAR, XNUR_L => XNUR, XALPHAS_L => XALPHAS, XNUS_L => XNUS, &
                                 XALPHAG_L => XALPHAG, XNUG_L => XNUG, XALPHAI_L => XALPHAI, XNUI_L => XNUI, &
                                 XRTMIN_L => XRTMIN, XALPHAC_L => XALPHAC, XNUC_L => XNUC
USE MODD_PARAM_LIMA_COLD,  ONLY: XDI_L => XDI, XLBEXI_L => XLBEXI, XLBI_L => XLBI, XAI_L => XAI, XBI_L => XBI, XC_I_L => XC_I, &
                                 XLBEXS_L => XLBEXS, XLBS_L => XLBS, XCCS_L => XCCS,                                           &
                                 XAS_L => XAS, XBS_L => XBS, XCXS_L => XCXS,                                                   &
                                 XLBDAS_MAX_L => XLBDAS_MAX, XLBDAS_MIN_L => XLBDAS_MIN,                                       &
                                 XNS_L => XNS, XTRANS_MP_GAMMAS_L=>XTRANS_MP_GAMMAS
USE MODD_PARAM_LIMA_MIXED, ONLY: XDG_L => XDG, XLBEXG_L => XLBEXG, XLBG_L => XLBG, XCCG_L => XCCG, &
                                 XAG_L => XAG, XBG_L => XBG, XCXG_L => XCXG, XCG_L => XCG
USE MODD_PARAM_LIMA_WARM,  ONLY: XLBEXR_L => XLBEXR, XLBR_L => XLBR, XBR_L => XBR, XAR_L => XAR, &
                                 XBC_L => XBC, XAC_L => XAC
USE MODD_PARAM_n,          ONLY: CCLOUD, CRAD, CSURF
USE MODD_PROFILER_n
USE MODD_RAIN_ICE_DESCR,   ONLY: XALPHAR_I => XALPHAR, XNUR_I => XNUR, XLBEXR_I => XLBEXR,                                 &
                                 XLBR_I => XLBR, XCCR_I => XCCR, XBR_I => XBR, XAR_I => XAR,                               &
                                 XALPHAC_I => XALPHAC, XNUC_I => XNUC,                                                     &
                                 XLBC_I => XLBC, XBC_I => XBC, XAC_I => XAC,                                               &
                                 XALPHAC2_I => XALPHAC2, XNUC2_I => XNUC2,                                                 &
                                 XALPHAS_I => XALPHAS, XNUS_I => XNUS, XLBEXS_I => XLBEXS,                                 &
                                 XLBS_I => XLBS, XCCS_I => XCCS, XAS_I => XAS, XBS_I => XBS, XCXS_I => XCXS,               &
                                 XALPHAG_I => XALPHAG, XNUG_I => XNUG, XDG_I => XDG, XLBEXG_I => XLBEXG,                   &
                                 XLBG_I => XLBG, XCCG_I => XCCG, XAG_I => XAG, XBG_I => XBG, XCXG_I => XCXG, XCG_I => XCG, &
                                 XALPHAI_I => XALPHAI, XNUI_I => XNUI, XDI_I => XDI, XLBEXI_I => XLBEXI,                   &
                                 XLBI_I => XLBI, XAI_I => XAI, XBI_I => XBI, XC_I_I => XC_I,                               &
                                 XNS_I => XNS, XRTMIN_I => XRTMIN, XCONC_LAND, XCONC_SEA,                                  &
                                 XLBDAS_MAX_I => XLBDAS_MAX, XLBDAS_MIN_I => XLBDAS_MIN,                                   &
                                 XTRANS_MP_GAMMAS_I => XTRANS_MP_GAMMAS
!
USE MODE_FGAU,             ONLY: GAULAG
USE MODE_FSCATTER,         ONLY: BHMIE, QEPSI, QEPSW, MG, MOMG
USE MODE_MSG
USE MODE_STATPROF_TOOLS,   ONLY: STATPROF_INSTANT, STATPROF_INTERP_2D, STATPROF_INTERP_3D, &
                                 STATPROF_INTERP_3D_U, STATPROF_INTERP_3D_V
!
USE MODI_GPS_ZENITH_GRID
USE MODI_WATER_SUM
!
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
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
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PCIT   ! ice concentration
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN)     :: PSEA   ! for radar
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!
INTEGER, PARAMETER :: JPTS_GAULAG = 9 ! number of points for Gauss-Laguerre quadrature
!
INTEGER :: IKB
INTEGER :: IKE
INTEGER :: IKU
!
!
REAL, DIMENSION(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3),SIZE(PSV,4))  :: ZWORK  
REAL, DIMENSION(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3),SIZE(PAER,4))  :: ZWORK2  
!
INTEGER :: IN     ! time index
INTEGER :: JSV    ! loop counter
INTEGER :: JK     ! loop
INTEGER :: JP     ! loop for profilers
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
REAL,DIMENSION(SIZE(PTH,1),SIZE(PTH,2))              :: ZZTD,ZZHD,ZZWD
REAL,DIMENSION(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3))  :: ZTEMP,ZTHV,ZTEMPV
REAL,DIMENSION(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3))  :: ZVISIGUL, ZVISIKUN
REAL ::  ZK1,ZK2,ZK3            ! k1, k2 and K3 atmospheric refractivity constants
REAL  :: ZRDSRV                 ! XRD/XRV
!
! specific to cloud radar
INTEGER                        :: JLOOP    ! loop counter
REAL, DIMENSION(SIZE(PR,3))    :: ZTEMPZ! vertical profile of temperature
REAL, DIMENSION(SIZE(PR,3))    :: ZRHODREFZ ! vertical profile of dry air density of the reference state
REAL, DIMENSION(SIZE(PR,3))    :: ZCIT     ! pristine ice concentration
REAL, DIMENSION(SIZE(PR,3))    :: ZCCI,ZCCR,ZCCC     ! ICE,RAIN CLOUD concentration (LIMA)
REAL, DIMENSION(SIZE(PR,3),SIZE(PR,4)+1) :: ZRZ  ! vertical profile of hydrometeor mixing ratios
REAL                           :: ZA, ZB, ZCC, ZCX, ZALPHA, ZNS, ZNU, ZLB, ZLBEX, ZRHOHYD   ! generic microphysical parameters
INTEGER                        :: JJ    ! loop counter for quadrature
COMPLEX                        :: QMW,QMI,QM,QB,QEPSIW,QEPSWI   ! dielectric parameter
REAL                           :: ZAETOT,ZAETMP,ZREFLOC,ZQSCA,ZQBACK,ZQEXT ! temporary scattering parameters
REAL,DIMENSION(:),ALLOCATABLE  :: ZAELOC,ZZMZ ! temporary arrays
REAL                           :: ZLBDA   ! slope distribution parameter
REAL                           :: ZDELTA_EQUIV ! mass-equivalent Gauss-Laguerre point
REAL                           :: ZFW ! liquid fraction
REAL                           :: ZFPW ! weight for mixed-phase reflectivity
REAL                           :: ZN   ! number concentration
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
!*      2.1  Indices
!            -------
!
IKU =   SIZE(PZ,3)     ! nombre de niveaux sur la verticale
IKB = JPVEXT+1
IKE = IKU-JPVEXT
!
!----------------------------------------------------------------------------
!
!*      3.4  instant of storage
!            ------------------
!
CALL  STATPROF_INSTANT( TPROFILERS_TIME, IN )
IF ( IN < 1 ) RETURN !No profiler storage at this time step
!
!----------------------------------------------------------------------------
!
!*      8.   DATA RECORDING
!            --------------
!
!PW: TODO: ne faire le calcul que si necessaire (presence de profileurs locaux,...)
ZTEMP(:,:,:)=PTH(:,:,:)*(PP(:,:,:)/ XP00) **(XRD/XCPD)
! Theta_v
ZTHV(:,:,:) = PTH(:,:,:) / (1.+WATER_SUM(PR(:,:,:,:)))*(1.+PR(:,:,:,1)/ZRDSRV)
! virtual temperature
ZTEMPV(:,:,:)=ZTHV(:,:,:)*(PP(:,:,:)/ XP00) **(XRD/XCPD)
CALL GPS_ZENITH_GRID(PR(:,:,:,1),ZTEMP,PP,ZZTD,ZZHD,ZZWD)

IF ( CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO' ) THEN
  ! Gultepe formulation
  ZVISIGUL(:,:,:) = 10E5 !default value
  WHERE ( (PR(:,:,:,2) /=0. ) .AND. (PSV(:,:,:,NSV_C2R2BEG+1) /=0. ) )
    ZVISIGUL(:,:,:) =1.002/(PR(:,:,:,2)*PRHODREF(:,:,:)*PSV(:,:,:,NSV_C2R2BEG+1))**0.6473
  END WHERE
END IF

IF ( CCLOUD /= 'NONE' .AND. CCLOUD /= 'REVE' ) THEN
  ! Kunkel formulation
  ZVISIKUN(:,:,:) = 10E5  !default value
  WHERE ( PR(:,:,:,2) /=0 )
    ZVISIKUN(:,:,:) =0.027/(10**(-8)+(PR(:,:,:,2)/(1+PR(:,:,:,2))*PRHODREF(:,:,:)*1000))**0.88
  END WHERE
END IF
!
PROFILER: DO JP = 1, NUMBPROFILER_LOC
  ZZ(:)                  = STATPROF_INTERP_3D( TPROFILERS(JP), PZ )
  ZRHOD(:)               = STATPROF_INTERP_3D( TPROFILERS(JP), PRHODREF )
  ZPRES(:)               = STATPROF_INTERP_3D( TPROFILERS(JP), PP )
  ZU_PROFILER(:)         = STATPROF_INTERP_3D_U( TPROFILERS(JP), PU )
  ZV_PROFILER(:)         = STATPROF_INTERP_3D_V( TPROFILERS(JP), PV )
  ZGAM                   = (XRPK * (TPROFILERS(JP)%XLON - XLON0) - XBETA)*(XPI/180.)
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
  XZS_GPS=TPROFILERS(JP)%XZ
  IF ( ABS( ZZ(IKB)-XZS_GPS ) < 150 ) THEN ! distance between real and model orography ok
    ZRV(:)                 = STATPROF_INTERP_3D( TPROFILERS(JP), PR(:,:,:,1) )
    ZT(:)                  = STATPROF_INTERP_3D( TPROFILERS(JP), ZTEMP )
    ZE(:)                  = ZPRES(:)*ZRV(:)/(ZRDSRV+ZRV(:))
    ZTV(:)                 = STATPROF_INTERP_3D( TPROFILERS(JP), ZTEMPV )
    ZZTD_PROFILER          = STATPROF_INTERP_2D( TPROFILERS(JP), ZZTD )
    ZZHD_PROFILER          = STATPROF_INTERP_2D( TPROFILERS(JP), ZZHD )
    ZZWD_PROFILER          = STATPROF_INTERP_2D( TPROFILERS(JP), ZZWD )
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
    TPROFILERS(JP)%XIWV(IN)= ZIWV
    TPROFILERS(JP)%XZTD(IN)= ZZTD_PROFILER
    TPROFILERS(JP)%XZWD(IN)= ZZWD_PROFILER
    TPROFILERS(JP)%XZHD(IN)= ZZHD_PROFILER
  ELSE
    CMNHMSG(1) = 'altitude of profiler ' // TRIM( TPROFILERS(JP)%CNAME ) // ' is too far from orography'
    CMNHMSG(2) = 'some variables are therefore not computed (IWV, ZTD, ZWD, ZHD)'
    CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'PROFILER_n', OLOCAL = .TRUE. )
    TPROFILERS(JP)%XIWV(IN)= XUNDEF
    TPROFILERS(JP)%XZTD(IN)= XUNDEF
    TPROFILERS(JP)%XZWD(IN)= XUNDEF
    TPROFILERS(JP)%XZHD(IN)= XUNDEF
  END IF
  TPROFILERS(JP)%XZON (IN,:) = ZU_PROFILER(:) * COS(ZGAM) + ZV_PROFILER(:) * SIN(ZGAM)
  TPROFILERS(JP)%XMER (IN,:) = - ZU_PROFILER(:) * SIN(ZGAM) + ZV_PROFILER(:) * COS(ZGAM)
  TPROFILERS(JP)%XFF  (IN,:) = ZFF(:)
  TPROFILERS(JP)%XDD  (IN,:) = ZDD(:)
  TPROFILERS(JP)%XW   (IN,:) = STATPROF_INTERP_3D( TPROFILERS(JP), PW )
  TPROFILERS(JP)%XTH  (IN,:) = STATPROF_INTERP_3D( TPROFILERS(JP), PTH )
  TPROFILERS(JP)%XTHV (IN,:) = STATPROF_INTERP_3D( TPROFILERS(JP), ZTHV )
  IF ( CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO' )  TPROFILERS(JP)%XVISIGUL(IN,:) = STATPROF_INTERP_3D( TPROFILERS(JP), ZVISIGUL )
  IF ( CCLOUD /= 'NONE' .AND. CCLOUD /= 'REVE' ) TPROFILERS(JP)%XVISIKUN(IN,:) = STATPROF_INTERP_3D( TPROFILERS(JP), ZVISIKUN )
  TPROFILERS(JP)%XZZ  (IN,:) = ZZ(:)
  TPROFILERS(JP)%XRHOD(IN,:) = ZRHOD(:)
  IF ( CCLOUD == 'ICE3' .OR. CCLOUD == 'ICE4' ) &
    TPROFILERS(JP)%XCIZ(IN,:) = STATPROF_INTERP_3D( TPROFILERS(JP), PCIT )
! add RARE
    ! initialization CRARE and CRARE_ATT + LWC and IWC
  TPROFILERS(JP)%XCRARE(IN,:) = 0.
  TPROFILERS(JP)%XCRARE_ATT(IN,:) = 0.
  TPROFILERS(JP)%XLWCZ  (IN,:) = 0.
  TPROFILERS(JP)%XIWCZ  (IN,:) = 0.
  IF (CCLOUD=="LIMA" .OR. CCLOUD=="ICE3" ) THEN ! only for ICE3 and LIMA
    TPROFILERS(JP)%XLWCZ  (IN,:) = STATPROF_INTERP_3D( TPROFILERS(JP), (PR(:,:,:,2)+PR(:,:,:,3))*PRHODREF(:,:,:) )
    TPROFILERS(JP)%XIWCZ  (IN,:) = STATPROF_INTERP_3D( TPROFILERS(JP), (PR(:,:,:,4)+PR(:,:,:,5)+PR(:,:,:,6))*PRHODREF(:,:,:) )
    ZTEMPZ(:)=STATPROF_INTERP_3D( TPROFILERS(JP), ZTEMP(:,:,:) )
    ZRHODREFZ(:)=STATPROF_INTERP_3D( TPROFILERS(JP), PRHODREF(:,:,:) )
    ZCIT(:)=STATPROF_INTERP_3D( TPROFILERS(JP), PCIT(:,:,:) )
    IF (CCLOUD=="LIMA") THEN
      ZCCI(:)=STATPROF_INTERP_3D( TPROFILERS(JP), PSV(:,:,:,NSV_LIMA_NI) )
      ZCCR(:)=STATPROF_INTERP_3D( TPROFILERS(JP), PSV(:,:,:,NSV_LIMA_NR) )
      ZCCC(:)=STATPROF_INTERP_3D( TPROFILERS(JP), PSV(:,:,:,NSV_LIMA_NC) )
    END IF
    DO JLOOP=3,6
      ZRZ(:,JLOOP)=STATPROF_INTERP_3D( TPROFILERS(JP), PR(:,:,:,JLOOP) )
    END DO
    IF (CSURF=="EXTE") THEN
      DO JK=1,IKU
        ZRZ(JK,2)=STATPROF_INTERP_2D( TPROFILERS(JP), PR(:,:,JK,2)*PSEA(:,:) )       ! becomes cloud mixing ratio over sea
        ZRZ(JK,7)=STATPROF_INTERP_2D( TPROFILERS(JP), PR(:,:,JK,2)*(1.-PSEA(:,:)) )  ! becomes cloud mixing ratio over land
      END DO
    ELSE
      ZRZ(:,2)=STATPROF_INTERP_3D( TPROFILERS(JP), PR(:,:,:,2) )
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
    END IF
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
        END IF
        IF (GCALC) THEN
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
              END IF
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
              END IF
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
              END IF
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
              END IF
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
              END IF
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
              END IF
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
              END IF
          END SELECT
          IF ( JLOOP ==  5 .AND. CCLOUD=='LIMA'.AND.LSNOW_T_L ) THEN
            IF (ZTEMPZ(JK)>XTT-10.) THEN
              ZLBDA = MAX(MIN(XLBDAS_MAX_L, 10**(14.554-0.0423*ZTEMPZ(JK))),XLBDAS_MIN_L)*XTRANS_MP_GAMMAS_L
            ELSE
              ZLBDA = MAX(MIN(XLBDAS_MAX_L, 10**(6.226-0.0106*ZTEMPZ(JK))),XLBDAS_MIN_L)*XTRANS_MP_GAMMAS_L
            END IF
            ZN=ZNS*ZRHODREFZ(JK)*ZRZ(JK,JLOOP)*ZLBDA**ZB
          ELSE IF (JLOOP.EQ.5 .AND. (CCLOUD=='ICE3'.AND.LSNOW_T_I) ) THEN
            IF (ZTEMPZ(JK)>XTT-10.) THEN
              ZLBDA = MAX(MIN(XLBDAS_MAX_I, 10**(14.554-0.0423*ZTEMPZ(JK))),XLBDAS_MIN_I)*XTRANS_MP_GAMMAS_I
            ELSE
              ZLBDA = MAX(MIN(XLBDAS_MAX_I, 10**(6.226-0.0106*ZTEMPZ(JK))),XLBDAS_MIN_I)*XTRANS_MP_GAMMAS_I
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
          TPROFILERS(JP)%XCRARE(IN,JK)=TPROFILERS(JP)%XCRARE(IN,JK)+ZREFLOC
          ZAELOC(JK)=ZAELOC(JK)+ZAETMP
        END IF
      END DO
    END DO
 ! apply attenuation
    ALLOCATE(ZZMZ(IKU))
    ZZMZ = ZZ(:) ! STATPROF_INTERP_3D( TPROFILERS(JP), ZZM(:,:,:) )
!        ZZMZ(1)=ZZM_STAT
    ! zenith
    ZAETOT=1.
    DO JK = 2,IKU
         ! attenuation from ZAELOC(JK) and ZAELOC(JK-1)
      ZAETOT=ZAETOT*EXP(-(ZAELOC(JK-1)+ZAELOC(JK))*(ZZMZ(JK)-ZZMZ(JK-1)))
      TPROFILERS(JP)%XCRARE_ATT(IN,JK)=TPROFILERS(JP)%XCRARE(IN,JK)*ZAETOT
    END DO
    DEALLOCATE(ZZMZ,ZAELOC)
     ! m^3 b mm^6/m^3 b dBZ
    WHERE(TPROFILERS(JP)%XCRARE(IN,:)>0)
      TPROFILERS(JP)%XCRARE(IN,:)=10.*LOG10(1.E18*TPROFILERS(JP)%XCRARE(IN,:))
    ELSEWHERE
      TPROFILERS(JP)%XCRARE(IN,:)=XUNDEF
    END WHERE
    WHERE(TPROFILERS(JP)%XCRARE_ATT(IN,:)>0)
      TPROFILERS(JP)%XCRARE_ATT(IN,:)=10.*LOG10(1.E18*TPROFILERS(JP)%XCRARE_ATT(IN,:))
    ELSEWHERE
      TPROFILERS(JP)%XCRARE_ATT(IN,:)=XUNDEF
    END WHERE
    DEALLOCATE(ZX,ZW,ZRTMIN)
  END IF ! end LOOP ICE3
! end add RARE
!!
  TPROFILERS(JP)%XP   (IN,:) = STATPROF_INTERP_3D( TPROFILERS(JP), PP )
  !
  DO JSV=1,SIZE(PR,4)
    TPROFILERS(JP)%XR   (IN,:,JSV) = STATPROF_INTERP_3D( TPROFILERS(JP), PR(:,:,:,JSV) )
  END DO
    ZWORK(:,:,:,:)=PSV(:,:,:,:)
    ZWORK(:,:,1,:)=PSV(:,:,2,:)
  DO JSV=1,SIZE(PSV,4)
    TPROFILERS(JP)%XSV  (IN,:,JSV) = STATPROF_INTERP_3D( TPROFILERS(JP), ZWORK(:,:,:,JSV) )
  END DO
  ZWORK2(:,:,:,:) = 0.
  DO JK=IKB,IKE
    IKRAD = JK - JPVEXT
    ZWORK2(:,:,JK,:)=PAER(:,:,IKRAD,:)
  END DO
  DO JSV=1,SIZE(PAER,4)
    TPROFILERS(JP)%XAER(IN,:,JSV) = STATPROF_INTERP_3D( TPROFILERS(JP), ZWORK2(:,:,:,JSV) )
  END DO
  IF (SIZE(PTKE)>0) TPROFILERS(JP)%XTKE  (IN,:) = STATPROF_INTERP_3D( TPROFILERS(JP), PTKE )
  !
  IF (LDIAG_IN_RUN) THEN
    TPROFILERS(JP)%XT2M   (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_T2M    )
    TPROFILERS(JP)%XQ2M   (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_Q2M    )
    TPROFILERS(JP)%XHU2M  (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_HU2M   )
    TPROFILERS(JP)%XZON10M(IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_ZON10M )
    TPROFILERS(JP)%XMER10M(IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_MER10M )
    TPROFILERS(JP)%XRN    (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_RN     )
    TPROFILERS(JP)%XH     (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_H      )
    TPROFILERS(JP)%XLE    (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_LE     )
    TPROFILERS(JP)%XLEI   (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_LEI    )
    TPROFILERS(JP)%XGFLUX (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_GFLUX  )
    IF (CRAD /= 'NONE') THEN
      TPROFILERS(JP)%XSWD   (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_SWD    )
      TPROFILERS(JP)%XSWU   (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_SWU    )
      TPROFILERS(JP)%XLWD   (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_LWD    )
      TPROFILERS(JP)%XLWU   (IN)     = STATPROF_INTERP_2D( TPROFILERS(JP), XCURRENT_LWU    )
    END IF
    TPROFILERS(JP)%XTKE_DISS(IN,:) = STATPROF_INTERP_3D( TPROFILERS(JP), XCURRENT_TKE_DISS )
  END IF
END DO PROFILER
!
!----------------------------------------------------------------------------
!
END SUBROUTINE PROFILER_n
