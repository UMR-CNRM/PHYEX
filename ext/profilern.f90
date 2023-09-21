!MNH_LIC Copyright 2002-2023 CNRS, Meteo-France and Universite Paul Sabatier
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
!  P. Wautelet 01/06/2023: deduplicate code => moved to modd/mode_sensors.f90
! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_ALLPROFILER_n,    ONLY: LDIAG_SURFRAD_PROF
USE MODD_CST,              ONLY: XCPD, XG, XP00, XPI, XRD, XRV
USE MODD_DIAG_IN_RUN,      ONLY: XCURRENT_TKE_DISS
USE MODD_GRID,             ONLY: XBETA, XLON0, XRPK
USE MODD_NSV,              ONLY: NSV_C2R2BEG, NSV_LIMA_NC, NSV_LIMA_NR, NSV_LIMA_NI
USE MODD_PARAMETERS,       ONLY: JPVEXT, XUNDEF
USE MODD_PARAM_n,          ONLY: CCLOUD, CRAD
USE MODD_PROFILER_n
!
USE MODE_FGAU,             ONLY: GAULAG
USE MODE_MSG
USE MODE_SENSOR,           ONLY: Sensor_rare_compute, Sensor_wc_compute
USE MODE_STATPROF_TOOLS,   ONLY: STATPROF_DIAG_SURFRAD
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
INTEGER :: IKB
INTEGER :: IKE
INTEGER :: IKU
!
!
REAL, DIMENSION(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3),SIZE(PSV,4))  :: ZWORK
REAL, DIMENSION(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3),SIZE(PAER,4)) :: ZWORK2
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
IF ( .NOT. TPROFILERS_TIME%STORESTEP_CHECK_AND_SET( IN ) ) RETURN !No profiler storage at this time step
!
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
  TPROFILERS(JP)%NSTORE_CUR = IN

  ZZ(:)                  = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PZ )
  ZRHOD(:)               = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PRHODREF )
  ZPRES(:)               = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PP )
  ZU_PROFILER(:)         = TPROFILERS(JP)%INTERP_HOR_FROM_UPOINT( PU )
  ZV_PROFILER(:)         = TPROFILERS(JP)%INTERP_HOR_FROM_VPOINT( PV )
  ZGAM                   = (XRPK * (TPROFILERS(JP)%XLON_CUR - XLON0) - XBETA)*(XPI/180.)
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
  XZS_GPS=TPROFILERS(JP)%XZ_CUR
  IF ( ABS( ZZ(IKB)-XZS_GPS ) < 150 ) THEN ! distance between real and model orography ok
    ZRV(:)                 = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PR(:,:,:,1) )
    ZT(:)                  = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( ZTEMP )
    ZE(:)                  = ZPRES(:)*ZRV(:)/(ZRDSRV+ZRV(:))
    ZTV(:)                 = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( ZTEMPV )
    ZZTD_PROFILER          = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( ZZTD )
    ZZHD_PROFILER          = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( ZZHD )
    ZZWD_PROFILER          = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( ZZWD )
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
  TPROFILERS(JP)%XZON (:,IN) = ZU_PROFILER(:) * COS(ZGAM) + ZV_PROFILER(:) * SIN(ZGAM)
  TPROFILERS(JP)%XMER (:,IN) = - ZU_PROFILER(:) * SIN(ZGAM) + ZV_PROFILER(:) * COS(ZGAM)
  TPROFILERS(JP)%XFF  (:,IN) = ZFF(:)
  TPROFILERS(JP)%XDD  (:,IN) = ZDD(:)
  TPROFILERS(JP)%XW   (:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PW )
  TPROFILERS(JP)%XTH  (:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PTH )
  TPROFILERS(JP)%XTHV (:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( ZTHV )
  IF ( CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO' ) &
    TPROFILERS(JP)%XVISIGUL(:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( ZVISIGUL )
  IF ( CCLOUD /= 'NONE' .AND. CCLOUD /= 'REVE' ) &
    TPROFILERS(JP)%XVISIKUN(:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( ZVISIKUN )
  TPROFILERS(JP)%XZZ  (:,IN) = ZZ(:)
  TPROFILERS(JP)%XRHOD(:,IN) = ZRHOD(:)
  IF (CCLOUD=="LIMA") THEN
    TPROFILERS(JP)%XCIZ(:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PSV(:,:,:,NSV_LIMA_NI) )
    TPROFILERS(JP)%XCCZ(:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PSV(:,:,:,NSV_LIMA_NC) )
    TPROFILERS(JP)%XCRZ(:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PSV(:,:,:,NSV_LIMA_NR) )
  ELSE IF ( CCLOUD=="ICE3" .OR. CCLOUD=="ICE4" ) THEN
    TPROFILERS(JP)%XCIZ(:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PCIT )
  END IF

  CALL Sensor_wc_compute(   TPROFILERS(JP), IN, PR, PRHODREF )
  CALL Sensor_rare_compute( TPROFILERS(JP), IN, PR, PSV, PRHODREF, PCIT, ZTEMP, ZZ, PSEA )
  !!
  TPROFILERS(JP)%XP   (:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PP )
  !
  DO JSV=1,SIZE(PR,4)
    TPROFILERS(JP)%XR   (:,IN,JSV) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PR(:,:,:,JSV) )
  END DO
    ZWORK(:,:,:,:)=PSV(:,:,:,:)
    ZWORK(:,:,1,:)=PSV(:,:,2,:)
  DO JSV=1,SIZE(PSV,4)
    TPROFILERS(JP)%XSV  (:,IN,JSV) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( ZWORK(:,:,:,JSV) )
  END DO
  ZWORK2(:,:,:,:) = 0.
  DO JK=IKB,IKE
    IKRAD = JK - JPVEXT
    ZWORK2(:,:,JK,:)=PAER(:,:,IKRAD,:)
  END DO
  DO JSV=1,SIZE(PAER,4)
    TPROFILERS(JP)%XAER(:,IN,JSV) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( ZWORK2(:,:,:,JSV) )
  END DO
  IF (SIZE(PTKE)>0) TPROFILERS(JP)%XTKE  (:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PTKE )

  ! XRHOD_SENSOR is not computed for profilers because not very useful
  ! If needed, the interpolation must also be done vertically
  ! (and therefore the vertical interpolation coefficients have to be computed)
  ! TPROFILERS(JP)%XRHOD_SENSOR(IN) = ...

  IF ( CRAD /= 'NONE' ) TPROFILERS(JP)%XTSRAD(IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( PTS )
  !
  IF ( LDIAG_SURFRAD_PROF ) CALL STATPROF_DIAG_SURFRAD(TPROFILERS(JP), IN )
  TPROFILERS(JP)%XTKE_DISS(:,IN) = TPROFILERS(JP)%INTERP_HOR_FROM_MASSPOINT( XCURRENT_TKE_DISS )
END DO PROFILER
!
!----------------------------------------------------------------------------
!
END SUBROUTINE PROFILER_n
