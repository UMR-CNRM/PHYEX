!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########
MODULE MODI_GROUND_PARAM_n
!     ##########
!
INTERFACE 
!
      SUBROUTINE GROUND_PARAM_n(D, PSFTH, PSFRV, PSFSV, PSFCO2, PSFU, PSFV, &
                                 PDIR_ALB, PSCA_ALB, PEMIS, PTSRAD, KTCOUNT, TPFILE )
!
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!* surface fluxes
!  --------------
!
USE MODD_IO,      ONLY: TFILEDATA
!
TYPE(DIMPHYEX_t),     INTENT(IN)   :: D
REAL, DIMENSION(:,:), INTENT(OUT) :: PSFTH ! surface flux of potential temperature (Km/s)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSFRV ! surface flux of water vapor           (m/s*kg/kg)
REAL, DIMENSION(:,:,:),INTENT(OUT):: PSFSV ! surface flux of scalar                (m/s*kg/kg)
                                           ! flux of chemical var.                 (ppp.m/s)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSFCO2! surface flux of CO2                   (m/s*kg/kg)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSFU  ! surface fluxes of horizontal   
REAL, DIMENSION(:,:), INTENT(OUT) :: PSFV  ! momentum in x and y directions        (m2/s2)
!
!* Radiative parameters
!  --------------------
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDIR_ALB  ! direct  albedo for each spectral band (-)
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each spectral band (-)
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEMIS     ! surface emissivity                    (-)
REAL, DIMENSION(:,:),   INTENT(OUT) :: PTSRAD    ! surface radiative temperature         (K)
!
INTEGER,                INTENT(IN)  :: KTCOUNT   ! temporal iteration count
TYPE(TFILEDATA),        INTENT(IN)  :: TPFILE    ! Synchronous output file
END SUBROUTINE GROUND_PARAM_n
!
END INTERFACE
!
END MODULE MODI_GROUND_PARAM_n
!
!     ######################################################################
      SUBROUTINE GROUND_PARAM_n(D, PSFTH, PSFRV, PSFSV, PSFCO2, PSFU, PSFV, &
                                 PDIR_ALB, PSCA_ALB, PEMIS, PTSRAD, KTCOUNT, TPFILE )
!     #######################################################################
!
!
!!****  *GROUND_PARAM*  
!!
!!    PURPOSE
!!    -------
!       Monitor to call the externalized surface
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    Noilhan and Planton (1989)
!!      
!!    AUTHOR
!!    ------
!!	S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/03/95 
!!      (J.Stein)   25/10/95  add the rain flux computation at the ground
!!                            and the lbc
!!      (J.Stein)   15/11/95  include the strong slopes cases
!!      (J.Stein)   06/02/96  bug correction for the precipitation flux writing 
!!      (J.Stein)   20/05/96  set the right IGRID value for the rain rate
!!      (J.Viviand) 04/02/97  add cold and convective precipitation rate
!!      (J.Stein)   22/06/97  use the absolute pressure    
!!      (V.Masson)  09/07/97  add directional z0 computations and RESA correction     
!!      (V.Masson)  13/02/98  merge the ISBA and TSZ0 routines,
!!                            rename the routine as a monitor, called by PHYS_PARAMn
!!                            add the town parameterization
!!                            recomputes z0 where snow is.
!!                            pack and unpack of 2D fields into 1D fields
!!      (V.Masson)  04/01/00  removes the TSZ0 case
!       (F.Solmon/V.Masson)   adapatation for patch approach
!                             modification of internal subroutine pack/ allocation in function 
!                              of patch indices
!                             calling of isba for each defined patch
!                             averaging of patch fluxes to get nat fluxes 
!       (P. Tulet/G.Guenais)  04/02/01  separation of vegetatives class
!                                           for friction velocity and 
!                                           aerodynamical resistance
!      (S Donnier)            09/12/02  add specific humidity at 2m for diagnostic
!      (V.Masson)             01/03/03  externalisation of the surface schemes!
!      (P.Tulet )             01/11/03  externalisation of the surface chemistry!
!!     (D.Gazen)              01/12/03  change emissions handling for surf. externalization
!!     (J.escobar)            18/10/2012 missing USE MODI_COUPLING_SURF_ATM_n & MODI_DIAG_SURF_ATM_n
!      (J.escobar)            02/2014 add Forefire coupling
!!     (G.Delautier)          06/2016 phasage surfex 8
!!     (B.Vie)                2016 LIMA
!!     (J.Pianezze)           08/2016 add send/recv oasis functions
!!      (M.Leriche)            24/03/16 remove flag for chemical surface fluxes
!!      (M.Leriche)           01/07/2017 Add DIAG chimical surface fluxes
!!  01/2018      (G.Delautier) SURFEX 8.1
!!                   02/2018 Q.Libois ECRAD
!!     (P.Wautelet) 28/03/2018 replace TEMPORAL_DIST by DATETIME_DISTANCE

!!     (V. Vionnet)           18/07/2017 add coupling for blowing snow module 
!!     (Bielli S.) 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  P. Wautelet 09/02/2022: bugfix: add missing XCURRENT_LEI computation
!  A. Costes      12/2021: Blaze Fire model
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
! 
#ifdef CPLOASIS
USE MODI_GET_HALO
USE MODI_MNH_OASIS_RECV
USE MODI_MNH_OASIS_SEND
USE MODD_SFX_OASIS, ONLY : LOASIS
USE MODD_DYN,       ONLY : XSEGLEN
USE MODD_DYN_n,     ONLY : DYN_MODEL
#endif
!
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_CST,        ONLY : XP00, XCPD, XRD, XRV,XRHOLW, XDAY, XPI, XLVTT, XMD, XAVOGADRO
USE MODD_DIMPHYEX,   ONLY : DIMPHYEX_t
USE MODD_PARAMETERS, ONLY : JPVEXT, XUNDEF
USE MODD_DYN_n,      ONLY : XTSTEP
USE MODD_CH_MNHC_n,  ONLY : LUSECHEM
USE MODD_CH_M9_n,   ONLY : CNAMES
USE MODD_FIELD_n,    ONLY : XUT, XVT, XWT, XTHT, XRT, XPABST, XSVT, XTKET, XZWS,&
XLSPHI, XBMAP, XFMR0, XFMRFA, XFMWF0, XFMR00, XFMIGNITION, XFMFUELTYPE,&
XFIRETAU, XFLUXPARAMH, XFLUXPARAMW, XFIRERW, XFMASE, XFMAWC, XFMWALKIG,&
XFMFLUXHDH, XFMFLUXHDW, XRTHS, XRRS, XFMHWS, XFMWINDU, XFMWINDV, XFMWINDW, XGRADLSPHIX, &
XGRADLSPHIY, XFIREWIND, XFMGRADOROX, XFMGRADOROY
USE MODD_METRICS_n,  ONLY : XDXX, XDYY, XDZZ
USE MODD_DIM_n,      ONLY : NKMAX
USE MODD_GRID_n,     ONLY : XLON, XZZ, XDIRCOSXW, XDIRCOSYW, XDIRCOSZW, &
                            XCOSSLOPE, XSINSLOPE, XZS
USE MODD_REF_n,      ONLY : XRHODREF,XRHODJ,XEXNREF
USE MODD_CONF_n,     ONLY : NRR
USE MODD_PARAM_n,    ONLY : CDCONV,CCLOUD, CRAD
USE MODD_PRECIP_n,   ONLY : XINPRC, XINPRR, XINPRS, XINPRG, XINPRH
USE MODD_DEEP_CONVECTION_n, ONLY : XPRCONV, XPRSCONV
USE MODD_CONF,       ONLY : LCARTESIAN, CPROGRAM
USE MODD_TIME_n,     ONLY : TDTCUR
USE MODD_RADIATIONS_n,   ONLY : XFLALWD, XCCO2, XTSIDER, &
                                XSW_BANDS, XDIRSRFSWD, XSCAFLASWD, &
                                XZENITH, XAZIM, XAER, XSWU, XLWU
USE MODD_NSV
USE MODD_GRID,       ONLY : XLON0, XRPK, XBETA
USE MODD_PARAM_ICE,  ONLY : LSEDIC
USE MODD_PARAM_C2R2, ONLY : LSEDC
USE MODD_DIAG_IN_RUN
USE MODD_DUST,       ONLY : LDUST 
USE MODD_SALT,       ONLY : LSALT
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
USE MODD_CH_AEROSOL, ONLY : LORILAM
USE MODD_CSTS_DUST,  ONLY : XMOLARWEIGHT_DUST
USE MODD_CSTS_SALT,  ONLY : XMOLARWEIGHT_SALT
USE MODD_CH_FLX_n, ONLY : XCHFLX
USE MODD_DIAG_FLAG, ONLY : LCHEMDIAG
!
USE MODI_NORMAL_INTERPOL
USE MODE_ROTATE_WIND, ONLY : ROTATE_WIND
USE MODI_SHUMAN
USE MODI_MNHGET_SURF_PARAM_n
USE MODI_COUPLING_SURF_ATM_n
USE MODI_DIAG_SURF_ATM_n
USE MODD_MNH_SURFEX_n
!
USE MODE_DATETIME
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
#ifdef MNH_FOREFIRE
!** MODULES FOR FOREFIRE **!
USE MODD_FOREFIRE
USE MODD_FOREFIRE_n
USE MODI_COUPLING_FOREFIRE_n
#endif
!
USE MODD_TIME_n
USE MODD_TIME
!
USE MODD_PARAM_LIMA, ONLY : MSEDC=>LSEDC
!
USE MODD_FIRE
USE MODD_FIELD
USE MODI_FIRE_MODEL
USE MODD_CONF, ONLY : NVERB, NHALO
USE MODE_MNH_TIMING, ONLY : SECOND_MNH2
USE MODE_MSG
USE MODD_IO,      ONLY: TFILEDATA
!
IMPLICIT NONE
!
!
!
!*      0.1    declarations of arguments
!
!* surface fluxes
!  --------------
!
TYPE(DIMPHYEX_t),     INTENT(IN)   :: D
REAL, DIMENSION(:,:), INTENT(OUT) :: PSFTH ! surface flux of potential temperature (Km/s)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSFRV ! surface flux of water vapor           (m/s*kg/kg)
REAL, DIMENSION(:,:,:),INTENT(OUT):: PSFSV ! surface flux of scalar                (m/s*kg/kg)
                                           ! flux of chemical var.                 (ppp.m/s)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSFCO2! surface flux of CO2                   (m/s*kg/kg)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSFU  ! surface fluxes of horizontal   
REAL, DIMENSION(:,:), INTENT(OUT) :: PSFV  ! momentum in x and y directions        (m2/s2)
!
!* Radiative parameters
!  --------------------
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDIR_ALB  ! direct  albedo for each spectral band (-)
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each spectral band (-)
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEMIS     ! surface emissivity                    (-)
REAL, DIMENSION(:,:),   INTENT(OUT) :: PTSRAD    ! surface radiative temperature         (K)
!
INTEGER,                INTENT(IN)  :: KTCOUNT   ! temporal iteration count
TYPE(TFILEDATA),        INTENT(IN)  :: TPFILE    ! Synchronous output file
!
!-------------------------------------------------------------------------------
!
!
!
!*      0.2    declarations of local variables
!              -------------------------------
!
!
!* Atmospheric variables
!  ---------------------
!
REAL, DIMENSION(:,:,:), ALLOCATABLE         :: ZRV    ! vapor mixing ratio
!
!            suffix 'A' stands for atmospheric variable at first model level
!
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZZREF  ! Forcing height
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZTA    ! Temperature
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZRVA   ! vapor mixing ratio
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZQA    ! humidity (kg/m3)
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZPA    ! Pressure
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZPS    ! Pressure
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZEXNA  ! Exner function
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZEXNS  ! Exner function
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZTHA   ! potential temperature
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZRAIN  ! liquid precipitation (kg/m2/s)
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZSNOW  ! solid precipitation  (kg/m2/s)
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZTSUN  ! solar time           (s since midnight)
!
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZUA    ! u component of the wind
!                                                       ! parallel to the orography
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZVA    ! v component of the wind
!                                                       ! parallel to the orography
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZU     ! zonal wind
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZV     ! meridian wind
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZWIND  ! wind parallel to the orography
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZRHOA  ! air density
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZDIR   ! wind direction (rad from N clockwise)
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZSFU   ! zonal momentum flux
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZSFV   ! meridian momentum flux
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZCO2   ! CO2 concentration (kg/kg)
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZALFA  ! angle between the wind
!                                                       ! and the x axis
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2),1):: ZU2D   ! u and v component of the
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2),1):: ZV2D   ! wind at mass point
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZSFTH  ! Turbulent flux of heat
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZSFTQ  ! Turbulent flux of water
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2))  :: ZSFCO2 ! Turbulent flux of CO2
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2),NSV):: ZSFTS! Turbulent flux of scalar
!
REAL, DIMENSION(SIZE(PSFTH,1),SIZE(PSFTH,2),NBLOWSNOW_2D)  :: ZBLOWSNOW_2D  ! 2D blowing snow variables
                                                                            ! after advection
                                  ! They refer to the 2D fields advected by MNH including:
                                  !             - total number concentration in Canopy
                                  !             - total mass concentration in Canopy
                                  !             - equivalent concentration in the saltation layer
!
!* Dimensions
!  ----------
!
INTEGER :: IIB  ! physical boundary
INTEGER :: IIE  ! physical boundary
INTEGER :: IJB  ! physical boundary
INTEGER :: IJE  ! physical boundary
INTEGER :: IKB  ! physical boundary
INTEGER :: IKE  ! physical boundary
INTEGER :: IKU  ! vertical array sizes
!
INTEGER :: JLAYER ! loop counter
INTEGER :: JSV    ! loop counter
INTEGER :: JI,JJ,JK ! loop index
!
INTEGER :: IDIM1 ! X physical dimension
INTEGER :: IDIM2 ! Y physical dimension
INTEGER :: IDIM1D! total physical dimension
INTEGER :: IKRAD
!
INTEGER :: KSV_SURF  ! Number of scalar variables sent to SURFEX
!
!* Arrays put in 1D vectors
!  ------------------------
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_TSUN     ! solar time
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_ZENITH   ! zenithal angle
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_AZIM     ! azimuthal angle
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_ZREF     ! forcing height
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_ZS       ! orography
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_U        ! zonal wind
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_V        ! meridian wind
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_QA       ! air humidity  (kg/m3)
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_TA       ! air temperature
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_RHOA     ! air density
REAL, DIMENSION(:,:), ALLOCATABLE :: ZP_SV       ! scalar at first atmospheric level
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_CO2      ! air CO2 concentration
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_RAIN     ! liquid precipitation
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_SNOW     ! solid precipitation
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_LW       ! incoming longwave
REAL, DIMENSION(:,:), ALLOCATABLE :: ZP_DIR_SW   ! direct incoming shortwave
REAL, DIMENSION(:,:), ALLOCATABLE :: ZP_SCA_SW   ! diffuse incoming shortwave
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_PS       ! surface pressure
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_PA       ! pressure at first atmospheric level
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_ZWS      ! significant wave height (m)

REAL, DIMENSION(:),   ALLOCATABLE :: ZP_SFTQ     ! water vapor flux
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_SFTH     ! potential temperature flux
REAL, DIMENSION(:,:), ALLOCATABLE :: ZP_SFTS     ! scalar flux
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_SFCO2    ! CO2 flux
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_SFU      ! zonal momentum flux
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_SFV      ! meridian momentum flux
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_TSRAD    ! radiative surface temperature
REAL, DIMENSION(:,:), ALLOCATABLE :: ZP_DIR_ALB  ! direct albedo
REAL, DIMENSION(:,:), ALLOCATABLE :: ZP_SCA_ALB  ! diffuse albedo
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_EMIS     ! emissivity

REAL, DIMENSION(:),   ALLOCATABLE :: ZP_TSURF
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_Z0
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_Z0H
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_QSURF

REAL, DIMENSION(:),   ALLOCATABLE :: ZP_PEW_A_COEF ! coefficients for
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_PEW_B_COEF ! implicit coupling
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_PET_A_COEF
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_PEQ_A_COEF
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_PET_B_COEF
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_PEQ_B_COEF
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_RN        ! net radiation           (W/m2)
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_H         ! sensible heat flux      (W/m2)
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_LE        ! Total latent heat flux  (W/m2)
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_LEI       ! Solid Latent heat flux  (W/m2)
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_GFLUX     ! ground flux             (W/m2)
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_T2M       ! Air temperature at 2 meters (K)
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_Q2M       ! Air humidity at 2 meters    (kg/kg)
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_HU2M      ! Air relative humidity at 2 meters (-)
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_ZON10M    ! zonal Wind at 10 meters     (m/s)
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_MER10M    ! meridian Wind at 10 meters  (m/s)
TYPE(LIST_ll), POINTER            :: TZFIELDSURF_ll    ! list of fields to exchange
INTEGER                           :: IINFO_ll       ! return code of parallel routine
!
!
CHARACTER(LEN=6), DIMENSION(:), ALLOCATABLE :: YSV_SURF ! name of the scalar variables
                                                        ! sent to SURFEX
!                                                        
REAL                              :: ZTIMEC
INTEGER           :: ILUOUT         ! logical unit
!
! Fire model
REAL, DIMENSION(2)                    :: ZFIRETIME1, ZFIRETIME2           ! CPU time for Blaze perf profiling
REAL, DIMENSION(2)                    :: ZGRADTIME1, ZGRADTIME2           ! CPU time for Blaze perf profiling
REAL, DIMENSION(2)                    :: ZPROPAGTIME1, ZPROPAGTIME2       ! CPU time for Blaze perf profiling
REAL, DIMENSION(2)                    :: ZFLUXTIME1, ZFLUXTIME2           ! CPU time for Blaze perf profiling
REAL, DIMENSION(2)                    :: ZROSWINDTIME1, ZROSWINDTIME2     ! CPU time for Blaze perf profiling
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZFIREFUELMAP                     ! Fuel map
CHARACTER(LEN=7)                      :: YFUELMAPFILE                     ! Fuel Map file name
TYPE(LIST_ll), POINTER                :: TZFIELDFIRE_ll                   ! list of fields to exchange

!-------------------------------------------------------------------------------
!
!
ILUOUT=TLUOUT%NLU
IKB= 1+JPVEXT
IKU=NKMAX + 2* JPVEXT
IKE=IKU-JPVEXT
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
PSFTH    = XUNDEF
PSFRV    = XUNDEF
PSFSV    = XUNDEF
PSFCO2   = XUNDEF
PSFU     = XUNDEF
PSFV     = XUNDEF
PDIR_ALB = XUNDEF
PSCA_ALB = XUNDEF
PEMIS    = XUNDEF
PTSRAD   = XUNDEF
!
!
!-------------------------------------------------------------------------------
!
!*       1.     CONVERSION OF THE ATMOSPHERIC VARIABLES
!               --------------------------------------- 
!
!        1.1    water vapor
!               -----------

!
ALLOCATE(ZRV(SIZE(PSFTH,1),SIZE(PSFTH,2),IKU))
!
IF(NRR>0) THEN
  ZRV(:,:,:)=XRT(:,:,:,1)
ELSE
  ZRV(:,:,:)=0.
END IF
!
!        1.2    Horizontal wind direction (rad from N clockwise)
!               -------------------------
!
ZU2D(:,:,:)=MXF(XUT(:,:,IKB:IKB))
ZV2D(:,:,:)=MYF(XVT(:,:,IKB:IKB))
!
!* angle between Y axis and wind (rad., clockwise)
!
ZALFA = 0.
WHERE(ZU2D(:,:,1)/=0. .OR. ZV2D(:,:,1)/=0.)
  ZALFA(:,:)=ATAN2(ZU2D(:,:,1),ZV2D(:,:,1))
END WHERE
WHERE(ZALFA(:,:)<0.) ZALFA(:,:) = ZALFA(:,:) + 2. * XPI
!
!* angle between North and wind (rad., clockwise)
!
IF (.NOT. LCARTESIAN) THEN
  ZDIR =   ( (XRPK*(XLON(:,:)-XLON0)) - XBETA ) * XPI/180.  + ZALFA
ELSE
  ZDIR =                              - XBETA   * XPI/180.  + ZALFA
END IF
!
!
!        1.3    Rotate the wind
!               ---------------
!
CALL ROTATE_WIND(D,XUT,XVT,XWT,           &
     XDIRCOSXW, XDIRCOSYW, XDIRCOSZW,   &
     XCOSSLOPE,XSINSLOPE,               &
     XDXX,XDYY,XDZZ,                    &
     ZUA,ZVA                            )

!
!        1.4    zonal and meridian components of the wind parallel to the slope
!               ---------------------------------------------------------------
!
ZWIND(:,:) = SQRT( ZUA**2 + ZVA**2 )
!
ZU(:,:) = ZWIND(:,:) * SIN(ZDIR)
ZV(:,:) = ZWIND(:,:) * COS(ZDIR)
!
!        1.5   Horizontal interpolation the thermodynamic fields
!              -------------------------------------------------
!
CALL NORMAL_INTERPOL(XTHT,ZRV,XPABST,                    &
     XDIRCOSXW, XDIRCOSYW, XDIRCOSZW,                    &
     XCOSSLOPE,XSINSLOPE,                                &
     XDXX,XDYY,XDZZ,                                     &
     ZTHA,ZRVA,ZEXNA                                     )
!
DEALLOCATE(ZRV)
!
!
!        1.6    Pressure and Exner function
!               ---------------------------
!
!
ZPA(:,:) = XP00 * ZEXNA(:,:) **(XCPD/XRD)
!
ZEXNS(:,:) = 0.5 * ( (XPABST(:,:,IKB-1)/XP00)**(XRD/XCPD)  &
                    +(XPABST(:,:,IKB  )/XP00)**(XRD/XCPD)  &
                   )
ZPS(:,:) = XP00 * ZEXNS(:,:) **(XCPD/XRD)
!
!        1.7    humidity in kg/m3 from the mixing ratio
!               ---------------------------------------
!
!
ZQA(:,:) = ZRVA(:,:) * XRHODREF(:,:,IKB)
!
!
!        1.8    Temperature from the potential temperature
!               ------------------------------------------
!
!
ZTA(:,:) = ZTHA(:,:) * ZEXNA(:,:)
!
!
!        1.9    Air density
!               -----------
!
ZRHOA(:,:) = ZPA(:,:)/(XRD * ZTA(:,:) * ((1. + (XRD/XRV)*ZRVA(:,:))/ &
                                         (1. + ZRVA(:,:))))
!
!
!        1.10   Precipitations
!               --------------
!
ZRAIN=0.
ZSNOW=0.
IF (NRR>2 .AND. SIZE(XINPRR)>0 ) THEN
  IF (( CCLOUD(1:3) == 'ICE'                                          .AND. LSEDIC) .OR. &
      ((CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'KHKO') .AND. LSEDC)  .OR. &
      ( CCLOUD=='LIMA'                                                .AND. MSEDC))   THEN
    ZRAIN = ZRAIN + XINPRR * XRHOLW + XINPRC * XRHOLW
  ELSE
    ZRAIN = ZRAIN + XINPRR * XRHOLW
  END IF
END IF
IF (CDCONV == 'KAFR')  THEN
  ZRAIN = ZRAIN + (XPRCONV - XPRSCONV)  * XRHOLW
  ZSNOW = ZSNOW + XPRSCONV * XRHOLW
END IF
IF( NRR    >= 5    .AND. SIZE(XINPRS)>0 ) ZSNOW = ZSNOW + XINPRS * XRHOLW
IF( NRR    >= 6    .AND. SIZE(XINPRG)>0 ) ZSNOW = ZSNOW + XINPRG * XRHOLW
IF( NRR    >= 7    .AND. SIZE(XINPRH)>0 ) ZSNOW = ZSNOW + XINPRH * XRHOLW
!
!
!        1.11   Solar time
!               ----------
!
IF (.NOT. LCARTESIAN) THEN
  ZTSUN(:,:) = MOD(TDTCUR%xtime -XTSIDER*3600. +XLON(:,:)*240., XDAY)
ELSE
  ZTSUN(:,:) = MOD(TDTCUR%xtime -XTSIDER*3600. +XLON0    *240., XDAY)
END IF
!
!        1.12   Forcing level
!               -------------
!
ZZREF(:,:) = 0.5*( XZZ(:,:,IKB+1)-XZZ(:,:,IKB) )*XDIRCOSZW(:,:)
!
!
!        1.13   CO2 concentration (kg/m3)
!               -----------------
!
ZCO2(:,:) = XCCO2 * XRHODREF(:,:,IKB)
!
!
!
!        1.14   Blowing snow scheme (optional)
!               -----------------
!
ZBLOWSNOW_2D=0.

IF(LBLOWSNOW) THEN
    KSV_SURF = NSV+NBLOWSNOW_2D ! When blowing snow scheme is used
                  ! NBLOWSN0W_2D variables are sent to SURFEX through ZP_SV.
                  ! They refer to the 2D fields advected by MNH including:
                  !             - total number concentration in Canopy
                  !             - total mass concentration in Canopy
                  !             - equivalent concentration in the saltation layer
    ! Initialize array of scalar to be sent to SURFEX including 2D blowing snow fields
    ALLOCATE(YSV_SURF(KSV_SURF))
    YSV_SURF(1:NSV)          = CSV(:)
    YSV_SURF(NSV+1:KSV_SURF) = YPBLOWSNOW_2D(:)


    DO JSV=1,NBLOWSNOW_2D 
       ZBLOWSNOW_2D(:,:,JSV) = XRSNWCANOS(:,:,JSV)*XTSTEP/XRHODJ(:,:,IKB)
    END DO

ELSE
    KSV_SURF = NSV
    ALLOCATE(YSV_SURF(KSV_SURF))
    YSV_SURF(:)     = CSV(:)
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.     Call to surface monitor with 2D variables
!               -----------------------------------------
!
!
! initial values:
!
IDIM1       = IIE-IIB+1
IDIM2       = IJE-IJB+1
IDIM1D      = IDIM1*IDIM2
!
!
! Transform 2D input fields into 1D:
!
CALL RESHAPE_SURF(IDIM1D)
!
! call to have the cumulated time since beginning of simulation
!
CALL DATETIME_DISTANCE(TDTSEG,TDTCUR,ZTIMEC)

#ifdef CPLOASIS
IF (LOASIS) THEN
  IF ( MOD(ZTIMEC,1.0) .LE. 1E-2 .OR. (1.0 - MOD(ZTIMEC,1.0)) .LE. 1E-2 ) THEN
    IF ( NINT(ZTIMEC-(XSEGLEN-DYN_MODEL(1)%XTSTEP)) .LT. 0 ) THEN
      WRITE(ILUOUT,*) '----------------------------'
      WRITE(ILUOUT,*) ' Reception des champs avec OASIS'
      WRITE(ILUOUT,*) 'NINT(ZTIMEC)=', NINT(ZTIMEC)
      CALL MNH_OASIS_RECV(CPROGRAM,IDIM1D,SIZE(XSW_BANDS),ZTIMEC+XTSTEP,XTSTEP,         &
                        ZP_ZENITH,XSW_BANDS                                         ,         &
                        ZP_TSRAD,ZP_DIR_ALB,ZP_SCA_ALB,ZP_EMIS,ZP_TSURF)
      WRITE(ILUOUT,*) '----------------------------'
    END IF
  END IF
END IF
#endif
!
! Call to surface schemes
!                       
CALL COUPLING_SURF_ATM_n(YSURF_CUR,'MESONH', 'E',ZTIMEC,                                            &
               XTSTEP, TDTCUR%nyear, TDTCUR%nmonth, TDTCUR%nday, TDTCUR%xtime,                      &
               IDIM1D,KSV_SURF,SIZE(XSW_BANDS),                                                     &
               ZP_TSUN, ZP_ZENITH,ZP_ZENITH, ZP_AZIM,                                               &
               ZP_ZREF, ZP_ZREF, ZP_ZS, ZP_U, ZP_V, ZP_QA, ZP_TA, ZP_RHOA, ZP_SV, ZP_CO2, YSV_SURF, &
               ZP_RAIN, ZP_SNOW, ZP_LW, ZP_DIR_SW, ZP_SCA_SW, XSW_BANDS, ZP_PS, ZP_PA,              &
               ZP_SFTQ, ZP_SFTH, ZP_SFTS, ZP_SFCO2, ZP_SFU, ZP_SFV,                                 &
               ZP_TSRAD, ZP_DIR_ALB, ZP_SCA_ALB, ZP_EMIS, ZP_TSURF, ZP_Z0, ZP_Z0H, ZP_QSURF,        &
               ZP_PEW_A_COEF, ZP_PEW_B_COEF,                                                        &
               ZP_PET_A_COEF, ZP_PEQ_A_COEF, ZP_PET_B_COEF, ZP_PEQ_B_COEF,ZP_ZWS,                   &
               'OK'                                                                                 )
!
#ifdef CPLOASIS
IF (LOASIS) THEN
  IF ( MOD(ZTIMEC,1.0) .LE. 1E-2 .OR. (1.0 - MOD(ZTIMEC,1.0)) .LE. 1E-2 ) THEN
    IF (NINT(ZTIMEC-(XSEGLEN-DYN_MODEL(1)%XTSTEP)) .LT. 0) THEN
      WRITE(ILUOUT,*) '----------------------------'
      WRITE(ILUOUT,*) ' Envoi des champs avec OASIS'
      WRITE(ILUOUT,*) 'NINT(ZTIMEC)=', NINT(ZTIMEC)
      CALL MNH_OASIS_SEND(CPROGRAM,IDIM1D,ZTIMEC+XTSTEP,XTSTEP)
      WRITE(ILUOUT,*) '----------------------------'
    END IF
  END IF
END IF
#endif
!
IF (CPROGRAM=='DIAG  ' .OR. LDIAG_IN_RUN) THEN
  CALL DIAG_SURF_ATM_n(YSURF_CUR,'MESONH')
  CALL  MNHGET_SURF_PARAM_n( PRN = ZP_RN,         PH = ZP_H,           PLE = ZP_LE,   PLEI = ZP_LEI,   &
                             PGFLUX = ZP_GFLUX,   PT2M = ZP_T2M,       PQ2M = ZP_Q2M, PHU2M = ZP_HU2M, &
                             PZON10M = ZP_ZON10M, PMER10M = ZP_MER10M                                  )
END IF
!
! Transform 1D output fields into 2D:
!
CALL UNSHAPE_SURF(IDIM1,IDIM2)
#ifdef MNH_FOREFIRE
!------------------------!
! COUPLING WITH FOREFIRE !
!------------------------!
	
IF ( LFOREFIRE ) THEN
	CALL FOREFIRE_DUMP_FIELDS_n(XUT, XVT, XWT, XSVT&
	           , XTHT, XRT(:,:,:,1), XPABST, XTKET&
	           , IDIM1+2, IDIM2+2, NKMAX+2)
END IF

IF ( FFCOUPLING ) THEN

	CALL SEND_GROUND_WIND_n(XUT, XVT, IKB, IINFO_ll)
	
	CALL FOREFIRE_RECEIVE_PARAL_n()
   
   CALL COUPLING_FOREFIRE_n(XTSTEP, ZSFTH, ZSFTQ, ZSFTS)
	
	CALL FOREFIRE_SEND_PARAL_n(IINFO_ll)
   
END IF

FF_TIME = FF_TIME + XTSTEP
#endif
!
! Friction of components along slope axes (U: largest local slope axis, V: zero slope axis)
!
! 
PSFU(:,:) = 0.
PSFV(:,:) = 0.
!
WHERE (ZSFU(:,:)/=XUNDEF .AND. ZWIND(:,:)>0.)
  PSFU(:,:) = - SQRT(ZSFU**2+ZSFV**2) * ZUA(:,:) / ZWIND(:,:) / XRHODREF(:,:,IKB)
  PSFV(:,:) = - SQRT(ZSFU**2+ZSFV**2) * ZVA(:,:) / ZWIND(:,:) / XRHODREF(:,:,IKB)
END WHERE
!

!*       2.1    Blaze Fire Model
!               ----------------
!
IF (LBLAZE) THEN
  ! get start time
  CALL SECOND_MNH2( ZFIRETIME1 )

  !*       2.1.1  Local variables allocation
  !               --------------------------
  !

  ! Parallel fuel
  NULLIFY(TZFIELDFIRE_ll)
  IF (KTCOUNT <= 1) THEN
    ! fuelmap
    SELECT CASE (CPROPAG_MODEL)
    CASE('SANTONI2011')
      !
      ALLOCATE( ZFIREFUELMAP(SIZE(XLSPHI,1), SIZE(XLSPHI,2), SIZE(XLSPHI,3), 22) );
      ! Parallel fuel
      CALL ADD4DFIELD_ll( TZFIELDFIRE_ll, ZFIREFUELMAP(:,:,:,1::22), 'MODEL_n::ZFIREFUELMAP' )
      ! Default value
      ZFIREFUELMAP(:,:,:,:) = 0.
    END SELECT

    !*       2.1.2  Read fuel map file
    !               ------------------
    !
    ! Fuel map file name
    YFUELMAPFILE = 'FuelMap'
    !
    CALL FIRE_READFUEL( TPFILE, ZFIREFUELMAP, XFMIGNITION, XFMWALKIG )

    !*       2.1.3  Ignition LS function with ignition map
    !               --------------------------------------
    !
    SELECT CASE (CFIRE_CPL_MODE)
    CASE('2WAYCPL', 'ATM2FIR')
      ! force ignition
      WHERE (XFMIGNITION <= TDTCUR%XTIME ) XLSPHI = 1.
      ! walking ignition
      CALL FIRE_LS_RECONSTRUCTION_FROM_BMAP( XLSPHI, XFMWALKIG, 0.)
      !
      !*       2.1.4  Update BMAP
      !               -----------
      !
      WHERE (XLSPHI >= .5 .AND. XBMAP < 0) XBMAP = TDTCUR%XTIME
      !
    CASE('FIR2ATM')
      CALL FIRE_READBMAP(TPFILE,XBMAP)
      
    END SELECT
    !
    !*       2.1.5  Compute R0, A, Wf0, R00
    !               -----------------------
    !
    SELECT CASE (CPROPAG_MODEL)
    CASE('SANTONI2011')
      CALL FIRE_NOWINDROS( ZFIREFUELMAP, XFMR0, XFMRFA, XFMWF0, XFMR00, XFMFUELTYPE, XFIRETAU, XFLUXPARAMH, &
                           XFLUXPARAMW, XFMASE, XFMAWC )
    END SELECT
    !
    !*       2.1.6  Compute orographic gradient
    !               ---------------------------
    CALL FIRE_GRAD_OROGRAPHY( XZS, XFMGRADOROX, XFMGRADOROY )
    !
    !*       2.1.7  Test halo size
    !               --------------
    IF (NHALO < 2 .AND. NFIRE_WENO_ORDER == 3) THEN
      WRITE(ILUOUT,'(A/A)') 'ERROR BLAZE-FIRE : WENO3 fire gradient calculation needs NHALO >= 2'
      !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','GROUND_PARAM_n','')
    ELSEIF (NHALO < 3 .AND. NFIRE_WENO_ORDER == 5) THEN
      WRITE(ILUOUT,'(A/A)') 'ERROR : WENO5 fire gradient calculation needs NHALO >= 3'
      !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','GROUND_PARAM_n','')
    END IF
    !
  END IF
  !
  !*       2.1.6  Compute grad of level set function phi
  !               --------------------------------------
  !
  SELECT CASE (CFIRE_CPL_MODE)
  CASE('2WAYCPL', 'ATM2FIR')
    ! get time 1
    CALL SECOND_MNH2( ZGRADTIME1 )
    CALL FIRE_GRADPHI( XLSPHI, XGRADLSPHIX, XGRADLSPHIY )

    ! get time 2
    CALL SECOND_MNH2( ZGRADTIME2 )
    XGRADPERF = XGRADPERF + ZGRADTIME2 - ZGRADTIME1
    !
    !*       2.1.7  Get horizontal wind speed projected on LS gradient direction
    !               ------------------------------------------------------------
    !
    CALL FIRE_GETWIND( XUT, XVT, XWT, XGRADLSPHIX, XGRADLSPHIY, XFIREWIND, KTCOUNT, XTSTEP, XFMGRADOROX, XFMGRADOROY )
    !
    !*       2.1.8  Compute ROS XFIRERW with wind
    !               -----------------------------
    !
    !
    SELECT CASE (CPROPAG_MODEL)
    CASE('SANTONI2011')
      CALL FIRE_RATEOFSPREAD( XFMFUELTYPE, XFMR0, XFMRFA, XFMWF0, XFMR00, XFIREWIND, XGRADLSPHIX, XGRADLSPHIY, &
                              XFMGRADOROX, XFMGRADOROY, XFIRERW )
    END SELECT
    CALL SECOND_MNH2( ZROSWINDTIME2 )
    XROSWINDPERF = XROSWINDPERF + ZROSWINDTIME2 - ZGRADTIME2
    !
    !*       2.1.8  Integrate model on atm time step to propagate
    !               ---------------------------------------------
    !
    SELECT CASE (CPROPAG_MODEL)
    CASE('SANTONI2011')
      CALL FIRE_PROPAGATE( XLSPHI, XBMAP, XFMIGNITION, XFMWALKIG, XGRADLSPHIX, XGRADLSPHIY, XTSTEP, XFIRERW )
    END SELECT
    CALL SECOND_MNH2( ZPROPAGTIME2 )
    XPROPAGPERF = XPROPAGPERF + ZPROPAGTIME2 - ZROSWINDTIME2
    !
  CASE('FIR2ATM')
    !
    CALL SECOND_MNH2( ZPROPAGTIME1 )
    CALL FIRE_LS_RECONSTRUCTION_FROM_BMAP( XLSPHI, XBMAP, XTSTEP )
    CALL SECOND_MNH2( ZPROPAGTIME2 )
    XPROPAGPERF = XPROPAGPERF + ZPROPAGTIME2 - ZPROPAGTIME1
    XGRADPERF(:) = 0.
    !
  END SELECT
  !
  !*       2.1.8  Compute fluxes
  !               --------------
  !
  SELECT CASE (CFIRE_CPL_MODE)
  CASE('2WAYCPL','FIR2ATM')
    CALL SECOND_MNH2( ZFLUXTIME1 )
    ! 2 way coupling
    CALL FIRE_HEATFLUXES( XLSPHI, XBMAP, XFIRETAU, XTSTEP, XFLUXPARAMH, XFLUXPARAMW, XFMFLUXHDH, XFMFLUXHDW, XFMASE, XFMAWC )
    ! vertical distribution of fire heat fluxes
    CALL FIRE_VERTICALFLUXDISTRIB( XFMFLUXHDH, XFMFLUXHDW, XRTHS, XRRS, ZSFTS, XEXNREF, XRHODJ, XRT, XRHODREF )
    !
    CALL SECOND_MNH2( ZFLUXTIME2 )
    XFLUXPERF = XFLUXPERF + ZFLUXTIME2 - ZFLUXTIME1
  CASE DEFAULT
    XFLUXPERF(:) = 0.
  END SELECT
  ! get end time
  CALL SECOND_MNH2( ZFIRETIME2 )
  ! add to Blaze time
  XFIREPERF = XFIREPERF + ZFIRETIME2 - ZFIRETIME1
END IF
!* conversion from H (W/m2) to w'Theta'
!
PSFTH(:,:) = ZSFTH(:,:) /  XCPD / XRHODREF(:,:,IKB)
!
!
!* conversion from water flux (kg/m2/s) to w'rv'
!
PSFRV(:,:) = ZSFTQ(:,:) / XRHODREF(:,:,IKB)
!
!
!* conversion from scalar flux (kg/m2/s) to w'rsv'
!
IF(NSV .GT. 0) THEN
   DO JSV=1,NSV
     PSFSV(:,:,JSV) = ZSFTS(:,:,JSV) / XRHODREF(:,:,IKB)
   END DO
END IF
!
!* conversion from chemistry flux (molec/m2/s) to (ppp.m.s-1)
!
IF (LUSECHEM) THEN
   DO JSV=NSV_CHEMBEG,NSV_CHEMEND
      PSFSV(:,:,JSV) = ZSFTS(:,:,JSV) * XMD / ( XAVOGADRO * XRHODREF(:,:,IKB))
      IF ((LCHEMDIAG).AND.(CPROGRAM == 'DIAG  ')) XCHFLX(:,:,JSV-NSV_CHEMBEG+1) = PSFSV(:,:,JSV)    
   END DO
ELSE
  PSFSV(:,:,NSV_CHEMBEG:NSV_CHEMEND) = 0.
END IF
!
!* conversion from dust flux (kg/m2/s) to (ppp.m.s-1)
!
IF (LDUST) THEN
  DO JSV=NSV_DSTBEG,NSV_DSTEND
     PSFSV(:,:,JSV) = ZSFTS(:,:,JSV) * XMD / (XMOLARWEIGHT_DUST * XRHODREF(:,:,IKB))
  END DO
ELSE
  PSFSV(:,:,NSV_DSTBEG:NSV_DSTEND) = 0.
END IF
!
!* conversion from sea salt flux (kg/m2/s) to (ppp.m.s-1)
!
IF (LSALT) THEN
  DO JSV=NSV_SLTBEG,NSV_SLTEND
     PSFSV(:,:,JSV) = ZSFTS(:,:,JSV) * XMD / (XMOLARWEIGHT_SALT * XRHODREF(:,:,IKB))
  END DO
ELSE
  PSFSV(:,:,NSV_SLTBEG:NSV_SLTEND) = 0.
END IF
!
!* conversion from aerosol flux (molec/m2/s) to (ppp.m.s-1)
!
IF (LORILAM) THEN
  DO JSV=NSV_AERBEG,NSV_AEREND
    PSFSV(:,:,JSV) = ZSFTS(:,:,JSV) * XMD / ( XAVOGADRO * XRHODREF(:,:,IKB))
  END DO
ELSE
  PSFSV(:,:,NSV_AERBEG:NSV_AEREND) = 0.
END IF
!
!* conversion from blowing snow flux (kg/m2/s) to [kg(snow)/kg(dry air).m.s-1]
!
IF (LBLOWSNOW) THEN
  DO JSV=NSV_SNWBEG,NSV_SNWEND
     PSFSV(:,:,JSV) = ZSFTS(:,:,JSV)/ (ZRHOA(:,:))
  END DO
  !* Update tendency for blowing snow 2D fields
  DO JSV=1,(NBLOWSNOW_2D)
     XRSNWCANOS(:,:,JSV) = ZBLOWSNOW_2D(:,:,JSV)*XRHODJ(:,:,IKB)/(XTSTEP*ZRHOA(:,:))
  END DO

ELSE
  PSFSV(:,:,NSV_SNWBEG:NSV_SNWEND) = 0.
END IF
!
!* conversion from CO2 flux (kg/m2/s) to w'CO2'
!
PSFCO2(:,:) = ZSFCO2(:,:) / XRHODREF(:,:,IKB)
!
!
!*  Diagnostics
!   -----------
!
!
IF (LDIAG_IN_RUN) THEN
  !
  XCURRENT_SFCO2(:,:) = ZSFCO2(:,:)
  XCURRENT_DSTAOD(:,:)=0.0
  XCURRENT_SLTAOD(:,:)=0.0
  IF (CRAD/='NONE') THEN
  XCURRENT_LWD  (:,:) = XFLALWD(:,:)
  XCURRENT_SWD  (:,:) = SUM(XDIRSRFSWD(:,:,:)+XSCAFLASWD(:,:,:),DIM=3)
  XCURRENT_LWU  (:,:) = XLWU(:,:,IKB) 
  XCURRENT_SWU  (:,:) = XSWU(:,:,IKB) 
  XCURRENT_SWDIR(:,:) = SUM(XDIRSRFSWD,DIM=3)
  XCURRENT_SWDIFF(:,:) = SUM(XSCAFLASWD(:,:,:),DIM=3) 
  DO JK=IKB,IKE
    IKRAD = JK - 1 
    DO JJ=IJB,IJE
       DO JI=IIB,IIE
         XCURRENT_DSTAOD(JI,JJ)=XCURRENT_DSTAOD(JI,JJ)+XAER(JI,JJ,IKRAD,3)
         XCURRENT_SLTAOD(JI,JJ)=XCURRENT_SLTAOD(JI,JJ)+XAER(JI,JJ,IKRAD,2)
       ENDDO
    ENDDO
  ENDDO
  END IF
!
  NULLIFY(TZFIELDSURF_ll)
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_RN,     'GROUND_PARAM_n::XCURRENT_RN'     )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_H,      'GROUND_PARAM_n::XCURRENT_H'      )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_LE,     'GROUND_PARAM_n::XCURRENT_LE'     )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_GFLUX,  'GROUND_PARAM_n::XCURRENT_GFLUX'  )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_SWD,    'GROUND_PARAM_n::XCURRENT_SWD'    )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_SWU,    'GROUND_PARAM_n::XCURRENT_SWU'    )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_LWD,    'GROUND_PARAM_n::XCURRENT_LWD'    )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_LWU,    'GROUND_PARAM_n::XCURRENT_LWU'    )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_T2M,    'GROUND_PARAM_n::XCURRENT_T2M'    )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_Q2M,    'GROUND_PARAM_n::XCURRENT_Q2M'    )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_HU2M,   'GROUND_PARAM_n::XCURRENT_HU2M'   )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_ZON10M, 'GROUND_PARAM_n::XCURRENT_ZON10M' )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_MER10M, 'GROUND_PARAM_n::XCURRENT_MER10M' )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_DSTAOD, 'GROUND_PARAM_n::XCURRENT_DSTAOD' )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_SLTAOD, 'GROUND_PARAM_n::XCURRENT_SLTAOD' )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_ZWS,    'GROUND_PARAM_n::XCURRENT_ZWS'    )
  CALL ADD2DFIELD_ll( TZFIELDSURF_ll,XCURRENT_SFCO2,  'GROUND_PARAM_n::XCURRENT_SFCO2'  )

  CALL UPDATE_HALO_ll(TZFIELDSURF_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDSURF_ll)
END IF
!
IF (LBLAZE) THEN
  IF (KTCOUNT <= 1) THEN
      DEALLOCATE(ZFIREFUELMAP)
  END IF
  CALL CLEANLIST_ll(TZFIELDFIRE_ll)
END IF
!==================================================================================
!
CONTAINS
!
!==================================================================================
!
SUBROUTINE RESHAPE_SURF(KDIM1D)
!
INTEGER, INTENT(IN)   :: KDIM1D
INTEGER, DIMENSION(1) :: ISHAPE_1
!
ISHAPE_1 = (/KDIM1D/)
!
ALLOCATE(ZP_TSUN    (KDIM1D))
ALLOCATE(ZP_ZENITH  (KDIM1D))
ALLOCATE(ZP_AZIM    (KDIM1D))
ALLOCATE(ZP_ZREF    (KDIM1D))
ALLOCATE(ZP_ZS      (KDIM1D))
ALLOCATE(ZP_U       (KDIM1D))
ALLOCATE(ZP_V       (KDIM1D))
ALLOCATE(ZP_QA      (KDIM1D))
ALLOCATE(ZP_TA      (KDIM1D))
ALLOCATE(ZP_RHOA    (KDIM1D))
ALLOCATE(ZP_SV      (KDIM1D,KSV_SURF))
ALLOCATE(ZP_CO2     (KDIM1D))
ALLOCATE(ZP_RAIN    (KDIM1D))
ALLOCATE(ZP_SNOW    (KDIM1D))
ALLOCATE(ZP_LW      (KDIM1D))
ALLOCATE(ZP_DIR_SW  (KDIM1D,SIZE(XDIRSRFSWD,3)))
ALLOCATE(ZP_SCA_SW  (KDIM1D,SIZE(XSCAFLASWD,3)))
ALLOCATE(ZP_PS      (KDIM1D))
ALLOCATE(ZP_PA      (KDIM1D))
ALLOCATE(ZP_ZWS     (KDIM1D))

ALLOCATE(ZP_SFTQ    (KDIM1D))
ALLOCATE(ZP_SFTH    (KDIM1D))
ALLOCATE(ZP_SFU     (KDIM1D))
ALLOCATE(ZP_SFV     (KDIM1D))
ALLOCATE(ZP_SFTS    (KDIM1D,KSV_SURF))
ALLOCATE(ZP_SFCO2   (KDIM1D))
ALLOCATE(ZP_TSRAD   (KDIM1D))
ALLOCATE(ZP_DIR_ALB (KDIM1D,SIZE(PDIR_ALB,3)))
ALLOCATE(ZP_SCA_ALB (KDIM1D,SIZE(PSCA_ALB,3)))
ALLOCATE(ZP_EMIS    (KDIM1D))
ALLOCATE(ZP_TSURF   (KDIM1D))
ALLOCATE(ZP_Z0      (KDIM1D))
ALLOCATE(ZP_Z0H     (KDIM1D))
ALLOCATE(ZP_QSURF   (KDIM1D))
ALLOCATE(ZP_RN      (KDIM1D))
ALLOCATE(ZP_H       (KDIM1D))
ALLOCATE(ZP_LE      (KDIM1D))
ALLOCATE(ZP_LEI     (KDIM1D))
ALLOCATE(ZP_GFLUX   (KDIM1D))
ALLOCATE(ZP_T2M     (KDIM1D))
ALLOCATE(ZP_Q2M     (KDIM1D))
ALLOCATE(ZP_HU2M    (KDIM1D))
ALLOCATE(ZP_ZON10M  (KDIM1D))
ALLOCATE(ZP_MER10M  (KDIM1D))

!* explicit coupling only
ALLOCATE(ZP_PEW_A_COEF  (KDIM1D))
ALLOCATE(ZP_PEW_B_COEF  (KDIM1D))
ALLOCATE(ZP_PET_A_COEF  (KDIM1D))
ALLOCATE(ZP_PEQ_A_COEF  (KDIM1D))
ALLOCATE(ZP_PET_B_COEF  (KDIM1D))
ALLOCATE(ZP_PEQ_B_COEF  (KDIM1D))

ZP_TSUN(:)        = RESHAPE(ZTSUN(IIB:IIE,IJB:IJE),        ISHAPE_1)
ZP_TA(:)          = RESHAPE(ZTA(IIB:IIE,IJB:IJE),          ISHAPE_1)
ZP_QA(:)          = RESHAPE(ZQA(IIB:IIE,IJB:IJE),          ISHAPE_1)
ZP_RHOA(:)        = RESHAPE(ZRHOA(IIB:IIE,IJB:IJE),        ISHAPE_1)
ZP_U(:)           = RESHAPE(ZU(IIB:IIE,IJB:IJE),           ISHAPE_1)
ZP_V(:)           = RESHAPE(ZV(IIB:IIE,IJB:IJE),           ISHAPE_1)
ZP_PS(:)          = RESHAPE(ZPS(IIB:IIE,IJB:IJE),          ISHAPE_1)
ZP_PA(:)          = RESHAPE(ZPA(IIB:IIE,IJB:IJE),          ISHAPE_1)
ZP_ZS(:)          = RESHAPE(XZS(IIB:IIE,IJB:IJE),          ISHAPE_1)
ZP_CO2(:)         = RESHAPE(ZCO2(IIB:IIE,IJB:IJE),         ISHAPE_1)
ZP_SNOW(:)        = RESHAPE(ZSNOW(IIB:IIE,IJB:IJE),        ISHAPE_1)
ZP_RAIN(:)        = RESHAPE(ZRAIN(IIB:IIE,IJB:IJE),        ISHAPE_1)
ZP_ZREF(:)        = RESHAPE(ZZREF(IIB:IIE,IJB:IJE),        ISHAPE_1)
ZP_ZWS(:)         = RESHAPE(XZWS(IIB:IIE,IJB:IJE),         ISHAPE_1)

DO JLAYER=1,NSV
  ZP_SV(:,JLAYER) = RESHAPE(XSVT(IIB:IIE,IJB:IJE,IKB,JLAYER), ISHAPE_1)
END DO
!
IF(LBLOWSNOW) THEN
  DO JLAYER=1,NBLOWSNOW_2D
      ZP_SV(:,NSV+JLAYER) = RESHAPE(ZBLOWSNOW_2D(IIB:IIE,IJB:IJE,JLAYER), ISHAPE_1)
  END DO
END IF
!
!chemical conversion : from part/part to molec./m3
DO JLAYER=NSV_CHEMBEG,NSV_CHEMEND
  ZP_SV(:,JLAYER) = ZP_SV(:,JLAYER) * XAVOGADRO * ZP_RHOA(:) / XMD
END DO
DO JLAYER=NSV_AERBEG,NSV_AEREND
  ZP_SV(:,JLAYER) = ZP_SV(:,JLAYER) * XAVOGADRO * ZP_RHOA(:) / XMD
END DO
!dust  conversion : from part/part to kg/m3
DO JLAYER=NSV_DSTBEG,NSV_DSTEND
 ZP_SV(:,JLAYER) = ZP_SV(:,JLAYER) *  XMOLARWEIGHT_DUST* ZP_RHOA(:) / XMD
END DO
!sea salt  conversion : from part/part to kg/m3
DO JLAYER=NSV_SLTBEG,NSV_SLTEND
 ZP_SV(:,JLAYER) = ZP_SV(:,JLAYER) *  XMOLARWEIGHT_SALT* ZP_RHOA(:) / XMD
END DO
!
!blowing snow conversion : from kg(snow)/kg(dry air) to kg(snow)/m3
DO JLAYER=NSV_SNWBEG,NSV_SNWEND
 ZP_SV(:,JLAYER) = ZP_SV(:,JLAYER) * ZP_RHOA(:)
END DO

IF(LBLOWSNOW) THEN ! Convert 2D blowing snow fields
                    ! from kg(snow)/kg(dry air) to kg(snow)/m3
  DO JLAYER=(NSV+1),KSV_SURF
     ZP_SV(:,JLAYER) = ZP_SV(:,JLAYER) * ZP_RHOA(:)
  END DO
END IF
!
ZP_ZENITH(:)      = RESHAPE(XZENITH(IIB:IIE,IJB:IJE),      ISHAPE_1)
ZP_AZIM  (:)      = RESHAPE(XAZIM  (IIB:IIE,IJB:IJE),      ISHAPE_1)
ZP_LW(:)          = RESHAPE(XFLALWD(IIB:IIE,IJB:IJE),      ISHAPE_1)
DO JLAYER=1,SIZE(XDIRSRFSWD,3)
  ZP_DIR_SW(:,JLAYER)  = RESHAPE(XDIRSRFSWD(IIB:IIE,IJB:IJE,JLAYER),  ISHAPE_1)
  ZP_SCA_SW(:,JLAYER)  = RESHAPE(XSCAFLASWD(IIB:IIE,IJB:IJE,JLAYER),  ISHAPE_1)
END DO
!
ZP_PEW_A_COEF = 0.
ZP_PEW_B_COEF = 0.
ZP_PET_A_COEF = 0.
ZP_PEQ_A_COEF = 0.
ZP_PET_B_COEF = 0.
ZP_PEQ_B_COEF = 0.
!
END SUBROUTINE RESHAPE_SURF
!================================================i=================================
SUBROUTINE UNSHAPE_SURF(KDIM1,KDIM2)
!
INTEGER, INTENT(IN)   :: KDIM1, KDIM2
INTEGER, DIMENSION(2) :: ISHAPE_2
!
ISHAPE_2 = (/KDIM1,KDIM2/)
!
! Arguments in call to surface:
!
ZSFTH = XUNDEF
ZSFTQ = XUNDEF
IF (NSV>0) ZSFTS = XUNDEF
ZSFCO2 = XUNDEF
ZSFU = XUNDEF
ZSFV = XUNDEF
!
ZSFTH   (IIB:IIE,IJB:IJE)       = RESHAPE(ZP_SFTH(:),       ISHAPE_2)
ZSFTQ   (IIB:IIE,IJB:IJE)       = RESHAPE(ZP_SFTQ(:),       ISHAPE_2)
DO JLAYER=1,SIZE(PSFSV,3)
  ZSFTS   (IIB:IIE,IJB:IJE,JLAYER) = RESHAPE(ZP_SFTS(:,JLAYER),  ISHAPE_2)
END DO
ZSFCO2  (IIB:IIE,IJB:IJE)       = RESHAPE(ZP_SFCO2(:),    ISHAPE_2)
ZSFU    (IIB:IIE,IJB:IJE)       = RESHAPE(ZP_SFU(:),      ISHAPE_2)
ZSFV    (IIB:IIE,IJB:IJE)       = RESHAPE(ZP_SFV(:),      ISHAPE_2)
DO JLAYER=1,SIZE(PEMIS,3)
    PEMIS   (IIB:IIE,IJB:IJE,JLAYER)     = RESHAPE(ZP_EMIS(:),     ISHAPE_2)
END DO
PTSRAD  (IIB:IIE,IJB:IJE)       = RESHAPE(ZP_TSRAD(:),    ISHAPE_2)
IF(LBLOWSNOW) THEN
  DO JLAYER=1,NBLOWSNOW_2D
    ZBLOWSNOW_2D(IIB:IIE,IJB:IJE,JLAYER)  =  RESHAPE(ZP_SFTS(:,NSV+JLAYER),  ISHAPE_2)
  END DO
END IF
!
IF (LDIAG_IN_RUN) THEN
  XCURRENT_RN      (IIB:IIE,IJB:IJE)  = RESHAPE(ZP_RN(:),     ISHAPE_2)
  XCURRENT_H       (IIB:IIE,IJB:IJE)  = RESHAPE(ZP_H (:),     ISHAPE_2)
  XCURRENT_LE      (IIB:IIE,IJB:IJE)  = RESHAPE(ZP_LE(:),     ISHAPE_2)
  XCURRENT_LEI     (IIB:IIE,IJB:IJE)  = RESHAPE(ZP_LEI(:),    ISHAPE_2)
  XCURRENT_GFLUX   (IIB:IIE,IJB:IJE)  = RESHAPE(ZP_GFLUX(:),  ISHAPE_2)
  XCURRENT_T2M     (IIB:IIE,IJB:IJE)  = RESHAPE(ZP_T2M(:),    ISHAPE_2)
  XCURRENT_Q2M     (IIB:IIE,IJB:IJE)  = RESHAPE(ZP_Q2M(:),    ISHAPE_2)
  XCURRENT_HU2M    (IIB:IIE,IJB:IJE)  = RESHAPE(ZP_HU2M(:),   ISHAPE_2)
  XCURRENT_ZON10M  (IIB:IIE,IJB:IJE)  = RESHAPE(ZP_ZON10M(:), ISHAPE_2)
  XCURRENT_MER10M  (IIB:IIE,IJB:IJE)  = RESHAPE(ZP_MER10M(:), ISHAPE_2)
  XCURRENT_ZWS  (IIB:IIE,IJB:IJE)  = RESHAPE(ZP_ZWS(:), ISHAPE_2)
ENDIF
!
DO JLAYER=1,SIZE(PDIR_ALB,3)
  PDIR_ALB(IIB:IIE,IJB:IJE,JLAYER) = RESHAPE(ZP_DIR_ALB(:,JLAYER), ISHAPE_2)
  PSCA_ALB(IIB:IIE,IJB:IJE,JLAYER) = RESHAPE(ZP_SCA_ALB(:,JLAYER), ISHAPE_2)
END DO
!
DEALLOCATE(ZP_TSUN    )
DEALLOCATE(ZP_ZENITH  )
DEALLOCATE(ZP_AZIM    )
DEALLOCATE(ZP_ZREF    )
DEALLOCATE(ZP_ZS      )
DEALLOCATE(ZP_U       )
DEALLOCATE(ZP_V       )
DEALLOCATE(ZP_QA      )
DEALLOCATE(ZP_TA      )
DEALLOCATE(ZP_RHOA    )
DEALLOCATE(ZP_SV      )
DEALLOCATE(ZP_CO2     )
DEALLOCATE(ZP_RAIN    )
DEALLOCATE(ZP_SNOW    )
DEALLOCATE(ZP_LW      )
DEALLOCATE(ZP_DIR_SW  )
DEALLOCATE(ZP_SCA_SW  )
DEALLOCATE(ZP_PS      )
DEALLOCATE(ZP_PA      )
DEALLOCATE(ZP_ZWS     )

DEALLOCATE(ZP_SFTQ    )
DEALLOCATE(ZP_SFTH    )
DEALLOCATE(ZP_SFTS    )
DEALLOCATE(ZP_SFCO2   )
DEALLOCATE(ZP_SFU     )
DEALLOCATE(ZP_SFV     )
DEALLOCATE(ZP_TSRAD   )
DEALLOCATE(ZP_DIR_ALB )
DEALLOCATE(ZP_SCA_ALB )
DEALLOCATE(ZP_EMIS    )
DEALLOCATE(ZP_RN      )
DEALLOCATE(ZP_H       )
DEALLOCATE(ZP_LE      )
DEALLOCATE(ZP_LEI     )
DEALLOCATE(ZP_GFLUX   )
DEALLOCATE(ZP_T2M     )
DEALLOCATE(ZP_Q2M     )
DEALLOCATE(ZP_HU2M     )
DEALLOCATE(ZP_ZON10M  )
DEALLOCATE(ZP_MER10M  )

DEALLOCATE(ZP_PEW_A_COEF  )
DEALLOCATE(ZP_PEW_B_COEF  )
DEALLOCATE(ZP_PET_A_COEF  )
DEALLOCATE(ZP_PEQ_A_COEF  )
DEALLOCATE(ZP_PET_B_COEF  )
DEALLOCATE(ZP_PEQ_B_COEF  )
!
END SUBROUTINE UNSHAPE_SURF
!==================================================================================
!
END SUBROUTINE GROUND_PARAM_n
