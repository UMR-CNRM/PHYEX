!MNH_LIC Copyright 1995-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 25/02/2019: split rain_ice (cleaner and easier to maintain/debug)
!  P. Wautelet 03/06/2019: remove PACK/UNPACK intrinsics (to get more performance and better OpenACC support)
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_WARM

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: RAIN_ICE_WARM

CONTAINS

SUBROUTINE RAIN_ICE_WARM(OMICRO, KMICRO, K1, K2, K3,                                                           &
                         PRHODREF, PRVT, PRCT, PRRT, PHLC_HCF, PHLC_LCF, PHLC_HRC, PHLC_LRC,                   &
                         PRHODJ, PPRES, PZT, PLBDAR, PLBDAR_RF, PLVFACT, PCJ, PKA, PDV, PRF, PCF, PTHT, PTHLT, &
                         PRVS, PRCS, PRRS, PTHS, PUSW, PEVAP3D)
!
!*      0. DECLARATIONS
!          ------------
!
use modd_budget,         only: lbudget_th, lbudget_rv, lbudget_rc, lbudget_rr, &
                               NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RR, &
                               tbudgets
USE MODD_CST,            only: XALPW, XBETAW, XCL, XCPV, XGAMW, XLVTT, XMD, XMV, XRV, XTT
USE MODD_PARAM_ICE_n,      only: CSUBG_RC_RR_ACCR, CSUBG_RR_EVAP
USE MODD_RAIN_ICE_DESCR_n, only: XCEXVT, XRTMIN
USE MODD_RAIN_ICE_PARAM_n, only: X0EVAR, X1EVAR, XCRIAUTC, XEX0EVAR, XEX1EVAR, XEXCACCR, XFCACCR, XTIMAUTC

use mode_budget,         only: Budget_store_add
#ifdef MNH_OPENACC
USE MODE_MNH_ZWORK,      ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
#endif
use mode_mppdb
use MODE_MSG
#if defined(MNH_BITREP) || defined(MNH_BITREP_OMP)
USE MODI_BITREP
#endif

#if defined(MNH_COMPILER_CCE) && defined(MNH_BITREP_OMP)
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,  DIMENSION(:,:,:), intent(in)    :: OMICRO   ! Test where to compute all processes
INTEGER,                    intent(in)    :: KMICRO
INTEGER,  DIMENSION(:),     intent(in)    :: K1
INTEGER,  DIMENSION(:),     intent(in)    :: K2
INTEGER,  DIMENSION(:),     intent(in)    :: K3
REAL,     DIMENSION(:),     intent(in)    :: PRHODREF ! RHO Dry REFerence
REAL,     DIMENSION(:),     intent(in)    :: PRVT     ! Water vapor m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRCT     ! Cloud water m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRRT     ! Rain water m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PHLC_HCF ! HLCLOUDS : fraction of High Cloud Fraction in grid
REAL,     DIMENSION(:),     intent(in)    :: PHLC_LCF ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
REAL,     DIMENSION(:),     intent(in)    :: PHLC_HRC ! HLCLOUDS : LWC that is High LWC in grid
REAL,     DIMENSION(:),     intent(in)    :: PHLC_LRC ! HLCLOUDS : LWC that is Low  LWC in grid
REAL,     DIMENSION(:),     intent(in)    :: PRHODJ   ! RHO times Jacobian
REAL,     DIMENSION(:),     intent(in)    :: PPRES    ! Pressure
REAL,     DIMENSION(:),     intent(in)    :: PZT      ! Temperature
REAL,     DIMENSION(:),     intent(in)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL,     DIMENSION(:),     intent(in)    :: PLBDAR_RF! Slope parameter of the raindrop  distribution
                                                      ! for the Rain Fraction part
REAL,     DIMENSION(:),     intent(in)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL,     DIMENSION(:),     intent(in)    :: PCJ      ! Function to compute the ventilation coefficient
REAL,     DIMENSION(:),     intent(in)    :: PKA      ! Thermal conductivity of the air
REAL,     DIMENSION(:),     intent(in)    :: PDV      ! Diffusivity of water vapor in the air
REAL,     DIMENSION(:),     intent(in)    :: PRF      ! Rain fraction
REAL,     DIMENSION(:),     intent(in)    :: PCF      ! Cloud fraction
REAL,     DIMENSION(:),     intent(in)    :: PTHT     ! Potential temperature
REAL,     DIMENSION(:),     intent(in)    :: PTHLT    ! Liquid potential temperature
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRVS     ! Water vapor m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRCS     ! Cloud water m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRRS     ! Rain water m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PTHS     ! Theta source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PUSW     ! Undersaturation over water
REAL,     DIMENSION(:,:,:), INTENT(INOUT) :: PEVAP3D  ! Rain evap profile
!
!*       0.2  declaration of local variables
!
INTEGER                            :: JL
#ifndef MNH_OPENACC
LOGICAL, DIMENSION(:), ALLOCATABLE :: GWORK
REAL,    DIMENSION(:), ALLOCATABLE :: ZZW  ! Work array
REAL,    DIMENSION(:), ALLOCATABLE :: ZZW2 ! Work array
REAL,    DIMENSION(:), ALLOCATABLE :: ZZW3 ! Work array
REAL,    DIMENSION(:), ALLOCATABLE :: ZZW4 ! Work array
#else
LOGICAL, DIMENSION(:), pointer, contiguous :: GWORK
REAL,    DIMENSION(:), pointer, contiguous :: ZZW  ! Work array
REAL,    DIMENSION(:), pointer, contiguous :: ZZW2 ! Work array
REAL,    DIMENSION(:), pointer, contiguous :: ZZW3 ! Work array
REAL,    DIMENSION(:), pointer, contiguous :: ZZW4 ! Work array
#endif
!
INTEGER :: JLU
!
LOGICAL :: GCSUBG_RR_EVAP ! temporary variable for OpenCC character limitation (Cray CCE)
!
!-------------------------------------------------------------------------------
!
! IN variables
!
!$acc data present( OMICRO, K1, K2, K3, PRHODREF, PRVT, PRCT, PRRT, &
!$acc &             PHLC_HCF, PHLC_LCF, PHLC_HRC, PHLC_LRC,         &
!$acc &             PRHODJ, PPRES, PZT, PLBDAR, PLBDAR_RF, PLVFACT, &
!$acc &             PCJ, PKA, PDV, PRF, PCF, PTHT,PTHLT,            &
!
! INOUT variables
!
!$acc &             PRVS, PRCS, PRRS, PTHS, PUSW, PEVAP3D )
!
! OUT variables
!
!NONE

IF (MPPDB_INITIALIZED) THEN
  !Check all IN arrays
  CALL MPPDB_CHECK(OMICRO,"RAIN_ICE_WARM beg:OMICRO")
  CALL MPPDB_CHECK(PRHODREF,"RAIN_ICE_WARM beg:PRHODREF")
  CALL MPPDB_CHECK(PRVT,"RAIN_ICE_WARM beg:PRVT")
  CALL MPPDB_CHECK(PRCT,"RAIN_ICE_WARM beg:PRCT")
  CALL MPPDB_CHECK(PRRT,"RAIN_ICE_WARM beg:PRRT")
  CALL MPPDB_CHECK(PHLC_HCF,"RAIN_ICE_WARM beg:PHLC_HCF")
  CALL MPPDB_CHECK(PHLC_LCF,"RAIN_ICE_WARM beg:PHLC_LCF")
  CALL MPPDB_CHECK(PHLC_HRC,"RAIN_ICE_WARM beg:PHLC_HRC")
  CALL MPPDB_CHECK(PHLC_LRC,"RAIN_ICE_WARM beg:PHLC_LRC")
  CALL MPPDB_CHECK(PRHODJ,"RAIN_ICE_WARM beg:PRHODJ")
  CALL MPPDB_CHECK(PPRES,"RAIN_ICE_WARM beg:PPRES")
  CALL MPPDB_CHECK(PZT,"RAIN_ICE_WARM beg:PZT")
  CALL MPPDB_CHECK(PLBDAR,"RAIN_ICE_WARM beg:PLBDAR")
  CALL MPPDB_CHECK(PLBDAR_RF,"RAIN_ICE_WARM beg:PLBDAR_RF")
  CALL MPPDB_CHECK(PLVFACT,"RAIN_ICE_WARM beg:PLVFACT")
  CALL MPPDB_CHECK(PCJ,"RAIN_ICE_WARM beg:PCJ")
  CALL MPPDB_CHECK(PKA,"RAIN_ICE_WARM beg:PKA")
  CALL MPPDB_CHECK(PDV,"RAIN_ICE_WARM beg:PDV")
  CALL MPPDB_CHECK(PRF,"RAIN_ICE_WARM beg:PRF")
  CALL MPPDB_CHECK(PCF,"RAIN_ICE_WARM beg:PCF")
  CALL MPPDB_CHECK(PTHT,"RAIN_ICE_WARM beg:PTHT")
  CALL MPPDB_CHECK(PTHLT,"RAIN_ICE_WARM beg:PTHLT")
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PRVS,"RAIN_ICE_WARM beg:PRVS")
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_WARM beg:PRCS")
  CALL MPPDB_CHECK(PRRS,"RAIN_ICE_WARM beg:PRRS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_WARM beg:PTHS")
  CALL MPPDB_CHECK(PUSW,"RAIN_ICE_WARM beg:PUSW")
  CALL MPPDB_CHECK(PEVAP3D,"RAIN_ICE_WARM beg:PEVAP3D")
END IF
!
JLU = size(PRHODREF)
!
#ifndef MNH_OPENACC
ALLOCATE( GWORK(JLU) )
ALLOCATE( ZZW  (JLU) )
ALLOCATE( ZZW2 (JLU) )
ALLOCATE( ZZW3 (JLU) )
ALLOCATE( ZZW4 (JLU) )
#else
!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_WARM' )

CALL MNH_MEM_GET( GWORK, JLU )
CALL MNH_MEM_GET( ZZW,   JLU )
CALL MNH_MEM_GET( ZZW2,  JLU )
CALL MNH_MEM_GET( ZZW3,  JLU )
CALL MNH_MEM_GET( ZZW4,  JLU )

!$acc data present( gwork, zzw, zzw2, zzw3, zzw4 )
#endif
!
!*       4.2    compute the autoconversion of r_c for r_r production: RCAUTR
!
!$acc kernels
!$mnh_do_concurrent(JL=1:JLU)
zzw(JL) = 0.
    GWORK(JL) = PRCS(JL)>0.0 .AND. PHLC_HCF(JL)>0.0
    IF( GWORK(JL) )THEN
      ZZW(JL) = XTIMAUTC*MAX( PHLC_HRC(JL)/PHLC_HCF(JL)  - XCRIAUTC/PRHODREF(JL),0.0)
      ZZW(JL) = MIN( PRCS(JL),PHLC_HCF(JL)*ZZW(JL))
      PRCS(JL) = PRCS(JL) - ZZW(JL)
      PRRS(JL) = PRRS(JL) + ZZW(JL)
    ENDIF
!$mnh_end_do()
!$acc end kernels

    if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'AUTO', &
                                             Unpack( -zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'AUTO', &
                                             Unpack(  zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
!
!*       4.3    compute the accretion of r_c for r_r production: RCACCR
!
    zzw(:) = 0.
    IF (CSUBG_RC_RR_ACCR=='NONE') THEN
!$acc kernels
      !CLoud water and rain are diluted over the grid box
      GWORK(:) = PRCT(:)>XRTMIN(2) .AND. PRRT(:)>XRTMIN(3) .AND. PRCS(:)>0.0
     !$mnh_do_concurrent ( JL=1:JLU )
        IF ( GWORK(JL) )  THEN
           ZZW(JL) = MIN( PRCS(JL), XFCACCR * PRCT(JL) &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)       
                * PLBDAR(JL)**XEXCACCR    &
                * PRHODREF(JL)**(-XCEXVT) )
#else
                * BR_POW(PLBDAR(JL),XEXCACCR)  &
                * BR_POW(PRHODREF(JL),-XCEXVT) )
#endif
           PRCS(JL) = PRCS(JL) - ZZW(JL)
           PRRS(JL) = PRRS(JL) + ZZW(JL)
        END IF
     !$mnh_end_do() ! CONCURRENT
!$acc end kernels

    ELSEIF (CSUBG_RC_RR_ACCR=='PRFR') THEN
#ifdef MNH_OPENACC
      CALL PRINT_MSG(NVERB_WARNING,'GEN','RAIN_ICE_WARM','OPENACC: CSUBG_RC_RR_ACCR=="PRFR" not yet tested')
#endif
!$acc kernels
      !Cloud water is concentrated over its fraction with possibly to parts with high and low content as set for autoconversion
      !Rain is concnetrated over its fraction
      !Rain in high content area fraction: PHLC_HCF
      !Rain in low content area fraction:
      ! if PRF<PCF (rain is entirely falling in cloud): PRF-PHLC_HCF
      ! if PRF>PCF (rain is falling in cloud and in clear sky): PCF-PHLC_HCF
      ! => min(PCF, PRF)-PHLC_HCF
      ZZW(:) = 0.
      GWORK(:) = PHLC_HRC(:)>XRTMIN(2) .AND. PRRT(:)>XRTMIN(3) .AND. PRCS(:)>0.0 .AND. PHLC_HCF(:)>0
     !$mnh_do_concurrent ( JL=1:JLU )
        IF ( GWORK(JL) ) THEN
           !Accretion due to rain falling in high cloud content
           ZZW(JL) = XFCACCR * ( PHLC_HRC(JL)/PHLC_HCF(JL) ) &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)      
                * PLBDAR_RF(JL)**XEXCACCR         &
                * PRHODREF(JL)**(-XCEXVT)         &
#else
                * BR_POW(PLBDAR_RF(JL),XEXCACCR)         &
                * BR_POW(PRHODREF(JL),-XCEXVT)           &
#endif
                * PHLC_HCF(JL)
        END IF
     !$mnh_end_do() ! CONCURRENT
     GWORK(:) = PHLC_LRC(:)>XRTMIN(2) .AND. PRRT(:)>XRTMIN(3) .AND. PRCS(:)>0.0 .AND. PHLC_LCF(:)>0
 !$mnh_do_concurrent ( JL=1:JLU )
    IF ( GWORK(JL) ) THEN     
        !We add acrretion due to rain falling in low cloud content
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
        ZZW(JL) = ZZW(JL) + XFCACCR * ( PHLC_LRC(JL)/PHLC_LCF(JL) ) &
                        * PLBDAR_RF(JL)**XEXCACCR                &
                        * PRHODREF(JL)**(-XCEXVT)                &
                        * (MIN(PCF(JL), PRF(JL))-PHLC_HCF(JL))
#else
        ZZW(JL) = ZZW(JL) + XFCACCR * ( PHLC_LRC(JL)/PHLC_LCF(JL) ) &
                        * BR_POW(PLBDAR_RF(JL),XEXCACCR)         &
                        * BR_POW(PRHODREF(JL),-XCEXVT)           &
                        * (MIN(PCF(JL), PRF(JL))-PHLC_HCF(JL))
#endif
      END IF
      ZZW(JL)=MIN(PRCS(JL), ZZW(JL))
      PRCS(JL) = PRCS(JL) - ZZW(JL)
      PRRS(JL) = PRRS(JL) + ZZW(JL)
   !$mnh_end_do() ! CONCURRENT
!$acc end kernels

    ELSE
      call Print_msg( NVERB_FATAL, 'GEN', 'RAIN_ICE_WARM', 'invalid CSUBG_RC_RR_ACCR value: '//Trim(csubg_rc_rr_accr) )
    ENDIF

    if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'ACCR', &
                                             Unpack( -zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'ACCR', &
                                             Unpack(  zzw(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
!
!*       4.4    compute the evaporation of r_r: RREVAV
!
!$acc kernels
    ZZW(:) = 0.0
!$acc end kernels

    IF (CSUBG_RR_EVAP=='NONE') THEN
!$acc kernels
       !Evaporation only when there's no cloud (RC must be 0)
       GWORK(:) = PRRT(:)>XRTMIN(3) .AND. PRCT(:)<=XRTMIN(2)
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)       
       !$mnh_expand_where(JL=1:JLU)
       WHERE( GWORK(:) )
          ZZW(:)  = EXP( XALPW - XBETAW/PZT(:) - XGAMW*ALOG(PZT(:) ) ) ! es_w
          PUSW(:) = 1.0 - PRVT(:)*( PPRES(:)-ZZW(:) ) / ( (XMV/XMD) * ZZW(:) )
                                                        ! Undersaturation over water
          ZZW(:) = ( XLVTT+(XCPV-XCL)*(PZT(:)-XTT) )**2 / ( PKA(:)*XRV*PZT(:)**2 ) &
               + ( XRV*PZT(:) ) / ( PDV(:)*ZZW(:) )
          ZZW(:) = MIN( PRRS(:),( MAX( 0.0,PUSW(:) )/(PRHODREF(:)*ZZW(:)) ) *      &
            ( X0EVAR*PLBDAR(:)**XEX0EVAR+X1EVAR*PCJ(:)*PLBDAR(:)**XEX1EVAR ) )
          PRRS(:) = PRRS(:) - ZZW(:)
          PRVS(:) = PRVS(:) + ZZW(:)
          PTHS(:) = PTHS(:) - ZZW(:)*PLVFACT(:)
       END WHERE
       !$mnh_end_expand_where()
#else
       !$mnh_do_concurrent ( JL=1:JLU )
          IF ( GWORK(JL) ) THEN
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)       
             ZZW(JL)  = EXP( XALPW - XBETAW/PZT(JL) - XGAMW*LOG(PZT(JL) ) ) ! es_w
#else
             ZZW(JL)  = BR_EXP( XALPW - XBETAW/PZT(JL) - XGAMW*BR_LOG(PZT(JL) ) ) ! es_w
#endif
             PUSW(JL) = 1.0 - PRVT(JL)*( PPRES(JL)-ZZW(JL) ) / ( (XMV/XMD) * ZZW(JL) )
             ! Undersaturation over water
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
             ZZW(JL) = ( XLVTT+(XCPV-XCL)*(PZT(JL)-XTT) )**2 / ( PKA(JL)*XRV*PZT(JL)**2 ) &
#else
             ZZW(JL) = BR_P2( XLVTT+(XCPV-XCL)*(PZT(JL)-XTT) ) / ( PKA(JL)*XRV*BR_P2(PZT(JL)) ) &
#endif
                  + ( XRV*PZT(JL) ) / ( PDV(JL)*ZZW(JL) )
             ZZW(JL) = MIN( PRRS(JL),( MAX( 0.0,PUSW(JL) )/(PRHODREF(JL)*ZZW(JL)) ) *            &
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)       
                  ( X0EVAR*PLBDAR(JL)**XEX0EVAR+X1EVAR*PCJ(JL)*PLBDAR(JL)**XEX1EVAR ) )
#else
                  ( X0EVAR*BR_POW(PLBDAR(JL),XEX0EVAR)+X1EVAR*PCJ(JL)*BR_POW(PLBDAR(JL),XEX1EVAR) ) )
#endif
             PRRS(JL) = PRRS(JL) - ZZW(JL)
             PRVS(JL) = PRVS(JL) + ZZW(JL)
             PTHS(JL) = PTHS(JL) - ZZW(JL)*PLVFACT(JL)
          END IF
       !$mnh_end_do() ! CONCURRENT
#endif       
!$acc end kernels
    ELSEIF (CSUBG_RR_EVAP=='CLFR' .OR. CSUBG_RR_EVAP=='PRFR') THEN
#ifdef MNH_OPENACC
      CALL PRINT_MSG(NVERB_WARNING,'GEN','RAIN_ICE_WARM','OPENACC: CSUBG_RR_EVAP=="CLFR" or "PRFR" not yet tested')
#endif
GCSUBG_RR_EVAP=.false.
IF (CSUBG_RR_EVAP=='CLFR') GCSUBG_RR_EVAP=.true.
!$acc kernels
      !Evaporation in clear sky part
      !With CLFR, rain is diluted over the grid box
      !With PRFR, rain is concentrated in its fraction
      !Use temperature and humidity in clear sky part like Bechtold et al. (1993)
      !IF (CSUBG_RR_EVAP=='CLFR') THEN
      IF (GCSUBG_RR_EVAP) THEN
        ZZW4(:)=1. !Precipitation fraction
        ZZW3(:)=PLBDAR(:)
      ELSE
        ZZW4(:)=PRF(:) !Precipitation fraction
        ZZW3(:)=PLBDAR_RF(:)
      ENDIF

      !ATTENTION
      !Il faudrait recalculer les variables PKA, PDV, PCJ en tenant compte de la température T^u
      !Ces variables devraient être sorties de rain_ice_slow et on mettrait le calcul de T^u, T^s
      !et plusieurs versions (comme actuellement, en ciel clair, en ciel nuageux) de PKA, PDV, PCJ dans rain_ice
      !On utiliserait la bonne version suivant l'option NONE, CLFR... dans l'évaporation et ailleurs
      !$mnh_do_concurrent ( JL=1:JLU )
        GWORK(JL) = PRRT(JL)>XRTMIN(3) .AND. ZZW4(JL)>PCF(JL)
        IF ( GWORK(JL) ) THEN
        ! outside the cloud (environment) the use of T^u (unsaturated) instead of T
        ! Bechtold et al. 1993
        !
        ! T^u = T_l = theta_l * (T/theta)
        ZZW2(JL) =  PTHLT(JL) * PZT(JL) / PTHT(JL)
        !
        ! es_w with new T^u
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
        ZZW(JL)  = EXP( XALPW - XBETAW/ZZW2(JL) - XGAMW*ALOG(ZZW2(JL) ) )
#else
        ZZW(JL)  = BR_EXP( XALPW - XBETAW/ZZW2(JL) - XGAMW*BR_LOG(ZZW2(JL) ) )
#endif
        !
        ! S, Undersaturation over water (with new theta^u)
        PUSW(JL) = 1.0 - PRVT(JL)*( PPRES(JL)-ZZW(JL) ) / ( (XMV/XMD) * ZZW(JL) )
        !
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
        ZZW(JL) = ( XLVTT+(XCPV-XCL)*(ZZW2(JL)-XTT) )**2 / ( PKA(JL)*XRV*ZZW2(JL)**2 ) &
               + ( XRV*ZZW2(JL) ) / ( PDV(JL)*ZZW(JL) )
        !
        ZZW(JL) = MAX( 0.0,PUSW(JL) )/(PRHODREF(JL)*ZZW(JL))  *      &
               ( X0EVAR*ZZW3(JL)**XEX0EVAR+X1EVAR*PCJ(JL)*ZZW3(JL)**XEX1EVAR )
#else
        ZZW(JL) = BR_P2( XLVTT+(XCPV-XCL)*(ZZW2(JL)-XTT) ) / ( PKA(JL)*XRV*BR_P2(ZZW2(JL)) ) &
               + ( XRV*ZZW2(JL) ) / ( PDV(JL)*ZZW(JL) )
        !
        ZZW(JL) = MAX( 0.0,PUSW(JL) )/(PRHODREF(JL)*ZZW(JL))  *      &
               ( X0EVAR*BR_POW(ZZW3(JL),XEX0EVAR)+X1EVAR*PCJ(JL)*BR_POW(ZZW3(JL),XEX1EVAR) )
#endif
        !
        ZZW(JL) = MIN( PRRS(JL),  ZZW(JL) *( ZZW4(JL) - PCF(JL) ) )
        !
        PRRS(JL) = PRRS(JL) - ZZW(JL)
        PRVS(JL) = PRVS(JL) + ZZW(JL)
        PTHS(JL) = PTHS(JL) - ZZW(JL)*PLVFACT(JL)
        END IF
     !$mnh_end_do() ! CONCURRENT
!$acc end kernels

    ELSE
      call Print_msg( NVERB_FATAL, 'GEN', 'RAIN_ICE_WARM', 'invalid CSUBG_RR_EVAP value: '//Trim( csubg_rr_evap ) )
    END IF

    if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'REVA', &
                                             Unpack( -zzw(:) * plvfact(:) * prhodj(:), mask = omicro(:,:,:), field = 0.) )
    if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'REVA', &
                                             Unpack(  zzw(:)              * prhodj(:), mask = omicro(:,:,:), field = 0.) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'REVA', &
                                             Unpack( -zzw(:)              * prhodj(:), mask = omicro(:,:,:), field = 0.) )

!$acc kernels
!$acc loop independent
    DO JL = 1, KMICRO
      PEVAP3D(K1(JL), K2(JL), K3(JL)) = ZZW( JL )
    END DO
!$acc end kernels
!
IF (MPPDB_INITIALIZED) THEN
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PRVS,"RAIN_ICE_WARM end:PRVS")
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_WARM end:PRCS")
  CALL MPPDB_CHECK(PRRS,"RAIN_ICE_WARM end:PRRS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_WARM end:PTHS")
  CALL MPPDB_CHECK(PUSW,"RAIN_ICE_WARM end:PUSW")
  CALL MPPDB_CHECK(PEVAP3D,"RAIN_ICE_WARM end:PEVAP3D")
END IF

!$acc end data

#ifdef MNH_OPENACC
!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'RAIN_ICE_WARM' )
#endif

!$acc end data

  END SUBROUTINE RAIN_ICE_WARM

END MODULE MODE_RAIN_ICE_WARM

