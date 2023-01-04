!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_WARM
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_WARM(CST, ICEP, ICED, KPROMA, KSIZE, LDSOFT, LDCOMPUTE, HSUBG_RC_RR_ACCR, HSUBG_RR_EVAP, &
                    &PRHODREF, PLVFACT, PT, PPRES, PTHT, &
                    &PLBDAR, PLBDAR_RF, PKA, PDV, PCJ, &
                    &PHLC_LCF, PHLC_HCF, PHLC_LRC, PHLC_HRC, &
                    &PCF, PRF, &
                    &PRVT, PRCT, PRRT, &
                    &PRCAUTR, PRCACCR, PRREVAV)
!!
!!**  PURPOSE
!!    -------
!!      Computes the warm process
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the plitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     R. El Khatib 24-Aug-2021 Optimizations
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: CST_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
!
USE MODE_MSG
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
INTEGER,                      INTENT(IN)    :: KPROMA, KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KPROMA),   INTENT(IN)    :: LDCOMPUTE
CHARACTER(LEN=80),            INTENT(IN)    :: HSUBG_RC_RR_ACCR ! subgrid rc-rr accretion
CHARACTER(LEN=80),            INTENT(IN)    :: HSUBG_RR_EVAP ! subgrid rr evaporation
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PTHT     ! Theta at time t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAR_RF!like PLBDAR but for the Rain Fraction part
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PHLC_HCF ! HLCLOUDS : fraction of High Cloud Fraction in grid
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PHLC_LCF ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PHLC_HRC ! HLCLOUDS : LWC that is High LWC in grid
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PHLC_LRC ! HLCLOUDS : LWC that is Low  LWC in grid
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PCF      ! Cloud fraction
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRF      ! Rain fraction
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRCAUTR   ! Autoconversion of r_c for r_r production
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRCACCR  ! Accretion of r_c for r_r production
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRREVAV  ! Evaporation of r_r
!
!*       0.2  declaration of local variables
!
REAL :: ZZW2, ZZW3, ZZW4
REAL, DIMENSION(KPROMA) :: ZUSW ! Undersaturation over water
REAL, DIMENSION(KPROMA) :: ZTHLT    ! Liquid potential temperature
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JL
LOGICAL :: LMASK, LMASK1, LMASK2
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_WARM', 0, ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
!*       4.2    compute the autoconversion of r_c for r_r production: RCAUTR
!
DO JL=1, KSIZE
#ifdef REPRO55
  IF(PHLC_HRC(JL)>ICED%XRTMIN(2) .AND. PHLC_HCF(JL)>1.E-20 .AND. LDCOMPUTE(JL)) THEN
#else
  IF(PHLC_HRC(JL)>ICED%XRTMIN(2) .AND. PHLC_HCF(JL)>0. .AND. LDCOMPUTE(JL)) THEN
#endif
    IF(.NOT. LDSOFT) THEN
#if defined(REPRO48) || defined(REPRO55)
      PRCAUTR(JL) = ICEP%XTIMAUTC*MAX(PHLC_HRC(JL)/PHLC_HCF(JL) - ICEP%XCRIAUTC/PRHODREF(JL), 0.0)
      PRCAUTR(JL) = PHLC_HCF(JL)*PRCAUTR(JL)
#else
      !HCF*autoconv(HRC/HCF) with simplification
      PRCAUTR(JL) = ICEP%XTIMAUTC*MAX(PHLC_HRC(JL) - PHLC_HCF(JL)*ICEP%XCRIAUTC/PRHODREF(JL), 0.0)
#endif
    ENDIF
  ELSE
    PRCAUTR(JL) = 0.
  ENDIF
ENDDO
!
!
!*       4.3    compute the accretion of r_c for r_r production: RCACCR
!
IF (HSUBG_RC_RR_ACCR=='NONE') THEN
  !CLoud water and rain are diluted over the grid box
  DO JL=1, KSIZE
    IF(PRCT(JL)>ICED%XRTMIN(2) .AND. PRRT(JL)>ICED%XRTMIN(3) .AND. LDCOMPUTE(JL)) THEN
      IF(.NOT. LDSOFT) THEN
        PRCACCR(JL) = ICEP%XFCACCR * PRCT(JL)                &
                    & * PLBDAR(JL)**ICEP%XEXCACCR    &
                    & * PRHODREF(JL)**(-ICED%XCEXVT)
      ENDIF
    ELSE
      PRCACCR(JL) = 0.
    ENDIF
  ENDDO

ELSEIF (HSUBG_RC_RR_ACCR=='PRFR') THEN
  !Cloud water is concentrated over its fraction with possibly to parts with high and low content as set for autoconversion
  !Rain is concentrated over its fraction
  !Rain in high content area fraction: PHLC_HCF
  !Rain in low content area fraction:
  ! if PRF<PCF (rain is entirely falling in cloud): PRF-PHLC_HCF
  ! if PRF>PCF (rain is falling in cloud and in clear sky): PCF-PHLC_HCF
  ! => min(PCF, PRF)-PHLC_HCF
  DO JL=1, KSIZE
    LMASK = PRCT(JL)>ICED%XRTMIN(2) .AND. PRRT(JL)>ICED%XRTMIN(3) .AND. LDCOMPUTE(JL)
#ifdef REPRO55
    LMASK1 = LMASK .AND. PHLC_HRC(JL)>ICED%XRTMIN(2) .AND. PHLC_HCF(JL)>1.E-20
#else
    LMASK1 = LMASK .AND. PHLC_HRC(JL)>ICED%XRTMIN(2) .AND. PHLC_HCF(JL)>0.
#endif
#ifdef REPRO48
    LMASK2 = LMASK .AND. PHLC_LRC(JL)>ICED%XRTMIN(2) .AND. PHLC_LCF(JL)>0.
#else
    LMASK2 = LMASK .AND. PHLC_LRC(JL)>ICED%XRTMIN(2) .AND. PHLC_LCF(JL)>1.E-20
#endif
    IF(LMASK1 .OR. LMASK2) THEN
      IF(.NOT. LDSOFT) THEN
        IF(LMASK1) THEN
          !Accretion due to rain falling in high cloud content
#if defined(REPRO48) || defined(REPRO55)
          PRCACCR(JL) = ICEP%XFCACCR * ( PHLC_HRC(JL)/PHLC_HCF(JL) )     &
                      &*PLBDAR_RF(JL)**ICEP%XEXCACCR &
                      &*PRHODREF(JL)**(-ICED%XCEXVT) &
                      &*PHLC_HCF(JL)
#else
          !HCF*accretion(HRC/HCF) with simplification
          PRCACCR(:) = ICEP%XFCACCR * PHLC_HRC(JL)     &
                      &*PLBDAR_RF(JL)**ICEP%XEXCACCR &
                      &*PRHODREF(JL)**(-ICED%XCEXVT)
#endif
        ELSE
          PRCACCR(JL)=0.
        ENDIF
        IF(LMASK2) THEN
          !We add acrretion due to rain falling in low cloud content
          PRCACCR(JL) = PRCACCR(JL) + ICEP%XFCACCR * ( PHLC_LRC(JL)/PHLC_LCF(JL) )     &
                      &*PLBDAR_RF(JL)**ICEP%XEXCACCR &
                      &*PRHODREF(JL)**(-ICED%XCEXVT) &
                      &*(MIN(PCF(JL), PRF(JL))-PHLC_HCF(JL))
        ENDIF
      ENDIF
    ELSE
      PRCACCR(JL)=0.
    ENDIF
  ENDDO
ELSE
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_WARM','wrong HSUBG_RC_RR_ACCR case')
ENDIF
!
!*       4.4    compute the evaporation of r_r: RREVAV
!
IF (HSUBG_RR_EVAP=='NONE') THEN
  DO JL=1, KSIZE
    IF(PRRT(JL)>ICED%XRTMIN(3) .AND. PRCT(JL)<=ICED%XRTMIN(2) .AND. LDCOMPUTE(JL)) THEN
      IF(.NOT. LDSOFT) THEN
        PRREVAV(JL) = EXP(CST%XALPW - CST%XBETAW/PT(JL) - CST%XGAMW*ALOG(PT(JL))) ! es_w
        ZUSW(JL) = 1. - PRVT(JL)*(PPRES(JL)-PRREVAV(JL)) / (CST%XEPSILO * PRREVAV(JL)) ! Undersaturation over water
        PRREVAV(JL) = (CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(JL)-CST%XTT) )**2 / (PKA(JL)*CST%XRV*PT(JL)**2) &
                    &+(CST%XRV*PT(JL)) / (PDV(JL)*PRREVAV(JL))
        PRREVAV(JL) = (MAX(0.,ZUSW(JL) )/(PRHODREF(JL)*PRREVAV(JL)) ) * &
                    & (ICEP%X0EVAR*PLBDAR(JL)**ICEP%XEX0EVAR+ICEP%X1EVAR*PCJ(JL)*PLBDAR(JL)**ICEP%XEX1EVAR)
      ENDIF
    ELSE
      PRREVAV(JL)=0.
    ENDIF
  ENDDO

ELSEIF (HSUBG_RR_EVAP=='CLFR' .OR. HSUBG_RR_EVAP=='PRFR') THEN
  !ATTENTION
  !Il faudrait recalculer les variables PKA, PDV, PCJ en tenant compte de la température T^u
  !Ces variables devraient être sorties de rain_ice_slow et on mettrait le calcul de T^u, T^s
  !et plusieurs versions (comme actuellement, en ciel clair, en ciel nuageux) de PKA, PDV, PCJ dans rain_ice
  !On utiliserait la bonne version suivant l'option NONE, CLFR... dans l'évaporation et ailleurs

  DO JL=1, KSIZE
    !Evaporation in clear sky part
    !With CLFR, rain is diluted over the grid box
    !With PRFR, rain is concentrated in its fraction
    !Use temperature and humidity in clear sky part like Bechtold et al. (1993)
    IF (HSUBG_RR_EVAP=='CLFR') THEN
      ZZW4=1. !Precipitation fraction
      ZZW3=PLBDAR(JL)
    ELSE
      ZZW4=PRF(JL) !Precipitation fraction
      ZZW3=PLBDAR_RF(JL)
    ENDIF

    IF(PRRT(JL)>ICED%XRTMIN(3) .AND. ZZW4>PCF(JL) .AND. LDCOMPUTE(JL)) THEN
      IF(.NOT. LDSOFT) THEN
        ! outside the cloud (environment) the use of T^u (unsaturated) instead of T
        ! Bechtold et al. 1993
        !
        ! T_l
        ZTHLT(JL) = PTHT(JL) - CST%XLVTT*PTHT(JL)/CST%XCPD/PT(JL)*PRCT(JL)
        !
        ! T^u = T_l = theta_l * (T/theta)
        ZZW2 =  ZTHLT(JL) * PT(JL) / PTHT(JL)
        !
        ! es_w with new T^u
        PRREVAV(JL)  = EXP(CST%XALPW - CST%XBETAW/ZZW2 - CST%XGAMW*ALOG(ZZW2))
        !
        ! S, Undersaturation over water (with new theta^u)
        ZUSW(JL) = 1.0 - PRVT(JL)*(PPRES(JL)-PRREVAV(JL)) / (CST%XEPSILO * PRREVAV(JL))
        !
        PRREVAV(JL) = (CST%XLVTT+(CST%XCPV-CST%XCL)*(ZZW2-CST%XTT))**2 / (PKA(JL)*CST%XRV*ZZW2**2) &
                    &+(CST%XRV*ZZW2) / (PDV(JL)*PRREVAV(JL))
        !
        PRREVAV(JL) = MAX(0., ZUSW(JL))/(PRHODREF(JL)*PRREVAV(JL))  *      &
                    & (ICEP%X0EVAR*ZZW3**ICEP%XEX0EVAR+ICEP%X1EVAR*PCJ(JL)*ZZW3**ICEP%XEX1EVAR)
        !
        PRREVAV(JL) = PRREVAV(JL)*(ZZW4-PCF(JL))
      ENDIF
    ELSE
      PRREVAV(JL)=0.
    ENDIF
  ENDDO

ELSE
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_WARM','wrong HSUBG_RR_EVAP case')
END IF
!
IF (LHOOK) CALL DR_HOOK('ICE4_WARM', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_WARM
END MODULE MODE_ICE4_WARM
