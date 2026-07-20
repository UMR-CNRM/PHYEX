!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_WARM
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_WARM(CST, PARAMI, ICEP, ICED, D, LDSOFT, LDCOMPUTE, HSUBG_RC_RR_ACCR, HSUBG_RR_EVAP, &
                    &PRHODREF, PLVFACT, PT, PPRES, PTHT, &
                    &PLBDAR, PLBDAR_RF, PKA, PDV, PCJ, &
                    &PHLC_LCF, PHLC_HCF, PHLC_LRC, PHLC_HRC, &
                    &PCF, PRF, &
                    &PRVT, PRCT, PRRT, PCONC, PACRF, &
                    &PRCAUTR, PRCACCR, PRREVAV)

!$ACDC singlecolumn

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
!      K.I Ivarsson 10-2018: Possibility to use Kogan autoconversion and alternative collision factor
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
!
USE MODE_MSG, ONLY: NVERB_FATAL, PRINT_MSG
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),        INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
TYPE(DIMPHYEX_t),             INTENT(IN)    :: D
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(D%NIJT),   INTENT(IN)    :: LDCOMPUTE
CHARACTER(LEN=80),            INTENT(IN)    :: HSUBG_RC_RR_ACCR ! subgrid rc-rr accretion
CHARACTER(LEN=80),            INTENT(IN)    :: HSUBG_RR_EVAP ! subgrid rr evaporation
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PTHT     ! Theta at time t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAR_RF!like PLBDAR but for the Rain Fraction part
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PHLC_HCF ! HLCLOUDS : fraction of High Cloud Fraction in grid
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PHLC_LCF ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PHLC_HRC ! HLCLOUDS : LWC that is High LWC in grid
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PHLC_LRC ! HLCLOUDS : LWC that is Low  LWC in grid
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCF      ! Cloud fraction
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRF      ! Rain fraction
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCONC    ! Cloud Droplet number concentration
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PACRF    ! Collision factor cloud droplet-rain
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRCAUTR  ! Autoconversion of r_c for r_r production
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRCACCR  ! Accretion of r_c for r_r production
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRREVAV  ! Evaporation of r_r
!
!*       0.2  declaration of local variables
!
REAL :: ZZW2, ZZW3, ZZW4
REAL, DIMENSION(D%NIJT) :: ZUSW ! Undersaturation over water
REAL, DIMENSION(D%NIJT) :: ZTHLT    ! Liquid potential temperature
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JIJ
LOGICAL :: LMASK, LMASK1, LMASK2
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_WARM', 0, ZHOOK_HANDLE)
!
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
#endif
!
!-------------------------------------------------------------------------------
!
!*       4.2    compute the autoconversion of r_c for r_r production: RCAUTR
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(LDCOMPUTE(JIJ)) THEN
    IF(PARAMI%LKOGAN .AND. .NOT. LDSOFT) THEN
      PRCAUTR(JIJ) = 0.
      IF (PRCT(JIJ) >  1.0E-8) THEN
        PRCAUTR(JIJ) =  1350.0 * ICEP%XFRMIN(10)* PCONC(JIJ)**(-1.79) * &
            &  (PRCT(JIJ)/(MAX(ICEP%XFRMIN(11),PCF(JIJ))))**2.47
        PRCAUTR(JIJ) = PRCAUTR(JIJ)*MAX(ICEP%XFRMIN(11),PCF(JIJ))
      ENDIF
    ELSE
      IF(PHLC_HRC(JIJ)>ICED%XRTMIN(2) .AND. PHLC_HCF(JIJ)>0. ) THEN
        IF (.NOT. LDSOFT)THEN
          !HCF*autoconv(HRC/HCF) with simplification
          PRCAUTR(JIJ) = ICEP%XTIMAUTC*MAX(PHLC_HRC(JIJ) - PHLC_HCF(JIJ)*ICEP%XCRIAUTC/PRHODREF(JIJ), 0.0)
        ENDIF
      ELSE
        PRCAUTR(JIJ) = 0.
      ENDIF
    ENDIF
  ELSE
    PRCAUTR(JIJ) = 0.
  ENDIF
ENDDO
!$mnh_end_do()

!
!
!*       4.3    compute the accretion of r_c for r_r production: RCACCR
!
IF (HSUBG_RC_RR_ACCR=='NONE') THEN
  !CLoud water and rain are diluted over the grid box

  !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
  DO JIJ=D%NIJB, D%NIJE
    IF(PRCT(JIJ)>ICED%XRTMIN(2) .AND. PRRT(JIJ)>ICED%XRTMIN(3) .AND. LDCOMPUTE(JIJ)) THEN
      IF(.NOT. LDSOFT) THEN
        PRCACCR(JIJ) = ICEP%XFCACCR * PRCT(JIJ) * PACRF(JIJ)   &
                    & * PLBDAR(JIJ)**ICEP%XEXCACCR    &
                    & * PRHODREF(JIJ)**(-ICED%XCEXVT)
      ENDIF
    ELSE
      PRCACCR(JIJ) = 0.
    ENDIF
  ENDDO
  !$mnh_end_do()

ELSEIF (HSUBG_RC_RR_ACCR=='PRFR') THEN
  !Cloud water is concentrated over its fraction with possibly to parts with high and low content as set for autoconversion
  !Rain is concentrated over its fraction
  !Rain in high content area fraction: PHLC_HCF
  !Rain in low content area fraction:
  ! if PRF<PCF (rain is entirely falling in cloud): PRF-PHLC_HCF
  ! if PRF>PCF (rain is falling in cloud and in clear sky): PCF-PHLC_HCF
  ! => min(PCF, PRF)-PHLC_HCF

  !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
  DO JIJ=D%NIJB, D%NIJE
    LMASK = PRCT(JIJ)>ICED%XRTMIN(2) .AND. PRRT(JIJ)>ICED%XRTMIN(3) .AND. LDCOMPUTE(JIJ)
    LMASK1 = LMASK .AND. PHLC_HRC(JIJ)>ICED%XRTMIN(2) .AND. PHLC_HCF(JIJ)>0.
    LMASK2 = LMASK .AND. PHLC_LRC(JIJ)>ICED%XRTMIN(2) .AND. PHLC_LCF(JIJ)>1.E-20
    IF(LMASK1 .OR. LMASK2) THEN
      IF(.NOT. LDSOFT) THEN
        IF(LMASK1) THEN
          !Accretion due to rain falling in high cloud content
          !HCF*accretion(HRC/HCF) with simplification
          PRCACCR(JIJ) = ICEP%XFCACCR * PHLC_HRC(JIJ)     &
                      &*PLBDAR_RF(JIJ)**ICEP%XEXCACCR &
                      &*PRHODREF(JIJ)**(-ICED%XCEXVT)
        ELSE
          PRCACCR(JIJ)=0.
        ENDIF
        IF(LMASK2) THEN
          !We add acrretion due to rain falling in low cloud content
          PRCACCR(JIJ) = PRCACCR(JIJ) + ICEP%XFCACCR * ( PHLC_LRC(JIJ)/PHLC_LCF(JIJ) )     &
                      &*PLBDAR_RF(JIJ)**ICEP%XEXCACCR &
                      &*PRHODREF(JIJ)**(-ICED%XCEXVT) &
                      &*(MIN(PCF(JIJ), PRF(JIJ))-PHLC_HCF(JIJ))
        ENDIF
      ENDIF
    ELSE
      PRCACCR(JIJ)=0.
    ENDIF
  ENDDO
  !$mnh_end_do()

ELSE
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_WARM','wrong HSUBG_RC_RR_ACCR case')
ENDIF
!
!*       4.4    compute the evaporation of r_r: RREVAV
!
IF (HSUBG_RR_EVAP=='NONE') THEN

  !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
  DO JIJ=D%NIJB, D%NIJE
    IF(PRRT(JIJ)>ICED%XRTMIN(3) .AND. PRCT(JIJ)<=ICED%XRTMIN(2) .AND. LDCOMPUTE(JIJ)) THEN
      IF(.NOT. LDSOFT) THEN
        PRREVAV(JIJ) = EXP(CST%XALPW - CST%XBETAW/PT(JIJ) - CST%XGAMW*LOG(PT(JIJ))) ! es_w
        ZUSW(JIJ) = 1. - PRVT(JIJ)*(PPRES(JIJ)-PRREVAV(JIJ)) / (CST%XEPSILO * PRREVAV(JIJ)) ! Undersaturation over water
        PRREVAV(JIJ) = (CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(JIJ)-CST%XTT) )**2 / (PKA(JIJ)*CST%XRV*PT(JIJ)**2) &
                    &+(CST%XRV*PT(JIJ)) / (PDV(JIJ)*PRREVAV(JIJ))
        PRREVAV(JIJ) = (MAX(0.,ZUSW(JIJ) )/(PRHODREF(JIJ)*PRREVAV(JIJ)) ) * &
                    & (ICEP%X0EVAR*PLBDAR(JIJ)**ICEP%XEX0EVAR+ICEP%X1EVAR*PCJ(JIJ)*PLBDAR(JIJ)**ICEP%XEX1EVAR)
      ENDIF
    ELSE
      PRREVAV(JIJ)=0.
    ENDIF
  ENDDO
  !$mnh_end_do()


ELSEIF (HSUBG_RR_EVAP=='CLFR' .OR. HSUBG_RR_EVAP=='PRFR') THEN
  !ATTENTION
  !Il faudrait recalculer les variables PKA, PDV, PCJ en tenant compte de la température T^u
  !Ces variables devraient être sorties de rain_ice_slow et on mettrait le calcul de T^u, T^s
  !et plusieurs versions (comme actuellement, en ciel clair, en ciel nuageux) de PKA, PDV, PCJ dans rain_ice
  !On utiliserait la bonne version suivant l'option NONE, CLFR... dans l'évaporation et ailleurs


  !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
  DO JIJ=D%NIJB, D%NIJE
    !Evaporation in clear sky part
    !With CLFR, rain is diluted over the grid box
    !With PRFR, rain is concentrated in its fraction
    !Use temperature and humidity in clear sky part like Bechtold et al. (1993)
    IF (HSUBG_RR_EVAP=='CLFR') THEN
      ZZW4=1. !Precipitation fraction
      ZZW3=PLBDAR(JIJ)
    ELSE
      ZZW4=PRF(JIJ) !Precipitation fraction
      ZZW3=PLBDAR_RF(JIJ)
    ENDIF

    IF(PRRT(JIJ)>ICED%XRTMIN(3) .AND. ZZW4>PCF(JIJ) .AND. LDCOMPUTE(JIJ)) THEN
      IF(.NOT. LDSOFT) THEN
        ! outside the cloud (environment) the use of T^u (unsaturated) instead of T
        ! Bechtold et al. 1993
        !
        ! T_l
        ZTHLT(JIJ) = PTHT(JIJ) - CST%XLVTT*PTHT(JIJ)/CST%XCPD/PT(JIJ)*PRCT(JIJ)
        !
        ! T^u = T_l = theta_l * (T/theta)
        ZZW2 =  ZTHLT(JIJ) * PT(JIJ) / PTHT(JIJ)
        !
        ! es_w with new T^u
        PRREVAV(JIJ)  = EXP(CST%XALPW - CST%XBETAW/ZZW2 - CST%XGAMW*LOG(ZZW2))
        !
        ! S, Undersaturation over water (with new theta^u)
        ZUSW(JIJ) = 1.0 - PRVT(JIJ)*(PPRES(JIJ)-PRREVAV(JIJ)) / (CST%XEPSILO * PRREVAV(JIJ))
        !
        PRREVAV(JIJ) = (CST%XLVTT+(CST%XCPV-CST%XCL)*(ZZW2-CST%XTT))**2 / (PKA(JIJ)*CST%XRV*ZZW2**2) &
                    &+(CST%XRV*ZZW2) / (PDV(JIJ)*PRREVAV(JIJ))
        !
        PRREVAV(JIJ) = MAX(0., ZUSW(JIJ))/(PRHODREF(JIJ)*PRREVAV(JIJ))  *      &
                    & (ICEP%X0EVAR*ZZW3**ICEP%XEX0EVAR+ICEP%X1EVAR*PCJ(JIJ)*ZZW3**ICEP%XEX1EVAR)
        !
        PRREVAV(JIJ) = PRREVAV(JIJ)*(ZZW4-PCF(JIJ))
      ENDIF
    ELSE
      PRREVAV(JIJ)=0.
    ENDIF
  ENDDO
  !$mnh_end_do()

ELSE
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_WARM','wrong HSUBG_RR_EVAP case')
END IF
!
IF (LHOOK) CALL DR_HOOK('ICE4_WARM', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_WARM
END MODULE MODE_ICE4_WARM
