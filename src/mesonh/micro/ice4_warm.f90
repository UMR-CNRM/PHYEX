!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_ICE4_WARM
INTERFACE
SUBROUTINE ICE4_WARM(KSIZE, LDSOFT, PCOMPUTE, HSUBG_RC_RR_ACCR, HSUBG_RR_EVAP, &
                     PRHODREF, PLVFACT, PT, PPRES, PTHT, &
                     PLBDAR, PLBDAR_RF, PKA, PDV, PCJ, &
                     PHLC_LCF, PHLC_HCF, PHLC_LRC, PHLC_HRC, &
                     PCF, PRF, &
                     PRVT, PRCT, PRRT, &
                     PRCAUTR, PRCACCR, PRREVAV, &
                     PA_TH, PA_RV, PA_RC, PA_RR)
IMPLICIT NONE
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
CHARACTER(len=80),            INTENT(IN)    :: HSUBG_RC_RR_ACCR ! subgrid rc-rr accretion
CHARACTER(len=80),            INTENT(IN)    :: HSUBG_RR_EVAP ! subgrid rr evaporation
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PTHT     ! Theta at time t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR_RF!like PLBDAR but for the Rain Fraction part
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_HCF ! HLCLOUDS : fraction of High Cloud Fraction in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_LCF ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_HRC ! HLCLOUDS : LWC that is High LWC in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_LRC ! HLCLOUDS : LWC that is Low  LWC in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCF      ! Cloud fraction
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRF      ! Rain fraction
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCAUTR   ! Autoconversion of r_c for r_r production
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCACCR  ! Accretion of r_c for r_r production
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRREVAV  ! Evaporation of r_r
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RV
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
END SUBROUTINE ICE4_WARM
END INTERFACE
END MODULE MODI_ICE4_WARM
SUBROUTINE ICE4_WARM(KSIZE, LDSOFT, PCOMPUTE, HSUBG_RC_RR_ACCR, HSUBG_RR_EVAP, &
                     PRHODREF, PLVFACT, PT, PPRES, PTHT, &
                     PLBDAR, PLBDAR_RF, PKA, PDV, PCJ, &
                     PHLC_LCF, PHLC_HCF, PHLC_LRC, PHLC_HRC, &
                     PCF, PRF, &
                     PRVT, PRCT, PRRT, &
                     PRCAUTR, PRCACCR, PRREVAV, &
                     PA_TH, PA_RV, PA_RC, PA_RR)
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
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: XALPW,XBETAW,XCL,XCPD,XCPV,XGAMW,XLVTT,XMD,XMV,XRV,XTT,XEPSILO
USE MODD_RAIN_ICE_DESCR, ONLY: XCEXVT,XRTMIN
USE MODD_RAIN_ICE_PARAM, ONLY: X0EVAR,X1EVAR,XCRIAUTC,XEX0EVAR,XEX1EVAR,XEXCACCR,XFCACCR,XTIMAUTC
!
USE MODE_MSG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
CHARACTER(len=80),            INTENT(IN)    :: HSUBG_RC_RR_ACCR ! subgrid rc-rr accretion
CHARACTER(len=80),            INTENT(IN)    :: HSUBG_RR_EVAP ! subgrid rr evaporation
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PTHT     ! Theta at time t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR_RF!like PLBDAR but for the Rain Fraction part
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_HCF ! HLCLOUDS : fraction of High Cloud Fraction in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_LCF ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_HRC ! HLCLOUDS : LWC that is High LWC in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PHLC_LRC ! HLCLOUDS : LWC that is Low  LWC in grid
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCF      ! Cloud fraction
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRF      ! Rain fraction
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCAUTR   ! Autoconversion of r_c for r_r production
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCACCR  ! Accretion of r_c for r_r production
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRREVAV  ! Evaporation of r_r
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RV
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KSIZE) :: ZZW2, ZZW3, ZZW4
REAL, DIMENSION(KSIZE) :: ZUSW ! Undersaturation over water
REAL, DIMENSION(KSIZE) :: ZTHLT    ! Liquid potential temperature
REAL, DIMENSION(KSIZE) :: ZMASK, ZMASK1, ZMASK2
INTEGER :: JL
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
!
!*       4.2    compute the autoconversion of r_c for r_r production: RCAUTR
!
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(2)-PHLC_HRC(JL))) * & ! PHLC_HRC(:)>XRTMIN(2)
           &MAX(0., -SIGN(1., 1.E-20-PHLC_HCF(JL))) * & ! PHLC_HCF(:) .GT. 0.
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRCAUTR(JL)=PRCAUTR(JL)*ZMASK(JL)
  ENDDO
ELSE
  PRCAUTR(:) = 0.
  WHERE(ZMASK(:)==1.)
    PRCAUTR(:) = XTIMAUTC*MAX(PHLC_HRC(:)/PHLC_HCF(:) - XCRIAUTC/PRHODREF(:), 0.0)
    PRCAUTR(:) = PHLC_HCF(:)*PRCAUTR(:)
  END WHERE
ENDIF
DO JL=1, KSIZE
  PA_RC(JL) = PA_RC(JL) - PRCAUTR(JL)
  PA_RR(JL) = PA_RR(JL) + PRCAUTR(JL)
ENDDO
!
!
!*       4.3    compute the accretion of r_c for r_r production: RCACCR
!
IF (HSUBG_RC_RR_ACCR=='NONE') THEN
  !CLoud water and rain are diluted over the grid box
  DO JL=1, KSIZE
    ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(2)-PRCT(JL))) * & ! PRCT(:)>XRTMIN(2)
             &MAX(0., -SIGN(1., XRTMIN(3)-PRRT(JL))) * & ! PRRT(:)>XRTMIN(3)
             &PCOMPUTE(JL)
  ENDDO
  IF(LDSOFT) THEN
    DO JL=1, KSIZE
      PRCACCR(JL)=PRCACCR(JL) * ZMASK(JL)
    ENDDO
  ELSE
    PRCACCR(:) = 0.
    WHERE(ZMASK(:)==1.)
      PRCACCR(:) = XFCACCR * PRCT(:)                &
                 * PLBDAR(:)**XEXCACCR    &
                 * PRHODREF(:)**(-XCEXVT)
    END WHERE
  ENDIF

ELSEIF (HSUBG_RC_RR_ACCR=='PRFR') THEN
  !Cloud water is concentrated over its fraction with possibly to parts with high and low content as set for autoconversion
  !Rain is concnetrated over its fraction
  !Rain in high content area fraction: PHLC_HCF
  !Rain in low content area fraction:
  ! if PRF<PCF (rain is entirely falling in cloud): PRF-PHLC_HCF
  ! if PRF>PCF (rain is falling in cloud and in clear sky): PCF-PHLC_HCF
  ! => min(PCF, PRF)-PHLC_HCF
  DO JL=1, KSIZE
    ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(2)-PRCT(JL))) * & ! PRCT(:)>XRTMIN(2)
             &MAX(0., -SIGN(1., XRTMIN(3)-PRRT(JL))) * & ! PRRT(:)>XRTMIN(3)
             &PCOMPUTE(JL)
    ZMASK1(JL)=ZMASK(JL) * &
              &MAX(0., -SIGN(1., XRTMIN(2)-PHLC_HRC(JL))) * & ! PHLC_HRC(:)>XRTMIN(2)
              &MAX(0., -SIGN(1., 1.E-20-PHLC_HCF(JL))) ! PHLC_HCF(:)>0.
    ZMASK2(JL)=ZMASK(JL) * &
              &MAX(0., -SIGN(1., XRTMIN(2)-PHLC_LRC(JL))) * & ! PHLC_LRC(:)>XRTMIN(2)
              &MAX(0., -SIGN(1., 1.E-20-PHLC_LCF(JL))) ! PHLC_LCF(:)>0.
  ENDDO
  IF(LDSOFT) THEN
    DO JL=1, KSIZE
      PRCACCR(JL)=PRCACCR(JL) * MIN(1., ZMASK1(JL)+ZMASK2(JL))
    ENDDO
  ELSE
    PRCACCR(:)=0.
    WHERE(ZMASK1(:)==1.)
      !Accretion due to rain falling in high cloud content
      PRCACCR(:) = XFCACCR * ( PHLC_HRC(:)/PHLC_HCF(:) )     &
             * PLBDAR_RF(:)**XEXCACCR &
             * PRHODREF(:)**(-XCEXVT) &
             * PHLC_HCF
    END WHERE
    WHERE(ZMASK2(:)==1.)
      !We add acrretion due to rain falling in low cloud content
      PRCACCR(:) = PRCACCR(:) + XFCACCR * ( PHLC_LRC(:)/PHLC_LCF(:) )     &
                      * PLBDAR_RF(:)**XEXCACCR &
                      * PRHODREF(:)**(-XCEXVT) &
                      * (MIN(PCF(:), PRF(:))-PHLC_HCF(:))
    END WHERE
  ENDIF
ELSE
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_WARM','wrong HSUBG_RC_RR_ACCR case')
ENDIF
DO JL=1, KSIZE
  PA_RC(JL) = PA_RC(JL) - PRCACCR(JL)
  PA_RR(JL) = PA_RR(JL) + PRCACCR(JL)
ENDDO
!
!*       4.4    compute the evaporation of r_r: RREVAV
!
IF (HSUBG_RR_EVAP=='NONE') THEN
  DO JL=1, KSIZE
    ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(3)-PRRT(JL))) * & ! PRRT(:)>XRTMIN(3)
             &MAX(0., SIGN(1., XRTMIN(2)-PRCT(JL))) * & ! PRCT(:)<=XRTMIN(2)
             &PCOMPUTE(JL)
  ENDDO
  IF(LDSOFT) THEN
    DO JL=1, KSIZE
      PRREVAV(JL)=PRREVAV(JL)*ZMASK(JL)
    ENDDO
  ELSE
    PRREVAV(:) = 0.
    !Evaporation only when there's no cloud (RC must be 0)
    WHERE(ZMASK(:)==1.)
      PRREVAV(:)  = EXP( XALPW - XBETAW/PT(:) - XGAMW*ALOG(PT(:) ) ) ! es_w
      ZUSW(:) = 1.0 - PRVT(:)*( PPRES(:)-PRREVAV(:) ) / ( XEPSILO * PRREVAV(:) )
                                                    ! Undersaturation over water
      PRREVAV(:) = ( XLVTT+(XCPV-XCL)*(PT(:)-XTT) )**2 / ( PKA(:)*XRV*PT(:)**2 ) &
           + ( XRV*PT(:) ) / ( PDV(:)*PRREVAV(:) )
      PRREVAV(:) = ( MAX( 0.0,ZUSW(:) )/(PRHODREF(:)*PRREVAV(:)) ) *      &
        ( X0EVAR*PLBDAR(:)**XEX0EVAR+X1EVAR*PCJ(:)*PLBDAR(:)**XEX1EVAR )
    END WHERE
  ENDIF

ELSEIF (HSUBG_RR_EVAP=='CLFR' .OR. HSUBG_RR_EVAP=='PRFR') THEN
  !Evaporation in clear sky part
  !With CLFR, rain is diluted over the grid box
  !With PRFR, rain is concentrated in its fraction
  !Use temperature and humidity in clear sky part like Bechtold et al. (1993)
  IF (HSUBG_RR_EVAP=='CLFR') THEN
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
  DO JL=1, KSIZE
    ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(3)-PRRT(JL))) * & ! PRRT(:)>XRTMIN(3)
             &MAX(0., -SIGN(1., PCF(JL)-ZZW4(JL))) * & ! ZZW4(:) > PCF(:)
             &PCOMPUTE(JL)
  ENDDO
  IF(LDSOFT) THEN
    DO JL=1, KSIZE
      PRREVAV(JL)=PRREVAV(JL)*ZMASK(JL)
    ENDDO
  ELSE
    PRREVAV(:) = 0.
    WHERE(ZMASK(:)==1)
      ! outside the cloud (environment) the use of T^u (unsaturated) instead of T
      ! Bechtold et al. 1993
      !
      ! T_l
      ZTHLT(:) = PTHT(:) - XLVTT*PTHT(:)/XCPD/PT(:)*PRCT(:)
      !
      ! T^u = T_l = theta_l * (T/theta)
      ZZW2(:) =  ZTHLT(:) * PT(:) / PTHT(:)
      !
      ! es_w with new T^u
      PRREVAV(:)  = EXP( XALPW - XBETAW/ZZW2(:) - XGAMW*ALOG(ZZW2(:) ) )
      !
      ! S, Undersaturation over water (with new theta^u)
      ZUSW(:) = 1.0 - PRVT(:)*( PPRES(:)-PRREVAV(:) ) / ( XEPSILO * PRREVAV(:) )
      !
      PRREVAV(:) = ( XLVTT+(XCPV-XCL)*(ZZW2(:)-XTT) )**2 / ( PKA(:)*XRV*ZZW2(:)**2 ) &
             + ( XRV*ZZW2(:) ) / ( PDV(:)*PRREVAV(:) )
      !
      PRREVAV(:) = MAX( 0.0,ZUSW(:) )/(PRHODREF(:)*PRREVAV(:))  *      &
             ( X0EVAR*ZZW3(:)**XEX0EVAR+X1EVAR*PCJ(:)*ZZW3(:)**XEX1EVAR )
      !
      PRREVAV(:) = PRREVAV(:)*(ZZW4(:)-PCF(:))
    END WHERE
  ENDIF

ELSE
  CALL PRINT_MSG(NVERB_FATAL,'GEN','ICE4_WARM','wrong HSUBG_RR_EVAP case')
END IF
DO JL=1, KSIZE
  PA_RR(JL) = PA_RR(JL) - PRREVAV(JL)
  PA_RV(JL) = PA_RV(JL) + PRREVAV(JL)
  PA_TH(JL) = PA_TH(JL) - PRREVAV(JL)*PLVFACT(JL)
ENDDO
!
!
END SUBROUTINE ICE4_WARM
