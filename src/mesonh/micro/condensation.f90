!MNH_LIC Copyright 2002-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
      MODULE MODI_CONDENSATION
!     ########################
!
INTERFACE
!
       SUBROUTINE CONDENSATION( KIU, KJU, KKU, KIB, KIE, KJB, KJE, KKB, KKE, KKL,&
          HFRAC_ICE, HCONDENS, HLAMBDA3, &
          PPABS, PZZ, PRHODREF, PT, PRV, PRC, PRI, PRS, PRG, PSIGS, PMFCONV, PCLDFR, PSIGRC, OUSERI,&
          OSIGMAS, PSIGQSAT, PLV, PLS, PCPH, PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF)
!
INTEGER,                      INTENT(IN)    :: KIU    ! horizontal dimension in x
INTEGER,                      INTENT(IN)    :: KJU    ! horizontal dimension in y
INTEGER,                      INTENT(IN)    :: KKU    ! vertical dimension
INTEGER,                      INTENT(IN)    :: KIB    ! value of the first point in x
INTEGER,                      INTENT(IN)    :: KIE    ! value of the last  point in x
INTEGER,                      INTENT(IN)    :: KJB    ! value of the first point in y
INTEGER,                      INTENT(IN)    :: KJE    ! value of the last  point in y
INTEGER,                      INTENT(IN)    :: KKB    ! value of the first point in z
INTEGER,                      INTENT(IN)    :: KKE    ! value of the last  point in z
INTEGER,                      INTENT(IN)    :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
CHARACTER(len=1),             INTENT(IN)    :: HFRAC_ICE
CHARACTER(len=4),             INTENT(IN)    :: HCONDENS
CHARACTER(len=*),             INTENT(IN)    :: HLAMBDA3 ! formulation for lambda3 coeff
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PPABS  ! pressure (Pa)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PZZ    ! height of model levels (m)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRHODREF
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PT     ! grid scale T  (K)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRC    ! grid scale r_c mixing ratio (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRI    ! grid scale r_i (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRS    ! grid scale mixing ration of snow (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRG    ! grid scale mixing ration of graupel (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PSIGS  ! Sigma_s from turbulence scheme
REAL, DIMENSION(:,:,:),       INTENT(IN)    :: PMFCONV! convective mass flux (kg /s m^2)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PCLDFR ! cloud fraction
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PSIGRC ! s r_c / sig_s^2
LOGICAL, INTENT(IN)                         :: OUSERI ! logical switch to compute both
                                                      ! liquid and solid condensate (OUSERI=.TRUE.)
                                                      ! or only solid condensate (OUSERI=.FALSE.)
LOGICAL, INTENT(IN)                         :: OSIGMAS! use present global Sigma_s values
                                                      ! or that from turbulence scheme
REAL, INTENT(IN)                            :: PSIGQSAT ! use an extra "qsat" variance contribution (OSIGMAS case)
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PLV
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PLS
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PCPH
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLC_HRC   !cloud water content in precipitating part
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLC_HCF   !precipitating part
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLI_HRI   !
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLI_HCF   !

END SUBROUTINE CONDENSATION
!
END INTERFACE
!
END MODULE MODI_CONDENSATION
!     ######spl
    SUBROUTINE CONDENSATION( KIU, KJU, KKU, KIB, KIE, KJB, KJE, KKB, KKE, KKL,         &
       HFRAC_ICE, HCONDENS, HLAMBDA3,                                                  &
       PPABS, PZZ, PRHODREF, PT, PRV, PRC, PRI, PRS, PRG, PSIGS, PMFCONV, PCLDFR, PSIGRC, OUSERI,&
       OSIGMAS, PSIGQSAT, PLV, PLS, PCPH, PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF )
!   ################################################################################
!
!!
!!    PURPOSE
!!    -------
!!**  Routine to diagnose cloud fraction, liquid and ice condensate mixing ratios
!!    and s'rl'/sigs^2
!!
!!
!!**  METHOD
!!    ------
!!    Based on the large-scale fields of temperature, water vapor, and possibly
!!    liquid and solid condensate, the conserved quantities r_t and h_l are constructed
!!    and then fractional cloudiness, liquid and solid condensate is diagnosed.
!!
!!    The total variance is parameterized as the sum of  stratiform/turbulent variance
!!    and a convective variance.
!!    The turbulent variance is parameterized as a function of first-order moments, and
!!    the convective variance is modelled as a function of the convective mass flux
!!    (units kg/s m^2) as provided by the  mass flux convection scheme.
!!
!!    Nota: if the host model does not use prognostic values for liquid and solid condensate
!!    or does not provide a convective mass flux, put all these values to zero.
!!
!!
!!    EXTERNAL
!!    --------
!!      INI_CST
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST       : contains physical constants
!!
!!    REFERENCE
!!    ---------
!!      Chaboureau J.P. and P. Bechtold (J. Atmos. Sci. 2002)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original: 31.1.2002
!!     modified : 21.3.2002
!!     S.Malardel : 05.2006 : Correction sur le calcul de la fonction de
!!                                         Bougeault F2
!!     W. de Rooy: 06-06-2010: Modification in the statistical cloud scheme
!!                             more specifically adding a variance term
!!                             following ideas of Lenderink & Siebesma 2002
!!                             and adding a height dependence
!!     S. Riette, 18 May 2010 : PSIGQSAT is added
!!     S. Riette, 11 Oct 2011 : MIN function in PDF for continuity
!!                              modification of minimum value for Rc+Ri to create cloud and minimum value for sigma
!!                              Use of guess point as a starting point instead of liquid point
!!                              Better computation of ZCPH and dRsat/dT
!!                              Set ZCOND to zero if PCLDFR==0
!!                              Safety limitation to .99*Pressure for saturation vapour pressure
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!      2015   C.Lac   Change min value of ZSIGMA to be in agreement with AROME
!!      2016   G.Delautier   Restore min value of ZSIGMA (instability)
!!      2016   S.Riette Change INQ1 
!!      2016-11 S. Riette: use HFRAC_ICE, output adjusted state
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_RAIN_ICE_PARAM, ONLY : XCRIAUTC, XCRIAUTI, XACRIAUTI, XBCRIAUTI
!
use mode_msg
!
USE MODI_COMPUTE_FRAC_ICE
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                      INTENT(IN)    :: KIU    ! horizontal dimension in x
INTEGER,                      INTENT(IN)    :: KJU    ! horizontal dimension in y
INTEGER,                      INTENT(IN)    :: KKU    ! vertical dimension
INTEGER,                      INTENT(IN)    :: KIB    ! value of the first point in x
INTEGER,                      INTENT(IN)    :: KIE    ! value of the last  point in x
INTEGER,                      INTENT(IN)    :: KJB    ! value of the first point in y
INTEGER,                      INTENT(IN)    :: KJE    ! value of the last  point in y
INTEGER,                      INTENT(IN)    :: KKB    ! value of the first point in z
INTEGER,                      INTENT(IN)    :: KKE    ! value of the last  point in z
INTEGER,                      INTENT(IN)    :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
CHARACTER(len=1),             INTENT(IN)    :: HFRAC_ICE
CHARACTER(len=4),             INTENT(IN)    :: HCONDENS
CHARACTER(len=*),             INTENT(IN)    :: HLAMBDA3 ! formulation for lambda3 coeff
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PPABS  ! pressure (Pa)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PZZ    ! height of model levels (m)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRHODREF
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PT     ! grid scale T  (K)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRC    ! grid scale r_c mixing ratio (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRI    ! grid scale r_i (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRS    ! grid scale mixing ration of snow (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRG    ! grid scale mixing ration of graupel (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PSIGS  ! Sigma_s from turbulence scheme
REAL, DIMENSION(:,:,:),       INTENT(IN)    :: PMFCONV! convective mass flux (kg /s m^2)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PCLDFR ! cloud fraction
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PSIGRC ! s r_c / sig_s^2

REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PLV    ! Latent heat L_v
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PLS    ! Latent heat L_s
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PCPH   ! Specific heat C_ph
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLC_HRC
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLC_HCF ! cloud fraction
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLI_HRI
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLI_HCF
LOGICAL, INTENT(IN)                         :: OUSERI ! logical switch to compute both
                                                      ! liquid and solid condensate (OUSERI=.TRUE.)
                                                      ! or only solid condensate (OUSERI=.FALSE.)
LOGICAL, INTENT(IN)                         :: OSIGMAS! use present global Sigma_s values
                                                      ! or that from turbulence scheme
REAL, INTENT(IN)                            :: PSIGQSAT ! use an extra "qsat" variance contribution (OSIGMAS case)

!
!
!*       0.2   Declarations of local variables :
!
INTEGER  :: JI, JJ, JK, JKP, JKM, IKTB, IKTE    ! loop index
REAL, DIMENSION(KIU,KJU,KKU) :: ZTLK, ZRT       ! work arrays for T_l and total water mixing ratio
REAL, DIMENSION(KIU,KJU,KKU) :: ZL              ! length scale
REAL, DIMENSION(KIU,KJU,KKU) :: ZFRAC           ! Ice fraction
REAL, DIMENSION(KIU,KJU,KKU) :: ZCRIAUTI        !
INTEGER, DIMENSION(KIU,KJU)  :: ITPL            ! top levels of troposphere 
REAL,    DIMENSION(KIU,KJU)  :: ZTMIN           ! minimum Temp. related to ITPL
!
REAL, DIMENSION(KIU,KJU,KKU) :: ZLV, ZLS, ZCPD
REAL, DIMENSION(KIU,KJU,KKU) :: ZCOND
REAL :: ZGCOND, ZSBAR, ZSBARC, ZQ1, ZAUTC, ZAUTI, ZGAUV, ZGAUC, ZGAUI, ZGAUTC, ZGAUTI   ! Used for integration in Gaussian Probability Density Function
REAL :: ZTEMP, ZPV, ZQSL, ZPIV, ZQSI, ZLVS ! thermodynamics
REAL :: ZLL, DZZ, ZZZ                           ! used for length scales 
REAL :: ZAH, ZA, ZB, ZSIGMA, ZDRW, ZDTL, ZSIG_CONV ! related to computation of Sig_s
REAL :: ZRCOLD, ZRIOLD
INTEGER  :: INQ1
REAL :: ZINC
!
!*       0.3  Definition of constants :
!
!-------------------------------------------------------------------------------
!
REAL,PARAMETER :: ZL0     = 600.        ! tropospheric length scale
REAL,PARAMETER :: ZCSIGMA = 0.2         ! constant in sigma_s parameterization
REAL,PARAMETER :: ZCSIG_CONV = 0.30E-2  ! scaling factor for ZSIG_CONV as function of mass flux
!

REAL, DIMENSION(-22:11),PARAMETER :: ZSRC_1D =(/                         &
       0.           ,  0.           ,  2.0094444E-04,   0.316670E-03,    &
       4.9965648E-04,  0.785956E-03 ,  1.2341294E-03,   0.193327E-02,    &
       3.0190963E-03,  0.470144E-02 ,  7.2950651E-03,   0.112759E-01,    &
       1.7350994E-02,  0.265640E-01 ,  4.0427860E-02,   0.610997E-01,    &
       9.1578111E-02,  0.135888E+00 ,  0.1991484    ,   0.230756E+00,    &
       0.2850565    ,  0.375050E+00 ,  0.5000000    ,   0.691489E+00,    &
       0.8413813    ,  0.933222E+00 ,  0.9772662    ,   0.993797E+00,    &
       0.9986521    ,  0.999768E+00 ,  0.9999684    ,   0.999997E+00,    &
       1.0000000    ,  1.000000     /)
!
!-------------------------------------------------------------------------------
!
!

IKTB=1+JPVEXT
IKTE=KKU-JPVEXT

PCLDFR(:,:,:) = 0. ! Initialize values
PSIGRC(:,:,:) = 0. ! Initialize values
!
!
!-------------------------------------------------------------------------------
! store total water mixing ratio
DO JK=IKTB,IKTE
  DO JJ=KJB,KJE
    DO JI=KIB,KIE
      ZRT(JI,JJ,JK)  = PRV(JI,JJ,JK) + PRC(JI,JJ,JK) + PRI(JI,JJ,JK)
    END DO
  END DO
END DO
!-------------------------------------------------------------------------------
! Preliminary calculations
! latent heat of vaporisation/sublimation
IF(PRESENT(PLV) .AND. PRESENT(PLS)) THEN
  ZLV(:,:,:)=PLV(:,:,:)
  ZLS(:,:,:)=PLS(:,:,:)
ELSE
  DO JK=IKTB,IKTE
    DO JJ=KJB,KJE
      DO JI=KIB,KIE
        ZTEMP  = PT(JI,JJ,JK)
        ! latent heat of vaporisation/sublimation
        ZLV(JI,JJ,JK) = XLVTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
        ZLS(JI,JJ,JK) = XLSTT + ( XCPV - XCI ) * ( ZTEMP - XTT )
      ENDDO
    ENDDO
  ENDDO
ENDIF
IF(PRESENT(PCPH)) THEN
  ZCPD(:,:,:)=PCPH(:,:,:)
ELSE
  DO JK=IKTB,IKTE
    DO JJ=KJB,KJE
      DO JI=KIB,KIE
        ZCPD(JI,JJ,JK) = XCPD + XCPV*PRV(JI,JJ,JK) + XCL*PRC(JI,JJ,JK) + XCI*PRI(JI,JJ,JK) + &
                                XCI*(PRS(JI,JJ,JK) + PRG(JI,JJ,JK) )
      ENDDO
    ENDDO
  ENDDO
ENDIF
!-------------------------------------------------------------------------------
! Preliminary calculations needed for computing the "turbulent part" of Sigma_s
IF ( .NOT. OSIGMAS ) THEN
  DO JK=IKTB,IKTE
    DO JJ=KJB,KJE
      DO JI=KIB,KIE
        ZTEMP  = PT(JI,JJ,JK)
        ! store temperature at saturation
        ZTLK(JI,JJ,JK) = ZTEMP - ZLV(JI,JJ,JK)*PRC(JI,JJ,JK)/ZCPD(JI,JJ,JK) &
                               - ZLS(JI,JJ,JK)*PRI(JI,JJ,JK)/ZCPD(JI,JJ,JK)
      END DO
    END DO
  END DO
  ! Determine tropopause/inversion  height from minimum temperature
  ITPL(:,:)  = KIB+1
  ZTMIN(:,:) = 400.
  DO JK = IKTB+1,IKTE-1
    DO JJ=KJB,KJE
      DO JI=KIB,KIE
        IF ( PT(JI,JJ,JK) < ZTMIN(JI,JJ) ) THEN
          ZTMIN(JI,JJ) = PT(JI,JJ,JK)
          ITPL(JI,JJ) = JK
        ENDIF
      END DO
    END DO
  END DO
  ! Set the mixing length scale
  ZL(:,:,KKB) = 20.
  DO JK = KKB+KKL,KKE,KKL
    DO JJ=KJB,KJE
      DO JI=KIB,KIE
        ! free troposphere
        ZL(JI,JJ,JK) = ZL0
        ZZZ =  PZZ(JI,JJ,JK) -  PZZ(JI,JJ,KKB)
        JKP = ITPL(JI,JJ)
        ! approximate length for boundary-layer
        IF ( ZL0 > ZZZ ) ZL(JI,JJ,JK) = ZZZ
        ! gradual decrease of length-scale near and above tropopause
        IF ( ZZZ > 0.9*(PZZ(JI,JJ,JKP)-PZZ(JI,JJ,KKB)) ) &
             ZL(JI,JJ,JK) = .6 * ZL(JI,JJ,JK-KKL)
      END DO
    END DO
  END DO
END IF
!-------------------------------------------------------------------------------
!
!
!Ice fraction
ZFRAC(:,:,:) = 0.
IF (OUSERI) THEN
  WHERE(PRC(:,:,:)+PRI(:,:,:) > 1.E-20)
    ZFRAC(:,:,:) = PRI(:,:,:) / (PRC(:,:,:)+PRI(:,:,:))
  ENDWHERE
  CALL COMPUTE_FRAC_ICE(HFRAC_ICE, ZFRAC, PT)
ENDIF
!
DO JK=IKTB,IKTE
  JKP=MAX(MIN(JK+KKL,IKTE),IKTB)
  JKM=MAX(MIN(JK-KKL,IKTE),IKTB)
  DO JJ=KJB,KJE
    DO JI=KIB,KIE
      ! latent heats
      ZTEMP  = PT(JI,JJ,JK)
      ! saturated water vapor mixing ratio over liquid water
      ZPV    = MIN(EXP( XALPW - XBETAW / ZTEMP - XGAMW * LOG( ZTEMP ) ), .99*PPABS(JI,JJ,JK))
      ZQSL   = XRD / XRV * ZPV / ( PPABS(JI,JJ,JK) - ZPV )

      ! saturated water vapor mixing ratio over ice
      ZPIV   = MIN(EXP( XALPI - XBETAI / ZTEMP - XGAMI * LOG( ZTEMP ) ), .99*PPABS(JI,JJ,JK))
      ZQSI   = XRD / XRV * ZPIV / ( PPABS(JI,JJ,JK) - ZPIV )

      ! interpolate between liquid and solid as function of temperature
      ZQSL = (1. - ZFRAC(JI,JJ,JK)) * ZQSL + ZFRAC(JI,JJ,JK) * ZQSI
      ZLVS = (1. - ZFRAC(JI,JJ,JK)) * ZLV(JI,JJ,JK) + &
             & ZFRAC(JI,JJ,JK)      * ZLS(JI,JJ,JK)

      ! coefficients a and b
      ZAH  = ZLVS * ZQSL / ( XRV * ZTEMP**2 ) * (XRV * ZQSL / XRD + 1.)
      ZA   = 1. / ( 1. + ZLVS/ZCPD(JI,JJ,JK) * ZAH )
      ZB   = ZAH * ZA

      ZSBAR = ZA * ( ZRT(JI,JJ,JK) - ZQSL + &
                   & ZAH * ZLVS * (PRC(JI,JJ,JK)+PRI(JI,JJ,JK)) / ZCPD(JI,JJ,JK))

      ! switch to take either present computed value of SIGMAS
      ! or that of Meso-NH turbulence scheme
      IF ( OSIGMAS ) THEN
        IF (PSIGQSAT/=0.) THEN
          ZSIGMA = SQRT((2*PSIGS(JI,JJ,JK))**2 + (PSIGQSAT*ZQSL*ZA)**2)
        ELSE
          ZSIGMA = 2*PSIGS(JI,JJ,JK)
        END IF
      ELSE
        ! parameterize Sigma_s with first_order closure
        DZZ    =  PZZ(JI,JJ,JKP) - PZZ(JI,JJ,JKM)
        ZDRW   =  ZRT(JI,JJ,JKP) - ZRT(JI,JJ,JKM)
        ZDTL   =  ZTLK(JI,JJ,JKP) - ZTLK(JI,JJ,JKM) + XG/ZCPD(JI,JJ,JK) * DZZ
        ZLL = ZL(JI,JJ,JK)
        ! standard deviation due to convection
        ZSIG_CONV =0.
        IF( SIZE(PMFCONV) /= 0) &
             ZSIG_CONV = ZCSIG_CONV * PMFCONV(JI,JJ,JK) / ZA
        ! zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
        ZSIGMA =  SQRT( MAX( 1.E-25, ZCSIGMA * ZCSIGMA * ZLL*ZLL/(DZZ*DZZ)*(&
             ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL + ZB*ZB*ZDTL*ZDTL) + &
             ZSIG_CONV * ZSIG_CONV ) )
      END IF
      ZSIGMA= MAX( 1.E-10, ZSIGMA )
!     ZSIGMA= MAX( 1.E-12, ZSIGMA )

      ! normalized saturation deficit
      ZQ1   = ZSBAR/ZSIGMA

      IF(HCONDENS == 'GAUS')THEN
        ! Gaussian Probability Density Function around ZQ1
        ! Computation of ZG and ZGAM(=erf(ZG))
        ZGCOND = -ZQ1/SQRT(2.)

        !Approximation of erf function for Gaussian distribution
        ZGAUV = 1 - SIGN(1., ZGCOND) * SQRT(1-EXP(-4*ZGCOND**2/XPI))

        !Computation Cloud Fraction
        PCLDFR(JI,JJ,JK) = MAX( 0., MIN(1.,0.5*ZGAUV))

        !Computation of condensate
        ZCOND(JI,JJ,JK) = (EXP(-ZGCOND**2)-ZGCOND*SQRT(XPI)*ZGAUV)*ZSIGMA/SQRT(2.*XPI)
        ZCOND(JI,JJ,JK) = MAX(ZCOND(JI,JJ,JK), 0.)

        PSIGRC(JI,JJ,JK) = PCLDFR(JI,JJ,JK)

        !Computation warm/cold Cloud Fraction and content in high water content part
        IF(PRESENT(PHLC_HCF) .AND. PRESENT(PHLC_HRC))THEN
          IF(1-ZFRAC(JI,JJ,JK) > 1.E-20)THEN
            ZAUTC = (ZSBAR - XCRIAUTC/(PRHODREF(JI,JJ,JK)*(1-ZFRAC(JI,JJ,JK))))/ZSIGMA
            ZGAUTC = -ZAUTC/SQRT(2.)
            !Approximation of erf function for Gaussian distribution
            ZGAUC = 1 - SIGN(1., ZGAUTC) * SQRT(1-EXP(-4*ZGAUTC**2/XPI))
            PHLC_HCF(JI,JJ,JK) = MAX( 0., MIN(1.,0.5*ZGAUC))
            PHLC_HRC(JI,JJ,JK) = (1-ZFRAC(JI,JJ,JK))*(EXP(-ZGAUTC**2)-ZGAUTC*SQRT(XPI)*ZGAUC)*ZSIGMA/SQRT(2.*XPI)
            PHLC_HRC(JI,JJ,JK) = PHLC_HRC(JI,JJ,JK) + XCRIAUTC/PRHODREF(JI,JJ,JK) * PHLC_HCF(JI,JJ,JK)
            PHLC_HRC(JI,JJ,JK) = MAX(PHLC_HRC(JI,JJ,JK), 0.)
          ELSE
            PHLC_HCF(JI,JJ,JK)=0.
            PHLC_HRC(JI,JJ,JK)=0.
          ENDIF
        ENDIF

        IF(PRESENT(PHLI_HCF) .AND. PRESENT(PHLI_HRI))THEN
          IF(ZFRAC(JI,JJ,JK) > 1.E-20)THEN
            ZCRIAUTI(JI,JJ,JK)=MIN(XCRIAUTI,10**(XACRIAUTI*(PT(JI,JJ,JK)-XTT)+XBCRIAUTI))
            ZAUTI = (ZSBAR - ZCRIAUTI(JI,JJ,JK)/ZFRAC(JI,JJ,JK))/ZSIGMA
            ZGAUTI = -ZAUTI/SQRT(2.)
            !Approximation of erf function for Gaussian distribution
            ZGAUI = 1 - SIGN(1., ZGAUTI) * SQRT(1-EXP(-4*ZGAUTI**2/XPI))
            PHLI_HCF(JI,JJ,JK) = MAX( 0., MIN(1.,0.5*ZGAUI))
            PHLI_HRI(JI,JJ,JK) = ZFRAC(JI,JJ,JK)*(EXP(-ZGAUTI**2)-ZGAUTI*SQRT(XPI)*ZGAUI)*ZSIGMA/SQRT(2.*XPI)
            PHLI_HRI(JI,JJ,JK) = PHLI_HRI(JI,JJ,JK) + ZCRIAUTI(JI,JJ,JK)*PHLI_HCF(JI,JJ,JK)
            PHLI_HRI(JI,JJ,JK) = MAX(PHLI_HRI(JI,JJ,JK), 0.)
          ELSE
            PHLI_HCF(JI,JJ,JK)=0.
            PHLI_HRI(JI,JJ,JK)=0.
          ENDIF
        ENDIF

      ELSEIF(HCONDENS == 'CB02')THEN
        !Cloud fraction
        PCLDFR(JI,JJ,JK) = MAX( 0., MIN(1.,0.5+0.36*ATAN(1.55*ZQ1)) )

        !Total condensate
        IF (ZQ1 > 0. .AND. ZQ1 <= 2) THEN
                ZCOND(JI,JJ,JK) = MIN(EXP(-1.)+.66*ZQ1+.086*ZQ1**2, 2.) ! We use the MIN function for continuity
        ELSE IF (ZQ1 > 2.) THEN
                ZCOND(JI,JJ,JK) = ZQ1
        ELSE
                ZCOND(JI,JJ,JK) = EXP( 1.2*ZQ1-1. )
        ENDIF
        ZCOND(JI,JJ,JK) = ZCOND(JI,JJ,JK) * ZSIGMA

        INQ1 = MIN( MAX(-22,FLOOR(MIN(100., MAX(-100., 2*ZQ1))) ), 10)  !inner min/max prevents sigfpe when 2*zq1 does not fit into an int
        ZINC = 2.*ZQ1 - INQ1

        PSIGRC(JI,JJ,JK) =  MIN(1.,(1.-ZINC)*ZSRC_1D(INQ1)+ZINC*ZSRC_1D(INQ1+1))

        IF(PRESENT(PHLC_HCF) .AND. PRESENT(PHLC_HRC))THEN
          PHLC_HCF(JI,JJ,JK)=0.
          PHLC_HRC(JI,JJ,JK)=0.
        ENDIF
        IF(PRESENT(PHLI_HCF) .AND. PRESENT(PHLI_HRI))THEN
          PHLI_HCF(JI,JJ,JK)=0.
          PHLI_HRI(JI,JJ,JK)=0.
        ENDIF
      ENDIF

      IF ( ZCOND(JI,JJ,JK) < 1.E-12 ) THEN
        ZCOND(JI,JJ,JK) = 0.
        PCLDFR(JI,JJ,JK) = 0.
      ENDIF
      IF (PCLDFR(JI,JJ,JK)==0.) THEN
        ZCOND(JI,JJ,JK)=0.
      ENDIF

      ZRCOLD=PRC(JI,JJ,JK)
      ZRIOLD=PRI(JI,JJ,JK)

      PRC(JI,JJ,JK) = (1.-ZFRAC(JI,JJ,JK)) * ZCOND(JI,JJ,JK) ! liquid condensate
      PRI(JI,JJ,JK) = ZFRAC(JI,JJ,JK) * ZCOND(JI,JJ,JK)   ! solid condensate
      
      PT(JI,JJ,JK) = PT(JI,JJ,JK) + ((PRC(JI,JJ,JK)-ZRCOLD)*ZLV(JI,JJ,JK) + &
                                    &(PRI(JI,JJ,JK)-ZRIOLD)*ZLS(JI,JJ,JK)   ) &
                                  & /ZCPD(JI,JJ,JK)
      PRV(JI,JJ,JK) = ZRT(JI,JJ,JK) - PRC(JI,JJ,JK) - PRI(JI,JJ,JK)

! s r_c/ sig_s^2
!    PSIGRC(JI,JJ,JK) = PCLDFR(JI,JJ,JK)  ! use simple Gaussian relation
!
!    multiply PSRCS by the lambda3 coefficient
!
!      PSIGRC(JI,JJ,JK) = 2.*PCLDFR(JI,JJ,JK) * MIN( 3. , MAX(1.,1.-ZQ1) )
! in the 3D case lambda_3 = 1.

      IF(HLAMBDA3=='CB')THEN
          PSIGRC(JI,JJ,JK) = PSIGRC(JI,JJ,JK)* MIN( 3. , MAX(1.,1.-ZQ1) )
      ELSEIF(HLAMBDA3=='NONE') THEN
      ELSE
        call Print_msg( NVERB_FATAL, 'GEN', 'CONDENSATION', 'invalid value for HLAMBDA3: ' // TRIM( HLAMBDA3 ) )
      ENDIF

    END DO
  END DO
END DO
!
END SUBROUTINE CONDENSATION
