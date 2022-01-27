!MNH_LIC Copyright 2002-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
    SUBROUTINE CONDENSATION(D, CST, ICEP, NEB, &
       HFRAC_ICE, HCONDENS, HLAMBDA3,                                                  &
       PPABS, PZZ, PRHODREF, PT, PRV_IN, PRV_OUT, PRC_IN, PRC_OUT, PRI_IN, PRI_OUT,    &
       PRS, PRG, PSIGS, LMFCONV, PMFCONV, PCLDFR, PSIGRC, OUSERI,                               &
       OSIGMAS, OCND2, PSIGQSAT,                                                       &
       PLV, PLS, PCPH,                                                                 &
       PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,                                         &
       PICE_CLD_WGT)
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
!!      2014-11 K.I Ivarsson add possibility to run with OCND2 option
!!      2016   S.Riette Change INQ1
!!      2016-11 S. Riette: use HFRAC_ICE, output adjusted state
!!      2020-12 U. Andrae : Introduce SPP for HARMONIE-AROME
!!     R. El Khatib 24-Aug-2021 Optimizations
!!      2021-01: SPP computations moved in aro_adjust (AROME/HARMONIE)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_CST,            ONLY: CST_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE MODD_NEB,            ONLY: NEB_t
USE MODE_TIWMX,          ONLY : ESATW, ESATI
USE MODE_ICECLOUD,       ONLY : ICECLOUD
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
TYPE(DIMPHYEX_t),             INTENT(IN)    :: D
TYPE(CST_t),                  INTENT(IN)    :: CST
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)    :: ICEP
TYPE(NEB_t),                  INTENT(IN)    :: NEB
CHARACTER(LEN=1),             INTENT(IN)    :: HFRAC_ICE
CHARACTER(LEN=4),             INTENT(IN)    :: HCONDENS
CHARACTER(LEN=*),             INTENT(IN)    :: HLAMBDA3 ! formulation for lambda3 coeff
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PPABS  ! pressure (Pa)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PZZ    ! height of model levels (m)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PRHODREF
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT) :: PT     ! grid scale T  (K)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PRV_IN ! grid scale water vapor mixing ratio (kg/kg) in input
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PRV_OUT! grid scale water vapor mixing ratio (kg/kg) in output
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PRC_IN ! grid scale r_c mixing ratio (kg/kg) in input
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PRC_OUT! grid scale r_c mixing ratio (kg/kg) in output
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PRI_IN ! grid scale r_i (kg/kg) in input
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PRI_OUT! grid scale r_i (kg/kg) in output
LOGICAL, INTENT(IN)                         :: OUSERI ! logical switch to compute both
                                                      ! liquid and solid condensate (OUSERI=.TRUE.)
                                                      ! or only solid condensate (OUSERI=.FALSE.)
LOGICAL, INTENT(IN)                         :: OSIGMAS! use present global Sigma_s values
                                                      ! or that from turbulence scheme
LOGICAL, INTENT(IN)                         :: OCND2  ! logical switch to sparate liquid and ice
                                                      ! more rigid (DEFALT value : .FALSE.)
REAL, DIMENSION(D%NIT,D%NJT),     INTENT(IN)    :: PSIGQSAT ! use an extra "qsat" variance contribution (OSIGMAS case)
                                                        ! multiplied by PSIGQSAT
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PRS    ! grid scale mixing ration of snow (kg/kg)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PRG    ! grid scale mixing ration of graupel (kg/kg)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PSIGS  ! Sigma_s from turbulence scheme
LOGICAL,                                                       INTENT(IN)    ::  LMFCONV ! =SIZE(PMFCONV)!=0
REAL, DIMENSION(MERGE(D%NIT,0,LMFCONV),&
                MERGE(D%NJT,0,LMFCONV),&
                MERGE(D%NKT,0,LMFCONV)),              INTENT(IN)    :: PMFCONV! convective mass flux (kg /s m^2)
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PCLDFR ! cloud fraction
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PSIGRC ! s r_c / sig_s^2

REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PLV    ! Latent heat L_v
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PLS    ! Latent heat L_s
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCPH   ! Specific heat C_ph
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLC_HRC
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLC_HCF ! cloud fraction
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLI_HRI
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLI_HCF
REAL, DIMENSION(D%NIT,D%NJT),   OPTIONAL, INTENT(IN)   :: PICE_CLD_WGT
!
!
!*       0.2   Declarations of local variables :
!
INTEGER  :: JI, JJ, JK, JKP, JKM    ! loop index
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZTLK, ZRT       ! work arrays for T_l and total water mixing ratio
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZL              ! length scale
INTEGER, DIMENSION(D%NIT,D%NJT)  :: ITPL            ! top levels of troposphere 
REAL,    DIMENSION(D%NIT,D%NJT)  :: ZTMIN           ! minimum Temp. related to ITPL
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZLV, ZLS, ZCPD
REAL :: ZGCOND, ZAUTC, ZAUTI, ZGAUV, ZGAUC, ZGAUI, ZGAUTC, ZGAUTI, ZCRIAUTI   ! Used for Gaussian PDF integration
REAL :: ZLVS                                      ! thermodynamics
REAL, DIMENSION(D%NIB:D%NIE) :: ZPV, ZPIV, ZQSL, ZQSI ! thermodynamics
REAL :: ZLL, DZZ, ZZZ                           ! used for length scales 
REAL :: ZAH, ZDRW, ZDTL, ZSIG_CONV                     ! related to computation of Sig_s
REAL, DIMENSION(D%NIB:D%NIE) :: ZA, ZB, ZSBAR, ZSIGMA, ZQ1 ! related to computation of Sig_s
REAL, DIMENSION(D%NIB:D%NIE) :: ZCOND
REAL, DIMENSION(D%NIB:D%NIE) :: ZFRAC           ! Ice fraction
INTEGER  :: INQ1
REAL :: ZINC
! related to OCND2 noise check :
REAL :: ZRSP,  ZRSW, ZRFRAC, ZRSDIF, ZRCOLD
! related to OCND2  ice cloud calulation :
REAL, DIMENSION(D%NIB:D%NIE) :: ESATW_T
REAL :: ZDUM1,ZDUM2,ZDUM3,ZDUM4,ZPRIFACT
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: TCLD
REAL :: ZDZ(D%NIB:D%NIE), &
        ZARDUM(D%NIE-D%NIB+1), ZCLDUM(D%NIE-D%NIB+1)
! end OCND2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER, DIMENSION(D%NIT) :: IERR
!
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
IF (LHOOK) CALL DR_HOOK('CONDENSATION',0,ZHOOK_HANDLE)
!
PCLDFR(:,:,:) = 0. ! Initialize values
PSIGRC(:,:,:) = 0. ! Initialize values
ZPRIFACT = 1.      ! Initialize value
ZCLDUM=-1.         ! Initialize value
IF(OCND2)ZPRIFACT = 0.
!
!
!-------------------------------------------------------------------------------
! store total water mixing ratio
DO JK=D%NKTB,D%NKTE
  DO JJ=D%NJB,D%NJE
    DO JI=D%NIB,D%NIE
      ZRT(JI,JJ,JK)  = PRV_IN(JI,JJ,JK) + PRC_IN(JI,JJ,JK) + PRI_IN(JI,JJ,JK)*ZPRIFACT
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
  DO JK=D%NKTB,D%NKTE
    DO JJ=D%NJB,D%NJE
      DO JI=D%NIB,D%NIE
        ! latent heat of vaporisation/sublimation
        ZLV(JI,JJ,JK) = CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(JI,JJ,JK) - CST%XTT )
        ZLS(JI,JJ,JK) = CST%XLSTT + ( CST%XCPV - CST%XCI ) * ( PT(JI,JJ,JK) - CST%XTT )
      ENDDO
    ENDDO
  ENDDO
ENDIF
IF(PRESENT(PCPH)) THEN
  ZCPD(:,:,:)=PCPH(:,:,:)
ELSE
  DO JK=D%NKTB,D%NKTE
    DO JJ=D%NJB,D%NJE
      DO JI=D%NIB,D%NIE
        ZCPD(JI,JJ,JK) = CST%XCPD + CST%XCPV*PRV_IN(JI,JJ,JK) + CST%XCL*PRC_IN(JI,JJ,JK) + CST%XCI*PRI_IN(JI,JJ,JK) + &
                                    CST%XCI*(PRS(JI,JJ,JK) + PRG(JI,JJ,JK) )
      ENDDO
    ENDDO
  ENDDO
ENDIF
! Preliminary calculations needed for computing the "turbulent part" of Sigma_s
IF ( .NOT. OSIGMAS ) THEN
  DO JK=D%NKTB,D%NKTE
    DO JJ=D%NJB,D%NJE
      DO JI=D%NIB,D%NIE
        ! store temperature at saturation
        ZTLK(JI,JJ,JK) = PT(JI,JJ,JK) - ZLV(JI,JJ,JK)*PRC_IN(JI,JJ,JK)/ZCPD(JI,JJ,JK) &
                                      - ZLS(JI,JJ,JK)*PRI_IN(JI,JJ,JK)/ZCPD(JI,JJ,JK)*ZPRIFACT
      END DO
    END DO
  END DO
  ! Determine tropopause/inversion  height from minimum temperature
#ifdef REPRO48
  ITPL(:,:)  = D%NIB+1
  !I (Sébastien Riette) don't understand why tropopause level is set
  !with the index of the second physical point on the horizontal (i.e. 2+JPHEXT)!!!
  !I assume it is a bug...
#else
  ITPL(:,:)  = D%NKB+D%NKL
#endif
  ZTMIN(:,:) = 400.
  DO JK = D%NKTB+1,D%NKTE-1
    DO JJ=D%NJB,D%NJE
      DO JI=D%NIB,D%NIE
        IF ( PT(JI,JJ,JK) < ZTMIN(JI,JJ) ) THEN
          ZTMIN(JI,JJ) = PT(JI,JJ,JK)
          ITPL(JI,JJ) = JK
        ENDIF
      END DO
    END DO
  END DO
  ! Set the mixing length scale
  ZL(:,:,D%NKB) = 20.
  DO JK = D%NKB+D%NKL,D%NKE,D%NKL
    DO JJ=D%NJB,D%NJE
      DO JI=D%NIB,D%NIE
        ! free troposphere
        ZL(JI,JJ,JK) = ZL0
        ZZZ =  PZZ(JI,JJ,JK) -  PZZ(JI,JJ,D%NKB)
        JKP = ITPL(JI,JJ)
        ! approximate length for boundary-layer
        IF ( ZL0 > ZZZ ) ZL(JI,JJ,JK) = ZZZ
        ! gradual decrease of length-scale near and above tropopause
        IF ( ZZZ > 0.9*(PZZ(JI,JJ,JKP)-PZZ(JI,JJ,D%NKB)) ) &
             ZL(JI,JJ,JK) = .6 * ZL(JI,JJ,JK-D%NKL)
      END DO
    END DO
  END DO
END IF
!-------------------------------------------------------------------------------
!
DO JK=D%NKTB,D%NKTE
  JKP=MAX(MIN(JK+D%NKL,D%NKTE),D%NKTB)
  JKM=MAX(MIN(JK-D%NKL,D%NKTE),D%NKTB)
  DO JJ=D%NJB,D%NJE
    IF (OCND2) THEN
       ZDZ(D%NIB:D%NIE) = PZZ(D%NIB:D%NIE,JJ,JKP) - &
                                            PZZ(D%NIB:D%NIE,JJ,JKP-D%NKL)
       CALL ICECLOUD(D%NIE-D%NIB+1,PPABS(D%NIB,JJ,JK),PZZ(D%NIB,JJ,JK),ZDZ(D%NIB), &
            & PT(D%NIB,JJ,JK),PRV_IN(D%NIB,JJ,JK),1.,-1., &
            & ZCLDUM,1.,TCLD(D%NIB,JJ,JK), &
            & ZARDUM,ZARDUM,ZARDUM,ZARDUM)
       ! latent heats
       ! saturated water vapor mixing ratio over liquid water and ice
       DO JI=D%NIB,D%NIE
         ESATW_T(JI)=ESATW(PT(JI,JJ,JK))
         ZPV(JI)  = MIN(ESATW_T(JI), .99*PPABS(JI,JJ,JK))
         ZPIV(JI) = MIN(ESATI(PT(JI,JJ,JK)), .99*PPABS(JI,JJ,JK))
       END DO
    ELSE
       ! latent heats
       ! saturated water vapor mixing ratio over liquid water and ice
      DO JI=D%NIB,D%NIE
        ZPV(JI)  = MIN(EXP( CST%XALPW - CST%XBETAW / PT(JI,JJ,JK) - CST%XGAMW * LOG( PT(JI,JJ,JK) ) ), .99*PPABS(JI,JJ,JK))
        ZPIV(JI) = MIN(EXP( CST%XALPI - CST%XBETAI / PT(JI,JJ,JK) - CST%XGAMI * LOG( PT(JI,JJ,JK) ) ), .99*PPABS(JI,JJ,JK))
      END DO
    ENDIF
    !Ice fraction
    ZFRAC(:) = 0.
    IF (OUSERI .AND. .NOT.OCND2) THEN
      DO JI=D%NIB,D%NIE
        IF (PRC_IN(JI,JJ,JK)+PRI_IN(JI,JJ,JK) > 1.E-20) THEN
          ZFRAC(JI) = PRI_IN(JI,JJ,JK) / (PRC_IN(JI,JJ,JK)+PRI_IN(JI,JJ,JK))
        ENDIF
      END DO
      CALL COMPUTE_FRAC_ICE(HFRAC_ICE, NEB, ZFRAC(:), PT(:,JJ,JK), IERR) !error code IERR cannot be checked here to not break vectorization
    ENDIF
    DO JI=D%NIB,D%NIE
      ZQSL(JI)   = CST%XRD / CST%XRV * ZPV(JI) / ( PPABS(JI,JJ,JK) - ZPV(JI) )
      ZQSI(JI)   = CST%XRD / CST%XRV * ZPIV(JI) / ( PPABS(JI,JJ,JK) - ZPIV(JI) )

      ! interpolate between liquid and solid as function of temperature
      ZQSL(JI) = (1. - ZFRAC(JI)) * ZQSL(JI) + ZFRAC(JI) * ZQSI(JI)
      ZLVS = (1. - ZFRAC(JI)) * ZLV(JI,JJ,JK) + &
             & ZFRAC(JI)      * ZLS(JI,JJ,JK)

      ! coefficients a and b
      ZAH  = ZLVS * ZQSL(JI) / ( CST%XRV * PT(JI,JJ,JK)**2 ) * (CST%XRV * ZQSL(JI) / CST%XRD + 1.)
      ZA(JI)   = 1. / ( 1. + ZLVS/ZCPD(JI,JJ,JK) * ZAH )
      ZB(JI)   = ZAH * ZA(JI)
      ZSBAR(JI) = ZA(JI) * ( ZRT(JI,JJ,JK) - ZQSL(JI) + &
                   & ZAH * ZLVS * (PRC_IN(JI,JJ,JK)+PRI_IN(JI,JJ,JK)*ZPRIFACT) / ZCPD(JI,JJ,JK))
    END DO
    ! switch to take either present computed value of SIGMAS
    ! or that of Meso-NH turbulence scheme
    IF ( OSIGMAS ) THEN
      DO JI=D%NIB,D%NIE
        IF (PSIGQSAT(JI,JJ)/=0.) THEN
          ZSIGMA(JI) = SQRT((2*PSIGS(JI,JJ,JK))**2 + (PSIGQSAT(JI,JJ)*ZQSL(JI)*ZA(JI))**2)
        ELSE
          ZSIGMA(JI) = 2*PSIGS(JI,JJ,JK)
        END IF
      END DO
    ELSE
      DO JI=D%NIB,D%NIE
        ! parameterize Sigma_s with first_order closure
        DZZ    =  PZZ(JI,JJ,JKP) - PZZ(JI,JJ,JKM)
        ZDRW   =  ZRT(JI,JJ,JKP) - ZRT(JI,JJ,JKM)
        ZDTL   =  ZTLK(JI,JJ,JKP) - ZTLK(JI,JJ,JKM) + CST%XG/ZCPD(JI,JJ,JK) * DZZ
        ZLL = ZL(JI,JJ,JK)
        ! standard deviation due to convection
        ZSIG_CONV =0.
        IF( SIZE(PMFCONV) /= 0) ZSIG_CONV = ZCSIG_CONV * PMFCONV(JI,JJ,JK) / ZA(JI)
        ! zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
        ZSIGMA(JI) =  SQRT( MAX( 1.E-25, ZCSIGMA * ZCSIGMA * ZLL*ZLL/(DZZ*DZZ)*(&
             ZA(JI)*ZA(JI)*ZDRW*ZDRW - 2.*ZA(JI)*ZB(JI)*ZDRW*ZDTL + ZB(JI)*ZB(JI)*ZDTL*ZDTL) + &
             ZSIG_CONV * ZSIG_CONV ) )
      END DO
    END IF
    DO JI=D%NIB,D%NIE
      ZSIGMA(JI)= MAX( 1.E-10, ZSIGMA(JI) )

      ! normalized saturation deficit
      ZQ1(JI)   = ZSBAR(JI)/ZSIGMA(JI)
    END DO
    IF(HCONDENS == 'GAUS') THEN
      DO JI=D%NIB,D%NIE
        ! Gaussian Probability Density Function around ZQ1
        ! Computation of ZG and ZGAM(=erf(ZG))
        ZGCOND = -ZQ1(JI)/SQRT(2.)

        !Approximation of erf function for Gaussian distribution
        ZGAUV = 1 - SIGN(1., ZGCOND) * SQRT(1-EXP(-4*ZGCOND**2/CST%XPI))

        !Computation Cloud Fraction
        PCLDFR(JI,JJ,JK) = MAX( 0., MIN(1.,0.5*ZGAUV))

        !Computation of condensate
        ZCOND(JI) = (EXP(-ZGCOND**2)-ZGCOND*SQRT(CST%XPI)*ZGAUV)*ZSIGMA(JI)/SQRT(2.*CST%XPI)
        ZCOND(JI) = MAX(ZCOND(JI), 0.)

        PSIGRC(JI,JJ,JK) = PCLDFR(JI,JJ,JK)
      END DO
      !Computation warm/cold Cloud Fraction and content in high water content part
      IF(PRESENT(PHLC_HCF) .AND. PRESENT(PHLC_HRC))THEN
        DO JI=D%NIB,D%NIE
          IF(1-ZFRAC(JI) > 1.E-20)THEN
            ZAUTC = (ZSBAR(JI) - ICEP%XCRIAUTC/(PRHODREF(JI,JJ,JK)*(1-ZFRAC(JI))))/ZSIGMA(JI)
            ZGAUTC = -ZAUTC/SQRT(2.)
            !Approximation of erf function for Gaussian distribution
            ZGAUC = 1 - SIGN(1., ZGAUTC) * SQRT(1-EXP(-4*ZGAUTC**2/CST%XPI))
            PHLC_HCF(JI,JJ,JK) = MAX( 0., MIN(1.,0.5*ZGAUC))
            PHLC_HRC(JI,JJ,JK) = (1-ZFRAC(JI))*(EXP(-ZGAUTC**2)-ZGAUTC*SQRT(CST%XPI)*ZGAUC)*ZSIGMA(JI)/SQRT(2.*CST%XPI)
            PHLC_HRC(JI,JJ,JK) = PHLC_HRC(JI,JJ,JK) + ICEP%XCRIAUTC/PRHODREF(JI,JJ,JK) * PHLC_HCF(JI,JJ,JK)
            PHLC_HRC(JI,JJ,JK) = MAX(PHLC_HRC(JI,JJ,JK), 0.)
          ELSE
            PHLC_HCF(JI,JJ,JK)=0.
            PHLC_HRC(JI,JJ,JK)=0.
          ENDIF
        END DO
      ENDIF

      IF(PRESENT(PHLI_HCF) .AND. PRESENT(PHLI_HRI))THEN
        DO JI=D%NIB,D%NIE
          IF(ZFRAC(JI) > 1.E-20)THEN
            ZCRIAUTI=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(PT(JI,JJ,JK)-CST%XTT)+ICEP%XBCRIAUTI))
            ZAUTI = (ZSBAR(JI) - ZCRIAUTI/ZFRAC(JI))/ZSIGMA(JI)
            ZGAUTI = -ZAUTI/SQRT(2.)
            !Approximation of erf function for Gaussian distribution
            ZGAUI = 1 - SIGN(1., ZGAUTI) * SQRT(1-EXP(-4*ZGAUTI**2/CST%XPI))
            PHLI_HCF(JI,JJ,JK) = MAX( 0., MIN(1.,0.5*ZGAUI))
            PHLI_HRI(JI,JJ,JK) = ZFRAC(JI)*(EXP(-ZGAUTI**2)-ZGAUTI*SQRT(CST%XPI)*ZGAUI)*ZSIGMA(JI)/SQRT(2.*CST%XPI)
            PHLI_HRI(JI,JJ,JK) = PHLI_HRI(JI,JJ,JK) + ZCRIAUTI*PHLI_HCF(JI,JJ,JK)
            PHLI_HRI(JI,JJ,JK) = MAX(PHLI_HRI(JI,JJ,JK), 0.)
          ELSE
            PHLI_HCF(JI,JJ,JK)=0.
            PHLI_HRI(JI,JJ,JK)=0.
          ENDIF
        END DO
      ENDIF

    ELSEIF(HCONDENS == 'CB02')THEN
      DO JI=D%NIB,D%NIE
        !Total condensate
        IF (ZQ1(JI) > 0. .AND. ZQ1(JI) <= 2) THEN
          ZCOND(JI) = MIN(EXP(-1.)+.66*ZQ1(JI)+.086*ZQ1(JI)**2, 2.) ! We use the MIN function for continuity
        ELSE IF (ZQ1(JI) > 2.) THEN
          ZCOND(JI) = ZQ1(JI)
        ELSE
          ZCOND(JI) = EXP( 1.2*ZQ1(JI)-1. )
        ENDIF
        ZCOND(JI) = ZCOND(JI) * ZSIGMA(JI)

        !Cloud fraction
        IF (ZCOND(JI) < 1.E-12) THEN
          PCLDFR(JI,JJ,JK) = 0.
        ELSE
          PCLDFR(JI,JJ,JK) = MAX( 0., MIN(1.,0.5+0.36*ATAN(1.55*ZQ1(JI))) )
        ENDIF
        IF (PCLDFR(JI,JJ,JK)==0.) THEN
          ZCOND(JI)=0.
        ENDIF

        INQ1 = MIN( MAX(-22,FLOOR(MIN(100., MAX(-100., 2*ZQ1(JI)))) ), 10)  !inner min/max prevents sigfpe when 2*zq1 does not fit into an int
        ZINC = 2.*ZQ1(JI) - INQ1

        PSIGRC(JI,JJ,JK) =  MIN(1.,(1.-ZINC)*ZSRC_1D(INQ1)+ZINC*ZSRC_1D(INQ1+1))
      END DO
      IF(PRESENT(PHLC_HCF) .AND. PRESENT(PHLC_HRC))THEN
        PHLC_HCF(:,JJ,JK)=0.
        PHLC_HRC(:,JJ,JK)=0.
      ENDIF
      IF(PRESENT(PHLI_HCF) .AND. PRESENT(PHLI_HRI))THEN
        PHLI_HCF(:,JJ,JK)=0.
        PHLI_HRI(:,JJ,JK)=0.
      ENDIF
    END IF !HCONDENS

    IF(.NOT. OCND2) THEN
      DO JI=D%NIB,D%NIE
        PRC_OUT(JI,JJ,JK) = (1.-ZFRAC(JI)) * ZCOND(JI) ! liquid condensate
        PRI_OUT(JI,JJ,JK) = ZFRAC(JI) * ZCOND(JI)   ! solid condensate
        PT(JI,JJ,JK) = PT(JI,JJ,JK) + ((PRC_OUT(JI,JJ,JK)-PRC_IN(JI,JJ,JK))*ZLV(JI,JJ,JK) + &
                                      &(PRI_OUT(JI,JJ,JK)-PRI_IN(JI,JJ,JK))*ZLS(JI,JJ,JK)   ) &
                                    & /ZCPD(JI,JJ,JK)
        PRV_OUT(JI,JJ,JK) = ZRT(JI,JJ,JK) - PRC_OUT(JI,JJ,JK) - PRI_OUT(JI,JJ,JK)*ZPRIFACT
      END DO
    ELSE
      DO JI=D%NIB,D%NIE
        PRC_OUT(JI,JJ,JK) = (1.-ZFRAC(JI)) * ZCOND(JI) ! liquid condensate
        !
!       This check is mainly for noise reduction :
!       -------------------------
        IF(ABS(PRC_IN(JI,JJ,JK)-PRC_OUT(JI,JJ,JK))>1.0E-12 .AND. ESATW_T(JI) < PPABS(JI,JJ,JK)*0.5)THEN
           ZRCOLD = PRC_OUT(JI,JJ,JK)
           ZRFRAC = PRV_IN(JI,JJ,JK) - ZCOND(JI) + PRC_OUT(JI,JJ,JK)
           IF( PRV_IN(JI,JJ,JK) < ZRSW )THEN ! sub - saturation over water:
              ! Avoid drying of cloudwater leading to supersaturation with
              ! respect to water
              ZRSDIF= MIN(0.,ZRSP-ZRFRAC)
           ELSE  ! super - saturation over water:
              ! Avoid depostition of water leading to sub-saturation with
              ! respect to water
              !            ZRSDIF= MAX(0.,ZRSP-ZRFRAC)
              ZRSDIF= MAX(0.,ZRSP*PCLDFR(JI,JJ,JK) - ZRFRAC) 
           ENDIF
           PRC_OUT(JI,JJ,JK) = ZCOND(JI)  - ZRSDIF
        ELSE
          ZRCOLD = PRC_IN(JI,JJ,JK)
        ENDIF
 !      end check 

 !      compute separate ice cloud:
        ZDUM1 = MIN(1.0,20.* PRC_OUT(JI,JJ,JK)*SQRT(ZDZ(JI))/ZQSL(JI)) ! clould liquid water 
                                                       ! factor 

        ZDUM3 = MAX(0.,TCLD(JI,JJ,JK)-PCLDFR(JI,JJ,JK)) ! pure ice cloud part

        IF (JK==D%NKTB) THEN
          ZDUM4 = PRI_IN(JI,JJ,JK)
        ELSE
          ZDUM4 = PRI_IN(JI,JJ,JK) + PRS(JI,JJ,JK)*0.5 + PRG(JI,JJ,JK)*0.25
        ENDIF

        ZDUM4 = MAX(0.,MIN(1.,PICE_CLD_WGT(JI,JJ)*ZDUM4*SQRT(ZDZ(JI))/ZQSI(JI))) ! clould ice+solid 
                                                           ! precip. water factor 

        ZDUM2 = (0.8*PCLDFR(JI,JJ,JK)+0.2)*MIN(1.,ZDUM1 + ZDUM4*PCLDFR(JI,JJ,JK))
        ! water cloud, use 'statistical' cloud, but reduce it in case of low liquid content

        PCLDFR(JI,JJ,JK) = MIN(1., ZDUM2 + (0.9*ZDUM3+0.1)*ZDUM4) ! Rad cloud
             ! Reduce ice cloud part in case of low ice water content
        PRI_OUT(JI,JJ,JK) = PRI_IN(JI,JJ,JK)
        PT(JI,JJ,JK) = PT(JI,JJ,JK) + ((PRC_OUT(JI,JJ,JK)-ZRCOLD)*ZLV(JI,JJ,JK) + &
                                      &(PRI_OUT(JI,JJ,JK)-PRI_IN(JI,JJ,JK))*ZLS(JI,JJ,JK)   ) &
                                    & /ZCPD(JI,JJ,JK)
        PRV_OUT(JI,JJ,JK) = ZRT(JI,JJ,JK) - PRC_OUT(JI,JJ,JK) - PRI_OUT(JI,JJ,JK)*ZPRIFACT
      END DO
    END IF ! End OCND2
    IF(HLAMBDA3=='CB')THEN
      DO JI=D%NIB,D%NIE
        ! s r_c/ sig_s^2
        !    PSIGRC(JI,JJ,JK) = PCLDFR(JI,JJ,JK)  ! use simple Gaussian relation
        !
        !    multiply PSRCS by the lambda3 coefficient
        !
        !      PSIGRC(JI,JJ,JK) = 2.*PCLDFR(JI,JJ,JK) * MIN( 3. , MAX(1.,1.-ZQ1(JI)) )
        ! in the 3D case lambda_3 = 1.

        PSIGRC(JI,JJ,JK) = PSIGRC(JI,JJ,JK)* MIN( 3. , MAX(1.,1.-ZQ1(JI)) )
      END DO
    END IF
  END DO
END DO
!
IF (LHOOK) CALL DR_HOOK('CONDENSATION',1,ZHOOK_HANDLE)
!
CONTAINS
INCLUDE "compute_frac_ice.func.h"
!
END SUBROUTINE CONDENSATION
