!MNH_LIC Copyright 2002-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
    SUBROUTINE CONDENSATION(D, CST, ICEP, NEBN, TURBN, &
                           &HFRAC_ICE, HCONDENS, HLAMBDA3,                                                  &
                           &PPABS, PZZ, PRHODREF, PT, PRV_IN, PRV_OUT, PRC_IN, PRC_OUT, PRI_IN, PRI_OUT,    &
                           &PRR, PRS, PRG, PSIGS, LMFCONV, PMFCONV, PCLDFR, PSIGRC, OUSERI,                 &
                           &OSIGMAS, OCND2,                                                                 &
                           &PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR, PSIGQSAT,                                 &
                           &PLV, PLS, PCPH,                                                                 &
                           &PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,                                         &
                           &PICE_CLD_WGT)
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
!!      2018-02 K.I Ivarsson: Some modificatons of OCND2 option, mainly for optimation - new outputs
!!      2019-06 W.C. de Rooy: Mods for new set up statistical cloud scheme
!!      2019-07 K.I.Ivarsson: Switch for height dependent VQSIGSAT: LHGT_QS
!!      2020-12 U. Andrae : Introduce SPP for HARMONIE-AROME
!!     R. El Khatib 24-Aug-2021 Optimizations
!!      2021-01: SPP computations moved in aro_adjust (AROME/HARMONIE)
!!      Jan 2025 A. Marcel: use ERF instead of approximate formulation
!!      Jan 2025: protection against too small values
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_CST,            ONLY: CST_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_NEB_n,          ONLY: NEB_t
USE MODD_TURB_n,     ONLY: TURB_t
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
TYPE(NEB_t),                  INTENT(IN)    :: NEBN
TYPE(TURB_t),                 INTENT(IN)    :: TURBN
CHARACTER(LEN=1),             INTENT(IN)    :: HFRAC_ICE
CHARACTER(LEN=4),             INTENT(IN)    :: HCONDENS
CHARACTER(LEN=4),             INTENT(IN)    :: HLAMBDA3 ! formulation for lambda3 coeff
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PPABS  ! pressure (Pa)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PZZ    ! height of model levels (m)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRHODREF
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PT     ! grid scale T  (K)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRV_IN ! grid scale water vapor mixing ratio (kg/kg) in input
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PRV_OUT! grid scale water vapor mixing ratio (kg/kg) in output
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRC_IN ! grid scale r_c mixing ratio (kg/kg) in input
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PRC_OUT! grid scale r_c mixing ratio (kg/kg) in output
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRI_IN ! grid scale r_i (kg/kg) in input
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PRI_OUT! grid scale r_i (kg/kg) in output
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRR    ! grid scale mixing ration of rain (kg/kg)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRS    ! grid scale mixing ration of snow (kg/kg)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRG    ! grid scale mixing ration of graupel (kg/kg)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PSIGS  ! Sigma_s from turbulence scheme
LOGICAL,                                                       INTENT(IN)    ::  LMFCONV ! =SIZE(PMFCONV)!=0
REAL, DIMENSION(MERGE(D%NIJT,0,LMFCONV),&
                MERGE(D%NKT,0,LMFCONV)),              INTENT(IN)    :: PMFCONV! convective mass flux (kg /s m^2)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PCLDFR ! cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PSIGRC ! s r_c / sig_s^2

LOGICAL, INTENT(IN)                         :: OUSERI ! logical switch to compute both
                                                      ! liquid and solid condensate (OUSERI=.TRUE.)
                                                      ! or only solid condensate (OUSERI=.FALSE.)
LOGICAL, INTENT(IN)                         :: OSIGMAS! use present global Sigma_s values
                                                      ! or that from turbulence scheme
LOGICAL, INTENT(IN)                         :: OCND2  ! logical switch to sparate liquid and ice
                                                      ! more rigid (DEFALT value : .FALSE.)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PICLDFR  ! ice cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PWCLDFR  ! water or mixed-phase cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PSSIO    ! Super-saturation with respect to ice in the  
                                                              ! supersaturated fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PSSIU    ! Sub-saturation with respect to ice in the  
                                                              ! subsaturated fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PIFR     ! Ratio cloud ice moist part
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PSIGQSAT ! use an extra "qsat" variance contribution (OSIGMAS case)
                                                              ! multiplied by PSIGQSAT

REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PLV    ! Latent heat L_v
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PLS    ! Latent heat L_s
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCPH   ! Specific heat C_ph
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLC_HRC
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLC_HCF ! cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLI_HRI
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLI_HCF
REAL, DIMENSION(D%NIJT),       OPTIONAL, INTENT(IN)    :: PICE_CLD_WGT
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JIJ, JK, JKP, JKM                    ! loop index
INTEGER :: IKTB, IKTE, IKB, IKE, IKL, IIJB, IIJE
REAL, DIMENSION(D%NIJT,D%NKT) :: ZTLK, ZRT     ! work arrays for T_l and total water mixing ratio
REAL, DIMENSION(D%NIJT,D%NKT) :: ZL            ! length scale
INTEGER, DIMENSION(D%NIJT)  :: ITPL            ! top levels of troposphere
REAL,    DIMENSION(D%NIJT)  :: ZTMIN           ! minimum Temp. related to ITPL
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZLV, ZLS, ZCPD
REAL, DIMENSION(D%NIJT,D%NKT) :: ZGCOND, ZAUTC, ZAUTI, ZGAUV, ZGAUC, ZGAUI, ZGAUTC, ZGAUTI, ZCRIAUTI   ! Used for Gaussian PDF integration
REAL, DIMENSION(D%NIJT,D%NKT) :: ZLVS                                      ! thermodynamics
REAL, DIMENSION(D%NIJT,D%NKT) :: ZPV, ZPIV, ZQSL, ZQSI ! thermodynamics
REAL, DIMENSION(D%NIJT,D%NKT) :: ZLL, DZZ, ZZZ                           ! used for length scales
REAL, DIMENSION(D%NIJT,D%NKT) :: ZAH, ZDRW, ZDTL, ZSIG_CONV                     ! related to computation of Sig_s
REAL, DIMENSION(D%NIJT,D%NKT) :: ZA, ZB, ZSBAR, ZSIGMA, ZQ1 ! related to computation of Sig_s
REAL, DIMENSION(D%NIJT,D%NKT) :: ZCOND
REAL, DIMENSION(D%NIJT,D%NKT) :: ZFRAC           ! Ice fraction
INTEGER, DIMENSION(D%NIJT,D%NKT)  :: INQ1
REAL, DIMENSION(D%NIJT,D%NKT) :: ZINC
! related to OCND2 noise check :
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRSP,  ZRSW, ZRFRAC, ZRSDIF, ZRCOLD
! related to OCND2  ice cloud calulation :
REAL, DIMENSION(D%NIJT,D%NKT) :: ESATW_T
REAL, DIMENSION(D%NIJT,D%NKT) :: ZDUM1,ZDUM2,ZDUM3,ZDUM4,ZPRIFACT,ZLWINC
REAL, DIMENSION(D%NIJT,D%NKT) :: ZDZ, ZARDUM, ZARDUM2, ZCLDINI
REAL, DIMENSION(D%NIJT,D%NKT) :: ZSIGQSAT, ZICE_CLD_WGT
! end OCND2

! LHGT_QS:
REAL, DIMENSION(D%NIJT,D%NKT) :: ZDZFACT,ZDZREF
! LHGT_QS END

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: IERR
INTEGER, DIMENSION(D%NKT) :: JKPK, JKMK
!
!
!*       0.3  Definition of constants :
!
!-------------------------------------------------------------------------------
!
REAL,PARAMETER :: ZL0     = 600.        ! tropospheric length scale
REAL,PARAMETER :: ZCSIGMA = 0.2         ! constant in sigma_s parameterization
REAL,PARAMETER :: ZCSIG_CONV = 0.30E-2  ! scaling factor for ZSIG_CONV(JIJ,JK) as function of mass flux
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
IKTB=D%NKTB
IKTE=D%NKTE
IKB=D%NKB
IKE=D%NKE
IKL=D%NKL
IIJB=D%NIJB
IIJE=D%NIJE
!

DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE
    ZSIGQSAT(JIJ,JK) = PSIGQSAT(JIJ)
  ENDDO
ENDDO

IF(PRESENT(PICE_CLD_WGT)) THEN

  DO JK=IKTB,IKTE
    DO JIJ=IIJB,IIJE
      ZICE_CLD_WGT(JIJ,JK) = PICE_CLD_WGT(JIJ)
    ENDDO
  ENDDO

END IF


PCLDFR(:,:) = 0. ! Initialize values
PSIGRC(:,:) = 0. ! Initialize values
PRV_OUT(:,:)= 0. ! Initialize values
PRC_OUT(:,:)= 0. ! Initialize values
PRI_OUT(:,:)= 0. ! Initialize values
ZPRIFACT(:,:) = 1.    ! Initialize value
ZARDUM2(:,:) = 0.  ! Initialize values
ZARDUM(:,:) = 0.  ! Initialize values
ZCLDINI(:,:) = -1. ! Dummy Initialized cloud input to icecloud routine
PIFR(:,:) = 10. ! ratio of cloud ice water mixing ratio wet to dry
           ! part of a gridbox
ZDZREF(:,:) = ICEP%XFRMIN(25) ! Thickness for unchanged vqsigsat (only used for LHGT_QS)
!
IF(OCND2) ZRSW(:,:) = -9999 ! This variable was never initialized !
!
IF(OCND2)ZPRIFACT(:,:) = 0.

!
!
!-------------------------------------------------------------------------------
! store total water mixing ratio


DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE
    ZRT(JIJ,JK)  = PRV_IN(JIJ,JK) + PRC_IN(JIJ,JK) + PRI_IN(JIJ,JK)*ZPRIFACT(JIJ,JK)
  END DO
END DO

!-------------------------------------------------------------------------------
! Preliminary calculations
! latent heat of vaporisation/sublimation
IF(PRESENT(PLV) .AND. PRESENT(PLS)) THEN

  ZLV(:,:)=PLV(:,:)
  ZLS(:,:)=PLS(:,:)

ELSE


  DO JK=IKTB,IKTE
    DO JIJ=IIJB,IIJE
      ! latent heat of vaporisation/sublimation
      ZLV(JIJ,JK) = CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(JIJ,JK) - CST%XTT )
      ZLS(JIJ,JK) = CST%XLSTT + ( CST%XCPV - CST%XCI ) * ( PT(JIJ,JK) - CST%XTT )
    ENDDO
  ENDDO

ENDIF
IF(PRESENT(PCPH)) THEN

  ZCPD(:,:)=PCPH(:,:)

ELSE


  DO JK=IKTB,IKTE
    DO JIJ=IIJB,IIJE
      ZCPD(JIJ,JK) = CST%XCPD + CST%XCPV*PRV_IN(JIJ,JK) + CST%XCL*PRC_IN(JIJ,JK) + CST%XCI*PRI_IN(JIJ,JK) + &
                                  CST%XCL*PRR(JIJ,JK) +  &
                                  CST%XCI*(PRS(JIJ,JK) + PRG(JIJ,JK) )
    ENDDO
  ENDDO

ENDIF
! Preliminary calculations needed for computing the "turbulent part" of Sigma_s
IF ( .NOT. OSIGMAS ) THEN


  DO JK=IKTB,IKTE
    DO JIJ=IIJB,IIJE
      ! store temperature at saturation
      ZTLK(JIJ,JK) = PT(JIJ,JK) - ZLV(JIJ,JK)*PRC_IN(JIJ,JK)/ZCPD(JIJ,JK) &
                                    - ZLS(JIJ,JK)*PRI_IN(JIJ,JK)/ZCPD(JIJ,JK)*ZPRIFACT(JIJ,JK)
    END DO
  END DO

  ! Determine tropopause/inversion  height from minimum temperature

  ITPL(:)  = IKB+IKL
  ZTMIN(:) = 400.



  DO JK = IKTB+1,IKTE-1
    DO JIJ=IIJB,IIJE
      IF ( PT(JIJ,JK) < ZTMIN(JIJ) ) THEN
        ZTMIN(JIJ) = PT(JIJ,JK)
        ITPL(JIJ) = JK
      ENDIF
    END DO
  END DO

  ! Set the mixing length scale

  ZL(:,IKB) = 20.

  DO JK = IKB+IKL,IKE,IKL
    DO JIJ=IIJB,IIJE
      ! free troposphere
      ZL(JIJ,JK) = ZL0
      ZZZ(JIJ,JK) =  PZZ(JIJ,JK) -  PZZ(JIJ,IKB)
      JKP = ITPL(JIJ)
      ! approximate length for boundary-layer
      IF ( ZL0 > ZZZ(JIJ,JK) ) ZL(JIJ,JK) = ZZZ(JIJ,JK)
      ! gradual decrease of length-scale near and above tropopause
      IF ( ZZZ(JIJ,JK) > 0.9*(PZZ(JIJ,JKP)-PZZ(JIJ,IKB)) ) &
           ZL(JIJ,JK) = .6 * ZL(JIJ,JK-IKL)
    END DO
  END DO

END IF
!-------------------------------------------------------------------------------
!


DO JK = IKTB, IKTE
 JKPK(JK) = MAX(MIN(JK+IKL, IKTE),IKTB)
 JKMK(JK) = MAX(MIN(JK-IKL, IKTE),IKTB)
END DO

!
!
IF (OCND2) THEN


  DO JK=IKTB,IKTE
    DO JIJ = IIJB, IIJE
      ZDZ(JIJ,JK) = PZZ(JIJ,JKPK(JK)) - PZZ(JIJ,JKPK(JK)-IKL)
    END DO
  END DO

  CALL ICECLOUD(D, CST, ICEP, PPABS(:,:),PZZ(:,:),ZDZ(:,:), &
       & PT(:,:),PRV_IN(:,:),1.,-1., &
       & ZCLDINI(:,:),PIFR(IIJB,1),PICLDFR(:,:), &
       & PSSIO(:,:),PSSIU(:,:),ZARDUM2(:,:),ZARDUM(:,:))
  ! latent heats
  ! saturated water vapor mixing ratio over liquid water and ice


  DO JK=IKTB,IKTE
    DO JIJ = IIJB, IIJE
      ESATW_T(JIJ,JK)=ESATW(ICEP%TIWMX, PT(JIJ,JK))
      ZPV(JIJ,JK)  = MIN(ESATW_T(JIJ,JK), .99*PPABS(JIJ,JK))
      ZPIV(JIJ,JK) = MIN(ESATI(ICEP%TIWMX, PT(JIJ,JK)), .99*PPABS(JIJ,JK))
    END DO
  END DO

END IF



DO JK=IKTB,IKTE
  DO JIJ = IIJB, IIJE
  IF (.NOT. OCND2) THEN
     ! latent heats
     ! saturated water vapor mixing ratio over liquid water and ice
      ZPV(JIJ,JK)  = MIN(EXP( CST%XALPW - CST%XBETAW / PT(JIJ,JK) - CST%XGAMW * LOG( PT(JIJ,JK) ) ), .99*PPABS(JIJ,JK))
      ZPIV(JIJ,JK) = MIN(EXP( CST%XALPI - CST%XBETAI / PT(JIJ,JK) - CST%XGAMI * LOG( PT(JIJ,JK) ) ), .99*PPABS(JIJ,JK))
  ENDIF
  !Ice fraction
  ZFRAC(JIJ,JK) = 0.
  IF (OUSERI .AND. .NOT.OCND2) THEN
      IF (PRC_IN(JIJ,JK)+PRI_IN(JIJ,JK) > 1.E-20) THEN
        ZFRAC(JIJ,JK) = PRI_IN(JIJ,JK) / (PRC_IN(JIJ,JK)+PRI_IN(JIJ,JK))
      ENDIF
    CALL COMPUTE_FRAC_ICE(CST, HFRAC_ICE, NEBN, ZFRAC(JIJ,JK), PT(JIJ,JK), IERR) !error code IERR cannot be checked here to not break vectorization
  ENDIF
    ZQSL(JIJ,JK)   = CST%XRD / CST%XRV * ZPV(JIJ,JK) / ( PPABS(JIJ,JK) - ZPV(JIJ,JK) )
    ZQSI(JIJ,JK)   = CST%XRD / CST%XRV * ZPIV(JIJ,JK) / ( PPABS(JIJ,JK) - ZPIV(JIJ,JK) )

    ! interpolate between liquid and solid as function of temperature
    ZQSL(JIJ,JK) = (1. - ZFRAC(JIJ,JK)) * ZQSL(JIJ,JK) + ZFRAC(JIJ,JK) * ZQSI(JIJ,JK)
    ZLVS(JIJ,JK) = (1. - ZFRAC(JIJ,JK)) * ZLV(JIJ,JK) + &
           & ZFRAC(JIJ,JK)      * ZLS(JIJ,JK)

    ! coefficients a and b
    ZAH(JIJ,JK)  = ZLVS(JIJ,JK) * ZQSL(JIJ,JK) / ( CST%XRV * PT(JIJ,JK)**2 ) * (CST%XRV * ZQSL(JIJ,JK) / CST%XRD + 1.)
    ZA(JIJ,JK)   = 1. / ( 1. + ZLVS(JIJ,JK)/ZCPD(JIJ,JK) * ZAH(JIJ,JK) )
    ZB(JIJ,JK)   = ZAH(JIJ,JK) * ZA(JIJ,JK)
    ZSBAR(JIJ,JK) = ZA(JIJ,JK) * ( ZRT(JIJ,JK) - ZQSL(JIJ,JK) + &
                 & ZAH(JIJ,JK) * ZLVS(JIJ,JK) * (PRC_IN(JIJ,JK)+PRI_IN(JIJ,JK)*ZPRIFACT(JIJ,JK)) / ZCPD(JIJ,JK))
  ! switch to take either present computed value of SIGMAS
  ! or that of Meso-NH turbulence scheme
  IF ( OSIGMAS ) THEN
      IF (ZSIGQSAT(JIJ,JK)/=0.) THEN
        ZDZFACT(JIJ,JK) = 1.
        IF(NEBN%LHGT_QS .AND. JK+1 <= IKTE)THEN
           ZDZFACT(JIJ,JK)= MAX(ICEP%XFRMIN(23),MIN(ICEP%XFRMIN(24),(PZZ(JIJ,JK) - PZZ(JIJ,JK+1))/ZDZREF(JIJ,JK)))
        ELSEIF(NEBN%LHGT_QS)THEN
           ZDZFACT(JIJ,JK)= MAX(ICEP%XFRMIN(23),MIN(ICEP%XFRMIN(24),((PZZ(JIJ,JK-1) - PZZ(JIJ,JK)))*0.8/ZDZREF(JIJ,JK)))
        ENDIF
        IF (NEBN%LSTATNW) THEN
          ZSIGMA(JIJ,JK) = SQRT((PSIGS(JIJ,JK))**2 + (ZSIGQSAT(JIJ,JK)*ZDZFACT(JIJ,JK)*ZQSL(JIJ,JK)*ZA(JIJ,JK))**2)
        ELSE
          ZSIGMA(JIJ,JK) = SQRT((2*PSIGS(JIJ,JK))**2 + (ZSIGQSAT(JIJ,JK)*ZQSL(JIJ,JK)*ZA(JIJ,JK))**2)
        ENDIF
      ELSE
        IF (NEBN%LSTATNW) THEN
          ZSIGMA(JIJ,JK) = PSIGS(JIJ,JK)
        ELSE
          ZSIGMA(JIJ,JK) = 2*PSIGS(JIJ,JK)
        ENDIF
      END IF
      IF (NEBN%LCONDBORN) THEN
        ZSIGMA(JIJ,JK)=MIN(ZSIGMA(JIJ,JK),ZRT(JIJ,JK)/2.) 
      ENDIF
  ELSE
      ! parameterize Sigma_s with first_order closure
      DZZ(JIJ,JK)    =  PZZ(JIJ,JKPK(JK)) - PZZ(JIJ,JKMK(JK))
      ZDRW(JIJ,JK)   =  ZRT(JIJ,JKPK(JK)) - ZRT(JIJ,JKMK(JK))
      ZDTL(JIJ,JK)   =  ZTLK(JIJ,JKPK(JK)) - ZTLK(JIJ,JKMK(JK)) + CST%XG/ZCPD(JIJ,JK) * DZZ(JIJ,JK)
      ZLL(JIJ,JK) = ZL(JIJ,JK)
      ! standard deviation due to convection
      ZSIG_CONV(JIJ,JK) =0.
      IF(LMFCONV) ZSIG_CONV(JIJ,JK) = ZCSIG_CONV * PMFCONV(JIJ,JK) / ZA(JIJ,JK)
      ! zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
      ZSIGMA(JIJ,JK) =  SQRT( MAX( 1.E-25, ZCSIGMA * ZCSIGMA * ZLL(JIJ,JK)*ZLL(JIJ,JK)/(DZZ(JIJ,JK)*DZZ(JIJ,JK))*(&
           ZA(JIJ,JK)*ZA(JIJ,JK)*ZDRW(JIJ,JK)*ZDRW(JIJ,JK) - 2.*ZA(JIJ,JK)*ZB(JIJ,JK)*ZDRW(JIJ,JK)*ZDTL(JIJ,JK) + ZB(JIJ,JK)*ZB(JIJ,JK)*ZDTL(JIJ,JK)*ZDTL(JIJ,JK)) + &
           ZSIG_CONV(JIJ,JK) * ZSIG_CONV(JIJ,JK) ) )
  END IF
    ZSIGMA(JIJ,JK)= MAX( 1.E-10, ZSIGMA(JIJ,JK) )

    ! normalized saturation deficit
    ZQ1(JIJ,JK)   = ZSBAR(JIJ,JK)/ZSIGMA(JIJ,JK)
  IF(HCONDENS == 'GAUS') THEN
      ! Gaussian Probability Density Function around ZQ1
      ! Computation of ZG and ZGAM(=erf(ZG))
      ZGCOND(JIJ,JK) = -ZQ1(JIJ,JK)/SQRT(2.)
      ZGAUV(JIJ,JK) = 1 + ERF(-ZGCOND(JIJ,JK))

      !Computation Cloud Fraction
      PCLDFR(JIJ,JK) = MAX( 0., MIN(1.,0.5*ZGAUV(JIJ,JK)))

      !Computation of condensate
      ZCOND(JIJ,JK) = (EXP(-ZGCOND(JIJ,JK)**2)-ZGCOND(JIJ,JK)*SQRT(CST%XPI)*ZGAUV(JIJ,JK))*ZSIGMA(JIJ,JK)/SQRT(2.*CST%XPI)
      ZCOND(JIJ,JK) = MAX(ZCOND(JIJ,JK), 0.)
      IF (NEBN%LCONDBORN) THEN
        ZCOND(JIJ,JK) = MIN(ZCOND(JIJ,JK),ZRT(JIJ,JK)-1.E-7)
      ENDIF

      IF (ZCOND(JIJ,JK) < 1.E-12 .OR. PCLDFR(JIJ,JK) == 0.) THEN
        PCLDFR(JIJ,JK)=0.
        ZCOND(JIJ,JK)=0.
      ENDIF

      PSIGRC(JIJ,JK) = PCLDFR(JIJ,JK)
    !Computation warm/cold Cloud Fraction and content in high water content part
    IF(PRESENT(PHLC_HCF) .AND. PRESENT(PHLC_HRC))THEN
        IF(1-ZFRAC(JIJ,JK) > 1.E-20)THEN
          ZAUTC(JIJ,JK) = (ZSBAR(JIJ,JK) - ICEP%XCRIAUTC/(PRHODREF(JIJ,JK)*(1-ZFRAC(JIJ,JK))))/ZSIGMA(JIJ,JK)
          ZGAUTC(JIJ,JK) = -ZAUTC(JIJ,JK)/SQRT(2.)
          ZGAUC(JIJ,JK) = 1 + ERF(-ZGAUTC(JIJ,JK))
          PHLC_HCF(JIJ,JK) = MAX( 0., MIN(1.,0.5*ZGAUC(JIJ,JK)))
          PHLC_HRC(JIJ,JK) = (1-ZFRAC(JIJ,JK))*(EXP(-ZGAUTC(JIJ,JK)**2)-ZGAUTC(JIJ,JK)*SQRT(CST%XPI)*ZGAUC(JIJ,JK))*ZSIGMA(JIJ,JK)/SQRT(2.*CST%XPI)
          PHLC_HRC(JIJ,JK) = PHLC_HRC(JIJ,JK) + ICEP%XCRIAUTC/PRHODREF(JIJ,JK) * PHLC_HCF(JIJ,JK)
          PHLC_HRC(JIJ,JK) = MAX(PHLC_HRC(JIJ,JK), 0.)
          IF(PHLC_HRC(JIJ,JK) < 1.E-12 .OR. PHLC_HCF(JIJ,JK) < 1.E-6) THEN
            PHLC_HRC(JIJ,JK)=0.
            PHLC_HCF(JIJ,JK)=0.
          ENDIF
        ELSE
          PHLC_HCF(JIJ,JK)=0.
          PHLC_HRC(JIJ,JK)=0.
        ENDIF
    ENDIF

    IF(PRESENT(PHLI_HCF) .AND. PRESENT(PHLI_HRI))THEN
        IF(ZFRAC(JIJ,JK) > 1.E-20)THEN
          ZCRIAUTI(JIJ,JK)=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(PT(JIJ,JK)-CST%XTT)+ICEP%XBCRIAUTI))
          ZAUTI(JIJ,JK) = (ZSBAR(JIJ,JK) - ZCRIAUTI(JIJ,JK)/ZFRAC(JIJ,JK))/ZSIGMA(JIJ,JK)
          ZGAUTI(JIJ,JK) = -ZAUTI(JIJ,JK)/SQRT(2.)
          ZGAUI(JIJ,JK) = 1 + ERF(-ZGAUTI(JIJ,JK))
          PHLI_HCF(JIJ,JK) = MAX( 0., MIN(1.,0.5*ZGAUI(JIJ,JK)))
          PHLI_HRI(JIJ,JK) = ZFRAC(JIJ,JK)*(EXP(-ZGAUTI(JIJ,JK)**2)-ZGAUTI(JIJ,JK)*SQRT(CST%XPI)*ZGAUI(JIJ,JK))*ZSIGMA(JIJ,JK)/SQRT(2.*CST%XPI)
          PHLI_HRI(JIJ,JK) = PHLI_HRI(JIJ,JK) + ZCRIAUTI(JIJ,JK)*PHLI_HCF(JIJ,JK)
          PHLI_HRI(JIJ,JK) = MAX(PHLI_HRI(JIJ,JK), 0.)
          IF(PHLI_HRI(JIJ,JK) < 1.E-12 .OR. PHLI_HCF(JIJ,JK) < 1.E-6) THEN
            PHLI_HRI(JIJ,JK)=0.
            PHLI_HCF(JIJ,JK)=0.
          ENDIF
        ELSE
          PHLI_HCF(JIJ,JK)=0.
          PHLI_HRI(JIJ,JK)=0.
        ENDIF
    ENDIF

  ELSEIF(HCONDENS == 'CB02')THEN
      !Total condensate
      IF (ZQ1(JIJ,JK) > 0. .AND. ZQ1(JIJ,JK) <= 2) THEN
        ZCOND(JIJ,JK) = MIN(EXP(-1.)+.66*ZQ1(JIJ,JK)+.086*ZQ1(JIJ,JK)**2, 2.) ! We use the MIN function for continuity
      ELSE IF (ZQ1(JIJ,JK) > 2.) THEN
        ZCOND(JIJ,JK) = ZQ1(JIJ,JK)
      ELSE
        ZCOND(JIJ,JK) = EXP( 1.2*ZQ1(JIJ,JK)-1. )
      ENDIF
      ZCOND(JIJ,JK) = ZCOND(JIJ,JK) * ZSIGMA(JIJ,JK)
      IF (NEBN%LCONDBORN) THEN
        ZCOND(JIJ,JK) = MIN(ZCOND(JIJ,JK),ZRT(JIJ,JK)-1.E-7)
      ENDIF

      !Cloud fraction
      IF (ZCOND(JIJ,JK) < 1.E-12) THEN
        PCLDFR(JIJ,JK) = 0.
      ELSE
        PCLDFR(JIJ,JK) = MAX( 0., MIN(1.,0.5+0.36*ATAN(1.55*ZQ1(JIJ,JK))) )
      ENDIF
      IF (PCLDFR(JIJ,JK)==0.) THEN
        ZCOND(JIJ,JK)=0.
      ENDIF

      INQ1(JIJ,JK) = MIN( MAX(-22,FLOOR(MIN(100., MAX(-100., 2*ZQ1(JIJ,JK)))) ), 10)  !inner min/max prevents sigfpe when 2*zq1 does not fit into an int
      ZINC(JIJ,JK) = 2.*ZQ1(JIJ,JK) - INQ1(JIJ,JK)

      PSIGRC(JIJ,JK) =  MIN(1.,(1.-ZINC(JIJ,JK))*ZSRC_1D(INQ1(JIJ,JK))+ZINC(JIJ,JK)*ZSRC_1D(INQ1(JIJ,JK)+1))
    IF(PRESENT(PHLC_HCF) .AND. PRESENT(PHLC_HRC))THEN
      PHLC_HCF(JIJ,JK)=0.
      PHLC_HRC(JIJ,JK)=0.
    ENDIF
    IF(PRESENT(PHLI_HCF) .AND. PRESENT(PHLI_HRI))THEN
      PHLI_HCF(JIJ,JK)=0.
      PHLI_HRI(JIJ,JK)=0.
    ENDIF
  END IF !HCONDENS

  IF(.NOT. OCND2) THEN
      PRC_OUT(JIJ,JK) = (1.-ZFRAC(JIJ,JK)) * ZCOND(JIJ,JK) ! liquid condensate
      PRI_OUT(JIJ,JK) = ZFRAC(JIJ,JK) * ZCOND(JIJ,JK)   ! solid condensate
      PT(JIJ,JK) = PT(JIJ,JK) + ((PRC_OUT(JIJ,JK)-PRC_IN(JIJ,JK))*ZLV(JIJ,JK) + &
                                    &(PRI_OUT(JIJ,JK)-PRI_IN(JIJ,JK))*ZLS(JIJ,JK)   ) &
                                  & /ZCPD(JIJ,JK)
      PRV_OUT(JIJ,JK) = ZRT(JIJ,JK) - PRC_OUT(JIJ,JK) - PRI_OUT(JIJ,JK)*ZPRIFACT(JIJ,JK)
  ELSE
      PRC_OUT(JIJ,JK) = (1.-ZFRAC(JIJ,JK)) * ZCOND(JIJ,JK) ! liquid condensate
      ZLWINC(JIJ,JK) = PRC_OUT(JIJ,JK) - PRC_IN(JIJ,JK)
      !
!     This check is mainly for noise reduction :
!     -------------------------
      IF(ABS(ZLWINC(JIJ,JK))>1.0E-12  .AND.  ESATW(ICEP%TIWMX, PT(JIJ,JK)) < PPABS(JIJ,JK)*0.5 )THEN
         ZRCOLD(JIJ,JK) = PRC_OUT(JIJ,JK)
         ZRFRAC(JIJ,JK) = PRV_IN(JIJ,JK) - ZLWINC(JIJ,JK)
         IF( PRV_IN(JIJ,JK) < ZRSW(JIJ,JK) )THEN ! sub - saturation over water:
            ! Avoid drying of cloudwater leading to supersaturation with
            ! respect to water
            ZRSDIF(JIJ,JK)= MIN(0.,ZRSP(JIJ,JK)-ZRFRAC(JIJ,JK))
         ELSE  ! super - saturation over water:
            ! Avoid deposition of water leading to sub-saturation with
            ! respect to water
            !            ZRSDIF(JIJ,JK)= MAX(0.,ZRSP(JIJ,JK)-ZRFRAC(JIJ,JK))
            ZRSDIF(JIJ,JK)= 0. ! t7
         ENDIF
         PRC_OUT(JIJ,JK) = ZCOND(JIJ,JK)  - ZRSDIF(JIJ,JK)
      ELSE
        ZRCOLD(JIJ,JK) = PRC_IN(JIJ,JK)
      ENDIF
 !   end check

 !    compute separate ice cloud:
      PWCLDFR(JIJ,JK) = PCLDFR(JIJ,JK)
      ZDUM1(JIJ,JK) = MIN(1.0,20.* PRC_OUT(JIJ,JK)*SQRT(ZDZ(JIJ,JK))/ZQSL(JIJ,JK)) ! cloud liquid water factor
      ZDUM3(JIJ,JK) = MAX(0.,PICLDFR(JIJ,JK)-PWCLDFR(JIJ,JK)) ! pure ice cloud part
      IF (JK==IKTB) THEN
        ZDUM4(JIJ,JK) = PRI_IN(JIJ,JK)
      ELSE
        ZDUM4(JIJ,JK) = PRI_IN(JIJ,JK) + PRS(JIJ,JK)*0.5 + PRG(JIJ,JK)*0.25
      ENDIF

      ZDUM4(JIJ,JK) = MAX(0.,MIN(1.,ZICE_CLD_WGT(JIJ,JK)*ZDUM4(JIJ,JK)*SQRT(ZDZ(JIJ,JK))/ZQSI(JIJ,JK))) ! clould ice+solid 
                                                         ! precip. water factor 

      ZDUM2(JIJ,JK) = (0.8*PCLDFR(JIJ,JK)+0.2)*MIN(1.,ZDUM1(JIJ,JK) + ZDUM4(JIJ,JK)*PCLDFR(JIJ,JK))
      ! water cloud, use 'statistical' cloud, but reduce it in case of low liquid content

      PCLDFR(JIJ,JK) = MIN(1., ZDUM2(JIJ,JK) + (0.5*ZDUM3(JIJ,JK)+0.5)*ZDUM4(JIJ,JK)) ! Rad cloud
           ! Reduce ice cloud part in case of low ice water content
      PRI_OUT(JIJ,JK) = PRI_IN(JIJ,JK)
      PT(JIJ,JK) = PT(JIJ,JK) + ((PRC_OUT(JIJ,JK)-ZRCOLD(JIJ,JK))*ZLV(JIJ,JK) + &
                                    &(PRI_OUT(JIJ,JK)-PRI_IN(JIJ,JK))*ZLS(JIJ,JK)   ) &
                                  & /ZCPD(JIJ,JK)
      PRV_OUT(JIJ,JK) = ZRT(JIJ,JK) - PRC_OUT(JIJ,JK) - PRI_OUT(JIJ,JK)*ZPRIFACT(JIJ,JK)
  END IF ! End OCND2
  IF(HLAMBDA3=='CB')THEN
      ! s r_c/ sig_s^2
      !    PSIGRC(JIJ,JK) = PCLDFR(JIJ,JK)  ! use simple Gaussian relation
      !
      !    multiply PSRCS by the lambda3 coefficient
      !
      !      PSIGRC(JIJ,JK) = 2.*PCLDFR(JIJ,JK) * MIN( 3. , MAX(1.,1.-ZQ1(JIJ,JK)) )
      ! in the 3D case lambda_3 = 1.

      PSIGRC(JIJ,JK) = PSIGRC(JIJ,JK)* MIN( 3. , MAX(1.,1.-ZQ1(JIJ,JK)) )
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
