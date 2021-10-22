!     ######spl
    SUBROUTINE CONDENSATION( KIU, KJU, KKU, KIB, KIE, KJB, KJE, KKB, KKE, KKL,         &
       HFRAC_ICE,                                                                      &
       PPABS, PZZ, PT, PRV, PRC, PRI, PRS, PRG, PSIGS, PMFCONV, PCLDFR, PSIGRC, OUSERI,&
       OSIGMAS, OCND2, PSIGQSAT, YSPP_PSIGQSAT,YSPP_ICE_CLD_WGT,                       &
       PLV, PLS, PCPH  )
    USE PARKIND1, ONLY : JPRB
    USE MODD_SPP_TYPE
    USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
!!      2016-11 S. Riette: use HFRAC_ICE, output adjusted state
!!      2020-12 U. Andrae : Introduce SPP for HARMONIE-AROME
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAMETERS
USE MODI_COMPUTE_FRAC_ICE
USE MODD_TIWMX
USE MODI_ICECLOUD
USE MODI_TIWMX
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
CHARACTER*1,                  INTENT(IN)    :: HFRAC_ICE
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PPABS  ! pressure (Pa)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PZZ    ! height of model levels (m)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PT     ! grid scale T  (K)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
LOGICAL, INTENT(IN)                         :: OUSERI ! logical switch to compute both
						      ! liquid and solid condensate (OUSERI=.TRUE.)
                                                      ! or only solid condensate (OUSERI=.FALSE.)
LOGICAL, INTENT(IN)                         :: OCND2  ! logical switch to sparate liquid and ice
                                                      ! more rigid (DEFALT value : .FALSE.)
LOGICAL, INTENT(IN)                         :: OSIGMAS! use present global Sigma_s values
                                                      ! or that from turbulence scheme
REAL, INTENT(IN)                            :: PSIGQSAT ! use an extra "qsat" variance contribution (OSIGMAS case)
                                                        ! multiplied by PSIGQSAT
TYPE(TSPP_CONFIG_MPA), INTENT(IN)           :: YSPP_PSIGQSAT,YSPP_ICE_CLD_WGT
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRC    ! grid scale r_c mixing ratio (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PRI    ! grid scale r_i (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRS    ! grid scale mixing ration of snow (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRG    ! grid scale mixing ration of graupel (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PSIGS  ! Sigma_s from turbulence scheme
REAL, DIMENSION(:,:,:),       INTENT(IN)    :: PMFCONV! convective mass flux (kg /s m^2)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PCLDFR ! cloud fraction
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PSIGRC ! s r_c / sig_s^2
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PLV
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PLS
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PCPH
!
!
!*       0.2   Declarations of local variables :
!
INTEGER  :: JI, JJ, JK, JKP, JKM, IKTB, IKTE    ! loop index
REAL, DIMENSION(KIU,KJU,KKU) :: ZTLK, ZRT       ! work arrays for T_l and total water mixing ratio
REAL, DIMENSION(KIU,KJU,KKU) :: ZL              ! length scale
REAL, DIMENSION(KIU,KJU,KKU) :: ZFRAC           ! Ice fraction
INTEGER, DIMENSION(KIU,KJU)  :: ITPL            ! top levels of troposphere 
REAL,    DIMENSION(KIU,KJU)  :: ZTMIN           ! minimum Temp. related to ITPL
!
REAL, DIMENSION(KIU,KJU,KKU) :: ZLV, ZLS, ZCPD
REAL :: ZTEMP, ZPV, ZQSL, ZPIV, ZQSI, ZCOND, ZLVS ! thermodynamics
REAL :: ZLL, DZZ, ZZZ                           ! used for length scales 
REAL :: ZAH, ZA, ZB, ZSBAR, ZQ1, ZSIGMA, ZDRW, ZDTL, ZSIG_CONV ! related to computation of Sig_s
REAL :: ZRIOLD

! related to OCND2 noise check :
REAL :: ZRCOLD,ZTSP,ZRSP,  ZRSW, ZRFRAC,ZRSDIF ,ZRSB
! related to OCND2  ice cloud calulation :
REAL :: ZDUM1,ZDUM2,ZDUM3,ZDUM4,ZPRIFACT
REAL, DIMENSION(KIU,KJU,KKU) :: TCLD
REAL :: ZDZ(KIB:KIE), ZARDUM(KIE-KIB+1),ZCLDUM(KIE-KIB+1)
REAL, DIMENSION(KIB:KIE) :: ZSIGQSAT,ZICE_CLD_WGT
! end OCND2

!
!*       0.3  Definition of constants :
!
!-------------------------------------------------------------------------------
!
REAL :: ZL0     = 600.        ! tropospheric length scale
REAL :: ZCSIGMA = 0.2         ! constant in sigma_s parameterization
REAL :: ZCSIG_CONV = 0.30E-2  ! scaling factor for ZSIG_CONV as function of mass flux
!

REAL, DIMENSION(-22:11) :: ZSRC_1D =(/                                   &
       0.           ,  0.           ,  2.0094444E-04,   0.316670E-03,    &
       4.9965648E-04,  0.785956E-03 ,  1.2341294E-03,   0.193327E-02,    &
       3.0190963E-03,  0.470144E-02 ,  7.2950651E-03,   0.112759E-01,    &
       1.7350994E-02,  0.265640E-01 ,  4.0427860E-02,   0.610997E-01,    &
       9.1578111E-02,  0.135888E+00 ,  0.1991484    ,   0.230756E+00,    &
       0.2850565    ,  0.375050E+00 ,  0.5000000    ,   0.691489E+00,    &
       0.8413813    ,  0.933222E+00 ,  0.9772662    ,   0.993797E+00,    &
       0.9986521    ,  0.999768E+00 ,  0.9999684    ,   0.999997E+00,    &
       1.0000000    ,  1.000000     /)
INTEGER  :: INQ1
REAL :: ZINC,ZMU,ZVAL
!
!-------------------------------------------------------------------------------
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CONDENSATION',0,ZHOOK_HANDLE)

IKTB=1+JPVEXT
IKTE=KKU-JPVEXT

PCLDFR(:,:,:) = 0. ! Initialize values
PSIGRC(:,:,:) = 0. ! Initialize values
ZPRIFACT = 1.      ! Initialize value
ZCLDUM=-1.         ! Initialize value
IF(OCND2)ZPRIFACT = 0.

!
! Perturb PSIGQSAT
!

IF (YSPP_PSIGQSAT%LPERT) THEN
  IF (YSPP_PSIGQSAT%LLNN_MEAN1.OR.YSPP_PSIGQSAT%LLNN_MEAN1_SELF) THEN
    ZMU = -0.5_JPRB * (YSPP_PSIGQSAT%CMPERT * YSPP_PSIGQSAT%SDEV)**2
  ELSE
    ZMU = 0._JPRB
  ENDIF
  DO JI=KIB,KIE
    ZVAL = &
      PSIGQSAT*EXP(ZMU+YSPP_PSIGQSAT%CMPERT*YSPP_PSIGQSAT%PGP2DSPP(JI)) 
    ZSIGQSAT(JI) = MAX(YSPP_PSIGQSAT%CLIP(1),MIN(ZVAL,YSPP_PSIGQSAT%CLIP(2)))
  ENDDO
ELSE
  ZSIGQSAT(KIB:KIE) = PSIGQSAT
ENDIF

!
! Perturb ICE_CLD_WGT
!

IF (YSPP_ICE_CLD_WGT%LPERT) THEN
  IF (YSPP_ICE_CLD_WGT%LLNN_MEAN1.OR.YSPP_ICE_CLD_WGT%LLNN_MEAN1_SELF) THEN
    ZMU = -0.5_JPRB * (YSPP_ICE_CLD_WGT%CMPERT * YSPP_ICE_CLD_WGT%SDEV)**2
  ELSE
    ZMU = 0._JPRB
  ENDIF
  DO JI=KIB,KIE
! Awaiting HARMONIE-AROME physics changes
!    ZVAL = &
!     XFRMIN(21)* EXP(ZMU+YSPP_ICE_CLD_WGT%CMPERT*YSPP_ICE_CLD_WGT%PGP2DSPP(JI)) 
     ZVAL = &
      1.5* EXP(ZMU+YSPP_ICE_CLD_WGT%CMPERT*YSPP_ICE_CLD_WGT%PGP2DSPP(JI)) 
    ZICE_CLD_WGT(JI) = &
    MAX(YSPP_ICE_CLD_WGT%CLIP(1),MIN(ZVAL,YSPP_ICE_CLD_WGT%CLIP(2)))
  ENDDO
ELSE
!  ZICE_CLD_WGT(KIB:KIE) = XFRMIN(21)
   ZICE_CLD_WGT(KIB:KIE) = 1.5
ENDIF

!
!
!-------------------------------------------------------------------------------
! store total water mixing ratio
DO JK=IKTB,IKTE
  DO JJ=KJB,KJE
    DO JI=KIB,KIE
      ZRT(JI,JJ,JK)  = PRV(JI,JJ,JK) + PRC(JI,JJ,JK) + PRI(JI,JJ,JK)*ZPRIFACT
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
! Preliminary calculations needed for computing the "turbulent part" of Sigma_s
IF ( .NOT. OSIGMAS ) THEN
  DO JK=IKTB,IKTE
    DO JJ=KJB,KJE
      DO JI=KIB,KIE
        ZTEMP  = PT(JI,JJ,JK)
        ! store temperature at saturation
        ZTLK(JI,JJ,JK) = ZTEMP - ZLV(JI,JJ,JK)*PRC(JI,JJ,JK)/ZCPD(JI,JJ,JK) &
                               - ZLS(JI,JJ,JK)*PRI(JI,JJ,JK)/ZCPD(JI,JJ,JK)*ZPRIFACT
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
WHERE(PRC(:,:,:)+PRI(:,:,:) > 1.E-20)
  ZFRAC(:,:,:) = PRI(:,:,:) / (PRC(:,:,:)+PRI(:,:,:))
ENDWHERE
CALL COMPUTE_FRAC_ICE(HFRAC_ICE, ZFRAC, PT)
IF(.NOT. OUSERI) ZFRAC(:,:,:)=0.
IF(OCND2) ZFRAC=0.
!
DO JK=IKTB,IKTE
  JKP=MAX(MIN(JK+KKL,IKTE),IKTB)
  JKM=MAX(MIN(JK-KKL,IKTE),IKTB) 
  DO JJ=KJB,KJE
    IF (OCND2) THEN
       ZDZ(KIB:KIE) = PZZ(KIB:KIE,JJ,JKP) - PZZ(KIB:KIE,JJ,JKP+1)
       CALL ICECLOUD(KIE-KIB+1,PPABS(KIB,JJ,JK),PZZ(KIB,JJ,JK),ZDZ(KIB), &
            & PT(KIB,JJ,JK),PRV(KIB,JJ,JK),1.,-1., &
            & ZCLDUM,1.,TCLD(KIB,JJ,JK), &
            & ZARDUM,ZARDUM,ZARDUM,ZARDUM)
    ENDIF
    DO JI=KIB,KIE
      ! latent heats
      ZTEMP  = PT(JI,JJ,JK)
      ! saturated water vapor mixing ratio over liquid water
      IF (OCND2) THEN
         ZPV    = MIN(ESATW(ZTEMP), .99*PPABS(JI,JJ,JK))
      ELSE
         ZPV    = MIN(EXP( XALPW - XBETAW / ZTEMP - XGAMW * LOG( ZTEMP ) ), .99*PPABS(JI,JJ,JK))
      ENDIF
      ZQSL   = XRD / XRV * ZPV / ( PPABS(JI,JJ,JK) - ZPV )

      ! saturated water vapor mixing ratio over ice
      IF (OCND2) THEN
         ZPIV   = MIN(ESATI(ZTEMP), .99*PPABS(JI,JJ,JK))
      ELSE
         ZPIV   = MIN(EXP( XALPI - XBETAI / ZTEMP - XGAMI * LOG( ZTEMP ) ), .99*PPABS(JI,JJ,JK))
      ENDIF
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
                   & ZAH * ZLVS * (PRC(JI,JJ,JK)+PRI(JI,JJ,JK)*ZPRIFACT) / ZCPD(JI,JJ,JK))

      ! switch to take either present computed value of SIGMAS
      ! or that of Meso-NH turbulence scheme
      IF ( OSIGMAS ) THEN
        IF (ZSIGQSAT(JI)/=0.) THEN
          ZSIGMA = SQRT((2*PSIGS(JI,JJ,JK))**2 + (ZSIGQSAT(JI)*ZQSL*ZA)**2)
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

      ! normalized saturation deficit
      ZQ1   = ZSBAR/ZSIGMA

      ! cloud fraction
      PCLDFR(JI,JJ,JK) = MAX( 0., MIN(1.,0.5+0.36*ATAN(1.55*ZQ1)) )

      ! total condensate
      IF (ZQ1 > 0. .AND. ZQ1 <= 2 ) THEN
        ZCOND = MIN(EXP(-1.)+.66*ZQ1+.086*ZQ1**2, 2.) ! We use the MIN function for continuity
      ELSE IF (ZQ1 > 2.) THEN
        ZCOND = ZQ1
      ELSE
        ZCOND = EXP( 1.2*ZQ1-1. )
      ENDIF
      ZCOND = ZCOND * ZSIGMA

      IF ( ZCOND < 1.E-12 ) THEN
        ZCOND = 0.
        PCLDFR(JI,JJ,JK) = 0.
      ENDIF
      IF (PCLDFR(JI,JJ,JK)==0.) THEN
        ZCOND=0.
      ENDIF

      ZRCOLD=PRC(JI,JJ,JK)
      ZRIOLD=PRI(JI,JJ,JK)

      IF(.NOT. OCND2) THEN
        PRC(JI,JJ,JK) = (1.-ZFRAC(JI,JJ,JK)) * ZCOND ! liquid condensate
        PRI(JI,JJ,JK) = ZFRAC(JI,JJ,JK) * ZCOND   ! solid condensate
      ELSE
         PRC(JI,JJ,JK) = (1.-ZFRAC(JI,JJ,JK)) * ZCOND ! liquid condensate
         !
!        This check is mainly for noise reduction :
!        -------------------------
         IF(ABS(ZRCOLD-PRC(JI,JJ,JK))>1.0E-12 .AND. ESATW(ZTEMP) < PPABS(JI,JJ,JK)*0.5)THEN
            ZRCOLD=PRC(JI,JJ,JK)
            ZTSP= TIWMX_TAB(PPABS(JI,JJ,JK),PT(JI,JJ,JK),PRV(JI,JJ,JK),0.,ZRSP,ZRSW,0.1)

            ZRFRAC =    PRV(JI,JJ,JK) - ZCOND + ZRCOLD
            IF( PRV(JI,JJ,JK) < ZRSW )THEN ! sub - saturation over water:
               ! Avoid drying of cloudwater leading to supersaturation with
               ! respect to water
! 
               ZRSDIF= MIN(0.,ZRSP-ZRFRAC)
            ELSE  ! super - saturation over water:
               ! Avoid depostition of water leading to sub-saturation with
               ! respect to water
               !            ZRSDIF= MAX(0.,ZRSP-ZRFRAC)
               ZRSDIF= MAX(0.,ZRSP*PCLDFR(JI,JJ,JK) - ZRFRAC) 
            ENDIF
            ZRSB= ZCOND  - ZRSDIF

            PRC(JI,JJ,JK) =  ZRSB
         ENDIF
 !       end check 

 !       compute separate ice cloud:
         ZZZ = ZDZ(JI)

         ZDUM1 = MIN(1.0,20.* PRC(JI,JJ,JK)*SQRT(ZZZ)/ZQSL) ! clould liquid water 
                                                        ! factor 
 
         ZDUM2 = (0.8*PCLDFR(JI,JJ,JK)+0.2)*MIN(1.,ZDUM1 + ZDUM4*PCLDFR(JI,JJ,JK))
         ! water cloud, use 'statistical' cloud, but reduce it in case of low liquid content

         ZDUM3 = MAX(0.,TCLD(JI,JJ,JK)-PCLDFR(JI,JJ,JK)) ! pure ice cloud part

         ZDUM4 = PRI(JI,JJ,JK) + PRS(JI,JJ,JK)*0.5 + PRG(JI,JJ,JK)*0.25 
         IF (JK==IKTB) ZDUM4 = PRI(JI,JJ,JK)

         ZDUM4 = MAX(0.,MIN(1.,ZICE_CLD_WGT(JI)*ZDUM4*SQRT(ZZZ)/ZQSI)) ! clould ice+solid 
                                                            ! precip. water factor 
         PCLDFR(JI,JJ,JK) = MIN(1., ZDUM2 + (0.9*ZDUM3+0.1)*ZDUM4) ! Rad cloud
              ! Reduce ice cloud part in case of low ice water content
      ENDIF ! End OCND2 
      
      PT(JI,JJ,JK) = PT(JI,JJ,JK) + ((PRC(JI,JJ,JK)-ZRCOLD)*ZLV(JI,JJ,JK) + &
                                    &(PRI(JI,JJ,JK)-ZRIOLD)*ZLS(JI,JJ,JK)   ) &
                                  & /ZCPD(JI,JJ,JK)
      PRV(JI,JJ,JK) = ZRT(JI,JJ,JK) - PRC(JI,JJ,JK) - PRI(JI,JJ,JK)*ZPRIFACT
     
! s r_c/ sig_s^2
!    PSIGRC(JI,JJ,JK) = PCLDFR(JI,JJ,JK)  ! use simple Gaussian relation
!
!    multiply PSRCS by the lambda3 coefficient
!
!      PSIGRC(JI,JJ,JK) = 2.*PCLDFR(JI,JJ,JK) * MIN( 3. , MAX(1.,1.-ZQ1) )
! in the 3D case lambda_3 = 1.
      INQ1 = MIN( MAX(-22,FLOOR(MIN(100., MAX(-100., 2*ZQ1))) ), 10) !inner min/max prevents sigfpe when 2*zq1 does not fit into an int
      ZINC = 2.*ZQ1 - INQ1

      PSIGRC(JI,JJ,JK) =  MIN(1.,(1.-ZINC)*ZSRC_1D(INQ1)+ZINC*ZSRC_1D(INQ1+1))

      PSIGRC(JI,JJ,JK) = PSIGRC(JI,JJ,JK)* MIN( 3. , MAX(1.,1.-ZQ1) )

    END DO
  END DO
END DO
!
IF (LHOOK) CALL DR_HOOK('CONDENSATION',1,ZHOOK_HANDLE)
END SUBROUTINE CONDENSATION
