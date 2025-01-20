!MNH_LIC Copyright 2011-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
     MODULE MODE_COMPUTE_MF_CLOUD_BIGAUS
!    ###################################
!
IMPLICIT NONE
CONTAINS
      SUBROUTINE COMPUTE_MF_CLOUD_BIGAUS(D, CST, PARAMMF, ICEP, KRR, &
                                  PRT_UP, PRV_UP, PRC_UP, PRI_UP, PTH_UP, PFRAC_UP, &
                                  PRTM, PTHM, PRM, &
                                  PRHODREF, PEXNM, PPABSM, &
                                  PRC_MF, PRI_MF, PCF_MF, PSIGMF, &
                                  PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF)
!     #################################################################
!!
!!****  *COMPUTE_MF_CLOUD_BIGAUS* -
!!       compute diagnostic subgrid cumulus cloud caracteristics with a statistical scheme
!!       based on a bi-gaussian PDF. In this routine, we only compute the shallow convection
!!       part of this bi-gaussian
!!
!!    PURPOSE
!!    -------
!!****  With this option, a formulation for the computation of the variance of the departure
!!      to saturation is proposed. This variance is used to compute the cloud fraction and
!!      the mean convective cloud content from the bi-gaussian PDF
!!
!
!!**  METHOD
!!    ------
!!      Updraft variables are used to diagnose the variance
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     REFERENCE
!!     ---------
!!
!!
!!     AUTHOR
!!     ------
!!     S. Riette moving of code previously in compute_mf_cloud code
!!
!!    MODIFICATIONS
!!    -------------
!!      Original 25 Aug 2011
!!      S. Riette Jan 2012: support for both order of vertical levels
!!      S. Riette Jun 2019: remove unused PRC_UP and PRI_UP, use SIGN in ERFC computation
!!      A. Marcel Jun 2025: new bi-gaussian cloud scheme implementation
!!      A. Marcel Jan 2025: bi-Gaussian PDF and associated subgrid precipitation
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST,             ONLY: CST_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
!
USE MODI_SHUMAN_MF, ONLY: MZF_MF, GZ_M_W_MF
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*                    0.1  Declaration of Arguments
!
TYPE(DIMPHYEX_t),                INTENT(IN)   :: D
TYPE(CST_t),                     INTENT(IN)   :: CST
TYPE(PARAM_MFSHALL_t),           INTENT(IN)   :: PARAMMF
TYPE(RAIN_ICE_PARAM_t),          INTENT(IN)   :: ICEP
INTEGER,                         INTENT(IN)   :: KRR                     ! number of moist var.
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PRT_UP, PRV_UP, PRC_UP, PRI_UP, PTH_UP, PFRAC_UP ! updraft characteritics
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PTHM, PRTM              ! env. var. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT, KRR),   INTENT(IN)   :: PRM                ! hydrometeors
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PRHODREF, PEXNM, PPABSM ! env density and pressure
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  :: PRC_MF, PRI_MF          ! cloud content
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  :: PCF_MF                  ! and cloud fraction for MF scheme
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  :: PSIGMF                  ! contribution to the env s variance
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  :: PHLC_HCF, PHLC_HRC, PHLI_HCF, PHLI_HRI ! low/high cloud diagnostics
!
!*                    0.1  Declaration of local variables
!
!
INTEGER                                    :: JK, JIJ  ! loop index
INTEGER                                    :: IIJE, IIJB, IKT ! bounds
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRT_UP_M,&
                               & ZRV_UP_M, ZRC_UP_M, ZRI_UP_M,&
                               & ZTH_UP_M,&
                               & ZFRAC_UP_M
REAL :: ZQ1, ZCPH, ZTEMP, ZTEMP_UP, ZLV, ZLS, ZLVS, ZQSL, ZAH, ZA, &
        ZSBAR_MEAN, ZSBAR_ENV, ZSIGS_UP, ZALPHA_UP_M, ZPV, ZQSI, &
        ZAUTC, ZGAUTC, ZGAUC, ZAUTI, ZGAUTI, ZGAUI, ZCRIAUTI, ZAUT, &
        ZCOND, ZGAM, ZFRAC_ICE_UP_M, ZFRAC_ICE, ZSBAR_UP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_BIGAUS',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
!Computation is done on mass points
!Interpolation on mass points
CALL MZF_MF(D, PRC_UP(:,:), ZRC_UP_M(:,:))
CALL MZF_MF(D, PRI_UP(:,:), ZRI_UP_M(:,:))
CALL MZF_MF(D, PRT_UP(:,:), ZRT_UP_M(:,:))
CALL MZF_MF(D, PRV_UP(:,:), ZRV_UP_M(:,:))
CALL MZF_MF(D, PFRAC_UP(:,:), ZFRAC_UP_M(:,:))
CALL MZF_MF(D, PTH_UP(:,:), ZTH_UP_M(:,:))

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    IF (ZFRAC_UP_M(JIJ, JK)>0.) THEN
      ! Ice franction in environment and in updraft
      ZFRAC_ICE_UP_M=0.
      IF (ZRC_UP_M(JIJ, JK)+ZRI_UP_M(JIJ, JK) > 1.E-20) THEN
        ZFRAC_ICE_UP_M = ZRI_UP_M(JIJ, JK) / (ZRC_UP_M(JIJ, JK)+ZRI_UP_M(JIJ, JK))
      ENDIF
      ZFRAC_ICE = 0.
      IF (PRM(JIJ, JK, 2)+PRM(JIJ, JK, 4) > 1.E-20) THEN
        ZFRAC_ICE = PRM(JIJ, JK, 4) / (PRM(JIJ, JK, 2)+PRM(JIJ, JK, 4))
      ENDIF

      !----------------------------------------------------------------------------
      !
      !*      1. Computation of the variance
      !          ------------------------------------------------
      !
      !
      !#### Computation of S in updraft
      ZCPH = CST%XCPD + CST%XCPV*ZRV_UP_M(JIJ, JK) + CST%XCL*ZRC_UP_M(JIJ, JK) + CST%XCI*ZRI_UP_M(JIJ, JK)
      ZTEMP_UP = ZTH_UP_M(JIJ, JK)*PEXNM(JIJ, JK)
      ZLV = CST%XLVTT + (CST%XCPV-CST%XCL) * (ZTEMP_UP-CST%XTT)
      ZLS = CST%XLSTT + (CST%XCPV-CST%XCI) * (ZTEMP_UP-CST%XTT)

      ZPV = MIN(EXP( CST%XALPW - CST%XBETAW / ZTEMP_UP - CST%XGAMW * LOG(ZTEMP_UP)), .99*PPABSM(JIJ, JK))
      ZQSL = CST%XRD / CST%XRV * ZPV / ( PPABSM(JIJ, JK) - ZPV )
      ZPV = MIN(EXP( CST%XALPI - CST%XBETAI / ZTEMP_UP - CST%XGAMI * LOG(ZTEMP_UP)), .99*PPABSM(JIJ, JK))
      ZQSI = CST%XRD / CST%XRV * ZPV / ( PPABSM(JIJ, JK) - ZPV )

      ZQSL = (1. - ZFRAC_ICE_UP_M) * ZQSL + ZFRAC_ICE_UP_M * ZQSI
      ZLVS = (1. - ZFRAC_ICE_UP_M) * ZLV + ZFRAC_ICE_UP_M * ZLS

      ZAH = ZLVS * ZQSL / ( CST%XRV * ZTEMP_UP**2 ) * (1. + (CST%XRV * ZQSL / CST%XRD))
      ZA = 1. / ( 1. + ZLVS/ZCPH * ZAH )

      ZSBAR_UP = ZA * (ZRT_UP_M(JIJ, JK) - ZQSL + ZAH * ZLVS * (ZRC_UP_M(JIJ, JK)+ZRI_UP_M(JIJ, JK)) / ZCPH)

      !#### Computation of S on mean grid
      ZCPH = CST%XCPD + CST%XCPV*PRM(JIJ, JK, 1) + CST%XCL*PRM(JIJ, JK, 2) + CST%XCI*PRM(JIJ, JK, 4)
      ZTEMP = PTHM(JIJ, JK)*PEXNM(JIJ, JK)
      ZLV = CST%XLVTT + (CST%XCPV-CST%XCL) * (ZTEMP-CST%XTT)
      ZLS = CST%XLSTT + (CST%XCPV-CST%XCI) * (ZTEMP-CST%XTT)

      ZPV = MIN(EXP( CST%XALPW - CST%XBETAW / ZTEMP - CST%XGAMW * LOG(ZTEMP)), .99*PPABSM(JIJ, JK))
      ZQSL = CST%XRD / CST%XRV * ZPV / ( PPABSM(JIJ, JK) - ZPV )
      ZPV = MIN(EXP( CST%XALPI - CST%XBETAI / ZTEMP - CST%XGAMI * LOG(ZTEMP)), .99*PPABSM(JIJ, JK))
      ZQSI = CST%XRD / CST%XRV * ZPV / ( PPABSM(JIJ, JK) - ZPV )

      ZQSL = (1. - ZFRAC_ICE) * ZQSL + ZFRAC_ICE * ZQSI
      ZLVS = (1. - ZFRAC_ICE) * ZLV + ZFRAC_ICE * ZLS

      ZAH = ZLVS * ZQSL / ( CST%XRV * ZTEMP**2 ) * (1. + (CST%XRV * ZQSL / CST%XRD))
      ZA = 1. / ( 1. + ZLVS/ZCPH * ZAH )

      ZSBAR_MEAN = ZA * (PRTM(JIJ, JK) - ZQSL + ZAH * ZLVS * (PRM(JIJ, JK, 2)+PRM(JIJ, JK, 4)) / ZCPH)

      !#### Computation of S on environment
      ZSBAR_ENV = (ZSBAR_MEAN-ZFRAC_UP_M(JIJ, JK)*ZSBAR_UP)/(1.-ZFRAC_UP_M(JIJ, JK))

      !#### Compute parametrized variances
      ! Contribution to the environment internal variance
      PSIGMF(JIJ, JK)=PARAMMF%XSIGMA_ENV*(1. -ZFRAC_UP_M(JIJ, JK))*(ZSBAR_ENV-ZSBAR_MEAN)**2

      ! Updraft internal variance
      ZSIGS_UP =PARAMMF%XSIGMA_MF*ZFRAC_UP_M(JIJ, JK)*(ZSBAR_UP-ZSBAR_MEAN)**2
      ZSIGS_UP = SQRT(MAX(ABS(ZSIGS_UP), 1.E-40))

      !*      2. PDF integration
      !          ------------------------------------------------
      !
      ZALPHA_UP_M = ZFRAC_UP_M(JIJ, JK)*PARAMMF%XALPHA_MF ! XALPHA_MF = 1 --> EDMF :: XALPHA_MF > 1 --> PARAM
      ZALPHA_UP_M = MAX(0., MIN(ZALPHA_UP_M, 1.))

      !The mean of the distribution is ZSBAR_UP
      !Computation of ZA and ZGAM (=efrc(ZA)) coefficient

      !Gaussian computation for S variable
      ZQ1 = ZSBAR_UP/ZSIGS_UP
      ZA = -ZQ1/SQRT(2.)

      ZGAM = 1 + ERF(-ZA)
      ZCOND = (EXP(-ZA**2)-ZA*SQRT(CST%XPI)*ZGAM)*ZSIGS_UP/SQRT(2.*CST%XPI) * ZALPHA_UP_M

      !Cloud fraction
      PCF_MF(JIJ, JK) = MAX(0., MIN(1., 0.5*ZGAM * ZALPHA_UP_M))
      !computation of condensate, then PRC and PRI
      ZCOND = MAX(ZCOND, 0.)
      PRC_MF(JIJ, JK) = (1.-ZFRAC_ICE_UP_M) * ZCOND
      PRI_MF(JIJ, JK) = ZFRAC_ICE_UP_M * ZCOND

      ! For low/high cloud diagnostics
      IF(1-ZFRAC_ICE_UP_M > 1.E-20)THEN
        ZAUTC = (ZSBAR_UP - ICEP%XCRIAUTC/(PRHODREF(JIJ, JK)*(1-ZFRAC_ICE_UP_M)))/ZSIGS_UP
        ZGAUTC = -ZAUTC/SQRT(2.)

        ZGAUC = 1 + ERF(-ZGAUTC)
        PHLC_HRC(JIJ, JK) = (1-ZFRAC_ICE_UP_M)*(EXP(-ZGAUTC**2)-ZGAUTC*SQRT(CST%XPI)*ZGAUC)*ZSIGS_UP/SQRT(2.*CST%XPI) &
         * ZALPHA_UP_M

        PHLC_HCF(JIJ, JK) = MAX( 0., MIN(1., 0.5*ZGAUC)) * ZALPHA_UP_M
        IF (PHLC_HCF(JIJ, JK)<1.E-6)THEN !Protection for subgrid rain using HSUBG_RC_RR_ACCR=='PRFR' (valid with RC)
          PHLC_HCF(JIJ, JK) = 0.
        ENDIF
        PHLC_HRC(JIJ, JK) = PHLC_HRC(JIJ, JK) + ICEP%XCRIAUTC/PRHODREF(JIJ, JK) * PHLC_HCF(JIJ, JK)
        PHLC_HRC(JIJ, JK) = MAX(PHLC_HRC(JIJ, JK), 0.)
      ELSE
        PHLC_HCF(JIJ, JK)=0.
        PHLC_HRC(JIJ, JK)=0.
      ENDIF

      IF(ZFRAC_ICE_UP_M > 1.E-20)THEN
        ZCRIAUTI = MIN(ICEP%XCRIAUTI, 10**(ICEP%XACRIAUTI*(ZTEMP_UP-CST%XTT)+ICEP%XBCRIAUTI))
        ZAUTI = (ZSBAR_UP - ZCRIAUTI/ZFRAC_ICE_UP_M)/ZSIGS_UP
        ZGAUTI = -ZAUTI/SQRT(2.)

        ZGAUI = 1 + ERF(-ZGAUTI)
        PHLI_HRI(JIJ, JK) = ZFRAC_ICE_UP_M*(EXP(-ZGAUTI**2)-ZGAUTI*SQRT(CST%XPI)*ZGAUI)*ZSIGS_UP/SQRT(2.*CST%XPI) &
        * ZALPHA_UP_M

        PHLI_HCF(JIJ, JK) = MAX( 0., MIN(1., 0.5*ZGAUI)) * ZALPHA_UP_M
        IF (PHLI_HCF(JIJ, JK)<1.E-6)THEN !Protection for subgrid rain using HSUBG_RC_RR_ACCR=='PRFR' (valid with RC)
          PHLI_HCF(JIJ, JK) = 0.
        ENDIF
        PHLI_HRI(JIJ, JK) = PHLI_HRI(JIJ, JK) + ZCRIAUTI*PHLI_HCF(JIJ, JK)
        PHLI_HRI(JIJ, JK) = MAX(PHLI_HRI(JIJ, JK), 0.)
      ELSE
        PHLI_HCF(JIJ, JK)=0.
        PHLI_HRI(JIJ, JK)=0.
      ENDIF
    ELSE
      PCF_MF(JIJ, JK) = 0.
      PRC_MF(JIJ, JK) = 0.
      PRI_MF(JIJ, JK) = 0.
      PHLC_HCF(JIJ, JK) = 0.
      PHLC_HRC(JIJ, JK) = 0.
      PHLI_HCF(JIJ, JK) = 0.
      PHLI_HRI(JIJ, JK) = 0.
    ENDIF
  ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_BIGAUS',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_MF_CLOUD_BIGAUS
END MODULE MODE_COMPUTE_MF_CLOUD_BIGAUS
