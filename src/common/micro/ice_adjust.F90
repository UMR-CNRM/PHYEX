!MNH_LIC Copyright 1996-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################################################################
      SUBROUTINE ICE_ADJUST (D, CST, ICEP, NEB, BUCONF, KRR,                   &
                             HFRAC_ICE, HCONDENS, HLAMBDA3,&
                             HBUNAME, OSUBG_COND, OSIGMAS, OCND2, HSUBG_MF_PDF,&
                             PTSTEP, PSIGQSAT,                                 &
                             PRHODJ, PEXNREF, PRHODREF, PSIGS, LMFCONV, PMFCONV,&
                             PPABST, PZZ,                                      &
                             PEXN, PCF_MF, PRC_MF, PRI_MF,                     &
                             PRV, PRC, PRVS, PRCS, PTH, PTHS, PSRCS, PCLDFR,   &
                             PRR, PRI, PRIS, PRS, PRG, PRH,                    &
                             POUT_RV, POUT_RC, POUT_RI, POUT_TH,               &
                             PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,           &
                             TBUDGETS, KBUDGETS,                               &
                             PICE_CLD_WGT)
!     #########################################################################
!
!!****  *ICE_ADJUST* -  compute the ajustment of water vapor in mixed-phase 
!!                      clouds
!!
!!    PURPOSE
!!    -------
!!    The purpose of this routine is to compute the fast microphysical sources
!!    through a saturation ajustement procedure in case of mixed-phase clouds.
!!
!!
!!**  METHOD
!!    ------
!!    Langlois, Tellus, 1973 for the cloudless version.
!!    When cloud water is taken into account, refer to book 1 of the
!!    documentation.
!!
!!     
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!         XP00               ! Reference pressure
!!         XMD,XMV            ! Molar mass of dry air and molar mass of vapor
!!         XRD,XRV            ! Gaz constant for dry air, gaz constant for vapor
!!         XCPD,XCPV          ! Cpd (dry air), Cpv (vapor)
!!         XCL                ! Cl (liquid)
!!         XCI                ! Ci (ice)
!!         XTT                ! Triple point temperature
!!         XLVTT              ! Vaporization heat constant
!!         XLSTT              ! Sublimation  heat constant
!!         XALPW,XBETAW,XGAMW ! Constants for saturation vapor over liquid
!!                            !  pressure  function 
!!         XALPI,XBETAI,XGAMI ! Constants for saturation vapor over ice
!!                            !  pressure  function 
!!      Module  MODD_CONF 
!!         CCONF
!!      Module MODD_BUDGET:
!!         NBUMOD 
!!         CBUTYPE
!!         LBU_RTH    
!!         LBU_RRV  
!!         LBU_RRC  
!!         LBU_RRI  
!!
!!
!!    REFERENCE
!!    ---------
!!      Book 1 and Book2 of documentation ( routine ICE_ADJUST )
!!      Langlois, Tellus, 1973
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty    * Laboratoire d'Aerologie*
!!   
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/12/96 
!!      M. Tomasini 27/11/00 Change CND and DEP fct of the T instead of rc and ri
!!                           Avoid the sub- and super-saturation before the ajustment
!!                           Avoid rc>0 if T<T00 before the ajustment
!!      P Bechtold 12/02/02  change subgrid condensation
!!      JP Pinty   29/11/02  add ICE2 and IC4 cases
!!      (P. Jabouille) 27/05/04 safety test for case where esw/i(T)> pabs (~Z>40km)
!!      J.Pergaud and S.Malardel Add EDKF case
!!      S. Riette ice for EDKF
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!      2016-07 S. Riette: adjustement is now realized on state variables (PRV, PRC, PRI, PTH)
!!                         whereas tendencies are still applied on S variables.
!!                         This modification allows to call ice_adjust on T variable
!!                         or to call it on S variables
!!      2016-11 S. Riette: all-or-nothing adjustment now uses condensation
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!!      2020-12 U. Andrae : Introduce SPP for HARMONIE-AROME
!!     R. El Khatib 24-Aug-2021 Optimizations
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_CST,        ONLY: CST_t
USE MODD_NEB,        ONLY: NEB_t
USE MODD_BUDGET,     ONLY: TBUDGETDATA, TBUDGETCONF_t, NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RI
USE MODD_RAIN_ICE_PARAM, ONLY : RAIN_ICE_PARAM_t
!
USE MODE_BUDGET,         ONLY: BUDGET_STORE_INIT, BUDGET_STORE_END
!
USE MODI_CONDENSATION
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments :
!
!
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(NEB_t),              INTENT(IN)    :: NEB
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER(LEN=1),         INTENT(IN)    :: HFRAC_ICE
CHARACTER(LEN=80),        INTENT(IN)    :: HCONDENS
CHARACTER(LEN=4),         INTENT(IN)    :: HLAMBDA3 ! formulation for lambda3 coeff
CHARACTER(LEN=4),         INTENT(IN)    :: HBUNAME  ! Name of the budget
LOGICAL,                  INTENT(IN)    :: OSUBG_COND ! Switch for Subgrid 
                                                    ! Condensation
LOGICAL                                 :: OSIGMAS  ! Switch for Sigma_s: 
                                                    ! use values computed in CONDENSATION
                                                    ! or that from turbulence scheme
LOGICAL                                 :: OCND2    ! logical switch to sparate liquid 
                                                    ! and ice
                                                    ! more rigid (DEFALT value : .FALSE.)
CHARACTER(LEN=80),        INTENT(IN)    :: HSUBG_MF_PDF
REAL,                     INTENT(IN)   :: PTSTEP    ! Double Time step
                                                    ! (single if cold start)
REAL, DIMENSION(D%NIT,D%NJT),                INTENT(IN)    :: PSIGQSAT  ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    ::  PEXNREF ! Reference Exner function
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    ::  PRHODREF
!
REAL, DIMENSION(MERGE(D%NIT,0,OSUBG_COND),&
                MERGE(D%NJT,0,OSUBG_COND),&
                MERGE(D%NKT,0,OSUBG_COND)),           INTENT(IN)    ::  PSIGS   ! Sigma_s at time t
LOGICAL,                                                       INTENT(IN)    ::  LMFCONV ! =SIZE(PMFCONV)!=0
REAL, DIMENSION(MERGE(D%NIT,0,LMFCONV),&
                MERGE(D%NJT,0,LMFCONV),&
                MERGE(D%NKT,0,LMFCONV)),              INTENT(IN)   ::  PMFCONV ! convective mass flux
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    ::  PPABST  ! Absolute Pressure at t        
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    ::  PZZ     ! height of model layer
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    ::  PEXN    ! Exner function
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PRI_MF! Convective Mass Flux ice mixing ratio
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    ::  PRV     ! Water vapor m.r. to adjust
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    ::  PRC     ! Cloud water m.r. to adjust
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PTH     ! Theta to adjust
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                                                        ! s'rc'/2Sigma_s2 at time t+1
                                                                                        ! multiplied by Lambda_3
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)   :: PCLDFR  ! Cloud fraction          
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT)::  PRIS ! Cloud ice  m.r. at t+1
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRR  ! Rain water m.r. to adjust
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRI  ! Cloud ice  m.r. to adjust
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRS  ! Aggregate  m.r. to adjust
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   ::  PRG  ! Graupel    m.r. to adjust
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)   ::  PRH  ! Hail       m.r. to adjust
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  POUT_RV ! Adjusted value
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  POUT_RC ! Adjusted value
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  POUT_RI ! Adjusted value
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  POUT_TH ! Adjusted value
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  PHLC_HRC
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  PHLC_HCF
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  PHLI_HRI
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  PHLI_HCF
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS),                        INTENT(INOUT) :: TBUDGETS
INTEGER,                                                       INTENT(IN)    :: KBUDGETS
REAL, DIMENSION(D%NIT,D%NJT),                OPTIONAL, INTENT(IN)   :: PICE_CLD_WGT
!
!*       0.2   Declarations of local variables :
!
!
REAL  :: ZW1,ZW2    ! intermediate fields
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) &
                         :: ZT,   &  ! adjusted temperature
                   ZRV, ZRC, ZRI, &  ! adjusted state
                            ZCPH, &  ! guess of the CPh for the mixing
                            ZLV,  &  ! guess of the Lv at t+1
                            ZLS      ! guess of the Ls at t+1
REAL :: ZCRIAUT, & ! Autoconversion thresholds
        ZHCF, ZHR
!
INTEGER             :: JITER,ITERMAX ! iterative loop for first order adjustment
INTEGER             :: JI, JJ, JK
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZSIGS, ZSRCS
REAL, DIMENSION(D%NIT,D%NJT) :: ZSIGQSAT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
IF (LHOOK) CALL DR_HOOK('ICE_ADJUST',0,ZHOOK_HANDLE)
!
ITERMAX=1
!
IF(BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT(TBUDGETS(NBUDGET_TH), TRIM(HBUNAME), PTHS(:, :, :)*PRHODJ(:, :, :))
IF(BUCONF%LBUDGET_RV) CALL BUDGET_STORE_INIT(TBUDGETS(NBUDGET_RV), TRIM(HBUNAME), PRVS(:, :, :)*PRHODJ(:, :, :))
IF(BUCONF%LBUDGET_RC) CALL BUDGET_STORE_INIT(TBUDGETS(NBUDGET_RC), TRIM(HBUNAME), PRCS(:, :, :)*PRHODJ(:, :, :))
IF(BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT(TBUDGETS(NBUDGET_RI), TRIM(HBUNAME), PRIS(:, :, :)*PRHODJ(:, :, :))
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE QUANTITIES WITH THE GUESS OF THE FUTURE INSTANT
!               -------------------------------------------------------
!
!
!    beginning of the iterative loop (to compute the adjusted state)
!
DO JITER =1,ITERMAX
  !
  !*       2.3    compute the latent heat of vaporization Lv(T*) at t+1
  !                   and the latent heat of sublimation  Ls(T*) at t+1
  !
  DO JK=D%NKTB,D%NKTE
    DO JJ=D%NJB,D%NJE
      DO JI=D%NIB,D%NIE
        IF (JITER==1) ZT(JI,JJ,JK) = PTH(JI,JJ,JK) * PEXN(JI,JJ,JK)
        ZLV(JI,JJ,JK) = CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZT(JI,JJ,JK) -CST%XTT )
        ZLS(JI,JJ,JK) = CST%XLSTT + ( CST%XCPV - CST%XCI ) * ( ZT(JI,JJ,JK) -CST%XTT )
      ENDDO
    ENDDO
  ENDDO
  !
  !*       2.4   Iterate
  !
  IF (JITER==1) THEN
    ! compute with input values
    CALL ITERATION(PRV,PRC,PRI,ZRV,ZRC,ZRI)
  ELSE
    ! compute with updated values
    CALL ITERATION(ZRV,ZRC,ZRI,ZRV,ZRC,ZRI)
  ENDIF
ENDDO         ! end of the iterative loop
!
!*       5.     COMPUTE THE SOURCES AND STORES THE CLOUD FRACTION
!               -------------------------------------------------
!
!
DO JK=D%NKTB,D%NKTE
  DO JJ=D%NJB,D%NJE
    DO JI=D%NIB,D%NIE
      !
      !*       5.0    compute the variation of mixing ratio
      !
                                                           !         Rc - Rc*
      ZW1 = (ZRC(JI,JJ,JK) - PRC(JI,JJ,JK)) / PTSTEP       ! Pcon = ----------
                                                           !         2 Delta t
      ZW2 = (ZRI(JI,JJ,JK) - PRI(JI,JJ,JK)) / PTSTEP       ! idem ZW1 but for Ri
      !
      !*       5.1    compute the sources
      !
      IF( ZW1 < 0.0 ) THEN
        ZW1 = MAX ( ZW1, -PRCS(JI,JJ,JK) )
      ELSE
        ZW1 = MIN ( ZW1,  PRVS(JI,JJ,JK) )
      ENDIF
      PRVS(JI,JJ,JK) = PRVS(JI,JJ,JK) - ZW1
      PRCS(JI,JJ,JK) = PRCS(JI,JJ,JK) + ZW1
      PTHS(JI,JJ,JK) = PTHS(JI,JJ,JK) +        &
                      ZW1 * ZLV(JI,JJ,JK) / (ZCPH(JI,JJ,JK) * PEXNREF(JI,JJ,JK))
      !
      IF( ZW2 < 0.0 ) THEN
        ZW2 = MAX ( ZW2, -PRIS(JI,JJ,JK) )
      ELSE
        ZW2 = MIN ( ZW2,  PRVS(JI,JJ,JK) )
      ENDIF
      PRVS(JI,JJ,JK) = PRVS(JI,JJ,JK) - ZW2
      PRIS(JI,JJ,JK) = PRIS(JI,JJ,JK) + ZW2
      PTHS(JI,JJ,JK) = PTHS(JI,JJ,JK) +        &
                    ZW2 * ZLS(JI,JJ,JK) / (ZCPH(JI,JJ,JK) * PEXNREF(JI,JJ,JK))
    ENDDO
    !
    !*       5.2    compute the cloud fraction PCLDFR
    !
    IF ( .NOT. OSUBG_COND ) THEN
      DO JI=D%NIB,D%NIE
        IF (PRCS(JI,JJ,JK) + PRIS(JI,JJ,JK) > 1.E-12 / PTSTEP) THEN
          PCLDFR(JI,JJ,JK)  = 1.
        ELSE
          PCLDFR(JI,JJ,JK)  = 0.
        ENDIF
        IF ( SIZE(PSRCS,3) /= 0 ) THEN
          PSRCS(JI,JJ,JK) = PCLDFR(JI,JJ,JK)
        END IF
      ENDDO
    ELSE !OSUBG_COND case
      DO JI=D%NIB,D%NIE
        !We limit PRC_MF+PRI_MF to PRVS*PTSTEP to avoid negative humidity
        ZW1=PRC_MF(JI,JJ,JK)/PTSTEP
        ZW2=PRI_MF(JI,JJ,JK)/PTSTEP
        IF(ZW1+ZW2>PRVS(JI,JJ,JK)) THEN
          ZW1=ZW1*PRVS(JI,JJ,JK)/(ZW1+ZW2)
          ZW2=PRVS(JI,JJ,JK)-ZW1
        ENDIF
        PCLDFR(JI,JJ,JK)=MIN(1.,PCLDFR(JI,JJ,JK)+PCF_MF(JI,JJ,JK))
        PRCS(JI,JJ,JK)=PRCS(JI,JJ,JK)+ZW1
        PRIS(JI,JJ,JK)=PRIS(JI,JJ,JK)+ZW2
        PRVS(JI,JJ,JK)=PRVS(JI,JJ,JK)-(ZW1+ZW2)
        PTHS(JI,JJ,JK) = PTHS(JI,JJ,JK) + &
                      (ZW1 * ZLV(JI,JJ,JK) + ZW2 * ZLS(JI,JJ,JK)) / ZCPH(JI,JJ,JK) / PEXNREF(JI,JJ,JK)
        !
        IF(PRESENT(PHLC_HRC) .AND. PRESENT(PHLC_HCF)) THEN
          ZCRIAUT=ICEP%XCRIAUTC/PRHODREF(JI,JJ,JK)
          IF(HSUBG_MF_PDF=='NONE')THEN
            IF(ZW1*PTSTEP>PCF_MF(JI,JJ,JK) * ZCRIAUT) THEN
              PHLC_HRC(JI,JJ,JK)=PHLC_HRC(JI,JJ,JK)+ZW1*PTSTEP
              PHLC_HCF(JI,JJ,JK)=MIN(1.,PHLC_HCF(JI,JJ,JK)+PCF_MF(JI,JJ,JK))
            ENDIF
          ELSEIF(HSUBG_MF_PDF=='TRIANGLE')THEN
            !ZHCF is the precipitating part of the *cloud* and not of the grid cell
            IF(ZW1*PTSTEP>PCF_MF(JI,JJ,JK)*ZCRIAUT) THEN
              ZHCF=1.-.5*(ZCRIAUT*PCF_MF(JI,JJ,JK) / MAX(1.E-20, ZW1*PTSTEP))**2
              ZHR=ZW1*PTSTEP-(ZCRIAUT*PCF_MF(JI,JJ,JK))**3 / &
                                          &(3*MAX(1.E-20, ZW1*PTSTEP)**2)
            ELSEIF(2.*ZW1*PTSTEP<=PCF_MF(JI,JJ,JK) * ZCRIAUT) THEN
              ZHCF=0.
              ZHR=0.
            ELSE
              ZHCF=(2.*ZW1*PTSTEP-ZCRIAUT*PCF_MF(JI,JJ,JK))**2 / &
                         &(2.*MAX(1.E-20, ZW1*PTSTEP)**2)
              ZHR=(4.*(ZW1*PTSTEP)**3-3.*ZW1*PTSTEP*(ZCRIAUT*PCF_MF(JI,JJ,JK))**2+&
                          (ZCRIAUT*PCF_MF(JI,JJ,JK))**3) / &
                        &(3*MAX(1.E-20, ZW1*PTSTEP)**2)
            ENDIF
            ZHCF=ZHCF*PCF_MF(JI,JJ,JK) !to retrieve the part of the grid cell
            PHLC_HCF(JI,JJ,JK)=MIN(1.,PHLC_HCF(JI,JJ,JK)+ZHCF) !total part of the grid cell that is precipitating
            PHLC_HRC(JI,JJ,JK)=PHLC_HRC(JI,JJ,JK)+ZHR
          ENDIF
        ENDIF
        IF(PRESENT(PHLI_HRI) .AND. PRESENT(PHLI_HCF)) THEN
          ZCRIAUT=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(ZT(JI,JJ,JK)-CST%XTT)+ICEP%XBCRIAUTI))
          IF(HSUBG_MF_PDF=='NONE')THEN
            IF(ZW2*PTSTEP>PCF_MF(JI,JJ,JK) * ZCRIAUT) THEN
              PHLI_HRI(JI,JJ,JK)=PHLI_HRI(JI,JJ,JK)+ZW2*PTSTEP
              PHLI_HCF(JI,JJ,JK)=MIN(1.,PHLI_HCF(JI,JJ,JK)+PCF_MF(JI,JJ,JK))
            ENDIF
          ELSEIF(HSUBG_MF_PDF=='TRIANGLE')THEN
            !ZHCF is the precipitating part of the *cloud* and not of the grid cell
            IF(ZW2*PTSTEP>PCF_MF(JI,JJ,JK)*ZCRIAUT) THEN
              ZHCF=1.-.5*(ZCRIAUT*PCF_MF(JI,JJ,JK) / (ZW2*PTSTEP))**2
              ZHR=ZW2*PTSTEP-(ZCRIAUT*PCF_MF(JI,JJ,JK))**3/(3*(ZW2*PTSTEP)**2)
            ELSEIF(2.*ZW2*PTSTEP<=PCF_MF(JI,JJ,JK) * ZCRIAUT) THEN
              ZHCF=0.
              ZHR=0.
            ELSE
              ZHCF=(2.*ZW2*PTSTEP-ZCRIAUT*PCF_MF(JI,JJ,JK))**2 / (2.*(ZW2*PTSTEP)**2)
              ZHR=(4.*(ZW2*PTSTEP)**3-3.*ZW2*PTSTEP*(ZCRIAUT*PCF_MF(JI,JJ,JK))**2+&
                          (ZCRIAUT*PCF_MF(JI,JJ,JK))**3)/(3*(ZW2*PTSTEP)**2)
            ENDIF
            ZHCF=ZHCF*PCF_MF(JI,JJ,JK) !to retrieve the part of the grid cell
            PHLI_HCF(JI,JJ,JK)=MIN(1.,PHLI_HCF(JI,JJ,JK)+ZHCF) !total part of the grid cell that is precipitating
            PHLI_HRI(JI,JJ,JK)=PHLI_HRI(JI,JJ,JK)+ZHR
          ENDIF
        ENDIF
      ENDDO
      !
      IF(PRESENT(POUT_RV) .OR. PRESENT(POUT_RC) .OR. &
        &PRESENT(POUT_RI) .OR. PRESENT(POUT_TH)) THEN
        DO JI=D%NIB,D%NIE
          ZW1=PRC_MF(JI,JJ,JK)
          ZW2=PRI_MF(JI,JJ,JK)
          IF(ZW1+ZW2>ZRV(JI,JJ,JK)) THEN
            ZW1=ZW1*ZRV(JI,JJ,JK)/(ZW1+ZW2)
            ZW2=ZRV(JI,JJ,JK)-ZW1
          ENDIF
          ZRC(JI,JJ,JK)=ZRC(JI,JJ,JK)+ZW1
          ZRI(JI,JJ,JK)=ZRI(JI,JJ,JK)+ZW2
          ZRV(JI,JJ,JK)=ZRV(JI,JJ,JK)-(ZW1+ZW2)
          ZT(JI,JJ,JK) = ZT(JI,JJ,JK) + &
                      (ZW1 * ZLV(JI,JJ,JK) + ZW2 * ZLS(JI,JJ,JK)) / ZCPH(JI,JJ,JK)
        ENDDO
      ENDIF
    ENDIF !OSUBG_COND
  ENDDO
ENDDO
!
IF(PRESENT(POUT_RV)) POUT_RV=ZRV
IF(PRESENT(POUT_RC)) POUT_RC=ZRC
IF(PRESENT(POUT_RI)) POUT_RI=ZRI
IF(PRESENT(POUT_TH)) POUT_TH=ZT / PEXN(:,:,:)
!
!
!*       6.  STORE THE BUDGET TERMS
!            ----------------------
!
IF(BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END(TBUDGETS(NBUDGET_TH), TRIM(HBUNAME), PTHS(:, :, :)*PRHODJ(:, :, :))
IF(BUCONF%LBUDGET_RV) CALL BUDGET_STORE_END(TBUDGETS(NBUDGET_RV), TRIM(HBUNAME), PRVS(:, :, :)*PRHODJ(:, :, :))
IF(BUCONF%LBUDGET_RC) CALL BUDGET_STORE_END(TBUDGETS(NBUDGET_RC), TRIM(HBUNAME), PRCS(:, :, :)*PRHODJ(:, :, :))
IF(BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END(TBUDGETS(NBUDGET_RI), TRIM(HBUNAME), PRIS(:, :, :)*PRHODJ(:, :, :))
!------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('ICE_ADJUST',1,ZHOOK_HANDLE)
!
CONTAINS
SUBROUTINE ITERATION(PRV_IN,PRC_IN,PRI_IN,PRV_OUT,PRC_OUT,PRI_OUT)

REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRV_IN ! Water vapor m.r. to adjust in input
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRC_IN ! Cloud water m.r. to adjust in input
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN) :: PRI_IN ! Cloud ice   m.r. to adjust in input
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PRV_OUT ! Water vapor m.r. to adjust in output
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PRC_OUT ! Cloud water m.r. to adjust in output
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PRI_OUT ! Cloud ice   m.r. to adjust in output
!
!*       2.4    compute the specific heat for moist air (Cph) at t+1
DO JK=D%NKTB,D%NKTE
  DO JJ=D%NJB,D%NJE
    DO JI=D%NIB,D%NIE
      SELECT CASE(KRR)
        CASE(7)
          ZCPH(JI,JJ,JK) = CST%XCPD + CST%XCPV * PRV_IN(JI,JJ,JK)                             &
                                    + CST%XCL  * (PRC_IN(JI,JJ,JK) + PRR(JI,JJ,JK))             &
                                    + CST%XCI  * (PRI_IN(JI,JJ,JK) + PRS(JI,JJ,JK) + PRG(JI,JJ,JK) + PRH(JI,JJ,JK))
        CASE(6)
          ZCPH(JI,JJ,JK) = CST%XCPD + CST%XCPV * PRV_IN(JI,JJ,JK)                             &
                                    + CST%XCL  * (PRC_IN(JI,JJ,JK) + PRR(JI,JJ,JK))             &
                                    + CST%XCI  * (PRI_IN(JI,JJ,JK) + PRS(JI,JJ,JK) + PRG(JI,JJ,JK))
        CASE(5)
          ZCPH(JI,JJ,JK) = CST%XCPD + CST%XCPV * PRV_IN(JI,JJ,JK)                             &
                                    + CST%XCL  * (PRC_IN(JI,JJ,JK) + PRR(JI,JJ,JK))             &
                                    + CST%XCI  * (PRI_IN(JI,JJ,JK) + PRS(JI,JJ,JK))
        CASE(3)
          ZCPH(JI,JJ,JK) = CST%XCPD + CST%XCPV * PRV_IN(JI,JJ,JK)               &
                                    + CST%XCL  * (PRC_IN(JI,JJ,JK) + PRR(JI,JJ,JK))
        CASE(2)
          ZCPH(JI,JJ,JK) = CST%XCPD + CST%XCPV * PRV_IN(JI,JJ,JK) &
                                    + CST%XCL  * PRC_IN(JI,JJ,JK)
      END SELECT
    ENDDO
  ENDDO
ENDDO
!
IF ( OSUBG_COND ) THEN
  !
  !*       3.     SUBGRID CONDENSATION SCHEME
  !               ---------------------------
  !
  !   PSRC= s'rci'/Sigma_s^2
  !   ZT is INOUT
  CALL CONDENSATION(D, CST, ICEP, NEB, &
       HFRAC_ICE, HCONDENS, HLAMBDA3,                                    &
       PPABST, PZZ, PRHODREF, ZT, PRV_IN, PRV_OUT, PRC_IN, PRC_OUT, PRI_IN, PRI_OUT, &
       PRS, PRG, PSIGS, LMFCONV, PMFCONV, PCLDFR, &
       PSRCS, .TRUE., OSIGMAS,                                           &
       OCND2, PSIGQSAT,                                                  &
       PLV=ZLV, PLS=ZLS, PCPH=ZCPH,                                      &
       PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF, PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF,&
       PICE_CLD_WGT=PICE_CLD_WGT)
ELSE
  !
  !*       4.     ALL OR NOTHING CONDENSATION SCHEME
  !                            FOR MIXED-PHASE CLOUD
  !               -----------------------------------------------
  !
  ZSIGS(:,:,:)=0.
  ZSIGQSAT(:,:)=0.
  !We use ZSRCS because in Méso-NH, PSRCS can be a zero-length array in this case
  !ZT is INOUT
  CALL CONDENSATION(D, CST, ICEP, NEB, &
       HFRAC_ICE, HCONDENS, HLAMBDA3,                                    &
       PPABST, PZZ, PRHODREF, ZT, PRV_IN, PRV_OUT, PRC_IN, PRC_OUT, PRI_IN, PRI_OUT, &
       PRS, PRG, ZSIGS, LMFCONV, PMFCONV, PCLDFR, &
       ZSRCS, .TRUE., OSIGMAS=.TRUE.,                                    &
       OCND2=OCND2, PSIGQSAT=ZSIGQSAT,                                   &
       PLV=ZLV, PLS=ZLS, PCPH=ZCPH,                                      &
       PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF, PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF,&
       PICE_CLD_WGT=PICE_CLD_WGT)
ENDIF

END SUBROUTINE ITERATION

END SUBROUTINE ICE_ADJUST 
