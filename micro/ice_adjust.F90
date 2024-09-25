!MNH_LIC Copyright 1996-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################################################################
      SUBROUTINE ICE_ADJUST (D, CST, ICEP, NEBN, TURBN, PARAMI, BUCONF, KRR,   &
                            &HBUNAME,                                          &
                            &PTSTEP, PSIGQSAT,                                 &
                            &PRHODJ, PEXNREF, PRHODREF, PSIGS, LMFCONV, PMFCONV,&
                            &PPABST, PZZ,                                      &
                            &PEXN, PCF_MF, PRC_MF, PRI_MF,                     &
                            &PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR,             &
                            &PRV, PRC, PRVS, PRCS, PTH, PTHS,                  &
                            &OCOMPUTE_SRC, PSRCS, PCLDFR,                      &
                            &PRR, PRI, PRIS, PRS, PRG, TBUDGETS, KBUDGETS,     &
                            &PICE_CLD_WGT,                                     &
                            &PRH,                                              &
                            &POUT_RV, POUT_RC, POUT_RI, POUT_TH,               &
                            &PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF)

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
!!      2018-02 K.I.Ivarsson : More outputs for OCND2 option
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!!      2020-12 U. Andrae : Introduce SPP for HARMONIE-AROME
!!     R. El Khatib 24-Aug-2021 Optimizations
!!     R. El Khatib 24-Oct-2023 Re-vectorize ;-)
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_CST,        ONLY: CST_t
USE MODD_NEB_n,      ONLY: NEB_t
USE MODD_TURB_n,         ONLY: TURB_t
USE MODD_PARAM_ICE_n,    ONLY: PARAM_ICE_t
USE MODD_BUDGET,     ONLY: TBUDGETDATA, TBUDGETCONF_t, NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RI
USE MODD_RAIN_ICE_PARAM_n, ONLY : RAIN_ICE_PARAM_t
!
USE MODE_BUDGET_PHY,         ONLY: BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY
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
TYPE(NEB_t),              INTENT(IN)    :: NEBN
TYPE(TURB_t),             INTENT(IN)    :: TURBN
TYPE(PARAM_ICE_t),        INTENT(IN)    :: PARAMI
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER(LEN=4),         INTENT(IN)    :: HBUNAME  ! Name of the budget
REAL,                     INTENT(IN)   :: PTSTEP    ! Double Time step
                                                    ! (single if cold start)
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PSIGQSAT  ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    ::  PEXNREF ! Reference Exner function
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    ::  PRHODREF
!
REAL, DIMENSION(MERGE(D%NIJT,0,NEBN%LSUBG_COND),&
                MERGE(D%NKT,0,NEBN%LSUBG_COND)),           INTENT(IN)    ::  PSIGS   ! Sigma_s at time t
LOGICAL,                                              INTENT(IN)    ::  LMFCONV ! =SIZE(PMFCONV)!=0
REAL, DIMENSION(MERGE(D%NIJT,0,LMFCONV),&
                MERGE(D%NKT,0,LMFCONV)),              INTENT(IN)   ::  PMFCONV ! convective mass flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    ::  PPABST  ! Absolute Pressure at t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    ::  PZZ     ! height of model layer
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    ::  PEXN    ! Exner function
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PCF_MF   ! Convective Mass Flux Cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRC_MF   ! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRI_MF   ! Convective Mass Flux ice mixing ratio
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRV     ! Water vapor m.r. to adjust
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRC     ! Cloud water m.r. to adjust
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PTH     ! Theta to adjust
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PTHS    ! Theta source
LOGICAL,                            INTENT(IN)    :: OCOMPUTE_SRC
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                                       ! s'rc'/2Sigma_s2 at time t+1
                                                                       ! multiplied by Lambda_3
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PCLDFR  ! Cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PICLDFR ! ice cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PWCLDFR ! water or mixed-phase cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PSSIO   ! Super-saturation with respect to ice in the  
                                                        ! supersaturated fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PSSIU   ! Sub-saturation with respect to ice in the  
                                                        ! subsaturated fraction 
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PIFR    ! Ratio cloud ice moist part to dry part
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT)::  PRIS ! Cloud ice  m.r. at t+1
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRR  ! Rain water m.r. to adjust
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRI  ! Cloud ice  m.r. to adjust
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRS  ! Aggregate  m.r. to adjust
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PRG  ! Graupel    m.r. to adjust
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS),       INTENT(INOUT)::  TBUDGETS
INTEGER,                                      INTENT(IN)   ::  KBUDGETS
REAL, DIMENSION(D%NIJT),       OPTIONAL, INTENT(IN)   ::  PICE_CLD_WGT
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)   ::  PRH  ! Hail       m.r. to adjust
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  POUT_RV ! Adjusted value
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  POUT_RC ! Adjusted value
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  POUT_RI ! Adjusted value
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  POUT_TH ! Adjusted value
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  PHLC_HRC
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  PHLC_HCF
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  PHLI_HRI
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)  ::  PHLI_HCF
!
!
!*       0.2   Declarations of local variables :
!
!
REAL  :: ZW1,ZW2    ! intermediate fields
REAL, DIMENSION(D%NIJT,D%NKT) &
                         :: ZT,   &  ! adjusted temperature
                   ZRV, ZRC, ZRI, &  ! adjusted state
                            ZCPH, &  ! guess of the CPh for the mixing
                            ZLV,  &  ! guess of the Lv at t+1
                            ZLS      ! guess of the Ls at t+1
REAL :: ZCRIAUT, & ! Autoconversion thresholds
        ZHCF, ZHR
!
INTEGER             :: JITER,ITERMAX ! iterative loop for first order adjustment
INTEGER             :: JIJ, JK
INTEGER :: IKTB, IKTE, IIJB, IIJE
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZSIGS, ZSRCS
REAL, DIMENSION(D%NIJT) :: ZSIGQSAT
LOGICAL :: LLNONE, LLTRIANGLE, LLHLC_H, LLHLI_H

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
IF (LHOOK) CALL DR_HOOK('ICE_ADJUST',0,ZHOOK_HANDLE)
!
IKTB=D%NKTB
IKTE=D%NKTE
IIJB=D%NIJB
IIJE=D%NIJE
!
ITERMAX=1
!
IF(BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), TRIM(HBUNAME), PTHS(:, :)*PRHODJ(:, :))
IF(BUCONF%LBUDGET_RV) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), TRIM(HBUNAME), PRVS(:, :)*PRHODJ(:, :))
IF(BUCONF%LBUDGET_RC) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), TRIM(HBUNAME), PRCS(:, :)*PRHODJ(:, :))
IF(BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), TRIM(HBUNAME), PRIS(:, :)*PRHODJ(:, :))
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
  DO JK=IKTB,IKTE
    DO JIJ=IIJB,IIJE
      IF (JITER==1) ZT(JIJ,JK) = PTH(JIJ,JK) * PEXN(JIJ,JK)
      ZLV(JIJ,JK) = CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZT(JIJ,JK) -CST%XTT )
      ZLS(JIJ,JK) = CST%XLSTT + ( CST%XCPV - CST%XCI ) * ( ZT(JIJ,JK) -CST%XTT )
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
DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE
    !
    !*       5.0    compute the variation of mixing ratio
    !
                                                         !         Rc - Rc*
    ZW1 = (ZRC(JIJ,JK) - PRC(JIJ,JK)) / PTSTEP       ! Pcon = ----------
                                                         !         2 Delta t
    ZW2 = (ZRI(JIJ,JK) - PRI(JIJ,JK)) / PTSTEP       ! idem ZW1 but for Ri
    !
    !*       5.1    compute the sources
    !
    IF( ZW1 < 0.0 ) THEN
      ZW1 = MAX ( ZW1, -PRCS(JIJ,JK) )
    ELSE
      ZW1 = MIN ( ZW1,  PRVS(JIJ,JK) )
    ENDIF
    PRVS(JIJ,JK) = PRVS(JIJ,JK) - ZW1
    PRCS(JIJ,JK) = PRCS(JIJ,JK) + ZW1
    PTHS(JIJ,JK) = PTHS(JIJ,JK) +        &
                    ZW1 * ZLV(JIJ,JK) / (ZCPH(JIJ,JK) * PEXNREF(JIJ,JK))
    !
    IF( ZW2 < 0.0 ) THEN
      ZW2 = MAX ( ZW2, -PRIS(JIJ,JK) )
    ELSE
      ZW2 = MIN ( ZW2,  PRVS(JIJ,JK) )
    ENDIF
    PRVS(JIJ,JK) = PRVS(JIJ,JK) - ZW2
    PRIS(JIJ,JK) = PRIS(JIJ,JK) + ZW2
    PTHS(JIJ,JK) = PTHS(JIJ,JK) +        &
                  ZW2 * ZLS(JIJ,JK) / (ZCPH(JIJ,JK) * PEXNREF(JIJ,JK))
  ENDDO
  !
  !*       5.2    compute the cloud fraction PCLDFR
  !
  IF ( .NOT. NEBN%LSUBG_COND ) THEN
    DO JIJ=IIJB,IIJE
      IF (PRCS(JIJ,JK) + PRIS(JIJ,JK) > 1.E-12 / PTSTEP) THEN
        PCLDFR(JIJ,JK)  = 1.
      ELSE
        PCLDFR(JIJ,JK)  = 0.
      ENDIF
      IF (OCOMPUTE_SRC) THEN
        PSRCS(JIJ,JK) = PCLDFR(JIJ,JK)
      END IF
    ENDDO
  ELSE !NEBN%LSUBG_COND case
    ! Tests on characters strings can break the vectorization, or at least they would
    ! slow down considerably the performance of a vector loop. One should use tests on
    ! reals, integers or booleans only. REK.
    LLNONE=PARAMI%CSUBG_MF_PDF=='NONE'
    LLTRIANGLE=PARAMI%CSUBG_MF_PDF=='TRIANGLE'
    LLHLC_H=PRESENT(PHLC_HRC).AND.PRESENT(PHLC_HCF)
    LLHLI_H=PRESENT(PHLI_HRI).AND.PRESENT(PHLI_HCF)
    DO JIJ=IIJB,IIJE
      !We limit PRC_MF+PRI_MF to PRVS*PTSTEP to avoid negative humidity
      ZW1=PRC_MF(JIJ,JK)/PTSTEP
      ZW2=PRI_MF(JIJ,JK)/PTSTEP
      IF(ZW1+ZW2>PRVS(JIJ,JK)) THEN
        ZW1=ZW1*PRVS(JIJ,JK)/(ZW1+ZW2)
        ZW2=PRVS(JIJ,JK)-ZW1
      ENDIF
      PCLDFR(JIJ,JK)=MIN(1.,PCLDFR(JIJ,JK)+PCF_MF(JIJ,JK))
      PRCS(JIJ,JK)=PRCS(JIJ,JK)+ZW1
      PRIS(JIJ,JK)=PRIS(JIJ,JK)+ZW2
      PRVS(JIJ,JK)=PRVS(JIJ,JK)-(ZW1+ZW2)
      PTHS(JIJ,JK) = PTHS(JIJ,JK) + &
                    (ZW1 * ZLV(JIJ,JK) + ZW2 * ZLS(JIJ,JK)) / ZCPH(JIJ,JK) / PEXNREF(JIJ,JK)
      !
      IF(LLHLC_H) THEN
        ZCRIAUT=ICEP%XCRIAUTC/PRHODREF(JIJ,JK)
        IF(LLNONE)THEN
          IF(ZW1*PTSTEP>PCF_MF(JIJ,JK) * ZCRIAUT) THEN
            PHLC_HRC(JIJ,JK)=PHLC_HRC(JIJ,JK)+ZW1*PTSTEP
            PHLC_HCF(JIJ,JK)=MIN(1.,PHLC_HCF(JIJ,JK)+PCF_MF(JIJ,JK))
          ENDIF
        ELSEIF(LLTRIANGLE)THEN
          !ZHCF is the precipitating part of the *cloud* and not of the grid cell
          IF(ZW1*PTSTEP>PCF_MF(JIJ,JK)*ZCRIAUT) THEN
            ZHCF=1.-.5*(ZCRIAUT*PCF_MF(JIJ,JK) / MAX(1.E-20, ZW1*PTSTEP))**2
            ZHR=ZW1*PTSTEP-(ZCRIAUT*PCF_MF(JIJ,JK))**3 / &
                                        &(3*MAX(1.E-20, ZW1*PTSTEP)**2)
          ELSEIF(2.*ZW1*PTSTEP<=PCF_MF(JIJ,JK) * ZCRIAUT) THEN
            ZHCF=0.
            ZHR=0.
          ELSE
            ZHCF=(2.*ZW1*PTSTEP-ZCRIAUT*PCF_MF(JIJ,JK))**2 / &
                       &(2.*MAX(1.E-20, ZW1*PTSTEP)**2)
            ZHR=(4.*(ZW1*PTSTEP)**3-3.*ZW1*PTSTEP*(ZCRIAUT*PCF_MF(JIJ,JK))**2+&
                        (ZCRIAUT*PCF_MF(JIJ,JK))**3) / &
                      &(3*MAX(1.E-20, ZW1*PTSTEP)**2)
          ENDIF
          ZHCF=ZHCF*PCF_MF(JIJ,JK) !to retrieve the part of the grid cell
          PHLC_HCF(JIJ,JK)=MIN(1.,PHLC_HCF(JIJ,JK)+ZHCF) !total part of the grid cell that is precipitating
          PHLC_HRC(JIJ,JK)=PHLC_HRC(JIJ,JK)+ZHR
        ENDIF
      ENDIF
      IF(LLHLI_H) THEN
        ZCRIAUT=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(ZT(JIJ,JK)-CST%XTT)+ICEP%XBCRIAUTI))
        IF(LLNONE)THEN
          IF(ZW2*PTSTEP>PCF_MF(JIJ,JK) * ZCRIAUT) THEN
            PHLI_HRI(JIJ,JK)=PHLI_HRI(JIJ,JK)+ZW2*PTSTEP
            PHLI_HCF(JIJ,JK)=MIN(1.,PHLI_HCF(JIJ,JK)+PCF_MF(JIJ,JK))
          ENDIF
        ELSEIF(LLTRIANGLE)THEN
          !ZHCF is the precipitating part of the *cloud* and not of the grid cell
          IF(ZW2*PTSTEP>PCF_MF(JIJ,JK)*ZCRIAUT) THEN
            ZHCF=1.-.5*(ZCRIAUT*PCF_MF(JIJ,JK) / (ZW2*PTSTEP))**2
            ZHR=ZW2*PTSTEP-(ZCRIAUT*PCF_MF(JIJ,JK))**3/(3*(ZW2*PTSTEP)**2)
          ELSEIF(2.*ZW2*PTSTEP<=PCF_MF(JIJ,JK) * ZCRIAUT) THEN
            ZHCF=0.
            ZHR=0.
          ELSE
            ZHCF=(2.*ZW2*PTSTEP-ZCRIAUT*PCF_MF(JIJ,JK))**2 / (2.*(ZW2*PTSTEP)**2)
            ZHR=(4.*(ZW2*PTSTEP)**3-3.*ZW2*PTSTEP*(ZCRIAUT*PCF_MF(JIJ,JK))**2+&
                        (ZCRIAUT*PCF_MF(JIJ,JK))**3)/(3*(ZW2*PTSTEP)**2)
          ENDIF
          ZHCF=ZHCF*PCF_MF(JIJ,JK) !to retrieve the part of the grid cell
          PHLI_HCF(JIJ,JK)=MIN(1.,PHLI_HCF(JIJ,JK)+ZHCF) !total part of the grid cell that is precipitating
          PHLI_HRI(JIJ,JK)=PHLI_HRI(JIJ,JK)+ZHR
        ENDIF
      ENDIF
    ENDDO
    !
    IF(PRESENT(POUT_RV) .OR. PRESENT(POUT_RC) .OR. &
      &PRESENT(POUT_RI) .OR. PRESENT(POUT_TH)) THEN
      DO JIJ=IIJB,IIJE
        ZW1=PRC_MF(JIJ,JK)
        ZW2=PRI_MF(JIJ,JK)
        IF(ZW1+ZW2>ZRV(JIJ,JK)) THEN
          ZW1=ZW1*ZRV(JIJ,JK)/(ZW1+ZW2)
          ZW2=ZRV(JIJ,JK)-ZW1
        ENDIF
        ZRC(JIJ,JK)=ZRC(JIJ,JK)+ZW1
        ZRI(JIJ,JK)=ZRI(JIJ,JK)+ZW2
        ZRV(JIJ,JK)=ZRV(JIJ,JK)-(ZW1+ZW2)
        ZT(JIJ,JK) = ZT(JIJ,JK) + &
                    (ZW1 * ZLV(JIJ,JK) + ZW2 * ZLS(JIJ,JK)) / ZCPH(JIJ,JK)
      ENDDO
    ENDIF
  ENDIF !NEBN%LSUBG_COND
ENDDO
!
IF(PRESENT(POUT_RV)) POUT_RV=ZRV
IF(PRESENT(POUT_RC)) POUT_RC=ZRC
IF(PRESENT(POUT_RI)) POUT_RI=ZRI
IF(PRESENT(POUT_TH)) POUT_TH=ZT / PEXN(:,:)
!
!
!*       6.  STORE THE BUDGET TERMS
!            ----------------------
!
IF(BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), TRIM(HBUNAME), PTHS(:, :)*PRHODJ(:, :))
IF(BUCONF%LBUDGET_RV) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), TRIM(HBUNAME), PRVS(:, :)*PRHODJ(:, :))
IF(BUCONF%LBUDGET_RC) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), TRIM(HBUNAME), PRCS(:, :)*PRHODJ(:, :))
IF(BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), TRIM(HBUNAME), PRIS(:, :)*PRHODJ(:, :))
!------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('ICE_ADJUST',1,ZHOOK_HANDLE)
!
CONTAINS
SUBROUTINE ITERATION(PRV_IN,PRC_IN,PRI_IN,PRV_OUT,PRC_OUT,PRI_OUT)

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRV_IN ! Water vapor m.r. to adjust in input
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRC_IN ! Cloud water m.r. to adjust in input
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRI_IN ! Cloud ice   m.r. to adjust in input
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PRV_OUT ! Water vapor m.r. to adjust in output
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PRC_OUT ! Cloud water m.r. to adjust in output
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PRI_OUT ! Cloud ice   m.r. to adjust in output
!
!*       2.4    compute the specific heat for moist air (Cph) at t+1
DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE
    SELECT CASE(KRR)
      CASE(7)
        ZCPH(JIJ,JK) = CST%XCPD + CST%XCPV * PRV_IN(JIJ,JK)                             &
                                + CST%XCL  * (PRC_IN(JIJ,JK) + PRR(JIJ,JK))             &
                                + CST%XCI  * (PRI_IN(JIJ,JK) + PRS(JIJ,JK) + PRG(JIJ,JK) + PRH(JIJ,JK))
      CASE(6)
        ZCPH(JIJ,JK) = CST%XCPD + CST%XCPV * PRV_IN(JIJ,JK)                             &
                                + CST%XCL  * (PRC_IN(JIJ,JK) + PRR(JIJ,JK))             &
                                + CST%XCI  * (PRI_IN(JIJ,JK) + PRS(JIJ,JK) + PRG(JIJ,JK))
      CASE(5)
        ZCPH(JIJ,JK) = CST%XCPD + CST%XCPV * PRV_IN(JIJ,JK)                             &
                                + CST%XCL  * (PRC_IN(JIJ,JK) + PRR(JIJ,JK))             &
                                + CST%XCI  * (PRI_IN(JIJ,JK) + PRS(JIJ,JK))
      CASE(3)
        ZCPH(JIJ,JK) = CST%XCPD + CST%XCPV * PRV_IN(JIJ,JK)               &
                                + CST%XCL  * (PRC_IN(JIJ,JK) + PRR(JIJ,JK))
      CASE(2)
        ZCPH(JIJ,JK) = CST%XCPD + CST%XCPV * PRV_IN(JIJ,JK) &
                                + CST%XCL  * PRC_IN(JIJ,JK)
    END SELECT
  ENDDO
ENDDO
!
IF ( NEBN%LSUBG_COND ) THEN
  !
  !*       3.     SUBGRID CONDENSATION SCHEME
  !               ---------------------------
  !
  !   PSRC= s'rci'/Sigma_s^2
  !   ZT is INOUT
  CALL CONDENSATION(D, CST, ICEP, NEBN, TURBN, &
       NEBN%CFRAC_ICE_ADJUST,NEBN%CCONDENS, NEBN%CLAMBDA3,                             &
       PPABST, PZZ, PRHODREF, ZT, PRV_IN, PRV_OUT, PRC_IN, PRC_OUT, PRI_IN, PRI_OUT, &
       PRR, PRS, PRG, PSIGS, LMFCONV, PMFCONV, PCLDFR, &
       PSRCS, .TRUE., NEBN%LSIGMAS, PARAMI%LOCND2,                       &
       PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR, PSIGQSAT,                   &
       PLV=ZLV, PLS=ZLS, PCPH=ZCPH,                                      &
       PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF, PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF,&
       PICE_CLD_WGT=PICE_CLD_WGT)
ELSE
  !
  !*       4.     ALL OR NOTHING CONDENSATION SCHEME
  !                            FOR MIXED-PHASE CLOUD
  !               -----------------------------------------------
  !
  ZSIGS(:,:)=0.
  ZSIGQSAT(:)=0.
  !We use ZSRCS because in Méso-NH, PSRCS can be a zero-length array in this case
  !ZT is INOUT
  CALL CONDENSATION(D, CST, ICEP, NEBN, TURBN, &
       NEBN%CFRAC_ICE_ADJUST,NEBN%CCONDENS, NEBN%CLAMBDA3,                             &
       PPABST, PZZ, PRHODREF, ZT, PRV_IN, PRV_OUT, PRC_IN, PRC_OUT, PRI_IN, PRI_OUT, &
       PRR, PRS, PRG, ZSIGS, LMFCONV, PMFCONV, PCLDFR, &
       ZSRCS, .TRUE., .TRUE., PARAMI%LOCND2,                             &
       PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR, ZSIGQSAT,                   &
       PLV=ZLV, PLS=ZLS, PCPH=ZCPH,                                      &
       PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF, PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF,&
       PICE_CLD_WGT=PICE_CLD_WGT)
ENDIF

END SUBROUTINE ITERATION

END SUBROUTINE ICE_ADJUST 
