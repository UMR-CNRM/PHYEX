!     ##########################################################################
      SUBROUTINE ICE_ADJUST (KKA, KKU, KKL, KRR, HFRAC_ICE,                    &
                             HBUNAME, OSUBG_COND, OSIGMAS, OCND2,              &
                             PTSTEP, PSIGQSAT,                                 &
                             PRHODJ, PEXNREF, PSIGS, PMFCONV, PPABST, PZZ,     &
                             PEXN, PCF_MF, PRC_MF, PRI_MF,                     &
                             PRV, PRC, PRVS, PRCS, PTH, PTHS, PSRCS, PCLDFR,   &
                             PRR, PRI, PRIS, PRS, PRG,                         &
                             PRH, POUT_RV, POUT_RC, POUT_RI, POUT_TH,          &
                             YSPP_PSIGQSAT, YSPP_ICE_CLD_WGT,                  &
                             YDDDH, YDLDDH, YDMDDH)
      USE PARKIND1, ONLY : JPRB
      USE MODD_SPP_TYPE
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      USE DDH_MIX , ONLY : TYP_DDH
      USE YOMLDDH , ONLY : TLDDH
      USE YOMMDDH , ONLY : TMDDH
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
!!         NBUPROCCTR 
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
!!      2016-07 S. Riette: adjustement is now realized on state variables (PRV, PRC, PRI, PTH)
!!                         whereas tendencies are still applied on S variables.
!!                         This modification allows to call ice_adjust on T variable
!!                         or to call it on S variables
!!      2016-11 S. Riette: all-or-nothing adjustment now uses condensation
!!      2020-12 U. Andrae : Introduce SPP for HARMONIE-AROME
!!     R. El Khatib 24-Aug-2021 Optimizations
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_CONF
USE MODD_BUDGET
!
USE MODI_CONDENSATION
USE MODI_BUDGET
USE MODE_FMWRIT
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN)    :: KKA  !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU  !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL  !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
CHARACTER*1,              INTENT(IN)    :: HFRAC_ICE
CHARACTER*4,              INTENT(IN)    :: HBUNAME  ! Name of the budget
LOGICAL,                  INTENT(IN)    :: OSUBG_COND ! Switch for Subgrid 
                                                    ! Condensation
LOGICAL                                 :: OSIGMAS  ! Switch for Sigma_s: 
                                                    ! use values computed in CONDENSATION
                                                    ! or that from turbulence scheme
LOGICAL                                 :: OCND2    ! logical switch to sparate liquid 
                                                    ! and ice
                                                    ! more rigid (DEFALT value : .FALSE.)
REAL,                     INTENT(IN)   :: PTSTEP    ! Double Time step
                                                    ! (single if cold start)
REAL,                     INTENT(IN)   :: PSIGQSAT  ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PEXNREF ! Reference Exner function
!
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PSIGS   ! Sigma_s at time t
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PMFCONV ! convective mass flux
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PPABST  ! Absolute Pressure at t        
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PZZ     ! height of model layer
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PEXN    ! Exner function
!
REAL, DIMENSION(:,:,:), CONTIGUOUS,     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
REAL, DIMENSION(:,:,:), CONTIGUOUS,     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:), CONTIGUOUS,     INTENT(IN)    :: PRI_MF! Convective Mass Flux ice mixing ratio
!
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PRV     ! Water vapor m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)   ::  PRC     ! Cloud water m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(IN)    :: PTH     ! Theta to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                   ! s'rc'/2Sigma_s2 at time t+1
                                                   ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:), CONTIGUOUS,   INTENT(OUT)   :: PCLDFR  ! Cloud fraction          
!
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(INOUT)::  PRIS ! Cloud ice  m.r. at t+1
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(IN)   ::  PRR  ! Rain water m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(IN)   ::  PRI  ! Cloud ice  m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(IN)   ::  PRS  ! Aggregate  m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS,  INTENT(IN)   ::  PRG  ! Graupel    m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(IN)   ::  PRH  ! Hail       m.r. to adjust
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  POUT_RV ! Adjusted value
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  POUT_RC ! Adjusted value
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  POUT_RI ! Adjusted value
REAL, DIMENSION(:,:,:), CONTIGUOUS, OPTIONAL, INTENT(OUT)  ::  POUT_TH ! Adjusted value
TYPE(TSPP_CONFIG_MPA),    INTENT(IN)   :: YSPP_PSIGQSAT,YSPP_ICE_CLD_WGT
TYPE(TYP_DDH)         , INTENT(INOUT) ::  YDDDH
TYPE(TLDDH)         , INTENT(IN) ::  YDLDDH
TYPE(TMDDH)         , INTENT(IN) ::  YDMDDH
!
!*       0.2   Declarations of local variables :
!
!
REAL  :: ZEPS  ! Mv/Md
REAL  :: ZT00,ZT0   ! Min and max temperature for the mixed phase liquid and solid water
                    ! for the coeff CND of the barycentric mixing ratio
REAL  :: ZW1,ZW2    ! intermediate fields
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) &
                         :: ZT,   &  ! adjusted temperature
                   ZRV, ZRC, ZRI, &  ! adjusted state
                            ZCPH, &  ! guess of the CPh for the mixing
                            ZLV,  &  ! guess of the Lv at t+1
                            ZLS,  &  ! guess of the Ls at t+1
                            ZH2O     ! Fraction of all types of water (kg/kg) (for OCND2)
!
INTEGER             :: IIU,IJU,IKU! dimensions of dummy arrays
INTEGER             :: IIB,IJB    ! Horz index values of the first inner mass points
INTEGER             :: IIE,IJE    ! Horz index values of the last inner mass points
INTEGER             :: IKB        ! K index value of the first inner mass point
INTEGER             :: IKE        ! K index value of the last inner mass point
INTEGER             :: JITER,ITERMAX ! iterative loop for first order adjustment
INTEGER             :: JI, JJ, JK
!
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) :: ZSIGS
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
IF (LHOOK) CALL DR_HOOK('ICE_ADJUST',0,ZHOOK_HANDLE)
IIU = SIZE(PEXNREF,1)
IJU = SIZE(PEXNREF,2)
IKU = SIZE(PEXNREF,3)
IIB = 1 + JPHEXT
IIE = IIU - JPHEXT
IJB = 1 + JPHEXT
IJE = IJU - JPHEXT
IKB=KKA+JPVEXT*KKL
IKE=KKU-JPVEXT*KKL
!
ITERMAX=1
!
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
  DO JK=1,IKU
    DO JJ=1,IJU
      DO JI=1,IIU
        IF (JITER==1) ZT(JI,JJ,JK) = PTH(JI,JJ,JK) * PEXN(JI,JJ,JK)
        ZLV(JI,JJ,JK) = XLVTT + ( XCPV - XCL ) * ( ZT(JI,JJ,JK) -XTT )
        ZLS(JI,JJ,JK) = XLSTT + ( XCPV - XCI ) * ( ZT(JI,JJ,JK) -XTT )
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
!       ZW8 = 4.*ZW2 ! + 0.5_JPRB*PRSS*PTSTEP + 0.25_JPRB*PRGS*PTSTEP ! include snow and graupel to cloudcover
!       ZW8(:,:,1) =  ZW2(:,:,1)                       ! but don't do this for the lowest model
!       ZW8(:,:,2) =  ZW2(:,:,2)                       ! level (Hail not concidered yet)
!
!*       5.0    compute the variation of mixing ratio
!
DO JK=1,IKU
  DO JJ=1,IJU
    DO JI=1,IIU
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
!
!*       5.2    compute the cloud fraction PCLDFR
!
    IF ( .NOT. OSUBG_COND ) THEN
      DO JI=1,IIU
        IF (PRCS(JI,JJ,JK) + PRIS(JI,JJ,JK) > 1.E-12 / PTSTEP) THEN
          PCLDFR(JI,JJ,JK)  = 1.
        ELSE
          PCLDFR(JI,JJ,JK)  = 0.
        ENDIF
        IF ( SIZE(PSRCS,3) /= 0 ) THEN
          PSRCS(JI,JJ,JK) = PCLDFR(JI,JJ,JK)
        END IF
      ENDDO
    ELSE
      DO JI=1,IIU
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
      ENDDO
      IF(PRESENT(POUT_RV) .OR. PRESENT(POUT_RC) .OR. &
          &PRESENT(POUT_RI) .OR. PRESENT(POUT_TH)) THEN
        DO JI=1,IIU
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
          IF(PRESENT(POUT_RV)) POUT_RV(JI,JJ,JK)=ZRV(JI,JJ,JK)
          IF(PRESENT(POUT_RC)) POUT_RC(JI,JJ,JK)=ZRC(JI,JJ,JK)
          IF(PRESENT(POUT_RI)) POUT_RI(JI,JJ,JK)=ZRI(JI,JJ,JK)
          IF(PRESENT(POUT_TH)) POUT_TH(JI,JJ,JK)=ZT(JI,JJ,JK) / PEXN(JI,JJ,JK)
        ENDDO
      ENDIF
    ENDIF
  ENDDO
ENDDO
!
!
!*       6.  STORE THE BUDGET TERMS
!            ----------------------
!
IF (LBUDGET_RV) CALL BUDGET (PRVS(:,:,:) * PRHODJ(:,:,:),6,HBUNAME//'_BU_RRV',YDDDH, YDLDDH, YDMDDH)
IF (LBUDGET_RC) CALL BUDGET (PRCS(:,:,:) * PRHODJ(:,:,:),7,HBUNAME//'_BU_RRC',YDDDH, YDLDDH, YDMDDH)
IF (LBUDGET_RI) CALL BUDGET (PRIS(:,:,:) * PRHODJ(:,:,:),9,HBUNAME//'_BU_RRI',YDDDH, YDLDDH, YDMDDH)
IF (LBUDGET_TH) CALL BUDGET (PTHS(:,:,:) * PRHODJ(:,:,:),4,HBUNAME//'_BU_RTH',YDDDH, YDLDDH, YDMDDH)
!
!------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('ICE_ADJUST',1,ZHOOK_HANDLE)

CONTAINS
SUBROUTINE ITERATION(PRV_IN,PRC_IN,PRI_IN,PRV_OUT,PRC_OUT,PRI_OUT)

REAL, DIMENSION(:,:,:), CONTIGUOUS, INTENT(IN) :: PRV_IN ! Water vapor m.r. to adjust in input
REAL, DIMENSION(:,:,:), CONTIGUOUS, INTENT(IN) :: PRC_IN ! Cloud water m.r. to adjust in input
REAL, DIMENSION(:,:,:), CONTIGUOUS, INTENT(IN) :: PRI_IN ! Cloud ice   m.r. to adjust in input
REAL, DIMENSION(:,:,:), CONTIGUOUS, INTENT(OUT) :: PRV_OUT ! Water vapor m.r. to adjust in output
REAL, DIMENSION(:,:,:), CONTIGUOUS, INTENT(OUT) :: PRC_OUT ! Cloud water m.r. to adjust in output
REAL, DIMENSION(:,:,:), CONTIGUOUS, INTENT(OUT) :: PRI_OUT ! Cloud ice   m.r. to adjust in output

!*       2.4    compute the specific heat for moist air (Cph) at t+1

SELECT CASE(KRR)
  CASE(7)
    ZCPH(:,:,:) = XCPD + XCPV * PRV_IN(:,:,:)                             &
                       + XCL  * (PRC_IN(:,:,:) + PRR(:,:,:))             &
                       + XCI  * (PRI_IN(:,:,:) + PRS(:,:,:) + PRG(:,:,:) + PRH(:,:,:))
  CASE(6)
    ZCPH(:,:,:) = XCPD + XCPV * PRV_IN(:,:,:)                             &
                       + XCL  * (PRC_IN(:,:,:) + PRR(:,:,:))             &
                       + XCI  * (PRI_IN(:,:,:) + PRS(:,:,:) + PRG(:,:,:))
  CASE(5)
    ZCPH(:,:,:) = XCPD + XCPV * PRV_IN(:,:,:)                             &
                       + XCL  * (PRC_IN(:,:,:) + PRR(:,:,:))             &
                       + XCI  * (PRI_IN(:,:,:) + PRS(:,:,:))
  CASE(3)
    ZCPH(:,:,:) = XCPD + XCPV * PRV_IN(:,:,:)               &
                       + XCL  * (PRC_IN(:,:,:) + PRR(:,:,:))
  CASE(2)
    ZCPH(:,:,:) = XCPD + XCPV * PRV_IN(:,:,:) &
                       + XCL  * PRC_IN(:,:,:)
END SELECT
!
IF ( OSUBG_COND ) THEN
  !
  !*       3.     SUBGRID CONDENSATION SCHEME
  !               ---------------------------
  !
  !   PSRC= s'rci'/Sigma_s^2
  !   ZT is INOUT
 
     CALL CONDENSATION(IIU, IJU, IKU, IIB, IIE, IJB, IJE, IKB, IKE, KKL, &
       HFRAC_ICE,                                                        &
       PPABST, PZZ, ZT, PRV_IN, PRV_OUT, PRC_IN, PRC_OUT, PRI_IN, PRI_OUT, PRS, PRG, PSIGS, PMFCONV, PCLDFR, &
       PSRCS, .TRUE., OSIGMAS,                                           &
       OCND2, PSIGQSAT, YSPP_PSIGQSAT,YSPP_ICE_CLD_WGT,                  &
       PLV=ZLV, PLS=ZLS, PCPH=ZCPH)
ELSE
  !
  !*       4.     ALL OR NOTHING CONDENSATION SCHEME
  !                            FOR MIXED-PHASE CLOUD
  !               -----------------------------------------------
     ZSIGS(:,:,:)=0.
     CALL CONDENSATION(IIU, IJU, IKU, IIB, IIE, IJB, IJE, IKB, IKE, KKL, &
       HFRAC_ICE,                                                        &
       PPABST, PZZ, ZT, PRV_IN, PRV_OUT, PRC_IN, PRC_OUT, PRI_IN, PRI_OUT, PRS, PRG, ZSIGS, PMFCONV, PCLDFR, &
       PSRCS, .TRUE., OSIGMAS=.TRUE.,                                    &
       OCND2=OCND2, PSIGQSAT=0.,                                         & 
       YSPP_PSIGQSAT=YSPP_PSIGQSAT,                                      &
       YSPP_ICE_CLD_WGT=YSPP_ICE_CLD_WGT,                                &
       PLV=ZLV, PLS=ZLS, PCPH=ZCPH)
ENDIF

END SUBROUTINE ITERATION

END SUBROUTINE ICE_ADJUST 
