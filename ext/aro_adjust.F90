!     ######spl
      SUBROUTINE  ARO_ADJUST(KLON,KIDIA,KFDIA,KLEV,  KRR,  &
                                  HFRAC_ICE, HCONDENS, HLAMBDA3, OSUBG_COND, &
                                  OSIGMAS, CMICRO, OCND2, LHGT_QS, HSUBG_MF_PDF, &
                                  PTSTEP, PSIGQSAT, &
                                  PZZF, PRHODJ, PEXNREF, PRHODREF,&
                                  PPABSM, PTHT, PRT, PSIGS, &
                                  PMFCONV, PRC_MF, PRI_MF, PCF_MF, &
                                  PTHS, PRS,  PSRCS, PCLDFR,&
                                  PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR, &
                                  PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,&
                                  YDDDH,YDLDDH,YDMDDH,&
                                  YSPP_PSIGQSAT,YSPP_ICE_CLD_WGT)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##########################################################################
!
!!****  * -  compute the  resolved clouds and precipitation
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the  microphysical sources
!!    related to the resolved clouds and precipitation
!!
!!
!!**  METHOD
!!    ------
!!      The main actions of this routine is to call the routines computing the
!!    microphysical sources. Before that:
!!        - it computes the real absolute pressure,
!!        - negative values of the current guess of all mixing ratio are removed.
!!          This is done by a global filling algorithm based on a multiplicative
!!          method (Rood, 1987), in order to conserved the total mass in the
!!          simulation domain.
!!        - Sources are transformed in physical tendencies, by removing the
!!          multiplicative term Rhod*J.
!!        - External points values are filled owing to the use of cyclic
!!          l.b.c., in order to performe computations on the full domain.
!!      After calling to microphysical routines, the physical tendencies are
!!    switched back to prognostic variables.
!!
!!
!!    EXTERNAL
!!    --------
!!      Subroutine FMLOOK: to recover the logical unit number linked to a FMfile
!!      Subroutine SLOW_TERMS: Computes the explicit microphysical sources
!!      Subroutine FAST_TERMS: Performs the saturation adjustment for l
!!      Subroutine RAIN_ICE  : Computes the explicit microphysical sources for i
!!      Subroutine ICE_ADJUST: Performs the saturation adjustment for i+l
!!      MIN_ll,SUM3D_ll : distributed functions equivalent to MIN and SUM
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains declarations of parameter variables
!!         JPHEXT       : Horizontal external points number
!!         JPVEXT       : Vertical external points number
!!      Module MODD_CST
!!          XP00               ! Reference pressure
!!          XRD                ! Gaz  constant for dry air
!!          XCPD               ! Cpd (dry air)
!!
!!    REFERENCE
!!    ---------
!!
!!      Documentation AROME
!!
!!    AUTHOR
!!    ------
!!    S.Malardel and Y.Seity
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/03/03
!!      T. Kovacic  11-05-05, Call to budgets for NEGA1_
!!      S. Riette ice for EDKF
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!      2016-11 S. Riette: new ice_adjust interface, add OLD3/OLD4 schemes
!!      2018-02 K.I Ivarsson : More outputs from OCND2 option
!!      2020-12 U. Andrae : Introduce SPP for HARMONIE-AROME
!!     R. El Khatib 24-Aug-2021 Optimizations
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY: CST
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM
USE MODD_NEB, ONLY: NEB
USE MODD_TURB_n, ONLY: TURBN
USE MODD_BUDGET, ONLY: TBUDGETDATA, NBUDGET_RI, TBUCONF
USE SPP_MOD_TYPE, ONLY : TSPP_CONFIG_TYPE, CLEAR_SPP_TYPE, APPLY_SPP
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
USE MODI_ICE_ADJUST
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
!
!
USE DDH_MIX , ONLY : TYP_DDH
USE YOMLDDH , ONLY : TLDDH
USE YOMMDDH , ONLY : TMDDH
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!

!
INTEGER,                  INTENT(IN)   :: KLON ! array length (NPROMA)
INTEGER,                  INTENT(IN)   :: KIDIA    !start index (=1)
INTEGER,                  INTENT(IN)   :: KFDIA    !end index (=KLON only if block is full)
INTEGER,                  INTENT(IN)   :: KLEV     !Number of vertical levels
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
CHARACTER*1,              INTENT(IN)   :: HFRAC_ICE
CHARACTER*80,             INTENT(IN)   :: HCONDENS
CHARACTER*4,              INTENT(IN)   :: HLAMBDA3 ! formulation for lambda3 coeff
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid Cond.
LOGICAL,                  INTENT(IN)   :: OSIGMAS  ! Switch for Sigma_s:
                                        ! use values computed in CONDENSATION
                                        ! or that from turbulence scheme
CHARACTER (LEN=4),        INTENT(IN)   :: CMICRO  ! Microphysics scheme
LOGICAL,                  INTENT(IN)   :: OCND2
LOGICAL,                  INTENT(IN)   :: LHGT_QS
CHARACTER*80,             INTENT(IN)   :: HSUBG_MF_PDF
REAL,                     INTENT(IN)   :: PTSTEP   ! Time step
REAL,                     INTENT(IN)   :: PSIGQSAT ! coeff applied to qsat variance contribution
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PZZF     ! Height (z)
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PEXNREF ! Reference Exner function
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PRHODREF
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PPABSM  ! abs. pressure at time t-dt
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT) :: PRT     ! Moist variables at time t
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PSIGS   ! Sigma_s at time t
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PMFCONV ! convective mass flux
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PRC_MF, PRI_MF, PCF_MF
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(INOUT) :: PTHS  ! Theta source
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT) :: PRS   ! Moist  variable sources
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(OUT)   :: PSRCS ! Second-order flux
                                                 ! s'rc'/2Sigma_s2 at time t+1
                                                 ! multiplied by Lambda_3
REAL, DIMENSION(KLON,1,KLEV), INTENT(INOUT)   :: PCLDFR! Cloud fraction
REAL, DIMENSION(KLON,1,KLEV),   INTENT(OUT)   :: PICLDFR ! ice cloud fraction
REAL, DIMENSION(KLON,1,KLEV),   INTENT(OUT)   :: PWCLDFR ! water or mixed-phase cloud fraction
REAL, DIMENSION(KLON,1,KLEV),   INTENT(OUT)   :: PSSIO   ! Super-saturation with respect to ice in the 
                                                         ! supersaturated fraction
REAL, DIMENSION(KLON,1,KLEV),   INTENT(OUT)   :: PSSIU   ! Sub-saturation with respect to ice in the 
                                                         ! subsaturated fraction 
REAL, DIMENSION(KLON,1,KLEV),   INTENT(OUT)   :: PIFR    ! Ratio cloud ice moist part to dry part
!
REAL, DIMENSION(KLON,1,KLEV), INTENT(OUT)   :: PHLC_HRC
REAL, DIMENSION(KLON,1,KLEV), INTENT(OUT)   :: PHLC_HCF
REAL, DIMENSION(KLON,1,KLEV), INTENT(OUT)   :: PHLI_HRI
REAL, DIMENSION(KLON,1,KLEV), INTENT(OUT)   :: PHLI_HCF
!
TYPE(TYP_DDH), INTENT(INOUT), TARGET :: YDDDH
TYPE(TLDDH), INTENT(IN), TARGET :: YDLDDH
TYPE(TMDDH), INTENT(IN), TARGET :: YDMDDH
!
TYPE(TSPP_CONFIG_TYPE), INTENT(INOUT) :: YSPP_PSIGQSAT,YSPP_ICE_CLD_WGT
!
!*       0.2   Declarations of local variables :

CHARACTER*4               :: HBUNAME  ! Name of the budget
!
INTEGER :: JRR           ! Loop index for the moist and scalar variables
INTEGER :: JLON, JLEV
REAL :: ZT, ZTWOTSTEP
REAL, DIMENSION(KLON) :: ZLV,ZLS,ZCPH
LOGICAL :: LL(KLON)
REAL, DIMENSION(KLON,1,KLEV,0:KRR) :: ZRS
REAL, DIMENSION(KLON,1,KLEV) :: ZZZ
                                    ! model layer height
REAL  :: ZMASSTOT                   ! total mass  for one water category
                                    ! including the negative values
REAL  :: ZMASSPOS                   ! total mass  for one water category
                                    ! after removing the negative values
REAL  :: ZRATIO                     ! ZMASSTOT / ZMASSCOR
REAL  :: ZCOR(KLON)                 ! for the correction of negative rv
!
REAL, DIMENSION(KLON,1) :: ZSIGQSAT, ZICE_CLD_WGT
TYPE(TBUDGETDATA), DIMENSION(NBUDGET_RI) :: YLBUDGET !NBUDGET_RI is the one with the highest number
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_ADJUST',0,ZHOOK_HANDLE)

!Dimensions
CALL FILL_DIMPHYEX(YLDIMPHYEX, KLON, 1, KLEV, 0, KFDIA)

!
! Apply SPP perturbations
!

IF (YSPP_PSIGQSAT%LPERT) THEN
 CALL APPLY_SPP(YSPP_PSIGQSAT,KLON,1,KLON,PSIGQSAT,ZSIGQSAT)
ELSE
 ZSIGQSAT(:,:) = PSIGQSAT
ENDIF

IF (YSPP_ICE_CLD_WGT%LPERT) THEN
 CALL APPLY_SPP(YSPP_ICE_CLD_WGT,KLON,1,KLON,RAIN_ICE_PARAM%XFRMIN(21),ZICE_CLD_WGT)
ELSE
 ZICE_CLD_WGT(:,:) = RAIN_ICE_PARAM%XFRMIN(21)
ENDIF

HBUNAME='DEPI'
!
!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
!
! personal comment:  tranfering these variables to the
!                    microphysical routines would save
!                    computing time

! Well, getting rid of array syntax already saves a lot ;-) .REK.
!
!
!
!*       3.     REMOVE NEGATIVE VALUES
!               ----------------------
!
!*       3.1    Non local correction for precipitating species (Rood 87)
!
IF (CMICRO == 'KESS' .OR. CMICRO == 'ICE3' .OR. CMICRO == 'ICE2' &
    .OR.  CMICRO == 'C2R2' .OR. CMICRO == 'C3R5'.OR. CMICRO == 'ICE4' &
    .OR. CMICRO == 'OLD3' .OR. CMICRO == 'OLD4') THEN

  DO JRR = 3,KRR
    SELECT CASE (JRR)
      CASE(3,5,6,7) ! rain, snow, graupel and hail

        IF ( MINVAL( PRS(KIDIA:KFDIA,:,:,JRR)) < 0.0 ) THEN
! For AROME, we cannot use MAX_ll so that according to JPP's advises
!  we only correct negative values but not the total mass
! compute the total water mass computation
!
!          ZMASSTOT = MAX( 0. , SUM( PRS(:,:,:,JRR) ))
!
! remove the negative values
!
          PRS(KIDIA:KFDIA,:,:,JRR) = MAX( 0., PRS(KIDIA:KFDIA,:,:,JRR) )
!
! compute the new total mass
!
!          ZMASSPOS = MAX(1.E-60,SUM( PRS(:,:,:,JRR) ))
!
! correct again in such a way to conserve the total mass
!
!          ZRATIO = ZMASSTOT / ZMASSPOS
!          PRS(:,:,:,JRR) = PRS(:,:,:,JRR) * ZRATIO

        END IF
    END SELECT
  END DO
END IF
!
!*       3.2    Adjustement for liquid and solid cloud
!

ZTWOTSTEP=2.*PTSTEP

SELECT CASE ( CMICRO )
!
!
  CASE('ICE2','ICE3','ICE4', 'OLD3', 'OLD4')

  DO JLEV=1,KLEV

    DO JLON=KIDIA,KFDIA
      ZT = PTHT(JLON,1,JLEV)*PEXNREF(JLON,1,JLEV)
      ZLV(JLON)=CST%XLVTT +(CST%XCPV-CST%XCL) *(ZT-CST%XTT)
      ZLS(JLON)=CST%XLSTT +(CST%XCPV-CST%XCI) *(ZT-CST%XTT)
      ZCPH(JLON)=CST%XCPD +CST%XCPV*2.*PTSTEP*PRS(JLON,1,JLEV,1)
    ENDDO

    DO JLON=KIDIA,KFDIA
      IF (PRS(JLON,1,JLEV,4) < 0.) THEN
        PRS(JLON,1,JLEV,1) = PRS(JLON,1,JLEV,1) + PRS(JLON,1,JLEV,4)
        PTHS(JLON,1,JLEV)  = PTHS(JLON,1,JLEV)  - PRS(JLON,1,JLEV,4) * ZLS(JLON) / ZCPH(JLON) / PEXNREF(JLON,1,JLEV)
        PRS(JLON,1,JLEV,4) = 0.
      ENDIF
    ENDDO
!
!   cloud
    DO JLON=KIDIA,KFDIA
      IF (PRS(JLON,1,JLEV,2) < 0.) THEN
        PRS(JLON,1,JLEV,1) = PRS(JLON,1,JLEV,1) + PRS(JLON,1,JLEV,2)
        PTHS(JLON,1,JLEV)  = PTHS(JLON,1,JLEV)  - PRS(JLON,1,JLEV,2) * ZLV(JLON) / ZCPH(JLON) / PEXNREF(JLON,1,JLEV)
        PRS(JLON,1,JLEV,2) = 0.
      ENDIF
    ENDDO
!
! if rc or ri are positive, we can correct negative rv
!   cloud
    DO JLON=KIDIA,KFDIA
      LL(JLON) = (PRS(JLON,1,JLEV,1) <0.) .AND. (PRS(JLON,1,JLEV,2)> 0.) 
      IF (LL(JLON)) THEN
#ifdef REPRO48
        ZCOR(JLON)=PRS(JLON,1,JLEV,2)
#else
        ZCOR(JLON)=MIN(-PRS(JLON,1,JLEV,1),PRS(JLON,1,JLEV,2))
#endif
      ENDIF
    ENDDO
    DO JLON=KIDIA,KFDIA
      IF (LL(JLON)) THEN
        PRS(JLON,1,JLEV,1) = PRS(JLON,1,JLEV,1) + ZCOR(JLON)
        PTHS(JLON,1,JLEV)  = PTHS(JLON,1,JLEV)  - ZCOR(JLON) * ZLV(JLON) / ZCPH(JLON) / PEXNREF(JLON,1,JLEV)
#ifdef REPRO48
        PRS(JLON,1,JLEV,2) = 0.
#else
        PRS(JLON,1,JLEV,2) = PRS(JLON,1,JLEV,2) - ZCOR(JLON)
#endif
      ENDIF
    ENDDO

!   ice
    IF (KRR > 3) THEN
      DO JLON=KIDIA,KFDIA
        LL(JLON) = (PRS(JLON,1,JLEV,1) < 0.).AND.(PRS(JLON,1,JLEV,4) > 0.)
        IF (LL(JLON)) THEN
          ZCOR(JLON)=MIN(-PRS(JLON,1,JLEV,1),PRS(JLON,1,JLEV,4))
        ENDIF
      ENDDO
      DO JLON=KIDIA,KFDIA
        IF (LL(JLON)) THEN
          PRS(JLON,1,JLEV,1) = PRS(JLON,1,JLEV,1) + ZCOR(JLON)
          PTHS(JLON,1,JLEV)  = PTHS(JLON,1,JLEV)  - ZCOR(JLON) * ZLS(JLON) / ZCPH(JLON) / PEXNREF(JLON,1,JLEV)
          PRS(JLON,1,JLEV,4) = PRS(JLON,1,JLEV,4) - ZCOR(JLON)
        ENDIF
      ENDDO
    ENDIF

  ENDDO ! JLEV
!
END SELECT
!
!
!*       3.3  STORE THE BUDGET TERMS
!            ----------------------
!
!IF (LBUDGET_RV) CALL BUDGET (PRS(:,:,:,1) * PRHODJ(:,:,:), 6,'NEGA_BU_RRV',YDDDH)
!IF (LBUDGET_RC) CALL BUDGET (PRS(:,:,:,2) * PRHODJ(:,:,:), 7,'NEGA_BU_RRC',YDDDH)
!IF (LBUDGET_RR) CALL BUDGET (PRS(:,:,:,3) * PRHODJ(:,:,:), 8,'NEGA_BU_RRR',YDDDH)
!IF (LBUDGET_RI) CALL BUDGET (PRS(:,:,:,4) * PRHODJ(:,:,:) ,9,'NEGA_BU_RRI',YDDDH)
!IF (LBUDGET_RS) CALL BUDGET (PRS(:,:,:,5) * PRHODJ(:,:,:),10,'NEGA_BU_RRS',YDDDH)
!IF (LBUDGET_RG) CALL BUDGET (PRS(:,:,:,6) * PRHODJ(:,:,:),11,'NEGA_BU_RRG',YDDDH)
!IF (LBUDGET_RH) CALL BUDGET (PRS(:,:,:,7) * PRHODJ(:,:,:),12,'NEGA_BU_RRH',YDDDH)
!IF (LBUDGET_TH) CALL BUDGET (PTHS(:,:,:)  * PRHODJ(:,:,:),4,'NEGA_BU_RTH',YDDDH)

DO JRR=1, NBUDGET_RI
  YLBUDGET(JRR)%NBUDGET=JRR
  YLBUDGET(JRR)%YDDDH=>YDDDH
  YLBUDGET(JRR)%YDLDDH=>YDLDDH
  YLBUDGET(JRR)%YDMDDH=>YDMDDH
ENDDO

!
!-------------------------------------------------------------------------------
!

!*       9.     MIXED-PHASE MICROPHYSICAL SCHEME (WITH 3 ICE SPECIES)
!               -----------------------------------------------------
!
DO JRR = 0,KRR
  IF (JRR==0) THEN
    ZRS(KIDIA:KFDIA,:,:,0)=PTHS(KIDIA:KFDIA,:,:)*2.*PTSTEP
  ELSE
    ZRS(KIDIA:KFDIA,:,:,JRR)=PRS(KIDIA:KFDIA,:,:,JRR)*2.*PTSTEP
  ENDIF
ENDDO
ZZZ(KIDIA:KFDIA,:,:) =  PZZF(KIDIA:KFDIA,:,:)
!
!*       9.2    Perform the saturation adjustment over cloud ice and cloud water
!
IF (KRR==6) THEN
  CALL ICE_ADJUST ( YLDIMPHYEX, CST=CST, ICEP=RAIN_ICE_PARAM, NEB=NEB, TURBN=TURBN, BUCONF=TBUCONF, KRR=KRR,&
    & HFRAC_ICE=HFRAC_ICE, HCONDENS=HCONDENS, HLAMBDA3=HLAMBDA3, HBUNAME=HBUNAME, &
    & OSUBG_COND=OSUBG_COND, OSIGMAS=OSIGMAS, &
    & OCND2=OCND2, LHGT_QS=LHGT_QS, HSUBG_MF_PDF=HSUBG_MF_PDF, &
    & PTSTEP=ZTWOTSTEP,PSIGQSAT=ZSIGQSAT, &
    & PRHODJ=PRHODJ ,PEXNREF=PEXNREF, PRHODREF=PRHODREF,   &
    & PSIGS=PSIGS, LMFCONV=.TRUE., PMFCONV=PMFCONV, PPABST=PPABSM, PZZ=ZZZ,    &
    & PEXN=PEXNREF, PCF_MF=PCF_MF,PRC_MF=PRC_MF,PRI_MF=PRI_MF, &
    & PICLDFR=PICLDFR, PWCLDFR=PWCLDFR, & 
    & PSSIO=PSSIO, PSSIU=PSSIU, PIFR=PIFR, &
    & PRV=ZRS(:,:,:,1), PRC=ZRS(:,:,:,2),  &
    & PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2), &
    & PTH=ZRS(:,:,:,0), PTHS=PTHS,OCOMPUTE_SRC=.TRUE.,PSRCS=PSRCS, PCLDFR=PCLDFR, &
    & PRR=ZRS(:,:,:,3), &
    & PRI=ZRS(:,:,:,4), PRIS=PRS(:,:,:,4), &
    & PRS=ZRS(:,:,:,5), &
    & PRG=ZRS(:,:,:,6), &
    & TBUDGETS=YLBUDGET, KBUDGETS=SIZE(YLBUDGET), &
    & PICE_CLD_WGT=ZICE_CLD_WGT(:,:), &
    & PHLC_HRC=PHLC_HRC(:,:,:), PHLC_HCF=PHLC_HCF(:,:,:), &
    & PHLI_HRI=PHLI_HRI(:,:,:), PHLI_HCF=PHLI_HCF(:,:,:))
ELSE
  CALL ICE_ADJUST ( YLDIMPHYEX, CST=CST, ICEP=RAIN_ICE_PARAM, NEB=NEB, TURBN=TURBN, BUCONF=TBUCONF, KRR=KRR,&
    & HFRAC_ICE=HFRAC_ICE, HCONDENS=HCONDENS, HLAMBDA3=HLAMBDA3, HBUNAME=HBUNAME,    &
    & OSUBG_COND=OSUBG_COND, OSIGMAS=OSIGMAS, &
    & OCND2=OCND2, LHGT_QS=LHGT_QS, HSUBG_MF_PDF=HSUBG_MF_PDF, &
    & PTSTEP=ZTWOTSTEP,PSIGQSAT=ZSIGQSAT, &
    & PRHODJ=PRHODJ ,PEXNREF=PEXNREF, PRHODREF=PRHODREF,   &
    & PSIGS=PSIGS, LMFCONV=.TRUE., PMFCONV=PMFCONV, PPABST=PPABSM, PZZ=ZZZ,    &
    & PEXN=PEXNREF, PCF_MF=PCF_MF,PRC_MF=PRC_MF,PRI_MF=PRI_MF, &
    & PICLDFR=PICLDFR, PWCLDFR=PWCLDFR, & 
    & PSSIO=PSSIO, PSSIU=PSSIU, PIFR=PIFR, &
    & PRV=ZRS(:,:,:,1), PRC=ZRS(:,:,:,2), &
    & PRVS=PRS(:,:,:,1), PRCS=PRS(:,:,:,2), &
    & PTH=ZRS(:,:,:,0), PTHS=PTHS,OCOMPUTE_SRC=.TRUE.,PSRCS=PSRCS, PCLDFR=PCLDFR, &
    & PRR=ZRS(:,:,:,3), &
    & PRI=ZRS(:,:,:,4), PRIS=PRS(:,:,:,4), &
    & PRS=ZRS(:,:,:,5), &
    & PRG=ZRS(:,:,:,6), &
    & TBUDGETS=YLBUDGET, KBUDGETS=SIZE(YLBUDGET), &
    & PICE_CLD_WGT=ZICE_CLD_WGT(:,:), &
    & PRH=ZRS(:,:,:,7), &
    & PHLC_HRC=PHLC_HRC(:,:,:), PHLC_HCF=PHLC_HCF(:,:,:), & 
    & PHLI_HRI=PHLI_HRI(:,:,:), PHLI_HCF=PHLI_HCF(:,:,:))
ENDIF

CALL CLEAR_SPP_TYPE(YSPP_PSIGQSAT)
CALL CLEAR_SPP_TYPE(YSPP_ICE_CLD_WGT)

!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_ADJUST',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_ADJUST
