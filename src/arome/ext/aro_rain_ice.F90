!     ######spl
      SUBROUTINE  ARO_RAIN_ICE(PHYEX, &
                             & KKA,KKU,KKL,KLON,KLEV, &
                             & KIDIA, KFDIA, KRR, &
                             & CMICRO, &
                             & PTSTEP, PDZZ, PRHODJ, PRHODREF, PEXNREF,&
                             & PPABSM, PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF, PTHT, PRT, PSIGS,PCLDFR, &
                             & PTHS, PRS, PEVAP,  &
                             & PCIT, PSEA, PTOWN,   &
                             & PICLDFR, PSSIO, PSSIU, PIFR,  &
                             & PINPRR,PINPRS,PINPRG,PINPRH,PFPR,     &
                             & OAERONRT, OAEIFN, PCLDROP, PIFNNC, &
                             & YDDDH, YDLDDH, YDMDDH, &
                             & YSPP_ICENU,YSPP_KGN_ACON,YSPP_KGN_SBGR)

      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
!!      15/05/05 T. Kovacic, budgets for negative correction
!!      29/09/08 Y. Seity, add PEVAP for chemistry
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!      2013-11, D. Degrauwe: Export upper-air precipitation fluxes PFPR.
!!      2013-11 S. Riette, subgrid precipitation
!!      2014-11 S. Riette, ICE3/ICE4 modified, old versions under OLD3/OLD4
!!      2014-11 S. Riette, ICE3/ICE4 modified, old versions under OLD3/OLD4
!!      2020-12 U. Andrae : Introduce SPP for HARMONIE-AROME
!!      2018-02 K.I: Ivarsson: More inputs to OCND2-option for saving computing time.
!!     R. El Khatib 24-Aug-2021 Specific cache-blocking factor for microphysics
!!      2021-01 D. Martin-Perez: n.r.t. aerosoles
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
USE MODD_BUDGET, ONLY: TBUDGETDATA_PTR, NBUDGET_RH
USE MODE_BUDGET_IAL, ONLY: TBUDGETDATA_IAL, BUDGET_DDH
USE MODD_LES, ONLY: TLES, LES_ASSOCIATE !only used by rain_ice_old
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
!
USE MODI_RAIN_ICE
!
USE SPP_MOD_TYPE, ONLY : TSPP_CONFIG_TYPE, APPLY_SPP
!
USE MODI_RAIN_ICE_OLD
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
TYPE(PHYEX_t),            INTENT(IN) :: PHYEX
INTEGER,                  INTENT(IN) :: KKA  !near ground array index
INTEGER,                  INTENT(IN) :: KKU  !uppest atmosphere array index
INTEGER,                  INTENT(IN) :: KKL  !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN) :: KLON     !NPROMA under CPG
INTEGER,                  INTENT(IN) :: KLEV     !Number of vertical levels
INTEGER,                  INTENT(IN) :: KIDIA
INTEGER,                  INTENT(IN) :: KFDIA
INTEGER,                  INTENT(IN) :: KRR      ! Number of moist variables
CHARACTER (LEN=4),        INTENT(IN) :: CMICRO  ! Microphysics scheme
REAL,                     INTENT(IN) :: PTSTEP   ! Time step
!
!
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PDZZ     ! Height (z)
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PRHODREF! Reference dry air density
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PEXNREF ! Reference Exner function
!
!
REAL, DIMENSION(KLON,KLEV),     INTENT(IN)    :: PPABSM  ! abs. pressure at time t-dt
REAL, DIMENSION(KLON,KLEV),     INTENT(INOUT) :: PHLC_HRC
REAL, DIMENSION(KLON,KLEV),     INTENT(INOUT) :: PHLC_HCF
REAL, DIMENSION(KLON,KLEV),     INTENT(INOUT) :: PHLI_HRI
REAL, DIMENSION(KLON,KLEV),     INTENT(INOUT) :: PHLI_HCF
REAL, DIMENSION(KLON,KLEV),     INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(KLON,KLEV,KRR), INTENT(INOUT) :: PRT   ! Moist variables at time t
REAL, DIMENSION(KLON,KLEV),     INTENT(IN)    :: PSIGS   ! Sigma_s at time t
REAL, DIMENSION(KLON,KLEV),     INTENT(IN)    :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(KLON,KLEV),     INTENT(INOUT) :: PTHS  ! Theta source
REAL, DIMENSION(KLON,KLEV,KRR), INTENT(INOUT) :: PRS   ! Moist  variable sources
REAL, DIMENSION(KLON,KLEV),     INTENT(INOUT) :: PEVAP ! Rain evap profile
!
REAL, DIMENSION(KLON,KLEV),     INTENT(INOUT) :: PCIT  ! Pristine ice number
                                                         ! concentration at time t
REAL, DIMENSION(KLON),          INTENT(IN)    :: PSEA  ! Land sea mask
REAL, DIMENSION(KLON),          INTENT(IN)    :: PTOWN  ! Town mask
! input from aro_adjust / condensation with OCND2, dummy if OCND2 = F
REAL, DIMENSION(KLON,KLEV),     INTENT(IN)    :: PICLDFR ! ice cloud fraction
REAL, DIMENSION(KLON,KLEV),     INTENT(IN)    :: PSSIO   ! Super-saturation with respect to ice in the  
                                                           ! supersaturated fraction
REAL, DIMENSION(KLON,KLEV),     INTENT(IN)    :: PSSIU   ! Sub-saturation with respect to ice in the  
                                                           ! subsaturated fraction 
REAL, DIMENSION(KLON,KLEV),     INTENT(INOUT) :: PIFR    ! Ratio cloud ice moist part to dry part 
!
! input from aro_adjust / condensation with OCND2 END.
!
REAL, DIMENSION(KLON),          INTENT(OUT)   :: PINPRR! Rain instant precip
REAL, DIMENSION(KLON),          INTENT(OUT)   :: PINPRS! Snow instant precip
REAL, DIMENSION(KLON),          INTENT(OUT)   :: PINPRG! Graupel instant precip
REAL, DIMENSION(KLON),          INTENT(OUT)   :: PINPRH! Hail instant precip
REAL, DIMENSION(KLON,KLEV,KRR), INTENT(INOUT) :: PFPR ! upper-air precip
!
LOGICAL,                        INTENT(IN) :: OAERONRT ! Switch for Near-Real-Time aerosol
LOGICAL,                        INTENT(IN) :: OAEIFN   ! Switch to activate ice nuclei
REAL, DIMENSION(KLON,KLEV),     INTENT(IN) :: PCLDROP  ! Cloud droplet number concentration
REAL, DIMENSION(KLON,KLEV),     INTENT(IN) :: PIFNNC   ! Ice Forming Nuclei number concentration
!
TYPE(TYP_DDH), INTENT(INOUT), TARGET :: YDDDH
TYPE(TLDDH),   INTENT(IN),    TARGET :: YDLDDH
TYPE(TMDDH),   INTENT(IN),    TARGET :: YDMDDH
!
TYPE(TSPP_CONFIG_TYPE), INTENT(INOUT) :: YSPP_ICENU 
TYPE(TSPP_CONFIG_TYPE), INTENT(INOUT) :: YSPP_KGN_ACON 
TYPE(TSPP_CONFIG_TYPE), INTENT(INOUT) :: YSPP_KGN_SBGR 
!
!
!*       0.2   Declarations of local variables :
INTEGER :: JRR           ! Loop index for the moist and scalar variables
INTEGER :: JLON, JLEV
!
!
!
REAL, DIMENSION(KLON,KLEV) :: ZT
REAL, DIMENSION(KLON,KLEV) :: ZLV
REAL, DIMENSION(KLON,KLEV) :: ZLS
REAL, DIMENSION(KLON,KLEV) :: ZCPH
REAL, DIMENSION(KLON,KLEV) :: ZCOR
REAL, DIMENSION(KLON)      :: ZINDEP     ! surf cloud deposition (already contained in sedimentation)
REAL, DIMENSION(KLON,KLEV) :: ZRAINFR
REAL, DIMENSION(KLON)      :: ZICENU
REAL, DIMENSION(KLON)      :: ZKGN_ACON
REAL, DIMENSION(KLON)      :: ZKGN_SBGR
REAL, DIMENSION(KLON)      :: ZINPRC    ! surf cloud sedimentation
                                    ! for the correction of negative rv
REAL  :: ZMASSTOT                   ! total mass  for one water category
                                    ! including the negative values
REAL  :: ZMASSPOS                   ! total mass  for one water category
                                    ! after removing the negative values
REAL  :: ZRATIO                     ! ZMASSTOT / ZMASSCOR
REAL  :: ZTWOTSTEP

REAL, DIMENSION(KLON,KLEV) :: ZWORK
TYPE(TBUDGETDATA_IAL), DIMENSION(NBUDGET_RH), TARGET :: YLBUDGET !NBUDGET_RH is the one with the highest number
TYPE(TBUDGETDATA_PTR), DIMENSION(NBUDGET_RH) :: YLBUDGET_PTR
REAL  :: ZTHVREFZIKB ! for electricity use only
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
LOGICAL, DIMENSION(KLON,KLEV) :: LLMICRO
INTEGER :: ISIZE
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_RAIN_ICE',0,ZHOOK_HANDLE)

!Dimensions
CALL FILL_DIMPHYEX(YLDIMPHYEX, KLON, 1, KLEV, 0, KIDIA, KFDIA)

LLMICRO(:,:) = .FALSE.

!LES init (only for rain_ice_old)
CALL LES_ASSOCIATE()
TLES%LLES=.FALSE.
TLES%LLES_CALL=.FALSE.

ZTWOTSTEP=2*PTSTEP
ZINPRC=0.
PINPRH=0.
IF (PHYEX%MISC%CELEC /='NONE') THEN
CALL ABOR1('ARO_RAIN_ICE : CELEC ELECTRICITY SCHEME NOT YET CORRECLY PLUGGED HERE')
! The following value of ZTHVREFZIKB must be removed from the electricity scheme or computed correctly here
ELSE
  ZTHVREFZIKB = 0. ! for electricity use only
END IF

!Mask to limit computation
IF ( KRR == 7 ) THEN
  IF (CMICRO /= 'ICE4') THEN
    CALL ABOR1('ARO_RAIN_ICE : KRR==7 NOT COMPATIBLE WITH CMICRO /= ICE4')
  ENDIF
END IF


!------------------------------------------------------------------------------
!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
!
!  complete the vertical boundaries
!
!
! personal comment:  tranfering these variables to the
!                    microphysical routines would save
!                    computing time
!
DO JLEV = 1,KLEV
  DO JLON = KIDIA,KFDIA
    ZT(JLON,JLEV) = PTHT(JLON,JLEV)*PEXNREF(JLON,JLEV)
  ENDDO
ENDDO

DO JLEV = 1,KLEV
  DO JLON = KIDIA,KFDIA
    ZLV(JLON,JLEV) = PHYEX%CST%XLVTT + (PHYEX%CST%XCPV - PHYEX%CST%XCL)*(ZT(JLON,JLEV) - PHYEX%CST%XTT)
  ENDDO
ENDDO

DO JLEV = 1,KLEV
  DO JLON = KIDIA,KFDIA
    ZLS(JLON,JLEV) = PHYEX%CST%XLSTT + (PHYEX%CST%XCPV - PHYEX%CST%XCI)*(ZT(JLON,JLEV) - PHYEX%CST%XTT)
  ENDDO
ENDDO

DO JLEV = 1,KLEV
  DO JLON = KIDIA,KFDIA
    ZCPH(JLON,JLEV) = PHYEX%CST%XCPD + PHYEX%CST%XCPV*2.*PTSTEP*PRS(JLON,JLEV,1)
  ENDDO
ENDDO
!
!
!*       3.     REMOVE NEGATIVE VALUES
!               ----------------------
!
!*       3.1    Non local correction for precipitating species (Rood 87)
!
IF (CMICRO == 'KESS' .OR. CMICRO == 'ICE3' .OR. CMICRO == 'ICE2' &
    .OR.  CMICRO == 'C2R2' .OR. CMICRO == 'C3R5'.OR. CMICRO == 'ICE4') THEN

  ! For AROME, we cannot use MAX_ll so that according to JPP's advises
  ! we only correct negative values but not the total mass
  ! compute the total water mass computation

  IF (KRR >= 3) THEN
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        PRS(JLON,JLEV,3) = MAX(0., PRS(JLON,JLEV,3))
      ENDDO
    ENDDO
  ENDIF

  IF (KRR >= 5) THEN
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        PRS(JLON,JLEV,5) = MAX(0., PRS(JLON,JLEV,5))
      ENDDO
    ENDDO
  ENDIF
  IF (KRR >= 6) THEN
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        PRS(JLON,JLEV,6) = MAX(0., PRS(JLON,JLEV,6))
      ENDDO
    ENDDO
  ENDIF
  IF (KRR >= 7) THEN
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        PRS(JLON,JLEV,7) = MAX(0., PRS(JLON,JLEV,7))
      ENDDO
    ENDDO
  ENDIF
END IF

!
!*       3.2    Adjustement for liquid and solid cloud
!

IF (CMICRO == 'ICE2' .OR. CMICRO == 'ICE3' .OR. CMICRO == 'ICE4') THEN

  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      IF (PRS(JLON,JLEV,4) < 0.) THEN
        PRS(JLON,JLEV,1) = PRS(JLON,JLEV,1) + PRS(JLON,JLEV,4)
        PTHS(JLON,JLEV) = PTHS(JLON,JLEV) &
                        & - PRS(JLON,JLEV,4) * ZLS(JLON,JLEV) / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
        PRS(JLON,JLEV,4) = 0.
      ENDIF
    ENDDO
  ENDDO

! cloud
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      IF (PRS(JLON,JLEV,2) < 0.) THEN
        PRS(JLON,JLEV,1) = PRS(JLON,JLEV,1) + PRS(JLON,JLEV,2)
        PTHS(JLON,JLEV) = PTHS(JLON,JLEV) &
                        & - PRS(JLON,JLEV,2) * ZLV(JLON,JLEV) / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
        PRS(JLON,JLEV,2) = 0.
      ENDIF
    ENDDO
  ENDDO

! if rc or ri are positive, we can correct negative rv
! cloud
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      IF ((PRS(JLON,JLEV,1) < 0.) .AND. (PRS(JLON,JLEV,2) > 0.)) THEN
        PRS(JLON,JLEV,1) = PRS(JLON,JLEV,1) + PRS(JLON,JLEV,2)
        PTHS(JLON,JLEV) = PTHS(JLON,JLEV) &
                        & - PRS(JLON,JLEV,2) * ZLV(JLON,JLEV) / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
        PRS(JLON,JLEV,2) = 0.
      ENDIF
    ENDDO
  ENDDO

! ice
  IF(KRR > 3) THEN
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        IF ((PRS(JLON,JLEV,1) < 0.) .AND. (PRS(JLON,JLEV,4) > 0.)) THEN
          ZCOR(JLON,JLEV)=MIN(-PRS(JLON,JLEV,1),PRS(JLON,JLEV,4))
          PRS(JLON,JLEV,1) = PRS(JLON,JLEV,1) + ZCOR(JLON,JLEV)
          PTHS(JLON,JLEV) = PTHS(JLON,JLEV) &
                          & - ZCOR(JLON,JLEV) * ZLS(JLON,JLEV) / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
          PRS(JLON,JLEV,4) = PRS(JLON,JLEV,4) -ZCOR(JLON,JLEV)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
ENDIF
!
!
!*       3.3  STORE THE BUDGET TERMS
!            ----------------------
IF (PHYEX%MISC%TBUCONF%LBUDGET_RV) THEN
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      ZWORK(JLON,JLEV) = PRS(JLON,JLEV,1) * PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK, 6,'NEGA_BU_RRV',YDDDH, YDLDDH, YDMDDH)
ENDIF

IF (PHYEX%MISC%TBUCONF%LBUDGET_RC) THEN
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      ZWORK(JLON,JLEV) = PRS(JLON,JLEV,2) * PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK, 7,'NEGA_BU_RRC',YDDDH, YDLDDH, YDMDDH)
ENDIF

IF (PHYEX%MISC%TBUCONF%LBUDGET_RR) THEN
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      ZWORK(JLON,JLEV) = PRS(JLON,JLEV,3) * PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK, 8,'NEGA_BU_RRR',YDDDH, YDLDDH, YDMDDH)
ENDIF

IF (PHYEX%MISC%TBUCONF%LBUDGET_RI) THEN
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      ZWORK(JLON,JLEV) = PRS(JLON,JLEV,4) * PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK, 9,'NEGA_BU_RRI',YDDDH, YDLDDH, YDMDDH)
ENDIF

IF (PHYEX%MISC%TBUCONF%LBUDGET_RS) THEN
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      ZWORK(JLON,JLEV) = PRS(JLON,JLEV,5) * PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK,10,'NEGA_BU_RRS',YDDDH, YDLDDH, YDMDDH)
ENDIF

IF (PHYEX%MISC%TBUCONF%LBUDGET_RG) THEN
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      ZWORK(JLON,JLEV) = PRS(JLON,JLEV,6) * PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK,11,'NEGA_BU_RRG',YDDDH, YDLDDH, YDMDDH)
ENDIF

IF (PHYEX%MISC%TBUCONF%LBUDGET_RH .AND. KRR == 7) THEN
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      ZWORK(JLON,JLEV) = PRS(JLON,JLEV,7) * PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK,12,'NEGA_BU_RRH',YDDDH, YDLDDH, YDMDDH)
ENDIF

IF (PHYEX%MISC%TBUCONF%LBUDGET_TH) THEN
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      ZWORK(JLON,JLEV) = PTHS(JLON,JLEV) * PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK, 4,'NEGA_BU_RTH',YDDDH, YDLDDH, YDMDDH)
ENDIF

DO JRR=1, NBUDGET_RH
  YLBUDGET(JRR)%NBUDGET=JRR
  YLBUDGET(JRR)%YDDDH=>YDDDH
  YLBUDGET(JRR)%YDLDDH=>YDLDDH
  YLBUDGET(JRR)%YDMDDH=>YDMDDH
  YLBUDGET_PTR(JRR)%PTR=>YLBUDGET(JRR)
ENDDO
!
!
!-------------------------------------------------------------------------------
!

!*       9.     MIXED-PHASE MICROPHYSICAL SCHEME (WITH 3 ICE SPECIES)
!               -----------------------------------------------------
!
!*                   Compute the explicit microphysical sources
!
!
!
IF (CMICRO=='ICE4' .AND. PHYEX%PARAM_ICEN%LRED) THEN
    CALL RAIN_ICE(  YLDIMPHYEX, PHYEX%CST, PHYEX%PARAM_ICEN, PHYEX%RAIN_ICE_PARAMN, &
                 &  PHYEX%RAIN_ICE_DESCRN, PHYEX%ELEC_PARAM, PHYEX%ELEC_DESCR, &
                 &  PHYEX%MISC%TBUCONF, OELEC=PHYEX%MISC%OELEC, OSEDIM_BEARD=PHYEX%MISC%OSEDIM_BEARD, &
                 &  PTHVREFZIKB=ZTHVREFZIKB, PTSTEP=ZTWOTSTEP, &
                 &  KRR=KRR, PEXN=PEXNREF,            &
                 &  PDZZ=PDZZ, PRHODJ=PRHODJ, PRHODREF=PRHODREF, PEXNREF=PEXNREF,&
                 &  PPABST=PPABSM, PCIT=PCIT, PCLDFR=PCLDFR,  &
                 &  PICLDFR=PICLDFR, &
                 &  PSSIO=PSSIO, PSSIU=PSSIU, PIFR=PIFR,   &
                 &  PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF, &
                 &  PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF, &
                 &  PTHT=PTHT,PRT=PRT, PTHS=PTHS, PRS=PRS, &
                 &  PINPRC=ZINPRC,PINPRR=PINPRR,PEVAP3D=PEVAP,&
                 &  PINPRS=PINPRS, PINPRG=PINPRG, PINDEP=ZINDEP, PRAINFR=ZRAINFR, &
                 &  PSIGS=PSIGS, &
                 &  TBUDGETS=YLBUDGET_PTR, KBUDGETS=SIZE(YLBUDGET_PTR), &
                 &  PSEA=PSEA, PTOWN=PTOWN, PCONC3D=PCLDROP, &
                 &  PINPRH=PINPRH, PFPR=PFPR)
ELSEIF (CMICRO=='ICE3' .AND. PHYEX%PARAM_ICEN%LRED) THEN
    CALL RAIN_ICE(  YLDIMPHYEX, PHYEX%CST, PHYEX%PARAM_ICEN, PHYEX%RAIN_ICE_PARAMN, &
                 &  PHYEX%RAIN_ICE_DESCRN,  PHYEX%ELEC_PARAM, PHYEX%ELEC_DESCR, &
                 &  PHYEX%MISC%TBUCONF, OELEC=PHYEX%MISC%OELEC, OSEDIM_BEARD=PHYEX%MISC%OSEDIM_BEARD, &
                 &  PTHVREFZIKB=ZTHVREFZIKB, PTSTEP=ZTWOTSTEP, &
                 &  KRR=KRR, PEXN=PEXNREF,            &
                 &  PDZZ=PDZZ, PRHODJ=PRHODJ, PRHODREF=PRHODREF,PEXNREF=PEXNREF,&
                 &  PPABST=PPABSM, PCIT=PCIT, PCLDFR=PCLDFR,  &
                 &  PICLDFR=PICLDFR, &
                 &  PSSIO=PSSIO, PSSIU=PSSIU, PIFR=PIFR,   &
                 &  PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF, &
                 &  PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF, &
                 &  PTHT=PTHT,PRT=PRT, PTHS=PTHS, PRS=PRS, &
                 &  PINPRC=ZINPRC,PINPRR=PINPRR,PEVAP3D=PEVAP,&
                 &  PINPRS=PINPRS, PINPRG=PINPRG, PINDEP=ZINDEP, PRAINFR=ZRAINFR, &
                 &  PSIGS=PSIGS, &
                 &  TBUDGETS=YLBUDGET_PTR, KBUDGETS=SIZE(YLBUDGET_PTR), &
                 &  PSEA=PSEA, PTOWN=PTOWN, PCONC3D=PCLDROP, PFPR=PFPR)
ELSEIF (CMICRO=='ICE4' .AND. .NOT. PHYEX%PARAM_ICEN%LRED) THEN
    IF (YSPP_ICENU%LPERT) THEN
     CALL APPLY_SPP(YSPP_ICENU,KLON,1,KLON,PHYEX%RAIN_ICE_PARAMN%XFRMIN(9),ZICENU)
    ELSE
     ZICENU(:) = PHYEX%RAIN_ICE_PARAMN%XFRMIN(9)
    ENDIF
   
    IF (YSPP_KGN_ACON%LPERT) THEN
     CALL APPLY_SPP(YSPP_KGN_ACON,KLON,1,KLON,PHYEX%RAIN_ICE_PARAMN%XFRMIN(10),ZKGN_ACON)
    ELSE
     ZKGN_ACON(:) = PHYEX%RAIN_ICE_PARAMN%XFRMIN(10)
    ENDIF
   
    IF (YSPP_KGN_SBGR%LPERT) THEN
     CALL APPLY_SPP(YSPP_KGN_SBGR,KLON,1,KLON,PHYEX%RAIN_ICE_PARAMN%XFRMIN(11),ZKGN_SBGR)
    ELSE
     ZKGN_SBGR(:) = PHYEX%RAIN_ICE_PARAMN%XFRMIN(11)
    ENDIF

    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        LLMICRO(JLON,JLEV) = PRT(JLON,JLEV,2)>PHYEX%RAIN_ICE_DESCRN%XRTMIN(2) &
                        & .OR. PRT(JLON,JLEV,3)>PHYEX%RAIN_ICE_DESCRN%XRTMIN(3) &
                        & .OR. PRT(JLON,JLEV,4)>PHYEX%RAIN_ICE_DESCRN%XRTMIN(4) &
                        & .OR. PRT(JLON,JLEV,5)>PHYEX%RAIN_ICE_DESCRN%XRTMIN(5) &
                        & .OR. PRT(JLON,JLEV,6)>PHYEX%RAIN_ICE_DESCRN%XRTMIN(6) &
                        & .OR. PRT(JLON,JLEV,7)>PHYEX%RAIN_ICE_DESCRN%XRTMIN(7)
      ENDDO
    ENDDO

    IF(PHYEX%PARAM_ICEN%LOCND2) THEN
      DO JLEV = 1,KLEV
        DO JLON = KIDIA,KFDIA
          LLMICRO(JLON,JLEV) = LLMICRO(JLON, JLEV) &
                          & .OR. PSSIO(JLON,JLEV) > PHYEX%RAIN_ICE_PARAMN%XFRMIN(12)
        ENDDO
      ENDDO
    ENDIF
    ISIZE=COUNT(LLMICRO)
    CALL RAIN_ICE_OLD(YLDIMPHYEX, PHYEX%CST, PHYEX%PARAM_ICEN, PHYEX%RAIN_ICE_PARAMN, &
                 &  PHYEX%RAIN_ICE_DESCRN, PHYEX%MISC%TBUCONF, TLES, &
                 &  OSEDIC=PHYEX%PARAM_ICEN%LSEDIC, OCND2=PHYEX%PARAM_ICEN%LOCND2, LKOGAN=PHYEX%PARAM_ICEN%LKOGAN, &
                 &  LMODICEDEP=PHYEX%PARAM_ICEN%LMODICEDEP, &
                 &  HSEDIM=PHYEX%PARAM_ICEN%CSEDIM, HSUBG_AUCV_RC=PHYEX%PARAM_ICEN%CSUBG_AUCV_RC, &
                 &  OWARM=PHYEX%PARAM_ICEN%LWARM,KKA=KKA,KKU=KKU,KKL=KKL,KSPLITR=PHYEX%CLOUDPARN%NSPLITR, &
                 &  PTSTEP=ZTWOTSTEP, KRR=KRR, KSIZE=ISIZE, GMICRO=LLMICRO, &
                 &  PDZZ=PDZZ, PRHODJ=PRHODJ, PRHODREF=PRHODREF, PEXNREF=PEXNREF,&
                 &  PPABST=PPABSM, PCIT=PCIT, PCLDFR=PCLDFR,  &
                 &  PICLDFR=PICLDFR, &
                 &  PSSIO=PSSIO, PSSIU=PSSIU, PIFR=PIFR, &
                 &  PTHT=PTHT,PRVT= PRT(:,:,1),PRCT= PRT(:,:,2), &
                 &  PRRT=PRT(:,:,3), &
                 &  PRIT=PRT(:,:,4), PRST=PRT(:,:,5), &
                 &  PRGT=PRT(:,:,6),       &
                 &  PTHS=PTHS, PRVS=PRS(:,:,1),PRCS=PRS(:,:,2),&
                 &  PRRS=PRS(:,:,3),&
                 &  PRIS=PRS(:,:,4),PRSS= PRS(:,:,5),PRGS= PRS(:,:,6),&
                 &  PINPRC=ZINPRC,PINPRR=PINPRR,PEVAP3D=PEVAP,&
                 &  PINPRS=PINPRS, PINPRG=PINPRG, &
                 &  PSIGS=PSIGS, PSEA=PSEA, PTOWN=PTOWN, &
                 &  OAERONRT=OAERONRT, OAEIFN=OAEIFN, PCLDROP=PCLDROP, PIFNNC=PIFNNC, &
                 &  TBUDGETS=YLBUDGET_PTR, KBUDGETS=SIZE(YLBUDGET_PTR), &
                 &  PRHT=PRT(:,:,7),&
                 &  PRHS=PRS(:,:,7), PINPRH=PINPRH, PFPR=PFPR, &
                 &  PICENU=ZICENU, &
                 &  PKGN_ACON=ZKGN_ACON, &
                 &  PKGN_SBGR=ZKGN_SBGR)
ELSE
    IF (YSPP_ICENU%LPERT) THEN
     CALL APPLY_SPP(YSPP_ICENU,KLON,1,KLON,PHYEX%RAIN_ICE_PARAMN%XFRMIN(9),ZICENU)
    ELSE
     ZICENU(:) = PHYEX%RAIN_ICE_PARAMN%XFRMIN(9)
    ENDIF
   
    IF (YSPP_KGN_ACON%LPERT) THEN
     CALL APPLY_SPP(YSPP_KGN_ACON,KLON,1,KLON,PHYEX%RAIN_ICE_PARAMN%XFRMIN(10),ZKGN_ACON)
    ELSE
     ZKGN_ACON(:) = PHYEX%RAIN_ICE_PARAMN%XFRMIN(10)
    ENDIF
   
    IF (YSPP_KGN_SBGR%LPERT) THEN
     CALL APPLY_SPP(YSPP_KGN_SBGR,KLON,1,KLON,PHYEX%RAIN_ICE_PARAMN%XFRMIN(11),ZKGN_SBGR)
    ELSE
     ZKGN_SBGR(:) = PHYEX%RAIN_ICE_PARAMN%XFRMIN(11)
    ENDIF
    IF(PHYEX%PARAM_ICEN%LOCND2) THEN
      DO JLEV = 1,KLEV
        DO JLON = KIDIA,KFDIA
          LLMICRO(JLON,JLEV) = PSSIO(JLON,JLEV)>PHYEX%RAIN_ICE_PARAMN%XFRMIN(12) .OR. &
                         PRT(JLON,JLEV,2)>PHYEX%RAIN_ICE_PARAMN%XFRMIN(13) .OR. &
                         PRT(JLON,JLEV,3)>PHYEX%RAIN_ICE_PARAMN%XFRMIN(13) .OR. &
                         PRT(JLON,JLEV,4)>PHYEX%RAIN_ICE_PARAMN%XFRMIN(13) .OR. &
                         PRT(JLON,JLEV,5)>PHYEX%RAIN_ICE_PARAMN%XFRMIN(13) .OR. &
                         PRT(JLON,JLEV,6)>PHYEX%RAIN_ICE_PARAMN%XFRMIN(13)
        END DO
      END DO
    ELSE
      DO JLEV = 1,KLEV
        DO JLON = KIDIA,KFDIA
          LLMICRO(JLON,JLEV)=PRT(JLON,JLEV,2)>PHYEX%RAIN_ICE_DESCRN%XRTMIN(2) .OR. &                                                                     
                         PRT(JLON,JLEV,3)>PHYEX%RAIN_ICE_DESCRN%XRTMIN(3) .OR. &                                                                     
                         PRT(JLON,JLEV,4)>PHYEX%RAIN_ICE_DESCRN%XRTMIN(4) .OR. &                                                                     
                         PRT(JLON,JLEV,5)>PHYEX%RAIN_ICE_DESCRN%XRTMIN(5) .OR. &                                                                     
                         PRT(JLON,JLEV,6)>PHYEX%RAIN_ICE_DESCRN%XRTMIN(6)
        END DO
      END DO
    ENDIF
    ISIZE=COUNT(LLMICRO)
    CALL RAIN_ICE_OLD(YLDIMPHYEX, PHYEX%CST, PHYEX%PARAM_ICEN, PHYEX%RAIN_ICE_PARAMN, &
                 &  PHYEX%RAIN_ICE_DESCRN, PHYEX%MISC%TBUCONF, TLES, &
                 &  OSEDIC=PHYEX%PARAM_ICEN%LSEDIC, OCND2=PHYEX%PARAM_ICEN%LOCND2, &
                 &  LKOGAN=PHYEX%PARAM_ICEN%LKOGAN, LMODICEDEP=PHYEX%PARAM_ICEN%LMODICEDEP, &
                 &  HSEDIM=PHYEX%PARAM_ICEN%CSEDIM, HSUBG_AUCV_RC=PHYEX%PARAM_ICEN%CSUBG_AUCV_RC, &
                 &  OWARM=PHYEX%PARAM_ICEN%LWARM,KKA=KKA,KKU=KKU,KKL=KKL,KSPLITR=PHYEX%CLOUDPARN%NSPLITR, &
                 &  PTSTEP=ZTWOTSTEP, KRR=KRR, KSIZE=ISIZE, GMICRO=LLMICRO, &
                 &  PDZZ=PDZZ, PRHODJ=PRHODJ, PRHODREF=PRHODREF, PEXNREF=PEXNREF,&
                 &  PPABST=PPABSM, PCIT=PCIT, PCLDFR=PCLDFR,  &
                 &  PICLDFR=PICLDFR, &
                 &  PSSIO=PSSIO, PSSIU=PSSIU, PIFR=PIFR, &
                 &  PTHT=PTHT,PRVT= PRT(:,:,1),PRCT= PRT(:,:,2), &
                 &  PRRT=PRT(:,:,3), &
                 &  PRIT=PRT(:,:,4), PRST=PRT(:,:,5), &
                 &  PRGT=PRT(:,:,6),       &
                 &  PTHS=PTHS, PRVS=PRS(:,:,1),PRCS=PRS(:,:,2),&
                 &  PRRS=PRS(:,:,3),&
                 &  PRIS=PRS(:,:,4),PRSS= PRS(:,:,5),PRGS= PRS(:,:,6),&
                 &  PINPRC=ZINPRC,PINPRR=PINPRR,PEVAP3D=PEVAP,&
                 &  PINPRS=PINPRS, PINPRG=PINPRG, &
                 &  PSIGS=PSIGS, PSEA=PSEA, PTOWN=PTOWN, &
                 &  OAERONRT=OAERONRT, OAEIFN=OAEIFN, PCLDROP=PCLDROP, PIFNNC=PIFNNC, &
                 &  TBUDGETS=YLBUDGET_PTR, KBUDGETS=SIZE(YLBUDGET_PTR), &
                 &  PFPR=PFPR, &
                 &  PICENU=ZICENU, &
                 &  PKGN_ACON=ZKGN_ACON, &
                 &  PKGN_SBGR=ZKGN_SBGR)
ENDIF

!add ZINPRC in PINPRR
DO JLON = KIDIA,KFDIA
  PINPRR(JLON) = PINPRR(JLON) + ZINPRC(JLON)
END DO
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_RAIN_ICE',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_RAIN_ICE
