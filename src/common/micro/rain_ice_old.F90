!     ######spl
      SUBROUTINE RAIN_ICE_OLD (D, CST, PARAMI, ICEP, ICED, BUCONF,                    &
                               OSEDIC, OCND2, LKOGAN, LMODICEDEP,                     &
                               HSEDIM, HSUBG_AUCV_RC, OWARM,                          &
                               KKA, KKU, KKL,                                         &
                               KSPLITR, PTSTEP, KRR, KSIZE, GMICRO,                   &
                               PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                               PICLDFR, PSSIO, PSSIU, PIFR,                           &
                               PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                    &
                               PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,        &
                               PINPRC, PINPRR, PEVAP3D,                               &
                               PINPRS, PINPRG, PSIGS, PSEA, PTOWN,                    &
                               TBUDGETS, KBUDGETS,                                    &
                               PICENU, PKGN_ACON, PKGN_SBGR,                          &
                               PRHT, PRHS, PINPRH, PFPR)

      USE PARKIND1,            ONLY: JPRB
      USE YOMHOOK,             ONLY: LHOOK, DR_HOOK, JPHOOK
      USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_T
      USE MODD_CST,            ONLY: CST_T
      USE MODD_PARAM_ICE_N,    ONLY: PARAM_ICE_t
      USE MODD_RAIN_ICE_PARAM_N, ONLY: RAIN_ICE_PARAM_T
      USE MODD_RAIN_ICE_DESCR_N, ONLY: RAIN_ICE_DESCR_T
!     ######################################################################
!
!!****  * -  compute the explicit microphysical sources
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the slow microphysical sources
!!    which can be computed explicitly
!!
!!
!!**  METHOD
!!    ------
!!      The autoconversion computation follows Kessler (1969).
!!      The sedimentation rate is computed with a time spliting technique and
!!    an upstream scheme, written as a difference of non-advective fluxes. This
!!    source term is added to the future instant ( split-implicit process ).
!!      The others microphysical processes are evaluated at the central instant
!!    (split-explicit process ): autoconversion, accretion and rain evaporation.
!!      These last 3 terms are bounded in order not to create negative values
!!    for the water species at the future instant.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!          JPHEXT       : Horizontal external points number
!!          JPVEXT       : Vertical external points number
!!      Module MODD_CONF :
!!          CCONF configuration of the model for the first time step
!!      Module MODD_CST
!!          XP00               ! Reference pressure
!!          XRD,XRV            ! Gaz  constant for dry air, vapor
!!          XMD,XMV            ! Molecular weight for dry air, vapor
!!          XCPD               ! Cpd (dry air)
!!          XCL                ! Cl (liquid)
!!          XCI                ! Ci (solid)
!!          XTT                ! Triple point temperature
!!          XLVTT              ! Vaporization heat constant
!!          XALPW,XBETAW,XGAMW ! Constants for saturation vapor pressure
!!                               function over liquid water
!!          XALPI,XBETAI,XGAMI ! Constants for saturation vapor pressure
!!                               function over solid ice
!!      Module MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask
!!                          'NONE'  ' for no budget
!!         NBUPROCCTR   : process counter used for each budget variable
!!         LBU_RTH      : logical for budget of RTH (potential temperature)
!!                        .TRUE. = budget of RTH
!!                        .FALSE. = no budget of RTH
!!         LBU_RRV      : logical for budget of RRV (water vapor)
!!                        .TRUE. = budget of RRV
!!                        .FALSE. = no budget of RRV
!!         LBU_RRC      : logical for budget of RRC (cloud water)
!!                        .TRUE. = budget of RRC
!!                        .FALSE. = no budget of RRC
!!         LBU_RRI      : logical for budget of RRI (cloud ice)
!!                        .TRUE. = budget of RRI
!!                        .FALSE. = no budget of RRI
!!         LBU_RRR      : logical for budget of RRR (rain water)
!!                        .TRUE. = budget of RRR
!!                        .FALSE. = no budget of RRR
!!         LBU_RRS      : logical for budget of RRS (aggregates)
!!                        .TRUE. = budget of RRS
!!                        .FALSE. = no budget of RRS
!!         LBU_RRG      : logical for budget of RRG (graupeln)
!!                        .TRUE. = budget of RRG
!!                        .FALSE. = no budget of RRG
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1 and Book2 of documentation ( routine RAIN_ICE )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/11/95
!!      (J.Viviand) 04/02/97  debug accumulated prcipitation & convert
!!                            precipitation rate in m/s
!!      (J.-P. Pinty) 17/02/97  add budget calls
!!      (J.-P. Pinty) 17/11/97  set ice sedim. for cirrus ice, reset RCHONI
!!                              and RRHONG, reverse order for DEALLOCATE
!!      (J.-P. Pinty) 11/02/98  correction of the air dynamical viscosity and
!!                              add advance of the budget calls
!!      (J.-P. Pinty) 18/05/98  correction of the air density in the RIAUTS
!!                              process
!!      (J.-P. Pinty) 18/11/98  split the main routine
!!      (V. Masson)   18/11/98  bug in IVEC1 and IVEC2 upper limits
!!      (J. Escobar & J.-P. Pinty)
!!                    11/12/98  contains and rewrite count+pack
!!      (J. Stein & J.-P. Pinty)
!!                    14/10/99  correction for very small RIT
!!      (J. Escobar & J.-P. Pinty)
!!                    24/07/00  correction for very samll m.r. in
!!                              the sedimentation subroutine
!!      (M. Tomasini) 11/05/01  Autoconversion of rc into rr modification to take
!!                              into account the subgrid variance
!!                              (cf Redelsperger & Sommeria JAS 86)
!!      (G. Molinie)  21/05/99  bug in RRCFRIG process, RHODREF**(-1) missing
!!                              in RSRIMCG
!!      (G. Molinie & J.-P. Pinty)
!!                    21/06/99  bug in RACCS process
!!      (P. Jabouille) 27/05/04 safety test for case where esw/i(T)> pabs (~Z>40km)
!!      (J-.P. Chaboureau) 12/02/05  temperature depending ice-to-snow autocon-
!                              version threshold (Chaboureau and Pinty GRL 2006)
!!      (J.-P. Pinty) 01/01/O1  add the hail category and correction of the
!!                              wet growth rate of the graupeln
!!      (S.Remy & C.Lac) 06/06 Add the cloud sedimentation
!!      (S.Remy & C.Lac) 06/06 Sedimentation becoming the last process
!!      to settle the precipitating species created during the current time step
!!      (S.Remy & C.Lac) 06/06 Modification of the algorithm of sedimentation
!!      to settle n times the precipitating species created during Dt/n instead
!!      of Dt
!!      (C.Lac) 11/06 Optimization of the sedimentation loop for NEC
!!      (J.Escobar) 18/01/2008 Parallel Bug in Budget when IMICRO >= 1
!!                  --> Path inhibit this test by IMICRO >= 0 allway true
!!      (Y.Seity) 03/2008 Add Statistic sedimentation
!!      (Y.Seity) 10/2009 Added condition for the raindrop accretion of the aggregates
!!         into graupeln process (5.2.6) to avoid negative graupel mixing ratio
!!      (V.Masson, C.Lac) 09/2010 Correction in split sedimentation for
!!                                reproducibility
!!      (S. Riette) Oct 2010 Better vectorisation of RAIN_ICE_SEDIMENTATION_STAT
!!      (Y. Seity), 02-2012  add possibility to run with reversed vertical levels
!!      (L. Bengtsson), 02-2013 Passing in land/sea mask and town fraction in
!!                      order to use different cloud droplet number conc. over
!!                      land, sea and urban areas in the cloud sedimentation.
!!      (D. Degrauwe), 2013-11: Export upper-air precipitation fluxes PFPR.
!!      (S. Riette) Nov 2013 Protection against null sigma
!!      (C. Abiven, Y. Léauté, V. Seigner, S. Riette) Phasing of Turner rain subgrid param
!!      (K.I Ivarsson 2014) OCND2-option, possible to use derived cloud dropelt conc for cloudphysics,
!!                      replace XMV/XMD with XEPSILO
!!      (K.I Ivarsson 2016) LKOGAN-option, possible to use Kogan autoconversion of liquid regardless of OCND2 option.
!!      (K.I Ivarsson 2018-02 New varibles for input/ output mainly for optimation. Some updates of OCND2 option.
!!                          Sedimented ice should be preciptation
!!      (U. Andrae Dec 2020) Introduce SPP for HARMONIE-AROME
!!      (C. Wittmann Jan 2021) Introduce sublimation factor tuning
!
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPVEXT
USE MODD_BUDGET,     ONLY: TBUDGETDATA, TBUDGETCONF_t, NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, &
                           NBUDGET_RI, NBUDGET_RR, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, &
                           LBU_ENABLE
USE MODD_LES,        ONLY: TLES
USE MODE_BUDGET_PHY, ONLY: BUDGET_STORE_ADD_PHY, BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY
USE MODI_GAMMA,      ONLY: GAMMA
USE MODE_TIWMX,      ONLY: ESATI, ESATW, AA2, BB3, AA2W, BB3W
USE MODE_TIWMX_TAB,  ONLY: TIWMX_TAB
!
USE MODE_RAIN_ICE_OLD_NUCLEATION,          ONLY: RAIN_ICE_OLD_NUCLEATION
USE MODE_RAIN_ICE_OLD_SEDIMENTATION_STAT,  ONLY: RAIN_ICE_OLD_SEDIMENTATION_STAT
USE MODE_RAIN_ICE_OLD_SEDIMENTATION_SPLIT, ONLY: RAIN_ICE_OLD_SEDIMENTATION_SPLIT
USE MODE_RAIN_ICE_OLD_SLOW,                ONLY: RAIN_ICE_OLD_SLOW
USE MODE_RAIN_ICE_OLD_WARM,                ONLY: RAIN_ICE_OLD_WARM
USE MODE_RAIN_ICE_OLD_FAST_RS,             ONLY: RAIN_ICE_OLD_FAST_RS
USE MODE_RAIN_ICE_OLD_FAST_RG,             ONLY: RAIN_ICE_OLD_FAST_RG
USE MODE_RAIN_ICE_OLD_FAST_RH,             ONLY: RAIN_ICE_OLD_FAST_RH
USE MODE_RAIN_ICE_OLD_FAST_RI,             ONLY: RAIN_ICE_OLD_FAST_RI


IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_T),       INTENT(IN) :: D
TYPE(CST_T),            INTENT(IN) :: CST 
TYPE(PARAM_ICE_t),      INTENT(IN) :: PARAMI
TYPE(RAIN_ICE_PARAM_t), INTENT(IN) :: ICEP
TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF

LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
LOGICAL,                  INTENT(IN)    :: OCND2  ! Logical switch to separate liquid and ice
LOGICAL,                  INTENT(IN)    :: LKOGAN ! Logical switch for using Kogan autoconversion of liquid.
LOGICAL,                  INTENT(IN)    :: LMODICEDEP ! Logical switch for alternative dep/evap of ice
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Switch for rc->rr Subgrid autoconversion
                                        ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step
                                      ! integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
!
LOGICAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: GMICRO    ! Layer thickness (m)

INTEGER, INTENT(IN) :: KSIZE
!
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PDZZ    ! Layer thickness (m)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PCLDFR  ! Cloud fraction
!
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PSIGS   ! Sigma_s at t
! input from aro_adjust / condensation with OCND2, dummy if OCND2 = F
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PICLDFR ! ice cloud fraction
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PSSIO   ! Super-saturation with respect to ice in the
                                                 ! supersaturated fraction
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PSSIU   ! Sub-saturation with respect to ice in the
                                                 ! subsaturated fraction
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PIFR    ! Ratio cloud ice moist part to dry part
! input from aro_adjust / condensation with OCND2 END.
!
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRC! Cloud instant precip
REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRR! Rain instant precip
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PEVAP3D! Rain evap profile
REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRS! Snow instant precip
REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRG! Graupel instant precip
REAL, DIMENSION(D%NIT),       INTENT(IN)  :: PSEA ! Sea Mask
REAL, DIMENSION(D%NIT),       INTENT(IN)  :: PTOWN! Fraction that is town
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
REAL, DIMENSION(D%NIT), INTENT(IN)            :: PICENU, PKGN_ACON, PKGN_SBGR
REAL, DIMENSION(D%NIT,D%NKT),   OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(D%NIT,D%NKT),   OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(D%NIT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(D%NIT,D%NKT,KRR), OPTIONAL, INTENT(OUT) :: PFPR    ! upper-air precipitation fluxes

!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK            ! Vertical loop index for the rain sedimentation
INTEGER :: JI            ! Loop index for the interpolation
INTEGER :: IKB           !
INTEGER :: IKE           !
!
INTEGER :: IMICRO ! Case number of sedimentation, T>0 (for HEN) and r_x>0 locations
REAL, DIMENSION(D%NIT,D%NKT) :: ZW        ! work array
REAL, DIMENSION(D%NIT)       :: ZCONC_TMP ! Weighted concentration
REAL, DIMENSION(D%NIT,D%NKT) :: ZT        ! Temperature
REAL, DIMENSION(D%NIT,D%NKT) :: ZRAY      ! Cloud Mean radius
REAL, DIMENSION(D%NIT,D%NKT) :: ZLBC      ! XLBC weighted by sea fraction
REAL, DIMENSION(D%NIT,D%NKT) :: ZFSEDC
REAL, DIMENSION(D%NIT,D%NKT) :: ZCONC3D   ! droplet concentration m-3
REAL, DIMENSION(D%NIT,D%NKT) :: ZZZZ      ! geometric height
REAL, DIMENSION(D%NIT,D%NKT) :: ZZZT      ! tempoary value for geometric height

!Diagnostics
REAL, DIMENSION(D%NIT,D%NKT) :: ZRAINFR

REAL, DIMENSION(KSIZE) :: ZRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE) :: ZRCT    ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE) :: ZRRT    ! Rain water m.r. at t
REAL, DIMENSION(KSIZE) :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE) :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE) :: ZRGT    ! Graupel m.r. at t
REAL, DIMENSION(KSIZE) :: ZRHT    ! Hail m.r. at t
REAL, DIMENSION(KSIZE) :: ZCIT    ! Pristine ice conc. at t
!
REAL, DIMENSION(KSIZE) :: ZRVS    ! Water vapor m.r. source
REAL, DIMENSION(KSIZE) :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(KSIZE) :: ZRRS    ! Rain water m.r. source
REAL, DIMENSION(KSIZE) :: ZRIS    ! Pristine ice m.r. source
REAL, DIMENSION(KSIZE) :: ZRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(KSIZE) :: ZRGS    ! Graupel m.r. source
REAL, DIMENSION(KSIZE) :: ZRHS    ! Hail m.r. source
REAL, DIMENSION(KSIZE) :: ZTHS    ! Theta source
REAL, DIMENSION(KSIZE) :: ZTHT    ! Potential temperature
REAL, DIMENSION(KSIZE) :: ZTHLT   ! Liquid potential temperature
!
REAL, DIMENSION(KSIZE) :: ZRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(KSIZE) :: ZRHODJ    ! RHO times Jacobian
REAL, DIMENSION(KSIZE) :: ZEXNREF   ! EXNer Pressure REFerence
REAL, DIMENSION(KSIZE) :: ZLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE) :: ZLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE) :: ZLBDAR    ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KSIZE) :: ZLBDAR_RF ! Slope parameter of the raindrop  distribution
                                    ! for the Rain Fraction part
REAL, DIMENSION(KSIZE) :: ZLBDAS    ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE) :: ZLBDAG    ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KSIZE) :: ZLBDAH    ! Slope parameter of the hail      distribution
REAL, DIMENSION(KSIZE) :: ZAI       ! Thermodynamical function
REAL, DIMENSION(KSIZE) :: ZCJ       ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE) :: ZKA       ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE) :: ZDV       ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE) :: ZSIGMA_RC ! Standard deviation of rc at time t
REAL, DIMENSION(KSIZE) :: ZCF       ! Cloud fraction
REAL, DIMENSION(KSIZE) :: ZRF       ! Rain fraction

!                          *******  used for logical switch OCND2 : *******
REAL, DIMENSION(KSIZE) :: ZSSIO      ! Super-saturation with respect to ice in the
                                     ! supersaturated fraction of gridbox
REAL, DIMENSION(KSIZE) :: ZSSIU      ! Sub-saturation with respect to ice in the
                                     ! sub-saturated fraction of gridbox
REAL, DIMENSION(KSIZE) :: ZW2D       ! Factor for subgridscale calculations
REAL, DIMENSION(KSIZE) :: ZXW2D      ! Ratio cloud ice moist part to dry part
REAL, DIMENSION(KSIZE) :: ZXW2D13    ! ZXW2D**0.333 or other expression for LMODICEDEP=T
REAL, DIMENSION(KSIZE) :: ZCOLF      ! collision factor cloud liquid to snow / graupel
REAL, DIMENSION(KSIZE) :: ZACRF      ! collision factor cloud liquid to rain
REAL, DIMENSION(KSIZE) :: ZCONCM     ! Same as ZCONC3D but GMICRO-array only and cm-3 instead of m-3

REAL, DIMENSION(KSIZE) :: ZHLC_HCF   ! HLCLOUDS : fraction of High Cloud Fraction in grid
REAL, DIMENSION(KSIZE) :: ZHLC_LCF   ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                     !    note that ZCF = ZHLC_HCF + ZHLC_LCF
REAL, DIMENSION(KSIZE) :: ZHLC_HRC   ! HLCLOUDS : LWC that is High LWC in grid
REAL, DIMENSION(KSIZE) :: ZHLC_LRC   ! HLCLOUDS : LWC that is Low  LWC in grid
                                     !    note that ZRC = ZHLC_HRC + ZHLC_LRC
REAL, DIMENSION(KSIZE) :: ZHLC_RCMAX ! HLCLOUDS : maximum value for RC in distribution
REAL, DIMENSION(KSIZE) :: ZRCRAUTC   ! RC value to begin rain formation =XCRIAUTC/RHODREF
!                          *******  end logical switch OCND2 *******

REAL, DIMENSION(KSIZE) :: ZHLC_HRCLOCAL ! HLCLOUDS : LWC that is High LWC local in HCF
REAL, DIMENSION(KSIZE) :: ZHLC_LRCLOCAL ! HLCLOUDS : LWC that is Low  LWC local in LCF
                                        !    note that ZRC/CF = ZHLC_HRCLOCAL+ ZHLC_LRCLOCAL
                                        !                     = ZHLC_HRC/HCF+ ZHLC_LRC/LCF

REAL, DIMENSION(KSIZE) :: ZAA2    ! Part of ZAI used for optimized code
REAL, DIMENSION(KSIZE) :: ZBB3    ! Part of ZAI used for optimized code
REAL, DIMENSION(KSIZE) :: ZAA2W   ! as ZAA2 but for liquid
REAL, DIMENSION(KSIZE) :: ZBB3W   ! as ZBB3 but for liquid
REAL, DIMENSION(KSIZE) :: ZTIW    ! Wet bulb temperature

REAL, DIMENSION(KSIZE) :: ZZT   ! Temperature
REAL, DIMENSION(KSIZE) :: ZPRES ! Pressure
REAL, DIMENSION(KSIZE) :: ZZW   ! Work array
REAL, DIMENSION(KSIZE) :: ZSSI  ! Supersaturation over ice

!                          *******  used for logical switch OCND2 : *******
REAL, DIMENSION(KSIZE) :: ZSIFRC ! subgridscale fraction with supersaturation with
REAL, DIMENSION(KSIZE) :: ZESI   ! saturation pressure over ice
REAL, DIMENSION(KSIZE) :: ZESW   ! saturation pressure over water
!                          *******  end logical switch OCND2 *******

!    *******  used for logical switch OCND2 : *******
REAL            :: ZRVSOLD, ZTSP
REAL            :: ZRSP, ZRISOLD, ZRSSOLD, ZRGSOLD ! Function,old ice
REAL            :: ZRISFRC, ZRSSFRC, ZRGSFRC, ZRFRAC, ZRSA, ZRSTS, ZRSB, ZRSDIF, ZR20
REAL            :: ZRSI, ZRSW, ZDICRIT, ZCITRED23, ZCITRED, ZRCW, ZVT, ZST
REAL            :: ZREDGR, ZREDSN    ! Possible reduction of the rate of graupel,snow growth
REAL            :: ZKVO  ! factor used for caluclate maximum mass in the ice
                         ! distubution.
!    *******  end logical switch OCND2 *******

! SPP arrays
REAL, DIMENSION(KSIZE) :: ZZKGN_ACON,ZZKGN_SBGR

     !internal fractions etc, finally saturation ratio over ice 'source' value

INTEGER, DIMENSION(D%NIT*D%NKT) :: I1, I2 ! Used to replace the COUNT
INTEGER                         :: JL       ! and PACK intrinsics
REAL :: ZCOEFFRCM
LOGICAL :: LCHECKNOISE ! Noise check on/off
LOGICAL :: LTIW   ! Use TIW for graupel melting ( set by XFRMIN(18) ~ 1)
!
REAL, DIMENSION(D%NIT,D%NKT) :: ZBU0
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.1     COMPUTE THE LOOP BOUNDS
!               -----------------------
!
IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD',0,ZHOOK_HANDLE)
LCHECKNOISE=.TRUE.
IKB=KKA+JPVEXT*KKL
IKE=KKU-JPVEXT*KKL

!
!*       1.2     COMPUTE SOME CONSTANT PARAMETERS
!

ZCITRED = 0.1     ! ratio of ice crystal concentration wet to dry
                  ! part of a gridbox
ZDICRIT = ICEP%XFRMIN(15)  ! Critical diameter of ice crystal to define what
                      ! is cloud ice and what is snow (m)

ZCITRED23 = ZCITRED**0.667
IF (LMODICEDEP) THEN
  ZCITRED = 1.
  !from spherical diameter to max diameter
  ZDICRIT = (700.*CST%XPI/ICED%XAI/6.)**(1./ICED%XBI)*ZDICRIT**(3./ICED%XBI) 
  ZCITRED23 = ZCITRED**(1.+ ICED%XLBEXI)
  ZKVO      = ((ICED%XALPHAI*ICED%XNUI + ICED%XBI -1.)/ICED%XALPHAI)**(1./ICED%XALPHAI)
  ZKVO =  ZKVO/ZDICRIT/ICEP%XFRMIN(14)
  PIFR = 1.
ENDIF

ZREDGR  = 1.      ! Tuning of the deposition of graupel, 1. is ref. value

ZREDSN  = 1.      ! Tuning of the deposition of snow, 1. is ref. value

IF(OCND2) THEN
   IF (.NOT. LMODICEDEP) THEN
      ZREDGR  = ICEP%XFRMIN(39)  ! Tuning factor, may be /= 1.
      ZREDSN  = ICEP%XFRMIN(40)  ! Tuning factor, may be /= 1.
   ENDIF
ENDIF

LTIW=.FALSE.
IF (NINT(ICEP%XFRMIN(18)) == 1) LTIW=.TRUE.

!
!*       1.3    COMPUTE THE DROPLET NUMBER CONCENTRATION
!   	        ----------------------------------------
!         (Do it already here, since also used with OCND2=T )

IF (OSEDIC.OR.OCND2) THEN
   ZRAY(:,:)   = 0.
   ZZZZ(:,D%NKTE)   = PDZZ(:,D%NKTE)*0.5
   ZZZT(:,D%NKTE)   = PDZZ(:,D%NKTE)
   IF (ICEP%XFRMIN(26)>0.001) THEN ! Use alternative concentration given by (XFRMIN(26)
      ZCONC_TMP(:) = ICEP%XFRMIN(26)
      DO JK=D%NKTB,D%NKTE
         ZLBC(:,JK)   = 0.5* (ICED%XLBC(2)+ICED%XLBC(1)) ! Assume "average" distr. func for simplicity
         ZFSEDC(:,JK) = 0.5* (ICEP%XFSEDC(2)+ICEP%XFSEDC(1))
         ZFSEDC(:,JK) = MAX(MIN(ICEP%XFSEDC(1),ICEP%XFSEDC(2)),ZFSEDC(:,JK))
         ZCONC3D(:,JK)= ZCONC_TMP(:)*PPABST(:,JK)/CST%XP00 ! Let it be diluted with decreasing pressure
         ZRAY(:,JK)   = 0.5*( 0.5*GAMMA(ICED%XNUC+1.0/ICED%XALPHAC)/(GAMMA(ICED%XNUC)) + &
           0.5*GAMMA(ICED%XNUC2+1.0/ICED%XALPHAC2)/(GAMMA(ICED%XNUC2)))
      ENDDO
   ELSE
     ZCONC_TMP(:)=PSEA(:)*ICED%XCONC_SEA+(1.-PSEA(:))*ICED%XCONC_LAND

     DO JK=D%NKTB,D%NKTE
        ZLBC(:,JK)   = PSEA(:)*ICED%XLBC(2)+(1.-PSEA(:))*ICED%XLBC(1)
        ZFSEDC(:,JK) = (PSEA(:)*ICEP%XFSEDC(2)+(1.-PSEA(:))*ICEP%XFSEDC(1))
        ZFSEDC(:,JK) = MAX(MIN(ICEP%XFSEDC(1),ICEP%XFSEDC(2)),ZFSEDC(:,JK))
        ZCONC3D(:,JK)= (1.-PTOWN(:))*ZCONC_TMP(:)+PTOWN(:)*ICED%XCONC_URBAN
        ZRAY(:,JK)   = 0.5*((1.-PSEA(:))*GAMMA(ICED%XNUC+1.0/ICED%XALPHAC)/(GAMMA(ICED%XNUC)) + &
           PSEA(:)*GAMMA(ICED%XNUC2+1.0/ICED%XALPHAC2)/(GAMMA(ICED%XNUC2)))
     ENDDO
   ENDIF

   ZCONC3D(:,D%NKTE)= ZCONC3D(:,D%NKTE)*MAX(0.001,ICEP%XFRMIN(22))
   ZRAY(:,:)      = MAX(1.,ZRAY(:,:))
   ZLBC(:,:)      = MAX(MIN(ICED%XLBC(1),ICED%XLBC(2)),ZLBC(:,:))

   DO JK=D%NKTE-1,D%NKTB,-1
     ZZZT(:,JK) = ZZZT(:,JK+1) + PDZZ(:,JK)
     ZZZZ(:,JK) = ZZZT(:,JK) - 0.5*PDZZ(:,JK)
   ENDDO
ENDIF

ZT(:,:) = PTHT(:,:) * ( PPABST(:,:) / CST%XP00 ) ** (CST%XRD/CST%XCPD)

IF(BUCONF%LBU_ENABLE) THEN
  IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'HENU', PTHS(:,:)*PRHODJ(:,:))
  IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'HENU', PRVS(:,:)*PRHODJ(:,:))
  IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'HENU', PRIS(:,:)*PRHODJ(:,:))
ENDIF
CALL RAIN_ICE_OLD_NUCLEATION(D, CST, ICEP, COUNT(ZT(D%NIB:D%NIE,D%NKTB:D%NKTE)<CST%XTT), &
                             OCND2, LMODICEDEP, KRR, PTSTEP, &
                             PTHT, PPABST, PEXNREF, PICLDFR, PRHODJ, PRHODREF, &
                             PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                             PTHS, PRVS, PRIS, PCIT, &
                             PICENU, ZT, ZZZZ, &
                             PRHT)
IF(BUCONF%LBU_ENABLE) THEN
  IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'HENU', PTHS(:,:)*PRHODJ(:,:))
  IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'HENU', PRVS(:,:)*PRHODJ(:,:))
  IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'HENU', PRIS(:,:)*PRHODJ(:,:))
ENDIF

IMICRO = 0
DO JK = 1, D%NKT
  DO JI = 1, D%NIT
    IF(GMICRO(JI,JK)) THEN
      IMICRO = IMICRO + 1
      I1(IMICRO) = JI
      I2(IMICRO) = JK
    END IF
  END DO
END DO

IF ( KSIZE >= 0 ) THEN

  IF (OCND2) THEN
     IF (LMODICEDEP) THEN
        DO JL=1,KSIZE
           ZXW2D(JL) = PIFR(I1(JL),I2(JL))
           ZXW2D13(JL)=ZXW2D(JL)**(-ICED%XLBEXI)
        ENDDO
     ELSE
        DO JL=1,KSIZE
           ZXW2D(JL) = PIFR(I1(JL),I2(JL))
           ZXW2D13(JL)=ZXW2D(JL)**0.333
        ENDDO
     ENDIF
  ENDIF

  DO JL=1,KSIZE
    ZRVT(JL) = PRVT(I1(JL),I2(JL))
    ZRCT(JL) = PRCT(I1(JL),I2(JL))
    ZRRT(JL) = PRRT(I1(JL),I2(JL))
    ZRIT(JL) = PRIT(I1(JL),I2(JL))
    ZRST(JL) = PRST(I1(JL),I2(JL))
    ZRGT(JL) = PRGT(I1(JL),I2(JL))

    IF ( KRR == 7 ) ZRHT(JL) = PRHT(I1(JL),I2(JL))
    ZCIT(JL) = PCIT(I1(JL),I2(JL))
    ZCF(JL) = PCLDFR(I1(JL),I2(JL))
    IF ( HSUBG_AUCV_RC == 'PDF ' .AND. PARAMI%CSUBG_PR_PDF == 'SIGM' ) THEN
      ZSIGMA_RC(JL) = PSIGS(I1(JL),I2(JL)) * 2.
    END IF

    ZRVS(JL) = PRVS(I1(JL),I2(JL))
    ZRCS(JL) = PRCS(I1(JL),I2(JL))
    ZRRS(JL) = PRRS(I1(JL),I2(JL))
    ZRIS(JL) = PRIS(I1(JL),I2(JL))
    ZRSS(JL) = PRSS(I1(JL),I2(JL))
    ZRGS(JL) = PRGS(I1(JL),I2(JL))
    IF ( KRR == 7 ) ZRHS(JL) = PRHS(I1(JL),I2(JL))
    ZTHS(JL) = PTHS(I1(JL),I2(JL))
!

    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL))
    ZZT(JL) = ZT(I1(JL),I2(JL))
    ZTHT(JL) = PTHT(I1(JL),I2(JL))
    ZTHLT(JL) = ZTHT(JL) - CST%XLVTT * ZTHT(JL) / CST%XCPD / ZZT(JL) * ZRCT(JL)
    ZPRES(JL) = PPABST(I1(JL),I2(JL))
    ZEXNREF(JL) = PEXNREF(I1(JL),I2(JL))
    ZCOLF(JL)=1. ! No change from orignal when  OCND2 = .FALSE.
    ZACRF(JL)=1. !    "      "       "            "       "
    ZCONCM(JL)=ZCONC3D(I1(JL),I2(JL))*0.000001 ! From m-3 to cm-3
    IF (LTIW) ZTIW(JL)=TIWMX_TAB(ZPRES(JL),ZZT(JL), ZRVS(JL)*PTSTEP,0._JPRB,ZRSP,ZRSW,0.1_JPRB)
    ZZKGN_ACON(JL)=PKGN_ACON(I1(JL))
    ZZKGN_SBGR(JL)=PKGN_SBGR(I1(JL))

    IF (OCND2) THEN
       ZESI(JL) = ESATI(ZZT(JL))
       ZESW(JL) = ESATW(ZZT(JL))
       ZAA2(JL) = AA2(ZZT(JL))
       ZBB3(JL) = BB3(ZZT(JL))
       ZAA2W(JL)= AA2W(ZZT(JL))
       ZBB3W(JL)= BB3W(ZZT(JL))
       ZSIFRC(JL) = PICLDFR(I1(JL),I2(JL))
       ZSSIO(JL) = PSSIO(I1(JL),I2(JL))
       ZSSIU(JL) = PSSIU(I1(JL),I2(JL))
       ZW2D(JL) = 1./(ZXW2D(JL)*ZSIFRC(JL) + 1. -ZSIFRC(JL))
       ZCOLF(JL)=0.00001
       ZACRF(JL)=0.00001
       IF(ZRCT(JL)>1.0E-10)THEN
          ! mean cloud droplet radius in cm
          ZRCW =  0.1*(0.75*ZRCT(JL)*ZRHODREF(JL)/(CST%XPI*ZCONCM(JL)))**0.333
          ! fall speed for mean cloud droplet with cloud droplet radius in cm/s
          IF(ZRCW < 0.0065 )THEN
             ZVT   =  1.19E6*ZRCW**2
          ELSE
             ZVT   =  8000.*ZRCW
          ENDIF
          ZVT = MIN(10.,ZVT)
          ZST = MAX(0.01,2.*(100.-ZVT)*ZVT/(CST%XG*10.))
          IF(ZST > 0.1) ZCOLF(JL) =  MAX(0.01,MIN(1.,0.939*ZST**2.657))
          IF( ZRRS(JL) > 1.0E-10 .AND. ZRCW >1.0E-5)THEN
            ZR20 = EXP(ZRCW*2000.)  ! This ZRCW is in cm . To convert to micro meter : x 10000
            ZACRF(JL)  = (ZR20 -1.)/(ZR20 +1.)             ! ZRCW is then multiplied with 0.2
          ENDIF

       ENDIF

    ENDIF

  ENDDO

  ZZW(:)  = ZEXNREF(:)*(CST%XCPD+CST%XCPV*ZRVT(:) + CST%XCL*(ZRCT(:)+ZRRT(:)) &
                                  + CST%XCI*(ZRIT(:)+ZRST(:) + ZRGT(:)) )
  ZLSFACT(:) = (CST%XLSTT + (CST%XCPV - CST%XCI)*(ZZT(:) - CST%XTT))/ZZW(:) ! L_s/(Pi_ref*C_ph)
  ZLVFACT(:) = (CST%XLVTT + (CST%XCPV - CST%XCL)*(ZZT(:) - CST%XTT))/ZZW(:) ! L_v/(Pi_ref*C_ph)

  IF(OCND2)THEN
    ZSSI(:) = ZRVT(:)*( ZPRES(:)-ZESI(:) ) / ( CST%XEPSILO * ZESI(:) ) - 1.0
  ELSE                                                  ! Supersaturation over ice
    ZZW(:) = EXP( CST%XALPI - CST%XBETAI/ZZT(:) - CST%XGAMI*ALOG(ZZT(:) ) )
    ZSSI(:) = ZRVT(:)*( ZPRES(:)-ZZW(:) ) / ( CST%XEPSILO * ZZW(:) ) - 1.0
                                                    ! Supersaturation over ice
  ENDIF

  IF (LBU_ENABLE .OR. TLES%LLES_CALL) THEN

    ZRHODJ(:) = PACK( PRHODJ(:,:),MASK=GMICRO(:,:) )

  END IF
!
  !Cloud water split between high and low content part is done here
  !according to autoconversion option
  ZRCRAUTC(:)   = ICEP%XCRIAUTC/ZRHODREF(:) ! Autoconversion rc threshold
  IF (HSUBG_AUCV_RC == 'NONE') THEN
    !Cloud water is entirely in low or high part
    WHERE (ZRCT(:) > ZRCRAUTC(:))
      ZHLC_HCF(:) = 1.
      ZHLC_LCF(:) = 0.0
      ZHLC_HRC(:) = ZRCT(:)
      ZHLC_LRC(:) = 0.0
      ZRF(:)      = 1.
    ELSEWHERE (ZRCT(:) > ICED%XRTMIN(2))
      ZHLC_HCF(:) = 0.0
      ZHLC_LCF(:) = 1.
      ZHLC_HRC(:) = 0.0
      ZHLC_LRC(:) = ZRCT(:)
      ZRF(:)      = 0.
    ELSEWHERE
      ZHLC_HCF(:) = 0.0
      ZHLC_LCF(:) = 0.0
      ZHLC_HRC(:) = 0.0
      ZHLC_LRC(:) = 0.0
      ZRF(:)      = 0.
    END WHERE

  ELSEIF (HSUBG_AUCV_RC == 'CLFR') THEN
    !Cloud water is only in the cloudy part and entirely in low or high part
    WHERE (ZCF(:) > 0. )
      WHERE (ZRCT(:)/ZCF(:) > ZRCRAUTC(:))
        ZHLC_HCF(:) = ZCF(:)
        ZHLC_LCF(:) = 0.0
        ZHLC_HRC(:) = ZRCT(:)
        ZHLC_LRC(:) = 0.0
        ZRF(:)      = ZCF(:)
      ELSEWHERE (ZRCT(:) > ICED%XRTMIN(2))
        ZHLC_HCF(:) = 0.0
        ZHLC_LCF(:) = ZCF(:)
        ZHLC_HRC(:) = 0.0
        ZHLC_LRC(:) = ZRCT(:)
        ZRF(:)      = 0.
      ELSEWHERE
        ZHLC_HCF(:) = 0.0
        ZHLC_LCF(:) = 0.0
        ZHLC_HRC(:) = 0.0
        ZHLC_LRC(:) = 0.0
        ZRF(:)      = 0.
      END WHERE
    ELSEWHERE
      ZHLC_HCF(:) = 0.0
      ZHLC_LCF(:) = 0.0
      ZHLC_HRC(:) = 0.0
      ZHLC_LRC(:) = 0.0
      ZRF(:)      = 0.
    END WHERE

  ELSEIF (HSUBG_AUCV_RC == 'PDF ') THEN
    !Cloud water is split between high and low part according to a PDF
    !    'HLCRECTPDF'    : rectangular PDF form
    !    'HLCTRIANGPDF'  : triangular PDF form
    !    'HLCQUADRAPDF'  : second order quadratic PDF form
    !    'HLCISOTRIPDF'  : isocele triangular PDF
    !    'SIGM'          : Redelsperger and Sommeria (1986)

    IF ( PARAMI%CSUBG_PR_PDF == 'SIGM' ) THEN
      ! Redelsperger and Sommeria (1986) but organised according to Turner (2011, 2012)
      WHERE ( ZRCT(:) > ZRCRAUTC(:) + ZSIGMA_RC(:))
        ZHLC_HCF(:) = 1.
        ZHLC_LCF(:) = 0.0
        ZHLC_HRC(:) = ZRCT(:)
        ZHLC_LRC(:) = 0.0
        ZRF(:)      = 1.
      ELSEWHERE ( ZRCT(:) >  ( ZRCRAUTC(:) - ZSIGMA_RC(:) ) .AND. &
                & ZRCT(:) <= ( ZRCRAUTC(:) + ZSIGMA_RC(:) )       )
        ZHLC_HCF(:) = (ZRCT(:)+ZSIGMA_RC(:)-ZRCRAUTC(:))/ &
                     &(2.*ZSIGMA_RC(:))
        ZHLC_LCF(:) = MAX(0., ZCF(:)-ZHLC_HCF(:))
        ZHLC_HRC(:) = (ZRCT(:)+ZSIGMA_RC(:)-ZRCRAUTC(:))* &
                     &(ZRCT(:)+ZSIGMA_RC(:)+ZRCRAUTC(:))/ &
                     &(4.*ZSIGMA_RC(:))
        ZHLC_LRC(:) = MAX(0., ZRCT(:)-ZHLC_HRC(:))
        ZRF(:)      = ZHLC_HCF(:)
      ELSEWHERE ( ZRCT(:)>ICED%XRTMIN(2) .AND. ZCF(:)>0. )
        ZHLC_HCF(:) = 0.0
        ZHLC_LCF(:) = ZCF(:)
        ZHLC_HRC(:) = 0.0
        ZHLC_LRC(:) = ZRCT(:)
        ZRF(:)      = 0.
      ELSEWHERE
        ZHLC_HCF(:) = 0.0
        ZHLC_LCF(:) = 0.0
        ZHLC_HRC(:) = 0.0
        ZHLC_LRC(:) = 0.0
        ZRF(:)      = 0.
      END WHERE

    ! Turner (2011, 2012)
    ELSEIF ( PARAMI%CSUBG_PR_PDF== 'HLCRECTPDF' .OR. PARAMI%CSUBG_PR_PDF == 'HLCISOTRIPDF' .OR. &
           & PARAMI%CSUBG_PR_PDF == 'HLCTRIANGPDF' .OR. PARAMI%CSUBG_PR_PDF == 'HLCQUADRAPDF' ) THEN
      ! Calculate maximum value r_cM from PDF forms
      IF ( PARAMI%CSUBG_PR_PDF == 'HLCRECTPDF' .OR. PARAMI%CSUBG_PR_PDF == 'HLCISOTRIPDF' ) THEN
        ZCOEFFRCM = 2.0
      ELSE IF ( PARAMI%CSUBG_PR_PDF == 'HLCTRIANGPDF' ) THEN
        ZCOEFFRCM = 3.0
      ELSE IF ( PARAMI%CSUBG_PR_PDF == 'HLCQUADRAPDF' ) THEN
        ZCOEFFRCM = 4.0
      END IF
      WHERE (ZRCT(:).GT.0. .AND. ZCF(:).GT.0.)
        ZHLC_RCMAX(:) = ZCOEFFRCM * ZRCT(:) / ZCF(:)
      END WHERE

      ! Split available water and cloud fraction in two parts
      ! Calculate local mean values int he low and high parts for the 3 PDF forms:
      IF ( PARAMI%CSUBG_PR_PDF == 'HLCRECTPDF' ) THEN
        WHERE (ZRCT(:).GT.0. .AND. ZCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
          ZHLC_LRCLOCAL(:) = 0.5*ZRCRAUTC(:)
          ZHLC_HRCLOCAL(:) = ( ZHLC_RCMAX(:) + ZRCRAUTC(:)) / 2.0
        END WHERE
      ELSE IF ( PARAMI%CSUBG_PR_PDF == 'HLCTRIANGPDF' ) THEN
        WHERE (ZRCT(:).GT.0. .AND. ZCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
          ZHLC_LRCLOCAL(:) = ( ZRCRAUTC(:) *(3.0 * ZHLC_RCMAX(:) - 2.0 * ZRCRAUTC(:) ) ) &
                          / (3.0 * (2.0 * ZHLC_RCMAX(:) - ZRCRAUTC(:)  ) )
          ZHLC_HRCLOCAL(:) = (ZHLC_RCMAX(:) + 2.0*ZRCRAUTC(:)) / 3.0
        END WHERE
      ELSE IF ( PARAMI%CSUBG_PR_PDF == 'HLCQUADRAPDF' ) THEN
        WHERE (ZRCT(:).GT.0. .AND. ZCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
          ZHLC_LRCLOCAL(:) = (3.0 *ZRCRAUTC(:)**3 - 8.0 *ZRCRAUTC(:)**2 * ZHLC_RCMAX(:) &
                          + 6.0*ZRCRAUTC(:) *ZHLC_RCMAX(:)**2 ) &
                          / &
                          (4.0* ZRCRAUTC(:)**2 -12.0*ZRCRAUTC(:) *ZHLC_RCMAX(:) &
                          + 12.0 * ZHLC_RCMAX(:)**2 )
          ZHLC_HRCLOCAL(:) =  (ZHLC_RCMAX(:) + 3.0*ZRCRAUTC(:)) / 4.0
        END WHERE
      ELSE IF ( PARAMI%CSUBG_PR_PDF == 'HLCISOTRIPDF' ) THEN
        WHERE (ZRCT(:).GT.0. .AND. ZCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
          WHERE ( (ZRCT(:) / ZCF(:)).LE.ZRCRAUTC(:) )
            ZHLC_LRCLOCAL(:) = ( (ZHLC_RCMAX(:))**3 &
                             - (12.0 * (ZHLC_RCMAX(:))*(ZRCRAUTC(:))**2) &
                             + (8.0 * ZRCRAUTC(:)**3) ) &
                             / ( (6.0 * (ZHLC_RCMAX(:))**2) &
                             - (24.0 * (ZHLC_RCMAX(:)) * ZRCRAUTC(:)) &
                             + (12.0 * ZRCRAUTC(:)**2) )
            ZHLC_HRCLOCAL(:) = ( ZHLC_RCMAX(:) + 2.0 * ZRCRAUTC(:) ) / 3.0
          ELSEWHERE
            ZHLC_LRCLOCAL(:) = (2.0/3.0) * ZRCRAUTC(:)
            ZHLC_HRCLOCAL(:) = (3.0*ZHLC_RCMAX(:)**3 - 8.0*ZRCRAUTC(:)**3) &
                             / (6.0 * ZHLC_RCMAX(:)**2 - 12.0*ZRCRAUTC(:)**2)
          END WHERE
        END WHERE
      END IF

      ! Compare r_cM  to r_cR to know if cloud water content is high enough to split in two parts or not
      WHERE (ZRCT(:).GT.0. .AND. ZCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
        ! Calculate final values for LCF and HCF:
        ZHLC_LCF(:) = ZCF(:) &
                       * ( ZHLC_HRCLOCAL - &
                       ( ZRCT(:) / ZCF(:) ) ) &
                       / (ZHLC_HRCLOCAL - ZHLC_LRCLOCAL)
        ZHLC_HCF(:) = MAX(0., ZCF(:) - ZHLC_LCF(:))
        !
        ! Calculate final values for LRC and HRC:
        ZHLC_LRC(:) = ZHLC_LRCLOCAL * ZHLC_LCF(:)
        ZHLC_HRC(:) = MAX(0., ZRCT(:) - ZHLC_LRC(:))
      ELSEWHERE (ZRCT(:).GT.0. .AND. ZCF(:).GT.0. .AND. ZHLC_RCMAX(:).LE.ZRCRAUTC(:))
        ! Put all available cloud water and his fraction in the low part
        ZHLC_LCF(:) = ZCF(:)
        ZHLC_HCF(:) = 0.0
        ZHLC_LRC(:) = ZRCT(:)
        ZHLC_HRC(:) = 0.0
      ELSEWHERE
        ZHLC_LCF(:) = 0.
        ZHLC_HCF(:) = 0.0
        ZHLC_LRC(:) = 0.
        ZHLC_HRC(:) = 0.0
      END WHERE

      ZRF(:)=ZHLC_HCF(:) !Precipitation fraction

    ELSE
      !wrong CSUBG_PR_PDF case
      CALL ABORT
      STOP 'wrong CSUBG_PR_PDF case'
    ENDIF
  ELSE
    !wrong HSUBG_AUCV_RC case
    CALL ABORT
    STOP 'wrong HSUBG_AUCV_RC case'
  ENDIF

  !Diagnostic of precipitation fraction
  ZW(:,:) = 0.
  ZRAINFR(:,:) = UNPACK( ZRF(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )

  DO JI = D%NIB,D%NIE
    ZRAINFR(JI,IKE)=0.
    DO JK=IKE-KKL, IKB, -KKL
      IF (PRRT(JI,JK) .GT. ICED%XRTMIN(3)) THEN

        ZRAINFR(JI,JK)=MAX(ZRAINFR(JI,JK), ZRAINFR(JI,JK+KKL))

        IF (ZRAINFR(JI,JK)==0) THEN
          ZRAINFR(JI,JK)=1.
        END IF
      ELSE
        ZRAINFR(JI,JK)=0.
      END IF
    END DO
  END DO

  DO JL=1,KSIZE
    ZRF(JL)=ZRAINFR(I1(JL),I2(JL))
  END DO
!
  CALL RAIN_ICE_OLD_SLOW(D, CST, ICED, ICEP, BUCONF, &
                         KSIZE, OCND2, LMODICEDEP, &
                         PTSTEP, ZREDSN, &
                         GMICRO, PRHODJ, PTHS, PRVS, &
                         ZRCT, ZRRT, ZRIT, ZRRS, &
                         ZRGS, ZRST, ZRGT, ZCIT, &
                         ZRHODREF, ZRHODJ, ZLBDAS, &
                         ZZT, ZLSFACT, ZLVFACT, ZPRES, ZSSI, &
                         ZRVS, ZRCS, ZRIS, ZRSS, ZTHS, &
                         ZLBDAG, ZKA, ZDV, &
                         ZAI, ZCJ, ZAA2, ZBB3, &
                         ZDICRIT, ZREDGR, ZKVO, &
                         TBUDGETS, KBUDGETS)
!
!-------------------------------------------------------------------------------
!
!
!*       3.     COMPUTES THE SLOW WARM PROCESS SOURCES
!               --------------------------------------
!
!*       3.1    compute the slope parameter Lbda_r
!
  !ZLBDAR will be used when we consider rain diluted over the grid box
  WHERE( ZRRT(:)>0.0 )
    ZLBDAR(:)  = ICED%XLBR*( ZRHODREF(:)*MAX( ZRRT(:),ICED%XRTMIN(3) ) )**ICED%XLBEXR
  END WHERE
  !ZLBDAR_RF will be used when we consider rain concentrated in its fraction
  WHERE( ZRRT(:)>0.0 .AND. ZRF(:)>0.0 )
    ZLBDAR_RF(:)  = ICED%XLBR*( ZRHODREF(:) *MAX( ZRRT(:)/ZRF(:), ICED%XRTMIN(3) ) )**ICED%XLBEXR
  ELSEWHERE
    ZLBDAR_RF(:)  = 0.
  END WHERE
!
  IF( OWARM ) THEN    !  Check if the formation of the raindrops by the slow
                      !  warm processes is allowed
    PEVAP3D(:,:)= 0.
    CALL RAIN_ICE_OLD_WARM(D, CST, PARAMI, ICEP, ICED, BUCONF, &
                           KSIZE, OCND2, LKOGAN, GMICRO, &
                           PRHODJ, PEVAP3D, PTHS, PRVS, &
                           ZRVT, ZRCT, ZRRT, ZRCS, ZRRS, ZTHS, &
                           ZRVS, ZTHT, ZTHLT, &
                           ZCJ, ZKA, ZCF, ZDV, ZRF, &
                           ZACRF, ZCONCM, &
                           ZRHODREF, ZRHODJ, ZLVFACT, ZLBDAR, ZLBDAR_RF, &
                           ZZKGN_ACON, ZZKGN_SBGR, &
                           ZHLC_HCF, ZHLC_LCF, ZHLC_HRC, ZHLC_LRC, &
                           ZAA2W, ZBB3W, &
                           ZZT, ZPRES, ZESW, &
                           TBUDGETS, KBUDGETS)
  END IF
!
!-------------------------------------------------------------------------------
!
!
!*       4.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_s
!               ----------------------------------------------
!
  CALL RAIN_ICE_OLD_FAST_RS(D, CST, ICEP, ICED, BUCONF, &
                            PTSTEP, KSIZE, KRR, GMICRO, &
                            PRHODJ, PTHS, &
                            ZRVT, ZRCT, ZRRT, ZRST, &
                            ZRRS, ZRCS, ZRSS, ZRGS, ZTHS, &
                            ZRHODREF, ZRHODJ, ZLSFACT, ZLVFACT, &
                            ZCJ, ZKA, ZDV, &
                            ZLBDAR, ZLBDAS, ZCOLF, ZPRES, ZZT, &
                            TBUDGETS, KBUDGETS)
!
!-------------------------------------------------------------------------------
!
!
!*       5.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_g
!               ----------------------------------------------
!

  CALL RAIN_ICE_OLD_FAST_RG(D, CST, ICEP, ICED, BUCONF, &
                            PTSTEP, KSIZE, KRR, &
                            OCND2, LTIW, GMICRO, &
                            PRHODJ, PTHS, &
                            ZRVT, ZRCT, ZRIT, ZRRT, ZRST, ZRGT, ZCIT, &
                            ZRIS, ZRRS, ZRCS, ZRSS, ZRGS, ZRHS, ZTHS, &
                            ZRHODREF, ZRHODJ, ZLSFACT, ZLVFACT, &
                            ZCJ, ZKA, ZDV, &
                            ZLBDAR, ZLBDAG, ZLBDAS, &
                            ZTIW, ZZT, ZPRES, &
                            TBUDGETS, KBUDGETS)
!
!-------------------------------------------------------------------------------
!
!
!*       6.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_h
!               ----------------------------------------------
!
 IF ( KRR == 7 ) THEN
  CALL RAIN_ICE_OLD_FAST_RH(D, CST, ICEP, ICED, BUCONF, &
                            KSIZE, KRR, &
                            GMICRO, &
                            PTHS, PRHODJ, &
                            ZRVT, ZRCT, ZRIT, ZRST, ZRGT, ZRHT, &
                            ZRIS, ZRRS, ZRCS, ZRSS, ZRGS, ZRHS, ZTHS, &
                            ZRHODREF, ZRHODJ, ZLSFACT, ZLVFACT, &
                            ZLBDAS, ZLBDAG, ZLBDAH, &
                            ZCJ, ZKA, ZDV, &
                            ZZT, ZPRES, &
                            TBUDGETS, KBUDGETS)
 END IF
!
!-------------------------------------------------------------------------------
!
!
!*       7.     COMPUTES SPECIFIC SOURCES OF THE WARM AND COLD CLOUDY SPECIES
!               -------------------------------------------------------------
!

  CALL RAIN_ICE_OLD_FAST_RI(D, CST, ICEP, ICED, BUCONF, &
                            PTSTEP, KSIZE, &
                            OCND2, LMODICEDEP, GMICRO, &
                            PRHODJ, PTHS, &
                            ZRIT, ZCIT, &
                            ZRVS, ZRCS, ZRIS, ZRSS, ZTHS, &
                            ZRHODREF, ZRHODJ, &
                            ZLSFACT, ZLVFACT, &
                            ZAI, ZCJ, &
                            ZSSIO, ZSSIU, ZW2D, ZXW2D13, &
                            ZZT, ZPRES, ZSSI, &
                            ZSIFRC, ZESI, &
                            ZCITRED, ZCITRED23, ZDICRIT, &
                            TBUDGETS, KBUDGETS)

  IF (OCND2.AND.LCHECKNOISE) THEN
!*       8     This check is mainly for noise reduction:
!              ----------------------------------------
! Do not override saturation point over ice for temperatures below freezing.
! If so, adjust total ice and then moisture  and temperature.

    DO JL=1,KSIZE
      ZRSA=ZRIS(JL)+ZRSS(JL) +ZRGS(JL) ! total solid
      ZRSTS=ZRIT(JL)+ZRST(JL) +ZRGT(JL) ! total solid timestep t
      IF(ZZT(JL)<CST%XTT .AND. ABS(ZRSA*PTSTEP-ZRSTS)> 1.0E-12 .AND. &
      &  ZESI(JL) < ZPRES(JL)*0.5 )THEN
        ZTSP = TIWMX_TAB(ZPRES(JL),ZZT(JL), ZRVS(JL)*PTSTEP,1._JPRB,ZRSP,ZRSI,0.1_JPRB)
        ZRVSOLD =ZRVS(JL)
        ZRISOLD =ZRIS(JL)
        ZRSSOLD =ZRSS(JL)
        ZRGSOLD =ZRGS(JL)
        !       Fractions of total solid for cloud ice, snow and graupel
        !       (hail not concidered yet):
        ZRISFRC = 1._JPRB !(Default)
        ZRSSFRC = 0._JPRB !(Default)
        ZRGSFRC = 0._JPRB !(Default)
        IF(ZRSA > 0._JPRB )THEN
          ZRISFRC = ZRISOLD/ZRSA
          ZRSSFRC = ZRSSOLD/ZRSA
          ZRGSFRC = ZRGSOLD/ZRSA
        ENDIF

        ZRSDIF =0._JPRB
        ZRFRAC=   ZRVS(JL)*PTSTEP - ZRSA*PTSTEP  +ZRSTS
        IF(ZRVS(JL)*PTSTEP < ZRSI )THEN ! sub - saturation over ice:
          ! Avoid drying of ice leading to supersaturation with
          ! respect to ice

          ! ZRFRAC should not exceed ZRSP, if so adjust
          ZRSDIF = MIN(0._JPRB,ZRSP-ZRFRAC)
        ELSE  ! super - saturation over ice:
          ! ZRFRAC should not go below ZRSP, if so adjust
!          ZRSDIF = MAX(0._JPRB,ZRSP-ZRFRAC)
        ENDIF
        ZRSB = ZRSA*PTSTEP  - ZRSDIF
        ZRVS(JL) = ZRVS(JL) -  (ZRSB/PTSTEP-ZRSA) ! total H2O should not change
        ZTHS(JL) = ZTHS(JL) +  (ZRSB/PTSTEP-ZRSA)*ZLSFACT(JL) ! total energy should not change

        ZRIS(JL) = ZRSB*ZRISFRC/PTSTEP ! individual fractions should not change
        ZRSS(JL) = ZRSB*ZRSSFRC/PTSTEP ! execpt for increase from no ice, when
        ZRGS(JL) = ZRSB*ZRGSFRC/PTSTEP ! new all becomes cloud ice only. (No precipipitation)

      ENDIF
    ENDDO
    ! End check

  ENDIF

!-------------------------------------------------------------------------------

  ZW(:,:) = PRVS(:,:)
  PRVS(:,:) = UNPACK( ZRVS(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
  ZW(:,:) = PRCS(:,:)
  PRCS(:,:) = UNPACK( ZRCS(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
  ZW(:,:) = PRRS(:,:)
  PRRS(:,:) = UNPACK( ZRRS(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
  ZW(:,:) = PRIS(:,:)
  PRIS(:,:) = UNPACK( ZRIS(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
  ZW(:,:) = PRSS(:,:)
  PRSS(:,:) = UNPACK( ZRSS(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
  ZW(:,:) = PRGS(:,:)
  PRGS(:,:) = UNPACK( ZRGS(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
  IF ( KRR == 7 ) THEN
    ZW(:,:) = PRHS(:,:)
    PRHS(:,:) = UNPACK( ZRHS(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
  END IF
  ZW(:,:) = PTHS(:,:)
  PTHS(:,:) = UNPACK( ZTHS(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
  ZW(:,:) = PCIT(:,:)
  PCIT(:,:) = UNPACK( ZCIT(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
!
  ZW(:,:) = ZRAINFR(:,:)
  ZRAINFR(:,:) = UNPACK( ZRF(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
!
  ELSE
!
! Advance the budget calls
!
 IF (BUCONF%LBU_ENABLE) THEN
! Reordered for compability with flexible structures like in AROME
 ZBU0(:,:)=0.
 ! rain_ice_slow
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'HON', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'HON', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'HON', ZBU0(:,:))
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'SFR', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RR), 'SFR', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RG), 'SFR', ZBU0(:,:))
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'DEPS', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RV), 'DEPS', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RS), 'DEPS', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'AGGS', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RS), 'AGGS', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'AUTS', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RS), 'AUTS', ZBU0(:,:))
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'DEPG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RV), 'DEPG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RG), 'DEPG', ZBU0(:,:))

 IF (OWARM) THEN ! rain_ice_warm
   IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'AUTO', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RR), 'AUTO', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'ACCR', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RR), 'ACCR', ZBU0(:,:))
   IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'REVA', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RV), 'REVA', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RR), 'REVA', ZBU0(:,:))
 ENDIF

 !rain_ice_fast_rs
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'RIM', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'RIM', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RS), 'RIM', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RG), 'RIM', ZBU0(:,:))
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'ACC', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RR), 'ACC', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RS), 'ACC', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RG), 'ACC', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RS), 'CMEL', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RG), 'CMEL', ZBU0(:,:))

 !rain_ice_fast_rg
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'CFRZ', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RR), 'CFRZ', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'CFRZ', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RG), 'CFRZ', ZBU0(:,:))
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'WETG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'WETG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RR), 'WETG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'WETG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RS), 'WETG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RG), 'WETG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RH), 'WETG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'DRYG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'DRYG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RR), 'DRYG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'DRYG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RS), 'DRYG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RG), 'DRYG', ZBU0(:,:))
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'GMLT', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RR), 'GMLT', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RG), 'GMLT', ZBU0(:,:))

 IF(KRR==7) THEN ! rain_ice_fast_rh
   IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'WETH', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'WETH', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RR), 'WETH', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'WETH', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RS), 'WETH', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RG), 'WETH', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RH), 'WETH', ZBU0(:,:))
   IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'HMLT', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RR), 'HMLT', ZBU0(:,:))
   IF (BUCONF%LBUDGET_RH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RH), 'HMLT', ZBU0(:,:))
 ENDIF

 !rain_ice_fast_ri
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'IMLT', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'IMLT', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'IMLT', ZBU0(:,:))
 IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_TH), 'BERFI', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RC), 'BERFI', ZBU0(:,:))
 IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_ADD_PHY(D, TBUDGETS(NBUDGET_RI), 'BERFI', ZBU0(:,:))

 ENDIF !BUCONF%LBU_ENABLE
END IF
!
!-------------------------------------------------------------------------------
!
!*       8.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!               -------------------------------------
!
!*       8.1    time splitting loop initialization
!
IF (HSEDIM == 'STAT') THEN
  IF(BUCONF%LBU_ENABLE) THEN
    IF (BUCONF%LBUDGET_RC .AND. OSEDIC) &
                    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC),'SEDI', PRCS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR),'SEDI', PRRS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI),'SEDI', PRIS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'SEDI', PRSS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'SEDI', PRGS(:,:)*PRHODJ(:,:))
    IF (KRR == 7 .AND. BUCONF%LBUDGET_RH) &
                    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RH), 'SEDI', PRHS(:,:)*PRHODJ(:,:))
  ENDIF

  CALL RAIN_ICE_OLD_SEDIMENTATION_STAT(D, CST, ICEP, ICED, &
                                       KRR, OSEDIC, PTSTEP, KKL, IKB, IKE, &
                                       PDZZ, PRHODJ, PRHODREF, PPABST, &
                                       PTHT, PRCT, PRRT, PRST, PRGT, &
                                       PRCS, PRRS, PRIS, PRSS, PRGS, &
                                       PINPRC, PINPRR, PINPRS, PINPRG, &
                                       ZRAY, ZLBC, ZFSEDC, ZCONC3D, &
                                       PRHT, PRHS, PINPRH, PFPR)

  IF(BUCONF%LBU_ENABLE) THEN
    IF (BUCONF%LBUDGET_RC .AND. OSEDIC) &
                    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC),'SEDI', PRCS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR),'SEDI', PRRS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI),'SEDI', PRIS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'SEDI', PRSS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'SEDI', PRGS(:,:)*PRHODJ(:,:))
    IF (KRR == 7 .AND. BUCONF%LBUDGET_RH) &
                    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RH), 'SEDI', PRHS(:,:)*PRHODJ(:,:))
  ENDIF

ELSEIF (HSEDIM == 'SPLI') THEN
  IF(BUCONF%LBU_ENABLE) THEN
    IF (BUCONF%LBUDGET_RC .AND. OSEDIC) &
                    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC),'SEDI', PRCS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR),'SEDI', PRRS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI),'SEDI', PRIS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'SEDI', PRSS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'SEDI', PRGS(:,:)*PRHODJ(:,:))
    IF (KRR == 7 .AND. BUCONF%LBUDGET_RH) &
                    CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RH), 'SEDI', PRHS(:,:)*PRHODJ(:,:))
  ENDIF

  CALL RAIN_ICE_OLD_SEDIMENTATION_SPLIT(D, CST, ICEP, ICED, KSIZE, &
                                        KRR, OSEDIC, PTSTEP, KKL, IKB, KSPLITR, &
                                        PDZZ, PRHODJ, PRHODREF, PPABST, &
                                        PTHT, PRCT, PRRT, PRST, PRGT, &
                                        PRCS, PRRS, PRIS, PRSS, PRGS, &
                                        PINPRC, PINPRR, PINPRS, PINPRG, &
                                        ZRAY, ZLBC, ZFSEDC, ZCONC3D, &
                                        PRHT, PRHS, PINPRH, PFPR)

  IF(BUCONF%LBU_ENABLE) THEN
    IF (BUCONF%LBUDGET_RC .AND. OSEDIC) &
                    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC),'SEDI', PRCS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR),'SEDI', PRRS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI),'SEDI', PRIS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'SEDI', PRSS(:,:)*PRHODJ(:,:))
    IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'SEDI', PRGS(:,:)*PRHODJ(:,:))
    IF (KRR == 7 .AND. BUCONF%LBUDGET_RH) &
                    CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RH), 'SEDI', PRHS(:,:)*PRHODJ(:,:))
  ENDIF
ELSE
  WRITE(*,*) ' STOP'
  WRITE(*,*) ' NO SEDIMENTATION SCHEME FOR HSEDIM=',HSEDIM
  CALL ABORT
  STOP
END IF

  !sedimentation of rain fraction
  DO JI = D%NIB,D%NIE
    ZRAINFR(JI,IKE)=0.
    DO JK=IKE-KKL, IKB, -KKL
      IF (PRRS(JI,JK)*PTSTEP .GT. ICED%XRTMIN(3)) THEN

        ZRAINFR(JI,JK)=MAX(ZRAINFR(JI,JK), ZRAINFR(JI,JK+KKL))

        IF (ZRAINFR(JI,JK)==0) THEN
          ZRAINFR(JI,JK)=1.
        END IF
      ELSE
        ZRAINFR(JI,JK)=0.
      END IF
    END DO
  END DO

  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD',1,ZHOOK_HANDLE)

END SUBROUTINE RAIN_ICE_OLD
