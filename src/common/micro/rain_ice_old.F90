!     ######spl
      SUBROUTINE RAIN_ICE_OLD (D, CST, PARAMI, ICEP, ICED,                         &
                            OSEDIC, OCND2, LKOGAN, LMODICEDEP,                     &
                            HSEDIM, HSUBG_AUCV_RC, OWARM,                          &
                            KKA,KKU,KKL,                                           &
                            KSPLITR, PTSTEP, KRR, KSIZE, GMICRO,                   &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                            PICLDFR, PSSIO, PSSIU, PIFR,                           &
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                    &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,        &
                            PINPRC, PINPRR, PEVAP3D,                               &
                            PINPRS, PINPRG, PSIGS, PSEA, PTOWN,                    &
                            YDDDH, YDLDDH, YDMDDH,                                 &
                            PICENU, PKGN_ACON, PKGN_SBGR,                          &
                            PRHT, PRHS, PINPRH, PFPR)

      USE PARKIND1,            ONLY: JPRB
      USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
      USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_T
      USE MODD_CST,            ONLY: CST_T
      USE MODD_PARAM_ICE,      ONLY: PARAM_ICE_t
      USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_T
      USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_T
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
USE MODD_BUDGET,     ONLY: LBU_ENABLE, LBUDGET_TH, LBUDGET_RG, LBUDGET_RR, LBUDGET_RC, &
                                       LBUDGET_RI, LBUDGET_RS, LBUDGET_RH, LBUDGET_RV
USE MODD_LES,        ONLY: TLES
USE MODE_BUDGET,     ONLY: BUDGET_DDH
USE MODI_GAMMA,      ONLY: GAMMA
USE MODE_TIWMX,      ONLY: ESATI, ESATW, AA2, BB3, AA2W, BB3W
USE MODE_ICECLOUD,   ONLY: ICECLOUD
USE MODE_TIWMX_TAB,  ONLY: TIWMX_TAB
USE DDH_MIX,         ONLY: TYP_DDH
USE YOMLDDH,         ONLY: TLDDH
USE YOMMDDH,         ONLY: TMDDH
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_T
!
USE MODE_RAIN_ICE_OLD_NUCLEATION,          ONLY: RAIN_ICE_OLD_NUCLEATION
USE MODE_RAIN_ICE_OLD_SEDIMENTATION_STAT,  ONLY: RAIN_ICE_OLD_SEDIMENTATION_STAT
USE MODE_RAIN_ICE_OLD_SEDIMENTATION_SPLIT, ONLY: RAIN_ICE_OLD_SEDIMENTATION_SPLIT
USE MODE_RAIN_ICE_OLD_SLOW,                ONLY: RAIN_ICE_OLD_SLOW
!
use iso_fortran_env, only: output_unit

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_T), INTENT(IN)       :: D
TYPE(CST_T), INTENT(IN)            :: CST 
TYPE(PARAM_ICE_t),      INTENT(IN) :: PARAMI
TYPE(RAIN_ICE_PARAM_T), INTENT(IN) :: ICEP
TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED

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
TYPE(TYP_DDH),        INTENT(INOUT)     :: YDDDH
TYPE(TLDDH),          INTENT(IN)        :: YDLDDH
TYPE(TMDDH),          INTENT(IN)        :: YDMDDH
REAL, DIMENSION(D%NIT), INTENT(IN)            :: PICENU, PKGN_ACON, PKGN_SBGR
REAL, DIMENSION(D%NIT,D%NKT),   OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(D%NIT,D%NKT),   OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(D%NIT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(D%NIT,D%NKT,KRR), OPTIONAL, INTENT(OUT) :: PFPR    ! upper-air precipitation fluxes

!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK            ! Vertical loop index for the rain sedimentation
INTEGER :: JJ            ! Loop index for the interpolation
INTEGER :: JI            ! Loop index for the interpolation
INTEGER :: IKB           !
INTEGER :: IKE           !
!
INTEGER :: IMICRO ! Case number of sedimentation, T>0 (for HEN) and r_x>0 locations
INTEGER :: IGRIM, IGACC, IGDRY ! Case number of riming, accretion and dry growth
                               ! locations
INTEGER :: IGWET, IHAIL   ! wet growth locations and case number
LOGICAL, DIMENSION(D%NIT,D%NKT) :: GNEGT  ! Test where to compute the HEN process
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1,IVEC2       ! Vectors of indices for
                                ! interpolations
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for
                                ! interpolations
REAL,    DIMENSION(D%NIT,D%NKT)     :: ZW ! work array
REAL,    DIMENSION(D%NIT)           :: ZCONC_TMP    ! Weighted concentration
REAL,    DIMENSION(D%NIT,D%NKT)     :: ZT ! Temperature
REAL,    DIMENSION(D%NIT,D%NKT)     :: ZRAY,   & ! Cloud Mean radius
                                       ZLBC,   & ! XLBC weighted by sea fraction
                                       ZFSEDC, &
                                       ZCONC3D,&  !  droplet concentration m-3
                                       ZZZZ,&! geometric height
                                       ZZZT  ! tempoary value for geometric height

!Diagnostics
REAL, DIMENSION(D%NIT,D%NKT) :: ZRAINFR,    &
                                ZHLC_HCF3D, & ! HLCLOUDS cloud fraction in high water content part
                                ZHLC_LCF3D, & ! HLCLOUDS cloud fraction in low water content part
                                ZHLC_HRC3D, & ! HLCLOUDS cloud water content in high water content
                                ZHLC_LRC3D    ! HLCLOUDS cloud water content in low water content
REAL, DIMENSION(KSIZE) :: ZRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRCT    ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE) :: ZRRT    ! Rain water m.r. at t
REAL, DIMENSION(KSIZE) :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE) :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE) :: ZRGT    ! Graupel m.r. at t
REAL, DIMENSION(KSIZE) :: ZRHT    ! Hail m.r. at t
REAL, DIMENSION(KSIZE) :: ZCIT    ! Pristine ice conc. at t
!
REAL, DIMENSION(KSIZE) :: ZRVS    ! Water vapor m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRRS    ! Rain water m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRGS    ! Graupel m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRHS    ! Hail m.r. source
REAL, DIMENSION(KSIZE) :: ZTHS    ! Theta source
REAL, DIMENSION(KSIZE) :: ZTHT    ! Potential temperature
REAL, DIMENSION(KSIZE) :: ZTHLT   ! Liquid potential temperature
!
REAL, DIMENSION(KSIZE) :: ZRHODREF, & ! RHO Dry REFerence
                          ZRHODJ,   & ! RHO times Jacobian
                          ZEXNREF,  & ! EXNer Pressure REFerence
                          ZZWC,     & ! Work array
                          ZZW2,     & ! Work array
                          ZZW3,     & ! Work array
                          ZZW4,     & ! Work array
                          ZLSFACT,  & ! L_s/(Pi_ref*C_ph)
                          ZLVFACT,  & ! L_v/(Pi_ref*C_ph)
                          ZLBDAR,   & ! Slope parameter of the raindrop  distribution
                          ZLBDAR_RF,& ! Slope parameter of the raindrop  distribution
                                      ! for the Rain Fraction part
                          ZLBDAS,   & ! Slope parameter of the aggregate distribution
                          ZLBDAG,   & ! Slope parameter of the graupel   distribution
                          ZLBDAH,   & ! Slope parameter of the hail      distribution
                          ZRDRYG,   & ! Dry growth rate of the graupeln
                          ZRWETG,   & ! Wet growth rate of the graupeln
                          ZAI,      & ! Thermodynamical function
                          ZCJ,      & ! Function to compute the ventilation coefficient
                          ZKA,      & ! Thermal conductivity of the air
                          ZDV,      & ! Diffusivity of water vapor in the air
                          ZSIGMA_RC,& ! Standard deviation of rc at time t
                          ZCF,      & ! Cloud fraction
                          ZRF,      & ! Rain fraction
!                          *******  used for logical switch OCND2 : *******
                          ZSSIO,    & ! Super-saturation with respect to ice in the
                                      ! supersaturated fraction of gridbox
                          ZSSIU,    & ! Sub-saturation with respect to ice in the
                                      ! sub-saturated fraction  of gridbox
                          ZW2D,     & ! Factor for subgridscale calculations
                          ZXW2D,    & ! Ratio cloud ice moist part to dry part
                          ZXW2D13,  & ! ZXW2D**0.333 or other expression for LMODICEDEP=T
                          ZCRYSHA,  & ! Ice crystal habit factor
                          ZCI2S,    & ! factor to turn cloud ice with few lagre crystals into snow
                          ZCOLF,    & ! collision factor cloud liquid to snow / graupel
                          ZACRF,    & ! collision factor cloud liquid to rain
                          ZCONCM,   & ! Same as ZCONC3D but GMICRO-array only and cm-3 instead of m-3

                          ZHLC_HCF, & ! HLCLOUDS : fraction of High Cloud Fraction in grid
                          ZHLC_LCF, & ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                      !    note that ZCF = ZHLC_HCF + ZHLC_LCF
                          ZHLC_HRC, & ! HLCLOUDS : LWC that is High LWC in grid
                          ZHLC_LRC, & ! HLCLOUDS : LWC that is Low  LWC in grid
                                      !    note that ZRC = ZHLC_HRC + ZHLC_LRC
                        ZHLC_RCMAX, & ! HLCLOUDS : maximum value for RC in distribution
                          ZRCRAUTC, & ! RC value to begin rain formation =XCRIAUTC/RHODREF
!                          *******  end logical switch OCND2 *******
                     ZHLC_HRCLOCAL, & ! HLCLOUDS : LWC that is High LWC local in HCF
                     ZHLC_LRCLOCAL, & ! HLCLOUDS : LWC that is Low  LWC local in LCF
                                      !    note that ZRC/CF = ZHLC_HRCLOCAL+ ZHLC_LRCLOCAL
                                      !                     = ZHLC_HRC/HCF+ ZHLC_LRC/LCF
                          ZAA2,     & ! Part of ZAI used for optimized code
                          ZBB3,     & ! Part of ZAI used for optimized code
                          ZAA2W,    & ! as ZAA2 but for liquid
                          ZBB3W,    & ! as ZBB3 but for liquid
                          ZTIW,     & ! Wet bulb temperature
                          ZARTMP      ! temporary work array

REAL, DIMENSION(:), ALLOCATABLE :: &
                  ZZT,      & ! Temperature
                  ZPRES,    & ! Pressure
                  ZZW,      & ! Work array
                  ZUSW,     & ! Undersaturation over water
                  ZSSI,     & ! Supersaturation over ice
!                          *******  used for logical switch OCND2 : *******
                  ZSIFRC,   & ! subgridscale fraction with supersaturation with
                              ! respect to ice.
                  ZESI,     & ! saturation pressure over ice
                  ZESW,     & ! saturation pressure over water
!                          *******  end logical switch OCND2 *******
                  ZWLBDC      ! Slope parameter of the droplet  distribution

REAL, DIMENSION(KSIZE, KRR) :: ZZW1 ! Work arrays
REAL            :: ZTIMAUTIC,XDUMMY6,XDUMMY7
REAL            :: ZINVTSTEP

!    *******  used for logical switch OCND2 : *******
REAL            :: ZRVSOLD,ZTSP,&
     &ZRSP,ZRISOLD,ZRSSOLD,ZRGSOLD,& ! Function,old ice
     &ZRISFRC,ZRSSFRC,ZRGSFRC,ZRFRAC,ZRSA,ZRSTS,ZRSB,ZRSDIF,ZR20,  &
     &ZRSI,ZRSW,ZTC,ZHU,ZTMP,ZQIMAX,ZDICRIT,ZCITRED23,ZCITRED,ZRCW,ZVT,ZST, &
     &ZREDGR,ZREDSN, &    ! Possible reduction of the rate of graupel,snow growth
           ZRSPO, &! Mixing ratio for saturation point for
                   ! supersaturated part of gridbox
           ZRSPU, &! Mixing ratio for saturation point for
                   ! subsaturated part of gridbox
            ZKVO, &! factor used for caluclate maximum mass in the ice
                   ! distubution.
         ZTIMESC   ! Timescale for conversion lagre ice crystals to snow.
                   ! distubution.
!    *******  end logical switch OCND2 *******

! SPP arrays
REAL, DIMENSION(:), ALLOCATABLE :: ZZICENU
REAL, DIMENSION(KSIZE) :: ZZKGN_ACON,ZZKGN_SBGR

!  Tuning of sublimation (neg. sublimation)
REAL            :: ZRDEPSRED, ZRDEPGRED

     !internal fractions etc, finally saturation ratio over ice 'source' value

REAL, DIMENSION(SIZE(ICED%XRTMIN))     :: ZRTMIN
! XRTMIN = Minimum value for the mixing ratio
! ZRTMIN = Minimum value for the source (tendency)
!
INTEGER , DIMENSION(SIZE(GMICRO)) :: I1, I2 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
CHARACTER (LEN=100) :: YCOMMENT   ! Comment string in LFIFM file
CHARACTER (LEN=16)  :: YRECFM     ! Name of the desired field in LFIFM file
REAL :: ZCOEFFRCM,ZMU
LOGICAL LTEST ! Only for test !
LOGICAL LCHECKNOISE ! Noise check on/off
LOGICAL LTIW   ! Use TIW for graupel melting ( set by XFRMIN(18) ~ 1)
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.1     COMPUTE THE LOOP BOUNDS
!               -----------------------
!
IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD',0,ZHOOK_HANDLE)
LTEST=.FALSE.
LCHECKNOISE=.TRUE.
IKB=KKA+JPVEXT*KKL
IKE=KKU-JPVEXT*KKL

!
!*       1.2     COMPUTE SOME CONSTANT PARAMETERS
!
ZINVTSTEP=1./PTSTEP


ZCITRED = 0.1     ! ratio of ice crystal concentration wet to dry
                  ! part of a gridbox
ZDICRIT = ICEP%XFRMIN(15)  ! Critical diameter of ice crystal to define what
                      ! is cloud ice and what is snow (m)

ZCITRED23 = ZCITRED**0.667
IF (LMODICEDEP) THEN
  ZCITRED = 1.
  ZTIMESC = ICEP%XFRMIN(14)
  !from spherical diameter to max diameter
  ZDICRIT = (700.*CST%XPI/ICED%XAI/6.)**(1./ICED%XBI)*ZDICRIT**(3./ICED%XBI) 
  ZCITRED23 = ZCITRED**(1.+ ICED%XLBEXI)
  ZKVO      = ((ICED%XALPHAI*ICED%XNUI + ICED%XBI -1.)/ICED%XALPHAI)**(1./ICED%XALPHAI)
  ZKVO =  ZKVO/ZDICRIT/ZTIMESC
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

ZRDEPSRED = ICEP%XRDEPSRED
ZRDEPGRED = ICEP%XRDEPGRED
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

CALL RAIN_ICE_OLD_NUCLEATION(D, CST, ICEP, COUNT(ZT(D%NIB:D%NIE,D%NKTB:D%NKTE)<CST%XTT), &
                             OCND2, LMODICEDEP, KRR, PTSTEP, &
                             PTHT, PPABST, PEXNREF, PICLDFR, PRHODJ, PRHODREF, &
                             PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                             PTHS, PRVS, PRIS, PCIT, &
                             PICENU, ZT, ZZZZ, &
                             PRHT)

IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'HENU_BU_RTH',YDDDH, YDLDDH, YDMDDH)
IF (LBUDGET_RV) CALL BUDGET_DDH (PRVS(:,:)*PRHODJ(:,:),6,'HENU_BU_RRV',YDDDH, YDLDDH, YDMDDH)
IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9,'HENU_BU_RRI',YDDDH, YDLDDH, YDMDDH)

CALL COUNTJV(IMICRO, GMICRO(:,:), I1(:), I2(:))

IF ( KSIZE >= 0 ) THEN
  ALLOCATE(ZRCT(KSIZE))
  ALLOCATE(ZRCS(KSIZE))
  ALLOCATE(ZRRS(KSIZE))
  ALLOCATE(ZRIS(KSIZE))
  ALLOCATE(ZRSS(KSIZE))
  ALLOCATE(ZRGS(KSIZE))
  IF ( KRR == 7 ) ALLOCATE(ZRHS(KSIZE))
  ALLOCATE(ZZT(KSIZE))
  ALLOCATE(ZPRES(KSIZE))
  IF (OCND2) THEN
     ALLOCATE(ZESI(KSIZE))
     ALLOCATE(ZESW(KSIZE))
     ALLOCATE(ZSIFRC(KSIZE))
  ENDIF

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

  ALLOCATE(ZZW(KSIZE))

    ZZW(:)  = ZEXNREF(:)*(CST%XCPD+CST%XCPV*ZRVT(:) + CST%XCL*(ZRCT(:)+ZRRT(:)) &
                                    + CST%XCI*(ZRIT(:)+ZRST(:) + ZRGT(:)) )
    ZLSFACT(:) = (CST%XLSTT + (CST%XCPV - CST%XCI)*(ZZT(:) - CST%XTT))/ZZW(:) ! L_s/(Pi_ref*C_ph)
    ZLVFACT(:) = (CST%XLVTT + (CST%XCPV - CST%XCL)*(ZZT(:) - CST%XTT))/ZZW(:) ! L_v/(Pi_ref*C_ph)

  ALLOCATE(ZUSW(KSIZE))
  ALLOCATE(ZSSI(KSIZE))

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
  CALL RAINFR_VERT(ZRAINFR(:,:), PRRT(:,:))
  DO JL=1,KSIZE
    ZRF(JL)=ZRAINFR(I1(JL),I2(JL))
  END DO
!
  CALL RAIN_ICE_OLD_SLOW(D, CST, ICED, ICEP, &
                         KSIZE, OCND2, LMODICEDEP, &
                         PTSTEP, ZREDSN, &
                         GMICRO, PRHODJ, PTHS, PRVS, &
                         ZRCT, ZRRT, ZRIT, ZRRS, &
                         ZRGS, ZRST, ZRGT, ZCIT, &
                         ZRHODREF, ZRHODJ, ZZW2, ZLBDAS, &
                         ZZT, ZLSFACT, ZLVFACT, ZPRES, ZSSI, &
                         ZRVS, ZRCS, ZRIS, ZRSS, ZTHS, &
                         ZLBDAG, ZKA, ZDV, &
                         ZAI, ZCJ, ZAA2, ZBB3, &
                         ZRDEPSRED, ZRDEPGRED, &
                         ZDICRIT, ZREDGR, ZKVO, &
                         YDDDH, YDLDDH, YDMDDH)
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
    CALL RAIN_ICE_WARM
  END IF
!
!-------------------------------------------------------------------------------
!
!
!*       4.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_s
!               ----------------------------------------------
!
  CALL RAIN_ICE_FAST_RS
!
!-------------------------------------------------------------------------------
!
!
!*       5.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_g
!               ----------------------------------------------
!

  CALL RAIN_ICE_FAST_RG
!
!-------------------------------------------------------------------------------
!
!
!*       6.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_h
!               ----------------------------------------------
!
 IF ( KRR == 7 ) THEN
  CALL RAIN_ICE_FAST_RH
 END IF
!
!-------------------------------------------------------------------------------
!
!
!*       7.     COMPUTES SPECIFIC SOURCES OF THE WARM AND COLD CLOUDY SPECIES
!               -------------------------------------------------------------
!

  CALL RAIN_ICE_FAST_RI

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
!
!
!-------------------------------------------------------------------------------
!
!
!
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
  ZW(:,:) = 0.
  ZHLC_HCF3D(:,:) = UNPACK( ZHLC_HCF(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
!
  ZW(:,:) = 0.
  ZHLC_LCF3D(:,:) = UNPACK( ZHLC_LCF(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
!
  ZW(:,:) = 0.
  ZHLC_HRC3D(:,:) = UNPACK( ZHLC_HRC(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
!
  ZW(:,:) = 0.
  ZHLC_LRC3D(:,:) = UNPACK( ZHLC_LRC(:),MASK=GMICRO(:,:),FIELD=ZW(:,:) )
!
!
  DEALLOCATE(ZSSI)
  DEALLOCATE(ZUSW)
  DEALLOCATE(ZZW)
  DEALLOCATE(ZPRES)
  DEALLOCATE(ZZT)
  IF ( KRR == 7 ) DEALLOCATE(ZRHS)
  DEALLOCATE(ZRGS)
  DEALLOCATE(ZRSS)
  DEALLOCATE(ZRIS)
  DEALLOCATE(ZRRS)
  DEALLOCATE(ZRCS)
  DEALLOCATE(ZRCT)
  IF (OCND2) THEN
     DEALLOCATE(ZESI)
     DEALLOCATE(ZESW)
     DEALLOCATE(ZSIFRC)
  ENDIF

!
  ELSE
!
! Advance the budget calls
!
! Reordered for compability with flexible structures like in AROME

 ! rain_ice_slow
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'HON_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:)*PRHODJ(:,:),7,'HON_BU_RRC',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9,'HON_BU_RRI',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'SFR_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8,'SFR_BU_RRR',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'SFR_BU_RRG',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'DEPS_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RV) CALL BUDGET_DDH (PRVS(:,:)*PRHODJ(:,:),6,'DEPS_BU_RRV',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RS) CALL BUDGET_DDH (PRSS(:,:)*PRHODJ(:,:),10,'DEPS_BU_RRS',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9,'AGGS_BU_RRI',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RS) CALL BUDGET_DDH (PRSS(:,:)*PRHODJ(:,:),10,'AGGS_BU_RRS',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9,'AUTS_BU_RRI',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RS) CALL BUDGET_DDH (PRSS(:,:)*PRHODJ(:,:),10,'AUTS_BU_RRS',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'DEPG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RV) CALL BUDGET_DDH (PRVS(:,:)*PRHODJ(:,:),6,'DEPG_BU_RRV',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'DEPG_BU_RRG',YDDDH, YDLDDH, YDMDDH)

 IF (OWARM) THEN ! rain_ice_warm
   IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:)*PRHODJ(:,:),7,'AUTO_BU_RRC',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8,'AUTO_BU_RRR',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:)*PRHODJ(:,:),7,'ACCR_BU_RRC',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8,'ACCR_BU_RRR',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'REVA_BU_RTH',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RV) CALL BUDGET_DDH (PRVS(:,:)*PRHODJ(:,:),6,'REVA_BU_RRV',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8,'REVA_BU_RRR',YDDDH, YDLDDH, YDMDDH)
 ENDIF

 !rain_ice_fast_rs
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'RIM_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:)*PRHODJ(:,:),7,'RIM_BU_RRC',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RS) CALL BUDGET_DDH (PRSS(:,:)*PRHODJ(:,:),10,'RIM_BU_RRS',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'RIM_BU_RRG',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'ACC_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8,'ACC_BU_RRR',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RS) CALL BUDGET_DDH (PRSS(:,:)*PRHODJ(:,:),10,'ACC_BU_RRS',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'ACC_BU_RRG',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RS) CALL BUDGET_DDH (PRSS(:,:)*PRHODJ(:,:),10,'CMEL_BU_RRS',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'CMEL_BU_RRG',YDDDH, YDLDDH, YDMDDH)

 !rain_ice_fast_rg
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'CFRZ_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8,'CFRZ_BU_RRR',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9,'CFRZ_BU_RRI',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'CFRZ_BU_RRG',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'WETG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:)*PRHODJ(:,:),7,'WETG_BU_RRC',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8,'WETG_BU_RRR',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9,'WETG_BU_RRI',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RS) CALL BUDGET_DDH (PRSS(:,:)*PRHODJ(:,:),10,'WETG_BU_RRS',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'WETG_BU_RRG',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RH) CALL BUDGET_DDH (PRHS(:,:)*PRHODJ(:,:),12,'WETG_BU_RRH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'DRYG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:)*PRHODJ(:,:),7,'DRYG_BU_RRC',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8,'DRYG_BU_RRR',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9,'DRYG_BU_RRI',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RS) CALL BUDGET_DDH (PRSS(:,:)*PRHODJ(:,:),10,'DRYG_BU_RRS',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'DRYG_BU_RRG',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'GMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8,'GMLT_BU_RRR',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'GMLT_BU_RRG',YDDDH, YDLDDH, YDMDDH)

 IF(KRR==7) THEN ! rain_ice_fast_rh
   IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'WETH_BU_RTH',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:)*PRHODJ(:,:),7,'WETH_BU_RRC',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8,'WETH_BU_RRR',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9,'WETH_BU_RRI',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RS) CALL BUDGET_DDH (PRSS(:,:)*PRHODJ(:,:),10,'WETH_BU_RRS',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'WETH_BU_RRG',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RH) CALL BUDGET_DDH (PRHS(:,:)*PRHODJ(:,:),12,'WETH_BU_RRH',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'HMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8,'HMLT_BU_RRR',YDDDH, YDLDDH, YDMDDH)
   IF (LBUDGET_RH) CALL BUDGET_DDH (PRHS(:,:)*PRHODJ(:,:),12,'HMLT_BU_RRH',YDDDH, YDLDDH, YDMDDH)
 ENDIF

 !rain_ice_fast_ri
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'IMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:)*PRHODJ(:,:),7,'IMLT_BU_RRC',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9,'IMLT_BU_RRI',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:)*PRHODJ(:,:),4,'BERFI_BU_RTH',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RC) CALL BUDGET_DDH (PRCS(:,:)*PRHODJ(:,:),7,'BERFI_BU_RRC',YDDDH, YDLDDH, YDMDDH)
 IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9,'BERFI_BU_RRI',YDDDH, YDLDDH, YDMDDH)

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
  CALL RAIN_ICE_OLD_SEDIMENTATION_STAT(D, CST, ICEP, ICED, &
                                       KRR, OSEDIC, PTSTEP, KKL, IKB, IKE, &
                                       PDZZ, PRHODJ, PRHODREF, PPABST, &
                                       PTHT, PRCT, PRRT, PRST, PRGT, &
                                       PRCS, PRRS, PRIS, PRSS, PRGS, &
                                       PINPRC, PINPRR, PINPRS, PINPRG, &
                                       ZRAY, ZLBC, ZFSEDC, ZCONC3D, &
                                       PRHT, PRHS, PINPRH, PFPR)

    IF (LBUDGET_RC .AND. OSEDIC) &
                    CALL BUDGET_DDH (PRCS(:,:)*PRHODJ(:,:),7 ,'SEDI_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8 ,'SEDI_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9 ,'SEDI_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RS) CALL BUDGET_DDH (PRSS(:,:)*PRHODJ(:,:),10,'SEDI_BU_RRS',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'SEDI_BU_RRG',YDDDH, YDLDDH, YDMDDH)
    IF ( KRR == 7 .AND. LBUDGET_RH) &
                    CALL BUDGET_DDH (PRHS(:,:)*PRHODJ(:,:),12,'SEDI_BU_RRH',YDDDH, YDLDDH, YDMDDH)

ELSEIF (HSEDIM == 'SPLI') THEN
  CALL RAIN_ICE_OLD_SEDIMENTATION_SPLIT(D, CST, ICEP, ICED, KSIZE, &
                                        KRR, OSEDIC, PTSTEP, KKL, IKB, KSPLITR, &
                                        PDZZ, PRHODJ, PRHODREF, PPABST, &
                                        PTHT, PRCT, PRRT, PRST, PRGT, &
                                        PRCS, PRRS, PRIS, PRSS, PRGS, &
                                        PINPRC, PINPRR, PINPRS, PINPRG, &
                                        ZRAY, ZLBC, ZFSEDC, ZCONC3D, &
                                        PRHT, PRHS, PINPRH, PFPR)

  IF (LBUDGET_RC .AND. OSEDIC) &
                  CALL BUDGET_DDH (PRCS(:,:)*PRHODJ(:,:),7 ,'SEDI_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET_DDH (PRRS(:,:)*PRHODJ(:,:),8 ,'SEDI_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET_DDH (PRIS(:,:)*PRHODJ(:,:),9 ,'SEDI_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET_DDH (PRSS(:,:)*PRHODJ(:,:),10,'SEDI_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET_DDH (PRGS(:,:)*PRHODJ(:,:),11,'SEDI_BU_RRG',YDDDH, YDLDDH, YDMDDH)
  IF ( KRR == 7 .AND. LBUDGET_RH) &
                  CALL BUDGET_DDH (PRHS(:,:)*PRHODJ(:,:),12,'SEDI_BU_RRH',YDDDH, YDLDDH, YDMDDH)

ELSE
  WRITE(*,*) ' STOP'
  WRITE(*,*) ' NO SEDIMENTATION SCHEME FOR HSEDIM=',HSEDIM
  CALL ABORT
  STOP
END IF
!sedimentation of rain fraction
CALL RAINFR_VERT(ZRAINFR, PRRS(:,:)*PTSTEP)

!
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD',1,ZHOOK_HANDLE)
CONTAINS
!
!
!-------------------------------------------------------------------------------
!
!
      REAL FUNCTION ICENUMBER2 (Q_ICE, T3D)

      IMPLICIT NONE
      REAL, PARAMETER:: ICE_DENSITY = 890.0
      REAL, PARAMETER:: PI = 3.1415926536
      INTEGER IDX_REI
      REAL CORR, REICE, DEICE, Q_ICE, T3D
      DOUBLE PRECISION LAMBDA

!+---+-----------------------------------------------------------------+
!..Table of lookup values of radiative effective radius of ice crystals
!.. as a function of Temperature from -94C to 0C.  Taken from WRF RRTMG
!.. radiation code where it is attributed to Jon Egill Kristjansson
!.. and coauthors.
!+---+-----------------------------------------------------------------+

      REAL RETAB(95)
      DATA RETAB /                                                      &
         5.92779, 6.26422, 6.61973, 6.99539, 7.39234,                   &
         7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,          &
         10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,          &
         15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,          &
         20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,          &
         27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,          &
         31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,          &
         34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,          &
         38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,          &
         42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,          &
         50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,          &
         65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,          &
         93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,          &
         124.954, 130.630, 136.457, 142.446, 148.608, 154.956,          &
         161.503, 168.262, 175.248, 182.473, 189.952, 197.699,          &
         205.728, 214.055, 222.694, 231.661, 240.971, 250.639/

!+---+-----------------------------------------------------------------+
!..From the model 3D temperature field, subtract 179K for which
!.. index value of retab as a start.  Value of corr is for
!.. interpolating between neighboring values in the table.
!+---+-----------------------------------------------------------------+

      IDX_REI = INT(T3D-179.)
      IDX_REI = MIN(MAX(IDX_REI,1),95)
      CORR = T3D - INT(T3D)
      REICE = RETAB(IDX_REI)*(1.-CORR) + RETAB(MIN(95,IDX_REI+1))*CORR
      DEICE = 2.*REICE * 1.E-6

!+---+-----------------------------------------------------------------+
!..Now we have the final radiative effective size of ice (as function
!.. of temperature only).  This size represents 3rd moment divided by
!.. second moment of the ice size distribution, so we can compute a
!.. number concentration from the mean size and mass mixing ratio.
!.. The mean (radiative effective) diameter is 3./Slope for an inverse
!.. exponential size distribution.  So, starting with slope, work
!.. backwords to get number concentration.
!+---+-----------------------------------------------------------------+

      LAMBDA = 3.0 / DEICE
      ICENUMBER2 = Q_ICE * LAMBDA*LAMBDA*LAMBDA / (PI*ICE_DENSITY)

!+---+-----------------------------------------------------------------+
!..Example1: Common ice size coming from Thompson scheme is about 30 microns.
!.. An example ice mixing ratio could be 0.001 g/kg for a temperature of -50C.
!.. Remember to convert both into MKS units.  This gives N_ice=357652 per kg.
!..Example2: Lower in atmosphere at T=-10C matching ~162 microns in retab,
!.. and assuming we have 0.1 g/kg mixing ratio, then N_ice=28122 per kg,
!.. which is 28 crystals per liter of air if the air density is 1.0.
!+---+-----------------------------------------------------------------+

      RETURN
      END

!
!
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE RAIN_ICE_SLOW
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
!
!*       3.2     compute the homogeneous nucleation source: RCHONI
!
REAL, DIMENSION(KSIZE) :: ZBFT ! Mean time for a pristine ice crystal to reach
                               ! size of an snow/graupel particle (ZDICRIT)
REAL, DIMENSION(KSIZE) :: ZCRIAUTI ! Snow-to-ice autoconversion thres.

  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SLOW',0,ZHOOK_HANDLE)
  ZZW(:) = 0.0

  WHERE( (ZZT(:)<CST%XTT-35.0) .AND. (ZRCT(:)>ICED%XRTMIN(2)) .AND. (ZRCS(:)>0.) )
    ZZW(:) = MIN( ZRCS(:),ICEP%XHON*ZRHODREF(:)*ZRCT(:)       &
                                 *EXP( ICEP%XALPHA3*(ZZT(:)-CST%XTT)-ICEP%XBETA3 ) )
    ZRIS(:) = ZRIS(:) + ZZW(:)
    ZRCS(:) = ZRCS(:) - ZZW(:)
    ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RCHONI))
  ENDWHERE

  IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),  &
                                   4,'HON_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),   &
                                   7,'HON_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),   &
                                   9,'HON_BU_RRI',YDDDH, YDLDDH, YDMDDH)

!*       3.3     compute the spontaneous freezing source: RRHONG

  ZZW(:) = 0.0
  WHERE( (ZZT(:)<CST%XTT-35.0) .AND. (ZRRT(:)>ICED%XRTMIN(3)) .AND. (ZRRS(:)>0.) )
    ZZW(:) = MIN( ZRRS(:),ZRRT(:)* ZINVTSTEP )
    ZRGS(:) = ZRGS(:) + ZZW(:)
    ZRRS(:) = ZRRS(:) - ZZW(:)
    ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RRHONG))
  ENDWHERE

  IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),  &
                                   4,'SFR_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),   &
                                   8,'SFR_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),   &
                                  11,'SFR_BU_RRG',YDDDH, YDLDDH, YDMDDH)

!*       3.4    compute the deposition, aggregation and autoconversion sources

  ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - CST%XTT )              ! k_a
  ZDV(:) = 0.211E-4 * (ZZT(:)/CST%XTT)**1.94 * (CST%XP00/ZPRES(:)) ! D_v
!
!*       3.4.1  compute the thermodynamical function A_i(T,P)
!*              and the c^prime_j (in the ventilation factor)
!
  IF(OCND2)THEN
     ZAI(:) = ZAA2(:) + ZBB3(:)*ZPRES(:)
  ELSE
     ZAI(:) = EXP( CST%XALPI - CST%XBETAI/ZZT(:) - CST%XGAMI*ALOG(ZZT(:) ) ) ! es_i
     ZAI(:) = ( CST%XLSTT + (CST%XCPV-CST%XCI)*(ZZT(:)-CST%XTT) )**2 / (ZKA(:)*CST%XRV*ZZT(:)**2) &
                                 + ( CST%XRV*ZZT(:) ) / (ZDV(:)*ZAI(:))
  ENDIF
  ZCJ(:) = ICEP%XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-CST%XTT) )
!
!*       3.4.2  compute the riming-conversion of r_c for r_i production: RCAUTI
!
!  ZZW(:) = 0.0
!  ZTIMAUTIC = SQRT XTIMAUTI*XTIMAUTC )
!  WHERE ( (ZRCT(:)>0.0) .AND. (ZRIT(:)>0.0) .AND. (ZRCS(:)>0.0) )
!    ZZW(:) = MIN( ZRCS(:),ZTIMAUTIC * MAX( SQRT( ZRIT(:)*ZRCT(:) ),0.0 ) )
!    ZRIS(:) = ZRIS(:) + ZZW(:)
!    ZRCS(:) = ZRCS(:) - ZZW(:)
!    ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RCAUTI))
!  END WHERE
!
!*       3.4.3  compute the deposition on r_s: RVDEPS
!
  WHERE ( ZRST(:)>0.0 )
    ZLBDAS(:)  = MIN( ICED%XLBDAS_MAX,                                           &
                      ICED%XLBS*( ZRHODREF(:)*MAX( ZRST(:),ICED%XRTMIN(5) ) )**ICED%XLBEXS )
  END WHERE
  ZZW(:) = 0.0

  IF(OCND2)THEN
     WHERE ( (ZRST(:)>ICED%XRTMIN(5)) .AND. (ZRSS(:)>0.0) )
        ZZW(:) = ( ZSSI(:)/(ZRHODREF(:)*ZAI(:)) ) *  &
             ( ICEP%X0DEPS*ZLBDAS(:)**ICEP%XEX0DEPS + ICEP%X1DEPS*ZCJ(:)*ZLBDAS(:)**ICEP%XEX1DEPS )
        ZZW(:) = MIN( ZRVS(:),MAX(-ZRSS(:),ZZW(:)))  ! Simpler
        ZZW(:) = ZZW(:)*ZREDSN ! Possible tuning by using ZREDSN /=  1
        ZRSS(:) = ZRSS(:) + ZZW(:)
        ZRVS(:) = ZRVS(:) - ZZW(:)
        ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
     END WHERE
  ELSE
     WHERE ( (ZRST(:)>ICED%XRTMIN(5)) .AND. (ZRSS(:)>0.0) )
        ZZW(:) = ( ZSSI(:)/(ZRHODREF(:)*ZAI(:)) ) *          &
             ( ICEP%X0DEPS*ZLBDAS(:)**ICEP%XEX0DEPS + ICEP%X1DEPS*ZCJ(:)*ZLBDAS(:)**ICEP%XEX1DEPS )
        ZZW(:) =         MIN( ZRVS(:),ZZW(:)      )*(0.5+SIGN(0.5,ZZW(:))) &
             - MIN( ZRSS(:),ABS(ZZW(:)) )*(0.5-SIGN(0.5,ZZW(:)))
        WHERE (ZZW(:) < 0.0 )
          ZZW(:) = ZZW(:) * ZRDEPSRED
        END WHERE
        ZRSS(:) = ZRSS(:) + ZZW(:)
        ZRVS(:) = ZRVS(:) - ZZW(:)
        ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
     END WHERE
  ENDIF

  IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),  &
                                   4,'DEPS_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RV) CALL BUDGET_DDH(UNPACK(ZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:),  &
                                   6,'DEPS_BU_RRV',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),   &
                                  10,'DEPS_BU_RRS',YDDDH, YDLDDH, YDMDDH)

!*       3.4.4  compute the aggregation on r_s: RIAGGS

  ZZW(:) = 0.0
  WHERE ( (ZRIT(:)>ICED%XRTMIN(4)) .AND. (ZRST(:)>ICED%XRTMIN(5)) .AND. (ZRIS(:)>0.0) )
    ZZW(:) = MIN( ZRIS(:),ICEP%XFIAGGS * EXP( ICEP%XCOLEXIS*(ZZT(:)-CST%XTT) ) &
                                  * ZRIT(:)                      &
                                  * ZLBDAS(:)**ICEP%XEXIAGGS     &
                                  * ZRHODREF(:)**(-ICED%XCEXVT)       )
    ZRSS(:)  = ZRSS(:)  + ZZW(:)
    ZRIS(:)  = ZRIS(:)  - ZZW(:)
  END WHERE

  IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   9,'AGGS_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                  10,'AGGS_BU_RRS',YDDDH, YDLDDH, YDMDDH)

!*       3.4.5  compute the autoconversion of r_i for r_s production: RIAUTS

  ZCRIAUTI(:)=MIN(ICEP%XCRIAUTI,10**(ICEP%XACRIAUTI*(ZZT(:)-CST%XTT)+ICEP%XBCRIAUTI))
  ZZW(:) = 0.0
  WHERE ( (ZRIT(:)>ICED%XRTMIN(4)) .AND. (ZRIS(:)>0.0) )
    ZZW(:) = MIN( ZRIS(:),ICEP%XTIMAUTI * EXP( ICEP%XTEXAUTI*(ZZT(:)-CST%XTT) ) &
                            * MAX( ZRIT(:)-ZCRIAUTI(:),0.0 ) )
    ZRSS(:)  = ZRSS(:)  + ZZW(:)
    ZRIS(:)  = ZRIS(:)  - ZZW(:)
  END WHERE

  IF (OCND2 .AND. .NOT. LMODICEDEP) THEN ! 3.4.5 B:
                ! Turn ice crystals lagrer than a precribed size into snow:
                ! (For the moment sperical ice crystals are assumed)

     WHERE (  (ZRIS(:)>0.0_JPRB) .AND.(ZSSI(:)>0.001_JPRB) )
        ZBFT(:) =   0.5_JPRB*87.5_JPRB*(ZDICRIT)**2*ZAI(:)/ ZSSI(:)
        ZBFT(:) =   PTSTEP/ MAX(PTSTEP,ZBFT(:)*2._JPRB)
        ZRSS(:) =   ZRSS(:)  + ZBFT(:)*ZRIS(:)
        ZRIS(:) =   ZRIS(:)  - ZBFT(:)*ZRIS(:)
     END WHERE

  ENDIF

  IF (OCND2 .AND. LMODICEDEP) THEN ! 3.4.5 B:
                ! Turn ice to snow if ice crystal distrubution is such that
                ! the ice crystal diameter for the (mass x N_i) maximum
                ! is lagrer than a precribed size.
                ! (ZDICRIT) The general gamma function is assumed
     DO JL=1,KSIZE
        ZZW2(JL) = &
        MAX(ZCIT(JL),ICENUMBER2(ZRIS(JL)*PTSTEP,ZZT(JL))*ZRHODREF(JL))
     ENDDO

    WHERE (  ZRIS(:)>ICEP%XFRMIN(13) .AND.ZCIT(:) > 0. )
           ! LAMBDA for ICE
           ZZW2(:) = MIN(1.E8,ICED%XLBI*( ZRHODREF(:)*ZRIS(:)* PTSTEP/ZZW2(:) )**ICED%XLBEXI) 
           ZBFT(:) = 1. - 0.5**( ZKVO /ZZW2(:))
           ZBFT(:) = MIN(0.9*ZRIS(:)*PTSTEP, ZBFT(:)*ZRIS(:)*PTSTEP)
           ZRSS(:) =   ZRSS(:)  + ZBFT(:)
           ZRIS(:) =   ZRIS(:)  - ZBFT(:)
     END WHERE

  ENDIF

  IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   9,'AUTS_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                  10,'AUTS_BU_RRS',YDDDH, YDLDDH, YDMDDH)

!*       3.4.6  compute the deposition on r_g: RVDEPG

  ZZW2(:) = 0.0
  IF (ICEP%XFRMIN(5)> 1.0E-12 .AND. ICEP%XFRMIN(6) > 0.01) &
        &        ZZW2(:) = MAX(0., MIN(1., (ICEP%XFRMIN(5) - ZRGS(:))/ICEP%XFRMIN(5)))* &
        & MAX(0.,MIN(1.,ZSSI(:)/ICEP%XFRMIN(6)))


  WHERE ( ZRGT(:)>0.0 )
    ZLBDAG(:)  = ICED%XLBG*( ZRHODREF(:)*MAX( ZRGT(:),ICED%XRTMIN(6) ) )**ICED%XLBEXG
  END WHERE
  ZZW(:) = 0.0
  WHERE ( (ZRGT(:)>ICED%XRTMIN(6)) .AND. (ZRGS(:)>0.0) )
    ZZW(:) = ( ZSSI(:)/(ZRHODREF(:)*ZAI(:)) ) *                               &
             ( ICEP%X0DEPG*ZLBDAG(:)**ICEP%XEX0DEPG + ICEP%X1DEPG*ZCJ(:)*ZLBDAG(:)**ICEP%XEX1DEPG )
    ZZW(:) =         MIN( ZRVS(:),ZZW(:)      )*(0.5+SIGN(0.5,ZZW(:))) &
                   - MIN( ZRGS(:),ABS(ZZW(:)) )*(0.5-SIGN(0.5,ZZW(:)))
    ZZW(:)  = ZZW(:)*ZREDGR
    WHERE (ZZW(:) < 0.0 )
      ZZW(:)  = ZZW(:) * ZRDEPGRED
    END WHERE
    ZRSS(:) = (ZZW(:) + ZRGS(:))* ZZW2(:) + ZRSS(:)
    ZRGS(:) = (ZZW(:) + ZRGS(:))*(1. - ZZW2(:))
    ZRVS(:) = ZRVS(:) - ZZW(:)
    ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
  END WHERE
  WHERE (ZZW(:) < 0.0 )
    ZZW(:)  = ZZW(:) * ZRDEPGRED
  END WHERE

  IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),   &
                                   4,'DEPG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RV) CALL BUDGET_DDH(UNPACK(ZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:),   &
                                   6,'DEPG_BU_RRV',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                  11,'DEPG_BU_RRG',YDDDH, YDLDDH, YDMDDH)

  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_SLOW',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_SLOW
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE RAIN_ICE_WARM
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!
!-------------------------------------------------------------------------------
!
!*       4.2    compute the autoconversion of r_c for r_r production: RCAUTR
!
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_WARM',0,ZHOOK_HANDLE)

    IF (LKOGAN) THEN
       WHERE( ZRCT(:) >  1.0E-8) ! Closely following Kogan autoconversion
          ZZW(:) = 1350.0*ZZKGN_ACON(:)* ZCONCM(:)**(-1.79) * &
         &  (ZRCT(:)/(MAX(ZZKGN_SBGR(:),ZCF(:))))**2.47
          ZZW(:) = ZZW(:)*MAX(ZZKGN_SBGR(:),ZCF(:))
          ZZW(:) = MIN( ZRCS(:),ZZW(:))
          ZRCS(:) = ZRCS(:) - ZZW(:)
          ZRRS(:) = ZRRS(:) + ZZW(:)
       END WHERE
    ELSE
       WHERE( ZRCS(:)>0.0 .AND. ZHLC_HCF(:).GT.0.0 )
          ZZW(:) = ICEP%XTIMAUTC*MAX( ZHLC_HRC(:)/ZHLC_HCF(:)  - ICEP%XCRIAUTC/ZRHODREF(:),0.0)
          ZZW(:) = MIN( ZRCS(:),ZHLC_HCF(:)*ZZW(:))
          ZRCS(:) = ZRCS(:) - ZZW(:)
          ZRRS(:) = ZRRS(:) + ZZW(:)
       END WHERE
    ENDIF

      IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                       7,'AUTO_BU_RRC',YDDDH, YDLDDH, YDMDDH)
      IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                       8,'AUTO_BU_RRR',YDDDH, YDLDDH, YDMDDH)

!*       4.3    compute the accretion of r_c for r_r production: RCACCR

    IF (PARAMI%CSUBG_RC_RR_ACCR=='NONE') THEN
      !CLoud water and rain are diluted over the grid box
      WHERE( ZRCT(:)>ICED%XRTMIN(2) .AND. ZRRT(:)>ICED%XRTMIN(3) .AND. ZRCS(:)>0.0 )
        ZZW(:) = MIN( ZRCS(:), ICEP%XFCACCR * ZRCT(:)*ZACRF(:) &
                 * ZLBDAR(:)**ICEP%XEXCACCR    &
                 * ZRHODREF(:)**(-ICED%XCEXVT) )
        ZRCS(:) = ZRCS(:) - ZZW(:)
        ZRRS(:) = ZRRS(:) + ZZW(:)
      END WHERE

    ELSEIF (PARAMI%CSUBG_RC_RR_ACCR=='PRFR') THEN
      !Cloud water is concentrated over its fraction with possibly to parts with high and low content as set for autoconversion
      !Rain is concnetrated over its fraction
      !Rain in high content area fraction: ZHLC_HCF
      !Rain in low content area fraction:
      ! if ZRF<ZCF (rain is entirely falling in cloud): ZRF-ZHLC_HCF
      ! if ZRF>ZCF (rain is falling in cloud and in clear sky): ZCF-ZHLC_HCF
      ! => min(ZCF, ZRF)-ZHLC_HCF
      ZZW(:) = 0.
      WHERE( ZHLC_HRC(:)>ICED%XRTMIN(2) .AND. ZRRT(:)>ICED%XRTMIN(3) .AND. ZRCS(:)>0.0 &
            .AND. ZHLC_HCF(:)>0 )
        !Accretion due to rain falling in high cloud content
        ZZW(:) = ICEP%XFCACCR * ( ZHLC_HRC(:)/ZHLC_HCF(:) ) &
               * ZLBDAR_RF(:)**ICEP%XEXCACCR &
               * ZRHODREF(:)**(-ICED%XCEXVT) &
               * ZHLC_HCF
      END WHERE
      WHERE( ZHLC_LRC(:)>ICED%XRTMIN(2) .AND. ZRRT(:)>ICED%XRTMIN(3) .AND. ZRCS(:)>0.0 &
            .AND. ZHLC_LCF(:)>0 )
        !We add acrretion due to rain falling in low cloud content
        ZZW(:) = ZZW(:) + ICEP%XFCACCR * ( ZHLC_LRC(:)/ZHLC_LCF(:) ) &
                        * ZLBDAR_RF(:)**ICEP%XEXCACCR &
                        * ZRHODREF(:)**(-ICED%XCEXVT) &
                        * (MIN(ZCF(:), ZRF(:))-ZHLC_HCF(:))
      END WHERE
      ZZW(:)=MIN(ZRCS(:), ZZW(:))
      ZRCS(:) = ZRCS(:) - ZZW(:)
      ZRRS(:) = ZRRS(:) + ZZW(:)

    ELSE
      !wrong CSUBG_RC_RR_ACCR case
      CALL ABORT
      STOP 'wrong CSUBG_RC_RR_ACCR case'
    ENDIF

    IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     7,'ACCR_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     8,'ACCR_BU_RRR',YDDDH, YDLDDH, YDMDDH)

!*       4.4    compute the evaporation of r_r: RREVAV

    ZZW(:) = 0.0

    IF (PARAMI%CSUBG_RR_EVAP=='NONE') THEN
      !Evaporation only when there's no cloud (RC must be 0)
       IF(OCND2)THEN
          WHERE( (ZRRT(:)>ICED%XRTMIN(3)) .AND. (ZRCT(:)<=ICED%XRTMIN(2)) )
             ZZW(:) = ZAA2W(:) + ZBB3W(:)*ZPRES(:)
             ZUSW(:) = 1.0 - ZRVT(:)*( ZPRES(:)-ZESW(:) ) / ( CST%XEPSILO * ZESW(:) )
                                                    ! Undersaturation over water
             ZZW(:) = MIN( ZRRS(:),( MAX( 0.0,ZUSW(:) )/(ZRHODREF(:)*ZZW(:)) ) *      &
                  ( ICEP%X0EVAR*ZLBDAR(:)**ICEP%XEX0EVAR+ICEP%X1EVAR*ZCJ(:)*ZLBDAR(:)**ICEP%XEX1EVAR))
             ZRRS(:) = ZRRS(:) - ZZW(:)
             ZRVS(:) = ZRVS(:) + ZZW(:)
             ZTHS(:) = ZTHS(:) - ZZW(:)*ZLVFACT(:)
          END WHERE
       ELSE
          WHERE( (ZRRT(:)>ICED%XRTMIN(3)) .AND. (ZRCT(:)<=ICED%XRTMIN(2)) )
             ZZW(:)  = EXP( CST%XALPW - CST%XBETAW/ZZT(:) - CST%XGAMW*ALOG(ZZT(:) ) ) ! es_w
             ZUSW(:) = 1.0 - ZRVT(:)*( ZPRES(:)-ZZW(:) ) / ( CST%XEPSILO * ZZW(:) )
                                                        ! Undersaturation over water
             ZZW(:) = (CST%XLVTT+(CST%XCPV-CST%XCL)*(ZZT(:)-CST%XTT))**2 / (ZKA(:)*CST%XRV*ZZT(:)**2) &
               + ( CST%XRV*ZZT(:) ) / ( ZDV(:)*ZZW(:) )
             ZZW(:) = MIN( ZRRS(:),( MAX( 0.0,ZUSW(:) )/(ZRHODREF(:)*ZZW(:)) ) *      &
            ( ICEP%X0EVAR*ZLBDAR(:)**ICEP%XEX0EVAR+ICEP%X1EVAR*ZCJ(:)*ZLBDAR(:)**ICEP%XEX1EVAR))
             ZRRS(:) = ZRRS(:) - ZZW(:)
             ZRVS(:) = ZRVS(:) + ZZW(:)
            ZTHS(:) = ZTHS(:) - ZZW(:)*ZLVFACT(:)
          END WHERE
       ENDIF
    ELSEIF (PARAMI%CSUBG_RR_EVAP=='CLFR' .OR. PARAMI%CSUBG_RR_EVAP=='PRFR') THEN
      !Evaporation in clear sky part
      !With CLFR, rain is diluted over the grid box
      !With PRFR, rain is concentrated in its fraction
      !Use temperature and humidity in clear sky part like Bechtold et al. (1993)
      IF (PARAMI%CSUBG_RR_EVAP=='CLFR') THEN
        ZZW4(:)=1. !Precipitation fraction
        ZZW3(:)=ZLBDAR(:)
      ELSE
        ZZW4(:)=ZRF(:) !Precipitation fraction
        ZZW3(:)=ZLBDAR_RF(:)
      ENDIF

      !ATTENTION
      !Il faudrait recalculer les variables ZKA, ZDV, ZCJ en tenant compte de la température T^u
      !Ces variables devraient être sorties de rain_ice_slow et on mettrait le calcul de T^u, T^s
      !et plusieurs versions (comme actuellement, en ciel clair, en ciel nuageux) de ZKA, ZDV, ZCJ dans rain_ice
      !On utiliserait la bonne version suivant l'option NONE, CLFR... dans l'évaporation et ailleurs

      IF(OCND2.AND.LTEST) THEN
         DO JK= 1,KSIZE
            IF((ZRRT(JK)>ICED%XRTMIN(3)) .AND. ( ZZW4(JK) > ZCF(JK)) )THEN
              ! outside the cloud (environment) the use of T^u (unsaturated) instead of T
              ! Bechtold et al. 1993
              !
              ! T^u = T_l = theta_l * (T/theta)
               ZZW2(JK) =  ZTHLT(JK) * ZZT(JK) / ZTHT(JK) ! ZZW2 = Temperature
               ZZW(JK) = AA2W(ZZW2(JK)) + BB3W(ZZW2(JK))*ZPRES(JK) ! ZZW = Droplet function
               ZARTMP(JK)= ESATW(ZZW2(JK)) ! saturation pressure, water
            ENDIF
         ENDDO
         WHERE(  (ZRRT(:)>ICED%XRTMIN(3)) .AND. ( ZZW4(:) > ZCF(:) ) )

          ! S, Undersaturation over water (with new theta^u)
            ZUSW(:) = 1.0 - ZRVT(:)*( ZPRES(:)-ZARTMP(:) ) / ( CST%XEPSILO * ZARTMP(:) )
          ! New ZCJ(:) for  T^u
            ZARTMP(:) = ICEP%XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZW2(:)-CST%XTT) ) 
            ZZW(:) = MAX( 0.0,ZUSW(:) )/(ZRHODREF(:)*ZZW(:))  *      &
               ( ICEP%X0EVAR*ZZW3(:)**ICEP%XEX0EVAR+ICEP%X1EVAR*ZARTMP(:)*ZZW3(:)**ICEP%XEX1EVAR )
        !
            ZZW(:) = MIN( ZRRS(:),  ZZW(:) *( ZZW4(:) - ZCF(:) ) )
        !
            ZRRS(:) = ZRRS(:) - ZZW(:)
            ZRVS(:) = ZRVS(:) + ZZW(:)
            ZTHS(:) = ZTHS(:) - ZZW(:)*ZLVFACT(:)
         END WHERE
      ELSE
      WHERE(  (ZRRT(:)>ICED%XRTMIN(3)) .AND. ( ZZW4(:) > ZCF(:) ) )
        ! outside the cloud (environment) the use of T^u (unsaturated) instead of T
        ! Bechtold et al. 1993
        !
        ! T^u = T_l = theta_l * (T/theta)
        ZZW2(:) =  ZTHLT(:) * ZZT(:) / ZTHT(:)
        !
        ! es_w with new T^u
        ZZW(:)  = EXP( CST%XALPW - CST%XBETAW/ZZW2(:) - CST%XGAMW*ALOG(ZZW2(:) ) )
        !
        ! S, Undersaturation over water (with new theta^u)
        ZUSW(:) = 1.0 - ZRVT(:)*( ZPRES(:)-ZZW(:) ) / ( CST%XEPSILO * ZZW(:) )
        !
        ZZW(:) = (CST%XLVTT+(CST%XCPV-CST%XCL)*(ZZW2(:)-CST%XTT))**2 / (ZKA(:)*CST%XRV*ZZW2(:)**2) &
               + (CST%XRV*ZZW2(:)) / (ZDV(:)*ZZW(:))
        !
        ZZW(:) = MAX( 0.0,ZUSW(:) )/(ZRHODREF(:)*ZZW(:))  *      &
               ( ICEP%X0EVAR*ZZW3(:)**ICEP%XEX0EVAR+ICEP%X1EVAR*ZCJ(:)*ZZW3(:)**ICEP%XEX1EVAR )
        !
        ZZW(:) = MIN( ZRRS(:),  ZZW(:) *( ZZW4(:) - ZCF(:) ) )
        !
        ZRRS(:) = ZRRS(:) - ZZW(:)
        ZRVS(:) = ZRVS(:) + ZZW(:)
        ZTHS(:) = ZTHS(:) - ZZW(:)*ZLVFACT(:)
      END WHERE
      ENDIF
    ELSE
      !wrong CSUBG_RR_EVAP case
      CALL ABORT
      STOP 'wrong CSUBG_RR_EVAP case'
    END IF

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),   &
                                     4,'REVA_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RV) CALL BUDGET_DDH(UNPACK(ZRVS(:),MASK=GMICRO(:,:),FIELD=PRVS)*PRHODJ(:,:),   &
                                     6,'REVA_BU_RRV',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     8,'REVA_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    ZW(:,:)=PEVAP3D(:,:)
    PEVAP3D(:,:)=UNPACK(ZZW(:),MASK=GMICRO(:,:),FIELD=ZW(:,:))
!
  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_WARM',1,ZHOOK_HANDLE)
  END SUBROUTINE RAIN_ICE_WARM
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE RAIN_ICE_FAST_RS
!
!*      0. DECLARATIONS
!          ------------
!
  IMPLICIT NONE

  LOGICAL, DIMENSION(KSIZE) :: GRIM ! Test where to compute riming
  LOGICAL, DIMENSION(KSIZE) :: GACC ! Test where to compute accretion
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       5.1    cloud droplet riming of the aggregates
!
  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RS',0,ZHOOK_HANDLE)
  ZZW1(:,:) = 0.0
!
  GRIM(:) = (ZRCT(:)>ICED%XRTMIN(2)) .AND. (ZRST(:)>ICED%XRTMIN(5)) .AND.            &
                                (ZRCS(:)>0.0) .AND. (ZZT(:)<CST%XTT)
  IGRIM = COUNT( GRIM(:) )
!
  IF( IGRIM>0 ) THEN
!
!        5.1.0  allocations
!
    ALLOCATE(ZVEC1(IGRIM))
    ALLOCATE(ZVEC2(IGRIM))
    ALLOCATE(IVEC1(IGRIM))
    ALLOCATE(IVEC2(IGRIM))
!
!        5.1.1  select the ZLBDAS
!
    ZVEC1(:) = PACK( ZLBDAS(:),MASK=GRIM(:) )
!
!        5.1.2  find the next lower indice for the ZLBDAS in the geometrical
!               set of Lbda_s used to tabulate some moments of the incomplete
!               gamma function
!
    ZVEC2(1:IGRIM) = MAX( 1.00001, MIN( FLOAT(ICEP%NGAMINC)-0.00001,           &
                          ICEP%XRIMINTP1 * LOG( ZVEC1(1:IGRIM) ) + ICEP%XRIMINTP2))
    IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
    ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - FLOAT( IVEC2(1:IGRIM) )
!
!        5.1.3  perform the linear interpolation of the normalized
!               "2+XDS"-moment of the incomplete gamma function
!
    ZVEC1(1:IGRIM) =   ICEP%XGAMINC_RIM1( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                     - ICEP%XGAMINC_RIM1( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 )
!
!        5.1.4  riming of the small sized aggregates
!
    WHERE ( GRIM(:) )
      ZZW1(:,1) = MIN( ZRCS(:),                                &
                     ICEP%XCRIMSS * ZZW(:) * ZRCT(:)*ZCOLF(:)       & ! RCRIMSS
                                      *   ZLBDAS(:)**ICEP%XEXCRIMSS &
                                      * ZRHODREF(:)**(-ICED%XCEXVT) )
      ZRCS(:) = ZRCS(:) - ZZW1(:,1)
      ZRSS(:) = ZRSS(:) + ZZW1(:,1)
      ZTHS(:) = ZTHS(:) + ZZW1(:,1)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RCRIMSS))
    END WHERE
!
!        5.1.5  perform the linear interpolation of the normalized
!               "XBS"-moment of the incomplete gamma function
!
    ZVEC1(1:IGRIM) =  ICEP%XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - ICEP%XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 )
!
!        5.1.6  riming-conversion of the large sized aggregates into graupeln
!
!
    WHERE ( GRIM(:) .AND. (ZRSS(:)>0.0) )
      ZZW1(:,2) = MIN( ZRCS(:),                     &
                   ICEP%XCRIMSG * ZRCT(:)*ZCOLF(:)       & ! RCRIMSG
                           *  ZLBDAS(:)**ICEP%XEXCRIMSG  &
                           * ZRHODREF(:)**(-ICED%XCEXVT) &
                           - ZZW1(:,1)              )
      ZZW1(:,3) = MIN( ZRSS(:),                         &
                       ICEP%XSRIMCG * ZLBDAS(:)**ICEP%XEXSRIMCG   & ! RSRIMCG
                               * (1.0 - ZZW(:) )/(PTSTEP*ZRHODREF(:)) )
      ZRCS(:) = ZRCS(:) - ZZW1(:,2)
      ZRSS(:) = ZRSS(:) - ZZW1(:,3)
      ZRGS(:) = ZRGS(:) + ZZW1(:,2)+ZZW1(:,3)
      ZTHS(:) = ZTHS(:) + ZZW1(:,2)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RCRIMSG))
    END WHERE
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
  END IF

  IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),   &
                                   4,'RIM_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   7,'RIM_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                  10,'RIM_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                  11,'RIM_BU_RRG',YDDDH, YDLDDH, YDMDDH)

!*       5.2    rain accretion onto the aggregates

  ZZW1(:,2:3) = 0.0
  GACC(:) = (ZRRT(:)>ICED%XRTMIN(3)) .AND. (ZRST(:)>ICED%XRTMIN(5)) .AND.            &
                            (ZRRS(:)>0.0) .AND. (ZZT(:)<CST%XTT)
  IGACC = COUNT( GACC(:) )
!
  IF( IGACC>0 ) THEN
!
!        5.2.0  allocations
!
    ALLOCATE(ZVEC1(IGACC))
    ALLOCATE(ZVEC2(IGACC))
    ALLOCATE(ZVEC3(IGACC))
    ALLOCATE(IVEC1(IGACC))
    ALLOCATE(IVEC2(IGACC))
!
!        5.2.1  select the (ZLBDAS,ZLBDAR) couplet
!
    ZVEC1(:) = PACK( ZLBDAS(:),MASK=GACC(:) )
    ZVEC2(:) = PACK( ZLBDAR(:),MASK=GACC(:) )
!
!        5.2.2  find the next lower indice for the ZLBDAS and for the ZLBDAR
!               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
!               tabulate the RACCSS-kernel
!
    ZVEC1(1:IGACC) = MAX( 1.00001, MIN( FLOAT(ICEP%NACCLBDAS)-0.00001,           &
                          ICEP%XACCINTP1S * LOG( ZVEC1(1:IGACC) ) + ICEP%XACCINTP2S ) )
    IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
    ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - FLOAT( IVEC1(1:IGACC) )
!
    ZVEC2(1:IGACC) = MAX( 1.00001, MIN( FLOAT(ICEP%NACCLBDAR)-0.00001,           &
                          ICEP%XACCINTP1R * LOG( ZVEC2(1:IGACC) ) + ICEP%XACCINTP2R ) )
    IVEC2(1:IGACC) = INT( ZVEC2(1:IGACC) )
    ZVEC2(1:IGACC) = ZVEC2(1:IGACC) - FLOAT( IVEC2(1:IGACC) )
!
!        5.2.3  perform the bilinear interpolation of the normalized
!               RACCSS-kernel
!
    DO JJ = 1,IGACC
      ZVEC3(JJ) =  (  ICEP%XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - ICEP%XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * ZVEC1(JJ) &
                 - (  ICEP%XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - ICEP%XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
!        5.2.4  raindrop accretion on the small sized aggregates
!
    WHERE ( GACC(:) )
      ZZW1(:,2) =                                            & !! coef of RRACCS
              ICEP%XFRACCSS*( ZLBDAS(:)**ICED%XCXS )*( ZRHODREF(:)**(-ICED%XCEXVT-1.) ) &
         *( ICEP%XLBRACCS1/((ZLBDAS(:)**2)               ) +                  &
            ICEP%XLBRACCS2/( ZLBDAS(:)    * ZLBDAR(:)    ) +                  &
            ICEP%XLBRACCS3/(               (ZLBDAR(:)**2)) )/ZLBDAR(:)**4
      ZZW1(:,4) = MIN( ZRRS(:),ZZW1(:,2)*ZZW(:) )           ! RRACCSS
      ZRRS(:) = ZRRS(:) - ZZW1(:,4)*ICEP%XFRMIN(7)
      ZRSS(:) = ZRSS(:) + ZZW1(:,4)*ICEP%XFRMIN(7)
      ZTHS(:) = ZTHS(:) + ZZW1(:,4)*(ZLSFACT(:)-ZLVFACT(:))*ICEP%XFRMIN(7) ! f(L_f*(RRACCSS))
    END WHERE
!
!        5.2.4b perform the bilinear interpolation of the normalized
!               RACCS-kernel
!
    DO JJ = 1,IGACC
      ZVEC3(JJ) =  (   ICEP%XKER_RACCS(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    -  ICEP%XKER_RACCS(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                                   * ZVEC2(JJ) &
                 - (   ICEP%XKER_RACCS(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    -  ICEP%XKER_RACCS(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                           * (ZVEC2(JJ) - 1.0)
    END DO
    ZZW1(:,2) = ZZW1(:,2)*UNPACK( VECTOR=ZVEC3(:),MASK=GACC(:),FIELD=0.0 )
                                                                       !! RRACCS!
!        5.2.5  perform the bilinear interpolation of the normalized
!               SACCRG-kernel
!
    DO JJ = 1,IGACC
      ZVEC3(JJ) =  (  ICEP%XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    - ICEP%XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                          * ZVEC2(JJ) &
                 - (  ICEP%XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    - ICEP%XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                          * (ZVEC2(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
!        5.2.6  raindrop accretion-conversion of the large sized aggregates
!               into graupeln
!
    WHERE ( GACC(:) .AND. (ZRSS(:)>0.0) )
      ZZW1(:,2) = MAX( MIN( ZRRS(:),ZZW1(:,2)-ZZW1(:,4) ),0.0 )       ! RRACCSG
    END WHERE
    WHERE ( GACC(:) .AND. (ZRSS(:)>0.0) .AND. ZZW1(:,2)>0.0 .AND. ZRSS(:)>ICEP%XFRMIN(1)/PTSTEP )
      ZZW1(:,3) = MIN( ZRSS(:),ICEP%XFSACCRG*ZZW(:)*                     & ! RSACCRG
            ( ZLBDAS(:)**(ICED%XCXS-ICED%XBS) )*( ZRHODREF(:)**(-ICED%XCEXVT-1.) ) &
           *( ICEP%XLBSACCR1/((ZLBDAR(:)**2)               ) +           &
              ICEP%XLBSACCR2/( ZLBDAR(:)    * ZLBDAS(:)    ) +           &
              ICEP%XLBSACCR3/(               (ZLBDAS(:)**2)) )/ZLBDAR(:) )
      ZRRS(:) = ZRRS(:) - ZZW1(:,2)
      ZRSS(:) = ZRSS(:) - ZZW1(:,3)
      ZRGS(:) = ZRGS(:) + ZZW1(:,2)+ZZW1(:,3)
      ZTHS(:) = ZTHS(:) + ZZW1(:,2)*(ZLSFACT(:)-ZLVFACT(:)) !
                               ! f(L_f*(RRACCSG))
    END WHERE
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
  END IF

  IF (LBUDGET_TH) CALL BUDGET_DDH (                                     &
               UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:), &
                                                             4,'ACC_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET_DDH (                                      &
                   UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                                             8,'ACC_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET_DDH (                                      &
                   UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                                            10,'ACC_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET_DDH (                                      &
                   UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                                            11,'ACC_BU_RRG',YDDDH, YDLDDH, YDMDDH)

!*       5.3    Conversion-Melting of the aggregates

  ZZW(:) = 0.0
  WHERE( (ZRST(:)>ICED%XRTMIN(5)) .AND. (ZRSS(:)>0.0) .AND. (ZZT(:)>CST%XTT) )
    ZZW(:) = ZRVT(:)*ZPRES(:)/(CST%XEPSILO+ZRVT(:)) ! Vapor pressure
    ZZW(:) =  ZKA(:)*(CST%XTT-ZZT(:)) +                                 &
               ( ZDV(:)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZZT(:) - CST%XTT )) &
                           *(CST%XESTT-ZZW(:))/(CST%XRV*ZZT(:))             )
!
! compute RSMLT
!
    ZZW(:)  = MIN( ZRSS(:), ICEP%XFSCVMG*MAX( 0.0,( -ZZW(:) *             &
                           ( ICEP%X0DEPS*       ZLBDAS(:)**ICEP%XEX0DEPS + &
                             ICEP%X1DEPS*ZCJ(:)*ZLBDAS(:)**ICEP%XEX1DEPS ) -   &
                                     ( ZZW1(:,1)+ZZW1(:,4) ) *       &
                              (ZRHODREF(:)*CST%XCL*(CST%XTT-ZZT(:)))) /    &
                                             ( ZRHODREF(:)*CST%XLMTT ) ) )
!
! note that RSCVMG = RSMLT*ICEP%XFSCVMG but no heat is exchanged (at the rate RSMLT)
! because the graupeln produced by this process are still icy!!!
!
    ZRSS(:) = ZRSS(:) - ZZW(:)
    ZRGS(:) = ZRGS(:) + ZZW(:)
  END WHERE
  IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                  10,'CMEL_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                  11,'CMEL_BU_RRG',YDDDH, YDLDDH, YDMDDH)

  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RS',1,ZHOOK_HANDLE)
  END SUBROUTINE RAIN_ICE_FAST_RS
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE RAIN_ICE_FAST_RG
!
!*      0. DECLARATIONS
!          ------------
!
  IMPLICIT NONE
!
  LOGICAL, DIMENSION(KSIZE) :: GDRY ! Test where to compute dry growth
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       6.1    rain contact freezing
!
  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RG',0,ZHOOK_HANDLE)
  ZZW1(:,3:4) = 0.0
  WHERE( (ZRIT(:)>ICED%XRTMIN(4) .AND. ZRIT(:)>ICEP%XFRMIN(2)) .AND. (ZRRT(:)>ICED%XRTMIN(3)) .AND. &
                             (ZRIS(:)>0.0) .AND. (ZRRS(:)>0.0) )
    ZZW1(:,3) = MIN( ZRIS(:),ICEP%XICFRR * ZRIT(:)                & ! RICFRRG
                                    * ZLBDAR(:)**ICEP%XEXICFRR    &
                                    * ZRHODREF(:)**(-ICED%XCEXVT) )
    ZZW1(:,4) = MIN( ZRRS(:),ICEP%XRCFRI * ZCIT(:)                & ! RRCFRIG
                                    * ZLBDAR(:)**ICEP%XEXRCFRI    &
                                    * ZRHODREF(:)**(-ICED%XCEXVT-1.) )
    ZRIS(:) = ZRIS(:) - ZZW1(:,3)
    ZRRS(:) = ZRRS(:) - ZZW1(:,4)
    ZRGS(:) = ZRGS(:) + ZZW1(:,3)+ZZW1(:,4)
    ZTHS(:) = ZTHS(:) + ZZW1(:,4)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*RRCFRIG)
  END WHERE
  IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:), &
                                   4,'CFRZ_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                   8,'CFRZ_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                   9,'CFRZ_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                  11,'CFRZ_BU_RRG',YDDDH, YDLDDH, YDMDDH)

!*       6.2    compute the Dry growth case

  ZZW1(:,:) = 0.0
  WHERE( (ZRGT(:)>ICED%XRTMIN(6)) .AND. ((ZRCT(:)>ICED%XRTMIN(2) .AND. ZRCS(:)>0.0)) )
    ZZW(:) = ZLBDAG(:)**(ICED%XCXG-ICED%XDG-2.0) * ZRHODREF(:)**(-ICED%XCEXVT)
    ZZW1(:,1) = MIN( ZRCS(:),ICEP%XFCDRYG * ZRCT(:) * ZZW(:) )             ! RCDRYG
  END WHERE
  WHERE( (ZRGT(:)>ICED%XRTMIN(6)) .AND. ((ZRIT(:)>ICED%XRTMIN(4) .AND. ZRIS(:)>0.0)) )
    ZZW(:) = ZLBDAG(:)**(ICED%XCXG-ICED%XDG-2.0) * ZRHODREF(:)**(-ICED%XCEXVT)
    ZZW1(:,2) = MIN( ZRIS(:),ICEP%XFIDRYG * EXP( ICEP%XCOLEXIG*(ZZT(:)-CST%XTT) ) &
                                     * ZRIT(:) * ZZW(:) )             ! RIDRYG
  END WHERE
!
!*       6.2.1  accretion of aggregates on the graupeln
!
  GDRY(:) = (ZRST(:)>ICED%XRTMIN(5)) .AND. (ZRGT(:)>ICED%XRTMIN(6)) .AND. (ZRSS(:)>0.0)
  IGDRY = COUNT( GDRY(:) )
!
  IF( IGDRY>0 ) THEN
!
!*       6.2.2  allocations
!
    ALLOCATE(ZVEC1(IGDRY))
    ALLOCATE(ZVEC2(IGDRY))
    ALLOCATE(ZVEC3(IGDRY))
    ALLOCATE(IVEC1(IGDRY))
    ALLOCATE(IVEC2(IGDRY))
!
!*       6.2.3  select the (ZLBDAG,ZLBDAS) couplet
!
    ZVEC1(:) = PACK( ZLBDAG(:),MASK=GDRY(:) )
    ZVEC2(:) = PACK( ZLBDAS(:),MASK=GDRY(:) )
!
!*       6.2.4  find the next lower indice for the ZLBDAG and for the ZLBDAS
!               in the geometrical set of (Lbda_g,Lbda_s) couplet use to
!               tabulate the SDRYG-kernel
!
    ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( FLOAT(ICEP%NDRYLBDAG)-0.00001,           &
                          ICEP%XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + ICEP%XDRYINTP2G ) )
    IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
    ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - FLOAT( IVEC1(1:IGDRY) )
!
    ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( FLOAT(ICEP%NDRYLBDAS)-0.00001,           &
                          ICEP%XDRYINTP1S * LOG( ZVEC2(1:IGDRY) ) + ICEP%XDRYINTP2S ) )
    IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
    ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - FLOAT( IVEC2(1:IGDRY) )
!
!*       6.2.5  perform the bilinear interpolation of the normalized
!               SDRYG-kernel
!
    DO JJ = 1,IGDRY
      ZVEC3(JJ) =  (  ICEP%XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - ICEP%XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * ZVEC1(JJ) &
                 - (  ICEP%XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - ICEP%XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
    IF (OCND2) THEN
      ZZW1(:,3) = 0.
    ELSE
      WHERE( GDRY(:) )
        ZZW1(:,3) = MIN( ZRSS(:),ICEP%XFSDRYG*ZZW(:)                         & ! RSDRYG
                                        * EXP( ICEP%XCOLEXSG*(ZZT(:)-CST%XTT) )  &
                      *( ZLBDAS(:)**(ICED%XCXS-ICED%XBS) )*( ZLBDAG(:)**ICED%XCXG )    &
                      *( ZRHODREF(:)**(-ICED%XCEXVT-1.) )                    &
                           *( ICEP%XLBSDRYG1/( ZLBDAG(:)**2              ) + &
                              ICEP%XLBSDRYG2/( ZLBDAG(:)   * ZLBDAS(:)   ) + &
                              ICEP%XLBSDRYG3/(               ZLBDAS(:)**2) ) )
      END WHERE
    ENDIF
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
  END IF
!
!*       6.2.6  accretion of raindrops on the graupeln
!
  GDRY(:) = (ZRRT(:)>ICED%XRTMIN(3)) .AND. (ZRGT(:)>ICED%XRTMIN(6)) .AND. (ZRRS(:)>0.0)
  IGDRY = COUNT( GDRY(:) )
!
  IF( IGDRY>0 ) THEN
!
!*       6.2.7  allocations
!
    ALLOCATE(ZVEC1(IGDRY))
    ALLOCATE(ZVEC2(IGDRY))
    ALLOCATE(ZVEC3(IGDRY))
    ALLOCATE(IVEC1(IGDRY))
    ALLOCATE(IVEC2(IGDRY))
!
!*       6.2.8  select the (ZLBDAG,ZLBDAR) couplet
!
    ZVEC1(:) = PACK( ZLBDAG(:),MASK=GDRY(:) )
    ZVEC2(:) = PACK( ZLBDAR(:),MASK=GDRY(:) )
!
!*       6.2.9  find the next lower indice for the ZLBDAG and for the ZLBDAR
!               in the geometrical set of (Lbda_g,Lbda_r) couplet use to
!               tabulate the RDRYG-kernel
!
    ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( FLOAT(ICEP%NDRYLBDAG)-0.00001,           &
                          ICEP%XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + ICEP%XDRYINTP2G ) )
    IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
    ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - FLOAT( IVEC1(1:IGDRY) )
!
    ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( FLOAT(ICEP%NDRYLBDAR)-0.00001,           &
                          ICEP%XDRYINTP1R * LOG( ZVEC2(1:IGDRY) ) + ICEP%XDRYINTP2R ) )
    IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
    ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - FLOAT( IVEC2(1:IGDRY) )
!
!*       6.2.10 perform the bilinear interpolation of the normalized
!               RDRYG-kernel
!
    DO JJ = 1,IGDRY
      ZVEC3(JJ) =  (  ICEP%XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - ICEP%XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                  * ZVEC1(JJ) &
                 - (  ICEP%XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - ICEP%XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
    WHERE( GDRY(:) )
      ZZW1(:,4) = MIN( ZRRS(:),ICEP%XFRDRYG*ZZW(:)                    & ! RRDRYG
                        *( ZLBDAR(:)**(-4) )*( ZLBDAG(:)**ICED%XCXG ) &
                               *( ZRHODREF(:)**(-ICED%XCEXVT-1.) )   &
                    *( ICEP%XLBRDRYG1/( ZLBDAG(:)**2              ) + &
                       ICEP%XLBRDRYG2/( ZLBDAG(:)   * ZLBDAR(:)   ) + &
                       ICEP%XLBRDRYG3/(               ZLBDAR(:)**2) ) )
    END WHERE
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
  END IF
!
  ZRDRYG(:) = ZZW1(:,1) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,4)
!
!*       6.3    compute the Wet growth case
!
  ZZW(:) = 0.0
  ZRWETG(:) = 0.0
  WHERE( ZRGT(:)>ICED%XRTMIN(6) )
    ZZW1(:,5) = MIN( ZRIS(:),                                    &
                ZZW1(:,2) / (ICEP%XCOLIG*EXP(ICEP%XCOLEXIG*(ZZT(:)-CST%XTT)) ) ) ! RIWETG
    ZZW1(:,6) = MIN( ZRSS(:),                                    &
                ZZW1(:,3) / (ICEP%XCOLSG*EXP(ICEP%XCOLEXSG*(ZZT(:)-CST%XTT)) ) ) ! RSWETG
!
    ZZW(:) = ZRVT(:)*ZPRES(:)/(CST%XEPSILO+ZRVT(:)) ! Vapor pressure
    ZZW(:) =   ZKA(:)*(CST%XTT-ZZT(:)) +                              &
             ( ZDV(:)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZZT(:) - CST%XTT )) &
                           *(CST%XESTT-ZZW(:))/(CST%XRV*ZZT(:))           )
!
! compute RWETG
!
    ZRWETG(:)=MAX( 0.0,                                               &
                 ( ZZW(:) * ( ICEP%X0DEPG*       ZLBDAG(:)**ICEP%XEX0DEPG + &
                              ICEP%X1DEPG*ZCJ(:)*ZLBDAG(:)**ICEP%XEX1DEPG ) + &
                 ( ZZW1(:,5)+ZZW1(:,6) ) *                            &
                 ( ZRHODREF(:)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-ZZT(:)))   ) ) / &
                            ( ZRHODREF(:)*(CST%XLMTT-CST%XCL*(CST%XTT-ZZT(:))) )   )
  END WHERE
!
!*       6.4    Select Wet or Dry case
!
   ZZW(:) = 0.0
  IF     ( KRR == 7 ) THEN
   WHERE( ZRGT(:)>ICED%XRTMIN(6) .AND. ZZT(:)<CST%XTT                            &
                                        .AND.                          & ! Wet
                              ZRDRYG(:)>=ZRWETG(:) .AND. ZRWETG(:)>0.0 ) ! case
     ZZW(:) = ZRWETG(:) - ZZW1(:,5) - ZZW1(:,6) ! RCWETG+RRWETG
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
    ZZW1(:,7) = MAX( 0.0,MIN( ZZW(:),ZRRS(:)+ZZW1(:,1) ) )
    ZUSW(:)   = ZZW1(:,7) / ZZW(:)
    ZZW1(:,5) = ZZW1(:,5)*ZUSW(:)
    ZZW1(:,6) = ZZW1(:,6)*ZUSW(:)
    ZRWETG(:) = ZZW1(:,7) + ZZW1(:,5) + ZZW1(:,6)
!
    ZRCS(:) = ZRCS(:) - ZZW1(:,1)
    ZRIS(:) = ZRIS(:) - ZZW1(:,5)
    ZRSS(:) = ZRSS(:) - ZZW1(:,6)
!
! assume a linear percent of conversion of graupel into hail
!
    ZRGS(:) = ZRGS(:) + ZRWETG(:)                     !     Wet growth
    ZZW(:)  = ZRGS(:)*ZRDRYG(:)/(ZRWETG(:)+ZRDRYG(:)) !        and
    ZRGS(:) = ZRGS(:) - ZZW(:)                        !   partial conversion
    ZRHS(:) = ZRHS(:) + ZZW(:)                        ! of the graupel into hail
!
    ZRRS(:) = MAX( 0.0,ZRRS(:) - ZZW1(:,7) + ZZW1(:,1) )
    ZTHS(:) = ZTHS(:) + ZZW1(:,7)*(ZLSFACT(:)-ZLVFACT(:))
                                                 ! f(L_f*(RCWETG+RRWETG))
   END WHERE
   ELSE IF( KRR == 6 ) THEN
     WHERE( ZRGT(:)>ICED%XRTMIN(6) .AND. ZRGT(:)>ICEP%XFRMIN(3) .AND.            &
            ZRIS(:)*PTSTEP>ICEP%XFRMIN(3) .AND. ZZT(:)<CST%XTT              &
                                        .AND.                          & ! Wet
                              ZRDRYG(:)>=ZRWETG(:) .AND. ZRWETG(:)>0.0 ) ! case
    ZZW(:)  = ZRWETG(:)
    ZRCS(:) = ZRCS(:) - ZZW1(:,1)
    ZRIS(:) = ZRIS(:) - ZZW1(:,5)
    ZRSS(:) = ZRSS(:) - ZZW1(:,6)
    ZRGS(:) = ZRGS(:) + ZZW(:)
!
    ZRRS(:) = ZRRS(:) - ZZW(:) + ZZW1(:,5) + ZZW1(:,6) + ZZW1(:,1)
    ZTHS(:) = ZTHS(:) + (ZZW(:)-ZZW1(:,5)-ZZW1(:,6))*(ZLSFACT(:)-ZLVFACT(:))
                                                 ! f(L_f*(RCWETG+RRWETG))
   END WHERE
 END IF
  IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:), &
                                   4,'WETG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   7,'WETG_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   8,'WETG_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   9,'WETG_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                  10,'WETG_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                  11,'WETG_BU_RRG',YDDDH, YDLDDH, YDMDDH)
  IF ( KRR == 7 ) THEN
    IF (LBUDGET_RH) CALL BUDGET_DDH(UNPACK(ZRHS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    12,'WETG_BU_RRH',YDDDH, YDLDDH, YDMDDH)
  END IF


  WHERE( ZRGT(:)>ICED%XRTMIN(6) .AND. ZRGT(:)>ICEP%XFRMIN(4) .AND.               &
         ZRIS(:)*PTSTEP>ICEP%XFRMIN(4) .AND. ZZT(:)<CST%XTT                 &
                                        .AND.                          &
                               ZRDRYG(:)<ZRWETG(:) .AND. ZRDRYG(:)>0.0 ) ! Dry
    ZRCS(:) = ZRCS(:) - ZZW1(:,1)
    ZRIS(:) = ZRIS(:) - ZZW1(:,2)
    ZRSS(:) = ZRSS(:) - ZZW1(:,3)
    ZRRS(:) = ZRRS(:) - ZZW1(:,4)
    ZRGS(:) = ZRGS(:) + ZRDRYG(:)
    ZTHS(:) = ZTHS(:) + (ZZW1(:,1)+ZZW1(:,4))*(ZLSFACT(:)-ZLVFACT(:)) !
                      ! f(L_f*(RCDRYG+RRDRYG))
  END WHERE
  IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),   &
                                   4,'DRYG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   7,'DRYG_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   8,'DRYG_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   9,'DRYG_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                  10,'DRYG_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                  11,'DRYG_BU_RRG',YDDDH, YDLDDH, YDMDDH)

!*       6.5    Melting of the graupeln

  ZZW(:) = 0.0
  IF (LTIW) THEN

    WHERE( (ZRGT(:)>ICED%XRTMIN(6)) .AND. (ZRGS(:)>0.0) .AND. (ZTIW(:)>CST%XTT) )
      ZZW(:) = ZRVT(:)*ZPRES(:)/(CST%XEPSILO+ZRVT(:)) ! Vapor pressure
      ZZW(:) =  ZKA(:)*(CST%XTT-ZTIW(:)) +                                 &
                 ( ZDV(:)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZTIW(:) - CST%XTT )) &
                             *(CST%XESTT-ZZW(:))/(CST%XRV*ZTIW(:))             )
!
! compute RGMLTR
!
      ZZW(:)  = ICEP%XFRMIN(8)*MIN( ZRGS(:), MAX( 0.0,( -ZZW(:) *           &
                             ( ICEP%X0DEPG*       ZLBDAG(:)**ICEP%XEX0DEPG + &
                               ICEP%X1DEPG*ZCJ(:)*ZLBDAG(:)**ICEP%XEX1DEPG ) - &
                                       ( ZZW1(:,1)+ZZW1(:,4) ) *       &
                                ( ZRHODREF(:)*CST%XCL*(CST%XTT-ZTIW(:))) ) /   &
                                               ( ZRHODREF(:)*CST%XLMTT)))


      ZRRS(:) = ZRRS(:) + ZZW(:)
      ZRGS(:) = ZRGS(:) - ZZW(:)
      ZTHS(:) = ZTHS(:) - ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(-RGMLTR))
    END WHERE
  ELSE

    WHERE( (ZRGT(:)>ICED%XRTMIN(6)) .AND. (ZRGS(:)>0.0) .AND. (ZZT(:)>CST%XTT) )
      ZZW(:) = ZRVT(:)*ZPRES(:)/(CST%XEPSILO+ZRVT(:)) ! Vapor pressure
      ZZW(:) =  ZKA(:)*(CST%XTT-ZZT(:)) +                                 &
                 ( ZDV(:)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZZT(:) - CST%XTT )) &
                             *(CST%XESTT-ZZW(:))/(CST%XRV*ZZT(:)))
!
! compute RGMLTR
!
      ZZW(:)  = ICEP%XFRMIN(8)*MIN( ZRGS(:), MAX( 0.0,( -ZZW(:) *           &
                             ( ICEP%X0DEPG*       ZLBDAG(:)**ICEP%XEX0DEPG +  &
                               ICEP%X1DEPG*ZCJ(:)*ZLBDAG(:)**ICEP%XEX1DEPG ) - &
                                       ( ZZW1(:,1)+ZZW1(:,4) ) *       &
                                ( ZRHODREF(:)*CST%XCL*(CST%XTT-ZZT(:))) ) /    &
                                               (ZRHODREF(:)*CST%XLMTT)))
      ZRRS(:) = ZRRS(:) + ZZW(:)
      ZRGS(:) = ZRGS(:) - ZZW(:)
      ZTHS(:) = ZTHS(:) - ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(-RGMLTR))
    END WHERE
  ENDIF

    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:), &
                                     4,'GMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     8,'GMLT_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    11,'GMLT_BU_RRG',YDDDH, YDLDDH, YDMDDH)
  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RG',1,ZHOOK_HANDLE)
!
  END SUBROUTINE RAIN_ICE_FAST_RG
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE RAIN_ICE_FAST_RH
!
!*      0. DECLARATIONS
!          ------------
!
  IMPLICIT NONE

  LOGICAL, DIMENSION(KSIZE) :: GWET  ! Test where to compute wet growth
  LOGICAL, DIMENSION(KSIZE) :: GHAIL ! Test where to compute hail growth

  REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RH',0,ZHOOK_HANDLE)

  GHAIL(:) = ZRHT(:)>ICED%XRTMIN(7)
  IHAIL = COUNT(GHAIL(:))
!
  IF( IHAIL>0 ) THEN
!
!*       7.2    compute the Wet growth of hail
!
    WHERE ( GHAIL(:) )
      ZLBDAH(:)  = ICED%XLBH*( ZRHODREF(:)*MAX( ZRHT(:),ICED%XRTMIN(7) ) )**ICED%XLBEXH
    END WHERE
!
    ZZW1(:,:) = 0.0
    WHERE( GHAIL(:) .AND. ((ZRCT(:)>ICED%XRTMIN(2) .AND. ZRCS(:)>0.0)) )
      ZZW(:) = ZLBDAH(:)**(ICED%XCXH-ICED%XDH-2.0) * ZRHODREF(:)**(-ICED%XCEXVT)
      ZZW1(:,1) = MIN( ZRCS(:),ICEP%XFWETH * ZRCT(:) * ZZW(:) )             ! RCWETH
    END WHERE
    WHERE( GHAIL(:) .AND. ((ZRIT(:)>ICED%XRTMIN(4) .AND. ZRIS(:)>0.0)) )
      ZZW(:) = ZLBDAH(:)**(ICED%XCXH-ICED%XDH-2.0) * ZRHODREF(:)**(-ICED%XCEXVT)
      ZZW1(:,2) = MIN( ZRIS(:),ICEP%XFWETH * ZRIT(:) * ZZW(:) )             ! RIWETH
    END WHERE
!
!*       7.2.1  accretion of aggregates on the hailstones
!
    GWET(:) = GHAIL(:) .AND. (ZRST(:)>ICED%XRTMIN(5) .AND. ZRSS(:)>0.0)
    IGWET = COUNT( GWET(:) )
!
    IF( IGWET>0 ) THEN
!
!*       7.2.2  allocations
!
      ALLOCATE(ZVEC1(IGWET))
      ALLOCATE(ZVEC2(IGWET))
      ALLOCATE(ZVEC3(IGWET))
      ALLOCATE(IVEC1(IGWET))
      ALLOCATE(IVEC2(IGWET))
!
!*       7.2.3  select the (ZLBDAH,ZLBDAS) couplet
!
      ZVEC1(:) = PACK( ZLBDAH(:),MASK=GWET(:) )
      ZVEC2(:) = PACK( ZLBDAS(:),MASK=GWET(:) )
!
!*       7.2.4  find the next lower indice for the ZLBDAG and for the ZLBDAS
!               in the geometrical set of (Lbda_h,Lbda_s) couplet use to
!               tabulate the SWETH-kernel
!
      ZVEC1(1:IGWET) = MAX( 1.00001, MIN( FLOAT(ICEP%NWETLBDAH)-0.00001,           &
                            ICEP%XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + ICEP%XWETINTP2H ) )
      IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
      ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - FLOAT( IVEC1(1:IGWET) )
!
      ZVEC2(1:IGWET) = MAX( 1.00001, MIN( FLOAT(ICEP%NWETLBDAS)-0.00001,           &
                            ICEP%XWETINTP1S * LOG( ZVEC2(1:IGWET) ) + ICEP%XWETINTP2S ) )
      IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
      ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - FLOAT( IVEC2(1:IGWET) )
!
!*       7.2.5  perform the bilinear interpolation of the normalized
!               SWETH-kernel
!
      DO JJ = 1,IGWET
        ZVEC3(JJ) = (  ICEP%XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - ICEP%XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                   * ZVEC1(JJ) &
                   - ( ICEP%XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - ICEP%XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
      END DO
      ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
      WHERE( GWET(:) )
        ZZW1(:,3) = MIN( ZRSS(:),ICEP%XFSWETH*ZZW(:)                       & ! RSWETH
                      *( ZLBDAS(:)**(ICED%XCXS-ICED%XBS) )*( ZLBDAH(:)**ICED%XCXH )  &
                         *( ZRHODREF(:)**(-ICED%XCEXVT-1.) )               &
                         *( ICEP%XLBSWETH1/( ZLBDAH(:)**2              ) + &
                            ICEP%XLBSWETH2/( ZLBDAH(:)   * ZLBDAS(:)   ) + &
                            ICEP%XLBSWETH3/(               ZLBDAS(:)**2) ) )
      END WHERE
      DEALLOCATE(IVEC2)
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC3)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(ZVEC1)
    END IF
!
!*       7.2.6  accretion of graupeln on the hailstones
!
    GWET(:) = GHAIL(:) .AND. (ZRGT(:)>ICED%XRTMIN(6) .AND. ZRGS(:)>0.0)
    IGWET = COUNT( GWET(:) )
!
    IF( IGWET>0 ) THEN
!
!*       7.2.7  allocations
!
      ALLOCATE(ZVEC1(IGWET))
      ALLOCATE(ZVEC2(IGWET))
      ALLOCATE(ZVEC3(IGWET))
      ALLOCATE(IVEC1(IGWET))
      ALLOCATE(IVEC2(IGWET))
!
!*       7.2.8  select the (ZLBDAH,ZLBDAG) couplet
!
      ZVEC1(:) = PACK( ZLBDAH(:),MASK=GWET(:) )
      ZVEC2(:) = PACK( ZLBDAG(:),MASK=GWET(:) )
!
!*       7.2.9  find the next lower indice for the ZLBDAH and for the ZLBDAG
!               in the geometrical set of (Lbda_h,Lbda_g) couplet use to
!               tabulate the GWETH-kernel
!
      ZVEC1(1:IGWET) = MAX( 1.00001, MIN( FLOAT(ICEP%NWETLBDAG)-0.00001,           &
                            ICEP%XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + ICEP%XWETINTP2H ) )
      IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
      ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - FLOAT( IVEC1(1:IGWET) )
!
      ZVEC2(1:IGWET) = MAX( 1.00001, MIN( FLOAT(ICEP%NWETLBDAG)-0.00001,           &
                            ICEP%XWETINTP1G * LOG( ZVEC2(1:IGWET) ) + ICEP%XWETINTP2G ) )
      IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
      ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - FLOAT( IVEC2(1:IGWET) )
!
!*       7.2.10 perform the bilinear interpolation of the normalized
!               GWETH-kernel
!
      DO JJ = 1,IGWET
        ZVEC3(JJ) = (  ICEP%XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - ICEP%XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                   * ZVEC1(JJ) &
                  - (  ICEP%XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - ICEP%XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
      END DO
      ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
      WHERE( GWET(:) )
        ZZW1(:,5) = MAX(MIN( ZRGS(:),ICEP%XFGWETH*ZZW(:)                       & ! RGWETH
                      *( ZLBDAG(:)**(ICED%XCXG-ICED%XBG) )*( ZLBDAH(:)**ICED%XCXH )  &
                         *( ZRHODREF(:)**(-ICED%XCEXVT-1.) )               &
                         *( ICEP%XLBGWETH1/( ZLBDAH(:)**2              ) + &
                            ICEP%XLBGWETH2/( ZLBDAH(:)   * ZLBDAG(:)   ) + &
                            ICEP%XLBGWETH3/(               ZLBDAG(:)**2) ) ),0. )
      END WHERE
      DEALLOCATE(IVEC2)
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC3)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(ZVEC1)
    END IF
!
!*       7.3    compute the Wet growth of hail
!
    ZZW(:) = 0.0
    WHERE( GHAIL(:) .AND. ZZT(:)<CST%XTT )
      ZZW(:) = ZRVT(:)*ZPRES(:)/(CST%XEPSILO+ZRVT(:)) ! Vapor pressure
      ZZW(:) = ZKA(:)*(CST%XTT-ZZT(:)) +                                 &
                ( ZDV(:)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZZT(:) - CST%XTT )) &
                            *(CST%XESTT-ZZW(:))/(CST%XRV*ZZT(:)))
!
! compute RWETH
!
      ZZW(:)  =  MAX(0.,  ( ZZW(:) * ( ICEP%X0DEPH*       ZLBDAH(:)**ICEP%XEX0DEPH + &
                                ICEP%X1DEPH*ZCJ(:)*ZLBDAH(:)**ICEP%XEX1DEPH ) + &
                   ( ZZW1(:,2)+ZZW1(:,3)+ZZW1(:,5) ) *                  &
                   ( ZRHODREF(:)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-ZZT(:))))) / &
                         ( ZRHODREF(:)*(CST%XLMTT-CST%XCL*(CST%XTT-ZZT(:))) ) )
!
      ZZW1(:,6) = MAX( ZZW(:) - ZZW1(:,2) - ZZW1(:,3) - ZZW1(:,5),0.) ! RCWETH+RRWETH
    END WHERE
    WHERE ( GHAIL(:) .AND. ZZT(:)<CST%XTT  .AND. ZZW1(:,6)/=0.)
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
      ZZW1(:,4) = MAX( 0.0,MIN( ZZW1(:,6),ZRRS(:)+ZZW1(:,1) ) )
      ZUSW(:)   = ZZW1(:,4) / ZZW1(:,6)
      ZZW1(:,2) = ZZW1(:,2)*ZUSW(:)
      ZZW1(:,3) = ZZW1(:,3)*ZUSW(:)
      ZZW1(:,5) = ZZW1(:,5)*ZUSW(:)
      ZZW(:)    = ZZW1(:,4) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,5)
!
!*       7.1.6  integrate the Wet growth of hail
!
      ZRCS(:) = ZRCS(:) - ZZW1(:,1)
      ZRIS(:) = ZRIS(:) - ZZW1(:,2)
      ZRSS(:) = ZRSS(:) - ZZW1(:,3)
      ZRGS(:) = ZRGS(:) - ZZW1(:,5)
      ZRHS(:) = ZRHS(:) + ZZW(:)
      ZRRS(:) = MAX( 0.0,ZRRS(:) - ZZW1(:,4) + ZZW1(:,1) )
      ZTHS(:) = ZTHS(:) + ZZW1(:,4)*(ZLSFACT(:)-ZLVFACT(:))
                           ! f(L_f*(RCWETH+RRWETH))
    END WHERE
  END IF
    IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:), &
                                     4,'WETH_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     7,'WETH_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     8,'WETH_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                     9,'WETH_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RS) CALL BUDGET_DDH(UNPACK(ZRSS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    10,'WETH_BU_RRS',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RG) CALL BUDGET_DDH(UNPACK(ZRGS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    11,'WETH_BU_RRG',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RH) CALL BUDGET_DDH(UNPACK(ZRHS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                    12,'WETH_BU_RRH',YDDDH, YDLDDH, YDMDDH)

!*       7.45   Conversion of the hailstones into graupel

  IF( IHAIL>0 ) THEN
!
!*       7.5    Melting of the hailstones
!
    ZZW(:) = 0.0
    WHERE( GHAIL(:) .AND. (ZRHS(:)>0.0) .AND. (ZZT(:)>CST%XTT) )
      ZZW(:) = ZRVT(:)*ZPRES(:)/(CST%XEPSILO+ZRVT(:)) ! Vapor pressure
      ZZW(:) = ZKA(:)*(CST%XTT-ZZT(:)) +                              &
             ( ZDV(:)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZZT(:) - CST%XTT )) &
                             *(CST%XESTT-ZZW(:))/(CST%XRV*ZZT(:)))
!
! compute RHMLTR
!
      ZZW(:)  = MIN( ZRHS(:), MAX( 0.0,( -ZZW(:) *                     &
                             ( ICEP%X0DEPH*       ZLBDAH(:)**ICEP%XEX0DEPH +     &
                               ICEP%X1DEPH*ZCJ(:)*ZLBDAH(:)**ICEP%XEX1DEPH ) -   &
                      ZZW1(:,6)*( ZRHODREF(:)*CST%XCL*(CST%XTT-ZZT(:))) ) /    &
                                               ( ZRHODREF(:)*CST%XLMTT)))
      ZRRS(:) = ZRRS(:) + ZZW(:)
      ZRHS(:) = ZRHS(:) - ZZW(:)
      ZTHS(:) = ZTHS(:) - ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(-RHMLTR))
    END WHERE
  END IF

  IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),&
                                   4,'HMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET_DDH(UNPACK(ZRRS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                   8,'HMLT_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RH) CALL BUDGET_DDH(UNPACK(ZRHS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0), &
                                  12,'HMLT_BU_RRH',YDDDH, YDLDDH, YDMDDH)

  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RH',1,ZHOOK_HANDLE)
  END SUBROUTINE RAIN_ICE_FAST_RH
!
!-------------------------------------------------------------------------------
!
!
  SUBROUTINE RAIN_ICE_FAST_RI
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
!*       7.1    cloud ice melting
!
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RI',0,ZHOOK_HANDLE)
  ZZW(:) = 0.0
  WHERE( (ZRIS(:)>0.0) .AND. (ZZT(:)>CST%XTT) )
    ZZW(:)  = ZRIS(:)
    ZRCS(:) = ZRCS(:) + ZRIS(:)
    ZTHS(:) = ZTHS(:) - ZRIS(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(-RIMLTC))
    ZRIS(:) = 0.0
    ZCIT(:) = 0.0
  END WHERE

  IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),   &
                                   4,'IMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   7,'IMLT_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   9,'IMLT_BU_RRI',YDDDH, YDLDDH, YDMDDH)

!*       7.2    Bergeron-Findeisen effect: RCBERI

  ZZW(:) = 0.0
  IF(OCND2)THEN

     ! Sub gridscale decomposition into a supersaturation part of the gridbox,
     ! ZSIFRC with a superaturation ZSSIO and a subsaturated part (1.- ZSIFRC)
     ! with a (negative) superaturation of ZSSIU

     IF (LMODICEDEP) THEN

       DO JL=1,KSIZE
         ZZW2(JL) = MAX(ZCIT(JL),ICENUMBER2(ZRIS(JL)*PTSTEP,ZZT(JL))* &
         ZRHODREF(JL))
       ENDDO

       WHERE( ZZW2(:)>0.0 .AND. ZESI(:) < ZPRES(:)*0.5)
          ZZW(:)= ICEP%X0DEPI/(ICED%XLBI*ZAI(:)) *(ZZW2(:)/ZRHODREF(:))**(1.+ICED%XLBEXI) * &
             & (PTSTEP*MAX(ICED%XRTMIN(4)/PTSTEP,ZRIS(:))*ZW2D(:) )**(-ICED%XLBEXI)
          ZZW(:)=  MAX(-ZRIS(:)*ZW2D(:)*(1.-ZSIFRC(:))+ZZW(:)*ZSSIO(:)* ZSIFRC(:)* ZXW2D13(:), &
        &  ZZW(:)* ( ZSSIO(:)* ZSIFRC(:)* ZXW2D13(:)  + ZCITRED23*ZSSIU(:)* (1.-ZSIFRC(:)) ))

          ZRIS(:) = ZRIS(:) + ZZW(:)
          ZRVS(:) = ZRVS(:) - ZZW(:)  ! Budget here: ! cloud ice + vapor = const
          ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:) ! f(L_f*(RCBERI))

       END WHERE

     ELSE

      DO JK=1,KSIZE

        ZTC =  MAX(-18.,MIN(-1.,ZZT(JK)-CST%XTT))
        ZHU =  MIN(0.15,MAX(0.,ZSSI(JK)))
        ZCRYSHA(JK)=1.1+ 3.*ZHU*(1.+ SIN(0.64*ZTC -1.3))
!       icedensity*4/3 *pi /8. =366.5 ; icedensity=700 kg/m3
        ZQIMAX = 366.5 * ZDICRIT**3 * ZCIT(JK)*ZCITRED/ZRHODREF(JK)
        ZCI2S(JK) = 0.
        IF(ZRIS(JK)*PTSTEP > 1.0e-12)THEN
            ZCI2S(JK)  =  ZRIS(JK)*(1. - MIN(1., 0.5*ZQIMAX /ZRIS(JK)/PTSTEP))* &
                &  (1.-ZSIFRC(JK))*ZW2D(JK)
        ENDIF

      ENDDO

      WHERE( ZCIT(:)>0.0 .AND. ZESI(:) < ZPRES(:)*0.5)
        ZZWC(:)=ZCRYSHA(:)*0.878/ZAI(:)*(ZCIT(:)/ZRHODREF(:))**0.667 &
             &*(MAX(ICED%XRTMIN(4)/PTSTEP,ZRIS(:))*PTSTEP*ZW2D(:))**0.333
!     Ice supersaturated part of grid box:
        WHERE( ZSSIO(:)>0. .AND. ZSIFRC(:) > 0.02_JPRB )
           ZZW(:)  = ZZWC(:)*ZXW2D13(:)*ZSSIO(:)
           ZRIS(:) = ZRIS(:) + ZZW(:)*ZSIFRC(:)
           ZRVS(:) = ZRVS(:) - ZZW(:)*ZSIFRC(:)  ! Budget here: ! cloud ice + vapor = const
           ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)*ZSIFRC(:) ! f(L_f*(RCBERI))
        END WHERE
!    Ice subsaturated part of grid box:
        WHERE(  ZSSIU(:)<0. .AND. ZSIFRC(:) <0.98_JPRB )
           ZRIS(:) =ZRIS(:) - ZCI2S(:)
           ZRSS(:) =ZRSS(:) + ZCI2S(:)
           ZZW(:)  =ZZWC(:)*ZCITRED23*ZSSIU(:)
           ZRIS(:) = ZRIS(:) + ZZW(:)*(1.-ZSIFRC(:))
           ZRVS(:) = ZRVS(:) - ZZW(:)*(1.-ZSIFRC(:))
           ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)*(1.-ZSIFRC(:))
        END WHERE
      END WHERE

     ENDIF

  ELSE ! End OCND2
  WHERE( (ZRCS(:)>0.0) .AND. (ZSSI(:)>0.0) .AND. &
         (ZRIT(:)>ICED%XRTMIN(4)) .AND. (ZCIT(:)>0.0)       )
    ZZW(:) = MIN(1.E8,ICED%XLBI*( ZRHODREF(:)*ZRIT(:)/ZCIT(:) )**ICED%XLBEXI) ! Lbda_i
    ZZW(:) = MIN( ZRCS(:),( ZSSI(:) / (ZRHODREF(:)*ZAI(:)) ) * ZCIT(:) * &
                  ( ICEP%X0DEPI/ZZW(:) + ICEP%X2DEPI*ZCJ(:)*ZCJ(:)/ZZW(:)**(ICED%XDI+2.0) ) )
    ZRCS(:) = ZRCS(:) - ZZW(:)
    ZRIS(:) = ZRIS(:) + ZZW(:)
    ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:))
  END WHERE
  ENDIF

  IF (LBUDGET_TH) CALL BUDGET_DDH(UNPACK(ZTHS(:),MASK=GMICRO(:,:),FIELD=PTHS)*PRHODJ(:,:),   &
                                   4,'BERFI_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET_DDH(UNPACK(ZRCS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   7,'BERFI_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET_DDH(UNPACK(ZRIS(:)*ZRHODJ(:),MASK=GMICRO(:,:),FIELD=0.0),    &
                                   9,'BERFI_BU_RRI',YDDDH, YDLDDH, YDMDDH)

  IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAIN_ICE_FAST_RI',1,ZHOOK_HANDLE)

  END SUBROUTINE RAIN_ICE_FAST_RI
!
SUBROUTINE RAINFR_VERT(ZPRFR, ZRR)

IMPLICIT NONE
REAL, DIMENSION(:,:), INTENT(OUT) :: ZPRFR !Precipitation fraction
REAL, DIMENSION(:,:), INTENT(IN)  :: ZRR !Rain field
!
!-------------------------------------------------------------------------------
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JI, JK
IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAINFR_VERT',0,ZHOOK_HANDLE)
!
  DO JI = D%NIB,D%NIE
    ZPRFR(JI,IKE)=0.
    DO JK=IKE-KKL, IKB, -KKL
       IF (ZRR(JI,JK) .GT. ICED%XRTMIN(3)) THEN
          ZPRFR(JI,JK)=MAX(ZPRFR(JI,JK),ZPRFR(JI,JK+KKL))
          IF (ZPRFR(JI,JK)==0) THEN
             ZPRFR(JI,JK)=1.
          END IF
       ELSE
          ZPRFR(JI,JK)=0.
       END IF
    END DO
  END DO
!
IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:RAINFR_VERT',1,ZHOOK_HANDLE)
!
END SUBROUTINE RAINFR_VERT
!
!
!-------------------------------------------------------------------------------
!
!
SUBROUTINE COUNTJV(IC, LTAB, I1, I2)

IMPLICIT NONE

INTEGER, INTENT(OUT) :: IC
LOGICAL, DIMENSION(:,:), INTENT(IN) :: LTAB ! Mask
INTEGER, DIMENSION(:), INTENT(OUT) :: I1,I2 ! Used to replace the COUNT and PACK
INTEGER :: JI,JK
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:COUNTJV',0,ZHOOK_HANDLE)

IC = 0
DO JK = 1,SIZE(LTAB,2)
  DO JI = 1,SIZE(LTAB,1)
    IF(LTAB(JI,JK)) THEN
      IC = IC +1
      I1(IC) = JI
      I2(IC) = JK
    END IF
  END DO
END DO

IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:COUNTJV',1,ZHOOK_HANDLE)

END SUBROUTINE COUNTJV

SUBROUTINE COUNTJV2(IC, LTAB, I1)
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
!
INTEGER, INTENT(OUT) :: IC
LOGICAL, DIMENSION(:), INTENT(IN) :: LTAB ! Mask
INTEGER, DIMENSION(:), INTENT(OUT) :: I1   ! Used to replace the COUNT and PACK
INTEGER :: JI
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:COUNTJV2',0,ZHOOK_HANDLE)

IC = 0
DO JI = 1,SIZE(LTAB,1)
  IF( LTAB(JI) ) THEN
    IC = IC +1
    I1(IC) = JI
  END IF
END DO
!
IF (LHOOK) CALL DR_HOOK('RAIN_ICE_OLD:COUNTJV2',1,ZHOOK_HANDLE)
END SUBROUTINE COUNTJV2
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE RAIN_ICE_OLD
