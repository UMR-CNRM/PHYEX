!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
      SUBROUTINE RAIN_ICE ( D, CST, PARAMI, ICEP, ICED, BUCONF,                   &
                            KPROMA, KSIZE,                                        &
                            OCND2,HSUBG_AUCV_RC, HSUBG_AUCV_RI,  &
                            PTSTEP, KRR, ODMICRO, PEXN,                           &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                            PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,               &
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                   &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                            PINPRC, PINPRR, PEVAP3D,                              &
                            PINPRS, PINPRG, PINDEP, PRAINFR, PSIGS,               &
                            TBUDGETS, KBUDGETS,                                   &
                            PSEA, PTOWN,                                          &
                            PRHT, PRHS, PINPRH, PFPR                              )
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
!!      (C. Lac) FIT temporal scheme : instant M removed
!!      (JP Pinty), 01-2014 : ICE4 : partial reconversion of hail to graupel
!!              July, 2015 (O.Nuissier/F.Duffourg) Add microphysics diagnostic for
!!                                      aircraft, ballon and profiler
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!      C.Lac : 10/2016 : add droplet deposition
!!      C.Lac : 01/2017 : correction on droplet deposition
!!      J.Escobar : 10/2017 : for real*4 , limit exp() in RAIN_ICE_SLOW with XMNH_HUGE_12_LOG
!!      (C. Abiven, Y. Léauté, V. Seigner, S. Riette) Phasing of Turner rain subgrid param
!!      (S. Riette) Source code split into several files
!!                  02/2019 C.Lac add rain fraction as an output field
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet 29/05/2019: remove PACK/UNPACK intrinsics (to get more performance and better OpenACC support)
!  P. Wautelet 17/01/2020: move Quicksort to tools.f90
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!  P. Wautelet 25/02/2020: bugfix: add missing budget: WETH_BU_RRG
!!     R. El Khatib 24-Aug-2021 Optimizations
!  J. Wurtz       03/2022: New snow characteristics with LSNOW_T
!-----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_BUDGET,         ONLY: TBUDGETDATA, TBUDGETCONF_t, NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, &
                               NBUDGET_RI, NBUDGET_RR, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE MODD_FIELDS_ADDRESS, ONLY : & ! common fields adress
      & ITH,     & ! Potential temperature
      & IRV,     & ! Water vapor
      & IRC,     & ! Cloud water
      & IRR,     & ! Rain water
      & IRI,     & ! Pristine ice
      & IRS,     & ! Snow/aggregate
      & IRG,     & ! Graupel
      & IRH,     & ! Hail
      & IBUNUM,       & ! Number of tendency terms
      & IBUNUM_MR,    & ! Number of tendency terms expressed as mixing ratio changes
      & IBUNUM_EXTRA, & ! Number of extra tendency terms
      & IRREVAV,      & ! Index for the evaporation tendency
      & IBUEXTRAIND     ! Index indirection

USE MODE_BUDGET,         ONLY: BUDGET_STORE_ADD_PHY, BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY
USE MODE_MSG,            ONLY: PRINT_MSG, NVERB_FATAL

USE MODE_ICE4_BUDGETS, ONLY: ICE4_BUDGETS
USE MODE_ICE4_RAINFR_VERT, ONLY: ICE4_RAINFR_VERT
USE MODE_ICE4_SEDIMENTATION, ONLY: ICE4_SEDIMENTATION
USE MODE_ICE4_TENDENCIES, ONLY: ICE4_TENDENCIES
USE MODE_ICE4_NUCLEATION, ONLY: ICE4_NUCLEATION
USE MODE_ICE4_CORRECT_NEGATIVITIES, ONLY: ICE4_CORRECT_NEGATIVITIES
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),        INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
INTEGER,                  INTENT(IN)    :: KPROMA ! cache-blocking factor for microphysic loop
INTEGER,                  INTENT(IN)    :: KSIZE
LOGICAL                                 :: OCND2  ! Logical switch to separate liquid and ice
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Kind of Subgrid autoconversion method
CHARACTER(LEN=80),        INTENT(IN)    :: HSUBG_AUCV_RI ! Kind of Subgrid autoconversion method
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: ODMICRO ! mask to limit computation
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PHLC_HRC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PHLC_HCF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PHLI_HRI
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PHLI_HCF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(D%NIJT), INTENT(OUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(D%NIJT), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(D%NIJT), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(D%NIJT), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(MERGE(D%NIJT, 0, PARAMI%LDEPOSC)), INTENT(OUT) :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) :: PRAINFR !Precipitation fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
REAL, DIMENSION(D%NIJT), OPTIONAL, INTENT(IN)        :: PSEA ! Sea Mask
REAL, DIMENSION(D%NIJT), OPTIONAL, INTENT(IN)        :: PTOWN! Fraction that is town
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL,  INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL,  INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(D%NIJT), OPTIONAL, INTENT(OUT)      :: PINPRH! Hail instant precip
REAL, DIMENSION(D%NIJT,D%NKT,KRR), OPTIONAL, INTENT(OUT)  :: PFPR ! upper-air precipitation fluxes
!
!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
INTEGER :: JIJ, JK
INTEGER :: IKTB, IKTE, IKB, IIJB, IIJE
INTEGER :: ISTIJ, ISTK
!
!Arrays for nucleation call outisde of ODMICRO points
REAL,    DIMENSION(D%NIJT, D%NKT) :: ZW ! work array
REAL,    DIMENSION(D%NIJT, D%NKT) :: ZT ! Temperature
REAL, DIMENSION(D%NIJT, D%NKT) :: ZZ_RVHENI_MR, & ! heterogeneous nucleation mixing ratio change
                                & ZZ_RVHENI       ! heterogeneous nucleation
REAL, DIMENSION(D%NIJT, D%NKT) :: ZZ_LVFACT, ZZ_LSFACT, ZZ_DIFF
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRCT    ! Cloud water m.r. source at t
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRRT    ! Rain water m.r. source at t
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRIT    ! Pristine ice m.r. source at t
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRST    ! Snow/aggregate m.r. source at t
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRGT    ! Graupel m.r. source at t
REAL, DIMENSION(D%NIJT,D%NKT) :: ZRHT    ! Hail m.r. source at t
REAL, DIMENSION(D%NIJT,D%NKT) :: ZCITOUT ! Output value for CIT

!
LOGICAL :: GEXT_TEND
LOGICAL :: LSOFT ! Must we really compute tendencies or only adjust them to new T variables
INTEGER :: INB_ITER_MAX ! Maximum number of iterations (with real tendencies computation)
REAL :: ZW0D
REAL :: ZTSTEP ! length of sub-timestep in case of time splitting
REAL :: ZINV_TSTEP ! Inverse ov PTSTEP
REAL :: ZTIME_THRESHOLD ! Time to reach threshold
!For total tendencies computation
REAL, DIMENSION(D%NIJT,D%NKT,0:7) :: ZWR
!
!Output packed total mixing ratio change (for budgets only)
REAL, DIMENSION(KSIZE, IBUNUM-IBUNUM_EXTRA) :: ZBU_PACK
!
!For packing
INTEGER :: IMICRO ! Case r_x>0 locations
INTEGER :: JL, JV, JJV
REAL, DIMENSION(KPROMA) :: ZTIME ! Current integration time (starts with 0 and ends with PTSTEP)
REAL, DIMENSION(KPROMA) :: &
                        & ZMAXTIME, & ! Time on which we can apply the current tendencies
                        & ZTIME_LASTCALL, &     ! Integration time when last tendecies call has been done
                        & ZSSI,     &
                        & ZCIT,     & ! Pristine ice conc. at t
                        & ZRHODREF, & ! RHO Dry REFerence
                        & ZZT,      & ! Temperature
                        & ZPRES,    & ! Pressure
                        & ZEXN,     & ! EXNer Pressure
                        & ZLSFACT,  & ! L_s/(Pi*C_ph)
                        & ZLVFACT,  & ! L_v/(Pi*C_ph)
                        & ZSIGMA_RC,& ! Standard deviation of rc at time t
                        & ZCF,      & ! Cloud fraction
                        & ZHLC_HCF, & ! HLCLOUDS : fraction of High Cloud Fraction in grid
                        & ZHLC_LCF, & ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                      !    note that ZCF = ZHLC_HCF + ZHLC_LCF
                        & ZHLC_HRC, & ! HLCLOUDS : LWC that is High LWC in grid
                        & ZHLC_LRC, & ! HLCLOUDS : LWC that is Low  LWC in grid
                                      !    note that ZRC = ZHLC_HRC + ZHLC_LRC
                        & ZHLI_HCF, &
                        & ZHLI_LCF, &
                        & ZHLI_HRI, &
                        & ZHLI_LRI
LOGICAL, DIMENSION(KPROMA) :: LLCOMPUTE ! .TRUE. or points where we must compute tendencies,
!
!Output packed tendencies (for budgets only)
REAL, DIMENSION(KPROMA, IBUNUM-IBUNUM_EXTRA) :: ZBU_SUM
REAL, DIMENSION(KPROMA, IBUNUM) :: ZBU_INST
!
!For mixing-ratio-splitting
LOGICAL :: LLCPZ0RT
REAL :: ZTIME_THRESHOLD1D(KPROMA) ! Time to reach threshold
REAL, DIMENSION(KPROMA, KRR) :: Z0RT ! Mixing-ratios at the beginig of the current loop
!
REAL, DIMENSION(KPROMA,0:7) :: &
                        & ZVART, & !Packed variables
                        & ZEXTPK, & !To take into acount external tendencies inside the splitting
                        & ZA, ZB
!
REAL, DIMENSION(KPROMA, 8) :: ZRS_TEND, ZRG_TEND
REAL, DIMENSION(KPROMA,10) :: ZRH_TEND

INTEGER, DIMENSION(KPROMA) :: &
                       & I1,I2, & ! Used to replace the COUNT and PACK intrinsics on variables
                       & IITER    ! Number of iterations done (with real tendencies computation)
INTEGER, DIMENSION(KSIZE) :: I1TOT, I2TOT ! Used to replace the COUNT and PACK intrinsics
!
REAL, DIMENSION(KPROMA) :: ZSUM2, ZMAXB
REAL :: ZDEVIDE, ZX, ZRICE
!
INTEGER :: IC, JMICRO
LOGICAL :: LLSIGMA_RC, LL_ANY_ITER, LL_AUCV_ADJU
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZW3D
LOGICAL, DIMENSION(D%NIJT,D%NKT) :: LLW3D
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RAIN_ICE', 0, ZHOOK_HANDLE)
!
IKTB=D%NKTB
IKTE=D%NKTE
IKB=D%NKB
IIJB=D%NIJB
IIJE=D%NIJE
!-------------------------------------------------------------------------------
!
IF(OCND2) THEN
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'RAIN_ICE', 'OCND2 OPTION NOT CODED IN THIS RAIN_ICE VERSION')
END IF
IF(KPROMA /= KSIZE) THEN
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'RAIN_ICE', 'For now, KPROMA must be equal to KSIZE, see code for explanation')
  ! Microphyscs was optimized by introducing chunks of KPROMA size
  ! Thus, in ice4_tendencies, the 1D array represent only a fraction of the points where microphisical species are present
  ! We cannot rebuild the entire 3D arrays in the subroutine, so we cannot call ice4_rainfr_vert in it
  ! A solution would be to suppress optimisation in this case by setting KPROMA=KSIZE in rain_ice
  ! Another solution would be to compute column by column?
  ! Another one would be to cut tendencies in 3 parts: before rainfr_vert, rainfr_vert, after rainfr_vert
ENDIF
!
!*       1.     COMPUTE THE LOOP BOUNDS
!               -----------------------
!
ZINV_TSTEP=1./PTSTEP
GEXT_TEND=.TRUE.
!
! LSFACT and LVFACT without exner
DO JK = IKTB,IKTE
  DO JIJ = IIJB,IIJE
    IF (KRR==7) THEN
      ZRICE=PRIT(JIJ,JK)+PRST(JIJ,JK)+PRGT(JIJ,JK)+PRHT(JIJ,JK)
    ELSE
      ZRICE=PRIT(JIJ,JK)+PRST(JIJ,JK)+PRGT(JIJ,JK)
    ENDIF
    ZDEVIDE = CST%XCPD + CST%XCPV*PRVT(JIJ,JK) + CST%XCL*(PRCT(JIJ,JK)+PRRT(JIJ,JK)) + CST%XCI*ZRICE
    ZT(JIJ,JK) = PTHT(JIJ,JK) * PEXN(JIJ,JK)
    ZZ_LSFACT(JIJ,JK)=(CST%XLSTT+(CST%XCPV-CST%XCI)*(ZT(JIJ,JK)-CST%XTT)) / ZDEVIDE
    ZZ_LVFACT(JIJ,JK)=(CST%XLVTT+(CST%XCPV-CST%XCL)*(ZT(JIJ,JK)-CST%XTT)) / ZDEVIDE
  ENDDO
ENDDO
!
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!               -------------------------------------
!
IF(.NOT. PARAMI%LSEDIM_AFTER) THEN
  IF(KRR==7) THEN
    CALL ICE4_SEDIMENTATION(D, CST, ICEP, ICED, PARAMI, BUCONF, &
                           &PTSTEP, KRR, PDZZ, &
                           &ZZ_LVFACT, ZZ_LSFACT, PRHODREF, PPABST, PTHT, ZT, PRHODJ, &
                           &PTHS, PRVS, PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                           &PINPRC, PINPRR, PINPRS, PINPRG, &
                           &TBUDGETS, KBUDGETS, &
                           &PSEA=PSEA, PTOWN=PTOWN, &
                           &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR)
  ELSE
    CALL ICE4_SEDIMENTATION(D, CST, ICEP, ICED, PARAMI, BUCONF, &
                           &PTSTEP, KRR, PDZZ, &
                           &ZZ_LVFACT, ZZ_LSFACT, PRHODREF, PPABST, PTHT, ZT, PRHODJ, &
                           &PTHS, PRVS, PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                           &PINPRC, PINPRR, PINPRS, PINPRG, &
                           &TBUDGETS, KBUDGETS, &
                           &PSEA=PSEA, PTOWN=PTOWN, &
                           &PFPR=PFPR)
  ENDIF
ENDIF
!

DO JK = IKTB,IKTE
  !Backup of T variables
  ZWR(:,JK,IRV)=PRVT(:,JK)
  ZWR(:,JK,IRC)=PRCT(:,JK)
  ZWR(:,JK,IRR)=PRRT(:,JK)
  ZWR(:,JK,IRI)=PRIT(:,JK)
  ZWR(:,JK,IRS)=PRST(:,JK)
  ZWR(:,JK,IRG)=PRGT(:,JK)
  IF (KRR==7) THEN
    ZWR(:,JK,IRH)=PRHT(:,JK)
  ELSE
    ZWR(:,JK,IRH)=0.
  ENDIF

  !Preset for output 3D variables
  IF(PARAMI%LWARM) THEN
    PEVAP3D(:,JK)=0.
  ENDIF
  PRAINFR(:,JK)=0.
#ifdef REPRO55
  ZCITOUT(:,JK)=PCIT(:,JK)
#else
  ZCITOUT(:,JK)=0. !We want 0 outside of IMICRO points
#endif
ENDDO

IF(BUCONF%LBU_ENABLE) THEN
  DO JV=1, IBUNUM-IBUNUM_EXTRA
    ZBU_PACK(:, JV)=0.
    ZBU_SUM(:, JV)=0.
  ENDDO
ENDIF

!-------------------------------------------------------------------------------
!  optimization by looking for locations where
!  the microphysical fields are larger than a minimal value only !!!
!
IF (KSIZE /= COUNT(ODMICRO(IIJB:IIJE,IKTB:IKTE))) THEN
    CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'RAIN_ICE', 'RAIN_ICE : KSIZE /= COUNT(ODMICRO)')
ENDIF

IF (KSIZE > 0) THEN

  !Maximum number of iterations
  !We only count real iterations (those for which we *compute* tendencies)
  INB_ITER_MAX=PARAMI%NMAXITER
  IF(PARAMI%XTSTEP_TS/=0.)THEN
    INB_ITER_MAX=MAX(1, INT(PTSTEP/PARAMI%XTSTEP_TS)) !At least the number of iterations needed for the time-splitting
    ZTSTEP=PTSTEP/INB_ITER_MAX
    INB_ITER_MAX=MAX(PARAMI%NMAXITER, INB_ITER_MAX) !For the case XMRSTEP/=0. at the same time
  ENDIF

!===============================================================================================================
! Cache-blocking loop :

  LLSIGMA_RC=(HSUBG_AUCV_RC=='PDF ' .AND. PARAMI%CSUBG_PR_PDF=='SIGM')
  LL_AUCV_ADJU=(HSUBG_AUCV_RC=='ADJU' .OR. HSUBG_AUCV_RI=='ADJU')

  ! starting indexes :
  IC=0
  ISTK=IKTB
  ISTIJ=IIJB

  DO JMICRO=1,KSIZE,KPROMA

    IMICRO=MIN(KPROMA,KSIZE-JMICRO+1)
!
!*       3.     PACKING
!               --------

    ! Setup packing parameters
    OUTER_LOOP: DO JK = ISTK, IKTE
      IF (ANY(ODMICRO(:,JK))) THEN
        DO JIJ = ISTIJ, IIJE
          IF (ODMICRO(JIJ,JK)) THEN
            IC=IC+1
            ! Initialization of variables in packed format :
            ZVART(IC, ITH)=PTHT(JIJ, JK)
            ZVART(IC, IRV)=PRVT(JIJ, JK)
            ZVART(IC, IRC)=PRCT(JIJ, JK)
            ZVART(IC, IRR)=PRRT(JIJ, JK)
            ZVART(IC, IRI)=PRIT(JIJ, JK)
            ZVART(IC, IRS)=PRST(JIJ, JK)
            ZVART(IC, IRG)=PRGT(JIJ, JK)
            IF (KRR==7) THEN
              ZVART(IC, IRH)=PRHT(JIJ, JK)
            ENDIF
            IF (GEXT_TEND) THEN
              !The th tendency is not related to a mixing ratio change, there is no exn/exnref issue here
              ZEXTPK(IC, ITH)=PTHS(JIJ, JK)
              ZEXTPK(IC, IRV)=PRVS(JIJ, JK)
              ZEXTPK(IC, IRC)=PRCS(JIJ, JK)
              ZEXTPK(IC, IRR)=PRRS(JIJ, JK)
              ZEXTPK(IC, IRI)=PRIS(JIJ, JK)
              ZEXTPK(IC, IRS)=PRSS(JIJ, JK)
              ZEXTPK(IC, IRG)=PRGS(JIJ, JK)
              IF (KRR==7) THEN
                ZEXTPK(IC, IRH)=PRHS(JIJ, JK)
              ENDIF
            ENDIF
            ZCIT       (IC)=PCIT    (JIJ, JK)
            ZCF        (IC)=PCLDFR  (JIJ, JK)
            ZRHODREF   (IC)=PRHODREF(JIJ, JK)
            ZPRES      (IC)=PPABST  (JIJ, JK)
            ZEXN       (IC)=PEXN    (JIJ, JK)
            IF(LLSIGMA_RC) THEN
              ZSIGMA_RC(IC)=PSIGS   (JIJ, JK)
            ENDIF
            IF (LL_AUCV_ADJU) THEN
              ZHLC_HCF(IC) = PHLC_HCF(JIJ, JK)
              ZHLC_HRC(IC) = PHLC_HRC(JIJ, JK)
              ZHLI_HCF(IC) = PHLI_HCF(JIJ, JK)
              ZHLI_HRI(IC) = PHLI_HRI(JIJ, JK)
            ENDIF
            ! Save indices for later usages:
            I1(IC) = JIJ
            I2(IC) = JK
            I1TOT(JMICRO+IC-1)=JIJ
            I2TOT(JMICRO+IC-1)=JK
            IF (IC==IMICRO) THEN
              ! the end of the chunk has been reached, then reset the starting index :
              ISTIJ=JIJ+1
              IF (ISTIJ <= IIJE) THEN
                ISTK=JK
              ELSE
                ! end of line, restart from 1 and increment upper loop
                ISTK=JK+1
                IF (ISTK > IKTE) THEN
                  ! end of line, restart from 1
                  ISTK=IKTB
                ENDIF
              ENDIF
              IC=0
              EXIT OUTER_LOOP
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      ! restart inner loop on JIJ :
      ISTIJ=IIJB
    ENDDO OUTER_LOOP

    IF (GEXT_TEND) THEN
      DO JV=0, KRR
        DO JL=1, IMICRO
          ZEXTPK(JL, JV)=ZEXTPK(JL, JV)-ZVART(JL, JV)*ZINV_TSTEP
        ENDDO
      ENDDO
    ENDIF
    IF (LLSIGMA_RC) THEN
      DO JL=1, IMICRO
        ZSIGMA_RC(JL)=ZSIGMA_RC(JL)*2.
      ENDDO 
    ENDIF
    IF (LL_AUCV_ADJU) THEN
      DO JL=1, IMICRO
        ZHLC_LRC(JL) = ZVART(JL, IRC) - ZHLC_HRC(JL)
        ZHLI_LRI(JL) = ZVART(JL, IRI) - ZHLI_HRI(JL)
        IF(ZVART(JL, IRC)>0.) THEN
          ZHLC_LCF(JL) = ZCF(JL)- ZHLC_HCF(JL)
        ELSE
          ZHLC_LCF(JL)=0.
        ENDIF
        IF(ZVART(JL, IRI)>0.) THEN
          ZHLI_LCF(JL) = ZCF(JL)- ZHLI_HCF(JL)
        ELSE
          ZHLI_LCF(JL)=0.
        ENDIF
      ENDDO
    ENDIF

!-------------------------------------------------------------------------------
!
!*       4.     LOOP
!               ----
!
    IITER(1:IMICRO)=0
    ZTIME(1:IMICRO)=0. ! Current integration time (all points may have a different integration time)

    DO WHILE(ANY(ZTIME(1:IMICRO)<PTSTEP)) ! Loop to *really* compute tendencies

      IF(PARAMI%XTSTEP_TS/=0.) THEN
        ! In this case we need to remember the time when tendencies were computed
        ! because when time has evolved more than a limit, we must re-compute tendencies
        ZTIME_LASTCALL(1:IMICRO)=ZTIME(1:IMICRO)
      ENDIF
      DO JL=1, IMICRO
        IF (ZTIME(JL) < PTSTEP) THEN
          LLCOMPUTE(JL)=.TRUE. ! Computation (.TRUE.) only for points for which integration time has not reached the timestep
          IITER(JL)=IITER(JL)+1
        ELSE
          LLCOMPUTE(JL)=.FALSE.
        ENDIF
      ENDDO
      LL_ANY_ITER=ANY(IITER(1:IMICRO) < INB_ITER_MAX)
      LLCPZ0RT=.TRUE.
      LSOFT=.FALSE. ! We *really* compute the tendencies

      DO WHILE(ANY(LLCOMPUTE(1:IMICRO))) ! Loop to adjust tendencies when we cross the 0°C or when a species disappears
!$OMP SIMD
        DO JL=1, IMICRO
          ZSUM2(JL)=SUM(ZVART(JL,IRI:KRR))
        ENDDO
        DO JL=1, IMICRO
          ZDEVIDE=(CST%XCPD + CST%XCPV*ZVART(JL, IRV) + CST%XCL*(ZVART(JL, IRC)+ZVART(JL, IRR)) + CST%XCI*ZSUM2(JL)) * ZEXN(JL)
          ZZT(JL) = ZVART(JL, ITH) * ZEXN(JL)
          ZLSFACT(JL)=(CST%XLSTT+(CST%XCPV-CST%XCI)*(ZZT(JL)-CST%XTT)) / ZDEVIDE
          ZLVFACT(JL)=(CST%XLVTT+(CST%XCPV-CST%XCL)*(ZZT(JL)-CST%XTT)) / ZDEVIDE
        ENDDO
        !
        !***       4.1 Tendencies computation
        !
        ! Tendencies are *really* computed when LSOFT==.FALSE. and only adjusted otherwise
    CALL ICE4_TENDENCIES(D, CST, PARAMI, ICEP, ICED, BUCONF, &
                        &KPROMA, IMICRO, &
                        &KRR, LSOFT, LLCOMPUTE, &
                        &HSUBG_AUCV_RC, HSUBG_AUCV_RI, &
                        &ZEXN, ZRHODREF, ZLVFACT, ZLSFACT, I1, I2, &
                        &ZPRES, ZCF, ZSIGMA_RC, &
                        &ZCIT, &
                        &ZZT, ZVART, &
                        &ZBU_INST, &
                        &ZRS_TEND, ZRG_TEND, ZRH_TEND, ZSSI, &
                        &ZA, ZB, &
                        &ZHLC_HCF, ZHLC_LCF, ZHLC_HRC, ZHLC_LRC, &
                        &ZHLI_HCF, ZHLI_LCF, ZHLI_HRI, ZHLI_LRI, PRAINFR)

        ! External tendencies
        IF(GEXT_TEND) THEN
          DO JV=0, KRR
            DO JL=1, IMICRO
              ZA(JL, JV) = ZA(JL, JV) + ZEXTPK(JL, JV)
            ENDDO
          ENDDO
        ENDIF
        !
        !***       4.2 Integration time
        !
        ! If we can, we shall use these tendencies until the end of the timestep
        DO JL=1, IMICRO
          IF(LLCOMPUTE(JL)) THEN
            ZMAXTIME(JL)=(PTSTEP-ZTIME(JL)) ! Remaining time until the end of the timestep
          ELSE
            ZMAXTIME(JL)=0.
          ENDIF
        ENDDO

        !We need to adjust tendencies when temperature reaches 0
        IF(PARAMI%LFEEDBACKT) THEN
          DO JL=1, IMICRO
            !Is ZB(:, ITH) enough to change temperature sign?
            ZX=CST%XTT/ZEXN(JL)
            IF ((ZVART(JL, ITH) - ZX) * (ZVART(JL, ITH) + ZB(JL, ITH) - ZX) < 0.) THEN
              ZMAXTIME(JL)=0.
            ENDIF
            !Can ZA(:, ITH) make temperature change of sign?
            IF (ABS(ZA(JL,ITH)) > 1.E-20 ) THEN
              ZTIME_THRESHOLD=(ZX - ZB(JL, ITH) - ZVART(JL, ITH))/ZA(JL, ITH)
              IF (ZTIME_THRESHOLD > 0.) THEN
                ZMAXTIME(JL)=MIN(ZMAXTIME(JL), ZTIME_THRESHOLD)
              ENDIF
            ENDIF
          ENDDO
        ENDIF

        !We need to adjust tendencies when a species disappears
        !When a species is missing, only the external tendencies can be negative (and we must keep track of it)
        DO JV=1, KRR
          DO JL=1, IMICRO
            IF (ZA(JL, JV) < -1.E-20 .AND. ZVART(JL, JV) > ICED%XRTMIN(JV)) THEN
              ZMAXTIME(JL)=MIN(ZMAXTIME(JL), -(ZB(JL, JV)+ZVART(JL, JV))/ZA(JL, JV))
            ENDIF
          ENDDO
        ENDDO

        !We stop when the end of the timestep is reached
        DO JL=1, IMICRO
          IF (ZTIME(JL)+ZMAXTIME(JL) >= PTSTEP) THEN
            LLCOMPUTE(JL)=.FALSE.
          ENDIF
        ENDDO
        !We must recompute tendencies when the end of the sub-timestep is reached
        IF (PARAMI%XTSTEP_TS/=0.) THEN
          DO JL=1, IMICRO
            IF ((IITER(JL) < INB_ITER_MAX) .AND. (ZTIME(JL)+ZMAXTIME(JL) > ZTIME_LASTCALL(JL)+ZTSTEP)) THEN
              ZMAXTIME(JL)=ZTIME_LASTCALL(JL)-ZTIME(JL)+ZTSTEP
              LLCOMPUTE(JL)=.FALSE.
            ENDIF
          ENDDO
        ENDIF

        !We must recompute tendencies when the maximum allowed change is reached
        !When a species is missing, only the external tendencies can be active and we do not want to recompute
        !the microphysical tendencies when external tendencies are negative (results won't change because species was already missing)
        IF (PARAMI%XMRSTEP/=0.) THEN
          IF (LL_ANY_ITER) THEN
            ! In this case we need to remember the initial mixing ratios used to compute the tendencies
            ! because when mixing ratio has evolved more than a threshold, we must re-compute tendencies
            ! Thus, at first iteration (ie when LLCPZ0RT=.TRUE.) we copy ZVART into Z0RT
            DO JV=1,KRR
              IF (LLCPZ0RT) Z0RT(1:IMICRO, JV)=ZVART(1:IMICRO, JV)
              DO JL=1, IMICRO
                IF (IITER(JL)<INB_ITER_MAX .AND. ABS(ZA(JL,JV))>1.E-20) THEN
                  ZTIME_THRESHOLD1D(JL)=(SIGN(1., ZA(JL, JV))*PARAMI%XMRSTEP+ &
                                        &Z0RT(JL, JV)-ZVART(JL, JV)-ZB(JL, JV))/ZA(JL, JV)
                ELSE
                  ZTIME_THRESHOLD1D(JL)=-1.
                ENDIF
              ENDDO
              DO JL=1, IMICRO
                IF (ZTIME_THRESHOLD1D(JL)>=0 .AND. ZTIME_THRESHOLD1D(JL)<ZMAXTIME(JL) .AND. &
                   &(ZVART(JL, JV)>ICED%XRTMIN(JV) .OR. ZA(JL, JV)>0.)) THEN
                  ZMAXTIME(JL)=MIN(ZMAXTIME(JL), ZTIME_THRESHOLD1D(JL))
                  LLCOMPUTE(JL)=.FALSE.
                ENDIF
              ENDDO
            ENDDO
            LLCPZ0RT=.FALSE.
!$OMP SIMD
            DO JL=1,IMICRO
              ZMAXB(JL)=MAXVAL(ABS(ZB(JL,1:KRR)))
            ENDDO
            DO JL=1, IMICRO
              IF (IITER(JL)<INB_ITER_MAX .AND. ZMAXB(JL)>PARAMI%XMRSTEP) THEN
                ZMAXTIME(JL)=0.
                LLCOMPUTE(JL)=.FALSE.
              ENDIF
            ENDDO
          ENDIF ! LL_ANY_ITER
        ENDIF ! XMRSTEP/=0.
        !
        !***       4.3 New values of variables for next iteration
        !
        DO JV=0, KRR
          DO JL=1, IMICRO
            ZVART(JL, JV)=ZVART(JL, JV)+ZA(JL, JV)*ZMAXTIME(JL)+ZB(JL, JV)
          ENDDO
        ENDDO
        DO JL=1, IMICRO
#ifdef REPRO55
          ZCIT(JL)=ZCIT(JL) * MAX(0., -SIGN(1., -ZVART(JL,IRI)))
#else
          IF (ZVART(JL,IRI)<=0.) ZCIT(JL) = 0.
#endif
          ZTIME(JL)=ZTIME(JL)+ZMAXTIME(JL)
        ENDDO

        !
        !***       4.4 Mixing ratio change due to each process
        !
        IF(BUCONF%LBU_ENABLE) THEN
          !Mixing ratio change due to a tendency
          DO JV=1, IBUNUM-IBUNUM_MR-IBUNUM_EXTRA
            DO JL=1, IMICRO
              ZBU_SUM(JL, JV) = ZBU_SUM(JL, JV) + ZBU_INST(JL, JV)*ZMAXTIME(JL)
            ENDDO
          ENDDO

          !Mixing ratio change due to a mixing ratio change
          DO JV=IBUNUM-IBUNUM_MR-IBUNUM_EXTRA+1, IBUNUM-IBUNUM_EXTRA
            DO JL=1, IMICRO
              ZBU_SUM(JL, JV) = ZBU_SUM(JL, JV) + ZBU_INST(JL, JV)
            ENDDO
          ENDDO

          !Extra contribution as a mixing ratio change
          DO JV=IBUNUM-IBUNUM_EXTRA+1, IBUNUM
            JJV=IBUEXTRAIND(JV)
            DO JL=1, IMICRO
              ZBU_SUM(JL, JJV) = ZBU_SUM(JL, JJV) + ZBU_INST(JL, JV)
            ENDDO
          ENDDO
        ENDIF
        !
        !***       4.5 Next loop
        !
        LSOFT=.TRUE. ! We try to adjust tendencies (inner while loop)
      ENDDO
    ENDDO

    IF(GEXT_TEND) THEN
      !Z..T variables contain the external tendency, we substract it
      DO JV=0, KRR
        DO JL=1, IMICRO
          ZVART(JL, JV) = ZVART(JL, JV) - ZEXTPK(JL, JV) * PTSTEP
        ENDDO
      ENDDO
    ENDIF

!-------------------------------------------------------------------------------
!
!*       5.     UNPACKING DIAGNOSTICS
!               ---------------------
!
    DO JL=1, IMICRO
      ZCITOUT  (I1(JL),I2(JL))=ZCIT   (JL)
      IF(PARAMI%LWARM) THEN
        PEVAP3D(I1(JL),I2(JL))=ZBU_INST(JL, IRREVAV)
      ENDIF
      ZWR(I1(JL),I2(JL),IRV)=ZVART(JL, IRV)
      ZWR(I1(JL),I2(JL),IRC)=ZVART(JL, IRC)
      ZWR(I1(JL),I2(JL),IRR)=ZVART(JL, IRR)
      ZWR(I1(JL),I2(JL),IRI)=ZVART(JL, IRI)
      ZWR(I1(JL),I2(JL),IRS)=ZVART(JL, IRS)
      ZWR(I1(JL),I2(JL),IRG)=ZVART(JL, IRG)
      IF (KRR==7) THEN
        ZWR(I1(JL),I2(JL),IRH)=ZVART(JL, IRH)
      ENDIF
    ENDDO
    IF(BUCONF%LBU_ENABLE) THEN
      DO JV=1, IBUNUM-IBUNUM_EXTRA
        DO JL=1, IMICRO
          ZBU_PACK(JMICRO+JL-1, JV) = ZBU_SUM(JL, JV)
        ENDDO
      ENDDO
    ENDIF


  ENDDO ! JMICRO
ENDIF ! KSIZE > 0
PCIT(:,:)=ZCITOUT(:,:)

!==========================================================================================================


!
!*       6.     COMPUTES THE SLOW COLD PROCESS SOURCES OUTSIDE OF ODMICRO POINTS
!               ----------------------------------------------------------------
!
LLW3D(:,:)=.FALSE.
DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE
    IF (.NOT. ODMICRO(JIJ, JK)) THEN
      LLW3D(JIJ, JK)=.TRUE.
      ZW3D(JIJ, JK)=ZZ_LSFACT(JIJ, JK)/PEXN(JIJ, JK)
    ELSE
      LLW3D(JIJ, JK)=.FALSE.
    ENDIF
  ENDDO
ENDDO
CALL ICE4_NUCLEATION(CST, PARAMI, ICEP, ICED, D%NIJT*D%NKT, LLW3D(:,:), &
                     PTHT(:, :), PPABST(:, :), PRHODREF(:, :), &                                       
                     PEXN(:, :), ZW3D(:, :), ZT(:, :), &                                                           
                     PRVT(:, :), &                                                                                 
                     PCIT(:, :), ZZ_RVHENI_MR(:, :))
!
!-------------------------------------------------------------------------------
!
!*       7.     TOTAL TENDENCIES
!               ----------------
!
!
!***     7.1    total tendencies limited by available species
!
DO JK = IKTB, IKTE
  DO CONCURRENT (JIJ=IIJB:IIJE)
    !LV/LS
    ZZ_LSFACT(JIJ,JK)=ZZ_LSFACT(JIJ,JK)/PEXNREF(JIJ,JK)
    ZZ_LVFACT(JIJ,JK)=ZZ_LVFACT(JIJ,JK)/PEXNREF(JIJ,JK)

    !Tendency dure to nucleation on non ODMICRO points
    ZZ_RVHENI(JIJ,JK) = MIN(PRVS(JIJ,JK), ZZ_RVHENI_MR(JIJ,JK)/PTSTEP)

    !Hydrometeor tendencies is the difference between old state and new state (can be negative)
    ZWR(JIJ,JK,IRV)=(ZWR(JIJ,JK,IRV)-PRVT(JIJ,JK))*ZINV_TSTEP
    ZWR(JIJ,JK,IRC)=(ZWR(JIJ,JK,IRC)-PRCT(JIJ,JK))*ZINV_TSTEP
    ZWR(JIJ,JK,IRR)=(ZWR(JIJ,JK,IRR)-PRRT(JIJ,JK))*ZINV_TSTEP
    ZWR(JIJ,JK,IRI)=(ZWR(JIJ,JK,IRI)-PRIT(JIJ,JK))*ZINV_TSTEP
    ZWR(JIJ,JK,IRS)=(ZWR(JIJ,JK,IRS)-PRST(JIJ,JK))*ZINV_TSTEP
    ZWR(JIJ,JK,IRG)=(ZWR(JIJ,JK,IRG)-PRGT(JIJ,JK))*ZINV_TSTEP
    IF(KRR==7) THEN
      ZWR(JIJ,JK,IRH)=(ZWR(JIJ,JK,IRH)-PRHT(JIJ,JK))*ZINV_TSTEP
    ENDIF

    !Theta tendency computed from hydrometeors tendencies
    ZWR(JIJ,JK, ITH) = (ZWR(JIJ,JK,IRC)+ZWR(JIJ,JK,IRR))*ZZ_LVFACT(JIJ,JK)+ &
                       & (ZWR(JIJ,JK,IRI)+ZWR(JIJ,JK,IRS)+ZWR(JIJ,JK,IRG)+ &
                       &  ZWR(JIJ,JK,IRH))*ZZ_LSFACT(JIJ,JK)

    !We apply these tendencies to the S variables
    !including the nucleation part
    PTHS(JIJ,JK) = PTHS(JIJ,JK) + ZWR(JIJ,JK,ITH)+ZZ_RVHENI(JIJ,JK)*ZZ_LSFACT(JIJ,JK)
    PRVS(JIJ,JK) = PRVS(JIJ,JK) + ZWR(JIJ,JK,IRV)-ZZ_RVHENI(JIJ,JK)
    PRCS(JIJ,JK) = PRCS(JIJ,JK) + ZWR(JIJ,JK,IRC)
    PRRS(JIJ,JK) = PRRS(JIJ,JK) + ZWR(JIJ,JK,IRR)
    PRIS(JIJ,JK) = PRIS(JIJ,JK) + ZWR(JIJ,JK,IRI)+ZZ_RVHENI(JIJ,JK)
    PRSS(JIJ,JK) = PRSS(JIJ,JK) + ZWR(JIJ,JK,IRS)
    PRGS(JIJ,JK) = PRGS(JIJ,JK) + ZWR(JIJ,JK,IRG)
    IF (KRR==7) THEN
      PRHS(JIJ,JK) = PRHS(JIJ,JK) + ZWR(JIJ,JK,IRH)
    ENDIF
  ENDDO
ENDDO

!
!***     7.2    LBU_ENABLE case
!
IF(BUCONF%LBU_ENABLE) THEN
  !Budgets for the different processes
  !They have been put in a subroutine because they perform pack/unpack
  CALL ICE4_BUDGETS(D, PARAMI, BUCONF, KSIZE, KPROMA, PTSTEP, KRR, I1TOT, I2TOT, &
                    ZZ_LVFACT, ZZ_LSFACT, PRHODJ, &
                    ZZ_RVHENI, ZBU_PACK, &
                    TBUDGETS, KBUDGETS)
  !
  !***     7.3    Final tendencies
  !
  IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'CORR', PTHS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'CORR', PRVS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'CORR', PRCS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RR), 'CORR', PRSS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'CORR', PRIS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RS), 'CORR', PRSS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RG), 'CORR', PRGS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RH .AND. KRR==7) CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RH), 'CORR', PRHS(:, :)*PRHODJ(:, :))
END IF

!NOTE:
!  This call cannot be moved before the preeceding budget calls because,
!  with AROME, the BUDGET_STORE_INIT does nothing. The equivalent is done only
!  once before the physics call and copies of the S variables evolve automatically
!  internally to the budget (DDH) machinery at each BUDGET_STORE_ADD and
!  BUDGET_STORE_END calls. Thus, the difference between the DDH internal version
!  of the S variables and the S variables used in the folowing BUDGET_STORE_END
!  call must only be due to the correction of negativities.
!
!We correct negativities with conservation
CALL ICE4_CORRECT_NEGATIVITIES(D, ICED, KRR, PRVS, PRCS, PRRS, &
                              &PRIS, PRSS, PRGS, &
                              &PTHS, ZZ_LVFACT, ZZ_LSFACT, PRHS)

IF (BUCONF%LBU_ENABLE) THEN
  IF (BUCONF%LBUDGET_TH) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'CORR', PTHS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RV) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'CORR', PRVS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RC) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'CORR', PRCS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RR) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RR), 'CORR', PRRS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RI) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'CORR', PRIS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RS) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RS), 'CORR', PRSS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RG) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RG), 'CORR', PRGS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RH .AND. KRR==7) CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RH), 'CORR', PRHS(:, :)*PRHODJ(:, :))
END IF
!
!-------------------------------------------------------------------------------
!
!*       8.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!               -------------------------------------
!
IF(PARAMI%LSEDIM_AFTER) THEN
  IF(KRR==7) THEN
    CALL ICE4_SEDIMENTATION(D, CST, ICEP, ICED, PARAMI, BUCONF, &
                           &PTSTEP, KRR, PDZZ, &
                           &ZZ_LVFACT, ZZ_LSFACT, PRHODREF, PPABST, PTHT, ZT, PRHODJ, &
                           &PTHS, PRVS, PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                           &PINPRC, PINPRR, PINPRS, PINPRG, &
                           &TBUDGETS, KBUDGETS, &
                           &PSEA=PSEA, PTOWN=PTOWN, &
                           &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR)
  ELSE
    CALL ICE4_SEDIMENTATION(D, CST, ICEP, ICED, PARAMI, BUCONF, &
                           &PTSTEP, KRR, PDZZ, &
                           &ZZ_LVFACT, ZZ_LSFACT, PRHODREF, PPABST, PTHT, ZT, PRHODJ, &
                           &PTHS, PRVS, PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                           &PINPRC, PINPRR, PINPRS, PINPRG, &
                           &TBUDGETS, KBUDGETS, &
                           &PSEA=PSEA, PTOWN=PTOWN, &
                           &PFPR=PFPR)
  ENDIF

  !"sedimentation" of rain fraction
  IF (PRESENT(PRHS)) THEN
    CALL ICE4_RAINFR_VERT(D, ICED, PRAINFR, PRRS(:,:)*PTSTEP, &
                       &PRSS(:,:)*PTSTEP, PRGS(:,:)*PTSTEP, PRHS(:,:)*PTSTEP)
  ELSE
    CALL ICE4_RAINFR_VERT(D, ICED, PRAINFR, PRRS(:,:)*PTSTEP, &
                       &PRSS(:,:)*PTSTEP, PRGS(:,:)*PTSTEP)
  ENDIF
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       9.     COMPUTE THE FOG DEPOSITION TERM 
!               -------------------------------------
!
IF (PARAMI%LDEPOSC) THEN !cloud water deposition on vegetation
  IF (BUCONF%LBU_ENABLE .AND. BUCONF%LBUDGET_RC) &
     & CALL BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'DEPO', PRCS(:, :)*PRHODJ(:, :))

  PINDEP(:)=0.
!DEC$ IVDEP
  DO JIJ = IIJB, IIJE
    PINDEP(JIJ) = PARAMI%XVDEPOSC * PRCT(JIJ, IKB) * PRHODREF(JIJ, IKB) / CST%XRHOLW
    PRCS(JIJ, IKB) = PRCS(JIJ, IKB) - PARAMI%XVDEPOSC * PRCT(JIJ, IKB) / PDZZ(JIJ, IKB)
    PINPRC(JIJ) = PINPRC(JIJ) + PINDEP(JIJ)
  ENDDO

  IF (BUCONF%LBU_ENABLE .AND. BUCONF%LBUDGET_RC) &
     & CALL BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'DEPO', PRCS(:, :)*PRHODJ(:, :))
ENDIF

IF (LHOOK) CALL DR_HOOK('RAIN_ICE', 1, ZHOOK_HANDLE)
!
END SUBROUTINE RAIN_ICE
