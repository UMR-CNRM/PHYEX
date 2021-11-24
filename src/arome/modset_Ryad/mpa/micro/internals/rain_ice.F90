!     ######spl
      SUBROUTINE RAIN_ICE ( KPROMA, KIT, KJT, KKT, KSIZE,                            &
                            OSEDIC, OCND2, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,&
                            PTSTEP, KRR, LDMICRO, PEXN,                           &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                   &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                            PINPRC, PINPRR, PEVAP3D,                              &
                            PINPRS, PINPRG, PSIGS, PSEA, PTOWN,                   &
                            PRHT, PRHS, PINPRH, PFPR,                             &
                            YDDDH, YDLDDH, YDMDDH                    )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
!!      (S. Riette) Source code split into several files
!!     R. El Khatib 24-Aug-2021 Optimizations
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM
USE MODD_PARAM_ICE
USE MODD_BUDGET
USE MODD_LES
USE MODD_FIELDS_ADDRESS, ONLY : & ! common fields adress
      & ITH,     & ! Potential temperature
      & IRV,     & ! Water vapor
      & IRC,     & ! Cloud water
      & IRR,     & ! Rain water
      & IRI,     & ! Pristine ice
      & IRS,     & ! Snow/aggregate
      & IRG,     & ! Graupel
      & IRH        ! Hail
USE MODI_BUDGET
USE MODI_ICE4_RAINFR_VERT
USE MODI_ICE4_SEDIMENTATION_STAT
USE MODI_ICE4_SEDIMENTATION_SPLIT
USE MODI_ICE4_NUCLEATION_WRAPPER
USE MODI_ICE4_TENDENCIES
USE DDH_MIX, ONLY : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
INTEGER,                  INTENT(IN)    :: KPROMA ! cache-blocking factor for microphysic loop
INTEGER,                  INTENT(IN)    :: KIT, KJT, KKT ! arrays size
INTEGER,                  INTENT(IN)    :: KSIZE
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
LOGICAL                                 :: OCND2  ! Logical switch to separate liquid and ice
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)   :: LDMICRO ! mask to limit computation
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PCLDFR  ! Convective Mass Flux Cloud fraction
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(KIT,KJT,KKT), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(IN)        :: PSEA ! Sea Mask
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(IN)        :: PTOWN! Fraction that is town
REAL, DIMENSION(KIT,KJT,KKT), OPTIONAL,  INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), OPTIONAL,  INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(OUT)      :: PINPRH! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)  :: PFPR ! upper-air precipitation fluxes
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!REAL(KIND=JPRB) :: ZHOOK_HANDLE1
!REAL(KIND=JPRB) :: ZHOOK_HANDLE2
!REAL(KIND=JPRB) :: ZHOOK_HANDLE3
!REAL(KIND=JPRB) :: ZHOOK_HANDLE4
!REAL(KIND=JPRB) :: ZHOOK_HANDLE5

INTEGER :: IIB           !  Define the domain where is
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           !
INTEGER :: IJE           !
INTEGER :: IKB, IKTB     !
INTEGER :: IKE, IKTE     !
!
INTEGER :: JI, JJ, JK
INTEGER :: ISTI, ISTJ, ISTK
!
!Arrays for nucleation call outisde of LDMICRO points
REAL,    DIMENSION(KIT, KJT, KKT) :: ZW ! work array
REAL,    DIMENSION(KIT, KJT, KKT) :: ZT ! Temperature
REAL, DIMENSION(KIT, KJT, KKT) :: &
                                  & ZZ_RVHENI_MR ! heterogeneous nucleation mixing ratio change
REAL, DIMENSION(KIT, KJT, KKT) :: ZZ_LVFACT, ZZ_LSFACT
REAL :: ZZ_RVHENI       ! heterogeneous nucleation
!
REAL, DIMENSION(KIT,KJT,KKT) :: ZRCT    ! Cloud water m.r. source at t
REAL, DIMENSION(KIT,KJT,KKT) :: ZRRT    ! Rain water m.r. source at t
REAL, DIMENSION(KIT,KJT,KKT) :: ZRIT    ! Pristine ice m.r. source at t
REAL, DIMENSION(KIT,KJT,KKT) :: ZRST    ! Snow/aggregate m.r. source at t
REAL, DIMENSION(KIT,KJT,KKT) :: ZRGT    ! Graupel m.r. source at t
REAL, DIMENSION(KIT,KJT,KKT) :: ZRHT    ! Hail m.r. source at t

!Diagnostics
LOGICAL :: LLRAIN_FRACTION=.FALSE. ! activate or rain fraction
REAL, DIMENSION(KIT, KJT, KKT) :: ZRAINFR3D
LOGICAL :: LLHLC=.FALSE. ! activate or not HLCLOUDS
REAL, DIMENSION(KIT, KJT, KKT) :: &
                                & ZHLC_HCF3D,& ! HLCLOUDS cloud fraction in high water content part
                                & ZHLC_LCF3D,& ! HLCLOUDS cloud fraction in low water content part
                                & ZHLC_HRC3D,& ! HLCLOUDS cloud water content in high water content
                                & ZHLC_LRC3D   ! HLCLOUDS cloud water content in low water content
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2)) :: ZINPRI ! Pristine ice instant precip
!
LOGICAL :: LEXT_TEND
LOGICAL :: LSOFT ! Must we really compute tendencies or only adjust them to new T variables
INTEGER :: INB_ITER_MAX ! Maximum number of iterations (with real tendencies computation)
REAL :: ZW1D
REAL :: ZTSTEP ! length of sub-timestep in case of time splitting
REAL :: ZINV_TSTEP ! Inverse ov PTSTEP
REAL :: ZTIME_THRESHOLD ! Time to reach threshold
!For total tendencies computation
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3),0:7) :: ZWR
!
!Output packed total mixing ratio change (for budgets only)
REAL, DIMENSION(KSIZE) :: ZTOT_RVHENI, & ! heterogeneous nucleation mixing ratio change
                        & ZTOT_RCHONI, & ! Homogeneous nucleation
                        & ZTOT_RRHONG, & ! Spontaneous freezing mixing ratio change
                        & ZTOT_RVDEPS, & ! Deposition on r_s,
                        & ZTOT_RIAGGS, & ! Aggregation on r_s
                        & ZTOT_RIAUTS, & ! Autoconversion of r_i for r_s production
                        & ZTOT_RVDEPG, & ! Deposition on r_g
                        & ZTOT_RCAUTR,  & ! Autoconversion of r_c for r_r production
                        & ZTOT_RCACCR, & ! Accretion of r_c for r_r production
                        & ZTOT_RREVAV, & ! Evaporation of r_r
                        & ZTOT_RCRIMSS, ZTOT_RCRIMSG, ZTOT_RSRIMCG, & ! Cloud droplet riming of the aggregates
                        & ZTOT_RIMLTC, & ! Cloud ice melting mixing ratio change
                        & ZTOT_RCBERI, & ! Bergeron-Findeisen effect
                        & ZTOT_RHMLTR, & ! Melting of the hailstones
                        & ZTOT_RSMLTG, & ! Conversion-Melting of the aggregates
                        & ZTOT_RCMLTSR, & ! Cloud droplet collection onto aggregates by positive temperature
                        & ZTOT_RRACCSS, ZTOT_RRACCSG, ZTOT_RSACCRG, & ! Rain accretion onto the aggregates
                        & ZTOT_RICFRRG, ZTOT_RRCFRIG, ZTOT_RICFRR, & ! Rain contact freezing
                        & ZTOT_RCWETG, ZTOT_RIWETG, ZTOT_RRWETG, ZTOT_RSWETG, &  ! Graupel wet growth
                        & ZTOT_RCDRYG, ZTOT_RIDRYG, ZTOT_RRDRYG, ZTOT_RSDRYG, &  ! Graupel dry growth
                        & ZTOT_RWETGH, & ! Conversion of graupel into hail
                        & ZTOT_RGMLTR, & ! Melting of the graupel
                        & ZTOT_RCWETH, ZTOT_RIWETH, ZTOT_RSWETH, ZTOT_RGWETH, ZTOT_RRWETH, & ! Dry growth of hailstone
                        & ZTOT_RCDRYH, ZTOT_RIDRYH, ZTOT_RSDRYH, ZTOT_RRDRYH, ZTOT_RGDRYH, & ! Wet growth of hailstone
                        & ZTOT_RDRYHG    ! Conversion of hailstone into graupel
!
!For packing
INTEGER :: IMICRO ! Case r_x>0 locations
INTEGER :: JL, JV
REAL,    DIMENSION(KPROMA) :: ZTIME ! Current integration time (starts with 0 and ends with PTSTEP)
REAL,    DIMENSION(KPROMA) :: &
                        & ZMAXTIME, & ! Time on which we can apply the current tendencies
                        & ZTIME_LASTCALL, &     ! Integration time when last tendecies call has been done
                        & ZCOMPUTE, & ! 1. for points where we must compute tendencies, 0. elsewhere
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
                        & ZRAINFR     ! rain fraction
!
!Output packed tendencies (for budgets only)
REAL, DIMENSION(KPROMA) :: ZRVHENI_MR, & ! heterogeneous nucleation mixing ratio change
                        & ZRCHONI, & ! Homogeneous nucleation
                        & ZRRHONG_MR, & ! Spontaneous freezing mixing ratio change
                        & ZRVDEPS, & ! Deposition on r_s,
                        & ZRIAGGS, & ! Aggregation on r_s
                        & ZRIAUTS, & ! Autoconversion of r_i for r_s production
                        & ZRVDEPG, & ! Deposition on r_g
                        & ZRCAUTR,  & ! Autoconversion of r_c for r_r production
                        & ZRCACCR, & ! Accretion of r_c for r_r production
                        & ZRREVAV, & ! Evaporation of r_r
                        & ZRIMLTC_MR, & ! Cloud ice melting mixing ratio change
                        & ZRCBERI, & ! Bergeron-Findeisen effect
                        & ZRHMLTR, & ! Melting of the hailstones
                        & ZRSMLTG, & ! Conversion-Melting of the aggregates
                        & ZRCMLTSR, & ! Cloud droplet collection onto aggregates by positive temperature
                        & ZRRACCSS, ZRRACCSG, ZRSACCRG, & ! Rain accretion onto the aggregates
                        & ZRCRIMSS, ZRCRIMSG, ZRSRIMCG, ZRSRIMCG_MR, & ! Cloud droplet riming of the aggregates
                        & ZRICFRRG, ZRRCFRIG, ZRICFRR, & ! Rain contact freezing
                        & ZRCWETG, ZRIWETG, ZRRWETG, ZRSWETG, &  ! Graupel wet growth
                        & ZRCDRYG, ZRIDRYG, ZRRDRYG, ZRSDRYG, &  ! Graupel dry growth
                        & ZRWETGH, & ! Conversion of graupel into hail
                        & ZRWETGH_MR, & ! Conversion of graupel into hail, mr change
                        & ZRGMLTR, & ! Melting of the graupel
                        & ZRCWETH, ZRIWETH, ZRSWETH, ZRGWETH, ZRRWETH, & ! Dry growth of hailstone
                        & ZRCDRYH, ZRIDRYH, ZRSDRYH, ZRRDRYH, ZRGDRYH, & ! Wet growth of hailstone
                        & ZRDRYHG    ! Conversion of hailstone into graupel
!
!For time- or mixing-ratio- splitting at the begining of the current loop
LOGICAL, DIMENSION(KRR) :: LLCPZ0T
REAL, DIMENSION(KPROMA) :: Z0T
!
REAL, DIMENSION(KPROMA,0:7) :: &
                        & ZVART, & !Packed variables
                        & ZEXTPK, & !To take into acount external tendencies inside the splitting
                        & ZA, ZB
!
REAL, DIMENSION(KPROMA, 8) :: ZRS_TEND, ZRG_TEND
REAL, DIMENSION(KPROMA,10) :: ZRH_TEND

INTEGER, DIMENSION(KPROMA) :: &
                       & I1,I2,I3, & ! Used to replace the COUNT and PACK intrinsics on variables
                       & IITER       ! Number of iterations done (with real tendencies computation)
!
REAL, DIMENSION(KPROMA) :: ZSUM2, ZMAXB, ZTHRESHOLD
REAL :: ZDEVIDE, ZX, ZSUM1
!
INTEGER :: IOFF, IC, JMICRO
LOGICAL :: LLSIGMA_RC, LL_ANY_ITER

#include "abor1.intfb.h"

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RAIN_ICE', 0, ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!IF (LHOOK) CALL DR_HOOK('RAIN_ICE:ANTE_MICRO', 0, ZHOOK_HANDLE1)
IF(OCND2) THEN
  WRITE(*,*) ' STOP'
  WRITE(*,*) ' OCND2 OPTION NOT CODED IN THIS RAIN_ICE VERSION'
  CALL ABORT
  STOP
END IF
!*       1.     COMPUTE THE LOOP BOUNDS
!               -----------------------
!
IIB=1+JPHEXT
IIE=SIZE(PDZZ,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PDZZ,2) - JPHEXT
IKB=KKA+JPVEXT*KKL
IKE=KKU-JPVEXT*KKL
IKTB=1+JPVEXT
IKTE=KKT-JPVEXT
!
ZINV_TSTEP=1./PTSTEP
LEXT_TEND=.TRUE.
!
! LSFACT and LVFACT without exner
DO JK = 1, KKT
  DO JJ = 1, KJT
    DO JI = 1, KIT
      IF (KRR==7) THEN
        ZSUM1=PRIT(JI,JJ,JK)+PRST(JI,JJ,JK)+PRGT(JI,JJ,JK)+PRHT(JI,JJ,JK)
      ELSE
        ZSUM1=PRIT(JI,JJ,JK)+PRST(JI,JJ,JK)+PRGT(JI,JJ,JK)
      ENDIF
      ZDEVIDE = XCPD + XCPV*PRVT(JI,JJ,JK) + XCL*(PRCT(JI,JJ,JK)+PRRT(JI,JJ,JK)) + XCI*ZSUM1
      ZT(JI,JJ,JK) = PTHT(JI,JJ,JK) * PEXN(JI,JJ,JK)
      ZZ_LSFACT(JI,JJ,JK)=(XLSTT+(XCPV-XCI)*(ZT(JI,JJ,JK)-XTT)) / ZDEVIDE
      ZZ_LVFACT(JI,JJ,JK)=(XLVTT+(XCPV-XCL)*(ZT(JI,JJ,JK)-XTT)) / ZDEVIDE
    ENDDO
  ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!               -------------------------------------
!
IF(.NOT. LSEDIM_AFTER) THEN
  !
  !*       2.1     sedimentation
  !
  IF(HSEDIM=='STAT') THEN
    IF(KRR==7) THEN
      DO JK = 1, KKT
        DO JJ = 1, KJT
          DO JI = 1, KIT
            ZRCT(JI,JJ,JK)=PRCS(JI,JJ,JK)*PTSTEP
            ZRRT(JI,JJ,JK)=PRRS(JI,JJ,JK)*PTSTEP
            ZRIT(JI,JJ,JK)=PRIS(JI,JJ,JK)*PTSTEP
            ZRST(JI,JJ,JK)=PRSS(JI,JJ,JK)*PTSTEP
            ZRGT(JI,JJ,JK)=PRGS(JI,JJ,JK)*PTSTEP
            ZRHT(JI,JJ,JK)=PRHS(JI,JJ,JK)*PTSTEP
          ENDDO
        ENDDO
      ENDDO
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, &
                                  &PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, ZRCT, PRRS, ZRRT, PRIS, ZRIT,&
                                  &PRSS, ZRST, PRGS, ZRGT,&
                                  &PINPRC, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA=PSEA, PTOWN=PTOWN, &
                                  &PINPRH=PINPRH, PRHT=ZRHT, PRHS=PRHS, PFPR=PFPR)
    ELSE
      DO JK = 1, KKT
        DO JJ = 1, KJT
          DO JI = 1, KIT
            ZRCT(JI,JJ,JK)=PRCS(JI,JJ,JK)*PTSTEP
            ZRRT(JI,JJ,JK)=PRRS(JI,JJ,JK)*PTSTEP
            ZRIT(JI,JJ,JK)=PRIS(JI,JJ,JK)*PTSTEP
            ZRST(JI,JJ,JK)=PRSS(JI,JJ,JK)*PTSTEP
            ZRGT(JI,JJ,JK)=PRGS(JI,JJ,JK)*PTSTEP
          ENDDO
        ENDDO
      ENDDO
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, &
                                  &PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, ZRCT, PRRS, ZRRT, PRIS, ZRIT,&
                                  &PRSS, ZRST, PRGS, ZRGT,&
                                  &PINPRC, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA=PSEA, PTOWN=PTOWN, &
                                  &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !No negativity correction here as we apply sedimentation on ZR.T variables
  ELSEIF(HSEDIM=='SPLI') THEN
    !SR: It *seems* that we must have two separate calls for ifort
    IF(KRR==7) THEN
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, &
                                   &PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA=PSEA, PTOWN=PTOWN, &
                                   &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR)
    ELSE
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, &
                                   &PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA=PSEA, PTOWN=PTOWN, &
                                   &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !We correct negativities with conservation
    !SPLI algorith uses a time-splitting. Inside the loop a temporary m.r. is used.
    !   It is initialized with the m.r. at T and is modified by two tendencies:
    !   sedimentation tendency and an external tendency which represents all other
    !   processes (mainly advection and microphysical processes). If both tendencies
    !   are negative, sedimentation can remove a species at a given sub-timestep. From
    !   this point sedimentation stops for the remaining sub-timesteps but the other tendency
    !   will be still active and will lead to negative values.
    !   We could prevent the algorithm to not consume too much a species, instead we apply
    !   a correction here.
    CALL CORRECT_NEGATIVITIES(KIT, KJT, KKT, KRR, PRVS, PRCS, PRRS, &
                             &PRIS, PRSS, PRGS, &
                             &PTHS, ZZ_LVFACT, ZZ_LSFACT, PRHS)
  ELSEIF(HSEDIM=='NONE') THEN
  ELSE
    WRITE(*,*) ' STOP'
    WRITE(*,*) ' NO SEDIMENTATION SCHEME FOR HSEDIM=', HSEDIM
    CALL ABORT
    STOP
  END IF
  !
  !*       2.2     budget storage
  !
  IF (LBUDGET_RC .AND. OSEDIC) &
                  CALL BUDGET (PRCS(:,:,:)*PRHODJ(:,:,:), 7 , 'SEDI_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:)*PRHODJ(:,:,:), 8 , 'SEDI_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET (PRIS(:,:,:)*PRHODJ(:,:,:), 9 , 'SEDI_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET (PRSS(:,:,:)*PRHODJ(:,:,:), 10, 'SEDI_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET (PRGS(:,:,:)*PRHODJ(:,:,:), 11, 'SEDI_BU_RRG',YDDDH, YDLDDH, YDMDDH)
  IF ( KRR == 7 .AND. LBUDGET_RH) &
                  CALL BUDGET (PRHS(:,:,:)*PRHODJ(:,:,:), 12, 'SEDI_BU_RRH',YDDDH, YDLDDH, YDMDDH)
ENDIF
!

! Bakup of tendencies + PEVAP3D preset :
DO JK = 1,KKT
  ZWR(:,:,JK,IRV)=PRVT(:,:,JK)
  ZWR(:,:,JK,IRC)=PRCT(:,:,JK)
  ZWR(:,:,JK,IRR)=PRRT(:,:,JK)
  ZWR(:,:,JK,IRI)=PRIT(:,:,JK)
  ZWR(:,:,JK,IRS)=PRST(:,:,JK)
  ZWR(:,:,JK,IRG)=PRGT(:,:,JK)
  IF (KRR==7) THEN
    ZWR(:,:,JK,IRH)=PRHT(:,:,JK)
  ELSE
    ZWR(:,:,JK,IRH)=0.
  ENDIF
  IF(OWARM) THEN
    PEVAP3D(:,:,JK)=0.
  ENDIF
  IF (KSIZE==0) THEN
    PCIT(:,:,JK) = 0.
  ENDIF
  IF (LLHLC) THEN
    ZHLC_HCF3D(:,:,JK)=0.
    ZHLC_LCF3D(:,:,JK)=0.
    ZHLC_HRC3D(:,:,JK)=0.
    ZHLC_LRC3D(:,:,JK)=0.
  ENDIF
ENDDO

IF(LBU_ENABLE) THEN
  ZTOT_RVHENI(:)=0.
  ZTOT_RCHONI(:)=0.
  ZTOT_RRHONG(:)=0.
  ZTOT_RVDEPS(:)=0.
  ZTOT_RIAGGS(:)=0.
  ZTOT_RIAUTS(:)=0.
  ZTOT_RVDEPG(:)=0.
  ZTOT_RCAUTR(:)=0.
  ZTOT_RCACCR(:)=0.
  ZTOT_RREVAV(:)=0.
  ZTOT_RCRIMSS(:)=0.
  ZTOT_RCRIMSG(:)=0.
  ZTOT_RSRIMCG(:)=0.
  ZTOT_RIMLTC(:)=0.
  ZTOT_RCBERI(:)=0.
  ZTOT_RHMLTR(:)=0.
  ZTOT_RSMLTG(:)=0.
  ZTOT_RCMLTSR(:)=0.
  ZTOT_RRACCSS(:)=0.
  ZTOT_RRACCSG(:)=0.
  ZTOT_RSACCRG(:)=0.
  ZTOT_RICFRRG(:)=0.
  ZTOT_RRCFRIG(:)=0.
  ZTOT_RICFRR(:)=0.
  ZTOT_RCWETG(:)=0.
  ZTOT_RIWETG(:)=0.
  ZTOT_RRWETG(:)=0.
  ZTOT_RSWETG(:)=0.
  ZTOT_RCDRYG(:)=0.
  ZTOT_RIDRYG(:)=0.
  ZTOT_RRDRYG(:)=0.
  ZTOT_RSDRYG(:)=0.
  ZTOT_RWETGH(:)=0.
  ZTOT_RGMLTR(:)=0.
  ZTOT_RCWETH(:)=0.
  ZTOT_RIWETH(:)=0.
  ZTOT_RSWETH(:)=0.
  ZTOT_RGWETH(:)=0.
  ZTOT_RRWETH(:)=0.
  ZTOT_RCDRYH(:)=0.
  ZTOT_RIDRYH(:)=0.
  ZTOT_RSDRYH(:)=0.
  ZTOT_RRDRYH(:)=0.
  ZTOT_RGDRYH(:)=0.
  ZTOT_RDRYHG(:)=0.
ENDIF

!IF (LHOOK) CALL DR_HOOK('RAIN_ICE:ANTE_MICRO', 1, ZHOOK_HANDLE1)
!-------------------------------------------------------------------------------
!  optimization by looking for locations where
!  the microphysical fields are larger than a minimal value only !!!
!
IF (KSIZE /= COUNT(LDMICRO)) CALL ABOR1('RAIN_ICE : KSIZE /= COUNT(LDMICRO)')

IF (KSIZE > 0) THEN

  !Maximum number of iterations
  !We only count real iterations (those for which we *compute* tendencies)
  INB_ITER_MAX=NMAXITER
  IF(XTSTEP_TS/=0.)THEN
    INB_ITER_MAX=MAX(1, INT(PTSTEP/XTSTEP_TS)) !At least the number of iterations needed for the time-splitting
    ZTSTEP=PTSTEP/INB_ITER_MAX
    INB_ITER_MAX=MAX(NMAXITER, INB_ITER_MAX) !Fot the case XMRSTEP/=0. at the same time
  ENDIF

!===============================================================================================================
! Cache-blocking loop :

  LLSIGMA_RC=(HSUBG_AUCV_RC=='PDF ' .AND. CSUBG_PR_PDF=='SIGM')

  ! starting indexes :
  IC=0
  ISTK=1
  ISTJ=1
  ISTI=1

  DO JMICRO=1,KSIZE,KPROMA

    IMICRO=MIN(KPROMA,KSIZE-JMICRO+1)
!
!*       3.     PACKING
!               --------
    !IF (LHOOK) CALL DR_HOOK('RAIN_ICE:PACK', 0, ZHOOK_HANDLE4)

    ! Setup packing parameters
    OUTER_LOOP: DO JK = ISTK, KKT
      DO JJ = ISTJ, KJT
        IF (ANY(LDMICRO(:,JJ,JK))) THEN
          DO JI = ISTI, KIT
            IF (LDMICRO(JI,JJ,JK)) THEN
              IC=IC+1
              ! Initialization of variables in packed format :
              ZVART(IC,ITH)=PTHT(JI,JJ,JK)
              ZVART(IC,IRV)=PRVT(JI,JJ,JK)
              ZVART(IC,IRC)=PRCT(JI,JJ,JK)
              ZVART(IC,IRR)=PRRT(JI,JJ,JK)
              ZVART(IC,IRI)=PRIT(JI,JJ,JK)
              ZVART(IC,IRS)=PRST(JI,JJ,JK)
              ZVART(IC,IRG)=PRGT(JI,JJ,JK)
              IF (KRR==7) THEN
                ZVART(IC,IRH)=PRHT(JI,JJ,JK)
              ENDIF
              IF (LEXT_TEND) THEN
                ZEXTPK(IC,ITH)=PTHS(JI,JJ,JK)
                ZEXTPK(IC,IRV)=PRVS(JI,JJ,JK)
                ZEXTPK(IC,IRC)=PRCS(JI,JJ,JK)
                ZEXTPK(IC,IRR)=PRRS(JI,JJ,JK)
                ZEXTPK(IC,IRI)=PRIS(JI,JJ,JK)
                ZEXTPK(IC,IRS)=PRSS(JI,JJ,JK)
                ZEXTPK(IC,IRG)=PRGS(JI,JJ,JK)
                IF (KRR==7) THEN
                  ZEXTPK(IC,IRH)=PRHS(JI,JJ,JK)
                ENDIF
                !The th tendency is not related to a mixing ratio change, there is no exn/exnref issue here
              ENDIF
              ZCIT       (IC)=PCIT    (JI,JJ,JK)
              ZCF        (IC)=PCLDFR  (JI,JJ,JK)
              ZRHODREF   (IC)=PRHODREF(JI,JJ,JK)
              ZPRES      (IC)=PPABST  (JI,JJ,JK)
              ZEXN       (IC)=PEXN    (JI,JJ,JK)
              IF(LLSIGMA_RC) THEN
                ZSIGMA_RC(IC)=PSIGS   (JI,JJ,JK)
              ENDIF
              ! Save indices for later usages:
              I1(IC) = JI
              I2(IC) = JJ
              I3(IC) = JK
              IF (IC==IMICRO) THEN
                ! the end of the chunk has been reached, then reset the starting index :
                ISTI=JI+1
                IF (ISTI <= KIT) THEN
                  ISTJ=JJ
                  ISTK=JK
                ELSE
                  ! end of line, restart from 1 and increment upper loop
                  ISTI=1
                  ISTJ=JJ+1
                  IF (ISTJ <= KJT) THEN
                    ISTK=JK
                  ELSE
                    ! end of line, restart from 1 and increment upper loop
                    ISTJ=1
                    ISTK=JK+1
                    IF (ISTK > KKT) THEN
                      ! end of line, restart from 1
                      ISTK=1
                    ENDIF
                  ENDIF
                ENDIF
                IC=0
                EXIT OUTER_LOOP
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        ! restart inner loop on JI :
        ISTI=1
      ENDDO
      ! restart inner loop on JJ :
      ISTJ=1
    ENDDO OUTER_LOOP

    IF (LEXT_TEND) THEN
      DO JL = 1,IMICRO
        ZEXTPK(JL,ITH)=ZEXTPK(JL,ITH)-ZVART(JL,ITH)*ZINV_TSTEP
        ZEXTPK(JL,IRV)=ZEXTPK(JL,IRV)-ZVART(JL,IRV)*ZINV_TSTEP
        ZEXTPK(JL,IRC)=ZEXTPK(JL,IRC)-ZVART(JL,IRC)*ZINV_TSTEP
        ZEXTPK(JL,IRR)=ZEXTPK(JL,IRR)-ZVART(JL,IRR)*ZINV_TSTEP
        ZEXTPK(JL,IRI)=ZEXTPK(JL,IRI)-ZVART(JL,IRI)*ZINV_TSTEP
        ZEXTPK(JL,IRS)=ZEXTPK(JL,IRS)-ZVART(JL,IRS)*ZINV_TSTEP
        ZEXTPK(JL,IRG)=ZEXTPK(JL,IRG)-ZVART(JL,IRG)*ZINV_TSTEP
        IF (KRR==7) THEN
          ZEXTPK(JL,IRH)=ZEXTPK(JL,IRH)-ZVART(JL,IRH)*ZINV_TSTEP
        ENDIF
      ENDDO 
    ENDIF
    IF (LLSIGMA_RC) THEN
      DO JL = 1,IMICRO
        ZSIGMA_RC(JL)=ZSIGMA_RC(JL)*2.
      ENDDO 
    ENDIF

    !IF (LHOOK) CALL DR_HOOK('RAIN_ICE:PACK', 1, ZHOOK_HANDLE4)

!-------------------------------------------------------------------------------
!
!*       4.     LOOP
!               ----
!
    !IF (LHOOK) CALL DR_HOOK('RAIN_ICE:IMICRO', 0, ZHOOK_HANDLE2)

    IITER(1:IMICRO)=0
    ZTIME(1:IMICRO)=0. ! Current integration time (all points may have a different integration time)

    DO WHILE(ANY(ZTIME(1:IMICRO)<PTSTEP)) ! Loop to *really* compute tendencies
   
      IF(XTSTEP_TS/=0.) THEN
        ! In this case we need to remember the time when tendencies were computed
        ! because when time has evolved more than a limit, we must re-compute tendencies
        ZTIME_LASTCALL(1:IMICRO)=ZTIME(1:IMICRO)
      ENDIF
      DO JL=1, IMICRO
        IF (ZTIME(JL) < PTSTEP) THEN
          ZCOMPUTE(JL)=1. ! Computation (1.) only for points for which integration time has not reached the timestep
          IITER(JL)=IITER(JL)+1
        ELSE
          ZCOMPUTE(JL)=0.
        ENDIF
      ENDDO
      LL_ANY_ITER=ANY(IITER(1:IMICRO) < INB_ITER_MAX)
      LLCPZ0T(:)=XMRSTEP/=0. .AND. LL_ANY_ITER
      LSOFT=.FALSE. ! We *really* compute the tendencies
  
      DO WHILE(ANY(ZCOMPUTE(1:IMICRO)==1.)) ! Loop to adjust tendencies when we cross the 0°C or when a species disappears

!$OMP SIMD
        DO JL=1, IMICRO
          ZSUM2(JL)=SUM(ZVART(JL,IRI:KRR))
        ENDDO
        DO JL=1, IMICRO
          ZDEVIDE=(XCPD + XCPV*ZVART(JL,IRV) + XCL*(ZVART(JL,IRC)+ZVART(JL,IRR)) + XCI*ZSUM2(JL)) * ZEXN(JL)
          ZZT(JL) = ZVART(JL,ITH) * ZEXN(JL)
          ZLSFACT(JL)=(XLSTT+(XCPV-XCI)*(ZZT(JL)-XTT)) / ZDEVIDE
          ZLVFACT(JL)=(XLVTT+(XCPV-XCL)*(ZZT(JL)-XTT)) / ZDEVIDE
        ENDDO
        !
        !***       4.1 Tendencies computation
        !
        ! Tendencies are *really* computed when LSOFT==.FALSE. and only adjusted otherwise
        CALL ICE4_TENDENCIES(KPROMA,IMICRO, IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, KKT, KKL, &
                            &KRR, LSOFT, ZCOMPUTE, &
                            &OWARM, CSUBG_RC_RR_ACCR, CSUBG_RR_EVAP, HSUBG_AUCV_RC, CSUBG_PR_PDF, &
                            &ZEXN, ZRHODREF, ZLVFACT, ZLSFACT, I1, I2, I3, &
                            &ZPRES, ZCF, ZSIGMA_RC, &
                            &ZCIT, &
                            &ZZT, ZVART, &
                            &PRRT, &
                            &ZRVHENI_MR, ZRRHONG_MR, ZRIMLTC_MR, ZRSRIMCG_MR, &
                            &ZRCHONI, ZRVDEPS, ZRIAGGS, ZRIAUTS, ZRVDEPG, &
                            &ZRCAUTR, ZRCACCR, ZRREVAV, &
                            &ZRCRIMSS, ZRCRIMSG, ZRSRIMCG, ZRRACCSS, ZRRACCSG, ZRSACCRG, ZRSMLTG, ZRCMLTSR, &
                            &ZRICFRRG, ZRRCFRIG, ZRICFRR, ZRCWETG, ZRIWETG, ZRRWETG, ZRSWETG, &
                            &ZRCDRYG, ZRIDRYG, ZRRDRYG, ZRSDRYG, ZRWETGH, ZRWETGH_MR, ZRGMLTR, &
                            &ZRCWETH, ZRIWETH, ZRSWETH, ZRGWETH, ZRRWETH, &
                            &ZRCDRYH, ZRIDRYH, ZRSDRYH, ZRRDRYH, ZRGDRYH, ZRDRYHG, ZRHMLTR, &
                            &ZRCBERI, &
                            &ZRS_TEND, ZRG_TEND, ZRH_TEND, ZSSI, &
                            &ZA,ZB, &
                            &ZHLC_HCF, ZHLC_LCF, ZHLC_HRC, ZHLC_LRC)

        ! External tendencies
        IF(LEXT_TEND) THEN
          DO JV=0,KRR
            DO JL=1, IMICRO
              ZA(JL,JV) = ZA(JL,JV) + ZEXTPK(JL,JV)
            ENDDO
          ENDDO
        ENDIF
        !
        !***       4.2 Integration time
        !
        ! If we can, we shall use these tendencies until the end of the timestep
        DO JL=1, IMICRO
          ZMAXTIME(JL)=ZCOMPUTE(JL) * (PTSTEP-ZTIME(JL)) ! Remaining time until the end of the timestep
        ENDDO

        !We need to adjust tendencies when temperature reaches 0
        IF(LFEEDBACKT) THEN
          DO JL=1, IMICRO
            !Is ZB(:,ITH) enough to change temperature sign?
            ZX=XTT/ZEXN(JL)
            IF ( (ZVART(JL,ITH) - ZX)*(ZVART(JL,ITH) + ZB(JL,ITH) - ZX) < 0.) THEN
              ZMAXTIME(JL)=0.
            ENDIF 
            !Can ZA(:,ITH) make temperature change of sign?
            IF (ABS(ZA(JL,ITH)) > 1.E-20 ) THEN
              ZTIME_THRESHOLD= (ZX - ZB(JL,ITH) - ZVART(JL,ITH))/ZA(JL,ITH)
              IF (ZTIME_THRESHOLD > 0.) THEN
                ZMAXTIME(JL)=MIN(ZMAXTIME(JL), ZTIME_THRESHOLD)
              ENDIF
            ENDIF
          ENDDO
        ENDIF

        !We need to adjust tendencies when a specy disappears
        !When a specy is missing, only the external tendencies can be negative (and we must keep track of it)
        DO JV=1,KRR
          DO JL=1, IMICRO
            IF ( ZA(JL,JV) < -1.E-20 .AND. ZVART(JL,JV) > XRTMIN(JV) ) THEN
              ZMAXTIME(JL)=MIN(ZMAXTIME(JL),-(ZB(JL,JV)+ZVART(JL,JV))/ZA(JL,JV))
            ENDIF
          ENDDO
        ENDDO

        !We must recompute tendencies when the end of the sub-timestep is reached
        DO JL=1, IMICRO
          !We stop when the end of the timestep is reached
          IF (ZTIME(JL)+ZMAXTIME(JL) >= PTSTEP) THEN
            ZCOMPUTE(JL)=0.
          ENDIF
        ENDDO
        IF (XTSTEP_TS/=0.) THEN
          DO JL=1, IMICRO
            IF ((IITER(JL) < INB_ITER_MAX).AND.(ZTIME(JL)+ZMAXTIME(JL) > ZTIME_LASTCALL(JL)+ZTSTEP)) THEN
              ZMAXTIME(JL)=ZTIME_LASTCALL(JL)-ZTIME(JL)+ZTSTEP
              ZCOMPUTE(JL)=0.
            ENDIF
          ENDDO
        ENDIF

        !We must recompute tendencies when the maximum allowed change is reached
        !When a specy is missing, only the external tendencies can be active and we do not want to recompute
        !the microphysical tendencies when external tendencies are negative (results won't change because specy was already missing)
        IF (XMRSTEP/=0.) THEN
          IF (LL_ANY_ITER) THEN
            ! In this case we need to remember the mixing ratios used to compute the tendencies
            ! because when mixing ratio has evolved more than a threshold, we must re-compute tendencies
            DO JV=1,KRR
              IF (LLCPZ0T(JV)) Z0T(:)=ZVART(:,JV)
              DO JL=1,IMICRO
                IF (IITER(JL)<INB_ITER_MAX .AND. ABS(ZA(JL,JV))>1.E-20) THEN
                  ZTHRESHOLD(JL)=(XMRSTEP+Z0T(JL)-ZVART(JL,JV)-ZB(JL,JV))/ZA(JL,JV)
                ELSE
                  ZTHRESHOLD(JL)=-1.
                ENDIF
              ENDDO
              DO JL=1,IMICRO
                IF ( ZTHRESHOLD(JL)>=0. .AND. ZTHRESHOLD(JL)<ZMAXTIME(JL) .AND. (ZVART(JL,JV)>XRTMIN(6) .OR. ZA(JL,JV)>0.) ) THEN
                  ZMAXTIME(JL)=MIN(ZMAXTIME(JL),ZTHRESHOLD(JL))
                  ZCOMPUTE(JL)=0.
                ENDIF
              ENDDO
            ENDDO
            LLCPZ0T(:)=.FALSE.
!$OMP SIMD
            DO JL=1,IMICRO
              ZMAXB(JL)=MAXVAL(ABS(ZB(JL,1:KRR)))
            ENDDO
            DO JL=1,IMICRO
              IF (IITER(JL)<INB_ITER_MAX .AND. ZMAXB(JL)>XMRSTEP) THEN
                ZMAXTIME(JL)=0.
                ZCOMPUTE(JL)=0.
              ENDIF
            ENDDO
          ENDIF ! LL_ANY_ITER
        ENDIF ! XMRSTEP/=0.
        !
        !***       4.3 New values of variables for next iteration
        !
        DO JV=0,KRR 
          DO JL=1, IMICRO
            ZVART(JL,JV)=ZVART(JL,JV)+ZA(JL,JV)*ZMAXTIME(JL)+ZB(JL,JV)
          ENDDO
        ENDDO
        DO JL=1, IMICRO
          IF (ZVART(JL,IRI)==0.) ZCIT(JL) = 0.
          ZTIME(JL)=ZTIME(JL)+ZMAXTIME(JL)
        ENDDO

        !
        !***       4.4 Mixing ratio change due to each process
        !
        IF(LBU_ENABLE) THEN
          DO JL=1, IMICRO
            ZTOT_RVHENI(JMICRO+JL-1)=ZTOT_RVHENI(JMICRO+JL-1)+ZRVHENI_MR(JL)
            ZTOT_RCHONI(JMICRO+JL-1)=ZTOT_RCHONI(JMICRO+JL-1)+ZRCHONI(JL)*ZMAXTIME(JL)
            ZTOT_RRHONG(JMICRO+JL-1)=ZTOT_RRHONG(JMICRO+JL-1)+ZRRHONG_MR(JL)
            ZTOT_RVDEPS(JMICRO+JL-1)=ZTOT_RVDEPS(JMICRO+JL-1)+ZRVDEPS(JL)*ZMAXTIME(JL)
            ZTOT_RIAGGS(JMICRO+JL-1)=ZTOT_RIAGGS(JMICRO+JL-1)+ZRIAGGS(JL)*ZMAXTIME(JL)
            ZTOT_RIAUTS(JMICRO+JL-1)=ZTOT_RIAUTS(JMICRO+JL-1)+ZRIAUTS(JL)*ZMAXTIME(JL)
            ZTOT_RVDEPG(JMICRO+JL-1)=ZTOT_RVDEPG(JMICRO+JL-1)+ZRVDEPG(JL)*ZMAXTIME(JL)
            ZTOT_RCAUTR(JMICRO+JL-1)=ZTOT_RCAUTR(JMICRO+JL-1)+ZRCAUTR(JL)*ZMAXTIME(JL)
            ZTOT_RCACCR(JMICRO+JL-1)=ZTOT_RCACCR(JMICRO+JL-1)+ZRCACCR(JL)*ZMAXTIME(JL)
            ZTOT_RREVAV(JMICRO+JL-1)=ZTOT_RREVAV(JMICRO+JL-1)+ZRREVAV(JL)*ZMAXTIME(JL)
            ZTOT_RCRIMSS(JMICRO+JL-1)=ZTOT_RCRIMSS(JMICRO+JL-1)+ZRCRIMSS(JL)*ZMAXTIME(JL)
            ZTOT_RCRIMSG(JMICRO+JL-1)=ZTOT_RCRIMSG(JMICRO+JL-1)+ZRCRIMSG(JL)*ZMAXTIME(JL)
            ZTOT_RSRIMCG(JMICRO+JL-1)=ZTOT_RSRIMCG(JMICRO+JL-1)+ZRSRIMCG(JL)*ZMAXTIME(JL)+ZRSRIMCG_MR(JL)
            ZTOT_RRACCSS(JMICRO+JL-1)=ZTOT_RRACCSS(JMICRO+JL-1)+ZRRACCSS(JL)*ZMAXTIME(JL)
            ZTOT_RRACCSG(JMICRO+JL-1)=ZTOT_RRACCSG(JMICRO+JL-1)+ZRRACCSG(JL)*ZMAXTIME(JL)
            ZTOT_RSACCRG(JMICRO+JL-1)=ZTOT_RSACCRG(JMICRO+JL-1)+ZRSACCRG(JL)*ZMAXTIME(JL)
            ZTOT_RSMLTG(JMICRO+JL-1)=ZTOT_RSMLTG(JMICRO+JL-1)+ZRSMLTG(JL)*ZMAXTIME(JL)
            ZTOT_RCMLTSR(JMICRO+JL-1)=ZTOT_RCMLTSR(JMICRO+JL-1)+ZRCMLTSR(JL)*ZMAXTIME(JL)
            ZTOT_RICFRRG(JMICRO+JL-1)=ZTOT_RICFRRG(JMICRO+JL-1)+ZRICFRRG(JL)*ZMAXTIME(JL)
            ZTOT_RRCFRIG(JMICRO+JL-1)=ZTOT_RRCFRIG(JMICRO+JL-1)+ZRRCFRIG(JL)*ZMAXTIME(JL)
            ZTOT_RICFRR(JMICRO+JL-1)=ZTOT_RICFRR(JMICRO+JL-1)+ZRICFRR(JL)*ZMAXTIME(JL)
            ZTOT_RCWETG(JMICRO+JL-1)=ZTOT_RCWETG(JMICRO+JL-1)+ZRCWETG(JL)*ZMAXTIME(JL)
            ZTOT_RIWETG(JMICRO+JL-1)=ZTOT_RIWETG(JMICRO+JL-1)+ZRIWETG(JL)*ZMAXTIME(JL)
            ZTOT_RRWETG(JMICRO+JL-1)=ZTOT_RRWETG(JMICRO+JL-1)+ZRRWETG(JL)*ZMAXTIME(JL)
            ZTOT_RSWETG(JMICRO+JL-1)=ZTOT_RSWETG(JMICRO+JL-1)+ZRSWETG(JL)*ZMAXTIME(JL)
            ZTOT_RWETGH(JMICRO+JL-1)=ZTOT_RWETGH(JMICRO+JL-1)+ZRWETGH(JL)*ZMAXTIME(JL)+ZRWETGH_MR(JL)
            ZTOT_RCDRYG(JMICRO+JL-1)=ZTOT_RCDRYG(JMICRO+JL-1)+ZRCDRYG(JL)*ZMAXTIME(JL)
            ZTOT_RIDRYG(JMICRO+JL-1)=ZTOT_RIDRYG(JMICRO+JL-1)+ZRIDRYG(JL)*ZMAXTIME(JL)
            ZTOT_RRDRYG(JMICRO+JL-1)=ZTOT_RRDRYG(JMICRO+JL-1)+ZRRDRYG(JL)*ZMAXTIME(JL)
            ZTOT_RSDRYG(JMICRO+JL-1)=ZTOT_RSDRYG(JMICRO+JL-1)+ZRSDRYG(JL)*ZMAXTIME(JL)
            ZTOT_RGMLTR(JMICRO+JL-1)=ZTOT_RGMLTR(JMICRO+JL-1)+ZRGMLTR(JL)*ZMAXTIME(JL)
            ZTOT_RCWETH(JMICRO+JL-1)=ZTOT_RCWETH(JMICRO+JL-1)+ZRCWETH(JL)*ZMAXTIME(JL)
            ZTOT_RIWETH(JMICRO+JL-1)=ZTOT_RIWETH(JMICRO+JL-1)+ZRIWETH(JL)*ZMAXTIME(JL)
            ZTOT_RSWETH(JMICRO+JL-1)=ZTOT_RSWETH(JMICRO+JL-1)+ZRSWETH(JL)*ZMAXTIME(JL)
            ZTOT_RGWETH(JMICRO+JL-1)=ZTOT_RGWETH(JMICRO+JL-1)+ZRGWETH(JL)*ZMAXTIME(JL)
            ZTOT_RRWETH(JMICRO+JL-1)=ZTOT_RRWETH(JMICRO+JL-1)+ZRRWETH(JL)*ZMAXTIME(JL)
            ZTOT_RCDRYH(JMICRO+JL-1)=ZTOT_RCDRYH(JMICRO+JL-1)+ZRCDRYH(JL)*ZMAXTIME(JL)
            ZTOT_RIDRYH(JMICRO+JL-1)=ZTOT_RIDRYH(JMICRO+JL-1)+ZRIDRYH(JL)*ZMAXTIME(JL)
            ZTOT_RSDRYH(JMICRO+JL-1)=ZTOT_RSDRYH(JMICRO+JL-1)+ZRSDRYH(JL)*ZMAXTIME(JL)
            ZTOT_RRDRYH(JMICRO+JL-1)=ZTOT_RRDRYH(JMICRO+JL-1)+ZRRDRYH(JL)*ZMAXTIME(JL)
            ZTOT_RGDRYH(JMICRO+JL-1)=ZTOT_RGDRYH(JMICRO+JL-1)+ZRGDRYH(JL)*ZMAXTIME(JL)
            ZTOT_RDRYHG(JMICRO+JL-1)=ZTOT_RDRYHG(JMICRO+JL-1)+ZRDRYHG(JL)*ZMAXTIME(JL)
            ZTOT_RHMLTR(JMICRO+JL-1)=ZTOT_RHMLTR(JMICRO+JL-1)+ZRHMLTR(JL)*ZMAXTIME(JL)
            ZTOT_RIMLTC(JMICRO+JL-1)=ZTOT_RIMLTC(JMICRO+JL-1)+ZRIMLTC_MR(JL)
            ZTOT_RCBERI(JMICRO+JL-1)=ZTOT_RCBERI(JMICRO+JL-1)+ZRCBERI(JL)*ZMAXTIME(JL)
          ENDDO
        ENDIF
        !
        !***       4.5 Next loop
        !
        LSOFT=.TRUE. ! We try to adjust tendencies (inner while loop)
      ENDDO
    ENDDO

    IF(LEXT_TEND) THEN
      !Z..T variables contain the external tendency, we substract it
      DO JV=0,KRR
        DO JL=1, IMICRO
          ZVART(JL,JV) = ZVART(JL,JV) - ZEXTPK(JL,JV) * PTSTEP
        ENDDO
      ENDDO
    ENDIF

    !IF (LHOOK) CALL DR_HOOK('RAIN_ICE:IMICRO', 1, ZHOOK_HANDLE2)
!-------------------------------------------------------------------------------
!
!*       5.     UNPACKING DIAGNOSTICS
!               ---------------------

    !IF (LHOOK) CALL DR_HOOK('RAIN_ICE:UNPACK', 0, ZHOOK_HANDLE5)
    !
    DO JL = 1,IMICRO
      PCIT      (I1(JL),I2(JL),I3(JL))=ZCIT    (JL)
      IF(OWARM) THEN
        PEVAP3D(I1(JL),I2(JL),I3(JL))=ZRREVAV(JL)
      ENDIF
      IF (LLHLC) THEN
        ZHLC_HCF3D(I1(JL),I2(JL),I3(JL))=ZHLC_HCF(JL)
        ZHLC_LCF3D(I1(JL),I2(JL),I3(JL))=ZHLC_LCF(JL)
        ZHLC_HRC3D(I1(JL),I2(JL),I3(JL))=ZHLC_HRC(JL)
        ZHLC_LRC3D(I1(JL),I2(JL),I3(JL))=ZHLC_LRC(JL)
      ENDIF
      !ZWR variables will contain the new S variables values
      ZWR(I1(JL),I2(JL),I3(JL),1)=ZVART(JL,1)
      ZWR(I1(JL),I2(JL),I3(JL),2)=ZVART(JL,2)
      ZWR(I1(JL),I2(JL),I3(JL),3)=ZVART(JL,3)
      ZWR(I1(JL),I2(JL),I3(JL),4)=ZVART(JL,4)
      ZWR(I1(JL),I2(JL),I3(JL),5)=ZVART(JL,5)
      ZWR(I1(JL),I2(JL),I3(JL),6)=ZVART(JL,6)
      IF (KRR==7) THEN
        ZWR(I1(JL),I2(JL),I3(JL),7)=ZVART(JL,7)
      ENDIF
    ENDDO
    !
    !IF (LHOOK) CALL DR_HOOK('RAIN_ICE:UNPACK', 1, ZHOOK_HANDLE5)

  ENDDO ! JMICRO
ENDIF ! KSIZE > 0

!==========================================================================================================

!IF (LHOOK) CALL DR_HOOK('RAIN_ICE:POST_MICRO', 0, ZHOOK_HANDLE3)

!*       6.     COMPUTES THE SLOW COLD PROCESS SOURCES OUTSIDE OF LDMICRO POINTS
!               ----------------------------------------------------------------
!
DO JK = 1, KKT
  DO JJ = 1, KJT
    DO JI = 1, KIT
      ZWR(JI,JJ,JK,IRV)=(ZWR(JI,JJ,JK,IRV)-PRVT(JI,JJ,JK))*ZINV_TSTEP
      ZWR(JI,JJ,JK,IRC)=(ZWR(JI,JJ,JK,IRC)-PRCT(JI,JJ,JK))*ZINV_TSTEP
      ZWR(JI,JJ,JK,IRR)=(ZWR(JI,JJ,JK,IRR)-PRRT(JI,JJ,JK))*ZINV_TSTEP
      ZWR(JI,JJ,JK,IRI)=(ZWR(JI,JJ,JK,IRI)-PRIT(JI,JJ,JK))*ZINV_TSTEP
      ZWR(JI,JJ,JK,IRS)=(ZWR(JI,JJ,JK,IRS)-PRST(JI,JJ,JK))*ZINV_TSTEP
      ZWR(JI,JJ,JK,IRG)=(ZWR(JI,JJ,JK,IRG)-PRGT(JI,JJ,JK))*ZINV_TSTEP
      IF(KRR==7) THEN
        ZWR(JI,JJ,JK,IRH)=(ZWR(JI,JJ,JK,IRH)-PRHT(JI,JJ,JK))*ZINV_TSTEP
      ENDIF
    ENDDO
  ENDDO
ENDDO
!
CALL ICE4_NUCLEATION_WRAPPER(KIT, KJT, KKT, .NOT. LDMICRO, &
                             PTHT, PPABST, PRHODREF, PEXN, ZZ_LSFACT/PEXN, ZT, &
                             PRVT, &
                             PCIT, ZZ_RVHENI_MR)
DO JK = 1, KKT
  DO JJ = 1, KJT
!DEC$ IVDEP
    DO JI = 1, KIT
      ZZ_LSFACT(JI,JJ,JK)=ZZ_LSFACT(JI,JJ,JK)/PEXNREF(JI,JJ,JK)
      ZZ_LVFACT(JI,JJ,JK)=ZZ_LVFACT(JI,JJ,JK)/PEXNREF(JI,JJ,JK)
      ZZ_RVHENI = MIN(PRVS(JI,JJ,JK), ZZ_RVHENI_MR(JI,JJ,JK)/PTSTEP)
      PRIS(JI,JJ,JK)=PRIS(JI,JJ,JK)+ZZ_RVHENI
      PRVS(JI,JJ,JK)=PRVS(JI,JJ,JK)-ZZ_RVHENI
      PTHS(JI,JJ,JK)=PTHS(JI,JJ,JK) + ZZ_RVHENI*ZZ_LSFACT(JI,JJ,JK)
!-------------------------------------------------------------------------------
!* 7.     TOTAL TENDENCIES
      ZWR(JI,JJ,JK,ITH) = (ZWR(JI,JJ,JK,IRC)+ZWR(JI,JJ,JK,IRR))*ZZ_LVFACT(JI,JJ,JK) &
                      & + (ZWR(JI,JJ,JK,IRI)+ZWR(JI,JJ,JK,IRS)+ZWR(JI,JJ,JK,IRG)+ZWR(JI,JJ,JK,IRH))*ZZ_LSFACT(JI,JJ,JK)
      ZWR(JI,JJ,JK,ITH) = PTHS(JI,JJ,JK) + ZWR(JI,JJ,JK,ITH)
      !We apply these tendencies to the S variables
      ZWR(JI,JJ,JK,IRV) = PRVS(JI,JJ,JK) + ZWR(JI,JJ,JK,IRV)
      ZWR(JI,JJ,JK,IRC) = PRCS(JI,JJ,JK) + ZWR(JI,JJ,JK,IRC)
      ZWR(JI,JJ,JK,IRR) = PRRS(JI,JJ,JK) + ZWR(JI,JJ,JK,IRR)
      ZWR(JI,JJ,JK,IRI) = PRIS(JI,JJ,JK) + ZWR(JI,JJ,JK,IRI)
      ZWR(JI,JJ,JK,IRS) = PRSS(JI,JJ,JK) + ZWR(JI,JJ,JK,IRS)
      ZWR(JI,JJ,JK,IRG) = PRGS(JI,JJ,JK) + ZWR(JI,JJ,JK,IRG)
      IF (KRR==7) THEN
        ZWR(JI,JJ,JK,IRH) = PRHS(JI,JJ,JK) + ZWR(JI,JJ,JK,IRH)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!We correct negativities with conservation
CALL CORRECT_NEGATIVITIES(KIT, KJT, KKT, KRR, ZWR(:,:,:,IRV), ZWR(:,:,:,IRC), ZWR(:,:,:,IRR), ZWR(:,:,:,IRI), &
                         &ZWR(:,:,:,IRS), ZWR(:,:,:,IRG), ZWR(:,:,:,ITH), ZZ_LVFACT, ZZ_LSFACT, ZWR(:,:,:,IRH))

!
!***     7.2    LBU_ENABLE case
!
IF(LBU_ENABLE) THEN
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RVHENI(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) + ZW(:,:,:)
  PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*ZZ_LSFACT(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'HENU_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RV) CALL BUDGET(PRVS(:,:,:)*PRHODJ(:,:,:), 6, 'HENU_BU_RRV',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'HENU_BU_RRI',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCHONI(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) + ZW(:,:,:)
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'HON_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'HON_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'HON_BU_RRI',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRHONG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'SFR_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'SFR_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'SFR_BU_RRG',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RVDEPS(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) + ZW(:,:,:)
  PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*ZZ_LSFACT(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'DEPS_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RV) CALL BUDGET(PRVS(:,:,:)*PRHODJ(:,:,:), 6, 'DEPS_BU_RRV',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'DEPS_BU_RRS',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RIAGGS(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) + ZW(:,:,:)
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'AGGS_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'AGGS_BU_RRS',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RIAUTS(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) + ZW(:,:,:)
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'AUTS_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'AUTS_BU_RRS',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RVDEPG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*ZZ_LSFACT(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'DEPG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RV) CALL BUDGET(PRVS(:,:,:)*PRHODJ(:,:,:), 6, 'DEPG_BU_RRV',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'DEPG_BU_RRG',YDDDH, YDLDDH, YDMDDH)

  IF(OWARM) THEN
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RCAUTR(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
    PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'AUTO_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'AUTO_BU_RRR',YDDDH, YDLDDH, YDMDDH)

    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RCACCR(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
    PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'ACCR_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'ACCR_BU_RRR',YDDDH, YDLDDH, YDMDDH)

    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RREVAV(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
    PRVS(:,:,:) = PRVS(:,:,:) + ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) - ZW(:,:,:)*ZZ_LVFACT(:,:,:)
    IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'REVA_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RV) CALL BUDGET(PRVS(:,:,:)*PRHODJ(:,:,:), 6, 'REVA_BU_RRV',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'REVA_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  ENDIF

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCRIMSS(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRSS(:,:,:) = PRSS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCRIMSG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RSRIMCG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'RIM_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'RIM_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'RIM_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'RIM_BU_RRG',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRACCSS(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PRSS(:,:,:) = PRSS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRACCSG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RSACCRG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'ACC_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'ACC_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'ACC_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'ACC_BU_RRG',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RSMLTG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCMLTSR, MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'CMEL_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'CMEL_BU_RRG',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'CMEL_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'CMEL_BU_RRR',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RICFRRG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRCFRIG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RICFRR(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'CFRZ_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'CFRZ_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'CFRZ_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'CFRZ_BU_RRG',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCWETG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRWETG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RIWETG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RSWETG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'WETG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'WETG_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'WETG_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'WETG_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'WETG_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'WETG_BU_RRG',YDDDH, YDLDDH, YDMDDH)

  IF(KRR==7) THEN
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RWETGH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRGS(:,:,:) = PRGS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'GHCV_BU_RRG',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'GHCV_BU_RRH',YDDDH, YDLDDH, YDMDDH)
  END IF

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCDRYG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRDRYG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RIDRYG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RSDRYG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'DRYG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'DRYG_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'DRYG_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'DRYG_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'DRYG_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'DRYG_BU_RRG',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RGMLTR(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) - ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'GMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'GMLT_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'GMLT_BU_RRG',YDDDH, YDLDDH, YDMDDH)

  IF(KRR==7) THEN
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RCWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RRWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RIWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RSWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RGWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRGS(:,:,:) = PRGS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'WETH_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'WETH_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'WETH_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'WETH_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'WETH_BU_RRS',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'WETH_BU_RRH',YDDDH, YDLDDH, YDMDDH)

    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RGWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRGS(:,:,:) = PRGS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'HGCV_BU_RRG',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'HGCV_BU_RRH',YDDDH, YDLDDH, YDMDDH)

    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RCDRYH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RRDRYH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RIDRYH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RSDRYH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RGDRYH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRGS(:,:,:) = PRGS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RDRYHG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRHS(:,:,:) = PRHS(:,:,:) - ZW(:,:,:)
    PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'DRYH_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'DRYH_BU_RRC',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'DRYH_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'DRYH_BU_RRI',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'DRYH_BU_RRS',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'DRYH_BU_RRG',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'DRYH_BU_RRH',YDDDH, YDLDDH, YDMDDH)

    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RHMLTR(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) - ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) - ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
    IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'HMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'HMLT_BU_RRR',YDDDH, YDLDDH, YDMDDH)
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'HMLT_BU_RRH',YDDDH, YDLDDH, YDMDDH)
  ENDIF

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RIMLTC(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  PRCS(:,:,:) = PRCS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) - ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'IMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'IMLT_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'IMLT_BU_RRI',YDDDH, YDLDDH, YDMDDH)

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCBERI(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRIS(:,:,:) = PRIS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'BERFI_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'BERFI_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'BERFI_BU_RRI',YDDDH, YDLDDH, YDMDDH)
ENDIF
!

!***     7.3    Final tendencies
!
DO JK = 1, KKT
  PTHS(:,:,JK) = ZWR(:,:,JK,ITH)
  PRVS(:,:,JK) = ZWR(:,:,JK,IRV)
  PRCS(:,:,JK) = ZWR(:,:,JK,IRC)
  PRRS(:,:,JK) = ZWR(:,:,JK,IRR)
  PRIS(:,:,JK) = ZWR(:,:,JK,IRI)
  PRSS(:,:,JK) = ZWR(:,:,JK,IRS)
  PRGS(:,:,JK) = ZWR(:,:,JK,IRG)
  IF (KRR==7) THEN
    PRHS(:,:,JK) = ZWR(:,:,JK,IRH)
  ENDIF
ENDDO

IF(LBU_ENABLE) THEN
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'CORR_BU_RTH',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RV) CALL BUDGET(PRVS(:,:,:)*PRHODJ(:,:,:), 6, 'CORR_BU_RRV',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'CORR_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'CORR_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'CORR_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'CORR_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'CORR_BU_RRG',YDDDH, YDLDDH, YDMDDH)
  IF (KRR==7) THEN
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'CORR_BU_RRH',YDDDH, YDLDDH, YDMDDH)
  ENDIF
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       8.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!               -------------------------------------
!
IF(LSEDIM_AFTER) THEN
  !
  !*       8.1     sedimentation
  !
  IF(HSEDIM=='STAT') THEN
    IF (KRR==7) THEN
      DO JK = 1, KKT
        DO JJ = 1, KJT
          DO JI = 1, KIT
            ZRCT(JI,JJ,JK)=PRCS(JI,JJ,JK)*PTSTEP
            ZRRT(JI,JJ,JK)=PRRS(JI,JJ,JK)*PTSTEP
            ZRIT(JI,JJ,JK)=PRIS(JI,JJ,JK)*PTSTEP
            ZRST(JI,JJ,JK)=PRSS(JI,JJ,JK)*PTSTEP
            ZRGT(JI,JJ,JK)=PRGS(JI,JJ,JK)*PTSTEP
            ZRHT(JI,JJ,JK)=PRHS(JI,JJ,JK)*PTSTEP
          ENDDO
        ENDDO
      ENDDO
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, &
                                  &PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, ZRCT, PRRS, ZRRT, PRIS, ZRIT,&
                                  &PRSS, ZRST, PRGS, ZRGT,&
                                  &PINPRC, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA=PSEA, PTOWN=PTOWN, &
                                  &PINPRH=PINPRH, PRHT=ZRHT, PRHS=PRHS, PFPR=PFPR)
    ELSE
      DO JK = 1, KKT
        DO JJ = 1, KJT
          DO JI = 1, KIT
            ZRCT(JI,JJ,JK)=PRCS(JI,JJ,JK)*PTSTEP
            ZRRT(JI,JJ,JK)=PRRS(JI,JJ,JK)*PTSTEP
            ZRIT(JI,JJ,JK)=PRIS(JI,JJ,JK)*PTSTEP
            ZRST(JI,JJ,JK)=PRSS(JI,JJ,JK)*PTSTEP
            ZRGT(JI,JJ,JK)=PRGS(JI,JJ,JK)*PTSTEP
          ENDDO
        ENDDO
      ENDDO
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, &
                                  &PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, ZRCT, PRRS, ZRRT, PRIS, ZRIT,&
                                  &PRSS, ZRST, PRGS, ZRGT,&
                                  &PINPRC, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA=PSEA, PTOWN=PTOWN, &
                                  &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !No negativity correction here as we apply sedimentation on ZR.T variables
  ELSEIF(HSEDIM=='SPLI') THEN
    !SR: It *seems* that we must have two separate calls for ifort
    IF(KRR==7) THEN
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, &
                                   &PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA=PSEA, PTOWN=PTOWN, &
                                   &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR)
    ELSE
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, &
                                   &PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA=PSEA, PTOWN=PTOWN, &
                                   &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !We correct negativities with conservation
    !SPLI algorith uses a time-splitting. Inside the loop a temporary m.r. is used.
    !   It is initialized with the m.r. at T and is modified by two tendencies:
    !   sedimentation tendency and an external tendency which represents all other
    !   processes (mainly advection and microphysical processes). If both tendencies
    !   are negative, sedimentation can remove a species at a given sub-timestep. From
    !   this point sedimentation stops for the remaining sub-timesteps but the other tendency
    !   will be still active and will lead to negative values.
    !   We could prevent the algorithm to not consume too much a species, instead we apply
    !   a correction here.
    CALL CORRECT_NEGATIVITIES(KIT, KJT, KKT, KRR, PRVS, PRCS, PRRS, &
                             &PRIS, PRSS, PRGS, &
                             &PTHS, ZZ_LVFACT, ZZ_LSFACT, PRHS)
  ELSE
    WRITE(*,*) ' STOP'
    WRITE(*,*) ' NO SEDIMENTATION SCHEME FOR HSEDIM=', HSEDIM
    CALL ABORT
    STOP
  END IF
  !
  !*       8.2     budget storage
  !
  IF (LBUDGET_RC .AND. OSEDIC) &
                  CALL BUDGET (PRCS(:,:,:)*PRHODJ(:,:,:), 7 , 'SEDI_BU_RRC',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:)*PRHODJ(:,:,:), 8 , 'SEDI_BU_RRR',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RI) CALL BUDGET (PRIS(:,:,:)*PRHODJ(:,:,:), 9 , 'SEDI_BU_RRI',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RS) CALL BUDGET (PRSS(:,:,:)*PRHODJ(:,:,:), 10, 'SEDI_BU_RRS',YDDDH, YDLDDH, YDMDDH)
  IF (LBUDGET_RG) CALL BUDGET (PRGS(:,:,:)*PRHODJ(:,:,:), 11, 'SEDI_BU_RRG',YDDDH, YDLDDH, YDMDDH)
  IF ( KRR == 7 .AND. LBUDGET_RH) &
                  CALL BUDGET (PRHS(:,:,:)*PRHODJ(:,:,:), 12, 'SEDI_BU_RRH',YDDDH, YDLDDH, YDMDDH)
  !
  !sedimentation of rain fraction
  IF (LLRAIN_FRACTION) THEN
    CALL ICE4_RAINFR_VERT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, KKT, KKL, ZRAINFR3D, PRRS(:,:,:)*PTSTEP)
  ENDIF

ENDIF
!
!-------------------------------------------------------------------------------
!
!*       9.     COMPUTE THE FOG DEPOSITION TERM 
!               -------------------------------------
!
IF (LDEPOSC) THEN !cloud water deposition on vegetation
  DO JJ = 1, KJT
!DEC$ IVDEP
    DO JI = 1, KIT
      PRCS(JI,JJ,IKB) = PRCS(JI,JJ,IKB) - XVDEPOSC * PRCT(JI,JJ,IKB) / PDZZ(JI,JJ,IKB)
      PINPRC(JI,JJ) = PINPRC(JI,JJ) + XVDEPOSC * PRCT(JI,JJ,IKB) * PRHODREF(JI,JJ,IKB)/XRHOLW
    ENDDO
  ENDDO

 IF (LBUDGET_RC) CALL BUDGET (PRCS(:,:,:)*PRHODJ(:,:,:),7,'DEPO_BU_RRC',YDDDH, YDLDDH, YDMDDH)
ENDIF

!IF (LHOOK) CALL DR_HOOK('RAIN_ICE:POST_MICRO', 1, ZHOOK_HANDLE3)

IF (LHOOK) CALL DR_HOOK('RAIN_ICE', 1, ZHOOK_HANDLE)
!
CONTAINS
  !
  SUBROUTINE CORRECT_NEGATIVITIES(KIT, KJT, KKT, KRR, PRV, PRC, PRR, &
                                 &PRI, PRS, PRG, &
                                 &PTH, PLVFACT, PLSFACT, PRH)
  !
  IMPLICIT NONE
  !
  INTEGER,                INTENT(IN)    :: KIT, KJT, KKT, KRR
  REAL, DIMENSION(KIT, KJT, KKT), INTENT(INOUT) :: PRV, PRC, PRR, PRI, PRS, PRG, PTH
  REAL, DIMENSION(KIT, KJT, KKT), INTENT(IN)    :: PLVFACT, PLSFACT
  REAL, DIMENSION(KIT, KJT, KKT), OPTIONAL, INTENT(INOUT) :: PRH
  !
  REAL :: ZW
  INTEGER :: JI, JJ, JK
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  !
  IF (LHOOK) CALL DR_HOOK('RAIN_ICE:CORRECT_NEGATIVITIES', 0, ZHOOK_HANDLE)
  !
  !We correct negativities with conservation
  ! 1) deal with negative values for mixing ratio, except for vapor
  DO JK = 1, KKT
    DO JJ = 1, KJT
      DO JI = 1, KIT
        ZW =PRC(JI,JJ,JK)-MAX(PRC(JI,JJ,JK), 0.)
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW*PLVFACT(JI,JJ,JK)
        PRC(JI,JJ,JK)=PRC(JI,JJ,JK)-ZW

        ZW =PRR(JI,JJ,JK)-MAX(PRR(JI,JJ,JK), 0.)
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW*PLVFACT(JI,JJ,JK)
        PRR(JI,JJ,JK)=PRR(JI,JJ,JK)-ZW

        ZW =PRI(JI,JJ,JK)-MAX(PRI(JI,JJ,JK), 0.)
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW*PLSFACT(JI,JJ,JK)
        PRI(JI,JJ,JK)=PRI(JI,JJ,JK)-ZW

        ZW =PRS(JI,JJ,JK)-MAX(PRS(JI,JJ,JK), 0.)
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW*PLSFACT(JI,JJ,JK)
        PRS(JI,JJ,JK)=PRS(JI,JJ,JK)-ZW

        ZW =PRG(JI,JJ,JK)-MAX(PRG(JI,JJ,JK), 0.)
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW*PLSFACT(JI,JJ,JK)
        PRG(JI,JJ,JK)=PRG(JI,JJ,JK)-ZW

        IF(KRR==7) THEN
          ZW =PRH(JI,JJ,JK)-MAX(PRH(JI,JJ,JK), 0.)
          PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW
          PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW*PLSFACT(JI,JJ,JK)
          PRH(JI,JJ,JK)=PRH(JI,JJ,JK)-ZW
        ENDIF

  ! 2) deal with negative vapor mixing ratio

        ! for rc and ri, we keep ice fraction constant
        ZW=MIN(1., MAX(XRTMIN(1)-PRV(JI,JJ,JK), 0.) / &
                            &MAX(PRC(JI,JJ,JK)+PRI(JI,JJ,JK), 1.E-20)) ! Proportion of rc+ri to convert into rv
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW* &
                     &(PRC(JI,JJ,JK)*PLVFACT(JI,JJ,JK)+PRI(JI,JJ,JK)*PLSFACT(JI,JJ,JK))
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW*(PRC(JI,JJ,JK)+PRI(JI,JJ,JK))
        PRC(JI,JJ,JK)=(1.-ZW)*PRC(JI,JJ,JK)
        PRI(JI,JJ,JK)=(1.-ZW)*PRI(JI,JJ,JK)

        ZW=MIN(MAX(PRR(JI,JJ,JK), 0.), &
                        &MAX(XRTMIN(1)-PRV(JI,JJ,JK), 0.)) ! Quantity of rr to convert into rv
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW
        PRR(JI,JJ,JK)=PRR(JI,JJ,JK)-ZW
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW*PLVFACT(JI,JJ,JK)

        ZW=MIN(MAX(PRS(JI,JJ,JK), 0.), &
                        &MAX(XRTMIN(1)-PRV(JI,JJ,JK), 0.)) ! Quantity of rs to convert into rv
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW
        PRS(JI,JJ,JK)=PRS(JI,JJ,JK)-ZW
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW*PLSFACT(JI,JJ,JK)

        ZW=MIN(MAX(PRG(JI,JJ,JK), 0.), &
                        &MAX(XRTMIN(1)-PRV(JI,JJ,JK), 0.)) ! Quantity of rg to convert into rv
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW
        PRG(JI,JJ,JK)=PRG(JI,JJ,JK)-ZW
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW*PLSFACT(JI,JJ,JK)

        IF(KRR==7) THEN
          ZW=MIN(MAX(PRH(JI,JJ,JK), 0.), &
                          &MAX(XRTMIN(1)-PRV(JI,JJ,JK), 0.)) ! Quantity of rh to convert into rv
          PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW
          PRH(JI,JJ,JK)=PRH(JI,JJ,JK)-ZW
          PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW*PLSFACT(JI,JJ,JK)
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  !
  IF (LHOOK) CALL DR_HOOK('RAIN_ICE:CORRECT_NEGATIVITIES', 1, ZHOOK_HANDLE)
  !
  END SUBROUTINE CORRECT_NEGATIVITIES

!
END SUBROUTINE RAIN_ICE
