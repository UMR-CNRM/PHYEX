!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
       MODULE MODI_RAIN_ICE_OLD
!      ####################
!
INTERFACE
      SUBROUTINE RAIN_ICE_OLD (D, OSEDIC,HSEDIM, HSUBG_AUCV, OWARM, KKA, KKU, KKL,      &
                            KSPLITR, PTSTEP, KRR,                            &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
                            PINPRC,PINPRR, PINPRR3D, PEVAP3D,           &
                            PINPRS, PINPRG, PSIGS, PINDEP, PRAINFR, PSEA, PTOWN,  &
                            PRHT, PRHS, PINPRH, PFPR                        )
!
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV ! Switch for rc->rr Subgrid autoconversion
                                        ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
!
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step 
                                      ! integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ     ! Layer thikness (m)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRC! Cloud instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRR! Rain instant precip
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEVAP3D! Rain evap profile
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG! Graupel instant precip
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRAINFR! Rain fraction
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(IN)    :: PSEA ! Sea Mask
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(IN)    :: PTOWN! Fraction that is town
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(INOUT) :: PINPRH! Hail instant precip
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)   :: PFPR ! upper-air precipitation fluxes
!
END SUBROUTINE RAIN_ICE_OLD
END INTERFACE
END MODULE MODI_RAIN_ICE_OLD
!     ######spl
      SUBROUTINE RAIN_ICE_OLD (D, OSEDIC,HSEDIM, HSUBG_AUCV, OWARM, KKA, KKU, KKL,      &
                            KSPLITR, PTSTEP, KRR,                            &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
                            PINPRC,PINPRR, PINPRR3D, PEVAP3D,           &
                            PINPRS, PINPRG, PSIGS, PINDEP, PRAINFR, PSEA, PTOWN,  &
                            PRHT, PRHS, PINPRH, PFPR                        )
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
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1 and Book2 of documentation ( routine RAIN_ICE_OLD )
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
!!      J.Escobar : 8/2018 : for real*4 , bis => limit exp() in RAIN_ICE_SLOW with XMNH_HUGE_12_LOG
!!      P.Wautelet 01/02/2019: add missing initialization for PFPR
!!                   02/2019 C.Lac add rain fraction as an output field
!  P. Wautelet 25/02/2019: split rain_ice (cleaner and easier to maintain/debug)
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet 29/05/2019: remove PACK/UNPACK intrinsics (to get more performance and better OpenACC support)
!  J. Escobar  09/07/2019: for reproductiblity MPPDB_CHECK, add missing LCHECK test in ZRHODJ de/allocate
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets (no more budget calls in this subroutine)
!-----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,         only: lbu_enable
use MODD_CONF,           only: LCHECK
use MODD_CST,            only: XCI, XCL, XCPD, XCPV, XLSTT, XLVTT, XTT, &
                               XALPI, XBETAI, XGAMI, XMD, XMV, XTT
use MODD_LES,            only: LLES_CALL
use MODD_PARAMETERS,     only: JPVEXT
use MODD_PARAM_ICE,      only: CSUBG_PR_PDF, LDEPOSC
use MODD_RAIN_ICE_DESCR, only: RAIN_ICE_DESCR, XLBEXR, XLBR, XRTMIN
use MODD_RAIN_ICE_PARAM, only: XCRIAUTC
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t

use MODE_MSG
use MODE_RAIN_ICE_FAST_RG,             only: RAIN_ICE_FAST_RG
use MODE_RAIN_ICE_FAST_RH,             only: RAIN_ICE_FAST_RH
use MODE_RAIN_ICE_FAST_RI,             only: RAIN_ICE_FAST_RI
use MODE_RAIN_ICE_FAST_RS,             only: RAIN_ICE_FAST_RS
use MODE_RAIN_ICE_NUCLEATION,          only: RAIN_ICE_NUCLEATION
use MODE_RAIN_ICE_SEDIMENTATION_SPLIT, only: RAIN_ICE_SEDIMENTATION_SPLIT
use MODE_RAIN_ICE_SEDIMENTATION_STAT,  only: RAIN_ICE_SEDIMENTATION_STAT
use MODE_RAIN_ICE_SLOW,                only: RAIN_ICE_SLOW
use MODE_RAIN_ICE_WARM,                only: RAIN_ICE_WARM
use mode_tools,                        only: Countjv
use mode_tools_ll,                     only: GET_INDICE_ll

USE MODE_ICE4_RAINFR_VERT
!

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV ! Switch for rc->rr Subgrid autoconversion
                                        ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
!
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step
                                      ! integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ     ! Layer thikness (m)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRC! Cloud instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRR! Rain instant precip
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEVAP3D! Rain evap profile
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG! Graupel instant precip
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PRAINFR! Rain fraction
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(IN)    :: PSEA ! Sea Mask
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(IN)    :: PTOWN! Fraction that is town
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(INOUT) :: PINPRH! Hail instant precip
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)   :: PFPR ! upper-air precipitation fluxes
!
!*       0.2   Declarations of local variables :
!
INTEGER                           :: IIB           !  Define the domain where is
INTEGER                           :: IIE           !  the microphysical sources have to be computed
INTEGER                           :: IIT           !
INTEGER                           :: IJB           !
INTEGER                           :: IJE           !
INTEGER                           :: IJT           !
INTEGER                           :: IKB,IKTB,IKT  !
INTEGER                           :: IKE,IKTE      !
!
INTEGER                           :: IMICRO
INTEGER, DIMENSION(SIZE(PEXNREF)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
LOGICAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) &
                                  :: GMICRO ! Test where to compute all processes
REAL                              :: ZINVTSTEP
REAL                              :: ZCOEFFRCM
REAL, DIMENSION(:), ALLOCATABLE   :: ZRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRRT    ! Rain water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRGT    ! Graupel m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZRHT    ! Hail m.r. at t
REAL, DIMENSION(:), ALLOCATABLE   :: ZCIT    ! Pristine ice conc. at t
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZRVS    ! Water vapor m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRRS    ! Rain water m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRGS    ! Graupel m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZRHS    ! Hail m.r. source
REAL, DIMENSION(:), ALLOCATABLE   :: ZTHS    ! Theta source
REAL, DIMENSION(:), ALLOCATABLE   :: ZTHT    ! Potential temperature
REAL, DIMENSION(:), ALLOCATABLE   :: ZTHLT   ! Liquid potential temperature
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZRHODREF, &      ! RHO Dry REFerence
                                     ZRHODJ,   &      ! RHO times Jacobian
                                     ZZT,      &      ! Temperature
                                     ZPRES,    &      ! Pressure
                                     ZEXNREF,  &      ! EXNer Pressure REFerence
                                     ZZW,      &      ! Work array
                                     ZLSFACT,  &      ! L_s/(Pi_ref*C_ph)
                                     ZLVFACT,  &      ! L_v/(Pi_ref*C_ph)
                                     ZUSW,     &      ! Undersaturation over water
                                     ZSSI,     &      ! Supersaturation over ice
                                     ZLBDAR,   &      ! Slope parameter of the raindrop  distribution
                                     ZLBDAR_RF,&      ! Slope parameter of the raindrop  distribution
                                                      ! for the Rain Fraction part
                                     ZLBDAS,   &      ! Slope parameter of the aggregate distribution
                                     ZLBDAG,   &      ! Slope parameter of the graupel   distribution
                                     ZLBDAH,   &      ! Slope parameter of the hail      distribution
                                     ZRDRYG,   &      ! Dry growth rate of the graupeln
                                     ZRWETG,   &      ! Wet growth rate of the graupeln
                                     ZAI,      &      ! Thermodynamical function
                                     ZCJ,      &      ! Function to compute the ventilation coefficient
                                     ZKA,      &      ! Thermal conductivity of the air
                                     ZDV,      &      ! Diffusivity of water vapor in the air
                                     ZSIGMA_RC,&      ! Standard deviation of rc at time t
                                     ZCF,      &      ! Cloud fraction
                                     ZRF,      &      ! Rain fraction
                                     ZHLC_HCF, &      ! HLCLOUDS : fraction of High Cloud Fraction in grid
                                     ZHLC_LCF, &      ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                                      !    note that ZCF = ZHLC_HCF + ZHLC_LCF
                                     ZHLC_HRC, &      ! HLCLOUDS : LWC that is High LWC in grid
                                     ZHLC_LRC, &      ! HLCLOUDS : LWC that is Low  LWC in grid
                                                      !    note that ZRC = ZHLC_HRC + ZHLC_LRC
                                     ZHLC_RCMAX, &    ! HLCLOUDS : maximum value for RC in distribution
                                     ZRCRAUTC, &      ! RC value to begin rain formation =XCRIAUTC/RHODREF
                                     ZHLC_HRCLOCAL, & ! HLCLOUDS : LWC that is High LWC local in HCF
                                     ZHLC_LRCLOCAL    ! HLCLOUDS : LWC that is Low  LWC local in LCF
                                                      !    note that ZRC/CF = ZHLC_HRCLOCAL+ ZHLC_LRCLOCAL
                                                      !                     = ZHLC_HRC/HCF+ ZHLC_LRC/LCF
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZW1 ! Work arrays
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
                                  :: ZW ! work array
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
                                  :: ZT ! Temperature
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!               -----------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IIT=SIZE(PDZZ,1)
IJT=SIZE(PDZZ,2)
IKB=KKA+JPVEXT*KKL
IKE=KKU-JPVEXT*KKL
IKT=SIZE(PDZZ,3)
IKTB=1+JPVEXT
IKTE=IKT-JPVEXT
!
!
ZINVTSTEP=1./PTSTEP
!
!
!*       2.     COMPUTES THE SLOW COLD PROCESS SOURCES
!               --------------------------------------
!
CALL RAIN_ICE_NUCLEATION(IIB, IIE, IJB, IJE, IKTB, IKTE,KRR,PTSTEP,&
     PTHT,PPABST,PRHODJ,PRHODREF,PRVT,PRCT,PRRT,PRIT,PRST,PRGT,&
     PCIT,PEXNREF,PTHS,PRVS,PRIS,ZT,PRHT)
!
!
!  optimization by looking for locations where
!  the microphysical fields are larger than a minimal value only !!!
!
GMICRO(:,:,:) = .FALSE.

 IF ( KRR == 7 ) THEN
  GMICRO(IIB:IIE,IJB:IJE,:) =                          &
                PRCT(IIB:IIE,IJB:IJE,:)>XRTMIN(2) .OR. &
                PRRT(IIB:IIE,IJB:IJE,:)>XRTMIN(3) .OR. &
                PRIT(IIB:IIE,IJB:IJE,:)>XRTMIN(4) .OR. &
                PRST(IIB:IIE,IJB:IJE,:)>XRTMIN(5) .OR. &
                PRGT(IIB:IIE,IJB:IJE,:)>XRTMIN(6) .OR. &
                PRHT(IIB:IIE,IJB:IJE,:)>XRTMIN(7)
 ELSE IF( KRR == 6 ) THEN
  GMICRO(IIB:IIE,IJB:IJE,:) =                          &
                PRCT(IIB:IIE,IJB:IJE,:)>XRTMIN(2) .OR. &
                PRRT(IIB:IIE,IJB:IJE,:)>XRTMIN(3) .OR. &
                PRIT(IIB:IIE,IJB:IJE,:)>XRTMIN(4) .OR. &
                PRST(IIB:IIE,IJB:IJE,:)>XRTMIN(5) .OR. &
                PRGT(IIB:IIE,IJB:IJE,:)>XRTMIN(6)
 END IF

IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
IF( IMICRO >= 0 ) THEN
  ALLOCATE(ZRVT(IMICRO))
  ALLOCATE(ZRCT(IMICRO))
  ALLOCATE(ZRRT(IMICRO))
  ALLOCATE(ZRIT(IMICRO))
  ALLOCATE(ZRST(IMICRO))
  ALLOCATE(ZRGT(IMICRO))
  IF ( KRR == 7 ) THEN
    ALLOCATE(ZRHT(IMICRO))
  ELSE
    ALLOCATE(ZRHT(0))
  END IF
  ALLOCATE(ZCIT(IMICRO))
  ALLOCATE(ZRVS(IMICRO))
  ALLOCATE(ZRCS(IMICRO))
  ALLOCATE(ZRRS(IMICRO))
  ALLOCATE(ZRIS(IMICRO))
  ALLOCATE(ZRSS(IMICRO))
  ALLOCATE(ZRGS(IMICRO))
  IF ( KRR == 7 ) THEN
    ALLOCATE(ZRHS(IMICRO))
  ELSE
    ALLOCATE(ZRHS(0))
  END IF
  ALLOCATE(ZTHS(IMICRO))
  ALLOCATE(ZTHT(IMICRO))
  ALLOCATE(ZTHLT(IMICRO))
  ALLOCATE(ZRHODREF(IMICRO))
  ALLOCATE(ZZT(IMICRO))
  ALLOCATE(ZPRES(IMICRO))
  ALLOCATE(ZEXNREF(IMICRO))
  ALLOCATE(ZSIGMA_RC(IMICRO))
  ALLOCATE(ZCF(IMICRO))
  ALLOCATE(ZRF(IMICRO))
  ALLOCATE(ZHLC_HCF(IMICRO))
  ALLOCATE(ZHLC_LCF(IMICRO))
  ALLOCATE(ZHLC_HRC(IMICRO))
  ALLOCATE(ZHLC_LRC(IMICRO))
  ALLOCATE(ZHLC_RCMAX(IMICRO))
  ALLOCATE(ZRCRAUTC(IMICRO))
  ALLOCATE(ZHLC_HRCLOCAL(IMICRO))
  ALLOCATE(ZHLC_LRCLOCAL(IMICRO))

  DO JL=1,IMICRO
    ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
    ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
    ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
    ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
    ZRST(JL) = PRST(I1(JL),I2(JL),I3(JL))
    ZRGT(JL) = PRGT(I1(JL),I2(JL),I3(JL))
    IF ( KRR == 7 ) ZRHT(JL) = PRHT(I1(JL),I2(JL),I3(JL))
    ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
    ZCF(JL) = PCLDFR(I1(JL),I2(JL),I3(JL))
    IF ( HSUBG_AUCV == 'PDF ' .AND. CSUBG_PR_PDF == 'SIGM' ) THEN
      ZSIGMA_RC(JL) = PSIGS(I1(JL),I2(JL),I3(JL)) * 2.
!     ZSIGMA_RC(JL) = MAX(PSIGS(I1(JL),I2(JL),I3(JL)) * 2., 1.E-12)
    END IF
    ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
    ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
    ZRRS(JL) = PRRS(I1(JL),I2(JL),I3(JL))
    ZRIS(JL) = PRIS(I1(JL),I2(JL),I3(JL))
    ZRSS(JL) = PRSS(I1(JL),I2(JL),I3(JL))
    ZRGS(JL) = PRGS(I1(JL),I2(JL),I3(JL))
    IF ( KRR == 7 ) ZRHS(JL) = PRHS(I1(JL),I2(JL),I3(JL))
    ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
!
    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
    ZZT(JL) = ZT(I1(JL),I2(JL),I3(JL))
    ZTHT(JL) = PTHT(I1(JL),I2(JL),I3(JL))
    ZTHLT(JL) = ZTHT(JL) - XLVTT * ZTHT(JL) / XCPD / ZZT(JL) * ZRCT(JL)
    ZPRES(JL) = PPABST(I1(JL),I2(JL),I3(JL))
    ZEXNREF(JL) = PEXNREF(I1(JL),I2(JL),I3(JL))
  ENDDO
  ALLOCATE(ZZW(IMICRO))
  ALLOCATE(ZLSFACT(IMICRO))
  ALLOCATE(ZLVFACT(IMICRO))
    ZZW(:)  = ZEXNREF(:)*( XCPD+XCPV*ZRVT(:)+XCL*(ZRCT(:)+ZRRT(:)) &
                                    +XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)) )
    ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZW(:) ! L_s/(Pi_ref*C_ph)
    ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZW(:) ! L_v/(Pi_ref*C_ph)
  ALLOCATE(ZUSW(IMICRO))
  ALLOCATE(ZSSI(IMICRO))
    ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) )
    ZSSI(:) = ZRVT(:)*( ZPRES(:)-ZZW(:) ) / ( (XMV/XMD) * ZZW(:) ) - 1.0
                                                      ! Supersaturation over ice
!
  ALLOCATE(ZLBDAR(IMICRO))
  ALLOCATE(ZLBDAR_RF(IMICRO))
  ALLOCATE(ZLBDAS(IMICRO))
  ALLOCATE(ZLBDAG(IMICRO))
  IF ( KRR == 7 ) THEN
    ALLOCATE(ZLBDAH(IMICRO))
  ELSE
    ALLOCATE(ZLBDAH(0))
  END IF
  ALLOCATE(ZRDRYG(IMICRO))
  ALLOCATE(ZRWETG(IMICRO))
  ALLOCATE(ZAI(IMICRO))
  ALLOCATE(ZCJ(IMICRO))
  ALLOCATE(ZKA(IMICRO))
  ALLOCATE(ZDV(IMICRO))
!
  IF ( KRR == 7 ) THEN
    ALLOCATE(ZZW1(IMICRO,7))
  ELSE IF( KRR == 6 ) THEN
    ALLOCATE(ZZW1(IMICRO,6))
  ENDIF
!
  IF (LBU_ENABLE .OR. LLES_CALL .OR. LCHECK ) THEN
    ALLOCATE(ZRHODJ(IMICRO))
    DO JL=1,IMICRO
      ZRHODJ(JL) = PRHODJ(I1(JL),I2(JL),I3(JL))
    END DO
  ELSE
    ALLOCATE(ZRHODJ(0))
  END IF
!

  !Cloud water split between high and low content part is done here
  !according to autoconversion option
  ZRCRAUTC(:)   = XCRIAUTC/ZRHODREF(:) ! Autoconversion rc threshold
  IF (HSUBG_AUCV == 'NONE') THEN
    !Cloud water is entirely in low or high part
    WHERE (ZRCT(:) > ZRCRAUTC(:))
      ZHLC_HCF(:) = 1.
      ZHLC_LCF(:) = 0.0
      ZHLC_HRC(:) = ZRCT(:)
      ZHLC_LRC(:) = 0.0
      ZRF(:)      = 1.
    ELSEWHERE (ZRCT(:) > XRTMIN(2))
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

  ELSEIF (HSUBG_AUCV == 'CLFR') THEN
    !Cloud water is only in the cloudy part and entirely in low or high part
      WHERE (ZCF(:) > 0. .AND. ZRCT(:) > ZRCRAUTC(:)*ZCF(:))
        ZHLC_HCF(:) = ZCF(:)
        ZHLC_LCF(:) = 0.0
        ZHLC_HRC(:) = ZRCT(:)
        ZHLC_LRC(:) = 0.0
        ZRF(:)      = ZCF(:)
      ELSEWHERE (ZCF(:) > 0. .AND. ZRCT(:) > XRTMIN(2))
        ZHLC_HCF(:) = 0.0
        ZHLC_LCF(:) = ZCF(:)
        ZHLC_HRC(:) = 0.0
        ZHLC_LRC(:) = ZRCT(:)
        ZRF(:)      = 0.
      ELSEWHERE (ZCF(:) > 0.)
        ZHLC_HCF(:) = 0.0
        ZHLC_LCF(:) = 0.0
        ZHLC_HRC(:) = 0.0
        ZHLC_LRC(:) = 0.0
        ZRF(:)      = 0.
      ELSEWHERE
        ZHLC_HCF(:) = 0.0
        ZHLC_LCF(:) = 0.0
        ZHLC_HRC(:) = 0.0
        ZHLC_LRC(:) = 0.0
        ZRF(:)      = 0.
      END WHERE

  ELSEIF (HSUBG_AUCV == 'PDF ') THEN
    !Cloud water is split between high and low part according to a PDF
    !    'HLCRECTPDF'    : rectangular PDF form
    !    'HLCTRIANGPDF'  : triangular PDF form
    !    'HLCQUADRAPDF'  : second order quadratic PDF form
    !    'HLCISOTRIPDF'  : isocele triangular PDF
    !    'SIGM'          : Redelsperger and Sommeria (1986)

    IF ( CSUBG_PR_PDF == 'SIGM' ) THEN
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
      ELSEWHERE ( ZRCT(:)>XRTMIN(2) .AND. ZCF(:)>0. )
        ZHLC_LCF(:) = 0.0
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
    ELSEIF ( CSUBG_PR_PDF== 'HLCRECTPDF' .OR. CSUBG_PR_PDF == 'HLCISOTRIPDF' .OR. &
           & CSUBG_PR_PDF == 'HLCTRIANGPDF' .OR. CSUBG_PR_PDF == 'HLCQUADRAPDF' ) THEN
      ! Calculate maximum value r_cM from PDF forms
      IF ( CSUBG_PR_PDF == 'HLCRECTPDF' .OR. CSUBG_PR_PDF == 'HLCISOTRIPDF' ) THEN
        ZCOEFFRCM = 2.0
      ELSE IF ( CSUBG_PR_PDF == 'HLCTRIANGPDF' ) THEN
        ZCOEFFRCM = 3.0
      ELSE IF ( CSUBG_PR_PDF == 'HLCQUADRAPDF' ) THEN
        ZCOEFFRCM = 4.0
      END IF
      WHERE (ZRCT(:).GT.0. .AND. ZCF(:).GT.0.)
        ZHLC_RCMAX(:) = ZCOEFFRCM * ZRCT(:) / ZCF(:)
      END WHERE

      ! Split available water and cloud fraction in two parts
      ! Calculate local mean values int he low and high parts for the 3 PDF forms:
      IF ( CSUBG_PR_PDF == 'HLCRECTPDF' ) THEN
        WHERE (ZRCT(:).GT.0. .AND. ZCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
          ZHLC_LRCLOCAL(:) = 0.5*ZRCRAUTC(:)
          ZHLC_HRCLOCAL(:) = ( ZHLC_RCMAX(:) + ZRCRAUTC(:)) / 2.0
        END WHERE
      ELSE IF ( CSUBG_PR_PDF == 'HLCTRIANGPDF' ) THEN
        WHERE (ZRCT(:).GT.0. .AND. ZCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
          ZHLC_LRCLOCAL(:) = ( ZRCRAUTC(:) *(3.0 * ZHLC_RCMAX(:) - 2.0 * ZRCRAUTC(:) ) ) &
                          / (3.0 * (2.0 * ZHLC_RCMAX(:) - ZRCRAUTC(:)  ) )
          ZHLC_HRCLOCAL(:) = (ZHLC_RCMAX(:) + 2.0*ZRCRAUTC(:)) / 3.0
        END WHERE
      ELSE IF ( CSUBG_PR_PDF == 'HLCQUADRAPDF' ) THEN
        WHERE (ZRCT(:).GT.0. .AND. ZCF(:).GT.0. .AND. ZHLC_RCMAX(:).GT.ZRCRAUTC(:))
          ZHLC_LRCLOCAL(:) = (3.0 *ZRCRAUTC(:)**3 - 8.0 *ZRCRAUTC(:)**2 * ZHLC_RCMAX(:) &
                          + 6.0*ZRCRAUTC(:) *ZHLC_RCMAX(:)**2 ) &
                          / &
                          (4.0* ZRCRAUTC(:)**2 -12.0*ZRCRAUTC(:) *ZHLC_RCMAX(:) &
                          + 12.0 * ZHLC_RCMAX(:)**2 )
          ZHLC_HRCLOCAL(:) =  (ZHLC_RCMAX(:) + 3.0*ZRCRAUTC(:)) / 4.0
        END WHERE
      ELSE IF ( CSUBG_PR_PDF == 'HLCISOTRIPDF' ) THEN
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
      WRITE(*,*) 'wrong CSUBG_PR_PDF case'
      CALL PRINT_MSG(NVERB_FATAL,'GEN','RAIN_ICE_OLD','')      
    ENDIF
  ELSE
    !wrong HSUBG_AUCV case
    WRITE(*,*)'wrong HSUBG_AUCV case'
    CALL PRINT_MSG(NVERB_FATAL,'GEN','RAIN_ICE_OLD','')  
  ENDIF

  !Diagnostic of precipitation fraction
  PRAINFR(:,:,:) = 0.
  DO JL=1,IMICRO
    PRAINFR(I1(JL),I2(JL),I3(JL)) = ZRF(JL)
  END DO
  CALL ICE4_RAINFR_VERT(D, RAIN_ICE_DESCR, PRAINFR, PRRT(:,:,:),      &
                         RESHAPE( SOURCE = [ ( 0., JL = 1, SIZE( PRSS ) ) ], SHAPE = SHAPE( PRSS ) ), &
                         RESHAPE( SOURCE = [ ( 0., JL = 1, SIZE( PRGS ) ) ], SHAPE = SHAPE( PRGS ) )  )
  DO JL=1,IMICRO
    ZRF(JL)=PRAINFR(I1(JL),I2(JL),I3(JL))
  END DO
!
  CALL RAIN_ICE_SLOW(GMICRO, ZINVTSTEP, ZRHODREF,                      &
                     ZRCT, ZRRT, ZRIT, ZRST, ZRGT, ZRHODJ, ZZT, ZPRES, &
                     ZLSFACT, ZLVFACT, ZSSI,                           &
                     ZRVS, ZRCS, ZRRS, ZRIS, ZRSS, ZRGS, ZTHS,         &
                     ZAI, ZCJ, ZKA, ZDV, ZLBDAS, ZLBDAG)
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
    ZLBDAR(:)  = XLBR*( ZRHODREF(:)*MAX( ZRRT(:),XRTMIN(3) ) )**XLBEXR
  END WHERE
  !ZLBDAR_RF will be used when we consider rain concentrated in its fraction
  WHERE( ZRRT(:)>0.0 .AND. ZRF(:)>0.0 )
    ZLBDAR_RF(:)  = XLBR*( ZRHODREF(:) *MAX( ZRRT(:)/ZRF(:)  , XRTMIN(3) ) )**XLBEXR
  ELSEWHERE
    ZLBDAR_RF(:)  = 0.
  END WHERE
!
  IF( OWARM ) THEN    !  Check if the formation of the raindrops by the slow
                      !  warm processes is allowed
    PEVAP3D(:,:,:)= 0.
    CALL RAIN_ICE_WARM(GMICRO, IMICRO, I1, I2, I3,                                                           &
                       ZRHODREF, ZRVT, ZRCT, ZRRT, ZHLC_HCF, ZHLC_LCF, ZHLC_HRC, ZHLC_LRC,                   &
                       ZRHODJ, ZPRES, ZZT, ZLBDAR, ZLBDAR_RF, ZLVFACT, ZCJ, ZKA, ZDV, ZRF, ZCF, ZTHT, ZTHLT, &
                       ZRVS, ZRCS, ZRRS, ZTHS, ZUSW, PEVAP3D)
  END IF
!
!-------------------------------------------------------------------------------
!
!
!*       4.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_s
!               ----------------------------------------------
!
  CALL RAIN_ICE_FAST_RS(PTSTEP, GMICRO, ZRHODREF, ZRVT, ZRCT, ZRRT, ZRST, ZRHODJ, ZPRES, ZZT, &
                        ZLBDAR, ZLBDAS, ZLSFACT, ZLVFACT, ZCJ, ZKA, ZDV, &
                        ZRCS, ZRRS, ZRSS, ZRGS, ZTHS)
!
!-------------------------------------------------------------------------------
!
!
!*       5.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_g
!               ----------------------------------------------
!
  CALL RAIN_ICE_FAST_RG(KRR, GMICRO, ZRHODREF, ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT, ZCIT, &
                        ZRHODJ, ZPRES, ZZT, ZLBDAR, ZLBDAS, ZLBDAG, ZLSFACT, ZLVFACT, &
                        ZCJ, ZKA, ZDV, &
                        ZRCS, ZRRS, ZRIS, ZRSS, ZRGS, ZRHS, ZTHS, &
                        ZUSW, ZRDRYG, ZRWETG)
!
!-------------------------------------------------------------------------------
!
!
!*       6.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_h
!               ----------------------------------------------
!
 IF ( KRR == 7 ) THEN
  CALL RAIN_ICE_FAST_RH(GMICRO, ZRHODREF, ZRVT, ZRCT, ZRIT, ZRST, ZRGT, ZRHT, ZRHODJ, ZPRES, &
                        ZZT, ZLBDAS, ZLBDAG, ZLBDAH, ZLSFACT, ZLVFACT, ZCJ, ZKA, ZDV, &
                        ZRCS, ZRRS, ZRIS, ZRSS, ZRGS, ZRHS, ZTHS, ZUSW)
 END IF
!
!-------------------------------------------------------------------------------
!
!
!*       7.     COMPUTES SPECIFIC SOURCES OF THE WARM AND COLD CLOUDY SPECIES
!               -------------------------------------------------------------
!
  CALL RAIN_ICE_FAST_RI(GMICRO, ZRHODREF, ZRIT, ZRHODJ, ZZT, ZSSI, ZLSFACT, ZLVFACT, &
                        ZAI, ZCJ, ZCIT, ZRCS, ZRIS, ZTHS)
!
!
!-------------------------------------------------------------------------------
!
!
!
  DO JL=1,IMICRO
    PRVS(I1(JL),I2(JL),I3(JL)) = ZRVS(JL)
    PRCS(I1(JL),I2(JL),I3(JL)) = ZRCS(JL)
    PRRS(I1(JL),I2(JL),I3(JL)) = ZRRS(JL)
    PRIS(I1(JL),I2(JL),I3(JL)) = ZRIS(JL)
    PRSS(I1(JL),I2(JL),I3(JL)) = ZRSS(JL)
    PRGS(I1(JL),I2(JL),I3(JL)) = ZRGS(JL)
    PTHS(I1(JL),I2(JL),I3(JL)) = ZTHS(JL)
    PCIT(I1(JL),I2(JL),I3(JL)) = ZCIT(JL)
    !
    PRAINFR(I1(JL),I2(JL),I3(JL)) = ZRF(JL)
  END DO
  IF ( KRR == 7 ) THEN
    DO JL=1,IMICRO
      PRHS(I1(JL),I2(JL),I3(JL)) = ZRHS(JL)
    END DO
  END IF
!
!
!
  DEALLOCATE(ZZW1)
  DEALLOCATE(ZDV)
  DEALLOCATE(ZCJ)
  DEALLOCATE(ZRDRYG)
  DEALLOCATE(ZRWETG)
  DEALLOCATE(ZLBDAG)
  DEALLOCATE(ZLBDAH)
  DEALLOCATE(ZLBDAS)
  DEALLOCATE(ZLBDAR)
  DEALLOCATE(ZLBDAR_RF)
  DEALLOCATE(ZSSI)
  DEALLOCATE(ZUSW)
  DEALLOCATE(ZLVFACT)
  DEALLOCATE(ZLSFACT)
  DEALLOCATE(ZZW)
  DEALLOCATE(ZEXNREF)
  DEALLOCATE(ZPRES)
  DEALLOCATE(ZRHODREF)
  DEALLOCATE(ZZT)
  IF(LBU_ENABLE .OR. LLES_CALL .OR. LCHECK ) DEALLOCATE(ZRHODJ)
  DEALLOCATE(ZTHS)
  DEALLOCATE(ZTHT)
  DEALLOCATE(ZTHLT)
  DEALLOCATE(ZRHS)
  DEALLOCATE(ZRGS)
  DEALLOCATE(ZRSS)
  DEALLOCATE(ZRIS)
  DEALLOCATE(ZRRS)
  DEALLOCATE(ZRCS)
  DEALLOCATE(ZRVS)
  DEALLOCATE(ZCIT)
  DEALLOCATE(ZRGT)
  DEALLOCATE(ZRHT)
  DEALLOCATE(ZRST)
  DEALLOCATE(ZRIT)
  DEALLOCATE(ZRRT)
  DEALLOCATE(ZAI)
  DEALLOCATE(ZRCT)
  DEALLOCATE(ZKA)
  DEALLOCATE(ZRVT)
  DEALLOCATE(ZSIGMA_RC)
  DEALLOCATE(ZCF)
  DEALLOCATE(ZRF)
  DEALLOCATE(ZHLC_HCF)
  DEALLOCATE(ZHLC_LCF)
  DEALLOCATE(ZHLC_HRC)
  DEALLOCATE(ZHLC_LRC)
  DEALLOCATE(ZHLC_RCMAX)
  DEALLOCATE(ZRCRAUTC)
  DEALLOCATE(ZHLC_HRCLOCAL)
  DEALLOCATE(ZHLC_LRCLOCAL)
END IF
!
!-------------------------------------------------------------------------------
!
!*       8.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!               -------------------------------------
!
!*       8.1    time splitting loop initialization
!
!
!
IF (HSEDIM == 'STAT') THEN
  CALL RAIN_ICE_SEDIMENTATION_STAT( IIB, IIE, IJB, IJE, IKB, IKE, IKTB, IKTE, IKT, KKL, KRR,                &
                                    PTSTEP, OSEDIC, PINPRC, PINDEP,                                         &
                                    PINPRR, PINPRS, PINPRG, PDZZ, PRHODREF, PPABST, PTHT, PRHODJ, PINPRR3D, &
                                    PRCS, PRCT, PRRS, PRRT, PRIS, PRSS, PRST, PRGS, PRGT,                   &
                                    PSEA, PTOWN, PINPRH, PRHS, PRHT, PFPR )
ELSEIF (HSEDIM == 'SPLI') THEN
  CALL RAIN_ICE_SEDIMENTATION_SPLIT(IIB, IIE, IJB, IJE, IKB, IKE, IKTB, IKTE, IKT, KKL,&
  KSPLITR,PTSTEP, &
  KRR,OSEDIC,LDEPOSC,PINPRC,PINDEP,PINPRR,PINPRS,PINPRG,PDZZ,PRHODREF,PPABST,PTHT,PRHODJ,&
      PINPRR3D,PRCS,PRCT,PRRS,PRRT,PRIS,PRIT,PRSS,PRST,PRGS,PRGT,PSEA,PTOWN,PINPRH,PRHS,PRHT,PFPR)
ELSE
  call Print_msg( NVERB_FATAL, 'GEN', 'RAIN_ICE_OLD', 'no sedimentation scheme for HSEDIM='//HSEDIM )
END IF
!sedimentation of rain fraction
CALL ICE4_RAINFR_VERT(D, RAIN_ICE_DESCR, PRAINFR, PRRS(:,:,:)*PTSTEP,  &
                      PRSS(:,:,:)*PTSTEP, PRGS(:,:,:)*PTSTEP)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE RAIN_ICE_OLD
