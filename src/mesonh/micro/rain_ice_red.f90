!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
       MODULE MODI_RAIN_ICE_RED
!      ########################
!
INTERFACE
      SUBROUTINE RAIN_ICE_RED ( KIT, KJT, KKT, KSIZE,                                 &
                            OSEDIC, HSEDIM, HSUBG_AUCV_RC, HSUBG_AUCV_RI, &
                            OWARM, KKA, KKU, KKL,   &
                            PTSTEP, KRR, ODMICRO, PEXN,             &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                            PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,&
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                   &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                            PINPRC,PINPRR, PINPRR3D, PEVAP3D,           &
                            PINPRS, PINPRG, PINDEP, PRAINFR, PSIGS, PSEA, PTOWN,  &
                            PRHT, PRHS, PINPRH, PFPR                              )
!
!
INTEGER,                  INTENT(IN)    :: KIT, KJT, KKT ! arrays size
INTEGER,                  INTENT(IN)    :: KSIZE
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Switch for rc->rr Subgrid autoconversion
                                        ! Kind of Subgrid autoconversion method
CHARACTER(LEN=80),        INTENT(IN)    :: HSUBG_AUCV_RI ! Switch for ri->rs Subgrid autoconversion
                                        ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
!
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)   :: ODMICRO ! mask to limit computation
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLC_HRC
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLC_HCF
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLI_HRI
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLI_HCF
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source

!
REAL, DIMENSION(:,:), INTENT(OUT)           :: PINPRC! Cloud instant precip
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(KIT,KJT,KKT), INTENT(OUT)   :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(KIT,KJT,KKT), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(:,:), INTENT(OUT)           :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(OUT) :: PRAINFR
REAL, DIMENSION(:,:,:),   INTENT(IN)        :: PSIGS   ! Sigma_s at t
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(IN)        :: PSEA ! Sea Mask
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(IN)        :: PTOWN! Fraction that is town
REAL, DIMENSION(KIT,KJT,KKT), OPTIONAL,  INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), OPTIONAL,  INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(OUT)      :: PINPRH! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)  :: PFPR ! upper-air precipitation fluxes
!
END SUBROUTINE RAIN_ICE_RED
END INTERFACE
END MODULE MODI_RAIN_ICE_RED
!     ######spl
      SUBROUTINE RAIN_ICE_RED ( KIT, KJT, KKT, KSIZE,                                 &
                            OSEDIC, HSEDIM, HSUBG_AUCV_RC, HSUBG_AUCV_RI,  &
                            OWARM,KKA,KKU,KKL,&
                            PTSTEP, KRR, ODMICRO, PEXN,                           &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                            PHLC_HRC, PHLC_HCF, PHLI_HRI,  PHLI_HCF,     &
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                   &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                            PINPRC,PINPRR, PINPRR3D, PEVAP3D,           &
                            PINPRS, PINPRG, PINDEP, PRAINFR, PSIGS, PSEA, PTOWN,  &
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
!-----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

use modd_budget,         only: lbu_enable,                                                                                     &
                               lbudget_th, lbudget_rv, lbudget_rc, lbudget_rr, lbudget_ri, lbudget_rs, lbudget_rg, lbudget_rh, &
                               NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, &
                               tbudgets
USE MODD_CST,            ONLY: XCI,XCL,XCPD,XCPV,XLSTT,XLVTT,XTT
USE MODD_PARAMETERS,     ONLY: JPVEXT,XUNDEF
USE MODD_PARAM_ICE,      ONLY: CSUBG_PR_PDF,CSUBG_RC_RR_ACCR,CSUBG_RR_EVAP,LDEPOSC,LFEEDBACKT,LSEDIM_AFTER, &
                               NMAXITER,XMRSTEP,XTSTEP_TS,XVDEPOSC
USE MODD_RAIN_ICE_DESCR, ONLY: XRTMIN
USE MODD_VAR_ll,         ONLY: IP

use mode_budget,         only: Budget_store_add, Budget_store_init, Budget_store_end
USE MODE_ll
USE MODE_MSG
use mode_tools,          only: Countjv

USE MODI_ICE4_NUCLEATION_WRAPPER
USE MODI_ICE4_RAINFR_VERT
USE MODI_ICE4_SEDIMENTATION_STAT
USE MODI_ICE4_SEDIMENTATION_SPLIT
USE MODI_ICE4_TENDENCIES

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
INTEGER,                  INTENT(IN)    :: KIT, KJT, KKT ! arrays size
INTEGER,                  INTENT(IN)    :: KSIZE
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Kind of Subgrid autoconversion method
CHARACTER(LEN=80),        INTENT(IN)    :: HSUBG_AUCV_RI ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)   :: ODMICRO ! mask to limit computation
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
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLC_HRC
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLC_HCF
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLI_HRI
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLI_HCF
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(:,:),     INTENT(OUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(KIT,KJT,KKT),INTENT(OUT)      :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(KIT,KJT,KKT), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(OUT)       :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(OUT) :: PRAINFR
REAL, DIMENSION(:,:,:),       INTENT(IN)    :: PSIGS   ! Sigma_s at t
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(IN)        :: PSEA ! Sea Mask
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(IN)        :: PTOWN! Fraction that is town
REAL, DIMENSION(KIT,KJT,KKT), OPTIONAL,  INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), OPTIONAL,  INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(OUT)      :: PINPRH! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)  :: PFPR ! upper-air precipitation fluxes
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB           !  Define the domain where is
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           !
INTEGER :: IJE           !
INTEGER :: IKB, IKTB     !
INTEGER :: IKE, IKTE     !
!
INTEGER :: JI, JJ, JK
!
!For packing
INTEGER :: IMICRO ! Case r_x>0 locations
INTEGER, DIMENSION(KSIZE) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                   :: JL       ! and PACK intrinsics
!
!Arrays for nucleation call outisde of LDMICRO points
REAL,    DIMENSION(KIT, KJT, KKT) :: ZW ! work array
REAL,    DIMENSION(KIT, KJT, KKT) :: ZT ! Temperature
REAL, DIMENSION(KIT, KJT, KKT) :: &
                                  & ZZ_RVHENI_MR, & ! heterogeneous nucleation mixing ratio change
                                  & ZZ_RVHENI       ! heterogeneous nucleation
real, dimension(:,:,:), allocatable :: zw1, zw2, zw3, zw4, zw5, zw6 !Work arrays
real, dimension(:,:,:), allocatable :: zz_diff
REAL, DIMENSION(KIT, KJT, KKT) :: ZZ_LVFACT, ZZ_LSFACT
!
!Diagnostics
REAL, DIMENSION(KIT, KJT, KKT) :: &
                                & ZHLC_HCF3D,& ! HLCLOUDS cloud fraction in high water content part
                                & ZHLC_LCF3D,& ! HLCLOUDS cloud fraction in low water content part
                                & ZHLC_HRC3D,& ! HLCLOUDS cloud water content in high water content
                                & ZHLC_LRC3D,& ! HLCLOUDS cloud water content in low water content
                                & ZHLI_HCF3D,& ! HLCLOUDS cloud fraction in high ice content part
                                & ZHLI_LCF3D,& ! HLCLOUDS cloud fraction in low ice content part
                                & ZHLI_HRI3D,& ! HLCLOUDS cloud water content in high ice content
                                & ZHLI_LRI3D   ! HLCLOUDS cloud water content in high ice content

REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2)) :: ZINPRI ! Pristine ice instant precip
!
!Packed variables
REAL, DIMENSION(KSIZE) :: ZRVT,     & ! Water vapor m.r. at t
                        & ZRCT,     & ! Cloud water m.r. at t
                        & ZRRT,     & ! Rain water m.r. at t
                        & ZRIT,     & ! Pristine ice m.r. at t
                        & ZRST,     & ! Snow/aggregate m.r. at t
                        & ZRGT,     & ! Graupel m.r. at t
                        & ZRHT,     & ! Hail m.r. at t
                        & ZCIT,     & ! Pristine ice conc. at t
                        & ZTHT,     & ! Potential temperature
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
                        & ZHLI_LRI, &
                        & ZFRAC
!
!Output packed tendencies (for budgets only)
REAL, DIMENSION(KSIZE) :: ZRVHENI_MR, & ! heterogeneous nucleation mixing ratio change
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
!For time- or mixing-ratio- splitting
REAL, DIMENSION(KSIZE) :: Z0RVT,     &   ! Water vapor m.r. at the beginig of the current loop
                        & Z0RCT,     &   ! Cloud water m.r. at the beginig of the current loop
                        & Z0RRT,     &   ! Rain water m.r. at the beginig of the current loop
                        & Z0RIT,     &   ! Pristine ice m.r. at the beginig of the current loop
                        & Z0RST,     &   ! Snow/aggregate m.r. at the beginig of the current loop
                        & Z0RGT,     &   ! Graupel m.r. at the beginig of the current loop
                        & Z0RHT,     &   ! Hail m.r. at the beginig of the current loop
                        & ZA_TH, ZA_RV, ZA_RC, ZA_RR, ZA_RI, ZA_RS, ZA_RG, ZA_RH, &
                        & ZB_TH, ZB_RV, ZB_RC, ZB_RR, ZB_RI, ZB_RS, ZB_RG, ZB_RH
!
!To take into acount external tendencies inside the splitting
REAL, DIMENSION(KSIZE) :: ZEXT_RV,   &   ! External tendencie for rv
                        & ZEXT_RC,   &   ! External tendencie for rc
                        & ZEXT_RR,   &   ! External tendencie for rr
                        & ZEXT_RI,   &   ! External tendencie for ri
                        & ZEXT_RS,   &   ! External tendencie for rs
                        & ZEXT_RG,   &   ! External tendencie for rg
                        & ZEXT_RH,   &   ! External tendencie for rh
                        & ZEXT_TH,   &   ! External tendencie for th
                        & ZEXT_WW        ! Working array
LOGICAL :: GEXT_TEND
!
INTEGER, DIMENSION(KSIZE) :: IITER ! Number of iterations done (with real tendencies computation)
INTEGER :: INB_ITER_MAX ! Maximum number of iterations (with real tendencies computation)
REAL, DIMENSION(KSIZE) :: ZTIME,    & ! Current integration time (starts with 0 and ends with PTSTEP)
                        & ZMAXTIME, & ! Time on which we can apply the current tendencies
                        & ZTIME_THRESHOLD, & ! Time to reach threshold
                        & ZTIME_LASTCALL     ! Integration time when last tendecies call has been done
REAL, DIMENSION(KSIZE) :: ZW1D
REAL, DIMENSION(KSIZE) :: ZCOMPUTE ! 1. for points where we must compute tendencies, 0. elsewhere
LOGICAL :: LSOFT ! Must we really compute tendencies or only adjust them to new T variables
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)):: GDEP
REAL :: ZTSTEP ! length of sub-timestep in case of time splitting
REAL :: ZINV_TSTEP ! Inverse ov PTSTEP
REAL, DIMENSION(KSIZE, 8) :: ZRS_TEND
REAL, DIMENSION(KSIZE, 8) :: ZRG_TEND
REAL, DIMENSION(KSIZE, 10) :: ZRH_TEND
REAL, DIMENSION(KSIZE) :: ZSSI
!
!For total tendencies computation
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: &
        &ZW_RVS, ZW_RCS, ZW_RRS, ZW_RIS, ZW_RSS, ZW_RGS, ZW_RHS, ZW_THS
!
!-------------------------------------------------------------------------------
if ( lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HENU', pths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'HENU', prvs(:, :, :) * prhodj(:, :, :) )
end if
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!               -----------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=KKA+JPVEXT*KKL
IKE=KKU-JPVEXT*KKL
IKTB=1+JPVEXT
IKTE=KKT-JPVEXT
!
ZINV_TSTEP=1./PTSTEP
GEXT_TEND=.TRUE.
!
! LSFACT and LVFACT without exner
IF(KRR==7) THEN
  DO JK = 1, KKT
    DO JJ = 1, KJT
      DO JI = 1, KIT
        ZT(JI,JJ,JK) = PTHT(JI,JJ,JK) * PEXN(JI,JJ,JK)
        ZZ_LSFACT(JI,JJ,JK)=(XLSTT+(XCPV-XCI)*(ZT(JI,JJ,JK)-XTT))   &
                         /( XCPD + XCPV*PRVT(JI,JJ,JK) + XCL*(PRCT(JI,JJ,JK)+PRRT(JI,JJ,JK))   &
                         + XCI*(PRIT(JI,JJ,JK)+PRST(JI,JJ,JK)+PRGT(JI,JJ,JK)+PRHT(JI,JJ,JK)))
        ZZ_LVFACT(JI,JJ,JK)=(XLVTT+(XCPV-XCL)*(ZT(JI,JJ,JK)-XTT))   &
                         /( XCPD + XCPV*PRVT(JI,JJ,JK) + XCL*(PRCT(JI,JJ,JK)+PRRT(JI,JJ,JK))   &
                         + XCI*(PRIT(JI,JJ,JK)+PRST(JI,JJ,JK)+PRGT(JI,JJ,JK)+PRHT(JI,JJ,JK)))
      ENDDO
    ENDDO
  ENDDO
ELSE
  DO JK = 1, KKT
    DO JJ = 1, KJT
      DO JI = 1, KIT
        ZT(JI,JJ,JK) = PTHT(JI,JJ,JK) * PEXN(JI,JJ,JK)
        ZZ_LSFACT(JI,JJ,JK)=(XLSTT+(XCPV-XCI)*(ZT(JI,JJ,JK)-XTT))   &
                         /( XCPD + XCPV*PRVT(JI,JJ,JK) + XCL*(PRCT(JI,JJ,JK)+PRRT(JI,JJ,JK))   &
                         + XCI*(PRIT(JI,JJ,JK)+PRST(JI,JJ,JK)+PRGT(JI,JJ,JK)))
        ZZ_LVFACT(JI,JJ,JK)=(XLVTT+(XCPV-XCL)*(ZT(JI,JJ,JK)-XTT))   &
                         /( XCPD + XCPV*PRVT(JI,JJ,JK) + XCL*(PRCT(JI,JJ,JK)+PRRT(JI,JJ,JK))   &
                         + XCI*(PRIT(JI,JJ,JK)+PRST(JI,JJ,JK)+PRGT(JI,JJ,JK)))
      ENDDO
    ENDDO
  ENDDO
ENDIF
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
  if ( lbudget_rc .and. osedic ) call Budget_store_init( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr  )             call Budget_store_init( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri  )             call Budget_store_init( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs  )             call Budget_store_init( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg  )             call Budget_store_init( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh  )             call Budget_store_init( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )

  !Init only if not osedic (to prevent crash with double init)
  !Remark: the 2 source terms SEDI and DEPO could be mixed and stored in the same source term (SEDI)
  !        if osedic=T and ldeposc=T (a warning is printed in ini_budget in that case)
  if ( lbudget_rc .and. ldeposc .and. .not.osedic ) &
    call Budget_store_init( tbudgets(NBUDGET_RC), 'DEPO', prcs(:, :, :) * prhodj(:, :, :) )

  IF(HSEDIM=='STAT') THEN
    !SR: It *seems* that we must have two separate calls for ifort
    IF(KRR==7) THEN
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCS*PTSTEP, PRRS, PRRS*PTSTEP, PRIS, PRIS*PTSTEP,&
                                  &PRSS, PRSS*PTSTEP, PRGS, PRGS*PTSTEP,&
                                  &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA=PSEA, PTOWN=PTOWN, &
                                  &PINPRH=PINPRH, PRHT=PRHS*PTSTEP, PRHS=PRHS, PFPR=PFPR)
    ELSE
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCS*PTSTEP, PRRS, PRRS*PTSTEP, PRIS, PRIS*PTSTEP,&
                                  &PRSS, PRSS*PTSTEP, PRGS, PRGS*PTSTEP,&
                                  &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA=PSEA, PTOWN=PTOWN, &
                                  &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !No negativity correction here as we apply sedimentation on PR.S*PTSTEP variables
  ELSEIF(HSEDIM=='SPLI') THEN
    !SR: It *seems* that we must have two separate calls for ifort
    IF(KRR==7) THEN
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA=PSEA, PTOWN=PTOWN, &
                                   &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR)
    ELSE
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA=PSEA, PTOWN=PTOWN, &
                                   &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !We correct negativities with conservation
    !SPLI algorith uses a time-splitting. Inside the loop a temporary m.r. is used.
    !   It is initialized with the m.r. at T and is modified by two tendencies:
    !   sedimentation tendency and an external tendency which represents all other
    !   processes (mainly advection and microphysical processes). If both tendencies
    !   are negative, sedimentation can remove a specie at a given sub-timestep. From
    !   this point sedimentation stops for the remaining sub-timesteps but the other tendency
    !   will be still active and will lead to negative values.
    !   We could prevent the algorithm to not consume too much a specie, instead we apply
    !   a correction here.
    CALL CORRECT_NEGATIVITIES(KIT, KJT, KKT, KRR, PRVS, PRCS, PRRS, &
                             &PRIS, PRSS, PRGS, &
                             &PTHS, ZZ_LVFACT, ZZ_LSFACT, PRHS)
  ELSEIF(HSEDIM=='NONE') THEN
  ELSE
    call Print_msg( NVERB_FATAL, 'GEN', 'RAIN_ICE_RED', 'no sedimentation scheme for HSEDIM='//HSEDIM )
  END IF
  !
  !*       2.2     budget storage
  !
  if ( lbudget_rc .and. osedic ) call Budget_store_end( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr  )             call Budget_store_end( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri  )             call Budget_store_end( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs  )             call Budget_store_end( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg  )             call Budget_store_end( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh  )             call Budget_store_end( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )

  !If osedic=T and ldeposc=T, DEPO is in fact mixed and stored with the SEDI source term
  !(a warning is printed in ini_budget in that case)
  if ( lbudget_rc .and. ldeposc .and. .not.osedic) &
    call Budget_store_end( tbudgets(NBUDGET_RC), 'DEPO', prcs(:, :, :) * prhodj(:, :, :) )
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.     PACKING
!               --------
!  optimization by looking for locations where
!  the microphysical fields are larger than a minimal value only !!!
!
IMICRO=0
IF(KSIZE/=0) IMICRO=COUNTJV(ODMICRO(:,:,:), I1(:), I2(:), I3(:))
!Packing
IF(IMICRO>0) THEN
  DO JL=1, IMICRO
    ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
    ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
    ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
    ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
    ZRST(JL) = PRST(I1(JL),I2(JL),I3(JL))
    ZRGT(JL) = PRGT(I1(JL),I2(JL),I3(JL))
    ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
    ZCF(JL) = PCLDFR(I1(JL),I2(JL),I3(JL))
    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
    ZTHT(JL) = PTHT(I1(JL),I2(JL),I3(JL))
    ZPRES(JL) = PPABST(I1(JL),I2(JL),I3(JL))
    ZEXN(JL) = PEXN(I1(JL),I2(JL),I3(JL))
    ZHLC_HCF(JL) = PHLC_HCF(I1(JL),I2(JL),I3(JL))
    ZHLC_HRC(JL) = PHLC_HRC(I1(JL),I2(JL),I3(JL))
    ZHLC_LRC(JL) = ZRCT(JL) - ZHLC_HRC(JL)
    ZHLI_HCF(JL) = PHLI_HCF(I1(JL),I2(JL),I3(JL))
    ZHLI_HRI(JL) = PHLI_HRI(I1(JL),I2(JL),I3(JL))
    ZHLI_LRI(JL) = ZRIT(JL) - ZHLI_HRI(JL)
    IF(ZRCT(JL)>0.) THEN
      ZHLC_LCF(JL) = ZCF(JL)- ZHLC_HCF(JL)
    ELSE
      ZHLC_LCF(JL)=0.
    ENDIF
    IF(ZRIT(JL)>0.) THEN
      ZHLI_LCF(JL) = ZCF(JL)- ZHLI_HCF(JL)
    ELSE
      ZHLI_LCF(JL)=0.
    ENDIF
  ENDDO
  IF(GEXT_TEND) THEN
    DO JL=1, IMICRO
      ZEXT_RV(JL) = PRVS(I1(JL),I2(JL),I3(JL)) - ZRVT(JL)*ZINV_TSTEP
      ZEXT_RC(JL) = PRCS(I1(JL),I2(JL),I3(JL)) - ZRCT(JL)*ZINV_TSTEP
      ZEXT_RR(JL) = PRRS(I1(JL),I2(JL),I3(JL)) - ZRRT(JL)*ZINV_TSTEP
      ZEXT_RI(JL) = PRIS(I1(JL),I2(JL),I3(JL)) - ZRIT(JL)*ZINV_TSTEP
      ZEXT_RS(JL) = PRSS(I1(JL),I2(JL),I3(JL)) - ZRST(JL)*ZINV_TSTEP
      ZEXT_RG(JL) = PRGS(I1(JL),I2(JL),I3(JL)) - ZRGT(JL)*ZINV_TSTEP
      ZEXT_TH(JL) = PTHS(I1(JL),I2(JL),I3(JL)) - ZTHT(JL)*ZINV_TSTEP
      !The th tendency is not related to a mixing ratio change, there is no exn/exnref issue here
    ENDDO
  ENDIF
  IF(HSUBG_AUCV_RC=='PDF ' .AND. CSUBG_PR_PDF=='SIGM') THEN
    DO JL=1, IMICRO
      ZSIGMA_RC(JL) = PSIGS(I1(JL),I2(JL),I3(JL))*2.
    ENDDO
  ENDIF
  IF(KRR==7) THEN
    DO JL=1, IMICRO
      ZRHT(JL) = PRHT(I1(JL),I2(JL),I3(JL))
    ENDDO
    IF(GEXT_TEND) THEN
      DO JL=1, IMICRO
        ZEXT_RH(JL) = PRHS(I1(JL),I2(JL),I3(JL)) - ZRHT(JL)*ZINV_TSTEP
      ENDDO
    ENDIF
  ELSE
    ZRHT(:)=0.
    IF(GEXT_TEND) ZEXT_RH(:)=0.
  ENDIF
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
ENDIF
!-------------------------------------------------------------------------------
!
!*       4.     LOOP
!               ----
!
!Maximum number of iterations
!We only count real iterations (those for which we *compute* tendencies)
INB_ITER_MAX=NMAXITER
IF(XTSTEP_TS/=0.)THEN
  INB_ITER_MAX=MAX(1, INT(PTSTEP/XTSTEP_TS)) !At least the number of iterations needed for the time-splitting
  ZTSTEP=PTSTEP/INB_ITER_MAX
  INB_ITER_MAX=MAX(NMAXITER, INB_ITER_MAX) !For the case XMRSTEP/=0. at the same time
ENDIF
IITER(:)=0
ZTIME(:)=0. ! Current integration time (all points may have a different integration time)
DO WHILE(ANY(ZTIME(:)<PTSTEP)) ! Loop to *really* compute tendencies
  IF(XMRSTEP/=0.) THEN
    ! In this case we need to remember the mixing ratios used to compute the tendencies
    ! because when mixing ratio has evolved more than a threshold, we must re-compute tendecies
    DO JL=1, IMICRO
      Z0RVT(JL)=ZRVT(JL)
      Z0RCT(JL)=ZRCT(JL)
      Z0RRT(JL)=ZRRT(JL)
      Z0RIT(JL)=ZRIT(JL)
      Z0RST(JL)=ZRST(JL)
      Z0RGT(JL)=ZRGT(JL)
      Z0RHT(JL)=ZRHT(JL)
    ENDDO
  ENDIF
  IF(XTSTEP_TS/=0.) THEN
    ! In this case we need to remember the time when tendencies were computed
    ! because when time has evolved more than a limit, we must re-compute tendecies
    ZTIME_LASTCALL(:)=ZTIME(:)
  ENDIF
  ZCOMPUTE(:)=MAX(0., -SIGN(1., ZTIME(:)-PTSTEP)) ! Compuation (1.) only for points for which integration time has not reached the timestep
  LSOFT=.FALSE. ! We *really* compute the tendencies
  IITER(:)=IITER(:)+INT(ZCOMPUTE(:))
  DO WHILE(SUM(ZCOMPUTE(:))>0.) ! Loop to adjust tendencies when we cross the 0°C or when a specie disappears
    IF(KRR==7) THEN
      DO JL=1, IMICRO
        ZZT(JL) = ZTHT(JL) * ZEXN(JL)
        ZLSFACT(JL)=(XLSTT+(XCPV-XCI)*(ZZT(JL)-XTT))   &
                   &/( (XCPD + XCPV*ZRVT(JL) + XCL*(ZRCT(JL)+ZRRT(JL))   &
                   &+ XCI*(ZRIT(JL)+ZRST(JL)+ZRGT(JL)+ZRHT(JL)))*ZEXN(JL) )
        ZLVFACT(JL)=(XLVTT+(XCPV-XCL)*(ZZT(JL)-XTT))   &
                   &/( (XCPD + XCPV*ZRVT(JL) + XCL*(ZRCT(JL)+ZRRT(JL))   &
                   &+ XCI*(ZRIT(JL)+ZRST(JL)+ZRGT(JL)+ZRHT(JL)))*ZEXN(JL) )
      ENDDO
    ELSE
      DO JL=1, IMICRO
        ZZT(JL) = ZTHT(JL) * ZEXN(JL)
        ZLSFACT(JL)=(XLSTT+(XCPV-XCI)*(ZZT(JL)-XTT))   &
                   &/( (XCPD + XCPV*ZRVT(JL) + XCL*(ZRCT(JL)+ZRRT(JL))   &
                   &+ XCI*(ZRIT(JL)+ZRST(JL)+ZRGT(JL)))*ZEXN(JL) )
        ZLVFACT(JL)=(XLVTT+(XCPV-XCL)*(ZZT(JL)-XTT))   &
                   &/( (XCPD + XCPV*ZRVT(JL) + XCL*(ZRCT(JL)+ZRRT(JL))   &
                   &+ XCI*(ZRIT(JL)+ZRST(JL)+ZRGT(JL)))*ZEXN(JL) )
      ENDDO
    ENDIF
    !
    !***       4.1 Tendecies computation
    !
    ! Tendencies are *really* computed when LSOFT==.FALSE. and only adjusted otherwise
    CALL ICE4_TENDENCIES(IMICRO, IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, KKT, KKL, &
                        &KRR, LSOFT, ZCOMPUTE, &
                        &OWARM, CSUBG_RC_RR_ACCR, CSUBG_RR_EVAP, &
                        &HSUBG_AUCV_RC, HSUBG_AUCV_RI, CSUBG_PR_PDF, &
                        &ZEXN, ZRHODREF, ZLVFACT, ZLSFACT, I1, I2, I3, &
                        &ZPRES, ZCF, ZSIGMA_RC,&
                        &ZCIT, &
                        &ZZT, ZTHT, &
                        &ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT, ZRHT, &
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
                        &ZA_TH, ZA_RV, ZA_RC, ZA_RR, ZA_RI, ZA_RS, ZA_RG, ZA_RH, &
                        &ZB_TH, ZB_RV, ZB_RC, ZB_RR, ZB_RI, ZB_RS, ZB_RG, ZB_RH, &
                        &ZHLC_HCF, ZHLC_LCF, ZHLC_HRC, ZHLC_LRC, &
                        &ZHLI_HCF, ZHLI_LCF, ZHLI_HRI, ZHLI_LRI, PRAINFR)
    ! External tendencies
    IF(GEXT_TEND) THEN
      DO JL=1, IMICRO
        ZA_TH(JL) = ZA_TH(JL) + ZEXT_TH(JL)
        ZA_RV(JL) = ZA_RV(JL) + ZEXT_RV(JL)
        ZA_RC(JL) = ZA_RC(JL) + ZEXT_RC(JL)
        ZA_RR(JL) = ZA_RR(JL) + ZEXT_RR(JL)
        ZA_RI(JL) = ZA_RI(JL) + ZEXT_RI(JL)
        ZA_RS(JL) = ZA_RS(JL) + ZEXT_RS(JL)
        ZA_RG(JL) = ZA_RG(JL) + ZEXT_RG(JL)
        ZA_RH(JL) = ZA_RH(JL) + ZEXT_RH(JL)
      ENDDO
    ENDIF
    !
    !***       4.2 Integration time
    !
    ! If we can, we will use these tendencies until the end of the timestep
    ZMAXTIME(:)=ZCOMPUTE(:) * (PTSTEP-ZTIME(:)) ! Remaining time until the end of the timestep

    !We need to adjust tendencies when temperature reaches 0
    IF(LFEEDBACKT) THEN
      DO JL=1, IMICRO
        !Is ZB_TH enough to change temperature sign?
        ZW1D(JL)=(ZTHT(JL) - XTT/ZEXN(JL)) * (ZTHT(JL) + ZB_TH(JL) - XTT/ZEXN(JL))
        ZMAXTIME(JL)=ZMAXTIME(JL)*MAX(0., SIGN(1., ZW1D(JL)))
        !Can ZA_TH make temperature change of sign?
        ZW1D(JL)=MAX(0., -SIGN(1., 1.E-20 - ABS(ZA_TH(JL)))) ! WHERE(ABS(ZA_TH(:))>1.E-20)
        ZTIME_THRESHOLD(JL)=(1. - ZW1D(JL))*(-1.) + &
                            ZW1D(JL) * &
                            (XTT/ZEXN(JL) - ZB_TH(JL) - ZTHT(JL))/ &
                            SIGN(MAX(ABS(ZA_TH(JL)), 1.E-20), ZA_TH(JL))
        ZW1D(JL)=MAX(0., -SIGN(1., 1.E-20 - ZTIME_THRESHOLD(JL))) ! WHERE(ZTIME_THRESHOLD(:)>1.E-20)
        ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                     ZW1D(JL) * MIN(ZMAXTIME(JL), ZTIME_THRESHOLD(JL))
      ENDDO
    ENDIF

    !We need to adjust tendencies when a specy disappears
    !When a species is missing, only the external tendencies can be negative (and we must keep track of it)
    DO JL=1, IMICRO
      ZW1D(JL)=MAX(0., -SIGN(1., ZA_RV(JL)+1.E-20)) * & ! WHERE(ZA_RV(:)<-1.E-20)
              &MAX(0., -SIGN(1., XRTMIN(1)-ZRVT(JL)))   ! WHERE(ZRVT(:)>XRTMIN(1))
      ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                  &ZW1D(JL) * MIN(ZMAXTIME(JL), -(ZB_RV(JL)+ZRVT(JL))/MIN(ZA_RV(JL), -1.E-20))

      ZW1D(JL)=MAX(0., -SIGN(1., ZA_RC(JL)+1.E-20)) * & ! WHERE(ZA_RC(:)<-1.E-20)
              &MAX(0., -SIGN(1., XRTMIN(2)-ZRCT(JL)))   ! WHERE(ZRCT(:)>XRTMIN(2))
      ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                  &ZW1D(JL) * MIN(ZMAXTIME(JL), -(ZB_RC(JL)+ZRCT(JL))/MIN(ZA_RC(JL), -1.E-20))

      ZW1D(JL)=MAX(0., -SIGN(1., ZA_RR(JL)+1.E-20)) * & ! WHERE(ZA_RR(:)<-1.E-20)
              &MAX(0., -SIGN(1., XRTMIN(3)-ZRRT(JL)))   ! WHERE(ZRRT(:)>XRTMIN(3))
      ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                  &ZW1D(JL) * MIN(ZMAXTIME(JL), -(ZB_RR(JL)+ZRRT(JL))/MIN(ZA_RR(JL), -1.E-20))

      ZW1D(JL)=MAX(0., -SIGN(1., ZA_RI(JL)+1.E-20)) * & ! WHERE(ZI_RV(:)<-1.E-20)
              &MAX(0., -SIGN(1., XRTMIN(4)-ZRIT(JL)))   ! WHERE(ZRIT(:)>XRTMIN(4))
      ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                  &ZW1D(JL) * MIN(ZMAXTIME(JL), -(ZB_RI(JL)+ZRIT(JL))/MIN(ZA_RI(JL), -1.E-20))

      ZW1D(JL)=MAX(0., -SIGN(1., ZA_RS(JL)+1.E-20)) * & ! WHERE(ZA_RS(:)<-1.E-20)
              &MAX(0., -SIGN(1., XRTMIN(5)-ZRST(JL)))   ! WHERE(ZRST(:)>XRTMIN(5))
      ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                  &ZW1D(JL) * MIN(ZMAXTIME(JL), -(ZB_RS(JL)+ZRST(JL))/MIN(ZA_RS(JL), -1.E-20))

      ZW1D(JL)=MAX(0., -SIGN(1., ZA_RG(JL)+1.E-20)) * & ! WHERE(ZA_RG(:)<-1.E-20)
              &MAX(0., -SIGN(1., XRTMIN(6)-ZRGT(JL)))   ! WHERE(ZRGT(:)>XRTMIN(6))
      ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                  &ZW1D(JL) * MIN(ZMAXTIME(JL), -(ZB_RG(JL)+ZRGT(JL))/MIN(ZA_RG(JL), -1.E-20))
    ENDDO

    IF(KRR==7) THEN
      DO JL=1, IMICRO
        ZW1D(JL)=MAX(0., -SIGN(1., ZA_RH(JL)+1.E-20)) * & ! WHERE(ZA_RH(:)<-1.E-20)
                &MAX(0., -SIGN(1., XRTMIN(7)-ZRHT(JL)))   ! WHERE(ZRHT(:)>XRTMIN(7))
        ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                    &ZW1D(JL) * MIN(ZMAXTIME(JL), -(ZB_RH(JL)+ZRHT(JL))/MIN(ZA_RH(JL), -1.E-20))
      ENDDO
    ENDIF

    !We stop when the end of the timestep is reached
    ZCOMPUTE(:)=ZCOMPUTE(:) * MAX(0., -SIGN(1., ZTIME(:)+ZMAXTIME(:)-PTSTEP))

    !We must recompute tendencies when the end of the sub-timestep is reached
    IF(XTSTEP_TS/=0.) THEN
      DO JL=1, IMICRO
        ZW1D(JL)=MAX(0., -SIGN(1., IITER(JL)-INB_ITER_MAX+0.)) * & ! WHERE(IITER(:)<INB_ITER_MAX)
                &MAX(0., -SIGN(1., ZTIME_LASTCALL(JL)+ZTSTEP-ZTIME(JL)-ZMAXTIME(JL))) ! WHERE(ZTIME(:)+ZMAXTIME(:)>ZTIME_LASTCALL(:)+ZTSTEP)
        ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                    &ZW1D(JL) * (ZTIME_LASTCALL(JL)-ZTIME(JL)+ZTSTEP)
        ZCOMPUTE(JL)=ZCOMPUTE(JL) * (1. - ZW1D(JL))
      ENDDO
    ENDIF

    !We must recompute tendencies when the maximum allowed change is reached
    !When a specy is missing, only the external tendencies can be active and we do not want to recompute
    !the microphysical tendencies when external tendencies are negative (results won't change because specy was already missing)
    IF(XMRSTEP/=0.) THEN
      DO JL=1, IMICRO
        ZW1D(JL)=MAX(0., -SIGN(1., IITER(JL)-INB_ITER_MAX+0.)) * & ! WHERE(IITER(:)<INB_ITER_MAX)
                &MAX(0., -SIGN(1., 1.E-20-ABS(ZA_RV(JL)))) ! WHERE(ABS(ZA_RV(:))>1.E-20)
        ZTIME_THRESHOLD(JL)=(1.-ZW1D(JL))*(-1.) + &
                           &ZW1D(JL)*(SIGN(1., ZA_RV(JL))*XMRSTEP+Z0RVT(JL)-ZRVT(JL)-ZB_RV(JL))/ &
                           &SIGN(MAX(ABS(ZA_RV(JL)), 1.E-20), ZA_RV(JL))
        ZW1D(JL)=MAX(0., SIGN(1., ZTIME_THRESHOLD(JL))) * & !WHERE(ZTIME_THRESHOLD(:)>=0.)
                &MAX(0., -SIGN(1., ZTIME_THRESHOLD(JL)-ZMAXTIME(JL))) * & !WHERE(ZTIME_THRESHOLD(:)<ZMAXTIME(:))
                &MIN(1., MAX(0., -SIGN(1., XRTMIN(6)-ZRVT(JL))) + & !WHERE(ZRVT(:)>XRTMIN(6)) .OR.
                        &MAX(0., -SIGN(1., -ZA_RV(JL))))            !WHERE(ZA_RV(:)>0.)
        ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                    &ZW1D(JL)*MIN(ZMAXTIME(JL), ZTIME_THRESHOLD(JL))
        ZCOMPUTE(JL)=ZCOMPUTE(JL) * (1. - ZW1D(JL))

        ZW1D(JL)=MAX(0., -SIGN(1., IITER(JL)-INB_ITER_MAX+0.)) * & ! WHERE(IITER(:)<INB_ITER_MAX)
                &MAX(0., -SIGN(1., 1.E-20-ABS(ZA_RC(JL)))) ! WHERE(ABS(ZA_RC(:))>1.E-20)
        ZTIME_THRESHOLD(JL)=(1.-ZW1D(JL))*(-1.) + &
                           &ZW1D(JL)*(SIGN(1., ZA_RC(JL))*XMRSTEP+Z0RCT(JL)-ZRCT(JL)-ZB_RC(JL))/ &
                           &SIGN(MAX(ABS(ZA_RC(JL)), 1.E-20), ZA_RC(JL))
        ZW1D(JL)=MAX(0., SIGN(1., ZTIME_THRESHOLD(JL))) * & !WHERE(ZTIME_THRESHOLD(:)>=0.)
                &MAX(0., -SIGN(1., ZTIME_THRESHOLD(JL)-ZMAXTIME(JL))) * & !WHERE(ZTIME_THRESHOLD(:)<ZMAXTIME(:))
                &MIN(1., MAX(0., -SIGN(1., XRTMIN(6)-ZRCT(JL))) + & !WHERE(ZRCT(:)>XRTMIN(6)) .OR.
                        &MAX(0., -SIGN(1., -ZA_RC(JL))))            !WHERE(ZA_RC(:)>0.)
        ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                    &ZW1D(JL)*MIN(ZMAXTIME(JL), ZTIME_THRESHOLD(JL))
        ZCOMPUTE(JL)=ZCOMPUTE(JL) * (1. - ZW1D(JL))

        ZW1D(JL)=MAX(0., -SIGN(1., IITER(JL)-INB_ITER_MAX+0.)) * & ! WHERE(IITER(:)<INB_ITER_MAX)
                &MAX(0., -SIGN(1., 1.E-20-ABS(ZA_RR(JL)))) ! WHERE(ABS(ZA_RR(:))>1.E-20)
        ZTIME_THRESHOLD(JL)=(1.-ZW1D(JL))*(-1.) + &
                           &ZW1D(JL)*(SIGN(1., ZA_RR(JL))*XMRSTEP+Z0RRT(JL)-ZRRT(JL)-ZB_RR(JL))/ &
                           &SIGN(MAX(ABS(ZA_RR(JL)), 1.E-20), ZA_RR(JL))
        ZW1D(JL)=MAX(0., SIGN(1., ZTIME_THRESHOLD(JL))) * & !WHERE(ZTIME_THRESHOLD(:)>=0.)
                &MAX(0., -SIGN(1., ZTIME_THRESHOLD(JL)-ZMAXTIME(JL))) * & !WHERE(ZTIME_THRESHOLD(:)<ZMAXTIME(:))
                &MIN(1., MAX(0., -SIGN(1., XRTMIN(6)-ZRRT(JL))) + & !WHERE(ZRRT(:)>XRTMIN(6)) .OR.
                        &MAX(0., -SIGN(1., -ZA_RR(JL))))            !WHERE(ZA_RR(:)>0.)
        ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                    &ZW1D(JL)*MIN(ZMAXTIME(JL), ZTIME_THRESHOLD(JL))
        ZCOMPUTE(JL)=ZCOMPUTE(JL) * (1. - ZW1D(JL))

        ZW1D(JL)=MAX(0., -SIGN(1., IITER(JL)-INB_ITER_MAX+0.)) * & ! WHERE(IITER(:)<INB_ITER_MAX)
                &MAX(0., -SIGN(1., 1.E-20-ABS(ZA_RI(JL)))) ! WHERE(ABS(ZA_RI(:))>1.E-20)
        ZTIME_THRESHOLD(JL)=(1.-ZW1D(JL))*(-1.) + &
                           &ZW1D(JL)*(SIGN(1., ZA_RI(JL))*XMRSTEP+Z0RIT(JL)-ZRIT(JL)-ZB_RI(JL))/ &
                           &SIGN(MAX(ABS(ZA_RI(JL)), 1.E-20), ZA_RI(JL))
        ZW1D(JL)=MAX(0., SIGN(1., ZTIME_THRESHOLD(JL))) * & !WHERE(ZTIME_THRESHOLD(:)>=0.)
                &MAX(0., -SIGN(1., ZTIME_THRESHOLD(JL)-ZMAXTIME(JL))) * & !WHERE(ZTIME_THRESHOLD(:)<ZMAXTIME(:))
                &MIN(1., MAX(0., -SIGN(1., XRTMIN(6)-ZRIT(JL))) + & !WHERE(ZRIT(:)>XRTMIN(6)) .OR.
                        &MAX(0., -SIGN(1., -ZA_RI(JL))))            !WHERE(ZA_RI(:)>0.)
        ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                    &ZW1D(JL)*MIN(ZMAXTIME(JL), ZTIME_THRESHOLD(JL))
        ZCOMPUTE(JL)=ZCOMPUTE(JL) * (1. - ZW1D(JL))

        ZW1D(JL)=MAX(0., -SIGN(1., IITER(JL)-INB_ITER_MAX+0.)) * & ! WHERE(IITER(:)<INB_ITER_MAX)
                &MAX(0., -SIGN(1., 1.E-20-ABS(ZA_RS(JL)))) ! WHERE(ABS(ZA_RS(:))>1.E-20)
        ZTIME_THRESHOLD(JL)=(1.-ZW1D(JL))*(-1.) + &
                           &ZW1D(JL)*(SIGN(1., ZA_RS(JL))*XMRSTEP+Z0RST(JL)-ZRST(JL)-ZB_RS(JL))/ &
                           &SIGN(MAX(ABS(ZA_RS(JL)), 1.E-20), ZA_RS(JL))
        ZW1D(JL)=MAX(0., SIGN(1., ZTIME_THRESHOLD(JL))) * & !WHERE(ZTIME_THRESHOLD(:)>=0.)
                &MAX(0., -SIGN(1., ZTIME_THRESHOLD(JL)-ZMAXTIME(JL))) * & !WHERE(ZTIME_THRESHOLD(:)<ZMAXTIME(:))
                &MIN(1., MAX(0., -SIGN(1., XRTMIN(6)-ZRST(JL))) + & !WHERE(ZRST(:)>XRTMIN(6)) .OR.
                        &MAX(0., -SIGN(1., -ZA_RS(JL))))            !WHERE(ZA_RS(:)>0.)
        ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                    &ZW1D(JL)*MIN(ZMAXTIME(JL), ZTIME_THRESHOLD(JL))
        ZCOMPUTE(JL)=ZCOMPUTE(JL) * (1. - ZW1D(JL))

        ZW1D(JL)=MAX(0., -SIGN(1., IITER(JL)-INB_ITER_MAX+0.)) * & ! WHERE(IITER(:)<INB_ITER_MAX)
                &MAX(0., -SIGN(1., 1.E-20-ABS(ZA_RG(JL)))) ! WHERE(ABS(ZA_RG(:))>1.E-20)
        ZTIME_THRESHOLD(JL)=(1.-ZW1D(JL))*(-1.) + &
                           &ZW1D(JL)*(SIGN(1., ZA_RG(JL))*XMRSTEP+Z0RGT(JL)-ZRGT(JL)-ZB_RG(JL))/ &
                           &SIGN(MAX(ABS(ZA_RG(JL)), 1.E-20), ZA_RG(JL))
        ZW1D(JL)=MAX(0., SIGN(1., ZTIME_THRESHOLD(JL))) * & !WHERE(ZTIME_THRESHOLD(:)>=0.)
                &MAX(0., -SIGN(1., ZTIME_THRESHOLD(JL)-ZMAXTIME(JL))) * & !WHERE(ZTIME_THRESHOLD(:)<ZMAXTIME(:))
                &MIN(1., MAX(0., -SIGN(1., XRTMIN(6)-ZRGT(JL))) + & !WHERE(ZRGT(:)>XRTMIN(6)) .OR.
                        &MAX(0., -SIGN(1., -ZA_RG(JL))))            !WHERE(ZA_RG(:)>0.)
        ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                    &ZW1D(JL)*MIN(ZMAXTIME(JL), ZTIME_THRESHOLD(JL))
        ZCOMPUTE(JL)=ZCOMPUTE(JL) * (1. - ZW1D(JL))
      ENDDO

      IF(KRR==7) THEN
        DO JL=1, IMICRO
          ZW1D(JL)=MAX(0., -SIGN(1., IITER(JL)-INB_ITER_MAX+0.)) * & ! WHERE(IITER(:)<INB_ITER_MAX)
                  &MAX(0., -SIGN(1., 1.E-20-ABS(ZA_RH(JL)))) ! WHERE(ABS(ZA_RH(:))>1.E-20)
          ZTIME_THRESHOLD(JL)=(1.-ZW1D(JL))*(-1.) + &
                             &ZW1D(JL)*(SIGN(1., ZA_RH(JL))*XMRSTEP+Z0RHT(JL)-ZRHT(JL)-ZB_RH(JL))/ &
                             &SIGN(MAX(ABS(ZA_RH(JL)), 1.E-20), ZA_RH(JL))
          ZW1D(JL)=MAX(0., SIGN(1., ZTIME_THRESHOLD(JL))) * & !WHERE(ZTIME_THRESHOLD(:)>=0.)
                  &MAX(0., -SIGN(1., ZTIME_THRESHOLD(JL)-ZMAXTIME(JL))) * & !WHERE(ZTIME_THRESHOLD(:)<ZMAXTIME(:))
                  &MIN(1., MAX(0., -SIGN(1., XRTMIN(6)-ZRHT(JL))) + & !WHERE(ZRHT(:)>XRTMIN(6)) .OR.
                          &MAX(0., -SIGN(1., -ZA_RH(JL))))            !WHERE(ZA_RH(:)>0.)
          ZMAXTIME(JL)=(1.-ZW1D(JL)) * ZMAXTIME(JL) + &
                      &ZW1D(JL)*MIN(ZMAXTIME(JL), ZTIME_THRESHOLD(JL))
          ZCOMPUTE(JL)=ZCOMPUTE(JL) * (1. - ZW1D(JL))
        ENDDO
      ENDIF

      DO JL=1, IMICRO
        ZW1D(JL)=MAX(ABS(ZB_RV(JL)), ABS(ZB_RC(JL)), ABS(ZB_RR(JL)), ABS(ZB_RI(JL)), &
                    &ABS(ZB_RS(JL)), ABS(ZB_RG(JL)), ABS(ZB_RH(JL)))
        ZW1D(JL)=MAX(0., -SIGN(1., IITER(JL)-INB_ITER_MAX+0.)) * & !WHERE(IITER(:)<INB_ITER_MAX)
                &MAX(0., -SIGN(1., XMRSTEP-ZW1D(JL))) !WHERE(ZW1D(:)>XMRSTEP)
        ZMAXTIME(JL)=(1.-ZW1D(JL))*ZMAXTIME(JL)
        ZCOMPUTE(JL)=ZCOMPUTE(JL) * (1. - ZW1D(JL))
      ENDDO
    ENDIF
    !
    !***       4.3 New values of variables for next iteration
    !
    DO JL=1, IMICRO
      ZTHT(JL)=ZTHT(JL)+ZA_TH(JL)*ZMAXTIME(JL)+ZB_TH(JL)
      ZRVT(JL)=ZRVT(JL)+ZA_RV(JL)*ZMAXTIME(JL)+ZB_RV(JL)
      ZRCT(JL)=ZRCT(JL)+ZA_RC(JL)*ZMAXTIME(JL)+ZB_RC(JL)
      ZRRT(JL)=ZRRT(JL)+ZA_RR(JL)*ZMAXTIME(JL)+ZB_RR(JL)
      ZRIT(JL)=ZRIT(JL)+ZA_RI(JL)*ZMAXTIME(JL)+ZB_RI(JL)
      ZRST(JL)=ZRST(JL)+ZA_RS(JL)*ZMAXTIME(JL)+ZB_RS(JL)
      ZRGT(JL)=ZRGT(JL)+ZA_RG(JL)*ZMAXTIME(JL)+ZB_RG(JL)
      ZCIT(JL)=ZCIT(JL) * MAX(0., -SIGN(1., -ZRIT(JL))) ! WHERE(ZRIT(:)==0.) ZCIT(:) = 0.
    ENDDO
    IF(KRR==7) ZRHT(:)=ZRHT(:)+ZA_RH(:)*ZMAXTIME(:)+ZB_RH(:)
    !
    !***       4.4 Mixing ratio change due to each process
    !
    IF(LBU_ENABLE) THEN
      ZTOT_RVHENI(:)= ZTOT_RVHENI(:) +ZRVHENI_MR(:)
      ZTOT_RCHONI(:)= ZTOT_RCHONI(:) +ZRCHONI(:) *ZMAXTIME(:)
      ZTOT_RRHONG(:)= ZTOT_RRHONG(:) +ZRRHONG_MR(:)
      ZTOT_RVDEPS(:)= ZTOT_RVDEPS(:) +ZRVDEPS(:) *ZMAXTIME(:)
      ZTOT_RIAGGS(:)= ZTOT_RIAGGS(:) +ZRIAGGS(:) *ZMAXTIME(:)
      ZTOT_RIAUTS(:)= ZTOT_RIAUTS(:) +ZRIAUTS(:) *ZMAXTIME(:)
      ZTOT_RVDEPG(:)= ZTOT_RVDEPG(:) +ZRVDEPG(:) *ZMAXTIME(:)
      ZTOT_RCAUTR(:)= ZTOT_RCAUTR(:) +ZRCAUTR(:) *ZMAXTIME(:)
      ZTOT_RCACCR(:)= ZTOT_RCACCR(:) +ZRCACCR(:) *ZMAXTIME(:)
      ZTOT_RREVAV(:)= ZTOT_RREVAV(:) +ZRREVAV(:) *ZMAXTIME(:)
      ZTOT_RCRIMSS(:)=ZTOT_RCRIMSS(:)+ZRCRIMSS(:)*ZMAXTIME(:)
      ZTOT_RCRIMSG(:)=ZTOT_RCRIMSG(:)+ZRCRIMSG(:)*ZMAXTIME(:)
      ZTOT_RSRIMCG(:)=ZTOT_RSRIMCG(:)+ZRSRIMCG(:)*ZMAXTIME(:)+ZRSRIMCG_MR(:)
      ZTOT_RRACCSS(:)=ZTOT_RRACCSS(:)+ZRRACCSS(:)*ZMAXTIME(:)
      ZTOT_RRACCSG(:)=ZTOT_RRACCSG(:)+ZRRACCSG(:)*ZMAXTIME(:)
      ZTOT_RSACCRG(:)=ZTOT_RSACCRG(:)+ZRSACCRG(:)*ZMAXTIME(:)
      ZTOT_RSMLTG(:)= ZTOT_RSMLTG(:) +ZRSMLTG(:) *ZMAXTIME(:)
      ZTOT_RCMLTSR(:)=ZTOT_RCMLTSR(:)+ZRCMLTSR(:) *ZMAXTIME(:)
      ZTOT_RICFRRG(:)=ZTOT_RICFRRG(:)+ZRICFRRG(:)*ZMAXTIME(:)
      ZTOT_RRCFRIG(:)=ZTOT_RRCFRIG(:)+ZRRCFRIG(:)*ZMAXTIME(:)
      ZTOT_RICFRR(:)= ZTOT_RICFRR(:) +ZRICFRR(:) *ZMAXTIME(:)
      ZTOT_RCWETG(:)= ZTOT_RCWETG(:) +ZRCWETG(:) *ZMAXTIME(:)
      ZTOT_RIWETG(:)= ZTOT_RIWETG(:) +ZRIWETG(:) *ZMAXTIME(:)
      ZTOT_RRWETG(:)= ZTOT_RRWETG(:) +ZRRWETG(:) *ZMAXTIME(:)
      ZTOT_RSWETG(:)= ZTOT_RSWETG(:) +ZRSWETG(:) *ZMAXTIME(:)
      ZTOT_RWETGH(:)= ZTOT_RWETGH(:) +ZRWETGH(:) *ZMAXTIME(:)+ZRWETGH_MR(:)
      ZTOT_RCDRYG(:)= ZTOT_RCDRYG(:) +ZRCDRYG(:) *ZMAXTIME(:)
      ZTOT_RIDRYG(:)= ZTOT_RIDRYG(:) +ZRIDRYG(:) *ZMAXTIME(:)
      ZTOT_RRDRYG(:)= ZTOT_RRDRYG(:) +ZRRDRYG(:) *ZMAXTIME(:)
      ZTOT_RSDRYG(:)= ZTOT_RSDRYG(:) +ZRSDRYG(:) *ZMAXTIME(:)
      ZTOT_RGMLTR(:)= ZTOT_RGMLTR(:) +ZRGMLTR(:) *ZMAXTIME(:)
      ZTOT_RCWETH(:)= ZTOT_RCWETH(:) +ZRCWETH(:) *ZMAXTIME(:)
      ZTOT_RIWETH(:)= ZTOT_RIWETH(:) +ZRIWETH(:) *ZMAXTIME(:)
      ZTOT_RSWETH(:)= ZTOT_RSWETH(:) +ZRSWETH(:) *ZMAXTIME(:)
      ZTOT_RGWETH(:)= ZTOT_RGWETH(:) +ZRGWETH(:) *ZMAXTIME(:)
      ZTOT_RRWETH(:)= ZTOT_RRWETH(:) +ZRRWETH(:) *ZMAXTIME(:)
      ZTOT_RCDRYH(:)= ZTOT_RCDRYH(:) +ZRCDRYH(:) *ZMAXTIME(:)
      ZTOT_RIDRYH(:)= ZTOT_RIDRYH(:) +ZRIDRYH(:) *ZMAXTIME(:)
      ZTOT_RSDRYH(:)= ZTOT_RSDRYH(:) +ZRSDRYH(:) *ZMAXTIME(:)
      ZTOT_RRDRYH(:)= ZTOT_RRDRYH(:) +ZRRDRYH(:) *ZMAXTIME(:)
      ZTOT_RGDRYH(:)= ZTOT_RGDRYH(:) +ZRGDRYH(:) *ZMAXTIME(:)
      ZTOT_RDRYHG(:)= ZTOT_RDRYHG(:) +ZRDRYHG(:) *ZMAXTIME(:)
      ZTOT_RHMLTR(:)= ZTOT_RHMLTR(:) +ZRHMLTR(:) *ZMAXTIME(:)
      ZTOT_RIMLTC(:)= ZTOT_RIMLTC(:) +ZRIMLTC_MR(:)
      ZTOT_RCBERI(:)= ZTOT_RCBERI(:) +ZRCBERI(:) *ZMAXTIME(:)
    ENDIF
    !
    !***       4.5 Next loop
    !
    LSOFT=.TRUE. ! We try to adjust tendencies (inner while loop)
    ZTIME(:)=ZTIME(:)+ZMAXTIME(:)
  ENDDO
ENDDO
!-------------------------------------------------------------------------------
!
!*       5.     UNPACKING DIAGNOSTICS
!               ---------------------
!
IF(IMICRO>0) THEN
  ZHLC_HCF3D(:,:,:)=0.
  ZHLC_LCF3D(:,:,:)=0.
  ZHLC_HRC3D(:,:,:)=0.
  ZHLC_LRC3D(:,:,:)=0.
  ZHLI_HCF3D(:,:,:)=0.
  ZHLI_LCF3D(:,:,:)=0.
  ZHLI_HRI3D(:,:,:)=0.
  ZHLI_LRI3D(:,:,:)=0.
  DO JL=1,IMICRO
    ZHLC_HCF3D(I1(JL), I2(JL), I3(JL)) = ZHLC_HCF(JL)
    ZHLC_LCF3D(I1(JL), I2(JL), I3(JL)) = ZHLC_LCF(JL)
    ZHLC_HRC3D(I1(JL), I2(JL), I3(JL)) = ZHLC_HRC(JL)
    ZHLC_LRC3D(I1(JL), I2(JL), I3(JL)) = ZHLC_LRC(JL)
    ZHLI_LCF3D(I1(JL), I2(JL), I3(JL)) = ZHLI_LCF(JL)
    ZHLI_HCF3D(I1(JL), I2(JL), I3(JL)) = ZHLI_HCF(JL)
    ZHLI_HRI3D(I1(JL), I2(JL), I3(JL)) = ZHLI_HRI(JL)
    ZHLI_LRI3D(I1(JL), I2(JL), I3(JL)) = ZHLI_LRI(JL)
    PCIT(I1(JL), I2(JL), I3(JL)) = ZCIT(JL)
  END DO
ELSE
  PRAINFR(:,:,:)=0.
  ZHLC_HCF3D(:,:,:)=0.
  ZHLC_LCF3D(:,:,:)=0.
  ZHLC_HRC3D(:,:,:)=0.
  ZHLC_LRC3D(:,:,:)=0.
  ZHLI_HCF3D(:,:,:)=0.
  ZHLI_LCF3D(:,:,:)=0.
  ZHLI_HRI3D(:,:,:)=0.
  ZHLI_LRI3D(:,:,:)=0.
  PCIT(:,:,:) = 0.
ENDIF
IF(OWARM) THEN
  PEVAP3D(:,:,:) = 0.
  DO JL=1,IMICRO
    PEVAP3D(I1(JL), I2(JL), I3(JL)) = ZRREVAV(JL)
  END DO
ENDIF
!
!
!*       6.     COMPUTES THE SLOW COLD PROCESS SOURCES OUTSIDE OF ODMICRO POINTS
!               ----------------------------------------------------------------
!
CALL ICE4_NUCLEATION_WRAPPER(KIT, KJT, KKT, .NOT. ODMICRO, &
                             PTHT, PPABST, PRHODREF, PEXN, ZZ_LSFACT/PEXN, ZT, &
                             PRVT, &
                             PCIT, ZZ_RVHENI_MR)
DO JK = 1, KKT
  DO JJ = 1, KJT
    DO JI = 1, KIT
      ZZ_LSFACT(JI,JJ,JK)=ZZ_LSFACT(JI,JJ,JK)/PEXNREF(JI,JJ,JK)
      ZZ_LVFACT(JI,JJ,JK)=ZZ_LVFACT(JI,JJ,JK)/PEXNREF(JI,JJ,JK)
      ZZ_RVHENI(JI,JJ,JK) = MIN(PRVS(JI,JJ,JK), ZZ_RVHENI_MR(JI,JJ,JK)/PTSTEP)
      PRIS(JI,JJ,JK)=PRIS(JI,JJ,JK)+ZZ_RVHENI(JI,JJ,JK)
      PRVS(JI,JJ,JK)=PRVS(JI,JJ,JK)-ZZ_RVHENI(JI,JJ,JK)
      PTHS(JI,JJ,JK)=PTHS(JI,JJ,JK) + ZZ_RVHENI(JI,JJ,JK)*ZZ_LSFACT(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
!
if ( lbu_enable ) then
  !Note: there is an other contribution for HENU later
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HENU', pths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'HENU', prvs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'HENU', zz_rvheni(:, :, :) * prhodj(:, :, :) )
end if
!-------------------------------------------------------------------------------
!
!*       7.     UNPACKING AND TOTAL TENDENCIES
!               ------------------------------
!
!
!***     7.1    total tendencies limited by available species
!
! ZW_??S variables will contain the new S variables values
!
IF(GEXT_TEND) THEN
  !Z..T variables contain the exeternal tendency, we substract it
  DO JL=1, IMICRO
    ZRVT(JL) = ZRVT(JL) - ZEXT_RV(JL) * PTSTEP
    ZRCT(JL) = ZRCT(JL) - ZEXT_RC(JL) * PTSTEP
    ZRRT(JL) = ZRRT(JL) - ZEXT_RR(JL) * PTSTEP
    ZRIT(JL) = ZRIT(JL) - ZEXT_RI(JL) * PTSTEP
    ZRST(JL) = ZRST(JL) - ZEXT_RS(JL) * PTSTEP
    ZRGT(JL) = ZRGT(JL) - ZEXT_RG(JL) * PTSTEP
    ZTHT(JL) = ZTHT(JL) - ZEXT_TH(JL) * PTSTEP
  ENDDO
  IF (KRR==7) ZRHT(:) = ZRHT(:) - ZEXT_RH(:) * PTSTEP
ENDIF
!Tendencies computed from difference between old state and new state (can be negative)
  ZW_RVS(:,:,:) = 0.
  ZW_RCS(:,:,:) = 0.
  ZW_RRS(:,:,:) = 0.
  ZW_RIS(:,:,:) = 0.
  ZW_RSS(:,:,:) = 0.
  ZW_RGS(:,:,:) = 0.
  ZW_RHS(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW_RVS(I1(JL), I2(JL), I3(JL)) = ( ZRVT(JL) - PRVT(I1(JL), I2(JL), I3(JL)) ) * ZINV_TSTEP
    ZW_RCS(I1(JL), I2(JL), I3(JL)) = ( ZRCT(JL) - PRCT(I1(JL), I2(JL), I3(JL)) ) * ZINV_TSTEP
    ZW_RRS(I1(JL), I2(JL), I3(JL)) = ( ZRRT(JL) - PRRT(I1(JL), I2(JL), I3(JL)) ) * ZINV_TSTEP
    ZW_RIS(I1(JL), I2(JL), I3(JL)) = ( ZRIT(JL) - PRIT(I1(JL), I2(JL), I3(JL)) ) * ZINV_TSTEP
    ZW_RSS(I1(JL), I2(JL), I3(JL)) = ( ZRST(JL) - PRST(I1(JL), I2(JL), I3(JL)) ) * ZINV_TSTEP
    ZW_RGS(I1(JL), I2(JL), I3(JL)) = ( ZRGT(JL) - PRGT(I1(JL), I2(JL), I3(JL)) ) * ZINV_TSTEP
  END DO
  IF(KRR==7) THEN
  DO JL=1,IMICRO
    ZW_RHS(I1(JL), I2(JL), I3(JL)) = ( ZRHT(JL) - PRHT(I1(JL), I2(JL), I3(JL)) ) * ZINV_TSTEP
  END DO
END IF
ZW_THS(:,:,:) = (ZW_RCS(:,:,:)+ZW_RRS(:,:,:)                            )*ZZ_LVFACT(:,:,:) + &
              & (ZW_RIS(:,:,:)+ZW_RSS(:,:,:)+ZW_RGS(:,:,:)+ZW_RHS(:,:,:))*ZZ_LSFACT(:,:,:)
!We apply these tendencies to the S variables
ZW_RVS(:,:,:) = PRVS(:,:,:) + ZW_RVS(:,:,:)
ZW_RCS(:,:,:) = PRCS(:,:,:) + ZW_RCS(:,:,:)
ZW_RRS(:,:,:) = PRRS(:,:,:) + ZW_RRS(:,:,:)
ZW_RIS(:,:,:) = PRIS(:,:,:) + ZW_RIS(:,:,:)
ZW_RSS(:,:,:) = PRSS(:,:,:) + ZW_RSS(:,:,:)
ZW_RGS(:,:,:) = PRGS(:,:,:) + ZW_RGS(:,:,:)
IF(KRR==7) ZW_RHS(:,:,:) = PRHS(:,:,:) + ZW_RHS(:,:,:)
ZW_THS(:,:,:) = PTHS(:,:,:) + ZW_THS(:,:,:)

if ( lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'CORR', zw_ths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'CORR', zw_rvs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'CORR', zw_rcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'CORR', zw_rrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'CORR', zw_ris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'CORR', zw_rss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'CORR', zw_rgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'CORR', zw_rhs(:, :, :) * prhodj(:, :, :) )
end if

!We correct negativities with conservation
CALL CORRECT_NEGATIVITIES(KIT, KJT, KKT, KRR, ZW_RVS, ZW_RCS, ZW_RRS, &
                         &ZW_RIS, ZW_RSS, ZW_RGS, &
                         &ZW_THS, ZZ_LVFACT, ZZ_LSFACT, ZW_RHS)

if ( lbu_enable ) then
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'CORR', zw_ths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'CORR', zw_rvs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'CORR', zw_rcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'CORR', zw_rrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'CORR', zw_ris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'CORR', zw_rss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'CORR', zw_rgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'CORR', zw_rhs(:, :, :) * prhodj(:, :, :) )
end if
!
!***     7.2    LBU_ENABLE case
!
IF(LBU_ENABLE) THEN
  allocate( zw1( size( pexnref, 1 ), size( pexnref, 2 ), size( pexnref, 3 ) ) )
  allocate( zw2( size( pexnref, 1 ), size( pexnref, 2 ), size( pexnref, 3 ) ) )
  allocate( zw3( size( pexnref, 1 ), size( pexnref, 2 ), size( pexnref, 3 ) ) )
  allocate( zw4( size( pexnref, 1 ), size( pexnref, 2 ), size( pexnref, 3 ) ) )
  if ( krr == 7 ) then
    allocate( zw5( size( pexnref, 1 ), size( pexnref, 2 ), size( pexnref, 3 ) ) )
    allocate( zw6( size( pexnref, 1 ), size( pexnref, 2 ), size( pexnref, 3 ) ) )
  end if

  if ( lbudget_th ) then
    allocate( zz_diff( size( zz_lsfact, 1 ), size( zz_lsfact, 2 ), size( zz_lsfact, 3 ) ) )
    zz_diff(:, :, :) = zz_lsfact(:, :, :) - zz_lvfact(:, :, :)
  end if

  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RVHENI(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'HENU',  zw(:, :, :) * zz_lsfact(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'HENU', -zw(:, :, :)                      * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'HENU',  zw(:, :, :)                      * prhodj(:, :, :) )

  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RCHONI(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'HON',  zw(:, :, :) * zz_diff(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'HON', -zw(:, :, :)                    * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'HON',  zw(:, :, :)                    * prhodj(:, :, :) )

  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RRHONG(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'SFR',  zw(:, :, :) * zz_diff(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'SFR', -zw(:, :, :)                    * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'SFR',  zw(:, :, :)                    * prhodj(:, :, :) )

  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RVDEPS(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'DEPS',  zw(:, :, :) * zz_lsfact(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'DEPS', -zw(:, :, :)                      * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'DEPS',  zw(:, :, :)                      * prhodj(:, :, :) )

  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RIAGGS(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'AGGS', -zw(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'AGGS',  zw(:, :, :) * prhodj(:, :, :) )

  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RIAUTS(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'AUTS', -zw(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'AUTS',  zw(:, :, :) * prhodj(:, :, :) )

  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RVDEPG(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'DEPG',  zw(:, :, :) * zz_lsfact(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'DEPG', -zw(:, :, :)                      * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'DEPG',  zw(:, :, :)                      * prhodj(:, :, :) )

  IF(OWARM) THEN
    ZW(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RCAUTR(JL) * ZINV_TSTEP
    END DO
    if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'AUTO', -zw(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'AUTO',  zw(:, :, :) * prhodj(:, :, :) )

    ZW(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RCACCR(JL) * ZINV_TSTEP
    END DO
    if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'ACCR', -zw(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'ACCR',  zw(:, :, :) * prhodj(:, :, :) )

    ZW(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RREVAV(JL) * ZINV_TSTEP
    END DO
    if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'REVA', -zw(:, :, :) * zz_lvfact(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'REVA',  zw(:, :, :)                      * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'REVA', -zw(:, :, :)                      * prhodj(:, :, :) )
  ENDIF

  ZW1(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW1(I1(JL), I2(JL), I3(JL)) = ZTOT_RCRIMSS(JL) * ZINV_TSTEP
  END DO
  ZW2(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW2(I1(JL), I2(JL), I3(JL)) = ZTOT_RCRIMSG(JL) * ZINV_TSTEP
  END DO
  ZW3(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW3(I1(JL), I2(JL), I3(JL)) = ZTOT_RSRIMCG(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) &
    call Budget_store_add( tbudgets(NBUDGET_TH), 'RIM', (  zw1(:, :, :) + zw2(:, :, :) ) * zz_diff(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'RIM', ( -zw1(:, :, :) - zw2(:, :, :) )    * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'RIM', (  zw1(:, :, :) - zw3(:, :, :) )    * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'RIM', (  zw2(:, :, :) + zw3(:, :, :) )    * prhodj(:, :, :) )

  ZW1(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW1(I1(JL), I2(JL), I3(JL)) = ZTOT_RRACCSS(JL) * ZINV_TSTEP
  END DO
  ZW2(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW2(I1(JL), I2(JL), I3(JL)) = ZTOT_RRACCSG(JL) * ZINV_TSTEP
  END DO
  ZW3(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW3(I1(JL), I2(JL), I3(JL)) = ZTOT_RSACCRG(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) &
    call Budget_store_add( tbudgets(NBUDGET_TH), 'ACC', (  zw1(:, :, :) + zw2(:, :, :) ) * zz_diff(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'ACC', ( -zw1(:, :, :) - zw2(:, :, :) )    * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'ACC', (  zw1(:, :, :) - zw3(:, :, :) )    * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'ACC', (  zw2(:, :, :) + zw3(:, :, :) )    * prhodj(:, :, :) )

  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RSMLTG(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'CMEL', -zw(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'CMEL',  zw(:, :, :) * prhodj(:, :, :) )
  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RCMLTSR(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'CMEL', -zw(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'CMEL',  zw(:, :, :) * prhodj(:, :, :) )

  ZW1(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW1(I1(JL), I2(JL), I3(JL)) = ZTOT_RICFRRG(JL) * ZINV_TSTEP
  END DO
  ZW2(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW2(I1(JL), I2(JL), I3(JL)) = ZTOT_RRCFRIG(JL) * ZINV_TSTEP
  END DO
  ZW3(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW3(I1(JL), I2(JL), I3(JL)) = ZTOT_RICFRR(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) &
    call Budget_store_add( tbudgets(NBUDGET_TH), 'CFRZ', zw2(:, :, :) * zz_diff(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'CFRZ', ( -zw2(:, :, :) + zw3(:, :, :) )    * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'CFRZ', ( -zw1(:, :, :) - zw3(:, :, :) )    * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'CFRZ', (  zw1(:, :, :) + zw2(:, :, :) )    * prhodj(:, :, :) )

  ZW1(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW1(I1(JL), I2(JL), I3(JL)) = ZTOT_RCWETG(JL) * ZINV_TSTEP
  END DO
  ZW2(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW2(I1(JL), I2(JL), I3(JL)) = ZTOT_RRWETG(JL) * ZINV_TSTEP
  END DO
  ZW3(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW3(I1(JL), I2(JL), I3(JL)) = ZTOT_RIWETG(JL) * ZINV_TSTEP
  END DO
  ZW4(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW4(I1(JL), I2(JL), I3(JL)) = ZTOT_RSWETG(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) &
    call Budget_store_add( tbudgets(NBUDGET_TH), 'WETG', (  zw1(:, :, :) + zw2(:, :, :) ) * zz_diff(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'WETG',   -zw1(:, :, :)                     * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'WETG',   -zw2(:, :, :)                     * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'WETG',   -zw3(:, :, :)                     * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'WETG',   -zw4(:, :, :)                     * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'WETG', (  zw1(:, :, :) + zw2(:, :, :)   &
                                                                          + zw3(:, :, :) + zw4(:, :, :) ) &
                                                                            * prhodj(:, :, :) )

  IF(KRR==7) THEN
    ZW(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RWETGH(JL) * ZINV_TSTEP
    END DO
    if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'GHCV', -zw(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_add( tbudgets(NBUDGET_RH), 'GHCV',  zw(:, :, :) * prhodj(:, :, :) )
  END IF

  ZW1(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW1(I1(JL), I2(JL), I3(JL)) = ZTOT_RCDRYG(JL) * ZINV_TSTEP
  END DO
  ZW2(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW2(I1(JL), I2(JL), I3(JL)) = ZTOT_RRDRYG(JL) * ZINV_TSTEP
  END DO
  ZW3(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW3(I1(JL), I2(JL), I3(JL)) = ZTOT_RIDRYG(JL) * ZINV_TSTEP
  END DO
  ZW4(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW4(I1(JL), I2(JL), I3(JL)) = ZTOT_RSDRYG(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) &
    call Budget_store_add( tbudgets(NBUDGET_TH), 'DRYG', (  zw1(:, :, :) + zw2(:, :, :) ) * zz_diff(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'DRYG',   -zw1(:, :, :)                     * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'DRYG',   -zw2(:, :, :)                     * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'DRYG',   -zw3(:, :, :)                     * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'DRYG',   -zw4(:, :, :)                     * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'DRYG', (  zw1(:, :, :) + zw2(:, :, :)   &
                                                                          + zw3(:, :, :) + zw4(:, :, :) ) &
                                                                            * prhodj(:, :, :) )

  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RGMLTR(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'GMLT', -zw(:, :, :) * zz_diff(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'GMLT',  zw(:, :, :)                    * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'GMLT', -zw(:, :, :)                    * prhodj(:, :, :) )

  IF(KRR==7) THEN
    ZW1(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW1(I1(JL), I2(JL), I3(JL)) = ZTOT_RCWETH(JL) * ZINV_TSTEP
    END DO
    ZW2(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW2(I1(JL), I2(JL), I3(JL)) = ZTOT_RRWETH(JL) * ZINV_TSTEP
    END DO
    ZW3(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW3(I1(JL), I2(JL), I3(JL)) = ZTOT_RIWETH(JL) * ZINV_TSTEP
    END DO
    ZW4(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW4(I1(JL), I2(JL), I3(JL)) = ZTOT_RSWETH(JL) * ZINV_TSTEP
    END DO
    ZW5(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW5(I1(JL), I2(JL), I3(JL)) = ZTOT_RGWETH(JL) * ZINV_TSTEP
    END DO
    if ( lbudget_th ) &
      call Budget_store_add( tbudgets(NBUDGET_TH), 'WETH', (  zw1(:, :, :) + zw2(:, :, :) ) * zz_diff(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'WETH',   -zw1(:, :, :)                     * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'WETH',   -zw2(:, :, :)                     * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'WETH',   -zw3(:, :, :)                     * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'WETH',   -zw4(:, :, :)                     * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'WETH',   -zw5(:, :, :)                     * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_add( tbudgets(NBUDGET_RH), 'WETH', (  zw1(:, :, :) + zw2(:, :, :) + zw3(:, :, :) &
                                                                              + zw4(:, :, :) + zw5(:, :, : ) )           &
                                                                              * prhodj(:, :, :) )

    ZW(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RGWETH(JL) * ZINV_TSTEP
    END DO
    if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'HGCV', -zw(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_add( tbudgets(NBUDGET_RH), 'HGCV',  zw(:, :, :) * prhodj(:, :, :) )

    ZW1(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW1(I1(JL), I2(JL), I3(JL)) = ZTOT_RCDRYH(JL) * ZINV_TSTEP
    END DO
    ZW2(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW2(I1(JL), I2(JL), I3(JL)) = ZTOT_RRDRYH(JL) * ZINV_TSTEP
    END DO
    ZW3(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW3(I1(JL), I2(JL), I3(JL)) = ZTOT_RIDRYH(JL) * ZINV_TSTEP
    END DO
    ZW4(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW4(I1(JL), I2(JL), I3(JL)) = ZTOT_RSDRYH(JL) * ZINV_TSTEP
    END DO
    ZW5(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW5(I1(JL), I2(JL), I3(JL)) = ZTOT_RGDRYH(JL) * ZINV_TSTEP
    END DO
    ZW6(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW6(I1(JL), I2(JL), I3(JL)) = ZTOT_RDRYHG(JL) * ZINV_TSTEP
    END DO
    if ( lbudget_th ) &
      call Budget_store_add( tbudgets(NBUDGET_TH), 'DRYH', (  zw1(:, :, :) + zw2(:, :, :) ) * zz_diff(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'DRYH',   -zw1(:, :, :)                     * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'DRYH',   -zw2(:, :, :)                     * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'DRYH',   -zw3(:, :, :)                     * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_add( tbudgets(NBUDGET_RS), 'DRYH',   -zw4(:, :, :)                     * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'DRYH', ( -zw5(:, :, :) + zw6(:, :, : ) )   * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_add( tbudgets(NBUDGET_RH), 'DRYH', (  zw1(:, :, :) + zw2(:, :, :) + zw3(:, :, :)   &
                                                                            + zw4(:, :, :) + zw5(:, :, : )- zw6(:, :, :) ) &
                                                                              * prhodj(:, :, :) )

    ZW(:,:,:) = 0.
    DO JL=1,IMICRO
      ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RHMLTR(JL) * ZINV_TSTEP
    END DO
    if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'HMLT', -zw(:, :, :) * zz_diff(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_add( tbudgets(NBUDGET_RR), 'HMLT',  zw(:, :, :)                    * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_add( tbudgets(NBUDGET_RH), 'HMLT', -zw(:, :, :)                    * prhodj(:, :, :) )
  ENDIF

  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RIMLTC(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'IMLT', -zw(:, :, :) * zz_diff(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'IMLT',  zw(:, :, :)                    * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'IMLT', -zw(:, :, :)                    * prhodj(:, :, :) )

  ZW(:,:,:) = 0.
  DO JL=1,IMICRO
    ZW(I1(JL), I2(JL), I3(JL)) = ZTOT_RCBERI(JL) * ZINV_TSTEP
  END DO
  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'BERFI',  zw(:, :, :) * zz_diff(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'BERFI', -zw(:, :, :)                    * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'BERFI',  zw(:, :, :)                    * prhodj(:, :, :) )

  deallocate( zw1, zw2, zw3, zw4 )
  if ( krr == 7 ) deallocate( zw5, zw6 )
  if ( lbudget_th ) deallocate( zz_diff )
ENDIF
!
!***     7.3    Final tendencies
!
DO JK = 1, KKT
  DO JJ = 1, KJT
    DO JI = 1, KIT
      PRVS(JI,JJ,JK) = ZW_RVS(JI,JJ,JK)
      PRCS(JI,JJ,JK) = ZW_RCS(JI,JJ,JK)
      PRRS(JI,JJ,JK) = ZW_RRS(JI,JJ,JK)
      PRIS(JI,JJ,JK) = ZW_RIS(JI,JJ,JK)
      PRSS(JI,JJ,JK) = ZW_RSS(JI,JJ,JK)
      PRGS(JI,JJ,JK) = ZW_RGS(JI,JJ,JK)
      PTHS(JI,JJ,JK) = ZW_THS(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
IF (KRR==7) PRHS(:,:,:) = ZW_RHS(:,:,:)
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
  if ( lbudget_rc .and. osedic ) call Budget_store_init( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr )              call Budget_store_init( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri )              call Budget_store_init( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs )              call Budget_store_init( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg )              call Budget_store_init( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh )              call Budget_store_init( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )

  !Init only if not osedic (to prevent crash with double init)
  !Remark: the 2 source terms SEDI and DEPO could be mixed and stored in the same source term (SEDI)
  !        if osedic=T and ldeposc=T (a warning is printed in ini_budget in that case)
  if ( lbudget_rc .and. ldeposc .and. .not.osedic ) &
    call Budget_store_init( tbudgets(NBUDGET_RC), 'DEPO', prcs(:, :, :) * prhodj(:, :, :) )

  IF(HSEDIM=='STAT') THEN
    !SR: It *seems* that we must have two separate calls for ifort
    IF(KRR==7) THEN
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCS*PTSTEP, PRRS, PRRS*PTSTEP, PRIS, PRIS*PTSTEP,&
                                  &PRSS, PRSS*PTSTEP, PRGS, PRGS*PTSTEP,&
                                  &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA=PSEA, PTOWN=PTOWN, &
                                  &PINPRH=PINPRH, PRHT=PRHS*PTSTEP, PRHS=PRHS, PFPR=PFPR)
    ELSE
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ,&
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCS*PTSTEP, PRRS, PRRS*PTSTEP, PRIS, PRIS*PTSTEP,&
                                  &PRSS, PRSS*PTSTEP, PRGS, PRGS*PTSTEP,&
                                  &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA=PSEA, PTOWN=PTOWN, &
                                  &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !No negativity correction here as we apply sedimentation on PR.S*PTSTEP variables
  ELSEIF(HSEDIM=='SPLI') THEN
    !SR: It *seems* that we must have two separate calls for ifort
    IF(KRR==7) THEN
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA=PSEA, PTOWN=PTOWN, &
                                   &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR)
    ELSE
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, IKTB, IKTE, KKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA=PSEA, PTOWN=PTOWN, &
                                   &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !We correct negativities with conservation
    !SPLI algorith uses a time-splitting. Inside the loop a temporary m.r. is used.
    !   It is initialized with the m.r. at T and is modified by two tendencies:
    !   sedimentation tendency and an external tendency which represents all other
    !   processes (mainly advection and microphysical processes). If both tendencies
    !   are negative, sedimentation can remove a specie at a given sub-timestep. From
    !   this point sedimentation stops for the remaining sub-timesteps but the other tendency
    !   will be still active and will lead to negative values.
    !   We could prevent the algorithm to not consume too much a specie, instead we apply
    !   a correction here.
    CALL CORRECT_NEGATIVITIES(KIT, KJT, KKT, KRR, PRVS, PRCS, PRRS, &
                             &PRIS, PRSS, PRGS, &
                             &PTHS, ZZ_LVFACT, ZZ_LSFACT, PRHS)
  ELSE
    call Print_msg( NVERB_FATAL, 'GEN', 'RAIN_ICE_RED', 'no sedimentation scheme for HSEDIM='//HSEDIM )
  END IF
  !
  !*       8.2     budget storage
  !
  if ( lbudget_rc .and. osedic ) call Budget_store_end( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr )              call Budget_store_end( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri )              call Budget_store_end( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs )              call Budget_store_end( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg )              call Budget_store_end( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh )              call Budget_store_end( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )

  !If osedic=T and ldeposc=T, DEPO is in fact mixed and stored with the SEDI source term
  !(a warning is printed in ini_budget in that case)
  if ( lbudget_rc .and. ldeposc .and. .not.osedic) &
    call Budget_store_end( tbudgets(NBUDGET_RC), 'DEPO', prcs(:, :, :) * prhodj(:, :, :) )

  !sedimentation of rain fraction
  IF (PRESENT(PRHS)) THEN
    CALL ICE4_RAINFR_VERT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, KKT, KKL, PRAINFR, PRRS(:,:,:)*PTSTEP, &
                       &PRSS(:,:,:)*PTSTEP, PRGS(:,:,:)*PTSTEP, PRHS(:,:,:)*PTSTEP)
  ELSE
    CALL ICE4_RAINFR_VERT(IIB, IIE, KIT, IJB, IJE, KJT, IKB, IKE, KKT, KKL, PRAINFR, PRRS(:,:,:)*PTSTEP, &
                       &PRSS(:,:,:)*PTSTEP, PRGS(:,:,:)*PTSTEP)
  ENDIF
ENDIF
!
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
  REAL, DIMENSION(KIT, KJT, KKT) :: ZW
  INTEGER :: JI, JJ, JK
  !
  !
  !We correct negativities with conservation
  ! 1) deal with negative values for mixing ratio, except for vapor
  DO JK = 1, KKT
    DO JJ = 1, KJT
      DO JI = 1, KIT
        ZW(JI,JJ,JK) =PRC(JI,JJ,JK)-MAX(PRC(JI,JJ,JK), 0.)
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW(JI,JJ,JK)
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW(JI,JJ,JK)*PLVFACT(JI,JJ,JK)
        PRC(JI,JJ,JK)=PRC(JI,JJ,JK)-ZW(JI,JJ,JK)

        ZW(JI,JJ,JK) =PRR(JI,JJ,JK)-MAX(PRR(JI,JJ,JK), 0.)
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW(JI,JJ,JK)
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW(JI,JJ,JK)*PLVFACT(JI,JJ,JK)
        PRR(JI,JJ,JK)=PRR(JI,JJ,JK)-ZW(JI,JJ,JK)

        ZW(JI,JJ,JK) =PRI(JI,JJ,JK)-MAX(PRI(JI,JJ,JK), 0.)
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW(JI,JJ,JK)
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW(JI,JJ,JK)*PLSFACT(JI,JJ,JK)
        PRI(JI,JJ,JK)=PRI(JI,JJ,JK)-ZW(JI,JJ,JK)

        ZW(JI,JJ,JK) =PRS(JI,JJ,JK)-MAX(PRS(JI,JJ,JK), 0.)
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW(JI,JJ,JK)
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW(JI,JJ,JK)*PLSFACT(JI,JJ,JK)
        PRS(JI,JJ,JK)=PRS(JI,JJ,JK)-ZW(JI,JJ,JK)

        ZW(JI,JJ,JK) =PRG(JI,JJ,JK)-MAX(PRG(JI,JJ,JK), 0.)
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW(JI,JJ,JK)
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW(JI,JJ,JK)*PLSFACT(JI,JJ,JK)
        PRG(JI,JJ,JK)=PRG(JI,JJ,JK)-ZW(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO

  IF(KRR==7) THEN
    DO JK = 1, KKT
      DO JJ = 1, KJT
        DO JI = 1, KIT
          ZW(JI,JJ,JK) =PRH(JI,JJ,JK)-MAX(PRH(JI,JJ,JK), 0.)
          PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW(JI,JJ,JK)
          PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW(JI,JJ,JK)*PLSFACT(JI,JJ,JK)
          PRH(JI,JJ,JK)=PRH(JI,JJ,JK)-ZW(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! 2) deal with negative vapor mixing ratio

  DO JK = 1, KKT
    DO JJ = 1, KJT
      DO JI = 1, KIT
        ! for rc and ri, we keep ice fraction constant
        ZW(JI,JJ,JK)=MIN(1., MAX(XRTMIN(1)-PRV(JI,JJ,JK), 0.) / &
                            &MAX(PRC(JI,JJ,JK)+PRI(JI,JJ,JK), 1.E-20)) ! Proportion of rc+ri to convert into rv
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW(JI,JJ,JK)* &
                     &(PRC(JI,JJ,JK)*PLVFACT(JI,JJ,JK)+PRI(JI,JJ,JK)*PLSFACT(JI,JJ,JK))
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW(JI,JJ,JK)*(PRC(JI,JJ,JK)+PRI(JI,JJ,JK))
        PRC(JI,JJ,JK)=(1.-ZW(JI,JJ,JK))*PRC(JI,JJ,JK)
        PRI(JI,JJ,JK)=(1.-ZW(JI,JJ,JK))*PRI(JI,JJ,JK)

        ZW(JI,JJ,JK)=MIN(MAX(PRR(JI,JJ,JK), 0.), &
                        &MAX(XRTMIN(1)-PRV(JI,JJ,JK), 0.)) ! Quantity of rr to convert into rv
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW(JI,JJ,JK)
        PRR(JI,JJ,JK)=PRR(JI,JJ,JK)-ZW(JI,JJ,JK)
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW(JI,JJ,JK)*PLVFACT(JI,JJ,JK)

        ZW(JI,JJ,JK)=MIN(MAX(PRS(JI,JJ,JK), 0.), &
                        &MAX(XRTMIN(1)-PRV(JI,JJ,JK), 0.)) ! Quantity of rs to convert into rv
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW(JI,JJ,JK)
        PRS(JI,JJ,JK)=PRS(JI,JJ,JK)-ZW(JI,JJ,JK)
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW(JI,JJ,JK)*PLSFACT(JI,JJ,JK)

        ZW(JI,JJ,JK)=MIN(MAX(PRG(JI,JJ,JK), 0.), &
                        &MAX(XRTMIN(1)-PRV(JI,JJ,JK), 0.)) ! Quantity of rg to convert into rv
        PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW(JI,JJ,JK)
        PRG(JI,JJ,JK)=PRG(JI,JJ,JK)-ZW(JI,JJ,JK)
        PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW(JI,JJ,JK)*PLSFACT(JI,JJ,JK)
      ENDDO
    ENDDO
  ENDDO

  IF(KRR==7) THEN
    DO JK = 1, KKT
      DO JJ = 1, KJT
        DO JI = 1, KIT
          ZW(JI,JJ,JK)=MIN(MAX(PRH(JI,JJ,JK), 0.), &
                          &MAX(XRTMIN(1)-PRV(JI,JJ,JK), 0.)) ! Quantity of rh to convert into rv
          PRV(JI,JJ,JK)=PRV(JI,JJ,JK)+ZW(JI,JJ,JK)
          PRH(JI,JJ,JK)=PRH(JI,JJ,JK)-ZW(JI,JJ,JK)
          PTH(JI,JJ,JK)=PTH(JI,JJ,JK)-ZW(JI,JJ,JK)*PLSFACT(JI,JJ,JK)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  !
  !
  END SUBROUTINE CORRECT_NEGATIVITIES

!
END SUBROUTINE RAIN_ICE_RED
