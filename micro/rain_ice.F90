!MNH_LIC Copyright 1995-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
      SUBROUTINE RAIN_ICE ( D, CST, PARAMI, ICEP, ICED, ELECP, ELECD, BUCONF,     &
                            OELEC, OSEDIM_BEARD, PTHVREFZIKB,                     &
                            PTSTEP, KRR, PEXN,                                    &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                            PICLDFR, PSSIO, PSSIU, PIFR,                 &
                            PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,               &
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                   &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                            PINPRC, PINPRR, PEVAP3D,                              &
                            PINPRS, PINPRG, PINDEP, PRAINFR, PSIGS,               &
                            TBUDGETS, KBUDGETS,                                   &
                            PQPIT, PQCT, PQRT, PQIT, PQST, PQGT, PQNIT,           &
                            PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS,           &                            
                            PEFIELDW, PLATHAM_IAGGS,                              &
                            PSEA, PTOWN, PCONC3D,                                 &
                            PRHT, PRHS, PINPRH, PFPR, PQHT, PQHS                  )
!     #############################################################################
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
!!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!!  P. Wautelet 29/05/2019: remove PACK/UNPACK intrinsics (to get more performance and better OpenACC support)
!!  P. Wautelet 17/01/2020: move Quicksort to tools.f90
!!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!!  P. Wautelet 25/02/2020: bugfix: add missing budget: WETH_BU_RRG
!!     R. El Khatib 24-Aug-2021 Optimizations
!  J. Wurtz       03/2022: New snow characteristics with LSNOW_T
!  C. Barthe      03/2023: Add call to cloud electrification
!  C. Barthe      06/2023: Add retroaction of electric field on IAGGS
!  C. Barthe      07/2023: use new data structures for electricity
!! S. Riette Sept 23: e from ice4_tendencies
!-----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_BUDGET,         ONLY: TBUDGETDATA_PTR, TBUDGETCONF_t, NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, &
                               NBUDGET_RI, NBUDGET_RR, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_ELEC_PARAM,     ONLY: ELEC_PARAM_t
USE MODD_ELEC_DESCR,     ONLY: ELEC_DESCR_t
USE MODD_FIELDS_ADDRESS

USE MODE_MSG,                ONLY: PRINT_MSG, NVERB_FATAL

USE MODE_ICE4_RAINFR_VERT,   ONLY: ICE4_RAINFR_VERT
USE MODE_ICE4_COMPUTE_PDF,   ONLY: ICE4_COMPUTE_PDF
USE MODE_ICE4_SEDIMENTATION, ONLY: ICE4_SEDIMENTATION
USE MODE_ICE4_PACK, ONLY: ICE4_PACK
USE MODE_ICE4_CORRECT_NEGATIVITIES, ONLY: ICE4_CORRECT_NEGATIVITIES
!
USE MODE_ELEC_TENDENCIES, ONLY : ELEC_TENDENCIES
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
TYPE(ELEC_PARAM_t),       INTENT(IN)    :: ELECP   ! electrical parameters
TYPE(ELEC_DESCR_t),       INTENT(IN)    :: ELECD   ! electrical descriptive csts
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
LOGICAL,                  INTENT(IN)    :: OELEC  ! Switch for cloud electricity
LOGICAL,                  INTENT(IN)    :: OSEDIM_BEARD  ! Switch for effect of electrical forces on sedim.
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
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
! OCND2
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PICLDFR ! Ice cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSSIO   ! Super-saturation with respect to ice in the supersaturated fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSSIU   ! Sub-saturation with respect to ice in the  subsaturated fraction 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PIFR    ! Ratio cloud ice moist part to dry part 

REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLC_HRC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLC_HCF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLI_HRI
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLI_HCF
!
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
REAL, DIMENSION(D%NIJT), INTENT(OUT) :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) :: PRAINFR !Precipitation fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
REAL, INTENT(IN)                :: PTHVREFZIKB ! Reference thv at IKB for electricity
!
TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
! scalar variables for cloud electricity
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(IN)    :: PQPIT  ! Positive ion  -
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQCT   ! Cloud droplet | 
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQRT   ! Rain          | electric
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQIT   ! Ice crystals  |  charge 
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQST   ! Snow          |   at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQGT   ! Graupel       |
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(IN)    :: PQNIT  ! Negative ion  -
!
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQPIS  ! Positive ion  -
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQCS   ! Cloud droplet | 
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQRS   ! Rain          | electric
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQIS   ! Ice crystals  |  charge 
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQSS   ! Snow          |  source
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQGS   ! Graupel       |
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQNIS  ! Negative ion  -
!
REAL, DIMENSION(MERGE(D%NIJT,0,OSEDIM_BEARD),MERGE(D%NKT,0,OSEDIM_BEARD)), OPTIONAL, INTENT(IN) :: PEFIELDW ! vertical electric field
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(IN)    :: PLATHAM_IAGGS ! E Function to simulate
                                                                                            ! enhancement of IAGGS
!
! optional variables
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL,  INTENT(IN)    :: PCONC3D    ! Cloud droplet number concentration
REAL, DIMENSION(D%NIJT,D%NKT),     OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(D%NIJT,D%NKT),     OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(D%NIJT,D%NKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQHT ! Hail electric charge at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQHS ! Hail electric charge source
!
!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
LOGICAL :: LLSEA_AND_TOWN, LLCONC
INTEGER :: JIJ, JK, JRR
INTEGER :: IKTB, IKTE, IKB, IIJB, IIJE, IIJT
!
LOGICAL, DIMENSION(D%NIJT,D%NKT) :: LLMICRO ! mask to limit computation
!Arrays for nucleation call outisde of LLMICRO points
REAL,    DIMENSION(D%NIJT, D%NKT) :: ZT ! Temperature
REAL, DIMENSION(D%NIJT, D%NKT) :: ZZ_RVHENI       ! heterogeneous nucleation
REAL, DIMENSION(D%NIJT, D%NKT) :: ZZ_LVFACT, ZZ_LSFACT
REAL, DIMENSION(D%NIJT, D%NKT) :: ZSIGMA_RC
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLC_LCF
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLC_LRC
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLI_LCF
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLI_LRI
!
REAL :: ZINV_TSTEP ! Inverse ov PTSTEP
!For total tendencies computation
REAL, DIMENSION(D%NIJT,D%NKT,0:7) :: ZWR
REAL, DIMENSION(KRR) :: ZICEDRTMIN
!
REAL :: ZDEVIDE, ZRICE, ZZSUM
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZW3D
REAL, DIMENSION(D%NIJT,D%NKT) :: ZZZZ, ZCONC3D
REAL, DIMENSION(D%NIJT) :: ZCONC_TMP

LOGICAL, DIMENSION(D%NIJT,D%NKT) :: LLW3D
REAL, DIMENSION(KRR) :: ZRSMIN
INTEGER :: ISIZE, IPROMA, IGPBLKS, ISIZE2
!
LOGICAL :: LSAVE_MICRO ! if true, microphysical tendencies are saved for cloud electricity
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC),MERGE(IBUNUM-IBUNUM_EXTRA,0,OELEC)) :: &
           ZMICRO_TEND ! Total mixing ratio change, used for electric charge tendencies
LOGICAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)) :: GMASK_ELEC
INTEGER :: IELEC ! nb of points where microphysical tendencies are not null
INTEGER :: JM    ! loop index
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RAIN_ICE', 0, ZHOOK_HANDLE)
!
!*       1.     GENERALITIES
!               ------------
!
IKTB=D%NKTB
IKTE=D%NKTE
IKB=D%NKB
IIJB=D%NIJB
IIJE=D%NIJE
IIJT=D%NIJT
ZICEDRTMIN(1:KRR)=ICED%XRTMIN(1:KRR)
!-------------------------------------------------------------------------------
!

ZINV_TSTEP=1./PTSTEP
!
LLSEA_AND_TOWN=PRESENT(PSEA).AND.PRESENT(PTOWN)
LLCONC=PRESENT(PCONC3D)
!
! LSFACT and LVFACT without exner, and LLMICRO
! LLMICRO is a mask with a True value on points where microphysics is active
DO JRR=1, KRR
  ZRSMIN(JRR) = ZICEDRTMIN(JRR) * ZINV_TSTEP
END DO
LLMICRO(:,:)=.FALSE.



DO JK = IKTB,IKTE
  DO JIJ = IIJB,IIJE
    !LSFACT and LVFACT
    IF (KRR==7) THEN
      ZRICE=PRIT(JIJ,JK)+PRST(JIJ,JK)+PRGT(JIJ,JK)+PRHT(JIJ,JK)
    ELSE
      ZRICE=PRIT(JIJ,JK)+PRST(JIJ,JK)+PRGT(JIJ,JK)
    ENDIF
    ZDEVIDE = CST%XCPD + CST%XCPV*PRVT(JIJ,JK) + CST%XCL*(PRCT(JIJ,JK)+PRRT(JIJ,JK)) + CST%XCI*ZRICE
    ZT(JIJ,JK) = PTHT(JIJ,JK) * PEXN(JIJ,JK)
    ZZ_LSFACT(JIJ,JK)=(CST%XLSTT+(CST%XCPV-CST%XCI)*(ZT(JIJ,JK)-CST%XTT)) / ZDEVIDE
    ZZ_LVFACT(JIJ,JK)=(CST%XLVTT+(CST%XCPV-CST%XCL)*(ZT(JIJ,JK)-CST%XTT)) / ZDEVIDE

    !LLMICRO
    IF ( PARAMI%LOCND2 ) THEN
      IF (KRR==7) THEN
        LLMICRO(JIJ,JK)=PSSIO(JIJ,JK)>ICEP%XFRMIN(12) .OR. &
                        PRCT(JIJ,JK)>ICEP%XFRMIN(13) .OR. &
                        PRRT(JIJ,JK)>ICEP%XFRMIN(13) .OR. &
                        PRIT(JIJ,JK)>ICEP%XFRMIN(13) .OR. &
                        PRST(JIJ,JK)>ICEP%XFRMIN(13) .OR. &
                        PRGT(JIJ,JK)>ICEP%XFRMIN(13) .OR. &
                        PRHT(JIJ,JK)>ICEP%XFRMIN(13)
      ELSE
        LLMICRO(JIJ,JK)=PSSIO(JIJ,JK)>ICEP%XFRMIN(12) .OR. &
                        PRCT(JIJ,JK)>ICEP%XFRMIN(13) .OR. &
                        PRRT(JIJ,JK)>ICEP%XFRMIN(13) .OR. &
                        PRIT(JIJ,JK)>ICEP%XFRMIN(13) .OR. &
                        PRST(JIJ,JK)>ICEP%XFRMIN(13) .OR. &
                        PRGT(JIJ,JK)>ICEP%XFRMIN(13)
      ENDIF
    ELSE
      IF (KRR==7) THEN
        LLMICRO(JIJ,JK)=PRCT(JIJ,JK)>ICED%XRTMIN(2) .OR. &
                        PRRT(JIJ,JK)>ICED%XRTMIN(3) .OR. &
                        PRIT(JIJ,JK)>ICED%XRTMIN(4) .OR. &
                        PRST(JIJ,JK)>ICED%XRTMIN(5) .OR. &
                        PRGT(JIJ,JK)>ICED%XRTMIN(6) .OR. &
                        PRHT(JIJ,JK)>ICED%XRTMIN(7)
      ELSE
        LLMICRO(JIJ,JK)=PRCT(JIJ,JK)>ICED%XRTMIN(2) .OR. &
                        PRRT(JIJ,JK)>ICED%XRTMIN(3) .OR. &
                        PRIT(JIJ,JK)>ICED%XRTMIN(4) .OR. &
                        PRST(JIJ,JK)>ICED%XRTMIN(5) .OR. &
                        PRGT(JIJ,JK)>ICED%XRTMIN(6)
      ENDIF
    ENDIF
  ENDDO
ENDDO
!        1.1    Compute cloud liquid number concentration
!               -----------------------------------------
! Model level height
DO JIJ=IIJB,IIJE
  ZZSUM = 0.
  DO JK=IKTE,IKTB, -1
    ZZSUM = ZZSUM + PDZZ(JIJ,JK)
    ZZZZ(JIJ,JK) = ZZSUM - PDZZ(JIJ,JK)*0.5
  ENDDO
ENDDO

ZCONC3D(:,:) = ICED%XCONC_LAND
IF( LLCONC .AND. PARAMI%LEXCLDROP) THEN
  DO JK=IKTB,IKTE
    DO JIJ=IIJB,IIJE
      ZCONC3D(JIJ,JK)=PCONC3D(JIJ,JK)
    ENDDO
  ENDDO
ELSE
  IF (PARAMI%LEXCLDROP) THEN
   CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'RAIN_ICE', 'WITH LEXCLDROP=TRUE CLOUD DROPLET FIELDS MUST BE PRESENT IN RAIN_ICE')
  END IF
  IF(LLSEA_AND_TOWN .AND. ICEP%XFRMIN(26)<0.001) THEN
    DO JK=IKTB,IKTE
      DO JIJ=IIJB,IIJE
        ZCONC_TMP(JIJ) = PSEA(JIJ)*ICED%XCONC_SEA+(1.-PSEA(JIJ))*ICED%XCONC_LAND
        ZCONC3D(JIJ,JK)= (1.-PTOWN(JIJ))*ZCONC_TMP(JIJ) +PTOWN(JIJ)*ICED%XCONC_URBAN
      ENDDO
    ENDDO
  ENDIF
  IF(ICEP%XFRMIN(26)>0.001)THEN
    DO JK=IKTB,IKTE
      DO JIJ=IIJB,IIJE
        ZCONC3D(JIJ,JK) = ICEP%XFRMIN(26)*PPABST(JIJ,JK)/CST%XP00
        IF(ICEP%XFRMIN(22)>0.001 .AND. LLSEA_AND_TOWN)THEN
          ZCONC_TMP(JIJ)  = MAX(0., MIN(1.,(ZZZZ(JIJ,JK)- ZZZZ(JIJ,IKTE))/&
              &(MAX(0.01,ICEP%XFRMIN(29)-ZZZZ(JIJ,IKTE)))))
          ZCONC3D(JIJ,JK) =  ZCONC3D(JIJ,JK)*(ICEP%XFRMIN(22)*(PSEA(JIJ)*ICEP%XFRMIN(30)+&
              &(1.-PSEA(JIJ)))*(1.-ZCONC_TMP(JIJ)) + ZCONC_TMP(JIJ))
        ENDIF
      ENDDO
    ENDDO
  ENDIF
ENDIF
!
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!               -------------------------------------
!
IF(.NOT. PARAMI%LSEDIM_AFTER) THEN
  CALL ICE4_SEDIMENTATION(D, CST, ICEP, ICED, PARAMI, ELECP, ELECD, BUCONF, &
                         &OELEC, OSEDIM_BEARD, PTSTEP, KRR, PDZZ, PTHVREFZIKB, &
                         &ZZ_LVFACT, ZZ_LSFACT, PRHODREF, PPABST, PTHT, ZT, PRHODJ, &
                         &PTHS, PRVS, PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                         &ZCONC3D, &
                         &PINPRC, PINPRR, PINPRS, PINPRG, &
                         &PQCT, PQRT, PQIT, PQST, PQGT, PQCS, PQRS, PQIS, PQSS, PQGS, PEFIELDW, &
                         &TBUDGETS, KBUDGETS, &
                         &PSEA=PSEA, PTOWN=PTOWN, &
                         &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR, &
                         &PQHT=PQHT, PQHS=PQHS)
ENDIF
!
!
!-------------------------------------------------------------------------------
!
!*       3.     INITIAL VALUES SAVING
!               ---------------------
!

DO JK = IKTB,IKTE
DO JIJ=1, D%NIJT
    !Copy of T variables to keep untouched the prognostic variables
    ZWR(JIJ, JK, ITH)=PTHT(JIJ, JK)
    ZWR(JIJ, JK, IRV)=PRVT(JIJ, JK)
    ZWR(JIJ, JK, IRC)=PRCT(JIJ, JK)
    ZWR(JIJ, JK, IRR)=PRRT(JIJ, JK)
    ZWR(JIJ, JK, IRI)=PRIT(JIJ, JK)
    ZWR(JIJ, JK, IRS)=PRST(JIJ, JK)
    ZWR(JIJ, JK, IRG)=PRGT(JIJ, JK)
    IF (KRR==7) THEN
      ZWR(JIJ, JK, IRH)=PRHT(JIJ, JK)
    ELSE
      ZWR(JIJ, JK, IRH)=0.
    ENDIF

    !Preset for output 3D variables
    IF(PARAMI%LWARM) THEN
      PEVAP3D(JIJ, JK)=0.
    ENDIF
    PRAINFR(JIJ, JK)=0.
END DO
ENDDO

!
!
!*       4.1    COMPUTES THE SLOW COLD PROCESS SOURCES OUTSIDE OF LLMICRO POINTS
!               -----------------------------------------------------------------
!
!The nucleation must be called everywhere
!This call is for points outside of the LLMICR mask, another call is coded in ice4_tendencies

LLW3D(:,:)=.FALSE.

DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE
    IF (.NOT. LLMICRO(JIJ, JK)) THEN
      LLW3D(JIJ, JK)=.TRUE.
      ZW3D(JIJ, JK)=ZZ_LSFACT(JIJ, JK)/PEXN(JIJ, JK)
      PCIT(JIJ,JK)=0. !ri=0 because where are in the not odmicro case
    ELSE
      LLW3D(JIJ, JK)=.FALSE.
    ENDIF
  ENDDO
ENDDO




DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE
    CALL ICE4_NUCLEATION(CST, PARAMI, ICEP, ICED, LLW3D(JIJ, JK), &
                         PTHT(JIJ, JK), PPABST(JIJ, JK), PRHODREF(JIJ, JK), &
                         PEXN(JIJ, JK), ZW3D(JIJ, JK), ZT(JIJ, JK), &
                         PRVT(JIJ, JK), PICLDFR(JIJ, JK), ZZZZ(JIJ, JK), &
                         PCIT(JIJ, JK), ZZ_RVHENI(JIJ, JK))
  ENDDO
ENDDO




DO JK = IKTB, IKTE
  DO JIJ=IIJB, IIJE
    ZZ_RVHENI(JIJ,JK) = MIN(PRVS(JIJ,JK), ZZ_RVHENI(JIJ,JK)/PTSTEP)
  ENDDO
ENDDO

!
!
!*       4.2    COMPUTES PRECIPITATION FRACTION
!               -------------------------------
!
!The ICE4_RAINFR_VERT call was previously in ice4_tendencies to be computed again at each iteration.
!The computation has been moved here to separate (for GPUs) the part of the code
!where column computation can occur (here, alongside with the sedimentation) and
!other routines where computation are only 0D (point by point).
!This is not completly exact but we can think that the precipitation fraction
!diagnostic does not evolve too much during a time-step.
!ICE4_RAINFR_VERT needs the output of ICE4_COMPUTE_PDF; thus this routine
!is called here but it's still called from within ice4_tendencies.
IF (PARAMI%CSUBG_RC_RR_ACCR=='PRFR' .OR. PARAMI%CSUBG_RR_EVAP=='PRFR') THEN
  IF (PARAMI%CSUBG_AUCV_RC=='PDF ' .AND. PARAMI%CSUBG_PR_PDF=='SIGM') THEN
    DO JK = IKTB, IKTE                                                                                                                  
      DO JIJ=IIJB, IIJE
        ZSIGMA_RC(JIJ, JK)=PSIGS(JIJ, JK)**2
      ENDDO
    ENDDO
  ENDIF
  IF (PARAMI%CSUBG_AUCV_RC=='ADJU' .OR. PARAMI%CSUBG_AUCV_RI=='ADJU') THEN
    DO JK = IKTB, IKTE                                                                                                                
      DO JIJ=IIJB, IIJE
        ZHLC_LRC(JIJ, JK) = ZWR(JIJ, JK, IRC) - PHLC_HRC(JIJ, JK)
        ZHLI_LRI(JIJ, JK) = ZWR(JIJ, JK, IRI) - PHLI_HRI(JIJ, JK)
        IF(ZWR(JIJ, JK, IRC)>0.) THEN
          ZHLC_LCF(JIJ, JK) = PCLDFR(JIJ, JK)- PHLC_HCF(JIJ, JK)
        ELSE
          ZHLC_LCF(JIJ, JK)=0.
        ENDIF
        IF(ZWR(JIJ, JK, IRI)>0.) THEN
          ZHLI_LCF(JIJ, JK) = PCLDFR(JIJ, JK)- PHLI_HCF(JIJ, JK)
        ELSE
          ZHLI_LCF(JIJ, JK)=0.
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  !We cannot use ZWR(:,IKTB:IKTE,IRC) which is not contiguous
  CALL ICE4_COMPUTE_PDF(CST, ICEP, ICED, IIJT*(IKTE-IKTB+1), PARAMI%CSUBG_AUCV_RC, PARAMI%CSUBG_AUCV_RI, PARAMI%CSUBG_PR_PDF,&
                        LLMICRO(:,IKTB:IKTE), PRHODREF(:,IKTB:IKTE), PRCT(:,IKTB:IKTE), PRIT(:,IKTB:IKTE), &
                        PCLDFR(:,IKTB:IKTE), ZT(:,IKTB:IKTE), ZSIGMA_RC(:,IKTB:IKTE), &
                        PHLC_HCF(:,IKTB:IKTE), ZHLC_LCF(:,IKTB:IKTE), PHLC_HRC(:,IKTB:IKTE), ZHLC_LRC(:,IKTB:IKTE), &
                        PHLI_HCF(:,IKTB:IKTE), ZHLI_LCF(:,IKTB:IKTE), PHLI_HRI(:,IKTB:IKTE), ZHLI_LRI(:,IKTB:IKTE), &
                        PRAINFR(:,IKTB:IKTE))
!CALL ICE4_COMPUTE_PDF2D(D, CST, ICEP, ICED, PARAMI%CSUBG_AUCV_RC, PARAMI%CSUBG_AUCV_RI, PARAMI%CSUBG_PR_PDF, &
!                        LLMICRO, PRHODREF, ZWR(:,:,IRC), ZWR(:,:,IRI), PCLDFR, ZT, ZSIGMA_RC,&
!                            PHLC_HCF, ZHLC_LCF, PHLC_HRC, ZHLC_LRC, &
!                            PHLI_HCF, ZHLI_LCF, PHLI_HRI, ZHLI_LRI, PRAINFR)
  IF (PRESENT(PRHS)) THEN
    CALL ICE4_RAINFR_VERT(D, ICED, PRAINFR, ZWR(:,:,IRR), &
                         &ZWR(:,:,IRS), ZWR(:,:,IRG), ZWR(:,:,IRH))
  ELSE
    CALL ICE4_RAINFR_VERT(D, ICED, PRAINFR, ZWR(:,:,IRR), &
                         &ZWR(:,:,IRS), ZWR(:,:,IRG)) 
  ENDIF
ELSE
  PRAINFR(:,:)=1.
ENDIF
!
!
!*       5.     TENDENCIES COMPUTATION
!               ----------------------
!
IF(PARAMI%LPACK_MICRO) THEN
  ISIZE=0
  DO JK=1,D%NKT
    DO JIJ=1,D%NIJT
      IF(LLMICRO(JIJ,JK)) ISIZE=ISIZE+1 ! Number of points with active microphysics
    END DO
  END DO
  !PARAMI%NPROMICRO is the requested size for cache_blocking loop
  !IPROMA is the effective size
  !This parameter must be computed here because it is used for array dimensioning in ice4_pack
  IF (PARAMI%NPROMICRO > 0 .AND. ISIZE > 0) THEN
    ! Cache-blocking is active
    ! number of chunks :
    IGPBLKS = (ISIZE-1)/MIN(PARAMI%NPROMICRO,ISIZE)+1
    ! Adjust IPROMA to limit the number of small chunks
    IPROMA=(ISIZE-1)/IGPBLKS+1
  ELSE
    IPROMA=ISIZE ! no cache-blocking
  ENDIF
  ISIZE2=IPROMA
ELSE
  ISIZE=D%NIJT*D%NKT
  IPROMA=0
  ISIZE2=ISIZE
ENDIF
!
!Microphysical tendencies must be saved for some physical parameterizations
IF (OELEC) THEN
  LSAVE_MICRO = .TRUE.
  ZMICRO_TEND(:,:,:) = 0.
ELSE
  LSAVE_MICRO = .FALSE.
END IF
!
!This part is put in another routine to separate pack/unpack operations from computations
CALL ICE4_PACK(D, CST, PARAMI, ICEP, ICED, BUCONF,                   &
               IPROMA, ISIZE, ISIZE2,                                &
               PTSTEP, KRR, LSAVE_MICRO, LLMICRO, OELEC, PEXN,       &
               PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,      &
               PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,               &
               PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,             &
               PEVAP3D,                                              &
               PRAINFR, PSIGS,                                       &
               ZZ_RVHENI, ZZ_LVFACT, ZZ_LSFACT,                      &
               PICLDFR, ZZZZ, ZCONC3D,                               &
               PSSIO, PSSIU, PIFR,                                   &
               ZWR,                                                  &
               TBUDGETS, KBUDGETS,                                   &
               ZMICRO_TEND, PLATHAM_IAGGS, PRHS                      )
!
!
!-------------------------------------------------------------------------------
!
!*       7.     CALL TO PHYSICAL PARAMETERIZATIONS CLOSELY LINKED TO MICROPHYSICS
!               -----------------------------------------------------------------
!
! Cloud electrification, water isotopes and aqueous chemistry need the mixing ratio tendencies
! to compute the evolution of electric charges, water isotopes and ...
!
!*       7.1    Cloud electrification
!
IF (OELEC) THEN
  DO JK = IKTB, IKTE
    DO JIJ = IIJB, IIJE
      DO JM = 1, IBUNUM-IBUNUM_EXTRA
        ZMICRO_TEND(JIJ,JK,JM) = ZMICRO_TEND(JIJ,JK,JM) * ZINV_TSTEP
        !
        ! transfer of electric charges occurs only where transfer of mass is non null
        GMASK_ELEC(JIJ,JK) = GMASK_ELEC(JIJ,JK) .OR. (ZMICRO_TEND(JIJ,JK,JM) .NE. 0.)
      END DO
    END DO
  END DO
  !
  IELEC = COUNT(GMASK_ELEC)
  ! 
  ! RVHENI : ajout de prvheni ?
  ! traitement des deux termes extra ? irwetgh_mr et irsrimcg_mr ?
  IF (KRR == 7) THEN
    CALL ELEC_TENDENCIES(D, CST, ICED, ICEP, ELECD, ELECP,                                    &
                         KRR, IELEC, PTSTEP, GMASK_ELEC,                                      &
                         BUCONF, TBUDGETS, KBUDGETS,                                          &
                         'ICE4', PTHVREFZIKB, PRHODREF, PRHODJ, ZT, PCIT,                     &
                         PRVT, PRCT, PRRT, PRIT, PRST, PRGT,                                  &
                         PQPIT, PQCT, PQRT, PQIT, PQST, PQGT, PQNIT,                          &
                         PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS,                          &
                         ZMICRO_TEND(:,:,IRVHENI_MR), ZMICRO_TEND(:,:,IRRHONG_MR),            &
                         ZMICRO_TEND(:,:,IRIMLTC_MR), ZMICRO_TEND(:,:,IRCHONI),               &
                         ZMICRO_TEND(:,:,IRVDEPS), ZMICRO_TEND(:,:,IRIAGGS),                  &
                         ZMICRO_TEND(:,:,IRIAUTS), ZMICRO_TEND(:,:,IRVDEPG),                  &
                         ZMICRO_TEND(:,:,IRCAUTR), ZMICRO_TEND(:,:,IRCACCR),                  &
                         ZMICRO_TEND(:,:,IRREVAV), ZMICRO_TEND(:,:,IRCRIMSS),                 &
                         ZMICRO_TEND(:,:,IRCRIMSG), ZMICRO_TEND(:,:,IRSRIMCG),                &
                         ZMICRO_TEND(:,:,IRRACCSS), ZMICRO_TEND(:,:,IRRACCSG),                &
                         ZMICRO_TEND(:,:,IRSACCRG), ZMICRO_TEND(:,:,IRSMLTG),                 &
                         ZMICRO_TEND(:,:,IRICFRRG), ZMICRO_TEND(:,:,IRRCFRIG),                &
                         ZMICRO_TEND(:,:,IRCWETG), ZMICRO_TEND(:,:,IRIWETG),                  &
                         ZMICRO_TEND(:,:,IRRWETG), ZMICRO_TEND(:,:,IRSWETG),                  &
                         ZMICRO_TEND(:,:,IRCDRYG), ZMICRO_TEND(:,:,IRIDRYG),                  &
                         ZMICRO_TEND(:,:,IRRDRYG), ZMICRO_TEND(:,:,IRSDRYG),                  &
                         ZMICRO_TEND(:,:,IRGMLTR), ZMICRO_TEND(:,:,IRCBERI),                  &
                         PRCMLTSR=ZMICRO_TEND(:,:,IRCMLTSR), PRICFRR=ZMICRO_TEND(:,:,IRICFRR),&
                         PRWETGH=ZMICRO_TEND(:,:,IRWETGH),                                    &
                         PRCWETH=ZMICRO_TEND(:,:,IRCWETH), PRIWETH=ZMICRO_TEND(:,:,IRIWETH),  &
                         PRSWETH=ZMICRO_TEND(:,:,IRSWETH),                                    &
                         PRGWETH=ZMICRO_TEND(:,:,IRGWETH), PRRWETH=ZMICRO_TEND(:,:,IRRWETH),  &
                         PRCDRYH=ZMICRO_TEND(:,:,IRCDRYH), PRIDRYH=ZMICRO_TEND(:,:,IRIDRYH),  &
                         PRSDRYH=ZMICRO_TEND(:,:,IRSDRYH),                                    &
                         PRRDRYH=ZMICRO_TEND(:,:,IRRDRYH), PRGDRYH=ZMICRO_TEND(:,:,IRGDRYH),  &
                         PRDRYHG=ZMICRO_TEND(:,:,IRDRYHG), PRHMLTR=ZMICRO_TEND(:,:,IRHMLTR),  &
                         PRHT=PRHT, PRHS=PRHS, PQHT=PQHT, PQHS=PQHS                       )
  ELSE
    CALL ELEC_TENDENCIES(D, CST, ICED, ICEP, ELECD, ELECP,                                   &
                         KRR, ISIZE, PTSTEP, LLMICRO,                                        &
                         BUCONF, TBUDGETS, KBUDGETS,                                         &
                         'ICE3', PTHVREFZIKB, PRHODREF, PRHODJ, ZT, PCIT,                    &
                         PRVT, PRCT, PRRT, PRIT, PRST, PRGT,                                 &
                         PQPIT, PQCT, PQRT, PQIT, PQST, PQGT, PQNIT,                         &
                         PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS,                         &
                         ZMICRO_TEND(:,:,IRVHENI_MR), ZMICRO_TEND(:,:,IRRHONG_MR),           &
                         ZMICRO_TEND(:,:,IRIMLTC_MR), ZMICRO_TEND(:,:,IRCHONI),              &
                         ZMICRO_TEND(:,:,IRVDEPS), ZMICRO_TEND(:,:,IRIAGGS),                 &
                         ZMICRO_TEND(:,:,IRIAUTS), ZMICRO_TEND(:,:,IRVDEPG),                 &
                         ZMICRO_TEND(:,:,IRCAUTR), ZMICRO_TEND(:,:,IRCACCR),                 &
                         ZMICRO_TEND(:,:,IRREVAV), ZMICRO_TEND(:,:,IRCRIMSS),                &
                         ZMICRO_TEND(:,:,IRCRIMSG), ZMICRO_TEND(:,:,IRSRIMCG),               &
                         ZMICRO_TEND(:,:,IRRACCSS), ZMICRO_TEND(:,:,IRRACCSG),               &
                         ZMICRO_TEND(:,:,IRSACCRG), ZMICRO_TEND(:,:,IRSMLTG),                &
                         ZMICRO_TEND(:,:,IRICFRRG), ZMICRO_TEND(:,:,IRRCFRIG),               &
                         ZMICRO_TEND(:,:,IRCWETG), ZMICRO_TEND(:,:,IRIWETG),                 &
                         ZMICRO_TEND(:,:,IRRWETG), ZMICRO_TEND(:,:,IRSWETG),                 &
                         ZMICRO_TEND(:,:,IRCDRYG), ZMICRO_TEND(:,:,IRIDRYG),                 &
                         ZMICRO_TEND(:,:,IRRDRYG), ZMICRO_TEND(:,:,IRSDRYG),                 &
                         ZMICRO_TEND(:,:,IRGMLTR), ZMICRO_TEND(:,:,IRCBERI),                 &
                         PRCMLTSR=ZMICRO_TEND(:,:,IRCMLTSR), PRICFRR=ZMICRO_TEND(:,:,IRICFRR))
  END IF
END IF
!
!
!*       7.2    Water isotopologues
!
!
!*       7.3    Aqueous chemistry
!
!
!-------------------------------------------------------------------------------
!
!*       8.     TOTAL TENDENCIES
!               ----------------
!
!
!***     8.1    total tendencies limited by available species
!


DO JK = IKTB, IKTE
  DO JIJ=IIJB, IIJE
    !LV/LS
    ZZ_LSFACT(JIJ,JK)=ZZ_LSFACT(JIJ,JK)/PEXNREF(JIJ,JK)
    ZZ_LVFACT(JIJ,JK)=ZZ_LVFACT(JIJ,JK)/PEXNREF(JIJ,JK)

    !Hydrometeor tendencies is the difference between new state and old state (can be negative)
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

!-------------------------------------------------------------------------------
!
!***     8.2    Negative corrections
!
!NOTE:
!  This call cannot be moved before the preeceding budget calls because,
!  with AROME, the BUDGET_STORE_INIT does nothing. The equivalent is done only
!  once before the physics call and copies of the S variables evolve automatically
!  internally to the budget (DDH) machinery at each BUDGET_STORE_ADD and
!  BUDGET_STORE_END calls. Thus, the difference between the DDH internal version
!  of the S variables and the S variables used in the folowing BUDGET_STORE_END
!  call must only be due to the correction of negativities.
!
IF(BUCONF%LBU_ENABLE) THEN
  IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'CORR', PTHS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RV) CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'CORR', PRVS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'CORR', PRCS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'CORR', PRRS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'CORR', PRIS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'CORR', PRSS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'CORR', PRGS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RH .AND. KRR==7) CALL TBUDGETS(NBUDGET_RH)%PTR%INIT_PHY(D, 'CORR', PRHS(:, :)*PRHODJ(:, :))
END IF
!++cb-- ajouter les bilans pour l'elec !!!

!We correct negativities with conservation
CALL ICE4_CORRECT_NEGATIVITIES(D, ICED, KRR, PRVS, PRCS, PRRS, &
                              &PRIS, PRSS, PRGS, &
                              &PTHS, ZZ_LVFACT, ZZ_LSFACT, PRHS)
!CALL ICE4_CORRECT_NEGATIVITIES(D, ICED, KRR, OELEC, PRVS, PRCS, PRRS, &
!                              &PRIS, PRSS, PRGS, &
!                              &PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS, &                              
!                              &PTHS, ZZ_LVFACT, ZZ_LSFACT, PRHS, PQHS)

IF (BUCONF%LBU_ENABLE) THEN
  IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'CORR', PTHS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RV) CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'CORR', PRVS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'CORR', PRCS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'CORR', PRRS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'CORR', PRIS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'CORR', PRSS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'CORR', PRGS(:, :)*PRHODJ(:, :))
  IF (BUCONF%LBUDGET_RH .AND. KRR==7) CALL TBUDGETS(NBUDGET_RH)%PTR%END_PHY(D, 'CORR', PRHS(:, :)*PRHODJ(:, :))
END IF
!
!-------------------------------------------------------------------------------
!
!*       9.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!               -------------------------------------
!
IF(PARAMI%LSEDIM_AFTER) THEN
  CALL ICE4_SEDIMENTATION(D, CST, ICEP, ICED, PARAMI, ELECP, ELECD, BUCONF, &
                         &OELEC, OSEDIM_BEARD, PTSTEP, KRR, PDZZ, PTHVREFZIKB, &
                         &ZZ_LVFACT, ZZ_LSFACT, PRHODREF, PPABST, PTHT, ZT, PRHODJ, &
                         &PTHS, PRVS, PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                         &ZCONC3D, &
                         &PINPRC, PINPRR, PINPRS, PINPRG, &
                         &PQCT, PQRT, PQIT, PQST, PQGT, PQCS, PQRS, PQIS, PQSS, PQGS, PEFIELDW, &
                         &TBUDGETS, KBUDGETS, &
                         &PSEA=PSEA, PTOWN=PTOWN, &
                         &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR, &
                         &PQHT=PQHT, PQHS=PQHS)
        
  !"sedimentation" of rain fraction


  DO JK = IKTB, IKTE
    DO JIJ=IIJB,IIJE
      ZWR(JIJ,JK,IRR)=PRRS(JIJ,JK)*PTSTEP
      ZWR(JIJ,JK,IRS)=PRSS(JIJ,JK)*PTSTEP
      ZWR(JIJ,JK,IRG)=PRGS(JIJ,JK)*PTSTEP
      IF(KRR==7) THEN
        ZWR(JIJ,JK,IRH)=PRHS(JIJ,JK)*PTSTEP
      ENDIF
    ENDDO
  ENDDO

  IF (PRESENT(PRHS)) THEN
    CALL ICE4_RAINFR_VERT(D, ICED, PRAINFR, ZWR(:,:,IRR), &
                         &ZWR(:,:,IRS), ZWR(:,:,IRG), ZWR(:,:,IRH))
  ELSE
    CALL ICE4_RAINFR_VERT(D, ICED, PRAINFR, ZWR(:,:,IRR), &
                         &ZWR(:,:,IRS), ZWR(:,:,IRG)) 
  ENDIF
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       10.    COMPUTE THE FOG DEPOSITION TERM 
!               -------------------------------------
!
IF (PARAMI%LDEPOSC) THEN !cloud water deposition on vegetation
  IF (BUCONF%LBU_ENABLE .AND. BUCONF%LBUDGET_RC) &
     & CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'DEPO', PRCS(:, :)*PRHODJ(:, :))

  PINDEP(:)=0.
!DEC$ IVDEP

  DO JIJ = IIJB, IIJE
    PINDEP(JIJ) = PARAMI%XVDEPOSC * PRCT(JIJ, IKB) * PRHODREF(JIJ, IKB) / CST%XRHOLW
    PRCS(JIJ, IKB) = PRCS(JIJ, IKB) - PARAMI%XVDEPOSC * PRCT(JIJ, IKB) / PDZZ(JIJ, IKB)
    PINPRC(JIJ) = PINPRC(JIJ) + PINDEP(JIJ)
  ENDDO


  IF (BUCONF%LBU_ENABLE .AND. BUCONF%LBUDGET_RC) &
     & CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'DEPO', PRCS(:, :)*PRHODJ(:, :))
ENDIF

IF (LHOOK) CALL DR_HOOK('RAIN_ICE', 1, ZHOOK_HANDLE)
!
CONTAINS
INCLUDE "ice4_nucleation.func.h"
END SUBROUTINE RAIN_ICE
