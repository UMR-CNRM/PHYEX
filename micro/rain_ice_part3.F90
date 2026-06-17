!MNH_LIC Copyright 1995-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
       SUBROUTINE RAIN_ICE_PART3 ( D, CST, PARAMI, ICEP, ICED, ELECP, ELECD, BUCONF,     &
                                   OELEC, OSEDIM_BEARD, PTHVREFZIKB,                     &
                                   PTSTEP, KRR,                                    &
                                   PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, &
                                   PTHT, PRT, PTHS, PRS,                                 &
                                   PINPRC, PINPRR,                              &
                                   PINPRS, PINPRG, PINDEP, PRAINFR,               &
                                   TBUDGETS, KBUDGETS,                                   &
                                   PQPIT, PQCT, PQRT, PQIT, PQST, PQGT, PQNIT,           &
                                   PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS,           &
                                   PEFIELDW,                              &
                                   PSEA, PTOWN,                                 &
                                   PINPRH, PFPR, PQHT, PQHS,                             &
                                   PT, PLVFACT, PLSFACT, PWR, PWTH, PCONC3D, PRVHENI,   &
                                   ZBUDGETS )
!     #############################################################################
!!!****  * -  compute the explicit microphysical sources
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
!!                              version threshold (Chaboureau and Pinty GRL 2006)
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
USE MODD_FIELDS_ADDRESS, ONLY: IBUNUM, IBUNUM_EXTRA, IRC, IRCACCR, IRCAUTR, IRCBERI, IRCDRYG, &
                               IRCDRYH, IRCHONI, IRCMLTSR, IRCRIMSG, IRCRIMSS, IRCWETG, IRCWETH, &
                               IRDRYHG, IRG, IRGDRYH, IRGMLTR, IRGWETH, IRH, IRHMLTR, IRI, IRIAGGS, &
                               IRIAUTS, IRICFRR, IRICFRRG, IRIDRYG, IRIDRYH, IRIMLTC_MR, IRIWETG, &
                               IRIWETH, IRR, IRRACCSG, IRRACCSS, IRRCFRIG, IRRDRYG, IRRDRYH, IRREVAV, &
                               IRRHONG_MR, IRRWETG, IRRWETH, IRS, IRSACCRG, IRSDRYG, IRSDRYH, IRSMLTG, &
                               IRSRIMCG, IRSWETG, IRSWETH, IRV, IRVDEPG, IRVDEPS, IRVHENI_MR, IRWETGH

USE MODE_ICE4_RAINFR_VERT,   ONLY: ICE4_RAINFR_VERT
USE MODE_ICE4_SEDIMENTATION, ONLY: ICE4_SEDIMENTATION
USE MODE_ICE4_BUDGETS,       ONLY: ICE4_BUDGETS
USE MODE_ICE4_CORRECT_NEGATIVITIES, ONLY: ICE4_CORRECT_NEGATIVITIES

USE MODE_ELEC_TENDENCIES, ONLY : ELEC_TENDENCIES

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
REAL,                     INTENT(IN)    :: PTHVREFZIKB ! Reference thv at IKB for electricity
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIJT,D%NKT,KRR),   INTENT(IN)    :: PRT    ! m.r. at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(D%NIJT,D%NKT,KRR),   INTENT(INOUT) :: PRS    ! m.r. source
!
REAL, DIMENSION(D%NIJT), INTENT(INOUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(D%NIJT), INTENT(INOUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(D%NIJT), INTENT(INOUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(D%NIJT), INTENT(INOUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(D%NIJT), INTENT(INOUT) :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRAINFR !Precipitation fraction
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
!
! optional variables
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(INOUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(D%NIJT,D%NKT,KRR), OPTIONAL, INTENT(INOUT)   :: PFPR    ! upper-air precipitation fluxes
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQHT ! Hail electric charge at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQHS ! Hail electric charge source
!
! Bridge variables from part1/part2
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PT
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PLVFACT
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PLSFACT
REAL, DIMENSION(D%NIJT,D%NKT,7), INTENT(INOUT) :: PWR
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PWTH
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PCONC3D
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRVHENI
!
! Budget array
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC .OR. BUCONF%LBU_ENABLE), &
                MERGE(D%NKT,0,OELEC .OR. BUCONF%LBU_ENABLE), &
                MERGE(IBUNUM-IBUNUM_EXTRA,0,OELEC .OR. BUCONF%LBU_ENABLE)), &
           INTENT(INOUT) :: ZBUDGETS
!
!*       0.2   Local variables
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JIJ, JK, JM
INTEGER :: IKTB, IKTE, IKB, IIJB, IIJE, IIJT
REAL :: ZINV_TSTEP
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWKBUD
LOGICAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)) :: GMASK_ELEC
INTEGER :: IELEC
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RAIN_ICE_PART3', 0, ZHOOK_HANDLE)
!
IKTB=D%NKTB
IKTE=D%NKTE
IKB=D%NKB
IIJB=D%NIJB
IIJE=D%NIJE
IIJT=D%NIJT
ZINV_TSTEP=1./PTSTEP
!
!
!*       5.     TENDENCIES COMPUTATION (budget part)
!               ----------------------
!
IF(BUCONF%LBU_ENABLE .OR. OELEC) THEN
  DO JK = IKTB, IKTE
    DO JIJ = IIJB, IIJE
      DO JM = 1, IBUNUM-IBUNUM_EXTRA
        ZBUDGETS(JIJ,JK,JM) = ZBUDGETS(JIJ,JK,JM) * ZINV_TSTEP
      ENDDO
    ENDDO
  ENDDO
ENDIF
!
IF(BUCONF%LBU_ENABLE) THEN
  CALL ICE4_BUDGETS(D, PARAMI, BUCONF, KRR, &
                    PLVFACT, PLSFACT, PRHODJ, PEXNREF, &
                    PRVHENI, ZBUDGETS, &
                    TBUDGETS, KBUDGETS)
ENDIF
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
  GMASK_ELEC(:,:)=.FALSE.
  DO JK = IKTB, IKTE
    DO JIJ = IIJB, IIJE
      DO JM = 1, IBUNUM-IBUNUM_EXTRA
        ! transfer of electric charges occurs only where transfer of mass is non null
        GMASK_ELEC(JIJ,JK) = GMASK_ELEC(JIJ,JK) .OR. (ZBUDGETS(JIJ,JK,JM) .NE. 0.)
      END DO
    END DO
  END DO
  !
  IELEC = COUNT(GMASK_ELEC)
  !
  IF (KRR == 7) THEN
    CALL ELEC_TENDENCIES(D, CST, ICED, ICEP, ELECD, ELECP,                                    &
                         KRR, IELEC, PTSTEP, GMASK_ELEC,                                      &
                         BUCONF, TBUDGETS, KBUDGETS,                                          &
                         'ICE4', PTHVREFZIKB, PRHODREF, PRHODJ, PT, PCIT,                     &
                         PRT(:,:,IRV), PRT(:,:,IRC), PRT(:,:,IRR), PRT(:,:,IRI), PRT(:,:,IRS), PRT(:,:,IRG), &
                         PQPIT, PQCT, PQRT, PQIT, PQST, PQGT, PQNIT,                          &
                         PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS,                          &
                         ZBUDGETS(:,:,IRVHENI_MR), ZBUDGETS(:,:,IRRHONG_MR),            &
                         ZBUDGETS(:,:,IRIMLTC_MR), ZBUDGETS(:,:,IRCHONI),               &
                         ZBUDGETS(:,:,IRVDEPS), ZBUDGETS(:,:,IRIAGGS),                  &
                         ZBUDGETS(:,:,IRIAUTS), ZBUDGETS(:,:,IRVDEPG),                  &
                         ZBUDGETS(:,:,IRCAUTR), ZBUDGETS(:,:,IRCACCR),                  &
                         ZBUDGETS(:,:,IRREVAV), ZBUDGETS(:,:,IRCRIMSS),                 &
                         ZBUDGETS(:,:,IRCRIMSG), ZBUDGETS(:,:,IRSRIMCG),                &
                         ZBUDGETS(:,:,IRRACCSS), ZBUDGETS(:,:,IRRACCSG),                &
                         ZBUDGETS(:,:,IRSACCRG), ZBUDGETS(:,:,IRSMLTG),                 &
                         ZBUDGETS(:,:,IRICFRRG), ZBUDGETS(:,:,IRRCFRIG),                &
                         ZBUDGETS(:,:,IRCWETG), ZBUDGETS(:,:,IRIWETG),                  &
                         ZBUDGETS(:,:,IRRWETG), ZBUDGETS(:,:,IRSWETG),                  &
                         ZBUDGETS(:,:,IRCDRYG), ZBUDGETS(:,:,IRIDRYG),                  &
                         ZBUDGETS(:,:,IRRDRYG), ZBUDGETS(:,:,IRSDRYG),                  &
                         ZBUDGETS(:,:,IRGMLTR), ZBUDGETS(:,:,IRCBERI),                  &
                         PRCMLTSR=ZBUDGETS(:,:,IRCMLTSR), PRICFRR=ZBUDGETS(:,:,IRICFRR),&
                         PRWETGH=ZBUDGETS(:,:,IRWETGH),                                    &
                         PRCWETH=ZBUDGETS(:,:,IRCWETH), PRIWETH=ZBUDGETS(:,:,IRIWETH),  &
                         PRSWETH=ZBUDGETS(:,:,IRSWETH),                                    &
                         PRGWETH=ZBUDGETS(:,:,IRGWETH), PRRWETH=ZBUDGETS(:,:,IRRWETH),  &
                         PRCDRYH=ZBUDGETS(:,:,IRCDRYH), PRIDRYH=ZBUDGETS(:,:,IRIDRYH),  &
                         PRSDRYH=ZBUDGETS(:,:,IRSDRYH),                                    &
                         PRRDRYH=ZBUDGETS(:,:,IRRDRYH), PRGDRYH=ZBUDGETS(:,:,IRGDRYH),  &
                         PRDRYHG=ZBUDGETS(:,:,IRDRYHG), PRHMLTR=ZBUDGETS(:,:,IRHMLTR),  &
                         PRHT=PRT(:,:,IRH), PQHT=PQHT, PQHS=PQHS                       )
  ELSE
    CALL ELEC_TENDENCIES(D, CST, ICED, ICEP, ELECD, ELECP,                                   &
                         KRR, IELEC, PTSTEP, GMASK_ELEC,                                     &
                         BUCONF, TBUDGETS, KBUDGETS,                                         &
                         'ICE3', PTHVREFZIKB, PRHODREF, PRHODJ, PT, PCIT,                    &
                         PRT(:,:,IRV), PRT(:,:,IRC), PRT(:,:,IRR), PRT(:,:,IRI), PRT(:,:,IRS), PRT(:,:,IRG), &
                         PQPIT, PQCT, PQRT, PQIT, PQST, PQGT, PQNIT,                         &
                         PQPIS, PQCS, PQRS, PQIS, PQSS, PQGS, PQNIS,                         &
                         ZBUDGETS(:,:,IRVHENI_MR), ZBUDGETS(:,:,IRRHONG_MR),           &
                         ZBUDGETS(:,:,IRIMLTC_MR), ZBUDGETS(:,:,IRCHONI),              &
                         ZBUDGETS(:,:,IRVDEPS), ZBUDGETS(:,:,IRIAGGS),                 &
                         ZBUDGETS(:,:,IRIAUTS), ZBUDGETS(:,:,IRVDEPG),                 &
                         ZBUDGETS(:,:,IRCAUTR), ZBUDGETS(:,:,IRCACCR),                 &
                         ZBUDGETS(:,:,IRREVAV), ZBUDGETS(:,:,IRCRIMSS),                &
                         ZBUDGETS(:,:,IRCRIMSG), ZBUDGETS(:,:,IRSRIMCG),               &
                         ZBUDGETS(:,:,IRRACCSS), ZBUDGETS(:,:,IRRACCSG),               &
                         ZBUDGETS(:,:,IRSACCRG), ZBUDGETS(:,:,IRSMLTG),                &
                         ZBUDGETS(:,:,IRICFRRG), ZBUDGETS(:,:,IRRCFRIG),               &
                         ZBUDGETS(:,:,IRCWETG), ZBUDGETS(:,:,IRIWETG),                 &
                         ZBUDGETS(:,:,IRRWETG), ZBUDGETS(:,:,IRSWETG),                 &
                         ZBUDGETS(:,:,IRCDRYG), ZBUDGETS(:,:,IRIDRYG),                 &
                         ZBUDGETS(:,:,IRRDRYG), ZBUDGETS(:,:,IRSDRYG),                 &
                         ZBUDGETS(:,:,IRGMLTR), ZBUDGETS(:,:,IRCBERI),                 &
                         PRCMLTSR=ZBUDGETS(:,:,IRCMLTSR), PRICFRR=ZBUDGETS(:,:,IRICFRR))
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
    PLSFACT(JIJ,JK)=PLSFACT(JIJ,JK)/PEXNREF(JIJ,JK)
    PLVFACT(JIJ,JK)=PLVFACT(JIJ,JK)/PEXNREF(JIJ,JK)

    !Hydrometeor tendencies is the difference between new state and old state (can be negative)
    PWR(JIJ,JK,IRV)=(PWR(JIJ,JK,IRV)-PRT(JIJ,JK,IRV))*ZINV_TSTEP
    PWR(JIJ,JK,IRC)=(PWR(JIJ,JK,IRC)-PRT(JIJ,JK,IRC))*ZINV_TSTEP
    PWR(JIJ,JK,IRR)=(PWR(JIJ,JK,IRR)-PRT(JIJ,JK,IRR))*ZINV_TSTEP
    PWR(JIJ,JK,IRI)=(PWR(JIJ,JK,IRI)-PRT(JIJ,JK,IRI))*ZINV_TSTEP
    PWR(JIJ,JK,IRS)=(PWR(JIJ,JK,IRS)-PRT(JIJ,JK,IRS))*ZINV_TSTEP
    PWR(JIJ,JK,IRG)=(PWR(JIJ,JK,IRG)-PRT(JIJ,JK,IRG))*ZINV_TSTEP
    IF(KRR==7) THEN
      PWR(JIJ,JK,IRH)=(PWR(JIJ,JK,IRH)-PRT(JIJ,JK,IRH))*ZINV_TSTEP
    ENDIF

    !Theta tendency computed from hydrometeors tendencies
    PWTH(JIJ,JK) = (PWR(JIJ,JK,IRC)+PWR(JIJ,JK,IRR))*PLVFACT(JIJ,JK)+ &
                 & (PWR(JIJ,JK,IRI)+PWR(JIJ,JK,IRS)+PWR(JIJ,JK,IRG)+ &
                 &  PWR(JIJ,JK,IRH))*PLSFACT(JIJ,JK)

    !We apply these tendencies to the S variables
    !including the nucleation part
    PTHS(JIJ,JK) = PTHS(JIJ,JK) + PWTH(JIJ,JK)+PRVHENI(JIJ,JK)*PLSFACT(JIJ,JK)
    PRS(JIJ,JK,IRV) = PRS(JIJ,JK,IRV) + PWR(JIJ,JK,IRV)-PRVHENI(JIJ,JK)
    PRS(JIJ,JK,IRC) = PRS(JIJ,JK,IRC) + PWR(JIJ,JK,IRC)
    PRS(JIJ,JK,IRR) = PRS(JIJ,JK,IRR) + PWR(JIJ,JK,IRR)
    PRS(JIJ,JK,IRI) = PRS(JIJ,JK,IRI) + PWR(JIJ,JK,IRI)+PRVHENI(JIJ,JK)
    PRS(JIJ,JK,IRS) = PRS(JIJ,JK,IRS) + PWR(JIJ,JK,IRS)
    PRS(JIJ,JK,IRG) = PRS(JIJ,JK,IRG) + PWR(JIJ,JK,IRG)
    IF (KRR==7) THEN
      PRS(JIJ,JK,IRH) = PRS(JIJ,JK,IRH) + PWR(JIJ,JK,IRH)
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
  IF (BUCONF%LBUDGET_TH) ZWKBUD(:, :) = PTHS(:, :)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RV) ZWKBUD(:, :) = PRS(:, :, IRV)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RV) CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RC) ZWKBUD(:, :) = PRS(:, :, IRC)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RR) ZWKBUD(:, :) = PRS(:, :, IRR)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%INIT_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RI) ZWKBUD(:, :) = PRS(:, :, IRI)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RS) ZWKBUD(:, :) = PRS(:, :, IRS)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%INIT_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RG) ZWKBUD(:, :) = PRS(:, :, IRG)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%INIT_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RH .AND. KRR==7) ZWKBUD(:, :) = PRS(:, :, IRH)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RH .AND. KRR==7) CALL TBUDGETS(NBUDGET_RH)%PTR%INIT_PHY(D, 'CORR', ZWKBUD)
END IF
!++cb-- ajouter les bilans pour l'elec !!!

!We correct negativities with conservation
CALL ICE4_CORRECT_NEGATIVITIES(D, ICED, KRR, PRS, &
                              &PTHS, PLVFACT, PLSFACT)

IF (BUCONF%LBU_ENABLE) THEN
  IF (BUCONF%LBUDGET_TH) ZWKBUD(:, :) = PTHS(:, :)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_TH) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RV) ZWKBUD(:, :) = PRS(:, :, IRV)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RV) CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RC) ZWKBUD(:, :) = PRS(:, :, IRC)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RC) CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RR) ZWKBUD(:, :) = PRS(:, :, IRR)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RR) CALL TBUDGETS(NBUDGET_RR)%PTR%END_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RI) ZWKBUD(:, :) = PRS(:, :, IRI)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RI) CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RS) ZWKBUD(:, :) = PRS(:, :, IRS)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RS) CALL TBUDGETS(NBUDGET_RS)%PTR%END_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RG) ZWKBUD(:, :) = PRS(:, :, IRG)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RG) CALL TBUDGETS(NBUDGET_RG)%PTR%END_PHY(D, 'CORR', ZWKBUD)
  IF (BUCONF%LBUDGET_RH .AND. KRR==7) ZWKBUD(:, :) = PRS(:, :, IRH)*PRHODJ(:, :)
  IF (BUCONF%LBUDGET_RH .AND. KRR==7) CALL TBUDGETS(NBUDGET_RH)%PTR%END_PHY(D, 'CORR', ZWKBUD)
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
                         &PLVFACT, PLSFACT, PRHODREF, PPABST, PTHT, PT, PRHODJ, &
                         &PTHS, PRT, PRS, PCONC3D, &
                         &PINPRC, PINPRR, PINPRS, PINPRG, &
                         &PQCT, PQRT, PQIT, PQST, PQGT, PQCS, PQRS, PQIS, PQSS, PQGS, PEFIELDW, &
                         &TBUDGETS, KBUDGETS, &
                         &PSEA=PSEA, PTOWN=PTOWN, &
                         &PINPRH=PINPRH, PFPR=PFPR, &
                         &PQHT=PQHT, PQHS=PQHS)
        
  !"sedimentation" of rain fraction


  DO JK = IKTB, IKTE
    DO JIJ=IIJB,IIJE
      PWR(JIJ,JK,IRR)=PRS(JIJ,JK,IRR)*PTSTEP
      PWR(JIJ,JK,IRS)=PRS(JIJ,JK,IRS)*PTSTEP
      PWR(JIJ,JK,IRG)=PRS(JIJ,JK,IRG)*PTSTEP
      IF(KRR==7) THEN
        PWR(JIJ,JK,IRH)=PRS(JIJ,JK,IRH)*PTSTEP
      ENDIF
    ENDDO
  ENDDO

  IF (KRR==7) THEN
    CALL ICE4_RAINFR_VERT(D, ICED, PRAINFR, PWR(:,:,IRR), &
                         &PWR(:,:,IRS), PWR(:,:,IRG), PWR(:,:,IRH))
  ELSE
    CALL ICE4_RAINFR_VERT(D, ICED, PRAINFR, PWR(:,:,IRR), &
                         &PWR(:,:,IRS), PWR(:,:,IRG)) 
  ENDIF
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       10.    COMPUTE THE FOG DEPOSITION TERM 
!               -------------------------------------
!
IF (PARAMI%LDEPOSC) THEN !cloud water deposition on vegetation
  IF (BUCONF%LBU_ENABLE .AND. BUCONF%LBUDGET_RC) ZWKBUD(:, :) = PRS(:, :, IRC)*PRHODJ(:, :)
  IF (BUCONF%LBU_ENABLE .AND. BUCONF%LBUDGET_RC) &
     & CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'DEPO', ZWKBUD)

  PINDEP(:)=0.
!DEC$ IVDEP

  DO JIJ = IIJB, IIJE
    PINDEP(JIJ) = PARAMI%XVDEPOSC * PRT(JIJ, IKB, IRC) * PRHODREF(JIJ, IKB) / CST%XRHOLW
    PRS(JIJ, IKB,IRC) = PRS(JIJ, IKB, IRC) - PARAMI%XVDEPOSC * PRT(JIJ, IKB, IRC) / PDZZ(JIJ, IKB)
    PINPRC(JIJ) = PINPRC(JIJ) + PINDEP(JIJ)
  ENDDO


  IF (BUCONF%LBU_ENABLE .AND. BUCONF%LBUDGET_RC) ZWKBUD(:, :) = PRS(:, :, IRC)*PRHODJ(:, :)
  IF (BUCONF%LBU_ENABLE .AND. BUCONF%LBUDGET_RC) &
     & CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'DEPO', ZWKBUD)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('RAIN_ICE_PART3', 1, ZHOOK_HANDLE)
!
END SUBROUTINE RAIN_ICE_PART3
