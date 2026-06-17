!MNH_LIC Copyright 1995-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
       SUBROUTINE RAIN_ICE_PART1 ( D, CST, PARAMI, ICEP, ICED, ELECP, ELECD, BUCONF,     &
                                   OELEC, OSEDIM_BEARD, PTHVREFZIKB,                     &
                                   PTSTEP, KRR, PEXN,                                    &
                                   PDZZ, PRHODJ, PRHODREF, PPABST, PCIT, PCLDFR,         &
                                   PICLDFR, PSSIO,                                        &
                                   PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,               &
                                   PTHT, PRT, PTHS, PRS,                                 &
                                   PINPRC, PINPRR, PEVAP3D,                              &
                                   PINPRS, PINPRG, PRAINFR, PSIGS,               &
                                   TBUDGETS, KBUDGETS,                                   &
                                   PQCT, PQRT, PQIT, PQST, PQGT,           &
                                   PQCS, PQRS, PQIS, PQSS, PQGS,           &
                                   PEFIELDW,                              &
                                   PSEA, PTOWN, PCONC3D,                                 &
                                   PINPRH, PFPR, PQHT, PQHS,                             &
                                   PT, PLVFACT, PLSFACT, PWR, PWTH, PZCONC3D, OMICRO,    &
                                   PRVHENI, PZZZZ )
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
USE MODD_BUDGET,         ONLY: TBUDGETDATA_PTR, TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_ELEC_PARAM,     ONLY: ELEC_PARAM_t
USE MODD_ELEC_DESCR,     ONLY: ELEC_DESCR_t
USE MODD_FIELDS_ADDRESS, ONLY: IRC, IRG, IRH, IRI, IRR, IRS, IRV

USE MODE_MSG,                ONLY: PRINT_MSG, NVERB_FATAL

USE MODE_ICE4_RAINFR_VERT,   ONLY: ICE4_RAINFR_VERT
USE MODE_ICE4_COMPUTE_PDF,   ONLY: ICE4_COMPUTE_PDF
USE MODE_ICE4_SEDIMENTATION, ONLY: ICE4_SEDIMENTATION

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
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PICLDFR ! Ice cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSSIO   ! Super-saturation with respect to ice in the supersaturated fraction
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLC_HRC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLC_HCF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLI_HRI
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLI_HCF
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIJT,D%NKT,KRR),   INTENT(IN)    :: PRT    ! m.r. at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(D%NIJT,D%NKT,KRR),   INTENT(INOUT) :: PRS    ! m.r. source
!
REAL, DIMENSION(D%NIJT), INTENT(OUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(D%NIJT), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(D%NIJT), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(D%NIJT), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT) :: PRAINFR !Precipitation fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
!
TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
! scalar variables for cloud electricity
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQCT   ! Cloud droplet |
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQRT   ! Rain          | electric
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQIT   ! Ice crystals  |  charge 
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQST   ! Snow          |   at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQGT   ! Graupel       |
!
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQCS   ! Cloud droplet |
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQRS   ! Rain          | electric
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQIS   ! Ice crystals  |  charge 
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQSS   ! Snow          |  source
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQGS   ! Graupel       |
!
REAL, DIMENSION(MERGE(D%NIJT,0,OSEDIM_BEARD),MERGE(D%NKT,0,OSEDIM_BEARD)), OPTIONAL, INTENT(IN) :: PEFIELDW ! vertical electric field
!
! optional variables
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PSEA    ! Sea Mask
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(IN)    :: PTOWN   ! Fraction that is town
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL,  INTENT(IN)    :: PCONC3D    ! Cloud droplet number concentration
REAL, DIMENSION(D%NIJT),           OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(D%NIJT,D%NKT,KRR), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQHT ! Hail electric charge at t
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), OPTIONAL, INTENT(INOUT) :: PQHS ! Hail electric charge source
!
! Bridge variables
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PT
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PLVFACT
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PLSFACT
REAL, DIMENSION(D%NIJT,D%NKT,7), INTENT(OUT) :: PWR
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PWTH
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PZCONC3D
LOGICAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: OMICRO
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PRVHENI
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PZZZZ
!
!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
LOGICAL :: LLSEA_AND_TOWN, LLCONC
INTEGER :: JIJ, JK
INTEGER :: IKTB, IKTE, IKB, IIJB, IIJE, IIJT
!
!Arrays for nucleation call outside of LLMICRO points
REAL, DIMENSION(D%NIJT, D%NKT) :: ZSIGMA_RC
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLC_LCF
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLC_LRC
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLI_LCF
REAL, DIMENSION(D%NIJT, D%NKT) :: ZHLI_LRI
!
REAL :: ZINV_TSTEP ! Inverse ov PTSTEP
REAL, DIMENSION(KRR) :: ZICEDRTMIN
!
REAL :: ZDEVIDE, ZRICE, ZZSUM
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZW3D
REAL, DIMENSION(D%NIJT) :: ZCONC_TMP

LOGICAL, DIMENSION(D%NIJT,D%NKT) :: LLW3D
INTEGER :: ISIZE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RAIN_ICE_PART1', 0, ZHOOK_HANDLE)
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
!$acc kernels
ZINV_TSTEP=1./PTSTEP
!
LLSEA_AND_TOWN=PRESENT(PSEA).AND.PRESENT(PTOWN)
LLCONC=PRESENT(PCONC3D)
!
! LSFACT and LVFACT without exner, and LLMICRO
! LLMICRO is a mask with a True value on points where microphysics is active
OMICRO(:,:)=.FALSE.
!$acc end kernels
!$acc kernels
!$acc loop independent collapse(2)
DO JK = IKTB,IKTE
  DO JIJ = IIJB,IIJE
    !LSFACT and LVFACT
    IF (KRR==7) THEN
      ZRICE=PRT(JIJ,JK,IRI)+PRT(JIJ,JK,IRS)+PRT(JIJ,JK,IRG)+PRT(JIJ,JK,IRH)
    ELSE
      ZRICE=PRT(JIJ,JK,IRI)+PRT(JIJ,JK,IRS)+PRT(JIJ,JK,IRG)
    ENDIF
    ZDEVIDE = CST%XCPD + CST%XCPV*PRT(JIJ,JK,IRV) + CST%XCL*(PRT(JIJ,JK,IRC)+PRT(JIJ,JK,IRR)) + CST%XCI*ZRICE
    PT(JIJ,JK) = PTHT(JIJ,JK) * PEXN(JIJ,JK)
    PLSFACT(JIJ,JK)=(CST%XLSTT+(CST%XCPV-CST%XCI)*(PT(JIJ,JK)-CST%XTT)) / ZDEVIDE
    PLVFACT(JIJ,JK)=(CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(JIJ,JK)-CST%XTT)) / ZDEVIDE

    !OMICRO
    IF ( PARAMI%LOCND2 ) THEN
      IF (KRR==7) THEN
        OMICRO(JIJ,JK)=PSSIO(JIJ,JK)>ICEP%XFRMIN(12) .OR. &
                        PRT(JIJ,JK,IRC)>ICEP%XFRMIN(13) .OR. &
                        PRT(JIJ,JK,IRR)>ICEP%XFRMIN(13) .OR. &
                        PRT(JIJ,JK,IRI)>ICEP%XFRMIN(13) .OR. &
                        PRT(JIJ,JK,IRS)>ICEP%XFRMIN(13) .OR. &
                        PRT(JIJ,JK,IRG)>ICEP%XFRMIN(13) .OR. &
                        PRT(JIJ,JK,IRH)>ICEP%XFRMIN(13)
      ELSE
        OMICRO(JIJ,JK)=PSSIO(JIJ,JK)>ICEP%XFRMIN(12) .OR. &
                        PRT(JIJ,JK,IRC)>ICEP%XFRMIN(13) .OR. &
                        PRT(JIJ,JK,IRR)>ICEP%XFRMIN(13) .OR. &
                        PRT(JIJ,JK,IRI)>ICEP%XFRMIN(13) .OR. &
                        PRT(JIJ,JK,IRS)>ICEP%XFRMIN(13) .OR. &
                        PRT(JIJ,JK,IRG)>ICEP%XFRMIN(13)
      ENDIF
    ELSE
      IF (KRR==7) THEN
        OMICRO(JIJ,JK)=PRT(JIJ,JK,IRC)>ICED%XRTMIN(IRC) .OR. &
                        PRT(JIJ,JK,IRR)>ICED%XRTMIN(IRR) .OR. &
                        PRT(JIJ,JK,IRI)>ICED%XRTMIN(IRI) .OR. &
                        PRT(JIJ,JK,IRS)>ICED%XRTMIN(IRS) .OR. &
                        PRT(JIJ,JK,IRG)>ICED%XRTMIN(IRG) .OR. &
                        PRT(JIJ,JK,IRH)>ICED%XRTMIN(IRH)
      ELSE
        OMICRO(JIJ,JK)=PRT(JIJ,JK,IRC)>ICED%XRTMIN(IRC) .OR. &
                        PRT(JIJ,JK,IRR)>ICED%XRTMIN(IRR) .OR. &
                        PRT(JIJ,JK,IRI)>ICED%XRTMIN(IRI) .OR. &
                        PRT(JIJ,JK,IRS)>ICED%XRTMIN(IRS) .OR. &
                        PRT(JIJ,JK,IRG)>ICED%XRTMIN(IRG)
      ENDIF
    ENDIF
  ENDDO
ENDDO
!$acc end kernels
!
!        1.1    Compute cloud liquid number concentration
!               -----------------------------------------
! Model level height
!$acc kernels
!$acc loop independent
DO JIJ=IIJB,IIJE
  ZZSUM = 0.
!$acc loop seq
  DO JK=IKTE,IKTB, -1
    ZZSUM = ZZSUM + PDZZ(JIJ,JK)
    PZZZZ(JIJ,JK) = ZZSUM - PDZZ(JIJ,JK)*0.5
  ENDDO
ENDDO
!$acc end kernels

!$acc kernels
PZCONC3D(:,:) = ICED%XCONC_LAND
IF( LLCONC .AND. PARAMI%LEXCLDROP) THEN
!$acc loop independent collapse(2)
  DO JK=IKTB,IKTE
    DO JIJ=IIJB,IIJE
      PZCONC3D(JIJ,JK)=PCONC3D(JIJ,JK)
    ENDDO
  ENDDO
ELSE
  IF (PARAMI%LEXCLDROP) THEN
   CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'RAIN_ICE_PART1', 'WITH LEXCLDROP=TRUE CLOUD DROPLET FIELDS MUST BE PRESENT IN RAIN_ICE')
  END IF
  IF(LLSEA_AND_TOWN .AND. ICEP%XFRMIN(26)<0.001) THEN
!$acc loop independent collapse(2)
    DO JK=IKTB,IKTE
      DO JIJ=IIJB,IIJE
        ZCONC_TMP(JIJ) = PSEA(JIJ)*ICED%XCONC_SEA+(1.-PSEA(JIJ))*ICED%XCONC_LAND
        PZCONC3D(JIJ,JK)= (1.-PTOWN(JIJ))*ZCONC_TMP(JIJ) +PTOWN(JIJ)*ICED%XCONC_URBAN
      ENDDO
    ENDDO
  ENDIF
  IF(ICEP%XFRMIN(26)>0.001)THEN
!$acc loop independent collapse(2)
    DO JK=IKTB,IKTE
      DO JIJ=IIJB,IIJE
        PZCONC3D(JIJ,JK) = ICEP%XFRMIN(26)*PPABST(JIJ,JK)/CST%XP00
        IF(ICEP%XFRMIN(22)>0.001 .AND. LLSEA_AND_TOWN)THEN
          ZCONC_TMP(JIJ)  = MAX(0., MIN(1.,(PZZZZ(JIJ,JK)- PZZZZ(JIJ,IKTE))/&
              &(MAX(0.01,ICEP%XFRMIN(29)-PZZZZ(JIJ,IKTE)))))
          PZCONC3D(JIJ,JK) =  PZCONC3D(JIJ,JK)*(ICEP%XFRMIN(22)*(PSEA(JIJ)*ICEP%XFRMIN(30)+&
              &(1.-PSEA(JIJ)))*(1.-ZCONC_TMP(JIJ)) + ZCONC_TMP(JIJ))
        ENDIF
      ENDDO
    ENDDO
  ENDIF
ENDIF
!$acc end kernels
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
                         &PLVFACT, PLSFACT, PRHODREF, PPABST, PTHT, PT, PRHODJ, &
                         &PTHS, PRT, PRS, PZCONC3D,&
                         &PINPRC, PINPRR, PINPRS, PINPRG, &
                         &PQCT, PQRT, PQIT, PQST, PQGT, PQCS, PQRS, PQIS, PQSS, PQGS, PEFIELDW, &
                         &TBUDGETS, KBUDGETS, &
                         &PSEA=PSEA, PTOWN=PTOWN, &
                         &PINPRH=PINPRH, PFPR=PFPR, &
                         &PQHT=PQHT, PQHS=PQHS)
ENDIF
!
!
!-------------------------------------------------------------------------------
!
!*       3.     INITIAL VALUES SAVING
!               ---------------------
!
!$acc kernels
DO JK = IKTB,IKTE
!$mnh_expand_array(JIJ=1:D%NIJT)
  !Copy of T variables to keep untouched the prognostic variables
  PWTH(:,JK)=PTHT(:,JK)
  PWR(:,JK,IRV)=PRT(:,JK,IRV)
  PWR(:,JK,IRC)=PRT(:,JK,IRC)
  PWR(:,JK,IRR)=PRT(:,JK,IRR)
  PWR(:,JK,IRI)=PRT(:,JK,IRI)
  PWR(:,JK,IRS)=PRT(:,JK,IRS)
  PWR(:,JK,IRG)=PRT(:,JK,IRG)
  IF (KRR==7) THEN
    PWR(:,JK,IRH)=PRT(:,JK,IRH)
  ELSE
    PWR(:,JK,IRH)=0.
  ENDIF

  !Preset for output 3D variables
  IF(PARAMI%LWARM) THEN
    PEVAP3D(:,JK)=0.
  ENDIF
  PRAINFR(:,JK)=0.
!$mnh_end_expand_array(JIJ=1:D%NIJT)
ENDDO
!$acc end kernels
!
!
!*       4.1    COMPUTES THE SLOW COLD PROCESS SOURCES OUTSIDE OF LLMICRO POINTS
!               -----------------------------------------------------------------
!
!The nucleation must be called everywhere
!This call is for points outside of the OMICRO mask, another call is coded in ice4_tendencies
!$acc kernels
LLW3D(:,:)=.FALSE.
!$acc loop independent collapse(2)
DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE
    IF (.NOT. OMICRO(JIJ, JK)) THEN
      LLW3D(JIJ, JK)=.TRUE.
      ZW3D(JIJ, JK)=PLSFACT(JIJ, JK)/PEXN(JIJ, JK)
      PCIT(JIJ,JK)=0. !ri=0 because where are in the not odmicro case
    ELSE
      LLW3D(JIJ, JK)=.FALSE.
    ENDIF
  ENDDO
ENDDO
!$acc end kernels

!$acc kernels
!$acc loop independent collapse(2)
DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE
!NEC$ noinline
    CALL ICE4_NUCLEATION(CST, PARAMI, ICEP, ICED, LLW3D(JIJ, JK), &
                         PTHT(JIJ, JK), PPABST(JIJ, JK), PRHODREF(JIJ, JK), &
                         PEXN(JIJ, JK), ZW3D(JIJ, JK), PT(JIJ, JK), &
                         PRT(JIJ, JK, IRV), PICLDFR(JIJ, JK), PZZZZ(JIJ, JK), &
                         PCIT(JIJ, JK), PRVHENI(JIJ, JK))
  ENDDO
ENDDO
!$acc end kernels

!$acc kernels
!$acc loop independent collapse(2)
DO JK = IKTB, IKTE
  DO JIJ=IIJB, IIJE
    PRVHENI(JIJ,JK) = MIN(PRS(JIJ,JK,IRV), PRVHENI(JIJ,JK)/PTSTEP)
  ENDDO
ENDDO
!$acc end kernels
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
        ZHLC_LRC(JIJ, JK) = PWR(JIJ, JK, IRC) - PHLC_HRC(JIJ, JK)
        ZHLI_LRI(JIJ, JK) = PWR(JIJ, JK, IRI) - PHLI_HRI(JIJ, JK)
        IF(PWR(JIJ, JK, IRC)>0.) THEN
          ZHLC_LCF(JIJ, JK) = PCLDFR(JIJ, JK)- PHLC_HCF(JIJ, JK)
        ELSE
          ZHLC_LCF(JIJ, JK)=0.
        ENDIF
        IF(PWR(JIJ, JK, IRI)>0.) THEN
          ZHLI_LCF(JIJ, JK) = PCLDFR(JIJ, JK)- PHLI_HCF(JIJ, JK)
        ELSE
          ZHLI_LCF(JIJ, JK)=0.
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  !We cannot use PWR(:,IKTB:IKTE,IRC) which is not contiguous
  ISIZE=IIJT*(IKTE-IKTB+1)
  CALL ICE4_COMPUTE_PDF(CST, ICEP, ICED, ISIZE, PARAMI%CSUBG_AUCV_RC, PARAMI%CSUBG_AUCV_RI, PARAMI%CSUBG_PR_PDF,&
                        OMICRO(:,IKTB:IKTE), PRHODREF(:,IKTB:IKTE), PRT(:,IKTB:IKTE,IRC), PRT(:,IKTB:IKTE,IRI), &
                        PCLDFR(:,IKTB:IKTE), PT(:,IKTB:IKTE), ZSIGMA_RC(:,IKTB:IKTE), &
                        PHLC_HCF(:,IKTB:IKTE), ZHLC_LCF(:,IKTB:IKTE), PHLC_HRC(:,IKTB:IKTE), ZHLC_LRC(:,IKTB:IKTE), &
                        PHLI_HCF(:,IKTB:IKTE), ZHLI_LCF(:,IKTB:IKTE), PHLI_HRI(:,IKTB:IKTE), ZHLI_LRI(:,IKTB:IKTE), &
                        PRAINFR(:,IKTB:IKTE))
  IF (KRR==7) THEN
    CALL ICE4_RAINFR_VERT(D, ICED, PRAINFR, PWR(:,:,IRR), &
                         &PWR(:,:,IRS), PWR(:,:,IRG), PWR(:,:,IRH))
  ELSE
    CALL ICE4_RAINFR_VERT(D, ICED, PRAINFR, PWR(:,:,IRR), &
                         &PWR(:,:,IRS), PWR(:,:,IRG))
  ENDIF
ELSE
  PRAINFR(:,:)=1.
ENDIF
!
IF (LHOOK) CALL DR_HOOK('RAIN_ICE_PART1', 1, ZHOOK_HANDLE)
!
CONTAINS
INCLUDE "ice4_nucleation.func.h"
END SUBROUTINE RAIN_ICE_PART1
