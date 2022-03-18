!OPTION! -pvctl noloopfusion
SUBROUTINE APL_AROME(YDGEOMETRY,YDSURF, YDCFU, YDXFU, YDMODEL, KBL, KGPCOMP, KIDIA , KFDIA , KLON ,&
 & KTDIA  , KLEV , KSTEP ,&
 & KMAXDRAFT, KSGST, KNFRRC, PDT, LDXFUMSE, PINDX, PINDY ,&
 & PGEMU,PGELAM,POROG,PGM,PMU0,PMU0LU,PMU0M,PMU0N,PCLON, PSLON,PVO3ABC,PLSM,&
 & PAESEA , PAELAN , PAESOO , PAEDES , PAESUL, PAEVOL,&
 & PGP2DSDT, PGP2DSPP, &
 !---------------------------------------------------------------------
 ! - INPUT A M
 & PAPHIM,PAPHIFM,PAPRSM, PAPRSFM, PRDELPM, PDELPM, PTM, PQVM ,&
 & PCPM    , PRM     ,PALPHM , PLNPRM,&
 & PQCM    , PQIM    ,PQRM   , PQSM, PQGM, PQHM,&
 
 & PLIMAM  , &
 & PTKEM   , PEFB1   ,PEFB2  , PEFB3,&
 & PSIGM,PSVM,&
 & PUM    , PVM, PWM, PEDR,&
 & PFORCEU,PFORCEV,PFORCET,PFORCEQ,&
 !---------------------------------------------------------------------
 !  - INOUT A S
 & PGPAR, PEMTD, PEMTU, PTRSO,&
 & PGDEOSI, PGUEOSI, PGMU0, PGMU0_MIN, PGMU0_MAX,&
 & PGDEOTI, PGDEOTI2, PGUEOTI, PGUEOTI2, PGEOLT, PGEOXT,&
 & PGRPROX, PGMIXP, PGFLUXC, PGRSURF,&
 & PTURB3D,&
 !  - OUT A S
 & PQLRAD, PQIRAD, PRH, PCLFS, PSIGS,&
 & PTENDT, PTENDR, PTENDU, PTENDV,PTENDW,&

 & PTENDLIMA, &
 & PTENDTKE, PTENDEFB1, PTENDEFB2, PTENDEFB3,&
 & PTENDEXT,PFRTH, PFRSO,PFRTHDS, PFRSODS, PFRSOPS, PFRSDNI,&
 & PFRSOPT, PFRTHC, PFRSOC, &
 !---------------------------------------------------------------------
 !  - IN FOR RADIATION IF NO SURFACE SCHEME
 & PALBIN,PEMIS,&
 !  - INOUT for easy diag
 & PEZDIAG,&
 !  - INOUT for CFU XFU
 & PCLCH,PCLCL,PCLCM,PCLCT,PFPLSL,PFPLSN,PFPLSG,PFPLSH,PSTRTU,PSTRTV,PFCS,PFCLL,&
 & PFCLN,PUCLS,PVCLS,PNUCLS,PNVCLS,PTCLS,PQCLS,PHUCLS,PUGST,PVGST,PFEVL,PFEVN, PPBLH,PSPSG,PSPSGR,&
 & PSDUR,PDIAGH,PFLASH,PSFORC,PTPWCLS,PDPRECIPS,PDPRECIPS2,PVISICLD,PVISIHYDRO,PMXCLWC,&
 ! daand: radflex
 & YDPROCSET ,YDDDH)

!**** *APL_AROME * - CALL OF PHYSICAL PARAMETERISATIONS FOR ALARO/AROME

!     Sujet.
!     ------
!     - APPEL DES SOUS-PROGRAMMES DE PARAMETRISATION

!**   Interface.
!     ----------
!        *CALL* *APL_AROME*

!-----------------------------------------------------------------------

! -   ARGUMENTS D'ENTREE.
! -   INPUT ARGUMENTS.
!     -------------------

! - NOM DES PARAMETRES DE DIMENSIONNEMENT DE LA PHYSIQUE.
! - DIMENSIONS.

! KBL  : NUMERO DE BLOC NPROMA
! KBL  : NPROMA-PACKETS NUMBER
! KGPCOMP : NOMBRE TOTAL DE POINTS DE GRILLE SUR LE DOMAINE
! KGPCOMP : TOTAL GRID POINTS NUMBER IN THE DOMAIN
! KIDIA, KFDIA : BORNES BOUCLES HORIZONTALES   (IST,IEND DANS CPG).
! KIDIA, KFDIA : START/END OF HORIZONTAL LOOP  (IST,IEND IN *CPG*).
! KLON : DIMENSION HORIZONTALE                 (NPROMA DANS CPG).
! KLON : HORIZONTAL DIMENSION                  (NPROMA IN *CPG*).
! KTDIA : DEBUT BOUCLE VERTICALE DANS LA PHYSIQUE.
! KTDIA : START OF THE VERTICAL LOOP IN THE PHYSICS (IF SOME LEVELS ARE
!                     SKIPPED AT THE TOP OF THE MODEL).
! KLEV : FIN BOUCLE VERTICE ET DIMENSION VERTICALE (NFLEVG DANS CPG).
! KLEV : END OF VERTICAL LOOP AND VERTICAL DIMENSION(NFLEVG IN *CPG*).
! KSTEP : TIME STEP NUMBER (starting with zero)
! KMAXDRAFT : MAX NUMBER OF DRAFTS (FOR DIMENSIONNING)
! KSGST : NUMBER OF SUBGRID SURFACE TEMPERATURES AND FLUXES (NTSSG IN *CPG*)
! KNFRRC : FREQUENCY FOR CLEAR SKY RADIATION CALCULATION
! PDT : TIME STEP (in s) 
! LDXFUMSE : T if CDCONF=X in order not to increment surfex timer in that case
!-----------------------------------------------------------------------
! PGEMU      : SINE OF GEOGRAPHICAL LATITUDE
! PGELAM     :  LONGITUDE
! POROG      : g * OROGRAPHY
! PGM        : MAP FACTOR (used in ALARO convection only)
! PMU0       : COSINUS LOCAL INSTANTANE DE L'ANGLE ZENITHAL SOLAIRE.
! PMU0LU     : COSINUS LOCAL INSTANTANE DE L'ANGLE ZENITHAL LUNAIRE.
! PMU0       : LOCAL COSINE OF INSTANTANEOUS SOLAR ZENITH ANGLE.
! PMU0M      : COSINUS LOCAL MOYEN DE L'ANGLE ZENITHAL.
! PMU0M      : LOCAL COSINE OF AVERAGED SOLAR ZENITH ANGLE.
! PMU0N      : COSINUS LOCAL AU PAS DE TEMPS SUIVANT DE L'ANGLE ZENITHAL SOLAIRE.
! PMU0N      : NEXT TIME STEP COSINUS LOCAL INSTANTANE DE L'ANGLE ZENITHAL SOLAIRE.
! PCLON      : cosine of geographical longitude.
! PSLON      : sine of geographical longitude.
! PVO3ABC    : OZONE COEFFICIENTS FOR ARPEGE PROFILES
! PLSM       : -ATMOSPHERIC MODEL- LAND-SEA MASK (! MAY BE DIFFERENT FROM 
!              THE SURFACE ONE) 
! PAESEA     : MARINE AEROSOLS (IF NVCLIA >= 4)
! PAELAN     : CONTINENTAL AEROSOLS (IF NVCLIA >= 4)
! PAESOO     : SOOT AEROSOLS (IF NVCLIA >= 4)
! PAEDES     : DESERT AEROSOLS (IF NVCLIA >= 4)
! PAESUL     : SULFATE AEROSOLS  (IF LAEROSUL=.T.)
! PAEVOL     : VOLCANO AEROSOLS  (IF LAEROVOL=.T.)
! PGP2DSDT   : STOCHASTIC PHYSICS PATTERNS

! FIELDS WITH SUBSCRIPT M FOR TIME T-DT IN 3TL OR T IN 2TL

! PAPHIM     : GEOPOTENTIAL ON HALF-LEVELS
! PAPHIFM    : GEOPOTENTIAL ON FULL-LEVELS
! PAPRSM     : PRESSURE ON HALF LEVELS
! PAPRSFM    : PRESSURE ON FULL LEVELS.
! PRDELPM    : INVERSE OF PDELP 
! PDELPM     : LAYER THICKNESS IN PRESSURE UNITS

! PTM        : TEMPERATURE.
! PQVM        : SPECIFIC HUMIDITY OF WATER VAPOR
! PCPM        : SPECIFIC HEAT AT CONSTANT PRESSURE FOR AIR
! PRM         : GAS CONSTANT FOR AIR
! PALPHM      : "alpha" on layers
! PLNPRM      : "delta" on layers

! PQCM        : SPECIFIC HUMIDITY OF CLOUD WATER
! PQIM        : SPECIFIC HUMIDITY OF ICE
! PQRM        : SPECIFIC HUMIDITY OF RAIN
! PQSM        : SPECIFIC HUMIDITY OF SNOW
! PQGM        : SPECIFIC HUMIDITY OF GRAUPEL
! PQHM        : SPECIFIC HUMIDITY OF HAIL
! PTKEM       : TURBULENT KINETIC ENERGY
! PSVM        : PASSIVE SCALARS
! PSIGM       : SIGMA FOR SUBGRIDCOND
! PUM         : ZONAL WIND
! PVM         : MERIDIAN WIND
! PWM         : VERTICAL VELOCITY (m/s)

!-----------------------------------------------------------------------
! - INOUT

! PGPAR       : BUFFER FOR 2D FIELDS - CONTAINS PRECIP, ALBEDO, EMISS, TS 
!             : SURFACE FLUXES
! PEMTD       : DOWNWARD LONGWAVE EMISSIVITY
! PEMTU       : UPWARD   LONGWAVE EMISSIVITY
! PTRSO       : SHORTWAVE TRANSMISSIVITY

! ACRANEB2 intermittency storage

! PGDEOSI   : DESCENDING INCREMENTAL OPTICAL DEPTHS, SOLAR
! PGUEOSI   : ASCENDING  INCREMENTAL OPTICAL DEPTHS, SOLAR
! PGMU0     : COSINE OF SOLAR ZENITH ANGLE, APPROXIMATE ACTUAL VALUE
! PGMU0_MIN : COSINE OF SOLAR ZENITH ANGLE, MIN VALUE
! PGMU0_MAX : COSINE OF SOLAR ZENITH ANGLE, MAX VALUE
! PGDEOTI     : DESCENDING INCREMENTAL OPTICAL DEPTHS, dB/dT(T0) WEIGHTS
! PGDEOTI2    : DESCENDING INCREMENTAL OPTICAL DEPTHS, B WEIGHTS WITH
!               LINEAR T_e CORRECTION
! PGUEOTI     : ASCENDING INCREMENTAL OPTICAL DEPTHS, dB/dT(T0) WEIGHTS
! PGUEOTI2    : ASCENDING INCREMENTAL OPTICAL DEPTHS, B WEIGHTS WITH
!               LINEAR T_e CORRECTION
! PGEOLT      : LOCAL OPTICAL DEPTHS, dB/dT(T0) WEIGHTS
! PGEOXT      : MAXIMUM OPTICAL DEPTHS FOR EBL-EAL, dB/dT(T0) WEIGHTS
! PGRPROX     : CORRECTION TERM FOR ADJACENT EXCHANGES
! PGMIXP      : NON-STATISTICAL WEIGHTS FOR BRACKETING
! PGFLUXC     : OUT OF BRACKET PART OF CLEARSKY EBL, RESP. EBL-EAL FLUX
! PGRSURF     : CORRECTIVE RATIO FOR SURFACE CTS CONTRIBUTION
! PTURB3D     : MATRICE DE GRADIENTS HORIZONTAUX
!-----------------------------------------------------------------------
! - OUTPUT (SUBSCRIPT S FOR T+DT)

! PCLFS       : CLOUD FRACTION
! PQLRAD      : SPECIFIC HUMIDITY OF CLOUD WATER FOR RTTOV
! PQIRAD      : SPECIFIC HUMIDITY OF ICE FOR RTTOV
! PSIGS       : SIGMA FOR SUBGRIDCOND
! PTENDT      : TEMPERATURE TENDENCY
! PTENDR      : HYDROMETEORE TENDENCIES
! PTENDU      : ZONAL WIND TENDENCY
! PTENDV      : MERIDIAN WIND TENDENCY
! PTENDW      : VERTICAL VELOCITY TENDENCY
! PTENDTKE    : TKE TENDENCY
! PTENDEXT    : PASSIVE SCALARS TENDENCY
! PFRTH       : LONGWAVE RADIATIVE FLUX
! PFRSO       : SHORTWAVE RADIATIVE FLUX
! PFRTHDS     : LONGWAVE DOWNWARD SURFACE RADIATIVE FLUX
! PFRSOPS     : SHORTWAVE DOWNWARD SURFACE RADIATIVE FLUX DIRECT
! PFRSDNI     : SHORTWAVE DIRECT NORMAL IRRADIANCE
! PFRSODS     : SHORTWAVE DOWNWARD SURFACE RADIATIVE FLUX GLOBAL
! PFRSOPT     : SHORTWAVE DOWNWARD TOP RADIATIVE FLUX DIRECT
! - 2D (0:1)
! PFRTHC      : LONGWAVE CLEAR SKY NET RADIATIVE FLUX
! PFRSOC      : SHORTWAVE CLEAR SKY NET RADIATIVE FLUX

! variables used in input for radiation in case no surface scheme is used 

! PALBIN     : MODEL SURFACE SHORTWAVE ALBEDO.
! PEMIS      : MODEL SURFACE LONGWAVE EMISSIVITY.

! Part of GFL strcture dedicated to easy diagnostics (to be used as a print...)
! PEZDIAG    : MULPITPLE ARRAY TO BE FILLED BY THE USER BY 3D FIELDS
!              (NGFL_EZDIAG ONES)
! output for CFU XFU
! PCLCH      : HIGH CLOUD COVER (DIAGNOSTIC).
! PCLCL      : LOW CLOUD COVER (DIAGNOSTIC).
! PCLCM      : MEDIUM CLOUD COVER (DIAGNOSTIC).
! PCLCT      : TOTAL CLOUD COVER (DIAGNOSTIC).
! PFPLSL     : RESOLVED PRECIPITATION AS RAIN.
! PFPLSN     : RESOLVED PRECIPITATION AS SNOW
! PFPLSG     : RESOLVED PRECIPITATION AS GRAUPEL
! PFPLSH     : RESOLVED PRECIPITATION AS HAIL
! PSTRTU     : TURBULENT FLUX OF MOMENTUM "U".
! PSTRTV     : TURBULENT FLUX OF MOMENTUM "V".
! PFCS       : SENSIBLE HEAT FLUX AT SURFACE LEVEL.
! PFCLL      : LATENT HEAT FLUX AT SURFACE LEVEL OVER WATER.
! PFCLN      : LATENT HEAT FLUX AT SURFACE LEVEL OVER SNOW.
! PUCLS      : SORTIE DIAGNOSTIQUE DU VENT EN X A HUV METEO.
! PUCLS      : U-COMPONENT OF WIND AT 10 METERS (DIAGNOSTIC).
! PVCLS      : SORTIE DIAGNOSTIQUE DU VENT EN Y A HUV METEO.
! PVCLS      : V-COMPONENT OF WIND AT 10 METERS (DIAGNOSTIC).
! PNUCLS     : SORTIE DIAGNOSTIQUE DU VENT NEUTRE EN X A HUV METEO.
! PNUCLS     : U-COMPONENT OF NEUTRAL WIND AT 10 METERS (DIAGNOSTIC).
! PNVCLS     : SORTIE DIAGNOSTIQUE DU VENT NEUTRE EN Y A HUV METEO.
! PNVCLS     : V-COMPONENT OF NEUTRAL WIND AT 10 METERS (DIAGNOSTIC).
! PTCLS      : SORTIE DIAGNOSTIQUE DE LA TEMPERATURE A HTQ METEO.
! PTCLS      : TEMPERATURE AT 2 METERS (DIAGNOSTIC).
! PQCLS      : SORTIE DIAGNOSTIQUE DE L'HUMIDITE SPECIFIQUE A HTQ METEO.
! PQCLS      : SPECIFIC HUMIDITY AT 2 METERS (DIAGNOSTIC).
! PHUCLS     : SORTIE DIAGNOSTIQUE DE L'HUMIDITE RELATIVE A HTQ METEO.
! PHUCLS     : RELATIVE HUMIDITY AT 2 METERS (DIAGNOSTIC).
! PUGST      : SORTIE DIAGNOSTIQUE DU VENT RAFALE EN X A HUV METEO.
! PUGST      : U-COMPONENT OF WIND GUST AT 10 METERS (DIAGNOSTIC).
! PVGST      : SORTIE DIAGNOSTIQUE DU VENT RAFALE EN Y A HUV METEO.
! PVGST      : V-COMPONENT OF WIND GUST AT 10 METERS (DIAGNOSTIC).
! PFEVL      : FLUX DE VAPEUR D'EAU SUR EAU LIQUIDE (OU SOL HUMIDE).
! PFEVL      : WATER VAPOUR FLUX OVER LIQUID WATER (OR WET SOIL)
! PFEVN      : FLUX DE VAPEUR D'EAU SUR NEIGE (OU GLACE) ET SOL GELE.
! PFEVN      : WATER VAPOUR FLUX OVER SNOW (OR ICE) AND FROZEN SOIL.
! PPBLH      : PSEUDO-HISTORICAL ARRAY FOR PBL HEIGHT
! PSPSG      : SNOW COVER
! PSPSGR     : SNOW DENSITY
! PSDUR      : SUNSHINE DURATION [s]
! PDIAGH     : HAIL DIAGNOSTIC
! PFLASH     : LIGHTNING DIAGNOSTICS
! PVISICLD   : VISIBILITY DUE TO CLOUD WATER AND CLOUD ICE
! PVISIHYDRO : VISIBILITY DUE TO RAIN AND SNOW AND GRAUPEL
! PMXCLWC    : CLOUD WATER LIQUID CONTENT AT HVISI METERS
! PDPRECIPS  : PRECIPITATION TYPE
! PDPRECIPS2 : PRECIPITATION TYPE FOR 2NDE PERIOD

!-----------------------------------------------------------------------

!     Externes.
!     ---------

!     Method
!     ------
!     - convert aladin variables into mesonh variables (level inversion 
!       and q to r, t to theta) 
!     - call mesoNH physics and ECMWF radiation scheme
!     - convert mesoNH tendencies to aladin tendencies

!     Auteur.
!     -------
!      S.Malardel et Y. Seity
!      10-03-03
!      big cleaning (18/06/04) S. Malardel and Y. Seity
!     externalisation of surface scheme call + small cleaning (20-07-04) Y.Seity
!     Modifications
!     -------------
!      G. Hello 04-02-06: Add the call of KFB-convection scheme 
!                         for future use in ALARO
!      T.Kovacic 04-05-05: Added ZCVTENDPR_ and ZCVTENDPRS_ 
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!      F.Bouyssel 04-05-05: New arguments in ACRADIN
!     Y. Seity 30-Sept-2005 Add MNH Chemistry scheme
!     R. Zaaboul 15-feb-2006 add surface scheme call
!     T.Kovacic  2006-03-23: calls to subroutines for budgets 
!                             and new arguments PFRTH and PFRSO
!     Y. Seity   2007-05-07: add CFU and XFU calculations 
!                           and call aro_ground_diag
!     S.Ivatek-S 2007-04-17: Over dimensioning of PGPAR by NGPAR+1 just 
!                            (KLON,NGPAR) is used boundary checking bf
!     T.Kovacic  2007-03-16: Fourth dim. in APFT
!     JJMorcrette, ECMWF, 20080325: dummy arguments for RADACT to allow for 
!                        using a new sulphate climatology in the ECMWF model
!     Y. Seity   2008-06-15: correct calculations of PFRTHDS, PFRSODS and PFCLL
!     Y. Seity   2008-09-29: phasing Chemistry corrections
!     O.Riviere  2008-10-01: introduction of new data flow for DDH in Arome
!     Y. Seity   2009-05-03: new version of EDKF and implementation of EDMF
!     Y. Seity   2009-10-03: add missed deallocations 
!     S. Riette  2009-03-25: Arguments modification for AROCLDIA to add HTKERAF
!     K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!     A. Alias   2009-09-01: Sulfate and Volcano aerosols added (call radaer)
!     S. Riette  2010-01-19: ZUM__, ZVM__ and ZDEPTH_HEIGHT_ are given
!                            ARO_GROUND_DIAG in 3D.                     
!     Y. Seity   2010-03-09: add PFEVN and PFEVL
!     Y. Bouteloup 2010-03-26 : Add PQLRAD et PQIRAD
!     Y. Seity : Test TKE > 0.
!     Y. Seity : Optimized version of EDKF + diag HCLS
!     Y. Seity : 2010-09 Save Ts at the end of apl_arome for ICMSH+0000
!     L. Bengtsson (2010): Introduce cloud diagnostics based on geop. 
!                               height (LWMOCLOUD), AND cloud-overlap assumptions 
!                               from C. Wittman 2009 (LACPANMX + WMXOV)
!     S. Riette: 2010-12 aro_ground_diag interface modified
!     Y. Seity: 2010-12 add hail diagnostic
!     R. El Khatib 30-Jun-2010 NEC directive noloopfusion to preserve critical regions
!     P.Marguinaud 2010-06-29 : KSURFEXCTL flag (disable SURFEX) 
!     2010-12    B. Decharme  : modify the radiative coupling with surfex (SW per band in ACRADIN and RADHEAT)
!     2011-02    A. Voldoire : add ZAERINDS to CALL RADAER and ACRADIN
!                              for sulfate indirect effect computation
!     2011-06: M. Jerczynski - some cleaning to meet norms
!     S. Riette: 2011-10 : Modifications for DUAL-MF scheme (according to Wim de Rooy's apl_arome version)
!                          Ice in EDKF
!     Y. Seity : 2012-03 : add LMNHLEV option to revert/or not arrays for MesoNH parameterisations
!     F. Bouttier: 2012-07 add SPPT stochastic physics
!     JJMorcrette, ECMWF, 20120815 additional dummy due to changes in RADACT
!     P. Marguinaud : 2012-09 : Add control threshold for orography
!     Y. Seity : 2013-01 Cleaning LMNHLEV and remove JPVEXT points
!     Y. Seity : 2013-02 Cleaning (add compute_neb)
!     L. Bengtsson: 2013-02: add LOLSMC and LOTOWNC options to compute (or not) cloud sedimentation
!                            using different cloud droplet number conc. depending on land/sea/town.
!     2013-11, D. Degrauwe: Introduction of radflex interface, export
!                           upper-air precipitation fluxes PFPR.
!     2013-11, J. Masek: Inclusion of ACRANEB2 radiation scheme.
!     S. Riette: 2013-11: subgrid precipitation
!     K. Yessad (July 2014): Move some variables.
!     2014-09, C. Wastl: Adaptations for orographic shadowing
!     2014-11, Y. Seity: add TKE budgets for DDH
!     2016-03, E. Bazile: Phasing MUSC for surf_ideal_flux
!     2016-04, J. Masek: LRNUEXP cloud overlap option (COMPUTE_NEB replaced
!                        by ACNPART), passing of sushine duration, fix of
!                        E. Gleeson for ACRANEB2 with SURFEX.
!     2016-09, J. Masek: Proper calculation of sunshine duration in ACRANEB2.
!     2016-10, P. Marguinaud : Port to single precision
!     S. Riette 2016-11: Changes in ICE3/ICE4
!     2018-09, E. Gleeson: Corrected misplaced arguments in ACRANEB2 call.
!     2019-09-24 J.M. Piriou arguments for convective gusts.
!     R. El Khatib 30-Oct-2018 substantial rewrite for optimization and coding standards respect.
!     2018-10, I. Etchevers : add Visibilities
!     2019-01, I. Etchevers, Y. Seity : add Precipitation Type
!     2019-09, J. Masek: Corrected dimensioning of dummy argument PGMU0.
!                        Modified call to ACRANEB2 (clearsky fluxes).
!     2019-10, I. Etchevers : Visibilities in ACVISIH, AROCLDIA=>ACCLDIA
!     2019-10, Y.Bouteloup and M. Bouzghaiam : Radiation modifications. Remove acradin.F90 direct
!              call to recmwf.F90 and add interface to ecrad (in recmwf !)
!     2020-10, J. Masek: Modified call to ACCLDIA.
!     2020-12, F. Meier add call to latent heat nudging if LNUDGLH is TRUE
!     2020-12, U. Andrae : Introduce SPP for HARMONIE-AROME
!     R. El Khatib 24-Aug-2021 NPROMICRO specific cache-blocking factor for microphysics
! End modifications
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMXFU             , ONLY : TXFU
USE YOMCFU             , ONLY : TCFU
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPIM     ,JPRB     ,JPRD
USE YOMHOOK            , ONLY : LHOOK    ,DR_HOOK

! AROME SPECIFIC


! OTHERS
USE YOESW      , ONLY : RSUN2
USE YOMCST     , ONLY : RG       ,RCPD     ,RD       ,RATM     ,RTT      ,&
          &             RCW      ,RCPV     ,RLVTT    ,RCS      ,RLSTT    ,RGAMW    ,&
          &             RBETW    ,RALPW    ,RGAMS    ,RBETS    ,RALPS    ,RGAMD    ,&
          &             RBETD    ,RALPD    ,RETV     ,RV       ,RKAPPA   ,RHOUR
USE YOMLUN     , ONLY : NULOUT
USE YOMCT0     , ONLY : LTWOTL, LSFORCS
USE YOMVERT    , ONLY : VP00
USE YOMRIP0    , ONLY : NINDAT
USE YOMNUDGLH , ONLY :  LNUDGLH, NSTARTNUDGLH, NSTOPNUDGLH, NINTNUDGLH, NTAUNUDGLH, &
          &             RAMPLIFY,RMAXNUDGLH,RMINNUDGLH,LNUDGLHCOMPT,NTIMESPLITNUDGLH
USE YOMNSV     , ONLY : NSV_CO2
USE DDH_MIX    , ONLY : ADD_FIELD_3D, NEW_ADD_FIELD_3D, TYP_DDH ! for new diag data flow
USE YOMSPSDT   , ONLY : YSPPT_CONFIG, YSPPT
USE SPP_MOD    , ONLY : YSPP_CONFIG, YSPP
USE YOMLSFORC  , ONLY : LMUSCLFA, NMUSCLFA, REMIS_FORC, RALB_FORC
! daand: radflex
USE INTFLEX_MOD, ONLY : LINTFLEX, LRADFLEX,&
                      & TYPE_INTPROC, TYPE_INTPROCSET,&
                      & NEWINTFIELD, NEWINTPROC
USE YOMMP0     , ONLY : MYPROC     

!     --------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TXFU)        ,INTENT(INOUT) :: YDXFU
TYPE(TCFU)        ,INTENT(INOUT) :: YDCFU
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KBL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPCOMP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMAXDRAFT
INTEGER(KIND=JPIM),INTENT(IN)    :: KSGST
INTEGER(KIND=JPIM),INTENT(IN)    :: KNFRRC
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
LOGICAL           ,INTENT(IN)    :: LDXFUMSE 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINDX(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINDY(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0LU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0M(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0N(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLON(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLON(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVO3ABC(KLON,3) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAESEA(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAELAN(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAESOO(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAEDES(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAESUL(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAEVOL(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSDT(KLON,YSPPT%YGPSDT(1)%NG2D,YSPPT%N2D)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGP2DSPP(KLON,YSPP%N2D)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIM(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSM(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSFM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELPM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELPM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQVM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCPM(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPHM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPRM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQCM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQIM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQRM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQGM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQHM(KLON,KLEV) 

REAL(KIND=JPRB)   ,INTENT(IN)    :: PLIMAM(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NLIMA)
REAL(KIND=JPRB)   ,INTENT(IN), TARGET :: PTKEM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEFB1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEFB2(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEFB3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSVM(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIGM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORCEU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORCEV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORCET(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORCEQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGPAR(KLON,YDMODEL%YRML_PHY_MF%YRPARAR%NGPAR+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEMTD(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEMTU(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTRSO(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGDEOSI(KLON,0:KLEV,2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGUEOSI(KLON,0:KLEV,2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMU0(KLON,0:YDMODEL%YRML_PHY_MF%YRPHY%NSORAYFR-1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMU0_MIN(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMU0_MAX(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGDEOTI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGDEOTI2(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGUEOTI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGUEOTI2(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGEOLT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGEOXT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGRPROX(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMIXP(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLUXC(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGRSURF(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTURB3D(KLON,YDMODEL%YRML_PHY_MF%YRARPHY%NGRADIENTS,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCLFS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQLRAD(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEDR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQIRAD(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSIGS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDR(KLON,KLEV,YDMODEL%YRML_PHY_MF%YRPARAR%NRR) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDW(KLON,KLEV) 
 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDLIMA(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NLIMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDTKE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRH(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENDEFB1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENDEFB2(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENDEFB3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDEXT(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBIN(KLON)
! daand: radflex; made target
REAL(KIND=JPRB)   ,INTENT(INOUT), TARGET :: PFRTH(KLON,0:KLEV,KSGST+1)
REAL(KIND=JPRB)   ,INTENT(OUT), TARGET :: PFRSO(KLON,0:KLEV,KSGST+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRTHDS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSOPS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSDNI(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSODS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSOPT(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSOC(KLON,0:1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRTHC(KLON,0:1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEZDIAG(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCH(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCM(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCT(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLSL(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLSN(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLSG(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLSH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRTU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRTV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCS(KLON,KSGST+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCLL(KLON,KSGST+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCLN(KLON,KSGST+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEVL(KLON,KSGST+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEVN(KLON,KSGST+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUCLS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVCLS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNUCLS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNVCLS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTCLS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQCLS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHUCLS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUGST(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVGST(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PPBLH(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPSG(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPSGR(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSDUR(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIAGH(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTPWCLS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLASH(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDPRECIPS(KLON,YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS%NDTPREC)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDPRECIPS2(KLON,YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS%NDTPREC2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVISICLD(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVISIHYDRO(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMXCLWC(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSFORC(KLON,YDMODEL%YRML_PHY_MF%YRPHYDS%NSFORC)
! daand: radflex
TYPE(TYPE_INTPROCSET), INTENT(INOUT) :: YDPROCSET
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH
!*
!     ------------------------------------------------------------------

! 3D arrays de reference dans mesoNH. En 1D, thetavref=thetavM, mais la question
! concernant la facon d initialiser cette variable dans le 3D reste ouverte (idem pour RHODREF)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                                READ ME, PLEASE !

!                         CODING CONVENTIONS FOR ARPEGE vs MNH PHYSICS

! The horizontal representation in MNH physics is such that 2 dimensions are needed, while a single dimension
! is used in ARPEGE and AROME. The remapping from 1 to 2 dimensions will be made inplicitly through the 
! subroutines interfaces (this is a fortran property). Therefore there is non need to add an explicit dimension
! sized to 1.

! Local 3D arrays with extra levels for Meso-NH turbulence scheme :
! - first dimension is KFDIA not KDLON in order to limit array copies
! - suffixed with two underscore to be easily identified
! These arrays are passed in argument as ZXXX__(:,1:KLEV) except for aro_turb_mnh where they are passed as ZXXX__.

! Local 3D arrays with regular number of levels for Meso-NH interfaces :
! - first dimension is KFDIA not KDLON in order to limit array copies due to non-contiguous data.
! - suffixed with one underscore to be easily identified.

! Local 4D arrays with regular number of levels for Meso-NH interfaces :
! - first dimension is KFDIA not KDLON in order to limit array copies
! - suffixed with one underscore to be easily identified

! Local 2D arrays for Meso-NH interfaces :
! - first dimension is KFDIA not KDLON as a convention as well (no arrays copies to fear)
! - suffixed with one underscore to be easily identified

! Other arrays, which can be dummy arguments, or local but used as argument to IFS/ARPEGE physics
! should remained dimensionned KDLON and should not be suffixed with undersores.

!                         DO NOT USE ARRAY SYNTAX FOR COMPUTATIONAL LOOPS !!

! - They make the code less performant because memory cache is poorly used
! - They can make the code even less readable if the indexes are removed

!                         AVOID ARRAYS COPIES, OR MAKE THEM FAST !

! - if you do need to initialize or copy an array, do it as follows with explicit array syntax in first dimension
! because the compiler will be able to use an optimized function to initialize/copy a segment of memory,
! and may be able to address simultaneously several cach lines :
! 1D array : 
!   Z(KIDIA:KFDIA)=value
! 2D arrays :
!  DO JLEV=1,KLEV
!     ZX(KIDIA:KFDIA,JLEV)=xval
!     ZY(KIDIA:KFDIA,JLEV)=yval
!  ENDDO

! - if you need the bakup of an array, use a swapp mechanism, as what is done here for instance for
! ZSVM_ and ZSVMIN_ : 
!  for the developer ZSVMIN_ is always the bakup and ZSVM_ is always the current array.
!  ZSVMSWAP_ or ZSVMSAVE_ are used for the swapp mechanism, they should never be used by the developer.

! - do not initialize if not necessary. To avoid useless initialization use the mechanizm 'INIT0' coded below :
! IF (INIT0 == 0) THEN
!   ZVALUE=HUGE(1._JPRB)
! ELSE
!   ZVALUE=default_value
! ENDIF
! IF (INIT0 >= 0) THEN
!   Z(:)=ZVALUE
! ENDIF
! INIT0= 0 : initialize to HUGE (testing/debugging)
! INIT0= 1 : initialize to realistic value (discouraged !)
! INIT0=-1 : no initialization (optimized code) - this is the default.

! - pointers can advantageously be used to avoid copies or to avoid the initialization to zero of a sum.
! Example :
!  DO JI=IFIRST,ILAST
!    IF (JI == IFIRST) THEN
!      ! Fill the sum at the first iteration
!      ZARG => ZSUM(:)
!    ELSE
!      ! increment
!      ZARG => ZINC(:)
!    ENDIF
!    CALL COMPUTE(ZARG)
!    IF (JI > IFIRST) THEN
!      ! Add increment
!      ZSUM(KIDIA:KFDIA)=ZSUM(KIDIA:KFDIA)+ZINC(KIDIA:KFDIA)
!    ENDIF
!  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(KIND=JPRB) :: ZRHODJM__(KFDIA,0:KLEV+1),      ZRHODREFM__(KFDIA,0:KLEV+1),   ZPABSM__(KFDIA,0:KLEV+1)
REAL(KIND=JPRB) :: ZUM__(KFDIA,0:KLEV+1),          ZVM__(KFDIA,0:KLEV+1),         ZTHM__(KFDIA,0:KLEV+1)
REAL(KIND=JPRB) :: ZUS__(KFDIA,0:KLEV+1),          ZVS__(KFDIA,0:KLEV+1),         ZWS__(KFDIA,0:KLEV+1)
REAL(KIND=JPRB) :: ZTKES_OUT__(KFDIA,0:KLEV+1),    ZMF_UP__(KFDIA,0:KLEV+1),      ZTHVREFM__(KFDIA,0:KLEV+1)  ! thetav de l etat
REAL(KIND=JPRB) :: ZTENDU_TURB__(KFDIA,0:KLEV+1),  ZTENDV_TURB__(KFDIA,0:KLEV+1), ZTENDTHL_TURB__(KFDIA,0:KLEV+1)
REAL(KIND=JPRB) :: ZTENDRT_TURB__(KFDIA,0:KLEV+1), ZTKEM__(KFDIA,0:KLEV+1),       ZSRCS__(KFDIA,0:KLEV+1) 
REAL(KIND=JPRB) :: ZHLC_HRC__(KFDIA,0:KLEV+1),     ZHLC_HCF__(KFDIA,0:KLEV+1), &
                 & ZHLI_HRI__(KFDIA,0:KLEV+1),     ZHLI_HCF__(KFDIA,0:KLEV+1)

REAL(KIND=JPRB) :: ZSIGS__(KFDIA,0:KLEV+1),        ZEDR__(KFDIA,0:KLEV+1)
! THE DDH budgets
REAL(KIND=JPRB) :: ZDP__(KFDIA,0:KLEV+1),          ZTP__(KFDIA,0:KLEV+1),         ZTPMF__(KFDIA,0:KLEV+1)
REAL(KIND=JPRB) :: ZTDIFF__(KFDIA,0:KLEV+1),       ZTDISS__(KFDIA,0:KLEV+1)
! length scales for momentum and heat for mnh level definitions in case LHARATU=TRUE
REAL(KIND=JPRB) :: ZLENGTHM__(KFDIA,0:KLEV+1), ZLENGTHH__(KFDIA,0:KLEV+1)

REAL(KIND=JPRB), POINTER :: ZTHS__(:,:)
! horizontal gradients and diagnostics
REAL(KIND=JPRB) :: ZTURB3D__(KFDIA,0:KLEV+1,YDMODEL%YRML_PHY_MF%YRARPHY%NGRADIENTS)
! WARNING ! Don't use ZTHSWAP__ or ZTHSAVE__ below because they may be swapped !
! Use only the pointer ZTHS__, and possibly ZTHSIN_ if you need the backup of input data.
REAL(KIND=JPRB), TARGET  :: ZTHSWAP__(KFDIA,0:KLEV+1),        ZTHSAVE__(KFDIA,0:KLEV+1)
REAL(KIND=JPRB), TARGET  :: ZFLXZTHVMF_SUM__(KFDIA,0:KLEV+1), ZWM__(KFDIA,0:KLEV+1)


! Updraft characteristics for Meso-NH world (input of ARO_SHALLOW_MF)
REAL(KIND=JPRB) :: ZTHETAL_UP_(KFDIA,KLEV), ZTHETAV_UP_(KFDIA,KLEV), ZZFRAC_UP_(KFDIA,KLEV)
REAL(KIND=JPRB) :: ZRT_UP_(KFDIA,KLEV),     ZRC_UP_(KFDIA,KLEV),     ZRI_UP_(KFDIA,KLEV)
REAL(KIND=JPRB) :: ZZU_UP_(KFDIA,KLEV),     ZZV_UP_(KFDIA,KLEV),     ZZW_UP_(KFDIA,KLEV)
REAL(KIND=JPRB) :: ZZRV_UP_(KFDIA,KLEV),    ZTKES_(KFDIA,KLEV),      ZZZ_(KFDIA,KLEV)
REAL(KIND=JPRB) :: ZDZZ_(KFDIA,KLEV),       ZZZ_F_(KFDIA,KLEV),      ZDZZ_F_(KFDIA,KLEV)
REAL(KIND=JPRB) :: ZCIT_(KFDIA,KLEV),       ZMFM_(KFDIA,KLEV),       ZEXNREFM_(KFDIA,KLEV)
REAL(KIND=JPRB) :: ZSIGM_(KFDIA,KLEV),      ZNEBMNH_(KFDIA,KLEV),    ZEVAP_(KFDIA,KLEV)
! additions for MF scheme (Pergaud et al)
REAL(KIND=JPRB) :: ZSIGMF_(KFDIA,KLEV),     ZRC_MF_(KFDIA,KLEV),     ZRI_MF_(KFDIA,KLEV)
REAL(KIND=JPRB) :: ZCF_MF_(KFDIA,KLEV),     ZAERD_(KFDIA,KLEV),      ZCVTENDT_(KFDIA,KLEV)
REAL(KIND=JPRB) :: ZCVTENDRV_(KFDIA,KLEV),  ZCVTENDRC_(KFDIA,KLEV),  ZCVTENDRI_(KFDIA,KLEV)
REAL(KIND=JPRB) :: ZMFS_(KFDIA,KLEV),       ZTHLS_(KFDIA,KLEV),      ZRTS_(KFDIA,KLEV)
REAL(KIND=JPRB) :: ZMFUS_(KFDIA,KLEV),      ZMFVS_(KFDIA,KLEV),      ZDEPTH_HEIGHT_(KFDIA,KLEV)

REAL(KIND=JPRB), TARGET :: ZFLXZTHVMF_(KFDIA,KLEV)
REAL(KIND=JPRB), POINTER :: ZARG_FLXZTHVMF_(:,:)


! WARNING ! Don't use ZRSWAP_ or ZRSAVE_ below because they may be swapped !
! Use only the pointer ZRS_, and possibly ZRSIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZRSIN_(:,:,:), ZRS_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZRSWAP_(KFDIA,KLEV,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
REAL(KIND=JPRB), TARGET :: ZRSAVE_(KFDIA,KLEV,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

REAL(KIND=JPRB), POINTER  :: ZPTRWNU_(:,:), ZTHSIN_(:,:)
REAL(KIND=JPRB), TARGET :: ZWNU_(KFDIA,KLEV)

! WARNING ! Don't use ZSVSWAP_ or ZSVSAVE_ below because they may be swapped !
! Use only the pointer ZSVS_, and possibly ZSVSIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZSVSIN_(:,:,:), ZSVS_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZSVSWAP_(KFDIA,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB), TARGET :: ZSVSAVE_(KFDIA,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

REAL(KIND=JPRB) :: ZSVXXX_(KFDIA,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

! WARNING ! Don't use ZSVMSWAP_ or ZSVMSAVE_ below because they may be swapped !
! Use only the pointer ZSVM_, and possibly ZSVMIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZSVMIN_(:,:,:), ZSVM_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZSVMSWAP_(KFDIA,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) 
REAL(KIND=JPRB), TARGET :: ZSVMSAVE_(KFDIA,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB) :: ZSVMB_(KFDIA,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

! WARNING ! Don't use ZLIMASWAP_ or ZLIMASAVE_ below because they may be swapped !
! Use only the pointer ZLIMAS_, and possibly ZLIMASIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZLIMAS_(:,:,:), ZLIMASIN_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZLIMASWAP_(KFDIA,KLEV,YDMODEL%YRML_GCONF%YGFL%NLIMA)
REAL(KIND=JPRB), TARGET :: ZLIMASAVE_(KFDIA,KLEV,YDMODEL%YRML_GCONF%YGFL%NLIMA)

REAL(KIND=JPRB) :: ZLIMAM_(KFDIA,KLEV,YDMODEL%YRML_GCONF%YGFL%NLIMA)
!INTEGER(KIND=JPIM) :: KSV_TURB !CPtoclean?
!CPtoclean REAL(KIND=JPRB) :: ZTURBM(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA)
!CPtoclean REAL(KIND=JPRB) :: ZTURBS(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA)
!not (yet ?) used. REK
!REAL(KIND=JPRB) :: ZSFTURB(KLON,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA) ! surf. flux of SV (=0)
!REAL(KIND=JPRB) :: ZTENDSV_TURB2(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA)  ! SV (=0)
REAL(KIND=JPRB) :: ZSFSVLIMA_(KFDIA,YDMODEL%YRML_GCONF%YGFL%NLIMA) ! surf. flux of LIMA vars
REAL(KIND=JPRB) :: ZTENDSV_TURBLIMA_(KFDIA,KLEV,YDMODEL%YRML_GCONF%YGFL%NLIMA) ! LIMA

! For radiation scheme
REAL(KIND=JPRB) :: ZRM_(KFDIA,KLEV,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
REAL(KIND=JPRB) :: ZPFPR_(KFDIA,KLEV,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

REAL(KIND=JPRB) :: ZPEZDIAG_(KFDIA,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG)

REAL(KIND=JPRB) :: ZSFSV_(KFDIA,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! surf. flux of scalars

! Single scattering albedo of dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZPIZA_DST_(KFDIA,KLEV,YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)
! Assymetry factor for dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZCGA_DST_(KFDIA,KLEV,YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)
! tau/tau_{550} dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZTAUREL_DST_(KFDIA,KLEV,YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)


! surface flux of theta and surface flux of vapor ; surface flux of CO2
REAL(KIND=JPRB) :: ZSFTH_(KFDIA),  ZSFRV_(KFDIA),          ZSFCO2_(KFDIA)
REAL(KIND=JPRB) :: ZACPRG_(KFDIA), ZINPRG_NOTINCR_(KFDIA), ZINPRG_(KFDIA)
REAL(KIND=JPRB) :: ZACPRR_(KFDIA), ZINPRR_NOTINCR_(KFDIA), ZINPRR_(KFDIA)
REAL(KIND=JPRB) :: ZACPRS_(KFDIA), ZINPRS_NOTINCR_(KFDIA), ZINPRS_(KFDIA)
REAL(KIND=JPRB) :: ZCFBTH_(KFDIA), ZINPRH_NOTINCR_(KFDIA), ZINPRH_(KFDIA)
REAL(KIND=JPRB) :: ZZS_(KFDIA),    ZSSO_STDEV_(KFDIA),     ZALB_UV_(KIDIA)
REAL(KIND=JPRB) :: ZLAT_(KIDIA),   ZLON_(KIDIA),           ZZENITH_(KIDIA)
REAL(KIND=JPRB) :: ZGZ0_(KFDIA),   ZGZ0H_(KFDIA),          ZTOWNS_(KFDIA) 
REAL(KIND=JPRB) :: ZCFAQ_(KFDIA),  ZCFBQ_(KFDIA),          ZCFATH_(KFDIA)
REAL(KIND=JPRB) :: ZCFAU_(KFDIA),  ZCFBU_(KFDIA),          ZCFBV_(KFDIA)
REAL(KIND=JPRB) :: ZBUDTH_ (KFDIA),ZBUDSO_(KFDIA),         ZFCLL_(KFDIA)
REAL(KIND=JPRB) :: ZCD_(KFDIA),    ZSEA_(KFDIA),           ZTOWN_(KFDIA)
REAL(KIND=JPRB) :: ZZTOP_(KFDIA),  ZCVTENDPR_(KFDIA),      ZCVTENDPRS_(KFDIA)
! surface flux of x and y component of wind. are they really necessary ? REK
REAL(KIND=JPRB) :: ZSFU_(KFDIA),   ZSFV_(KFDIA)


! Arpege-style dimensionning :
! --------------------------

!Variables used in case LHARATU=TRUE
! length scales for momentum and heat and TKE
REAL(KIND=JPRB) :: ZLENGTH_M(KLON,KLEV),ZLENGTH_H(KLON,KLEV)
REAL(KIND=JPRB) :: ZTKEEDMF(KLON,KLEV)
REAL(KIND=JPRB) :: ZTKEEDMFS(KLON,KLEV)

REAL(KIND=JPRB) :: ZEMIS (KLON)
REAL(KIND=JPRB) :: ZQICE(KLON,KLEV), ZQLIQ(KLON,KLEV)

REAL(KIND=JPRB) :: ZAER(KLON,KLEV,6)
REAL(KIND=JPRB) :: ZAERINDS(KLON,KLEV)

REAL(KIND=JPRB) :: ZQSAT(KLON,KLEV)

REAL(KIND=JPRB) :: ZFRSOFS(KLON)
REAL(KIND=JPRB) :: ZLH(KLON,KLEV), ZLSCPE(KLON,KLEV), ZGEOSLC(KLON,KLEV)
REAL(KIND=JPRB) :: ZQDM(KLON,KLEV), ZQV(KLON,KLEV)
REAL(KIND=JPRB) :: ZQCO2(KLON,KLEV)
REAL(KIND=JPRB), TARGET :: ZTKEM4SLDDH(KLON,KLEV)
REAL(KIND=JPRB), POINTER :: ZTKEM(:,:)
REAL(KIND=JPRB) :: ZQW(KLON,KLEV), ZTW(KLON,KLEV)

REAL(KIND=JPRB) :: ZTENT(KLON,KLEV)
REAL(KIND=JPRB) :: ZTENDT(KLON,KLEV) ! array to save heating profile for LHN
REAL(KIND=JPRB) :: ZMAXTEND,ZMINTEND
REAL(KIND=JPRB) :: ZDZZ(KLON,KLEV)
REAL(KIND=JPRB) :: ZTPW(KLON,KLEV)

! POUR GROUND 
REAL(KIND=JPRB) :: ZZS_FSWDIR(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZZS_FSWDIF(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZTRSODIR (KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZTRSODIF (KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)

REAL(KIND=JPRB) :: ZTRSOD(KLON)
REAL(KIND=JPRB) :: ZSUDU(KLON), ZSDUR(KLON), ZDSRP(KLON)
REAL(KIND=JPRB) :: ZCEMTR(KLON,2), ZCTRSO(KLON,2)

REAL(KIND=JPRB) :: ZALBD(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZALBP(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZALBD1(KLON), ZALBP1(KLON)
REAL(KIND=JPRB) :: ZAPHIM(KLON,0:KLEV), ZAPHIFM(KLON,KLEV)

REAL(KIND=JPRB) :: ZTM(KLON,KLEV), ZQVM(KLON,KLEV), ZQIM(KLON,KLEV)
REAL(KIND=JPRB) :: ZQCM(KLON,KLEV),ZQHM(KLON,KLEV), ZQHGM(KLON,KLEV)
REAL(KIND=JPRB) :: ZQRM(KLON,KLEV), ZQSM(KLON,KLEV), ZQGM(KLON,KLEV)
REAL(KIND=JPRB) :: ZUPGENL(KLON,KLEV)
REAL(KIND=JPRB) :: ZUPGENN(KLON,KLEV)
REAL(KIND=JPRB) :: ZCLFR(KLON)

REAL(KIND=JPRB) :: ZCPM(KLON,KLEV), ZRHM(KLON,KLEV)

! Variables concerning updraft rain/snow for EDMF
REAL(KIND=JPRB) :: ZTENDTUP(KLON,KLEV), ZTENDQVUP(KLON,KLEV)

! specific to new data flow for diagnostics
REAL(KIND=JPRB) :: ZTENDTBAK(KLON,KLEV), ZTENDRBAK(KLON,KLEV,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
REAL(KIND=JPRB) :: ZTMPAF(KLON,KLEV)

! daand: radflex
REAL(KIND=JPRB)  :: ZFPR(KLON,0:KLEV,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

! Target should not be necessary. REK
REAL(KIND=JPRB), TARGET :: ZCON1(KLON,KLEV)
REAL(KIND=JPRB), TARGET :: ZCON2(KLON,KLEV)
REAL(KIND=JPRB), TARGET :: ZCON3(KLON,KLEV)

! Ajout pour MF Dual Scheme (KNMI et al)
! Updraft characteristics in Arpege/IFS world
REAL(KIND=JPRB) :: ZMF_UP(KLON,0:KLEV,KMAXDRAFT)
REAL(KIND=JPRB) :: ZTHETAL_UP(KLON,0:KLEV,KMAXDRAFT)
REAL(KIND=JPRB) :: ZQT_UP(KLON,0:KLEV,KMAXDRAFT)
REAL(KIND=JPRB) :: ZTHTV_UP(KLON,0:KLEV,KMAXDRAFT)
REAL(KIND=JPRB) :: ZQC_UP(KLON,0:KLEV,KMAXDRAFT)
REAL(KIND=JPRB) :: ZQI_UP(KLON,0:KLEV,KMAXDRAFT)
REAL(KIND=JPRB) :: ZU_UP(KLON,0:KLEV,KMAXDRAFT)
REAL(KIND=JPRB) :: ZV_UP(KLON,0:KLEV,KMAXDRAFT)
REAL(KIND=JPRB) :: ZTSURF(KLON), ZTN(KLON), ZQS(KLON)

REAL(KIND=JPRB) :: ZFRSOLU(KLON), ZFRSODS(KLON)
REAL(KIND=JPRB) :: ZFSDNN(KLON), ZFSDNV(KLON)

REAL(KIND=JPRB) :: ZSURFPREP(KLON), ZSURFSNOW(KLON)

REAL(KIND=JPRB) :: ZQO3(KLON,0:KLEV)

REAL(KIND=JPRB) :: ZZS_FTH_(KLON), ZZS_FRV_(KLON), ZZS_FU_(KLON), ZZS_FV_(KLON)

! Surface forcing arrays for MUSC
REAL(KIND=JPRB) :: ZRHODREFM(KLON), ZTHETAS(KLON)

! ACRANEB2 local variables
REAL(KIND=JPRB) :: ZNEB0    (KLON,KLEV)  ! protected cloud fractions
REAL(KIND=JPRB) :: ZCLCT_RAD(KLON)       ! total cloud cover for radiation
REAL(KIND=JPRB) :: ZDECRD   (KLON)       ! decorrelation depth

! Stochastic physics pattern & dummy tendencies for calling sppten
! Bof. REK
REAL(KIND=JPRB) :: ZDUMMY(KLON,KLEV)
REAL(KIND=JPRB) :: ZDUMMY1(KLON)

REAL(KIND=JPRB) :: ZROZ(KLON,KLEV)

! Can we remove ? REK
REAL(KIND=JPRB) :: ZEPSM(0,0,0) ! Dissipation of TKE (eps) at time t-dt
REAL(KIND=JPRB) :: ZEPSS(0,0,0) ! Dissipation of TKE at time t+dt


!    Integers
INTEGER(KIND=JPIM) :: JLEV, JLON, JRR, JGFL, JGR, ISPLITR
INTEGER(KIND=JPIM) :: IJN  ! max. number of day/night slices within NRPOMA
INTEGER(KIND=JPIM) :: IKL  !ordering of vert levels 1:MNH -1:AROME
INTEGER(KIND=JPIM) :: IOFF_MFSHAL, IEZDIAG_CHEM
INTEGER(KIND=JPIM) :: IKA,IKB,IKU,IKT,IKTE,IKTB ! vertical points as in mpa
INTEGER(KIND=JPIM) :: JSG, JK, JR, JSW
INTEGER(KIND=JPIM) :: IDRAFT,JDRAFT,INDRAFT
INTEGER(KIND=JPIM) :: ISURFEX
INTEGER(KIND=JPIM) :: IDAY,IYEAR,IMONTH

INTEGER(KIND=JPIM) :: INIT0 ! Kind of safety/debugging initialization :
                            ! 0 = initialize to HUGE (debugging)
                            ! 1 = initialize to realistic value (discouraged)
                            ! -1 = no initialization (optimized code) - this is the default.

INTEGER(KIND=JPIM) :: ICLPH(KLON)             !PBL top level
INTEGER(KIND=JPIM) :: JLHSTEP,ISTEP

!       Real
REAL(KIND=JPRB) :: ZRHO
REAL(KIND=JPRB) :: ZAEO, ZAEN, ZSALBCOR
REAL(KIND=JPRB) :: ZDT, ZDT2, ZINVDT, ZINVG, ZRSCP, ZINVATM, Z_WMAX, Z_WMIN
 ! pas de temps pour la surface externalise
REAL(KIND=JPRB) :: ZDTMSE,ZRHGMT,ZSTATI
REAL(KIND=JPRB) :: ZDELTA
REAL(KIND=JPRB) :: ZEPSNEB

! default values for initialization :
REAL(KIND=JPRB) :: ZVALUE, ZVALUE_ONE, ZVALUE_T, ZVALUE_P, ZVALUE_L, ZVALUE_EPSILON

REAL(KIND=JPRB) :: ZVETAH(0:KLEV)

!       Boolean
LOGICAL :: LLMSE, LLMSE_PARAM, LLMSE_DIAG
LOGICAL :: LLAROME
LOGICAL :: LLRAD
LOGICAL :: LLSWAP_THS, LLSWAP_RS, LLSWAP_SVS, LLSWAP_SVM, LLSWAP_LIMAS ! logical to swap or not pointers in and out
LOGICAL :: LLHN(KLON,KLEV)
LOGICAL :: LNUDGLHNREAD

!       Characters
CHARACTER(LEN=11) :: CLNAME
CHARACTER(LEN=2),DIMENSION(7):: CLVARNAME=(/"QV","QL","QR","QI","QS","QG","QH"/)

! daand: radflex
REAL(KIND=JPRB), POINTER :: ZFRSO(:,:), ZFRTH(:,:)
TYPE(TYPE_INTPROC), POINTER :: YLRADPROC
REAL(KIND=JPRB)   :: ZCAPE(KLON), ZDCAPE(KLON)

!
! Phaser team note from CY43T1:
! there was a USE MODD_CTURB for accessing XTKEMIN here, but that created a forbidden
! dependence of APL_AROME (in "ifsarp") to the Méso-NH/Arome interfaces (in "mpa").
! There should be no USE MODD_* in APL_*.
! We decided to change the variable here to a local one, with the classical initial value for TKE.
!
REAL(KIND=JPRB), PARAMETER :: PPTKEMIN = 1.E-6


! Perturbed radiation-cloud interaction coef
REAL(KIND=JPRB), DIMENSION (KLON) :: ZRADGR,ZRADSN
REAL(KIND=JPRB) :: ZMU,ZVAL
INTEGER(KIND=JPIM) :: JKO,JKE

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!     --------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "recmwf.intfb.h"
#include "acraneb2.intfb.h"
#include "actqsat.intfb.h"
#include "acnpart.intfb.h"
#include "bri2acconv.intfb.h"
#include "gpgeo.intfb.h"
#include "gprcp.intfb.h"
#include "radheat.intfb.h"
#include "suozon.intfb.h"
#include "radaer.intfb.h"
#include "radozc.intfb.h"
#include "accldia.intfb.h"
#include "vdfhghthl.intfb.h"
#include "sppten.intfb.h"
#include "surf_ideal_flux.intfb.h"
#include "ecr1d.intfb.h"
#include "apl_arome2intflex.intfb.h"

#include "aro_rain_ice.h"
#include "nudglhprecip.intfb.h"
#include "nudglh.intfb.h"
#include "nudglhclimprof.intfb.h"
#include "nudglhprep.intfb.h"
#include "aro_turb_mnh.h"
#include "aro_adjust.h"
#include "aro_mnhc.h"
#include "aro_mnhdust.h"
#include "aro_startbu.h"
#include "aro_convbu.h"
#include "aro_ground_param.h"
#include "aro_ground_diag.h"
#include "aro_shallow_mf.h"
#include "aro_rainaero.h"
#include "aro_lima.h"
#include "diagflash.intfb.h"
#include "dprecips.intfb.h"
#include "ppwetpoint.intfb.h"
#include "acvisih.intfb.h"
#include "aro_ground_diag_2isba.h"

!     ------------------------------------------------------------------

#include "fcttrm.func.h"

!     --------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('APL_AROME',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE,  YDSTA=>YDGEOMETRY%YRSTA, YDLAP=>YDGEOMETRY%YRLAP, &
 & YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  YDVSPLIP=>YDGEOMETRY%YRVSPLIP, YDVSLETA=>YDGEOMETRY%YRVSLETA, &
 & YDHSLMER=>YDGEOMETRY%YRHSLMER,  YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,  YDSPGEOM=>YDGEOMETRY%YSPGEOM, &
 & YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, &
 & YLDDH=>YDMODEL%YRML_DIAG%YRLDDH,YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH,  &
 & YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI,YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY,  &
 & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,YDPHY3=>YDMODEL%YRML_PHY_MF%YRPHY3,YDPHYDS=>YDMODEL%YRML_PHY_MF%YRPHYDS,  &
 & YDPHY0=>YDMODEL%YRML_PHY_MF%YRPHY0, YDVISI=>YDMODEL%YRML_PHY_MF%YRPHY%YRDVISI,&
 & YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,&
 & YGFL=>YDMODEL%YRML_GCONF%YGFL,YDMSE=>YDMODEL%YRML_PHY_MF%YRMSE, &
 & YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR,YDPRECIPS=>YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS, &
 & YDSTOPH=>YDMODEL%YRML_PHY_STOCH%YRSTOPH)

ASSOCIATE(MINPRR=>YDPARAR%MINPRR, MINPRS=>YDPARAR%MINPRS, MVQS=>YDPARAR%MVQS, &
 & MINPRG=>YDPARAR%MINPRG, LOTOWNC=>YDPARAR%LOTOWNC, LFPREC3D=>YDPARAR%LFPREC3D, &
 & NGPAR=>YDPARAR%NGPAR, CSUBG_PR_PDF=>YDPARAR%CSUBG_PR_PDF, &
 & NRRI=>YDPARAR%NRRI, NRRL=>YDPARAR%NRRL, CSUBG_AUCV_RC=>YDPARAR%CSUBG_AUCV_RC, &
 & CSUBG_AUCV_RI=>YDPARAR%CSUBG_AUCV_RI, CCONDENS=>YDPARAR%CCONDENS, &
 & CSUBG_MF_PDF=>YDPARAR%CSUBG_MF_PDF, &
 & LTOTPREC=>YDPARAR%LTOTPREC, CSUBG_RC_RR_ACCR=>YDPARAR%CSUBG_RC_RR_ACCR, CSUBG_RR_EVAP=>YDPARAR%CSUBG_RR_EVAP, &
 & NPRINTFR=>YDPARAR%NPRINTFR, CMF_CLOUD=>YDPARAR%CMF_CLOUD, &
 & MALBDIR=>YDPARAR%MALBDIR, NSWB_MNH=>YDPARAR%NSWB_MNH, &
 & XSW_BANDS=>YDPARAR%XSW_BANDS, MACPRG=>YDPARAR%MACPRG, MSWDIR=>YDPARAR%MSWDIR, &
 & LMIXUV=>YDPARAR%LMIXUV, MSWDIF=>YDPARAR%MSWDIF, LOLSMC=>YDPARAR%LOLSMC, &
 & NDIAGWMAX=>YDPARAR%NDIAGWMAX, MACPRS=>YDPARAR%MACPRS, MACPRR=>YDPARAR%MACPRR, &
 & LSQUALL=>YDPARAR%LSQUALL, VSIGQSAT=>YDPARAR%VSIGQSAT, &
 & MALBSCA=>YDPARAR%MALBSCA,&
 & RADSN=>YDPARAR%RADSN, LOSEDIC=>YDPARAR%LOSEDIC, LDIAGWMAX=>YDPARAR%LDIAGWMAX, &
 & CSEDIM=>YDPARAR%CSEDIM, CLAMBDA3=>YDPARAR%CLAMBDA3, &
 & NPTP=>YDPARAR%NPTP, NSPLITR=>YDPARAR%NSPLITR, NSPLITG=>YDPARAR%NSPLITG, NSV=>YDPARAR%NSV, &
 & CFRAC_ICE_SHALLOW_MF=>YDPARAR%CFRAC_ICE_SHALLOW_MF, CFRAC_ICE_ADJUST=>YDPARAR%CFRAC_ICE_ADJUST, &
 & MVTS=>YDPARAR%MVTS, NREFROI2=>YDPARAR%NREFROI2, NREFROI1=>YDPARAR%NREFROI1, &
 & MVEMIS=>YDPARAR%MVEMIS, LOWARM=>YDPARAR%LOWARM, LOCND2=>YDPARAR%LOCND2, &
 & LGRSN=>YDPARAR%LGRSN, LOSIGMAS=>YDPARAR%LOSIGMAS, NRR=>YDPARAR%NRR, &
 & LOSUBG_COND=>YDPARAR%LOSUBG_COND, RADGR=>YDPARAR%RADGR, &
 & CMF_UPDRAFT=>YDPARAR%CMF_UPDRAFT, LHARATU=>YDPARAR%LHARATU, &
 & XMINLM=>YDPHY0%XMINLM, XMAXLM=>YDPHY0%XMAXLM, AERCS1=>YDPHY0%AERCS1, &
 & AERCS3=>YDPHY0%AERCS3, AERCS5=>YDPHY0%AERCS5, &
 & RDECRD1=>YDPHY0%RDECRD1, RDECRD2=>YDPHY0%RDECRD2, &
 & RDECRD3=>YDPHY0%RDECRD3, RDECRD4=>YDPHY0%RDECRD4, &
 & LMPA=>YDARPHY%LMPA, LUSECHEM=>YDARPHY%LUSECHEM, LKFBCONV=>YDARPHY%LKFBCONV, &
 & LMFSHAL=>YDARPHY%LMFSHAL, LMICRO=>YDARPHY%LMICRO, &
 & CCOUPLING=>YDARPHY%CCOUPLING, LTURB=>YDARPHY%LTURB, LGRADHPHY=>YDARPHY%LGRADHPHY, &
 & NSURFEX_ITER=>YDARPHY%NSURFEX_ITER, LRDUST=>YDARPHY%LRDUST, &
 & NGRADIENTS=>YDARPHY%NGRADIENTS, &
 & LRDEPOS=>YDARPHY%LRDEPOS, LSURFEX_CRITICAL=>YDARPHY%LSURFEX_CRITICAL, &
 & LRCO2=>YDARPHY%LRCO2, LMSE=>YDARPHY%LMSE, &
 & LSURFEX_KFROM=>YDARPHY%LSURFEX_KFROM, &
 & NSURFEXCTL=>YDMSE%NSURFEXCTL, XZSEPS=>YDMSE%XZSEPS, &
 & NDLUNG=>YDDIM%NDLUNG, NDGUNG=>YDDIM%NDGUNG, NPROMA=>YDDIM%NPROMA, &
 & NDLUXG=>YDDIM%NDLUXG, NDGUXG=>YDDIM%NDGUXG, &
 & NSFORC=>YDPHYDS%NSFORC, &
 & NGFL_EXT=>YGFL%NGFL_EXT, YLRAD=>YGFL%YLRAD, YIRAD=>YGFL%YIRAD, &
 & NGFL_EZDIAG=>YGFL%NGFL_EZDIAG, &
 & NLIMA=>YGFL%NLIMA, CMICRO=>YDPARAR%CMICRO,NPROMICRO=>YDPARAR%NPROMICRO, &
 & YSD_VAD=>YDSURF%YSD_VAD, &
 & QCO2=>YDPHY3%QCO2, &
 & NTEND_DIAG_POS=>YDPHY%NTEND_DIAG_POS, NTEND_DIAG_FREQ_RESET=>YDPHY%NTEND_DIAG_FREQ_RESET, &
 & NRAY=>YDPHY%NRAY, LRAYFM=>YDPHY%LRAYFM, &
 & LO3ABC=>YDPHY%LO3ABC, LRAY=>YDPHY%LRAY, &
 & LAEROVOL=>YDPHY%LAEROVOL, LRSTAER=>YDPHY%LRSTAER, LRNUEXP=>YDPHY%LRNUEXP, &
 & AMAGSTOPH_CASBS=> YDSTOPH%AMAGSTOPH_CASBS, LFORCENL=>YDSTOPH%LFORCENL, &
 & NFORCESTART=>YDSTOPH%NFORCESTART, NFORCEEND=>YDSTOPH%NFORCEEND, &
 & NTRADI=>YDTOPH%NTRADI, NTQSAT=>YDTOPH%NTQSAT, NTNEBU=>YDTOPH%NTNEBU,&
 & NAER=>YDERAD%NAER, &
 & LHLRADUPD=>YDPHY%LHLRADUPD, &
 & TSPHY=>YDPHY2%TSPHY,&
 & NMODE=>YDERAD%NMODE, &
 & NOZOCL=>YDERAD%NOZOCL, &
 & NRADFR=>YDERAD%NRADFR, &
 & NSW=>YDERAD%NSW, &
 & RCARDI=>YDERDI%RCARDI, &
 & LFLEXDIA=>YLDDH%LFLEXDIA, LDDH_OMP=>YLDDH%LDDH_OMP, LRSLDDH=>YLDDH%LRSLDDH, &
 & RDECLI=>YDRIP%RDECLI, &
 & RCODEC=>YDRIP%RCODEC, &
 & RHGMT=>YDRIP%RHGMT, &
 & RSIDEC=>YDRIP%RSIDEC, &
 & RSOVR=>YDRIP%RSOVR, &
 & RSTATI=>YDRIP%RSTATI, &
 & TSTEP=>YDRIP%TSTEP, &
 & STPREH=>YDSTA%STPREH, &
 & LXXDIAGH=>YDXFU%LXXDIAGH,&
 & LFLASH =>YDCFU%LFLASH,&
 & LDPRECIPS=>YDPHY%LDPRECIPS,LDPRECIPS2=>YDPHY%LDPRECIPS2,&
 & NDTPREC=>YDPRECIPS%NDTPREC,NDTPREC2=>YDPRECIPS%NDTPREC2,&
 & NDTPRECCUR=>YDPRECIPS%NDTPRECCUR,NDTPRECCUR2=>YDPRECIPS%NDTPRECCUR2,&
 & NGPTOT=>YDGEM%NGPTOT,&
 & NGPBLKS=>YDDIM%NGPBLKS)
!     --------------------------------------------------------------------------

!    ------------------------------------------------------------------
!     1 - Initialisations
!    - --------------------------------------------------------------------

INIT0=-1

IF (INIT0 == 0) THEN
  ZVALUE=HUGE(1._JPRB)
  ZVALUE_ONE=HUGE(1._JPRB)
  ZVALUE_T=HUGE(1._JPRB)
  ZVALUE_P=HUGE(1._JPRB)
  ZVALUE_L=HUGE(1._JPRB)
  ZVALUE_EPSILON=HUGE(1._JPRB)
ELSE
  ZVALUE=0._JPRB
  ZVALUE_ONE=1._JPRB
  ZVALUE_T=293._JPRB
  ZVALUE_P=101325._JPRB
  ZVALUE_L=0.01_JPRB
  ZVALUE_EPSILON=1E-12_JPRB
ENDIF

LLSWAP_THS=.TRUE. ! it can be as well true as false, actually : this is just to start up the swapp process
LLSWAP_RS=.TRUE. ! it can be as well true as false, actually : this is just to start up the swapp process
LLSWAP_SVS=.TRUE. ! it can be as well true as false, actually : this is just to start up the swapp process
LLSWAP_SVM=.TRUE. ! it can be as well true as false, actually : this is just to start up the swapp process
LLSWAP_LIMAS=.TRUE. ! it can be as well true as false, actually : this is just to start up the swapp process

!        1.0 numerical safety

IF (JPRD == JPRB) THEN
  ZEPSNEB=1.E-12
ELSE
  ZEPSNEB=1.E-06
ENDIF

NSV=0

!         1.3 time step initialisation
!             the mesoNH physics (turb and microphysics) is written 
!             for leap frog scheme
!             !!! be carefull for 2TL or 3TL 

IF (LTWOTL) THEN
  ZDT=PDT/2._JPRB
ELSE
  IF (KSTEP/=0) THEN
    ZDT=PDT/2._JPRB
  ELSE
    ZDT=PDT
  ENDIF
ENDIF

ZINVDT=1/PDT

ZINVG=1._JPRB/RG 

! initialisation de ZDTMSE
IF (LDXFUMSE) THEN
  ZDTMSE=0.01_JPRB
  ZSTATI=REAL(RSTATI,JPRB)-ZDTMSE/2._JPRB
ELSE
  ZDTMSE=PDT
  ZSTATI=REAL(RSTATI,JPRB)
ENDIF

IF(LTWOTL) THEN
  ZRHGMT=REAL(RHGMT,JPRB)-ZDTMSE/2._JPRB
ELSE
  ZRHGMT=REAL(RHGMT,JPRB)
ENDIF


LLMSE=LMSE.AND.(NSURFEXCTL >= 2)
LLMSE_PARAM=LLMSE
LLMSE_DIAG=LLMSE.AND.(NSURFEXCTL >= 3)


!  Vertical points
IKA=KLEV
IKB=KLEV
IKU=1
IKT=KLEV
IKTE=KLEV
IKTB=1
IKL=-1

!  SETUP

IF (INIT0 >= 0) THEN

  ZUM__(:,:)=ZVALUE
  ZUS__(:,:)=ZVALUE
  ZVM__(:,:)=ZVALUE
  ZVS__(:,:)=ZVALUE
  ZWM__(:,:)=ZVALUE
  ZTHSWAP__(:,:)=ZVALUE
  ZTHSAVE__(:,:)=ZVALUE
  ZSRCS__(:,:)=ZVALUE
  ZSIGS__(:,:)=ZVALUE
  ZTHM__(:,:)=ZVALUE_T
  ZRHODREFM__(:,:)=ZVALUE_ONE
  ZPABSM__(:,:)=ZVALUE_P
  ZTURB3D__(:,:,:)=ZVALUE
  ZLENGTHM__(:,:)=ZVALUE_L
  ZLENGTHH__(:,:)=ZVALUE_L

  ZZZ_(:,:)=ZVALUE
  ZMFM_(:,:)=ZVALUE
  ZSIGM_(:,:)=ZVALUE
  ZNEBMNH_(:,:)=ZVALUE
  ZEVAP_(:,:)=ZVALUE

  ZRSWAP_(:,:,:)=ZVALUE
  ZRSAVE_(:,:,:)=ZVALUE

  ZRM_(:,:,:)=ZVALUE

  ZLIMASWAP_(:,:,:)=ZVALUE
  ZLIMASAVE_(:,:,:)=ZVALUE

  ZLIMAM_(:,:,:)=ZVALUE

  ZFLXZTHVMF_(:,:)=ZVALUE
  ZSIGMF_(:,:)=ZVALUE
  ZRC_MF_(:,:)=ZVALUE
  ZRI_MF_(:,:)=ZVALUE
  ZCF_MF_(:,:)=ZVALUE

  ZSVSWAP_(:,:,:)=ZVALUE
  ZSVSAVE_(:,:,:)=ZVALUE

  ZSVMSWAP_(:,:,:)=ZVALUE
  ZSVMSAVE_(:,:,:)=ZVALUE
  ZSVMB_(:,:)=ZVALUE

  ZPIZA_DST_(:,:,:)  = ZVALUE
  ZCGA_DST_(:,:,:)   = ZVALUE
  ZTAUREL_DST_(:,:,:) = ZVALUE_EPSILON

  ZAERD_(:,:)=ZVALUE

  ZMFS_(:,:)=ZVALUE

  ZFPR(:,:,:)=ZVALUE

  IF(LKFBCONV) THEN
    ZCVTENDRV_(:,:)=ZVALUE
    ZCVTENDRC_(:,:)=ZVALUE
    ZCVTENDRI_(:,:)=ZVALUE
    ZCVTENDT_(:,:)=ZVALUE
  ENDIF

  ZACPRG_(:)=ZVALUE
  ZINPRR_NOTINCR_(:)=ZVALUE
  ZINPRS_NOTINCR_(:)=ZVALUE
  ZINPRG_NOTINCR_(:)=ZVALUE
  ZINPRH_NOTINCR_(:)=ZVALUE

  ZTHLS_(:,:)=ZVALUE
  ZMFUS_(:,:)=ZVALUE
  ZMFVS_(:,:)=ZVALUE

  ZTKEEDMF(:,:)=ZVALUE

  ZSFTH_(:)=ZVALUE
  ZSFRV_(:)=ZVALUE
  ZSFU_(:)=ZVALUE
  ZSFV_(:)=ZVALUE

  ZSFCO2_(:)=ZVALUE
  ZEMIS(:)=ZVALUE_ONE

  ZQICE(:,:)=ZVALUE
  ZQLIQ(:,:)=ZVALUE
  ZQO3(:,:)=ZVALUE

  ZAER(:,:,:)=ZVALUE
  ZAERINDS(:,:)=ZVALUE

  ZQSAT(:,:)=ZVALUE
  ZLH(:,:)=ZVALUE
  ZLSCPE(:,:)=ZVALUE
  ZGEOSLC(:,:)=ZVALUE

  ZFRSOFS(:)=ZVALUE

  ZQW(:,:)=ZVALUE
  PRH(:,:)=ZVALUE
  ZTW(:,:)=ZVALUE

  ZTRSODIF(:,:)=ZVALUE
  ZTRSODIR(:,:)=ZVALUE
  ZZS_FSWDIR(:,:)=ZVALUE
  ZZS_FSWDIF(:,:)=ZVALUE

  ZSDUR(:)=ZVALUE
  ZDSRP(:)=ZVALUE
  ZALBD(:,:)=ZVALUE
  ZALBP(:,:)=ZVALUE
  ZALBD1(:)=ZVALUE
  ZALBP1(:)=ZVALUE

  ZQV(:,:)=ZVALUE
  ZSFSV_(:,:)=ZVALUE

  PFRSOC(:,:)=ZVALUE
  PFRTHC(:,:)=ZVALUE
  PFRSOPS(:)=ZVALUE
  PDIAGH(:)=ZVALUE

  ZTENDSV_TURBLIMA_(:,:,:)=ZVALUE

ENDIF

!  INITIALIZE (CUMULATED) TENDENCIES

DO JLEV=1,KLEV
  PTENDT(KIDIA:KFDIA,JLEV)=0.0_JPRB
  PTENDU(KIDIA:KFDIA,JLEV)=0.0_JPRB
  PTENDV(KIDIA:KFDIA,JLEV)=0.0_JPRB
  PTENDW(KIDIA:KFDIA,JLEV)=0.0_JPRB
  PTENDTKE(KIDIA:KFDIA,JLEV)=0.0_JPRB
ENDDO
DO JRR=1,NRR
  DO JLEV=1,KLEV
    PTENDR(KIDIA:KFDIA,JLEV,JRR)=0.0_JPRB
  ENDDO
ENDDO
DO JGFL=1,NGFL_EXT
  DO JLEV=1,KLEV
    PTENDEXT(KIDIA:KFDIA,JLEV,JGFL)=0.0_JPRB
  ENDDO
ENDDO
DO JGFL=1,NLIMA
  DO JLEV=1,KLEV
    PTENDLIMA(KIDIA:KFDIA,JLEV,JGFL)=0.0_JPRB
  ENDDO
ENDDO

!  INITIALIZE CUMULATED STUFF

! Small array, OK. REK
ZINPRH_(KIDIA:KFDIA)=0._JPRB
ZINPRR_(KIDIA:KFDIA)=0._JPRB
ZACPRR_(KIDIA:KFDIA)=0._JPRB
ZINPRS_(KIDIA:KFDIA)=0._JPRB
ZACPRS_(KIDIA:KFDIA)=0._JPRB
ZINPRG_(KIDIA:KFDIA)=0._JPRB


DO JLEV = 1,KLEV
  DO JLON = KIDIA,KFDIA
    ZZZ_F_(JLON,JLEV)=PAPHIFM(JLON,JLEV)*ZINVG
    ZTENDT(JLON,JLEV)=0._JPRB
  ENDDO
ENDDO

! adhoc solution to avoid negative tke values
! when SL advective ddh is activated 
IF (LRSLDDH) THEN
  DO JLEV=1, KLEV
    DO JLON = KIDIA, KFDIA
      ZTKEM4SLDDH(JLON,JLEV)=MAX(PTKEM(JLON,JLEV),PPTKEMIN)
    ENDDO
  ENDDO
  ZTKEM => ZTKEM4SLDDH(:,:)
ELSE
  ZTKEM => PTKEM(:,:)
  !test TKE > 0.
  IF (MINVAL(ZTKEM(KIDIA:KFDIA,1:KLEV)) <= 0._JPRB) THEN
    CALL ABOR1('TKE < 0 under APL_AROME check YTKE_NL%NREQIN')
  ENDIF
ENDIF


!test invalid combinations
IF (LHARATU .AND. CMF_UPDRAFT == 'EDKF') THEN
  CALL ABOR1('Combination LHARATU and EDKF not valid!')
ENDIF

!initialisation of first useful field for EZDIAG use in Chemistry/Dust
IOFF_MFSHAL=1
IF(LFPREC3D) IOFF_MFSHAL=2

!    ------------------------------------------------------------------
!     2 - CHANGEMENTS DE VARIABLES ET INVERSION DES NIVEAUX 
!     --------------------------------------------------------------------

IF (LMICRO.OR.LTURB.OR.LLMSE.OR.LKFBCONV) THEN

  ZRSCP=RD/RCPD
  ZINVATM=1/RATM
  !initialisation de ZZZ_
  DO JLEV = 1,KLEV
   !initialisation de qdm (utile localement pour calculer rho  
   !et convertir q en r 
    IF (NRR==7) THEN
      DO JLON = KIDIA,KFDIA
        ZQDM(JLON,JLEV)=1._JPRB-PQVM(JLON,JLEV)-PQCM(JLON,JLEV)-PQRM(JLON,JLEV)&
         & -PQIM(JLON,JLEV)-PQSM(JLON,JLEV)-PQGM(JLON,JLEV)-PQHM(JLON,JLEV)   
      ENDDO
    ELSE
      DO JLON = KIDIA,KFDIA
        ZQDM(JLON,JLEV)=1._JPRB-PQVM(JLON,JLEV)-PQCM(JLON,JLEV)-PQRM(JLON,JLEV)&
         & -PQIM(JLON,JLEV)-PQSM(JLON,JLEV)-PQGM(JLON,JLEV)
      ENDDO
    ENDIF
    DO JLON = KIDIA,KFDIA
   !initialisation de ZRHODREFM__ (=qd*zrho) 
      ZRHO=PAPRSFM(JLON,JLEV)/(PRM(JLON,JLEV)*PTM(JLON,JLEV))
      ZRHODREFM__(JLON,JLEV)=ZRHO*ZQDM(JLON,JLEV)
      ZRHODJM__(JLON,JLEV)=PDELPM(JLON,JLEV)*ZINVG
  !initialisation de ZEXNREFM_
      ZEXNREFM_(JLON,JLEV)=(PAPRSFM(JLON,JLEV)*ZINVATM)**(ZRSCP)
 ! vent horizontal et TKE
      ZPABSM__(JLON,JLEV)=PAPRSFM(JLON,JLEV)
      ZUM__(JLON,JLEV)= PUM(JLON,JLEV)
      ZVM__(JLON,JLEV)= PVM(JLON,JLEV)
      ZWM__(JLON,JLEV)= PWM(JLON,JLEV)
      ZTKEM__(JLON,JLEV)= ZTKEM(JLON,JLEV)
      ZZZ_(JLON,JLEV)=PAPHIM(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO
  !initialise sigma for subgrid condensation coming
  !from previous time step turbulence scheme
  IF (LOSIGMAS) THEN
    DO JLEV = 1, KLEV 
      ZSIGM_(KIDIA:KFDIA,JLEV)= PSIGM(KIDIA:KFDIA,JLEV)
    ENDDO
  ENDIF
  !initialise convective mas flux for subgrid condensation coming 
  !from previous time step convection scheme
  IF ( LKFBCONV.AND.LOSUBG_COND.AND..NOT.LOSIGMAS) THEN
    DO JLEV = 1, KLEV 
      ZMFM_(KIDIA:KFDIA,JLEV)=PSIGM(KIDIA:KFDIA,JLEV) 
    ENDDO
  ENDIF
!!! initialisation des variables d etat MNH �t

  !initialisation de ZRM_ pour les hydrometeores (ri=qi/qd)
  DO JLEV = 1, KLEV
    DO JLON = KIDIA,KFDIA
      ZTHM__(JLON,JLEV)=PTM(JLON,JLEV)/ZEXNREFM_(JLON,JLEV)
      ZRM_(JLON,JLEV,1)=PQVM(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,2)=PQCM(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,3)=PQRM(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,4)=PQIM(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,5)=PQSM(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,6)=PQGM(JLON,JLEV)/ZQDM(JLON,JLEV)  
    ENDDO
  ENDDO
  
  IF (NRR==7) THEN
    DO JLEV = 1, KLEV
      DO JLON = KIDIA,KFDIA
        ZRM_(JLON,JLEV,7)=PQHM(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ENDDO
    ENDDO
  ENDIF
   
  IF (NRR==6) THEN
    !initialisation de ZTHVREFM__
    DO JLEV = 1, KLEV 
      DO JLON = KIDIA,KFDIA
        ZTHVREFM__(JLON,JLEV)=ZTHM__(JLON,JLEV)*&
         & (1._JPRB+ZRM_(JLON,JLEV,1)*(RV/RD))/&
         & (1._JPRB+ZRM_(JLON,JLEV,1)+ZRM_(JLON,JLEV,2) +&
         & ZRM_(JLON,JLEV,3)+ZRM_(JLON,JLEV,4)+&
         & ZRM_(JLON,JLEV,5)+ZRM_(JLON,JLEV,6))
      ENDDO
    ENDDO
  ELSEIF (NRR==7) THEN
    DO JLEV = 1, KLEV 
      DO JLON = KIDIA,KFDIA
        ZTHVREFM__(JLON,JLEV)=ZTHM__(JLON,JLEV)*&
         & (1._JPRB+ZRM_(JLON,JLEV,1)*(RV/RD))/&
         & (1._JPRB+ZRM_(JLON,JLEV,1)+ZRM_(JLON,JLEV,2) +&
         & ZRM_(JLON,JLEV,3)+ZRM_(JLON,JLEV,4)+&
         & ZRM_(JLON,JLEV,5)+ZRM_(JLON,JLEV,6)+&
         &  ZRM_(JLON,JLEV,7) )
      ENDDO
    ENDDO
  ENDIF

!!! initialisation des variables d etat MNH a t+dt
!!! division pas le pas de temps
!!!(la multiplication par rhodj est faite plus tard, si necessaire, 
!!! suivant les parametrisations)   

  ! initialise pointers :
  CALL SWAP_THS
  ! vent horizontal
  DO JLEV = 1, KLEV 
    DO JLON=KIDIA,KFDIA
      ZUS__(JLON,JLEV)= PUM(JLON,JLEV)*ZINVDT
      ZVS__(JLON,JLEV)= PVM(JLON,JLEV)*ZINVDT
      ZWS__(JLON,JLEV)= PWM(JLON,JLEV)*ZINVDT
      ZTKES_(JLON,JLEV)= ZTKEM(JLON,JLEV)*ZINVDT
      ZTHS__(JLON,JLEV)=ZTHM__(JLON,JLEV)*ZINVDT
    ENDDO
  ENDDO

  !initialisation de ZRS_ pour les hydrometeores
  ! initialise pointers :
  CALL SWAP_RS
  DO JRR=1,NRR 
    DO JLEV = 1, KLEV
      DO JLON=KIDIA,KFDIA
        ZRS_(JLON,JLEV,JRR)=ZRM_(JLON,JLEV,JRR)*ZINVDT
      ENDDO
    ENDDO
  ENDDO

!!! Initialisations temporaires d'arguments non-utilises
  !initialisation de ZCIT_
  ZCIT_(KIDIA:KFDIA,1:IKT)=0.0_JPRB
  
  !initialisation des tableaux de precipitations inst. and cumulated 
  !and surface fluxes for turbulence
  IF (LLMSE.OR.LSFORCS) THEN
    ZACPRR_(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MACPRR)
    ZACPRS_(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MACPRS)
    ZACPRG_(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MACPRG)
    ZINPRR_NOTINCR_(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MINPRR)
    ZINPRS_NOTINCR_(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MINPRS)
    ZINPRG_NOTINCR_(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MINPRG)
  ENDIF

  !initialisation des scalaires passifs
  ! initialise pointers :
  CALL SWAP_SVM
  CALL SWAP_SVS
  DO JGFL=1,NGFL_EXT
    DO JLEV = 1, KLEV
      DO JLON=KIDIA,KFDIA
        ZSVM_(JLON,JLEV,JGFL)=PSVM(JLON,JLEV,JGFL)
        ZSVS_(JLON,JLEV,JGFL)=PSVM(JLON,JLEV,JGFL)*ZINVDT
      ENDDO
    ENDDO
  ENDDO

  !initialisation des concentrations LIMA
  ! initialise pointers :
  CALL SWAP_LIMAS
  DO JGFL=1,NLIMA
    DO JLEV = 1, KLEV
      DO JLON=KIDIA,KFDIA
        ZLIMAM_(JLON,JLEV,JGFL)=PLIMAM(JLON,JLEV,JGFL)
        ZLIMAS_(JLON,JLEV,JGFL)=PLIMAM(JLON,JLEV,JGFL)*ZINVDT
      ENDDO
    ENDDO
  ENDDO

ENDIF

! daand: radflex
ZFRSO => PFRSO(:,:,1)
ZFRTH => PFRTH(:,:,1)

!    ------------------------------------------------------------------
!     3 - PRINTS FOR DIAGNOSTICS IF NEEDED
!    ------------------------------------------------------------------
IF (LDIAGWMAX) THEN
  IF (MOD(KSTEP+1,NDIAGWMAX)==0) THEN
  ! calcul de wmax
    DO JLEV = 1 , KLEV
      Z_WMAX=0._JPRB
      Z_WMIN=0._JPRB
      DO JLON=KIDIA,KFDIA
        IF (PWM(JLON,JLEV)>Z_WMAX) THEN
          Z_WMAX=PWM(JLON,JLEV)
        ENDIF
        IF (PWM(JLON,JLEV)<Z_WMIN) THEN
          Z_WMIN=PWM(JLON,JLEV)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
ENDIF

IF (LFLEXDIA) THEN
  !save tendencies
  ZTENDTBAK(KIDIA:KFDIA,1:KLEV)=PTENDT(KIDIA:KFDIA,1:KLEV)
  DO JR=1,NRR
    ZTENDRBAK(KIDIA:KFDIA,1:KLEV,JR)=PTENDR(KIDIA:KFDIA,1:KLEV,JR)
  ENDDO
ENDIF


!    ------------------------------------------------------------------
!     4 - ADJUSTMENT (CALLED IF THE MICROPHYSICS IS SWITCH ON) 
!    ------------------------------------------------------------------

IF (LMICRO) THEN

  ! Swap pointers because input values of THS and RS should be saved
  CALL SWAP_THS
  CALL SWAP_RS

  IF (LMFSHAL .AND. (CMF_CLOUD=='DIRE'.OR.CMF_CLOUD=='BIGA')) THEN
    IOFF_MFSHAL=IOFF_MFSHAL+3
    IF (KSTEP==0) THEN
      DO JLEV = 1, KLEV 
        ZRC_MF_(KIDIA:KFDIA,JLEV)=0._JPRB
        ZRI_MF_(KIDIA:KFDIA,JLEV)=0._JPRB
        ZCF_MF_(KIDIA:KFDIA,JLEV)=0._JPRB
      ENDDO
    ELSE
      DO JLEV = 1, KLEV 
        ZRC_MF_(KIDIA:KFDIA,JLEV)=PEZDIAG(KIDIA:KFDIA,JLEV,1)
        ZRI_MF_(KIDIA:KFDIA,JLEV)=PEZDIAG(KIDIA:KFDIA,JLEV,3)
        ZCF_MF_(KIDIA:KFDIA,JLEV)=PEZDIAG(KIDIA:KFDIA,JLEV,2)
      ENDDO
    ENDIF
    PEZDIAG(KIDIA:KFDIA,1:KLEV,1:3)=0._JPRB
  ENDIF

  IF (MOD(KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'avant aro_adjust sous apl_arome'
    WRITE(NULOUT,*)'JLEV   ZZZ_F_   RHODJM   EXNREFM       PABSM       THM      SIGM         MFM    '
    DO JLEV=1,KLEV+1 
      WRITE(NULOUT,'(I2,X,7F10.3)') JLEV,ZZZ_F_(NPTP,JLEV),ZRHODJM__(NPTP,JLEV),&
       & ZEXNREFM_(NPTP,JLEV), ZPABSM__(NPTP,JLEV), ZTHM__(NPTP,JLEV), ZSIGM_(NPTP,JLEV), ZMFM_(NPTP,JLEV)
    ENDDO 
    WRITE (NULOUT,*)'JLEV  rhoQv  rhoQc   rhoQr   rhoQi   rhoQs   rhoQg'
    DO JLEV=1,KLEV
      WRITE(NULOUT,'(I2,X,6E11.4)') JLEV,ZRM_(NPTP,JLEV,1),&
       & ZRM_(NPTP,JLEV,2), ZRM_(NPTP,JLEV,3),ZRM_(NPTP,JLEV,4),ZRM_(NPTP,JLEV,5), ZRM_(NPTP,JLEV,6)  
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZRC_MF_  ZRI_MF_  ZCF_MF_ ZTHSIN_  ZSRCS__ ZNEBMNH_'
    DO JLEV=1,KLEV
      WRITE(NULOUT,'(I2,X,6E11.4)') JLEV,ZRC_MF_(NPTP,JLEV),&
       & ZRI_MF_(NPTP,JLEV),ZCF_MF_(NPTP,JLEV), ZTHSIN_(NPTP,JLEV),ZSRCS__(NPTP,JLEV), ZNEBMNH_(NPTP,JLEV)
    ENDDO
  ENDIF

  ! for now a copy is needed (see below, inside). I don't like than :-( REK
  ZTHS__(KIDIA:KFDIA,1:KLEV)=ZTHSIN_(KIDIA:KFDIA,1:KLEV)
  ZRS_(KIDIA:KFDIA,1:KLEV,1:NRR)=ZRSIN_(KIDIA:KFDIA,1:KLEV,1:NRR)

  IF (CMICRO == 'LIMA') THEN

    CALL SWAP_LIMAS
    ! for now a copy is needed (see below, inside). I don't like than :-( REK
    ZLIMAS_(KIDIA:KFDIA,1:KLEV,1:NLIMA)=ZLIMASIN_(KIDIA:KFDIA,1:KLEV,1:NLIMA)

    CALL ARO_ADJUST_LIMA (KLEV,IKU,IKL,KFDIA,KLEV,NRR,NLIMA,KSTEP+1,&
     & LOSUBG_COND, LOSIGMAS, LOCND2, &
     & ZDT,VSIGQSAT,ZZZ_F_,&
     & ZRHODJM__(:,1:KLEV),&
     & ZRHODREFM__(:,1:KLEV),&
     & ZEXNREFM_,&
     & ZPABSM__(:,1:KLEV),&
     & ZTHM__(:,1:KLEV),&
     & ZRM_,&
     & ZLIMAM_,&
     & ZSIGM_,&
     & ZMFM_,ZRC_MF_,&
     & ZRI_MF_,ZCF_MF_,&
     & ZTHS__(:,1:KLEV),ZRS_,&
     & ZLIMAS_,&
     & ZSRCS__(:,1:KLEV),ZNEBMNH_,&
     & YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH)
  ELSE
!    CALL ARO_ADJUST (KLON,KIDIA,KFDIA,KLEV,NRR,& !this is the target version
    CALL ARO_ADJUST (KFDIA,KIDIA,KFDIA,KLEV,NRR,&
     & NGFL_EZDIAG, &
     & CFRAC_ICE_ADJUST, CCONDENS, CLAMBDA3, LOSUBG_COND, &
     & LOSIGMAS, CMICRO, LOCND2, CSUBG_MF_PDF, &
     & ZDT,VSIGQSAT,ZZZ_F_,&
     & ZRHODJM__(:,1:KLEV),&
     & ZEXNREFM_,&
     & ZRHODREFM__(:,1:KLEV),&
     & ZPABSM__(:,1:KLEV),&
     & ZTHM__(:,1:KLEV),&
     & ZRM_,ZSIGM_,&
     & ZMFM_,ZRC_MF_,&
     & ZRI_MF_,ZCF_MF_,&
     & ZTHS__(:,1:KLEV),ZRS_,&
     & ZSRCS__(:,1:KLEV),ZNEBMNH_,&
     & ZHLC_HRC__(:,1:KLEV), ZHLC_HCF__(:,1:KLEV),&
     & ZHLI_HRI__(:,1:KLEV), ZHLI_HCF__(:,1:KLEV),&
     & PGP2DSPP,PEZDIAG,&
     & YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH)

  ENDIF
  
  IF (MOD(KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'apres aro_adjust sous apl_arome'
    WRITE(NULOUT,*)'JLEV   ZZZ_F_   RHODJM   EXNREFM       PABSM       THM      SIGM         MFM    '   
    DO JLEV=1,KLEV+1 
      WRITE(NULOUT,'(I2,X,7F10.3)') JLEV,ZZZ_F_(NPTP,JLEV),ZRHODJM__(NPTP,JLEV),&
       & ZEXNREFM_(NPTP,JLEV), ZPABSM__(NPTP,JLEV), ZTHM__(NPTP,JLEV), ZSIGM_(NPTP,JLEV), ZMFM_(NPTP,JLEV)
    ENDDO 
    WRITE (NULOUT,*)'JLEV  rhoQv  rhoQc   rhoQr   rhoQi   rhoQs   rhoQg'
    DO JLEV=1,KLEV
      WRITE(NULOUT,'(I2,X,6E11.4)') JLEV,ZRS_(NPTP,JLEV,1),&
       & ZRS_(NPTP,JLEV,2), ZRS_(NPTP,JLEV,3),ZRS_(NPTP,JLEV,4),ZRS_(NPTP,JLEV,5), ZRS_(NPTP,JLEV,6)  
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZRC_MF_  ZRI_MF_  ZCF_MF_ ZTHS__ ZSRCS__ ZNEBMNH_'
    DO JLEV=1,KLEV
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,ZRC_MF_(NPTP,JLEV),&
       & ZRI_MF_(NPTP,JLEV),ZCF_MF_(NPTP,JLEV), ZTHS__(NPTP,JLEV),ZSRCS__(NPTP,JLEV), ZNEBMNH_(NPTP,JLEV)
    ENDDO
  ENDIF

  DO JLEV=1,KLEV
    PCLFS(KIDIA:KFDIA,JLEV)=ZNEBMNH_(KIDIA:KFDIA,JLEV) 
  ENDDO

  !adjusted zthm and zrm
  DO JLEV = 1, KLEV
    DO JLON = KIDIA,KFDIA
      ZTHM__(JLON,JLEV)=ZTHS__(JLON,JLEV)*PDT
    ENDDO
  ENDDO

  DO JRR=1,NRR
    DO JLEV = 1, KLEV
      DO JLON = KIDIA,KFDIA
        ZRM_(JLON,JLEV,JRR)=ZRS_(JLON,JLEV,JRR)*PDT
      ENDDO
    ENDDO
  ENDDO


  !initialisation de qdm utile pour 
  !convertir tendance de r en tendance de q 
  IF (NRR==6) THEN
    DO JLEV=1,KLEV
      DO JLON= KIDIA, KFDIA
        ZQDM(JLON,JLEV)=1._JPRB/(1._JPRB+ZRM_(JLON,JLEV,1)+&
        &ZRM_(JLON,JLEV,2)+ZRM_(JLON,JLEV,3)+ZRM_(JLON,JLEV,4)+ZRM_(JLON,JLEV,5)+&
        &ZRM_(JLON,JLEV,6) )
      ENDDO
    ENDDO
  ELSEIF (NRR==7) THEN
    DO JLEV=1,KLEV
      DO JLON= KIDIA, KFDIA
        ZQDM(JLON,JLEV)=1._JPRB/(1._JPRB+ZRM_(JLON,JLEV,1)+&
        &ZRM_(JLON,JLEV,2)+ZRM_(JLON,JLEV,3)+ZRM_(JLON,JLEV,4)+ZRM_(JLON,JLEV,5)+&
        &ZRM_(JLON,JLEV,6)+ZRM_(JLON,JLEV,7) )
      ENDDO
    ENDDO
  ENDIF 
  !reinitialisation des qi
  DO JLEV = 1, KLEV
    DO JLON = KIDIA,KFDIA
      ZQVM(JLON,JLEV)=ZRM_(JLON,JLEV,1)*ZQDM(JLON,JLEV)
      ZQCM(JLON,JLEV)=ZRM_(JLON,JLEV,2)*ZQDM(JLON,JLEV)
      ZQRM(JLON,JLEV)=ZRM_(JLON,JLEV,3)*ZQDM(JLON,JLEV)
      ZQIM(JLON,JLEV)=ZRM_(JLON,JLEV,4)*ZQDM(JLON,JLEV)
      ZQSM(JLON,JLEV)=ZRM_(JLON,JLEV,5)*ZQDM(JLON,JLEV)
      ZQGM(JLON,JLEV)=ZRM_(JLON,JLEV,6)*ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  IF (NRR==7) THEN
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        ZQHM(JLON,JLEV)=ZRM_(JLON,JLEV,7)*ZQDM(JLON,JLEV)
      ENDDO
    ENDDO
  ELSE
    ZQHM(KIDIA:KFDIA,1:KLEV)=0._JPRB
  ENDIF

  ! Tendances des variables LIMA
  DO JGFL=1,NLIMA
    DO JLEV = 1, KLEV
      ! Réinitialisation des variables LIMA
      DO JLON=KIDIA,KFDIA
        ZLIMAM_(JLON,JLEV,JGFL)=ZLIMAS_(JLON,JLEV,JGFL)*PDT
        PTENDLIMA(JLON,JLEV,JGFL)=PTENDLIMA(JLON,JLEV,JGFL)+(ZLIMAS_(JLON,JLEV,JGFL)-ZLIMASIN_(JLON,JLEV,JGFL))  
      ENDDO
    ENDDO
  ENDDO

  !modif de R et CP
  ZQHGM(KIDIA:KFDIA,:)=ZQHM(KIDIA:KFDIA,:)+ZQGM(KIDIA:KFDIA,:)
  CALL GPRCP(KLON,KIDIA,KFDIA,KLEV,PQ=ZQVM,PQI=ZQIM,PQL=ZQCM,PQR=ZQRM,PQS=ZQSM,PQG=ZQHGM,PCP=ZCPM,PR=ZRHM)  

  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      ZTM(JLON,JLEV)=ZTHM__(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
  !reinitialisation de ZRHODREFM__ (=qd*zrho)
      ZRHO=PAPRSFM(JLON,JLEV)/(ZRHM(JLON,JLEV)*ZTM(JLON,JLEV))
      ZRHODREFM__(JLON,JLEV)=ZRHO*ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  !geopotentiel calculation
 
  ZAPHIM(KIDIA:KFDIA,KLEV)=PAPHIM(KIDIA:KFDIA,KLEV)
  CALL GPGEO(KLON,KIDIA,KFDIA,KLEV,ZAPHIM,ZAPHIFM,ZTM,ZRHM,PLNPRM,PALPHM,YDGEOMETRY%YRVERT_GEOM)
   
  !calcul de l'altitude
  DO JLEV = 1, KLEV
    DO JLON = KIDIA,KFDIA
      ZZZ_(JLON,JLEV)=ZAPHIM(JLON,JLEV)*ZINVG
      !initialisation de ZZZ_F_
      ZZZ_F_(JLON,JLEV)=ZAPHIFM(JLON,JLEV)*ZINVG
      ! tendency of T
      PTENDT(JLON,JLEV)=PTENDT(JLON,JLEV)+(ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV))*ZEXNREFM_(JLON,JLEV)
      ZTENDT(JLON,JLEV)=ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV)
    ENDDO
  ENDDO
  !inversion niveaux tendances des ri et conversion en qi en multipliant par qd 
  DO JR=1,NRR
    DO JLEV=1,KLEV
      DO JLON=KIDIA,KFDIA
        PTENDR(JLON,JLEV,JR)=PTENDR(JLON,JLEV,JR)+(ZRS_(JLON,JLEV,JR)-ZRSIN_(JLON,JLEV,JR))*ZQDM(JLON,JLEV)
      ENDDO
    ENDDO
  ENDDO
  !initialisation de ZDZZ_
  DO JLON = KIDIA,KFDIA
    ZDZZ_(JLON,1)=ZAPHIM(JLON,0)*ZINVG-ZZZ_(JLON,1)
  ENDDO
  DO JLEV = 2, KLEV
    DO JLON = KIDIA,KFDIA
      ZDZZ_(JLON,JLEV)=ZZZ_(JLON,JLEV+IKL)-ZZZ_(JLON,JLEV)
    ENDDO
  ENDDO

ELSE

  ZTM (KIDIA:KFDIA,1:KLEV)=PTM(KIDIA:KFDIA,1:KLEV)
  ZRHM(KIDIA:KFDIA,1:KLEV)=PRM(KIDIA:KFDIA,1:KLEV)
  ZQVM(KIDIA:KFDIA,1:KLEV)=PQVM(KIDIA:KFDIA,1:KLEV)
  ZQIM(KIDIA:KFDIA,1:KLEV)=PQIM(KIDIA:KFDIA,1:KLEV)
  ZQCM(KIDIA:KFDIA,1:KLEV)=PQCM(KIDIA:KFDIA,1:KLEV)
  ZQRM(KIDIA:KFDIA,1:KLEV)=PQRM(KIDIA:KFDIA,1:KLEV)
  ZQSM(KIDIA:KFDIA,1:KLEV)=PQSM(KIDIA:KFDIA,1:KLEV)
  ZQGM(KIDIA:KFDIA,1:KLEV)=PQGM(KIDIA:KFDIA,1:KLEV)
  IF (NRR==7) THEN
    ZQHM(KIDIA:KFDIA,1:KLEV)=PQHM(KIDIA:KFDIA,1:KLEV)
  ELSE
    ZQHM(KIDIA:KFDIA,1:KLEV)=0._JPRB
  ENDIF
  ZCPM(KIDIA:KFDIA,KTDIA:KLEV)=PCPM(KIDIA:KFDIA,KTDIA:KLEV)
  ZAPHIM(KIDIA:KFDIA,0:KLEV)=PAPHIM(KIDIA:KFDIA,0:KLEV)
  ZAPHIFM(KIDIA:KFDIA,1:KLEV)=PAPHIFM(KIDIA:KFDIA,1:KLEV)
  ZZZ_(KIDIA:KFDIA,1:KLEV)=PAPHIM(KIDIA:KFDIA,1:KLEV)*ZINVG
  !initialisation of PCLFS outside LMICRO to be zero in case LMICRO=F
  PCLFS(KIDIA:KFDIA,1:KLEV)=0._JPRB

ENDIF ! ADJUSTMENT LMICRO

IF (LFLEXDIA) THEN
  DO JLEV = 1, KLEV
    DO JLON=KIDIA,KFDIA 
      ZTMPAF(JLON,JLEV)=(PTENDT(JLON,JLEV)-ZTENDTBAK(JLON,JLEV))*PDELPM(JLON,JLEV)*ZINVG*ZCPM(JLON,JLEV)
    ENDDO
  ENDDO
  IF (LDDH_OMP) THEN
    CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,'TCTADJU',YDDDH)
  ELSE
    CALL ADD_FIELD_3D(YLDDH,ZTMPAF,'TCTADJU','T','ARP',.TRUE.,.TRUE.)
  ENDIF
  DO JR=1,NRR
    CLNAME='T'//CLVARNAME(JR)//'ADJU'
    DO JLEV = 1, KLEV
      DO JLON=KIDIA,KFDIA 
        ZTMPAF(JLON,JLEV)=(PTENDR(JLON,JLEV,JR)-ZTENDRBAK(JLON,JLEV,JR))*PDELPM(JLON,JLEV)*ZINVG
      ENDDO
    ENDDO
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,CLNAME,YDDDH)
    ELSE
      CALL ADD_FIELD_3D(YLDDH,ZTMPAF,CLNAME,'T','ARP',.TRUE.,.TRUE.)
    ENDIF
    DO JLEV = 1, KLEV
      DO JLON=KIDIA,KFDIA 
        ZTMPAF(JLON,JLEV)=PCLFS(JLON,JLEV)*PDELPM(JLON,JLEV)
      ENDDO
    ENDDO
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,'VNT',YDDDH)
    ELSE
      CALL ADD_FIELD_3D(YLDDH,ZTMPAF,'VNT','V','ARP',.TRUE.,.TRUE.)
    ENDIF
  ENDDO
! specific to new data flow for diagnostics
  IF (LDDH_OMP) THEN
    DO JLEV = 1, KLEV
      ZCON1(KIDIA:KFDIA,JLEV) = 1.0_JPRB
      ZCON2(KIDIA:KFDIA,JLEV) = ZQDM(KIDIA:KFDIA,JLEV)
    ENDDO
    DO JLEV = 1, KLEV
      DO JLON=KIDIA,KFDIA 
        ZCON3(JLON,JLEV) = PCPM(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
      ENDDO
    ENDDO
    ! missing interface !!! REK
    CALL ARO_SUINTBUDGET_OMP(KLON,KLEV,KSTEP,ZCON1,ZCON2,ZCON3,YDDDH)
  ELSE
    ! missing interface !!! REK
    CALL ARO_SUINTBUDGET(KLON,KLEV,KSTEP,KFDIA,ZQDM,ZEXNREFM_,PCPM)
  ENDIF
ENDIF


DO JLEV = 1, KLEV-1
  DO JLON = KIDIA,KFDIA
    ZDZZ_F_(JLON,JLEV)=ZZZ_F_(JLON,JLEV)-ZZZ_F_(JLON,JLEV-IKL)
  ENDDO
ENDDO
DO JLON = KIDIA,KFDIA
  ZDZZ_F_(JLON,KLEV)=ZZZ_F_(JLON,KLEV)-POROG(JLON)*ZINVG
ENDDO


!     --------------------------------------------------------------------
!     5 - COMPUTE DUST PROPERTIES FOR RADIATION IF LRDUST=T
!     --------------------------------------------------------------------
IF (LRDUST) THEN
  PEZDIAG(KIDIA:KFDIA,1:KLEV,IOFF_MFSHAL:NGFL_EZDIAG)=0.0_JPRB
  ! input dust scalar concentration in ppp from
  IEZDIAG_CHEM=NGFL_EZDIAG-IOFF_MFSHAL+1
  ! Swapp because IN and OUT will be needed simultaneously
  CALL SWAP_SVM
  ! input dust scalar concentration in ppp from
  CALL ARO_MNHDUST (IKL,KFDIA,KLEV,NGFL_EXT, PDT,ZSVMIN_,ZZZ_,ZDZZ_,&
   & ZPABSM__(:,1:KLEV),ZTHM__(:,1:KLEV),ZRHODREFM__(:,1:KLEV),&
   & NSWB_MNH,KSTEP+1,ZSVM_,ZPIZA_DST_,ZCGA_DST_,ZTAUREL_DST_,ZAERD_,IEZDIAG_CHEM,&
   & ZPEZDIAG_(:,:,IOFF_MFSHAL:NGFL_EZDIAG)           )
  PEZDIAG(KIDIA:KFDIA,1:KLEV,IOFF_MFSHAL:NGFL_EZDIAG)=ZPEZDIAG_(KIDIA:KFDIA,1:KLEV,IOFF_MFSHAL:NGFL_EZDIAG)
! return to tendency
  DO JGFL=1, NGFL_EXT
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        PTENDEXT(JLON,JLEV,JGFL)=PTENDEXT(JLON,JLEV,JGFL)+(ZSVM_(JLON,JLEV,JGFL)-ZSVMIN_(JLON,JLEV,JGFL))*ZINVDT
      ENDDO
    ENDDO
  ENDDO
ENDIF ! LRDUST

IF (LSFORCS) THEN   ! <== Surface forcing for MUSC

  ZTSURF(KIDIA:KFDIA) = PTM(KIDIA:KFDIA,KLEV)
  ZTN(KIDIA:KFDIA)    = PTM(KIDIA:KFDIA,KLEV)
  ZQS(KIDIA:KFDIA)    = PQVM(KIDIA:KFDIA,KLEV)
  DO JLON=KIDIA,KFDIA
    ZRHODREFM(JLON) = PAPRSFM(JLON,KLEV)/(PTM(JLON,KLEV)*PRM(JLON,KLEV))
    ZTHETAS(JLON)   = ZTSURF(JLON)*(RATM/PAPRSM(JLON,KLEV))**RKAPPA
  ENDDO

  LLAROME=.TRUE.
  CALL SURF_IDEAL_FLUX(YDRIP,YDPHY0,YDPHYDS, LLAROME, KIDIA , KFDIA  , KLON, PAPHIFM(:,KLEV), &
   & ZRHODREFM, PSFORC,ZTN,ZTSURF,PLSM,PQVM(:,KLEV), PUM(:,KLEV), PVM(:,KLEV), ZTHETAS, &
   & ZSFTH_, ZSFRV_, ZSFU_, ZSFV_)

!* Compute PBL-diagnostics
   
   ZCAPE(:)=0._JPRB
   ZDCAPE(:)=0._JPRB
  CALL ACCLDIA(YDXFU,YDPHY,YDMODEL%YRML_PHY_MF%YRPHY2,YDTOPH,KIDIA,KFDIA,KLON,KLEV,PUCLS,&
   & PVCLS,PUM(:,1:KLEV),PVM(:,1:KLEV),ZCAPE,ZDCAPE,ZTKEM(:,1:KLEV),PAPHIFM(:,1:KLEV),POROG,&
   & PUGST,PVGST,PPBLH,ICLPH)

  PPBLH(KIDIA:KFDIA)=MIN(XMAXLM,MAX(XMINLM,PPBLH(KIDIA:KFDIA))) 

ENDIF    ! <== End of surface forcing for MUSC

!     --------------------------------------------------------------------
!     6 - RADIATION LRAYFM (IFS) or LRAY (ACRANEB2) 
!     --------------------------------------------------------------------
IF (LRAYFM.OR.LRAY) THEN
  ! prepare some input for both radiation schemes at every time step

  ! test de coherence sur le nombre de bandes spectrales entre ce qui sort de
  ! la surface et ce qu'attend le rayonnement
  IF( NSWB_MNH /= NSW) THEN
    CALL ABOR1 (' NSWB_MNH must be equal to NSW !')
  ENDIF

  ! compute saturated specific humidity
  CALL ACTQSAT ( YDPHY,KIDIA,KFDIA,KLON,NTQSAT,KLEV, PAPRSFM, ZCPM, ZQVM, ZTM,&
   & ZGEOSLC, ZLH, ZLSCPE, ZQSAT, ZQW, PRH, ZTW)  

  IF (YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RADGR) THEN
   IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RADGR) THEN
    ZMU = -0.5_JPRB * (YSPP_CONFIG%CMPERT_RADGR * YSPP_CONFIG%SDEV)**2
   ELSE
    ZMU = 0._JPRB
   ENDIF
   DO JLON=KIDIA,KFDIA
    ZVAL = RADGR*EXP(ZMU+YSPP_CONFIG%CMPERT_RADGR*PGP2DSPP(JLON,YSPP%MP_RADGR))
    ZRADGR(JLON) = MAX(YSPP_CONFIG%CLIP_RADGR(1),MIN(ZVAL,YSPP_CONFIG%CLIP_RADGR(2)))
   ENDDO
   IF ( YSPP_CONFIG%IEZDIAG_POS > 0 ) THEN
    JKO=2*YSPP%MP_RADGR-1
    JKE=2*YSPP%MP_RADGR
    DO JLON=KIDIA,KFDIA
     PEZDIAG(JLON,JKO,YSPP_CONFIG%IEZDIAG_POS) = PGP2DSPP(JLON,YSPP%MP_RADGR)
     PEZDIAG(JLON,JKE,YSPP_CONFIG%IEZDIAG_POS) = ZRADGR(JLON)
    ENDDO
   ENDIF
  ELSE
   DO JLON=KIDIA,KFDIA
    ZRADGR(JLON) = RADGR
   ENDDO
  ENDIF

  IF (YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RADSN) THEN
   IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RADSN) THEN
    ZMU = -0.5_JPRB * (YSPP_CONFIG%CMPERT_RADSN * YSPP_CONFIG%SDEV)**2
   ELSE
    ZMU = 0._JPRB
   ENDIF
   DO JLON=KIDIA,KFDIA
    ZVAL = RADSN*EXP(ZMU+YSPP_CONFIG%CMPERT_RADSN*PGP2DSPP(JLON,YSPP%MP_RADSN))
    ZRADSN(JLON) = MAX(YSPP_CONFIG%CLIP_RADSN(1),MIN(ZVAL,YSPP_CONFIG%CLIP_RADSN(2)))
   ENDDO
   IF ( YSPP_CONFIG%IEZDIAG_POS > 0 ) THEN
    JKO=2*YSPP%MP_RADSN-1
    JKE=2*YSPP%MP_RADSN
    DO JLON=KIDIA,KFDIA
     PEZDIAG(JLON,JKO,YSPP_CONFIG%IEZDIAG_POS) = PGP2DSPP(JLON,YSPP%MP_RADSN)
     PEZDIAG(JLON,JKE,YSPP_CONFIG%IEZDIAG_POS) = ZRADSN(JLON)
    ENDDO
   ENDIF
  ELSE
   DO JLON=KIDIA,KFDIA
    ZRADSN(JLON) = RADSN
   ENDDO
  ENDIF
  ! initialisation des humidite (dans le rayonnement, l'eau liquide nuageuse 
  ! et la glace sont donne par des hu par rapport au gaz.
  ! (qi/qa+qv pour ice par ex. C'est donc different de ri)
  DO JLEV=KTDIA,KLEV
    DO JLON=KIDIA,KFDIA 
       ZQICE(JLON,JLEV)= MAX(0.0_JPRB,&
        & (ZQIM(JLON,JLEV) + ZQSM(JLON,JLEV)*ZRADSN(JLON) + ZQGM(JLON,JLEV)*ZRADGR(JLON))/&
        & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
        & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
       ZQLIQ(JLON,JLEV)=MAX(0.0_JPRB, ZQCM(JLON,JLEV)/&
        & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
        & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
       ZQV(JLON,JLEV)=MAX(0.0_JPRB, ZQVM(JLON,JLEV)/&
        & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
        & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
    ENDDO
  ENDDO

  ! store cloud water content for RTTOV
  IF (YIRAD%LGP) PQIRAD(KIDIA:KFDIA,:) = ZQICE(KIDIA:KFDIA,:)
  IF (YLRAD%LGP) PQLRAD(KIDIA:KFDIA,:) = ZQLIQ(KIDIA:KFDIA,:)

  ! Hannu Savijarvi diffuse -> direct albedo correction from hlradia,
  ! Assuming that SURFEX does not make difference between  
  ! dir/dif albedo as surfex/SURFEX/albedo_from_nir_vis.F90 defines
  ! PSCA_ALB(:,:) = PDIR_ALB(:,:)
  
! Albedo dans les intervalles, direct (parallel) et diffus (diffuse).
  IF (NSW==6.OR.NSW==1) THEN
    IF (LLMSE) THEN
      DO JSW=1,NSW
        ZALBP(KIDIA:KFDIA,JSW)=PGPAR(KIDIA:KFDIA,MALBDIR-1+JSW)
        ZALBD(KIDIA:KFDIA,JSW)=PGPAR(KIDIA:KFDIA,MALBSCA-1+JSW)
        IF (LHLRADUPD) THEN
          DO JLON=KIDIA,KFDIA
            ZSALBCOR=0.2_JPRB/(1._JPRB+PMU0(JLON))-0.12_JPRB
            ZALBP(JLON,JSW)=ZALBD(JLON,JSW)+ZSALBCOR
          ENDDO
        ENDIF
      ENDDO
    ELSEIF (LSFORCS) THEN
      DO JSW=1,NSW
        ZALBP(KIDIA:KFDIA,JSW)=RALB_FORC
        ZALBD(KIDIA:KFDIA,JSW)=RALB_FORC
!  direct>diffuse correction might be applied to RALB_FORC,too:
!              ZALBP(JLON,JSW)=RALB_FORC+ZSALBCOR
      ENDDO
    ELSE
     !pour pouvoir tourner sans la surface
      DO JSW=1,NSW
        ZALBP(KIDIA:KFDIA,JSW)=PALBIN(KIDIA:KFDIA)
        ZALBD(KIDIA:KFDIA,JSW)=PALBIN(KIDIA:KFDIA)
!              ZALBP(JLON,JSW)=PALBIN(JLON)+ZSALBCOR
      ENDDO
    ENDIF

  ! Spectral average albedo done with RSUN2 weights, 
  ! to be applied for HLRADIA, ACRANEB2 which use a single solar spectral band
    IF (LHLRADUPD) THEN
      ZALBP1(KIDIA:KFDIA)=0._JPRB
      ZALBD1(KIDIA:KFDIA)=0._JPRB
      DO JSW=1,NSW
        DO JLON=KIDIA,KFDIA
          ZALBP1(JLON)=ZALBP1(JLON)+RSUN2(JSW)*ZALBP(JLON,JSW)
          ZALBD1(JLON)=ZALBD1(JLON)+RSUN2(JSW)*ZALBD(JLON,JSW)
        ENDDO
      ENDDO
    ELSE
       ZALBP1(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MALBDIR)
       ZALBD1(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MALBSCA)
    ENDIF

  ELSE

    CALL ABOR1 ('ALBEDO FOR NSW/= 1 or 6 not defined in apl_arome')

  ENDIF

  ! all albedo operations

  IF (LLMSE) THEN
    ZEMIS(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MVEMIS)
    ZTSURF(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MVTS)
    ! protection for E Zone, Where surface scheme send back EMIS and T =0
    ! the protection in aro_ground_paramn is not sufficient !!! WHY ??
    DO JLON=KIDIA,KFDIA
      IF (ZEMIS(JLON)==0._JPRB) THEN
        ZEMIS(JLON)=1.0_JPRB
        ZTSURF(JLON)=288.0_JPRB
      ENDIF
    ENDDO
  ELSEIF (LSFORCS) THEN
    ZEMIS(KIDIA:KFDIA)=REMIS_FORC
  ELSE
    ZEMIS(KIDIA:KFDIA)=0.5_JPRB ! value 0.5 is suspicious
    ZTSURF(KIDIA:KFDIA)=ZTM(KIDIA:KFDIA,KLEV)
  ENDIF !LLMSE EMIS
  
  ! initialization of aerosols
  IF (LRAY.AND.NRAY == 2.AND.LRADFLEX.AND.LRSTAER) THEN

    ! old ("standard") aerosols for LRAY only
    ZAER(KIDIA:KFDIA,1:KTDIA-1,1)=0._JPRB
    DO JLEV=KTDIA-1,KLEV
      ZVETAH(JLEV)=STPREH(JLEV)/VP00
    ENDDO
    ZAEO=AERCS1*ZVETAH(KTDIA-1)+AERCS3*ZVETAH(KTDIA-1)**3+AERCS5*ZVETAH(KTDIA-1)**5
    DO JLEV=KTDIA,KLEV
      ZAEN=AERCS1*ZVETAH(JLEV)+AERCS3*ZVETAH(JLEV)**3+AERCS5*ZVETAH(JLEV)**5
      ZAER(KIDIA:KFDIA,JLEV,1)=ZAEN-ZAEO
      ZAEO=ZAEN
    ENDDO
    ZAER(KIDIA:KFDIA,1:KLEV,2:6)=0._JPRB
  
  ELSE
    
    IF (NAER >= 1 ) THEN
      IF(YSD_VAD%NUMFLDS >= 4) THEN
        CALL RADAER ( YDMODEL%YRML_PHY_RAD%YREAERD,YDERAD,YDPHY, KIDIA , KFDIA , KLON  , KLEV,&
         & PAPRSM,PAPRSFM,ZTM,ZTSURF,PAESEA,PAELAN,PAESOO,PAEDES,PAESUL,PAEVOL,ZAER,ZAERINDS)
      ELSE
        WRITE(NULOUT,*) 'YSD_VAD%NUMFLDS SHOULD BE >= 4, IT IS: ',YSD_VAD%NUMFLDS
        CALL ABOR1('APL_AROME: PB AEROSOLS!')
        ! NB : this abort excludes the use of radact. REK.
      ENDIF
    ENDIF

    IF (LRDUST) THEN
      ! We use the extinction coefficient explicitly solved by ARO_MNHDUST
      ZAER(KIDIA:KFDIA,1:KLEV,3) = ZAERD_(KIDIA:KFDIA,1:KLEV)
    ENDIF

  ENDIF
  ! end of old or new aerosols

  ! initialization of ozone
  IF (NOZOCL == 1) THEN
    ! as in IFS
    CALL RADOZC(KIDIA,KFDIA,KLON,KLEV,1,KLON,0,PAPRSM,PGEMU,ZROZ)
    DO JK=1,KLEV
      DO JLON=KIDIA,KFDIA
        ZQO3(JLON,JK)=ZROZ(JLON,JK)/PDELPM(JLON,JK)
      ENDDO
    ENDDO
  ELSEIF (NOZOCL == 2) THEN
    ! as in ARPEGE (from clim profiles)
    CALL SUOZON(KIDIA,KFDIA,KLON,KLEV,ZQO3,.FALSE.,PAPRSM,PRDELPM,LO3ABC,PVO3ABC)
  ENDIF

ELSE

  DO JSW=1,NSW
    ZALBP(KIDIA:KFDIA,JSW)=0._JPRB 
    ZALBD(KIDIA:KFDIA,JSW)=0._JPRB
  ENDDO

ENDIF
 !of preparation of input for LRAYFM, LRAY at every time step
 
 IF (LRAYFM) THEN
   ! Intermittent call to radiation interface
   IF (MOD(KSTEP,NRADFR) == 0) THEN 
     CALL RECMWF (YDGEOMETRY%YRDIMV, YDMODEL,         &
     & KIDIA   , KFDIA , KLON     , KLEV     ,           & 
     & ZALBD   , ZALBP , PAPRSM   , PAPRSFM  ,           &
     & PCLFS   , ZQO3  , ZAER     , PDELPM   , ZEMIS   , &
     & PMU0M   , ZQV   , ZQSAT    , ZQICE    , ZQLIQ   , &
     & ZQSM    , ZQRM  , PLSM     , ZTM      , ZTSURF  , &
     & PGP2DSPP, PEZDIAG,                                &
     & PEMTD   , PEMTU , PTRSO    , PFRTHC   , PFRTH   , & 
     & PFRSOC  , PFRSO , ZZS_FSWDIR, ZZS_FSWDIF, ZFSDNN  , &
     & ZFSDNV  , ZCTRSO, ZCEMTR   , ZTRSOD   , ZTRSODIR, &
     & ZTRSODIF, ZPIZA_DST_,ZCGA_DST_,ZTAUREL_DST_,ZAERINDS,&
     & PGELAM , PGEMU ,PGPAR , &
     & PMU0LU , ZALBD1 , ZFRSOLU)
   ELSE
     IF (LLMSE) THEN 
       DO JSW=1,NSW 
         ZTRSODIR(KIDIA:KFDIA,JSW)=PGPAR(KIDIA:KFDIA,MSWDIR+JSW-1)
         ZTRSODIF(KIDIA:KFDIA,JSW)=PGPAR(KIDIA:KFDIA,MSWDIF+JSW-1)
       ENDDO 
     ENDIF
     ZCTRSO(:,:)=0._JPRB
   ENDIF
   ! daand: radflex
   IF (LRADFLEX) THEN
     YLRADPROC => NEWINTPROC(YDPROCSET,'Radiation')
     ZFRSO => NEWINTFIELD(YLRADPROC,KLON,KLEV,'FRSO','H','F')
     ZFRTH => NEWINTFIELD(YLRADPROC,KLON,KLEV,'FRTH','H','F')
   ENDIF

    DO JLEV=1,KLEV
      ZTENT(KIDIA:KFDIA,JLEV)=0.0_JPRB
    ENDDO

   ZSUDU(KIDIA:KFDIA)=0.0_JPRB

   CALL RADHEAT&
   & (  YDERAD,YDERDI,YDMODEL%YRML_PHY_MF, KIDIA  , KFDIA  , KLON   , KLEV,&
   & PAPRSM  , ZEMIS  , PEMTD  , PMU0, ZQVM,&
   & ZTENT  , PTRSO  , ZTRSOD , ZTSURF   , PDT,&
   & ZTRSODIR,ZTRSODIF, ZALBD , ZALBP,&
   ! daand: radflex; replaced PFRSO and PRFTH by pointers
   & ZFRSO  , ZFRTH  , PFRSODS, PFRTHDS, ZCEMTR , ZCTRSO , PFRSOC , PFRTHC,&
   & ZSUDU  , ZSDUR  , ZDSRP  , ZZS_FSWDIR , ZZS_FSWDIF    ,&
   & PFRSOPS, ZFRSOFS, PFRSOPT )

  ! daand: radflex
  IF (LRADFLEX) THEN
    ! store for further calculations and diagnostics
    ! warning : pointers. REK
    PFRSO(:,:,1)=ZFRSO
    PFRTH(:,:,1)=ZFRTH
  ELSE
  ! daand: if LRADFLEX, the contribution to temperature is done by
  ! cptend_flex/cputqy
    ! update temperature tendency by radiative contribution
    DO JLEV=1,KLEV
      DO JLON = KIDIA, KFDIA
        PTENDT(JLON,JLEV)=PTENDT(JLON,JLEV)+ZTENT(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF

  DO JLON = KIDIA, KFDIA
    ! update sunshine duration [s]
    PSDUR(JLON)=PSDUR(JLON)+ZSDUR(JLON)*TSTEP
    ! Estimate of the direct normal irradiance, with securities
    IF (PMU0(JLON) > 3.0E-02_JPRB) THEN
      PFRSDNI(JLON)=MAX(0.0_JPRB,PFRSOPS(JLON)/PMU0(JLON))
    ELSE
      PFRSDNI(JLON)=MAX(0.0_JPRB,PFRSOPS(JLON))
    ENDIF
  ENDDO 

  IF( MOD(KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'sous apl_arome apres rayonnement ZTENT=',ZTENT(NPTP,30:41)
    IF (LLMSE) THEN
      DO JSW=1, NSW
        WRITE(NULOUT,*)'ZSFSWDIR ZSFSWDIF ZFSDNN ZFSDNV PFRSO',&
         & ZZS_FSWDIR(NPTP,JSW),ZZS_FSWDIF(NPTP,JSW),ZFSDNN(NPTP), ZFSDNV(NPTP),PFRSO(NPTP,KLEV,1)
        WRITE(NULOUT,*)'ZALBD ZALBP',ZALBD(NPTP,JSW),ZALBP(NPTP,JSW)
      ENDDO
    ENDIF
    WRITE(NULOUT,*)ZFSDNN(NPTP),ZFSDNV(NPTP)
    WRITE (NULOUT,*)'TSURF EMIS ZFRTH',ZTSURF(NPTP),ZEMIS(NPTP),PFRTHDS(NPTP)
  ENDIF

  IF (LFLEXDIA) THEN
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,PFRSO(:,:,1),'FCTRAYSO',YDDDH)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,PFRTH(:,:,1),'FCTRAYTH',YDDDH)
    ELSE
      CALL ADD_FIELD_3D(YLDDH,PFRSO(:,:,1),'FCTRAYSO','F','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_3D(YLDDH,PFRTH(:,:,1),'FCTRAYTH','F','ARP',.TRUE.,.TRUE.)
    ENDIF
  ENDIF

ELSE

  PFRSOC(KIDIA:KFDIA,0:1)=0.0_JPRB
  PFRTHC(KIDIA:KFDIA,0:1)=0.0_JPRB

ENDIF  ! LRAYFM

!     ------------------------------------------------------------------
!     NEBULOSITE (CONVECTIVE+STRATIFORME) A TROIS NIVEAUX.
!     DIAGNOSTIC OF THREE LEVELS (CONVECTIVE+STRATIFORM) CLOUDINESS.

! protect cloudiness from being 0 or 1 (needed for ACRANEB2 and ACNPART)
DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZNEB0(JLON,JLEV)=MAX(ZEPSNEB,MIN(1._JPRB-ZEPSNEB,PCLFS(JLON,JLEV)))
  ENDDO
ENDDO

! decorrelation depth for cloud overlaps

IF (LRNUEXP) THEN
  DO JLON=KIDIA,KFDIA
    ZDECRD(JLON)=RDECRD1+RDECRD2*EXP(-((ASIN(PGEMU(JLON))-RDECRD3*RDECLI)/RDECRD4)**2)
  ENDDO
ENDIF

! calculate high, medium, low and total cloud cover
CALL ACNPART(YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTNEBU,KLEV,&
 & PAPHIM,PAPHIFM,PAPRSFM,ZDECRD,ZNEB0,&
 & PCLCH,PCLCM,PCLCL,PCLCT,ZCLCT_RAD)

IF (LRAY.AND.NRAY == 2.AND.LRADFLEX) THEN

  ! -------------------------
  ! ACRANEB2 radiation scheme
  ! -------------------------

!+++ The next input preparations are redundant:

  ! initialization of cloud ice, cloud liquid and specific humidity
  ! (with respect to moist air, i.e. excluding hydrometeors)
  DO JLEV=KTDIA,KLEV
    DO JLON=KIDIA,KFDIA
      ZQICE(JLON,JLEV)=MAX(0.0_JPRB, ZQIM(JLON,JLEV)/&
       & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
       & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
      ZQLIQ(JLON,JLEV)=MAX(0.0_JPRB, ZQCM(JLON,JLEV)/&
       & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
       & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
      ZQV(JLON,JLEV)=MAX(0.0_JPRB, ZQVM(JLON,JLEV)/&
       & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
       & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
    ENDDO
  ENDDO

  ! store cloud water content for RTTOV
  IF (YIRAD%LGP) PQIRAD(KIDIA:KFDIA,:) = ZQICE(KIDIA:KFDIA,:)
  IF (YLRAD%LGP) PQLRAD(KIDIA:KFDIA,:) = ZQLIQ(KIDIA:KFDIA,:)

  ! initialization of ozone
  IF (NOZOCL == 1) THEN
    ! as in IFS
    CALL RADOZC(KIDIA,KFDIA,KLON,KLEV,1,KLON,0,PAPRSM,PGEMU,ZROZ)
    DO JK=1,KLEV
      DO JLON=KIDIA,KFDIA
        ZQO3(JLON,JK)=ZROZ(JLON,JK)/PDELPM(JLON,JK)
      ENDDO
    ENDDO
  ELSEIF (NOZOCL == 2) THEN
    ! as in ARPEGE (from clim profiles)
    CALL SUOZON(KIDIA,KFDIA,KLON,KLEV,ZQO3,.FALSE.,PAPRSM,PRDELPM,LO3ABC,PVO3ABC)
  ENDIF

  ! initialization of aerosols
  IF (LRAY.AND.NRAY == 2.AND.LRADFLEX.AND.LRSTAER) THEN

    ZAER(KIDIA:KFDIA,1:KTDIA-1,1)=0._JPRB
    ! old ("standard") aerosols
    DO JLEV=KTDIA-1,KLEV
      ZVETAH(JLEV)=STPREH(JLEV)/VP00
    ENDDO
    ZAEO=AERCS1*ZVETAH(KTDIA-1)+AERCS3*ZVETAH(KTDIA-1)**3+AERCS5*ZVETAH(KTDIA-1)**5
    DO JLEV=KTDIA,KLEV
      ZAEN=AERCS1*ZVETAH(JLEV)+AERCS3*ZVETAH(JLEV)**3+AERCS5*ZVETAH(JLEV)**5
      ZAER(KIDIA:KFDIA,JLEV,1)=ZAEN-ZAEO
      ZAEO=ZAEN
    ENDDO
    ZAER(KIDIA:KFDIA,1:KLEV,2:6)=0._JPRB
  
  ELSE

    IF (NAER >= 1) THEN
      IF (YSD_VAD%NUMFLDS >= 4) THEN
        ! initialisation of aerosols as in ARPEGE (from clim files)
        CALL RADAER (YDMODEL%YRML_PHY_RAD%YREAERD,YDERAD,YDPHY, KIDIA , KFDIA , KLON  , KLEV,&
         & PAPRSM,PAPRSFM,ZTM,ZTSURF,PAESEA,PAELAN,PAESOO,PAEDES,PAESUL,PAEVOL,ZAER,ZAERINDS)
      ELSE
        CALL ABOR1('APL_AROME: PB AEROSOLS!')
        ! NB : this abort excludes the use of radact. REK.
      ENDIF
    ENDIF

    IF (LRDUST) THEN
      ! We uses the extinction coefficient explicitely solved by ARO_MNHDUST
      ZAER(KIDIA:KFDIA,1:KLEV,3) = ZAERD_(KIDIA:KFDIA,1:KLEV)
    ENDIF

  ENDIF ! (LRAY.AND.NRAY == 2.AND.LRADFLEX.AND.LRSTAER) THEN

  ! get diffuse and direct surface albedo, emissivity and temperature
  IF (.NOT.LHLRADUPD) THEN
    ZALBD1(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MALBSCA)
    ZALBP1(KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MALBDIR)
  ENDIF
  ZEMIS  (KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MVEMIS)
  ZTSURF (KIDIA:KFDIA)=PGPAR(KIDIA:KFDIA,MVTS)
  DO JLON=KIDIA,KFDIA
    ! protection of E-zone (not to have zero emissivity and T_surf there)
    IF (ZEMIS(JLON) == 0._JPRB) THEN
      ZEMIS (JLON)=  1._JPRB
      ZTSURF(JLON)=288._JPRB
    ENDIF
  ENDDO

!+++ End of redundant input preparations for ACRANEB

  ! initialization of CO2(+), differs from IFS radiation scheme!
  ZQCO2(KIDIA:KFDIA,1:KLEV)=QCO2

  ! daand: radflex
  YLRADPROC => NEWINTPROC(YDPROCSET,'Radiation')
  ZFRSO => NEWINTFIELD(YLRADPROC,KLON,KLEV, 'FRSO','H','F')
  ZFRTH => NEWINTFIELD(YLRADPROC,KLON,KLEV, 'FRTH','H','F')

  ! call radiation scheme
  IJN=KLON
  CALL ACRANEB2(YDERDI,YDRIP,YDMODEL%YRML_PHY_MF,&
   & KIDIA,KFDIA,KLON,NTRADI,KLEV,IJN,KSTEP,KNFRRC,&
   & PAPRSM,PAPRSFM,PCPM,PRM,PDELPM,ZNEB0,&
   & ZQV,ZQCO2,ZQICE,ZQLIQ,ZQO3,PTM,&
   & ZALBD1,ZALBP1,ZEMIS,PGELAM,PGEMU,PMU0,PMU0LU,ZTSURF, &
   & ZDECRD,ZCLCT_RAD,&
   & PGDEOSI,PGUEOSI,PGMU0,PGMU0_MIN,PGMU0_MAX,&
   & PGDEOTI,PGDEOTI2,PGUEOTI,PGUEOTI2,PGEOLT,PGEOXT,&
   & PGRPROX,PGMIXP,PGFLUXC,PGRSURF,PSDUR,&
   & ZFRSO,ZFRTH,&
   & PFRSOC,PFRTHC,ZFRSODS,PFRSOPS,ZFRSOLU,PFRTHDS,ZAER)

  ! daand: radflex
  ! store for further calculations and diagnostics
  ! warning : pointers. REK
  PFRSO(:,:,1)=ZFRSO
  PFRTH(:,:,1)=ZFRTH

  ! extract surface fluxes
  DO JLON=KIDIA,KFDIA
    PFRSODS(JLON)=ZFRSODS(JLON)+PFRSOPS(JLON)  ! downward surface sw flux
  ENDDO

  IF (LLMSE) THEN
    IF (LHLRADUPD) THEN
      DO JSW = 1,NSW
        DO JLON=KIDIA,KFDIA
          ZZS_FSWDIR(JLON,JSW) = PFRSOPS(JLON)*RSUN2(JSW)
          ZZS_FSWDIF(JLON,JSW) = ZFRSODS(JLON)*RSUN2(JSW)
         ENDDO
      ENDDO
    ELSE
      ZZS_FSWDIR(KIDIA:KFDIA,1)=PFRSOPS(KIDIA:KFDIA) ! direct surface swdn flux
      ZZS_FSWDIF(KIDIA:KFDIA,1)=ZFRSODS(KIDIA:KFDIA) ! diffuse surface swdn flux
    ENDIF
  ENDIF

  ! Estimate of the direct normal irradiance, with securities
  PFRSDNI(KIDIA:KFDIA)=PFRSOPS(KIDIA:KFDIA)
  DO JLON = KIDIA, KFDIA
    IF (PMU0(JLON) > 3.0E-02_JPRB) THEN
      PFRSDNI(JLON)=PFRSOPS(JLON)/PMU0(JLON)
    ENDIF
  ENDDO
  DO JLON = KIDIA, KFDIA
    PFRSDNI(JLON)=MAX(0.0_JPRB,PFRSDNI(JLON))
  ENDDO

  IF (LFLEXDIA) THEN
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,PFRSO(:,:,1),'FCTRAYSO',YDDDH)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,PFRTH(:,:,1),'FCTRAYSO',YDDDH)
    ELSE
      CALL ADD_FIELD_3D(YLDDH,PFRSO(:,:,1),'FCTRAYSO','F','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_3D(YLDDH,PFRTH(:,:,1),'FCTRAYTH','F','ARP',.TRUE.,.TRUE.)
    ENDIF
  ENDIF

ENDIF

IF (.NOT.(LRAY.AND.NRAY == 2.AND.LRADFLEX).AND..NOT.LRAYFM) THEN
  DO JSW = 1,NSW
    ZZS_FSWDIR(KIDIA:KFDIA,JSW) = 0._JPRB
    ZZS_FSWDIF(KIDIA:KFDIA,JSW) = 0._JPRB
  ENDDO
  PFRSOPS(KIDIA:KFDIA)=0._JPRB
ENDIF

IF (LFLEXDIA) THEN
  CALL ARO_STARTBU( KIDIA, KFDIA, KLEV, NRR,NGFL_EXT,ZRHODJM__(:,1:KLEV),&
   & ZUS__(:,1:KLEV), ZVS__(:,1:KLEV), ZWS__(:,1:KLEV), ZTHS__(:,1:KLEV), &
   & ZRS_, ZTKES_, YDDDH,YDMODEL%YRML_DIAG%YRLDDH,YDMODEL%YRML_DIAG%YRMDDH)
ENDIF


!    ------------------------------------------------------------------
!     7 - CONVECTION. 
!     --------------------------------------------------------------------

IF(LKFBCONV) THEN

  ! No swapp needed becaus IN and OUT are not needed simultaneously

  CALL BRI2ACCONV(YDMODEL%YRML_PHY_MF,YDGEOMETRY%YREGEO,KIDIA,KFDIA,KFDIA,KLEV,PGM(KIDIA:KFDIA),&
   & ZPABSM__(:,1:KLEV),ZZZ_F_, ZTM(KIDIA:KFDIA,:), ZRM_(:,:,1),ZRM_(:,:,2), ZRM_(:,:,4), &
   & ZRHODREFM__(:,1:KLEV), ZUM__(:,1:KLEV),ZVM__(:,1:KLEV), ZWM__(:,1:KLEV),ZMFS_,&
   & ZCVTENDT_, ZCVTENDRV_,ZCVTENDRC_, ZCVTENDRI_,ZCVTENDPR_, ZCVTENDPRS_ &
   & )

  IF(MOD(KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)"Pluie conv au sol", ZCVTENDPR_(NPTP), &
     & MAXVAL(ZCVTENDPR_(KIDIA:KFDIA)) ,MINVAL(ZCVTENDPR_(KIDIA:KFDIA))
  ENDIF

  DO JLEV = 1,KLEV
    DO JLON = KIDIA, KFDIA
      PTENDT(JLON,JLEV)=PTENDT(JLON,JLEV) + ZCVTENDT_(JLON,JLEV)
      PTENDR(JLON,JLEV,1) = PTENDR(JLON,JLEV,1) + ZCVTENDRV_(JLON,JLEV)*ZQDM(JLON,JLEV)  
      PTENDR(JLON,JLEV,2) = PTENDR(JLON,JLEV,2) + ZCVTENDRC_(JLON,JLEV)*ZQDM(JLON,JLEV)  
      PTENDR(JLON,JLEV,4) = PTENDR(JLON,JLEV,4) + ZCVTENDRI_(JLON,JLEV)*ZQDM(JLON,JLEV)  
      ZRS_(JLON,JLEV,1)=ZRS_(JLON,JLEV,1)+ZCVTENDRV_(JLON,JLEV)
      ZRS_(JLON,JLEV,2)=ZRS_(JLON,JLEV,2)+ZCVTENDRC_(JLON,JLEV)
      ZRS_(JLON,JLEV,4)=ZRS_(JLON,JLEV,4)+ZCVTENDRI_(JLON,JLEV)
      ZTHS__(JLON,JLEV)=ZTHS__(JLON,JLEV)+ZCVTENDT_(JLON,JLEV)*(RATM/PAPRSFM(JLON,JLEV))**(RD/RCPD)  
    ENDDO
  ENDDO
  DO JLON =KIDIA, KFDIA
    ZINPRR_(JLON)=ZINPRR_(JLON)+ZCVTENDPR_(JLON)-ZCVTENDPRS_(JLON)
    ZACPRR_(JLON)=ZACPRR_(JLON)+(ZCVTENDPR_(JLON)-ZCVTENDPRS_(JLON))*PDT
    ZINPRS_(JLON)=ZINPRS_(JLON)+ZCVTENDPRS_(JLON)
    ZACPRS_(JLON)=ZACPRS_(JLON)+ZCVTENDPRS_(JLON)*PDT
  ENDDO
  ! avance temporelle et inversion niveau pour ZMFS_
  ! on utilise PSIGS pour le flux de masse pour la condensation sous maille 
  ! car PSIGS n est utilise que si LOSIGMAS=T
  IF (LOSUBG_COND.AND..NOT.LOSIGMAS) THEN
    PSIGS(KIDIA:KFDIA,1:KLEV)=ZMFS_(KIDIA:KFDIA,1:KLEV)
  ENDIF
  IF(MOD(KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)"aps CONV, TENRV, TENRC, TENRI"
    DO JLEV=1,KLEV
      WRITE(NULOUT,*)PTENDR(NPTP,JLEV,1),PTENDR(NPTP,JLEV,2),PTENDR(NPTP,JLEV,4)
    ENDDO
  ENDIF
  CALL ARO_CONVBU(KFDIA,KLEV,NRR,ZRHODJM__(:,1:KLEV),ZRS_,ZTHS__(:,1:KLEV), &
   & YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH)

ENDIF

!    ------------------------------------------------------------------
!     8 - SURFACE. 
!     --------------------------------------------------------------------

IF (LLMSE) THEN
! A loop around SURFEX in order to test OpenMP

  SURFEX_LOOP : DO ISURFEX = 1, NSURFEX_ITER

! Initialisations 

  DO JLON=KIDIA,KFDIA
    ZZS_(JLON)=POROG(JLON)*ZINVG 
  ENDDO
  DO JLEV = 1,KLEV
    DO JLON=KIDIA,KFDIA
      ZDEPTH_HEIGHT_(JLON,JLEV)=ZZZ_F_(JLON,JLEV)-ZZS_(JLON)  
    ENDDO
  ENDDO
  IF (MINVAL(ZDEPTH_HEIGHT_(KIDIA:KFDIA,IKB)) <= 0._JPRB) THEN
    DO JLON=KIDIA,KFDIA
      IF (ZDEPTH_HEIGHT_(JLON,IKB) <= 0._JPRB) THEN
        WRITE (NULOUT,*)'sous apl_arome pb height en', JLON,ZAPHIFM(JLON,KLEV),POROG(JLON)
      ENDIF
    ENDDO
  ENDIF
  ! Can't use a section of pointer. An explicit copy shows, by the way, that a copy is needed
  ! because data is not contiguous. REK
  ZSVMB_(KIDIA:KFDIA,1:NGFL_EXT)=ZSVM_(KIDIA:KFDIA,IKB,1:NGFL_EXT)

  IF (LSURFEX_CRITICAL) THEN

!$OMP CRITICAL (ARO_GROUND_PARAM_LOCK)

    IF (LLMSE_PARAM) THEN
      CALL ARO_GROUND_PARAM( KBL,KGPCOMP,&
       & KFDIA,KIDIA,KFDIA,KSTEP,&
       & NRR,NSW,NGFL_EXT,NDGUNG, NDGUXG, NDLUNG, NDLUXG,LSURFEX_KFROM,&
       & LMPA,CCOUPLING,LDXFUMSE,&
       & NINDAT,ZRHGMT,ZSTATI,RSOVR,RCODEC,RSIDEC,&
       & PINDX(KIDIA:KFDIA),PINDY(KIDIA:KFDIA),&
       & ZUM__(:,IKB),&
       & ZVM__(:,IKB),&
       & ZTM(KIDIA:KFDIA,KLEV),ZRM_(:,IKB,1),&
       & ZSVMB_,&
       & RCARDI,ZRHODREFM__(:,IKB),&
       & ZPABSM__(:,IKB),PAPRSM(KIDIA:KFDIA,KLEV),&
       & ZDTMSE,ZDEPTH_HEIGHT_(:,IKB),ZZS_, XZSEPS,&
       & PMU0(KIDIA:KFDIA),PMU0N(KIDIA:KFDIA),PGELAM(KIDIA:KFDIA),&
       & PGEMU(KIDIA:KFDIA),XSW_BANDS,&
       & ZINPRR_NOTINCR_,ZINPRS_NOTINCR_,&
       & ZINPRG_NOTINCR_,&
       & PFRTHDS(KIDIA:KFDIA),ZZS_FSWDIF(KIDIA:KFDIA,1:NSW),&
       & ZZS_FSWDIR(KIDIA:KFDIA,1:NSW),&
       & ZCFAQ_, ZCFATH_, ZCFAU_,ZCFBQ_, ZCFBTH_, ZCFBU_,ZCFBV_,&
       & ZSFTH_,ZSFRV_,&
       & ZSFSV_,ZSFCO2_,&
       & ZSFU_,ZSFV_,&
       & ZALBP(KIDIA:KFDIA,1:NSW),ZALBD(KIDIA:KFDIA,1:NSW),&
       & ZEMIS(KIDIA:KFDIA),ZTSURF(KIDIA:KFDIA),PFRTH(KIDIA:KFDIA,KLEV,1))

    ENDIF

    IF (LRCO2) THEN
      ZSFSV_(KIDIA:KFDIA,NSV_CO2)= ZSFCO2_(KIDIA:KFDIA)
!print*,' FLUX CO2 =', MINVAL(ZSFSV_(KIDIA:KFDIA,NSV_CO2)),&
!                    & MAXVAL(ZSFSV_(KIDIA:KFDIA,NSV_CO2))
    ENDIF

!!!!! TEST DDH ATTENTION
!ZSFRV_(KIDIA:KFDIA) = 0._JPRB

    IF (LLMSE_DIAG) THEN

      CALL ARO_GROUND_DIAG( KBL, KGPCOMP,&
       & KFDIA,KIDIA,KFDIA,KLEV, IKL,&
       & NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM,&
       & ZZS_,ZSFRV_,&
       & ZUM__(:,IKTB:IKTE),&
       & ZVM__(:,IKTB:IKTE),&
       & ZDEPTH_HEIGHT_(:,IKTB:IKTE),&
       & PFRTH(KIDIA:KFDIA,KLEV,1),PFRSO(KIDIA:KFDIA,KLEV,1),&
       & PINDX(KIDIA:KFDIA),PINDY(KIDIA:KFDIA),&
       & ZQS(KIDIA:KFDIA),ZGZ0_,ZGZ0H_,&
       & PTCLS(KIDIA:KFDIA),PQCLS(KIDIA:KFDIA),PHUCLS(KIDIA:KFDIA),&
       & PUCLS(KIDIA:KFDIA),PVCLS(KIDIA:KFDIA),&
       & PNUCLS(KIDIA:KFDIA),PNVCLS(KIDIA:KFDIA),&
       & PFCLL(KIDIA:KFDIA,1),PFCLN(KIDIA:KFDIA,1),&
       & PFEVL(KIDIA:KFDIA,1),PFEVN(KIDIA:KFDIA,1),&
       & ZSSO_STDEV_, PSPSG(KIDIA:KFDIA),&
       & ZBUDTH_, ZBUDSO_,&
       & ZFCLL_, ZTOWNS_,&
       & ZCD_                         )
      CALL ARO_GROUND_DIAG_2ISBA( KBL, KGPCOMP, &
       & KFDIA, KIDIA, KFDIA, &
       & NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, &
       & PINDX(KIDIA:KFDIA), PINDY(KIDIA:KFDIA), &
       & PLSM, ZDUMMY1, ZDUMMY1, ZDUMMY1, ZTSURF(KIDIA:KFDIA), PSPSG(KIDIA:KFDIA), &
       & ZDUMMY1(KIDIA:KFDIA), ZDUMMY1(KIDIA:KFDIA), ZDUMMY1(KIDIA:KFDIA), &
       & ZDUMMY1(KIDIA:KFDIA), ZDUMMY1(KIDIA:KFDIA), ZDUMMY1(KIDIA:KFDIA), &
       & ZDUMMY1(KIDIA:KFDIA), PSPSGR(KIDIA:KFDIA), ZDUMMY1 )

    ENDIF
 
!$OMP END CRITICAL (ARO_GROUND_PARAM_LOCK)
  ELSE

    IF (LLMSE_PARAM) THEN

      CALL ARO_GROUND_PARAM( KBL,KGPCOMP,&
       & KFDIA,KIDIA,KFDIA,KSTEP,&
       & NRR,NSW,NGFL_EXT,NDGUNG, NDGUXG, NDLUNG, NDLUXG,LSURFEX_KFROM,&
       & LMPA,CCOUPLING,LDXFUMSE,&
       & NINDAT,ZRHGMT,ZSTATI,RSOVR,RCODEC,RSIDEC,&
       & PINDX(KIDIA:KFDIA),PINDY(KIDIA:KFDIA),&
       & ZUM__(:,IKB),&
       & ZVM__(:,IKB),&
       & ZTM(KIDIA:KFDIA,KLEV),ZRM_(:,IKB,1),&
       & ZSVMB_,&
       & RCARDI,ZRHODREFM__(:,IKB),&
       & ZPABSM__(:,IKB),PAPRSM(KIDIA:KFDIA,KLEV),&
       & ZDTMSE,ZDEPTH_HEIGHT_(:,IKB),ZZS_, XZSEPS,&
       & PMU0(KIDIA:KFDIA),PMU0N(KIDIA:KFDIA),PGELAM(KIDIA:KFDIA),&
       & PGEMU(KIDIA:KFDIA),XSW_BANDS,&
       & ZINPRR_NOTINCR_,ZINPRS_NOTINCR_,&
       & ZINPRG_NOTINCR_,&
       & PFRTHDS(KIDIA:KFDIA),ZZS_FSWDIF(KIDIA:KFDIA,1:NSW),&
       & ZZS_FSWDIR(KIDIA:KFDIA,1:NSW),&
       & ZCFAQ_, ZCFATH_, ZCFAU_,ZCFBQ_, ZCFBTH_, ZCFBU_,ZCFBV_,&
       & ZSFTH_,ZSFRV_,&
       & ZSFSV_,ZSFCO2_,&
       & ZSFU_,ZSFV_,&
       & ZALBP(KIDIA:KFDIA,1:NSW),ZALBD(KIDIA:KFDIA,1:NSW),&
       & ZEMIS(KIDIA:KFDIA),ZTSURF(KIDIA:KFDIA),PFRTH(KIDIA:KFDIA,KLEV,1))

    ENDIF

    IF (LRCO2) THEN
      ZSFSV_(KIDIA:KFDIA,NSV_CO2)= ZSFCO2_(KIDIA:KFDIA)
!print*,' FLUX CO2 =', MINVAL(ZSFSV_(KIDIA:KFDIA,NSV_CO2)),&
!                    & MAXVAL(ZSFSV_(KIDIA:KFDIA,NSV_CO2))
    ENDIF

!!!!! TEST DDH ATTENTION
!ZSFRV_(KIDIA:KFDIA) = 0._JPRB

    IF (LLMSE_DIAG) THEN

      CALL ARO_GROUND_DIAG( KBL, KGPCOMP,&
       & KFDIA,KIDIA,KFDIA,KLEV, IKL,&
       & NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM,&
       & ZZS_,ZSFRV_,&
       & ZUM__(:,IKTB:IKTE),&
       & ZVM__(:,IKTB:IKTE),&
       & ZDEPTH_HEIGHT_(:,IKTB:IKTE),&
       & PFRTH(KIDIA:KFDIA,KLEV,1),PFRSO(KIDIA:KFDIA,KLEV,1),&
       & PINDX(KIDIA:KFDIA),PINDY(KIDIA:KFDIA),&
       & ZQS(KIDIA:KFDIA),ZGZ0_,ZGZ0H_,&
       & PTCLS(KIDIA:KFDIA),PQCLS(KIDIA:KFDIA),PHUCLS(KIDIA:KFDIA),&
       & PUCLS(KIDIA:KFDIA),PVCLS(KIDIA:KFDIA),&
       & PNUCLS(KIDIA:KFDIA),PNVCLS(KIDIA:KFDIA),&
       & PFCLL(KIDIA:KFDIA,1),PFCLN(KIDIA:KFDIA,1),&
       & PFEVL(KIDIA:KFDIA,1),PFEVN(KIDIA:KFDIA,1),&
       & ZSSO_STDEV_, PSPSG(KIDIA:KFDIA),&
       & ZBUDTH_, ZBUDSO_,&
       & ZFCLL_, ZTOWNS_,&
       & ZCD_                         )

      CALL ARO_GROUND_DIAG_2ISBA( KBL, KGPCOMP, &
       & KFDIA, KIDIA, KFDIA, &
       & NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, &
       & PINDX(KIDIA:KFDIA), PINDY(KIDIA:KFDIA), &
       & PLSM, ZDUMMY1, ZDUMMY1, ZDUMMY1, ZTSURF(KIDIA:KFDIA), PSPSG(KIDIA:KFDIA), &
       & ZDUMMY1(KIDIA:KFDIA), ZDUMMY1(KIDIA:KFDIA), ZDUMMY1(KIDIA:KFDIA), &
       & ZDUMMY1(KIDIA:KFDIA), ZDUMMY1(KIDIA:KFDIA), ZDUMMY1(KIDIA:KFDIA), &
       & ZDUMMY1(KIDIA:KFDIA), PSPSGR(KIDIA:KFDIA), ZDUMMY1 )
 
    ENDIF

  ENDIF

  ENDDO SURFEX_LOOP

!* Compute PBL-diagnostics

   ZCAPE(:)=0._JPRB
   ZDCAPE(:)=0._JPRB   
   CALL ACCLDIA(YDXFU,YDPHY,YDMODEL%YRML_PHY_MF%YRPHY2,YDTOPH, KIDIA,KFDIA,KLON,KLEV,PUCLS,&
    & PVCLS,PUM(:,1:KLEV),PVM(:,1:KLEV), ZCAPE,ZDCAPE,ZTKEM(:,1:KLEV),PAPHIFM(:,1:KLEV),POROG,&
    & PUGST,PVGST,PPBLH,ICLPH)

   PPBLH(KIDIA:KFDIA)=MIN(XMAXLM,MAX(XMINLM,PPBLH(KIDIA:KFDIA))) 

   CALL ACVISIH(YDVISI,KIDIA,KFDIA,KLON,KTDIA,KLEV,PAPHIM,PAPHIFM,PAPRSFM,&
    & ZTM,ZRHM,ZQCM,ZQIM,ZQRM,ZQSM,ZQGM,PVISICLD, PVISIHYDRO,PMXCLWC)

ELSE

  ZSFSV_(KIDIA:KFDIA,:)=0._JPRB

ENDIF    !  <== End block "IF (LMSE)"

!*            IDEALIZED TURBULENT SURFACE FLUXES FOR SQUALL LINE CASE
!                --------------------------------------------------------

IF (LSQUALL.AND.LTURB) THEN
  ! on n'a besoin que d'un flux sur V (U est nul). 
  DO JLON=KIDIA,KFDIA
    IF (ABS(ZVM__(JLON,IKB)) <= 1.E-12) THEN
      ZSFV_(JLON)=0._JPRB
    ELSE
      ZSFV_(JLON)=-(ZVM__(JLON,IKB))**2 *&
       & (0.4_JPRB  /(LOG(ZZZ_F_(JLON,IKB)/0.2_JPRB) ) )**2&
       & *ZVM__(JLON,IKB)/ABS(ZVM__(JLON,IKB))  
    ENDIF
  ENDDO
ENDIF

!    ------------------------------------------------------------------
!    9.  Shallow Mass Flux Mixing
!    ------------------------------------------------------------------


IF (LMFSHAL) THEN
  IF (CMF_UPDRAFT=='DUAL') THEN
    ! Updraft computation from EDMF/ECMWF dual proposal
    ! Version May 2007
    !
    ! The following routine  are using arrays with the vertical Arpege/IFS fashion (as in the radiation scheme)

    IDRAFT = 2 ! beginning of the loop for MF tendency equation
               ! only 2 and 3 are used for tendency computation in ARO_SHALLOW_MF
    INDRAFT=3   ! 1 for test, 2 for dry, 3 for wet

    IF (KMAXDRAFT < INDRAFT) THEN
      CALL ABOR1('APL_AROME : KMAXDRAFT TOO SMALL !')
    ENDIF

    DO JLON = KIDIA, KFDIA
      ZZS_FTH_(JLON)=-1._JPRB*ZSFTH_(JLON)*(PAPRSM(JLON,KLEV)*ZINVATM)**(ZRSCP)
      ZZS_FRV_(JLON)=-1._JPRB*ZSFRV_(JLON)
    ENDDO
    ZZS_FU_(KIDIA:KFDIA)=ZSFU_(KIDIA:KFDIA)
    ZZS_FV_(KIDIA:KFDIA)=ZSFV_(KIDIA:KFDIA)

    !  IF LHARATU=TRUE then TKE at t-dt is needed as input for vdfexcuhl so fill ZTKEEDMF with t-1 value  from PTKEM

    IF (LHARATU) THEN
      DO JLEV=1,KLEV
        ZTKEEDMF(KIDIA:KFDIA,JLEV)=ZTKEM(KIDIA:KFDIA,JLEV)
        ZLENGTH_M(KIDIA:KFDIA,JLEV)=0.01_JPRB
        ZLENGTH_H(KIDIA:KFDIA,JLEV)=0.01_JPRB
      ENDDO
      IF (MAXVAL(ZTKEM(KIDIA:KFDIA,1:KLEV)) > 3300._JPRB) THEN
        DO JLEV=1, KLEV
          DO JLON = KIDIA, KFDIA
            IF (ZTKEM(JLON,JLEV) > 3300._JPRB) THEN
              WRITE (NULOUT,*) 'TKE > 3300 ! '
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    CALL VDFHGHTHL(YDMODEL%YRML_PHY_G%YRVDF,YDMODEL%YRML_PHY_SLIN%YREPHLI, &
     & YDMODEL%YRML_PHY_EC%YRECUMF,YDMODEL%YRML_PHY_EC%YREPHY,YDPARAR, &
     & KSTEP,KIDIA,KFDIA,KLON,KLEV,INDRAFT,&
     & PDT, PUM, PVM,&
     & ZTM,ZQVM,ZQCM,ZQIM,PCLFS,&
     & PAPRSM, PAPRSFM, ZAPHIFM,ZAPHIM,&
     & ZZS_FTH_,ZZS_FRV_,ZZS_FU_,ZZS_FV_,&
     & ZMF_UP,ZTHETAL_UP,ZQT_UP,ZTHTV_UP,ZQC_UP,ZQI_UP,&
     & ZU_UP, ZV_UP,&
     & PGP2DSPP, NGFL_EZDIAG, PEZDIAG, &
     & ZTENDQVUP,ZTENDTUP,ZSURFPREP,ZSURFSNOW, &
     & ZUPGENL,ZUPGENN, ZCLFR, &
     & ZLENGTH_M, ZLENGTH_H, ZTKEEDMF)

    !  tendtup, tendqvup  tendencies for non-conserved AROME
    !  variables due to updraft precipitation/snow (and its evaporation)
    DO JLEV = 2 ,KLEV
      DO JLON = KIDIA,KFDIA
        PTENDT(JLON,JLEV)=PTENDT(JLON,JLEV) + ZTENDTUP(JLON,JLEV)
        PTENDR(JLON,JLEV,1)=PTENDR(JLON,JLEV,1) + ZTENDQVUP(JLON,JLEV)
      ENDDO
    ENDDO
    
    IF(LTOTPREC)THEN
      !Add rain and snow tendencies from the sub-grid scheme to tendencies and sources,
      !at all vertical levels, instead of diagnosing only surface precip. 
      ZSURFPREP(KIDIA:KFDIA)=0.0_JPRB
      ZSURFSNOW(KIDIA:KFDIA)=0.0_JPRB
      DO JLEV= 1, KLEV
        DO JLON = KIDIA, KFDIA
          !Add rain and snow to sources:
          ZRS_(JLON,JLEV,3)=ZRS_(JLON,JLEV,3)+ZUPGENL(JLON,JLEV)
          ZRS_(JLON,JLEV,5)=ZRS_(JLON,JLEV,5)+ZUPGENN(JLON,JLEV)
          ZTHS__(JLON,JLEV)=ZTHS__(JLON,JLEV)+ZTENDTUP(JLON,JLEV)*(RATM/&
           & PAPRSFM(JLON,JLEV))**(RD/RCPD)
          !Update rain/snow tendencies:
          PTENDR(JLON,JLEV,3)=PTENDR(JLON,JLEV,3)+ZUPGENL(JLON,JLEV)
          PTENDR(JLON,JLEV,5)=PTENDR(JLON,JLEV,5)+ZUPGENN(JLON,JLEV)
        ENDDO
      ENDDO 
    ENDIF

  ELSE
    IDRAFT=3 ! only a wet updraft
    INDRAFT=1
    ZSURFPREP(KIDIA:KFDIA)=0._JPRB
    ZSURFSNOW(KIDIA:KFDIA)=0._JPRB
  ENDIF

  DO JDRAFT=IDRAFT,3

    ! No need to swapp because IN and OUT are never needed simultaneously

    !!! Call mass fluxes computations
    ! If CMF_UPDRAFT='DUAL', the updraft characteritics are already computed and will be passed as inputs of SHALLOW_MF
    ! if not, they will be computed in SHALLOW_MF itself (from Méso-NH type routines)

    ! JDRAFT=2 : dry updraft
    ! JDRAFT=3 : wet updraft

    IF (CMF_UPDRAFT=='DUAL') THEN
      ! Goes from one of the updraft from the IFS level world to the Méso-NH level world
      ! go from q to r)
      DO JLEV = 1,KLEV
        DO JLON = KIDIA,KFDIA
          ZMF_UP__(JLON,JLEV) = ZMF_UP(JLON,JLEV,JDRAFT)
          ZZU_UP_(JLON,JLEV) = ZU_UP(JLON,JLEV,JDRAFT)
          ZZV_UP_(JLON,JLEV) = ZV_UP(JLON,JLEV,JDRAFT)
          ZTHETAL_UP_(JLON,JLEV) = ZTHETAL_UP(JLON,JLEV,JDRAFT)
          ZTHETAV_UP_(JLON,JLEV) = ZTHTV_UP(JLON,JLEV,JDRAFT)
          ZRT_UP_(JLON,JLEV)  = ZQT_UP(JLON,JLEV,JDRAFT)/(1.-ZQT_UP(JLON,JLEV,JDRAFT))
          ZRC_UP_(JLON,JLEV)  = ZQC_UP(JLON,JLEV,JDRAFT)/(1.-ZQT_UP(JLON,JLEV,JDRAFT))
          ZRI_UP_(JLON,JLEV)  = ZQI_UP(JLON,JLEV,JDRAFT)/(1.-ZQT_UP(JLON,JLEV,JDRAFT))
        ENDDO
      ENDDO
      ZZW_UP_(KIDIA:KFDIA,1:IKT)=0._JPRB
      ZZFRAC_UP_(KIDIA:KFDIA,1:IKT)=0._JPRB
      IF (LHARATU) THEN
        DO JLEV = 1,KLEV
          DO JLON = KIDIA,KFDIA
            ZLENGTHM__(JLON,JLEV) = MAX(0.01_JPRB,ZLENGTH_M(JLON,JLEV))
            ZLENGTHH__(JLON,JLEV) = MAX(0.01_JPRB,ZLENGTH_H(JLON,JLEV))
            ! TKE should be bigger than a minimum value:
            ZTKEEDMFS(JLON,JLEV) = MAX(ZTKEEDMF(JLON,JLEV),PPTKEMIN)*ZINVDT
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    IF(MOD(KSTEP+1,NPRINTFR)==0) THEN
      WRITE(NULOUT,*)"apres surface zsfth zsfrv",ZSFTH_(NPTP),ZSFRV_(NPTP)
    ENDIF

    DO JLEV = 1, KLEV 
      ZRC_MF_(KIDIA:KFDIA,JLEV)=0._JPRB
      ZRI_MF_(KIDIA:KFDIA,JLEV)=0._JPRB
      ZCF_MF_(KIDIA:KFDIA,JLEV)=0._JPRB
    ENDDO

    IF (JDRAFT == IDRAFT) THEN
      ! Fill the sum at the first iteration
      ZARG_FLXZTHVMF_ => ZFLXZTHVMF_SUM__(:,1:KLEV)
    ELSE
      ! increment
      ZARG_FLXZTHVMF_ => ZFLXZTHVMF_(:,1:KLEV)
    ENDIF

    CALL ARO_SHALLOW_MF (KKL=IKL, KLON=KFDIA,KLEV=KLEV,KRR=NRR,KRRL=NRRL,&
     & KRRI=NRRI,KSV=NGFL_EXT,HMF_UPDRAFT=CMF_UPDRAFT, HMF_CLOUD=CMF_CLOUD,&
     & HFRAC_ICE=CFRAC_ICE_SHALLOW_MF,&
     & OMIXUV=LMIXUV, ONOMIXLG=.FALSE.,KSV_LGBEG=0,KSV_LGEND=0,&
     & KTCOUNT=KSTEP+1, PTSTEP=ZDT,&
     & PZZ=ZZZ_,PZZF=ZZZ_F_,&
     & PDZZF=ZDZZ_F_,&
     & PRHODJ=ZRHODJM__(:,1:KLEV),&
     & PRHODREF=ZRHODREFM__(:,1:KLEV),&
     & PPABSM=ZPABSM__(:,1:KLEV),&
     & PEXNM=ZEXNREFM_,&
     & PSFTH=ZSFTH_,PSFRV=ZSFRV_,&
     & PTHM=ZTHM__(:,1:KLEV),PRM=ZRM_,&
     & PUM=ZUM__(:,1:KLEV),PVM=ZVM__(:,1:KLEV),&
     & PTKEM=ZTKEM__(:,1:KLEV),PSVM=ZSVM_,&
     & PDUDT_MF=ZMFUS_,PDVDT_MF=ZMFVS_,&
     & PDTHLDT_MF=ZTHLS_,PDRTDT_MF=ZRTS_,&
     & PDSVDT_MF=ZSVXXX_,&
     & PSIGMF=ZSIGMF_,PRC_MF=ZRC_MF_,&
     & PRI_MF=ZRI_MF_,&
     & PCF_MF=ZCF_MF_,PFLXZTHVMF=ZARG_FLXZTHVMF_,&
     & PTHL_UP=ZTHETAL_UP_,PRT_UP= ZRT_UP_,&
     & PRV_UP=ZZRV_UP_,&
     & PRC_UP=ZRC_UP_,PRI_UP=ZRI_UP_,&
     & PU_UP=ZZU_UP_,PV_UP=ZZV_UP_,&
     & PTHV_UP=ZTHETAV_UP_,PW_UP=ZZW_UP_,&
     & PFRAC_UP=ZZFRAC_UP_,PEMF=ZMF_UP__(:,1:KLEV))

    IF (JDRAFT > IDRAFT) THEN
      ! Add increment
      ZFLXZTHVMF_SUM__(KIDIA:KFDIA,1:KLEV)=ZFLXZTHVMF_SUM__(KIDIA:KFDIA,1:KLEV)+ZFLXZTHVMF_(KIDIA:KFDIA,1:KLEV)
    ENDIF

    ! traitement des sorties pour repasser dans le monde Aladin

    IF ((CMF_CLOUD=='DIRE'.OR.CMF_CLOUD=='BIGA').AND.JDRAFT==3) THEN
      ! sauvegarde pour le schema de nuage
      DO JLEV = 1,KLEV
        PEZDIAG(KIDIA:KFDIA,JLEV,1)=ZRC_MF_(KIDIA:KFDIA,JLEV)
        PEZDIAG(KIDIA:KFDIA,JLEV,3)=ZRI_MF_(KIDIA:KFDIA,JLEV)
        PEZDIAG(KIDIA:KFDIA,JLEV,2)=ZCF_MF_(KIDIA:KFDIA,JLEV)
      ENDDO
    ENDIF
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        ZUS__(JLON,JLEV)=ZUS__(JLON,JLEV)+ZMFUS_(JLON,JLEV)
        ZVS__(JLON,JLEV)=ZVS__(JLON,JLEV)+ZMFVS_(JLON,JLEV)
        ZTHS__(JLON,JLEV)=ZTHS__(JLON,JLEV)+ZTHLS_(JLON,JLEV)
        ZRS_(JLON,JLEV,1)=ZRS_(JLON,JLEV,1)+ZRTS_(JLON,JLEV)
        !calcul de tendance et inversion des niveaux pour le vent horizontal
        PTENDU(JLON,JLEV)=PTENDU(JLON,JLEV)+ZMFUS_(JLON,JLEV)
        PTENDV(JLON,JLEV)=PTENDV(JLON,JLEV)+ZMFVS_(JLON,JLEV)
        !conversion de la tendance de theta en tendance de T et inversion niveau
        PTENDT(JLON,JLEV)=PTENDT(JLON,JLEV)+ZTHLS_(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
        ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+ZTHLS_(JLON,JLEV)
        !inversion niveaux tendances des ri et conversion en qi en multipliant par qd
        PTENDR(JLON,JLEV,1) = PTENDR(JLON,JLEV,1)+ZRTS_(JLON,JLEV)*ZQDM(JLON,JLEV)
      ENDDO
    ENDDO

  ENDDO ! JDRAFT

ENDIF ! LMFSHAL


!    ------------------------------------------------------------------
!     10 - TURBULENCE.
!     --------------------------------------------------------------------

IF (LTURB) THEN

  ! Swapp because IN and OUT might be needed simultaneously (though commented out)
  CALL SWAP_LIMAS

  ! Swapp because IN and OUT will be needed simultaneously
  CALL SWAP_SVS
  ! well let's keep the copy, though for now  OUT=IN anyway.
  IF (NGFL_EXT /=0 ) THEN
    ZSVSIN_(KIDIA:KFDIA,1:KLEV,1:NGFL_EXT)=ZSVS_(KIDIA:KFDIA,1:KLEV,1:NGFL_EXT)
  ENDIF

  !prints
  IF(MOD(KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome U'
    WRITE(NULOUT,*)MAXVAL(ZUM__(:,IKB)), MINVAL(ZUM__(:,IKB))
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome V'
    WRITE(NULOUT,*)MAXVAL(ZVM__(:,IKB)), MINVAL(ZVM__(:,IKB))
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome W'
    WRITE(NULOUT,*)MAXVAL(ZWM__(:,IKB)), MINVAL(ZWM__(:,IKB))
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome TKE'
    WRITE(NULOUT,*)MAXVAL(ZTKEM__(:,IKB)), MINVAL(ZTKEM__(:,IKB))
    DO JLEV=1,KLEV
      WRITE(NULOUT,*)JLEV,ZUM__(NPTP,JLEV),ZVM__(NPTP,JLEV),ZWM__(NPTP,JLEV),ZTKEM__(NPTP,JLEV)
    ENDDO
    WRITE(NULOUT,*)'u v w tke a S'
    DO JLEV=1,KLEV
      WRITE(NULOUT,*)JLEV,ZUS__(NPTP,JLEV),ZVS__(NPTP,JLEV),ZWS__(NPTP,JLEV),ZTKES_(NPTP,JLEV)
    ENDDO
    WRITE(NULOUT,*)'ZTHS__ avant turb'
    DO JLEV=1,KLEV
      WRITE(NULOUT,*)JLEV,ZTHS__(NPTP,JLEV)
    ENDDO
  ENDIF

!!$
!!$! Allocation des variables SV (NGFL_EXT + NLIMA)
!!$  KSV_TURB=NGFL_EXT+NLIMA
!!$!
!!$  IF (NGFL_EXT/=0) THEN
!!$     DO JGFL=1,NGFL_EXT
!!$        DO JLON=KIDIA,KFDIA
!!$           ZSFTURB(JLON,JGFL)=ZSFSV_(JLON,JGFL)
!!$           DO JLEV = 1, KLEV
!!$              ZTURBM(JLON,JLEV,JGFL)=ZSVM_(JLON,1,JLEV,JGFL)
!!$              ZTURBS(JLON,JLEV,JGFL)=ZSVSIN_(JLON,1,JLEV,JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF
!!$!
!!$  IF (NLIMA/=0) THEN
!!$     DO JGFL=1,NLIMA
!!$        DO JLON=KIDIA,KFDIA
!!$           ZSFTURB(JLON,NGFL_EXT+JGFL)=0.
!!$           DO JLEV = 1, KLEV
!!$              ZTURBM(JLON,JLEV,NGFL_EXT+JGFL)=ZLIMAM_(JLON,JLEV,JGFL)
!!$              ZTURBS(JLON,JLEV,NGFL_EXT+JGFL)=ZLIMASIN_(JLON,JLEV,JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF

  ! Input variable indeed. REK
  ZSFSVLIMA_(KIDIA:KFDIA,1:NLIMA)=0._JPRB

  ! 10.2 calcul TURB
  ZZTOP_(KIDIA:KFDIA)=ZAPHIM(KIDIA:KFDIA,0)*ZINVG

  IF (LGRADHPHY) THEN
  !   
    DO JLEV = 1,KLEV
      DO JGR=1,NGRADIENTS
        ZTURB3D__(KIDIA:KFDIA,JLEV,JGR)=PTURB3D(KIDIA:KFDIA,JGR,JLEV)
      ENDDO
    ENDDO
  
  ENDIF

! Appel avec les arguments modifiés pour variables LIMA :
! KSV_TURB, ZSFTURB, ZTURBM, ZTURBS, ZTENDSV_TURB
  CALL ARO_TURB_MNH(KKA=IKA,KKU=IKU,KKL=IKL,KLON=KFDIA,KLEV=KLEV,&
   & KRR=NRR, KRRL=NRRL,KRRI= NRRI,&
   & KSV=NLIMA,KTCOUNT=KSTEP+1,KGRADIENTS=NGRADIENTS,&
   &LDHARATU=LHARATU,CMICRO=CMICRO,PTSTEP=ZDT,&
   & PZZ=ZZZ_,PZZF=ZZZ_F_,&
   & PZZTOP=ZZTOP_,PRHODJ=ZRHODJM__,PTHVREF=ZTHVREFM__,&
   & PRHODREF=ZRHODREFM__,HINST_SFU='M',&
   & HMF_UPDRAFT=CMF_UPDRAFT,&
   & PSFTH=ZSFTH_,PSFRV=ZSFRV_,&
   & PSFSV=ZSFSVLIMA_,PSFU=ZSFU_,&
   & PSFV=ZSFV_,PPABSM=ZPABSM__,&
   & PUM=ZUM__,PVM=ZVM__,PWM=ZWM__,PTKEM=ZTKEM__,PEPSM=ZEPSM,&
   & PSVM=ZLIMAM_,&
   & PSRCM=ZSRCS__,PTHM=ZTHM__,&
   & PRM=ZRM_,&
   & PRUS=ZUS__,PRVS=ZVS__,PRWS=ZWS__,PRTHS=ZTHS__,&
   & PRRS=ZRS_,&
   & PRSVSIN=ZLIMASIN_,PRSVS=ZLIMAS_,&
   & PRTKES=ZTKES_,PRTKES_OUT=ZTKES_OUT__,&
   & PREPSS=ZEPSS,PHGRAD=ZTURB3D__,PSIGS=ZSIGS__,&
   & OSUBG_COND=LOSUBG_COND,&
   & PFLXZTHVMF=ZFLXZTHVMF_SUM__,&
   & PLENGTHM=ZLENGTHM__,PLENGTHH=ZLENGTHH__,MFMOIST=ZMF_UP__,PDRUS_TURB=ZTENDU_TURB__, &
   & PDRVS_TURB=ZTENDV_TURB__,PDRTHLS_TURB=ZTENDTHL_TURB__,PDRRTS_TURB=ZTENDRT_TURB__,&
   & PDRSVS_TURB=ZTENDSV_TURBLIMA_,&
   & PDP=ZDP__, PTP=ZTP__,PTPMF=ZTPMF__,PTDIFF=ZTDIFF__,PTDISS=ZTDISS__,PEDR=ZEDR__,YDDDH=YDDDH,&
   & YDLDDH=YDMODEL%YRML_DIAG%YRLDDH,YDMDDH=YDMODEL%YRML_DIAG%YRMDDH)


! Séparation des variables SV (NGFL_EXT + NLIMA)
!!$  IF (NGFL_EXT/=0) THEN
!!$     DO JGFL=1,NGFL_EXT
!!$        DO JLON=KIDIA,KFDIA
!!$           ZSFSV_(JLON,JGFL)=ZSFTURB(JLON,JGFL)
!!$           DO JLEV = 1, KLEV
!!$              ZSVM_(JLON,1,JLEV,JGFL)=ZTURBM(JLON,JLEV,JGFL)
!!$              ZSVS_(JLON,1,JLEV,JGFL)=ZTURBS(JLON,JLEV,JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF
!!$!
!!$  IF (NLIMA/=0) THEN
!!$     DO JGFL=1,NLIMA
!!$        DO JLON=KIDIA,KFDIA
!!$           DO JLEV = 1, KLEV
!!$              ZLIMAM_(JLON,JLEV,JGFL)=ZTURBM(JLON,JLEV,NGFL_EXT+JGFL)
!!$              ZLIMAS_(JLON,JLEV,JGFL)=ZTURBS(JLON,JLEV,NGFL_EXT+JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF


  DO JLEV = 1 , KLEV
     PEDR(KIDIA:KFDIA,JLEV)=ZEDR__(KIDIA:KFDIA,JLEV)
  ENDDO
   
  IF (LFLEXDIA) THEN
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        ZDP__(JLON,JLEV)=ZDP__(JLON,JLEV)*PDELPM(JLON,JLEV)*ZINVG
        ZTP__(JLON,JLEV)=(ZTP__(JLON,JLEV)-ZTPMF__(JLON,JLEV))*PDELPM(JLON,JLEV)*ZINVG
        ZTPMF__(JLON,JLEV)=ZTPMF__(JLON,JLEV)*PDELPM(JLON,JLEV)*ZINVG
        ZTDIFF__(JLON,JLEV)=ZTDIFF__(JLON,JLEV)*PDELPM(JLON,JLEV)*ZINVG
        ZTDISS__(JLON,JLEV)=ZTDISS__(JLON,JLEV)*PDELPM(JLON,JLEV)*ZINVG
      ENDDO
    ENDDO
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZDP__(:,1:KLEV),'TKEPRDY',YDDDH)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTP__(:,1:KLEV),'TKEPRTH',YDDDH)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTPMF__(:,1:KLEV),'TKEPRTHMF',YDDDH)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTDIFF__(:,1:KLEV),'TKEDIFF',YDDDH)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTDISS__(:,1:KLEV),'TKEDISS',YDDDH)
    ELSE
      CALL ADD_FIELD_3D(YLDDH,ZDP__(:,1:KLEV),'TKEPRDY','T','ARO',.TRUE.,.TRUE.)
      CALL ADD_FIELD_3D(YLDDH,ZTP__(:,1:KLEV),'TKEPRTH','T','ARO',.TRUE.,.TRUE.)
      CALL ADD_FIELD_3D(YLDDH,ZTPMF__(:,1:KLEV),'TKEPRTHMF','T','ARO',.TRUE.,.TRUE.)
      CALL ADD_FIELD_3D(YLDDH,ZTDIFF__(:,1:KLEV),'TKEDIFF','T','ARO',.TRUE.,.TRUE.)
      CALL ADD_FIELD_3D(YLDDH,ZTDISS__(:,1:KLEV),'TKEDISS','T','ARO',.TRUE.,.TRUE.)
    ENDIF
  ENDIF 

  IF(MOD(KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'u v w a S apres turb'
    DO JLEV=1,KLEV
      WRITE(NULOUT,*)JLEV,ZUS__(NPTP,JLEV),ZVS__(NPTP,JLEV),ZWS__(NPTP,JLEV),ZTKES_OUT__(NPTP,JLEV)
    ENDDO
    WRITE(NULOUT,*)'THS TKES SIGS apres turb'
    DO JLEV=1,KLEV
      WRITE(NULOUT,*)JLEV,ZTHS__(NPTP,JLEV),ZTKES_OUT__(NPTP,JLEV),ZSIGS__(NPTP,JLEV)
    ENDDO
  ENDIF

  ! avance temporelle et inversion niveau pour ZSIGS__
  IF (LOSUBG_COND .AND. LOSIGMAS) THEN
    IF (CMF_CLOUD=='DIRE'.OR.CMF_CLOUD=='BIGA'.OR.CMF_CLOUD=='NONE') THEN
      DO JLEV = 1,KLEV
        PSIGS(KIDIA:KFDIA,JLEV)=ZSIGS__(KIDIA:KFDIA,JLEV)
      ENDDO
    ELSEIF (CMF_CLOUD=='STAT') THEN
      DO JLEV = 1,KLEV
        DO JLON = KIDIA,KFDIA
          PSIGS(JLON,JLEV)=SQRT(ZSIGS__(JLON,JLEV)**2+ZSIGMF_(JLON,JLEV)**2 )
        ENDDO
      ENDDO
    ENDIF
  ENDIF


  !10.3. traitement des sorties pour repasser dans le monde Aladin
  !calcul de tendance et inversion des niveaux pour le vent horizontal et la TKE

  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      PTENDU(JLON,JLEV)=PTENDU(JLON,JLEV)+ZTENDU_TURB__(JLON,JLEV)
      PTENDV(JLON,JLEV)=PTENDV(JLON,JLEV)+ZTENDV_TURB__(JLON,JLEV)
      ! for the moment, turbulence do not compute w tendency:
      PTENDW(JLON,JLEV)=0.0_JPRB
      ! PTENDW(JLON,JLEV)+(ZWS__(JLON,JLEV)-&
      ! & ZWS_AVE(JLON,1,JLEV))
      !conversion de la tendance de theta en tendance de T et inversion niveau
      PTENDT(JLON,JLEV)=PTENDT(JLON,JLEV)+ZTENDTHL_TURB__(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
  !inversion niveaux tendances des rv et conversion en qv en multipliant par qd
      PTENDR(JLON,JLEV,1)= PTENDR(JLON,JLEV,1)+ZTENDRT_TURB__(JLON,JLEV)*ZQDM(JLON,JLEV)
    ENDDO
  ENDDO


  IF (LHARATU) THEN
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
          PTENDTKE(JLON,JLEV)=PTENDTKE(JLON,JLEV)+(ZTKEEDMFS(JLON,JLEV)-ZTKES_(JLON,JLEV))
      ENDDO
    ENDDO
  ELSE
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
         PTENDTKE(JLON,JLEV)=PTENDTKE(JLON,JLEV)+(ZTKES_OUT__(JLON,JLEV)-ZTKES_(JLON,JLEV))
      ENDDO
    ENDDO
  ENDIF

  DO JGFL=1,NGFL_EXT
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        PTENDEXT(JLON,JLEV,JGFL)=PTENDEXT(JLON,JLEV,JGFL)+(ZSVS_(JLON,JLEV,JGFL)-ZSVSIN_(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO

! Tendances LIMA
  DO JGFL=1,NLIMA
    DO JLEV = 1, KLEV
      DO JLON=KIDIA,KFDIA
        PTENDLIMA(JLON,JLEV,JGFL)=PTENDLIMA(JLON,JLEV,JGFL)+ZTENDSV_TURBLIMA_(JLON,JLEV,JGFL)
!        PTENDLIMA(JLON,JLEV,:)=PTENDLIMA(JLON,JLEV,:)+ (ZLIMAS_(JLON,JLEV,:)-ZLIMASIN_(JLON,JLEV,:))
      ENDDO
    ENDDO
  ENDDO

ENDIF
!    ------------------------------------------------------------------
!     11 - MICROPHYSIQUE. 
!     --------------------------------------------------------------------

IF (LMICRO) THEN

  ! Swap pointers because input values of THS and RS should be saved
  CALL SWAP_THS
  CALL SWAP_RS
  CALL SWAP_LIMAS

  ! for now a copy is needed (see below, inside). I don't like than :-( REK
  ZTHS__(KIDIA:KFDIA,1:KLEV)=ZTHSIN_(KIDIA:KFDIA,1:KLEV)
  ZRS_(KIDIA:KFDIA,1:KLEV,1:NRR)=ZRSIN_(KIDIA:KFDIA,1:KLEV,1:NRR)
  ! for now a copy is needed (see below, inside). I don't like than :-( REK
  ZLIMAS_(KIDIA:KFDIA,1:KLEV,1:NLIMA)=ZLIMASIN_(KIDIA:KFDIA,1:KLEV,1:NLIMA)
     
  !prints
  IF (MOD(KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'avant rain_ice sous apl_arome'
    WRITE(NULOUT,*)'JLEV   ZZZ_F_      ZZZ_      ZRHODREF',&
     & '    ZRHODJ      ZPABSM__        ZTHSIN_       ZTHM__      '   
    DO JLEV=1,KLEV+1 
      WRITE(NULOUT,'(I2,X,7F10.3)')JLEV,ZZZ_F_(NPTP,JLEV),ZZZ_(NPTP,JLEV), ZRHODREFM__(NPTP,JLEV),&
       & ZRHODJM__(NPTP,JLEV), ZPABSM__(NPTP,JLEV), ZTHSIN_(NPTP,JLEV), ZTHM__(NPTP,JLEV)  
    ENDDO 
    WRITE(NULOUT,*)'JLEV        PDELPM        ZPABSM__         ZEXNREF','          ZSIGS__'
    DO JLEV=2,KLEV
      WRITE(NULOUT,'(I2,X,4f10.3)')JLEV, PDELPM(NPTP,JLEV),&
       & ZPABSM__(NPTP,JLEV),ZEXNREFM_(NPTP,JLEV),ZSIGS__(NPTP,JLEV)  
    ENDDO
    WRITE(NULOUT,*)'JLEV    PTM       PRM          PCPM'
    DO JLEV=1,KLEV
      WRITE(NULOUT,'(I2,X,3f10.3)')JLEV,ZTM(NPTP,KLEV+1-JLEV), ZRHM(NPTP,KLEV+1-JLEV) ,ZCPM(NPTP,KLEV+1-JLEV)  
    ENDDO
    WRITE (NULOUT,*)'JLEV  rhoQv  rhoQc   rhoQr   rhoQi   rhoQs   rhoQg'
    DO JLEV=1,KLEV
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,ZRM_(NPTP,JLEV,1), ZRM_(NPTP,JLEV,2),&
       & ZRM_(NPTP,JLEV,3),ZRM_(NPTP,JLEV,4),ZRM_(NPTP,JLEV,5), ZRM_(NPTP,JLEV,6)  
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZRSQv  ZRSQc   ZRSQr   ZRSQi   ZRSQs   ZRSQg'
    DO JLEV=1,KLEV
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,ZRS_(NPTP,JLEV,1), ZRS_(NPTP,JLEV,2),&
       & ZRSIN_(NPTP,JLEV,3),ZRSIN_(NPTP,JLEV,4),ZRSIN_(NPTP,JLEV,5), ZRSIN_(NPTP,JLEV,6)  
    ENDDO
    WRITE(NULOUT,*)'ZDT=',ZDT
    WRITE(NULOUT,*)'NRR and co',NRR,KSTEP+1,NSPLITR,LOSUBG_COND, LOSIGMAS, CSUBG_AUCV_RC,LOWARM  
  ENDIF
  

  ZSEA_(KIDIA:KFDIA)=0.0_JPRB
  IF (LOLSMC) THEN
    DO JLON = KIDIA, KFDIA
      IF (PLSM(JLON) < 0.5) THEN
        ZSEA_(JLON) = 1.0_JPRB
      ENDIF
    ENDDO
  ENDIF
         
  IF (LOTOWNC) THEN
    ZTOWN_(KIDIA:KFDIA) = ZTOWNS_(KIDIA:KFDIA)
  ELSE
    ZTOWN_(KIDIA:KFDIA)=0.0_JPRB
  ENDIF  

  IF (CMICRO == 'LIMA') THEN

    IF (LTURB) THEN
      DO JLON=KIDIA,KFDIA
        DO JLEV=1,KLEV
          ZWNU_(JLON,JLEV) = ZWM__(JLON,JLEV) + 0.66*SQRT(ZTKEM__(JLON,JLEV))
        ENDDO
      ENDDO
      ZPTRWNU_ => ZWNU_(1:KFDIA,1:KLEV)
    ELSE
      ZPTRWNU_ => ZWM__(1:KFDIA,1:KLEV)
    ENDIF
    CALL ARO_LIMA(KLEV,IKU,IKL,KFDIA,KLEV,NRR,NLIMA,KSTEP+1,NSPLITR,NSPLITG, ZDT,ZDZZ_ ,&
     & ZRHODJM__(:,1:KLEV),ZRHODREFM__(:,1:KLEV), ZEXNREFM_, ZPABSM__(:,1:KLEV), ZPTRWNU_, &
     & ZTHM__(:,1:KLEV),ZRM_, ZLIMAM_,  ZTHS__(:,1:KLEV),ZRS_,  ZLIMAS_,  ZEVAP_, &
     & ZINPRR_NOTINCR_, ZINPRS_NOTINCR_, ZINPRG_NOTINCR_, ZINPRH_NOTINCR_,ZPFPR_,&
     & YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH )
  ELSE
    !CALL ARO_RAIN_ICE (NPROMICRO,KLEV,IKU,IKL,KLON,KLEV,KFDIA,NRR,KSTEP+1,NSPLITR,NGFL_EZDIAG,&      !this is the target version
    CALL ARO_RAIN_ICE (NPROMICRO,KLEV,IKU,IKL,KFDIA,KLEV,KFDIA,NRR,KSTEP+1,NSPLITR,NGFL_EZDIAG,&
     & LOSUBG_COND, CSUBG_AUCV_RC, CSUBG_AUCV_RI, LOSEDIC,CSEDIM, CMICRO, ZDT,ZDZZ_ ,&
     & ZRHODJM__(:,1:KLEV),ZRHODREFM__(:,1:KLEV), ZEXNREFM_, ZPABSM__(:,1:KLEV),&
     & ZHLC_HRC__(:,1:KLEV), ZHLC_HCF__(:,1:KLEV),&
     & ZHLI_HRI__(:,1:KLEV), ZHLI_HCF__(:,1:KLEV),&
     & ZTHM__(:,1:KLEV),ZRM_, ZSIGS__(:,1:KLEV), ZNEBMNH_, ZTHS__(:,1:KLEV),ZRS_,&
     & ZEVAP_, ZCIT_,LOWARM,ZSEA_,ZTOWN_, LOCND2,LGRSN,&
     & ZINPRR_NOTINCR_, ZINPRS_NOTINCR_, ZINPRG_NOTINCR_, ZINPRH_NOTINCR_,ZPFPR_,&
     & PGP2DSPP,PEZDIAG, &
     & YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH)
  ENDIF

  DO JLON=KIDIA,KFDIA
    ZINPRR_(JLON)=ZINPRR_(JLON)+ZINPRR_NOTINCR_(JLON)
    ZINPRS_(JLON)=ZINPRS_(JLON)+ZINPRS_NOTINCR_(JLON)
    ZINPRG_(JLON)=ZINPRG_(JLON)+ZINPRG_NOTINCR_(JLON)
    ZINPRH_(JLON)=ZINPRH_(JLON)+ZINPRH_NOTINCR_(JLON)
  ENDDO

  !conversion de la tendance de theta en tendance de T et inversion niveau
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      PTENDT(JLON,JLEV)= PTENDT(JLON,JLEV)+(ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV))*ZEXNREFM_(JLON,JLEV)  
      ZTENDT(JLON,JLEV)= ZTENDT(JLON,JLEV)+ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV)
    ENDDO
  ENDDO
  
  !inversion niveaux tendances des ri et conversion en qi en multipliant par qd
  DO JR=1,NRR
    DO JLEV=1,KLEV
      DO JLON=KIDIA,KFDIA
        PTENDR(JLON,JLEV,JR)=PTENDR(JLON,JLEV,JR)+(ZRS_(JLON,JLEV,JR)-ZRSIN_(JLON,JLEV,JR))*ZQDM(JLON,JLEV)  
      ENDDO
    ENDDO
  ENDDO

  ! Tendances des variables LIMA
  DO JGFL=1,NLIMA
    DO JLEV = 1, KLEV
      DO JLON=KIDIA,KFDIA
        PTENDLIMA(JLON,JLEV,JGFL)=PTENDLIMA(JLON,JLEV,JGFL)+(ZLIMAS_(JLON,JLEV,JGFL)-ZLIMASIN_(JLON,JLEV,JGFL))  
      ENDDO
    ENDDO
  ENDDO

  IF (LINTFLEX) THEN
    !inversion of levels of upper-air precipitation
    DO JR=2,NRR ! no precip for qv
      ZFPR(KIDIA:KFDIA,0,JR)=0._JPRB  ! zero precip at top of atmosphere
      DO JLEV=1,KLEV
        ZFPR(KIDIA:KFDIA,JLEV,JR)=ZPFPR_(KIDIA:KFDIA,JLEV,JR)
      ENDDO
    ENDDO
  ENDIF

  !store cumulative 3D precipitations for mocage      
  IF (LFPREC3D) THEN
    DO JR=2,NRR ! no precip for qv
      DO JLEV=1,KLEV
        DO JLON=KIDIA,KFDIA
          PEZDIAG(JLON,JLEV,4)=PEZDIAG(JLON,JLEV,4)+ZPFPR_(JLON,JLEV,JR)*1000._JPRB*PDT
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  !prints    
  IF(MOD(KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'PTENDT en sortie de rain_ice'
    WRITE(NULOUT,*)'ZTHS__ en sortie de rain_ice'
    DO JLEV=1,KLEV
      WRITE(NULOUT,*)PTENDT(NPTP,JLEV),ZTHS__(NPTP,JLEV)
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZTENDQv  ZTZNDQc   ZTENDQr   ZTENDQi' ,'ZTENDQs   ZTENDQg'  
    DO JLEV=1,KLEV
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,PTENDR(NPTP,JLEV,1),PTENDR(NPTP,JLEV,2),&
       & PTENDR(NPTP,JLEV,3),PTENDR(NPTP,JLEV,4),PTENDR(NPTP,JLEV,5),PTENDR(NPTP,JLEV,6)  
    ENDDO
    WRITE (NULOUT,*) 'ZSRCS__ et ZNEBMNH_',MAXVAL(ZSRCS__),MAXVAL(ZNEBMNH_) 
  ENDIF
  
  IF (LRDEPOS) THEN
    ISPLITR=NSPLITR
    ! Swapp because IN and OUT will be needed simultaneously
    CALL SWAP_SVM
    CALL ARO_RAINAERO(KFDIA,KLEV,NGFL_EXT,NRR, PDT,ZSVMIN_,ZZZ_,&
     & ZPABSM__(:,1:KLEV),ZTHM__(:,1:KLEV),ZRHODREFM__(:,1:KLEV),&
     & KSTEP+1,ZRM_,ZEVAP_, ISPLITR, ZSVM_           )
    ! return to tendency
    DO JGFL=1,NGFL_EXT
      DO JLEV = 1,KLEV
        DO JLON = KIDIA,KFDIA
          PTENDEXT(JLON,JLEV,JGFL)=PTENDEXT(JLON,JLEV,JGFL)+(ZSVM_(JLON,JLEV,JGFL)-ZSVMIN_(JLON,JLEV,JGFL))*ZINVDT
        ENDDO
      ENDDO
    ENDDO
  ENDIF ! LRDEPOS

ENDIF ! LMICRO

! start LHN F.Meier 2020 ******

LNUDGLHNREAD=.TRUE.
IF(MYPROC==1.AND.KSTEP==1.AND.LNUDGLH)THEN
  CALL NUDGLHCLIMPROF(KLEV,LNUDGLHNREAD)
ENDIF
! save accumulated precipitation for LHN
IF (LNUDGLH.AND.KSTEP == NSTARTNUDGLH.AND.NSTARTNUDGLH > 0) THEN
  !IF(MYPROC==1) WRITE(NULOUT,*)'save precip for LHN - STEP:',KSTEP, &
  !  & 'NUDGINGINT:',NINTNUDGLH,'NSTARTNUDGLH:',NSTARTNUDGLH
  CALL NUDGLHPRECIP(KIDIA,KFDIA,KLON,ZACPRR_(KIDIA:KFDIA),&
         & ZACPRS_(KIDIA:KFDIA),ZACPRG_(KIDIA:KFDIA),KBL)
ENDIF
ISTEP=KSTEP-NSTARTNUDGLH
! if LNUDGLH and KSTEP in nudging interval
IF (LNUDGLH.AND.KSTEP > NSTARTNUDGLH.AND.KSTEP <= NSTOPNUDGLH) THEN
  ! safe LH profile for step before LHN step
  LLHN=.FALSE.
  IF(MOD(ISTEP+1,NINTNUDGLH)==0) THEN
    CALL NUDGLHPREP(KIDIA,KFDIA,KLON,KLEV,ZTENDT,KBL)
  ENDIF
  ! LHN step
  IF(MOD(ISTEP,NINTNUDGLH)==0) THEN
    !IF(MYPROC==1) WRITE(NULOUT,*)'LH nudging applied - STEP:',KSTEP, &
    !  & 'NUDGINGINT:',NINTNUDGLH
    ! get index for correctly reading observation from array
    ! first two indices are reserved for other LHN stuff
    JLHSTEP=NINT(1.0_JPRB*ISTEP/(NTIMESPLITNUDGLH*NINTNUDGLH))+2
    !IF(MYPROC==1) WRITE(NULOUT,*)'observation array:',JLHSTEP
    ! call nudging routine to modify LHN profile where necessary
    CALL NUDGLH(NGPTOT,NPROMA,NGPBLKS,KIDIA,KFDIA,KLON,KLEV,ZACPRR_(KIDIA:KFDIA),&
       & ZACPRS_(KIDIA:KFDIA),ZACPRG_(KIDIA:KFDIA),&
       & ZTENDT(KIDIA:KFDIA,1:KLEV),JLHSTEP,KBL,&
       & ZEXNREFM_(KIDIA:KFDIA,1:KLEV),.TRUE.,&
       & LLHN(KIDIA:KFDIA,1:KLEV),ZPABSM__(KIDIA:KFDIA,1:KLEV),ZDT,&
       & ZTHM__(KIDIA:KFDIA,1:KLEV),ZRM_(KIDIA:KFDIA,1:KLEV,:),&
       & ZQDM(KIDIA:KFDIA,1:KLEV),PTENDR(KIDIA:KFDIA,1:KLEV,1),NRR,LNUDGLHNREAD)
    !IF(MYPROC==1) WRITE(NULOUT,*)'calling LH successful - convert TH to T and
    !add temperature tendency'
    ! add LHN tendency to physics tendency, limit LHN tendency
    ZMAXTEND=0.0_JPRB
    ZMINTEND=0.0_JPRB
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        IF(LLHN(JLON,JLEV))THEN
          ZTENDT(JLON,JLEV)=MAX(ZTENDT(JLON,JLEV),RMINNUDGLH)
          ZTENDT(JLON,JLEV)=MIN(ZTENDT(JLON,JLEV),RMAXNUDGLH)
          PTENDT(JLON,JLEV)=PTENDT(JLON,JLEV)+ZTENDT(JLON,JLEV)*&
           & RAMPLIFY*ZEXNREFM_(JLON,JLEV)
          ZMINTEND=MIN(ZTENDT(JLON,JLEV),ZMINTEND)
          ZMAXTEND=MAX(ZTENDT(JLON,JLEV),ZMAXTEND)
          ! keep RH constant if LNUDGLHCOMPT=T
          IF(PTM(JLON,JLEV)>0.01_JPRB.AND.LNUDGLHCOMPT)THEN
            PTENDR(JLON,JLEV,1)=PTENDR(JLON,JLEV,1)+RLVTT/RV/((PTM(JLON,JLEV))**2._JPRB)* &
            & ZTENDT(JLON,JLEV)*RAMPLIFY*ZEXNREFM_(JLON,JLEV)*ZQSAT(JLON,JLEV)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !IF(ZMINTEND<-0.01)WRITE(*,*)'ZMINTEND',ZMINTEND
    !IF(ZMAXTEND>0.01)WRITE(*,*)'ZMAXTEND',ZMAXTEND
    ! write LH profiles to array to save it for next time step
    CALL NUDGLHPREP(KIDIA,KFDIA,KLON,KLEV,ZTENDT,KBL)
    IF(MYPROC==1) WRITE(NULOUT,*)'calling LH successful finished'
    ! use LHN factor again on following time steps depending on NTAUNUDGLH
  ELSEIF(MOD(ISTEP,NINTNUDGLH)<NTAUNUDGLH.AND.MOD(ISTEP,NINTNUDGLH)>0 &
      & .AND.ISTEP>NINTNUDGLH) THEN
    IF(MYPROC==1)THEN
      WRITE(NULOUT,*)'LH nudging applied-STEP:',KSTEP,'NUDGINGINT:',NINTNUDGLH
      WRITE(NULOUT,*)'NTAUNUDGLH:',NTAUNUDGLH
    ENDIF
    ! get index for reading correctly most recent obs
    JLHSTEP=2+NINT(1.0_JPRB*(ISTEP-MOD(ISTEP,NINTNUDGLH))/(NTIMESPLITNUDGLH*NINTNUDGLH))
    !IF(MYPROC==1) WRITE(NULOUT,*)'observation array:',JLHSTEP
    ! call nudging routine to modify LHN profile where necessary
    ! LHN factor is not recalculated but might be damped by RDAMPNUDGLH
    CALL NUDGLH(NGPTOT,NPROMA,NGPBLKS,KIDIA,KFDIA,KLON,KLEV,ZACPRR_(KIDIA:KFDIA),&
        & ZACPRS_(KIDIA:KFDIA),ZACPRG_(KIDIA:KFDIA),&
        & ZTENDT(KIDIA:KFDIA,1:KLEV),JLHSTEP,KBL,&
        & ZEXNREFM_(KIDIA:KFDIA,1:KLEV),.FALSE.,&
        & LLHN(KIDIA:KFDIA,1:KLEV),ZPABSM__(KIDIA:KFDIA,1:KLEV),ZDT,&
        & ZTHM__(KIDIA:KFDIA,1:KLEV),ZRM_(KIDIA:KFDIA,1:KLEV,:),&
        & ZQDM(KIDIA:KFDIA,1:KLEV),PTENDR(KIDIA:KFDIA,1:KLEV,1),NRR,LNUDGLHNREAD)
    !IF(MYPROC==1) WRITE(NULOUT,*)'calling LH successful - convert TH to T and
    !add temperature tendency'
    ! add LHN tendency to physics tendency, limit LHN tendency
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        IF(LLHN(JLON,JLEV))THEN
          ZTENDT(JLON,JLEV)=MAX(ZTENDT(JLON,JLEV),RMINNUDGLH)
          ZTENDT(JLON,JLEV)=MIN(ZTENDT(JLON,JLEV),RMAXNUDGLH)
          PTENDT(JLON,JLEV)= PTENDT(JLON,JLEV)+ZTENDT(JLON,JLEV)*&
           & RAMPLIFY*ZEXNREFM_(JLON,JLEV)
          ZMINTEND=MIN(ZTENDT(JLON,JLEV),ZMINTEND)
          ZMAXTEND=MAX(ZTENDT(JLON,JLEV),ZMAXTEND)
          ! keep RH constant if LNUDGLHCOMPT=T
          IF(PTM(JLON,JLEV)>0.01_JPRB.AND.LNUDGLHCOMPT)THEN
             PTENDR(JLON,JLEV,1)=PTENDR(JLON,JLEV,1)+RLVTT/RV/((PTM(JLON,JLEV))**2._JPRB)*&
              & ZTENDT(JLON,JLEV)*RAMPLIFY*ZEXNREFM_(JLON,JLEV)*ZQSAT(JLON,JLEV)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !IF(ZMAXTEND>0.01) WRITE(*,*)'ZMAXTEND',ZMAXTEND
    !IF(ZMINTEND<-0.01) WRITE(*,*)'ZMINTEND',ZMINTEND
    ! write LHN profiles to array for next timestep
    CALL NUDGLHPREP(KIDIA,KFDIA,KLON,KLEV,ZTENDT,KBL)
    IF(MYPROC==1) WRITE(NULOUT,*)'calling LH successful finished'
  ENDIF
ENDIF
! **end latent heat nudging***********

    
!    ------------------------------------------------------------------
!     11 - SAVE FIELDS FOR EXT. SURFACE.
!     --------------------------------------------------------------------
!    Cette partie n'est plus necessaire apres branchement de la physique 
!    de surface sous apl_arome

!    ------------------------------------------------------------------
!     12 - CALL CHEMICAL SCHEME.
!     --------------------------------------------------------------------
IF (LUSECHEM) THEN

  ! ANNEE
  IYEAR = NINDAT / 10000
  ! MOIS
  IMONTH = (NINDAT - 10000*IYEAR ) / 100
  ! JOUR DU MOIS
  IDAY = NINDAT - 10000*IYEAR - 100*IMONTH

  DO JLON = KIDIA,KFDIA
    ZLAT_(JLON) = 180. * ASIN(PGEMU(JLON)) / (2.*ASIN(1.))
    ZLON_(JLON) = 180. * PGELAM(JLON) / (2.*ASIN(1.))
    ZZENITH_(JLON) = ACOS( PMU0(JLON) )
    ZZS_(JLON)=POROG(JLON)/RG
    ZALB_UV_(JLON)=ZALBP(JLON,1)
  ENDDO

  ! Swapp because IN and OUT will be needed simultaneously
  CALL SWAP_SVS

  DO JGFL=1,NGFL_EXT
    DO JLEV=1,KLEV
      DO JLON= KIDIA,KFDIA
        ! modify input
        ZSVSIN_(JLON,JLEV,JGFL)=MAX(0.0_JPRB, ZSVSIN_(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO
  IEZDIAG_CHEM=NGFL_EZDIAG-IOFF_MFSHAL+1
  CALL ARO_MNHC(ZSVSIN_,&
   & ZRHODREFM__(:,1:KLEV),PDT, ZTHM__(:,1:KLEV), ZPABSM__(:,1:KLEV), ZRM_, ZLAT_, ZLON_, ZALB_UV_, &
   & ZZS_, ZZENITH_,ZZZ_, IYEAR,IMONTH,IDAY, REAL(RHGMT,JPRB)+PDT/2._JPRB,&
   & KFDIA,KLEV,NGFL_EXT, NRR, KSTEP+1,NULOUT,IEZDIAG_CHEM, ZPEZDIAG_(:,:,IOFF_MFSHAL:NGFL_EZDIAG),ZSVS_ )
 
  PEZDIAG(KIDIA:KFDIA,1:KLEV,IOFF_MFSHAL:NGFL_EZDIAG)=ZPEZDIAG_(KIDIA:KFDIA,1:KLEV,IOFF_MFSHAL:NGFL_EZDIAG)

  !inversion niveau de la tendance des scalaires passifs
  DO JGFL=1,NGFL_EXT
    DO JLEV = 1,KLEV
      DO JLON = KIDIA,KFDIA
        PTENDEXT(JLON,JLEV,JGFL)=PTENDEXT(JLON,JLEV,JGFL)+(ZSVS_(JLON,JLEV,JGFL)-ZSVSIN_(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO

ENDIF ! LUSECHEM

!    ------------------------------------------------------------------
!     13 - STOCHASTIC PHYSICS : PERTURB TENDENCIES
!     -----------------------------------------------------------------

IF(YSPPT_CONFIG%LSPSDT) THEN

  ZDUMMY(KIDIA:KFDIA,1:KLEV)=0.0_JPRB               ! Dummy nonphys tendency for compatibility with ecmwf stochphy
  CALL SPPTEN (YDMODEL%YRML_PHY_EC%YRECLDP,YGFL, &
   & KIDIA,KFDIA,KLON,KLEV,1,PDT,         &  ! In: block indices, physicstimestep
   & PTSL=PTM,PQSL=PQVM, PA=PCLFS, PAP=PAPRSFM, PAPH=PAPRSM,   &  ! In: (T,Q,cloud) forsupersatcheck, Pfull, Phalf
   & PDYN_U=ZDUMMY,PDYN_V=ZDUMMY,PDYN_T=ZDUMMY,PDYN_Q=ZDUMMY,  &  ! In: dummy nonphys tendencies
   & PUNP_U=PTENDU,PUNP_V=PTENDV,PUNP_T=PTENDT,PUNP_Q=PTENDR(:,:,1), &  ! In: (u,v,t,qv) tendencies to perturb
   & PMULNOISE=PGP2DSDT(1,1,1),                           &  ! In: stochphy 3D random multiplicative pattern (less one)
   & PTENU=PTENDU,PTENV=PTENDV,PTENT=PTENDT,PTENQ=PTENDR(:,:,1) )    ! Out: (u,v,t,qv) total perturbed tendencies
ENDIF

IF(LFORCENL.AND.(KSTEP*(TSPHY/RHOUR)>=NFORCESTART).AND.&
              & (KSTEP*(TSPHY/RHOUR)<=NFORCEEND)) THEN
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      PTENDU(JLON,JLEV)=PTENDU(JLON,JLEV)+AMAGSTOPH_CASBS*PFORCEU(JLON,JLEV)
      PTENDV(JLON,JLEV)=PTENDV(JLON,JLEV)+AMAGSTOPH_CASBS*PFORCEV(JLON,JLEV)
      PTENDT(JLON,JLEV)=PTENDT(JLON,JLEV)+AMAGSTOPH_CASBS*PFORCET(JLON,JLEV)
      PTENDR(JLON,JLEV,1)=PTENDR(JLON,JLEV,1)+AMAGSTOPH_CASBS*PFORCEQ(JLON,JLEV)
    ENDDO
  ENDDO
ENDIF

!    ------------------------------------------------------------------
!     14 - FINAL CALCULATIONS.
!     --------------------------------------------------------------------

!forcage pour declencher la ligne de grain 
IF (LSQUALL) THEN
  IF (LTWOTL) THEN
    ZDT2=2*ZDT
  ELSE
    ZDT2=ZDT
  ENDIF
  IF((KSTEP+1)*ZDT2 < 600._JPRB) THEN
    WRITE(NULOUT, *)'refroidissement impose de',NREFROI1,' a ',NREFROI2
    DO JLEV=KLEV,KLEV-20,-1
      PTENDT(NREFROI1:NREFROI2,JLEV)=-0.01_JPRB
    ENDDO
  ENDIF
ENDIF


!ecriture du buffer
IF(LLMSE.OR.LSFORCS) THEN
  DO JLON = KIDIA,KFDIA
    PGPAR(JLON,MINPRR)=ZINPRR_(JLON)+ZSURFPREP(JLON)/1000._JPRB
    PGPAR(JLON,MINPRS)=ZINPRS_(JLON)+ZSURFSNOW(JLON)/1000._JPRB
    PGPAR(JLON,MINPRG)=ZINPRG_(JLON)+ZINPRH_(JLON)
    PGPAR(JLON,MACPRR)=PGPAR(JLON,MACPRR)+(ZINPRR_(JLON)+ZSURFPREP(JLON)/1000._JPRB)*PDT
    PGPAR(JLON,MACPRS)=PGPAR(JLON,MACPRS)+(ZINPRS_(JLON)+ZSURFSNOW(JLON)/1000._JPRB)*PDT
    PGPAR(JLON,MACPRG)=PGPAR(JLON,MACPRG)+(ZINPRG_(JLON)+ZINPRH_(JLON))*PDT
  ENDDO
  PGPAR(KIDIA:KFDIA,MVTS)=ZTSURF(KIDIA:KFDIA)
  PGPAR(KIDIA:KFDIA,MVEMIS)=ZEMIS(KIDIA:KFDIA)
  PGPAR(KIDIA:KFDIA,MVQS)=ZQS(KIDIA:KFDIA)
  DO JSW=1,NSW
    PGPAR(KIDIA:KFDIA,MALBDIR-1+JSW)=ZALBP(KIDIA:KFDIA,JSW)
    PGPAR(KIDIA:KFDIA,MALBSCA-1+JSW)=ZALBD(KIDIA:KFDIA,JSW)
  ENDDO
ENDIF

IF (LMUSCLFA) CALL ECR1D(NMUSCLFA,'PCLCT_apl',PCLCT,1,KLON)
! initialisations for CFU for Rainfalls
DO JLEV = 0,KLEV
  DO JLON = KIDIA,KFDIA
    ! conversion from m/s in mm/s
    PFPLSL(JLON,JLEV)= ZINPRR_(JLON)*1000._JPRB+ZSURFPREP(JLON)
    PFPLSN(JLON,JLEV)= ZINPRS_(JLON)*1000._JPRB+ZSURFSNOW(JLON)
    PFPLSG(JLON,JLEV)= ZINPRG_(JLON)*1000._JPRB
    PFPLSH(JLON,JLEV)= ZINPRH_(JLON)*1000._JPRB
    ! conversion in correct Unit for BADP (same as ALADIN)
    PSTRTU(JLON,JLEV)= ZSFU_(JLON)*ZRHODREFM__(JLON,IKB) 
    PSTRTV(JLON,JLEV)= ZSFV_(JLON)*ZRHODREFM__(JLON,IKB) 
  ENDDO
ENDDO
!Hail diagnostic
PDIAGH(KIDIA:KFDIA)=0._JPRB
IF (LXXDIAGH) THEN
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      PDIAGH(JLON)=PDIAGH(JLON)+ZQGM(JLON,JLEV)*PDELPM(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO
ENDIF
! lightening density
IF (LFLASH) THEN
  IF (KSTEP==0) PFLASH=0._JPRB

  CALL DIAGFLASH(YDCFU,KIDIA,KFDIA,KLON,KLEV,KSTEP,&
    &ZQCM,ZQIM,ZQRM,ZQSM,ZQGM,ZQHM,&
    &PDELPM,ZTM,PWM,PLSM,PFLASH)
ENDIF
!!! modif pour LMSE non activee
IF (LLMSE) THEN
  DO JLEV=1,KSGST+1
    DO JLON = KIDIA,KFDIA
      PFCS(JLON,JLEV)=-ZSFTH_(JLON)*ZRHODREFM__(JLON,IKB)*RCPD
      PFCLL(JLON,JLEV) = PFCLL(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
      PFCLN(JLON,JLEV) = PFCLN(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
      PFEVL(JLON,JLEV) = PFEVL(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
      PFEVN(JLON,JLEV) = PFEVN(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
    ENDDO
  ENDDO
ENDIF
IF (LSFORCS) THEN
  DO JLEV=1,KSGST+1
    DO JLON = KIDIA,KFDIA
      PFCS(JLON,JLEV)=-ZSFTH_(JLON)*ZRHODREFM__(JLON,IKB)*RCPD
      ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZTSURF(JLON)))
      PFCLL(JLON,JLEV)=-ZSFRV_(JLON)*ZRHODREFM__(JLON,IKB)* FOLH (ZTSURF(JLON),0._JPRB)*(1.0_JPRB-ZDELTA)
      PFCLN(JLON,JLEV)=-ZSFRV_(JLON)*ZRHODREFM__(JLON,IKB)* FOLH (ZTSURF(JLON),0._JPRB)*ZDELTA
    ENDDO
  ENDDO
ENDIF  

DO JSG  = 1, KSGST+1
  DO JLEV = 0, KLEV
    DO JLON = KIDIA, KFDIA
      PFRTH(JLON,JLEV,JSG)=PFRTH(JLON,JLEV,JSG)+ZBUDTH_(JLON)
    ENDDO
  ENDDO
ENDDO

! daand: radflex
IF (LINTFLEX) THEN
  ! account for radiation separately
  LLRAD=.NOT.LRADFLEX
    
  CALL APL_AROME2INTFLEX(YGFL,YDPARAR,YDPHY,KLON,KIDIA,KFDIA,KLEV, PDT,&
   & PRDELPM, PUM, PVM, PTM, PGPAR(1,MVTS), PCPM,&
   & ZFPR,&! precipitation fluxes
   & LLRAD, PFRTH, PFRSO,&! radiative fluxes
   & PTENDU, PTENDV, PTENDT,&! momentum and temperature tendencies
   & PTENDR, PTENDTKE, PTENDEXT,&! total gfl tendencies
   & YDPROCSET)
ENDIF


! Precipitation Type

! Compute wet-bulb temperature at 2 meters (suppose homogeneity of qv/ql/qi )
!ZPCLS(KIDIA:KFDIA)=PAPRSM(KIDIA:KFDIA,KLEV)-2._JPRB/ZZZF(KIDIA:KFDIA,1,KLEV)*&
!                 &(PAPRSM(KIDIA:KFDIA,KLEV)-PAPRSFM(KIDIA:KFDIA,KLEV))

CALL PPWETPOINT(YDPHY,KIDIA,KFDIA,KLON,PAPRSM(:,KLEV),PTCLS,&
  & PQCLS,ZQCM(:,KLEV),ZQIM(:,KLEV),PTPWCLS)

IF (LDPRECIPS.OR.LDPRECIPS2) THEN

  !initialisation de ZDZZ
  DO JLON = KIDIA,KFDIA
    ZDZZ(JLON,1)=ZAPHIM(JLON,0)*ZINVG-ZZZ_(JLON,1)
  ENDDO
  DO JLEV = 2, KLEV
    DO JLON = KIDIA,KFDIA
      ZDZZ(JLON,JLEV)=ZZZ_(JLON,JLEV+IKL)-ZZZ_(JLON,JLEV)
    ENDDO
  ENDDO


  ! Compute wet-bulb temperature
  DO JLEV=1,KLEV
      CALL PPWETPOINT(YDPHY,KIDIA,KFDIA,KLON,PAPRSFM(:,JLEV),ZTM(:,JLEV),&
       & ZQVM(:,JLEV),ZQCM(:,JLEV),ZQIM(:,JLEV),ZTPW(:,JLEV))
  ENDDO

  IF (LDPRECIPS) THEN
   ! Defined precipitation type 
   !
   NDTPRECCUR=INT(MOD(ZSTATI/TSTEP,REAL(NDTPREC)))+1_JPIM
   !PDPRECIPS(:,NDTPRECCUR)=HUGE(1._JPRB)
   PDPRECIPS(:,NDTPRECCUR)=0._JPRB

   !WRITE(NULOUT,*)'sous apl_arome NDTPRECCUR=',NDTPRECCUR,NDTPREC
   CALL DPRECIPS(YDPRECIPS,KIDIA,KFDIA,KLON,KLEV,POROG,PTPWCLS,PDIAGH,PAPHIFM,&
      & ZDZZ,ZTPW,ZQCM,PFPLSL(:,KLEV),PFPLSN(:,KLEV),PFPLSG(:,KLEV),PDPRECIPS(:,NDTPRECCUR))
  ENDIF

  IF (LDPRECIPS2) THEN

   !Idem for an other time step and an other period
   NDTPRECCUR2=INT(MOD(ZSTATI/TSTEP,REAL(NDTPREC2)))+1_JPIM
   PDPRECIPS2(:,NDTPRECCUR2)=0._JPRB

   CALL DPRECIPS(YDPRECIPS,KIDIA,KFDIA,KLON,KLEV,POROG,PTPWCLS,PDIAGH,PAPHIFM,&
      & ZDZZ,ZTPW,ZQCM,PFPLSL(:,KLEV),PFPLSN(:,KLEV),PFPLSG(:,KLEV),PDPRECIPS2(:,NDTPRECCUR2))

  ENDIF

ENDIF



!     --------------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('APL_AROME',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE SWAP_THS
IF (LLSWAP_THS) THEN
  ZTHSIN_ => ZTHSAVE__(:,1:KLEV)
  ZTHS__ => ZTHSWAP__
ELSE
  ZTHSIN_ => ZTHSWAP__(:,1:KLEV)
  ZTHS__ => ZTHSAVE__
ENDIF
LLSWAP_THS=.NOT.LLSWAP_THS
END SUBROUTINE SWAP_THS

SUBROUTINE SWAP_RS
IF (LLSWAP_RS) THEN
  ZRSIN_ => ZRSAVE_
  ZRS_   => ZRSWAP_
ELSE
  ZRSIN_ => ZRSWAP_
  ZRS_   => ZRSAVE_
ENDIF
LLSWAP_RS=.NOT.LLSWAP_RS
END SUBROUTINE SWAP_RS

SUBROUTINE SWAP_SVS
IF (LLSWAP_SVS) THEN
  ZSVSIN_ => ZSVSAVE_
  ZSVS_   => ZSVSWAP_
ELSE
  ZSVSIN_ => ZSVSWAP_
  ZSVS_   => ZSVSAVE_
ENDIF
LLSWAP_SVS=.NOT.LLSWAP_SVS
END SUBROUTINE SWAP_SVS

SUBROUTINE SWAP_SVM
IF (LLSWAP_SVM) THEN
  ZSVMIN_ => ZSVMSAVE_
  ZSVM_   => ZSVMSWAP_
ELSE
  ZSVMIN_ => ZSVMSWAP_
  ZSVM_   => ZSVMSAVE_
ENDIF
LLSWAP_SVM=.NOT.LLSWAP_SVM
END SUBROUTINE SWAP_SVM

SUBROUTINE SWAP_LIMAS
IF (LLSWAP_LIMAS) THEN
  ZLIMASIN_ => ZLIMASAVE_
  ZLIMAS_   => ZLIMASWAP_
ELSE
  ZLIMASIN_ => ZLIMASWAP_
  ZLIMAS_   => ZLIMASAVE_
ENDIF
LLSWAP_LIMAS=.NOT.LLSWAP_LIMAS
END SUBROUTINE SWAP_LIMAS

END SUBROUTINE APL_AROME
