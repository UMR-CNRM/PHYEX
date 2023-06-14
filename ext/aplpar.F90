#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE APLPAR(YDCST, YDMF_PHYS_BASE_STATE, YDMF_PHYS_NEXT_STATE, YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS,  &
& YDCPG_MISC, YDCPG_GPAR, YDCPG_PHY0, YDMF_PHYS, YDCPG_DYN0, YDMF_PHYS_SURF, YDCPG_SL1, YDCPG_SL2, &
& YDVARS, YDGMV, YDSURF, YDCFU, YDXFU, YDMODEL, PGFL, PGMVT1, PGFLT1, PTRAJ_PHYS, &
& YDDDH)

!**** *APLPAR * - APPEL DES PARAMETRISATIONS PHYSIQUES.

!     Sujet.
!     ------
!     - APPEL DES SOUS-PROGRAMMES DE PARAMETRISATION
!       INTERFACE AVEC LES PARAMETRISATIONS PHYSIQUES (IALPP).
!     - CALL THE SUBROUTINES OF THE E.C.M.W.F. PHYSICS PACKAGE.

!**   Interface.
!     ----------
!        *CALL* *APLPAR*

!-----------------------------------------------------------------------

! - 2D (1:KLEV) .

! PGFL       : GFL FIELDS
! PKOZO      : CHAMPS POUR LA PHOTOCHIMIE DE L'OZONE (KVCLIS CHAMPS).
! PKOZO      : FIELDS FOR PHOTOCHEMISTERY OF OZONE   (KVCLIS FIELDS).

! PGPAR       : BUFFER FOR 2D FIELDS - CONTAINS PRECIP, ALBEDO, EMISS, TS
!             : SURFACE FLUXES
! - INPUT/OUTPUT 1D
! YDDDH      : DDH superstructure

!-----------------------------------------------------------------------

!     Externes.
!     ---------

!     Methode.
!     --------
!     - TERMINE LES INITIALISATIONS.
!     - APPELLE LES SS-PRGMS TAMPONS SUIVANT LA LOGIQUE TROUVEE
!        DANS /YOMPHY/. EUX MEMES VONT DECLARER LES TABLEAUX DE TRAVAIL
!        ET APPELER LES PARAMETRISATIONS ELLES MEMES.
!     - FINISH UP THE INITIALIZATION.
!     - CALL THE BUFFER SUBROUTINES FOLLOWING /YOEPHY/ REQUIREMENTS
!        WHICH IN TURN CALL THE ACTUAL PHYSICS SUBROUTINES
!        (THIS LAST POINT NOT PARTIALLY DONE)

!     Auteur.
!     -------
!     90-09-28: A. Joly, *CNRM*.

!     Modifications.
!     --------------
!     2007-02-01 M.Janousek : Introduction of 3MT routines
!     2007-02-19 R.Brozkova : Cleaning obsolet features (LSRCON, LSRCONT, LNEBT,
!                             pre-ISBA, modularisation and racionalisation.
!     2007-04-17 S.Ivatek-S : Over dimensioning of PGPAR (KLON,NGPAR+1) is used
!                             boundary checking bf
!     2007-05-10 E. Bazile  : Introduction of the AROME shallow convection (LCVPPKF)
!     2007-03-21 A. Alias   : Modifications for SURFEX (IGFL_EXT)
!     2007-05-07 F. Bouyssel: Several modifications for SURFEX
!     2007-05-07 F. Bouyssel: New argument in ACCOEFK
!     2007-06-27 A. Alias   : Use NGFL_EXT instead of IGFL_EXT
!     2008-02-18 F. Bouyssel: New acdifv1 & acdifv2 & arp_ground_param
!     2008-02-21 E. Bazile  : Cleaning for the call of the AROME shallow convection (LCVPPKF)
!     4-Mar-2008 Y. Seity : Cleaning IR and WV similated sat radiances
!                            (replaced by Fullpos Calculations)
!     2008-03-14 Y. Bouteloup: Store diffusion coefficients from non-linear model
!     2007-10-11 A. Alias   : New Call to ACHMT/ACNEBR/ACPBLH (P. Marquet/JF. Gueremy)
!     2008-02-01 P. Marquet : modify ZALBD/ZALBP and PFRSODS if LRAYFM15 (idem V4)
!     2008-03-26 F. Bouyssel: Intrduction of LACDIFUS
!     2008-04-28 E. Bazile  : Introduction of ZPROTH_CVPP for the TKE scheme
!     2008-05-09, J.F. Gueremy : Flux MEMO sur mer (ACFLUSO/LFLUSO) +
!            and  P. Marquet   : ZCEROV as new argument in ACDIFUS
!     2008-06-01 F. Bouyssel: Interface of radozc (ECMWF ozone)
!     2008-09-02 F. Vana  : Better split of ACPTKE and ACDIFV1 code
!     2008-10-01 F. Bouyssel: Call of radozcmf instead of radozc (ECMWF ozone)
!     2008-10-05 E. Bazile : Computation of the PBL height from the TKE.
!     03-Oct-2008 J. Masek    parameters for NER statistical model via namelist
!     2009-Jan-21 F. Vana : new mixing lengths for pTKE + few fixes
!     2008-11    C. Payan: Neutral Wind (new arg in the call of ACHMT)
!     2008-11-15 F. Bouyssel: Correction of negative humidities
!     2009-05-01 Y. Bouteloup : Store radiative cloud water and ice in GFL (YIRAD and YLRAD)
!     2009-05-25 F. Bouyssel: Cleaning
!     K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!     2009-08-07 A. Alias   : LCALLSFX introduced to call only once SURFEX at NSTEP=0
!                             Computation of ZRTI (INVERSE DE R*T) added (after ACHMTLS)-not done
!                             Negatives humidity correction PFCQVNG in acdifus (J.F. Gueremy)
!                             add ZAESUL/ZAEVOL to CALL RADAER
!     2009-09-21  D. Banciu: complete the cascade within 3MT frame;
!            prepare the environment for Rash and Kristjansson condensation (RK) scheme
!            remove some arguments of ACPUM and ACUPD (LUDEN option was removed)
!     12-Oct-2009 F. Vana : optimization + update of mixing lengths for p/eTKE
!     2009-10-15 Y. Bouteloup : Store radiative cloud water and ice in GFL (YIRAD and YLRAD)
!     2009-10-23 O.Riviere Intro. of LGWDSPNL for GWD in simpl. phys.
!     2010-01-21 S. Riette: PU, PV and ZDEPTH_HEIGHT in 3D for ARO_GROUND_DIAG
!     2010-05-11 F. Bouyssel : Use of PINDX, PINDY
!     2010-06-20 Y. Seity : Use AROCLDIA to compute PBLH
!     2010-10    A. Alias Compute Sunshine duration
!     2010-10    A. Alias modify ZALBD/ZALBP and PEMIS if LRAYFM for CLIMAT (JF Gueremy)
!     2010-08-10 P.marguinaud : Cleaning
!     2011-01-10 F. Bouyssel: Intro. of LADJCLD and some cleaning.
!     2010-12-01 E. Bazile: TKE1 for AROCLDIA and contributions terms of the
!          TKE equations for DDH.
!     2010-12 S. Riette: aro_ground_diag interface modified to add snow cover
!     2010-12    B. Decharme  : modify the radiative coupling with surfex (SW per band in ACRADIN and RADHEAT)
!     2011-02    A. Alias     : Computation of ZRTI (INVERSE DE R*T) added (after ACHMTLS)
!     2011-02    A. Voldoire : add ZAERINDS to CALL RADAER and ACRADIN
!                              for sulfate indirect effect computation
!     L. Bengtsson-Sedlar & F. Vana 18-Feb-2011 : CA scheme for convection
!     I. Bastak-Duran, F. Vana & R. Brozkova  16-Mar-2011: TOUCANS, version 0
!     2011-02-01 M. Mokhtari: Several modifications for aplpar and introduction of the key LMDUST
!                             (treatment of the desert aerosols)
!     2011-02-24 Y. Bouteloup : EDKF + Surface forcing for MUSC
!     2011-03-26 F. Bouyssel: Intro. of PSPSG (snow cover with surfex)
!     2011-09-07 J.M. Piriou: PCMT convection scheme.
!     2011-11-17 J.F. Gueremy: ZQLI_CVP diagnostic convective water content
!     2011-06: M. Jerczynski - some cleaning to meet norms
!     2011-11-29 K-I. Ivarsson, L. Bengtsson: RK-scheme modifications   
!     26-Jan-2012: F. Vana + I. Bastak-Duran - TOUCANS update + bugfixes 
!     2012-04-24 F. Bouyssel: Bug correction on surface water fluxes with surfex
!     2012-06-09 M. Mile: Bug correction for undefined z0;z0h at 0th step CALL ARO_GROUND_DIAG_Z0
!     2012-09-11 : P.Marguinaud : Add control threshold for
!     2013-06-17 J.M. Piriou: evaporation for PCMT scheme.
!     T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     2013-11-08 Y. Bouteloup New version of ACDIFV1 and ACDIFV2 for "full implicit PMMC09 scheme"
!     F. Vana  28-Nov-2013 : Redesigned trajectory handling.
!     2013-11, J. Masek: Introduction of ACRANEB2 scheme, externalized
!                        computation of direct albedo for ACRANEB/ACRANEB2.
!                        Phasing to cy40t1.
!     K. Yessad (July 2014): Move some variables.
!     2014-09, C. Wastl: Adaptations for orographic shadowing
!     2014-10, R. Brozkova: phasing TOUCANS.
!     2016-03, E. Bazile: phasing MUSC for surf_ideal_flux
!     2016-03, L. Gerard: LNSDO AND LCVCSD
!     2016-04, J. Masek: Exponential-random cloud overlap with variable
!                        decorrelation depth.
!     2016-09, J. Masek: Proper diagnostics of sunshine duration in ACRANEB2.
!     2016-09, M. Mokhtari & A. Ambar: preliminary calculation for passive scalar
!     2016-10, P. Marguinaud : Port to single precision
!     K. Yessad (Dec 2016): Prune obsolete options.
!     2016-06, F.Taillefer: add of aro_ground_diag_2isba call
!     2017-09, Y.Bouteloup: Phased Francoise's modification on cy45
!     2017-09, J. Masek: Fix for protected convective cloudiness,
!                        shifted dimensioning of PGMU0.
!     R. El Khatib 05-Feb-2018 fix bounds violations
!     2018-09, F. Duruisseau: Add PQRCONV1 and PQSCONV1 out (BAYRAD)
!     2018-09, D. St-Martin: Add non-orographic GWD scheme (ACNORGWD)
!     2018-09, R. Roehrig: add ACTKE input/output (ZQLC/ZQIC and ZKQROV/ZKQLROV) (from JF Guérémy)
!     2018-09, M. Michou: Add call to chem_main to activate ARPEGE-Climat chemistry scheme
!     2018-09, R. Brozkova: Fixes in thermodynamic adjustment - deep convective
!                           condensates protection. Passing of diagnostic hail.
!     2018-09, J. Masek: Calculation of snow fractions over bare ground and
!                        vegetation moved to ACSOL (case LVGSN=T). Coding of
!                        ALARO-1 fixes for LZ0HSREL=T. Diagnostics of global
!                        normal irradiance and mean radiant temperature.
!     2018-11, J.M. Piriou: correct 2010 historical bug about adding cloud sedimentation to resolved surface precipitation.
!     R. Hogan     24-Jan-2019 Removed radiation scheme from cycle 15
!     R. El Khatib 30-Apr-2019 fix uninitialized variable
!     2018-10, I. Etchevers : add Visibilities
!     2019-01, I. Etchevers, Y. Seity : add Precipitation Type
!     2019-05, J.M. Piriou: LCVRESDYN + LADJCLD.
!     2019-09, M. Hrastinski: Dataflow for TKE and TTE terms in ALARO DDH (PFTCNS).
!     2019-09, L. Gerard: Modified call to ACNSDO.
!     2019-09, R. Brozkova: Introduction of new NDIFFNEB options.
!     2019-09, J. Masek: Introduction of ETKE_MIN, efficient ACRANEB2 clearsky
!                        computations.
!    2019-10, I. Etchevers : Visibilities in ACVISIH, AROCLDIA=>ACCLDIA
!    2019-10, Y.Bouteloup : New anti-GPS in accvimp.F90 
!    2019-10, Y.Bouteloup and M. Bouzghaiam : Radiation modifications. Remove of FMR15, remove acradin.F90 direct
!                   call to recmwf.F90 and add interface to ecrad (in recmwf !)
!    2020-07, J.M. Piriou and O. Jaron: interface Tiedtke scheme with lightning flash density.
!    2020-10, J. Masek : modified call to ACCLDIA
!    2020-10, M. Hrastinski: Reorganized computation of the moist gustiness
!             correction. Modified call of ACMRIP and ACMIXELEN subroutines.
!    2020-11, Y.Bouteloup : Interface to IFS deep convection scheme under LCVTDK key
!    2020-12, U.Andrae : Introduce SPP for HARMONIE-AROME
!    2021-04, J.M. Piriou: interface PCMT with lightning flash density.
!    2021-09, J.M. Piriou: initialize LLLAND for PCMT lightning computations.
!    2021-    R.Brozkova,D.Nemec: interfaced subroutine DIAGFLASH for ALARO 
!                                 bf for ALARO graupel + removed ZMELNET/ZMELGET
!     13-Jul-2022 R. El Khatib Fix initialization of ZQGM

! End Modifications
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!     ******************************************************************
!     ****** IDIOSYNCRASIES *** IDIOSYNCRASIES *** IDIOSYNCRASIES ******
!     ******************************************************************
!     ***  HEALTH WARNING:                                           ***
!     ***  ===============                                           ***
!     ***  NOTE THAT WITHIN THE E.C.M.W.F. PHYSICS HALF-LEVELS       ***
!     ***  ARE INDEXED FROM 1 TO NFLEVG+1 WHILE THEY ARE BETWEEN     ***
!     ***  0 AND NFLEVG IN THE REST OF THE MODEL. THE CHANGE IS TAKEN***
!     ***  CARE OF IN THE CALL TO THE VARIOUS SUBROUTINES OF THE     ***
!     ***  PHYSICS PACKAGE                                           ***
!     ***                                                            ***
!     ***    THIS IS SUPPOSED TO BE A "TEMPORARY" FEATURE TO BE      ***
!     ***    STRAIGHTENED OUT IN THE "NEAR" FUTURE                   ***
!     ******************************************************************

!     ------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_MISC_TYPE, CPG_DYN_TYPE, CPG_PHY_TYPE, &
                              & CPG_SL1_TYPE, CPG_SL2_TYPE, CPG_GPAR_TYPE
USE CPG_OPTS_TYPE_MOD   , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMXFU             , ONLY : TXFU
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMVERT            , ONLY : VP00
USE YOMCST             , ONLY : TCST
USE YOMRIP0            , ONLY : NINDAT
USE DDH_MIX            , ONLY : TYP_DDH
USE YOMLUN             , ONLY : NULOUT
USE YOMLSFORC          , ONLY : LMUSCLFA,NMUSCLFA
USE YOMTRAJ            , ONLY : TRAJ_PHYS_TYPE, LPRTTRAJ
USE YOMCFU             , ONLY : TCFU !!! for parameters of FLASH
!USE SPP_MOD            , ONLY : YSPP, YSPP_CONFIG
USE MF_PHYS_BASE_STATE_TYPE_MOD &
                     & , ONLY : MF_PHYS_BASE_STATE_TYPE
USE MF_PHYS_NEXT_STATE_TYPE_MOD &
                     & , ONLY : MF_PHYS_NEXT_STATE_TYPE


USE YOMGMV             , ONLY : TGMV
USE SC2PRG_MOD         , ONLY : SC2PRG

USE YOMCT0             , ONLY : LCALLSFX ,LSFORCS, LELAM, LAROME, LCORWAT
USE YOMNUD             , ONLY : NFNUDG   ,LNUDG 
USE YOMSNU             , ONLY : XPNUDG
USE YOMSCM             , ONLY : LGSCM
USE YOMCHET            , ONLY : GCHETN
USE YOMDYNCORE         , ONLY : RPLDARE, RPLRG

USE INTFLEX_MOD        , ONLY : LINTFLEX, TYPE_INTPROCSET, NEWINTPROCSET, CLEANINTPROCSET
!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST),                    INTENT(IN)    :: YDCST
TYPE (MF_PHYS_BASE_STATE_TYPE), INTENT(IN)    :: YDMF_PHYS_BASE_STATE
TYPE (MF_PHYS_NEXT_STATE_TYPE), INTENT(INOUT) :: YDMF_PHYS_NEXT_STATE
TYPE(GEOMETRY),                 INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE),            INTENT(IN)    :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),            INTENT(IN)    :: YDCPG_OPTS
TYPE(CPG_MISC_TYPE),            INTENT(INOUT) :: YDCPG_MISC
TYPE(CPG_GPAR_TYPE),            INTENT(INOUT) :: YDCPG_GPAR
TYPE(CPG_PHY_TYPE),             INTENT(IN)    :: YDCPG_PHY0
TYPE(MF_PHYS_TYPE),             INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),             INTENT(IN)    :: YDCPG_DYN0
TYPE(MF_PHYS_SURF_TYPE),        INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(CPG_SL1_TYPE),             INTENT(INOUT) :: YDCPG_SL1
TYPE(CPG_SL2_TYPE),             INTENT(INOUT) :: YDCPG_SL2
TYPE(FIELD_VARIABLES),          INTENT(INOUT) :: YDVARS
TYPE(TGMV),                     INTENT(IN)    :: YDGMV
TYPE(TSURF),                    INTENT(IN)    :: YDSURF
TYPE(TCFU),                     INTENT(IN)    :: YDCFU
TYPE(TXFU),                     INTENT(IN)    :: YDXFU
TYPE(MODEL),                    INTENT(IN)    :: YDMODEL

 
REAL(KIND=JPRB),                INTENT(INOUT) :: PGFL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB),                INTENT(INOUT) :: PGMVT1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB),                INTENT(INOUT) :: PGFLT1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
 

TYPE (TRAJ_PHYS_TYPE), INTENT(INOUT) :: PTRAJ_PHYS
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH

!     ------------------------------------------------------------------
LOGICAL :: LL_SAVE_PHSURF

INTEGER(KIND=JPIM) :: IFIELDSS

INTEGER(KIND=JPIM) :: INSTEP_DEB,INSTEP_FIN
INTEGER(KIND=JPIM) :: JROF, JSPP

!     --- UPPER AIR PHYSICAL TENDENCIES.
REAL(KIND=JPRB) :: ZTENDH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)     ! Enthalpy tendency.
REAL(KIND=JPRB) :: ZTENDQ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)     ! Moisture tendency.
REAL(KIND=JPRB) :: ZTENDPTKE(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)  ! Pseudo progn. TKE

! GFL tendencies for APL_AROME (assumes YDMODEL%YRML_GCONF%YGFL%NUMFLDS>=YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
! for now, use Jovi's trick :
REAL(KIND=JPRB) :: ZTENDGFL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)   ! GFL tendencies

!     --- UPPER AIR PHYSICAL TENDENCIES FOR AROME.
!       (the previous one are not used in AROME)
REAL(KIND=JPRB) :: ZTENDD (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)        ! d  tendency
REAL(KIND=JPRB) :: ZTENDEXT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)      ! GFL EXTRA tendency
REAL(KIND=JPRB) :: ZTENDEXT_DEP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)  ! GFL EXTRA tendency
REAL(KIND=JPRB) :: ZDIFEXT(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)     ! Extra-GFL fluxes.

REAL(KIND=JPRB) :: ZTENDU (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! U tendency without deep convection contribution
REAL(KIND=JPRB) :: ZTENDV (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! V tendency without deep convection contribution

!     --- RADIATION COEFFICIENTS FOR SIMPLIFIED PHYSICS IN GRID-POINT ---
REAL(KIND=JPRB) :: ZAC(YDCPG_OPTS%KLON,(YDCPG_OPTS%KFLEVG+1)*(YDCPG_OPTS%KFLEVG+1))   ! Curtis matrix.
REAL(KIND=JPRB) :: ZAC_HC(YDCPG_OPTS%KFLEVG+1,YDCPG_OPTS%KFLEVG+1)           ! horizontally-constant field for ZAC.



! required for INTFLEX
TYPE(TYPE_INTPROCSET) :: YLPROCSET

! SPP
REAL(KIND=JPRB) :: ZGP2DSPP(YDCPG_OPTS%KLON,YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL)


REAL(KIND=JPRB), POINTER :: ZPTENDEFB11(:,:), ZPTENDEFB21(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDEFB31(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDG1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDICONV1(:,:), ZPTENDI1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDLCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZP1EZDIAG(:,:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDQ1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDRCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDR1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDSCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDS1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDTKE1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDL1(:,:)

!     ------------------------------------------------------------------
!        ATTENTION SI KVCLIG < 7 LES CHAMPS SUIVANTS NE SONT
!        PAS REELLEMENT ALLOUES EN MEMOIRE.
!*
!     ------------------------------------------------------------------
!     DECLARATION DES TABLEAUX LOCAUX-GLOBAUX DES PARAMETRISATIONS
INTEGER(KIND=JPIM) :: INLAB(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), INLAB_CVPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZNLAB(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZNLABCVP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

!     ------------------------------------------------------------------
!     ARE DIMENSIONNED 0:KLEV ONLY IN ORDER TO KEEP IN MIND
!     THAT THEY ARE COMPUTED AT "HALF LEVELS".
!     THEY ARE USED HOWEVER FROM 1 TO KLEV.

REAL(KIND=JPRB) :: ZXTROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZXUROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZXPTKEROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZRRCOR(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZMRIPP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZMRIFPP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZBNEBCVPP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZBNEBQ(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZMN2PP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZMN2_ES(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZMN2_EQ(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZMN2_DS(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZMN2_DQ(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
! ZMRIMC     : M(C) from P. Marquet's moist Ri computation - for TKE correction after TOMs
! ZMRICTERM  : Rv/R.F(C)-1/M(C).T/Tv from P. Marquet's moist Ri computation - for TKE correction after TOMs
REAL(KIND=JPRB) :: ZMRIMC(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZMRICTERM(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTSTAR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZTSTAR2(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) !diag 
REAL(KIND=JPRB) :: ZTSTARQ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZTSTAR2Q(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)!diag
REAL(KIND=JPRB) :: ZTAU_TKE(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)!DISSIPATION TIME SCALE TAU  -FOR TOM's CALCULATION
REAL(KIND=JPRB) :: ZF_EPS(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)   !  Conversion function lm-L
REAL(KIND=JPRB) :: ZFUN_TTE(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)   !  Function in computation of tte_tilde
REAL(KIND=JPRB) :: ZKTROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZKUROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZNBVNO(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZKQROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZKQLROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZKNROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZFHORM(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZFHORH(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! arrays for 3D turb
! ZFMTKE - F_m function for static K_m computation
! ZFHTKE - F_h function for static K_h computation
REAL(KIND=JPRB) :: ZFMTKE(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZFTTKE(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZRHS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
! ZAUTKE - alpha_u for dry AF scheme
! ZATTKE - alpha_theta for dry AF scheme
REAL(KIND=JPRB) :: ZAUTKE(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZATTKE(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
!ZTH_FUN, ZWW_FUN - T_h, A_h, F_ww - stability functions for TOMs par.
REAL(KIND=JPRB) :: ZTH_FUN(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZWW_FUN(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
! ZFMGST - stability function F_m for moist gustiness correction
REAL(KIND=JPRB) :: ZFMGST(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)


!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZATSLC(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZNEBS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQLIS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZNEBS0(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQLIS0(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZNEBC0(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)   !Nebulosite convective radiative
REAL(KIND=JPRB) :: ZNEBDIFF(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) !Nebulosite: calcul de la diffusion
REAL(KIND=JPRB) :: ZNEBCH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)   !Nebulosite convective condensation
REAL(KIND=JPRB) :: ZUNEBH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)   !Nebulosite convective histo
REAL(KIND=JPRB) :: ZDETFI(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)   !fraction of instantaneous detrained air
REAL(KIND=JPRB) :: ZFPCOR(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZFHP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZLMT(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZZLMT(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZLMU(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZLMU2(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZLMT2(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! temporary storage of lm,lh
REAL(KIND=JPRB) :: ZLML(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! TKE type mixing length
REAL(KIND=JPRB) :: ZLMLTILD(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! 'STATIC' TKE type mixing length
REAL(KIND=JPRB) :: ZOME(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)   ! updraught envt vert vel*dt
REAL(KIND=JPRB) :: ZFALLR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! fall velocity of rain
REAL(KIND=JPRB) :: ZFALLS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! fall velocity of snow
REAL(KIND=JPRB) :: ZFALLG(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! fall velocity of graupel
REAL(KIND=JPRB) :: ZICEFR1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)! Resolved Condensate ice fraction
REAL(KIND=JPRB) :: ZRHCRI(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! Smith scheme critical RH
REAL(KIND=JPRB) :: ZRHDFDA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)! RK scheme change in RH over cloud
REAL(KIND=JPRB) :: ZLHS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)   ! Sublimation latent heat
REAL(KIND=JPRB) :: ZLHV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)   ! Evaporation latent heat
REAL(KIND=JPRB) :: ZLH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! Temporar storage for updated PLH
REAL(KIND=JPRB) :: ZLSCPE(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! Temporar storage for updated PLSCPE
REAL(KIND=JPRB) :: ZQSAT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)  ! Temporar storage for updated PQSAT
REAL(KIND=JPRB) :: ZQSATS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! QSAT of resolved cond./evap. scheme
REAL(KIND=JPRB) :: ZQW(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! Temporar storage for updated PQW
REAL(KIND=JPRB) :: ZRH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! Temporar storage for updated PRH
REAL(KIND=JPRB) :: ZTW(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! Temporar storage for updated PTW)
REAL(KIND=JPRB) :: ZDQ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! Saturation departure for a given thermodynamic state
REAL(KIND=JPRB) :: ZDQM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)   ! maximum saturation departure
REAL(KIND=JPRB) :: ZPOID(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)  ! DP/(RG*DT) FOR A GIVEN LEVEL AND A GIVEN TIME STEP.
REAL(KIND=JPRB) :: ZIPOI(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)  ! INVERSE OF ZPOID.

REAL(KIND=JPRB) :: ZQU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! Updraught Specific moisture
REAL(KIND=JPRB) :: ZTU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! Updraught Temperature
REAL(KIND=JPRB) :: ZUU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! Updraught zonal wind
REAL(KIND=JPRB) :: ZVU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! Updraught merid. wind

REAL(KIND=JPRB) :: ZTMIC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! Temperature for microphysics
REAL(KIND=JPRB) :: ZQMIC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! Specific moisture for microphysics

REAL(KIND=JPRB) :: ZT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)     ! updated temperature T for cascading parameterization
REAL(KIND=JPRB) :: ZTCORR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! temperature corr. for convective cloud
REAL(KIND=JPRB) :: ZU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)     ! updated zonal velocity
REAL(KIND=JPRB) :: ZV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)     ! updated meridional velocity

REAL(KIND=JPRB) :: ZQV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) vapour
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQI(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) cloud ice
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) cloud liquid
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) rain
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) snow
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQG(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) graupel
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) hail
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZCP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! new cp for turbulent diffusion
REAL(KIND=JPRB) :: ZQT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENHA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENQVA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZFCQVNG(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! correction flux increment for neg vapour
REAL(KIND=JPRB) :: ZFCQING(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! correction flux increment for neg ice
REAL(KIND=JPRB) :: ZFCQLNG(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! correction flux increment for neg liquid water
REAL(KIND=JPRB) :: ZFPLSL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! total liquid water flux: diff+sedi+rain
REAL(KIND=JPRB) :: ZFPLSN(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! total solid water flux: diff+sedi+snow
REAL(KIND=JPRB) :: ZFCQL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)   ! condensation flux(liquid)
REAL(KIND=JPRB) :: ZFCQI(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)   ! condensation flux(ice)
REAL(KIND=JPRB) :: ZDIFCQD(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! downdraft flux of specific humidity
REAL(KIND=JPRB) :: ZDIFCQLD(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! downdraft flux of liquid water
REAL(KIND=JPRB) :: ZDIFCQID(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! downdraft flux of  solid water
REAL(KIND=JPRB) :: ZSEDIQL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! sedimentation flux of cloud liquid water
REAL(KIND=JPRB) :: ZSEDIQI(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! sedimentation flux of cloud ice water
REAL(KIND=JPRB) :: ZDIFCSD(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! downdraft entalphy flux
REAL(KIND=JPRB) :: ZSTRCUD(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! change in horizontal mom.
REAL(KIND=JPRB) :: ZSTRCVD(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! change in horizontal mom.
REAL(KIND=JPRB) :: ZRCVOTT(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! degree of inhomogeneity in precips.
REAL(KIND=JPRB) :: ZSIGPC(YDCPG_OPTS%KLON)         ! Convective precipit mesh fraction
REAL(KIND=JPRB) :: ZSIGP(YDCPG_OPTS%KLON)         ! Precipitation mesh fraction
REAL(KIND=JPRB) :: ZAUXPRC(YDCPG_OPTS%KLON)        ! Precipitation auxilary
REAL(KIND=JPRB) :: ZDIFCVPPQ(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de CVPP (KFB or EDKF) sur Qv
REAL(KIND=JPRB) :: ZDIFCVPPS(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de CVPP (KFB or EDKF) sur CpT
REAL(KIND=JPRB) :: ZDIFCVTH(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de CV sur Theta air sec
REAL(KIND=JPRB) :: ZDIFCVPPU(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de CVPP (EDKF) sur U
REAL(KIND=JPRB) :: ZDIFCVPPV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de CVPP (EDKF) sur V

REAL(KIND=JPRB) :: ZEDMFQ(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Mass flux part of EDMF flux for Qv
REAL(KIND=JPRB) :: ZEDMFS(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Mass flux part of EDMF flux for CpT
REAL(KIND=JPRB) :: ZEDMFU(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Mass flux part of EDMF flux for U
REAL(KIND=JPRB) :: ZEDMFV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Mass flux part of EDMF flux for V
REAL(KIND=JPRB) :: ZMF_UP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Mass flux for implicit formulation of EDMF equation (LEDMFI)
REAL(KIND=JPRB) :: ZMU(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de masse (updraft) pour XIOS output
REAL(KIND=JPRB) :: ZMD(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de masse (downdraft) pour XIOS output

REAL(KIND=JPRB) :: ZCONDCVPPL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de condensation liquide du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZCONDCVPPI(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de condensation glace du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZPRODTH_CVPP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de production thermique de TKE du a CVPP(KFB)
REAL(KIND=JPRB) :: ZDTRAD(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! radiation contribution to T tendency
REAL(KIND=JPRB) :: ZDQVDIFF(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! turtb.diff contribution to Qv tendency
REAL(KIND=JPRB) :: ZRKQCTEND(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! Qc input for RK condensation scheme
REAL(KIND=JPRB) :: ZRKQVTEND(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! Qv input for RK condensation scheme
REAL(KIND=JPRB) :: ZRKTTEND(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! T input for RK condensation scheme
REAL(KIND=JPRB) :: ZDQV, ZDQI, ZDQL, ZDQR, ZDQS, ZDQC, ZGDT, ZGDTI,&
                 & ZQV0, ZQX0, ZQX1,&
                 & ZCONVC, ZTOTC,ZDTURDIFF
REAL(KIND=JPRB) :: ZTMPPRODTH(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! temporary array

!!for BAYRAD
REAL(KIND=JPRB) :: ZDE2MR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! temporary array for conversion of density to mixing ratio
!-----------------------------------------------------------------

! - 2D (0:KLEV) .

! ZKTROV     : COEFFICIENT D'ECHANGE VERTICAL DE T ET Q EN KG/(M*M*S).
! ZKUROV     : COEFFICIENT D'ECHANGE VERTICAL DE U ET V EN KG/(M*M*S).
! ZKNROV     : COEFFICIENT D'ECHANGE VERTICAL NEUTRE EN KG/(M*M*S).
! ZNBVNO     : FREQUENCE DE BRUNT-VAISALA DIVISEE PAR G FOIS LA DENSITE.

! ZNEBS      : NEBULOSITE STRATIFORME (SCHEMA STATISTIQUE DE NUAGES).
!            : STRATIFORM CLOUDINESS (STATISTICAL CLOUD SCHEME)
! ZQLIS      : QUANTITE D'EAU LIQUIDE STRATIFORME (SCHEMA STATISTIQUE).
!            : STRATIFORM LIQUID WATER (STATISTICAL CLOUD SCHEME)


! - 2D (1:KLEV) .

! INLAB      : INDICE D'INSTABILITE CONVECTIVE.

INTEGER(KIND=JPIM) :: INND(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZXDROV(YDCPG_OPTS%KLON),ZXHROV(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZUGST(YDCPG_OPTS%KLON),ZVGST(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZCDROV(YDCPG_OPTS%KLON),ZCHROV(YDCPG_OPTS%KLON),ZDQSTS(YDCPG_OPTS%KLON),ZGWDCS(YDCPG_OPTS%KLON),&
 & ZHQ(YDCPG_OPTS%KLON),ZHU(YDCPG_OPTS%KLON),ZHTR(YDCPG_OPTS%KLON),ZCDNH(YDCPG_OPTS%KLON),ZMOD(YDCPG_OPTS%KLON),&
 & ZRTI(YDCPG_OPTS%KLON),ZDPHI(YDCPG_OPTS%KLON),ZPRS(YDCPG_OPTS%KLON),ZSTAB(YDCPG_OPTS%KLON),ZTAUX(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZWFC(YDCPG_OPTS%KLON),ZWPMX(YDCPG_OPTS%KLON),ZWLMX(YDCPG_OPTS%KLON),ZWSEQ(YDCPG_OPTS%KLON),&
 & ZWSMX(YDCPG_OPTS%KLON),ZWWILT(YDCPG_OPTS%KLON),&
 & ZC3(YDCPG_OPTS%KLON),ZCG(YDCPG_OPTS%KLON),ZCN(YDCPG_OPTS%KLON),&
 & ZNEIJG(YDCPG_OPTS%KLON),ZNEIJV(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZPCLS(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZPREN(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZFRSODS(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZCD(YDCPG_OPTS%KLON)


! - 1D (DIAGNOSTIQUE) .

! ZCDROV     : PCD RENORME EN DENSITE FOIS VITESSE.
! ZCDNH      : COEFFICIENT NEUTRE D'ECHANGE EN SURFACE POUR LA CHALEUR.
! ZCDNH      : EXCHANGE COEFF. AT SURFACE LEVEL IN NEUTRAL CONDITIONS FOR HEAT.
! ZCG        : COEFFICIENT THERMIQUE DU SOL NU.
! ZCG        : THERMICAL COEFFICIENT OF BARE GROUND.
! ZCN        : COEFFICIENT THERMIQUE DE LA NEIGE.
! ZCN        : THERMICAL COEFFICIENT OF SNOW.
! ZCHROV     : PCH RENORME EN DENSITE FOIS VITESSE.
! ZC3        : COEFFICIENT UTILE POUR LE CALCUL DU DRAINAGE
! ZDQSTS     : DERIVEE DE PQSATS PAR RAPPORT A LA TEMPERATURE.
! ZGWDCS     : VARIABLE DE SURFACE POUR LE DRAG OROGRAPHIQUE (RHO*N0/G).
! ZHQ        : POIDS DE L'HUMIDITE DE L'AIR DANS L'HUMIDITE DE SURFACE.
! ZHTR       : RESISTANCE A LA TRANSPIRATION DU COUVERT VEGETAL.
! ZHTR       : FOLIAGE TRANSPIRATION RESISTANCE.
! ZHU        : POIDS DE L'HUMIDITE SATURANTE DANS L'HUMIDITE DE SURFACE.
! ZNEIJG     : FRACTION DE NEIGE RECOUVRANT LE SOL.
! ZNEIJV     : FRACTION DE NEIGE RECOUVRANT LA VEGETATION.
! ZRTI       : INVERSE DE R*T.
! ZDPHI      : EPAISSEUR EN GEOPOTENTIEL DU NIVEAU DE SURFACE.
! ZPRS       : CONSTANTE DES GAZ POUR L'AIR AU SOL.
! ZSTAB      : INDICE DE STABILITE A LA SURFACE.
! INND       : INDICE DE PRECIPITATIONS CONVECTIVES.
! ZWFC       : TENEUR EN EAU A LA CAPACITE AUX CHAMPS.
! ZWFC       : FIELD CAPACITY WATER CONTENT.
! ZWPMX      : TENEUR EN EAU MAXIMALE DU RESERVOIR PROFOND.
! ZWPMX      : MAXIMUM WATER CONTENT OF THE DEEP WATER-TANK.
! ZWLMX      : TENEUR EN EAU MAXIMALE DU RESERVOIR D'INTERCEPTION.
! ZWLMX      : MAXIMUM WATER CONTENT OF THE INTERCEPTION WATER-TANK.
! ZWSEQ      : TENEUR EN EAU A L'EQUILIBRE (EQUILIBRE ENTRE GRAVITE ET
!              CAPILLARITE) EN SURFACE.
! ZWSEQ      : SURFACE WATER CONTENT FOR THE BALANCE BETWEEN GRAVITY
!              AND CAPILLARITY
! ZWSMX      : TENEUR EN EAU MAXIMALE DU RESERVOIR SUPERFICIEL.
! ZWSMX      : MAXIMUM WATER CONTENT FOR THE SUPERFICIAL WATER-TANK.
! ZWWILT     : TENEUR EN EAU AU POINT DE FLETRISSEMENT.
! ZWWILT     : WATER CONTENT AT THE WILTING POINT.
! ZSSO_STDEV : OROGRAPHY STANDARD DEVIATION
! ZTWSNOW    : SNOW COVER FROM SURFEX
! ZTOWNS     : FRACTION OF TOWN FROM SURFEX
REAL(KIND=JPRB) :: ZDAER(YDCPG_OPTS%KFLEVG), ZBLH(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZQO3(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZAER(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,6)
REAL(KIND=JPRB) :: ZAERO(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,12)
REAL(KIND=JPRB) :: ZAERD(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZAERINDS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQCO2(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZROZ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZDPHIV(YDCPG_OPTS%KLON),ZDPHIT(YDCPG_OPTS%KLON)

REAL(KIND=JPRB) :: ZMAN(0:YDCPG_OPTS%KFLEVG), ZMAK(0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZSSO_STDEV(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZTWSNOW(YDCPG_OPTS%KLON),ZTOWNS(YDCPG_OPTS%KLON)
! - (PROFILS DIAGNOSTIQUES)
! ZDAER     : EPAISSEUR OPTIQUE DES AEROSOLS STANDARDS DANS LA COUCHE
! ZDAER     : OPTICAL DEPTH OF STANDARD AEROSOLS IN THE LAYER.
!   ZQO3    : RAPPORT DE MELANGE MASSIQUE D'OZONE
!        (0):     "         "    MOYEN AU-DESSUS DU MODELE
!   ZQO3    : OZONE MIXING RATIO (MASS).
!        (0):AVERAGED-ABOVE      "         " .
!   ZQCO2   : RAPPORT MASSIQUE LOCAL DU CO2.
!   ZQCO2   : CO2 MIXING RATIO (MASS).

! IJN        : DIMENSION TABLEAUX ETENDUS POUR CYCLE DIURNE RAYONNEMENT
!              IJN AU PLUSL A KLON

!* INPUT ARGUMENTS FOR ACRADIN ( RAYT ECMWF POUR CLIMAT )

!            1-D ARRAYS
!            ----------

REAL(KIND=JPRB) :: ZTRSOD(YDCPG_OPTS%KLON)

!            2-D ARRAYS
!            ----------

!* OUTPUT ARGUMENTS FOR THE ECMWF PHYSICS

!            0.2  LOCAL ARRAYS FOR ECMWF PHYSICS PACKAGE
!                 --------------------------------------

REAL(KIND=JPRB) :: ZCEMTR(YDCPG_OPTS%KLON,0:1) , ZCTRSO(YDCPG_OPTS%KLON,0:1)
REAL(KIND=JPRB) :: ZALBD(YDCPG_OPTS%KLON,YDCPG_OPTS%KSW), ZALBP(YDCPG_OPTS%KLON,YDCPG_OPTS%KSW),&
                 & ZALB(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZSFSWDIR (YDCPG_OPTS%KLON,YDCPG_OPTS%KSW), ZSFSWDIF (YDCPG_OPTS%KLON,YDCPG_OPTS%KSW)
REAL(KIND=JPRB) :: ZTRSODIR (YDCPG_OPTS%KLON,YDCPG_OPTS%KSW), ZTRSODIF (YDCPG_OPTS%KLON,YDCPG_OPTS%KSW)
REAL(KIND=JPRB) :: ZFSDNN(YDCPG_OPTS%KLON),ZFSDNV(YDCPG_OPTS%KLON)

!            1-D ARRAYS
!            ----------

REAL(KIND=JPRB) :: ZSUDU(YDCPG_OPTS%KLON) , ZDSRP(YDCPG_OPTS%KLON) , ZSDUR(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZTHETAVS(YDCPG_OPTS%KLON), ZTHETAS(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZAESEA(YDCPG_OPTS%KLON), ZAELAN(YDCPG_OPTS%KLON), ZAESOO(YDCPG_OPTS%KLON), ZAEDES(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZAESUL(YDCPG_OPTS%KLON), ZAEVOL(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZMERL(YDCPG_OPTS%KLON)

!            2-D ARRAYS
!            ----------

REAL(KIND=JPRB) :: ZTENT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZGEOSLC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTHETAV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZDUM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

!            LOCAL ARRAYS FOR TKE
! ZCOEFN : COEFFICIENT STATISTIQUE POUR LES FLUX D'EAUX CONDENSEES.

REAL(KIND=JPRB) :: ZCOEFN(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

!            LOCAL ARRAYS FOR ACVPPKF
REAL(KIND=JPRB) :: ZQLI_CVPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZNEB_CVPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
!            LOCAL ARRAYS FOR EDKF
REAL(KIND=JPRB) :: ZIMPL

!           2-D ARRAY FOR SIMPL.RADIATION SCHEME

REAL(KIND=JPRB) :: ZZNEB(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

INTEGER(KIND=JPIM) :: IJN, IMLTYPE, JCHA, JLEV, JLON, JSG, JGFL, JDIAG, IAERO
INTEGER(KIND=JPIM) :: ILONMNH, IKRR, ISPLITR ! useful size of klon arrays for mesonh physics
LOGICAL :: LLCLS, LLHMT, LLMAF, LLREDPR

REAL(KIND=JPRB) :: ZAEN, ZAEO, ZALBV, ZCARDI, ZEPS0, ZEPSNEB, ZEPSO3
REAL(KIND=JPRB) :: ZALBPMER

!            2-D ARRAYS

! ZNEBC      : NEBULOSITE  CONVECTIVE A L'ECHELLE DE LA MAILLE.
! ZQLIC      : EAU LIQUIDE CONVECTIVE A L'ECHELLE DE LA MAILLE.
! ZQCL       : CONDENSAT STRATIFORME LIQUIDE
! ZQCI       : CONDENSAT STRATIFORME SOLIDE
! ZFHEVPPC   : FLUX DE CHALEUR DU A L'EVAPORATION DES PREC. CONVECTIVES.
! ZFHMLTSC   : FLUX DE CHALEUR DU A LA FONTE/GEL DES PREC. CONVECTIVES.
! ZFPEVPPC   : EVAPORATION DES PREC. CONVECTIVES.
! ICIS       : INDICE DE NIVEAU D'INSTABILITE SECHE.

REAL(KIND=JPRB) :: ZNEBC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQLIC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQCL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQCI(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZFHEVPPC(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZFHMLTSC(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)&
 & ,ZFPEVPPC(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
INTEGER(KIND=JPIM) :: ICIS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)


!           SURFEX local VARIABLES
!           ----------------------------

! Implicit coupling coefficients
INTEGER(KIND=JPIM) :: IRR
REAL(KIND=JPRB) :: ZDTMSE,ZRHGMT,ZSTATI,ZADTMS
REAL(KIND=JPRB) :: ZCFAQ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZCFAS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),&
 & ZCFATH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZCFAU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZCFBQ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),&
 & ZCFBS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZCFBTH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),&
 & ZCFBU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZCFBV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),&
 & ZCFBU_G(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZCFBV_G(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),&
 & ZCFBS_G(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZCFBQ_G(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZDSE(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZFEV(YDCPG_OPTS%KLON),ZFMDU(YDCPG_OPTS%KLON),ZFMDV(YDCPG_OPTS%KLON),ZFEVS(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZSRAIN(YDCPG_OPTS%KLON), ZSSNOW(YDCPG_OPTS%KLON), ZSGROUPEL(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZSFCO2(YDCPG_OPTS%KLON), ZRHODREFM(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZDEPTH_HEIGHT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZZS(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZTSN(YDCPG_OPTS%KLON),ZTN(YDCPG_OPTS%KLON)
REAL(KIND=JPRB)   :: ZBUDTH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)   :: ZBUDSO (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)   :: ZFCLL  (YDCPG_OPTS%KLON)
! FOR Hv
REAL(KIND=JPRB)   :: ZHV2(YDCPG_OPTS%KLON)
! FOR DUST
REAL(KIND=JPRB), DIMENSION (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ::  ZQDM
REAL(KIND=JPRB) :: ZCFASV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZCFBSV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZSMOOTRAC(1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZINVDT, ZINVG, ZRSCP, ZINVATM

REAL(KIND=JPRB),  DIMENSION(YDCPG_OPTS%KLON,1,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT):: ZZI_SVM
REAL(KIND=JPRB),  DIMENSION (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG):: ZZI_PEZDIAG
REAL(KIND=JPRB), DIMENSION (:,:,:), ALLOCATABLE  :: ZSVM, ZPSV
REAL(KIND=JPRB), DIMENSION (:,:),   ALLOCATABLE  :: ZSFSV  ! passifs scalaires surf flux
! TRAITEMENT DES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZTM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), DIMENSION (YDCPG_OPTS%KLON,1,YDCPG_OPTS%KFLEVG) :: ZZZ,ZDZZ,ZZI_PABSM, ZZI_THM,&
                & ZZI_EXNREFM, ZZI_RHODREFM,ZEVAP,ZZDEP,ZZI_RHO
REAL(KIND=JPRB) :: ZZI_APHI(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), DIMENSION (YDCPG_OPTS%KLON,1,YDCPG_OPTS%KFLEVG,6) :: ZZI_RM
!            3-D ARRAYS
REAL(KIND=JPRB), DIMENSION(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KSW):: ZPIZA_DST !Single scattering
                                             ! albedo of dust (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KSW):: ZCGA_DST  !Assymetry factor
                                             ! for dust (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KSW):: ZTAUREL_DST !tau/tau_{550}
                                             !dust (points,lev,wvl)

!         ACFLUSO (ECUME) local variable
!-------------------------------------------
REAL(KIND=JPRB) :: ZCE(YDCPG_OPTS%KLON), ZCEROV(YDCPG_OPTS%KLON), ZCRTI(YDCPG_OPTS%KLON)


!        New ACDIFV1 local variable
!--------------------------------------------
REAL(KIND=JPRB)   :: ZXURO(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZXQRO(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZXTRO(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)

!        New ACNORGWD local variables
!--------------------------------------------
REAL(KIND=JPRB) :: ZD_U(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZD_V(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: Z_PP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) 
REAL(KIND=JPRB) :: Z_UU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) 
REAL(KIND=JPRB) :: Z_VV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: Z_TT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: Z_VO(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)  
REAL(KIND=JPRB) :: ZFLX_LOTT_GWU(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG), ZFLX_LOTT_GWV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZPRECGWD(YDCPG_OPTS%KLON)

!    TKE+ for ACCLDIA
REAL(KIND=JPRB)   :: ZTKE1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZTPRDY(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

!   For ACVISIH
REAL(KIND=JPRB)   :: ZQGM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)


!        New ARP_GROUND_PARAM local variable
!------------------------------------------------

REAL(KIND=JPRB)   :: ZALPHA1(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZCOEFA (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZLVT   (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZQICE  (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZDIFWQ (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)   :: ZDIFWS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)   :: ZSC_FEVI (YDCPG_OPTS%KLON),ZSC_FEVN(YDCPG_OPTS%KLON),ZSC_FCLL(YDCPG_OPTS%KLON),ZSC_FCLN(YDCPG_OPTS%KLON)

!           TRAJECTORY (For diffusion !) local VARIABLES
!           ----------------------------
REAL(KIND=JPRB) :: ZCDROV_SAVE(YDCPG_OPTS%KLON),ZCHROV_SAVE(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZKTROV_SAVE(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZKUROV_SAVE(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTRAJGWD(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) !Traj buffer saved for TL/AD (YDMODEL%YRML_PHY_MF%YRSIMPHL%LGWDSPNL)

REAL(KIND=JPRB)    :: ZRVMD,ZDELTA
LOGICAL            :: LLAERO, LLAROME, LLCALLRAD
REAL(KIND=JPRB)    :: ZAIPCMT(YDCPG_OPTS%KLON) ! Activity Index of PCMT: 1. if PCMT is active, 0. else case.
REAL(KIND=JPRB)    :: ZALF_CAPE(YDCPG_OPTS%KLON)
REAL(KIND=JPRB)    :: ZALF_CVGQ(YDCPG_OPTS%KLON)
REAL(KIND=JPRB)    :: ZQIC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQLC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQRC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQSC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQVI
REAL(KIND=JPRB)    :: ZQLI_CVP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)    :: ZQC_DET_PCMT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)    :: ZFPLS(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZFPLC(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZFPL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)    :: ZCSGC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)    :: ZTZER
CHARACTER(LEN=200) :: CLERR

! Tracers: prognostique aerosols, passive scalars...
INTEGER(KIND=JPIM) :: INBTRA
INTEGER(KIND=JPIM) :: INBTRA_DEP
REAL(KIND=JPRB), ALLOCATABLE :: ZSTRCTRA(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZTRA(:,:,:)

!        New chemistry local variables
!-------------------------------------
INTEGER(KIND=JPIM)               :: IFLDX, IFLDX2, ILEVX 
INTEGER(KIND=JPIM), ALLOCATABLE  :: INDCHEM(:), IGPLAT(:)
REAL(KIND=JPRB), ALLOCATABLE     :: ZSD_XA(:,:,:), ZSD_X2(:,:), ZTENGFL(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE     :: ZCFLX(:,:), ZCFLXO(:,:), ZCHEMDV(:,:), ZAEROP(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE     :: ZTENC(:,:,:), ZDELP(:,:), ZWND(:), ZDUMMY1(:,:), ZGELAT(:)
REAL(KIND=JPRB), ALLOCATABLE     :: ZNEEFLX(:), ZCHEM2AER(:,:,:)
REAL(KIND=JPRB) :: ZDCAPE(YDCPG_OPTS%KLON) ! Descending CAPE for gusts.

! ACRANEB/ACRANEB2 local variables
! --------------------------------
REAL(KIND=JPRB) :: ZLAMB           ! proportion of Lambertian scattering
REAL(KIND=JPRB) :: ZALBDIR  (YDCPG_OPTS%KLON) ! direct (parallel) surface albedo
REAL(KIND=JPRB) :: ZCLCT_RAD(YDCPG_OPTS%KLON) ! total cloud cover for radiation
REAL(KIND=JPRB) :: ZDECRD   (YDCPG_OPTS%KLON) ! decorrelation depth for cloud overlaps
REAL(KIND=JPRB) :: ZDECRD_MF(YDCPG_OPTS%KLON) ! decorrelation depth for cloud overlaps
                                   ! in microphysics

REAL(KIND=JPRB) :: ZDQG

! IFS deep convection scheme local variables
! --------------------------------
INTEGER(KIND=JPIM) :: ITOPC(YDCPG_OPTS%KLON),IBASC(YDCPG_OPTS%KLON),ITYPE(YDCPG_OPTS%KLON),ISPPN2D
INTEGER(KIND=JPIM) :: ICBOT(YDCPG_OPTS%KLON),ICTOP(YDCPG_OPTS%KLON),IBOTSC(YDCPG_OPTS%KLON)
INTEGER(KIND=JPIM) :: ICBOT_LIG(YDCPG_OPTS%KLON),ICTOP_LIG(YDCPG_OPTS%KLON)
LOGICAL :: LLDSLPHY,LLPTQ,LLLAND(YDCPG_OPTS%KLON),LLCUM(YDCPG_OPTS%KLON),LLSC(YDCPG_OPTS%KLON),LLSHCV(YDCPG_OPTS%KLON),LLLINOX(YDCPG_OPTS%KLON)
LOGICAL :: LLCUM_LIG(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZLIGH_CTG(YDCPG_OPTS%KLON),ZCTOPH(YDCPG_OPTS%KLON),ZPRECMX(YDCPG_OPTS%KLON),ZICE(YDCPG_OPTS%KLON),ZCDEPTH(YDCPG_OPTS%KLON),ZWMFU(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZVERVEL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZGEOM1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZGEOMH(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENQ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZTENU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZTENV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENTA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZTENQA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZLCRIT_AER(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZSNDE(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,2)
REAL(KIND=JPRB) :: ZCUCONVCA(YDCPG_OPTS%KLON),ZGAW(YDCPG_OPTS%KLON),ZVDIFTS,ZDXTDK(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZLU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZLUDE(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZMFU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZLISUM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZMFD(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZWMEAN(YDCPG_OPTS%KLON),ZACPR(YDCPG_OPTS%KLON),ZDIFF(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZVDISCU(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZMFUDE_RATE(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZMFDDE_RATE(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZFHPCL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZFHPCN(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZFCQLF(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZFCQLI(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZFRSO(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZFRTH(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZGP2DSPPA(YDCPG_OPTS%KLON,1),ZLUDELI(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,4),ZLRAIN(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZRSUD(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,2)
REAL(KIND=JPRB), ALLOCATABLE :: ZCEN(:,:,:),ZSCAV(:)

! Precipitation type diagnostics
!--------------------------------
REAL(KIND=JPRB)   :: ZFPLSG(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZENTCH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

INTEGER (KIND=JPIM)  :: IMOC_CLPH (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZADJ_DTAJU (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZADJ_TAUX (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZBAY_QRCONV (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZBAY_QSCONV (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZDSA_C1 (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZDSA_C2 (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZDSA_CPS (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZDSA_LHS (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZDSA_RS (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZFLU_CDN (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZFLU_CD (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZFLU_CH (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZFLU_EMIS (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZFLU_FEVI (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KTSSG+1)
REAL (KIND=JPRB)     :: ZFLU_NEIJ (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZFLU_QS1 (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZFLU_QSATS (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZFLU_QSAT (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZFLU_VEG (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZKUR_KTROV_H (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZKUR_KUROV_H (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZMSC_FHP (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZMSC_FRMQ (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZMSC_LH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZMSC_LSCPE (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZMSC_QW (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZMSC_TW (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FEFB1 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FEFB2 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FEFB3 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FPLCH (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FPLSH (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FPLSN (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FP (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FTKEI (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FTKE (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_COR (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB3C (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB3N (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB4C (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB4N (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB6C (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB6N (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT1C (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT1N (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT2C (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT2N (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT3C (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT3N (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT4C (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT4N (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT5C (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT5N (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZTDS_TDALBNS (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDRHONS (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDSNS (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDTP (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%YRSURF_DIMS%YSP_SBD%NLEVS)
REAL (KIND=JPRB)     :: ZTDS_TDTS (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDWL (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDWPI (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDWP (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDWSI (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDWS (YDCPG_OPTS%KLON)

REAL (KIND=JPRB)     :: ZPRC_DPRECIPS2 (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%NDTPREC2)
REAL (KIND=JPRB)     :: ZPRC_DPRECIPS (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%NDTPREC)
REAL (KIND=JPRB)     :: ZRDG_CVGQ (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDG_LCVQ (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZRDG_MMU0 (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZRDG_MU0LU (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZRDG_MU0M (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZRDG_MU0N (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZRDG_MU0 (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZSAV_DDAL (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZSAV_DDOM (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZSAV_ENTCH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZSAV_FHPS (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZSAV_GZ0F (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZSAV_GZ0HF (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZSAV_HV (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZSAV_PBLH (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZSAV_QSH  (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZSAV_UDAL (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZSAV_UDGRO (YDCPG_OPTS%KLON)
REAL (KIND=JPRB)     :: ZSAV_UDOM (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB)     :: ZSAV_UNEBH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_1

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "acaa1.intfb.h"
#include "acajucv.intfb.h"
#include "accdev.intfb.h"
#include "accldia.intfb.h"
#include "acclph.intfb.h"
#include "accoefk.intfb.h"
#include "accsu.intfb.h"
#include "accvimpgy.intfb.h"
#include "accvimp.intfb.h"
#include "accvimp_v3.intfb.h"
#include "accvud.intfb.h"
#include "acdayd.intfb.h"
#include "acdifoz.intfb.h"
#include "acdifus.intfb.h"
#include "acdifv1.intfb.h"
#include "acdifv2.intfb.h"
#include "acdifv3.intfb.h"
#include "acdnshf.intfb.h"
#include "acdrac.intfb.h"
#include "acdrag.intfb.h"
#include "acdrme.intfb.h"
#include "acdrov.intfb.h"
#include "acevadcape.intfb.h"
#include "acfluso.intfb.h"
#include "achmt.intfb.h"
#include "achmtls.intfb.h"
#include "acmixelen.intfb.h"
#include "acmixlentm.intfb.h"
#include "acmixlenz.intfb.h"
#include "acmodo.intfb.h"
#include "acmrip.intfb.h"
#include "acmris.intfb.h"
#include "acmriss.intfb.h"
#include "acnebc.intfb.h"
#include "acnebcond.intfb.h"
#include "acnebn.intfb.h"
#include "acnebr.intfb.h"
#include "acnorgwd.intfb.h"
#include "acnpart.intfb.h"
#include "acnsdo.intfb.h"
#include "acozone.intfb.h"
#include "acpblh.intfb.h"
#include "acpblhtm.intfb.h"
#include "acpcmt.intfb.h"
#include "acpluie.intfb.h"
#include "acpluis.intfb.h"
#include "acpluiz.intfb.h"
#include "acptke.intfb.h"
#include "acradcoef.intfb.h"
#include "acraneb2.intfb.h"
#include "acraneb.intfb.h"
#include "acrso.intfb.h"
#include "acsol.intfb.h"
#include "actkecoefkh.intfb.h"
#include "actkecoefk.intfb.h"
#include "actkehmt.intfb.h"
#include "actke.intfb.h"
#include "actkezotls.intfb.h"
#include "actqsat.intfb.h"
#include "actqsats.intfb.h"
#include "acupd.intfb.h"
#include "acupm.intfb.h"
#include "acuptq.intfb.h"
#include "acupu.intfb.h"
#include "acveg.intfb.h"
#include "acvisih.intfb.h"
#include "acvppkf.intfb.h"
#include "aplmphys.intfb.h"
#include "aplpar2intflex.intfb.h"
#include "aplpar_init.intfb.h"
#include "aro_ground_diag_2isba.h"
#include "aro_ground_diag.h"
#include "aro_ground_diag_z0.h"
#include "aro_ground_param.h"
#include "aro_mnhdust.h"
#include "arp_ground_param.intfb.h"
#include "checkmv.intfb.h"
!include "chem_main.intfb.h"
#include "cpchet.intfb.h"
#include "cpmvvps.intfb.h"
#include "cpnudg.intfb.h"
#include "cpozo.intfb.h"
#include "cpphinp.intfb.h"
#include "cpqsol.intfb.h"
#include "cptend_flex.intfb.h"
#include "cptend_new.intfb.h"
#include "cptends.intfb.h"
#include "cputqy_aplpar_expl.intfb.h"
#include "cputqy_aplpar_loop.intfb.h"
#include "cpwts.intfb.h"
#include "cucalln_mf.intfb.h"
!#include "culight.intfb.h"
#include "dprecips.intfb.h"
#include "mean_rad_temp.intfb.h"
#include "mf_phys_bayrad.intfb.h"
#include "mf_phys_corwat.intfb.h"
#include "mf_phys_cvv.intfb.h"
#include "mf_phys_fpl_part1.intfb.h"
#include "mf_phys_fpl_part2.intfb.h"
#include "mf_phys_mocon.intfb.h"
#include "mf_phys_precips.intfb.h"
#include "mf_phys_save_phsurf_part1.intfb.h"
#include "mf_phys_save_phsurf_part2.intfb.h"
#include "mf_phys_transfer.intfb.h"
#include "ppwetpoint.intfb.h"
#include "profilechet.intfb.h"
#include "qngcor.intfb.h"
#include "radaer.intfb.h"
#include "radheat.intfb.h"
#include "radozcmf.intfb.h"
#include "recmwf.intfb.h"
#include "suozon.intfb.h"
#include "surf_ideal_flux.intfb.h"
#include "writephysio.intfb.h"
#include "wrphtrajm.intfb.h"
#include "wrphtrajtm_nl.intfb.h"
#include "wrradcoef.intfb.h"
#include "wrscmr.intfb.h"
#include "aplpar_flexdia.intfb.h"
#include "checknan.intfb.h"
#include "aclight.intfb.h"
#include "diagflash.intfb.h"

!     ------------------------------------------------------------------

#include "fcttrm.func.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('APLPAR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDVAB=>YDGEOMETRY%YRVAB, YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,      &
& YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1, YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2, YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH,           &
& YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL, YDRIP=>YDMODEL%YRML_GCONF%YRRIP, YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH,                  &
& YDRCOEF=>YDMODEL%YRML_PHY_RAD%YRRCOEF, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,             &
& YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, YGFL=>YDMODEL%YRML_GCONF%YGFL,                        &
& YDEPHY=> YDMODEL%YRML_PHY_EC%YREPHY, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR, YDPRECIPS=>YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS, &
& YDGEM=>YDGEOMETRY%YRGEM, YDSTA=>YDGEOMETRY%YRSTA, YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI, YDMCC=>YDMODEL%YRML_AOC%YRMCC,       &
& YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD, YDPHY3=>YDMODEL%YRML_PHY_MF%YRPHY3, YDPHY1=>YDMODEL%YRML_PHY_MF%YRPHY1,                &
& YDPHY0=>YDMODEL%YRML_PHY_MF%YRPHY0, YDNORGWD=>YDMODEL%YRML_PHY_MF%YRNORGWD, YDMSE=>YDMODEL%YRML_PHY_MF%YRMSE,               &
& YDPHYDS=>YDMODEL%YRML_PHY_MF%YRPHYDS, YDDYNA=>YDMODEL%YRML_DYN%YRDYNA,                                                      &
& YDSPP_CONFIG=>YDMODEL%YRML_GCONF%YRSPP_CONFIG, YDECUMF=>YDMODEL%YRML_PHY_EC%YRECUMF)



ASSOCIATE(CMF_UPDRAFT=>YDPARAR%PHYEX%PARAM_MFSHALLN%CMF_UPDRAFT, TSPHY=>YDPHY2%TSPHY, NTSSG=>YDDPHY%NTSSG, LMDUST=>YDARPHY%LMDUST,              &
& LMSE=>YDARPHY%LMSE, YI=>YGFL%YI, YEZDIAG=>YGFL%YEZDIAG, YL     =>YGFL%YL, YEXT=>YGFL%YEXT, YG=>YGFL%YG,                                 &
& YQ=>YGFL%YQ, YR=>YGFL%YR, YSCONV=>YGFL%YSCONV, YS=>YGFL%YS, YEFB3=>YGFL%YEFB3, YEFB2=>YGFL%YEFB2, YEFB1=>YGFL%YEFB1,                    &
& LCHEM_ARPCLIM=>YDMODEL%YRML_CHEM%YRCHEM%LCHEM_ARPCLIM, NGFL_EXT=>YGFL%NGFL_EXT, YTKE=>YGFL%YTKE, YLCONV=>YGFL%YLCONV,                   &
& YRCONV=>YGFL%YRCONV, YICONV=>YGFL%YICONV, YSP_SBD=>YDSURF%YSP_SBD, LTRAJPS=>YDSIMPHL%LTRAJPS, LNEBR=>YDPHY%LNEBR,                       &
& LNEBN=>YDPHY%LNEBN, LSTRAPRO=>YDPHY%LSTRAPRO, LPTKE=> YDPHY%LPTKE, NDPSFI=>YDPHY%NDPSFI, LOZONE=>YDPHY%LOZONE,                          &
& L3MT=>YDPHY%L3MT, LGPCMT=>YDPHY%LGPCMT, LAJUCV=>YDPHY%LAJUCV, LCVPGY=>YDPHY%LCVPGY, LRRGUST=>YDPHY%LRRGUST,                             &
& LEDR=>YDPHY%LEDR, NTAJUC=> YDTOPH%NTAJUC, NTPLUI=>YDTOPH%NTPLUI, LDPRECIPS=>YDPHY%LDPRECIPS, LDPRECIPS2=>YDPHY%LDPRECIPS2,              &
& LRCOEF    =>YDRCOEF%LRCOEF, NG3SR=>YDRCOEF%NG3SR, XMINLM=>YDPHY0%XMINLM, RTCAPE=>YDPHY0%RTCAPE, GRSO=>YDPHY0%GRSO,                      &
& GCVTSMO=>YDPHY0%GCVTSMO, XKLM=>YDPHY0%XKLM, GAEPS=>YDPHY0%GAEPS, AERCS1=>YDPHY0%AERCS1, AERCS3=>YDPHY0%AERCS3,                          &
& AERCS5=>YDPHY0%AERCS5, HUTIL2=>YDPHY0%HUTIL2, HUTIL1=>YDPHY0%HUTIL1, XMAXLM=>YDPHY0%XMAXLM, TEQK=>YDPHY0%TEQK,                          &
& HUCOE=>YDPHY0%HUCOE, LCVNHD=>YDPHY0%LCVNHD, TEQC=>YDPHY0%TEQC, UHDIFV=>YDPHY0%UHDIFV, HUTIL=>YDPHY0%HUTIL,                              &
& NPCLO1=>YDPHY0%NPCLO1, NPCLO2=>YDPHY0%NPCLO2, RDECRD=>YDPHY0%RDECRD, RDECRD1=>YDPHY0%RDECRD1, RDECRD2=>YDPHY0%RDECRD2,                  &
& RDECRD3=>YDPHY0%RDECRD3, RDECRD4=>YDPHY0%RDECRD4, ETKE_MIN=>YDPHY0%ETKE_MIN, HSOLIWR=>YDPHY1%HSOLIWR,                                   &
& ALCRIN=>YDPHY1%ALCRIN, ALBMED=>YDPHY1%ALBMED, WSMX=>YDPHY1%WSMX, LALBMERCLIM=>YDPHY1%LALBMERCLIM, HSOLIT0=>YDPHY1%HSOLIT0,              &
& HSOL=>YDPHY1%HSOL, WPMX=>YDPHY1%WPMX, EMCRIN=>YDPHY1%EMCRIN, EMMMER=>YDPHY1%EMMMER, TMERGL=>YDPHY1%TMERGL,                              &
& EMMGLA=>YDPHY1%EMMGLA, LRAFTKE=>YDPHY2%LRAFTKE, LRAFTUR=>YDPHY2%LRAFTUR, HVCLS=>YDPHY2%HVCLS, HTCLS=>YDPHY2%HTCLS,                      &
& FSM_HH=>YDPHY3%FSM_HH, FSM_GG=>YDPHY3%FSM_GG, FSM_FF=>YDPHY3%FSM_FF, FSM_EE=>YDPHY3%FSM_EE, FSM_II=>YDPHY3%FSM_II,                      &
& FSM_CC=>YDPHY3%FSM_CC, FSM_DD=>YDPHY3%FSM_DD, RLAMB_WATER=>YDPHY3%RLAMB_WATER, RII0=>YDPHY3%RII0, QCO2=>YDPHY3%QCO2,                    &
& RLAMB_SOLID=>YDPHY3%RLAMB_SOLID, NDLUNG=>YDDIM%NDLUNG, NDGUNG=>YDDIM%NDGUNG, NDLUXG=>YDDIM%NDLUXG,                                      &
& NDGUXG=>YDDIM%NDGUXG, LRDEPOS=>YDARPHY%LRDEPOS, LMPA=>YDARPHY%LMPA, CCOUPLING=>YDARPHY%CCOUPLING, LSURFEX_KFROM=>YDARPHY%LSURFEX_KFROM, &
& YA=>YGFL%YA, NGFL_EZDIAG=>YGFL%NGFL_EZDIAG, YFQTUR=>YGFL%YFQTUR, YFSTUR=>YGFL%YFSTUR, YIRAD=>YGFL%YIRAD,                                &
& YLRAD=>YGFL%YLRAD, XZSEPS=>YDMSE%XZSEPS, LVDIFSPNL=>YDSIMPHL%LVDIFSPNL, LGWDSPNL=>YDSIMPHL%LGWDSPNL,                                    &
& LRAYSP=>YDSIMPHL%LRAYSP, LSTRA=>YDPHY%LSTRA, LAEROSOO=>YDPHY%LAEROSOO, LCDDPRO=>YDPHY%LCDDPRO, LRKCDEV=>YDPHY%LRKCDEV,                  &
& LHUCN=>YDPHY%LHUCN, LCOEFK_TOMS=>YDPHY%LCOEFK_TOMS, LVDIF=>YDPHY%LVDIF, LRRMES=>YDPHY%LRRMES, LCVRA=>YDPHY%LCVRA,                       &
& LCVTDK=>YDPHY%LCVTDK, LCOEFK_RIS=>YDPHY%LCOEFK_RIS, LAEROLAN=>YDPHY%LAEROLAN, LNEBECT=>YDPHY%LNEBECT,                                   &
& LAERODES=>YDPHY%LAERODES, LNEWSTAT=>YDPHY%LNEWSTAT, LTHERMO=>YDPHY%LTHERMO, LO3FL=>YDPHY%LO3FL, LPHSPSH=>YDPHY%LPHSPSH,                 &
& LSNV=>YDPHY%LSNV, LECSHAL=>YDPHY%LECSHAL, LECT=>YDPHY%LECT, LDIFCONS=>YDPHY%LDIFCONS, LNODIFQC=>YDPHY%LNODIFQC,                         &
& LAEROVOL=>YDPHY%LAEROVOL, LRSTAER=>YDPHY%LRSTAER, NCALLRAD=>YDPHY%NCALLRAD, LNCVPGY=>YDPHY%LNCVPGY,                                     &
& LRAYLU=>YDPHY%LRAYLU, LAEROSUL=>YDPHY%LAEROSUL, LO3ABC=>YDPHY%LO3ABC, LSTRAS=>YDPHY%LSTRAS, LCOEFK_PTTE=>YDPHY%LCOEFK_PTTE,             &
& LSFHYD=>YDPHY%LSFHYD, LAEROSEA=>YDPHY%LAEROSEA, NDIFFNEB=>YDPHY%NDIFFNEB, LEDKF=>YDPHY%LEDKF, LRCVOTT=>YDPHY%LRCVOTT,                   &
& LMPHYS=>YDPHY%LMPHYS, LZ0HSREL=>YDPHY%LZ0HSREL, LCAMOD=>YDPHY%LCAMOD, LCOMOD=>YDPHY%LCOMOD, LCONDWT=>YDPHY%LCONDWT,                     &
& LCVCSD=>YDPHY%LCVCSD, LNSDO=>YDPHY%LNSDO, LUDEVOL=>YDPHY%LUDEVOL, LRAY=>YDPHY%LRAY, LGWD=>YDPHY%LGWD,                                   &
& LCOEFKTKE=>YDPHY%LCOEFKTKE, LRAYFM=>YDPHY%LRAYFM, LECDEEP=>YDPHY%LECDEEP, LCVGQD=>YDPHY%LCVGQD, LGWDC=>YDPHY%LGWDC,                     &
& LNORGWD=>YDPHY%LNORGWD, LFLUSO=>YDPHY%LFLUSO, LNEBCO=>YDPHY%LNEBCO, LNEBCV=>YDPHY%LNEBCV, LSOLV=>YDPHY%LSOLV,                           &
& LCOEFKSURF=>YDPHY%LCOEFKSURF, LRNUEXP=>YDPHY%LRNUEXP, NRAY=>YDPHY%NRAY, LDAYD=>YDPHY%LDAYD, LEDMFI=>YDPHY%LEDMFI,                       &
& LFPCOR=>YDPHY%LFPCOR, LRPROX=>YDPHY%LRPROX, LPROCLD=>YDPHY%LPROCLD, LACDIFUS=>YDPHY%LACDIFUS, LCAPE=>YDPHY%LCAPE,                       &
& LCVRAV3=>YDPHY%LCVRAV3, LHMTO=>YDPHY%LHMTO, LVGSN=>YDPHY%LVGSN, LCVPPKF=>YDPHY%LCVPPKF, LCVPRO=>YDPHY%LCVPRO,                           &
& LADJCLD=>YDPHY%LADJCLD, LGRAPRO=>YDPHY%LGRAPRO, RDECLI=>YDRIP%RDECLI, CMF_CLOUD=>YDPARAR%PHYEX%PARAM_MFSHALLN%CMF_CLOUD,                      &
& XSW_BANDS=>YDPARAR%XSW_BANDS, NSWB_MNH=>YDPARAR%NSWB_MNH, LMIXUV=>YDPARAR%PHYEX%PARAM_MFSHALLN%LMIXUV, NTRADI=>YDTOPH%NTRADI,                 &
& NTNEBU=>YDTOPH%NTNEBU, NTDIFU=>YDTOPH%NTDIFU, NTOZON=>YDTOPH%NTOZON, NTDRME=>YDTOPH%NTDRME, NTCVIM=>YDTOPH%NTCVIM,                      &
& NTCOET=>YDTOPH%NTCOET, NTDRAG=>YDTOPH%NTDRAG, RMESOQ=>YDTOPH%RMESOQ, RMESOT=>YDTOPH%RMESOT, RMESOU=>YDTOPH%RMESOU,                      &
& NTQSAT=>YDTOPH%NTQSAT, NTCOEF=>YDTOPH%NTCOEF, NAER=>YDERAD%NAER, NOZOCL=>YDERAD%NOZOCL, NRADFR=>YDERAD%NRADFR,                          &
& NAERMACC=>YDERAD%NAERMACC, &
& NSW=>YDERAD%NSW, RCARDI=>YDERDI%RCARDI, RSUNDUR=>YDERDI%RSUNDUR, LXVISI=>YDXFU%LXVISI, LXVISI2=>YDXFU%LXVISI2,                          &
& LMCC03=>YDMCC%LMCC03, NSTOP=>YDRIP%NSTOP, RCODEC=>YDRIP%RCODEC, RHGMT=>YDRIP%RHGMT, RSIDEC=>YDRIP%RSIDEC,                               &
& RSOVR=>YDRIP%RSOVR, RSTATI=>YDRIP%RSTATI, TSTEP=>YDRIP%TSTEP, NDTPREC=>YDPHY%YRDPRECIPS%NDTPREC, NDTPREC2=>YDPHY%YRDPRECIPS%NDTPREC2,   &
& STPRE=>YDSTA%STPRE, STPREH=>YDSTA%STPREH, STTEM=>YDSTA%STTEM, NORGWD_NNOVERDIF=>YDNORGWD%NORGWD_NNOVERDIF,                              &
& LGCHECKMV=>YDPHY%LGCHECKMV, NAERO=>YGFL%NAERO, NCHEM=>YDMODEL%YRML_GCONF%YGFL%NCHEM, NACTAERO=>YGFL%NACTAERO,                           &
 & RG=>YDCST%RG, RSIGMA=>YDCST%RSIGMA, RCPV=>YDCST%RCPV, RETV=>YDCST%RETV, &
 & RCW=>YDCST%RCW, RCS=>YDCST%RCS, RLVTT=>YDCST%RLVTT, RLSTT=>YDCST%RLSTT, &
 & RTT=>YDCST%RTT, RALPW=>YDCST%RALPW, RBETW=>YDCST%RBETW, RGAMW=>YDCST%RGAMW, &
 & RALPS=>YDCST%RALPS, RBETS=>YDCST%RBETS, RGAMS=>YDCST%RGAMS, RALPD=>YDCST%RALPD, &
 & RBETD=>YDCST%RBETD, RGAMD=>YDCST%RGAMD, RCPD=>YDCST%RCPD, RATM=>YDCST%RATM, &
 & RKAPPA=>YDCST%RKAPPA, RV=>YDCST%RV, RD=>YDCST%RD, &
& LXMRT=>YDXFU%LXMRT, RDELXN=>YDGEM%RDELXN, LFLASH =>YDCFU%LFLASH, CGMIXLEN=>YDMODEL%YRML_PHY_MF%YRPHY%CGMIXLEN                            )

CALL SC2PRG(YEFB1%MP1,  ZTENDGFL, ZPTENDEFB11)
CALL SC2PRG(YEFB2%MP1,  ZTENDGFL, ZPTENDEFB21)
CALL SC2PRG(YEFB3%MP1,  ZTENDGFL, ZPTENDEFB31)
CALL SC2PRG(YG%MP1,     ZTENDGFL, ZPTENDG1)
CALL SC2PRG(YICONV%MP1, ZTENDGFL, ZPTENDICONV1)
CALL SC2PRG(YI%MP1,     ZTENDGFL, ZPTENDI1)
CALL SC2PRG(YLCONV%MP1, ZTENDGFL, ZPTENDLCONV1)
CALL SC2PRG(YL%MP1,     ZTENDGFL, ZPTENDL1)
CALL SC2PRG(YQ%MP1,     ZTENDGFL, ZPTENDQ1)
CALL SC2PRG(YRCONV%MP1, ZTENDGFL, ZPTENDRCONV1)
CALL SC2PRG(YR%MP1,     ZTENDGFL, ZPTENDR1)
CALL SC2PRG(YSCONV%MP1, ZTENDGFL, ZPTENDSCONV1)
CALL SC2PRG(YS%MP1,     ZTENDGFL, ZPTENDS1)
CALL SC2PRG(YTKE%MP1,   ZTENDGFL, ZPTENDTKE1)

CALL SC2PRG(1, YEZDIAG(:)%MP, YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG, PGFL, ZP1EZDIAG)

!     ------------------------------------------------------------------

!        0.    constructor for procset
IF (LINTFLEX) YLPROCSET=NEWINTPROCSET()

!        1.    Preliminary calculations necessary
!              for all types of physics.
!              ------------------------------------


INSTEP_DEB=1
INSTEP_FIN=1

! SPP
! am:phasing 49 : dead piece of code ?
!IF ( YDSPP_CONFIG%LSPP ) THEN
! DO JSPP=1,YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL
!   ZGP2DSPP(:,JSPP) = YSPP%GP_ARP(JSPP)%GP2D(:,1,YDCPG_BNDS%KBL)
! ENDDO
!ENDIF

CALL CPPHINP(YDGEOMETRY, YDMODEL, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GELAM%T0,  &
& YDVARS%U%T0, YDVARS%V%T0, YDVARS%Q%T0, YDVARS%Q%DL, YDVARS%Q%DM, YDVARS%CVGQ%DL, YDVARS%CVGQ%DM, YDCPG_PHY0%XYB%RDELP,  &
& YDCPG_DYN0%CTY%EVEL, YDVARS%CVGQ%T0, ZRDG_MU0, ZRDG_MU0LU, ZRDG_MU0M, ZRDG_MU0N, ZRDG_CVGQ)
ZRDG_LCVQ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=ZRDG_CVGQ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)

DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
  ZFLU_QSATS(JROF)=0.0_JPRB
ENDDO

CALL MF_PHYS_FPL_PART1 (YDCPG_BNDS, YDCPG_OPTS, ZPFL_FPLCH, ZPFL_FPLSH, YDVARS%CPF%T0, YDVARS%SPF%T0, &
& YDMODEL)


! * In some cases, some pseudo-historic surface buffers (like z0) should
!   not be modified between the entrance and the output of APLPAR
!   (this is the case for example if LDCONFX=T).
!   For the time being, we must save:
!   - HV (group VV) : resistance to evapotranspiration
!   - Z0F (group VD): gravity * surface roughness length
!   - Z0H (group VV): gravity * roughness length for heat
!   - PBLH (group VH): PBL height
!   - SPSH (group VH):
!   - QSH (group VH):

LL_SAVE_PHSURF=YDCPG_OPTS%LCONFX
IF (LL_SAVE_PHSURF) THEN
  CALL MF_PHYS_SAVE_PHSURF_PART1 (YDCPG_BNDS, YDCPG_OPTS, ZSAV_DDAL, ZSAV_DDOM, ZSAV_ENTCH,                           &
  & ZSAV_FHPS, ZSAV_GZ0F, ZSAV_GZ0HF, ZSAV_HV, ZSAV_PBLH, ZSAV_QSH, ZSAV_UDAL, ZSAV_UDGRO, ZSAV_UDOM,                 &
  & ZSAV_UNEBH, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VH%PPBLH, YDMF_PHYS_SURF%GSD_VH%PQSH,                  &
  & YDMF_PHYS_SURF%GSD_VH%PSPSH, YDMF_PHYS_SURF%GSD_VK%PUDGRO, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VV%PZ0H, &
  & YDVARS%DAL%T0, YDVARS%DOM%T0, YDVARS%UAL%T0, YDVARS%UEN%T0, YDVARS%UNEBH%T0, YDVARS%UOM%T0,                       &
  & YDMODEL)
ENDIF


IF (LMDUST) THEN
  ZDIFEXT (:,:,:) = 0.0_JPRB
ENDIF

CALL APLPAR_INIT (YDCPG_OPTS%LAROME, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, NTSSG, &
& YSP_SBD%NLEVS, YDMF_PHYS_SURF%GSD_VF%PVEG, ZMSC_FRMQ, ZDSA_CPS, ZDSA_LHS, ZDSA_RS, ZMSC_LH, ZMSC_LSCPE, ZFLU_QSAT,&
& ZMSC_QW, ZMSC_TW, ZFLU_CD, ZFLU_CDN, ZFLU_CH, ZDSA_C1, ZDSA_C2, ZFLU_EMIS, ZFLU_FEVI, ZPFL_FTKE,                  &
& ZPFL_FTKEI, ZPFL_FEFB1, ZPFL_FEFB2, ZPFL_FEFB3, ZFLU_NEIJ, ZFLU_VEG, ZFLU_QSATS, IMOC_CLPH)

!*       2.    Complete physics.
!              -----------------

!        2.2  Complete physics.
!             -----------------

! PAS DE TEMPS DE LA PHYSIQUE (/YOMPHY2/)
! Dans le cas des iterations de l'initialisation par modes normaux,
! le pas de temps pour la physique ne peux pas etre nul pour APLPAR et
! CPATY (par contre c'est bien PDTPHY qui est passe en argument aux autres
! sous-prog. de la physique). Ceci est du a l'impossibilite de prendre en
! compte des flux qui deviennent infinis pour TSPHY=0 (flux de masse du au
! reajustement des sursaturations par exemple...). Mais les tendances phys.
! sont bien nulles dans le cas de la configuration 'E' (Modes Normaux).
! PHYSICS TIME STEP (/YOMPHY2/)
! In case of normal mode initialisation iterations, the physics time
! step cannot be zero for APLPAR and CPATY (nevertheless it is PDTPHY
! which is passed to other physics subroutines). This is due to the
! impossibility to take into account fluxes which are infinite for TSPHY=0
! (e.g.: mass flux due to oversaturation...).

! CALL PARAMETERISATIONS

DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
  YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)=REAL(NINT(YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)),JPRB)
ENDDO

IF (YDDYNA%LTWOTL) THEN
  IF (LAJUCV) THEN
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZADJ_TAUX(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)
      ENDDO
    ENDDO
    CALL ACAJUCV(YDMODEL%YRML_PHY_MF%YRPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,     &
    & NTPLUI, YDCPG_OPTS%KFLEVG, NTAJUC, YDCPG_PHY0%PREHYD, YDCPG_PHY0%XYB%ALPH, YDCPG_PHY0%XYB%DELP, &
    & YDCPG_PHY0%XYB%LNPR, YDVARS%T%T0)
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZADJ_DTAJU(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)-ZADJ_TAUX(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
ELSE
   ! IF (LAJUCV) THEN
   !   missing code under LAJUCV for leap-frog schemes.
   ! ENDIF
ENDIF


!
!-------------------------------------------------
! Check magnitude of model variables.
!-------------------------------------------------
!
IF(LGCHECKMV) CALL CHECKMV(YDCPG_OPTS%NINDAT, YDCST, YDRIP, YDPHY0, YDPHY2, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, &
              & YDCPG_OPTS%KFLEVG, YDCPG_OPTS%NSTEP, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,              &
              & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDVARS%GEOMETRY%GELAM%T0,                   &
              & YDVARS%GEOMETRY%GEMU%T0, ZRDG_MU0, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_BASE_STATE%T,                                   &
              & YDMF_PHYS_BASE_STATE%Q, YDMF_PHYS_BASE_STATE%YGSP_RR%T                      )
!     ------------------------------------------------------------------

LLREDPR=LCVCSD
ZRVMD=RV-RD
! SURFEX  and passive scalar
IF (YDCPG_OPTS%LCONFX) THEN
  ZDTMSE=0.01_JPRB
  ZSTATI=RSTATI-ZDTMSE*.5_JPRB
  ZADTMS=0._JPRB
ELSE
  ZDTMSE=TSPHY
  ZSTATI=RSTATI
  ZADTMS=ZDTMSE
ENDIF
ZRHGMT=REAL(RHGMT,JPRB)
ZAIPCMT(:)=0._JPRB

ALLOCATE(ZSVM   (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,NGFL_EXT))
ALLOCATE(ZSFSV  (YDCPG_OPTS%KLON,NGFL_EXT))
ALLOCATE(ZPSV   (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,NGFL_EXT))

!     ------------------------------------------------------------------
!     1.- INITIALISATIONS COMPLEMENTAIRES
!     -----------------------------------
IJN = YDCPG_OPTS%KLON

!     CALCUL FIN DE IJN
!      IJN = 1
!      ZMU0 = MAX(0.,-SIGN(1.,-PMU0(KIDIA)))
!      DO 1  JLON = KIDIA+1, KFDIA
!      ZMUN = MAX(0.,-SIGN(1.,-PMU0(JLON)))
!      IJN = IJN + IABS(NINT(ZMUN)-NINT(ZMU0))
!      ZMU0 = ZMUN
!   1  CONTINUE

!*        1.0 DECORRELATION DEPTH FOR CLOUD OVERLAPS
IF ( RDECRD <= 0._JPRB .OR. LRNUEXP ) THEN
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZDECRD(JLON)=RDECRD1+RDECRD2* &
     & EXP(-((ASIN(YDVARS%GEOMETRY%GEMU%T0(JLON))-RDECRD3*RDECLI)/RDECRD4)**2)
  ENDDO
ENDIF
IF ( RDECRD <= 0._JPRB ) THEN 
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZDECRD_MF(JLON)=ZDECRD(JLON)
  ENDDO
ELSE
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZDECRD_MF(JLON)=RDECRD
  ENDDO
ENDIF

!*        1.1 INITIALISATION DE L'OZONE

IF (LMPHYS) THEN
  IF (LOZONE) THEN

    ! L'ozone est initialise par PO3 (ozone GFL).
    ! Ozone is computed from PO3 (GFL ozone).

    ZEPSO3=1.E-11_JPRB
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZQO3(JLON,0)=1.E-9_JPRB
    ENDDO
    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZQO3(JLON,JLEV) = MAX(ZEPSO3,YDMF_PHYS_BASE_STATE%O3(JLON,JLEV))
      ENDDO
    ENDDO

  ELSEIF(YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS == 1) THEN
    ZEPSO3=1.E-11_JPRB
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZQO3(JLON,0)=1.E-9_JPRB
    ENDDO
    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZQO3(JLON,JLEV) = YDCPG_MISC%KOZO(JLON,JLEV,1)
      ENDDO
    ENDDO

  ELSEIF ((LO3FL).AND.(NOZOCL == 1).AND.(LRAYFM)) THEN
    IF (MOD(YDCPG_OPTS%NSTEP,NRADFR) == 0) THEN
      CALL RADOZCMF(YDMODEL%YRCST,YDRIP%YREOZOC, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
      & YDVARS%GEOMETRY%GEMU%T0, ZROZ)
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZQO3(JLON,0)=1.E-9_JPRB
      ENDDO
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZQO3(JLON,JLEV)=ZROZ(JLON,JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ELSE
    CALL SUOZON(YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, ZQO3,         &
    & .FALSE., YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, LO3ABC, &
    & YDMF_PHYS_SURF%GSD_VC%PA, YDMF_PHYS_SURF%GSD_VC%PB, YDMF_PHYS_SURF%GSD_VC%PC)
  ENDIF

!     GAZ CARBONIQUE.

  DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZQCO2(JLON,JLEV)=QCO2
    ENDDO
  ENDDO


  ! INITIALISATION DE LA COORDONNEE ETA.
  ! INITIALISATION DE LA COORDONNEE ETA.

!     EPAISSEUR STD AEROSOLS
  ZAEO = AERCS1*YDSTA%SVETAH(YDCPG_OPTS%KTDIA-1) + AERCS3*YDSTA%SVETAH(YDCPG_OPTS%KTDIA-1)**3&
   & + AERCS5*YDSTA%SVETAH(YDCPG_OPTS%KTDIA-1)**5
  DO JLEV = YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG
    ZAEN = AERCS1*YDSTA%SVETAH(JLEV) + AERCS3*YDSTA%SVETAH(JLEV)**3&
     & + AERCS5*YDSTA%SVETAH(JLEV)**5
    ZDAER(JLEV) = ZAEN - ZAEO
    ZAEO = ZAEN
  ENDDO

  IF ( LNEWSTAT ) THEN

    DO JLEV = YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG-1
      ZMAN(JLEV) = FSM_CC * TANH(FSM_DD*YDSTA%SVETAH(JLEV))
      ZMAK(JLEV) = FSM_EE * YDSTA%SVETAH(JLEV)**FSM_FF +&
       & FSM_GG * (1-YDSTA%SVETAH(JLEV))**FSM_HH + FSM_II
    ENDDO

!     MATHEMATICAL FILTER FOR THE EDGES

    IF ( LRPROX ) THEN
      DO JLEV = YDCPG_OPTS%KTDIA, YDCPG_OPTS%KTDIA+3
        ZMAK(JLEV) = ZMAK(JLEV) / 2**(YDCPG_OPTS%KTDIA+4-JLEV)
      ENDDO
      DO JLEV = YDCPG_OPTS%KFLEVG-4, YDCPG_OPTS%KFLEVG-1
        ZMAK(JLEV) = ZMAK(JLEV) / 2**(5-YDCPG_OPTS%KFLEVG+JLEV)
      ENDDO
    ENDIF

  ELSE
    DO JLEV = YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG-1
      ZMAN(JLEV) = 0.3_JPRB*YDSTA%SVETAH(JLEV)
      ZMAK(JLEV) = 0.1_JPRB
    ENDDO
  ENDIF   !LNEWSTAT

!3MT
!     INCREMENTAL CORRECTION FLUX FOR NEGAVTIVE HUMIDITY VALUES

  ZFCQVNG(:,:)=0.0_JPRB
  ZFCQING(:,:)=0.0_JPRB
  ZFCQLNG(:,:)=0.0_JPRB

  ZPRODTH_CVPP(:,:)=0.0_JPRB

  ZGDT=RG*TSPHY
  ZGDTI=1.0_JPRB/ZGDT

  DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA

!     LAYER WEIGHTS

      ZPOID(JLON,JLEV)=YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZGDTI
      ZIPOI(JLON,JLEV)=YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)*ZGDT

!     CALCULATION OF LATENT HEATS

      ZLHV(JLON,JLEV)=FOLH(YDMF_PHYS_BASE_STATE%T(JLON,JLEV),0.0_JPRB)
      ZLHS(JLON,JLEV)=FOLH(YDMF_PHYS_BASE_STATE%T(JLON,JLEV),1.0_JPRB)

    ENDDO
  ENDDO

!    ------------------------------------------------------------------
!     PROGNOSTIC GEMS/MACC AEROSOLS - INITIAL COMPUTATIONS
!     IMPORTANT for IFS: Tracer order is : CO2 - other tracers - react Gases - Aerosol - extra GFL
!    ------------------------------------------------------------------

  ! Preliminary for prog. aerosol or extra gfl
  INBTRA=0
  IF (LMDUST.AND.(NGFL_EXT/=0)) INBTRA=NGFL_EXT
  IF (NAERO>0)                  INBTRA=NAERO           ! the two cases exclude each other
  ALLOCATE(ZSTRCTRA(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,INBTRA)) ! to cover both prog aero and extra gfl cases
  ALLOCATE(ZTRA    (YDCPG_OPTS%KLON,  YDCPG_OPTS%KFLEVG,INBTRA))
  IF (INBTRA > 0) THEN
    ZSTRCTRA(:,:,:) = 0._JPRB
    ZTRA(:,:,:)     = 0._JPRB
  ENDIF 
  IF(INBTRA == 0) THEN
    INBTRA_DEP=0
  ELSE
    INBTRA_DEP=1
  ENDIF

!    ------------------------------------------------------------------
!      - CHANGEMENTS DE VARIABLES ET INVERSION DES NIVEAUX
!          POUR LE TRAITEMENT DES SCALAIRES PASSIFS
!     --------------------------------------------------------------------

  ZSFSV=0.0_JPRB    ! surf. flux of scalars
  IF (LMDUST.AND.(NGFL_EXT/=0)) THEN
!  SIZE OF ARRAY FOR MSE
  ILONMNH=YDCPG_BNDS%KFDIA-YDCPG_BNDS%KIDIA+1
  ZINVDT=1/YDCPG_OPTS%ZDTPHY
  ZINVG=1/RG
  ZZI_APHI=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI
  ZTM=YDMF_PHYS_BASE_STATE%T
!SETUP
  ZZI_SVM=0.0_JPRB
  ZZI_PEZDIAG=0.0_JPRB
  ZZI_PABSM=101325._JPRB
  ZZI_RHODREFM=1.0_JPRB
  ZZI_RHO=1.0_JPRB
  ZZZ=0.0_JPRB
  ZAERD=0.0_JPRB
  ZP1EZDIAG=0.0_JPRB

!Initialisation des scalaires passifs pour aro_ground_param
  DO JGFL=1,NGFL_EXT
    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
       ZSVM(JLON,JLEV,JGFL)=MAX(YDMF_PHYS_BASE_STATE%P1EXT(JLON,JLEV,JGFL),0.0_JPRB)
      ENDDO
    ENDDO
  ENDDO

  !initialisation de ZZZ
  DO JLEV = 1,YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZZZ(JLON,1,JLEV)=ZZI_APHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO

  !initialisation de ZDZZ
  DO JLEV = 2, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZDZZ(JLON,1,JLEV)=ZZZ(JLON,1,JLEV-1)-ZZZ(JLON,1,JLEV)
    ENDDO
  ENDDO
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZDZZ(JLON,1,1)=ZZI_APHI(JLON,0)*ZINVG-ZZZ(JLON,1,1)
  ENDDO


!Initialisation de ZZI_RHODREFM
  DO JLEV = 1 , YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZZI_RHODREFM(JLON,1,JLEV)=YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)/&
       & (YDMF_PHYS_BASE_STATE%T(JLON,JLEV)*YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(JLON,JLEV))
     ENDDO
  ENDDO

!Initialisation de ZZI_PABSM
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    DO JLEV = 1 , YDCPG_OPTS%KFLEVG
      ZZI_PABSM(JLON,1,JLEV)=YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)
    ENDDO
  ENDDO

!Initialisation de ZZI_EXNREFM
  ZRSCP=RD/RCPD
  ZINVATM=1/RATM
  DO JLEV = 1 , YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZZI_EXNREFM(JLON,1,JLEV)=(YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)*ZINVATM)**(ZRSCP)
    ENDDO
  ENDDO

!Initialisation de ZZI_THM
  DO JLEV = 1 , YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZZI_THM(JLON,1,JLEV)=ZTM(JLON,JLEV)/ZZI_EXNREFM(JLON,1,JLEV)
    ENDDO
  ENDDO

!Initialisation des scalaires passifs pour aro_mnhdust (inversion des niveaux)
  DO JGFL=1,NGFL_EXT
    DO JLON=1,YDCPG_OPTS%KLON
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        ZZI_SVM(JLON,1,JLEV,JGFL)=ZSVM(JLON,JLEV,JGFL)
      ENDDO
    ENDDO
  ENDDO

  ENDIF  ! ENDIF  (LMDUST & NGFL_EXT)
ENDIF  !LMPHYS
!**
!     ------------------------------------------------------------------
!     2.- MISES A ZERO DE SECURITE EN CAS DE NON-APPEL DES PARAMETRIS.
!     ------------------------------------------------------------------
ZEPS0=1.E-12_JPRB
ZEPSNEB=1.E-10_JPRB

! To profitize from the vectorization collapsing the (:,:) form is preferable.
! (Even better would be to completely avoid any useless initialization.)

! arrays dimensioned from 0:KLEV (half level quantities)
ZFPCOR  (:,:) = 0.0_JPRB
ZFHP    (:,:) = 0.0_JPRB
ZXTROV  (:,:) = 1.0_JPRB
ZXUROV  (:,:) = 1.0_JPRB
ZLMT    (:,:) = 0.0_JPRB
ZZLMT   (:,:) = 0.0_JPRB
ZLMU    (:,:) = 0.0_JPRB
ZLMU2   (:,:) = 0.0_JPRB
ZLMT2   (:,:) = 0.0_JPRB
ZKTROV  (:,:) = 0.0_JPRB
ZKQROV  (:,:) = 0.0_JPRB
ZKQLROV (:,:) = 0.0_JPRB
ZKUROV  (:,:) = 0.0_JPRB
ZFHEVPPC(:,:) = 0.0_JPRB
ZFHMLTSC(:,:) = 0.0_JPRB
ZFPEVPPC(:,:) = 0.0_JPRB
ZFCQL   (:,:) = 0.0_JPRB
ZFCQI   (:,:) = 0.0_JPRB
ZDIFCVPPQ (:,:) = 0.0_JPRB
ZDIFCVPPS (:,:) = 0.0_JPRB
ZDIFCVTH (:,:) = 0.0_JPRB
ZDIFCVPPU (:,:) = 0.0_JPRB
ZDIFCVPPV (:,:) = 0.0_JPRB
ZCONDCVPPL(:,:) = 0.0_JPRB
ZCONDCVPPI(:,:) = 0.0_JPRB
ZSEDIQL(:,:) = 0.0_JPRB
ZSEDIQI(:,:) = 0.0_JPRB

ZXURO   (:,:) = 0.0_JPRB
ZXQRO   (:,:) = 0.0_JPRB
ZXTRO   (:,:) = 0.0_JPRB

ZALPHA1 (:,:) = 0.0_JPRB
ZCOEFA  (:,:) = 0.0_JPRB
ZLVT    (:,:) = 0.0_JPRB
ZQICE   (:,:) = 0.0_JPRB

ZF_EPS (:,:) = 1.0_JPRB
ZFUN_TTE (:,:) = 1.0_JPRB
ZMRIPP (:,:) = 1.E-12_JPRB
ZMRIMC  (:,:) = 1.0_JPRB
ZMRICTERM (:,:) = 1.0_JPRB
ZRRCOR (:,:) = 1.0_JPRB
ZTAU_TKE (:,:) = 0.0_JPRB
ZTH_FUN (:,:) = 1.0_JPRB
ZMRIFPP (:,:) = 1.E-12_JPRB
ZMN2PP (:,:) = 1.E-12_JPRB
ZMN2_ES (:,:) = 1.0_JPRB
ZMN2_EQ (:,:) = 1.0_JPRB
ZMN2_DS (:,:) = 1.0_JPRB
ZMN2_DQ (:,:) = 1.0_JPRB
ZTSTAR (:,:) = 1.E-12_JPRB
ZTSTAR2 (:,:) = 1.E-12_JPRB
ZTSTARQ (:,:) = 1.E-12_JPRB
ZTSTAR2Q (:,:) = 1.E-12_JPRB
ZFMGST (:,:) = 1.0_JPRB
ZFMTKE (:,:) = 1.0_JPRB
ZFTTKE (:,:) = 1.0_JPRB
ZAUTKE (:,:) = 1.0_JPRB
ZATTKE (:,:) = 1.0_JPRB
ZTH_FUN(:,:) = 1.0_JPRB
ZWW_FUN(:,:) = 1.0_JPRB
ZBNEBCVPP(:,:) = 0.0_JPRB
ZBNEBQ(:,:)   = 0.0_JPRB
ZRHS(:,:)     = 0.0_JPRB
ZLML(:,:)     = 1.0_JPRB
ZLMLTILD(:,:) = 1.0_JPRB

ZDIFWQ  (:) = 0.0_JPRB
ZDIFWS  (:) = 0.0_JPRB
ZSC_FEVI(:) = 1.0_JPRB     
ZSC_FEVN(:) = 1.0_JPRB     
ZSC_FCLL(:) = 1.0_JPRB     
ZSC_FCLN(:) = 1.0_JPRB
ZCDNH(:)    = 1.0_JPRB

! arrays dimensioned from 1:KLEV (full level quantities)
ZTENT   (:,:) = 0.0_JPRB
ZNEBS   (:,:) = ZEPS0
ZNEBC   (:,:) = ZEPS0
ZNEBS0  (:,:) = ZEPS0
ZNEBC0  (:,:) = ZEPS0
ZNEBCH  (:,:) = 0.0_JPRB
ZUNEBH  (:,:) = 0.0_JPRB
ZDETFI (:,:) = 0.0_JPRB
ZNEBDIFF(:,:) = 0.0_JPRB
ZQLIS   (:,:) = 0.0_JPRB
ZQLIS0  (:,:) = 0.0_JPRB
ZCFATH  (:,:) = 0.0_JPRB
ZCFAU   (:,:) = 0.0_JPRB
ZCFBTH  (:,:) = 0.0_JPRB
ZCFBU   (:,:) = 0.0_JPRB
ZCFBV   (:,:) = 0.0_JPRB
ZQLIC   (:,:) = 0.0_JPRB
INLAB   (:,:) = 0
INLAB_CVPP(:,:) = 0
ICIS    (:,:) = 1
ZQLI_CVPP(:,:) = 0.0_JPRB
ZNEB_CVPP(:,:) = ZEPS0

ZEDMFQ  (:,:)  = 0.0_JPRB
ZEDMFS  (:,:)  = 0.0_JPRB
ZEDMFU  (:,:)  = 0.0_JPRB
ZEDMFV  (:,:)  = 0.0_JPRB
ZMF_UP  (:,: ) = 0.0_JPRB
ZQLI_CVP(:,:) = 0.0_JPRB
ZQC_DET_PCMT(:,:) = 0.0_JPRB
ZTENHA(:,:)   = 0.0_JPRB
ZTENQVA(:,:)  = 0.0_JPRB
ZRHDFDA(:,:)  = 0.0_JPRB
ZQIC   (:,:)  = 0.0_JPRB
ZQLC   (:,:)  = 0.0_JPRB
ZQRC   (:,:)  = 0.0_JPRB
ZQSC   (:,:)  = 0.0_JPRB
ZQG    (:,:)  = 0.0_JPRB
ZQH    (:,:)  = 0.0_JPRB

!  ---------------------------------------------------------------------
!  Correction of negative advected humidity and precipitation values
!  ---------------------------------------------------------------------

IF (LCONDWT) THEN
  IF (L3MT .OR. LSTRAPRO .OR. LPROCLD) THEN
    IF (LGRAPRO) THEN
!cdir unroll=8
      DO JLEV = 1, YDCPG_OPTS%KFLEVG
        DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZQI(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%I(JLON,JLEV))
          ZQL(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%L(JLON,JLEV))
          ZQR(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%R(JLON,JLEV))
          ZQS(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%S(JLON,JLEV))
          ZQG(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%G(JLON,JLEV))
  ! CORRECTION OF NEGATIVE ADVECTED VALUES:
  !                         VAPOUR PUT IN PFCQVNG
  !                         LIQUID,ICE PUT IN PFCQL/ING
  !                         FOR RAIN/SNOW/GRAUPEL PUT IN PFCQR/S/GNG

          ZDQI=ZQI(JLON,JLEV)-YDMF_PHYS_BASE_STATE%I(JLON,JLEV)
          ZDQL=ZQL(JLON,JLEV)-YDMF_PHYS_BASE_STATE%L(JLON,JLEV)
          ZDQR=ZQR(JLON,JLEV)-YDMF_PHYS_BASE_STATE%R(JLON,JLEV)
          ZDQS=ZQS(JLON,JLEV)-YDMF_PHYS_BASE_STATE%S(JLON,JLEV)
          ZDQG=ZQG(JLON,JLEV)-YDMF_PHYS_BASE_STATE%G(JLON,JLEV)
          ZDQC=ZDQI+ZDQL+ZDQR+ZDQS+ZDQG

          ZQV0=YDMF_PHYS_BASE_STATE%Q(JLON,JLEV)-ZIPOI(JLON,JLEV)*(0.0_JPRB- YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1) &
          & -YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1) &
          & -YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQGNG(JLON,JLEV-1))
          ZQV(JLON,JLEV)=MAX(0.0_JPRB,ZQV0-ZDQC)
          ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-YDMF_PHYS_BASE_STATE%Q(JLON,JLEV)

          YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQRNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1)-ZDQR*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQSNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1)-ZDQS*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQGNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQGNG(JLON,JLEV-1)-ZDQG*ZPOID(JLON,JLEV)
        ENDDO
      ENDDO
    ELSE
      DO JLEV = 1, YDCPG_OPTS%KFLEVG
        DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZQI(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%I(JLON,JLEV))
          ZQL(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%L(JLON,JLEV))
          ZQR(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%R(JLON,JLEV))
          ZQS(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%S(JLON,JLEV))

    ! CORRECTION OF NEGATIVE ADVECTED VALUES:
    !                         VAPOUR PUT IN PFCQVNG
    !                         LIQUID,ICE PUT IN PFCQL/ING
    !                         FOR RAIN/SNOW PUT IN PFCQR/SNG

          ZDQI=ZQI(JLON,JLEV)-YDMF_PHYS_BASE_STATE%I(JLON,JLEV)
          ZDQL=ZQL(JLON,JLEV)-YDMF_PHYS_BASE_STATE%L(JLON,JLEV)
          ZDQR=ZQR(JLON,JLEV)-YDMF_PHYS_BASE_STATE%R(JLON,JLEV)
          ZDQS=ZQS(JLON,JLEV)-YDMF_PHYS_BASE_STATE%S(JLON,JLEV)
          ZDQC=ZDQI+ZDQL+ZDQR+ZDQS

          ZQV0=YDMF_PHYS_BASE_STATE%Q(JLON,JLEV)-ZIPOI(JLON,JLEV)*(0.0_JPRB- YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1)&
          & -YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1)&
          & -YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1))
          ZQV(JLON,JLEV)=MAX(0.0_JPRB,ZQV0-ZDQC)
          ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-YDMF_PHYS_BASE_STATE%Q(JLON,JLEV)

          YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQRNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1)-ZDQR*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQSNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1)-ZDQS*ZPOID(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

    IF (LGPCMT) THEN
      DO JLEV = 1, YDCPG_OPTS%KFLEVG
        DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZTZER=MAX(0.0_JPRB,-SIGN(1.0_JPRB,-(YDVARS%ICONV%T0(JLON,JLEV)-ZEPS0)+0.0_JPRB))
          ZQIC(JLON,JLEV)=ZTZER*YDVARS%ICONV%T0(JLON,JLEV)
          ZDQI=ZQIC(JLON,JLEV)-YDVARS%ICONV%T0(JLON,JLEV)

          ZTZER=MAX(0.0_JPRB,-SIGN(1.0_JPRB,-(YDVARS%LCONV%T0(JLON,JLEV)-ZEPS0)+0.0_JPRB))
          ZQLC(JLON,JLEV)=ZTZER*YDVARS%LCONV%T0(JLON,JLEV)
          ZDQL=ZQLC(JLON,JLEV)-YDVARS%LCONV%T0(JLON,JLEV)

          ZTZER=MAX(0.0_JPRB,-SIGN(1.0_JPRB,-(YDVARS%RCONV%T0(JLON,JLEV)-ZEPS0)+0.0_JPRB))
          ZQRC(JLON,JLEV)=ZTZER*YDVARS%RCONV%T0(JLON,JLEV)
          ZDQR=ZQRC(JLON,JLEV)-YDVARS%RCONV%T0(JLON,JLEV)

          ZTZER=MAX(0.0_JPRB,-SIGN(1.0_JPRB,-(YDVARS%SCONV%T0(JLON,JLEV)-ZEPS0)+0.0_JPRB))
          ZQSC(JLON,JLEV)=ZTZER*YDVARS%SCONV%T0(JLON,JLEV)
          ZDQS=ZQSC(JLON,JLEV)-YDVARS%SCONV%T0(JLON,JLEV)

          ZDQC=ZDQI+ZDQL+ZDQR+ZDQS

          ZQV0=ZQV(JLON,JLEV)-ZIPOI(JLON,JLEV)*(0.0_JPRB-ZFCQVNG(JLON,JLEV-1)&
          & -YDMF_PHYS%OUT%FCNEGQIC(JLON,JLEV-1)-YDMF_PHYS%OUT%FCNEGQLC(JLON,JLEV-1)&
          & -YDMF_PHYS%OUT%FCNEGQRC(JLON,JLEV-1)-YDMF_PHYS%OUT%FCNEGQSC(JLON,JLEV-1))
          ZQVI=MAX(0.0_JPRB,ZQV0-ZDQC)
          ZDQV=ZQVI-ZQV(JLON,JLEV)
          ZQV(JLON,JLEV)=ZQVI

          ZFCQVNG(JLON,JLEV)=ZFCQVNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV)+ZFCQVNG(JLON,JLEV)
          YDMF_PHYS%OUT%FCNEGQIC(JLON,JLEV)=YDMF_PHYS%OUT%FCNEGQIC(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCNEGQLC(JLON,JLEV)=YDMF_PHYS%OUT%FCNEGQLC(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCNEGQRC(JLON,JLEV)=YDMF_PHYS%OUT%FCNEGQRC(JLON,JLEV-1)-ZDQR*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCNEGQSC(JLON,JLEV)=YDMF_PHYS%OUT%FCNEGQSC(JLON,JLEV-1)-ZDQS*ZPOID(JLON,JLEV)
        ENDDO
      ENDDO
      ZFCQVNG(:,:)=0.0_JPRB
    ENDIF
  ELSE
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZQI(JLON,JLEV)=YDMF_PHYS_BASE_STATE%I(JLON,JLEV)
        ZQL(JLON,JLEV)=YDMF_PHYS_BASE_STATE%L(JLON,JLEV)
        ZQR(JLON,JLEV)=YDMF_PHYS_BASE_STATE%R(JLON,JLEV)
        ZQS(JLON,JLEV)=YDMF_PHYS_BASE_STATE%S(JLON,JLEV)
        IF (LGRAPRO) THEN
          ZQG(JLON,JLEV)=YDMF_PHYS_BASE_STATE%G(JLON,JLEV)
        ENDIF
        ZQV(JLON,JLEV)=YDMF_PHYS_BASE_STATE%Q(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF
ELSE
  DO JLEV = 1, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZQI(JLON,JLEV)=0.0_JPRB
      ZQL(JLON,JLEV)=0.0_JPRB
      ZQV(JLON,JLEV)=YDMF_PHYS_BASE_STATE%Q(JLON,JLEV)
    ENDDO
  ENDDO
ENDIF ! LCONDWT

DO JCHA = 1, 6
  DO JLEV = YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG
    DO JLON = 1, YDCPG_OPTS%KLON
      ZAER(JLON,JLEV,JCHA)=0.0_JPRB
    ENDDO
  ENDDO
ENDDO

DO JLEV = YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG
  DO JLON = 1, YDCPG_OPTS%KLON
    ZAERINDS(JLON,JLEV)=0.0_JPRB
  ENDDO
ENDDO

DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
  ZCEMTR  (JLON,0) = 0.0_JPRB
  ZCEMTR  (JLON,1) = 0.0_JPRB
  ZCTRSO  (JLON,0) = 0.0_JPRB
  ZCTRSO  (JLON,1) = 0.0_JPRB
  ZTRSOD  (JLON)   = 0.0_JPRB
  ZSUDU   (JLON)   = 0.0_JPRB
  ZXDROV  (JLON)   = 1.0_JPRB
  ZXHROV  (JLON)   = 1.0_JPRB
  ZTAUX   (JLON)   = ZEPS0
  IMOC_CLPH   (JLON)   = YDCPG_OPTS%KFLEVG
ENDDO
DO JSG = 1, NSW
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZALBD   (JLON,JSG) = 0.0_JPRB
    ZALBP   (JLON,JSG) = 0.0_JPRB
    ZSFSWDIF(JLON,JSG) = 0.0_JPRB
    ZSFSWDIR(JLON,JSG) = 0.0_JPRB
    ZTRSODIF(JLON,JSG) = 0.0_JPRB
    ZTRSODIR(JLON,JSG) = 0.0_JPRB
  ENDDO
ENDDO

!  -------------------------------------------------------
!  Security values for pseudo-historical arrays at KSTEP=0
!  -------------------------------------------------------

IF((LNEBR.OR.(TRIM(CGMIXLEN) == 'TM')&
       & .OR.(TRIM(CGMIXLEN) == 'TMC')).AND.YDCPG_OPTS%NSTEP == 0) THEN
!DEC$ IVDEP
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    YDMF_PHYS_SURF%GSD_VH%PPBLH(JLON)=(XMINLM+XMAXLM)*0.5_JPRB
  ENDDO
ENDIF
IF((LNEBCO.OR.LGWDC).AND.YDCPG_OPTS%NSTEP == 0) THEN
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    YDMF_PHYS_SURF%GSD_VH%PTCCH(JLON)=0.0_JPRB
    YDMF_PHYS_SURF%GSD_VH%PSCCH(JLON)=0.0_JPRB
    YDMF_PHYS_SURF%GSD_VH%PBCCH(JLON)=0.0_JPRB
  ENDDO
ENDIF

IF((LNEBN.OR.LNEBR.OR.LRRGUST).AND.YDCPG_OPTS%NSTEP == 0) THEN
  DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZPFL_FPLCH(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF
IF(LRRGUST.AND.YDCPG_OPTS%NSTEP == 0) THEN
  DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZPFL_FPLSH(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF
IF(LCVPGY.AND.YDCPG_OPTS%NSTEP == 0) THEN
  DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS_BASE_STATE%CVV(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF
IF(LPHSPSH.AND.YDCPG_OPTS%NSTEP == 0) THEN
  YDMF_PHYS_SURF%GSD_VH%PSPSH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
ENDIF
IF(LCOEFKTKE.AND.YDCPG_OPTS%NSTEP == 0) THEN
  YDMF_PHYS_SURF%GSD_VH%PQSH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_BASE_STATE%Q(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
ENDIF
IF(LCVCSD.AND.LUDEVOL.AND.YDCPG_OPTS%NSTEP==0) THEN
  YDMF_PHYS_SURF%GSD_VK%PUDGRO(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
ENDIF
IF(LRKCDEV.AND.YDCPG_OPTS%NSTEP == 0) THEN
  DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
     DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZDTRAD(JLON,JLEV)=0.0_JPRB
        ZRKTTEND(JLON,JLEV)=0.0_JPRB
        ZRKQVTEND(JLON,JLEV)=0.0_JPRB
        ZRKQCTEND(JLON,JLEV)=0.0_JPRB
        ZDQVDIFF(JLON,JLEV)=0.0_JPRB
        YDVARS%RKTH%T0(JLON,JLEV) = 0.0_JPRB
        YDVARS%RKTQV%T0(JLON,JLEV)= 0.0_JPRB
        YDVARS%RKTQC%T0(JLON,JLEV)= 0.0_JPRB
     ENDDO
  ENDDO
ENDIF
IF (L3MT) THEN
  IF (LCVPRO) THEN
    YDVARS%UNEBH%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=MAX(0._JPRB,YDVARS%UNEBH%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:))
    ZUNEBH(:,:)=MIN(1.0_JPRB-ZEPS0,YDVARS%UNEBH%T0(:,:)+MAX(0._JPRB,YDVARS%UAL%T0(:,:)))
  ELSE
    ZUNEBH(:,:)=YDVARS%UNEBH%T0(:,:)
  ENDIF
ENDIF

!     ------------------------------------------------------------------
!  The LMPHYS and LEPHYS keys should in the following be
! dispached to individual parametrizations.
IF(LMPHYS) THEN
!*
!     ------------------------------------------------------------------
!     4.- CALCULS THERMODYNAMIQUES
!     ----------------------------
  IF ( LTHERMO ) THEN
    CALL ACTQSAT (YDCST, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTQSAT, YDCPG_OPTS%KFLEVG,  &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, YDMF_PHYS_BASE_STATE%T, &
    & ZGEOSLC, ZMSC_LH, ZMSC_LSCPE, ZFLU_QSAT, ZMSC_QW, YDCPG_MISC%RH, ZMSC_TW)
  ENDIF

!*
!     ------------------------------------------------------------------
!     4.BIS. COEFFICIENTS THERMO-HYDRIQUES DU SOL
!     -------------------------------------------

  IF (LSFORCS) THEN  ! Surface forcing for 1D model MUSC
     DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZTSN(JLON)=YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)
        ZTN(JLON) =YDMF_PHYS_BASE_STATE%T(JLON,YDCPG_OPTS%KFLEVG)
     ENDDO
     DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZRHODREFM(JLON) = YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,YDCPG_OPTS%KFLEVG)/(YDMF_PHYS_BASE_STATE%T(JLON,YDCPG_OPTS%KFLEVG)*YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(JLON,YDCPG_OPTS%KFLEVG))
        ZTHETAS(JLON)   = ZTSN(JLON)*(RATM/YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(JLON,YDCPG_OPTS%KFLEVG))**RKAPPA
     ENDDO
     LLAROME = .FALSE.
     CALL SURF_IDEAL_FLUX(YDRIP, YDPHY0, YDPHYDS, LLAROME, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                           &
     & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(:, YDCPG_OPTS%KFLEVG), ZRHODREFM, YDMF_PHYS_SURF%GSD_SFO%PGROUP,                                &
     & ZTN, ZTSN, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_BASE_STATE%Q(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS_BASE_STATE%U(:, YDCPG_OPTS%KFLEVG), &
     & YDMF_PHYS_BASE_STATE%V(:, YDCPG_OPTS%KFLEVG), ZTHETAS, YDMF_PHYS%OUT%FCS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1),                    &
     & ZFEV, ZFMDU, ZFMDV)
     DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)=ZTSN(JLON)
     ENDDO
     DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)))
        ! To be equivalent to surfex forcing
        YDMF_PHYS%OUT%FEVL(JLON,1)=ZFEV(JLON)*(1.0_JPRB-ZDELTA)
        YDMF_PHYS%OUT%FEVN(JLON,1)=ZFEV(JLON)*ZDELTA
        YDMF_PHYS%OUT%FCLL(JLON,1)=YDMF_PHYS%OUT%FEVL(JLON,1)*ZMSC_LH(JLON,YDCPG_OPTS%KFLEVG)
        YDMF_PHYS%OUT%FCLN(JLON,1)=YDMF_PHYS%OUT%FEVN(JLON,1)*ZMSC_LH(JLON,YDCPG_OPTS%KFLEVG)
        ZDSA_LHS(JLON)=FOLH(YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON),0.0_JPRB)
     ENDDO
  ENDIF  ! End of surface forcing for 1D model MUSC

  IF ( .NOT.LMSE ) THEN
    IF ( LSOLV ) THEN
      LLHMT=.FALSE.
      CALL ACSOL (YDCPG_OPTS%YRCLI, YDCST, YDPHY, YDPHY1, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDMF_PHYS_SURF%GSD_VV%PARG,  &
      & YDMF_PHYS_SURF%GSD_VV%PD2, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF,                            &
      & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS_SURF%GSD_VV%PLAI, YDMF_PHYS_BASE_STATE%YGSP_SG%A,                        &
      & YDMF_PHYS_BASE_STATE%YGSP_SG%R, YDMF_PHYS_SURF%GSD_VV%PSAB, YDMF_PHYS_BASE_STATE%YGSP_SG%F, YDMF_PHYS_BASE_STATE%YGSP_RR%T,                 &
      & YDMF_PHYS_SURF%GSD_VF%PVEG, YDMF_PHYS_BASE_STATE%YGSP_SB%Q, YDMF_PHYS_BASE_STATE%YGSP_SB%TL,                                                &
      & YDMF_PHYS_BASE_STATE%YGSP_RR%W, YDMF_PHYS_BASE_STATE%YGSP_RR%IC, LLHMT, ZDSA_C1, ZDSA_C2,                                                   &
      & ZC3, ZCG, ZCN, YDMF_PHYS%OUT%CT, ZNEIJG, ZNEIJV, ZWFC, ZWPMX, ZWSEQ, ZWSMX, ZWWILT)
    ELSE

!            INITIALISATION DE L'INERTIE THERMIQUE DU SOL.

!DEC$ IVDEP
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDMF_PHYS%OUT%CT(JLON)=HSOL /(&
         & 1.0_JPRB+HSOLIWR*(YDMF_PHYS_BASE_STATE%YGSP_RR%W(JLON)+YDMF_PHYS_BASE_STATE%YGSP_SB%Q(JLON,1))/(WSMX+WPMX)&
         & *EXP(-0.5_JPRB*(HSOLIT0*(YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)-RTT))**2))
      ENDDO
    ENDIF
  ENDIF

!*
!     ------------------------------------------------------------------
!     4.TER.- INITIALISATIONS LIEES AU SCHEMA DE SURFACE EXTERNALISE
!     ------------------------------------------------------------------

  IF (LMSE) THEN

!     INITIALISATION DU SCHEMA DE SURFACE EXTERNALISE ET DES
!     VARIABLES PSEUDO-HISTORIQUES ASSOCIEES

    IF ( (NSWB_MNH /= NSW) .AND. LRAYFM ) CALL ABOR1('APLPAR: NSWB_MNH not = NSW')

    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%GZ0     (JLON) = YDCPG_GPAR%GZ0(JLON)
      YDMF_PHYS%OUT%GZ0H    (JLON) = YDCPG_GPAR%GZ0H(JLON)
      ZFLU_EMIS    (JLON) = YDCPG_GPAR%VEMIS(JLON)
      YDCPG_MISC%QS      (JLON) = YDCPG_GPAR%VQS(JLON)
      YDMF_PHYS_BASE_STATE%YGSP_RR%T      (JLON) = YDCPG_GPAR%VTS(JLON)
      ZSRAIN   (JLON) = YDCPG_GPAR%RAIN(JLON)
      ZSSNOW   (JLON) = YDCPG_GPAR%SNOW(JLON)
      ZSGROUPEL(JLON) = 0._JPRB
      ZTSN     (JLON) = YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)
    ENDDO

    IF ( LRAY ) THEN
      ! ACRANEB/ACRANEB2 radiation => one solar band
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDMF_PHYS%OUT%ALB   (JLON) = YDCPG_GPAR%ALBSCA(JLON,1)
        ZALBDIR(JLON) = YDCPG_GPAR%ALBDIR(JLON,1)
      ENDDO
    ELSE
      ! FMR  radiation => NSW solar bands
      DO JSG=1,NSW
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZALBP(JLON,JSG) = YDCPG_GPAR%ALBDIR(JLON,JSG)
          ZALBD(JLON,JSG) = YDCPG_GPAR%ALBSCA(JLON,JSG)
        ENDDO
      ENDDO
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDMF_PHYS%OUT%ALB(JLON)=0.0_JPRB
      ENDDO
      DO JSG=1,NSW
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS%OUT%ALB(JLON)+0.5*(ZALBP(JLON,JSG)+ZALBD(JLON,JSG))
        ENDDO
      ENDDO
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS%OUT%ALB(JLON)/REAL(NSW,JPRB)
      ENDDO
    ENDIF

    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      IF (ZFLU_EMIS(JLON)==0._JPRB) THEN
        ZFLU_EMIS(JLON) = 0.99_JPRB
        YDMF_PHYS_BASE_STATE%YGSP_RR%T  (JLON) = 288.0_JPRB
        YDMF_PHYS%OUT%ALB (JLON) = 0.1_JPRB
      ENDIF
    ENDDO

  ENDIF  ! LMSE
!
! Define z0;z0h if it's necessary
!
  IF (LMSE.AND.(.NOT.LCOEFKTKE).AND.(.NOT.LCOEFK_TOMS).AND.YDCPG_OPTS%NSTEP == 0) THEN
  CALL ARO_GROUND_DIAG_Z0( YDCPG_BNDS%KBL, YDCPG_OPTS%KGPCOMP, YDCPG_BNDS%KFDIA-YDCPG_BNDS%KIDIA+1, YDCPG_BNDS%KIDIA,                 &
  & YDCPG_BNDS%KFDIA, NDGUNG, NDGUXG, NDLUNG, NDLUXG, YDVARS%GEOMETRY%RINDX%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                    &
  & YDVARS%GEOMETRY%RINDY%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), LSURFEX_KFROM, YDMF_PHYS%OUT%GZ0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), &
  & YDMF_PHYS%OUT%GZ0H(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA))
  ENDIF
!*
!     ------------------------------------------------------------------
!     5.- STRUCTURE ET CHAMPS DANS LA COUCHE LIMITE DE SURFACE
!     ------------------------------------------------------------------

!        INITIALISATION DES HAUTEURS "METEO".

  IF ( LHMTO ) THEN
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZDPHIV(JLON)=RG*HVCLS
      ZDPHIT(JLON)=RG*HTCLS
    ENDDO
  ENDIF

  IF ( LVDIF.OR.LHMTO.OR.LGWD ) THEN
    LLCLS=LGWD.OR.LVDIF
    LLHMT=LHMTO
    IF (LMSE) THEN
      IF(LCOEFKTKE.AND.LCOEFKSURF) THEN

       IF (YDCPG_OPTS%NSTEP == 0) THEN 
         CALL ARO_GROUND_DIAG_Z0( YDCPG_BNDS%KBL, YDCPG_OPTS%KGPCOMP, YDCPG_BNDS%KFDIA-YDCPG_BNDS%KIDIA+1,                                   &
         & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, NDGUNG, NDGUXG, NDLUNG, NDLUXG, YDVARS%GEOMETRY%RINDX%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),  &
         & YDVARS%GEOMETRY%RINDY%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), LSURFEX_KFROM, YDMF_PHYS%OUT%GZ0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), &
         & YDMF_PHYS%OUT%GZ0H(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA))
        ENDIF


        CALL ACTKEZOTLS ( YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,             &
        & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
        & YDMF_PHYS_BASE_STATE%T, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, YDMF_PHYS_BASE_STATE%YGSP_RR%T,              &
        & YDCPG_MISC%QS, ZFLU_CDN, ZCDNH, ZDSA_CPS, ZRTI, ZDSA_RS)
         
      ELSE
        CALL ACHMTLS (YDCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,    &
        & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,               &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%T,                        &
        & YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDCPG_MISC%QS,                             &
        & ZDPHIT, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, LLCLS, ZNBVNO, ZMRIPP, ZDSA_CPS, ZGWDCS,                                    &
        & ZDSA_LHS, ZPCLS, ZFLU_CD, ZFLU_CDN)
!       Computation of ZRTI
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZDPHI(JLON)=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(JLON,YDCPG_OPTS%KFLEVG)-YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,YDCPG_OPTS%KFLEVG)
          ZPRS(JLON)=RD+ZRVMD*YDCPG_MISC%QS(JLON)
          ZRTI(JLON)=2.0_JPRB/(YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(JLON,YDCPG_OPTS%KFLEVG)*YDMF_PHYS_BASE_STATE%T(JLON,YDCPG_OPTS%KFLEVG)+RKAPPA*ZDPHI(JLON)&
           & +ZPRS(JLON)*YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON))
        ENDDO
      ENDIF
    ELSE
      IF (LCOEFKSURF) THEN
        CALL ACTKEHMT ( YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDSURF%YSD_VVD%NUMFLDS>=8.AND.LSOLV, &
        & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,               &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,     &
        & YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZPFL_FPLSH,                                        &
        & ZPFL_FPLCH, ZDPHIT, ZDPHIV, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF,          &
        & YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VF%PLSM, ZNEIJG, ZNEIJV, YDMF_PHYS_BASE_STATE%YGSP_SG%F,                     &
        & YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VF%PVEG, ZWFC, YDMF_PHYS_BASE_STATE%YGSP_RR%W,                          &
        & YDMF_PHYS_BASE_STATE%YGSP_RR%IC, LLCLS, LLHMT, ZNBVNO, ZMRIPP, ZFLU_CD, ZFLU_CDN, ZCDROV,                                  &
        & ZFLU_CH, ZCHROV, ZDSA_CPS, ZDQSTS, ZGWDCS, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, ZHQ, ZHU,                                &
        & ZFLU_NEIJ, YDMF_PHYS%OUT%QCLS, YDCPG_MISC%QS, ZFLU_QSATS, YDMF_PHYS%OUT%RHCLS, ZDSA_RS,                                    &
        & ZRTI, ZSTAB, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%NUCLS,                              &
        & YDMF_PHYS%OUT%NVCLS, ZPCLS, ZFLU_VEG, ZXDROV, ZXHROV, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST                               &
        &                                                                                                                            )
      ELSE
        CALL ACHMT (YDCPG_OPTS%YRCLI, YDMODEL%YRML_PHY_MF%YRPHY, YDMODEL%YRML_PHY_MF%YRPHY0, YDMODEL%YRML_PHY_MF%YRPHY1, YDMODEL%YRML_PHY_MF%YRPHY2, &
        & YDCST, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDSURF%YSD_VVD%NUMFLDS>=8.AND.LSOLV,                &
        & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                               &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,                     &
        & YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZPFL_FPLSH,                                                        &
        & ZPFL_FPLCH, ZDPHIT, ZDPHIV, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF,                          &
        & YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VF%PLSM, ZNEIJG, ZNEIJV, YDMF_PHYS_BASE_STATE%YGSP_SG%F,                                     &
        & YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VF%PVEG, ZWFC, YDMF_PHYS_BASE_STATE%YGSP_RR%W,                                          &
        & YDMF_PHYS_BASE_STATE%YGSP_RR%IC, LLCLS, LLHMT, ZNBVNO, ZMRIPP, ZFLU_CD, ZFLU_CDN, ZCDROV,                                                  &
        & ZFLU_CH, ZCHROV, ZDSA_CPS, ZDQSTS, ZGWDCS, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, ZHQ, ZHU,                                                &
        & ZFLU_NEIJ, YDMF_PHYS%OUT%QCLS, YDCPG_MISC%QS, ZFLU_QSATS, YDMF_PHYS%OUT%RHCLS, ZDSA_RS,                                                    &
        & ZRTI, ZSTAB, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%NUCLS,                                              &
        & YDMF_PHYS%OUT%NVCLS, ZPCLS, ZFLU_VEG, ZXDROV, ZXHROV, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST                                               &
        &                                                                                                                                            )
      ENDIF
    ENDIF

    IF (LPTKE) THEN
      YDMF_PHYS_BASE_STATE%TKE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG) = MAX(YDMF_PHYS_BASE_STATE%TKE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG),ETKE_MIN)
    ENDIF
    IF (LCOEFK_PTTE) THEN
      YDVARS%TTE%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG) = MAX(YDVARS%TTE%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG),ETKE_MIN)
    ENDIF

    IF(LCOEFKTKE) THEN
      ZCP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG) = RCPD*(1.0_JPRB+(RCPV/RCPD-1.0_JPRB)*(&
        & ZQV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)+ZQI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)+&
        & ZQL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)))
    ELSE
      ZCP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG) = YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
    ENDIF

    IF(LCOEFK_RIS .AND. LCOEFKTKE) THEN
      !  computation of Ri*,Ri** for mixing lenth computation
      CALL ACMRISS ( YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                             &
      & NTCOEF, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,                  &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, &
      & ZQV, ZQL, ZQI, ZFLU_QSAT, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U,     &
      & YDMF_PHYS_BASE_STATE%V, ZMSC_LSCPE, YDMF_PHYS%OUT%GZ0, ZMN2PP, ZMRIPP)
    ENDIF

    ! COMPUTATION OF mixing lengths from  Ri*,Ri** - FIRST GUES for moist AF

    !---------------------------------------------------
    ! COMPUTATION OF 'DRY' mixing lengths : lm_d lh_d
    ! COMPUTATION OF ZPBLH - PBL HEIGHT

    IF (CGMIXLEN == 'Z'  .OR. &
     &  CGMIXLEN == 'EL0'.OR. &
     &  CGMIXLEN == 'EL1'.OR. &
     &  CGMIXLEN == 'EL2'.OR. &
     &  CGMIXLEN == 'AY' .OR. &
     &  CGMIXLEN == 'AYC'.AND.(.NOT.LECT)) THEN
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZTHETAV(JLON,JLEV)=YDMF_PHYS_BASE_STATE%T(JLON,JLEV)*(1.0_JPRB+RETV*ZQV(JLON,JLEV))&
           & *(RATM/YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,JLEV))**RKAPPA  
        ENDDO
      ENDDO
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZTHETAVS(JLON)=YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)*(1.0_JPRB+RETV*YDCPG_MISC%QS(JLON))&
         & *(RATM/YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(JLON,YDCPG_OPTS%KFLEVG))**RKAPPA  
      ENDDO
      CALL ACCLPH (YDCST, YDPHY0, YDPHY2, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA, &
      & YDCPG_OPTS%KFLEVG, ZTHETAV, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,           &
      & YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZTHETAVS, IMOC_CLPH, YDMF_PHYS%OUT%CLPH, YDMF_PHYS%OUT%VEIN, &
      & ZUGST, ZVGST)
      ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      IF (.NOT.LRAFTUR) THEN
        YDMF_PHYS%OUT%UGST(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZUGST(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
        YDMF_PHYS%OUT%VGST(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZVGST(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      ENDIF
    ELSE
      ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VH%PPBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
    ENDIF

  ENDIF ! end of LVDIF or LHMTO or LGWD

  IF ( (LVDIF.OR.LGWD).AND.(.NOT.(LNEBR.OR.LECT)) ) THEN

    IF(TRIM(CGMIXLEN) == 'Z') THEN
      !-------------------------------------------------
      ! "z dependent" mixing length.
      !-------------------------------------------------
      ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=1.0_JPRB/UHDIFV
      CALL ACMIXLENZ ( YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, 1, YDCPG_OPTS%KFLEVG, &
      & .FALSE., YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZBLH, YDMF_PHYS%OUT%GZ0, &
      & YDMF_PHYS%OUT%GZ0H, ZLMU, ZLMT)

    ELSEIF((TRIM(CGMIXLEN) == 'TMC').OR.(TRIM(CGMIXLEN) == 'AYC')) THEN
      !     Cubique du climat
      ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)))
      CALL ACMIXLENTM ( YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,               &
      & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, &
      & ZBLH, ZLMU, ZLMT)

    ELSEIF(TRIM(CGMIXLEN) == 'TM') THEN
      !     Ancienne formulation pour Lm
      ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)))*XKLM
      CALL ACMIXLENZ ( YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, 1, YDCPG_OPTS%KFLEVG, &
      & .FALSE., YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZBLH, YDMF_PHYS%OUT%GZ0, &
      & YDMF_PHYS%OUT%GZ0H, ZLMU, ZLMT)

      !     Cubique du climat pour Lh
      ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)/XKLM
      CALL ACMIXLENTM ( YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,               &
      & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, &
      & ZBLH, ZZLMT, ZLMT                                                                                             &
      &                                                                                                                                                                                          )

    ELSEIF(TRIM(CGMIXLEN) == 'AY') THEN
      !     new Ayotte-Tudor ZBLH & mixing length
      CALL ACMIXLENZ ( YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, 1, YDCPG_OPTS%KFLEVG, &
      & .FALSE., YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZBLH, YDMF_PHYS%OUT%GZ0, &
      & YDMF_PHYS%OUT%GZ0H, ZLMU, ZLMT)

    ELSEIF((CGMIXLEN(1:2) == 'EL').AND.LPTKE) THEN
      !     e-type mixing length converted to Prandtl type
      CALL ACMIXLENZ(YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, 1, YDCPG_OPTS%KFLEVG,  &
      & .TRUE., YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZBLH, YDMF_PHYS%OUT%GZ0, &
      & YDMF_PHYS%OUT%GZ0H, ZLMU, ZLMT)
      ZLMU2(:,:)=ZLMU(:,:)
      ZLMT2(:,:)=ZLMT(:,:)
      IF     (CGMIXLEN == 'EL0') THEN
        IMLTYPE=0
        ! to have identical mixing length like in pTKE
        ! e-type mixing length converted to Prandtl type
        CALL ACMIXLENZ(YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, 1, YDCPG_OPTS%KFLEVG,   &
        & .FALSE., YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZBLH, YDMF_PHYS%OUT%GZ0, &
        & YDMF_PHYS%OUT%GZ0H, ZLMU, ZLMT)
        ZLMU2(:,:)=ZLMU(:,:)
        ZLMT2(:,:)=ZLMT(:,:)
      ELSEIF (CGMIXLEN == 'EL1') THEN
        IMLTYPE=1
      ELSEIF (CGMIXLEN == 'EL2') THEN
        IMLTYPE=2
      ELSE
        CLERR='APLPAR: UNEXPECTED VALUE FOR CGMIXLEN: '//TRIM(CGMIXLEN)
        CALL ABOR1(CLERR)
      ENDIF

      IF( LCOEFK_RIS) THEN
        LLMAF=.TRUE.
        CALL ACMIXELEN(YGFL, YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,                        &
        & YDCPG_OPTS%KFLEVG, YDCPG_OPTS%NSTEP, IMLTYPE, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,            &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%T,                            &
        & ZQV, ZQL, ZQI, YDMF_PHYS_BASE_STATE%R, YDMF_PHYS_BASE_STATE%S, YDMF_PHYS_BASE_STATE%TKE, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,   &
        & YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, &
        & ZMN2PP, ZFMGST, ZPFL_FPLSH, ZPFL_FPLCH, YDMF_PHYS%OUT%GZ0, YDMF_PHYS_BASE_STATE%YGSP_RR%T,                                      &
        & ZFLU_CDN, ZBLH, ZLMU, ZLMT, YDVARS%MXL%T0, ZLML, ZLMLTILD, ZRRCOR, LLMAF)
      ENDIF

    ELSE
      CLERR='APLPAR: UNEXPECTED VALUE FOR CGMIXLEN: '//TRIM(CGMIXLEN)
      CALL ABOR1(CLERR)
    ENDIF

    IF(LCOEFKTKE) THEN

      ! ------------------------------------------------------------- 
      ! COMPUTATION OF Ri', NCVPP AND COEFFICIENT FOR MOIST GUSTINESS
      ! ------------------------------------------------------------- 
      IF(LCOEFK_RIS) THEN
        CALL ACMRIS ( YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                              &
        & NTCOEF, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,                  &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZMSC_LSCPE, ZQV, ZQL, ZQI, ZFLU_QSAT, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
        & YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZLMU, ZLMT, YDMF_PHYS%OUT%GZ0,             &
        & ZMRIPP)
      ENDIF

      CALL ACMRIP(YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTCOEF, YDCPG_OPTS%KFLEVG,                 &
      & YDCPG_OPTS%NSTEP, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, ZQV, ZQL, ZQI, ZCP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH,                             &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZFLU_QSAT, ZMSC_QW, ZMSC_TW, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,                      &
      & YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDVARS%FQTUR%T0, YDVARS%FSTUR%T0,                      &
      & YDVARS%SHTUR%T0, YDMF_PHYS_BASE_STATE%TKE, YDVARS%TTE%T0, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VH%PQSH,          &
      & ZDSA_RS, ZDSA_CPS, ZRTI, YDMF_PHYS%OUT%GZ0, LLCLS, ZMRIPP, ZMRIFPP, ZBNEBCVPP, ZBNEBQ,                                         &
      & ZNBVNO, ZFMTKE, ZFTTKE, ZF_EPS, ZFUN_TTE, ZAUTKE, ZATTKE, ZFHORM, ZFHORH, ZTH_FUN, ZWW_FUN,                                    &
      & ZMRIMC, ZMRICTERM, ZMN2PP, ZMN2_ES, ZMN2_EQ, ZMN2_DS, ZMN2_DQ, ZFMGST)

    ENDIF ! LCOEFTKE

    ! FINISHING MIXING LENGTH COMPUTATION
    IF((CGMIXLEN(1:2) == 'EL').AND.LPTKE) THEN
      ZLMU(:,:)=ZLMU2(:,:)
      ZLMT(:,:)=ZLMT2(:,:)
      IF     (CGMIXLEN == 'EL0') THEN
        IMLTYPE=0
      ELSEIF (CGMIXLEN == 'EL1') THEN
        IMLTYPE=1
      ELSEIF (CGMIXLEN == 'EL2') THEN
        IMLTYPE=2
      ENDIF
      LLMAF=.FALSE.
      CALL ACMIXELEN(YGFL, YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,                        &
      & YDCPG_OPTS%KFLEVG, YDCPG_OPTS%NSTEP, IMLTYPE, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,            &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%T,                            &
      & ZQV, ZQL, ZQI, YDMF_PHYS_BASE_STATE%R, YDMF_PHYS_BASE_STATE%S, YDMF_PHYS_BASE_STATE%TKE, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,   &
      & YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, &
      & ZMN2PP, ZFMGST, ZPFL_FPLSH, ZPFL_FPLCH, YDMF_PHYS%OUT%GZ0, YDMF_PHYS_BASE_STATE%YGSP_RR%T,                                      &
      & ZFLU_CDN, ZBLH, ZLMU, ZLMT, YDVARS%MXL%T0, ZLML, ZLMLTILD, ZRRCOR, LLMAF)
    ENDIF

  ENDIF ! (LVDIF or LGWD) and( not(LNEBR or LECT))

  IF ( LVDIF.OR.LHMTO.OR.LGWD ) THEN

    IF (LFLUSO.AND.(.NOT.LMSE)) THEN
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZCRTI(JLON) = 1.0_JPRB/(YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)*ZDSA_RS(JLON))
      ENDDO
      CALL ACFLUSO (YDCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,      &
      & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,   &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%Q, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, &
      & YDMF_PHYS_BASE_STATE%V, ZDPHIT, ZDPHIV, YDMF_PHYS%OUT%GZ0, YDMF_PHYS_SURF%GSD_VF%PLSM,                         &
      & ZFLU_QSATS, ZCRTI, YDMF_PHYS_BASE_STATE%YGSP_RR%T, LLHMT, ZFLU_CD, ZFLU_CDN, ZCDROV, ZCE,                      &
      & ZCEROV, ZFLU_CH, ZCHROV, YDMF_PHYS%OUT%QCLS, YDMF_PHYS%OUT%RHCLS, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%UCLS,      &
      & YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST)
    ELSE
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZCE   (JLON) = ZFLU_CH   (JLON)
        ZCEROV(JLON) = ZCHROV(JLON)
      ENDDO
    ENDIF

  ENDIF ! LVDIF or LHMTO or LGWD

  IF (LRAFTKE) THEN
    YDMF_PHYS%OUT%CAPE(:)=0._JPRB
    ZDCAPE(:)=0._JPRB
    CALL ACCLDIA(YDCST, YDCPG_OPTS%LXCLP, YDCPG_OPTS%LXTGST, YDCPG_OPTS%LXXGST, YDPHY, YDPHY2, YDTOPH, YDCPG_BNDS%KIDIA,      &
    & YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS_BASE_STATE%U,           &
    & YDMF_PHYS_BASE_STATE%V, YDMF_PHYS%OUT%CAPE, ZDCAPE, YDMF_PHYS_BASE_STATE%TKE, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,               &
    & YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST, ZBLH, IMOC_CLPH)
  ENDIF

!**
!     ------------------------------------------------------------------
!     6.- TURBULENCE: COEFFICIENTS D'ECHANGE
!     ------------------------------------------------------------------

  IF ( (LVDIF.OR.LGWD).AND.(.NOT.(LNEBR.OR.LECT)) ) THEN
     !-------------------------------------------------
     ! Compute diffusion coefficients.
     !-------------------------------------------------

    IF(LCOEFKTKE) THEN
      CALL ACTKECOEFK ( YDPHY0, YDPHY2, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTCOEF,                                    &
      & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
      & ZFMTKE, ZFTTKE, ZAUTKE, ZATTKE, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%T,                                    &
      & YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZLMU, ZLMT, YDMF_PHYS%OUT%GZ0, ZKTROV,                                          &
      & ZKUROV, ZKNROV, ZXTROV, ZXUROV, ZXPTKEROV)
    ELSE
      CALL ACCOEFK ( YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                             &
      & NTCOEF, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,                  &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZMSC_LSCPE, ZQV, ZQL, ZQI, ZFLU_QSAT, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
      & YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZPFL_FPLSH,                                &
      & ZPFL_FPLCH, ZLMU, ZLMT, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, ZBLH, ZKTROV, ZKUROV, ZKNROV,                       &
      & ZNBVNO, ZXTROV, ZXUROV, ZXPTKEROV)
    ENDIF
  ENDIF

  IF (LCVPPKF) THEN
    CALL ACVPPKF(YDCST,YDMODEL%YRML_PHY_MF, YDCPG_BNDS, YDCPG_OPTS, NTCVIM,                                               &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,  &
    & YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%T, ZQV, ZQL, ZQI, YDMF_PHYS_BASE_STATE%U,                 &
    & YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%TKE,                             &
    & ZDIFCVPPQ, ZDIFCVPPS, ZCONDCVPPL, ZCONDCVPPI, ZPRODTH_CVPP, INLAB_CVPP, ZQLI_CVPP, ZNEB_CVPP,                       &
    & INND)
  ENDIF

     !-------------------------------------------------
     !  Call to EDKF
     !-------------------------------------------------
  IF (LEDKF) THEN
    IF (LEDMFI) THEN
       ZIMPL=0._JPRB
    ELSE  
       ZIMPL=1._JPRB 
    ENDIF
    CALL ABOR1('APLPAR: CODE MUST BE UPDATED, IMPL_MF IS NOW SET IN NAMELIST')
      
    IF (YDCPG_OPTS%NSTEP == 0) YDMF_PHYS_SURF%GSD_SFL%PGROUP(:,:) = 0.0_JPRB
    CALL ARP_SHALLOW_MF( YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG,              &
    & TSPHY, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,        &
    & YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, CMF_UPDRAFT, CMF_CLOUD, LMIXUV, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V,     &
    & YDMF_PHYS_BASE_STATE%T, ZQV, ZQL, ZQI, ZQR, ZQS, YDMF_PHYS_BASE_STATE%TKE, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,         &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZEDMFQ, ZEDMFS, ZEDMFU, ZEDMFV, YDMF_PHYS_SURF%GSD_SFL%PGROUP(:, 1),              &
    & YDMF_PHYS_SURF%GSD_SFL%PGROUP(:, 2), ZPRODTH_CVPP, ZQLI_CVPP, ZNEB_CVPP, INLAB_CVPP, ZMF_UP)

  ENDIF

  IF ( LVDIF.AND.LECT ) THEN
    IF ( LCONDWT ) THEN
       YDCPG_MISC%QICE(:,:)= ZQI(:,:)
       YDCPG_MISC%QLI(:,:) = ZQL(:,:)
    ELSE
       YDCPG_MISC%QICE(:,:)= 0.0_JPRB
       YDCPG_MISC%QLI(:,:) = 0.0_JPRB
    ENDIF

! Computation of the 2 znlab used in acbl89
    IF (.NOT. LECSHAL) INLAB_CVPP(:,:) = 0
    IF (LECDEEP) THEN
       ZNLABCVP(:,:) = 1.0_JPRB
    ELSE
       ZNLABCVP(:,:) = 0.0_JPRB
    ENDIF
    IF(LNEBN.OR.LNEBR.OR.LRRGUST) THEN
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZNLABCVP(JLON,JLEV) = ZNLABCVP(JLON,JLEV)&
         & *MAX(0.0_JPRB,SIGN(1.0_JPRB,ZPFL_FPLCH(JLON,JLEV)-ZPFL_FPLCH(JLON,JLEV-1)-YDPHY0%REPS))
          ZNLAB(JLON,JLEV) = REAL(INLAB_CVPP(JLON,JLEV),JPRB)
        ENDDO
      ENDDO
    ENDIF
    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZNLAB(JLON,JLEV) = REAL(INLAB_CVPP(JLON,JLEV),JPRB)
      ENDDO
    ENDDO

    CALL ACTKE (YDCST, YDLDDH, YDMODEL%YRML_DIAG%YRMDDH, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,        &
    & YDCPG_OPTS%KLON, NTCOEF, NTCOET, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,       &
    & YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V,               &
    & ZQV, ZQIC, ZQLC, ZMSC_LSCPE, ZFLU_CD, ZFLU_CH, YDMF_PHYS%OUT%GZ0, YDMF_PHYS_BASE_STATE%YGSP_RR%T,                          &
    & YDCPG_MISC%QS, YDCPG_MISC%QICE, YDCPG_MISC%QLI, YDMF_PHYS_BASE_STATE%TKE, ZPRODTH_CVPP, ZNLAB,                             &
    & ZNLABCVP, ZKTROV, ZKQROV, ZKQLROV, ZKUROV, ZXTROV, ZXUROV, ZNBVNO, ZNEBS, ZQLIS, ZNEBS0,                                   &
    & ZQLIS0, ZCOEFN, ZPFL_FTKE, ZPFL_FTKEI, ZTKE1, ZTPRDY, YDMF_PHYS%OUT%EDR, YDDDH)
    YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=MIN(XMAXLM,MAX(XMINLM,YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)))
  ENDIF


     !-------------------------------------------------
     ! Store diffusion coefficients in trajectory in temporary variables
     ! before final writing.
     !-------------------------------------------------
  ZKTROV_SAVE(:,:)=ZKTROV(:,:)
  ZKUROV_SAVE(:,:)=ZKUROV(:,:)
  ZCDROV_SAVE(:)=ZCDROV(:)
  ZCHROV_SAVE(:)=ZCHROV(:)

!**
!     ------------------------------------------------------------------
!     7.- RAYONNEMENT
!     ----------------
!     --------------------------------------------------------------------
!      - COMPUTE DUST PROPERTIES FOR RADIATION IF LMDUST=T
!     --------------------------------------------------------------------
  IF (LMDUST.AND.(NGFL_EXT/=0)) THEN
! input dust scalar concentration in ppp from

  CALL ARO_MNHDUST (1, ILONMNH, YDCPG_OPTS%KFLEVG, NGFL_EXT, YDCPG_OPTS%ZDTPHY, ZZI_SVM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :, 1:NGFL_EXT),                   &
  & ZZZ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :), ZDZZ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :), ZZI_PABSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :), &
  & ZZI_THM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :), ZZI_RHODREFM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :),                                         &
  & NSWB_MNH, YDCPG_OPTS%NSTEP+1, ZZI_SVM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :, 1:NGFL_EXT), ZPIZA_DST(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :),  &
  & ZCGA_DST(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :), ZTAUREL_DST(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :),                                         &
  & ZAERD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :), NGFL_EZDIAG, ZZI_PEZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :)                                   &
  &                                                                                                                                                                                                                              )

  ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,:)=ZZI_PEZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,:)

! return to aladin environment (inversion des niveaux)
   DO JGFL=1,NGFL_EXT
     DO JLON=1,YDCPG_OPTS%KLON
       DO JLEV=1,YDCPG_OPTS%KFLEVG
         ZSVM(JLON,JLEV,JGFL)=ZZI_SVM(JLON,1,JLEV,JGFL)
       ENDDO
     ENDDO
   ENDDO
  ENDIF

!      7.1 Albedo et emissivite en presence de neige
!          Albedo and emissivity with snow

  IF (.NOT.LMSE) THEN
!DEC$ IVDEP
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      IF (LSNV) THEN
        IF ((YDMF_PHYS_SURF%GSD_VF%PVEG(JLON) < 0.01_JPRB).OR.(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON) >= 0.60_JPRB)) THEN
          ZALBV=0.0_JPRB
          YDMF_PHYS_SURF%GSD_VF%PALBSF(JLON)=YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)
        ELSE
          ZALBV=(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)-(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*YDMF_PHYS_SURF%GSD_VF%PALBSF(JLON))/YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)
        ENDIF
        YDMF_PHYS%OUT%ALB(JLON)= (1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*(1.0_JPRB-ZNEIJG(JLON)) *&
         & YDMF_PHYS_SURF%GSD_VF%PALBSF(JLON)&
         & + (1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*ZNEIJG(JLON) *&
         & MAX(YDMF_PHYS_SURF%GSD_VF%PALBSF(JLON),YDMF_PHYS_BASE_STATE%YGSP_SG%A(JLON,1))&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*ZNEIJV(JLON) *&
         & MAX(ZALBV,YDMF_PHYS_BASE_STATE%YGSP_SG%A(JLON,1))&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*(1.0_JPRB-ZNEIJV(JLON)) * ZALBV
        ZFLU_EMIS(JLON)= (1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*(1.0_JPRB-ZNEIJG(JLON)) *&
         & YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)&
         & + (1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*ZNEIJG(JLON) * EMCRIN&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*ZNEIJV(JLON) * EMCRIN&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*(1.0_JPRB-ZNEIJV(JLON)) * YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)
      ELSE
        IF (LVGSN) THEN
          IF (LZ0HSREL.AND.LCOEFKSURF) THEN
            ! new treatment, PNEIJ is gridbox snow fraction
            YDMF_PHYS%OUT%ALB(JLON)=(1.0_JPRB-ZFLU_VEG(JLON)-ZFLU_NEIJ(JLON))*YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)+ &
             & ZFLU_VEG(JLON)*YDMF_PHYS_SURF%GSD_VV%PALV(JLON)+ZFLU_NEIJ(JLON)*YDMF_PHYS_BASE_STATE%YGSP_SG%A(JLON,1)
          ELSE
            ! old treatment, PNEIJ is snow fraction for bare ground
            YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)-ZFLU_NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)- &
             & YDMF_PHYS_BASE_STATE%YGSP_SG%A(JLON,1))+(ZFLU_NEIJ(JLON)-ZNEIJV(JLON))*     &
             & YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*(YDMF_PHYS_SURF%GSD_VV%PALV(JLON)-YDMF_PHYS_BASE_STATE%YGSP_SG%A(JLON,1))
          ENDIF

          YDMF_PHYS%OUT%ALB(JLON)=MIN(ABS(YDMF_PHYS_SURF%GSD_VV%PIVEG(JLON)-2._JPRB),1.0_JPRB) * YDMF_PHYS%OUT%ALB(JLON) +(&
           & 1.0_JPRB-MIN(ABS(YDMF_PHYS_SURF%GSD_VV%PIVEG(JLON)-2._JPRB),1.0_JPRB))&
           & * MAX(ALCRIN,YDMF_PHYS%OUT%ALB(JLON))
          YDMF_PHYS_SURF%GSP_SG%PT_T1(JLON,1)=YDMF_PHYS%OUT%ALB(JLON)

          ZFLU_EMIS(JLON)=YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-ZFLU_NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-EMCRIN)

        ELSE
          YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)-ZFLU_NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)&
           & -MAX(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON),ALCRIN))
          ZFLU_EMIS(JLON)=YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-ZFLU_NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-EMCRIN)
        ENDIF
      ENDIF
    ENDDO

    IF (LRAYFM) THEN
      ! diffuse and direct (parallel) albedo in NSW solar intervals
      IF (LALBMERCLIM) THEN
        DO JSG=1,NSW
          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            ZALBD(JLON,JSG)=YDMF_PHYS%OUT%ALB(JLON)
            ZALBPMER=(1.0_JPRB+&
             & 0.5_JPRB*ZRDG_MU0M(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))/   (&
             & 1.0_JPRB+ZRDG_MU0M(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))**2
            ZALBP(JLON,JSG)=YDMF_PHYS%OUT%ALB(JLON)*          YDMF_PHYS_SURF%GSD_VF%PLSM(JLON)+&
                          & ZALBPMER  *(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PLSM(JLON))
          ENDDO
        ENDDO
      ELSE
!DEC$ IVDEP
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZMERL(JLON)=(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PLSM(JLON))*(1.0_JPRB&
            & -MAX(0.0_JPRB,SIGN(1.0_JPRB,TMERGL-YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON))))
          ZFLU_EMIS(JLON)=ZFLU_EMIS(JLON)*YDMF_PHYS_SURF%GSD_VF%PLSM(JLON)+ZMERL(JLON)*EMMMER&
            & +(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PLSM(JLON))*(1.0_JPRB-ZMERL(JLON))*EMMGLA
        ENDDO
        DO JSG=1,NSW
          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            ZALBD(JLON,JSG)=(1.0_JPRB-ZMERL(JLON))*YDMF_PHYS%OUT%ALB(JLON)+&
                                    & ZMERL(JLON) *ALBMED
            ZALBP(JLON,JSG)=(1.0_JPRB-ZMERL(JLON))*YDMF_PHYS%OUT%ALB(JLON)+&
             & ZMERL(JLON) *&
             & MAX(0.037_JPRB/(1.1_JPRB*ZRDG_MU0(JLON)**1.4_JPRB+0.15_JPRB),ZEPS0)
          ENDDO
        ENDDO
      ENDIF
    ELSEIF (LRAY) THEN
      ! direct (parallel) albedo for ACRANEB/ACRANEB2, Geleyn's formula
      ! with given proportion of Lambertian scattering
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        IF ( YDMF_PHYS_SURF%GSD_VF%PLSM(JLON) < 0.5_JPRB .AND. YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON) >= TMERGL ) THEN
          ZLAMB=RLAMB_WATER  ! water surface (open sea)
        ELSE
          ZLAMB=RLAMB_SOLID  ! solid surface (frozen sea or land)
        ENDIF
        ZALBDIR(JLON)=ZLAMB*YDMF_PHYS%OUT%ALB(JLON)+(1._JPRB-ZLAMB)*(1._JPRB+&
         & 0.5_JPRB*ZRDG_MU0(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))/   (&
         & 1.0_JPRB+ZRDG_MU0(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))**2
      ENDDO
    ENDIF

  ENDIF  ! .NOT.LMSE

  ! Appel de la routine d'aerosols

  LLAERO=LAEROSEA.AND.LAEROLAN.AND.LAEROSOO.AND.LAERODES

  IF    (   (LRAYFM.AND.(MOD(YDCPG_OPTS%NSTEP,NRADFR) == 0)) &
  & .OR.  ( (LRAY.OR.LRAYSP).AND.(.NOT.LRSTAER)) ) THEN

    IF (LLAERO) THEN
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZAESEA(JLON) = YDMF_PHYS_SURF%GSD_VA%PSEA(JLON)
        ZAELAN(JLON) = YDMF_PHYS_SURF%GSD_VA%PLAN(JLON)
        ZAESOO(JLON) = YDMF_PHYS_SURF%GSD_VA%PSOO(JLON)
        ZAEDES(JLON) = YDMF_PHYS_SURF%GSD_VA%PDES(JLON)
      ENDDO
    ELSE
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZAESEA(JLON) = 0.0_JPRB
        ZAELAN(JLON) = 0.0_JPRB
        ZAESOO(JLON) = 0.0_JPRB
        ZAEDES(JLON) = 0.0_JPRB
      ENDDO
    ENDIF
    IF (LAEROSUL) THEN
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZAESUL(JLON) = YDMF_PHYS_SURF%GSD_VA%PSUL(JLON)
      ENDDO
    ELSE
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZAESUL(JLON) = 0.0_JPRB
      ENDDO
    ENDIF
    IF (LAEROVOL) THEN
       DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZAEVOL(JLON) = YDMF_PHYS_SURF%GSD_VA%PVOL(JLON)
      ENDDO
    ELSE
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZAEVOL(JLON) = 0.0_JPRB
      ENDDO
    ENDIF

    IF ( ( (LRAYFM.AND.NAER /= 0) .OR.LRAY.OR.LRAYSP).AND.LLAERO )  THEN
      CALL RADAER ( YDMODEL%YRML_PHY_RAD%YREAERD, YDERAD, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,                     &
      & YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
      & YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%YGSP_RR%T, ZAESEA, ZAELAN, ZAESOO, ZAEDES,                          &
      & ZAESUL, ZAEVOL, ZAER, ZAERINDS                                                     )
    ENDIF

  ELSEIF ( (LRAY.OR.LRAYSP).AND.LRSTAER ) THEN

    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZAER(JLON,JLEV,1)=ZDAER(JLEV)
        ZAER(JLON,JLEV,2:6)=0._JPRB
      ENDDO
    ENDDO

  ENDIF ! FOR AEROSOLS

! We uses the extinction coefficient explicitely solved by ARO_MNHDUST
  IF (LMDUST.AND.(NGFL_EXT/=0)) THEN
         ZAER(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:,3) = ZAERD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
  ENDIF


!      7.2 Flux radiatifs par ciel clair (Code Geleyn)
!          Clear sky radiative fluxes    (Geleyn's scheme)

  ! separate clearsky call is kept only for old ACRANEB; for ACRANEB2
  ! duplicit calculation of gaseous transmissions is avoided
  IF (LRAY.AND.NRAY == 1.AND.YDCFU%NFRRC /= 0) THEN
    IF (MOD(YDCPG_OPTS%NSTEP,YDCFU%NFRRC) == 0) THEN
      CALL ACRANEB(YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                       &
      & NTRADI, YDCPG_OPTS%KFLEVG, IJN, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,         &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%NEB, ZQV, ZQCO2, YDCPG_MISC%QICE, YDCPG_MISC%QLI,              &
      & ZQO3, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS%OUT%ALB, ZALBDIR, ZFLU_EMIS, ZRDG_MU0, YDVARS%GEOMETRY%GEMU%T0,           &
      & YDVARS%GEOMETRY%GELAM%T0, ZRDG_MU0LU, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH,     &
      & ZFRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOLU, YDMF_PHYS%OUT%FRTHDS, ZAER, ZMAK, ZMAN)
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDMF_PHYS%OUT%FRSOC(JLON,0)=YDMF_PHYS%OUT%FRSO(JLON,NTRADI-1,1)
        YDMF_PHYS%OUT%FRSOC(JLON,1)=YDMF_PHYS%OUT%FRSO(JLON,YDCPG_OPTS%KFLEVG,1)
        YDMF_PHYS%OUT%FRTHC(JLON,0)=YDMF_PHYS%OUT%FRTH(JLON,NTRADI-1,1)
        YDMF_PHYS%OUT%FRTHC(JLON,1)=YDMF_PHYS%OUT%FRTH(JLON,YDCPG_OPTS%KFLEVG,1)
      ENDDO
    ENDIF
  ENDIF

!      7.3 Nebulosite et Convection
!          Cloud cover and Convection
!      7.3.1 Shallow + Deep convection

  IF (LCVPGY) THEN
    ! Le schema de convection de J. F. Gueremy
    IF (LCONDWT) THEN
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZQCL(JLON,JLEV)=ZQL(JLON,JLEV)
          ZQCI(JLON,JLEV)=ZQI(JLON,JLEV)
        ENDDO
      ENDDO
    ELSE
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZQCL(JLON,JLEV)=0.0_JPRB
          ZQCI(JLON,JLEV)=0.0_JPRB
        ENDDO
      ENDDO
    ENDIF
    CALL ACCVIMPGY ( YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                             &
    & NTCVIM, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,               &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP,   &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZMSC_LH, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZQV,                        &
    & ZQCI, ZQCL, ZQLIS, ZFLU_QSAT, ZMSC_QW, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, &
    & YDMF_PHYS_BASE_STATE%T, ZMSC_TW, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZDSA_CPS, YDVARS%GEOMETRY%GM%T0,    &
    & YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCQN,                     &
    & YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,           &
    & YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, ZFHMLTSC, ZFHEVPPC, ZFPEVPPC, YDMF_PHYS%OUT%FPLCL,                     &
    & YDMF_PHYS%OUT%FPLCN, ZNEBC, ZQLIC, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, ICIS, INLAB, INND,                      &
    & YDMF_PHYS_BASE_STATE%CVV)

    DO JLEV = YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
        YDMF_PHYS%OUT%DIFCS(JLON,JLEV) = YDMF_PHYS%OUT%DIFCS(JLON,JLEV) - ZFHEVPPC(JLON,JLEV)&
                                            & - ZFHMLTSC(JLON,JLEV)
        YDMF_PHYS%OUT%DIFCQ(JLON,JLEV) = YDMF_PHYS%OUT%DIFCQ(JLON,JLEV) + YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV)&
                                            & + YDMF_PHYS%OUT%FPEVPCN(JLON,JLEV)
      ENDDO
    ENDDO
    YDMF_PHYS%OUT%FPEVPCL=0._JPRB
    YDMF_PHYS%OUT%FPEVPCN=0._JPRB
    IF (LGRAPRO) THEN
      YDMF_PHYS%OUT%FPEVPCG=0._JPRB
    ENDIF
    ! Prise en compte des nuages convectifs diagnostiques sortant d'ACMTUD
    IF (LNCVPGY) THEN
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZNEBC0(JLON,JLEV)=ZNEBC(JLON,JLEV)
          ZQLI_CVP(JLON,JLEV)=ZQLIC(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

    ! Annulation possible des flux convectifs pour les eaux condensees.
    IF (.FALSE.) THEN
      DO JLEV=0,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDMF_PHYS%OUT%DIFCQL(JLON,JLEV)=0.0_JPRB
          YDMF_PHYS%OUT%DIFCQN(JLON,JLEV)=0.0_JPRB
        ENDDO
      ENDDO
    ENDIF
  ENDIF !  (LCVPGY)

  IF ( LCONDWT.AND.(.NOT.LNEBECT)) THEN

    IF(LCVPRO.AND.LNEBCV) THEN
! convective cloudiness in case we need protection of convective cloud water.
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZNEBCH(JLON,JLEV)=ZUNEBH(JLON,JLEV)
        ENDDO
      ENDDO
      IF (LCVCSD) THEN
        DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            ZNEBC0(JLON,JLEV)=MAX(ZEPSNEB,&
                            & MIN(1._JPRB-ZEPSNEB,ZUNEBH(JLON,JLEV)))
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    CALL ACNEBCOND (YDCST, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,            &
    & NTPLUI, YDCPG_OPTS%KFLEVG, LLREDPR, YDSTA, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,        &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,        &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%RH, ZBLH, ZQV, ZQI, ZQL, ZMSC_QW, YDMF_PHYS_BASE_STATE%T,             &
    & ZNEBCH, YDVARS%GEOMETRY%GM%T0, YDMF_PHYS_BASE_STATE%YGSP_RR%T, ZQLIS, ZNEBS, ZRHCRI, ZRH,                                &
    & ZQSATS, ZICEFR1, ZQLIS0, ZNEBS0)
     
    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'ZRHCRI', ZRHCRI, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG) 
    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'ZQLIS0', ZQLIS0, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'ZQLIS', ZQLIS, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

    IF(LRKCDEV) THEN
! Rash-Kristiansson cloud water scheme - second part.
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ! analytical solution dRH/dCLOUD:                                                
          ZRHDFDA(JLON,JLEV)=2._JPRB*(1._JPRB - ZNEBS(JLON,JLEV))&
               & *(1.0_JPRB-ZRHCRI(JLON,JLEV))
        ENDDO
      ENDDO
    ENDIF


  ENDIF ! LCONDWT .AND. .NOT.LNEBECT

  !-------------------------------------------------
  ! PCMT convection scheme.
  !-------------------------------------------------
  IF(LGPCMT) THEN
    ZSMOOTRAC(1:INBTRA) = GCVTSMO ! as usual 

    IF(LEDMFI) THEN
      CALL ACPCMT(YDCST, YDGEM, YDGEOMETRY%YRDIM, YDGEOMETRY%YREGEO, YDLDDH, YDMODEL%YRML_DIAG%YRMDDH,                      &
      & YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTCVIM, YDCPG_OPTS%KFLEVG,                 &
      & INBTRA, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,      &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDCPG_DYN0%CTY%VVEL(:, 1:),                    &
      & YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZMSC_LH, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, YDMF_PHYS_BASE_STATE%Q, ZQI, ZQL, YDMF_PHYS_BASE_STATE%R,                           &
      & YDMF_PHYS_BASE_STATE%S, ZQLIS, ZFLU_QSAT, ZMSC_QW, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                                 &
      & ZRDG_CVGQ, YDMF_PHYS_BASE_STATE%T, ZMSC_TW, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%TKE,       &
      & ZTRA(:, :, INBTRA_DEP:INBTRA), ZSMOOTRAC(1:INBTRA), ZQLC, ZQIC, ZQRC, ZQSC, YDVARS%GEOMETRY%GM%T0,                          &
      & YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDCPG_MISC%QS, ZFLU_NEIJ, YDMF_PHYS_SURF%GSD_VF%PLSM, YDSTA,                                &
      & ZEDMFQ, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC,                           &
      & ZEDMFS, ZDIFCVTH, ZTMPPRODTH, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, ZEDMFU, ZEDMFV, ZSTRCTRA(:, :, INBTRA_DEP:INBTRA),  &
      & YDMF_PHYS%OUT%FIMCC, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN,                &
      & ZMF_UP, ZMU, ZMD, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC,                   &
      & YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC, YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC,         &
      & YDMF_PHYS%OUT%FCNEGQSC, INLAB, ZNEBC0, ZQLI_CVP, ZTU, ZQU, ZQC_DET_PCMT, ZCSGC, ZENTCH, INND,                               &
      & YDMF_PHYS%OUT%CAPE, ZAIPCMT, ZALF_CAPE, ZALF_CVGQ, YDVARS%UAL%T0, YDVARS%UOM%T0, YDVARS%DAL%T0,                             &
      & YDVARS%DOM%T0, YDDDH)
    ELSE
      CALL ACPCMT(YDCST, YDGEM, YDGEOMETRY%YRDIM, YDGEOMETRY%YREGEO, YDLDDH, YDMODEL%YRML_DIAG%YRMDDH,                      &
      & YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTCVIM, YDCPG_OPTS%KFLEVG,                 &
      & INBTRA, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,      &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDCPG_DYN0%CTY%VVEL(:, 1:),                    &
      & YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZMSC_LH, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, YDMF_PHYS_BASE_STATE%Q, ZQI, ZQL, YDMF_PHYS_BASE_STATE%R,                           &
      & YDMF_PHYS_BASE_STATE%S, ZQLIS, ZFLU_QSAT, ZMSC_QW, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                                 &
      & ZRDG_CVGQ, YDMF_PHYS_BASE_STATE%T, ZMSC_TW, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%TKE,       &
      & ZTRA(:, :, INBTRA_DEP:INBTRA), ZSMOOTRAC(1:INBTRA), ZQLC, ZQIC, ZQRC, ZQSC, YDVARS%GEOMETRY%GM%T0,                          &
      & YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDCPG_MISC%QS, ZFLU_NEIJ, YDMF_PHYS_SURF%GSD_VF%PLSM, YDSTA,                                &
      & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC,              &
      & YDMF_PHYS%OUT%DIFCS, ZDIFCVTH, ZTMPPRODTH, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%STRCU,                   &
      & YDMF_PHYS%OUT%STRCV, ZSTRCTRA(:, :, INBTRA_DEP:INBTRA), YDMF_PHYS%OUT%FIMCC, YDMF_PHYS%OUT%FPEVPCL,                         &
      & YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, ZMF_UP, ZMU, ZMD, YDMF_PHYS%OUT%FPFPCL,                    &
      & YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC,               &
      & YDMF_PHYS%OUT%FCNEGQLC, YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC, YDMF_PHYS%OUT%FCNEGQSC,                             &
      & INLAB, ZNEBC0, ZQLI_CVP, ZTU, ZQU, ZQC_DET_PCMT, ZCSGC, ZENTCH, INND, YDMF_PHYS%OUT%CAPE,                                   &
      & ZAIPCMT, ZALF_CAPE, ZALF_CVGQ, YDVARS%UAL%T0, YDVARS%UOM%T0, YDVARS%DAL%T0, YDVARS%DOM%T0,                                  &
      & YDDDH)
    ENDIF
    IF(LFLASH) THEN
      ! Lightning flashes: interface between PCMT and CULIGHT input data.
      CALL ACLIGHT(YDMODEL%YRML_PHY_MF,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,YDCPG_OPTS%KLON,NTCVIM,YDCPG_OPTS%KFLEVG,&
                   & ZQLC,ZQIC,ZLU,ICBOT_LIG,ICTOP_LIG,LLCUM_LIG,YDMF_PHYS%OUT%CAPE)
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZMFU(JLON,JLEV)=ZMF_UP(JLON,JLEV)
        ENDDO
      ENDDO
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        LLLAND(JLON)=YDMF_PHYS_SURF%GSD_VF%PLSM(JLON) > 0.5_JPRB
      ENDDO
    ENDIF
  ENDIF

!         Appel du calcul de nebulosite.

  IF(LNEBN) THEN
    CALL ACNEBN (YDCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTNEBU, YDCPG_OPTS%KFLEVG,  &
    & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZQL, ZQI,                                        &
    & ZFLU_QSAT, YDMF_PHYS_BASE_STATE%T, ZPFL_FPLCH, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, ZUNEBH, ZNEBS0, ZQLIS0, ZQLI_CVP, ZNEB_CVPP, ZQLI_CVPP,                                  &
    & ZAIPCMT, YDCPG_MISC%NEB, ZNEBC0, YDCPG_MISC%QICE, YDCPG_MISC%QLI, YDSTA)
    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
!DEC$ IVDEP
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDCPG_MISC%NEB(JLON,JLEV)=YDCPG_MISC%NEB(JLON,JLEV)*GAEPS
      ENDDO
    ENDDO
  ENDIF

  IF ( LNEBR ) THEN
    CALL ACNEBR ( YDERAD, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,              &
    & NTCOEF, NTNEBU, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,  &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, ZMSC_LSCPE, ZQV,              &
    & ZFLU_QSAT, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U,            &
    & YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,   &
    & YDMF_PHYS%OUT%GZ0, ZPFL_FPLCH, YDMF_PHYS_SURF%GSD_VH%PPBLH, YDCPG_MISC%QS, YDMF_PHYS_BASE_STATE%YGSP_RR%T, &
    & ZKTROV, ZKUROV, ZNBVNO, YDCPG_MISC%NEB, ZNEBS, YDCPG_MISC%QICE, YDCPG_MISC%QLI, ZQLIS)
  ENDIF

!         Diagnostique de nebulosite partielle.
  CALL ACNPART(YDCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTNEBU, YDCPG_OPTS%KFLEVG, &
  & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,                  &
  & ZDECRD, YDCPG_MISC%NEB, YDMF_PHYS%OUT%CLCH, YDMF_PHYS%OUT%CLCM, YDMF_PHYS%OUT%CLCL, YDCPG_MISC%CLCT,                           &
  & ZCLCT_RAD, PCLCC=YDMF_PHYS%OUT%CLCC, PNEBC=ZNEBC0, PTOPC=YDMF_PHYS%OUT%CTOP)

!     7.3.5 Computation of the equivalent coefficients for simplified
!           radiation scheme

  IF ( LRAYSP .AND.(YDCPG_OPTS%NSTEP == 1).AND.LRCOEF ) THEN
    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZZNEB(JLON,JLEV)=0.0_JPRB
      ENDDO
    ENDDO

    CALL ACRADCOEF ( YDRCOEF, YDPHY3, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA, &
    & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,        &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZZNEB, ZQV, ZQCO2, YDCPG_MISC%QICE, YDCPG_MISC%QLI,            &
    & ZQO3, YDMF_PHYS_BASE_STATE%T, ZRDG_MMU0, YDCPG_MISC%DHSF, ZFLU_EMIS, ZAER, ZAC, ZAC_HC,                &
    & ZRDT_COR, ZRDT_RAB3C, ZRDT_RAB3N, ZRDT_RAB4C, ZRDT_RAB4N, ZRDT_RAB6C, ZRDT_RAB6N, ZRDT_RAT1C,          &
    & ZRDT_RAT1N, ZRDT_RAT2C, ZRDT_RAT2N, ZRDT_RAT3C, ZRDT_RAT3N, ZRDT_RAT4C, ZRDT_RAT4N, ZRDT_RAT5C,        &
    & ZRDT_RAT5N)
  ENDIF

!     7.3.6 Module chimique - Chemistry module

  IF (LCHEM_ARPCLIM) THEN ! at this stage call when ARPEGE-Climat chemistry only

     ! initialisation below needs to be refined later for more general use
     IFLDX  = 1_JPIM
     IFLDX2 = 1_JPIM
     ILEVX  = 1_JPIM
     ALLOCATE(ZSD_XA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,IFLDX), ZSD_X2(YDCPG_OPTS%KLON,IFLDX2))
     ALLOCATE(INDCHEM(NCHEM),IGPLAT(YDCPG_OPTS%KLON))
     ALLOCATE(ZCFLX(YDCPG_OPTS%KLON,NCHEM), ZCFLXO(YDCPG_OPTS%KLON,NCHEM), ZCHEMDV(YDCPG_OPTS%KLON,0)) ! no species with dry deposition in ARPCLIM
     ALLOCATE(ZAEROP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,NACTAERO),ZTENC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,NCHEM))
     ALLOCATE(ZDELP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZWND(YDCPG_OPTS%KLON),ZDUMMY1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZGELAT(YDCPG_OPTS%KLON))
     ALLOCATE(ZNEEFLX(YDCPG_OPTS%KLON),ZCHEM2AER(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,6))
     INDCHEM(:)     = 1_JPIM ! we should have here the indexes of the chemical species in the YCHEM array, to be implemented later
     IGPLAT (:)     = 1_JPIM
     ZSD_XA (:,:,:) = 0._JPRB
     ZSD_X2 (:,:)   = 0._JPRB
     DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
       DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
         ZDELP(JLON,JLEV) = YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(JLON,JLEV) - YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(JLON,JLEV-1)
       ENDDO
     ENDDO
     ZWND  (:)      = 0._JPRB  ! not used in  ARPEGE-Climat chemistry
     ZDUMMY1 (:,:)  = 1.0E-18_JPRB
     ZGELAT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA) = ASIN(YDVARS%GEOMETRY%GEMU%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA))
     ZCFLX (:,:)        = 0._JPRB
     ZCFLXO (:,:)       = 0._JPRB
     ZCHEMDV (:,:)      = 0._JPRB
     ZAEROP (:,:,:)     = 0._JPRB  ! no interaction aerosol/chemistry in ARPCLIM
     ZTENGFL(:,:,:)     = 0._JPRB  ! only used later for diagnostics
     ZTENC (:,:,:)      = 0._JPRB  ! not used in  ARPEGE-Climat chemistry
     ZNEEFLX (:)        = 0._JPRB  ! not used in  ARPEGE-Climat chemistry
     ZCHEM2AER (:,:,:)  = 0._JPRB  ! not used in  ARPEGE-Climat chemistry
     
!#if false
! !   YDVAB YDDIMV not used for ARPEGE-Climat chemistry
!      CALL CHEM_MAIN &
!    &( YDVAB, YDDIMV, YDMODEL, KIDIA , KFDIA , KLON , KLEV, KVCLIS, NCHEM, INDCHEM,&
!    &  TSPHY , IGPLAT,  IFLDX  , IFLDX2 , ILEVX,&
!    &  ZSD_XA , ZSD_X2,  ZDELP, PAPRS, PAPRSF, PAPHI, PQ, PT,&
!    &  ZDUMMY1, ZDUMMY1, ZDUMMY1, ZDUMMY1,  PNEB,& 
!    &  PFPLCL, PFPLCN, PFPLSL, PFPLSN, ZDUMMY1,& 
!    &  PALB, ZWND, PLSM,& 
!    &  PMU0, ZGELAT, PGELAM, PGEMU, PKOZO, ZCFLX, ZCFLXO, ZCHEMDV, PGFL,&
!    &  ZAEROP, ZTENGFL, PCHEM, ZTENC, ZNEEFLX, ZCHEM2AER )
!#endif     
  ENDIF

  
!      7.4 Rayonnement Geleyn
!          Geleyn's radiation

  IF ( LRAY ) THEN
    SELECT CASE (NRAY)
      CASE(1)
        CALL ACRANEB(YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                       &
        & NTRADI, YDCPG_OPTS%KFLEVG, IJN, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,         &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%NEB, ZQV, ZQCO2, YDCPG_MISC%QICE, YDCPG_MISC%QLI,              &
        & ZQO3, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS%OUT%ALB, ZALBDIR, ZFLU_EMIS, ZRDG_MU0, YDVARS%GEOMETRY%GEMU%T0,           &
        & YDVARS%GEOMETRY%GELAM%T0, ZRDG_MU0LU, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH,     &
        & ZFRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOLU, YDMF_PHYS%OUT%FRTHDS, ZAER, ZMAK,                            &
        & ZMAN)

        ! update sunshine duration [s]
!DEC$ IVDEP
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          IF ( YDMF_PHYS%OUT%FRSOPS(JLON) > RSUNDUR*ZRDG_MU0(JLON) ) THEN
            !YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)=YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)+TSTEP
            YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)=YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)+ZADTMS ! fix stepx case
          ENDIF
        ENDDO
      CASE(2)
        CALL ACRANEB2(YDERDI, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,              &
        & NTRADI, YDCPG_OPTS%KFLEVG, IJN, YDCPG_OPTS%NSTEP, YDCFU%NFRRC, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,              &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%NEB, ZQV, ZQCO2, YDCPG_MISC%QICE, YDCPG_MISC%QLI,              &
        & ZQO3, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS%OUT%ALB, ZALBDIR, ZFLU_EMIS, YDVARS%GEOMETRY%GELAM%T0,                    &
        & YDVARS%GEOMETRY%GEMU%T0, ZRDG_MU0, ZRDG_MU0LU, YDMF_PHYS_BASE_STATE%YGSP_RR%T, ZDECRD,                            &
        & ZCLCT_RAD, YDMF_PHYS%OPT%GDEOSI, YDMF_PHYS%OPT%GUEOSI, YDMF_PHYS%OPT%GMU0, YDMF_PHYS%OPT%GMU0_MIN,                &
        & YDMF_PHYS%OPT%GMU0_MAX, YDMF_PHYS%OPT%GDEOTI, YDMF_PHYS%OPT%GDEOTI2, YDMF_PHYS%OPT%GUEOTI,                        &
        & YDMF_PHYS%OPT%GUEOTI2, YDMF_PHYS%OPT%GEOLT, YDMF_PHYS%OPT%GEOXT, YDMF_PHYS%OPT%GRPROX, YDMF_PHYS%OPT%GMIXP,       &
        & YDMF_PHYS%OPT%GFLUXC, YDMF_PHYS%OPT%GRSURF, YDMF_PHYS_SURF%GSD_VD%PSUND, YDMF_PHYS%OUT%FRSO,                      &
        & YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRTHC, ZFRSODS, YDMF_PHYS%OUT%FRSOPS,                      &
        & YDMF_PHYS%OUT%FRSOLU, YDMF_PHYS%OUT%FRTHDS, ZAER)
    ENDSELECT

    ! sum downward diffuse and direct solar radiation at surface
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%FRSODS(JLON)=ZFRSODS(JLON)+YDMF_PHYS%OUT%FRSOPS(JLON)
    ENDDO

!      7.5 Rayonnement Morcrette
!          Morcrette's radiation

  ELSEIF ( LRAYFM ) THEN

    LLCALLRAD=(MOD(YDCPG_OPTS%NSTEP,NRADFR) == 0 )
!    IF (NCALLRAD==1) ! <== not yet
    IF (NCALLRAD==2) LLCALLRAD=(LLCALLRAD.AND.(YDCPG_OPTS%NSTEP<=NSTOP-1))
!    IF (NCALLRAD==3) ! <== not yet
    IAERO=SIZE(ZAERO,3)
    ! ---- Intermittent call to radiation scheme
    IF (LLCALLRAD) THEN
      CALL RECMWF(YDGEOMETRY%YRDIMV, YDMODEL, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,    &
      & YDCPG_OPTS%KSW, &
      & NOZOCL, NAERMACC, IAERO, &
      & ZALBD, ZALBP, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,       &
      & YDCPG_MISC%NEB, ZQO3, ZDUM,ZDUM,ZDUM, &
      & ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,&
      & ZAER, ZAERO, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZFLU_EMIS, &
      & ZRDG_MU0M, ZQV, ZFLU_QSAT, YDCPG_MISC%QICE, YDCPG_MISC%QLI, &
      & ZQS, ZQR, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%YGSP_RR%T, &
      & YDMF_PHYS%RAD%EMTD, YDMF_PHYS%RAD%EMTU, YDMF_PHYS%RAD%TRSW, YDMF_PHYS%OUT%FRTHC, &
      & YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRSO, ZSFSWDIR, ZSFSWDIF, ZFSDNN,                         &
      & ZFSDNV, ZCTRSO, ZCEMTR, ZTRSOD, ZTRSODIR, &
      & ZTRSODIF, ZPIZA_DST, ZCGA_DST, ZTAUREL_DST,                            &
      & ZAERINDS, YDVARS%GEOMETRY%GELAM%T0, YDVARS%GEOMETRY%GEMU%T0, YDCPG_GPAR%SWDIR, YDCPG_GPAR%SWDIF,                 &
      & ZRDG_MU0LU, YDMF_PHYS%OUT%ALB, YDMF_PHYS%RAD%RMOON)
    ELSE
      IF (LMSE) THEN
      DO JSG=1,NSW
        ZTRSODIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSG)=YDCPG_GPAR%SWDIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSG)
        ZTRSODIF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSG)=YDCPG_GPAR%SWDIF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSG)
      ENDDO
      ENDIF
    ENDIF

    IF (LRAYLU)  YDMF_PHYS%OUT%FRSOLU(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS%RAD%RMOON(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)

    ! ---- Flux update and radiative heating rates
    CALL RADHEAT (YDMODEL%YRCST,  YDMODEL%YRML_PHY_EC%YRTHF, YDERAD, YDERDI, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, &
    & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, ZFLU_EMIS, YDMF_PHYS%RAD%EMTD,                &
    & ZRDG_MU0, ZQV, ZTENT, YDMF_PHYS%RAD%TRSW, ZTRSOD, YDMF_PHYS_BASE_STATE%YGSP_RR%T, TSPHY,               &
    & ZTRSODIR, ZTRSODIF, ZALBD, ZALBP, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRSODS,        &
    & YDMF_PHYS%OUT%FRTHDS, ZCEMTR, ZCTRSO, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRTHC, ZSUDU, ZSDUR,          &
    & ZDSRP, ZSFSWDIR, ZSFSWDIF, YDMF_PHYS%OUT%FRSOPS, ZFRSODS, YDMF_PHYS%OUT%FRSOPT )


    ! ---- Take into account day duration depending on altitude.
    IF(LDAYD) CALL ACDAYD(YDCST, YDRIP, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,    &
              & YDCPG_OPTS%KTDIA, NTSSG, YDVARS%GEOMETRY%GEMU%T0, ZRDG_MU0, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,            &
              & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS%OUT%FRSO)

    ! ---- Correct solar absorption as a function of pmu0.
    IF(GRSO < 1._JPRB) CALL ACRSO(YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,        &
                       & YDCPG_OPTS%KTDIA, NTSSG, YDVARS%GEOMETRY%GEMU%T0, ZRDG_MU0, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, &
                       & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS%OUT%FRSO)

    IF(.NOT.LMSE) THEN
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZALB(JLON)=0.0_JPRB
        DO JSG=1,NSW
          ZALB(JLON)=ZALB(JLON)+0.5_JPRB*(ZALBD(JLON,JSG)+ZALBP(JLON,JSG))
        ENDDO
        ZALB(JLON)=ZALB(JLON)/FLOAT(NSW)
        YDMF_PHYS%OUT%FRSODS(JLON)=YDMF_PHYS%OUT%FRSO(JLON,YDCPG_OPTS%KFLEVG,1)/(1.0_JPRB-ZALB(JLON))
      ENDDO
    ENDIF
    ! Compute Sunshine Duration (in seconds)
!DEC$ IVDEP
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      IF(YDMF_PHYS%OUT%FRSODS(JLON) >= RSUNDUR) THEN
        !YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)=YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)+1.0_JPRB*TSTEP
        YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)=YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)+1.0_JPRB*ZADTMS ! fix stepx case
      ENDIF
    ENDDO

  ENDIF

  IF (LRAY.OR.LRAYFM ) THEN

    ! Direct normal irradiance with securities
    DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%FRSDNI(JLON)=YDMF_PHYS%OUT%FRSOPS(JLON)
      IF (ZRDG_MU0(JLON) > 3.0E-02_JPRB) THEN
        YDMF_PHYS%OUT%FRSDNI(JLON)=YDMF_PHYS%OUT%FRSOPS(JLON)/ZRDG_MU0(JLON)
      ENDIF
      YDMF_PHYS%OUT%FRSDNI(JLON)=MAX(0.0_JPRB,YDMF_PHYS%OUT%FRSDNI(JLON))
    ENDDO
  ENDIF

  IF (LRAY.OR.LRAYFM) THEN

    ! global normal irradiance
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%FRSGNI(JLON)=YDMF_PHYS%OUT%FRSDNI(JLON)+0.5_JPRB*(                         &
       & (1.0_JPRB+ZRDG_MU0(JLON))*(YDMF_PHYS%OUT%FRSODS(JLON)-YDMF_PHYS%OUT%FRSOPS(JLON)       )+ &
       & (1.0_JPRB-ZRDG_MU0(JLON))*(YDMF_PHYS%OUT%FRSODS(JLON)-YDMF_PHYS%OUT%FRSO  (JLON,YDCPG_OPTS%KFLEVG,1)))
    ENDDO

    ! mean radiant temperature
    IF (LXMRT) THEN
      CALL MEAN_RAD_TEMP(YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, ZRDG_MU0, YDMF_PHYS%OUT%FRSO(:, YDCPG_OPTS%KFLEVG, 1), &
      & YDMF_PHYS%OUT%FRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRTH(:, YDCPG_OPTS%KFLEVG, 1), YDMF_PHYS%OUT%FRTHDS,               &
      & YDMF_PHYS%OUT%MRT) 
    ENDIF

  ENDIF

!*
!     ------------------------------------------------------------------

!     7.BIS. BILAN HYDRIQUE DU SOL
!     ----------------------------
!     CALCUL DES RESISTANCES A L'EVAPOTRANSPIRATION HV ET
!     A LA TRANSPIRATION
!     ------------------------------------------------------------------
!     HTR DU COUVERT VEGETAL
!     ----------------------

  IF (LSOLV.AND.(.NOT.LMSE)) THEN
    CALL ACVEG ( YDPHY, YDPHY1, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,                  &
    & YDMF_PHYS%OUT%FRSO, ZQV, ZFLU_QSAT, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_SURF%GSD_VV%PD2, YDMF_PHYS_SURF%GSD_VF%PLSM, &
    & YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS_SURF%GSD_VV%PLAI, ZFLU_NEIJ, ZFLU_VEG, YDMF_PHYS_SURF%GSD_VV%PRSMIN,        &
    & ZCHROV, ZGWDCS, ZWFC, YDMF_PHYS_BASE_STATE%YGSP_RR%FC, ZWWILT, YDMF_PHYS_BASE_STATE%YGSP_SB%Q,                     &
    & ZFLU_QSATS, ZHQ, ZHTR, ZHU, YDMF_PHYS_SURF%GSD_VV%PHV, ZWLMX)
  ENDIF

!     ------------------------------------------------------------------
!     8.- DIFFUSION VERTICALE TURBULENTE
!     ----------------------------------
  IF ( LVDIF ) THEN

!   Sauvegarde temporaire de l'ancien acdifus pour les besoins du Climat
    IF ( LACDIFUS ) THEN

      IF(NDIFFNEB == 1) THEN
        ZNEBDIFF(:,:)=ZNEBS(:,:)
      ELSEIF(NDIFFNEB == 2) THEN
        ZNEBDIFF(:,:)=YDCPG_MISC%NEB(:,:)
      ELSEIF(NDIFFNEB == 3) THEN
        ZNEBDIFF(:,:)=ZNEBS(:,:)+(1.0_JPRB-ZNEBS(:,:))*ZNEBCH(:,:)
      ENDIF

      CALL ACDIFUS ( YDMCC, YGFL, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                &
      & NTDIFU, YDCPG_OPTS%KFLEVG, YSP_SBD%NLEVS, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,   &
      & YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS%OUT%FRSO,                  &
      & ZKTROV, ZKUROV, ZKNROV, ZQV, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_BASE_STATE%T,                      &
      & YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZXTROV, ZXUROV, ZXPTKEROV, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, ZQL, ZQI,                               &
      & ZNEBDIFF, YDMF_PHYS_BASE_STATE%TKE, ZLMU, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZFLU_CD,                         &
      & ZFLU_CDN, ZCDROV, ZCHROV, ZCEROV, ZDSA_CPS, YDMF_PHYS%OUT%CT, ZDQSTS, ZFLU_EMIS, YDMF_PHYS_SURF%GSD_VH%PSPSH,      &
      & ZHQ, ZHTR, ZHU, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG,                &
      & ZFLU_NEIJ, YDCPG_MISC%QS, ZFLU_QSATS, YDMF_PHYS_BASE_STATE%YGSP_SB%T, YDMF_PHYS_BASE_STATE%YGSP_RR%T,              &
      & ZFLU_VEG, ZXDROV, ZXHROV, YDMF_PHYS_BASE_STATE%YGSP_RR%W, YDMF_PHYS_BASE_STATE%YGSP_RR%IC, YDMF_PHYS%OUT%DIFTQ,    &
      & YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV, YDMF_PHYS%OUT%DIFTQL,          &
      & YDMF_PHYS%OUT%DIFTQN, ZTENDPTKE, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN,                      &
      & YDMF_PHYS%OUT%FCS, ZFLU_FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FTR,       &
      & ZDSA_LHS, YDMF_PHYS%OUT%FRTH, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H)

    ELSE

      IF ( LPTKE ) THEN
        IF ( LMSE.AND.LCALLSFX ) THEN       
          IF (YDCPG_OPTS%NSTEP == 0) THEN
            ZFLU_CD(:)=ZFLU_CDN(:)   ! very first approximation
          ELSE
            ZFLU_CD(:)=MAX(YDCPG_GPAR%CD(:),ZEPS0)
          ENDIF
        ENDIF
        CALL ACPTKE(YGFL, YDLDDH, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                             &
        & YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, YDCPG_OPTS%NSTEP, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, ZKNROV, YDMF_PHYS_BASE_STATE%T,              &
        & YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%TKE,                         &
        & ZLMU, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZF_EPS, ZFUN_TTE, ZMRIPP, ZMRIFPP, ZRHS,                                        &
        & ZMN2_ES, ZMN2_EQ, ZRRCOR, YDVARS%FQTUR%T0, YDVARS%FSTUR%T0, YDVARS%SHTUR%T0, ZFLU_CD, YDMF_PHYS%OUT%GZ0,                      &
        & ZKUROV, ZKTROV, ZXPTKEROV, ZLMLTILD, YDVARS%MXL%T0, ZLML, ZTH_FUN, ZWW_FUN, ZTENDPTKE, YDVARS%TTE%T0,                         &
        & YDCPG_MISC%FTCNS)
      ENDIF

      CALL ACDIFV1 (YDCST, YGFL, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                   &
      & YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,              &
      & ZCP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZKTROV, ZKQROV, ZKUROV, ZQV, ZQL, ZQI, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, &
      & YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZSVM, ZXTROV, ZXUROV,                                &
      & ZXDROV, ZXHROV, ZEDMFS, ZEDMFQ, ZEDMFU, ZEDMFV, ZMF_UP, ZXURO, ZXQRO, ZXTRO, ZCFAQ, ZCFAS,                                   &
      & ZCFATH, ZCFAU, ZCFASV, ZCFBQ, ZCFBS, ZCFBTH, ZCFBU, ZCFBV, ZCFBSV, ZDSE, ZQT)   
        

      IF ( LMSE.AND.LCALLSFX ) THEN

        IF (LRAYFM) THEN
          ZCARDI=RCARDI
        ELSEIF (LRAY) THEN
          ZCARDI=QCO2
          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            ZSFSWDIF(JLON,1)=ZFRSODS(JLON)
            ZSFSWDIR(JLON,1)=YDMF_PHYS%OUT%FRSOPS(JLON)
          ENDDO
        ENDIF

        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZRHODREFM(JLON)=YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,YDCPG_OPTS%KFLEVG)/(YDMF_PHYS_BASE_STATE%T(JLON,YDCPG_OPTS%KFLEVG)*YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(JLON,YDCPG_OPTS%KFLEVG))
          ZDEPTH_HEIGHT(JLON,:)=(YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(JLON,:)-YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,YDCPG_OPTS%KFLEVG))/RG
          ZZS(JLON)=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,YDCPG_OPTS%KFLEVG)/RG
        ENDDO

        IRR=2

        CALL ARO_GROUND_PARAM( YDCPG_BNDS%KBL, YDCPG_OPTS%KGPCOMP, YDCPG_BNDS%KFDIA-YDCPG_BNDS%KIDIA+1,                                                                        &
        & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%NSTEP, IRR, NSW, NGFL_EXT, NDGUNG, NDGUXG,                                                                            &
        & NDLUNG, NDLUXG, LSURFEX_KFROM, LMPA, CCOUPLING, YDCPG_OPTS%LCONFX, NINDAT, ZRHGMT, ZSTATI, RSOVR,                                                                              &
        & RCODEC, RSIDEC, YDVARS%GEOMETRY%RINDX%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDVARS%GEOMETRY%RINDY%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                            &
        & YDMF_PHYS_BASE_STATE%U(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG),                                                                      &
        & YDMF_PHYS_BASE_STATE%V(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG),                                                                      &
        & YDMF_PHYS_BASE_STATE%T(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG),                                                                      &
        & YDMF_PHYS_BASE_STATE%Q(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG),                                                                      &
        & ZSVM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG, 1:NGFL_EXT),                                                                            &
        & ZCARDI, ZRHODREFM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG), &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG),                                                        &
        & ZDTMSE, ZDEPTH_HEIGHT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG), ZZS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                                 &
        & XZSEPS, ZRDG_MU0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZRDG_MU0N(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                                                   &
        & YDVARS%GEOMETRY%GELAM%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDVARS%GEOMETRY%GEMU%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                             &
        & XSW_BANDS, ZSRAIN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZSSNOW(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                                                     &
        & ZSGROUPEL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%FRTHDS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                                               &
        & ZSFSWDIF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:NSW), ZSFSWDIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:NSW),                                                              &
        & ZCFAQ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG), ZCFATH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG),       &
        & ZCFAU(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG), ZCFBQ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG),        &
        & ZCFBTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG), ZCFBU(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG),       &
        & ZCFBV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG:YDCPG_OPTS%KFLEVG), YDMF_PHYS%OUT%FCS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1),                              &
        & ZFEV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZSFSV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:NGFL_EXT),                                                                       &
        & ZSFCO2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZFMDU(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZFMDV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                       &
        & ZALBP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:NSW), ZALBD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:NSW),                                                                    &
        & ZFLU_EMIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZTSN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%FRTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, 1)   &
        &                               )               !orographic shadowing

! Opposite water vapor flux
        ZFEVS(:)=-ZFEV(:)

        CALL ARO_GROUND_DIAG( YDCPG_BNDS%KBL, YDCPG_OPTS%KGPCOMP, YDCPG_BNDS%KFDIA-YDCPG_BNDS%KIDIA+1,                                                                    &
        & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, -1, NDGUNG, NDGUXG, NDLUNG, NDLUXG,                                                                      &
        & LSURFEX_KFROM, ZZS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZFEVS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                                                &
        & YDMF_PHYS_BASE_STATE%U(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), YDMF_PHYS_BASE_STATE%V(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), &
        & ZDEPTH_HEIGHT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), YDMF_PHYS%OUT%FRTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, 1),             &
        & YDMF_PHYS%OUT%FRSO(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, 1), YDVARS%GEOMETRY%RINDX%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                       &
        & YDVARS%GEOMETRY%RINDY%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDCPG_MISC%QS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                                  &
        & YDMF_PHYS%OUT%GZ0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%GZ0H(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                                    &
        & YDMF_PHYS%OUT%TCLS (YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%QCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                                  &
        & YDMF_PHYS%OUT%RHCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%UCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                                  &
        & YDMF_PHYS%OUT%VCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%NUCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                                  &
        & YDMF_PHYS%OUT%NVCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%FCLL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1),                                               &
        & YDMF_PHYS%OUT%FCLN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1), YDMF_PHYS%OUT%FEVL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1),                                             &
        & YDMF_PHYS%OUT%FEVN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1), ZSSO_STDEV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                                        &
        & ZTWSNOW(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZBUDTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZBUDSO(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                               &
        & ZFCLL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZTOWNS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZCD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                    &
        & YDMF_PHYS%OUT%SIC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA))

        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDCPG_GPAR%GZ0(JLON)   = YDMF_PHYS%OUT%GZ0 (JLON)
          YDCPG_GPAR%GZ0H(JLON)  = YDMF_PHYS%OUT%GZ0H(JLON)
          YDCPG_GPAR%VEMIS(JLON) = ZFLU_EMIS(JLON)
          YDCPG_GPAR%VQS(JLON)   = YDCPG_MISC%QS  (JLON)
          YDCPG_GPAR%CD(JLON)    = ZCD  (JLON)
        ENDDO

        DO JSG=1,NSW
          YDCPG_GPAR%ALBSCA(:,JSG) = ZALBD(:,JSG)
          YDCPG_GPAR%ALBDIR(:,JSG) = ZALBP(:,JSG)
        ENDDO

        DO JSG  = 1, NTSSG+1
          DO JLEV = 0, YDCPG_OPTS%KFLEVG
            DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
              YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)=YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)+ZBUDTH(JLON)
            ENDDO
          ENDDO
        ENDDO

! calculation of variables for the old "ISBA" atmosphere scheme
        IF (.NOT. LELAM) THEN
          CALL ARO_GROUND_DIAG_2ISBA( YDCPG_BNDS%KBL, YDCPG_OPTS%KGPCOMP, YDCPG_BNDS%KFDIA-YDCPG_BNDS%KIDIA+1,                                              &
          & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, YDVARS%GEOMETRY%RINDX%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), &
          & YDVARS%GEOMETRY%RINDY%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS_SURF%GSD_VF%PLSM,                                                        &
          & YDMF_PHYS_SURF%GSD_VV%PARG, YDMF_PHYS_SURF%GSD_VV%PSAB, YDMF_PHYS_SURF%GSD_VV%PD2, ZTSN,                                                        &
          & ZTWSNOW, YDMF_PHYS_SURF%GSP_SB%PT_T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1), YDMF_PHYS_SURF%GSP_RR%PW_T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),     &
          & YDMF_PHYS_SURF%GSP_SB%PQ_T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1), YDMF_PHYS_SURF%GSP_RR%PIC_T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),             &
          & YDMF_PHYS_SURF%GSP_SB%PTL_T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1), YDMF_PHYS_SURF%GSP_RR%PFC_T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),            &
          & YDMF_PHYS_SURF%GSP_SG%PA_T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1), YDMF_PHYS_SURF%GSP_SG%PR_T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1),           &
          & ZHV2 )
          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            YDMF_PHYS_SURF%GSD_VV%PHV(JLON) = ZHV2(JLON)
          ENDDO
        ENDIF

        IF (LDIFCONS.AND..NOT.LNODIFQC) THEN
          CALL ACAA1 (YDCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA, &
          & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, ZCOEFN, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP,               &
          & ZQL, ZQI, YDMF_PHYS_BASE_STATE%T, ZALPHA1, ZCOEFA, ZLVT, ZQICE)
        ENDIF

      ELSEIF (.NOT.LMSE) THEN

        IF(NDIFFNEB == 5) THEN
          ZNEBDIFF(:,:)=ZNEBS(:,:)
        ELSEIF(NDIFFNEB == 2) THEN
          ZNEBDIFF(:,:)=YDCPG_MISC%NEB(:,:)
        ELSEIF(NDIFFNEB == 3) THEN
          ZNEBDIFF(:,:)=ZNEBS(:,:)+(1.0_JPRB-ZNEBS(:,:))*ZNEBCH(:,:)
        ENDIF

        IF (.NOT. LSFORCS) THEN
          IF (LEDMFI) THEN
            ZCFBS_G(:,:) = ZCFBS(:,:)
            ZCFBQ_G(:,:) = ZCFBQ(:,:) 
            ZCFBU_G(:,:) = ZCFBU(:,:)
            ZCFBV_G(:,:) = ZCFBV(:,:)
          ENDIF   
          CALL ARP_GROUND_PARAM (YDCST,  YDMCC, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, &
          & YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, YSP_SBD%NLEVS, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI,                                &
          & ZCP, YDMF_PHYS%OUT%FRSO, ZQV, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_BASE_STATE%T,                        &
          & YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZXURO, ZXQRO, ZXTRO, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,            &
          & ZQL, ZQI, ZNEBDIFF, ZFLU_CD, ZFLU_CDN, ZCDROV, ZCHROV, ZCEROV, ZDSA_CPS, YDMF_PHYS%OUT%CT,                            &
          & ZDQSTS, ZFLU_EMIS, YDMF_PHYS_SURF%GSD_VH%PSPSH, ZHQ, ZHTR, ZHU, YDMF_PHYS_SURF%GSD_VV%PHV,                            &
          & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG, ZFLU_NEIJ, YDCPG_MISC%QS,                                    &
          & ZFLU_QSATS, YDMF_PHYS_BASE_STATE%YGSP_SB%T, YDMF_PHYS_BASE_STATE%YGSP_RR%T, ZFLU_VEG,                                 &
          & ZXDROV, ZXHROV, YDMF_PHYS_BASE_STATE%YGSP_RR%W, YDMF_PHYS_BASE_STATE%YGSP_RR%IC, ZDSE,                                &
          & ZCFAS, ZCFAU, ZCFBS, ZCFBU, ZCFBV, ZCFBQ, ZCOEFA, ZALPHA1, ZLVT, ZQICE, ZDIFWQ, ZDIFWS,                               &
          & ZFMDU, ZFMDV, ZSC_FEVI, ZSC_FEVN, ZSC_FCLL, ZSC_FCLN, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL,                        &
          & YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS, ZFLU_FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN,                             &
          & YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FTR, ZDSA_LHS, YDMF_PHYS_SURF%GSD_VH%PQSH, YDMF_PHYS%OUT%FRTH,                      &
          & YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H)

          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            YDMF_PHYS%OUT%DIFTQ(JLON,YDCPG_OPTS%KFLEVG)=ZDIFWQ(JLON)
            YDMF_PHYS%OUT%DIFTS(JLON,YDCPG_OPTS%KFLEVG)=ZDIFWS(JLON)
          ENDDO

          IF (LEDMFI) THEN
            ZCFBQ(:,:)  = ZCFBQ_G(:,:)
            ZCFBS(:,:)  = ZCFBS_G(:,:) 
            ZCFBU(:,YDCPG_OPTS%KFLEVG)  = ZCFBU_G(:,YDCPG_OPTS%KFLEVG)
            ZCFBV(:,YDCPG_OPTS%KFLEVG)  = ZCFBV_G(:,YDCPG_OPTS%KFLEVG) 
          ENDIF   

        ENDIF   ! <== LSFORCS

      ELSEIF (.NOT. LCALLSFX) THEN
        YDMF_PHYS%OUT%FCS(:,:) = 0.0_JPRB
        ZFEV(:)   = 0.0_JPRB

      ENDIF  !LMSE.AND.LCALLSFX

      CALL ACDIFV2 (LSFORCS, YDCST,  YGFL, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,&
      & YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, ZCFAQ, ZCFAS,                            &
      & ZCFAU, ZCFASV, ZCFBQ, ZCFBS, ZCFBU, ZCFBV, ZCFBSV, ZKTROV, ZKQROV, ZKQLROV, ZKUROV, ZDSE,                           &
      & ZQT, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZPOID, YDMF_PHYS_BASE_STATE%T, ZQL, ZQI,                       &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, ZCOEFA, ZALPHA1, ZLVT, ZQICE, ZSFSV, YDMF_PHYS%OUT%FCS,                    &
      & ZFEV, ZFMDU, ZFMDV, ZTSN, ZXHROV, ZDIFEXT, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%STRTU,           &
      & YDMF_PHYS%OUT%STRTV, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS,          &
      & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, YDVARS%SHTUR%T0          )

      IF (LEDKF) THEN
        IF (LSFORCS) THEN
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
           YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,1) = YDMF_PHYS%OUT%FCS(JLON,1)
           YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,2) = ZFEV(JLON)
        ENDDO
        ELSE
          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA             
            YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,1) = YDMF_PHYS%OUT%DIFTS(JLON,YDCPG_OPTS%KFLEVG)
            YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,2) = YDMF_PHYS%OUT%DIFTQ(JLON,YDCPG_OPTS%KFLEVG)
          ENDDO
        ENDIF   
      ENDIF

      IF ( LMSE ) THEN
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDCPG_GPAR%VTS(JLON) = ZTSN (JLON)
          YDMF_PHYS_SURF%GSP_SG%PF_T1(JLON,1) = ZTWSNOW (JLON)
        ENDDO
      ENDIF

      ! First compute horizontal exchange coefficients for momentum:
      !  (there's mo TOMs contribution, thus has to be done at latest here)
      IF (YDDYNA%L3DTURB) THEN
        CALL ACTKECOEFKH(YDRIP, YDMODEL%YRML_PHY_MF, YDDYNA, YDGEOMETRY%YREGEO, &
        & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,                      &
        & YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%TKE, ZTENDPTKE,                             &
        & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%T,                         &
        & YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%DIV,      &
        & YDMF_PHYS_BASE_STATE%VOR, YDVARS%U%DL, YDVARS%V%DL, YDMF_PHYS_BASE_STATE%YCPG_PHY%W, YDMF_PHYS_BASE_STATE%YCPG_PHY%WL, &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%WM, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,           &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZLML, ZFHORM, ZFHORH, ZKUROV, ZRTI, ZFLU_CD,                                   &
        & ZCDROV, LPTKE, ZKUR_KUROV_H, ZKUR_KTROV_H, ZRHS)
      ENDIF

      IF (LCOEFKTKE) THEN
        IF (NDIFFNEB == 1) THEN
          ZCOEFA(:,:) = ZBNEBQ(:,:)
        ELSEIF (NDIFFNEB == 4) THEN
          ZCOEFA(:,:) = ZBNEBCVPP(:,:)
        ENDIF
        CALL ACDIFV3 ( YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,       &
        & YDCPG_OPTS%KFLEVG, YDCPG_OPTS%NSTEP, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI,  &
        & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZCP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZKTROV,                       &
        & ZXTROV, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                          &
        & ZCHROV, ZXHROV, ZCOEFA, ZQV, ZQL, ZQI, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%T, &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZLML,                          &
        & ZTH_FUN, ZWW_FUN, ZF_EPS, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS_BASE_STATE%TKE,                  &
        & ZTENDPTKE, ZMN2_ES, ZMN2_EQ, ZMN2_DS, ZMN2_DQ, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTQN,                     &
        & ZDIFWQ, ZDIFWS, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS, ZFLU_FEVI, YDMF_PHYS%OUT%FEVL,      &
        & YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%GZ0, ZRTI, ZSC_FEVI, ZSC_FEVN, ZSC_FCLL, ZSC_FCLN,                           &
        & ZTSTAR, ZTSTAR2, ZTSTARQ, ZTSTAR2Q)

        !store fluxes and shear term
        !GFL fields are on full levels, fluxes on half levels
        IF (YFQTUR%LGP.AND.YFSTUR%LGP) THEN
          DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
            DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
              YDVARS%FQTUR%T0(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)+YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)&
               &                                 +YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)
              YDVARS%FSTUR%T0(JLON,JLEV)=YDMF_PHYS%OUT%DIFTS(JLON,JLEV)
            ENDDO
          ENDDO
        ENDIF
      ENDIF ! LCOEFKTKE

      ! Now the heat coefficient can be completed by TKE+ containing
      !  the TOMs contribution. 
      IF (YDDYNA%L3DTURB) THEN
        CALL ACTKECOEFKH(YDRIP, YDMODEL%YRML_PHY_MF, YDDYNA, YDGEOMETRY%YREGEO, &
        & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,                      &
        & YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%TKE, ZTENDPTKE,                             &
        & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%T,                         &
        & YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%DIV,      &
        & YDMF_PHYS_BASE_STATE%VOR, YDVARS%U%DL, YDVARS%V%DL, YDMF_PHYS_BASE_STATE%YCPG_PHY%W, YDMF_PHYS_BASE_STATE%YCPG_PHY%WL, &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%WM, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,           &
        & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZLML, ZFHORM, ZFHORH, ZKUROV, ZRTI, ZFLU_CD,                                   &
        & ZCDROV, LPTKE, ZKUR_KUROV_H, ZKUR_KTROV_H, ZRHS)
      ENDIF

    ENDIF

!-----------------------------------------------------------------------------
!   THE DEEP CONVECTION WILL SEE THE SHALLOW PART FROM KFB AS IT IS WITH LOUIS
!   SCHEME AND THE MODIFIED RI
!----------------------------------------------------------------------------
    IF (LCVPPKF.OR.(LEDKF .AND. .NOT. LEDMFI)) THEN
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDMF_PHYS%OUT%DIFTQ   (JLON,JLEV) = YDMF_PHYS%OUT%DIFTQ(JLON,JLEV) + ZDIFCVPPQ(JLON,JLEV)
          YDMF_PHYS%OUT%DIFTS   (JLON,JLEV) = YDMF_PHYS%OUT%DIFTS(JLON,JLEV) + ZDIFCVPPS(JLON,JLEV)
          YDMF_PHYS%OUT%STRTU   (JLON,JLEV) = YDMF_PHYS%OUT%STRTU(JLON,JLEV) + ZDIFCVPPU(JLON,JLEV)
          YDMF_PHYS%OUT%STRTV   (JLON,JLEV) = YDMF_PHYS%OUT%STRTV(JLON,JLEV) + ZDIFCVPPV(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

    IF (L3MT) THEN

    ! ------------------------------------------------------------------
    ! UPDATE TEMPERATURE, LIQUID WATER AND ICE BY THE TURBULENT FLUXES
    ! INCLUDING CORRECTION OF NEGATIVE VALUES OF WATER SPECIES
    ! SAVE THE INCREMENTAL PART DUE TO THE TURBULENT DIFFUSION PROCESSES
    ! ------------------------------------------------------------------

!   setup of auxiliary variables for water vapour updating for RK condensation scheme

    IF(LRKCDEV) THEN
      ZDTURDIFF=1._JPRB
      IF (LCVGQD) ZDTURDIFF=0._JPRB
    ENDIF

!cdir unroll=8
    DO JLEV = YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
!DEC$ IVDEP
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA

        ZT(JLON,JLEV)=YDMF_PHYS_BASE_STATE%T(JLON,JLEV)-ZIPOI(JLON,JLEV)/YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)&
         & *(YDMF_PHYS%OUT%DIFTS(JLON,JLEV)-YDMF_PHYS%OUT%DIFTS(JLON,JLEV-1))
        ZQX0=ZQI(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQN(JLON,JLEV-1))
        ZQI(JLON,JLEV)=MAX(0.0_JPRB,ZQX1)
        ZDQI=MAX(0.0_JPRB,ZQX1)-ZQX1
        ZFCQING(JLON,JLEV)=ZFCQING(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)+ZFCQING(JLON,JLEV)
        ZQX0=ZQL(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQL(JLON,JLEV-1))
        ZQL(JLON,JLEV)=MAX(0.0_JPRB,ZQX1)
        ZDQL= MAX(0.0_JPRB,ZQX1)-ZQX1
        ZFCQLNG(JLON,JLEV)=ZFCQLNG(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)+ZFCQLNG(JLON,JLEV)
        ZDQC=ZDQI+ZDQL

        ZQX0=ZQV(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQ(JLON,JLEV-1))
        ZQV0=ZQX1-ZIPOI(JLON,JLEV)*(0.0_JPRB&
          & -ZFCQVNG(JLON,JLEV-1)-ZFCQING(JLON,JLEV-1)&
          & -ZFCQLNG(JLON,JLEV-1))
        IF(LCVGQD) THEN
          ZQV(JLON,JLEV)=MAX(0.0_JPRB,ZQV0-ZDQC)
        ENDIF
        ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-ZQX1
        ZDQVDIFF(JLON,JLEV)=ZDQV+ZQX1-ZQX0
        ZFCQVNG(JLON,JLEV)=ZFCQVNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV)+ZFCQVNG(JLON,JLEV)
        IF(LRKCDEV.AND.YDCPG_OPTS%NSTEP>0) THEN
          ZDTRAD(JLON,JLEV)=ZIPOI(JLON,JLEV)/YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)* (&
            & YDMF_PHYS%OUT%FRSO(JLON,JLEV-1,1)-YDMF_PHYS%OUT%FRSO(JLON,JLEV,1)&
            & +YDMF_PHYS%OUT%FRTH(JLON,JLEV-1,1)-YDMF_PHYS%OUT%FRTH(JLON,JLEV,1) )
          ZRKTTEND(JLON,JLEV)=(ZT(JLON,JLEV)-YDVARS%RKTH%T0(JLON,JLEV)&
             & +ZDTRAD(JLON,JLEV))/TSPHY
          ZRKQVTEND(JLON,JLEV)=(MAX(0._JPRB,ZQV(JLON,JLEV)&
             & +ZDQVDIFF(JLON,JLEV)*ZDTURDIFF)&
             & -MAX(0._JPRB,YDVARS%RKTQV%T0(JLON,JLEV)))/TSPHY
          ZRKQCTEND(JLON,JLEV)= (MAX(0._JPRB,ZQI(JLON,JLEV)+ZQL(JLON,JLEV))&
             & -MAX(0._JPRB,YDVARS%RKTQC%T0(JLON,JLEV)))/TSPHY
        ENDIF
      ENDDO
    ENDDO

    ELSEIF (LSTRAPRO) THEN

    ! ------------------------------------------------------------------
    ! UPDATE THE CORRECTION FLUXES FOR NEGATIVE VALUES OF WATER SPECIES
    ! SAVE THE INCREMENTAL PART DUE TO THE TURBULENT DIFFUSION PROCESSES
    ! ------------------------------------------------------------------
    DO JLEV = YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA

        ZT (JLON,JLEV)= YDMF_PHYS_BASE_STATE%T(JLON,JLEV)
        ZQX0=ZQI(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQN(JLON,JLEV-1))
        ZDQI= MAX(0.0_JPRB,ZQX1)-ZQX1
        ZFCQING(JLON,JLEV)=ZFCQING(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)+ZFCQING(JLON,JLEV)
        ZQX0=ZQL(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQL(JLON,JLEV-1))
        ZDQL= MAX(0.0_JPRB,ZQX1)-ZQX1
        ZFCQLNG(JLON,JLEV)=ZFCQLNG(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)+ZFCQLNG(JLON,JLEV)
        ZDQC=ZDQI+ZDQL

        ZQX0=ZQV(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQ(JLON,JLEV-1))
        ZQV0=ZQX1-ZIPOI(JLON,JLEV)*(0.0_JPRB&
          & -ZFCQVNG(JLON,JLEV-1)-ZFCQING(JLON,JLEV-1)&
          & -ZFCQLNG(JLON,JLEV-1))
        ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-ZQX1
        ZFCQVNG(JLON,JLEV)=ZFCQVNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV)+ZFCQVNG(JLON,JLEV)

      ENDDO
    ENDDO

    ENDIF ! 3MT || LSTRAPRO

    ! ------------------------------------------------------
    ! DIAGNOSTIC DE LA HAUTEUR DE COUCHE LIMITE, SELON TROEN ET MAHRT,
    ! POUR USAGE AU PAS DE TEMPS SUIVANT.
    ! ------------------------------------------------------

    IF((TRIM(CGMIXLEN) == 'TM').OR.(TRIM(CGMIXLEN) == 'TMC')) THEN
      ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=MIN(XMAXLM,MAX(XMINLM,YDMF_PHYS_SURF%GSD_VH%PPBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)))
      CALL ACPBLHTM ( YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTDIFU, YDCPG_OPTS%KFLEVG,              &
      & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,     &
      & YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS%OUT%STRTU,                   &
      & YDMF_PHYS%OUT%STRTV, ZPCLS, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%QCLS, YDMF_PHYS%OUT%FCS, YDMF_PHYS%OUT%FCLL,         &
      & YDMF_PHYS%OUT%FCLN, ZDSA_LHS, YDMF_PHYS%OUT%GZ0, IMOC_CLPH, ZBLH, ZUGST, ZVGST)
      YDMF_PHYS_SURF%GSD_VH%PPBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      IF (.NOT.LRAFTUR) THEN
        YDMF_PHYS%OUT%UGST(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZUGST(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
        YDMF_PHYS%OUT%VGST(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZVGST(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      ENDIF
    ENDIF

    IF(LNEBCO.AND.LNEBR)THEN
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZPREN(JLON)=YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(JLON,YDCPG_OPTS%KFLEVG)
      ENDDO
      CALL ACPBLH ( YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTDIFU, YDCPG_OPTS%KFLEVG,                 &
      & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP,        &
      & ZQV, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
      & YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV, ZPREN, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDCPG_MISC%QS,                     &
      & YDMF_PHYS%OUT%FCS, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, ZDSA_LHS, ZRTI, YDMF_PHYS_SURF%GSD_VH%PPBLH              &
      &                                                                                                  )
    ENDIF
    ! ------------------------------------------------------------------
    ! UPDATE PASSIFS SCALAIRS DUE TO THE TURBULENT DIFFUSION PROCESSES
    ! ------------------------------------------------------------------
    IF (LMDUST.AND.(NGFL_EXT/=0)) THEN
      DO JGFL=1,NGFL_EXT
        DO JLON=1,YDCPG_OPTS%KLON
          DO JLEV=1,YDCPG_OPTS%KFLEVG
            ZSVM(JLON,JLEV,JGFL)=ZSVM(JLON,JLEV,JGFL)-ZIPOI(JLON,JLEV)&
          & *(ZDIFEXT(JLON,JLEV,JGFL)-ZDIFEXT(JLON,JLEV-1,JGFL))
            ZSVM(JLON,JLEV,JGFL)=MAX(ZSVM(JLON,JLEV,JGFL),0.0_JPRB)
          ENDDO
        ENDDO
      ENDDO
    ENDIF  ! LMDUST
  ENDIF ! LVDIF

!     DIAGNOSTIC SUPPLEMENTAIRE FLUX DE RAYONNEMENT (PFRTHDS ET PFRSOPT)
!     ADDITIONAL DIAGNOSTICS OF RADIATIVE FLUXES (PFRTHDS AND PFRSOPT)

!DEC$ IVDEP
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    IF(.NOT.LMSE)YDMF_PHYS%OUT%FRTHDS(JLON)=YDMF_PHYS%OUT%FRTH(JLON,YDCPG_OPTS%KFLEVG,1)/ZFLU_EMIS(JLON)+RSIGMA*YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)**4
    IF(ZRDG_MU0(JLON) <= 0.0_JPRB) THEN


      YDMF_PHYS%OUT%FRSOPT(JLON)=0.0_JPRB
    ELSE
      YDMF_PHYS%OUT%FRSOPT(JLON)=RII0*ZRDG_MU0(JLON)
    ENDIF
  ENDDO

! ADDITIONAL DIAGNOSTIC OF THE DERIVATIVE OF THE NON SOLAR SURFACE
! HEAT FLUX WITH RESPECT TO SURFACE TEMPERATURE (PDERNSHF)

  IF(LMCC03)THEN
    CALL ACDNSHF(YDCST, YDPHY, YDPHY1, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,          &
    & YDCPG_OPTS%KFLEVG, ZFLU_EMIS, YDMF_PHYS_SURF%GSD_VF%PLSM, ZFLU_NEIJ, ZQV, YDCPG_MISC%QS, YDMF_PHYS_BASE_STATE%YGSP_RR%T, &
    & ZCHROV, ZDQSTS, YDMF_PHYS%OUT%DRNSHF)
  ENDIF

!*
!     ------------------------------------------------------------------
!     9.- TRAINEE DES ONDES DE GRAVITE INDUITES PAR LE RELIEF
!     ------------------------------------------------------------------
  IF ( LGWD ) THEN
    CALL ACDRAG (YDCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTDRAG, YDCPG_OPTS%KFLEVG, &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,           &
    & ZNBVNO, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V,                               &
    & YDVARS%GEOMETRY%RCORI%T0, YDMF_PHYS_SURF%GSD_VF%PGETRL, ZGWDCS, YDMF_PHYS_SURF%GSD_VF%PVRLAN, YDMF_PHYS_SURF%GSD_VF%PVRLDI,    &
    & YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV, ZTRAJGWD)
  ENDIF
 ! SAVE FOR TL/NL COEFS FROM VERT. DIFF AND GWD

  IF (LTRAJPS.AND.(LVDIFSPNL.OR.LGWDSPNL)) THEN


     IF(.NOT.LVDIFSPNL) THEN
       ZKTROV_SAVE(:,:)=0.0_JPRB
       ZKUROV_SAVE(:,:)=0.0_JPRB
       ZCDROV_SAVE(:)=0.0_JPRB
       ZCHROV_SAVE(:)=0.0_JPRB
     ENDIF
     IF(.NOT. LGWDSPNL) THEN
       ZTRAJGWD(:,:)=0.0_JPRB
     ENDIF

     CALL WRPHTRAJTM_NL(YDGEOMETRY, YDSIMPHL, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, PTRAJ_PHYS,                        &
     & ZKTROV_SAVE, ZKUROV_SAVE, ZCDROV_SAVE, ZCHROV_SAVE, ZTRAJGWD, YDMF_PHYS_BASE_STATE%L, YDMF_PHYS_BASE_STATE%I, &
     & YDMF_PHYS_BASE_STATE%R, YDMF_PHYS_BASE_STATE%S, ZQLIS, ZNEBS)
  ENDIF

!     ------------------------------------------------------------------
!     10.- PRECIPITATIONS STRATIFORMES.
!     ---------------------------------

  IF ( LSTRA.AND.(.NOT.LSTRAPRO) ) THEN
    CALL ACPLUIE ( YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTPLUI,             &
    & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,           &
    & ZQV, ZMSC_QW, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FCSQL,     &
    & YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN&
    & )
  ENDIF

  IF ( LSTRAS ) THEN
    CALL ACPLUIS (YDCST,  YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTPLUI, &
    & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,             &
    & ZNEBS, ZQV, ZQLIS, ZMSC_QW, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS%OUT%FCSQL, &
    & YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN,  &
    & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN       )
  ENDIF

  IF (L3MT.OR.LSTRAPRO) THEN

    IF (LCOEFKTKE.OR.LLREDPR) THEN
      CALL ACTQSATS ( YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTQSAT, YDCPG_OPTS%KFLEVG, &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZT,                        &
      & ZGEOSLC, ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)
      CALL ACNEBCOND (YDCST, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,        &
      & NTPLUI, YDCPG_OPTS%KFLEVG, LLREDPR, YDSTA, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,    &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,    &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%RH, ZBLH, ZQV, ZQI, ZQL, ZQW, ZT, ZNEBCH,                         &
      & YDVARS%GEOMETRY%GM%T0, YDMF_PHYS_BASE_STATE%YGSP_RR%T, ZQLIS, ZNEBS, ZRHCRI, ZRH, ZQSATS,                            &
      & ZICEFR1, ZQLIS0, ZNEBS0)
      CALL ACCDEV (YDCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTPLUI,          &
      & YDCPG_OPTS%KFLEVG, YDSTA, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                   &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
      & ZQV, ZQI, ZQL, ZQS, ZQR, ZQG, ZQW, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, ZT, ZMSC_TW, ZNEBS,                          &
      & ZRHCRI, ZICEFR1, ZQSATS, ZNEBCH, ZRHDFDA, ZRKTTEND, ZRKQVTEND, ZRKQCTEND, ZQLI_CVPP,                                 &
      & ZNEB_CVPP, YDMF_PHYS_SURF%GSD_VF%PLSM, ZFLU_NEIJ, YDVARS%GEOMETRY%GM%T0, YDMF_PHYS_BASE_STATE%YGSP_RR%T,             &
      & ZDECRD_MF, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FCSQL,                    &
      & YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPLSL,       &
      & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG, ZSEDIQL, ZSEDIQI, YDMF_PHYS%OUT%DIAGH)
    ELSE
      CALL ACCDEV (YDCST,  YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTPLUI,         &
      & YDCPG_OPTS%KFLEVG, YDSTA, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                   &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
      & ZQV, ZQI, ZQL, ZQS, ZQR, ZQG, ZMSC_QW, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%T,                  &
      & ZMSC_TW, ZNEBS, ZRHCRI, ZICEFR1, ZQSATS, ZNEBCH, ZRHDFDA, ZRKTTEND, ZRKQVTEND, ZRKQCTEND,                            &
      & ZQLI_CVPP, ZNEB_CVPP, YDMF_PHYS_SURF%GSD_VF%PLSM, ZFLU_NEIJ, YDVARS%GEOMETRY%GM%T0, YDMF_PHYS_BASE_STATE%YGSP_RR%T,  &
      & ZDECRD_MF, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FCSQL,                    &
      & YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPLSL,       &
      & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG, ZSEDIQL, ZSEDIQI, YDMF_PHYS%OUT%DIAGH)
    ENDIF

  ENDIF

  IF (L3MT) THEN


    !  ----------------------------------------------------------------
    !   UPDATE HUMIDITY VARIABLES BY RESOLVED CONDENSATION FLUXES
    !  ----------------------------------------------------------------
!cdir unroll=8
    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZDQL=ZIPOI(JLON,JLEV)* (&
          & YDMF_PHYS%OUT%FCSQL(JLON,JLEV)-YDMF_PHYS%OUT%FCSQL(JLON,JLEV-1))
        ZQL(JLON,JLEV)=ZQL(JLON,JLEV)+ZDQL
        ZDQI=ZIPOI(JLON,JLEV)* (&
          & YDMF_PHYS%OUT%FCSQN(JLON,JLEV)-YDMF_PHYS%OUT%FCSQN(JLON,JLEV-1))
        ZQI(JLON,JLEV)=ZQI(JLON,JLEV)+ZDQI
        ZQV(JLON,JLEV)=ZQV(JLON,JLEV)-ZDQI-ZDQL

        ZT(JLON,JLEV)=ZT(JLON,JLEV)+(ZLHV(JLON,JLEV)*ZDQL&
                   & +ZLHS(JLON,JLEV)*ZDQI)/YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)
      ENDDO
    ENDDO

!rb : Store values for next time-step: all cases
    IF(LRKCDEV) THEN
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDVARS%RKTH%T0(JLON,JLEV)= ZT(JLON,JLEV)+ZDTRAD(JLON,JLEV)
          YDVARS%RKTQV%T0(JLON,JLEV)=ZQV(JLON,JLEV)+ZDQVDIFF(JLON,JLEV)*ZDTURDIFF
          YDVARS%RKTQC%T0(JLON,JLEV)= ZQI(JLON,JLEV)+ZQL(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ENDIF ! L3MT

!     11.- PRECIPITATIONS SOUS-MAILLES.
!     ---------------------------------

  IF (LVDIF.AND.(LCVRA.OR.L3MT.OR.LCVRAV3)) THEN
!           LA TENDANCE DYNAMIQUE DE Q EST MULTIPLIEE PAR UN
!           FACTEUR CORRECTIF FONCTION DE LA RESOLUTION LOCALE
!           DU MODELE PUIS AUGMENTEE DE LA CONTRIBUTION
!           DE L EVAPORATION DU SOL.
    IF(LCVCSD.OR..NOT.LCVGQD) THEN
      IF (LCOMOD) THEN
          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            ZMOD(JLON)=1._JPRB/(1.0_JPRB+YDVARS%GEOMETRY%GM%T0(JLON)*TEQK)
         ENDDO
      ELSE
         ZMOD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=1._JPRB
      ENDIF
!cdir unroll=8
      DO JLEV = YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG
        DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
          ZRDG_CVGQ(JLON,JLEV) = ZRDG_CVGQ(JLON,JLEV)*ZMOD(JLON) &
           & -RG*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)&
           & *(YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)-YDMF_PHYS%OUT%DIFTQ(JLON,JLEV-1)&
           & +ZFCQVNG(JLON,JLEV)-ZFCQVNG(JLON,JLEV-1))
        ENDDO
      ENDDO
    ENDIF
    IF (LCAPE) THEN
      IF (LCAMOD) THEN
         DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            ZTAUX(JLON)=RTCAPE*(YDVARS%GEOMETRY%GM%T0(JLON)*TEQC)
         ENDDO
      ELSE
         ZTAUX=RTCAPE
      ENDIF ! lcamod
    ENDIF
  ENDIF


  IF (L3MT) THEN
    ! - TEMPORAIRES
    ! ZSIGPC: CONVECTIVE FRACTION OF PRECIPITATION FLUX, USED A POSTERIORI 
    ! ZSIGP : PRECIPITATION MESH FRACTION
    ZSIGPC=0.0_JPRB
    ZSIGP=1.0_JPRB  ! used to limit/scale dd area
  IF (LCVPRO) THEN
  !  -------------------------
  !  UPDRAUGHT CONTRIBUTION
  !  -------------------------

    CALL ACTQSATS ( YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTQSAT, YDCPG_OPTS%KFLEVG, &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZT, ZGEOSLC,               &
    & ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)

    IF (.NOT.LCVCSD) THEN
      CALL ACCVUD ( YDMODEL%YRML_PHY_EC%YRECUCONVCA, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,             &
      & YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH,                     &
      & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,      &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
      & ZRDG_CVGQ, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR,                        &
      & ZQV, ZQI, ZQL, ZQSAT, ZQW, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, ZT, ZTW, YDMF_PHYS_BASE_STATE%U,              &
      & YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%VOR, YDVARS%DAL%T0, YDMF_PHYS_BASE_STATE%YGSP_RR%T,                  &
      & ZTAUX, YDVARS%GEOMETRY%RCORI%T0, YDVARS%GEOMETRY%GM%T0, ZDECRD_MF, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQL,     &
      & YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%STRCU,         &
      & YDMF_PHYS%OUT%STRCV, ZATSLC, INLAB, YDVARS%UNEBH%T0, ZDETFI, ZTU, ZQU, ZUU, ZVU, INND, YDMF_PHYS%OUT%CAPE,        &
      & YDMF_PHYS%OUT%CUCONVCA, YDMF_PHYS%OUT%NLCONVCA, ZGEOSLC, YDVARS%UEN%T0, YDVARS%UAL%T0, YDVARS%UOM%T0              &
      & )

    ! Initialize temperature correction
      ZTCORR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=ZTU(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)-ZT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
    ! Initialize Environment horizontal velocity
      ZU(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=YDMF_PHYS_BASE_STATE%U(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
      ZV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=YDMF_PHYS_BASE_STATE%V(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
    ELSE
! PUT IN ZTCORR THE CONDENSATE USED IN TRIGGERING
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
       DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZTCORR(JLON,JLEV)=MAX(0._JPRB,ZQI(JLON,JLEV)+ZQL(JLON,JLEV) -&
          &  MAX(0._JPRB,YDMF_PHYS_BASE_STATE%I(JLON,JLEV)+YDMF_PHYS_BASE_STATE%L(JLON,JLEV)) )
       ENDDO
      ENDDO
! PASS T9 CONDENSATES TO ACCSU
      CALL ACCSU ( YDMODEL%YRML_PHY_EC%YRECUCONVCA, YDMODEL%YRML_PHY_EC%YRECUMF, YDMODEL%YRML_PHY_MF,                                     &
      & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, &
      & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                      &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,                 &
      & ZRDG_CVGQ, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR,                                        &
      & ZQV, ZTCORR, YDMF_PHYS_BASE_STATE%I, YDMF_PHYS_BASE_STATE%L, ZQSAT, ZQW, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                 &
      & ZT, ZTW, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%VOR, YDCPG_DYN0%CTY%VVEL(:, 1:),                    &
      & YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS_BASE_STATE%YGSP_RR%T, ZTAUX, YDVARS%GEOMETRY%RCORI%T0,                        &
      & ZDECRD_MF, YDMF_PHYS_SURF%GSD_VK%PUDGRO, YDMF_PHYS%OUT%CUCONVCA, YDMF_PHYS%OUT%NLCONVCA, YDMF_PHYS%OUT%DIFCQ,                     &
      & YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN,                        &
      & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, ZATSLC, INLAB, YDVARS%UNEBH%T0, ZDETFI, ZTU, ZQU,                                       &
      & ZUU, ZVU, INND, YDMF_PHYS%OUT%CAPE, ZGEOSLC, YDVARS%UEN%T0, YDVARS%UAL%T0, YDVARS%UOM%T0)

    ! Initialize temperature correction
      ZTCORR(:,:)=ZTU(:,:)
    ! Initialize Environment horizontal velocity
      ZU(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=YDMF_PHYS_BASE_STATE%U(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
      ZV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=YDMF_PHYS_BASE_STATE%V(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
    ENDIF


    CALL ACUPU(YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,                             &
    & YDCPG_OPTS%KFLEVG, YDVARS%UAL%T0, YDVARS%UNEBH%T0, ZDETFI, ZPOID, ZIPOI, ZLHV, ZLHS, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, &
    & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL,                 &
    & YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, ZSIGP, ZSIGPC, ZNEBS, ZQI, ZQL,                             &
    & ZQV, ZT, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQLNG)

!rb store "-B" state: no protection of Ncv
    IF(LRKCDEV.AND.(.NOT.LNEBCV)) THEN
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDVARS%RKTH%T0(JLON,JLEV)=YDVARS%RKTH%T0(JLON,JLEV)-ZT(JLON,JLEV)
          YDVARS%RKTQV%T0(JLON,JLEV)=YDVARS%RKTQV%T0(JLON,JLEV)-ZQV(JLON,JLEV)
          YDVARS%RKTQC%T0(JLON,JLEV)=YDVARS%RKTQC%T0(JLON,JLEV)-ZQI(JLON,JLEV)&
        & -ZQL(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ENDIF !LCVPRO


  !  ---------------------------------------------------------
  !  MICROPHYSICS (AUTOCONVERSION, COLLECTION AND EVAPORATION)
  !        OF TOTAL CONDENSATE (RESOLVED + SUB-GRID)
  !  ---------------------------------------------------------

    ZAUXPRC(:)=0._JPRB
    ZRCVOTT(:,:)=0._JPRB

    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA

    !    COMPUTE TOTAL CONDENSATION FLUX
    !    -------------------------------
        ZFCQL(JLON,JLEV)=YDMF_PHYS%OUT%FCSQL(JLON,JLEV)+YDMF_PHYS%OUT%FCCQL(JLON,JLEV)
        ZFCQI(JLON,JLEV)=YDMF_PHYS%OUT%FCSQN(JLON,JLEV)+YDMF_PHYS%OUT%FCCQN(JLON,JLEV)
        IF (LRCVOTT) THEN
          ZCONVC=MAX(0.0_JPRB, YDMF_PHYS%OUT%FCCQL(JLON,JLEV)+YDMF_PHYS%OUT%FCCQN(JLON,JLEV))
          ZAUXPRC(JLON)=ZAUXPRC(JLON)+&
         & MAX(0.0_JPRB,YDMF_PHYS%OUT%FCSQL(JLON,JLEV)+YDMF_PHYS%OUT%FCSQN(JLON,JLEV)&
         & -YDMF_PHYS%OUT%FCSQL(JLON,JLEV-1)-YDMF_PHYS%OUT%FCSQN(JLON,JLEV-1))
          ZTOTC=ZCONVC+ZAUXPRC(JLON)
          ZRCVOTT(JLON,JLEV)=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZTOTC-ZEPS0))&
         & *ZCONVC/MAX(ZTOTC,ZEPS0)
        ENDIF
      ENDDO
    ENDDO

    CALL ACTQSATS ( YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTQSAT, YDCPG_OPTS%KFLEVG, &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZT, ZGEOSLC,               &
    & ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)

    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZDQ(JLON,JLEV)=ZQW(JLON,JLEV)&
         & -(ZQV(JLON,JLEV)+ZQI(JLON,JLEV)+ZQL(JLON,JLEV))
        ZDQM(JLON,JLEV)=ZQW(JLON,JLEV)
      ENDDO
    ENDDO

    CALL APLMPHYS( YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,          &
    & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,                   &
    & YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZQI, ZQL, ZQS, ZQR, ZQG, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,          &
    & ZT, ZIPOI, ZDQ, ZDQM, ZLHS, ZLHV, ZNEBS, ZPOID, ZRCVOTT, ZTCORR, YDMF_PHYS_SURF%GSD_VF%PLSM,                      &
    & ZFLU_NEIJ, YDMF_PHYS_BASE_STATE%YGSP_RR%T, ZDECRD_MF, ZFCQL, ZFCQI, ZFALLR, ZFALLS, ZFALLG, YDMF_PHYS%OUT%FPFPSL, &
    & YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG,  &
    & YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG, ZSEDIQL, ZSEDIQI, YDMF_PHYS%OUT%DIAGH)

    !  ------------------------------------------------------------
    !   UPDATE AFTER MICROPHYSICS - 3MT
    !  ------------------------------------------------------------

   CALL ACUPM( YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,             &
   & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZLHS, ZLHV, YDCPG_DYN0%CTY%EVEL, ZIPOI,                  &
   & ZPOID, YDVARS%UAL%T0, YDVARS%UOM%T0, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG,         &
   & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,       &
   & YDMF_PHYS%OUT%FPLSG, ZFCQL, ZFCQI, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, ZOME, ZFHP, ZQV,                     &
   & ZQL, ZQI, ZQR, ZQS, ZQG, ZT, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FCQRNG, &
   & YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG)

  IF (LCDDPRO) THEN
  !  -------------------------
  !  DOWNDRAUGHT CONTRIBUTION
  !  -------------------------

    ZDIFCQD(:,:)=0.0_JPRB
    ZDIFCQLD(:,:)=0.0_JPRB
    ZDIFCQID(:,:)=0.0_JPRB
    ZDIFCSD(:,:)=0.0_JPRB
    ZSTRCUD(:,:)=0.0_JPRB
    ZSTRCVD(:,:)=0.0_JPRB

    CALL ACTQSATS ( YDPHY, YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTQSAT, YDCPG_OPTS%KFLEVG, &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZT, ZGEOSLC,               &
    & ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)

    IF (LNSDO) THEN
      CALL ACNSDO(YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,                &
      & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,                       &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZQV, ZQI, ZQL, ZQR, ZQS, ZQW, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,       &
      & ZSIGP, ZT, ZTW, ZU, ZV, YDCPG_DYN0%CTY%VVEL(:, 1:), ZATSLC, ZGEOSLC, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,       &
      & YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, ZDIFCQD, ZDIFCQLD, ZDIFCQID, ZDIFCSD, YDMF_PHYS%OUT%FPEVPCL,           &
      & YDMF_PHYS%OUT%FPEVPCN, ZSTRCUD, ZSTRCVD, YDVARS%DAL%T0, YDVARS%DOM%T0)
    ELSE
      CALL ACMODO(YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,              &
      & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,                     &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, &
      & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZQV, ZQI,                          &
      & ZQL, ZQW, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                            &
      & ZSIGP, ZT, ZTW, ZU, ZV, YDCPG_DYN0%CTY%EVEL, ZOME, ZATSLC, ZGEOSLC, ZFHP, YDMF_PHYS%OUT%FPLSL,                     &
      & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG, ZDIFCQD, ZDIFCQLD, ZDIFCQID, ZDIFCSD, YDMF_PHYS%OUT%FPEVPCL,             &
      & YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPEVPCG, ZSTRCUD, ZSTRCVD, YDVARS%DAL%T0, YDVARS%DOM%T0                       &
      &                                                                                                                                                               )
    ENDIF
  !  ---------------------------------------------
  !  UPDATE VARIABLES  BY DOWNDRAUGHT CONTRIBUTION
  !  ---------------------------------------------


    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA

    !   UPDATE CONVECTIVE DIFFUSION AND EVAPORATION FLUXES
    !   --------------------------------------------------

        YDMF_PHYS%OUT%DIFCS(JLON,JLEV) =YDMF_PHYS%OUT%DIFCS(JLON,JLEV) +ZDIFCSD(JLON,JLEV)
        YDMF_PHYS%OUT%DIFCQ(JLON,JLEV) =YDMF_PHYS%OUT%DIFCQ(JLON,JLEV) +ZDIFCQD(JLON,JLEV)
        YDMF_PHYS%OUT%DIFCQL(JLON,JLEV)=YDMF_PHYS%OUT%DIFCQL(JLON,JLEV)+ZDIFCQLD(JLON,JLEV)
        YDMF_PHYS%OUT%DIFCQN(JLON,JLEV)=YDMF_PHYS%OUT%DIFCQN(JLON,JLEV)+ZDIFCQID(JLON,JLEV)
        YDMF_PHYS%OUT%STRCU(JLON,JLEV)=YDMF_PHYS%OUT%STRCU(JLON,JLEV)+ZSTRCUD(JLON,JLEV)
        YDMF_PHYS%OUT%STRCV(JLON,JLEV)=YDMF_PHYS%OUT%STRCV(JLON,JLEV)+ZSTRCVD(JLON,JLEV)
      ENDDO
    ENDDO

    CALL ACUPD(YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,               &
    & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP,                    &
    & ZLHS, ZLHV, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, ZIPOI, ZPOID, ZFALLR, ZFALLS, ZFALLG, YDMF_PHYS%OUT%FPEVPSL,      &
    & YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPEVPCG, &
    & ZFHP, ZDIFCSD, ZDIFCQD, ZDIFCQLD, ZDIFCQID, ZQV, ZQL, ZQI, ZQR, ZQS, ZQG, ZT, YDVARS%UEN%T0, YDMF_PHYS%OUT%FCQNG,  &
    & YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG,      &
    & YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG)

!rb  store "+C": case with no protection of Ncv
    IF((LRKCDEV).AND.(.NOT.LNEBCV)) THEN
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDVARS%RKTH%T0(JLON,JLEV)=YDVARS%RKTH%T0(JLON,JLEV)+ZT(JLON,JLEV)
          YDVARS%RKTQV%T0(JLON,JLEV)=YDVARS%RKTQV%T0(JLON,JLEV)+ZQV(JLON,JLEV)
          YDVARS%RKTQC%T0(JLON,JLEV)=YDVARS%RKTQC%T0(JLON,JLEV)&
          & +ZQI(JLON,JLEV)+ZQL(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  !  PARTITION CONVECTIVE/STRATIFORM PRECIPITATION
  ! ---------------------------------------------
    DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
         YDMF_PHYS%OUT%FPLCL(JLON,JLEV)=ZSIGPC(JLON)*YDMF_PHYS%OUT%FPLSL(JLON,JLEV)
         YDMF_PHYS%OUT%FPLSL(JLON,JLEV)=YDMF_PHYS%OUT%FPLSL(JLON,JLEV)-YDMF_PHYS%OUT%FPLCL(JLON,JLEV)
         YDMF_PHYS%OUT%FPLCN(JLON,JLEV)=ZSIGPC(JLON)*YDMF_PHYS%OUT%FPLSN(JLON,JLEV)
         YDMF_PHYS%OUT%FPLSN(JLON,JLEV)=YDMF_PHYS%OUT%FPLSN(JLON,JLEV)-YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
      ENDDO
    ENDDO
    IF (LGRAPRO) THEN
      DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
           YDMF_PHYS%OUT%FPLCG(JLON,YDCPG_OPTS%KFLEVG)=0.0_JPRB !ZSIGPC(JLON)*PFPLSG(JLON,JLEV)
           YDMF_PHYS%OUT%FPLSG(JLON,JLEV)=YDMF_PHYS%OUT%FPLSG(JLON,JLEV)-YDMF_PHYS%OUT%FPLCG(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF
  ENDIF ! LCDDPRO

  ENDIF ! L3MT

! -----------------------------------------------------

  IF ( LCVRA ) THEN

    IF (LSTRAPRO) THEN
      DO JLEV = YDCPG_OPTS%KTDIA-1, YDCPG_OPTS%KFLEVG
        DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
          YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)+ZFCQVNG(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

    CALL ACCVIMP ( YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTCVIM,                                    &
    & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP,                &
    & ZRDG_CVGQ, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZQV,                                   &
    & ZFLU_QSAT, ZMSC_QW, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                                 &
    & YDMF_PHYS_BASE_STATE%T, ZMSC_TW, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%VOR,                        &
    & YDCPG_DYN0%CTY%VVEL(:, 1:), YDVARS%GEOMETRY%GM%T0, YDMF_PHYS_BASE_STATE%YGSP_RR%T, ZTAUX, YDVARS%GEOMETRY%RCORI%T0,               &
    & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FPLCL,                          &
    & YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN,                    &
    & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, ZFPCOR, INLAB, YDMF_PHYS%OUT%CAPE, INND, YDMF_PHYS%OUT%DIFTQ,                           &
    & YDMF_PHYS%OUT%DIFTS, ZGEOSLC, YDVARS%GEOMETRY%GEMU%T0)

    IF (LSTRAPRO) THEN
      DO JLEV = YDCPG_OPTS%KTDIA-1, YDCPG_OPTS%KFLEVG
        DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
          YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)-ZFCQVNG(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ELSEIF (LCVRAV3) THEN

    CALL ACCVIMP_V3 (YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                                         &
    & NTCVIM, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI,                            &
    & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,                 &
    & YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZRDG_CVGQ, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, &
    & ZQV, ZMSC_QW, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_BASE_STATE%T,              &
    & ZMSC_TW, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VH%PPBLH,            &
    & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FPLCL,                         &
    & YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN,                   &
    & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, INLAB, INND, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS                                  &
    &                                                                                                                                                    )

  ELSEIF (LCVTDK) THEN     ! <== IFS deep convection scheme

    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZGEOM1(JLON,JLEV)     = YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(JLON,JLEV)-YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,YDCPG_OPTS%KFLEVG)
        ZVERVEL(JLON,JLEV)    = YDCPG_DYN0%CTY%VVEL(JLON,JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)
        ZLISUM(JLON,JLEV)     = 0.0_JPRB
        ZLCRIT_AER(JLON,JLEV) = 5.E-4_JPRB
      ENDDO
    ENDDO
    DO JLEV=0,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZGEOMH(JLON,JLEV)=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,JLEV)-YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,YDCPG_OPTS%KFLEVG)
        ZFCQLF(JLON,JLEV)=0.0_JPRB
        ZFCQLI(JLON,JLEV)=0.0_JPRB
      ENDDO
    ENDDO
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
       LLLAND(JLON)    = (YDMF_PHYS_SURF%GSD_VF%PLSM(JLON)>0.5_JPRB)
       LLSHCV(JLON)    = .FALSE.
       ZCUCONVCA(JLON) = 0.0_JPRB
       ZACPR(JLON)     = 0.0_JPRB
       ZDXTDK(JLON)    = RDELXN/YDVARS%GEOMETRY%GM%T0(JLON)
    ENDDO
    LLPTQ = .FALSE.
    CALL ACUPTQ (YDCST,  YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, LLPTQ, YDMF_PHYS%OUT%FRSO, &
    & YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS,                       &
    & YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPEVPCL,                    &
    & YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                   &
    & YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%Q, YDMF_PHYS_BASE_STATE%YGSP_RR%T, ZTENHA,                                       &
    & ZTENQVA )
    
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA 
         ZTENT(JLON,JLEV) = ZTENHA(JLON,JLEV)/YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)
         ZTENQ(JLON,JLEV) = ZTENQVA(JLON,JLEV)
         ZTENU(JLON,JLEV) = 0.0_JPRB
         ZTENV(JLON,JLEV) = 0.0_JPRB
         ZTENTA(JLON,JLEV) = 0.0_JPRB
         ZTENQA(JLON,JLEV) = 0.0_JPRB
      ENDDO
    ENDDO
    LLDSLPHY=.TRUE.
    ZVDIFTS = 0._JPRB
    ISPPN2D = 0
         
    CALL CUCALLN_MF (RPLDARE, RPLRG, YDCPG_OPTS%NSTEP, YDMODEL%YRML_PHY_EC%YRTHF, YDCST, YDMODEL%YRML_PHY_RAD%YRERAD, &
    & YDMODEL%YRML_PHY_SLIN, YDMODEL%YRML_PHY_EC, YDMODEL%YRML_GCONF%YGFL,                                                    &
    & YDMODEL%YRML_CHEM%YRCHEM, YDMODEL%YRML_GCONF%YRSPP_CONFIG, &
    & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, 0, YDCPG_OPTS%KFLEVG, ZDXTDK, ISPPN2D,                             &
    & YDMODEL%YRML_PHY_MF%YRPHY%YRCAPE%LMCAPEA, LLLAND, LLDSLPHY, TSPHY, ZVDIFTS, YDMF_PHYS_BASE_STATE%T, ZQV, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, &
    & ZLISUM, ZVERVEL, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                        &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, ZGEOM1, ZGEOMH, YDVARS%GEOMETRY%GM%T0,     &
    & ZCUCONVCA, ZGP2DSPPA, ZTENT, ZTENQ, ZTENU, ZTENV, ZTENTA, ZTENQA, ZACPR, ITOPC, IBASC, ITYPE, ICBOT, ICTOP,             &
    & IBOTSC, LLCUM, LLSC, ICBOT_LIG, ICTOP_LIG, LLCUM_LIG, &
    & LLSHCV, ZLCRIT_AER, ZLU, ZLUDE, ZLUDELI, ZSNDE, ZMFU, ZMFD, YDMF_PHYS%OUT%DIFCQ,                   &
    & YDMF_PHYS%OUT%DIFCS, ZFHPCL, ZFHPCN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, ZLRAIN, ZRSUD, YDMF_PHYS%OUT%STRCU,      &
    & YDMF_PHYS%OUT%STRCV, ZFCQLF, ZFCQLI, ZMFUDE_RATE, ZMFDDE_RATE, YDMF_PHYS%OUT%CAPE, ZWMEAN, ZVDISCU,                     &
    & ZDIFF, 0, ZCEN, ZTENC, ZSCAV)
    DO JLEV=0,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA 
         ZFPCOR  (JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
         YDMF_PHYS%OUT%FCCQL  (JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)
         YDMF_PHYS%OUT%FCCQN  (JLON,JLEV)=YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
         YDMF_PHYS%OUT%FPFPCL (JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)
         YDMF_PHYS%OUT%FPFPCN (JLON,JLEV)=YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
         YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV)=0._JPRB
         YDMF_PHYS%OUT%FPEVPCN(JLON,JLEV)=0._JPRB
       ENDDO
    ENDDO

  ENDIF
  IF(LFLASH .AND. (LCVTDK .OR. LGPCMT)) THEN
    ! LIGHTNING PARAMETERIZATION. OUTPUT IS PFLASH, TOTAL LIGHTNING FLASH RATES.
    ZGAW(:)=0._JPRB
    CALL ABOR1("Call to CULIGHT is not phased wrt 48R1 (ZQPFROZ, CHARGE)")
    !CALL CULIGHT (RPLDARE, RPLRG, YDMODEL%YRML_PHY_EC%YRTHF, YDCST, YDEPHY, YGFL, YDECUMF,&
    !& YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,                                                     &
    !& ZGAW, ZGAW, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, &
    !& YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, LLLAND, YDMF_PHYS_BASE_STATE%T, ZLU, ZMFU, YDMF_PHYS%OUT%CAPE,                          &
    !& YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, ZQPFROZ, LLCUM_LIG, ICBOT_LIG, ICTOP_LIG, LLLINOX, YDMF_PHYS%OUT%FLASH,           &
    !& ZLIGH_CTG, ZCTOPH, ZPRECMX, ZICE, ZCDEPTH, ZWMFU, YDMF_PHYS%OUT%CHARGE)
    ! LIGHTNING FLASH RATES ARE CONVERTED IN fl/km2/s BEFORE ENTERING CFU TIME ACCUMULATION.
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%FLASH(JLON)=YDMF_PHYS%OUT%FLASH(JLON)/86400._JPRB
    ENDDO
  ELSEIF(LFLASH) THEN
    DO JLON=YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%FLASH(JLON)=0._JPRB
    ENDDO
  ENDIF

  IF ( LADJCLD ) THEN
    ZFRSO(:,:) = 0.0_JPRB
    ZFRTH(:,:) = 0.0_JPRB
    LLPTQ = .TRUE.
    CALL ACUPTQ (YDCST,  YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, LLPTQ,   &
    & ZFRSO, ZFRTH, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS,           &
    & YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPEVPCL,  &
    & YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, &
    & YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%Q, YDMF_PHYS_BASE_STATE%YGSP_RR%T, ZTENHA,                     &
    & ZTENQVA )
  ENDIF

  IF ( LPROCLD ) THEN
    IF(LGPCMT) THEN
      ! Microphysics occurs in a shell between updraft and its resolved environment.
      ZNEBS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)=1._JPRB-(1._JPRB-ZNEBS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG))&
        & *(1._JPRB-ZCSGC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG))
      ZTMIC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%T(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)
      ZQMIC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)=ZCSGC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)&
        & *ZFLU_QSAT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)&
        & +(1._JPRB-ZCSGC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG))*ZQV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)
      ZQLIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)=MAX(ZQLIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG),&
        & ZCSGC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)&
        & *(ZQL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)+ZQI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG))&
        & +(1._JPRB-ZCSGC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG))*ZQLIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG))
      ZQC_DET_PCMT(:,:)=0._JPRB
    ELSE
      ! Microphysics occurs in the resolved state.
      ZTMIC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%T(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)
      ZQMIC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)=ZQV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)
    ENDIF

    CALL ACPLUIZ (YDCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTPLUI,                        &
    & YDCPG_OPTS%KFLEVG, ZTMIC, ZQMIC, ZQL, ZQI, ZQR, ZQS, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
    & YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, ZNEBS, ZQLIS,                                            &
    & ZNEB_CVPP, ZQLI_CVPP, ZQC_DET_PCMT, ZTENHA, ZTENQVA, LADJCLD, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI,                                    &
    & YDMF_PHYS_BASE_STATE%YGSP_RR%T, ZFLU_NEIJ, YDMF_PHYS_SURF%GSD_VF%PLSM, YDVARS%GEOMETRY%GM%T0,                                       &
    & YDSTA, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPEVPSL,                   &
    & YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, ZSEDIQL, ZSEDIQI )

    IF(LGPCMT.AND.LCVNHD) THEN
      ! Evaporation processes within convective environment
      ! will be used by the convection, to compute their feedback on convective updrafts. 
      ! This information is put here in PDDAL.
      YDVARS%DAL%T0(:,:)=0._JPRB
      DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDVARS%DAL%T0(JLON,JLEV)=( &
            & +YDMF_PHYS%OUT%FPEVPSL(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPSL(JLON,JLEV-1)   &
            & +YDMF_PHYS%OUT%FPEVPSN(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPSN(JLON,JLEV-1)   &
            & +YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV-1)   &
            & +YDMF_PHYS%OUT%FPEVPCN(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPCN(JLON,JLEV-1)) &
            & /RG*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)
        ENDDO
      ENDDO
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'PFPEVPSL0', YDMF_PHYS%OUT%FPEVPSL, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG+1&
                   &             )
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'PFPEVPSN0', YDMF_PHYS%OUT%FPEVPSN, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG+1&
                   &             )
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'PFPEVPCL0', YDMF_PHYS%OUT%FPEVPCL, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG+1&
                   &             )
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'PFPEVPCN0', YDMF_PHYS%OUT%FPEVPCN, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG+1&
                   &             )
    ENDIF
  ENDIF

  IF ( LNEBCO ) THEN
    CALL ACNEBC ( YDPHY0, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTNEBU, YDCPG_OPTS%KFLEVG,              &
    & YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, INLAB, INND, YDMF_PHYS_SURF%GSD_VH%PTCCH, YDMF_PHYS_SURF%GSD_VH%PSCCH, &
    & YDMF_PHYS_SURF%GSD_VH%PBCCH)
  ENDIF
  IF ( LGWDC ) THEN
    CALL ACDRAC ( YDPHY0, YDPHY2, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTDRAG, YDCPG_OPTS%KFLEVG,    &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, ZNBVNO, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_BASE_STATE%U, &
    & YDMF_PHYS_BASE_STATE%V, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS_SURF%GSD_VH%PSCCH,                 &
    & YDMF_PHYS_SURF%GSD_VH%PBCCH, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV)
  ENDIF

  ZFLX_LOTT_GWU(:, :) = 0.0_JPRB
  ZFLX_LOTT_GWV(:, :) = 0.0_JPRB

  IF ( LNORGWD ) THEN

    ! Inversion du sens des niveaux
    DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
      DO JLEV = 1, YDCPG_OPTS%KFLEVG
        Z_PP(JLON, JLEV) = YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON, YDCPG_OPTS%KFLEVG - JLEV + 1)
        Z_UU(JLON, JLEV) = YDMF_PHYS_BASE_STATE%U(JLON, YDCPG_OPTS%KFLEVG - JLEV + 1)
        Z_VV(JLON, JLEV) = YDMF_PHYS_BASE_STATE%V(JLON, YDCPG_OPTS%KFLEVG - JLEV + 1)
        Z_TT(JLON, JLEV) = YDMF_PHYS_BASE_STATE%T(JLON, YDCPG_OPTS%KFLEVG - JLEV + 1)
        Z_VO(JLON, JLEV) = YDMF_PHYS_BASE_STATE%VOR(JLON, YDCPG_OPTS%KFLEVG - JLEV + 1)
        ZD_U(JLON, JLEV) = TSPHY * YDMF_PHYS_BASE_STATE%P1NOGW(JLON, YDCPG_OPTS%KFLEVG - JLEV + 1)
        ZD_V(JLON, JLEV) = TSPHY * YDMF_PHYS_BASE_STATE%P2NOGW(JLON, YDCPG_OPTS%KFLEVG - JLEV + 1)
      ENDDO
    ENDDO

    ZPRECGWD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA) = MAX(0.0_JPRB,                  &
      & YDMF_PHYS%OUT%FPLCL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG) + YDMF_PHYS%OUT%FPLCN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG))
    
    CALL ACNORGWD(YDNORGWD, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, &
    & TSPHY, Z_PP, YDVARS%GEOMETRY%GEMU%T0, Z_TT, Z_UU, Z_VV, Z_VO, ZPRECGWD, ZD_U, ZD_V)
    
    DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
      DO JLEV = 1, YDCPG_OPTS%KFLEVG
        YDMF_PHYS_BASE_STATE%P1NOGW(JLON, JLEV) = ZD_U(JLON, YDCPG_OPTS%KFLEVG - JLEV + 1) / TSPHY
        YDMF_PHYS_BASE_STATE%P2NOGW(JLON, JLEV) = ZD_V(JLON, YDCPG_OPTS%KFLEVG - JLEV + 1) / TSPHY
      ENDDO
    ENDDO

    !-- CALCUL DU FLUX, PAR INTEGRATION DE LA TENDANCE DE HAUT EN BAS.
    !   LES FLUX SONT SUPPOSES NULS AU PREMIER NIVEAU (1 = TOP) DE CALCUL.
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
       DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
          ZFLX_LOTT_GWU(JLON, JLEV) = ZFLX_LOTT_GWU(JLON, JLEV - 1) - YDMF_PHYS_BASE_STATE%P1NOGW(JLON, JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON, JLEV)/RG
          ZFLX_LOTT_GWV(JLON, JLEV) = ZFLX_LOTT_GWV(JLON, JLEV - 1) - YDMF_PHYS_BASE_STATE%P2NOGW(JLON, JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON, JLEV)/RG
          YDMF_PHYS%OUT%STRCU(JLON, JLEV) = YDMF_PHYS%OUT%STRCU(JLON, JLEV) + ZFLX_LOTT_GWU(JLON, JLEV)
          YDMF_PHYS%OUT%STRCV(JLON, JLEV) = YDMF_PHYS%OUT%STRCV(JLON, JLEV) + ZFLX_LOTT_GWV(JLON, JLEV)          
    !+ TODOLATER +!      PSTRNORGWDU(JLON, JLEV) = PSTRNORGWDU(JLON, JLEV) + ZFLX_LOTT_GWU(JLON, JLEV)
    !+ TODOLATER +!      PSTRNORGWDV(JLON, JLEV) = PSTRNORGWDV(JLON, JLEV) + ZFLX_LOTT_GWV(JLON, JLEV)          
       ENDDO
    ENDDO

    ! NO VERTICAL DIFFUSION IN THE MIDDLE ATMOSPHERE
    DO JLEV = 1, NORGWD_NNOVERDIF
       DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
          YDMF_PHYS%OUT%STRTU(JLON, JLEV) = 0.0
          YDMF_PHYS%OUT%STRTV(JLON, JLEV) = 0.0          
       ENDDO
    ENDDO

    !+ TODOLATER +! IF (LNOWINDTEND) THEN
    !+ TODOLATER +!  DO JLEV = 0, KLEV
    !+ TODOLATER +!    DO JLON = KIDIA, KFDIA
    !+ TODOLATER +!      PSTRCU(JLON,JLEV) = 0._JPRB
    !+ TODOLATER +!      PSTRCV(JLON,JLEV) = 0._JPRB
    !+ TODOLATER +!    ENDDO
    !+ TODOLATER +!  ENDDO
    !+ TODOLATER +! ENDIF

  ENDIF

!*
!     ------------------------------------------------------------------
!         SAUVEGARDE DES FLUX DE PRECIPITATION CONVECTIVE ET STRATIFORME.

  IF(LNEBN.OR.LNEBR.OR.LRRGUST) THEN
    IF (LFPCOR) THEN
      IF (LCDDPRO) THEN
        DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            ZPFL_FPLCH(JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
          ENDDO
        ENDDO
      ELSE
        DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            ZPFL_FPLCH(JLON,JLEV)=ZFPCOR(JLON,JLEV)
          ENDDO
        ENDDO
      ENDIF
    ELSE
      DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZPFL_FPLCH(JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  IF(LRRGUST) THEN
    DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZPFL_FPLSH(JLON,JLEV)=YDMF_PHYS%OUT%FPLSL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSG(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! UPDATE TRANSPORT FLUXES DUE TO SEDIMENTATION OF CLOUDS.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)+ZSEDIQL(JLON,JLEV)
      YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)+ZSEDIQI(JLON,JLEV)
      ZFPLSL (JLON,JLEV)=YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSL (JLON,JLEV)
      ZFPLSN (JLON,JLEV)=YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN (JLON,JLEV)
    ENDDO
  ENDDO
  IF (LGRAPRO) THEN
    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZFPLSN (JLON,JLEV)=ZFPLSN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSG(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF

  ! - - - - - - - - - - - - - - - - -
  ! CORRECT NEGATIVE WATER CONTENTS.
  ! - - - - - - - - - - - - - - - - -

  IF ( LCONDWT.AND.LPROCLD.AND..NOT.LGPCMT ) THEN
    CALL QNGCOR (YDCST, YDPHY2, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTPLUI, YDCPG_OPTS%KFLEVG,&
    & ZQV, ZQL, ZQI, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN,               &
    & YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%FPEVPSL,    &
    & YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, &
    & YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL,       &
    & YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FCQNG )

    ! Comment the following lines, which generate negative QR/QS values in CPTEND_NEW.
    !DO JLON=KIDIA,KFDIA
    !  PFPLSL(JLON,KLEV)=PFPLSL(JLON,KLEV)+PDIFTQL(JLON,KLEV)
    !  PFPLSN(JLON,KLEV)=PFPLSN(JLON,KLEV)+PDIFTQI(JLON,KLEV)
    !ENDDO

  ENDIF

!*
!     ------------------------------------------------------------------
!     12. - BILAN HYDRIQUE DU SOL
!     ---------------------------
  IF ( LMSE ) THEN

    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDCPG_GPAR%RAIN(JLON)=ZFPLSL(JLON,YDCPG_OPTS%KFLEVG)+YDMF_PHYS%OUT%FPLCL(JLON,YDCPG_OPTS%KFLEVG)
      YDCPG_GPAR%SNOW(JLON)=ZFPLSN(JLON,YDCPG_OPTS%KFLEVG)+YDMF_PHYS%OUT%FPLCN(JLON,YDCPG_OPTS%KFLEVG)
    ENDDO

  ELSE

    IF ( LSFHYD.AND.LSOLV ) THEN
      CALL ACDROV (YDCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,&
      & YSP_SBD%NLEVS, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, ZFPLSL, ZFPLSN, YDMF_PHYS%OUT%FRSO,                          &
      & YDMF_PHYS%OUT%FRTH, ZDSA_C1, ZDSA_C2, ZC3, ZCN, YDMF_PHYS%OUT%CT, YDMF_PHYS_SURF%GSD_VV%PD2,                          &
      & YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FTR, YDMF_PHYS_SURF%GSD_VV%PLAI, ZFLU_NEIJ, ZFLU_VEG, ZWFC,                         &
      & ZWPMX, YDMF_PHYS_BASE_STATE%YGSP_RR%FC, ZWLMX, ZWSEQ, ZWSMX, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL,                 &
      & YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS, ZFLU_FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS_SURF%GSD_VF%PLSM, &
      & YDMF_PHYS_BASE_STATE%YGSP_SG%F, YDMF_PHYS_BASE_STATE%YGSP_SB%T, YDMF_PHYS_BASE_STATE%YGSP_RR%T,                       &
      & YDMF_PHYS_BASE_STATE%YGSP_SB%Q, YDMF_PHYS_BASE_STATE%YGSP_SB%TL, YDMF_PHYS_BASE_STATE%YGSP_RR%W,                      &
      & YDMF_PHYS_BASE_STATE%YGSP_RR%IC, YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS, YDMF_PHYS%OUT%FLWSP,                        &
      & YDMF_PHYS%OUT%FONTE, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%RUISS)
    ENDIF

  ENDIF

!*
!-  --------------------------------------------------------------------
!     13.- DRAG MESOSPHERIQUE POUR UN MODELE POSSEDANT DES NIVEAUX
!               AU-DESSUS DE 50 KM (I.E. DANS LA MESOSPHERE)
!     ------------------------------------------------------------------
  IF ( LRRMES ) THEN
    CALL ACDRME (YDCST, YDSTA, YDPHY2, YDTOPH, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,    &
    & NTDRME, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
    & YDMF_PHYS_BASE_STATE%T, ZQV, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, YDMF_PHYS%OUT%FRMH,         &
    & ZMSC_FRMQ, YDMF_PHYS%OUT%STRMU, YDMF_PHYS%OUT%STRMV)
  ENDIF

!*
!     ------------------------------------------------------------------
!     14.- FLUX PHOTO-CHIMIQUE D'OZONE
!     ------------------------------------------------------------------
  IF ( LOZONE ) THEN
    CALL ACOZONE ( YDPHY2, YDTOPH, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTOZON, YDCPG_OPTS%KFLEVG,   &
    & YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%KOZO, ZQO3(1, 1), YDMF_PHYS_BASE_STATE%T,                   &
    & ZRDG_MU0, YDMF_PHYS%OUT%FCHOZ)

!           DIFFUSION TURBULENTE/DEPOT SEC DE L'OZONE

    CALL ACDIFOZ ( YDRIP, YDPHY2, YDMODEL%YRML_PHY_MF%YRVDOZ, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,               &
    & YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,        &
    & YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FRSO, ZKTROV, ZQO3(1, 1), YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, &
    & YDMF_PHYS_BASE_STATE%T, ZXTROV, ZFLU_NEIJ, YDVARS%GEOMETRY%GEMU%T0, YDMF_PHYS_SURF%GSD_VV%PIVEG,                           &
    & YDMF_PHYS%OUT%FCHOZ                                                                                                        &
    &                                                                                                       )
  ENDIF
!*
!     ------------------------------------------------------------------
!     15.- CALCUL DES FLUX D'ENTHALPIE ET DE CHALEUR SENSIBLE LIES AUX
!          PRECIPITATIONS EN FONCTION DES FLUX DE PRECIPITATION
!          ET DE CONDENSATION.
!     ------------------------------------------------------------------

  ! STORE THE PSEUDO-HISTORIC SURFACE PRECIPITATION SENSIBLE HEAT FLUX
  ! -------------------------------------------------------------------

  IF (LPHSPSH) THEN
    DO JLON=YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
       YDMF_PHYS_SURF%GSD_VH%PSPSH(JLON)=(YDMF_PHYS_BASE_STATE%T(JLON,YDCPG_OPTS%KFLEVG)-YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)) * ((&
          & RCW-RCPD)*(YDMF_PHYS%OUT%FPLSL(JLON,YDCPG_OPTS%KFLEVG)+YDMF_PHYS%OUT%FPLCL(JLON,YDCPG_OPTS%KFLEVG))&
          & +(RCS-RCPD)*(YDMF_PHYS%OUT%FPLSN(JLON,YDCPG_OPTS%KFLEVG)+YDMF_PHYS%OUT%FPLCN(JLON,YDCPG_OPTS%KFLEVG)) )
    ENDDO
    IF (LGRAPRO) THEN
      DO JLON=YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
        YDMF_PHYS_SURF%GSD_VH%PSPSH(JLON)=YDMF_PHYS_SURF%GSD_VH%PSPSH(JLON)+(YDMF_PHYS_BASE_STATE%T(JLON,YDCPG_OPTS%KFLEVG)-YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)) * (&
        &(RCS-RCPD)*(YDMF_PHYS%OUT%FPLSG(JLON,YDCPG_OPTS%KFLEVG)+YDMF_PHYS%OUT%FPLCG(JLON,YDCPG_OPTS%KFLEVG)))
      ENDDO
    ENDIF
  ENDIF !LPHSPSH
     
!*
!     ------------------------------------------------------------------
!     16.- CALCUL DE LA DENSITE DES FOUDRES (DIAGNOSTIQUE).
!     ------------------------------------------------------------------

  IF (LFLASH) THEN
! tbd
! about vertical velocity:
! PWW is the resolved w in m/s given by the dynamics.
! In case moist deep convection is not fully resolved, we should consider
! computing w from the mass flux. See transport in ACCVUD.
!
    CALL DIAGFLASH(YDCFU,YDMODEL%YRML_PHY_MF,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,&
         & YDCPG_OPTS%NSTEP,&
         & ZQL,ZQI,YDMF_PHYS_BASE_STATE%R,YDMF_PHYS_BASE_STATE%S,ZQG,ZQH,YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,&
         & YDMF_PHYS_BASE_STATE%T,YDMF_PHYS_BASE_STATE%YCPG_PHY%W,YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,&
         & YDVARS%UAL%T0,YDVARS%UOM%T0,YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,&
         & YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS%OUT%FLASH)
  ENDIF

! Store radiative cloudiness in GFL structure for ISP, Historical files or PostProcessing
  IF (YIRAD%LGP) YDVARS%IRAD%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = YDCPG_MISC%QICE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
  IF (YLRAD%LGP) YDVARS%LRAD%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = YDCPG_MISC%QLI (YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
  IF (YA%LGP)    YDVARS%A%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = YDCPG_MISC%NEB (YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)


!*
!     -----------------------------------------------------------------------
!     16.- DDH FLEXIBLES POUR LES CHAMPS DE SURFACE VARIABLES/FLUX/TENDANCES.
!     -----------------------------------------------------------------------

  IF (YDMODEL%YRML_DIAG%YRLDDH%LFLEXDIA) THEN
    CALL APLPAR_FLEXDIA (YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, YDMF_PHYS, YDMF_PHYS_SURF, YDMODEL, &
    & YDDDH, YDMF_PHYS_BASE_STATE)
  ENDIF

  IF (LECT.AND.LRAFTKE) THEN
     ! DCAPE due to precipitation evaporation.
     CALL ACEVADCAPE(YDMODEL%YRML_PHY_MF%YRPHY2, YDCST, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMF_PHYS%OUT%FPLSL, &
     & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                      &
     & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP,                      &
     & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, ZDCAPE)
     ! Gusts.
     CALL ACCLDIA(YDCST, YDCPG_OPTS%LXCLP, YDCPG_OPTS%LXTGST, YDCPG_OPTS%LXXGST, YDPHY, YDPHY2, YDTOPH,        &
     & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, &
     & YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V,  YDMF_PHYS%OUT%CAPE, ZDCAPE, ZTKE1,                             &
     & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST, ZBLH, IMOC_CLPH)
     YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)))
  ENDIF

  IF (LXVISI.OR.LXVISI2) THEN
     ZQGM(:,:)=ZEPSNEB
     CALL ACVISIH(YDCST, YDMODEL%YRML_PHY_MF%YRPHY%YRDVISI, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,  &
     & YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,        &
     & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,                &
     & YDCPG_MISC%QLI, YDCPG_MISC%QICE, ZQR, ZQS, ZQGM, YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD,                     &
     & YDMF_PHYS%OUT%MXCLWC    )
  ENDIF

  !-------------------------------------------------
  ! Check if fluxes are different from NaN.
  !-------------------------------------------------
  IF (YDPHY%LGCHECKNAN) THEN
    CALL CHECKNAN(YDRIP,YDPHY0,YDPHY2, YDCPG_OPTS%NINDAT, &
    & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDCPG_OPTS%NSTEP, &
    & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, &
    & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
    & YDVARS%GEOMETRY%GELAM%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDVARS%GEOMETRY%GEMU%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), &
    & ZRDG_MU0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS_SURF%GSD_VF%PLSM, &
    & YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%Q, ZFLU_QSAT, YDMF_PHYS_BASE_STATE%YGSP_RR%T, &
    & YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC, YDMF_PHYS%OUT%FIMCC, &
    & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, &
    & YDMF_PHYS%OUT%FCNEGQLC, YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC, YDMF_PHYS%OUT%FCNEGQSC, &
    & ZPFL_FPLCH, &
    & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, &
    & YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS, &
    & YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FCQNG, &
    & YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, &
    & YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG, &
    & YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG, &
    & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV, &
    & YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV, YDMF_PHYS%OUT%STRMU, YDMF_PHYS%OUT%STRMV, &
    & YDMF_PHYS%OUT%FRMH, ZMSC_FRMQ, YDMF_PHYS%OUT%FCHOZ, &
    & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, &
    & YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN)
  ENDIF
ENDIF !LMPHYS
!----------------------------------------------------
!  CALCUL DE DEPOT HUMIDE POUR LES AEROSOLS DESERTIQUES
!----------------------------------------------------
IF (LMDUST.AND.(NGFL_EXT/=0).AND.LRDEPOS) THEN

IKRR=6
ISPLITR=1
ZEVAP=0.0_JPRB
ZZDEP=0.0_JPRB

  DO JLEV = 1 , YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZZDEP(JLON,1,JLEV)=ZZI_APHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO


  DO JLEV=1,YDCPG_OPTS%KFLEVG
    DO JLON= YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
      ZQDM(JLON,JLEV)=1._JPRB-ZQV(JLON,JLEV)-ZQL(JLON,JLEV)-ZQR(JLON,JLEV) &
       & -ZQI(JLON,JLEV)-ZQS(JLON,JLEV)
    ENDDO
  ENDDO

  DO JLEV = 1, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZZI_RM(JLON,1,JLEV,2)=ZQL(JLON,JLEV)&
       & /ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  DO JLEV = 1, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZZI_RM(JLON,1,JLEV,3)=ZQR(JLON,JLEV)&
       & /ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  DO JLEV = 1, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZEVAP(JLON,1,JLEV)=YDMF_PHYS%OUT%FPEVPSL(JLON,JLEV)&
       & +YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV)
    ENDDO
  ENDDO

  DO JGFL=1,NGFL_EXT
    DO JLON=1,YDCPG_OPTS%KLON
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        ZZI_SVM(JLON,1,JLEV,JGFL)=MAX(0._JPRB,ZSVM(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO

    CALL ARO_WETDEP(ILONMNH, YDCPG_OPTS%KFLEVG, NGFL_EXT, IKRR, YDCPG_OPTS%ZDTPHY, ZZI_SVM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :, 1:NGFL_EXT), &
    & ZZDEP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :), ZZI_PABSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :),                             &
    & ZZI_THM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :), ZZI_RHODREFM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :),                        &
    & YDCPG_OPTS%NSTEP+1, ZZI_RM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :, :), ZEVAP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :, :),         &
    & ISPLITR               )
! return to tendency
   DO JGFL=1,NGFL_EXT
     DO JLON=1,YDCPG_OPTS%KLON
       DO JLEV=1,YDCPG_OPTS%KFLEVG
         ZPSV(JLON,JLEV,JGFL)=ZZI_SVM(JLON,1,JLEV,JGFL)
       ENDDO
     ENDDO
   ENDDO

      DO JGFL=1,NGFL_EXT
        DO JLON=1,YDCPG_OPTS%KLON
          DO JLEV=1,YDCPG_OPTS%KFLEVG
            ZTENDEXT_DEP(JLON,JLEV,JGFL)=(ZPSV(JLON,JLEV,JGFL)-ZSVM(JLON,JLEV,JGFL))*ZINVDT
          ENDDO
        ENDDO
      ENDDO

ENDIF ! ENDIF  (LMDUST & NGFL_EXT & LRDEPOS)

IF(LMUSCLFA) THEN
  DO JLEV=0,YDCPG_OPTS%KFLEVG
    DO JLON=1,YDCPG_OPTS%KLON
      ZFPLS(JLON,JLEV)=YDMF_PHYS%OUT%FPLSL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN(JLON,JLEV)
      ZFPLC(JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
      ZFPL (JLON,JLEV)=ZFPLC (JLON,JLEV)+ZFPLS (JLON,JLEV)
    ENDDO
  ENDDO
  CALL WRSCMR(NMUSCLFA, 'ZFPLS', ZFPLS, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG+1)
  CALL WRSCMR(NMUSCLFA, 'ZFPLC', ZFPLC, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG+1)
  CALL WRSCMR(NMUSCLFA, 'ZFPL', ZFPL, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG+1)
ENDIF

!     ------------------------------------------------------------------

! BAYRAD
! Compute convective hydrometeors mixing ratio from diagnistic fluxes
!--------------------------------------------------------------------
IF ( (.NOT. LGPCMT) ) THEN

   ! Convert from flux [kg/m2/s] to density [kg/m3] using old RTTOV-SCATT
   !  a    b         ! RR = a * LWC^b, [RR]=mm/h, [LWC]=g/m^3
   ! 20.89 1.15      ! rain
   ! 29.51 1.10      ! snow

   ZBAY_QRCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) =  0.0_JPRB
   ZBAY_QSCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) =  0.0_JPRB


   DO JLEV=1,YDCPG_OPTS%KFLEVG
     DO JLON= YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
       ZBAY_QRCONV(JLON,JLEV)  = 0.001_JPRB * ( (ABS(YDMF_PHYS%OUT%FPLCL(JLON,JLEV)) * YDPHY0%RCOEFRAIN(1)) ** YDPHY0%RCOEFRAIN(2) )
       ZBAY_QSCONV(JLON,JLEV)  = 0.001_JPRB * ( (ABS(YDMF_PHYS%OUT%FPLCN(JLON,JLEV)) * YDPHY0%RCOEFSNOW(1)) ** YDPHY0%RCOEFSNOW(2) )    
     ENDDO
   ENDDO


   ! Convert density [kg/m3] to mixing ratio [kg/kg]
   ! R_dry (dry air constant)

   ZDE2MR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)  = RD * YDMF_PHYS_BASE_STATE%T(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) / YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
   ZBAY_QRCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = ZBAY_QRCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) * ZDE2MR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
   ZBAY_QSCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = ZBAY_QSCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) * ZDE2MR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)

ENDIF


! Precipitation Type

! Compute wet-bulb temperature at 2 meters (suppose homogeneity of qv/ql/qi )
!ZPCLS(KIDIA:KFDIA)=PAPRS(KIDIA:KFDIA,KLEV)-2._JPRB/ZZZF(KIDIA:KFDIA,1,KLEV)*&
!                 &(PAPRS(KIDIA:KFDIA,KLEV)-PAPRSF(KIDIA:KFDIA,KLEV))
CALL PPWETPOINT(YDCST, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, ZPCLS, YDMF_PHYS%OUT%TCLS,                       &
& YDMF_PHYS_BASE_STATE%Q(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS_BASE_STATE%L(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS_BASE_STATE%I(:, YDCPG_OPTS%KFLEVG), &
& YDMF_PHYS%OUT%TPWCLS                                                                                                                      &
&                                                                                                                                                                                                                        )


IF (LDPRECIPS .OR. LDPRECIPS2) THEN
   ! Defined precipitation type
   !
   ! Compute wet-bulb temperature
  DO JLEV=1,YDCPG_OPTS%KFLEVG
    CALL PPWETPOINT(YDCST, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(:, JLEV), &
    & YDMF_PHYS_BASE_STATE%T(:, JLEV), YDMF_PHYS_BASE_STATE%Q(:, JLEV), YDMF_PHYS_BASE_STATE%L(:, JLEV),                                       &
    & YDMF_PHYS_BASE_STATE%I(:, JLEV), ZTW(:, JLEV))
  ENDDO

  DO JLON=1,YDCPG_OPTS%KLON
      ZFPLS(JLON,YDCPG_OPTS%KFLEVG)=YDMF_PHYS%OUT%FPLCN(JLON,YDCPG_OPTS%KFLEVG)+YDMF_PHYS%OUT%FPLSN(JLON,YDCPG_OPTS%KFLEVG)
      ZFPLC(JLON,YDCPG_OPTS%KFLEVG)=YDMF_PHYS%OUT%FPLCL(JLON,YDCPG_OPTS%KFLEVG)+YDMF_PHYS%OUT%FPLSL(JLON,YDCPG_OPTS%KFLEVG)
      ZFPLSG(JLON,YDCPG_OPTS%KFLEVG)=0._JPRB
  ENDDO

  !initialisation de ZZZ
  DO JLEV = 1,YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZZZ(JLON,1,JLEV)=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO

  !initialisation de ZDZZ
  DO JLEV = 2, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZDZZ(JLON,1,JLEV)=ZZZ(JLON,1,JLEV-1)-ZZZ(JLON,1,JLEV)
    ENDDO
  ENDDO
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZDZZ(JLON,1,1)=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,0)*ZINVG-ZZZ(JLON,1,1)
  ENDDO

  IF (LDPRECIPS) THEN

   CALL DPRECIPS (YDCST, YDPHY%YRDPRECIPS, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,      &
   & YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,                    &
   & ZDZZ, ZTW, YDMF_PHYS_BASE_STATE%L, ZFPLC(:, YDCPG_OPTS%KFLEVG), ZFPLS(:, YDCPG_OPTS%KFLEVG), ZFPLSG(:, YDCPG_OPTS%KFLEVG), &
   & ZPRC_DPRECIPS(:, YDCPG_OPTS%NDTPRECCUR))

  ENDIF

  IF (LDPRECIPS2) THEN

 !Idem for an other time step and an other period

   CALL DPRECIPS(YDCST, YDPHY%YRDPRECIPS, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,       &
   & YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,                    &
   & ZDZZ, ZTW, YDMF_PHYS_BASE_STATE%L, ZFPLC(:, YDCPG_OPTS%KFLEVG), ZFPLS(:, YDCPG_OPTS%KFLEVG), ZFPLSG(:, YDCPG_OPTS%KFLEVG), &
   & ZPRC_DPRECIPS2(:, YDCPG_OPTS%NDTPRECCUR2))

  ENDIF
ENDIF  

IF (LAJUCV) THEN
  DO JLEV=1,YDCPG_OPTS%KFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDVARS%T%T0(JROF,JLEV)=ZADJ_TAUX(JROF,JLEV)
    ENDDO
  ENDDO
  DO JLEV=1,YDCPG_OPTS%KFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS_NEXT_STATE%T%P (JROF, JLEV) = YDMF_PHYS_NEXT_STATE%T%P (JROF, JLEV) + ZADJ_DTAJU(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!    convert to flexible interface structure
IF (LINTFLEX) THEN
  CALL APLPAR2INTFLEX(YGFL, YDPHY, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG,              &
  & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS, ZDIFEXT, YDMF_PHYS%OUT%DIFTQ, &
  & YDMF_PHYS%OUT    %DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN,      &
  & YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLCL,            &
  & YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN,    &
  & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FCQLNG,       &
  & YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FRMH,          &
  & ZMSC_FRMQ, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, YDMF_PHYS%OUT%STRDU,   &
  & YDMF_PHYS%OUT%STRDV, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV, YDMF_PHYS%OUT%STRMU, YDMF_PHYS%OUT%STRMV,            &
  & YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC, YDMF_PHYS%OUT%FIMCC, YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC,      &
  & YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC, YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC, &
  & YDMF_PHYS%OUT%FCNEGQSC, ZPFL_FTKE, ZTENDPTKE, ZTENDEXT, ZTENDEXT_DEP, YLPROCSET )
ENDIF

!        2.3  Computes MOCON in the CLP.
!             --------------------------
CALL MF_PHYS_MOCON (YDCPG_BNDS, YDCPG_OPTS, ZRDG_LCVQ, IMOC_CLPH, YDMF_PHYS, YDMF_PHYS_BASE_STATE)

! Store surface water flux P and E for water conservation
IF (LCORWAT) THEN
  CALL MF_PHYS_CORWAT (YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%FPLCL,                &
  & YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS_SURF%GSD_VH%PEVA, YDMF_PHYS_SURF%GSD_VH%PPCL, &
  & YDMF_PHYS_SURF%GSD_VH%PPCN, YDMF_PHYS_SURF%GSD_VH%PPSL, YDMF_PHYS_SURF%GSD_VH%PPSN)
ENDIF

!        2.4  Stores radiation coefficients.
!             ------------------------------

! * writes grid-point transmission coefficients for simplified physics.


IF (LRCOEF.AND.(YDCPG_OPTS%NSTEP == 1)) THEN
    IFIELDSS=NG3SR*YDCPG_OPTS%KFLEVG
    CALL WRRADCOEF(YDGEOMETRY, YDRCOEF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_BNDS%KSTGLO,        &
    & IFIELDSS, ZRDT_COR, ZRDT_RAB3C, ZRDT_RAB3N, ZRDT_RAB4C, ZRDT_RAB4N, ZRDT_RAB6C, ZRDT_RAB6N,     &
    & ZRDT_RAT1C, ZRDT_RAT1N, ZRDT_RAT2C, ZRDT_RAT2N, ZRDT_RAT3C, ZRDT_RAT3N, ZRDT_RAT4C, ZRDT_RAT4N, &
    & ZRDT_RAT5C, ZRDT_RAT5N, ZAC_HC)
ENDIF

!       2.5   Ozone
!             -----


IF (LOZONE) THEN
  ! * Caution: this part has not been yet validated relative
  !   to the GFL implementation, and LOZONE (the setup of
  !   which has not yet been updated) can be true only if
  !   the GFL ozone is activated as a prognostic and advected
  !   variable.
  CALL CPOZO (YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, YDCPG_OPTS%ZDTPHY, YDMF_PHYS%OUT%FCHOZ, &
  & YDMF_PHYS_NEXT_STATE%O3%P (:, 1:YDCPG_OPTS%KFLEVG), YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP)
ENDIF
  
!        2.5.1 Chemical species   
!              ----------------
IF (LCHEM_ARPCLIM) THEN
   ! Processes described in my_phys ARPEGE-Climat 6.3 : to be added later here
   ! Modify also calls in CPTEND_NEW, etc.. as done ARPEGE-Climat 6.3.   
ENDIF   

!        2.6   surface specific humidity necessary to compute the vertical
!              advection of q in the case "delta m=1" (unlagged physics only).
!              ---------------------------------------------------------------

IF (NDPSFI == 1) THEN
  CALL CPQSOL(YDCST, YDGEOMETRY%YRDIMV, YDPHY, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_PHY0%PREHYD, &
  & YDMF_PHYS_SURF%GSP_RR%PT_T0, YDCPG_MISC%QS, ZFLU_QSATS, YDCPG_MISC%QSOL)
ENDIF

!        2.7  Computation of tendencies T,u,v and Q.
!             --------------------------------------

! Set GFL tendencies to 0

ZTENDGFL(:,:,:) = 0.0_JPRB

! * CPTEND+CPUTQY = Old( CPATY + CPDUP + CPDTHP )
! Calcul des tendances de T , U et de Q et modifications
! eventuelles de W et de OMEGA/P

IF (LINTFLEX.AND.(.NOT.YDCPG_OPTS%LCONFX)) THEN
  CALL CPTEND_FLEX( YDLDDH, YDMDDH, YGFL, YDPHY, YDDYNA, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,        &
  & YDCPG_OPTS%KFLEVG, YDVARS%GEOMETRY%GNORDL%T0, YDVARS%GEOMETRY%GNORDM%T0, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
  & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%U,           &
  & YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%YGSP_RR%T, PGFL, YLPROCSET,                 &
  & YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDH, ZTENDGFL, YDMF_PHYS%OUT%FHSCL, YDMF_PHYS%OUT%FHSCN,            &
  & YDMF_PHYS%OUT%FHSSL, YDMF_PHYS%OUT%FHSSN, YDMF_PHYS%OUT%FHPCL, YDMF_PHYS%OUT%FHPCN, YDMF_PHYS%OUT%FHPSL,         &
  & YDMF_PHYS%OUT%FHPSN, PFHP=ZMSC_FHP, PFP=ZPFL_FP, PFEPFP=YDMF_PHYS%OUT%FEPFP, PFCMPCQ=YDMF_PHYS%OUT%FCMPCQ,       &
  & PFCMPSN=YDMF_PHYS%OUT%FCMPSN, PFCMPSL=YDMF_PHYS%OUT%FCMPSL, YDDDH=YDDDH)
ELSE
  CALL CPTEND_NEW(YDCST,  YDMODEL, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, YDVARS%GEOMETRY%GNORDL%T0, &
  & YDVARS%GEOMETRY%GNORDM%T0, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS,                          &
  & ZDIFEXT, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCCQL,                       &
  & YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,                                  &
  & YDMF_PHYS%OUT%FPLSG, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG, YDMF_PHYS%OUT%FPEVPSL,                                &
  & YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPEVPCG,                        &
  & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,                             &
  & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG,                             &
  & YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FRMH, ZMSC_FRMQ, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%STRCU,                          &
  & YDMF_PHYS%OUT%STRCV, YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV,                                  &
  & YDMF_PHYS%OUT%STRMU, YDMF_PHYS%OUT%STRMV, YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC, YDMF_PHYS%OUT%FIMCC,                              &
  & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,                           &
  & YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC, YDMF_PHYS%OUT%FCNEGQSC, ZPFL_FTKE, ZPFL_FTKEI,                                            &
  & ZPFL_FEFB1, ZPFL_FEFB2, ZPFL_FEFB3, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                      &
  & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%U,                                         &
  & YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%Q, YDMF_PHYS_BASE_STATE%I, YDMF_PHYS_BASE_STATE%L,                   &
  & YDVARS%LCONV%T0, YDVARS%ICONV%T0, YDVARS%RCONV%T0, YDVARS%SCONV%T0, YDMF_PHYS_BASE_STATE%R, YDMF_PHYS_BASE_STATE%S,                       &
  & YDMF_PHYS_BASE_STATE%G, ZDSA_CPS, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS%OUT%FHSCL, YDMF_PHYS%OUT%FHSCN,                               &
  & YDMF_PHYS%OUT%FHSSL, YDMF_PHYS%OUT%FHSSN, YDMF_PHYS%OUT%FHSSG, YDMF_PHYS%OUT%FHPCL, YDMF_PHYS%OUT%FHPCN,                                  &
  & YDMF_PHYS%OUT%FHPCG, YDMF_PHYS%OUT%FHPSL, YDMF_PHYS%OUT%FHPSN, YDMF_PHYS%OUT%FHPSG, ZMSC_FHP,                                             &
  & ZPFL_FP, YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, YDMF_PHYS%OUT%FCMPSL, YDMF_PHYS%OUT%TENDU,                      &
  & YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDH, ZPTENDQ1, ZPTENDI1, ZPTENDL1, ZPTENDLCONV1,                                                  &
  & ZPTENDICONV1, ZPTENDRCONV1, ZPTENDSCONV1, ZPTENDR1, ZPTENDS1, ZPTENDG1, ZPTENDTKE1, ZPTENDEFB11,                                          &
  & ZPTENDEFB21, ZPTENDEFB31, ZTENDEXT, YDDDH)
ENDIF

IF (YDDYNA%LTWOTL) THEN

ELSE
    
  IF ( L3MT.OR.LSTRAPRO.OR.(NDPSFI==1)) THEN
!     PFEPFP was ZFEPFP in CPTEND_NEW, before, ZFEPFP still in CPFHPFS
    DO JLEV= 0, YDCPG_OPTS%KFLEVG 
      DO JROF = 1, YDCPG_OPTS%KLON
        YDMF_PHYS%OUT%FEPFP(JROF,JLEV) = 0.0_JPRB
        YDMF_PHYS%OUT%FCMPCQ(JROF,JLEV) = 0.0_JPRB
        YDMF_PHYS%OUT%FCMPSN(JROF,JLEV) = 0.0_JPRB
        YDMF_PHYS%OUT%FCMPSL(JROF,JLEV) = 0.0_JPRB
      ENDDO
    ENDDO
  ENDIF

ENDIF 



!        2.7.1  Diagnostics on physical tendencies
!               ----------------------------------

IF (.NOT.YDCPG_OPTS%LCONFX) THEN
  IF ((GCHETN%LFREQD).OR.(GCHETN%LCOORD).OR.(GCHETN%LPROFV)) THEN
    CALL CPCHET (YDMF_PHYS, YDMF_PHYS_BASE_STATE, YDCPG_MISC, YDRIP, YDPHY, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, &
    & YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, YDCPG_OPTS%NSTEP, ZMSC_FRMQ, ZDSA_CPS, ZTENDH, ZPTENDQ1,            &
    & ZPTENDI1, ZPTENDL1, ZPTENDR1, ZPTENDS1, YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GELAM%T0)
  ENDIF
  
  IF (GCHETN%LPROFV)&
   & CALL PROFILECHET(YDGEOMETRY, YDCPG_MISC, YDMF_PHYS, ZRDG_CVGQ, ZRDG_LCVQ, ZRDG_MU0, YDCPG_DYN0,                   &
     & YDMF_PHYS_SURF, YDVARS, YDSURF, YDDPHY, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KFDIA, YDVARS%GEOMETRY%GELAM%T0, &
     & YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GM%T0, YDVARS%GEOMETRY%OROG%T0, YDVARS%GEOMETRY%RCORI%T0,              &
     & YDVARS%GEOMETRY%RATATH%T0, YDVARS%GEOMETRY%RATATX%T0)

ENDIF

!        2.8  Modification of vertical velocities
!             by some physics output when required.
!             -------------------------------------


! * MODIFICATION DE LA VITESSE VERTICALE ET DE LA TENDANCE DE
! PRESSION DE SURFACE SI NDPSFI=1 ( MASSE VARIABLE ).
! Ajout de la physique dans l'equation de continuite/Add physics
! in continuity equation.

IF (NDPSFI == 1) THEN
  CALL CPMVVPS(YDCST, YDVAB, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, YDCPG_OPTS%ZDTPHY,  &
  & ZPFL_FP, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN,                 &
  & YDCPG_DYN0%CTY%EVEL, YDCPG_DYN0%CTY%PSDVBC, YDMF_PHYS_NEXT_STATE%SP%P)
ENDIF

!        2.9  Computation of evolution of T, u, v and Q.
!             ------------------------------------------

!  ALARO does not respect the coding rules, tendency of pseudo-TKE is computed in APLPAR and not
!  in CPTEND_NEW. To use the new version of cputqy it is then necessary to write it in GFL tendencies array.
! This memory transfer is not necessary, please respect coding rules to avoid it.

! Not necessary for intflex: already done in aplpar2intflex
IF (.NOT.(LINTFLEX.AND.(.NOT.YDCPG_OPTS%LCONFX))) THEN
  IF (LPTKE) THEN
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZPTENDTKE1(JROF,JLEV) = ZTENDPTKE(JROF,JLEV)
      ENDDO
    ENDDO    
  ENDIF
  ! Extra-GFL
  IF(LMDUST.AND.(NGFL_EXT/=0)) THEN
    DO JGFL=1, NGFL_EXT
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZTENDGFL(JROF,JLEV,YEXT(JGFL)%MP1) = ZTENDEXT(JROF,JLEV,JGFL)+&! turbulent tendency
                                             & ZTENDEXT_DEP(JROF,JLEV,JGFL) ! moist tendency
        ENDDO
      ENDDO 
    ENDDO   
  ENDIF
ENDIF

! ky: non-zero option not yet coded for the time being.
ZTENDD=0.0_JPRB

! Calcul de T , Q et du Vent a l'instant 1

IF (LHOOK) CALL DR_HOOK ('CPUTQY',0,ZHOOK_HANDLE_1)

CALL CPUTQY_APLPAR_EXPL(YDCST, YDCPG_BNDS, YDCPG_OPTS, YDDYNA, YDMF_PHYS_NEXT_STATE, YDMF_PHYS_BASE_STATE, YDVARS, &
& YDPHY, YDCPG_OPTS%ZDTPHY, ZTENDH, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDD,              &
& ZPTENDEFB11, ZPTENDEFB21, ZPTENDEFB31, ZPTENDG1, ZPTENDICONV1, ZPTENDI1, ZPTENDLCONV1, ZPTENDL1,                 &
& ZPTENDQ1, ZPTENDRCONV1, ZPTENDR1, ZPTENDSCONV1, ZPTENDS1, ZPTENDTKE1, YDMF_PHYS%OUT%FDIS)

CALL CPUTQY_APLPAR_LOOP(YDMODEL%YRML_DYN%YRDYN, YDMODEL%YRML_DYN%YRDYNA, YDGEOMETRY%YRDIMV, YDGMV, YGFL, YDPTRSLB1, YDPHY, &
& YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, &
& YDCPG_OPTS%ZDTPHY, ZTENDGFL, YDCPG_SL1%ZVIEW, PGMVT1,                &
& PGFLT1)

IF (LHOOK) CALL DR_HOOK ('CPUTQY',1,ZHOOK_HANDLE_1)

CALL MF_PHYS_FPL_PART2 (YDCPG_BNDS, YDCPG_OPTS, ZPFL_FPLCH, ZPFL_FPLSH, YDVARS%CPF%T1, YDVARS%SPF%T1, &
& YDMODEL)

!       2.9b Prognostic convection etc.
!            --------------------------

! TRANSFER NOT ADVECTED VARIABLES INTO PGFLT1
CALL MF_PHYS_TRANSFER (YDCPG_BNDS, YDCPG_OPTS, YDVARS, YDMODEL%YRML_PHY_MF%YRPHY, YDMODEL%YRML_GCONF%YGFL)

!        2.10  Surface variables.
!              ------------------

IF ((.NOT.LSFORCS)) THEN
  
  IF (.NOT.LMSE) THEN
    DO JLEV=0,YDCPG_OPTS%KFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZPFL_FPLSN(JROF,JLEV)=YDMF_PHYS%OUT%FPLSN(JROF,JLEV)+YDMF_PHYS%OUT%FPLSG(JROF,JLEV)
      ENDDO
    ENDDO
    CALL CPTENDS(YDCST,  YDMODEL%YRML_PHY_MF, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, &
    & YSP_SBD%NLEVS, YDCPG_OPTS%ZDTPHY, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLCN,                      &
    & ZPFL_FPLSN, YDMF_PHYS  %OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS_SURF%GSP_SG%PA_T1, YDMF_PHYS%OUT%CT,                  &
    & ZDSA_C1, ZDSA_C2, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL, YDMF_PHYS  %OUT%FCLN, YDMF_PHYS%OUT%FCS,                   &
    & ZFLU_FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT  %FEVN, YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS,     &
    & YDMF_PHYS%OUT%FLWSP, YDMF_PHYS%OUT%FONTE, YDMF_PHYS%OUT%FTR, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSP_SG%PR_T1, &
    & YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISS, YDMF_PHYS_SURF%GSP_SG%PF_T1,                           &
    & ZFLU_VEG, ZTDS_TDTS, ZTDS_TDTP, ZTDS_TDWS, ZTDS_TDWSI, ZTDS_TDWP, ZTDS_TDWPI, ZTDS_TDWL,                              &
    & ZTDS_TDSNS, ZTDS_TDALBNS, ZTDS_TDRHONS)  

    CALL CPWTS(YDCST, YDCPG_OPTS, YDMODEL%YRML_AOC%YRMCC, YDPHY, YDMODEL%YRML_PHY_MF%YRPHY1, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, &
    & YDCPG_BNDS%KFDIA, YSP_SBD%NLEVS, YDCPG_OPTS%ZDTPHY, ZTDS_TDTS, ZTDS_TDTP, ZTDS_TDWS, ZTDS_TDWSI,                                  &
    & ZTDS_TDWP, ZTDS_TDWPI, ZTDS_TDWL, ZTDS_TDSNS, ZTDS_TDALBNS, ZTDS_TDRHONS, YDMF_PHYS_SURF%GSD_VP%PTPC,                             &
    & YDMF_PHYS_SURF%GSD_VP%PWPC, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS_SURF%GSP_RR%PT_T1,                 &
    & YDMF_PHYS_SURF%GSP_SB%PT_T1, YDMF_PHYS_SURF%GSP_RR%PW_T1, YDMF_PHYS_SURF%GSP_RR%PIC_T1, YDMF_PHYS_SURF%GSP_SB%PQ_T1,              &
    & YDMF_PHYS_SURF%GSP_SB%PTL_T1, YDMF_PHYS_SURF%GSP_RR%PFC_T1, YDMF_PHYS_SURF%GSP_SG%PF_T1, YDMF_PHYS_SURF%GSP_SG%PA_T1,             &
    & YDMF_PHYS_SURF%GSP_SG%PR_T1                              )  
  ELSE
    IF (YDCPG_OPTS%LCONFX) THEN
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDMF_PHYS_SURF%GSP_RR%PT_T0(JROF)=YDCPG_GPAR%VTS(JROF)
      ENDDO
    ELSE
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDMF_PHYS_SURF%GSP_RR%PT_T1(JROF)=YDCPG_GPAR%VTS(JROF)
      ENDDO
    ENDIF
  ENDIF
  IF(LNUDG)THEN
    CALL CPNUDG ( YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, NFNUDG, YDCPG_OPTS%KFLEVG, YDCPG_BNDS%KBL,         &
    & XPNUDG, YDMF_PHYS_SURF%GSD_VF%PNUDM, YDMF_PHYS_SURF%GSP_RR%PT_T1, YDMF_PHYS_SURF%GSP_RR%PW_T1,                      &
    & YDMF_PHYS_SURF%GSP_SB%PQ_T1, YDMF_PHYS_SURF%GSP_SG%PF_T1, YDMF_PHYS_NEXT_STATE%T%P (:, 1:YDCPG_OPTS%KFLEVG),          &
    & YDMF_PHYS_NEXT_STATE%Q%P (:, 1:YDCPG_OPTS%KFLEVG), YDMF_PHYS_NEXT_STATE%U%P (:, 1:YDCPG_OPTS%KFLEVG),                   &
    & YDMF_PHYS_NEXT_STATE%V%P (:, 1:YDCPG_OPTS%KFLEVG), YDMF_PHYS_NEXT_STATE%SP%P, YDVARS%T%T0, YDVARS%Q%T0,                 &
    & YDVARS%U%T0, YDVARS%V%T0, YDCPG_PHY0%PREHYD(:, YDCPG_OPTS%KFLEVG), YDVARS%GEOMETRY%GM%T0, YDMF_PHYS_SURF%GSD_VF%PLSM&
    &             )
  ENDIF
ENDIF


IF(YDMODEL%YRML_PHY_MF%YRPHY%LCVPGY) THEN
  CALL MF_PHYS_CVV (YDCPG_BNDS, YDCPG_OPTS, YDVARS%CVV%T0, YDVARS%CVV%T1)
ENDIF

!        3.3  Store the model trajectory at t-dt (leap-frog) or t (sl2tl).
!             ------------------------------------------------------------

IF (LTRAJPS) THEN
  PTRAJ_PHYS%PQSSMF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_MISC%QS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
  PTRAJ_PHYS%PTSMF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA) =YDMF_PHYS_BASE_STATE%YGSP_RR%T(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
  PTRAJ_PHYS%PSNSMF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_BASE_STATE%YGSP_SG%F(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)

  IF (.NOT. YDDYNA%LTWOTL) THEN
    CALL WRPHTRAJM(YDGEOMETRY, YDSIMPHL, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, PTRAJ_PHYS, YDVARS%U%T9, &
    & YDVARS%V%T9, YDVARS%T%T9, YDVARS%Q%T9, YDVARS%L%T9, YDVARS%I%T9, YDVARS%SP%T9)  
  ENDIF

  IF (LPRTTRAJ.AND.PTRAJ_PHYS%LASTCHUNK) WRITE(NULOUT,*)'GREPTRAJ STORE TRAJ_PHYS in APLPAR'
ENDIF

!     ------------------------------------------------------------------

!*       5.    Final calculations.
!              -------------------

! * Restore the initial value of some pseudo-historical surface buffers
!   if relevant.
IF (LL_SAVE_PHSURF) THEN
  CALL MF_PHYS_SAVE_PHSURF_PART2 (YDCPG_BNDS, YDCPG_OPTS, ZSAV_DDAL, ZSAV_DDOM, ZSAV_ENTCH,                           &
  & ZSAV_FHPS, ZSAV_GZ0F, ZSAV_GZ0HF, ZSAV_HV, ZSAV_PBLH, ZSAV_QSH, ZSAV_UDAL, ZSAV_UDGRO, ZSAV_UDOM,                 &
  & ZSAV_UNEBH, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VH%PPBLH, YDMF_PHYS_SURF%GSD_VH%PQSH,                  &
  & YDMF_PHYS_SURF%GSD_VH%PSPSH, YDMF_PHYS_SURF%GSD_VK%PUDGRO, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VV%PZ0H, &
  & YDVARS%DAL%T0, YDVARS%DOM%T0, YDVARS%UAL%T0, YDVARS%UEN%T0, YDVARS%UNEBH%T0, YDVARS%UOM%T0,                       &
  & YDMODEL)
ENDIF

! Store horizontal exchange coefficients (3D turbulence) to SL2 buffers
IF (YDDYNA%L3DTURB) THEN
  DO JLEV=1,YDCPG_OPTS%KFLEVG
    YDCPG_SL2%KAPPAM (YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, JLEV) = ZKUR_KUROV_H(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
    YDCPG_SL2%KAPPAH (YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, JLEV) = ZKUR_KTROV_H(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
  ENDDO
ENDIF

CALL MF_PHYS_BAYRAD (YDCPG_BNDS, YDCPG_OPTS, ZBAY_QRCONV, ZBAY_QSCONV, YDVARS%RCONV%T1, YDVARS%SCONV%T1, &
& YDMODEL)


!-------------------------------------------------
! Extract Single Column Model profiles from 3D run or 
! write LFA file for MUSC (1D model)
!-------------------------------------------------
IF(LGSCM.OR.LMUSCLFA) THEN
  IF (LAROME) THEN
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDCPG_MISC%NEB(JROF,JLEV)=YDVARS%A%T1(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
  CALL WRITEPHYSIO(YDGEOMETRY, YDCPG_MISC, YDCPG_PHY0, YDMF_PHYS, YDCPG_DYN0, YDMF_PHYS_SURF, YDVARS,         &
  & YDSURF, YDDPHY, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KFDIA, YDCPG_BNDS%KIDIA, YDCPG_OPTS%KGL1,          &
  & YDCPG_OPTS%KGL2, YDCPG_BNDS%KSTGLO, YDCPG_OPTS%NSTEP, NTSSG, YSP_SBD%NLEVS, YDVARS%GEOMETRY%GELAM%T0,     &
  & YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GM%T0, YDVARS%GEOMETRY%OROG%T0, YDVARS%GEOMETRY%RCORI%T0,        &
  & YDVARS%GEOMETRY%RATATH%T0, YDVARS%GEOMETRY%RATATX%T0, YDVARS%GEOMETRY%GECLO%T0, YDVARS%GEOMETRY%GESLO%T0, &
  & ZRDG_CVGQ, ZRDG_LCVQ, ZRDG_MU0, ZDSA_C1, ZDSA_C2, ZDSA_CPS, ZDSA_LHS, ZDSA_RS, ZFLU_CD, ZFLU_CDN,         &
  & ZFLU_CH, ZFLU_EMIS, ZFLU_FEVI, ZFLU_NEIJ, ZFLU_QSAT, ZFLU_QSATS, ZFLU_VEG, IMOC_CLPH, ZMSC_FRMQ,          &
  & ZMSC_LH, ZMSC_LSCPE, ZMSC_QW, ZMSC_TW, ZPFL_FEFB1, ZPFL_FEFB2, ZPFL_FEFB3, ZPFL_FPLCH,                    &
  & ZPFL_FPLSH, ZPFL_FTKE  )
ENDIF

IF (LEDR) THEN
  YDMF_PHYS_SURF%GSD_DI%PXEDR(:,:)=YDMF_PHYS%OUT%EDR(:,:)
ENDIF

CALL MF_PHYS_PRECIPS (YDCPG_BNDS, YDCPG_OPTS, ZPRC_DPRECIPS, ZPRC_DPRECIPS2, YDMF_PHYS_SURF%GSD_XP%PPRECIP, &
& YDMF_PHYS_SURF%GSD_XP2%PPRECIP2, YDMODEL)

!     ------------------------------------------------------------------

!       6. destructor for procset
IF (LINTFLEX) CALL CLEANINTPROCSET(YLPROCSET)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('APLPAR',1,ZHOOK_HANDLE)
END SUBROUTINE APLPAR
