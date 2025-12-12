#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE APL_AROME(YDCST, YDMF_PHYS_BASE_STATE, YDMF_PHYS_NEXT_STATE, YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, &
& YDFORCESPPT, YDCPG_MISC, YDGPAR, YDCPG_PHY0, YDMF_PHYS, YDRADF, YDHGRAD, YDCPG_DYN0, YDMF_PHYS_SURF, &
& YDVARS, YDGEOMVARS, YDSURF, YDCFU, YDXFU, YDMODEL, PGFL, PGP2DSDT, YDDDH)

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
!
!   YDCST:                Model constants
!   YDMF_PHYS_BASE_STATE:
!   YDMF_PHYS_NEXT_STATE:
!   YDGEOMETRY:           Geomtery settings from IFS model
!   YDCPG_BNDS:           Array dimensions and bounds
!   YDCPG_OPTS:           Constant input variables
!   YDCPG_MISC:
!   YDGPAR:
!   YDCPG_PHY0:
!   YDMF_PHYS:            Physics fields
!   YDHGRAD:
!   YDCPG_DYN0:           Dynamic fields
!   YDMF_PHYS_SURF:       Surface fields
!   YDVARS:               Persistent fields
!   YDSURF:               Various surface fields and, declarations, and definitions
!   YDCFU:                Controls for cumulated fluxes
!   YDXFU:                Controls for instantaneous fluxes
!   YDMODEL:              Model configuration settings
!   PGFL:                 GFL fields (Should be removed)
!   PGP2DSDT:             SPPT perturbation pattern
!   YDDDH:                DDH fields and settings
!
!-----------------------------------------------------------------------

!     Externes.
!     ---------

!     Method
!     ------
!     - convert aladin variables into mesonh variables (level inversion
!       and q to r, t to theta)
!     - call mesoNH physics and ECMWF radiation scheme
!     - convert mesoNH tendencies to aladin tendencies
!
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
!                            (YDCPG_OPTS%KLON,NGPAR) is used boundary checking bf
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
!     K.I Ivarsson 2018-02 : Some new variables for microphysics
!     2018-09, E. Gleeson: Corrected misplaced arguments in ACRANEB2 call.
!     2019-09-24 J.M. Piriou arguments for convective gusts.
!     R. El Khatib 30-Oct-2018 substantial rewrite for optimization and coding standards respect.
!     2018-10, I. Etchevers : add Visibilities
!     2019-01, I. Etchevers, Y. Seity : add Precipitation Type
!     2019-06, W. de Rooy: Modifications for new set-up statistical cloud scheme (LSTATNW)
!     2019-09, J. Masek: Corrected dimensioning of dummy argument PGMU0.
!                        Modified call to ACRANEB2 (clearsky fluxes).
!     2019-10, I. Etchevers : Visibilities in ACVISIH, AROCLDIA=>ACCLDIA
!     2019-10, Y.Bouteloup and M. Bouzghaiam : Radiation modifications. Remove acradin.F90 direct
!              call to recmwf.F90 and add interface to ecrad (in recmwf !)
!     2020-10, J. Masek: Modified call to ACCLDIA.
!     2020-12, F. Meier add call to latent heat nudging if LNUDGLH is TRUE
!     2020-12, U. Andrae : Introduce SPP for HARMONIE-AROME
!     2021-11, R. Honnert: Add Thermal vertical velocity
!     2021-12, F. Meier add call to windfarm parametrisation if LWINDFARM is TRUE
!     R. El Khatib 08-Jul-2022 Contribution to the encapsulation of YOMCST and YOETHF
!     2022-04-13 J.M. Piriou : Sun eclipses (solar obscuration ZRDG_SOLO).
!     2023-09 R. H. Myhre : Major cleanup and split routine into multiple subroutines
!     2024-11 R. H. Myhre : New refactoring cycle. All array dimensions lead with KLON and SWAP removed
!     A. Marcel Jan 2025: EDMF contribution to dynamic TKE production
!     2025-02 M. Mokhtari: Fixed confusion between LIMA and DUST in TURB
! End modifications
!-------------------------------------------------------------------------------

USE PARKIND1,                      ONLY: JPIM, JPRB, JPRD
USE YOMHOOK,                       ONLY: LHOOK, DR_HOOK, JPHOOK

USE GEOMETRY_MOD,                  ONLY: GEOMETRY
USE MF_PHYS_TYPE_MOD,              ONLY: MF_PHYS_TYPE
USE CPG_TYPE_MOD,                  ONLY: CPG_MISC_TYPE, CPG_DYN_TYPE, CPG_PHY_TYPE
USE YOMGPAR,                       ONLY: TGPAR
USE CPG_OPTS_TYPE_MOD,             ONLY: CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,      ONLY: MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD,           ONLY: FIELD_VARIABLES
USE YOMGEOMVARS,                   ONLY: TGEOMVARS         
USE SURFACE_FIELDS_MIX,            ONLY: TSURF
USE YOMXFU,                        ONLY: TXFU
USE YOMCFU,                        ONLY: TCFU
USE TYPE_MODEL,                    ONLY: MODEL
USE YOMCST,                        ONLY: TCST
USE DDH_MIX,                       ONLY: TYP_DDH ! for new data flow
USE SPP_MOD_TYPE,                  ONLY: UA_PHYS_SPP_VARS
USE MF_PHYS_BASE_STATE_TYPE_MOD,   ONLY: MF_PHYS_BASE_STATE_TYPE
USE MF_PHYS_NEXT_STATE_TYPE_MOD,   ONLY: MF_PHYS_NEXT_STATE_TYPE
USE SC2PRG_MOD,                    ONLY: SC2PRG
USE YOMDGRADIENT_TYPE_MOD,         ONLY: GRADIENT_TYPE
USE YOMFORCESPPT,                  ONLY: TFORCESPPT
USE YOMRADF,                       ONLY: TRADF

!     -------------------------------------------------------------------------

IMPLICIT NONE


TYPE(TCST)                    ,INTENT(IN)     :: YDCST
TYPE(MF_PHYS_BASE_STATE_TYPE) ,INTENT(IN)     :: YDMF_PHYS_BASE_STATE
TYPE(MF_PHYS_NEXT_STATE_TYPE) ,INTENT(INOUT)  :: YDMF_PHYS_NEXT_STATE
TYPE(GEOMETRY)                ,INTENT(IN)     :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE)           ,INTENT(IN)     :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE)           ,INTENT(IN)     :: YDCPG_OPTS
TYPE(TFORCESPPT)              ,INTENT(IN)     :: YDFORCESPPT
TYPE(CPG_MISC_TYPE)           ,INTENT(INOUT)  :: YDCPG_MISC
TYPE(TGPAR)                   ,INTENT(INOUT)  :: YDGPAR
TYPE(CPG_PHY_TYPE)            ,INTENT(IN)     :: YDCPG_PHY0
TYPE(MF_PHYS_TYPE)            ,INTENT(INOUT)  :: YDMF_PHYS
TYPE(TRADF)                   ,INTENT(INOUT)  :: YDRADF
TYPE(GRADIENT_TYPE)           ,INTENT(IN)     :: YDHGRAD
TYPE(CPG_DYN_TYPE)            ,INTENT(IN)     :: YDCPG_DYN0
TYPE(MF_PHYS_SURF_TYPE)       ,INTENT(INOUT)  :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES)         ,INTENT(INOUT)  :: YDVARS
TYPE(TGEOMVARS)               ,INTENT(IN)     :: YDGEOMVARS
TYPE(TSURF)                   ,INTENT(IN)     :: YDSURF
TYPE(TCFU)                    ,INTENT(IN)     :: YDCFU
TYPE(TXFU)                    ,INTENT(IN)     :: YDXFU
TYPE(MODEL)                   ,INTENT(IN)     :: YDMODEL
REAL(KIND=JPRB)               ,INTENT(INOUT)  :: PGFL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)               ,INTENT(IN)     :: PGP2DSDT(YDCPG_OPTS%KLON,YDMODEL%YRML_SPPT%YGPSDT(1)%NG2D,YDMODEL%YRML_SPPT%N2D)
TYPE(TYP_DDH)                 ,INTENT(INOUT)  :: YDDDH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                                READ ME, PLEASE !
!
!                         CODING CONVENTIONS FOR ARPEGE vs MNH PHYSICS
!
! Local 3D arrays with extra levels for Meso-NH turbulence scheme :
! - suffixed with two underscore to be easily identified
!
! Local 3D arrays with regular number of levels for Meso-NH interfaces :
! - suffixed with one underscore to be easily identified.
!
! Local 4D arrays with regular number of levels for Meso-NH interfaces :
! - suffixed with one underscore to be easily identified
!
! Local 2D arrays for Meso-NH interfaces :
! - suffixed with one underscore to be easily identified
!
! Other arrays, which can be dummy arguments, or local but used as argument to IFS/ARPEGE physics
! must be dimensioned KLON and should not be suffixed with undersores.
!
!                         DO NOT USE ARRAY SYNTAX FOR COMPUTATIONAL LOOPS !!
!
! - Explicit loops are needed for code transformation
! - They make the code less performant because memory cache is poorly used
! - They can make the code even less readable if the indexes are removed
!
!                         AVOID ARRAYS COPIES, OR MAKE THEM FAST !
!
! - if you do need to initialize or copy an array, do it as follows with explicit array syntax in first dimension
! because the compiler will be able to use an optimized function to initialize/copy a segment of memory,
! and may be able to address simultaneously several cach lines :
! 1D array :
!   Z(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=value
! 2D arrays :
!  DO JLEV=1,YDCPG_OPTS%KFLEVG
!     ZX(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=xval
!     ZY(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=yval
!  ENDDO
!
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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(KIND=JPRB) :: ZRHODJM__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZRHODREFM__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZPABSM__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)

REAL(KIND=JPRB) :: ZUM__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZVM__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZWM__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZUS__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZVS__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZWS__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)

REAL(KIND=JPRB) :: ZTHM__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)

REAL(KIND=JPRB) :: ZTKES_OUT__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZMF_UP__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZTHVREFM__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)  ! thetav de l etat

REAL(KIND=JPRB) :: ZTENDU_TURB__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZTENDV_TURB__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)

REAL(KIND=JPRB) :: ZTKEM__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZSRCS__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZSIGS__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)

! length scales for momentum and heat for mnh level definitions in case LHARATU=TRUE
REAL(KIND=JPRB) :: ZLENGTHM__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZLENGTHH__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)

REAL(KIND=JPRB) :: ZFLXZTHVMF_SUM__(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZFLXZUMF_SUM__(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZFLXZVMF_SUM__(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG+1)

! Subgrid autoconversions
REAL(KIND=JPRB) :: ZHLC_HRC_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZHLC_HCF_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZHLI_HRI_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZHLI_HCF_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

! Updraft characteristics for Meso-NH world (input of ARO_SHALLOW_MF)
REAL(KIND=JPRB) :: ZTKES_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZZZ_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZDZZ_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZZZ_F_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZEXNREFM_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZNEBMNH_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

! additions for future ice cloud fraction and precipitation fraction
REAL(KIND=JPRB) :: ZICEFR_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZPRCFR_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

! additions for MF scheme (Pergaud et al)
REAL(KIND=JPRB) :: ZSIGMF_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZAERD_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZDTHRAD_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZTHS__(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZRS_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
REAL(KIND=JPRB) :: ZSVM_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB) :: ZLIMAS_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NLIMA)
REAL(KIND=JPRB) :: ZSVS_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

REAL(KIND=JPRB) :: ZLIMAM_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NLIMA)

! NRT aerosol
! Aerosol species (MR:Mass mixing ratio, NC:number condentration)
REAL(KIND=JPRB) :: ZAEROM_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NAERO)
! Cloud condensation nuclei, ice nuclei and cloud droplets
REAL(KIND=JPRB) :: ZIFN_NC_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZCLDROP_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZRM_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

! Single scattering albedo of dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZPIZA_DST_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)
! Assymetry factor for dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZCGA_DST_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)
! tau/tau_{550} dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZTAUREL_DST_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)

! surface flux of theta and surface flux of vapor ; surface flux of CO2
REAL(KIND=JPRB) :: ZSFTH_(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZSFRV_(YDCPG_OPTS%KLON)

REAL(KIND=JPRB) :: ZACPRG_(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZINPRG_NOTINCR_(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZINPRG_(YDCPG_OPTS%KLON)

REAL(KIND=JPRB) :: ZACPRR_(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZINPRR_NOTINCR_(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZINPRR_(YDCPG_OPTS%KLON)

REAL(KIND=JPRB) :: ZACPRS_(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZINPRS_NOTINCR_(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZINPRS_(YDCPG_OPTS%KLON)

REAL(KIND=JPRB) :: ZINPRH_(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZTOWNS_(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZBUDTH_(YDCPG_OPTS%KLON)

REAL(KIND=JPRB) :: ZSFSV_(YDCPG_OPTS%KLON, YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! surf. flux of scalars

! surface flux of x and y component of wind. are they really necessary ? REK
REAL(KIND=JPRB) :: ZSFU_(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZSFV_(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZZS_(YDCPG_OPTS%KLON)

! local areas for OCND2 option:
! ZICLDFR = ice cloud fraction , ZWCLDFR = water or mixed-phase cloud fraction,
! ZSSIO = Super-saturation with respect to ice in ZICLDFR ,
! ZSSIU = Sub-saturation with respect to ice outside ZICLDFR,
! ZIFR = variable used for calulation of subgridscale ice
! Meso-NH world
REAL(KIND=JPRB) :: ZICLDFR_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZWCLDFR_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZSSIO_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZSSIU_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZIFR_(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

! Arpege-style dimensionning :
! --------------------------

!Variables used in case LHARATU=TRUE
! length scales for momentum and heat and TKE
REAL(KIND=JPRB) :: ZTKEEDMFS(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZEMIS(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZQSAT(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQDM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTKEM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDTT(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG) ! array to save heating profile for LHN

! POUR GROUND
REAL(KIND=JPRB) :: ZZS_FSWDIR(YDCPG_OPTS%KLON, YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZZS_FSWDIF(YDCPG_OPTS%KLON, YDMODEL%YRML_PHY_RAD%YRERAD%NSW)

REAL(KIND=JPRB) :: ZALBD(YDCPG_OPTS%KLON, YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZALBP(YDCPG_OPTS%KLON, YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZAPHIM(YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZAPHIFM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZTM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQVM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQIM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQCM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQHM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQRM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQSM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQGM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZCPM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZRHM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

! Ajout pour MF Dual Scheme (KNMI et al)
! Updraft characteristics in Arpege/IFS world
REAL(KIND=JPRB) :: ZTSURF(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZQS(YDCPG_OPTS%KLON)

REAL(KIND=JPRB) :: ZSURFPREP(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZSURFSNOW(YDCPG_OPTS%KLON)

!    Integers
INTEGER(KIND=JPIM) :: JLEV, JLON
INTEGER(KIND=JPIM) :: JSPP, JSW

INTEGER(KIND=JPIM), PARAMETER :: INIT0 = -1 ! Kind of safety/debugging initialization :
                                            !  0 = initialize to HUGE (debugging)
                                            !  1 = initialize to realistic value (discouraged)
                                            ! -1 = no initialization (optimized code) - this is the default.

REAL(KIND=JPRB) :: ZDT, ZINVG
 ! pas de temps pour la surface externalise
REAL(KIND=JPRB) :: ZADTMS

! default values for initialization :
REAL(KIND=JPRB) :: ZVALUE, ZVALUE_ONE, ZVALUE_T, ZVALUE_P, ZVALUE_L, ZVALUE_EPSILON

!       Boolean
LOGICAL :: LLMSE

TYPE(UA_PHYS_SPP_VARS) :: ZSPP_UA

!     --- UPPER AIR PHYSICAL TENDENCIES FOR AROME.
!       (the previous one are not used in AROME)
REAL(KIND=JPRB) :: ZTENDT (YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG) ! temperature tendency
REAL(KIND=JPRB) :: ZTENDW (YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG) ! W  tendency

!     ---FOR AROME PHYSICS  ---

REAL (KIND=JPRB) :: ZSAV_GZ0F(YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_UDOM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZRDG_LCVQ(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZSAV_DDAL(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZRDG_MU0M(YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZRDG_CVGQ(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZRDG_CVGT(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZRDG_MU0N(YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_UDAL(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZSAV_UDGRO(YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_HV(YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_QSH(YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_DDOM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZSAV_UNEBH(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZSAV_MO(YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZRDG_MU0(YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZRDG_SOLO(YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_ENTCH(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZRDG_MU0LU(YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_PBLH(YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_GZ0HF(YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_FHPS(YDCPG_OPTS%KLON)

! SPP
REAL(KIND=JPRB) :: ZGP2DSPP(YDCPG_OPTS%KLON,YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL)

REAL(KIND=JPRB) :: ZTENDRA(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

REAL(KIND=JPRB) :: ZTENDLIMA(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NLIMA)
REAL(KIND=JPRB) :: ZTENDAERO(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NAERO) ! 3D N.R.T. AEROSOL MMR TENDENCIES
REAL(KIND=JPRB) :: ZTENDTKE(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDEXT(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

REAL(KIND=JPRB), POINTER :: ZP1EZDIAG(:,:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

#include "cpphinp.intfb.h"
#include "writephysio.intfb.h"
#include "mf_phys_save_phsurf_part1.intfb.h"
#include "mf_phys_save_phsurf_part2.intfb.h"
#include "mf_phys_transfer.intfb.h"
#include "transfer_cloud_fields.intfb.h"
#include "ecr1d.intfb.h"
#include "aro_startbu.h"
#include "checkmv.intfb.h"
#include "checknanaro.intfb.h"

#include "apl_arome_init.intfb.h"
#include "apl_arome_adjust.intfb.h"
#include "apl_arome_dust.intfb.h"
#include "apl_arome_surface_forcing.intfb.h"
#include "apl_arome_radiation.intfb.h"
#include "apl_arome_convection.intfb.h"
#include "apl_arome_surface.intfb.h"
#include "apl_arome_shallow.intfb.h"
#include "apl_arome_turbulence.intfb.h"
#include "apl_arome_windfarm.intfb.h"
#include "apl_arome_micro.intfb.h"
#include "apl_arome_nudge.intfb.h"
#include "apl_arome_sppt.intfb.h"
#include "apl_arome_final.intfb.h"
#include "apl_arome_ddh.intfb.h"
#include "apl_arome_tw_tendency.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('APL_AROME',0,ZHOOK_HANDLE)

ASSOCIATE(KLON => YDCPG_OPTS%KLON, KIDIA => YDCPG_BNDS%KIDIA, KFDIA => YDCPG_BNDS%KFDIA, KFLEVG => YDCPG_OPTS%KFLEVG)

CALL SC2PRG(1, YDMODEL%YRML_GCONF%YGFL%YEZDIAG(:)%MP, YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG, PGFL, ZP1EZDIAG)

!     ------------------------------------------------------------------
!        1.    Preliminary calculations necessary
!              for all types of physics.
!              ------------------------------------


! SPP
IF ( YDMODEL%YRML_GCONF%YRSPP_CONFIG%LSPP ) THEN
  DO JSPP = 1, YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL
    DO JLON = KIDIA, KFDIA
      ZGP2DSPP(JLON, JSPP) = YDMODEL%YRML_SPP%GP_ARP(JSPP)%GP2D(JLON, 1, YDCPG_BNDS%KBL)
    ENDDO
  ENDDO
ENDIF

IF (.NOT. YDMODEL%YRML_PHY_EC%YREPHY%LAGPHY) THEN

  CALL CPPHINP (YDGEOMETRY, YDMODEL, KIDIA, KFDIA, YDGEOMVARS%GEMU, YDGEOMVARS%GELAM, YDVARS%U%T0, YDVARS%V%T0, &
    & YDVARS%T%DL, YDVARS%T%DM, YDVARS%Q%T0, YDVARS%Q%DL, YDVARS%Q%DM, YDVARS%CVGQ%DL, YDVARS%CVGQ%DM, YDCPG_PHY0%XYB%RDELP, &
    & YDCPG_DYN0%CTY%EVEL, YDVARS%CVGQ%T0, ZRDG_MU0, ZRDG_SOLO, ZRDG_MU0LU, ZRDG_MU0M, ZRDG_MU0N, ZRDG_CVGQ, ZRDG_CVGT)

  ZRDG_LCVQ(:,:) = ZRDG_CVGQ(:,:)

ENDIF

! * In some cases, some pseudo-historic surface buffers (like z0) should
!   not be modified between the entrance and the output of APL_AROME
!   (this is the case for example if LDCONFX=T).
!   For the time being, we must save:
!   - HV (group VV) : resistance to evapotranspiration
!   - Z0F (group VD): gravity * surface roughness length
!   - Z0H (group VV): gravity * roughness length for heat
!   - PBLH (group VH): PBL height
!   - SPSH (group VH):
!   - QSH (group VH):


IF (.NOT. YDMODEL%YRML_PHY_EC%YREPHY%LAGPHY .AND. YDCPG_OPTS%LCONFX) THEN
  CALL MF_PHYS_SAVE_PHSURF_PART1 (YDCPG_BNDS, YDCPG_OPTS, ZSAV_DDAL, ZSAV_DDOM, ZSAV_ENTCH,                              &
  & ZSAV_FHPS, ZSAV_GZ0F, ZSAV_GZ0HF, ZSAV_HV, ZSAV_PBLH, ZSAV_QSH, ZSAV_UDAL, ZSAV_UDGRO,                               &
  & ZSAV_UDOM, ZSAV_UNEBH, ZSAV_MO, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VH%PPBLH, YDMF_PHYS_SURF%GSD_VH%PQSH, &
  & YDMF_PHYS_SURF%GSD_VH%PSPSH, YDMF_PHYS_SURF%GSD_VK%PUDGRO, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VV%PZ0H,    &
  & YDVARS%DAL%T0, YDVARS%DOM%T0, YDVARS%UAL%T0, YDVARS%UEN%T0, YDVARS%UNEBH%T0, YDVARS%UOM%T0, YDVARS%MO%T0,            &
  & YDMODEL)
ENDIF

CALL MF_PHYS_TRANSFER(YDCPG_BNDS, YDCPG_OPTS, YDVARS, YDMODEL%YRML_PHY_MF%YRPHY, YDMODEL%YRML_GCONF%YGFL)

!    ---------------------------------------------------------------------
!     0 - Check magnitude of model variables.
!    ---------------------------------------------------------------------

IF (YDMODEL%YRML_PHY_MF%YRPHY%LGCHECKMV) THEN

  CALL CHECKMV (YDCPG_OPTS%NINDAT, YDCST, YDMODEL%YRML_GCONF%YRRIP, YDMODEL%YRML_PHY_MF%YRPHY0, &
    & YDMODEL%YRML_PHY_MF%YRPHY2, KIDIA, KFDIA, KLON, KFLEVG, YDCPG_OPTS%NSTEP, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, &
    & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF, &
    & YDGEOMVARS%GELAM, YDGEOMVARS%GEMU, ZRDG_MU0, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_BASE_STATE%T%P, &
    & YDMF_PHYS_BASE_STATE%Q%P, YDGPAR%VTS)

ENDIF

!    ---------------------------------------------------------------------
!     1 - Initialisations
!    ---------------------------------------------------------------------

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

!         1.3 time step initialisation
!             the mesoNH physics (turb and microphysics) is written
!             for leap frog scheme
!             !!! be carefull for 2TL or 3TL

IF (YDMODEL%YRML_DYN%YRDYNA%LTWOTL) THEN
  ZDT=YDCPG_OPTS%ZDTPHY/2._JPRB
ELSE
  IF (YDCPG_OPTS%NSTEP/=0) THEN
    ZDT=YDCPG_OPTS%ZDTPHY/2._JPRB
  ELSE
    ZDT=YDCPG_OPTS%ZDTPHY
  ENDIF
ENDIF

ZINVG=1._JPRB/YDCST%RG

IF (YDCPG_OPTS%LCONFX) THEN
  ZADTMS = 0._JPRB
ELSE
  ZADTMS = YDCPG_OPTS%ZDTPHY
ENDIF


LLMSE = YDMODEL%YRML_PHY_MF%YRARPHY%LMSE .AND. (YDMODEL%YRML_PHY_MF%YRMSE%NSURFEXCTL >= 2)

!  SETUP

IF (INIT0 >= 0) THEN

  ZRHODJM__(:,:)        = ZVALUE
  ZRHODREFM__(:,:)      = ZVALUE_ONE
  ZPABSM__(:,:)         = ZVALUE_P

  ZUM__(:,:)            = ZVALUE
  ZVM__(:,:)            = ZVALUE
  ZWM__(:,:)            = ZVALUE
  ZUS__(:,:)            = ZVALUE
  ZVS__(:,:)            = ZVALUE
  ZWS__(:,:)            = ZVALUE

  ZTHM__(:,:)           = ZVALUE_T

  ZTKES_OUT__(:,:)      = ZVALUE
  ZMF_UP__(:,:)         = ZVALUE
  ZTHVREFM__(:,:)       = ZVALUE

  ZTENDU_TURB__(:,:)    = ZVALUE
  ZTENDV_TURB__(:,:)    = ZVALUE

  ZTKEM__(:,:)          = ZVALUE
  ZSRCS__(:,:)          = ZVALUE
  ZSIGS__(:,:)          = ZVALUE

  ZLENGTHM__(:,:)       = ZVALUE_L
  ZLENGTHH__(:,:)       = ZVALUE_L

  ZFLXZTHVMF_SUM__(:,:) = ZVALUE
  ZFLXZUMF_SUM__(:,:) = ZVALUE
  ZFLXZVMF_SUM__(:,:) = ZVALUE

  ZHLC_HRC_(:,:)        = ZVALUE
  ZHLC_HCF_(:,:)        = ZVALUE
  ZHLI_HRI_(:,:)        = ZVALUE
  ZHLI_HCF_(:,:)        = ZVALUE

  ZTKES_(:,:)           = ZVALUE
  ZZZ_(:,:)             = ZVALUE
  ZDZZ_(:,:)            = ZVALUE
  ZZZ_F_(:,:)           = ZVALUE

  ZEXNREFM_(:,:)        = ZVALUE
  ZNEBMNH_(:,:)         = ZVALUE

  ZICEFR_(:,:)          = ZVALUE
  ZPRCFR_(:,:)          = ZVALUE

  ZSIGMF_(:,:)          = ZVALUE
  ZAERD_(:,:)           = ZVALUE

  ZDTHRAD_(:,:)         = ZVALUE

  ZTHS__(:,:)           = ZVALUE
  ZRS_(:,:,:)           = ZVALUE
  ZSVM_(:,:,:)          = ZVALUE
  ZLIMAS_(:,:,:)        = ZVALUE
  ZSVS_(:,:,:)          = ZVALUE

  ZLIMAM_(:,:,:)        = ZVALUE

  ZRM_(:,:,:)           = ZVALUE

  ZPIZA_DST_(:,:,:)     = ZVALUE
  ZCGA_DST_(:,:,:)      = ZVALUE
  ZTAUREL_DST_(:,:,:)   = ZVALUE_EPSILON

  ZSFTH_(:)             = ZVALUE
  ZSFRV_(:)             = ZVALUE

  ZACPRG_(:)            = ZVALUE
  ZINPRG_NOTINCR_(:)    = ZVALUE
  ZINPRG_(:)            = ZVALUE

  ZACPRR_(:)            = ZVALUE
  ZINPRR_NOTINCR_(:)    = ZVALUE
  ZINPRR_(:)            = ZVALUE

  ZACPRS_(:)            = ZVALUE
  ZINPRS_NOTINCR_(:)    = ZVALUE
  ZINPRS_(:)            = ZVALUE

  ZINPRH_(:)            = ZVALUE
  ZTOWNS_(:)            = ZVALUE
  ZBUDTH_(:)            = ZVALUE

  ZSFU_(:)              = ZVALUE
  ZSFV_(:)              = ZVALUE

  ZICLDFR_(:,:)         = ZVALUE
  ZWCLDFR_(:,:)         = ZVALUE
  ZSSIO_(:,:)           = ZVALUE
  ZSSIU_(:,:)           = ZVALUE
  ZIFR_(:,:)            = ZVALUE

  ZTKEEDMFS(:,:)        = ZVALUE
  ZEMIS(:)              = ZVALUE_ONE
  ZQSAT(:,:)            = ZVALUE
  ZQDM(:,:)             = ZVALUE
  ZTKEM(:,:)            = ZVALUE
  ZTENDTT(:,:)          = ZVALUE

  ZZS_FSWDIR(:,:)       = ZVALUE
  ZZS_FSWDIF(:,:)       = ZVALUE

  ZALBD(:,:)            = ZVALUE
  ZALBP(:,:)            = ZVALUE
  ZAPHIM(:,:)           = ZVALUE

  ZTM(:,:)              = ZVALUE
  ZQVM(:,:)             = ZVALUE
  ZQIM(:,:)             = ZVALUE
  ZQCM(:,:)             = ZVALUE
  ZQHM(:,:)             = ZVALUE
  ZQRM(:,:)             = ZVALUE
  ZQSM(:,:)             = ZVALUE
  ZQGM(:,:)             = ZVALUE

  ZCPM(:,:)             = ZVALUE
  ZRHM(:,:)             = ZVALUE

  ZTSURF(:)             = ZVALUE
  ZQS(:)                = ZVALUE

  ZSURFPREP(:)          = ZVALUE
  ZSURFSNOW(:)          = ZVALUE

  YDMF_PHYS%FRSOC(:,:) = ZVALUE
  YDMF_PHYS%FRTHC(:,:) = ZVALUE
  YDMF_PHYS%FRSOPS(:)  = ZVALUE
  YDMF_PHYS%DIAGH(:)   = ZVALUE

ENDIF

!  INITIALIZE (CUMULATED) TENDENCIES

ZTENDT(:,:) = 0.0_JPRB
YDMF_PHYS%TENDU(:,:) = 0.0_JPRB
YDMF_PHYS%TENDV(:,:) = 0.0_JPRB
ZTENDW(:,:) = 0.0_JPRB

ZTENDRA(:,:,:)   = 0.0_JPRB
ZTENDLIMA(:,:,:) = 0.0_JPRB
ZTENDTKE(:,:)    = 0.0_JPRB
ZTENDEXT(:,:,:)  = 0.0_JPRB

!  INITIALIZE CUMULATED STUFF

! Small array, OK. REK
ZINPRH_(:) = 0._JPRB
ZINPRR_(:) = 0._JPRB
ZACPRR_(:) = 0._JPRB
ZINPRS_(:) = 0._JPRB
ZACPRS_(:) = 0._JPRB
ZINPRG_(:) = 0._JPRB


DO JLEV = 1,KFLEVG
  DO JLON = KIDIA,KFDIA
    ZZZ_F_(JLON, JLEV) = YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(JLON, JLEV)*ZINVG
  ENDDO
ENDDO
ZTENDTT(:,:) = 0._JPRB

! 1.5 SPP settings
IF (YDMODEL%YRML_GCONF%YRSPP_CONFIG%LSPP) THEN
  CALL ZSPP_UA%SET(KLON,KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG, &
   &  KIDIA,KFDIA,YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL, &
   &  ZGP2DSPP,ZP1EZDIAG, &
   &  YDMODEL%YRML_GCONF%YRSPP_CONFIG)
ENDIF

!    ---------------------------------------------------------------------
!     2 - CHANGEMENTS DE VARIABLES ET INVERSION DES NIVEAUX
!    ---------------------------------------------------------------------

IF (YDMODEL%YRML_PHY_MF%YRARPHY%LMICRO .OR. YDMODEL%YRML_PHY_MF%YRARPHY%LTURB .OR. &
    & LLMSE .OR. YDMODEL%YRML_PHY_MF%YRARPHY%LKFBCONV) THEN

  CALL APL_AROME_INIT(YDCST, YDMODEL, YDCPG_BNDS, YDCPG_OPTS, YDGPAR, &
                    & YDMF_PHYS_BASE_STATE, &
                    & LLMSE, YDCPG_OPTS%LSFORCS, &
                    & ZQDM, ZRHODREFM__, ZRHODJM__, ZEXNREFM_, ZPABSM__, &
                    & ZUM__, ZVM__, ZWM__, ZUS__, ZVS__, ZWS__, &
                    & ZTKEM, ZTKEM__, ZTKES_, &
                    & ZZZ_, ZRM_, ZTHM__, ZTHS__, &
                    & ZTHVREFM__, ZRS_, &
                    & ZACPRR_, ZACPRS_, ZACPRG_, &
                    & ZINPRR_NOTINCR_, ZINPRS_NOTINCR_, ZINPRG_NOTINCR_, &
                    & ZSVM_, ZSVS_, ZLIMAM_, ZLIMAS_, ZAEROM_)

ENDIF

!    ---------------------------------------------------------------------
!     4 - ADJUSTMENT (CALLED IF THE MICROPHYSICS IS SWITCH ON)
!    ---------------------------------------------------------------------


CALL APL_AROME_ADJUST(YDCST, YDMODEL, YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, YDGEOMETRY, YDVARS, &
                    & YDGEOMVARS, YDMF_PHYS, YDMF_PHYS_BASE_STATE, &
                    & ZSPP_UA, YDDDH, &
                    & ZRHODJM__, ZEXNREFM_, &
                    & ZPABSM__, ZWM__, ZTKEM__, &
                    & ZRHODREFM__, &
                    & ZTHM__, ZRM_, ZLIMAM_, ZTENDRA, &
                    & ZTHS__(:, 1:KFLEVG), ZRS_, ZLIMAS_, &
                    & ZTENDLIMA, ZZZ_F_, ZTENDT, &
                    & ZSRCS__, ZNEBMNH_, &
                    & ZICLDFR_, ZWCLDFR_, ZSSIO_, ZSSIU_, ZIFR_, &
                    & ZQDM, ZQVM, ZQCM, ZQRM, ZQIM, ZQSM, ZQGM, ZQHM, &
                    & ZCPM, ZRHM, ZTM, &
                    & ZAPHIM, ZAPHIFM, &
                    & ZZZ_, ZTENDTT, ZDZZ_, &
                    & ZDTHRAD_, ZICEFR_, ZPRCFR_, &
                    & ZHLC_HRC_, ZHLC_HCF_, ZHLI_HRI_, ZHLI_HCF_)


IF (YDMODEL%YRML_DIAG%YRLDDH%LFLEXDIA) THEN

  CALL ARO_STARTBU(KIDIA, KFDIA, KLON, KFLEVG, &
                 & YDMODEL%YRML_PHY_MF%YRPARAR%NRR, YDMODEL%YRML_GCONF%YGFL%NGFL_EXT, ZRHODJM__(:, 1:KFLEVG), &
                 & ZUS__(:, 1:KFLEVG), ZVS__(:, 1:KFLEVG), &
                 & ZWS__(:, 1:KFLEVG), ZTHS__(:, 1:KFLEVG), &
                 & ZRS_, ZTKES_, YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH)

ENDIF

!     --------------------------------------------------------------------
!      CALCULATE NUMBER OF DROPLETS FROM NRT AEROSOLS
!      FOR MICROPHYSICS AND RADIATION
!     --------------------------------------------------------------------

IF (YDMODEL%YRML_PHY_MF%YRPHY%LAERONRT) THEN
  ! Call function to get the number of condensation nuclei
  CALL ARO_CCN(YDMODEL%YRML_PHY_MF%YRPARAR%PHYEX, &
             & YDMODEL%YRML_PHY_MF%YRPARAR%NRR, &
             & KIDIA, KFDIA, &
             & KLON, KFLEVG, ZDT,&
             & ZRHODREFM__(:, 1:KFLEVG), &
             & ZEXNREFM_, &
             & ZPABSM__(:, 1:KFLEVG), &
             & ZDZZ_, &
             & ZNEBMNH_, &
             & ZTHM__(:, 1:KFLEVG), &
             & ZRM_, &
             & ZWM__(:, 1:KFLEVG), &
             & YDMODEL%YRML_GCONF%YGFL%NAERO, &
             & ZAEROM_, &
             & ZCLDROP_, &
             & ZIFN_NC_)
ENDIF

!    ---------------------------------------------------------------------
!     5 - COMPUTE DUST PROPERTIES FOR RADIATION IF LRDUST=T
!    ---------------------------------------------------------------------

IF (YDMODEL%YRML_PHY_MF%YRARPHY%LRDUST) THEN

  CALL APL_AROME_DUST(YDMODEL, YDCPG_BNDS, YDCPG_OPTS, YDVARS, &
                    & YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG - 2, &
                    & ZZZ_, ZDZZ_, ZPABSM__, &
                    & ZTHM__, ZRHODREFM__, &
                    & ZAERD_, ZTENDEXT, &
                    & ZSVM_, &
                    & ZPIZA_DST_, ZCGA_DST_, ZTAUREL_DST_)

ENDIF ! LRDUST

IF (YDCPG_OPTS%LSFORCS) THEN   ! <== Surface forcing for MUSC

  CALL APL_AROME_SURFACE_FORCING(YDCST, YDMODEL, YDMF_PHYS, YDMF_PHYS_BASE_STATE, YDMF_PHYS_SURF, YDXFU, &
                               & YDCPG_BNDS, YDCPG_OPTS, YDVARS, &
                               & YDGEOMVARS, ZTKEM, &
                               & ZTSURF, &
                               & ZSFTH_, ZSFRV_, ZSFU_, ZSFV_, ZQS)

ENDIF    ! <== End of surface forcing for MUSC

!    ---------------------------------------------------------------------
!     6 - RADIATION
!    ---------------------------------------------------------------------

IF (YDMODEL%YRML_PHY_MF%YRPHY%LRAYFM) THEN

  CALL APL_AROME_RADIATION(YDCST, YDMODEL, YDMF_PHYS_BASE_STATE, YDMF_PHYS_SURF, YDMF_PHYS, YDRADF, YDGEOMETRY, &
                         & YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, YDGPAR, YDVARS, YDGEOMVARS, YDSURF, &
                         & ZSPP_UA, YDDDH, &
                         & LLMSE, &
                         & ZADTMS, &
                         & ZCPM, ZQVM, ZTM, &
                         & ZQIM, ZQSM, ZQGM, ZQCM, ZQRM, ZQHM, &
                         & ZEXNREFM_, ZRDG_MU0, ZRDG_SOLO, ZAERD_, &
                         & ZRDG_MU0M, ZRDG_MU0LU, ZCLDROP_, &
                         & ZPIZA_DST_, ZCGA_DST_, ZTAUREL_DST_, ZAEROM_, &
                         & ZTENDT, &
                         & ZEMIS, ZTSURF, &
                         & ZQSAT, ZALBP, ZALBD, &
                         & ZZS_FSWDIR, ZZS_FSWDIF, &
                         & ZDTHRAD_)
ELSE

  DO JSW=1, YDMODEL%YRML_PHY_RAD%YRERAD%NSW
    DO JLON = KIDIA, KFDIA
      ZALBP(JLON, JSW) = 0._JPRB
      ZALBD(JLON, JSW) = 0._JPRB
    ENDDO
  ENDDO

ENDIF

!    ---------------------------------------------------------------------
!     7 - CONVECTION.
!    ---------------------------------------------------------------------

IF(YDMODEL%YRML_PHY_MF%YRARPHY%LKFBCONV) THEN

  CALL APL_AROME_CONVECTION(YDCST, YDMODEL, YDMF_PHYS_BASE_STATE, &
                          & YDCPG_BNDS, YDCPG_OPTS, YDGEOMETRY, YDVARS, &
                          & YDGEOMVARS, ZPABSM__, ZZZ_F_, ZTM, ZRM_, ZRHODREFM__, ZUM__, ZVM__, ZWM__, &
                          & ZTENDT, &
                          & ZTENDRA(:,:,YDMODEL%YRML_PHY_MF%YRPARAR%YRTENDRA%Q), &
                          & ZTENDRA(:,:,YDMODEL%YRML_PHY_MF%YRPARAR%YRTENDRA%L), &
                          & ZTENDRA(:,:,YDMODEL%YRML_PHY_MF%YRPARAR%YRTENDRA%I), &
                          & ZRS_, ZQDM, ZTHS__, &
                          & ZINPRR_, ZACPRR_, ZINPRS_, ZACPRS_, ZRHODJM__, &
                          & YDDDH)

ENDIF

!    ---------------------------------------------------------------------
!     8 - SURFACE.
!    ---------------------------------------------------------------------

CALL APL_AROME_SURFACE (YDCST, YDMODEL, YDMF_PHYS, YDMF_PHYS_BASE_STATE, YDMF_PHYS_SURF, YDCPG_BNDS, YDCPG_OPTS, &
  & YDGEOMETRY, YDVARS, YDGEOMVARS, YDXFU, YDDDH, ZZZ_F_, ZSVM_, ZUM__, ZVM__, ZTM, ZRM_, ZRHODREFM__, ZPABSM__, ZRDG_MU0, &
  & ZRDG_MU0N, ZINPRR_NOTINCR_, ZINPRS_NOTINCR_, ZINPRG_NOTINCR_, ZTKEM, ZALBD, ZALBP, ZEMIS, ZTSURF, ZZS_FSWDIR, &
  & ZZS_FSWDIF, ZSFTH_, ZSFRV_, ZSFU_, ZSFV_, YDCPG_OPTS%NINDAT, LLMSE, ZQS, ZBUDTH_, ZTOWNS_, ZSFSV_, ZZS_)

!    ---------------------------------------------------------------------
!     9 - Shallow Mass Flux Mixing
!    ---------------------------------------------------------------------

IF (YDMODEL%YRML_PHY_MF%YRARPHY%LMFSHAL) THEN

  CALL APL_AROME_SHALLOW(YDCST, YDMODEL, YDMF_PHYS_BASE_STATE, YDGEOMETRY, YDMF_PHYS, &
                       & YDCPG_BNDS, YDCPG_OPTS, YDVARS, YDGEOMVARS, ZSPP_UA, YDDDH, &
                       & ZSFTH_, ZSFRV_, ZSFU_, ZSFV_, ZEXNREFM_, ZWCLDFR_, ZTKEM__, &
                       & ZTM, ZQVM, ZQCM, ZQIM, ZQDM, &
                       & ZRM_, ZSVM_, &
                       & ZAPHIM, ZAPHIFM, &
                       & ZZZ_, ZZZ_F_, &
                       & ZRHODJM__, ZRHODREFM__, ZPABSM__, &
                       & ZTHM__, ZUM__, ZVM__, &
                       & ZRS_, ZTHS__, ZTKES_, &
                       & ZTENDTKE, ZTENDT, ZTENDTT, &
                       & ZTENDRA(:,:,YDMODEL%YRML_PHY_MF%YRPARAR%YRTENDRA%Q), &
                       & ZTENDRA(:,:,YDMODEL%YRML_PHY_MF%YRPARAR%YRTENDRA%R), &
                       & ZTENDRA(:,:,YDMODEL%YRML_PHY_MF%YRPARAR%YRTENDRA%S), &
                       & ZUS__, ZVS__, ZTKEM, &
                       & ZMF_UP__, &
                       & ZSURFPREP, ZSURFSNOW, &
                       & ZLENGTHM__, ZLENGTHH__, ZTKEEDMFS, &
                       & ZFLXZTHVMF_SUM__, ZFLXZUMF_SUM__, ZFLXZVMF_SUM__, ZSIGMF_)

ENDIF

!    ---------------------------------------------------------------------
!     10 - TURBULENCE.
!    ---------------------------------------------------------------------

IF (YDMODEL%YRML_PHY_MF%YRARPHY%LTURB) THEN

  CALL APL_AROME_TURBULENCE(YDCST, YDMODEL, YDGEOMETRY, YDMF_PHYS_BASE_STATE, YDMF_PHYS, &
                          & YDCPG_BNDS, YDCPG_OPTS, YDVARS, YDGEOMVARS, YDHGRAD, YDXFU, ZSPP_UA, YDDDH, &
                          & ZAPHIM, ZAPHIFM, ZZZ_, ZZZ_F_, &
                          & ZLIMAM_, ZSVM_, ZTKES_, ZSIGMF_, &
                          & ZEXNREFM_, ZQDM, ZTKEEDMFS, ZZS_, &
                          & ZRHODJM__, ZTHVREFM__, ZRHODREFM__, &
                          & ZSFTH_, ZSFRV_, ZSFU_, ZSFV_, &
                          & ZPABSM__, ZUM__, ZVM__, ZWM__, &
                          & ZTKEM__, ZSRCS__, ZTHM__, &
                          & ZRM_, ZUS__, ZVS__, ZWS__, &
                          & ZTHS__, ZRS_, &
                          & ZFLXZTHVMF_SUM__, ZFLXZUMF_SUM__, ZFLXZVMF_SUM__, ZLENGTHM__, ZLENGTHH__,&
                          & ZMF_UP__, ZTENDT, &
                          & ZTENDRA(:,:,YDMODEL%YRML_PHY_MF%YRPARAR%YRTENDRA%Q), &
                          & ZTENDTKE, ZTENDLIMA, ZTENDEXT, ZSFSV_, &
                          & ZLIMAS_, ZSVS_, ZTKES_OUT__, ZSIGS__, &
                          & ZTENDU_TURB__, ZTENDV_TURB__)

ENDIF

IF (YDMODEL%YRML_PHY_MF%YRPHY%LWINDFARM) THEN

    CALL APL_AROME_WINDFARM(YDGEOMETRY, &
                          & YDMODEL, &
                          & YDCPG_BNDS, &
                          & YDCPG_OPTS, &
                          & YDMF_PHYS_BASE_STATE, &
                          & YDVARS, &
                          & YDGEOMVARS, &
                          & YDMF_PHYS, &
                          & ZTENDTKE)
ENDIF

!    ---------------------------------------------------------------------
!     11 - MICROPHYSIQUE.
!    ---------------------------------------------------------------------

IF (YDMODEL%YRML_PHY_MF%YRARPHY%LMICRO) THEN

  CALL APL_AROME_MICRO(YDMODEL, YDCPG_BNDS, YDCPG_OPTS, &
                     & YDMF_PHYS, YDMF_PHYS_SURF, &
                     & ZSPP_UA, YDDDH, &
                     & ZTOWNS_, &
                     & ZWM__, ZTKEM__, &
                     & ZDT, ZDZZ_, ZRHODJM__, ZRHODREFM__, &
                     & ZEXNREFM_, ZPABSM__, ZTHM__, &
                     & ZDTHRAD_, ZSIGS__, &
                     & ZICLDFR_, ZWCLDFR_, &
                     & ZSSIO_, ZSSIU_, ZIFR_, &
                     & ZQDM, ZZZ_, &
                     & ZIFN_NC_, ZCLDROP_, ZAEROM_, &
                     & ZRM_, ZLIMAM_, &
                     & ZHLC_HRC_, ZHLC_HCF_, ZHLI_HRI_, ZHLI_HCF_, &
                     & ZICEFR_, ZPRCFR_, &
                     & ZINPRG_NOTINCR_, ZINPRR_NOTINCR_, ZINPRS_NOTINCR_, &
                     & ZINPRG_, ZINPRR_, ZINPRS_, ZINPRH_, &
                     & ZNEBMNH_, ZTENDT, ZTENDTT, &
                     & ZTENDRA, ZTENDLIMA, ZTENDEXT, &
                     & ZTHS__(:, 1:KFLEVG), ZRS_, ZLIMAS_, ZSVM_, &
                     & ZTENDAERO)
ENDIF ! LMICRO

CALL APL_AROME_NUDGE(YDCST, YDMODEL, YDMF_PHYS_BASE_STATE, &
                   & YDCPG_BNDS, YDCPG_OPTS, YDGEOMETRY, &
                   & ZDT, &
                   & ZEXNREFM_, ZPABSM__, ZTHM__, ZRM_, ZQDM, ZQSAT, &
                   & ZTENDT, ZTENDTT, &
                   & ZTENDRA(:,:,YDMODEL%YRML_PHY_MF%YRPARAR%YRTENDRA%Q), &
                   & ZACPRR_, ZACPRS_, ZACPRG_)

!    ---------------------------------------------------------------------
!     12 - STOCHASTIC PHYSICS : PERTURB TENDENCIES
!    ---------------------------------------------------------------------

CALL APL_AROME_SPPT(YDCST, YDMODEL, YDCPG_BNDS, YDCPG_OPTS, YDVARS, &
                  & YDFORCESPPT, YDMF_PHYS, YDMF_PHYS_BASE_STATE, &
                  & PGP2DSDT(:,1,1), &
                  & ZTENDT, &
                  & ZTENDRA(:,:,YDMODEL%YRML_PHY_MF%YRPARAR%YRTENDRA%Q))

!    ---------------------------------------------------------------------
!     13 - FINAL CALCULATIONS.
!    ---------------------------------------------------------------------

CALL APL_AROME_FINAL(YDCST, YDMODEL, YDCPG_BNDS, YDCPG_OPTS, YDVARS, &
                   & YDGEOMVARS, YDCFU, YDXFU, YDMF_PHYS_BASE_STATE, YDMF_PHYS_SURF, &
                   & YDGPAR, YDMF_PHYS, &
                   & YDCPG_OPTS%LSFORCS, &
                   & ZDT, &
                   & ZINPRR_, ZINPRS_, ZINPRG_, ZINPRH_, &
                   & ZTSURF, ZEMIS, ZQS, &
                   & ZALBD, ZALBP, &
                   & ZSURFPREP, ZSURFSNOW, &
                   & ZSFU_, ZSFV_, &
                   & ZRHODREFM__, ZTM, ZSFTH_, ZSFRV_, &
                   & ZQCM, ZQIM, ZQRM, ZQSM, ZQGM, ZQHM, ZQVM, &
                   & ZRHM, ZBUDTH_, ZAPHIM, ZAPHIFM, ZZZ_, &
                   & ZTENDT)

IF (YDMODEL%YRML_PHY_FORCING%LMUSCLFA) THEN
  CALL ECR1D(YDMODEL%YRML_PHY_FORCING%NMUSCLFA, 'PCLCT_apl', YDCPG_MISC%CLCT, 1, KLON)
ENDIF

IF (YDMODEL%YRML_DIAG%YRLDDH%LFLEXDIA) THEN

  CALL APL_AROME_DDH(YDCPG_BNDS, YDCPG_OPTS, YDMODEL, YDMF_PHYS, YDMF_PHYS_SURF, YDVARS, YDGEOMVARS, YDDDH, &
                   & ZTSURF)

ENDIF

!      4.2  COMPUTE THE PHYS. TENDENCY FOR "T" AND "w"
!           ------------------------------------------

CALL APL_AROME_TW_TENDENCY(YDCST, YDMODEL, YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_DYN0, YDVARS, &
                         & YDMF_PHYS_BASE_STATE, YDMF_PHYS, &
                         & YDMF_PHYS_NEXT_STATE, &
                         & ZTENDT, ZTENDW, &
                         & ZTENDRA, &
                         & ZTENDLIMA, ZTENDAERO, ZTENDTKE, ZTENDEXT)

!     ------------------------------------------------------------------
!*       5.    Final calculations.
!              -------------------

! * Restore the initial value of some pseudo-historical surface buffers
!   if relevant.
IF (.NOT. YDMODEL%YRML_PHY_EC%YREPHY%LAGPHY .AND. YDCPG_OPTS%LCONFX) THEN
    CALL MF_PHYS_SAVE_PHSURF_PART2 (YDCPG_BNDS, YDCPG_OPTS, ZSAV_DDAL, ZSAV_DDOM, ZSAV_ENTCH,                              &
    & ZSAV_FHPS, ZSAV_GZ0F, ZSAV_GZ0HF, ZSAV_HV, ZSAV_PBLH, ZSAV_QSH, ZSAV_UDAL, ZSAV_UDGRO,                               &
    & ZSAV_UDOM, ZSAV_UNEBH, ZSAV_MO, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VH%PPBLH, YDMF_PHYS_SURF%GSD_VH%PQSH, &
    & YDMF_PHYS_SURF%GSD_VH%PSPSH, YDMF_PHYS_SURF%GSD_VK%PUDGRO, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VV%PZ0H,    &
    & YDVARS%DAL%T0, YDVARS%DOM%T0, YDVARS%UAL%T0, YDVARS%UEN%T0, YDVARS%UNEBH%T0, YDVARS%UOM%T0, YDVARS%MO%T0,            &
    & YDMODEL)

  ! make radiative cloudiness available for postprocessing (in zero X step)
  CALL TRANSFER_CLOUD_FIELDS(YDMODEL, YDVARS)
ENDIF

!-------------------------------------------------
! Extract Single Column Model profiles from 3D run or
! write LFA file for MUSC (1D model)
!-------------------------------------------------

IF(YDMODEL%YRML_PHY_SCM%LGSCM .OR. YDMODEL%YRML_PHY_FORCING%LMUSCLFA) THEN
  IF (YDCPG_OPTS%LAROME) THEN
    DO JLEV=1,KFLEVG
      DO JLON=KIDIA,KFDIA
        YDCPG_MISC%NEB(JLON,JLEV)=YDVARS%A%T1(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF

  CALL WRITEPHYSIO (YDGEOMETRY, YDCPG_MISC, YDCPG_PHY0, YDMF_PHYS, YDRADF, YDCPG_DYN0, YDMF_PHYS_SURF, YDVARS, YDSURF, &
    & YDMODEL%YRML_PHY_G%YRDPHY, YDMODEL%YRML_GCONF%YRRIP, YDMODEL%YRML_PHY_MF, YDMODEL%YRML_PHY_FORCING, &
    & YDMODEL%YRML_PHY_SCM, KFDIA, KIDIA, YDCPG_OPTS%KGL1, YDCPG_OPTS%KGL2, YDCPG_OPTS%NSTEP, &
    & YDMODEL%YRML_PHY_G%YRDPHY%NTSSG, YDSURF%YSP_SBD%NLEVS, YDGEOMVARS%GELAM, YDGEOMVARS%GEMU, YDGEOMVARS%GM, &
    & YDGEOMVARS%OROG, YDGEOMVARS%RCORI, YDGEOMVARS%RATATH, YDGEOMVARS%RATATX, YDGEOMVARS%GECLO, YDGEOMVARS%GESLO, ZRDG_CVGQ, &
    & ZRDG_LCVQ, ZRDG_MU0)

ENDIF

IF (YDMODEL%YRML_PHY_MF%YRPHY%LEDR) THEN
  YDMF_PHYS_SURF%GSD_DI%PXEDR(:,:) = YDMF_PHYS%EDR(:,:)
ENDIF

IF (YDMODEL%YRML_PHY_MF%YRPHY%LGCHECKNAN) THEN

  CALL CHECKNANARO (YDMODEL%YRML_GCONF%YRRIP, YDMODEL%YRML_PHY_MF%YRPHY0, YDMODEL%YRML_PHY_MF%YRPHY2, KIDIA, KFDIA, KLON, &
    & YDMODEL%YRML_PHY_MF%YRPARAR%NRR, KFLEVG, YDCPG_OPTS%NSTEP, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, &
    & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF, &
    & YDGEOMVARS%GELAM, YDGEOMVARS%GEMU, ZRDG_MU0, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_BASE_STATE%T%P, &
    & YDMF_PHYS_BASE_STATE%Q%P, ZQSAT, ZTSURF, ZTENDT, ZTENDRA, ZTENDW, ZTENDTKE, YDMF_PHYS%FRTH, YDMF_PHYS%FRSO, &
    & YDMF_PHYS%FRTHDS, YDMF_PHYS%FRSODS, YDMF_PHYS%CLCH, YDMF_PHYS%CLCM, YDMF_PHYS%CLCL, &
    & YDCPG_MISC%CLCT, YDMF_PHYS%FPLSL, YDMF_PHYS%FPLSN, YDMF_PHYS%FPLSG, YDMF_PHYS%STRTU, &
    & YDMF_PHYS%STRTV, YDMF_PHYS%FCS, YDMF_PHYS%FCLL, YDMF_PHYS%FCLN, YDMF_PHYS%FEVL, YDMF_PHYS%FEVN)

ENDIF

! Clear SPP
IF (YDMODEL%YRML_GCONF%YRSPP_CONFIG%LSPP) CALL ZSPP_UA%CLEAR()

!     ------------------------------------------------------------------

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('APL_AROME',1,ZHOOK_HANDLE)

END SUBROUTINE APL_AROME

