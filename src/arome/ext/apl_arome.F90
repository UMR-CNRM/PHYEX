#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE APL_AROME(YDCST, YDMF_PHYS_BASE_STATE, YDMF_PHYS_NEXT_STATE, YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, &
& YDCPG_MISC, YDCPG_GPAR, YDCPG_PHY0, YDMF_PHYS, YDCPG_DYN0, YDMF_PHYS_SURF, &
& YDCPG_SL1, YDVARS, YDGMV, YDSURF, YDCFU, YDXFU, YDMODEL, PGFL, PGP2DSDT, PGMVT1,  &
& PGFLT1, PTRAJ_PHYS, YDDDH)

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

! KMAXDRAFT : MAX NUMBER OF DRAFTS (FOR DIMENSIONNING)
! KSGST : NUMBER OF SUBGRID SURFACE TEMPERATURES AND FLUXES (NTSSG IN *CPG*)
! KNFRRC : FREQUENCY FOR CLEAR SKY RADIATION CALCULATION
! PDT : TIME STEP (in s)
! LDXFUMSE : T if CDCONF=X in order not to increment surfex timer in that case
!-----------------------------------------------------------------------
! YDVARS%GEOMETRY%GEMU%T0      : SINE OF GEOGRAPHICAL LATITUDE
! PGELAM     :  LONGITUDE
! POROG      : g * OROGRAPHY
! PGM        : MAP FACTOR (used in ALARO convection only)
! PCLON      : cosine of geographical longitude.
! PSLON      : sine of geographical longitude.
! PGP2DSDT   : STOCHASTIC PHYSICS PATTERNS

! FIELDS WITH SUBSCRIPT M FOR TIME T-DT IN 3TL OR T IN 2TL

! PDELPM     : LAYER THICKNESS IN PRESSURE UNITS

! PTM        : TEMPERATURE.
! PCPM        : SPECIFIC HEAT AT CONSTANT PRESSURE FOR AIR
! PRM         : GAS CONSTANT FOR AIR

! PTKEM       : TURBULENT KINETIC ENERGY
! PSVM        : PASSIVE SCALARS
! PUM         : ZONAL WIND
! PVM         : MERIDIAN WIND
! PWM         : VERTICAL VELOCITY (m/s)

!-----------------------------------------------------------------------
! - INOUT

! PGPAR       : BUFFER FOR 2D FIELDS - CONTAINS PRECIP, ALBEDO, EMISS, TS
!             : SURFACE FLUXES

! ACRANEB2 intermittency storage

!               LINEAR T_e CORRECTION
!               LINEAR T_e CORRECTION
!-----------------------------------------------------------------------
! - OUTPUT (SUBSCRIPT S FOR T+DT)

! PSIGS       : SIGMA FOR SUBGRIDCOND
! PTENDT      : TEMPERATURE TENDENCY
! PTENDR      : HYDROMETEORE TENDENCIES
! PTENDW      : VERTICAL VELOCITY TENDENCY
! PTENDTKE    : TKE TENDENCY
! PTENDEXT    : PASSIVE SCALARS TENDENCY
! PFRSO       : SHORTWAVE RADIATIVE FLUX
! - 2D (0:1)

! variables used in input for radiation in case no surface scheme is used

! PSIC       : MODEL SEA ICE CONCENTRATION

! Part of GFL strcture dedicated to easy diagnostics (to be used as a print...)
! PEZDIAG    : MULPITPLE ARRAY TO BE FILLED BY THE USER BY 3D FIELDS
!              (NGFL_EZDIAG ONES)
! output for CFU XFU

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
!     2021-12, F. Meier add call to windfarm parametrisation if LWINDFARM is TRUE
!     R. El Khatib 08-Jul-2022 Contribution to the encapsulation of YOMCST and YOETHF
! End modifications
!-------------------------------------------------------------------------------


USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_MISC_TYPE, CPG_GPAR_TYPE, CPG_SL1_TYPE, CPG_DYN_TYPE, CPG_PHY_TYPE
USE CPG_OPTS_TYPE_MOD   , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMXFU             , ONLY : TXFU
USE YOMCFU             , ONLY : TCFU
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPIM     ,JPRB     ,JPRD
USE YOMHOOK            , ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMCST     , ONLY : TCST
USE YOMLUN     , ONLY : NULOUT
USE YOMCT0     , ONLY : LSFORCS, LAROME
USE YOMVERT    , ONLY : VP00
USE YOMRIP0    , ONLY : NINDAT
USE YOMNUDGLH , ONLY :  LNUDGLH, NSTARTNUDGLH, NSTOPNUDGLH, NINTNUDGLH, NTAUNUDGLH, &
          &             RAMPLIFY,RMAXNUDGLH,RMINNUDGLH,LNUDGLHCOMPT,NTIMESPLITNUDGLH
USE YOMNSV     , ONLY : NSV_CO2
USE DDH_MIX    , ONLY : NEW_ADD_FIELD_3D, NEW_ADD_FIELD_2D,&
                    & NTOTSVAR, NTOTSURF, NTOTSVFS, TYP_DDH ! for new data flow
!USE SPP_MOD    , ONLY : YSPP_CONFIG, YSPP
USE SPP_MOD_TYPE, ONLY : ALL_SPP_VARS, SET_ALL_SPP, CLEAR_ALL_SPP, APPLY_SPP
USE YOMLSFORC  , ONLY : LMUSCLFA, NMUSCLFA, REMIS_FORC, RALB_FORC
USE INTFLEX_MOD, ONLY : LINTFLEX, LRADFLEX,&
                      & TYPE_INTPROC, TYPE_INTPROCSET,&
                      & NEWINTFIELD, NEWINTPROC, NEWINTPROCSET, CLEANINTPROCSET
USE YOMGFL     , ONLY : GFL_WKA, GFL_WKA2
USE YOMMP0     , ONLY : MYPROC     
USE MF_PHYS_BASE_STATE_TYPE_MOD &
             & , ONLY : MF_PHYS_BASE_STATE_TYPE
USE MF_PHYS_NEXT_STATE_TYPE_MOD &
             & , ONLY : MF_PHYS_NEXT_STATE_TYPE
USE YOMGMV             , ONLY : TGMV
USE SC2PRG_MOD         , ONLY : SC2PRG

USE YOMSCM             , ONLY : LGSCM
USE YOMTRAJ            , ONLY : TRAJ_PHYS_TYPE

!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST),                     INTENT(IN)    :: YDCST
TYPE (MF_PHYS_BASE_STATE_TYPE), INTENT(IN)    :: YDMF_PHYS_BASE_STATE
TYPE (MF_PHYS_NEXT_STATE_TYPE), INTENT(INOUT) :: YDMF_PHYS_NEXT_STATE
TYPE(GEOMETRY),                 INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE),            INTENT(IN)    :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),            INTENT(IN)    :: YDCPG_OPTS
TYPE(CPG_MISC_TYPE),            INTENT(INOUT) :: YDCPG_MISC
TYPE(CPG_GPAR_TYPE),            INTENT(INOUT) :: YDCPG_GPAR
TYPE(CPG_PHY_TYPE),             INTENT(IN)    :: YDCPG_PHY0
TYPE(MF_PHYS_TYPE),             INTENT(IN)    :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),             INTENT(IN)    :: YDCPG_DYN0
TYPE(MF_PHYS_SURF_TYPE),        INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(CPG_SL1_TYPE),             INTENT(INOUT) :: YDCPG_SL1
TYPE(FIELD_VARIABLES),          INTENT(INOUT) :: YDVARS
TYPE(TGMV),                     INTENT(IN)    :: YDGMV
TYPE(TSURF),                    INTENT(IN)    :: YDSURF
TYPE(TCFU),                     INTENT(IN)    :: YDCFU
TYPE(TXFU),                     INTENT(IN)    :: YDXFU
TYPE(MODEL),                    INTENT(IN)    :: YDMODEL

 
REAL(KIND=JPRB),                INTENT(INOUT) :: PGFL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB),                INTENT(IN)    :: PGP2DSDT(YDCPG_OPTS%KLON,YDMODEL%YRML_SPPT%YGPSDT(1)%NG2D,YDMODEL%YRML_SPPT%N2D)
REAL(KIND=JPRB),                INTENT(INOUT) :: PGMVT1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB),                INTENT(INOUT) :: PGFLT1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
TYPE (TRAJ_PHYS_TYPE),          INTENT(INOUT) :: PTRAJ_PHYS

TYPE(TYP_DDH),                  INTENT(INOUT) :: YDDDH

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
! These arrays are passed in argument as ZXXX__(:,1:YDCPG_OPTS%KFLEVG) except for aro_turb_mnh where they are passed as ZXXX__.

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
!   Z(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=value
! 2D arrays :
!  DO JLEV=1,YDCPG_OPTS%KFLEVG
!     ZX(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=xval
!     ZY(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=yval
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
!      ZSUM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZSUM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)+ZINC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
!    ENDIF
!  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(KIND=JPRB) :: ZRHODJM__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),      ZRHODREFM__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),   ZPABSM__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZUM__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),          ZVM__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),         ZTHM__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZUS__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),          ZVS__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),         ZWS__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZTKES_OUT__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),    ZMF_UP__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),      ZTHVREFM__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1)  ! thetav de l etat
REAL(KIND=JPRB) :: ZTENDU_TURB__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),  ZTENDV_TURB__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1), ZTENDTHL_TURB__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZTENDRT_TURB__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1), ZTKEM__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),       ZSRCS__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1) 
REAL(KIND=JPRB) :: ZSIGS__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),        ZEDR__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1)
! THE DDH budgets
REAL(KIND=JPRB) :: ZDP__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),          ZTP__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),         ZTPMF__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZTDIFF__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),       ZTDISS__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1)
! length scales for momentum and heat for mnh level definitions in case LHARATU=TRUE
REAL(KIND=JPRB) :: ZLENGTHM__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1), ZLENGTHH__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1)

REAL(KIND=JPRB), POINTER :: ZTHS__(:,:)
! horizontal gradients and diagnostics
REAL(KIND=JPRB) :: ZTURB3D__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1,YDMODEL%YRML_PHY_MF%YRARPHY%NGRADIENTS)
! WARNING ! Don't use ZTHSWAP__ or ZTHSAVE__ below because they may be swapped !
! Use only the pointer ZTHS__, and possibly ZTHSIN_ if you need the backup of input data.
REAL(KIND=JPRB), TARGET  :: ZTHSWAP__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1),        ZTHSAVE__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB), TARGET  :: ZFLXZTHVMF_SUM__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1), ZWM__(YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG+1)


! Updraft characteristics for Meso-NH world (input of ARO_SHALLOW_MF)
REAL(KIND=JPRB) :: ZTHETAL_UP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG), ZTHETAV_UP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG), ZZFRAC_UP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZRT_UP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),     ZRC_UP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),     ZRI_UP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZZU_UP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),     ZZV_UP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),     ZZW_UP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZZRV_UP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),    ZTKES_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),      ZZZ_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZDZZ_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),       ZZZ_F_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),      ZDZZ_F_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZCIT_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),       ZMFM_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),       ZEXNREFM_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZSIGM_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),      ZNEBMNH_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),    ZEVAP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
! additions for future ice cloud fraction and precipitation fraction
REAL(KIND=JPRB) :: ZICEFR_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),     ZPRCFR_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
! additions for MF scheme (Pergaud et al)
REAL(KIND=JPRB) :: ZSIGMF_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),     ZRC_MF_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),     ZRI_MF_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZCF_MF_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),     ZAERD_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),      ZCVTENDT_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZCVTENDRV_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),  ZCVTENDRC_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),  ZCVTENDRI_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZMFS_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),       ZTHLS_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),      ZRTS_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZMFUS_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),      ZMFVS_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),      ZDEPTH_HEIGHT_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZDTHRAD_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB), TARGET :: ZFLXZTHVMF_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), POINTER :: ZARG_FLXZTHVMF_(:,:)

! Subgrid autoconversions
REAL(KIND=JPRB) :: ZHLC_HRC_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),     ZHLC_HCF_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG), &
                   ZHLI_HRI_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),     ZHLI_HCF_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)


! WARNING ! Don't use ZRSWAP_ or ZRSAVE_ below because they may be swapped !
! Use only the pointer ZRS_, and possibly ZRSIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZRSIN_(:,:,:), ZRS_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZRSWAP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
REAL(KIND=JPRB), TARGET :: ZRSAVE_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

REAL(KIND=JPRB), POINTER  :: ZPTRWNU_(:,:), ZTHSIN_(:,:)
REAL(KIND=JPRB), TARGET :: ZWNU_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)

! WARNING ! Don't use ZSVSWAP_ or ZSVSAVE_ below because they may be swapped !
! Use only the pointer ZSVS_, and possibly ZSVSIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZSVSIN_(:,:,:), ZSVS_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZSVSWAP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB), TARGET :: ZSVSAVE_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

REAL(KIND=JPRB) :: ZSVXXX_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

! WARNING ! Don't use ZSVMSWAP_ or ZSVMSAVE_ below because they may be swapped !
! Use only the pointer ZSVM_, and possibly ZSVMIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZSVMIN_(:,:,:), ZSVM_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZSVMSWAP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) 
REAL(KIND=JPRB), TARGET :: ZSVMSAVE_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB) :: ZSVMB_(YDCPG_BNDS%KFDIA,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

! WARNING ! Don't use ZLIMASWAP_ or ZLIMASAVE_ below because they may be swapped !
! Use only the pointer ZLIMAS_, and possibly ZLIMASIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZLIMAS_(:,:,:), ZLIMASIN_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZLIMASWAP_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NLIMA)
REAL(KIND=JPRB), TARGET :: ZLIMASAVE_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NLIMA)

REAL(KIND=JPRB) :: ZLIMAM_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NLIMA)
!INTEGER(KIND=JPIM) :: KSV_TURB !CPtoclean?
!CPtoclean REAL(KIND=JPRB) :: ZTURBM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA)
!CPtoclean REAL(KIND=JPRB) :: ZTURBS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA)
!not (yet ?) used. REK
!REAL(KIND=JPRB) :: ZSFTURB(YDCPG_OPTS%KLON,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA) ! surf. flux of SV (=0)
!REAL(KIND=JPRB) :: ZTENDSV_TURB2(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA)  ! SV (=0)
REAL(KIND=JPRB) :: ZSFSVLIMA_(YDCPG_BNDS%KFDIA,YDMODEL%YRML_GCONF%YGFL%NLIMA) ! surf. flux of LIMA vars
REAL(KIND=JPRB) :: ZTENDSV_TURBLIMA_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NLIMA) ! LIMA

! For radiation scheme
REAL(KIND=JPRB) :: ZRM_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
REAL(KIND=JPRB) :: ZPFPR_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

REAL(KIND=JPRB) :: ZPEZDIAG_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG)

REAL(KIND=JPRB) :: ZSFSV_(YDCPG_BNDS%KFDIA,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! surf. flux of scalars
REAL(KIND=JPRD) :: ZGEMU_D(YDCPG_OPTS%KLON)  ! double precision version of YDVARS%GEOMETRY%GEMU%T0, for RADACT

! Single scattering albedo of dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZPIZA_DST_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)
! Assymetry factor for dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZCGA_DST_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)
! tau/tau_{550} dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZTAUREL_DST_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)


! surface flux of theta and surface flux of vapor ; surface flux of CO2
REAL(KIND=JPRB) :: ZSFTH_(YDCPG_BNDS%KFDIA),  ZSFRV_(YDCPG_BNDS%KFDIA),          ZSFCO2_(YDCPG_BNDS%KFDIA)
REAL(KIND=JPRB) :: ZACPRG_(YDCPG_BNDS%KFDIA), ZINPRG_NOTINCR_(YDCPG_BNDS%KFDIA), ZINPRG_(YDCPG_BNDS%KFDIA)
REAL(KIND=JPRB) :: ZACPRR_(YDCPG_BNDS%KFDIA), ZINPRR_NOTINCR_(YDCPG_BNDS%KFDIA), ZINPRR_(YDCPG_BNDS%KFDIA)
REAL(KIND=JPRB) :: ZACPRS_(YDCPG_BNDS%KFDIA), ZINPRS_NOTINCR_(YDCPG_BNDS%KFDIA), ZINPRS_(YDCPG_BNDS%KFDIA)
REAL(KIND=JPRB) :: ZCFBTH_(YDCPG_BNDS%KFDIA), ZINPRH_NOTINCR_(YDCPG_BNDS%KFDIA), ZINPRH_(YDCPG_BNDS%KFDIA)
REAL(KIND=JPRB) :: ZZS_(YDCPG_BNDS%KFDIA),    ZSSO_STDEV_(YDCPG_BNDS%KFDIA),     ZALB_UV_(YDCPG_BNDS%KIDIA)
REAL(KIND=JPRB) :: ZLAT_(YDCPG_BNDS%KIDIA),   ZLON_(YDCPG_BNDS%KIDIA),           ZZENITH_(YDCPG_BNDS%KIDIA)
REAL(KIND=JPRB) :: ZGZ0_(YDCPG_BNDS%KFDIA),   ZGZ0H_(YDCPG_BNDS%KFDIA),          ZTOWNS_(YDCPG_BNDS%KFDIA) 
REAL(KIND=JPRB) :: ZCFAQ_(YDCPG_BNDS%KFDIA),  ZCFBQ_(YDCPG_BNDS%KFDIA),          ZCFATH_(YDCPG_BNDS%KFDIA)
REAL(KIND=JPRB) :: ZCFAU_(YDCPG_BNDS%KFDIA),  ZCFBU_(YDCPG_BNDS%KFDIA),          ZCFBV_(YDCPG_BNDS%KFDIA)
REAL(KIND=JPRB) :: ZBUDTH_ (YDCPG_BNDS%KFDIA),ZBUDSO_(YDCPG_BNDS%KFDIA),         ZFCLL_(YDCPG_BNDS%KFDIA)
REAL(KIND=JPRB) :: ZCD_(YDCPG_BNDS%KFDIA),    ZSEA_(YDCPG_BNDS%KFDIA),           ZTOWN_(YDCPG_BNDS%KFDIA)
REAL(KIND=JPRB) :: ZZTOP_(YDCPG_BNDS%KFDIA),  ZCVTENDPR_(YDCPG_BNDS%KFDIA),      ZCVTENDPRS_(YDCPG_BNDS%KFDIA)
! surface flux of x and y component of wind. are they really necessary ? REK
REAL(KIND=JPRB) :: ZSFU_(YDCPG_BNDS%KFDIA),   ZSFV_(YDCPG_BNDS%KFDIA)

! local areas for OCND2 option:
! ZICLDFR = ice cloud fraction , ZWCLDFR = water or mixed-phase cloud fraction,
! ZSSIO = Super-saturation with respect to ice in ZICLDFR , 
! ZSSIU = Sub-saturation with respect to ice outside ZICLDFR,
! ZIFR = variable used for calulation of subgridscale ice
! Meso-NH world
REAL(KIND=JPRB) :: ZICLDFR_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),    ZWCLDFR_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZSSIO_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),      ZSSIU_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),      ZIFR_(YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)

! Arpege-style dimensionning :
! --------------------------

!Variables used in case LHARATU=TRUE
! length scales for momentum and heat and TKE
REAL(KIND=JPRB) :: ZLENGTH_M(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZLENGTH_H(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTKEEDMF(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTKEEDMFS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZEMIS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZTMP2(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZTMP(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZQICE(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQLIQ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZAER(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,6)
REAL(KIND=JPRB) :: ZAERINDS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZRAER(YDCPG_OPTS%KLON,6,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZAERO(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,12)


REAL(KIND=JPRB) :: ZQSAT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZFRSOFS(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZLH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZLSCPE(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZGEOSLC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQDM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZQCO2(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQCH4(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQN2O(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQNO2(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQC11(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQC12(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQC22(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQCL4(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZCHTIX(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG+1), ZCAPH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG+1), ZTH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG+1)
REAL(KIND=JPRB) :: ZDUM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZGELAM(YDCPG_OPTS%KLON)

REAL(KIND=JPRB) :: ZTKEM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQW(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZTW(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZTENT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDTT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! array to save heating profile for LHN
REAL(KIND=JPRB) :: ZMAXTEND,ZMINTEND
REAL(KIND=JPRB) :: ZDZZ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTPW(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

! POUR GROUND 
REAL(KIND=JPRB) :: ZZS_FSWDIR(YDCPG_OPTS%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZZS_FSWDIF(YDCPG_OPTS%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZTRSODIR (YDCPG_OPTS%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZTRSODIF (YDCPG_OPTS%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZTP2(YDCPG_OPTS%KLON), ZWS2(YDCPG_OPTS%KLON), ZWP2(YDCPG_OPTS%KLON), ZWSI2(YDCPG_OPTS%KLON), ZWPI2(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZWR2(YDCPG_OPTS%KLON), ZSNA2(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZTRSOD(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZSUDU(YDCPG_OPTS%KLON), ZSDUR(YDCPG_OPTS%KLON), ZDSRP(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZCEMTR(YDCPG_OPTS%KLON,2), ZCTRSO(YDCPG_OPTS%KLON,2)

REAL(KIND=JPRB) :: ZALBD(YDCPG_OPTS%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZALBP(YDCPG_OPTS%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZALBD1(YDCPG_OPTS%KLON), ZALBP1(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZAPHIM(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG), ZAPHIFM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZTM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQVM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQIM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQCM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQHM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQHGM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZQRM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQSM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZQGM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZUPGENL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZUPGENN(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZCLFR(YDCPG_OPTS%KLON)

REAL(KIND=JPRB) :: ZCPM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZRHM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

! Variables concerning updraft rain/snow for EDMF
REAL(KIND=JPRB) :: ZTENDTUP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZTENDQVUP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

! specific to new data flow for diagnostics
REAL(KIND=JPRB) :: ZTENDTBAK(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZTENDRBAK(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
REAL(KIND=JPRB) :: ZTMPAF(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

! daand: radflex
REAL(KIND=JPRB)  :: ZFPR(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

! Target should not be necessary. REK
REAL(KIND=JPRB), TARGET :: ZCON1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), TARGET :: ZCON2(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB), TARGET :: ZCON3(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

! Ajout pour MF Dual Scheme (KNMI et al)
! Updraft characteristics in Arpege/IFS world
REAL(KIND=JPRB) :: ZMF_UP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KMAXDRAFT)
REAL(KIND=JPRB) :: ZTHETAL_UP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KMAXDRAFT)
REAL(KIND=JPRB) :: ZQT_UP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KMAXDRAFT)
REAL(KIND=JPRB) :: ZTHTV_UP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KMAXDRAFT)
REAL(KIND=JPRB) :: ZQC_UP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KMAXDRAFT)
REAL(KIND=JPRB) :: ZQI_UP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KMAXDRAFT)
REAL(KIND=JPRB) :: ZU_UP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KMAXDRAFT)
REAL(KIND=JPRB) :: ZV_UP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KMAXDRAFT)
REAL(KIND=JPRB) :: ZTSURF(YDCPG_OPTS%KLON), ZTN(YDCPG_OPTS%KLON), ZQS(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZZEXNREFM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG), ZZWCLDFR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZFRSOLU(YDCPG_OPTS%KLON), ZFRSODS(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZFSDNN(YDCPG_OPTS%KLON), ZFSDNV(YDCPG_OPTS%KLON)

REAL(KIND=JPRB) :: ZSURFPREP(YDCPG_OPTS%KLON), ZSURFSNOW(YDCPG_OPTS%KLON)

REAL(KIND=JPRB) :: ZQO3(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZZS_FTH_(YDCPG_OPTS%KLON), ZZS_FRV_(YDCPG_OPTS%KLON), ZZS_FU_(YDCPG_OPTS%KLON), ZZS_FV_(YDCPG_OPTS%KLON)

! Surface forcing arrays for MUSC
REAL(KIND=JPRB) :: ZRHODREFM(YDCPG_OPTS%KLON), ZTHETAS(YDCPG_OPTS%KLON)

! ACRANEB2 local variables
REAL(KIND=JPRB) :: ZNEB0    (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)  ! protected cloud fractions
REAL(KIND=JPRB) :: ZCLCT_RAD(YDCPG_OPTS%KLON)       ! total cloud cover for radiation
REAL(KIND=JPRB) :: ZDECRD   (YDCPG_OPTS%KLON)       ! decorrelation depth

! Stochastic physics pattern & dummy tendencies for calling sppten
! Bof. REK
REAL(KIND=JPRB) :: ZMULNOISE(YDCPG_OPTS%KLON,1)
REAL(KIND=JPRB) :: ZDUMMY(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZDUMMY1(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: PTENDENCYU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,1),PTENDENCYV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,1)
REAL(KIND=JPRB) :: PTENDENCYT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,1),PTENDENCYQ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,1)

REAL(KIND=JPRB) :: ZROZ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

! Can we remove ? REK
REAL(KIND=JPRB) :: ZEPSM(0,0,0) ! Dissipation of TKE (eps) at time t-dt
REAL(KIND=JPRB) :: ZEPSS(0,0,0) ! Dissipation of TKE at time t+dt


!    Integers
INTEGER(KIND=JPIM) :: JLEV, JLON, JRR, JGFL, JGR
INTEGER(KIND=JPIM) :: IJN  ! max. number of day/night slices within NRPOMA
INTEGER(KIND=JPIM) :: IKL  !ordering of vert levels 1:MNH -1:AROME
INTEGER(KIND=JPIM) :: IOFF_MFSHAL, IEZDIAG_CHEM
INTEGER(KIND=JPIM) :: IKA,IKB,IKU,IKT,IKTE,IKTB ! vertical points as in mpa
INTEGER(KIND=JPIM) :: JSG, JK, JR, JSW, JAE
INTEGER(KIND=JPIM) :: IDRAFT,JDRAFT,INDRAFT
INTEGER(KIND=JPIM) :: ISURFEX
INTEGER(KIND=JPIM) :: IDAY,IYEAR,IMONTH,IAERO

INTEGER(KIND=JPIM) :: INIT0 ! Kind of safety/debugging initialization :
                            ! 0 = initialize to HUGE (debugging)
                            ! 1 = initialize to realistic value (discouraged)
                            ! -1 = no initialization (optimized code) - this is the default.

INTEGER(KIND=JPIM) :: ICLPH(YDCPG_OPTS%KLON)             !PBL top level
INTEGER(KIND=JPIM) :: JLHSTEP,ISTEP

!       Real
REAL(KIND=JPRB) :: ZRHO
REAL(KIND=JPRB) :: ZAEO, ZAEN, ZSALBCOR
REAL(KIND=JPRB) :: ZDT, ZDT2, ZINVDT, ZINVG, ZRSCP, ZINVATM, Z_WMAX, Z_WMIN
 ! pas de temps pour la surface externalise
REAL(KIND=JPRB) :: ZDTMSE,ZRHGMT,ZSTATI,ZADTMS
REAL(KIND=JPRB) :: ZDELTA
REAL(KIND=JPRB) :: ZEPSNEB

! default values for initialization :
REAL(KIND=JPRB) :: ZVALUE, ZVALUE_ONE, ZVALUE_T, ZVALUE_P, ZVALUE_L, ZVALUE_EPSILON

REAL(KIND=JPRB) :: ZVETAH(0:YDCPG_OPTS%KFLEVG)

!       Boolean
LOGICAL :: LLMSE, LLMSE_PARAM, LLMSE_DIAG
LOGICAL :: LLAROME
LOGICAL :: LLRAD
LOGICAL :: LLSWAP_THS, LLSWAP_RS, LLSWAP_SVS, LLSWAP_SVM, LLSWAP_LIMAS ! logical to swap or not pointers in and out
LOGICAL :: LLHN(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
LOGICAL :: LNUDGLHNREAD
LOGICAL :: LLIMAINIT

!       Characters
CHARACTER(LEN=11) :: CLNAME
CHARACTER(LEN=2),DIMENSION(7):: CLVARNAME=(/"QV","QL","QR","QI","QS","QG","QH"/)

! daand: radflex
REAL(KIND=JPRB), POINTER :: ZFRSO(:,:), ZFRTH(:,:)
TYPE(TYPE_INTPROC), POINTER :: YLRADPROC
REAL(KIND=JPRB)   :: ZCAPE(YDCPG_OPTS%KLON), ZDCAPE(YDCPG_OPTS%KLON)

!
! Phaser team note from CY43T1:
! there was a USE MODD_CTURB for accessing XTKEMIN here, but that created a forbidden
! dependence of APL_AROME (in "ifsarp") to the Méso-NH/Arome interfaces (in "mpa").
! There should be no USE MODD_* in APL_*.
! We decided to change the variable here to a local one, with the classical initial value for TKE.
!
REAL(KIND=JPRB), PARAMETER :: PPTKEMIN = 1.E-6


! Perturbed radiation-cloud interaction coef
REAL(KIND=JPRB), DIMENSION (YDCPG_OPTS%KLON) :: ZRADGR,ZRADSN

TYPE(ALL_SPP_VARS) :: ZSPP_ALL

!     ------------------------------------------------------------------
LOGICAL :: LLDIAB
LOGICAL :: LL_SAVE_PHSURF
LOGICAL :: LLXFUMSE

INTEGER(KIND=JPIM) :: IPTREXT,IEFB1,IEFB2,IEFB3
INTEGER(KIND=JPIM) :: IPTR(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
INTEGER(KIND=JPIM) :: IPTRLIMA
INTEGER(KIND=JPIM) :: IRR ! pointer of 1st hydrometeors in ZTENDGFLR
INTEGER(KIND=JPIM) :: IPTRTKE ! pointer of TKE in ZTENDGFLR

INTEGER(KIND=JPIM) :: IPGFL(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)

INTEGER(KIND=JPIM) :: INSTEP_DEB,INSTEP_FIN
INTEGER(KIND=JPIM) :: JROF, JSPP
INTEGER(KIND=JPIM) :: ISLB1U9  ,ISLB1V9  ,ISLB1T9  ,ISLB1GFL9, ISLB1VD9

!     --- UPPER AIR PHYSICAL TENDENCIES.
REAL(KIND=JPRB) :: ZTENDH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)     ! Enthalpy tendency.
REAL(KIND=JPRB) :: ZTENDPTKE(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)  ! Pseudo progn. TKE

! GFL tendencies for APL_AROME (assumes YDMODEL%YRML_GCONF%YGFL%NUMFLDS>=YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
! for now, use Jovi's trick :
REAL(KIND=JPRB), TARGET :: ZTENDGFLR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,0:YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
REAL(KIND=JPRB) :: ZTENDGFL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)   ! GFL tendencies

!     --- UPPER AIR PHYSICAL TENDENCIES FOR AROME.
!       (the previous one are not used in AROME)
REAL(KIND=JPRB) :: ZTENDT (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)        ! temperature tendency
REAL(KIND=JPRB) :: ZTENDW (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)        ! W  tendency
REAL(KIND=JPRB) :: ZTENDD (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)        ! d  tendency

REAL(KIND=JPRB) :: ZTENDU (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! U tendency without deep convection contribution
REAL(KIND=JPRB) :: ZTENDV (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! V tendency without deep convection contribution


!     ---FOR AROME PHYSICS  ---
REAL(KIND=JPRB) :: ZGWT1(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)     ! vertical velocity calculated by cputqy_arome before convertion in d
REAL(KIND=JPRB) :: ZTT1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)       ! Temperature at t1

! ZRTT1: appropriate version of R*T at t1 for gnhgw2svd
!  Version of R must be consistent with definition of vertical divergence.
REAL(KIND=JPRB) :: ZRTT1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)


REAL (KIND=JPRB) :: ZSAV_GZ0F (YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_UDOM (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZPRC_DPRECIPS2 (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%NDTPREC2)
REAL (KIND=JPRB) :: ZRDG_LCVQ (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZSAV_DDAL (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZPRC_DPRECIPS (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%NDTPREC)
REAL (KIND=JPRB) :: ZRDG_MU0M (YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZRDG_CVGQ (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZRDG_MU0N (YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_UDAL (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZSAV_UDGRO (YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_HV (YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_QSH  (YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_DDOM (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZSAV_UNEBH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZRDG_MU0 (YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_ENTCH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL (KIND=JPRB) :: ZRDG_MU0LU (YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_PBLH (YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_GZ0HF (YDCPG_OPTS%KLON)
REAL (KIND=JPRB) :: ZSAV_FHPS (YDCPG_OPTS%KLON)

! required for INTFLEX
TYPE(TYPE_INTPROCSET) :: YLPROCSET

! SPP
REAL(KIND=JPRB) :: ZGP2DSPP(YDCPG_OPTS%KLON,YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL)

REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTEND_Q    (:,:)  
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTEND_L    (:,:)  
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTEND_R    (:,:)  
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTEND_I    (:,:)  
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTEND_S    (:,:)  
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTEND_G    (:,:)  
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTEND_H    (:,:)  
REAL(KIND=JPRB), TARGET :: ZDUM2 (1,1)

REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDRA   (:,:,:) 

REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDLIMA (:,:,:) 
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDTKE  (:,:) 
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEFB1 (:,:) 
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEFB2 (:,:)
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEFB3 (:,:)
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEXT  (:,:,:)

REAL(KIND=JPRB), POINTER :: ZP1EZDIAG(:,:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "cpphinp.intfb.h"
#include "cptend_flex.intfb.h"
#include "cputqy_arome_expl.intfb.h"
#include "cputqy_arome_loop.intfb.h"
#include "cputqy.intfb.h"
#include "cp_ptrslb1.intfb.h"
#include "gnhgw2svdarome.intfb.h"
#include "writephysio.intfb.h"
#include "mf_phys_save_phsurf_part1.intfb.h"
#include "mf_phys_save_phsurf_part2.intfb.h"
#include "mf_phys_transfer.intfb.h"
#include "mf_phys_precips.intfb.h"
#include "apl_arome_calc_iptr.intfb.h"
#include "apl_arome_calc_ipgfl.intfb.h"
#include "abor1.intfb.h"
#include "recmwf.intfb.h"
#include "acraneb2.intfb.h"
#include "actqsat.intfb.h"
#include "acnpart.intfb.h"
#include "bri2acconv.intfb.h"
#include "gpgeo.intfb.h"
#include "gprcp_qlirsg.intfb.h"
#include "radheat.intfb.h"
#include "radghg.intfb.h"
#include "suozon.intfb.h"
#include "radaer.intfb.h"
#include "radact.intfb.h"
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
#include "aro_windfarm.intfb.h"
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


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('APL_AROME',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDGEM=>YDGEOMETRY%YRGEM, YDSTA=>YDGEOMETRY%YRSTA, YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,     &
& YLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH, YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI,                &
& YDRIP=>YDMODEL%YRML_GCONF%YRRIP, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,              &
& YDPHY3=>YDMODEL%YRML_PHY_MF%YRPHY3, YDPHYDS=>YDMODEL%YRML_PHY_MF%YRPHYDS, YDPHY0=>YDMODEL%YRML_PHY_MF%YRPHY0,            &
& YDVISI=>YDMODEL%YRML_PHY_MF%YRPHY%YRDVISI, YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, YGFL=>YDMODEL%YRML_GCONF%YGFL,            &
& YDMSE=>YDMODEL%YRML_PHY_MF%YRMSE, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR, YDPRECIPS=>YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS, &
& YDSTOPH=>YDMODEL%YRML_PHY_STOCH%YRSTOPH, YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1, YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH,        &
& YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY, YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YDEPHY=> YDMODEL%YRML_PHY_EC%YREPHY,                &
& YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, YDSPPT_CONFIG=>YDMODEL%YRML_GCONF%YRSPPT_CONFIG, YDSPPT=>YDMODEL%YRML_SPPT,             &
& YDSPP_CONFIG=>YDMODEL%YRML_GCONF%YRSPP_CONFIG,YDSPP=>YDMODEL%YRML_SPP,YDCVER=>YDGEOMETRY%YRCVER)

ASSOCIATE(MINPRR=>YDPARAR%MINPRR, MINPRS=>YDPARAR%MINPRS, MVQS=>YDPARAR%MVQS, MINPRG=>YDPARAR%MINPRG,                                 &
& LOTOWNC=>YDPARAR%LOTOWNC, LFPREC3D=>YDPARAR%LFPREC3D, NRRI=>YDPARAR%NRRI, NRRL=>YDPARAR%NRRL, &
& LTOTPREC=>YDPARAR%LTOTPREC, NPRINTFR=>YDPARAR%NPRINTFR, MALBDIR=>YDPARAR%MALBDIR,                     &
& NSWB_MNH=>YDPARAR%NSWB_MNH, XSW_BANDS=>YDPARAR%XSW_BANDS, MACPRG=>YDPARAR%MACPRG, MSWDIR=>YDPARAR%MSWDIR,                           &
& MSWDIF=>YDPARAR%MSWDIF, LOLSMC=>YDPARAR%LOLSMC, NDIAGWMAX=>YDPARAR%NDIAGWMAX,                               &
& MACPRS=>YDPARAR%MACPRS, MACPRR=>YDPARAR%MACPRR, LSQUALL=>YDPARAR%LSQUALL, &
& MALBSCA=>YDPARAR%MALBSCA, RADSN=>YDPARAR%RADSN, LDIAGWMAX=>YDPARAR%LDIAGWMAX,                                                       &
& NPTP=>YDPARAR%NPTP,                                          &
& NREFROI2=>YDPARAR%NREFROI2, NREFROI1=>YDPARAR%NREFROI1, MVEMIS=>YDPARAR%MVEMIS,                                                 &
& NRR=>YDPARAR%NRR, &
& RADGR=>YDPARAR%RADGR, XMINLM=>YDPHY0%XMINLM,                            &
& XMAXLM=>YDPHY0%XMAXLM, AERCS1=>YDPHY0%AERCS1, AERCS3=>YDPHY0%AERCS3, AERCS5=>YDPHY0%AERCS5, RDECRD1=>YDPHY0%RDECRD1,                &
& RDECRD2=>YDPHY0%RDECRD2, RDECRD3=>YDPHY0%RDECRD3, RDECRD4=>YDPHY0%RDECRD4, LMPA=>YDARPHY%LMPA, LUSECHEM=>YDARPHY%LUSECHEM,          &
& LKFBCONV=>YDARPHY%LKFBCONV, LMFSHAL=>YDARPHY%LMFSHAL, LMICRO=>YDARPHY%LMICRO, CCOUPLING=>YDARPHY%CCOUPLING,                         &
& LTURB=>YDARPHY%LTURB, LGRADHPHY=>YDARPHY%LGRADHPHY, LRDUST=>YDARPHY%LRDUST,                     &
& NGRADIENTS=>YDARPHY%NGRADIENTS, LRDEPOS=>YDARPHY%LRDEPOS,                               &
& LRCO2=>YDARPHY%LRCO2, LMSE=>YDARPHY%LMSE, LSURFEX_KFROM=>YDARPHY%LSURFEX_KFROM, NSURFEXCTL=>YDMSE%NSURFEXCTL,                       &
& XZSEPS=>YDMSE%XZSEPS, NDLUNG=>YDDIM%NDLUNG, NDGUNG=>YDDIM%NDGUNG, NPROMA=>YDDIM%NPROMA, NDLUXG=>YDDIM%NDLUXG,                       &
& NDGUXG=>YDDIM%NDGUXG, NGFL_EXT=>YGFL%NGFL_EXT, YLRAD=>YGFL%YLRAD, YIRAD=>YGFL%YIRAD, NGFL_EZDIAG=>YGFL%NGFL_EZDIAG,                 &
& NLIMA=>YGFL%NLIMA, CMICRO=>YDPARAR%CMICRO,                                                                                          &
& PHYEX=>YDPARAR%PHYEX,                                                                                                               &
& YSD_VAD=>YDSURF%YSD_VAD, QCO2=>YDPHY3%QCO2, NRAY=>YDPHY%NRAY,                            &
& LRAYFM=>YDPHY%LRAYFM, LO3ABC=>YDPHY%LO3ABC, LRAY=>YDPHY%LRAY, LRSTAER=>YDPHY%LRSTAER, LRNUEXP=>YDPHY%LRNUEXP,                       &
& AMAGSTOPH_CASBS=> YDSTOPH%AMAGSTOPH_CASBS, LFORCENL=>YDSTOPH%LFORCENL, NFORCESTART=>YDSTOPH%NFORCESTART,                            &
& NFORCEEND=>YDSTOPH%NFORCEEND, NTRADI=>YDTOPH%NTRADI, NTQSAT=>YDTOPH%NTQSAT, NTNEBU=>YDTOPH%NTNEBU, NAERMACC=>YDERAD%NAERMACC,       &
& NAER=>YDERAD%NAER, LHLRADUPD=>YDPHY%LHLRADUPD, TSPHY=>YDPHY2%TSPHY, NOZOCL=>YDERAD%NOZOCL, NRADFR=>YDERAD%NRADFR,                   &
& NSW=>YDERAD%NSW, RCARDI=>YDERDI%RCARDI, LFLEXDIA=>YLDDH%LFLEXDIA, LDDH_OMP=>YLDDH%LDDH_OMP, LRSLDDH=>YLDDH%LRSLDDH,                 &
& RDECLI=>YDRIP%RDECLI, RCODEC=>YDRIP%RCODEC, RHGMT=>YDRIP%RHGMT, RSIDEC=>YDRIP%RSIDEC, RSOVR=>YDRIP%RSOVR,                           &
& RSTATI=>YDRIP%RSTATI, TSTEP=>YDRIP%TSTEP, STPREH=>YDSTA%STPREH, LXXDIAGH=>YDXFU%LXXDIAGH, LFLASH =>YDCFU%LFLASH,                    &
& LDPRECIPS=>YDPHY%LDPRECIPS, LDPRECIPS2=>YDPHY%LDPRECIPS2, NDTPREC=>YDPRECIPS%NDTPREC, NDTPREC2=>YDPRECIPS%NDTPREC2,                 &
& NGPTOT=>YDGEM%NGPTOT, NGPBLKS=>YDDIM%NGPBLKS, NTSSG=>YDDPHY%NTSSG, YEZDIAG=>YGFL%YEZDIAG, YEXT=>YGFL%YEXT,                          &
& YNOGW=>YGFL%YNOGW, YCHEM=>YGFL%YCHEM, YSP_SBD=>YDSURF%YSP_SBD, LEDR=>YDPHY%LEDR, LAGPHY=>YDEPHY%LAGPHY,                             &
& YLIMA=>YGFL%YLIMA, LSPSDT => YDSPPT_CONFIG%LSPSDT, LKOGAN=>YDPARAR%LKOGAN,                                                          &
& LMODICEDEP=>YDPARAR%LMODICEDEP, LWINDFARM=>YDPHY%LWINDFARM,                                                                         &
 & RG=>YDCST%RG, RCPD=>YDCST%RCPD, RATM=>YDCST%RATM, RTT=>YDCST%RTT, RPI=>YDCST%RPI, &
 & RCW=>YDCST%RCW, RCPV=>YDCST%RCPV, RLVTT=>YDCST%RLVTT, RCS=>YDCST%RCS, RLSTT=>YDCST%RLSTT, &
 & RGAMW=>YDCST%RGAMW, RBETW=>YDCST%RBETW, RALPW=>YDCST%RALPW, RGAMS=>YDCST%RGAMS, &
 & RBETS=>YDCST%RBETS, RALPS=>YDCST%RALPS, RGAMD=>YDCST%RGAMD, RBETD=>YDCST%RBETD, &
 & RALPD=>YDCST%RALPD, RETV=>YDCST%RETV, RKAPPA=>YDCST%RKAPPA, RHOUR=>YDCST%RHOUR, RV=>YDCST%RV, RD=>YDCST%RD, &
& LTOTPRECL=>YDPARAR%LTOTPRECL, RSUN2=>YDMODEL%YRML_PHY_RAD%YRESWRT%RSUN2,&
& LSLAG=>YDDYNA%LSLAG, LGWADV=>YDDYNA%LGWADV, L_RDRY_VD=>YDDYNA%L_RDRY_VD,&
& LVERTFE=>YDCVER%LVERTFE, LVFE_GWMPA=>YDCVER%LVFE_GWMPA)

CALL SC2PRG(1, YEZDIAG(:)%MP, YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG, PGFL, ZP1EZDIAG)

!     ------------------------------------------------------------------

!        0.    constructor for procset
IF (LINTFLEX) YLPROCSET=NEWINTPROCSET()

!        1.    Preliminary calculations necessary
!              for all types of physics.
!              ------------------------------------

INSTEP_DEB=1
INSTEP_FIN=1

! initialisation for surfex if XFU
LLXFUMSE=.FALSE.
IF (YDCPG_OPTS%LCONFX) THEN
  LLXFUMSE=.TRUE.
ENDIF

! SPP 
IF ( YDSPP_CONFIG%LSPP ) THEN
 DO JSPP=1,YDMODEL%YRML_GCONF%YRSPP_CONFIG%SM%NRFTOTAL
   ZGP2DSPP(:,JSPP) = YDSPP%GP_ARP(JSPP)%GP2D(:,1,YDCPG_BNDS%KBL)
 ENDDO
ENDIF

! Complete physics is called.
LLDIAB=(.NOT.LAGPHY)

IF (LLDIAB) THEN
  CALL CPPHINP(YDGEOMETRY, YDMODEL, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDVARS%GEOMETRY%GEMU%T0,                                &
  & YDVARS%GEOMETRY%GELAM%T0, YDVARS%U%T0, YDVARS%V%T0, YDVARS%Q%T0, YDVARS%Q%DL, YDVARS%Q%DM, YDVARS%CVGQ%DL, YDVARS%CVGQ%DM,  &
  & YDCPG_PHY0%XYB%RDELP, YDCPG_DYN0%CTY%EVEL, YDVARS%CVGQ%T0, ZRDG_MU0, ZRDG_MU0LU, ZRDG_MU0M,                                 &
  & ZRDG_MU0N, ZRDG_CVGQ)
  ZRDG_LCVQ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=ZRDG_CVGQ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
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

LL_SAVE_PHSURF = .FALSE.

IF (LLDIAB) THEN
  LL_SAVE_PHSURF=YDCPG_OPTS%LCONFX
  IF (LL_SAVE_PHSURF) THEN
    CALL MF_PHYS_SAVE_PHSURF_PART1 (YDCPG_BNDS, YDCPG_OPTS, ZSAV_DDAL, ZSAV_DDOM, ZSAV_ENTCH,                           &
    & ZSAV_FHPS, ZSAV_GZ0F, ZSAV_GZ0HF, ZSAV_HV, ZSAV_PBLH, ZSAV_QSH, ZSAV_UDAL, ZSAV_UDGRO,                            &
    & ZSAV_UDOM, ZSAV_UNEBH, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VH%PPBLH, YDMF_PHYS_SURF%GSD_VH%PQSH,       &
    & YDMF_PHYS_SURF%GSD_VH%PSPSH, YDMF_PHYS_SURF%GSD_VK%PUDGRO, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VV%PZ0H, &
    & YDVARS%DAL%T0, YDVARS%DOM%T0, YDVARS%UAL%T0, YDVARS%UEN%T0, YDVARS%UNEBH%T0, YDVARS%UOM%T0,                       &
    & YDMODEL)
  ENDIF
ENDIF


CALL APL_AROME_CALC_IPGFL (YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDMODEL, IPGFL)

CALL MF_PHYS_TRANSFER (YDCPG_BNDS, YDCPG_OPTS, YDVARS, YDMODEL%YRML_PHY_MF%YRPHY, YDMODEL%YRML_GCONF%YGFL)

CALL APL_AROME_CALC_IPTR (YDMODEL, IEFB1, IEFB2, IEFB3, IPTR, IPTREXT, IPTRLIMA, IPTRTKE, IRR)

! If an incorrect address is used, then the initialization below will detect it :
ZTENDGFLR(:,:,0)=HUGE(1._JPRB)

ZTEND_Q    => ZTENDGFLR (:,:,IRR+0) ! 1 Q
ZTEND_L    => ZTENDGFLR (:,:,IRR+1) ! 2 L
ZTEND_R    => ZTENDGFLR (:,:,IRR+2) ! 3 R
ZTEND_I    => ZTENDGFLR (:,:,IRR+3) ! 4 I
ZTEND_S    => ZTENDGFLR (:,:,IRR+4) ! 5 S
ZTEND_G    => ZTENDGFLR (:,:,IRR+5) ! 6 G
IF (YDMODEL%YRML_PHY_MF%YRPARAR%NRR == 7) THEN
  ZTEND_H    => ZTENDGFLR (:,:,IRR+6) ! 7 H
ELSE
  ZTEND_H    => ZDUM2
ENDIF


ZTENDRA   => ZTENDGFLR (:, :, IRR:IRR+YDMODEL%YRML_PHY_MF%YRPARAR%NRR-1)
ZTENDLIMA => ZTENDGFLR (:, :, IPTRLIMA:IPTRLIMA+YDMODEL%YRML_GCONF%YGFL%NLIMA-1)
ZTENDTKE  => ZTENDGFLR (:, :, IPTRTKE)
ZTENDEFB1 => ZTENDGFLR (:, :, IEFB1)
ZTENDEFB2 => ZTENDGFLR (:, :, IEFB2)
ZTENDEFB3 => ZTENDGFLR (:, :, IEFB3)
ZTENDEXT  => ZTENDGFLR (:, :, IPTREXT:IPTREXT+YDMODEL%YRML_GCONF%YGFL%NGFL_EXT-1)






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

!         1.3 time step initialisation
!             the mesoNH physics (turb and microphysics) is written
!             for leap frog scheme
!             !!! be carefull for 2TL or 3TL

IF (YDDYNA%LTWOTL) THEN
  ZDT=YDCPG_OPTS%ZDTPHY/2._JPRB
ELSE
  IF (YDCPG_OPTS%NSTEP/=0) THEN
    ZDT=YDCPG_OPTS%ZDTPHY/2._JPRB
  ELSE
    ZDT=YDCPG_OPTS%ZDTPHY
  ENDIF
ENDIF

ZINVDT=1/YDCPG_OPTS%ZDTPHY

ZINVG=1._JPRB/RG 

!set concentration for LIMA
LLIMAINIT=.FALSE.
IF (YDCPG_OPTS%NSTEP==0 .AND. CMICRO=='LIMA') THEN
  LLIMAINIT=.TRUE.
  ZP1EZDIAG(:,:,1)=0._JPRB
  ZP1EZDIAG(:,:,2)=0._JPRB
  ZP1EZDIAG(:,:,3)=0._JPRB
  ZP1EZDIAG(:,:,4)=0._JPRB
  ZP1EZDIAG(:,:,5)=0._JPRB
ENDIF

! initialisation de ZDTMSE
IF (LLXFUMSE) THEN
  ZDTMSE=0.01_JPRB
  ZSTATI=REAL(RSTATI,JPRB)-ZDTMSE*.5_JPRB
  ZADTMS=0._JPRB
ELSE
  ZDTMSE=YDCPG_OPTS%ZDTPHY
  ZSTATI=REAL(RSTATI,JPRB)
  ZADTMS=ZDTMSE
ENDIF

IF(YDDYNA%LTWOTL) THEN
  ZRHGMT=REAL(RHGMT,JPRB)-ZDTMSE*.5_JPRB
ELSE
  ZRHGMT=REAL(RHGMT,JPRB)
ENDIF


LLMSE=LMSE.AND.(NSURFEXCTL >= 2)
LLMSE_PARAM=LLMSE
LLMSE_DIAG=LLMSE.AND.(NSURFEXCTL >= 3)


!  Vertical points
IKA=YDCPG_OPTS%KFLEVG
IKB=YDCPG_OPTS%KFLEVG
IKU=1
IKT=YDCPG_OPTS%KFLEVG
IKTE=YDCPG_OPTS%KFLEVG
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
  ZICEFR_(:,:)=ZVALUE
  ZPRCFR_(:,:)=ZVALUE
  ZICLDFR_(:,:)=ZVALUE
  ZWCLDFR_(:,:)=ZVALUE
  ZSSIO_(:,:)=ZVALUE
  ZSSIU_(:,:)=ZVALUE
  ZIFR_(:,:)=ZVALUE
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
  ZDTHRAD_(:,:)=ZVALUE

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
  YDCPG_MISC%RH(:,:)=ZVALUE
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

  YDMF_PHYS%OUT%FRSOC(:,:)=ZVALUE
  YDMF_PHYS%OUT%FRTHC(:,:)=ZVALUE
  YDMF_PHYS%OUT%FRSOPS(:)=ZVALUE
  YDMF_PHYS%OUT%DIAGH(:)=ZVALUE

  ZTENDSV_TURBLIMA_(:,:,:)=ZVALUE

ENDIF

!  INITIALIZE (CUMULATED) TENDENCIES

DO JLEV=1,YDCPG_OPTS%KFLEVG
  ZTENDT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.0_JPRB
  YDMF_PHYS%OUT%TENDU(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.0_JPRB
  YDMF_PHYS%OUT%TENDV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.0_JPRB
  ZTENDW(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.0_JPRB
  ZTENDTKE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.0_JPRB
ENDDO
DO JRR=1,NRR
  DO JLEV=1,YDCPG_OPTS%KFLEVG
    ZTENDRA(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,JRR)=0.0_JPRB
  ENDDO
ENDDO
DO JGFL=1,NGFL_EXT
  DO JLEV=1,YDCPG_OPTS%KFLEVG
    ZTENDEXT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,JGFL)=0.0_JPRB
  ENDDO
ENDDO
DO JGFL=1,NLIMA
  DO JLEV=1,YDCPG_OPTS%KFLEVG
    ZTENDLIMA(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,JGFL)=0.0_JPRB
  ENDDO
ENDDO

!  INITIALIZE CUMULATED STUFF

! Small array, OK. REK
ZINPRH_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
ZINPRR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
ZACPRR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
ZINPRS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
ZACPRS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
ZINPRG_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB


DO JLEV = 1,YDCPG_OPTS%KFLEVG
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZZZ_F_(JLON,JLEV)=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(JLON,JLEV)*ZINVG
    ZTENDTT(JLON,JLEV)=0._JPRB
  ENDDO
ENDDO

! adhoc solution to avoid negative tke values
DO JLEV=1, YDCPG_OPTS%KFLEVG
  DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
    ZTKEM(JLON,JLEV)=MAX(YDMF_PHYS_BASE_STATE%TKE(JLON,JLEV),PPTKEMIN)
  ENDDO
ENDDO

!initialisation of first useful field for EZDIAG use in Chemistry/Dust
IOFF_MFSHAL=1
IF(LFPREC3D) IOFF_MFSHAL=2

! 1.5 SPP settings
IF (YDSPP_CONFIG%LSPP) THEN
  CALL SET_ALL_SPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,NGFL_EZDIAG, &
   &  YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,YDSPP_CONFIG%SM%NRFTOTAL, &
   &  ZGP2DSPP,ZP1EZDIAG, &
   &  YDSPP_CONFIG,ZSPP_ALL)
ENDIF

!    ------------------------------------------------------------------
!     2 - CHANGEMENTS DE VARIABLES ET INVERSION DES NIVEAUX
!     --------------------------------------------------------------------

IF (LMICRO.OR.LTURB.OR.LLMSE.OR.LKFBCONV) THEN

  ZRSCP=RD/RCPD
  ZINVATM=1/RATM
  !initialisation de ZZZ_
  DO JLEV = 1,YDCPG_OPTS%KFLEVG
   !initialisation de qdm (utile localement pour calculer rho
   !et convertir q en r 
    IF (NRR==7) THEN
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZQDM(JLON,JLEV)=1._JPRB-YDMF_PHYS_BASE_STATE%Q(JLON,JLEV)-YDMF_PHYS_BASE_STATE%L(JLON,JLEV)-YDMF_PHYS_BASE_STATE%R(JLON,JLEV)&
         & -YDMF_PHYS_BASE_STATE%I(JLON,JLEV)-YDMF_PHYS_BASE_STATE%S(JLON,JLEV)-YDMF_PHYS_BASE_STATE%G(JLON,JLEV)-YDMF_PHYS_BASE_STATE%H(JLON,JLEV)   
      ENDDO
    ELSE
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZQDM(JLON,JLEV)=1._JPRB-YDMF_PHYS_BASE_STATE%Q(JLON,JLEV)-YDMF_PHYS_BASE_STATE%L(JLON,JLEV)-YDMF_PHYS_BASE_STATE%R(JLON,JLEV)&
         & -YDMF_PHYS_BASE_STATE%I(JLON,JLEV)-YDMF_PHYS_BASE_STATE%S(JLON,JLEV)-YDMF_PHYS_BASE_STATE%G(JLON,JLEV)
      ENDDO
    ENDIF
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
   !initialisation de ZRHODREFM__ (=qd*zrho)
      ZRHO=YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(JLON,JLEV)/(YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(JLON,JLEV)*YDMF_PHYS_BASE_STATE%T(JLON,JLEV))
      ZRHODREFM__(JLON,JLEV)=ZRHO*ZQDM(JLON,JLEV)
      ZRHODJM__(JLON,JLEV)=YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
  !initialisation de ZEXNREFM_
      ZEXNREFM_(JLON,JLEV)=(YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(JLON,JLEV)*ZINVATM)**(ZRSCP)
 ! vent horizontal et TKE
      ZPABSM__(JLON,JLEV)=YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(JLON,JLEV)
      ZUM__(JLON,JLEV)= YDMF_PHYS_BASE_STATE%U(JLON,JLEV)
      ZVM__(JLON,JLEV)= YDMF_PHYS_BASE_STATE%V(JLON,JLEV)
      ZWM__(JLON,JLEV)= YDMF_PHYS_BASE_STATE%YCPG_PHY%W(JLON,JLEV)
      ZTKEM__(JLON,JLEV)= ZTKEM(JLON,JLEV)
      ZZZ_(JLON,JLEV)=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO
  !initialise sigma for subgrid condensation coming
  !from previous time step turbulence scheme
  IF (PHYEX%NEBN%LSIGMAS) THEN
    DO JLEV = 1, YDCPG_OPTS%KFLEVG 
      ZSIGM_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)= YDMF_PHYS_BASE_STATE%SRC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
    ENDDO
  ENDIF
  !initialise convective mas flux for subgrid condensation coming
  !from previous time step convection scheme
  IF (PHYEX%NEBN%LSUBG_COND.AND..NOT.PHYEX%NEBN%LSIGMAS) THEN
    IF (LKFBCONV) THEN
      DO JLEV = 1, YDCPG_OPTS%KFLEVG
        ZMFM_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=YDMF_PHYS_BASE_STATE%SRC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
      ENDDO
    ELSE
      DO JLEV = 1, YDCPG_OPTS%KFLEVG
        ZMFM_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.
      ENDDO
    ENDIF
  ENDIF
!!! initialisation des variables d etat MNH

  !initialisation de ZRM_ pour les hydrometeores (ri=qi/qd)
  DO JLEV = 1, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZTHM__(JLON,JLEV)=YDMF_PHYS_BASE_STATE%T(JLON,JLEV)/ZEXNREFM_(JLON,JLEV)
      ZRM_(JLON,JLEV,1)=YDMF_PHYS_BASE_STATE%Q(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,2)=YDMF_PHYS_BASE_STATE%L(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,3)=YDMF_PHYS_BASE_STATE%R(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,4)=YDMF_PHYS_BASE_STATE%I(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,5)=YDMF_PHYS_BASE_STATE%S(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,6)=YDMF_PHYS_BASE_STATE%G(JLON,JLEV)/ZQDM(JLON,JLEV)  
    ENDDO
  ENDDO

  IF (NRR==7) THEN
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZRM_(JLON,JLEV,7)=YDMF_PHYS_BASE_STATE%H(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ENDDO
    ENDDO
  ENDIF

  IF (NRR==6) THEN
    !initialisation de ZTHVREFM__
    DO JLEV = 1, YDCPG_OPTS%KFLEVG 
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZTHVREFM__(JLON,JLEV)=ZTHM__(JLON,JLEV)*&
         & (1._JPRB+ZRM_(JLON,JLEV,1)*(RV/RD))/&
         & (1._JPRB+ZRM_(JLON,JLEV,1)+ZRM_(JLON,JLEV,2) +&
         & ZRM_(JLON,JLEV,3)+ZRM_(JLON,JLEV,4)+&
         & ZRM_(JLON,JLEV,5)+ZRM_(JLON,JLEV,6))
      ENDDO
    ENDDO
  ELSEIF (NRR==7) THEN
    DO JLEV = 1, YDCPG_OPTS%KFLEVG 
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
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
  DO JLEV = 1, YDCPG_OPTS%KFLEVG 
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZUS__(JLON,JLEV)= YDMF_PHYS_BASE_STATE%U(JLON,JLEV)*ZINVDT
      ZVS__(JLON,JLEV)= YDMF_PHYS_BASE_STATE%V(JLON,JLEV)*ZINVDT
      ZWS__(JLON,JLEV)= YDMF_PHYS_BASE_STATE%YCPG_PHY%W(JLON,JLEV)*ZINVDT
      ZTKES_(JLON,JLEV)= ZTKEM(JLON,JLEV)*ZINVDT
      ZTHS__(JLON,JLEV)=ZTHM__(JLON,JLEV)*ZINVDT
    ENDDO
  ENDDO

  !initialisation de ZRS_ pour les hydrometeores
  ! initialise pointers :
  CALL SWAP_RS
  DO JRR=1,NRR 
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZRS_(JLON,JLEV,JRR)=ZRM_(JLON,JLEV,JRR)*ZINVDT
      ENDDO
    ENDDO
  ENDDO

!!! Initialisations temporaires d'arguments non-utilises
  !initialisation de ZCIT_
  ZCIT_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:IKT)=0.0_JPRB
  
  !initialisation des tableaux de precipitations inst. and cumulated
  !and surface fluxes for turbulence
  IF (LLMSE.OR.LSFORCS) THEN
    ZACPRR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%ACPRR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
    ZACPRS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%ACPRS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
    ZACPRG_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%ACPRG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
    ZINPRR_NOTINCR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%INPRR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
    ZINPRS_NOTINCR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%INPRS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
    ZINPRG_NOTINCR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%INPRG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
  ENDIF

  !initialisation des scalaires passifs
  ! initialise pointers :
  CALL SWAP_SVM
  CALL SWAP_SVS
  DO JGFL=1,NGFL_EXT
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZSVM_(JLON,JLEV,JGFL)=YDMF_PHYS_BASE_STATE%P1EXT(JLON,JLEV,JGFL)
        ZSVS_(JLON,JLEV,JGFL)=YDMF_PHYS_BASE_STATE%P1EXT(JLON,JLEV,JGFL)*ZINVDT
      ENDDO
    ENDDO
  ENDDO

  !initialisation des concentrations LIMA
  ! initialise pointers :
  CALL SWAP_LIMAS
  DO JGFL=1,NLIMA
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZLIMAM_(JLON,JLEV,JGFL)=YDMF_PHYS_BASE_STATE%P1LIMA(JLON,JLEV,JGFL)
        ZLIMAS_(JLON,JLEV,JGFL)=YDMF_PHYS_BASE_STATE%P1LIMA(JLON,JLEV,JGFL)*ZINVDT
      ENDDO
    ENDDO
  ENDDO

  !initialisation de ZZI_THRAD
  IF (CMICRO=='LIMA') THEN
     IF (YDCPG_OPTS%NSTEP==0) THEN
        DO JLEV = 1, YDCPG_OPTS%KFLEVG 
           ZDTHRAD_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
        ENDDO
     ELSE
        DO JLEV = 1, YDCPG_OPTS%KFLEVG 
           ZDTHRAD_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,5)
        ENDDO
     ENDIF
  ENDIF

ENDIF

! daand: radflex
ZFRSO => YDMF_PHYS%OUT%FRSO(:,:,1)
ZFRTH => YDMF_PHYS%OUT%FRTH(:,:,1)

!    ------------------------------------------------------------------
!     3 - PRINTS FOR DIAGNOSTICS IF NEEDED
!    ------------------------------------------------------------------
IF (LDIAGWMAX) THEN
  IF (MOD(YDCPG_OPTS%NSTEP+1,NDIAGWMAX)==0) THEN
  ! calcul de wmax
    DO JLEV = 1 , YDCPG_OPTS%KFLEVG
      Z_WMAX=0._JPRB
      Z_WMIN=0._JPRB
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        IF (YDMF_PHYS_BASE_STATE%YCPG_PHY%W(JLON,JLEV)>Z_WMAX) THEN
          Z_WMAX=YDMF_PHYS_BASE_STATE%YCPG_PHY%W(JLON,JLEV)
        ENDIF
        IF (YDMF_PHYS_BASE_STATE%YCPG_PHY%W(JLON,JLEV)<Z_WMIN) THEN
          Z_WMIN=YDMF_PHYS_BASE_STATE%YCPG_PHY%W(JLON,JLEV)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
ENDIF

IF (LFLEXDIA) THEN
  !save tendencies
  ZTENDTBAK(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=ZTENDT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  DO JR=1,NRR
    ZTENDRBAK(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,JR)=ZTENDRA(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,JR)
  ENDDO
ENDIF


!    ------------------------------------------------------------------
!     4 - ADJUSTMENT (CALLED IF THE MICROPHYSICS IS SWITCH ON)
!    ------------------------------------------------------------------

IF (LMICRO) THEN

  ! Swap pointers because input values of THS and RS should be saved
  CALL SWAP_THS
  CALL SWAP_RS

  IF (LMFSHAL .AND. (PHYEX%PARAM_MFSHALLN%CMF_CLOUD=='DIRE'.OR.PHYEX%PARAM_MFSHALLN%CMF_CLOUD=='BIGA')) THEN
    IOFF_MFSHAL=IOFF_MFSHAL+3
    IF (YDCPG_OPTS%NSTEP==0) THEN
      DO JLEV = 1, YDCPG_OPTS%KFLEVG 
        ZRC_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
        ZRI_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
        ZCF_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
      ENDDO
    ELSE
      DO JLEV = 1, YDCPG_OPTS%KFLEVG 
        ZRC_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,1)
        ZRI_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,3)
        ZCF_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,2)
      ENDDO
    ENDIF
    ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,1:3)=0._JPRB
  ELSE
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      ZRC_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
      ZRI_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
      ZCF_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
    ENDDO
  ENDIF

  IF (MOD(YDCPG_OPTS%NSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'avant aro_adjust sous apl_arome'
    WRITE(NULOUT,*)'JLEV   ZZZ_F_   RHODJM   EXNREFM       PABSM       THM      SIGM         MFM    '
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,7F10.3)') JLEV,ZZZ_F_(NPTP,JLEV),ZRHODJM__(NPTP,JLEV),&
       & ZEXNREFM_(NPTP,JLEV), ZPABSM__(NPTP,JLEV), ZTHM__(NPTP,JLEV), ZSIGM_(NPTP,JLEV), ZMFM_(NPTP,JLEV)
    ENDDO 
    WRITE (NULOUT,*)'JLEV  rhoQv  rhoQc   rhoQr   rhoQi   rhoQs   rhoQg'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)') JLEV,ZRM_(NPTP,JLEV,1),&
       & ZRM_(NPTP,JLEV,2), ZRM_(NPTP,JLEV,3),ZRM_(NPTP,JLEV,4),ZRM_(NPTP,JLEV,5), ZRM_(NPTP,JLEV,6)
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZRC_MF_  ZRI_MF_  ZCF_MF_ ZTHSIN_  ZSRCS__ ZNEBMNH_'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)') JLEV,ZRC_MF_(NPTP,JLEV),&
       & ZRI_MF_(NPTP,JLEV),ZCF_MF_(NPTP,JLEV), ZTHSIN_(NPTP,JLEV),ZSRCS__(NPTP,JLEV), ZNEBMNH_(NPTP,JLEV)
    ENDDO
  ENDIF

  ! for now a copy is needed (see below, inside). I don't like than :-( REK
  ZTHS__(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=ZTHSIN_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ZRS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,1:NRR)=ZRSIN_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,1:NRR)

  IF (CMICRO == 'LIMA') THEN

    IF (LTURB) THEN
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        DO JLEV=1,YDCPG_OPTS%KFLEVG
          ZWNU_(JLON,JLEV) = ZWM__(JLON,JLEV) + 0.66*SQRT(ZTKEM__(JLON,JLEV))
        ENDDO
      ENDDO
      ZPTRWNU_ => ZWNU_(1:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
    ELSE
      ZPTRWNU_ => ZWM__(1:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
    ENDIF

    CALL SWAP_LIMAS
    ! for now a copy is needed (see below, inside). I don't like than :-( REK
    ZLIMAS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,1:NLIMA)=ZLIMASIN_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,1:NLIMA)

    CALL ARO_ADJUST_LIMA (PHYEX, &
    & YDCPG_OPTS%KFLEVG, IKU, IKL, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, YDCPG_BNDS%KFDIA, NRR,        &
    & NLIMA, YDCPG_OPTS%NSTEP+1, PHYEX%NEBN%LSUBG_COND, PHYEX%NEBN%LSIGMAS, ZDT, PHYEX%NEBN%VSIGQSAT, ZZZ_F_, ZRHODJM__(:, 1:YDCPG_OPTS%KFLEVG), &
    & ZRHODREFM__(:, 1:YDCPG_OPTS%KFLEVG), ZEXNREFM_, ZPABSM__(:, 1:YDCPG_OPTS%KFLEVG), ZTHM__(:, 1:YDCPG_OPTS%KFLEVG),   &
    & ZRM_, ZLIMAM_, ZSIGM_, ZPTRWNU_, ZDTHRAD_, ZMFM_, ZRC_MF_, ZRI_MF_, ZCF_MF_, ZTHS__(:, 1:YDCPG_OPTS%KFLEVG), ZRS_,                      &
    & ZLIMAS_, ZSRCS__(:, 1:YDCPG_OPTS%KFLEVG), ZNEBMNH_, ZICEFR_, ZPRCFR_, YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH,       &
    & LLIMAINIT                              )
  ELSE

    CALL ARO_ADJUST (PHYEX, &
    & YDCPG_BNDS%KFDIA, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, NRR, &
    & CMICRO, &
    & ZDT, ZZZ_F_, ZRHODJM__(:, 1:YDCPG_OPTS%KFLEVG),                                                                 &
    & ZEXNREFM_, ZRHODREFM__(:, 1:YDCPG_OPTS%KFLEVG), ZPABSM__(:, 1:YDCPG_OPTS%KFLEVG),                                         &
    & ZTHM__(:, 1:YDCPG_OPTS%KFLEVG), ZRM_, ZSIGM_,                                                                             &
    & ZMFM_, ZRC_MF_, ZRI_MF_, ZCF_MF_, ZTHS__(:, 1:YDCPG_OPTS%KFLEVG), ZRS_, ZSRCS__(:, 1:YDCPG_OPTS%KFLEVG),                  &
    & ZNEBMNH_, &
    & ZICLDFR_,ZWCLDFR_,ZSSIO_,ZSSIU_,ZIFR_,&
    & ZHLC_HRC_, ZHLC_HCF_, ZHLI_HRI_, ZHLI_HCF_,                                                                               &
    & YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH,&
    & ZSPP_ALL%YSPP_PSIGQSAT,ZSPP_ALL%YSPP_ICE_CLD_WGT)

  ENDIF
  
  IF (MOD(YDCPG_OPTS%NSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'apres aro_adjust sous apl_arome'
    WRITE(NULOUT,*)'JLEV   ZZZ_F_   RHODJM   EXNREFM       PABSM       THM      SIGM         MFM    '
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,7F10.3)') JLEV,ZZZ_F_(NPTP,JLEV),ZRHODJM__(NPTP,JLEV),&
       & ZEXNREFM_(NPTP,JLEV), ZPABSM__(NPTP,JLEV), ZTHM__(NPTP,JLEV), ZSIGM_(NPTP,JLEV), ZMFM_(NPTP,JLEV)
    ENDDO 
    WRITE (NULOUT,*)'JLEV  rhoQv  rhoQc   rhoQr   rhoQi   rhoQs   rhoQg'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)') JLEV,ZRS_(NPTP,JLEV,1),&
       & ZRS_(NPTP,JLEV,2), ZRS_(NPTP,JLEV,3),ZRS_(NPTP,JLEV,4),ZRS_(NPTP,JLEV,5), ZRS_(NPTP,JLEV,6)
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZRC_MF_  ZRI_MF_  ZCF_MF_ ZTHS__ ZSRCS__ ZNEBMNH_'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,ZRC_MF_(NPTP,JLEV),&
       & ZRI_MF_(NPTP,JLEV),ZCF_MF_(NPTP,JLEV), ZTHS__(NPTP,JLEV),ZSRCS__(NPTP,JLEV), ZNEBMNH_(NPTP,JLEV)
    ENDDO
  ENDIF

  DO JLEV=1,YDCPG_OPTS%KFLEVG
    YDVARS%A%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZNEBMNH_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV) 
  ENDDO

  !adjusted zthm and zrm
  DO JLEV = 1, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZTHM__(JLON,JLEV)=ZTHS__(JLON,JLEV)*YDCPG_OPTS%ZDTPHY
    ENDDO
  ENDDO

  DO JRR=1,NRR
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZRM_(JLON,JLEV,JRR)=ZRS_(JLON,JLEV,JRR)*YDCPG_OPTS%ZDTPHY
      ENDDO
    ENDDO
  ENDDO

  !initialisation de qdm utile pour
  !convertir tendance de r en tendance de q
  IF (NRR==6) THEN
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JLON= YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
        ZQDM(JLON,JLEV)=1._JPRB/(1._JPRB+ZRM_(JLON,JLEV,1)+&
        &ZRM_(JLON,JLEV,2)+ZRM_(JLON,JLEV,3)+ZRM_(JLON,JLEV,4)+ZRM_(JLON,JLEV,5)+&
        &ZRM_(JLON,JLEV,6) )
      ENDDO
    ENDDO
  ELSEIF (NRR==7) THEN
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JLON= YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
        ZQDM(JLON,JLEV)=1._JPRB/(1._JPRB+ZRM_(JLON,JLEV,1)+&
        &ZRM_(JLON,JLEV,2)+ZRM_(JLON,JLEV,3)+ZRM_(JLON,JLEV,4)+ZRM_(JLON,JLEV,5)+&
        &ZRM_(JLON,JLEV,6)+ZRM_(JLON,JLEV,7) )
      ENDDO
    ENDDO
  ENDIF 
  !reinitialisation des qi
  DO JLEV = 1, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZQVM(JLON,JLEV)=ZRM_(JLON,JLEV,1)*ZQDM(JLON,JLEV)
      ZQCM(JLON,JLEV)=ZRM_(JLON,JLEV,2)*ZQDM(JLON,JLEV)
      ZQRM(JLON,JLEV)=ZRM_(JLON,JLEV,3)*ZQDM(JLON,JLEV)
      ZQIM(JLON,JLEV)=ZRM_(JLON,JLEV,4)*ZQDM(JLON,JLEV)
      ZQSM(JLON,JLEV)=ZRM_(JLON,JLEV,5)*ZQDM(JLON,JLEV)
      ZQGM(JLON,JLEV)=ZRM_(JLON,JLEV,6)*ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  IF (NRR==7) THEN
    DO JLEV = 1,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZQHM(JLON,JLEV)=ZRM_(JLON,JLEV,7)*ZQDM(JLON,JLEV)
      ENDDO
    ENDDO
  ELSE
    ZQHM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=0._JPRB
  ENDIF

  ! Tendances des variables LIMA
  DO JGFL=1,NLIMA
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      ! Réinitialisation des variables LIMA
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZLIMAM_(JLON,JLEV,JGFL)=ZLIMAS_(JLON,JLEV,JGFL)*YDCPG_OPTS%ZDTPHY
        ZTENDLIMA(JLON,JLEV,JGFL)=ZTENDLIMA(JLON,JLEV,JGFL)+(ZLIMAS_(JLON,JLEV,JGFL)-ZLIMASIN_(JLON,JLEV,JGFL))  
      ENDDO
    ENDDO
  ENDDO

  !modif de R et CP
  ZQHGM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=ZQHM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)+ZQGM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
  CALL GPRCP_QLIRSG(YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, PQ=ZQVM, PQI=ZQIM, &
  & PQL=ZQCM, PQR=ZQRM, PQS=ZQSM, PQG=ZQHGM, PCP=ZCPM, PR=ZRHM)  

  DO JLEV = 1,YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZTM(JLON,JLEV)=ZTHM__(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
  !reinitialisation de ZRHODREFM__ (=qd*zrho)
      ZRHO=YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(JLON,JLEV)/(ZRHM(JLON,JLEV)*ZTM(JLON,JLEV))
      ZRHODREFM__(JLON,JLEV)=ZRHO*ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  !geopotentiel calculation
 
  ZAPHIM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
  CALL GPGEO(YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, ZAPHIM, ZAPHIFM,                &
  & ZTM, ZRHM, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDGEOMETRY%YRVERT_GEOM&
  &             )
   
  !calcul de l'altitude
  DO JLEV = 1, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZZZ_(JLON,JLEV)=ZAPHIM(JLON,JLEV)*ZINVG
      !initialisation de ZZZ_F_
      ZZZ_F_(JLON,JLEV)=ZAPHIFM(JLON,JLEV)*ZINVG
      ! tendency of T
      ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+(ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV))*ZEXNREFM_(JLON,JLEV)
      ZTENDTT(JLON,JLEV)=ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV)
    ENDDO
  ENDDO
  !inversion niveaux tendances des ri et conversion en qi en multipliant par qd
  DO JR=1,NRR
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZTENDRA(JLON,JLEV,JR)=ZTENDRA(JLON,JLEV,JR)+(ZRS_(JLON,JLEV,JR)-ZRSIN_(JLON,JLEV,JR))*ZQDM(JLON,JLEV)
      ENDDO
    ENDDO
  ENDDO
  !initialisation de ZDZZ_
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZDZZ_(JLON,1)=ZAPHIM(JLON,0)*ZINVG-ZZZ_(JLON,1)
  ENDDO
  DO JLEV = 2, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZDZZ_(JLON,JLEV)=ZZZ_(JLON,JLEV+IKL)-ZZZ_(JLON,JLEV)
    ENDDO
  ENDDO

ELSE

  ZTM (YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%T(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ZRHM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ZQVM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%Q(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ZQIM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%I(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ZQCM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%L(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ZQRM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%R(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ZQSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%S(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ZQGM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%G(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  IF (NRR==7) THEN
    ZQHM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%H(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ELSE
    ZQHM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=0._JPRB
  ENDIF
  ZCPM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KTDIA:YDCPG_OPTS%KFLEVG)
  ZAPHIM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,0:YDCPG_OPTS%KFLEVG)
  ZAPHIFM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ZZZ_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)*ZINVG
  !initialisation of PCLFS outside LMICRO to be zero in case LMICRO=F
  YDVARS%A%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=0._JPRB

ENDIF ! ADJUSTMENT LMICRO

!     ------------------------------------------------------------------
!     NEBULOSITE (CONVECTIVE+STRATIFORME) A TROIS NIVEAUX.
!     DIAGNOSTIC OF THREE LEVELS (CONVECTIVE+STRATIFORM) CLOUDINESS.

! protect cloudiness from being 0 or 1  (needed for ACRANEB2 and ACNPART)
DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZNEB0(JLON,JLEV)=MAX(ZEPSNEB,MIN(1._JPRB-ZEPSNEB,YDVARS%A%T1(JLON,JLEV)))
  ENDDO
ENDDO

! decorrelation depth for cloud overlaps

IF (LRNUEXP) THEN
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZDECRD(JLON)=RDECRD1+RDECRD2*EXP(-((ASIN(YDVARS%GEOMETRY%GEMU%T0(JLON))-RDECRD3*RDECLI)/RDECRD4)**2)
  ENDDO
ENDIF

! calculate high, medium, low and total cloud cover
CALL ACNPART(YDCST, YDMODEL%YRML_PHY_MF,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,YDCPG_OPTS%KLON,NTNEBU,YDCPG_OPTS%KFLEVG,&
 & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI,YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF,ZDECRD,ZNEB0,&
 & YDMF_PHYS%OUT%CLCH, YDMF_PHYS%OUT%CLCM, YDMF_PHYS%OUT%CLCL, YDCPG_MISC%CLCT, ZCLCT_RAD)

IF (LFLEXDIA) THEN
  DO JLEV = 1, YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA 
      ZTMPAF(JLON,JLEV)=(ZTENDT(JLON,JLEV)-ZTENDTBAK(JLON,JLEV))*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG*ZCPM(JLON,JLEV)
    ENDDO
  ENDDO
  CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,'TCTADJU',YDDDH)
  DO JR=1,NRR
    CLNAME='T'//CLVARNAME(JR)//'ADJU'
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA 
        ZTMPAF(JLON,JLEV)=(ZTENDRA(JLON,JLEV,JR)-ZTENDRA(JLON,JLEV,JR))*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
      ENDDO
    ENDDO
    CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,CLNAME,YDDDH)
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA 
        ZTMPAF(JLON,JLEV)=YDVARS%A%T1(JLON,JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)
      ENDDO
    ENDDO
    CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,'VNT',YDDDH)
  ENDDO
! specific to new data flow for diagnostics
  DO JLEV = 1, YDCPG_OPTS%KFLEVG
      ZCON1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV) = 1.0_JPRB
      ZCON2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV) = ZQDM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
  ENDDO
  DO JLEV = 1, YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA 
        ZCON3(JLON,JLEV) = YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
    ENDDO
  ENDDO
  ! missing interface !!! REK
    CALL ARO_SUINTBUDGET_OMP(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDCPG_OPTS%NSTEP, ZCON1, ZCON2, &
    & ZCON3, YDDDH)

ENDIF


DO JLEV = 1, YDCPG_OPTS%KFLEVG-1
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZDZZ_F_(JLON,JLEV)=ZZZ_F_(JLON,JLEV)-ZZZ_F_(JLON,JLEV-IKL)
  ENDDO
ENDDO
DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
  ZDZZ_F_(JLON,YDCPG_OPTS%KFLEVG)=ZZZ_F_(JLON,YDCPG_OPTS%KFLEVG)-YDVARS%GEOMETRY%OROG%T0(JLON)*ZINVG
ENDDO


!     --------------------------------------------------------------------
!     5 - COMPUTE DUST PROPERTIES FOR RADIATION IF LRDUST=T
!     --------------------------------------------------------------------
IF (LRDUST) THEN
  ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,IOFF_MFSHAL:NGFL_EZDIAG)=0.0_JPRB
  ! input dust scalar concentration in ppp from
  IEZDIAG_CHEM=NGFL_EZDIAG-IOFF_MFSHAL+1
  ! Swapp because IN and OUT will be needed simultaneously
  CALL SWAP_SVM
  ! input dust scalar concentration in ppp from
  CALL ARO_MNHDUST (IKL, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, NGFL_EXT, YDCPG_OPTS%ZDTPHY, ZSVMIN_, ZZZ_, ZDZZ_,                                       &
  & ZPABSM__(:, 1:YDCPG_OPTS%KFLEVG), ZTHM__(:, 1:YDCPG_OPTS%KFLEVG), ZRHODREFM__(:, 1:YDCPG_OPTS%KFLEVG),                                  &
  & NSWB_MNH, YDCPG_OPTS%NSTEP+1, ZSVM_, ZPIZA_DST_, ZCGA_DST_, ZTAUREL_DST_, ZAERD_, IEZDIAG_CHEM, ZPEZDIAG_(:, :, IOFF_MFSHAL:NGFL_EZDIAG)&
  &                                                                                                                                                                  )
  ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,IOFF_MFSHAL:NGFL_EZDIAG)=ZPEZDIAG_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,IOFF_MFSHAL:NGFL_EZDIAG)
! return to tendency
  DO JGFL=1, NGFL_EXT
    DO JLEV = 1,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZTENDEXT(JLON,JLEV,JGFL)=ZTENDEXT(JLON,JLEV,JGFL)+(ZSVM_(JLON,JLEV,JGFL)-ZSVMIN_(JLON,JLEV,JGFL))*ZINVDT
      ENDDO
    ENDDO
  ENDDO
ENDIF ! LRDUST

IF (LSFORCS) THEN   ! <== Surface forcing for MUSC

  ZTSURF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA) = YDMF_PHYS_BASE_STATE%T(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
  ZTN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)    = YDMF_PHYS_BASE_STATE%T(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
  ZQS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)    = YDMF_PHYS_BASE_STATE%Q(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZRHODREFM(JLON) = YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(JLON,YDCPG_OPTS%KFLEVG)/(YDMF_PHYS_BASE_STATE%T(JLON,YDCPG_OPTS%KFLEVG)*YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(JLON,YDCPG_OPTS%KFLEVG))
    ZTHETAS(JLON)   = ZTSURF(JLON)*(RATM/YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(JLON,YDCPG_OPTS%KFLEVG))**RKAPPA
  ENDDO

  LLAROME=.TRUE.
  CALL SURF_IDEAL_FLUX(YDRIP, YDPHY0, YDPHYDS, LLAROME, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,                             &
  & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(:, YDCPG_OPTS%KFLEVG), ZRHODREFM, YDMF_PHYS_SURF%GSD_SFO%PGROUP,                                  &
  & ZTN, ZTSURF, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_BASE_STATE%Q(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS_BASE_STATE%U(:, YDCPG_OPTS%KFLEVG), &
  & YDMF_PHYS_BASE_STATE%V(:, YDCPG_OPTS%KFLEVG), ZTHETAS, ZSFTH_, ZSFRV_, ZSFU_, ZSFV_)

!* Compute PBL-diagnostics
   
   ZCAPE(:)=0._JPRB
   ZDCAPE(:)=0._JPRB
  CALL ACCLDIA(YDCST, YDCPG_OPTS%LXCLP, YDCPG_OPTS%LXTGST, YDCPG_OPTS%LXXGST, YDPHY, YDMODEL%YRML_PHY_MF%YRPHY2, YDTOPH, YDCPG_BNDS%KIDIA,    &
  & YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS_BASE_STATE%U(:, 1:YDCPG_OPTS%KFLEVG),     &
  & YDMF_PHYS_BASE_STATE%V(:, 1:YDCPG_OPTS%KFLEVG), ZCAPE, ZDCAPE, ZTKEM(:, 1:YDCPG_OPTS%KFLEVG), YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(:, 1:YDCPG_OPTS%KFLEVG), &
  & YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST, YDMF_PHYS%OUT%CLPH, ICLPH)

  YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=MIN(XMAXLM,MAX(XMINLM,YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA))) 

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
  CALL ACTQSAT (YDCST, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTQSAT, YDCPG_OPTS%KFLEVG,   &
                & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF, ZCPM, ZQVM, ZTM, ZGEOSLC, ZLH, ZLSCPE, ZQSAT, ZQW, YDCPG_MISC%RH, &
                & ZTW)  

  IF (ZSPP_ALL%YSPP_RADGR%LPERT) THEN
   CALL APPLY_SPP(ZSPP_ALL%YSPP_RADGR, &
                & YDCPG_OPTS%KLON,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA, &
                & RADGR,ZRADGR)
  ELSE
   DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZRADGR(JLON) = RADGR
   ENDDO
  ENDIF

  IF (ZSPP_ALL%YSPP_RADSN%LPERT) THEN
   CALL APPLY_SPP(ZSPP_ALL%YSPP_RADSN, &
                & YDCPG_OPTS%KLON,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA, &
                & RADSN,ZRADSN)
  ELSE
   DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZRADSN(JLON) = RADSN
   ENDDO
  ENDIF

  ! initialisation des humidite (dans le rayonnement, l'eau liquide nuageuse
  ! et la glace sont donne par des hu par rapport au gaz.
  ! (qi/qa+qv pour ice par ex. C'est donc different de ri)
  DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA 
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
  IF (YIRAD%LGP) YDVARS%IRAD%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = ZQICE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
  IF (YLRAD%LGP) YDVARS%LRAD%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = ZQLIQ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)

  ! Hannu Savijarvi diffuse -> direct albedo correction from hlradia,
  ! Assuming that SURFEX does not make difference between
  ! dir/dif albedo as surfex/SURFEX/albedo_from_nir_vis.F90 defines
  ! PSCA_ALB(:,:) = PDIR_ALB(:,:)

! Albedo dans les intervalles, direct (parallel) et diffus (diffuse).
  IF (NSW==6.OR.NSW==1) THEN
    IF (LLMSE) THEN
      DO JSW=1,NSW
        ZALBP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=YDCPG_GPAR%ALBDIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)
        ZALBD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=YDCPG_GPAR%ALBSCA(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)
        IF (LHLRADUPD) THEN
          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            ZSALBCOR=0.2_JPRB/(1._JPRB+ZRDG_MU0(JLON))-0.12_JPRB
            ZALBP(JLON,JSW)=ZALBD(JLON,JSW)+ZSALBCOR
          ENDDO
        ENDIF
      ENDDO
    ELSEIF (LSFORCS) THEN
      DO JSW=1,NSW
        ZALBP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=RALB_FORC
        ZALBD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=RALB_FORC
!  direct>diffuse correction might be applied to RALB_FORC,too:
!              ZALBP(JLON,JSW)=RALB_FORC+ZSALBCOR
      ENDDO
    ELSE
     !pour pouvoir tourner sans la surface
      DO JSW=1,NSW
        ZALBP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=YDMF_PHYS_SURF%GSD_VF%PALBF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
        ZALBD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=YDMF_PHYS_SURF%GSD_VF%PALBF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
!              ZALBP(JLON,JSW)=PALBIN(JLON)+ZSALBCOR
      ENDDO
    ENDIF

  ! Spectral average albedo done with RSUN2 weights,
  ! to be applied for HLRADIA, ACRANEB2 which use a single solar spectral band
    IF (LHLRADUPD) THEN
      ZALBP1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
      ZALBD1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
      DO JSW=1,NSW
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZALBP1(JLON)=ZALBP1(JLON)+RSUN2(JSW)*ZALBP(JLON,JSW)
          ZALBD1(JLON)=ZALBD1(JLON)+RSUN2(JSW)*ZALBD(JLON,JSW)
        ENDDO
      ENDDO
    ELSE
       ZALBP1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%ALBDIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)
       ZALBD1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%ALBSCA(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)
    ENDIF
  ELSE
    CALL ABOR1 ('ALBEDO FOR NSW/= 1 or 6 not defined in apl_arome')
  ENDIF

  ! all albedo operations

  IF (LLMSE) THEN
    ZEMIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%VEMIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
    ZTSURF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%VTS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
    ! protection for E Zone, Where surface scheme send back EMIS and T =0
    ! the protection in aro_ground_paramn is not sufficient !!! WHY ??
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      IF (ZEMIS(JLON)==0._JPRB) THEN
        ZEMIS(JLON)=1.0_JPRB
        ZTSURF(JLON)=288.0_JPRB
      ENDIF
    ENDDO
  ELSEIF (LSFORCS) THEN
    ZEMIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=REMIS_FORC
  ELSE
    ZEMIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.5_JPRB ! value 0.5 is suspicious
    ZTSURF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZTM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
  ENDIF !LLMSE EMIS
  
  ! initialization of aerosols
  IF (LRAY.AND.NRAY == 2.AND.LRADFLEX.AND.LRSTAER) THEN

    ! old ("standard") aerosols for LRAY only
    ZAER(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KTDIA-1,1)=0._JPRB
    DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
      ZVETAH(JLEV)=STPREH(JLEV)/VP00
    ENDDO
    ZAEO=AERCS1*ZVETAH(YDCPG_OPTS%KTDIA-1)+AERCS3*ZVETAH(YDCPG_OPTS%KTDIA-1)**3+AERCS5*ZVETAH(YDCPG_OPTS%KTDIA-1)**5
    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      ZAEN=AERCS1*ZVETAH(JLEV)+AERCS3*ZVETAH(JLEV)**3+AERCS5*ZVETAH(JLEV)**5
      ZAER(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,1)=ZAEN-ZAEO
      ZAEO=ZAEN
    ENDDO
    ZAER(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,2:6)=0._JPRB
  
  ELSE
    
    IF (NAER >= 1 .AND. NAERMACC == 0) THEN
      IF(YSD_VAD%NUMFLDS >= 4) THEN
        CALL RADAER ( YDMODEL%YRML_PHY_RAD%YREAERD, YDERAD, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,               &
        & YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF, &
        & ZTM, ZTSURF, YDMF_PHYS_SURF%GSD_VA%PSEA, YDMF_PHYS_SURF%GSD_VA%PLAN, YDMF_PHYS_SURF%GSD_VA%PSOO,           &
        & YDMF_PHYS_SURF%GSD_VA%PDES, YDMF_PHYS_SURF%GSD_VA%PSUL, YDMF_PHYS_SURF%GSD_VA%PVOL, ZAER,                  &
        & ZAERINDS)
      ELSE
        WRITE(NULOUT,*) 'YSD_VAD%NUMFLDS SHOULD BE >= 4, IT IS: ',YSD_VAD%NUMFLDS
        CALL ABOR1('APL_AROME: PB AEROSOLS!')
        ! NB : this abort excludes the use of radact. REK.
      ENDIF
    ELSE
      !AEROSOLS from MACC (NAERMACC=1)
      ZDUM=1._JPRB
      ! in E Zone, there are YDVARS%GEOMETRY%GEMU%T0 < 0.
      ZGELAM=YDVARS%GEOMETRY%GELAM%T0
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        IF (ZGELAM(JLON)<0._JPRB) THEN
              ZGELAM(JLON)=ZGELAM(JLON)+2*RPI
        ENDIF
      ENDDO
      ! Init ZCHTIX
      ! Warning YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE is 0:YDCPG_OPTS%KFLEVG here whereas in radintg it is
      ! 1:YDCPG_OPTS%KFLEVG+1
      DO JK=2,YDCPG_OPTS%KFLEVG
        ZCAPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK)=YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK-1)
        ZTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK)=(ZTM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK-1)*YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK-1)&
      & *(YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK)-YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK-1))&
      & +ZTM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK)*YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK)*(YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK-1)&
     & -YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK-1)))&
     &*(1.0_JPRB/(YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK-1)*(YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK)-YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK-1))))
        ZCHTIX(1:YDCPG_OPTS%KLON,JK)=ZTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JK)
      ENDDO
      ! QUANTITIES AT BOUNDARIES
      ZCAPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)=YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,0)
      ZCAPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG+1)=YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)
      ZTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)=ZTM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)-YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)*(ZTM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)&
       & -ZTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,2))/(YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)-YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1))
      ZTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG+1)=ZTSURF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      ZCHTIX(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG+1)=ZTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG+1)
      ZCHTIX(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)=ZTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)
      IAERO=SIZE(ZAERO,3)
      ZGEMU_D=REAL(YDVARS%GEOMETRY%GEMU%T0,JPRD)
      CALL RADACT(YDMODEL%YRML_PHY_RAD%YREAERD,YDERAD,YDMODEL%YRML_PHY_AER%YREAERSNK,YDRIP, YDSPP_CONFIG,&
         & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,&
         & 1 , YDCPG_OPTS%KLON , YDCPG_OPTS%KLON  ,0 , 1,&
         & YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE , ZGELAM, ZGEMU_D, YDVARS%GEOMETRY%GECLO%T0, YDVARS%GEOMETRY%GESLO%T0, ZCHTIX,&
         & ZQVM   , ZQSAT  , ZDUM  ,&
         & ZRAER  , ZAERO,ZROZ  )
      DO JAE=1,6
        DO JK=1,YDCPG_OPTS%KFLEVG
          DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            ZAER(JLON,JK,JAE)=ZRAER(JLON,JAE,JK)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (LRDUST) THEN
      ! We use the extinction coefficient explicitly solved by ARO_MNHDUST
      ZAER(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,3) = ZAERD_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
    ENDIF

  ENDIF
  ! end of old or new aerosols

  ! initialization of ozone
  IF (NOZOCL == 1) THEN
    ! as in IFS
    CALL RADOZC(YDRIP%YREOZOC,YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, 1, YDCPG_OPTS%KLON, &
    & 0, YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, YDVARS%GEOMETRY%GEMU%T0, ZROZ)
    DO JK=1,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZQO3(JLON,JK)=ZROZ(JLON,JK)/YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JK)
      ENDDO
    ENDDO
  ELSEIF (NOZOCL == 2) THEN
    ! as in ARPEGE (from clim profiles)
    CALL SUOZON(YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, ZQO3,                                  &
    & .FALSE., YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, LO3ABC, YDMF_PHYS_SURF%GSD_VC%PA,   &
    & YDMF_PHYS_SURF%GSD_VC%PB, YDMF_PHYS_SURF%GSD_VC%PC)
  ENDIF
  IF (NOZOCL==3.OR.NOZOCL==4) THEN ! Clims MACC
     CALL RADGHG (YDERAD,YDRIP,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KLON, &
     & YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE,YDVARS%GEOMETRY%GEMU%T0,&
     & ZQCO2, ZQCH4, ZQN2O, ZQNO2, ZQC11, ZQC12, ZROZ, ZQC22, ZQCL4 )
    DO JK=1,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZQO3(JLON,JK)=ZROZ(JLON,JK)
      ENDDO
    ENDDO
  ENDIF


ELSE

  DO JSW=1,NSW
    ZALBP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=0._JPRB 
    ZALBD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=0._JPRB
  ENDDO

ENDIF
 !of preparation of input for LRAYFM, LRAY at every time step

 IF (LRAYFM) THEN
   ! Intermittent call to radiation interface
   IF (MOD(YDCPG_OPTS%NSTEP,NRADFR) == 0) THEN 
     CALL RECMWF (YDGEOMETRY%YRDIMV, YDMODEL, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, &
     & YDCPG_OPTS%KSW, &
     & NOZOCL   ,NAERMACC, IAERO, &
     & ZALBD, ZALBP, YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF,           &
     & YDVARS%A%T1, ZQO3,ZQCO2 , ZQCH4    , ZQN2O   , &
     & ZQNO2   , ZQC11   , ZQC12 , ZQC22    , ZQCL4   , &
     & ZAER, ZAERO, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZEMIS, &
     & ZRDG_MU0M, ZQV, ZQSAT,                 &
     & ZQICE, ZQLIQ, ZQSM, ZQRM, YDMF_PHYS_SURF%GSD_VF%PLSM, ZTM, ZTSURF, YDMF_PHYS%RAD%EMTD, YDMF_PHYS%RAD%EMTU,     &
     & YDMF_PHYS%RAD%TRSW, YDMF_PHYS%OUT%FRTHC, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRSO,          &
     & ZZS_FSWDIR, ZZS_FSWDIF, ZFSDNN, ZFSDNV, ZCTRSO, ZCEMTR, ZTRSOD, ZTRSODIR, ZTRSODIF,                            &
     & ZPIZA_DST_, ZCGA_DST_, ZTAUREL_DST_, ZAERINDS, YDVARS%GEOMETRY%GELAM%T0, YDVARS%GEOMETRY%GEMU%T0,              &
     & YDCPG_GPAR%SWDIR, YDCPG_GPAR%SWDIF, ZRDG_MU0LU, ZALBD1, ZFRSOLU, &
     & YSPP_RSWINHF=ZSPP_ALL%YSPP_RSWINHF, YSPP_RLWINHF=ZSPP_ALL%YSPP_RLWINHF)
   ELSE
     IF (LLMSE) THEN
       DO JSW=1,NSW
         ZTRSODIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=YDCPG_GPAR%SWDIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)
         ZTRSODIF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=YDCPG_GPAR%SWDIF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)
       ENDDO
     ENDIF
     ZCTRSO(:,:)=0._JPRB
   ENDIF
   ! daand: radflex
   IF (LRADFLEX) THEN
     YLRADPROC => NEWINTPROC(YLPROCSET,'Radiation')
     ZFRSO => NEWINTFIELD(YLRADPROC,YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,'FRSO','H','F')
     ZFRTH => NEWINTFIELD(YLRADPROC,YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,'FRTH','H','F')
   ENDIF

    DO JLEV=1,YDCPG_OPTS%KFLEVG
      ZTENT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.0_JPRB
    ENDDO

   ZSUDU(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.0_JPRB

   CALL RADHEAT  (YDMODEL%YRCST,  YDMODEL%YRML_PHY_EC%YRTHF, YDERAD, YDERDI, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,     &
   & YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, ZEMIS, YDMF_PHYS%RAD%EMTD, ZRDG_MU0, ZQVM,             &
   & ZTENT, YDMF_PHYS%RAD%TRSW, ZTRSOD, ZTSURF, YDCPG_OPTS%ZDTPHY, ZTRSODIR, ZTRSODIF, ZALBD, ZALBP, ZFRSO,                  &
   & ZFRTH, YDMF_PHYS%OUT%FRSODS, YDMF_PHYS%OUT%FRTHDS, ZCEMTR, ZCTRSO, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRTHC, &
   & ZSUDU, ZSDUR, ZDSRP, ZZS_FSWDIR, ZZS_FSWDIF, YDMF_PHYS%OUT%FRSOPS, ZFRSOFS, YDMF_PHYS%OUT%FRSOPT             &
   &                                                                                                                        )

  ! daand: radflex
  IF (LRADFLEX) THEN
    ! store for further calculations and diagnostics
    ! warning : pointers. REK
    YDMF_PHYS%OUT%FRSO(:,:,1)=ZFRSO
    YDMF_PHYS%OUT%FRTH(:,:,1)=ZFRTH
  ELSE
  ! daand: if LRADFLEX, the contribution to temperature is done by
  ! cptend_flex/cputqy
    ! update temperature tendency by radiative contribution
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
        ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+ZTENT(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF

  !initialisation de ZZI_THRAD
  IF (CMICRO=='LIMA') THEN
     DO JLEV = 1, YDCPG_OPTS%KFLEVG
        DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
           ZDTHRAD_(JLON,JLEV)=ZTENT(JLON,JLEV)/ZEXNREFM_(JLON,JLEV)
        END DO
        ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,5)=ZDTHRAD_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
     ENDDO
  ENDIF
     
  DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
    ! update sunshine duration [s]
    !YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)=YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)+ZSDUR(JLON)*TSTEP
    YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)=YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)+ZSDUR(JLON)*ZADTMS ! fix stepx case
    ! Estimate of the direct normal irradiance, with securities
    IF (ZRDG_MU0(JLON) > 3.0E-02_JPRB) THEN
      YDMF_PHYS%OUT%FRSDNI(JLON)=MAX(0.0_JPRB,YDMF_PHYS%OUT%FRSOPS(JLON)/ZRDG_MU0(JLON))
    ELSE
      YDMF_PHYS%OUT%FRSDNI(JLON)=MAX(0.0_JPRB,YDMF_PHYS%OUT%FRSOPS(JLON))
    ENDIF
  ENDDO 

  IF( MOD(YDCPG_OPTS%NSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'sous apl_arome apres rayonnement ZTENT=',ZTENT(NPTP,30:41)
    IF (LLMSE) THEN
      DO JSW=1, NSW
        WRITE(NULOUT,*)'ZSFSWDIR ZSFSWDIF ZFSDNN ZFSDNV PFRSO',&
         & ZZS_FSWDIR(NPTP,JSW),ZZS_FSWDIF(NPTP,JSW),ZFSDNN(NPTP), ZFSDNV(NPTP),YDMF_PHYS%OUT%FRSO(NPTP,YDCPG_OPTS%KFLEVG,1)
        WRITE(NULOUT,*)'ZALBD ZALBP',ZALBD(NPTP,JSW),ZALBP(NPTP,JSW)
      ENDDO
    ENDIF
    WRITE(NULOUT,*)ZFSDNN(NPTP),ZFSDNV(NPTP)
    WRITE (NULOUT,*)'TSURF EMIS ZFRTH',ZTSURF(NPTP),ZEMIS(NPTP),YDMF_PHYS%OUT%FRTHDS(NPTP)
  ENDIF

  IF (LFLEXDIA) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, YDMF_PHYS%OUT%FRSO(:, :, 1), 'FCTRAYSO', YDDDH&
      &             )
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, YDMF_PHYS%OUT%FRTH(:, :, 1), 'FCTRAYTH', YDDDH&
      &             )
  ENDIF

ELSE

  YDMF_PHYS%OUT%FRSOC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,0:1)=0.0_JPRB
  YDMF_PHYS%OUT%FRTHC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,0:1)=0.0_JPRB

ENDIF  ! LRAYFM


IF (LRAY.AND.NRAY == 2.AND.LRADFLEX) THEN

  ! -------------------------
  ! ACRANEB2 radiation scheme
  ! -------------------------

!+++ The next input preparations are redundant:

  ! initialization of cloud ice, cloud liquid and specific humidity
  ! (with respect to moist air, i.e. excluding hydrometeors)
  DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
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
  IF (YIRAD%LGP) YDVARS%IRAD%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = ZQICE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
  IF (YLRAD%LGP) YDVARS%LRAD%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = ZQLIQ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)

  ! initialization of ozone
  IF (NOZOCL == 1) THEN
    ! as in IFS
    CALL RADOZC(YDRIP%YREOZOC,YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, 1, YDCPG_OPTS%KLON, &
    & 0, YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, YDVARS%GEOMETRY%GEMU%T0, ZROZ)
    DO JK=1,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZQO3(JLON,JK)=ZROZ(JLON,JK)/YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JK)
      ENDDO
    ENDDO
  ELSEIF (NOZOCL == 2) THEN
    ! as in ARPEGE (from clim profiles)
    CALL SUOZON(YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, ZQO3,                                  &
    & .FALSE., YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, LO3ABC, YDMF_PHYS_SURF%GSD_VC%PA,   &
    & YDMF_PHYS_SURF%GSD_VC%PB, YDMF_PHYS_SURF%GSD_VC%PC)
  ENDIF
  IF (NOZOCL==3.OR.NOZOCL==4) THEN ! Clims MACC
     CALL RADGHG (YDERAD,YDRIP,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KLON, &
     & YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE,YDVARS%GEOMETRY%GEMU%T0,&
     & ZQCO2, ZQCH4, ZQN2O, ZQNO2, ZQC11, ZQC12, ZROZ, ZQC22, ZQCL4 )
    DO JK=1,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZQO3(JLON,JK)=ZROZ(JLON,JK)
      ENDDO
    ENDDO
  ENDIF


  ! initialization of aerosols
  IF (LRAY.AND.NRAY == 2.AND.LRADFLEX.AND.LRSTAER) THEN

    ZAER(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KTDIA-1,1)=0._JPRB
    ! old ("standard") aerosols
    DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
      ZVETAH(JLEV)=STPREH(JLEV)/VP00
    ENDDO
    ZAEO=AERCS1*ZVETAH(YDCPG_OPTS%KTDIA-1)+AERCS3*ZVETAH(YDCPG_OPTS%KTDIA-1)**3+AERCS5*ZVETAH(YDCPG_OPTS%KTDIA-1)**5
    DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
      ZAEN=AERCS1*ZVETAH(JLEV)+AERCS3*ZVETAH(JLEV)**3+AERCS5*ZVETAH(JLEV)**5
      ZAER(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,1)=ZAEN-ZAEO
      ZAEO=ZAEN
    ENDDO
    ZAER(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,2:6)=0._JPRB
  
  ELSE

    IF (NAER >= 1) THEN
      IF (YSD_VAD%NUMFLDS >= 4) THEN
        ! initialisation of aerosols as in ARPEGE (from clim files)
        CALL RADAER (YDMODEL%YRML_PHY_RAD%YREAERD, YDERAD, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,                &
        & YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF, &
        & ZTM, ZTSURF, YDMF_PHYS_SURF%GSD_VA%PSEA, YDMF_PHYS_SURF%GSD_VA%PLAN, YDMF_PHYS_SURF%GSD_VA%PSOO,           &
        & YDMF_PHYS_SURF%GSD_VA%PDES, YDMF_PHYS_SURF%GSD_VA%PSUL, YDMF_PHYS_SURF%GSD_VA%PVOL, ZAER,                  &
        & ZAERINDS)
      ELSE
        CALL ABOR1('APL_AROME: PB AEROSOLS!')
        ! NB : this abort excludes the use of radact. REK.
      ENDIF
    ENDIF

    IF (LRDUST) THEN
      ! We uses the extinction coefficient explicitely solved by ARO_MNHDUST
      ZAER(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,3) = ZAERD_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
    ENDIF

  ENDIF ! (LRAY.AND.NRAY == 2.AND.LRADFLEX.AND.LRSTAER) THEN

  ! get diffuse and direct surface albedo, emissivity and temperature
  IF (.NOT.LHLRADUPD) THEN
    ZALBD1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%ALBSCA(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)
    ZALBP1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%ALBDIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)
  ENDIF
  ZEMIS  (YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%VEMIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
  ZTSURF (YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDCPG_GPAR%VTS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ! protection of E-zone (not to have zero emissivity and T_surf there)
    IF (ZEMIS(JLON) == 0._JPRB) THEN
      ZEMIS (JLON)=  1._JPRB
      ZTSURF(JLON)=288._JPRB
    ENDIF
  ENDDO

!+++ End of redundant input preparations for ACRANEB

  ! initialization of CO2(+), differs from IFS radiation scheme!
  ZQCO2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=QCO2

  ! daand: radflex
  YLRADPROC => NEWINTPROC(YLPROCSET,'Radiation')
  ZFRSO => NEWINTFIELD(YLRADPROC,YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG, 'FRSO','H','F')
  ZFRTH => NEWINTFIELD(YLRADPROC,YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG, 'FRTH','H','F')

  ! call radiation scheme
  IJN=YDCPG_OPTS%KLON
  CALL ACRANEB2(YDERDI, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,               &
  & NTRADI, YDCPG_OPTS%KFLEVG, IJN, YDCPG_OPTS%NSTEP, YDCFU%NFRRC, YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE,                  &
  & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,     &
  & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZNEB0, ZQV, ZQCO2, ZQICE, ZQLIQ, ZQO3, YDMF_PHYS_BASE_STATE%T,             &
  & ZALBD1, ZALBP1, ZEMIS, YDVARS%GEOMETRY%GELAM%T0, YDVARS%GEOMETRY%GEMU%T0, ZRDG_MU0, ZRDG_MU0LU,                    &
  & ZTSURF, ZDECRD, ZCLCT_RAD, YDMF_PHYS%OPT%GDEOSI, YDMF_PHYS%OPT%GUEOSI, YDMF_PHYS%OPT%GMU0, YDMF_PHYS%OPT%GMU0_MIN, &
  & YDMF_PHYS%OPT%GMU0_MAX, YDMF_PHYS%OPT%GDEOTI, YDMF_PHYS%OPT%GDEOTI2, YDMF_PHYS%OPT%GUEOTI, YDMF_PHYS%OPT%GUEOTI2,  &
  & YDMF_PHYS%OPT%GEOLT, YDMF_PHYS%OPT%GEOXT, YDMF_PHYS%OPT%GRPROX, YDMF_PHYS%OPT%GMIXP, YDMF_PHYS%OPT%GFLUXC,         &
  & YDMF_PHYS%OPT%GRSURF, YDMF_PHYS_SURF%GSD_VD%PSUND, ZFRSO, ZFRTH, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRTHC,         &
  & ZFRSODS, YDMF_PHYS%OUT%FRSOPS, ZFRSOLU, YDMF_PHYS%OUT%FRTHDS, ZAER)

  ! daand: radflex
  ! store for further calculations and diagnostics
  ! warning : pointers. REK
  YDMF_PHYS%OUT%FRSO(:,:,1)=ZFRSO
  YDMF_PHYS%OUT%FRTH(:,:,1)=ZFRTH

  ! extract surface fluxes
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    YDMF_PHYS%OUT%FRSODS(JLON)=ZFRSODS(JLON)+YDMF_PHYS%OUT%FRSOPS(JLON)  ! downward surface sw flux
  ENDDO

  IF (LLMSE) THEN
    IF (LHLRADUPD) THEN
      DO JSW = 1,NSW
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZZS_FSWDIR(JLON,JSW) = YDMF_PHYS%OUT%FRSOPS(JLON)*RSUN2(JSW)
          ZZS_FSWDIF(JLON,JSW) = ZFRSODS(JLON)*RSUN2(JSW)
         ENDDO
      ENDDO
    ELSE
      ZZS_FSWDIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)=YDMF_PHYS%OUT%FRSOPS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA) ! direct surface swdn flux
      ZZS_FSWDIF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)=ZFRSODS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA) ! diffuse surface swdn flux
    ENDIF
  ENDIF

  ! Estimate of the direct normal irradiance, with securities
  YDMF_PHYS%OUT%FRSDNI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS%OUT%FRSOPS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
  DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
    IF (ZRDG_MU0(JLON) > 3.0E-02_JPRB) THEN
      YDMF_PHYS%OUT%FRSDNI(JLON)=YDMF_PHYS%OUT%FRSOPS(JLON)/ZRDG_MU0(JLON)
    ENDIF
  ENDDO
  DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
    YDMF_PHYS%OUT%FRSDNI(JLON)=MAX(0.0_JPRB,YDMF_PHYS%OUT%FRSDNI(JLON))
  ENDDO

  IF (LFLEXDIA) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, YDMF_PHYS%OUT%FRSO(:, :, 1), 'FCTRAYSO', YDDDH&
      &             )
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, YDMF_PHYS%OUT%FRTH(:, :, 1), 'FCTRAYSO', YDDDH&
      &             )
  ENDIF

ENDIF

IF (.NOT.(LRAY.AND.NRAY == 2.AND.LRADFLEX).AND..NOT.LRAYFM) THEN
  DO JSW = 1,NSW
    ZZS_FSWDIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW) = 0._JPRB
    ZZS_FSWDIF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW) = 0._JPRB
  ENDDO
  YDMF_PHYS%OUT%FRSOPS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
ENDIF

IF (LFLEXDIA) THEN
  CALL ARO_STARTBU( YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, NRR, NGFL_EXT, ZRHODJM__(:, 1:YDCPG_OPTS%KFLEVG),     &
  & ZUS__(:, 1:YDCPG_OPTS%KFLEVG), ZVS__(:, 1:YDCPG_OPTS%KFLEVG), ZWS__(:, 1:YDCPG_OPTS%KFLEVG), ZTHS__(:, 1:YDCPG_OPTS%KFLEVG), &
  & ZRS_, ZTKES_, YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH)
ENDIF


!    ------------------------------------------------------------------
!     7 - CONVECTION.
!     --------------------------------------------------------------------

IF(LKFBCONV) THEN

  ! No swapp needed becaus IN and OUT are not needed simultaneously

  CALL BRI2ACCONV(YDMODEL%YRML_PHY_MF, YDGEOMETRY%YREGEO, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_BNDS%KFDIA,                          &
  & YDCPG_OPTS%KFLEVG, YDVARS%GEOMETRY%GM%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZPABSM__(:, 1:YDCPG_OPTS%KFLEVG),                       &
  & ZZZ_F_, ZTM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, :), ZRM_(:, :, 1), ZRM_(:, :, 2), ZRM_(:, :, 4), ZRHODREFM__(:, 1:YDCPG_OPTS%KFLEVG), &
  & ZUM__(:, 1:YDCPG_OPTS%KFLEVG), ZVM__(:, 1:YDCPG_OPTS%KFLEVG), ZWM__(:, 1:YDCPG_OPTS%KFLEVG),                                         &
  & ZMFS_, ZCVTENDT_, ZCVTENDRV_, ZCVTENDRC_, ZCVTENDRI_, ZCVTENDPR_, ZCVTENDPRS_   )

  IF(MOD(YDCPG_OPTS%NSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)"Pluie conv au sol", ZCVTENDPR_(NPTP), &
     & MAXVAL(ZCVTENDPR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)) ,MINVAL(ZCVTENDPR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA))
  ENDIF

  DO JLEV = 1,YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
      ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV) + ZCVTENDT_(JLON,JLEV)
      ZTEND_Q(JLON,JLEV) = ZTEND_Q(JLON,JLEV) + ZCVTENDRV_(JLON,JLEV)*ZQDM(JLON,JLEV)  
      ZTEND_L(JLON,JLEV) = ZTEND_L(JLON,JLEV) + ZCVTENDRC_(JLON,JLEV)*ZQDM(JLON,JLEV)  
      ZTEND_I(JLON,JLEV) = ZTEND_I(JLON,JLEV) + ZCVTENDRI_(JLON,JLEV)*ZQDM(JLON,JLEV)  
      ZRS_(JLON,JLEV,1)=ZRS_(JLON,JLEV,1)+ZCVTENDRV_(JLON,JLEV)
      ZRS_(JLON,JLEV,2)=ZRS_(JLON,JLEV,2)+ZCVTENDRC_(JLON,JLEV)
      ZRS_(JLON,JLEV,4)=ZRS_(JLON,JLEV,4)+ZCVTENDRI_(JLON,JLEV)
      ZTHS__(JLON,JLEV)=ZTHS__(JLON,JLEV)+ZCVTENDT_(JLON,JLEV)*(RATM/YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(JLON,JLEV))**(RD/RCPD)  
    ENDDO
  ENDDO
  DO JLON =YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
    ZINPRR_(JLON)=ZINPRR_(JLON)+ZCVTENDPR_(JLON)-ZCVTENDPRS_(JLON)
    ZACPRR_(JLON)=ZACPRR_(JLON)+(ZCVTENDPR_(JLON)-ZCVTENDPRS_(JLON))*YDCPG_OPTS%ZDTPHY
    ZINPRS_(JLON)=ZINPRS_(JLON)+ZCVTENDPRS_(JLON)
    ZACPRS_(JLON)=ZACPRS_(JLON)+ZCVTENDPRS_(JLON)*YDCPG_OPTS%ZDTPHY
  ENDDO
  ! avance temporelle et inversion niveau pour ZMFS_
  ! on utilise PSIGS pour le flux de masse pour la condensation sous maille
  ! car PSIGS n est utilise que si LOSIGMAS=T
  IF (PHYEX%NEBN%LSUBG_COND.AND..NOT.PHYEX%NEBN%LSIGMAS) THEN
    YDVARS%SRC%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=ZMFS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ENDIF
  IF(MOD(YDCPG_OPTS%NSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)"aps CONV, TENRV, TENRC, TENRI"
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,*)ZTEND_Q(NPTP,JLEV),ZTEND_L(NPTP,JLEV),ZTEND_I(NPTP,JLEV)
    ENDDO
  ENDIF
  CALL ARO_CONVBU(YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, NRR, ZRHODJM__(:, 1:YDCPG_OPTS%KFLEVG), ZRS_, &
  & ZTHS__(:, 1:YDCPG_OPTS%KFLEVG), YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH)

ENDIF

!    ------------------------------------------------------------------
!     8 - SURFACE.
!     --------------------------------------------------------------------

IF (LLMSE) THEN
! Initialisations

  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZZS_(JLON)=YDVARS%GEOMETRY%OROG%T0(JLON)*ZINVG 
  ENDDO
  DO JLEV = 1,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZDEPTH_HEIGHT_(JLON,JLEV)=ZZZ_F_(JLON,JLEV)-ZZS_(JLON)
    ENDDO
  ENDDO
  IF (MINVAL(ZDEPTH_HEIGHT_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,IKB)) <= 0._JPRB) THEN
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      IF (ZDEPTH_HEIGHT_(JLON,IKB) <= 0._JPRB) THEN
        WRITE (NULOUT,*)'sous apl_arome pb height en', JLON,ZAPHIFM(JLON,YDCPG_OPTS%KFLEVG),YDVARS%GEOMETRY%OROG%T0(JLON)
      ENDIF
    ENDDO
  ENDIF
  ! Can't use a section of pointer. An explicit copy shows, by the way, that a copy is needed
  ! because data is not contiguous. REK
  ZSVMB_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NGFL_EXT)=ZSVM_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,IKB,1:NGFL_EXT)

  IF (LLMSE_PARAM) THEN

      CALL ARO_GROUND_PARAM( YDCPG_BNDS%KBL, YDCPG_OPTS%KGPCOMP, YDCPG_BNDS%KFDIA, YDCPG_BNDS%KIDIA,                                                                     &
      & YDCPG_BNDS%KFDIA, YDCPG_OPTS%NSTEP, NRR, NSW, NGFL_EXT, NDGUNG, NDGUXG, NDLUNG, NDLUXG,                                                                          &
      & LSURFEX_KFROM, LMPA, CCOUPLING, LLXFUMSE, NINDAT, ZRHGMT, ZSTATI, RSOVR, RCODEC, RSIDEC, YDVARS%GEOMETRY%RINDX%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),            &
      & YDVARS%GEOMETRY%RINDY%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZUM__(:, IKB), ZVM__(:, IKB), ZTM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG),            &
      & ZRM_(:, IKB, 1), ZSVMB_, RCARDI, ZRHODREFM__(:, IKB), ZPABSM__(:, IKB), YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG), &
      & ZDTMSE, ZDEPTH_HEIGHT_(:, IKB), ZZS_, XZSEPS, ZRDG_MU0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZRDG_MU0N(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                         &
      & YDVARS%GEOMETRY%GELAM%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDVARS%GEOMETRY%GEMU%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                         &
      & XSW_BANDS, ZINPRR_NOTINCR_, ZINPRS_NOTINCR_, ZINPRG_NOTINCR_, YDMF_PHYS%OUT%FRTHDS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                           &
      & ZZS_FSWDIF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:NSW), ZZS_FSWDIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:NSW),                                                      &
      & ZCFAQ_, ZCFATH_, ZCFAU_, ZCFBQ_, ZCFBTH_, ZCFBU_, ZCFBV_, ZSFTH_, ZSFRV_, ZSFSV_, ZSFCO2_,                                                                       &
      & ZSFU_, ZSFV_, ZALBP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:NSW), ZALBD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:NSW),                                                  &
      & ZEMIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZTSURF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%FRTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, 1) &
      &     )

  ENDIF

  IF (LRCO2) THEN
      ZSFSV_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,NSV_CO2)= ZSFCO2_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
!print*,' FLUX CO2 =', MINVAL(ZSFSV_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,NSV_CO2)),&
!                    & MAXVAL(ZSFSV_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,NSV_CO2))
  ENDIF

!!!!! TEST DDH ATTENTION
!ZSFRV_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA) = 0._JPRB

  IF (LLMSE_DIAG) THEN

      CALL ARO_GROUND_DIAG( YDCPG_BNDS%KBL, YDCPG_OPTS%KGPCOMP, YDCPG_BNDS%KFDIA, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,                             &
      & YDCPG_OPTS%KFLEVG, IKL, NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, ZZS_, ZSFRV_, ZUM__(:, IKTB:IKTE),                                 &
      & ZVM__(:, IKTB:IKTE), ZDEPTH_HEIGHT_(:, IKTB:IKTE), YDMF_PHYS%OUT%FRTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, 1),           &
      & YDMF_PHYS%OUT%FRSO(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, 1), YDVARS%GEOMETRY%RINDX%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), &
      & YDVARS%GEOMETRY%RINDY%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZQS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                                      &
      & ZGZ0_, ZGZ0H_, YDMF_PHYS%OUT%TCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%QCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),              &
      & YDMF_PHYS%OUT%RHCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%UCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                            &
      & YDMF_PHYS%OUT%VCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%NUCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                            &
      & YDMF_PHYS%OUT%NVCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS%OUT%FCLL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1),                         &
      & YDMF_PHYS%OUT%FCLN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1), YDMF_PHYS%OUT%FEVL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1),                       &
      & YDMF_PHYS%OUT%FEVN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1), ZSSO_STDEV_, YDMF_PHYS_SURF%GSP_SG%PF_T1,                                       &
      & ZBUDTH_, ZBUDSO_, ZFCLL_, ZTOWNS_, ZCD_, YDMF_PHYS%OUT%SIC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA))

      CALL ARO_GROUND_DIAG_2ISBA( YDCPG_BNDS%KBL, YDCPG_OPTS%KGPCOMP, YDCPG_BNDS%KFDIA, YDCPG_BNDS%KIDIA,                                   &
      & YDCPG_BNDS%KFDIA, NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, YDVARS%GEOMETRY%RINDX%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),       &
      & YDVARS%GEOMETRY%RINDY%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS_SURF%GSD_VF%PLSM,                                            &
      & ZDUMMY1, ZDUMMY1, ZDUMMY1, ZTSURF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS_SURF%GSP_SG%PF_T1,                                  &
      & ZTP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZWS2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZWP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), &
      & ZWSI2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZWPI2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZWR2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), &
      & ZDUMMY1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDMF_PHYS_SURF%GSP_SG%PR_T1, ZDUMMY1 )
 
  ENDIF


!* Compute PBL-diagnostics

   ZCAPE(:)=0._JPRB
   ZDCAPE(:)=0._JPRB
   CALL ACCLDIA(YDCST, YDCPG_OPTS%LXCLP, YDCPG_OPTS%LXTGST, YDCPG_OPTS%LXXGST, YDPHY, YDMODEL%YRML_PHY_MF%YRPHY2, YDTOPH, YDCPG_BNDS%KIDIA,            &
   & YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS_BASE_STATE%U(:, 1:YDCPG_OPTS%KFLEVG),             &
   & YDMF_PHYS_BASE_STATE%V(:, 1:YDCPG_OPTS%KFLEVG), ZCAPE, ZDCAPE, ZTKEM(:, 1:YDCPG_OPTS%KFLEVG), YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(:, 1:YDCPG_OPTS%KFLEVG), &
   & YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST, YDMF_PHYS%OUT%CLPH, ICLPH)

   YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=MIN(XMAXLM,MAX(XMINLM,YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA))) 

   CALL ACVISIH(YDCST, YDVISI, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, &
   & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF,                  &
   & ZTM, ZRHM, ZQCM, ZQIM, ZQRM, ZQSM, ZQGM, YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD, YDMF_PHYS%OUT%MXCLWC                 &
   &                         )

ELSE

  ZSFSV_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)=0._JPRB

ENDIF    !  <== End block "IF (LMSE)"

!*            IDEALIZED TURBULENT SURFACE FLUXES FOR SQUALL LINE CASE
!                --------------------------------------------------------

IF (LSQUALL.AND.LTURB) THEN
  ! on n'a besoin que d'un flux sur V (U est nul).
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
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
  IF (PHYEX%PARAM_MFSHALLN%CMF_UPDRAFT=='DUAL') THEN
    ! Updraft computation from EDMF/ECMWF dual proposal
    ! Version May 2007
    !
    ! The following routine  are using arrays with the vertical Arpege/IFS fashion (as in the radiation scheme)

    IDRAFT = 2 ! beginning of the loop for MF tendency equation
               ! only 2 and 3 are used for tendency computation in ARO_SHALLOW_MF
    INDRAFT=3   ! 1 for test, 2 for dry, 3 for wet

    IF (YDCPG_OPTS%KMAXDRAFT < INDRAFT) THEN
      CALL ABOR1('APL_AROME : KMAXDRAFT TOO SMALL !')
    ENDIF

    DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
      ZZS_FTH_(JLON)=-1._JPRB*ZSFTH_(JLON)*(YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(JLON,YDCPG_OPTS%KFLEVG)*ZINVATM)**(ZRSCP)
      ZZS_FRV_(JLON)=-1._JPRB*ZSFRV_(JLON)
    ENDDO
    ZZS_FU_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZSFU_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
    ZZS_FV_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZSFV_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      ZZEXNREFM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZEXNREFM_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
    ENDDO

    ! If OCND2 and non-LIMA microphysics use water or mixed-phase cloud fraction for VDFHGHTHL, else use cloud fraction YDVARS%A%T1
    IF (PHYEX%PARAM_ICEN%LOCND2.AND.LMICRO.AND.CMICRO/='LIMA') THEN
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        ZZWCLDFR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV) = ZWCLDFR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
      ENDDO
    ELSE
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        ZZWCLDFR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV) = YDVARS%A%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
      ENDDO
    ENDIF

    !  IF LHARATU=TRUE then TKE at t-dt is needed as input for vdfexcuhl so fill ZTKEEDMF with t-1 value  from PTKEM

    IF (PHYEX%TURBN%LHARAT) THEN
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        ZTKEEDMF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZTKEM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
        ZLENGTH_M(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.01_JPRB
        ZLENGTH_H(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.01_JPRB
      ENDDO
      IF (MAXVAL(ZTKEM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)) > 3300._JPRB) THEN
        DO JLEV=1, YDCPG_OPTS%KFLEVG
          DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
            IF (ZTKEM(JLON,JLEV) > 3300._JPRB) THEN
              WRITE (NULOUT,*) 'TKE > 3300 ! '
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    CALL VDFHGHTHL(YDMODEL%YRML_PHY_G%YRVDF, YDMODEL%YRML_PHY_SLIN%YREPHLI, YDMODEL%YRML_PHY_EC%YRECUMF,          &
    & YDMODEL%YRML_PHY_EC%YREPHY, YDPARAR, YDCPG_OPTS%NSTEP, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, &
    & YDCPG_OPTS%KFLEVG, INDRAFT, YDCPG_OPTS%ZDTPHY, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, ZTM, ZQVM,              &
    & ZQCM, ZQIM, ZZWCLDFR, YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF,             &
    & ZAPHIFM, ZAPHIM, ZZEXNREFM, ZZS_FTH_, ZZS_FRV_, ZZS_FU_, ZZS_FV_, ZMF_UP, ZTHETAL_UP, ZQT_UP, ZTHTV_UP,                &
    & ZQC_UP, ZQI_UP, ZU_UP, ZV_UP, &
    & ZSPP_ALL%YSPP_CLDDPTH,ZSPP_ALL%YSPP_CLDDPTHDP, &
    & ZSPP_ALL%YSPP_RFAC_TWOC,ZSPP_ALL%YSPP_RZC_H,ZSPP_ALL%YSPP_RZL_INF, &
    & ZTENDQVUP, ZTENDTUP, ZSURFPREP,             &
    & ZSURFSNOW, ZUPGENL, ZUPGENN, ZCLFR, ZLENGTH_M, ZLENGTH_H, ZTKEEDMF)


    !  tendtup, tendqvup  tendencies for non-conserved AROME
    !  variables due to updraft precipitation/snow (and its evaporation)
    DO JLEV = 2 ,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV) + ZTENDTUP(JLON,JLEV)
        ZTEND_Q(JLON,JLEV)=ZTEND_Q(JLON,JLEV) + ZTENDQVUP(JLON,JLEV)
      ENDDO
    ENDDO

    IF (LTOTPREC.OR.LTOTPRECL) THEN
      !Add rain and snow tendencies from the sub-grid scheme to tendencies and sources,
      !at all vertical levels, instead of diagnosing only surface precip.
      ZSURFPREP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.0_JPRB
      ZSURFSNOW(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.0_JPRB
      DO JLEV= 1, YDCPG_OPTS%KFLEVG
        DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
          !Add rain and snow to sources:
          ZRS_(JLON,JLEV,3)=ZRS_(JLON,JLEV,3)+ZUPGENL(JLON,JLEV)
          ZRS_(JLON,JLEV,5)=ZRS_(JLON,JLEV,5)+ZUPGENN(JLON,JLEV)
          ZTHS__(JLON,JLEV)=ZTHS__(JLON,JLEV)+ZTENDTUP(JLON,JLEV)*(RATM/&
           & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(JLON,JLEV))**(RD/RCPD)
          !Update rain/snow tendencies:
          ZTEND_R(JLON,JLEV)=ZTEND_R(JLON,JLEV)+ZUPGENL(JLON,JLEV)
          ZTEND_S(JLON,JLEV)=ZTEND_S(JLON,JLEV)+ZUPGENN(JLON,JLEV)
        ENDDO
      ENDDO 
    ENDIF

  ELSE
    IDRAFT=3 ! only a wet updraft
    INDRAFT=1
    ZSURFPREP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
    ZSURFSNOW(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
  ENDIF

  DO JDRAFT=IDRAFT,3

    ! No need to swapp because IN and OUT are never needed simultaneously

    !!! Call mass fluxes computations
    ! If CMF_UPDRAFT='DUAL', the updraft characteritics are already computed and will be passed as inputs of SHALLOW_MF
    ! if not, they will be computed in SHALLOW_MF itself (from Méso-NH type routines)

    ! JDRAFT=2 : dry updraft
    ! JDRAFT=3 : wet updraft

    IF (PHYEX%PARAM_MFSHALLN%CMF_UPDRAFT=='DUAL') THEN
      ! Goes from one of the updraft from the IFS level world to the Méso-NH level world
      ! go from q to r)
      DO JLEV = 1,YDCPG_OPTS%KFLEVG
        DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
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
      ZZW_UP_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:IKT)=0._JPRB
      ZZFRAC_UP_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:IKT)=0._JPRB
      IF (PHYEX%TURBN%LHARAT) THEN
        DO JLEV = 1,YDCPG_OPTS%KFLEVG
          DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            ZLENGTHM__(JLON,JLEV) = MAX(0.01_JPRB,ZLENGTH_M(JLON,JLEV))
            ZLENGTHH__(JLON,JLEV) = MAX(0.01_JPRB,ZLENGTH_H(JLON,JLEV))
            ! TKE should be bigger than a minimum value:
            ZTKEEDMFS(JLON,JLEV) = MAX(ZTKEEDMF(JLON,JLEV),PPTKEMIN)*ZINVDT
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    IF(MOD(YDCPG_OPTS%NSTEP+1,NPRINTFR)==0) THEN
      WRITE(NULOUT,*)"apres surface zsfth zsfrv",ZSFTH_(NPTP),ZSFRV_(NPTP)
    ENDIF

    DO JLEV = 1, YDCPG_OPTS%KFLEVG 
      ZRC_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
      ZRI_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
      ZCF_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
    ENDDO

    IF (JDRAFT == IDRAFT) THEN
      ! Fill the sum at the first iteration
      ZARG_FLXZTHVMF_ => ZFLXZTHVMF_SUM__(:,1:YDCPG_OPTS%KFLEVG)
    ELSE
      ! increment
      ZARG_FLXZTHVMF_ => ZFLXZTHVMF_(:,1:YDCPG_OPTS%KFLEVG)
    ENDIF

    CALL ARO_SHALLOW_MF (PHYEX, &
    & KKL=IKL, KLON=YDCPG_BNDS%KFDIA, KLEV=YDCPG_OPTS%KFLEVG, KFDIA=YDCPG_BNDS%KFDIA, KRR=NRR, KRRL=NRRL,                         &
    & KRRI=NRRI, KSV=NGFL_EXT, &
    & KSV_LGBEG=0, KSV_LGEND=0, PTSTEP=ZDT, &
    & PDX=YDGEOMETRY%YREGEO%EDELX, PDY=YDGEOMETRY%YREGEO%EDELY,                                                                     &
    & PZZ=ZZZ_, PZZF=ZZZ_F_, PDZZF=ZDZZ_F_, PRHODJ=ZRHODJM__(:, 1:YDCPG_OPTS%KFLEVG), PRHODREF=ZRHODREFM__(:, 1:YDCPG_OPTS%KFLEVG), &
    & PPABSM=ZPABSM__(:, 1:YDCPG_OPTS%KFLEVG), PEXNM=ZEXNREFM_, PSFTH=ZSFTH_, PSFRV=ZSFRV_, PTHM=ZTHM__(:, 1:YDCPG_OPTS%KFLEVG),    &
    & PRM=ZRM_, PUM=ZUM__(:, 1:YDCPG_OPTS%KFLEVG), PVM=ZVM__(:, 1:YDCPG_OPTS%KFLEVG), PTKEM=ZTKEM__(:, 1:YDCPG_OPTS%KFLEVG),        &
    & PSVM=ZSVM_, PDUDT_MF=ZMFUS_, PDVDT_MF=ZMFVS_, PDTHLDT_MF=ZTHLS_, PDRTDT_MF=ZRTS_, PDSVDT_MF=ZSVXXX_,                          &
    & PSIGMF=ZSIGMF_, PRC_MF=ZRC_MF_, PRI_MF=ZRI_MF_, PCF_MF=ZCF_MF_, PFLXZTHVMF=ZARG_FLXZTHVMF_, PTHL_UP=ZTHETAL_UP_,              &
    & PRT_UP= ZRT_UP_, PRV_UP=ZZRV_UP_, PRC_UP=ZRC_UP_, PRI_UP=ZRI_UP_, PU_UP=ZZU_UP_, PV_UP=ZZV_UP_,                               &
    & PTHV_UP=ZTHETAV_UP_, PW_UP=ZZW_UP_, PFRAC_UP=ZZFRAC_UP_, PEMF=ZMF_UP__(:, 1:YDCPG_OPTS%KFLEVG),                               &
    & YDDDH=YDDDH, YDLDDH=YDMODEL%YRML_DIAG%YRLDDH, YDMDDH=YDMODEL%YRML_DIAG%YRMDDH)

    !wc No variance due to dry updraft yet.
    ! Putting ZSIGMF to 0 for dry updraft might be obsolete 
    IF (PHYEX%NEBN%LSTATNW) THEN
      IF (JDRAFT .EQ. 2) THEN
        DO JLEV = 1,YDCPG_OPTS%KFLEVG
          ZSIGMF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
        ENDDO
      ENDIF
    ENDIF

    IF (JDRAFT > IDRAFT) THEN
      ! Add increment
      ZFLXZTHVMF_SUM__(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=ZFLXZTHVMF_SUM__(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)+ZFLXZTHVMF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
    ENDIF

    ! traitement des sorties pour repasser dans le monde Aladin

    IF ((PHYEX%PARAM_MFSHALLN%CMF_CLOUD=='DIRE'.OR.PHYEX%PARAM_MFSHALLN%CMF_CLOUD=='BIGA').AND.JDRAFT==3) THEN
      ! sauvegarde pour le schema de nuage
      DO JLEV = 1,YDCPG_OPTS%KFLEVG
        ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,1)=ZRC_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
        ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,3)=ZRI_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
        ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,2)=ZCF_MF_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
      ENDDO
    ENDIF
    DO JLEV = 1,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZUS__(JLON,JLEV)=ZUS__(JLON,JLEV)+ZMFUS_(JLON,JLEV)
        ZVS__(JLON,JLEV)=ZVS__(JLON,JLEV)+ZMFVS_(JLON,JLEV)
        ZTHS__(JLON,JLEV)=ZTHS__(JLON,JLEV)+ZTHLS_(JLON,JLEV)
        ZRS_(JLON,JLEV,1)=ZRS_(JLON,JLEV,1)+ZRTS_(JLON,JLEV)
        !calcul de tendance et inversion des niveaux pour le vent horizontal
        YDMF_PHYS%OUT%TENDU(JLON,JLEV)=YDMF_PHYS%OUT%TENDU(JLON,JLEV)+ZMFUS_(JLON,JLEV)
        YDMF_PHYS%OUT%TENDV(JLON,JLEV)=YDMF_PHYS%OUT%TENDV(JLON,JLEV)+ZMFVS_(JLON,JLEV)
        !conversion de la tendance de theta en tendance de T et inversion niveau
        ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+ZTHLS_(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
        ZTENDTT(JLON,JLEV)=ZTENDTT(JLON,JLEV)+ZTHLS_(JLON,JLEV)
        !inversion niveaux tendances des ri et conversion en qi en multipliant par qd
        ZTEND_Q(JLON,JLEV) = ZTEND_Q(JLON,JLEV)+ZRTS_(JLON,JLEV)*ZQDM(JLON,JLEV)
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
    ZSVSIN_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,1:NGFL_EXT)=ZSVS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,1:NGFL_EXT)
  ENDIF

  !prints
  IF(MOD(YDCPG_OPTS%NSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome U'
    WRITE(NULOUT,*)MAXVAL(ZUM__(:,IKB)), MINVAL(ZUM__(:,IKB))
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome V'
    WRITE(NULOUT,*)MAXVAL(ZVM__(:,IKB)), MINVAL(ZVM__(:,IKB))
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome W'
    WRITE(NULOUT,*)MAXVAL(ZWM__(:,IKB)), MINVAL(ZWM__(:,IKB))
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome TKE'
    WRITE(NULOUT,*)MAXVAL(ZTKEM__(:,IKB)), MINVAL(ZTKEM__(:,IKB))
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,*)JLEV,ZUM__(NPTP,JLEV),ZVM__(NPTP,JLEV),ZWM__(NPTP,JLEV),ZTKEM__(NPTP,JLEV)
    ENDDO
    WRITE(NULOUT,*)'u v w tke a S'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,*)JLEV,ZUS__(NPTP,JLEV),ZVS__(NPTP,JLEV),ZWS__(NPTP,JLEV),ZTKES_(NPTP,JLEV)
    ENDDO
    WRITE(NULOUT,*)'ZTHS__ avant turb'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,*)JLEV,ZTHS__(NPTP,JLEV)
    ENDDO
  ENDIF

!!$
!!$! Allocation des variables SV (NGFL_EXT + NLIMA)
!!$  KSV_TURB=NGFL_EXT+NLIMA
!!$!
!!$  IF (NGFL_EXT/=0) THEN
!!$     DO JGFL=1,NGFL_EXT
!!$        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
!!$           ZSFTURB(JLON,JGFL)=ZSFSV_(JLON,JGFL)
!!$           DO JLEV = 1, YDCPG_OPTS%KFLEVG
!!$              ZTURBM(JLON,JLEV,JGFL)=ZSVM_(JLON,1,JLEV,JGFL)
!!$              ZTURBS(JLON,JLEV,JGFL)=ZSVSIN_(JLON,1,JLEV,JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF
!!$!
!!$  IF (NLIMA/=0) THEN
!!$     DO JGFL=1,NLIMA
!!$        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
!!$           ZSFTURB(JLON,NGFL_EXT+JGFL)=0.
!!$           DO JLEV = 1, YDCPG_OPTS%KFLEVG
!!$              ZTURBM(JLON,JLEV,NGFL_EXT+JGFL)=ZLIMAM_(JLON,JLEV,JGFL)
!!$              ZTURBS(JLON,JLEV,NGFL_EXT+JGFL)=ZLIMASIN_(JLON,JLEV,JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF

  ! Input variable indeed. REK
  ZSFSVLIMA_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NLIMA)=0._JPRB

  ! 10.2 calcul TURB
  ZZTOP_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZAPHIM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,0)*ZINVG

  IF (LGRADHPHY) THEN
  !
  DO JGR=1,NGRADIENTS  
    DO JLEV = 1,YDCPG_OPTS%KFLEVG
        ZTURB3D__(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,JGR)=YDMF_PHYS%GRA%G(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,JGR)
      ENDDO
    ENDDO
  
  ENDIF

! Appel avec les arguments modifiés pour variables LIMA :
! KSV_TURB, ZSFTURB, ZTURBM, ZTURBS, ZTENDSV_TURB
  CALL ARO_TURB_MNH(PHYEX, &
  & KKA=IKA, KKU=IKU, KKL=IKL, KLON=YDCPG_BNDS%KFDIA, KLEV=YDCPG_OPTS%KFLEVG, KRR=NRR, &
  & KRRL=NRRL, KRRI= NRRI, KSV=NLIMA, KGRADIENTS=NGRADIENTS, &
  & CMICRO=CMICRO, &
  & PTSTEP=ZDT, PZZ=ZZZ_, PZZF=ZZZ_F_, PZZTOP=ZZTOP_, PRHODJ=ZRHODJM__, PTHVREF=ZTHVREFM__,                          &
  & PSFTH=ZSFTH_, PSFRV=ZSFRV_, PSFSV=ZSFSVLIMA_, PSFU=ZSFU_,                                                        &
  & PSFV=ZSFV_, PPABSM=ZPABSM__, PUM=ZUM__, PVM=ZVM__, PWM=ZWM__, PTKEM=ZTKEM__, PEPSM=ZEPSM, PSVM=ZLIMAM_,          &
  & PSRCM=ZSRCS__, PTHM=ZTHM__, PRM=ZRM_, PRUS=ZUS__, PRVS=ZVS__, PRWS=ZWS__, PRTHS=ZTHS__, PRRS=ZRS_,               &
  & PRSVSIN=ZLIMASIN_, PRSVS=ZLIMAS_, PRTKES=ZTKES_, PRTKES_OUT=ZTKES_OUT__, PREPSS=ZEPSS, PHGRAD=ZTURB3D__,         &
  & PSIGS=ZSIGS__, PFLXZTHVMF=ZFLXZTHVMF_SUM__, PLENGTHM=ZLENGTHM__, PLENGTHH=ZLENGTHH__,    &
  & MFMOIST=ZMF_UP__, PDRUS_TURB=ZTENDU_TURB__, PDRVS_TURB=ZTENDV_TURB__, PDRTHLS_TURB=ZTENDTHL_TURB__,              &
  & PDRRTS_TURB=ZTENDRT_TURB__, PDRSVS_TURB=ZTENDSV_TURBLIMA_, PDP=ZDP__, PTP=ZTP__, PTPMF=ZTPMF__, PTDIFF=ZTDIFF__, &
  & PTDISS=ZTDISS__, PEDR=ZEDR__, YDDDH=YDDDH, YDLDDH=YDMODEL%YRML_DIAG%YRLDDH, YDMDDH=YDMODEL%YRML_DIAG%YRMDDH      &
  &                                                 )


! Séparation des variables SV (NGFL_EXT + NLIMA)
!!$  IF (NGFL_EXT/=0) THEN
!!$     DO JGFL=1,NGFL_EXT
!!$        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
!!$           ZSFSV_(JLON,JGFL)=ZSFTURB(JLON,JGFL)
!!$           DO JLEV = 1, YDCPG_OPTS%KFLEVG
!!$              ZSVM_(JLON,1,JLEV,JGFL)=ZTURBM(JLON,JLEV,JGFL)
!!$              ZSVS_(JLON,1,JLEV,JGFL)=ZTURBS(JLON,JLEV,JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF
!!$!
!!$  IF (NLIMA/=0) THEN
!!$     DO JGFL=1,NLIMA
!!$        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
!!$           DO JLEV = 1, YDCPG_OPTS%KFLEVG
!!$              ZLIMAM_(JLON,JLEV,JGFL)=ZTURBM(JLON,JLEV,NGFL_EXT+JGFL)
!!$              ZLIMAS_(JLON,JLEV,JGFL)=ZTURBS(JLON,JLEV,NGFL_EXT+JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF


  DO JLEV = 1 , YDCPG_OPTS%KFLEVG
     YDMF_PHYS%OUT%EDR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZEDR__(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
  ENDDO
   
  IF (LFLEXDIA) THEN
    DO JLEV = 1,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZDP__(JLON,JLEV)=ZDP__(JLON,JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
        ZTP__(JLON,JLEV)=(ZTP__(JLON,JLEV)-ZTPMF__(JLON,JLEV))*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
        ZTPMF__(JLON,JLEV)=ZTPMF__(JLON,JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
        ZTDIFF__(JLON,JLEV)=ZTDIFF__(JLON,JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
        ZTDISS__(JLON,JLEV)=ZTDISS__(JLON,JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
      ENDDO
    ENDDO
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZDP__(:, 1:YDCPG_OPTS%KFLEVG), 'TKEPRDY', &
      & YDDDH            )
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTP__(:, 1:YDCPG_OPTS%KFLEVG), 'TKEPRTH', &
      & YDDDH            )
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTPMF__(:, 1:YDCPG_OPTS%KFLEVG), 'TKEPRTHMF', &
      & YDDDH)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTDIFF__(:, 1:YDCPG_OPTS%KFLEVG), 'TKEDIFF', &
      & YDDDH)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTDISS__(:, 1:YDCPG_OPTS%KFLEVG), 'TKEDISS', &
      & YDDDH)

  ENDIF 

  IF(MOD(YDCPG_OPTS%NSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'u v w a S apres turb'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,*)JLEV,ZUS__(NPTP,JLEV),ZVS__(NPTP,JLEV),ZWS__(NPTP,JLEV),ZTKES_OUT__(NPTP,JLEV)
    ENDDO
    WRITE(NULOUT,*)'THS TKES SIGS apres turb'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,*)JLEV,ZTHS__(NPTP,JLEV),ZTKES_OUT__(NPTP,JLEV),ZSIGS__(NPTP,JLEV)
    ENDDO
  ENDIF

  ! avance temporelle et inversion niveau pour ZSIGS__
  IF (PHYEX%NEBN%LSUBG_COND .AND. PHYEX%NEBN%LSIGMAS) THEN
    IF (PHYEX%PARAM_MFSHALLN%CMF_CLOUD=='DIRE'.OR.PHYEX%PARAM_MFSHALLN%CMF_CLOUD=='BIGA'.OR. &
       &PHYEX%PARAM_MFSHALLN%CMF_CLOUD=='NONE') THEN
      DO JLEV = 1,YDCPG_OPTS%KFLEVG
        YDVARS%SRC%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZSIGS__(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
      ENDDO
    ELSEIF (PHYEX%PARAM_MFSHALLN%CMF_CLOUD=='STAT') THEN
      DO JLEV = 1,YDCPG_OPTS%KFLEVG
        DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDVARS%SRC%T1(JLON,JLEV)=SQRT(ZSIGS__(JLON,JLEV)**2+ZSIGMF_(JLON,JLEV)**2 )
        ENDDO
      ENDDO
    ENDIF
  ENDIF


  !10.3. traitement des sorties pour repasser dans le monde Aladin
  !calcul de tendance et inversion des niveaux pour le vent horizontal et la TKE

  DO JLEV = 1,YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%TENDU(JLON,JLEV)=YDMF_PHYS%OUT%TENDU(JLON,JLEV)+ZTENDU_TURB__(JLON,JLEV)
      YDMF_PHYS%OUT%TENDV(JLON,JLEV)=YDMF_PHYS%OUT%TENDV(JLON,JLEV)+ZTENDV_TURB__(JLON,JLEV)
      ! for the moment, turbulence do not compute w tendency:
      ZTENDW(JLON,JLEV)=0.0_JPRB
      ! PTENDW(JLON,JLEV)+(ZWS__(JLON,JLEV)-&
      ! & ZWS_AVE(JLON,1,JLEV))
      !conversion de la tendance de theta en tendance de T et inversion niveau
      ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+ZTENDTHL_TURB__(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
  !inversion niveaux tendances des rv et conversion en qv en multipliant par qd
      ZTEND_Q(JLON,JLEV)= ZTEND_Q(JLON,JLEV)+ZTENDRT_TURB__(JLON,JLEV)*ZQDM(JLON,JLEV)
    ENDDO
  ENDDO


  IF (PHYEX%TURBN%LHARAT) THEN
    DO JLEV = 1,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZTENDTKE(JLON,JLEV)=ZTENDTKE(JLON,JLEV)+(ZTKEEDMFS(JLON,JLEV)-ZTKES_(JLON,JLEV))
      ENDDO
    ENDDO
  ELSE
    DO JLEV = 1,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
         ZTENDTKE(JLON,JLEV)=ZTENDTKE(JLON,JLEV)+(ZTKES_OUT__(JLON,JLEV)-ZTKES_(JLON,JLEV))
      ENDDO
    ENDDO
  ENDIF

  DO JGFL=1,NGFL_EXT
    DO JLEV = 1,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZTENDEXT(JLON,JLEV,JGFL)=ZTENDEXT(JLON,JLEV,JGFL)+(ZSVS_(JLON,JLEV,JGFL)-ZSVSIN_(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO

! Tendances LIMA
  DO JGFL=1,NLIMA
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZTENDLIMA(JLON,JLEV,JGFL)=ZTENDLIMA(JLON,JLEV,JGFL)+ZTENDSV_TURBLIMA_(JLON,JLEV,JGFL)
!        PTENDLIMA(JLON,JLEV,:)=PTENDLIMA(JLON,JLEV,:)+ (ZLIMAS_(JLON,JLEV,:)-ZLIMASIN_(JLON,JLEV,:))
      ENDDO
    ENDDO
  ENDDO

ENDIF
IF(LWINDFARM)THEN
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZP1EZDIAG(JLON,1,5)=0.0_JPRB
    IF(GFL_WKA(JLON,YDCPG_BNDS%KBL)>0.0)THEN
      JGFL=INT(GFL_WKA2(1,6,YDCPG_BNDS%KBL))*1_JPIM
      CALL ARO_WINDFARM(YDGEOMETRY,YDMODEL,YDVARS%GEOMETRY%OROG%T0(JLON), &
          & YDMF_PHYS_BASE_STATE%U(JLON,1:YDCPG_OPTS%KFLEVG),YDMF_PHYS_BASE_STATE%V(JLON,1:YDCPG_OPTS%KFLEVG),&
          & YDMF_PHYS_BASE_STATE%YCPG_PHY%W(JLON,1:YDCPG_OPTS%KFLEVG),ZTKEM(JLON,1:YDCPG_OPTS%KFLEVG),&
          & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,0:YDCPG_OPTS%KFLEVG),&
          & YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(JLON,1:YDCPG_OPTS%KFLEVG),ZTENDU_TURB__(JLON,1:YDCPG_OPTS%KFLEVG),&
          & ZTENDV_TURB__(JLON,1:YDCPG_OPTS%KFLEVG),ZTKES_OUT__(JLON,1:YDCPG_OPTS%KFLEVG),&
          & GFL_WKA(JLON,YDCPG_BNDS%KBL),GFL_WKA2(1:JGFL,1:6,YDCPG_BNDS%KBL),YDCPG_OPTS%KFLEVG,&
          & YDCPG_OPTS%NSTEP,JGFL,ZP1EZDIAG(JLON,1,5))
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        YDMF_PHYS%OUT%TENDU(JLON,JLEV)=YDMF_PHYS%OUT%TENDU(JLON,JLEV)+ZTENDU_TURB__(JLON,JLEV)
        YDMF_PHYS%OUT%TENDV(JLON,JLEV)=YDMF_PHYS%OUT%TENDV(JLON,JLEV)+ZTENDV_TURB__(JLON,JLEV)
        ZTENDTKE(JLON,JLEV)=ZTENDTKE(JLON,JLEV)+ZTKES_OUT__(JLON,JLEV)
        !WRITE(*,*) 'windfarm tendencies dU/dt,dV/dt,dTKE/dt'
        !WRITE(*,*)ZTENDU_TURB(JLON,1,JLEV),ZTENDV_TURB(JLON,1,JLEV),&
        !        & ZZI_TKES(JLON,1,JLEV)
      ENDDO
    ENDIF
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
  ZTHS__(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=ZTHSIN_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ZRS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,1:NRR)=ZRSIN_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,1:NRR)
  ! for now a copy is needed (see below, inside). I don't like than :-( REK
  ZLIMAS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,1:NLIMA)=ZLIMASIN_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,1:NLIMA)
     
  !prints
  IF (MOD(YDCPG_OPTS%NSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'avant rain_ice sous apl_arome'
    WRITE(NULOUT,*)'JLEV   ZZZ_F_      ZZZ_      ZRHODREF',&
     & '    ZRHODJ      ZPABSM__        ZTHSIN_       ZTHM__      '
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,7F10.3)')JLEV,ZZZ_F_(NPTP,JLEV),ZZZ_(NPTP,JLEV), ZRHODREFM__(NPTP,JLEV),&
       & ZRHODJM__(NPTP,JLEV), ZPABSM__(NPTP,JLEV), ZTHSIN_(NPTP,JLEV), ZTHM__(NPTP,JLEV)  
    ENDDO 
    WRITE(NULOUT,*)'JLEV        PDELPM        ZPABSM__         ZEXNREF','          ZSIGS__'
    DO JLEV=2,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,4f10.3)')JLEV, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(NPTP,JLEV),&
       & ZPABSM__(NPTP,JLEV),ZEXNREFM_(NPTP,JLEV),ZSIGS__(NPTP,JLEV)
    ENDDO
    WRITE(NULOUT,*)'JLEV    PTM       PRM          PCPM'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,3f10.3)')JLEV,ZTM(NPTP,YDCPG_OPTS%KFLEVG+1-JLEV), ZRHM(NPTP,YDCPG_OPTS%KFLEVG+1-JLEV) ,ZCPM(NPTP,YDCPG_OPTS%KFLEVG+1-JLEV)  
    ENDDO
    WRITE (NULOUT,*)'JLEV  rhoQv  rhoQc   rhoQr   rhoQi   rhoQs   rhoQg'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,ZRM_(NPTP,JLEV,1), ZRM_(NPTP,JLEV,2),&
       & ZRM_(NPTP,JLEV,3),ZRM_(NPTP,JLEV,4),ZRM_(NPTP,JLEV,5), ZRM_(NPTP,JLEV,6)
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZRSQv  ZRSQc   ZRSQr   ZRSQi   ZRSQs   ZRSQg'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,ZRS_(NPTP,JLEV,1), ZRS_(NPTP,JLEV,2),&
       & ZRSIN_(NPTP,JLEV,3),ZRSIN_(NPTP,JLEV,4),ZRSIN_(NPTP,JLEV,5), ZRSIN_(NPTP,JLEV,6)
    ENDDO
    WRITE(NULOUT,*)'ZDT=',ZDT
    WRITE(NULOUT,*)'NRR and co',NRR,YDCPG_OPTS%NSTEP+1
  ENDIF

  ZSEA_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.0_JPRB
  IF (LOLSMC) THEN
    DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
      IF (YDMF_PHYS_SURF%GSD_VF%PLSM(JLON) < 0.5) THEN
        ZSEA_(JLON) = 1.0_JPRB
      ENDIF
    ENDDO
  ENDIF

  IF (LOTOWNC) THEN
    ZTOWN_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA) = ZTOWNS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
  ELSE
    ZTOWN_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.0_JPRB
  ENDIF

  IF (CMICRO == 'LIMA') THEN

    IF (LTURB) THEN
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        DO JLEV=1,YDCPG_OPTS%KFLEVG
          ZWNU_(JLON,JLEV) = ZWM__(JLON,JLEV) + 0.66*SQRT(ZTKEM__(JLON,JLEV))
        ENDDO
      ENDDO
      ZPTRWNU_ => ZWNU_(1:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
    ELSE
      ZPTRWNU_ => ZWM__(1:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
    ENDIF
    CALL ARO_LIMA(PHYEX, YDCPG_OPTS%KFLEVG, IKU, IKL, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG,YDCPG_BNDS%KFDIA,NRR, NLIMA, &
    & ZDT, ZDZZ_, ZRHODJM__(:, 1:YDCPG_OPTS%KFLEVG), ZRHODREFM__(:, 1:YDCPG_OPTS%KFLEVG),               &
    & ZEXNREFM_, ZPABSM__(:, 1:YDCPG_OPTS%KFLEVG), ZPTRWNU_, ZDTHRAD_, ZTHM__(:, 1:YDCPG_OPTS%KFLEVG), ZRM_,                     &
    & ZLIMAM_, ZTHS__(:, 1:YDCPG_OPTS%KFLEVG), ZRS_, ZLIMAS_, ZEVAP_, ZINPRR_NOTINCR_,                                    &
    & ZINPRS_NOTINCR_, ZINPRG_NOTINCR_, ZINPRH_NOTINCR_, ZPFPR_, ZNEBMNH_, ZICEFR_, ZPRCFR_, YDDDH, YDMODEL%YRML_DIAG%YRLDDH,               &
    & YDMODEL%YRML_DIAG%YRMDDH)
  ELSE
    CALL ARO_RAIN_ICE (PHYEX, &
    & YDCPG_OPTS%KFLEVG,IKU,IKL,YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,YDCPG_BNDS%KFDIA,NRR,                     &
    & CMICRO, ZDT, ZDZZ_, &
    & ZRHODJM__(:, 1:YDCPG_OPTS%KFLEVG), ZRHODREFM__(:, 1:YDCPG_OPTS%KFLEVG), ZEXNREFM_, ZPABSM__(:, 1:YDCPG_OPTS%KFLEVG),      &
    & ZHLC_HRC_, ZHLC_HCF_, ZHLI_HRI_, ZHLI_HCF_,                                                                               &
    & ZTHM__(:, 1:YDCPG_OPTS%KFLEVG),                                                                                           &
    & ZRM_, ZSIGS__(:, 1:YDCPG_OPTS%KFLEVG), ZNEBMNH_, ZTHS__(:, 1:YDCPG_OPTS%KFLEVG), ZRS_, ZEVAP_,                            &
    & ZCIT_, ZSEA_, ZTOWN_, &
    & ZICLDFR_, ZWCLDFR_, ZSSIO_, ZSSIU_, ZIFR_, &
    & LKOGAN, LMODICEDEP,&
    & ZINPRR_NOTINCR_, ZINPRS_NOTINCR_, ZINPRG_NOTINCR_,                           &
    & ZINPRH_NOTINCR_, ZPFPR_, &
    & YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH,&
    & ZSPP_ALL%YSPP_ICENU,ZSPP_ALL%YSPP_KGN_ACON,ZSPP_ALL%YSPP_KGN_SBGR)
  ENDIF

  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZINPRR_(JLON)=ZINPRR_(JLON)+ZINPRR_NOTINCR_(JLON)
    ZINPRS_(JLON)=ZINPRS_(JLON)+ZINPRS_NOTINCR_(JLON)
    ZINPRG_(JLON)=ZINPRG_(JLON)+ZINPRG_NOTINCR_(JLON)
    ZINPRH_(JLON)=ZINPRH_(JLON)+ZINPRH_NOTINCR_(JLON)
  ENDDO

  !conversion de la tendance de theta en tendance de T et inversion niveau
  DO JLEV = 1,YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZTENDT(JLON,JLEV)= ZTENDT(JLON,JLEV)+(ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV))*ZEXNREFM_(JLON,JLEV)  
      ZTENDTT(JLON,JLEV)= ZTENDTT(JLON,JLEV)+ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV)
    ENDDO
  ENDDO

  !inversion niveaux tendances des ri et conversion en qi en multipliant par qd
  DO JR=1,NRR
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZTENDRA(JLON,JLEV,JR)=ZTENDRA(JLON,JLEV,JR)+(ZRS_(JLON,JLEV,JR)-ZRSIN_(JLON,JLEV,JR))*ZQDM(JLON,JLEV)  
      ENDDO
    ENDDO
  ENDDO

  ! Tendances des variables LIMA
  DO JGFL=1,NLIMA
    DO JLEV = 1, YDCPG_OPTS%KFLEVG
      DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZTENDLIMA(JLON,JLEV,JGFL)=ZTENDLIMA(JLON,JLEV,JGFL)+(ZLIMAS_(JLON,JLEV,JGFL)-ZLIMASIN_(JLON,JLEV,JGFL))  
      ENDDO
    ENDDO
  ENDDO

  IF (LINTFLEX) THEN
    !inversion of levels of upper-air precipitation
    DO JR=2,NRR ! no precip for qv
      ZFPR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,0,JR)=0._JPRB  ! zero precip at top of atmosphere
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        ZFPR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,JR)=ZPFPR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,JR)
      ENDDO
    ENDDO
  ENDIF

  !store cumulative 3D precipitations for mocage
  IF (LFPREC3D) THEN
    DO JR=2,NRR ! no precip for qv
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZP1EZDIAG(JLON,JLEV,4)=ZP1EZDIAG(JLON,JLEV,4)+ZPFPR_(JLON,JLEV,JR)*1000._JPRB*YDCPG_OPTS%ZDTPHY
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  !prints
  IF(MOD(YDCPG_OPTS%NSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'PTENDT en sortie de rain_ice'
    WRITE(NULOUT,*)'ZTHS__ en sortie de rain_ice'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,*)ZTENDT(NPTP,JLEV),ZTHS__(NPTP,JLEV)
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZTENDQv  ZTZNDQc   ZTENDQr   ZTENDQi' ,'ZTENDQs   ZTENDQg'
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,ZTEND_Q(NPTP,JLEV),ZTEND_L(NPTP,JLEV),&
       & ZTEND_R(NPTP,JLEV),ZTEND_I(NPTP,JLEV),ZTEND_S(NPTP,JLEV),ZTEND_G(NPTP,JLEV)  
    ENDDO
    WRITE (NULOUT,*) 'ZSRCS__ et ZNEBMNH_',MAXVAL(ZSRCS__),MAXVAL(ZNEBMNH_)
  ENDIF

  IF (LRDEPOS) THEN
    ! Swapp because IN and OUT will be needed simultaneously
    CALL SWAP_SVM
    CALL ARO_RAINAERO(YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, NGFL_EXT, NRR, YDCPG_OPTS%ZDTPHY, ZSVMIN_, ZZZ_, ZPABSM__(:, 1:YDCPG_OPTS%KFLEVG), &
    & ZTHM__(:, 1:YDCPG_OPTS%KFLEVG), ZRHODREFM__(:, 1:YDCPG_OPTS%KFLEVG), YDCPG_OPTS%NSTEP+1, ZRM_,                               &
    & ZEVAP_, YDPARAR%PHYEX%CLOUDPARN%NSPLITR, ZSVM_           )
    ! return to tendency
    DO JGFL=1,NGFL_EXT
      DO JLEV = 1,YDCPG_OPTS%KFLEVG
        DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZTENDEXT(JLON,JLEV,JGFL)=ZTENDEXT(JLON,JLEV,JGFL)+(ZSVM_(JLON,JLEV,JGFL)-ZSVMIN_(JLON,JLEV,JGFL))*ZINVDT
        ENDDO
      ENDDO
    ENDDO
  ENDIF ! LRDEPOS

ENDIF ! LMICRO

! start LHN F.Meier 2020 ******

LNUDGLHNREAD=.TRUE.
IF(MYPROC==1.AND.YDCPG_OPTS%NSTEP==1.AND.LNUDGLH)THEN
  CALL NUDGLHCLIMPROF(YDCPG_OPTS%KFLEVG, LNUDGLHNREAD)
ENDIF
! save accumulated precipitation for LHN
IF (LNUDGLH.AND.YDCPG_OPTS%NSTEP == NSTARTNUDGLH.AND.NSTARTNUDGLH > 0) THEN
  !IF(MYPROC==1) WRITE(NULOUT,*)'save precip for LHN - STEP:',KSTEP, &
  !  & 'NUDGINGINT:',NINTNUDGLH,'NSTARTNUDGLH:',NSTARTNUDGLH
  CALL NUDGLHPRECIP(YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, ZACPRR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), &
  & ZACPRS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZACPRG_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), YDCPG_BNDS%KBL           &
  &                                                                               )
ENDIF
ISTEP=YDCPG_OPTS%NSTEP-NSTARTNUDGLH
! if LNUDGLH and KSTEP in nudging interval
IF (LNUDGLH.AND.YDCPG_OPTS%NSTEP > NSTARTNUDGLH.AND.YDCPG_OPTS%NSTEP <= NSTOPNUDGLH) THEN
  ! safe LH profile for step before LHN step
  LLHN=.FALSE.
  IF(MOD(ISTEP+1,NINTNUDGLH)==0) THEN
    CALL NUDGLHPREP(YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, ZTENDTT, &
    & YDCPG_BNDS%KBL            )
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
    CALL NUDGLH(NGPTOT, NPROMA, NGPBLKS, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,                                                   &
    & ZACPRR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZACPRS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZACPRG_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                          &
    & ZTENDTT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), JLHSTEP, YDCPG_BNDS%KBL, ZEXNREFM_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), &
    & .TRUE., LLHN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), ZPABSM__(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG),                      &
    & ZDT, ZTHM__(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), ZRM_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG, :),                        &
    & ZQDM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), ZTEND_Q(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG),                               &
    & NRR, LNUDGLHNREAD)
    !IF(MYPROC==1) WRITE(NULOUT,*)'calling LH successful - convert TH to T and
    !add temperature tendency'
    ! add LHN tendency to physics tendency, limit LHN tendency
    ZMAXTEND=0.0_JPRB
    ZMINTEND=0.0_JPRB
    DO JLEV = 1,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        IF(LLHN(JLON,JLEV))THEN
          ZTENDTT(JLON,JLEV)=MAX(ZTENDTT(JLON,JLEV),RMINNUDGLH)
          ZTENDTT(JLON,JLEV)=MIN(ZTENDTT(JLON,JLEV),RMAXNUDGLH)
          ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+ZTENDTT(JLON,JLEV)*&
           & RAMPLIFY*ZEXNREFM_(JLON,JLEV)
          ZMINTEND=MIN(ZTENDTT(JLON,JLEV),ZMINTEND)
          ZMAXTEND=MAX(ZTENDTT(JLON,JLEV),ZMAXTEND)
          ! keep RH constant if LNUDGLHCOMPT=T
          IF(YDMF_PHYS_BASE_STATE%T(JLON,JLEV)>0.01_JPRB.AND.LNUDGLHCOMPT)THEN
            ZTEND_Q(JLON,JLEV)=ZTEND_Q(JLON,JLEV)+RLVTT/RV/((YDMF_PHYS_BASE_STATE%T(JLON,JLEV))**2._JPRB)* &
            & ZTENDTT(JLON,JLEV)*RAMPLIFY*ZEXNREFM_(JLON,JLEV)*ZQSAT(JLON,JLEV)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !IF(ZMINTEND<-0.01)WRITE(*,*)'ZMINTEND',ZMINTEND
    !IF(ZMAXTEND>0.01)WRITE(*,*)'ZMAXTEND',ZMAXTEND
    ! write LH profiles to array to save it for next time step
    CALL NUDGLHPREP(YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, ZTENDTT, &
    & YDCPG_BNDS%KBL            )
    IF(MYPROC==1) WRITE(NULOUT,*)'calling LH successful finished'
    ! use LHN factor again on following time steps depending on NTAUNUDGLH
  ELSEIF(MOD(ISTEP,NINTNUDGLH)<NTAUNUDGLH.AND.MOD(ISTEP,NINTNUDGLH)>0 &
      & .AND.ISTEP>NINTNUDGLH) THEN
    IF(MYPROC==1)THEN
      WRITE(NULOUT,*)'LH nudging applied-STEP:',YDCPG_OPTS%NSTEP,'NUDGINGINT:',NINTNUDGLH
      WRITE(NULOUT,*)'NTAUNUDGLH:',NTAUNUDGLH
    ENDIF
    ! get index for reading correctly most recent obs
    JLHSTEP=2+NINT(1.0_JPRB*(ISTEP-MOD(ISTEP,NINTNUDGLH))/(NTIMESPLITNUDGLH*NINTNUDGLH))
    !IF(MYPROC==1) WRITE(NULOUT,*)'observation array:',JLHSTEP
    ! call nudging routine to modify LHN profile where necessary
    ! LHN factor is not recalculated but might be damped by RDAMPNUDGLH
    CALL NUDGLH(NGPTOT, NPROMA, NGPBLKS, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,                                                   &
    & ZACPRR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZACPRS_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA), ZACPRG_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),                          &
    & ZTENDTT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), JLHSTEP, YDCPG_BNDS%KBL, ZEXNREFM_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), &
    & .FALSE., LLHN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), ZPABSM__(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG),                     &
    & ZDT, ZTHM__(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), ZRM_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG, :),                        &
    & ZQDM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG), ZTEND_Q(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA, 1:YDCPG_OPTS%KFLEVG),                               &
    & NRR, LNUDGLHNREAD)
    !IF(MYPROC==1) WRITE(NULOUT,*)'calling LH successful - convert TH to T and
    !add temperature tendency'
    ! add LHN tendency to physics tendency, limit LHN tendency
    DO JLEV = 1,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        IF(LLHN(JLON,JLEV))THEN
          ZTENDTT(JLON,JLEV)=MAX(ZTENDTT(JLON,JLEV),RMINNUDGLH)
          ZTENDTT(JLON,JLEV)=MIN(ZTENDTT(JLON,JLEV),RMAXNUDGLH)
          ZTENDT(JLON,JLEV)= ZTENDT(JLON,JLEV)+ZTENDTT(JLON,JLEV)*&
           & RAMPLIFY*ZEXNREFM_(JLON,JLEV)
          ZMINTEND=MIN(ZTENDTT(JLON,JLEV),ZMINTEND)
          ZMAXTEND=MAX(ZTENDTT(JLON,JLEV),ZMAXTEND)
          ! keep RH constant if LNUDGLHCOMPT=T
          IF(YDMF_PHYS_BASE_STATE%T(JLON,JLEV)>0.01_JPRB.AND.LNUDGLHCOMPT)THEN
             ZTEND_Q(JLON,JLEV)=ZTEND_Q(JLON,JLEV)+RLVTT/RV/((YDMF_PHYS_BASE_STATE%T(JLON,JLEV))**2._JPRB)*&
              & ZTENDTT(JLON,JLEV)*RAMPLIFY*ZEXNREFM_(JLON,JLEV)*ZQSAT(JLON,JLEV)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !IF(ZMAXTEND>0.01) WRITE(*,*)'ZMAXTEND',ZMAXTEND
    !IF(ZMINTEND<-0.01) WRITE(*,*)'ZMINTEND',ZMINTEND
    ! write LHN profiles to array for next timestep
    CALL NUDGLHPREP(YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, ZTENDTT, &
    & YDCPG_BNDS%KBL            )
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

  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZLAT_(JLON) = 180. * ASIN(YDVARS%GEOMETRY%GEMU%T0(JLON)) / (2.*ASIN(1.))
    ZLON_(JLON) = 180. * YDVARS%GEOMETRY%GELAM%T0(JLON) / (2.*ASIN(1.))
    ZZENITH_(JLON) = ACOS( ZRDG_MU0(JLON) )
    ZZS_(JLON)=YDVARS%GEOMETRY%OROG%T0(JLON)/RG
    ZALB_UV_(JLON)=ZALBP(JLON,1)
  ENDDO

  ! Swapp because IN and OUT will be needed simultaneously
  CALL SWAP_SVS

  DO JGFL=1,NGFL_EXT
    DO JLEV=1,YDCPG_OPTS%KFLEVG
      DO JLON= YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ! modify input
        ZSVSIN_(JLON,JLEV,JGFL)=MAX(0.0_JPRB, ZSVSIN_(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO
  IEZDIAG_CHEM=NGFL_EZDIAG-IOFF_MFSHAL+1
  CALL ARO_MNHC(ZSVSIN_, ZRHODREFM__(:, 1:YDCPG_OPTS%KFLEVG), YDCPG_OPTS%ZDTPHY, ZTHM__(:, 1:YDCPG_OPTS%KFLEVG), &
  & ZPABSM__(:, 1:YDCPG_OPTS%KFLEVG), ZRM_, ZLAT_, ZLON_, ZALB_UV_, ZZS_, ZZENITH_, ZZZ_, IYEAR,      &
  & IMONTH, IDAY, REAL(RHGMT, JPRB)+YDCPG_OPTS%ZDTPHY/2._JPRB, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, NGFL_EXT,    &
  & NRR, YDCPG_OPTS%NSTEP+1, NULOUT, IEZDIAG_CHEM, ZPEZDIAG_(:, :, IOFF_MFSHAL:NGFL_EZDIAG), ZSVS_ )
 
  ZP1EZDIAG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,IOFF_MFSHAL:NGFL_EZDIAG)=ZPEZDIAG_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG,IOFF_MFSHAL:NGFL_EZDIAG)

  !inversion niveau de la tendance des scalaires passifs
  DO JGFL=1,NGFL_EXT
    DO JLEV = 1,YDCPG_OPTS%KFLEVG
      DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZTENDEXT(JLON,JLEV,JGFL)=ZTENDEXT(JLON,JLEV,JGFL)+(ZSVS_(JLON,JLEV,JGFL)-ZSVSIN_(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO

ENDIF ! LUSECHEM

!    ------------------------------------------------------------------
!     13 - STOCHASTIC PHYSICS : PERTURB TENDENCIES
!     -----------------------------------------------------------------

IF(LSPSDT) THEN

  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZMULNOISE(JLON,1)=PGP2DSDT(JLON,1,1)       ! Use a single 2D pattern for all levels
  ENDDO

  ZDUMMY(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=0.0_JPRB               ! Dummy nonphys tendency for compatibility with ecmwf stochphy
  CALL SPPTEN (YDMODEL, YGFL, &
  & YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, 1, YDCPG_OPTS%ZDTPHY, &  ! In: block indices, physicstimestep
  & PTSL=YDMF_PHYS_BASE_STATE%T, PQSL=YDMF_PHYS_BASE_STATE%Q, PA=YDVARS%A%T1,       &  ! In: (T,Q,cloud) forsupersatcheck
  & PAP=YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF, PAPH=YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE, &  ! In: Pfull, Phalf
  & PDYN_U=ZDUMMY, PDYN_V=ZDUMMY, PDYN_T=ZDUMMY, PDYN_Q=ZDUMMY,                     &  ! In: dummy nonphys tendencies
  & PUNP_U=ZDUMMY, PUNP_V=ZDUMMY, PUNP_T=ZDUMMY, PUNP_Q=ZDUMMY,                     &  ! In: (u,v,t,qv) tendencies to perturb
  & PPHY_U=YDMF_PHYS%OUT%TENDU, PPHY_V=YDMF_PHYS%OUT%TENDV, PPHY_T=ZTENDT, PPHY_Q=ZTEND_Q,     &  ! In: (u,v,t,qv) tendencies to perturb
  & PMULNOISE=ZMULNOISE,                                                            &  ! In: stochphy 3D random multiplicative pattern (less one)  
  & PTENU=YDMF_PHYS%OUT%TENDU, PTENV=YDMF_PHYS%OUT%TENDV, PTENT=ZTENDT, PTENQ=ZTEND_Q )  ! Out: (u,v,t,qv) total perturbed tendencies

ENDIF

IF(LFORCENL.AND.(YDCPG_OPTS%NSTEP*(TSPHY/RHOUR)>=NFORCESTART).AND.&
              & (YDCPG_OPTS%NSTEP*(TSPHY/RHOUR)<=NFORCEEND)) THEN
  DO JLEV=1,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%TENDU(JLON,JLEV)=YDMF_PHYS%OUT%TENDU(JLON,JLEV)+AMAGSTOPH_CASBS*YDMF_PHYS%FOR%U(JLON,JLEV)
      YDMF_PHYS%OUT%TENDV(JLON,JLEV)=YDMF_PHYS%OUT%TENDV(JLON,JLEV)+AMAGSTOPH_CASBS*YDMF_PHYS%FOR%V(JLON,JLEV)
      ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+AMAGSTOPH_CASBS*YDMF_PHYS%FOR%T(JLON,JLEV)
      ZTEND_Q(JLON,JLEV)=ZTEND_Q(JLON,JLEV)+AMAGSTOPH_CASBS*YDMF_PHYS%FOR%Q(JLON,JLEV)
    ENDDO
  ENDDO
ENDIF

!    ------------------------------------------------------------------
!     14 - FINAL CALCULATIONS.
!     --------------------------------------------------------------------

!forcage pour declencher la ligne de grain
IF (LSQUALL) THEN
  IF (YDDYNA%LTWOTL) THEN
    ZDT2=2*ZDT
  ELSE
    ZDT2=ZDT
  ENDIF
  IF((YDCPG_OPTS%NSTEP+1)*ZDT2 < 600._JPRB) THEN
    WRITE(NULOUT, *)'refroidissement impose de',NREFROI1,' a ',NREFROI2
    DO JLEV=YDCPG_OPTS%KFLEVG,YDCPG_OPTS%KFLEVG-20,-1
      ZTENDT(NREFROI1:NREFROI2,JLEV)=-0.01_JPRB
    ENDDO
  ENDIF
ENDIF


!ecriture du buffer
IF(LLMSE.OR.LSFORCS) THEN
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    YDCPG_GPAR%INPRR(JLON)=ZINPRR_(JLON)+ZSURFPREP(JLON)/1000._JPRB
    YDCPG_GPAR%INPRS(JLON)=ZINPRS_(JLON)+ZSURFSNOW(JLON)/1000._JPRB
    YDCPG_GPAR%INPRG(JLON)=ZINPRG_(JLON)+ZINPRH_(JLON)
    YDCPG_GPAR%ACPRR(JLON)=YDCPG_GPAR%ACPRR(JLON)+(ZINPRR_(JLON)+ZSURFPREP(JLON)/1000._JPRB)*YDCPG_OPTS%ZDTPHY
    YDCPG_GPAR%ACPRS(JLON)=YDCPG_GPAR%ACPRS(JLON)+(ZINPRS_(JLON)+ZSURFSNOW(JLON)/1000._JPRB)*YDCPG_OPTS%ZDTPHY
    YDCPG_GPAR%ACPRG(JLON)=YDCPG_GPAR%ACPRG(JLON)+(ZINPRG_(JLON)+ZINPRH_(JLON))*YDCPG_OPTS%ZDTPHY
  ENDDO
  YDCPG_GPAR%VTS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZTSURF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
  YDCPG_GPAR%VEMIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZEMIS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
  YDCPG_GPAR%VQS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=ZQS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
  DO JSW=1,NSW
    YDCPG_GPAR%ALBDIR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=ZALBP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)
    YDCPG_GPAR%ALBSCA(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)=ZALBD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JSW)
  ENDDO
ENDIF

IF (LMUSCLFA) CALL ECR1D(NMUSCLFA, 'PCLCT_apl', YDCPG_MISC%CLCT, 1, YDCPG_OPTS%KLON)
! initialisations for CFU for Rainfalls
DO JLEV = 0,YDCPG_OPTS%KFLEVG
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ! conversion from m/s in mm/s
    YDMF_PHYS%OUT%FPLSL(JLON,JLEV)= ZINPRR_(JLON)*1000._JPRB+ZSURFPREP(JLON)
    YDMF_PHYS%OUT%FPLSN(JLON,JLEV)= ZINPRS_(JLON)*1000._JPRB+ZSURFSNOW(JLON)
    YDMF_PHYS%OUT%FPLSG(JLON,JLEV)= ZINPRG_(JLON)*1000._JPRB
    YDMF_PHYS%OUT%FPLSH(JLON,JLEV)= ZINPRH_(JLON)*1000._JPRB
    ! conversion in correct Unit for BADP (same as ALADIN)
    YDMF_PHYS%OUT%STRTU(JLON,JLEV)= ZSFU_(JLON)*ZRHODREFM__(JLON,IKB) 
    YDMF_PHYS%OUT%STRTV(JLON,JLEV)= ZSFV_(JLON)*ZRHODREFM__(JLON,IKB) 
  ENDDO
ENDDO
!Hail diagnostic
YDMF_PHYS%OUT%DIAGH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
IF (LXXDIAGH) THEN
  DO JLEV=1,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%DIAGH(JLON)=YDMF_PHYS%OUT%DIAGH(JLON)+ZQGM(JLON,JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO
ENDIF
! lightening density
IF (LFLASH) THEN
  IF (YDCPG_OPTS%NSTEP==0) YDMF_PHYS%OUT%FLASH=0._JPRB

  CALL DIAGFLASH(YDCFU,YDMODEL%YRML_PHY_MF,YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDCPG_OPTS%NSTEP,&
    &ZQCM,ZQIM,ZQRM,ZQSM,ZQGM,ZQHM,YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,ZTM,YDMF_PHYS_BASE_STATE%YCPG_PHY%W,ZDUMMY,ZDUMMY,ZDUMMY,ZDUMMY,&
    &YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS%OUT%FLASH)
ENDIF
!!! modif pour LMSE non activee
IF (LLMSE) THEN
  DO JLEV=1,NTSSG+1
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%FCS(JLON,JLEV)=-ZSFTH_(JLON)*ZRHODREFM__(JLON,IKB)*RCPD
      YDMF_PHYS%OUT%FCLL(JLON,JLEV) = YDMF_PHYS%OUT%FCLL(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
      YDMF_PHYS%OUT%FCLN(JLON,JLEV) = YDMF_PHYS%OUT%FCLN(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
      YDMF_PHYS%OUT%FEVL(JLON,JLEV) = YDMF_PHYS%OUT%FEVL(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
      YDMF_PHYS%OUT%FEVN(JLON,JLEV) = YDMF_PHYS%OUT%FEVN(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
    ENDDO
  ENDDO
ENDIF
IF (LSFORCS) THEN
  DO JLEV=1,NTSSG+1
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%FCS(JLON,JLEV)=-ZSFTH_(JLON)*ZRHODREFM__(JLON,IKB)*RCPD
      ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZTSURF(JLON)))
      YDMF_PHYS%OUT%FCLL(JLON,JLEV)=-ZSFRV_(JLON)*ZRHODREFM__(JLON,IKB)* FOLH (ZTSURF(JLON),0._JPRB)*(1.0_JPRB-ZDELTA)
      YDMF_PHYS%OUT%FCLN(JLON,JLEV)=-ZSFRV_(JLON)*ZRHODREFM__(JLON,IKB)* FOLH (ZTSURF(JLON),0._JPRB)*ZDELTA
    ENDDO
  ENDDO
ENDIF

DO JSG  = 1, NTSSG+1
  DO JLEV = 0, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)=YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)+ZBUDTH_(JLON)
    ENDDO
  ENDDO
ENDDO

IF (LFLEXDIA) THEN
      !   surface variables
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS_SURF%GSD_VF%PLSM,'SVLSM',YDDDH,CDTYPE='V')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTSURF,'SVTS',YDDDH,CDTYPE='V')
      !am: FIXME: issue of shape of PF_T1 (2D) when arg should be 1D
      !CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS_SURF%GSP_SG%PF_T1,'SVWN',YDDDH,CDTYPE='V')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%TCLS,'SVTCLS',YDDDH,CDTYPE='V')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%QCLS,'SVQCLS',YDDDH,CDTYPE='V')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%RHCLS,'SVHUCLS',YDDDH,CDTYPE='V')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%UCLS,'SVUCLS',YDDDH,CDTYPE='V')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%VCLS,'SVVCLS',YDDDH,CDTYPE='V')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%CLPH,'SVPBLH',YDDDH,CDTYPE='V')
      ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*ZWS2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SVWS',YDDDH,CDTYPE='V')
      ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*ZWP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SVWP',YDDDH,CDTYPE='V')
      ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*ZWSI2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SVWIS',YDDDH,CDTYPE='V')
      ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*ZWPI2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SVWIP',YDDDH,CDTYPE='V')
      ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*ZTP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SVTP',YDDDH,CDTYPE='V')
      !am: FIXME: issue of shape of PF_T1 (2D) when arg should be 1D
      !ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*YDMF_PHYS_SURF%GSP_SG%PF_T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDVARS%GEOMETRY%OROG%T0,'SVOROG',YDDDH,CDTYPE='V')
      ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*YDMF_PHYS%OUT%FEVL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SFEVAPLIQ',YDDDH,CDTYPE='F')
      ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*YDMF_PHYS%OUT%FEVN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SFEVAPNEG',YDDDH,CDTYPE='F')
      ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*YDMF_PHYS%OUT%FCLL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SFCHLATLI',YDDDH,CDTYPE='F')
      ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*YDMF_PHYS%OUT%FCLN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SFCHLATNE',YDDDH,CDTYPE='F')
      ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*YDMF_PHYS%OUT%FCS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SFCHSENS',YDDDH,CDTYPE='F')
      !set to 0._JPRB
      ! ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0._JPRB
      !CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SFFONTESL',YDDDH,CDTYPE='F')
      !CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SFLIQSNPL',YDDDH,CDTYPE='F')
      !CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SFFONTESN',YDDDH,CDTYPE='F')

      ! surface radiation
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%FRSO(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,1),'SFRAYSO',YDDDH,CDTYPE='F')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%FRTH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,1),'SFRAYTH',YDDDH,CDTYPE='F')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%FRSODS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),'SFRAYSODS',YDDDH,CDTYPE='F')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%FRTHDS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA),'SFRAYTHDS',YDDDH,CDTYPE='F')
      ZTMP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=1_JPRB-YDMF_PHYS%OUT%FRSO(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG,1)/MAX(ZEPSNEB,YDMF_PHYS%OUT%FRSODS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA))
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP,'SVALB',YDDDH,CDTYPE='V')
      ! surface precipitations
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%FPLSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),'SFPRELIGE',YDDDH,CDTYPE='F')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%FPLSN(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),'SFPRENEGE',YDDDH,CDTYPE='F')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%FPLSG(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),'SFPREGRPL',YDDDH,CDTYPE='F')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%FPLSH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),'SFPREHAIL',YDDDH,CDTYPE='F')
      ! surface wind stress
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%STRTU(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),'SFUUTUR',YDDDH,CDTYPE='F')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%STRTV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG),'SFVVTUR',YDDDH,CDTYPE='F')

     ! WRITE(NULOUT,*) 'LFLEXDIA ARPEGE WITH NTOTSURF = ',NTOTSURF,&
     !    & ' AND NTOTSVAR = ',NTOTSVAR, ' AND NTOTSVFS = ',NTOTSVFS
     ! 3D Variables :
      IF (LINTFLEX) THEN
        ZTMP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,0)=0._JPRB
        DO JLEV=1,YDCPG_OPTS%KFLEVG
          ZTMP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZPFPR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,2)+ZPFPR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,4)
        ENDDO
        CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP2(:,:),'FQTPRECISTL',YDDDH)
        DO JLEV=1,YDCPG_OPTS%KFLEVG                                                                                                   
          ZTMP2(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=ZPFPR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,4)+ZPFPR_(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV,5)
        ENDDO
        CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMP2(:,:),'FQTPRECISTN',YDDDH)
      ENDIF

ENDIF


! daand: radflex
IF (LINTFLEX) THEN
  ! account for radiation separately
  LLRAD=.NOT.LRADFLEX
    
  CALL APL_AROME2INTFLEX(YGFL, YDPARAR, YDPHY, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, &
  & YDCPG_OPTS%ZDTPHY, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V,                   &
  & YDMF_PHYS_BASE_STATE%T, YDCPG_GPAR%VTS, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZFPR, LLRAD, YDMF_PHYS%OUT%FRTH,     &
  & YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDT, ZTENDRA, ZTENDTKE,                           &
  & ZTENDEXT, YLPROCSET)
ENDIF


! Precipitation Type

! Compute wet-bulb temperature at 2 meters (suppose homogeneity of qv/ql/qi )
!ZPCLS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)-2._JPRB/ZZZF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1,YDCPG_OPTS%KFLEVG)*&
!                 &(YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG)-YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,YDCPG_OPTS%KFLEVG))

CALL PPWETPOINT(YDCST, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDMF_PHYS_BASE_STATE%YCPG_PHY%PRE(:, YDCPG_OPTS%KFLEVG), &
& YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%QCLS, ZQCM(:, YDCPG_OPTS%KFLEVG), ZQIM(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS%OUT%TPWCLS               &
&                                                             )

IF (LDPRECIPS.OR.LDPRECIPS2) THEN

  !initialisation de ZDZZ
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZDZZ(JLON,1)=ZAPHIM(JLON,0)*ZINVG-ZZZ_(JLON,1)
  ENDDO
  DO JLEV = 2, YDCPG_OPTS%KFLEVG
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZDZZ(JLON,JLEV)=ZZZ_(JLON,JLEV+IKL)-ZZZ_(JLON,JLEV)
    ENDDO
  ENDDO


  ! Compute wet-bulb temperature
  DO JLEV=1,YDCPG_OPTS%KFLEVG
      CALL PPWETPOINT(YDCST, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(:, JLEV), &
      & ZTM(:, JLEV), ZQVM(:, JLEV), ZQCM(:, JLEV), ZQIM(:, JLEV), ZTPW(:, JLEV))
  ENDDO

  IF (LDPRECIPS) THEN
   ! Defined precipitation type 
   !
   ZPRC_DPRECIPS(:,YDCPG_OPTS%NDTPRECCUR)=0._JPRB

   CALL DPRECIPS(YDCST, YDPRECIPS, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDVARS%GEOMETRY%OROG%T0, &
   & YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZDZZ, ZTPW, ZQCM,                                       &
   & YDMF_PHYS%OUT%FPLSL(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS%OUT%FPLSN(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS%OUT%FPLSG(:, YDCPG_OPTS%KFLEVG),       &
   & ZPRC_DPRECIPS(:, YDCPG_OPTS%NDTPRECCUR)                                                                                                &
   &                                                                                             )
  ENDIF

  IF (LDPRECIPS2) THEN

   !Idem for an other time step and an other period
   ZPRC_DPRECIPS2(:,YDCPG_OPTS%NDTPRECCUR2)=0._JPRB

   CALL DPRECIPS(YDCST, YDPRECIPS, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDVARS%GEOMETRY%OROG%T0,  &
   & YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZDZZ, ZTPW, ZQCM,                                        &
   & YDMF_PHYS%OUT%FPLSL(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS%OUT%FPLSN(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS%OUT%FPLSG(:, YDCPG_OPTS%KFLEVG),        &
   & ZPRC_DPRECIPS2(:, YDCPG_OPTS%NDTPRECCUR2)                                                                                               &
   &                                                                                 )

  ENDIF

ENDIF

!Save surface temperature
IF (LMSE.OR.LSFORCS) THEN
  IF (LLXFUMSE) THEN
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS_SURF%GSP_RR%PT_T0(JROF)=YDCPG_GPAR%VTS(JROF)
    ENDDO
  ELSE
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS_SURF%GSP_RR%PT_T1(JROF)=YDCPG_GPAR%VTS(JROF)
    ENDDO
  ENDIF 
ENDIF  
!      4.2  COMPUTE THE PHYS. TENDENCY FOR "T" AND "w"
!           ------------------------------------------

IF (LVERTFE.AND.LVFE_GWMPA) THEN
  ! * case LVFE_GWMPA not yet coded.
  !   (in this case ZGWT1 must be computed at full levels and
  !   not at half levels)
  CALL ABOR1(' APL_AROME: case LVFE_GWMPA not yet coded if LMPA=T!')
ENDIF

! * compute ZTT1:
IF (LSLAG.AND.YDDYNA%LTWOTL) THEN
  DO JLEV=1,YDCPG_OPTS%KFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZTT1(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)+YDCPG_OPTS%ZDTPHY*ZTENDT(JROF,JLEV)
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,YDCPG_OPTS%KFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZTT1(JROF,JLEV)=YDVARS%T%T9(JROF,JLEV)+YDCPG_OPTS%ZDTPHY*ZTENDT(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

! * compute ZGWT1 = tendency of gw:
IF (YDDYNA%LNHDYN) THEN
  ! Valid for LVFE_GWMPA=F only; ZGWT1 assumed to be half level values.
  DO JLEV=1,YDCPG_OPTS%KFLEVG-1
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZGWT1(JROF,JLEV)=0.5_JPRB*RG*(ZTENDW(JROF,JLEV)+ZTENDW(JROF,JLEV+1))
    ENDDO
  ENDDO
  DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZGWT1(JROF,YDCPG_OPTS%KFLEVG)=0.0_JPRB
    ZGWT1(JROF,0)=0.0_JPRB
  ENDDO
ENDIF

! * convert gw tendency in d tendency:
IF(YDDYNA%LNHDYN) THEN

  IF (LGWADV) THEN
    ZTENDD(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=ZGWT1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
  ELSE

    ! * Provide the appropriate version of (RT) at t+dt for GNHGW2SVDAROME:
    IF (L_RDRY_VD) THEN
      ! Use Rd because "dver" is currently defined with Rd.
      ZRTT1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=RD*ZTT1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)
    ELSE
      ! Use "moist R" because "dver" is defined with "moist R".
      ! Unfortunately, R(t+dt) is not yet available there, use R(t) instead.
      ! "Moist R" tendency is neglected in the below call to GNHGW2SVDAROME.
      DO JLEV=1,YDCPG_OPTS%KFLEVG
        DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZRTT1(JROF,JLEV)=YDCPG_DYN0%RCP%R(JROF,JLEV)*ZTT1(JROF,JLEV)
        ENDDO
      ENDDO
    ENDIF

    ! * Do conversion:
    CALL GNHGW2SVDAROME(YDGEOMETRY, YDDYNA%LNHEE, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
                      & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZRTT1, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF, ZGWT1, ZTENDD)  

  ENDIF
ELSE
  ZTENDD=0.0_JPRB
ENDIF

!      4.3  PUT THE TENDENCIES IN PB1/GFLT1/GMVT1.
!           --------------------------------------


IF (LSLAG) CALL CP_PTRSLB1(YDMODEL%YRML_DYN%YRDYN, YDPTRSLB1, ISLB1U9, ISLB1V9, ISLB1T9, ISLB1VD9, &
           & ISLB1GFL9)
IF ( LINTFLEX ) THEN

  ! Set GFL tendencies to 0
  ZTENDGFL(:,:,:) = 0.0_JPRB

  CALL CPTEND_FLEX( YDLDDH, YDMDDH, YGFL, YDPHY, YDDYNA, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,                &
  & YDCPG_OPTS%KFLEVG, YDVARS%GEOMETRY%GNORDL%T0, YDVARS%GEOMETRY%GNORDM%T0, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
  & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%U,           &
  & YDMF_PHYS_BASE_STATE%V, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%YGSP_RR%T, PGFL, YLPROCSET,                 &
  & YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDH, ZTENDGFL, YDMF_PHYS%OUT%FHSCL, YDMF_PHYS%OUT%FHSCN,            &
  & YDMF_PHYS%OUT%FHSSL, YDMF_PHYS%OUT%FHSSN, YDMF_PHYS%OUT%FHPCL, YDMF_PHYS%OUT%FHPCN, YDMF_PHYS%OUT%FHPSL,         &
  & YDMF_PHYS%OUT%FHPSN, PFEPFP =YDMF_PHYS%OUT%FEPFP, PFCMPCQ=YDMF_PHYS%OUT%FCMPCQ, PFCMPSN=YDMF_PHYS%OUT%FCMPSN,    &
  & PFCMPSL=YDMF_PHYS%OUT%FCMPSL, YDDDH=YDDDH )
  
  CALL CPUTQY(YDGEOMETRY%YRDIMV, YDGMV, YGFL, YDPTRSLB1, YDPHY, YDDYNA, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,  &
  & YDCPG_OPTS%KFLEVG, YDCPG_OPTS%ZDTPHY, IPGFL, ISLB1T9, ISLB1U9, ISLB1V9, ISLB1VD9, ISLB1GFL9, ZTENDH, ZTENDT,                 &
  & YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDD, ZTENDGFL, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, &
  & YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V,   &
  & YDCPG_SL1%ZVIEW, PGMVT1, PGFLT1, YDMF_PHYS%OUT%FDIS)    
   
ELSE

  ! start ZTENDGFLR at 1 because it is dimensionned (:,:,0:n)
  CALL CPUTQY_AROME_EXPL (YDMF_PHYS_NEXT_STATE, YDVARS, YDMODEL, YDGEOMETRY%YRDIMV, YDCPG_BNDS,  &
  & YDCPG_OPTS, YDCPG_OPTS%ZDTPHY, ZTEND_Q, ZTEND_L, ZTEND_R, ZTEND_I, ZTEND_S, ZTEND_G, ZTEND_H, ZTENDTKE, &
  & ZTENDT, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDD)
  CALL CPUTQY_AROME_LOOP (YDMODEL, YDGEOMETRY%YRDIMV, YDGMV, YDCPG_BNDS, YDCPG_OPTS, YDCPG_OPTS%ZDTPHY, IPGFL, &
  & IPTR, ZTENDGFLR(:, :, 1:), YDCPG_SL1%ZVIEW, PGMVT1, PGFLT1)
ENDIF
  

!     ------------------------------------------------------------------ 
!*       5.    Final calculations.
!              -------------------

! * Restore the initial value of some pseudo-historical surface buffers
!   if relevant.
IF (LLDIAB) THEN
  IF (LL_SAVE_PHSURF) THEN
    CALL MF_PHYS_SAVE_PHSURF_PART2 (YDCPG_BNDS, YDCPG_OPTS, ZSAV_DDAL, ZSAV_DDOM, ZSAV_ENTCH,                           &
    & ZSAV_FHPS, ZSAV_GZ0F, ZSAV_GZ0HF, ZSAV_HV, ZSAV_PBLH, ZSAV_QSH, ZSAV_UDAL, ZSAV_UDGRO,                            &
    & ZSAV_UDOM, ZSAV_UNEBH, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VH%PPBLH, YDMF_PHYS_SURF%GSD_VH%PQSH,       &
    & YDMF_PHYS_SURF%GSD_VH%PSPSH, YDMF_PHYS_SURF%GSD_VK%PUDGRO, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VV%PZ0H, &
    & YDVARS%DAL%T0, YDVARS%DOM%T0, YDVARS%UAL%T0, YDVARS%UEN%T0, YDVARS%UNEBH%T0, YDVARS%UOM%T0,                       &
    & YDMODEL)
  ENDIF
ENDIF

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
  & ZRDG_CVGQ, ZRDG_LCVQ, ZRDG_MU0               )
ENDIF

IF (LEDR) THEN
  YDMF_PHYS_SURF%GSD_DI%PXEDR(:,:)=YDMF_PHYS%OUT%EDR(:,:)
ENDIF

CALL MF_PHYS_PRECIPS (YDCPG_BNDS, YDCPG_OPTS, ZPRC_DPRECIPS, ZPRC_DPRECIPS2, YDMF_PHYS_SURF%GSD_XP%PPRECIP, &
& YDMF_PHYS_SURF%GSD_XP2%PPRECIP2, YDMODEL)

!       6. destructor for procset
IF (LINTFLEX) CALL CLEANINTPROCSET(YLPROCSET)

! Clear SPP
IF (YDSPP_CONFIG%LSPP) CALL CLEAR_ALL_SPP(ZSPP_ALL)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('APL_AROME',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE SWAP_THS
IF (LLSWAP_THS) THEN
  ZTHSIN_ => ZTHSAVE__(:,1:YDCPG_OPTS%KFLEVG)
  ZTHS__ => ZTHSWAP__
ELSE
  ZTHSIN_ => ZTHSWAP__(:,1:YDCPG_OPTS%KFLEVG)
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
