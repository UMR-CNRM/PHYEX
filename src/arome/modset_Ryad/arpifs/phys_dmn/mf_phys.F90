#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE MF_PHYS(YDGEOMETRY,YDGMV,YDSURF,YDCFU,YDXFU,YDMODEL,&
 !---------------------------------------------------------------------
 ! - INPUT and INOUT.
 & KBL,KGPCOMP,KST,KEND,KGL1,KGL2,KSTGLO,&
 & LDCONFX,PDTPHY,&
 & KIBL,POROGL,POROGM,&
 & PCUCONVCA,PNLCONVCA,&
 & PGMV,PGMVS,PGFL,PWT0,PWT0L,PWT0M,&
 & PRCP0,PHI0,PHIF0,PRE0,PRE0F,PREPHY0,PREPHY0F,PXYB0,&
 & PWT9,PRCP9,PHI9,PHIF9,PRE9,PRE9F,PREPHY9,PREPHY9F,PXYB9,&
 & PKOZO,PGP2DSDT,PGRADH_PHY,&
 & PFORCEU,PFORCEV,PFORCET,PFORCEQ,&
 & PCTY0,&
 & PB1,PB2,PGMVT1,PGMVT1S,PGFLT1,&
 & PSP_SB,PSP_SG,PSP_RR,&
 & PSD_VF,PSD_VP,PSD_VV,PSD_VH,PSD_VK,PSD_VA,PSD_VC,PSD_DI,PSD_VD,PSD_SFL,&
 & PSD_SFO,PSD_XP,PSD_XP2,&
 & PEMTD,PEMTU,PTRSW,PRMOON,PGPAR,&
 & PGDEOSI,PGUEOSI,PGMU0,PGMU0_MIN,PGMU0_MAX,&
 & PGDEOTI,PGDEOTI2,PGUEOTI,PGUEOTI2,PGEOLT,PGEOXT,&
 & PGRPROX,PGMIXP,PGFLUXC,PGRSURF,&
 !---------------------------------------------------------------------
 ! - OUTPUT .
 & PDHSF,PALBDG,PCAPE,PCTOP,PCLCC,PCLCH,PCLCL,PCLCM,PCLCT,PCLPH,PVEIN,PCT,&
 & PDIFCQ,PDIFCQI,PDIFCQL,PDIFCS,PDIFTQ,PDIFTQI,PDIFTQL,PDIFTS,&
 & PFCCQL,PFCCQN,PFCHOZ,PFCHSP,PFCLL,PFCLN,PFCQING,PFCQLNG,PFCQNG,PFCS,&
 & PFCSQL,PFCSQN,PFEVL,PFEVN,PFEVV,PFGEL,PFGELS,PFLWSP,PFONTE,&
 & PFPLCL,PFPLCN,PFPLCG,PFPLCHL,PFPLSL,PFPLSN,PFPLSG,PFPLSHL,&
 & PMRT,PFRMH,PFRSO,PFRSOC,PFRSODS,PFRSOLU,PFRSGNI,&
 & PFRSDNI,PFRSOPS,PFRSOPT,PFRTH,PFRTHC,PFRTHDS,PFTR,PGZ0,PGZ0H,PNEB,&
 & PQCLS,PQICE,PQLI,PQS,&
 & PRH,PRHCLS,PRUISL,PRUISP,PRUISS,&
 & PSTRCU,PSTRCV,PSTRDU,PSTRDV,PSTRMU,PSTRMV,PSTRTU,PSTRTV,&
 & PDIFCQLC,PDIFCQIC,PFIMCC,&
 & PFEDQLC,PFEDQIC,PFEDQRC,PFEDQSC,PFCNEGQLC,PFCNEGQIC,PFCNEGQRC,PFCNEGQSC,&
 & PUCLS,PVCLS,PNUCLS,PNVCLS,PTCLS,PUGST,PVGST,&
 & PMOCON,PFDIS,&
 & PFHPCL,PFHPCN,PFHPCG,PFHPSL,PFHPSN,PFHPSG,PFHSCL,PFHSCN,PFHSSL,PFHSSN,PFHSSG,& 
 & PFEPFP, PFCMPCQ, PFCMPSN, PFCMPSL,&
 & PTENDU,PTENDV,&
 & PQSOL,&
 & PFPFPSL,PFPFPSN,PFPFPSG,PFPFPCL,PFPFPCN,PFPEVPSL,PFPEVPSN,PFPEVPSG,&
 & PFPEVPCL,PFPEVPCN,PFPEVPCG,PDERNSHF,PFCQRNG, PFCQSNG,PFCQGNG,&
 & PDIAGH,PTPWCLS,PVISICLD, PVISIHYDRO,PMXCLWC,&
 & PFLASH,PTRAJ_PHYS,YDDDH,PFTCNS)

!**** *MF_PHYS* METEO-FRANCE PHYSICS.

!     Purpose.
!     --------
!         Call METEO-FRANCE physics and physical tendencies.

!**   Interface.
!     ----------
!        *CALL* *MF_PHYS(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KBL       : NPROMA-packets number
!        KGPCOMP   : total number of grid points in the domain
!        KST       : first element of work.
!        KEND      : last element of work.
!        KGL1,KGL2 : first and last latitude of computations.
!                    - bounds NGPTOT-packets for DM calculations.
!        KSTGLO    : global offset.
!        LDCONFX   : (see in CPG)
!        PDTPHY    : timestep used in the physics.
!        KIBL      : index into YRCSGEOM/YRGSGEOM types in YDGEOMETRY
!        POROGL,POROGM: components of grad(orography).
!        PCUCONVCA : CA array for interaction with the physics
!        PNLCONVCA : CA array for interaction with the physics
!        PGMV      : GMV at time t and t-dt.
!        PGMVS     : GMVS at time t and t-dt.
!        PGFL      : GFL at time t and t-dt.
!        PWT0      : w-wind time t.
!        PWT0L     : zonal derivative of w-wind at time t.
!        PWT0M     : merid derivative of w-wind at time t.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PHI0      : geopotential height at half levels at time t.
!        PHIF0     : geopotential height at full levels at time t.
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE0F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PREPHY0   : input pressure "pre" for AROME at half levels at time t.
!        PREPHY0F  : input pressure "pre" for AROME at full levels at time t.
!        PXYB0     : contains pressure depth, "delta", "alpha" at time t.
!        PRKQVH    : Rasch-Kristjansson scheme - water vapour tendency
!        PRKQCH    : Rasch-Kristjansson scheme - condensates tendency
!        PWT9      : Vertical wind time t-dt.
!        PRCP9     : contains "cp", "R" and "Kap=R/Cp" at t-dt.
!        PHI9      : geopotential height at half levels at time t-dt.
!        PHIF9     : geopotential height at full levels at time t-dt.
!        PRE9      : hydrostatic pressure "prehyd" at half levels at time t-dt.
!        PRE9F     : hydrostatic pressure "prehyd" at full levels at time t-dt.
!        PREPHY9   : input pressure "pre" for AROME at half levels at time t-dt.
!        PREPHY9F  : input pressure "pre" for AROME at full levels at time t-dt.
!        PXYB9     : contains pressure depth, "delta", "alpha" at time t-dt.
!        PKOZO     : fields for photochemistery of ozon.
!        PGP2DSDT  : stochastic physics random pattern.
!        PGRADH_PHY: horizontal gradients for physics

!     INPUT/OUTPUT:
!     -------------
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        PB1       : "SLB1"-buffer, used for interpolations in the SL scheme.
!        PB2       : "SLB2"-buffer.
!        PGFLT1    : GFL t+dt
!        PGPAR     : surface fields for AROME.
!        PGDEOSI   : DESCENDING INCREMENTAL OPTICAL DEPTHS, SOLAR
!        PGUEOSI   : ASCENDING  INCREMENTAL OPTICAL DEPTHS, SOLAR
!        PGMU0     : COSINE OF SOLAR ZENITH ANGLE, APPROXIMATE ACTUAL VALUE
!        PGMU0_MIN : COSINE OF SOLAR ZENITH ANGLE, MIN VALUE
!        PGMU0_MAX : COSINE OF SOLAR ZENITH ANGLE, MAX VALUE
!        PGDEOTI   : descending incremental optical depths, dB/dT(T0) weights
!        PGDEOTI2  : descending incremental optical depths, B weights with
!                    linear T_e correction
!        PGUEOTI   : ascending incremental optical depths, dB/dT(T0) weights
!        PGUEOTI2  : ascending incremental optical depths, B weights with
!                    linear T_e correction
!        PGEOLT    : local optical depths, dB/dT(T0) weights
!        PGEOXT    : maximum optical depths for EBL-EAL, dB/dT(T0) weights
!        PGRPROX   : correction term for adjacent exchanges
!        PGMIXP    : non-statistical weights for bracketing
!        PGFLUXC   : out of bracket part of clearsky EBL, resp. EBL-EAL flux
!        PGRSURF   : corrective ratio for surface cts contribution

!     OUTPUT:
!     -------
!        PDHSF     : distribution of horizontal mean weights used for
!                    simplified radiation scheme.
!        ---------------------- output of aplpar ------------------------------
!        PALBDG      : modele surface shortwave albedo (diagnostic).
!        PCAPE     : CAPE.
!        PCTOP     : top of convective nebulosity (diagnostic).
!        PCLCC     : convective cloud cover (diagnostic).
!        PCLCH     : high cloud cover (diagnostic).
!        PCLCL     : low cloud cover (diagnostic).
!        PCLCM     : medium cloud cover (diagnostic).
!        PCLCT     : total cloud cover (diagnostic).
!        PCLPH     : height (in meters) of the PBL.
!        PVEIN     : ventilation index in the PBL.
!        PCT       : thermical coefficient of soil-vegetation middle.
!        PDIFCQ    : convective flux of specific humidity (not rain/snow).
!        PDIFCQI   : convective flux of solid water (not rain/snow).
!        PDIFCQL   : convective flux of liquid water (not rain/snow).
!        PDIFCS    : convective flux of enthalpy (not rain/snow).
!        PDIFTQ    : turbulent flux (inc. "q" negative) of specific humidity.
!        PDIFTQI   : turbulent flux (inc. "q" negative) of solid water.
!        PDIFTQL   : turbulent flux (inc. "q" negative) of liquid water.
!        PDIFTS    : turbulent flux of enthalpy (or dry static energy).
!        PFCCQL    : convective condensation flux for liquid water.
!        PFCCQN    : convective condensation flux for ice.
!        PFCHOZ    : ozon photo-chemical flux.
!        PFPFPSL   : flux of liquid resol. precipitation: the generation term. 
!        PFPFPSN   : flux of solid resolved precipitation: the generation term. 
!        PFPFPCL   : flux of liquid conv. precipitation: the generation term. 
!        PFPFPCN   : flux of solid conv. precipitation: the generation term. 
!        PFPEVPSL  : resolved precipitation flux due to evaporation.
!        PFPEVPSN  : resolved precipitation flux due to sublimation.
!        PFPEVPCL  : convective precipitation flux due to evaporation.
!        PFPEVPCN  : convective precipitation flux due to sublimation.
!        PFCHSP    : heat flux from surface to deep soil.
!        PFCLL     : latent heat flux over liquid water (or wet soil).
!        PFCLN     : latent heat flux over snow (or ice).
!        PFCQING   : pseudo-flux of ice to correct for "qi"<0.
!        PFCQLNG   : pseudo-flux of liquid water to correct for "ql"<0.
!        PFCQNG    : pseudo-flux of water to correct for Q<0.
!        PFCS      : sensible heat flux at surface level.
!        PFCSQL    : stratiform condensation flux for liquid water.
!        PFCSQN    : stratiform condensation flux for ice.
!        PFEVL     : water vapour flux over liquid water (or wet soil).
!        PFEVN     : water vapour flux over snow (or ice) and frozen soil.
!        PFEVV     : evapotranspiration flux.
!        PFGEL     : freezing flux of soil water.
!        PFGELS    : freezing flux of soil water at surface level.
!        PFLWSP    : water flux from surface to deep soil.
!        PFONTE    : water flux corresponding to surface snow melt.
!        PFPLCL    : convective precipitation as rain.
!        PFPLCN    : convective precipitation as snow.
!        PFPLCG    : convective precipitation as graupel.
!        PFPLCHL   : convective precipitation as hail.
!        PFPLSL    : stratiform precipitation as rain.
!        PFPLSN    : stratiform precipitation as snow.
!        PFPLSG    : stratiform precipitation as graupel.
!        PFPLSHL   : stratiform precipitation as hail.
!        PMRT      : mean radiant temperature.
!        PFRMH     : mesospheric enthalpy flux.
!        PFRSO     : shortwave radiative flux.
!        PFRSOC    : shortwave clear sky radiative flux.
!        PFRSODS   : surface downwards solar flux.
!        PFRSOLU   : downward lunar flux at surface.
!        PFRSGNI   : Global normal irradiance
!        PFRSDNI   : Direct normal irradiance
!        PFRSOPS   : surface parallel solar flux.
!        PFRSOPT   : top parallel solar flux.
!        PFRTH     : longwave radiative flux.
!        PFRTHC    : longwave clear sky radiative flux.
!        PFRTHDS   : surface downwards IR flux.
!        PFTR      : transpiration flux.
!        PGZ0      : g*roughness length (current).
!        PGZ0H     : current g*thermal roughness length (if KVCLIV >=8).
!        PNEB      : fractional cloudiness for radiation.
!        PQCLS     : specific humidity at 2 meters (diagnostic).
!        PQICE     : specific humidity of solid water for radiation.
!        PQLI      : specific humidity of liquid water for radiation.
!        PQS       : specific humidity at surface level.
!        PRH       : relative humidity.
!        PRHCLS    : relative humidity at 2 meters (diagnostic).
!        PRUISL    : run-off flux out the interception water-tank.
!        PRUISP    : run-off flux in soil.
!        PRUISS    : run-off flux at surface level.
!        PSTRCU    : convective flux of momentum "U".
!        PSTRCV    : convective flux of momentum "V".
!        PSTRDU    : gravity wave drag flux "U".
!        PSTRDV    : gravity wave drag flux "V".
!        PSTRMU    : mesospheric flux for "U"-momentum.
!        PSTRMV    : mesospheric flux for "V"-momentum.
!        PSTRTU    : turbulent flux of momentum "U".
!        PSTRTV    : turbulent flux of momentum "V".
!        PDIFCQLC to PFCNEGQSC:
!        PUCLS     : U-component of wind at 10 meters (diagnostic).
!        PVCLS     : V-component of wind at 10 meters (diagnostic).
!        PNUCLS    : U-component of neutral wind at 10 meters (diagnostic).
!        PNVCLS    : V-component of neutral wind at 10 meters (diagnostic).
!        PTCLS     : temperature at 2 meters (diagnostic).
!        PTPWCLS   : wet-bulb temperature at 2 meters (diagnostic)
!        PUGST     : U-component of gusts (diagnostic).
!        PVGST     : V-component of gusts (diagnostic).
!        PDERNSHF  : derivative of the non solar surface with respect to Tsurf
!        ---------------------- end of output of aplpar -----------------------
!        PMOCON    : moisture convergence.
!        PFDIS     : enthalpy flux due to dissipation of kinetic energy. 
!        PFHPCL    : liquid water convective condensation enthalpy flux.
!        PFHPCN    : snow convective condensation enthalpy flux.
!        PFHPSL    : liquid water stratiform condensation enthalpy flux.
!        PFHPSN    : snow stratiform condensation enthalpy flux.
!        PFHSCL    : sensible heat flux due to liquid convective precipitations
!        PFHSCN    : sensible heat flux due to snow convective precipitations.
!        PFHSSL    : sensible heat flux due to liquid stratiform precipitations
!        PFHSSN    : sensible heat flux due to snow stratiform precipitations.
!        PTENDU    : "U"-wind tendency due to physics.
!        PTENDV    : "V"-wind tendency due to physics.
!        PQSOL     : surface specific humidity used in case "delta m=1".
!        PFCQRNG   : pseudo-flux of rain to correct for Q<0
!        PFCQSNG   : pseudo-flux of snow to correct for Q<0
!        PDIAGH    : Add Hail diagnostic PDIAGH (AROME)
!        PFLASH    : Add lightening density (fl/ km2 /s )
!        PVISICLD  : Visibility due to ice and/or water cloud
!        PVISIHYDRO : Vsibility due to precipitations(rain, graupel, snow)
!        PMXCLWC   : Cloud Water Liquid Content at HVISI meters

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!       2000-12-04: F. Bouyssel & J.M. Piriou

!     Modifications.
!     --------------
!       04-Mar-2009 A.Alias : call CPTEND/INITAPLPAR modified to add
!                         Humidity Mesopheric flux (ZFRMQ).
!                     and IVCLIA removed and call to CPNUDG modified as
!                         Nuding mask is now in SD_VF group
!                         call HL_APLPAR modified to add PFCQNG for acdifus
!                         call APL_AROME modified to add Sulfate/Volcano aerosols for radaer
!       K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!       2009-10-15 Y. Bouteloup : Store radiative cloud water and ice in GFL (YIRAD and YLRAD)
!       F. Vana   15-Oct-2009 : NSPLTHOI option
!       K.Yessad (Feb 2010): use YM_RADTC and RFORADTC
!       2010-03-26 Y. Bouteloup : Store radiative cloud water and ice in GFL (AROME case)
!       2010-04-26 Y. Bouteloup : Only one call to cputqy, cputqys and cputqy_arome
!            This need the use of ZTENDGFL as argument of cptend, cptend_new and apl_arome.
!       2010-05-11 F. Bouyssel : Use of PINDX, PINDY
!       2010-05-28 C. Geijo    : Fix error in IPTR array element referencing 
!       2010-06-21 O.Riviere/F. Bouyssel : Fix to have Ts evolving in Fa files with Surfex
!       Dec 2010 A.Alias   : ZMU0N added to call CPPHINP/APLPAR/APL_AROME/HL_APLPAR
!                            CALL to CPNUDG with or with LMSE (A.Voldoire)
!       K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!       L. Bengtsson-Sedlar & F. Vana 18-Feb-2011 : CA scheme for convection
!       F. Vana   22-Feb-2011 : 3D turbulence
!       2011-02-01 M. Mokhtari: Add LMDUST and PEXTT9 and PEXTT0 IN APLPAR
!                             (treatment of the desert aerosols) 
!       2011-03 A.Alias  : new argument to  for sunshine hours YSD_VD%YSUND 
!                      CPNUDG if LMSE=.T. or LMSE=.F. (bugfix)
!                      debug ozone GFL (IPO3) (D. St-Martin)
!                      Humidity Mesopheric flux (ZFRMQ) added in CPTEND_NEW
!       F.Bouyssel (26-03-2011): Fix to have Snow in hist file with surfex
!       2011-06: M. Jerczynski - some cleaning to meet norms
!       E. Bazile 2011-08-26 : Output for MUSC 1D with LFA files with WRITEPHYSIO
!         used previously for extracting profiles from 3D (now also available for AROME).
!       K. Yessad (Dec 2011): use YDOROG, YDGSGEOM and YDCSGEOM.
!       2011-11-21 JF Gueremy : dry convective adjustment (LAJUCV)
!       F. Vana  26-Jan-2012 : historic Qs for TOM's BBC.
!       F.Bouttier Jul 2012: stochastic physics for AROME
!       Z. SASSI  : 07-Mar-2013   INITIALIZING THE WEIGHT VECTORS PDHSF(NPROMA)
!       [DISTRIBUTION OF HORIZONTAL MEANS WEIGHTS]
!       F. Vana  28-Nov-2013 : Redesigned trajectory handling
!       T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!       2013-11, D. Degrauwe: Flexible interface CPTEND_FLEX.
!       2013-11, J. Masek: Passing intermittency arrays for ACRANEB2.
!       K. Yessad (July 2014): Move some variables.
!       2016-04, J. Masek: Passing sunshine duration to APL_AROME.
!       2016-09, M. Mokhtari & A. Ambar: replacement of ZEXT and ZEZDIAG by PGFL
!                                        in aplpar.F90 argument.
!       2016-10, P. Marguinaud : Port to single precision
!       K. Yessad (Dec 2016): Prune obsolete options.
!       K. Yessad (June 2017): Introduce NHQE model.
!       2017-09, J. Masek: Shifted dimensioning of PGMU0.
!       K. Yessad (Feb 2018): remove deep-layer formulations.
!       K. Yessad (Apr 2018): introduce key L_RDRY_VD (ensure consistent definition of "dver" everywhere).
!       2018-09, F. Duruisseau: add rconv and sconv in gfl for bayrad
!       2018-09, R. Brozkova: Passing of diagnostic hail, global normal
!         irradiance and mean radiant temperature from APLPAR.
!       2018-09, D. St-Martin : add NOGWD inputs in aplpar
!       2018-09, M. Michou : add ARPEGE-Climat chemistry call in aplpar  
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!       2019-05, I. Etchevers : add visibilities and precipitation type
!   R. El Khatib 27-02-2019 memory bandwidth savings.
!   R. El Khatib 30-Oct-2018 IMAXDRAFT
!       2019-09, M. Hrastinski: Dataflow for TKE and TTE terms in ALARO DDH (PFTCNS).
!       2019-09, J. Masek: Modified call to APL_AROME (added argument NFRRC).
!       2019-12, Y. Bouteloup: Introduction of ZTENDU and ZTENDV for computation of ZDEC in cputqy
!                diferent from PTENDU and PTENDV in the case of the use of Tiedtke scheme to avoid double counting
!       2020-12, U. Andrae : Introduce SPP for HARMONIE-AROME
!       2021-01, R. Brozkova: ALARO graupel fix.
!     R. El Khatib 24-Aug-2021 Fix potentially non-associated pointers
! End Modifications
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGMV             , ONLY : TGMV
USE YOMCFU             , ONLY : TCFU
USE YOMXFU             , ONLY : TXFU
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK
USE SC2PRG_MOD         , ONLY : SC2PRG

USE YOMCT0             , ONLY : LSLAG, LTWOTL, LNHDYN, LAROME, LSFORCS, LNHQE,LCORWAT
USE YOMCT3             , ONLY : NSTEP 
USE YOMCVER            , ONLY : LVERTFE  ,LVFE_GWMPA 
USE YOMDYNA            , ONLY : LGWADV, L3DTURB, L_RDRY_VD
USE YOMNUD             , ONLY : NFNUDG   ,LNUDG 
USE YOMSNU             , ONLY : XPNUDG
USE MODULE_RADTC_MIX   , ONLY : YM_RADTC
USE YOMSCM             , ONLY : LGSCM
USE YOMCST             , ONLY : RG, RD
USE YOMCHET            , ONLY : GCHETN
USE INTDYN_MOD         , ONLY : YYTCTY0  ,YYTRCP0  ,YYTRCP9  ,YYTXYB0_PHY,YYTXYB9_PHY
USE YOMLSFORC          , ONLY : LMUSCLFA
USE YOMSPSDT           , ONLY : YSPPT
USE YOMTRAJ            , ONLY : TRAJ_PHYS_TYPE, LPRTTRAJ
USE YOMLUN             , ONLY : NULOUT

USE DDH_MIX            , ONLY : TYP_DDH
USE INTFLEX_MOD        , ONLY : LINTFLEX, TYPE_INTPROCSET, NEWINTPROCSET, CLEANINTPROCSET
USE SPP_MOD , ONLY : YSPP_CONFIG,YSPP
!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TCFU)        ,INTENT(INOUT) :: YDCFU
TYPE(TXFU)        ,INTENT(INOUT) :: YDXFU
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KBL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPCOMP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGL1 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGL2 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO 
LOGICAL           ,INTENT(IN)    :: LDCONFX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTPHY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCUCONVCA(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNLCONVCA(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRCP0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP0%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHI0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIF0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE0F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREPHY0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREPHY0F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXYB0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB0_PHY%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRCP9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP9%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHI9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIF9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE9F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREPHY9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREPHY9F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXYB9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB9_PHY%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKOZO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG*YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSDT(YDGEOMETRY%YRDIM%NPROMA,YSPPT%YGPSDT(1)%NG2D)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORCEU(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORCEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORCET(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORCEQ(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGRADH_PHY(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_PHY_MF%YRARPHY%NGRADIENTS,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCTY0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY0%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSP_SB(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSP_SBD%NLEVS,YDSURF%YSP_SBD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSP_SG(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSP_SGD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSP_RR(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSP_RRD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VF(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_VFD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VP(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_VPD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VV(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_VVD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VH(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_VHD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VK(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_VKD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VA(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_VAD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VC(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_VCD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_DI(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_DID%NLEVS,YDSURF%YSD_DID%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_XP(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_XPD%NLEVS,YDSURF%YSD_XPD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_XP2(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_XP2D%NLEVS,YDSURF%YSD_XP2D%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VD(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_VDD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_SFL(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_SFLD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_SFO(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_SFOD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEMTD(YDGEOMETRY%YRDIM%NPROMA,1:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEMTU(YDGEOMETRY%YRDIM%NPROMA,1:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTRSW(YDGEOMETRY%YRDIM%NPROMA,1:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRMOON(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGPAR(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_PHY_MF%YRPARAR%NGPAR+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDHSF(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALBDG(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCAPE(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCTOP(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCC(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCH(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCL(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCM(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCT(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLPH(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVEIN(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCT(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFCQ(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFCQI(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFCQL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFCS(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFTQ(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFTQI(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFTQL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFTS(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCCQL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCCQN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCHOZ(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPFPSL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPFPSN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPFPSG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPFPCL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPFPCN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPEVPSL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPEVPSN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPEVPSG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPEVPCL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPEVPCN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPEVPCG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCHSP(YDGEOMETRY%YRDIM%NPROMM,YDSURF%YSP_SBD%NLEVS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCLL(YDGEOMETRY%YRDIM%NPROMM,YDMODEL%YRML_PHY_G%YRDPHY%NTSSG+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCLN(YDGEOMETRY%YRDIM%NPROMM,YDMODEL%YRML_PHY_G%YRDPHY%NTSSG+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQING(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCQLNG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCQNG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCS(YDGEOMETRY%YRDIM%NPROMM,YDMODEL%YRML_PHY_G%YRDPHY%NTSSG+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCSQL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCSQN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEVL(YDGEOMETRY%YRDIM%NPROMM,YDMODEL%YRML_PHY_G%YRDPHY%NTSSG+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEVN(YDGEOMETRY%YRDIM%NPROMM,YDMODEL%YRML_PHY_G%YRDPHY%NTSSG+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEVV(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFGEL(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFGELS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFLWSP(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFONTE(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLCL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLCN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLCG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLCHL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLSL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLSN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLSG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLSHL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMRT(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFRMH(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFCQLC(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFCQIC(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFIMCC (YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEDQLC(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEDQIC(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEDQRC(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEDQSC(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCNEGQLC(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCNEGQIC(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCNEGQRC(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFCNEGQSC(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFRSO(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_G%YRDPHY%NTSSG+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSOC(YDGEOMETRY%YRDIM%NPROMM,0:1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSODS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSOLU(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSGNI(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSDNI(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSOPS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRSOPT(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFRTH(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_G%YRDPHY%NTSSG+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRTHC(YDGEOMETRY%YRDIM%NPROMM,0:1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRTHDS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFTR(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGZ0(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGZ0H(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNEB(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQCLS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQICE(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQLI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQS(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRHCLS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRUISL(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRUISP(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRUISS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRCU(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRCV(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRDU(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRDV(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRMU(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRMV(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRTU(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSTRTV(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUCLS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVCLS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNUCLS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNVCLS(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTCLS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUGST(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVGST(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMOCON(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFDIS(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPCL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPCN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPCG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHSCL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHSCN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHSSL(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHSSN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHSSG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFEPFP(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCMPCQ(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCMPSN(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCMPSL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENDU(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENDV(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSOL(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDERNSHF(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQRNG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQSNG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQGNG(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIAGH(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLASH(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGDEOSI(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG,2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGUEOSI(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG,2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMU0(YDGEOMETRY%YRDIM%NPROMM,0:YDMODEL%YRML_PHY_MF%YRPHY%NSORAYFR-1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMU0_MIN(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMU0_MAX(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGDEOTI(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGDEOTI2(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGUEOTI(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGUEOTI2(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGEOLT(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGEOXT(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGRPROX(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMIXP(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLUXC(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGRSURF(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVISICLD(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFTCNS(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,6)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVISIHYDRO(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMXCLWC(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTPWCLS(YDGEOMETRY%YRDIM%NPROMM)

TYPE (TRAJ_PHYS_TYPE), INTENT(INOUT) :: PTRAJ_PHYS
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH

!     ------------------------------------------------------------------
LOGICAL :: LLDIAB
LOGICAL :: LL_SAVE_PHSURF
LOGICAL :: LLXFUMSE

INTEGER(KIND=JPIM) :: IFIELDSS
INTEGER(KIND=JPIM) :: IBLK
INTEGER(KIND=JPIM) :: IPQ,IPO3,ITDIA,IPTREXT,IPTR_CONT,IEFB1,IEFB2,IEFB3
INTEGER(KIND=JPIM) :: IPGFL(YDMODEL%YRML_GCONF%YGFL%NUMFLDS),IPTR(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)

INTEGER(KIND=JPIM) :: INSTEP_DEB,INSTEP_FIN
INTEGER(KIND=JPIM) :: JLEV, JGFL
INTEGER(KIND=JPIM) :: JROF
INTEGER(KIND=JPIM) :: ISLB1U9  ,ISLB1V9  ,ISLB1T9  ,ISLB1GFL9, ISLB1VD9

INTEGER(KIND=JPIM) :: ICLPH(YDGEOMETRY%YRDIM%NPROMA)          ! cf. KCLPH in APLPAR.

!     --- UPPER AIR PHYSICAL TENDENCIES.
REAL(KIND=JPRB) :: ZTENDH(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! Enthalpy tendency.
REAL(KIND=JPRB) :: ZTENDQ(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! Moisture tendency.
REAL(KIND=JPRB) :: ZTENDPTKE(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)  ! Pseudo progn. TKE

! GFL tendencies for APL_AROME (assumes YDMODEL%YRML_GCONF%YGFL%NUMFLDS>=YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
! for now, use Jovi's trick :
REAL(KIND=JPRB) :: ZTENDGFLR(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG,0:YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
REAL(KIND=JPRB) :: ZTENDGFL(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)   ! GFL tendencies

REAL(KIND=JPRB) :: ZUDOM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZUDAL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDDOM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZDDAL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUNEBH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZENTCH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     --- UPPER AIR PHYSICAL TENDENCIES FOR AROME.
!       (the previous one are not used in AROME)
REAL(KIND=JPRB) :: ZTENDT (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)        ! temperature tendency
REAL(KIND=JPRB) :: ZTENDW (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)        ! W  tendency
REAL(KIND=JPRB) :: ZTENDD (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)        ! d  tendency
REAL(KIND=JPRB) :: ZTENDEXT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)      ! GFL EXTRA tendency
REAL(KIND=JPRB) :: ZTENDEXT_DEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)  ! GFL EXTRA tendency
REAL(KIND=JPRB) :: ZDIFEXT(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)     ! Extra-GFL fluxes.

REAL(KIND=JPRB) :: ZTENDU (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! U tendency without deep convection contribution
REAL(KIND=JPRB) :: ZTENDV (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! V tendency without deep convection contribution
!     --- SURFACE AND DEEP RESERVOIR PHYSICAL TENDENCIES.
REAL(KIND=JPRB) :: ZTDTS(YDGEOMETRY%YRDIM%NPROMM)             ! Surface temperature tendency.
REAL(KIND=JPRB) :: ZTDTP(YDGEOMETRY%YRDIM%NPROMM,YDSURF%YSP_SBD%NLEVS)        ! Deep temperature tendency.
REAL(KIND=JPRB) :: ZTDWS(YDGEOMETRY%YRDIM%NPROMM)             ! Surface water res. tendency.
REAL(KIND=JPRB) :: ZTDWP(YDGEOMETRY%YRDIM%NPROMM)             ! Deep water res. tendency.
REAL(KIND=JPRB) :: ZTDWL(YDGEOMETRY%YRDIM%NPROMM)             ! Interception water res. tendency
REAL(KIND=JPRB) :: ZTDSNS(YDGEOMETRY%YRDIM%NPROMM)            ! Snow res. tendency. 
REAL(KIND=JPRB) :: ZTDWPI(YDGEOMETRY%YRDIM%NPROMM)            ! Deep ice res. tendency.
REAL(KIND=JPRB) :: ZTDWSI(YDGEOMETRY%YRDIM%NPROMM)            ! Surface ice res. tendency.
REAL(KIND=JPRB) :: ZTDALBNS(YDGEOMETRY%YRDIM%NPROMA)          ! Snow albedo tendency.
REAL(KIND=JPRB) :: ZTDRHONS(YDGEOMETRY%YRDIM%NPROMA)          ! Snow density tendency.

!     --- FLUXES FROM PARAMETERISATIONS AND TENDENCIES COMP.
REAL(KIND=JPRB) :: ZFP(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)      ! Total rainfall flux.
REAL(KIND=JPRB) :: ZFPLCH(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)   ! convective precipitation flux,
                                             ! local
REAL(KIND=JPRB) :: ZFPLSH(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)   ! stratiform precipitation flux,
                                             ! local
REAL(KIND=JPRB) :: ZFPLSN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)   ! all solid stratiform precipitation flux,
                                             ! local
REAL(KIND=JPRB) :: ZFTKE(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! TKE flux.
REAL(KIND=JPRB) :: ZFTKEI(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! TKE flux.
REAL(KIND=JPRB) :: ZFEFB1(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! EFB1 flux.
REAL(KIND=JPRB) :: ZFEFB2(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! EFB2 flux.
REAL(KIND=JPRB) :: ZFEFB3(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! EFB3 flux.

!     --- Fields of potential use for HIRLAM .
REAL(KIND=JPRB) :: ZCVGQL(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! local array for cloud water tendency 
REAL(KIND=JPRB) :: ZCVGQI(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! local array for cloud ice tendency
REAL(KIND=JPRB) :: ZCVGT(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)      ! local array for temperature tendency

!     --- MISCELLANEOUS PARAMETERISATIONS, 2D ARRAYS ---
REAL(KIND=JPRB) :: ZCVGQ(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)      ! convergence of humidity
                                             ! ("Kuo" condition).
REAL(KIND=JPRB) :: ZFHP(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)     ! Total enthalpy flux
                                             ! + sensible heat flux.
REAL(KIND=JPRB) :: ZLCVQ(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)      ! limited physical contribution to
                                             ! the diagnostic of moisture
                                             ! convergence (remove
                                             ! large scale precip. contrib.).
REAL(KIND=JPRB) :: ZLSCPE(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! cf. PLSCPE in APLPAR.
REAL(KIND=JPRB) :: ZLH(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)        ! cf. PLH in APLPAR.
REAL(KIND=JPRB) :: ZQW(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)        ! cf. PQW in APLPAR.
REAL(KIND=JPRB) :: ZTW(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)        ! cf. PTW in APLPAR.

REAL(KIND=JPRB) :: ZFRMQ(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! cf. MESOSPHERIC humidity flux in APLPAR.

!     --- 1D DIAGNOSTIC FIELDS, SURFACE FLUXES
REAL(KIND=JPRB) :: ZCD(YDGEOMETRY%YRDIM%NPROMM)               ! cf. PCD in APLPAR.
REAL(KIND=JPRB) :: ZCDN(YDGEOMETRY%YRDIM%NPROMM)              ! cf. PCDN in APLPAR.
REAL(KIND=JPRB) :: ZCH(YDGEOMETRY%YRDIM%NPROMM)               ! cf. PCH in APLPAR.
REAL(KIND=JPRB) :: ZEMIS(YDGEOMETRY%YRDIM%NPROMM)             ! cf. PEMIS in APLPAR.
REAL(KIND=JPRB) :: ZFEVI(YDGEOMETRY%YRDIM%NPROMM,YDMODEL%YRML_PHY_G%YRDPHY%NTSSG+1)     ! cf. PFEVI in APLPAR.
REAL(KIND=JPRB) :: ZNEIJ(YDGEOMETRY%YRDIM%NPROMM)             ! cf. PNEIJ in APLPAR.
REAL(KIND=JPRB) :: ZVEG(YDGEOMETRY%YRDIM%NPROMM)              ! cf. PVEG in APLPAR.
REAL(KIND=JPRB) :: ZQSAT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)      ! specific humidity at saturation.
REAL(KIND=JPRB) :: ZQSATS(YDGEOMETRY%YRDIM%NPROMA)            ! cf. PQSATS in APLPAR.
REAL(KIND=JPRB) :: ZQS1(YDGEOMETRY%YRDIM%NPROMM)              ! cf. PQS1 in APLPARS.

!     --- DIAGNOSTIC FIELDS STATE OF SURFACE AIR ---
REAL(KIND=JPRB) :: ZC1(YDGEOMETRY%YRDIM%NPROMM)               ! cf. PC1 in APLPAR.
REAL(KIND=JPRB) :: ZC2(YDGEOMETRY%YRDIM%NPROMM)               ! cf. PC2 in APLPAR.
REAL(KIND=JPRB) :: ZCPS(YDGEOMETRY%YRDIM%NPROMM)              ! cf. PCPS in APLPAR.
REAL(KIND=JPRB) :: ZLHS(YDGEOMETRY%YRDIM%NPROMM)              ! cf. PLHS in APLPAR.
REAL(KIND=JPRB) :: ZRS(YDGEOMETRY%YRDIM%NPROMM)               ! cf. PRS in APLPAR.

!     --- RADIATION COEFFICIENTS FOR SIMPLIFIED PHYSICS IN GRID-POINT ---
REAL(KIND=JPRB) :: ZAC(YDGEOMETRY%YRDIM%NPROMM,(YDGEOMETRY%YRDIMV%NFLEVG+1)*(YDGEOMETRY%YRDIMV%NFLEVG+1))   ! Curtis matrix.
REAL(KIND=JPRB) :: ZAC_HC(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1)           ! horizontally-constant field for ZAC.
REAL(KIND=JPRB) :: ZRADTC(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG,YM_RADTC%NDIM) ! others than Curtis matrix.

!     --- GEOMETRY FOR RADIATION ---
REAL(KIND=JPRB) :: ZMMU0(YDGEOMETRY%YRDIM%NPROMA)  ! mean solar angle (for given day for
                                  ! simpl. rad. scheme)
REAL(KIND=JPRB) :: ZMU0(YDGEOMETRY%YRDIM%NPROMA)   ! local cosine of instantaneous solar zenith
                                  ! angle.
REAL(KIND=JPRB) :: ZMU0LU(YDGEOMETRY%YRDIM%NPROMM) ! local cosine of instantaneous lunar zenith
                                  ! angle.
REAL(KIND=JPRB) :: ZMU0M(YDGEOMETRY%YRDIM%NPROMA)  ! local cosine of averaged solar zenith angle
REAL(KIND=JPRB) :: ZMU0N(YDGEOMETRY%YRDIM%NPROMA)  ! same as ZMU0 for next time step (used for YDMODEL%YRML_PHY_MF%YRARPHY%LMSE)    

!     ---FOR AROME PHYSICS  ---
REAL(KIND=JPRB) :: ZDT   !pour cputqy_arome, a changer peut etre plus tard...
REAL(KIND=JPRB) :: ZGWT1(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)     ! vertical velocity calculated by 
                                              ! cputqy_arome before convertion in d
REAL(KIND=JPRB) :: ZTT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)       ! Temperature at t1

! ZRTT1: appropriate version of R*T at t1 for gnhgw2svd
!  Version of R must be consistent with definition of vertical divergence.
REAL(KIND=JPRB) :: ZRTT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     --- Buffers to save the initial value of   ---
!     --- some pseudo-historical surface buffers ---
REAL(KIND=JPRB) :: ZHV(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZGZ0F(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZGZ0HF(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZPBLH(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZFHPS(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZQSH (YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZUDGRO(YDGEOMETRY%YRDIM%NPROMM)

!     --- FOR BAYRAD ALLSKY FRAMEWORK ---
REAL(KIND=JPRB) :: ZQRCONV(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZQSCONV(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)

!  Horizontal exchange coefficients for 3D turbulence
REAL(KIND=JPRB) :: ZKUROV_H(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZKTROV_H(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!  Empty arrays for 3TL scheme (to be coded later)
!!later REAL(KIND=JPRB) :: ZDIVT9(NPROMA,NFLEVG)
!!later REAL(KIND=JPRB) :: ZUT9L(NPROMA,NFLEVG)
!!later REAL(KIND=JPRB) :: ZVT9L(NPROMA,NFLEVG)
REAL(KIND=JPRB) :: ZWT9L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWT9M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDPRECIPS(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS%NDTPREC)
REAL(KIND=JPRB) :: ZDPRECIPS2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS%NDTPREC2)

REAL(KIND=JPRB) :: ZTAUX(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 &  ZDTAJU(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZEDR(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

! To save Tt for NHQE model
REAL(KIND=JPRB) :: ZTT0_SAVE(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTT0L_SAVE(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTT0M_SAVE(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTT9_SAVE(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)

! required for INTFLEX
TYPE(TYPE_INTPROCSET) :: YLPROCSET

REAL(KIND=JPRB), TARGET :: ZAESUL_NULL(YDGEOMETRY%YRDIM%NPROMA), ZVAVOL_NULL(YDGEOMETRY%YRDIM%NPROMA)
! SPP
REAL(KIND=JPRB) :: ZGP2DSPP(YDGEOMETRY%YRDIM%NPROMA,YSPP%N2D)

INTEGER(KIND=JPIM) :: IMAXDRAFT

INTEGER(KIND=JPIM) :: IPTRLIMA
INTEGER(KIND=JPIM) :: IRR ! pointer of 1st hydrometeors in ZTENDGFLR
INTEGER(KIND=JPIM) :: IPTRTKE ! pointer of TKE in ZTENDGFLR

REAL(KIND=JPRB), POINTER :: ZT0DIV(:,:), ZT0SP(:), ZT0SPD(:,:), ZT0SPDL(:,:)
REAL(KIND=JPRB), POINTER :: ZT0SPDM(:,:), ZT0SPL(:), ZT0SPM(:), ZT0T(:,:)
REAL(KIND=JPRB), POINTER :: ZT0TL(:,:), ZT0TM(:,:), ZT0U(:,:), ZT0UL(:,:)
REAL(KIND=JPRB), POINTER :: ZT0V(:,:), ZT0VL(:,:), ZT0VOR(:,:), ZT1SP(:)
REAL(KIND=JPRB), POINTER :: ZT9DIV(:,:), ZT9SP(:), ZT9SPD(:,:), ZT9T(:,:)
REAL(KIND=JPRB), POINTER :: ZT9U(:,:), ZT9UL(:,:), ZT9V(:,:), ZT9VL(:,:)
REAL(KIND=JPRB), POINTER :: ZT9VOR(:,:)
REAL(KIND=JPRB), POINTER :: ZDIXEDR(:,:), ZSFLSFL1(:), ZSFOSFO1(:), ZVADES(:)
REAL(KIND=JPRB), POINTER :: ZVALAN(:), ZVASEA(:), ZVASOO(:), ZVASUL(:)
REAL(KIND=JPRB), POINTER :: ZVAVOL(:), ZVCVC1(:), ZVDSUND(:), ZVFALBF(:)
REAL(KIND=JPRB), POINTER :: ZVFALBSF(:), ZVFEMISF(:), ZVFGETRL(:), ZVFLSM(:)
REAL(KIND=JPRB), POINTER :: ZVFNUDM(:), ZVFVEG(:), ZVFVRLAN(:), ZVFVRLDI(:)
REAL(KIND=JPRB), POINTER :: ZVFZ0F(:), ZVFZ0RLF(:), ZVHBCCH(:), ZVHPBLH(:)
REAL(KIND=JPRB), POINTER :: ZVHQSH(:), ZVHSCCH(:), ZVHSPSH(:), ZVHTCCH(:)
REAL(KIND=JPRB), POINTER :: ZVKUDGRO(:), ZVPTPC(:), ZVPWPC(:), ZVVALV(:)
REAL(KIND=JPRB), POINTER :: ZVVARG(:), ZVVD2(:), ZVVHV(:), ZVVIVEG(:)
REAL(KIND=JPRB), POINTER :: ZVVLAI(:), ZVVRSMIN(:), ZVVSAB(:), ZVVZ0H(:)
REAL(KIND=JPRB), POINTER :: ZRRFC0(:), ZRRFC1(:), ZRRFC9(:), ZRRIC0(:)
REAL(KIND=JPRB), POINTER :: ZRRIC1(:), ZRRIC9(:), ZRRT0(:), ZRRT1(:)
REAL(KIND=JPRB), POINTER :: ZRRT9(:), ZRRW0(:), ZRRW1(:), ZRRW9(:)
REAL(KIND=JPRB), POINTER :: ZSBQ0(:,:), ZSBQ1(:,:), ZSBQ9(:,:), ZSBTL0(:,:)
REAL(KIND=JPRB), POINTER :: ZSBTL1(:,:), ZSBTL9(:,:), ZSBT0(:,:), ZSBT1(:,:)
REAL(KIND=JPRB), POINTER :: ZSBT9(:,:), ZSGA0(:), ZSGA1(:), ZSGA9(:)
REAL(KIND=JPRB), POINTER :: ZSGF0(:), ZSGF1(:), ZSGF9(:), ZSGR0(:)
REAL(KIND=JPRB), POINTER :: ZSGR1(:), ZSGR9(:), ZSGT1(:)
REAL(KIND=JPRB), POINTER :: ZMCOR(:,:), ZMRAB3C(:,:), ZMRAB3N(:,:), ZMRAB4C(:,:)
REAL(KIND=JPRB), POINTER :: ZMRAB4N(:,:), ZMRAB6C(:,:), ZMRAB6N(:,:), ZMRAT1C(:,:)
REAL(KIND=JPRB), POINTER :: ZMRAT1N(:,:), ZMRAT2C(:,:), ZMRAT2N(:,:), ZMRAT3C(:,:)
REAL(KIND=JPRB), POINTER :: ZMRAT3N(:,:), ZMRAT4C(:,:), ZMRAT4N(:,:), ZMRAT5C(:,:)
REAL(KIND=JPRB), POINTER :: ZMRAT5N(:,:)
REAL(KIND=JPRB), POINTER :: ZPA1(:,:), ZP1CHEM(:,:), ZP1CHEM9(:,:), ZPCPF(:,:)
REAL(KIND=JPRB), POINTER :: ZPCPF1(:,:), ZPCVGQ(:,:), ZPCVGQL(:,:), ZPCVGQM(:,:)
REAL(KIND=JPRB), POINTER :: ZPCVV(:,:), ZPCVV1(:,:), ZPCVV9(:,:), ZPDAL(:,:)
REAL(KIND=JPRB), POINTER :: ZPDAL1(:,:), ZPDOM(:,:), ZPDOM1(:,:), ZPEFB1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDEFB11(:,:), ZPEFB19(:,:), ZPEFB2(:,:), ZPTENDEFB21(:,:)
REAL(KIND=JPRB), POINTER :: ZPEFB29(:,:), ZPEFB3(:,:), ZPTENDEFB31(:,:), ZPEFB39(:,:)
REAL(KIND=JPRB), POINTER :: ZP1EXT(:,:), ZP1EXT9(:,:), ZP1EZDIAG(:,:), ZPFQTUR(:,:)
REAL(KIND=JPRB), POINTER :: ZPFQTUR1(:,:), ZPFSTUR(:,:), ZPFSTUR1(:,:), ZPG(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDG1(:,:), ZPG9(:,:), ZPH(:,:), ZPH9(:,:)
REAL(KIND=JPRB), POINTER :: ZPICONV(:,:), ZPTENDICONV1(:,:), ZPI(:,:), ZPTENDI1(:,:)
REAL(KIND=JPRB), POINTER :: ZPI9(:,:), ZPIRAD1(:,:), ZPLCONV(:,:), ZPTENDLCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZP1LIMA(:,:), ZP1LIMA9(:,:), ZPL(:,:), ZPTENDL1(:,:)
REAL(KIND=JPRB), POINTER :: ZPL9(:,:), ZPLRAD1(:,:), ZPMXL(:,:), ZPMXL1(:,:)
REAL(KIND=JPRB), POINTER :: ZP1NOGW(:,:), ZP1NOGW9(:,:), ZP2NOGW(:,:), ZP2NOGW9(:,:)
REAL(KIND=JPRB), POINTER :: ZPO3(:,:), ZPO31(:,:), ZPO39(:,:), ZPQ(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDQ1(:,:), ZPQ9(:,:), ZPQL(:,:), ZPQM(:,:)
REAL(KIND=JPRB), POINTER :: ZPRCONV(:,:), ZPRCONV1(:,:), ZPTENDRCONV1(:,:), ZPRKTH(:,:)
REAL(KIND=JPRB), POINTER :: ZPRKTH1(:,:), ZPRKTQC(:,:), ZPRKTQC1(:,:), ZPRKTQV(:,:)
REAL(KIND=JPRB), POINTER :: ZPRKTQV1(:,:), ZPR(:,:), ZPTENDR1(:,:), ZPR9(:,:)
REAL(KIND=JPRB), POINTER :: ZPSCONV(:,:), ZPSCONV1(:,:), ZPTENDSCONV1(:,:), ZPSHTUR(:,:)
REAL(KIND=JPRB), POINTER :: ZPSHTUR1(:,:), ZPS(:,:), ZPTENDS1(:,:), ZPS9(:,:)
REAL(KIND=JPRB), POINTER :: ZPSPF(:,:), ZPSPF1(:,:), ZPSRC(:,:), ZPSRC1(:,:)
REAL(KIND=JPRB), POINTER :: ZPSRC9(:,:), ZPTKE(:,:), ZPTENDTKE1(:,:), ZPTKE9(:,:)
REAL(KIND=JPRB), POINTER :: ZPTTE(:,:), ZPTTE1(:,:), ZPUAL(:,:), ZPUAL1(:,:)
REAL(KIND=JPRB), POINTER :: ZPUEN(:,:), ZPUEN1(:,:), ZPUNEBH(:,:), ZPUNEBH1(:,:)
REAL(KIND=JPRB), POINTER :: ZPUOM(:,:), ZPUOM1(:,:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "apl_arome.intfb.h"
#include "aplpar.intfb.h"
#include "aplpar2intflex.intfb.h"
#include "aplpars.intfb.h"
#include "aplpassh.intfb.h"
#include "cpchet.intfb.h"
#include "cpmvvps.intfb.h"
#include "cpnudg.intfb.h"
#include "cpozo.intfb.h"
#include "cpphinp.intfb.h"
#include "cppsolan.intfb.h"
#include "cpqsol.intfb.h"
#include "cptend_new.intfb.h"
#include "cptend_flex.intfb.h"
#include "cptends.intfb.h"
#include "cptendsm.intfb.h"
#include "cputqy_arome.intfb.h"
#include "cputqy.intfb.h"
#include "cputqys.intfb.h"
#include "cpwts.intfb.h"
#include "cp_ptrslb1.intfb.h"
#include "gnhgw2svdarome.intfb.h"
#include "initaplpar.intfb.h"
#include "profilechet.intfb.h"
#include "rdradcoef.intfb.h"
#include "writephysio.intfb.h"
#include "wrphtrajm.intfb.h"
#include "wrradcoef.intfb.h"
#include "acajucv.intfb.h"
#include "gnhqe_conv_tempe.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MF_PHYS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDVAB=>YDGEOMETRY%YRVAB, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL),  YDOROG=>YDGEOMETRY%YROROG(KIBL), &
 & YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1,YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2, &
 & YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH,YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL,YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
 & YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH,YDRCOEF=>YDMODEL%YRML_PHY_RAD%YRRCOEF, &
 & YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY,YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, &
 & YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,  &
 & YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR, YDPRECIPS=>YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS, &
 & YDSTOPH=>YDMODEL%YRML_PHY_STOCH%YRSTOPH)

ASSOCIATE(MVTS=>YDPARAR%MVTS, NRR=>YDPARAR%NRR, NGPAR=>YDPARAR%NGPAR, CMF_UPDRAFT=>YDPARAR%CMF_UPDRAFT,&
 & TSPHY=>YDPHY2%TSPHY, &
 & NPROMA=>YDDIM%NPROMA, NPROMM=>YDDIM%NPROMM, &
 & NTSSG=>YDDPHY%NTSSG, NVCLIS=>YDDPHY%NVCLIS, &
 & YT1=>YDGMV%YT1, YT0=>YDGMV%YT0, YT9=>YDGMV%YT9, &
 & LMDUST=>YDARPHY%LMDUST, LMPA=>YDARPHY%LMPA, LMSE=>YDARPHY%LMSE, LMFSHAL=>YDARPHY%LMFSHAL,&
 & YI=>YGFL%YI, YH=>YGFL%YH, YEZDIAG=>YGFL%YEZDIAG, YL=>YGFL%YL, &
 & YEXT=>YGFL%YEXT, YA=>YGFL%YA, YSRC=>YGFL%YSRC, YSPF=>YGFL%YSPF, &
 & YUEN=>YGFL%YUEN, YG=>YGFL%YG, YCVGQ=>YGFL%YCVGQ, YMXL=>YGFL%YMXL, &
 & YQ=>YGFL%YQ, YCPF=>YGFL%YCPF, YR=>YGFL%YR, YSCONV=>YGFL%YSCONV, YS=>YGFL%YS, &
 & YEFB3=>YGFL%YEFB3, YEFB2=>YGFL%YEFB2, YEFB1=>YGFL%YEFB1, YRKTH=>YGFL%YRKTH, &
 & YDOM=>YGFL%YDOM, YFQTUR=>YGFL%YFQTUR, YFSTUR=>YGFL%YFSTUR, &
 & YUNEBH=>YGFL%YUNEBH, YUOM=>YGFL%YUOM, YCVV=>YGFL%YCVV, YO3=>YGFL%YO3, YNOGW=>YGFL%YNOGW, &
 & LCHEM_ARPCLIM=>YDMODEL%YRML_CHEM%YRCHEM%LCHEM_ARPCLIM,&
 & YCHEM=>YGFL%YCHEM, NGFL_EXT=>YGFL%NGFL_EXT, YRKTQC=>YGFL%YRKTQC, YTKE=>YGFL%YTKE, &
 & YIRAD=>YGFL%YIRAD, YDAL=>YGFL%YDAL, YCOMP=>YGFL%YCOMP, &
 & NGFL_EZDIAG=>YGFL%NGFL_EZDIAG, YRKTQV=>YGFL%YRKTQV, YUAL=>YGFL%YUAL, &
 & YTTE=>YGFL%YTTE, YLCONV=>YGFL%YLCONV, YSHTUR=>YGFL%YSHTUR, &
 & YRCONV=>YGFL%YRCONV, YICONV=>YGFL%YICONV, YLRAD=>YGFL%YLRAD, &
 & YSP_SBD=>YDSURF%YSP_SBD, YSD_SFLD=>YDSURF%YSD_SFLD, YSD_VDD=>YDSURF%YSD_VDD, &
 & YSP_RR=>YDSURF%YSP_RR, YSD_VFD=>YDSURF%YSD_VFD, YSD_VP=>YDSURF%YSD_VP, &
 & YSP_RRD=>YDSURF%YSP_RRD, YSD_VC=>YDSURF%YSD_VC, YSD_VA=>YDSURF%YSD_VA, &
 & YSD_VHD=>YDSURF%YSD_VHD, YSD_VKD=>YDSURF%YSD_VKD, &
 & YSD_VF=>YDSURF%YSD_VF, YSD_VD=>YDSURF%YSD_VD, YSD_VH=>YDSURF%YSD_VH, &
 & YSD_VK=>YDSURF%YSD_VK, YSD_SFOD=>YDSURF%YSD_SFOD, YSP_SG=>YDSURF%YSP_SG, &
 & YSP_SB=>YDSURF%YSP_SB, YSP_SGD=>YDSURF%YSP_SGD, YSD_DI=>YDSURF%YSD_DI, &
 & YSD_VV=>YDSURF%YSD_VV, YSD_VCD=>YDSURF%YSD_VCD, YSD_SFO=>YDSURF%YSD_SFO, &
 & YSD_VVD=>YDSURF%YSD_VVD, YSD_SFL=>YDSURF%YSD_SFL, YSD_VPD=>YDSURF%YSD_VPD, &
 & YSD_VAD=>YDSURF%YSD_VAD, YSD_XP=>YDSURF%YSD_XP,YSD_XP2=>YDSURF%YSD_XP2,&
 & NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, NFLSUL=>YDDIMV%NFLSUL, &
 & LTRAJPS=>YDSIMPHL%LTRAJPS, LSIMPH=>YDSIMPHL%LSIMPH, LRAYSP=>YDSIMPHL%LRAYSP, &
 & LNEBR=>YDPHY%LNEBR, LCDDPRO=>YDPHY%LCDDPRO, LNEBN=>YDPHY%LNEBN, &
 & LMPHYS=>YDPHY%LMPHYS, NSORAYFR=>YDPHY%NSORAYFR, &
 & LCVCSD=>YDPHY%LCVCSD, LCVPRO=>YDPHY%LCVPRO, LNSDO=>YDPHY%LNSDO,&
 & LSTRAPRO=>YDPHY%LSTRAPRO, LUDEVOL=>YDPHY%LUDEVOL, LPTKE=>YDPHY%LPTKE, &
 & NDPSFI=>YDPHY%NDPSFI, LOZONE=>YDPHY%LOZONE, L3MT=>YDPHY%L3MT, &
 & LGPCMT=>YDPHY%LGPCMT, LAJUCV=>YDPHY%LAJUCV, LCVPGY=>YDPHY%LCVPGY, &
 & LRRGUST=>YDPHY%LRRGUST, LEDR=>YDPHY%LEDR, &
 & NTAJUC=>YDTOPH%NTAJUC, NTPLUI=>YDTOPH%NTPLUI, NDIM=>YGFL%NDIM, &
 & NDIM1=>YGFL%NDIM1, NUMFLDS=>YGFL%NUMFLDS, LAGPHY=>YDEPHY%LAGPHY, &
 & LEPHYS=>YDEPHY%LEPHYS, NDIMGMV=>YDGMV%NDIMGMV, &
 & LDPRECIPS=>YDPHY%LDPRECIPS, LDPRECIPS2=>YDPHY%LDPRECIPS2, &
 & NDTPREC=>YDPRECIPS%NDTPREC, NDTPREC2=>YDPRECIPS%NDTPREC2, &
 & NDTPRECCUR=>YDPRECIPS%NDTPRECCUR, NDTPRECCUR2=>YDPRECIPS%NDTPRECCUR2,&
 & NDIMGMVS=>YDGMV%NDIMGMVS, &
 & LSDDH=>YDLDDH%LSDDH, &
 & HDSF=>YDMDDH%HDSF, &
 & MSLB1SP9=>YDPTRSLB1%MSLB1SP9, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB2KAPPAH=>YDPTRSLB2%MSLB2KAPPAH, &
 & MSLB2KAPPAM=>YDPTRSLB2%MSLB2KAPPAM, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & LRCOEF=>YDRCOEF%LRCOEF, &
 & LTLADDIA=>YDRCOEF%LTLADDIA, &
 & NG3SR=>YDRCOEF%NG3SR, &
 & NSTOP=>YDRIP%NSTOP, NLIMA=>YGFL%NLIMA, YLIMA=>YGFL%YLIMA)

CALL SC2PRG(YA%MP1    ,PGFLT1   ,ZPA1)
CALL SC2PRG(1,YCHEM(:)%MP,PGFL  ,ZP1CHEM)
CALL SC2PRG(1,YCHEM(:)%MP9,PGFL ,ZP1CHEM9)
CALL SC2PRG(YCPF%MP   ,PGFL     ,ZPCPF)
CALL SC2PRG(YCPF%MP1  ,PGFLT1   ,ZPCPF1)
CALL SC2PRG(YCVGQ%MP  ,PGFL     ,ZPCVGQ)
CALL SC2PRG(YCVGQ%MPL ,PGFL     ,ZPCVGQL)
CALL SC2PRG(YCVGQ%MPM ,PGFL     ,ZPCVGQM)
CALL SC2PRG(YCVV%MP   ,PGFL     ,ZPCVV)
CALL SC2PRG(YCVV%MP1  ,PGFLT1   ,ZPCVV1)
CALL SC2PRG(YCVV%MP9  ,PGFL     ,ZPCVV9)
CALL SC2PRG(YDAL%MP   ,PGFL     ,ZPDAL)
CALL SC2PRG(YDAL%MP1  ,PGFLT1   ,ZPDAL1)
CALL SC2PRG(YDOM%MP   ,PGFL     ,ZPDOM)
CALL SC2PRG(YDOM%MP1  ,PGFLT1   ,ZPDOM1)
CALL SC2PRG(YEFB1%MP  ,PGFL     ,ZPEFB1)
CALL SC2PRG(YEFB1%MP1 ,ZTENDGFL ,ZPTENDEFB11)
CALL SC2PRG(YEFB1%MP9 ,PGFL     ,ZPEFB19)
CALL SC2PRG(YEFB2%MP  ,PGFL     ,ZPEFB2)
CALL SC2PRG(YEFB2%MP1 ,ZTENDGFL ,ZPTENDEFB21)
CALL SC2PRG(YEFB2%MP9 ,PGFL     ,ZPEFB29)
CALL SC2PRG(YEFB3%MP  ,PGFL     ,ZPEFB3)
CALL SC2PRG(YEFB3%MP1 ,ZTENDGFL ,ZPTENDEFB31)
CALL SC2PRG(YEFB3%MP9 ,PGFL     ,ZPEFB39)
CALL SC2PRG(1,YEXT(:)%MP,PGFL   ,ZP1EXT)
CALL SC2PRG(1,YEXT(:)%MP9,PGFL  ,ZP1EXT9)
CALL SC2PRG(1,YEZDIAG(:)%MP,PGFL,ZP1EZDIAG)
CALL SC2PRG(YFQTUR%MP ,PGFL     ,ZPFQTUR)
CALL SC2PRG(YFQTUR%MP1,PGFLT1   ,ZPFQTUR1)
CALL SC2PRG(YFSTUR%MP ,PGFL     ,ZPFSTUR)
CALL SC2PRG(YFSTUR%MP1,PGFLT1   ,ZPFSTUR1)
CALL SC2PRG(YG%MP     ,PGFL     ,ZPG)
CALL SC2PRG(YG%MP1    ,ZTENDGFL ,ZPTENDG1)
CALL SC2PRG(YG%MP9    ,PGFL     ,ZPG9)
CALL SC2PRG(YH%MP     ,PGFL     ,ZPH)
CALL SC2PRG(YH%MP9    ,PGFL     ,ZPH9)
CALL SC2PRG(YICONV%MP ,PGFL     ,ZPICONV)
CALL SC2PRG(YICONV%MP1,ZTENDGFL ,ZPTENDICONV1)
CALL SC2PRG(YI%MP     ,PGFL     ,ZPI)
CALL SC2PRG(YI%MP1    ,ZTENDGFL ,ZPTENDI1)
CALL SC2PRG(YI%MP9    ,PGFL     ,ZPI9)
CALL SC2PRG(YIRAD%MP1 ,PGFLT1   ,ZPIRAD1)
CALL SC2PRG(YLCONV%MP ,PGFL     ,ZPLCONV)
CALL SC2PRG(YLCONV%MP1,ZTENDGFL ,ZPTENDLCONV1)
CALL SC2PRG(1,YLIMA(:)%MP,PGFL  ,ZP1LIMA)
CALL SC2PRG(1,YLIMA(:)%MP9,PGFL ,ZP1LIMA9)
CALL SC2PRG(YL%MP     ,PGFL     ,ZPL)
CALL SC2PRG(YL%MP1    ,ZTENDGFL ,ZPTENDL1)
CALL SC2PRG(YL%MP9    ,PGFL     ,ZPL9)
CALL SC2PRG(YLRAD%MP1 ,PGFLT1   ,ZPLRAD1)
CALL SC2PRG(YMXL%MP   ,PGFL     ,ZPMXL)
CALL SC2PRG(YMXL%MP1  ,PGFLT1   ,ZPMXL1)
CALL SC2PRG(1,YNOGW(:)%MP,PGFL  ,ZP1NOGW)
CALL SC2PRG(1,YNOGW(:)%MP9,PGFL ,ZP1NOGW9)
CALL SC2PRG(2,YNOGW(:)%MP,PGFL  ,ZP2NOGW)
CALL SC2PRG(2,YNOGW(:)%MP9,PGFL ,ZP2NOGW9)
CALL SC2PRG(YO3%MP    ,PGFL     ,ZPO3)
CALL SC2PRG(YO3%MP1   ,PGFLT1   ,ZPO31)
CALL SC2PRG(YO3%MP9   ,PGFL     ,ZPO39)
CALL SC2PRG(YQ%MP     ,PGFL     ,ZPQ)
CALL SC2PRG(YQ%MP1    ,ZTENDGFL ,ZPTENDQ1)
CALL SC2PRG(YQ%MP9    ,PGFL     ,ZPQ9)
CALL SC2PRG(YQ%MPL    ,PGFL     ,ZPQL)
CALL SC2PRG(YQ%MPM    ,PGFL     ,ZPQM)
CALL SC2PRG(YRCONV%MP ,PGFL     ,ZPRCONV)
CALL SC2PRG(YRCONV%MP1,PGFLT1   ,ZPRCONV1)
CALL SC2PRG(YRCONV%MP1,ZTENDGFL ,ZPTENDRCONV1)
CALL SC2PRG(YRKTH%MP  ,PGFL     ,ZPRKTH)
CALL SC2PRG(YRKTH%MP1 ,PGFLT1   ,ZPRKTH1)
CALL SC2PRG(YRKTQC%MP ,PGFL     ,ZPRKTQC)
CALL SC2PRG(YRKTQC%MP1,PGFLT1   ,ZPRKTQC1)
CALL SC2PRG(YRKTQV%MP ,PGFL     ,ZPRKTQV)
CALL SC2PRG(YRKTQV%MP1,PGFLT1   ,ZPRKTQV1)
CALL SC2PRG(YR%MP     ,PGFL     ,ZPR)
CALL SC2PRG(YR%MP1    ,ZTENDGFL ,ZPTENDR1)
CALL SC2PRG(YR%MP9    ,PGFL     ,ZPR9)
CALL SC2PRG(YSCONV%MP ,PGFL     ,ZPSCONV)
CALL SC2PRG(YSCONV%MP1,PGFLT1   ,ZPSCONV1)
CALL SC2PRG(YSCONV%MP1,ZTENDGFL ,ZPTENDSCONV1)
CALL SC2PRG(YSHTUR%MP ,PGFL     ,ZPSHTUR)
CALL SC2PRG(YSHTUR%MP1,PGFLT1   ,ZPSHTUR1)
CALL SC2PRG(YS%MP     ,PGFL     ,ZPS)
CALL SC2PRG(YS%MP1    ,ZTENDGFL ,ZPTENDS1)
CALL SC2PRG(YS%MP9    ,PGFL     ,ZPS9)
CALL SC2PRG(YSPF%MP   ,PGFL     ,ZPSPF)
CALL SC2PRG(YSPF%MP1  ,PGFLT1   ,ZPSPF1)
CALL SC2PRG(YSRC%MP   ,PGFL     ,ZPSRC)
CALL SC2PRG(YSRC%MP1  ,PGFLT1   ,ZPSRC1)
CALL SC2PRG(YSRC%MP9  ,PGFL     ,ZPSRC9)
CALL SC2PRG(YTKE%MP   ,PGFL     ,ZPTKE)
CALL SC2PRG(YTKE%MP1  ,ZTENDGFL ,ZPTENDTKE1)
CALL SC2PRG(YTKE%MP9  ,PGFL     ,ZPTKE9)
CALL SC2PRG(YTTE%MP   ,PGFL     ,ZPTTE)
CALL SC2PRG(YTTE%MP1  ,PGFLT1   ,ZPTTE1)
CALL SC2PRG(YUAL%MP   ,PGFL     ,ZPUAL)
CALL SC2PRG(YUAL%MP1  ,PGFLT1   ,ZPUAL1)
CALL SC2PRG(YUEN%MP   ,PGFL     ,ZPUEN)
CALL SC2PRG(YUEN%MP1  ,PGFLT1   ,ZPUEN1)
CALL SC2PRG(YUNEBH%MP ,PGFL     ,ZPUNEBH)
CALL SC2PRG(YUNEBH%MP1,PGFLT1   ,ZPUNEBH1)
CALL SC2PRG(YUOM%MP   ,PGFL     ,ZPUOM)
CALL SC2PRG(YUOM%MP1  ,PGFLT1   ,ZPUOM1)
CALL SC2PRG(YM_RADTC%MCOR,ZRADTC   ,ZMCOR)
CALL SC2PRG(YM_RADTC%MRAB3C,ZRADTC ,ZMRAB3C)
CALL SC2PRG(YM_RADTC%MRAB3N,ZRADTC ,ZMRAB3N)
CALL SC2PRG(YM_RADTC%MRAB4C,ZRADTC ,ZMRAB4C)
CALL SC2PRG(YM_RADTC%MRAB4N,ZRADTC ,ZMRAB4N)
CALL SC2PRG(YM_RADTC%MRAB6C,ZRADTC ,ZMRAB6C)
CALL SC2PRG(YM_RADTC%MRAB6N,ZRADTC ,ZMRAB6N)
CALL SC2PRG(YM_RADTC%MRAT1C,ZRADTC ,ZMRAT1C)
CALL SC2PRG(YM_RADTC%MRAT1N,ZRADTC ,ZMRAT1N)
CALL SC2PRG(YM_RADTC%MRAT2C,ZRADTC ,ZMRAT2C)
CALL SC2PRG(YM_RADTC%MRAT2N,ZRADTC ,ZMRAT2N)
CALL SC2PRG(YM_RADTC%MRAT3C,ZRADTC ,ZMRAT3C)
CALL SC2PRG(YM_RADTC%MRAT3N,ZRADTC ,ZMRAT3N)
CALL SC2PRG(YM_RADTC%MRAT4C,ZRADTC ,ZMRAT4C)
CALL SC2PRG(YM_RADTC%MRAT4N,ZRADTC ,ZMRAT4N)
CALL SC2PRG(YM_RADTC%MRAT5C,ZRADTC ,ZMRAT5C)
CALL SC2PRG(YM_RADTC%MRAT5N,ZRADTC ,ZMRAT5N)
CALL SC2PRG(YSD_DI%YXEDR%MP,PSD_DI   ,ZDIXEDR)
CALL SC2PRG(1,YSD_SFL%YSFL(:)%MP,PSD_SFL  ,ZSFLSFL1)
CALL SC2PRG(1,YSD_SFO%YSFO(:)%MP,PSD_SFO  ,ZSFOSFO1)
CALL SC2PRG(YSD_VA%YDES%MP,PSD_VA  ,ZVADES)
CALL SC2PRG(YSD_VA%YLAN%MP,PSD_VA  ,ZVALAN)
CALL SC2PRG(YSD_VA%YSEA%MP,PSD_VA  ,ZVASEA)
CALL SC2PRG(YSD_VA%YSOO%MP,PSD_VA  ,ZVASOO)
CALL SC2PRG(YSD_VA%YSUL%MP,PSD_VA  ,ZVASUL)
CALL SC2PRG(YSD_VA%YVOL%MP,PSD_VA  ,ZVAVOL)
CALL SC2PRG(1,YSD_VC%YVC(:)%MP,PSD_VC,ZVCVC1)
CALL SC2PRG(YSD_VD%YSUND%MP,PSD_VD ,ZVDSUND)
CALL SC2PRG(YSD_VF%YALBF%MP,PSD_VF ,ZVFALBF)
CALL SC2PRG(YSD_VF%YALBSF%MP,PSD_VF,ZVFALBSF)
CALL SC2PRG(YSD_VF%YEMISF%MP,PSD_VF,ZVFEMISF)
CALL SC2PRG(YSD_VF%YGETRL%MP,PSD_VF,ZVFGETRL)
CALL SC2PRG(YSD_VF%YLSM%MP,PSD_VF  ,ZVFLSM)
CALL SC2PRG(YSD_VF%YNUDM%MP,PSD_VF ,ZVFNUDM)
CALL SC2PRG(YSD_VF%YVEG%MP,PSD_VF  ,ZVFVEG)
CALL SC2PRG(YSD_VF%YVRLAN%MP,PSD_VF,ZVFVRLAN)
CALL SC2PRG(YSD_VF%YVRLDI%MP,PSD_VF,ZVFVRLDI)
CALL SC2PRG(YSD_VF%YZ0F%MP,PSD_VF  ,ZVFZ0F)
CALL SC2PRG(YSD_VF%YZ0RLF%MP,PSD_VF,ZVFZ0RLF)
CALL SC2PRG(YSD_VH%YBCCH%MP,PSD_VH ,ZVHBCCH)
CALL SC2PRG(YSD_VH%YPBLH%MP,PSD_VH ,ZVHPBLH)
CALL SC2PRG(YSD_VH%YQSH%MP,PSD_VH  ,ZVHQSH)
CALL SC2PRG(YSD_VH%YSCCH%MP,PSD_VH ,ZVHSCCH)
CALL SC2PRG(YSD_VH%YSPSH%MP,PSD_VH ,ZVHSPSH)
CALL SC2PRG(YSD_VH%YTCCH%MP,PSD_VH ,ZVHTCCH)
CALL SC2PRG(YSD_VK%YUDGRO%MP,PSD_VK,ZVKUDGRO)
CALL SC2PRG(YSD_VP%YTPC%MP,PSD_VP  ,ZVPTPC)
CALL SC2PRG(YSD_VP%YWPC%MP,PSD_VP  ,ZVPWPC)
CALL SC2PRG(YSD_VV%YALV%MP,PSD_VV  ,ZVVALV)
CALL SC2PRG(YSD_VV%YARG%MP,PSD_VV  ,ZVVARG)
CALL SC2PRG(YSD_VV%YD2%MP,PSD_VV   ,ZVVD2)
CALL SC2PRG(YSD_VV%YHV%MP,PSD_VV   ,ZVVHV)
CALL SC2PRG(YSD_VV%YIVEG%MP,PSD_VV ,ZVVIVEG)
CALL SC2PRG(YSD_VV%YLAI%MP,PSD_VV  ,ZVVLAI)
CALL SC2PRG(YSD_VV%YRSMIN%MP,PSD_VV,ZVVRSMIN)
CALL SC2PRG(YSD_VV%YSAB%MP,PSD_VV  ,ZVVSAB)
CALL SC2PRG(YSD_VV%YZ0H%MP,PSD_VV  ,ZVVZ0H)
CALL SC2PRG(YSP_RR%YFC%MP0,PSP_RR  ,ZRRFC0)
CALL SC2PRG(YSP_RR%YFC%MP1,PSP_RR  ,ZRRFC1)
CALL SC2PRG(YSP_RR%YFC%MP9,PSP_RR  ,ZRRFC9)
CALL SC2PRG(YSP_RR%YIC%MP0,PSP_RR  ,ZRRIC0)
CALL SC2PRG(YSP_RR%YIC%MP1,PSP_RR  ,ZRRIC1)
CALL SC2PRG(YSP_RR%YIC%MP9,PSP_RR  ,ZRRIC9)
CALL SC2PRG(YSP_RR%YT%MP0,PSP_RR   ,ZRRT0)
CALL SC2PRG(YSP_RR%YT%MP1,PSP_RR   ,ZRRT1)
CALL SC2PRG(YSP_RR%YT%MP9,PSP_RR   ,ZRRT9)
CALL SC2PRG(YSP_RR%YW%MP0,PSP_RR   ,ZRRW0)
CALL SC2PRG(YSP_RR%YW%MP1,PSP_RR   ,ZRRW1)
CALL SC2PRG(YSP_RR%YW%MP9,PSP_RR   ,ZRRW9)
CALL SC2PRG(YSP_SB%YQ%MP0,PSP_SB   ,ZSBQ0)
CALL SC2PRG(YSP_SB%YQ%MP1,PSP_SB   ,ZSBQ1)
CALL SC2PRG(YSP_SB%YQ%MP9,PSP_SB   ,ZSBQ9)
CALL SC2PRG(YSP_SB%YTL%MP0,PSP_SB  ,ZSBTL0)
CALL SC2PRG(YSP_SB%YTL%MP1,PSP_SB  ,ZSBTL1)
CALL SC2PRG(YSP_SB%YTL%MP9,PSP_SB  ,ZSBTL9)
CALL SC2PRG(YSP_SB%YT%MP0,PSP_SB   ,ZSBT0)
CALL SC2PRG(YSP_SB%YT%MP1,PSP_SB   ,ZSBT1)
CALL SC2PRG(YSP_SB%YT%MP9,PSP_SB   ,ZSBT9)
CALL SC2PRG(YSP_SG%YA%MP0,PSP_SG   ,ZSGA0)
CALL SC2PRG(YSP_SG%YA%MP1,PSP_SG   ,ZSGA1)
CALL SC2PRG(YSP_SG%YA%MP9,PSP_SG   ,ZSGA9)
CALL SC2PRG(YSP_SG%YF%MP0,PSP_SG   ,ZSGF0)
CALL SC2PRG(YSP_SG%YF%MP1,PSP_SG   ,ZSGF1)
CALL SC2PRG(YSP_SG%YF%MP9,PSP_SG   ,ZSGF9)
CALL SC2PRG(YSP_SG%YR%MP0,PSP_SG   ,ZSGR0)
CALL SC2PRG(YSP_SG%YR%MP1,PSP_SG   ,ZSGR1)
CALL SC2PRG(YSP_SG%YR%MP9,PSP_SG   ,ZSGR9)
CALL SC2PRG(YSP_SG%YT%MP1,PSP_SG   ,ZSGT1)
CALL SC2PRG(YT0%MDIV  ,PGMV    ,ZT0DIV)
CALL SC2PRG(YT0%MSP   ,PGMVS   ,ZT0SP)
CALL SC2PRG(YT0%MSPD  ,PGMV    ,ZT0SPD)
CALL SC2PRG(YT0%MSPDL ,PGMV    ,ZT0SPDL)
CALL SC2PRG(YT0%MSPDM ,PGMV    ,ZT0SPDM)
CALL SC2PRG(YT0%MSPL  ,PGMVS   ,ZT0SPL)
CALL SC2PRG(YT0%MSPM  ,PGMVS   ,ZT0SPM)
CALL SC2PRG(YT0%MT    ,PGMV    ,ZT0T)
CALL SC2PRG(YT0%MTL   ,PGMV    ,ZT0TL)
CALL SC2PRG(YT0%MTM   ,PGMV    ,ZT0TM)
CALL SC2PRG(YT0%MU    ,PGMV    ,ZT0U)
CALL SC2PRG(YT0%MUL   ,PGMV    ,ZT0UL)
CALL SC2PRG(YT0%MV    ,PGMV    ,ZT0V)
CALL SC2PRG(YT0%MVL   ,PGMV    ,ZT0VL)
CALL SC2PRG(YT0%MVOR  ,PGMV    ,ZT0VOR)
CALL SC2PRG(YT1%MSP   ,PGMVT1S ,ZT1SP)
CALL SC2PRG(YT9%MDIV  ,PGMV    ,ZT9DIV)
CALL SC2PRG(YT9%MSP   ,PGMVS   ,ZT9SP)
CALL SC2PRG(YT9%MSPD  ,PGMV    ,ZT9SPD)
CALL SC2PRG(YT9%MT    ,PGMV    ,ZT9T)
CALL SC2PRG(YT9%MU    ,PGMV    ,ZT9U)
CALL SC2PRG(YT9%MUL   ,PGMV    ,ZT9UL)
CALL SC2PRG(YT9%MV    ,PGMV    ,ZT9V)
CALL SC2PRG(YT9%MVL   ,PGMV    ,ZT9VL)
CALL SC2PRG(YT9%MVOR  ,PGMV    ,ZT9VOR)
!     ------------------------------------------------------------------

!        0.    constructor for procset
IF (LINTFLEX) YLPROCSET=NEWINTPROCSET()

!        1.    Preliminary calculations necessary
!              for all types of physics.
!              ------------------------------------

IF (LSLAG) CALL CP_PTRSLB1(YDMODEL%YRML_DYN%YRDYN,YDPTRSLB1,ISLB1U9,ISLB1V9,ISLB1T9,ISLB1VD9,ISLB1GFL9)

IBLK=(KSTGLO-1)/NPROMA+1
INSTEP_DEB=1
INSTEP_FIN=1

! initialisation for surfex if XFU
LLXFUMSE=.FALSE.
IF (LDCONFX) THEN
  LLXFUMSE=.TRUE.
ENDIF

! SPP 
IF ( YSPP_CONFIG%LSPP ) THEN
 DO JROF=1,YSPP%N2D
   ZGP2DSPP(:,JROF) = YSPP%GP_ARP(JROF)%GP2D(:,1,KIBL)
 ENDDO
ENDIF

! Complete physics is called.
LLDIAB=(LMPHYS.OR.LEPHYS).AND.(.NOT.LAGPHY)

! In the NHQE model, MF_PHYS enters with Tt and grad(Tt), where Tt = T * exp(-(R/cp) log(pre/prehyd)).
! But calculations of MF_PHYS must use T and grad(T).
! So we do a conversion Tt -> T.
IF (LNHQE) THEN
  ! Valid for NPDVAR=2 only.
  ! At instant t (with the derivatives):
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZTT0_SAVE(JROF,JLEV)=ZT0T(JROF,JLEV)
      ZTT0L_SAVE(JROF,JLEV)=ZT0TL(JROF,JLEV)
      ZTT0M_SAVE(JROF,JLEV)=ZT0TM(JROF,JLEV)
    ENDDO
  ENDDO
  CALL GNHQE_CONV_TEMPE(YDGEOMETRY,.TRUE.,YDMODEL%YRML_GCONF%YGFL%NDIM,KST,KEND,&
   & ZT0SPD,ZT0SP,ZT0T,&
   & KGFLTYP=0,PGFL=PGFL,KDDER=2,PQCHAL=ZT0SPDL,PQCHAM=ZT0SPDM,&
   & PTL=ZT0TL,PTM=ZT0TM)
  ! At instant t-dt for leap-frog advections (without the derivatives):
  IF (.NOT.LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZTT9_SAVE(JROF,JLEV)=ZT9T(JROF,JLEV)
      ENDDO
    ENDDO
    CALL GNHQE_CONV_TEMPE(YDGEOMETRY,.TRUE.,YDMODEL%YRML_GCONF%YGFL%NDIM,KST,KEND,&
     & ZT9SPD,ZT9SP,ZT9T,&
     & KGFLTYP=9,PGFL=PGFL)
  ENDIF
ENDIF

IF (LRAYSP.AND.(NSTEP >= INSTEP_DEB .AND. NSTEP <= INSTEP_FIN)) THEN  
  CALL CPPSOLAN(YDGEOMETRY%YRDIM,KST,KEND,YDGSGEOM%GEMU,YDGSGEOM%GELAM,ZMMU0)
  IF (.NOT.LSDDH) THEN
!-----------INITIALIZING THE WEIGHT VECTORS-------------
DO JROF=KST,KEND
  PDHSF(JROF)=HDSF(JROF+KSTGLO-1)
ENDDO
!-------------------------------------------------------
  ENDIF
ENDIF

IF (LLDIAB.OR.LSIMPH) THEN
  CALL CPPHINP(YDGEOMETRY,YDMODEL,KST,KEND,&
   & YDGSGEOM%GEMU,YDGSGEOM%GELAM,&
   & ZT0U,ZT0V,&
   & ZPQ,ZPQL,ZPQM,ZPCVGQL,ZPCVGQM,&
   & PXYB0(1,1,YYTXYB0_PHY%M_RDELP),PCTY0(1,0,YYTCTY0%M_EVEL),ZPCVGQ,&
   & ZMU0,ZMU0LU,ZMU0M,ZMU0N,ZCVGQ)
  ZLCVQ(KST:KEND,1:NFLEVG)=ZCVGQ(KST:KEND,1:NFLEVG)
ENDIF

DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    ZCVGQL(JROF,JLEV)=0._JPRB
    ZCVGQI(JROF,JLEV)=0._JPRB
    ZCVGT (JROF,JLEV)=0._JPRB
  ENDDO
ENDDO

DO JROF=KST,KEND
  ZQSATS(JROF)=0.0_JPRB
  ZFPLCH(JROF,0)=0.0_JPRB
  ZFPLSH(JROF,0)=0.0_JPRB
ENDDO

IF(LNEBN.OR.LNEBR.OR.LRRGUST) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZFPLCH(JROF,JLEV)=ZPCPF(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF
IF(LRRGUST) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZFPLSH(JROF,JLEV)=ZPSPF(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

! * In some cases, some pseudo-historic surface buffers (like z0) should
!   not be modified between the entrance and the output of MF_PHYS
!   (this is the case for example if LDCONFX=T).
!   For the time being, we must save:
!   - HV (group VV) : resistance to evapotranspiration
!   - Z0F (group VD): gravity * surface roughness length
!   - Z0H (group VV): gravity * roughness length for heat
!   - PBLH (group VH): PBL height
!   - SPSH (group VH):
!   - QSH (group VH):
LL_SAVE_PHSURF=LLDIAB.AND.LDCONFX
IF (LL_SAVE_PHSURF) THEN
  IF(YSD_VV%YHV%LSET) ZHV(1:NPROMA)=ZVVHV(1:NPROMA)
  IF(YSD_VF%YZ0F%LSET) ZGZ0F(1:NPROMA)=ZVFZ0F(1:NPROMA)
  IF(YSD_VV%YZ0H%LSET) ZGZ0HF(1:NPROMA)=ZVVZ0H(1:NPROMA)
  IF(YSD_VH%YPBLH%LSET) ZPBLH(1:NPROMA)=ZVHPBLH(1:NPROMA)
  IF(YSD_VH%YSPSH%LSET) ZFHPS(1:NPROMA)=ZVHSPSH(1:NPROMA)
  IF(YSD_VH%YQSH%LSET) ZQSH(1:NPROMA)=ZVHQSH(1:NPROMA)
  IF(YSD_VK%YUDGRO%LSET) ZUDGRO(1:NPROMA)=ZVKUDGRO(1:NPROMA)
  IF(LCVPRO.OR.LGPCMT) THEN
    ZUDAL(:,:)=ZPUAL(:,:)
    ZUDOM(:,:)=ZPUOM(:,:)
    IF(LCDDPRO) THEN
      ZDDAL(:,:)=ZPDAL(:,:)
      ZDDOM(:,:)=ZPDOM(:,:)
    ENDIF
  ENDIF
  IF(YUNEBH%LACTIVE) ZUNEBH(:,:)=ZPUNEBH(:,:)
  IF(YUEN%LACTIVE)   ZENTCH(:,:)=ZPUEN(:,:)
ENDIF

IF (.NOT.ASSOCIATED(ZVASUL)) THEN
  ZAESUL_NULL(:)=0._JPRB
  ZVASUL=>ZAESUL_NULL
ENDIF
IF (.NOT.ASSOCIATED(ZVAVOL)) THEN
  ZVAVOL_NULL(:)=0._JPRB
  ZVAVOL=>ZVAVOL_NULL
ENDIF

CALL INITAPLPAR ( YGFL, YDARPHY, KST, KEND, NPROMA, NFLEVG, NTSSG, YSP_SBD%NLEVS,&
  & ZVFVEG ,&
  & PDIFCQ , PDIFCQI, PDIFCQL, PDIFCS , ZDIFEXT, PDIFTQ , PDIFTQI, PDIFTQL,&
  & PDIFTS , PFCCQL , PFCCQN , PFCSQL , PFCSQN , PFCQNG , PFCQING,&
  & PFCQLNG, PFCQRNG, PFCQSNG,PFCQGNG,&
  & PFPLCL , PFPLCN , PFPLCG , PFPLCHL, PFPLSL , PFPLSN ,PFPLSG ,&
  & PFPLSHL, PFRSO  , PFRSOC ,&
  & PFRTH  , PFRTHC , PSTRCU , PSTRCV , PSTRDU , PSTRDV , PSTRTU ,&
  & PSTRTV , PSTRMU , PSTRMV , PFRMH  , ZFRMQ  ,&
  & PDIFCQLC,PDIFCQIC,PFIMCC,&
  & PFEDQLC, PFEDQIC, PFEDQRC, PFEDQSC, PFCNEGQLC,PFCNEGQIC,PFCNEGQRC,PFCNEGQSC,&
  & PFCHOZ , ZCPS   , ZLHS   ,&
  & ZRS    , ZLH    , ZLSCPE , PNEB   , PQICE  , PQLI   , ZQSAT  ,&
  & ZQW    , PRH    , ZTW    , PALBDG , ZCD    , ZCDN   , ZCH    ,&
  & ZC1    , ZC2    , PCT    , ZEMIS  , PFCHSP , PFCLL  , PFCLN  ,&
  & PFCS   , ZFEVI  , PFEVL  , PFEVN  , PFEVV  , PFLASH , PFTR   , PFLWSP ,&
  & PFONTE , PFGEL  , PFGELS ,&
  & PFRSGNI, PFRSDNI, PFRSODS, PFRSOPS, PFRSOPT, PFRSOLU, PFRTHDS,&
  & PFPFPSL, PFPFPSN,PFPFPSG, PFPFPCL, PFPFPCN,PFPEVPSL,PFPEVPSN,PFPEVPSG,PFPEVPCL,&
  & PFPEVPCN,PFPEVPCG,ZFTKE  , ZFTKEI , ZFEFB1 , ZFEFB2   , ZFEFB3,&
  & PGZ0   , PGZ0H  , ZNEIJ  , ZVEG   , PQS    , ZQSATS , PRUISL ,&
  & PRUISP , PRUISS , PUCLS  , PVCLS  , PNUCLS , PNVCLS , PTCLS  , PMRT,&
  & PQCLS  , PRHCLS , PCLCT  , PCLCH  , PCLCM  , PCLCL  , PCLCC  ,&
  & PCAPE  , PCTOP  , ICLPH  , PCLPH  , PVEIN  , PUGST  , PVGST  ,&
  & PDIAGH , ZEDR, PVISICLD, PVISIHYDRO, PMXCLWC)

!*       2.    Complete physics.
!              -----------------

!        2.2  Complete physics.
!             -----------------

IF (LLDIAB.AND.(.NOT.LMPA)) THEN

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

  TSPHY = MAX(PDTPHY,1.0_JPRB)

  ! CALL PARAMETERISATIONS

  DO JROF=KST,KEND
    ZVFLSM(JROF)=REAL(NINT(ZVFLSM(JROF)),JPRB)
  ENDDO
  
    IF (LTWOTL) THEN

      IF (LAJUCV) THEN
        DO JLEV=1,NFLEVG
          DO JROF=KST,KEND
            ZTAUX(JROF,JLEV)=ZT0T(JROF,JLEV)
          ENDDO
        ENDDO
        CALL ACAJUCV(YDMODEL%YRML_PHY_MF%YRPHY0,KST,KEND,NPROMA,NTPLUI,NFLEVG,NTAJUC,&
         & PRE0,PXYB0(1,1,YYTXYB0_PHY%M_ALPH),PXYB0(1,1,YYTXYB0_PHY%M_DELP),&
         & PXYB0(1,1,YYTXYB0_PHY%M_LNPR),ZT0T)
        DO JLEV=1,NFLEVG
          DO JROF=KST,KEND
            ZDTAJU(JROF,JLEV)=ZT0T(JROF,JLEV)-ZTAUX(JROF,JLEV)
          ENDDO
        ENDDO
      ENDIF

      ITDIA=1_JPIM
      CALL APLPAR(YDGEOMETRY,YDSURF, YDXFU, YDCFU, YDMODEL, KST    , KEND   , NPROMA ,&
       & ITDIA  , NFLEVG  , KSTGLO,&
       & NVCLIS , YSD_VVD%NUMFLDS ,&
       & NSTEP  ,&
       & NTSSG  , YSP_SBD%NLEVS   ,&
       & KBL    , KGPCOMP, YDCFU%NFRRC, PDTPHY,YDCSGEOM%RINDX,YDCSGEOM%RINDY, LLXFUMSE,&
       & PHI0  , PRE0  , PHIF0 , PRE0F ,PXYB0(1,1,YYTXYB0_PHY%M_ALPH),ZVVARG,&
       & ZVVD2    ,&
       & PXYB0(1,1,YYTXYB0_PHY%M_DELP),ZVVIVEG,ZVVLAI,&
       & PXYB0(1,1,YYTXYB0_PHY%M_LNPR),PXYB0(1,1,YYTXYB0_PHY%M_RDELP),&
       & ZVVRSMIN , ZVVSAB ,&
       & ZVVZ0H , ZVASEA , ZVALAN ,&
       & ZVASOO , ZVADES , ZVASUL ,&
       & ZVAVOL ,YDGSGEOM%RCORI, ZP1EXT,&
       & ZT0U,ZT0V,ZT0T,&
       & ZPQ,ZPI,ZPL,&
       & ZPLCONV,ZPICONV,ZPRCONV,ZPSCONV,&
       & ZPS,ZPR,ZPG,ZPTKE,&
       & ZPEFB1,ZPEFB2,ZPEFB3,&
       & ZPCVV,ZPO3,ZP1CHEM,ZP1NOGW,ZP2NOGW, PGFL, &
       & ZT0VOR,&
       & PRCP0(1,1,YYTRCP0%M_CP), ZCVGQ  ,PRCP0(1,1,YYTRCP0%M_R), PKOZO  , ZFPLCH , ZFPLSH ,&
       & ZPDAL,ZPDOM,ZPUEN,ZPUAL,&
       & ZPUOM,ZPUNEBH, PCTY0(1,0,YYTCTY0%M_EVEL),&
       & ZPRKTH,ZPRKTQV,ZPRKTQC, ZPTTE,&
       & ZPMXL,ZPSHTUR,ZPFQTUR, ZPFSTUR,&
       & ZVHTCCH ,  ZVHSCCH , ZVHBCCH ,&
       & ZVHSPSH ,ZVHPBLH ,&
       & ZVHQSH ,ZVKUDGRO ,&
       & PGPAR , PCUCONVCA, PNLCONVCA,&
       & ZSGF0 , ZSGA0, ZSGR0,&
       & ZSBT0 , ZRRT0 , ZRRFC0 ,&
       & ZSBQ0 , ZSBTL0, ZRRW0  ,&
       & ZRRIC0,&
       & ZVVHV    , PCTY0(1,1,YYTCTY0%M_VVEL),&
       & PEMTD , PEMTU ,PTRSW ,&
       & ZVVALV  ,&
       & ZVFALBF , ZVFALBSF , ZVFEMISF ,&
       & ZVFGETRL , ZVFLSM , ZVFVEG ,&
       & ZVFZ0F , ZVFZ0RLF,&
       & ZVFVRLAN , ZVFVRLDI , ZVCVC1  ,&
       & ZSFLSFL1 , ZSFOSFO1 ,&
       & PRMOON ,&
       & ZMU0   , ZMU0LU , ZMU0M  ,ZMU0N,YDGSGEOM%GELAM,YDGSGEOM%GEMU,YDGSGEOM%GM,&
       & ZAC    , ZAC_HC , ZMCOR , ZMMU0  , PDHSF  ,&
       & ZMRAB3C,ZMRAB3N,&
       & ZMRAB4C,ZMRAB4N,&
       & ZMRAB6C,ZMRAB6N,&
       & ZMRAT1C,ZMRAT1N,&
       & ZMRAT2C,ZMRAT2N,&
       & ZMRAT3C,ZMRAT3N,&
       & ZMRAT4C,ZMRAT4N,&
       & ZMRAT5C,ZMRAT5N,&
       & YDOROG%OROG,&
       & PWT0,ZT0DIV,ZT0UL,ZT0VL,PWT0L,PWT0M,&
       & PGDEOSI,PGUEOSI,PGMU0,PGMU0_MIN,PGMU0_MAX,&
       & PGDEOTI,PGDEOTI2,PGUEOTI,PGUEOTI2,PGEOLT,PGEOXT,&
       & PGRPROX,PGMIXP,PGFLUXC,PGRSURF,&
       & PDIFCQ , PDIFCQI, PDIFCQL, PDIFCS , ZDIFEXT, PDIFTQ , PDIFTQI,&
       & PDIFTQL,&
       & PDIFTS , PFCCQL , PFCCQN , PFCSQL , PFCSQN , PFCQNG ,&
       & PFCQING, PFCQLNG, PFCQRNG, PFCQSNG, PFCQGNG,&
       & PFPLCL , PFPLCN , PFPLCG , PFPLSL , PFPLSN , PFPLSG , PFRSO  ,&
       & PFRSOC ,&
       & PFRTH  , PFRTHC , PSTRCU , PSTRCV , PSTRDU , PSTRDV ,&
       & PSTRTU ,&
       & PSTRTV , PSTRMU , PSTRMV , PFRMH  , ZFRMQ  ,&
       & PDIFCQLC,PDIFCQIC,PFIMCC,&
       & PFEDQLC, PFEDQIC, PFEDQRC, PFEDQSC, PFCNEGQLC,PFCNEGQIC,PFCNEGQRC,PFCNEGQSC,&
       & PFCHOZ ,&
       & ZCPS   , ZLHS   ,&
       & ZRS    , ZLH    , ZLSCPE , PNEB   , PQICE  , PQLI   ,&
       & ZPA1,ZPLRAD1,ZPIRAD1,&
       & ZQRCONV, ZQSCONV,&
       & ZQSAT  ,&
       & ZQW    , PRH    , ZTW    , PALBDG   , ZCD    , ZCDN   ,&
       & ZCH    ,&
       & ZC1    , ZC2    , PCT    , ZEMIS  , PFCHSP , PFCLL  ,&
       & PFCLN  ,&
       & PFCS   , ZFEVI  , PFEVL  , PFEVN  , PFEVV  , PFLASH,PFTR   ,&
       & PFLWSP , PFONTE , PSP_SG(1,YSP_SG%YT%MP1), PFGEL  , PFGELS ,&
       & PFRSGNI, PFRSDNI, PFRSODS, PFRSOPS, PFRSOPT, PFRSOLU, PFRTHDS,&
       & PFPFPSL, PFPFPSN, PFPFPSG, PFPFPCL, PFPFPCN, PFPEVPSL, PFPEVPSN,PFPEVPSG, PFPEVPCL,&
       & PFPEVPCN,PFPEVPCG,ZFTKE  , ZFTKEI,  ZFEFB1 , ZFEFB2 , ZFEFB3  ,PGZ0   , PGZ0H   ,&
       & ZNEIJ  , ZVEG   , PQS,&
       & ZQSATS , PRUISL , PRUISP , PRUISS ,&
       & PUCLS  , PVCLS  , PNUCLS , PNVCLS , PTCLS  , PMRT   ,&
       & PQCLS  , PRHCLS , PCLCT  , PCLCH  , PCLCM  , PCLCL  ,&
       & PCLCC  , PCAPE  , PCTOP  , ICLPH  , PCLPH  , PVEIN  , PUGST  , PVGST  , PDIAGH, &
       & ZP1EZDIAG, PTPWCLS,ZDPRECIPS,ZDPRECIPS2,PVISICLD,PVISIHYDRO,PMXCLWC,&
       & ZSGF1, ZTENDPTKE, ZKUROV_H, ZKTROV_H,&
       & PDERNSHF, ZTENDEXT_DEP, ZVDSUND,PTRAJ_PHYS,&
       & ZEDR, YDDDH, ZSBT1, ZRRW1, &
       & ZSBQ1,ZRRIC1,ZSBTL1,&
       & ZRRFC1,ZSGA1,ZSGR1,PFTCNS,&
       & ZGP2DSPP)

      IF (LAJUCV) THEN
        DO JLEV=1,NFLEVG
          DO JROF=KST,KEND
            ZT0T(JROF,JLEV)=ZTAUX(JROF,JLEV)
          ENDDO
        ENDDO
        DO JLEV=1,NFLEVG
          DO JROF=KST,KEND
            PB1(JROF,ISLB1T9+JLEV-NFLSA)=PB1(JROF,ISLB1T9+JLEV-NFLSA)+ZDTAJU(JROF,JLEV)
          ENDDO
        ENDDO
      ENDIF

    ELSE

      ! IF (LAJUCV) THEN
      !   missing code under LAJUCV for leap-frog schemes.
      ! ENDIF

      ITDIA=1_JPIM
      CALL APLPAR(YDGEOMETRY,YDSURF,  YDXFU, YDCFU,YDMODEL,  KST,    KEND   , NPROMA ,&
       & ITDIA  , NFLEVG  , KSTGLO ,&
       & NVCLIS , YSD_VVD%NUMFLDS ,&
       & NSTEP  ,&
       & NTSSG  , YSP_SBD%NLEVS   ,&
       & KBL    , KGPCOMP, YDCFU%NFRRC, PDTPHY,YDCSGEOM%RINDX,YDCSGEOM%RINDY, LLXFUMSE,&
       & PHI9  , PRE9  , PHIF9 , PRE9F , PXYB9(1,1,YYTXYB9_PHY%M_ALPH) , ZVVARG   ,&
       & ZVVD2    ,&
       & PXYB9(1,1,YYTXYB9_PHY%M_DELP),ZVVIVEG,ZVVLAI,&
       & PXYB9(1,1,YYTXYB9_PHY%M_LNPR),PXYB9(1,1,YYTXYB9_PHY%M_RDELP),&
       & ZVVRSMIN , ZVVSAB ,&
       & ZVVZ0H , ZVASEA , ZVALAN ,&
       & ZVASOO , ZVADES , ZVASUL ,&
       & ZVAVOL,YDGSGEOM%RCORI, ZP1EXT9,&
       & ZT9U,ZT9V,ZT9T,&
       & ZPQ9,ZPI9,ZPL9,&
       & ZPLCONV,ZPICONV,ZPRCONV,ZPSCONV,&
       & ZPS9,ZPR9,ZPG9,ZPTKE9,&
       & ZPEFB19,ZPEFB29,ZPEFB39,&
       & ZPCVV9,ZPO39,ZP1CHEM9,ZP1NOGW9,ZP2NOGW9, PGFL, &
       & ZT9VOR,&
       & PRCP9(1,1,YYTRCP9%M_CP), ZCVGQ  ,PRCP9(1,1,YYTRCP9%M_R), PKOZO  , ZFPLCH, ZFPLSH  ,&
       & ZPDAL,ZPDOM,ZPUEN,ZPUAL,&
       & ZPUOM,ZPUNEBH, PCTY0(1,0,YYTCTY0%M_EVEL),&
       & ZPRKTH,ZPRKTQV,ZPRKTQC, ZPTTE,&
       & ZPMXL,ZPSHTUR,ZPFQTUR,ZPFSTUR,&
       & ZVHTCCH ,  ZVHSCCH , ZVHBCCH ,&
       & ZVHSPSH , ZVHPBLH ,&
       & ZVHQSH ,ZVKUDGRO ,&
       & PGPAR , PCUCONVCA, PNLCONVCA,&
       & ZSGF9 , ZSGA9, ZSGR9,&
       & ZSBT9 , ZRRT9 , ZRRFC9 ,&
       & ZSBQ9 , ZSBTL9, ZRRW9  ,&
       & ZRRIC9,&
       & ZVVHV    , PCTY0(1,1,YYTCTY0%M_VVEL),&
       & PEMTD , PEMTU ,PTRSW ,&
       & ZVVALV  ,&
       & ZVFALBF , ZVFALBSF , ZVFEMISF ,&
       & ZVFGETRL , ZVFLSM , ZVFVEG ,&
       & ZVFZ0F , ZVFZ0RLF,&
       & ZVFVRLAN , ZVFVRLDI , ZVCVC1  ,&
       & ZSFLSFL1 , ZSFOSFO1 , PRMOON ,&
       & ZMU0   , ZMU0LU , ZMU0M  , ZMU0N,YDGSGEOM%GELAM,YDGSGEOM%GEMU,YDGSGEOM%GM,&
       & ZAC    , ZAC_HC , ZMCOR , ZMMU0  , PDHSF  ,&
       & ZMRAB3C,ZMRAB3N,&
       & ZMRAB4C,ZMRAB4N,&
       & ZMRAB6C,ZMRAB6N,&
       & ZMRAT1C,ZMRAT1N,&
       & ZMRAT2C,ZMRAT2N,&
       & ZMRAT3C,ZMRAT3N,&
       & ZMRAT4C,ZMRAT4N,&
       & ZMRAT5C,ZMRAT5N,&
       & YDOROG%OROG,&
       & PWT9,ZT9DIV,ZT9UL,ZT9VL,ZWT9L,ZWT9M,&
       & PGDEOSI,PGUEOSI,PGMU0,PGMU0_MIN,PGMU0_MAX,&
       & PGDEOTI,PGDEOTI2,PGUEOTI,PGUEOTI2,PGEOLT,PGEOXT,&
       & PGRPROX,PGMIXP,PGFLUXC,PGRSURF,&
       & PDIFCQ , PDIFCQI, PDIFCQL, PDIFCS , ZDIFEXT, PDIFTQ , PDIFTQI,&
       & PDIFTQL,&
       & PDIFTS , PFCCQL , PFCCQN , PFCSQL , PFCSQN , PFCQNG ,&
       & PFCQING, PFCQLNG, PFCQRNG, PFCQSNG,PFCQGNG,&
       & PFPLCL , PFPLCN , PFPLCG , PFPLSL , PFPLSN , PFPLSG , PFRSO  ,&
       & PFRSOC ,&
       & PFRTH  , PFRTHC , PSTRCU , PSTRCV , PSTRDU , PSTRDV ,&
       & PSTRTU ,&
       & PSTRTV , PSTRMU , PSTRMV , PFRMH  , ZFRMQ  ,&
       & PDIFCQLC,PDIFCQIC,PFIMCC,&
       & PFEDQLC, PFEDQIC, PFEDQRC, PFEDQSC, PFCNEGQLC,PFCNEGQIC,PFCNEGQRC,PFCNEGQSC,&
       & PFCHOZ ,&
       & ZCPS   , ZLHS   ,&
       & ZRS    , ZLH    , ZLSCPE , PNEB   , PQICE  , PQLI   ,&
       & ZPA1,ZPLRAD1,ZPIRAD1,&
       & ZQRCONV, ZQSCONV,&
       & ZQSAT  ,&
       & ZQW    , PRH    , ZTW    , PALBDG   , ZCD    , ZCDN   ,&
       & ZCH    ,&
       & ZC1    , ZC2    , PCT    , ZEMIS  , PFCHSP , PFCLL  ,&
       & PFCLN  ,&
       & PFCS   , ZFEVI  , PFEVL  , PFEVN  , PFEVV  , PFLASH, PFTR   ,&
       & PFLWSP , PFONTE , ZSGT1, PFGEL  , PFGELS ,&
       & PFRSGNI, PFRSDNI, PFRSODS, PFRSOPS, PFRSOPT, PFRSOLU, PFRTHDS,&
       & PFPFPSL, PFPFPSN, PFPFPSG, PFPFPCL, PFPFPCN, PFPEVPSL,PFPEVPSN,PFPEVPSG, PFPEVPCL,&
       & PFPEVPCN,PFPEVPCG,ZFTKE  , ZFTKEI,  ZFEFB1 , ZFEFB2 , ZFEFB3  ,PGZ0   , PGZ0H   ,&
       & ZNEIJ  , ZVEG   , PQS,&
       & ZQSATS , PRUISL , PRUISP , PRUISS ,&
       & PUCLS  , PVCLS  , PNUCLS , PNVCLS , PTCLS  , PMRT   ,&
       & PQCLS  , PRHCLS , PCLCT  , PCLCH  , PCLCM  , PCLCL  ,&
       & PCLCC  , PCAPE  , PCTOP  , ICLPH  , PCLPH  , PVEIN  , PUGST  , PVGST  , PDIAGH, &
       & ZP1EZDIAG, PTPWCLS,ZDPRECIPS,ZDPRECIPS2,PVISICLD,PVISIHYDRO,PMXCLWC,&
       & ZSGF1, ZTENDPTKE, ZKUROV_H, ZKTROV_H ,&
       & PDERNSHF, ZTENDEXT_DEP, ZVDSUND,PTRAJ_PHYS,&
       & ZEDR,YDDDH, ZSBT1, ZRRW1, &
       & ZSBQ1,ZRRIC1,ZSBTL1,&
       & ZRRFC1,ZSGA1,ZSGR1,PFTCNS,&
       & ZGP2DSPP)

      ! IF (LAJUCV) THEN
      !   missing code under LAJUCV for leap-frog schemes.
      ! ENDIF

    ENDIF

  !    convert to flexible interface structure
  IF (LINTFLEX) THEN
    CALL APLPAR2INTFLEX(YGFL,YDPHY,NPROMA, KST, KEND, NFLEVG,&
      & PDIFCQ , PDIFCQI, PDIFCQL, PDIFCS ,ZDIFEXT,&
      & PDIFTQ , PDIFTQI, PDIFTQL, PDIFTS ,&
      & PFCCQL , PFCCQN , PFCSQL , PFCSQN ,&
      & PFPLSL , PFPLSN , PFPLCL,  PFPLCN,&
      & PFPEVPSL,PFPEVPSN,PFPEVPCL,PFPEVPCN,&
      & PFPFPSL, PFPFPSN, PFPFPCL ,PFPFPCN ,&
      & PFCQLNG, PFCQING, PFCQRNG, PFCQSNG,&
      & PFCQNG , PFRMH  , ZFRMQ  , PFRSO  , PFRTH  ,&
      & PSTRCU , PSTRCV , PSTRDU , PSTRDV ,&
      & PSTRTU , PSTRTV , PSTRMU , PSTRMV ,&
      & PDIFCQLC,PDIFCQIC,PFIMCC,&
      & PFEDQLC, PFEDQIC, PFEDQRC, PFEDQSC, PFCNEGQLC,PFCNEGQIC,PFCNEGQRC,PFCNEGQSC,&
      & ZFTKE,&
      & ZTENDPTKE, ZTENDEXT, ZTENDEXT_DEP,&
      & YLPROCSET )
  ENDIF

ENDIF ! LLDIAB.AND..NOT.LMPA

!        2.3  Computes MOCON in the CLP.
!             --------------------------
PMOCON(KST:KEND)=0.0_JPRB
IF (LLDIAB.AND.(.NOT.LMPA)) THEN

  IF (LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF = KST, KEND
        PMOCON(JROF) = PMOCON(JROF)+(ZLCVQ(JROF,JLEV)-ZPQ(JROF,JLEV)*&
         & ZT0DIV(JROF,JLEV))*PXYB0(JROF,JLEV,YYTXYB0_PHY%M_DELP)&
         & *MAX(0,SIGN(1,JLEV-ICLPH(JROF)))  
      ENDDO
    ENDDO
    DO JROF = KST, KEND
      PMOCON(JROF) = PMOCON(JROF)/(PRE0(JROF,NFLEVG)-PRE0(JROF,ICLPH(JROF)-1))
    ENDDO
  ELSE
    DO JLEV=1,NFLEVG
      DO JROF = KST, KEND
        PMOCON(JROF) = PMOCON(JROF)+(ZLCVQ(JROF,JLEV)-ZPQ9(JROF,JLEV)*&
         & ZT0DIV(JROF,JLEV))*PXYB9(JROF,JLEV,YYTXYB9_PHY%M_DELP)&
         & *MAX(0,SIGN(1,JLEV-ICLPH(JROF)))  
      ENDDO
    ENDDO
    DO JROF = KST, KEND
      PMOCON(JROF) = PMOCON(JROF)/(PRE9(JROF,NFLEVG)-PRE9(JROF,ICLPH(JROF)-1))
    ENDDO
  ENDIF
ENDIF

! Store surface water flux P and E for water conservation
IF (LCORWAT) THEN
  PSD_VH(KST:KEND,YSD_VH%YPSL%MP) = PFPLSL(KST:KEND,NFLEVG)
  PSD_VH(KST:KEND,YSD_VH%YPCL%MP) = PFPLCL(KST:KEND,NFLEVG)
  PSD_VH(KST:KEND,YSD_VH%YPSN%MP) = PFPLSN(KST:KEND,NFLEVG)
  PSD_VH(KST:KEND,YSD_VH%YPCN%MP) = PFPLCN(KST:KEND,NFLEVG)
  PSD_VH(KST:KEND,YSD_VH%YEVA%MP) = PFEVN(KST:KEND,1)-PFEVL(KST:KEND,1)
ENDIF

!        2.4  Stores radiation coefficients.
!             ------------------------------

! * writes grid-point transmission coefficients for simplified physics.

IF (LRCOEF.AND.(NSTEP == 1).AND.LLDIAB.AND.(.NOT.LMPA)) THEN
  IFIELDSS=NG3SR*NFLEVG
  CALL WRRADCOEF(YDGEOMETRY,YDRCOEF,KST,KEND,KSTGLO,IFIELDSS,ZRADTC,ZAC_HC)
ENDIF

!       2.5   Ozone
!             -----

IF (LLDIAB.AND.LOZONE.AND.(.NOT.LMPA)) THEN
  ! * Caution: this part has not been yet validated relative
  !   to the GFL implementation, and LOZONE (the setup of
  !   which has not yet been updated) can be true only if
  !   the GFL ozone is activated as a prognostic and advected
  !   variable.
  IPO3=(YO3%MP_SL1-1)*(NFLEVG+2*NFLSUL)
  IF (LSLAG) THEN
    IF (LTWOTL) THEN
      CALL CPOZO (NPROMA,KST,KEND,NFLEVG,PDTPHY,PFCHOZ,&
       & PB1(1,ISLB1GFL9+IPO3+1-NFLSA),PXYB0(1,1,YYTXYB0_PHY%M_RDELP))  
    ELSE
      CALL CPOZO (NPROMA,KST,KEND,NFLEVG,PDTPHY,PFCHOZ,&
       & PB1(1,ISLB1GFL9+IPO3+1-NFLSA),PXYB9(1,1,YYTXYB9_PHY%M_RDELP))  
    ENDIF
  ELSE
    CALL CPOZO (NPROMA,KST,KEND,NFLEVG,PDTPHY,PFCHOZ,&
     & ZPO31,PXYB9(1,1,YYTXYB9_PHY%M_RDELP))
  ENDIF
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

IF (LLDIAB.AND.(NDPSFI == 1)) THEN
  CALL CPQSOL(YDGEOMETRY%YRDIMV,YDPHY,NPROMA,KST,KEND,PRE0,ZRRT0,PQS,ZQSATS,PQSOL)
ENDIF

!        2.7  Computation of tendencies T,u,v and Q.
!             --------------------------------------

! Set GFL tendencies to 0

ZTENDGFL(:,:,:) = 0.0_JPRB

IF (LLDIAB.AND.(.NOT.LSIMPH).AND.(.NOT.LMPA)) THEN

  TSPHY = MAX(PDTPHY,1.0_JPRB)

  ! * CPTEND+CPUTQY = Old( CPATY + CPDUP + CPDTHP )
  ! Calcul des tendances de T , U et de Q et modifications
  ! eventuelles de W et de OMEGA/P

  IF (LTWOTL) THEN

    IF (LINTFLEX.AND.(.NOT.LDCONFX)) THEN
      CALL CPTEND_FLEX( YDLDDH,YDMDDH,YGFL,YDPHY,NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & PXYB0(1,1,YYTXYB0_PHY%M_DELP) ,&
       & PXYB0(1,1,YYTXYB0_PHY%M_RDELP), PRCP0(1,1,YYTRCP0%M_CP),&
       & ZT0U,ZT0V,ZT0T,ZRRT0,&
       & PGFL,&
       & YLPROCSET,&
       & PTENDU , PTENDV , ZTENDH , ZTENDGFL,&
       & PFHSCL ,PFHSCN,PFHSSL,PFHSSN,&
       & PFHPCL ,PFHPCN,PFHPSL,PFHPSN,&
       & ZFHP   ,ZFP   ,  PFEPFP, PFCMPCQ, PFCMPSN, PFCMPSL,YDDDH )
    ELSE
      CALL CPTEND_NEW( YDMODEL, NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & PDIFCQ , PDIFCQI, PDIFCQL, PDIFCS ,ZDIFEXT,&
       & PDIFTQ , PDIFTQI, PDIFTQL, PDIFTS ,&
       & PFCCQL , PFCCQN , PFCSQL , PFCSQN ,&
       & PFPLSL , PFPLSN , PFPLSG , PFPLCL,  PFPLCN, PFPLCG,&
       & PFPEVPSL,PFPEVPSN,PFPEVPSG,PFPEVPCL,PFPEVPCN,PFPEVPCG,&
       & PFPFPSL, PFPFPSN, PFPFPSG, PFPFPCL ,PFPFPCN ,&
       & PFCQLNG, PFCQING, PFCQRNG, PFCQSNG, PFCQGNG ,&
       & PFCQNG , PFRMH  , ZFRMQ  , PFRSO  , PFRTH  ,&
       & PSTRCU , PSTRCV , PSTRDU , PSTRDV ,&
       & PSTRTU , PSTRTV , PSTRMU , PSTRMV ,&
       & PDIFCQLC,PDIFCQIC,PFIMCC,&
       & PFEDQLC, PFEDQIC, PFEDQRC, PFEDQSC, PFCNEGQLC,PFCNEGQIC,PFCNEGQRC,PFCNEGQSC,&
       & ZFTKE  , ZFTKEI,  ZFEFB1 , ZFEFB2 , ZFEFB3 ,PXYB0(1,1,YYTXYB0_PHY%M_DELP) ,&
       & PXYB0(1,1,YYTXYB0_PHY%M_RDELP), PHIF0  , PRCP0(1,1,YYTRCP0%M_CP),&
       & ZT0U,ZT0V,ZT0T,&
       & ZPQ,ZPI,ZPL,&
       & ZPLCONV,ZPICONV,ZPRCONV,ZPSCONV,&
       & ZPR,ZPS,ZPG,&
       & ZCPS   , ZRRT0  ,&
       & PFHSCL ,PFHSCN,PFHSSL,PFHSSN,PFHSSG,&
       & PFHPCL ,PFHPCN,PFHPCG,PFHPSL,PFHPSN,PFHPSG,&
       & ZFHP   ,ZFP   ,  PFEPFP, PFCMPCQ, PFCMPSN, PFCMPSL,&
       & PTENDU , PTENDV , ZTENDU, ZTENDV, ZTENDH ,&
       & ZPTENDQ1,ZPTENDI1,ZPTENDL1,&
       & ZPTENDLCONV1,ZPTENDICONV1,&
       & ZPTENDRCONV1,ZPTENDSCONV1,&
       & ZPTENDR1,ZPTENDS1,ZPTENDG1,ZPTENDTKE1,&
       & ZPTENDEFB11,ZPTENDEFB21,ZPTENDEFB31,&
       & ZTENDEXT,YDDDH)
    ENDIF
  ELSE
    IF (LINTFLEX.AND.(.NOT.LDCONFX)) THEN
      CALL CPTEND_FLEX( YDLDDH,YDMDDH,YGFL,YDPHY,NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & PXYB9(1,1,YYTXYB9_PHY%M_DELP) ,&
       & PXYB9(1,1,YYTXYB9_PHY%M_RDELP), PRCP9(1,1,YYTRCP9%M_CP),&
       & ZT9U,ZT9V,ZT9T,ZRRT9,&
       & PGFL,&
       & YLPROCSET,&
       & PTENDU , PTENDV , ZTENDH , ZTENDGFL,&
       & PFHSCL ,PFHSCN,PFHSSL,PFHSSN,&
       & PFHPCL ,PFHPCN,PFHPSL,PFHPSN,&
       & ZFHP   ,ZFP   ,  PFEPFP, PFCMPCQ, PFCMPSN, PFCMPSL,YDDDH )
    ELSE
      CALL CPTEND_NEW( YDMODEL, NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & PDIFCQ , PDIFCQI, PDIFCQL, PDIFCS ,ZDIFEXT,&
       & PDIFTQ , PDIFTQI, PDIFTQL, PDIFTS ,&
       & PFCCQL , PFCCQN , PFCSQL , PFCSQN ,&
       & PFPLSL , PFPLSN , PFPLSG , PFPLCL,  PFPLCN, PFPLCG,&
      & PFPEVPSL,PFPEVPSN,PFPEVPSG,PFPEVPCL,PFPEVPCN,PFPEVPCG,&
       & PFPFPSL, PFPFPSN, PFPFPSG, PFPFPCL ,PFPFPCN,&
       & PFCQLNG, PFCQING, PFCQRNG, PFCQSNG, PFCQGNG,&
       & PFCQNG , PFRMH  , ZFRMQ  , PFRSO  , PFRTH  ,&
       & PSTRCU , PSTRCV , PSTRDU , PSTRDV ,&
       & PSTRTU , PSTRTV , PSTRMU , PSTRMV ,&
       & PDIFCQLC,PDIFCQIC,PFIMCC,&
       & PFEDQLC, PFEDQIC, PFEDQRC, PFEDQSC, PFCNEGQLC,PFCNEGQIC,PFCNEGQRC,PFCNEGQSC,&
       & ZFTKE  , ZFTKEI,  ZFEFB1 , ZFEFB2 , ZFEFB3 , PXYB9(1,1,YYTXYB9_PHY%M_DELP) ,&
       & PXYB9(1,1,YYTXYB9_PHY%M_RDELP), PHIF9  , PRCP9(1,1,YYTRCP9%M_CP),&
       & ZT9U,ZT9V,ZT9T,&
       & ZPQ9,ZPI9,ZPL9,&
       & ZPLCONV,ZPICONV,ZPRCONV,ZPSCONV,&
       & ZPR9,ZPS9,ZPG9,&
       & ZCPS   , ZRRT9,&
       & PFHSCL ,PFHSCN,PFHSSL,PFHSSN,PFHSSG,&
       & PFHPCL ,PFHPCN,PFHPCG,PFHPSL,PFHPSN,PFHPSG,&
       & ZFHP   ,ZFP   ,  PFEPFP, PFCMPCQ, PFCMPSN, PFCMPSL,&
       & PTENDU , PTENDV , ZTENDU, ZTENDV, ZTENDH ,&
       & ZPTENDQ1,ZPTENDI1,ZPTENDL1,&
       & ZPTENDLCONV1,ZPTENDICONV1,&
       & ZPTENDRCONV1,ZPTENDSCONV1,&
       & ZPTENDR1,ZPTENDS1,ZPTENDG1,ZPTENDTKE1,&
       & ZPTENDEFB11,ZPTENDEFB21,ZPTENDEFB31,&
       & ZTENDEXT,YDDDH)
    ENDIF
    
    IF ( L3MT.OR.LSTRAPRO.OR.(NDPSFI==1)) THEN
!     PFEPFP was ZFEPFP in CPTEND_NEW, before, ZFEPFP still in CPFHPFS
      DO JLEV= 0, NFLEVG 
        DO JROF = 1, NPROMA
          PFEPFP(JROF,JLEV) = 0.0_JPRB
          PFCMPCQ(JROF,JLEV) = 0.0_JPRB
          PFCMPSN(JROF,JLEV) = 0.0_JPRB
          PFCMPSL(JROF,JLEV) = 0.0_JPRB
        ENDDO
      ENDDO
    ENDIF

  ENDIF ! LTWOTL



!        2.7.1  Diagnostics on physical tendencies
!               ----------------------------------

  IF (.NOT.LDCONFX) THEN
    IF ((GCHETN%LFREQD).OR.(GCHETN%LCOORD).OR.(GCHETN%LPROFV)) THEN
      IF (LTWOTL) THEN
        CALL CPCHET(  YDRIP,YDPHY, NPROMA, KST, KEND, NFLEVG, NSTEP,&
        & PFHSCL , PFHSCN , PFHSSL , PFHSSN ,&
        & PFHPCL , PFHPCN , PFHPSL , PFHPSN ,&
        & PDIFCQ , PDIFCQI, PDIFCQL, PDIFCS ,&
        & PDIFTQ , PDIFTQI, PDIFTQL, PDIFTS ,&
        & PFCCQL , PFCCQN , PFCSQL , PFCSQN ,&
        & PFPLCL , PFPLCN , PFPLSL , PFPLSN ,&
        & PFPFPSL, PFPFPSN, PFPFPCL, PFPFPCN,&
        & PFPEVPSL,PFPEVPSN,PFPEVPCL,PFPEVPCN,&
        & PFRMH  , ZFRMQ  , PFRSO  , PFRTH  ,&
        & PSTRCU , PSTRCV , PSTRDU , PSTRDV ,&
        & PSTRTU , PSTRTV , PSTRMU , PSTRMV ,&
        & PXYB0(1,1,YYTXYB0_PHY%M_RDELP), PRCP0(1,1,YYTRCP0%M_CP),&
        & ZT0T,ZPQ,ZPI,ZPL,&
        & ZCPS   , ZRRT0 , PQS  ,&
        & PTENDU , PTENDV , ZTENDH ,&
        & ZPTENDQ1,ZPTENDI1,ZPTENDL1,&
        & ZPTENDR1,ZPTENDS1,&
        & YDGSGEOM%GEMU,YDGSGEOM%GELAM, PHI0   , PHIF0  )
      ELSE
        CALL CPCHET(  YDRIP,YDPHY, NPROMA, KST, KEND, NFLEVG, NSTEP,&
        & PFHSCL , PFHSCN , PFHSSL , PFHSSN ,&
        & PFHPCL , PFHPCN , PFHPSL , PFHPSN ,&
        & PDIFCQ , PDIFCQI, PDIFCQL, PDIFCS ,&
        & PDIFTQ , PDIFTQI, PDIFTQL, PDIFTS ,&
        & PFCCQL , PFCCQN , PFCSQL , PFCSQN ,&
        & PFPLCL , PFPLCN , PFPLSL , PFPLSN ,&
        & PFPFPSL, PFPFPSN, PFPFPCL, PFPFPCN,&
        & PFPEVPSL,PFPEVPSN,PFPEVPCL,PFPEVPCN,&
        & PFRMH  , ZFRMQ  , PFRSO  , PFRTH  ,&
        & PSTRCU , PSTRCV , PSTRDU , PSTRDV ,&
        & PSTRTU , PSTRTV , PSTRMU , PSTRMV ,&
        & PXYB9(1,1,YYTXYB9_PHY%M_RDELP), PRCP9(1,1,YYTRCP9%M_CP),&
        & ZT9T,ZPQ9,ZPI9,ZPL9,&
        & ZCPS   , ZRRT9  , PQS    ,&
        & PTENDU , PTENDV , ZTENDH ,&
        & ZPTENDQ1,ZPTENDI1,ZPTENDL1,&
        & ZPTENDR1,ZPTENDS1,&
        & YDGSGEOM%GEMU,YDGSGEOM%GELAM, PHI9   , PHIF9  )
      ENDIF
    ENDIF
    
    IF (GCHETN%LPROFV)&
     & CALL PROFILECHET(YDGEOMETRY,YDSURF,YDDPHY,YDRIP,YDMODEL%YRML_PHY_MF,&
     & KEND,&
     & YDGSGEOM%GELAM,YDGSGEOM%GEMU,&
     & YDGSGEOM%GM,ZMU0,YDOROG%OROG,POROGL,POROGM,YDGSGEOM%RCORI,YDCSGEOM%RATATH,YDCSGEOM%RATATX,&
     & ZVFLSM, ZVVARG, ZVVSAB, ZVFALBF,&
     & ZVFALBSF, ZVFEMISF, ZVVD2, ZVVIVEG,&
     & ZVVLAI, PCT, ZVFZ0F, ZVVZ0H, ZVFZ0RLF,&
     & ZVFGETRL, ZVFVRLAN, ZVFVRLDI, ZVVRSMIN,&
     & ZVFVEG, ZVVHV, ZVASEA, ZVALAN,&
     & ZVASOO, ZVADES, ZVCVC1,&
     & ZT0SP,ZT0SPL,ZT0SPM,&
     & ZT0T,ZT0TL,ZT0TM,&
     & ZPQ,ZPQL,ZPQM,&
     & ZT0U,ZT0V,ZT0VOR,ZT0DIV, ZCVGQ, ZLCVQ,&
     & ZRRT9, ZSBT9, ZRRFC9, ZRRW9,&
     & ZRRIC9, ZSBQ9, ZSBTL9, ZSGF9,&
     & PDIFCQ, PDIFCQI, PDIFCQL, PDIFCS, PDIFTQ, PDIFTQI, PDIFTQL,&
     & PDIFTS, PFCCQL, PFCCQN, PFCSQL, PFCSQN, PFCQNG, PFCQING,&
     & PFCQLNG, PFPLCL, PFPLCN, PFPLSL, PFPLSN, PFRSO, PFRTH,&
     & PSTRCU, PSTRCV, PSTRDU, PSTRDV, PSTRTU, PSTRTV, PSTRMU, PSTRMV,&
     & PFRMH, PFCHOZ, PNEB, PQICE, PQLI, PRH,&
     & PFCS, PFCLL, PFCLN, PFEVL, PFEVN, PFEVV, PFTR, PFLWSP, PFONTE,&
     & PFGEL, PFGELS, PFRSODS, PFRSOPS, PFRSOPT, PFRTHDS,&
     & PQS, PRUISL, PRUISP, PRUISS,&
     & PUCLS, PVCLS, PTCLS, PQCLS, PRHCLS,&
     & PCLCT, PCLCH, PCLCM, PCLCL, PCLCC,&
     & ZFPLCH, ZFPLSH,&
     & ZVHTCCH,  ZVHSCCH, ZVHBCCH, ZVHPBLH )

  ENDIF

ENDIF

!        2.8  Modification of vertical velocities
!             by some physics output when required.
!             -------------------------------------

IF (LLDIAB.AND.(.NOT.LSIMPH)) THEN

  ! * MODIFICATION DE LA VITESSE VERTICALE ET DE LA TENDANCE DE
  ! PRESSION DE SURFACE SI NDPSFI=1 ( MASSE VARIABLE ).
  ! Ajout de la physique dans l'equation de continuite/Add physics
  ! in continuity equation.

  IF (NDPSFI == 1) THEN
    IF (LSLAG .AND. LTWOTL) THEN
      CALL CPMVVPS(YDVAB,NPROMA,KST,KEND,NFLEVG,PDTPHY,&
       & ZFP,PRE0(1,NFLEVG),PFEVL,PFEVN,&
       & PCTY0(1,0,YYTCTY0%M_EVEL),PCTY0(1,0,YYTCTY0%M_PSDVBC),PB1(1,MSLB1SP9))
    ELSEIF (LSLAG .AND. (.NOT.LTWOTL)) THEN
      CALL CPMVVPS(YDVAB,NPROMA,KST,KEND,NFLEVG,PDTPHY,&
       & ZFP,PRE9(1,NFLEVG),PFEVL,PFEVN,&
       & PCTY0(1,0,YYTCTY0%M_EVEL),PCTY0(1,0,YYTCTY0%M_PSDVBC),PB1(1,MSLB1SP9))
    ELSE
      CALL CPMVVPS(YDVAB,NPROMA,KST,KEND,NFLEVG,PDTPHY,&
       & ZFP,PRE9(1,NFLEVG),PFEVL,PFEVN,&
       & PCTY0(1,0,YYTCTY0%M_EVEL),PCTY0(1,0,YYTCTY0%M_PSDVBC),ZT1SP)
    ENDIF
  ENDIF

ENDIF

!        2.9  Computation of evolution of T, u, v and Q.
!             ------------------------------------------

! * Calculation of IPGFL, since the old pointers
!   MSLB1[X]9 (=MSLB1GFL9+IP[X]) do not exist any longer in PTRSLB1.

! usefull pointer for new version of cputqy

DO JGFL=1,NUMFLDS
  IF ((YCOMP(JGFL)%MP1 > 0) .AND. (YCOMP(JGFL)%MP_SL1 > 0)) THEN
     IPGFL(YCOMP(JGFL)%MP1) = (YCOMP(JGFL)%MP_SL1-1)*(NFLEVG+2*NFLSUL)
  ENDIF   
ENDDO  

!  ALARO does not respect the coding rules, tendency of pseudo-TKE is computed in APLPAR and not
!  in CPTEND_NEW. To use the new version of cputqy it is then necessary to write it in GFL tendencies array.
! This memory transfer is not necessary, please respect coding rules to avoid it.

! Not necessary for intflex: already done in aplpar2intflex
IF (.NOT.(LINTFLEX.AND.(.NOT.LDCONFX))) THEN
  IF (LPTKE) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZPTENDTKE1(JROF,JLEV) = ZTENDPTKE(JROF,JLEV)
      ENDDO
    ENDDO    
  ENDIF
  ! Extra-GFL
  IF(LMDUST.AND.(NGFL_EXT/=0)) THEN
    DO JGFL=1, NGFL_EXT
      DO JLEV=1,NFLEVG
        DO JROF=KST,KEND
          ZTENDGFL(JROF,JLEV,YEXT(JGFL)%MP1) = ZTENDEXT(JROF,JLEV,JGFL)+&! turbulent tendency
                                             & ZTENDEXT_DEP(JROF,JLEV,JGFL) ! moist tendency
! moist tendency
! moist tendency
! moist tendency
! moist tendency
! moist tendency
! moist tendency
! moist tendency
! moist tendency
        ENDDO
      ENDDO 
    ENDDO   
  ENDIF
ENDIF

! ky: non-zero option not yet coded for the time being.
ZTENDD=0.0_JPRB

IF (LLDIAB.AND.(.NOT.LSIMPH).AND.(.NOT.LMPA)) THEN
  ! Calcul de T , Q et du Vent a l'instant 1

  CALL CPUTQY(YDGEOMETRY%YRDIMV,YDGMV,YGFL,YDPTRSLB1,YDPHY,NPROMA,KST,KEND,NFLEVG,PDTPHY,IPGFL,&
   & ISLB1T9,ISLB1U9,ISLB1V9,ISLB1VD9,ISLB1GFL9,&
   & ZTENDH, ZTENDT, PTENDU, PTENDV, ZTENDU, ZTENDV, ZTENDD, ZTENDGFL,&
   & PRCP0(1,1,YYTRCP0%M_CP),PXYB0(1,1,YYTXYB0_PHY%M_DELP),ZT0T,ZT0U,ZT0V,&
   & PRCP9(1,1,YYTRCP9%M_CP),PXYB9(1,1,YYTXYB9_PHY%M_DELP),ZT9T,ZT9U,ZT9V,&
   & PB1, PGMVT1, PGFLT1,&
   & PFDIS)

ENDIF

!        2.9a Evolution of precipitation fluxes
!             ------------------------------------------

IF(LNEBN.OR.LNEBR.OR.LRRGUST) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZPCPF1(JROF,JLEV)=ZFPLCH(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF
IF(LRRGUST) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZPSPF1(JROF,JLEV)=ZFPLSH(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!       2.9b Prognostic convection etc.
!            --------------------------

! TRANSFER NOT ADVECTED VARIABLES INTO PGFLT1
IF (LCVPRO.OR.LGPCMT) THEN
   IF (.NOT.YUAL%LADV) THEN
      ZPUAL1(KST:KEND,1:NFLEVG)=ZPUAL(KST:KEND,1:NFLEVG)
      ZPUOM1(KST:KEND,1:NFLEVG)=ZPUOM(KST:KEND,1:NFLEVG)
   ENDIF
   IF (LCDDPRO) THEN
    IF (.NOT.YDAL%LADV) THEN
      ZPDAL1(KST:KEND,1:NFLEVG)=ZPDAL(KST:KEND,1:NFLEVG)
      ZPDOM1(KST:KEND,1:NFLEVG)=ZPDOM(KST:KEND,1:NFLEVG)
    ENDIF
   ENDIF
ENDIF
IF (YUNEBH%LACTIVE.AND..NOT.YUNEBH%LADV)THEN
   ZPUNEBH1(KST:KEND,1:NFLEVG)=ZPUNEBH(KST:KEND,1:NFLEVG)
ENDIF
IF (YUEN%LACTIVE.AND..NOT.YUEN%LADV) THEN
   ZPUEN1(KST:KEND,1:NFLEVG)=ZPUEN(KST:KEND,1:NFLEVG)
ENDIF
IF (YTTE%LACTIVE.AND..NOT.YTTE%LADV) THEN
   ZPTTE1(KST:KEND,1:NFLEVG)=ZPTTE(KST:KEND,1:NFLEVG)
ENDIF
IF (YMXL%LACTIVE.AND..NOT.YMXL%LADV) THEN
   ZPMXL1(KST:KEND,1:NFLEVG)=ZPMXL(KST:KEND,1:NFLEVG)
ENDIF
IF (YSHTUR%LACTIVE.AND..NOT.YSHTUR%LADV) THEN
   ZPSHTUR1(KST:KEND,1:NFLEVG)=ZPSHTUR(KST:KEND,1:NFLEVG)
ENDIF
IF (YFQTUR%LACTIVE.AND..NOT.YFQTUR%LADV) THEN
   ZPFQTUR1(KST:KEND,1:NFLEVG)=ZPFQTUR(KST:KEND,1:NFLEVG)
ENDIF
IF (YFSTUR%LACTIVE.AND..NOT.YFSTUR%LADV) THEN
   ZPFSTUR1(KST:KEND,1:NFLEVG)=ZPFSTUR(KST:KEND,1:NFLEVG)
ENDIF
IF (YRKTH%LACTIVE)THEN
   ZPRKTH1(KST:KEND,1:NFLEVG)=ZPRKTH(KST:KEND,1:NFLEVG)
ENDIF
IF (YRKTQV%LACTIVE)THEN
   ZPRKTQV1(KST:KEND,1:NFLEVG)=ZPRKTQV(KST:KEND,1:NFLEVG)
ENDIF
IF (YRKTQC%LACTIVE)THEN
   ZPRKTQC1(KST:KEND,1:NFLEVG)=ZPRKTQC(KST:KEND,1:NFLEVG)
ENDIF

!        2.10  Surface variables.
!              ------------------

IF (LLDIAB.AND.LMPHYS.AND.(.NOT.LMPA).AND.(.NOT.LSFORCS)) THEN
  
  IF (.NOT.LMSE) THEN
    DO JLEV=0,NFLEVG
      DO JROF=KST,KEND
        ZFPLSN(JROF,JLEV)=PFPLSN(JROF,JLEV)+PFPLSG(JROF,JLEV)
      ENDDO
    ENDDO
    CALL CPTENDS( YDMODEL%YRML_PHY_MF, NPROMA, KST, KEND, NFLEVG, YSP_SBD%NLEVS, PDTPHY,&
     & PFPLCL, PFPLSL, PFPLCN, ZFPLSN,&
     & PFRSO, PFRTH,&
     & ZSGA1,PCT, ZC1, ZC2,&
     & PFCHSP, PFCLL, PFCLN, PFCS,&
     & ZFEVI,PFEVL, PFEVN,&
     & PFEVV, PFGEL, PFGELS, PFLWSP, PFONTE, PFTR,&
     & ZVFLSM, ZSGR1,&
     & PRUISL, PRUISP, PRUISS, ZSGF1, ZVEG,&
     & ZTDTS, ZTDTP, ZTDWS, ZTDWSI, ZTDWP, ZTDWPI, ZTDWL,&
     & ZTDSNS, ZTDALBNS, ZTDRHONS)  

    CALL CPWTS(YDSURF, YDMODEL%YRML_AOC%YRMCC,YDPHY,YDMODEL%YRML_PHY_MF%YRPHY1, NPROMA, KST,KEND, YSP_SBD%NLEVS, PDTPHY,&
     & ZTDTS, ZTDTP, ZTDWS, ZTDWSI, ZTDWP, ZTDWPI, ZTDWL,&
     & ZTDSNS, ZTDALBNS, ZTDRHONS,&
     & ZVPTPC,ZVPWPC,ZVFLSM,ZVVIVEG,&
     & ZRRT1,ZSBT1,ZRRW1,&
     & ZRRIC1,&
     & ZSBQ1,ZSBTL1,ZRRFC1,&
     & ZSGF1,ZSGA1,ZSGR1)  
  ELSE
    IF (LLXFUMSE) THEN
      DO JROF=KST,KEND
        ZRRT0(JROF)=PGPAR(JROF,MVTS)
      ENDDO
    ELSE
      DO JROF=KST,KEND
        ZRRT1(JROF)=PGPAR(JROF,MVTS)
      ENDDO
    ENDIF
  ENDIF
  IF(LNUDG)THEN

    ! * Calculation of IPQ since the old pointers
    !   MSLB1[X]9 (=MSLB1GFL9+IP[X]) do not exist any longer in PTRSLB1.
    IPQ=(YQ%MP_SL1-1)*(NFLEVG+2*NFLSUL)

    CALL CPNUDG ( NPROMA, KST,KEND, NFNUDG, NFLEVG, IBLK,&
     & XPNUDG,&
     & ZVFNUDM,&
     & ZRRT1,ZRRW1,&
     & ZSBQ1,ZSGF1,&
     & PB1(1,ISLB1T9+1-NFLSA),PB1(1,ISLB1GFL9+IPQ+1-NFLSA),&
     & PB1(1,ISLB1U9+1-NFLSA),PB1(1,ISLB1V9+1-NFLSA),&
     & PB1(1,MSLB1SP9),&
     & ZT0T,ZPQ,ZT0U,&
     & ZT0V,PRE0(1,NFLEVG),YDGSGEOM%GM,ZVFLSM)
  ENDIF
ENDIF

!        2.11 Evolution of CVV (GY)
!             ---------------------

IF(LCVPGY) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZPCVV1(JROF,JLEV)=ZPCVV(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!    -------------------------------------------------------------------

!*       3.    Simplified physics.
!              -------------------

!        3.1  Preliminary calculations necessary for simplified physics.
!             ----------------------------------------------------------

IF (LSIMPH) THEN

  ! * read grid-point transmission coefficients for simplified physics.
  IF (LRAYSP.AND.LRCOEF.AND.(NSTEP > 1.OR.LTLADDIA)) THEN
    IFIELDSS=NG3SR*NFLEVG
    CALL RDRADCOEF(YDGEOMETRY,YDRCOEF,KST,KEND,KSTGLO,IFIELDSS,ZRADTC,ZAC_HC)
  ENDIF

ENDIF

!        3.2  Simplified physics.
!             -------------------

IF (LSIMPH) THEN

  TSPHY = MAX(PDTPHY,1.0_JPRB)

  IF (LTWOTL) THEN

    CALL APLPASSH (YDPHY,YDMODEL%YRML_PHY_MF%YRPHY1,KST,KEND,NPROMA,NFLEVG,&
     & PRE0,ZPQ,&
     & ZSGF1,ZRRT1,ZRRW1,&
     & ZVFLSM,ZVFVEG,&
     & ZQS1)  

    IF (.NOT.LLDIAB) THEN
      CALL APLPASSH (YDPHY,YDMODEL%YRML_PHY_MF%YRPHY1,KST,KEND,NPROMA,NFLEVG,&
       & PRE0,ZPQ,&
       & ZSGF0,ZRRT0,ZRRW0,&
       & ZVFLSM,ZVFVEG,&
       & PQS)  
    ENDIF

    CALL APLPARS(YDGEOMETRY,YDRCOEF,YDMODEL%YRML_PHY_MF,KST,KEND,NPROMA,1,NFLEVG,NTSSG,NSTEP,&
     & PHI0,PRE0,PHIF0,PRE0F,PXYB0(1,1,YYTXYB0_PHY%M_DELP),PXYB0(1,1,YYTXYB0_PHY%M_RDELP),&
     & ZT0U,ZT0V,ZT0T,ZPQ,PRCP0(1,1,YYTRCP0%M_CP),ZCVGQ,&
     & ZSGF0,ZRRT0,PQS,&
     & ZRRT1,ZQS1,&
     & ZVFGETRL,ZVFLSM,&
     & ZVFZ0F,ZVFVRLAN,ZVFVRLDI,&
     & ZVFALBF,ZMU0,YDGSGEOM%GM,&
     & ZAC_HC,ZMCOR,&
     & ZMRAB3C,ZMRAB3N,&
     & ZMRAB4C,ZMRAB4N,&
     & ZMRAB6C,ZMRAB6N,&
     & ZMRAT1C,ZMRAT1N,&
     & ZMRAT2C,ZMRAT2N,&
     & ZMRAT3C,ZMRAT3N,&
     & ZMRAT4C,ZMRAT4N,&
     & ZMRAT5C,ZMRAT5N,&
     & PDIFCQ,PDIFCS,PDIFTQ,PDIFTS,&
     & PFCCQL,PFCCQN,PFCSQL,PFCSQN,&
     & PFPLCL,PFPLCN,PFPLSL,PFPLSN,PFRSO,PFRTH,&
     & PSTRCU,PSTRCV,PSTRDU,PSTRDV,PSTRTU,PSTRTV,&
     & PSTRMU,PSTRMV,PFRMH)

  ELSE

    CALL APLPASSH (YDPHY,YDMODEL%YRML_PHY_MF%YRPHY1,KST,KEND,NPROMA,NFLEVG,&
     & PRE0,ZPQ,&
     & ZSGF1,ZRRT1,ZRRW1,&
     & ZVFLSM,ZVFVEG,&
     & ZQS1)  

    IF (.NOT.LLDIAB) THEN
      CALL APLPASSH (YDPHY,YDMODEL%YRML_PHY_MF%YRPHY1,KST,KEND,NPROMA,NFLEVG,&
       & PRE9,ZPQ9,&
       & ZSGF0,ZRRT9,ZRRW9,&
       & ZVFLSM,ZVFVEG,&
       & PQS)  
    ENDIF

    CALL APLPARS(YDGEOMETRY,YDRCOEF,YDMODEL%YRML_PHY_MF,KST,KEND,NPROMA,1,NFLEVG,NTSSG,NSTEP,&
     & PHI9,PRE9,PHIF9,PRE9F,PXYB9(1,1,YYTXYB9_PHY%M_DELP),PXYB9(1,1,YYTXYB9_PHY%M_RDELP),&
     & ZT9U,ZT9V,ZT9T,ZPQ9,PRCP9(1,1,YYTRCP9%M_CP),ZCVGQ,&
     & ZSGF0,ZRRT9,PQS,&
     & ZRRT1,ZQS1,&
     & ZVFGETRL,ZVFLSM,&
     & ZVFZ0F,ZVFVRLAN,ZVFVRLDI,&
     & ZVFALBF,ZMU0,YDGSGEOM%GM,&
     & ZAC_HC,ZMCOR,&
     & ZMRAB3C,ZMRAB3N,&
     & ZMRAB4C,ZMRAB4N,&
     & ZMRAB6C,ZMRAB6N,&
     & ZMRAT1C,ZMRAT1N,&
     & ZMRAT2C,ZMRAT2N,&
     & ZMRAT3C,ZMRAT3N,&
     & ZMRAT4C,ZMRAT4N,&
     & ZMRAT5C,ZMRAT5N,&
     & PDIFCQ,PDIFCS,PDIFTQ,PDIFTS,&
     & PFCCQL,PFCCQN,PFCSQL,PFCSQN,&
     & PFPLCL,PFPLCN,PFPLSL,PFPLSN,PFRSO,PFRTH,&
     & PSTRCU,PSTRCV,PSTRDU,PSTRDV,PSTRTU,PSTRTV,&
     & PSTRMU,PSTRMV,PFRMH)

  ENDIF

ENDIF

!        3.3  Store the model trajectory at t-dt (leap-frog) or t (sl2tl).
!             ------------------------------------------------------------

IF (LTRAJPS) THEN
  IF (LTWOTL) THEN
    PTRAJ_PHYS%PQSSMF(KST:KEND)=PQS(KST:KEND)
    PTRAJ_PHYS%PTSMF(KST:KEND) =ZRRT0(KST:KEND)
    PTRAJ_PHYS%PSNSMF(KST:KEND)=ZSGF0(KST:KEND)
  ELSE
    CALL WRPHTRAJM(YDGEOMETRY,YDSIMPHL,KST,KEND,PTRAJ_PHYS,&
     & ZT9U,ZT9V,ZT9T,&
     & ZPQ9,ZPL9,ZPI9,ZT9SP)  

    PTRAJ_PHYS%PQSSMF(KST:KEND)=PQS(KST:KEND)
    PTRAJ_PHYS%PTSMF(KST:KEND) =ZRRT9(KST:KEND)
    PTRAJ_PHYS%PSNSMF(KST:KEND)=ZSGF9(KST:KEND)
  ENDIF
  IF (LPRTTRAJ.AND.PTRAJ_PHYS%LASTCHUNK) WRITE(NULOUT,*)'GREPTRAJ STORE TRAJ_PHYS in MF_PHYS'
ENDIF

!        3.4  Computation of tendencies T,u,v and Q.
!             --------------------------------------

IF (LSIMPH) THEN

  TSPHY = MAX(PDTPHY,1.0_JPRB)

  IF (LTWOTL) THEN
    CALL CPTENDSM (YDPHY,NPROMA,KST,KEND,NFLEVG,&
     & PDIFCQ,PDIFCS,PDIFTQ,PDIFTS,&
     & PFCCQL,PFCCQN,PFCSQL,PFCSQN,&
     & PFPLCL,PFPLCN,PFPLSL,PFPLSN,&
     & PFRSO,PFRTH,&
     & PSTRCU,PSTRCV,PSTRDU,PSTRDV,&
     & PSTRTU,PSTRTV,&
     & PSTRMU,PSTRMV,PFRMH,&
     & PXYB0(1,1,YYTXYB0_PHY%M_RDELP),PHIF0,&
     & ZT0U,ZT0V,ZT0T,ZPQ,&
     & PQS,ZRRT0,YDOROG%OROG,&
     & PTENDU,PTENDV,ZTENDH,ZTENDQ,&
     & PFHPCL,PFHPCN,PFHPSL,PFHPSN,&
     & PFHSCL,PFHSCN,PFHSSL,PFHSSN)  
  ELSE
    CALL CPTENDSM (YDPHY,NPROMA,KST,KEND,NFLEVG,&
     & PDIFCQ,PDIFCS,PDIFTQ,PDIFTS,&
     & PFCCQL,PFCCQN,PFCSQL,PFCSQN,&
     & PFPLCL,PFPLCN,PFPLSL,PFPLSN,&
     & PFRSO,PFRTH,&
     & PSTRCU,PSTRCV,PSTRDU,PSTRDV,&
     & PSTRTU,PSTRTV,&
     & PSTRMU,PSTRMV,PFRMH,&
     & PXYB9(1,1,YYTXYB9_PHY%M_RDELP),PHIF9,&
     & ZT9U,ZT9V,ZT9T,ZPQ9,&
     & PQS,ZRRT9,YDOROG%OROG,&
     & PTENDU,PTENDV,ZTENDH,ZTENDQ,&
     & PFHPCL,PFHPCN,PFHPSL,PFHPSN,&
     & PFHSCL,PFHSCN,PFHSSL,PFHSSN)  
  ENDIF

ENDIF

!        3.5  Computation of evolution of T, u, v and Q.
!             ------------------------------------------

IF (LSIMPH) THEN

  CALL CPUTQYS(YDGEOMETRY%YRDIMV,YDGMV,YGFL,YDPTRSLB1,&
   & YDSTOPH, YDPHY2, &
   & NPROMA,KST,KEND,NFLEVG,PDTPHY,IPGFL,&
   & ISLB1U9,ISLB1V9,ISLB1T9,ISLB1GFL9,&
   & ZTENDH, ZTENDQ, PTENDU, PTENDV,&
   & PFORCEU, PFORCEV, PFORCET, PFORCEQ, &
   & PRCP0(1,1,YYTRCP0%M_CP),PXYB0(1,1,YYTXYB0_PHY%M_DELP),ZT0T,ZT0U,ZT0V,&
   & PRCP9(1,1,YYTRCP9%M_CP),PXYB9(1,1,YYTXYB9_PHY%M_DELP),ZT9T,ZT9U,ZT9V,&
   & PB1, PGMVT1, PGFLT1,&
   & PFDIS)

ENDIF

!     ------------------------------------------------------------------

!*       4.    AROME  physics.
!              ---------------

IF (LMPA) THEN

  !      4.1  CALL APL_AROME
  !           --------------

  IPTR(:) = 0 ! means no fields defined in ZGFLTENDR at start ; > 0 means defined.

  ! * ZTENDR     
  IRR=1
  IPTR_CONT = IRR
  IF (YQ%LACTIVE) THEN
    IPTR(YQ%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YL%LACTIVE) THEN
    IPTR(YL%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YR%LACTIVE) THEN
    IPTR(YR%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YI%LACTIVE) THEN
    IPTR(YI%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YS%LACTIVE) THEN
    IPTR(YS%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YG%LACTIVE) THEN
    IPTR(YG%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YH%LACTIVE) THEN
    IPTR(YH%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   

  ! * ZTENDTKE   
  IF (YTKE%LACTIVE) THEN
    IPTR(YTKE%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
    IPTRTKE=IPTR(YTKE%MP1)
  ELSE
    IPTRTKE=0
  ENDIF
  ! * ZTENDEFB1 ZTENDEFB2 ZTENDEFB3
  IF (YEFB1%LACTIVE) THEN
    IPTR(YEFB1%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
    IEFB1 = IPTR(YEFB1%MP1)
  ELSE
    IEFB1 = 0
  ENDIF   
  IF (YEFB2%LACTIVE) THEN
    IPTR(YEFB2%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
    IEFB2 = IPTR(YEFB2%MP1)
  ELSE
    IEFB2 = 0
  ENDIF
  IF (YEFB3%LACTIVE) THEN   
    IPTR(YEFB3%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
    IEFB3 = IPTR(YEFB3%MP1)
  ELSE
    IEFB3 = 0
  ENDIF

  ! * ZTENDEXT
  IF (NGFL_EXT > 0) THEN
    DO JGFL=1,NGFL_EXT
      IF (YEXT(JGFL)%LACTIVE) THEN   
        IPTR(YEXT(JGFL)%MP1)=IPTR_CONT
        IPTR_CONT = IPTR_CONT+1   
      ENDIF
    ENDDO  
    IPTREXT=IPTR(YEXT(1)%MP1)  
  ELSE
    IPTREXT=0
  ENDIF

  ! * LIMA
  IF (NLIMA > 0) THEN
    DO JGFL=1,NLIMA
      IF (YLIMA(JGFL)%LACTIVE) THEN   
        IPTR(YLIMA(JGFL)%MP1)=IPTR_CONT
        IPTR_CONT = IPTR_CONT+1   
      ENDIF
    ENDDO  
    IPTRLIMA=IPTR(YLIMA(1)%MP1)  
  ELSE
    IPTRLIMA=0
  ENDIF
 

  ! If an incorrect address is used, then the initialization below will detect it :
  ZTENDGFLR(:,:,0)=HUGE(1._JPRB)
  
  IF (LMFSHAL .AND. CMF_UPDRAFT=='DUAL') THEN
    IMAXDRAFT=3
  ELSE
    IMAXDRAFT=0
  ENDIF

  IF(LTWOTL) THEN
    CALL APL_AROME(YDGEOMETRY,YDSURF,YDCFU,YDXFU,YDMODEL, KBL, KGPCOMP, KST  , KEND   , NPROMA ,&
     & 1      , NFLEVG , NSTEP  ,&
     & IMAXDRAFT, NTSSG, YDCFU%NFRRC, PDTPHY, LLXFUMSE,YDCSGEOM%RINDX,YDCSGEOM%RINDY,&
     & YDGSGEOM%GEMU,YDGSGEOM%GELAM,YDOROG%OROG,YDGSGEOM%GM,&
     & ZMU0,ZMU0LU,ZMU0M,ZMU0N,&
     & YDGSGEOM%GECLO,YDGSGEOM%GESLO,ZVCVC1,ZVFLSM,&
     & ZVASEA , ZVALAN ,ZVASOO , ZVADES ,&
     & ZVASUL , ZVAVOL ,&
     & PGP2DSDT, ZGP2DSPP, &
     & PHI0,PHIF0,PREPHY0,PREPHY0F,PXYB0(1,1,YYTXYB0_PHY%M_RDELP),PXYB0(1,1,YYTXYB0_PHY%M_DELP),&
     & ZT0T,ZPQ,&
     & PRCP0(1,1,YYTRCP0%M_CP),PRCP0(1,1,YYTRCP0%M_R),PXYB0(1,1,YYTXYB0_PHY%M_ALPH),PXYB0(1,1,YYTXYB0_PHY%M_LNPR),&
     & ZPL,ZPI,ZPR,ZPS,&
     & ZPG,ZPH,ZP1LIMA,ZPTKE,&
     & ZPEFB1,ZPEFB2,ZPEFB3,&
     & ZPSRC,ZP1EXT,&
     & ZT0U,ZT0V, PWT0,ZEDR,&
     & PFORCEU,PFORCEV,PFORCET,PFORCEQ,&
     & PGPAR,PEMTD,PEMTU,PTRSW,&
     & PGDEOSI,PGUEOSI,PGMU0,PGMU0_MIN,PGMU0_MAX,&
     & PGDEOTI,PGDEOTI2,PGUEOTI,PGUEOTI2,PGEOLT,PGEOXT,&
     & PGRPROX,PGMIXP,PGFLUXC,PGRSURF,&
     & PGRADH_PHY,&
     ! outputs
     & ZPLRAD1,ZPIRAD1,PRH,ZPA1,ZPSRC1,&
     & ZTENDT, ZTENDGFLR(1,1,IRR), PTENDU, PTENDV,ZTENDW, ZTENDGFLR(1,1,IPTRLIMA), ZTENDGFLR(1,1,IPTRTKE),&
     & ZTENDGFLR(1,1,IEFB1),ZTENDGFLR(1,1,IEFB2),ZTENDGFLR(1,1,IEFB3),&
     & ZTENDGFLR(1,1,IPTREXT),&
     & PFRTH, PFRSO, PFRTHDS, PFRSODS, PFRSOPS, PFRSDNI, PFRSOPT, PFRTHC, PFRSOC,&
     & ZVFALBF, ZVFEMISF,ZP1EZDIAG,&
     & PCLCH,PCLCL,PCLCM,PCLCT,PFPLSL,PFPLSN,PFPLSG,PFPLSHL,PSTRTU,PSTRTV,&
     & PFCS,PFCLL,PFCLN,PUCLS,PVCLS,PNUCLS,PNVCLS,PTCLS,PQCLS,PRHCLS,PUGST,PVGST,&
     & PFEVL,PFEVN,PCLPH,ZSGF1,ZSGR1,ZVDSUND,&
     & PDIAGH,PFLASH,ZSFOSFO1,PTPWCLS, ZDPRECIPS, ZDPRECIPS2,&
     & PVISICLD, PVISIHYDRO, PMXCLWC,YLPROCSET,YDDDH )
  ELSE    
    CALL APL_AROME(YDGEOMETRY,YDSURF,YDCFU,YDXFU,YDMODEL, KBL, KGPCOMP, KST  , KEND   , NPROMA ,&
     & 1      , NFLEVG , NSTEP  ,&
     & IMAXDRAFT, NTSSG, YDCFU%NFRRC, PDTPHY, LLXFUMSE,YDCSGEOM%RINDX,YDCSGEOM%RINDY,&
     & YDGSGEOM%GEMU,YDGSGEOM%GELAM,YDOROG%OROG,YDGSGEOM%GM,&
     & ZMU0,ZMU0LU,ZMU0M,ZMU0N,&
     & YDGSGEOM%GECLO,YDGSGEOM%GESLO,ZVCVC1,ZVFLSM,&
     & ZVASEA , ZVALAN ,ZVASOO , ZVADES ,&
     & ZVASUL , ZVAVOL ,&
     & PGP2DSDT, ZGP2DSPP, &
     & PHI9,PHIF9,PREPHY9,PREPHY9F,PXYB9(1,1,YYTXYB9_PHY%M_RDELP),PXYB9(1,1,YYTXYB9_PHY%M_DELP),&
     & ZT9T,ZPQ9,&
     & PRCP9(1,1,YYTRCP9%M_CP),PRCP9(1,1,YYTRCP9%M_R),PXYB9(1,1,YYTXYB9_PHY%M_ALPH),PXYB9(1,1,YYTXYB9_PHY%M_LNPR),&
     & ZPL9,ZPI9,ZPR9,ZPS9,&
     & ZPG9,ZPH9,ZP1LIMA9,ZPTKE9,&
     & ZPEFB19,ZPEFB29,ZPEFB39,&
     & ZPSRC9,ZP1EXT9,&
     & ZT9U,ZT9V, PWT9,ZEDR ,&
     & PFORCEU,PFORCEV,PFORCET,PFORCEQ,&
     & PGPAR,PEMTD,PEMTU,PTRSW,&
     & PGDEOSI,PGUEOSI,PGMU0,PGMU0_MIN,PGMU0_MAX,&
     & PGDEOTI,PGDEOTI2,PGUEOTI,PGUEOTI2,PGEOLT,PGEOXT,&
     & PGRPROX,PGMIXP,PGFLUXC,PGRSURF,&
     & PGRADH_PHY,&
     ! outputs
     & ZPLRAD1,ZPIRAD1,PRH,ZPA1,ZPSRC1,&
     & ZTENDT, ZTENDGFLR(1,1,IRR), PTENDU, PTENDV,ZTENDW,  ZTENDGFLR(1,1,IPTRLIMA), ZTENDGFLR(1,1,IPTRTKE),&
     & ZTENDGFLR(1,1,IEFB1),ZTENDGFLR(1,1,IEFB2),ZTENDGFLR(1,1,IEFB3),&
     & ZTENDGFLR(1,1,IPTREXT),&
     & PFRTH, PFRSO, PFRTHDS, PFRSODS, PFRSOPS, PFRSDNI, PFRSOPT, PFRTHC, PFRSOC,&
     & ZVFALBF, ZVFEMISF,ZP1EZDIAG,&
     & PCLCH,PCLCL,PCLCM,PCLCT,PFPLSL,PFPLSN,PFPLSG,PFPLSHL,PSTRTU,PSTRTV,&
     & PFCS,PFCLL,PFCLN,PUCLS,PVCLS,PNUCLS,PNVCLS,PTCLS,PQCLS,PRHCLS,PUGST,PVGST,&
     & PFEVL,PFEVN,PCLPH,ZSGF1,ZSGR1,ZVDSUND,&
     & PDIAGH,PFLASH,ZSFOSFO1, PTPWCLS, ZDPRECIPS, ZDPRECIPS2,&
     & PVISICLD, PVISIHYDRO, PMXCLWC,YLPROCSET,YDDDH )
  ENDIF 
  !Save surface temperature
  IF (LMSE.OR.LSFORCS) THEN
    IF (LLXFUMSE) THEN
      DO JROF=KST,KEND
        ZRRT0(JROF)=PGPAR(JROF,MVTS)
      ENDDO
    ELSE
      DO JROF=KST,KEND
        ZRRT1(JROF)=PGPAR(JROF,MVTS)
      ENDDO
    ENDIF 
  ENDIF  
  !      4.2  COMPUTE THE PHYS. TENDENCY FOR "T" AND "w"
  !           ------------------------------------------


  IF (LVERTFE.AND.LVFE_GWMPA) THEN
    ! * case LVFE_GWMPA not yet coded.
    !   (in this case ZGWT1 must be computed at full levels and
    !   not at half levels)
    CALL ABOR1(' MF_PHYS: case LVFE_GWMPA not yet coded if LMPA=T!')
  ENDIF

  ZDT = PDTPHY

  ZTENDD=0.0_JPRB

  ! * compute ZTT1:
  IF (LSLAG.AND.LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZTT1(JROF,JLEV)=ZT0T(JROF,JLEV)+ZDT*ZTENDT(JROF,JLEV)
      ENDDO
    ENDDO
  ELSE
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZTT1(JROF,JLEV)=ZT9T(JROF,JLEV)+ZDT*ZTENDT(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF

  ! * compute ZGWT1 = tendency of gw:
  IF (LNHDYN) THEN
    ! Valid for LVFE_GWMPA=F only; ZGWT1 assumed to be half level values.
    DO JLEV=1,NFLEVG-1
      DO JROF=KST,KEND
        ZGWT1(JROF,JLEV)=0.5_JPRB*RG*(ZTENDW(JROF,JLEV)+ZTENDW(JROF,JLEV+1))
      ENDDO
    ENDDO
    DO JROF=KST,KEND
      ZGWT1(JROF,NFLEVG)=0.0_JPRB
      ZGWT1(JROF,0)=0.0_JPRB
    ENDDO
  ENDIF

  ! * convert gw tendency in d tendency:
  IF(LNHDYN) THEN

    IF (LGWADV) THEN
      ZTENDD(KST:KEND,1:NFLEVG)=ZGWT1(KST:KEND,1:NFLEVG)
    ELSE

      ! * Provide the appropriate version of (RT) at t+dt for GNHGW2SVDAROME:
      IF (L_RDRY_VD) THEN
        ! Use Rd because "dver" is currently defined with Rd.
        ZRTT1(KST:KEND,1:NFLEVG)=RD*ZTT1(KST:KEND,1:NFLEVG)
      ELSE
        ! Use "moist R" because "dver" is defined with "moist R".
        ! Unfortunately, R(t+dt) is not yet available there, use R(t) instead.
        ! "Moist R" tendency is neglected in the below call to GNHGW2SVDAROME.
        DO JLEV=1,NFLEVG
          DO JROF=KST,KEND
            ZRTT1(JROF,JLEV)=PRCP0(JROF,JLEV,YYTRCP0%M_R)*ZTT1(JROF,JLEV)
          ENDDO
        ENDDO
      ENDIF

      ! * Do conversion:
      IF (LSLAG.AND.LTWOTL) THEN
        CALL GNHGW2SVDAROME(YDGEOMETRY,KST,KEND,PRE0F,PXYB0(1,1,YYTXYB0_PHY%M_LNPR),ZRTT1,PREPHY0F,ZGWT1,&
         & ZTENDD)  
      ELSE
        CALL GNHGW2SVDAROME(YDGEOMETRY,KST,KEND,PRE9F,PXYB9(1,1,YYTXYB9_PHY%M_LNPR),ZRTT1,PREPHY9F,ZGWT1,&
         & ZTENDD)  
      ENDIF

    ENDIF
  ELSE
    ZTENDD=0.0_JPRB
  ENDIF

  !      4.3  PUT THE TENDENCIES IN PB1/GFLT1/GMVT1.
  !           --------------------------------------


  IF ( LINTFLEX ) THEN
    IF (LTWOTL) THEN
      CALL CPTEND_FLEX( YDLDDH,YDMDDH,YGFL,YDPHY,NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & PXYB0(1,1,YYTXYB0_PHY%M_DELP) ,&
       & PXYB0(1,1,YYTXYB0_PHY%M_RDELP), PRCP0(1,1,YYTRCP0%M_CP),&
       & ZT0U,ZT0V,ZT0T,ZRRT0,&
       & PGFL,&
       & YLPROCSET,&
       & PTENDU , PTENDV , ZTENDH , ZTENDGFL,&
       & PFHSCL ,PFHSCN,PFHSSL,PFHSSN,&
       & PFHPCL ,PFHPCN,PFHPSL,PFHPSN,&
       & ZFHP   ,ZFP   ,  PFEPFP, PFCMPCQ, PFCMPSN, PFCMPSL,YDDDH )
    ELSE
      CALL CPTEND_FLEX( YDLDDH,YDMDDH,YGFL,YDPHY,NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & PXYB9(1,1,YYTXYB9_PHY%M_DELP) ,&
       & PXYB9(1,1,YYTXYB9_PHY%M_RDELP), PRCP9(1,1,YYTRCP9%M_CP),&
       & ZT9U,ZT9V,ZT9T,ZRRT9,&
       & PGFL,&
       & YLPROCSET,&
       & PTENDU , PTENDV , ZTENDH , ZTENDGFL,&
       & PFHSCL ,PFHSCN,PFHSSL,PFHSSN,&
       & PFHPCL ,PFHPCN,PFHPSL,PFHPSN,&
       & ZFHP   ,ZFP   ,  PFEPFP, PFCMPCQ, PFCMPSN, PFCMPSL,YDDDH )      
    ENDIF
    
    CALL CPUTQY(YDGEOMETRY%YRDIMV,YDGMV,YGFL,YDPTRSLB1,YDPHY,NPROMA,KST,KEND,NFLEVG,PDTPHY,IPGFL,&
     & ISLB1T9,ISLB1U9,ISLB1V9,ISLB1VD9,ISLB1GFL9,&
     & ZTENDH, ZTENDT, PTENDU, PTENDV, ZTENDU, ZTENDV, ZTENDD, ZTENDGFL,&
     & PRCP0(1,1,YYTRCP0%M_CP),PXYB0(1,1,YYTXYB0_PHY%M_DELP),ZT0T,ZT0U,ZT0V,&
     & PRCP9(1,1,YYTRCP9%M_CP),PXYB9(1,1,YYTXYB9_PHY%M_DELP),ZT9T,ZT9U,ZT9V,&
     & PB1, PGMVT1, PGFLT1,&
     & PFDIS)    
     
  ELSE
  
    ! start ZTENDGFLR at 1 because it is dimensionned (:,:,0:n)
    CALL CPUTQY_AROME(YDGEOMETRY%YRDIMV,YDGMV,YGFL,YDPTRSLB1,NPROMA,KST,KEND,NFLEVG,ZDT,IPGFL,IPTR,&
       & ISLB1U9,ISLB1V9,ISLB1T9,ISLB1GFL9,ISLB1VD9 ,&
       & ZTENDT, ZTENDGFLR(:,:,1), PTENDU, PTENDV, ZTENDD ,&
       & PB1, PGMVT1, PGFLT1)
  ENDIF
  
ENDIF

!     ------------------------------------------------------------------

!*       5.    Final calculations.
!              -------------------

! * Restore the initial value of some pseudo-historical surface buffers
!   if relevant.
IF (LL_SAVE_PHSURF) THEN
  IF(YSD_VV%YHV%LSET) ZVVHV(1:NPROMA)=ZHV(1:NPROMA)
  IF(YSD_VF%YZ0F%LSET) ZVFZ0F(1:NPROMA)=ZGZ0F(1:NPROMA)
  IF(YSD_VV%YZ0H%LSET) ZVVZ0H(1:NPROMA)=ZGZ0HF(1:NPROMA)
  IF(YSD_VH%YPBLH%LSET) ZVHPBLH(1:NPROMA)=ZPBLH(1:NPROMA)
  IF(YSD_VH%YSPSH%LSET) ZVHSPSH(1:NPROMA)=ZFHPS(1:NPROMA)
  IF(YSD_VH%YQSH%LSET)  ZVHQSH(1:NPROMA)=ZQSH(1:NPROMA)
  IF(YSD_VK%YUDGRO%LSET) ZVKUDGRO(1:NPROMA)=ZUDGRO(1:NPROMA)
  IF(LCVPRO.OR.LGPCMT) THEN
    ZPUAL(:,:)=ZUDAL(:,:)
    ZPUOM(:,:)=ZUDOM(:,:)
    IF(LCDDPRO) THEN
      ZPDAL(:,:)=ZDDAL(:,:)
      ZPDOM(:,:)=ZDDOM(:,:)
    ENDIF
  ENDIF
  IF(YUNEBH%LACTIVE) ZPUNEBH(:,:)=ZUNEBH(:,:)
  IF(YUEN%LACTIVE)   ZPUEN(:,:)=ZENTCH(:,:)
ENDIF

! Store horizontal exchange coefficients (3D turbulence) to SL2 buffers
IF (L3DTURB) THEN
  DO JLEV=1,NFLEVG
    PB2(KST:KEND,MSLB2KAPPAM+JLEV-1)=ZKUROV_H(KST:KEND,JLEV)
    PB2(KST:KEND,MSLB2KAPPAH+JLEV-1)=ZKTROV_H(KST:KEND,JLEV)
  ENDDO
ENDIF

!--------------------------------------------------------------------
! BAYRAD
! Fill convective hydrometeors mixing ratio in GFL
!--------------------------------------------------------------------
IF((.NOT.LGPCMT).AND.(.NOT.LAROME)) THEN
   IF(YRCONV%LACTIVE) THEN
     ZPRCONV1(:,:) = ZQRCONV(:,:)
   ENDIF
   IF(YSCONV%LACTIVE) THEN
     ZPSCONV1(:,:) = ZQSCONV(:,:)
   ENDIF
ENDIF



!-------------------------------------------------
! Extract Single Column Model profiles from 3D run or 
! write LFA file for MUSC (1D model)
!-------------------------------------------------
IF(LGSCM.OR.LMUSCLFA) THEN
  IF (LAROME) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        PNEB(JROF,JLEV)=ZPA1(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
  CALL WRITEPHYSIO(YDGEOMETRY,YDSURF,YDDPHY,YDRIP,YDMODEL%YRML_PHY_MF,&
       & KEND,&
       & KST, KGL1, KGL2, KSTGLO,&
       & NSTEP  , NTSSG  , YSP_SBD%NLEVS   ,&
       & YDGSGEOM%GELAM,YDGSGEOM%GEMU,&
       & YDGSGEOM%GM, ZMU0,YDOROG%OROG,POROGL,POROGM,YDGSGEOM%RCORI,YDCSGEOM%RATATH,YDCSGEOM%RATATX,&
       & PHI0  , PRE0  , PHIF0 , PRE0F , PXYB0(1,1,YYTXYB0_PHY%M_ALPH), PXYB0(1,1,YYTXYB0_PHY%M_DELP),&
       & PXYB0(1,1,YYTXYB0_PHY%M_LNPR), PXYB0(1,1,YYTXYB0_PHY%M_RDELP),&
       & ZVFLSM, ZVVARG, ZVVSAB,&
       & ZVVD2, ZVVIVEG, ZVVLAI, PCT,&
       & ZVVALV, PALBDG, ZVFALBF, ZVFALBSF, ZSGT1,&
       & ZVFEMISF,&
       & ZVFZ0F, ZVVZ0H, ZVFZ0RLF,&
       & ZVFGETRL, ZVFVRLAN, ZVFVRLDI, ZVVRSMIN,&
       & ZVFVEG, ZVVHV, ZVASEA, ZVALAN,&
       & ZVASOO, ZVADES, ZVCVC1,&
       & ZT0SP,ZT0SPL,ZT0SPM,&
       & ZT0T,ZT0TL,ZT0TM,&
       & ZPQ,ZPQL,ZPQM,&
       & ZPI,ZPL, ZPS,ZPR,ZPG,ZPTKE,&
       & ZPEFB1,ZPEFB2,ZPEFB3,&
       & PRCP0(1,1,YYTRCP0%M_CP), PRCP0(1,1,YYTRCP0%M_R),&
       & ZT0U,ZT0V,ZT0VOR,ZT0DIV, ZCVGQ, ZLCVQ,&
       & ZRRT9, ZSBT9, ZRRFC9, ZRRW9,&
       & ZRRIC9, ZSBQ9, ZSBTL9, ZSGF9,&
       & ZSGA0, ZSGR0, PCTY0(1,1,YYTCTY0%M_VVEL),&
       & PEMTD , PEMTU ,PTRSW ,YDGSGEOM%GECLO,YDGSGEOM%GESLO,&
       & PDIFCQ, PDIFCQI, PDIFCQL, PDIFCS, PDIFTQ, PDIFTQI, PDIFTQL,&
       & PDIFTS, PFCCQL, PFCCQN, PFCSQL, PFCSQN, PFCQNG, PFCQING,&
       & PFCQLNG, PFCQRNG, PFCQSNG, PFPLCL, PFPLCN, PFPLSL, PFPLSN, PFPLSG, PFPLSHL,&
       & PFPFPSL, PFPFPSN, PFPFPCL, PFPFPCN, PFPEVPSL,PFPEVPSN,PFPEVPCL, PFPEVPCN,&
       & ZFTKE, ZFEFB1, ZFEFB2, ZFEFB3, PFRSO, PFRSOC, PFRTH, PFRTHC, PFRSOLU,&
       & PSTRCU, PSTRCV, PSTRDU, PSTRDV, PSTRTU, PSTRTV, PSTRMU, PSTRMV,&
       & PFRMH, ZFRMQ, PFCHOZ, PNEB, PQICE, PQLI, PRH,&
       & PFCS, PFCLL, PFCLN, PFEVL, PFEVN, ZFEVI, PFEVV, PFTR, PFLWSP, PFONTE,&
       & PFGEL, PFGELS, PFCHSP, PFRSODS, PFRSOPS, PFRSOPT, PFRTHDS,&
       & ZCD, ZCDN, ZCH, ZC1, ZC2, ZEMIS, PGZ0   , PGZ0H  , ZNEIJ  , ZVEG,&
       & ZCPS, ZLHS, ZRS, ZLH, ZLSCPE, ZQSAT, ZQW, ZTW,&
       & PQS, ZQSATS, PRUISL, PRUISP, PRUISS,&
       & PUCLS, PVCLS, PTCLS, PQCLS, PRHCLS,&
       & PCLCT, PCLCH, PCLCM, PCLCL, PCLCC,&
       & PCAPE  , PCTOP, ICLPH  , PCLPH  , PUGST  , PVGST,&
       & ZFPLCH, ZFPLSH,&
       & ZVHTCCH,  ZVHSCCH, ZVHBCCH, ZVHPBLH )
ENDIF

IF (LEDR) THEN
  ZDIXEDR(:,:)=ZEDR(:,:)
ENDIF

IF (LDPRECIPS) THEN
  PSD_XP(KST:KEND,NDTPRECCUR,YSD_XP%YPRECIP%MP)=ZDPRECIPS(KST:KEND,NDTPRECCUR)
ENDIF

IF (LDPRECIPS2) THEN
  PSD_XP2(KST:KEND,NDTPRECCUR2,YSD_XP2%YPRECIP2%MP)=ZDPRECIPS2(KST:KEND,NDTPRECCUR2)
ENDIF

! Restore Tt and grad(Tt) for NHQE model.
IF (LNHQE) THEN
  ! At instant t (with the derivatives):
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZT0T(JROF,JLEV)=ZTT0_SAVE(JROF,JLEV)
      ZT0TL(JROF,JLEV)=ZTT0L_SAVE(JROF,JLEV)
      ZT0TM(JROF,JLEV)=ZTT0M_SAVE(JROF,JLEV)
    ENDDO
  ENDDO
  ! At instant t-dt for leap-frog advections (without the derivatives):
  IF (.NOT.LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZT9T(JROF,JLEV)=ZTT9_SAVE(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!       6. destructor for procset
IF (LINTFLEX) CALL CLEANINTPROCSET(YLPROCSET)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('MF_PHYS',1,ZHOOK_HANDLE)
END SUBROUTINE MF_PHYS
