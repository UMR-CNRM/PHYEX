#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE CPG_GP(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_TND, YDCPG_MISC, YDCPG_GPAR, YDCPG_DYN0, &
& YDCPG_DYN9, YDMF_PHYS_SURF, YDVARS, YDMODEL, YDCPG_SL2, YDPGTERM, YDDDH)

!**** *CPG_GP* - Grid point calculations:
!                initial part of not lagged grid-point calculations.

!     Purpose.
!     --------
!           Grid point calculations:
!           initial part of not lagged grid-point calculations.
!           - get data in buffers.
!           - multiply p-order horizontal derivatives by M**p.
!           - second part of the temporal filter.
!           - grid-point calculations for nudging.
!           - calls some GP... routines to initialise some auxiliary variables.
!           - sets-up and PB2.

!           Abbreviation "vwv" stands for "vertical wind variable".

!**   Interface.
!     ----------
!        *CALL* *CPG_GP(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        LD_DFISTEP   : 'D' -> DFI computations
!        KST       : first element of work.
!        KEND      : last element of work.
!        KSTC,KENDC: the same as KST,KEND but including zone C for ALADIN.
!        KBL       : block number.
!        KSTGLO    : global offset.
!        LDLFSTEP  : .T.: first time-step?
!        LDLDIAB   : .T. if complete physics is activated and predictor step.
!        PDT       : For a leap-frog scheme (three time level scheme):
!                     'dt' at the first time-step, '2 dt' otherwise.
!                    For a 2TL SL scheme: timestep 'dt'.
!        PTE       : 1. or 0. according to different configurations.
!        KIBL      : index into YRGSGEOM/YRCSGEOM types in YDGEOMETRY

!     INPUT/OUTPUT:
!     -------------
!        PGFL      : unified_treatment grid-point fields at t
!        PGMV      : upper air GMV variables at time t and t-dt.
!        PGMVS     : surface GMV variables at time t and t-dt.
!        YDPGTERM  : pressure gradient quantities, only used if LPGFSF=T
!        PGMVTNDSL : GMV(t+dt,F)-GMV(t or t-dt,O) for DDH
!        PGFLTNDSL : GFL(t+dt,F)-GFL(t or t-dt,O) for DDH

!     OUTPUT:
!     -------
!        POROGL    : zonal component of "grad(surf orography)"
!        POROGM    : meridian component of "grad(surf orography)"
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE0L     : zonal component of "grad prehyds" at t.
!        PRE0M     : meridian component of "grad prehyds" at t.
!        PRE0F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PNHPRE0F  : "pre" at full levels (time t).
!        PNHPRE0H  : "pre" at half levels (time t).
!        PXYB0     : contains pressure depth, "delta", "alpha" at t.
!        PUVH0     : horizontal wind at time t at half levels.
!        PHIF0     : geopotential height "gz" at t at full levels.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        PGWFT0    : [Gw] at full layers at t.
!        PKENE0    : kinetic energy at t.
!        PRE9      : hydrostatic pressure "prehyd" at half levels at t-dt.
!        PRE9F     : hydrostatic pressure "prehyd" at full levels at time t-dt.
!        PNHPRE9F  : "pre" at full levels (time t-dt).
!        PNHPRE9H  : "pre" at half levels (time t-dt).
!        PXYB9     : contains pressure depth, "delta", "alpha" at t-dt.
!        PHI9      : geopotential height "gz" at t-dt at half levels.
!        PHIF9     : geopotential height "gz" at t-dt at full levels.
!        PRCP9     : contains "cp", "R" and "Kap=R/Cp" at t-dt.
!        PGWFT9    : [Gw] at full levels at t-dt.
!        PGPAR     : surface fields for AROME.
!        PB2       : "SLB2" buffer.
!        PGMVT1    : upper air GMV variables at t+dt.
!        PGMVT1S   : surface GMV variables at t+dt.
!        PGFLT1    : GFL variables at t+dt.
!        PKOZO     : fields for photochemistery of ozon.
!        PQS       : specific humidity at surface level.
!        PQICE     : specific humidity of solid water for radiation.
!        PQLI      : specific humidity of liquid water for radiation.
!        PQRAIN    : specific humidity of rain for radiation.
!        PQSNOW    : specific humidity of snow for radiation.
!        PQGRAUPEL : specific humidity of graupel for radiation.
!        PATND     : adiabatic Lagrangian tendencies.
!        PDBBC     : [D (Gw)_surf / Dt]_adiab.
!        PRDPHI    : HYD: not used.
!                    NHEE: contains pre/(R T prehyd [Delta log(prehyd)]) at t.
!                    NHQE: contains 1/(R Tt [Delta log(prehyd)]) at t.
!                    "R" is the version of R (may be Rdry or Rmoist) used in the definition of vertical divergence "dver".
!        PGWT0     : [Gw] at t (LGWADV=T only).
!        PGWT9     : [Gw] at t-dt (LGWADV=T only).
!                    for PGWT0 and PGWT9:
!                    * half level 0 to nflevg-1 values if LVFE_GW=F.
!                    * full level 1 to nflevg values if LVFE_GW=T.
!        PGWS      : [Gw]_surf at t (LRDBBC only).
!        PGWFL     : zonal comp grad(Gw) at full level at t.
!        PGWFM     : meridian comp grad(Gw) at full level at t.
!        PNHXT0    : term 'NHX' at t.
!        PNHXT9    : term 'NHX' at t-dt.
!        PNHYT0    : term 'NHY' at t.    (NHY = gW - gw)
!        PNHYT9    : term 'NHY' at t-dt. (NHY = gW - gw)
!        PQCHA0L   : zonal comp grad(log(pre/prehyd)).
!        PQCHA0M   : merid comp grad(log(pre/prehyd)).
!        PEXTRA    : additional quantity for diagnostics.
!        YDDDH     : diagnostic superstructure

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation vol 2 ch 1 and vol 3 ch 6

!     Author.
!     -------
!        K. YESSAD, after parts 1, 2 and 3 of old CPG.
!        Original : 16-08-2001

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Sep 2008): prune enhanced diffusion (lfrein).
!   F. Bouyssel (Nov 2008): Removal of LPROCLD protection
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (March 2009): correct false comments for LRWSDLG=T
!   K. Yessad (Aug 2009): remove LPC_OLD in TL and AD codes.
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP.
!   K. Yessad (Nov 2009): prune lpc_old.
!   K. Yessad (Dec 2009): LRWSDLW,LRWSDLR,LRWSDLG=T,T,T in NH model for LGWADV=F.
!   A. Alias  (Mar 2011) No nudging of SST when LSME=.T.
!   K. Yessad (Nov 2011): various contributions.
!   N. Wedi   (Nov 2011): add LGRADSP
!   M. Ahlgrimm  31-Oct-2011 add rain, snow and PEXTRA to DDH output
!   K. Yessad (Nov 2012): simplify testings.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!   K. Yessad (July 2014): Move some variables.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   R. Roehrig (Sept 2018): add argument to CP_FORCING + add possibility to
!                           impose (time-evolving) surface pressure in MUSC
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!   F. Voitus (Dec 2019): NVDVAR=5.
!   P. Smolikova (Sep 2020): Remove obsolete calculations in hyd.model.
!   I. Polichtchouk & F. Vana (Oct-2020): LPGFSF with LPGREUSE
!   I. Polichtchouk (Jul 2021): Add LSACC option
!   D. Nemec (July 2021): Add PQGRAUPEL.
! End Modifications
!     ------------------------------------------------------------------

USE TYPE_MODEL              , ONLY : MODEL
USE GEOMETRY_MOD            , ONLY : GEOMETRY
USE CPG_TYPE_MOD            , ONLY : CPG_DYN_TYPE, CPG_GPAR_TYPE, CPG_MISC_TYPE, CPG_SL2_TYPE, CPG_TND_TYPE
USE CPG_OPTS_TYPE_MOD       , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD, ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD     , ONLY : FIELD_VARIABLES
USE DDH_MIX                 , ONLY : TYP_DDH
USE PARKIND1                , ONLY : JPIM, JPRB
USE YOMHOOK                 , ONLY : DR_HOOK, LHOOK, JPHOOK
USE INTDYN_MOD              , ONLY : TPG_TYPE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)          ,INTENT(IN)               :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE)     ,INTENT(IN)               :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE)     ,INTENT(IN)               :: YDCPG_OPTS
TYPE(CPG_TND_TYPE)      ,INTENT(INOUT)            :: YDCPG_TND
TYPE(CPG_MISC_TYPE)     ,INTENT(INOUT)            :: YDCPG_MISC
TYPE(CPG_GPAR_TYPE)     ,INTENT(INOUT)            :: YDCPG_GPAR
TYPE(CPG_DYN_TYPE)      ,INTENT(INOUT)            :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE)      ,INTENT(INOUT)            :: YDCPG_DYN9
TYPE(MF_PHYS_SURF_TYPE) ,INTENT(INOUT)            :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES)   ,INTENT(INOUT)            :: YDVARS
TYPE(MODEL)             ,INTENT(IN)               :: YDMODEL
TYPE(CPG_SL2_TYPE)      ,INTENT(INOUT)            :: YDCPG_SL2
TYPE(TPG_TYPE)          ,INTENT(INOUT)            :: YDPGTERM
TYPE(TYP_DDH)           ,INTENT(INOUT)            :: YDDDH

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "cp_forcing.intfb.h"
#include "cp_forcing_ps.intfb.h"
#include "gpinislb_part1_expl.intfb.h"
#include "gpinislb_part2_expl.intfb.h"
#include "gpinislb_part3_expl.intfb.h"
#include "gpinozst.intfb.h"
#include "gpmpfc_expl.intfb.h"
#include "gpnspng_expl.intfb.h"
#include "gp_spv.intfb.h"
#include "gptf2_expl_2tl.intfb.h"
#include "gptf2_expl_3tl_part1.intfb.h"
#include "gptf2_expl_3tl_part2.intfb.h"
#include "updsst.intfb.h"
#include "cpg_gp_hyd.intfb.h"
#include "cpg_gp_nhee.intfb.h"
#include "cpg_gp_nhqe.intfb.h"
#include "cpg_gp_sacc.intfb.h"

!     ------------------------------------------------------------------

REAL (KIND=JPRB) :: Z_GDW_T0 (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL (KIND=JPRB) :: Z_GWHT_T0 (YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL (KIND=JPRB) :: Z_OROGLL_T0 (YDGEOMETRY%YRDIM%NPROMA)
REAL (KIND=JPRB) :: Z_OROGLM_T0 (YDGEOMETRY%YRDIM%NPROMA)
REAL (KIND=JPRB) :: Z_OROGMM_T0 (YDGEOMETRY%YRDIM%NPROMA)

INTEGER(KIND=JPIM) :: JROF,IFLAG, JGFL
LOGICAL :: LLSTR, LLGPXX, LLUVH
REAL(KIND=JPRB) :: ZEPS


REAL(KIND=JPRB) :: ZGM2

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_GP', 0, ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR, &
& YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY, YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, YGFL=>YDMODEL%YRML_GCONF%YGFL,                      &
 & YDCVER=>YDGEOMETRY%YRCVER, &
& YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, YDDYN=>YDMODEL%YRML_DYN%YRDYN,                     &
& YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL, YDPHLC=>YDMODEL%YRML_PHY_SLIN%YRPHLC, YDDYNA=>YDMODEL%YRML_DYN%YRDYNA)

ASSOCIATE(YQ=>YGFL%YQ, NPROMA=>YDDIM%NPROMA, NFLEVG=>YDDIMV%NFLEVG, NTOZ1D=>YDDPHY%NTOZ1D, NTOZ2D=>YDDPHY%NTOZ2D, &
& NTOZ3D=>YDDPHY%NTOZ3D, NCURRENT_ITER=>YDDYN%NCURRENT_ITER, LEO3CH=>YDEPHY%LEO3CH, RSTRET=>YDGEM%RSTRET,         &
& LSPHLC=>YDPHLC%LSPHLC, LSIMPH=>YDSIMPHL%LSIMPH, LMSE=>YDARPHY%LMSE, LMPA=>YDARPHY%LMPA, LMPHYS=>YDPHY%LMPHYS    &
& )


!     ------------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)
LLSTR=(ABS(RSTRET-1.0_JPRB)>ZEPS)

!     ------------------------------------------------------------------

!*       1.    INTERFACE TO GLOBAL ARRAYS/WORK FILES.
!              --------------------------------------

!*     1.3   SPONGE AT THE TOP OF THE MODEL (ACADEMIC SIMULATIONS).

! new sponge for 2D and 3D models (grid-point GFL only).
IF (YDMODEL%YRML_DYN%YRSPNG%LNSPONGE) THEN
  CALL GPNSPNG_EXPL (YDMODEL%YRML_DYN%YRSPNG, NPROMA, NFLEVG, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, &
  & YDVARS)
ENDIF

!*     1.4   PART 2 OF TIME FILTER.

IF(NCURRENT_ITER == 0) THEN
  IF (YDDYNA%LTWOTL) THEN
    CALL GPTF2_EXPL_2TL (YDGEOMETRY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%LFSTEP, P0U=YDVARS%U%T0, &
    & P0V=YDVARS%V%T0, P9U=YDVARS%U%T9, P9V=YDVARS%V%T9)
  ELSE
    CALL GPTF2_EXPL_3TL_PART1 (YDGEOMETRY, YDMODEL%YRML_GCONF, YDDYN, YDDYNA, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,      &
    & YDCPG_OPTS%LFSTEP, P0SP=YDVARS%SP%T0, P0SPL=YDVARS%SP%DL,         &
    & P0SPM=YDVARS%SP%DM, P9SP=YDVARS%SP%T9, P9SPL=YDVARS%SP%DL9, P9SPM=YDVARS%SP%DM9, P0DIV=YDVARS%DIV%T0,            &
    & P0NHX=YDVARS%NHX%T0, P0SPD=YDVARS%SPD%T0, P0SPDL=YDVARS%SPD%DL, P0SPDM=YDVARS%SPD%DM, P0SVD=YDVARS%SVD%T0,       &
    & P0SVDL=YDVARS%SVD%DL, P0SVDM=YDVARS%SVD%DM, P0T=YDVARS%T%T0, P0TL=YDVARS%T%DL, P0TM=YDVARS%T%DM,                 &
    & P0U=YDVARS%U%T0, P0V=YDVARS%V%T0, P9DIV=YDVARS%DIV%T9, P9NHX=YDVARS%NHX%T9, P9SPD=YDVARS%SPD%T9,                 &
    & P9SPDL=YDVARS%SPD%DL9, P9SPDM=YDVARS%SPD%DM9, P9SVD=YDVARS%SVD%T9, P9SVDL=YDVARS%SVD%DL9, P9SVDM=YDVARS%SVD%DM9, &
    & P9T=YDVARS%T%T9, P9TL=YDVARS%T%DL9, P9TM=YDVARS%T%DM9, P9U=YDVARS%U%T9, P9V=YDVARS%V%T9)
    CALL GPTF2_EXPL_3TL_PART2 (YDGEOMETRY, YDDYN, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%LFSTEP, &
    & YDVARS)
  ENDIF
ENDIF

!*     1.5   MAP FACTOR AND SURFACE PRESSURE VARIABLES.

YDCPG_DYN0%OROGL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDVARS%GEOMETRY%OROGL%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*YDVARS%GEOMETRY%GM%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
YDCPG_DYN0%OROGM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDVARS%GEOMETRY%OROGM%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*YDVARS%GEOMETRY%GM%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
IF(YDDYNA%LNHDYN) THEN
  DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZGM2 = YDVARS%GEOMETRY%GM%T0(JROF)**2
    Z_OROGLL_T0(JROF)=YDVARS%GEOMETRY%OROGLL%T0(JROF)*ZGM2
    Z_OROGMM_T0(JROF)=YDVARS%GEOMETRY%OROGMM%T0(JROF)*ZGM2
    Z_OROGLM_T0(JROF)=YDVARS%GEOMETRY%OROGLM%T0(JROF)*ZGM2
  ENDDO
ENDIF

IF(LLSTR.AND.(.NOT.YDDYNA%LPGFSF)) THEN
  IFLAG=0
  CALL GPMPFC_EXPL(YDVARS, YDMODEL%YRML_GCONF, YDDYN, YDDYNA, NPROMA, NFLEVG, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, IFLAG, YDVARS%GEOMETRY%GM%T0)
ENDIF

CALL GP_SPV(YDGEOMETRY, YDDYN, YDDYNA, YDSIMPHL,  .FALSE., YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDVARS%SP%T0, &
& YDVARS%SP%DL, YDVARS%SP%DM, YDVARS%SP%T9,  YDVARS%SP%DL9, YDVARS%SP%DM9, YDCPG_DYN0%PRE, YDCPG_DYN0%PREL,  &
& YDCPG_DYN0%PREM, YDCPG_DYN9%PRE, YDCPG_DYN9%PREL,  YDCPG_DYN9%PREM)

IF (YDCPG_OPTS%LSFORC .AND. YDCPG_OPTS%LSPS_FRC) THEN
  CALL CP_FORCING_PS(YDGEOMETRY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_DYN0%PRE)
ENDIF

IF (NCURRENT_ITER == 0) THEN

!*     1.6   SURFACE VARIABLES.

  IF (YDCPG_OPTS%LDIAB.OR.LSPHLC.OR.LSIMPH) THEN
    YDCPG_MISC%QS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.0_JPRB
    IF(YDCPG_OPTS%LFSTEP .AND. .NOT. YDDYNA%LTWOTL) THEN
      CALL YDMF_PHYS_SURF%SET9TO0 ()
    ENDIF
  ENDIF

!*     1.8   OZON PHYSICO-CHEMICAL PROPERTIES AND SUBGRID SURFACE TEMPERATURE.

  IF (YDCPG_OPTS%LDIAB) THEN
    IF(NTOZ3D > 0.OR.NTOZ2D > 0.OR.NTOZ1D > 0.AND.(.NOT.LEO3CH)) THEN
      CALL GPINOZST(YDGEOMETRY, YDMODEL%YRML_CHEM%YROZO, YDDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, &
      & YDCPG_BNDS%KSTGLO, YDCPG_MISC%KOZO)
    ENDIF
  ENDIF

!*     1.9  NUDGING.

 IF(YDCPG_OPTS%LNUDG.AND.(.NOT.LMSE)) THEN
    CALL UPDSST(YDMODEL%YRML_AOC%YRMCC, YDMODEL%YRML_PHY_MF%YRPHY1, NPROMA, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,                &
    & YDCPG_OPTS%NFNUDG, YDCPG_OPTS%YRSURF_DIMS%YSP_SBD%NLEVS, YDCPG_BNDS%KBL, YDCPG_OPTS%XPNUDG, YDMF_PHYS_SURF%GSP_RR%PT_T0, &
    & YDMF_PHYS_SURF%GSP_SB%PT_T0, YDMF_PHYS_SURF%GSP_RR%PW_T0, YDMF_PHYS_SURF%GSP_SB%PQ_T0, YDMF_PHYS_SURF%GSP_SG%PF_T0,      &
    & YDMF_PHYS_SURF%GSP_SG%PA_T0, YDMF_PHYS_SURF%GSP_SG%PR_T0, YDMF_PHYS_SURF%GSP_RR%PIC_T0, YDMF_PHYS_SURF%GSP_SB%PTL_T0,    &
    & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VF%PALBF, YDMF_PHYS_SURF%GSD_VF%PEMISF, YDMF_PHYS_SURF%GSD_VF%PZ0F,       &
    & YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VV%PD2, YDMF_PHYS_SURF%GSD_VV%PIVEG)
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       2.    INTERFACE TO GLOBAL ARRAYS/WORK FILES.
!              --------------------------------------

!     ------------------------------------------------------------------

!*       3.    INITIALISE AUXILIARY VARIABLES.
!              -------------------------------

IF (YDDYNA%LNHQE) THEN
  ! * NHQE model:
  LLGPXX=.NOT.(YDDYNA%LSLAG.AND.(YDDYNA%NVDVAR == 4).AND.(YDDYNA%ND4SYS==2))
  CALL CPG_GP_NHQE(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDVARS, YDCPG_DYN0, YDCPG_DYN9, YDMODEL, &
  & LLGPXX, YDCPG_OPTS%LDIAB, LMPA, Z_OROGLL_T0, Z_OROGMM_T0, Z_OROGLM_T0, Z_GDW_T0, Z_GWHT_T0, &
  & YDCPG_TND)

ELSEIF (YDDYNA%LNHEE) THEN
  ! * NHEE model:
  LLGPXX=.NOT.(YDDYNA%LSLAG.AND.(YDDYNA%NVDVAR == 4 .OR. YDDYNA%NVDVAR == 5).AND.(YDDYNA%ND4SYS==2))
  CALL CPG_GP_NHEE(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDVARS, YDCPG_DYN0, YDCPG_DYN9, YDMODEL, &
  & LLGPXX, YDCPG_OPTS%LDIAB, LMPA, Z_OROGLL_T0, Z_OROGMM_T0, Z_OROGLM_T0, Z_GDW_T0, Z_GWHT_T0, &
  & YDCPG_TND)

ELSEIF (YDDYNA%LSACC) THEN
  ! * Shallow-atmosphere model with complete Coriolis of Tort & Dubos (2013)
  ! * LLPGX=t is needed because we need [gw] term 
  LLGPXX=.TRUE.
  CALL CPG_GP_SACC(YDGEOMETRY,YDMODEL,&
   !---------------------------------------------------------------------
   ! - INPUT .
   & YDCPG_BNDS, YDCPG_OPTS, YDCPG_TND, &
   & LLGPXX,YDCPG_OPTS%LDIAB,LMPA,&
   !---------------------------------------------------------------------
   ! - INPUT/OUTPUT .
   & YDCPG_DYN0, YDCPG_DYN9, &
   & YDVARS,&
   !---------------------------------------------------------------------
   ! - OUTPUT .
   & YDPGTERM%PHI0H,YDPGTERM%PHI0FL,YDPGTERM%PHI0FM,&
   & Z_GWHT_T0,&
   & YDPGTERM%RT0L,YDPGTERM%RT0M)
ELSE
  ! * Hydrostatic model:
  IF (ASSOCIATED(YDCPG_DYN0%NHX)) YDCPG_DYN0%NHX(:,:)=0.0_JPRB
  IF (ASSOCIATED(YDCPG_DYN0%GWFT)) YDCPG_DYN0%GWFT(:,:)=0.0_JPRB
  IF (.NOT. YDDYNA%LTWOTL) THEN
    IF (ASSOCIATED(YDCPG_DYN9%GWFT)) YDCPG_DYN9%GWFT(:,:)=0.0_JPRB
    IF (ASSOCIATED(YDCPG_DYN9%NHX)) YDCPG_DYN9%NHX(:,:)=0.0_JPRB
  ENDIF
  IF(YDCVER%LVERTFE.AND..NOT.YDCVER%LVFE_ECMWF .OR. .NOT.((LMPHYS.OR.LSIMPH) .AND. NCURRENT_ITER == 0)) THEN
    LLGPXX=.FALSE.
  ELSE
    LLGPXX=.TRUE.
  ENDIF
  LLUVH=LLGPXX.OR.YDMODEL%YRML_DYN%YRDYNA%LSLHD
  CALL CPG_GP_HYD(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDVARS, YDCPG_DYN0, YDCPG_DYN9, YDMODEL, LLUVH, &
  & LLGPXX, YDCPG_OPTS%LDIAB, LMPA, Z_OROGLL_T0, Z_OROGMM_T0, Z_OROGLM_T0, Z_GDW_T0, Z_GWHT_T0,       &
  & YDPGTERM, YDCPG_TND)
ENDIF

!     ------------------------------------------------------------------

!*       4.    TIME DIMENSION.
!              ---------------

CALL GPINISLB_PART1_EXPL (YDGEOMETRY, YGFL, YDDYN, YDDYNA, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%ZTE, &
& YDVARS%U%T9, YDVARS%V%T9,  YDVARS%T%T9, YDVARS%SPD%T9, YDVARS%SVD%T9, YDCPG_DYN9%NHX, YDVARS%SP%T9,          &
& YDVARS%U%T1, YDVARS%V%T1,  YDVARS%T%T1, YDVARS%SPD%T1, YDVARS%SVD%T1, YDVARS%NHX%T1, YDVARS%SP%T1,           &
& YDCPG_MISC%QICE, YDCPG_MISC%QLI,  YDCPG_MISC%QRAIN, YDCPG_MISC%QSNOW, YDCPG_MISC%QGRAUPEL, PL0=YDVARS%L%T0,  &
& PL9=YDVARS%L%T9, PI0=YDVARS%I%T0, PI9=YDVARS%I%T9,  PR0=YDVARS%R%T0, PR9=YDVARS%R%T9, PS0=YDVARS%S%T0,       &
& PS9=YDVARS%S%T9, PG0=YDVARS%G%T0, PG9=YDVARS%G%T9) 

CALL GPINISLB_PART2_EXPL (YDDYN, YDDYNA, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%ZTE, YDVARS)

CALL GPINISLB_PART3_EXPL (YDGEOMETRY, YDDYNA, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_DYN0%CTY%VVEL (:, 1:),   &
& YDCPG_DYN0%PREF, YDCPG_SL2%VVEL (:, 1:),  YDCPG_SL2%GWF (:, 1:), YDCPG_SL2%GDW (:, 1:), YDCPG_SL2%GWS (:, 1:), &
& YDCPG_DYN0%GWFT, Z_GDW_T0, Z_GWHT_T0(:, NFLEVG)    )


IF (NCURRENT_ITER == 0) THEN
  IF (YDCPG_OPTS%LDIAB.OR.LSIMPH) THEN
    CALL YDMF_PHYS_SURF%SET1TO9 ()
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       5.    TENDENCIES FOR 1D MODEL.
!              ------------------------

! Add large scale forcing for SCUM (Single Column Unified Model)
! - incrementation of PATND for GMV variables.
! - fills ZATND_Q for humidity.

IF (YDCPG_OPTS%LSFORC) THEN
  IF (NCURRENT_ITER > 0 .AND. YDMODEL%YRML_DYN%YRDYNA%LPC_FULL) THEN
    ! this call to CP_FORCING has not been well checked for the corrector step,
    !  but the corrector step cannot work for GFL variables.
    CALL ABOR1('CPG_GP: part 5')
  ENDIF

  IF (YQ%MP /= YQ%MP1) THEN
    CALL ABOR1('CPG_GP: CP_FORCING NEEDS DOUBLE CHECK ON Q POINTER. REK')
  ENDIF
  CALL CP_FORCING(YDGEOMETRY, YDMODEL, YDCPG_OPTS, YDCPG_BNDS, YDCPG_DYN0, YDCPG_TND, YDVARS, YDDDH)
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_GP', 1, ZHOOK_HANDLE)
END SUBROUTINE CPG_GP

