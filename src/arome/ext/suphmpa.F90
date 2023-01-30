SUBROUTINE SUPHMPA(YDGEOMETRY,YDLDDH,YDML_GCONF,YDML_PHY_MF,KULOUT)

!**** *SUPHMPA*   - Initialize common meso_NH MODD_ used in physics for AROME

!     Purpose.
!     --------
!           Initialize MODD_PARAMETERS, MODD_CST, MODD_CONF,
!           MODD_RAIN_ICE_DESCR, MODD_RAIN_ICE_PARAM, MODD_BUDGET  
!           parameters used in meso_NH Physics and aladin/meso_NH physics
!           interface 

!**   Interface.
!     ----------
!        *CALL* *SUPHMPA(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY2

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation AROME 

!     Author.
!     -------
!       R. Zaaboul
!        Original : 28-Feb-2006

!     Modifications.
!     --------------
!           E. BAZILE : 01-09-2006 : Modified for LCVPPKF.
!        Y. Seity : 28-March-2007 Move chemistry setup under suphmse       
!        O.Riviere: 01/10/2008 removal of call to now obsolete aro_iniapft
!        S.Riette: 24 Aug 2011 add call to AROINI_NEB
!        Y.Seity: 9 Feb 2014 : add autoconversion setup (*CRIAUT*)
!        Y.Seity: 12 Nov  2014 : add test on NGFL_EZDIAG
!        S. Riette (Jan 2015): new ICE3 and ICE4 parameters with new aroini_micro interface
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOMLUN    ,ONLY : NULNAM

USE YOMLDDH   , ONLY : TLDDH
USE YOMCT0 ,ONLY : LTWOTL, LELAM

USE MODD_BUDGET, ONLY : TBUCONF_ASSOCIATE
USE MODI_INI_PHYEX, ONLY: INI_PHYEX

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TLDDH)       ,INTENT(INOUT) :: YDLDDH
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT

!     ------------------------------------------------------------------


REAL(KIND=JPRB) :: ZTSTEP, ZDZMIN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
CHARACTER(LEN=6) :: CPROGRAM

LOGICAL :: LLNOTMAP
!     ------------------------------------------------------------------

#include "sucvmnh.intfb.h"
#include "aroini_cstmnh.h"
#include "aroini_budget.h"
#include "aroini_turb.h"
#include "abor1.intfb.h"
#include "aroini_mfshal.h"

#include "aroini_micro_lima.h"

IF (LHOOK) CALL DR_HOOK('SUPHMPA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
  & YDPHY=>YDML_PHY_MF%YRPHY,YDRIP=>YDML_GCONF%YRRIP,YDARPHY=>YDML_PHY_MF%YRARPHY,YDPARAR=>YDML_PHY_MF%YRPARAR)
ASSOCIATE(XDETR_LUP=>YDPARAR%XDETR_LUP, XCMF=>YDPARAR%XCMF, &
 & XBDETR=>YDPARAR%XBDETR, XLINI=>YDPARAR%XLINI, XABUO=>YDPARAR%XABUO, &
 & XLAMBDA=>YDPARAR%XLAMBDA, &
 & XKCF_MF=>YDPARAR%XKCF_MF, XALP_PERT=>YDPARAR%XALP_PERT, &
 & NSPLITR=>YDPARAR%NSPLITR, NSPLITG=>YDPARAR%NSPLITG, XSIGMA_MF=>YDPARAR%XSIGMA_MF, XA1=>YDPARAR%XA1, &
 & CMICRO=>YDPARAR%CMICRO, XENTR_DRY=>YDPARAR%XENTR_DRY, &
 & XENTR_MF=>YDPARAR%XENTR_MF, NSV=>YDPARAR%NSV, &
 & XFRAC_UP_MAX=>YDPARAR%XFRAC_UP_MAX, XB=>YDPARAR%XB, XC=>YDPARAR%XC, &
 & XTAUSIGMF=>YDPARAR%XTAUSIGMF, &
 & XDETR_DRY=>YDPARAR%XDETR_DRY, XR=>YDPARAR%XR, &
 & XBENTR=>YDPARAR%XBENTR, XBETA1=>YDPARAR%XBETA1, &
 & LAROBU_ENABLE=>YDPARAR%LAROBU_ENABLE, &
 & XKRC_MF=>YDPARAR%XKRC_MF, XALPHA_MF=>YDPARAR%XALPHA_MF, &
 & XPRES_UV=>YDPARAR%XPRES_UV, NRR=>YDPARAR%NRR, XCRAD_MF=>YDPARAR%XCRAD_MF, &
 & CMF_UPDRAFT=>YDPARAR%CMF_UPDRAFT, LHARATU=>YDPARAR%LHARATU, LOSUBG_COND=>YDPARAR%LOSUBG_COND,&
 & LSTATNW=>YDPARAR%LSTATNW, &
 & LMPA=>YDARPHY%LMPA, LKFBCONV=>YDARPHY%LKFBCONV, LMFSHAL=>YDARPHY%LMFSHAL, &
 & LGRADHPHY=>YDARPHY%LGRADHPHY, &
 & NPROMA=>YDDIM%NPROMA, &
 & LEDKF=>YDPHY%LEDKF, LCVPPKF=>YDPHY%LCVPPKF, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LSDDH=>YDLDDH%LSDDH, TSTEP=>YDRIP%TSTEP, &
 & LMIXUV=>YDPARAR%LMIXUV,CMF_CLOUD=>YDPARAR%CMF_CLOUD,CCONDENS=>YDPARAR%CCONDENS,&
 & CLAMBDA3=>YDPARAR%CLAMBDA3,CSUBG_MF_PDF=>YDPARAR%CSUBG_MF_PDF,LSIGMAS=>YDPARAR%LOSIGMAS,&
 & RADGR=>YDPARAR%RADGR, RADSN=>YDPARAR%RADSN,&
 & PARAM_ICE=>YDPARAR%PARAM_ICE, RAIN_ICE_DESCR=>YDPARAR%RAIN_ICE_DESCR, RAIN_ICE_PARAM=>YDPARAR%RAIN_ICE_PARAM, &
 & CLOUDPARN=>YDPARAR%CLOUDPARN)
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!       1. Initialisation of MesoNH constantes

IF (LELAM) THEN
  LLNOTMAP=.NOT.YDGEOMETRY%YREGEO%LMAP
ELSE
  LLNOTMAP=.TRUE.
ENDIF
CALL AROINI_CSTMNH (KULOUT,LTWOTL,LLNOTMAP)

!       2. Initialisation for microphysics scheme
CPROGRAM='AROME'
ZTSTEP=TSTEP
ZDZMIN=20.
CALL INI_PHYEX(CPROGRAM, NULNAM, .TRUE., KULOUT, 0, 1, &
               ZTSTEP, ZDZMIN, &
               CMICRO, &
               PARAM_ICE_INOUT=PARAM_ICE, RAIN_ICE_DESCR_INOUT=RAIN_ICE_DESCR, &
               RAIN_ICE_PARAM_INOUT=RAIN_ICE_PARAM, CLOUDPARN_INOUT=CLOUDPARN)

! Ensure consistency
IF (.NOT. PARAM_ICE%LOCND2) THEN
   RADGR=0._JPRB
   RADSN=0._JPRB
ENDIF

IF (CMICRO == 'LIMA') THEN
  CALL AROINI_MICRO_LIMA (KULOUT,NULNAM,ZTSTEP,CMICRO,NSPLITR,NSPLITG)
ENDIF

!       3. Initialisation of Budget

LAROBU_ENABLE=LMPA.AND.LSDDH
CALL AROINI_BUDGET(LAROBU_ENABLE)


!       4. Initialisation of Turbulence scheme

CALL AROINI_TURB(XLINI,LHARATU,LSTATNW,LOSUBG_COND,CCONDENS,CLAMBDA3,CSUBG_MF_PDF,LSIGMAS)

!       5. Initialisation of Mass Flux Shallow convection scheme

IF(LMFSHAL.OR.LEDKF) CALL AROINI_MFSHAL(XALP_PERT,XABUO,XBENTR,XBDETR,XCMF,XENTR_MF,XCRAD_MF,XENTR_DRY,&
 &          XDETR_DRY,XDETR_LUP,XKCF_MF,XKRC_MF,XTAUSIGMF,XPRES_UV,XFRAC_UP_MAX,&
 &          XALPHA_MF,XSIGMA_MF,XA1,XB,XC,XBETA1,XR,XLAMBDA,CMF_UPDRAFT,CMF_CLOUD,LMIXUV)

IF (LMFSHAL.AND.YDML_GCONF%YGFL%NGFL_EZDIAG < 3) THEN
  CALL ABOR1 ("With LMFSHAL NGFL_EZDIAG should be >= 3 !")
ENDIF

!       6. Initialisation of Convection scheme

IF(LKFBCONV.OR.LCVPPKF) THEN
  CALL SUCVMNH(YDML_PHY_MF,KULOUT)
ENDIF

!       7. Initialisation of nebulosity computation

CALL AROINI_NEB

!       8. Initialisation of The Horizontal Gradient on Z levels for 3D turbulence 
!       Quand il y aura des initialisations 
IF (LGRADHPHY .AND. .NOT. LELAM) THEN
  CALL ABOR1 ("With LGRADHPHY, LELAM should be TRUE !")
ENDIF

! -----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPHMPA',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHMPA
