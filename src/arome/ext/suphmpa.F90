SUBROUTINE SUPHMPA(YDGEOMETRY,YDLDDH,YDML_GCONF,YDDYNA,YDML_PHY_MF,KULOUT)

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
!        Y.Seity: 7 Jul  2023 : remove test on NGFL_EZDIAG (as these GFL fields
!        are no longer used, replaced by IMF LMF and AMF)
!        R. Honnert (sept 2023) : GOGER and LEONARD key in INI_PHYEX 
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULNAM

USE YOMDYNA  , ONLY  : TDYNA
USE YOMLDDH   , ONLY : TLDDH
USE YOMCT0 ,ONLY : LELAM

USE MODD_BUDGET, ONLY : TBUCONF_ASSOCIATE, TBUCONF
USE MODI_INI_PHYEX, ONLY: INI_PHYEX
USE MODD_IO, ONLY : TFILEDATA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TLDDH)       ,INTENT(INOUT) :: YDLDDH
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TDYNA), INTENT(IN)          :: YDDYNA
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT
INTEGER, SAVE                    :: PHYEX_NUM = 0
TYPE(TFILEDATA) :: TPFILE

!     ------------------------------------------------------------------


REAL(KIND=JPRB) :: ZTSTEP, ZDZMIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
CHARACTER(LEN=4) :: CSCONV

LOGICAL :: LLNOTMAP
!     ------------------------------------------------------------------

#include "sucvmnh.intfb.h"
#include "aroini_conf.h"
#include "aroini_budget.h"
#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('SUPHMPA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
  & YDPHY=>YDML_PHY_MF%YRPHY,YDRIP=>YDML_GCONF%YRRIP,YDARPHY=>YDML_PHY_MF%YRARPHY,YDPARAR=>YDML_PHY_MF%YRPARAR)
ASSOCIATE(CMICRO=>YDPARAR%CMICRO, CTURB=>YDPARAR%CTURB, NSV=>YDPARAR%NSV, &
 & LAROBU_ENABLE=>YDPARAR%LAROBU_ENABLE, &
 & NRR=>YDPARAR%NRR,&
 & LMPA=>YDARPHY%LMPA, LKFBCONV=>YDARPHY%LKFBCONV, LMFSHAL=>YDARPHY%LMFSHAL, &
 & LLEONARD=>YDARPHY%LLEONARD,LGOGER=>YDARPHY%LGOGER,&
 & NPROMA=>YDDIM%NPROMA, &
 & LEDKF=>YDPHY%LEDKF, LCVPPKF=>YDPHY%LCVPPKF, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LSDDH=>YDLDDH%LSDDH, TSTEP=>YDRIP%TSTEP, &
 & RADGR=>YDPARAR%RADGR, RADSN=>YDPARAR%RADSN,&
 & PHYEX=>YDPARAR%PHYEX)
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!       1. Initialisation of MesoNH constantes
!          Initialisation for microphysics scheme
!          Initialisation of Mass Flux Shallow convection scheme
!
ZTSTEP=TSTEP
ZDZMIN=20.
IF(LMFSHAL.OR.LEDKF) THEN
  CSCONV='EDKF'
ELSE
  CSCONV='NONE'
ENDIF

!          Initialisation of Turbulence scheme
!          Test of constistancy between NAMPARAR and NAMARPHY
!       

!INCREMENTING KTO
PHYEX_NUM=PHYEX_NUM+1
TPFILE%NLU = NULNAM

CALL INI_PHYEX(HPROGRAM=PHYEX%MISC%CPROGRAM, TPFILE=TPFILE, &
             & LDNEEDNAM=.TRUE., KLUOUT=KULOUT, KFROM=0, KTO=PHYEX_NUM, &
             & PTSTEP=ZTSTEP, PDZMIN=ZDZMIN, &
             & CMICRO=CMICRO, CSCONV=CSCONV, CTURB=CTURB, &
             & KPRINT=2, &
             & PHYEX_OUT=PHYEX)

! Ensure consistency
IF (.NOT. PHYEX%TURBN%LGOGER .EQV. LGOGER .OR. .NOT. PHYEX%TURBN%LLEONARD .EQV. LLEONARD) THEN
  CALL ABOR1 ("LLEONARD or LGOGER are not consistent in NAM_TURBn and NAMARPHY.")
ENDIF
IF (.NOT. PHYEX%PARAM_ICEN%LOCND2) THEN
   RADGR=0._JPRB
   RADSN=0._JPRB
ENDIF
IF (LELAM) THEN
  LLNOTMAP=.NOT.YDGEOMETRY%YREGEO%LMAP
ELSE
  LLNOTMAP=.TRUE.
ENDIF
CALL AROINI_CONF (KULOUT,YDDYNA%LTWOTL,LLNOTMAP) !only needed for LFLAT key used in turbulence
                                                 !call must be suppressed once turbulence source code will be updated

!       3. Initialisation of Budget

LAROBU_ENABLE=LMPA.AND.LSDDH
CALL AROINI_BUDGET(LAROBU_ENABLE)

!       4. PHYEX%MISC
PHYEX%MISC%TBUCONF = TBUCONF
IF (LELAM) THEN
  PHYEX%MISC%OFLAT=.NOT.YDGEOMETRY%YREGEO%LMAP
ELSE
  PHYEX%MISC%OFLAT=.TRUE.
ENDIF

IF (PHYEX%TURBN%LHARAT .AND. PHYEX%PARAM_MFSHALLN%CMF_UPDRAFT == 'EDKF') THEN
  CALL ABOR1('Combination LHARATU and EDKF not valid!')
ENDIF

IF (PHYEX%TURBN%CTURBDIM == '3DIM') THEN
  CALL ABOR1('TURBDIM cannot be 3DIM with AROME')
ENDIF

!       5. Initialisation of Convection scheme

IF(LKFBCONV.OR.LCVPPKF) THEN
  CALL SUCVMNH(YDML_PHY_MF,KULOUT)
ENDIF


! -----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPHMPA',1,ZHOOK_HANDLE)
END SUBROUTINE SUPHMPA
