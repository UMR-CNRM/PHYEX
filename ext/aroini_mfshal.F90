SUBROUTINE AROINI_MFSHAL(PALP_PERT,PABUO,PBENTR,PBDETR,PCMF,PENTR_MF,PCRAD_MF,PENTR_DRY,&
 &          PDETR_DRY,PDETR_LUP,PKCF_MF,PKRC_MF,PTAUSIGMF,PPRES_UV,PFRAC_UP_MAX,&
 &          PALPHA_MF,PSIGMA_MF,PA1,PB,PC,PBETA1,PR,PLAMBDA,HMF_UPDRAFT,HMF_CLOUD,OMIXUV)
 
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!**** *AROINI_MFSHAL*   - 
!     Purpose.
!     --------
!           Call Méso-NH routine INI_CMFSHALL 
!              (setup of constants for Mass Flux scheme of Pergaud et al)
!
!**   Interface.
!     ----------
!        *CALL* *AROINI_MFSHAL

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        

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
!        S. Malardel

!     Modifications.
!     --------------
!        Original : 07-10-30
!     ------------------------------------------------------------------

USE MODI_INI_CMFSHALL
USE MODD_PARAM_MFSHALL_n, ONLY: LMIXUV, CMF_UPDRAFT, CMF_CLOUD

IMPLICIT NONE

REAL,   INTENT(IN)   :: PALP_PERT
REAL,   INTENT(IN)   :: PABUO
REAL,   INTENT(IN)   :: PBENTR
REAL,   INTENT(IN)   :: PBDETR
REAL,   INTENT(IN)   :: PCMF
REAL,   INTENT(IN)   :: PENTR_MF
REAL,   INTENT(IN)   :: PCRAD_MF
REAL,   INTENT(IN)   :: PENTR_DRY
REAL,   INTENT(IN)   :: PDETR_DRY
REAL,   INTENT(IN)   :: PDETR_LUP
REAL,   INTENT(IN)   :: PKCF_MF
REAL,   INTENT(IN)   :: PKRC_MF
REAL,   INTENT(IN)   :: PTAUSIGMF
REAL,   INTENT(IN)   :: PPRES_UV
REAL,   INTENT(IN)   :: PFRAC_UP_MAX
REAL,   INTENT(IN)   :: PALPHA_MF
REAL,   INTENT(IN)   :: PSIGMA_MF
REAL,   INTENT(IN)   :: PA1
REAL,   INTENT(IN)   :: PB
REAL,   INTENT(IN)   :: PC
REAL,   INTENT(IN)   :: PBETA1   
REAL,   INTENT(IN)   :: PR   
REAL,   INTENT(IN)   :: PLAMBDA
CHARACTER (LEN=4), INTENT(IN) :: HMF_UPDRAFT
CHARACTER (LEN=4), INTENT(IN) :: HMF_CLOUD
LOGICAL, INTENT(IN) :: OMIXUV

!     ------------------------------------------------------------------

!         1. Set implicit default values for MODD_CMFSHALL

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AROINI_MFSHAL',0,ZHOOK_HANDLE)
CALL INI_CMFSHALL(PALP_PERT,PABUO,PBENTR,PBDETR,PCMF,PENTR_MF,PCRAD_MF,PENTR_DRY,&
 &          PDETR_DRY,PDETR_LUP,PKCF_MF,PKRC_MF,PTAUSIGMF,PPRES_UV,PFRAC_UP_MAX,&
 &          PALPHA_MF,PSIGMA_MF,PA1,PB,PC,PBETA1,PR,PLAMBDA)
!
LMIXUV=OMIXUV
CMF_UPDRAFT=HMF_UPDRAFT
CMF_CLOUD=HMF_CLOUD
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AROINI_MFSHAL',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE AROINI_MFSHAL
