MODULE SPP_MOD_TYPE

 USE MODE_MSG,            ONLY: PRINT_MSG, NVERB_FATAL

 TYPE TSPP_CONFIG_TYPE

  INTEGER :: MP_SELF=-1

  LOGICAL :: LPERT=.FALSE.                
  LOGICAL :: LPRINT=.TRUE.                
  LOGICAL :: LLNN_MEAN1=.FALSE.
  LOGICAL :: LPERT_UNIFORM=.FALSE.

  REAL :: CMPERT
  REAL :: UNIFORM_OFFSET
  REAL :: SDEV
  REAL :: CLIP(2)
  REAL, POINTER :: PGP2DSPP(:) => NULL(), &
                   PTRNDIAG(:) => NULL()

  CHARACTER(LEN=20) :: CTAG = '#'

 END TYPE TSPP_CONFIG_TYPE

 TYPE ALL_SPP_VARS

  ! Gather all parameter holders for convenience

  TYPE(TSPP_CONFIG_TYPE) :: YSPP_RADGR,YSPP_RADSN, &
  YSPP_CLDDPTH,YSPP_CLDDPTHDP, &
  YSPP_RFAC_TWOC,YSPP_RZC_H,YSPP_RZL_INF, &
  YSPP_PSIGQSAT,YSPP_ICE_CLD_WGT, &
  YSPP_RSWINHF,YSPP_RLWINHF, &
  YSPP_ICENU,YSPP_KGN_ACON,YSPP_KGN_SBGR

 END TYPE ALL_SPP_VARS

 CONTAINS

 !
 !-----------------------------------------------------------------------
 !

 SUBROUTINE CLEAR_SPP_TYPE(TSPP)
  IMPLICIT NONE
  TYPE(TSPP_CONFIG_TYPE), INTENT(INOUT) :: TSPP
 END SUBROUTINE CLEAR_SPP_TYPE

 !
 !-----------------------------------------------------------------------
 !

 SUBROUTINE SET_SPP_TYPE(TSPP,CTAG,LLNN_MEAN1_SELF, &
                         LPERT_UNIFORM, &
                         CMPERT,UNIFORM_OFFSET,SDEV,CLIP,MP_SELF, &
                         KLON,KLEV,NEZDIAG, &
                         KSTA,KEND, &
                         PGP2DSPP,PEZDIAG)
  IMPLICIT NONE
  TYPE(TSPP_CONFIG_TYPE),  INTENT(INOUT) :: TSPP
  CHARACTER(LEN=*),        INTENT(IN   ) :: CTAG
  LOGICAL,                 INTENT(IN   ) :: LLNN_MEAN1_SELF,LPERT_UNIFORM
  REAL           ,         INTENT(IN   ) :: CMPERT,UNIFORM_OFFSET,SDEV,CLIP(2)
  INTEGER           ,      INTENT(IN   ) :: MP_SELF,KLON,KLEV,KSTA,KEND, &
                                            NEZDIAG
  REAL           , TARGET, INTENT(IN   ) :: PGP2DSPP(KLON,0)
  REAL           , TARGET, INTENT(INOUT) :: PEZDIAG(KLON,KLEV,NEZDIAG)
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'SPP_MOD_TYPE', 'SET_SPP_TYPE is not implemented in PHYEX')
 END SUBROUTINE SET_SPP_TYPE

 !
 !-----------------------------------------------------------------------
 !

 SUBROUTINE APPLY_SPP(TSPP, &
                      KLON,KSTA,KEND, &
                      PREFVAL,PFIELD)
  IMPLICIT NONE
  TYPE(TSPP_CONFIG_TYPE),    INTENT(INOUT) :: TSPP
  INTEGER           ,        INTENT(IN   ) :: KLON,KSTA,KEND
  REAL           ,           INTENT(IN   ) :: PREFVAL
  REAL           ,           INTENT(INOUT) :: PFIELD(KLON)
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'SPP_MOD_TYPE', 'APPLY_SPP is not implemented in PHYEX')
 END SUBROUTINE APPLY_SPP

 !
 !-----------------------------------------------------------------------
 !

 SUBROUTINE DIA_SPP(TSPP,KSTA,KEND)
  IMPLICIT NONE
  TYPE(TSPP_CONFIG_TYPE),    INTENT(IN) :: TSPP
  INTEGER           ,        INTENT(IN) :: KSTA,KEND
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'SPP_MOD_TYPE', 'DIA_SPP is not implemented in PHYEX')
 END SUBROUTINE DIA_SPP

 !
 !-----------------------------------------------------------------------
 !
 SUBROUTINE SET_ALL_SPP(KLON,KLEV,NGFL_EZDIAG, &
  KIDIA,KFDIA,PGP2DSPP,PEZDIAG,YSPP_ALL)
  IMPLICIT NONE
  INTEGER           ,         INTENT(IN   ) :: KLON,KLEV,NGFL_EZDIAG,KIDIA,KFDIA
  REAL           , TARGET,    INTENT(IN   ) :: PGP2DSPP(KLON,0)
  REAL           ,            INTENT(INOUT) :: PEZDIAG(KLON,KLEV,NGFL_EZDIAG)
  TYPE(ALL_SPP_VARS),         INTENT(INOUT) :: YSPP_ALL
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'SPP_MOD_TYPE', 'SET_ALL_SPP is not implemented in PHYEX')
 END SUBROUTINE SET_ALL_SPP

 !
 !-----------------------------------------------------------------------
 !

 SUBROUTINE CLEAR_ALL_SPP(YSPP_ALL)
  IMPLICIT NONE
  TYPE(ALL_SPP_VARS),         INTENT(INOUT) :: YSPP_ALL
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'SPP_MOD_TYPE', 'CLEAR_ALL_SPP is not implemented in PHYEX')
 END SUBROUTINE CLEAR_ALL_SPP

END MODULE SPP_MOD_TYPE


