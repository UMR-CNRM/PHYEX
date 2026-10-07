!     ######spl

MODULE MODE_INI_ICET_TABLE_EFSW
  IMPLICIT NONE
CONTAINS

  SUBROUTINE INI_ICET_TABLE_EFSW(ICE_T_PARAMETERS, RAIN_ICE_DESCRN)

!     ###########################################################
!
!!****  *INI_ICET_TABLE_EFSW * - initialize the t_Efsw table for ICE-T
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the t_Efsw table for
!!    the ICE-T implementation.
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XRHOLW               ! Liquid water density
!!        XPI
!!      Module MODD_RAIN_ICE_DESCR
!!        XAS
!!        XBS
!!        XCS
!!        XDS
!!      Module MODD_RAIN_ICE_PARAM
!!        NBC
!!        NBS
!!        XAM_R
!!        XD0S
!!        XITDC
!!        XITDS
!!        XOBMR
!!        XT_EFSW
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      B.J.K. Engdahl   * MetNo *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    2022
!!      B. van Ulft 01/11/23 moved to separate routine
!!      Rolf H. Myhre   2024 Adapt for CY49T1
!!
!-------------------------------------------------------------------------------

!*       0.    DECLARATIONS
!              ------------

    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_CST, ONLY: XPI, XRHOLW
    USE MODD_PRECISION, ONLY: MNHREAL64
    USE MODD_ICET_PARAM, ONLY: ICET_PARAM_t
    USE MODD_ICET_PARAM, ONLY: NBC, NBS, XFV_S, XD0S
    USE MODD_RAIN_ICE_DESCR_N, ONLY: RAIN_ICE_DESCR_t

    IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

    TYPE(ICET_PARAM_t), INTENT(INOUT)    :: ICE_T_PARAMETERS
    TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: RAIN_ICE_DESCRN

!*       0.2   Declarations of local variables :

    REAL(KIND=MNHREAL64) :: ZDS_M, ZVTS, ZVTC, ZSTOKES, ZREYNOLDS, ZEF_SW
    REAL(KIND=MNHREAL64) :: ZP, ZYC0, ZF, ZG, ZH, ZZ, ZK0
    INTEGER :: JI, JJ

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------

    IF (LHOOK) CALL DR_HOOK('INI_ICET_TABLE_EFSW',0,ZHOOK_HANDLE)

!-------------------------------------------------------------------------------

    DO JJ = 1, NBC
      ZVTC = 1.19D4 * (1.0D4*ICE_T_PARAMETERS%XITDC(JJ)*ICE_T_PARAMETERS%XITDC(JJ)*0.25D0)
      DO JI = 1, NBS
        ZVTS = RAIN_ICE_DESCRN%XCS*ICE_T_PARAMETERS%XITDS(JI)**RAIN_ICE_DESCRN%XDS &
           & *EXP(-XFV_S*ICE_T_PARAMETERS%XITDS(JI)) - ZVTC

        ZDS_M = (RAIN_ICE_DESCRN%XAS &
            &   *ICE_T_PARAMETERS%XITDS(JI)**RAIN_ICE_DESCRN%XBS &
            &   /ICE_T_PARAMETERS%XAM_R)**ICE_T_PARAMETERS%XOBMR

        ZP = ICE_T_PARAMETERS%XITDC(JJ)/ZDS_M

        IF (ZP.gt.0.25 .OR. ICE_T_PARAMETERS%XITDS(JI).LT.XD0S .OR. ICE_T_PARAMETERS%XITDC(JJ).LT.6.E-6 &
           .OR. ZVTS.LT.1.E-3) THEN
          ICE_T_PARAMETERS%XT_EFSW(JI,JJ) = 0.0
        ELSE
          ZSTOKES = ICE_T_PARAMETERS%XITDC(JJ)*ICE_T_PARAMETERS%XITDC(JJ)*ZVTS*XRHOLW/(9.*1.718E-5*ZDS_M)
          ZREYNOLDS = 9.*ZSTOKES/(ZP*ZP*XRHOLW)

          ZF = LOG(ZREYNOLDS)
          ZG = -0.1007D0 - 0.358D0*ZF + 0.0261D0*ZF*ZF
          ZK0 = EXP(ZG)
          ZZ = LOG(ZSTOKES/(ZK0+1.D-15))
          ZH = 0.1465D0 + 1.302D0*ZZ - 0.607D0*ZZ*ZZ + 0.293D0*ZZ*ZZ*ZZ
          ZYC0 = 2.0D0/XPI * ATAN(ZH)
          ZEF_SW = (ZYC0+ZP)*(ZYC0+ZP) / ((1.+ZP)*(1.+ZP))

          ICE_T_PARAMETERS%XT_EFSW(JI,JJ) = MAX(0.0, MIN(REAL(ZEF_SW), 0.95))
        ENDIF
      ENDDO
    ENDDO

!-------------------------------------------------------------------------------

    IF (LHOOK) CALL DR_HOOK('INI_ICET_TABLE_EFSW',1,ZHOOK_HANDLE)

!-------------------------------------------------------------------------------

  END SUBROUTINE INI_ICET_TABLE_EFSW
END MODULE MODE_INI_ICET_TABLE_EFSW
