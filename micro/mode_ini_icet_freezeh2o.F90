!     ######spl

MODULE MODE_INI_ICET_FREEZEH2O
  IMPLICIT NONE
CONTAINS

  SUBROUTINE INI_ICET_FREEZEH2O(ICE_T_PARAMETERS, RAIN_ICE_DESCRN)

!     ###########################################################
!
!!****  *INI_ICET_FREEZEH2O * - initialize the tpi_qcfz table for ICE-T
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the tpi_qcfz table for
!!    the ICE-T implementation.
!!
!!**  METHOD
!!    ------
!!      This is a literal adaptation of Bigg (1954) probability of drops of
!!    a particular volume freezing.  Given this probability, simply freeze
!!    the proportion of drops summing their masses.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XRHOLW               ! Liquid water density
!!      Module MODD_RAIN_ICE_DESCR
!!        XAR
!!        XBR
!!      Module MODD_RAIN_ICE_PARAM
!!        NBC
!!        NBR
!!        NTB_IN
!!        NTB_R
!!        NTB_R1
!!        XITDR
!!        XITDC
!!        XMU_R
!!        XN0R_EXP
!!        XOBMR
!!        XTPI_QCFZ
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
    USE MODD_CST, ONLY: XRHOLW
    USE MODD_PRECISION, ONLY: MNHREAL64
    USE MODD_RAIN_ICE_DESCR_N, ONLY: RAIN_ICE_DESCR_t
    USE MODD_ICET_PARAM, ONLY: ICET_PARAM_t
    USE MODD_ICET_PARAM, ONLY: XNT_IN, NTb_R1, NTB_R, XN0R_EXP, XR_R, XMU_R,NTB_C, XR_C
    USE MODD_ICET_PARAM, ONLY: NTB_IN, NTB_C, NBR, NBC

    IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

    TYPE(ICET_PARAM_t), INTENT(INOUT)    :: ICE_T_PARAMETERS
    TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: RAIN_ICE_DESCRN

!*       0.2   Declarations of local variables :

    INTEGER :: JI, JJ, JK, JM, JN, JN2
    INTEGER :: INU_C
    REAL(KIND=MNHREAL64), DIMENSION(NBR) :: ZN_R, ZMASSR
    REAL(KIND=MNHREAL64), DIMENSION(NBC) :: ZN_C, ZMASSC
    REAL(KIND=MNHREAL64) :: ZSUM1, ZSUM2, ZSUMN1, ZSUMN2, &
                        ZPROB, ZVOL, ZTEXP, ZORHO_W, &
                        ZLAM_EXP, ZLAM_R, ZN0_R, ZLAM_C, ZN0_C
    REAL :: ZT_ADJUST

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------

    IF (LHOOK) CALL DR_HOOK('INI_ICET_FREEZEH2O',0,ZHOOK_HANDLE)

!-------------------------------------------------------------------------------

    ZORHO_W = 1./XRHOLW
    DO JN2 = 1, NBR
      ZMASSR(JN2) = RAIN_ICE_DESCRN%XAR*ICE_T_PARAMETERS%XITDR(JN2)**RAIN_ICE_DESCRN%XBR
    ENDDO
    DO JN = 1, NBC
      ZMASSC(JN) = RAIN_ICE_DESCRN%XAR*ICE_T_PARAMETERS%XITDC(JN)**RAIN_ICE_DESCRN%XBR
    ENDDO

    !..Freeze water (smallest drops become cloud ice, otherwise graupel).
    DO JM = 1, NTB_IN
      ZT_ADJUST = MAX(-3.0, MIN(3.0 - LOG10(XNT_IN(JM)), 3.0))
      DO JK = 1, 45
        ZTEXP = EXP( DFLOAT(JK) - ZT_ADJUST*1.0D0 ) - 1.0D0
        DO JJ = 1, NTB_R1
          DO JI = 1, NTB_R
            ZLAM_EXP = (XN0R_EXP(JJ)*RAIN_ICE_DESCRN%XAR*ICE_T_PARAMETERS%XCR_GM(1)/XR_R(JI))**ICE_T_PARAMETERS%XORE1
            ZLAM_R = ZLAM_EXP &
                 & * (ICE_T_PARAMETERS%XCR_GM(3)*ICE_T_PARAMETERS%XORG2*ICE_T_PARAMETERS%XORG1)**ICE_T_PARAMETERS%XOBMR
            ZN0_R = XN0R_EXP(JJ)/(ICE_T_PARAMETERS%XCR_GM(2)*ZLAM_EXP) * ZLAM_R**ICE_T_PARAMETERS%XCR_EX(2)
            ZSUM1 = 0.0D0
            ZSUM2 = 0.0D0
            ZSUMN1 = 0.0D0
            ZSUMN2 = 0.0D0
            DO JN2 = NBR, 1, -1
              ZN_R(JN2) = ZN0_R*ICE_T_PARAMETERS%XITDR(JN2)**XMU_R*EXP(-ZLAM_R*ICE_T_PARAMETERS%XITDR(JN2))&
                      & * ICE_T_PARAMETERS%XITDTR(JN2)
              ZVOL = ZMASSR(JN2)*ZORHO_W
              ZPROB = 1.0D0 - EXP(-120.0D0*ZVOL*5.2D-4 * ZTEXP)
              IF (ZMASSR(JN2) .LT. ICE_T_PARAMETERS%XM0G) THEN
                ZSUMN1 = ZSUMN1 + ZPROB*ZN_R(JN2)
                ZSUM1 = ZSUM1 + ZPROB*ZN_R(JN2)*ZMASSR(JN2)
              ELSE
                ZSUMN2 = ZSUMN2 + ZPROB*ZN_R(JN2)
                ZSUM2 = ZSUM2 + ZPROB*ZN_R(JN2)*ZMASSR(JN2)
              ENDIF
              IF ((ZSUM1+ZSUM2) .GE. XR_R(JI)) EXIT
            ENDDO
            ICE_T_PARAMETERS%XTPI_QRFZ(JI,JJ,JK,JM) = ZSUM1
            ICE_T_PARAMETERS%XTPG_QRFZ(JI,JJ,JK,JM) = ZSUM2
          ENDDO
        ENDDO

        DO JJ = 1, NBC
          INU_C = MIN(15, NINT(1000.E6/ICE_T_PARAMETERS%XT_NC(JJ)) + 2)
          DO JI = 1, NTB_C
            ZLAM_C = (ICE_T_PARAMETERS%XT_NC(JJ) * RAIN_ICE_DESCRN%XAR * ICE_T_PARAMETERS%XCC_GM(2,INU_C)&
                 & * ICE_T_PARAMETERS%XOCG1(INU_C) / XR_C(JI))**ICE_T_PARAMETERS%XOBMR
            ZN0_C = ICE_T_PARAMETERS%XT_NC(JJ) * ICE_T_PARAMETERS%XOCG1(INU_C) * ZLAM_C**ICE_T_PARAMETERS%XCC_EX(1,INU_C)
            ZSUM1 = 0.0d0
            ZSUMN2 = 0.0d0
            DO JN = NBC, 1, -1
              ZVOL = ZMASSC(JN)*ZORHO_W
              ZPROB = 1.0D0 - EXP(-120.0D0*ZVOL*5.2D-4 * ZTEXP)
              ZN_C(JN) = ZN0_C*ICE_T_PARAMETERS%XITDC(JN)**INU_C&
                     & * EXP(-ZLAM_C*ICE_T_PARAMETERS%XITDC(JN))*ICE_T_PARAMETERS%XITDTC(JN)
              ZSUM1 = ZSUM1 + ZPROB*ZN_C(JN)*ZMASSC(JN)
              IF (ZSUM1 .GE. XR_C(JI)) EXIT
            ENDDO
            ICE_T_PARAMETERS%XTPI_QCFZ(JI,JJ,JK,JM) = ZSUM1
          ENDDO
        ENDDO
      ENDDO
    ENDDO

!-------------------------------------------------------------------------------

    IF (LHOOK) CALL DR_HOOK('INI_ICET_FREEZEH2O',1,ZHOOK_HANDLE)

!-------------------------------------------------------------------------------

  END SUBROUTINE INI_ICET_FREEZEH2O
END MODULE MODE_INI_ICET_FREEZEH2O
