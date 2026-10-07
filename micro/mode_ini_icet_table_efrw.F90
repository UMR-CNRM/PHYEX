!     ######spl

MODULE MODE_INI_ICET_TABLE_EFRW
  IMPLICIT NONE
CONTAINS

  SUBROUTINE INI_ICET_TABLE_EFRW(ICE_T_PARAMETERS)

!     ###########################################################
!
!!****  *INI_ICET_TABLE_EFRW * - initialize the t_Efrw table for ICE-T
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the t_Efrw table for
!!    the ICE-T implementation.
!!
!!**  METHOD
!!    ------
!!      Variable collision efficiency for rain collecting cloud water using
!!    method of Beard and Grover, 1974 if a/A less than 0.25; otherwise
!!    uses polynomials to get close match of Pruppacher & Klett Fig 14-9.
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
!!      Module MODD_RAIN_ICE_PARAM
!!        XITDC
!!        XITDR
!!        XT_EFRW
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
!!      Original        2022
!!      B. van Ulft     01/11/23 moved to separate routine
!!      Rolf H. Myhre   2024 Adapt for CY49T1
!!
!-------------------------------------------------------------------------------

!*       0.    DECLARATIONS
!              ------------

    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_CST, ONLY: XPI, XRHOLW
    USE MODD_PRECISION, ONLY: MNHREAL64
    USE MODD_ICET_PARAM, ONLY: ICET_PARAM_t
    USE MODD_ICET_PARAM, ONLY: NBC, NBR

    IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

    TYPE(ICET_PARAM_t), INTENT(INOUT) :: ICE_T_PARAMETERS

!*       0.2   Declarations of local variables :

    REAL(KIND=MNHREAL64) :: ZVTR, ZSTOKES, ZREYNOLDS, ZEF_RW
    REAL(KIND=MNHREAL64) :: ZP, ZYC0, ZF, ZG, ZH, ZZ, ZK0, ZX
    INTEGER :: JI, JJ

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------

    IF (LHOOK) CALL DR_HOOK('INI_ICET_TABLE_EFRW',0,ZHOOK_HANDLE)

!-------------------------------------------------------------------------------

    DO JJ = 1, NBC
      DO JI = 1, NBR
        ZEF_RW = 0.0
        ZP = ICE_T_PARAMETERS%XITDC(JJ)/ICE_T_PARAMETERS%XITDR(JI)
        IF (ICE_T_PARAMETERS%XITDR(JI).LT.50.E-6 .OR. ICE_T_PARAMETERS%XITDC(JJ).LT.3.E-6) THEN
          ICE_T_PARAMETERS%XT_EFRW(JI,JJ) = 0.0
        ELSEIF (ZP.gt.0.25) THEN
          ZX = ICE_T_PARAMETERS%XITDC(JJ)*1.D6
          IF (ICE_T_PARAMETERS%XITDR(JI) .LT. 75.E-6) THEN
            ZEF_RW = 0.026794*ZX - 0.20604
          ELSEIF (ICE_T_PARAMETERS%XITDR(JI) .LT. 125.E-6) THEN
            ZEF_RW = -0.00066842*ZX*ZX + 0.061542*ZX - 0.37089
          ELSEIF (ICE_T_PARAMETERS%XITDR(JI) .LT. 175.E-6) THEN
            ZEF_RW = 4.091e-06*ZX*ZX*ZX*ZX - 0.00030908*ZX*ZX*ZX              &
                   + 0.0066237*ZX*ZX - 0.0013687*ZX - 0.073022
          ELSEIF (ICE_T_PARAMETERS%XITDR(JI) .LT. 250.E-6) THEN
            ZEF_RW = 9.6719e-5*ZX*ZX*ZX - 0.0068901*ZX*ZX + 0.17305*ZX        &
                   - 0.65988
          ELSEIF (ICE_T_PARAMETERS%XITDR(JI) .LT. 350.E-6) THEN
            ZEF_RW = 9.0488e-5*ZX*ZX*ZX - 0.006585*ZX*ZX + 0.16606*ZX         &
                   - 0.56125
          ELSE
            ZEF_RW = 0.00010721*ZX*ZX*ZX - 0.0072962*ZX*ZX + 0.1704*ZX        &
                   - 0.46929
          ENDIF
        ELSE
          ZVTR = -0.1021 &
             & + 4.932E3*ICE_T_PARAMETERS%XITDR(JI) &
             & - 0.9551E6*ICE_T_PARAMETERS%XITDR(JI) &
             &           *ICE_T_PARAMETERS%XITDR(JI) &
             & + 0.07934E9*ICE_T_PARAMETERS%XITDR(JI) &
             &            *ICE_T_PARAMETERS%XITDR(JI) &
             &            *ICE_T_PARAMETERS%XITDR(JI) &
             & - 0.002362E12*ICE_T_PARAMETERS%XITDR(JI) &
             &              *ICE_T_PARAMETERS%XITDR(JI) &
             &              *ICE_T_PARAMETERS%XITDR(JI) &
             &              *ICE_T_PARAMETERS%XITDR(JI)

          ZSTOKES = ICE_T_PARAMETERS%XITDC(JJ)*ICE_T_PARAMETERS%XITDC(JJ)*ZVTR*XRHOLW/(9.*1.718E-5*ICE_T_PARAMETERS%XITDR(JI))
          ZREYNOLDS = 9.*ZSTOKES/(ZP*ZP*XRHOLW)

          ZF = LOG(ZREYNOLDS)
          ZG = -0.1007D0 - 0.358D0*ZF + 0.0261D0*ZF*ZF
          ZK0 = EXP(ZG)
          ZZ = LOG(ZSTOKES/(ZK0+1.D-15))
          ZH = 0.1465D0 + 1.302D0*ZZ - 0.607D0*ZZ*ZZ + 0.293D0*ZZ*ZZ*ZZ
          ZYC0 = 2.0D0/XPI * ATAN(ZH)
          ZEF_RW = (ZYC0+ZP)*(ZYC0+ZP) / ((1.+ZP)*(1.+ZP))

        ENDIF

        ICE_T_PARAMETERS%XT_EFRW(JI,JJ) = MAX(0.0, MIN(REAL(ZEF_RW), 0.95))

      ENDDO
    ENDDO

!-------------------------------------------------------------------------------

    IF (LHOOK) CALL DR_HOOK('INI_ICET_TABLE_EFRW',1,ZHOOK_HANDLE)

!-------------------------------------------------------------------------------

  END SUBROUTINE INI_ICET_TABLE_EFRW
END MODULE MODE_INI_ICET_TABLE_EFRW
