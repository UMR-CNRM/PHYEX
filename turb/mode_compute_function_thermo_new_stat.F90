!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_COMPUTE_FUNCTION_THERMO_NEW_STAT
IMPLICIT NONE
CONTAINS
SUBROUTINE COMPUTE_FUNCTION_THERMO_NEW_STAT (D, CST, PALP, PBETA, PGAM, PLTT, PC, PT, PEXN, PCP, &
                                             & PLOCPEXN, PAMOIST, PATHETA, PPABST)
    !     ########################################################################
    !!
    !!****  *COMPUTE_FUNCTION_THERMO* routine to compute several thermo functions
    !
    !!    AUTHOR
    !!    ------
    !!
    !!     JP Pinty      *LA*
    !!
    !!    MODIFICATIONS
    !!    -------------
    !!      Original   24/02/03
    !!     Modified: Wim de Rooy 06-02-2019
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !              ------------
    !
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MODD_CST, ONLY: CST_t
    USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t

    IMPLICIT NONE
    !
    !*       0.1   Declarations of dummy arguments
    !
    TYPE(DIMPHYEX_t), INTENT(IN) :: D  ! PHYEX variables dimensions structure
    REAL, INTENT(IN) :: PALP, PBETA, PGAM, PLTT, PC
    REAL, INTENT(IN), DIMENSION(D%NIJT, D%NKT) :: PT, PEXN, PCP
    !
    REAL, INTENT(OUT), DIMENSION(D%NIJT, D%NKT) :: PLOCPEXN
    REAL, INTENT(OUT), DIMENSION(D%NIJT, D%NKT) :: PAMOIST, PATHETA
    TYPE(CST_t), INTENT(IN) :: CST
    REAL, INTENT(IN) :: PPABST(D%NIJT, D%NKT)

    REAL :: ZDRVSATDT(D%NIJT, D%NKT)
    REAL :: ZRVSAT(D%NIJT, D%NKT)
    REAL :: ZEPS
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE2
    INTEGER :: JK, JIJ, IIJB, IIJE, IKT
    !
    !-------------------------------------------------------------------------------
    !
    IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO_NEW_STAT', 0, ZHOOK_HANDLE2)

    IKT = D%NKT
    IIJE = D%NIJE
    IIJB = D%NIJB
    ZEPS = CST%XMV / CST%XMD
    !
    !*       1.1 Lv/Cph at  t
    !

DO JK=1, IKT
      DO JIJ=IIJB, IIJE
        PLOCPEXN(JIJ, JK) = (PLTT + (CST%XCPV - PC)*(PT(JIJ, JK) - CST%XTT)) / PCP(JIJ, JK)
        !
        !*      1.2 Saturation vapor pressure at t
        !
        ZRVSAT(JIJ, JK) = EXP(PALP - PBETA / PT(JIJ, JK) - PGAM*LOG(PT(JIJ, JK)))
        !
        !*      1.3 saturation  mixing ratio at t
        !
        ZRVSAT(JIJ, JK) = ZRVSAT(JIJ, JK)*ZEPS / (PPABST(JIJ, JK) - ZRVSAT(JIJ, JK))
        !
        !*      1.4 compute the saturation mixing ratio derivative (rvs')
        !
        ZDRVSATDT(JIJ, JK) = (PBETA / PT(JIJ, JK) - PGAM) / PT(JIJ, JK)*ZRVSAT(JIJ, JK)*(1. +  &
        & ZRVSAT(JIJ, JK) / ZEPS)
        !
        !*      1.5 compute Amoist
        !
        PAMOIST(JIJ, JK) = 1.0 / (1.0 + ZDRVSATDT(JIJ, JK)*PLOCPEXN(JIJ, JK))
        !
        !*      1.6 compute Atheta
        !
        PATHETA(JIJ, JK) = PAMOIST(JIJ, JK)*PEXN(JIJ, JK)*ZDRVSATDT(JIJ, JK)
        !
        !*      1.7 Lv/Cph/Exner at t-1
        !
        PLOCPEXN(JIJ, JK) = PLOCPEXN(JIJ, JK) / PEXN(JIJ, JK)
  END DO
    END DO

    !
    IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO_NEW_STAT', 1, ZHOOK_HANDLE2)
  END SUBROUTINE COMPUTE_FUNCTION_THERMO_NEW_STAT
END MODULE MODE_COMPUTE_FUNCTION_THERMO_NEW_STAT
