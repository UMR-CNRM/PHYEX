!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_COMPUTE_FUNCTION_THERMO
IMPLICIT NONE
CONTAINS
SUBROUTINE COMPUTE_FUNCTION_THERMO (D, CST, PALP, PBETA, PGAM, PLTT, PC, PT, PEXN, PCP, &
                                   & PLOCPEXN, PAMOIST, PATHETA, PRT, PPABST, KRR)
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
    INTEGER, INTENT(IN) :: KRR
    REAL, INTENT(IN) :: PALP, PBETA, PGAM, PLTT, PC
    REAL, INTENT(IN), DIMENSION(D%NIJT, D%NKT) :: PT, PEXN, PCP
    !
    REAL, INTENT(OUT), DIMENSION(D%NIJT, D%NKT) :: PLOCPEXN
    REAL, INTENT(OUT), DIMENSION(D%NIJT, D%NKT) :: PAMOIST, PATHETA
    REAL, INTENT(INOUT) :: PRT(D%NIJT, D%NKT, KRR)
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
    IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO', 0, ZHOOK_HANDLE2)

    IKT = D%NKT
    IIJE = D%NIJE
    IIJB = D%NIJB
    ZEPS = CST%XMV / CST%XMD
    !
    !*       1.1 Lv/Cph at  t
    !
    ! present(ZRVSAT,ZDRVSATDT) ! present(PLOCPEXN) ! present(ZDRVSATDT)
!$acc kernels present_cr( PLOCPEXN )
!$mnh_expand_array( JIJ=IIJB:IIJE,JK=1:IKT )
    PLOCPEXN(:, :) = (PLTT + (CST%XCPV - PC)*(PT(:, :) - CST%XTT)) / PCP(:, :)
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
!$acc kernels present_cr( ZRVSAT,ZDRVSATDT )
!$mnh_expand_array( JIJ=IIJB:IIJE,JK=1:IKT )
    !*      1.2 Saturation vapor pressure at t
    !
    ZRVSAT(:, :) = EXP(PALP - PBETA / PT(:, :) - PGAM*LOG(PT(:, :)))
    !
    !*      1.3 saturation  mixing ratio at t
    !
    !YS Added protection (AROME 2024-03-12 crashs)
    ZRVSAT(:, :) = ZRVSAT(:, :)*ZEPS / MAX(1.E-3, PPABST(:, :) - ZRVSAT(:, :))
    !
    !*      1.4 compute the saturation mixing ratio derivative (rvs')
    !
    ZDRVSATDT(:, :) = (PBETA / PT(:, :) - PGAM) / PT(:, :)*ZRVSAT(:, :)*(1. +  &
    & ZRVSAT(:, :) / ZEPS)
    !
    !*      1.5 compute Amoist
    !
    PAMOIST(:, :) = 0.5 / (1.0 + ZDRVSATDT(:, :)*PLOCPEXN(:, :))
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    !*      1.6 compute Atheta
    !
!$acc kernels
!$mnh_expand_array( JIJ=IIJB:IIJE,JK=1:IKT )
    PATHETA(:, :) = PAMOIST(:, :)*PEXN(:, :)*((ZRVSAT(:, :) &
    & - PRT(:, :, 1))*PLOCPEXN(:, :) / (1. + ZDRVSATDT(:, :)*& 
    & PLOCPEXN(:, :))*(ZRVSAT(:, :)*(1. + ZRVSAT(:, :) / ZEPS)*&
    & (-2.*PBETA / PT(:, :) + PGAM) / PT(:, :)**2 +  &
    & ZDRVSATDT(:, :)*(1. + 2.*ZRVSAT(:, :) / ZEPS)*(PBETA / PT(:, :) - PGAM) / &
    & PT(:,:)) - ZDRVSATDT(:, :))
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    !*      1.7 Lv/Cph/Exner at t-1
    !
!$acc kernels
!$mnh_expand_array( JIJ=IIJB:IIJE,JK=1:IKT )
    PLOCPEXN(:, :) = PLOCPEXN(:, :) / PEXN(:, :)
!$mnh_end_expand_array ( JIJ=IIJB:IIJE,JK=1:IKT )
!$acc end kernels
    !
    IF (LHOOK) CALL DR_HOOK('TURB:COMPUTE_FUNCTION_THERMO', 1, ZHOOK_HANDLE2)
  END SUBROUTINE COMPUTE_FUNCTION_THERMO
END MODULE MODE_COMPUTE_FUNCTION_THERMO


