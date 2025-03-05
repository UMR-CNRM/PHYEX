!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_DROPLETS_AUTOCONVERSION
  IMPLICIT NONE
CONTAINS
!     ##########################################################################
  SUBROUTINE LIMA_DROPLETS_AUTOCONVERSION (CST, LIMAP, LIMAW, KSIZE, ODCOMPUTE,               &
                                           PRHODREF,                       &
                                           PRCT, PCCT, PLBDC, PLBDR,       &
                                           P_RC_AUTO, P_CC_AUTO, P_CR_AUTO )
!     ##########################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the autoconversion of cloud droplets
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018 
!!      B. Vie 02/03/2020 : missing CC process
!       Delbeke/Vie     03/2022 : KHKO option
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(CST_T),INTENT(IN)::CST
INTEGER,              INTENT(IN)    :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCCT    ! Cloud water conc. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDR   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_AUTO
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CC_AUTO
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CR_AUTO
!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL, DIMENSION(SIZE(PRCT)) :: ZW1, ZW2, ZW3 ! work arrays
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPLETS_AUTOCONVERSION', 0, ZHOOK_HANDLE)
P_RC_AUTO(:) = 0.0
P_CC_AUTO(:) = 0.0
P_CR_AUTO(:) = 0.0
!
ZW3(:) = 0.0
ZW2(:) = 0.0
ZW1(:) = 0.0
!
IF (LIMAP%NMOM_C.EQ.1 .AND. LIMAP%LKESSLERAC) THEN
   P_RC_AUTO(:) = - 1.E-3 * MAX ( PRCT(:) - 0.5E-3 / PRHODREF(:), 0. )
ELSE IF (LIMAP%LKHKO) THEN
!
!        1. Autoconversion of cloud droplets
!        -----------------------------------
!
   WHERE ( PRCT(:)>LIMAP%XRTMIN(2) .AND. PCCT(:)>LIMAP%XCTMIN(2) .AND. ODCOMPUTE(:) )
!
      ZW1(:)= 1350.0 * PRCT(:)**(2.47) * (PCCT(:)* PRHODREF(:)/1.0E6)**(-1.79) ! ZCCT in cm-3         
!
      P_RC_AUTO(:) = - ZW1(:)
!
      ZW2(:) = ZW1(:) * 3./(4.*CST%XPI*CST%XRHOLW*(LIMAW%XR0)**(3.))
      P_CR_AUTO(:) = ZW2(:)
!
      ZW3(:) = - ZW1(:) * PCCT(:) / PRCT(:)
      P_CC_AUTO(:) = ZW3(:)
!
   END WHERE
!
ELSE
!
!        2. Autoconversion of cloud droplets (Berry-Reinhardt parameterization)
!        ----------------------------------------------------------------------
!
   WHERE( PRCT(:)>LIMAP%XRTMIN(2) .AND. PCCT(:)>LIMAP%XCTMIN(2) .AND. PLBDC(:)>0. .AND. ODCOMPUTE(:) )
      ZW2(:) = MAX( 0.0, &
                     LIMAW%XLAUTR*PRHODREF(:)*PRCT(:)*(LIMAW%XAUTO1/MIN(PLBDC(:),1.E9)**4-LIMAW%XLAUTR_THRESHOLD) ) ! L 
!
      ZW3(:) = MAX( 0.0, &
                     LIMAW%XITAUTR*ZW2(:)*PRCT(:)*(LIMAW%XAUTO2/PLBDC(:)-LIMAW%XITAUTR_THRESHOLD) ) ! L/tau
!
      P_RC_AUTO(:) = - ZW3(:)
!
      ZW1(:) = MIN( MIN( 1.2E4, &
                         (LIMAW%XACCR4/PLBDC(:)-LIMAW%XACCR5)/LIMAW%XACCR3 ), &
                         PLBDR(:)/LIMAW%XACCR1 )         ! D**-1 threshold diameter for 
                                                   ! switching the autoconversion regimes
                                                   ! min (80 microns, D_h, D_r)
      ZW3(:) = ZW3(:) * MAX( 0.0,ZW1(:) )**3 / LIMAW%XAC 
!
      P_CC_AUTO(:) = -ZW3(:)
      P_CR_AUTO(:) = ZW3(:)
!
   END WHERE
!
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPLETS_AUTOCONVERSION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_DROPLETS_AUTOCONVERSION
END MODULE MODE_LIMA_DROPLETS_AUTOCONVERSION
