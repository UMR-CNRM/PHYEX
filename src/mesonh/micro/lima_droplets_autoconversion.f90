!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!      #################################
       MODULE MODI_LIMA_DROPLETS_AUTOCONVERSION
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_DROPLETS_AUTOCONVERSION (LDCOMPUTE,                      &
                                            PRHODREF,                       &
                                            PRCT, PCCT, PLBDC, PLBDR,       &
                                            P_RC_AUTO, P_CC_AUTO, P_CR_AUTO )
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! Cloud water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR   ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RC_AUTO
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CC_AUTO
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CR_AUTO
!
END SUBROUTINE LIMA_DROPLETS_AUTOCONVERSION
END INTERFACE
END MODULE MODI_LIMA_DROPLETS_AUTOCONVERSION
!
!     ##########################################################################
      SUBROUTINE LIMA_DROPLETS_AUTOCONVERSION (LDCOMPUTE,                      &
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
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, LKHKO, LKESSLERAC, NMOM_C
USE MODD_PARAM_LIMA_WARM, ONLY : XLAUTR, XAUTO1, XLAUTR_THRESHOLD, &
                                 XITAUTR, XAUTO2, XITAUTR_THRESHOLD, &
                                 XACCR4, XACCR5, XACCR3, XACCR1, XAC, XR0
USE MODD_CST,             ONLY : XPI, XRHOLW
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT    ! Cloud water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR   ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RC_AUTO
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CC_AUTO
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CR_AUTO
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRCT)) :: ZW1, ZW2, ZW3 ! work arrays
!
!-------------------------------------------------------------------------------
!
P_RC_AUTO(:) = 0.0
P_CC_AUTO(:) = 0.0
P_CR_AUTO(:) = 0.0
!
ZW3(:) = 0.0
ZW2(:) = 0.0
ZW1(:) = 0.0
!
IF (NMOM_C.EQ.1 .AND. LKESSLERAC) THEN
   P_RC_AUTO(:) = - 1.E-3 * MAX ( PRCT(:) - 0.5E-3 / PRHODREF(:), 0. )
ELSE IF (LKHKO) THEN
!
!        1. Autoconversion of cloud droplets (Berry-Reinhardt parameterization)
!   	 ----------------------------------------------------------------------
!
   WHERE ( PRCT(:)>XRTMIN(2) .AND. PCCT(:)>XCTMIN(2) .AND. LDCOMPUTE(:) )
!
      ZW1(:)= 1350.0 * PRCT(:)**(2.47) * (PCCT(:)* PRHODREF(:)/1.0E6)**(-1.79) ! ZCCT in cm-3         
!
      P_RC_AUTO(:) = - ZW1(:)
!
      ZW2(:) = ZW1(:) * 3./(4.*XPI*XRHOLW*(XR0)**(3.))
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
!   	 ----------------------------------------------------------------------
!
   WHERE( PRCT(:)>XRTMIN(2) .AND. PCCT(:)>XCTMIN(2) .AND. PLBDC(:)>0. .AND. LDCOMPUTE(:) )
      ZW2(:) = MAX( 0.0, &
                     XLAUTR*PRHODREF(:)*PRCT(:)*(XAUTO1/min(PLBDC(:),1.e9)**4-XLAUTR_THRESHOLD) ) ! L 
!
      ZW3(:) = MAX( 0.0, &
                     XITAUTR*ZW2(:)*PRCT(:)*(XAUTO2/PLBDC(:)-XITAUTR_THRESHOLD) ) ! L/tau
!
      P_RC_AUTO(:) = - ZW3(:)
!
      ZW1(:) = MIN( MIN( 1.2E4, &
                         (XACCR4/PLBDC(:)-XACCR5)/XACCR3 ), &
                         PLBDR(:)/XACCR1 )         ! D**-1 threshold diameter for 
                                                   ! switching the autoconversion regimes
                                                   ! min (80 microns, D_h, D_r)
      ZW3(:) = ZW3(:) * MAX( 0.0,ZW1(:) )**3 / XAC 
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
END SUBROUTINE LIMA_DROPLETS_AUTOCONVERSION
