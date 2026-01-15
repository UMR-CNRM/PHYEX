!MNH_LIC Copyright 2013-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_LIMA_DROPLETS_HOM_FREEZING
  IMPLICIT NONE
CONTAINS
!     ##########################################################################
  SUBROUTINE LIMA_DROPLETS_HOM_FREEZING (CST, LIMAP, LIMAC, KSIZE, PTSTEP,  ODCOMPUTE,        &
                                         PT, PLVFACT, PLSFACT,             &
                                         PRCT, PCCT, PLBDC,                &
                                         P_TH_HONC, P_RC_HONC, P_CC_HONC   )
!     ##########################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the cloud droplets homogeneous freezing rate
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T 
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK 
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(CST_T),INTENT(IN)::CST
INTEGER,              INTENT(IN)    :: KSIZE
REAL,                 INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT        ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLVFACT   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRCT      ! Cloud water m.r. at t 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCCT      ! Cloud water C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDC     ! Cloud water lambda
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_HONC
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_HONC
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CC_HONC
!
!*       0.2   Declarations of local variables :
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL, DIMENSION(SIZE(PT)) ::  ZZW, ZZX, ZZY, ZTCELSIUS
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Cloud droplets homogeneous freezing
!               -----------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPLETS_HOM_FREEZING', 0, ZHOOK_HANDLE)
P_TH_HONC(:) = 0.
P_RC_HONC(:) = 0.
P_CC_HONC(:) = 0.
!
WHERE ( PT(:)<CST%XTT-35.0 .AND. PCCT(:)>LIMAP%XCTMIN(2) .AND. PRCT(:)>LIMAP%XRTMIN(2) )
   ZTCELSIUS(:) = PT(:)-CST%XTT                                    ! T [°C]
   !
   ZZW(:) = 0.0
   ZZX(:) = 0.0
   ZZY(:) = 0.0

   ZZX(:) = 1.0 / ( 1.0 + (LIMAC%XC_HONC/PLBDC(:))*PTSTEP*           &
           EXP( LIMAC%XTEXP1_HONC + ZTCELSIUS(:)*(                   &
           LIMAC%XTEXP2_HONC + ZTCELSIUS(:)*(                        &
           LIMAC%XTEXP3_HONC + ZTCELSIUS(:)*(                        &
           LIMAC%XTEXP4_HONC + ZTCELSIUS(:)*LIMAC%XTEXP5_HONC))) ) )**LIMAP%XNUC
!
   ZZW(:) = PCCT(:) * (1.0 - ZZX(:))                    ! CCHONI
   ZZY(:) = PRCT(:) * (1.0 - ZZX(:))                    ! RCHONI
!
   P_RC_HONC(:) = - ZZY(:)/PTSTEP
   P_CC_HONC(:) = - ZZW(:)/PTSTEP
!   P_TH_HONC(:) = P_RC_HONC(:) * (PLSFACT(:)-PLVFACT(:))
!
END WHERE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPLETS_HOM_FREEZING', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_DROPLETS_HOM_FREEZING
END MODULE MODE_LIMA_DROPLETS_HOM_FREEZING
