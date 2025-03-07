!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_RAIN_EVAPORATION
  IMPLICIT NONE
CONTAINS
!     ###############################################################################
  SUBROUTINE LIMA_RAIN_EVAPORATION (CST, LIMAP, LIMAW, KSIZE, PTSTEP, ODCOMPUTE,                   &
                                    PRHODREF, PT, PLV, PLVFACT, PEVSAT, PRVSAT, &
                                    PRVT, PRCT, PRRT, PCRT, PLBDR,              &
                                    P_TH_EVAP, P_RR_EVAP, P_CR_EVAP,            &
                                    PEVAP3D                                     )
!     ###############################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the raindrop evaporation
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!!
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
REAL,                 INTENT(IN)    :: PTSTEP     ! Time step
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE  !
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT         ! Temperature
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLV        ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLVFACT    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PEVSAT     !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRVSAT     !
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCRT       ! Rain water conc at t 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDR      ! Lambda(rain)
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_EVAP
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RR_EVAP
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CR_EVAP
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: PEVAP3D    ! Rain evap profile
!
!*       0.1   Declarations of local variables :
!
! 
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GEVAP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL, DIMENSION(SIZE(PRHODREF))    :: ZZW1, ZZW2
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PREPARE COMPUTATIONS - PACK
!               ---------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_RAIN_EVAPORATION', 0, ZHOOK_HANDLE)
P_TH_EVAP(:) = 0.
P_RR_EVAP(:) = 0.
P_CR_EVAP(:) = 0.
!
ZZW1(:) = 0.
ZZW2(:) = 0.
!
GEVAP(:) = .FALSE.
GEVAP(:) = ODCOMPUTE(:)      .AND. &
           PRRT(:)>LIMAP%XRTMIN(3) .AND. &
           PRVT(:)<PRVSAT(:) .AND. &
           PCRT(:)>LIMAP%XCTMIN(3)
!
!
!
IF (LIMAP%LKHKO) THEN

   ZZW1(:) = MAX((1.0 - PRVT(:)/ZZW1(:)),0.0)  ! Subsaturation
    
   ZZW2(:) = 1. / ( CST%XRHOLW*((((PLV(:)/PT(:))**2)/(LIMAW%XTHCO*CST%XRV)) +          & ! G
        (CST%XRV*PT(:))/(LIMAW%XDIVA*PEVSAT(:))))

   ZZW2(:) = 3.0 * LIMAW%XCEVAP * ZZW2(:) * (4.*CST%XPI*CST%XRHOLW/(3.))**(2./3.) *    &
        (PRRT(:))**(1./3.) * (PCRT(:))**(2./3.) * ZZW1(:)                            
   P_RR_EVAP(:) = - ZZW2(:)

   ZZW2(:) = ZZW2(:) * PCRT(:)/PRRT(:)
   P_CR_EVAP = - ZZW2(:)

ELSE

   WHERE ( GEVAP )
!
      ZZW1(:) = MAX((1.0 - PRVT(:)/PRVSAT(:)),0.0)  ! Subsaturation
!
! Compute the function G(T)
!
      ZZW2(:) = 1. / ( CST%XRHOLW*((((PLV(:)/PT(:))**2)/(LIMAW%XTHCO*CST%XRV)) + (CST%XRV*PT(:))/(LIMAW%XDIVA*PEVSAT(:))))
!
! Compute the evaporation tendency
!
      ZZW2(:) = ZZW2(:) * ZZW1(:) * PRRT(:) *        &
           (LIMAW%X0EVAR * PLBDR(:)**LIMAW%XEX0EVAR + LIMAW%X1EVAR * PRHODREF(:)**LIMAW%XEX2EVAR * PLBDR(:)**LIMAW%XEX1EVAR)
      ZZW2(:) = MAX(ZZW2(:),0.0)
!
      P_RR_EVAP(:) = - ZZW2(:)
!   P_TH_EVAP(:) = P_RR_EVAP(:) * PLVFACT(:)
!   PEVAP3D(:) = - P_RR_EVAP(:)
!
   END WHERE
!
END IF
!-----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_RAIN_EVAPORATION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_RAIN_EVAPORATION
END MODULE MODE_LIMA_RAIN_EVAPORATION
