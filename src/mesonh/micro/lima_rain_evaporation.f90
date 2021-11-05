!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ##########################
       MODULE MODI_LIMA_RAIN_EVAPORATION
!      ##########################
!
INTERFACE
      SUBROUTINE LIMA_RAIN_EVAPORATION (PTSTEP, LDCOMPUTE,                          &
                                        PRHODREF, PT, PLV, PLVFACT, PEVSAT, PRVSAT, &
                                        PRVT, PRCT, PRRT, PLBDR,                    &
                                        P_TH_EVAP, P_RR_EVAP,                       &
                                        PEVAP3D                                     )
!
REAL,                 INTENT(IN)    :: PTSTEP     ! Time step
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE  !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:),   INTENT(IN)    :: PT         ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PLV        ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PEVSAT     !
REAL, DIMENSION(:),   INTENT(IN)    :: PRVSAT     !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR     ! Lambda(rain)
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_EVAP
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RR_EVAP
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PEVAP3D    ! Rain evap profile
!
END SUBROUTINE LIMA_RAIN_EVAPORATION
END INTERFACE
END MODULE MODI_LIMA_RAIN_EVAPORATION
!     ###############################################################################
      SUBROUTINE LIMA_RAIN_EVAPORATION (PTSTEP, LDCOMPUTE,                          &
                                        PRHODREF, PT, PLV, PLVFACT, PEVSAT, PRVSAT, &
                                        PRVT, PRCT, PRRT, PLBDR,                    &
                                        P_TH_EVAP, P_RR_EVAP,                       &
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,             ONLY : XRHOLW, XRV
USE MODD_PARAM_LIMA,      ONLY : XRTMIN
USE MODD_PARAM_LIMA_WARM, ONLY : X0EVAR, XEX0EVAR, X1EVAR, XEX2EVAR, XEX1EVAR, XTHCO, XDIVA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                 INTENT(IN)    :: PTSTEP     ! Time step
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE  !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF   ! Reference density
REAL, DIMENSION(:),   INTENT(IN)    :: PT         ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PLV        ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT    !
REAL, DIMENSION(:),   INTENT(IN)    :: PEVSAT     !
REAL, DIMENSION(:),   INTENT(IN)    :: PRVSAT     !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR     ! Lambda(rain)
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_EVAP
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RR_EVAP
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PEVAP3D    ! Rain evap profile
!
!*       0.1   Declarations of local variables :
!
! 
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GEVAP
REAL, DIMENSION(SIZE(PRHODREF))    :: ZZW1, ZZW2
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PREPARE COMPUTATIONS - PACK
!   	        ---------------------------
!
P_TH_EVAP(:) = 0.
P_RR_EVAP(:) = 0.
!
GEVAP(:) = .FALSE.
GEVAP(:) = LDCOMPUTE(:)      .AND. &
           PRRT(:)>XRTMIN(3) .AND. &
           PRVT(:)<PRVSAT(:)
!
WHERE ( GEVAP )
!
!-------------------------------------------------------------------------------
!
!
!*       2. compute the evaporation of rain drops    
!   	 ----------------------------------------
!
!
   ZZW1(:) = MAX((1.0 - PRVT(:)/PRVSAT(:)),0.0)  ! Subsaturation
!
! Compute the function G(T)
!
   ZZW2(:) = 1. / ( XRHOLW*((((PLV(:)/PT(:))**2)/(XTHCO*XRV)) +          & ! G
          (XRV*PT(:))/(XDIVA*PEVSAT(:))))
!
! Compute the evaporation tendency
!
   ZZW2(:) = ZZW2(:) * ZZW1(:) * PRRT(:) *        &
        (X0EVAR * PLBDR(:)**XEX0EVAR + X1EVAR * PRHODREF(:)**XEX2EVAR * PLBDR(:)**XEX1EVAR)
   ZZW2(:) = MAX(ZZW2(:),0.0)
!
   P_RR_EVAP(:) = - ZZW2(:)
!   P_TH_EVAP(:) = P_RR_EVAP(:) * PLVFACT(:)
!   PEVAP3D(:) = - P_RR_EVAP(:)
!
END WHERE
!
!-----------------------------------------------------------------------------
!
END SUBROUTINE LIMA_RAIN_EVAPORATION
