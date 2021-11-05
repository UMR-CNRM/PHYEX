!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!      ###############################
       MODULE MODI_LIMA_INST_PROCS
!      ###############################
!
INTERFACE
   SUBROUTINE LIMA_INST_PROCS (PTSTEP, LDCOMPUTE,                        &
                               PEXNREF, PPABST,                          &
                               PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                               PCCT, PCRT, PCIT,                         &
                               PINT,                                     &
                               P_CR_BRKU,                                & ! spontaneous break up of drops (BRKU) : Nr
                               P_TH_HONR, P_RR_HONR, P_CR_HONR,          & ! rain drops homogeneous freezing (HONR) : rr, Nr, rg=-rr, th
                               P_TH_IMLT, P_RC_IMLT, P_CC_IMLT,          & ! ice melting (IMLT) : rc, Nc, ri=-rc, Ni=-Nc, th, IFNF, IFNA
                               PB_TH, PB_RV, PB_RC, PB_RR, PB_RI, PB_RG, &
                               PB_CC, PB_CR, PB_CI,                      &
                               PB_IFNN,                                  &
                               PCF1D, PIF1D, PPF1D                       )
!
REAL,                 INTENT(IN)    :: PTSTEP     ! Time step
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(:),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PTHT       ! Theta at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT       ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRST       ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT       ! Rain water m.r. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT       ! Cloud water conc. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT       ! Rain water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT       ! Prinstine ice conc. at t
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PINT       ! IFN C. activated at t
!
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_CR_BRKU  ! Concentration change (#/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_TH_HONR  ! 
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_RR_HONR  ! mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_CR_HONR  ! Concentration change (#/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_TH_IMLT  ! 
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_RC_IMLT  ! mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_CC_IMLT  ! Concentration change (#/kg)
!
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_TH      ! Cumulated theta change
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_RV      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_RC      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_RR      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_RI      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_RG      ! Cumulated mr change (kg/kg)
!
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_CC      ! Cumulated concentration change (#/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_CR      ! Cumulated concentration change (#/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_CI      ! Cumulated concentration change (#/kg)
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PB_IFNN    ! Cumulated concentration change (#/kg)
!
REAL, DIMENSION(:)  , INTENT(INOUT) :: PCF1D      ! Liquid cloud fraction
REAL, DIMENSION(:)  , INTENT(INOUT) :: PIF1D      ! Ice cloud fraction
REAL, DIMENSION(:)  , INTENT(INOUT) :: PPF1D      ! Precipitation fraction
!
   END SUBROUTINE LIMA_INST_PROCS
END INTERFACE
END MODULE MODI_LIMA_INST_PROCS
!
!
!     ###########################################################################
SUBROUTINE LIMA_INST_PROCS (PTSTEP, LDCOMPUTE,                                  &
                            PEXNREF, PPABST,                                    &
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,           &
                            PCCT, PCRT, PCIT,                                   &
                            PINT,                                               &
                            P_CR_BRKU,                                          & ! spontaneous break up of drops (BRKU) : Nr
                            P_TH_HONR, P_RR_HONR, P_CR_HONR,                    & ! rain drops homogeneous freezing (HONR) : rr, Nr, rg=-rr, th
                            P_TH_IMLT, P_RC_IMLT, P_CC_IMLT,                    & ! ice melting (IMLT) : rc, Nc, ri=-rc, Ni=-Nc, th, IFNF, IFNA
                            PB_TH, PB_RV, PB_RC, PB_RR, PB_RI, PB_RG,           &
                            PB_CC, PB_CR, PB_CI,                                &
                            PB_IFNN,                                            &
                            PCF1D, PIF1D, PPF1D                                 )
!     ###########################################################################
!
!!    PURPOSE
!!    -------
!!      Compute sources of instantaneous microphysical processes for the
!!    time-split version of LIMA
!!
!!    AUTHOR
!!    ------
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!!
!-------------------------------------------------------------------------------
!
!
USE MODD_PARAM_LIMA, ONLY : LCOLD, LWARM, LRAIN
!
USE MODI_LIMA_DROPS_BREAK_UP
USE MODI_LIMA_DROPS_HOM_FREEZING
USE MODI_LIMA_ICE_MELTING

IMPLICIT NONE



REAL,                 INTENT(IN)    :: PTSTEP     ! Time step
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(:),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PTHT       ! Theta at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT       ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRST       ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT       ! Rain water m.r. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT       ! Cloud water conc. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT       ! Rain water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT       ! Prinstine ice conc. at t
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PINT       ! IFN C. activated at t
!
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_CR_BRKU  ! Concentration change (#/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_TH_HONR  ! 
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_RR_HONR  ! mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_CR_HONR  ! Concentration change (#/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_TH_IMLT  ! 
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_RC_IMLT  ! mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: P_CC_IMLT  ! Concentration change (#/kg)
!
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_TH      ! Cumulated theta change
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_RV      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_RC      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_RR      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_RI      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_RG      ! Cumulated mr change (kg/kg)
!
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_CC      ! Cumulated concentration change (#/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_CR      ! Cumulated concentration change (#/kg)
REAL, DIMENSION(:)  , INTENT(INOUT) :: PB_CI      ! Cumulated concentration change (#/kg)
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PB_IFNN    ! Cumulated concentration change (#/kg)
!
REAL, DIMENSION(:)  , INTENT(INOUT) :: PCF1D      ! Liquid cloud fraction
REAL, DIMENSION(:)  , INTENT(INOUT) :: PIF1D      ! Ice cloud fraction
REAL, DIMENSION(:)  , INTENT(INOUT) :: PPF1D      ! Precipitation fraction
!
!-------------------------------------------------------------------------------
!
IF (LWARM .AND. LRAIN) THEN
   CALL LIMA_DROPS_BREAK_UP (LDCOMPUTE,    & ! no dependance on CF, IF or PF
                             PCRT, PRRT,   &
                             P_CR_BRKU,    &
                             PB_CR         )
END IF
!
!-------------------------------------------------------------------------------
!
IF (LCOLD .AND. LWARM .AND. LRAIN) THEN
   CALL LIMA_DROPS_HOM_FREEZING (PTSTEP, LDCOMPUTE,                        & ! no dependance on CF, IF or PF
                                 PEXNREF, PPABST,                          &
                                 PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                 PCRT,                                     &
                                 P_TH_HONR, P_RR_HONR, P_CR_HONR,          &
                                 PB_TH, PB_RR, PB_CR, PB_RG                )
END IF
!
!-------------------------------------------------------------------------------
!
IF (LCOLD .AND. LWARM) THEN
   CALL LIMA_ICE_MELTING (PTSTEP, LDCOMPUTE,                        & ! no dependance on CF, IF or PF
                          PEXNREF, PPABST,                          & ! but ice fraction becomes cloud fraction
                          PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, & ! -> where ?
                          PCIT, PINT,                               &
                          P_TH_IMLT, P_RC_IMLT, P_CC_IMLT,          &
                          PB_TH, PB_RC, PB_CC, PB_RI, PB_CI, PB_IFNN)
   !
   !PCF1D(:)=MAX(PCF1D(:),PIF1D(:))
   !PIF1D(:)=0.
   !
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_INST_PROCS
