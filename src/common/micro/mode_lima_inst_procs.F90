!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_INST_PROCS
  IMPLICIT NONE
CONTAINS
!     ###########################################################################
  SUBROUTINE LIMA_INST_PROCS (CST, LIMAP, LIMAW, KSIZE, PTSTEP, ODCOMPUTE,                           &
                              PEXNREF, PPABST,                                    &
                              PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,           &
                              PCCT, PCRT, PCIT, PCIT_SHAPE,                       &
                              PINT,                                               &
                              P_CR_BRKU,                                          & ! spontaneous break up of drops (BRKU) : Nr
                              P_TH_HONR, P_RR_HONR, P_CR_HONR,                    & ! rain drops homogeneous freezing (HONR) : rr, Nr, rg=-rr, th
                              P_TH_IMLT, P_RC_IMLT, P_CC_IMLT, P_SHCI_IMLT,       & ! ice melting (IMLT) : rc, Nc, ri=-rc, Ni=-Nc, th, IFNF, IFNA
                              PB_TH, PB_RV, PB_RC, PB_RR, PB_RI, PB_RG,           &
                              PB_CC, PB_CR, PB_CI, PB_CI_SHAPE, PB_CG,            &
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
!
USE MODE_LIMA_DROPS_BREAK_UP, ONLY: LIMA_DROPS_BREAK_UP
USE MODE_LIMA_DROPS_HOM_FREEZING, ONLY: LIMA_DROPS_HOM_FREEZING
USE MODE_LIMA_ICE_MELTING, ONLY: LIMA_ICE_MELTING
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(CST_T),INTENT(IN)::CST
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
INTEGER,              INTENT(IN)    :: KSIZE
REAL,                 INTENT(IN)    :: PTSTEP     ! Time step
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PEXNREF    ! Reference Exner function
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PPABST     ! abs. pressure at time t
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PTHT       ! Theta at t 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRVT       ! Water vapor m.r. at t 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRIT       ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRST       ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRGT       ! Rain water m.r. at t
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCCT       ! Cloud water conc. at t 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCRT       ! Rain water conc. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCIT       ! Prinstine ice conc. at t
REAL, DIMENSION(KSIZE, LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN) :: PCIT_SHAPE  ! 
!
REAL, DIMENSION(KSIZE,LIMAP%NMOD_IFN), INTENT(IN)    :: PINT       ! IFN C. activated at t
!
REAL, DIMENSION(KSIZE)  , INTENT(OUT) :: P_CR_BRKU  ! Concentration change (#/kg)
REAL, DIMENSION(KSIZE)  , INTENT(OUT) :: P_TH_HONR  ! 
REAL, DIMENSION(KSIZE)  , INTENT(OUT) :: P_RR_HONR  ! mr change (kg/kg)
REAL, DIMENSION(KSIZE)  , INTENT(OUT) :: P_CR_HONR  ! Concentration change (#/kg)
REAL, DIMENSION(KSIZE)  , INTENT(OUT) :: P_TH_IMLT  ! 
REAL, DIMENSION(KSIZE)  , INTENT(OUT) :: P_RC_IMLT  ! mr change (kg/kg)
REAL, DIMENSION(KSIZE)  , INTENT(OUT) :: P_CC_IMLT  ! Concentration change (#/kg)
REAL, DIMENSION(KSIZE, LIMAP%NNB_CRYSTAL_SHAPE), INTENT(OUT) :: P_SHCI_IMLT  ! 
!
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PB_TH      ! Cumulated theta change
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PB_RV      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PB_RC      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PB_RR      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PB_RI      ! Cumulated mr change (kg/kg)
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PB_RG      ! Cumulated mr change (kg/kg)
!
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PB_CC      ! Cumulated concentration change (#/kg)
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PB_CR      ! Cumulated concentration change (#/kg)
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PB_CI      ! Cumulated concentration change (#/kg)
REAL, DIMENSION(KSIZE, LIMAP%NNB_CRYSTAL_SHAPE), INTENT(INOUT) :: PB_CI_SHAPE  ! 
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PB_CG      ! Cumulated concentration change (#/kg)
!
REAL, DIMENSION(KSIZE,LIMAP%NMOD_IFN), INTENT(INOUT) :: PB_IFNN    ! Cumulated concentration change (#/kg)
!
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PCF1D      ! Liquid cloud fraction
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PIF1D      ! Ice cloud fraction
REAL, DIMENSION(KSIZE)  , INTENT(INOUT) :: PPF1D      ! Precipitation fraction
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_INST_PROCS', 0, ZHOOK_HANDLE)
IF (LIMAP%NMOM_R.GE.2) THEN
   CALL LIMA_DROPS_BREAK_UP (LIMAP, LIMAW, KSIZE, ODCOMPUTE, & ! no dependance on CF, IF or PF
                             PCRT, PRRT,       &
                             P_CR_BRKU,        &
                             PB_CR             )
END IF
!
!-------------------------------------------------------------------------------
!
IF (LIMAP%NMOM_G.GE.1 .AND. LIMAP%NMOM_R.GE.1) THEN
   CALL LIMA_DROPS_HOM_FREEZING (CST, LIMAP, KSIZE, PTSTEP, ODCOMPUTE,                 & ! no dependance on CF, IF or PF
                                 PEXNREF, PPABST,                          &
                                 PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                 PCRT,                                     &
                                 P_TH_HONR, P_RR_HONR, P_CR_HONR,          &
                                 PB_TH, PB_RR, PB_CR, PB_RG, PB_CG         )
END IF
!
!-------------------------------------------------------------------------------
!
IF (LIMAP%NMOM_C.GE.1 .AND. LIMAP%NMOM_I.GE.1) THEN
   CALL LIMA_ICE_MELTING (CST, LIMAP, KSIZE, PTSTEP, ODCOMPUTE,                 & ! no dependance on CF, IF or PF
                          PEXNREF, PPABST,                          & ! but ice fraction becomes cloud fraction
                          PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, & ! -> where ?
                          PCIT, PCIT_SHAPE, PINT,                               &
                          P_TH_IMLT, P_RC_IMLT, P_CC_IMLT, P_SHCI_IMLT,         &
                          PB_TH, PB_RC, PB_CC, PB_RI, PB_CI, PB_CI_SHAPE, PB_IFNN)
   !
   !PCF1D(:)=MAX(PCF1D(:),PIF1D(:))
   !PIF1D(:)=0.
   !
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_INST_PROCS', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_INST_PROCS
END MODULE MODE_LIMA_INST_PROCS
