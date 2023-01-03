!      ###############################
       MODULE MODI_LIMA_INST_PROCS
!      ###############################
!
INTERFACE
   SUBROUTINE LIMA_INST_PROCS (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,             &
                            PEXNREF, PPABST,           &
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,           &
                            PCCT, PCRT, PCIT,                                   &
                            PINT,                 &
                            P_CR_BRKU,                                          & ! spontaneous break up of drops (BRKU) : Nr
                            P_TH_HONR, P_RR_HONR, P_CR_HONR,                    & ! rain drops homogeneous freezing (HONR) : rr, Nr, rg=-rr, th
                            P_TH_IMLT, P_RC_IMLT, P_CC_IMLT,                    & ! ice melting (IMLT) : rc, Nc, ri=-rc, Ni=-Nc, th, IFNF, IFNA
                            PB_TH, PB_RV, PB_RC, PB_RR, PB_RI, PB_RG,           &
                            PB_CC, PB_CR, PB_CI,                                &
                            PB_IFNN)
!
REAL,                 INTENT(IN)    :: PTSTEP     ! Double Time step
CHARACTER(LEN=*),     INTENT(IN)    :: HFMFILE    ! Name of the output FM-file
LOGICAL,              INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
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
   END SUBROUTINE LIMA_INST_PROCS
END INTERFACE
END MODULE MODI_LIMA_INST_PROCS
!     #############################################################################
SUBROUTINE LIMA_INST_PROCS (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,             &
                            PEXNREF, PPABST,           &
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,           &
                            PCCT, PCRT, PCIT,                                   &
                            PINT,                 &
                            P_CR_BRKU,                                          & ! spontaneous break up of drops (BRKU) : Nr
                            P_TH_HONR, P_RR_HONR, P_CR_HONR,                    & ! rain drops homogeneous freezing (HONR) : rr, Nr, rg=-rr, th
                            P_TH_IMLT, P_RC_IMLT, P_CC_IMLT,                    & ! ice melting (IMLT) : rc, Nc, ri=-rc, Ni=-Nc, th, IFNF, IFNA
                            PB_TH, PB_RV, PB_RC, PB_RR, PB_RI, PB_RG,           &
                            PB_CC, PB_CR, PB_CI,                                &
                            PB_IFNN)
!     #############################################################################
!
USE MODD_PARAM_LIMA, ONLY : LCOLD_LIMA, LNUCL_LIMA, LMEYERS_LIMA, LSNOW_LIMA, LWARM_LIMA, LACTI_LIMA, LRAIN_LIMA, LHHONI_LIMA, NMOD_CCN, NMOD_IFN
!
USE MODI_LIMA_DROPS_BREAK_UP
USE MODI_LIMA_DROPS_HOM_FREEZING
USE MODI_LIMA_ICE_MELTING

IMPLICIT NONE



REAL,                 INTENT(IN)    :: PTSTEP     ! Double Time step
CHARACTER(LEN=*),     INTENT(IN)    :: HFMFILE    ! Name of the output FM-file
LOGICAL,              INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
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
!-------------------------------------------------------------------------------
!
IF (LWARM_LIMA .AND. LRAIN_LIMA) THEN
   CALL LIMA_DROPS_BREAK_UP (HFMFILE, OCLOSE_OUT, LDCOMPUTE,    &
                             PCRT, PRRT,                        &
                             P_CR_BRKU,                         &
                             PB_CR                              )
END IF
!
!-------------------------------------------------------------------------------
!
IF (LCOLD_LIMA .AND. LWARM_LIMA .AND. LRAIN_LIMA) THEN
   CALL LIMA_DROPS_HOM_FREEZING (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,   &
                                 PEXNREF, PPABST,                          &
                                 PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                 PCRT,                                     &
                                 P_TH_HONR, P_RR_HONR, P_CR_HONR,          &
                                 PB_TH, PB_RR, PB_CR, PB_RG                )
END IF
!
!-------------------------------------------------------------------------------
!
IF (LCOLD_LIMA .AND. LWARM_LIMA) THEN
   CALL LIMA_ICE_MELTING (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,   &
                          PEXNREF, PPABST,                          &
                          PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                          PCIT, PINT,                               &
                          P_TH_IMLT, P_RC_IMLT, P_CC_IMLT,          &
                          PB_TH, PB_RC, PB_CC, PB_RI, PB_CI, PB_IFNN)
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_INST_PROCS
