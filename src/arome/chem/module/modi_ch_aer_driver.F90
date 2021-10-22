!     ######spl
     MODULE MODI_CH_AER_DRIVER
!    ######################### 
! 
INTERFACE
! 
SUBROUTINE CH_AER_DRIVER(PM, PSIG0, PRG0, PN0, PCTOTG, PCTOTA,&
                           PCCTOT,  PDTACT, PSEDA,&
                           PMU, PLAMBDA, PRHOP0, POM, PSO4RAT, &
                           PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PMASK,&
                           PTIME, PSOLORG)
IMPLICIT NONE
REAL,                                   INTENT(IN)    :: PDTACT, PTIME
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PRHOP0, POM
REAL,                 DIMENSION(:),     INTENT(INOUT) :: PLAMBDA, PMU, PSO4RAT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSIG0, PRG0, PN0
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSEDA
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PMASK
REAL,                 DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL,                 DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC
END SUBROUTINE CH_AER_DRIVER
! 
END INTERFACE
! 
END MODULE MODI_CH_AER_DRIVER
