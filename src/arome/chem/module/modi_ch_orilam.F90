!     ######spl
     MODULE MODI_CH_ORILAM
!!   ######################### 
!!
INTERFACE
!!
SUBROUTINE CH_ORILAM(PAERO, PCHEM, PM, PSIG0, PRG0, PN0, PCTOTG, PCTOTA,&
                           PCCTOT, PDTACT, PSEDA,&
                           PMU, PLAMBDA, PRHOP0, POM, PSO4RAT,  &
                           PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PFRAC, PMI,&
                           PTIME, GSCHEME, PSOLORG)
IMPLICIT NONE
REAL,                                   INTENT(IN)    :: PDTACT, PTIME
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PRHOP0, POM
REAL,                 DIMENSION(:),     INTENT(INOUT) :: PLAMBDA, PMU, PSO4RAT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSIG0, PRG0, PN0
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSEDA
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCHEM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PAERO
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PFRAC
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PMI
REAL,                 DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG
REAL,                 DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC
CHARACTER(LEN=10),                      INTENT(IN)    :: GSCHEME


END SUBROUTINE CH_ORILAM
!!
END INTERFACE
!!
END MODULE MODI_CH_ORILAM
