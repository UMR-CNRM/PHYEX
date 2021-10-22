!     ######spl
     MODULE MODI_SEDIM_DUST
!!   ##############################
!!
INTERFACE
!
SUBROUTINE SEDIM_DUST(  &
     PTHT               & !I [K] theta
     ,PDTMONITOR        & !I Time step
     ,PRHODREF          & !I [kg/m3] air density
     ,PPABST            & !I [Pa] pressure
     ,PZZ               & !I [m] height of layers
     ,PSVT              & !IO [scalar variable, ppp] sea salt concentration
     )

IMPLICIT NONE

REAL,  INTENT(IN) :: PDTMONITOR
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSVT   !scalar variable 
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PTHT,PRHODREF, PZZ
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PPABST


END SUBROUTINE SEDIM_DUST
!!
END INTERFACE
!!
END MODULE MODI_SEDIM_DUST
