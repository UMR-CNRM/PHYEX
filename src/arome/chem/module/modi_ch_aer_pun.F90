!     ######spl
     MODULE MODI_CH_AER_PUN
!!   #######################
!!
INTERFACE
!!
SUBROUTINE CH_AER_PUN(   &
     PTOT                & !I [ug/m3] aerosol phase concentrations
     ,PTOTG              & !I [ug/m3] gas phase concentrations
     ,PTEMP              & !I [K] Temperature
     ,PRH                & !I [0-1] relative humidity
     ,PLWC               & !I [ug/m3] liquid water content
     ,PPROTON            & !I [mole/g_{water}] H+ concentration
     ,PTOTNEW            & !O [ug/m3] new aerosol concentration
     ,PTOTGNEW           & !O [ug/m3] new aerosol concentration
     )

IMPLICIT NONE

REAL, DIMENSION(:,:),   INTENT(IN)    :: PTOT     !I [ug/m3] aerosol phase concentrations
REAL, DIMENSION(:,:),   INTENT(IN)    :: PTOTG    !I [ug/m3] gas phase concentrations
REAL, DIMENSION(:),     INTENT(IN)    :: PTEMP    !I [K] temperature
REAL, DIMENSION(:),     INTENT(IN)    :: PRH      !I [0-1] relative humidity
REAL, DIMENSION(:),     INTENT(INOUT) :: PLWC     !I [ug/m3] liquid water content
REAL, DIMENSION(:),     INTENT(INOUT) :: PPROTON  !I [mole/g_{water}] H+ concentrations
REAL, DIMENSION(:,:),   INTENT(OUT)   :: PTOTNEW  !O [ug/m3] aerosol phase concentrations
REAL, DIMENSION(:,:),   INTENT(OUT)   :: PTOTGNEW !O [ug/m3] gas phase concentrations

END SUBROUTINE CH_AER_PUN
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_PUN
