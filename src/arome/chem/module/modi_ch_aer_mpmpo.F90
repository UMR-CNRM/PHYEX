!     ######spl
     MODULE MODI_CH_AER_MPMPO
!!   #######################
!!
INTERFACE
!!
SUBROUTINE CH_AER_MPMPO( &
     PTOT                & !I [ug/m3] aerosol phase concentration
     ,PTOTG              & !I/O [ug/m3] gas phase concentrations
     ,PTEMP              & !I [K] temperature
     ,PRH                & !I [0-1] relative humidty
     ,PLWC               & !I [ug/m3] liquid water content
     ,PPROTON            & !I [mole/g_{water}]
     ,PTOTNEW            & !O [ug/m3] new aerosol concentration
     ,PTOTGNEW           & !O [ug/m3] new gas phase concentration
     ,PSOLORG            & !IO [%] Solubility of SOA (fraction)
     )

IMPLICIT NONE

REAL, DIMENSION(:,:),   INTENT(IN)    :: PTOT  ![ug/m3] aerosol conc
REAL, DIMENSION(:,:),   INTENT(IN)    :: PTOTG ![ug/m3] gas concentration
REAL, DIMENSION(:),     INTENT(IN)    :: PTEMP, PRH ![K/-] temp and RH
REAL, DIMENSION(:),     INTENT(IN)    :: PLWC  ![ug/m3] liquid water
REAL, DIMENSION(:),     INTENT(IN)    :: PPROTON ![mole/g_{water]] H+
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PTOTNEW ![ug/m3] new aerosol conc
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PTOTGNEW ![ug/m3] new gas conc
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG ! [%] Solubility of SOA (fraction)


END SUBROUTINE CH_AER_MPMPO
!
END INTERFACE
!
END MODULE MODI_CH_AER_MPMPO
