!     ######spl
      MODULE MODI_CH_AQUA
!!    ###################
!!
INTERFACE
!!
!!    ###########################################################
      SUBROUTINE CH_AQUA (TEMP, P, LWC,                         &
                          H2O2, O3, SO2, CO2, HNO3, H2SO4, NH3, &
                          PFRAC, PPH, POX                       )
!!    ###########################################################
IMPLICIT NONE
REAL, INTENT(IN) :: TEMP ! air temperature
REAL, INTENT(IN) :: P    ! pressure (in atm)
REAL, INTENT(IN) :: LWC  ! cloud liquid water mixing ratio (kg/kg)
REAL, INTENT(IN) :: H2O2, O3, SO2, CO2, HNO3, H2SO4, NH3
                         ! total concentration of these chemical species
REAL, INTENT(INOUT)               :: PPH        ! pH of the droplets
REAL, INTENT(OUT), DIMENSION(3,4) :: PFRAC      ! fraction of each dissolved
                         ! species in the four states
REAL, INTENT(OUT), DIMENSION(2)   :: POX        ! SO2 oxidation by O3
                                                ! and H2O2
END SUBROUTINE CH_AQUA
!!
END INTERFACE
!!
END MODULE MODI_CH_AQUA
