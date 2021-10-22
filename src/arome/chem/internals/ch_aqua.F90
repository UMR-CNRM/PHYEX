!     ######spl
      SUBROUTINE CH_AQUA (TEMP, P, LWC,                         &
                          H2O2, O3, SO2, CO2, HNO3, H2SO4, NH3, &
                          PFRAC, PPH, POX                       )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ###########################################################
!!
!!****  *CH_AQUA*
!!
!!    PURPOSE
!!    -------
!!
!!    empty box for ch_aqua.f90
!!
!!    EXPLICIT ARGUMENTS
!     ------------------
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
!=======================================================================
!=======================================================================
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('CH_AQUA',0,ZHOOK_HANDLE)
  PFRAC(:,1)   = 1.0
  PFRAC(:,2:4) = 0.0
  POX(:)       = 0.0
  IF (LHOOK) CALL DR_HOOK('CH_AQUA',1,ZHOOK_HANDLE)
  RETURN
END SUBROUTINE CH_AQUA
