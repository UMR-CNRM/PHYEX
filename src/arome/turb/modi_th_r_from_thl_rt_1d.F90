!     ######spl
      MODULE MODI_TH_R_FROM_THL_RT_1D
!     ###############################
!
      INTERFACE
      SUBROUTINE TH_R_FROM_THL_RT_1D(HFRAC_ICE,PFRAC_ICE,PP,             &
                                  PTHL, PRT, PTH, PRV, PRL, PRI,         &
                                  PRSATW, PRSATI, PRR, PRS, PRG, PRH     )
CHARACTER*1       , INTENT(IN)    :: HFRAC_ICE
REAL, DIMENSION(:), INTENT(INOUT) :: PFRAC_ICE
REAL, DIMENSION(:), INTENT(IN) :: PP     ! Pressure
REAL, DIMENSION(:), INTENT(IN) :: PTHL   ! Liquid pot. temp.
REAL, DIMENSION(:), INTENT(IN) :: PRT    ! Total mixing ratios 
REAL, DIMENSION(:),OPTIONAL,INTENT(IN) :: PRR, PRS, PRG, PRH
REAL, DIMENSION(:), INTENT(OUT):: PTH    ! Potential temp.
REAL, DIMENSION(:), INTENT(OUT):: PRV    ! vapor mixing ratio
REAL, DIMENSION(:), INTENT(INOUT):: PRL  ! cloud mixing ratio
REAL, DIMENSION(:), INTENT(INOUT):: PRI  ! ice   mixing ratio
REAL, DIMENSION(:), INTENT(OUT)  :: PRSATW ! estimated mixing ration at saturation over water
REAL, DIMENSION(:), INTENT(OUT)  :: PRSATI ! estimated mixing ration at saturation over ice

      END SUBROUTINE TH_R_FROM_THL_RT_1D
      END INTERFACE
      END MODULE MODI_TH_R_FROM_THL_RT_1D
