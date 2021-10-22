!     ######spl
      MODULE MODI_CH_MONITOR_n
!!    ########################
!!
!
INTERFACE
!!
SUBROUTINE CH_MONITOR_n(KTCOUNT,PTSTEP, KLUOUT, KVERB)
IMPLICIT NONE
INTEGER, INTENT(IN) :: KTCOUNT    ! iteration count
REAL,  INTENT(IN)   ::  PTSTEP    ! Double timestep except 
                                  ! for the first time step (single one)
INTEGER, INTENT(IN) :: KLUOUT     ! unit for output listing count
INTEGER, INTENT(IN) :: KVERB      ! verbosity level
END SUBROUTINE CH_MONITOR_n
!!
END INTERFACE
!!
END MODULE MODI_CH_MONITOR_n
