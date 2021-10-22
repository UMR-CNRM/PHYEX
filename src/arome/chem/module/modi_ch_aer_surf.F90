!     ######spl
     MODULE MODI_CH_AER_SURF
!!   ########################
!!
INTERFACE
!!
SUBROUTINE CH_AER_SURF(PM, PRG, PSIG, PSURF)

IMPLICIT NONE
REAL,   DIMENSION(:,:), INTENT(IN) :: PM      ! moments
REAL,   DIMENSION(:,:), INTENT(IN) :: PRG     ! radius
REAL,   DIMENSION(:,:), INTENT(IN) :: PSIG    ! dispersion
REAL,   DIMENSION(:,:), INTENT(OUT) :: PSURF  ! aerosol surface
END SUBROUTINE CH_AER_SURF
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_SURF
