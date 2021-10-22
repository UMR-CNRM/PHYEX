!     ######spl
     MODULE MODI_COMPUTE_FUNCTION_THERMO_MF
!    ######################################
!
INTERFACE
      
!     #################################################################
      SUBROUTINE COMPUTE_FUNCTION_THERMO_MF( KRR,KRRL,KRRI,                  &
                                       PTH, PR, PEXN, PFRAC_ICE, PPABS,      &
                                       PT, PAMOIST,PATHETA                   )
!     #################################################################

!*               1.1  Declaration of Arguments
!

INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.

REAL, DIMENSION(:,:), INTENT(IN)   :: PTH      ! theta
REAL, DIMENSION(:,:,:), INTENT(IN) :: PR       ! water species
REAL, DIMENSION(:,:)  , INTENT(IN) :: PPABS,PEXN    ! pressure, Exner funct.
REAL, DIMENSION(:,:)  , INTENT(IN) :: PFRAC_ICE     ! ice fraction

REAL, DIMENSION(:,:), INTENT(OUT)   :: PT      ! temperature

REAL, DIMENSION(:,:), INTENT(OUT)  ::  PAMOIST,PATHETA
!
END SUBROUTINE COMPUTE_FUNCTION_THERMO_MF

END INTERFACE
!
END MODULE MODI_COMPUTE_FUNCTION_THERMO_MF
