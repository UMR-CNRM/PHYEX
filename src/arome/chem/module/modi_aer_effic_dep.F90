!     ######spl
      MODULE MODI_AER_EFFIC_DEP
!!    ########################
!!
!
INTERFACE
!!
SUBROUTINE AER_EFFIC_DEP(PRG,PVGG,  & !aerosol radius/fall speed (m/s)
                PRHODREF,           & !Air density     
                PMUW, PMU,          & !mu water/air
                PDPG, PEFC,         & !diffusivity, efficiency
                PRRS,               & ! Rain water m.r. at time 
                KMODE,              & ! Number of aerosol modes
                PTEMP, PCOR,        & ! air temp, cunningham corr factor
                PDENSITY_AER )        ! aerosol density
!
IMPLICIT NONE
REAL, DIMENSION(:,:), INTENT(IN) ::  PRG,  PVGG
REAL, DIMENSION(:),   INTENT(IN) ::  PRHODREF
REAL, DIMENSION(:,:), INTENT(IN) ::  PDPG
REAL, DIMENSION(:),   INTENT(IN) ::  PMU, PMUW
REAL, DIMENSION(:,:), INTENT(INOUT) :: PEFC
REAL, DIMENSION(:),   INTENT(IN)    :: PRRS
REAL, DIMENSION(:),   INTENT(IN)    :: PTEMP
REAL, DIMENSION(:,:), INTENT(IN)    :: PCOR
INTEGER, INTENT(IN)                 :: KMODE
REAL, INTENT(IN)                    :: PDENSITY_AER


END SUBROUTINE  AER_EFFIC_DEP 
!!
END INTERFACE
END MODULE MODI_AER_EFFIC_DEP
