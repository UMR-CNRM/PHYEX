!     ######spl
      MODULE MODI_AER_WET_DEP
!!    ################################
!!
!
INTERFACE
!!
SUBROUTINE AER_WET_DEP         (KSPLITR, PTSTEP, PZZ, PRHODREF,       &
                                PRCT, PRRT,                           &
                                PRCS, PRRS,  PSVT, PTHT,              &
                                PPABST, PRGAER, PEVAP3D, KMODE,       &
                                PDENSITY_AER, PMASSMIN, PSEA, PTOWN,  &
                                PCCT, PCRT )
!
IMPLICIT NONE
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step 
                                                   ! integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference [kg/m3] air density
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT    ! Tracer m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCS    ! Cloud water conc derived from source term
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRS    ! Rain water conc derifed from source term
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEVAP3D ! Instantaneous 3D Rain Evaporation flux (KG/KG/S)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT       !Potential temp
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST     ! [Pa] pressure
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRGAER     ! Aerosol radius (um)
INTEGER,                  INTENT(IN)    :: KMODE      ! Nb aerosols mode
REAL,                     INTENT(IN)    :: PDENSITY_AER ! Begin Index for aerosol in cloud
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PMASSMIN   ! Aerosol mass minimum value
REAL, DIMENSION(:,:),OPTIONAL, INTENT(IN)   :: PSEA  ! Sea mask
REAL, DIMENSION(:,:),OPTIONAL, INTENT(IN)   :: PTOWN ! Town mask
REAL, DIMENSION(:,:,:),OPTIONAL, INTENT(IN) :: PCCT   ! Cloud water concentration
REAL, DIMENSION(:,:,:),OPTIONAL, INTENT(IN) :: PCRT   ! Rain water concentration
!
END SUBROUTINE  AER_WET_DEP 
!!
END INTERFACE
END MODULE MODI_AER_WET_DEP
