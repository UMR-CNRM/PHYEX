!     ######spl
     MODULE MODI_TURB_VER_SV_FLUX 
!    ####################
!
INTERFACE 
!
      SUBROUTINE TURB_VER_SV_FLUX(KKA,KKU,KKL,                      &
                      OCLOSE_OUT,OTURB_FLX,HTURBDIM,                &
                      PIMPL,PEXPL,                                  &
                      PTSTEP,                                       &
                      HFMFILE,HLUOUT,                               &
                      PDZZ,PDIRCOSZW,                               &
                      PRHODJ,PWM,                                   &
                      PSFSVM,PSFSVP,                                &
                      PSVM,                                         &
                      PTKEM,PLM,MFMOIST,PPSI_SV,                    &
                      PRSVS,PWSV                                    )
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=AR
LOGICAL,                INTENT(IN)   ::  OCLOSE_OUT   ! switch for syncronous
                                                      ! file opening       
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
CHARACTER*4,            INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP       ! Double Time Step
CHARACTER(LEN=*),       INTENT(IN)   ::  HFMFILE      ! Name of the output
                                                      ! FM-file 
CHARACTER(LEN=*),       INTENT(IN)   ::  HLUOUT       ! Output-listing name for
                                                      ! model n
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDZZ 
                                                      ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  MFMOIST      ! moist mass flux dual scheme

!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVM       ! t - deltat 
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVP       ! t + deltat 
!
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PWM          ! vertical wind
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PPSI_SV      ! Inv.Turb.Sch.for scalars
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS
                            ! cumulated sources for the prognostic variables
REAL, DIMENSION(:,:,:,:), INTENT(OUT)  :: PWSV        ! scalar flux

!
!
END SUBROUTINE TURB_VER_SV_FLUX
!
END INTERFACE
!
END MODULE MODI_TURB_VER_SV_FLUX
