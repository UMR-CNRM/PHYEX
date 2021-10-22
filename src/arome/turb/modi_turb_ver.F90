!     ######spl
     MODULE MODI_TURB_VER 
!    ####################
!
INTERFACE 
!
      SUBROUTINE TURB_VER(KKA,KKU,KKL,KRR,KRRL,KRRI,                &
                      OCLOSE_OUT,OTURB_FLX,                         &
                      HTURBDIM,HTOM,PIMPL,PEXPL,                    & 
                      PTSTEP_UVW,PTSTEP_MET, PTSTEP_SV,             &
                      HFMFILE,HLUOUT,                               &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PDIRCOSZW,PZZ,       &
                      PCOSSLOPE,PSINSLOPE,                          &
                      PRHODJ,PTHVREF,                               &
                      PSFTHM,PSFRM,PSFSVM,PSFTHP,PSFRP,PSFSVP,      &
                      PCDUEFF,PTAU11M,PTAU12M,PTAU33M,              &
                      PUM,PVM,PWM,PUSLOPEM,PVSLOPEM,PTHLM,PRM,PSVM, &
                      PTKEM,PLM,PLENGTHM,PLENGTHH,PLEPS,MFMOIST,    &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PFRAC_ICE,    &
                      PFWTH,PFWR,PFTH2,PFR2,PFTHR,PBL_DEPTH,        &
                      PSBL_DEPTH,PLMO,                              &
                      PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,             &
                      PDP,PTP,PSIGS,PWTH,PWRC,PWSV          )

!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=ARO
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
LOGICAL,                INTENT(IN)   ::  OCLOSE_OUT   ! switch for syncronous
                                                      ! file opening       
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
CHARACTER*4,            INTENT(IN)   ::  HTURBDIM     ! dimensionality of the 
                                                      ! turbulence scheme
CHARACTER*4,            INTENT(IN)   ::  HTOM         ! type of Third Order Moment
REAL,                   INTENT(IN)   ::  PIMPL, PEXPL ! Coef. for temporal disc.
REAL,                   INTENT(IN)   ::  PTSTEP_UVW   ! Dynamical timestep 
REAL,                   INTENT(IN)   ::  PTSTEP_MET   ! Timestep for meteorological variables                        
REAL,                   INTENT(IN)   ::  PTSTEP_SV    ! Timestep for tracer variables
CHARACTER(LEN=*),       INTENT(IN)   ::  HFMFILE      ! Name of the output
                                                      ! FM-file 
CHARACTER(LEN=*),       INTENT(IN)   ::  HLUOUT       ! Output-listing name for
                                                      ! model n
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                      ! Metric coefficients
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDIRCOSZW    ! Director Cosinus of the
                                                      ! normal to the ground surface
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PZZ          ! altitudes at flux points
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCOSSLOPE    ! cosinus of the angle 
                                      ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSINSLOPE    ! sinus of the angle 
                                      ! between i and the slope vector
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRHODJ       ! dry density * grid volum
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  MFMOIST      ! moist mass flux dual scheme

REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHVREF      ! ref. state Virtual 
                                                      ! Potential Temperature 
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSFTHM,PSFRM ! surface fluxes at time
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVM       ! t - deltat 
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSFTHP,PSFRP ! surface fluxes at time
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSFSVP       ! t + deltat 
!
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCDUEFF     ! Cd * || u || at time t
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU11M      ! <uu> in the axes linked 
       ! to the maximum slope direction and the surface normal and the binormal 
       ! at time t - dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU12M      ! <uv> in the same axes
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTAU33M      ! <ww> in the same axes
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PUM,PVM,PWM,PTHLM 
  ! Wind and potential temperature at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRM          ! Mixing ratios 
                                                      ! at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PUSLOPEM     ! wind component along the 
                                     ! maximum slope direction
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PVSLOPEM     ! wind component along the 
                                     ! direction normal to the maximum slope one
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLENGTHM     ! Turb. mixing length momentum   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLENGTHH     ! Turb. mixing length heat/moisture  

REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLEPS        ! dissipative length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exner at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFRAC_ICE    ! ri fraction of rc+ri
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFWTH        ! d(w'2th' )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFWR         ! d(w'2r'  )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFTH2        ! d(w'th'2 )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFR2         ! d(w'r'2  )/dz
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PFTHR        ! d(w'th'r')/dz
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PBL_DEPTH    ! BL depth
REAL, DIMENSION(:,:),   INTENT(INOUT)::  PSBL_DEPTH   ! SBL depth
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PLMO         ! Monin-Obukhov length
!
REAL, DIMENSION(:,:,:), INTENT(INOUT)   ::  PRUS, PRVS, PRWS, PRTHLS
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS,PRRS
                            ! cumulated sources for the prognostic variables
!
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PDP,PTP   ! Dynamic and thermal
                                                   ! TKE production terms
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PSIGS     ! Vert. part of Sigma_s at t
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWTH       ! heat flux
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWRC       ! cloud water flux
REAL, DIMENSION(:,:,:,:),INTENT(OUT) :: PWSV       ! scalar flux

!
!
END SUBROUTINE TURB_VER
!
END INTERFACE
!
END MODULE MODI_TURB_VER
