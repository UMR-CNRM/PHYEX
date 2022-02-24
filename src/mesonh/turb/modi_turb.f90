!     ######spl
     MODULE MODI_TURB  
!    ################ 
!
INTERFACE
!
!
      SUBROUTINE TURB(KKA,KKU,KKL,KMI,KRR,KRRL,KRRI,HLBCX,HLBCY,      &
              & KSPLIT,KMODEL_CL,                                     &
              & OTURB_FLX,OTURB_DIAG,OSUBG_COND,ORMC01,               &
              & HTURBDIM,HTURBLEN,HTOM,HTURBLEN_CL,HCLOUD,PIMPL,      &
              & PTSTEP,TPFILE,PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,           &
              & PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,PCOSSLOPE,PSINSLOPE,    &
              & PRHODJ,PTHVREF,                                       &
              & PSFTH,PSFRV,PSFSV,PSFU,PSFV,                          &
              & PPABST,PUT,PVT,PWT,PTKET,PSVT,PSRCT,                  &
              & PLENGTHM,PLENGTHH,MFMOIST,                            &
              & PBL_DEPTH, PSBL_DEPTH,                                &
              & PCEI,PCEI_MIN,PCEI_MAX,PCOEF_AMPL_SAT,                &
              & PTHLT,PRT,                                            &
              & PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,PRTKES,PRTKEMS,PSIGS,&
              & PFLXZTHVMF,PWTH,PWRC,PWSV,PDYP,PTHP,PTDIFF,PTDISS,    &
              & TBUDGETS, KBUDGETS                                   ,&
              & PTR,PDISS,PEDR,PLEM                            )

USE MODD_BUDGET, ONLY : TBUDGETDATA
USE MODD_IO, ONLY : TFILEDATA
!
!
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=ARO
INTEGER,                INTENT(IN)   :: KMI           ! model index number  
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
CHARACTER(LEN=*),DIMENSION(:),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
INTEGER,                INTENT(IN)   :: KSPLIT        ! number of time-splitting
INTEGER,                INTENT(IN)   :: KMODEL_CL     ! model number for cloud mixing length
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                 ! diagnostic fields in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OSUBG_COND   ! switch for SUBGrid 
                                 ! CONDensation
LOGICAL,                INTENT(IN)   ::  ORMC01       ! switch for RMC01 lengths in SBL
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBDIM     ! dimensionality of the
                                                      ! turbulence scheme
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBLEN     ! kind of mixing length
CHARACTER(LEN=4),       INTENT(IN)   ::  HTOM         ! kind of Third Order Moment
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBLEN_CL  ! kind of cloud mixing length
REAL,                   INTENT(IN)   ::  PIMPL        ! degree of implicitness
CHARACTER (LEN=4),      INTENT(IN)   ::  HCLOUD       ! Kind of microphysical scheme
REAL,                   INTENT(IN)   ::  PTSTEP       ! timestep 
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(:,:,:), INTENT(IN)   :: PDXX,PDYY,PDZZ,PDZX,PDZY
                                        ! metric coefficients
REAL, DIMENSION(:,:,:), INTENT(IN)   :: PZZ       !  physical distance 
! between 2 succesive grid points along the K direction
REAL, DIMENSION(:,:),   INTENT(IN)      ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PRHODJ    ! dry density * Grid size
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  MFMOIST ! moist mass flux dual scheme
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PTHVREF   ! Virtual Potential
                                        ! Temperature of the reference state
!
REAL, DIMENSION(:,:),   INTENT(IN)      ::  PSFTH,PSFRV,   &
! normal surface fluxes of theta and Rv 
                                            PSFU,PSFV
! normal surface fluxes of (u,v) parallel to the orography
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PSFSV
! normal surface fluxes of Scalar var. 
!
!    prognostic variables at t- deltat
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PPABST      ! Pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PUT,PVT,PWT ! wind components
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PTKET       ! TKE
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVT        ! passive scal. var.
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PSRCT       ! Second-order flux
                      ! s'rc'/2Sigma_s2 at time t-1 multiplied by Lambda_3
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PBL_DEPTH  ! BL height for TOMS
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PSBL_DEPTH ! SBL depth for RMC01
!
!    variables for cloud mixing length
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PCEI ! Cloud Entrainment instability
                                                 ! index to emphasize localy 
                                                 ! turbulent fluxes
REAL, INTENT(IN)      ::  PCEI_MIN ! minimum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCEI_MAX ! maximum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCOEF_AMPL_SAT ! saturation of the amplification coefficient
!
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PTHLT       ! conservative pot. temp.
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRT         ! water var.  where 
                             ! PRT(:,:,:,1) is the conservative mixing ratio        
!
! sources of momentum, conservative potential temperature, Turb. Kin. Energy, 
! TKE dissipation
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRUS,PRVS,PRWS,PRTHLS,PRTKES
! Source terms for all water kinds, PRRS(:,:,:,1) is used for the conservative
! mixing ratio
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRTKEMS
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRRS 
! Source terms for all passive scalar variables
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS
! Sigma_s at time t+1 : square root of the variance of the deviation to the 
! saturation
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PSIGS
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PFLXZTHVMF 
!                                           MF contribution for vert. turb. transport
!                                           used in the buoy. prod. of TKE
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWTH       ! heat flux
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWRC       ! cloud water flux
REAL, DIMENSION(:,:,:,:),INTENT(OUT) :: PWSV       ! scalar flux
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PDYP  ! Dynamical production of TKE
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PTHP  ! Thermal production of TKE
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PTDIFF   ! Transport production of TKE
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PTDISS ! Dissipation of TKE
!
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
! length scale from vdfexcu
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLENGTHM, PLENGTHH
!
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL  :: PTR   ! Transport production of TKE
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL  :: PDISS ! Dissipation of TKE
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL  :: PEDR  ! EDR
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL  :: PLEM  ! Mixing length
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TURB
!
END INTERFACE
!
END MODULE MODI_TURB
