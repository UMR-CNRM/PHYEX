!     ######spl
     MODULE MODI_TURB  
!    ################ 
!
INTERFACE
!
      SUBROUTINE TURB(KKA, KKU, KKL, KMI,KRR,KRRL,KRRI,HLBCX,HLBCY,   &
              & KSPLIT,KMODEL_CL, &
              & OCLOSE_OUT,OTURB_FLX,OTURB_DIAG,OSUBG_COND,ORMC01,    &
              & HTURBDIM,HTURBLEN,HTOM,HTURBLEN_CL,HINST_SFU,         &
              & HMF_UPDRAFT,PIMPL,PTSTEP_UVW, PTSTEP_MET,PTSTEP_SV,   &
              & HFMFILE,HLUOUT,PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,          &
              & PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,PCOSSLOPE,PSINSLOPE,    &
              & PRHODJ,PTHVREF,PRHODREF,                              &
              & PSFTH,PSFRV,PSFSV,PSFU,PSFV,                          &
              & PPABST,PUT,PVT,PWT,PTKET,PSVT,PSRCT,                  &
              & PLENGTHM,PLENGTHH,MFMOIST,                            &
              & PBL_DEPTH, PSBL_DEPTH,                                &
              & PCEI,PCEI_MIN,PCEI_MAX,PCOEF_AMPL_SAT,    &
              & PTHLT,PRT,                                            &
              & PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,PRTKES,              &
              & PHGRAD,PSIGS,                                         &
              & PDRUS_TURB,PDRVS_TURB,                                &
              & PDRTHLS_TURB,PDRRTS_TURB,PDRSVS_TURB,                 &
              & PFLXZTHVMF,PWTH,PWRC,PWSV,PDP,PTP,PTPMF,PTDIFF,PTDISS,&
              & YDDDH,YDLDDH,YDMDDH,                                  &
              & TBUDGETS, KBUDGETS,                                   &
              & PTR,PDISS,PEDR,PLEM                                   ) 
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
USE MODD_BUDGET, ONLY : TBUDGETDATA
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=AR
INTEGER,                INTENT(IN)   :: KMI           ! model index number  
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
CHARACTER(LEN=*),DIMENSION(:),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
CHARACTER(LEN=4),INTENT(IN)          :: HMF_UPDRAFT   ! Type of mass flux
INTEGER,                INTENT(IN)   :: KSPLIT        ! number of time-splitting
INTEGER,                INTENT(IN)   :: KMODEL_CL     ! model number for cloud mixing length
LOGICAL,                INTENT(IN)   ::  OCLOSE_OUT   ! switch for syncronous
                                                      ! file opening
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                 ! diagnostic fields in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OSUBG_COND   ! switch for SUBGrid 
                                 ! CONDensation
LOGICAL,                INTENT(IN)   ::  ORMC01       ! switch for RMC01 lengths in SBL
CHARACTER*4           , INTENT(IN)      ::  HTURBDIM  ! dimensionality of the 
                                 ! turbulence scheme
CHARACTER*4           , INTENT(IN)   ::  HTURBLEN     ! kind of mixing length
CHARACTER*4           , INTENT(IN)   ::  HTOM         ! kind of Third Order Moment
CHARACTER*4           , INTENT(IN)   ::  HTURBLEN_CL  ! kind of cloud mixing length
CHARACTER*1           , INTENT(IN)   ::  HINST_SFU    ! temporal location of the
                                                      ! surface friction flux
REAL,                   INTENT(IN)   ::  PIMPL        ! degree of implicitness
REAL,                   INTENT(IN)   ::  PTSTEP_UVW   ! Dynamical timestep 
REAL,                   INTENT(IN)   ::  PTSTEP_MET   ! Timestep for meteorological variables                        
REAL,                   INTENT(IN)   ::  PTSTEP_SV    ! Timestep for tracer variables
CHARACTER(LEN=*),       INTENT(IN)   ::  HFMFILE      ! Name of the output
                                                      ! FM-file
CHARACTER(LEN=*),       INTENT(IN)   ::  HLUOUT       ! Output-listing name for
                                                      ! model n
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
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  MFMOIST   ! Moist mass flux DUal scheme

REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PTHVREF   ! Virtual Potential
                                        ! Temperature of the reference state
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PRHODREF  ! dry density of the 
                                        ! reference state
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
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PBL_DEPTH  ! BL depth for TOMS
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PSBL_DEPTH ! SBL depth for RMC01
!
!    variables for cloud mixing length
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PCEI ! Cloud Entrainment instability
                                                 ! index to emphasize localy 
                                                 ! turbulent fluxes
REAL, INTENT(IN)      ::  PCEI_MIN ! minimum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCEI_MAX ! maximum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCOEF_AMPL_SAT ! saturation of the amplification coefficient
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PTHLT       ! conservative pot. temp.
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRT         ! water var.  where 
                             ! PRM(:,:,:,1) is the conservative mixing ratio        
!
! sources of momentum, conservative potential temperature, Turb. Kin. Energy, 
! TKE dissipation
REAL, DIMENSION(:,:,:),   INTENT(INOUT) ::  PRUS,PRVS,PRWS,PRTHLS,PRTKES
! Source terms for all water kinds, PRRS(:,:,:,1) is used for the conservative
! mixing ratio
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRRS 
! Source terms for all passive scalar variables
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS
! Sigma_s at time t+1 : square root of the variance of the deviation to the 
! saturation 
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PHGRAD
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PSIGS
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PDRUS_TURB   ! evolution of rhoJ*U   by turbulence only
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PDRVS_TURB   ! evolution of rhoJ*V   by turbulence only
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PDRTHLS_TURB ! evolution of rhoJ*thl by turbulence only
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PDRRTS_TURB  ! evolution of rhoJ*rt  by turbulence only
REAL, DIMENSION(:,:,:,:), INTENT(OUT)   ::  PDRSVS_TURB  ! evolution of rhoJ*Sv  by turbulence only
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PFLXZTHVMF 
!                                           MF contribution for vert. turb. transport
!                                           used in the buoy. prod. of TKE
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWTH       ! heat flux
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWRC       ! cloud water flux
REAL, DIMENSION(:,:,:,:),INTENT(OUT) :: PWSV       ! scalar flux
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PDP        ! Dynamic TKE production
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PTP        ! Thermal TKE production
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PTPMF      ! Thermal TKE production
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PTDIFF     ! Diffusion TKE term
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PTDISS     ! Dissipation TKE term
!
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PLENGTHM 
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PLENGTHH 
!
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH),   INTENT(IN)    :: YDLDDH
TYPE(TMDDH),   INTENT(IN)    :: YDMDDH
!
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
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
