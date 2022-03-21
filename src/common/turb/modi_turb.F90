!     ######spl
     MODULE MODI_TURB  
!    ################ 
!
INTERFACE
!
      SUBROUTINE TURB(CST,CSTURB,BUCONF,TURBN,              &
              & KKA, KKU, KKL, KMI,KRR,KRRL,KRRI,HLBCX,HLBCY,         &
              & KSPLIT,KMODEL_CL,KSV,KSV_LGBEG,KSV_LGEND,             &
              & HPROGRAM, O2D, ONOMIXLG, OFLAT,                       &
              & OLES_CALL,OCOUPLES,OBLOWSNOW,                         &
              & OTURB_FLX,OTURB_DIAG,OSUBG_COND,                      &
              & ORMC01,OOCEAN,ODEEPOC,OHARAT,                         &
              & HTURBDIM,HTURBLEN,HTOM,HTURBLEN_CL,HCLOUD,            &
              & PIMPL,PTSTEP,TPFILE,                                  &
              & PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,                         &
              & PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,PCOSSLOPE,PSINSLOPE,    &
              & PRHODJ,PTHVREF,                                       &
              & PSFTH,PSFRV,PSFSV,PSFU,PSFV,                          &
              & PPABST,PUT,PVT,PWT,PTKET,PSVT,PSRCT,                  &
              & PLENGTHM,PLENGTHH,MFMOIST,                            &
              & PBL_DEPTH, PSBL_DEPTH,                                &
              & PCEI,PCEI_MIN,PCEI_MAX,PCOEF_AMPL_SAT,                &
              & PTHLT,PRT,                                            &
              & PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,PRTKES,              &
              & PSIGS,                                                &
              & PFLXZTHVMF,PWTH,PWRC,PWSV,PDP,PTP,PTDIFF,PTDISS,      &
              & TBUDGETS, KBUDGETS,                                   &
              & PEDR,PLEM,PRTKEMS,PTPMF,                              &
              & PDRUS_TURB,PDRVS_TURB,                                &
              & PDRTHLS_TURB,PDRRTS_TURB,PDRSVS_TURB                  ) 
!
USE MODD_BUDGET, ONLY : TBUDGETDATA,TBUDGETCONF_t
USE MODD_IO, ONLY : TFILEDATA
USE MODD_CST, ONLY: CST_t
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_TURB_n, ONLY: TURB_t
!
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
TYPE(TBUDGETCONF_t),    INTENT(IN)   :: BUCONF
TYPE(TURB_t),           INTENT(IN)   :: TURBN
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=AR
INTEGER,                INTENT(IN)   :: KMI           ! model index number  
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV, KSV_LGBEG, KSV_LGEND ! number of scalar variables
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
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  ODEEPOC      ! activates sfc forcing for ideal ocean deep conv
LOGICAL,                INTENT(IN)   ::  OHARAT
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
LOGICAL,                INTENT(IN)   ::  OLES_CALL    ! compute the LES diagnostics at current time-step
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES version 
LOGICAL,                INTENT(IN)   ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
CHARACTER(LEN=4)       , INTENT(IN)   ::  HTURBDIM  ! dimensionality of the 
                                 ! turbulence scheme
CHARACTER(LEN=4)      , INTENT(IN)   ::  HTURBLEN     ! kind of mixing length
CHARACTER(LEN=4)      , INTENT(IN)   ::  HTOM         ! kind of Third Order Moment
CHARACTER(LEN=4)      , INTENT(IN)   ::  HTURBLEN_CL  ! kind of cloud mixing length
REAL,                   INTENT(IN)   ::  PIMPL        ! degree of implicitness
CHARACTER (LEN=4),      INTENT(IN)   ::  HCLOUD       ! Kind of microphysical scheme
REAL,                   INTENT(IN)   ::  PTSTEP       ! Timestep
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file for MesoNH
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
REAL, DIMENSION(:,:,:),   INTENT(IN),OPTIONAL    ::  PRTKEMS
! Source terms for all passive scalar variables
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) ::  PRSVS
! Sigma_s at time t+1 : square root of the variance of the deviation to the 
! saturation 
REAL, DIMENSION(:,:,:), INTENT(OUT)     ::  PSIGS
REAL, DIMENSION(:,:,:), INTENT(OUT),OPTIONAL     ::  PDRUS_TURB   ! evolution of rhoJ*U   by turbulence only
REAL, DIMENSION(:,:,:), INTENT(OUT),OPTIONAL     ::  PDRVS_TURB   ! evolution of rhoJ*V   by turbulence only
REAL, DIMENSION(:,:,:), INTENT(OUT),OPTIONAL     ::  PDRTHLS_TURB ! evolution of rhoJ*thl by turbulence only
REAL, DIMENSION(:,:,:), INTENT(OUT),OPTIONAL     ::  PDRRTS_TURB  ! evolution of rhoJ*rt  by turbulence only
REAL, DIMENSION(:,:,:,:), INTENT(OUT),OPTIONAL   ::  PDRSVS_TURB  ! evolution of rhoJ*Sv  by turbulence only
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PFLXZTHVMF 
!                                           MF contribution for vert. turb. transport
!                                           used in the buoy. prod. of TKE
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWTH       ! heat flux
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PWRC       ! cloud water flux
REAL, DIMENSION(:,:,:,:),INTENT(OUT) :: PWSV       ! scalar flux
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PDP        ! Dynamic TKE production
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PTP        ! Thermal TKE production
REAL, DIMENSION(:,:,:), INTENT(OUT),OPTIONAL  :: PTPMF      ! Thermal TKE production
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PTDIFF     ! Diffusion TKE term
REAL, DIMENSION(:,:,:), INTENT(OUT)  :: PTDISS     ! Dissipation TKE term
!
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PLENGTHM 
REAL, DIMENSION(:,:,:), INTENT(IN)      ::  PLENGTHH 
!
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! CPROGRAM is the program currently running (modd_conf)
LOGICAL, INTENT(IN) :: ONOMIXLG          ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL, INTENT(IN) :: O2D               ! Logical for 2D model version (modd_conf)
!
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
