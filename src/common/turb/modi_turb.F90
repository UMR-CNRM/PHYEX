!     ######spl
     MODULE MODI_TURB  
!    ################ 
!
IMPLICIT NONE
INTERFACE
!
      SUBROUTINE TURB(CST,CSTURB,BUCONF,TURBN,NEBN,D,TLES,            &
              & KRR,KRRL,KRRI,HLBCX,HLBCY,KGRADIENTS,KHALO,           &
              & KSPLIT,OCLOUDMODIFLM,KSV,KSV_LGBEG,KSV_LGEND,         &
              & KSV_LIMA_NR, KSV_LIMA_NS, KSV_LIMA_NG, KSV_LIMA_NH,   &
              & O2D,ONOMIXLG,OFLAT,OCOUPLES,OBLOWSNOW,OIBM,OFLYER,    &
              & OCOMPUTE_SRC, PRSNOW,                                 &
              & OOCEAN,ODEEPOC,ODIAG_IN_RUN,                          &
              & HTURBLEN_CL,HCLOUD,HELEC,                             &
              & PTSTEP,TPFILE,                                        &
              & PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,                         &
              & PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,PCOSSLOPE,PSINSLOPE,    &
              & PRHODJ,PTHVREF,PHGRAD,PZS,                            &
              & PSFTH,PSFRV,PSFSV,PSFU,PSFV,                          &
              & PPABST,PUT,PVT,PWT,PTKET,PSVT,PSRCT,                  &
              & PLENGTHM,PLENGTHH,MFMOIST,                            &
              & PBL_DEPTH,PSBL_DEPTH,                                 &
              & PCEI,PCEI_MIN,PCEI_MAX,PCOEF_AMPL_SAT,                &
              & PTHLT,PRT,                                            &
              & PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,PRTKES,              &
              & PSIGS,                                                &
              & PFLXZTHVMF,PWTH,PWRC,PWSV,PDP,PTP,PTDIFF,PTDISS,      &
              & TBUDGETS, KBUDGETS,                                   &
              & PEDR,PLEM,PRTKEMS,PTPMF,                              &
              & PDRUS_TURB,PDRVS_TURB,                                &
              & PDRTHLS_TURB,PDRRTS_TURB,PDRSVS_TURB,PTR,PDISS,       &
              & PIBM_LS, PIBM_XMUT,                                   &
              & PCURRENT_TKE_DISS, PSSTFL, PSSTFL_C, PSSRFL_C,        &
              & PSSUFL_C, PSSVFL_C,PSSUFL,PSSVFL                      )
!
USE MODD_BUDGET, ONLY : TBUDGETDATA,TBUDGETCONF_t
USE MODD_IO, ONLY : TFILEDATA
USE MODD_CST, ONLY: CST_t
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_TURB_n, ONLY: TURB_t
USE MODD_NEB_n, ONLY: NEB_t
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_LES, ONLY: TLES_t
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D             ! PHYEX variables dimensions structure
TYPE(CST_t),            INTENT(IN)   :: CST           ! modd_cst general constant structure
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB        ! modd_csturb turb constant structure
TYPE(TBUDGETCONF_t),    INTENT(IN)   :: BUCONF        ! budget structure
TYPE(TURB_t),           INTENT(IN)   :: TURBN         ! modn_turbn (turb namelist) structure
TYPE(NEB_t),            INTENT(IN)   :: NEBN          ! modd_nebn structure
TYPE(TLES_t),           INTENT(INOUT)   :: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KGRADIENTS    ! Number of stored horizontal gradients
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV, KSV_LGBEG, KSV_LGEND ! number of scalar variables
INTEGER,                INTENT(IN)   :: KSV_LIMA_NR,KSV_LIMA_NS,KSV_LIMA_NG,KSV_LIMA_NH
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
INTEGER,                INTENT(IN)   :: KSPLIT        ! number of time-splitting
LOGICAL,                INTENT(IN)   :: OCLOUDMODIFLM ! cloud mixing length modifs
INTEGER,                INTENT(IN)   ::  KHALO        ! Size of the halo for parallel distribution
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  ODEEPOC      ! activates sfc forcing for ideal ocean deep conv
LOGICAL,                INTENT(IN)   ::  OFLYER       ! MesoNH flyer diagnostic
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES version 
LOGICAL,                INTENT(IN)   ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
LOGICAL,                INTENT(IN)   ::  ODIAG_IN_RUN ! switch to activate online diagnostics (mesonh)
LOGICAL,                INTENT(IN)   ::  OIBM         ! switch to modity mixing length near building with IBM
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBLEN_CL  ! kind of cloud mixing length
CHARACTER (LEN=4),      INTENT(IN)   ::  HCLOUD       ! Kind of microphysical scheme
CHARACTER (LEN=4),      INTENT(IN)   ::  HELEC        ! Kind of cloud electricity scheme
REAL,                   INTENT(IN)   ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
REAL,                   INTENT(IN)   ::  PTSTEP       ! timestep
TYPE(TFILEDATA),        INTENT(IN)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PDXX,PDYY,PDZZ,PDZX,PDZY
                                        ! metric coefficients
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PZZ       !  physical distance
! between 2 succesive grid points along the K direction
REAL, DIMENSION(D%NIJT),   INTENT(IN)      ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PZS ! orography (for LEONARD terms)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  PRHODJ    ! dry density * Grid size
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  MFMOIST ! moist mass flux dual scheme
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  PTHVREF   ! Virtual Potential
                                        ! Temperature of the reference state
REAL, DIMENSION(D%NIJT,D%NKT,KGRADIENTS),   INTENT(IN) ::  PHGRAD      ! horizontal gradients
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)      ::  PSFTH,PSFRV,   &
! normal surface fluxes of theta and Rv
                                            PSFU,PSFV
! normal surface fluxes of (u,v) parallel to the orography
REAL, DIMENSION(D%NIJT,KSV), INTENT(IN)      ::  PSFSV
! normal surface fluxes of Scalar var.
!
!    prognostic variables at t- deltat
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PPABST      ! Pressure at time t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PUT,PVT,PWT ! wind components
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PTKET       ! TKE
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVT        ! passive scal. var.
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)),   INTENT(IN) ::  PSRCT       ! Second-order flux
                      ! s'rc'/2Sigma_s2 at time t-1 multiplied by Lambda_3
REAL, DIMENSION(MERGE(D%NIJT,0,TURBN%CTOM=='TM06')),INTENT(INOUT) :: PBL_DEPTH  ! BL height for TOMS
REAL, DIMENSION(MERGE(D%NIJT,0,TURBN%LRMC01)),INTENT(INOUT) :: PSBL_DEPTH ! SBL depth for RMC01
!
!    variables for cloud mixing length
REAL, DIMENSION(MERGE(D%NIJT,0,OCLOUDMODIFLM),&
                MERGE(D%NKT,0,OCLOUDMODIFLM)),INTENT(IN)      ::  PCEI
                                                 ! Cloud Entrainment instability
                                                 ! index to emphasize localy
                                                 ! turbulent fluxes
REAL, INTENT(IN)      ::  PCEI_MIN ! minimum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCEI_MAX ! maximum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCOEF_AMPL_SAT ! saturation of the amplification coefficient
!
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PTHLT       ! conservative pot. temp.
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT) ::  PRT         ! water var.  where
                             ! PRT(:,:,:,1) is the conservative mixing ratio
!
! sources of momentum, conservative potential temperature, Turb. Kin. Energy,
! TKE dissipation
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PRUS,PRVS,PRWS,PRTHLS,PRTKES
! Source terms for all water kinds, PRRS(:,:,:,1) is used for the conservative
! mixing ratio
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN),OPTIONAL    ::  PRTKEMS
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT) ::  PRRS
! Source terms for all passive scalar variables
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(INOUT) ::  PRSVS
! Sigma_s at time t+1 : square root of the variance of the deviation to the
! saturation
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(OUT)     ::  PSIGS
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRUS_TURB   ! evolution of rhoJ*U   by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRVS_TURB   ! evolution of rhoJ*V   by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRTHLS_TURB ! evolution of rhoJ*thl by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRRTS_TURB  ! evolution of rhoJ*rt  by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT),OPTIONAL ::  PDRSVS_TURB  ! evolution of rhoJ*Sv  by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  PFLXZTHVMF
!                                           MF contribution for vert. turb. transport
!                                           used in the buoy. prod. of TKE
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PWTH       ! heat flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PWRC       ! cloud water flux
REAL, DIMENSION(D%NIJT,D%NKT,KSV),INTENT(OUT) :: PWSV       ! scalar flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTP        ! Thermal TKE production
                                                   ! MassFlux + turb
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL  :: PTPMF      ! Thermal TKE production
                                                   ! MassFlux Only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PDP        ! Dynamic TKE production
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTDIFF     ! Diffusion TKE term
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTDISS     ! Dissipation TKE term
!
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
LOGICAL, INTENT(IN) :: ONOMIXLG          ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL, INTENT(IN) :: O2D               ! Logical for 2D model version (modd_conf)
!
! length scale from vdfexcu
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PLENGTHM, PLENGTHH
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT), OPTIONAL  :: PEDR  ! EDR
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT), OPTIONAL  :: PLEM  ! Mixing length
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT), OPTIONAL  ::  PTR          ! Transport prod. of TKE
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT), OPTIONAL  ::  PDISS        ! Dissipation of TKE
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT), OPTIONAL  ::  PCURRENT_TKE_DISS ! if ODIAG_IN_RUN in mesonh
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSTFL        ! Time evol Flux of T at sea surface (LOCEAN)
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSTFL_C  ! O-A interface flux for theta(LOCEAN and LCOUPLES)
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSRFL_C  ! O-A interface flux for vapor (LOCEAN and LCOUPLES) 
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSUFL_C        ! Time evol Flux of U at sea surface (LOCEAN)
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSVFL_C  !
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSUFL   
REAL, DIMENSION(D%NIJT), INTENT(IN),OPTIONAL   ::  PSSVFL  !
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT), OPTIONAL :: PIBM_XMUT ! IBM turbulent viscosity
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN), OPTIONAL  :: PIBM_LS ! IBM Level-set function
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TURB
!
END INTERFACE
!
END MODULE MODI_TURB
