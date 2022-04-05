!     ######spl
     MODULE MODI_TURB  
!    ################ 
!
INTERFACE
!
      SUBROUTINE TURB(CST,CSTURB,BUCONF,TURBN,D,              &
              & KMI,KRR,KRRL,KRRI,HLBCX,HLBCY,         &
              & KSPLIT,KMODEL_CL,KSV,KSV_LGBEG,KSV_LGEND,             &
              & HPROGRAM, O2D, ONOMIXLG, OFLAT,                       &
              & OLES_CALL,OCOUPLES,OBLOWSNOW,                         &
              & OTURB_FLX,OTURB_DIAG,OSUBG_COND,OCOMPUTE_SRC,         &
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
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
TYPE(TBUDGETCONF_t),    INTENT(IN)   :: BUCONF
TYPE(TURB_t),           INTENT(IN)   :: TURBN
INTEGER,                INTENT(IN)   :: KMI           ! model index number  
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV, KSV_LGBEG, KSV_LGEND ! number of scalar variables
CHARACTER(LEN=*),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
INTEGER,                INTENT(IN)   :: KSPLIT        ! number of time-splitting
INTEGER,                INTENT(IN)   :: KMODEL_CL     ! model number for cloud mixing length
LOGICAL,                INTENT(IN)   ::  OTURB_FLX    ! switch to write the
                                 ! turbulent fluxes in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                 ! diagnostic fields in the syncronous FM-file
LOGICAL,                INTENT(IN)   ::  OSUBG_COND   ! switch for SUBGrid CONDensation
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables
LOGICAL,                INTENT(IN)   ::  ORMC01       ! switch for RMC01 lengths in SBL
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  ODEEPOC      ! activates sfc forcing for ideal ocean deep conv
LOGICAL,                INTENT(IN)   ::  OHARAT       ! switch for LHARATU from AROME
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
LOGICAL,                INTENT(IN)   ::  OLES_CALL    ! compute the LES diagnostics at current time-step
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES version 
LOGICAL,                INTENT(IN)   ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   :: PDXX,PDYY,PDZZ,PDZX,PDZY
                                        ! metric coefficients
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   :: PZZ       !  physical distance 
! between 2 succesive grid points along the K direction
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)      ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)      ::  PRHODJ    ! dry density * Grid size
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)      ::  MFMOIST ! moist mass flux dual scheme
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)      ::  PTHVREF   ! Virtual Potential
                                        ! Temperature of the reference state
!
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)      ::  PSFTH,PSFRV,   &
! normal surface fluxes of theta and Rv 
                                            PSFU,PSFV
! normal surface fluxes of (u,v) parallel to the orography
REAL, DIMENSION(D%NIT,D%NJT,KSV), INTENT(IN)      ::  PSFSV
! normal surface fluxes of Scalar var. 
!
!    prognostic variables at t- deltat
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) ::  PPABST      ! Pressure at time t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) ::  PUT,PVT,PWT ! wind components
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) ::  PTKET       ! TKE
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(IN) ::  PSVT        ! passive scal. var.
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN) ::  PSRCT       ! Second-order flux
                      ! s'rc'/2Sigma_s2 at time t-1 multiplied by Lambda_3
REAL, DIMENSION(MERGE(D%NIT,0,HTOM=='TMO6'),&
                MERGE(D%NJT,0,HTOM=='TMO6')),INTENT(INOUT) :: PBL_DEPTH  ! BL height for TOMS
REAL, DIMENSION(MERGE(D%NIT,0,ORMC01),&
                MERGE(D%NJT,0,ORMC01)),INTENT(INOUT) :: PSBL_DEPTH ! SBL depth for RMC01
!
!    variables for cloud mixing length
REAL, DIMENSION(MERGE(D%NIT,0,KMODEL_CL==KMI .AND. HTURBLEN_CL/='NONE'),&
                MERGE(D%NJT,0,KMODEL_CL==KMI .AND. HTURBLEN_CL/='NONE'),&
                MERGE(D%NJT,0,KMODEL_CL==KMI .AND. HTURBLEN_CL/='NONE')),INTENT(IN)      ::  PCEI 
                                                 ! Cloud Entrainment instability
                                                 ! index to emphasize localy 
                                                 ! turbulent fluxes
REAL, INTENT(IN)      ::  PCEI_MIN ! minimum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCEI_MAX ! maximum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCOEF_AMPL_SAT ! saturation of the amplification coefficient
!
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PTHLT       ! conservative pot. temp.
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(INOUT) ::  PRT         ! water var.  where 
                             ! PRT(:,:,:,1) is the conservative mixing ratio        
!
! sources of momentum, conservative potential temperature, Turb. Kin. Energy, 
! TKE dissipation
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(INOUT) ::  PRUS,PRVS,PRWS,PRTHLS,PRTKES
! Source terms for all water kinds, PRRS(:,:,:,1) is used for the conservative
! mixing ratio
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),   INTENT(IN),OPTIONAL    ::  PRTKEMS
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(INOUT) ::  PRRS
! Source terms for all passive scalar variables
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(INOUT) ::  PRSVS
! Sigma_s at time t+1 : square root of the variance of the deviation to the 
! saturation
REAL, DIMENSION(MERGE(D%NIT,0,OCOMPUTE_SRC),&
                MERGE(D%NIT,0,OCOMPUTE_SRC),&
                MERGE(D%NIT,0,OCOMPUTE_SRC)), INTENT(OUT)     ::  PSIGS
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRUS_TURB   ! evolution of rhoJ*U   by turbulence only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRVS_TURB   ! evolution of rhoJ*V   by turbulence only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRTHLS_TURB ! evolution of rhoJ*thl by turbulence only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRRTS_TURB  ! evolution of rhoJ*rt  by turbulence only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(OUT),OPTIONAL ::  PDRSVS_TURB  ! evolution of rhoJ*Sv  by turbulence only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)      ::  PFLXZTHVMF 
!                                           MF contribution for vert. turb. transport
!                                           used in the buoy. prod. of TKE
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PWTH       ! heat flux
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PWRC       ! cloud water flux
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV),INTENT(OUT) :: PWSV       ! scalar flux
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PTP        ! Thermal TKE production
                                                   ! MassFlux + turb
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT),OPTIONAL  :: PTPMF      ! Thermal TKE production
                                                   ! MassFlux Only
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PDP        ! Dynamic TKE production
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PTDIFF     ! Diffusion TKE term
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: PTDISS     ! Dissipation TKE term
!
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! CPROGRAM is the program currently running (modd_conf)
LOGICAL, INTENT(IN) :: ONOMIXLG          ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL, INTENT(IN) :: O2D               ! Logical for 2D model version (modd_conf)
!
! length scale from vdfexcu
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PLENGTHM, PLENGTHH
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT), OPTIONAL  :: PEDR  ! EDR
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT), OPTIONAL  :: PLEM  ! Mixing length
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TURB
!
END INTERFACE
!
END MODULE MODI_TURB
