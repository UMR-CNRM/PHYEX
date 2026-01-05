MODULE phyex_bridge
    USE ISO_C_BINDING
    ! Import the original routines and required modules
    USE MODI_ICE_ADJUST, ONLY : ICE_ADJUST
    USE MODI_RAIN_ICE, ONLY : RAIN_ICE
    USE MODI_SHALLOW_CONVECTION, ONLY : SHALLOW_CONVECTION
    USE PARKIND1, ONLY : JPIM, JPRB
    USE MODD_DIMPHYEX, ONLY : DIMPHYEX_t
    USE MODD_CST, ONLY : CST_t, CST
    USE MODD_RAIN_ICE_PARAM_n
    USE MODD_RAIN_ICE_DESCR_n
    USE MODD_PARAM_ICE_n
    USE MODD_NEB_n, ONLY : NEB_t
    USE MODD_TURB_n, ONLY : TURB_t
    USE MODD_BUDGET, ONLY : TBUDGETCONF_t, TBUDGETDATA_PTR
    USE MODD_CONVPAR, ONLY : CONVPAR_t
    USE MODD_CONVPAR_SHAL, ONLY : CONVPAR_SHAL
    USE MODD_NSV, ONLY : NSV_t
    USE MODE_INI_CST, ONLY : INI_CST
    USE MODE_INI_RAIN_ICE, ONLY : INI_RAIN_ICE

    IMPLICIT NONE

CONTAINS

    ! C-callable wrapper for ICE_ADJUST
    SUBROUTINE c_ice_adjust_wrap(                                          &
        nlon, nlev, krr, timestep,                                         &
        ptr_sigqsat, ptr_pabs, ptr_sigs, ptr_th, ptr_exn, ptr_exn_ref,    &
        ptr_rho_dry_ref, ptr_rv, ptr_rc, ptr_ri, ptr_rr, ptr_rs, ptr_rg,  &
        ptr_cf_mf, ptr_rc_mf, ptr_ri_mf,                                   &
        ptr_rvs, ptr_rcs, ptr_ris, ptr_ths,                                &
        ptr_cldfr, ptr_icldfr, ptr_wcldfr                                  &
    ) BIND(C, name="c_ice_adjust")
    
        ! C-compatible arguments (using C_FLOAT for single precision)
        INTEGER(C_INT), VALUE, INTENT(IN) :: nlon, nlev, krr
        REAL(C_FLOAT), VALUE, INTENT(IN) :: timestep
        
        ! C pointers for input arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_sigqsat      ! 1D: (nlon)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_pabs         ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_sigs         ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_th           ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_exn          ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_exn_ref      ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rho_dry_ref  ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rv           ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rc           ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ri           ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rr           ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rs           ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rg           ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_cf_mf        ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rc_mf        ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ri_mf        ! 2D: (nlon, nlev)
        
        ! C pointers for input/output tendency arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rvs          ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rcs          ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ris          ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ths          ! 2D: (nlon, nlev)
        
        ! C pointers for output arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_cldfr        ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_icldfr       ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_wcldfr       ! 2D: (nlon, nlev)

        ! Fortran pointers to map C data (using C_FLOAT for single precision)
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:) :: f_sigqsat
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_pabs, f_sigs, f_th, f_exn, f_exn_ref
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_rho_dry_ref, f_rv, f_rc, f_ri
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_rr, f_rs, f_rg
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_cf_mf, f_rc_mf, f_ri_mf
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_rvs, f_rcs, f_ris, f_ths
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_cldfr, f_icldfr, f_wcldfr
        
        ! Local variables for PHYEX structures
        TYPE(DIMPHYEX_t) :: D
        TYPE(RAIN_ICE_PARAM_t) :: ICEP
        TYPE(NEB_t) :: NEBN
        TYPE(TURB_t) :: TURBN
        TYPE(PARAM_ICE_t) :: PARAMI
        TYPE(TBUDGETCONF_t) :: BUCONF
        TYPE(TBUDGETDATA_PTR), DIMENSION(0) :: TBUDGETS
        
        ! Additional required arrays (using C_FLOAT)
        REAL(KIND=C_FLOAT), ALLOCATABLE, DIMENSION(:,:) :: PRHODJ, PZZ
        REAL(KIND=C_FLOAT), ALLOCATABLE, DIMENSION(:,:) :: PMFCONV
        REAL(KIND=C_FLOAT), ALLOCATABLE, DIMENSION(:,:) :: PWEIGHT_MF_CLOUD
        REAL(KIND=C_FLOAT), ALLOCATABLE, DIMENSION(:,:) :: PSSIO, PSSIU, PIFR
        REAL(KIND=C_FLOAT), ALLOCATABLE, DIMENSION(:,:) :: PSRCS
        LOGICAL :: LMFCONV, OCOMPUTE_SRC
        
        ! Convert C pointers to Fortran arrays
        CALL C_F_POINTER(ptr_sigqsat, f_sigqsat, [nlon])
        CALL C_F_POINTER(ptr_pabs, f_pabs, [nlon, nlev])
        CALL C_F_POINTER(ptr_sigs, f_sigs, [nlon, nlev])
        CALL C_F_POINTER(ptr_th, f_th, [nlon, nlev])
        CALL C_F_POINTER(ptr_exn, f_exn, [nlon, nlev])
        CALL C_F_POINTER(ptr_exn_ref, f_exn_ref, [nlon, nlev])
        CALL C_F_POINTER(ptr_rho_dry_ref, f_rho_dry_ref, [nlon, nlev])
        CALL C_F_POINTER(ptr_rv, f_rv, [nlon, nlev])
        CALL C_F_POINTER(ptr_rc, f_rc, [nlon, nlev])
        CALL C_F_POINTER(ptr_ri, f_ri, [nlon, nlev])
        CALL C_F_POINTER(ptr_rr, f_rr, [nlon, nlev])
        CALL C_F_POINTER(ptr_rs, f_rs, [nlon, nlev])
        CALL C_F_POINTER(ptr_rg, f_rg, [nlon, nlev])
        CALL C_F_POINTER(ptr_cf_mf, f_cf_mf, [nlon, nlev])
        CALL C_F_POINTER(ptr_rc_mf, f_rc_mf, [nlon, nlev])
        CALL C_F_POINTER(ptr_ri_mf, f_ri_mf, [nlon, nlev])
        CALL C_F_POINTER(ptr_rvs, f_rvs, [nlon, nlev])
        CALL C_F_POINTER(ptr_rcs, f_rcs, [nlon, nlev])
        CALL C_F_POINTER(ptr_ris, f_ris, [nlon, nlev])
        CALL C_F_POINTER(ptr_ths, f_ths, [nlon, nlev])
        CALL C_F_POINTER(ptr_cldfr, f_cldfr, [nlon, nlev])
        CALL C_F_POINTER(ptr_icldfr, f_icldfr, [nlon, nlev])
        CALL C_F_POINTER(ptr_wcldfr, f_wcldfr, [nlon, nlev])
        
        ! Initialize DIMPHYEX structure
        D%NIT = nlon
        D%NIB = 1
        D%NIE = nlon
        D%NJT = 1
        D%NJB = 1
        D%NJE = 1
        D%NKT = nlev
        D%NKL = 1        ! Ground to space ordering
        D%NKA = 1
        D%NKU = nlev
        D%NKB = 1
        D%NKE = nlev
        D%NKTB = 1
        D%NKTE = nlev
        D%NIBC = 1
        D%NJBC = 1
        D%NIEC = nlon
        D%NJEC = 1
        D%NIJT = nlon
        D%NIJB = 1
        D%NIJE = nlon
        D%NKLES = nlev
        D%NLESMASK = 0
        D%NLES_TIMES = 0
        
        ! Initialize physical constants (uses global CST module)
        CALL INI_CST()
        
        ! Initialize NEBN (nebulosity/cloud parameters) - AROME defaults
        NEBN%LSUBG_COND = .FALSE.     ! No subgrid condensation
        NEBN%LSIGMAS = .TRUE.         ! Use sigma_s
        NEBN%CFRAC_ICE_ADJUST = 'S'   ! Standard
        NEBN%CCONDENS = 'CB02'        ! Condensation scheme
        NEBN%CLAMBDA3 = 'CB'          ! Lambda3 formulation
        
        ! Initialize PARAMI (ice parameters)
        PARAMI%CSUBG_MF_PDF = 'NONE'  ! No mass flux PDF
        PARAMI%LOCND2 = .FALSE.       ! OCND2 option
        
        ! Initialize budget configuration (disable budgets)
        BUCONF%LBUDGET_TH = .FALSE.
        BUCONF%LBUDGET_RV = .FALSE.
        BUCONF%LBUDGET_RC = .FALSE.
        BUCONF%LBUDGET_RI = .FALSE.
        
        ! Allocate and initialize additional required arrays
        ALLOCATE(PRHODJ(nlon, nlev))
        ALLOCATE(PZZ(nlon, nlev))
        ALLOCATE(PMFCONV(nlon, nlev))
        ALLOCATE(PWEIGHT_MF_CLOUD(nlon, nlev))
        ALLOCATE(PSSIO(nlon, nlev))
        ALLOCATE(PSSIU(nlon, nlev))
        ALLOCATE(PIFR(nlon, nlev))
        ALLOCATE(PSRCS(nlon, nlev))
        
        ! Compute PRHODJ from density and assume unit Jacobian
        PRHODJ = f_rho_dry_ref
        
        ! Set height field (simplified - could be passed as parameter)
        PZZ = 0.0_C_FLOAT
        
        ! Initialize mass flux arrays
        PMFCONV = 0.0_C_FLOAT
        PWEIGHT_MF_CLOUD = 0.0_C_FLOAT
        LMFCONV = .FALSE.
        
        ! Initialize output arrays
        PSSIO = 0.0_C_FLOAT
        PSSIU = 0.0_C_FLOAT
        PIFR = 0.0_C_FLOAT
        PSRCS = 0.0_C_FLOAT
        OCOMPUTE_SRC = .FALSE.
        
        ! Call the actual ICE_ADJUST routine (using global CST module)
        CALL ICE_ADJUST(                                                   &
            D, CST, ICEP, NEBN, TURBN, PARAMI, BUCONF, krr,                &
            'BRID',                                                        &
            timestep, f_sigqsat,                                           &
            PRHODJ, f_exn_ref, f_rho_dry_ref, f_sigs, LMFCONV, PMFCONV,   &
            f_pabs, PZZ,                                                   &
            f_exn, f_cf_mf, f_rc_mf, f_ri_mf, PWEIGHT_MF_CLOUD,            &
            f_icldfr, f_wcldfr, PSSIO, PSSIU, PIFR,                        &
            f_rv, f_rc, f_rvs, f_rcs, f_th, f_ths,                         &
            OCOMPUTE_SRC, PSRCS, f_cldfr,                                  &
            f_rr, f_ri, f_ris, f_rs, f_rg, TBUDGETS, 0                     &
        )
        
        ! Cleanup
        DEALLOCATE(PRHODJ, PZZ, PMFCONV, PWEIGHT_MF_CLOUD)
        DEALLOCATE(PSSIO, PSSIU, PIFR, PSRCS)

    END SUBROUTINE c_ice_adjust_wrap

    ! C-callable initialization for RAIN_ICE
    SUBROUTINE c_ini_rain_ice_wrap(timestep, dzmin, krr, hcloud) BIND(C, name="c_ini_rain_ice")
        REAL(C_FLOAT), VALUE, INTENT(IN) :: timestep
        REAL(C_FLOAT), VALUE, INTENT(IN) :: dzmin
        INTEGER(C_INT), VALUE, INTENT(IN) :: krr
        CHARACTER(KIND=C_CHAR), DIMENSION(4), INTENT(IN) :: hcloud

        CHARACTER(LEN=4) :: f_hcloud
        INTEGER :: i, ksplitr
        INTEGER, PARAMETER :: KMODEL = 1  ! Model number (default: 1)

        DO i = 1, 4
            f_hcloud(i:i) = hcloud(i)
        END DO

        ! Initialize physical constants
        CALL INI_CST()

        ! Set up module pointers for model 1
        CALL PARAM_ICE_GOTO_MODEL(0, KMODEL)
        CALL RAIN_ICE_PARAM_GOTO_MODEL(0, KMODEL)
        CALL RAIN_ICE_DESCR_GOTO_MODEL(0, KMODEL)

        ! Initialize RAIN_ICE microphysics scheme
        CALL INI_RAIN_ICE("AROME ", KMODEL, timestep, dzmin, ksplitr, f_hcloud)
    END SUBROUTINE c_ini_rain_ice_wrap

    ! C-callable wrapper for RAIN_ICE
    SUBROUTINE c_rain_ice_wrap(                                            &
        nlon, nlev, krr, timestep,                                         &
        ptr_exn, ptr_dzz, ptr_rhodj, ptr_rhodref, ptr_exnref, ptr_pabs,   &
        ptr_cldfr, ptr_icldfr, ptr_ssio, ptr_ssiu, ptr_ifr,               &
        ptr_tht, ptr_rvt, ptr_rct, ptr_rrt, ptr_rit, ptr_rst, ptr_rgt,    &
        ptr_sigs,                                                          &
        ptr_cit,                                                           &
        ptr_hlc_hrc, ptr_hlc_hcf, ptr_hli_hri, ptr_hli_hcf,               &
        ptr_ths, ptr_rvs, ptr_rcs, ptr_rrs, ptr_ris, ptr_rss, ptr_rgs,    &
        ptr_evap3d, ptr_rainfr,                                            &
        ptr_inprc, ptr_inprr, ptr_inprs, ptr_inprg, ptr_indep,            &
        ptr_rain_ice_param, ptr_rain_ice_descr                             &
    ) BIND(C, name="c_rain_ice")
    
        ! C-compatible arguments (using C_FLOAT for single precision)
        INTEGER(C_INT), VALUE, INTENT(IN) :: nlon, nlev, krr
        REAL(C_FLOAT), VALUE, INTENT(IN) :: timestep
        
        ! C pointers for 2D input arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_exn, ptr_dzz, ptr_rhodj
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rhodref, ptr_exnref, ptr_pabs
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_cldfr, ptr_icldfr
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ssio, ptr_ssiu, ptr_ifr
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_tht, ptr_rvt, ptr_rct
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rrt, ptr_rit, ptr_rst, ptr_rgt
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_sigs
        
        ! C pointers for 2D input/output arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_cit
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_hlc_hrc, ptr_hlc_hcf
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_hli_hri, ptr_hli_hcf
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ths, ptr_rvs, ptr_rcs
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rrs, ptr_ris, ptr_rss, ptr_rgs
        
        ! C pointers for 2D output arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_evap3d, ptr_rainfr
        
        ! C pointers for 1D output arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_inprc, ptr_inprr
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_inprs, ptr_inprg, ptr_indep

        ! C pointers for parameter structures
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rain_ice_param, ptr_rain_ice_descr

        ! Fortran pointers to map C data (using C_FLOAT for single precision)
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_exn, f_dzz, f_rhodj
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_rhodref, f_exnref, f_pabs
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_cldfr, f_icldfr
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_ssio, f_ssiu, f_ifr
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_tht, f_rvt, f_rct
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_rrt, f_rit, f_rst, f_rgt
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_sigs
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_cit
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_hlc_hrc, f_hlc_hcf
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_hli_hri, f_hli_hcf
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_ths, f_rvs, f_rcs
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_rrs, f_ris, f_rss, f_rgs
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_evap3d, f_rainfr
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:) :: f_inprc, f_inprr
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:) :: f_inprs, f_inprg, f_indep

        ! Local variables for PHYEX structures
        TYPE(DIMPHYEX_t) :: D
        TYPE(PARAM_ICE_t) :: PARAMI
        TYPE(TBUDGETCONF_t) :: BUCONF
        TYPE(TBUDGETDATA_PTR), DIMENSION(0) :: TBUDGETS

        ! Fortran pointers for parameter structures (passed from Python)
        TYPE(RAIN_ICE_PARAM_t), POINTER :: rain_ice_param_local
        TYPE(RAIN_ICE_DESCR_t), POINTER :: rain_ice_descr_local

        ! Convert C pointers to Fortran structure pointers
        CALL C_F_POINTER(ptr_rain_ice_param, rain_ice_param_local)
        CALL C_F_POINTER(ptr_rain_ice_descr, rain_ice_descr_local)

        ! Convert C pointers to Fortran arrays
        CALL C_F_POINTER(ptr_exn, f_exn, [nlon, nlev])
        CALL C_F_POINTER(ptr_dzz, f_dzz, [nlon, nlev])
        CALL C_F_POINTER(ptr_rhodj, f_rhodj, [nlon, nlev])
        CALL C_F_POINTER(ptr_rhodref, f_rhodref, [nlon, nlev])
        CALL C_F_POINTER(ptr_exnref, f_exnref, [nlon, nlev])
        CALL C_F_POINTER(ptr_pabs, f_pabs, [nlon, nlev])
        CALL C_F_POINTER(ptr_cldfr, f_cldfr, [nlon, nlev])
        CALL C_F_POINTER(ptr_icldfr, f_icldfr, [nlon, nlev])
        CALL C_F_POINTER(ptr_ssio, f_ssio, [nlon, nlev])
        CALL C_F_POINTER(ptr_ssiu, f_ssiu, [nlon, nlev])
        CALL C_F_POINTER(ptr_ifr, f_ifr, [nlon, nlev])
        CALL C_F_POINTER(ptr_tht, f_tht, [nlon, nlev])
        CALL C_F_POINTER(ptr_rvt, f_rvt, [nlon, nlev])
        CALL C_F_POINTER(ptr_rct, f_rct, [nlon, nlev])
        CALL C_F_POINTER(ptr_rrt, f_rrt, [nlon, nlev])
        CALL C_F_POINTER(ptr_rit, f_rit, [nlon, nlev])
        CALL C_F_POINTER(ptr_rst, f_rst, [nlon, nlev])
        CALL C_F_POINTER(ptr_rgt, f_rgt, [nlon, nlev])
        CALL C_F_POINTER(ptr_sigs, f_sigs, [nlon, nlev])
        CALL C_F_POINTER(ptr_cit, f_cit, [nlon, nlev])
        CALL C_F_POINTER(ptr_hlc_hrc, f_hlc_hrc, [nlon, nlev])
        CALL C_F_POINTER(ptr_hlc_hcf, f_hlc_hcf, [nlon, nlev])
        CALL C_F_POINTER(ptr_hli_hri, f_hli_hri, [nlon, nlev])
        CALL C_F_POINTER(ptr_hli_hcf, f_hli_hcf, [nlon, nlev])
        CALL C_F_POINTER(ptr_ths, f_ths, [nlon, nlev])
        CALL C_F_POINTER(ptr_rvs, f_rvs, [nlon, nlev])
        CALL C_F_POINTER(ptr_rcs, f_rcs, [nlon, nlev])
        CALL C_F_POINTER(ptr_rrs, f_rrs, [nlon, nlev])
        CALL C_F_POINTER(ptr_ris, f_ris, [nlon, nlev])
        CALL C_F_POINTER(ptr_rss, f_rss, [nlon, nlev])
        CALL C_F_POINTER(ptr_rgs, f_rgs, [nlon, nlev])
        CALL C_F_POINTER(ptr_evap3d, f_evap3d, [nlon, nlev])
        CALL C_F_POINTER(ptr_rainfr, f_rainfr, [nlon, nlev])
        CALL C_F_POINTER(ptr_inprc, f_inprc, [nlon])
        CALL C_F_POINTER(ptr_inprr, f_inprr, [nlon])
        CALL C_F_POINTER(ptr_inprs, f_inprs, [nlon])
        CALL C_F_POINTER(ptr_inprg, f_inprg, [nlon])
        CALL C_F_POINTER(ptr_indep, f_indep, [nlon])
        
        ! Initialize DIMPHYEX structure
        D%NIT = nlon
        D%NIB = 1
        D%NIE = nlon
        D%NJT = 1
        D%NJB = 1
        D%NJE = 1
        D%NKT = nlev
        D%NKL = 1
        D%NKA = 1
        D%NKU = nlev
        D%NKB = 1
        D%NKE = nlev
        D%NKTB = 1
        D%NKTE = nlev
        D%NIBC = 1
        D%NJBC = 1
        D%NIEC = nlon
        D%NJEC = 1
        D%NIJT = nlon
        D%NIJB = 1
        D%NIJE = nlon
        D%NKLES = nlev
        D%NLESMASK = 0
        D%NLES_TIMES = 0
        
        ! Initialize physical constants (uses global CST module)
        CALL INI_CST()
        
        ! Initialize PARAMI (microphysics parameters) - same as ICE_ADJUST
        PARAMI%CSUBG_AUCV_RC = 'NONE'
        PARAMI%CSUBG_AUCV_RI = 'NONE'
        PARAMI%CSUBG_PR_PDF = 'SIGM'
        PARAMI%CSUBG_RC_RR_ACCR = 'NONE'
        PARAMI%CSUBG_RR_EVAP = 'NONE'
        PARAMI%LOCND2 = .FALSE.
        PARAMI%LSEDIM_AFTER = .FALSE.
        PARAMI%LWARM = .TRUE.
        PARAMI%LPACK_MICRO = .FALSE.
        PARAMI%NPROMICRO = 0
        PARAMI%LEXCLDROP = .FALSE.
        
        ! Initialize ICEP and ICED (will use default values)
        ! These structures are complex and would need proper initialization
        ! For now, we rely on Fortran's default initialization
        
        ! Initialize budget configuration (disabled)
        BUCONF%LBU_ENABLE = .FALSE.
        BUCONF%LBUDGET_TH = .FALSE.
        BUCONF%LBUDGET_RV = .FALSE.
        BUCONF%LBUDGET_RC = .FALSE.
        BUCONF%LBUDGET_RI = .FALSE.
        BUCONF%LBUDGET_RR = .FALSE.
        BUCONF%LBUDGET_RS = .FALSE.
        BUCONF%LBUDGET_RG = .FALSE.
        BUCONF%LBUDGET_RH = .FALSE.
        
        ! Call the actual RAIN_ICE routine with locally passed structures
        CALL RAIN_ICE(                                                     &
            D, CST, PARAMI, rain_ice_param_local, rain_ice_descr_local, BUCONF, &
            timestep, krr, f_exn,                                          &
            f_dzz, f_rhodj, f_rhodref, f_exnref, f_pabs, f_cit, f_cldfr,   &
            f_icldfr, f_ssio, f_ssiu, f_ifr,                               &
            f_hlc_hrc, f_hlc_hcf, f_hli_hri, f_hli_hcf,                    &
            f_tht, f_rvt, f_rct, f_rrt, f_rit, f_rst,                      &
            f_rgt, f_ths, f_rvs, f_rcs, f_rrs, f_ris, f_rss, f_rgs,        &
            f_inprc, f_inprr, f_evap3d,                                    &
            f_inprs, f_inprg, f_indep, f_rainfr, f_sigs,                   &
            TBUDGETS, 0                                                    &
        )

    END SUBROUTINE c_rain_ice_wrap

    ! C-callable wrapper for SHALLOW_CONVECTION
    SUBROUTINE c_shallow_convection_wrap(                                     &
        nlon, nlev, kice, kbdia, ktdia,                                       &
        osettadj_int, ptadjs, och1conv_int, kch1,                             &
        ptr_ppabst, ptr_pzz, ptr_ptkecls, ptr_ptt, ptr_prvt, ptr_prct,       &
        ptr_prit, ptr_pwt, ptr_ptten, ptr_prvten, ptr_prcten, ptr_priten,    &
        ptr_kcltop, ptr_kclbas, ptr_pumf, ptr_pch1, ptr_pch1ten              &
    ) BIND(C, name="c_shallow_convection")

        ! C-compatible arguments
        INTEGER(C_INT), VALUE, INTENT(IN) :: nlon, nlev, kice, kbdia, ktdia
        INTEGER(C_INT), VALUE, INTENT(IN) :: osettadj_int, och1conv_int, kch1
        REAL(C_FLOAT), VALUE, INTENT(IN) :: ptadjs

        ! C pointers for 1D input arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ptkecls     ! 1D: (nlon)

        ! C pointers for 2D input arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ppabst      ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_pzz         ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ptt         ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_prvt        ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_prct        ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_prit        ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_pwt         ! 2D: (nlon, nlev)

        ! C pointers for 2D input/output arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ptten       ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_prvten      ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_prcten      ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_priten      ! 2D: (nlon, nlev)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_pumf        ! 2D: (nlon, nlev)

        ! C pointers for 1D input/output arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_kcltop      ! 1D: (nlon)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_kclbas      ! 1D: (nlon)

        ! C pointers for 3D chemical tracer arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_pch1        ! 3D: (nlon, nlev, kch1)
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_pch1ten     ! 3D: (nlon, nlev, kch1)

        ! Fortran pointers to map C data
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:) :: f_ptkecls
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_ppabst, f_pzz, f_ptt
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_prvt, f_prct, f_prit, f_pwt
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_ptten, f_prvten, f_prcten, f_priten
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_pumf
        INTEGER(KIND=C_INT), POINTER, DIMENSION(:) :: f_kcltop, f_kclbas
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:,:) :: f_pch1, f_pch1ten

        ! Local variables for PHYEX structures
        TYPE(DIMPHYEX_t) :: D
        TYPE(NSV_t) :: NSV
        TYPE(CONVPAR_t) :: CONVPAR
        TYPE(CONVPAR_SHAL) :: CVP_SHAL
        LOGICAL :: LOSETTADJ, LOCH1CONV

        ! Convert C integers to Fortran logicals
        LOSETTADJ = (osettadj_int /= 0)
        LOCH1CONV = (och1conv_int /= 0)

        ! Convert C pointers to Fortran arrays
        CALL C_F_POINTER(ptr_ptkecls, f_ptkecls, [nlon])
        CALL C_F_POINTER(ptr_ppabst, f_ppabst, [nlon, nlev])
        CALL C_F_POINTER(ptr_pzz, f_pzz, [nlon, nlev])
        CALL C_F_POINTER(ptr_ptt, f_ptt, [nlon, nlev])
        CALL C_F_POINTER(ptr_prvt, f_prvt, [nlon, nlev])
        CALL C_F_POINTER(ptr_prct, f_prct, [nlon, nlev])
        CALL C_F_POINTER(ptr_prit, f_prit, [nlon, nlev])
        CALL C_F_POINTER(ptr_pwt, f_pwt, [nlon, nlev])
        CALL C_F_POINTER(ptr_ptten, f_ptten, [nlon, nlev])
        CALL C_F_POINTER(ptr_prvten, f_prvten, [nlon, nlev])
        CALL C_F_POINTER(ptr_prcten, f_prcten, [nlon, nlev])
        CALL C_F_POINTER(ptr_priten, f_priten, [nlon, nlev])
        CALL C_F_POINTER(ptr_kcltop, f_kcltop, [nlon])
        CALL C_F_POINTER(ptr_kclbas, f_kclbas, [nlon])
        CALL C_F_POINTER(ptr_pumf, f_pumf, [nlon, nlev])
        CALL C_F_POINTER(ptr_pch1, f_pch1, [nlon, nlev, kch1])
        CALL C_F_POINTER(ptr_pch1ten, f_pch1ten, [nlon, nlev, kch1])

        ! Initialize DIMPHYEX structure
        D%NIT = nlon
        D%NIB = 1
        D%NIE = nlon
        D%NJT = 1
        D%NJB = 1
        D%NJE = 1
        D%NKT = nlev
        D%NKL = 1
        D%NKA = 1
        D%NKU = nlev
        D%NKB = 1
        D%NKE = nlev
        D%NKTB = 1
        D%NKTE = nlev
        D%NIBC = 1
        D%NJBC = 1
        D%NIEC = nlon
        D%NJEC = 1
        D%NIJT = nlon
        D%NIJB = 1
        D%NIJE = nlon
        D%NKLES = nlev
        D%NLESMASK = 0
        D%NLES_TIMES = 0

        ! Initialize NSV structure (tracers)
        NSV%NSV_USER = 0
        NSV%NSV_C2R2BEG = 0
        NSV%NSV_C2R2END = 0
        NSV%NSV_C1R3BEG = 0
        NSV%NSV_C1R3END = 0
        NSV%NSV_ELECBEG = 0
        NSV%NSV_ELECEND = 0
        NSV%NSV_LNOXBEG = 0
        NSV%NSV_LNOXEND = 0
        NSV%NSV_DSTBEG = 0
        NSV%NSV_DSTEND = 0
        NSV%NSV_SLTBEG = 0
        NSV%NSV_SLTEND = 0
        NSV%NSV_PPBEG = 0
        NSV%NSV_PPEND = 0
        NSV%NSV_CSBEG = 0
        NSV%NSV_CSEND = 0
        NSV%NSV_AERBEG = 0
        NSV%NSV_AEREND = 0
        NSV%NSV_SNWBEG = 0
        NSV%NSV_SNWEND = 0
        NSV%NSV_CHEMBEG = 0
        NSV%NSV_CHEMEND = 0

        ! Initialize CONVPAR structure (deep convection parameters)
        CONVPAR%XA25 = 625.0E6_C_FLOAT     ! Reference grid area (25km)^2
        CONVPAR%XCRAD = 1500.0_C_FLOAT     ! Cloud radius (m)
        CONVPAR%XCDEPTH = 3000.0_C_FLOAT   ! Minimum necessary cloud depth
        CONVPAR%XENTR = 0.03_C_FLOAT       ! Entrainment constant
        CONVPAR%XZLCL = 3500.0_C_FLOAT     ! Max LCL height
        CONVPAR%XZPBL = 6000.0_C_FLOAT     ! Minimum PBL height
        CONVPAR%XWTRIG = 6.0_C_FLOAT       ! Trigger vertical velocity
        CONVPAR%XNHGAM = 1.3333_C_FLOAT    ! Non-hydrostatic pressure factor
        CONVPAR%XTFRZ1 = 268.16_C_FLOAT    ! Freezing interval begin
        CONVPAR%XTFRZ2 = 248.16_C_FLOAT    ! Freezing interval end
        CONVPAR%XRHDBC = 0.9_C_FLOAT       ! Relative humidity below cloud
        CONVPAR%XRCONV = 0.015_C_FLOAT     ! Precipitation conversion constant
        CONVPAR%XSTABT = 0.75_C_FLOAT      ! Stability in fractional time integration
        CONVPAR%XSTABC = 0.95_C_FLOAT      ! Stability in CAPE adjustment
        CONVPAR%XUSRDPTH = 16500.0_C_FLOAT ! Pressure thickness for updraft moisture
        CONVPAR%XMELDPTH = 10000.0_C_FLOAT ! Layer for precipitation melt
        CONVPAR%XUVDP = 0.7_C_FLOAT        ! Pressure perturbation in momentum transport

        ! Initialize CVP_SHAL structure (shallow convection parameters)
        CVP_SHAL%XA25 = 625.0E6_C_FLOAT       ! Reference grid area
        CVP_SHAL%XCRAD = 1500.0_C_FLOAT       ! Cloud radius
        CVP_SHAL%XCTIME_SHAL = 10800.0_C_FLOAT ! Convective adjustment time
        CVP_SHAL%XCDEPTH = 2500.0_C_FLOAT     ! Minimum cloud depth
        CVP_SHAL%XCDEPTH_D = 3000.0_C_FLOAT   ! Maximum cloud thickness
        CVP_SHAL%XDTPERT = 1.0_C_FLOAT        ! Temperature perturbation at LCL
        CVP_SHAL%XATPERT = 0.0_C_FLOAT        ! Parameter for temp perturbation
        CVP_SHAL%XBTPERT = 0.0_C_FLOAT        ! Parameter for temp perturbation
        CVP_SHAL%XENTR = 0.03_C_FLOAT         ! Entrainment constant
        CVP_SHAL%XZLCL = 3500.0_C_FLOAT       ! Max LCL height
        CVP_SHAL%XZPBL = 6000.0_C_FLOAT       ! Minimum PBL height
        CVP_SHAL%XWTRIG = 6.0_C_FLOAT         ! Trigger vertical velocity
        CVP_SHAL%XNHGAM = 1.3333_C_FLOAT      ! Non-hydrostatic pressure factor
        CVP_SHAL%XTFRZ1 = 268.16_C_FLOAT      ! Freezing interval begin
        CVP_SHAL%XTFRZ2 = 248.16_C_FLOAT      ! Freezing interval end
        CVP_SHAL%XSTABT = 0.75_C_FLOAT        ! Stability factor
        CVP_SHAL%XSTABC = 0.95_C_FLOAT        ! Stability in CAPE adjustment
        CVP_SHAL%XAW = 1.0_C_FLOAT            ! WLCL parameter A
        CVP_SHAL%XBW = 0.0_C_FLOAT            ! WLCL parameter B
        CVP_SHAL%LLSMOOTH = .TRUE.            ! Smoothing flag

        ! Initialize physical constants
        CALL INI_CST()

        ! Call the actual SHALLOW_CONVECTION routine
        CALL SHALLOW_CONVECTION(                                               &
            CVP_SHAL, CST, D, NSV, CONVPAR, kbdia, ktdia,                      &
            kice, LOSETTADJ, ptadjs, f_ppabst, f_pzz,                          &
            f_ptkecls, f_ptt, f_prvt, f_prct, f_prit, f_pwt,                   &
            f_ptten, f_prvten, f_prcten, f_priten,                             &
            f_kcltop, f_kclbas, f_pumf, LOCH1CONV, kch1,                       &
            f_pch1, f_pch1ten                                                  &
        )

    END SUBROUTINE c_shallow_convection_wrap

END MODULE phyex_bridge