MODULE phyex_bridge_acc
    USE ISO_C_BINDING
    ! Import the OpenACC GPU-accelerated routines
    USE MODI_ICE_ADJUST
    USE MODI_RAIN_ICE
    USE PARKIND1, ONLY : JPIM, JPRB
    USE MODD_DIMPHYEX, ONLY : DIMPHYEX_t
    USE MODD_CST, ONLY : CST_t, CST
    USE MODD_RAIN_ICE_PARAM_n
    USE MODD_RAIN_ICE_DESCR_n
    USE MODD_PARAM_ICE_n
    USE MODD_NEB_n, ONLY : NEB_t
    USE MODD_TURB_n, ONLY : TURB_t
    USE MODD_BUDGET, ONLY : TBUDGETCONF_t, TBUDGETDATA_PTR
    USE MODE_INI_CST, ONLY : INI_CST
    USE MODE_INI_RAIN_ICE, ONLY : INI_RAIN_ICE

    IMPLICIT NONE

CONTAINS

    !===========================================================================
    ! C-callable GPU wrapper for ICE_ADJUST_ACC
    !
    ! This subroutine expects ALL input arrays to be GPU device pointers.
    ! Use CuPy or other GPU array library to allocate and manage GPU memory.
    !
    ! Example Python/CuPy usage:
    !   import cupy as cp
    !   th_gpu = cp.asarray(th_cpu, dtype=np.float32)
    !   # Pass th_gpu.data.ptr to this function
    !===========================================================================
    SUBROUTINE c_ice_adjust_acc_wrap(                                          &
        nlon, nlev, krr, timestep,                                             &
        ptr_sigqsat, ptr_pabs, ptr_sigs, ptr_th, ptr_exn, ptr_exn_ref,        &
        ptr_rho_dry_ref, ptr_rv, ptr_rc, ptr_ri, ptr_rr, ptr_rs, ptr_rg,      &
        ptr_cf_mf, ptr_rc_mf, ptr_ri_mf,                                       &
        ptr_rvs, ptr_rcs, ptr_ris, ptr_ths,                                    &
        ptr_cldfr, ptr_icldfr, ptr_wcldfr                                      &
    ) BIND(C, name="c_ice_adjust_acc")

        !-----------------------------------------------------------------------
        ! Arguments (C-compatible)
        !-----------------------------------------------------------------------
        INTEGER(C_INT), VALUE, INTENT(IN) :: nlon, nlev, krr
        REAL(C_FLOAT), VALUE, INTENT(IN) :: timestep

        ! GPU device pointers for input arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_sigqsat      ! 1D: (nlon) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_pabs         ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_sigs         ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_th           ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_exn          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_exn_ref      ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rho_dry_ref  ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rv           ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rc           ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ri           ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rr           ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rs           ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rg           ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_cf_mf        ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rc_mf        ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ri_mf        ! 2D: (nlon, nlev) - GPU

        ! GPU device pointers for input/output tendency arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rvs          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rcs          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ris          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ths          ! 2D: (nlon, nlev) - GPU

        ! GPU device pointers for output arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_cldfr        ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_icldfr       ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_wcldfr       ! 2D: (nlon, nlev) - GPU

        !-----------------------------------------------------------------------
        ! Fortran pointers to GPU data (using C_FLOAT for single precision)
        !-----------------------------------------------------------------------
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:) :: f_sigqsat
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_pabs, f_sigs, f_th, f_exn, f_exn_ref
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_rho_dry_ref, f_rv, f_rc, f_ri
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_rr, f_rs, f_rg
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_cf_mf, f_rc_mf, f_ri_mf
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_rvs, f_rcs, f_ris, f_ths
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_cldfr, f_icldfr, f_wcldfr

        !-----------------------------------------------------------------------
        ! Local variables for PHYEX structures
        !-----------------------------------------------------------------------
        TYPE(DIMPHYEX_t) :: D
        TYPE(RAIN_ICE_PARAM_t) :: ICEP
        TYPE(NEB_t) :: NEBN
        TYPE(TURB_t) :: TURBN
        TYPE(PARAM_ICE_t) :: PARAMI
        TYPE(TBUDGETCONF_t) :: BUCONF
        TYPE(TBUDGETDATA_PTR), DIMENSION(0) :: TBUDGETS

        !-----------------------------------------------------------------------
        ! Additional required arrays (allocated on GPU)
        !-----------------------------------------------------------------------
        REAL(KIND=C_FLOAT), ALLOCATABLE, DIMENSION(:,:) :: PRHODJ, PZZ
        REAL(KIND=C_FLOAT), ALLOCATABLE, DIMENSION(:,:) :: PMFCONV
        REAL(KIND=C_FLOAT), ALLOCATABLE, DIMENSION(:,:) :: PWEIGHT_MF_CLOUD
        REAL(KIND=C_FLOAT), ALLOCATABLE, DIMENSION(:,:) :: PSSIO, PSSIU, PIFR
        REAL(KIND=C_FLOAT), ALLOCATABLE, DIMENSION(:,:) :: PSRCS
        LOGICAL :: LMFCONV, OCOMPUTE_SRC

        !-----------------------------------------------------------------------
        ! Convert C pointers to Fortran pointers
        !-----------------------------------------------------------------------
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

        !-----------------------------------------------------------------------
        ! Initialize DIMPHYEX structure
        !-----------------------------------------------------------------------
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

        !-----------------------------------------------------------------------
        ! Initialize physical constants (uses global CST module)
        !-----------------------------------------------------------------------
        CALL INI_CST()

        !-----------------------------------------------------------------------
        ! Initialize NEBN (nebulosity/cloud parameters) - AROME defaults
        !-----------------------------------------------------------------------
        NEBN%LSUBG_COND = .FALSE.     ! No subgrid condensation
        NEBN%LSIGMAS = .TRUE.         ! Use sigma_s
        NEBN%CFRAC_ICE_ADJUST = 'S'   ! Standard
        NEBN%CCONDENS = 'CB02'        ! Condensation scheme
        NEBN%CLAMBDA3 = 'CB'          ! Lambda3 formulation

        !-----------------------------------------------------------------------
        ! Initialize PARAMI (ice parameters)
        !-----------------------------------------------------------------------
        PARAMI%CSUBG_MF_PDF = 'NONE'  ! No mass flux PDF
        PARAMI%LOCND2 = .FALSE.       ! OCND2 option

        !-----------------------------------------------------------------------
        ! Initialize budget configuration (disable budgets)
        !-----------------------------------------------------------------------
        BUCONF%LBUDGET_TH = .FALSE.
        BUCONF%LBUDGET_RV = .FALSE.
        BUCONF%LBUDGET_RC = .FALSE.
        BUCONF%LBUDGET_RI = .FALSE.

        !-----------------------------------------------------------------------
        ! Allocate additional required arrays ON GPU
        !-----------------------------------------------------------------------
        ALLOCATE(PRHODJ(nlon, nlev))
        ALLOCATE(PZZ(nlon, nlev))
        ALLOCATE(PMFCONV(nlon, nlev))
        ALLOCATE(PWEIGHT_MF_CLOUD(nlon, nlev))
        ALLOCATE(PSSIO(nlon, nlev))
        ALLOCATE(PSSIU(nlon, nlev))
        ALLOCATE(PIFR(nlon, nlev))
        ALLOCATE(PSRCS(nlon, nlev))

        !-----------------------------------------------------------------------
        ! Initialize arrays (on GPU)
        !-----------------------------------------------------------------------
        !$acc data create(PRHODJ, PZZ, PMFCONV, PWEIGHT_MF_CLOUD, PSSIO, PSSIU, PIFR, PSRCS) &
        !$acc&     deviceptr(f_sigqsat, f_pabs, f_sigs, f_th, f_exn, f_exn_ref, f_rho_dry_ref) &
        !$acc&     deviceptr(f_rv, f_rc, f_ri, f_rr, f_rs, f_rg) &
        !$acc&     deviceptr(f_cf_mf, f_rc_mf, f_ri_mf) &
        !$acc&     deviceptr(f_rvs, f_rcs, f_ris, f_ths) &
        !$acc&     deviceptr(f_cldfr, f_icldfr, f_wcldfr)

        ! Initialize on GPU
        !$acc kernels
        PRHODJ = f_rho_dry_ref
        PZZ = 0.0_C_FLOAT
        PMFCONV = 0.0_C_FLOAT
        PWEIGHT_MF_CLOUD = 0.0_C_FLOAT
        PSSIO = 0.0_C_FLOAT
        PSSIU = 0.0_C_FLOAT
        PIFR = 0.0_C_FLOAT
        PSRCS = 0.0_C_FLOAT
        !$acc end kernels

        LMFCONV = .FALSE.
        OCOMPUTE_SRC = .FALSE.

        !-----------------------------------------------------------------------
        ! Call the GPU-accelerated ICE_ADJUST routine
        !-----------------------------------------------------------------------
        CALL ICE_ADJUST(                                                       &
            D, CST, ICEP, NEBN, TURBN, PARAMI, BUCONF, krr,                    &
            'BRID',                                                            &
            timestep, f_sigqsat,                                               &
            PRHODJ, f_exn_ref, f_rho_dry_ref, f_sigs, LMFCONV, PMFCONV,       &
            f_pabs, PZZ,                                                       &
            f_exn, f_cf_mf, f_rc_mf, f_ri_mf, PWEIGHT_MF_CLOUD,               &
            f_icldfr, f_wcldfr, PSSIO, PSSIU, PIFR,                           &
            f_rv, f_rc, f_rvs, f_rcs, f_th, f_ths,                            &
            OCOMPUTE_SRC, PSRCS, f_cldfr,                                     &
            f_rr, f_ri, f_ris, f_rs, f_rg, TBUDGETS, 0                        &
        )

        !$acc end data

        !-----------------------------------------------------------------------
        ! Cleanup
        !-----------------------------------------------------------------------
        DEALLOCATE(PRHODJ, PZZ, PMFCONV, PWEIGHT_MF_CLOUD)
        DEALLOCATE(PSSIO, PSSIU, PIFR, PSRCS)

    END SUBROUTINE c_ice_adjust_acc_wrap

    !===========================================================================
    ! C-callable GPU wrapper for RAIN_ICE_ACC
    !
    ! This subroutine expects ALL input arrays to be GPU device pointers.
    ! Use CuPy or other GPU array library to allocate and manage GPU memory.
    !
    ! Example Python/CuPy usage:
    !   import cupy as cp
    !   tht_gpu = cp.asarray(tht_cpu, dtype=np.float32)
    !   # Pass tht_gpu.data.ptr to this function
    !===========================================================================
    SUBROUTINE c_rain_ice_acc_wrap(                                            &
        nlon, nlev, krr, timestep,                                             &
        ptr_exn, ptr_dzz, ptr_rhodj, ptr_rhodref, ptr_exnref, ptr_pabs,       &
        ptr_cit, ptr_cldfr, ptr_icldfr, ptr_ssio, ptr_ssiu, ptr_ifr,          &
        ptr_hlc_hrc, ptr_hlc_hcf, ptr_hli_hri, ptr_hli_hcf,                   &
        ptr_tht, ptr_rvt, ptr_rct, ptr_rrt, ptr_rit, ptr_rst, ptr_rgt,        &
        ptr_ths, ptr_rvs, ptr_rcs, ptr_rrs, ptr_ris, ptr_rss, ptr_rgs,        &
        ptr_inprc, ptr_inprr, ptr_evap3d, ptr_inprs, ptr_inprg, ptr_indep,    &
        ptr_rainfr, ptr_sigs, ptr_sea, ptr_town, ptr_conc3d,                  &
        ptr_rht, ptr_rhs, ptr_inprh, ptr_fpr                                  &
    ) BIND(C, name="c_rain_ice_acc")

        !-----------------------------------------------------------------------
        ! Arguments (C-compatible)
        !-----------------------------------------------------------------------
        INTEGER(C_INT), VALUE, INTENT(IN) :: nlon, nlev, krr
        REAL(C_FLOAT), VALUE, INTENT(IN) :: timestep

        ! GPU device pointers for input arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_exn          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_dzz          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rhodj        ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rhodref      ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_exnref       ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_pabs         ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_cit          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_cldfr        ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_icldfr       ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ssio         ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ssiu         ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ifr          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_hlc_hrc      ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_hlc_hcf      ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_hli_hri      ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_hli_hcf      ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_tht          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rvt          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rct          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rrt          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rit          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rst          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rgt          ! 2D: (nlon, nlev) - GPU

        ! GPU device pointers for input/output tendency arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ths          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rvs          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rcs          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rrs          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_ris          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rss          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rgs          ! 2D: (nlon, nlev) - GPU

        ! GPU device pointers for output arrays
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_inprc        ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_inprr        ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_evap3d       ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_inprs        ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_inprg        ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_indep        ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rainfr       ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_sigs         ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_sea          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_town         ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_conc3d       ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rht          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_rhs          ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_inprh        ! 2D: (nlon, nlev) - GPU
        TYPE(C_PTR), VALUE, INTENT(IN) :: ptr_fpr          ! 3D: (nlon, nlev, krr) - GPU

        !-----------------------------------------------------------------------
        ! Fortran pointers to GPU data (using C_FLOAT for single precision)
        !-----------------------------------------------------------------------
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_exn, f_dzz, f_rhodj, f_rhodref
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_exnref, f_pabs, f_cit, f_cldfr
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_icldfr, f_ssio, f_ssiu, f_ifr
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_hlc_hrc, f_hlc_hcf
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_hli_hri, f_hli_hcf
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_tht, f_rvt, f_rct, f_rrt
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_rit, f_rst, f_rgt
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_ths, f_rvs, f_rcs, f_rrs
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_ris, f_rss, f_rgs
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_inprc, f_inprr, f_evap3d
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_inprs, f_inprg, f_indep
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_rainfr, f_sigs, f_sea, f_town
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:) :: f_conc3d, f_rht, f_rhs, f_inprh
        REAL(KIND=C_FLOAT), POINTER, DIMENSION(:,:,:) :: f_fpr

        !-----------------------------------------------------------------------
        ! Local variables for PHYEX structures
        !-----------------------------------------------------------------------
        TYPE(DIMPHYEX_t) :: D
        TYPE(RAIN_ICE_PARAM_t) :: ICEP
        TYPE(RAIN_ICE_DESCR_t) :: ICED
        TYPE(PARAM_ICE_t) :: PARAMI
        TYPE(TBUDGETCONF_t) :: BUCONF
        TYPE(TBUDGETDATA_PTR), DIMENSION(0) :: TBUDGETS

        !-----------------------------------------------------------------------
        ! Convert C pointers to Fortran pointers
        !-----------------------------------------------------------------------
        CALL C_F_POINTER(ptr_exn, f_exn, [nlon, nlev])
        CALL C_F_POINTER(ptr_dzz, f_dzz, [nlon, nlev])
        CALL C_F_POINTER(ptr_rhodj, f_rhodj, [nlon, nlev])
        CALL C_F_POINTER(ptr_rhodref, f_rhodref, [nlon, nlev])
        CALL C_F_POINTER(ptr_exnref, f_exnref, [nlon, nlev])
        CALL C_F_POINTER(ptr_pabs, f_pabs, [nlon, nlev])
        CALL C_F_POINTER(ptr_cit, f_cit, [nlon, nlev])
        CALL C_F_POINTER(ptr_cldfr, f_cldfr, [nlon, nlev])
        CALL C_F_POINTER(ptr_icldfr, f_icldfr, [nlon, nlev])
        CALL C_F_POINTER(ptr_ssio, f_ssio, [nlon, nlev])
        CALL C_F_POINTER(ptr_ssiu, f_ssiu, [nlon, nlev])
        CALL C_F_POINTER(ptr_ifr, f_ifr, [nlon, nlev])
        CALL C_F_POINTER(ptr_hlc_hrc, f_hlc_hrc, [nlon, nlev])
        CALL C_F_POINTER(ptr_hlc_hcf, f_hlc_hcf, [nlon, nlev])
        CALL C_F_POINTER(ptr_hli_hri, f_hli_hri, [nlon, nlev])
        CALL C_F_POINTER(ptr_hli_hcf, f_hli_hcf, [nlon, nlev])
        CALL C_F_POINTER(ptr_tht, f_tht, [nlon, nlev])
        CALL C_F_POINTER(ptr_rvt, f_rvt, [nlon, nlev])
        CALL C_F_POINTER(ptr_rct, f_rct, [nlon, nlev])
        CALL C_F_POINTER(ptr_rrt, f_rrt, [nlon, nlev])
        CALL C_F_POINTER(ptr_rit, f_rit, [nlon, nlev])
        CALL C_F_POINTER(ptr_rst, f_rst, [nlon, nlev])
        CALL C_F_POINTER(ptr_rgt, f_rgt, [nlon, nlev])
        CALL C_F_POINTER(ptr_ths, f_ths, [nlon, nlev])
        CALL C_F_POINTER(ptr_rvs, f_rvs, [nlon, nlev])
        CALL C_F_POINTER(ptr_rcs, f_rcs, [nlon, nlev])
        CALL C_F_POINTER(ptr_rrs, f_rrs, [nlon, nlev])
        CALL C_F_POINTER(ptr_ris, f_ris, [nlon, nlev])
        CALL C_F_POINTER(ptr_rss, f_rss, [nlon, nlev])
        CALL C_F_POINTER(ptr_rgs, f_rgs, [nlon, nlev])
        CALL C_F_POINTER(ptr_inprc, f_inprc, [nlon, nlev])
        CALL C_F_POINTER(ptr_inprr, f_inprr, [nlon, nlev])
        CALL C_F_POINTER(ptr_evap3d, f_evap3d, [nlon, nlev])
        CALL C_F_POINTER(ptr_inprs, f_inprs, [nlon, nlev])
        CALL C_F_POINTER(ptr_inprg, f_inprg, [nlon, nlev])
        CALL C_F_POINTER(ptr_indep, f_indep, [nlon, nlev])
        CALL C_F_POINTER(ptr_rainfr, f_rainfr, [nlon, nlev])
        CALL C_F_POINTER(ptr_sigs, f_sigs, [nlon, nlev])
        CALL C_F_POINTER(ptr_sea, f_sea, [nlon, nlev])
        CALL C_F_POINTER(ptr_town, f_town, [nlon, nlev])
        CALL C_F_POINTER(ptr_conc3d, f_conc3d, [nlon, nlev])
        CALL C_F_POINTER(ptr_rht, f_rht, [nlon, nlev])
        CALL C_F_POINTER(ptr_rhs, f_rhs, [nlon, nlev])
        CALL C_F_POINTER(ptr_inprh, f_inprh, [nlon, nlev])
        CALL C_F_POINTER(ptr_fpr, f_fpr, [nlon, nlev, krr])

        !-----------------------------------------------------------------------
        ! Initialize DIMPHYEX structure
        !-----------------------------------------------------------------------
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

        !-----------------------------------------------------------------------
        ! Initialize physical constants (uses global CST module)
        !-----------------------------------------------------------------------
        CALL INI_CST()

        !-----------------------------------------------------------------------
        ! Initialize RAIN_ICE parameters
        !-----------------------------------------------------------------------
        CALL INI_RAIN_ICE(ICEP, ICED, PARAMI)

        !-----------------------------------------------------------------------
        ! Initialize budget configuration (disable budgets)
        !-----------------------------------------------------------------------
        BUCONF%LBUDGET_TH = .FALSE.
        BUCONF%LBUDGET_RV = .FALSE.
        BUCONF%LBUDGET_RC = .FALSE.
        BUCONF%LBUDGET_RR = .FALSE.
        BUCONF%LBUDGET_RI = .FALSE.
        BUCONF%LBUDGET_RS = .FALSE.
        BUCONF%LBUDGET_RG = .FALSE.

        !-----------------------------------------------------------------------
        ! Call the GPU-accelerated RAIN_ICE routine
        !-----------------------------------------------------------------------
        !$acc data deviceptr(f_exn, f_dzz, f_rhodj, f_rhodref, f_exnref, f_pabs) &
        !$acc&     deviceptr(f_cit, f_cldfr, f_icldfr, f_ssio, f_ssiu, f_ifr) &
        !$acc&     deviceptr(f_hlc_hrc, f_hlc_hcf, f_hli_hri, f_hli_hcf) &
        !$acc&     deviceptr(f_tht, f_rvt, f_rct, f_rrt, f_rit, f_rst, f_rgt) &
        !$acc&     deviceptr(f_ths, f_rvs, f_rcs, f_rrs, f_ris, f_rss, f_rgs) &
        !$acc&     deviceptr(f_inprc, f_inprr, f_evap3d, f_inprs, f_inprg, f_indep) &
        !$acc&     deviceptr(f_rainfr, f_sigs, f_sea, f_town, f_conc3d) &
        !$acc&     deviceptr(f_rht, f_rhs, f_inprh, f_fpr)

        CALL RAIN_ICE(                                                         &
            D, CST, PARAMI, ICEP, ICED, BUCONF,                                &
            timestep, krr, f_exn,                                              &
            f_dzz, f_rhodj, f_rhodref, f_exnref, f_pabs, f_cit, f_cldfr,      &
            f_icldfr, f_ssio, f_ssiu, f_ifr,                                   &
            f_hlc_hrc, f_hlc_hcf, f_hli_hri, f_hli_hcf,                        &
            f_tht, f_rvt, f_rct, f_rrt, f_rit, f_rst,                          &
            f_rgt, f_ths, f_rvs, f_rcs, f_rrs, f_ris, f_rss, f_rgs,            &
            f_inprc, f_inprr, f_evap3d,                                        &
            f_inprs, f_inprg, f_indep, f_rainfr, f_sigs,                       &
            TBUDGETS, 0,                                                       &
            f_sea, f_town, f_conc3d,                                           &
            f_rht, f_rhs, f_inprh, f_fpr                                       &
        )

        !$acc end data

    END SUBROUTINE c_rain_ice_acc_wrap

END MODULE phyex_bridge_acc
