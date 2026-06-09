MODULE physiqex_mod
   IMPLICIT NONE
   CONTAINS

   ! ------------------------------------------------------------------------------------------------
   !
   ! Main authors : F. Hourdin, S. Riette, Q. Libois, A. Idelkadi
   ! --------------
   !
   ! Object : A version of the LMDZ physics that include PHYEX 
   ! -------  Tu be used in particular to be coupled to non hydrostatic version of the dynamical
   !          core dynamico or its limited area version.
   !
   !
   ! Development steps :
   ! -------------------
   !
   ! physiqex_mod.F90 has been introducued in phylmd as an empty driver for the physics (with only
   !    a newtonian relaxation for temperature and drag at the surface).
   ! This driver has been used to introduce the PHYEX package, an externailzed physics package coming
   !    from MesoNH, representing turbulence, convection, clouds and precipitations.
   ! This work has been done for a large part in Dephy sessions (Oleron, 2022, Fréjus 2023, Lège 2024).
   !
   ! Fréjus 2023 : first version runing in 1D for Arm cumulus
   ! Lège 2024   : first versions with a generic soil model (developped at this occasion) and ECrad
   ! Paris 2026  : Updates following the rewriting of the LMDZ physical model and the interface 
   !               between LMDZ and the radiative transfer codes
   !
   ! To do list
   ! ---------
   ! To do list pour un branchement plus propre
   ! * PHYEX considère toutes les espèces microphysiques et la tke comme pronostiques, des termes de tendances 
   !   sont calculés. L'avance temporelle est faite en fin de pas de temps mais pourrait être déplacée au dessus si
   !   les tendances de ces variables passaient par l'interface
   ! * PHYEX a besoin de variables avec mémoire d'un pas de temps à l'autre (en plus des variables pronostiques
   !   décrites juste au dessus). Ce sont les tableaux ALLOCATABLE. Ces variables pourraient devenir
   !   des traceurs.
   ! * La variable ZDZMIN est ici calculée avec un MINVAL qui conduira à des résultats différents
   !   lorsque la distribution (noeuds/procs) sera différente. Cette variable est utile pour déterminer
   !   un critère CFL pour la sédimentation des précipitations. Il suffit donc d'avoir une valeur
   !   approchante.
   ! * Certains ALLOCATABLE sont en klev, d'autres en klev+2. Il est possible de changer ceci mais il faut gérer
   !   les recopies pour que les params voient des tableaux en klev+2 (en effectuant des recopies des
   !   niveaux extrêmes comme fait pour le vent par exemple)
   ! * L'eau tombée en surface (précipitations, sédimentation du nuage, terme de dépôt) se trouve dans
   !   les variables ZINPRC, ZINPRR, ZINPRS, ZINPRG
   !
   ! * Interface avec le rayonnement (A. Idelkadi) :
   !   - gestion des pas de temps (itap, itaprad)
   !   - albedos à passer au rayonnement (alb_dir/alb_dif)
   !   - temperature au sol 
   !   - ozone 
   !   - orbite 
   !   - proprietes optiques nuages:
   !       epaisseur optique, emissivite : cldfrarad, cldemirad, cldtaurad , non necessaire pour ecrad
   !       rayons effectifs : cas ecrad => calcules en utilisant la routine de l'IFS en fonct des contenus en eau 
   !       cloud liq/Ice water mixing contenent (kg/kg) ?
   !   - thermo : zqsat
   !       zqsat specific humidity at liq saturation => 
   !       calcul de Fred ou a  calculer dans ecrad enf fct de t et p
   !       flwc, fiwc : qx( 2) et qx(  3) 
   !   - aerosols ???
   !
   ! * Declarations (A. Abderrahamne) :
   !   - verifier les modules (supprimer les use non necessaires / rajouter ONLY)
   !   - declarations a revoir (IN/OUT/LOCALES)
   !-------------------------------------------------------------------------------------------------


   SUBROUTINE physiqex (nlon,nlev, &
     &            debut,lafin,pdtphys_, &
     &            paprs,pplay,pphi,pphis,presnivs, &
     &            u,v,rot,t,qx, &
     &            flxmass_w, &
     &            d_u, d_v, d_t, d_qx, d_ps)

      USE dimphy, only : klon,klev
      USE infotrac_phy, only : nqtot
      USE geometry_mod, only : latitude, cell_area, latitude_deg, longitude_deg
!      USE comcstphy, only : rg
      USE ioipsl, only : ymds2ju
      !USE phys_state_var_mod, only : phys_state_var_init
      USE phys_state_var_mod
      USE phyetat0_mod, only: phyetat0
      USE output_physiqex_mod, ONLY: output_physiqex
      USE yomcst_mod_h
      USE clesphys_mod_h
      USE flux_arp_mod_h
      USE time_phylmdz_mod, only: current_time, itau_phy, pdtphys, raz_date, update_time
      USE limit_read_mod, ONLY : init_limit_read
      USE conf_phys_m, only: conf_phys
      USE ozonecm_m, only: ozonecm
      USE aero_mod
      USE phys_cal_mod, only: year_len, mth_len, days_elapsed, jh_1jan, &
         year_cur, mth_cur,jD_cur, jH_cur, jD_ref, day_cur, hour, calend
      USE lmdz_call_radiation, only: call_radiation
      USE MODD_CST
      USE phys_output_var_mod, ONLY : phys_output_var_init, &
         cloud_cover_sw, ZFLUX_DIR, ZFLUX_DIR_CLEAR, ZFLUX_DIR_SUN
      USE write_field_phy
      USE phyaqua_mod, only: zenang_an
!      
      ! PHYEX internal modules
      USE MODE_INIT_PHYEX, ONLY: INIT_PHYEX, FILL_DIMPHYEX
      USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
      USE MODD_PHYEX,      ONLY: PHYEX_t
      USE MODD_MISC_LMDZ,  ONLY: MISC_LMDZ_t
      USE MODI_INI_PHYEX,  ONLY: INI_PHYEX
      USE MODI_ICE_ADJUST, ONLY: ICE_ADJUST
      USE MODI_RAIN_ICE, ONLY: RAIN_ICE
      USE MODI_TURB
      USE MODI_SHALLOW_MF

      USE MODD_CST
      USE ioipsl_getin_p_mod, ONLY : getin_p
      USE generic_soil_ini, ONLY : soil_ini
      USE generic_soil_conduct, ONLY : soil_conduct
      USE generic_soil_fluxes, ONLY : soil_fluxes

      USE orbite_mod, ONLY: angle, orbite, zenang
      USE conf_phys_m, ONLY : conf_phys
      ! TBD 2026/03/13 : les variables ci dessous ont été passées en module.
      !                  elles foivent disparaitre de physiq_mod.F90
      USE conf_phys_m, ONLY : ok_journe, ok_mensuel, ok_instan, ok_hf, ok_LES
      USE conf_phys_m, ONLY : ok_volcan, flag_volc_surfstrat, iflag_radia, facttemps, fact_cldcon, iflag_cld_th
      USE conf_phys_m , ONLY : ok_ade, ok_aie, ok_alw, ok_cdnc, bl95_b0, bl95_b1
      USE conf_phys_m , ONLY : aerosol_couple
      USE conf_phys_m , ONLY : chemistry_couple
      USE conf_phys_m , ONLY : flag_aerosol, read_climoz
      USE conf_phys_m , ONLY : flag_bc_internal_mixture
      USE conf_phys_m , ONLY : solarlong0,alp_offset,flag_aer_feedback,flag_aerosol_strat

      IMPLICIT none
!
! Routine argument:
!
! --- variables input ---------------------------------------------------------
      INTEGER, INTENT(in) :: nlon ! number of atmospheric colums
      INTEGER, INTENT(in) :: nlev ! number of vertical levels (should be =klev)
      LOGICAL, INTENT(in) :: debut ! signals first call to physics
      LOGICAL, INTENT(in) :: lafin ! signals last call to physics
      REAL, INTENT(in) :: pdtphys_ ! physics time step (s)
      REAL, INTENT(in) :: paprs(klon,klev+1) ! interlayer pressure (Pa)
      REAL, INTENT(in) :: pplay(klon,klev) ! mid-layer pressure (Pa)
      REAL, INTENT(in) :: pphi(klon,klev) ! geopotential at mid-layer
      REAL, INTENT(in) :: pphis(klon) ! surface geopotential
      REAL, INTENT(in) :: presnivs(klev) ! pseudo-pressure (Pa) of mid-layers
      REAL, INTENT(in) :: u(klon,klev) ! eastward zonal wind (m/s)
      REAL, INTENT(in) :: v(klon,klev) ! northward meridional wind (m/s)
      REAL, INTENT(in) :: rot(klon,klev) ! northward meridional wind (m/s)
      REAL, INTENT(in) :: t(klon,klev) ! temperature (K)
      REAL, INTENT(in) :: qx(klon,klev,nqtot) ! tracers (.../kg_air)
      REAL, INTENT(in) :: flxmass_w(klon,klev) ! vertical mass flux
!
! --- variables output --------------------------------------------------
      REAL, INTENT(out) :: d_u(klon,klev) ! physics tendency on u (m/s/s)
      REAL, INTENT(out) :: d_v(klon,klev) ! physics tendency on v (m/s/s)
      REAL, INTENT(out) :: d_t(klon,klev) ! physics tendency on t (K/s)
      REAL, INTENT(out) :: d_qx(klon,klev,nqtot) ! physics tendency on tracers
      REAL, INTENT(out) :: d_ps(klon) ! physics tendency on surface pressure
!
! --- variables locales -----------------------------------------------      
      REAL :: d_qr(klon, klev), d_qs(klon, klev), d_qg(klon, klev) ! tendency for rain, snow, graupel
      REAL :: d_tke(klon, klev)

      ! --- A. Idelkadi 04-2024 a nettoyer ------------------
      CHARACTER (LEN=20) :: modname='physiq_mod'
      CHARACTER*80 abort_message
      ! Declarations pas de temps
      INTEGER, SAVE :: itap=0         ! compteur pour la physique
      !$OMP THREADPRIVATE(itap)
      INTEGER lmt_pas
      SAVE lmt_pas                ! frequence de mise a jour
      !$OMP THREADPRIVATE(lmt_pas)
      REAL zdtime, zdtime1, zdtime2, zlongi
      REAL zzzo
      !REAL :: ro3i ! 0<=ro3i<=360 ; required time index in NetCDF file for
                 ! the ozone fields, old method.
      ! Declarations pour le rayonnement
      ! Le rayonnement n'est pas calcule tous les pas, il faut donc
      ! sauvegarder les sorties du rayonnement
      ! Surface
      REAL zxtsol(KLON)
      INTEGER itaprad
      SAVE itaprad
      !$OMP THREADPRIVATE(itaprad)
      ! Date equinoxe de PRINTemps
      !LOGICAL, PARAMETER ::  ok_rad=.true. ! a nettoyer
      ! 6 bandes SW : argument non utilise dans radlwsw a nettoyer ?
      REAL,DIMENSION(6), PARAMETER :: SFRWL = (/ 1.28432794E-03, 0.12304168, 0.33106142, &
                                         & 0.32870591, 0.18568763, 3.02191470E-02 /)
      LOGICAL, PARAMETER :: new_orbit = .TRUE.
      INTEGER, PARAMETER :: mth_eq=3, day_eq=21
      REAL :: day_since_equinox, jD_eq
      !REAL zdtime, zdtime1, zdtime2, zlongi
      REAL dist, rmu0(klon), fract(klon), zrmu0(klon), zfract(klon)
      CHARACTER(len=512) :: namelist_ecrad_file
      ! 
      REAL, DIMENSION(klon, klev) :: zqsat
      REAL, DIMENSION(klon, klev) :: ref_liq, ref_ice, ref_liq_pi, ref_ice_pi ! rayons effectifs des gout
      REAL, DIMENSION(klon, klev) :: cldtaurad   ! epaisseur optique
      REAL, DIMENSION(klon, klev) :: cldtaupirad ! epaisseur optique
      REAL, DIMENSION(klon, klev) :: cldemirad   ! emissivite pour
      REAL, DIMENSION(klon, klev) :: cldfrarad   ! fraction nuageuse
      ! aerosols AI attention a les declarer ailleurs dans phylmd : phys_output ou phy_stats et non phys_local
      REAL, DIMENSION(klon,klev+1) :: ZLWFT0_i, ZSWFT0_i,ZFLDN0, ZFLUP0, ZFSDN0, ZFSUP0
      REAL, DIMENSION(klon)        :: oplwad_aero, sollwad_aero,&
                                      toplwai_aero, sollwai_aero, &
                                      toplwad0_aero, sollwad0_aero
      REAL, DIMENSION(klon,9)      :: solsw_aero, solsw0_aero,topsw_aero, topsw0_aero
      REAL, DIMENSION(klon)        :: &
                                     topswcf_aero, solswcf_aero, &
                                     solswad_aero, &
                                     solswad0_aero, &
                                     topswai_aero, solswai_aero, &
                                     toplwad_aero, topswad0_aero, &
                                     topswad_aero
      REAL, DIMENSION(klon,klev,naero_tot) :: m_allaer

      REAL,DIMENSION(klon,klev)   :: d_t_swr, d_t_sw0, d_t_lwr, d_t_lw0
      !
      EXTERNAL suphel ! AI attention 
      !
      ! --- Fin A. Idelkadi ----------------

      INTEGER        length
      PARAMETER    ( length = 100 )
      REAL tabcntr0( length       )
      INTEGER, PARAMETER :: longcles=20
      REAL, SAVE :: clesphy0(longcles)
      !$OMP THREADPRIVATE(clesphy0)

! ---------------  Saved variables ----------------------------------------------------
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: PSIGS !variance of s
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: PCF_MF, PRC_MF, PRI_MF, PWEIGHT_MF_CLOUD !shallow convection cloud
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: PHLC_HCF_MF, PHLC_HRC_MF, PHLI_HCF_MF, PHLI_HRI_MF
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: ZQR, ZQS, ZQG !rain, snow, graupel specifiq contents
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE ::  PTKEM       ! TKE
      TYPE(DIMPHYEX_t),SAVE    :: D
      TYPE(PHYEX_t), SAVE      :: PHYEX
      TYPE(MISC_LMDZ_t), TARGET, SAVE :: MISC_LMDZ
      !
      INTEGER, PARAMETER       :: KRR=6, KRRL=2, KRRI=3, KSV=0
      INTEGER                  :: JRR
      ! Time-State variables and Sources variables
      REAL, DIMENSION(klon,klev+2)     ::  ZUT, ZVT ! wind component on klev+2
      REAL, DIMENSION(klon,klev+2)     ::  ZPABST   ! absolute pressure
      REAL, DIMENSION(klon,klev+2,KRR) ::  ZRX,ZRXS,ZRXS0 ! qx and source of q from LMDZ to rt
      REAL, DIMENSION(klon,klev+2)     ::  ZRHOD    ! rho_dry
      REAL, DIMENSION(klon,klev+2,0)   ::  ZSVT     ! passive scal. var.
      REAL, DIMENSION(klon,klev+2)     :: PWT ! vertical wind velocity (used only for diagnostic)
      REAL, DIMENSION(klon,klev+2) ::  ZRUS,ZRVS,ZRWS,ZRTHS,ZRTKES ! sources of momentum, conservative potential temperature, Turb. Kin. Energy
      REAL, DIMENSION(KLON,KLEV+2) :: ZTHETAS !tendency
      REAL, DIMENSION(KLON,KLEV+2) :: ZTHETAS0
      REAL, DIMENSION(KLON,KLEV+2) :: ZTKES
      REAL, DIMENSION(KLON,KLEV+2) :: ZTKES0
      ! Source terms for all water kinds, PRRS(:,:,:,1) is used for the conservative ! mixing ratio
      REAL, DIMENSION(klon,klev+2,KRR) ::  ZRRS
      REAL, DIMENSION(klon,klev+2)   ::  PTHVREF  ! Virtual Potential Temperature of the reference state
      ! Adjustment variables
      REAL, DIMENSION(KLON,KLEV+2) :: zqdm !1-qt=1/(1+rt)
      REAL, DIMENSION(KLON,KLEV+2) :: zqt !mixing ratio and total specifiq content
      REAL, DIMENSION(KLON,KLEV+2) :: ZTHETA, ZEXN !theta and exner function
      REAL, DIMENSION(KLON,KLEV+2) :: zz_mass !altitude above ground of mass points
      REAL, DIMENSION(KLON,KLEV+2) :: zz_flux !altitude above ground of flux points
      REAL, DIMENSION(KLON,KLEV+2) :: zdzm !distance between two consecutive mass points (expressed on flux points)
      REAL, DIMENSION(KLON,KLEV+2) :: zdzf !distance between two consecutive flux points (expressed on mass points)
      REAL, DIMENSION(KLON) :: zs !surface orography
      REAL, DIMENSION(KLON) :: zsigqsat !coeff for the extra term for s variance
      REAL, DIMENSION(0) :: ZMFCONV
      REAL, DIMENSION(KLON,KLEV+2) :: ZICLDFR, ZWCLDFR, ZSSIO, ZSSIU, ZIFR, ZICE_CLD_WGT !used only in HIRLAM config
      REAL, DIMENSION(KLON,KLEV+2) :: ZSRC ! bar(s'rc')
      REAL, DIMENSION(KLON,KLEV+2) :: ZCLDFR !cloud fraction
      REAL, DIMENSION(KLON,KLEV+2) :: ZHLC_HRC, ZHLC_HCF, ZHLI_HRI, ZHLI_HCF !subgrid autoconversion
      ! Rain_ice variables
      REAL :: ZDZMIN
      REAL, DIMENSION(klon, klev+2) :: ZCIT !ice concentration
      REAL, DIMENSION(klon) :: ZINPRC, ZINPRR, ZINPRS, ZINPRG !precipitation flux at ground
      REAL, DIMENSION(klon, klev+2) :: ZEVAP3D !evaporation (diag)
      REAL, DIMENSION(klon, klev+2) :: ZRAINFR
      REAL, DIMENSION(klon) :: ZINDEP !deposition flux, already contained in ZINPRC
      REAL, DIMENSION(klon)         :: ZSEA, ZTOWN !sea and town fractions in the frid cell
      ! Turbulence variables
      REAL, DIMENSION(klon,klev+2)  ::  PDXX,PDYY,PDZX,PDZY ! metric coefficients
      REAL, DIMENSION(klon,klev+2)  ::  PZZ       !  physical distance between 2 succesive grid points along the K direction
      REAL, DIMENSION(klon)         ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW ! Director Cosinus along x, y and z directions at surface w-point
      REAL, DIMENSION(klon)         ::  PCOSSLOPE       ! cosinus of the angle between i and the slope vector
      REAL, DIMENSION(klon)         ::  PSINSLOPE       ! sinus of the angle   between i and the slope vector
      REAL, DIMENSION(klon,klev+2)  ::  PRHODJ   ! dry density * Grid size
      REAL, DIMENSION(0,0)          ::  MFMOIST  ! moist mass flux dual scheme
      REAL, DIMENSION(0,0,0)        ::  PHGRADLEO   ! horizontal gradients
      REAL, DIMENSION(0,0,0)        ::  PHGRADGOG
      REAL, DIMENSION(klon)         ::  PSFTH,PSFRV,PSFU,PSFV ! normal surface fluxes of theta, Rv, (u,v) parallel to the orography
      REAL, DIMENSION(klon,0)       ::  PSFSV ! normal surface fluxes of Scalar var. KSV=0
      REAL, DIMENSION(klon)         ::  ZBL_DEPTH  ! BL height for TOMS
      REAL, DIMENSION(klon)         ::  ZSBL_DEPTH ! SBL depth for RMC01
      REAL, DIMENSION(klon,klev+2)  ::  ZCEI ! Cloud Entrainment instability index to emphasize localy turbulent fluxes
      REAL, DIMENSION(klon,klev+2,0)::  ZRSVS ! Source terms for all passive scalar variables
      REAL, DIMENSION(klon,klev+2)  ::  ZFLXZTHVMF ! MF contribution for vert. turb. transport used in the buoy. prod. of TKE
      REAL, DIMENSION(klon,klev+2)  ::  ZWTH       ! heat flux
      REAL, DIMENSION(klon,klev+2)  ::  ZWRC       ! cloud water flux
      REAL, DIMENSION(klon,klev+2,0)::  ZWSV       ! scalar flux
      REAL, DIMENSION(klon,klev+2)  ::  ZTP        ! Thermal TKE production (MassFlux + turb)
      REAL, DIMENSION(klon,klev+2)  ::  ZDP        ! Dynamic TKE production
      REAL, DIMENSION(klon,klev+2)  ::  ZTDIFF     ! Diffusion TKE term
      REAL, DIMENSION(klon,klev+2)  ::  ZTDISS     ! Dissipation TKE term
      REAL, DIMENSION(0,0)          ::  PLENGTHM, PLENGTHH ! length scale from vdfexcu (HARMONIE-AROME)
      REAL :: ZTHVREFZIKB ! for electricity scheme interface
      REAL, DIMENSION(klon,klev+2)  ::  ZDXX,ZDYY,ZDZX,ZDZY,ZZZ
      REAL, DIMENSION(klon)         ::  ZDIRCOSXW,ZDIRCOSYW,ZDIRCOSZW,ZCOSSLOPE,ZSINSLOPE
      REAL, DIMENSION(klon)         ::  ZSEA_UCU, ZSEA_VCU ! oceanic surface current
      ! Shallow variables 
      REAL, DIMENSION(klon,klev+2) ::  PDUDT_MF     ! tendency of U   by massflux scheme
      REAL, DIMENSION(klon,klev+2) ::  PDVDT_MF     ! tendency of V   by massflux scheme
      REAL, DIMENSION(klon,klev+2) ::  PDTHLDT_MF   ! tendency of thl by massflux scheme
      REAL, DIMENSION(klon,klev+2) ::  PDRTDT_MF    ! tendency of rt  by massflux scheme
      REAL, DIMENSION(klon,klev+2,KSV)::  PDSVDT_MF    ! tendency of Sv  by massflux scheme
      REAL, DIMENSION(klon,klev+2) ::  PDTKEDT_MF   ! tendency of TKE by massflux scheme
      REAL, DIMENSION(klon,klev+2) ::  PSIGMF
      REAL, DIMENSION(klon,klev+2) ::  ZFLXZTHMF
      REAL, DIMENSION(klon,klev+2) ::  ZFLXZRMF
      REAL, DIMENSION(klon,klev+2) ::  ZFLXZUMF
      REAL, DIMENSION(klon,klev+2) ::  ZFLXZVMF
      REAL, DIMENSION(klon,klev+2) ::  ZFLXZTKEMF
      REAL, DIMENSION(klon,klev+2) ::  PTHL_UP   ! Thl updraft characteristics
      REAL, DIMENSION(klon,klev+2) ::  PRT_UP    ! Rt  updraft characteristics
      REAL, DIMENSION(klon,klev+2) ::  PRV_UP    ! Vapor updraft characteristics
      REAL, DIMENSION(klon,klev+2) ::  PU_UP     ! U wind updraft characteristics
      REAL, DIMENSION(klon,klev+2) ::  PV_UP     ! V wind updraft characteristics
      REAL, DIMENSION(klon,klev+2) ::  PTKE_UP   ! TKE updraft characteristics
      REAL, DIMENSION(klon,klev+2) ::  PRC_UP    ! cloud content updraft characteristics
      REAL, DIMENSION(klon,klev+2) ::  PRI_UP    ! ice content   updraft characteristics
      REAL, DIMENSION(klon,klev+2) ::  PTHV_UP   ! Thv   updraft characteristics
      REAL, DIMENSION(klon,klev+2) ::  PW_UP     ! vertical speed updraft characteristics
      REAL, DIMENSION(klon,klev+2) ::  PFRAC_UP  ! updraft fraction
      REAL, DIMENSION(klon,klev+2) ::  PEMF      ! updraft mass flux
      REAL, DIMENSION(klon,klev+2) ::  PDETR     ! updraft detrainment
      REAL, DIMENSION(klon,klev+2) ::  PENTR     ! updraft entrainment
      INTEGER,DIMENSION(klon) ::IKLCL,IKETL,IKCTL ! level of LCL,ETL and CTL
      ! Values after saturation adjustement
      REAL, DIMENSION(klon,klev) :: t_adj ! Adjusted temperature
      REAL, DIMENSION(klon,klev) :: qv_adj, ql_adj, qi_adj, qr_adj, qs_adj, qg_adj !specific contents after adjustement
      !
      REAL :: temp_newton(klon,klev)
      INTEGER :: k
      LOGICAL, SAVE :: first=.true.
      !$OMP THREADPRIVATE(first)
      !REAL,save :: rg=9.81
      !!$OMP THREADPRIVATE(rg)

      ! ----------------- For I/Os ---------------------------------
      INTEGER :: itau0
      REAL, DIMENSION(nlon) :: swdnsrf,lwdnsrf,lwupsrf,sensible,latent,qsats,capcal,fluxgrd, fluxu,fluxv,fluxt,fluxq
      INTEGER, PARAMETER :: nsoil=11
      REAL, DIMENSION(nsoil) :: zout
      REAL, DIMENSION(nlon,nlev) :: tsoil_out

      INTEGER :: isoil
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: tsrf
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: tsoil
      !$OMP THREADPRIVATE(tsrf,tsoil)

      REAL :: thermal_inertia0,z0m0,z0h0,z0q0,betasoil

      REAL, SAVE :: stefan, fsw0, timesec
      !$OMP THREADPRIVATE(stefan, fsw0, timesec)
      INTEGER, SAVE :: iecri=1,iflag_surface=0
      !$OMP THREADPRIVATE(iecri,iflag_surface)
      REAL, SAVE :: zjulian
      !$OMP THREADPRIVATE(zjulian)

! =======================================================================================================
! Initialisations de la physique au premier pas de temps
!------------------------------------------------------------
      PRINT*,'Debut physiqex',debut
      ! AI
      pdtphys=pdtphys_
      CALL update_time(pdtphys)
      phys_tstep=NINT(pdtphys)

      ! initializations
      IF (debut) THEN ! Things to do only for the first call to physics 
        PRINT*,'Debut physiqex IN'
        PRINT*,'NOUVEAU PHYSIQEX 222222222222222222222'
        CALL suphel ! initialiser constantes et parametres phys.
        ! Read parameters and flags in files .def  
        CALL conf_phys

        ! IF (.NOT. create_etat0_limit) 
        CALL init_limit_read(days_elapsed)

        ! load initial conditions for physics (including the grid)
        ! AI
        !call phys_state_var_init(1) ! some initializations, required before calling phyetat0
        call phys_state_var_init ! some initializations, required before calling phyetat0
        call phyetat0("startphy.nc", clesphy0, tabcntr0)
        CALL phys_output_var_init
    
        ! AI
        ! Init increments et pas
        !itap    = 0
        itaprad = 0
        lmt_pas = NINT(86400./phys_tstep * 1.0)   ! tous les jours
        WRITE(*,*)'La frequence de lecture surface est de ',  &
                    lmt_pas
        IF (MOD(NINT(86400./phys_tstep),nbapp_rad).EQ.0) THEN
                  radpas = NINT( 86400./phys_tstep)/nbapp_rad
        ELSE
                WRITE(*,*) 'le nombre de pas de temps physique doit etre un ', &
                'multiple de nbapp_rad'
                WRITE(*,*) 'changer nbapp_rad ou alors commenter ce test ', &
                'mais 1+1<>2'
                abort_message='nbre de pas de temps physique n est pas multiple ' &
                // 'de nbapp_rad'
                CALL abort_physic(modname,abort_message,1)
        ENDIF

        ! Initialize outputs:
        itau0=0
        ! compute zjulian for annee0=1979 and month=1 dayref=1 and hour=0.0
        !CALL ymds2ju(annee0, month, dayref, hour, zjulian)
        !call ymds2ju(1979, 1, 1, 0.0, zjulian)
        !PRINT*,'zjulian=',zjulian
  
        ZDZMIN=MINVAL((pphi(:,2:) - pphi(:,1:klev-1))/9.81)
        CALL INIT_PHYEX(pdtphys, ZDZMIN, PHYEX, MISC_LMDZ)
        CALL FILL_DIMPHYEX(KLON, KLEV, D)

        ! Variables SAVEd
        ALLOCATE(PTKEM(klon,klev+2))
        ALLOCATE(PSIGS(klon,klev+2))
        ALLOCATE(PCF_MF(klon,klev+2))
        ALLOCATE(PWEIGHT_MF_CLOUD(klon,klev+2))
        ALLOCATE(PRC_MF(klon,klev+2))
        ALLOCATE(PRI_MF(klon,klev+2))
        ALLOCATE(PHLC_HCF_MF(klon,klev+2))
        ALLOCATE(PHLC_HRC_MF(klon,klev+2))
        ALLOCATE(PHLI_HCF_MF(klon,klev+2))
        ALLOCATE(PHLI_HRI_MF(klon,klev+2))
        ALLOCATE(ZQR(klon, klev))
        ALLOCATE(ZQS(klon, klev))
        ALLOCATE(ZQG(klon, klev))
        PSIGS=0.
        PCF_MF=0.
        PWEIGHT_MF_CLOUD=0.
        PRC_MF=0.
        PRI_MF=0.
        PHLC_HCF_MF=0.
        PHLC_HRC_MF=0.
        PHLI_HCF_MF=0.
        PHLI_HRI_MF=0.
        ZQR=0.
        ZQS=0.
        ZQG=0.
        PTKEM(:,:) = PHYEX%TURBN%XTKEMIN ! TODO: init from TKE at stationnary state

        ! ----------------------------------------------------------------
        ! Initialisation du sol generique
        ! ----------------------------------------------------------------
        allocate(tsrf(nlon))
        allocate(tsoil(nlon,nsoil))
        thermal_inertia0=2000.
        z0m0=0.05
        z0h0=0.05
        z0q0=0.05
        betasoil=0.2
        call soil_ini(klon,nsoil,PHYEX%CST%XCPD,PHYEX%CST%XLVTT,thermal_inertia0,z0m0,z0h0,z0q0,betasoil)
        !
#ifndef CPP_IOIPSL_NO_OUTPUT
       ! Initialize IOIPSL output file
#endif
        ! AI Partie Surface depplacee
        ! compute tendencies to return to the dynamics:
        ! "friction" on the first layer
        !   tsrf(1:nlon)=t(1:nlon,1)
        !   do isoil=1,nsoil
        !      tsoil(1:nlon,isoil)=tsrf(1:nlon)
        !      !tsoil(1:nlon,isoil)=302.
        !   enddo
        !   rpi=2.*asin(1.)
        !   timesec=5*3600
        !   stefan=5.67e-8
        !   fsw0=900.
        !   call getin_p('iecri',iecri)
        !   call getin_p('iflag_surface',iflag_surface)

      ENDIF ! of IF (debut)

      !------------------------------------------------------------
      ! Initialisations a chaque pas de temps
      !------------------------------------------------------------

      ! --- set all tendencies to zero
      d_u(1:klon,1:klev)=0.
      d_v(1:klon,1:klev)=0.
      d_t(1:klon,1:klev)=0.
      d_qx(1:klon,1:klev,1:nqtot)=0.
      d_qr(1:klon,1:klev)=0.
      d_qs(1:klon,1:klev)=0.
      d_qg(1:klon,1:klev)=0.
      d_ps(1:klon)=0.
      d_tke(1:klon,1:klev)=0.
      !
      ! A. Idelkadi attention
      d_t_swr = 0.
      d_t_sw0 = 0.
      d_t_lwr = 0.
      d_t_lw0 = 0.
      !
      ZDXX(:,:) = 0.
      ZDYY(:,:) = 0.
      ZDZX(:,:) = 0.
      ZDZY(:,:) = 0.
      ZDIRCOSXW(:) = 1.
      ZDIRCOSYW(:) = 1.
      ZDIRCOSZW(:) = 1.
      ZCOSSLOPE(:) = 0.
      ZSINSLOPE(:) = 1.
      PHGRADLEO(:,:,:) = 0.
      PHGRADGOG(:,:,:) = 0.
      ZSEA_UCU(:) = 0.
      ZSEA_VCU(:) = 0.
      ZBL_DEPTH(:) = 0. ! needed only with LRMC01 key (correction in the surface boundary layer)
      ZSBL_DEPTH(:) = 0.
      ZCEI(:,:) = 0.  ! needed only IF HTURBLEN_CL /= 'NONE' modification of mixing lengh inside clouds
      ZSVT(:,:,:) = 0.
      PWT(:,:) = 0.
      ZUT(:,2:klev+1) = u(:,:)
      ZVT(:,2:klev+1) = v(:,:)
      !
      !------------------------------------------------------------
      ! Conversions and extra levels
      !------------------------------------------------------------
      !TODO check in Meso-NH how values are extrapolated outside of the physical domain
      zqt(:,2:klev+1) = qx(:,:,1) + qx(:,:,2) + qx(:,:,3)
      zqt(:,1)=0.
      zqt(:,klev+2)=0.
      zqdm(:,:)=1.-zqt(:,:) !equal to 1/(1+rt)
      ZRX(:,2:klev+1,1) = qx(:,:,1) / zqdm(:,2:klev+1)
      ZRX(:,2:klev+1,2) = qx(:,:,2) / zqdm(:,2:klev+1)
      ZRX(:,2:klev+1,4) = qx(:,:,3) / zqdm(:,2:klev+1)
      ZRX(:,2:klev+1,3) = ZQR(:,:)
      ZRX(:,2:klev+1,5) = ZQS(:,:)
      ZRX(:,2:klev+1,6) = ZQG(:,:)
      DO JRR=1,KRR
        CALL VERTICAL_EXTEND(ZRX(:,:,JRR),klev)
      END DO
      !
      ZEXN(:,2:klev+1) = (pplay(:,:) / PHYEX%CST%XP00) ** (PHYEX%CST%XRD/PHYEX%CST%XCPD)
      ZTHETA(:,2:klev+1) = t(:,:) / ZEXN(:,2:klev+1)
      CALL VERTICAL_EXTEND(ZEXN,klev)
      CALL VERTICAL_EXTEND(ZTHETA,klev)

      !TODO check in Meso-NH how zz_mass and zz_flux are initialized outside of the physical domain
      zs(:) = pphis(:)/PHYEX%CST%XG
      zz_mass(:,2:klev+1) = pphi(:,:) / PHYEX%CST%XG
      zz_mass(:,1) = 2*zs-zz_mass(:,2)
      zz_mass(:,klev+2)=2.*zz_mass(:,klev+1)-zz_mass(:,klev)

      do k=2, klev+2
        zz_flux(:,k)=(zz_mass(:,k-1)+zz_mass(:,k))/2.
      enddo
      zz_flux(:,1)=2*zz_mass(:,1)-zz_flux(:,2)

      !zdzf is the distance between two consecutive flux points (expressed on mass points)
      do k=1,klev+1
        zdzf(:,k)=zz_flux(:,k+1)-zz_flux(:,k)
      enddo
      zdzf(:,klev+2)=(zz_mass(:,klev+2)-zz_flux(:,klev+2))*2.

      !zdzm distance between two consecutive mass points (expressed on flux points)
      do k=2,klev+2
        zdzm(:,k)=zz_mass(:,k)-zz_mass(:,k-1)
      enddo
      zdzm(:,1)=(zz_mass(:,1)-zz_flux(:,1))*2.

      ZPABST(:,2:klev+1) = pplay(:,:)
      ZRHOD(:,2:klev+1)=ZPABST(:,2:klev+1)/(t*(PHYEX%CST%XRD+ZRX(:,2:klev+1,1)*PHYEX%CST%XRV))
      DO k=2,klev+1
        PRHODJ(:,k) = ZRHOD(:,k) * (zdzf(:,k)*cell_area(:))
      END DO
      PTHVREF(:,:) = ZTHETA(:,:) * (1. + PHYEX%CST%XRV/PHYEX%CST%XRD * ZRX(:,:,1)) * ZQDM(:,:)

      CALL VERTICAL_EXTEND(ZPABST,klev)
      CALL VERTICAL_EXTEND(PRHODJ,klev)
      CALL VERTICAL_EXTEND(ZRHOD,klev)
      CALL VERTICAL_EXTEND(ZUT,klev)
      CALL VERTICAL_EXTEND(ZVT,klev)
 
      !------------------------------------------------------------
      ! Tendencies
      !------------------------------------------------------------
      !For Meso-NH, initialia values for the tendencies are filled with
      !a pseudo-tendecy computed by dividing the state variable by the time step
      !This mechanism enables the possibility for the different parametrisations
      !to guess the value at the end of the time step.
      !For the wind components, we could do the same way but it is not needed
      !as the parametrisations don't use the S varaible to guess the futur value of the wind.
      ZRXS(:,:,:) = ZRX(:,:,:)/pdtphys
      ZTHETAS(:,:)=ZTHETA(:,:)/pdtphys
      ZTKES(:,:)=PTKEM(:,:)/pdtphys
      !To compute the actual tendency, we SAVE the initial values of these variables
      ZRXS0(:,:,:) = ZRXS(:,:,:)
      ZTHETAS0(:,:)=ZTHETAS(:,:)
      ZTKES0(:,:)=ZTKES(:,:)

! ========================================================================================================
! Adjustment
!---------------------------------------------------------------------------------------------------------
      !
      ZSRC(:,:) = 0.
      ZSIGQSAT=PHYEX%NEBN%VSIGQSAT
      CALL ICE_ADJUST (D, PHYEX%CST, PHYEX%RAIN_ICE_PARAMN, PHYEX%NEBN, PHYEX%TURBN, PHYEX%PARAM_ICEN,    &
                &MISC_LMDZ%TBUCONF, KRR,                                                           &
                &'ADJU',                                                                            &
                &pdtphys, ZSIGQSAT,                                                                 &
                &PRHODJ, ZEXN, ZRHOD, PSIGS, .FALSE., zmfconv,                                      &
                &ZPABST, ZZ_MASS,                                                                   &
                &ZEXN, PCF_MF, PRC_MF, PRI_MF, PWEIGHT_MF_CLOUD,                                    &
                &ZICLDFR, ZWCLDFR, ZSSIO, ZSSIU, ZIFR,                                              &
                &ZRX(:,:,1), ZRX(:,:,2), ZRXS(:,:,1), ZRXS(:,:,2), ZTHETA, ZTHETAS,                 &
                &MISC_LMDZ%COMPUTE_SRC, ZSRC, ZCLDFR,                                              &
                &ZRX(:,:,3), ZRX(:,:,4), ZRXS(:,:,4), ZRX(:,:,5), ZRX(:,:,6),                       &
                &MISC_LMDZ%YLBUDGET, MISC_LMDZ%NBUDGET,                                           &
                &ZICE_CLD_WGT,                                                                      &
                &PHLC_HRC=ZHLC_HRC, PHLC_HCF=ZHLC_HCF, PHLI_HRI=ZHLI_HRI, PHLI_HCF=ZHLI_HCF,        &
                &PHLC_HRC_MF=PHLC_HRC_MF, PHLC_HCF_MF=PHLC_HCF_MF, PHLI_HRI_MF=PHLI_HRI_MF, PHLI_HCF_MF=PHLI_HCF_MF)
      !
      !Variables are updated with their adjusted values (to be used by the other parametrisations)
      ZTHETA(:,:)=ZTHETAS(:,:)*pdtphys
      ZRX(:,:,:)=ZRXS(:,:,:)*pdtphys
      t_adj=ZTHETA(:,2:klev+1)*ZEXN(:,2:klev+1)
      qv_adj=ZRX(:,2:klev+1,1)*zqdm(:,2:klev+1)
      ql_adj=ZRX(:,2:klev+1,2)*zqdm(:,2:klev+1)
      qi_adj=ZRX(:,2:klev+1,4)*zqdm(:,2:klev+1)
      qr_adj=ZRX(:,2:klev+1,3)*zqdm(:,2:klev+1)
      qs_adj=ZRX(:,2:klev+1,5)*zqdm(:,2:klev+1)
      qg_adj=ZRX(:,2:klev+1,6)*zqdm(:,2:klev+1)
      ZRHOD(:,2:klev+1)=ZPABST(:,2:klev+1)/(t*(PHYEX%CST%XRD+ZRX(:,2:klev+1,1)*PHYEX%CST%XRV))
      CALL VERTICAL_EXTEND(ZRHOD,klev)
      !

! ==========================================================================================================
! Radiation
!-----------------------------------------------------------------------------------------------------------
      !
      !A. Idelkadi a nettoyer apres
      ! calcul de qsat Fred mai 2024
      zqsat(:,:) = EXP( XALPW - XBETAW / t(:,:) - XGAMW * ALOG( t(:,:) ) )
      zqsat(:,:) = XRD / XRV * zqsat(:,:) / ( pplay(:,:) - zqsat(:,:))
      !
      IF (iflag_radia.ge.1) THEN
        !
        !ECRAD can be plugged here and can use the adjusted values for temperature and hydrometeors
        !Cloud fraction is available in the ZCLDFR variable
        !
        ! AI attention : dans un 1er temps read_climoz=0
        ! Update ozone IF day change
        !PRINT*,'solarlong0, itap, lmt_pas, days_elapsed : ', solarlong0, itap, lmt_pas, days_elapsed
        !stop
        IF (MOD(itap-1,lmt_pas) == 0) THEN
        !       IF (read_climoz <= 0) THEN
            ! Once per day, update ozone from Royer:
          IF (solarlong0<-999.) THEN
             ! Generic case with evolvoing season
             zzzo=real(days_elapsed+1)
          ELSE IF (abs(solarlong0-1000.)<1.e-4) THEN
             ! Particular case with annual mean insolation
             zzzo=real(90) ! could be revisited
             IF (read_climoz/=-1) THEN
                abort_message ='read_climoz=-1 is recommended when ' &
                               // 'solarlong0=1000.'
                CALL abort_physic (modname,abort_message,1)
             ENDIF
          ELSE
             ! Case where the season is imposed with solarlong0
             zzzo=real(90) ! could be revisited
          ENDIF
          zzzo=real(days_elapsed+1)
          wo(:,:,1)=ozonecm(latitude_deg, paprs,rjour=zzzo)
          !PRINT*,'wo = ', wo
        !       ELSE
        !          !--- ro3i = elapsed days number since current year 1st january, 0h
        !          ro3i=days_elapsed+jh_cur-jh_1jan
        !          !--- scaling for old style files (360 records)
        !          IF(SIZE(time_climoz)==360.AND..NOT.ok_daily_climoz) ro3i=ro3i*360./year_len
        !          IF(adjust_tropopause) THEN
        !             CALL regr_pr_time_av(ncid_climoz, vars_climoz(1:read_climoz),   &
        !                      ro3i, 'C', press_cen_climoz, pplay, wo, paprs(:,1),    &
        !                      time_climoz ,  longitude_deg,   latitude_deg,          &
        !                      dyn_tropopause(t_seri, ztsol, paprs, pplay, rot))
        !          ELSE
        !             CALL regr_pr_time_av(ncid_climoz,  vars_climoz(1:read_climoz),  &
        !                      ro3i, 'C', press_cen_climoz, pplay, wo, paprs(:,1),    &
        !                      time_climoz )
        !          ENDIF
        !          ! Convert from mole fraction of ozone to column density of ozone in a
        !          ! cell, in kDU:
        !          FORALL (l = 1: read_climoz) wo(:, :, l) = wo(:, :, l) * rmo3 / rmd &
        !               * zmasse / dobson_u / 1e3
        !          ! (By regridding ozone values for LMDZ only once a day, we
        !          ! have already neglected the variation of pressure in one
        !          ! day. So do not recompute "wo" at each time step even if
        !          ! "zmasse" changes a little.)
        !       ENDIF
        ENDIF
        !=========================================================================
        ! Calculs de l'orbite.
        ! Necessaires pour le rayonnement et la surface (calcul de l'albedo).
        ! doit donc etre plac\'e avant radlwsw et pbl_surface
        !CALL ymds2ju(year_cur, mth_eq, day_eq,0., jD_eq)
        CALL ymds2ju(year_cur, mth_eq, day_eq,0., jD_eq)
        PRINT*,'mth_eq, day_eq, jD_eq, jH_cur =', mth_eq, day_eq, jD_eq, jH_cur
        day_since_equinox = (jD_cur + jH_cur) - jD_eq
        !choix entre calcul de la longitude solaire vraie ou valeur fixee = solarlong0
        IF (solarlong0<-999.) THEN
          IF (new_orbit) THEN
             ! calcul selon la routine utilisee pour les planetes
             CALL solarlong(day_since_equinox, zlongi, dist)
             PRINT*,'day_since_equinox, zlongi, dist : ', day_since_equinox, zlongi, dist
          ELSE
            ! calcul selon la routine utilisee pour l'AR4
            CALL orbite(REAL(days_elapsed+1),zlongi,dist)
          ENDIF
        ELSE
          zlongi=solarlong0  ! longitude solaire vraie
          dist=1.            ! distance au soleil / moyenne
        ENDIF
        !IF (prt_level.ge.1) write(lunout,*)'Longitude solaire ',zlongi,solarlong0,dist

        ! AI a refaire en definissant les itap, itaprad    
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ============================
        ! Pour une solarlong0=1000., on calcule un ensoleillement moyen sur
        ! l'annee a partir d'une formule analytique.
        ! Cet ensoleillement est sym\'etrique autour de l'\'equateur et
        ! non nul aux poles.
        IF (abs(solarlong0-1000.)<1.e-4) THEN
           CALL zenang_an(iflag_cycle_diurne.GE.1,jH_cur, &
                latitude_deg,longitude_deg,rmu0,fract)
           swradcorr(:) = 1.0
           ! JrNt(:) = 1.0
           zrmu0(:) = rmu0(:)
        ELSE
           ! recode par Olivier Boucher en sept 2015
           SELECT CASE (iflag_cycle_diurne)
           CASE(0)
                !  Sans cycle diurne
                CALL angle(zlongi, latitude_deg, fract, rmu0)
                swradcorr = 1.0
                ! JrNt = 1.0
                zrmu0 = rmu0
           CASE(1)
                !  Avec cycle diurne sans application des poids
                !  bit comparable a l ancienne formulation cycle_diurne=true
                !  on integre entre gmtime et gmtime+radpas
                zdtime=phys_tstep*REAL(radpas) ! pas de temps du rayonnement (s)
                CALL zenang(zlongi,jH_cur,0.0,zdtime, &
                     latitude_deg,longitude_deg,rmu0,fract)
                zrmu0 = rmu0
                swradcorr = 1.0
                ! Calcul du flag jour-nuit
                ! JrNt = 0.0
                !WHERE (fract.GT.0.0) JrNt = 1.0
           CASE(2)
                !  Avec cycle diurne sans application des poids
                !  On integre entre gmtime-pdtphys et gmtime+pdtphys*(radpas-1)
                !  Comme cette routine est appele a tous les pas de temps de
                !  la physique meme si le rayonnement n'est pas appele je
                !  remonte en arriere les radpas-1 pas de temps
                !  suivant. Petite ruse avec MOD pour prendre en compte le
                !  premier pas de temps de la physique pendant lequel
                !  itaprad=0
                zdtime1=phys_tstep*REAL(-MOD(itaprad,radpas)-1)
                zdtime2=phys_tstep*REAL(radpas-MOD(itaprad,radpas)-1)
                CALL zenang(zlongi,jH_cur,zdtime1,zdtime2, &
                        latitude_deg,longitude_deg,rmu0,fract)
                !
                ! Calcul des poids
                !
                zdtime1=-phys_tstep !--on corrige le rayonnement pour representer le
                zdtime2=0.0    !--pas de temps de la physique qui se termine
                CALL zenang(zlongi,jH_cur,zdtime1,zdtime2, &
                        latitude_deg,longitude_deg,zrmu0,zfract)
                swradcorr = 0.0
                WHERE (rmu0.GE.1.e-10 .OR. fract.GE.1.e-10) &
                        swradcorr=zfract/fract*zrmu0/rmu0
                        ! Calcul du flag jour-nuit
                        !JrNt = 0.0
                        !WHERE (zfract.GT.0.0) JrNt = 1.0
          END SELECT
        ENDIF
        !PRINT*,'fract = ',fract
        !sza_o = ACOS (rmu0) *180./pi

        ! Calcul de l'ensoleillement :
        IF (MOD(itaprad,radpas).EQ.0) THEN
           namelist_ecrad_file='namelist_ecrad'
           !
           ! Q0 gestion des pas de temps :
           ! itap, itaprad ?????
           ! 
           ! Q1 champs lies a la surface ? sorties appel au modele de surface
           ! zxtsol temp au sol : dvrait etre recuperee de la partie surface
           ! albsol_dir, albsol_dif aussi
           albsol_dir = 1. ! AI attention
           albsol_dif = 1. ! AI attention
           zxtsol(:) = t(:,1) ! AI attention
           !
           ! Q2 ozone ? lmdz ok a completer
           ! wo : pour l'instant appel ozonecm
           !
           ! Q3 orbite ? lmdz ok a comleter
           ! dist, rmu0 et fract
           ! pour l'instant cas solarlong0<-999.
           ! comment introduire les pas de temps rayonnement
           !  
           ! Q3 proprietes optiques des nuages ? call cloud_optics_prop.F90 ? non necessaires pour ecrad
           ! epaisseur optique, emissivite : cldfrarad, cldemirad, cldtaurad , non necessaire pour ecrad
           ! rayons effectifs : ref_liq, ref_ice ?????
           ! cloud liq/Ice water mixing contenent (kg/kg) ????
           ! les rayons effectifs seront calcules en utilisant la routine de l'IFS
           ! en fonct des contenus en eau, .....
           cldemirad = 1.
           cldtaurad = 5.
           cldtaupirad = 5.
           ref_liq = 1.
           ref_ice = 1.
           ref_liq_pi = 1.
           ref_ice_pi = 1.
           !
           ! Q4 aerosols ?
           ! simplifier l'ecriture : mettre tous les flags et parametres dans aero_mod.F90 ?
           ! ok_ade.OR.flag_aerosol_strat.GT.0, ok_aie,  ok_volcan, flag_volc_surfstrat, &
           ! flag_aerosol, flag_aerosol_strat, flag_aer_feedback, &
           ! tau_aero, piz_aero, cg_aero, &
           ! tau_aero_sw_rrtm, piz_aero_sw_rrtm, cg_aero_sw_rrtm, &
           ! tau_aero_lw_rrtm, &
           !
           ! Q5 thermo : zqsat
           ! zqsat specific humidity at liq saturation => 
           ! calcul de Fred ou a  calculer dans ecrad enf fct de t et p
           ! flwc, fiwc : qx( 2) et qx(  3) 

           PRINT*,'Entree radlwsw a itap : ', itap
           CALL call_radiation &
                  (debut, dist, rmu0, fract,  &
                & paprs, pplay,zxtsol,SFRWL,albsol_dir, albsol_dif,  &
                & t,qx(:,:,1),wo, ZCLDFR, cldemirad, cldtaurad, &
                & tau_aero, piz_aero, cg_aero, &
                & tau_aero_sw_rrtm, piz_aero_sw_rrtm, cg_aero_sw_rrtm, &
                & tau_aero_lw_rrtm, &
                & cldtaupirad, m_allaer, &
                & zqsat, qx(:,:,2), qx(:,:,2), &
                & ref_liq, ref_ice, ref_liq_pi, ref_ice_pi, &
                & namelist_ecrad_file, &
                & heat,heat0,cool,cool0,albpla, heat_volc,cool_volc, &
                & topsw,toplw,solsw,solswfdiff,sollw, &
                & sollwdown, topsw0,toplw0,solsw0,sollw0, &
                & lwdnc0, lwdn0, lwdn, lwupc0, lwup0, lwup,  &
                & lwtoa0b, lwtoab , &
                & swdnc0, swdn0, swdn, swupc0, swup0, swup, &
                & topswad_aero, solswad_aero, topswai_aero, solswai_aero, &
                & topswad0_aero, solswad0_aero, topsw_aero, topsw0_aero, &
                & solsw_aero, solsw0_aero, topswcf_aero, solswcf_aero, &
                & toplwad_aero, sollwad_aero, toplwai_aero, sollwai_aero, &
                & toplwad0_aero, sollwad0_aero,&
                & ZLWFT0_i, ZFLDN0, ZFLUP0, ZSWFT0_i, ZFSDN0, ZFSUP0, &
                & ZFLUX_DIR, ZFLUX_DIR_CLEAR, ZFLUX_DIR_SUN, &
                & cloud_cover_sw)

           itaprad = 0
        ENDIF ! IF (MOD(itaprad,radpas).EQ.0)    
        itaprad = itaprad + 1

        IF (iflag_radia.eq.0) THEN
           heat=0.
           cool=0.
           sollw=0.   ! MPL 01032011
           solsw=0.
           radsol=0.
           swup=0.    ! MPL 27102011 pour les fichiers AMMA_profiles et AMMA_scalars
           swup0=0.
           lwup=0.
           lwup0=0.
           lwdn=0.
           lwdn0=0.
        ENDIF
        ! Ajouter la tendance des rayonnements (tous les pas)
        ! avec une correction pour le cycle diurne dans le SW
        !
        DO k=1, klev
        ! AI attention 
        !d_t_swr(:,k)=swradcorr(:)*heat(:,k)*phys_tstep/RDAY
        !d_t_sw0(:,k)=swradcorr(:)*heat0(:,k)*phys_tstep/RDAY
        !d_t_lwr(:,k)=-cool(:,k)*phys_tstep/RDAY
        !d_t_lw0(:,k)=-cool0(:,k)*phys_tstep/RDAY
           d_t_swr(:,k)=swradcorr(:)*heat(:,k)/RDAY
           d_t_sw0(:,k)=swradcorr(:)*heat0(:,k)/RDAY
           d_t_lwr(:,k)=-cool(:,k)/RDAY
           d_t_lw0(:,k)=-cool0(:,k)/RDAY
        ENDDO

        !t_seri(:,:) = t_seri(:,:) + d_t_swr(:,:) + d_t_lwr(:,:)
        d_t = d_t + d_t_swr + d_t_lwr

        !CALL writefield_phy('d_t_swr',d_t_swr,klev)
        !CALL writefield_phy('d_t_lwr',d_t_lwr,klev)
        !CALL writefield_phy('temp',t,klev)

      ENDIF   ! iflag_radia a nettoyer
      !
      ! AI passer peut-etre la surface avant radiation
      !    albedo, temp surf, ... pour le rayonnement
      !
      !ECRAD can be plugged here and can use the adjusted values for temperature and hydrometeors
      !Cloud fraction is available in the ZCLDFR variable

!
! ==============================================================================================
! Surface
!------------------------------------------------------------
      !
      ! A. Idelkadi 
      ! compute tendencies to return to the dynamics:
      ! "friction" on the first layer
      IF (debut) THEN
        tsrf(1:nlon)=t(1:nlon,1)
        do isoil=1,nsoil
           tsoil(1:nlon,isoil)=tsrf(1:nlon)
           !tsoil(1:nlon,isoil)=302.
        enddo
        rpi=2.*asin(1.)
        timesec=5*3600
        stefan=5.67e-8
        fsw0=900.
        call getin_p('iecri',iecri)
        call getin_p('iflag_surface',iflag_surface)
      ENDIF

      timesec=timesec+pdtphys

      ! Computation of qsats. Done even it is not used for ouptuts
      ! A. Idelkadi
      qsats(:) = EXP( XALPW - XBETAW / t(:,1) - XGAMW * ALOG( t(:,1) ) )
      qsats(:) = XRD / XRV * qsats(:) / ( pplay(:,1) - qsats(:) )
      !qsats(:) = zqsat(:,1)

      IF (iflag_surface == 0 ) THEN
          ! Latent and sensible heat fluxes imposed in 1D (exemple Arm cumulus)
          ! radiative fluxes set to 0.
          ! capcal and fluxgrd=0. replace generic_soil and allow to compute a surface temperature
          ! A case close to Rico can be obtained by imposing fsens=-5E-3 et flat=-6E-5
          capcal(:)=1.e6
          fluxgrd(:)=0.
          sensible(:)=-fsens
          latent(:)=-flat
          swdnsrf(:)=0.
          lwdnsrf(:)=0.
          fluxt(1:klon)=sensible/PHYEX%CST%XCPD
          fluxq(1:klon)=latent(:)/PHYEX%CST%XLVTT
          fluxu(1:klon)=-0.01*u(1:klon,1)
          fluxv(1:klon)=-0.01*v(1:klon,1)
      ELSE
      ! A. Idelkadi
         IF (iflag_radia.eq.0) THEN    
            ! Using the generic soil model
            IF ( nlon == 1 ) THEN
                ! Imposed radiative fluxes to simulate a 1D case close to Arm Cumulus
                swdnsrf(:)=fsw0*max(cos(2*rpi*(timesec-43200.)/(1.15*86400.)),0.)
                lwdnsrf(:)=400.
            ELSE
                ! Model forced in latitude by an idelaized radiative flux that depends on latitude
                ! only. Used for 3D tests without radiation.
                swdnsrf(:)=600*cos(latitude(:))+100
                lwdnsrf(:)=0.
            ENDIF
         ELSE
                swdnsrf(:)=solsw(:)
                lwdnsrf(:)=sollwdown(:)
         ENDIF     
         ! updating the soil temperature (below the surface) from the surface temperature
         ! en returning the surface conductive flux and a heat capacity (accounting for the
         ! sentivity of the conductive flux to the surface temperature.
         call soil_conduct(klon,nsoil,debut,pdtphys,tsrf,tsoil,capcal,fluxgrd,zout)

         ! tsoil_out = tsoil with the vertical dimension of atmosphere outputs

         ! Computing the turbulent fluxes Fqx = rho w'x' = rho (kappa/ln(z/z0x))^2 [ x_surface - x ]. Upward
         call soil_fluxes(klon,tsrf,qsats,pphi(:,1)/rg,ZRHOD(:,2),u(:,1),v(:,1),t(:,1),qx(:,1,1),fluxu,fluxv,fluxt,fluxq)
         sensible(:)=fluxt(:)*PHYEX%CST%XCPD
         latent(:)=fluxq(:)*PHYEX%CST%XLVTT
   ENDIF
   tsoil_out=0. ; tsoil_out(:,1:min(nsoil,nlev))=tsoil(:,1:min(nsoil,nlev))

   ! Computing the surface temperature
   lwupsrf=stefan*tsrf(:)**4
   tsrf(:)=tsrf(:)+pdtphys*(-latent-sensible+lwdnsrf+swdnsrf-lwupsrf+fluxgrd)/capcal(:)

   PSFTH(:) = fluxt(:)
   PSFRV(:) = fluxq(:)
   PSFU(:) = fluxu(:)
   PSFV(:) = fluxv(:)
   PSFSV(:,:) = 0.

!
! ==========================================================================================================
! Shallow convection
!------------------------------------------------------------
!
   CALL SHALLOW_MF(D, PHYEX%CST, PHYEX%NEBN, PHYEX%PARAM_MFSHALLN, PHYEX%TURBN, PHYEX%CSTURB, PHYEX%RAIN_ICE_PARAMN, &
                &KRR=KRR, KRRL=KRRL, KRRI=KRRI, KSV=KSV,                                                        &
                &ONOMIXLG=MISC_LMDZ%ONOMIXLG,KSV_LGBEG=MISC_LMDZ%KSV_LGBEG,KSV_LGEND=MISC_LMDZ%KSV_LGEND,    &
                &PTSTEP=pdtphys,                                                                                &
                &PDZZ=zdzm(:,:),PZZ=zz_mass(:,:),                                                               &
                &PRHODJ=PRHODJ(:,:),PRHODREF=ZRHOD(:,:),                                                        &
                &PPABSM=ZPABST(:,:),PEXNM=ZEXN(:,:),                                                            &
                &PSFTH=PSFTH(:),PSFRV=PSFRV(:),                                                                 &
                &PTHM=ZTHETA(:,:),PRM=ZRX(:,:,:),PUM=ZUT(:,:),PVM=ZVT(:,:),                                     &
                &PTKEM=PTKEM(:,:),PSVM=ZSVT(:,:,:),                                                             &
                &PDUDT_MF=PDUDT_MF(:,:),PDVDT_MF=PDVDT_MF(:,:), PDTKEDT_MF=PDTKEDT_MF(:,:),                     &
                &PDTHLDT_MF=PDTHLDT_MF(:,:),PDRTDT_MF=PDRTDT_MF(:,:),PDSVDT_MF=PDSVDT_MF(:,:,:),                &
                &PSIGMF=PSIGMF(:,:),PRC_MF=PRC_MF(:,:),PRI_MF=PRI_MF(:,:),PCF_MF=PCF_MF(:,:),                   &
                &PHLC_HRC=PHLC_HRC_MF(:,:), PHLC_HCF=PHLC_HCF_MF(:,:), PHLI_HRI=PHLI_HRI_MF(:,:), PHLI_HCF=PHLI_HCF_MF(:,:), &
                &PWEIGHT_MF_CLOUD=PWEIGHT_MF_CLOUD(:,:), PFLXZTHVMF=ZFLXZTHVMF(:,:), &
                &PFLXZTHMF=ZFLXZTHMF(:,:),PFLXZRMF=ZFLXZRMF(:,:),PFLXZUMF=ZFLXZUMF(:,:),PFLXZVMF=ZFLXZVMF(:,:), &
                &PFLXZTKEMF=ZFLXZTKEMF(:,:), &
                &PTHL_UP=PTHL_UP(:,:),PRT_UP=PRT_UP(:,:),PRV_UP=PRV_UP(:,:),                                    &
                &PRC_UP=PRC_UP(:,:),PRI_UP=PRI_UP(:,:),                                                         &
                &PU_UP=PU_UP(:,:), PV_UP=PV_UP(:,:), PTKE_UP=PTKE_UP(:,:), PTHV_UP=PTHV_UP(:,:), PW_UP=PW_UP(:,:),                    &
                &PFRAC_UP=PFRAC_UP(:,:),PEMF=PEMF(:,:),PDETR=PDETR(:,:),PENTR=PENTR(:,:),                       &
                &KKLCL=IKLCL(:),KKETL=IKETL(:),KKCTL=IKCTL(:),PDX=1000.0,PDY=1000.0,KBUDGETS=MISC_LMDZ%NBUDGET )

   ! Add tendencies of shallow to total physics tendency
   d_u(:,1:klev) = d_u(:,1:klev) + PDUDT_MF(:,2:klev+1)
   d_v(:,1:klev) = d_v(:,1:klev) + PDVDT_MF(:,2:klev+1) 
   ZRXS(:,:,1)=ZRXS(:,:,1)+PDRTDT_MF(:,:)
   ZTHETAS(:,:)=ZTHETAS(:,:)+PDTHLDT_MF(:,:)
   ZTKES(:,:) = ZTKES(:,:)+PDTKEDT_MF(:,:)
   ! TODO add SV tendencies

!
!==================================================================================================================
! Turbulence
!------------------------------------------------------------------------------------------------------------------
   ! out tendencies
   ZRUS(:,:) = 0.
   ZRVS(:,:) = 0.
   ZRWS(:,:) = 0.
   ZRSVS(:,:,:) = 0.
   ZRTKES(:,:) =  ZTKES(:,:) * PRHODJ(:,:)
   DO JRR=1, KRR
      ZRRS(:,:,JRR) = ZRXS(:,:,JRR) * PRHODJ(:,:)
   ENDDO
   ZRTHS(:,:) = ZTHETAS(:,:) * PRHODJ(:,:)
   CALL TURB(PHYEX%CST, PHYEX%CSTURB, MISC_LMDZ%TBUCONF, PHYEX%TURBN, PHYEX%NEBN, D, MISC_LMDZ%TLES,               &
        & KRR, KRRL, KRRI, MISC_LMDZ%HLBCX, MISC_LMDZ%HLBCY, MISC_LMDZ%KGRADIENTSLEO, &
        & MISC_LMDZ%KGRADIENTSGOG, MISC_LMDZ%KHALO,                &
        & PHYEX%TURBN%NTURBSPLIT, PHYEX%TURBN%LCLOUDMODIFLM, KSV, MISC_LMDZ%KSV_LGBEG, MISC_LMDZ%KSV_LGEND,          &
        & MISC_LMDZ%KSV_LIMA_NR, MISC_LMDZ%KSV_LIMA_NS, MISC_LMDZ%KSV_LIMA_NG, MISC_LMDZ%KSV_LIMA_NH,              &
        & MISC_LMDZ%O2D, MISC_LMDZ%ONOMIXLG, MISC_LMDZ%OFLAT, MISC_LMDZ%OCOUPLES,                                  &
        & MISC_LMDZ%OBLOWSNOW,MISC_LMDZ%OIBM,                                                                        &
        & MISC_LMDZ%OFLYER, MISC_LMDZ%COMPUTE_SRC, MISC_LMDZ%PRSNOW,                                                &
        & MISC_LMDZ%OOCEAN, MISC_LMDZ%ODEEPOC, MISC_LMDZ%ODIAG_IN_RUN,                                              &
        & PHYEX%TURBN%CTURBLEN_CLOUD, MISC_LMDZ%CMICRO, MISC_LMDZ%CELEC,                                             &
        & pdtphys,MISC_LMDZ%ZTFILE,                                                                 &
        & ZDXX(:,:),ZDYY(:,:),zdzm(:,:),                                                             &
        & ZDZX(:,:),ZDZY(:,:),zz_flux(:,:),                                                          &
        & ZDIRCOSXW(:),ZDIRCOSYW(:),ZDIRCOSZW(:),ZCOSSLOPE(:),ZSINSLOPE(:),                          &
        & PRHODJ(:,:),PTHVREF(:,:), PHGRADLEO(:,:,:), PHGRADGOG(:,:,:), zs(:),                       &
        & ZSEA_UCU(:),ZSEA_VCU(:),                                                                   &
        & PSFTH(:),PSFRV(:),PSFSV(:,:),PSFU(:),PSFV(:),                                              &
        & ZPABST(:,:),ZUT(:,:),ZVT(:,:),PWT(:,:),PTKEM(:,:),ZSVT(:,:,:),ZSRC(:,:),                   &
        & PLENGTHM(:,:),PLENGTHH(:,:),MFMOIST(:,:),                                                  &
        & ZBL_DEPTH(:),ZSBL_DEPTH(:),                                                                &
        & ZCEI(:,:), PHYEX%TURBN%XCEI_MIN, PHYEX%TURBN%XCEI_MAX, PHYEX%TURBN%XCOEF_AMPL_SAT,         &
        & ZTHETA(:,:),ZRX(:,:,:),                                                                    &
        & ZRUS(:,:),ZRVS(:,:),ZRWS(:,:),ZRTHS(:,:),ZRRS(:,:,:),ZRSVS(:,:,:),ZRTKES(:,:),             &
        & PSIGS(:,:),                                                                                &
        & ZFLXZTHVMF(:,:),ZFLXZUMF(:,:), ZFLXZVMF,                                                   &
        & ZWTH(:,:),ZWRC(:,:),ZWSV(:,:,:),ZDP(:,:),ZTP(:,:),ZTDIFF(:,:),ZTDISS(:,:),                 &
        & MISC_LMDZ%YLBUDGET, MISC_LMDZ%NBUDGET                                                    )
   DO JRR=1, KRR
      ZRXS(:,:,JRR) = ZRRS(:,:,JRR) / PRHODJ(:,:)
   ENDDO
   ZTHETAS(:,:) = ZRTHS(:,:) / PRHODJ(:,:)
   ZTKES(:,:) = ZRTKES(:,:) / PRHODJ(:,:)
   ! Add tendencies of turb to total physics tendency
   d_u(:,1:klev) = d_u(:,1:klev) + ZRUS(:,2:klev+1)/PRHODJ(:,2:klev+1)
   d_v(:,1:klev) = d_v(:,1:klev) + ZRVS(:,2:klev+1)/PRHODJ(:,2:klev+1)
   IF(PHYEX%PARAM_MFSHALLN%CMF_CLOUD=='STAT') THEN
      PSIGS(:,:)=SQRT(PSIGS(:,:)**2 + PSIGMF(:,:)**2)
   ENDIF

! ==================================================================================================================
! Microphysics
!------------------------------------------------------------
   ZSEA=1.
   ZTOWN=0.
   ZCIT=0.
   ZTHVREFZIKB=0
   CALL RAIN_ICE (D, PHYEX%CST, PHYEX%PARAM_ICEN, PHYEX%RAIN_ICE_PARAMN, PHYEX%RAIN_ICE_DESCRN,                      &
               PHYEX%ELEC_PARAM, PHYEX%ELEC_DESCR, MISC_LMDZ%TBUCONF,                                            &
               MISC_LMDZ%OELEC, MISC_LMDZ%OSEDIM_BEARD,                                                         &
               ZTHVREFZIKB,                                                                                       &
               pdtphys, KRR, ZEXN,                                                                                &
               zdzf, PRHODJ, ZRHOD, ZEXN, ZPABST, ZCIT, ZCLDFR,                                                   &
               ZICLDFR, ZSSIO, ZSSIU, ZIFR,                                                                       &
               ZHLC_HRC, ZHLC_HCF, ZHLI_HRI, ZHLI_HCF,                                                            &
               ztheta, ZRX, zthetas, ZRXS, &
               ZINPRC, ZINPRR, ZEVAP3D,                                                                           &
               ZINPRS, ZINPRG, ZINDEP, ZRAINFR, PSIGS,                                                            &
               MISC_LMDZ%YLBUDGET, MISC_LMDZ%NBUDGET,                                                           &
               ZSEA, ZTOWN                                                                                        )

! ====================================================================================================================
! Tendencies and time evolution (values for next time step)
!------------------------------------------------------------
   ! Tendencies, mixing ratio -> specific
   d_qx(:,1:klev,1)=d_qx(:,1:klev,1) + (ZRXS(:,2:klev+1,1)-ZRXS0(:,2:klev+1,1))*ZQDM(:,2:klev+1)
   d_qx(:,1:klev,2)=d_qx(:,1:klev,2) + (ZRXS(:,2:klev+1,2)-ZRXS0(:,2:klev+1,2))*ZQDM(:,2:klev+1)
   d_qx(:,1:klev,3)=d_qx(:,1:klev,3) + (ZRXS(:,2:klev+1,4)-ZRXS0(:,2:klev+1,4))*ZQDM(:,2:klev+1)
   d_qr(:,1:klev)=d_qr(:,1:klev) + (ZRXS(:,2:klev+1,3)-ZRXS0(:,2:klev+1,3))*ZQDM(:,2:klev+1)
   d_qs(:,1:klev)=d_qs(:,1:klev) + (ZRXS(:,2:klev+1,5)-ZRXS0(:,2:klev+1,5))*ZQDM(:,2:klev+1)
   d_qg(:,1:klev)=d_qg(:,1:klev) + (ZRXS(:,2:klev+1,6)-ZRXS0(:,2:klev+1,6))*ZQDM(:,2:klev+1)
   ! Tendency, theta -> T
   d_t(:,1:klev)=d_t(:,1:klev) + (zthetas(:,2:klev+1)-zthetas0(:,2:klev+1))*ZEXN(:,2:klev+1)
   ! TKE
   d_tke(:,1:klev)=d_tke(:,1:klev) + (ZTKES(:,2:klev+1) - ZTKES0(:,2:klev+1))

   !Time evolution
   ZQR(:,:)=ZQR(:,:)+d_qr(:,:)*pdtphys
   ZQS(:,:)=ZQS(:,:)+d_qs(:,:)*pdtphys
   ZQG(:,:)=ZQG(:,:)+d_qg(:,:)*pdtphys
   PTKEM(:,2:klev+1)=PTKEM(:,2:klev+1)+d_tke(:,:)*pdtphys

!
! ===================================================================================================================
! Entrees sorties
!------------------------------------------------------------

   itap=itap+1
   !zjulian=zjulian+pdtphys/86400.
   !PRINT*,'avant mod(itap,iecri) , zjulian itap',zjulian,itap,iecri
   IF ( mod(itap,iecri) == 1 .or. iecri == 1 ) THEN
        PRINT*,'ECRITURE ITAP= AAA',itap
        call output_physiqex(debut,zjulian,pdtphys_,presnivs,paprs,u,v,t,qx,ZCLDFR,ZQR,ZQS,ZQG,PTKEM,ZTHETA)
  call iophys_ecrit('dqmf',klev,'MF Hum tend',' ',PDRTDT_MF(:,2:klev+1))
  call iophys_ecrit('dtmf',klev,'MF temp tend',' ',PDTHLDT_MF(:,2:klev+1))
        call iophys_ecrit('tsrf',       1,'tsrf',       ' ',tsrf)
        call iophys_ecrit('capcal',     1,'capcal',     ' ',capcal)
        call iophys_ecrit('swdnsrf',    1,'swdnsrf',    ' ',swdnsrf)
        call iophys_ecrit('lwdnsrf',    1,'lwdnsrf',    ' ',lwdnsrf)
        call iophys_ecrit('lwupsrf',    1,'lwupsrf',    ' ',lwupsrf)
        call iophys_ecrit('sensible',   1,'sensible',   ' ',sensible)
        call iophys_ecrit('latent',     1,'latent',     ' ',latent)
        call iophys_ecrit('qsats',      1,'qsats',      ' ',qsats)
        call iophys_ecrit('tsoil',   nlev,'tsoi',      ' ',tsoil_out)
        !AI
        IF (iflag_radia.ge.1) THEN
                call iophys_ecrit('lwdn0',       1,'lwdn0',       ' ',lwdn0)
                call iophys_ecrit('lwdn',       1,'lwdn',       ' ',lwdn)
                call iophys_ecrit('lwup0',       1,'lwup0',       ' ',lwup0)
                call iophys_ecrit('lwup',       1,'lwup',       ' ',lwup)
                call iophys_ecrit('swdn0',       1,'swdn0',       ' ',swdn0)
                call iophys_ecrit('swdn',       1,'swdn',       ' ',swdn)
                call iophys_ecrit('swup0',       1,'swup0',       ' ',swup0)
                call iophys_ecrit('swup',       1,'swup',       ' ',swup)
        ENDIF
    ENDIF

    ! IF lastcall, THEN it is time to write "restartphy.nc" file
    IF (lafin) THEN
        call phyredem("restartphy.nc")
    ENDIF


  end subroutine physiqex
  !
  SUBROUTINE VERTICAL_EXTEND(PX,KLEV)

  ! fill extra vetical levels to fit MNH interface

        REAL, DIMENSION(:,:),   INTENT(INOUT)   :: PX
        INTEGER, INTENT(IN) :: KLEV
        PX(:,1     )= PX(:,2)
        PX(:,KLEV+2)= PX(:,KLEV+1)
  END SUBROUTINE VERTICAL_EXTEND

END MODULE physiqex_mod
