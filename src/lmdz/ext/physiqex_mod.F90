MODULE physiqex_mod

IMPLICIT NONE

CONTAINS


!TODO list pour un branchement plus propre
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
! * Certains allocatable sont en klev, d'autres en klev+2. Il est possible de changer ceci mais il faut gérer
!   les recopies pour que les params voient des tableaux en klev+2 (en effectuant des recopies des
!   niveaux extrêmes comme fait pour le vent par exemple)

!TODO à faire avant d'historiser dans LMDZ:
! * vérifier si rhodj ne devrait pas être calculé avec un rho humide ici (arome?)
! * brancher le sigma du schéma STAT
! * utiliser les variables ajustées














      SUBROUTINE physiqex (nlon,nlev, &
     &            debut,lafin,pdtphys, &
     &            paprs,pplay,pphi,pphis,presnivs, &
     &            u,v,rot,t,qx, &
     &            flxmass_w, &
     &            d_u, d_v, d_t, d_qx, d_ps)

      USE dimphy, only : klon,klev
      USE infotrac_phy, only : nqtot
      USE geometry_mod, only : latitude, cell_area
!      USE comcstphy, only : rg
      USE ioipsl, only : ymds2ju
      USE phys_state_var_mod, only : phys_state_var_init
      USE phyetat0_mod, only: phyetat0
      USE output_physiqex_mod, ONLY: output_physiqex
      ! PHYEX internal modules
      USE MODE_INIT_PHYEX, ONLY: INIT_PHYEX, FILL_DIMPHYEX
      USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
      USE MODD_PHYEX,      ONLY: PHYEX_t
      USE MODI_INI_PHYEX,  ONLY: INI_PHYEX
      USE MODI_ICE_ADJUST, ONLY: ICE_ADJUST
      USE MODI_RAIN_ICE, ONLY: RAIN_ICE
      USE MODI_TURB
      USE MODI_SHALLOW_MF
      IMPLICIT none
!
! Routine argument:
!
      integer,intent(in) :: nlon ! number of atmospheric colums
      integer,intent(in) :: nlev ! number of vertical levels (should be =klev)
      logical,intent(in) :: debut ! signals first call to physics
      logical,intent(in) :: lafin ! signals last call to physics
      real,intent(in) :: pdtphys ! physics time step (s)
      real,intent(in) :: paprs(klon,klev+1) ! interlayer pressure (Pa)
      real,intent(in) :: pplay(klon,klev) ! mid-layer pressure (Pa)
      real,intent(in) :: pphi(klon,klev) ! geopotential at mid-layer
      real,intent(in) :: pphis(klon) ! surface geopotential
      real,intent(in) :: presnivs(klev) ! pseudo-pressure (Pa) of mid-layers
      real,intent(in) :: u(klon,klev) ! eastward zonal wind (m/s)
      real,intent(in) :: v(klon,klev) ! northward meridional wind (m/s)
      real,intent(in) :: rot(klon,klev) ! northward meridional wind (m/s)
      real,intent(in) :: t(klon,klev) ! temperature (K)
      real,intent(in) :: qx(klon,klev,nqtot) ! tracers (.../kg_air)
      real,intent(in) :: flxmass_w(klon,klev) ! vertical mass flux
      real,intent(out) :: d_u(klon,klev) ! physics tendency on u (m/s/s)
      real,intent(out) :: d_v(klon,klev) ! physics tendency on v (m/s/s)
      real,intent(out) :: d_t(klon,klev) ! physics tendency on t (K/s)
      real,intent(out) :: d_qx(klon,klev,nqtot) ! physics tendency on tracers
      real,intent(out) :: d_ps(klon) ! physics tendency on surface pressure
      real :: d_qr(klon, klev), d_qs(klon, klev), d_qg(klon, klev) ! tendency for rain, snow, graupel
      real :: d_tke(klon, klev)

!    include "clesphys.h"
    include "flux_arp.h"

    INTEGER        length
    PARAMETER    ( length = 100 )
    REAL tabcntr0( length       )
    INTEGER, PARAMETER :: longcles=20
    REAL, SAVE :: clesphy0(longcles)
    !$OMP THREADPRIVATE(clesphy0)

! Saved variables
REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: PSIGS !variance of s
REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: PCF_MF, PRC_MF, PRI_MF !shallow convection cloud
REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: ZQR, ZQS, ZQG !rain, snow, graupel specifiq contents
REAL, DIMENSION(:,:), ALLOCATABLE, SAVE ::  PTKEM       ! TKE
TYPE(DIMPHYEX_t),SAVE    :: D
TYPE(PHYEX_t), SAVE      :: PHYEX
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
real :: ZDZMIN
real, dimension(klon, klev+2) :: ZCIT !ice concentration
real, dimension(klon) :: ZINPRC, ZINPRR, ZINPRS, ZINPRG !precipitation flux at ground
real, dimension(klon, klev+2) :: ZEVAP3D !evaporation (diag)
real, dimension(klon, klev+2) :: ZRAINFR
real, dimension(klon) :: ZINDEP !deposition flux, already contained in ZINPRC
real, dimension(klon)         :: ZSEA, ZTOWN !sea and town fractions in the frid cell
! Turbulence variables
REAL, DIMENSION(klon,klev+2)  ::  PDXX,PDYY,PDZX,PDZY ! metric coefficients
REAL, DIMENSION(klon,klev+2)  ::  PZZ       !  physical distance between 2 succesive grid points along the K direction
REAL, DIMENSION(klon)         ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW ! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(klon)         ::  PCOSSLOPE       ! cosinus of the angle between i and the slope vector
REAL, DIMENSION(klon)         ::  PSINSLOPE       ! sinus of the angle   between i and the slope vector
REAL, DIMENSION(klon,klev+2)  ::  PRHODJ   ! dry density * Grid size
REAL, DIMENSION(0,0)          ::  MFMOIST  ! moist mass flux dual scheme
REAL, DIMENSION(0,0,0)        ::  PHGRAD      ! horizontal gradients
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
REAL, DIMENSION(klon,klev+2)  ::  ZDXX,ZDYY,ZDZX,ZDZY,ZZZ
REAL, DIMENSION(klon)         ::  ZDIRCOSXW,ZDIRCOSYW,ZDIRCOSZW,ZCOSSLOPE,ZSINSLOPE
! Shallow variables 
REAL, DIMENSION(klon,klev+2) ::  PDUDT_MF     ! tendency of U   by massflux scheme
REAL, DIMENSION(klon,klev+2) ::  PDVDT_MF     ! tendency of V   by massflux scheme
REAL, DIMENSION(klon,klev+2) ::  PDTHLDT_MF   ! tendency of thl by massflux scheme
REAL, DIMENSION(klon,klev+2) ::  PDRTDT_MF    ! tendency of rt  by massflux scheme
REAL, DIMENSION(klon,klev+2,KSV)::  PDSVDT_MF    ! tendency of Sv  by massflux scheme
REAL, DIMENSION(klon,klev+2) ::  PSIGMF
REAL, DIMENSION(klon,klev+2) ::  ZFLXZTHMF
REAL, DIMENSION(klon,klev+2) ::  ZFLXZRMF
REAL, DIMENSION(klon,klev+2) ::  ZFLXZUMF
REAL, DIMENSION(klon,klev+2) ::  ZFLXZVMF
REAL, DIMENSION(klon,klev+2) ::  PTHL_UP   ! Thl updraft characteristics
REAL, DIMENSION(klon,klev+2) ::  PRT_UP    ! Rt  updraft characteristics
REAL, DIMENSION(klon,klev+2) ::  PRV_UP    ! Vapor updraft characteristics
REAL, DIMENSION(klon,klev+2) ::  PU_UP     ! U wind updraft characteristics
REAL, DIMENSION(klon,klev+2) ::  PV_UP     ! V wind updraft characteristics
REAL, DIMENSION(klon,klev+2) ::  PRC_UP    ! cloud content updraft characteristics
REAL, DIMENSION(klon,klev+2) ::  PRI_UP    ! ice content   updraft characteristics
REAL, DIMENSION(klon,klev+2) ::  PTHV_UP   ! Thv   updraft characteristics
REAL, DIMENSION(klon,klev+2) ::  PW_UP     ! vertical speed updraft characteristics
REAL, DIMENSION(klon,klev+2) ::  PFRAC_UP  ! updraft fraction
REAL, DIMENSION(klon,klev+2) ::  PEMF      ! updraft mass flux
REAL, DIMENSION(klon,klev+2) ::  PDETR     ! updraft detrainment
REAL, DIMENSION(klon,klev+2) ::  PENTR     ! updraft entrainment
INTEGER,DIMENSION(klon) ::IKLCL,IKETL,IKCTL ! level of LCL,ETL and CTL
!
real :: temp_newton(klon,klev)
integer :: k
logical, save :: first=.true.
!$OMP THREADPRIVATE(first)
real,save :: rg=9.81
!$OMP THREADPRIVATE(rg)

! For I/Os
integer :: itau0
real :: zjulian


!------------------------------------------------------------
! Initialisations de la physique au premier pas de temps
!------------------------------------------------------------

print*,'Debut physiqex',debut
! initializations
if (debut) then ! Things to do only for the first call to physics 
  print*,'Debut physiqex IN'
  
  ! load initial conditions for physics (including the grid)
    call phys_state_var_init(1) ! some initializations, required before calling phyetat0
    call phyetat0("startphy.nc", clesphy0, tabcntr0)
  
  ! Initialize outputs:
    itau0=0
    ! compute zjulian for annee0=1979 and month=1 dayref=1 and hour=0.0
    !CALL ymds2ju(annee0, month, dayref, hour, zjulian)
    call ymds2ju(1979, 1, 1, 0.0, zjulian)
  
  ZDZMIN=MINVAL((pphi(:,2:) - pphi(:,1:klev-1))/9.81)
  CALL INIT_PHYEX(pdtphys, ZDZMIN, PHYEX)
  CALL FILL_DIMPHYEX(KLON, KLEV, D)

  !Update default values
  PHYEX%NEBN%LSUBG_COND = .TRUE. 
  PHYEX%PARAM_ICEN%CSUBG_AUCV_RC='PDF'
  !  
  ! Variables saved
  ALLOCATE(PTKEM(klon,klev+2))
  ALLOCATE(PSIGS(klon,klev+2))
  ALLOCATE(PCF_MF(klon,klev+2))
  ALLOCATE(PRC_MF(klon,klev+2))
  ALLOCATE(PRI_MF(klon,klev+2))
  ALLOCATE(ZQR(klon, klev))
  ALLOCATE(ZQS(klon, klev))
  ALLOCATE(ZQG(klon, klev))
  PSIGS=0.
  PCF_MF=0.
  PRC_MF=0.
  PRI_MF=0.
  ZQR=0.
  ZQS=0.
  ZQG=0.
  PTKEM(:,:) = PHYEX%TURBN%XTKEMIN ! TODO: init from TKE at stationnary state
!
#ifndef CPP_IOIPSL_NO_OUTPUT
  ! Initialize IOIPSL output file
#endif
endif ! of if (debut)

!------------------------------------------------------------
! Initialisations a chaque pas de temps
!------------------------------------------------------------

! set all tendencies to zero
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
ZDXX(:,:) = 0.
ZDYY(:,:) = 0.
ZDZX(:,:) = 0.
ZDZY(:,:) = 0.
ZDIRCOSXW(:) = 1.
ZDIRCOSYW(:) = 1.
ZDIRCOSZW(:) = 1.
ZCOSSLOPE(:) = 0.
ZSINSLOPE(:) = 1.
PHGRAD(:,:,:) = 0.
ZBL_DEPTH(:) = 0. ! needed only with LRMC01 key (correction in the surface boundary layer)
ZSBL_DEPTH(:) = 0.
ZCEI(:,:) = 0.  ! needed only if HTURBLEN_CL /= 'NONE' modification of mixing lengh inside clouds
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
!To compute the actual tendecy, we save the initial values of these variables
ZRXS0(:,:,:) = ZRXS(:,:,:)
ZTHETAS0=ZTHETAS
ZTKES0(:,:)=ZTKES(:,:)
!------------------------------------------------------------
! Adjustment
!------------------------------------------------------------
!
ZSRC(:,:) = 0.
ZSIGQSAT=PHYEX%NEBN%VSIGQSAT
CALL ICE_ADJUST (D, PHYEX%CST, PHYEX%RAIN_ICE_PARAMN, PHYEX%NEBN, PHYEX%TURBN, PHYEX%PARAM_ICEN,    &
                &PHYEX%MISC%TBUCONF, KRR,                                                           &
                &'ADJU',                                                                            &
                &pdtphys, ZSIGQSAT,                                                                 &
                &PRHODJ, ZEXN, ZRHOD, PSIGS, .FALSE., zmfconv,                                      &
                &ZPABST, ZZ_MASS,                                                                   &
                &ZEXN, PCF_MF, PRC_MF, PRI_MF,                                                      &
                &ZICLDFR, ZWCLDFR, ZSSIO, ZSSIU, ZIFR,                                              &
                &ZRX(:,:,1), ZRX(:,:,2), ZRXS(:,:,1), ZRXS(:,:,2), ZTHETA, ZTHETAS,                 &
                &PHYEX%MISC%COMPUTE_SRC, ZSRC, ZCLDFR,                                              &
                &ZRX(:,:,3), ZRX(:,:,4), ZRXS(:,:,4), ZRX(:,:,5), ZRX(:,:,6),                       &
                &PHYEX%MISC%YLBUDGET, PHYEX%MISC%NBUDGET,                                           &
                &ZICE_CLD_WGT,                                                                      &
                !&POUT_RV, POUT_RC, POUT_RI, POUT_TH,               &
                &ZHLC_HRC, ZHLC_HCF, ZHLI_HRI, ZHLI_HCF                                             )
!
!------------------------------------------------------------
! Surface
!------------------------------------------------------------
!
! compute tendencies to return to the dynamics:
! "friction" on the first layer
d_u(1:klon,1)=d_u(1:klon,1)-u(1:klon,1)/86400.
d_v(1:klon,1)=d_v(1:klon,1)-v(1:klon,1)/86400.
!
! Flux RICO
PSFTH(:) = 5E-3 ! RICO
PSFRV(:) = 6E-5 ! RICO
! Flux ARMCU
!PSFTH(:) = -fsens/1000.
!PSFRV(:) = -flat/(2.5e6)
!
PSFSV(:,:) = 0.
PSFU(:) = 0.
PSFV(:) = 0.
!
!TODO PSIGMF option STAT
!------------------------------------------------------------
! Shallow convection
!------------------------------------------------------------
!
  CALL SHALLOW_MF(D, PHYEX%CST, PHYEX%NEBN, PHYEX%PARAM_MFSHALLN, PHYEX%TURBN, PHYEX%CSTURB,         &
     &KRR=KRR, KRRL=KRRL, KRRI=KRRI, KSV=KSV,                                                        &
     &ONOMIXLG=PHYEX%MISC%ONOMIXLG,KSV_LGBEG=PHYEX%MISC%KSV_LGBEG,KSV_LGEND=PHYEX%MISC%KSV_LGEND,    &
     &PTSTEP=pdtphys,                                                                                &
     &PDZZ=zdzm(:,:),PZZ=zz_mass(:,:),                                                               &
     &PRHODJ=PRHODJ(:,:),PRHODREF=ZRHOD(:,:),                                                        &
     &PPABSM=ZPABST(:,:),PEXNM=ZEXN(:,:),                                                            &
     &PSFTH=PSFTH(:),PSFRV=PSFRV(:),                                                                 &
     &PTHM=ZTHETA(:,:),PRM=ZRX(:,:,:),PUM=ZUT(:,:),PVM=ZVT(:,:),                                     &
     &PTKEM=PTKEM(:,:),PSVM=ZSVT(:,:,:),                                                             &
     &PDUDT_MF=PDUDT_MF(:,:),PDVDT_MF=PDVDT_MF(:,:),                                                 &
     &PDTHLDT_MF=PDTHLDT_MF(:,:),PDRTDT_MF=PDRTDT_MF(:,:),PDSVDT_MF=PDSVDT_MF(:,:,:),                &
     &PSIGMF=PSIGMF(:,:),PRC_MF=PRC_MF(:,:),PRI_MF=PRI_MF(:,:),PCF_MF=PCF_MF(:,:),                   &
     &PFLXZTHVMF=ZFLXZTHVMF(:,:),                                                                    &
     &PFLXZTHMF=ZFLXZTHMF(:,:),PFLXZRMF=ZFLXZRMF(:,:),PFLXZUMF=ZFLXZUMF(:,:),PFLXZVMF=ZFLXZVMF(:,:), &
     &PTHL_UP=PTHL_UP(:,:),PRT_UP=PRT_UP(:,:),PRV_UP=PRV_UP(:,:),                                    &
     &PRC_UP=PRC_UP(:,:),PRI_UP=PRI_UP(:,:),                                                         &
     &PU_UP=PU_UP(:,:), PV_UP=PV_UP(:,:), PTHV_UP=PTHV_UP(:,:), PW_UP=PW_UP(:,:),                    &
     &PFRAC_UP=PFRAC_UP(:,:),PEMF=PEMF(:,:),PDETR=PDETR(:,:),PENTR=PENTR(:,:),                       &
     &KKLCL=IKLCL(:),KKETL=IKETL(:),KKCTL=IKCTL(:),PDX=1000.0,PDY=1000.0,KBUDGETS=PHYEX%MISC%NBUDGET )

! Add tendencies of shallow to total physics tendency
d_u(:,1:klev) = d_u(:,1:klev) + PDUDT_MF(:,2:klev+1)
d_v(:,1:klev) = d_v(:,1:klev) + PDVDT_MF(:,2:klev+1) 
ZRXS(:,:,1)=ZRXS(:,:,1)+PDRTDT_MF(:,:)
ZTHETAS(:,:)=ZTHETAS(:,:)+PDTHLDT_MF(:,:)
! TODO add SV tendencies
!
!------------------------------------------------------------
! Turbulence
!------------------------------------------------------------
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
CALL TURB(PHYEX%CST, PHYEX%CSTURB, PHYEX%MISC%TBUCONF, PHYEX%TURBN, PHYEX%NEBN, D, PHYEX%MISC%TLES,               &
   & PHYEX%MISC%KMI, KRR, KRRL, KRRI, PHYEX%MISC%HLBCX, PHYEX%MISC%HLBCY, PHYEX%MISC%KGRADIENTS, PHYEX%MISC%KHALO,&
   & PHYEX%MISC%KSPLIT,PHYEX%MISC%KMI, KSV, PHYEX%MISC%KSV_LGBEG, PHYEX%MISC%KSV_LGEND,                           &
   & PHYEX%MISC%CPROGRAM,                                                                                         &
   & PHYEX%MISC%KSV_LIMA_NR, PHYEX%MISC%KSV_LIMA_NS, PHYEX%MISC%KSV_LIMA_NG, PHYEX%MISC%KSV_LIMA_NH,              &
   & PHYEX%MISC%O2D, PHYEX%MISC%ONOMIXLG, PHYEX%MISC%OFLAT, PHYEX%MISC%OCOUPLES,                                  &
   &  PHYEX%MISC%OBLOWSNOW,PHYEX%MISC%OIBM,                                                                       &
   & PHYEX%MISC%OFLYER, PHYEX%MISC%COMPUTE_SRC, PHYEX%MISC%PRSNOW,                                                &
   & PHYEX%MISC%OOCEAN, PHYEX%MISC%ODEEPOC, PHYEX%MISC%ODIAG_IN_RUN,                                              &
   & PHYEX%MISC%HTURBLEN_CL,PHYEX%MISC%CMICRO,                                                                    &
   & pdtphys,PHYEX%MISC%ZTFILE,                                                                 &
   & ZDXX(:,:),ZDYY(:,:),zdzm(:,:),                                                             &
   & ZDZX(:,:),ZDZY(:,:),zz_flux(:,:),                                                          &
   & ZDIRCOSXW(:),ZDIRCOSYW(:),ZDIRCOSZW(:),ZCOSSLOPE(:),ZSINSLOPE(:),                          &
   & PRHODJ(:,:),PTHVREF(:,:), PHGRAD(:,:,:), zs(:),                                            &
   & PSFTH(:),PSFRV(:),PSFSV(:,:),PSFU(:),PSFV(:),                                              &
   & ZPABST(:,:),ZUT(:,:),ZVT(:,:),PWT(:,:),PTKEM(:,:),ZSVT(:,:,:),ZSRC(:,:),                   &
   & PLENGTHM(:,:),PLENGTHH(:,:),MFMOIST(:,:),                                                  &
   & ZBL_DEPTH(:),ZSBL_DEPTH(:),                                                                &
   & ZCEI(:,:),PHYEX%MISC%ZCEI_MIN,PHYEX%MISC%ZCEI_MAX,PHYEX%MISC%ZCOEF_AMPL_SAT,               &
   & ZTHETA(:,:),ZRX(:,:,:),                                                                    &
   & ZRUS(:,:),ZRVS(:,:),ZRWS(:,:),ZRTHS(:,:),ZRRS(:,:,:),ZRSVS(:,:,:),ZRTKES(:,:),             &
   & PSIGS(:,:),                                                                                &
   & ZFLXZTHVMF(:,:),ZWTH(:,:),ZWRC(:,:),ZWSV(:,:,:),ZDP(:,:),ZTP(:,:),ZTDIFF(:,:),ZTDISS(:,:), &
   & PHYEX%MISC%YLBUDGET, PHYEX%MISC%NBUDGET                                                    )
DO JRR=1, KRR
  ZRXS(:,:,JRR) = ZRRS(:,:,JRR) / PRHODJ(:,:)
ENDDO
ZTHETAS(:,:) = ZRTHS(:,:) / PRHODJ(:,:)
ZTKES(:,:) = ZRTKES(:,:) / PRHODJ(:,:)
! Add tendencies of turb to total physics tendency
d_u(:,1:klev) = d_u(:,1:klev) + ZRUS(:,2:klev+1)/PRHODJ(:,2:klev+1)
d_v(:,1:klev) = d_v(:,1:klev) + ZRVS(:,2:klev+1)/PRHODJ(:,2:klev+1)
!------------------------------------------------------------
! Microphysics
!------------------------------------------------------------
ZSEA=1.
ZTOWN=0.
ZCIT=0.
CALL RAIN_ICE (D, PHYEX%CST, PHYEX%PARAM_ICEN, PHYEX%RAIN_ICE_PARAMN, PHYEX%RAIN_ICE_DESCRN, PHYEX%MISC%TBUCONF,  &
               pdtphys, KRR, ZEXN,                                                                                &
               zdzf, PRHODJ, ZRHOD, ZEXN, ZPABST, ZCIT, ZCLDFR,                                                   &
               ZHLC_HRC, ZHLC_HCF, ZHLI_HRI, ZHLI_HCF,                                                            &
               ztheta, ZRX(:,:,1), ZRX(:,:,2), ZRX(:,:,3), ZRX(:,:,4), ZRX(:,:,5),                                &
               ZRX(:,:,6), zthetas, ZRXS(:,:,1), ZRXS(:,:,2), ZRXS(:,:,3), ZRXS(:,:,4), ZRXS(:,:,5), ZRXS(:,:,6), &
               ZINPRC, ZINPRR, ZEVAP3D,                                                                           &
               ZINPRS, ZINPRG, ZINDEP, ZRAINFR, PSIGS,                                                            &
               PHYEX%MISC%YLBUDGET, PHYEX%MISC%NBUDGET,                                                           &
               ZSEA, ZTOWN                                                                                        )

!------------------------------------------------------------
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
!------------------------------------------------------------
! Entrees sorties
!------------------------------------------------------------

call output_physiqex(debut,zjulian,pdtphys,presnivs,paprs,u,v,t,qx,ZCLDFR,ZQR,ZQS,ZQG,PTKEM)

! if lastcall, then it is time to write "restartphy.nc" file
if (lafin) then
  call phyredem("restartphy.nc")
endif


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
