MODULE physiqex_mod

IMPLICIT NONE

CONTAINS

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
      USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
      USE MODD_PHYEX,      ONLY: PHYEX_t
      USE MODI_INI_PHYEX,  ONLY: INI_PHYEX
      USE MODI_ICE_ADJUST, ONLY: ICE_ADJUST
      USE MODD_BUDGET, ONLY: TBUCONF_ASSOCIATE, NBUDGET_RH, TBUCONF, LBU_ENABLE, &
                           & LBUDGET_U, LBUDGET_V, LBUDGET_W, LBUDGET_TH, &
                           & LBUDGET_TKE, LBUDGET_RV, LBUDGET_RC, LBUDGET_RR, &
                           & LBUDGET_RI, LBUDGET_RS, LBUDGET_RG, LBUDGET_RH, LBUDGET_SV, &
                           & TBUDGETDATA
      USE MODD_IO,         ONLY: TFILEDATA
      USE MODD_LES,        ONLY: TLES_t
      USE MODD_PARAMETERS, ONLY: JPVEXT_TURB
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

!    include "clesphys.h"
     include "flux_arp.h"



    INTEGER        length
    PARAMETER    ( length = 100 )
    REAL tabcntr0( length       )
    INTEGER, PARAMETER :: longcles=20
    REAL, SAVE :: clesphy0(longcles)
    TYPE(TLES_t)             :: TLES
    TYPE(DIMPHYEX_t),SAVE    :: D
    TYPE(PHYEX_t), SAVE      :: PHYEX
    real, dimension(klon, klev + 2) :: zqdm !1-qt=1/(1+rt)
    real, dimension(klon, klev + 2) :: zqt !mixing ratio and total specifiq content
    real, dimension(klon, klev + 2) :: ztheta, zexn !theta and exner function
    real, dimension(klon, klev + 2) :: zthetas !tendency
    real, dimension(klon, klev + 2) :: zthetas0
    real, dimension(klon, klev + 2) :: zz_mass !altitude above ground of mass points
    real, dimension(klon, klev + 2) :: zz_flux !altitude above ground of flux points
    real, dimension(klon, klev + 2) :: zdzm !distance between two consecutive mass points (expressed on flux points)
    real, dimension(klon, klev + 2) :: zdzf !distance between two consecutive flux points (expressed on mass points)
    real, dimension(klon) :: zs !surface orography
    real, dimension(klon) :: zsigqsat !coeff for the extra term for s variance
    real, dimension(0) :: zmfconv
    real, dimension(klon, klev+2) :: ZICLDFR, ZWCLDFR, ZSSIO, ZSSIU, ZIFR, ZICE_CLD_WGT !used only in HIRLAM config
    real, dimension(klon, klev+2) :: ZSRC ! bar(s'rc')
    real, dimension(klon, klev+2) :: ZCLDFR !cloud fraction
    real, dimension(klon, klev+2) :: ZHLC_HRC, ZHLC_HCF, ZHLI_HRI, ZHLI_HCF !subgrid autoconversion

    TYPE(TBUDGETDATA), DIMENSION(NBUDGET_RH), SAVE :: YLBUDGET

    !$OMP THREADPRIVATE(clesphy0)
! PHYEX variables
  REAL :: PTSTEP
  INTEGER, PARAMETER       :: KRR=6, KRRL=2, KRRI=3, KSV=0
  CHARACTER(LEN=4)         :: HBUNAME
  LOGICAL                  :: LMFCONV
  LOGICAL                  :: OCOMPUTE_SRC
  LOGICAL                  :: ONOMIXLG
  INTEGER                  :: KSV_LGBEG, KSV_LGEND
  REAL                     :: PDX, PDY
  INTEGER                  :: KMI, KSPLIT, KGRADIENTS, KHALO
  CHARACTER(LEN=4),DIMENSION(2)  :: HLBCX, HLBCY
  CHARACTER(LEN=6)         :: CPROGRAM
  INTEGER                  :: KSV_LIMA_NR, KSV_LIMA_NS, KSV_LIMA_NG, KSV_LIMA_NH
  LOGICAL                  :: O2D, OFLAT, OCOUPLES, OBLOWSNOW, OOCEAN, ODEEPOC
  LOGICAL                  :: OIBM, OFLYER
  TYPE(TFILEDATA)          :: ZTFILE
  REAL                     :: ZCEI_MAX, ZCEI_MIN, ZCOEF_AMPL_SAT
  REAL                     :: PRSNOW
  LOGICAL                  :: ODIAG_IN_RUN
  CHARACTER(LEN=4)         :: HTURBLEN_CL
  CHARACTER(LEN=4)         :: CMICRO
  INTEGER                  :: JRR
!
REAL, DIMENSION(klon,klev+2)   :: PDXX,PDYY,PDZX,PDZY
                                        ! metric coefficients
REAL, DIMENSION(klon,klev+2)   :: PZZ       !  physical distance
! between 2 succesive grid points along the K direction
REAL, DIMENSION(klon)      ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(klon)   ::  PCOSSLOPE       ! cosinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(klon)   ::  PSINSLOPE       ! sinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(klon,klev+2)      ::  PRHODJ    ! dry density * Grid size
REAL, DIMENSION(klon,klev+2)      ::  MFMOIST ! moist mass flux dual scheme
REAL, DIMENSION(klon,klev+2)      ::  PTHVREF   ! Virtual Potential
                                        ! Temperature of the reference state
REAL, DIMENSION(klon,klev+2,0) ::  PHGRAD      ! horizontal gradients
!
REAL, DIMENSION(klon)      ::  PSFTH,PSFRV,   &
! normal surface fluxes of theta and Rv
                                            PSFU,PSFV
! normal surface fluxes of (u,v) parallel to the orography
REAL, DIMENSION(klon,0)      ::  PSFSV
! normal surface fluxes of Scalar var.
!
!    prognostic variables at t- deltat
REAL, DIMENSION(:,:), ALLOCATABLE, SAVE ::  PTKEM       ! TKE
REAL, DIMENSION(klon,klev+2,0) ::  ZSVT        ! passive scal. var.
REAL, DIMENSION(klon) :: ZBL_DEPTH  ! BL height for TOMS
REAL, DIMENSION(klon) :: ZSBL_DEPTH ! SBL depth for RMC01
!
!    variables for cloud mixing length
REAL, DIMENSION(klon,klev+2)      ::  ZCEI
                                                 ! Cloud Entrainment instability
                                                 ! index to emphasize localy
                                                 ! turbulent fluxes
!
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(klon,klev+2) ::  ZUT, ZVT ! wind component on klev+2
REAL, DIMENSION(klon,klev+2) ::  ZPABST       ! absolute pressure
REAL, DIMENSION(klon,klev+2,KRR) ::  ZRX,ZRXS,ZRXS0 ! qx and source of q from LMDZ to rt
REAL, DIMENSION(klon,klev+2) :: PWT ! vertical wind velocity (used only for diagnostic)
!
! sources of momentum, conservative potential temperature, Turb. Kin. Energy,
! TKE dissipation
REAL, DIMENSION(klon,klev+2) ::  ZRUS,ZRVS,ZRWS,ZRTHS,ZRTKES
! Source terms for all water kinds, PRRS(:,:,:,1) is used for the conservative
! mixing ratio
REAL, DIMENSION(klon,klev+2,KRR) ::  ZRRS
! Source terms for all passive scalar variables
REAL, DIMENSION(klon,klev+2,0) ::  ZRSVS
! Sigma_s at time t+1 : square root of the variance of the deviation to the
! saturation
REAL, DIMENSION(klon,klev+2)     ::  ZFLXZTHVMF
!                                           MF contribution for vert. turb. transport
!                                           used in the buoy. prod. of TKE
REAL, DIMENSION(klon,klev+2)  :: ZWTH       ! heat flux
REAL, DIMENSION(klon,klev+2)  :: ZWRC       ! cloud water flux
REAL, DIMENSION(klon,klev+2,0):: ZWSV       ! scalar flux
REAL, DIMENSION(klon,klev+2)  :: ZTP        ! Thermal TKE production
                                                   ! MassFlux + turb
REAL, DIMENSION(klon,klev+2)  :: ZDP        ! Dynamic TKE production
REAL, DIMENSION(klon,klev+2)  :: ZTDIFF     ! Diffusion TKE term
REAL, DIMENSION(klon,klev+2)  :: ZTDISS     ! Dissipation TKE term
!
! length scale from vdfexcu
REAL, DIMENSION(klon,klev+2)    :: PLENGTHM, PLENGTHH
!
!
REAL, DIMENSION(klon,klev+2) :: ZDXX,ZDYY,ZDZX,ZDZY,ZZZ
REAL, DIMENSION(klon)      :: ZDIRCOSXW,ZDIRCOSYW,ZDIRCOSZW,ZCOSSLOPE,ZSINSLOPE
!
REAL, DIMENSION(klon,klev+2) ::  ZRHOD !rho_dry
!
REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: PSIGS !variance of s
REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: PCF_MF, PRC_MF, PRI_MF !shallow convection cloud
!
! Shallow specific variables 
REAL, DIMENSION(klon,klev+2)      ::  PRHODREF    ! dry density reference profile
REAL, DIMENSION(klon,klev+2)::  PDUDT_MF     ! tendency of U   by massflux scheme
REAL, DIMENSION(klon,klev+2)::  PDVDT_MF     ! tendency of V   by massflux scheme
REAL, DIMENSION(klon,klev+2)::  PDTHLDT_MF   ! tendency of thl by massflux scheme
REAL, DIMENSION(klon,klev+2)::  PDRTDT_MF    ! tendency of rt  by massflux scheme
REAL, DIMENSION(klon,klev+2,KSV)::  PDSVDT_MF    ! tendency of Sv  by massflux scheme

REAL, DIMENSION(klon,klev+2)     ::  PSIGMF
REAL, DIMENSION(klon,klev+2)     ::  ZFLXZTHMF
REAL, DIMENSION(klon,klev+2)     ::  ZFLXZRMF
REAL, DIMENSION(klon,klev+2)     ::  ZFLXZUMF
REAL, DIMENSION(klon,klev+2)     ::  ZFLXZVMF
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
  
  ! Initialize PHYEX
  CALL INI_PHYEX(HPROGRAM='AROME ', KUNITNML=0, LDNEEDNAM=.TRUE., &
                &KLUOUT=20, KFROM=0, KTO=1, &
                &PTSTEP=pdtphys, PDZMIN=999., &
                &CMICRO='ICE3', CSCONV='EDKF', CTURB='TKEL', &
                &LDDEFAULTVAL=.TRUE., LDREADNAM=.FALSE., LDCHECK=.FALSE., &
                &KPRINT=0, LDINIT=.TRUE., &
                &PHYEX_OUT=PHYEX)

  PHYEX%NEBN%LSUBG_COND = .TRUE. 

  D%NIT  = klon
  D%NIB  = 1
  D%NIE  = klon
  D%NJT  = 1
  D%NJB  = 1
  D%NJE  = 1
  D%NIJT = D%NIT * D%NJT
  D%NIJB = 1
  D%NIJE = 1
  D%NKL  = 1
  D%NKT  = klev+2
  D%NKA  = 1
  D%NKU  = klev+2
  D%NKB  = 1+JPVEXT_TURB
  D%NKE  = D%NKT - JPVEXT_TURB
  D%NKTB = 1+JPVEXT_TURB
  D%NKTE = D%NKT - JPVEXT_TURB
  D%NIBC = 1
  D%NJBC = 1
  D%NIEC = D%NIE
  D%NJEC = D%NJT
  
  !Budgets
  CALL TBUCONF_ASSOCIATE
  DO K=1, NBUDGET_RH
    YLBUDGET(K)%NBUDGET=K
  ENDDO
  LBU_ENABLE=.FALSE.
  LBUDGET_U=.FALSE.
  LBUDGET_V=.FALSE.
  LBUDGET_W=.FALSE.
  LBUDGET_TH=.FALSE.
  LBUDGET_TKE=.FALSE.
  LBUDGET_RV=.FALSE.
  LBUDGET_RC=.FALSE.
  LBUDGET_RR=.FALSE.
  LBUDGET_RI=.FALSE.
  LBUDGET_RS=.FALSE.
  LBUDGET_RG=.FALSE.
  LBUDGET_RH=.FALSE.
  LBUDGET_SV=.FALSE.

! TKE is not advected yet (TODO), for now, init at TKEMIN
ALLOCATE(PTKEM(klon,klev+2))
PTKEM(:,:) = PHYEX%TURBN%XTKEMIN

#ifndef CPP_IOIPSL_NO_OUTPUT
  ! Initialize IOIPSL output file
#endif

  ALLOCATE(PSIGS(klon,klev+2))
  ALLOCATE(PCF_MF(klon,klev+2))
  ALLOCATE(PRC_MF(klon,klev+2))
  ALLOCATE(PRI_MF(klon,klev+2))
  PSIGS=0.
  PCF_MF=0.
  PRC_MF=0.
  PRI_MF=0.
endif ! of if (debut)

!------------------------------------------------------------
! Initialisations a chaque pas de temps
!------------------------------------------------------------

! set all tendencies to zero
d_u(1:klon,1:klev)=0.
d_v(1:klon,1:klev)=0.
d_t(1:klon,1:klev)=0.
d_qx(1:klon,1:klev,1:nqtot)=0.
d_ps(1:klon)=0.

!------------------------------------------------------------
! Conversions and extra levels
!------------------------------------------------------------
!TODO check in Meso-NH how values are extrapolated outside of the physical domain
zqt(:,2:klev+1) = qx(:,:,1) + qx(:,:,2) + qx(:,:,3)
zqt(:,1)=0.
zqt(:,klev+2)=0.
zqdm(:,:)=1.-zqt(:,:)
ZRX(:,2:klev+1,1) = qx(:,:,1) / zqdm(:,2:klev+1)
ZRX(:,2:klev+1,2) = qx(:,:,2) / zqdm(:,2:klev+1)
ZRX(:,2:klev+1,4) = qx(:,:,3) / zqdm(:,2:klev+1)

ZRX(:,:,5:KRR) = 0. ! TODO init of snow, graupel, hail
ZRX(:,:,3) = 0. ! TODO init of rain
!
!TODO add hydrometeors
ZRX(:,1,1)=0.
ZRX(:,klev+2,1)=0.
ZRX(:,1,2)=0.
ZRX(:,klev+2,2)=0.
ZRX(:,1,4)=0.
ZRX(:,klev+2,4)=0.

ZRXS(:,:,:) = ZRX(:,:,:)/pdtphys

zexn(:,2:klev+1) = (pplay(:,:) / PHYEX%CST%XP00) ** (PHYEX%CST%XRD/PHYEX%CST%XCPD)
ztheta(:,2:klev+1) = t(:,:) / zexn(:,2:klev+1)
CALL VERTICAL_EXTEND(zexn,klev)
CALL VERTICAL_EXTEND(ztheta,klev)
zthetas(:,:)=ztheta(:,:)/pdtphys

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
DO k=2,klev+1
  PRHODJ(:,k) = pplay(:,k-1) / (t(:,k-1) * 287.0) * (zdzf(:,k)*cell_area(:))
END DO
PRHODREF(:,2:klev+1) = pplay(:,:) / (t(:,:) * 287.0)
CALL VERTICAL_EXTEND(ZPABST,klev)
CALL VERTICAL_EXTEND(PRHODJ,klev)
CALL VERTICAL_EXTEND(PRHODREF,klev)

ZRHOD(:,2:klev+1)=ZPABST(:,2:klev+1)/(t*(PHYEX%CST%XRD+ZRX(:,2:klev+1,1)*PHYEX%CST%XRV))
CALL VERTICAL_EXTEND(ZRHOD,klev)

!
!------------------------------------------------------------
! Adjustment
!------------------------------------------------------------
!
ZRXS0(:,:,:) = ZRXS(:,:,:)
ZSRC(:,:) = 0.
ZTHETAS0=ZTHETAS
ZSIGQSAT=PHYEX%NEBN%VSIGQSAT
CALL ICE_ADJUST (D, PHYEX%CST, PHYEX%RAIN_ICE_PARAMN, PHYEX%NEBN, PHYEX%TURBN, PHYEX%PARAM_ICEN, TBUCONF, KRR,   &
                &'ADJU',                                          &
                &pdtphys, ZSIGQSAT,                                 &
                &PRHODJ, zexn, ZRHOD, PSIGS, .FALSE., zmfconv,&
                &ZPABST, ZZ_MASS,                                      &
                &zexn, PCF_MF, PRC_MF, PRI_MF,                     &
                &ZICLDFR, ZWCLDFR, ZSSIO, ZSSIU, ZIFR,             &
                &ZRX(:,:,1), ZRX(:,:,2), ZRXS(:,:,1), ZRXS(:,:,2), ztheta, ZTHETAS,                  &
                &.TRUE., ZSRC, ZCLDFR,                      &
                &ZRX(:,:,3), ZRX(:,:,4), ZRXS(:,:,4), ZRX(:,:,5), ZRX(:,:,6), YLBUDGET, NBUDGET_RH,     &
                &ZICE_CLD_WGT,                                     &
                !&POUT_RV, POUT_RC, POUT_RI, POUT_TH,               &
                &ZHLC_HRC, ZHLC_HCF, ZHLI_HRI, ZHLI_HCF)
! Tendencies, mixing ratio -> specific
d_qx(:,1:klev,1)=d_qx(:,1:klev,1) + (ZRXS(:,2:klev+1,1)-ZRXS0(:,2:klev+1,1))*ZQDM(:,2:klev+1)
d_qx(:,1:klev,2)=d_qx(:,1:klev,2) + (ZRXS(:,2:klev+1,2)-ZRXS0(:,2:klev+1,2))*ZQDM(:,2:klev+1)
d_qx(:,1:klev,3)=d_qx(:,1:klev,3) + (ZRXS(:,2:klev+1,4)-ZRXS0(:,2:klev+1,4))*ZQDM(:,2:klev+1)
d_t(:,1:klev)=d_t(:,1:klev) + (zthetas(:,2:klev+1)-zthetas0(:,2:klev+1))*zexn(:,2:klev+1)
!
! compute tendencies to return to the dynamics:
! "friction" on the first layer
d_u(1:klon,1)=-u(1:klon,1)/86400.
d_v(1:klon,1)=-v(1:klon,1)/86400.
! newtonian relaxation towards temp_newton()
!do k=1,klev
!  temp_newton(1:klon,k)=280.+cos(latitude(1:klon))*40.-pphi(1:klon,k)/rg*6.e-3
!  d_t(1:klon,k)=d_t(1:klon,k) + (temp_newton(1:klon,k)-t(1:klon,k))/1.e5
!enddo


KSV_LGBEG = 0
KSV_LGEND = 0
ONOMIXLG=.FALSE.
KMI = 1
KGRADIENTS =0
HLBCX(:)=(/'CYCL','CYCL'/)
HLBCY(:)=(/'CYCL','CYCL'/)
KSPLIT = 1
KHALO=1
CPROGRAM='AROME '
O2D=.FALSE.
OFLAT=.FALSE.
OCOUPLES=.FALSE.
OBLOWSNOW=.FALSE.
OCOMPUTE_SRC=.TRUE.
OOCEAN=.FALSE.
ODEEPOC=.FALSE.
ZTFILE%LOPENED=.FALSE.
ZCEI_MAX=1.0
ZCEI_MIN=0.0
ZCOEF_AMPL_SAT=0.0
KSV_LIMA_NR=0
KSV_LIMA_NS=0
KSV_LIMA_NG=0
KSV_LIMA_NH=0
OIBM=.FALSE.
OFLYER=.FALSE.
PRSNOW=1.0
ODIAG_IN_RUN=.FALSE.
HTURBLEN_CL='NONE'
CMICRO='ICE3'
TLES%LLES=.FALSE.
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
!
! Flux surface RICO
PSFTH(:) = 5E-3 ! RICO
PSFRV(:) = 6E-5 ! RICO

! ARMCU
!PSFTH(:) = -fsens/1000.
!PSFRV(:) = -flat/(2.5e6)
!
PSFSV(:,:) = 0.
PSFU(:) = 0.
PSFV(:) = 0.
!
ZSVT(:,:,:) = 0.
PWT(:,:) = 0.
! needed only with LRMC01 key (correction in the surface boundary layer)
ZBL_DEPTH(:) = 0.
ZSBL_DEPTH(:) = 0.
! needed only if HTURBLEN_CL /= 'NONE' modification of mixing lengh inside clouds
ZCEI(:,:) = 0.
!
ZUT(:,2:klev+1) = u(:,:)
ZVT(:,2:klev+1) = v(:,:)
PTHVREF(:,:) = ztheta(:,:) ! profil de theta_v a calculer TODO 
!
CALL VERTICAL_EXTEND(ZUT,klev)
CALL VERTICAL_EXTEND(ZVT,klev)
CALL VERTICAL_EXTEND(PTHVREF,klev)
DO JRR=1,KRR
  CALL VERTICAL_EXTEND(ZRX(:,:,JRR),klev)
END DO 
!
!TODO PSIGMF option STAT
!------------------------------------------------------------
! Shallow convection
!------------------------------------------------------------
!
  CALL SHALLOW_MF(D, PHYEX%CST, PHYEX%NEBN, PHYEX%PARAM_MFSHALLN, PHYEX%TURBN, PHYEX%CSTURB,                    &
     &KRR=KRR, KRRL=KRRL, KRRI=KRRI, KSV=KSV,              &
     &ONOMIXLG=ONOMIXLG,KSV_LGBEG=KSV_LGBEG,KSV_LGEND=KSV_LGEND,      &
     &PTSTEP=pdtphys, &
     &PDZZ=zdzm(:,:),PZZ=zz_mass(:,:),                                                                 &
     &PRHODJ=PRHODJ(:,:),PRHODREF=PRHODREF(:,:),                                                    &
     &PPABSM=ZPABST(:,:),PEXNM=zexn(:,:),                                                          &
     &PSFTH=PSFTH(:),PSFRV=PSFRV(:),                                                            &
     &PTHM=ztheta(:,:),PRM=ZRX(:,:,:),PUM=ZUT(:,:),PVM=ZVT(:,:),&
     &PTKEM=PTKEM(:,:),PSVM=ZSVT(:,:,:),                            &
     &PDUDT_MF=PDUDT_MF(:,:),PDVDT_MF=PDVDT_MF(:,:),                                                &
     &PDTHLDT_MF=PDTHLDT_MF(:,:),PDRTDT_MF=PDRTDT_MF(:,:),PDSVDT_MF=PDSVDT_MF(:,:,:),                      &
     &PSIGMF=PSIGMF(:,:),PRC_MF=PRC_MF(:,:),PRI_MF=PRI_MF(:,:),PCF_MF=PCF_MF(:,:),&
     &PFLXZTHVMF=ZFLXZTHVMF(:,:),      &
     &PFLXZTHMF=ZFLXZTHMF(:,:),PFLXZRMF=ZFLXZRMF(:,:),PFLXZUMF=ZFLXZUMF(:,:),PFLXZVMF=ZFLXZVMF(:,:),     &
     &PTHL_UP=PTHL_UP(:,:),PRT_UP=PRT_UP(:,:),PRV_UP=PRV_UP(:,:),&
     &PRC_UP=PRC_UP(:,:),PRI_UP=PRI_UP(:,:),            &
     &PU_UP=PU_UP(:,:), PV_UP=PV_UP(:,:), PTHV_UP=PTHV_UP(:,:), PW_UP=PW_UP(:,:),                        &
     &PFRAC_UP=PFRAC_UP(:,:),PEMF=PEMF(:,:),PDETR=PDETR(:,:),PENTR=PENTR(:,:),                           &
     &KKLCL=IKLCL(:),KKETL=IKETL(:),KKCTL=IKCTL(:),PDX=1.0,PDY=1.0,KBUDGETS=NBUDGET_RH )

! Add tendencies of shallow to total physics tendency
d_u(:,1:klev) = d_u(:,1:klev) + PDUDT_MF(:,2:klev+1)
d_v(:,1:klev) = d_v(:,1:klev) + PDVDT_MF(:,2:klev+1) 
d_t(:,1:klev) = d_t(:,1:klev) + PDTHLDT_MF(:,2:klev+1)*zexn(:,2:klev+1) !TODO Theta_l en theta ?
d_qx(:,1:klev,1)=d_qx(:,1:klev,1) + PDRTDT_MF(:,2:klev+1)*zqdm(:,2:klev+1)
! TODO add SV tendencies
!
!------------------------------------------------------------
! Turbulence
!------------------------------------------------------------
! out tendencies
ZRUS(:,:) = 0.
ZRVS(:,:) = 0.
ZRWS(:,:) = 0.
ZRTHS(:,:) = 0.
ZRRS(:,:,:) = 0.
ZRSVS(:,:,:) = 0.
ZRTKES(:,:) = 0.
CALL TURB(PHYEX%CST, PHYEX%CSTURB, TBUCONF, PHYEX%TURBN, PHYEX%NEBN, D, TLES,&
   & KMI, KRR, KRRL, KRRI, HLBCX, HLBCY, KGRADIENTS, KHALO,&
   & KSPLIT,KMI, KSV, KSV_LGBEG, KSV_LGEND, &
   & CPROGRAM, &
   & KSV_LIMA_NR, KSV_LIMA_NS, KSV_LIMA_NG, KSV_LIMA_NH,&
   & O2D, ONOMIXLG, OFLAT, OCOUPLES, OBLOWSNOW,OIBM,&
   & OFLYER, OCOMPUTE_SRC, PRSNOW, &
   & OOCEAN, ODEEPOC, ODIAG_IN_RUN,   &
   & HTURBLEN_CL,CMICRO,           &
   & pdtphys,ZTFILE,                                       &
   & ZDXX(:,:),ZDYY(:,:),zdzm(:,:),      &
   & ZDZX(:,:),ZDZY(:,:),zz_flux(:,:),       &
   & ZDIRCOSXW(:),ZDIRCOSYW(:),ZDIRCOSZW(:),ZCOSSLOPE(:),ZSINSLOPE(:),    &
   & PRHODJ(:,:),PTHVREF(:,:), PHGRAD(:,:,:), zs(:),    &
   & PSFTH(:),PSFRV(:),PSFSV(:,:),PSFU(:),PSFV(:), &
   & ZPABST(:,:),ZUT(:,:),ZVT(:,:),PWT(:,:),PTKEM(:,:),ZSVT(:,:,:),ZSRC(:,:),&
   & PLENGTHM(:,:),PLENGTHH(:,:),MFMOIST(:,:),                            &
   & ZBL_DEPTH(:),ZSBL_DEPTH(:),                                 &
   & ZCEI(:,:),ZCEI_MIN,ZCEI_MAX,ZCOEF_AMPL_SAT,    &
   & ztheta(:,:),ZRX(:,:,:), &
   & ZRUS(:,:),ZRVS(:,:),ZRWS(:,:),ZRTHS(:,:),ZRRS(:,:,:),ZRSVS(:,:,:),ZRTKES(:,:),         &
   & PSIGS(:,:),                                         &
   & ZFLXZTHVMF(:,:),ZWTH(:,:),ZWRC(:,:),ZWSV(:,:,:),ZDP(:,:),ZTP(:,:),ZTDIFF(:,:),ZTDISS(:,:),&
   & YLBUDGET, NBUDGET_RH)

! Add tendencies of turb to total physics tendency
d_u(:,1:klev) = d_u(:,1:klev) + ZRUS(:,2:klev+1)/PRHODJ(:,2:klev+1)
d_v(:,1:klev) = d_v(:,1:klev) + ZRVS(:,2:klev+1)/PRHODJ(:,2:klev+1)
d_t(:,1:klev) = d_t(:,1:klev) + ZRTHS(:,2:klev+1)*zexn(:,2:klev+1)/PRHODJ(:,2:klev+1)
!
d_qx(:,1:klev,1)=d_qx(:,1:klev,1) + ZRRS(:,2:klev+1,1)/PRHODJ(:,2:klev+1)*zqdm(:,2:klev+1) !(ZRXS(:,2:klev+1,1)-ZRXS0(:,2:klev+1,1))*ZQDM(:,2:klev+1)
d_qx(:,1:klev,2)=d_qx(:,1:klev,2) + ZRRS(:,2:klev+1,2)/PRHODJ(:,2:klev+1)*zqdm(:,2:klev+1) !(ZRXS(:,2:klev+1,2)-ZRXS0(:,2:klev+1,2))*ZQDM(:,2:klev+1)
d_qx(:,1:klev,3)=d_qx(:,1:klev,3) + ZRRS(:,2:klev+1,4)/PRHODJ(:,2:klev+1)*zqdm(:,2:klev+1) !(ZRXS(:,2:klev+1,4)-ZRXS0(:,2:klev+1,4))*ZQDM(:,2:klev+1)
!
PTKEM(:,:) = PTKEM(:,:) + ZRTKES(:,:)/PRHODJ(:,:)*pdtphys
!
!------------------------------------------------------------
! Entrees sorties
!------------------------------------------------------------

call output_physiqex(debut,zjulian,pdtphys,presnivs,paprs,u,v,t,qx,ZCLDFR)

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
