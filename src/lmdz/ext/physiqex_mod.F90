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
      USE MODD_BUDGET, ONLY: NBUDGET_RH, TBUDGETDATA, TBUDGETCONF_t
      USE MODD_IO,         ONLY: TFILEDATA
      USE MODD_LES,        ONLY: TLES_t
      USE MODD_PARAMETERS, ONLY: JPVEXT_TURB
      USE MODI_TURB
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
      real,intent(in) :: flxmass_w(klon,klev+2) ! vertical mass flux
      real,intent(out) :: d_u(klon,klev+2) ! physics tendency on u (m/s/s)
      real,intent(out) :: d_v(klon,klev+2) ! physics tendency on v (m/s/s)
      real,intent(out) :: d_t(klon,klev+2) ! physics tendency on t (K/s)
      real,intent(out) :: d_qx(klon,klev+2,nqtot) ! physics tendency on tracers
      real,intent(out) :: d_ps(klon) ! physics tendency on surface pressure

!    include "clesphys.h"
    INTEGER        length
    PARAMETER    ( length = 100 )
    REAL tabcntr0( length       )
    INTEGER, PARAMETER :: longcles=20
    REAL, SAVE :: clesphy0(longcles)
    TYPE(TLES_t)             :: TLES
    TYPE(DIMPHYEX_t),SAVE    :: D
    TYPE(PHYEX_t), SAVE      :: PHYEX
    real, dimension(klon, klev + 2) :: zrv, zrc, zri, zrr, zrs, zrg, zqt, zqdm !mixing ratio and total specifiq content
    real, dimension(klon, klev + 2) :: ztheta, zexn !theta and exner function
    real, dimension(klon, klev + 2) :: zz_mass !altitude above ground of mass points
    real, dimension(klon, klev + 2) :: zz_flux !altitude above ground of flux points
    real, dimension(klon, klev + 2) :: zdzm !distance between two consecutive mass points (expressed on flux points)
    real, dimension(klon, klev + 2) :: zdzf !distance between two consecutive flux points (expressed on mass points)
    real, dimension(klon) :: zs !surface orography

    !$OMP THREADPRIVATE(clesphy0)
! PHYEX variables
  REAL :: PTSTEP
  INTEGER, PARAMETER       :: KRR= 6
  CHARACTER(LEN=4)         :: HBUNAME
  LOGICAL                  :: LMFCONV
  INTEGER                  :: KRRL, KRRI, KSV
  LOGICAL                  :: OCOMPUTE_SRC
  TYPE(TBUDGETDATA), DIMENSION(NBUDGET_RH) :: YLBUDGET
  TYPE(TBUDGETCONF_t)      :: TBUCONF
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
REAL, DIMENSION(klon,klev+2) ::  ZSRCT       ! Second-order flux
                      ! s'rc'/2Sigma_s2 at time t-1 multiplied by Lambda_3
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
REAL, DIMENSION(klon,klev+2,KRR) ::  ZRX ! qx from LMDZ
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
REAL, DIMENSION(klon,klev+2)     ::  ZSIGS
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

real :: temp_newton(klon,klev+2)
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

! TKE is not advected yet (TODO), for now, init at TKEMIN
ALLOCATE(PTKEM(klon,klev+2))
PTKEM(:,:) = PHYEX%TURBN%XTKEMIN

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
d_ps(1:klon)=0.

!------------------------------------------------------------
! Conversions
!------------------------------------------------------------
!TODO check in Meso-NH how values are extrapolated outside of the physical domain
zqt(:,2:klev+1) = qx(:,:,1) + qx(:,:,2) + qx(:,:,3)
zrv(:,2:klev+1) = qx(:,:,1) / (1 - zqt(:,2:klev+1))
zrc(:,2:klev+1) = qx(:,:,2) / (1 - zqt(:,2:klev+1))
zri(:,2:klev+1) = qx(:,:,3) / (1 - zqt(:,2:klev+1))

ZRX(:,:,1) = zrv(:,:)
ZRX(:,:,2) = zrc(:,:)
ZRX(:,:,3) = zri(:,:)
ZRX(:,:,4:KRR) = 0. ! TODO init of rain, graupel, hail

!zrr =
!zrs =
!zrg =
zqt(:,1)=0.
zqt(:,klev+2)=0.
zrv(:,1)=0.
zrv(:,klev+2)=0.
zrc(:,1)=0.
zrc(:,klev+2)=0.
zri(:,1)=0.
zri(:,klev+2)=0.
zqdm=1-zqt

zexn(:,2:klev+1) = (pplay / PHYEX%CST%XP00) ** (PHYEX%CST%XRD/PHYEX%CST%XCPD)
ztheta(:,2:klev+1) = t / zexn(:,2:klev+1)
CALL VERTICAL_EXTEND(zexn,klev)
CALL VERTICAL_EXTEND(ztheta,klev)

!TODO check in Meso-NH how zz_mass and zz_flux are initialized outside of the physical domain
zs = pphis/PHYEX%CST%XG
zz_mass(:,2:klev+1) = pphi / PHYEX%CST%XG
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

!------------------------------------------------------------
! Calculs
!------------------------------------------------------------

! compute tendencies to return to the dynamics:
! "friction" on the first layer
d_u(1:klon,1)=-u(1:klon,1)/86400.
d_v(1:klon,1)=-v(1:klon,1)/86400.
! newtonian relaxation towards temp_newton()
do k=1,klev
  temp_newton(1:klon,k)=280.+cos(latitude(1:klon))*40.-pphi(1:klon,k)/rg*6.e-3
  d_t(1:klon,k)=(temp_newton(1:klon,k)-t(1:klon,k))/1.e5
enddo


KSV_LGBEG = 0
KSV_LGEND = 0
ONOMIXLG=.FALSE.
KRRL = 2
KRRI = 3
KSV = 0
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
PSFSV(:,:) = 0.
PSFU(:) = 0.
PSFV(:) = 0.
!
ZSVT(:,:,:) = 0.
ZSRCT(:,:) =0.
PWT(:,:) = 0.
! needed only with LRMC01 key (correction in the surface boundary layer)
ZBL_DEPTH(:) = 0.
ZSBL_DEPTH(:) = 0.
! needed only if HTURBLEN_CL /= 'NONE' modification of mixing lengh inside clouds
ZCEI(:,:) = 0.
! out tendencies
ZRUS(:,:) = 0.
ZRVS(:,:) = 0.
ZRWS(:,:) = 0.
ZRTHS(:,:) = 0.
ZRRS(:,:,:) = 0.
ZRSVS(:,:,:) = 0.
ZRTKES(:,:) = 0.
!
ZUT(:,2:klev+1) = u(:,:)
ZVT(:,2:klev+1) = v(:,:)
ZPABST(:,2:klev+1) = pplay(:,:)
DO k=2,klev+1
  PRHODJ(:,k) = pplay(:,k-1) / (t(:,k-1) * 287.0) * (zdzf(:,k)*cell_area(:))
END DO
PTHVREF(:,:) = ztheta(:,:) ! profil de theta TODO ?
!ZFLXZTHVMF(:,2:klev+1) = flxmass_w(:,:) TODO : flxmass_w deja en klev+2 ???
!ZFLXZTHVMF(:,:) = flxmass_w(:,:) ! TODO to uncomment once a mass-flux scheme is plugged
ZFLXZTHVMF(:,:) = 0.
!
CALL VERTICAL_EXTEND(ZUT,klev)
CALL VERTICAL_EXTEND(ZVT,klev)
CALL VERTICAL_EXTEND(ZPABST,klev)
CALL VERTICAL_EXTEND(PTHVREF,klev)
CALL VERTICAL_EXTEND(PRHODJ,klev)
DO JRR=1,KRR
  CALL VERTICAL_EXTEND(ZRX(:,:,JRR),klev)
END DO 
!
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
   & ZPABST(:,:),ZUT(:,:),ZVT(:,:),PWT(:,:),PTKEM(:,:),ZSVT(:,:,:),ZSRCT(:,:),&
   & PLENGTHM(:,:),PLENGTHH(:,:),MFMOIST(:,:),                            &
   & ZBL_DEPTH(:),ZSBL_DEPTH(:),                                 &
   & ZCEI(:,:),ZCEI_MIN,ZCEI_MAX,ZCOEF_AMPL_SAT,    &
   & ztheta(:,:),ZRX(:,:,:), &
   & ZRUS(:,:),ZRVS(:,:),ZRWS(:,:),ZRTHS(:,:),ZRRS(:,:,:),ZRSVS(:,:,:),ZRTKES(:,:),         &
   & ZSIGS(:,:),                                         &
   & ZFLXZTHVMF(:,:),ZWTH(:,:),ZWRC(:,:),ZWSV(:,:,:),ZDP(:,:),ZTP(:,:),ZTDIFF(:,:),ZTDISS(:,:),&
   & YLBUDGET, NBUDGET_RH)

! Add tendencies of turb to total physics tendency
d_u(:,1:klev) = d_u(:,1:klev) + ZRUS(:,2:klev+1)/PRHODJ(:,2:klev+1)
d_v(:,1:klev) = d_v(:,1:klev) + ZRVS(:,2:klev+1)/PRHODJ(:,2:klev+1)
d_t(:,1:klev) = d_t(:,1:klev) + ZRTHS(:,2:klev+1)*zexn(:,2:klev+1)/PRHODJ(:,2:klev+1)
DO JRR=1,3 !3 for qv, ql, qi
  d_qx(:,1:klev,JRR)=d_qx(:,1:klev,JRR) + ZRRS(:,2:klev+1,JRR)/PRHODJ(:,2:klev+1)*zqdm(:,2:klev+1)
END DO
d_ps(1:klon)=0.

!------------------------------------------------------------
! Entrees sorties
!------------------------------------------------------------

call output_physiqex(debut,zjulian,pdtphys,presnivs,paprs,u,v,t)

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
