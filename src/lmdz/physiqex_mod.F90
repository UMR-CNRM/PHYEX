! $Id: physiq.F 1565 2011-08-31 12:53:29Z jghattas $
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
      USE geometry_mod, only : latitude
!      USE comcstphy, only : rg
      USE ioipsl, only : ymds2ju
      USE phys_state_var_mod, only : phys_state_var_init
      USE phyetat0_mod, only: phyetat0
      USE output_physiqex_mod, ONLY: output_physiqex

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
    INTEGER        length
    PARAMETER    ( length = 100 )
    REAL tabcntr0( length       )
    INTEGER, PARAMETER :: longcles=20
    REAL, SAVE :: clesphy0(longcles)
    !$OMP THREADPRIVATE(clesphy0)


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


!------------------------------------------------------------
! Entrees sorties
!------------------------------------------------------------


call output_physiqex(debut,zjulian,pdtphys,presnivs,paprs,u,v,t)

! if lastcall, then it is time to write "restartphy.nc" file
if (lafin) then
  call phyredem("restartphy.nc")
endif


end subroutine physiqex

END MODULE physiqex_mod
