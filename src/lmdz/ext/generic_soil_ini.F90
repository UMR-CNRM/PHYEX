module generic_soil_ini


! -------------------------------------------------------------------
!
! Object : initializing the "generic" soil model
!    Inheritating from a model developed in F Hourdin PhD (1992)
!    Tu be used in the PHYEX version of LMDZ
!
! Autor : Frédéric Hourdin
!
! First version of the generic soil model : Lege, Mai 2024
! -------------------------------------------------------------------

real, save   :: soil_min_period=2000., soil_dalpha=2.
!$OMP THREADPRIVATE(soil_min_period,soil_dalpha,nsoilmx)

real, allocatable, save :: zb(:),dz2(:),zc(:,:),zd(:,:)
real, save :: lambda
!$OMP THREADPRIVATE(zb,dz2,zc,zd, lambda)

real, save, allocatable, dimension(:) :: thermal_inertia,z0m,z0h,z0q,beta
!$OMP THREADPRIVATE(thermal_inertia,z0m,z0h,z0q,beta)

real, save :: cpd,lvtt,karman
!$OMP THREADPRIVATE(cpd,lvtt)

contains
  
    subroutine soil_ini(ngrid,nsoil,cpd_,lvtt_,thermal_inertia0,z0m0,z0h0,z0q0,beta0)

    USE ioipsl_getin_p_mod, ONLY : getin_p
    implicit none

    integer, intent(in) :: ngrid,nsoil
    real, intent(in) :: cpd_,lvtt_
    real, intent(in) :: thermal_inertia0,z0m0,z0h0,z0q0,beta0

    real :: fz,rk,fz1,rk1,rk2
    fz(rk)=fz1*(soil_dalpha**rk-1.)/(soil_dalpha-1.)

    integer :: jk

    cpd=cpd_
    lvtt=lvtt_
    karman=0.4

    !CALL getin_p('nsoilmx',nsoilmx)
    !CALL getin_p('soil_min_period',soil_min_period)
    !CALL getin_p('soil_dalpha',soil_dalpha)

    allocate(thermal_inertia(ngrid))
    allocate(z0m(ngrid),z0h(ngrid),z0q(ngrid),beta(ngrid))

    thermal_inertia=thermal_inertia0
    z0m(:)=z0m0
    z0h(:)=z0h0
    z0q(:)=z0q0
    beta(:)=beta0

! ---------------------------------------------------
! Specification of the soil vertical discretization
! ---------------------------------------------------

    allocate(zb(nsoil+1),dz2(nsoil))
    allocate(zc(ngrid,nsoil),zd(ngrid,nsoil))
!   la premiere couche represente un dixieme de cycle diurne
    fz1=sqrt(soil_min_period/3.14)

    ! dz2 delta_z des couches
    do jk=1,nsoil
       rk1=jk
       rk2=jk-1
       dz2(jk)=fz(rk1)-fz(rk2)
    enddo

    ! zb 1/delta_z entre les demi couches
    do jk=1,nsoil-1
       rk1=jk+.5
       rk2=jk-.5
       zb(jk+1)=1./(fz(rk1)-fz(rk2))
    enddo

    ! dz de la première demi couche
    lambda=fz(.5)*zb(2)
    zb(1)=1./fz(.5) ! pour le noveau schéma implicit
    PRINT*,'full layers, intermediate layers (secoonds)'

    ! Diagnostics
    do jk=1,nsoil
       rk=jk
       rk1=jk+.5
       rk2=jk-.5
       PRINT*,fz(rk1)*fz(rk2)*3.14,fz(rk)*fz(rk)*3.14
    enddo

    !do jk=1,nsoil
    !   print*,' dans soil ',jk
    !   rk=jk
    !   zout(jk)=fz(rk)
    !enddo

    end subroutine soil_ini

end module generic_soil_ini
