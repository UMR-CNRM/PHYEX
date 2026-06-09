module generic_soil_fluxes
contains

subroutine soil_fluxes(ngrid,tsrf,qsats,z1,rho1,u1,v1,t1,q1,fluxu,fluxv,fluxt,fluxq)
!=======================================================================
!
! Object : computing turbulent fluxes for the "generic" soil model
!    Inheritating from a model developed in F Hourdin PhD (1992)
!    Tu be used in the PHYEX version of LMDZ
!
! Autor :  Frederic Hourdin     04/06/2024
!
! First version of the generic soil model : Lege, Mai 2024
!
!=======================================================================
!   declarations:
!   -------------

use generic_soil_ini, only : karman,cpd,lvtt,beta,z0m,z0h,z0q

implicit none

integer, intent(in) :: ngrid
real, intent(in), dimension(ngrid) :: tsrf,qsats,z1,rho1,u1,v1,t1,q1
real, intent(out), dimension(ngrid) :: fluxu,fluxv,fluxt,fluxq

real, dimension(ngrid) :: cdm,cdh,cdq,rhov

!-----------------------------------------------------------------------
!  arguments
!  ---------


!-----------------------------------------------------------------------

cdm(:)=(karman/log(z1(:)/z0m(:)))**2
cdh(:)=(karman/log(z1(:)/z0h(:)))**2
cdq(:)=(karman/log(z1(:)/z0q(:)))**2
rhov(:)=rho1(:)*sqrt(u1(:)**2+v1(:)**2)

fluxt(:)=rhov(:)*cdh(:)*(tsrf(:)-t1(:))
   ! En l'absence de calcul de qsat, on suppose que LE=rho Cd V q1 * ( 1 -RH )
   !print*,'qsats',qsats,qx(:,1,1)
fluxq(:)=beta(:)*rhov(:)*cdq(:)*(qsats(:)-q1(:))

fluxu(:)=-rhov(:)*cdm(:)*u1(:)
fluxv(:)=-rhov(:)*cdm(:)*v1(:)


return
end subroutine soil_fluxes
end module generic_soil_fluxes
