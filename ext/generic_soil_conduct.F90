      MODULE generic_soil_conduct
      CONTAINS

      SUBROUTINE soil_conduct(ngrid,nsoil,firstcall,ptimestep,ptsrf,ptsoil,pcapcal,pfluxgrd,zout)

!=======================================================================
!
! Object : computing thermal conduction for the "generic" soil model
!    Inheritating from a model developed in F Hourdin PhD (1992)
!    Tu be used in the PHYEX version of LMDZ
!
! Autor : Frédéric Hourdin
!
! First version of the generic soil model : Lege, Mai 2024
!
!   Original version : Frederic Hourdin     30/01/92
!   -------
!
!   object:  computation of : the soil temperature evolution
!   -------                   the surfacic heat capacity "Capcal"
!                             the surface conduction flux pcapcal
!
!
!   Method: implicit time integration
!   -------
!   Consecutive ground temperatures are related by:
!           T(k) = C(k) + D(k)*T(k-1)  (1)
!   the coefficients C and D are computed at the t-dt time-step.
!   Routine structure:
!   1)new temperatures are computed  using (1)
!   2)C and D coefficients are computed from the new temperature
!     profile for the t+dt time-step
!   3)the coefficients A and B are computed where the diffusive
!     fluxes at the t+dt time-step is given by
!            Fdiff = A + B Ts(t+dt)
!     or     Fdiff = F0 + Capcal (Ts(t+dt)-Ts(t))/dt
!            with F0 = A + B (Ts(t))
!                 Capcal = B*dt
!           
! Conductivité thermique λ (W/mK), Capacité thermique volumétrique ρC (MJ/m³K)
!   Argile/limon sec	0.4	1.5 – 1.6 
! λs  = 1 W m−1 K−1 et rho Cs = 1 × 10^6 J m−3 K−1
! I = sqrt ( lambda rho C )
! z = z' * sqrt( lambda / rho C )
!
!   Interface:
!   ----------
!
!   Arguments:
!   ----------
!   ngird               number of grid-points
!   ptimestep              physical timestep (s)
!   pto(ngrid,nsoil)     temperature at time-step t (K)
!   ptn(ngrid,nsoil)     temperature at time step t+dt (K)
!   pcapcal(ngrid)      specific heat (W*m-2*s*K-1)
!   pfluxgrd(ngrid)      surface diffusive flux from ground (Wm-2)
!   
!=======================================================================
!   declarations:
!   -------------

use generic_soil_ini, only : lambda,zb,dz2,zc,zd,thermal_inertia

implicit none
!-----------------------------------------------------------------------
!  arguments
!  ---------

      logical, intent(in) :: firstcall
      integer, intent(in) :: ngrid,nsoil
      real, intent(in) :: ptimestep
      real, dimension(ngrid), intent(in) :: ptsrf
      real, dimension(ngrid), intent(out) :: pcapcal,pfluxgrd
      real, dimension(ngrid,nsoil), intent(inout) :: ptsoil
      real, dimension(nsoil),intent(out) :: zout

      integer :: iflag_srf, iflag_impl


!-----------------------------------------------------------------------
!  local arrays
!  ------------

      integer :: ig,jk
      real, dimension(ngrid,nsoil+1) :: zflux
      real, dimension(ngrid) :: z1
      real, dimension(nsoil) :: za

!   local saved variables:
!   ----------------------

!-----------------------------------------------------------------------
!   Depthts:
!   --------

      iflag_srf=0

!=========== Fin des initialisations =================================

  do jk=1,nsoil
     za(jk)=dz2(jk)/ptimestep
  enddo

  iflag_impl=1

!=========== Version explicite =======================================
!  Necessite des pas de temps ultra petits. Sans doute faux

  if ( iflag_impl == 0 ) then

      do ig=1,ngrid
         zflux(ig,1)=thermal_inertia(ig)*zb(1)*(ptsoil(ig,1)-ptsrf(ig))
         zflux(ig,nsoil+1)=0.
         pfluxgrd(ig)=zflux(ig,1)
         pcapcal(ig)=thermal_inertia(ig)*zb(1)*ptimestep
         Print*,'TSOIL pcapcal',pcapcal(ig)
         print*,'TSOIL zflux',1,zflux(ig,1),za(1)
      enddo

      do jk=2,nsoil
         do ig=1,ngrid
            zflux(ig,jk)=zb(jk)*(ptsoil(ig,jk)-ptsoil(ig,jk-1))
            print*,'TSOIL zflux',jk,zflux(ig,jk),za(jk)
         enddo
      enddo

      do jk=1,nsoil
         do ig=1,ngrid
            ptsoil(ig,jk)=ptsoil(ig,jk)+(zflux(ig,jk+1)-zflux(ig,jk))*ptimestep/za(jk)
            print*,'TSOIL ptsoil',jk,ptsoil(ig,jk)
         enddo
      enddo
  
  else

!=========== Version implicite =======================================
!  Necessite des pas de temps ultra petits. Sans doute faux

      if ( .not. firstcall ) then
!-----------------------------------------------------------------------
!   Computation of the soil temperatures using the Cgrd and Dgrd
!  coefficient computed at the previous time-step:
!  -----------------------------------------------

!    surface temperature
         if (iflag_srf==0) then
            do ig=1,ngrid
               ptsoil(ig,1)=(lambda*zc(ig,2)+ptsrf(ig))/(lambda*(1.-zd(ig,2))+1.)
            enddo
         else
            do ig=1,ngrid
               ptsoil(ig,1)=zc(ig,1)+zd(ig,1)*ptsrf(ig)
            enddo
         endif

!   other temperatures
         do jk=2,nsoil
            do ig=1,ngrid
               ptsoil(ig,jk)=zc(ig,jk)+zd(ig,jk)*ptsoil(ig,jk-1)
            enddo
         enddo

!        print*,'T fond ',maxval(abs(ptsoil(:,nsoil)-200.))

      endif

!-----------------------------------------------------------------------
!   Computation of the Cgrd and Dgrd coefficient for the next step:
!   ---------------------------------------------------------------

      do ig=1,ngrid
         z1(ig)=za(nsoil)+zb(nsoil)
         zc(ig,nsoil)=za(nsoil)*ptsoil(ig,nsoil)/z1(ig)
         zd(ig,nsoil)=zb(nsoil)/z1(ig)
      enddo

      ! On descends la boucle en jk=1 mais on peut s'arrêter à
      ! jk=2 pour iflag_srf=0
      do jk=nsoil-1,1,-1
         do ig=1,ngrid
            z1(ig)=1./(za(jk)+zb(jk)+zb(jk+1)*(1.-zd(ig,jk+1)))
            zc(ig,jk)=(ptsoil(ig,jk)*za(jk)+zb(jk+1)*zc(ig,jk+1))*z1(ig)
            zd(ig,jk)=zb(jk)*z1(ig)
         enddo
      enddo

!-----------------------------------------------------------------------
!   computation of the surface diffusive flux from ground and
!   calorific capacity of the ground:
!   ---------------------------------

      do ig=1,ngrid
         pfluxgrd(ig)=thermal_inertia(ig)*zb(2)*(zc(ig,2)+(zd(ig,2)-1.)*ptsoil(ig,1))
         pcapcal(ig)=thermal_inertia(ig)*(dz2(1)+ptimestep*(1.-zd(ig,2))*zb(2))
         z1(ig)=lambda*(1.-zd(ig,2))+1.
         pcapcal(ig)=pcapcal(ig)/z1(ig)
         ! Ecriture avec za :
         !pcapcal(ig)=thermal_inertia(ig)*ptimestep*(za(1)+1.-zd(ig,2))*zb(2)) / ( lambda*(1-zd(ig,2)) + 1. )
         pfluxgrd(ig)=pfluxgrd(ig)+pcapcal(ig)*(ptsoil(ig,1)*z1(ig)-lambda*zc(ig,2)-ptsrf(ig))/ptimestep
      enddo

  endif

return
end subroutine soil_conduct

END MODULE generic_soil_conduct

