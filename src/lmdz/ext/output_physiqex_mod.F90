MODULE output_physiqex_mod



CONTAINS 

SUBROUTINE output_physiqex(debut,zjulian,pdtphys,presnivs,paprs,u,v,t,qx,cf,zqr,zqs,zqg,ptke,theta)

      USE dimphy, only : klon,klev
      USE iophy, only : histbeg_phy,histwrite_phy
      USE ioipsl, only : histvert,histdef,histend,histsync
      USE mod_phys_lmdz_para, only : jj_nb
      USE ioipsl_getin_p_mod, ONLY : getin_p
      USE mod_grid_phy_lmdz, ONLY: nbp_lon,nbp_lat
      USE iophy, ONLY : init_iophy_new
      USE geometry_mod, ONLY: latitude_deg, longitude_deg
      USE infotrac_phy, only : nqtot



implicit none
logical, intent(in) :: debut
real, intent(in) :: pdtphys,zjulian
real,intent(in) :: presnivs(klev) ! pseudo-pressure (Pa) of mid-layers
real,intent(in) :: u(klon,klev) ! eastward zonal wind (m/s)
real,intent(in) :: v(klon,klev) ! northward meridional wind (m/s)
real,intent(in) :: t(klon,klev) ! temperature (K)
real,intent(in) :: theta(klon,klev) ! temperature (K)
real,intent(in) :: paprs(klon,klev+1) ! interlayer pressure (Pa)
real,intent(in) :: qx(klon,klev,nqtot) !tracers
real,intent(in) :: cf(klon,klev) !cloud fraction
real,intent(in) :: zqr(klon,klev) !rain specifiq content
real,intent(in) :: zqs(klon,klev) !snow specifiq content
real,intent(in) :: zqg(klon,klev) !graupel specifiq content
real,intent(in) :: ptke(klon,klev) !tke

real :: t_ops ! frequency of the IOIPSL operations (eg average over...)
real :: t_wrt ! frequency of the IOIPSL outputs
integer :: zvertid ! vertical coordinate ID
real :: dtime

integer,save :: iwrite_phys=1 ! output every iwrite_phys physics step
!$OMP THREADPRIVATE(iwrite_phys)
integer :: nhori ! horizontal coordinate ID
integer,save :: nid_hist ! output file ID
!$OMP THREADPRIVATE(nid_hist)
integer, save :: itau=0
!$OMP THREADPRIVATE(itau)

integer, save :: ioex=1


print*,'nnnnnnn ',nid_hist,debut,itau

if(debut)then

   call getin_p("iwrite_phys",iwrite_phys)

   !$OMP MASTER
   CALL iophys_ini(pdtphys)
   !$OMP END MASTER
   !$OMP BARRIER

endif


itau=itau+1

if (modulo(itau,iwrite_phys)==0) then
     call iophys_ecrit('temp',klev,'Temperature','K',t)
     call iophys_ecrit('u',klev,'zonal wind','m/s',u)
     call iophys_ecrit('v',klev,'meridinal wind','m/s',v)
     call iophys_ecrit('ps',1,'Surface pressure','Pa',paprs(:,1))
     call iophys_ecrit('qv',klev,'Water vapor specifiq content', 'kg/kg', qx(:,:,1))
     call iophys_ecrit('qc',klev,'Cloud liquid water specifiq content', 'kg/kg', qx(:,:,2))
     call iophys_ecrit('qi',klev,'Cloud solid water specifiq content', 'kg/kg', qx(:,:,3))
     call iophys_ecrit('CF',klev,'Cloud fraction', '0-1', cf)
     call iophys_ecrit('qr',klev,'Rain specifiq content', 'kg/kg', zqr)
     call iophys_ecrit('qs',klev,'Snow specifiq content', 'kg/kg', zqs)
     call iophys_ecrit('qg',klev,'Graupel specifiq content', 'kg/kg', zqg)
     call iophys_ecrit('TKE',klev,'TKE', 'm2/s2', ptke)
     call iophys_ecrit('theta',klev,'Temperature potentielle', 'K', theta)
endif


END SUBROUTINE output_physiqex
END MODULE output_physiqex_mod
