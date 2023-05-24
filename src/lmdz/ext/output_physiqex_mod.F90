MODULE output_physiqex_mod



CONTAINS 

SUBROUTINE output_physiqex(debut,zjulian,pdtphys,presnivs,paprs,u,v,t,qx)

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
real,intent(in) :: paprs(klon,klev+1) ! interlayer pressure (Pa)
real,intent(in) :: qx(klon,klev,nqtot) !tracers

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
   call getin_p("ioex",ioex)

   if ( ioex == 1 ) then
      CALL iophys_ini(pdtphys)
   else if ( ioex == 2 ) then

      CALL init_iophy_new(latitude_deg,longitude_deg)
      dtime=pdtphys
      itau=0
      call histbeg_phy("histins.nc",itau,zjulian,dtime,nhori,nid_hist)
      print*,'NNNNNNN ',nid_hist,debut
      print*,'NNNNNNN OK0'
      t_ops=pdtphys*iwrite_phys ! frequency of the IOIPSL operation
      t_wrt=pdtphys*iwrite_phys ! frequency of the outputs in the file
      print*,'NNNNNNN OK1'

   !$OMP MASTER

#ifndef CPP_IOIPSL_NO_OUTPUT 
       ! IOIPSL
       ! define vertical coordinate
       call histvert(nid_hist,"presnivs","Vertical levels","Pa",klev, &
                     presnivs,zvertid,'down')
       ! define variables which will be written in "histins.nc" file
       call histdef(nid_hist,'Temp','Atmospheric temperature','K', &
                    nbp_lon,jj_nb,nhori,klev,1,klev,zvertid,32, &
                    'inst(X)',t_ops,t_wrt)
       print*,'NNNNNNN OK2a',nid_hist,t_ops,t_wrt
       call histdef(nid_hist,'u','Eastward Zonal Wind','m/s', &
                    nbp_lon,jj_nb,nhori,klev,1,klev,zvertid,32, &
                    'inst(X)',t_ops,t_wrt)
       print*,'NNNNNNN OK2b',nid_hist,t_ops,t_wrt
       call histdef(nid_hist,'v','Northward Meridional Wind','m/s', &
                    nbp_lon,jj_nb,nhori,klev,1,klev,zvertid,32, &
                    'inst(X)',t_ops,t_wrt)
       print*,'NNNNNNN OK2c',nid_hist,t_ops,t_wrt
       call histdef(nid_hist,'ps','Surface Pressure','Pa', &
                    nbp_lon,jj_nb,nhori,1,1,1,zvertid,32, &
                    'inst(X)',t_ops,t_wrt)


       call histdef(nid_hist,'qv','Water vapor specifiq content', 'kg/kg', &
                    nbp_lon,jj_nb,nhori,klev,1,klev,zvertid,32, &
                    'inst(X)',t_ops,t_wrt)
       call histdef(nid_hist,'qc','Cloud liquid water specifiq content', 'kg/kg', &
                    nbp_lon,jj_nb,nhori,klev,1,klev,zvertid,32, &
                    'inst(X)',t_ops,t_wrt)
       call histdef(nid_hist,'qi','Cloud solid water specifiq content', 'kg/kg', &
                    nbp_lon,jj_nb,nhori,klev,1,klev,zvertid,32, &
                    'inst(X)',t_ops,t_wrt)

       ! end definition sequence
       print*,'NNNNNNN OK2',nid_hist,t_ops,t_wrt
       call histend(nid_hist)
       print*,'NNNNNNN OK3'
#endif

#ifdef CPP_XIOS
      !XIOS
          ! Declare available vertical axes to be used in output files:    
          CALL wxios_add_vaxis("presnivs", klev, presnivs)
      
          ! Declare calendar and time step
          CALL wxios_set_cal(dtime,"earth_360d",1,1,1,0.0,1,1,1,0.0)
          
          !Finalize the context:
          CALL wxios_closedef()
#endif

      !$OMP END MASTER
      !$OMP BARRIER
   endif

endif


itau=itau+1

! write some outputs:
! IOIPSL
#ifndef CPP_IOIPSL_NO_OUTPUT 
if (modulo(itau,iwrite_phys)==0) then
  if ( ioex == 1 ) then
     call iophys_ecrit('temp',klev,'Temperature','K',t)
     call iophys_ecrit('u',klev,'zonal wind','m/s',t)
     call iophys_ecrit('v',klev,'meridinal wind','m/s',t)
     call iophys_ecrit('ps',1,'Surface pressure','Pa',paprs(:,1))
     call iophys_ecrit('qv',klev,'Water vapor specifiq content', 'kg/kg', qx(:,:,1))
     call iophys_ecrit('qc',klev,'Cloud liquid water specifiq content', 'kg/kg', qx(:,:,2))
     call iophys_ecrit('qi',klev,'Cloud solid water specifiq content', 'kg/kg', qx(:,:,3))
  else if ( ioex == 2 ) then
     call histwrite_phy(nid_hist,.false.,"Temp",itau,t)
     call histwrite_phy(nid_hist,.false.,"u",itau,u)
     call histwrite_phy(nid_hist,.false.,"v",itau,v)
     call histwrite_phy(nid_hist,.false.,"ps",itau,paprs(:,1))
     call histwrite_phy(nid_hist,.false.,"qv",itau,qx(:,:,1))
     call histwrite_phy(nid_hist,.false.,"qc",itau,qx(:,:,2))
     call histwrite_phy(nid_hist,.false.,"qi",itau,qx(:,:,3))
     !$OMP MASTER
     CALL histsync(nid_hist)
     !$OMP END MASTER
  endif
endif
#endif

!XIOS
#ifdef CPP_XIOS
   !$OMP MASTER
       !Increment XIOS time
       CALL xios_update_calendar(itau)
   !$OMP END MASTER
   !$OMP BARRIER
   
       !Send fields to XIOS: (NB these fields must also be defined as
       ! <field id="..." /> in iodef.xml to be correctly used
       CALL histwrite_phy("Temp",t)
       CALL histwrite_phy("temp_newton",temp_newton)
       CALL histwrite_phy("u",u)
       CALL histwrite_phy("v",v)
       CALL histwrite_phy("ps",paprs(:,1))
       CALL histwrite_phy("qv",qx(:,:,1))
       CALL histwrite_phy("qc",qx(:,:,2))
       CALL histwrite_phy("qi",qx(:,:,3))
#endif


END SUBROUTINE output_physiqex
END MODULE output_physiqex_mod
