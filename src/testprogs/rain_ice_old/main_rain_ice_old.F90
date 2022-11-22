program main_rain_ice_old

  use xrd_getoptions, only: initoptions, getoption
  use getdata_rain_ice_old_mod, only: getdata_rain_ice_old

  use modi_rain_ice_old

  use yomhook, only: lhook, dr_hook
  use parkind1, only: jprb, jpim

  use ddh_mix, only: typ_ddh
  use yomlddh, only: tlddh
  use yommddh, only: tmddh
  use modd_dimphyex, only: dimphyex_t

  use iso_fortran_env, only: output_unit

  implicit none

  integer :: n_gp_blocks, &
             n_proma, &
             n_levels

  real, allocatable, dimension(:,:,:)   :: pdzz
  real, allocatable, dimension(:,:,:)   :: prhodj
  real, allocatable, dimension(:,:,:)   :: prhodref
  real, allocatable, dimension(:,:,:)   :: pexnref
  real, allocatable, dimension(:,:,:)   :: ppabsm
  real, allocatable, dimension(:,:,:)   :: pcit, pcit_out
  real, allocatable, dimension(:,:,:)   :: pcldfr
  real, allocatable, dimension(:,:,:)   :: ptht
  real, allocatable, dimension(:,:,:,:) :: prt
  real, allocatable, dimension(:,:,:)   :: pths, pths_out
  real, allocatable, dimension(:,:,:,:) :: prs, prs_out
  real, allocatable, dimension(:,:,:)   :: psigs
  real, allocatable, dimension(:,:)     :: psea
  real, allocatable, dimension(:,:)     :: ptown

  real, allocatable, dimension(:,:)     :: zinprc, zinprc_out
  real, allocatable, dimension(:,:)     :: pinprr, pinprr_out
  real, allocatable, dimension(:,:,:)   :: pevap, pevap_out
  real, allocatable, dimension(:,:)     :: pinprs, pinprs_out
  real, allocatable, dimension(:,:)     :: pinprg, pinprg_out
  real, allocatable, dimension(:,:,:,:) :: pfpr, pfpr_out

  real, allocatable, dimension(:,:)     :: pinprh, pinprh_out

  !Dummies for now
  !Spp stuff
  real, allocatable, dimension(:,:)   :: picenu, pkgn_acon, pkgn_sbgr
  !Ocnd2 stuff
  real, allocatable, dimension(:,:,:) :: picldfr ! Ice cloud fraction
  real, allocatable, dimension(:,:,:) :: pifr    ! Ratio cloud ice moist part to dry part
  real, allocatable, dimension(:,:,:) :: pssio   ! Super-saturation with respect to ice in the supersaturated fraction
  real, allocatable, dimension(:,:,:) :: pssiu   ! Sub-saturation with respect to ice in the subsaturated fraction

  type(dimphyex_t) :: D

  integer :: counter, c_rate
  logical :: l_verbose

  integer :: kka
  integer :: kku
  integer :: kkl
  integer :: krr
  integer :: ksplitr

  logical :: osedic
  logical :: ocnd2
  logical :: lkogan
  logical :: lmodicedep
  character(len=4) :: c_sedim
  character(len=4) :: c_micro
  character(len=4) :: csubg_aucv_rc
  logical :: owarm

  real    :: ptstep

  integer :: i

  type(typ_ddh) :: ydddh
  type(tlddh)   :: ydlddh
  type(tmddh)   :: ydmddh

  real(kind=jprb) :: time_start_real, time_end_real
  real(kind=jprb) :: time_start_cpu, time_end_cpu

  interface

    subroutine init_rain_ice_old(kulout, l_warm, c_sedim, csubg_aucv_rc)

      implicit none

      integer, intent (in)            :: kulout
      logical, intent (in)            :: l_warm
      character(4), intent (in)       :: c_sedim
      character(len=4), intent(inout) :: csubg_aucv_rc

    end subroutine init_rain_ice_old

  end interface


  n_gp_blocks = 150
  n_proma = 32
  n_levels = 90

  call initoptions()

  call getoption ("--blocks", n_gp_blocks)
  call getoption ("--nproma", n_proma)
  call getoption ("--nflevg", n_levels)
  call getoption ("--verbose", l_verbose)

  write(output_unit, *) 'n_gp_blocks: ', n_gp_blocks
  write(output_unit, *) 'n_proma:     ', n_proma
  write(output_unit, *) 'n_levels:    ', n_levels
  write(output_unit, *) 'total:       ', n_levels*n_proma*n_gp_blocks

  call getdata_rain_ice_old(n_proma, n_gp_blocks, n_levels, &
                            pdzz, prhodj, prhodref, &
                            pexnref, ppabsm, &
                            pcit, pcit_out, &
                            pcldfr, &
                            picldfr, pssio, pssiu, pifr,  &
                            ptht, prt, pths, pths_out, &
                            prs, prs_out, &
                            psigs, psea, ptown,     &
                            zinprc, zinprc_out, &
                            pinprr, pinprr_out, &
                            pevap, pevap_out,        &
                            pinprs, pinprs_out, &
                            pinprg, pinprg_out,      &
                            pinprh, pinprh_out,      &
                            picenu, pkgn_acon, pkgn_sbgr, &
                            pfpr, pfpr_out, l_verbose)

  osedic = .false.
  ocnd2 = .true.
  lkogan = .false.
  lmodicedep = .false.
  c_sedim = 'STAT'
  csubg_aucv_rc = 'PDF'
  owarm = .true.

  ptstep = 25.0000000000000

  kka = 1
  kku = n_levels
  kkl = -1
  krr = 6
  ksplitr = 2

  D%nit  = n_proma
  D%nib  = 1
  D%nie  = n_proma
  D%njt  = 1
  D%njb  = 1
  D%nje  = 1
  D%nijt = D%nit * D%njt
  D%nijb = 1
  D%nije = n_proma
  D%nkl  = -1
  D%nkt  = n_levels
  D%nka  = n_levels
  D%nku  = 1
  D%nkb  = n_levels
  D%nke  = 1
  D%nktb = 1
  D%nkte = n_levels

  call init_rain_ice_old(20, owarm, c_sedim, csubg_aucv_rc)

  call cpu_time(time_start_cpu)
  call system_clock(count=counter, count_rate=c_rate)
  time_start_real = real(counter,8)/c_rate

  do i = 1, n_gp_blocks

    call rain_ice_old(D, osedic=osedic, ocnd2=ocnd2,                                    &
                      lkogan=lkogan, lmodicedep=lmodicedep,                             &
                      hsedim=c_sedim, hsubg_aucv_rc=csubg_aucv_rc, owarm=owarm,         &
                      kka=kka, kku=kku, kkl=kkl,                                        &
                      ksplitr=ksplitr, ptstep=2*ptstep, krr=krr,                        &
                      pdzz=pdzz(:,:,i), prhodj=prhodj(:,:,i), prhodref=prhodref(:,:,i), &
                      pexnref=pexnref(:,:,i), ppabst=ppabsm(:,:,i),                     &
                      pcit=pcit(:,:,i), pcldfr=pcldfr(:,:,i),                           &
                      picldfr=picldfr(:,:,i), pssio=pssio(:,:,i), pssiu=pssiu(:,:,i),   &
                      pifr=pifr(:,:,i),                                                 &
                      ptht=ptht(:,:,i),                                                 &
                      prvt=prt(:,:,1,i), prct=prt(:,:,2,i), prrt=prt(:,:,3,i),          &
                      prit=prt(:,:,4,i), prst=prt(:,:,5,i), prgt=prt(:,:,6,i),          &
                      pths=pths(:,:,i),                                                 &
                      prvs=prs(:,:,1,i), prcs=prs(:,:,2,i), prrs=prs(:,:,3,i),          &
                      pris=prs(:,:,4,i), prss=prs(:,:,5,i), prgs=prs(:,:,6,i),          &
                      pinprc=zinprc(:,i), pinprr=pinprr(:,i), pevap3d=pevap(:,:,i),     &
                      pinprs=pinprs(:,i), pinprg=pinprg(:,i), psigs=psigs(:,:,i),       &
                      psea=psea(:,i), ptown=ptown(:,i),                                 &
                      ydddh=ydddh, ydlddh=ydlddh, ydmddh=ydmddh,                        &
                      picenu=picenu(:,i),                                               &
                      pkgn_acon=pkgn_acon(:,i), pkgn_sbgr=pkgn_sbgr(:,i),               &
                      pfpr=pfpr(:,:,:,i))

  enddo

  call cpu_time(time_end_cpu)
  call system_clock(count=counter, count_rate=c_rate)
  time_end_real = real(counter,8)/c_rate

  write (output_unit, '(a11,f8.2,a)') 'real time: ', time_end_real - time_start_real,' s'

  write (output_unit, '(a11,f8.2,a)') 'cpu time: ', time_end_cpu - time_start_cpu,' s'

end program

subroutine init_rain_ice_old(kulout, l_warm, c_sedim, csubg_aucv_rc)

  use modd_rain_ice_param, only: rain_ice_param_associate
  use modd_rain_ice_descr, only: rain_ice_descr_associate
  use modd_param_ice

  use modd_ref
  use modi_ini_rain_ice

  use modi_ini_cst
  use modd_budget

  use iso_fortran_env, only: output_unit

  implicit none

  integer, intent (in)            :: kulout
  logical, intent (in)            :: l_warm
  character(4), intent (in)       :: c_sedim
  character(len=4), intent(inout) :: csubg_aucv_rc
  character(len=4) :: c_micro

  call ini_cst

  call param_ice_associate

  call tbuconf_associate

  lbu_enable=.false.
  lbudget_u=.false.
  lbudget_v=.false.
  lbudget_w=.false.
  lbudget_th=.false.
  lbudget_tke=.false.
  lbudget_rv=.false.
  lbudget_rc=.false.
  lbudget_rr=.false.
  lbudget_ri=.false.
  lbudget_rs=.false.
  lbudget_rg=.false.
  lbudget_rh=.false.
  lbudget_sv=.false.

  ! 1. set implicit default values for modd_param_ice
  lwarm = l_warm
  cpristine_ice = 'PLAT'
  csedim = c_sedim
  csubg_aucv_rc = 'PDF'
  csubg_rc_rr_accr = 'NONE'
  csubg_rr_evap = 'NONE'
  csubg_pr_pdf = 'SIGM'
  c_micro = 'ICE3'

  ! snow riming
  csnowriming='m90 '
  xfracm90=0.1 ! fraction used for the murakami 1990 formulation

  lfeedbackt=.true. ! when .true. feed back on temperature is taken into account
  levlimit=.true.   ! when .true. water vapour pressure is limited by saturation
  lnullwetg=.true.  ! when .true. graupel wet growth is activated with null rate (to allow water shedding)
  lwetgpost=.true.  ! when .true. graupel wet growth is activated with positive temperature (to allow water shedding)
  lnullweth=.true.  ! same as lnullwetg but for hail
  lwethpost=.true.  ! same as lwetgpost but for hail
  lconvhg=.true. ! true to allow the conversion from hail to graupel
  lcrflimit=.true. !true to limit rain contact freezing to possible heat exchange
  cfrac_ice_adjust='S' ! ice/liquid partition rule to use in adjustment
  cfrac_ice_shallow_mf='S' ! ice/liquid partition rule to use in shallow_mf
  lsedim_after=.false. ! sedimentation done after microphysics
  xsplit_maxcfl=0.8
  ldeposc=.false.  ! water deposition on vegetation
  xvdeposc=0.02    ! deposition speed (2 cm.s-1)

  ! 2. set implicit default values for modd_rain_ice_descr and modd_rain_ice_param

  xthvrefz=300.

  call ini_rain_ice(kulout, c_micro)

end subroutine init_rain_ice_old

