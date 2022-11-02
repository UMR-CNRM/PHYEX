PROGRAM MAIN_RAIN_ICE_OLD

use xrd_getoptions, only: initoptions, getoption
use getdata_rain_ice_old_mod, only: getdata_rain_ice_old

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use ddh_mix, only  : typ_ddh
use yomlddh, only  : tlddh
use yommddh, only  : tmddh

use iso_fortran_env, only: output_unit

implicit none

integer :: n_gp_blocks, &
           n_proma, &
           n_levels

real, allocatable, dimension(:,:,:)     :: pexnref
real, allocatable, dimension(:,:,:,:)   :: pdzz
real, allocatable, dimension(:,:,:,:)   :: prhodj
real, allocatable, dimension(:,:,:,:)   :: prhodref
real, allocatable, dimension(:,:,:)     :: pexnref2
real, allocatable, dimension(:,:,:,:)   :: ppabsm
real, allocatable, dimension(:,:,:,:)   :: pcit
real, allocatable, dimension(:,:,:,:)   :: pcldfr
real, allocatable, dimension(:,:,:,:)   :: phlc_hrc
real, allocatable, dimension(:,:,:,:)   :: phlc_hcf
real, allocatable, dimension(:,:,:,:)   :: phli_hri
real, allocatable, dimension(:,:,:,:)   :: phli_hcf
real, allocatable, dimension(:,:,:,:)   :: ptht
real, allocatable, dimension(:,:,:,:)   :: prt
real, allocatable, dimension(:,:,:,:)   :: pths
real, allocatable, dimension(:,:,:,:)   :: prs
real, allocatable, dimension(:,:,:,:)   :: psigs
real, allocatable, dimension(:,:,:)     :: psea
real, allocatable, dimension(:,:,:)     :: ptown
real, allocatable, dimension(:,:,:,:)   :: pcit_out
real, allocatable, dimension(:,:,:,:,:) :: prs_out

real, allocatable, dimension(:,:,:)     :: zinprc, zinprc_out
real, allocatable, dimension(:,:,:)     :: pinprr, pinprr_out
real, allocatable, dimension(:,:,:,:)   :: pevap, pevap_out
real, allocatable, dimension(:,:,:)     :: pinprs, pinprs_out
real, allocatable, dimension(:,:,:)     :: pinprg, pinprg_out
real, allocatable, dimension(:,:,:)     :: zindep, zindep_out
real, allocatable, dimension(:,:,:,:)   :: zrainfr, zrainfr_out
real, allocatable, dimension(:,:,:,:,:) :: pfpr, pfpr_out

logical, allocatable, dimension(:,:,:)     :: l_micro

!Dummies for now
!SPP stuff
real, allocatable, dimension(:,:)   :: picenu, pkgn_acon, pkgn_sbgr
!OCND2 stuff
real, allocatable, dimension(:,:,:) :: picldfr ! Ice cloud fraction
real, allocatable, dimension(:,:,:) :: pwcldfr ! Water or mixed-phase cloud fraction
real, allocatable, dimension(:,:,:) :: pifr    ! Ratio cloud ice moist part to dry part
real, allocatable, dimension(:,:,:) :: pssio   ! Super-saturation with respect to ice in the supersaturated fraction
real, allocatable, dimension(:,:,:) :: pssiu   ! Sub-saturation with respect to ice in the subsaturated fraction

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
logical :: csedim
logical :: csubg_aucv_rc
logical :: owarm

real    :: ptstep

integer :: i_block

type(typ_ddh) :: ydddh
type(tlddh)   :: ydlddh
type(tmddh)   :: ydmddh

real(kind=jprb) :: time_start_real, time_end_real
real(kind=jprb) :: time_start_cpu, time_end_cpu

n_gp_blocks = 150
n_proma = 32
n_levels = 90

call initoptions()

call getoption ("--blocks", n_gp_blocks)
call getoption ("--nproma", n_proma)
call getoption ("--nflevg", n_levels)
call getoption ("--verbose", l_verbose)

write(output_unit, *) n_gp_blocks
write(output_unit, *) n_proma
write(output_unit, *) n_levels

call getdata_rain_ice_old(n_proma, n_gp_blocks, n_levels, l_micro, pexnref, pdzz, prhodj, prhodref, &
                          pexnref2, ppabsm, pcit, pcldfr, phlc_hrc, phlc_hcf, phli_hri, phli_hcf, ptht, prt, pths, &
                          prs, psigs, psea, ptown, pcit_out, prs_out, zinprc, zinprc_out, pinprr, pinprr_out, pevap, pevap_out, &
                          pinprs, pinprs_out, pinprg, pinprg_out, zindep, zindep_out, zrainfr, zrainfr_out, pfpr, pfpr_out, l_verbose)

!Dummies
allocate(picenu(n_proma, 1))
allocate(pkgn_acon(n_proma, 1))
allocate(pkgn_sbgr(n_proma, 1))

allocate(picldfr(n_proma, 1, n_levels))

picenu = 1.0
pkgn_acon = 1.0
pkgn_sbgr = 1.0

osedic = .false.
ocnd2 = .true.
lkogan = .false.
lmodicedep = .false.
csedim = .false.
csubg_aucv_rc = .false.
owarm = .true.

ptstep = 25.0000000000000

kka = 1
kku = n_levels
kkl = -1
krr = 6
ksplitr = 2

call cpu_time(time_start_cpu)
call system_clock(count=counter, count_rate=c_rate)
time_start_real = real(counter,8)/c_rate

do i_block = 1, n_gp_blocks

    CALL RAIN_ICE_OLD(OSEDIC=OSEDIC, OCND2=OCND2, LKOGAN=LKOGAN, LMODICEDEP=LMODICEDEP, &
                      HSEDIM=CSEDIM, HSUBG_AUCV_RC=CSUBG_AUCV_RC, &
                      OWARM=OWARM,KKA=KKA,KKU=KKU,KKL=KKL,KSPLITR=KSPLITR, &
                      PTSTEP=2*PTSTEP, KRR=KRR, &
                      PDZZ=PDZZ, PRHODJ=PRHODJ, PRHODREF=PRHODREF, PEXNREF=PEXNREF, &
                      PPABST=PPABSM, PCIT=PCIT, PCLDFR=PCLDFR, &
                      PICLDFR=PICLDFR, PWCLDFR=PWCLDFR, &
                      PSSIO=PSSIO, PSSIU=PSSIU, PIFR=PIFR, &
                      PTHT=PTHT,PRVT=PRT(:,:,:,1),PRCT= PRT(:,:,:,2), &
                      PRRT=PRT(:,:,:,3), &
                      PRIT=PRT(:,:,:,4), PRST=PRT(:,:,:,5), &
                      PRGT=PRT(:,:,:,6), &
                      PTHS=PTHS, PRVS=PRS(:,:,:,1),PRCS=PRS(:,:,:,2), &
                      PRRS=PRS(:,:,:,3), &
                      PRIS=PRS(:,:,:,4),PRSS= PRS(:,:,:,5),PRGS= PRS(:,:,:,6), &
                      PINPRC=ZINPRC,PINPRR=PINPRR,PEVAP3D=PEVAP, &
                      PINPRS=PINPRS, PINPRG=PINPRG, &
                      PSIGS=PSIGS, PSEA=PSEA, PTOWN=PTOWN, &
                      YDDDH=YDDDH,YDLDDH=YDLDDH,YDMDDH=YDMDDH, &
                      PICENU=PICENU, PKGN_ACON=PKGN_ACON, PKGN_SBGR=PKGN_SBGR, &
                      PFPR=PFPR(:,:,:,:,i_block))

enddo

call cpu_time(time_end_cpu)
call system_clock(count=counter, count_rate=c_rate)
time_end_real = real(counter,8)/c_rate


write (output_unit, '(a11,f8.2,a)') 'real time: ', time_end_real - time_start_real,' s'

write (output_unit, '(a11,f8.2,a)') 'cpu time: ', time_end_cpu - time_start_cpu,' s'

END PROGRAM

