PROGRAM MAIN_RAIN_ICE_OLD

use xrd_getoptions
use yomhook, only : lhook, dr_hook
use parkind1, only : jprb, jpim

use iso_fortran_env, only: output_unit

implicit none

integer :: n_gp_blocks, &
           n_proma, &
           n_lev

integer :: counter, c_rate

real(kind=jprb) :: time_start_real, time_end_real
real(kind=jprb) :: time_start_cpu, time_end_cpu

n_gp_blocks = 150
n_proma = 32
n_lev = 90

call initoptions ()

call getoption ("--blocks", n_gp_blocks)
call getoption ("--nproma", n_proma)
call getoption ("--nflevg", n_lev)

write(output_unit, *) n_gp_blocks
write(output_unit, *) n_proma
write(output_unit, *) n_lev

call cpu_time(time_start_cpu)
call system_clock(count=counter, count_rate=c_rate)
time_start_real = real(counter,8)/c_rate

call cpu_time(time_end_cpu)
call system_clock(count=counter, count_rate=c_rate)
time_end_real = real(counter,8)/c_rate


write (output_unit, '(a11,f8.2,a)') 'real time: ', time_end_real - time_start_real,' s'

write (output_unit, '(a11,f8.2,a)') 'cpu time: ', time_end_cpu - time_start_cpu,' s'

END PROGRAM

