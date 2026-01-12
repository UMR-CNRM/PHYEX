!MNH_LIC Copyright 2020-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Author: P. Wautelet 25/06/2020 (deduplication of code from advection_metsv, resolved_cloud and turb)
! Modifications:
!  P. Wautelet 30/06/2020: remove non-local corrections in resolved_cloud for NEGA => new local corrections here
!  J. Escobar  21/07/2020: bug <-> array of size(:,:,:,0) => return if krr=0
!  P. Wautelet 10/02/2021: budgets: add missing sources for NSV_C2R2BEG+3 budget
!  C. Barthe   03/02/2022: add corrections for electric charges
!  C. Barthe   24/01/2024: add corrections for concentration of ice crystal with 
!                          different habits in LIMA
!-----------------------------------------------------------------
module mode_sources_neg_correct

implicit none

private

public :: Sources_neg_correct,Sources_neg_correct_phy

contains

subroutine Sources_neg_correct_phy(D, KSV, hcloud, helec, hbudname, KRR, ptstep, ppabst, ptht, prt, prths, prrs, prsvs, prhodj)
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),            INTENT(IN) :: D
INTEGER,                     INTENT(IN) :: KSV
character(len=*),            intent(in)           :: hcloud   ! Kind of cloud parameterization
character(len=*),            intent(in)           :: helec    ! Kind of cloud electricity parameterization
character(len=*),            intent(in)           :: hbudname ! Budget name
integer,                     intent(in)           :: KRR      ! Number of moist variables
real,                        intent(in)           :: ptstep   ! Timestep
real, dimension(D%NIT,D%NJT,D%NKT),    intent(in)           :: ppabst   ! Absolute pressure at time t
real, dimension(D%NIT,D%NJT,D%NKT),    intent(in)           :: ptht     ! Theta at time t
real, dimension(D%NIT,D%NJT,D%NKT, KRR), intent(in)         :: prt      ! Moist variables at time t
real, dimension(D%NIT,D%NJT,D%NKT),    intent(inout)        :: prths    ! Source terms
real, dimension(D%NIT,D%NJT,D%NKT, KRR), intent(inout)      :: prrs     ! Source terms
real, dimension(D%NIT,D%NJT,D%NKT, KSV), intent(inout)      :: prsvs    ! Source terms
real, dimension(D%NIT,D%NJT,D%NKT),    intent(in), optional :: prhodj   ! Dry density * jacobian
!
CALL SOURCES_NEG_CORRECT(HCLOUD, HELEC, 'NETUR',KRR,PTSTEP,PPABST,PTHT,PRT,PRTHS,PRRS,PRSVS)
!
end subroutine Sources_neg_correct_phy
!
subroutine Sources_neg_correct( hcloud, helec, hbudname, krr, ptstep, ppabst, ptht, prt, prths, prrs, prsvs, prhodj )

use modd_budget,     only: lbudget_th, lbudget_rv, lbudget_rc, lbudget_rr, lbudget_ri, &
                           lbudget_rs, lbudget_rg, lbudget_rh, lbudget_sv,             &
                           NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, &
                           NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, NBUDGET_SV1,            &
                           tbudgets
use modd_cst,        only: cst_xci => xci , cst_xcl => xcl , cst_xcpd => xcpd , cst_xcpv => xcpv , cst_xlstt => xlstt , &
                           cst_xlvtt => xlvtt , cst_xp00 => xp00 , cst_xrd => xrd, cst_xtt => xtt
use modd_elec_descr, only: xrtmin_elec, xecharge
use modd_nsv,        only: nsv_c2r2beg, nsv_c2r2end, nsv_lima_beg, nsv_lima_end, nsv_lima_nc, nsv_lima_nr, &
                           nsv_lima_ni, nsv_lima_ns, nsv_lima_ng, nsv_lima_nh,                             &
                           nsv_elecbeg, nsv_elecend
use modd_param_lima, only: lspro_lima => lspro, &
                           xctmin_lima => xctmin, xrtmin_lima => xrtmin, &
                           lcrystal_shape, nnb_crystal_shape  !++cb-- 24/01/24

use mode_budget,         only: Budget_store_init, Budget_store_end
#ifdef MNH_OPENACC
use mode_mnh_zwork,  only: Mnh_mem_get, Mnh_mem_position_pin, Mnh_mem_release
#endif
use mode_mppdb
use mode_msg

#if defined(MNH_BITREP) || defined(MNH_BITREP_OMP)
use modi_bitrep
#endif

implicit none

character(len=*),            intent(in)           :: hcloud   ! Kind of cloud parameterization
character(len=*),            intent(in)           :: helec    ! Kind of cloud electricity parameterization
character(len=*),            intent(in)           :: hbudname ! Budget name
integer,                     intent(in)           :: krr      ! Number of moist variables
real,                        intent(in)           :: ptstep   ! Timestep
real, dimension(:, :, :),    intent(in)           :: ppabst   ! Absolute pressure at time t
real, dimension(:, :, :),    intent(in)           :: ptht     ! Theta at time t
real, dimension(:, :, :, :), intent(in)           :: prt      ! Moist variables at time t
real, dimension(:, :, :),    intent(inout)        :: prths    ! Source terms
real, dimension(:, :, :, :), intent(inout)        :: prrs     ! Source terms
real, dimension(:, :, :, :), intent(inout)        :: prsvs    ! Source terms
real, dimension(:, :, :),    intent(in), optional :: prhodj   ! Dry density * jacobian

integer :: jiu, jju, jku
integer :: ji, jj, jk
integer :: jr
integer :: jrmax
integer :: jsv
integer :: isv_lima_end
integer :: jsh  ! loop index for ice crystal shapes  !++cb--
real, dimension(:, :, :), pointer, contiguous :: zt, zexn, zlv, zls, zcph, zcor
REAL, DIMENSION(:,:,:), pointer , contiguous :: ZTEMP_BUD
logical, dimension(:, :, :), pointer, contiguous :: zmask
real, dimension(:, :, :), pointer, contiguous :: zadd, zion_number !++cb--
real, dimension(size(prths,1),size(prths,2),size(prths,3)) :: zni_tot ! total concentration of ice crystals !++cb--
logical :: GELEC_ELE , GELEC_ELE4

REAL :: xci, xcl, xcpd, xcpv, xlstt, xlvtt, xp00, xrd, xtt

if ( krr == 0 ) return

GELEC_ELE  = ( helec(1:3) == 'ELE' )
GELEC_ELE4 = ( helec == 'ELE4' )

zcor => null()

xci = cst_xci
xcl = cst_xcl
xcpd = cst_xcpd
xcpv = cst_xcpv
xlstt = cst_xlstt
xlvtt = cst_xlvtt
xp00 = cst_xp00
xrd = cst_xrd
xtt = cst_xtt

jiu = Size(prths, 1 )
jju = Size(prths, 2 )
jku = Size(prths, 3 )

#ifndef MNH_OPENACC
allocate( zt  ( jiu, jju, jku ) )
allocate( zexn( jiu, jju, jku ) )
allocate( zlv ( jiu, jju, jku ) )
allocate( zcph( jiu, jju, jku ) )
if ( hcloud == 'LIMA') allocate (zmask(jiu, jju, jku ) )
if ( hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. hcloud == 'LIMA' ) then
  allocate( zls( jiu, jju, jku ) )
  if ( krr > 3 ) then
    allocate( zcor( jiu, jju, jku ) )
  end if
else
  allocate( zls(0, 0, 0) )
end if
if ( .not. Associated( zcor ) ) Allocate( zcor(0, 0, 0) )
ALLOCATE(ZTEMP_BUD( jiu, jju, jku ))
#else
!Pin positions in the pools of MNH memory
call Mnh_mem_position_pin( 'Sources_neg_correct' )

call Mnh_mem_get( zt,   jiu, jju, jku )
call Mnh_mem_get( zexn, jiu, jju, jku )
call Mnh_mem_get( zlv,  jiu, jju, jku )
call Mnh_mem_get( zcph, jiu, jju, jku )
if ( hcloud == 'LIMA') call Mnh_mem_get( zmask, jiu, jju, jku )
if ( hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. hcloud == 'LIMA' ) then
  call Mnh_mem_get( zls, jiu, jju, jku )
  if ( krr > 3 ) then
    call Mnh_mem_get( zcor, jiu, jju, jku )
  else
    call Mnh_mem_get( zcor, 0, 0, 0 )
  end if
else
  call Mnh_mem_get( zls,  0, 0, 0 )
  call Mnh_mem_get( zcor, 0, 0, 0 )
end if
call MNH_MEM_GET( ZTEMP_BUD, jiu, jju, jku )
if ( GELEC_ELE ) then
  call Mnh_mem_get( zadd, jiu, jju, jku )
  call Mnh_mem_get( zion_number, jiu, jju, jku )
else
  call Mnh_mem_get( zadd, 0, 0, 0 )
  call Mnh_mem_get( zion_number, 0, 0, 0 )
end if
#endif

!$acc data present( ppabst, ptht, prt, prths, prrs, prsvs, prhodj )

if ( mppdb_initialized ) then
  !Check all IN arrays
  call Mppdb_check( ppabst, "Sources_neg_correct beg:ppabst")
  call Mppdb_check( ptht,   "Sources_neg_correct beg:ptht")
  call Mppdb_check( prt,    "Sources_neg_correct beg:prt")
  if ( Present( prhodj ) ) call Mppdb_check( prhodj, "Sources_neg_correct beg:prhodj")
  !Check all INOUT arrays
  call Mppdb_check( prths, "Sources_neg_correct beg:prths")
  call Mppdb_check( prrs,  "Sources_neg_correct beg:prrs")
  call Mppdb_check( prsvs, "Sources_neg_correct beg:prsvs")
end if

if ( hbudname /= 'NEADV' .and. hbudname /= 'NECON' .and. hbudname /= 'NEGA' .and. hbudname /= 'NETUR' ) &
  call Print_msg( NVERB_WARNING, 'GEN', 'Sources_neg_correct', 'budget '//hbudname//' not yet tested' )

if ( hcloud == 'LIMA' ) then
  ! The negativity correction does not apply to the SPRO (supersaturation) variable which may be naturally negative
  if ( lspro_lima ) then
    isv_lima_end = nsv_lima_end - 1
  else
    isv_lima_end = nsv_lima_end
  end if
end if

if ( hbudname /= 'NECON' .and. hbudname /= 'NEGA' ) then
  if ( hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. &
       hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), Trim( hbudname ), prths(:, :, :) )
    if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), Trim( hbudname ), prrs (:, :, :, 1) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), Trim( hbudname ), prrs (:, :, :, 2) )
    if ( lbudget_rr .and.                                                                                   &
       (   hbudname /= 'NETUR' .or.                                                                         &
         ( hbudname == 'NETUR' .and. ( hcloud == 'C2R2' .or. hcloud == 'KHKO' .or. hcloud == 'LIMA' ) ) ) ) &
                    call Budget_store_init( tbudgets(NBUDGET_RR), Trim( hbudname ), prrs (:, :, :, 3) )
        IF (lbudget_ri .and.                                                                                    &
       (   hbudname /= 'NETUR' .or.                                                                         &
         ( hbudname == 'NETUR' .and. ( hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. hcloud == 'LIMA' ) ) ) ) &
                    call Budget_store_init( tbudgets(NBUDGET_RI), Trim( hbudname ), prrs (:, :, :, 4) )
    if ( lbudget_rs .and. hbudname /= 'NETUR' ) call Budget_store_init( tbudgets(NBUDGET_RS), Trim( hbudname ), prrs (:, :, :, 5) )
    if ( lbudget_rg .and. hbudname /= 'NETUR' ) call Budget_store_init( tbudgets(NBUDGET_RG), Trim( hbudname ), prrs (:, :, :, 6) )
    if ( lbudget_rh .and. hbudname /= 'NETUR' ) call Budget_store_init( tbudgets(NBUDGET_RH), Trim( hbudname ), prrs (:, :, :, 7) )
  end if

  if ( lbudget_sv .and. ( hcloud == 'C2R2' .or. hcloud == 'KHKO' ) ) then
    do ji = nsv_c2r2beg, nsv_c2r2end
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) )
    end do
  end if
  if ( lbudget_sv .and. hcloud == 'LIMA' ) then
    do ji = nsv_lima_beg, isv_lima_end
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) )
    end do
  end if
  if ( lbudget_sv .and. GELEC_ELE ) then
    do ji = nsv_elecbeg, nsv_elecend
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) )
    end do
  end if
else !NECON + NEGA
  if ( .not. present( prhodj ) ) &
    call Print_msg( NVERB_FATAL, 'GEN', 'Sources_neg_correct', 'optional argument prhodj not present' )

  if ( hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. &
       hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' ) then
     if ( lbudget_th) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) = prths(:, :, :)    * prhodj(:, :, :)
        !$acc end kernels        
        call Budget_store_init( tbudgets(NBUDGET_TH), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rv) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 1)  * prhodj(:, :, :)
        !$acc end kernels                
        call Budget_store_init( tbudgets(NBUDGET_RV), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rc) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 2)  * prhodj(:, :, :)
        !$acc end kernels                
        call Budget_store_init( tbudgets(NBUDGET_RC), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rr) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 3)  * prhodj(:, :, :)
        !$acc end kernels                        
        call Budget_store_init( tbudgets(NBUDGET_RR), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_ri) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 4)  * prhodj(:, :, :)
        !$acc end kernels                                
        call Budget_store_init( tbudgets(NBUDGET_RI), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rs) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 5)  * prhodj(:, :, :)
        !$acc end kernels         
        call Budget_store_init( tbudgets(NBUDGET_RS), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rg) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 6)  * prhodj(:, :, :)
        !$acc end kernels                 
        call Budget_store_init( tbudgets(NBUDGET_RG), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rh) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 7)  * prhodj(:, :, :)
        !$acc end kernels                         
        call Budget_store_init( tbudgets(NBUDGET_RH), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
  end if

  if ( lbudget_sv .and. ( hcloud == 'C2R2' .or. hcloud == 'KHKO' ) ) then
    do ji = nsv_c2r2beg, nsv_c2r2end
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) = prsvs(:, :, :, ji) * prhodj(:, :, :)
        !$acc end kernels        
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
    end do
  end if
  if ( lbudget_sv .and. hcloud == 'LIMA' ) then
    do ji = nsv_lima_beg, isv_lima_end
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) = prsvs(:, :, :, ji) * prhodj(:, :, :)
        !$acc end kernels        
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
    end do
  end if
  if ( lbudget_sv .and. GELEC_ELE ) then
    do ji = nsv_elecbeg, nsv_elecend
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) * prhodj(:, :, :) )
    end do
  end if
end if

!$acc data present( zt, zexn, zlv, zcph, zls, zcor )
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
!$acc kernels present_cr(zexn,zt,zlv)
zexn(:, :, :) = ( ppabst(:, :, :) / xp00 ) ** (xrd / xcpd )
!$acc end kernels
#else
!$acc kernels
!$acc_nv loop collapse(3) independent
!$acc_cr loop independent
DO CONCURRENT(jk=1:jku,jj=1:jju,ji=1:jiu)
!DO CONCURRENT(jk=1:1,jj=1:1,ji=1:jiu*jju*jku)
 zexn(ji,jj,jk) = Br_pow( ppabst(ji,jj,jk) / xp00,  xrd / xcpd )
ENDDO
!$acc end kernels
#endif
!$acc kernels present_cr(zexn,zt,zlv)
zt  (:, :, :) = ptht(:, :, :) * zexn(:, :, :)
zlv (:, :, :) = xlvtt + ( xcpv - xcl ) * ( zt(:, :, :) - xtt )
!$acc end kernels
if ( hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. hcloud == 'LIMA' ) then
!$acc kernels present_cr(zls)
  zls(:, :, :) = xlstt + ( xcpv - xci ) * ( zt(:, :, :) - xtt )
!$acc end kernels
end if
!$acc kernels
zcph(:, :, :) = xcpd + xcpv * prt(:, :, :, 1)
!$acc end kernels

#ifndef MNH_OPENACC
!++cb++
if ( GELEC_ELE ) then
  allocate( zadd( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
  allocate( zion_number( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
end if

deallocate( zt )
#endif

CLOUD: select case ( hcloud )
  case ( 'KESS' )
    jrmax = Size( prrs, 4 )
    do jr = 2, jrmax
!PW: kernels directive inside do loop on jr because compiler bug... (NVHPC 21.7)
!$acc kernels present_cr(zexn,zcph,zlv)
      !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
      WHERE ( prrs(:, :, :, jr) < 0. )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, jr)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, jr) * zlv(:, :, :) /  &
           ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, jr) = 0.
      END WHERE
      !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
!$acc end kernels
    end do

!$acc kernels present_cr(zexn,zcph,zlv)
    !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
    WHERE ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
      prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
      prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
           ( zcph(:, :, :) * zexn(:, :, :) )
      prrs(:, :, :, 2) = 0.
    END WHERE
    !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
!$acc end kernels


  case( 'ICE3', 'ICE4' )
    if ( hbudname == 'NETUR' ) then
      jrmax = 4
    else
      jrmax = Size( prrs, 4 )
    end if
!$acc kernels
    do jr = 4, jrmax
      if ( GELEC_ELE ) then
      !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
        WHERE ( prrs(:, :, :, jr) < 0.)
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, jr)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, jr) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, jr) = 0.
          !
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+jr-1)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+jr-1))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+jr-1) = 0.0
        END WHERE
      !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
      else
      !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
        WHERE ( prrs(:, :, :, jr) < 0.)
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, jr)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, jr) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, jr) = 0.
        END WHERE
      !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
      end if
    end do
!$acc end kernels
!
!   cloud
    if ( hbudname == 'NETUR' ) then
      jrmax = 2
    else
      jrmax = 3
    end if
!$acc kernels
    do jr = 2, jrmax
      if ( GELEC_ELE ) then
        !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
        WHERE ( prrs(:, :, :, jr) < 0.)
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, jr)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, jr) * zlv(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, jr) = 0.
          !
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+jr-1)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+jr-1))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+jr-1) = 0.0
        END WHERE
        !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
      else
        !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
        WHERE ( prrs(:, :, :, jr) < 0.)
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, jr)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, jr) * zlv(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, jr) = 0.
        END WHERE
        !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
      end if
    end do
!$acc end kernels
!
! if rc or ri are positive, we can correct negative rv
    if ( GELEC_ELE ) then
!   cloud
!$acc kernels
      !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
      WHERE ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 2) = 0.
        !
        zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+1)) / xecharge
        zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+1))
        prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                   zadd(:,:,:) * zion_number(:,:,:)
        prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                  (1. - zadd(:,:,:)) * zion_number(:,:,:)
        prsvs(:,:,:,nsv_elecbeg+1) = 0.0   
      END WHERE
    !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
!   ice
      if ( krr > 3 ) then
        !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
        WHERE ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 4) > 0. )
          zcor(:, :, :) = Min( -prrs(:, :, :, 1), prrs(:, :, :, 4) )
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + zcor(:, :, :)
          prths(:, :, :) = prths(:, :, :) - zcor(:, :, :) * zls(:, :, :) /  &
               ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, 4) = prrs(:, :, :, 4) - zcor(:, :, :)
          !
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+3)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+3))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                     zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+3) = 0.0
        END WHERE
        !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
      end if
!$acc end kernels
    else
!$acc kernels
!   cloud
      !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
      WHERE ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 2) = 0.
      END WHERE
      !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
!   ice
      if ( krr > 3 ) then
        !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
        WHERE ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 4) > 0. )
          zcor(:, :, :) = Min( -prrs(:, :, :, 1), prrs(:, :, :, 4) )
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + zcor(:, :, :)
          prths(:, :, :) = prths(:, :, :) - zcor(:, :, :) * zls(:, :, :) /  &
               ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, 4) = prrs(:, :, :, 4) - zcor(:, :, :)
        END WHERE
        !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
      end if
!$acc end kernels
    end if
!
!
  case( 'C2R2', 'KHKO' )
!$acc kernels
    !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
    WHERE ( prrs(:, :, :, 2) < 0. .or. prsvs(:, :, :, nsv_c2r2beg + 1) < 0. )
      prsvs(:, :, :, nsv_c2r2beg) = 0.
    END WHERE
    !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
!$acc end kernels
    do jsv = 2, 3
!$acc kernels present_cr(zexn,zcph,zlv)
!PW: kernels directive inside do loop on jr because compiler bug... (NVHPC 21.7)
      !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
      WHERE ( prrs(:, :, :, jsv) < 0. .or. prsvs(:, :, :, nsv_c2r2beg - 1 + jsv) < 0. )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, jsv)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, jsv) * zlv(:, :, :) /  &
                ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, jsv)  = 0.
        prsvs(:, :, :, nsv_c2r2beg - 1 + jsv) = 0.
      END WHERE
      !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
!$acc end kernels
    end do
!$acc kernels present_cr(zexn,zcph,zlv)
    !$mnh_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
    WHERE ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
      prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
      prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
           ( zcph(:, :, :) * zexn(:, :, :) )
      prrs(:, :, :, 2) = 0.
      prsvs(:, :, :, nsv_c2r2beg + 1) = 0.
    END WHERE
    !$mnh_end_expand_where(ji = 1 : jiu, jj = 1 : jju, jk = 1 : jku)
!$acc end kernels
!
!
  case( 'LIMA' )
! Correction WHERE rc<0 or Nc<0
    if ( krr.GE.2 ) then
      zmask(:,:,:)=(prrs(:, :, :, 2) < xrtmin_lima(2) / ptstep)
      if (nsv_lima_nc.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_nc) < xctmin_lima(2) / ptstep )
      if ( GELEC_ELE4 ) then
        WHERE ( zmask(:,:,:) )
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
                   ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, 2)  = 0.
          !
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+1)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+1))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+1) = 0.0
        END WHERE
        WHERE ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, 2) = 0.
          !
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+1)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+1))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+1) = 0.0
        END WHERE
      else
        WHERE ( zmask(:,:,:) )
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
                   ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, 2)  = 0.
        END WHERE
        WHERE ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, 2) = 0.
        END WHERE
      end if
      if (nsv_lima_nc.gt.0) then
         WHERE (prrs(:, :, :, 2) == 0.)  prsvs(:, :, :, nsv_lima_nc) = 0.
      end if
    end if
!
! Correction WHERE rr<0 or Nr<0
    if ( krr.GE.3 .and. hbudname.ne.'NETUR' ) then
      zmask(:,:,:)=(prrs(:, :, :, 3) < xrtmin_lima(3) / ptstep)
      if (nsv_lima_nr.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_nr) < xctmin_lima(3) / ptstep )
      WHERE ( zmask(:,:,:) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 3)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 3) * zlv(:, :, :) /  &
               ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 3)  = 0.
      END WHERE
      if ( GELEC_ELE4 ) then
        WHERE ( zmask(:,:,:) )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+2)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+2))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+2) = 0.0
        END WHERE
      end if
      if (nsv_lima_nr.gt.0) then
        WHERE (prrs(:, :, :, 3) == 0.)  prsvs(:, :, :, nsv_lima_nr) = 0.
      end if
    end if
!
! Correction WHERE ri<0 or Ni<0
    if ( krr.GE.4 ) then
      zmask(:,:,:)=(prrs(:, :, :, 4) < xrtmin_lima(4) / ptstep)
!++cb++ 24/01/24
      zni_tot(:, :, :) = 0.
      if (lcrystal_shape) then
        ! compute the total ice crystal concentration
        do jsh = 1, nnb_crystal_shape
          zni_tot(:, :, :) = zni_tot(:, :, :) + prsvs(:, :, :, nsv_lima_ni+jsh-1)
        end do
      else
        zni_tot(:, :, :) = prsvs(:, :,: , nsv_lima_ni)
      end if 
!      if (nsv_lima_ni.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_ni) < xctmin_lima(4) / ptstep)
      if (nsv_lima_ni.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. zni_tot(:, :, :) < xctmin_lima(4) / ptstep)
!--cb--
      WHERE ( zmask(:,:,:) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 4)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 4) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 4)  = 0.
      END WHERE
      if ( GELEC_ELE4 ) then
        WHERE ( zmask(:,:,:) )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+3)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+3))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+3) = 0.0
        END WHERE
        WHERE ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 4) > 0. )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+3)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+3))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+3) = 0.0
        END WHERE
      end if
      WHERE ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 4) > 0. )
        zcor(:, :, :) = Min( -prrs(:, :, :, 1), prrs(:, :, :, 4) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + zcor(:, :, :)
        prths(:, :, :) = prths(:, :, :) - zcor(:, :, :) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 4) = prrs(:, :, :, 4) - zcor(:, :, :)
      END WHERE
      if (nsv_lima_ni.gt.0) then
!++cb++ 24/01/24
!         WHERE (prrs(:, :, :, 4) == 0.)  prsvs(:, :, :, nsv_lima_ni) = 0.
        if (.not. lcrystal_shape) then
          WHERE (prrs(:, :, :, 4) == 0.)  prsvs(:, :, :, nsv_lima_ni) = 0.
        else
          do jsh = 1, nnb_crystal_shape
            WHERE (prrs(:, :, :, 4) == 0.)  prsvs(:, :, :, nsv_lima_ni+jsh-1) = 0.
          end do
        end if
!--cb--
      end if
    end if
!
! Snow     
    if ( krr.GE.5 .and. hbudname.ne.'NETUR' ) then
      zmask(:,:,:)=(prrs(:, :, :, 5) < xrtmin_lima(5) / ptstep)
      if (nsv_lima_ns.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_ns) < xctmin_lima(5) / ptstep )
      WHERE ( zmask(:,:,:) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 5)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 5) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 5)  = 0.
      END WHERE
      if ( GELEC_ELE4 ) then
        WHERE ( zmask(:,:,:) )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+4)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+4))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+4) = 0.0
        END WHERE
      end if
      if (nsv_lima_ns.gt.0) then
        WHERE (prrs(:, :, :, 5) == 0.)  prsvs(:, :, :, nsv_lima_ns) = 0.
      end if
    end if
!
! Graupel
    if ( krr.GE.6 .and. hbudname.ne.'NETUR' ) then
      zmask(:,:,:)=(prrs(:, :, :, 6) < xrtmin_lima(6) / ptstep)
      if (nsv_lima_ng.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_ng) < xctmin_lima(6) / ptstep )
      WHERE ( zmask(:,:,:) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 6)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 6) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 6)  = 0.
      END WHERE
      if ( GELEC_ELE4 ) then
        WHERE ( zmask(:,:,:) )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+5)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+5))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+5) = 0.0
        END WHERE
      end if
      if (nsv_lima_ng.gt.0) then
        WHERE (prrs(:, :, :, 6) == 0.)  prsvs(:, :, :, nsv_lima_ng) = 0.
      end if
    end if
!
! Hail
    if ( krr.GE.7 .and. hbudname.ne.'NETUR' ) then
      zmask(:,:,:)=(prrs(:, :, :, 7) < xrtmin_lima(7) / ptstep)
      if (nsv_lima_nh.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_nh) < xctmin_lima(7) / ptstep )
      WHERE ( zmask(:,:,:) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 7)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 7) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 7)  = 0.
      END WHERE
      if ( GELEC_ELE4 ) then
        WHERE ( zmask(:,:,:) )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+6)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+6))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+6) = 0.0
        END WHERE
      end if
      if (nsv_lima_nh.gt.0) then
        WHERE (prrs(:, :, :, 7) == 0.)  prsvs(:, :, :, nsv_lima_nh) = 0.
      end if
   end if
!
   prsvs(:, :, :, nsv_lima_beg : isv_lima_end) = Max( 0.0, prsvs(:, :, :, nsv_lima_beg : isv_lima_end) )
end select CLOUD
!
! cascade the electric charge in the absence of hydrometeor
if ( GELEC_ELE ) then
  do jr = krr, 5, -1
    WHERE(prrs(:,:,:,jr) < xrtmin_elec(jr))
      prsvs(:,:,:,nsv_elecbeg-2+jr) = prsvs(:,:,:,nsv_elecbeg-2+jr) + &
                                      prsvs(:,:,:,nsv_elecbeg-1+jr)
      prsvs(:,:,:,nsv_elecbeg-1+jr) = 0.0
    END WHERE
  end do
  jr = 3
  WHERE(prrs(:,:,:,jr) < xrtmin_elec(jr))
    prsvs(:,:,:,nsv_elecbeg-2+jr) = prsvs(:,:,:,nsv_elecbeg-2+jr) + &
                                    prsvs(:,:,:,nsv_elecbeg-1+jr)
    prsvs(:,:,:,nsv_elecbeg-1+jr) = 0.0
  END WHERE
  do jr = 4, 2, -2
    WHERE(prrs(:,:,:,jr) < xrtmin_elec(jr))
      zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg-1+jr)) / xecharge
      zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg-1+jr))
      prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                 zadd(:,:,:) * zion_number(:,:,:)
      prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                (1. - zadd(:,:,:)) * zion_number(:,:,:)
      prsvs(:,:,:,nsv_elecbeg-1+jr) = 0.0
    END WHERE
  end do            
end if
!

if ( hbudname /= 'NECON' .and. hbudname /= 'NEGA' ) then
  if ( hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. &
       hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), Trim( hbudname ), prths(:, :, :) )
    if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), Trim( hbudname ), prrs (:, :, :, 1) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), Trim( hbudname ), prrs (:, :, :, 2) )
    if ( lbudget_rr .and.                                                                                   &
       (   hbudname /= 'NETUR' .or.                                                                         &
         ( hbudname == 'NETUR' .and. ( hcloud == 'C2R2' .or. hcloud == 'KHKO' .or. hcloud == 'LIMA' ) ) ) ) &
                    call Budget_store_end( tbudgets(NBUDGET_RR), Trim( hbudname ), prrs (:, :, :, 3) )
        IF (lbudget_ri .and.                                                                                    &
       (   hbudname /= 'NETUR' .or.                                                                         &
         ( hbudname == 'NETUR' .and. ( hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. hcloud == 'LIMA' ) ) ) ) &
                    call Budget_store_end( tbudgets(NBUDGET_RI), Trim( hbudname ), prrs (:, :, :, 4) )
    if ( lbudget_rs .and. hbudname /= 'NETUR' ) call Budget_store_end( tbudgets(NBUDGET_RS), Trim( hbudname ), prrs (:, :, :, 5) )
    if ( lbudget_rg .and. hbudname /= 'NETUR' ) call Budget_store_end( tbudgets(NBUDGET_RG), Trim( hbudname ), prrs (:, :, :, 6) )
    if ( lbudget_rh .and. hbudname /= 'NETUR' ) call Budget_store_end( tbudgets(NBUDGET_RH), Trim( hbudname ), prrs (:, :, :, 7) )
  end if

  if ( lbudget_sv .and. ( hcloud == 'C2R2' .or. hcloud == 'KHKO' ) ) then
    do ji = nsv_c2r2beg, nsv_c2r2end
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) )
    end do
  end if
  if ( lbudget_sv .and. hcloud == 'LIMA' ) then
    do ji = nsv_lima_beg, isv_lima_end
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) )
    end do
  end if
  if ( lbudget_sv .and. GELEC_ELE ) then
    do ji = nsv_elecbeg, nsv_elecend
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) )
    end do
  end if
else !NECON + NEGA
  if ( hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. &
       hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' ) then
     if ( lbudget_th) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) = prths(:, :, :)    * prhodj(:, :, :)
        !$acc end kernels
        call Budget_store_end( tbudgets(NBUDGET_TH), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rv) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 1)  * prhodj(:, :, :)
        !$acc end kernels        
        call Budget_store_end( tbudgets(NBUDGET_RV), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rc) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 2)  * prhodj(:, :, :)
        !$acc end kernels                
        call Budget_store_end( tbudgets(NBUDGET_RC), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rr) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 3)  * prhodj(:, :, :)
        !$acc end kernels                        
        call Budget_store_end( tbudgets(NBUDGET_RR), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_ri) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 4)  * prhodj(:, :, :)
        !$acc end kernels                                
        call Budget_store_end( tbudgets(NBUDGET_RI), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rs) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 5)  * prhodj(:, :, :)
        !$acc end kernels         
        call Budget_store_end( tbudgets(NBUDGET_RS), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rg) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 6)  * prhodj(:, :, :)
        !$acc end kernels                 
        call Budget_store_end( tbudgets(NBUDGET_RG), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
     if ( lbudget_rh) then
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) =  prrs (:, :, :, 7)  * prhodj(:, :, :)
        !$acc end kernels                         
        call Budget_store_end( tbudgets(NBUDGET_RH), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
     end if
  end if

  if ( lbudget_sv .and. ( hcloud == 'C2R2' .or. hcloud == 'KHKO' ) ) then
     do ji = nsv_c2r2beg, nsv_c2r2end
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) = prsvs(:, :, :, ji) * prhodj(:, :, :)
        !$acc end kernels        
        call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), ZTEMP_BUD(:,:,:) )
    end do
  end if
  if ( lbudget_sv .and. hcloud == 'LIMA' ) then
     do ji = nsv_lima_beg, isv_lima_end
        !$acc kernels present(ZTEMP_BUD)
        ZTEMP_BUD(:,:,:) = prsvs(:, :, :, ji) * prhodj(:, :, :)
        !$acc end kernels                
        call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ),  ZTEMP_BUD(:,:,:) )
    end do
  end if
  if ( lbudget_sv .and. GELEC_ELE ) then
    do ji = nsv_elecbeg, nsv_elecend
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji)  * prhodj(:, :, :) )
    end do
  end if
end if


#ifndef MNH_OPENACC
deallocate( zexn, zlv, zcph, zls, ZTEMP_BUD )
if ( hcloud == 'LIMA' ) deallocate( zmask )
IF ( hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. hcloud == 'LIMA' .and. KRR > 3) then
 DEALLOCATE(zcor)
END IF
IF ( GELEC_ELE ) THEN
  DEALLOCATE(zadd)
  DEALLOCATE(zion_number)
END IF
#else
!$acc end data

!Release all memory allocated with Mnh_mem_get calls since last call to Mnh_mem_position_pin
call Mnh_mem_release( 'Sources_neg_correct' )
#endif

if ( mppdb_initialized ) then
  !Check all INOUT arrays
  call Mppdb_check( prths, "Sources_neg_correct end:prths")
  call Mppdb_check( prrs,  "Sources_neg_correct end:prrs")
  call Mppdb_check( prsvs, "Sources_neg_correct end:prsvs")
end if

!$acc end data

end subroutine Sources_neg_correct

end module mode_sources_neg_correct
