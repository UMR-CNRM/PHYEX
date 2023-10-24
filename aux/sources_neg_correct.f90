!MNH_LIC Copyright 2020-2021 CNRS, Meteo-France and Universite Paul Sabatier
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
use modd_cst,        only: xci, xcl, xcpd, xcpv, xlstt, xlvtt, xp00, xrd, xtt
use modd_nsv,        only: nsv_c2r2beg, nsv_c2r2end, nsv_lima_beg, nsv_lima_end, nsv_lima_nc, nsv_lima_nr,&
                           nsv_lima_ni, nsv_lima_ns, nsv_lima_ng, nsv_lima_nh, &
                           nsv_elecbeg, nsv_elecend
use modd_param_lima, only: lspro_lima => lspro, &
                           xctmin_lima => xctmin, xrtmin_lima => xrtmin
use modd_elec_descr, only: xrtmin_elec, xecharge                   

use mode_budget,         only: Budget_store_init, Budget_store_end
use mode_msg

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

integer :: ji, jj, jk
integer :: jr
integer :: jrmax
integer :: jsv
integer :: isv_lima_end
real, dimension(:, :, :), allocatable :: zt, zexn, zlv, zls, zcph, zcor
logical, dimension(:, :, :), allocatable :: zmask
real, dimension(:, :, :), allocatable :: zadd, zion_number !++cb--

if ( krr == 0 ) return

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
  if ( lbudget_sv .and. helec(1:3) == 'ELE' ) then
    do ji = nsv_elecbeg, nsv_elecend
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) )
    end do
  end if
else !NECON + NEGA
  if ( .not. present( prhodj ) ) &
    call Print_msg( NVERB_FATAL, 'GEN', 'Sources_neg_correct', 'optional argument prhodj not present' )

  if ( hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. &
       hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' ) then
    if ( lbudget_th) call Budget_store_init( tbudgets(NBUDGET_TH), Trim( hbudname ), prths(:, :, :)    * prhodj(:, :, :) )
    if ( lbudget_rv) call Budget_store_init( tbudgets(NBUDGET_RV), Trim( hbudname ), prrs (:, :, :, 1) * prhodj(:, :, :) )
    if ( lbudget_rc) call Budget_store_init( tbudgets(NBUDGET_RC), Trim( hbudname ), prrs (:, :, :, 2) * prhodj(:, :, :) )
    if ( lbudget_rr) call Budget_store_init( tbudgets(NBUDGET_RR), Trim( hbudname ), prrs (:, :, :, 3) * prhodj(:, :, :) )
    if ( lbudget_ri) call Budget_store_init( tbudgets(NBUDGET_RI), Trim( hbudname ), prrs (:, :, :, 4) * prhodj(:, :, :) )
    if ( lbudget_rs) call Budget_store_init( tbudgets(NBUDGET_RS), Trim( hbudname ), prrs (:, :, :, 5) * prhodj(:, :, :) )
    if ( lbudget_rg) call Budget_store_init( tbudgets(NBUDGET_RG), Trim( hbudname ), prrs (:, :, :, 6) * prhodj(:, :, :) )
    if ( lbudget_rh) call Budget_store_init( tbudgets(NBUDGET_RH), Trim( hbudname ), prrs (:, :, :, 7) * prhodj(:, :, :) )
  end if

  if ( lbudget_sv .and. ( hcloud == 'C2R2' .or. hcloud == 'KHKO' ) ) then
    do ji = nsv_c2r2beg, nsv_c2r2end
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) * prhodj(:, :, :) )
    end do
  end if
  if ( lbudget_sv .and. hcloud == 'LIMA' ) then
    do ji = nsv_lima_beg, isv_lima_end
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) * prhodj(:, :, :) )
    end do
  end if
  if ( lbudget_sv .and. helec(1:3) == 'ELE' ) then
    do ji = nsv_elecbeg, nsv_elecend
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) * prhodj(:, :, :) )
    end do
  end if
end if

allocate( zt  ( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
allocate( zexn( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
allocate( zlv ( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
allocate( zcph( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )

zexn(:, :, :) = ( ppabst(:, :, :) / xp00 ) ** (xrd / xcpd )
zt  (:, :, :) = ptht(:, :, :) * zexn(:, :, :)
zlv (:, :, :) = xlvtt + ( xcpv - xcl ) * ( zt(:, :, :) - xtt )
if ( hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. hcloud == 'LIMA' ) then
  allocate( zls( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
  zls(:, :, :) = xlstt + ( xcpv - xci ) * ( zt(:, :, :) - xtt )
end if
zcph(:, :, :) = xcpd + xcpv * prt(:, :, :, 1)

!++cb++
if ( helec(1:3) == 'ELE' ) then
  allocate( zadd( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
  allocate( zion_number( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
end if
!--cb--

deallocate( zt )

CLOUD: select case ( hcloud )
  case ( 'KESS' )
    jrmax = Size( prrs, 4 )
    do jr = 2, jrmax
      where ( prrs(:, :, :, jr) < 0. )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, jr)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, jr) * zlv(:, :, :) /  &
           ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, jr) = 0.
      end where
    end do

    where ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
      prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
      prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
           ( zcph(:, :, :) * zexn(:, :, :) )
      prrs(:, :, :, 2) = 0.
    end where


  case( 'ICE3', 'ICE4' )
    if ( hbudname == 'NETUR' ) then
      jrmax = 4
    else
      jrmax = Size( prrs, 4 )
    end if
    do jr = 4, jrmax
      if ( helec(1:3) == 'ELE' ) then
        where ( prrs(:, :, :, jr) < 0.)
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
        end where
      else
        where ( prrs(:, :, :, jr) < 0.)
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, jr)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, jr) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, jr) = 0.
        end where
      end if
    end do
!
!   cloud
    if ( hbudname == 'NETUR' ) then
      jrmax = 2
    else
      jrmax = 3
    end if
    do jr = 2, jrmax
      if ( helec(1:3) == 'ELE' ) then
        where ( prrs(:, :, :, jr) < 0.)
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
        end where
      else
        where ( prrs(:, :, :, jr) < 0.)
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, jr)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, jr) * zlv(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, jr) = 0.
        end where
      end if
    end do
!
! if rc or ri are positive, we can correct negative rv
    if ( helec(1:3) == 'ELE' ) then
!   cloud
      where ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
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
      end where
!   ice
      if ( krr > 3 ) then
        allocate( zcor( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
        where ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 4) > 0. )
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
        end where
        deallocate(zcor)
      end if
    else
!   cloud
      where ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 2) = 0.
      end where
!   ice
      if ( krr > 3 ) then
        allocate( zcor( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
        where ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 4) > 0. )
          zcor(:, :, :) = Min( -prrs(:, :, :, 1), prrs(:, :, :, 4) )
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + zcor(:, :, :)
          prths(:, :, :) = prths(:, :, :) - zcor(:, :, :) * zls(:, :, :) /  &
               ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, 4) = prrs(:, :, :, 4) - zcor(:, :, :)
        end where
        deallocate(zcor)
      end if
    end if
!
!++cb++ 08/06/23 deplace a la fin pour traiter aussi le cas de lima
! cascade the electric charge in the absence of hydrometeor
!    if ( helec(1:3) == 'ELE' ) then
!      do jr = krr, 5, -1
!        where(prrs(:,:,:,jr) < xrtmin_elec(jr))
!          prsvs(:,:,:,nsv_elecbeg-2+jr) = prsvs(:,:,:,nsv_elecbeg-2+jr) + &
!                                          prsvs(:,:,:,nsv_elecbeg-1+jr)
!          prsvs(:,:,:,nsv_elecbeg-1+jr) = 0.0
!        end where
!      end do
!      jr = 3
!      where(prrs(:,:,:,jr) < xrtmin_elec(jr))
!        prsvs(:,:,:,nsv_elecbeg-2+jr) = prsvs(:,:,:,nsv_elecbeg-2+jr) + &
!                                        prsvs(:,:,:,nsv_elecbeg-1+jr)
!        prsvs(:,:,:,nsv_elecbeg-1+jr) = 0.0
!      end where
!      do jr = 4, 2, -2
!        where(prrs(:,:,:,jr) < xrtmin_elec(jr))
!          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg-1+jr)) / xecharge
!          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg-1+jr))
!          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
!                                     zadd(:,:,:) * zion_number(:,:,:)
!          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
!                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
!          prsvs(:,:,:,nsv_elecbeg-1+jr) = 0.0
!        end where
!      end do            
!    end if
!--cb--
!
  case( 'C2R2', 'KHKO' )
    where ( prrs(:, :, :, 2) < 0. .or. prsvs(:, :, :, nsv_c2r2beg + 1) < 0. )
      prsvs(:, :, :, nsv_c2r2beg) = 0.
    end where
    do jsv = 2, 3
      where ( prrs(:, :, :, jsv) < 0. .or. prsvs(:, :, :, nsv_c2r2beg - 1 + jsv) < 0. )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, jsv)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, jsv) * zlv(:, :, :) /  &
                ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, jsv)  = 0.
        prsvs(:, :, :, nsv_c2r2beg - 1 + jsv) = 0.
      end where
    end do
    where ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
      prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
      prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
           ( zcph(:, :, :) * zexn(:, :, :) )
      prrs(:, :, :, 2) = 0.
      prsvs(:, :, :, nsv_c2r2beg + 1) = 0.
    end where
!
!
  case( 'LIMA' )
    allocate( zmask  ( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
!
! Correction where rc<0 or Nc<0
    if ( krr.GE.2 ) then
      zmask(:,:,:)=(prrs(:, :, :, 2) < xrtmin_lima(2) / ptstep)
      if (nsv_lima_nc.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_nc) < xctmin_lima(2) / ptstep )
      if ( helec == 'ELE4' ) then
        where ( zmask(:,:,:) )
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
        end where
        where ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
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
        end where
      else
        where ( zmask(:,:,:) )
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
                   ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, 2)  = 0.
        end where
        where ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 2) > 0. )
          prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 2)
          prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 2) * zlv(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
          prrs(:, :, :, 2) = 0.
        end where
      end if
      if (nsv_lima_nc.gt.0) then
         where (prrs(:, :, :, 2) == 0.)  prsvs(:, :, :, nsv_lima_nc) = 0.
      end if
    end if
!
! Correction where rr<0 or Nr<0
    if ( krr.GE.3 .and. hbudname.ne.'NETUR' ) then
      zmask(:,:,:)=(prrs(:, :, :, 3) < xrtmin_lima(3) / ptstep)
      if (nsv_lima_nr.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_nr) < xctmin_lima(3) / ptstep )
      where ( zmask(:,:,:) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 3)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 3) * zlv(:, :, :) /  &
               ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 3)  = 0.
      end where
      if ( helec == 'ELE4' ) then
        where ( zmask(:,:,:) )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+2)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+2))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+2) = 0.0
        end where
      end if
      if (nsv_lima_nr.gt.0) then
        where (prrs(:, :, :, 3) == 0.)  prsvs(:, :, :, nsv_lima_nr) = 0.
      end if
    end if
!
! Correction where ri<0 or Ni<0
    if ( krr.GE.4 ) then
      zmask(:,:,:)=(prrs(:, :, :, 4) < xrtmin_lima(4) / ptstep)
      if (nsv_lima_ni.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_ni) < xctmin_lima(4) / ptstep)
      where ( zmask(:,:,:) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 4)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 4) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 4)  = 0.
      end where
      if ( helec == 'ELE4' ) then
        where ( zmask(:,:,:) )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+3)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+3))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+3) = 0.0
        end where
        where ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 4) > 0. )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+3)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+3))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+3) = 0.0
        end where
      end if
      allocate( zcor( Size( prths, 1 ), Size( prths, 2 ), Size( prths, 3 ) ) )
      where ( prrs(:, :, :, 1) < 0. .and. prrs(:, :, :, 4) > 0. )
        zcor(:, :, :) = Min( -prrs(:, :, :, 1), prrs(:, :, :, 4) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + zcor(:, :, :)
        prths(:, :, :) = prths(:, :, :) - zcor(:, :, :) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 4) = prrs(:, :, :, 4) - zcor(:, :, :)
      end where
      deallocate( zcor )
      if (nsv_lima_ni.gt.0) then
         where (prrs(:, :, :, 4) == 0.)  prsvs(:, :, :, nsv_lima_ni) = 0.
      end if
    end if
!
! Snow     
    if ( krr.GE.5 .and. hbudname.ne.'NETUR' ) then
      zmask(:,:,:)=(prrs(:, :, :, 5) < xrtmin_lima(5) / ptstep)
      if (nsv_lima_ns.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_ns) < xctmin_lima(5) / ptstep )
      where ( zmask(:,:,:) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 5)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 5) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 5)  = 0.
      end where
      if ( helec == 'ELE4' ) then
        where ( zmask(:,:,:) )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+4)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+4))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+4) = 0.0
        end where
      end if
      if (nsv_lima_ns.gt.0) then
        where (prrs(:, :, :, 5) == 0.)  prsvs(:, :, :, nsv_lima_ns) = 0.
      end if
    end if
!
! Graupel
    if ( krr.GE.6 .and. hbudname.ne.'NETUR' ) then
      zmask(:,:,:)=(prrs(:, :, :, 6) < xrtmin_lima(6) / ptstep)
      if (nsv_lima_ng.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_ng) < xctmin_lima(6) / ptstep )
      where ( zmask(:,:,:) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 6)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 6) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 6)  = 0.
      end where
      if ( helec == 'ELE4' ) then
        where ( zmask(:,:,:) )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+5)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+5))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+5) = 0.0
        end where
      end if
      if (nsv_lima_ng.gt.0) then
        where (prrs(:, :, :, 6) == 0.)  prsvs(:, :, :, nsv_lima_ng) = 0.
      end if
    end if
!
! Hail
    if ( krr.GE.7 .and. hbudname.ne.'NETUR' ) then
      zmask(:,:,:)=(prrs(:, :, :, 7) < xrtmin_lima(7) / ptstep)
      if (nsv_lima_nh.gt.0) zmask(:,:,:)=(zmask(:,:,:) .or. prsvs(:, :, :, nsv_lima_nh) < xctmin_lima(7) / ptstep )
      where ( zmask(:,:,:) )
        prrs(:, :, :, 1) = prrs(:, :, :, 1) + prrs(:, :, :, 7)
        prths(:, :, :) = prths(:, :, :) - prrs(:, :, :, 7) * zls(:, :, :) /  &
             ( zcph(:, :, :) * zexn(:, :, :) )
        prrs(:, :, :, 7)  = 0.
      end where
      if ( helec == 'ELE4' ) then
        where ( zmask(:,:,:) )
          zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg+6)) / xecharge
          zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg+6))
          prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                    zadd(:,:,:) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                    (1. - zadd(:,:,:)) * zion_number(:,:,:)
          prsvs(:,:,:,nsv_elecbeg+6) = 0.0
        end where
      end if
      if (nsv_lima_nh.gt.0) then
        where (prrs(:, :, :, 7) == 0.)  prsvs(:, :, :, nsv_lima_nh) = 0.
      end if
   end if
!
   prsvs(:, :, :, nsv_lima_beg : isv_lima_end) = Max( 0.0, prsvs(:, :, :, nsv_lima_beg : isv_lima_end) )
   deallocate(zmask)
end select CLOUD
!
! cascade the electric charge in the absence of hydrometeor
if ( helec(1:3) == 'ELE' ) then
  do jr = krr, 5, -1
    where(prrs(:,:,:,jr) < xrtmin_elec(jr))
      prsvs(:,:,:,nsv_elecbeg-2+jr) = prsvs(:,:,:,nsv_elecbeg-2+jr) + &
                                      prsvs(:,:,:,nsv_elecbeg-1+jr)
      prsvs(:,:,:,nsv_elecbeg-1+jr) = 0.0
    end where
  end do
  jr = 3
  where(prrs(:,:,:,jr) < xrtmin_elec(jr))
    prsvs(:,:,:,nsv_elecbeg-2+jr) = prsvs(:,:,:,nsv_elecbeg-2+jr) + &
                                    prsvs(:,:,:,nsv_elecbeg-1+jr)
    prsvs(:,:,:,nsv_elecbeg-1+jr) = 0.0
  end where
  do jr = 4, 2, -2
    where(prrs(:,:,:,jr) < xrtmin_elec(jr))
      zion_number(:,:,:) = abs(prsvs(:,:,:,nsv_elecbeg-1+jr)) / xecharge
      zadd(:,:,:) = 0.5 + sign(0.5, prsvs(:,:,:,nsv_elecbeg-1+jr))
      prsvs(:,:,:,nsv_elecbeg) = prsvs(:,:,:,nsv_elecbeg) +  &
                                 zadd(:,:,:) * zion_number(:,:,:)
      prsvs(:,:,:,nsv_elecend) = prsvs(:,:,:,nsv_elecend) +  &
                                (1. - zadd(:,:,:)) * zion_number(:,:,:)
      prsvs(:,:,:,nsv_elecbeg-1+jr) = 0.0
    end where
  end do            
end if
!
if (allocated(zion_number))  deallocate( zion_number )
if (allocated(zadd))         deallocate( zadd )
if (allocated(zls))          deallocate( zls )
deallocate( zexn )
deallocate( zlv )
deallocate( zcph )


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
  if ( lbudget_sv .and. helec(1:3) == 'ELE' ) then
    do ji = nsv_elecbeg, nsv_elecend
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) )
    end do
  end if
else !NECON + NEGA
  if ( hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. &
       hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' ) then
    if ( lbudget_th) call Budget_store_end( tbudgets(NBUDGET_TH), Trim( hbudname ), prths(:, :, :)    * prhodj(:, :, :) )
    if ( lbudget_rv) call Budget_store_end( tbudgets(NBUDGET_RV), Trim( hbudname ), prrs (:, :, :, 1) * prhodj(:, :, :) )
    if ( lbudget_rc) call Budget_store_end( tbudgets(NBUDGET_RC), Trim( hbudname ), prrs (:, :, :, 2) * prhodj(:, :, :) )
    if ( lbudget_rr) call Budget_store_end( tbudgets(NBUDGET_RR), Trim( hbudname ), prrs (:, :, :, 3) * prhodj(:, :, :) )
    if ( lbudget_ri) call Budget_store_end( tbudgets(NBUDGET_RI), Trim( hbudname ), prrs (:, :, :, 4) * prhodj(:, :, :) )
    if ( lbudget_rs) call Budget_store_end( tbudgets(NBUDGET_RS), Trim( hbudname ), prrs (:, :, :, 5) * prhodj(:, :, :) )
    if ( lbudget_rg) call Budget_store_end( tbudgets(NBUDGET_RG), Trim( hbudname ), prrs (:, :, :, 6) * prhodj(:, :, :) )
    if ( lbudget_rh) call Budget_store_end( tbudgets(NBUDGET_RH), Trim( hbudname ), prrs (:, :, :, 7) * prhodj(:, :, :) )
  end if

  if ( lbudget_sv .and. ( hcloud == 'C2R2' .or. hcloud == 'KHKO' ) ) then
    do ji = nsv_c2r2beg, nsv_c2r2end
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) * prhodj(:, :, :) )
    end do
  end if
  if ( lbudget_sv .and. hcloud == 'LIMA' ) then
    do ji = nsv_lima_beg, isv_lima_end
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji) * prhodj(:, :, :) )
    end do
  end if
  if ( lbudget_sv .and. helec(1:3) == 'ELE' ) then
    do ji = nsv_elecbeg, nsv_elecend
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + ji), Trim( hbudname ), prsvs(:, :, :, ji)  * prhodj(:, :, :) )
    end do
  end if
end if

end subroutine Sources_neg_correct

end module mode_sources_neg_correct
