!MNH_LIC Copyright 2019-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
module modi_tools
  ! /!\ converte this module MODE -> in MODI + subroutines/functions outside module
  ! /!\ to avoid very long compilation time if implementation change in subroutines
  ! /!\ This 'empty' module is here to avoid 'automatic' generation of wrong interface modi_quisort + modi_upcase
  ! /!\ in futur version rename all 'use mode_tools' -> 'use modi_tools'
  ! /!\ and change the module name below  mode_tools -> modi_tools
end module modi_tools
!################
module mode_tools
!################
!
!    Purpose
!    -------
!
!     The Purpose of this module is to provide useful tools for MesoNH
!
!    Author
!    ------
!     P. Wautelet 14/02/2019
!
! Modifications:
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet 05/06/2019: add Countjv_device
!  P. Wautelet 20/06/2019: add Countjv1d, Countjv1d_device and Countjv2d_device subroutines
!  P. Wautelet 17/01/2020: move Quicksort to tools.f90
!  J. Escobar  25/08/2021: nvhpc21.X bug on 'atomic' in host -> switch to version without atomic on HOST    

implicit none

private

public :: Quicksort
interface
   recursive subroutine Quicksort( ka, kbeg, kend, kpos )
     integer, dimension(:),           intent(inout) :: ka
     integer,                         intent(in)    :: kbeg, kend
     integer, dimension(:), optional, intent(inout) :: kpos
   end subroutine Quicksort
end interface

public :: Upcase
interface
   function Upcase(hstring)
     character(len=*), intent(in) :: hstring
     character(len=len(hstring))  :: upcase
   end function Upcase
end interface

public :: Countjv
interface Countjv
   function Countjv1d(ltab,i1) result(ic)
     logical, dimension(:), intent(in)  :: ltab ! Mask
     integer, dimension(:), intent(out) :: i1   ! Positions of elements with 'true' value
     integer                            :: ic   ! Total number of 'true' values
   end function Countjv1d
   function Countjv2d(ltab,i1,i2) result(ic)
     logical, dimension(:,:), intent(in)  :: ltab   ! Mask
     integer, dimension(:),   intent(out) :: i1, i2 ! Positions of elements with 'true' value
     integer                              :: ic     ! Total number of 'true' values
   end function Countjv2d
   function Countjv3d(ltab,i1,i2,i3) result(ic)
     logical, dimension(:,:,:), intent(in)  :: ltab       ! Mask
     integer, dimension(:),     intent(out) :: i1, i2, i3 ! Positions of elements with 'true' value
     integer                                :: ic         ! Total number of 'true' values
   end function Countjv3d
end interface Countjv

#ifdef MNH_OPENACC
public :: Countjv_device
interface Countjv_device
   subroutine Countjv1d_device(ltab, i1,ic)
     logical, dimension(:), intent(in)  :: ltab ! Mask
     integer, dimension(:), intent(out) :: i1   ! Positions of elements with 'true' value
     integer,               intent(out) :: ic   ! Total number of 'true' values
   end subroutine Countjv1d_device
   subroutine Countjv2d_device(ltab, i1, i2, ic)
     logical, dimension(:,:), intent(in)  :: ltab   ! Mask
     integer, dimension(:),   intent(out) :: i1, i2 ! Positions of elements with 'true' value
     integer,                 intent(out) :: ic     ! Total number of 'true' values
   end subroutine Countjv2d_device
   subroutine Countjv3d_device(ltab, i1, i2, i3, ic)
     logical, dimension(:,:,:), intent(in)  :: ltab       ! Mask
     integer, dimension(:),     intent(out) :: i1, i2, i3 ! Positions of elements with 'true' value
     integer,                   intent(out) :: ic         ! Total number of 'true' values
   end subroutine Countjv3d_device
end interface
#endif

end module mode_tools

function Countjv1d(ltab,i1) result(ic)
  logical, dimension(:), intent(in)  :: ltab ! Mask
  integer, dimension(:), intent(out) :: i1   ! Positions of elements with 'true' value
  integer                            :: ic   ! Total number of 'true' values

  integer :: ji

  ic = 0

  do ji = 1, size( ltab, 1 )
    if ( ltab(ji ) ) then
      ic = ic +1
      i1(ic) = ji
    end if
  end do
end function Countjv1d

function Countjv2d(ltab,i1,i2) result(ic)
  logical, dimension(:,:), intent(in)  :: ltab   ! Mask
  integer, dimension(:),   intent(out) :: i1, i2 ! Positions of elements with 'true' value
  integer                              :: ic     ! Total number of 'true' values

  integer :: ji, jj

  ic = 0

  do jj = 1, size( ltab, 2 )
    do ji = 1, size( ltab, 1 )
      if ( ltab(ji, jj ) ) then
        ic = ic +1
        i1(ic) = ji
        i2(ic) = jj
      end if
    end do
  end do
end function Countjv2d

function Countjv3d(ltab,i1,i2,i3) result(ic)
  logical, dimension(:,:,:), intent(in)  :: ltab       ! Mask
  integer, dimension(:),     intent(out) :: i1, i2, i3 ! Positions of elements with 'true' value
  integer                                :: ic         ! Total number of 'true' values

  integer :: ji, jj, jk

  ic = 0

  do jk = 1, size( ltab, 3 )
    do jj = 1, size( ltab, 2 )
      do ji = 1, size( ltab, 1 )
        if ( ltab(ji, jj, jk ) ) then
          ic = ic +1
          i1(ic) = ji
          i2(ic) = jj
          i3(ic) = jk
        end if
      end do
    end do
  end do
end function Countjv3d

#ifdef MNH_OPENACC
subroutine Countjv1d_device(ltab, i1,ic)

  use mode_mppdb, only: mppdb_initialized
#ifndef _FAKEOPENACC
  use MODE_OPENACC_SET_DEVICE, only : mnh_idevice_type_current, acc_device_nvidia, acc_device_host
#endif

  logical, dimension(:), intent(in)  :: ltab ! Mask
  integer, dimension(:), intent(out) :: i1   ! Positions of elements with 'true' value
  integer,               intent(out) :: ic   ! Total number of 'true' values

  integer :: idx
  integer :: ji
  logical :: latomic


#ifndef _FAKEOPENACC
latomic = ( (.not. mppdb_initialized ) .and. (mnh_idevice_type_current .ne. acc_device_host ) )
#else
latomic = (.not. mppdb_initialized )
#endif
if (latomic) then
ic = 0
   
! acc kernels
   
!To allow comparisons... (i1 is not fully used)
!Can be removed in production
!   i1(:) = -999
  
!Warning: if "independent" is set, content of i1, i2 and i3 can vary between 2
! different runs of this subroutine BUT final result should be the same
!Comment the following line + atomic directives to have consistent values for debugging
!Warning: huge impact on performance

  do ji = 1, size( ltab, 1 )
    if ( ltab(ji ) ) then

      ic  = ic +1
      idx = ic

      i1(idx) = ji
    end if
  end do
! acc end parallel

else
   
ic = 0
   


!To allow comparisons... (i1 is not fully used)
!!$  i1(:) = -999

  do ji = 1, size( ltab, 1 )
    if ( ltab(ji ) ) then
      ic  = ic +1
      idx = ic
      i1(idx) = ji
    end if
  end do


end if



end subroutine Countjv1d_device

subroutine Countjv2d_device(ltab, i1, i2, ic)

  use mode_mppdb, only: mppdb_initialized
#ifndef _FAKEOPENACC
  use MODE_OPENACC_SET_DEVICE, only : mnh_idevice_type_current, acc_device_nvidia, acc_device_host
#endif

  logical, dimension(:,:), intent(in)  :: ltab   ! Mask
  integer, dimension(:),   intent(out) :: i1, i2 ! Positions of elements with 'true' value
  integer,                 intent(out) :: ic     ! Total number of 'true' values

  integer :: idx
  integer :: ji, jj
  logical :: LATOMIC


#ifndef _FAKEOPENACC
latomic = ( (.not. mppdb_initialized ) .and. (mnh_idevice_type_current .ne. acc_device_host ) )
#else
latomic = (.not. mppdb_initialized )
#endif
if (latomic) then
ic = 0   
     
! acc kernels

!To allow comparisons... (i1/i2 are not fully used)
!Can be removed in production
!   i1(:) = -999
!   i2(:) = -999

!Warning: if "independent" is set, content of i1, i2 and i3 can vary between 2
! different runs of this subroutine BUT final result should be the same
!Comment the following line + atomic directives to have consistent values for debugging
!Warning: huge impact on performance

  do jj = 1, size( ltab, 2 )
    do ji = 1, size( ltab, 1 )
      if ( ltab(ji, jj ) ) then

        ic  = ic +1
        idx = ic

        i1(idx) = ji
        i2(idx) = jj
      end if
    end do
  end do
! acc end parallel

else

ic = 0
   


!To allow comparisons... (i1/i2 are not fully used)
  i1(:) = -999
  i2(:) = -999

  do jj = 1, size( ltab, 2 )
    do ji = 1, size( ltab, 1 )
      if ( ltab(ji, jj ) ) then
        ic  = ic +1
        idx = ic
        i1(idx) = ji
        i2(idx) = jj
      end if
    end do
  end do


end if



end subroutine Countjv2d_device

subroutine Countjv3d_device(ltab, i1, i2, i3, ic)
  use mode_mppdb, only: mppdb_initialized
#ifndef _FAKEOPENACC
  use MODE_OPENACC_SET_DEVICE, only : mnh_idevice_type_current, acc_device_nvidia, acc_device_host
#endif

  logical, dimension(:,:,:), intent(in)  :: ltab       ! Mask
  integer, dimension(:),     intent(out) :: i1, i2, i3 ! Positions of elements with 'true' value
  integer,                   intent(out) :: ic         ! Total number of 'true' values

  integer :: idx
  integer :: ji, jj, jk
  logical :: latomic


#ifndef _FAKEOPENACC
latomic = ( (.not. mppdb_initialized ) .and. (mnh_idevice_type_current .ne. acc_device_host ) )
#else
latomic = (.not. mppdb_initialized )
#endif
if (latomic) then
ic = 0
   
! acc kernels

!To allow comparisons... (i1/i2/i3 are not fully used)
!Can be removed in production
!   i1(:) = -999
!   i2(:) = -999
!   i3(:) = -999

!Warning: if "independent" is set, content of i1, i2 and i3 can vary between 2
! different runs of this subroutine BUT final result should be the same
!Comment the following line + atomic directives to have consistent values for debugging
!Warning: huge impact on performance

  do jk = 1, size( ltab, 3 )
    do jj = 1, size( ltab, 2 )
       do ji = 1, size( ltab, 1 )
        if ( ltab(ji, jj, jk ) ) then

          ic  = ic +1
          idx = ic

          i1(idx) = ji
          i2(idx) = jj
          i3(idx) = jk
        end if
      end do
    end do
  end do
! acc end parellel

else

ic = 0  



!To allow comparisons... (i1/i2/i3 are not fully used)
!!$  i1(:) = -999
!!$  i2(:) = -999
!!$  i3(:) = -999

  do jk = 1, size( ltab, 3 )
    do jj = 1, size( ltab, 2 )
      do ji = 1, size( ltab, 1 )
        if ( ltab(ji, jj, jk ) ) then
          ic  = ic +1
          idx = ic
          i1(idx) = ji
          i2(idx) = jj
          i3(idx) = jk
        end if
      end do
    end do
  end do


end if



end subroutine Countjv3d_device
#endif

recursive subroutine Quicksort( ka, kbeg, kend, kpos )
  integer, dimension(:),           intent(inout) :: ka
  integer,                         intent(in)    :: kbeg, kend
  integer, dimension(:), optional, intent(inout) :: kpos

  integer :: ji, jj
  integer :: itmp, itmp2, ival

  ival = ka( ( kbeg + kend ) / 2 )
  ji = kbeg
  jj = kend
  do
     do while ( ka(ji) < ival )
        ji = ji + 1
     end do
     do while ( ival < ka(jj) )
        jj = jj - 1
     end do
     if ( ji >= jj ) exit

     itmp = ka(ji)
     ka(ji) = ka(jj)
     ka(jj) = itmp

     if ( present( kpos ) ) then
      itmp2 = kpos(ji)
      kpos(ji) = kpos(jj)
      kpos(jj) = itmp2
     end if

     ji=ji+1
     jj=jj-1
  end do
  if ( kbeg   < ji - 1 ) call Quicksort( ka, kbeg,   ji - 1, kpos )
  if ( jj + 1 < kend   ) call Quicksort( ka, jj + 1, kend,   kpos )
end subroutine Quicksort

function Upcase(hstring)
  character(len=*), intent(in) :: hstring
  character(len=len(hstring))  :: upcase

  integer :: jc
  integer, parameter :: iamin = iachar("a")
  integer, parameter :: iamaj = iachar("A")

  do jc = 1,len(hstring)
    if ( hstring(jc:jc) >= "a" .and. hstring(jc:jc) <= "z" ) then
      upcase(jc:jc) = achar( iachar( hstring(jc:jc) ) - iamin + iamaj )
    else
      upcase(jc:jc) = hstring(jc:jc)
    end if
  end do
end function Upcase
