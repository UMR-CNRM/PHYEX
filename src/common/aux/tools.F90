!MNH_LIC Copyright 2019-2025 CNRS, Meteo-France and Universite Paul Sabatier
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
  implicit none
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
  implicit none
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
  implicit none
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
