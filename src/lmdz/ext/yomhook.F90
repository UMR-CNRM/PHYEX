! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! DR_HOOK is a profiling and debugging system for the IFS, and should
! be called at the beginning and end of each subroutine. This is a
! dummy implementation for offline packages.

module yomhook

  use parkind1, only : jprb

  public

  save

  logical :: lhook = .false.

  integer, parameter :: jphook=jprb

contains

  subroutine dr_hook(proc_name, iswitch, proc_key)

    use parkind1, only : jprb

    character(len=*), intent(in)    :: proc_name
    integer,          intent(in)    :: iswitch
    real(jprb),       intent(inout) :: proc_key
    ! Do nothing!

  end subroutine dr_hook

end module yomhook
