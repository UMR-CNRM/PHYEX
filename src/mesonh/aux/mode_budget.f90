!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications
!  P. Wautelet 28/01/2020: new subroutines: Budget_store_init, Budget_store_end and Budget_source_id_find in new module mode_budget
!  P. Wautelet 17/08/2020: treat LES budgets correctly
!  P. Wautelet 05/03/2021: measure cpu_time for budgets
!-----------------------------------------------------------------

!#################
module mode_budget
!#################

use modd_budget,     only: cbutype, nbutime, tbudgetdata, xtime_bu, xtime_bu_process
use modd_les_budget, only: xtime_les_bu, xtime_les_bu_process

use modi_cart_compress, only: Cart_compress
use modi_mask_compress, only: Mask_compress
use modi_second_mnh,    only: Second_mnh

use mode_msg

implicit none

private

public :: Budget_store_init
public :: Budget_store_end
public :: Budget_store_add

real :: ztime1, ztime2

contains

subroutine Budget_store_init( tpbudget, hsource, pvars )
  use modd_les, only: lles_call

  type(tbudgetdata),      intent(inout) :: tpbudget ! Budget datastructure
  character(len=*),       intent(in)    :: hsource  ! Name of the source term
  real, dimension(:,:,:), intent(in)    :: pvars    ! Current value to be stored

  integer :: iid ! Reference number of the current source term

  call Print_msg( NVERB_DEBUG, 'BUD', 'Budget_store_init', trim( tpbudget%cname )//':'//trim( hsource ) )

  if ( lles_call ) then
    call Second_mnh( ztime1 )

    if ( allocated( tpbudget%xtmplesstore ) ) then
      call Print_msg( NVERB_ERROR, 'BUD', 'Budget_store_init', 'xtmplesstore already allocated' )
    else
      allocate( tpbudget%xtmplesstore( Size( pvars, 1 ), Size( pvars, 2 ), Size ( pvars, 3 )  ) )
    end if
    tpbudget%xtmplesstore(:, :, :) = pvars(:, :, :)

    tpbudget%clessource = hsource

    call Second_mnh( ztime2 )
    xtime_les_bu         = xtime_les_bu         + ztime2 - ztime1
    xtime_les_bu_process = xtime_les_bu_process + ztime2 - ztime1
  end if

  ! Nothing else to do if budgets are not enabled
  if ( .not. tpbudget%lenabled ) return

  call Second_mnh( ztime1 )

  call Budget_source_id_find( tpbudget, hsource, iid )

  if ( tpbudget%ntmpstoresource /= 0 ) then
    cmnhmsg(1) = 'ntmpstoresource already set (previous call to '//'Budget_store_end missing?)'
    cmnhmsg(2) = 'Set for:    ' // Trim( tpbudget%cname ) // ':' // Trim( tpbudget%tsources(tpbudget%ntmpstoresource)%cmnhname )
    cmnhmsg(3) = 'Working on: ' // Trim( tpbudget%cname ) // ':' // Trim( hsource )
    call Print_msg( NVERB_ERROR, 'BUD', 'Budget_store_init' )
  end if

  if ( tpbudget%tsources(iid)%ldonotinit ) then
    ! If ldonotinit is set, this subroutine should not be called
    call Print_msg( NVERB_ERROR, 'BUD', 'Budget_store_init', 'should not be called for ' &
                    //trim( tpbudget%cname )//':'//trim( hsource ) )
    return
  end if

  if ( tpbudget%tsources(iid)%lenabled ) then
    if ( tpbudget%ntmpstoresource /= 0 ) then
      call Print_msg( NVERB_ERROR, 'BUD', 'Budget_store_init', 'xtmpstore already used by ' &
                      //trim( tpbudget%tsources(tpbudget%ntmpstoresource)%cmnhname ) )
      return
    end if

    tpbudget%ntmpstoresource = iid

    !Store data into the budget temporary array
    !This value will be subtracted from the next one (in Budget_store_end) to get the evolution of the array between the 2 calls
    if ( cbutype == 'CART' ) then
      tpbudget%xtmpstore(:, :, :) = Cart_compress( pvars(:, :, :) )
    else if ( cbutype == 'MASK' ) then
      tpbudget%xtmpstore(:, nbutime, :) = Mask_compress( pvars(:, :, :) )
    else
      call Print_msg( NVERB_ERROR, 'BUD', 'Budget_store_init', 'unknown cbutype: '//trim( cbutype ) )
    end if
  end if

  call Second_mnh( ztime2 )
  xtime_bu         = xtime_bu         + ztime2 - ztime1
  xtime_bu_process = xtime_bu_process + ztime2 - ztime1

  end subroutine Budget_store_init


subroutine Budget_store_end( tpbudget, hsource, pvars )
  use modd_les, only: lles_call

  use modi_les_budget, only: Les_budget

  type(tbudgetdata),      intent(inout) :: tpbudget ! Budget datastructure
  character(len=*),       intent(in) :: hsource     ! Name of the source term
  real, dimension(:,:,:), intent(in) :: pvars       ! Current value to be stored

  integer :: iid    ! Reference number of the current source term
  integer :: igroup ! Number of the group where to store the source term
  real, dimension(:,:,:), allocatable :: zvars_add

  call Print_msg( NVERB_DEBUG, 'BUD', 'Budget_store_end', trim( tpbudget%cname )//':'//trim( hsource ) )

  if ( lles_call ) then
    if ( hsource /= tpbudget%clessource ) &
      call Print_msg( NVERB_FATAL, 'BUD', 'Budget_store_end', 'hsource not the same as in Budget_store_init (' &
                      // Trim( hsource ) // ' / ' // Trim( tpbudget%clessource ) // ')' )

    tpbudget%clessource = 'reset'

    if ( allocated( tpbudget%xtmplesstore ) ) then
      ! Do the call to Les_budget with oadd=.true.
      ! This is necessary when the call to Budget_store_init was done with pvars not strictly
      ! equal to the source term
      Allocate( zvars_add( Size( pvars, 1 ), Size( pvars, 2 ), Size ( pvars, 3 ) ) )
      zvars_add(:, :, :) = pvars(:, :, :) - tpbudget%xtmplesstore(:, :, :)
      call Les_budget( zvars_add, tpbudget%nid, hsource, oadd = .true. )
      Deallocate( zvars_add )
      Deallocate( tpbudget%xtmplesstore )
    else
      call Les_budget( pvars, tpbudget%nid, hsource, oadd = .false. )
    end if
  end if

  ! Nothing to do if budgets are not enabled
  if ( .not. tpbudget%lenabled ) return

  call Second_mnh( ztime1 )

  call Budget_source_id_find( tpbudget, hsource, iid )

  if ( tpbudget%tsources(iid)%lenabled ) then
    if ( iid /= tpbudget%ntmpstoresource .and. .not.tpbudget%tsources(iid)%ldonotinit ) then
      if ( tpbudget%ntmpstoresource == 0 ) then
        call Print_msg( NVERB_ERROR, 'BUD', 'Budget_store_end', 'ntmpstoresource not set for ' &
                        //trim( tpbudget%tsources(iid)%cmnhname ) )
      else
        call Print_msg( NVERB_ERROR, 'BUD', 'Budget_store_end', 'xtmpstore used by an other source: '    &
                        //trim( tpbudget%tsources(tpbudget%ntmpstoresource)%cmnhname )//', expected: '   &
                        //trim( tpbudget%tsources(iid)%cmnhname ) )
      end if
    end if

    !Store data into the budget array
    !The values are computed by the difference between the values stored in the temporary array (filled in Budget_store_init)
    !and the current values added to the already stored ones.
    !Except if ldonotinit is true. In that case, overwrite the array.
    igroup = tpbudget%tsources(iid)%ngroup
    if ( cbutype == 'CART' ) then
      if ( tpbudget%tsources(iid)%ldonotinit ) then
        if ( tpbudget%tsources(iid)%loverwrite ) then
          tpbudget%tgroups(igroup)%xdata(:, :, :) =   Cart_compress( pvars(:, :, :) )
        else
          tpbudget%tgroups(igroup)%xdata(:, :, :) =   tpbudget%tgroups(igroup)%xdata(:, :, :) &
                                                    + Cart_compress( pvars(:, :, :) )
        end if
      else
        if ( tpbudget%tsources(iid)%loverwrite ) then
          tpbudget%tgroups(igroup)%xdata(:, :, :) =   Cart_compress( pvars(:, :, :) )          &
                                                    - tpbudget%xtmpstore(:, :, :)
        else
          tpbudget%tgroups(igroup)%xdata(:, :, :) =   tpbudget%tgroups(igroup)%xdata(:, :, :) &
                                                     + Cart_compress( pvars(:, :, :) )          &
                                                     - tpbudget%xtmpstore(:, :, :)
        end if
      end if
    else if ( cbutype == 'MASK' ) then
      if ( tpbudget%tsources(iid)%ldonotinit ) then
        if ( tpbudget%tsources(iid)%loverwrite ) then
          tpbudget%tgroups(igroup)%xdata(:, nbutime, :) =   Mask_compress( pvars(:, :, :) )
        else
          tpbudget%tgroups(igroup)%xdata(:, nbutime, :) =   tpbudget%tgroups(igroup)%xdata(:, nbutime, :) &
                                                           + Mask_compress( pvars(:, :, :) )
        end if
      else
        if ( tpbudget%tsources(iid)%loverwrite ) then
          tpbudget%tgroups(igroup)%xdata(:, nbutime, :) =   Mask_compress( pvars(:, :, :) )   &
                                                          - tpbudget%xtmpstore(:, nbutime, :)
        else
          tpbudget%tgroups(igroup)%xdata(:, nbutime, :) =   tpbudget%tgroups(igroup)%xdata(:, nbutime, :) &
                                                          + Mask_compress( pvars(:, :, :) )                &
                                                          - tpbudget%xtmpstore(:, nbutime, :)
        end if
      end if
    else
      call Print_msg( NVERB_ERROR, 'BUD', 'Budget_store_end', 'unknown cbutype: '//trim( cbutype ) )
    end if

    ! Release the budget temporary array
    tpbudget%ntmpstoresource = 0
  end if

  call Second_mnh( ztime2 )
  xtime_bu         = xtime_bu         + ztime2 - ztime1
  xtime_bu_process = xtime_bu_process + ztime2 - ztime1

end subroutine Budget_store_end


subroutine Budget_store_add( tpbudget, hsource, pvars )
  use modd_les, only: lles_call

  use modi_les_budget, only: Les_budget

  type(tbudgetdata),      intent(inout) :: tpbudget ! Budget datastructure
  character(len=*),       intent(in) :: hsource     ! Name of the source term
  real, dimension(:,:,:), intent(in) :: pvars       ! Current value to be stored

  integer :: iid    ! Reference number of the current source term
  integer :: igroup ! Number of the group where to store the source term

  call Print_msg( NVERB_DEBUG, 'BUD', 'Budget_store_add', trim( tpbudget%cname )//':'//trim( hsource ) )

  if ( tpbudget%ntmpstoresource /= 0 ) &
    call Print_msg( NVERB_ERROR, 'BUD', 'Budget_store_add', 'inside a Budget_store_init/Budget_store_end zone' )

  if ( lles_call ) call Les_budget( pvars, tpbudget%nid, hsource, oadd = .true. )

  ! Nothing to do if budgets are not enabled
  if ( .not. tpbudget%lenabled ) return

  call Second_mnh( ztime1 )

  call Budget_source_id_find( tpbudget, hsource, iid )

  if ( tpbudget%tsources(iid)%lenabled ) then
    if ( tpbudget%tsources(iid)%loverwrite ) &
      call Print_msg( NVERB_ERROR, 'BUD', 'Budget_store_add', 'loverwrite=.true. is not allowed' )

    !Store data into the budget array
    igroup = tpbudget%tsources(iid)%ngroup
    if ( cbutype == 'CART' ) then
      tpbudget%tgroups(igroup)%xdata(:, :, :) =   tpbudget%tgroups(igroup)%xdata(:, :, :) &
                                                + Cart_compress( pvars(:, :, :) )
    else if ( cbutype == 'MASK' ) then
      tpbudget%tgroups(igroup)%xdata(:, nbutime, :) =   tpbudget%tgroups(igroup)%xdata(:, nbutime, :) &
                                                      + Mask_compress( pvars(:, :, :) )
    else
      call Print_msg( NVERB_ERROR, 'BUD', 'Budget_store_add', 'unknown cbutype: '//trim( cbutype ) )
    end if
  end if

  call Second_mnh( ztime2 )
  xtime_bu         = xtime_bu         + ztime2 - ztime1
  xtime_bu_process = xtime_bu_process + ztime2 - ztime1

end subroutine Budget_store_add


subroutine Budget_source_id_find( tpbudget, hsource, kid )
  type(tbudgetdata), intent(in)  :: tpbudget ! Budget datastructure
  character(len=*),  intent(in)  :: hsource  ! Name of the source term
  integer,           intent(out) :: kid      ! Reference number of the current source term

  integer :: iid
  integer :: ji

  call Print_msg( NVERB_DEBUG, 'BUD', 'Budget_source_id_find', trim( tpbudget%cname )//':'//trim( hsource ) )

  iid = 0
  do ji = 1, tpbudget%nsources
    if ( trim( hsource ) == trim( tpbudget%tsources(ji)%cmnhname ) ) then
      iid = ji
      exit
    end if
  end do

  if ( iid > 0 ) then
    call Print_msg( NVERB_DEBUG, 'BUD', 'Budget_source_id_find', trim( tpbudget%cname )//':'//trim( hsource )//' found' )
  else
    !Search also in the non-available source term list
    do ji = tpbudget%nsources + 1, tpbudget%nsourcesmax
      if ( trim( hsource ) == trim( tpbudget%tsources(ji)%cmnhname ) ) then
        iid = ji
        exit
      end if
    end do

    if ( iid == 0 ) then
      call Print_msg( NVERB_ERROR, 'BUD', 'Budget_source_id_find', trim( tpbudget%cname )//':'//trim( hsource )//' not found' )
    else
      cmnhmsg(1) = Trim( tpbudget%cname ) // ':' // Trim( hsource ) // ' found'
      cmnhmsg(2) = 'in non-available source term list.'
      cmnhmsg(3) = 'Check availability condition in Ini_budget.'
      call Print_msg( NVERB_ERROR, 'BUD', 'Budget_source_id_find' )
    end if
  end if

  kid = iid
end subroutine Budget_source_id_find

end module mode_budget
