!     ######spl
     MODULE MODI_eqsam_v03d_sub
!!   ########################
!!
INTERFACE
!!
SUBROUTINE eqsam_v03d_sub(yi,yo,nca,nco,iopt,loop,imax)
IMPLICIT NONE
integer,intent(in)                            :: nca,nco,imax,loop
integer,intent(inout)                         :: iopt
real,dimension(:,:),intent(in)           :: yi
real,dimension(:,:),intent(inout)        :: yo
END SUBROUTINE eqsam_v03d_sub
!!
END INTERFACE
!!
END MODULE MODI_eqsam_v03d_sub
