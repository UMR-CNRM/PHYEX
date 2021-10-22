!     ######spl
      MODULE MODI_RRCOLSS
!     ###################
!
INTERFACE
!
      SUBROUTINE RRCOLSS( KND, PALPHAS, PNUS, PALPHAR, PNUR,                 &
                         PESR, PEXMASSR, PFALLS, PEXFALLS, PFALLR, PEXFALLR, &
                         PLBDASMAX, PLBDARMAX, PLBDASMIN, PLBDARMIN,         &
                         PDINFTY, PRRCOLSS, PAG, PBS, PAS                    )
!
INTEGER, INTENT(IN) :: KND    ! Number of discrete size intervals in DS and DR  
!
REAL, INTENT(IN) :: PALPHAS   ! First shape parameter of the aggregates 
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUS      ! Second shape parameter of the aggregates
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PALPHAR   ! First shape parameter of the rain  
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUR      ! Second shape parameter of the rain 
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PESR      ! Efficiency of aggregates collecting rain 
REAL, INTENT(IN) :: PEXMASSR  ! Mass exponent of rain 
REAL, INTENT(IN) :: PFALLS    ! Fall speed constant of aggregates
REAL, INTENT(IN) :: PEXFALLS  ! Fall speed exponent of aggregates
REAL, INTENT(IN) :: PFALLR    ! Fall speed constant of rain 
REAL, INTENT(IN) :: PEXFALLR  ! Fall speed exponent of rain 
REAL, INTENT(IN) :: PLBDASMAX ! Maximun slope of size distribution of aggregates
REAL, INTENT(IN) :: PLBDARMAX ! Maximun slope of size distribution of rain 
REAL, INTENT(IN) :: PLBDASMIN ! Minimun slope of size distribution of aggregates
REAL, INTENT(IN) :: PLBDARMIN ! Minimun slope of size distribution of rain 
REAL, INTENT(IN) :: PDINFTY   ! Factor to define the largest diameter up to
                              ! which the diameter integration is performed
REAL, INTENT(IN) :: PAG, PBS, PAS
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRRCOLSS! Scaled fall speed difference in
                                               ! the mass collection kernel as a
                                               ! function of LAMBDAX and LAMBDAZ
!
      END SUBROUTINE RRCOLSS
!
END INTERFACE
!
      END MODULE MODI_RRCOLSS
