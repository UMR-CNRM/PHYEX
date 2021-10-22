!     ######spl
      MODULE MODI_RZCOLX
!     ##################
!
INTERFACE
!
      SUBROUTINE RZCOLX( KND, PALPHAX, PNUX, PALPHAZ, PNUZ,                  &
                         PEXZ, PEXMASSZ, PFALLX, PEXFALLX, PFALLZ, PEXFALLZ, &
                         PLBDAXMAX, PLBDAZMAX, PLBDAXMIN, PLBDAZMIN,         &
                         PDINFTY, PRZCOLX                                    )
!
INTEGER, INTENT(IN) :: KND    ! Number of discrete size intervals in DX and DZ
!
!
REAL, INTENT(IN) :: PALPHAX   ! First shape parameter of the specy X
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUX      ! Second shape parameter of the specy X
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PALPHAZ   ! First shape parameter of the specy Z
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PNUZ      ! Second shape parameter of the specy Z
                              ! size distribution (generalized gamma law)
REAL, INTENT(IN) :: PEXZ      ! Efficiency of specy X collecting specy Z
REAL, INTENT(IN) :: PEXMASSZ  ! Mass exponent of specy Z
REAL, INTENT(IN) :: PFALLX    ! Fall speed constant of specy X
REAL, INTENT(IN) :: PEXFALLX  ! Fall speed exponent of specy X
REAL, INTENT(IN) :: PFALLZ    ! Fall speed constant of specy Z
REAL, INTENT(IN) :: PEXFALLZ  ! Fall speed exponent of specy Z
REAL, INTENT(IN) :: PLBDAXMAX ! Maximun slope of size distribution of specy X
REAL, INTENT(IN) :: PLBDAZMAX ! Maximun slope of size distribution of specy Z
REAL, INTENT(IN) :: PLBDAXMIN ! Minimun slope of size distribution of specy X
REAL, INTENT(IN) :: PLBDAZMIN ! Minimun slope of size distribution of specy Z
REAL, INTENT(IN) :: PDINFTY   ! Factor to define the largest diameter up to
                              ! which the diameter integration is performed
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRZCOLX ! Scaled fall speed difference in
                                               ! the mass collection kernel as a
                                               ! function of LAMBDAX and LAMBDAZ
!
      END SUBROUTINE RZCOLX
!
END INTERFACE
!
      END MODULE MODI_RZCOLX
