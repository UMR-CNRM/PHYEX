!     ######spl
      MODULE MODD_CH_METEO
!!    ######################
!!
!!*** *MODD_CH_METEO*
!!
!!    PURPOSE
!!    -------
!       contains the meteovariables that will be updated by ch_update_meteo
!!
!!**  AUTHOR
!!    ------
!!    K. Suhre     *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 20/04/99
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
!
IMPLICIT NONE
SAVE
!
INTEGER                           :: NMETEORECS
                    ! number of records
!
INTEGER                           :: NMETEORECACT = 1
                    ! actual record (used in temporal interpolation)
!
REAL, DIMENSION(:),   ALLOCATABLE :: XMETEOTIME
                    ! the time of the individual records
!
REAL, DIMENSION(:,:), ALLOCATABLE :: XMETEODATA
                    ! the meteodata (the first dimension is the number of
                    ! elements NMETEOVARS, the second dimension is the number
                    ! of records)
!
END MODULE MODD_CH_METEO
