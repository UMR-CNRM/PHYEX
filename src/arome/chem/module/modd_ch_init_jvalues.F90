!     ######spl
      MODULE MODD_CH_INIT_JVALUES
!!    ###########################
!!
!!*** *MODD_CH_INIT_JVALUES*
!!
!!    PURPOSE 
!!    -------
!!    Store J values variables common to all models 
!!     (i.e. before spatial and temporal interpolation)
!!    XJDATA is calculated at the first call of model 1. 
!!    XJDATA is calculated for a discrete number 
!!     of solar zenith angle, altitude and albedo.
!
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
SAVE
!
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: XJDATA 
INTEGER                               :: NSZA_INCR = 99 + 1
REAL, ALLOCATABLE, DIMENSION(:)       :: XSZA_JVAL
INTEGER, PARAMETER                    :: NZZ_JVAL = 30 + 1
REAL, ALLOCATABLE, DIMENSION(:)       :: XZZ_JVAL
INTEGER, PARAMETER                    :: JPJVMAX = 21    
INTEGER                               :: NBALB = 10   
!
END MODULE MODD_CH_INIT_JVALUES
