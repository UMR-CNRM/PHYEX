!     ######spl
      MODULE MODI_MPDATA_SCALAR
!     #########################
INTERFACE
      SUBROUTINE MPDATA_SCALAR   ( KLITER, HLBCX, HLBCY, KSV,             &
                                   PTSTEP, PRHODJ, PSVM, PSVT,            &
                                   PRUCT, PRVCT, PRWCT, PRSVS             )                
!
INTEGER,                  INTENT(IN)    :: KLITER  ! Number of iterations MPDATA 
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KSV  ! Number of Scalar Variables
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVM
                                                ! Variables at t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT,PRVCT,PRWCT ! Contravariants
                                                ! components of the momentum
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT
                                                ! Variables at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS
                                                ! Sources terms
END SUBROUTINE MPDATA_SCALAR 
!
END INTERFACE
!
END MODULE MODI_MPDATA_SCALAR 
