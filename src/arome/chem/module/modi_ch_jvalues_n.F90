!     ######spl
       MODULE MODI_CH_JVALUES_n
!      #############################################
!
       INTERFACE
       SUBROUTINE CH_JVALUES_n(KVECNPT, KIINDEX, KJINDEX, KMODELLEVEL, KRATES, PJVALUES, PRATES)
       IMPLICIT NONE
       INTEGER, DIMENSION(:), INTENT(IN)  :: KMODELLEVEL
       INTEGER,               INTENT(IN)  :: KVECNPT 
       INTEGER, DIMENSION(:), INTENT(IN)  :: KIINDEX,KJINDEX 
                                                     ! current model grid point
       INTEGER,               INTENT(IN)  :: KRATES          
                                                     ! dimension of PRATES
       REAL, DIMENSION(KVECNPT,KRATES), INTENT(OUT) :: PRATES          
                                                     ! photolysis rates 
       REAL,DIMENSION(:,:,:,:), INTENT(IN) :: PJVALUES    ! Tuv coefficients

       END SUBROUTINE CH_JVALUES_n
       END INTERFACE
       END MODULE MODI_CH_JVALUES_n
