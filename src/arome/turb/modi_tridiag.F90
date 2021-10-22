!     ######spl
      MODULE MODI_TRIDIAG
!     ###################
INTERFACE
!
       SUBROUTINE TRIDIAG(KKA,KKU,KKL,PVARM,PA,PTSTEP,PEXPL,PIMPL, &
                                  PRHODJ,PSOURCE,PVARP )
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=AR
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PVARM       ! variable at t-1  
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PA          ! upper diag. elements
REAL,                      INTENT(IN)  :: PTSTEP      ! Double time step
REAL,                      INTENT(IN)  :: PEXPL,PIMPL ! weights of the temporal scheme
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PRHODJ      ! (dry rho)*J
REAL, DIMENSION(:,:,:),    INTENT(IN)  :: PSOURCE     ! source term of PVAR    
!
REAL, DIMENSION(:,:,:),    INTENT(OUT) :: PVARP       ! variable at t+1        
!
END SUBROUTINE TRIDIAG
!
END INTERFACE
!
END MODULE MODI_TRIDIAG 
