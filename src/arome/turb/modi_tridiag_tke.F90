!     ######spl
      MODULE MODI_TRIDIAG_TKE
!     ##########################
INTERFACE
!
       SUBROUTINE TRIDIAG_TKE(KKA,KKU,KKL,PVARM,PA,PTSTEP,PEXPL,PIMPL, &
                                  PRHODJ,PSOURCE,PDIAG,PVARP )
!
INTEGER,                INTENT(IN)   :: KKA     !near ground array index  
INTEGER,                INTENT(IN)   :: KKU     !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL     !vert. levels type 1=MNH -1=ARO
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PVARM       ! variable at t-1  
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PA          ! upper diag. elements
REAL,                   INTENT(IN)  :: PTSTEP      ! Double time step
REAL,                   INTENT(IN)  :: PEXPL,PIMPL ! weights of the temporal scheme
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHODJ      ! (dry rho)*J
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSOURCE     ! source term of PVAR    
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDIAG       ! diagonal term linked to
                                                      ! the implicit dissipation
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PVARP       ! variable at t+1        
!
END SUBROUTINE TRIDIAG_TKE
!
END INTERFACE
!
END MODULE MODI_TRIDIAG_TKE
