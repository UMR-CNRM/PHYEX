!     ######spl
     MODULE MODI_TM06  
!    ################ 
!
INTERFACE
!
      SUBROUTINE TM06(KKA,KKU,KKL,PTHVREF,PBL_DEPTH,PZZ,PSFTH,PMWTH,PMTH2)
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=ARO
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHVREF    ! reference potential temperature
REAL, DIMENSION(:,:),   INTENT(IN) :: PBL_DEPTH ! boundary layer height
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ        ! altitude of flux levels
REAL, DIMENSION(:,:),   INTENT(IN) :: PSFTH      ! surface heat flux
REAL, DIMENSION(:,:,:), INTENT(OUT):: PMWTH      ! w'2th'
REAL, DIMENSION(:,:,:), INTENT(OUT):: PMTH2      ! w'th'2
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TM06
!
END INTERFACE
!
END MODULE MODI_TM06
