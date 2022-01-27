!     ######spl
     MODULE MODI_TM06_H  
!    ################ 
!
INTERFACE
!
      SUBROUTINE TM06_H(KKB,KKTB,KKTE,PTSTEP,PZZ,PFLXZ,PBL_DEPTH)
!
INTEGER,                INTENT(IN)    :: KKB       ! index of 1st physical level
                                                   ! close to ground 
INTEGER,                INTENT(IN)    :: KKTB      ! first physical level in k
INTEGER,                INTENT(IN)    :: KKTE      ! last physical level in k
REAL,                   INTENT(IN)    :: PTSTEP    ! Double time step
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PZZ       ! altitude of flux levels
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFLXZ     ! heat flux
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PBL_DEPTH ! boundary layer height
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TM06_H
!
END INTERFACE
!
END MODULE MODI_TM06_H
