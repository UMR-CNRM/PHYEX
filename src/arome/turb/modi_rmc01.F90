!     ######spl
      MODULE MODI_RMC01
!     ################
INTERFACE
      SUBROUTINE RMC01(HTURBLEN,KKA,KKU,KKL,PZZ,PDXX,PDYY,PDZZ,PDIRCOSZW, &
                       PSBL_DEPTH, PLMO, PLK, PLEPS                )
!
CHARACTER(LEN=4),         INTENT(IN)    :: HTURBLEN ! type of mixing length
INTEGER,                  INTENT(IN)    :: KKA           !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU           !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL           !vert. levels type 1=MNH -1=ARO
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ   ! altitude of flux points
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDXX  ! width of grid mesh (X dir)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDYY  ! width of grid mesh (Y dir)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ  ! width of vert. layers
REAL, DIMENSION(:,:),     INTENT(IN)    :: PDIRCOSZW ! Director Cosinus 
REAL, DIMENSION(:,:),     INTENT(IN)    :: PSBL_DEPTH! SBL depth
REAL, DIMENSION(:,:),     INTENT(IN)    :: PLMO  ! Monin Obuhkov length
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PLK   ! Mixing length
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PLEPS ! Dissipative length

END SUBROUTINE RMC01
END INTERFACE
END MODULE MODI_RMC01
