!     ######spl
     MODULE MODI_UPDRAFT_SOPE
!    #################################
!
INTERFACE
!
      SUBROUTINE UPDRAFT_SOPE(KRR,KRRL,KRRI,OMIXUV,                & 
                         PZZ,PDZZ,PSFTH,PSFRV,PPABSM,PRHODREF,     &
                         PTKEM,PTHM,PRM,PTHLM,PRTM,PUM,PVM,PSVM,   &
                         PTHL_UP,PRT_UP,PRV_UP,PU_UP,PV_UP,PSV_UP, &    
                         PRC_UP,PRI_UP,PTHV_UP,PW_UP,PFRAC_UP,PEMF,&
                         PDETR,PENTR,KKLCL,KKETL,KKCTL )
!
!
!
INTEGER,                INTENT(IN)   :: KRR       ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL      ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI      ! number of ice water var.
LOGICAL,                INTENT(IN)   :: OMIXUV    ! True if mixing of momentum

REAL, DIMENSION(:,:,:), INTENT(IN)   :: PZZ       !  Height at the flux point

REAL, DIMENSION(:,:,:), INTENT(IN)   :: PDZZ      ! depth between mass levels  

REAL, DIMENSION(:,:),   INTENT(IN)   ::  PSFTH,PSFRV
                                                  ! normal surface fluxes of theta,rv
!
!    prognostic variables at t- deltat
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PPABSM      ! Pressure at time t-1
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PRHODREF    ! dry density of the
                                                     ! reference state
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PTKEM       ! TKE
REAL, DIMENSION(:,:,:),   INTENT(IN) ::  PUM,PVM     ! momentum
!
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM        ! Scalar variables 
!
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PTHM       ! pot. temp. = PTHLM in turb.f90
REAL, DIMENSION(:,:,:,:), INTENT(IN)   ::  PRM        ! water species 
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PTHLM,PRTM !cons. var.       
REAL, DIMENSION(:,:,:),   INTENT(OUT)  ::  PTHL_UP,PRT_UP        ! updraft properties
REAL, DIMENSION(:,:,:),   INTENT(OUT)  ::  PRV_UP,PRC_UP,PRI_UP,&!Thl,Rt,Rv,Rc,Ri
                                           PW_UP,PFRAC_UP,PEMF, &!w,Updraft Fraction, Mass Flux
                                           PDETR,PENTR,PTHV_UP, &!entrainment, detrainment, ThV
                                           PU_UP, PV_UP          !updraft wind component
REAL, DIMENSION(:,:,:,:), INTENT(OUT)  ::  PSV_UP                ! updraft scalar variables                                             
INTEGER, DIMENSION(:,:),  INTENT(OUT)  ::  KKLCL,KKETL,KKCTL     !index for LCL,ETL,CTL                                      
!
!
END SUBROUTINE UPDRAFT_SOPE

END INTERFACE
!
END MODULE MODI_UPDRAFT_SOPE
