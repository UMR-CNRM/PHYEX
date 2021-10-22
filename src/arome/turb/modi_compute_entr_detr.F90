!     ######spl
     MODULE MODI_COMPUTE_ENTR_DETR
!    ##############################
!
INTERFACE
!
      SUBROUTINE COMPUTE_ENTR_DETR(KK,KKB,KKE,KKL,OTEST,OTESTLCL,HFRAC_ICE, &
                          PFRAC_ICE,PRHODREF,PPRE_MINUS_HALF,&
                          PPRE_PLUS_HALF,PZZ,PDZZ,&
                          PTHVM,PTHLM,PRTM,PW_UP2,PTH_UP,&
                          PTHL_UP,PRT_UP,PLUP,&
                          PRC_UP,PRI_UP,PTHV_UP,&
                          PRSAT_UP,PRC_MIX,PRI_MIX,      &
                          PENTR,PDETR,PENTR_CLD,PDETR_CLD,&
                          PBUO_INTEG_DRY,PBUO_INTEG_CLD,&
                          PPART_DRY)

!
!
!
INTEGER,                INTENT(IN)   :: KK          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
LOGICAL,DIMENSION(:),   INTENT(IN)   :: OTEST
LOGICAL,DIMENSION(:),   INTENT(IN)   :: OTESTLCL !test of condensation 
CHARACTER*1,            INTENT(IN)   :: HFRAC_ICE 
REAL, DIMENSION(:)  ,INTENT(IN)     :: PFRAC_ICE
 
!
!    prognostic variables at t- deltat
REAL, DIMENSION(:),   INTENT(IN) ::  PRHODREF  !rhodref
REAL, DIMENSION(:),   INTENT(IN) ::  PPRE_MINUS_HALF ! Pressure at flux level KK
REAL, DIMENSION(:),   INTENT(IN) ::  PPRE_PLUS_HALF ! Pressure at flux level KK+KKL
REAL, DIMENSION(:,:), INTENT(IN) ::  PZZ       !  Height at the flux point
REAL, DIMENSION(:,:), INTENT(IN) ::  PDZZ      !   metrics coefficient
REAL, DIMENSION(:,:), INTENT(IN) ::  PTHVM     ! ThetaV environment 
!
!
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(:,:), INTENT(IN)     ::  PTHLM        ! Thetal
REAL, DIMENSION(:,:), INTENT(IN)     ::  PRTM         ! total mixing ratio
REAL, DIMENSION(:,:), INTENT(IN)     ::  PW_UP2       ! Vertical velocity^2
REAL, DIMENSION(:),   INTENT(IN)     ::  PTH_UP,PTHL_UP,PRT_UP  ! updraft properties
REAL, DIMENSION(:),   INTENT(IN)     ::  PLUP         ! LUP compute from the ground
REAL, DIMENSION(:),   INTENT(IN)     ::  PRC_UP,PRI_UP   ! Updraft cloud content
REAL, DIMENSION(:),   INTENT(IN)     ::  PTHV_UP ! Thetav of updraft
REAL, DIMENSION(:),   INTENT(IN)     ::  PRSAT_UP ! Mixing ratio at saturation in updraft
REAL, DIMENSION(:),   INTENT(INOUT)  ::  PRC_MIX, PRI_MIX      ! Mixture cloud content
REAL, DIMENSION(:),   INTENT(OUT)    ::  PENTR        ! Mass flux entrainment of the updraft
REAL, DIMENSION(:),   INTENT(OUT)    ::  PDETR        ! Mass flux detrainment of the updraft
REAL, DIMENSION(:),   INTENT(OUT)    ::  PENTR_CLD ! Mass flux entrainment of the updraft in cloudy part
REAL, DIMENSION(:),   INTENT(OUT)    ::  PDETR_CLD ! Mass flux detrainment of the updraft in cloudy part
REAL, DIMENSION(:),   INTENT(OUT)    ::  PBUO_INTEG_DRY,PBUO_INTEG_CLD   ! Integrated Buoyancy
REAL, DIMENSION(:),   INTENT(OUT)    ::  PPART_DRY
!
!
END SUBROUTINE COMPUTE_ENTR_DETR

END INTERFACE
!
END MODULE MODI_COMPUTE_ENTR_DETR
