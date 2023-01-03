!     ######spl
       MODULE MODI_RAIN_ICE_OLD
!      ####################
!
INTERFACE
      SUBROUTINE RAIN_ICE_OLD ( OSEDIC, OCND2, LKOGAN, LMODICEDEP, HSEDIM, HSUBG_AUCV_RC, OWARM,    &
                            KKA, KKU, KKL,                                        &
                            KSPLITR, PTSTEP, KRR,                            &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                            PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR,                 &
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
                            PINPRC, PINPRR, PEVAP3D,                    &
                            PINPRS, PINPRG, PSIGS, PSEA, PTOWN,                   &
                            YDDDH, YDLDDH, YDMDDH, &
                            PICENU, PKGN_ACON, PKGN_SBGR, &
                            PRHT, PRHS, PINPRH, PFPR)
!
USE DDH_MIX, ONLY : TYP_DDH
USE YOMLDDH, ONLY : TLDDH
USE YOMMDDH, ONLY : TMDDH
!
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
LOGICAL,                  INTENT(IN)    :: OCND2  ! Logical switch to separate liquid and ice
LOGICAL,                  INTENT(IN)    :: LKOGAN ! Logical switch for using Kogan autoconversion of liquid.
LOGICAL,                  INTENT(IN)    :: LMODICEDEP ! Logical switch for alternative dep/evap of ice
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Switch for rc->rr Subgrid autoconversion
                                        ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
!
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step
                                      ! integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ    ! Layer thickness (m)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
! input from aro_adjust / condensation with OCND2, dummy if OCND2 = F
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PICLDFR ! ice cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PWCLDFR ! water or mixed-phase cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PSSIO   ! Super-saturation with respect to ice in the  
                                                   ! supersaturated fraction
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PSSIU   ! Sub-saturation with respect to ice in the  
                                                   ! subsaturated fraction 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PIFR    ! Ratio cloud ice moist part to dry part 
! input from aro_adjust / condensation with OCND2 END.
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRC! Cloud instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRR! Rain instant precip
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PEVAP3D! Rain evap profile
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRS! Snow instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRG! Graupel instant precip
REAL, DIMENSION(:,:), INTENT(IN)        :: PSEA ! Sea Mask
REAL, DIMENSION(:,:), INTENT(IN)        :: PTOWN! Fraction that is town
TYPE(TYP_DDH),        INTENT(INOUT)     :: YDDDH
TYPE(TLDDH),          INTENT(IN)        :: YDLDDH
TYPE(TMDDH),          INTENT(IN)        :: YDMDDH
REAL, DIMENSION(:,:), INTENT(IN)        :: PICENU, PKGN_ACON, PKGN_SBGR
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(:,:,:),   OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(:,:),     OPTIONAL, INTENT(INOUT) :: PINPRH  ! Hail instant precip
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)   :: PFPR    ! upper-air precipitation fluxes
!
END SUBROUTINE RAIN_ICE_OLD
END INTERFACE
END MODULE MODI_RAIN_ICE_OLD
