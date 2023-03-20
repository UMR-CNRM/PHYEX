!     ######spl
       MODULE MODI_RAIN_ICE_OLD
!      ####################
!
INTERFACE
      SUBROUTINE RAIN_ICE_OLD (D, CST, PARAMI, ICEP, ICED, BUCONF,                    &
                               OSEDIC, OCND2, LKOGAN, LMODICEDEP,                     &
                               HSEDIM, HSUBG_AUCV_RC, OWARM,                          &
                               KKA, KKU, KKL,                                         &
                               KSPLITR, PTSTEP, KRR, KSIZE, GMICRO,                   &
                               PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                               PICLDFR, PSSIO, PSSIU, PIFR,                           &
                               PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                    &
                               PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,        &
                               PINPRC, PINPRR, PEVAP3D,                               &
                               PINPRS, PINPRG, PSIGS, PSEA, PTOWN,                    &
                               TBUDGETS, KBUDGETS,                                    &
                               PICENU, PKGN_ACON, PKGN_SBGR,                          &
                               PRHT, PRHS, PINPRH, PFPR)
!
USE MODD_BUDGET,         ONLY: TBUDGETDATA, TBUDGETCONF_t
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_CST, ONLY: CST_T
USE MODD_PARAM_ICE_n,    ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_T
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_T
!
TYPE(DIMPHYEX_T),       INTENT(IN) :: D
TYPE(CST_T),            INTENT(IN) :: CST 
TYPE(PARAM_ICE_t),      INTENT(IN) :: PARAMI
TYPE(RAIN_ICE_PARAM_t), INTENT(IN) :: ICEP
TYPE(RAIN_ICE_DESCR_t), INTENT(IN) :: ICED
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF

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
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step
                                      ! integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
!
LOGICAL, DIMENSION(D%NIT,D%NKT), INTENT(IN) :: GMICRO    ! Layer thickness (m)

INTEGER, INTENT(IN) :: KSIZE
!
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PDZZ    ! Layer thickness (m)
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PCLDFR  ! Cloud fraction
!
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PSIGS   ! Sigma_s at t
! input from aro_adjust / condensation with OCND2, dummy if OCND2 = F
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PICLDFR ! ice cloud fraction
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PSSIO   ! Super-saturation with respect to ice in the
                                                 ! supersaturated fraction
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)    :: PSSIU   ! Sub-saturation with respect to ice in the
                                                 ! subsaturated fraction
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PIFR    ! Ratio cloud ice moist part to dry part
! input from aro_adjust / condensation with OCND2 END.
!
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(D%NIT,D%NKT), INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRC! Cloud instant precip
REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRR! Rain instant precip
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT) :: PEVAP3D! Rain evap profile
REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRS! Snow instant precip
REAL, DIMENSION(D%NIT),       INTENT(OUT) :: PINPRG! Graupel instant precip
REAL, DIMENSION(D%NIT),       INTENT(IN)  :: PSEA ! Sea Mask
REAL, DIMENSION(D%NIT),       INTENT(IN)  :: PTOWN! Fraction that is town
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
REAL, DIMENSION(D%NIT), INTENT(IN)            :: PICENU, PKGN_ACON, PKGN_SBGR
REAL, DIMENSION(D%NIT,D%NKT),   OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(D%NIT,D%NKT),   OPTIONAL, INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(D%NIT),         OPTIONAL, INTENT(OUT)   :: PINPRH  ! Hail instant precip
REAL, DIMENSION(D%NIT,D%NKT,KRR), OPTIONAL, INTENT(OUT) :: PFPR    ! upper-air precipitation fluxes
!
END SUBROUTINE RAIN_ICE_OLD
END INTERFACE
END MODULE MODI_RAIN_ICE_OLD
