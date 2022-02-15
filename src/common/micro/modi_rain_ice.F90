!     ######spl
       MODULE MODI_RAIN_ICE
!      ####################
!
INTERFACE
      SUBROUTINE RAIN_ICE ( KPROMA, KIT, KJT, KKT, KSIZE,                                 &
                            OSEDIC,OCND2, HSEDIM, HSUBG_AUCV_RC, HSUBG_AUCV_RI,   &
                            OWARM, KKA, KKU, KKL,   &
                            PTSTEP, KRR, LDMICRO, PEXN,             &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                            PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,               &
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                   &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                            PINPRC, PINPRR, PEVAP3D,                    &
                            PINPRS, PINPRG, PINDEP, PRAINFR, PSIGS, PSEA, PTOWN,            &
                            PRHT, PRHS, PINPRH, PFPR, &
                            TBUDGETS, KBUDGETS)
!
USE MODD_BUDGET, ONLY : TBUDGETDATA
USE MODD_PARAM_ICE,      ONLY: LDEPOSC
!
INTEGER,                  INTENT(IN)    :: KPROMA ! cache-blocking factor for microphysic loop
INTEGER,                  INTENT(IN)    :: KIT, KJT, KKT ! arrays size
INTEGER,                  INTENT(IN)    :: KSIZE
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
LOGICAL                                 :: OCND2  ! Logical switch to separate liquid and ice
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Switch for rc->rr Subgrid autoconversion
                                        ! Kind of Subgrid autoconversion method
CHARACTER(LEN=80),        INTENT(IN)    :: HSUBG_AUCV_RI ! Switch for ri->rs Subgrid autoconversion
                                        ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
!
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)   :: LDMICRO ! mask to limit computation
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLC_HRC
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLC_HCF
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLI_HRI
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PHLI_HCF
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
!
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source

!
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(KIT,KJT,KKT), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(KIT,KJT), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(MERGE(KIT, 0, LDEPOSC), MERGE(KJT, 0, LDEPOSC)),     INTENT(OUT)       :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(OUT) :: PRAINFR !Precipitation fraction
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(IN)        :: PSEA ! Sea Mask
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(IN)        :: PTOWN! Fraction that is town
REAL, DIMENSION(KIT,KJT,KKT), OPTIONAL,  INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(KIT,KJT,KKT), OPTIONAL,  INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(KIT,KJT), OPTIONAL, INTENT(OUT)      :: PINPRH! Hail instant precip
REAL, DIMENSION(KIT,KJT,KKT,KRR), OPTIONAL, INTENT(OUT)  :: PFPR ! upper-air precipitation fluxes
!
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
END SUBROUTINE RAIN_ICE
END INTERFACE
END MODULE MODI_RAIN_ICE
