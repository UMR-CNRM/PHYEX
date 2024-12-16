!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################################################################
SUBROUTINE LIMA_ADJUST_SPLIT(D, CST, BUCONF, TBUDGETS, KBUDGETS,                &
                             KRR, KMI, HCONDENS, HLAMBDA3,                      &
                             OSUBG_COND, OSIGMAS, PTSTEP, PSIGQSAT,             &
                             PRHODREF, PRHODJ, PEXNREF, PSIGS, LMFCONV, PMFCONV,&
                             PPABST, PPABSTT, PZZ, ODTHRAD, PDTHRAD, PW_NU,     &
                             PRT, PRS, PSVT, PSVS,                              &
                             HACTCCN, PAERO,PSOLORG, PMI,                       &
                             PTHS, OCOMPUTE_SRC, PSRCS, PCLDFR, PICEFR,         &
                             PRC_MF, PRI_MF, PCF_MF)
!     ###########################################################################
!
!!****  *MIMA_ADJUST* -  compute the fast microphysical sources 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the fast microphysical sources
!!      through an explict scheme and a saturation ajustement procedure.
!!
!!
!!**  METHOD
!!    ------
!!      Reisin et al.,    1996 for the explicit scheme when ice is present
!!      Langlois, Tellus, 1973 for the implict adjustment for the cloud water
!!      (refer also to book 1 of the documentation).
!!
!!      Computations are done separately for three cases :
!!        - ri>0 and rc=0
!!        - rc>0 and ri=0
!!        - ri>0 and rc>0
!!
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!         XP00               ! Reference pressure
!!         XMD,XMV            ! Molar mass of dry air and molar mass of vapor
!!         XRD,XRV            ! Gaz constant for dry air, gaz constant for vapor
!!         XCPD,XCPV          ! Cpd (dry air), Cpv (vapor)
!!         XCL                ! Cl (liquid)
!!         XTT                ! Triple point temperature
!!         XLVTT              ! Vaporization heat constant
!!         XALPW,XBETAW,XGAMW ! Constants for saturation vapor 
!!                            !  pressure  function 
!!      Module  MODD_CONF 
!!         CCONF
!!      Module MODD_BUDGET:
!!         NBUMOD 
!!         CBUTYPE
!!         LBU_RTH    
!!         LBU_RRV  
!!         LBU_RRC  
!!      Module MODD_LES : NCTR_LES,LTURB_LES,NMODNBR_LES
!!                        XNA declaration (cloud fraction as global var)
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 1 and Book2 of documentation ( routine FAST_TERMS )
!!      Langlois, Tellus, 1973
!!
!!    AUTHOR
!!    ------
!!      E. Richard       * Laboratoire d'Aerologie*
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             06/2021 forked from lima_adjust.f90 
!  P. Wautelet 23/07/2021: replace non-standard FLOAT function by REAL function
!  B. Vie         03/2022: Add option for 1-moment pristine ice
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET,   ONLY: TBUDGETDATA, TBUDGETCONF_t, NBUDGET_TH, NBUDGET_RV, &
                         NBUDGET_RC, NBUDGET_RI, NBUDGET_RV, NBUDGET_SV1, NBUMOD        
USE MODD_CST,            ONLY: CST_t
USE MODD_CH_AEROSOL,      ONLY: NSP,NCARB,NSOA
!USE MODD_CONF
!use modd_field,            only: TFIELDDATA, TYPEREAL
!USE MODD_IO,               ONLY: TFILEDATA
!USE MODD_LUNIT_n,          ONLY: TLUOUT
USE MODD_NSV
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_COLD
USE MODD_PARAM_LIMA_MIXED
USE MODD_PARAM_LIMA_WARM
USE MODD_RAIN_ICE_PARAM_n,   ONLY: RAIN_ICE_PARAMN
USE MODD_NEB_n,            ONLY: NEBN
USE MODD_TURB_n,           ONLY: TURBN
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_t
!
USE MODE_BUDGET_PHY,       ONLY: BUDGET_STORE_INIT_PHY, BUDGET_STORE_END_PHY
!USE MODE_IO_FIELD_WRITE,   only: IO_Field_write
use mode_msg
!
USE MODI_CONDENSATION
USE MODE_LIMA_CCN_ACTIVATION, ONLY: LIMA_CCN_ACTIVATION
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
TYPE(DIMPHYEX_t),         INTENT(IN)   :: D
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
TYPE(TBUDGETDATA), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
INTEGER,                  INTENT(IN)   :: KRR        ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KMI        ! Model index 
CHARACTER(len=80),        INTENT(IN)   :: HCONDENS
CHARACTER(len=4),         INTENT(IN)   :: HLAMBDA3   ! formulation for lambda3 coeff
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid
                                                     ! Condensation
LOGICAL,                  INTENT(IN)   :: OSIGMAS    ! Switch for Sigma_s:
                                                     ! use values computed in CONDENSATION
                                                     ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PTSTEP     ! Time step
REAL,                     INTENT(IN)   :: PSIGQSAT   ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PRHODREF  ! Dry density of the 
                                                                   ! reference state
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PRHODJ    ! Dry density * Jacobian
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PEXNREF   ! Reference Exner function
REAL, DIMENSION(MERGE(D%NIT,0,NEBN%LSUBG_COND), &
                MERGE(D%NJT,0,NEBN%LSUBG_COND), &
                MERGE(D%NKT,0,NEBN%LSUBG_COND)),   INTENT(IN)   ::  PSIGS     ! Sigma_s at time t
LOGICAL,                                  INTENT(IN)    ::  LMFCONV ! T to use PMFCONV
REAL, DIMENSION(MERGE(D%NIT,0,LMFCONV), &
                MERGE(D%NJT,0,LMFCONV), &
                MERGE(D%NKT,0,LMFCONV)),   INTENT(IN)   ::  PMFCONV   ! 
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PPABST    ! Absolute Pressure at t     
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PPABSTT   ! Absolute Pressure at t+dt     
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)   ::  PZZ       !     
LOGICAL,                                INTENT(IN)   :: ODTHRAD    ! Use radiative temperature tendency
REAL, DIMENSION(MERGE(D%NIT,0,ODTHRAD), &
                MERGE(D%NJT,0,ODTHRAD), &
                MERGE(D%NKT,0,ODTHRAD)),   INTENT(IN) :: PDTHRAD   ! Radiative temperature tendency
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(IN)    :: PW_NU     ! updraft velocity used for
REAL, DIMENSION(D%NIT, D%NJT, D%NKT ,NSV), INTENT(INOUT) :: PAERO    ! Aerosol concentration
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, 10),  INTENT(IN)    :: PSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, NSP+NCARB+NSOA), INTENT(IN)    :: PMI
CHARACTER(LEN=4),         INTENT(IN)    :: HACTCCN  ! kind of CCN activation
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, KRR), INTENT(IN)    :: PRT       ! m.r. at t
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, KRR), INTENT(INOUT) :: PRS       ! m.r. source
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, NSV), INTENT(IN)    :: PSVT ! Concentrations at time t
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT, NSV), INTENT(INOUT) :: PSVS ! Concentration sources
!
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(INOUT) :: PTHS      ! Theta source
!
LOGICAL,                                      INTENT(IN)    :: OCOMPUTE_SRC ! T to comput PSRCS
REAL, DIMENSION(MERGE(D%NIT,0,OCOMPUTE_SRC), &
                MERGE(D%NJT,0,OCOMPUTE_SRC), &
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(OUT)   :: PSRCS     ! Second-order flux
                                                                         ! s'rc'/2Sigma_s2 at time t+1
                                                                         ! multiplied by Lambda_3
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(INOUT)   :: PCLDFR    ! Cloud fraction          
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),   INTENT(INOUT)   :: PICEFR    ! Cloud fraction          
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),     INTENT(IN)    :: PRI_MF! Convective Mass Flux ice mixing ratio
REAL, DIMENSION(D%NIT, D%NJT, D%NKT),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
  !
  !
  !*       0.2   Declarations of local variables :
  !
  ! 3D Microphysical variables
  REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) &
       :: PTHT,        &
       PRVT,        & ! Water vapor m.r. at t
       PRCT,        & ! Cloud water m.r. at t
       PRRT,        & ! Rain water m.r. at t
       PRIT,        & ! Cloud ice  m.r. at t
                                !
       PRVS,        & ! Water vapor m.r. source
       PRCS,        & ! Cloud water m.r. source
       PRRS,        & ! Rain water m.r. source
       PRIS,        & ! Cloud ice  m.r. source
       PRSS,        & ! Aggregate  m.r. source
       PRGS,        & ! Graupel    m.r. source
       PRHS,        & ! Hail       m.r. source
                                !
       PCCT,        & ! Cloud water conc. at t
       PCIT,        & ! Cloud ice   conc. at t
                                !
       PCCS,        & ! Cloud water C. source
       PMAS           ! Mass of scavenged AP
  !
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE &
       :: PNFS,        & ! Free      CCN C. source
       PNAS,        & ! Activated CCN C. source
       PNFT,        & ! Free      CCN C.
       PNAT           ! Activated CCN C.
  !
  !
  !
  REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) &
       :: ZEXNS,&      ! guess of the Exner function at t+1
       ZT, ZT2,  &      ! guess of the temperature at t+1
       ZCPH, &      ! guess of the CPh for the mixing
       ZW,   &
       ZW1,  &
       ZW2,  &
       ZLV,  &      ! guess of the Lv at t+1
       ZLS,  &      ! guess of the Ls at t+1
       ZMASK,&
       ZRV, ZRV2,ZRV_IN,  &
       ZRC, ZRC2,ZRC_IN,  &
       ZRI, ZRI_IN,  &
       Z_SIGS, Z_SRCS, &
       ZW_MF, &
       ZCND, ZS, ZVEC1, ZDUM
  REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2)) :: ZSIGQSAT2D
  !
  INTEGER, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) :: IVEC1
  !
  INTEGER                           :: ISIZE
  LOGICAL                           :: G_SIGMAS, GUSERI
  REAL, DIMENSION(:), ALLOCATABLE   :: ZRTMIN
  REAL, DIMENSION(:), ALLOCATABLE   :: ZCTMIN
  !
  integer :: idx
  integer :: JI, JJ, JK, jl
  INTEGER                           :: JMOD
  !
  INTEGER :: ISV_LIMA_NC
  INTEGER :: ISV_LIMA_CCN_FREE
  INTEGER :: ISV_LIMA_CCN_ACTI
  INTEGER :: ISV_LIMA_SCAVMASS
  !
  !-------------------------------------------------------------------------------
  !
  !*       1.     PRELIMINARIES
  !               -------------
  !
  ISV_LIMA_NC       = NSV_LIMA_NC       - NSV_LIMA_BEG + 1
  ISV_LIMA_CCN_FREE = NSV_LIMA_CCN_FREE - NSV_LIMA_BEG + 1
  ISV_LIMA_CCN_ACTI = NSV_LIMA_CCN_ACTI - NSV_LIMA_BEG + 1
  ISV_LIMA_SCAVMASS = NSV_LIMA_SCAVMASS - NSV_LIMA_BEG + 1
  !
  ISIZE = SIZE(XRTMIN)
  ALLOCATE(ZRTMIN(ISIZE))
  ZRTMIN(:) = XRTMIN(:) / PTSTEP
  ISIZE = SIZE(XCTMIN)
  ALLOCATE(ZCTMIN(ISIZE))
  ZCTMIN(:) = XCTMIN(:) / PTSTEP
  !
  ! Prepare 3D water mixing ratios
  !
  PTHT = PTHS*PTSTEP
  !
  PRVT(:,:,:) = PRS(:,:,:,1)*PTSTEP
  PRVS(:,:,:) = PRS(:,:,:,1)
  !
  PRCT(:,:,:) = 0.
  PRCS(:,:,:) = 0.
  PRRT(:,:,:) = 0.
  PRRS(:,:,:) = 0.
  PRIT(:,:,:) = 0.
  PRIS(:,:,:) = 0.
  PRSS(:,:,:) = 0.
  PRGS(:,:,:) = 0.
  PRHS(:,:,:) = 0.
  !
  IF ( KRR .GE. 2 ) PRCT(:,:,:) = PRS(:,:,:,2)*PTSTEP
  IF ( KRR .GE. 2 ) PRCS(:,:,:) = PRS(:,:,:,2)
  IF ( KRR .GE. 3 ) PRRT(:,:,:) = PRT(:,:,:,3) 
  IF ( KRR .GE. 3 ) PRRS(:,:,:) = PRS(:,:,:,3)
  IF ( KRR .GE. 4 ) PRIT(:,:,:) = PRT(:,:,:,4)
  IF ( KRR .GE. 4 ) PRIS(:,:,:) = PRS(:,:,:,4) 
  IF ( KRR .GE. 5 ) PRSS(:,:,:) = PRS(:,:,:,5) 
  IF ( KRR .GE. 6 ) PRGS(:,:,:) = PRS(:,:,:,6)
  IF ( KRR .GE. 7 ) PRHS(:,:,:) = PRS(:,:,:,7)
  !
  ! Prepare 3D number concentrations
  !
  PCCT(:,:,:) = 0.
  PCCS(:,:,:) = 0.
  !
  IF ( NMOM_C.GE.2 ) PCCT(:,:,:) = PSVS(:,:,:,ISV_LIMA_NC)*PTSTEP
  IF ( NMOM_C.GE.2 ) PCCS(:,:,:) = PSVS(:,:,:,ISV_LIMA_NC)
  !
  IF ( LSCAV .AND. LAERO_MASS ) PMAS(:,:,:) = PSVS(:,:,:,ISV_LIMA_SCAVMASS)
  ! 
  IF ( NMOM_C.GE.1 .AND. NMOD_CCN.GE.1 ) THEN
     ALLOCATE( PNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
     ALLOCATE( PNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
     ALLOCATE( PNFT(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
     ALLOCATE( PNAT(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
     PNFS(:,:,:,:) = PSVS(:,:,:,ISV_LIMA_CCN_FREE:ISV_LIMA_CCN_FREE+NMOD_CCN-1)
     PNAS(:,:,:,:) = PSVS(:,:,:,ISV_LIMA_CCN_ACTI:ISV_LIMA_CCN_ACTI+NMOD_CCN-1)
     PNFT(:,:,:,:) = PSVS(:,:,:,ISV_LIMA_CCN_FREE:ISV_LIMA_CCN_FREE+NMOD_CCN-1)*PTSTEP
     PNAT(:,:,:,:) = PSVS(:,:,:,ISV_LIMA_CCN_ACTI:ISV_LIMA_CCN_ACTI+NMOD_CCN-1)*PTSTEP
  END IF
  !
  ! Initialize budgets
  !
  if ( nbumod == kmi .and. BUCONF%lbu_enable ) then
     if ( BUCONF%lbudget_th ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'CEDS', pths(:, :, :) * prhodj(:, :, :) )
     if ( BUCONF%lbudget_rv ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'CEDS', prvs(:, :, :) * prhodj(:, :, :) )
     if ( BUCONF%lbudget_rc ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'CEDS', prcs(:, :, :) * prhodj(:, :, :) )
     if ( BUCONF%lbudget_ri ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RI), 'CEDS', pris(:, :, :) * prhodj(:, :, :) )
     if ( BUCONF%lbudget_sv ) then
        if ( nmom_c.ge.2) &
             call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nc), 'CEDS', pccs(:, :, :) * prhodj(:, :, :) )
        if ( lscav .and. laero_mass ) &
             call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_scavmass), 'CEDS', pmas(:, :, :) * prhodj(:, :, :) )
        if ( nmom_c.ge.1 ) then
           do jl = 1, nmod_ccn
              idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_free - 1 + jl
              call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'CEDS', pnfs(:, :, :, jl) * prhodj(:, :, :) )
              idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_acti - 1 + jl
              call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'CEDS', pnas(:, :, :, jl) * prhodj(:, :, :) )
           end do
        end if
     end if
  end if
  !
  !-------------------------------------------------------------------------------
  !
  !*       2.     CONDENSATION
  !               ------------
  !
  WHERE ( PRVS(:,:,:)+PRCS(:,:,:)+PRIS(:,:,:) < 0.)
     PRVS(:,:,:) = -  PRCS(:,:,:) - PRIS(:,:,:)
  END WHERE
  !
  ZEXNS(:,:,:) = ( PPABSTT(:,:,:) / CST%XP00 ) ** (CST%XRD/CST%XCPD)  
  ZT(:,:,:) = ( PTHS(:,:,:) * PTSTEP ) * ZEXNS(:,:,:)
  ZT2(:,:,:) = ZT(:,:,:)
  ZCPH(:,:,:) = CST%XCPD + CST%XCPV  * PTSTEP *   PRVS(:,:,:)                             &
       + CST%XCL * PTSTEP * ( PRCS(:,:,:) + PRRS(:,:,:) )             &
       + CST%XCI * PTSTEP * ( PRIS(:,:,:) + PRSS(:,:,:) + PRGS(:,:,:) + PRHS(:,:,:) )
  ZLV(:,:,:) = CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZT(:,:,:) -CST%XTT )
  ZLS(:,:,:) = CST%XLSTT + ( CST%XCPV - CST%XCI ) * ( ZT(:,:,:) -CST%XTT )
  !
  IF (LADJ) THEN
     ZRV_IN=PRVS*PTSTEP
     ZRC_IN=PRCS*PTSTEP
     IF (NMOM_I.EQ.1) THEN
        ZRI_IN=PRIS*PTSTEP
        GUSERI=.TRUE.
     ELSE
        ZRI_IN=0.
        GUSERI=.FALSE.
     END IF
     IF (OSUBG_COND) THEN
        Z_SIGS=PSIGS
        G_SIGMAS=OSIGMAS
        ZSIGQSAT2D(:,:)=PSIGQSAT
     ELSE
        Z_SIGS=0.
        G_SIGMAS=.TRUE.
        ZSIGQSAT2D(:,:)=0.
     END IF
     CALL CONDENSATION(D, CST, RAIN_ICE_PARAMN, NEBN, TURBN,                    &
          'S', HCONDENS, HLAMBDA3,                                             &
          PPABST, PZZ, PRHODREF, ZT, ZRV_IN, ZRV, ZRC_IN, ZRC, ZRI_IN, ZRI,   &
          PRRS*PTSTEP,PRSS*PTSTEP, PRGS*PTSTEP, &
          Z_SIGS, .FALSE., PMFCONV, PCLDFR, Z_SRCS, GUSERI, G_SIGMAS, .FALSE., &
          ZDUM, ZDUM, ZDUM, ZDUM, ZDUM,              &
          ZSIGQSAT2D, PLV=ZLV, PLS=ZLS, PCPH=ZCPH )
     IF (NMOM_C.GE.1) THEN
        ZW1(:,:,:) = (ZRC(:,:,:) - PRCS(:,:,:)*PTSTEP) / PTSTEP ! Pcon = Delta_rc / Delta_t
        WHERE( ZW1(:,:,:) < 0.0 )
           ZW1(:,:,:) = MAX ( ZW1(:,:,:), -PRCS(:,:,:) )
        ELSEWHERE
           ZW1(:,:,:) = MIN ( ZW1(:,:,:),  PRVS(:,:,:) )
        END WHERE
        PRVS(:,:,:) = PRVS(:,:,:) - ZW1(:,:,:)
        PRCS(:,:,:) = PRCS(:,:,:) + ZW1(:,:,:)
        PTHS(:,:,:) = PTHS(:,:,:) +        &
             ZW1(:,:,:) * ZLV(:,:,:) / (ZCPH(:,:,:) * PEXNREF(:,:,:))
     END IF
     IF (NMOM_I.EQ.1) THEN
        PICEFR(:,:,:)=PCLDFR(:,:,:)
        ZW2(:,:,:) = (ZRI(:,:,:) - PRIS(:,:,:)*PTSTEP) / PTSTEP ! Pdep = Delta_ri / Delta_t
        !
        WHERE( ZW2(:,:,:) < 0.0 )
           ZW2(:,:,:) = MAX ( ZW2(:,:,:), -PRIS(:,:,:) )
        ELSEWHERE
           ZW2(:,:,:) = MIN ( ZW2(:,:,:),  PRVS(:,:,:) )
        END WHERE
        PRVS(:,:,:) = PRVS(:,:,:) - ZW2(:,:,:)
        PRIS(:,:,:) = PRIS(:,:,:) + ZW2(:,:,:)
        PTHS(:,:,:) = PTHS(:,:,:) +        &
             ZW2(:,:,:) * ZLS(:,:,:) / (ZCPH(:,:,:) * PEXNREF(:,:,:))
     END IF
  ELSE
     DO JI=1,SIZE(PRCS,1)
        DO JJ=1,SIZE(PRCS,2)
           DO JK=1,SIZE(PRCS,3)
              IF (PRCS(JI,JJ,JK).GE.XRTMIN(2) .AND. PCCS(JI,JJ,JK).GE.XCTMIN(2)) THEN
                 ZVEC1(JI,JJ,JK) = MAX( 1.0001, MIN( FLOAT(NAHEN)-0.0001, XAHENINTP1 * ZT(JI,JJ,JK) + XAHENINTP2 ) )
                 IVEC1(JI,JJ,JK) = INT( ZVEC1(JI,JJ,JK) )
                 ZVEC1(JI,JJ,JK) = ZVEC1(JI,JJ,JK) - FLOAT( IVEC1(JI,JJ,JK) )
                 ZW(JI,JJ,JK)=EXP( CST%XALPW - CST%XBETAW/ZT(JI,JJ,JK) - CST%XGAMW*ALOG(ZT(JI,JJ,JK) ) ) ! es_w
                 ZW(JI,JJ,JK)=CST%XMV / CST%XMD * ZW(JI,JJ,JK) / ( PPABST(JI,JJ,JK)-ZW(JI,JJ,JK) ) 
                 ZS(JI,JJ,JK) = PRVS(JI,JJ,JK)*PTSTEP / ZW(JI,JJ,JK) - 1.
                 ZW(JI,JJ,JK) = PCCS(JI,JJ,JK)*PTSTEP/(XLBC*PCCS(JI,JJ,JK)/PRCS(JI,JJ,JK))**XLBEXC
                 ZW2(JI,JJ,JK) = XAHENG3(IVEC1(JI,JJ,JK)+1)*ZVEC1(JI,JJ,JK)-XAHENG3(IVEC1(JI,JJ,JK))*(ZVEC1(JI,JJ,JK)-1.)
                 ZCND(JI,JJ,JK) = 2.*3.14*1000.*ZW2(JI,JJ,JK)*ZS(JI,JJ,JK)*ZW(JI,JJ,JK)
                 IF(ZCND(JI,JJ,JK).LE.0.) THEN
                    ZCND(JI,JJ,JK) = MAX ( ZCND(JI,JJ,JK), -PRCS(JI,JJ,JK) )
                 ELSE
                    ZCND(JI,JJ,JK) = MIN ( ZCND(JI,JJ,JK),  PRVS(JI,JJ,JK) )
                 END IF
                 PRVS(JI,JJ,JK) = PRVS(JI,JJ,JK) - ZCND(JI,JJ,JK)
                 PRCS(JI,JJ,JK) = PRCS(JI,JJ,JK) + ZCND(JI,JJ,JK)
                 PTHS(JI,JJ,JK) = PTHS(JI,JJ,JK) + ZCND(JI,JJ,JK) * ZLV(JI,JJ,JK) / (ZCPH(JI,JJ,JK) * PEXNREF(JI,JJ,JK))
              END IF
           END DO
        END DO
     END DO
  END IF
  !
  IF (OSUBG_COND .AND. NMOM_C.GE.2 .AND. LACTI) THEN
     PSRCS=Z_SRCS
     ZW_MF=0.
     ZRV2=PRVT
     ZRC2=PRCT
     CALL LIMA_CCN_ACTIVATION (CST,                          &
          PRHODREF, PEXNREF, PPABST, ZT2, PDTHRAD, PW_NU+ZW_MF, &
          PAERO,PSOLORG, PMI,  HACTCCN,                         &
          PTHT, ZRV2, ZRC2, PCCT, PRRT, PNFT, PNAT,             &
          PCLDFR                                                )      
  END IF
  !
  !-------------------------------------------------------------------------------
  !
  !*       3.     CLOUD FROM MASS-FLUX SCHEME
  !               ---------------------------
  !
  IF ( .NOT. OSUBG_COND ) THEN
     WHERE (PRCS(:,:,:) + PRIS(:,:,:) > 1.E-12 / PTSTEP)
        PCLDFR(:,:,:)  = 1.
     ELSEWHERE
        PCLDFR(:,:,:)  = 0. 
     ENDWHERE
     IF ( SIZE(PSRCS,3) /= 0 ) THEN
        PSRCS(:,:,:)  = PCLDFR(:,:,:)
     END IF
  ELSE
     ! We limit PRC_MF+PRI_MF to PRVS*PTSTEP to avoid negative humidity
     ZW1(:,:,:)=PRC_MF(:,:,:)/PTSTEP
     ZW2(:,:,:)=0.
     IF (NMOM_I.EQ.1) ZW2(:,:,:)=PRI_MF(:,:,:)/PTSTEP
     WHERE(ZW1(:,:,:)+ZW2(:,:,:)>PRVS(:,:,:))
        ZW1(:,:,:)=ZW1(:,:,:)*PRVS(:,:,:)/(ZW1(:,:,:)+ZW2(:,:,:))
        ZW2(:,:,:)=PRVS(:,:,:)-ZW1(:,:,:)
     ENDWHERE
     !
     PCLDFR(:,:,:) = MIN(1.,PCLDFR(:,:,:)+PCF_MF(:,:,:))
     PRVS(:,:,:)   = PRVS(:,:,:) - ZW1(:,:,:) -ZW2(:,:,:)
     PRCS(:,:,:)   = PRCS(:,:,:) + ZW1(:,:,:)
     IF (NMOM_I.EQ.1) PRIS(:,:,:)   = PRIS(:,:,:) + ZW2(:,:,:)
     IF (NMOM_C.GE.2) PCCS(:,:,:)   = PCCT(:,:,:) / PTSTEP
     IF (NMOD_CCN.GE.1) PNFS(:,:,:,:) = PNFT(:,:,:,:) / PTSTEP
     IF (NMOD_CCN.GE.1) PNAS(:,:,:,:) = PNAT(:,:,:,:) / PTSTEP
     PTHS(:,:,:)   = PTHS(:,:,:) + &
          (ZW1(:,:,:) * ZLV(:,:,:) + ZW2 * ZLS(:,:,:)) / ZCPH(:,:,:)     &
          /  PEXNREF(:,:,:)
  END IF
  !
  !-------------------------------------------------------------------------------
  !
  !*       4.     REMOVE SMALL NUMBERS OF DROPLETS
  !               --------------------------------
  !
  IF (NMOM_C .GE. 2) THEN
     ZMASK(:,:,:) = 0.0
     ZW(:,:,:) = 0.
     WHERE (PRCS(:,:,:) <= ZRTMIN(2) .OR. PCCS(:,:,:) <=0.) 
        PRVS(:,:,:) = PRVS(:,:,:) + PRCS(:,:,:) 
        PTHS(:,:,:) = PTHS(:,:,:) - PRCS(:,:,:)*ZLV(:,:,:)/(ZCPH(:,:,:)*ZEXNS(:,:,:))
        PRCS(:,:,:) = 0.0
        ZW(:,:,:)   = MAX(PCCS(:,:,:),0.)
        PCCS(:,:,:) = 0.0
     END WHERE
     !
     ZW1(:,:,:) = 0.
     IF (NMOD_CCN.GE.1) ZW1(:,:,:) = SUM(PNAS,DIM=4)
     ZW (:,:,:) = MIN( ZW(:,:,:), ZW1(:,:,:) )
     ZW2(:,:,:) = 0.
     WHERE ( ZW(:,:,:) > 0. )
        ZMASK(:,:,:) = 1.0
        ZW2(:,:,:) = ZW(:,:,:) / ZW1(:,:,:)
     ENDWHERE
     DO JMOD = 1, NMOD_CCN
        PNFS(:,:,:,JMOD) = PNFS(:,:,:,JMOD) +                           &
             ZMASK(:,:,:) * PNAS(:,:,:,JMOD) * ZW2(:,:,:)
        PNAS(:,:,:,JMOD) = PNAS(:,:,:,JMOD) -                           &
             ZMASK(:,:,:) * PNAS(:,:,:,JMOD) * ZW2(:,:,:)
        PNAS(:,:,:,JMOD) = MAX( 0.0 , PNAS(:,:,:,JMOD) )
     ENDDO
     IF (LSCAV .AND. LAERO_MASS) PMAS(:,:,:) = PMAS(:,:,:) * (1-ZMASK(:,:,:))
  END IF
  !
  !-------------------------------------------------------------------------------
  !
  !*       5.     SAVE CHANGES & CLEANING
  !               -----------------------
  !
  ! 3D water mixing ratios
  PRS(:,:,:,1) = PRVS(:,:,:)
  IF ( KRR .GE. 2 ) PRS(:,:,:,2) = PRCS(:,:,:)
  IF ( KRR .GE. 3 ) PRS(:,:,:,3) = PRRS(:,:,:)
  IF ( KRR .GE. 4 ) PRS(:,:,:,4) = PRIS(:,:,:)
  IF ( KRR .GE. 5 ) PRS(:,:,:,5) = PRSS(:,:,:)
  IF ( KRR .GE. 6 ) PRS(:,:,:,6) = PRGS(:,:,:)
  IF ( KRR .GE. 7 ) PRS(:,:,:,7) = PRHS(:,:,:)
  !
  ! 3D number concentrations
  IF ( NMOM_C.GE.2 ) PSVS(:,:,:,ISV_LIMA_NC) = PCCS(:,:,:)
  IF ( LSCAV .AND. LAERO_MASS ) PSVS(:,:,:,ISV_LIMA_SCAVMASS) = PMAS(:,:,:)
  IF ( NMOD_CCN.GE.1 ) THEN
     PSVS(:,:,:,ISV_LIMA_CCN_FREE:ISV_LIMA_CCN_FREE+NMOD_CCN-1) = PNFS(:,:,:,:)
     PSVS(:,:,:,ISV_LIMA_CCN_ACTI:ISV_LIMA_CCN_ACTI+NMOD_CCN-1) = PNAS(:,:,:,:)
  END IF
  !
  ! budgets
  if ( nbumod == kmi .and. BUCONF%lbu_enable ) then
     if ( BUCONF%lbudget_th ) call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_TH), 'CEDS', pths(:, :, :) * prhodj(:, :, :) )
     if ( BUCONF%lbudget_rv ) call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RV), 'CEDS', prvs(:, :, :) * prhodj(:, :, :) )
     if ( BUCONF%lbudget_rc ) call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RC), 'CEDS', prcs(:, :, :) * prhodj(:, :, :) )
     if ( BUCONF%lbudget_ri ) call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_RI), 'CEDS', pris(:, :, :) * prhodj(:, :, :) )
     if ( BUCONF%lbudget_sv ) then
        if ( nmom_c.ge.2) &
             call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_nc), 'CEDS', pccs(:, :, :) * prhodj(:, :, :) )
        if ( lscav .and. laero_mass ) &
             call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_scavmass), 'CEDS', pmas(:, :, :) * prhodj(:, :, :) )
        if ( nmom_c.ge.1 ) then
           do jl = 1, nmod_ccn
              idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_free - 1 + jl
              call BUDGET_STORE_END_PHY(D, TBUDGETS(idx), 'CEDS', pnfs(:, :, :, jl) * prhodj(:, :, :) )
              idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_acti - 1 + jl
              call BUDGET_STORE_END_PHY(D, TBUDGETS(idx), 'CEDS', pnas(:, :, :, jl) * prhodj(:, :, :) )
           end do
        end if
     end if
  end if
  !
  DEALLOCATE(ZRTMIN)
  DEALLOCATE(ZCTMIN)
  IF (ALLOCATED(PNFS)) DEALLOCATE(PNFS)
  IF (ALLOCATED(PNAS)) DEALLOCATE(PNAS)
  IF (ALLOCATED(PNFT)) DEALLOCATE(PNFT)
  IF (ALLOCATED(PNAT)) DEALLOCATE(PNAT)
  !
  !------------------------------------------------------------------------------
  !
END SUBROUTINE LIMA_ADJUST_SPLIT
