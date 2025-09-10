!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################################################################
SUBROUTINE LIMA_ADJUST_SPLIT(LIMAP, LIMAW, TNSV, D, CST, NEBN, TURBN, BUCONF, TBUDGETS, KBUDGETS,  &
                             KRR, HCONDENS, HLAMBDA3,                           &
                             KCARB, KSOA, KSP, ODUST, OSALT, OORILAM,           &
                             OSUBG_COND, OSIGMAS, PTSTEP, PSIGQSAT,             &
                             PRHODREF, PRHODJ, PEXNREF, PSIGS, OMFCONV, PMFCONV,&
                             PPABST, PZZ, ODTHRAD, PDTHRAD, PW_NU,              &
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
USE MODD_BUDGET,   ONLY: TBUDGETCONF_T, NBUDGET_TH, NBUDGET_RV,& 
                         NBUDGET_RC, NBUDGET_RI, NBUDGET_RV, NBUDGET_SV1, NBUMOD, TBUDGETDATA_PTR
USE MODD_CST,            ONLY: CST_T
!USE MODD_CONF
!use modd_field,            only: TFIELDDATA, TYPEREAL
!USE MODD_IO,               ONLY: TFILEDATA
!USE MODD_LUNIT_n,          ONLY: TLUOUT
USE MODD_NSV,             ONLY: NSV_T
USE MODD_PARAM_LIMA, ONLY: PARAM_LIMA_T
USE MODD_PARAM_LIMA_WARM, ONLY : PARAM_LIMA_WARM_T
USE MODD_RAIN_ICE_PARAM_N,   ONLY: RAIN_ICE_PARAMN
USE MODD_NEB_N,            ONLY: NEB_T
USE MODD_TURB_N,           ONLY: TURB_T
USE MODD_DIMPHYEX,         ONLY: DIMPHYEX_T
!
!USE MODE_IO_FIELD_WRITE,   only: IO_Field_write
USE MODE_MSG
!
USE MODI_CONDENSATION
USE MODE_LIMA_CCN_ACTIVATION, ONLY: LIMA_CCN_ACTIVATION
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(NSV_T),              INTENT(IN)    :: TNSV
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
TYPE(CST_T),              INTENT(IN)    :: CST
TYPE(NEB_T),              INTENT(IN)    :: NEBN
TYPE(TURB_T),             INTENT(IN)    :: TURBN
TYPE(TBUDGETCONF_T),      INTENT(IN)    :: BUCONF
TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
INTEGER,                  INTENT(IN)   :: KRR        ! Number of moist variables
CHARACTER(LEN=80),        INTENT(IN)   :: HCONDENS
CHARACTER(LEN=4),         INTENT(IN)   :: HLAMBDA3   ! formulation for lambda3 coeff
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid
                                                     ! Condensation
LOGICAL,                  INTENT(IN)   :: OSIGMAS    ! Switch for Sigma_s:
                                                     ! use values computed in CONDENSATION
                                                     ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PTSTEP     ! Time step
REAL, DIMENSION(D%NIJT),  INTENT(IN)   :: PSIGQSAT   ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)   ::  PRHODREF  ! Dry density of the 
                                                             ! reference state
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)   ::  PRHODJ    ! Dry density * Jacobian
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)   ::  PEXNREF   ! Reference Exner function
REAL, DIMENSION(MERGE(D%NIJT,0,NEBN%LSUBG_COND), &
                MERGE(D%NKT,0,NEBN%LSUBG_COND)),   INTENT(IN)   ::  PSIGS     ! Sigma_s at time t
LOGICAL,                                  INTENT(IN)    ::  OMFCONV ! T to use PMFCONV
REAL, DIMENSION(MERGE(D%NIJT,0,OMFCONV), &
                MERGE(D%NKT,0,OMFCONV)),   INTENT(IN)   ::  PMFCONV   ! 
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)   ::  PPABST    ! Absolute Pressure at t     
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)   ::  PZZ       !     
LOGICAL,                                INTENT(IN)   :: ODTHRAD    ! Use radiative temperature tendency
REAL, DIMENSION(MERGE(D%NIJT,0,ODTHRAD), &
                MERGE(D%NKT,0,ODTHRAD)),   INTENT(IN) :: PDTHRAD   ! Radiative temperature tendency
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(IN)    :: PW_NU     ! updraft velocity used for
INTEGER,                  INTENT(IN)    :: KCARB, KSOA, KSP ! for array size declarations
LOGICAL,                  INTENT(IN)    :: ODUST, OSALT, OORILAM
REAL, DIMENSION(D%NIJT, D%NKT ,TNSV%NSV), INTENT(INOUT) :: PAERO    ! Aerosol concentration
REAL, DIMENSION(D%NIJT, D%NKT, 10),  INTENT(IN)    :: PSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(D%NIJT, D%NKT, KSP+KCARB+KSOA), INTENT(IN)    :: PMI
CHARACTER(LEN=4),         INTENT(IN)    :: HACTCCN  ! kind of CCN activation
!
REAL, DIMENSION(D%NIJT, D%NKT, KRR), INTENT(IN)    :: PRT       ! m.r. at t
!
REAL, DIMENSION(D%NIJT, D%NKT, KRR), INTENT(INOUT) :: PRS       ! m.r. source
!
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), INTENT(IN)    :: PSVT ! Concentrations at time t
!
REAL, DIMENSION(D%NIJT, D%NKT, TNSV%NSV), INTENT(INOUT) :: PSVS ! Concentration sources
!
REAL, DIMENSION(D%NIJT, D%NKT),   INTENT(INOUT) :: PTHS      ! Theta source
!
LOGICAL,                                      INTENT(IN)    :: OCOMPUTE_SRC ! T to comput PSRCS
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC), &
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(OUT)   :: PSRCS     ! Second-order flux
                                                                         ! s'rc'/2Sigma_s2 at time t+1
                                                                         ! multiplied by Lambda_3
REAL, DIMENSION(D%NIJT, D%NKT),     INTENT(OUT) :: PCLDFR    ! Cloud fraction          
REAL, DIMENSION(D%NIJT, D%NKT),     INTENT(OUT) :: PICEFR    ! Cloud fraction          
REAL, DIMENSION(D%NIJT, D%NKT),     INTENT(IN)  :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(D%NIJT, D%NKT),     INTENT(IN)  :: PRI_MF! Convective Mass Flux ice mixing ratio
REAL, DIMENSION(D%NIJT, D%NKT),     INTENT(IN)  :: PCF_MF! Convective Mass Flux Cloud fraction 
!
!
!*       0.2   Declarations of local variables :
!
! 3D Microphysical variables
REAL, DIMENSION(D%NIJT,D%NKT) &
                         :: ZTHT,        &
                            ZRVT,        & ! Water vapor m.r. at t
                            ZRCT,        & ! Cloud water m.r. at t
                            ZRRT,        & ! Rain water m.r. at t
                            ZRIT,        & ! Cloud ice  m.r. at t
!
                            ZRVS,        & ! Water vapor m.r. source
                            ZRCS,        & ! Cloud water m.r. source
                            ZRRS,        & ! Rain water m.r. source
                            ZRIS,        & ! Cloud ice  m.r. source
                            ZRSS,        & ! Aggregate  m.r. source
                            ZRGS,        & ! Graupel    m.r. source
                            ZRHS,        & ! Hail       m.r. source
!
                            ZCCT,        & ! Cloud water conc. at t
                            ZCIT,        & ! Cloud ice   conc. at t
!
                            ZCCS,        & ! Cloud water C. source
                            ZMAS           ! Mass of scavenged AP
!
REAL, DIMENSION(:,:,:), ALLOCATABLE &
                         :: ZNFS,        & ! Free      CCN C. source
                            ZNAS,        & ! Activated CCN C. source
                            ZNFT,        & ! Free      CCN C.
                            ZNAT           ! Activated CCN C.
!
REAL, DIMENSION(D%NIJT,D%NKT) &
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
!
INTEGER, DIMENSION(D%NIJT,D%NKT) :: IVEC1
INTEGER                           :: ISIZE
LOGICAL                           :: G_SIGMAS, GUSERI
REAL, DIMENSION(:), ALLOCATABLE   :: ZRTMIN
REAL, DIMENSION(:), ALLOCATABLE   :: ZCTMIN
!
INTEGER :: IDX
INTEGER :: II, IK, IL
INTEGER                           :: IMOD
!
INTEGER :: ISV_LIMA_NC
INTEGER :: ISV_LIMA_CCN_FREE
INTEGER :: ISV_LIMA_CCN_ACTI
INTEGER :: ISV_LIMA_SCAVMASS
INTEGER :: ISV_LIMA_NI
INTEGER :: ISV_LIMA_IFN_FREE
INTEGER :: ISV_LIMA_IFN_NUCL
INTEGER :: ISV_LIMA_IMM_NUCL
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_ADJUST_SPLIT', 0, ZHOOK_HANDLE)
!
ISV_LIMA_NC       = TNSV%NSV_LIMA_NC       - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_CCN_FREE = TNSV%NSV_LIMA_CCN_FREE - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_CCN_ACTI = TNSV%NSV_LIMA_CCN_ACTI - TNSV%NSV_LIMA_BEG + 1
ISV_LIMA_SCAVMASS = TNSV%NSV_LIMA_SCAVMASS - TNSV%NSV_LIMA_BEG + 1
!
ISIZE = SIZE(LIMAP%XRTMIN)
ALLOCATE(ZRTMIN(ISIZE))
ZRTMIN(:) = LIMAP%XRTMIN(:) / PTSTEP
ISIZE = SIZE(LIMAP%XCTMIN)
ALLOCATE(ZCTMIN(ISIZE))
ZCTMIN(:) = LIMAP%XCTMIN(:) / PTSTEP
!
! Prepare 3D water mixing ratios
!
ZTHT = PTHS*PTSTEP
!
ZRVT(:,:) = PRS(:,:,1)*PTSTEP
ZRVS(:,:) = PRS(:,:,1)
!
ZRCT(:,:) = 0.
ZRCS(:,:) = 0.
ZRRT(:,:) = 0.
ZRRS(:,:) = 0.
ZRIT(:,:) = 0.
ZRIS(:,:) = 0.
ZRSS(:,:) = 0.
ZRGS(:,:) = 0.
ZRHS(:,:) = 0.
!
IF ( KRR .GE. 2 ) ZRCT(:,:) = PRS(:,:,2)*PTSTEP
IF ( KRR .GE. 2 ) ZRCS(:,:) = PRS(:,:,2)
IF ( KRR .GE. 3 ) ZRRT(:,:) = PRT(:,:,3) 
IF ( KRR .GE. 3 ) ZRRS(:,:) = PRS(:,:,3)
IF ( KRR .GE. 4 ) ZRIT(:,:) = PRT(:,:,4)
IF ( KRR .GE. 4 ) ZRIS(:,:) = PRS(:,:,4) 
IF ( KRR .GE. 5 ) ZRSS(:,:) = PRS(:,:,5) 
IF ( KRR .GE. 6 ) ZRGS(:,:) = PRS(:,:,6)
IF ( KRR .GE. 7 ) ZRHS(:,:) = PRS(:,:,7)
!
! Prepare 3D number concentrations
ZCCT(:,:) = 0.
ZCCS(:,:) = 0.
!
IF ( LIMAP%NMOM_C.GE.2 ) ZCCT(:,:) = PSVS(:,:,ISV_LIMA_NC)*PTSTEP
IF ( LIMAP%NMOM_C.GE.2 ) ZCCS(:,:) = PSVS(:,:,ISV_LIMA_NC)
!
IF ( LIMAP%LSCAV .AND. LIMAP%LAERO_MASS ) ZMAS(:,:) = PSVS(:,:,ISV_LIMA_SCAVMASS)
! 
IF ( LIMAP%NMOM_C.GE.1 .AND. LIMAP%NMOD_CCN.GE.1 ) THEN
   ALLOCATE( ZNFS(D%NIJT,D%NKT,LIMAP%NMOD_CCN) )
   ALLOCATE( ZNAS(D%NIJT,D%NKT,LIMAP%NMOD_CCN) )
   ALLOCATE( ZNFT(D%NIJT,D%NKT,LIMAP%NMOD_CCN) )
   ALLOCATE( ZNAT(D%NIJT,D%NKT,LIMAP%NMOD_CCN) )
   ZNFS(:,:,:) = PSVS(:,:,ISV_LIMA_CCN_FREE:ISV_LIMA_CCN_FREE+LIMAP%NMOD_CCN-1)
   ZNAS(:,:,:) = PSVS(:,:,ISV_LIMA_CCN_ACTI:ISV_LIMA_CCN_ACTI+LIMAP%NMOD_CCN-1)
   ZNFT(:,:,:) = PSVS(:,:,ISV_LIMA_CCN_FREE:ISV_LIMA_CCN_FREE+LIMAP%NMOD_CCN-1)*PTSTEP
   ZNAT(:,:,:) = PSVS(:,:,ISV_LIMA_CCN_ACTI:ISV_LIMA_CCN_ACTI+LIMAP%NMOD_CCN-1)*PTSTEP
END IF
!
! Initialize budgets
!
IF ( BUCONF%LBU_ENABLE ) then
  IF ( BUCONF%LBUDGET_TH ) CALL TBUDGETS(NBUDGET_TH)%PTR%INIT_PHY(D, 'CEDS', PTHS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RV ) CALL TBUDGETS(NBUDGET_RV)%PTR%INIT_PHY(D, 'CEDS', ZRVS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RC ) CALL TBUDGETS(NBUDGET_RC)%PTR%INIT_PHY(D, 'CEDS', ZRCS(:,:) * PRHODJ(:,:) )
  !Remark: ZRIS is not modified but source term kept for better coherence with lima_adjust and lima_notadjust
  IF ( BUCONF%LBUDGET_RI ) CALL TBUDGETS(NBUDGET_RI)%PTR%INIT_PHY(D, 'CEDS', ZRIS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_SV ) then
    IF ( LIMAP%NMOM_C.GE.2) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NC)%PTR%INIT_PHY(D, 'CEDS', ZCCS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%LSCAV .AND. LIMAP%LAERO_MASS ) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_SCAVMASS)%PTR%INIT_PHY(D, 'CEDS', ZMAS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%NMOM_C.GE.1 ) then
      DO IL = 1, LIMAP%NMOD_CCN
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_FREE - 1 + IL
        CALL TBUDGETS(IDX)%PTR%INIT_PHY(D, 'CEDS', ZNFS(:,:, IL) * PRHODJ(:,:) )
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_ACTI - 1 + IL
        CALL TBUDGETS(IDX)%PTR%INIT_PHY(D, 'CEDS', ZNAS(:,:, IL) * PRHODJ(:,:) )
      END DO
    END IF
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!
!*       2.     CONDENSATION
!               ------------
!
!
WHERE ( ZRVS(:,:)+ZRCS(:,:)+ZRIS(:,:) < 0.)
  ZRVS(:,:) = -  ZRCS(:,:) - ZRIS(:,:)
END WHERE
!
ZEXNS(:,:) = ( PPABST(:,:) / CST%XP00 ) ** (CST%XRD/CST%XCPD)  
!
ZT(:,:) = ( PTHS(:,:) * PTSTEP ) * ZEXNS(:,:)
ZT2(:,:) = ZT(:,:)
ZCPH(:,:) = CST%XCPD + CST%XCPV  *PTSTEP*   ZRVS(:,:)                     &
     + CST%XCL *PTSTEP* ( ZRCS(:,:) + ZRRS(:,:) )                         &
     + CST%XCI *PTSTEP* ( ZRIS(:,:) + ZRSS(:,:) + ZRGS(:,:) + ZRHS(:,:) )
ZLV(:,:) = CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZT(:,:) -CST%XTT )
ZLS(:,:) = CST%XLSTT + ( CST%XCPV - CST%XCI ) * ( ZT(:,:) -CST%XTT )
!
IF (LIMAP%LADJ) THEN
   ZRV_IN=ZRVS*PTSTEP
   ZRC_IN=ZRCS*PTSTEP
   IF (LIMAP%NMOM_I.EQ.1) THEN
      ZRI_IN=ZRIS*PTSTEP
      GUSERI=.TRUE.
   ELSE
      ZRI_IN=0.
      GUSERI=.FALSE.
   END IF
   IF (OSUBG_COND) THEN
      Z_SIGS=PSIGS
      G_SIGMAS=OSIGMAS
   ELSE
      Z_SIGS=0.
      G_SIGMAS=.TRUE.
   END IF
   !
   CALL CONDENSATION(D, CST, RAIN_ICE_PARAMN, NEBN, TURBN,                   &
        'S', HCONDENS, HLAMBDA3,                                             &
        PPABST, PZZ, PRHODREF, ZT, ZRV_IN, ZRV, ZRC_IN, ZRC, ZRI_IN, ZRI,    &
        ZRRS*PTSTEP,ZRSS*PTSTEP, ZRGS*PTSTEP,                                &
        Z_SIGS, .FALSE., PMFCONV, PCLDFR, Z_SRCS, GUSERI, G_SIGMAS, .FALSE., &
        ZDUM, ZDUM, ZDUM, ZDUM, ZDUM,                                        &
        PSIGQSAT, PLV=ZLV, PLS=ZLS, PCPH=ZCPH )
   !
   IF (LIMAP%NMOM_C.GE.1) THEN
      ZW1(:,:) = (ZRC(:,:) - ZRCS(:,:)*PTSTEP) / PTSTEP
      WHERE( ZW1(:,:) < 0.0 )
         ZW1(:,:) = MAX ( ZW1(:,:), -ZRCS(:,:) )
      ELSEWHERE
         ZW1(:,:) = MIN ( ZW1(:,:),  ZRVS(:,:) )
      END WHERE
      ZRVS(:,:) = ZRVS(:,:) - ZW1(:,:)
      ZRCS(:,:) = ZRCS(:,:) + ZW1(:,:)
      PTHS(:,:) = PTHS(:,:) +        &
           ZW1(:,:) * ZLV(:,:) / (ZCPH(:,:) * PEXNREF(:,:))
   END IF
   !
   IF (LIMAP%NMOM_I.EQ.1) THEN
      PICEFR(:,:)=PCLDFR(:,:)
      ZW2(:,:) = (ZRI(:,:) - ZRIS(:,:)*PTSTEP) / PTSTEP ! idem ZW1 but for Ri
      !
      WHERE( ZW2(:,:) < 0.0 )
         ZW2(:,:) = MAX ( ZW2(:,:), -ZRIS(:,:) )
      ELSEWHERE
         ZW2(:,:) = MIN ( ZW2(:,:),  ZRVS(:,:) )
      END WHERE
      ZRVS(:,:) = ZRVS(:,:) - ZW2(:,:)
      ZRIS(:,:) = ZRIS(:,:) + ZW2(:,:)
      PTHS(:,:) = PTHS(:,:) +        &
           ZW2(:,:) * ZLS(:,:) / (ZCPH(:,:) * PEXNREF(:,:))
   END IF
   !
ELSE
   DO II=1,SIZE(ZRCS,1)
      DO IK=1,SIZE(ZRCS,2)
         IF (ZRCS(II,IK).GE.LIMAP%XRTMIN(2) .AND. ZCCS(II,IK).GE.LIMAP%XCTMIN(2)) THEN
            ZVEC1(II,IK) = MAX( 1.0001, MIN( FLOAT(LIMAW%NAHEN)-0.0001, LIMAW%XAHENINTP1 * ZT(II,IK) + LIMAW%XAHENINTP2 ) )
            IVEC1(II,IK) = INT( ZVEC1(II,IK) )
            ZVEC1(II,IK) = ZVEC1(II,IK) - FLOAT( IVEC1(II,IK) )
            ZW(II,IK) = EXP( CST%XALPW - CST%XBETAW/ZT(II,IK) - CST%XGAMW*ALOG(ZT(II,IK) ) ) ! es_w
            ZW(II,IK) = CST%XMV / CST%XMD * ZW(II,IK) / ( PPABST(II,IK)-ZW(II,IK) ) 
            ZS(II,IK) = ZRVS(II,IK)*PTSTEP / ZW(II,IK) - 1.
            ZW(II,IK) = ZCCS(II,IK)*PTSTEP/(LIMAW%XLBC*ZCCS(II,IK)/ZRCS(II,IK))**LIMAW%XLBEXC
            ZW2(II,IK) = LIMAW%XAHENG3(IVEC1(II,IK)+1)*ZVEC1(II,IK)-LIMAW%XAHENG3(IVEC1(II,IK))*(ZVEC1(II,IK)-1.)
            ZCND(II,IK) = 2.*3.14*1000.*ZW2(II,IK)*ZS(II,IK)*ZW(II,IK)
            IF(ZCND(II,IK).LE.0.) THEN
               ZCND(II,IK) = MAX ( ZCND(II,IK), -ZRCS(II,IK) )
            ELSE
               ZCND(II,IK) = MIN ( ZCND(II,IK),  ZRVS(II,IK) )
            END IF
            ZRVS(II,IK) = ZRVS(II,IK) - ZCND(II,IK)
            ZRCS(II,IK) = ZRCS(II,IK) + ZCND(II,IK)
            PTHS(II,IK) = PTHS(II,IK) + ZCND(II,IK) * ZLV(II,IK) / (ZCPH(II,IK) * PEXNREF(II,IK))
         END IF
      END DO
   END DO
END IF
!
IF (OSUBG_COND) THEN
   PSRCS=Z_SRCS
ENDIF
IF (OSUBG_COND .AND. LIMAP%NMOM_C.GE.2 .AND. LIMAP%LACTI) THEN
   ZW_MF=0.
   ZRV2=ZRVT
   ZRC2=ZRCT
   CALL LIMA_CCN_ACTIVATION (LIMAP, LIMAW, TNSV, D, CST, NEBN,&
        KCARB, KSOA, KSP, ODUST, OSALT, OORILAM,              &
        PRHODREF, PEXNREF, PPABST, ZT2, PDTHRAD, PW_NU+ZW_MF, &
        PAERO,PSOLORG, PMI,  HACTCCN,                         &
        ZTHT, ZRV2, ZRC2, ZCCT, ZRRT, ZNFT, ZNAT,             &
        PCLDFR                                                )      
END IF
!
!-------------------------------------------------------------------------------
!
!*       3.     CLOUD FROM MASS-FLUX SCHEME
!               ---------------------------
!
IF ( .NOT. OSUBG_COND ) THEN
   WHERE (ZRCS(:,:) + ZRIS(:,:) > 1.E-12 / PTSTEP)
      PCLDFR(:,:)  = 1.
   ELSEWHERE
      PCLDFR(:,:)  = 0. 
   ENDWHERE
   IF ( SIZE(PSRCS,2) /= 0 ) THEN
      WHERE (ZRCS(:,:) + ZRIS(:,:) > 1.E-12 / PTSTEP)
         PSRCS(:,:)  = 1.
      ELSEWHERE
         PSRCS(:,:)  = 0.
      ENDWHERE
   END IF
ELSE
   ! We limit PRC_MF+PRI_MF to ZRVS*PTSTEP to avoid negative humidity
   ZW1(:,:)=PRC_MF(:,:)/PTSTEP
   ZW2(:,:)=0.
   IF (LIMAP%NMOM_I.EQ.1) ZW2(:,:)=PRI_MF(:,:)/PTSTEP
   WHERE(ZW1(:,:)+ZW2(:,:)>ZRVS(:,:))
      ZW1(:,:)=ZW1(:,:)*ZRVS(:,:)/(ZW1(:,:)+ZW2(:,:))
      ZW2(:,:)=ZRVS(:,:)-ZW1(:,:)
   ENDWHERE
   ! Compute CF and update rc, ri from MF scheme
   PCLDFR(:,:) = MIN(1.,PCLDFR(:,:)+PCF_MF(:,:))
   ZRVS(:,:)   = ZRVS(:,:) - ZW1(:,:) -ZW2(:,:)
   ZRCS(:,:)   = ZRCS(:,:) + ZW1(:,:)
   IF (LIMAP%NMOM_I.EQ.1) ZRIS(:,:)   = ZRIS(:,:) + ZW2(:,:)
   IF (LIMAP%NMOM_C.GE.1) ZCCS(:,:)   = ZCCT(:,:) / PTSTEP
   IF (LIMAP%NMOD_CCN.GE.1) ZNFS(:,:,:) = ZNFT(:,:,:) / PTSTEP
   IF (LIMAP%NMOD_CCN.GE.1) ZNAS(:,:,:) = ZNAT(:,:,:) / PTSTEP
   PTHS(:,:)   = PTHS(:,:) + &
        (ZW1(:,:) * ZLV(:,:) + ZW2(:,:) * ZLS(:,:)) / ZCPH(:,:)     &
        /  PEXNREF(:,:)
END IF
!
!-------------------------------------------------------------------------------
!
!*       4.     REMOVE SMALL NUMBERS OF DROPLETS
!               --------------------------------
!
IF (LIMAP%NMOM_C .GE. 2) THEN
   ZMASK(:,:) = 0.0
   ZW(:,:) = 0.
   WHERE (ZRCS(:,:) <= ZRTMIN(2) .OR. ZCCS(:,:) <=0.) 
      ZRVS(:,:) = ZRVS(:,:) + ZRCS(:,:) 
      PTHS(:,:) = PTHS(:,:) - ZRCS(:,:)*ZLV(:,:)/(ZCPH(:,:)*ZEXNS(:,:))
      ZRCS(:,:) = 0.0
      ZW(:,:)   = MAX(ZCCS(:,:),0.)
      ZCCS(:,:) = 0.0
      PCLDFR(:,:) = 0.
   END WHERE
   !
   ZW1(:,:) = 0.
   IF (LIMAP%NMOD_CCN.GE.1) ZW1(:,:) = SUM(ZNAS,DIM=3)
   ZW (:,:) = MIN( ZW(:,:), ZW1(:,:) )
   ZW2(:,:) = 0.
   WHERE ( ZW(:,:) > 0. )
      ZMASK(:,:) = 1.0
      ZW2(:,:) = ZW(:,:) / ZW1(:,:)
   ENDWHERE
   !
   IF (LIMAP%NMOD_CCN.GE.1) THEN
      DO IMOD = 1, LIMAP%NMOD_CCN
         ZNFS(:,:,IMOD) = ZNFS(:,:,IMOD) +                           &
              ZMASK(:,:) * ZNAS(:,:,IMOD) * ZW2(:,:)
         ZNAS(:,:,IMOD) = ZNAS(:,:,IMOD) -                           &
              ZMASK(:,:) * ZNAS(:,:,IMOD) * ZW2(:,:)
         ZNAS(:,:,IMOD) = MAX( 0.0 , ZNAS(:,:,IMOD) )
      ENDDO
   END IF
   IF (LIMAP%LSCAV .AND. LIMAP%LAERO_MASS) ZMAS(:,:) = ZMAS(:,:) * (1-ZMASK(:,:))
END IF
!
!-------------------------------------------------------------------------------
!
!*       5.     SAVE CHANGES & CLEANING
!               -----------------------
!
! 3D water mixing ratios
PRS(:,:,1) = ZRVS(:,:)
IF ( KRR .GE. 2 ) PRS(:,:,2) = ZRCS(:,:)
IF ( KRR .GE. 3 ) PRS(:,:,3) = ZRRS(:,:)
IF ( KRR .GE. 4 ) PRS(:,:,4) = ZRIS(:,:)
IF ( KRR .GE. 5 ) PRS(:,:,5) = ZRSS(:,:)
IF ( KRR .GE. 6 ) PRS(:,:,6) = ZRGS(:,:)
IF ( KRR .GE. 7 ) PRS(:,:,7) = ZRHS(:,:)
!
! Prepare 3D number concentrations
!
IF ( LIMAP%NMOM_C.GE.2 ) PSVS(:,:,ISV_LIMA_NC) = ZCCS(:,:)
IF ( LIMAP%LSCAV .AND. LIMAP%LAERO_MASS ) PSVS(:,:,ISV_LIMA_SCAVMASS) = ZMAS(:,:)
IF ( LIMAP%NMOM_C.GE.1 .AND. LIMAP%NMOD_CCN.GE.1 ) THEN
   PSVS(:,:,ISV_LIMA_CCN_FREE:ISV_LIMA_CCN_FREE+LIMAP%NMOD_CCN-1) = ZNFS(:,:,:)
   PSVS(:,:,ISV_LIMA_CCN_ACTI:ISV_LIMA_CCN_ACTI+LIMAP%NMOD_CCN-1) = ZNAS(:,:,:)
END IF
!
! Initialize budgets
!
IF ( BUCONF%LBU_ENABLE ) THEN
  IF ( BUCONF%LBUDGET_TH ) CALL TBUDGETS(NBUDGET_TH)%PTR%END_PHY(D, 'CEDS', PTHS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RV ) CALL TBUDGETS(NBUDGET_RV)%PTR%END_PHY(D, 'CEDS', ZRVS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RC ) CALL TBUDGETS(NBUDGET_RC)%PTR%END_PHY(D, 'CEDS', ZRCS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_RI ) CALL TBUDGETS(NBUDGET_RI)%PTR%END_PHY(D, 'CEDS', ZRIS(:,:) * PRHODJ(:,:) )
  IF ( BUCONF%LBUDGET_SV ) THEN
    IF ( LIMAP%NMOM_C.GE.2) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_NC)%PTR%END_PHY(D, 'CEDS', ZCCS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%LSCAV .AND. LIMAP%LAERO_MASS ) &
      CALL TBUDGETS(NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_SCAVMASS)%PTR%END_PHY(D, 'CEDS', ZMAS(:,:) * PRHODJ(:,:) )
    IF ( LIMAP%NMOM_C.GE.1 ) THEN
      DO IL = 1, LIMAP%NMOD_CCN
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_FREE - 1 + IL
        CALL TBUDGETS(IDX)%PTR%END_PHY(D, 'CEDS', ZNFS(:,:, IL) * PRHODJ(:,:) )
        IDX = NBUDGET_SV1 - 1 + TNSV%NSV_LIMA_CCN_ACTI - 1 + IL
        CALL TBUDGETS(IDX)%PTR%END_PHY(D, 'CEDS', ZNAS(:,:, IL) * PRHODJ(:,:) )
      END DO
    END IF
  END IF
END IF
!
DEALLOCATE(ZRTMIN)
DEALLOCATE(ZCTMIN)
IF (ALLOCATED(ZNFS)) DEALLOCATE(ZNFS)
IF (ALLOCATED(ZNAS)) DEALLOCATE(ZNAS)
IF (ALLOCATED(ZNFT)) DEALLOCATE(ZNFT)
IF (ALLOCATED(ZNAT)) DEALLOCATE(ZNAT)
!
!------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_ADJUST_SPLIT', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_ADJUST_SPLIT
