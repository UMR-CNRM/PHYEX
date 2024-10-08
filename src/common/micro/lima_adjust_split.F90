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
                            PRST,        & ! Aggregate  m.r. at t
                            PRGT,        & ! Graupel    m.r. at t
                            PRHT,        & ! Hail       m.r. at t
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
!                             PIFS,        & ! Free      IFN C. source
!                             PINS,        & ! Nucleated IFN C. source
!                             PNIS           ! Acti. IMM. nuclei C. source
!
!
!
REAL                     :: ZEPS         ! Mv/Md
REAL                     :: ZDT          ! Time increment (2*Delta t or Delta t if cold start)
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
!INTEGER                  :: IRESP      ! Return code of FM routines
INTEGER                  :: JITER,ITERMAX  ! iterative loop for first order adjustment
!INTEGER                  :: ILUOUT     ! Logical unit of output listing 
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
!!$TYPE(TFIELDMETADATA)     :: TZFIELD
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
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
ISV_LIMA_NC       = NSV_LIMA_NC       - NSV_LIMA_BEG + 1
ISV_LIMA_CCN_FREE = NSV_LIMA_CCN_FREE - NSV_LIMA_BEG + 1
ISV_LIMA_CCN_ACTI = NSV_LIMA_CCN_ACTI - NSV_LIMA_BEG + 1
ISV_LIMA_SCAVMASS = NSV_LIMA_SCAVMASS - NSV_LIMA_BEG + 1
ISV_LIMA_NI       = NSV_LIMA_NI       - NSV_LIMA_BEG + 1
ISV_LIMA_IFN_FREE = NSV_LIMA_IFN_FREE - NSV_LIMA_BEG + 1
ISV_LIMA_IFN_NUCL = NSV_LIMA_IFN_NUCL - NSV_LIMA_BEG + 1
ISV_LIMA_IMM_NUCL = NSV_LIMA_IMM_NUCL - NSV_LIMA_BEG + 1
!
!ILUOUT = TLUOUT%NLU
!
ZEPS= CST%XMV / CST%XMD
!
IF (OSUBG_COND) THEN
  ITERMAX=1
ELSE
  ITERMAX=1
END IF
!
ZDT = PTSTEP
!
ISIZE = SIZE(XRTMIN)
ALLOCATE(ZRTMIN(ISIZE))
ZRTMIN(:) = XRTMIN(:) / ZDT
ISIZE = SIZE(XCTMIN)
ALLOCATE(ZCTMIN(ISIZE))
ZCTMIN(:) = XCTMIN(:) / ZDT
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
PRST(:,:,:) = 0.
PRSS(:,:,:) = 0.
PRGT(:,:,:) = 0.
PRGS(:,:,:) = 0.
PRHT(:,:,:) = 0.
PRHS(:,:,:) = 0.
!
IF ( KRR .GE. 2 ) PRCT(:,:,:) = PRS(:,:,:,2)*PTSTEP
IF ( KRR .GE. 2 ) PRCS(:,:,:) = PRS(:,:,:,2)
IF ( KRR .GE. 3 ) PRRT(:,:,:) = PRT(:,:,:,3) 
IF ( KRR .GE. 3 ) PRRS(:,:,:) = PRS(:,:,:,3)
IF ( KRR .GE. 4 ) PRIT(:,:,:) = PRT(:,:,:,4)
IF ( KRR .GE. 4 ) PRIS(:,:,:) = PRS(:,:,:,4) 
IF ( KRR .GE. 5 ) PRST(:,:,:) = PRT(:,:,:,5) 
IF ( KRR .GE. 5 ) PRSS(:,:,:) = PRS(:,:,:,5) 
IF ( KRR .GE. 6 ) PRGT(:,:,:) = PRT(:,:,:,6)
IF ( KRR .GE. 6 ) PRGS(:,:,:) = PRS(:,:,:,6)
IF ( KRR .GE. 7 ) PRHT(:,:,:) = PRT(:,:,:,7)
IF ( KRR .GE. 7 ) PRHS(:,:,:) = PRS(:,:,:,7)
!
! Prepare 3D number concentrations
PCCT(:,:,:) = 0.
PCIT(:,:,:) = 0.
PCCS(:,:,:) = 0.
! PCIS(:,:,:) = 0.
!
IF ( NMOM_C.GE.2 ) PCCT(:,:,:) = PSVS(:,:,:,ISV_LIMA_NC)*PTSTEP
IF ( NMOM_I.GE.2 ) PCIT(:,:,:) = PSVT(:,:,:,ISV_LIMA_NI)
!
IF ( NMOM_C.GE.2 ) PCCS(:,:,:) = PSVS(:,:,:,ISV_LIMA_NC)
! IF (  NMOM_I.GE.2 ) PCIS(:,:,:) = PSVS(:,:,:,ISV_LIMA_NI)
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
! IF ( NMOM_I.GE.1 .AND. NMOD_IFN .GE. 1 ) THEN
!    ALLOCATE( PIFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
!    ALLOCATE( PINS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
!    PIFS(:,:,:,:) = PSVS(:,:,:,ISV_LIMA_IFN_FREE:ISV_LIMA_IFN_FREE+NMOD_IFN-1)
!    PINS(:,:,:,:) = PSVS(:,:,:,ISV_LIMA_IFN_NUCL:ISV_LIMA_IFN_NUCL+NMOD_IFN-1)
! END IF
!
! IF ( NMOD_IMM .GE. 1 ) THEN
!    ALLOCATE( PNIS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IMM) )
!    PNIS(:,:,:,:) = PSVS(:,:,:,ISV_LIMA_IMM_NUCL:ISV_LIMA_IMM_NUCL+NMOD_IMM-1)
! END IF
!
!
if ( nbumod == kmi .and. BUCONF%lbu_enable ) then
  if ( BUCONF%lbudget_th ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_TH), 'CEDS', pths(:, :, :) * prhodj(:, :, :) )
  if ( BUCONF%lbudget_rv ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RV), 'CEDS', prvs(:, :, :) * prhodj(:, :, :) )
  if ( BUCONF%lbudget_rc ) call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_RC), 'CEDS', prcs(:, :, :) * prhodj(:, :, :) )
  !Remark: PRIS is not modified but source term kept for better coherence with lima_adjust and lima_notadjust
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
!     if ( nmom_i.ge.2 ) then
!       call BUDGET_STORE_INIT_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CEDS', pcis(:, :, :) * prhodj(:, :, :) )
!       do jl = 1, nmod_ifn
!         idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_free - 1 + jl
!         call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'CEDS', pifs(:, :, :, jl) * prhodj(:, :, :) )
!         idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl - 1 + jl
!         call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'CEDS', pins(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!       do jl = 1, nmod_imm
!         idx = NBUDGET_SV1 - 1 + nsv_lima_imm_nucl - 1 + jl
!         call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'CEDS', pnis(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!     end if
  end if
end if
!
!-------------------------------------------------------------------------------
!
!
!*       2.     COMPUTE QUANTITIES WITH THE GUESS OF THE FUTURE INSTANT
!               -------------------------------------------------------
!
!*       2.1    remove negative non-precipitating negative water
!               ------------------------------------------------
!
!IF (ANY(PRVS(:,:,:)+PRCS(:,:,:)+PRIS(:,:,:) < 0.) .AND. NVERB>5) THEN
!  WRITE(ILUOUT,*) 'LIMA_ADJUST:  negative values of total water (reset to zero)'
!  WRITE(ILUOUT,*) '  location of minimum PRVS+PRCS+PRIS:',MINLOC(PRVS+PRCS+PRIS)
!  WRITE(ILUOUT,*) '  value of minimum    PRVS+PRCS+PRIS:',MINVAL(PRVS+PRCS+PRIS)
!END IF
!
WHERE ( PRVS(:,:,:)+PRCS(:,:,:)+PRIS(:,:,:) < 0.)
  PRVS(:,:,:) = -  PRCS(:,:,:) - PRIS(:,:,:)
END WHERE
!
!*       2.2    estimate the Exner function at t+1
!
ZEXNS(:,:,:) = ( PPABSTT(:,:,:) / CST%XP00 ) ** (CST%XRD/CST%XCPD)  
!
!    beginning of the iterative loop
!
DO JITER =1,ITERMAX
!
!*       2.3    compute the intermediate temperature at t+1, T*
!  
   ZT(:,:,:) = ( PTHS(:,:,:) * ZDT ) * ZEXNS(:,:,:)
   ZT2(:,:,:) = ZT(:,:,:)
!
!*       2.4    compute the specific heat for moist air (Cph) at t+1
!
   ZCPH(:,:,:) = CST%XCPD + CST%XCPV  *ZDT*   PRVS(:,:,:)                             &
                      + CST%XCL *ZDT* ( PRCS(:,:,:) + PRRS(:,:,:) )             &
                      + CST%XCI *ZDT* ( PRIS(:,:,:) + PRSS(:,:,:) + PRGS(:,:,:) + PRHS(:,:,:) )
!
!*       2.5    compute the latent heat of vaporization Lv(T*) at t+1
!               and of sublimation Ls(T*) at t+1
!
   ZLV(:,:,:) = CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( ZT(:,:,:) -CST%XTT )
   ZLS(:,:,:) = CST%XLSTT + ( CST%XCPV - CST%XCI ) * ( ZT(:,:,:) -CST%XTT )
!
!
!-------------------------------------------------------------------------------
!
!*       3.     FIRST ORDER SUBGRID CONDENSATION SCHEME
!               ---------------------------------------
!
   ZRV_IN=PRVS*PTSTEP
   ZRV2=PRVT
   ZRC_IN=PRCS*PTSTEP
   ZRC2=PRCT
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

   IF (LADJ) THEN
      CALL CONDENSATION(D, CST, RAIN_ICE_PARAMN, NEBN, TURBN,                    &
          'S', HCONDENS, HLAMBDA3,                                             &
           PPABST, PZZ, PRHODREF, ZT, ZRV_IN, ZRV, ZRC_IN, ZRC, ZRI_IN, ZRI,   &
           PRRS*PTSTEP,PRSS*PTSTEP, PRGS*PTSTEP, &
           Z_SIGS, .FALSE., PMFCONV, PCLDFR, Z_SRCS, GUSERI, G_SIGMAS, .FALSE., &
           ZDUM, ZDUM, ZDUM, ZDUM, ZDUM,              &
           ZSIGQSAT2D, PLV=ZLV, PLS=ZLS, PCPH=ZCPH )
   END IF
   IF (OSUBG_COND .AND. NMOM_C.GE.2 .AND. LACTI) THEN
      PSRCS=Z_SRCS
      ZW_MF=0.
      CALL LIMA_CCN_ACTIVATION (CST,                          &
           PRHODREF, PEXNREF, PPABST, ZT2, PDTHRAD, PW_NU+ZW_MF, &
           PTHT, ZRV2, ZRC2, PCCT, PRRT, PNFT, PNAT,             &
           PCLDFR                                                )      
   END IF

END DO
!
!*       5.1    compute the sources
!
IF (LADJ) THEN
                                                           !         Rc - Rc*
   ZW1(:,:,:) = (ZRC(:,:,:) - PRCS(:,:,:)*PTSTEP) / PTSTEP ! Pcon = ----------
                                                           !         2 Delta t
   WHERE( ZW1(:,:,:) < 0.0 )
      ZW1(:,:,:) = MAX ( ZW1(:,:,:), -PRCS(:,:,:) )
   ELSEWHERE
      ZW1(:,:,:) = MIN ( ZW1(:,:,:),  PRVS(:,:,:) )
   END WHERE
   PRVS(:,:,:) = PRVS(:,:,:) - ZW1(:,:,:)
   PRCS(:,:,:) = PRCS(:,:,:) + ZW1(:,:,:)
   PTHS(:,:,:) = PTHS(:,:,:) +        &
              ZW1(:,:,:) * ZLV(:,:,:) / (ZCPH(:,:,:) * PEXNREF(:,:,:))
ELSE
   DO JI=1,SIZE(PRCS,1)
      DO JJ=1,SIZE(PRCS,2)
         DO JK=1,SIZE(PRCS,3)
            IF (PRCS(JI,JJ,JK).GE.XRTMIN(2) .AND. PCCS(JI,JJ,JK).GE.XCTMIN(2)) THEN
               ZVEC1(JI,JJ,JK) = MAX( 1.0001, MIN( FLOAT(NAHEN)-0.0001, XAHENINTP1 * ZT(JI,JJ,JK) + XAHENINTP2 ) )
               IVEC1(JI,JJ,JK) = INT( ZVEC1(JI,JJ,JK) )
               ZVEC1(JI,JJ,JK) = ZVEC1(JI,JJ,JK) - FLOAT( IVEC1(JI,JJ,JK) )
               ZW(JI,JJ,JK)=EXP( CST%XALPW - CST%XBETAW/ZT(JI,JJ,JK) - CST%XGAMW*ALOG(ZT(JI,JJ,JK) ) ) ! es_w
               ZW(JI,JJ,JK)=ZEPS*ZW(JI,JJ,JK) / ( PPABST(JI,JJ,JK)-ZW(JI,JJ,JK) ) 
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
IF (NMOM_I.EQ.1) THEN
   ZW2(:,:,:) = (ZRI(:,:,:) - PRIS(:,:,:)*PTSTEP) / PTSTEP ! idem ZW1 but for Ri
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
!
!*       5.2    compute the cloud fraction PCLDFR
!
IF ( .NOT. OSUBG_COND ) THEN
  WHERE (PRCS(:,:,:) + PRIS(:,:,:) > 1.E-12 / PTSTEP)
    PCLDFR(:,:,:)  = 1.
  ELSEWHERE
    PCLDFR(:,:,:)  = 0. 
  ENDWHERE 
  IF ( SIZE(PSRCS,3) /= 0 ) THEN
     WHERE (PRCS(:,:,:) + PRIS(:,:,:) > 1.E-12 / ZDT)
        PSRCS(:,:,:)  = 1.
     ELSEWHERE
        PSRCS(:,:,:)  = 0.
     ENDWHERE
  END IF
ELSE
! We limit PRC_MF+PRI_MF to PRVS*PTSTEP to avoid negative humidity
   ZW1(:,:,:)=PRC_MF(:,:,:)/PTSTEP
   IF (NMOM_I.EQ.1) THEN
      ZW2(:,:,:)=PRI_MF(:,:,:)/PTSTEP
   ELSE
      ZW2(:,:,:)=0.
   END IF
   WHERE(ZW1(:,:,:)+ZW2(:,:,:)>PRVS(:,:,:))
      ZW1(:,:,:)=ZW1(:,:,:)*PRVS(:,:,:)/(ZW1(:,:,:)+ZW2(:,:,:))
      ZW2(:,:,:)=PRVS(:,:,:)-ZW1(:,:,:)
   ENDWHERE
! Compute CF and update rc, ri from MF scheme
   PRVS(:,:,:)   = PRVS(:,:,:) - ZW1(:,:,:) -ZW2(:,:,:)
   PRCS(:,:,:)   = PRCS(:,:,:) + ZW1(:,:,:)
   PRIS(:,:,:)   = PRIS(:,:,:) + ZW2(:,:,:)
   PCCS(:,:,:)   = PCCT(:,:,:) / PTSTEP
   PNFS(:,:,:,:) = PNFT(:,:,:,:) / PTSTEP
   PNAS(:,:,:,:) = PNAT(:,:,:,:) / PTSTEP
   PTHS(:,:,:)   = PTHS(:,:,:) + &
                   (ZW1(:,:,:) * ZLV(:,:,:) + ZW2 * ZLS(:,:,:)) / ZCPH(:,:,:)     &
                   /  PEXNREF(:,:,:)
END IF
!
! Remove cloud droplets if there are few
!
ZMASK(:,:,:) = 0.0
ZW(:,:,:) = 0.
IF (NMOM_C .GE. 2) THEN
   WHERE (PRCS(:,:,:) <= ZRTMIN(2) .OR. PCCS(:,:,:) <=0.) 
      PRVS(:,:,:) = PRVS(:,:,:) + PRCS(:,:,:) 
      PTHS(:,:,:) = PTHS(:,:,:) - PRCS(:,:,:)*ZLV(:,:,:)/(ZCPH(:,:,:)*ZEXNS(:,:,:))
      PRCS(:,:,:) = 0.0
      ZW(:,:,:)   = MAX(PCCS(:,:,:),0.)
      PCCS(:,:,:) = 0.0
   END WHERE
END IF
!
ZW1(:,:,:) = 0.
IF (NMOM_C.GE.1 .AND. NMOD_CCN.GE.1) ZW1(:,:,:) = SUM(PNAS,DIM=4)
ZW (:,:,:) = MIN( ZW(:,:,:), ZW1(:,:,:) )
ZW2(:,:,:) = 0.
WHERE ( ZW(:,:,:) > 0. )
   ZMASK(:,:,:) = 1.0
   ZW2(:,:,:) = ZW(:,:,:) / ZW1(:,:,:)
ENDWHERE
!
IF (NMOM_C.GE.1 .AND. NMOD_CCN.GE.1) THEN
   DO JMOD = 1, NMOD_CCN
      PNFS(:,:,:,JMOD) = PNFS(:,:,:,JMOD) +                           &
           ZMASK(:,:,:) * PNAS(:,:,:,JMOD) * ZW2(:,:,:)
      PNAS(:,:,:,JMOD) = PNAS(:,:,:,JMOD) -                           &
           ZMASK(:,:,:) * PNAS(:,:,:,JMOD) * ZW2(:,:,:)
      PNAS(:,:,:,JMOD) = MAX( 0.0 , PNAS(:,:,:,JMOD) )
   ENDDO
END IF
!
IF (LSCAV .AND. LAERO_MASS) PMAS(:,:,:) = PMAS(:,:,:) * (1-ZMASK(:,:,:))
!
!
!
PICEFR(:,:,:)=0.
IF (NMOM_I.EQ.1) THEN
   WHERE(PICEFR(:,:,:)<1.E-10 .AND. PRIT(:,:,:)>XRTMIN(4) .AND. PCIT(:,:,:)>XCTMIN(4)) PICEFR(:,:,:)=PCLDFR(:,:,:)
ELSE
   WHERE(PICEFR(:,:,:)<1.E-10 .AND. PRIT(:,:,:)>XRTMIN(4) .AND. PCIT(:,:,:)>XCTMIN(4)) PICEFR(:,:,:)=1.
END IF
!
!!$IF ( tpfile%lopened ) THEN
!!$   TZFIELD = TFIELDMETADATA(   &
!!$     CMNHNAME   = 'NEB',       &
!!$     CSTDNAME   = '',          &
!!$     CLONGNAME  = 'NEB',       &
!!$     CUNITS     = '1',         &
!!$     CDIR       = 'XY',        &
!!$     CCOMMENT   = 'X_Y_Z_NEB', &
!!$     NGRID      = 1,           &
!!$     NTYPE      = TYPEREAL,    &
!!$     NDIMS      = 3,           &
!!$     LTIMEDEP   = .TRUE.       )
!!$   CALL IO_Field_write(TPFILE,TZFIELD,PCLDFR)
!!$END IF
!
!
!*       6.  SAVE CHANGES IN PRS AND PSVS
!            ----------------------------
!
! Prepare 3D water mixing ratios
PRS(:,:,:,1) = PRVS(:,:,:)
IF ( KRR .GE. 2 ) PRS(:,:,:,2) = PRCS(:,:,:)
IF ( KRR .GE. 3 ) PRS(:,:,:,3) = PRRS(:,:,:)
IF ( KRR .GE. 4 ) PRS(:,:,:,4) = PRIS(:,:,:)
IF ( KRR .GE. 5 ) PRS(:,:,:,5) = PRSS(:,:,:)
IF ( KRR .GE. 6 ) PRS(:,:,:,6) = PRGS(:,:,:)
IF ( KRR .GE. 7 ) PRS(:,:,:,7) = PRHS(:,:,:)
!
! Prepare 3D number concentrations
!
IF ( NMOM_C.GE.2 ) PSVS(:,:,:,ISV_LIMA_NC) = PCCS(:,:,:)
! IF ( NMOM_I.GE.2 ) PSVS(:,:,:,ISV_LIMA_NI) = PCIS(:,:,:)
!
IF ( LSCAV .AND. LAERO_MASS ) PSVS(:,:,:,ISV_LIMA_SCAVMASS) = PMAS(:,:,:)
! 
IF ( NMOM_C.GE.1 .AND. NMOD_CCN.GE.1 ) THEN
   PSVS(:,:,:,ISV_LIMA_CCN_FREE:ISV_LIMA_CCN_FREE+NMOD_CCN-1) = PNFS(:,:,:,:)
   PSVS(:,:,:,ISV_LIMA_CCN_ACTI:ISV_LIMA_CCN_ACTI+NMOD_CCN-1) = PNAS(:,:,:,:)
END IF
!
! IF ( NMOM_I.GE.1 .AND. NMOD_IFN .GE. 1 ) THEN
!    PSVS(:,:,:,ISV_LIMA_IFN_FREE:ISV_LIMA_IFN_FREE+NMOD_IFN-1) = PIFS(:,:,:,:)
!    PSVS(:,:,:,ISV_LIMA_IFN_NUCL:ISV_LIMA_IFN_NUCL+NMOD_IFN-1) = PINS(:,:,:,:)
! END IF
!
! IF ( NMOM_I.GE.1 .AND. NMOD_IMM .GE. 1 ) THEN
!    PSVS(:,:,:,ISV_LIMA_IMM_NUCL:ISV_LIMA_IMM_NUCL+NMOD_IMM-1) = PNIS(:,:,:,:)
! END IF
!
! write SSI in LFI
!
!!$IF ( tpfile%lopened ) THEN
!!$   ZT(:,:,:) = ( PTHS(:,:,:) * ZDT ) * ZEXNS(:,:,:)
!!$   ZW(:,:,:) = EXP( XALPI - XBETAI/ZT(:,:,:) - XGAMI*ALOG(ZT(:,:,:) ) )
!!$   ZW1(:,:,:)= PPABSTT(:,:,:)
!!$   ZW(:,:,:) = PRVT(:,:,:)*( ZW1(:,:,:)-ZW(:,:,:) ) / ( (XMV/XMD) * ZW(:,:,:) ) - 1.0
!!$   
!!$   TZFIELD = TFIELDMETADATA(   &
!!$     CMNHNAME   = 'SSI',       &
!!$     CSTDNAME   = '',          &
!!$     CLONGNAME  = 'SSI',       &
!!$     CUNITS     = '',          &
!!$     CDIR       = 'XY',        &
!!$     CCOMMENT   = 'X_Y_Z_SSI', &
!!$     NGRID      = 1,           &
!!$     NTYPE      = TYPEREAL,    &
!!$     NDIMS      = 3,           &
!!$     LTIMEDEP   = .TRUE.       )
!!$   CALL IO_Field_write(TPFILE,TZFIELD,ZW)
!!$END IF
!
!
!*       7.  STORE THE BUDGET TERMS
!            ----------------------
!
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
!     if ( nmom_i.ge.2 ) then
!       call BUDGET_STORE_END_PHY(D, TBUDGETS(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CEDS', pcis(:, :, :) * prhodj(:, :, :) )
!       do jl = 1, nmod_ifn
!         idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_free - 1 + jl
!         call BUDGET_STORE_END_PHY(D, TBUDGETS(idx), 'CEDS', pifs(:, :, :, jl) * prhodj(:, :, :) )
!         idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl - 1 + jl
!         call BUDGET_STORE_END_PHY(D, TBUDGETS(idx), 'CEDS', pins(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!       do jl = 1, nmod_imm
!         idx = NBUDGET_SV1 - 1 + nsv_lima_imm_nucl - 1 + jl
!         call BUDGET_STORE_INIT_PHY(D, TBUDGETS(idx), 'CEDS', pnis(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!     end if
  end if
end if
!++cb++
DEALLOCATE(ZRTMIN)
DEALLOCATE(ZCTMIN)
IF (ALLOCATED(PNFS)) DEALLOCATE(PNFS)
IF (ALLOCATED(PNAS)) DEALLOCATE(PNAS)
IF (ALLOCATED(PNFT)) DEALLOCATE(PNFT)
IF (ALLOCATED(PNAT)) DEALLOCATE(PNAT)
! IF (ALLOCATED(PIFS)) DEALLOCATE(PIFS)
! IF (ALLOCATED(PINS)) DEALLOCATE(PINS)
! IF (ALLOCATED(PNIS)) DEALLOCATE(PNIS)
!--cb--
!
!------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_ADJUST_SPLIT
