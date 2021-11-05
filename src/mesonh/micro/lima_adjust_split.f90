!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################
      MODULE MODI_LIMA_ADJUST_SPLIT
!     #############################
!
INTERFACE
!
      SUBROUTINE LIMA_ADJUST_SPLIT(KRR, KMI, TPFILE, HCONDENS, HLAMBDA3,        &
                             OSUBG_COND, OSIGMAS, PTSTEP, PSIGQSAT,             &
                             PRHODREF, PRHODJ, PEXNREF, PPABSM, PSIGS, PMFCONV, &
                             PPABST, PZZ, PDTHRAD, PW_NU,                       &
                             PRT, PRS, PSVT, PSVS,                              &
                             PTHS, PSRCS, PCLDFR, PRC_MF, PCF_MF                )
!
USE MODD_IO,    ONLY: TFILEDATA
USE MODD_NSV,   only: NSV_LIMA_BEG
!
INTEGER,                  INTENT(IN)   :: KRR        ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KMI        ! Model index 
TYPE(TFILEDATA),          INTENT(IN)   :: TPFILE     ! Output file
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
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODREF  ! Dry density of the 
                                                     ! reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ    ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF   ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABSM    ! Absolute Pressure at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PSIGS     ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PMFCONV   ! 
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST    ! Absolute Pressure at t     
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PZZ       !     
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PDTHRAD   ! Radiative temperature tendency
REAL, DIMENSION(:,:,:),   INTENT(IN)   :: PW_NU     ! updraft velocity used for
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT       ! m.r. at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS       ! m.r. source
!
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(IN)    :: PSVT ! Concentrations at time t
!
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(INOUT) :: PSVS ! Concentration sources
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS      ! Theta source
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS     ! Second-order flux
                                                     ! s'rc'/2Sigma_s2 at time t+1
                                                     ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   :: PCLDFR    ! Cloud fraction          
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
!
END SUBROUTINE LIMA_ADJUST_SPLIT
!
END INTERFACE
!
END MODULE MODI_LIMA_ADJUST_SPLIT
!
!     ###########################################################################
      SUBROUTINE LIMA_ADJUST_SPLIT(KRR, KMI, TPFILE, HCONDENS, HLAMBDA3,        &
                             OSUBG_COND, OSIGMAS, PTSTEP, PSIGQSAT,             &
                             PRHODREF, PRHODJ, PEXNREF, PPABSM, PSIGS, PMFCONV, &
                             PPABST, PZZ, PDTHRAD, PW_NU,                       &
                             PRT, PRS, PSVT, PSVS,                              &
                             PTHS, PSRCS, PCLDFR, PRC_MF, PCF_MF                )
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,           only: lbu_enable, nbumod,                                          &
                                 lbudget_th, lbudget_rv, lbudget_rc, lbudget_ri, lbudget_sv,  &
                                 NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RI, NBUDGET_SV1, &
                                 tbudgets
USE MODD_CONF
USE MODD_CST
use modd_field,            only: TFIELDDATA, TYPEREAL
USE MODD_IO,               ONLY: TFILEDATA
USE MODD_LUNIT_n,          ONLY: TLUOUT
USE MODD_NSV
USE MODD_PARAMETERS
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_COLD
USE MODD_PARAM_LIMA_MIXED
USE MODD_PARAM_LIMA_WARM
!
use mode_budget,           only: Budget_store_init, Budget_store_end
USE MODE_IO_FIELD_WRITE,   only: IO_Field_write
use mode_msg
use mode_tools,            only: Countjv
!
USE MODI_CONDENS
USE MODI_CONDENSATION
USE MODI_LIMA_FUNCTIONS
USE MODI_LIMA_CCN_ACTIVATION
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN)   :: KRR        ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KMI        ! Model index 
TYPE(TFILEDATA),          INTENT(IN)   :: TPFILE     ! Output file
CHARACTER(len=80),        INTENT(IN)    :: HCONDENS
CHARACTER(len=4),         INTENT(IN)    :: HLAMBDA3  ! formulation for lambda3 coeff
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid 
                                                     ! Condensation
LOGICAL,                  INTENT(IN)   :: OSIGMAS    ! Switch for Sigma_s: 
                                                     ! use values computed in CONDENSATION
                                                     ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PTSTEP     ! Time step          
REAL,                     INTENT(IN)   :: PSIGQSAT   ! coeff applied to qsat variance contribution
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODREF  ! Dry density of the 
                                                     ! reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ    ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF   ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABSM    ! Absolute Pressure at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PSIGS     ! Sigma_s at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PMFCONV   ! 
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST    ! Absolute Pressure at t     
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PZZ       !     
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDTHRAD   ! Radiative temperature tendency
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU     ! updraft velocity used for
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT       ! m.r. at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS       ! m.r. source
!
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(IN)    :: PSVT ! Concentrations at time t
!
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(INOUT) :: PSVS ! Concentration sources
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS      ! Theta source
!
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS     ! Second-order flux
                                                     ! s'rc'/2Sigma_s2 at time t+1
                                                     ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   :: PCLDFR    ! Cloud fraction          
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PRC_MF! Convective Mass Flux liquid mixing ratio
REAL, DIMENSION(:,:,:),     INTENT(IN)    :: PCF_MF! Convective Mass Flux Cloud fraction 
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
!
                            PRVS,        & ! Water vapor m.r. source
                            PRCS,        & ! Cloud water m.r. source
                            PRRS,        & ! Rain water m.r. source
                            PRIS,        & ! Cloud ice  m.r. source
                            PRSS,        & ! Aggregate  m.r. source
                            PRGS,        & ! Graupel    m.r. source
!
                            PCCT,        & ! Cloud water conc. at t
                            PCIT,        & ! Cloud ice   conc. at t
!
                            PCCS,        & ! Cloud water C. source
                            PMAS,        & ! Mass of scavenged AP
                            PCIS           ! Ice crystal C. source
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
                            ZRV, ZRV2,  &
                            ZRC, ZRC2,  &
                            ZRI,  &
                            ZSIGS, &
                            ZW_MF
LOGICAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) &
                         :: GMICRO ! Test where to compute cond/dep proc.
INTEGER                  :: IMICRO
REAL, DIMENSION(:), ALLOCATABLE &
                         :: ZRVT, ZRCT, ZRIT, ZRVS, ZRCS, ZRIS, ZTHS,        &
                            ZCCT, ZCIT, ZCCS, ZCIS,                          &
                            ZRHODREF, ZZT, ZPRES, ZEXNREF, ZZCPH,            &
                            ZZW, ZLVFACT, ZLSFACT,                           &
                            ZRVSATW, ZRVSATI, ZRVSATW_PRIME, ZRVSATI_PRIME,  &
                            ZAW, ZAI, ZCJ, ZKA, ZDV, ZITW, ZITI, ZAWW, ZAIW, &
                            ZAWI, ZAII, ZFACT, ZDELTW,                       &
                            ZDELTI, ZDELT1, ZDELT2, ZCND, ZDEP, ZS, ZVEC1, ZZW2
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1
!
INTEGER                  :: IRESP      ! Return code of FM routines
INTEGER                  :: IIU,IJU,IKU! dimensions of dummy arrays
INTEGER                  :: IKB        ! K index value of the first inner mass point
INTEGER                  :: IKE        ! K index value of the last inner mass point
INTEGER                  :: IIB,IJB    ! Horz index values of the first inner mass points
INTEGER                  :: IIE,IJE    ! Horz index values of the last inner mass points
INTEGER                  :: JITER,ITERMAX  ! iterative loop for first order adjustment
INTEGER                  :: ILUOUT     ! Logical unit of output listing 
!
INTEGER                           :: ISIZE
REAL, DIMENSION(:), ALLOCATABLE   :: ZRTMIN
REAL, DIMENSION(:), ALLOCATABLE   :: ZCTMIN
!
integer :: idx
INTEGER , DIMENSION(SIZE(GMICRO)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
INTEGER                           :: JMOD, JMOD_IFN, JMOD_IMM
!
INTEGER , DIMENSION(3) :: BV
TYPE(TFIELDDATA)  :: TZFIELD
!
!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
ILUOUT = TLUOUT%NLU
!
IIU = SIZE(PEXNREF,1)
IJU = SIZE(PEXNREF,2)
IKU = SIZE(PEXNREF,3)
IIB = 1 + JPHEXT
IIE = SIZE(PRHODJ,1) - JPHEXT
IJB = 1 + JPHEXT
IJE = SIZE(PRHODJ,2) - JPHEXT
IKB = 1 + JPVEXT
IKE = SIZE(PRHODJ,3) - JPVEXT
!
ZEPS= XMV / XMD
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
!
! Prepare 3D number concentrations
PCCT(:,:,:) = 0.
PCIT(:,:,:) = 0.
PCCS(:,:,:) = 0.
! PCIS(:,:,:) = 0.
!
IF ( LWARM ) PCCT(:,:,:) = PSVS(:,:,:,NSV_LIMA_NC)*PTSTEP
IF ( LCOLD ) PCIT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NI)
!
IF ( LWARM ) PCCS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NC)
! IF ( LCOLD ) PCIS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NI)
!
IF ( LSCAV .AND. LAERO_MASS ) PMAS(:,:,:) = PSVS(:,:,:,NSV_LIMA_SCAVMASS)
! 
IF ( LWARM .AND. NMOD_CCN.GE.1 ) THEN
   ALLOCATE( PNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ALLOCATE( PNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ALLOCATE( PNFT(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ALLOCATE( PNAT(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   PNFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1)
   PNAS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1)
   PNFT(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1)*PTSTEP
   PNAT(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1)*PTSTEP
END IF
!
! IF ( LCOLD .AND. NMOD_IFN .GE. 1 ) THEN
!    ALLOCATE( PIFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
!    ALLOCATE( PINS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
!    PIFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1)
!    PINS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1)
! END IF
!
! IF ( NMOD_IMM .GE. 1 ) THEN
!    ALLOCATE( PNIS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IMM) )
!    PNIS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1)
! END IF
!
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'CEDS', pths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'CEDS', prvs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'CEDS', prcs(:, :, :) * prhodj(:, :, :) )
  !Remark: PRIS is not modified but source term kept for better coherence with lima_adjust and lima_notadjust
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'CEDS', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if ( lwarm ) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'CEDS', pccs(:, :, :) * prhodj(:, :, :) )
    if ( lscav .and. laero_mass ) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_scavmass), 'CEDS', pmas(:, :, :) * prhodj(:, :, :) )
    if ( lwarm ) then
      do jl = 1, nmod_ccn
        idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_free - 1 + jl
        call Budget_store_init( tbudgets(idx), 'CEDS', pnfs(:, :, :, jl) * prhodj(:, :, :) )
        idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_acti - 1 + jl
        call Budget_store_init( tbudgets(idx), 'CEDS', pnas(:, :, :, jl) * prhodj(:, :, :) )
      end do
    end if
!     if ( lcold ) then
!       call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CEDS', pcis(:, :, :) * prhodj(:, :, :) )
!       do jl = 1, nmod_ifn
!         idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_free - 1 + jl
!         call Budget_store_init( tbudgets(idx), 'CEDS', pifs(:, :, :, jl) * prhodj(:, :, :) )
!         idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl - 1 + jl
!         call Budget_store_init( tbudgets(idx), 'CEDS', pins(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!       do jl = 1, nmod_imm
!         idx = NBUDGET_SV1 - 1 + nsv_lima_imm_nucl - 1 + jl
!         call Budget_store_init( tbudgets(idx), 'CEDS', pnis(:, :, :, jl) * prhodj(:, :, :) )
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
IF (ANY(PRVS(:,:,:)+PRCS(:,:,:)+PRIS(:,:,:) < 0.) .AND. NVERB>5) THEN
  WRITE(ILUOUT,*) 'LIMA_ADJUST:  negative values of total water (reset to zero)'
  WRITE(ILUOUT,*) '  location of minimum PRVS+PRCS+PRIS:',MINLOC(PRVS+PRCS+PRIS)
  WRITE(ILUOUT,*) '  value of minimum    PRVS+PRCS+PRIS:',MINVAL(PRVS+PRCS+PRIS)
END IF
!
WHERE ( PRVS(:,:,:)+PRCS(:,:,:)+PRIS(:,:,:) < 0.)
  PRVS(:,:,:) = -  PRCS(:,:,:) - PRIS(:,:,:)
END WHERE
!
!*       2.2    estimate the Exner function at t+1
!
ZEXNS(:,:,:) = ( (2. * PPABST(:,:,:) - PPABSM(:,:,:)) / XP00 ) ** (XRD/XCPD)  
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
   ZCPH(:,:,:) = XCPD + XCPV  *ZDT*   PRVS(:,:,:)                             &
                      + XCL   *ZDT* ( PRCS(:,:,:) + PRRS(:,:,:) )             &
                      + XCI   *ZDT* ( PRIS(:,:,:) + PRSS(:,:,:) + PRGS(:,:,:) )
!
!*       2.5    compute the latent heat of vaporization Lv(T*) at t+1
!               and of sublimation Ls(T*) at t+1
!
   ZLV(:,:,:) = XLVTT + ( XCPV - XCL ) * ( ZT(:,:,:) -XTT )
   ZLS(:,:,:) = XLSTT + ( XCPV - XCI ) * ( ZT(:,:,:) -XTT )
!
!
!-------------------------------------------------------------------------------
!
!*       3.     FIRST ORDER SUBGRID CONDENSATION SCHEME
!               ---------------------------------------
!
   IF ( OSUBG_COND ) THEN
     !
      ZRV=PRVS*PTSTEP
      ZRC=PRCS*PTSTEP
      ZRV2=PRVT
      ZRC2=PRCT
      ZRI=0.
      ZSIGS=PSIGS
      CALL CONDENSATION(IIU, IJU, IKU, IIB, IIE, IJB, IJE, IKB, IKE, 1, 'S',   &
           HCONDENS, HLAMBDA3, &
           PPABST, PZZ, PRHODREF, ZT, ZRV, ZRC, ZRI, PRSS*PTSTEP, PRGS*PTSTEP, &
           ZSIGS, PMFCONV, PCLDFR, PSRCS, .FALSE., OSIGMAS, &
           PSIGQSAT, PLV=ZLV, PLS=ZLS, PCPH=ZCPH )
      PCLDFR(:,:,:) = MIN(PCLDFR(:,:,:) + PCF_MF(:,:,:) , 1.)
      ZRV(:,:,:) = ZRV(:,:,:) - MAX(MIN(PRC_MF(:,:,:), ZRV(:,:,:)),0.)
      ZRC(:,:,:) = ZRC(:,:,:) + MAX(MIN(PRC_MF(:,:,:), ZRV(:,:,:)),0.)
      ZW_MF=0.
      CALL LIMA_CCN_ACTIVATION (TPFILE,                         &
           PRHODREF, PEXNREF, PPABST, ZT2, PDTHRAD, PW_NU+ZW_MF, &
           PTHT, ZRV2, ZRC2, PCCT, PRRT, PNFT, PNAT,              &
           PCLDFR                                               )
!
   ELSE
!
!-------------------------------------------------------------------------------
!
!
!
!*              FULLY IMPLICIT CONDENSATION SCHEME
!               ---------------------------------
! 
!*              select cases where r_c>0
! 
!
      GMICRO(:,:,:) = .FALSE.
      GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =( PRCS(IIB:IIE,IJB:IJE,IKB:IKE)>0. .AND.        &
                                         PCCS(IIB:IIE,IJB:IJE,IKB:IKE)>0.      )
      IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
      IF( IMICRO >= 1 ) THEN
         ALLOCATE(ZRVT(IMICRO))
         ALLOCATE(ZRCT(IMICRO))
!
         ALLOCATE(ZRVS(IMICRO))
         ALLOCATE(ZRCS(IMICRO))
         ALLOCATE(ZCCS(IMICRO))
         ALLOCATE(ZTHS(IMICRO))
!
         ALLOCATE(ZRHODREF(IMICRO))
         ALLOCATE(ZZT(IMICRO))
         ALLOCATE(ZPRES(IMICRO))
         ALLOCATE(ZEXNREF(IMICRO))
         ALLOCATE(ZZCPH(IMICRO))
         DO JL=1,IMICRO
            ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
            ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
            !
            ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
            ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
            ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
            ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
            !
            ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
            ZZT(JL) = ZT(I1(JL),I2(JL),I3(JL))
            ZPRES(JL) = 2.0*PPABST(I1(JL),I2(JL),I3(JL))-PPABSM(I1(JL),I2(JL),I3(JL))
            ZEXNREF(JL) = PEXNREF(I1(JL),I2(JL),I3(JL))
            ZZCPH(JL) = ZCPH(I1(JL),I2(JL),I3(JL))
         ENDDO
         ALLOCATE(ZZW(IMICRO))
         ALLOCATE(ZLVFACT(IMICRO))
         ALLOCATE(ZRVSATW(IMICRO))
         ALLOCATE(ZCND(IMICRO))
         ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZCPH(:) ! L_v/C_ph
         ZZW(:) = EXP( XALPW - XBETAW/ZZT(:) - XGAMW*ALOG(ZZT(:) ) ) ! es_w
         ZRVSATW(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_sw

         IF (LADJ) THEN
            ALLOCATE(ZRVSATW_PRIME(IMICRO))
            ALLOCATE(ZAWW(IMICRO))
            ALLOCATE(ZDELT1(IMICRO))
            ALLOCATE(ZDELT2(IMICRO))
            ZRVSATW_PRIME(:) = (( XBETAW/ZZT(:) - XGAMW ) / ZZT(:))  &  ! r'_sw
                               * ZRVSATW(:) * ( 1. + ZRVSATW(:)/ZEPS )
            ZAWW(:) = 1.0 + ZRVSATW_PRIME(:)*ZLVFACT(:)
            ZDELT2(:) = (ZRVSATW_PRIME(:)*ZLVFACT(:)/ZAWW(:)) *                     &
                        ( ((-2.*XBETAW+XGAMW*ZZT(:))/(XBETAW-XGAMW*ZZT(:))          &
                        + (XBETAW/ZZT(:)-XGAMW)*(1.0+2.0*ZRVSATW(:)/ZEPS))/ZZT(:) )
            ZDELT1(:) = (ZLVFACT(:)/ZAWW(:)) * ( ZRVSATW(:) - ZRVS(:)*ZDT )
            ZCND(:) = - ZDELT1(:)*( 1.0 + 0.5*ZDELT1(:)*ZDELT2(:) ) / (ZLVFACT(:)*ZDT)
            DEALLOCATE(ZRVSATW_PRIME)
            DEALLOCATE(ZAWW)
            DEALLOCATE(ZDELT1)
            DEALLOCATE(ZDELT2)
         ELSE
            ALLOCATE(ZS(IMICRO))
            ALLOCATE(ZZW2(IMICRO))
            ALLOCATE(ZVEC1(IMICRO))
            ALLOCATE(IVEC1(IMICRO))
            ZVEC1(:) = MAX( 1.0001, MIN( FLOAT(NAHEN)-0.0001, XAHENINTP1 * ZZT(:) + XAHENINTP2 ) )
            IVEC1(:) = INT( ZVEC1(:) )
            ZVEC1(:) = ZVEC1(:) - FLOAT( IVEC1(:) )
            ZS(:) = ZRVS(:)*PTSTEP / ZRVSATW(:) - 1.
            ZZW(:) = ZCCS(:)*PTSTEP/(XLBC*ZCCS(:)/ZRCS(:))**XLBEXC
            ZZW2(:) = XAHENG3(IVEC1(:)+1)*ZVEC1(:)-XAHENG3(IVEC1(:))*(ZVEC1(:)-1.)
            ZCND(:) = 2.*3.14*1000.*ZZW2(:)*ZS(:)*ZZW(:)
            DEALLOCATE(ZS)
            DEALLOCATE(ZZW2)
            DEALLOCATE(ZVEC1)
            DEALLOCATE(IVEC1)
         END IF
!
!
! Integration
!
         WHERE( ZCND(:) < 0.0 )
            ZCND(:) = MAX ( ZCND(:), -ZRCS(:) )
         ELSEWHERE
            ZCND(:) = MIN ( ZCND(:),  ZRVS(:) )
         END WHERE
         ZRVS(:) = ZRVS(:) - ZCND(:)
         ZRCS(:) = ZRCS(:) + ZCND(:)
         ZTHS(:) = ZTHS(:) + ZCND(:) * ZLVFACT(:) / ZEXNREF(:)
!
         ZW(:,:,:) = PRVS(:,:,:)
         PRVS(:,:,:) = UNPACK( ZRVS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
         ZW(:,:,:) = PRCS(:,:,:)
         PRCS(:,:,:) = UNPACK( ZRCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
         ZW(:,:,:) = PTHS(:,:,:)
         PTHS(:,:,:) = UNPACK( ZTHS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
         DEALLOCATE(ZRVT)
         DEALLOCATE(ZRCT)
         DEALLOCATE(ZRVS)
         DEALLOCATE(ZRCS)
         DEALLOCATE(ZTHS)
         DEALLOCATE(ZRHODREF)
         DEALLOCATE(ZZT)
         DEALLOCATE(ZPRES)
         DEALLOCATE(ZEXNREF)
         DEALLOCATE(ZZCPH)
         DEALLOCATE(ZZW)
         DEALLOCATE(ZLVFACT)
         DEALLOCATE(ZRVSATW)
         DEALLOCATE(ZCND)
      END IF ! IMICRO
!
   END IF ! end of adjustment procedure (test on OSUBG_COND)
!
! Remove cloud droplets if there are few

   ZMASK(:,:,:) = 0.0
   ZW(:,:,:) = 0.
   WHERE (PRCS(:,:,:) <= ZRTMIN(2) .OR. PCCS(:,:,:) <= ZCTMIN(2)) 
      PRVS(:,:,:) = PRVS(:,:,:) + PRCS(:,:,:) 
      PTHS(:,:,:) = PTHS(:,:,:) - PRCS(:,:,:)*ZLV(:,:,:)/(ZCPH(:,:,:)*ZEXNS(:,:,:))
      PRCS(:,:,:) = 0.0
      ZW(:,:,:)   = MAX(PCCS(:,:,:),0.)
      PCCS(:,:,:) = 0.0
   END WHERE
!
   ZW1(:,:,:) = 0.
   IF (LWARM .AND. NMOD_CCN.GE.1) ZW1(:,:,:) = SUM(PNAS,DIM=4)
   ZW (:,:,:) = MIN( ZW(:,:,:), ZW1(:,:,:) )
   ZW2(:,:,:) = 0.
   WHERE ( ZW(:,:,:) > 0. )
      ZMASK(:,:,:) = 1.0
      ZW2(:,:,:) = ZW(:,:,:) / ZW1(:,:,:)
   ENDWHERE
!
   IF (LWARM .AND. NMOD_CCN.GE.1) THEN
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
END DO !  end of the iterative loop
!
!
!*       5.2    compute the cloud fraction PCLDFR (binary !!!!!!!)
!
IF ( .NOT. OSUBG_COND ) THEN
   WHERE (PRCS(:,:,:) + PRIS(:,:,:) + PRSS(:,:,:) > 1.E-12 / ZDT)
      PCLDFR(:,:,:)  = 1.
   ELSEWHERE
      PCLDFR(:,:,:)  = 0.
   ENDWHERE
END IF
!
IF ( SIZE(PSRCS,3) /= 0 ) THEN
   WHERE (PRCS(:,:,:) + PRIS(:,:,:) > 1.E-12 / ZDT)
      PSRCS(:,:,:)  = 1.
   ELSEWHERE
      PSRCS(:,:,:)  = 0.
   ENDWHERE
END IF
!
IF ( OSUBG_COND ) THEN
   !
   ! Mixing ratio change (cloud liquid water)
   !
   ZW1(:,:,:) = (ZRC(:,:,:) - PRCS(:,:,:)*PTSTEP) / PTSTEP
   WHERE( ZW1(:,:,:) < 0.0 )
      ZW1(:,:,:) = MAX ( ZW1(:,:,:), -PRCS(:,:,:) )
   ELSEWHERE
      ZW1(:,:,:) = MIN ( ZW1(:,:,:),  PRVS(:,:,:) )
   END WHERE

   WHERE (PCCT(:,:,:) < PCLDFR(:,:,:)*XCTMIN(2) .OR. ZRC(:,:,:)<PCLDFR(:,:,:)*XRTMIN(2))
      ZW1=-PRCS
      PCCS=0.
      PCLDFR=0.
   END WHERE
   
   PRVS(:,:,:)   = PRVS(:,:,:) - ZW1(:,:,:)
   PRCS(:,:,:)   = PRCS(:,:,:) + ZW1(:,:,:)
   PCCS(:,:,:)   = PCCT(:,:,:) / PTSTEP
   PNFS(:,:,:,:) = PNFT(:,:,:,:) / PTSTEP
   PNAS(:,:,:,:) = PNAT(:,:,:,:) / PTSTEP
   PTHS(:,:,:)   = PTHS(:,:,:) +        &
                   ZW1(:,:,:) * ZLV(:,:,:) / (ZCPH(:,:,:) * PEXNREF(:,:,:))
END IF ! fin test OSUBG_COND

IF ( tpfile%lopened ) THEN
   TZFIELD%CMNHNAME   = 'NEB'
   TZFIELD%CSTDNAME   = ''
   TZFIELD%CLONGNAME  = 'NEB'
   TZFIELD%CUNITS     = '1'
   TZFIELD%CDIR       = 'XY'
   TZFIELD%CCOMMENT   = 'X_Y_Z_NEB'
   TZFIELD%NGRID      = 1
   TZFIELD%NTYPE      = TYPEREAL
   TZFIELD%NDIMS      = 3
   TZFIELD%LTIMEDEP   = .TRUE.
   CALL IO_Field_write(TPFILE,TZFIELD,PCLDFR)
END IF
!
!
!*       6.  SAVE CHANGES IN PRS AND PSVS
!            ----------------------------
!
!
! Prepare 3D water mixing ratios
PRS(:,:,:,1) = PRVS(:,:,:)
IF ( KRR .GE. 2 ) PRS(:,:,:,2) = PRCS(:,:,:)
IF ( KRR .GE. 3 ) PRS(:,:,:,3) = PRRS(:,:,:)
IF ( KRR .GE. 4 ) PRS(:,:,:,4) = PRIS(:,:,:)
IF ( KRR .GE. 5 ) PRS(:,:,:,5) = PRSS(:,:,:)
IF ( KRR .GE. 6 ) PRS(:,:,:,6) = PRGS(:,:,:)
!
! Prepare 3D number concentrations
!
IF ( LWARM ) PSVS(:,:,:,NSV_LIMA_NC) = PCCS(:,:,:)
! IF ( LCOLD ) PSVS(:,:,:,NSV_LIMA_NI) = PCIS(:,:,:)
!
IF ( LSCAV .AND. LAERO_MASS ) PSVS(:,:,:,NSV_LIMA_SCAVMASS) = PMAS(:,:,:)
! 
IF ( LWARM .AND. NMOD_CCN .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1) = PNFS(:,:,:,:)
   PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1) = PNAS(:,:,:,:)
END IF
!
! IF ( LCOLD .AND. NMOD_IFN .GE. 1 ) THEN
!    PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1) = PIFS(:,:,:,:)
!    PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1) = PINS(:,:,:,:)
! END IF
!
! IF ( LCOLD .AND. NMOD_IMM .GE. 1 ) THEN
!    PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1) = PNIS(:,:,:,:)
! END IF
!
! write SSI in LFI
!
IF ( tpfile%lopened ) THEN
   ZT(:,:,:) = ( PTHS(:,:,:) * ZDT ) * ZEXNS(:,:,:)
   ZW(:,:,:) = EXP( XALPI - XBETAI/ZT(:,:,:) - XGAMI*ALOG(ZT(:,:,:) ) )
   ZW1(:,:,:)= 2.0*PPABST(:,:,:)-PPABSM(:,:,:)
   ZW(:,:,:) = PRVT(:,:,:)*( ZW1(:,:,:)-ZW(:,:,:) ) / ( (XMV/XMD) * ZW(:,:,:) ) - 1.0
   
   TZFIELD%CMNHNAME   = 'SSI'
   TZFIELD%CSTDNAME   = ''
   TZFIELD%CLONGNAME  = 'SSI'
   TZFIELD%CUNITS     = ''
   TZFIELD%CDIR       = 'XY'
   TZFIELD%CCOMMENT   = 'X_Y_Z_SSI'
   TZFIELD%NGRID      = 1
   TZFIELD%NTYPE      = TYPEREAL
   TZFIELD%NDIMS      = 3
   TZFIELD%LTIMEDEP   = .TRUE.
   CALL IO_Field_write(TPFILE,TZFIELD,ZW)
END IF
!
!
!*       7.  STORE THE BUDGET TERMS
!            ----------------------
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'CEDS', pths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'CEDS', prvs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'CEDS', prcs(:, :, :) * prhodj(:, :, :) )
  !Remark: PRIS is not modified but source term kept for better coherence with lima_adjust and lima_notadjust
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'CEDS', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if ( lwarm ) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'CEDS', pccs(:, :, :) * prhodj(:, :, :) )
    if ( lscav .and. laero_mass ) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_scavmass), 'CEDS', pmas(:, :, :) * prhodj(:, :, :) )
    if ( lwarm ) then
      do jl = 1, nmod_ccn
        idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_free - 1 + jl
        call Budget_store_end( tbudgets(idx), 'CEDS', pnfs(:, :, :, jl) * prhodj(:, :, :) )
        idx = NBUDGET_SV1 - 1 + nsv_lima_ccn_acti - 1 + jl
        call Budget_store_end( tbudgets(idx), 'CEDS', pnas(:, :, :, jl) * prhodj(:, :, :) )
      end do
    end if
!     if ( lcold ) then
!       call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CEDS', pcis(:, :, :) * prhodj(:, :, :) )
!       do jl = 1, nmod_ifn
!         idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_free - 1 + jl
!         call Budget_store_end( tbudgets(idx), 'CEDS', pifs(:, :, :, jl) * prhodj(:, :, :) )
!         idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl - 1 + jl
!         call Budget_store_end( tbudgets(idx), 'CEDS', pins(:, :, :, jl) * prhodj(:, :, :) )
!       end do
!       do jl = 1, nmod_imm
!         idx = NBUDGET_SV1 - 1 + nsv_lima_imm_nucl - 1 + jl
!         call Budget_store_init( tbudgets(idx), 'CEDS', pnis(:, :, :, jl) * prhodj(:, :, :) )
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
