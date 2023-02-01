!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_LIMA_ADJUST
!     #######################
!
INTERFACE
!
      SUBROUTINE LIMA_ADJUST(KRR, KMI, TPFILE,                    &
                             OSUBG_COND, PTSTEP,                  &
                             PRHODREF, PRHODJ, PEXNREF,           &
                             PPABST, PPABSTT,                     &
                             PRT, PRS, PSVT, PSVS,                &
                             PTHS, PSRCS, PCLDFR, PICEFR, PRAINFR )
!
USE MODD_IO,  ONLY: TFILEDATA
USE MODD_NSV, only: NSV_LIMA_BEG
!
INTEGER,                  INTENT(IN)   :: KRR        ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KMI        ! Model index 
TYPE(TFILEDATA),          INTENT(IN)   :: TPFILE     ! Output file
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid 
                                                     ! Condensation
REAL,                     INTENT(IN)   :: PTSTEP     ! Time step          
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODREF  ! Dry density of the 
                                                     ! reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ    ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF   ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST    ! Absolute Pressure at t     
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABSTT   ! Absolute Pressure at t+dt     
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
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   :: PICEFR    ! Cloud fraction          
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   :: PRAINFR   ! Cloud fraction          
!
END SUBROUTINE LIMA_ADJUST
!
END INTERFACE
!
END MODULE MODI_LIMA_ADJUST
!
!     #############################################################
      SUBROUTINE LIMA_ADJUST(KRR, KMI, TPFILE,                    &
                             OSUBG_COND, PTSTEP,                  &
                             PRHODREF, PRHODJ, PEXNREF,           &
                             PPABST, PPABSTT,                     &
                             PRT, PRS, PSVT, PSVS,                &
                             PTHS, PSRCS, PCLDFR, PICEFR, PRAINFR )
!     #############################################################
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
!!      Original             ??/??/13 
!!      C. Barthe  * LACy*   jan. 2014  add budgets
!!      JP Chaboureau *LA*   March 2014  fix the calculation of icy cloud fraction
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet    03/2020: use the new data structures and subroutines for budgets
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet 28/05/2020: bugfix: correct array start for PSVT and PSVS
!  P. Wautelet 01/02/2021: bugfix: add missing CEDS source terms for SV budgets
!  B. Vie         06/2020: fix PSRCS
!  P. Wautelet 23/07/2021: replace non-standard FLOAT function by REAL function
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
use modd_field,            only: TFIELDMETADATA, TYPEREAL
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
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid 
                                                     ! Condensation
REAL,                     INTENT(IN)   :: PTSTEP     ! Time step          
!
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODREF  ! Dry density of the 
                                                     ! reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PRHODJ    ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PEXNREF   ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABST    ! Absolute Pressure at t     
REAL, DIMENSION(:,:,:),   INTENT(IN)   ::  PPABSTT   ! Absolute Pressure at t+dt     
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
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   :: PICEFR    ! Cloud fraction          
REAL, DIMENSION(:,:,:),   INTENT(INOUT)   :: PRAINFR   ! Cloud fraction          
!
!
!*       0.2   Declarations of local variables :
!
! 3D Microphysical variables
REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) &
                         :: PRVT,        & ! Water vapor m.r. at t
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
                            PIFS,        & ! Free      IFN C. source 
                            PINS,        & ! Nucleated IFN C. source
                            PNIS           ! Acti. IMM. nuclei C. source
!
!
!
REAL                     :: ZEPS         ! Mv/Md
REAL                     :: ZDT          ! Time increment (2*Delta t or Delta t if cold start)
REAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) &
                         :: ZEXNS,&      ! guess of the Exner function at t+1
                            ZT,   &      ! guess of the temperature at t+1
                            ZCPH, &      ! guess of the CPh for the mixing
                            ZW,   &
                            ZW1,  &
                            ZW2,  &
                            ZLV,  &      ! guess of the Lv at t+1
                            ZLS,  &      ! guess of the Ls at t+1
                            ZMASK
LOGICAL, DIMENSION(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) &
                         :: GMICRO, GMICRO_RI, GMICRO_RC ! Test where to compute cond/dep proc.
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
TYPE(TFIELDMETADATA)  :: TZFIELD
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
  ITERMAX=2
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
PRVT(:,:,:) = PRT(:,:,:,1)
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
IF ( KRR .GE. 2 ) PRCT(:,:,:) = PRT(:,:,:,2)
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
PCIS(:,:,:) = 0.
!
IF ( LWARM ) PCCT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NC)
IF ( LCOLD ) PCIT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NI)
!
IF ( LWARM ) PCCS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NC)
IF ( LCOLD ) PCIS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NI)
!
IF ( LSCAV .AND. LAERO_MASS ) PMAS(:,:,:) = PSVS(:,:,:,NSV_LIMA_SCAVMASS)
! 
IF ( LWARM .AND. NMOD_CCN.GE.1 ) THEN
   ALLOCATE( PNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ALLOCATE( PNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   PNFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1)
   PNAS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1)
END IF
!
IF ( LCOLD .AND. NMOD_IFN .GE. 1 ) THEN
   ALLOCATE( PIFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
   ALLOCATE( PINS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
   PIFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1)
   PINS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1)
END IF
!
IF ( NMOD_IMM .GE. 1 ) THEN
   ALLOCATE( PNIS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IMM) )
   PNIS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1)
END IF

if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'CEDS', pths(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'CEDS', prvs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'CEDS', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'CEDS', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if ( lwarm .and. nmom_c.ge.2) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc),       'CEDS', pccs(:, :, :) * prhodj(:, :, :) )
    if ( lcold .and. nmom_i.ge.2) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni),       'CEDS', pcis(:, :, :) * prhodj(:, :, :) )
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
    if ( lcold ) then
      do jl = 1, nmod_ifn
        idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_free - 1 + jl
        call Budget_store_init( tbudgets(idx), 'CEDS', pifs(:, :, :, jl) * prhodj(:, :, :) )
        idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl - 1 + jl
        call Budget_store_init( tbudgets(idx), 'CEDS', pins(:, :, :, jl) * prhodj(:, :, :) )
      end do
      do jl = 1, nmod_imm
        idx = NBUDGET_SV1 - 1 + nsv_lima_imm_nucl - 1 + jl
        call Budget_store_init( tbudgets(idx), 'CEDS', pnis(:, :, :, jl) * prhodj(:, :, :) )
      end do
    end if
  end if
end if

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
ZEXNS(:,:,:) = ( PPABSTT(:,:,:) / XP00 ) ** (XRD/XCPD)  
!
!    beginning of the iterative loop
!
DO JITER =1,ITERMAX
!
!*       2.3    compute the intermediate temperature at t+1, T*
!  
  ZT(:,:,:) = ( PTHS(:,:,:) * ZDT ) * ZEXNS(:,:,:)
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
    call Print_msg( NVERB_FATAL, 'GEN', 'LIMA_ADJUST', 'OSUBG_COND=.true. not yet developed' )
  ELSE
!
!-------------------------------------------------------------------------------
!
!
!*       4.     FULLY EXPLICIT SCHEME FROM TZIVION et al. (1989)
!               -----------------------------------------------
! 
!*              select cases where r_i>0 and r_c=0
! 
GMICRO(:,:,:) = .FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                                         &
                         (PRIS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(4) .AND.        &
                          PCIS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(4)      )       &
             .AND. .NOT. (PRCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(2) .AND.        &
                          PCCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(2)      )
GMICRO_RI(:,:,:) = GMICRO(:,:,:)
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
IF( IMICRO >= 1 ) THEN
   ALLOCATE(ZRVT(IMICRO))
   ALLOCATE(ZRIT(IMICRO))
   ALLOCATE(ZCIT(IMICRO))
!
   ALLOCATE(ZRVS(IMICRO))
   ALLOCATE(ZRIS(IMICRO))
   ALLOCATE(ZCIS(IMICRO))  !!!BVIE!!!
   ALLOCATE(ZTHS(IMICRO))
!
   ALLOCATE(ZRHODREF(IMICRO))
   ALLOCATE(ZZT(IMICRO))
   ALLOCATE(ZPRES(IMICRO))
   ALLOCATE(ZEXNREF(IMICRO))
   ALLOCATE(ZZCPH(IMICRO))
   DO JL=1,IMICRO
      ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
      ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
      ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
!
      ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
      ZRIS(JL) = PRIS(I1(JL),I2(JL),I3(JL))
      ZCIS(JL) = PCIS(I1(JL),I2(JL),I3(JL)) !!!BVIE!!!
      ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
!
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
      ZZT(JL) = ZT(I1(JL),I2(JL),I3(JL))
      ZPRES(JL) = PPABSTT(I1(JL),I2(JL),I3(JL))
      ZEXNREF(JL) = PEXNREF(I1(JL),I2(JL),I3(JL))
      ZZCPH(JL) = ZCPH(I1(JL),I2(JL),I3(JL))
   ENDDO
   ALLOCATE(ZZW(IMICRO))
   ALLOCATE(ZLSFACT(IMICRO))
   ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZCPH(:) ! L_s/C_ph
   ALLOCATE(ZRVSATI(IMICRO))
   ALLOCATE(ZRVSATI_PRIME(IMICRO))
   ALLOCATE(ZDELTI(IMICRO))
   ALLOCATE(ZAI(IMICRO))
   ALLOCATE(ZCJ(IMICRO))
   ALLOCATE(ZKA(IMICRO))
   ALLOCATE(ZDV(IMICRO))
   ALLOCATE(ZITI(IMICRO))
!
   ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
   ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
   ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
!
   ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
   ZRVSATI_PRIME(:) = (( XBETAI/ZZT(:) - XGAMI ) / ZZT(:))  &  ! r'_si
                       * ZRVSATI(:) * ( 1. + ZRVSATI(:)/ZEPS )
!
   ZDELTI(:) = ZRVS(:)*ZDT - ZRVSATI(:)
   ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                  + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
   ZZW(:) = MIN(1.E8,( XLBI* MAX(ZCIT(:),XCTMIN(4))                       &
                           /(MAX(ZRIT(:),XRTMIN(4))) )**XLBEXI)
                                                                  ! Lbda_I
   ZITI(:) = ZCIT(:) * (X0DEPI/ZZW(:) + X2DEPI*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDI+2.0)) &
                     / (ZRVSATI(:)*ZAI(:))
!
   ALLOCATE(ZAII(IMICRO))
   ALLOCATE(ZDEP(IMICRO))
!
   ZAII(:) = 1.0 + ZRVSATI_PRIME(:)*ZLSFACT(:)
   ZDEP(:) = 0.0
!
   ZZW(:)  = ZAII(:)*ZITI(:)*ZDT ! R*delta_T
   WHERE( ZZW(:)<1.0E-2 )
      ZDEP(:) = ZITI(:)*ZDELTI(:)*(1.0 - (ZZW(:)/2.0)*(1.0-ZZW(:)/3.0))
   ELSEWHERE          
      ZDEP(:) = ZITI(:)*ZDELTI(:)*(1.0 - EXP(-ZZW(:)))/ZZW(:)
   END WHERE
!
! Integration
!
   WHERE( ZDEP(:) < 0.0 )
      ZDEP(:) = MAX ( ZDEP(:), -ZRIS(:) )
   ELSEWHERE
      ZDEP(:) = MIN ( ZDEP(:),  ZRVS(:) )
!      ZDEP(:) = MIN ( ZDEP(:),  ZCIS(:)*5.E-10 ) !!!BVIE!!!
   END WHERE
   WHERE( ZRIS(:) < ZRTMIN(4) )
      ZDEP(:) = 0.0
   END WHERE
   ZRVS(:) = ZRVS(:) - ZDEP(:)
   ZRIS(:) = ZRIS(:) + ZDEP(:)
   ZTHS(:) = ZTHS(:) + ZDEP(:) * ZLSFACT(:) / ZEXNREF(:)
!
!  Implicit ice crystal sublimation if ice saturated conditions are not met
!
   ZZT(:) = ( ZTHS(:) * ZDT ) * ( ZPRES(:) / XP00 ) ** (XRD/XCPD)
   ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
   WHERE( ZRVS(:)*ZDT<ZRVSATI(:) )
      ZZW(:)  = ZRVS(:) + ZRIS(:)
      ZRVS(:) = MIN( ZZW(:),ZRVSATI(:)/ZDT )
      ZTHS(:) = ZTHS(:) + ( MAX( 0.0,ZZW(:)-ZRVS(:) )-ZRIS(:) ) &
                          * ZLSFACT(:) / ZEXNREF(:)
      ZRIS(:) = MAX( 0.0,ZZW(:)-ZRVS(:) )
   END WHERE
!
!
   ZW(:,:,:) = PRVS(:,:,:)
   PRVS(:,:,:) = UNPACK( ZRVS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRIS(:,:,:)
   PRIS(:,:,:) = UNPACK( ZRIS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PTHS(:,:,:)
   PTHS(:,:,:) = UNPACK( ZTHS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
   DEALLOCATE(ZRVT)
   DEALLOCATE(ZRIT)
   DEALLOCATE(ZCIT)
   DEALLOCATE(ZRVS)
   DEALLOCATE(ZRIS)
   DEALLOCATE(ZCIS) !!!BVIE!!!
   DEALLOCATE(ZTHS)
   DEALLOCATE(ZRHODREF)
   DEALLOCATE(ZZT)
   DEALLOCATE(ZPRES)
   DEALLOCATE(ZEXNREF)
   DEALLOCATE(ZZCPH)
   DEALLOCATE(ZZW)
   DEALLOCATE(ZLSFACT)
   DEALLOCATE(ZRVSATI)
   DEALLOCATE(ZRVSATI_PRIME)
   DEALLOCATE(ZDELTI)
   DEALLOCATE(ZAI)
   DEALLOCATE(ZCJ)
   DEALLOCATE(ZKA)
   DEALLOCATE(ZDV)
   DEALLOCATE(ZITI)
   DEALLOCATE(ZAII)
   DEALLOCATE(ZDEP)
END IF ! IMICRO
!
!
!-------------------------------------------------------------------------------
!
!
!*       5.     FULLY IMPLICIT CONDENSATION SCHEME
!               ---------------------------------
! 
!*              select cases where r_c>0 and r_i=0
! 
!
GMICRO(:,:,:) = .FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =( PRCS(IIB:IIE,IJB:IJE,IKB:IKE)>0. .AND.        &
                                   PCCS(IIB:IIE,IJB:IJE,IKB:IKE)>0.      ) .AND. &
                                   .NOT.GMICRO_RI(IIB:IIE,IJB:IJE,IKB:IKE)
GMICRO_RC(:,:,:) = GMICRO(:,:,:)
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
      ZPRES(JL) = PPABSTT(I1(JL),I2(JL),I3(JL))
      ZEXNREF(JL) = PEXNREF(I1(JL),I2(JL),I3(JL))
      ZZCPH(JL) = ZCPH(I1(JL),I2(JL),I3(JL))
   ENDDO
   ALLOCATE(ZZW(IMICRO))
   ALLOCATE(ZLVFACT(IMICRO))
   ALLOCATE(ZCND(IMICRO))
   ALLOCATE(ZRVSATW(IMICRO))
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
      ZVEC1(:) = MAX( 1.0001, MIN( REAL(NAHEN)-0.0001, XAHENINTP1 * ZZT(:) + XAHENINTP2 ) )
      IVEC1(:) = INT( ZVEC1(:) )
      ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
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
   DEALLOCATE(ZCCS)
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
!
!-------------------------------------------------------------------------------
!
!
!*       6.     IMPLICIT-EXPLICIT SCHEME USING REISIN et al. (1996)
!               ---------------------------------------------------
! 
!*              select cases where r_i>0 and r_c>0 (supercooled water)
! 
!
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                                      &
           .NOT. GMICRO_RI(IIB:IIE,IJB:IJE,IKB:IKE)                    &
     .AND. .NOT. GMICRO_RC(IIB:IIE,IJB:IJE,IKB:IKE)                    &
           .AND. ( PRIS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(4) .AND.       &
                   PCIS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(4)       )     &
           .AND. ( PRCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZRTMIN(2) .AND.       &
                   PCCS(IIB:IIE,IJB:IJE,IKB:IKE)>ZCTMIN(2)       )
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
IF( IMICRO >= 1 ) THEN
   ALLOCATE(ZRVT(IMICRO))
   ALLOCATE(ZRCT(IMICRO))
   ALLOCATE(ZRIT(IMICRO))
   ALLOCATE(ZCCT(IMICRO))
   ALLOCATE(ZCIT(IMICRO))
!
   ALLOCATE(ZRVS(IMICRO))
   ALLOCATE(ZRCS(IMICRO))
   ALLOCATE(ZRIS(IMICRO))
   ALLOCATE(ZCCS(IMICRO))
   ALLOCATE(ZCIS(IMICRO))
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
      ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
      ZCCT(JL) = PCCT(I1(JL),I2(JL),I3(JL))
      ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
!
      ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
      ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
      ZRIS(JL) = PRIS(I1(JL),I2(JL),I3(JL))
      ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
      ZCIS(JL) = PCIS(I1(JL),I2(JL),I3(JL))
      ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
!
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
      ZZT(JL) = ZT(I1(JL),I2(JL),I3(JL))
      ZPRES(JL) = PPABSTT(I1(JL),I2(JL),I3(JL))
      ZEXNREF(JL) = PEXNREF(I1(JL),I2(JL),I3(JL))
      ZZCPH(JL) = ZCPH(I1(JL),I2(JL),I3(JL))
   ENDDO
   ALLOCATE(ZZW(IMICRO))
   ALLOCATE(ZLVFACT(IMICRO))
   ALLOCATE(ZLSFACT(IMICRO))
   ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZCPH(:) ! L_v/C_ph
   ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZCPH(:) ! L_s/C_ph
   ALLOCATE(ZRVSATW(IMICRO))
   ALLOCATE(ZRVSATI(IMICRO))
   ALLOCATE(ZRVSATW_PRIME(IMICRO))
   ALLOCATE(ZRVSATI_PRIME(IMICRO))
   ALLOCATE(ZDELTW(IMICRO))
   ALLOCATE(ZDELTI(IMICRO))
   ALLOCATE(ZAW(IMICRO))
   ALLOCATE(ZAI(IMICRO))
   ALLOCATE(ZCJ(IMICRO))
   ALLOCATE(ZKA(IMICRO))
   ALLOCATE(ZDV(IMICRO))
   ALLOCATE(ZITW(IMICRO))
   ALLOCATE(ZITI(IMICRO))
!
   ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
   ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
   ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
!
!*       6.2    implicit adjustment at water saturation
!
   ZZW(:) = EXP( XALPW - XBETAW/ZZT(:) - XGAMW*ALOG(ZZT(:) ) ) ! es_w
   ZRVSATW(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_sw
   ZRVSATW_PRIME(:) = (( XBETAW/ZZT(:) - XGAMW ) / ZZT(:))  &  ! r'_sw
                      * ZRVSATW(:) * ( 1. + ZRVSATW(:)/ZEPS )
   ZDELTW(:) = ABS( ZRVS(:)*ZDT - ZRVSATW(:) )
   ZAW(:) = ( XLSTT + (XCPV-XCL)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                  + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
   ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
   ZRVSATI_PRIME(:) = (( XBETAI/ZZT(:) - XGAMI ) / ZZT(:))  &  ! r'_si
                      * ZRVSATI(:) * ( 1. + ZRVSATI(:)/ZEPS )
   ZDELTI(:) = ABS( ZRVS(:)*ZDT - ZRVSATI(:) )
   ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                  + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
!
   ZZW(:) = MIN(1.E8,( XLBC* MAX(ZCCT(:),XCTMIN(2))                       &
                           /(MAX(ZRCT(:),XRTMIN(2))) )**XLBEXC)
                                                                  ! Lbda_c
   ZITW(:) = ZCCT(:) * (X0CNDC/ZZW(:) + X2CNDC*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDC+2.0)) &
                     / (ZRVSATW(:)*ZAW(:))
   ZZW(:) = MIN(1.E8,( XLBI* MAX(ZCIT(:),XCTMIN(4))                       &
                     /(MAX(ZRIT(:),XRTMIN(4))) )**XLBEXI)
                                                                  ! Lbda_I
   ZITI(:) = ZCIT(:) * (X0DEPI/ZZW(:) + X2DEPI*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDI+2.0)) &
                     / (ZRVSATI(:)*ZAI(:))
!
   ALLOCATE(ZAWW(IMICRO))
   ALLOCATE(ZAIW(IMICRO))
   ALLOCATE(ZAWI(IMICRO))
   ALLOCATE(ZAII(IMICRO))
!
   ALLOCATE(ZFACT(IMICRO))
   ALLOCATE(ZDELT1(IMICRO))
   ALLOCATE(ZDELT2(IMICRO))
!
   ZAII(:)  = ZITI(:)*ZDELTI(:)
   WHERE( ZAII(:)<1.0E-15 )
      ZFACT(:) = ZLVFACT(:)
   ELSEWHERE          
      ZFACT(:) = (ZLVFACT(:)*ZITW(:)*ZDELTW(:)+ZLSFACT(:)*ZITI(:)*ZDELTI(:)) &
                        / (ZITW(:)*ZDELTW(:)+ZITI(:)*ZDELTI(:))
   END WHERE
   ZAWW(:) = 1.0 + ZRVSATW_PRIME(:)*ZFACT(:)
!
   ZDELT2(:) = (ZRVSATW_PRIME(:)*ZFACT(:)/ZAWW(:)) *                       &
               ( ((-2.*XBETAW+XGAMW*ZZT(:))/(XBETAW-XGAMW*ZZT(:))          &
                 + (XBETAW/ZZT(:)-XGAMW)*(1.0+2.0*ZRVSATW(:)/ZEPS))/ZZT(:) )
   ZDELT1(:) = (ZFACT(:)/ZAWW(:)) * ( ZRVSATW(:) - ZRVS(:)*ZDT )
!
   ALLOCATE(ZCND(IMICRO))
   ALLOCATE(ZDEP(IMICRO))
   ZCND(:) = 0.0
   ZDEP(:) = 0.0
!
   ZZW(:) =  - ZDELT1(:)*( 1.0 + 0.5*ZDELT1(:)*ZDELT2(:) ) / (ZFACT(:)*ZDT) 
   WHERE( ZAII(:)<1.0E-15 )
      ZCND(:) = ZZW(:)
      ZDEP(:) = 0.0
   ELSEWHERE          
      ZCND(:) = ZZW(:)*ZITW(:)*ZDELTW(:) / (ZITW(:)*ZDELTW(:)+ZITI(:)*ZDELTI(:))
      ZDEP(:) = ZZW(:)*ZITI(:)*ZDELTI(:) / (ZITW(:)*ZDELTW(:)+ZITI(:)*ZDELTI(:))
   END WHERE
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
   WHERE( ZDEP(:) < 0.0 )
      ZDEP(:) = MAX ( ZDEP(:), -ZRIS(:) )
   ELSEWHERE
      ZDEP(:) = MIN ( ZDEP(:),  ZRVS(:) )
   END WHERE
   ZRVS(:) = ZRVS(:) - ZDEP(:)
   ZRIS(:) = ZRIS(:) + ZDEP(:)
   ZTHS(:) = ZTHS(:) + ZDEP(:) * ZLSFACT(:) / ZEXNREF(:)
!
!*       6.3    explicit integration of the final eva/dep rates
!
   ZZT(:) = ( ZTHS(:) * ZDT ) * ( ZPRES(:) / XP00 ) ** (XRD/XCPD)
   ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
!
!  If Si < 0, implicit adjustment to Si=0 using ice only
!
   WHERE( ZRVS(:)*ZDT<ZRVSATI(:) )
      ZZW(:)  = ZRVS(:) + ZRIS(:)
      ZRVS(:) = MIN( ZZW(:),ZRVSATI(:)/ZDT )
      ZTHS(:) = ZTHS(:) + ( MAX( 0.0,ZZW(:)-ZRVS(:) )-ZRIS(:) ) &
                          * ZLSFACT(:) / ZEXNREF(:)
      ZRIS(:) = MAX( 0.0,ZZW(:)-ZRVS(:) )
   END WHERE
!
!  Following the previous adjustment, the real procedure begins
!
   ZZT(:) = ( ZTHS(:) * ZDT ) * ( ZPRES(:) / XP00 ) ** (XRD/XCPD)
!
   ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZCPH(:) ! L_v/C_ph
   ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZCPH(:) ! L_s/C_ph
!
   ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
   ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
   ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
!
   ZZW(:) = EXP( XALPW - XBETAW/ZZT(:) - XGAMW*ALOG(ZZT(:) ) ) ! es_w
   ZRVSATW(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_sw
   ZRVSATW_PRIME(:) = (( XBETAW/ZZT(:) - XGAMW ) / ZZT(:))  &  ! r'_sw
                      * ZRVSATW(:) * ( 1. + ZRVSATW(:)/ZEPS )
   ZDELTW(:) = ZRVS(:)*ZDT - ZRVSATW(:)
   ZAW(:) = ( XLSTT + (XCPV-XCL)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                  + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
!
   ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
   ZRVSATI_PRIME(:) = (( XBETAI/ZZT(:) - XGAMI ) / ZZT(:))  &  ! r'_si
                      * ZRVSATI(:) * ( 1. + ZRVSATI(:)/ZEPS )
   ZDELTI(:) = ZRVS(:)*ZDT - ZRVSATI(:)
   ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                 + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
!                    
   ZZW(:) = MIN(1.E8,( XLBC* MAX(ZCCS(:),ZCTMIN(2))                       &
                           /(MAX(ZRCS(:),ZRTMIN(2))) )**XLBEXC)
                                                                  ! Lbda_c
   ZITW(:) = ZCCT(:) * (X0CNDC/ZZW(:) + X2CNDC*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDC+2.0)) &
                     / (ZRVSATW(:)*ZAW(:))
   ZZW(:) = MIN(1.E8,( XLBI* MAX(ZCIS(:),ZCTMIN(4))                       &
                           /(MAX(ZRIS(:),ZRTMIN(4))) )**XLBEXI)
                                                                  ! Lbda_I
   ZITI(:) = ZCIT(:) * (X0DEPI/ZZW(:) + X2DEPI*ZCJ(:)*ZCJ(:)/ZZW(:)**(XDI+2.0)) &
                     / (ZRVSATI(:)*ZAI(:))
!
   ZAWW(:) = 1.0 + ZRVSATW_PRIME(:)*ZLVFACT(:)
   ZAIW(:) = 1.0 + ZRVSATI_PRIME(:)*ZLVFACT(:)
   ZAWI(:) = 1.0 + ZRVSATW_PRIME(:)*ZLSFACT(:)
   ZAII(:) = 1.0 + ZRVSATI_PRIME(:)*ZLSFACT(:)
!
   ZCND(:) = 0.0      
   ZDEP(:) = 0.0
   ZZW(:) = ZAWW(:)*ZITW(:) + ZAII(:)*ZITI(:) ! R
   WHERE( ZZW(:)<1.0E-2 )
      ZFACT(:) = ZDT*(0.5 - (ZZW(:)*ZDT)/6.0)
   ELSEWHERE          
      ZFACT(:) = (1.0/ZZW(:))*(1.0-(1.0-EXP(-ZZW(:)*ZDT))/(ZZW(:)*ZDT))
   END WHERE
   ZCND(:) = ZITW(:)*(ZDELTW(:)-( ZAWW(:)*ZITW(:)*ZDELTW(:)           &
                                + ZAWI(:)*ZITI(:)*ZDELTI(:) )*ZFACT(:))
   ZDEP(:) = ZITI(:)*(ZDELTI(:)-( ZAIW(:)*ZITW(:)*ZDELTW(:)           &
                                + ZAII(:)*ZITI(:)*ZDELTI(:) )*ZFACT(:))
!                    
! Integration        
!           
   WHERE( ZCND(:) < 0.0 )
      ZCND(:) = MAX ( ZCND(:), -ZRCS(:) )
   ELSEWHERE          
      ZCND(:) = MIN ( ZCND(:),  ZRVS(:) )
   END WHERE
   WHERE( ZRCS(:) < ZRTMIN(2) )
      ZCND(:) = 0.0
   END WHERE
   ZRVS(:) = ZRVS(:) - ZCND(:)
   ZRCS(:) = ZRCS(:) + ZCND(:)
   ZTHS(:) = ZTHS(:) + ZCND(:) * ZLVFACT(:) / ZEXNREF(:)
!                    
   WHERE( ZDEP(:) < 0.0 )
      ZDEP(:) = MAX ( ZDEP(:), -ZRIS(:) )
   ELSEWHERE          
      ZDEP(:) = MIN ( ZDEP(:),  ZRVS(:) )
   END WHERE
   WHERE( ZRIS(:) < ZRTMIN(4) )
      ZDEP(:) = 0.0
   END WHERE
   ZRVS(:) = ZRVS(:) - ZDEP(:)
   ZRIS(:) = ZRIS(:) + ZDEP(:)
   ZTHS(:) = ZTHS(:) + ZDEP(:) * ZLSFACT(:) / ZEXNREF(:)
!
!  Implicit ice crystal sublimation if ice saturated conditions are not met
!
   ZZT(:) = ( ZTHS(:) * ZDT ) * ( ZPRES(:) / XP00 ) ** (XRD/XCPD)
   ZLVFACT(:) = (XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZCPH(:) ! L_v/C_ph
   ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZCPH(:) ! L_s/C_ph   
   ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
   ZRVSATI(:) = ZEPS*ZZW(:) / ( ZPRES(:)-ZZW(:) )              ! r_si
   WHERE( ZRVS(:)*ZDT<ZRVSATI(:) )
      ZZW(:)  = ZRVS(:) + ZRIS(:)
      ZRVS(:) = MIN( ZZW(:),ZRVSATI(:)/ZDT )
      ZTHS(:) = ZTHS(:) + ( MAX( 0.0,ZZW(:)-ZRVS(:) )-ZRIS(:) ) &
                                            * ZLSFACT(:) / ZEXNREF(:)
      ZRIS(:) = MAX( 0.0,ZZW(:)-ZRVS(:) )
   END WHERE
!
!                    
!                    
   ZW(:,:,:) = PRVS(:,:,:)
   PRVS(:,:,:) = UNPACK( ZRVS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRCS(:,:,:)
   PRCS(:,:,:) = UNPACK( ZRCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRIS(:,:,:)
   PRIS(:,:,:) = UNPACK( ZRIS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PTHS(:,:,:)
   PTHS(:,:,:) = UNPACK( ZTHS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
   DEALLOCATE(ZRVT)
   DEALLOCATE(ZRCT)
   DEALLOCATE(ZRIT)
   DEALLOCATE(ZCCT)
   DEALLOCATE(ZCIT)
   DEALLOCATE(ZRVS)
   DEALLOCATE(ZRCS)
   DEALLOCATE(ZRIS)
   DEALLOCATE(ZCCS)
   DEALLOCATE(ZCIS)
   DEALLOCATE(ZTHS)
   DEALLOCATE(ZRHODREF)
   DEALLOCATE(ZZT)
   DEALLOCATE(ZPRES)
   DEALLOCATE(ZEXNREF)
   DEALLOCATE(ZZCPH)
   DEALLOCATE(ZZW)
   DEALLOCATE(ZLVFACT)
   DEALLOCATE(ZLSFACT)
   DEALLOCATE(ZRVSATW)
   DEALLOCATE(ZRVSATI)
   DEALLOCATE(ZRVSATW_PRIME)
   DEALLOCATE(ZRVSATI_PRIME)
   DEALLOCATE(ZDELTW)
   DEALLOCATE(ZDELTI)
   DEALLOCATE(ZAW)
   DEALLOCATE(ZAI)
   DEALLOCATE(ZCJ)
   DEALLOCATE(ZKA)
   DEALLOCATE(ZDV)
   DEALLOCATE(ZITW)
   DEALLOCATE(ZITI)
   DEALLOCATE(ZAWW)
   DEALLOCATE(ZAIW)
   DEALLOCATE(ZAWI)
   DEALLOCATE(ZAII)
   DEALLOCATE(ZFACT)
   DEALLOCATE(ZDELT1)
   DEALLOCATE(ZDELT2)
   DEALLOCATE(ZCND)
   DEALLOCATE(ZDEP)
END IF ! IMICRO
!
END IF ! OSUBG_COND
!
! full sublimation of the cloud ice crystals if there are few
!
IF ( .NOT. OSUBG_COND ) THEN

ZMASK(:,:,:) = 0.0
ZW(:,:,:) = 0.
WHERE (PRIS(:,:,:) <= ZRTMIN(4) .OR. PCIS(:,:,:) <= ZCTMIN(4)) 
   PRVS(:,:,:) = PRVS(:,:,:) + PRIS(:,:,:) 
   PTHS(:,:,:) = PTHS(:,:,:) - PRIS(:,:,:)*ZLS(:,:,:)/(ZCPH(:,:,:)*ZEXNS(:,:,:))
   PRIS(:,:,:) = 0.0
   ZW(:,:,:)   = MAX(PCIS(:,:,:),0.)
   PCIS(:,:,:) = 0.0
END WHERE
!
IF (LCOLD .AND. (NMOD_IFN .GE. 1 .OR. NMOD_IMM .GE. 1)) THEN
   ZW1(:,:,:) = 0.
   IF (NMOD_IFN .GE. 1) ZW1(:,:,:) = ZW1(:,:,:) + SUM(PINS,DIM=4)
   IF (NMOD_IMM .GE. 1) ZW1(:,:,:) = ZW1(:,:,:) + SUM(PNIS,DIM=4)
   ZW (:,:,:) = MIN( ZW(:,:,:), ZW1(:,:,:) )
   ZW2(:,:,:) = 0.
   WHERE ( ZW(:,:,:) > 0. )
      ZMASK(:,:,:) = 1.0
      ZW2(:,:,:) = ZW(:,:,:) / ZW1(:,:,:)
   ENDWHERE
END IF
!
IF (LCOLD .AND. NMOD_IFN.GE.1) THEN
   DO JMOD_IFN = 1, NMOD_IFN
      PIFS(:,:,:,JMOD_IFN) = PIFS(:,:,:,JMOD_IFN) +                    &
           ZMASK(:,:,:) * PINS(:,:,:,JMOD_IFN) * ZW2(:,:,:)
      PINS(:,:,:,JMOD_IFN) = PINS(:,:,:,JMOD_IFN) -                    &
           ZMASK(:,:,:) * PINS(:,:,:,JMOD_IFN) * ZW2(:,:,:)
      PINS(:,:,:,JMOD_IFN) = MAX( 0.0 , PINS(:,:,:,JMOD_IFN) )
   ENDDO
END IF
!
IF (LCOLD .AND. NMOD_IMM.GE.1) THEN
   JMOD_IMM = 0
   DO JMOD = 1, NMOD_CCN
      IF (NIMM(JMOD) == 1) THEN 
         JMOD_IMM = JMOD_IMM + 1 
         PNAS(:,:,:,JMOD)     = PNAS(:,:,:,JMOD) +                     &
              ZMASK(:,:,:) * PNIS(:,:,:,JMOD_IMM) * ZW2(:,:,:)
         PNIS(:,:,:,JMOD_IMM) = PNIS(:,:,:,JMOD_IMM) -                 &
              ZMASK(:,:,:) * PNIS(:,:,:,JMOD_IMM) * ZW2(:,:,:)
         PNIS(:,:,:,JMOD_IMM) = MAX( 0.0 , PNIS(:,:,:,JMOD_IMM) )
      END IF
   ENDDO
END IF
!
! complete evaporation of the cloud droplets if there are few
!
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
!  end of the iterative loop
!
END IF ! .NOT.OSUBG_COND

END DO
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
   WHERE (PRIS(:,:,:) > 1.E-12 / ZDT)
      PICEFR(:,:,:)  = 1.
   ELSEWHERE
      PICEFR(:,:,:)  = 0.
   ENDWHERE
   WHERE (PRRS(:,:,:)+PRSS(:,:,:)+PRGS(:,:,:) > 1.E-12 / ZDT)
      PRAINFR(:,:,:)  = 1.
   ELSEWHERE
      PRAINFR(:,:,:)  = 0.
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
IF ( tpfile%lopened ) THEN
  TZFIELD = TFIELDMETADATA(   &
    CMNHNAME   = 'NEB',       &
    CSTDNAME   = '',          &
    CLONGNAME  = 'NEB',       &
    CUNITS     = '1',         &
    CDIR       = 'XY',        &
    CCOMMENT   = 'X_Y_Z_NEB', &
    NGRID      = 1,           &
    NTYPE      = TYPEREAL,    &
    NDIMS      = 3,           &
    LTIMEDEP   = .TRUE.       )
  CALL IO_Field_write(TPFILE,TZFIELD,ZW)
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
IF ( LCOLD ) PSVS(:,:,:,NSV_LIMA_NI) = PCIS(:,:,:)
!
IF ( LSCAV .AND. LAERO_MASS ) PSVS(:,:,:,NSV_LIMA_SCAVMASS) = PMAS(:,:,:)
! 
IF ( LWARM .AND. NMOD_CCN .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1) = PNFS(:,:,:,:)
   PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1) = PNAS(:,:,:,:)
END IF
!
IF ( LCOLD .AND. NMOD_IFN .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1) = PIFS(:,:,:,:)
   PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1) = PINS(:,:,:,:)
END IF
!
IF ( LCOLD .AND. NMOD_IMM .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1) = PNIS(:,:,:,:)
END IF
!
! write SSI in LFI
!
IF ( tpfile%lopened ) THEN
  ZT(:,:,:) = ( PTHS(:,:,:) * ZDT ) * ZEXNS(:,:,:)
  ZW(:,:,:) = EXP( XALPI - XBETAI/ZT(:,:,:) - XGAMI*ALOG(ZT(:,:,:) ) )
  ZW1(:,:,:)= PPABSTT(:,:,:)
  ZW(:,:,:) = PRVT(:,:,:)*( ZW1(:,:,:)-ZW(:,:,:) ) / ( (XMV/XMD) * ZW(:,:,:) ) - 1.0

  TZFIELD = TFIELDMETADATA(   &
    CMNHNAME   = 'SSI',       &
    CSTDNAME   = '',          &
    CLONGNAME  = 'SSI',       &
    CUNITS     = '',          &
    CDIR       = 'XY',        &
    CCOMMENT   = 'X_Y_Z_SSI', &
    NGRID      = 1,           &
    NTYPE      = TYPEREAL,    &
    NDIMS      = 3,           &
    LTIMEDEP   = .TRUE.       )
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
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'CEDS', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    if ( lwarm .and. nmom_c.ge.2) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc),       'CEDS', pccs(:, :, :) * prhodj(:, :, :) )
    if ( lcold .and. nmom_i.ge.2) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni),       'CEDS', pcis(:, :, :) * prhodj(:, :, :) )
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
    if ( lcold ) then
      do jl = 1, nmod_ifn
        idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_free - 1 + jl
        call Budget_store_end( tbudgets(idx), 'CEDS', pifs(:, :, :, jl) * prhodj(:, :, :) )
        idx = NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl - 1 + jl
        call Budget_store_end( tbudgets(idx), 'CEDS', pins(:, :, :, jl) * prhodj(:, :, :) )
      end do
      do jl = 1, nmod_imm
        idx = NBUDGET_SV1 - 1 + nsv_lima_imm_nucl - 1 + jl
        call Budget_store_end( tbudgets(idx), 'CEDS', pnis(:, :, :, jl) * prhodj(:, :, :) )
      end do
    end if
  end if
end if
!++cb++
DEALLOCATE(ZRTMIN)
DEALLOCATE(ZCTMIN)
IF (ALLOCATED(PNFS)) DEALLOCATE(PNFS)
IF (ALLOCATED(PNAS)) DEALLOCATE(PNAS)
IF (ALLOCATED(PIFS)) DEALLOCATE(PIFS)
IF (ALLOCATED(PINS)) DEALLOCATE(PINS)
IF (ALLOCATED(PNIS)) DEALLOCATE(PNIS)
!--cb--
!
!------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_ADJUST
