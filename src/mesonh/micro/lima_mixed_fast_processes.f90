!MNH_LIC Copyright 2013-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #####################################
       MODULE MODI_LIMA_MIXED_FAST_PROCESSES
!      #####################################
!
INTERFACE
      SUBROUTINE LIMA_MIXED_FAST_PROCESSES (PRHODREF, PZT, PPRES, PTSTEP,                   &
                                            PLSFACT, PLVFACT, PKA, PDV, PCJ,                &
                                            PRVT1D, PRCT1D, PRRT1D, PRIT1D, PRST1D, PRGT1D, &
                                            PRHT1D, PCCT1D, PCRT1D, PCIT1D, PCST1D, PCGT1D, PCHT1D, &
                                            PRCS1D, PRRS1D, PRIS1D, PRSS1D, PRGS1D, PRHS1D, &
                                            PTHS1D, PCCS1D, PCRS1D, PCIS1D, PCSS1D, PCGS1D, PCHS1D, &
                                            PLBDAC, PLBDAR, PLBDAI, PLBDAS, PLBDAG, PLBDAH, &
                                            PRHODJ1D, GMICRO, PRHODJ, KMI, PTHS,            &
                                            PRCS, PRRS, PRIS, PRSS, PRGS, PRHS,             &
                                            PCCS, PCRS, PCIS, PCSS, PCGS, PCHS              )
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(:),   INTENT(IN)    :: PZT       ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PPRES     ! Pressure
REAL,                 INTENT(IN)    :: PTSTEP    ! Time step          
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PKA       ! Thermal conductivity of the air
REAL, DIMENSION(:),   INTENT(IN)    :: PDV       ! Diffusivity of water vapor in the air
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ       ! Ventilation coefficient ?
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT1D    ! Water vapor m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT1D    ! Cloud water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT1D    ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT1D    ! Pristine ice m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRST1D    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT1D    ! Graupel m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRHT1D    ! Hail m.r. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT1D    ! Cloud water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT1D    ! Rain water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT1D    ! Pristine ice conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCST1D    ! Snow/aggregate conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCGT1D    ! Graupel conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCHT1D    ! Hail conc. at t
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PRCS1D    ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRRS1D    ! Rain water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRIS1D    ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRSS1D    ! Snow/aggregate m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRGS1D    ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRHS1D    ! Hail m.r. source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PTHS1D    ! Theta source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PCCS1D    ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCRS1D    ! Rain water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCIS1D    ! Pristine ice conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCSS1D    ! Snow/aggregate conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCGS1D    ! Graupel conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCHS1D    ! Hail conc. source
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAC  ! Slope param of the cloud droplet distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAR  ! Slope param of the raindrop  distr
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAI  ! Slope param of the ice distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAS  ! Slope param of the aggregate distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAG  ! Slope param of the graupel distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAH  ! Slope param of the hail distr.
!
! used for budget storage
REAL,    DIMENSION(:),     INTENT(IN) :: PRHODJ1D
LOGICAL, DIMENSION(:,:,:), INTENT(IN) :: GMICRO 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHODJ
INTEGER,                   INTENT(IN) :: KMI 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PTHS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRRS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRIS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRSS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRGS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCRS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCIS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCSS 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCGS 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCHS
!
END SUBROUTINE LIMA_MIXED_FAST_PROCESSES
END INTERFACE
END MODULE MODI_LIMA_MIXED_FAST_PROCESSES
!
!     ###############################################################################
      SUBROUTINE LIMA_MIXED_FAST_PROCESSES (PRHODREF, PZT, PPRES, PTSTEP,                   &
                                            PLSFACT, PLVFACT, PKA, PDV, PCJ,                &
                                            PRVT1D, PRCT1D, PRRT1D, PRIT1D, PRST1D, PRGT1D, &
                                            PRHT1D, PCCT1D, PCRT1D, PCIT1D, PCST1D, PCGT1D, PCHT1D, &
                                            PRCS1D, PRRS1D, PRIS1D, PRSS1D, PRGS1D, PRHS1D, &
                                            PTHS1D, PCCS1D, PCRS1D, PCIS1D, PCSS1D, PCGS1D, PCHS1D, &
                                            PLBDAC, PLBDAR, PLBDAI, PLBDAS, PLBDAG, PLBDAH, &
                                            PRHODJ1D, GMICRO, PRHODJ, KMI, PTHS,            &
                                            PRCS, PRRS, PRIS, PRSS, PRGS, PRHS,             &
                                            PCCS, PCRS, PCIS, PCSS, PCGS, PCHS              )
!     ###############################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the mixed-phase 
!!    fast processes :
!!      
!!      - Fast RS processes :
!!          - Cloud droplet riming of the aggregates
!!          - Hallett-Mossop ice multiplication process due to snow riming
!!          - Rain accretion onto the aggregates
!!          - Conversion-Melting of the aggregates
!!
!!      - Fast RG processes :
!!          - Rain contact freezing
!!          - Wet/Dry growth of the graupel
!!          - Hallett-Mossop ice multiplication process due to graupel riming
!!          - Melting of the graupeln
!!
!!
!!**  METHOD
!!    ------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Most of the parameterizations come from the ICE3 scheme, described in
!!    the MESO-NH scientific documentation.
!!
!!      Cohard, J.-M. and J.-P. Pinty, 2000: A comprehensive two-moment warm 
!!      microphysical bulk scheme. 
!!        Part I: Description and tests
!!        Part II: 2D experiments with a non-hydrostatic model
!!      Accepted for publication in Quart. J. Roy. Meteor. Soc. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      C. Barthe  * LACy *  jan. 2014    add budgets
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet    03/2020: use the new data structures and subroutines for budgets
!  C. Barthe   14/03/2022: - add CIBU (from T. Hoarau's work) and RDSF (from J.P. Pinty's work)
!                          - change the name of some arguments to match the DOCTOR norm
!                          - change conditions for HMG to occur
!  J. Wurtz       03/2022: new snow characteristics
!  M. Taufour     07/2022: add concentration for snow, graupel, hail
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,           only: lbu_enable, nbumod,                                                                              &
                                 lbudget_th, lbudget_rc, lbudget_rr, lbudget_ri, lbudget_rs, lbudget_rg, lbudget_rh, lbudget_sv,  &
                                 NBUDGET_TH, NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, NBUDGET_SV1, &
                                 tbudgets
USE MODD_CST
USE MODD_NSV
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_COLD
USE MODD_PARAM_LIMA_MIXED
USE MODD_PARAM_LIMA_WARM, ONLY : XBR, XDR

use mode_budget,           only: Budget_store_init, Budget_store_end

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(:),   INTENT(IN)    :: PZT       ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PPRES     ! Pressure
REAL,                 INTENT(IN)    :: PTSTEP    ! Time step          
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: PKA       ! Thermal conductivity of the air
REAL, DIMENSION(:),   INTENT(IN)    :: PDV       ! Diffusivity of water vapor in the air
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ       ! Ventilation coefficient ?
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT1D    ! Water vapor m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT1D    ! Cloud water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT1D    ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT1D    ! Pristine ice m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRST1D    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT1D    ! Graupel m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRHT1D    ! Hail m.r. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT1D    ! Cloud water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT1D    ! Rain water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT1D    ! Pristine ice conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCST1D    ! Snow/aggregate conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCGT1D    ! Graupel conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCHT1D    ! Hail conc. at t
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PRCS1D    ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRRS1D    ! Rain water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRIS1D    ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRSS1D    ! Snow/aggregate m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRGS1D    ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PRHS1D    ! Hail m.r. source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PTHS1D    ! Theta source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PCCS1D    ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCRS1D    ! Rain water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCIS1D    ! Pristine ice conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCSS1D    ! Snow/aggregate conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCGS1D    ! Graupel conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: PCHS1D    ! Hail conc. source
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAC  ! Slope param of the cloud droplet distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAR  ! Slope param of the raindrop  distr
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAI  ! Slope param of the ice distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAS  ! Slope param of the aggregate distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAG  ! Slope param of the graupel distr.
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDAH  ! Slope param of the hail distr.
!
! used for budget storage
REAL,    DIMENSION(:),     INTENT(IN) :: PRHODJ1D
LOGICAL, DIMENSION(:,:,:), INTENT(IN) :: GMICRO 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHODJ
INTEGER,                   INTENT(IN) :: KMI 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PTHS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRRS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRIS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRSS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRGS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCRS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCIS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCSS 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCGS 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCHS
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PZT)) :: GRIM, GACC, GDRY, GWET, GHAIL ! Test where to compute
INTEGER :: IGRIM, IGACC, IGDRY, IGWET, IHAIL
INTEGER :: JJ
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1,IVEC2        ! Vectors of indices
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1,ZVEC2, ZVEC3 ! Work vectors
REAL,    DIMENSION(SIZE(PZT))  :: ZZW, ZZX, ZZNW
REAL,    DIMENSION(SIZE(PRIT1D)) :: ZZDI, ZZDC
REAL,    DIMENSION(SIZE(PZT))  :: ZRDRYG, ZRWETG, ZNWETG
REAL,    DIMENSION(SIZE(PZT),7)  :: ZZW1, ZZNW1 
REAL :: NHAIL
REAL :: ZTHRH, ZTHRC
!
! Variables for CIBU
LOGICAL, DIMENSION(SIZE(PZT)) :: GCIBU ! Test where to compute collision process
LOGICAL, SAVE                 :: GFIRSTCALL = .TRUE. ! control switch for the first call
!
INTEGER                            :: ICIBU
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC2_S1,IVEC2_S2         ! Snow indice vector
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC2_G                   ! Graupel indice vector
INTEGER, PARAMETER                 :: I_SEED_PARAM = 26032012
INTEGER, DIMENSION(:), ALLOCATABLE :: I_SEED
INTEGER                            :: NI_SEED
!
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_S, ZVEC1_SW, ZVEC1_S1, ZVEC1_S2,  & ! Work vectors
                                      ZVEC1_S3, ZVEC1_S4,           &
                                      ZVEC1_S11, ZVEC1_S12,         & ! for snow
                                      ZVEC1_S21, ZVEC1_S22,         &
                                      ZVEC1_S31, ZVEC1_S32,         &
                                      ZVEC1_S41, ZVEC1_S42,         &
                                      ZVEC2_S1, ZVEC2_S2
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_G, ZVEC1_G1, ZVEC1_G2, & ! Work vectors
                                      ZVEC2_G                        ! for graupel
REAL,    DIMENSION(:), ALLOCATABLE :: ZINTG_SNOW_1, & ! incomplete gamma function
                                      ZINTG_SNOW_2, & ! for snow
                                      ZINTG_SNOW_3, &
                                      ZINTG_SNOW_4
REAL,    DIMENSION(:), ALLOCATABLE :: ZINTG_GRAUPEL_1, &  ! incomplete gamma
                                      ZINTG_GRAUPEL_2     ! function for graupel
REAL,    DIMENSION(:), ALLOCATABLE :: ZNI_CIBU, ZRI_CIBU  ! CIBU rates
REAL,    DIMENSION(:), ALLOCATABLE :: ZFRAGMENTS, ZHARVEST, ZFRAG_CIBU
REAL                               :: ZFACT1_XNDEBRIS, ZFACT2_XNDEBRIS
!
LOGICAL, DIMENSION(SIZE(PZT))      :: GRDSF              ! Test where to compute collision process
INTEGER :: IRDSF
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_R            ! Work vectors for rain
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1_R1           ! Work vectors for rain
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC2_R            ! Work vectors for rain
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC2_R            ! Rain indice vector
REAL,    DIMENSION(:), ALLOCATABLE :: ZINTG_RAIN         ! incomplete gamma function for rain
REAL,    DIMENSION(:), ALLOCATABLE :: ZNI_RDSF,ZRI_RDSF  ! RDSF rates
!
REAL,    DIMENSION(:),   ALLOCATABLE :: ZAUX     ! used to distribute
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZFACT    ! the total concentration in each shape
REAL,    DIMENSION(:),   ALLOCATABLE :: ZONEOVER_VAR ! for optimization
LOGICAL :: M2_ICE
!
!
!-------------------------------------------------------------------------------
!
M2_ICE = NMOM_S.GE.2 .AND. NMOM_G.GE.2
IF (LHAIL) M2_ICE = M2_ICE .AND. NMOM_H.GE.2
!
!                         #################
!                         FAST RS PROCESSES
!                         #################
!
SNOW: IF (LSNOW) THEN
!
!
!*       1.1  Cloud droplet riming of the aggregates  
!        -------------------------------------------
!
ZZW1(:,:) = 0.0
!
GRIM(:) = (PRCT1D(:)>XRTMIN(2)) .AND. (PRST1D(:)>XRTMIN(5)) .AND. (PRCS1D(:)>XRTMIN(2)/PTSTEP) .AND. (PZT(:)<XTT)
IF (NMOM_S.GE.2) GRIM(:) = GRIM(:) .AND. (PCST1D(:)>XCTMIN(5))
IGRIM = COUNT( GRIM(:) )
!
IF( IGRIM>0 ) THEN
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'RIM', &
                                            Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'RIM', &
                                            Unpack( prcs1d(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'RIM', &
                                            Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'RIM', &
                                            Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'RIM', &
                                            Unpack( pccs1d(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv .and. NMOM_S.GE.2 ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'RIM', &
                                            Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv .and. NMOM_G.GE.2 ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'RIM', &
                                            Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) )
  end if
!
!        1.1.0  allocations
!
  ALLOCATE(ZVEC1(IGRIM))
  ALLOCATE(ZVEC2(IGRIM))
  ALLOCATE(IVEC1(IGRIM))
  ALLOCATE(IVEC2(IGRIM))
!
!        1.1.1  select the PLBDAS
!
  ZVEC1(:) = PACK( PLBDAS(:),MASK=GRIM(:) )
!
!        1.1.2  find the next lower indice for the PLBDAS in the geometrical
!               set of Lbda_s used to tabulate some moments of the incomplete 
!               gamma function
!
  ZVEC2(1:IGRIM) = MAX( 1.0001, MIN( REAL(NGAMINC)-0.0001,           &
                         XRIMINTP1 * LOG( ZVEC1(1:IGRIM) ) + XRIMINTP2 ) )
  IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
  ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - REAL( IVEC2(1:IGRIM) )
!
!        1.1.3  perform the linear interpolation of the normalized
!               "2+XDS"-moment of the incomplete gamma function
!
  ZVEC1(1:IGRIM) =   XGAMINC_RIM1( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                   - XGAMINC_RIM1( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
  ZZW(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 )
!
!        1.1.4  riming of the small sized aggregates
!
  WHERE ( GRIM(:) )
     ZZW1(:,1) = MIN( PRCS1D(:),                          &
  	           XCRIMSS * ZZW(:) * PRCT1D(:) * PCST1D(:)      & 
  	                            *   PLBDAS(:)**XEXCRIMSS &
               * (1+(XFVELOS/PLBDAS(:))**XALPHAS)**(-XNUS+XEXCRIMSS/XALPHAS) &
    			            * PRHODREF(:)**(-XCEXVT) )
     PRCS1D(:) = PRCS1D(:) - ZZW1(:,1)
     PRSS1D(:) = PRSS1D(:) + ZZW1(:,1)
     PTHS1D(:) = PTHS1D(:) + ZZW1(:,1)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(RCRIMSS))
!
     PCCS1D(:) = MAX( PCCS1D(:)-ZZW1(:,1)*(PCCT1D(:)/PRCT1D(:)),0.0 ) ! Lambda_c**3
  END WHERE
!
!        1.1.5  perform the linear interpolation of the normalized
!               "XBS"-moment of the incomplete gamma function
!
  ZVEC1(1:IGRIM) =  XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                  - XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
  ZZW(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 )
!
!        1.1.6  riming-conversion of the large sized aggregates into graupeln
!
  WHERE ( GRIM(:) .AND. (PRSS1D(:)>XRTMIN(5)/PTSTEP) .AND. (PCSS1D(:)>XCTMIN(5)/PTSTEP))
     ZZW1(:,2) = MIN( PRCS1D(:),                 &
    	           XCRIMSG * PRCT1D(:) * PCST1D(:)      & ! RCRIMSG
    	           *  PLBDAS(:)**XEXCRIMSG*(1+(XFVELOS/PLBDAS(:))**XALPHAS)**(-XNUS+XEXCRIMSG/XALPHAS)  &
  	                   * PRHODREF(:)**(-XCEXVT) &
    		           - ZZW1(:,1)              )
     ZZW1(:,3) = MIN( PRSS1D(:),  PCST1D(:) *          &
                       XSRIMCG * PLBDAS(:)**XEXSRIMCG   & ! RSRIMCG 
   	                       * (1.0 - ZZW(:) )/PTSTEP)
     PRCS1D(:) = PRCS1D(:) - ZZW1(:,2)
     PRSS1D(:) = PRSS1D(:) - ZZW1(:,3)
     PRGS1D(:) = PRGS1D(:) + ZZW1(:,2) + ZZW1(:,3)
     PTHS1D(:) = PTHS1D(:) + ZZW1(:,2)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(RCRIMSG))
!
     PCCS1D(:) = MAX( PCCS1D(:)-ZZW1(:,2)*(PCCT1D(:)/PRCT1D(:)),0.0 ) ! Lambda_c**3
     PCSS1D(:) = MAX( PCSS1D(:)-ZZW1(:,3)*(PCST1D(:)/PRST1D(:)),0.0 )
     PCGS1D(:) = MAX( PCGS1D(:)+ZZW1(:,3)*(PCST1D(:)/PRST1D(:)),0.0 ) !
  END WHERE
  DEALLOCATE(IVEC2)
  DEALLOCATE(IVEC1)
  DEALLOCATE(ZVEC2)
  DEALLOCATE(ZVEC1)
  !
  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'RIM', &
                                           Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'RIM', &
                                           Unpack( prcs1d(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'RIM', &
                                           Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'RIM', &
                                           Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'RIM', &
                                           Unpack( pccs1d(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv .and. NMOM_S.GE.2) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'RIM', &
                                           Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv .and. NMOM_G.GE.2) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'RIM', &
                                           Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) )
  end if
END IF
!
!
!*       1.2  Hallett-Mossop ice multiplication process due to snow riming  
!        -----------------------------------------------------------------
!
GRIM(:) = (PZT(:)<XHMTMAX) .AND. (PZT(:)>XHMTMIN)                          &
                           .AND. (PRST1D(:)>XRTMIN(5)) .AND. (PRCT1D(:)>XRTMIN(2))
IGRIM = COUNT( GRIM(:) )
IF( IGRIM>0 ) THEN
  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'HMS', &
                                            Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'HMS', &
                                            Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HMS', &
                                            Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if

  ALLOCATE(ZVEC1(IGRIM))
  ALLOCATE(ZVEC2(IGRIM))
  ALLOCATE(IVEC2(IGRIM))
!
  ZVEC1(:) = PACK( PLBDAC(:),MASK=GRIM(:) )
  ZVEC2(1:IGRIM) = MAX( 1.0001, MIN( REAL(NGAMINC)-0.0001,           &
                        XHMLINTP1 * LOG( ZVEC1(1:IGRIM) ) + XHMLINTP2 ) )
  IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
  ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - REAL( IVEC2(1:IGRIM) )
  ZVEC1(1:IGRIM) =   XGAMINC_HMC( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                   - XGAMINC_HMC( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
  ZZX(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 ) ! Large droplets
!
  WHERE ( GRIM(:) .AND. ZZX(:)<0.99 )
    ZZW1(:,5) = (ZZW1(:,1)+ZZW1(:,2))*(PCCT1D(:)/PRCT1D(:))*(1.0-ZZX(:))* & 
                                                           XHM_FACTS* &
         MAX( 0.0, MIN( (PZT(:)-XHMTMIN)/3.0,(XHMTMAX-PZT(:))/2.0 ) ) ! CCHMSI
    PCIS1D(:) = PCIS1D(:) + ZZW1(:,5)
!
    ZZW1(:,6) = ZZW1(:,5) * XMNU0                                     ! RCHMSI
    PRIS1D(:) = PRIS1D(:) + ZZW1(:,6)
    PRSS1D(:) = PRSS1D(:) - ZZW1(:,6)
  END WHERE
  DEALLOCATE(IVEC2)
  DEALLOCATE(ZVEC2)
  DEALLOCATE(ZVEC1)
  !
  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'HMS', &
                                          Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'HMS', &
                                          Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HMS', &
                                          Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if
END IF
!
!
!*      1.3  Ice multiplication process due to ice-ice collisions
!       ---------------------------------------------------------
!
GCIBU(:) = LCIBU .AND. (PRST1D(:)>XRTMIN(5)) .AND. (PRGT1D(:)>XRTMIN(6))
ICIBU    = COUNT( GCIBU(:) )
!
IF (ICIBU > 0) THEN
!
  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'CIBU', &
                                            Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'CIBU', &
                                            Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CIBU', &
                                            Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if
!
!       1.3.0 randomization of XNDEBRIS_CIBU values
!
  IF (GFIRSTCALL) THEN
    CALL RANDOM_SEED(SIZE=NI_SEED) ! get size of seed
    ALLOCATE(I_SEED(NI_SEED))
    I_SEED(:) = I_SEED_PARAM !
    CALL RANDOM_SEED(PUT=I_SEED)
    GFIRSTCALL = .FALSE.
  END IF
!
  ALLOCATE(ZFRAGMENTS(ICIBU))
!
  IF (XNDEBRIS_CIBU >= 0.0) THEN
    ZFRAGMENTS(:) = XNDEBRIS_CIBU
  ELSE
!
! Mantissa gives the mean value (randomization around 10**MANTISSA)
! First digit after the comma provides the full range around 10**MANTISSA
!
    ALLOCATE(ZHARVEST(ICIBU))
!
    ZFACT1_XNDEBRIS = AINT(XNDEBRIS_CIBU)
    ZFACT2_XNDEBRIS = ABS(ANINT(10.0*(XNDEBRIS_CIBU - ZFACT1_XNDEBRIS)))
!
    CALL RANDOM_NUMBER(ZHARVEST(:))
!
    ZFRAGMENTS(:) = 10.0**(ZFACT2_XNDEBRIS*ZHARVEST(:) + ZFACT1_XNDEBRIS)
!
    DEALLOCATE(ZHARVEST)
!
! ZFRAGMENTS is a random variable containing the number of fragments per collision
! For XNDEBRIS_CIBU=-1.2345  => ZFRAGMENTS(:) = 10.0**(2.0*RANDOM_NUMBER(ZHARVEST(:)) - 1.0)
! and ZFRAGMENTS=[0.1, 10.0] centered around 1.0
!
  END IF
!
!
!       1.3.1 To compute the partial integration of snow gamma function
!
!       1.3.1.0 allocations
!
  ALLOCATE(ZVEC1_S(ICIBU))
  ALLOCATE(ZVEC1_SW(ICIBU))
  ALLOCATE(ZVEC1_S1(ICIBU))
  ALLOCATE(ZVEC1_S2(ICIBU))
  ALLOCATE(ZVEC1_S3(ICIBU))
  ALLOCATE(ZVEC1_S4(ICIBU))
  ALLOCATE(ZVEC1_S11(ICIBU))
  ALLOCATE(ZVEC1_S12(ICIBU))
  ALLOCATE(ZVEC1_S21(ICIBU))
  ALLOCATE(ZVEC1_S22(ICIBU))
  ALLOCATE(ZVEC1_S31(ICIBU))
  ALLOCATE(ZVEC1_S32(ICIBU))
  ALLOCATE(ZVEC1_S41(ICIBU))
  ALLOCATE(ZVEC1_S42(ICIBU))
  ALLOCATE(ZVEC2_S1(ICIBU))
  ALLOCATE(IVEC2_S1(ICIBU))
  ALLOCATE(ZVEC2_S2(ICIBU))
  ALLOCATE(IVEC2_S2(ICIBU))
!
!
!       1.3.1.1 select the PLBDAS
!
  ZVEC1_S(:) = PACK( PLBDAS(:),MASK=GCIBU(:) )
  ZVEC1_SW(:)= ( XFVELOS**XALPHAS + ZVEC1_S(:)**XALPHAS ) ** (1./XALPHAS) ! modified equivalent lambda 
!
!
!       1.3.1.2 find the next lower indice for the PLBDAS in the
!               geometrical set of Lbda_s used to tabulate some moments of the
!               incomplete gamma function, for boundary 1 (0.2 mm)
!
  ZVEC2_S1(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_S  &
                      * LOG( ZVEC1_S(1:ICIBU) ) + XCIBUINTP1_S  ) )
  IVEC2_S1(1:ICIBU) = INT( ZVEC2_S1(1:ICIBU) )
  ZVEC2_S1(1:ICIBU) = ZVEC2_S1(1:ICIBU) - FLOAT( IVEC2_S1(1:ICIBU) )
!
!
!       1.3.1.3 find the next lower indice for the PLBDAS in the
!               geometrical set of Lbda_s used to tabulate some moments of the
!               incomplete gamma function, for boundary 2 (1 mm)
!
  ZVEC2_S2(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_S  &
                      * LOG( ZVEC1_S(1:ICIBU) ) + XCIBUINTP2_S  ) )
  IVEC2_S2(1:ICIBU) = INT( ZVEC2_S2(1:ICIBU) )
  ZVEC2_S2(1:ICIBU) = ZVEC2_S2(1:ICIBU) - FLOAT( IVEC2_S2(1:ICIBU) )
!
!
!       1.3.1.4 perform the linear interpolation of the
!               normalized "0"-moment of the incomplete gamma function
!
! For lower boundary (0.2 mm)
  ZVEC1_S11(1:ICIBU) = XGAMINC_CIBU_S(1,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                     - XGAMINC_CIBU_S(1,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
  ZVEC1_S12(1:ICIBU) = XGAMINC_CIBU_S(1,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                     - XGAMINC_CIBU_S(1,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! Computation of spectrum from 0.2 mm to 1 mm
  ZVEC1_S1(1:ICIBU) = ZVEC1_S12(1:ICIBU) - ZVEC1_S11(1:ICIBU)
!
!
!       1.3.1.5 perform the linear interpolation of the
!               normalized "XBS"-moment of the incomplete gamma function
!
! For lower boundary (0.2 mm)
  ZVEC1_S31(1:ICIBU) = XGAMINC_CIBU_S(3,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                     - XGAMINC_CIBU_S(3,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
  ZVEC1_S32(1:ICIBU) = XGAMINC_CIBU_S(3,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                     - XGAMINC_CIBU_S(3,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! From 0.2 mm to 1 mm we need
  ZVEC1_S3(1:ICIBU) = XMOMGS_CIBU_2 * (ZVEC1_S32(1:ICIBU) - ZVEC1_S31(1:ICIBU))
!
!
!       1.3.1.2 find the next lower indice for the PLBDS in the
!               geometrical set of Lbda_s used to tabulate some moments of the
!               incomplete gamma function, for boundary 1 (0.2 mm) for modified lambda (Wurtz snow fall speed)
!
  ZVEC2_S1(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_S  &
                      * LOG( ZVEC1_SW(1:ICIBU) ) + XCIBUINTP1_S  ) )
  IVEC2_S1(1:ICIBU) = INT( ZVEC2_S1(1:ICIBU) )
  ZVEC2_S1(1:ICIBU) = ZVEC2_S1(1:ICIBU) - FLOAT( IVEC2_S1(1:ICIBU) )
!
!
!       1.3.1.3 find the next lower indice for the PLBDS in the
!               geometrical set of Lbda_s used to tabulate some moments of the
!               incomplete gamma function, for boundary 2 (1 mm) for modified lambda (Wurtz snow fall speed)
!
  ZVEC2_S2(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_S  &
                      * LOG( ZVEC1_SW(1:ICIBU) ) + XCIBUINTP2_S  ) )
  IVEC2_S2(1:ICIBU) = INT( ZVEC2_S2(1:ICIBU) )
  ZVEC2_S2(1:ICIBU) = ZVEC2_S2(1:ICIBU) - FLOAT( IVEC2_S2(1:ICIBU) )
!
!
!       1.3.1.5 perform the linear interpolation of the
!               normalized "XDS"-moment of the incomplete gamma function
!
! For lower boundary (0.2 mm)
  ZVEC1_S21(1:ICIBU) = XGAMINC_CIBU_S(2,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                     - XGAMINC_CIBU_S(2,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
  ZVEC1_S22(1:ICIBU) = XGAMINC_CIBU_S(2,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                     - XGAMINC_CIBU_S(2,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! From 0.2 mm to 1 mm we need
  ZVEC1_S2(1:ICIBU) = XMOMGS_CIBU_1 * (ZVEC1_S22(1:ICIBU) - ZVEC1_S21(1:ICIBU))
!
!
!       1.3.1.6 perform the linear interpolation of the
!               normalized "XBS+XDS"-moment of the incomplete gamma function
!
! For lower boundary (0.2 mm)
  ZVEC1_S41(1:ICIBU) = XGAMINC_CIBU_S(4,IVEC2_S1(1:ICIBU)+1) *  ZVEC2_S1(1:ICIBU) &
                     - XGAMINC_CIBU_S(4,IVEC2_S1(1:ICIBU))   * (ZVEC2_S1(1:ICIBU)-1.0)
!
! For upper boundary (1 mm)
  ZVEC1_S42(1:ICIBU) = XGAMINC_CIBU_S(4,IVEC2_S2(1:ICIBU)+1) *  ZVEC2_S2(1:ICIBU) &
                     - XGAMINC_CIBU_S(4,IVEC2_S2(1:ICIBU))   * (ZVEC2_S2(1:ICIBU)-1.0)
!
! From 0.2 mm to 1 mm we need
  ZVEC1_S4(1:ICIBU) = XMOMGS_CIBU_3 * (ZVEC1_S42(1:ICIBU) - ZVEC1_S41(1:ICIBU))
!
  ALLOCATE(ZINTG_SNOW_1(SIZE(PZT)))
  ALLOCATE(ZINTG_SNOW_2(SIZE(PZT)))
  ALLOCATE(ZINTG_SNOW_3(SIZE(PZT)))
  ALLOCATE(ZINTG_SNOW_4(SIZE(PZT)))
!
  ZINTG_SNOW_1(:) = UNPACK ( VECTOR=ZVEC1_S1(:),MASK=GCIBU,FIELD=0.0 )
  ZINTG_SNOW_2(:) = UNPACK ( VECTOR=ZVEC1_S2(:),MASK=GCIBU,FIELD=0.0 )
  ZINTG_SNOW_3(:) = UNPACK ( VECTOR=ZVEC1_S3(:),MASK=GCIBU,FIELD=0.0 )
  ZINTG_SNOW_4(:) = UNPACK ( VECTOR=ZVEC1_S4(:),MASK=GCIBU,FIELD=0.0 )
!
!
!       1.3.2 Compute the partial integration of graupel gamma function
!
!       1.3.2.0 allocations
!
  ALLOCATE(ZVEC1_G(ICIBU))
  ALLOCATE(ZVEC1_G1(ICIBU))
  ALLOCATE(ZVEC1_G2(ICIBU))
  ALLOCATE(ZVEC2_G(ICIBU))
  ALLOCATE(IVEC2_G(ICIBU))
!
!
!       1.3.2.1 select the PLBDAG
!
  ZVEC1_G(:) = PACK( PLBDAG(:),MASK=GCIBU(:) )
!
!
!       1.3.2.2 find the next lower indice for the PLBDAG in the
!               geometrical set of Lbda_g used to tabulate some moments of the
!               incomplete gamma function, for the "2mm" boundary
!
  ZVEC2_G(1:ICIBU) = MAX( 1.0001, MIN( FLOAT(NGAMINC)-0.0001,XCIBUINTP_G  &
                     * LOG( ZVEC1_G(1:ICIBU) ) + XCIBUINTP1_G  ) )
  IVEC2_G(1:ICIBU) = INT( ZVEC2_G(1:ICIBU) )
  ZVEC2_G(1:ICIBU) = ZVEC2_G(1:ICIBU) - FLOAT( IVEC2_G(1:ICIBU) )
!
!
!       1.3.2.3 perform the linear interpolation of the
!               normalized "2+XDG"-moment of the incomplete gamma function
!
  ZVEC1_G1(1:ICIBU) = XGAMINC_CIBU_G(1,IVEC2_G(1:ICIBU)+1) *  ZVEC2_G(1:ICIBU)    &
                    - XGAMINC_CIBU_G(1,IVEC2_G(1:ICIBU))   * (ZVEC2_G(1:ICIBU)-1.0)
!
! From 2 mm to infinity we need
  ZVEC1_G1(1:ICIBU) = XMOMGG_CIBU_1 * (1.0 - ZVEC1_G1(1:ICIBU))
!
!
!       1.3.2.4 perform the linear interpolation of the
!               normalized "2.0"-moment of the incomplete gamma function
!
  ZVEC1_G2(1:ICIBU) = XGAMINC_CIBU_G(2,IVEC2_G(1:ICIBU)+1) *  ZVEC2_G(1:ICIBU)    &
                    - XGAMINC_CIBU_G(2,IVEC2_G(1:ICIBU))   * (ZVEC2_G(1:ICIBU)-1.0)
!
! From 2 mm to infinity we need
  ZVEC1_G2(1:ICIBU) = XMOMGG_CIBU_2 * (1.0 - ZVEC1_G2(1:ICIBU))
!
!
  ALLOCATE(ZINTG_GRAUPEL_1(SIZE(PZT)))
  ALLOCATE(ZINTG_GRAUPEL_2(SIZE(PZT)))
!
  ZINTG_GRAUPEL_1(:) = UNPACK ( VECTOR=ZVEC1_G1(:),MASK=GCIBU,FIELD=0.0 )
  ZINTG_GRAUPEL_2(:) = UNPACK ( VECTOR=ZVEC1_G2(:),MASK=GCIBU,FIELD=0.0 )
!
!
!        1.3.3 To compute final "CIBU" contributions
!
  ALLOCATE(ZNI_CIBU(SIZE(PZT)))
  ALLOCATE(ZFRAG_CIBU(SIZE(PZT)))
!
  ZFRAG_CIBU(:) = UNPACK ( VECTOR=ZFRAGMENTS(:),MASK=GCIBU,FIELD=0.0 )
  ZNI_CIBU(:) = ZFRAG_CIBU(:) * (XFACTOR_CIBU_NI * PCST1D(:) * PCGT1D(:) / (PRHODREF(:)**XCEXVT)) * &
                (XCG * ZINTG_GRAUPEL_1(:) * ZINTG_SNOW_1(:) *                                               &
                 PLBDAG(:)**(-(XDG+2.0))                                             &
               - XCS * ZINTG_GRAUPEL_2(:) * ZINTG_SNOW_2(:) *                                               &
                 PLBDAS(:)**(-XDS) * PLBDAG(:)**(-2.0) *                                            &
                 (1+(XFVELOS/PLBDAS(:))**XALPHAS)**(-XNUS-XDS/XALPHAS) )

  PCIS1D(:) = PCIS1D(:) + MAX(ZNI_CIBU(:), 0.)
!
  DEALLOCATE(ZFRAG_CIBU)
  DEALLOCATE(ZFRAGMENTS)
!
! Max value of rs removed by CIBU
  ALLOCATE(ZRI_CIBU(SIZE(PZT)))
  ZRI_CIBU(:) = (XFACTOR_CIBU_RI * PCST1D(:) * PCGT1D(:) / (PRHODREF(:)**XCEXVT)) * &
                 (XCG * ZINTG_GRAUPEL_1(:) * ZINTG_SNOW_3(:) *                              &
                  PLBDAS(:)**(-XBS) * PLBDAG(:)**(-(XDG+2.0))                                               &
                - XCS * ZINTG_GRAUPEL_2(:) * ZINTG_SNOW_4(:) *                              &
                  PLBDAS(:)**(-XBS-XDS) * PLBDAG(:)**(-2.0) *                               &
                  (1+(XFVELOS/PLBDAS(:))**XALPHAS)**(-XNUS-(XBS+XDS)/XALPHAS) )
!
! The value of rs removed by CIBU is determined by the mean mass of pristine ice
  WHERE( PRIT1D(:)>XRTMIN(4) .AND. PCIT1D(:)>XCTMIN(4) )
    ZRI_CIBU(:) = MIN( ZRI_CIBU(:), PRSS1D(:), ZNI_CIBU(:)*PRIT1D(:)/PCIT1D(:) )
  ELSEWHERE
    ZRI_CIBU(:) = MIN( ZRI_CIBU(:), PRSS1D(:), MAX( ZNI_CIBU(:)*XMNU0,XRTMIN(4) ) )
  END WHERE
!
  PRIS1D(:) = PRIS1D(:) + MAX(ZRI_CIBU(:), 0.)   !
  PRSS1D(:) = PRSS1D(:) - MAX(ZRI_CIBU(:), 0.)   !
!
  DEALLOCATE(ZVEC1_S)
  DEALLOCATE(ZVEC1_SW)
  DEALLOCATE(ZVEC1_S1)
  DEALLOCATE(ZVEC1_S2)
  DEALLOCATE(ZVEC1_S3)
  DEALLOCATE(ZVEC1_S4)
  DEALLOCATE(ZVEC1_S11)
  DEALLOCATE(ZVEC1_S12)
  DEALLOCATE(ZVEC1_S21)
  DEALLOCATE(ZVEC1_S22)
  DEALLOCATE(ZVEC1_S31)
  DEALLOCATE(ZVEC1_S32)
  DEALLOCATE(ZVEC1_S41)
  DEALLOCATE(ZVEC1_S42)
  DEALLOCATE(ZVEC2_S1)
  DEALLOCATE(IVEC2_S1)
  DEALLOCATE(ZVEC2_S2)
  DEALLOCATE(IVEC2_S2)
  DEALLOCATE(ZVEC1_G)
  DEALLOCATE(ZVEC1_G1)
  DEALLOCATE(ZVEC1_G2)
  DEALLOCATE(ZVEC2_G)
  DEALLOCATE(IVEC2_G)
  DEALLOCATE(ZINTG_SNOW_1)
  DEALLOCATE(ZINTG_SNOW_2)
  DEALLOCATE(ZINTG_SNOW_3)
  DEALLOCATE(ZINTG_SNOW_4)
  DEALLOCATE(ZINTG_GRAUPEL_1)
  DEALLOCATE(ZINTG_GRAUPEL_2)
  DEALLOCATE(ZNI_CIBU)
  DEALLOCATE(ZRI_CIBU)
  !
  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'CIBU', &
                                          Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'CIBU', &
                                          Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CIBU', &
                                          Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if
END IF
!
!
!*       1.4  Rain accretion onto the aggregates  
!        ---------------------------------------
!
!
ZZW1(:,2:3) = 0.0
GACC(:) = (PRRT1D(:)>XRTMIN(3)) .AND. (PRST1D(:)>XRTMIN(5)) .AND. (PRRS1D(:)>XRTMIN(3)/PTSTEP) .AND. (PZT(:)<XTT)
IGACC = COUNT( GACC(:) )
!
IF( IGACC>0 .AND. LRAIN) THEN
  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'ACC', &
                                            Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'ACC', &
                                            Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'ACC', &
                                            Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'ACC', &
                                            Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'ACC', &
                                            Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv .and. NMOM_S.GE.2) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'ACC', &
                                            Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv .and. NMOM_G.GE.2) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'ACC', &
                                            Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) )
  end if
!
!        1.4.0  allocations
!
  ALLOCATE(ZVEC1(IGACC))
  ALLOCATE(ZVEC2(IGACC))
  ALLOCATE(ZVEC3(IGACC))
  ALLOCATE(IVEC1(IGACC))
  ALLOCATE(IVEC2(IGACC))
!
!        1.4.1  select the (PLBDAS,PLBDAR) couplet
!
  if (M2_ICE) then
     ZVEC1(:) = PACK( MAX(MIN(PLBDAS(:),5.E5),5.E1),MASK=GACC(:) )
  else
     ZVEC1(:) = PACK( PLBDAS(:),MASK=GACC(:) )
  end if
  ZVEC2(:) = PACK( PLBDAR(:),MASK=GACC(:) )
!
!        1.4.2  find the next lower indice for the PLBDAS and for the PLBDAR
!               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
!               tabulate the RACCSS-kernel
!
  ZVEC1(1:IGACC) = MAX( 1.0001, MIN( REAL(NACCLBDAS)-0.0001,           &
                        XACCINTP1S * LOG( ZVEC1(1:IGACC) ) + XACCINTP2S ) )
  IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
  ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - REAL( IVEC1(1:IGACC) )
!
  ZVEC2(1:IGACC) = MAX( 1.0001, MIN( REAL(NACCLBDAR)-0.0001,           &
                        XACCINTP1R * LOG( ZVEC2(1:IGACC) ) + XACCINTP2R ) )
  IVEC2(1:IGACC) = INT( ZVEC2(1:IGACC) )
  ZVEC2(1:IGACC) = ZVEC2(1:IGACC) - REAL( IVEC2(1:IGACC) )
!
!        1.4.3  perform the bilinear interpolation of the normalized
!               RACCSS-kernel
!
  DO JJ = 1,IGACC
     ZVEC3(JJ) =  (  XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                   - XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                          *  ZVEC1(JJ)          &
                - (  XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                   - XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
  END DO
  ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
  DEALLOCATE(ZVEC3)
  ALLOCATE(ZVEC3(IGACC))  
!
!                                                                         
!        1.4.3b  perform the bilinear interpolation of the normalized
!               RACCSS-kernel FOR CONCENTRATION
!
  DO JJ = 1,IGACC
     ZVEC3(JJ) =  (  XKER_N_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_N_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                				 	           * ZVEC1(JJ) &
                 - (  XKER_N_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_N_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
  	                    			             * (ZVEC1(JJ) - 1.0)
  END DO
  ZZNW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
!        1.4.4  raindrop accretion on the small sized aggregates
!
  WHERE ( GACC(:) )
     ZZW1(:,2) = PCRT1D(:) *                                      & !! coef of RRACCS 
                  XFRACCSS*( PCST1D(:) )*( PRHODREF(:)**(-XCEXVT+1.) ) &
             *( XLBRACCS1/((PLBDAS(:)**2)               ) +                  &
                XLBRACCS2/( PLBDAS(:)    * PLBDAR(:)    ) +                  &
                XLBRACCS3/(               (PLBDAR(:)**2)) )/PLBDAR(:)**3
!                                                                                
     ZZNW1(:,2) = PCRT1D(:) *                                           & !! coef of CRACCS
                  XFNRACCSS*( PCST1D(:) )*( PRHODREF(:)**(-XCEXVT+1.) ) &
             *( XLBNRACCS1/((PLBDAS(:)**2)               ) +                  &
               XLBNRACCS2/( PLBDAS(:)    * PLBDAR(:)    ) +                  &
                XLBNRACCS3/(               (PLBDAR(:)**2)) )
     ZZW1(:,4) = MIN( PRRS1D(:),ZZW1(:,2)*ZZW(:) )           ! RRACCSS
     PRRS1D(:) = PRRS1D(:) - ZZW1(:,4)
     PRSS1D(:) = PRSS1D(:) + ZZW1(:,4)
     PTHS1D(:) = PTHS1D(:) + ZZW1(:,4)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(RRACCSS))
!
     ZZNW1(:,4) = MIN( PCRS1D(:),ZZNW1(:,2)*ZZNW(:) )           ! NRACCSS      
     PCRS1D(:) = PCRS1D(:)-ZZNW1(:,4)
  END WHERE
!
!        1.4.4b perform the bilinear interpolation of the normalized
!               RACCS-kernel
!
  DO JJ = 1,IGACC
    ZVEC3(JJ) =  (   XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                  -  XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                         *  ZVEC1(JJ)          &
               - (   XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                  -  XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
  END DO
  ZZW1(:,2) = ZZW1(:,2)*UNPACK( VECTOR=ZVEC3(:),MASK=GACC(:),FIELD=0.0 ) !! RRACCS
!                                                                                 
!        1.3.4b2 perform the bilinear interpolation of the normalized
!               RACCS-kernel
!
  DO JJ = 1,IGACC
     ZVEC3(JJ) =  (   XKER_N_RACCS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                         -  XKER_N_RACCS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                         * ZVEC1(JJ) &
                         - (   XKER_N_RACCS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                         -  XKER_N_RACCS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                           * (ZVEC1(JJ) - 1.0)
  END DO
  ZZNW1(:,2) = ZZNW1(:,2)*UNPACK( VECTOR=ZVEC3(:),MASK=GACC(:),FIELD=0.0 ) !! NRACCS
!
!        1.4.5  perform the bilinear interpolation of the normalized
!               SACCRG-kernel
!
  DO JJ = 1,IGACC
    ZVEC3(JJ) =  (  XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)+1) *  ZVEC1(JJ)          &
                  - XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)  ) * (ZVEC1(JJ) - 1.0) ) &
                                                         *  ZVEC2(JJ)          &
               - (  XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)+1) *  ZVEC1(JJ)          &
                  - XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)  ) * (ZVEC1(JJ) - 1.0) ) &
                                                         * (ZVEC2(JJ) - 1.0)
  END DO
  ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
!        1.3.5b  perform the bilinear interpolation of the normalized
!               SACCRG-kernel
!
  DO JJ = 1,IGACC
     ZVEC3(JJ) =  (  XKER_N_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                        - XKER_N_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
      			 	                                   * ZVEC2(JJ) &
                     - (  XKER_N_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                        - XKER_N_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
			                                     * (ZVEC2(JJ) - 1.0)
  END DO
  ZZNW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
!        1.4.6  raindrop accretion-conversion of the large sized aggregates
!               into graupeln
!
  WHERE ( GACC(:) .AND. (PRSS1D(:)>XRTMIN(5)/PTSTEP)  .AND. (PCSS1D(:)>XCTMIN(5)/PTSTEP) )
     ZZW1(:,2) = MAX( MIN( PRRS1D(:),ZZW1(:,2)-ZZW1(:,4) ) , 0. )      ! RRACCSG
     ZZNW1(:,2) = MAX( MIN( PCRS1D(:),ZZNW1(:,2)-ZZNW1(:,4) ) , 0. )   ! NRACCSG  
     ZZW1(:,3) = MIN( PRSS1D(:),PCRT1D(:)*XFSACCRG*ZZW(:)* PCST1D(:) *           & ! RSACCRG 
                PLBDAS(:)**(-XBS) * ( PRHODREF(:)**(-XCEXVT+1.) )     & 
                *( XLBSACCR1/((PLBDAR(:)**2)               ) +           &
                  XLBSACCR2/( PLBDAR(:)    * PLBDAS(:)    ) +           &
                  XLBSACCR3/(               (PLBDAS(:)**2)) ) )
     ZZNW1(:,3) = MIN( PCSS1D(:),PCRT1D(:)*XFNSACCRG*ZZNW(:)* PCST1D(:) *        & ! NSACCRG 
                                      ( PRHODREF(:)**(-XCEXVT+1.) )     &            
               *( XLBNSACCR1/((PLBDAR(:)**2)               ) +          &            
                  XLBNSACCR2/( PLBDAR(:)    * PLBDAS(:)    ) +          &            
                  XLBNSACCR3/(               (PLBDAS(:)**2)) ) )           
     PRRS1D(:) = PRRS1D(:) - ZZW1(:,2)
     PRSS1D(:) = PRSS1D(:) - ZZW1(:,3)
     PRGS1D(:) = PRGS1D(:) + ZZW1(:,2)+ZZW1(:,3)
     PTHS1D(:) = PTHS1D(:) + ZZW1(:,2)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(RRACCSG))
!
     PCRS1D(:) = PCRS1D(:)-ZZNW1(:,2) !                                   
     PCSS1D(:) = PCSS1D(:)-ZZNW1(:,3)
     PCGS1D(:) = PCGS1D(:)+ZZNW1(:,3)
  END WHERE
   DEALLOCATE(IVEC2)
   DEALLOCATE(IVEC1)
   DEALLOCATE(ZVEC3)
   DEALLOCATE(ZVEC2)
   DEALLOCATE(ZVEC1)
  !
  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'ACC', &
                                           Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'ACC', &
                                           Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'ACC', &
                                           Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'ACC', &
                                           Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'ACC', &
                                           Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv .and. NMOM_S.GE.2) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'ACC', &
                                           Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv .and. NMOM_G.GE.2) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'ACC', &
                                           Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) )
  end if
END IF
!
!
!*       1.5  Conversion-Melting of the aggregates
!        -----------------------------------------
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'CMEL', &
                                          Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'CMEL', &
                                          Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv .and. NMOM_S.GE.2 .and. NMOM_G.GE.2 ) then
                   call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'CMEL', &
                                          Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'CMEL', &
                                          Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) ) 
  end if
end if
!
ZZW(:) = 0.0
WHERE( (PRST1D(:)>XRTMIN(5)) .AND. (PRSS1D(:)>XRTMIN(5)/PTSTEP) .AND. (PZT(:)>XTT) )
   ZZW(:) = PRVT1D(:)*PPRES(:)/((XMV/XMD)+PRVT1D(:)) ! Vapor pressure
   ZZW(:) =  PKA(:)*(XTT-PZT(:)) +                                 &
              ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT )) &
                          *(XESTT-ZZW(:))/(XRV*PZT(:))             )
!
! compute RSMLT
!
   ZZW(:)  = MIN( PRSS1D(:), XFSCVMG*MAX( 0.0,( -ZZW(:) * PCST1D(:) *   & 
                          ( X0DEPS*       PLBDAS(:)**XEX0DEPS +     &
                            X1DEPS*PCJ(:)*PLBDAS(:)**XEX1DEPS *     &
                                  (1+0.5*(XFVELOS/PLBDAS(:))**XALPHAS)**(-XNUS+XEX1DEPS/XALPHAS) ) -   &
                                    ( ZZW1(:,1)+ZZW1(:,4) ) *       &
                             ( PRHODREF(:)*XCL*(XTT-PZT(:))) ) /    &
                                            ( PRHODREF(:)*XLMTT ) ) )
   PRSS1D(:) = PRSS1D(:) - ZZW(:)
   PRGS1D(:) = PRGS1D(:) + ZZW(:)
   PCSS1D(:) = MAX( PCSS1D(:) - ZZW(:)*(MAX(PCST1D(:),XCTMIN(5))/MAX(PRST1D(:),XRTMIN(5))), 0.0 )
   PCGS1D(:) = MAX( PCGS1D(:) + ZZW(:)*(MAX(PCST1D(:),XCTMIN(5))/MAX(PRST1D(:),XRTMIN(5))), 0.0 )
END WHERE
!
! Budget storage
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'CMEL', &
                                         Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'CMEL', &
                                         Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv .and. NMOM_S.GE.2 .AND. NMOM_G.GE.2) then
                   call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'CMEL', &
                                         Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
                   call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'CMEL', &
                                         Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) ) 
  end if
end if

END IF SNOW
!
!------------------------------------------------------------------------------
!
!                         #################
!                         FAST RG PROCESSES
!                         #################
!
!
!*       2.1  Rain contact freezing  
!        --------------------------
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'CFRZ', &
                                         Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'CFRZ', &
                                         Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'CFRZ', &
                                         Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'CFRZ', &
                                         Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'CFRZ', &
                                         Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CFRZ', &
                                         Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv .and. NMOM_G.GE.2 ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'CFRZ', &
                                         Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) )
end if

ZZW1(:,3:4) = 0.0
WHERE( (PRIT1D(:)>XRTMIN(4)) .AND. (PRRT1D(:)>XRTMIN(3)) .AND. (PRIS1D(:)>XRTMIN(4)/PTSTEP) .AND. (PRRS1D(:)>XRTMIN(3)/PTSTEP) )
   ZZW1(:,3) = MIN( PRIS1D(:),XICFRR * PRIT1D(:) * PCRT1D(:)          & ! RICFRRG
                                   * PLBDAR(:)**XEXICFRR        &
                                   * PRHODREF(:)**(-XCEXVT-1.0) )
!
   ZZW1(:,4) = MIN( PRRS1D(:),XRCFRI * PCIT1D(:) * PCRT1D(:)          & ! RRCFRIG
                                   * PLBDAR(:)**XEXRCFRI        &
                                   * PRHODREF(:)**(-XCEXVT-2.0) )
   PRIS1D(:) = PRIS1D(:) - ZZW1(:,3)
   PRRS1D(:) = PRRS1D(:) - ZZW1(:,4)
   PRGS1D(:) = PRGS1D(:) + ZZW1(:,3)+ZZW1(:,4)
   PTHS1D(:) = PTHS1D(:) + ZZW1(:,4)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*RRCFRIG)
!
   PCIS1D(:) = MAX( PCIS1D(:)-ZZW1(:,3)*(PCIT1D(:)/PRIT1D(:)),0.0 )     ! CICFRRG
   PCRS1D(:) = MAX( PCRS1D(:)-ZZW1(:,4)*(PCRT1D(:)/PRRT1D(:)),0.0 )     ! CRCFRIG
   PCGS1D(:) = PCGS1D(:)+ZZW1(:,3)*(PCIT1D(:)/PRIT1D(:)) 
END WHERE
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'CFRZ', &
                                         Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'CFRZ', &
                                         Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'CFRZ', &
                                         Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'CFRZ', &
                                         Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'CFRZ', &
                                         Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CFRZ', &
                                         Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv .and. NMOM_G.GE.2 ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'CFRZ', &
                                         Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) )
end if
!
!
!*       2.2  Ice multiplication process following rain contact freezing
!        ---------------------------------------------------------------
!
GRDSF(:) = LRDSF .AND. (PRIT1D(:)>0.0) .AND. (PRRT1D(:)>0.0) .AND. &
                       (PRIS1D(:)>0.0) .AND. (PRRS1D(:)>0.0)
IRDSF    = COUNT( GRDSF(:) )
!
IF (IRDSF > 0) THEN
!
  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'RDSF', &
                                            Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'RDSF', &
                                            Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'RDSF', &
                                            Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if
!
  ALLOCATE(ZVEC1_R(IRDSF))
  ALLOCATE(ZVEC1_R1(IRDSF))
  ALLOCATE(ZVEC2_R(IRDSF))
  ALLOCATE(IVEC2_R(IRDSF))
!
!*       2.2.1  select the ZLBDAR
!
  ZVEC1_R(:) = PACK( PLBDAR(:),MASK=GRDSF(:) )
!
!*       2.2.2  find the next lower indice for the ZLBDAR in the
!               geometrical set of Lbda_r used to tabulate some moments of the
!               incomplete gamma function, for the lower boundary (0.1 mm)
!
  ZVEC2_R(1:IRDSF) = MAX( 1.00001, MIN( FLOAT(NGAMINC)-0.00001,XRDSFINTP_R  &
                        * LOG( ZVEC1_R(1:IRDSF) ) + XRDSFINTP1_R  ) )
  IVEC2_R(1:IRDSF) = INT( ZVEC2_R(1:IRDSF) )
  ZVEC2_R(1:IRDSF) = ZVEC2_R(1:IRDSF) - FLOAT( IVEC2_R(1:IRDSF) )
!
!*       2.2.3  perform the linear interpolation of the
!               normalized "2+XDR"-moment of the incomplete gamma function
!
  ZVEC1_R1(1:IRDSF) = XGAMINC_RDSF_R(IVEC2_R(1:IRDSF)+1) *  ZVEC2_R(1:IRDSF)    &
                    - XGAMINC_RDSF_R(IVEC2_R(1:IRDSF))   * (ZVEC2_R(1:IRDSF) - 1.0)
!
!  From 0.1 mm to infinity we need
  ZVEC1_R1(1:IRDSF) = XMOMGR_RDSF * (1.0 - ZVEC1_R1(1:IRDSF))
!
  ALLOCATE(ZINTG_RAIN(SIZE(PZT)))
  ZINTG_RAIN(:) = UNPACK ( VECTOR=ZVEC1_R1(:),MASK=GRDSF,FIELD=0.0 )
!
!*       2.2.4  To compute final "RDSF" contributions
!
  ALLOCATE(ZNI_RDSF(SIZE(PZT)))
  ZNI_RDSF(:) = (XFACTOR_RDSF_NI / (PRHODREF(:)**(XCEXVT-1.0))) * (  &
                 PCIT1D(:) * PCRT1D(:) * ZINTG_RAIN(:) * PLBDAR(:)**(-(XDR+6.0)) )
!
  PCIS1D(:) = PCIS1D(:) + ZNI_RDSF(:)
!
! The value of rg removed by RDSF is determined by the mean mass of pristine ice
  ALLOCATE(ZRI_RDSF(SIZE(PZT)))
  ZRI_RDSF(:) = MIN( PRGS1D(:), MAX( ZNI_RDSF(:)*XMNU0,XRTMIN(5) ) )
!
  PRIS1D(:) = PRIS1D(:) + ZRI_RDSF(:)
  PRGS1D(:) = PRGS1D(:) - ZRI_RDSF(:)
!
  DEALLOCATE(ZINTG_RAIN)
  DEALLOCATE(ZVEC1_R)
  DEALLOCATE(ZVEC1_R1)
  DEALLOCATE(ZVEC2_R)
  DEALLOCATE(IVEC2_R)
  DEALLOCATE(ZNI_RDSF)
  DEALLOCATE(ZRI_RDSF)
  !
  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'RDSF', &
                                            Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RG), 'RDSF', &
                                            Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'RDSF', &
                                            Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if
ENDIF
!
!
!*       2.3  Compute the Dry growth case
!        --------------------------------
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'WETG', &
                                          Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'WETG', &
                                          Unpack( prcs1d(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'WETG', &
                                          Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'WETG', &
                                          Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'WETG', &
                                          Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'WETG', &
                                          Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'WETG', &
                                          Unpack( prhs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'WETG', &
                                          Unpack( pccs1d(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'WETG', &
                                          Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'WETG', &
                                          Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
       if(M2_ICE) then
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'WETG', &
                                         Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'WETG', &
                                         Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) )
                  if (LHAIL) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nh), 'WETG', &
                                         Unpack( pchs1d(:), mask = gmicro(:, :, :), field = pchs(:, :, :) ) * prhodj(:, :, :) )
       end if
  end if
end if
!
ZZW1(:,:) = 0.0
ZZNW1(:,:) = 0.0   
WHERE( ((PRCT1D(:)>XRTMIN(2)) .AND. (PRGT1D(:)>XRTMIN(6)) .AND. (PRCS1D(:)>XRTMIN(2)/PTSTEP))  .OR. &
          ((PRIT1D(:)>XRTMIN(4)) .AND. (PRGT1D(:)>XRTMIN(6)) .AND. (PRIS1D(:)>XRTMIN(4)/PTSTEP))      )
   ZZW(:) = PLBDAG(:)**(-XDG-2.0) * PRHODREF(:)**(-XCEXVT) * PCGT1D(:) 
   ZZW1(:,1) = MIN( PRCS1D(:),XFCDRYG * PRCT1D(:) * ZZW(:) )             ! RCDRYG
   ZZNW1(:,1) =  MIN( PCCS1D(:),ZZW1(:,1) * PCCT1D(:) / MAX(PRCT1D(:),XRTMIN(2)) )
   ZZW1(:,2) = MIN( PRIS1D(:),XFIDRYG * EXP( XCOLEXIG*(PZT(:)-XTT) ) &
                                    * PRIT1D(:) * ZZW(:) )             ! RIDRYG
   ZZNW1(:,2) =  MIN( PCIS1D(:),ZZW1(:,2) * PCIT1D(:) / MAX(PRIT1D(:),XRTMIN(4)) )
END WHERE
!
!*       2.3.1  accretion of aggregates on the graupeln
!        ----------------------------------------------
!
GDRY(:) = (PRST1D(:)>XRTMIN(5)) .AND. (PRGT1D(:)>XRTMIN(6)) .AND. (PRSS1D(:)>XRTMIN(5)/PTSTEP)
if (M2_ICE) GDRY(:) = GDRY(:) .AND. (PCST1D(:)>XCTMIN(5)) .AND. PCSS1D(:)>XCTMIN(5)/PTSTEP
IGDRY = COUNT( GDRY(:) )
!
IF( IGDRY>0 ) THEN
!
!*       2.3.2  allocations
!
  ALLOCATE(ZVEC1(IGDRY))
  ALLOCATE(ZVEC2(IGDRY))
  ALLOCATE(ZVEC3(IGDRY))
  ALLOCATE(IVEC1(IGDRY))
  ALLOCATE(IVEC2(IGDRY))
!
!*       2.3.3  select the (PLBDAG,PLBDAS) couplet
!
  ZVEC1(:) = PACK( PLBDAG(:),MASK=GDRY(:) )
  ZVEC2(:) = PACK( PLBDAS(:),MASK=GDRY(:) )
!
!*       2.3.4  find the next lower indice for the PLBDAG and for the PLBDAS
!               in the geometrical set of (Lbda_g,Lbda_s) couplet use to
!               tabulate the SDRYG-kernel
!
  ZVEC1(1:IGDRY) = MAX( 1.0001, MIN( REAL(NDRYLBDAG)-0.0001,           &
                        XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + XDRYINTP2G ) )
  IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
  ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - REAL( IVEC1(1:IGDRY) )
!
  ZVEC2(1:IGDRY) = MAX( 1.0001, MIN( REAL(NDRYLBDAS)-0.0001,           &
                        XDRYINTP1S * LOG( ZVEC2(1:IGDRY) ) + XDRYINTP2S ) )
  IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
  ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
!
!*       2.3.5  perform the bilinear interpolation of the normalized
!               SDRYG-kernel
!
  DO JJ = 1,IGDRY
    ZVEC3(JJ) =  (  XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                  - XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                        *  ZVEC1(JJ)          &
               - (  XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                  - XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                        * (ZVEC1(JJ) - 1.0)
  END DO
  ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
  WHERE( GDRY(:) )
     ZZW1(:,3) = MIN( PRSS1D(:),XFSDRYG*ZZW(:)                         & ! RSDRYG
                                      * EXP( XCOLEXSG*(PZT(:)-XTT) )  &
                    *( PRST1D(:) )*( PCGT1D(:) )      &  
                    *( PRHODREF(:)**(-XCEXVT+1.) )                    &
                         *( XLBSDRYG1/( PLBDAG(:)**2              ) + &
                            XLBSDRYG2/( PLBDAG(:)   * PLBDAS(:)   ) + &
                            XLBSDRYG3/(               PLBDAS(:)**2) ) )
  END WHERE
!
!                                                                           
!*       2.2.5b  perform the bilinear interpolation of the normalized
!               SDRYG-kernel FOR CONCENTRATION
!
  DO JJ = 1,IGDRY
     ZVEC3(JJ) =  (  XKER_N_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                        - XKER_N_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
      	    		 	                                  * ZVEC1(JJ) &
                     - (  XKER_N_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                        - XKER_N_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                        * (ZVEC1(JJ) - 1.0)
  END DO
  ZZNW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
  WHERE( GDRY(:)   )
     ZZNW1(:,3) = MIN( PCSS1D(:),XFNSDRYG*ZZNW(:)                       & ! NSDRYG
                                   * EXP( XCOLEXSG*(PZT(:)-XTT) )  &
                 *( PCST1D(:)                     )*( PCGT1D(:) )      &  
                 *( PRHODREF(:)**(-XCEXVT+1.) )                    &
                      *( XLBNSDRYG1/( PLBDAG(:)**2             ) + &
                         XLBNSDRYG2/( PLBDAG(:)   * PLBDAS(:)  ) + &
                         XLBNSDRYG3/(               PLBDAS(:)**2)) )
!                                                                        
  END WHERE
  DEALLOCATE(IVEC2)
  DEALLOCATE(IVEC1)
  DEALLOCATE(ZVEC3)
  DEALLOCATE(ZVEC2)
  DEALLOCATE(ZVEC1)
END IF
!
!*       2.3.6  accretion of raindrops on the graupeln
!        ---------------------------------------------
!
GDRY(:) = (PRRT1D(:)>XRTMIN(3)) .AND. (PRGT1D(:)>XRTMIN(6)) .AND. (PRRS1D(:)>XRTMIN(3))
IGDRY = COUNT( GDRY(:) )
!
IF( IGDRY>0 ) THEN
!
!*       2.3.7  allocations
!
  ALLOCATE(ZVEC1(IGDRY))
  ALLOCATE(ZVEC2(IGDRY))
  ALLOCATE(ZVEC3(IGDRY))
  ALLOCATE(IVEC1(IGDRY))
  ALLOCATE(IVEC2(IGDRY))
!
!*       2.3.8  select the (PLBDAG,PLBDAR) couplet
!
  ZVEC1(:) = PACK( PLBDAG(:),MASK=GDRY(:) )
  ZVEC2(:) = PACK( PLBDAR(:),MASK=GDRY(:) )
!
!*       2.3.9  find the next lower indice for the PLBDAG and for the PLBDAR
!               in the geometrical set of (Lbda_g,Lbda_r) couplet use to
!               tabulate the RDRYG-kernel
!
  ZVEC1(1:IGDRY) = MAX( 1.0001, MIN( REAL(NDRYLBDAG)-0.0001,           &
                        XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + XDRYINTP2G ) )
  IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
  ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - REAL( IVEC1(1:IGDRY) )
!
  ZVEC2(1:IGDRY) = MAX( 1.0001, MIN( REAL(NDRYLBDAR)-0.0001,           &
                        XDRYINTP1R * LOG( ZVEC2(1:IGDRY) ) + XDRYINTP2R ) )
  IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
  ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
!
!*       2.3.10 perform the bilinear interpolation of the normalized
!               RDRYG-kernel
!
  DO JJ = 1,IGDRY
    ZVEC3(JJ) =  (  XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                  - XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                        *  ZVEC1(JJ)          &
               - (  XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                  - XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                        * (ZVEC1(JJ) - 1.0)
  END DO
  ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
  WHERE( GDRY(:) )
     ZZW1(:,4) = MIN( PRRS1D(:),XFRDRYG*ZZW(:) * PRRT1D(:) * PCGT1D(:) & !
                                *( PRHODREF(:)**(-XCEXVT+1.) )   &
                    *( XLBRDRYG1/( PLBDAG(:)**2              ) + &
                       XLBRDRYG2/( PLBDAG(:)   * PLBDAR(:)   ) + &
                       XLBRDRYG3/(               PLBDAR(:)**2) ) )    
  END WHERE
!
!                                                                          
!*       2.2.10b perform the bilinear interpolation of the normalized
!               RDRYG-kernel FOR CONCENTRATION
!
  DO JJ = 1,IGDRY
     ZVEC3(JJ) =  (  XKER_N_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                        - XKER_N_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                         			 	                  * ZVEC1(JJ) &
                    - (  XKER_N_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                        - XKER_N_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                     			     * (ZVEC1(JJ) - 1.0)
  END DO
  ZZNW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
  WHERE( GDRY(:)  )
     ZZNW1(:,4) = MIN( PCRS1D(:),XFNRDRYG*ZZNW(:) * PCRT1D(:) * PCGT1D(:) & ! NRDRYG
                                *( PRHODREF(:)**(-XCEXVT+1.) )      &
                    *( XLBNRDRYG1/( PLBDAG(:)**2              ) +   &
                       XLBNRDRYG2/( PLBDAG(:)   * PLBDAR(:)   ) +   &
                       XLBNRDRYG3/(               PLBDAR(:)**2) ) )
  END WHERE
  DEALLOCATE(IVEC2)
  DEALLOCATE(IVEC1)
  DEALLOCATE(ZVEC3)
  DEALLOCATE(ZVEC2)
  DEALLOCATE(ZVEC1)
END IF
!
ZRDRYG(:) = ZZW1(:,1) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,4)
!
!
!*       2.4  Compute the Wet growth case
!        --------------------------------
!
ZZW(:) = 0.0
ZRWETG(:) = 0.0
WHERE( PRGT1D(:)>XRTMIN(6) )
   ZZW1(:,5) = MIN( PRIS1D(:),                                    &
               ZZW1(:,2) / (XCOLIG*EXP(XCOLEXIG*(PZT(:)-XTT)) ) ) ! RIWETG
   ZZNW1(:,5) = MIN( PCIS1D(:),                      &  
               ZZNW1(:,2) / (XCOLIG*EXP(XCOLEXIG*(PZT(:)-XTT)) ) )  ! NIWETG
   ZZW1(:,6) = MIN( PRSS1D(:),                                    &
               ZZW1(:,3) / (XCOLSG*EXP(XCOLEXSG*(PZT(:)-XTT)) ) ) ! RSWETG
!
!
   ZZNW1(:,6) = MIN( PCSS1D(:),                      &       
               ZZNW1(:,3) / (XCOLSG*EXP(XCOLEXSG*(PZT(:)-XTT)) ) )  ! NSWETG
!
   ZZW(:) = PRVT1D(:)*PPRES(:)/((XMV/XMD)+PRVT1D(:)) ! Vapor pressure
   ZZW(:) =  PKA(:)*(XTT-PZT(:)) +                                  &
                ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT ))   &
                           *(XESTT-ZZW(:))/(XRV*PZT(:))             )
!  
! compute RWETG
!
   ZRWETG(:)  = MAX( 0.0,                                               &
                   ( ZZW(:) * PCGT1D(:) * ( X0DEPG* PLBDAG(:)**XEX0DEPG + &
                    X1DEPG*PCJ(:)*PLBDAG(:)**XEX1DEPG ) +   &
                   ( ZZW1(:,5)+ZZW1(:,6) ) *                            &
                   ( PRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-PZT(:)))   ) ) / &
                              ( PRHODREF(:)*(XLMTT-XCL*(XTT-PZT(:)))  ) )
  !We must agregate, at least, the cold species
   ZRWETG(:)=MAX(ZRWETG(:), ZZW1(:,5)+ZZW1(:,6))
END WHERE
!
!
!*       2.5  Select Wet or Dry case
!        ---------------------------
!
! Wet case and partial conversion to hail
!
ZZW(:) = 0.0
NHAIL = 0.
IF (LHAIL) NHAIL = 1.
DO JJ=1, SIZE(PRGT1D)
   IF ( PRGT1D(JJ)>XRTMIN(6) .AND. PZT(JJ)<XTT .AND. &
        (ZRDRYG(JJ)-ZZW1(JJ,2)-ZZW1(JJ,3))>(ZRWETG(JJ)-ZZW1(JJ,5)-ZZW1(JJ,6)) .AND. (ZRWETG(JJ)-ZZW1(JJ,5)-ZZW1(JJ,6))>0.0 ) THEN
      !
      ZZW(JJ) = ZRWETG(JJ) - ZZW1(JJ,5) - ZZW1(JJ,6) ! RCWETG+RRWETG
!   
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!   
      ZZW1(JJ,7) = MAX( 0.0,MIN( ZZW(JJ),PRRS1D(JJ)+ZZW1(JJ,1) ) )
      ZZX(JJ)    = ZZW1(JJ,7) / ZZW(JJ)
      ZZW1(JJ,5) = ZZW1(JJ,5)*ZZX(JJ)
      ZZNW1(:,5) = ZZNW1(JJ,5)*ZZX(JJ)
      ZZW1(JJ,6) = ZZW1(JJ,6)*ZZX(JJ)
      ZZNW1(JJ,6) = ZZNW1(JJ,6)*ZZX(JJ)
      ZRWETG(JJ) = ZZW1(JJ,7) + ZZW1(JJ,5) + ZZW1(JJ,6)
!   
      PRCS1D(JJ) = PRCS1D(JJ) - ZZW1(JJ,1)
      PRIS1D(JJ) = PRIS1D(JJ) - ZZW1(JJ,5)
      PRSS1D(JJ) = PRSS1D(JJ) - ZZW1(JJ,6)
!
! assume a linear percent of conversion of graupel into hail
!
      PRGS1D(JJ) = PRGS1D(JJ) + ZRWETG(JJ)
      ZZW(JJ)  = PRGS1D(JJ)*ZRDRYG(JJ)*NHAIL/(ZRWETG(JJ)+ZRDRYG(JJ)) 
      PRGS1D(JJ) = PRGS1D(JJ) - ZZW(JJ)                        
      PRHS1D(JJ) = PRHS1D(JJ) + ZZW(JJ)
      PRRS1D(JJ) = MAX( 0.0,PRRS1D(JJ) - ZZW1(JJ,7) + ZZW1(JJ,1) )
      PTHS1D(JJ) = PTHS1D(JJ) + ZZW1(JJ,7) * (PLSFACT(JJ) - PLVFACT(JJ))
                                                ! f(L_f*(RCWETG+RRWETG))
!
      PCCS1D(JJ) = MAX( PCCS1D(JJ)-ZZW1(JJ,1)*(PCCT1D(JJ)/MAX(PRCT1D(JJ),XRTMIN(2))),0.0 )
      PCIS1D(JJ) = MAX( PCIS1D(JJ)-ZZNW1(JJ,5),0.0 )
      PCRS1D(JJ) = MAX( PCRS1D(JJ)-MAX( ZZW1(JJ,7)-ZZW1(JJ,1),0.0 )                 &
           *(PCRT1D(JJ)/MAX(PRRT1D(JJ),XRTMIN(3))),0.0 )
      PCSS1D(JJ) = MAX( PCSS1D(JJ)-ZZNW1(JJ,6),0.0 )
      ZZNW(JJ)  = PCGS1D(JJ)*ZRDRYG(JJ)*NHAIL/(ZRWETG(JJ)+ZRDRYG(JJ))
      PCGS1D(JJ) = MAX( PCGS1D(JJ)-ZZNW(JJ),0.0 )
      PCHS1D(JJ) = MAX( PCHS1D(JJ)+ZZNW(JJ),0.0 )
   END IF
END DO
!
! Budget storage
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'WETG', &
                                         Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'WETG', &
                                         Unpack( prcs1d(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'WETG', &
                                         Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'WETG', &
                                         Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'WETG', &
                                         Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'WETG', &
                                         Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'WETG', &
                                         Unpack( prhs1d(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'WETG', &
                                         Unpack( pccs1d(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'WETG', &
                                         Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'WETG', &
                                         Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
       if(M2_ICE) then
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'WETG', &
                                         Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'WETG', &
                                         Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) )
                  if (LHAIL) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nh), 'WETG', &
                                         Unpack( pchs1d(:), mask = gmicro(:, :, :), field = pchs(:, :, :) ) * prhodj(:, :, :) )
       end if
  end if
end if
!
! Dry case
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'DRYG', &
                                          Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'DRYG', &
                                          Unpack( prcs1d(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'DRYG', &
                                          Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'DRYG', &
                                          Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'DRYG', &
                                          Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'DRYG', &
                                          Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'DRYG', &
                                          Unpack( pccs1d(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'DRYG', &
                                          Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'DRYG', &
                                          Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
    if (M2_ICE)   call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'DRYG', &
                                          Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
  end if
end if
!
WHERE( PRGT1D(:)>XRTMIN(6) .AND. PZT(:)<XTT                              &
     .AND. (ZRDRYG(:)-ZZW1(:,2)-ZZW1(:,3))<(ZRWETG(:)-ZZW1(:,5)-ZZW1(:,6)) .AND. ZRDRYG(:)>0.0 ) ! case
   PRCS1D(:) = PRCS1D(:) - ZZW1(:,1)
   PRIS1D(:) = PRIS1D(:) - ZZW1(:,2)
   PRSS1D(:) = PRSS1D(:) - ZZW1(:,3)
   PRRS1D(:) = PRRS1D(:) - ZZW1(:,4)
   PRGS1D(:) = PRGS1D(:) + ZRDRYG(:)
   PTHS1D(:) = PTHS1D(:) + (ZZW1(:,1)+ZZW1(:,4))*(PLSFACT(:)-PLVFACT(:)) !
  						        ! f(L_f*(RCDRYG+RRDRYG))
!
   PCCS1D(:) = MAX( PCCS1D(:)-ZZNW1(:,1),0.0 )                                    
   PCIS1D(:) = MAX( PCIS1D(:)-ZZNW1(:,2),0.0 )                                     
   PCRS1D(:) = MAX( PCRS1D(:)-ZZNW1(:,4),0.0 )                                    
   PCSS1D(:) = MAX( PCSS1D(:)-ZZNW1(:,3),0.0 )
END WHERE
!
! Budget storage
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'DRYG', &
                                         Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'DRYG', &
                                         Unpack( prcs1d(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'DRYG', &
                                         Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'DRYG', &
                                         Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'DRYG', &
                                         Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'DRYG', &
                                         Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'DRYG', &
                                         Unpack( pccs1d(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'DRYG', &
                                         Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'DRYG', &
                                         Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
    if (M2_ICE) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'DRYG', &
                                         Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
  end if
end if
!
!
!*       2.6  Hallett-Mossop ice multiplication process due to graupel riming
!        --------------------------------------------------------------------
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'HMG', &
                                          Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'HMG', &
                                          Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HMG', &
                                          Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
end if

GDRY(:) = (PZT(:)<XHMTMAX) .AND. (PZT(:)>XHMTMIN)            &
     .AND. (PRGT1D(:)>XRTMIN(6)) .AND. (PRCT1D(:)>XRTMIN(2)) &
     .AND. (ZRDRYG(:)-ZZW1(:,2)-ZZW1(:,3))<(ZRWETG(:)-ZZW1(:,5)-ZZW1(:,6))
IGDRY = COUNT( GDRY(:) )
IF( IGDRY>0 ) THEN
  ALLOCATE(ZVEC1(IGDRY))
  ALLOCATE(ZVEC2(IGDRY))
  ALLOCATE(IVEC2(IGDRY))
!
  ZVEC1(:) = PACK( PLBDAC(:),MASK=GDRY(:) )
  ZVEC2(1:IGDRY) = MAX( 1.0001, MIN( REAL(NGAMINC)-0.0001,           &
                        XHMLINTP1 * LOG( ZVEC1(1:IGDRY) ) + XHMLINTP2 ) )
  IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
  ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
  ZVEC1(1:IGDRY) =   XGAMINC_HMC( IVEC2(1:IGDRY)+1 )* ZVEC2(1:IGDRY)      &
                   - XGAMINC_HMC( IVEC2(1:IGDRY)   )*(ZVEC2(1:IGDRY) - 1.0)
  ZZX(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GDRY,FIELD=0.0 ) ! Large droplets
!
  WHERE ( GDRY(:) .AND. ZZX(:)<0.99 ) ! Dry case
    ZZW1(:,5) = ZZW1(:,1)*(PCCT1D(:)/PRCT1D(:))*(1.0-ZZX(:))*XHM_FACTG*  &
         MAX( 0.0, MIN( (PZT(:)-XHMTMIN)/3.0,(XHMTMAX-PZT(:))/2.0 ) ) ! CCHMGI
    PCIS1D(:) = PCIS1D(:) + ZZW1(:,5)
!
    ZZW1(:,6) = ZZW1(:,5) * XMNU0                                     ! RCHMGI
    PRIS1D(:) = PRIS1D(:) + ZZW1(:,6)
    PRGS1D(:) = PRGS1D(:) - ZZW1(:,6)
  END WHERE
  DEALLOCATE(IVEC2)
  DEALLOCATE(ZVEC2)
  DEALLOCATE(ZVEC1)
END IF
!
! Budget storage
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'HMG', &
                                         Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'HMG', &
                                         Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HMG', &
                                         Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
end if
!
!*       2.7  Melting of the graupeln
!        ----------------------------
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'GMLT', &
                                          Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'GMLT', &
                                          Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'GMLT', &
                                          Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'GMLT', &
                                          Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv .and. M2_ICE) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'GMLT', &
                                          Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) )              
end if

ZZW(:) = 0.0
   WHERE( (PRGT1D(:)>XRTMIN(6)) .AND. (PRGS1D(:)>XRTMIN(6)/PTSTEP) .AND. (PZT(:)>XTT) )
      ZZW(:) = PRVT1D(:)*PPRES(:)/((XMV/XMD)+PRVT1D(:)) ! Vapor pressure
      ZZW(:) =  PKA(:)*(XTT-PZT(:)) +                                 &
              ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT )) &
                          *(XESTT-ZZW(:))/(XRV*PZT(:))             )
!
! compute RGMLTR
!
      ZZW(:)  = MIN( PRGS1D(:), MAX( 0.0,( -ZZW(:) * PCGT1D(:) *           & 
                          ( X0DEPG*       PLBDAG(:)**XEX0DEPG +     &
                            X1DEPG*PCJ(:)*PLBDAG(:)**XEX1DEPG ) -   &
                                    ( ZZW1(:,1)+ZZW1(:,4) ) *       &
                             ( PRHODREF(:)*XCL*(XTT-PZT(:))) ) /    &
                                            ( PRHODREF(:)*XLMTT ) ) )
      PRRS1D(:) = PRRS1D(:) + ZZW(:)
      PRGS1D(:) = PRGS1D(:) - ZZW(:)
      PTHS1D(:) = PTHS1D(:) - ZZW(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(-RGMLTR))
!
      PCRS1D(:) = PCRS1D(:) + ZZW(:)*5.0E6  ! obtained after averaging
                                        ! Dshed=1mm and 500 microns
      PCGS1D(:) = MAX( PCGS1D(:) - ZZW(:)*5.0E6, 0.0)
   END WHERE    
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'GMLT', &
                                         Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'GMLT', &
                                         Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'GMLT', &
                                         Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'GMLT', &
                                         Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv .and. M2_ICE) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'GMLT', &
                                          Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) )
end if
!
!
!------------------------------------------------------------------------------
!
!                         #################
!                         FAST RH PROCESSES
!                         #################
!
!
HAIL: IF (LHAIL) THEN
!
GHAIL(:) = PRHT1D(:)>XRTMIN(7)
IHAIL = COUNT(GHAIL(:))
!
IF( IHAIL>0 ) THEN
!
!*       3.1 Wet growth of hail 
!        ----------------------------
!
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'WETH', &
                                            Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'WETH', &
                                            Unpack( prcs1d(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'WETH', &
                                            Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'WETH', &
                                            Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'WETH', &
                                            Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'WETH', &
                                            Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'WETH', &
                                            Unpack( prhs1d(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
                    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'WETH', &
                                            Unpack( pccs1d(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'WETH', &
                                            Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'WETH', &
                                            Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
      if (M2_ICE) then
                    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'WETH', &
                                           Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
                    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'WETH', &
                                           Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) )                                     
      end if
    end if
  end if

  ZZW1(:,:) = 0.0
  ZZNW1(:,:) = 0.0
  WHERE( GHAIL(:) .AND. ( (PRCT1D(:)>XRTMIN(2) .AND. PRCS1D(:)>XRTMIN(2)/PTSTEP ) .OR. &
                              (PRIT1D(:)>XRTMIN(4) .AND. PRIS1D(:)>XRTMIN(4)/PTSTEP)  )    )   
     ZZW(:) = PCHT1D(:) * PLBDAH(:)**(-XDH-2.0) * PRHODREF(:)**(1-XCEXVT)
     ZZW1(:,1) = MIN( PRCS1D(:),XFWETH * PRCT1D(:) * ZZW(:) )             ! RCWETH
     ZZNW1(:,1) = MIN( PCCS1D(:),XFWETH * PCCT1D(:) * ZZW(:) ) !ZZW1(:,1) * ZCCT(:) / MAX(ZRCT(:),XRTMIN(2)) )            ! NCWETH 
     ZZW1(:,2) = MIN( PRIS1D(:),XFWETH * PRIT1D(:) * ZZW(:) )             ! RIWETH
     ZZNW1(:,2) = MIN( PCIS1D(:),XFWETH * PCIT1D(:) * ZZW(:) ) !ZZW1(:,2) * ZCIT(:) / MAX(ZRIT(:),XRTMIN(4)) )      ! NIWETH 
  END WHERE
!
!*       3.1.1  accretion of aggregates on the hailstones
!        ------------------------------------------------
!
  GWET(:) = GHAIL(:) .AND. (PRST1D(:)>XRTMIN(5) .AND. PRSS1D(:)>XRTMIN(5)/PTSTEP)
  IGWET = COUNT( GWET(:) )
!
  IF( IGWET>0 ) THEN
!
!*       3.1.2  allocations
!
    ALLOCATE(ZVEC1(IGWET))
    ALLOCATE(ZVEC2(IGWET))
    ALLOCATE(ZVEC3(IGWET))
    ALLOCATE(IVEC1(IGWET))
    ALLOCATE(IVEC2(IGWET))
!
!*       3.1.3  select the (PLBDAH,PLBDAS) couplet
!
    ZVEC1(:) = PACK( PLBDAH(:),MASK=GWET(:) )
    ZVEC2(:) = PACK( PLBDAS(:),MASK=GWET(:) )
!
!*       3.1.4  find the next lower indice for the PLBDAG and for the PLBDAS
!               in the geometrical set of (Lbda_h,Lbda_s) couplet use to
!               tabulate the SWETH-kernel
!
    ZVEC1(1:IGWET) = MAX( 1.0001, MIN( REAL(NWETLBDAH)-0.0001,           &
                          XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + XWETINTP2H ) )
    IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
    ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - REAL( IVEC1(1:IGWET) )
!
    ZVEC2(1:IGWET) = MAX( 1.0001, MIN( REAL(NWETLBDAS)-0.0001,           &
                          XWETINTP1S * LOG( ZVEC2(1:IGWET) ) + XWETINTP2S ) )
    IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
    ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - REAL( IVEC2(1:IGWET) )
!
!*       3.1.5  perform the bilinear interpolation of the normalized
!               SWETH-kernel
!
    DO JJ = 1,IGWET
      ZVEC3(JJ) = (  XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                   - XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                         *  ZVEC1(JJ)          &
                 - ( XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                   - XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
    WHERE( GWET(:) )
       ZZW1(:,3) = MIN( PRSS1D(:),XFSWETH*ZZW(:) * PRST1D(:) * PCHT1D(:)  & 
            *( PRHODREF(:)**(-XCEXVT+1.) )               &
            *( XLBSWETH1/( PLBDAH(:)**2              ) + &
            XLBSWETH2/( PLBDAH(:)   * PLBDAS(:)   ) + &
            XLBSWETH3/(               PLBDAS(:)**2) ) )   
    END WHERE
!
!*       3.1.5b  perform the bilinear interpolation of the normalized
!               SWETH-kernel FOR CONCENTRATION
!
    DO JJ = 1,IGWET
       ZVEC3(JJ) = (  XKER_N_SWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_N_SWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
        		 	                                   * ZVEC1(JJ) &
                   - ( XKER_N_SWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_N_SWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
         		                                     * (ZVEC1(JJ) - 1.0)
    END DO
    ZZNW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
    WHERE( GWET(:) )
       ZZNW1(:,3) = MIN( PCSS1D(:),XFNSWETH*ZZNW(:) * PCST1D(:)            & ! NSWETH
                                              *( PCHT1D(:) )              & 
       	                 *( PRHODREF(:)**(-XCEXVT+1.) )               &
                         *( XLBNSWETH1/( PLBDAH(:)**2              ) + &
                            XLBNSWETH2/( PLBDAH(:)   * PLBDAS(:)   ) + &
                            XLBNSWETH3/(               PLBDAS(:)**2) ) )
    END WHERE
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
  END IF
!
!*       3.1.6  accretion of graupeln on the hailstones
!        ----------------------------------------------
!
  GWET(:) = GHAIL(:) .AND. (PRGT1D(:)>XRTMIN(6) .AND. PRGS1D(:)>XRTMIN(6)/PTSTEP)
  IGWET = COUNT( GWET(:) )
!
  IF( IGWET>0 ) THEN
!
!*       3.1.7  allocations
!
    ALLOCATE(ZVEC1(IGWET))
    ALLOCATE(ZVEC2(IGWET))
    ALLOCATE(ZVEC3(IGWET))
    ALLOCATE(IVEC1(IGWET))
    ALLOCATE(IVEC2(IGWET))
!
!*       3.1.8  select the (PLBDAH,PLBDAG) couplet
!
    ZVEC1(:) = PACK( PLBDAH(:),MASK=GWET(:) )
    ZVEC2(:) = PACK( PLBDAG(:),MASK=GWET(:) )
!
!*       3.1.9  find the next lower indice for the PLBDAH and for the PLBDAG
!               in the geometrical set of (Lbda_h,Lbda_g) couplet use to
!               tabulate the GWETH-kernel
!
    ZVEC1(1:IGWET) = MAX( 1.0001, MIN( REAL(NWETLBDAG)-0.0001,           &
                          XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + XWETINTP2H ) )
    IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
    ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - REAL( IVEC1(1:IGWET) )
!
    ZVEC2(1:IGWET) = MAX( 1.0001, MIN( REAL(NWETLBDAG)-0.0001,           &
                          XWETINTP1G * LOG( ZVEC2(1:IGWET) ) + XWETINTP2G ) )
    IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
    ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - REAL( IVEC2(1:IGWET) )
!
!*       3.1.10 perform the bilinear interpolation of the normalized
!               GWETH-kernel
!
    DO JJ = 1,IGWET
      ZVEC3(JJ) = (  XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                   - XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                         *  ZVEC1(JJ)          &
                - (  XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)+1) *  ZVEC2(JJ)          &
                   - XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)  ) * (ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
    WHERE( GWET(:) )
       ZZW1(:,5) = MAX(MIN( PRGS1D(:),XFGWETH*ZZW(:) * PRGT1D(:) * PCHT1D(:) &  
            *( PRHODREF(:)**(-XCEXVT+1.) )               &
            *( XLBGWETH1/( PLBDAH(:)**2              ) + &
            XLBGWETH2/( PLBDAH(:)   * PLBDAG(:)   ) + &
            XLBGWETH3/(               PLBDAG(:)**2) ) ),0. )
    END WHERE
!*       3.1.10b perform the bilinear interpolation of the normalized
!               GWETH-kernel FOR CONCENTRATION
!
    DO JJ = 1,IGWET
       ZVEC3(JJ) = (  XKER_N_GWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                        - XKER_N_GWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                          			 	                   * ZVEC1(JJ) &
                     - (  XKER_N_GWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                        - XKER_N_GWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                    			     * (ZVEC1(JJ) - 1.0)
    END DO
    ZZNW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
    WHERE( GWET(:) )
       ZZNW1(:,5) = MAX(MIN( PCGS1D(:),XFNGWETH*ZZNW(:) * PCGT1D(:)       & ! NGWETH
                                                   *( PCHT1D(:) )        &
                         *( PRHODREF(:)**(-XCEXVT+1.) )               &
                         *( XLBNGWETH1/( PLBDAH(:)**2              ) + &
                            XLBNGWETH2/( PLBDAH(:)   * PLBDAG(:)   ) + &
                            XLBNGWETH3/(               PLBDAG(:)**2) ) ),0. )
    END WHERE
! 
    DEALLOCATE(IVEC2)
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC3)
    DEALLOCATE(ZVEC2)
    DEALLOCATE(ZVEC1)
  END IF
!
!
!*       3.2    compute the Wet growth of hail
!        -------------------------------------
!
  ZZW(:) = 0.0
  WHERE( GHAIL(:) .AND. PZT(:)<XTT )
     ZZW(:) = PRVT1D(:)*PPRES(:)/((XMV/XMD)+PRVT1D(:)) ! Vapor pressure
     ZZW(:) = PKA(:)*(XTT-PZT(:)) +                                 &
                    ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT )) &
                                *(XESTT-ZZW(:))/(XRV*PZT(:))             )
!
! compute RWETH
!
     ZZW(:)  =  MAX(0.,  ( ZZW(:) * PCHT1D(:) * ( X0DEPH*       PLBDAH(:)**XEX0DEPH +     &  
                                    X1DEPH*PCJ(:)*PLBDAH(:)**XEX1DEPH ) +   &
                       ( ZZW1(:,2)+ZZW1(:,3)+ZZW1(:,5) ) *                  &
                       ( PRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-PZT(:)))   ) ) / &
                             ( PRHODREF(:)*(XLMTT-XCL*(XTT-PZT(:))) ) )            
!
     ZZW1(:,6) = MAX( ZZW(:) - ZZW1(:,2) - ZZW1(:,3) - ZZW1(:,5),0.) !RCWETH+RRWETH
  END WHERE
   !
  WHERE ( GHAIL(:) .AND. PZT(:)<XTT  .AND. ZZW1(:,6)/=0.)
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
     ZZW1(:,4) = MAX( 0.0,MIN( ZZW1(:,6),PRRS1D(:)+ZZW1(:,1) ) )
     ZZX(:)    = ZZW1(:,4) / ZZW1(:,6)
     ZZW1(:,2) = ZZW1(:,2) * ZZX(:)
     ZZNW1(:,2) = ZZNW1(:,2)*ZZX(:) 
     ZZW1(:,3) = ZZW1(:,3) * ZZX(:)
     ZZNW1(:,3) = ZZNW1(:,3)*ZZX(:)
     ZZW1(:,5) = ZZW1(:,5) * ZZX(:)
     ZZNW1(:,5) = ZZNW1(:,5)*ZZX(:)  
     ZZW(:)    = ZZW1(:,4) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,5)
!
!*       3.2.1  integrate the Wet growth of hail
!
     PRCS1D(:) = PRCS1D(:) - ZZW1(:,1)
     PRIS1D(:) = PRIS1D(:) - ZZW1(:,2)
     PRSS1D(:) = PRSS1D(:) - ZZW1(:,3)
     PRGS1D(:) = PRGS1D(:) - ZZW1(:,5)
     PRHS1D(:) = PRHS1D(:) + ZZW(:)
     PRRS1D(:) = MAX( 0.0,PRRS1D(:) - ZZW1(:,4) + ZZW1(:,1) )
     PTHS1D(:) = PTHS1D(:) + ZZW1(:,4) * (PLSFACT(:) - PLVFACT(:)) 
                                                ! f(L_f*(RCWETH+RRWETH))
!
     PCRS1D(:) = MAX( PCRS1D(:)-MAX( ZZW1(:,4)-ZZW1(:,1),0.0 )                 &
                                       *(PCRT1D(:)/MAX(PRRT1D(:),XRTMIN(3))),0.0 )
     PCSS1D(:) = MAX( PCSS1D(:)-ZZNW1(:,3),0.0 ) 
     PCGS1D(:) = MAX( PCGS1D(:)-ZZNW1(:,5),0.0 )   
     PCCS1D(:) = MAX( PCCS1D(:)-ZZNW1(:,1),0.0 )       
     PCIS1D(:) = MAX( PCIS1D(:)-ZZNW1(:,2),0.0 )
  END WHERE
!
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'WETH', &
                                           Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'WETH', &
                                           Unpack( prcs1d(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'WETH', &
                                           Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'WETH', &
                                           Unpack( pris1d(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'WETH', &
                                           Unpack( prss1d(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'WETH', &
                                           Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'WETH', &
                                           Unpack( prhs1d(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
                    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'WETH', &
                                           Unpack( pccs1d(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'WETH', &
                                           Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'WETH', &
                                           Unpack( pcis1d(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
      if (M2_ICE) then
                    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'WETH', &
                                           Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
                    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'WETH', &
                                           Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) ) 
      end if
    end if
  end if
END IF ! IHAIL>0
!
! Partial reconversion of hail to graupel when rc and rh are small    
!
!
!*       3.3   Conversion of the hailstones into graupel
!        -----------------------------------------------
!
IF ( IHAIL>0 ) THEN
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'COHG', &
                                            Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'COHG', &
                                            Unpack( prhs1d(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
    if (lbudget_sv .and. M2_ICE) then
                      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'COHG', &
                                           Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) ) 
                      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nh), 'COHG', &
                                           Unpack( pchs1d(:), mask = gmicro(:, :, :), field = pchs(:, :, :) ) * prhodj(:, :, :) )  
    end if
  end if
!
  ZTHRH = 0.01E-3
  ZTHRC = 0.001E-3
  ZZW(:) = 0.0
!
  WHERE( PRHT1D(:)<ZTHRH .AND. PRCT1D(:)<ZTHRC .AND. PZT(:)<XTT )
    ZZW(:) = MIN( 1.0,MAX( 0.0,1.0-(PRCT1D(:)/ZTHRC) ) )
!
! assume a linear percent conversion rate of hail into graupel
!
    ZZW(:)  = PRHS1D(:) * ZZW(:)
    PRGS1D(:) = PRGS1D(:) + ZZW(:)                      !   partial conversion
    PRHS1D(:) = PRHS1D(:) - ZZW(:)                      ! of hail into graupel
  END WHERE
    if (M2_ICE) then
       WHERE( PRHT1D(:)<ZTHRH .AND. PRCT1D(:)<ZTHRC .AND. PZT(:)<XTT )
          PCGS1D(:) = PCGS1D(:) + ZZW(:)* PCHS1D(:)/MAX(PRHS1D(:),XRTMIN(7)) 
          PCHS1D(:) = MAX( PCHS1D(:) - ZZW(:)* PCHS1D(:)/MAX(PRHS1D(:),XRTMIN(7)), 0.0 ) 
       END WHERE
    end if
!
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'COHG', &
                                           Unpack( prgs1d(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'COHG', &
                                           Unpack( prhs1d(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
    if (lbudget_sv .and. M2_ICE) then
                      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ng), 'COHG', &
                                           Unpack( pcgs1d(:), mask = gmicro(:, :, :), field = pcgs(:, :, :) ) * prhodj(:, :, :) ) 
                      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nh), 'COHG', &
                                           Unpack( pchs1d(:), mask = gmicro(:, :, :), field = pchs(:, :, :) ) * prhodj(:, :, :) )                                   
    end if
  end if
END IF
!
!
!*       3.4    Melting of the hailstones
!
IF ( IHAIL>0 ) THEN
   if ( nbumod == kmi .and. lbu_enable ) then
      if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HMLT', &
           Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
      if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'HMLT', &
           Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
      if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'HMLT', &
           Unpack( prhs1d(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
      if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'HMLT', &
           Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
      if ( lbudget_sv .and. M2_ICE ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nh), 'HMLT', &
           Unpack( pchs1d(:), mask = gmicro(:, :, :), field = pchs(:, :, :) ) * prhodj(:, :, :) )                                    
   end if
!
   ZZW(:) = 0.0
   WHERE( GHAIL(:) .AND. (PRHS1D(:)>XRTMIN(7)/PTSTEP) .AND. (PRHT1D(:)>XRTMIN(7)) .AND. (PZT(:)>XTT) )
      ZZW(:) = PRVT1D(:)*PPRES(:)/((XMV/XMD)+PRVT1D(:)) ! Vapor pressure
      ZZW(:) = PKA(:)*(XTT-PZT(:)) +                              &
              ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PZT(:) - XTT )) &
              *(XESTT-ZZW(:))/(XRV*PZT(:))         )
!
! compute RHMLTR
!
      ZZW(:)  = MIN( PRHS1D(:), MAX( 0.0,( -ZZW(:) * PCHT1D(:) *           &
              ( X0DEPH*       PLBDAH(:)**XEX0DEPH +     &
              X1DEPH*PCJ(:)*PLBDAH(:)**XEX1DEPH ) -   &
              ZZW1(:,6)*( PRHODREF(:)*XCL*(XTT-PZT(:))) ) /    &
              ( PRHODREF(:)*XLMTT ) ) )  
      PRRS1D(:) = PRRS1D(:) + ZZW(:)
      PRHS1D(:) = PRHS1D(:) - ZZW(:)
      PTHS1D(:) = PTHS1D(:) - ZZW(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(-RHMLTR))
!
      PCRS1D(:) = MAX( PCRS1D(:) + ZZW(:)*(PCHT1D(:)/PRHT1D(:)),0.0 )
      PCHS1D(:) = MAX( PCHS1D(:) - ZZW(:)*(PCHT1D(:)/PRHT1D(:)),0.0 )
!
   END WHERE
!
   if ( nbumod == kmi .and. lbu_enable ) then
      if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HMLT', &
           Unpack( pths1d(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
      if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'HMLT', &
           Unpack( prrs1d(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
      if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'HMLT', &
           Unpack( prhs1d(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
      if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'HMLT', &
           Unpack( pcrs1d(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
      if ( lbudget_sv .and. M2_ICE ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nh), 'HMLT', &
           Unpack( pchs1d(:), mask = gmicro(:, :, :), field = pchs(:, :, :) ) * prhodj(:, :, :) )
   end if
END IF
!
END IF HAIL
!
!
!
if( M2_ICE) then
   if ( nbumod == kmi .and. lbu_enable ) then
      if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'SSC', &
           Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
   end if
!*       1.3NEW  Aggregates self-collection  
!        -------------------------------
!
   ZZNW1(:,:) = 0.0
!
   GACC(:) = PCST1D(:)>XCTMIN(5) .AND. PRST1D(:)>XRTMIN(5) .AND. PZT(:)<XTT
   IGACC = COUNT( GACC(:) )
!
   IF( IGACC>0 ) THEN
!
!        1.3N.0  allocations
!
      ALLOCATE(ZVEC1(IGACC))
      ALLOCATE(IVEC1(IGACC))
!
!        1.3N.1  select the (ZLBDAS,ZLBDAS) couplet
!
      ZVEC1(:) = PACK( PLBDAS(:),MASK=GACC(:) )
!
!        1.3N.2  find the next lower indice for the ZLBDAS and for the ZLBDAS
!               in the geometrical set of (Lbda_s,Lbda_s) couplet use to
!               tabulate the SACCS-kernel
!
      ZVEC1(1:IGACC) = MAX( 1.0001, MIN( FLOAT(NSCLBDAS)-0.0001,           &
           XSCINTP1S * LOG( ZVEC1(1:IGACC) ) + XSCINTP2S ) )
      IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
      ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - FLOAT( IVEC1(1:IGACC) )
!
!        1.3N.3 perform the bilinear interpolation of the normalized
!               SSCS-kernel
!
      ALLOCATE(ZVEC3(IGACC))
      DO JJ = 1,IGACC
         ZVEC3(JJ) =  (   XKER_N_SSCS(IVEC1(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    -  XKER_N_SSCS(IVEC1(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                                         * ZVEC1(JJ) &
                 - (   XKER_N_SSCS(IVEC1(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    -  XKER_N_SSCS(IVEC1(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                           * (ZVEC1(JJ) - 1.0)
      END DO
      ZZNW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC(:),FIELD=0.0 ) !! NSACCS
      DEALLOCATE(ZVEC3)
!
      WHERE( GACC(:) )               
!          ZZNW1(:,5) = XSSC * (PRST1D(:) * PRHODREF(:))**EXPRS * PCST1D(:)**EXPNS      
         ZZNW1(:,5) = MIN( PCSS1D(:),XFNSSCS*ZZNW(:)                       & ! NSSCS
                                      * EXP( XCOLEXSS*(PZT(:)-XTT) )  &
                    *( PCST1D(:)                     )**( 2 )           &  
                    *( PRHODREF(:)**(-XCEXVT-1.) )                    &
                         *( XLBNSSCS1/( PLBDAS(:)**2             ) +  &
                            XLBNSSCS2/( PLBDAS(:)**2 ) )            )
!                                                                           
         PCSS1D(:) = MAX( PCSS1D(:)-ZZNW1(:,5),0.0 )    
      END WHERE
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC1)
   END IF
   if ( nbumod == kmi .and. lbu_enable ) then
      if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'SSC', &
           Unpack( pcss1d(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
   end if
end if
!
!------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_MIXED_FAST_PROCESSES
