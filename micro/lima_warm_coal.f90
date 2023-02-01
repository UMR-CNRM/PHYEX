!MNH_LIC Copyright 2013-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ##########################
       MODULE MODI_LIMA_WARM_COAL
!      ##########################
!
INTERFACE
      SUBROUTINE LIMA_WARM_COAL (PTSTEP, KMI,                                &
                                 PRHODREF, ZWLBDC3, ZWLBDC, ZWLBDR3, ZWLBDR, &
                                 PRCT, PRRT, PCCT, PCRT,                     &
                                 PRCS, PRRS, PCCS, PCRS,                     &
                                 PRHODJ                              )
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step
                                                      ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDC3    ! Lambda(cloud) **3
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDC     ! Lambda(cloud)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDR3    ! Lambda(rain) **3
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDR     ! Lambda(rain)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCCT       ! Cloud water C. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT       ! Rain water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS       ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS       ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS       ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS       ! Rain water C. source
!
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHODJ
!
      END SUBROUTINE LIMA_WARM_COAL
END INTERFACE
END MODULE MODI_LIMA_WARM_COAL
!     #############################################################################
      SUBROUTINE LIMA_WARM_COAL (PTSTEP, KMI,                                &
                                 PRHODREF, ZWLBDC3, ZWLBDC, ZWLBDR3, ZWLBDR, &
                                 PRCT, PRRT, PCCT, PCRT,                     &
                                 PRCS, PRRS, PCCS, PCRS,                     &
                                 PRHODJ                              )
!     #############################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the microphysical sources:
!!    nucleation, sedimentation, autoconversion, accretion, self-collection 
!!    and vaporisation which are parameterized according to Cohard and Pinty 
!!    QJRMS, 2000
!!
!!
!!**  METHOD
!!    ------
!!      Assuming a generalized gamma distribution law for the cloud droplets 
!!    and the raindrops, the zeroth and third order moments tendencies 
!!    are evaluated for all the coalescence terms by integrating the
!!    Stochastic Collection Equation. As autoconversion is a process that 
!!    cannot be resolved analytically, the Berry-Reinhardt parameterisation 
!!    is employed with modifications to initiate the raindrop spectrum mode.
!!     
!!    Computation steps :
!!      1- Check where computations are necessary, pack variables
!!      2- Self collection of cloud droplets
!!      3- Autoconversion of cloud droplets (Berry-Reinhardt parameterization)
!!      4- Accretion sources
!!      5- Self collection - Coalescence/Break-up
!!      6- Unpack variables, clean
!!
!!
!!    REFERENCE
!!    ---------
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
!!      C. Barthe  * LACy *  jan. 2014   add budgets
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets (no more budget calls in this subroutine)
!       Delbeke/Vie     03/2022 : KHKO option
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,          only: lbudget_rc, lbudget_rr, lbudget_sv, NBUDGET_RC, NBUDGET_RR, NBUDGET_SV1, tbudgets
USE MODD_CST,             ONLY: XPI, XRHOLW
USE MODD_NSV,             ONLY: NSV_LIMA_NC, NSV_LIMA_NR
USE MODD_PARAMETERS,      ONLY: JPHEXT, JPVEXT
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_WARM

use mode_budget,          only: Budget_store_init, Budget_store_end
use mode_tools,           only: Countjv

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                     INTENT(IN)    :: PTSTEP     ! Double Time step
                                                      ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KMI        ! Model index 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF   ! Reference density
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDC3    ! Lambda(cloud) **3
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDC     ! Lambda(cloud)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDR3    ! Lambda(rain) **3
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: ZWLBDR     ! Lambda(rain)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT       ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT       ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCCT       ! Cloud water C. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT       ! Rain water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS       ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS       ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS       ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS       ! Rain water C. source
!
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHODJ
!
!*       0.1   Declarations of local variables :
!
! Packing variables
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: GMICRO 
INTEGER :: IMICRO
INTEGER , DIMENSION(SIZE(GMICRO))   :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                             :: JL       ! and PACK intrinsics 
!
! Packed micophysical variables
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCCT    ! cloud conc. at t
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCRT    ! rain conc. at t
!
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZRRS    ! Rain water m.r. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCCS    ! cloud conc. source
REAL, DIMENSION(:)  , ALLOCATABLE :: ZCRS    ! rain conc. source
!
! Other packed variables
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZRHODREF ! RHO Dry REFerence
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZLBDC3
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZLBDC
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZLBDR3
REAL, DIMENSION(:)  , ALLOCATABLE  :: ZLBDR
!
! Work arrays
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) :: ZW
!
REAL, DIMENSION(:), ALLOCATABLE    :: ZZW1, ZZW2, ZZW3, ZZW4, ZSCBU
LOGICAL, DIMENSION(:), ALLOCATABLE :: GSELF,               &
                                      GACCR,               &
                                      GSCBU,               &
                                      GENABLE_ACCR_SCBU
! 
!
INTEGER :: ISELF, IACCR, ISCBU
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE        ! Physical domain
!
!-------------------------------------------------------------------------------
!
!
!*       1.     PREPARE COMPUTATIONS - PACK
!   	        ---------------------------
!
!
IIB=1+JPHEXT
IIE=SIZE(PRHODREF,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PRHODREF,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PRHODREF,3) - JPVEXT
!
GMICRO(:,:,:) = .FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                  &
     PRCT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(2) .OR.  &
     PRRT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(3)
!
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
!
IF( IMICRO >= 0 ) THEN
   ALLOCATE(ZRCT(IMICRO))
   ALLOCATE(ZRRT(IMICRO))
   ALLOCATE(ZCCT(IMICRO))
   ALLOCATE(ZCRT(IMICRO))
!
   ALLOCATE(ZRCS(IMICRO))
   ALLOCATE(ZRRS(IMICRO))
   ALLOCATE(ZCCS(IMICRO))
   ALLOCATE(ZCRS(IMICRO))
!
   ALLOCATE(ZLBDC(IMICRO)) 
   ALLOCATE(ZLBDC3(IMICRO))
   ALLOCATE(ZLBDR(IMICRO)) 
   ALLOCATE(ZLBDR3(IMICRO))
! 
   ALLOCATE(ZRHODREF(IMICRO))
   DO JL=1,IMICRO
      ZCCT(JL) = PCCT(I1(JL),I2(JL),I3(JL))
      ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
      ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
      ZCRT(JL) = PCRT(I1(JL),I2(JL),I3(JL))
      ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
      ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
      ZRRS(JL) = PRRS(I1(JL),I2(JL),I3(JL))
      ZCRS(JL) = PCRS(I1(JL),I2(JL),I3(JL))
      ZLBDR(JL) = ZWLBDR(I1(JL),I2(JL),I3(JL))
      ZLBDR3(JL) = ZWLBDR3(I1(JL),I2(JL),I3(JL))
      ZLBDC(JL) = ZWLBDC(I1(JL),I2(JL),I3(JL))
      ZLBDC3(JL) = ZWLBDC3(I1(JL),I2(JL),I3(JL))
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
   END DO
! 
   ALLOCATE(GSELF(IMICRO))
   ALLOCATE(GACCR(IMICRO))
   ALLOCATE(GSCBU(IMICRO))
   ALLOCATE(ZZW1(IMICRO))
   ALLOCATE(ZZW2(IMICRO))
   ALLOCATE(ZZW3(IMICRO))
!
!
!-------------------------------------------------------------------------------
!
IF (LRAIN) THEN
!
!*       2. Self-collection of cloud droplets    
!   	 ------------------------------------
!
!
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'SELF', pccs(:, :, :)  * prhodj(:, :, :) )

   GSELF(:) = ZCCT(:)>XCTMIN(2)
   ISELF = COUNT(GSELF(:))
   IF( ISELF>=0 .AND. .NOT.LKHKO) THEN
      ZZW1(:) = XSELFC*(ZCCT(:)/ZLBDC3(:))**2 * ZRHODREF(:) ! analytical integration
      WHERE( GSELF(:) )
         ZCCS(:) = ZCCS(:) - MIN( ZCCS(:),ZZW1(:) )
      END WHERE
   END IF

  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'SELF', &
                                           Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )

!-------------------------------------------------------------------------------
!
!
!*       3. Autoconversion of cloud droplets (Berry-Reinhardt parameterization)
!   	 ----------------------------------------------------------------------
!
!
!
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'AUTO', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'AUTO', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'AUTO', pccs(:, :, :) * prhodj(:, :, :) )
    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'AUTO', pcrs(:, :, :) * prhodj(:, :, :) )
  end if

   ZZW2(:) = 0.0
   ZZW1(:) = 0.0
   IF (LKHKO) THEN
      WHERE ( ZRCT(:) .GT. XRTMIN(2) .AND. ZCCT(:) .GT. XCTMIN(2)                 &
            .AND. (ZRCS(:) .GT. 0.0) .AND. (ZCCS(:) .GT. 0.0))
!
         ZZW1(:)= 1350.0 * ZRCT(:)**(2.47) * (ZCCT(:)/1.0E6)**(-1.79) ! ZCCT in cm-3         
         ZZW1(:) = min (ZRCS(:), ZZW1(:))
         ZRCS(:) = ZRCS(:) - ZZW1(:)
         ZRRS(:) = ZRRS(:) + ZZW1(:)
!
         ZCRS(:) = ZCRS(:) + ZZW1(:) * 3. * ZRHODREF(:)/(4.*XPI*XRHOLW*(XR0)**(3.))
!
         ZZW1(:) = min ( ZCCS(:),ZZW1(:) * ZCCT(:) / ZRCT(:))
         ZCCS(:) = ZCCS(:) - ZZW1(:)
!
      END WHERE
   ELSE
      WHERE( ZRCT(:)>XRTMIN(2) )
         ZZW2(:) = MAX( 0.0,XLAUTR*ZRHODREF(:)*ZRCT(:)*             &
                           (XAUTO1/min(ZLBDC(:),1.e9)**4-XLAUTR_THRESHOLD) ) ! L 
!
         ZZW3(:) = MIN( ZRCS(:), MAX( 0.0,XITAUTR*ZZW2(:)*ZRCT(:)*  &
                           (XAUTO2/ZLBDC(:)-XITAUTR_THRESHOLD) ) ) ! L/tau
!
         ZRCS(:) = ZRCS(:) - ZZW3(:)
         ZRRS(:) = ZRRS(:) + ZZW3(:)
!
         ZZW1(:) = MIN( MIN( 1.2E4,(XACCR4/ZLBDC(:)-XACCR5)/XACCR3),   &
                           ZLBDR(:)/XACCR1 ) ! D**-1 threshold diameter for 
                                             ! switching the autoconversion regimes
                                             ! min (80 microns, D_h, D_r)
         ZZW3(:) = ZZW3(:) * MAX( 0.0,ZZW1(:) )**3 / XAC 
         ZCRS(:) = ZCRS(:) + ZZW3(:)
      END WHERE
   END IF
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'AUTO', &
                                           Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'AUTO', &
                                           Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
    !This budget is = 0 for nsv_lima_nc => not necessary to call it (ZCCS is not modified in this part)
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'AUTO', &
                           Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'AUTO', &
                           Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
  end if

!-------------------------------------------------------------------------------
!
!
!*       4. Accretion sources
!   	 --------------------
!
!
   GACCR(:) = ZRRT(:)>XRTMIN(3) .AND. ZCRT(:)>XCTMIN(3)
   IACCR = COUNT(GACCR(:))
   IF( IACCR >= 0 ) THEN
      ALLOCATE(ZZW4(IMICRO)); ZZW4(:) = XACCR1/ZLBDR(:)
      ALLOCATE(GENABLE_ACCR_SCBU(IMICRO))
      GENABLE_ACCR_SCBU(:) = ZRRT(:)>1.2*ZZW2(:)/ZRHODREF(:) .OR.           &
                       ZZW4(:)>=MAX( XACCR2,XACCR3/(XACCR4/ZLBDC(:)-XACCR5) )
      GACCR(:) = GACCR(:) .AND. ZRCT(:)>XRTMIN(2) .AND. ZCCT(:)>XCTMIN(2) .AND. GENABLE_ACCR_SCBU(:)
   END IF
!
   IACCR = COUNT(GACCR(:))
   IF( IACCR >= 0 ) THEN
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'ACCR', &
                                              Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'ACCR', &
                                              Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'ACCR', &
                                              Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
    IF (LKHKO) THEN
       WHERE ( (ZRCT(:) .GT. XRTMIN(2)) .AND. (ZRRT(:) .GT. XRTMIN(3))                 &
               .AND. (ZRCS(:) .GT. 0.0) .AND. (ZCCS(:) .GT. 0.0))
          ZZW1(:) = 67.0 * ( ZRCT(:) * ZRRT(:) )**1.15
          ZZW1(:) = MIN (ZRCS(:),ZZW1(:))
          ZRCS(:) = ZRCS(:) - ZZW1(:)
          ZRRS(:) = ZRRS(:) + ZZW1(:)
!
          ZZW1(:) = MIN (ZCCS(:),ZZW1(:) * ZCCT(:) / ZRCT(:))
          ZCCS(:) = ZCCS(:) - ZZW1(:)
!
       END WHERE
    ELSE
       WHERE( GACCR(:).AND.(ZZW4(:)>1.E-4) ) ! Accretion for D>100 10-6 m
          ZZW3(:) = ZLBDC3(:) / ZLBDR3(:)
          ZZW1(:) = ( ZCCT(:)*ZCRT(:) / ZLBDC3(:) )*ZRHODREF(:)
          ZZW2(:) = MIN( ZZW1(:)*(XACCR_CLARGE1+XACCR_CLARGE2*ZZW3(:)),ZCCS(:) )
          ZCCS(:) = ZCCS(:) - ZZW2(:)
!
          ZZW1(:) = ( ZZW1(:) / ZLBDC3(:) )
          ZZW2(:) = MIN( ZZW1(:)*(XACCR_RLARGE1+XACCR_RLARGE2*ZZW3(:)),ZRCS(:) )
          ZRCS(:) = ZRCS(:) - ZZW2(:)
          ZRRS(:) = ZRRS(:) + ZZW2(:)
       END WHERE
       WHERE( GACCR(:).AND.(ZZW4(:)<=1.E-4) ) ! Accretion for D<100 10-6 m
          ZZW3(:) = MIN(ZLBDC3(:) / ZLBDR3(:), 1.E8)
          ZZW1(:) = ( ZCCT(:)*ZCRT(:) / ZLBDC3(:) )*ZRHODREF(:)
          ZZW1(:) = ZZW1(:) / ZLBDC3(:)
          ZZW3(:) = ZZW3(:)**2
          ZZW2(:) = MIN( ZZW1(:)*(XACCR_CSMALL1+XACCR_CSMALL2*ZZW3(:)),ZCCS(:) )
          ZCCS(:) = ZCCS(:) - ZZW2(:)
!
          ZZW1(:) = ZZW1(:) / ZLBDC3(:)
          ZZW2(:) = MIN( ZZW1(:)*(XACCR_RSMALL1+XACCR_RSMALL2*ZZW3(:))             &
                                                          ,ZRCS(:) )
          ZRCS(:) = ZRCS(:) - ZZW2(:)
          ZRRS(:) = ZRRS(:) + ZZW2(:)
       END WHERE
    END IF

    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'ACCR', &
                                             Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'ACCR', &
                                             Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'ACCR', &
                                             Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
   END IF
!-------------------------------------------------------------------------------
!
!
!*       5. Self collection - Coalescence/Break-up
!   	 -----------------------------------------
!
!
   IF( IACCR >= 0 ) THEN
      GSCBU(:) = ZCRT(:)>XCTMIN(3) .AND. GENABLE_ACCR_SCBU(:)
      ISCBU = COUNT(GSCBU(:))
   ELSE
      ISCBU = 0.0
   END IF
   IF( ISCBU>0 .AND. .NOT.LKHKO) THEN
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'SCBU', &
                                              Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
!
!*       5.1  efficiencies
!
      IF (.NOT.ALLOCATED(ZZW4)) ALLOCATE(ZZW4(IMICRO))
      ZZW4(:)  = XACCR1 / ZLBDR(:)                ! Mean diameter
      ALLOCATE(ZSCBU(IMICRO))
      ZSCBU(:) = 1.0
      WHERE (ZZW4(:)>=XSCBU_EFF1 .AND. GSCBU(:))   ZSCBU(:) = &  ! Coalescence
                            EXP(XSCBUEXP1*(ZZW4(:)-XSCBU_EFF1))  ! efficiency
      WHERE (ZZW4(:)>=XSCBU_EFF2) ZSCBU(:) = 0.0  ! Break-up
!
!*       5.2  integration
!
      ZZW1(:) = 0.0
      ZZW2(:) = 0.0
      ZZW3(:) = 0.0
      ZZW4(:) = XACCR1 / ZLBDR(:)                 ! Mean volume drop diameter
      WHERE (GSCBU(:).AND.(ZZW4(:)>1.E-4))              ! analytical integration
         ZZW1(:) = XSCBU2 * ZCRT(:)**2 / ZLBDR3(:)   ! D>100 10-6 m
         ZZW3(:) = ZZW1(:)*ZSCBU(:)
      END WHERE
      WHERE (GSCBU(:).AND.(ZZW4(:)<=1.E-4))
         ZZW2(:) = XSCBU3 *(ZCRT(:) / ZLBDR3(:))**2  ! D<100 10-6 m
         ZZW3(:) = ZZW2(:)
      END WHERE
      ZCRS(:) = ZCRS(:) - MIN( ZCRS(:),ZZW3(:) * ZRHODREF(:) )
      DEALLOCATE(ZSCBU)
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'SCBU', &
                                             Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
   END IF
END IF ! LRAIN
!
!
!-------------------------------------------------------------------------------
!
!
!*       6. Unpack and clean
!   	 -------------------
!
!
   ZW(:,:,:) = PRCS(:,:,:)
   PRCS(:,:,:) = UNPACK( ZRCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PRRS(:,:,:)
   PRRS(:,:,:) = UNPACK( ZRRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PCCS(:,:,:)
   PCCS(:,:,:) = UNPACK( ZCCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
   ZW(:,:,:) = PCRS(:,:,:)
   PCRS(:,:,:) = UNPACK( ZCRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
   DEALLOCATE(ZRCT)
   DEALLOCATE(ZRRT)
   DEALLOCATE(ZCCT)
   DEALLOCATE(ZCRT)
   DEALLOCATE(ZRCS)
   DEALLOCATE(ZRRS)
   DEALLOCATE(ZCRS)
   DEALLOCATE(ZCCS)
   DEALLOCATE(ZRHODREF) 
   DEALLOCATE(GSELF)
   DEALLOCATE(GACCR)
   DEALLOCATE(GSCBU)
   IF( ALLOCATED(GENABLE_ACCR_SCBU) ) DEALLOCATE(GENABLE_ACCR_SCBU)
   DEALLOCATE(ZZW1)
   DEALLOCATE(ZZW2)
   DEALLOCATE(ZZW3)
   IF( ALLOCATED(ZZW4) ) DEALLOCATE(ZZW4)
   DEALLOCATE(ZLBDR3)
   DEALLOCATE(ZLBDC3)
   DEALLOCATE(ZLBDR)
   DEALLOCATE(ZLBDC)
END IF ! IMICRO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_WARM_COAL
