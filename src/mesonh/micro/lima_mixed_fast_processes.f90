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
      SUBROUTINE LIMA_MIXED_FAST_PROCESSES (ZRHODREF, ZZT, ZPRES, PTSTEP,           &
                                            ZLSFACT, ZLVFACT, ZKA, ZDV, ZCJ,        &
                                            ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT,     &
                                            ZRHT, ZCCT, ZCRT, ZCIT,                 &
                                            ZRCS, ZRRS, ZRIS, ZRSS, ZRGS, ZRHS,     &
                                            ZTHS, ZCCS, ZCRS, ZCIS,                 &
                                            ZLBDAC, ZLBDAR, ZLBDAS, ZLBDAG, ZLBDAH, &
                                            PRHODJ1D, GMICRO, PRHODJ, KMI, PTHS,    &
                                            PRCS, PRRS, PRIS, PRSS, PRGS, PRHS,     &
                                            PCCS, PCRS, PCIS                        )
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(:),   INTENT(IN)    :: ZZT       ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: ZPRES     ! Pressure
REAL,                 INTENT(IN)    :: PTSTEP    ! Time step          
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZKA       ! Thermal conductivity of the air
REAL, DIMENSION(:),   INTENT(IN)    :: ZDV       ! Diffusivity of water vapor in the air
REAL, DIMENSION(:),   INTENT(IN)    :: ZCJ       ! Ventilation coefficient ?
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRRT    ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRGT    ! Graupel m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRHT    ! Hail m.r. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZCCT    ! Cloud water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZCRT    ! Rain water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZCIT    ! Pristine ice conc. at t
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRRS    ! Rain water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRGS    ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRHS    ! Hail m.r. source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZTHS    ! Theta source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCCS    ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCRS    ! Rain water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCIS    ! Pristine ice conc. source
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAC  ! Slope param of the cloud droplet distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAR  ! Slope param of the raindrop  distr
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAS  ! Slope param of the aggregate distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAG  ! Slope param of the graupel distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAH  ! Slope param of the hail distr.
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
!
END SUBROUTINE LIMA_MIXED_FAST_PROCESSES
END INTERFACE
END MODULE MODI_LIMA_MIXED_FAST_PROCESSES
!
!     #######################################################################
      SUBROUTINE LIMA_MIXED_FAST_PROCESSES (ZRHODREF, ZZT, ZPRES, PTSTEP,           &
                                            ZLSFACT, ZLVFACT, ZKA, ZDV, ZCJ,        &
                                            ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT,     &
                                            ZRHT, ZCCT, ZCRT, ZCIT,                 &
                                            ZRCS, ZRRS, ZRIS, ZRSS, ZRGS, ZRHS,     &
                                            ZTHS, ZCCS, ZCRS, ZCIS,                 &
                                            ZLBDAC, ZLBDAR, ZLBDAS, ZLBDAG, ZLBDAH, &
                                            PRHODJ1D, GMICRO, PRHODJ, KMI, PTHS,    &
                                            PRCS, PRRS, PRIS, PRSS, PRGS, PRHS,     &
                                            PCCS, PCRS, PCIS                        )
!     #######################################################################
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

use mode_budget,           only: Budget_store_init, Budget_store_end

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(:),   INTENT(IN)    :: ZZT       ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: ZPRES     ! Pressure
REAL,                 INTENT(IN)    :: PTSTEP    ! Time step          
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZKA       ! Thermal conductivity of the air
REAL, DIMENSION(:),   INTENT(IN)    :: ZDV       ! Diffusivity of water vapor in the air
REAL, DIMENSION(:),   INTENT(IN)    :: ZCJ       ! Ventilation coefficient ?
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRRT    ! Rain water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRHT    ! Hail m.r. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZCCT    ! Cloud water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZCRT    ! Rain water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZCIT    ! Pristine ice conc. at t
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRRS    ! Rain water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRGS    ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRHS    ! Hail m.r. source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZTHS    ! Theta source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCCS    ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCRS    ! Rain water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCIS    ! Pristine ice conc. source
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAC  ! Slope param of the cloud droplet distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAR  ! Slope param of the raindrop  distr
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAS  ! Slope param of the aggregate distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAG  ! Slope param of the graupel distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAH  ! Slope param of the hail distr.
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

!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(ZZT)) :: GRIM, GACC, GDRY, GWET, GHAIL ! Test where to compute
INTEGER :: IGRIM, IGACC, IGDRY, IGWET, IHAIL
INTEGER :: JJ
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1,IVEC2        ! Vectors of indices
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1,ZVEC2, ZVEC3 ! Work vectors
REAL,    DIMENSION(SIZE(ZZT))  :: ZZW, ZZX      
REAL,    DIMENSION(SIZE(ZZT))  :: ZRDRYG, ZRWETG   
REAL,    DIMENSION(SIZE(ZZT),7)  :: ZZW1 
REAL :: NHAIL
REAL :: ZTHRH, ZTHRC
!
!-------------------------------------------------------------------------------
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
GRIM(:) = (ZRCT(:)>XRTMIN(2)) .AND. (ZRST(:)>XRTMIN(5)) .AND. (ZRCS(:)>XRTMIN(2)/PTSTEP) .AND. (ZZT(:)<XTT)
IGRIM = COUNT( GRIM(:) )
!
IF( IGRIM>0 ) THEN
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'RIM', &
                                            Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'RIM', &
                                            Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'RIM', &
                                            Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'RIM', &
                                            Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'RIM', &
                                            Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
  end if
!
!        1.1.0  allocations
!
   ALLOCATE(ZVEC1(IGRIM))
   ALLOCATE(ZVEC2(IGRIM))
   ALLOCATE(IVEC1(IGRIM))
   ALLOCATE(IVEC2(IGRIM))
!
!        1.1.1  select the ZLBDAS
!
   ZVEC1(:) = PACK( ZLBDAS(:),MASK=GRIM(:) )
!
!        1.1.2  find the next lower indice for the ZLBDAS in the geometrical
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
      ZZW1(:,1) = MIN( ZRCS(:),                              &
  	           XCRIMSS * ZZW(:) * ZRCT(:)                & ! RCRIMSS
  	                            *   ZLBDAS(:)**XEXCRIMSS &
    			            * ZRHODREF(:)**(-XCEXVT) )
      ZRCS(:) = ZRCS(:) - ZZW1(:,1)
      ZRSS(:) = ZRSS(:) + ZZW1(:,1)
      ZTHS(:) = ZTHS(:) + ZZW1(:,1)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RCRIMSS))
!
      ZCCS(:) = MAX( ZCCS(:)-ZZW1(:,1)*(ZCCT(:)/ZRCT(:)),0.0 ) ! Lambda_c**3
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
!
   WHERE ( GRIM(:) .AND. (ZRSS(:)>XRTMIN(5)/PTSTEP) )
      ZZW1(:,2) = MIN( ZRCS(:),                     &
    	           XCRIMSG * ZRCT(:)                & ! RCRIMSG
    	                   *  ZLBDAS(:)**XEXCRIMSG  &
  	                   * ZRHODREF(:)**(-XCEXVT) &
    		           - ZZW1(:,1)              )
      ZZW1(:,3) = MIN( ZRSS(:),                         &
                       XSRIMCG * ZLBDAS(:)**XEXSRIMCG   & ! RSRIMCG
   	                       * (1.0 - ZZW(:) )/(PTSTEP*ZRHODREF(:)))
      ZRCS(:) = ZRCS(:) - ZZW1(:,2)
      ZRSS(:) = ZRSS(:) - ZZW1(:,3)
      ZRGS(:) = ZRGS(:) + ZZW1(:,2) + ZZW1(:,3)
      ZTHS(:) = ZTHS(:) + ZZW1(:,2)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RCRIMSG))
!
      ZCCS(:) = MAX( ZCCS(:)-ZZW1(:,2)*(ZCCT(:)/ZRCT(:)),0.0 ) ! Lambda_c**3
   END WHERE
   DEALLOCATE(IVEC2)
   DEALLOCATE(IVEC1)
   DEALLOCATE(ZVEC2)
   DEALLOCATE(ZVEC1)

  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'RIM', &
                                           Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'RIM', &
                                           Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'RIM', &
                                           Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'RIM', &
                                           Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'RIM', &
                                           Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
  end if
END IF
!
!*       1.2  Hallett-Mossop ice multiplication process due to snow riming  
!        -----------------------------------------------------------------
!
!
GRIM(:) = (ZZT(:)<XHMTMAX) .AND. (ZZT(:)>XHMTMIN)                          &
                           .AND. (ZRST(:)>XRTMIN(5)) .AND. (ZRCT(:)>XRTMIN(2))
IGRIM = COUNT( GRIM(:) )
IF( IGRIM>0 ) THEN
  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'HMS', &
                                            Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'HMS', &
                                            Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HMS', &
                                            Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if

   ALLOCATE(ZVEC1(IGRIM))
   ALLOCATE(ZVEC2(IGRIM))
   ALLOCATE(IVEC2(IGRIM))
!
   ZVEC1(:) = PACK( ZLBDAC(:),MASK=GRIM(:) )
   ZVEC2(1:IGRIM) = MAX( 1.0001, MIN( REAL(NGAMINC)-0.0001,           &
                         XHMLINTP1 * LOG( ZVEC1(1:IGRIM) ) + XHMLINTP2 ) )
   IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
   ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - REAL( IVEC2(1:IGRIM) )
   ZVEC1(1:IGRIM) =   XGAMINC_HMC( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - XGAMINC_HMC( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
   ZZX(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 ) ! Large droplets
!
   WHERE ( GRIM(:) .AND. ZZX(:)<0.99 )
      ZZW1(:,5) = (ZZW1(:,1)+ZZW1(:,2))*(ZCCT(:)/ZRCT(:))*(1.0-ZZX(:))* & 
                                                             XHM_FACTS* &
           MAX( 0.0, MIN( (ZZT(:)-XHMTMIN)/3.0,(XHMTMAX-ZZT(:))/2.0 ) ) ! CCHMSI
      ZCIS(:) = ZCIS(:) + ZZW1(:,5)
!
      ZZW1(:,6) = ZZW1(:,5) * XMNU0                                     ! RCHMSI
      ZRIS(:) = ZRIS(:) + ZZW1(:,6)
      ZRSS(:) = ZRSS(:) - ZZW1(:,6)
   END WHERE
   DEALLOCATE(IVEC2)
   DEALLOCATE(ZVEC2)
   DEALLOCATE(ZVEC1)

  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'HMS', &
                                          Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'HMS', &
                                          Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HMS', &
                                          Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if
END IF
!
!
!*       1.3  Rain accretion onto the aggregates  
!        ---------------------------------------
!
!
ZZW1(:,2:3) = 0.0
GACC(:) = (ZRRT(:)>XRTMIN(3)) .AND. (ZRST(:)>XRTMIN(5)) .AND. (ZRRS(:)>XRTMIN(3)/PTSTEP) .AND. (ZZT(:)<XTT)
IGACC = COUNT( GACC(:) )
!
IF( IGACC>0 .AND. LRAIN) THEN
  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'ACC', &
                                            Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'ACC', &
                                            Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'ACC', &
                                            Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'ACC', &
                                            Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'ACC', &
                                            Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
  end if
!
!        1.3.0  allocations
!
   ALLOCATE(ZVEC1(IGACC))
   ALLOCATE(ZVEC2(IGACC))
   ALLOCATE(ZVEC3(IGACC))
   ALLOCATE(IVEC1(IGACC))
   ALLOCATE(IVEC2(IGACC))
!
!        1.3.1  select the (ZLBDAS,ZLBDAR) couplet
!
   ZVEC1(:) = PACK( ZLBDAS(:),MASK=GACC(:) )
   ZVEC2(:) = PACK( ZLBDAR(:),MASK=GACC(:) )
!
!        1.3.2  find the next lower indice for the ZLBDAS and for the ZLBDAR
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
!        1.3.3  perform the bilinear interpolation of the normalized
!               RACCSS-kernel
!
   DO JJ = 1,IGACC
      ZVEC3(JJ) =  (  XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                				 	           * ZVEC1(JJ) &
                 - (  XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
  	                    			             * (ZVEC1(JJ) - 1.0)
   END DO
   ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
!        1.3.4  raindrop accretion on the small sized aggregates
!
   WHERE ( GACC(:) )
      ZZW1(:,2) = ZCRT(:) *                                           & !! coef of RRACCS
              XFRACCSS*( ZLBDAS(:)**XCXS )*( ZRHODREF(:)**(-XCEXVT-1.) ) &
         *( XLBRACCS1/((ZLBDAS(:)**2)               ) +                  &
            XLBRACCS2/( ZLBDAS(:)    * ZLBDAR(:)    ) +                  &
            XLBRACCS3/(               (ZLBDAR(:)**2)) )/ZLBDAR(:)**3
      ZZW1(:,4) = MIN( ZRRS(:),ZZW1(:,2)*ZZW(:) )           ! RRACCSS
      ZRRS(:) = ZRRS(:) - ZZW1(:,4)
      ZRSS(:) = ZRSS(:) + ZZW1(:,4)
      ZTHS(:) = ZTHS(:) + ZZW1(:,4)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RRACCSS))
!
      ZCRS(:) = MAX( ZCRS(:)-ZZW1(:,4)*(ZCRT(:)/ZRRT(:)),0.0 ) ! Lambda_r**3 
   END WHERE
!
!        1.3.4b perform the bilinear interpolation of the normalized
!               RACCS-kernel
!
   DO JJ = 1,IGACC
      ZVEC3(JJ) =  (   XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    -  XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                         * ZVEC1(JJ) &
                 - (   XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    -  XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                           * (ZVEC1(JJ) - 1.0)
   END DO
   ZZW1(:,2) = ZZW1(:,2)*UNPACK( VECTOR=ZVEC3(:),MASK=GACC(:),FIELD=0.0 ) !! RRACCS
!
!        1.3.5  perform the bilinear interpolation of the normalized
!               SACCRG-kernel
!
   DO JJ = 1,IGACC
      ZVEC3(JJ) =  (  XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    - XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
      			 	                                   * ZVEC2(JJ) &
                 - (  XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                    - XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
			                                     * (ZVEC2(JJ) - 1.0)
   END DO
   ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
!        1.3.6  raindrop accretion-conversion of the large sized aggregates
!               into graupeln
!
   WHERE ( GACC(:) .AND. (ZRSS(:)>XRTMIN(5)/PTSTEP) )
      ZZW1(:,2) = MAX( MIN( ZRRS(:),ZZW1(:,2)-ZZW1(:,4) ) , 0. )      ! RRACCSG
      ZZW1(:,3) = MIN( ZRSS(:),XFSACCRG*ZZW(:)*                     & ! RSACCRG
            ( ZLBDAS(:)**(XCXS-XBS) )*( ZRHODREF(:)**(-XCEXVT-1.) ) &
           *( XLBSACCR1/((ZLBDAR(:)**2)               ) +           &
              XLBSACCR2/( ZLBDAR(:)    * ZLBDAS(:)    ) +           &
              XLBSACCR3/(               (ZLBDAS(:)**2)) ) )
      ZRRS(:) = ZRRS(:) - ZZW1(:,2)
      ZRSS(:) = ZRSS(:) - ZZW1(:,3)
      ZRGS(:) = ZRGS(:) + ZZW1(:,2)+ZZW1(:,3)
      ZTHS(:) = ZTHS(:) + ZZW1(:,2)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RRACCSG))
!
      ZCRS(:) = MAX( ZCRS(:)-ZZW1(:,2)*(ZCRT(:)/ZRRT(:)),0.0 ) ! Lambda_r**3 
   END WHERE
   DEALLOCATE(IVEC2)
   DEALLOCATE(IVEC1)
   DEALLOCATE(ZVEC3)
   DEALLOCATE(ZVEC2)
   DEALLOCATE(ZVEC1)

  ! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'ACC', &
                                           Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'ACC', &
                                           Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'ACC', &
                                           Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'ACC', &
                                           Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'ACC', &
                                           Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
  end if
END IF
!
!*       1.4  Conversion-Melting of the aggregates
!        -----------------------------------------
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'CMEL', &
                                          Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'CMEL', &
                                          Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
end if

ZZW(:) = 0.0
WHERE( (ZRST(:)>XRTMIN(5)) .AND. (ZRSS(:)>XRTMIN(5)/PTSTEP) .AND. (ZZT(:)>XTT) )
   ZZW(:) = ZRVT(:)*ZPRES(:)/((XMV/XMD)+ZRVT(:)) ! Vapor pressure
   ZZW(:) =  ZKA(:)*(XTT-ZZT(:)) +                                 &
              ( ZDV(:)*(XLVTT + ( XCPV - XCL ) * ( ZZT(:) - XTT )) &
                          *(XESTT-ZZW(:))/(XRV*ZZT(:))             )
!
! compute RSMLT
!
   ZZW(:)  = MIN( ZRSS(:), XFSCVMG*MAX( 0.0,( -ZZW(:) *             &
                          ( X0DEPS*       ZLBDAS(:)**XEX0DEPS +     &
                            X1DEPS*ZCJ(:)*ZLBDAS(:)**XEX1DEPS ) -   &
                                    ( ZZW1(:,1)+ZZW1(:,4) ) *       &
                             ( ZRHODREF(:)*XCL*(XTT-ZZT(:))) ) /    &
                                            ( ZRHODREF(:)*XLMTT ) ) )
!
! note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
! because the graupeln produced by this process are still icy!!!
!
   ZRSS(:) = ZRSS(:) - ZZW(:)
   ZRGS(:) = ZRGS(:) + ZZW(:)
END WHERE
!
! Budget storage
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'CMEL', &
                                         Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'CMEL', &
                                         Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
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
                                         Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'CFRZ', &
                                         Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'CFRZ', &
                                         Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'CFRZ', &
                                         Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'CFRZ', &
                                         Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CFRZ', &
                                         Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
end if

ZZW1(:,3:4) = 0.0
WHERE( (ZRIT(:)>XRTMIN(4)) .AND. (ZRRT(:)>XRTMIN(3)) .AND. (ZRIS(:)>XRTMIN(4)/PTSTEP) .AND. (ZRRS(:)>XRTMIN(3)/PTSTEP) )
   ZZW1(:,3) = MIN( ZRIS(:),XICFRR * ZRIT(:) * ZCRT(:)          & ! RICFRRG
                                   * ZLBDAR(:)**XEXICFRR        &
                                   * ZRHODREF(:)**(-XCEXVT-1.0) )
!
   ZZW1(:,4) = MIN( ZRRS(:),XRCFRI * ZCIT(:) * ZCRT(:)          & ! RRCFRIG
                                   * ZLBDAR(:)**XEXRCFRI        &
                                   * ZRHODREF(:)**(-XCEXVT-2.0) )
   ZRIS(:) = ZRIS(:) - ZZW1(:,3)
   ZRRS(:) = ZRRS(:) - ZZW1(:,4)
   ZRGS(:) = ZRGS(:) + ZZW1(:,3)+ZZW1(:,4)
   ZTHS(:) = ZTHS(:) + ZZW1(:,4)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*RRCFRIG)
!
   ZCIS(:) = MAX( ZCIS(:)-ZZW1(:,3)*(ZCIT(:)/ZRIT(:)),0.0 )     ! CICFRRG
   ZCRS(:) = MAX( ZCRS(:)-ZZW1(:,4)*(ZCRT(:)/ZRRT(:)),0.0 )     ! CRCFRIG
END WHERE

if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'CFRZ', &
                                         Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'CFRZ', &
                                         Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'CFRZ', &
                                         Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'CFRZ', &
                                         Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'CFRZ', &
                                         Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CFRZ', &
                                         Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
end if
!
!*       2.2  Compute the Dry growth case
!        --------------------------------
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'WETG', &
                                          Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'WETG', &
                                          Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'WETG', &
                                          Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'WETG', &
                                          Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'WETG', &
                                          Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'WETG', &
                                          Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'WETG', &
                                          Unpack( zrhs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'WETG', &
                                          Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'WETG', &
                                          Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'WETG', &
                                          Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if
end if
!
ZZW1(:,:) = 0.0
WHERE( ((ZRCT(:)>XRTMIN(2)) .AND. (ZRGT(:)>XRTMIN(6)) .AND. (ZRCS(:)>XRTMIN(2)/PTSTEP)) .OR. &
       ((ZRIT(:)>XRTMIN(4)) .AND. (ZRGT(:)>XRTMIN(6)) .AND. (ZRIS(:)>XRTMIN(4)/PTSTEP))      )
   ZZW(:) = ZLBDAG(:)**(XCXG-XDG-2.0) * ZRHODREF(:)**(-XCEXVT)
   ZZW1(:,1) = MIN( ZRCS(:),XFCDRYG * ZRCT(:) * ZZW(:) )             ! RCDRYG
   ZZW1(:,2) = MIN( ZRIS(:),XFIDRYG * EXP( XCOLEXIG*(ZZT(:)-XTT) ) &
                                    * ZRIT(:) * ZZW(:) )             ! RIDRYG
END WHERE
!
!*       2.2.1  accretion of aggregates on the graupeln
!        ----------------------------------------------
!
GDRY(:) = (ZRST(:)>XRTMIN(5)) .AND. (ZRGT(:)>XRTMIN(6)) .AND. (ZRSS(:)>XRTMIN(5)/PTSTEP)
IGDRY = COUNT( GDRY(:) )
!
IF( IGDRY>0 ) THEN
!
!*       2.2.2  allocations
!
   ALLOCATE(ZVEC1(IGDRY))
   ALLOCATE(ZVEC2(IGDRY))
   ALLOCATE(ZVEC3(IGDRY))
   ALLOCATE(IVEC1(IGDRY))
   ALLOCATE(IVEC2(IGDRY))
!
!*       2.2.3  select the (ZLBDAG,ZLBDAS) couplet
!
   ZVEC1(:) = PACK( ZLBDAG(:),MASK=GDRY(:) )
   ZVEC2(:) = PACK( ZLBDAS(:),MASK=GDRY(:) )
!
!*       2.2.4  find the next lower indice for the ZLBDAG and for the ZLBDAS
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
!*       2.2.5  perform the bilinear interpolation of the normalized
!               SDRYG-kernel
!
   DO JJ = 1,IGDRY
      ZVEC3(JJ) =  (  XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
      			 	                                  * ZVEC1(JJ) &
                 - (  XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
       			                                    * (ZVEC1(JJ) - 1.0)
   END DO
   ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
   WHERE( GDRY(:) )
      ZZW1(:,3) = MIN( ZRSS(:),XFSDRYG*ZZW(:)                         & ! RSDRYG
                                      * EXP( XCOLEXSG*(ZZT(:)-XTT) )  &
                    *( ZLBDAS(:)**(XCXS-XBS) )*( ZLBDAG(:)**XCXG )    &
                    *( ZRHODREF(:)**(-XCEXVT-1.) )                    &
                         *( XLBSDRYG1/( ZLBDAG(:)**2              ) + &
                            XLBSDRYG2/( ZLBDAG(:)   * ZLBDAS(:)   ) + &
                            XLBSDRYG3/(               ZLBDAS(:)**2) ) )
   END WHERE
   DEALLOCATE(IVEC2)
   DEALLOCATE(IVEC1)
   DEALLOCATE(ZVEC3)
   DEALLOCATE(ZVEC2)
   DEALLOCATE(ZVEC1)
END IF
!
!*       2.2.6  accretion of raindrops on the graupeln
!        ---------------------------------------------
!
GDRY(:) = (ZRRT(:)>XRTMIN(3)) .AND. (ZRGT(:)>XRTMIN(6)) .AND. (ZRRS(:)>XRTMIN(3))
IGDRY = COUNT( GDRY(:) )
!
IF( IGDRY>0 ) THEN
!
!*       2.2.7  allocations
!
   ALLOCATE(ZVEC1(IGDRY))
   ALLOCATE(ZVEC2(IGDRY))
   ALLOCATE(ZVEC3(IGDRY))
   ALLOCATE(IVEC1(IGDRY))
   ALLOCATE(IVEC2(IGDRY))
!
!*       2.2.8  select the (ZLBDAG,ZLBDAR) couplet
!
   ZVEC1(:) = PACK( ZLBDAG(:),MASK=GDRY(:) )
   ZVEC2(:) = PACK( ZLBDAR(:),MASK=GDRY(:) )
!
!*       2.2.9  find the next lower indice for the ZLBDAG and for the ZLBDAR
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
!*       2.2.10 perform the bilinear interpolation of the normalized
!               RDRYG-kernel
!
   DO JJ = 1,IGDRY
      ZVEC3(JJ) =  (  XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                     			 	                  * ZVEC1(JJ) &
                 - (  XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                 			     * (ZVEC1(JJ) - 1.0)
   END DO
   ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
   WHERE( GDRY(:) )
      ZZW1(:,4) = MIN( ZRRS(:),XFRDRYG*ZZW(:) * ZCRT(:)                   & ! RRDRYG
                        *( ZLBDAR(:)**(-3) )*( ZLBDAG(:)**XCXG ) &
                                *( ZRHODREF(:)**(-XCEXVT-1.) )   &
                    *( XLBRDRYG1/( ZLBDAG(:)**2              ) + &
                       XLBRDRYG2/( ZLBDAG(:)   * ZLBDAR(:)   ) + &
                       XLBRDRYG3/(               ZLBDAR(:)**2) ) )
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
!*       2.3  Compute the Wet growth case
!        --------------------------------
!
!
ZZW(:) = 0.0
ZRWETG(:) = 0.0
WHERE( ZRGT(:)>XRTMIN(6) )
   ZZW1(:,5) = MIN( ZRIS(:),                                    &
               ZZW1(:,2) / (XCOLIG*EXP(XCOLEXIG*(ZZT(:)-XTT)) ) ) ! RIWETG
   ZZW1(:,6) = MIN( ZRSS(:),                                    &
               ZZW1(:,3) / (XCOLSG*EXP(XCOLEXSG*(ZZT(:)-XTT)) ) ) ! RSWETG
!
   ZZW(:) = ZRVT(:)*ZPRES(:)/((XMV/XMD)+ZRVT(:)) ! Vapor pressure
   ZZW(:) =  ZKA(:)*(XTT-ZZT(:)) +                                  &
             ( ZDV(:)*(XLVTT + ( XCPV - XCL ) * ( ZZT(:) - XTT ))   &
                           *(XESTT-ZZW(:))/(XRV*ZZT(:))             )
!
! compute RWETG
!
   ZRWETG(:)  = MAX( 0.0,                                               &
                   ( ZZW(:) * ( X0DEPG*       ZLBDAG(:)**XEX0DEPG +     &
                                X1DEPG*ZCJ(:)*ZLBDAG(:)**XEX1DEPG ) +   &
                   ( ZZW1(:,5)+ZZW1(:,6) ) *                            &
                   ( ZRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-ZZT(:)))   ) ) / &
                                   ( ZRHODREF(:)*(XLMTT-XCL*(XTT-ZZT(:))) )   )
END WHERE
!
!
!*       2.4  Select Wet or Dry case
!        ---------------------------
!
!
! Wet case and partial conversion to hail
!
ZZW(:) = 0.0
NHAIL = 0.
IF (LHAIL) NHAIL = 1. 
WHERE( ZRGT(:)>XRTMIN(6) .AND. ZZT(:)<XTT                               &
                         .AND. ZRDRYG(:)>=ZRWETG(:) .AND. ZRWETG(:)>0.0 ) 
!   
   ZZW(:) = ZRWETG(:) - ZZW1(:,5) - ZZW1(:,6) ! RCWETG+RRWETG
!   
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!   
   ZZW1(:,7) = MAX( 0.0,MIN( ZZW(:),ZRRS(:)+ZZW1(:,1) ) )
   ZZX(:)    = ZZW1(:,7) / ZZW(:)
   ZZW1(:,5) = ZZW1(:,5)*ZZX(:)
   ZZW1(:,6) = ZZW1(:,6)*ZZX(:)
   ZRWETG(:) = ZZW1(:,7) + ZZW1(:,5) + ZZW1(:,6)
!   
   ZRCS(:) = ZRCS(:) - ZZW1(:,1)
   ZRIS(:) = ZRIS(:) - ZZW1(:,5)
   ZRSS(:) = ZRSS(:) - ZZW1(:,6)
!
! assume a linear percent of conversion of graupel into hail
!
   ZRGS(:) = ZRGS(:) + ZRWETG(:)
   ZZW(:)  = ZRGS(:)*ZRDRYG(:)*NHAIL/(ZRWETG(:)+ZRDRYG(:)) 
   ZRGS(:) = ZRGS(:) - ZZW(:)                        
   ZRHS(:) = ZRHS(:) + ZZW(:)
   ZRRS(:) = MAX( 0.0,ZRRS(:) - ZZW1(:,7) + ZZW1(:,1) )
   ZTHS(:) = ZTHS(:) + ZZW1(:,7)*(ZLSFACT(:)-ZLVFACT(:))
                                                 ! f(L_f*(RCWETG+RRWETG))
!
   ZCCS(:) = MAX( ZCCS(:)-ZZW1(:,1)*(ZCCT(:)/MAX(ZRCT(:),XRTMIN(2))),0.0 )
   ZCIS(:) = MAX( ZCIS(:)-ZZW1(:,5)*(ZCIT(:)/MAX(ZRIT(:),XRTMIN(4))),0.0 )
   ZCRS(:) = MAX( ZCRS(:)-MAX( ZZW1(:,7)-ZZW1(:,1),0.0 )                 &
				   *(ZCRT(:)/MAX(ZRRT(:),XRTMIN(3))),0.0 )
END WHERE
!
! Budget storage
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'WETG', &
                                         Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'WETG', &
                                         Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'WETG', &
                                         Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'WETG', &
                                         Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'WETG', &
                                         Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'WETG', &
                                         Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'WETG', &
                                         Unpack( zrhs(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'WETG', &
                                         Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'WETG', &
                                         Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'WETG', &
                                         Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if
end if
!
! Dry case
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'DRYG', &
                                          Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'DRYG', &
                                          Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'DRYG', &
                                          Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'DRYG', &
                                          Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'DRYG', &
                                          Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'DRYG', &
                                          Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'DRYG', &
                                          Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'DRYG', &
                                          Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'DRYG', &
                                          Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if
end if

WHERE( ZRGT(:)>XRTMIN(6) .AND. ZZT(:)<XTT                              &
                         .AND. ZRDRYG(:)<ZRWETG(:) .AND. ZRDRYG(:)>0.0 ) ! case
   ZRCS(:) = ZRCS(:) - ZZW1(:,1)
   ZRIS(:) = ZRIS(:) - ZZW1(:,2)
   ZRSS(:) = ZRSS(:) - ZZW1(:,3)
   ZRRS(:) = ZRRS(:) - ZZW1(:,4)
   ZRGS(:) = ZRGS(:) + ZRDRYG(:)
   ZTHS(:) = ZTHS(:) + (ZZW1(:,1)+ZZW1(:,4))*(ZLSFACT(:)-ZLVFACT(:)) !
  						        ! f(L_f*(RCDRYG+RRDRYG))
!
   ZCCS(:) = MAX( ZCCS(:)-ZZW1(:,1)*(ZCCT(:)/MAX(ZRCT(:),XRTMIN(2))),0.0 )
   ZCIS(:) = MAX( ZCIS(:)-ZZW1(:,2)*(ZCIT(:)/MAX(ZRIT(:),XRTMIN(4))),0.0 )
   ZCRS(:) = MAX( ZCRS(:)-ZZW1(:,4)*(ZCRT(:)/MAX(ZRRT(:),XRTMIN(3))),0.0 ) 
                                                         ! Approximate rates
END WHERE
!
! Budget storage
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'DRYG', &
                                         Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'DRYG', &
                                         Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'DRYG', &
                                         Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'DRYG', &
                                         Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'DRYG', &
                                         Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'DRYG', &
                                         Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) then
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'DRYG', &
                                         Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'DRYG', &
                                         Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'DRYG', &
                                         Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
  end if
end if
!
!
!*       2.5  Hallett-Mossop ice multiplication process due to graupel riming
!        --------------------------------------------------------------------
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'HMG', &
                                          Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'HMG', &
                                          Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HMG', &
                                          Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
end if

GDRY(:) = (ZZT(:)<XHMTMAX) .AND. (ZZT(:)>XHMTMIN)    .AND. (ZRDRYG(:)<ZZW(:))&
                           .AND. (ZRGT(:)>XRTMIN(6)) .AND. (ZRCT(:)>XRTMIN(2))
IGDRY = COUNT( GDRY(:) )
IF( IGDRY>0 ) THEN
   ALLOCATE(ZVEC1(IGDRY))
   ALLOCATE(ZVEC2(IGDRY))
   ALLOCATE(IVEC2(IGDRY))
!
   ZVEC1(:) = PACK( ZLBDAC(:),MASK=GDRY(:) )
   ZVEC2(1:IGDRY) = MAX( 1.0001, MIN( REAL(NGAMINC)-0.0001,           &
                         XHMLINTP1 * LOG( ZVEC1(1:IGDRY) ) + XHMLINTP2 ) )
   IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
   ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
   ZVEC1(1:IGDRY) =   XGAMINC_HMC( IVEC2(1:IGDRY)+1 )* ZVEC2(1:IGDRY)      &
                    - XGAMINC_HMC( IVEC2(1:IGDRY)   )*(ZVEC2(1:IGDRY) - 1.0)
   ZZX(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GDRY,FIELD=0.0 ) ! Large droplets
!
   WHERE ( GDRY(:) .AND. ZZX(:)<0.99 ) ! Dry case
      ZZW1(:,5) = ZZW1(:,1)*(ZCCT(:)/ZRCT(:))*(1.0-ZZX(:))*XHM_FACTG*  &
           MAX( 0.0, MIN( (ZZT(:)-XHMTMIN)/3.0,(XHMTMAX-ZZT(:))/2.0 ) ) ! CCHMGI
      ZCIS(:) = ZCIS(:) + ZZW1(:,5)
!
      ZZW1(:,6) = ZZW1(:,5) * XMNU0                                     ! RCHMGI
      ZRIS(:) = ZRIS(:) + ZZW1(:,6)
      ZRGS(:) = ZRGS(:) - ZZW1(:,6)
   END WHERE
   DEALLOCATE(IVEC2)
   DEALLOCATE(ZVEC2)
   DEALLOCATE(ZVEC1)
END IF
!
! Budget storage
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'HMG', &
                                         Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'HMG', &
                                         Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'HMG', &
                                         Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
end if
!
!*       2.6  Melting of the graupeln
!        ----------------------------
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'GMLT', &
                                          Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'GMLT', &
                                          Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'GMLT', &
                                          Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'GMLT', &
                                          Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
end if

ZZW(:) = 0.0
WHERE( (ZRGT(:)>XRTMIN(6)) .AND. (ZRGS(:)>XRTMIN(6)/PTSTEP) .AND. (ZZT(:)>XTT) )
   ZZW(:) = ZRVT(:)*ZPRES(:)/((XMV/XMD)+ZRVT(:)) ! Vapor pressure
   ZZW(:) =  ZKA(:)*(XTT-ZZT(:)) +                                 &
              ( ZDV(:)*(XLVTT + ( XCPV - XCL ) * ( ZZT(:) - XTT )) &
                          *(XESTT-ZZW(:))/(XRV*ZZT(:))             )
!
! compute RGMLTR
!
   ZZW(:)  = MIN( ZRGS(:), MAX( 0.0,( -ZZW(:) *                     &
                          ( X0DEPG*       ZLBDAG(:)**XEX0DEPG +     &
                            X1DEPG*ZCJ(:)*ZLBDAG(:)**XEX1DEPG ) -   &
                                    ( ZZW1(:,1)+ZZW1(:,4) ) *       &
                             ( ZRHODREF(:)*XCL*(XTT-ZZT(:))) ) /    &
                                            ( ZRHODREF(:)*XLMTT ) ) )
   ZRRS(:) = ZRRS(:) + ZZW(:)
   ZRGS(:) = ZRGS(:) - ZZW(:)
   ZTHS(:) = ZTHS(:) - ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(-RGMLTR))
!
!   ZCRS(:) = MAX( ZCRS(:) + ZZW(:)*(XCCG*ZLBDAG(:)**XCXG/ZRGT(:)),0.0 )
   ZCRS(:) = ZCRS(:) + ZZW(:)*5.0E6  ! obtained after averaging
                                     ! Dshed=1mm and 500 microns
END WHERE
!
if ( nbumod == kmi .and. lbu_enable ) then
  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'GMLT', &
                                         Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'GMLT', &
                                         Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'GMLT', &
                                         Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'GMLT', &
                                         Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
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
GHAIL(:) = ZRHT(:)>XRTMIN(7)
IHAIL = COUNT(GHAIL(:))
!
IF( IHAIL>0 ) THEN
!
!*       3.1 Wet growth of hail 
!        ----------------------------
!
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'WETH', &
                                            Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'WETH', &
                                            Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'WETH', &
                                            Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'WETH', &
                                            Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'WETH', &
                                            Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_init( tbudgets(NBUDGET_RG), 'WETH', &
                                            Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'WETH', &
                                            Unpack( zrhs(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
                    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'WETH', &
                                            Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'WETH', &
                                            Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                    call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'WETH', &
                                            Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
    end if
  end if

   ZZW1(:,:) = 0.0
   WHERE( GHAIL(:) .AND. ( (ZRCT(:)>XRTMIN(2) .AND. ZRCS(:)>XRTMIN(2)/PTSTEP) .OR. &
                           (ZRIT(:)>XRTMIN(4) .AND. ZRIS(:)>XRTMIN(4)/PTSTEP) )    )    
      ZZW(:) = ZLBDAH(:)**(XCXH-XDH-2.0) * ZRHODREF(:)**(-XCEXVT)
      ZZW1(:,1) = MIN( ZRCS(:),XFWETH * ZRCT(:) * ZZW(:) )             ! RCWETH
      ZZW1(:,2) = MIN( ZRIS(:),XFWETH * ZRIT(:) * ZZW(:) )             ! RIWETH
   END WHERE
!
!*       3.1.1  accretion of aggregates on the hailstones
!        ------------------------------------------------
!
   GWET(:) = GHAIL(:) .AND. (ZRST(:)>XRTMIN(5) .AND. ZRSS(:)>XRTMIN(5)/PTSTEP)
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
!*       3.1.3  select the (ZLBDAH,ZLBDAS) couplet
!
      ZVEC1(:) = PACK( ZLBDAH(:),MASK=GWET(:) )
      ZVEC2(:) = PACK( ZLBDAS(:),MASK=GWET(:) )
!
!*       3.1.4  find the next lower indice for the ZLBDAG and for the ZLBDAS
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
        ZVEC3(JJ) = (  XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
        		 	                                   * ZVEC1(JJ) &
                   - ( XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
         		                                     * (ZVEC1(JJ) - 1.0)
      END DO
      ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
      WHERE( GWET(:) )
        ZZW1(:,3) = MIN( ZRSS(:),XFSWETH*ZZW(:)                       & ! RSWETH
                      *( ZLBDAS(:)**(XCXS-XBS) )*( ZLBDAH(:)**XCXH )  &
       	                 *( ZRHODREF(:)**(-XCEXVT-1.) )               &
                         *( XLBSWETH1/( ZLBDAH(:)**2              ) + &
                            XLBSWETH2/( ZLBDAH(:)   * ZLBDAS(:)   ) + &
                            XLBSWETH3/(               ZLBDAS(:)**2) ) )
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
    GWET(:) = GHAIL(:) .AND. (ZRGT(:)>XRTMIN(6) .AND. ZRGS(:)>XRTMIN(6)/PTSTEP)
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
!*       3.1.8  select the (ZLBDAH,ZLBDAG) couplet
!
      ZVEC1(:) = PACK( ZLBDAH(:),MASK=GWET(:) )
      ZVEC2(:) = PACK( ZLBDAG(:),MASK=GWET(:) )
!
!*       3.1.9  find the next lower indice for the ZLBDAH and for the ZLBDAG
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
        ZVEC3(JJ) = (  XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                       			 	                   * ZVEC1(JJ) &
                  - (  XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                     - XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                   			     * (ZVEC1(JJ) - 1.0)
      END DO
      ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GWET,FIELD=0.0 )
!
      WHERE( GWET(:) )
        ZZW1(:,5) = MAX(MIN( ZRGS(:),XFGWETH*ZZW(:)                       & ! RGWETH
                      *( ZLBDAG(:)**(XCXG-XBG) )*( ZLBDAH(:)**XCXH )  &
                         *( ZRHODREF(:)**(-XCEXVT-1.) )               &
                         *( XLBGWETH1/( ZLBDAH(:)**2              ) + &
                            XLBGWETH2/( ZLBDAH(:)   * ZLBDAG(:)   ) + &
                            XLBGWETH3/(               ZLBDAG(:)**2) ) ),0. )
      END WHERE
      DEALLOCATE(IVEC2)
      DEALLOCATE(IVEC1)
      DEALLOCATE(ZVEC3)
      DEALLOCATE(ZVEC2)
      DEALLOCATE(ZVEC1)
   END IF
!
!*       3.2    compute the Wet growth of hail
!        -------------------------------------
!
    ZZW(:) = 0.0
    WHERE( GHAIL(:) .AND. ZZT(:)<XTT )
       ZZW(:) = ZRVT(:)*ZPRES(:)/((XMV/XMD)+ZRVT(:)) ! Vapor pressure
       ZZW(:) = ZKA(:)*(XTT-ZZT(:)) +                                 &
                 ( ZDV(:)*(XLVTT + ( XCPV - XCL ) * ( ZZT(:) - XTT )) &
                             *(XESTT-ZZW(:))/(XRV*ZZT(:))             )
!
! compute RWETH
!
       ZZW(:)  =  MAX(0.,  ( ZZW(:) * ( X0DEPH*       ZLBDAH(:)**XEX0DEPH +     &
                                 X1DEPH*ZCJ(:)*ZLBDAH(:)**XEX1DEPH ) +   &
                    ( ZZW1(:,2)+ZZW1(:,3)+ZZW1(:,5) ) *                  &
                    ( ZRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-ZZT(:)))   ) ) / &
                          ( ZRHODREF(:)*(XLMTT-XCL*(XTT-ZZT(:))) ) )
!
       ZZW1(:,6) = MAX( ZZW(:) - ZZW1(:,2) - ZZW1(:,3) - ZZW1(:,5),0.) ! RCWETH+RRWETH
    END WHERE
    WHERE ( GHAIL(:) .AND. ZZT(:)<XTT  .AND. ZZW1(:,6)/=0.)
!
! limitation of the available rainwater mixing ratio (RRWETH < RRS !)
!
       ZZW1(:,4) = MAX( 0.0,MIN( ZZW1(:,6),ZRRS(:)+ZZW1(:,1) ) )
       ZZX(:)    = ZZW1(:,4) / ZZW1(:,6)
       ZZW1(:,2) = ZZW1(:,2)*ZZX(:)
       ZZW1(:,3) = ZZW1(:,3)*ZZX(:)
       ZZW1(:,5) = ZZW1(:,5)*ZZX(:)
       ZZW(:)    = ZZW1(:,4) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,5)
!
!*       3.2.1  integrate the Wet growth of hail
!
       ZRCS(:) = ZRCS(:) - ZZW1(:,1)
       ZRIS(:) = ZRIS(:) - ZZW1(:,2)
       ZRSS(:) = ZRSS(:) - ZZW1(:,3)
       ZRGS(:) = ZRGS(:) - ZZW1(:,5)
       ZRHS(:) = ZRHS(:) + ZZW(:)
       ZRRS(:) = MAX( 0.0,ZRRS(:) - ZZW1(:,4) + ZZW1(:,1) )
       ZTHS(:) = ZTHS(:) + ZZW1(:,4)*(ZLSFACT(:)-ZLVFACT(:)) 
       						        ! f(L_f*(RCWETH+RRWETH))
!
       ZCCS(:) = MAX( ZCCS(:)-ZZW1(:,1)*(ZCCT(:)/MAX(ZRCT(:),XRTMIN(2))),0.0 )
       ZCIS(:) = MAX( ZCIS(:)-ZZW1(:,2)*(ZCIT(:)/MAX(ZRIT(:),XRTMIN(4))),0.0 )
       ZCRS(:) = MAX( ZCRS(:)-MAX( ZZW1(:,4)-ZZW1(:,1),0.0 )                 &
                                       *(ZCRT(:)/MAX(ZRRT(:),XRTMIN(3))),0.0 )
    END WHERE

  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'WETH', &
                                           Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'WETH', &
                                           Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'WETH', &
                                           Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'WETH', &
                                           Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'WETH', &
                                           Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'WETH', &
                                           Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'WETH', &
                                           Unpack( zrhs(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
                    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'WETH', &
                                           Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
                    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'WETH', &
                                           Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
                    call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'WETH', &
                                           Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
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
                                            Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'COHG', &
                                            Unpack( zrhs(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
  end if

    ZTHRH=0.01E-3
    ZTHRC=0.001E-3
    ZZW(:) = 0.0
    WHERE( ZRHT(:)<ZTHRH .AND. ZRCT(:)<ZTHRC .AND. ZZT(:)<XTT )
       ZZW(:) = MIN( 1.0,MAX( 0.0,1.0-(ZRCT(:)/ZTHRC) ) )
!
! assume a linear percent conversion rate of hail into graupel
!
       ZZW(:)  = ZRHS(:)*ZZW(:)
       ZRGS(:) = ZRGS(:) + ZZW(:)                      !   partial conversion
       ZRHS(:) = ZRHS(:) - ZZW(:)                      ! of hail into graupel
!
    END WHERE

  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_rg ) call Budget_store_end( tbudgets(NBUDGET_RG), 'COHG', &
                                           Unpack( zrgs(:), mask = gmicro(:, :, :), field = prgs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'COHG', &
                                           Unpack( zrhs(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
  end if
END IF
!
!*       3.4    Melting of the hailstones
!
IF ( IHAIL>0 ) THEN
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HMLT', &
                                            Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'HMLT', &
                                            Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_init( tbudgets(NBUDGET_RH), 'HMLT', &
                                            Unpack( zrhs(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'HMLT', &
                                            Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
  end if

    ZZW(:) = 0.0
    WHERE( GHAIL(:) .AND. (ZRHS(:)>XRTMIN(7)/PTSTEP) .AND. (ZRHT(:)>XRTMIN(7)) .AND. (ZZT(:)>XTT) )
       ZZW(:) = ZRVT(:)*ZPRES(:)/((XMV/XMD)+ZRVT(:)) ! Vapor pressure
       ZZW(:) = ZKA(:)*(XTT-ZZT(:)) +                              &
            ( ZDV(:)*(XLVTT + ( XCPV - XCL ) * ( ZZT(:) - XTT )) &
            *(XESTT-ZZW(:))/(XRV*ZZT(:))         )
!
! compute RHMLTR
!
       ZZW(:)  = MIN( ZRHS(:), MAX( 0.0,( -ZZW(:) *                     &
                              ( X0DEPH*       ZLBDAH(:)**XEX0DEPH +     &
                                X1DEPH*ZCJ(:)*ZLBDAH(:)**XEX1DEPH ) -   &
                       ZZW1(:,6)*( ZRHODREF(:)*XCL*(XTT-ZZT(:))) ) /    &
                                                ( ZRHODREF(:)*XLMTT ) ) )
       ZRRS(:) = ZRRS(:) + ZZW(:)
       ZRHS(:) = ZRHS(:) - ZZW(:)
       ZTHS(:) = ZTHS(:) - ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(-RHMLTR))
!
       ZCRS(:) = MAX( ZCRS(:) + ZZW(:)*(XCCH*ZLBDAH(:)**XCXH/ZRHT(:)),0.0 )
!
    END WHERE

  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HMLT', &
                                           Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'HMLT', &
                                           Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rh ) call Budget_store_end( tbudgets(NBUDGET_RH), 'HMLT', &
                                           Unpack( zrhs(:), mask = gmicro(:, :, :), field = prhs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nr), 'HMLT', &
                                           Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
  end if
END IF

END IF HAIL
!
!------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_MIXED_FAST_PROCESSES
