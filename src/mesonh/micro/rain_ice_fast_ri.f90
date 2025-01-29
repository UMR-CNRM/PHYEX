!MNH_LIC Copyright 1995-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 25/02/2019: split rain_ice (cleaner and easier to maintain/debug)
!  P. Wautelet 05/06/2019: optimisations
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_FAST_RI

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: RAIN_ICE_FAST_RI

CONTAINS

SUBROUTINE RAIN_ICE_FAST_RI(OMICRO, PRHODREF, PRIT, PRHODJ, PZT, PSSI, PLSFACT, PLVFACT, &
                            PAI, PCJ, PCIT, PRCS, PRIS, PTHS)
!
!*      0. DECLARATIONS
!          ------------
!
use modd_budget,         only: lbudget_th, lbudget_rc, lbudget_ri, &
                               NBUDGET_TH, NBUDGET_RC, NBUDGET_RI, &
                               tbudgets
USE MODD_CST,            only: XTT
USE MODD_RAIN_ICE_DESCR_n, only: XDI, XLBEXI, XLBI, XRTMIN
USE MODD_RAIN_ICE_PARAM_n, only: X0DEPI, X2DEPI

use mode_budget,         only: Budget_store_add, Budget_store_end, Budget_store_init
#ifdef MNH_OPENACC
USE MODE_MNH_ZWORK, ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
#endif
use mode_mppdb

#if defined(MNH_BITREP) || defined(MNH_BITREP_OMP)
USE MODI_BITREP
#endif
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,  DIMENSION(:,:,:), intent(in)    :: OMICRO   ! Test where to compute all processes
REAL,     DIMENSION(:),     intent(in)    :: PRHODREF ! RHO Dry REFerence
REAL,     DIMENSION(:),     intent(in)    :: PRIT     ! Pristine ice m.r. at t
REAL,     DIMENSION(:),     intent(in)    :: PRHODJ   ! RHO times Jacobian
REAL,     DIMENSION(:),     intent(in)    :: PZT      ! Temperature
REAL,     DIMENSION(:),     intent(in)    :: PSSI     ! Supersaturation over ice
REAL,     DIMENSION(:),     intent(in)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL,     DIMENSION(:),     intent(in)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL,     DIMENSION(:),     intent(in)    :: PAI      ! Thermodynamical function
REAL,     DIMENSION(:),     intent(in)    :: PCJ      ! Function to compute the ventilation coefficient
REAL,     DIMENSION(:),     intent(inout) :: PCIT     ! Pristine ice conc. at t
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRCS     ! Cloud water m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PRIS     ! Pristine ice m.r. source
REAL,     DIMENSION(:),     INTENT(INOUT) :: PTHS     ! Theta source
!
!*       0.2  declaration of local variables
!
#ifndef MNH_OPENACC
LOGICAL, DIMENSION(:), ALLOCATABLE :: GWORK
REAL,    DIMENSION(:), ALLOCATABLE :: ZZW  ! Work array
REAL,    DIMENSION(:), ALLOCATABLE :: ZLBEXI
#else
LOGICAL, DIMENSION(:), pointer, contiguous :: GWORK
REAL,    DIMENSION(:), pointer, contiguous :: ZZW  ! Work array
REAL,    DIMENSION(:), pointer, contiguous :: ZLBEXI
#endif
!
INTEGER :: JL,JLU
!-------------------------------------------------------------------------------
!
! IN variables
!
!$acc data present( OMICRO, PRHODREF, PRIT,              &
!$acc &             PRHODJ, PZT, PSSI, PLSFACT, PLVFACT, &
!$acc &             PAI, PCJ,                            &
!
! INOUT variables
!
!$acc &             PCIT, PRCS, PRIS, PTHS )
!
! OUT variables
!
!NONE

IF (MPPDB_INITIALIZED) THEN
  !Check all IN arrays
  CALL MPPDB_CHECK(OMICRO,"RAIN_ICE_FAST_RI beg:OMICRO")
  CALL MPPDB_CHECK(PRHODREF,"RAIN_ICE_FAST_RI beg:PRHODREF")
  CALL MPPDB_CHECK(PRIT,"RAIN_ICE_FAST_RI beg:PRIT")
  CALL MPPDB_CHECK(PRHODJ,"RAIN_ICE_FAST_RI beg:PRHODJ")
  CALL MPPDB_CHECK(PZT,"RAIN_ICE_FAST_RI beg:PZT")
  CALL MPPDB_CHECK(PSSI,"RAIN_ICE_FAST_RI beg:PSSI")
  CALL MPPDB_CHECK(PLSFACT,"RAIN_ICE_FAST_RI beg:PLSFACT")
  CALL MPPDB_CHECK(PLVFACT,"RAIN_ICE_FAST_RI beg:PLVFACT")
  CALL MPPDB_CHECK(PAI,"RAIN_ICE_FAST_RI beg:PAI")
  CALL MPPDB_CHECK(PCJ,"RAIN_ICE_FAST_RI beg:PCJ")
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PCIT,"RAIN_ICE_FAST_RI beg:PCIT")
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_FAST_RI beg:PRCS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_FAST_RI beg:PRIS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_FAST_RI beg:PTHS")
END IF
!
JLU = size(PRHODREF)
!
#ifndef MNH_OPENACC
ALLOCATE( GWORK(size(PRHODREF)) )
ALLOCATE( ZZW  (size(PRHODREF)) )
ALLOCATE( ZLBEXI   (size(PRHODREF)) )
#else
!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_FAST_RI' )

CALL MNH_MEM_GET( GWORK,  SIZE(PRHODREF) )
CALL MNH_MEM_GET( ZZW,    SIZE(PRHODREF) )
CALL MNH_MEM_GET( ZLBEXI, SIZE(PRHODREF) )

!$acc data present( GWORK, ZZW , ZLBEXI )
#endif

!
!*       7.1    cloud ice melting
!
  if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'IMLT', Unpack ( pths(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'IMLT', Unpack ( prcs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'IMLT', Unpack ( pris(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )

!$acc kernels
  GWORK(:) = PRIS(:)>0.0 .AND. PZT(:)>XTT
 !$mnh_expand_where(JL=1:JLU)
  WHERE( GWORK(:) )
    PRCS(:) = PRCS(:) + PRIS(:)
    PTHS(:) = PTHS(:) - PRIS(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(-RIMLTC))
    PRIS(:) = 0.0
    PCIT(:) = 0.0
  END WHERE
 !$mnh_end_expand_where()
!$acc end kernels

  if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'IMLT', Unpack ( pths(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'IMLT', Unpack ( prcs(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'IMLT', Unpack ( pris(:) * prhodj(:), &
                                           mask = omicro(:,:,:), field = 0. ) )
!
!*       7.2    Bergeron-Findeisen effect: RCBERI
!
!$acc kernels
  zzw(:) = 0.
  GWORK(:) = PRCS(:)>0.0 .AND. PSSI(:)>0.0 .AND. PRIT(:)>XRTMIN(4) .AND. PCIT(:)>0.0
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
 !$mnh_expand_where(JL=1:JLU)
  WHERE( GWORK(:) )
    ZZW(:) = MIN(1.E8,XLBI*( PRHODREF(:)*PRIT(:)/PCIT(:) )**XLBEXI) ! Lbda_i
    ZZW(:) = MIN( PRCS(:),( PSSI(:) / (PRHODREF(:)*PAI(:)) ) * PCIT(:) * &
                  ( X0DEPI/ZZW(:) + X2DEPI*PCJ(:)*PCJ(:)/ZZW(:)**(XDI+2.0) ) )
    PRCS(:) = PRCS(:) - ZZW(:)
    PRIS(:) = PRIS(:) + ZZW(:)
    PTHS(:) = PTHS(:) + ZZW(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(RCBERI))
 END WHERE
 !$mnh_end_expand_where()
#else

!!$ Le DO concurrent n'est pas bit-reproductible BUG NVHPC 20.7
 !$mnh_do_concurrent( JL=1:JLU )
    ZLBEXI(JL) = XLBEXI
    IF ( GWORK(JL)  ) THEN
       ZZW(JL) = MIN(1.E8,XLBI*BR_POW( PRHODREF(JL)*PRIT(JL)/PCIT(JL), ZLBEXI(JL) ) ) ! Lbda_i
       ZZW(JL) = MIN( PRCS(JL),( PSSI(JL) / (PRHODREF(JL)*PAI(JL)) ) * PCIT(JL) * &
            ( X0DEPI/ZZW(JL) + X2DEPI*PCJ(JL)*PCJ(JL)/BR_POW(ZZW(JL),XDI+2.0) ) )
       PRCS(JL) = PRCS(JL) - ZZW(JL)
       PRIS(JL) = PRIS(JL) + ZZW(JL)
       PTHS(JL) = PTHS(JL) + ZZW(JL)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RCBERI))
    END IF
 !$mnh_end_do() ! CONCURRENT
 
!!! WHERE( GWORK(:) )
   !!!! ZLBEXI(:) = XLBEXI
    !!!!ZZW(:) = MIN(1.E8,XLBI*BR_POW( PRHODREF(:)*PRIT(:)/PCIT(:), ZLBEXI(:) ) ) ! Lbda_i
    !!!!ZZW(:) = MIN( PRCS(:),( PSSI(:) / (PRHODREF(:)*PAI(:)) ) * PCIT(:) * &
        !!! ( X0DEPI/ZZW(:) + X2DEPI*PCJ(:)*PCJ(:)/BR_POW(ZZW(:),XDI+2.0) ) )
    !!!PRCS(:) = PRCS(:) - ZZW(:)
    !!!!PRIS(:) = PRIS(:) + ZZW(:)
    !!!!PTHS(:) = PTHS(:) + ZZW(:)*(PLSFACT(:)-PLVFACT(:)) ! f(L_f*(RCBERI))
 !!!!END WHERE
 
#endif 
!$acc end kernels

  if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'BERFI', Unpack (  zzw(:) * ( plsfact(:) - plvfact(:) ) &
                                                                     * prhodj(:), mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_rc ) call Budget_store_add( tbudgets(NBUDGET_RC), 'BERFI', Unpack ( -zzw(:) * prhodj(:), &
                                                                                  mask = omicro(:,:,:), field = 0. ) )
  if ( lbudget_ri ) call Budget_store_add( tbudgets(NBUDGET_RI), 'BERFI', Unpack (  zzw(:) * prhodj(:), &
                                                                                  mask = omicro(:,:,:), field = 0. ) )
IF (MPPDB_INITIALIZED) THEN
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PCIT,"RAIN_ICE_FAST_RI end:PCIT")
  CALL MPPDB_CHECK(PRCS,"RAIN_ICE_FAST_RI end:PRCS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_FAST_RI end:PRIS")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_FAST_RI end:PTHS")
END IF

!$acc end data

#ifdef MNH_OPENACC
!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'RAIN_ICE_FAST_RI' )
#endif

!$acc end data

END SUBROUTINE RAIN_ICE_FAST_RI

END MODULE MODE_RAIN_ICE_FAST_RI
