!MNH_LIC Copyright 1995-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 25/02/2019: split rain_ice (cleaner and easier to maintain/debug)
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet 29/05/2019: remove PACK/UNPACK intrinsics (to get more performance and better OpenACC support)
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!-----------------------------------------------------------------
MODULE MODE_RAIN_ICE_NUCLEATION

  IMPLICIT NONE

  PRIVATE

  PUBLIC RAIN_ICE_NUCLEATION

CONTAINS

SUBROUTINE RAIN_ICE_NUCLEATION(KIB, KIE, KJB, KJE, KKTB, KKTE,KRR,PTSTEP,&
     PTHT,PPABST,PRHODJ,PRHODREF,PRVT,PRCT,PRRT,PRIT,PRST,PRGT,&
     PCIT,PEXNREF,PTHS,PRVS,PRIS,PT,PRHT)
!
!*      0. DECLARATIONS
!          ------------
!
use modd_budget,          only: lbudget_th, lbudget_rv, lbudget_ri, &
                                NBUDGET_TH, NBUDGET_RV, NBUDGET_RI, &
                                tbudgets
USE MODD_CST,             only: XALPI, XALPW, XBETAI, XBETAW, XCI, XCL, XCPD, XCPV, XGAMI, XGAMW, &
                                XLSTT, XMD, XMV, XP00, XRD, XTT
USE MODD_RAIN_ICE_PARAM_n,  only: XALPHA1, XALPHA2, XBETA1, XBETA2, XMNU0, XNU10, XNU20

use mode_budget,          only: Budget_store_init, Budget_store_end
use mode_mppdb
#ifndef MNH_OPENACC
use mode_tools,           only: Countjv
#else
use mode_tools,           only: Countjv_device
#endif

#if defined(MNH_BITREP) || defined(MNH_BITREP_OMP)
use modi_bitrep
#endif

#ifdef MNH_OPENACC
USE MODE_MNH_ZWORK,       ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
#endif

#if defined(MNH_COMPILER_CCE) && defined(MNH_BITREP_OMP)
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KIB, KIE, KJB, KJE, KKTB, KKTE
INTEGER,                          INTENT(IN)    :: KRR     ! Number of moist variable
REAL,                             INTENT(IN)    :: PTSTEP  ! Double Time step
                                                           ! (single if cold start)
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PTHT    ! Theta at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(:,:,:),           INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:),           INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),           INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),           INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),           INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),           INTENT(OUT)   :: PT      ! Temperature
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN)    :: PRHT    ! Hail m.r. at t
!
!*       0.2  declaration of local variables
!
INTEGER                                :: INEGT
INTEGER                                :: JL       ! and PACK intrinsics
INTEGER, DIMENSION(:),     POINTER, CONTIGUOUS :: I1,I2,I3 ! Used to replace the COUNT
LOGICAL, DIMENSION(:,:,:), POINTER, CONTIGUOUS :: GNEGT    ! Test where to compute the HEN process
REAL                                   :: ZZWMAX
REAL,    DIMENSION(:), POINTER, CONTIGUOUS     :: ZRVT     ! Water vapor m.r. at t
REAL,    DIMENSION(:), POINTER, CONTIGUOUS     :: ZCIT     ! Pristine ice conc. at t
REAL,    DIMENSION(:), POINTER, CONTIGUOUS     :: ZZT,   & ! Temperature
                                          ZPRES, & ! Pressure
                                          ZZW,   & ! Work array
                                          ZUSW,  & ! Undersaturation over water
                                          ZSSI     ! Supersaturation over ice
REAL,    DIMENSION(:,:,:), POINTER, CONTIGUOUS :: ZW      ! work array

INTEGER :: JIU,JJU,JKU, JIJKU

INTEGER :: JI,JJ,JK
!
!-------------------------------------------------------------------------------
!
! IN variables
!
!$acc data present( PTHT, PPABST, PRHODJ, PRHODREF, PRVT, PRCT, &
!$acc &             PRRT, PRIT, PRST, PRGT, PEXNREF, PRHT,      &
!
! INOUT variables
!
!$acc &             PCIT, PTHS, PRVS, PRIS,                     &
!
! OUT variables
!
!$acc &             PT )

IF (MPPDB_INITIALIZED) THEN
  !Check all IN arrays
  CALL MPPDB_CHECK(PTHT,"RAIN_ICE_NUCLEATION beg:PTHT")
  CALL MPPDB_CHECK(PPABST,"RAIN_ICE_NUCLEATION beg:PPABST")
  CALL MPPDB_CHECK(PRHODJ,"RAIN_ICE_NUCLEATION beg:PRHODJ")
  CALL MPPDB_CHECK(PRHODREF,"RAIN_ICE_NUCLEATION beg:PRHODREF")
  CALL MPPDB_CHECK(PRVT,"RAIN_ICE_NUCLEATION beg:PRVT")
  CALL MPPDB_CHECK(PRCT,"RAIN_ICE_NUCLEATION beg:PRCT")
  CALL MPPDB_CHECK(PRRT,"RAIN_ICE_NUCLEATION beg:PRRT")
  CALL MPPDB_CHECK(PRIT,"RAIN_ICE_NUCLEATION beg:PRIT")
  CALL MPPDB_CHECK(PRST,"RAIN_ICE_NUCLEATION beg:PRST")
  CALL MPPDB_CHECK(PRGT,"RAIN_ICE_NUCLEATION beg:PRGT")
  CALL MPPDB_CHECK(PEXNREF,"RAIN_ICE_NUCLEATION beg:PEXNREF")
  IF(PRESENT(PRHT)) CALL MPPDB_CHECK(PRHT,"RAIN_ICE_NUCLEATION beg:PRHT")
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PCIT,"RAIN_ICE_NUCLEATION beg:PCIT")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_NUCLEATION beg:PTHS")
  CALL MPPDB_CHECK(PRVS,"RAIN_ICE_NUCLEATION beg:PRVS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_NUCLEATION beg:PRIS")
END IF

JIU =  size(PEXNREF, 1 )
JJU =  size(PEXNREF, 2 )
JKU =  size(PEXNREF, 3 )
JIJKU = JIU * JJU * JKU

#ifndef MNH_OPENACC
allocate( i1( JIJKU ) )
allocate( i2( JIJKU ) )
allocate( i3( JIJKU ) )
allocate( gnegt( JIU,JJU,JKU ) )
allocate( zw   ( JIU,JJU,JKU ) )
#else
!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_NUCLEATION 1' )

CALL MNH_MEM_GET( i1, JIJKU )
CALL MNH_MEM_GET( i2, JIJKU )
CALL MNH_MEM_GET( i3, JIJKU )
CALL MNH_MEM_GET( gnegt, JIU,JJU,JKU )
CALL MNH_MEM_GET( zw   , JIU,JJU,JKU )
#endif

!$acc data present( i1, i2, i3, gnegt, zw )

if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HIN', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'HIN', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'HIN', pris(:, :, :) * prhodj(:, :, :) )
!
!  compute the temperature and the pressure
!
!$acc kernels
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
PT(:,:,:) = PTHT(:,:,:) * ( PPABST(:,:,:) / XP00 ) ** ( XRD / XCPD )
#else
!$mnh_do_concurrent ( JI=1:JIU,JJ=1:JJU,JK=1:JKU)
   PT(JI,JJ,JK) = PTHT(JI,JJ,JK) * BR_POW( PPABST(JI,JJ,JK) / XP00, XRD / XCPD )
!$mnh_end_do()
#endif
!
!$acc end kernels

!  optimization by looking for locations where
!  the temperature is negative only !!!
!$acc kernels present_cr(GNEGT)
GNEGT(:,:,:) = .FALSE.
GNEGT(KIB:KIE,KJB:KJE,KKTB:KKTE) = PT(KIB:KIE,KJB:KJE,KKTB:KKTE)<XTT
!$acc end kernels
#ifndef MNH_OPENACC
INEGT = COUNTJV( GNEGT(:,:,:),I1(:),I2(:),I3(:))
#else
CALL COUNTJV_DEVICE(GNEGT(:,:,:),I1(:),I2(:),I3(:),INEGT)
#endif
IF( INEGT >= 1 ) THEN
#ifndef MNH_OPENACC   
  ALLOCATE(ZRVT(INEGT)) 
  ALLOCATE(ZCIT(INEGT)) 
  ALLOCATE(ZZT(INEGT))  
  ALLOCATE(ZPRES(INEGT))
  ALLOCATE(ZZW(INEGT))
  ALLOCATE(ZUSW(INEGT))
  ALLOCATE(ZSSI(INEGT))
#else
  !Pin positions in the pools of MNH memory
  CALL MNH_MEM_POSITION_PIN( 'RAIN_ICE_NUCLEATION 2' )

  CALL MNH_MEM_GET( ZRVT,  INEGT )
  CALL MNH_MEM_GET( ZCIT,  INEGT )
  CALL MNH_MEM_GET( ZZT,   INEGT )
  CALL MNH_MEM_GET( ZPRES, INEGT )
  CALL MNH_MEM_GET( ZZW,   INEGT )
  CALL MNH_MEM_GET( ZUSW,  INEGT )
  CALL MNH_MEM_GET( ZSSI,  INEGT )
#endif
  
!$acc data present( zrvt, zcit, zzt, zpres, zzw, zusw, zssi )

  !$acc kernels 
  !$mnh_do_concurrent ( JL=1:INEGT )
    ZRVT(JL)  = PRVT  (I1(JL),I2(JL),I3(JL))
    ZCIT(JL)  = PCIT  (I1(JL),I2(JL),I3(JL))
    ZZT(JL)   = PT    (I1(JL),I2(JL),I3(JL))
    ZPRES(JL) = PPABST(I1(JL),I2(JL),I3(JL))
  !$mnh_end_do()
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
  !$mnh_do_concurrent ( JL=1:INEGT )
    ZZW(JL) = EXP( XALPI - XBETAI/ZZT(JL) - XGAMI*LOG(ZZT(JL) ) )      ! es_i
  !$mnh_end_do()
#else
  !$mnh_do_concurrent ( JL=1:INEGT )
    ZZW(JL) = BR_EXP( XALPI - XBETAI/ZZT(JL) - XGAMI*BR_LOG(ZZT(JL) ) )      ! es_i
  !$mnh_end_do()
#endif
  !$mnh_do_concurrent ( JL=1:INEGT )
    ZZW(JL) = MIN(ZPRES(JL)/2., ZZW(JL))             ! safety limitation
    ZSSI(JL) = ZRVT(JL)*( ZPRES(JL)-ZZW(JL) ) / ( (XMV/XMD) * ZZW(JL) ) - 1.0
                         ! Supersaturation over ice
  !$mnh_end_do()
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
  !$mnh_do_concurrent ( JL=1:INEGT )
    ZUSW(JL) = EXP( XALPW - XBETAW/ZZT(JL) - XGAMW*LOG(ZZT(JL) ) )     ! es_w
    ZUSW(JL) = MIN(ZPRES(JL)/2.,ZUSW(JL))            ! safety limitation
    ZUSW(JL) = ( ZUSW(JL)/ZZW(JL) )*( (ZPRES(JL)-ZZW(JL))/(ZPRES(JL)-ZUSW(JL)) ) - 1.0
                       ! Supersaturation of saturated water vapor over ice
  !$mnh_end_do()
#else
  !$mnh_do_concurrent ( JL=1:INEGT )
    ZUSW(JL) = BR_EXP( XALPW - XBETAW/ZZT(JL) - XGAMW*BR_LOG(ZZT(JL) ) )     ! es_w
    ZUSW(JL) = MIN(ZPRES(JL)/2.,ZUSW(JL))            ! safety limitation
    ZUSW(JL) = ( ZUSW(JL)/ZZW(JL) )*( (ZPRES(JL)-ZZW(JL))/(ZPRES(JL)-ZUSW(JL)) ) - 1.0
                       ! Supersaturation of saturated water vapor over ice
  !$mnh_end_do()
#endif  
!
!*       3.1     compute the heterogeneous nucleation source: RVHENI
!
!*       3.1.1   compute the cloud ice concentration
!
!$mnh_do_concurrent ( JL=1:INEGT )
  ZZW(JL) = 0.0
  ZSSI(JL) = MIN( ZSSI(JL), ZUSW(JL) ) ! limitation of SSi according to SSw=0
  IF ( (ZZT(JL)<XTT-5.0) .AND. (ZSSI(JL)>0.0) ) THEN
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
    ZZW(JL) = XNU20 * EXP( XALPHA2*ZSSI(JL)-XBETA2 )
#else
    ZZW(JL) = XNU20 * BR_EXP( XALPHA2*ZSSI(JL)-XBETA2 )
#endif
 END IF
!$mnh_end_do()
#if !defined(MNH_BITREP) && !defined(MNH_BITREP_OMP)
 !$mnh_do_concurrent ( JL=1:INEGT )
    IF ( (ZZT(JL)<=XTT-2.0) .AND. (ZZT(JL)>=XTT-5.0) .AND. (ZSSI(JL)>0.0) ) THEN
       ZZW(JL) = MAX( XNU20 * EXP( -XBETA2 ),XNU10 * EXP( -XBETA1*(ZZT(JL)-XTT) ) * &
            ( ZSSI(JL)/ZUSW(JL) )**XALPHA1 )
    END IF
 !$mnh_end_do() ! CONCURRENT
#else
 !$mnh_do_concurrent ( JL=1:INEGT )
    IF ( (ZZT(JL)<=XTT-2.0) .AND. (ZZT(JL)>=XTT-5.0) .AND. (ZSSI(JL)>0.0) ) THEN
       ZZW(JL) = MAX( XNU20 * BR_EXP( -XBETA2 ),XNU10 * BR_EXP( -XBETA1*(ZZT(JL)-XTT) ) * &
            BR_POW( ZSSI(JL)/ZUSW(JL),XALPHA1 ) )
    END IF
 !$mnh_end_do() ! CONCURRENT
#endif
 ! WARNING COMPILER BUG NVHPC20.X/3 <-> if array syntaxe ZZW(1:INEGT) = ZZW(1:INEGT)
 !$mnh_do_concurrent ( JL=1:INEGT )
    ZZW(JL) = ZZW(JL) - ZCIT(JL)
 !$mnh_end_do()
#ifndef MNH_COMPILER_NVHPC
 ZZWMAX = MAXVAL(ZZW(1:INEGT))
!$acc end kernels
#else
!$acc end kernels
 ZZWMAX = 0.0
!$acc parallel reduction(max:ZZWMAX)
 !$mnh_do_concurrent( JL=1:INEGT)
  ZZWMAX = MAX(ZZWMAX,ZZW(JL))
 !$mnh_end_do()
!$acc end parallel 
#endif 


  IF( ZZWMAX > 0.0 ) THEN
  !$acc kernels     
!
!*       3.1.2   update the r_i and r_v mixing ratios
     !
  !$mnh_do_concurrent ( JL=1:INEGT )
     ZZW(JL) = MIN( ZZW(JL),50.E3 ) ! limitation provisoire a 50 l^-1
  !$mnh_end_do()
    ZW(:,:,:) = 0.0
    !$mnh_do_concurrent ( JL=1:INEGT )
      ZW(I1(JL), I2(JL), I3(JL)) = ZZW( JL )
    !$mnh_end_do()
    ZW(:,:,:) = MAX( ZW(:,:,:) ,0.0 ) *XMNU0/(PRHODREF(:,:,:)*PTSTEP)
    PRIS(:,:,:) = PRIS(:,:,:) + ZW(:,:,:)
    PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
    IF ( KRR == 7 ) THEN
      PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(XLSTT+(XCPV-XCI)*(PT(:,:,:)-XTT))   &
                   /( (XCPD + XCPV*PRVT(:,:,:) + XCL*(PRCT(:,:,:)+PRRT(:,:,:))   &
       + XCI*(PRIT(:,:,:)+PRST(:,:,:)+PRGT(:,:,:)+PRHT(:,:,:)))*PEXNREF(:,:,:) )
    ELSE IF( KRR == 6 ) THEN
      PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(XLSTT+(XCPV-XCI)*(PT(:,:,:)-XTT))   &
                   /( (XCPD + XCPV*PRVT(:,:,:) + XCL*(PRCT(:,:,:)+PRRT(:,:,:))   &
                   + XCI*(PRIT(:,:,:)+PRST(:,:,:)+PRGT(:,:,:)))*PEXNREF(:,:,:) )
    END IF
                                 ! f(L_s*(RVHENI))
    !$mnh_do_concurrent ( JL=1:INEGT )
       ZZW(JL) = MAX( ZZW(JL)+ZCIT(JL),ZCIT(JL) )
    !$mnh_end_do()
    PCIT(:,:,:) = MAX( PCIT(:,:,:), 0.0 )
    !$mnh_do_concurrent ( JL=1:INEGT )
      PCIT(I1(JL), I2(JL), I3(JL)) = MAX( ZZW( JL ), PCIT(I1(JL), I2(JL), I3(JL)), 0.0 )
   !$mnh_end_do()
  !$acc end kernels
END IF

!$acc end data
  
#ifndef MNH_OPENACC
  DEALLOCATE(ZSSI)
  DEALLOCATE(ZUSW)
  DEALLOCATE(ZZW)
  DEALLOCATE(ZPRES)
  DEALLOCATE(ZZT)
  DEALLOCATE(ZCIT)
  DEALLOCATE(ZRVT)
#else
  !Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
  CALL MNH_MEM_RELEASE( 'RAIN_ICE_NUCLEATION 2' )
#endif
END IF
!
!*       3.1.3   budget storage
!
if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HIN', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'HIN', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'HIN', pris(:, :, :) * prhodj(:, :, :) )
!
IF (MPPDB_INITIALIZED) THEN
  !Check all INOUT arrays
  CALL MPPDB_CHECK(PCIT,"RAIN_ICE_NUCLEATION end:PCIT")
  CALL MPPDB_CHECK(PTHS,"RAIN_ICE_NUCLEATION end:PTHS")
  CALL MPPDB_CHECK(PRVS,"RAIN_ICE_NUCLEATION end:PRVS")
  CALL MPPDB_CHECK(PRIS,"RAIN_ICE_NUCLEATION end:PRIS")
  !Check all OUT arrays
  CALL MPPDB_CHECK(PT,"RAIN_ICE_NUCLEATION end:PT")
END IF

#ifndef MNH_OPENACC
deallocate (i1, i2, i3, gnegt, zw )
#else
!Release all memory allocated with MNH_MEM_GET calls since beginning of subroutine
CALL MNH_MEM_RELEASE( 'RAIN_ICE_NUCLEATION 1' )
#endif
!$acc end data

!$acc end data

END SUBROUTINE RAIN_ICE_NUCLEATION

END MODULE MODE_RAIN_ICE_NUCLEATION
