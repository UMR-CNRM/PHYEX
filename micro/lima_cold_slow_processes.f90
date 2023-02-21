!MNH_LIC Copyright 2013-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #####################
       MODULE MODI_LIMA_COLD_SLOW_PROCESSES
!      #####################
!
INTERFACE
      SUBROUTINE LIMA_COLD_SLOW_PROCESSES (PTSTEP, KMI, PZZ, PRHODJ,                 &
                                           PRHODREF, PEXNREF, PPABST,                &
                                           PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                           PTHS, PRVS, PRIS, PRSS,                   &
                                           PCIT, PCIS, PCST, PCSS                    )
!
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT    ! Ice crystal C. at t
REAL, DIMENSION(:,:,:),OPTIONAL,   INTENT(IN)    :: PCST    ! Ice crystal C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS    ! Ice crystal C. source
REAL, DIMENSION(:,:,:),OPTIONAL,   INTENT(INOUT) :: PCSS    ! Snow/aggregates C. source  
!
END SUBROUTINE LIMA_COLD_SLOW_PROCESSES
END INTERFACE
END MODULE MODI_LIMA_COLD_SLOW_PROCESSES
!
!     ################################################################################
      SUBROUTINE LIMA_COLD_SLOW_PROCESSES (PTSTEP, KMI, PZZ, PRHODJ,                 &
                                           PRHODREF, PEXNREF, PPABST,                &
                                           PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                           PTHS, PRVS, PRIS, PRSS,                   &
                                           PCIT, PCIS, PCST, PCSS                    )
!     ################################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the microphysical sources
!!    for slow cold processes :
!!      - conversion of snow to ice
!!      - deposition of vapor on snow
!!      - conversion of ice to snow (Harrington 1995)
!!      - aggregation of ice on snow
!!
!!
!!    REFERENCE
!!    ---------
!!
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
!  P. Wautelet    03/2020: use the new data structures and subroutines for budgets
!  J. Wurtz       03/2022: new snow characteristics
!  M. Taufour     07/2022: add concentration for snow, graupel, hail
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,          only: lbu_enable, nbumod,                                          &
                                lbudget_th, lbudget_rv, lbudget_ri, lbudget_rs, lbudget_sv,  &
                                NBUDGET_TH, NBUDGET_RV, NBUDGET_RI, NBUDGET_RS, NBUDGET_SV1, &
                                tbudgets
USE MODD_CST,             ONLY: XP00, XRD, XRV, XMV, XMD, XCPD, XCPV,        &
                                XCL, XCI, XTT, XLSTT, XALPI, XBETAI, XGAMI
USE MODD_NSV,             ONLY: NSV_LIMA_NI, NSV_LIMA_NS
USE MODD_PARAMETERS,      ONLY: JPHEXT, JPVEXT
USE MODD_PARAM_LIMA,      ONLY: LSNOW_T, XRTMIN, XCTMIN,              &
                                XALPHAI, XALPHAS, XNUI, XNUS, NMOM_S
USE MODD_PARAM_LIMA_COLD, ONLY: XLBI, XLBEXI, XLBS, XLBEXS, XNS, XBI, XCXS, XCCS, &
                                XLBDAS_MAX, XDSCNVI_LIM, XLBDASCNVI_MAX,     &
                                XC0DEPSI, XC1DEPSI, XR0DEPSI, XR1DEPSI,      &
                                XSCFAC, X1DEPS, X0DEPS, XEX1DEPS, XEX0DEPS,  &
                                XDICNVS_LIM, XLBDAICNVS_LIM,                 &
                                XC0DEPIS, XC1DEPIS, XR0DEPIS, XR1DEPIS,      &
                                XCOLEXIS, XAGGS_CLARGE1, XAGGS_CLARGE2,      &
                                XAGGS_RLARGE1, XAGGS_RLARGE2, XBS,           &
                                XLBDAS_MIN,XFVELOS,XTRANS_MP_GAMMAS

use mode_budget,          only: Budget_store_init, Budget_store_end
use mode_tools,           only: Countjv

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT    ! Ice crystal C. at t
REAL, DIMENSION(:,:,:),OPTIONAL,   INTENT(IN)    :: PCST    ! Snow/aggregates C. at t 
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIS    ! Ice crystal C. source
REAL, DIMENSION(:,:,:),OPTIONAL,   INTENT(INOUT) :: PCSS    ! Snow/aggregates C. source  
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
			  :: GMICRO ! Computations only where necessary
INTEGER :: IMICRO
INTEGER , DIMENSION(SIZE(GMICRO)) :: I1,I2,I3 ! Used to replace PACK
INTEGER                           :: JL       ! and PACK intrinsics
!
REAL, DIMENSION(:), ALLOCATABLE :: ZRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRRT    ! Rain water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(:), ALLOCATABLE :: ZCIT    ! Pristine ice conc. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZCST    ! Snow/aggregate conc. at t !
!
REAL, DIMENSION(:), ALLOCATABLE :: ZRVS    ! Water vapor m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRSS    ! Snow/aggregate m.r. source
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZTHS    ! Theta source
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZCIS    ! Pristine ice conc. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCSS    ! Snow/aggregates conc. source 
!
REAL, DIMENSION(:), ALLOCATABLE &
                   :: ZRHODREF, & ! RHO Dry REFerence
                      ZRHODJ,   & ! RHO times Jacobian
                      ZZT,      & ! Temperature
                      ZPRES,    & ! Pressure
                      ZEXNREF,  & ! EXNer Pressure REFerence
                      ZZW,      & ! Work array
                      ZZX,      & ! Work array
                      ZLSFACT,  & ! L_s/(Pi_ref*C_ph)
                      ZSSI,     & ! Supersaturation over ice
                      ZLBDAI,   & ! Slope parameter of the ice crystal distr.
                      ZLBDAS,   & ! Slope parameter of the aggregate distr.
                      ZAI,      & ! Thermodynamical function
                      ZCJ,      & ! used to compute the ventilation coefficient
                      ZKA,      & ! Thermal conductivity of the air
                      ZDV,      & ! Diffusivity of water vapor in the air
                      ZVISCA      ! Viscosity of air
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZW1      ! Work arrays
!
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZT, ZW ! Temperature
!
INTEGER :: IIB, IIE, IJB, IJE, IKB, IKE        ! Physical domain
!
REAL,    DIMENSION(:),   ALLOCATABLE :: ZRTMIN, ZCTMIN
!
!-------------------------------------------------------------------------------
!
! Physical domain
!
IIB=1+JPHEXT
IIE=SIZE(PZZ,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PZZ,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
!
! Physical limitations
!
ALLOCATE(ZRTMIN(SIZE(XRTMIN)))
ALLOCATE(ZCTMIN(SIZE(XCTMIN)))
ZRTMIN(:) = XRTMIN(:) / PTSTEP
ZCTMIN(:) = XCTMIN(:) / PTSTEP
!
! Temperature
ZT(:,:,:)  = PTHT(:,:,:) * ( PPABST(:,:,:)/XP00 ) ** (XRD/XCPD)
!
! Looking for regions where computations are necessary
!
GMICRO(:,:,:) = .FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                            &
     PRIT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(4) .OR. &
     PRST(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(5)
!
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
!
IF( IMICRO >= 1 ) THEN
!
!------------------------------------------------------------------------------
!
!
!*       1.    Optimization : packing variables
!              --------------------------------
!
!
!
   ALLOCATE(ZRVT(IMICRO)) 
   ALLOCATE(ZRCT(IMICRO)) 
   ALLOCATE(ZRRT(IMICRO)) 
   ALLOCATE(ZRIT(IMICRO)) 
   ALLOCATE(ZRST(IMICRO)) 
   ALLOCATE(ZRGT(IMICRO)) 
!
   ALLOCATE(ZCIT(IMICRO)) 
!
   ALLOCATE(ZRVS(IMICRO))  
   ALLOCATE(ZRIS(IMICRO))
   ALLOCATE(ZRSS(IMICRO))
!
   ALLOCATE(ZTHS(IMICRO))
!
   ALLOCATE(ZCIS(IMICRO))
   ALLOCATE(ZCST(IMICRO))
   ALLOCATE(ZCSS(IMICRO)) 
! 
   ALLOCATE(ZRHODREF(IMICRO)) 
   ALLOCATE(ZZT(IMICRO)) 
   ALLOCATE(ZPRES(IMICRO)) 
   ALLOCATE(ZEXNREF(IMICRO))
   DO JL=1,IMICRO   
      ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
      ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
      ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
      ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
      ZRST(JL) = PRST(I1(JL),I2(JL),I3(JL))
      ZRGT(JL) = PRGT(I1(JL),I2(JL),I3(JL))
!
      ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
!
      ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
      ZRIS(JL) = PRIS(I1(JL),I2(JL),I3(JL))
      ZRSS(JL) = PRSS(I1(JL),I2(JL),I3(JL))
!
      ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
!
      ZCIS(JL) = PCIS(I1(JL),I2(JL),I3(JL))
      if (NMOM_S.GE.2) then
              ZCST(JL) = PCST(I1(JL),I2(JL),I3(JL)) 
              ZCSS(JL) = PCSS(I1(JL),I2(JL),I3(JL)) 
      end if
!
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
      ZZT(JL)      = ZT(I1(JL),I2(JL),I3(JL))
      ZPRES(JL)    = PPABST(I1(JL),I2(JL),I3(JL))
      ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
   ENDDO
!
   IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
      ALLOCATE(ZRHODJ(IMICRO))
      ZRHODJ(:) = PACK( PRHODJ(:,:,:),MASK=GMICRO(:,:,:) )
   END IF
!
!
!------------------------------------------------------------------------------
!
!
!*       2.    Microphysical computations
!              --------------------------
! 
!
   ALLOCATE(ZZW(IMICRO))
   ALLOCATE(ZZX(IMICRO))
   ALLOCATE(ZLSFACT(IMICRO))
   ALLOCATE(ZSSI(IMICRO))
   ALLOCATE(ZLBDAI(IMICRO)) 
   ALLOCATE(ZLBDAS(IMICRO))
   ALLOCATE(ZAI(IMICRO))
   ALLOCATE(ZCJ(IMICRO))
   ALLOCATE(ZKA(IMICRO))
   ALLOCATE(ZDV(IMICRO))
   ALLOCATE(ZZW1(IMICRO,7))
!
! Preliminary computations
!
   ZZW(:)  = ZEXNREF(:)*( XCPD+XCPV*ZRVT(:)+XCL*(ZRCT(:)+ZRRT(:)) &
                                   +XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)) )
!
   ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZW(:) ! L_s/(Pi_ref*C_ph)
!
   ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) )
   ZSSI(:) = ZRVT(:)*( ZPRES(:)-ZZW(:) ) / ( (XMV/XMD) * ZZW(:) ) - 1.0
                                                       ! Supersaturation over ice
! Distribution parameters for ice and snow
   ZLBDAI(:)  = 1.E10
   WHERE (ZRIT(:)>XRTMIN(4) .AND. ZCIT(:)>XCTMIN(4))
      ZLBDAI(:) = ( XLBI*ZCIT(:) / ZRIT(:) )**XLBEXI
   END WHERE
   ZLBDAS(:)  = 1.E10
   IF (LSNOW_T .AND. NMOM_S.EQ.1) THEN
      WHERE(ZZT(:)>263.15 .AND. ZRST(:)>XRTMIN(5)) 
         ZLBDAS(:) = MAX(MIN(XLBDAS_MAX, 10**(14.554-0.0423*ZZT(:))),XLBDAS_MIN)
      END WHERE
      WHERE(ZZT(:)<=263.15 .AND. ZRST(:)>XRTMIN(5)) 
         ZLBDAS(:) = MAX(MIN(XLBDAS_MAX, 10**(6.226-0.0106*ZZT(:))),XLBDAS_MIN)
      END WHERE
      ZLBDAS(:) = ZLBDAS(:) * XTRANS_MP_GAMMAS
      ZCST(:) = (XNS*ZRST(:)*ZLBDAS(:)**XBS)
      ZCSS(:) = (XNS*ZRSS(:)*ZLBDAS(:)**XBS)
   ELSE IF (NMOM_S.GE.2) THEN
	WHERE (ZRST(:)>XRTMIN(5) .AND. ZCST(:)>XCTMIN(5))
                ZLBDAS(:) = ( XLBS*ZCST(:) / ZRST(:) )**XLBEXS 
        END WHERE
   ELSE
      WHERE (ZRST(:)>XRTMIN(5) )
         ZLBDAS(:) = MAX(MIN(XLBDAS_MAX,XLBS*( ZRHODREF(:)*ZRST(:) )**XLBEXS),XLBDAS_MIN)
      END WHERE
      ZCST(:) = XCCS*ZLBDAS(:)**XCXS / ZRHODREF(:)
      ZCSS(:) = XCCS*ZLBDAS(:)**XCXS / ZRHODREF(:)
   END IF
!
   ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
   ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
!
! Thermodynamical function ZAI = A_i(T,P)
   ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) )**2 / (ZKA(:)*XRV*ZZT(:)**2) &
                                         + ( XRV*ZZT(:) ) / (ZDV(:)*ZZW(:))
! ZCJ = c^prime_j/c_i (in the ventilation factor) ( c_i from v(D)=c_i*D^(d_i) )
   ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
!
!
!
!
!*       2.1    Conversion of snow to r_i: RSCNVI
!        ----------------------------------------
!
      if ( nbumod == kmi .and. lbu_enable ) then
        if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'CNVI', pris(:, :, :) * prhodj(:, :, :) )
        if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'CNVI', prss(:, :, :) * prhodj(:, :, :) )
        if ( lbudget_sv ) &
          call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CNVI', pcis(:, :, :) * prhodj(:, :, :) )
        if ( lbudget_sv .AND. NMOM_S.GE.2) &
          call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'CNVI', pcss(:, :, :) * prhodj(:, :, :) )
      end if
!
      ZZW(:) = 0.0
      WHERE ( ZLBDAS(:)<XLBDASCNVI_MAX .AND. (ZRST(:)>XRTMIN(5)) .AND.(ZCST(:)>XCTMIN(5))   &
                                         .AND. (ZSSI(:)<0.0)       )
         ZZW(:) = (ZLBDAS(:)*XDSCNVI_LIM)**(XALPHAS)
         ZZX(:) = ( -ZSSI(:)/ZAI(:) ) * (ZCST(:)) * (ZZW(:)**XNUS) * EXP(-ZZW(:))
!
         ZZW(:) = MIN( ( XR0DEPSI+XR1DEPSI*ZCJ(:) )*ZZX(:),ZRSS(:) )
         ZRIS(:) = ZRIS(:) + ZZW(:)
         ZRSS(:) = ZRSS(:) - ZZW(:)
!
         ZZW(:) = MIN(ZZW(:)*( XC0DEPSI+XC1DEPSI*ZCJ(:) )/( XR0DEPSI+XR1DEPSI*ZCJ(:)),ZCSS(:) )
         ZCIS(:) = ZCIS(:) + ZZW(:)
         ZCSS(:) = ZCSS(:) - ZZW(:)
      END WHERE
!
! Budget storage
      if ( nbumod == kmi .and. lbu_enable ) then
        if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'CNVI', &
                                               Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'CNVI', &
                                               Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CNVI', &
                                               Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_sv .AND. NMOM_S.GE.2 ) &
            call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'CNVI', &
                                               Unpack( zcss(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
      end if
!
!
!*       2.2    Deposition of water vapor on r_s: RVDEPS
!        -----------------------------------------------
!
      if ( nbumod == kmi .and. lbu_enable ) then
        if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'DEPS', pths(:, :, :) * prhodj(:, :, :) )
        if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'DEPS', prvs(:, :, :) * prhodj(:, :, :) )
        if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'DEPS', &
                                               Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
      end if

      ZZW(:) = 0.0
         WHERE ( (ZRST(:)>XRTMIN(5)) .AND. (ZRSS(:)>ZRTMIN(5)) )
            ZZW(:) = ( ZSSI(:)/(ZRHODREF(:)*ZAI(:)) ) * ZCST(:) *  &
                      ( X0DEPS*ZLBDAS(:)**XEX0DEPS +               &
                        X1DEPS*ZCJ(:)*ZLBDAS(:)**XEX1DEPS *        &
                             (1+0.5*(XFVELOS/ZLBDAS(:))**XALPHAS)**(-XNUS+XEX1DEPS/XALPHAS) )
            ZZW(:) =    MIN( ZRVS(:),ZZW(:)      )*(0.5+SIGN(0.5,ZZW(:))) &
                      - MIN( ZRSS(:),ABS(ZZW(:)) )*(0.5-SIGN(0.5,ZZW(:)))
            ZRSS(:) = ZRSS(:) + ZZW(:)
            ZRVS(:) = ZRVS(:) - ZZW(:)
            ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
         END WHERE
! Budget storage
      if ( nbumod == kmi .and. lbu_enable ) then
        if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'DEPS', &
                                               Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'DEPS', &
                                               Unpack( zrvs(:), mask = gmicro(:, :, :), field = prvs(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'DEPS', &
                                               Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
      end if
!
!*       2.3    Conversion of pristine ice to r_s: RICNVS
!        ------------------------------------------------
!
      if ( nbumod == kmi .and. lbu_enable ) then
        if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'CNVS', &
                                               Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'CNVS', &
                                               Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CNVS', &
                                               Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_sv .AND. NMOM_S.GE.2 ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'CNVS', &
                                               Unpack( zcss(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
      end if

      ZZW(:) = 0.0
      WHERE ( (ZLBDAI(:)<XLBDAICNVS_LIM) .AND. (ZCIT(:)>XCTMIN(4)) &
                                            .AND. (ZSSI(:)>0.0)       )
         ZZW(:) = (ZLBDAI(:)*XDICNVS_LIM)**(XALPHAI)
         ZZX(:) = ( ZSSI(:)/ZAI(:) )*ZCIT(:) * (ZZW(:)**XNUI) *EXP(-ZZW(:))
!
! Correction BVIE
!         ZZW(:) = MAX( MIN( ( XR0DEPIS + XR1DEPIS*ZCJ(:) )*ZZX(:)/ZRHODREF(:) &
         ZZW(:) = MAX( MIN( ( XR0DEPIS + XR1DEPIS*ZCJ(:) )*ZZX(:) , ZRIS(:) ) , 0. )
         ZRIS(:) = ZRIS(:) - ZZW(:)
         ZRSS(:) = ZRSS(:) + ZZW(:)
!
         ZZW(:) = MIN( ZZW(:)*(( XC0DEPIS+XC1DEPIS*ZCJ(:) )                   &
                                /( XR0DEPIS+XR1DEPIS*ZCJ(:) )),ZCIS(:) )
         ZCIS(:) = ZCIS(:) - ZZW(:)
         ZCSS(:) = ZCSS(:) + ZZW(:)
         END WHERE
!
! Budget storage
      if ( nbumod == kmi .and. lbu_enable ) then
        if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'CNVS', &
                                               Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'CNVS', &
                                               Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'CNVS', &
                                               Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_sv .AND. NMOM_S.GE.2 ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ns), 'CNVS', &
                                               Unpack( zcss(:), mask = gmicro(:, :, :), field = pcss(:, :, :) ) * prhodj(:, :, :) )
      end if
!
!
!*       2.4    Aggregation of r_i on r_s: CIAGGS and RIAGGS
!        ---------------------------------------------------
!
      if ( nbumod == kmi .and. lbu_enable ) then
        if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'AGGS', &
                                               Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_rs ) call Budget_store_init( tbudgets(NBUDGET_RS), 'AGGS', &
                                               Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'AGGS', &
                                               Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
      end if
!
      WHERE ( (ZRIT(:)>XRTMIN(4)) .AND. (ZRST(:)>XRTMIN(5)) .AND. (ZRIS(:)>ZRTMIN(4)) &
                                                               .AND. (ZCIS(:)>ZCTMIN(4)) )
         ZZW1(:,3) = (ZLBDAI(:) / ZLBDAS(:))**3
         ZZW1(:,1) = (ZCIT(:)*ZCST(:)*EXP( XCOLEXIS*(ZZT(:)-XTT) )) / (ZLBDAI(:)**3)
         ZZW1(:,2) = MIN( ZZW1(:,1)*(XAGGS_CLARGE1+XAGGS_CLARGE2*ZZW1(:,3)),ZCIS(:) )
         ZCIS(:) = ZCIS(:) - ZZW1(:,2)
!
         ZZW1(:,1) = ZZW1(:,1) / ZLBDAI(:)**XBI
         ZZW1(:,2) = MIN( ZZW1(:,1)*(XAGGS_RLARGE1+XAGGS_RLARGE2*ZZW1(:,3)),ZRIS(:) )
         ZRIS(:) = ZRIS(:) - ZZW1(:,2)
         ZRSS(:) = ZRSS(:) + ZZW1(:,2)
      END WHERE
! Budget storage
      if ( nbumod == kmi .and. lbu_enable ) then
        if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'AGGS', &
                                               Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_rs ) call Budget_store_end( tbudgets(NBUDGET_RS), 'AGGS', &
                                               Unpack( zrss(:), mask = gmicro(:, :, :), field = prss(:, :, :) ) * prhodj(:, :, :) )
        if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'AGGS', &
                                               Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
      end if
!------------------------------------------------------------------------------
!
!
!*       3.    Unpacking & Deallocating
!              ------------------------
!
!
  ZW(:,:,:) = PRVS(:,:,:)
  PRVS(:,:,:) = UNPACK( ZRVS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PRIS(:,:,:)
  PRIS(:,:,:) = UNPACK( ZRIS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PRSS(:,:,:)
  PRSS(:,:,:) = UNPACK( ZRSS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
  ZW(:,:,:) = PCIS(:,:,:)
  PCIS(:,:,:) = UNPACK( ZCIS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  if (NMOM_S.GE.2) then
    ZW(:,:,:) = PCSS(:,:,:)
    PCSS(:,:,:) = UNPACK( ZCSS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  end if

!
  ZW(:,:,:) = PTHS(:,:,:)
  PTHS(:,:,:) = UNPACK( ZTHS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
  DEALLOCATE(ZRVT) 
  DEALLOCATE(ZRCT) 
  DEALLOCATE(ZRRT) 
  DEALLOCATE(ZRIT) 
  DEALLOCATE(ZRST) 
  DEALLOCATE(ZRGT) 
  DEALLOCATE(ZCIT) 
  DEALLOCATE(ZRVS)  
  DEALLOCATE(ZRIS)
  DEALLOCATE(ZRSS)
  DEALLOCATE(ZTHS)
  DEALLOCATE(ZCIS)
  DEALLOCATE(ZCSS)  
  DEALLOCATE(ZCST)  
  DEALLOCATE(ZRHODREF) 
  DEALLOCATE(ZZT) 
  DEALLOCATE(ZPRES) 
  DEALLOCATE(ZEXNREF)
  DEALLOCATE(ZZW)
  DEALLOCATE(ZZX)
  DEALLOCATE(ZLSFACT)
  DEALLOCATE(ZSSI)
  DEALLOCATE(ZLBDAI) 
  DEALLOCATE(ZLBDAS)
  DEALLOCATE(ZAI)
  DEALLOCATE(ZCJ)
  DEALLOCATE(ZKA)
  DEALLOCATE(ZDV)
  DEALLOCATE(ZZW1)
  IF (NBUMOD==KMI .AND. LBU_ENABLE) DEALLOCATE(ZRHODJ)
!
END IF
!
!++cb++
DEALLOCATE(ZRTMIN)
DEALLOCATE(ZCTMIN)
!--cb--
!
END SUBROUTINE LIMA_COLD_SLOW_PROCESSES

