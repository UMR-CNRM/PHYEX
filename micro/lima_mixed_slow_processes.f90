!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #####################################
       MODULE MODI_LIMA_MIXED_SLOW_PROCESSES
!      #####################################
!
INTERFACE
      SUBROUTINE LIMA_MIXED_SLOW_PROCESSES(ZRHODREF, ZZT, ZSSI, PTSTEP,        &
                                           ZLSFACT, ZLVFACT, ZAI, ZCJ,         &
                                           ZRGT, ZRHT, ZCIT, ZCGT, ZCHT,       &
                                           ZRVS, ZRCS, ZRIS, ZRGS, ZRHS, ZTHS, &
                                           ZCCS, ZCIS, ZCGS, ZIFS, ZINS,       &
                                           ZLBDAI, ZLBDAG, ZLBDAH,             &
                                           PRHODJ1D, GMICRO, PRHODJ, KMI,      &
                                           PTHS, PRVS, PRCS, PRIS, PRGS, PRHS, &
                                           PCCS, PCIS                    )
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(:),   INTENT(IN)    :: ZZT       ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: ZSSI      ! Supersaturation over ice
REAL,                 INTENT(IN)    :: PTSTEP    ! Time-step
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZAI       ! Thermodynamical function
REAL, DIMENSION(:),   INTENT(IN)    :: ZCJ       ! for the ventilation coefficient
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRGT      ! Graupel/hail m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRHT      ! Hail m.r. at t         
REAL, DIMENSION(:),   INTENT(IN)    :: ZCIT      ! Pristine ice conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZCGT      ! Graupel conc. at t  
REAL, DIMENSION(:),   INTENT(IN)    :: ZCHT      ! Hail conc. at t   
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRVS      ! Water vapor m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRCS      ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRIS      ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRGS      ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRHS      ! Graupel/hail m.r. source 
REAL, DIMENSION(:),   INTENT(INOUT) :: ZTHS      ! Theta source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCCS      ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCIS      ! Pristine ice conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCGS      ! Graupel conc. source 
REAL, DIMENSION(:,:), INTENT(INOUT) :: ZIFS      ! Free Ice nuclei conc. source
REAL, DIMENSION(:,:), INTENT(INOUT) :: ZINS      ! Nucleated Ice nuclei conc. source 
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAI  ! Slope parameter of the ice crystal distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAG  ! Slope parameter of the graupel distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAH  ! Slope parameter of the hail distr. 
!
! used for budget storage
REAL,    DIMENSION(:),     INTENT(IN) :: PRHODJ1D
LOGICAL, DIMENSION(:,:,:), INTENT(IN) :: GMICRO 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHODJ
INTEGER,                   INTENT(IN) :: KMI 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PTHS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRVS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRIS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRGS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCIS
!
END SUBROUTINE LIMA_MIXED_SLOW_PROCESSES
END INTERFACE
END MODULE MODI_LIMA_MIXED_SLOW_PROCESSES
!
!     #######################################################################
      SUBROUTINE LIMA_MIXED_SLOW_PROCESSES(ZRHODREF, ZZT, ZSSI, PTSTEP,        &
                                           ZLSFACT, ZLVFACT, ZAI, ZCJ,         &
                                           ZRGT, ZRHT, ZCIT, ZCGT, ZCHT,       &
                                           ZRVS, ZRCS, ZRIS, ZRGS, ZRHS, ZTHS, &
                                           ZCCS, ZCIS, ZCGS, ZIFS, ZINS,       &
                                           ZLBDAI, ZLBDAG, ZLBDAH,             &
                                           PRHODJ1D, GMICRO, PRHODJ, KMI,      &
                                           PTHS, PRVS, PRCS, PRIS, PRGS, PRHS, &
                                           PCCS, PCIS                    )
!     #######################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the mixed-phase 
!!    slow processes : 
!!
!!      Deposition of water vapor on graupeln
!!      Cloud ice Melting
!!      Bergeron-Findeisen effect
!!
!!**  METHOD
!!    ------
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
!!      C. Barthe  * LACy *  jan. 2014   add budgets
!  P. Wautelet    03/2020: use the new data structures and subroutines for budgets
!  P. Wautelet 02/02/2021: budgets: add missing source terms for SV budgets in LIMA
!  M. Taufour     07/2022: add concentration for snow, graupel, hail
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,           only: lbu_enable, nbumod,                                                                  &
                                 lbudget_th, lbudget_rv, lbudget_rc, lbudget_rc, lbudget_ri, lbudget_rg, lbudget_rh, lbudget_sv,  &
                                 NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RC, NBUDGET_RI, NBUDGET_RG, NBUDGET_RH, NBUDGET_SV1, &
                                 tbudgets
USE MODD_CST,              ONLY : XTT, XALPI, XBETAI, XGAMI,          &
                                       XALPW, XBETAW, XGAMW
USE MODD_NSV
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCTMIN, NMOD_IFN, NMOM_S, NMOM_G, NMOM_H
USE MODD_PARAM_LIMA_COLD,  ONLY : XDI, X0DEPI, X2DEPI, XSCFAC
USE MODD_PARAM_LIMA_MIXED, ONLY : XLBG, XLBEXG, XLBDAG_MAX, XCCG, XCXG, &
                                  X0DEPG, XEX0DEPG, X1DEPG, XEX1DEPG,   &
                                  X0DEPH, XEX0DEPH, X1DEPH, XEX1DEPH  
use mode_budget,           only: Budget_store_add, Budget_store_init, Budget_store_end

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(:),   INTENT(IN)    :: ZZT       ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: ZSSI      ! Supersaturation over ice
REAL,                 INTENT(IN)    :: PTSTEP    ! Time-step
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZAI       ! Thermodynamical function
REAL, DIMENSION(:),   INTENT(IN)    :: ZCJ       ! for the ventilation coefficient
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRGT      ! Graupel/hail m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZRHT      ! Hail m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: ZCIT      ! Pristine ice conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZCGT      ! Graupel conc. at t     
REAL, DIMENSION(:),   INTENT(IN)    :: ZCHT      ! hail conc. at t  
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRVS      ! Water vapor m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRCS      ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRIS      ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRGS      ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRHS      ! hail m.r. source 
REAL, DIMENSION(:),   INTENT(INOUT) :: ZTHS      ! Theta source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCCS      ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCIS      ! Pristine ice conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCGS      ! Graupel conc. source 
REAL, DIMENSION(:,:), INTENT(INOUT) :: ZIFS      ! Free Ice nuclei conc. source
REAL, DIMENSION(:,:), INTENT(INOUT) :: ZINS      ! Nucleated Ice nuclei conc. source 
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAI  ! Slope parameter of the ice crystal distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAG  ! Slope parameter of the graupel distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAH  ! Slope parameter of the hail distr. 
!
! used for budget storage
REAL,    DIMENSION(:),     INTENT(IN) :: PRHODJ1D
LOGICAL, DIMENSION(:,:,:), INTENT(IN) :: GMICRO 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHODJ
INTEGER,                   INTENT(IN) :: KMI
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PTHS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRVS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRIS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRGS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHS 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCIS
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(ZZT)) :: ZZW, ZMASK    ! Work vectors
!
INTEGER :: JMOD_IFN
!
!-------------------------------------------------------------------------------
!
!*       1    Deposition of water vapor on r_g: RVDEPG
!        ---------------------------------------------
!
!
IF (NMOM_S.GE.1) THEN
   ZZW(:) = 0.0
   WHERE ( (ZRGT(:)>XRTMIN(6)) .AND. (ZRGS(:)>XRTMIN(6)/PTSTEP) ) 
      ZZW(:) = ( ZSSI(:)/ZAI(:)/ZRHODREF(:) ) *  ZCGT(:) *               &
           ( X0DEPG*ZLBDAG(:)**XEX0DEPG + X1DEPG*ZCJ(:)*ZLBDAG(:)**XEX1DEPG )
      ZZW(:) =         MIN( ZRVS(:),ZZW(:)      )*(0.5+SIGN(0.5,ZZW(:))) &
                     - MIN( ZRGS(:),ABS(ZZW(:)) )*(0.5-SIGN(0.5,ZZW(:)))
      ZRGS(:) = ZRGS(:) + ZZW(:)
      ZRVS(:) = ZRVS(:) - ZZW(:)
      ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
   END WHERE
!
! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'DEPG', &
                                           Unpack(  zzw(:) * zlsfact(:) * prhodj1d(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'DEPG', &
                                           Unpack( -zzw(:)              * prhodj1d(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rg ) call Budget_store_add( tbudgets(NBUDGET_RG), 'DEPG', &
                                           Unpack(  zzw(:)              * prhodj1d(:), mask = gmicro(:, :, :), field = 0. ) )
  end if
END IF
!
!                                                                          
!*       1.0  Deposition of water vapor on r_h: RVDEPH
!        ---------------------------------------------
!
!
IF (NMOM_H.GE.2) THEN
   ZZW(:) = 0.0
   WHERE ( (ZRHT(:)>XRTMIN(7)) .AND. (ZRHS(:)>XRTMIN(7)/PTSTEP) )
      ZZW(:) = ( ZSSI(:)/(ZAI(:)) ) *  ZCHT(:) *                                 &   
               ( X0DEPH*ZLBDAH(:)**XEX0DEPH + X1DEPH*ZCJ(:)*ZLBDAH(:)**XEX1DEPH )
      ZZW(:) =         MIN( ZRVS(:),ZZW(:)      )*(0.5+SIGN(0.5,ZZW(:))) &
                     - MIN( ZRHS(:),ABS(ZZW(:)) )*(0.5-SIGN(0.5,ZZW(:)))
      ZRHS(:) = ZRHS(:) + ZZW(:)
      ZRVS(:) = ZRVS(:) - ZZW(:)
      ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
   END WHERE
!
! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_add( tbudgets(NBUDGET_TH), 'DEPH', &
                                           Unpack(  zzw(:) * zlsfact(:) * prhodj1d(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rv ) call Budget_store_add( tbudgets(NBUDGET_RV), 'DEPH', &
                                           Unpack( -zzw(:)              * prhodj1d(:), mask = gmicro(:, :, :), field = 0. ) )
    if ( lbudget_rh ) call Budget_store_add( tbudgets(NBUDGET_RH), 'DEPH', &
                                           Unpack(  zzw(:)              * prhodj1d(:), mask = gmicro(:, :, :), field = 0. ) )
  end if
END IF
!       
!
!*       2    cloud ice Melting: RIMLTC and CIMLTC
!        -----------------------------------------
!
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'IMLT', &
                                           Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'IMLT', prcs(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'IMLT', pris(:, :, :) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'IMLT', pccs(:, :, :) * prhodj(:, :, :) )
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'IMLT', pcis(:, :, :) * prhodj(:, :, :) )
!      do jmod_ifn = 1,nmod_ifn
!        call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl + jmod_ifn - 1), 'IMLT', &
!                                pins(:, :, :, jmod_ifn) * prhodj(:, :, :) )
!      enddo
    end if
  end if

   ZMASK(:) = 1.0
   WHERE( (ZRIS(:)>XRTMIN(4)/PTSTEP) .AND. (ZZT(:)>XTT) )
      ZRCS(:) = ZRCS(:) + ZRIS(:)
      ZTHS(:) = ZTHS(:) - ZRIS(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(-RIMLTC))
      ZRIS(:) = 0.0
!
      ZCCS(:) = ZCCS(:) + ZCIS(:)
      ZCIS(:) = 0.0
      ZMASK(:)= 0.0
   END WHERE
   DO JMOD_IFN = 1,NMOD_IFN
! Correction BVIE aerosols not released but in droplets
!      ZIFS(:,JMOD_IFN) = ZIFS(:,JMOD_IFN) + ZINS(:,JMOD_IFN)*(1.-ZMASK(:)) 
!      ZINS(:,JMOD_IFN) = ZINS(:,JMOD_IFN) * ZMASK(:)
   ENDDO
!
! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'IMLT', &
                                           Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'IMLT', &
                                           Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'IMLT', &
                                           Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_sv ) then
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_nc), 'IMLT', &
                                           Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'IMLT', &
                                           Unpack( zcis(:), mask = gmicro(:, :, :), field = pcis(:, :, :) ) * prhodj(:, :, :) )
!      do jmod_ifn = 1,nmod_ifn
!        call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ifn_nucl + jmod_ifn - 1), 'IMLT', &
!                          Unpack( zins(:, jmod_ifn), mask = gmicro(:, :, :), field = pins(:, :, :, jmod_ifn) ) * prhodj(:, :, :) )
!      enddo
    end if
  end if
!
!*       3    Bergeron-Findeisen effect: RCBERI
!        --------------------------------------
!
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'BERFI', &
                                           Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'BERFI', &
                                           Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_init( tbudgets(NBUDGET_RI), 'BERFI', &
                                           Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  end if

   ZZW(:) = 0.0
   WHERE( (ZRCS(:)>XRTMIN(2)/PTSTEP) .AND. (ZRIS(:)>XRTMIN(4)/PTSTEP) .AND. (ZCIT(:)>XCTMIN(4)) )
      ZZW(:) = EXP( (XALPW-XALPI) - (XBETAW-XBETAI)/ZZT(:)          &
                                  - (XGAMW-XGAMI)*ALOG(ZZT(:)) ) -1.0 
                                      ! supersaturation of saturated water over ice
      ZZW(:) = MIN( ZRCS(:),( ZZW(:) / ZAI(:) ) * ZCIT(:) *        &
                    ( X0DEPI/ZLBDAI(:)+X2DEPI*ZCJ(:)*ZCJ(:)/ZLBDAI(:)**(XDI+2.0) ) )
      ZRCS(:) = ZRCS(:) - ZZW(:)
      ZRIS(:) = ZRIS(:) + ZZW(:)
      ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RCBERI))
   END WHERE
!
! Budget storage
  if ( nbumod == kmi .and. lbu_enable ) then
    if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'BERFI', &
                                           Unpack( zths(:), mask = gmicro(:, :, :), field = pths(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'BERFI', &
                                           Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
    if ( lbudget_ri ) call Budget_store_end( tbudgets(NBUDGET_RI), 'BERFI', &
                                           Unpack( zris(:), mask = gmicro(:, :, :), field = pris(:, :, :) ) * prhodj(:, :, :) )
  end if
!------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_MIXED_SLOW_PROCESSES
