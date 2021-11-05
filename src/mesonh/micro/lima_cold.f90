!MNH_LIC Copyright 2013-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #####################
       MODULE MODI_LIMA_COLD
!      #####################
!
INTERFACE
      SUBROUTINE LIMA_COLD (OSEDI, OHHONI, KSPLITG, PTSTEP, KMI,           &
                           KRR, PZZ, PRHODJ,                               &
                           PRHODREF, PEXNREF, PPABST, PW_NU,               &
                           PTHM, PPABSM,                                   &
                           PTHT, PRT, PSVT,                                &
                           PTHS, PRS, PSVS,                                &
                           PINPRS, PINPRG, PINPRH)
!
USE MODD_NSV,   only: NSV_LIMA_BEG
!
LOGICAL,                  INTENT(IN)    :: OSEDI   ! switch to activate the
                                                   ! cloud ice sedimentation
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
INTEGER,                  INTENT(IN)    :: KSPLITG ! Number of small time step 
                                                   ! for ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
!
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM    ! Theta at time t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABSM  ! abs. pressure at time t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT     ! m.r. at t 
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(IN)    :: PSVT ! Concentrations at time t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS     ! m.r. source
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(INOUT) :: PSVS ! Concentration sources
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS  ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRH  ! Hail instant precip
!
END SUBROUTINE LIMA_COLD
END INTERFACE
END MODULE MODI_LIMA_COLD
!
!     ######################################################################
      SUBROUTINE LIMA_COLD (OSEDI, OHHONI, KSPLITG, PTSTEP, KMI,           &
                           KRR, PZZ, PRHODJ,                               &
                           PRHODREF, PEXNREF, PPABST, PW_NU,               &
                           PTHM, PPABSM,                                   &
                           PTHT, PRT, PSVT,                                &
                           PTHS, PRS, PSVS,                                &
                           PINPRS, PINPRG, PINPRH)
!     ######################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the cold-phase 
!!    microphysical sources involving only primary ice and snow, except for 
!!    the sedimentation which also includes graupelns, and the homogeneous
!!    freezing of CCNs, cloud droplets and raindrops.
!!
!!
!!**  METHOD
!!    ------
!!      The nucleation of IFN is parameterized following either Meyers (1992)
!!    or Phillips (2008, 2013).
!!
!!      The sedimentation rates are computed with a time spliting technique: 
!!    an upstream scheme, written as a difference of non-advective fluxes. 
!!    This source term is added to the next coming time step (split-implicit 
!!    process).
!!
!!
!!    REFERENCES
!!    ----------
!!
!!      Cohard, J.-M. and J.-P. Pinty, 2000: A comprehensive two-moment warm 
!!      microphysical bulk scheme. 
!!        Part I: Description and tests
!!        Part II: 2D experiments with a non-hydrostatic model
!!      Accepted for publication in Quart. J. Roy. Meteor. Soc. 
!!
!!      Phillips et al., 2008: An empirical parameterization of heterogeneous
!!        ice nucleation for multiple chemical species of aerosols, J. Atmos. Sci. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets (no more budget calls in this subroutine)
!  P. Wautelet 28/05/2020: bugfix: correct array start for PSVT and PSVS
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------

use modd_budget,     only: lbu_enable,                                                  &
                           lbudget_ri, lbudget_rs, lbudget_rg, lbudget_rh, lbudget_sv,  &
                           NBUDGET_RI, NBUDGET_RS, NBUDGET_RG, NBUDGET_RH, NBUDGET_SV1, &
                           tbudgets
USE MODD_NSV
USE MODD_PARAM_LIMA

use mode_budget,          only: Budget_store_init, Budget_store_end

USE MODI_LIMA_COLD_HOM_NUCL
USE MODI_LIMA_COLD_SEDIMENTATION
USE MODI_LIMA_COLD_SLOW_PROCESSES
USE MODI_LIMA_MEYERS
USE MODI_LIMA_PHILLIPS

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL,                  INTENT(IN)    :: OSEDI   ! switch to activate the 
                                                   ! cloud ice sedimentation
LOGICAL,                  INTENT(IN)    :: OHHONI  ! enable haze freezing
INTEGER,                  INTENT(IN)    :: KSPLITG ! Number of small time step 
                                                   ! for ice sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
!
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variables
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM    ! Theta at time t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABSM  ! abs. pressure at time t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT     ! m.r. at t 
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(IN)    :: PSVT ! Concentrations at time t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS     ! m.r. source
REAL, DIMENSION(:,:,:,NSV_LIMA_BEG:), INTENT(INOUT) :: PSVS ! Concentration sources
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRS  ! Snow instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRG  ! Graupel instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRH  ! Hail instant precip
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))  &
                                    :: PRVT,    & ! Water vapor m.r. at t 
                                       PRCT,    & ! Cloud water m.r. at t 
                                       PRRT,    & ! Rain water m.r. at t 
                                       PRIT,    & ! Cloud ice m.r. at t 
                                       PRST,    & ! Snow/aggregate m.r. at t 
                                       PRGT,    & ! Graupel m.r. at t 
                                       PRHT,    & ! Graupel m.r. at t 
                                       !
                                       PRVS,    & ! Water vapor m.r. source
                                       PRCS,    & ! Cloud water m.r. source
                                       PRRS,    & ! Rain water m.r. source
                                       PRIS,    & ! Pristine ice m.r. source
                                       PRSS,    & ! Snow/aggregate m.r. source
                                       PRGS,    & ! Graupel/hail m.r. source
                                       PRHS,    & ! Graupel/hail m.r. source
                                       !
                                       PCCT,    & ! Cloud water C. at t
                                       PCRT,    & ! Rain water C. at t
                                       PCIT,    & ! Ice crystal C. at t
                                       !
                                       PCCS,    & ! Cloud water C. source
                                       PCRS,    & ! Rain water C. source
                                       PCIS       ! Ice crystal C. source
!
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PNFS     ! CCN C. available source
                                                  !used as Free ice nuclei for
                                                  !HOMOGENEOUS nucleation of haze
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PNAS     ! Cloud  C. nuclei C. source
                                                  !used as Free ice nuclei for
                                                  !IMMERSION freezing
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PIFS     ! Free ice nuclei C. source 
                                                  !for DEPOSITION and CONTACT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PINS     ! Activated ice nuclei C. source
                                                  !for DEPOSITION and CONTACT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PNIS     ! Activated ice nuclei C. source
                                                  !for IMMERSION
REAL, DIMENSION(:,:,:), ALLOCATABLE :: PNHS     ! Haze homogeneous activation
!
!-------------------------------------------------------------------------------
!
!
!*       0.     3D MICROPHYSCAL VARIABLES
!	        -------------------------
!
!
! Prepare 3D water mixing ratios
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
PRHT(:,:,:) = 0.
PRHS(:,:,:) = 0.
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
IF ( KRR .GE. 7 ) PRHT(:,:,:) = PRT(:,:,:,7)
IF ( KRR .GE. 7 ) PRHS(:,:,:) = PRS(:,:,:,7)
!
! Prepare 3D number concentrations
PCCT(:,:,:) = 0.
PCRT(:,:,:) = 0.
PCIT(:,:,:) = 0.
PCCS(:,:,:) = 0.
PCRS(:,:,:) = 0.
PCIS(:,:,:) = 0.
!
IF ( LWARM ) PCCT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NC)
IF ( LWARM .AND. LRAIN ) PCRT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NR)
IF ( LCOLD ) PCIT(:,:,:) = PSVT(:,:,:,NSV_LIMA_NI)
!
IF ( LWARM ) PCCS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NC)
IF ( LWARM .AND. LRAIN ) PCRS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NR)
IF ( LCOLD ) PCIS(:,:,:) = PSVS(:,:,:,NSV_LIMA_NI)
!
IF ( NMOD_CCN .GE. 1 ) THEN
   ALLOCATE( PNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   ALLOCATE( PNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_CCN) )
   PNFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1)
   PNAS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1)
ELSE
   ALLOCATE( PNFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   ALLOCATE( PNAS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   PNFS(:,:,:,:) = 0.
   PNAS(:,:,:,:) = 0.
END IF
!
IF ( NMOD_IFN .GE. 1 ) THEN
   ALLOCATE( PIFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
   ALLOCATE( PINS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IFN) )
   PIFS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1)
   PINS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1)
ELSE
   ALLOCATE( PIFS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   ALLOCATE( PINS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   PIFS(:,:,:,:) = 0.
   PINS(:,:,:,:) = 0.
END IF
!
IF ( NMOD_IMM .GE. 1 ) THEN
   ALLOCATE( PNIS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),NMOD_IMM) )
   PNIS(:,:,:,:) = PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1)
ELSE
   ALLOCATE( PNIS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3),1) )
   PNIS(:,:,:,:) = 0.0
END IF
!
IF ( OHHONI ) THEN
   ALLOCATE( PNHS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) )
   PNHS(:,:,:) = PSVS(:,:,:,NSV_LIMA_HOM_HAZE)
ELSE
   ALLOCATE( PNHS(SIZE(PRHODJ,1),SIZE(PRHODJ,2),SIZE(PRHODJ,3)) )
   PNHS(:,:,:) = 0.0
END IF
!
!-------------------------------------------------------------------------------
!
!
!*       1.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!	        -------------------------------------
!
if ( lbu_enable ) then
  if ( lbudget_ri .and. osedi ) call Budget_store_init( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs .and. lsnow ) call Budget_store_init( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg .and. lsnow ) call Budget_store_init( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh .and. lhail ) call Budget_store_init( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv .and. osedi ) &
      call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'SEDI', pcis(:, :, :) * prhodj(:, :, :) )
end if

CALL LIMA_COLD_SEDIMENTATION (OSEDI, KSPLITG, PTSTEP, KMI,     &
                              PZZ, PRHODJ, PRHODREF,           &
                              PRIT, PCIT,                      &
                              PRIS, PRSS, PRGS, PRHS, PCIS,    &
                              PINPRS, PINPRG,&
                              PINPRH                  )

if ( lbu_enable ) then
  if ( lbudget_ri .and. osedi ) call Budget_store_end( tbudgets(NBUDGET_RI), 'SEDI', pris(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rs .and. lsnow ) call Budget_store_end( tbudgets(NBUDGET_RS), 'SEDI', prss(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rg .and. lsnow ) call Budget_store_end( tbudgets(NBUDGET_RG), 'SEDI', prgs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rh .and. lhail ) call Budget_store_end( tbudgets(NBUDGET_RH), 'SEDI', prhs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv .and. osedi ) &
      call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_lima_ni), 'SEDI', pcis(:, :, :) * prhodj(:, :, :) )
end if
!-------------------------------------------------------------------------------
!
!
!               COMPUTE THE NUCLEATION PROCESS SOURCES
!   	        --------------------------------------
!
IF (LNUCL) THEN
!
   IF ( LMEYERS ) THEN
      PIFS(:,:,:,:) = 0.0
      PNIS(:,:,:,:) = 0.0
      CALL LIMA_MEYERS (OHHONI, PTSTEP, KMI,                            &
                        PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST,         &
                        PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PCCT, &
                        PTHS, PRVS, PRCS, PRIS,                         &
                        PCCS, PCIS, PINS )
   ELSE
      CALL LIMA_PHILLIPS (OHHONI, PTSTEP, KMI,                      &
                          PZZ, PRHODJ, PRHODREF, PEXNREF, PPABST,   &
                          PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                          PTHS, PRVS, PRCS, PRIS,                   &
                          PCIT, PCCS, PCIS,                         &
                          PNAS, PIFS, PINS, PNIS   )
   END IF
!
   IF (LWARM .OR. (LHHONI .AND. NMOD_CCN.GE.1)) THEN
      CALL LIMA_COLD_HOM_NUCL (OHHONI, PTSTEP, KMI,                         &
                               PZZ, PRHODJ,                                 &
                               PRHODREF, PEXNREF, PPABST, PW_NU,            &
                               PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT,    &
                               PTHS, PRVS, PRCS, PRRS, PRIS, PRGS,          &
                               PCCT,                                        &
                               PCCS, PCRS, PNFS, PCIS, PNHS                 )
   END IF
!
END IF
!
!------------------------------------------------------------------------------
!
!
!*       4.    SLOW PROCESSES: depositions, aggregation
!              ----------------------------------------
!
IF (LSNOW) THEN
!
   CALL LIMA_COLD_SLOW_PROCESSES(PTSTEP, KMI, PZZ, PRHODJ,                 &
                                 PRHODREF, PEXNREF, PPABST,                &
                                 PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                 PTHS, PRVS, PRIS, PRSS,                   &
                                 PCIT, PCIS                                )
!
END IF
!
!------------------------------------------------------------------------------
!
!
!*       4.    REPORT 3D MICROPHYSICAL VARIABLES IN PRS AND PSVS
!              -------------------------------------------------
!
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
IF ( LWARM ) PSVS(:,:,:,NSV_LIMA_NC) = PCCS(:,:,:)
IF ( LWARM .AND. LRAIN ) PSVS(:,:,:,NSV_LIMA_NR) = PCRS(:,:,:)
IF ( LCOLD ) PSVS(:,:,:,NSV_LIMA_NI) = PCIS(:,:,:)
!
IF ( NMOD_CCN .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_CCN_FREE:NSV_LIMA_CCN_FREE+NMOD_CCN-1) = PNFS(:,:,:,:)
   PSVS(:,:,:,NSV_LIMA_CCN_ACTI:NSV_LIMA_CCN_ACTI+NMOD_CCN-1) = PNAS(:,:,:,:)
END IF
!
IF ( NMOD_IFN .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_IFN_FREE:NSV_LIMA_IFN_FREE+NMOD_IFN-1) = PIFS(:,:,:,:)
   PSVS(:,:,:,NSV_LIMA_IFN_NUCL:NSV_LIMA_IFN_NUCL+NMOD_IFN-1) = PINS(:,:,:,:)
END IF
!
IF ( NMOD_IMM .GE. 1 ) THEN
   PSVS(:,:,:,NSV_LIMA_IMM_NUCL:NSV_LIMA_IMM_NUCL+NMOD_IMM-1) = PNIS(:,:,:,:)
END IF

IF ( OHHONI ) PSVS(:,:,:,NSV_LIMA_HOM_HAZE) = PNHS(:,:,:)
!
!++cb++
IF (ALLOCATED(PNFS)) DEALLOCATE(PNFS)
IF (ALLOCATED(PNAS)) DEALLOCATE(PNAS)
IF (ALLOCATED(PIFS)) DEALLOCATE(PIFS)
IF (ALLOCATED(PINS)) DEALLOCATE(PINS)
IF (ALLOCATED(PNIS)) DEALLOCATE(PNIS)
IF (ALLOCATED(PNHS)) DEALLOCATE(PNHS)
!--cb--
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_COLD
