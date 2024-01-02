!MNH_LIC Copyright 1996-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ######################
       MODULE MODI_RAIN_C2R2_KHKO
!      ######################
!
IMPLICIT NONE
INTERFACE
      SUBROUTINE RAIN_C2R2_KHKO(HCLOUD,OACTIT, OSEDC, ORAIN, KSPLITR, PTSTEP,   &
                            KMI,TPFILE,                                         &
                            PZZ, PRHODJ,                                        &
                            PRHODREF, PEXNREF,                                  &
                            PPABST, PTHT, PRVT, PRCT,                           &
                            PRRT, PTHM, PRCM, PPABSM,                           &
                            PW_NU,PDTHRAD, PTHS, PRVS, PRCS, PRRS,              &
                            PCNT, PCCT, PCRT, PCNS, PCCS, PCRS,                 &
                            PINPRC, PINPRR, PINPRR3D, PEVAP3D,PAEROT,           &
                            PSOLORG, PMI, HACTCCN,                              &
                            PINDEP, PSUPSAT, PNACT                      )
!
USE MODD_IO, ONLY: TFILEDATA
IMPLICIT NONE
!
CHARACTER(LEN=*),         INTENT(IN)    :: HCLOUD   !  kind of cloud

LOGICAL,                  INTENT(IN)    :: OACTIT ! Switch to activate the
                                                  ! activation by radiative
                                                  ! tendency
LOGICAL,                  INTENT(IN)    :: OSEDC   ! switch to activate the 
                                                   ! cloud droplet sedimentation
LOGICAL,                  INTENT(IN)    :: ORAIN   ! switch to activate the 
                                                   ! rain formation by coalescence
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step 
                                      ! integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP ! Time step :XTSTEP in namelist
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE   ! Output file
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM    ! Theta at time t-Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABSM   ! Pressure time t-Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCM    ! Cloud water m.r. at time t-Dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDTHRAD ! THeta RADiative Tendancy
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCNT    ! Water vapor C. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRT    ! Rain water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCNS    ! Water vapor C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS    ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS    ! Rain water C. source
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRR  ! Rain instant precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEVAP3D! Rain evap profile
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PAEROT  ! Aerosol concentration
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PMI
CHARACTER(LEN=4),         INTENT(IN)    :: HACTCCN  ! kind of CCN activation scheme
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSUPSAT  !sursat
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNACT   !concentrtaion d'aérosols activés au temps t
!
END SUBROUTINE RAIN_C2R2_KHKO
END INTERFACE
END MODULE MODI_RAIN_C2R2_KHKO
!     ######################################################################
      SUBROUTINE RAIN_C2R2_KHKO (HCLOUD,OACTIT, OSEDC, ORAIN, KSPLITR, PTSTEP,  &
                            KMI, TPFILE, PZZ, PRHODJ,                           &
                            PRHODREF, PEXNREF,                                  &
                            PPABST, PTHT, PRVT,  PRCT,                          &
                            PRRT, PTHM, PRCM, PPABSM,                           &
                            PW_NU,PDTHRAD, PTHS, PRVS, PRCS, PRRS,              &
                            PCNT, PCCT,  PCRT, PCNS, PCCS, PCRS,                &
                            PINPRC, PINPRR, PINPRR3D, PEVAP3D,PAEROT,           &
                            PSOLORG, PMI, HACTCCN,                              &
                            PINDEP, PSUPSAT, PNACT                      )
!     ######################################################################
!
!!****  * -  compute the explicit microphysical sources of cloud water and
!!           rain water concentrations and mixing ratios
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the microphysical sources
!!    for the two schemes C2R2 and KHKO
!!    For C2R2 the microphysical sources are :                    
!!    nucleation, sedimentation, autoconversion, accretion, self-collection 
!!    and vaporisation which are parameterized according to Cohard and Pinty 
!!    QJRMS, 2000
!!    For KHKO the microphysical sources are : 
!!    drizzle drops sedimentation, autoconversion, accretion and vaporisation   
!!    which are parameterized according to Khairoutdinov and Kogan 2000, 
!!    nucleation and cloud droplets sedimentation which are parameterized  
!!    according to Cohard and Pinty QJRMS, 2000
!!
!!**  METHOD
!!    ------
!!      The activation of CCN is checked for quasi-saturated air parcels 
!!    to update the cloud droplet number concentration. Then assuming a 
!!    generalized gamma distribution law for the cloud droplets and the 
!!    raindrops, the zeroth and third order moments tendencies are evaluated
!!    for all the coalescence terms by integrating the Stochastic Collection 
!!    Equation. As autoconversion is a process that cannot be resolved 
!!    analytically, the Berry-Reinhardt parameterisation is employed with
!!    modifications to initiate the raindrop spectrum mode. The integration
!!    of the raindrop evaporation of the raindrops below clouds is 
!!    straightformward.
!!
!!      The sedimentation rates are computed with a time spliting technique: 
!!    an upstream scheme, written as a difference of non-advective fluxes. 
!!    This source term is added to the next coming time step (split-implicit 
!!    process).
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!          JPHEXT       : Horizontal external points number
!!          JPVEXT       : Vertical external points number
!!      Module MODD_CONF :
!!          CCONF configuration of the model for the first time step
!!
!!      Module MODD_CST     
!!          XP00               ! Reference pressure
!!          XRD,XRV            ! Gaz  constant for dry air, vapor
!!          XMD,XMV            ! Molecular weight for dry air, vapor
!!          XCPD               ! Cpd (dry air)
!!          XCL                ! Cl (liquid)
!!          XTT                ! Triple point temperature
!!          XLVTT              ! Vaporization heat constant
!!          XALPW,XBETAW,XGAMW ! Constants for saturation vapor pressure
!!                             ! function over liquid water
!!      Module MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask 
!!                          'NONE'  ' for no budget
!!         LBU_RTH      : logical for budget of RTH (potential temperature)
!!                        .TRUE. = budget of RTH        
!!                        .FALSE. = no budget of RTH
!!         LBU_RRV      : logical for budget of RRV (water vapor)
!!                        .TRUE. = budget of RRV 
!!                        .FALSE. = no budget of RRV 
!!         LBU_RRC      : logical for budget of RRC (cloud water)
!!                        .TRUE. = budget of RRC 
!!                        .FALSE. = no budget of RRC 
!!         LBU_RRR      : logical for budget of RRR (rain water)
!!                        .TRUE. = budget of RRR 
!!                        .FALSE. = no budget of RRR 
!!
!!    REFERENCE
!!    ---------
!!
!!      Cohard, J.-M. and J.-P. Pinty, 2000: A comprehensive two-moment warm 
!!      microphysical bulk scheme. 
!!        Part I: Description and tests
!!        Part II: 2D experiments with a non-hydrostatic model
!!      Accepted for publication in Quart. J. Roy. Meteor. Soc. 
!!      M. Khairoutdinov and Y. Kogan,"A new Cloud Physics Parametererization
!!      in a Large-Eddy Simulation Model of Marine Stratocumulus"
!!      Mon. Weather Rev.,128, 229-243-2000
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      O. Geoffroy  * CNRM Meteo-France* : 07/2006
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             31/12/96 
!!      Jean-Pierre PINTY     7/ 7/00  Code cleaning
!!      Jean-Pierre PINTY    27/ 5/01  Review of rain transfer to cloud droplets
!!                                     in the case of strong evaporation
!!      C.Lac                11/09     Distinction of the TSTEPs
!!      C.Lac, V.Masson      09/10     Corrections in sedimentation and
!!                                     evaporation for reproducibility
!!      C.Lac                06/14     C2R2_SEDIMENTATION replaced by
!!                                     KHKO_SEDIMENTATION because of instability
!!      G.Tanguy             07/14     FUSION C2R2 and KHKO
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      J.Escobar : 07/10/2015 , Bug in parallel run , => comment test on INUCT>1 containing GET_HALO  
!!      M.Mazoyer : 04/2016 : Temperature radiative tendency used for  
!!                            activation by cooling (OACTIT : mis en commentaires)
!!      M.Mazoyer : 04/2016 : Add supersaturation diagnostics
!!      C.Lac     : 07/2016 : Add droplet deposition
!!      C.Lac     : 01/2017 : Correction on droplet deposition
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet    02/2020: use the new data structures and subroutines for budgets
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
use modd_budget,               only: lbudget_th, lbudget_rv, lbudget_rc, lbudget_rr, lbudget_sv,  &
                                     NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RR, NBUDGET_SV1, &
                                     tbudgets
USE MODD_CH_AEROSOL
USE MODD_CONF
USE MODD_CST
USE MODD_DUST
use modd_field,                only: tfieldmetadata, TYPEREAL
USE MODD_IO,                   ONLY: TFILEDATA
USE MODD_NSV,                  ONLY : NSV_C2R2BEG
USE MODD_PARAM_C2R2
USE MODD_PARAMETERS
USE MODD_RAIN_C2R2_DESCR
USE MODD_RAIN_C2R2_KHKO_PARAM
USE MODD_SALT

use mode_budget,               only: Budget_store_init, Budget_store_end
USE MODE_IO_FIELD_WRITE,       only: IO_Field_write
USE MODE_ll
use mode_tools,                only: Countjv

USE MODI_GAMMA

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
!
CHARACTER(LEN=*),         INTENT(IN)    :: HCLOUD   !  kind of cloud

LOGICAL,                  INTENT(IN)    :: OACTIT ! Switch to activate the
                                                  ! activation by radiative
                                                  ! tendency
LOGICAL,                  INTENT(IN)    :: OSEDC   ! switch to activate the 
                                                   ! cloud droplet sedimentation
LOGICAL,                  INTENT(IN)    :: ORAIN   ! switch to activate the 
                                                   ! rain formation by coalescence
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step 
                                      ! integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP ! Time step :XTSTEP in namelist
INTEGER,                  INTENT(IN)    :: KMI     ! Model index 
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE   ! Output file
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM    ! Theta at time t-Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABSM  ! Pressure time t-Dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCM    ! Cloud water m.r. at time t-Dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PW_NU   ! updraft velocity used for
                                                   ! the nucleation param.
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDTHRAD ! THeta RADiative Tendancy
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCNT    ! Water vapor C. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRT    ! Rain water C. at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCNS    ! Water vapor C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS    ! Cloud water C. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS    ! Rain water C. source
!
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRC  ! Cloud instant precip
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(:,:),     INTENT(INOUT) :: PINPRR  ! Rain instant precip
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEVAP3D! Rain evap profile
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PAEROT  ! Aerosol concentration
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PMI
CHARACTER(LEN=4),         INTENT(IN)    :: HACTCCN  ! kind of CCN activation scheme
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSUPSAT  !sursat
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNACT   !concentrtaion d'aérosols activés au temps t
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK            ! Vertical loop index for the rain sedimentation 
INTEGER :: JKBIS         ! For spectrum
INTEGER :: JN            ! Temporal loop index for the rain sedimentation
INTEGER :: IIB           !  Define the domain where is 
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           ! 
INTEGER :: IJE           !
INTEGER :: IKB           ! 
INTEGER :: IKE           !
INTEGER :: ISIZE         !
!
REAL    :: ZTSPLITR      ! Small time step for rain sedimentation
REAL    :: ZEPS          ! molar mass ratio
!
LOGICAL :: GBU           ! General condition prior calling any BUDGET routine
!
!
INTEGER :: ISEDIM, INUCT, & ! Case number of sedimentation, nucleation,
           IMICRO, IEVAP, & ! coalescence and rain_evaporation locations
           ISELF, IACCR, ISCBU
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
            :: GSEDIM ! Test where to compute the SED processes
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)):: GDEP
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
            :: GNUCT  ! Test where to compute the HEN process
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
            :: GMICRO ! Test where to compute coalescence proc.
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
            :: GEVAP  ! Test where to compute rain_evap. proc.
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1             ! Vectors of indices for
                                                        ! interpolations
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1             ! Work vectors for 
                                                        ! interpolations
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZW ! work array
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZWSEDR, ZWSEDC, &! sedimentation fluxes
                                     ZZW1LOG,         &  ! cloud sedimentation speed
                                     ZWSEDLOGR, ZWSEDLOGC ! sedimentation fluxes
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZT, ZTM,ZTDT, ZDRC ! Temperature
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZZA,ZCHEN
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZRVSAT
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZLV !latent heat of vaporization
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZWLBDR,ZWLBDR3,ZWLBDC,ZWLBDC3, & 
                                     ZWLBDA,  &  !libre parcours moyen
                                     ZRAY,    &  ! Mean volumic radius
                                     ZCC         ! Terminal vertical velocity
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZMVRR,ZVRR,ZVCR
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZPRCT, ZPCCT, ZPRRT, ZPCRT 
                                           ! For split sedimentation
REAL,    DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))   &
                                  :: ZMVRC !Cloud water mean volumic radius
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
     :: ZPRRS,ZPCRS   ! Rain and cloud source for sedim
REAL, DIMENSION(:), ALLOCATABLE :: ZRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:), ALLOCATABLE :: ZRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:), ALLOCATABLE :: ZRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:), ALLOCATABLE :: ZCNT    ! nucleus conc. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZCCT    ! cloud conc. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZCRT    ! rain conc. at t
!
REAL, DIMENSION(:), ALLOCATABLE :: ZRVS    ! Water vapor m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRCS    ! Cloud water m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZRRS    ! Rain water m.r. source
REAL, DIMENSION(:), ALLOCATABLE :: ZCNS    ! nucleus conc. source
REAL, DIMENSION(:), ALLOCATABLE :: ZCCS    ! cloud conc. source
REAL, DIMENSION(:), ALLOCATABLE :: ZTCC    ! Corrective factor for Terminal velocity
REAL, DIMENSION(:), ALLOCATABLE :: ZCRS    ! rain conc. source
REAL, DIMENSION(:), ALLOCATABLE :: ZTHS    ! Theta source
!-------------------------------------------------------------------------------
! Modification of XCHEN according to theta vertical gradient (J. Rangonio)
!REAL, DIMENSION(:), ALLOCATABLE :: ZCHEN_TMP
!REAL, DIMENSION(:), ALLOCATABLE :: ZCONC_CCN
!-------------------------------------------------------------------------------
!
REAL,    DIMENSION(:), ALLOCATABLE :: ZZVRR    !terminal velocity for drop concentration
REAL,    DIMENSION(:), ALLOCATABLE :: ZZVCR    !erminal velocity for rain water
!
LOGICAL, DIMENSION(:), ALLOCATABLE :: GSELF(:), GACCR(:), GSCBU(:)
LOGICAL, DIMENSION(:), ALLOCATABLE :: GENABLE_ACCR_SCBU(:)
REAL, DIMENSION(:), ALLOCATABLE :: &
                  ZRHODREF, & ! RHO Dry REFerence
                  ZZT,      & ! Temperature
                  ZTDTBIS, &  ! dT/dt
                  ZEXNREF,  & ! EXNer Pressure REFerence
                  ZZW1, ZZW2, ZZW3, ZZW4, ZZW5, & ! Work array
                  ZZLV,     & ! Latent heat of vaporization at T
                  ZSMAX,    & ! Maximum supersaturation
                  ZSCBU,    & ! optimisation mask
                  ZLBDC, ZLBDR,           & ! Lambda parameter
                  ZLBDC3, ZLBDR3,         & ! Lambda**3
                  ZKA,      & ! Thermal conductivity of the air
                  ZDV,      & ! Diffusivity of water vapor in the air
                  ZPABST, ZNCN, ZMCN
REAL, DIMENSION(:), ALLOCATABLE    :: ZDG3
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZAERO, ZAEROS, ZSOLORG, ZMI
REAL  :: ZFACT, JSV, ZMU, ZALPHA

REAL, DIMENSION(:), ALLOCATABLE    :: ZRTMIN
REAL, DIMENSION(:), ALLOCATABLE    :: ZCTMIN
REAL :: ZTMP
TYPE(TFIELDMETADATA) :: TZFIELD
!
!
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE SLOPE PARAMETERS ZLBDC,ZLBDR
!   	        ----------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
!
ISIZE = SIZE(XRTMIN)
ISIZE = SIZE(XCTMIN)
ALLOCATE(ZCTMIN(ISIZE))
ALLOCATE(ZRTMIN(ISIZE))
ZRTMIN(:) = XRTMIN(:) / PTSTEP
ZCTMIN(:) = XCTMIN(:) / PTSTEP
!
ZWLBDC3(:,:,:) = 1.E30
ZWLBDC(:,:,:)  = 1.E10
!
WHERE (PRCT(:,:,:)>XRTMIN(2) .AND. PCCT(:,:,:)>XCTMIN(2))
  ZWLBDC3(:,:,:) = XLBC * PCCT(:,:,:) / (PRHODREF(:,:,:) * PRCT(:,:,:))
  ZWLBDC(:,:,:)  = ZWLBDC3(:,:,:)**XLBEXC
END WHERE
!
IF (HCLOUD=='C2R2'.OR. HCLOUD=='C3R5' ) THEN
  ZWLBDR3(:,:,:) = 1.E30
  ZWLBDR(:,:,:)  = 1.E10
  WHERE (PRRT(:,:,:)>XRTMIN(3) .AND. PCRT(:,:,:)>XCTMIN(3))
    ZWLBDR3(:,:,:) = XLBR * PCRT(:,:,:) / (PRHODREF(:,:,:) * PRRT(:,:,:))
    ZWLBDR(:,:,:)  = ZWLBDR3(:,:,:)**XLBEXR
  END WHERE
ENDIF
!
ZT(:,:,:)  = PTHT(:,:,:) * (PPABST(:,:,:)/XP00)**(XRD/XCPD)

!
!*       2.     COMPUTES THE NUCLEATION PROCESS SOURCES
!   	        --------------------------------------
!
IF ((HACTCCN == 'ABRK').AND.((LORILAM).OR.(LDUST).OR.(LSALT))) THEN
  CALL AER_NUCLEATION
ELSE
  IF (.NOT. LSUPSAT) THEN
    CALL C2R2_KHKO_NUCLEATION
  ELSE
    ZEPS= XMV / XMD
    ZT(:,:,:)  = PTHT(:,:,:) * (PPABST(:,:,:)/XP00)**(XRD/XCPD)
!
    ZRVSAT(:,:,:) = ZEPS / (PPABST(:,:,:) * &
                    EXP(-XALPW+XBETAW/ZT(:,:,:)+XGAMW*LOG(ZT(:,:,:))) - 1.0)
  ENDIF
ENDIF
!
!------------------------------------------------------------------------------
!
!*       3.    COALESCENCE PROCESSES
!              ---------------------
!
IF (ORAIN) THEN
!
!  optimization by looking for locations where
!  the microphysical fields are larger than a minimal value only !!!
!
  IF (HCLOUD=='C2R2'.OR. HCLOUD=='C3R5') THEN
    CALL C2R2_COALESCENCE
  ELSE ! KHKO
    CALL KHKO_COALESCENCE
  ENDIF
!
!-------------------------------------------------------------------------------
!
!        4.    EVAPORATION OF RAINDROPS
!              ------------------------
!
  CALL C2R2_KHKO_EVAPORATION
!
!-------------------------------------------------------------------------------
!
!        5.    SPONTANEOUS BREAK-UP (NUMERICAL FILTER)
!              --------------------
!
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2), 'BRKU', pcrs(:, :, :) * prhodj(:, :, :) )

  ZWLBDR(:,:,:) = 1.E10
  WHERE (PRRS(:,:,:)>0.0.AND.PCRS(:,:,:)>0.0 )
    ZWLBDR3(:,:,:) = XLBR * PCRS(:,:,:) / (PRHODREF(:,:,:) * PRRS(:,:,:))
    ZWLBDR(:,:,:)  = ZWLBDR3(:,:,:)**XLBEXR
  END WHERE
  WHERE (ZWLBDR(:,:,:)<(XACCR1/XSPONBUD1))
    PCRS(:,:,:) = PCRS(:,:,:)*MAX((1.+XSPONCOEF2*(XACCR1/ZWLBDR(:,:,:)-XSPONBUD1)**2),&
                                                 (XACCR1/ZWLBDR(:,:,:)/XSPONBUD3)**3)
  END WHERE

  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2), 'BRKU', pcrs(:, :, :) * prhodj(:, :, :) )
ENDIF
!-------------------------------------------------------------------------------
!*       6.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!	        -------------------------------------
!
!*       6.1 Calculation of the mean volumic radius (ZRAY) and
!             the terminal vertical velocity ZCC for precipitating clouds
!
ZTSPLITR = PTSTEP / REAL(KSPLITR)       ! Small time step
!
!
!*       6.2    compute the sedimentation velocities for rain
!   	        --------------------------------------------
!
ZMVRR(:,:,:) = 0.
ZVRR(:,:,:) = 0.
ZVCR(:,:,:) = 0.
WHERE (PCRT(:,:,:) > XCTMIN(3) .and. PRRT(:,:,:)>XRTMIN(3) )
  ZMVRR(:,:,:) = ((3. * PRHODREF(:,:,:)*PRRT(:,:,:))/   &
                  (4. * XPI *XRHOLW*PCRT(:,:,:)))**0.333 ! in m
  ZVRR(:,:,:) = 0.012 * 1.0E6 * ZMVRR(:,:,:) - 0.2 ! velocity for mixing ratio
  ZVCR(:,:,:) = 0.007 * 1.0E6 * ZMVRR(:,:,:) - 0.1 ! velocity for concentration
END WHERE
WHERE (ZVRR(:,:,:) .lt. 0.0 .OR. ZVCR(:,:,:) .lt. 0.0)
  ZVRR(:,:,:) = 0.0
  ZVCR(:,:,:) = 0.0
END WHERE
!
CALL C2R2_KHKO_SEDIMENTATION
! 
DEALLOCATE(ZRTMIN)
DEALLOCATE(ZCTMIN)            
!
!------------------------------------------------------------------------------
!


!------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
 SUBROUTINE C2R2_KHKO_NUCLEATION
!
!*      0. DECLARATIONS
!          ------------
!JUAN
USE MODI_GET_HALO
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZTCELSIUS
INTEGER , DIMENSION(SIZE(GNUCT))  :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
INTEGER                           :: J1
!
!-------------------------------------------------------------------------------

if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HENU', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'HENU', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'HENU', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_sv ) then
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg    ), 'HENU', pcns(:, :, :) * prhodj(:, :, :) )
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'HENU', pccs(:, :, :) * prhodj(:, :, :) )
end if

! Modification of XCHEN according to theta vertical gradient (J. Rangonio)
!ZZA(:,:,2) = 1.
!DO JK=IKB,IKE-1
! WHERE (PZZ(:,:,JK+1) >= XAERHEIGHT) 
!  ZZA(:,:,JK+1) = ZZA(:,:,JK)
! ELSEWHERE
!  ZZA(:,:,JK+1) = ZZA(:,:,JK)* &
!   EXP(MIN(0.,-XAERDIFF*(PTHT(:,:,JK+1)-PTHT(:,:,JK))/(PZZ(:,:,JK+1)-PZZ(:,:,JK))))
! END WHERE
! ZCHEN(:,:,JK) = XCHEN*ZZA(:,:,JK)
!END DO
!ZCHEN(:,:,IKE) = ZCHEN(:,:,IKE-1)
!!
!!
! IF ( tpfile%lopened ) THEN
!   TZFIELD = TFIELDMETADATA(     &
!     CMNHNAME   = 'ZCHEN',       &
!     CSTDNAME   = '',            &
!     CLONGNAME  = 'ZCHEN',       &
!     CUNITS     = '',            &
!     CDIR       = 'XY',          &
!     CCOMMENT   = 'X_Y_Z_ZCHEN', &
!     NGRID      = 1,             &
!     NTYPE      = TYPEREAL,      &
!     NDIMS      = 3,             &
!     LTIMEDEP   = .TRUE.         )
!   CALL IO_Field_write(TPFILE,TZFIELD,ZCHEN)
! END IF
!
!-------------------------------------------------------------------------------
!
!  compute the saturation vapor mixing ratio  and
!          the radiative tendency and
!          the latent heat of vaporization Lv(T) and 
!          the specific heat for moist air Cph
!  
ZEPS= XMV / XMD
ZRVSAT(:,:,:) = ZEPS / (PPABST(:,:,:) * &
                   EXP(-XALPW+XBETAW/ZT(:,:,:)+XGAMW*ALOG(ZT(:,:,:))) - 1.0)
ZZW1LOG(:,:,:)= 0. ! supersaturation
ZTDT(:,:,:)   = 0.
ZDRC(:,:,:)   = 0.
IF (OACTIT) THEN
  ZTM(:,:,:)    = PTHM(:,:,:) * (PPABSM(:,:,:)/XP00)**(XRD/XCPD)
  ZTDT(:,:,:)   = (ZT(:,:,:)-ZTM(:,:,:))/PTSTEP                              ! dT/dt
  ZDRC(:,:,:)   = (PRCT(:,:,:)-PRCM(:,:,:))/PTSTEP                           ! drc/dt
! Modif M.Mazoyer
! ZTDT(:,:,:)   =  PDTHRAD(:,:,:)*(PPABST(:,:,:)/XP00)**(XRD/XCPD)
END IF
!
!  optimization by looking for locations where
!  the updraft velocity is positive!!!
!
GNUCT(:,:,:) = .FALSE.
IF( OACTIT ) THEN
 GNUCT(IIB:IIE,IJB:IJE,IKB:IKE) = (PW_NU(IIB:IIE,IJB:IJE,IKB:IKE)>XWMIN .OR. &
 ZTDT(IIB:IIE,IJB:IJE,IKB:IKE)<XTMIN)  .AND.   &
 PRVT(IIB:IIE,IJB:IJE,IKB:IKE)>(0.98*ZRVSAT(IIB:IIE,IJB:IJE,IKB:IKE))
ELSE
 GNUCT(IIB:IIE,IJB:IJE,IKB:IKE) = PW_NU(IIB:IIE,IJB:IJE,IKB:IKE)>XWMIN .AND.   &
           PRVT(IIB:IIE,IJB:IJE,IKB:IKE)>(0.98*ZRVSAT(IIB:IIE,IJB:IJE,IKB:IKE))
END IF
INUCT = COUNTJV( GNUCT(:,:,:),I1(:),I2(:),I3(:))

! IF( INUCT >= 1 ) THEN
  ALLOCATE(ZRVT(INUCT))
  ALLOCATE(ZRCT(INUCT))
  ALLOCATE(ZRRT(INUCT))
  ALLOCATE(ZCNS(INUCT))
  ALLOCATE(ZCCS(INUCT))
  ALLOCATE(ZZT(INUCT)) 
  ALLOCATE(ZTDTBIS(INUCT))
  ALLOCATE(ZZW1(INUCT))
  ALLOCATE(ZZW2(INUCT))
  ALLOCATE(ZZW3(INUCT))
  ALLOCATE(ZZW4(INUCT))
  ALLOCATE(ZZW5(INUCT))
  ALLOCATE(ZVEC1(INUCT))
  ALLOCATE(IVEC1(INUCT))
  ALLOCATE(ZRHODREF(INUCT)) 
  ALLOCATE(ZEXNREF(INUCT)) 
!-------------------------------------------------------------------------------
! Modification of XCHEN according to theta vertical gradient (J. Rangonio)
!  ALLOCATE(ZCHEN_TMP(INUCT))
!  ALLOCATE(ZCONC_CCN(INUCT)) 
!-------------------------------------------------------------------------------
  DO JL=1,INUCT
    ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
    ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
    ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
    ZCNS(JL) = PCNS(I1(JL),I2(JL),I3(JL))
    ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
    ZZT(JL)  = ZT(I1(JL),I2(JL),I3(JL))
    ZZW1(JL) = ZRVSAT(I1(JL),I2(JL),I3(JL))
    ZZW2(JL) = PW_NU(I1(JL),I2(JL),I3(JL))
    ZTDTBIS(JL) = ZTDT(I1(JL),I2(JL),I3(JL))
    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
    ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
!-------------------------------------------------------------------------------
! Modification of XCHEN according to theta vertical gradient (J. Rangonio)
!  ZCHEN_TMP(JL)= ZCHEN(I1(JL),I2(JL),I3(JL))
!-------------------------------------------------------------------------------
    
  ENDDO
!-------------------------------------------------------------------------------
! Modification of XCHEN according to theta vertical gradient (J. Rangonio)
! ZCONC_CCN(:)=XCONC_CCN*ZCHEN_TMP(:)/XCHEN
!-------------------------------------------------------------------------------
  ZZW1(:) = 1.0/ZEPS + 1.0/ZZW1(:)                                   &
          + (((XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/ZZT(:))**2)/(XCPD*XRV) ! Psi2
!
!*       3.1     compute the heterogeneous nucleation source: RVHENC, CVHENC
!
!*       3.1.1   compute the constant term (ZZW3)
!
  ZVEC1(:) = MAX( 1.00001, MIN( REAL(NAHEN)-0.00001, &
                  XAHENINTP1 * ZZT(:) + XAHENINTP2 )  )
  IVEC1(:) = INT( ZVEC1(:) )
  ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
  ALLOCATE(ZSMAX(INUCT))
!
!
  IF( HPARAM_CCN == 'TFH' ) THEN
    ZZW2(:) = 100.*ZZW2(:) ! FH is in CGS units
    ALLOCATE(ZTCELSIUS(INUCT)); ZTCELSIUS(:) = ZZT(:) - XTT
    ZZW3(:) =   XAHENF( IVEC1(:)+1 )* ZVEC1(:)      &
              - XAHENF( IVEC1(:)   )*(ZVEC1(:) - 1.0)      ! Cste*(Psi1/Gr)
    ZZW3(:) = ZZW3(:)/ZZW2(:)**(XWCOEF_F1+ZTCELSIUS(:)*  &
                     (XWCOEF_F2+XWCOEF_F3*ZTCELSIUS(:)))
    ZZW3(:) = (ZZW3(:)/ZZW1(:)) * ZZW2(:) * ZRHODREF(:) ! R.H.S. of
                                                        ! Eq. (12) in FH92
!
!*       3.1.1.1   compute the maximum fo supersaturation
!
    ZSMAX(:) = ZZW3(:)**(1.0/(XKHEN+1.0)) ! first estimate (y_bar=0)
!
! 4 iterations to estimate S_max for the TFH parameterization
!
    ZZW1(:) =   XAHENY( IVEC1(:)+1 )* ZVEC1(:)      &
              - XAHENY( IVEC1(:)   )*(ZVEC1(:) - 1.0)     ! y_bar
    ZZW1(:) = ZZW1(:)*ZZW2(:)** (XWCOEF_Y1+ZTCELSIUS(:)*  &
                      (XWCOEF_Y2+XWCOEF_Y3*ZTCELSIUS(:)))
    DO J1 = 1,4
      ZSMAX(:) = (ZZW1(:)*ZSMAX(:)**XKHEN + ZSMAX(:))**(1.0/(XKHEN+1.0))
    END DO
    DEALLOCATE(ZTCELSIUS)
    ZZW3(:) = 1.0
  ELSE
    IF (OACTIT) THEN
      ZZW4(:)=XPSI1( IVEC1(:)+1)*ZZW2(:)+XPSI3(IVEC1(:)+1)*ZTDTBIS(:)
      ZZW5(:)=XPSI1( IVEC1(:))*ZZW2(:)+XPSI3(IVEC1(:))*ZTDTBIS(:)
! Modif M.Mazoyer
!     ZZW4(:) =0.0
!     ZZW5(:) =0.0
!     WHERE  (ZZW2(:)>= XWMIN .AND. ZTDTBIS(:) < XTMIN )
!     ZZW4(:)=XPSI1( IVEC1(:)+1)*ZZW2(:)+XPSI3(IVEC1(:)+1)*ZTDTBIS(:)
!     ZZW5(:)=XPSI1( IVEC1(:))*ZZW2(:)+XPSI3(IVEC1(:))*ZTDTBIS(:)
!     ELSEWHERE  (ZZW2(:)< XWMIN .AND. ZTDTBIS(:) < XTMIN )
!     ZZW4(:)=XPSI3(IVEC1(:)+1)*ZTDTBIS(:)
!     ZZW5(:)=XPSI3(IVEC1(:))*ZTDTBIS(:)
!     ELSEWHERE  (ZZW2(:)< XWMIN .AND. ZTDTBIS(:) >= XTMIN  )  
!     ZZW4(:)=0.0
!     ZZW5(:)=0.0
!     ELSEWHERE  (ZZW2(:)>= XWMIN .AND. ZTDTBIS(:) >= XTMIN  ) 
!     ZZW4(:)=XPSI1( IVEC1(:)+1)*ZZW2(:)
!     ZZW5(:)=XPSI1( IVEC1(:))*ZZW2(:)
!     END WHERE
      WHERE (ZZW4(:) < 0. .OR. ZZW5(:) < 0.)
             ZZW4(:) = 0.
             ZZW5(:) = 0.
      END WHERE
      ZZW3(:) = XCHEN*XAHENG(IVEC1(:)+1)*(ZZW4(:)**1.5)*ZVEC1(:)/XCHEN      &
              - XCHEN*XAHENG( IVEC1(:))*(ZZW5(:)**1.5)*(ZVEC1(:) - 1.0)/XCHEN
                       ! Cste*((Psi1*w+Psi3*dT/dt)/(G))**1.5
!-------------------------------------------------------------------------------
! Modification of XCHEN according to theta vertical gradient (J. Rangonio)
!      ZZW3(:) = XCHEN*XAHENG(IVEC1(:)+1)*(ZZW4(:)**1.5)*ZVEC1(:)/ZCHEN_TMP(:)      &
!              - XCHEN*XAHENG( IVEC1(:))*(ZZW5(:)**1.5)*(ZVEC1(:) - 1.0)/ZCHEN_TMP(:)  
!                       ! Cste*((Psi1*w+Psi3*dT/dt)/(G))**1.5
!-------------------------------------------------------------------------------
    ELSE
      ZZW3(:) = XAHENG( IVEC1(:)+1)*((XPSI1( IVEC1(:)+1)*ZZW2(:))**1.5)* ZVEC1(:)      &
              - XAHENG( IVEC1(:))*((XPSI1(IVEC1(:))*ZZW2(:))**1.5)*(ZVEC1(:) - 1.0)
    END IF
    ZZW5(:) = 1.
    ZZW3(:) = (ZZW3(:)/ZZW1(:))*ZRHODREF(:) ! R.H.S. of
                                                      ! Eq 9 of CPB 98 
    WHERE (ZZW3(:) == 0.)
            ZZW5(:)= -1.
    END WHERE
!
!*       3.1.2.1   compute the maximum fo supersaturation
!
    ZSMAX(:) = ZZW3(:)**(1.0/(XKHEN+2.0)) ! Smax has no unit
!
! 4 iterations to estimate S_max for the CPB98 parameterization
!
    IF( HPARAM_CCN == 'CPB' ) THEN
      DO J1 = 1,4
       WHERE (ZZW5(:) > 0.)
        ZVEC1(:) = MAX( 1.00001, MIN( REAL(NHYP)-0.00001,      &
                        XHYPINTP1*LOG(ZSMAX(:))+XHYPINTP2 ) )
        IVEC1(:) = INT( ZVEC1(:) )
        ZVEC1(:) = ZVEC1(:) - REAL( IVEC1(:) )
        ZZW2(:)  =   XHYPF32( IVEC1(:)+1 )* ZVEC1(:)      &
                   - XHYPF32( IVEC1(:)   )*(ZVEC1(:) - 1.0)
        ZSMAX(:) = (ZZW3(:)/ZZW2(:))**(1.0/(XKHEN+2.0))
       ELSEWHERE
         ZSMAX(:)=0.
       END WHERE
      END DO
!
!*       3.2    compute the nucleus source
!
! ZSMAX(:) is used in percent in the nucleation formula
!
      ZZW3(:) =   XHYPF12( IVEC1(:)+1 )* ZVEC1(:)      &
                - XHYPF12( IVEC1(:)   )*(ZVEC1(:) - 1.0)
    ELSE
      ZZW3(:) = 1.0
    END IF
  END IF
  ZZW1LOG(:,:,:) = UNPACK( 100*ZSMAX(:),MASK=GNUCT(:,:,:),FIELD=0.0 )
  PSUPSAT(:,:,:) = 0.0
  PSUPSAT(:,:,:) = ZZW1LOG(:,:,:)
!
! the CCN spectra formula uses ZSMAX in percent
!
  IF (XCONC_CCN > 0) THEN
    ZZW1(:) = MIN( XCONC_CCN,XCHEN * (100.0*ZSMAX(:))**XKHEN * ZZW3(:) ) / PTSTEP
  ELSE
    ZZW1(:) = XCHEN * (100.0*ZSMAX(:))**XKHEN * ZZW3(:) / PTSTEP
  ENDIF
!-------------------------------------------------------------------------------
! Modification of XCHEN according to theta vertical gradient (J. Rangonio)
!  IF (XCONC_CCN > 0.) THEN
!    ZZW1(:) = MIN( ZCONC_CCN(:),ZCHEN_TMP(:) * (100.0*ZSMAX(:))**XKHEN * ZZW3(:) ) / PTSTEP
!  ELSE
!    ZZW1(:) = ZCHEN_TMP(:) * (100.0*ZSMAX(:))**XKHEN * ZZW3(:) / PTSTEP
!  ENDIF
!-------------------------------------------------------------------------------
  ZW(:,:,:)   = PCNS(:,:,:)
  PCNS(:,:,:) = UNPACK( MAX( ZZW1(:),ZCNS(:) ),MASK=GNUCT(:,:,:), &
                                                 FIELD=ZW(:,:,:)  )
!
  DEALLOCATE(IVEC1)
  DEALLOCATE(ZVEC1)
!
!*       3.3    compute the cloud water concentration and mixing ratio sources
!
  ZZW2(:) = MAX( (ZZW1(:)-ZCNS(:)),0.0 )
  
  PNACT(:,:,:) = 0.0
  PNACT(:,:,:) = UNPACK(ZZW2(:)*PTSTEP,MASK=GNUCT(:,:,:),FIELD=0.)
  
  ZZW1(:)=0.
  WHERE (ZZW5(:) > 0.)
    ZZW1(:) = MIN( XCSTDCRIT * ZZW2(:) / ( ((ZZT(:)*ZSMAX(:))**3.)*ZRHODREF(:) ),&
                1.E-5 )
  END WHERE
  CALL GET_HALO(PRVS)
  ZW(:,:,:) = MIN( UNPACK( ZZW1(:),MASK=GNUCT(:,:,:),FIELD=0.0 ),PRVS(:,:,:) )
!
  PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
  PRCS(:,:,:) = PRCS(:,:,:) + ZW(:,:,:) 
  ZW(:,:,:) = ZW(:,:,:)*(XLVTT+(XCPV-XCL)*(ZT(:,:,:)-XTT))/                      &
                   (PEXNREF(:,:,:)*( XCPD+XCPV*PRVT(:,:,:)+XCL*(PRCT(:,:,:)+PRRT(:,:,:))))
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)
!JUAN
  CALL GET_HALO(PTHS)
  CALL GET_HALO(PRCS)
!  
  ZW(:,:,:)   = PCCS(:,:,:)
  PCCS(:,:,:) = UNPACK( ZZW2(:)+ZCCS(:),MASK=GNUCT(:,:,:),FIELD=ZW(:,:,:) )
!
!
  DEALLOCATE(ZRVT)
  DEALLOCATE(ZRCT)
  DEALLOCATE(ZRRT)
  DEALLOCATE(ZCNS)
  DEALLOCATE(ZCCS)
  DEALLOCATE(ZZT)
  DEALLOCATE(ZSMAX)
  DEALLOCATE(ZZW1)
  DEALLOCATE(ZZW2)
  DEALLOCATE(ZZW3)
  DEALLOCATE(ZZW4)
  DEALLOCATE(ZZW5)
  DEALLOCATE(ZTDTBIS)
  DEALLOCATE(ZRHODREF)
  DEALLOCATE(ZEXNREF)
!-------------------------------------------------------------------------------
! Modification of XCHEN according to theta vertical gradient (J. Rangonio)
!  DEALLOCATE(ZCHEN_TMP)
!  DEALLOCATE(ZCONC_CCN)
!-------------------------------------------------------------------------------
! END IF
!                      
IF ( tpfile%lopened ) THEN
  TZFIELD = TFIELDMETADATA(    &
    CMNHNAME   = 'SMAX',       &
    CSTDNAME   = '',           &
    CLONGNAME  = 'SMAX',       &
    CUNITS     = '1',          &
    CDIR       = 'XY',         &
    CCOMMENT   = 'X_Y_Z_SMAX', &
    NGRID      = 1,            &
    NTYPE      = TYPEREAL,     &
    NDIMS      = 3,            &
    LTIMEDEP   = .TRUE.        )
  CALL IO_Field_write(TPFILE,TZFIELD,ZZW1LOG)
END IF
!
!*       3.4   budget storage
!
if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HENU', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'HENU', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'HENU', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_sv ) then
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg    ), 'HENU', pcns(:, :, :) * prhodj(:, :, :) )
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'HENU', pccs(:, :, :) * prhodj(:, :, :) )
end if

END SUBROUTINE C2R2_KHKO_NUCLEATION
!
!-------------------------------------------------------------------------------
!
 SUBROUTINE AER_NUCLEATION
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_NSV
USE MODE_AERO_PSD
USE MODI_CH_AER_ACTIVATION

IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZTCELSIUS
INTEGER , DIMENSION(SIZE(GNUCT))  :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
INTEGER                           :: J1
INTEGER                           :: JSV
!
!-------------------------------------------------------------------------------

if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH), 'HENU', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV), 'HENU', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'HENU', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_sv ) then
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg    ), 'HENU', pcns(:, :, :) * prhodj(:, :, :) )
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'HENU', pccs(:, :, :) * prhodj(:, :, :) )
end if
!
!  compute the saturation vapor mixing ratio  
!          the radiative tendency                    
!
ZEPS= XMV / XMD
!
!
ZRVSAT(:,:,:) = ZEPS / (PPABST(:,:,:) * &
                   EXP(-XALPW+XBETAW/ZT(:,:,:)+XGAMW*ALOG(ZT(:,:,:))) - 1.0)
ZZW1LOG(:,:,:)= 0. ! supersaturation
ZTDT(:,:,:)   = 0.
ZDRC(:,:,:)   = 0.
IF (OACTIT) THEN
  ZTM(:,:,:)    = PTHM(:,:,:) * (PPABSM(:,:,:)/XP00)**(XRD/XCPD)
  ZTDT(:,:,:)   = (ZT(:,:,:)-ZTM(:,:,:))/PTSTEP                              ! dT/dt
  ZDRC(:,:,:)   = (PRCT(:,:,:)-PRCM(:,:,:))/PTSTEP                           ! drc/dt
  ZTDT(:,:,:)   = MIN(0.,ZTDT(:,:,:)+(XG*PW_NU(:,:,:))/XCPD- &
  (XLVTT+(XCPV-XCL)*(ZT(:,:,:)-XTT))*ZDRC(:,:,:)/XCPD)
! Modif M.Mazoyer
! ZTDT(:,:,:)   =  PDTHRAD(:,:,:)*(PPABST(:,:,:)/XP00)**(XRD/XCPD)
END IF
!
!  optimization by looking for locations where
!  the updraft velocity is positive!!!
!
GNUCT(:,:,:) = .FALSE.
IF( OACTIT ) THEN
 GNUCT(IIB:IIE,IJB:IJE,IKB:IKE) = (PW_NU(IIB:IIE,IJB:IJE,IKB:IKE)>XWMIN .OR. &
 ZTDT(IIB:IIE,IJB:IJE,IKB:IKE)<XTMIN)  .AND.   &
 PRVT(IIB:IIE,IJB:IJE,IKB:IKE)>(0.98*ZRVSAT(IIB:IIE,IJB:IJE,IKB:IKE))
ELSE
 GNUCT(IIB:IIE,IJB:IJE,IKB:IKE) = PW_NU(IIB:IIE,IJB:IJE,IKB:IKE)>XWMIN .AND.   &
           PRVT(IIB:IIE,IJB:IJE,IKB:IKE)>(0.98*ZRVSAT(IIB:IIE,IJB:IJE,IKB:IKE))
END IF
!
INUCT = COUNTJV(GNUCT(:,:,:),I1(:),I2(:),I3(:))

IF( INUCT >= 1 ) THEN
  ALLOCATE(ZRVT(INUCT))
  ALLOCATE(ZRCT(INUCT))
  ALLOCATE(ZRRT(INUCT))
  ALLOCATE(ZZT(INUCT))
  ALLOCATE(ZTDTBIS(INUCT)) 
  ALLOCATE(ZZW1(INUCT))
  ALLOCATE(ZZW2(INUCT))
  ALLOCATE(ZZW3(INUCT))
  ALLOCATE(ZZW4(INUCT))
  ALLOCATE(ZZW5(INUCT))
  ALLOCATE(ZDG3(INUCT))
  ALLOCATE(ZCCS(INUCT))
  ALLOCATE(ZCNS(INUCT))
  ALLOCATE(ZRHODREF(INUCT)) 
  ALLOCATE(ZEXNREF(INUCT)) 
  ALLOCATE(ZPABST(INUCT)) 
  ALLOCATE(ZNCN(INUCT)) 
  ALLOCATE(ZMCN(INUCT)) 
  ALLOCATE(ZSMAX(INUCT))
  ALLOCATE(ZAERO(INUCT,SIZE(PAEROT,4)))
  ALLOCATE(ZSOLORG(INUCT,SIZE(PSOLORG,4))) 
  ALLOCATE(ZMI(INUCT,SIZE(PMI,4)))
  ALLOCATE(ZLBDC3(INUCT)) 

  DO JL=1,INUCT
    ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
    ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
    ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
    ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
    ZCNS(JL) = PCNS(I1(JL),I2(JL),I3(JL))
    ZZT(JL)  = ZT(I1(JL),I2(JL),I3(JL))
    ZZW1(JL) = ZRVSAT(I1(JL),I2(JL),I3(JL))
    ZZW2(JL) = PW_NU(I1(JL),I2(JL),I3(JL))
    ZTDTBIS(JL) = ZTDT(I1(JL),I2(JL),I3(JL))
    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
    ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
    ZPABST(JL)   = PPABST(I1(JL),I2(JL),I3(JL))
    ZAERO(JL,:)   = PAEROT(I1(JL),I2(JL),I3(JL),:)
    ZLBDC3(JL) = ZWLBDC3(I1(JL),I2(JL),I3(JL))

  ENDDO
!
  ZSMAX(:) = 0.
  IF (LORILAM) THEN
    DO JL=1,INUCT
      ZSOLORG(JL,:) = PSOLORG(I1(JL),I2(JL),I3(JL),:)
      ZMI(JL,:) = PMI(I1(JL),I2(JL),I3(JL),:)
    ENDDO
  ELSE
    ZSOLORG(:,:) = 0.
    ZMI(:,:) = 0.
  END IF

  CALL CH_AER_ACTIVATION(ZAERO, ZZT, ZZW2, ZTDTBIS, ZRHODREF, ZPABST,&
                       ZNCN, ZMCN, ZSOLORG, ZMI, ZSMAX)

! Nb de goutelettes activées

!test
  ZZW1(:) = MAX(ZNCN(:)/PTSTEP - ZCNS(:), 0.)
!
  ZW(:,:,:) = UNPACK( ZZW1(:),MASK=GNUCT(:,:,:),FIELD=0.0 ) 
  PCNS(:,:,:) = PCNS(:,:,:) + ZW(:,:,:)
!
! Modification reservoir eau (gaz et liquide)
!
! valeur de petites goutelettes type brouillard (test)
!  ZALPHA=0.8
!  ZMU=3.
!  ZDG3(:) =  1./ZLBDC3(:) * GAMMA(ZMU + 3./ZALPHA) / GAMMA(ZMU) ! integrated cubic diameter
!  ZZW2(:) = ZZW1(:) + ZCCS(:)
!  ZZW1(:) = XPI/6. * ZDG3(:)**3 * (ZZW1(:)) * 1000. /  ZRHODREF(:)
!  !
!  ZW(:,:,:) = MIN( UNPACK( ZZW1(:),MASK=GNUCT(:,:,:),FIELD=0.0 ),PRVS(:,:,:) )
!  !
!
!  PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
!  PRCS(:,:,:) = PRCS(:,:,:) + ZW(:,:,:) 
!  !
!  ! Modification temperature (diabatisme)
!  ZZW1(:) = ZZW1(:)*(XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))/                      &
!                      (ZEXNREF(:)*( XCPD+XCPV*ZRVT(:)+XCL*(ZRCT(:)+ZRRT(:))))
!  !
!  ZW(:,:,:) = MIN( UNPACK( ZZW1(:),MASK=GNUCT(:,:,:),FIELD=0.0 ),PRVS(:,:,:) )
!  !
!  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)
!  !
!  ! Modification gouttes nuages
!  ZW(:,:,:)   = PCCS(:,:,:)
!  PCCS(:,:,:) = UNPACK(ZZW2(:),MASK=GNUCT(:,:,:),FIELD=ZW(:,:,:))
  ZZW2(:) = MAX(ZNCN(:)/PTSTEP - ZCNS(:), 0.)
  ZZW1(:)=0.
  WHERE(ZZW2(:).gt.0.0)
    ZZW1(:)=MIN(XCSTDCRIT * ZZW2(:) / ( ((ZZT(:)*ZSMAX(:))**3.)&
    *ZRHODREF(:) ) , 1.E-5 )   
  END WHERE
  ZW(:,:,:) = MIN( UNPACK( ZZW1(:),MASK=GNUCT(:,:,:),FIELD=0.0 ),PRVS(:,:,:) )
  PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
  PRCS(:,:,:) = PRCS(:,:,:) + ZW(:,:,:) 
!
  DEALLOCATE(ZRVT)
  DEALLOCATE(ZRCT)
  DEALLOCATE(ZRRT)
  DEALLOCATE(ZZT)
  DEALLOCATE(ZTDTBIS) 
  DEALLOCATE(ZZW1)
  DEALLOCATE(ZZW2)
  DEALLOCATE(ZZW3)
  DEALLOCATE(ZZW4)
  DEALLOCATE(ZZW5)
  DEALLOCATE(ZDG3)
  DEALLOCATE(ZCCS)
  DEALLOCATE(ZCNS)
  DEALLOCATE(ZRHODREF) 
  DEALLOCATE(ZEXNREF) 
  DEALLOCATE(ZPABST) 
  DEALLOCATE(ZNCN) 
  DEALLOCATE(ZMCN) 
  DEALLOCATE(ZAERO) 
  DEALLOCATE(ZSMAX)
  DEALLOCATE(ZSOLORG) 
  DEALLOCATE(ZMI)
  DEALLOCATE(ZLBDC3) 

END IF
!
!
!*             budget storage
!
!
if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH), 'HENU', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV), 'HENU', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'HENU', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_sv ) then
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg    ), 'HENU', pcns(:, :, :) * prhodj(:, :, :) )
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'HENU', pccs(:, :, :) * prhodj(:, :, :) )
end if

  END SUBROUTINE AER_NUCLEATION
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE C2R2_COALESCENCE
!
!
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
INTEGER , DIMENSION(SIZE(GNUCT))  :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
!
!-------------------------------------------------------------------------------
!
!
GMICRO(:,:,:) = .FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =               &
  PRCT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(2) .OR.  &
  PRRT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(3)
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
IF( IMICRO >= 1 ) THEN
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
!
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
!*       4.1   Self-collection of cloud droplets
!
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'SELF', pccs(:, :, :) * prhodj(:, :, :) )

  GSELF(:) = ZCCT(:)>XCTMIN(2)
  ISELF = COUNT(GSELF(:))
  IF( ISELF>0 ) THEN
    ZZW1(:) = XSELFC*(ZCCT(:)/ZLBDC3(:))**2 ! analytical integration
      WHERE( GSELF(:) )
        ZCCS(:) = ZCCS(:) - MIN( ZCCS(:),ZZW1(:) )
      END WHERE
  END IF

  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'SELF', &
                                           Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
!
!*       4.2   Autoconversion of cloud droplets
!              using a Berry-Reinhardt parameterization
!
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC),                        'AUTO', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RR),                        'AUTO', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2), 'AUTO', pcrs(:, :, :) * prhodj(:, :, :) )

  ZZW2(:) = 0.0
  ZZW1(:) = 0.0
  WHERE( ZRCT(:)>XRTMIN(2) )
    ZZW2(:) = MAX( 0.0,XLAUTR*ZRHODREF(:)*ZRCT(:)*             &
                       (XAUTO1/ZLBDC(:)**4-XLAUTR_THRESHOLD) ) ! L 
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
    ZZW3(:) = ZZW3(:) * ZRHODREF(:)**2 * MAX( 0.0,ZZW1(:) )**3 / XAC 
    ZCRS(:) = ZCRS(:) + ZZW3(:)
  END WHERE

  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'AUTO', &
                                           Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'AUTO', &
                                           Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2), 'AUTO', &
                                           Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )
!
!
!*       4.3    Accretion sources
!
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'ACCR', &
                                            Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR), 'ACCR', &
                                            Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'ACCR', &
                                            Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
!
!*       4.31   test the criterium Df>Dh or Nr>Nrm
!
  GACCR(:) = ZRRT(:)>XRTMIN(3) .AND. ZCRT(:)>XCTMIN(3)
  IACCR = COUNT(GACCR(:))
  IF( IACCR>0 ) THEN
    ALLOCATE(ZZW4(IMICRO)); ZZW4(:) = XACCR1/ZLBDR(:)
    ALLOCATE(GENABLE_ACCR_SCBU(IMICRO))
    GENABLE_ACCR_SCBU(:) = ZRRT(:)>1.2*ZZW2(:)/ZRHODREF(:) .OR.           &
                     ZZW4(:)>=MAX( XACCR2,XACCR3/(XACCR4/ZLBDC(:)-XACCR5) )
    GACCR(:) = GACCR(:) .AND. ZRCT(:)>XRTMIN(2) .AND. GENABLE_ACCR_SCBU(:)
  END IF
!
  IACCR = COUNT(GACCR(:))
  IF( IACCR>0 ) THEN
    WHERE( GACCR(:).AND.(ZZW4(:)>1.E-4) ) ! Accretion for D>100 10-6 m
      ZZW3(:) = ZLBDC3(:) / ZLBDR3(:)
      ZZW1(:) = ZCCT(:)*ZCRT(:) / ZLBDC3(:)
      ZZW2(:) = MIN( ZZW1(:)*(XACCR_CLARGE1+XACCR_CLARGE2*ZZW3(:)),ZCCS(:) )
      ZCCS(:) = ZCCS(:) - ZZW2(:)
!
      ZZW1(:) = ZZW1(:) / ZLBDC3(:)
      ZZW2(:) = MIN( ZZW1(:)*(XACCR_RLARGE1+XACCR_RLARGE2*ZZW3(:))             &
                                                          /ZRHODREF(:),ZRCS(:) )
      ZRCS(:) = ZRCS(:) - ZZW2(:)
      ZRRS(:) = ZRRS(:) + ZZW2(:)
    END WHERE
    WHERE( GACCR(:).AND.(ZZW4(:)<=1.E-4) ) ! Accretion for D<100 10-6 m
      ZZW3(:) = ZLBDC3(:) / ZLBDR3(:)
      ZZW1(:) = ZCCT(:)*ZCRT(:) / ZLBDC3(:)**2
      ZZW3(:) = ZZW3(:)**2
      ZZW2(:) = MIN( ZZW1(:)*(XACCR_CSMALL1+XACCR_CSMALL2*ZZW3(:)),ZCCS(:) )
      ZCCS(:) = ZCCS(:) - ZZW2(:)
!
      ZZW1(:) = ZZW1(:) / ZLBDC3(:)
      ZZW2(:) = MIN( ZZW1(:)*(XACCR_RSMALL1+XACCR_RSMALL2*ZZW3(:))             &
                                                          /ZRHODREF(:),ZRCS(:) )
      ZRCS(:) = ZRCS(:) - ZZW2(:)
      ZRRS(:) = ZRRS(:) + ZZW2(:)
    END WHERE
  END IF

  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'ACCR', &
                                           Unpack( zrcs(:), mask = gmicro(:, :, :), field = prcs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR), 'ACCR', &
                                           Unpack( zrrs(:), mask = gmicro(:, :, :), field = prrs(:, :, :) ) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'ACCR', &
                                           Unpack( zccs(:), mask = gmicro(:, :, :), field = pccs(:, :, :) ) * prhodj(:, :, :) )
!
!*       4.4   Self collection - Coalescence/Break-up
!
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2), 'SCBU', &
                                            Unpack( zcrs(:), mask = gmicro(:, :, :), field = pcrs(:, :, :) ) * prhodj(:, :, :) )

  IF( IACCR>0 ) THEN
    GSCBU(:) = ZCRT(:)>XCTMIN(3) .AND. GENABLE_ACCR_SCBU(:)
    ISCBU = COUNT(GSCBU(:))
  ELSE
    ISCBU = 0.0
  END IF
  IF( ISCBU>0 ) THEN
!
!*       4.41  efficiencies
!
    IF (.NOT.ALLOCATED(ZZW4)) ALLOCATE(ZZW4(IMICRO))
    ZZW4(:)  = XACCR1 / ZLBDR(:)                ! Mean diameter
    ALLOCATE(ZSCBU(IMICRO))
    ZSCBU(:) = 1.0
    WHERE (ZZW4(:)>=XSCBU_EFF1 .AND. GSCBU(:))   ZSCBU(:) = &  ! Coalescence
                          EXP(XSCBUEXP1*(ZZW4(:)-XSCBU_EFF1))  ! efficiency
    WHERE (ZZW4(:)>=XSCBU_EFF2) ZSCBU(:) = 0.0  ! Break-up
!
!*       4.42  integration
!
    ZZW1(:) = 0.0
    ZZW2(:) = 0.0
    ZZW3(:) = 0.0
    ZZW4(:) = XACCR1 / ZLBDR(:)                 ! Mean volume drop diameter
    WHERE (GSCBU(:).AND.(ZZW4(:)>1.E-4))              ! analytical integration
      ZZW1(:) = XSCBU2 *  ZCRT(:)**2 / ZLBDR3(:)      ! D>100 10-6 m
      ZZW3(:) = ZZW1(:)*ZSCBU(:)
    END WHERE
    WHERE (GSCBU(:).AND.(ZZW4(:)<=1.E-4))
      ZZW2(:) = XSCBU3 * (ZCRT(:)    / ZLBDR3(:))**2  ! D<100 10-6 m
      ZZW3(:) = ZZW2(:)
    END WHERE
    ZCRS(:) = ZCRS(:) - MIN( ZCRS(:),ZZW3(:) )
    DEALLOCATE(ZSCBU)
  END IF
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

  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2), 'SCBU', pcrs(:, :, :) * prhodj(:, :, :) )

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
  IF( ALLOCATED(IVEC1) ) THEN
    DEALLOCATE(IVEC1)
    DEALLOCATE(ZVEC1)
  END IF
END IF
!
  END SUBROUTINE C2R2_COALESCENCE
!
!-------------------------------------------------------------------------------
!
SUBROUTINE KHKO_COALESCENCE
!
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
INTEGER , DIMENSION(SIZE(GNUCT))  :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
!
!-------------------------------------------------------------------------------
!
!
GMICRO(:,:,:) = .FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =               &
  PRCT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(2) .OR.  &
  PRRT(IIB:IIE,IJB:IJE,IKB:IKE)>XRTMIN(3)
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
IF( IMICRO >= 1 ) THEN
  ALLOCATE(ZRCT(IMICRO))
  ALLOCATE(ZRRT(IMICRO))
  ALLOCATE(ZCCT(IMICRO))
!
  ALLOCATE(ZRCS(IMICRO))
  ALLOCATE(ZRRS(IMICRO))
  ALLOCATE(ZCCS(IMICRO))
  ALLOCATE(ZCRS(IMICRO))
!
  ALLOCATE(ZRHODREF(IMICRO))
!
  DO JL=1,IMICRO
    ZCCT(JL) = PCCT(I1(JL),I2(JL),I3(JL))
    ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
    ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
    ZCCS(JL) = PCCS(I1(JL),I2(JL),I3(JL))
    ZRCS(JL) = PRCS(I1(JL),I2(JL),I3(JL))
    ZRRS(JL) = PRRS(I1(JL),I2(JL),I3(JL))
    ZCRS(JL) = PCRS(I1(JL),I2(JL),I3(JL))
    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
  END DO
!
  ALLOCATE(ZZW1(IMICRO))
!
!*       4.1.1   autoconversion
!
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'SELF', pccs(:, :, :) * prhodj(:, :, :) )

  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC),                        'AUTO', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR),                        'AUTO', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2), 'AUTO', pcrs(:, :, :) * prhodj(:, :, :) )

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
!
  ZW(:,:,:) = PRCS(:,:,:)
  PRCS(:,:,:) = UNPACK( ZRCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PCCS(:,:,:)
  PCCS(:,:,:) = UNPACK( ZCCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PRRS(:,:,:)
  PRRS(:,:,:) = UNPACK( ZRRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PCRS(:,:,:)
  PCRS(:,:,:) = UNPACK( ZCRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
!*       4.1.2   budget storage
!
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'SELF', pccs(:, :, :) * prhodj(:, :, :) )

  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC),                        'AUTO', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR),                        'AUTO', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2), 'AUTO', pcrs(:, :, :) * prhodj(:, :, :) )
!
!*       4.2.1    Accretion sources
!
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC),                        'ACCR', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR),                        'ACCR', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'ACCR', pccs(:, :, :) * prhodj(:, :, :) )

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
!
  ZW(:,:,:) = PRCS(:,:,:)
  PRCS(:,:,:) = UNPACK( ZRCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PCCS(:,:,:)
  PCCS(:,:,:) = UNPACK( ZCCS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PRRS(:,:,:)
  PRRS(:,:,:) = UNPACK( ZRRS(:),MASK=GMICRO(:,:,:),FIELD=ZW(:,:,:) )
!
  DEALLOCATE(ZRCT)
  DEALLOCATE(ZRRT)
  DEALLOCATE(ZCCT)
  DEALLOCATE(ZRCS)
  DEALLOCATE(ZRRS)
  DEALLOCATE(ZCRS)
  DEALLOCATE(ZCCS)
  DEALLOCATE(ZRHODREF) 
  DEALLOCATE(ZZW1)
!
!*       4.2.2   budget storage
!
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC),                        'ACCR', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR),                        'ACCR', prrs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'ACCR', pccs(:, :, :) * prhodj(:, :, :) )
END IF
!
  END SUBROUTINE KHKO_COALESCENCE
!
!------------------------------------------------------------------------------
!
 SUBROUTINE C2R2_KHKO_EVAPORATION
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
INTEGER , DIMENSION(SIZE(GNUCT))  :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
!
!-------------------------------------------------------------------------------

if ( lbudget_th ) call Budget_store_init( tbudgets(NBUDGET_TH),                        'REVA', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_init( tbudgets(NBUDGET_RV),                        'REVA', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rr ) call Budget_store_init( tbudgets(NBUDGET_RR),                        'REVA', prrs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2), 'REVA', pcrs(:, :, :) * prhodj(:, :, :) )
!
!  optimization by looking for locations where
!  the raindrop mixing ratio is non-zero
!
ZW(:,:,:) = 0.0
ZLV(:,:,:) = XLVTT + (XCPV-XCL)*(ZT(:,:,:)-XTT)   !!!latent heat of vaporization
!
GEVAP(:,:,:) = .FALSE.
IF (HCLOUD=='C2R2'.OR. HCLOUD=='C3R5') THEN
  GEVAP(IIB:IIE,IJB:IJE,IKB:IKE) =                              &
    PRRS(IIB:IIE,IJB:IJE,IKB:IKE)> 0.0 .AND.                    &
    PRVT(IIB:IIE,IJB:IJE,IKB:IKE)<ZRVSAT(IIB:IIE,IJB:IJE,IKB:IKE)
ELSE ! KHKO
  GEVAP(IIB:IIE,IJB:IJE,IKB:IKE) =                              &
     PRRS(IIB:IIE,IJB:IJE,IKB:IKE)> 0.0 .AND.                    &
     PCRS(IIB:IIE,IJB:IJE,IKB:IKE)> 0.0 .AND.                    &
     PRRT(IIB:IIE,IJB:IJE,IKB:IKE)> 0.0 .AND.                    &
     PCRT(IIB:IIE,IJB:IJE,IKB:IKE)> 0.0 .AND.                    &
     PRVT(IIB:IIE,IJB:IJE,IKB:IKE)<ZRVSAT(IIB:IIE,IJB:IJE,IKB:IKE)
ENDIF
IEVAP = COUNTJV( GEVAP(:,:,:),I1(:),I2(:),I3(:))

IF( IEVAP >= 1 ) THEN
  ALLOCATE(ZRVT(IEVAP))
  ALLOCATE(ZRCT(IEVAP))
  ALLOCATE(ZRRT(IEVAP))
  ALLOCATE(ZCRT(IEVAP))
  ALLOCATE(ZRVS(IEVAP))
  ALLOCATE(ZRRS(IEVAP))
  ALLOCATE(ZCRS(IEVAP))
  ALLOCATE(ZTHS(IEVAP))
  ALLOCATE(ZLBDR(IEVAP))
  ALLOCATE(ZRHODREF(IEVAP))
  ALLOCATE(ZEXNREF(IEVAP))
  ALLOCATE(ZZT(IEVAP))  
  ALLOCATE(ZZLV(IEVAP)) 
  ALLOCATE(ZZW1(IEVAP))
  ALLOCATE(ZZW2(IEVAP))
  ALLOCATE(ZZW3(IEVAP))
!
  DO JL=1,IEVAP
    ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
    ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
    ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
    ZCRT(JL) = PCRT(I1(JL),I2(JL),I3(JL))
    ZRRS(JL) = PRRS(I1(JL),I2(JL),I3(JL))
    ZRVS(JL) = PRVS(I1(JL),I2(JL),I3(JL))
    ZTHS(JL) = PTHS(I1(JL),I2(JL),I3(JL))
    ZCRS(JL) = PCRS(I1(JL),I2(JL),I3(JL))
    ZZT(JL) = ZT(I1(JL),I2(JL),I3(JL))
    ZZW1(JL) = ZRVSAT(I1(JL),I2(JL),I3(JL))
    ZLBDR(JL) = ZWLBDR(I1(JL),I2(JL),I3(JL))
    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
    ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
    ZZLV(JL)  = ZLV(I1(JL),I2(JL),I3(JL))
  END DO
!
!*       5.1  Compute the intermediate supersaturation mixing ratio
!
  ZZW3(:) = MAX((1.0 - ZRVT(:)/ZZW1(:)),0.0)  ! Subsaturation
!
!*       5.2  Compute the function G(T)
!
  ZZW2(:) = 1. / ( XRHOLW*((((ZZLV(:)/ZZT(:))**2)/(XTHCO*XRV)) +          & ! G
          (XRV*ZZT(:))/(XDIVA*EXP(XALPW-XBETAW/ZZT(:)-XGAMW*ALOG(ZZT(:))))))
!
!*       5.3  Compute the evaporation tendency
!
  IF (HCLOUD =='C2R2'.OR. HCLOUD=='C3R5') THEN
    ZZW2(:) = MIN( ZZW2(:) * ZZW3(:) * ZRRT(:) *        &
                (X0EVAR*ZLBDR(:)**XEX0EVAR + X1EVAR*ZRHODREF(:)**XEX2EVAR* &
                 ZLBDR(:)**XEX1EVAR),ZRRS(:) )
    ZZW2(:) = MAX(ZZW2(:),0.0)
  ELSE
    ZZW2(:) = 3.0 * XCEVAP * ZZW2(:) * (4.*XPI*XRHOLW/(3.*ZRHODREF(:)))**(2./3.) *    &
                               (ZRRT(:))**(1./3.) * (ZCRT(:))**(2./3.) * ZZW3(:)
    ZZW2(:) = MIN(ZZW2(:),ZRRS(:))        
  ENDIF
!
!*       5.4  Adjust sources
!
  ZRVS(:) = ZRVS(:) + ZZW2(:)
  ZRRS(:) = ZRRS(:) - ZZW2(:)
  ZTHS(:) = ZTHS(:) - ZZW2(:) * ZZLV(:) /                                        &
                    ( ZEXNREF(:)*(XCPD + XCPV*ZRVT(:) + XCL*(ZRCT(:) + ZRRT(:)) ) )
!
  ZW(:,:,:) = PRVS(:,:,:)
  PRVS(:,:,:) = UNPACK( ZRVS(:),MASK=GEVAP(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PRRS(:,:,:)
  PRRS(:,:,:) = UNPACK( ZRRS(:),MASK=GEVAP(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PTHS(:,:,:)
  PTHS(:,:,:) = UNPACK( ZTHS(:),MASK=GEVAP(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:)= PEVAP3D(:,:,:)
  PEVAP3D(:,:,:) = UNPACK( ZZW2(:),MASK=GEVAP(:,:,:),FIELD=ZW(:,:,:) )
! 
  IF (HCLOUD == 'KHKO') THEN
    ZZW2(:) = MIN(ZZW2(:) * ZCRT(:)/ZRRT(:),ZCRS(:))
    ZCRS(:) = ZCRS(:) - ZZW2(:)
    ZW(:,:,:) = PCRS(:,:,:)
    PCRS(:,:,:) = UNPACK( ZCRS(:),MASK=GEVAP(:,:,:),FIELD=ZW(:,:,:) )
  ENDIF
  DEALLOCATE(ZRCT)
  DEALLOCATE(ZRRT)
  DEALLOCATE(ZRVT)
  DEALLOCATE(ZCRT)
  DEALLOCATE(ZRVS)
  DEALLOCATE(ZRRS)
  DEALLOCATE(ZTHS)
  DEALLOCATE(ZCRS)  
  DEALLOCATE(ZZLV)
  DEALLOCATE(ZZT)
  DEALLOCATE(ZRHODREF)
  DEALLOCATE(ZEXNREF)
  DEALLOCATE(ZZW1)
  DEALLOCATE(ZZW2)
  DEALLOCATE(ZZW3)
  DEALLOCATE(ZLBDR)
!
END IF

IF (HCLOUD == 'C2R2'.OR. HCLOUD=='C3R5') THEN
!*       5.5  Update Nr if:  80 microns < Dr < D_h
!
  GEVAP(:,:,:) = PRRS(:,:,:)>ZRTMIN(3) .AND. PCRS(:,:,:)>ZCTMIN(3) .AND. &
                 PRCS(:,:,:)>ZRTMIN(2) .AND. PCCS(:,:,:)>ZCTMIN(2)
  WHERE (GEVAP(:,:,:))
    ZWLBDR3(:,:,:) = XLBR * PCRS(:,:,:) / (PRHODREF(:,:,:) * PRRS(:,:,:))
    ZWLBDR(:,:,:)  = ZWLBDR3(:,:,:)**XLBEXR
!
    ZWLBDC3(:,:,:) = XLBC * PCCS(:,:,:) / (PRHODREF(:,:,:) * PRCS(:,:,:))
    ZWLBDC(:,:,:)  = ZWLBDC3(:,:,:)**XLBEXC
    ZWLBDC3(:,:,:) = (XACCR1/XACCR3)*(XACCR4/ZWLBDC(:,:,:)-XACCR5)
                                                                    ! "Lambda_h"
  END WHERE
!
  GMICRO(:,:,:) = GEVAP(:,:,:) .AND. ZWLBDR(:,:,:)>ZWLBDC3(:,:,:)
                          ! the raindrops are too small, that is lower than D_h
  ZFACT = 1.2E4*XACCR1
  WHERE (GMICRO(:,:,:))
    ZWLBDC(:,:,:) = XLBR / MIN( ZFACT,ZWLBDC3(:,:,:) )**3
    ZW(:,:,:) = MIN( MAX(                                                      &
                   (PRHODREF(:,:,:)*PRRS(:,:,:) - ZWLBDC(:,:,:)*PCRS(:,:,:)) / &
                   (PRHODREF(:,:,:)*PRCS(:,:,:)/PCCS(:,:,:) - ZWLBDC(:,:,:)) , &
                    0.0 ),PCRS(:,:,:),                                         &
                          PCCS(:,:,:)*PRRS(:,:,:)/(PRHODREF(:,:,:)*PRCS(:,:,:)))
!
! Compute the percent (=1 if (ZWLBDR/XACCR1) >= 1.2E4
! of transfer with    (=0 if (ZWLBDR/XACCR1) <= (XACCR4/ZWLBDC-XACCR5)/XACCR3
!
    ZW(:,:,:) = ZW(:,:,:)*( (MIN(ZWLBDR(:,:,:),1.2E4*XACCR1)-ZWLBDC3(:,:,:)) / &
                            (                  1.2E4*XACCR1 -ZWLBDC3(:,:,:))   )
!
    ZWLBDC(:,:,:) = PCCS(:,:,:)      !temporary storage
    PCCS(:,:,:) = PCCS(:,:,:)+ZW(:,:,:)
    PCRS(:,:,:) = PCRS(:,:,:)-ZW(:,:,:)
    ZW(:,:,:)     = ZW(:,:,:) * (PRHODREF(:,:,:)*PRCS(:,:,:)/ZWLBDC(:,:,:))
    PRCS(:,:,:) = PRCS(:,:,:)+ZW(:,:,:)

    PRRS(:,:,:) = PRRS(:,:,:)-ZW(:,:,:)
  END WHERE
!
  GEVAP(:,:,:) = PRRS(:,:,:)<ZRTMIN(3) .OR. PCRS(:,:,:)<ZCTMIN(3)
  WHERE (GEVAP(:,:,:))
    PCRS(:,:,:) = 0.0
    PRRS(:,:,:) = 0.0
  END WHERE
!

ELSE ! KHKO
!*             correct negative values for rain
!   	      --------------------------------
!
  WHERE (PRRS(:,:,:)<0.) 
    PRCS(:,:,:) = PRCS(:,:,:)+PRRS(:,:,:)
    PRRS(:,:,:) = 0.
    PCRS(:,:,:) = 0.
  END WHERE
!
!*           REMOVES NON-PHYSICAL LOW VALUES
  GEVAP(:,:,:) = PRRS(:,:,:)<ZRTMIN(3) .AND. PCRS(:,:,:)< ZCTMIN(3)
  WHERE (GEVAP(:,:,:))
    PRVS(:,:,:) = PRVS(:,:,:) + PRRS(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) - PRRS(:,:,:) * ZLV(:,:,:) /                         &
     ( PEXNREF(:,:,:)*(XCPD + XCPV*PRVT(:,:,:) + XCL*(PRCT(:,:,:) + PRRT(:,:,:)) ) )
    PCRS(:,:,:) = 0.0
    PRRS(:,:,:) = 0.0
  END WHERE
ENDIF

if ( lbudget_th ) call Budget_store_end( tbudgets(NBUDGET_TH),                        'REVA', pths(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rv ) call Budget_store_end( tbudgets(NBUDGET_RV),                        'REVA', prvs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rr ) call Budget_store_end( tbudgets(NBUDGET_RR),                        'REVA', prrs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2), 'REVA', pcrs(:, :, :) * prhodj(:, :, :) )

  END SUBROUTINE C2R2_KHKO_EVAPORATION
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE C2R2_KHKO_SEDIMENTATION
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
!
INTEGER , DIMENSION(SIZE(GSEDIM)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics 
!
!-------------------------------------------------------------------------------

if ( lbudget_rc .and. osedc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rr             ) call Budget_store_init( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_sv ) then
  if ( osedc ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'SEDI', pccs(:, :, :) * prhodj(:, :, :) )
  call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2), 'SEDI', pcrs(:, :, :) * prhodj(:, :, :) )
end if
!
!*       2.1    compute the fluxes  
!
!  optimization by looking for locations where
!  the precipitating fields are larger than a minimal value only !!!
!
IF (OSEDC) PINPRC (:,:) = 0.
IF (LDEPOC) PINDEP (:,:) = 0.
!
DO JN = 1 , KSPLITR
  GSEDIM(:,:,:) = .FALSE.
  IF( OSEDC ) THEN
    GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) =                               &
            PRCT(IIB:IIE,IJB:IJE,IKB:IKE)/PTSTEP>ZRTMIN(2) .OR.     &
            (PRRT(IIB:IIE,IJB:IJE,IKB:IKE)/PTSTEP>ZRTMIN(3) .AND.   &
            PCRT(IIB:IIE,IJB:IJE,IKB:IKE)/PTSTEP>ZCTMIN(3))
  ELSE
    GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) =                               &
            PRRT(IIB:IIE,IJB:IJE,IKB:IKE)/PTSTEP>ZRTMIN(3) .AND.    &
            PCRT(IIB:IIE,IJB:IJE,IKB:IKE)/PTSTEP>ZCTMIN(3)
  END IF
!
  ISEDIM = COUNTJV( GSEDIM(:,:,:),I1(:),I2(:),I3(:))
!
    IF( JN==1 ) THEN
      IF( OSEDC ) THEN
        ZPCCT(:,:,:) = PCCT(:,:,:)
        ZPRCT(:,:,:) = PRCT(:,:,:)
        PCCS(:,:,:) = PCCS(:,:,:) * PTSTEP - PCCT(:,:,:)
        PRCS(:,:,:) = PRCS(:,:,:) * PTSTEP - PRCT(:,:,:)
      END IF
      ZPCRT(:,:,:) = PCRT(:,:,:)
      ZPRRT(:,:,:) = PRRT(:,:,:)
      PCRS(:,:,:) = PCRS(:,:,:) * PTSTEP - PCRT(:,:,:)
      PRRS(:,:,:) = PRRS(:,:,:) * PTSTEP - PRRT(:,:,:)
      DO JK = IKB , IKE
        ZW(:,:,JK) = ZTSPLITR/(PZZ(:,:,JK+1) -PZZ(:,:,JK))
      END DO
    END IF
!
  ZWSEDR(:,:,:) = 0.0
  ZWSEDC(:,:,:) = 0.0
!
  IF( ISEDIM >= 1 ) THEN
!
    ALLOCATE(ZRHODREF(ISEDIM))
    DO JL = 1,ISEDIM
      ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
    END DO
!
    ALLOCATE(ZZW1(ISEDIM)) 
    ALLOCATE(ZZW2(ISEDIM)) 
    ALLOCATE(ZZW3(ISEDIM)) 
!
!*       2.21   for cloud
!
    ZZW1(:) = 0.0
    ZZW2(:) = 0.0
    ZZW3(:) = 0.0
!
    IF( OSEDC.AND.MAXVAL(PRCS(:,:,:))>0.0 ) THEN
      ALLOCATE(ZRCT(ISEDIM))
      ALLOCATE(ZCCT(ISEDIM))
      ALLOCATE(ZLBDC(ISEDIM))
      DO JL = 1,ISEDIM
        ZRCT(JL) = ZPRCT(I1(JL),I2(JL),I3(JL))
        ZCCT(JL) = ZPCCT(I1(JL),I2(JL),I3(JL))
        ZLBDC(JL) = ZWLBDC(I1(JL),I2(JL),I3(JL))
      END DO
      WHERE( ZRCT(:)>XRTMIN(2) )
        ZZW3(:) = ZRHODREF(:)**(-XCEXVT) * ZLBDC(:)**(-XDC)
        ZZW1(:) = XFSEDRC * ZRCT(:) * ZZW3(:) * ZRHODREF(:)
        ZZW2(:) = XFSEDCC * ZCCT(:) * ZZW3(:)
      END WHERE
      ZWSEDR(:,:,:) = UNPACK( ZZW1(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
      ZWSEDC(:,:,:) = UNPACK( ZZW2(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
      DEALLOCATE(ZRCT)
      DEALLOCATE(ZCCT)
      DEALLOCATE(ZLBDC)
    END IF
!             
   END IF
!        
   IF( OSEDC ) THEN
     DO JK = IKB , IKE
        ZPRCT(:,:,JK) = ZPRCT(:,:,JK) + ZW(:,:,JK)*    &
                      (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
        ZPCCT(:,:,JK) = ZPCCT(:,:,JK) + ZW(:,:,JK)*    &
                      (ZWSEDC(:,:,JK+1)-ZWSEDC(:,:,JK))
     END DO
!             
     IF( JN.EQ.1 ) THEN
        PINPRC(:,:) = ZWSEDR(:,:,IKB)/XRHOLW                           ! in m/s
     END IF
   END IF
!
!*       2.22   for drizzle
!
  ZWSEDR(:,:,:) = 0.0
  ZWSEDC(:,:,:) = 0.0
  IF( ISEDIM >= 1 ) THEN
    ZZW1(:) = 0.0
    ZZW2(:) = 0.0
!
    IF( MAXVAL(PRRS(:,:,:))>0.0 ) THEN
      ALLOCATE(ZRRT(ISEDIM)) 
      ALLOCATE(ZCRT(ISEDIM))
      ALLOCATE(ZZVRR(ISEDIM))
      ALLOCATE(ZZVCR(ISEDIM))
      DO JL = 1,ISEDIM
        ZRRT(JL) = ZPRRT(I1(JL),I2(JL),I3(JL))
        ZCRT(JL) = ZPCRT(I1(JL),I2(JL),I3(JL))
        ZZVRR(JL) = ZVRR(I1(JL),I2(JL),I3(JL))
        ZZVCR(JL) = ZVCR(I1(JL),I2(JL),I3(JL))
      END DO
      WHERE (ZRRT(:)>XRTMIN(3) )
        ZZW1(:) = ZZVRR(:) * ZRRT(:) * ZRHODREF(:)
        ZZW2(:) = ZZVCR(:) * ZCRT(:)
      END WHERE
      ZWSEDR(:,:,:) = UNPACK( ZZW1(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
      ZWSEDC(:,:,:) = UNPACK( ZZW2(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
!
      DEALLOCATE(ZRRT)
      DEALLOCATE(ZCRT)
      DEALLOCATE(ZZVRR)
      DEALLOCATE(ZZVCR)
!
    END IF
!
    DEALLOCATE(ZRHODREF)
    DEALLOCATE(ZZW1)
    DEALLOCATE(ZZW2)
    DEALLOCATE(ZZW3)
!
  END IF
!
!*       2.3     update the rain tendency
!
  DO JK = IKB , IKE
        ZPRRT(:,:,JK) = ZPRRT(:,:,JK) + ZW(:,:,JK)* &
                      (ZWSEDR(:,:,JK+1)-ZWSEDR(:,:,JK))/PRHODREF(:,:,JK)
        ZPCRT(:,:,JK) = ZPCRT(:,:,JK) + ZW(:,:,JK)* &
                      (ZWSEDC(:,:,JK+1)-ZWSEDC(:,:,JK))
  END DO
!
!*       2.4     compute the explicit accumulated precipitations
!         
  IF( JN.EQ.1 ) THEN
      PINPRR(:,:) = ZWSEDR(:,:,IKB)/XRHOLW                           ! in m/s
      PINPRR3D(:,:,:) = ZWSEDR(:,:,:)/XRHOLW                           ! in m/s
  END IF
!
  IF( JN==KSPLITR ) THEN
      IF( OSEDC ) THEN
        PRCS(:,:,:) = ( PRCS(:,:,:) + ZPRCT(:,:,:) ) / PTSTEP
        PCCS(:,:,:) = ( PCCS(:,:,:) + ZPCCT(:,:,:) ) / PTSTEP
      END IF
      PRRS(:,:,:) = ( PRRS(:,:,:) + ZPRRT(:,:,:) ) / PTSTEP
      PCRS(:,:,:) = ( PCRS(:,:,:) + ZPCRT(:,:,:) ) / PTSTEP
  END IF
!   
 IF ( OSEDC .AND. tpfile%lopened ) THEN
  TZFIELD = TFIELDMETADATA(        &
    CMNHNAME   = 'SEDFLUXC',       &
    CSTDNAME   = '',               &
    CLONGNAME  = 'SEDFLUXC',       &
    CUNITS     = '',               &
    CDIR       = 'XY',             &
    CCOMMENT   = 'X_Y_Z_SEDFLUXC', &
    NGRID      = 1,                &
    NTYPE      = TYPEREAL,         &
    NDIMS      = 3,                &
    LTIMEDEP   = .TRUE.            )
  CALL IO_Field_write(TPFILE,TZFIELD,ZWSEDC)
  !
  TZFIELD = TFIELDMETADATA(        &
    CMNHNAME   = 'SEDFLUXR',       &
    CSTDNAME   = '',               &
    CLONGNAME  = 'SEDFLUXR',       &
    CUNITS     = '',               &
    CDIR       = 'XY',             &
    CCOMMENT   = 'X_Y_Z_SEDFLUXR', &
    NGRID      = 1,                &
    NTYPE      = TYPEREAL,         &
    NDIMS      = 3,                &
    LTIMEDEP   = .TRUE.            )
  CALL IO_Field_write(TPFILE,TZFIELD,ZWSEDR)
 END IF
END DO
!
!*       2.5     budget storage
!
if ( lbudget_rc .and. osedc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'SEDI', prcs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_rr             ) call Budget_store_end( tbudgets(NBUDGET_RR), 'SEDI', prrs(:, :, :) * prhodj(:, :, :) )
if ( lbudget_sv ) then
  if ( osedc ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'SEDI', pccs(:, :, :) * prhodj(:, :, :) )
  call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 2),              'SEDI', pcrs(:, :, :) * prhodj(:, :, :) )
end if
!
!*       2.6  DROPLET DEPOSITION AT THE 1ST LEVEL ABOVE GROUND
!
IF (LDEPOC) THEN
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC),                        'DEPO', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_init( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'DEPO', pccs(:, :, :) * prhodj(:, :, :) )

  GDEP(:,:) = .FALSE.
  GDEP(IIB:IIE,IJB:IJE) =    PRCS(IIB:IIE,IJB:IJE,2) >0 .AND. &
                                  PCCS(IIB:IIE,IJB:IJE,2) >0
  WHERE (GDEP)
     PRCS(:,:,2) = PRCS(:,:,2) - XVDEPOC * PRCT(:,:,2) / ( PZZ(:,:,3) - PZZ(:,:,2))
     PCCS(:,:,2) = PCCS(:,:,2) - XVDEPOC * PCCT(:,:,2) / ( PZZ(:,:,3) - PZZ(:,:,2))
     PINPRC(:,:) = PINPRC(:,:) + XVDEPOC * PRCT(:,:,2) * PRHODREF(:,:,2) /XRHOLW             
     PINDEP(:,:) = XVDEPOC * PRCT(:,:,2) * PRHODREF(:,:,2) /XRHOLW
  END WHERE

  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC),                        'DEPO', prcs(:, :, :) * prhodj(:, :, :) )
  if ( lbudget_sv ) call Budget_store_end( tbudgets(NBUDGET_SV1 - 1 + nsv_c2r2beg + 1), 'DEPO', pccs(:, :, :) * prhodj(:, :, :) )
END IF

  END SUBROUTINE C2R2_KHKO_SEDIMENTATION
!-------------------------------------------------------------------------------
!
END SUBROUTINE RAIN_C2R2_KHKO
